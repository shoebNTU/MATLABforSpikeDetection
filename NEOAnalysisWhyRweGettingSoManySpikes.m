%% this code is being written to analyze NEO, as in why is our NEO based implementation with quantized information
% is yielding spikes at almost all instances..there will be two NEO based
% spikes returned -- first one without any approximation and the other one
% with approximation..Lets see what is the difference..
%initial loading
clear all;
close all;
clc;
laptop = 1;

if (laptop)
    loadFrom = 'sdshaikh';
else
    loadFrom = 'my';
end

date = 151026;
session = 6;

%loading data
load(['C:\Users\' loadFrom '\Google Drive\ELMforChairM\MATLABInterfaceCodes\ExtractedRawData\' num2str(date) ...
    'Sess' num2str(session) '.mat']);

sr = (40/14/2)*10^6/10/11;
[b,a]=ellip(2,0.1,40,[300 3000]*2/sr); %filter design

% in I2R's code, data is zero meaned for every channel, filtered before threshold based spike detection
% in real time, zero mean is not feasible -- so ignoring it for now
%%
%threshold calculation block time is "35 secs"
channelArray = 1:100;
ReadTime = 35;
dataThr = data(:,1:ceil(ReadTime*sr)); %selecting the columns corresponding to the read-time (~ 10% of the total time)

Thr = zeros(length(channelArray),1); %threshold array dimensioned as (channelNos. x 1)
ThrNormto255 = zeros(length(channelArray),1);
for channelNo = channelArray
    dataThrFil = round(filter(b,a,round(255*double(dataThr(channelNo,:))/512))); %filtering the data
    %scale the data, the min should be 0, max - 255
%     dataThrFil = round(dataThr(channelNo,:));
    RangeOfChannel = max(dataThrFil) - min(dataThrFil); %Range of Channel
    dataThrFilNormto255 = round((dataThrFil/RangeOfChannel)*255 - (min(dataThrFil)/RangeOfChannel)*255);%eliminating -ve and range:0to255
    dataThrFilNEO = NEOfn(dataThrFil);
    dataThrFilNEONormto255 = NEOfn(dataThrFilNormto255);
    Thr(channelNo,1) = median((dataThrFilNEO));
    ThrNormto255(channelNo,1) = mean(dataThrFilNEONormto255);
end

CNEO = 8;
ThrwithC = round(CNEO*Thr); %multiplication with C
ThrwithCNormto255 = round(CNEO*ThrNormto255);
thchk = [ThrwithC ThrwithCNormto255];
%%
channelOfPlot = 1;
NoOfDataPoints = 78000:79000;
plot(NoOfDataPoints,dataThr(1,NoOfDataPoints));
hold on;
%%
plot(NoOfDataPoints,dataThrFil(1,NoOfDataPoints));
hold on;
plot(NoOfDataPoints,dataThrFilNEO(channelOfPlot,NoOfDataPoints));
%%

plot(NoOfDataPoints,dataThrFilNEO(channelOfPlot,NoOfDataPoints));
% axis([1 1000 -500 500]);
hold on;
plot(NoOfDataPoints,dataThrFilNormto255(channelOfPlot,NoOfDataPoints));
hold on;
plot(NoOfDataPoints,ThrwithC(channelOfPlot)*ones(1,length(NoOfDataPoints)));
hold on;
plot(NoOfDataPoints,ThrwithCNormto255(channelOfPlot)*ones(1,length(NoOfDataPoints)));

%% After passing the threshold, need to pass the data for the entire session
% send the 8 bit normalized data
dataFilt = zeros(length(channelArray),size(data,2));
% maxChannel = [];
for channelNo = channelArray
    filteredData = (filter(b,a,255*double((data(channelNo,:)))/512)); %filter the data
    
%     RangeOfChannel = max(filteredData) - min(filteredData);
%     filteredDataNormto255 = (filteredData/RangeOfChannel)*255 - (min(filteredData)/RangeOfChannel)*255;
    dataFilt(channelNo,:) = round(filteredData);
%     maxChannel = [maxChannel; max()];
%     filteredData(find(filteredData == 254)) = 253; %randomly setting 254 as the stop bit and replacing all 254's by 253
end

% clear data;
% got the 8 bit filtered, normalized data
%%
%which instances to pass:

% the no. of instances will be determined from t_0's,90's..from these t's
% we will get the instances of time of interest, meaning 500 ms before the
% corresponding time -- just to rewind the timings are obtained from the
% positions of the joystick with std deviation less than 200 corresponding
% to the recorded instances of time..check training_allThr fn for further
% info. Over here, directly loading the training file with timing info

%load the training file to get the timings
load(['C:\Users\' loadFrom '\Google Drive\ELMforChairM\cartman spike data (NG)\Training_' ...%ModifiedCodesforELM .. cartman spike data (NG)
    num2str(date) '_Session' num2str(session) '_allThr.mat'])
t_0=ModelParam.t_0;
t_90=ModelParam.t_90;
t_180=ModelParam.t_180;
t_stop=ModelParam.t_stop;
t=sortrows([t_0;t_90;t_180;t_stop]); %t has all the timing instances (in ascending order) where we have a target

tend = round(t(3)*sr);
tbegin =  round((t(3)-0.5)*sr);
noOfInstances = tbegin:(tbegin+50);

%% set random threshold here and check if the microcontroller's output is correct..Once this is done, pass different thresholds...
% for different channels and cross-check, by receiving spikes

%%
%validating in MATLAB - NEO spikes
dataSelect = dataFilt(channelArray,noOfInstances(1):(noOfInstances(end)));
zeroPad = zeros(length(channelArray),1); %this will be used it introduce delays

dataPlus1 = [zeroPad dataSelect(:,1:size(dataSelect,2)-1)];
dataMinus1 = [zeroPad dataSelect(:,3:size(dataSelect,2)) zeroPad];

NEOOutput = power(dataSelect,2) - dataPlus1.*dataMinus1;

% generate a spike array similar to the one received from MSP430 in order
% to compare the results
NEOOutputMSPForm = zeros(size(NEOOutput,1),size(NEOOutput,2)-2); %ignore first and last columns

for channelNo = channelArray
    indexNEOgreaterthanThr = find(NEOOutput(channelNo,2:(end-1))>ThrwithC(channelNo,1));
    NEOOutputMSPForm(channelNo,indexNEOgreaterthanThr) = 1;
end
%% this code was to check if the channel nos. are noted down correctly in the address pointer.
% chkStoredChannelNos = [];
% for jspk = 1:size(NEOOutputMSPForm,2)
%     for ispk = 1:length(channelArray)
%         if(NEOOutputMSPForm(ispk,jspk))
%             chkStoredChannelNos = [chkStoredChannelNos (ispk-1)];
%             break;
%         end
%     end
% end
%%
if(isequal(NEOOutputMSPForm,spkRxdTemporal))
    print = 'Outputs Match'
end

