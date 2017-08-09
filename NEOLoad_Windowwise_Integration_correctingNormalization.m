%% this code does NEO based spike detection
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
dataThr = data(:,1:ceil(ReadTime*sr)); %selecting the columns corresponding to the read-time (~ 10% of the total time = 35s)

Thr = zeros(length(channelArray),1); %threshold array dimensioned as (channelNos. x 1)
for channelNo = channelArray
    dataThrNorm = 255*double(dataThr(channelNo,:))/512; %normalizing before filtering
    dataThrFil = round(filter(b,a,dataThrNorm)); %filtering the data
    dataThrFilNEO = NEOfn(dataThrFil);
    Thr(channelNo,1) = mean((dataThrFilNEO));
end

CNEO = 8;
ThrwithC = round(CNEO*Thr); %multiplication with C
%%
% send the threshold%
delete(instrfindall);
uart=serial('COM7','BaudRate',115200 ,'DataBits',8,'StopBits',1); %
fopen(uart);

% ThrwithC(1,1) = -1;

for thresSend = 1:length(ThrwithC)
    fwrite(uart,ThrwithC(thresSend,1),'int8');
end

fclose(uart);
clear uart;
%% After passing the threshold, need to pass the data for the entire session
% send the 8 bit normalized data
dataFiltNormto255 = zeros(length(channelArray),size(data,2));
for channelNo = channelArray
    filteredData =round(255*filter(b,a,double((data(channelNo,:))))/512); %filter the data
    dataFiltNormto255(channelNo,:) = filteredData;
end
%%
dataFiltNormto255 (dataFiltNormto255 < -127) = -127; %stopbit is automatically taken care of here
dataFiltNormto255 (dataFiltNormto255 > 127) = 127;

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

%% Communication with MSP430
% copying the UART block as it is -
delete(instrfindall);
uart=serial('COM7','BaudRate',115200 ,'DataBits',8,'StopBits',1); %
% fopen(uart);
spkRxdTemporal = [];
fopen(uart);
for  jtime = 1:(length(noOfInstances)-2) %here the timing instances corresponding to each window can be set
    
    for ichannel = channelArray
        if(jtime == 1)
            pause(0.005);
            fwrite(uart,dataFiltNormto255(ichannel,noOfInstances(1):noOfInstances(3)),'int8');  %dataOneChNormto255(1,1:3)
        else
            pause(0.005);
            fwrite(uart, dataFiltNormto255(ichannel,noOfInstances(2+jtime)),'int8');
        end
    end
    
    %     spkRxd = [];
    %     %     for readChannel = channelArray
    %     %         spkRxd  = [spkRxd; fread(uart,1)];
    %     %     end
    %     spkRxd = fread(uart,length(channelArray));
    %     spkRxdTemporal = [spkRxdTemporal spkRxd];
    
end
%
pause(1);
fwrite(uart, -128,'int8'); %as of now the stop bit is 254
% fclose(uart);
%before just closing everything off write like a reset value
% fclose(uart);
% Read the output here..in the above section -
% pause(5);
freqHiddenNeuronCum=[];
noOfOutputNeurons=128;
countAcc = [];

for noOfInstancesRead=1 %there should be a 512 constraint here
    b=fread(uart,2*noOfOutputNeurons);
    
    k=1;
    for(i=1:2:2*noOfOutputNeurons-1)
        %if(b(i,1)==1)
        b(i,1)=256*b(i,1);
        %end
        c(k,1)=b(i,1)+b(i+1,1);
        k=k+1;
    end
    countAcc=[countAcc; (fliplr(c))'];
    % freqHiddenNeuron = (1/10)*fliplr(c);
end
fclose(uart);
%%
%validating in MATLAB - NEO spikes

dataSelect = dataFiltNormto255(channelArray,noOfInstances(1):(noOfInstances(end)));
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
%% this part is for generating timingcycles(timerB input values) and channel addresses
%These can be passed serially to emulate the same NEO based spike detection
%
NEOOutputMSPFormSum = sum(NEOOutputMSPForm,1);
cycleTrackerInit = find(NEOOutputMSPFormSum>0);
cycleTrackerInit = cycleTrackerInit';
cycleTracker = cycleTrackerInit - [0; cycleTrackerInit(1:(end-1))]; %finding the timing cycles in the incremental difference form
%going for address detection
channelAddress = [];
for channelTracker = cycleTrackerInit'
channelAddressChk = find(NEOOutputMSPForm(:,channelTracker)>0); %as we sweep through the timing instances of ...
% spike occurrence, we check which addresses correspond to the spike
channelAddress = [channelAddress; channelAddressChk(1)-1]; %we consider only the first while discarding the rest
end

%%
if(isequal(NEOOutputMSPForm,spkRxdTemporal))
    print = 'Outputs Match'
end

