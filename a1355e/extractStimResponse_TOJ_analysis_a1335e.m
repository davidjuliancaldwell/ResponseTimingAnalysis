% this is from my z_constants
close all;clear all;clc
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\TOJ';
sid = SIDS{5};

for s = 1:2
    % load in data
    folder_data = strcat(DATA_DIR);
    if s == 1
        load(fullfile(folder_data,'TOJ-1.mat'))
        block = '1';
    elseif s == 2
        load(fullfile(folder_data,'TOJ-2.mat'))
        block = '2';
        length1st = size(tact,1);
    end
    
    % load in data of interest
    if s == 1
        stim = Stim.data;
    else
        stim = [stim; Stim.data];
    end
    fs_stim = Stim.info.SamplingRateHz;
    
    clear Stim
    
    % the eco data is crashing it right now
    % eco1 = ECO1.data;
    % eco2 = ECO2.data;
    % eco3 = ECO3.data;
    fs_data = ECO1.info.SamplingRateHz;
    %
    %
    clear ECO1 ECO2 ECO3
    % data = [eco1 eco2 eco3];
    % clear  eco1 eco2 eco3
    
    if s ==1
        tact = Tact.data;
    else
        tact = [tact; Tact.data];
    end
    fs_tact = Tact.info.SamplingRateHz;
    clear Tact
end
samplesOfPulse = round(2*fs_stim/1e3);
% build a burst table with the timing of stimuli
preStim = round(fs_stim*1);
[trainTimes] = find(tact(:,8)==1)- preStim;
logicalMask = logical(ones(size(trainTimes)));
logicalMask(14) = 0;
trainTimes = trainTimes(logicalMask);
starts = trainTimes;
ends = trainTimes + samplesOfPulse;
delay = round(0.2867*fs_stim/1e3);
% extract data
% try and account for delay for the stim times
stimTimes = starts+delay;
trainTimes=stimTimes;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * fs_data); % pre time in sec
postsamps = round(1 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after

% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);

% get that button press
tact(tact(:,2) >= 0.009,2) = 0.009;
[buttonPksTemp,buttonLocsTemp] = findpeaks(tact(:,2),fs_tact,'minpeakdistance',2,'Minpeakheight',0.008);
tact(tact(:,2) < 0.009,2) = 0;
tact(:,2) = tact(:,2)*1000;

%%
figure
hold on
t = [0:length(tact(:,1))-1]/fs_stim;
plot(t,stim(:,1),'linewidth',2)
plot(t,tact(:,1),'linewidth',2)
plot(t,tact(:,4),'linewidth',2)
plot(t,tact(:,2),'linewidth',2)
legend({'stim','tactor','audio train','button'})
%vline(trainTimes/fs_stim)

% QUANTIFY RXN TIME TO CORTICAL STIM
sampsEnd = round(3*fs_stim);

epochedTactor = squeeze(getEpochSignal(tact(:,1),trainTimes,trainTimes+sampsEnd));
epochedAudio = squeeze(getEpochSignal(tact(:,4),trainTimes,trainTimes+sampsEnd));
epochedStim = squeeze(getEpochSignal(stim(:,1),trainTimes,trainTimes+sampsEnd));
epochedButton = squeeze(getEpochSignal(tact(:,2),trainTimes,trainTimes+sampsEnd)); 

t = [-preStim:sampsEnd-preStim-1]/fs_tact;
t_samps = [-preStim:sampsEnd-preStim-1];
% epoched button press
%%
numTrials = size(epochedAudio,2);
[p,n]=numSubplots(numTrials);
figure
whichPerceived = {'stim','same','stim','stim','tactor','same','tactor','stim','stim','stim','tactor','same','stim',...
    'stim','stim','tactor','stim','tactor','stim','tactor','stim'};
whichPerceiveMoresame = {'stim','same','stim','stim','tactor','same','tactor','stim','same','same','tactor','same','stim',...
    'stim','same','tactor','stim','tactor','stim','tactor','stim'};

for i = 1:numTrials
    subplot(p(1),p(2),i)
    hold on
    plot(t,epochedStim(:,i),'linewidth',2)
    plot(t,epochedTactor(:,i),'linewidth',2)
    plot(t,epochedAudio(:,i),'linewidth',2)
    plot(t,epochedButton(:,i),'linewidth',2)
    title(whichPerceived{i})
    xlabel('time (s)')
    ylim([-6 6])
end
legend({'stim','tactor','audio train','button'})

%%
for j = 1:numTrials
    
    [buttonPksTemp,buttonLocsTemp] = findpeaks((epochedButton(t>0,j)),t(t>0),'NPeaks',1,'Minpeakheight',3);
    [buttonPksTempSamps,buttonLocsTempSamps] = findpeaks((epochedButton(t_samps>0,j)),t_samps(t_samps>0),'NPeaks',1,'Minpeakheight',0.008);
    
    sprintf(['button ' num2str(buttonLocsTemp)])
    sprintf(['button ' num2str(buttonLocsTempSamps)])
    
    [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),t,'NPeaks',1,'Minpeakheight',1);
    [tactorPksTempSamps,tactorLocsTempSamps] = findpeaks((epochedTactor(:,j)),t_samps,'NPeaks',1,'Minpeakheight',1);
    
    sprintf(['tactor ' num2str(tactorLocsTemp)])
    sprintf(['tactor ' num2str(tactorLocsTempSamps)])
    
    [stimPksTemp,stimLocsTemp] = findpeaks((epochedStim(:,j)),t,'NPeaks',1,'Minpeakheight',1);
    [stimPksTempSamps,stimLocsTempSamps] = findpeaks((epochedStim(:,j)),t_samps,'NPeaks',1,'Minpeakheight',1);
    
    sprintf(['stim train ' num2str(stimLocsTemp)])
    sprintf(['tactor ' num2str(stimLocsTempSamps)])    
    
    if isempty(buttonPksTemp)
        buttonPksTemp = NaN;
        buttonLocsTemp = NaN;
        buttonPksTempSamps = NaN;
        buttonLocsTempSamps = NaN;
    end
    
    if isempty(tactorPksTemp)
        tactorPksTemp = NaN;
        tactorLocsTemp = NaN;
        tactorPksTempSamps = NaN;
        tactorLocsTempSamps = NaN;
    end
    
    if isempty(stimPksTemp)
       stimPksTemp = NaN;
        stimLocsTemp = NaN;
        stimPksTempSamps = NaN;
        stimLocsTempSamps = NaN;
    end
    
    buttonPksVec(j) = buttonPksTemp;
    buttonLocsVec(j) = buttonLocsTemp;
    
    % do samples too
    buttonPksVecSamps(j) = buttonPksTempSamps;
    buttonLocsVecSamps(j) = buttonLocsTempSamps;
    
    tactorPksVec(j) = tactorPksTemp;
    tactorLocsVec(j) = tactorLocsTemp;
    
    tactorPksVecSamps(j) = tactorPksTempSamps;
    tactorLocsVecSamps(j) = tactorLocsTempSamps;
    
    stimPksVec(j) = stimPksTemp;
    stimLocsVec(j) = stimLocsTemp;
    
    stimPksVecSamps(j) = stimPksTempSamps;
    stimLocsVecSamps(j) = stimLocsTempSamps;
    
end

%%
tactorStimDiff = tactorLocsVec - stimLocsVec;

responseTimes = buttonLocsVec - min(tactorLocsVec,stimLocsVec);

saveIt = 1;

if saveIt
    current_direc = pwd;
    
    save(fullfile(current_direc, [sid '_TOJ_matlab.mat']),'tactorStimDiff','responseTimes','t','trainTimes',...,
        'fs_stim','epochedButton','epochedTactor','epochedAudio','epochedStim','sampsEnd');    
end
