%% load in subject

% this is from my z_constants

sid = SIDS{3};

% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, '693ffd'))
    folder_data = strcat(DATA_DIR,'\693ffd');
    
    if s == 1
        load(fullfile(folder_data,'ReactionTime_693ffd-2.mat'))
        block = '1';
    elseif s == 2
        load(fullfile(folder_data,'ReactionTime_693ffd-4.mat'))
        block = '2';
    end
    
end


%% load in data of interest

stim = Stim.data;
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

sing = Sing.data;
fs_sing = Sing.info.SamplingRateHz;

clear Sing

tact = Tact.data;
fs_tact = Tact.info.SamplingRateHz;
clear Tact

valu = Valu.data;
fs_valu = Valu.info.SamplingRateHz;

clear Valu

%% figure out stim times
% vector of condition type - for first subject, looks like condition type
% is what was used , rather than test_condition
if s == 2
    condType= dlmread('rxnTime_condition_2_modified_9_23_2016.txt');
    train = dlmread('rxnTime_stimTrainDelivery_2.txt');
else
    condType= dlmread('rxnTime_condition_1_modified_9_23_2016.txt');
    train = dlmread('rxnTime_stimTrainDelivery_1.txt');
end

% stim cue from file
stimFromFile = tact(:,3);

% button press, start being onset of stimulation marker

trainTimes = find(stimFromFile~=0);

% shrink condition type to be 140
% as the experiment went S1 (110 stims), tactor (30 stims) , then break
condType = condType(1:140);

% pick condition type where stimulation was delivered

uniqueCond = unique(condType);
trainTimesTotal = {};

for i = 1:length(uniqueCond)
    trainTimesTotal{i} = trainTimes(condType==uniqueCond(i));
end

% just in case - keep these around
trainTimesCond_noStim = trainTimes(condType==0);
trainTimesCond_offTarget = trainTimes(condType==1);
trainTimesCond_tactor = trainTimes(condType==-1);
trainTimesCond_100 = trainTimes(condType==2);
trainTimesCond_200 = trainTimes(condType==3);
trainTimesCond_400 = trainTimes(condType==4);
trainTimesCond_800 = trainTimes(condType==5);

%% plot stim
%
figure
hold on
for i = 1:size(stim,2)
    
    t = (0:length(stim)-1)/fs_stim;
    subplot(3,2,i)
    plot(t*1e3,stim(:,i))
    title(sprintf('Channel %d',i))
    
end

xlabel('Time (ms)')
ylabel('Amplitude (V)')

subtitle('Stimulation Channels')



%% Sing looks like the wave to be delivered, with amplitude in uA


% 1st stim channel
Sing1 = sing(:,1);
% 2nd stim channel
Sing2 = sing(:,4);

samplesOfPulse = round(2*fs_stim/1e3);

% build a burst table with the timing of stimuli
bursts = [];
bursts(1,:) = condType;
bursts(2,:) = trainTimes;
bursts(3,:) = trainTimes + samplesOfPulse;

stims1 = squeeze(getEpochSignal(Sing1,(bursts(2,condType>=2)-1),(bursts(3,condType>=2))+1));
t = (0:size(stims1,1)-1)/fs_sing;
t = t*1e3;
figure
plot(t,stims1,'b','linewidth',2)
xlabel('Time (ms)');
ylabel('Current to be delivered (mA)')
ylim([(min(stims1(:))-100) (max(stims1(:))+100)])
title('Current to be delivered for all trials  on 1st channel')

stims2 = squeeze(getEpochSignal(Sing2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stims1,1)-1)/fs_sing;
t = t*1e3;
figure
plot(t,stims2,'b','linewidth',2)
xlabel('Time (ms)');
ylabel('Current to be delivered (mA)')
ylim([(min(stims2(:))-100) (max(stims2(:))+100)])
title('Current to be delivered for all trials  on 2nd channel')

% delay loks to be 0.2867 ms from below.

%% Plot stims with info from above

% 1st stimulation channel
stim1 = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,condType>=2)-1),(bursts(3,condType>=2))+1));
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim1Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Finding the delay between current output and stim delivery - 1st stim channel')

% hold on

plot(t,stims1)

% get the delay in stim times

delay = round(0.2867*fs_stim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType>=2)-1+delay;
stimTimesEnd = bursts(3,condType>=2)-1+delay;
stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim1Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Stim voltage monitoring with delay added in - 1st stim channel ')

% 2nd stimulation channel
stim2= stim(:,4);
stim2Epoched = squeeze(getEpochSignal(stim2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stim2Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim2Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Finding the delay between current output and stim delivery - 2nd stim channel')

% hold on

plot(t,stims2)

% get the delay in stim times

delay = round(0.2867*fs_stim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType==1)-1+delay;
stimTimesEnd = bursts(3,condType==1)-1+delay;
stim2Epoched = squeeze(getEpochSignal(stim2,stimTimesBegin,stimTimesEnd+5));
t = (0:size(stim2Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim2Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Stim voltage monitoring with delay added in - 2nd stim channel')

%% extract data
% try and account for delay for the stim times
stimTimes = bursts(2,:)+delay;
trainTimes=stimTimes;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * fs_data); % pre time in sec
postsamps = round(1 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after


% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);

%% look at tactor

tactorData = tact(:,1);
t_tact = (0:length(tactorData)-1)/fs_tact;
figure
plot(t_tact,tactorData);

title('tactor data')

% look at button press

buttonData = tact(:,2);
t_button = (0:length(buttonData)-1)/fs_tact;
figure
plot(t_button,buttonData);

title('button data')

% look at stim from file saved

t_stimFile = (0:length(stim)-1)/fs_tact;
figure
plot(t_stimFile,stimFromFile);
title('stim from file')

% look at all 3 + stim waveform

figure
ax1 = subplot(5,1,1);
plot(t_tact,tactorData)
title('tactor data')

ax2 = subplot(5,1,2);
plot(t_button,buttonData);
title('button data')

ax3 = subplot(5,1,3);
plot(t_stimFile,stimFromFile);
title('stim from file')

% assuming stim1 here is the channel where stim was being delivered
ax4 = subplot(5,1,4);
plot(t_stimFile,stim1);
title('S1 Stim Channel')

%
ax5 = subplot(5,1,5);
plot(t_stimFile,stim2);
title('Off Target Stim Channel')

%link axis
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

%% for 1st subject - only look at parts where t > 50 for sensory stim, t > 12 for tactor  (t_begin)

% t_begin = 10 ;
%
% buttonData = buttonData(t_button>t_begin);
% tactorData = tactorData(t_tact>t_begin);
% stimFromFile = stimFromFile(t_stimFile>t_begin);
% stim1 = stim1(t_stimFile>t_begin);
% stim2 = stim2(t_stimFile>t_begin);
%
% t_button = t_button(t_button>t_begin);
% t_tact = t_tact(t_tact>t_begin);
% t_stimFile = t_stimFile(t_stimFile>t_begin);

% set respLo and respHi, which are the values which for less or greater
% than rxn times aren't analyzed

% started with respLo - 0.150, try 0.100 
respLo = 0.100;
respHi = 1;

%% 9-13-2016 - start quantifying data

% find button data peaks

% set above certain threshold to 0.009
buttonDataClip = buttonData;
buttonDataClip(buttonData >= 0.009) = 0.009;
[buttonPks,buttonLocs] = findpeaks(buttonDataClip,t_button,'MinpeakDistance',2,'Minpeakheight',0.008);

figure
findpeaks(buttonDataClip,t_button,'MinpeakDistance',2,'Minpeakheight',0.008);
hold on
plot(t_stimFile,stimFromFile,'g');
plot(t_stimFile,stim1,'r');
%plot(t_stimFile,stim2,'m');

legend({'Button Data','Button Press Onset Peaks','Stimulation Times From File','S1 stim output','Off Target Stim Output'})


% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)
%% QUANTIFY RXN TIME TO CORTICAL STIM

sampsEnd = round(3.5*fs_stim);

% epoched button press

epochedButton = {};

% the last trial of the epoched button press for the null condition is
% clipped - so omit that one

trainTimesTotal{2}(end) = [];

for i = 1:length(uniqueCond)
    epochedButton{i} = squeeze(getEpochSignal(buttonDataClip,trainTimesTotal{i},(trainTimesTotal{i} + sampsEnd)));
end

epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesTotal{1},trainTimesTotal{1} + sampsEnd));

t_epoch = [0:size(epochedTactor,1)-1]/fs_stim;
t_epoch_samps = [0:size(epochedTactor,1)-1];


buttonPks = {};
buttonLocs = {};

buttonPksSamps = {};
buttonLocsSamps = {};

buttonPksTempVec = [];
buttonLocsTempVec = [];

buttonPksTempVecSamps = [];
buttonLocsTempVecSamps = [];
%%
for i = 1:length(uniqueCond)
    
    % for stimulation condititions
    if uniqueCond(i)~=-1
        for j = 1:length(trainTimesTotal{i})
            [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton{i}(:,j),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
            [buttonPksTempSamps,buttonLocsTempSamps] = findpeaks(epochedButton{i}(:,j),t_epoch_samps,'NPeaks',1,'Minpeakheight',0.008); % get sample number DJC 10-12-2017

            if isempty(buttonPksTemp)
                buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
                buttonPksTempSamps = NaN;
                buttonLocsTempSamps = NaN;
            end
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
                        
            % do samples too 
            buttonPksTempVecSamps(j) = buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
        end
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        buttonPksSamps{i} = buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
        % for tactor target condition
    elseif uniqueCond(i)==-1
        for j = 1:length(trainTimesTotal{i})
            
            [buttonPksTemp,buttonLocsTemp] = findpeaks((epochedButton{i}(t_epoch>1,j)),t_epoch(t_epoch>1),'NPeaks',1,'Minpeakheight',0.008);
            [buttonPksTempSamps,buttonLocsTempSamps] = findpeaks((epochedButton{i}(t_epoch_samps>24415,j)),t_epoch_samps(t_epoch_samps>24415),'NPeaks',1,'Minpeakheight',0.008);

            sprintf(['button ' num2str(buttonLocsTemp)])
            sprintf(['button ' num2str(buttonLocsTempSamps)])

            [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),t_epoch,'NPeaks',1,'Minpeakheight',1);
            [tactorPksTempSamps,tactorLocsTempSamps] = findpeaks((epochedTactor(:,j)),t_epoch_samps,'NPeaks',1,'Minpeakheight',1);

            sprintf(['tactor ' num2str(tactorLocsTemp)])
            sprintf(['tactor ' num2str(tactorLocsTempSamps)])

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
            
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
                        
            % do samples too 
            buttonPksTempVecSamps(j) = buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
            
            tactorPksVec(j) = tactorPksTemp;
            tactorLocsVec(j) = tactorLocsTemp;
                       
            tactorPksVecSamps(j) = tactorPksTempSamps;
            tactorLocsVecSamps(j) = tactorLocsTempSamps;
        end
        
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        buttonPksSamps{i} = buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
    end
end

%%
% calculate differences - MAKE SURE YOU ONLY DO THIS ONCE

buttonTactDiff = buttonLocs{uniqueCond==-1} - tactorLocsVec;
buttonTactDiffSamps = buttonLocsSamps{uniqueCond==-1} - tactorLocsVecSamps;


buttonLocs{uniqueCond==-1} = buttonTactDiff;
buttonLocsSamps{uniqueCond==-1} = buttonTactDiffSamps;

%% save it

current_direc = pwd;

save(fullfile(current_direc, [sid '_compareResponse_block_' block '.mat']),'buttonTactDiffSamps','buttonLocsSamps','s','block','sid','buttonLocs','tactorLocsVec','t_epoch','stimTimes','fs_stim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');

clearvars -except buttonTactDiffSamps buttonLocSamps s buttonLocs block t_epoch stimTimes fs_stim epochedButton tactorLocsVec epochedTactor condType uniqueCond respLo respHi SIDS DATA_DIR sid

close all

