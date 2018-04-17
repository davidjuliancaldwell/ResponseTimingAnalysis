%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks
addpath('./scripts')

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming';
sid = SIDS{5};

plotIt = 1;

% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, 'a1355e'))
    folder_data = strcat(DATA_DIR,'\2fd831');
    
    if s == 1
        load(fullfile(DATA_DIR,'ResponseTiming-3.mat'))
        block = '3';
    elseif s == 2
        load(fullfile(DATA_DIR,'ResponseTiming-4.mat'))
        block = '4';
    end
    
end

%% load in data of interest

stim = Stim.data;
fs_stim = Stim.info.SamplingRateHz;

clear Stim

fs_data = ECO1.info.SamplingRateHz;

clear ECO1 ECO2 ECO3

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
% is what was used , rather than test_condition, BOTH blocks used the same
% file

condType = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_condition_primingPilot.txt');
primedOption = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_primedOption_primingPilot.txt');
train = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_stimTrainDelivery_primingPilot.txt');
% stim cue from file
stimFromFile = tact(:,3);

%
trainTimes = find(stimFromFile~=0);

%%
% shrink condition type to be 140
% as the experiment went S1 (110 stims), tactor (30 stims) , then break

% for this subject, on the 1st/2nd block, seems to only be 139 trials
% for the 1st block, this dropped a no stim


% pick condition type where stimulation was delivered

uniqueCond = unique(primedOption);
trainTimesTotal = {};

for i = 1:length(uniqueCond)
    trainTimesTotal{i} = trainTimes(primedOption==uniqueCond(i));
end

% just in case - keep these around
trainTimesCondNoPrime = trainTimes(primedOption==0);
trainTimesCondPrime = trainTimes(primedOption==1);


%% plot stim
%
if plotIt
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
    
end
stim1 = stim(:,1);

%% Sing looks like the wave to be delivered, with amplitude in uA


% 1st stim channel
Sing1 = sing(:,1);
% 2nd stim channel
Sing2 = sing(:,4);

samplesOfPulse = round(2*fs_stim/1e3);

% build a burst table with the timing of stimuli
bursts = [];
bursts(1,:) = primedOption;
bursts(2,:) = trainTimes;
bursts(3,:) = trainTimes + samplesOfPulse;

% delay loks to be 0.2867 ms from below.

% get the delay in stim times

delay = round(0.2867*fs_stim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType==1)-1+delay;
stimTimesEnd = bursts(3,condType==1)-1+delay;

% extract data
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

ax1 = subplot(3,1,1);
plot(t_button,buttonData);
title('button data')

ax2 = subplot(3,1,2);
plot(t_stimFile,stimFromFile);
title('stim from file')

% assuming stim1 here is the channel where stim was being delivered
ax3 = subplot(3,1,3);
plot(t_stimFile,stim1);
title('S1 Stim Channel')

%link axis
linkaxes([ax1,ax2,ax3],'x')

respLo = 0.100;
respHi = 1;

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

legend({'Button Data','Button Press Onset Peaks','Stimulation Times From File','S1 stim output'})

% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)
%% QUANTIFY RXN TIME TO CORTICAL STIM

sampsEnd = round(3.5*fs_stim);

% epoched button press

epochedButton = {};

% the last trial of the epoched button press for the null condition is
% clipped - so omit that one

for i = 1:length(uniqueCond)
    epochedButton{i} = squeeze(getEpochSignal(buttonDataClip,trainTimesTotal{i},(trainTimesTotal{i} + sampsEnd)));
end

t_epoch = [0:size(epochedButton{1},1)-1]/fs_stim;
t_epoch_samps = [0:size(epochedButton{1},1)-1];


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
    
    
end

%%
% 9-13-2016 - script to compare response times once they've been calculated

respLo = 0.100;
respHi = 1;

for i = 1:length(uniqueCond)
    
    trim = buttonLocs{i};
    trim = trim(trim>respLo & trim<respHi);
    zTrim = zscore(trim);
    buttonLocsThresh{i} = 1e3.*trim(abs(zTrim)<3);
    
end


%% Histogram

% set number of bins
nbins = 15;

% make cell array for legends
uniqueCondText = cellstr(num2str(uniqueCond));
uniqueCondText{1} = 'no Prime';
uniqueCondText{2} = 'Prime';

%individual  histogram of each condition type
figure

for i = 1:length(uniqueCond)
    subplot(length(uniqueCond),1,i)
    a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
    title(['Histogram of reaction times for condition ' uniqueCondText{i} ])
    xlabel('Time (ms)')
    ylabel('Count')
    xlim([0 1000])
    a.BinWidth = 10;
    a.Normalization = 'probability';
end
subtitle([' block ', block])

% overall histogram

colormap lines;
cmap = colormap ;

figure
hold on
leg = {};
for i = 1:length(uniqueCond)
    a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
    leg{end+1} = uniqueCondText{i};
    a.FaceColor = cmap(i,:);
    a.BinWidth = 10;
    a.Normalization = 'probability';
    xlim([0 800])
    
end
xlabel('Time (ms)')
ylabel('Count')
title(['Histogram of response Times for block ', block])

legend(leg)

% histogram just for tactor and 100,200,400,800


%BOX PLOT

% change colormap to matlab default lines

combinedInfo = [];
groups = [];
colors = [];
leg = {};
temp = cmap(1,:);
cmap(1,:) = cmap(2,:);
cmap(2,:) = temp;

keeps = [1 2];

j = length(keeps);
k = 1;

for i = keeps
    combinedInfo = [buttonLocsThresh{i} combinedInfo ];
    groups = [(j).*ones(length(buttonLocsThresh{i}),1)' groups];
    colors(i,:) = cmap(i,:);
    leg{end+1} = uniqueCondText{i};
    j = j - 1;
    k = k + 1;
end

figure
prettybox(combinedInfo,groups,colors,1,true)
fig1 = gca;

ylim([0 1000])
fig1.XTick = [];
legend(findobj(gca,'Tag','Box'),leg)
ylabel('Response times (ms)')
title(['Reaction Times for block ',block])

%% statistics! kruskal wallis to start

groupsKW = [];
for i = keeps
    groupsKW = [cellstr(repmat(uniqueCondText{i},[length(buttonLocsThresh{i}),1])); groupsKW(:)];
end

[p,table,stats] = kruskalwallis(combinedInfo,groupsKW);
[c,m,h,nms] = multcompare(stats);

[nms num2cell(m)]

c((c(:,6)<0.05),[1 2 6])


%% save it

current_direc = pwd;

save(fullfile(current_direc, [sid '_Priming_block_' block '.mat']),'buttonLocsSamps',...
    's','block','sid','primedOption','buttonLocs','t_epoch','stimTimes','fs_stim','epochedButton',...
    'uniqueCond', 'respLo','respHi');


