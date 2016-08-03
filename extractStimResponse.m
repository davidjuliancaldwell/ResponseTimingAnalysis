%% starting with subject a

%% initialize output and meta dir
% clear workspace
close all; clear all; clc

% set input output working directories - for David's PC right now 
% Z_ConstantsStimResponse;

% add path for scripts to work with data tanks
addpath('./scripts')

% subject directory, change as needed
% SUB_DIR = fullfile(myGetenv('subject_dir')); - for David's PC right now 

% data directory

%PUT PATH TO DATA DIRECTORY WITH CONVERTED DATA FILES
DATA_DIR = 'C:\Users\djcald\Data\ConvertedTDTfiles';

SIDS = {'acabb1'};

%% load in subject

% this is from my z_constants

sid = SIDS{1};

% ui box for input
list_str = {'sensory stimulation','tactor stimulation','off target stimulation'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, 'acabb1'))
    folder_data = strcat(DATA_DIR,'\acabb1');
    
    if s == 1
        load(fullfile(folder_data,'ReactionTime-1.mat'))
    elseif s == 2
        load(fullfile(folder_data,'ReactionTime-3.mat'))
        
    elseif s == 3
        load(fullfile(folder_data,'ReactionTime-4.mat'))
        
    end
    
end

%% load in data of interest 

stim = Stim.data;
fs_stim = Stim.info.SamplingRateHz;
data = Wave.data;
fs_data = Wave.info.SamplingRateHz;
sing = Sing.data;
fs_sing = Sing.info.SamplingRateHz;

tact = Tact.data;
fs_tact = Tact.info.SamplingRateHz;

valu = Valu.data;
fs_valu = Valu.info.SamplingRateHz;

%% plot stim
%
figure
hold on
for i = 1:size(stim,2)

    t = (0:length(stim)-1)/fs_stim;
    subplot(2,2,i)
    plot(t*1e3,stim(:,i))
    title(sprintf('Channel %d',i))


end


xlabel('Time (ms)')
ylabel('Amplitude (V)')

subtitle('Stimulation Channels')



%% Sing looks like the wave to be delivered, with amplitude in uA
% Try working from this

% build a burst table with the timing of stimuli
bursts = [];

Sing1 = sing(:,1);

samplesOfPulse = round(2*fs_stim/1e3);



% trying something like A_BuildStimTables from BetaStim


Sing1Mask = Sing1~=0;
dmode = diff([0 Sing1Mask' 0 ]);


dmode(end-1) = dmode(end);


bursts(2,:) = find(dmode==1);
bursts(3,:) = find(dmode==-1);

stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stims,1)-1)/fs_sing;
t = t*1e3;
figure
plot(t,stims,'b','linewidth',2)
xlabel('Time (ms)');
ylabel('Current to be delivered (mA)')
ylim([(min(stims(:))-100) (max(stims(:))+100)])
title('Current to be delivered for all trials')

% delay loks to be 0.2867 ms from below.

%% Plot stims with info from above

stim1 = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,:)-1),(bursts(3,:))+1));
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim1Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Finding the delay between current output and stim delivery')

% hold on

plot(t,stims)

% get the delay in stim times

delay = round(0.2867*fs_stim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,:)-1+delay;
stimTimesEnd = bursts(3,:)-1+delay;
stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd+5));
t = (0:size(stim1Epoched,1)-1)/fs_stim;
t = t*1e3;
figure
plot(t,stim1Epoched)
xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Stim voltage monitoring with delay added in')



%% extract data

% try and account for delay for the stim times
stimTimes = bursts(2,:)-1+delay;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * fs_data); % pre time in sec
postsamps = round(1 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after


% sampling rate conversion between stim and data
fac = fs_stim/fs_data;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);

%% look at tactor data 
% 
% figure
% hold on
% for i = 1:size(tact,2)
% 
%     t = (0:length(tact)-1)/fs_tact;
%     subplot(2,2,i)
%     plot(t*1e3,tact(:,i))
%     title(sprintf('Channel %d',i))
% 
% 
% end
% 
% 
% xlabel('Time (ms)')
% ylabel('Amplitude (V)')
% 
% subtitle('tact')

%% look at tactor 

tactorData = tact(:,1);
t_tact = (0:length(tactorData)-1)/fs_tact;
figure
plot(t_tact,tactorData);

title('tactor data')



%% look at button press 

buttonData = tact(:,2);
t_button = (0:length(buttonData)-1)/fs_tact;
figure
plot(t_button,buttonData);

title('button data')

% make masked button
% threshold of 0.009
buttonMask = buttonData > 0.009;

%% look at stim from file saved

stimFromFile = tact(:,3);
t_stimFile = (0:length(stim)-1)/fs_tact;
figure
plot(t_stimFile,stimFromFile);
title('stim from file')

%% look at all 3 + stim waveform 

figure
ax1 = subplot(4,1,1);
plot(t_tact,tactorData)
title('tactor data')

ax2 = subplot(4,1,2);
plot(t_button,buttonData);
title('button data')

ax3 = subplot(4,1,3);
plot(t_stimFile,stimFromFile);
title('stim from file')


% assuming stim1 here is the channel where stim was being delivered
ax4 = subplot(4,1,4);
plot(t_stimFile,stim1);

%link axis
linkaxes([ax1,ax2,ax3,ax4],'x')

%% look at valu
figure
hold on
for i = 1:size(valu,2)

    t = (0:length(valu)-1)/fs_stim;
    subplot(5,1,i)
    plot(t*1e3,valu(:,i))
    title(sprintf('Channel %d',i))


end


xlabel('Time (ms)')
ylabel('Amplitude (V)')

subtitle('valu')


