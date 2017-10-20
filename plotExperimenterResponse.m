%% 7-11-2017 - script to plot expirementer's response times for Kurt

close all;clear all;clc

Z_ConstantsStimResponse;

% add path for scripts to work with data tanks
addpath('./scripts')
%addpath('./scripts/JennysConversionScripts')

% subject directory, change as needed
% SUB_DIR = fullfile(myGetenv('subject_dir')); - for David's PC right now

% data directory

%PUT PATH TO DATA DIRECTORY WITH CONVERTED DATA FILES

% DJC Desktop
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';

%%
% first subject - this is where we had to touch on the beep, really a test
% of my reaction time

sid = SIDS{1};
load([sid,'_compareResponse.mat']);
tactorTotal = [tactorLocsVecTact(tactorLocsVecTact>0.1 & tactorLocsVecTact<0.8)];
nbins = [15];
histogram(1e3.*tactorTotal,nbins)
title(["Distribution  of Experimenter's Reaction Times"])
xlabel('Time (ms)')
ylabel('Count')
xlim([0 1000])
set(gca,'FontSize',14)

%% look at jitter around the 1 s trigger for the tactor

tactorTotal = [];

% iterate through the rest
for i = 2:length(SIDS)
    sid = SIDS{i};
    for block = 1:2
        load([sid,'_compareResponse_block_',num2str(block),'.mat'])
        tactorTotal = [tactorTotal tactorLocsVec];
    end
end

% plot it
nbins = 15;
histogram(1e3.*tactorTotal-1000,nbins)
title({"Distribution of experimenter's reaction time","centered around desired stimulus onset"})
xlabel('Time (ms)')
ylabel('Count')
set(gca,'FontSize',14)