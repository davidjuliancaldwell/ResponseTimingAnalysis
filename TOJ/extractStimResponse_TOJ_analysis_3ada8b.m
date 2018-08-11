% this is from my z_constants
close all;clear all;clc
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\3ada8b\data\d10\MATLAB_conversions\3ada8b_TOJ';
sid = SIDS{6};

folder_data = strcat(DATA_DIR);
load(fullfile(folder_data,'TOJ-1.mat'))
load(fullfile(folder_data,'TOJ-1_TOJ.mat'))
%%
% load in data of interest
stim = Stim.data;

fsStim = Stim.info.SamplingRateHz;

clear Stim

clear ECO1 ECO2 ECO3

tact = Tact.data;

fsTact = Tact.info.SamplingRateHz;
clear Tact

plotIt = 1;
%%
% get that button press - threshold it
tact(tact(:,2) >= 0.009,2) = 0.009;

% check to makes ure it all is working
%[buttonPksTemp,buttonLocsTemp] = findpeaks(tact(:,2),fsTact,'minpeakdistance',2,'Minpeakheight',0.008);

tact(tact(:,2) < 0.009,2) = 0;
tact(:,2) = tact(:,2)*1000;

%plot  inputs continuously through time
if plotIt
    plot_TOJ_concurrently(stim,tact,fsStim)
end
%%
% QUANTIFY RXN TIME TO CORTICAL STIM
[stimTimes,trainTimes] = extract_stimulation_times_TOJ_readIn_v2(tact,fsStim,[]);
%%
[epochedTactor,epochedAudio,epochedStim,epochedButton,t,tSamps] = extract_epochs_TOJ(stim,tact,trainTimes,fsStim);
% epoched button press

%%
% map 0 to tactor first, map 1 to stim first
stimFirst = ismember(feltFirstVec,'s');
tactorFirst = ismember(feltFirstVec,'t');
whichPerceived = zeros(size(feltFirstVec));
whichPerceived(stimFirst) = 1;
whichPerceived(tactorFirst) = 0;
%%
numTrials = size(epochedAudio,2);

if plotIt
    plot_epochs_TOJ(epochedStim,epochedTactor,epochedAudio,epochedButton,whichPerceived,t,numTrials)
end


%%

[tactorLocsVec,stimLocsVec,buttonLocsVec,tactorStimDiff,responseTimes] = get_response_timing_segs_TOJ(epochedButton,epochedTactor,epochedStim,t,tSamps,numTrials);


%%
saveIt = 0;

if saveIt
    current_direc = pwd;
    
    save(fullfile(current_direc, [sid '_TOJ_matlab.mat']),'tactorStimDiff','responseTimes','t','trainTimes',...,
        'fsStim','epochedButton','epochedTactor','epochedAudio','epochedStim','whichPerceived');
end
