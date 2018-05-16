%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming';
sid = SIDS{5};

load(fullfile(DATA_DIR,'ResponseTiming-1.mat'))
block = '1';

plotIt = 1;
%% load in data of interest

[stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);

clear Stim Tact Sing

%%

% figure out stim times
% vector of condition type - for first subject, looks like condition type
% is what was used , rather than test_condition, BOTH blocks used the same
% file
condType= dlmread('C:\Users\djcald.CSENETID\SharedCode\StimulationResponseTimingAnalysis\2fd831\rxnTime_condition_1_modified_9_23_2016.txt');
train = dlmread('C:\Users\djcald.CSENETID\SharedCode\StimulationResponseTimingAnalysis\2fd831\rxnTime_stimTrainDelivery_1.txt');

% for this subject, on the 1st/2nd block, seems to only be 139 trials
% for the 1st block, this dropped a no stim
condType = condType(1:139);

[trainTimesTotal,stimFromFile,trainTimes,condType,uniqueCond] = extract_stimulation_times(tact,condType);

%% extract stimulus data, find delay, and get timing of stimuli

[bursts,delay] = extract_stimulus_delivery(stim,sing,condType,trainTimes,fsStim,fsSing,plotIt);

%% extract data
% try and account for delay for the stim times
stimTimes = bursts(2,:)+delay;
trainTimes=stimTimes;

%% look at all simultaneously

tactorData = tact(:,1);
buttonData = tact(:,2);

analyze_all_inputs_simultaneously(tactorData,buttonData,stim,stimFromFile,fsTact)

%% 
respLo = 0.150;
respHi = 1;

%% quantifying data

[buttonLocs,buttonLocsSamps,tactorLocsVec,tactorLocsVecSamps,tEpoch] = get_response_timing_segs(tactorData,uniqueCond,stim,buttonData,stimFromFile,fsStim,fsTact,trainTimesTotal,plotIt);
%
%% save it
saveIt = 0;
if saveIt
    current_direc = pwd;
    
    save(fullfile(current_direc, [sid '_compareResponse_block_' block '.mat']),'buttonTactDiffSamps','buttonLocsSamps','block','sid','buttonLocs','tactorLocsVec','tEpoch','stimTimes','fsStim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');
    
    clearvars -except buttonTactDiffSamps buttonLocSamps s buttonLocs block tEpoch stimTimes fsStim epochedButton tactorLocsVec epochedTactor condType uniqueCond respLo respHi SIDS DATA_DIR sid
    
    close all
    
end