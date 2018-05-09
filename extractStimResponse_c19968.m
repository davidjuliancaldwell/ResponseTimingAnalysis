%% load in subject

% this is from my z_constants
Z_ConstantsStimResponse;

sid = SIDS{2};
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';

% ui box for input
list_str = {'1st block','2nd block','1st block with no tactor','2nd block with no tactor'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, 'c19968'))
    folder_data = strcat(DATA_DIR,'\c19968');
    
    if s == 1
        load(fullfile(folder_data,'ReactionTime_c19968-7.mat'))
        block = '1';
    elseif s == 2
        load(fullfile(folder_data,'ReactionTime_c19968-11.mat'))
        block = '2';
    elseif s == 3
        load(fullfile(folder_data,'ReactionTime_c19968-3.mat'))
        block = '1_cort_stimOnly';
    elseif s == 4
        load(fullfile(folder_data,'ReactionTime_c19968-4.mat'))
        block = '1_cort_stimOnly';
    end
    
end

plotIt = 1;

%% load in data of interest

[stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);

clear Stim Tact Sing

%% figure out stim times
% vector of condition type - for first subject, looks like condition type
% is what was used , rather than test_condition
if s == 2
    condType= dlmread('rxnTime_condition_2.txt');
    train = dlmread('rxnTime_stimTrainDelivery_2.txt');
else
    condType= dlmread('rxnTime_condition_1.txt');
    train = dlmread('rxnTime_stimTrainDelivery_1.txt');
end

% shrink condition type to be 140
% as the experiment went S1 (110 stims), tactor (30 stims) , then break
condType = condType(1:140);

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

current_direc = pwd;

%save(fullfile(current_direc, [sid '_compareResponse_block_' block '.mat']),'buttonTactDiffSamps','buttonLocsSamps','s','block','sid','buttonLocs','tactorLocsVec','t_epoch','stimTimes','fs_stim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');

clearvars -except buttonTactDiffSamps buttonLocSamps s buttonLocs block tEpoch stimTimes fs_stim epochedButton tactorLocsVec epochedTactor condType uniqueCond respLo respHi SIDS DATA_DIR sid

%close all


