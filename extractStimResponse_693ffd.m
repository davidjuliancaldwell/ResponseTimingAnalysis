%% load in subject

% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
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

plotIt = 1;


%% load in data of interest

[stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);

clear Stim Tact Sing

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

%save(fullfile(current_direc, [sid '_compareResponse_block_' block '.mat']),'buttonTactDiffSamps','buttonLocsSamps','s','block','sid','buttonLocs','tactorLocsVec','tEpoch','stimTimes','fsStim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');

%clearvars -except buttonTactDiffSamps buttonLocSamps s buttonLocs block tEpoch stimTimes fsStim epochedButton tactorLocsVec epochedTactor condType uniqueCond respLo respHi SIDS DATA_DIR sid

%close all


