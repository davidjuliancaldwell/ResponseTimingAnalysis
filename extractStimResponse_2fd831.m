%% load in subject

% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
sid = SIDS{4};

% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, '2fd831'))
    folder_data = strcat(DATA_DIR,'\2fd831');
    
    if s == 1
        load(fullfile(folder_data,'ResponseTiming-1.mat'))
        block = '1';
    elseif s == 2
        load(fullfile(folder_data,'ResponseTiming-2.mat'))
        block = '2';
    end
    
end

plotIt = 1;

%% load in data of interest

[stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);

clear Stim Tact Sing


%% figure out stim times
% vector of condition type
% BOTH blocks used the same file
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

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * fsData); % pre time in sec
postsamps = round(1 * fsData); % post time in sec, % modified DJC to look at up to 300 ms after

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


