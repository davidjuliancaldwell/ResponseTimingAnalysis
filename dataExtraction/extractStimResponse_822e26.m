%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;

subjdir = getenv('SUBJECT_DIR');
DATA_DIR = fullfile(subjdir,'\ConvertedTDTfiles');

sid = SIDS{7};
%
% % ui box for input
% list_str = {'1st block','2nd block'};
%
% [s,v] = listdlg('PromptString','Pick experiment',...
%     'SelectionMode','single',...
%     'ListString',list_str);

for s = 2:2
    % load in data
    if (strcmp(sid, '822e26'))
        folder_data = strcat(DATA_DIR,'\',sid);
        
        if s == 1
            load(fullfile(folder_data,'responseTiming-2.mat'))
            block = '1';
        elseif s == 2
            load(fullfile(folder_data,'responseTiming-5.mat'))
            block = '2';
        end
        
    end
    plotIt = 1;
    %% load in data of interest
    
    [stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);
    
    clear Stim Tact Sing
    
    %% figure out stim times
    
    folder  = 'C:\Users\david\Data\Subjects\822e26\822e26_responseTiming';
    primedOption = dlmread(fullfile(folder,'rxnTime_primingOption_1_priming_5_cond.txt'));
    condType= dlmread(fullfile(folder,'rxnTime_condition_1_priming_5_cond.txt'));
    train = dlmread(fullfile(folder,'rxnTime_stimTrainDelivery_1_priming_5_cond.txt'));
    
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
    
    %% quantifying data
    
    [buttonLocs,buttonLocsSamps,tactorLocsVec,tactorLocsVecSamps,tEpoch,epochedButton,epochedTactor,buttonTactDiffSamps] = get_response_timing_segs_newTactor(tactorData,uniqueCond,stim,buttonData,stimFromFile,fsStim,fsTact,trainTimesTotal,plotIt);
    %% get ISI info
    
    %[ISICellSamps,ISICellSeconds,ISICondBefore,ISICellSampsNoNuOt,ISICellSecondsNoNuOt,ISIcondBeforeNoNuOt] = get_ISI(condType,uniqueCond,tactorLocsVecSamps,stimFromFile,fsStim,trainTimesTotal,trainTimes);
    %% look at RT vs ISI
    
    % [mdl,mdlNoNuOt] = compare_resp_times_ISI(uniqueCond,buttonLocs,ISICellSecondsNoNuOt,ISICellSeconds);
    %% save it
    saveIt = 1;
    
    if saveIt
        current_direc = pwd;
        
        %save(fullfile(current_direc, [sid '_compareResponse_block_tactorSub' block '.mat']),'buttonLocsSamps','block','sid','buttonLocs','tactorLocsVec','tEpoch','stimTimes','fsStim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');
        save(fullfile(current_direc, [sid '_compareResponse_block_' block '_changePts_noDelay.mat']),'buttonTactDiffSamps','buttonLocsSamps','s','block','sid','buttonLocs','tactorLocsVec','tEpoch','stimTimes','fsStim','epochedButton','epochedTactor','condType','uniqueCond', 'respLo','respHi');
        
        clearvars -except buttonTactDiffSamps buttonLocSamps s buttonLocs block tEpoch stimTimes fsStim epochedButton tactorLocsVec epochedTactor condType uniqueCond respLo respHi SIDS DATA_DIR sid
        
        close all
        
    end
    
end