%% look at the neural data in response to cortical stimulation  

% load in subject

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

%% neural data 

raw_eco = Wave.data;
eco_fs = Wave.info.SamplingRateHz;


%% get train times

% vector of condition type - for first subject, looks like condition type
% is what was used , rather than test_condition

condType= dlmread('condition.txt');
condTestType = dlmread('test_condition.txt');
train = dlmread('testTrain.txt');


%% cortical brain data 
if s == 1  
    
end 

%% tactor brain data