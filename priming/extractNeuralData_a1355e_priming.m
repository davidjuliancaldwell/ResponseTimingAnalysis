%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming';
sid = SIDS{5};

sVec = [1,2];

for s = sVec
    % load in data
    if (strcmp(sid, 'a1355e'))
        folderData = strcat(DATA_DIR,'\2fd831');
        
        if s == 1
            load(fullfile(DATA_DIR,'ResponseTiming-3.mat'))
            block = '3';
        elseif s == 2
            load(fullfile(DATA_DIR,'ResponseTiming-4.mat'))
            block = '4';
        end
        
    end
    
    
    %% neural data
    
    % the eco data is crashing it right now
    clearvars -except ECO1 ECO2 ECO3 Tact sid block s DATA_DIR s sVec folderData 
    eco1 = ECO1.data;
    fsData = ECO1.info.SamplingRateHz;
    ecoFs = fsData;
    clear ECO1
    eco2 = ECO2.data;
    clear ECO2
    
    eco3 = ECO3.data;
    clear ECO3
    
    data = [eco1 eco2 eco3];
    clearvars eco1 eco2 eco3
    
    % only 64 channels grid
    
    % data = data(:,1:64);
    
    %
    load([sid,'_Priming_block_',block,'.mat'])
    
    % get train times
    
    % look at stim from file saved (this is the sample where things were
    % delivered
    tact = Tact.data;
    fs_tact = Tact.info.SamplingRateHz;
    stimFromFile = tact(:,3);
    
    % get stimulation times of delivery
    trainTimes = find(stimFromFile~=0);
    
    % convert sample times for eco
    
    convertSamps = fs_tact/fsData;
    
    trainTimesConvert = round(stimTimes/convertSamps);
    
    trainTimesCell = {};
    trainTimesCellThresh = {};
    
    for i = 1:length(uniqueCond)
        
        trainTimesCell{i} = trainTimesConvert(primedOption==uniqueCond(i));
        trim = buttonLocs{i};
        trim = trim(trim>respLo & trim<respHi);
        zTrim = zscore(trim);
        if ~isempty(trainTimesCell{i}) % check to make sure not indexing empty cell
            %trainTimesCellThresh{i} = trainTimesCell{i}(abs(zTrim)<3); % z score
            % buttonLocsThresh = buttLocs{i}(abs(zTrim)<3);
            
            trainTimesCellThresh{i} = trainTimesCell{i};% no zscore
            buttonLocsThresh{i} = buttonLocs{i};% no zscore
            
        end
    end
    
    % ARTIFACT
    
    for i = 1:length(uniqueCond)
               condIntAns = uniqueCond(i);
        postStim = 2000;
        sampsPostStim = round(postStim/1e3*ecoFs);
        
        preStim = 1000;
        sampsPreStim = round(preStim/1e3*ecoFs);
        
        epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{i})-sampsPreStim,(trainTimesCellThresh{i}+ sampsPostStim)));
        
        response = buttonLocsThresh{i};
        
        epochedCortEco_cell{i} = epochedCortEco;
        
    end
    tEpoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;
    current_direc = pwd;

    save(fullfile(current_direc, [sid '_priming_neural_block_' num2str(s) '.mat']),'-v7.3','epochedCortEco_cell','fsData','tEpoch');
    
    
end