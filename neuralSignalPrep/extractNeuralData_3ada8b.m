%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\3ada8b\data\d9\MATLAB_conversions\3ada8b_ResponseTiming';
sid = SIDS{6};

% THIS SCRIPT MULTIPLIES BY 4 FOR THE ECOG SIGNAL


for s = 1:2
    % load in data
    if (strcmp(sid, '3ada8b'))
        
        if s == 1
            load(fullfile(DATA_DIR,'responseTiming-1.mat'))
            block = '1';
        elseif s == 2
            load(fullfile(DATA_DIR,'responseTiming-2.mat'))
            block = '2';
        end
        
    end
    
    plotIt = 1;
    
    %% neural data
    
    % the eco data is crashing it right now
    clearvars -except ECO1 ECO2 ECO3 Tact block s DATA_DIR s sVec folderData sid
    eco1 = ECO1.data;
    fsData = ECO1.info.SamplingRateHz;
    ecoFs = fsData;
    clear ECO1
    eco2 = ECO2.data;
    clear ECO2
    
    eco3 = ECO3.data;
    clear ECO3
    
    
    ECoG = 4*[eco1 eco2 eco3];
    clearvars eco1 eco2 eco3
    
    % only 64 channels grid
    
    % get rid of bad channels
    bads = [79:128];
    goodVec = logical(ones(size(ECoG,2),1));
    goodVec(bads) = 0;
    ECoG = ECoG(:,goodVec);
    
    %
    load([sid,'_compareResponse_block_',block,'.mat'])
    
    % get train times
    
    % look at stim from file saved (this is the sample where things were
    % delivered
    tact = Tact.data;
    fsTact = Tact.info.SamplingRateHz;
    stimFromFile = tact(:,3);
    
    % get stimulation times of delivery
    trainTimes = find(stimFromFile~=0);
    
    % convert sample times for eco
    
    convertSamps = fsTact/fsData;
    trainTimesConvert = round(stimTimes/convertSamps);
    stimChans = [3 4 24 32];
    
    convertSamps = fsTact/fsData;
    
    trainTimesConvert = round(stimTimes/convertSamps);
    
    trainTimesCell = {};
    trainTimesCellThresh = {};
    
    for i = 1:length(uniqueCond)
        
        trainTimesCell{i} = trainTimesConvert(condType==uniqueCond(i));
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
    
    %% ARTIFACT
    epochedCortEco_cell = {};
    for i = 1:length(uniqueCond)
        
        condIntAns = uniqueCond(i);
        condIntAns
        % ARTIFACT
        
        if (condIntAns >=0)
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*ecoFs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*ecoFs);
            
            epochedCortEco = squeeze(getEpochSignal(ECoG,(trainTimesCellThresh{i})-sampsPreStim,(trainTimesCellThresh{i}+ sampsPostStim)));
            response = buttonLocsThresh{i};
        elseif condIntAns == -1
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*ecoFs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*ecoFs);
            
            %response = buttonLocsThresh{condInt} + tactorLocsVec;
            response = buttonLocsThresh{i};
            responseSamps = round(tactorLocsVec*ecoFs);
            
            % 11-3-2017 - account for nan's in response_samps vector
            response_mask = (~isnan(responseSamps));
            
            adjustTact = 0;
            if adjustTact  == 1
                responseSamps = responseSamps - (ecoFs*9/1e3);
            end
            
            
            epochedCortEco = squeeze(getEpochSignal(ECoG,(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)-sampsPreStim),(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)+ sampsPostStim)));
            tactData = decimate(Tact.data,2)';
            epochedTactorNew = squeeze(getEpochSignal(tactData,(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)-sampsPreStim),(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)+ sampsPostStim)));
        end
        
        epochedCortEco_cell{i} = epochedCortEco;
    end
    
    tEpoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;
    current_direc = pwd;
    
    save(fullfile(current_direc, [sid '_priming_neural_block_' num2str(s) '.mat']),'-v7.3','epochedCortEco_cell','fsData','tEpoch');
end
