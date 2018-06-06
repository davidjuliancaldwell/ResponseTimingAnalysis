%% look at the neural data in response to cortical stimulation

% this data is DC coupled
%%
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks
addpath('./scripts')

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
sid = SIDS{2};

sVec = [1,2];

for s = sVec
    
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
    
    %% neural data
    clearvars -except ECO1 ECO2 ECO3 Tact sid block s DATA_DIR s sVec folder_data epochedCortEco_cell
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
    
    % DC coupled - get rid of end
    if s ==1
        data = data(1:5339085,:);
    elseif s == 2
        data = data(1:5389295,:);
    end
    
    % only 80 channels
    
    data = data(:,1:80);
    
   
    % subtract mean? Or line fit?
    order_poly = 10;
    data = polyfitSubtract(data,order_poly);
    
    plot(data(:,10));
    hold on
    
    clearvars data_int f_y p s mu
    %%
    load([sid,'_compareResponse_block_tactorSub',block,'.mat'])
    
    %% get train times
    
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

    for i = 1:length(uniqueCond)
        
        condIntAns = uniqueCond(i);
        %% ARTIFACT
        
        if (condIntAns == 100 || condIntAns == 200 || condIntAns == 400 || condIntAns == 800)
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*ecoFs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*ecoFs);
            
            epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{i})-sampsPreStim,(trainTimesCellThresh{i}+ sampsPostStim)));
            response = buttonLocsThresh{i};
        end
        
        if condIntAns == -1
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*ecoFs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*ecoFs);
            
            %response = buttonLocsThresh{condInt} + tactorLocsVec;
            response = buttonLocsThresh{i};
            responseSamps = round(tactorLocsVec*ecoFs);
            
            adjustTact = 1;
            if adjustTact  == 1
                responseSamps = responseSamps - (ecoFs*9/1e3);
            end
            
            epochedCortEco = squeeze(getEpochSignal(data,((trainTimesCellThresh{i}+responseSamps)-sampsPreStim),((trainTimesCellThresh{i}+responseSamps)+ sampsPostStim)));
            
            
        end

        epochedCortEco_cell{s}{i} = epochedCortEco;  
    end
    
    tEpoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;
end

current_direc = pwd;

save(fullfile(current_direc, [sid 'pooledData_tactorSub_10ms.mat']),'-v7.3','epochedCortEco_cell','fsData','tEpoch');

return
