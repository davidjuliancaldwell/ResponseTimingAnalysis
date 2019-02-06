%% look at the neural data in response to cortical stimulation

Z_ConstantsStimResponse;

dataDirTotal = fullfile(DATA_DIR,'ConvertedTDTfiles'); 
sid = SIDS{3};

s_vec = [1,2];

for s = s_vec
    % load in data
    if (strcmp(sid, '693ffd'))
        folder_data = strcat(dataDirTotal,'\693ffd');
        
        if s == 1
            load(fullfile(folder_data,'ReactionTime_693ffd-2.mat'))
            block = '1';
        elseif s == 2
            load(fullfile(folder_data,'ReactionTime_693ffd-4.mat'))
            block = '2';             
        end       
    end
    % On-target stim: 20, 29 active
    % Off-target stim: 57 active -58  ground
    
    %% neural data
    
    % the eco data is crashing it right now
    clearvars -except ECO1 ECO2 ECO3 Tact sid block s DATA_DIR dataDirTotal s s_vec folder_data epochedCortEco_cell
    eco1 = ECO1.data;
    fs_data = ECO1.info.SamplingRateHz;
    fsData = fs_data;
    ecoFs = fs_data;
    clear ECO1
    eco2 = ECO2.data;
    clear ECO2
    
    eco3 = ECO3.data;

    data = [eco1 eco2 eco3];
    clearvars eco1 eco2 eco3
    
    % git rid of the bit where it jumps down at the end
    if s==1
        data = data(1:5441534,:);
    elseif s == 2
        data = data(1:5456873,:);
    end
    
    % only 64 channels, and 56 and up look terrible
    
    data = data(:,1:56);
    
    
    %%
    % subtract mean? Or line fit?
    % subtract mean? Or line fit?
    order_poly = 10;
    data = polyfitSubtract(data,order_poly);
    
    % for i = 1:size(data,2)
    %
    %     data_int = data(:,i);
    %     [p,s,mu] = polyfit((1:numel(data_int))',data_int,10);
    %     f_y = polyval(p,(1:numel(data_int))',[],mu);
    %
    %     data(:,i) = data(:,i) - f_y;
    % end
    
    plot(data(:,55));
    hold on
    
    clearvars data_int f_y p s mu
    %%
  %  load([sid,'_compareResponse_block_tactorSub',block,'.mat'])
    load([sid,'_compareResponse_block_',block,'_changePts_noDelay.mat'])
    
    %% get train times
    
    % look at stim from file saved (this is the sample where things were
    % delivered
    tact = Tact.data;
    fs_tact = Tact.info.SamplingRateHz;
    stimFromFile = tact(:,3);
    
    % get stimulation times of delivery
    trainTimes = find(stimFromFile~=0);
    
    % convert sample times for eco
    convertSamps = fs_tact/fs_data;
    
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
        
        if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
            
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
            
                   adjustTact = 0;
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

%save(fullfile(current_direc, [sid 'pooledData_tactorSub_10ms.mat']),'-v7.3','epochedCortEco_cell','fs_data','t_epoch');
save(fullfile(current_direc, [sid 'pooledData_changePts_noDelay.mat']),'-v7.3','epochedCortEco_cell','fsData','tEpoch');

return

