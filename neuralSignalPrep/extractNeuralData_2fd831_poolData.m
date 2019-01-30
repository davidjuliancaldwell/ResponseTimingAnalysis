%% look at the neural data in response to cortical stimulation

%%
Z_ConstantsStimResponse;

dataDirTotal = fullfile(DATA_DIR,'ConvertedTDTfiles');

sid = SIDS{4};
s_vec = [1,2];

for s = s_vec
    % load in data
    if (strcmp(sid, '2fd831'))
        folder_data = strcat(dataDirTotal,'\2fd831');
        
        if s == 1
            load(fullfile(folder_data,'ResponseTiming-1.mat'))
            block = '1';
        elseif s == 2
            load(fullfile(folder_data,'ResponseTiming-2.mat'))
            block = '2';
        end
    end
    
    %% neural data
    
    % the eco data is crashing it right now
    clearvars -except ECO1 ECO2 ECO3 Tact sid block s dataDirTotal DATA_DIR s s_vec folder_data epochedCortEco_cell
    eco1 = ECO1.data;
    fs_data = ECO1.info.SamplingRateHz;
    fsData = fs_data;
    
    eco_fs = fs_data;
    ecosFs = eco_fs;
    clear ECO1
    eco2 = ECO2.data;
    clear ECO2
    
    eco3 = ECO3.data;
    clear ECO3
    
    data = [eco1 eco2 eco3];
    clearvars eco1 eco2 eco3
    
    % git rid of the bit where it jumps down at the end
    %%
    if s==1
        data = data(1:5919739,:);
    elseif s == 2
        data = data(1:6989821,:);
    end
    
    % only 64 channels grid
    
    data = data(:,1:64);
    %%
    
    % figure
    % for i = 1:size(data,2)
    %     plot(data(:,i))
    %         title(num2str(i))
    %
    %     pause(1)
    % end
    
    %%
    % subtract mean? Or line fit?
    order_poly = 10;
    data = polyfitSubtract(data,order_poly);
    
    %%
    clearvars data_int f_y p s mu
    %%
    %   load([sid,'_compareResponse_block_tactorSub',block,'.mat'])
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
    %%
    for i = 1:length(uniqueCond)
        
        condIntAns = uniqueCond(i);
        %% ARTIFACT
        
        if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*eco_fs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*eco_fs);
            
            epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{i})-sampsPreStim,(trainTimesCellThresh{i}+ sampsPostStim)));
            response = buttonLocsThresh{i};
        end
        
        if condIntAns == -1
            
            postStim = 2000;
            sampsPostStim = round(postStim/1e3*eco_fs);
            
            preStim = 1000;
            sampsPreStim = round(preStim/1e3*eco_fs);
            
            %response = buttonLocsThresh{condInt} + tactorLocsVec;
            response = buttonLocsThresh{i};
            responseSamps = round(tactorLocsVec*eco_fs);
            
            adjustTact = 0;
            if adjustTact  == 1
                responseSamps = responseSamps - (ecoFs*9/1e3);
            end         
            
            % 11-3-2017 - account for nan's in response_samps vector
            response_mask = (~isnan(responseSamps));
            
            epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)-sampsPreStim),(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)+ sampsPostStim)));
            
        end
        
        epochedCortEco_cell{s}{i} = epochedCortEco;
    end
    t_epoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;
    
end

current_direc = pwd;

%save(fullfile(current_direc, [sid 'pooledData_tactorSub_10ms.mat']),'-v7.3','epochedCortEco_cell','fs_data','t_epoch');
save(fullfile(current_direc, [sid 'pooledData_changePts_noDelay.mat']),'-v7.3','epochedCortEco_cell','fsData','tEpoch');

return

