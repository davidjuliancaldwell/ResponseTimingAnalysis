%% look at the neural data in response to cortical stimulation

% this data is DC coupled
%%
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks
addpath('./scripts')

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
sid = SIDS{2};

s_vec = [1,2];

for s = s_vec
    
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
    clearvars -except ECO1 ECO2 ECO3 Tact sid block s DATA_DIR s s_vec folder_data epochedCortEco_cell
    eco1 = ECO1.data;
    fs_data = ECO1.info.SamplingRateHz;
    eco_fs = fs_data;
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
    load([sid,'_compareResponse_block_',block,'.mat'])
    
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
        
        if (condIntAns == 100 || condIntAns == 200 || condIntAns == 400 || condIntAns == 800)
            
            post_stim = 2000;
            samps_post_stim = round(post_stim/1e3*eco_fs);
            
            pre_stim = 1000;
            samps_pre_stim = round(pre_stim/1e3*eco_fs);
            
            epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{i})-samps_pre_stim,(trainTimesCellThresh{i}+ samps_post_stim)));
            response = buttonLocsThresh{i};
        end
        
        if condIntAns == -1
            
            post_stim = 2000;
            samps_post_stim = round(post_stim/1e3*eco_fs);
            
            pre_stim = 1000;
            samps_pre_stim = round(pre_stim/1e3*eco_fs);
            
            %response = buttonLocsThresh{condInt} + tactorLocsVec;
            response = buttonLocsThresh{i};
            response_samps = round(tactorLocsVec*eco_fs);
            epochedCortEco = squeeze(getEpochSignal(data,((trainTimesCellThresh{i}+response_samps)-samps_pre_stim),((trainTimesCellThresh{i}+response_samps)+ samps_post_stim)));
            
            
        end

        epochedCortEco_cell{s}{i} = epochedCortEco;  
    end
    
    t_epoch = (-samps_pre_stim:samps_post_stim-1)/eco_fs;
end

current_direc = pwd;

save(fullfile(current_direc, [sid 'pooledData.mat']),'-v7.3','epochedCortEco_cell','fs_data','t_epoch');

return
