%% look at the neural data in response to cortical stimulation


Z_ConstantsStimResponse;
% add path for scripts to work with data tanks
addpath('./scripts')
%addpath('./scripts/JennysConversionScripts')

% subject directory, change as needed
% SUB_DIR = fullfile(myGetenv('subject_dir')); - for David's PC right now

% data directory

%PUT PATH TO DATA DIRECTORY WITH CONVERTED DATA FILES

% DJC Desktop
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
sid = SIDS{3};

s_vec = [1,2];

for s = s_vec
    
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
    
    % On-target stim: 20, 29 active
    % Off-target stim: 57 active -58  ground
    
    
    %% neural data
    
    % the eco data is crashing it right now
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
        
        if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
            
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

%save(fullfile(current_direc, [sid 'pooledData.mat']),'-v7.3','epochedCortEco_cell','fs_data','t_epoch');

return

