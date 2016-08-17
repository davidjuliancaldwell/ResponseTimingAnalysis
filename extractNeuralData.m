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

% look at stim from file saved (this is the sample where things were
% delivered
tact = Tact.data;
fs_stim = Stim.info.SamplingRateHz;
fs_tact = Tact.info.SamplingRateHz;
stimFromFile = tact(:,3);

% get stimulation times of delivery
trainTimes = find(stimFromFile~=0);

% shrink condition type to be 120
condType = condType(1:120);

% convert sample times for eco

convertSamps = fs_stim/eco_fs;

trainTimes = round(trainTimes/convertSamps);

% get values for z-scoring, in order to filter later if desired

zTact = zscore(tactorLocsVecTactTrim);
zDiff = zscore(buttonTactDiffTrim);
zCort = zscore(buttonLocsVecCortTrim);

tactor = 1e3.*tactorLocsVecTactTrim(abs(zTact)<3);
difference = 1e3.*buttonTactDiffTrim(abs(zDiff)<3);
cort = 1e3.*buttonLocsVecCortTrim(abs(zCort)<3);


%% cortical brain data
if s == 1
    
    % pick condition type where stimulation was delivered
    if s == 1
        trainTimesCond1 = trainTimes(condType==0);
    elseif s == 2
        trainTimesCond1 = trainTimes(condType==0 | condType==1);
    end
    
    % where to begin plotting with artifact
    artifact_end = round(0.21*eco_fs);
    
    %where to end plotting
    sampsEnd = round(2*eco_fs);
    
    % only pick train times which met certain threshold
    trainTimesScreen = trainTimesCond1(buttonLocsVecCort>respLo & buttonLocsVecCort<respHi);
    trainTimesScreen = trainTimesScreen(abs(zCort)< 3 );
    
    
    % epoched button press
    epochedCortEco = squeeze(getEpochSignal(Wave.data,(trainTimesScreen + artifact_end),(trainTimesScreen + sampsEnd)));
    
    
    figure
    for chan = 1:size(epochedCortEco,2)
        % channel of interest, plot mean
        t_epoch = [0:size(epochedCortEco,1)-1]/eco_fs;
        
        exampChan = mean(squeeze(epochedCortEco(:,chan,:)),2);
        plot(1e3*t_epoch,exampChan);
        xlabel('time (ms)')
        ylabel('Voltage (V)')
        title(['Raw data for Channel ', num2str(chan)])
        vline(mean(cort))
%         pause(1)
    end
    
    % trial by trial notch and extract power 
    
    % from quick screen
    
    trialHG = zeros(size(epochedCortEco,1),size(epochedCortEco,2),size(epochedCortEco,3));
    trialBeta = zeros(size(epochedCortEco,1),size(epochedCortEco,2),size(epochedCortEco,3));

    for i = 1:size(epochedCortEco,3)
   % notch filter to eliminate line noise
    sig = squeeze(epochedCortEco(:,:,i));
    
    fprintf('notch filtering\n');
    sig = notch(sig, [60 120 180], eco_fs, 4);
    
    % extract HG and Beta power bands
    fprintf('extracting HG power\n');
    logHGPower = log(hilbAmp(sig, [70 200], eco_fs).^2);
    fprintf('extracting Beta power\n');
    logBetaPower = log(hilbAmp(sig, [12 30], eco_fs).^2);
    
    trialHG(:,:,i) = logHGPower;
    trialBeta(:,:,i) = logBetaPower;
    end
    
        % sort by reaction time 
    [sorted,indexes] = sort(cort);
    
    trialHGsort = trialHG(:,:,indexes);
        trialBetasort = trialBeta(:,:,indexes);

    
    % plot example channel
    
    chan = 54;
    
    figure
    imagesc(1e3*t_epoch,trial,squeeze(trialHGsort(:,chan,:))')
    axis xy
    xlabel('Time (ms)')
    ylabel('Trial')
    title(['High Gamma power in channel ', num2str(chan)])
    
    figure
    imagesc(1e3*t_epoch,trial,squeeze(trialBetasort(:,chan,:))')
    axis xy
    xlabel('Time (ms)')
    ylabel('Trial')
    title(['Beta power in channel ', num2str(chan)])
    

    
    
    %     % from Stavros code - filter ?
    %
    %     % assume stimulation is over by 210 ms or so, so select segments where
    %     % t>210
    %
    %     % attempt to fit exponential decay to post-stim segment
    %     %pp = smooth(c2plot(supsam(2)+1:end),smoothfw);
    %     pp = exampChan((t_epoch>210):end);
    %     [f,gof] = fit([1:length(pp)]',pp','exp2'); % 2-term exponential
    %     if gof.adjrsquare>0.8
    %         % subtract exp function if R2>0.8
    %         ppc = pp' - f([1:length(pp)]');
    %     else
    %         ppc = pp';
    %     end
    %
    %     % smooth post-stim segment using Savitzky-Golay filtering
    %     smppc = sgolayfilt(ppc,5,81);
    %
    %     % end of stavros code
    
elseif s ==2
    
    
    % pick condition type where stimulation was delivered
    if s == 1
        trainTimesCond1 = trainTimes(condType==0);
    elseif s == 2
        trainTimesCond1 = trainTimes(condType==0 | condType==1);
    end
    
    sampsEnd = round(2*fs_stim);
    
    % epoched button press
    epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEnd)));
    
    figure
    t_epoch = [0:size(epochedButton,1)-1]/fs_stim;
    plot(t_epoch,epochedButton);
end

%% tactor brain data