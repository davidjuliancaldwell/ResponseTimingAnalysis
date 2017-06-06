%% look at the neural data in response to cortical stimulation


% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

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
clearvars -except ECO1 ECO2 ECO3 Tact sid block s
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

respLo = 0.1;
respHi = 1;

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


%% cortical brain data
% ui box for input
prompt = {'Channel of interest?','Trim ends?','condition'};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'19','y','4'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

% options for this subject are -1,0,1,2,3,4,5
chanInt = str2num(answer{1});
trimEnds = answer{2};
condIntAns = str2num(answer{3});
condInt = find(uniqueCond==condIntAns);
% where to begin plotting with artifact
%artifact_end = round(0.05*eco_fs);
artifact_end = 0;

%where to end plotting
sampsEnd = round(2*eco_fs);

%presamps - where to begin looking for "rest" period (500 ms before?)
presamps = round(0.5*eco_fs);

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
if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
    
    %stim_train_length = condIntAns;
    stim_train_length = 2000;
    post_stim = 2000;
    samps_post_stim = round(stim_train_length/1e3*eco_fs);
    
    pre_stim = 1000;
    samps_pre_stim = round(pre_stim/1e3*eco_fs);
    
    epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{condInt})-samps_pre_stim,(trainTimesCellThresh{condInt}+ samps_post_stim)));
    
    response = buttonLocsThresh{condInt};
end

if condIntAns == -1
    
    post_stim = 2000;
    samps_post_stim = round(post_stim/1e3*eco_fs);
    
    pre_stim = 1000;
    samps_pre_stim = round(pre_stim/1e3*eco_fs);
    
    %response = buttonLocsThresh{condInt} + tactorLocsVec;
    response = buttonLocsThresh{condInt};
    response_samps = round(tactorLocsVec*eco_fs);
    epochedCortEco = squeeze(getEpochSignal(data,((trainTimesCellThresh{condInt}+response_samps)-samps_pre_stim),((trainTimesCellThresh{condInt}+response_samps)+ samps_post_stim)));
    
    
end

%t_epoch = [1:size(epochedCortEco,1)]/eco_fs;
t_epoch = (-samps_pre_stim:samps_post_stim-1)/eco_fs;

exampChan = mean(squeeze(epochedCortEco(:,chanInt,:)),2);

figure
plot(1e3*t_epoch,exampChan);

clear exampChan

%%
% Process the signal with ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
    
    if s ==1
        scale_factor = 50;
        numComponentsSearch = 10;
    elseif s == 2
        scale_factor = 50;
        numComponentsSearch = 10;
        
        %      scale_factor = 1000;
        %      numComponentsSearch = 20;
        
    end
    plotIt = true;
    %stimChans = [20 29 57 58]; already excluded 57 and 58 from analysis up
    %top
    
    stimChans = [20 29];
    meanSub = 1;
    %
    % [subtracted_sig_matrixS_I, subtracted_sig_cellS_I,recon_artifact_matrix,recon_artifact,t] = ...
    %     ica_artifact_remove_train(t_epoch,epochedCortEco,stimChans,eco_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub);
    
    orderPoly = 6;
    [processedSig,~,~,~,t] = ...
        ica_artifact_remove_train(t_epoch,epochedCortEco,stimChans,eco_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub,orderPoly);
    
    stimTime = zeros(size(subtracted_sig_matrixS_I,3));
elseif (condIntAns == -1)
    processedSig = [];
    meanSub = 1;
    
    if meanSub == 1
        for i = 1:size(epochedCortEco,2)
            for j = 1:size(epochedCortEco,3)
                data_int_temp = squeeze(epochedCortEco(:,i,j));
                [p,s,mu] = polyfit((1:numel(data_int_temp))',data_int_temp,6);
                f_y = polyval(p,(1:numel(data_int_temp))',[],mu);
                
                % subtract poly fit
                processedSig(:,i,j) = data_int_temp - f_y;
                
            end
            
        end
    else
        processedSig = epochedCortEco;
    end
    
    %subtracted_sig_matrixS_I = epochedCortEco;
    %stimTime = 1e3*tactorLocsVec; %
    stimTime = zeros(size(processedSig,3)); % it is centered around zero now 
    
    % exclude bad channel
    
    stimChans = [1 9 24 29 32];
    lnFreqs = [60 120 180 240 300 360 420 480 540];
    order = 3;
    processedSig = notch(processedSig,lnFreqs,eco_fs,order);
    
end

return

%% PROCESS THE DATA
% process the wavelet using morlet process and PLV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial by trial wavelet decomp, PLV

%%%%%% PLV
freq_range = [8 12];
[plv] = plvWrapper(processedSig,eco_fs,freq_range,stimChans);
%%
%%%%%%% wavelet
time_res = 0.050; % 50 ms bins

[powerout,f_morlet,t_morlet,~] = waveletWrapper(processedSig,eco_fs,time_res,stimChans);

t_morlet = linspace(-pre_stim,post_stim,length(t_morlet))/1e3;

%% Visualize wavelets

% example wavelet decomp

chanInt = 22;


for i = 1:size(powerout,4)
    figure;
    subplot(3,1,1)
    imagesc(1e3*t_morlet,f_morlet,powerout(:,:,chanInt,i));
    axis xy;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    title(['Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(stimTime(i),'r','stim')
    vline(1e3*response(i),'g','response')
    
    %figure;
    h1 = subplot(3,1,2)
    plot(1e3*t_epoch,processedSig(:,chanInt,i))
    vline(stimTime(i),'r','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Processed Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylim_h1 = ylim;
    
    h2 = subplot(3,1,3)
    plot(1e3*t_epoch,epochedCortEco(:,chanInt,i))
    vline(stimTime(i),'r','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Raw Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylim(ylim_h1)
    
    linkaxes([h1,h2],'xy')
    
end
% sort by reaction time


%% Visualize PLV

% chan 1 is the lower valued chan, so e.g., 17, 20
chan1 = 32;
chan2 = 40;
figure;

% probably want to discard the number of samples for the order of the
% filter. So for alpha

desired_f = 10;
period = 1/desired_f;
time_4oscil = period*4; % time total in seconds
order = round(time_4oscil*eco_fs);
samps_discard = order;


plot(t_epoch, squeeze(plv(:, chan1, chan2)));
ylim([0 1])
vline(0)
xlabel('Time (s)');
ylabel('Plase Locking Value');
title(['PLV between Channel ' num2str(chan1) ' and ' num2str(chan2)])
%%
%%%%%%% find max PLV in time window 

t1 = 0;
t2 =  0.2;

t_sub = t_epoch((t_epoch>t1 & t_epoch<t2));
plv_sub = plv((t_epoch>t1 & t_epoch<t2),:,:);

figure
max_plv = squeeze(max(plv_sub,[],1));
imagesc(max_plv)
h = colorbar
h.Limits = [0 1]
xlabel('Channel')
ylabel('Channel')
title(['Max Phase Locking Value between T  = ' num2str(t1) ' seconds and T = ' num2str(t2) ' seconds'])


return



% Below here is old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


% plot channel of interest
figure

% channel of interest, plot mean
t_epoch = [artifact_end:artifact_end+size(epochedCortEco,1)-1]/eco_fs;
t_epochPre = [-presamps:0-1]/eco_fs;

exampChanPost = mean(squeeze(epochedCortEco(:,chanInt,:)),2)-mean(exampChanPost);
exampChanPre = mean(squeeze(epochedPreStim(:,chanInt,:)),2)-mean(exampChanPre);


subplot(2,1,1)
plot(1e3*t_epoch,exampChanPost);
xlabel('time (ms)')
ylabel('Voltage (V)')
title(['Post stim Raw data for Channel ', num2str(chanInt)])

subplot(2,1,2)
plot(1e3*t_epochPre,exampChanPre);
xlabel('time (ms)')
ylabel('Voltage (V)')
title(['Pre Stim Raw data for Channel ', num2str(chanInt)])
%         pause(1)

% trial by trial notch and extract power
%%
% from quick screen

% trials
trialHG = zeros(size(epochedCortEco,1),size(epochedCortEco,2),size(epochedCortEco,3));
trialBeta = zeros(size(epochedCortEco,1),size(epochedCortEco,2),size(epochedCortEco,3));

% rest
restHG = zeros(size(epochedPreStim,1),size(epochedPreStim,2),size(epochedPreStim,3));
restBeta = zeros(size(epochedPreStim,1),size(epochedPreStim,2),size(epochedPreStim,3));

% for each trial, run the functions below - only want to filter pre and
% post
for i = 1:size(epochedCortEco,3)
    % notch filter to eliminate line noise
    sigT = squeeze(epochedCortEco(:,:,i));
    sigR = squeeze(epochedPreStim(:,:,i));
    
    sigT = notch(sigT, [60 120 180], eco_fs, 4);
    sigR = notch(sigR, [60 120 180], eco_fs, 4);
    
    % extract HG and Beta power bands
    logHGPowerT = log(hilbAmp(sigT, [70 200], eco_fs).^2);
    logBetaPowerT = log(hilbAmp(sigT, [12 30], eco_fs).^2);
    
    logHGPowerR = log(hilbAmp(sigR, [70 200], eco_fs).^2);
    logBetaPowerR = log(hilbAmp(sigR, [12 30], eco_fs).^2);
    
    trialHG(:,:,i) = logHGPowerT;
    trialBeta(:,:,i) = logBetaPowerT;
    
    
    restHG(:,:,i) = logHGPowerR;
    restBeta(:,:,i) = logBetaPowerR;
end

% part for getting rid of funny filtering at beginning and end
t_min = 0.2; % in seconds
t_max = 1.5; % in seconds

t_minPre = -0.4;
t_maxPre = -0.1;

if strcmp(trimEnds,'y')
    
    %temporary ones
    
    t_epochT = t_epoch(t_epoch>t_min & t_epoch<t_max);
    t_epochPreT = t_epochPre(t_epochPre>t_minPre & t_epochPre<t_maxPre);
    trialHGT = trialHG(t_epoch>t_min & t_epoch<t_max,:,:);
    trialBetaT = trialBeta(t_epoch>t_min & t_epoch<t_max,:,:);
    restHGT = restHG(t_epochPre>t_minPre & t_epochPre<t_maxPre,:,:);
    restBetaT = restBeta(t_epochPre>t_minPre & t_epochPre<t_maxPre,:,:);
    
    clear t_epoch t_epochPre trialHG trialBeta restHG restBeta
    
    t_epoch = t_epochT;
    t_epochPre = t_epochPreT;
    trialHG = trialHGT;
    trialBeta = trialBetaT;
    restHG = restHGT;
    restBeta = restBetaT;
    
end
%%
% sort by reaction time
[sorted,indexes] = sort(cort);

trialHGsort = trialHG(:,:,indexes);
trialBetasort = trialBeta(:,:,indexes);

restHGsort = restHG(:,:,indexes);
restBetasort = restBeta(:,:,indexes);


% find significant differences for all channels

restHG_ave = squeeze(mean(restHGsort,1));
restBeta_ave = squeeze(mean(restBetasort,1));
trialHG_ave = squeeze(mean(trialHGsort,1));
trialBeta_ave = squeeze(mean(trialBetasort,1));

numChans = size(raw_eco,2);
chans = 1:numChans;

ptarg = 0.05 / numChans;

HGSigs = ttest2(restHG_ave, trialHG_ave, ptarg, 'r', 'unequal', 2);
BetaSigs = ttest2(restBeta_ave, trialBeta_ave, ptarg, 'r', 'unequal', 2);

HGSigs = HGSigs == 1; % make boolean
BetaSigs = BetaSigs == 1; % make boolean


HGRSAs = signedSquaredXCorrValue(restHG_ave, trialHG_ave, 2);
BetaRSAs = signedSquaredXCorrValue(restBeta_ave, trialBeta_ave, 2);

%HG figure
figure
plot(chans, HGRSAs);
hold on;
plot(chans(HGSigs), HGRSAs(HGSigs), '*');

xlabel('channel number');
ylabel('R^2');
title('Aggregated HG Response');
legend('aggregate activity');

% Beta
figure;
plot(chans, BetaRSAs);
hold on;
plot(chans(BetaSigs), BetaRSAs(BetaSigs), '*');

xlabel('channel number');
ylabel('R^2');
title('Aggregated Beta Response');
legend('aggregate activity');


% plot example channel
%% plot channel of interest


trial = 1:length(sorted);

figure
imagesc(1e3*t_epoch,trial,squeeze(trialHGsort(:,chanInt,:))')
axis xy
xlabel('Time (ms)')
ylabel('Trial')
title(['High Gamma power in channel ', num2str(chanInt)])

figure
imagesc(1e3*t_epoch,trial,squeeze(trialBetasort(:,chanInt,:))')
axis xy
xlabel('Time (ms)')
ylabel('Trial')
title(['Beta power in channel ', num2str(chanInt)])

%%


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

%% time frequency wavelet

fw = 1:1:200;


[C_post, ~, C_totPost, ~] = time_frequency_wavelet(squeeze(epochedCortEco(t_epoch>t_min & t_epoch<t_max,chanInt,:)), fw, eco_fs, 1, 1, 'CPUtest');
[C_pre, ~, C_totPre, ~] = time_frequency_wavelet(squeeze(epochedPreStim(t_epochPre>t_minPre &t_epochPre<t_maxPre,chanInt,:)), fw, eco_fs, 1, 1, 'CPUtest');
C_norm = normalize_data(C_post',C_pre');

figure
imagesc(1e3*t_epoch,fw,C_norm);
set_colormap_threshold(gcf, [-1 1], [-7 7], [1 1 1]);
colorbar
axis xy
xlabel('time (ms)');
ylabel('frequency (hz)');
title(['Normalized wavelet data for Channel ', num2str(chanInt)])


figure
imagesc(1e3*t_epoch,fw,C_post');
%set_colormap_threshold(gcf, [-1 1], [-7 7], [1 1 1]);
colorbar
axis xy
xlabel('time (ms)');
ylabel('frequency (hz)');
title(['Normalized wavelet data for Channel ', num2str(chanInt)])



figure
imagesc(1e3*t_epochPre,fw,C_pre');
%set_colormap_threshold(gcf, [-1 1], [-7 7], [1 1 1]);
colorbar
axis xy
xlabel('time (ms)');
ylabel('frequency (hz)');
title(['Normalized wavelet data for Channel ', num2str(chanInt)])




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

%% tactor brain data