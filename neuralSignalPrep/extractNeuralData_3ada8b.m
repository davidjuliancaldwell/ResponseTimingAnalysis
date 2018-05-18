%% load in subject
%close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\3ada8b\data\d9\MATLAB_conversions\3ada8b_ResponseTiming';
sid = SIDS{6};

% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, '3ada8b'))
    folder_data = strcat(DATA_DIR,'\2fd831');
    
    if s == 1
        load(fullfile(folder_data,'responseTiming-1.mat'))
        block = '1';
    elseif s == 2
        load(fullfile(folder_data,'responseTiming-2.mat'))
        block = '2';
    end
    
end

plotIt = 1;

%% neural data

% the eco data is crashing it right now
clearvars -except ECO1 ECO2 ECO3 Tact sid block s
eco1 = ECO1.data;
fsData = ECO1.info.SamplingRateHz;
ecoFs = fsData;
clear ECO1
eco2 = ECO2.data;
clear ECO2

eco3 = ECO3.data;
clear ECO3


ECoG = [eco1(1:end-1,:) eco2(1:end-1,:) eco3];
clearvars eco1 eco2 eco3

% only 64 channels grid

% get rid of bad channels
bads = [79:98];
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
stimChans = [1 2 16 24];

%% cortical brain data
%
% where to begin plotting with artifact
%artifact_end = round(0.05*eco_fs);
artifact_end = 0;

%where to end plotting
sampsEnd = round(2*ecoFs);

%presamps - where to begin looking for "rest" period (500 ms before?)
presamps = round(0.5*ecoFs);

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
    
    if (condIntAns == 0 || condIntAns == 1 || condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
        
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
        epochedCortEco = squeeze(getEpochSignal(ECoG,(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)-sampsPreStim),(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)+ sampsPostStim)));
        tactData = decimate(Tact.data,2)';
        epochedTactorNew = squeeze(getEpochSignal(tactData,(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)-sampsPreStim),(trainTimesCellThresh{i}(response_mask)+responseSamps(response_mask)+ sampsPostStim)));
    end
    
    epochedCortEco_cell{i} = epochedCortEco;
end
tEpoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;

%%
 % ui box for input
prompt = {'Channel of interest?','condition'};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'8','-1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

% options for this subject are 0 1
chanInt = str2num(answer{1});
condIntAns = str2num(answer{2});
condInt = find(uniqueCond==condIntAns);

dataInt = epochedCortEco_cell{condInt};

exampChan = mean(squeeze(dataInt(:,chanInt,:)),2);

figure
plot(1e3*tEpoch,exampChan);
xlim([-1000 2000])
ylim([-10e-5 10e-5])
title(['Subject ' sid ' Channel ' num2str(chanInt) ' Condition ' num2str(condIntAns)])
%clear exampChan

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
    
    meanSub = 1;
    %
    % [subtracted_sig_matrixS_I, subtracted_sig_cellS_I,recon_artifact_matrix,recon_artifact,t] = ...
    %     ica_artifact_remove_train(t_epoch,epochedCortEco,stimChans,eco_fs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub);
    
    orderPoly = 6;
    [processedSig,~,~,~,t] = ...
        ica_artifact_remove_train(tEpoch,epochedCortEco,stimChans,ecoFs,scale_factor,numComponentsSearch,plotIt,chanInt,meanSub,orderPoly);
    
    stimTime = zeros(size(subtracted_sig_matrixS_I,3));
elseif (condIntAns == -1)
    processedSig = [];
    meanSub = 0;
    
    if meanSub == 1
        for i = 1:size(dataInt,2)
            for j = 1:size(dataInt,3)
                data_int_temp = squeeze(dataInt(:,i,j));
                [p,s,mu] = polyfit((1:numel(data_int_temp))',data_int_temp,6);
                f_y = polyval(p,(1:numel(data_int_temp))',[],mu);
                
                % subtract poly fit
                processedSig(:,i,j) = data_int_temp - f_y;
                
            end
            
        end
    else
        processedSig = dataInt;
    end
    
    %subtracted_sig_matrixS_I = epochedCortEco;
    %stimTime = 1e3*tactorLocsVec; %
    stimTime = zeros(size(processedSig,3)); % it is centered around zero now
    
    % exclude bad channel
    
    %stimChans = [1 9 24 29 32]; this might be from the other subject???
    %2fd831
    %lnFreqs = [60 120 180 240 300 360 420 480 540];
    %order = 3;
   % processedSig = notch(processedSig,lnFreqs,ecoFs,order);
    
end
%%
chanIntList = [7 8 15 23 31 32 39 40 47 48];

for ind = chanIntList
    tEpoch = (-sampsPreStim:sampsPostStim-1)/ecoFs;
    
    exampChan = mean(squeeze(processedSig(:,ind,:)),2);
    
    figure
    subplot(2,1,1)
    plot(1e3*tEpoch,exampChan);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    subplot(2,1,2)
    exampChan = mean(squeeze(dataInt(:,ind,:)),2);
    plot(1e3*tEpoch,exampChan);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    title(['Raw Signal Average - Channel ' num2str(ind)])
    
    
    clear exampChan
end

% 9-28-2017 - save intermediate data for plotting response map

%%
sig = processedSig;
avgResponse = mean(sig,3);

smallMultiples_responseTiming(avgResponse,tEpoch,'type1',stimChans,'type2',0,'average',1)

%save([sid '_block' num2str(block) '_postICAProcessData'],processsedSig,'-v7.3')

%# sort by rxn time

% load subject data, need sid still
%load([sid,'_compareResponse_block_',block,'.mat'])

% for i = 1:length(uniqueCond)
%     % 12-10-2016
%     respLo = 0.150;
%     respHi = 1;
%
%
%     trim = buttonLocs{i};
%     trim = trim(trim>respLo & trim<respHi);
%     zTrim = zscore(trim);
%     buttonLocsThresh{i} = 1e3.*trim(abs(zTrim)<3);
%     %buttonLocsThresh{i} = 1e3.*trim;
%
% end

%%
rxnTimes = buttonLocs{condInt};
[sorted,indexes] = sort(rxnTimes);
sortedSig = sig(:,:,indexes);

% calculate first response time, shift others based from there
sortedBasedOffFirst = round((sorted - sorted(1))*ecoFs);
sigShifted = sortedSig;

num_non_nan = sum(~isnan(rxnTimes));

for i = 2:num_non_nan
    sigShifted(:,:,i) = circshift(sortedSig(:,:,i),sortedBasedOffFirst(i),1);
end

avgResponse = mean(sigShifted,3);
tShift = tEpoch - sorted(1)/ecoFs;
smallMultiples_responseTiming(avgResponse,tShift,'type1',stimChans,'type2',0,'average',1)


%% PROCESS THE DATA
% process the wavelet using morlet process and PLV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial by trial wavelet decomp, PLV

%%%%%% PLV
freqRange = [8 12];
[plv] = plvWrapper(processedSig,ecoFs,freqRange,stimChans);

%%%%%%% wavelet
timeRes = 0.050; % 50 ms bins

[powerout,fMorlet,tMorlet,~] = waveletWrapper(processedSig,ecoFs,timeRes,stimChans);

tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;

%% Visualize wavelets

% example wavelet decomp
%trialInt = 20; % which channel to check out
chanInt = 15;


response = buttonLocs{condInt};

for i = 1:size(powerout,4)
    totalFig = figure;
    totalFig.Units = 'inches';
    totalFig.Position = [12.1806 3.4931 6.0833 7.8056];
    subplot(3,1,1);
    imagesc(1e3*tMorlet,fMorlet,powerout(:,:,chanInt,i));
    axis xy;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    title(['Wavelet decomposition Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(stimTime(i),'r','stim')
    vline(1e3*response(i),'g','response')
    xlim([-200 1000]);
    set(gca,'fontsize',14)
    
    
    
    %figure;
    h1 = subplot(3,1,2);
    plot(1e3*tEpoch,1e6*processedSig(:,chanInt,i))
    vline(stimTime(i),'r','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Processed Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylims = [-(max(abs(1e6*processedSig(:,chanInt,i))) + 10) (max(abs(1e6*processedSig(:,chanInt,i))) + 10)];
    ylim(ylims);
    ylim_h1 = ylims;
    xlim([-200 1000]);
    set(gca,'fontsize',14)
    
    
    h2 = subplot(3,1,3);
    plot(1e3*tEpoch,1e6*dataInt(:,chanInt,i))
    vline(stimTime(i),'r','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Raw Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylim(ylim_h1);
    xlim([-200 1000]);
    set(gca,'fontsize',14);
    
    
    
    linkaxes([h1,h2],'xy');
    
end
%%
% plot average spectrogram
avgPower = mean(powerout,4);
smallMultiples_responseTiming_spectrogram(avgPower,tMorlet,fMorlet,'type1',stimChans,'type2',0,'average',1)

% sort by reaction time


%% Visualize PLV

% chan 1 is the lower valued chan, so e.g., 17, 20
chan1 = 34;
chan2 = 42;
figure;

% probably want to discard the number of samples for the order of the
% filter. So for alpha

desiredF = 10;
period = 1/desiredF;
time4oscil = period*4; % time total in seconds
order = round(time4oscil*ecoFs);
samps_discard = order;


plot(tEpoch, squeeze(plv(:, chan1, chan2)));
ylim([0 1])
vline(0)
xlabel('Time (s)');
ylabel('Plase Locking Value');
title(['PLV between Channel ' num2str(chan1) ' and ' num2str(chan2)])

t1 = 0;
t2 =  0.2;

tSub = tEpoch((tEpoch>t1 & tEpoch<t2));
plvSub = plv((tEpoch>t1 & tEpoch<t2),:,:);

figure
maxPlv = squeeze(max(plvSub,[],1));
imagesc(maxPlv)
h = colorbar;
h.Limits = [0 1];
xlabel('Channel')
ylabel('Channel')
title(['Max Phase Locking Value between T  = ' num2str(t1) ' seconds and T = ' num2str(t2) ' seconds'])


return

% Below here is old
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Try wavelet denoising

nc = 15;
br = 8;
er = 15;
plotIt = 1;
ec = 10;
et = 10;


processedSig = wavelet_denoise(epochedCortEco,'numComponents',nc,'beginRecon',br,'endRecon',er,...
    'plotIt',plotIt,'exampChan',ec,'exampTrial',et);

stimTime = zeros(size(processedSig,3));
stimChans = [9 17 50 58 ];

% median subtract
processedSig_medianS = processedSig - repmat(mean(processedSig,2),[1,size(processedSig,2),1]);


%%
% plot channel of interest
figure

% channel of interest, plot mean
t_epoch = [artifact_end:artifact_end+size(epochedCortEco,1)-1]/fs_data;
t_epochPre = [-presamps:0-1]/fs_data;

exampChanPost = mean(squeeze(epochedCortEco(:,chanInt,:)),2);
exampChanPre = mean(squeeze(epochedPreStim(:,chanInt,:)),2);

%exampChanPost = mean(squeeze(epochedCortEco(:,chanInt,:)),2)-mean(exampChanPost);
%exampChanPre = mean(squeeze(epochedPreStim(:,chanInt,:)),2)-mean(exampChanPre);


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
    
    sigT = notch(sigT, [60 120 180], fs_data, 4);
    sigR = notch(sigR, [60 120 180], fs_data, 4);
    
    % extract HG and Beta power bands
    logHGPowerT = log(hilbAmp(sigT, [70 200], fs_data).^2);
    logBetaPowerT = log(hilbAmp(sigT, [12 30], fs_data).^2);
    
    logHGPowerR = log(hilbAmp(sigR, [70 200], fs_data).^2);
    logBetaPowerR = log(hilbAmp(sigR, [12 30], fs_data).^2);
    
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

fw = 1:3:200;

[C_post, ~, C_totPost, ~] = time_frequency_wavelet(squeeze(epochedCortEco(t_epoch>t_min & t_epoch<t_max,chanInt,:)), fw, fs_data, 1, 1, 'CPUtest');
[C_pre, ~, C_totPre, ~] = time_frequency_wavelet(squeeze(epochedPreStim(t_epochPre>t_minPre &t_epochPre<t_maxPre,chanInt,:)), fw, fs_data, 1, 1, 'CPUtest');
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
