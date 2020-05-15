% this is from my z_constants
close all;clear all;clc
Z_ConstantsStimResponse;

subjdir = getenv('SUBJECT_DIR');
sid = SIDS{5};
DATA_DIR = fullfile(subjdir,sid,'\data\d7\Converted_Matlab\TOJ');

for s = 1:2
    % load in data
    folder_data = strcat(DATA_DIR);
    if s == 1
        load(fullfile(folder_data,'TOJ-1.mat'))
        block = '1';
        ECoG = cat(2,ECO1.data,ECO2.data,ECO3.data);
        %stim = Stim.data;
        %tact = Tact.data;
        
        
        clearvars ECO1 ECO2 ECO3 Stim tact
    elseif s == 2
        load(fullfile(folder_data,'TOJ-2.mat'))
        block = '2';
        ECoG = 4*[ECoG; cat(2,ECO1.data,ECO2.data,ECO3.data)];
        %stim = [stim; Stim.data];
        fsData = ECO1.info.SamplingRateHz;
       % fsStim = Stim.info.SamplingRateHz;
       % fsTact = Tact.info.SamplingRateHz;
        clearvars ECO1 ECO2 ECO3 Stim Tact
        
        % get rid of bad channels
        bads = [79:98];
        goodVec = logical(ones(size(ECoG,2),1));
        goodVec(bads) = 0;
        ECoG = ECoG(:,goodVec);
        
        % convert sampling rate
        
       % fac = fsTact/fsData;
        fac = 2;
        % stim chans, 16/24
        
        stimChans = [16 24];
        
    end
    
end
load(fullfile('a1355e_TOJ_matlab.mat'));

%% define what to epoch around for centering on stim trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trainTimes = round(trainTimes/fac) + round(1*fsData); % convert from the tactor sampling rate to the eco sampling rate, add one second to center it around the stim train

preTime = 1000; % ms
postTime = 2000; % ms
preSamps = round(preTime*fsData/1e3); % convert time to samps
postSamps = round(postTime*fsData/1e3); % convert time to samps
tEpoch = [-preSamps:postSamps-1]/fsData;

% get signal epochs
epochedECoG = getEpochSignal(ECoG,trainTimes-preSamps,trainTimes+postSamps); % break up the ECoG into chunks

%% behavioral data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delays manually used
delays = [2700,2550,2500,2600,2700,2550,2700,2550,2550,2700,2800,2550,2400,2400,2800,2850,2400,2850,2400,2580,2400]-2700;

% tactor = 0
% stim = 1
% same = 2

firstFeel = [1,2,1,1,0,2,0,1,1,1,0,2,1,1,1,0,1,0,1,0,1]; % which did he feel first
firstFeel([8,9,14]+1) = 2;
% ones with slight edge for stim put back as same

% process artifacts centered on stim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trainDuration = [0 500]; % this is how long the stimulation train was
xlims = [-200 2000]; % these are the x limits to visualize in plots
chanIntList = [6 7 8 15 22 23 32];;
% these are the channels of interest to visualize in closer detail
minDuration = 0.5; % minimum duration of artifact in ms


%% parameters
type = 'dictionary';

useFixedEnd = 0;
fixedDistance = 4;
plotIt = 0;
pre = 0.8; % started with 1
post = 1; % started with 0.2
% 2.8, 1, 0.5 was 3/19/2018

% these are the metrics used if the dictionary method is selected. The
% options are 'eucl', 'cosine', 'corr', for either euclidean distance,
% cosine similarity, or correlation for clustering and template matching.
distanceMetricDbscan = 'eucl';
distanceMetricSigMatch = 'corr';
amntPreAverage = 3;
normalize = 'preAverage';
%normalize = 'firstSamp';

onsetThreshold = 1.5;
recoverExp = 0;
threshVoltageCut = 75;
threshDiffCut = 75;
expThreshVoltageCut = 95;
expThreshDiffCut = 95;
bracketRange = [-6:2];
chanInt = 23;
minPts = 3;
minClustSize = 5;
outlierThresh = 0.95;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(epochedECoG,'type',type,...
    'fs',fsData,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,...
    'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
    'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
    'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,...
    'minDuration',minDuration,'bracketRange',bracketRange,'threshVoltageCut',threshVoltageCut,...
    'threshDiffCut',threshDiffCut,'expThreshVoltageCut',expThreshVoltageCut,...
    'expThreshDiffCut',expThreshDiffCut,'onsetThreshold',onsetThreshold,'chanInt',chanInt,...
    'minPts',minPts,'minClustSize',minClustSize,'outlierThresh',outlierThresh);

%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(processedSig,epochedECoG,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')


%% PROCESS THE DATA
% process the wavelet using morlet process and PLV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial by trial wavelet decomp, PLV

%%%%%% PLV
%freqRange = [8 12];
%[plv] = plvWrapper(processedSig,fsData,freqRange,stimChans);

% %%%%%%% wavelet
% timeRes = 0.050; % 50 ms bins
% 
% [powerout,fMorlet,tMorlet,~] = waveletWrapper(processedSig,fsData,timeRes,stimChans);
% 
% tMorlet = linspace(-preTime,postTime,length(tMorlet))/1e3;
% 
% 
% % normalize data
% dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
% %
% [normalizedData] = normalize_spectrogram_wavelet(dataRef,powerout);


%%%%%% wavelet
fprintf(['-------Beginning wavelet analysis-------- \n'])

timeRes = 0.01; % 10 ms bins

% [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
[powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
%
fprintf(['-------Ending wavelet analysis-------- \n'])


%% rereference
rerefMode = 'mean';
badChannels = stimChans;
processedSigReref = rereference_CAR_median(processedSig,rerefMode,badChannels);

%%
individual = 0;
average = 1;
%chanIntList = 3;
trainDuration = [];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-1 1];
vizFunc.small_multiples_time_series(processedSigReref,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)

%
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,epochedECoG,chanInt,individual,average,xlims)
end

%%
% stimTime = 0;
% % chanIntList = chanInt;
% for chanInt = chanIntList
%     visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
%         tEpoch,epochedECoG,chanInt,stimTime,response,individual,average)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HGPowerWavelet = squeeze(mean(squeeze(normalizedData(fMorlet < 150 & fMorlet > 70,:,:,:)),1));
%
vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);

%% process artifacts focused on tactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circular shift already processed signal

sigShifted = nan(size(processedSig));
for i = 1:length(tactorStimDiff)
    if ~isnan(tactorStimDiff(i))
        sigShifted(:,:,i) = circshift(processedSig(:,:,i),-round(fsData*tactorStimDiff(i)),1);
    end
end

avgResponseShift = nanmean(sigShifted,3);
avgResponse = nanmean(processedSig,3);

%% Visualize PLV

% chan 1 is the lower valued chan, so e.g., 17, 20
chan1 = 15;
chan2 = 23;
figure;

% probably want to discard the number of samples for the order of the
% filter. So for alpha

desiredF = 10;
period = 1/desiredF;
time4oscil = period*4; % time total in seconds
order = round(time4oscil*fsData);
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

% Below here is tactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% process artifacts focused on tactor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%stimTactDiffSamps = 

for i = 1:length(stimTactDiff)
    sigShifted(:,:,i) = circshift(sortedSig(:,:,i),sortedBasedOffFirst(i),1);
end

avgResponse = mean(sigShifted,3);
tShift = tEpoch - sorted(1)/fsData;
smallMultiples_responseTiming(avgResponse,tShift,'type1',stimChans,'type2',0,'average',1)



%%
trainTimesTactor = trainTimes + round(tactorStimDiff*fsData)';

preTime = 1000; % ms
postTime = 2000; % ms
preSamps = round(preTime*fsData/1e3); % convert time to samps
postSamps = round(postTime*fsData/1e3); % convert time to samps
tEpoch = [-preSamps:postSamps-1]/fsData;

% get signal epochs
epochedECoGTact = getEpochSignal(ECoG,trainTimesTactor-preSamps,trainTimesTactor+postSamps); % break up the ECoG into chunks


%% process artifacts centered on tactor 

trainDuration = [0 500]; % this is how long the stimulation train was
xlims = [-200 2000]; % these are the x limits to visualize in plots
chanIntList = [7 8 15 23 31 32 39 40 47 48];
%chanIntList = [8 31 32];
% these are the channels of interest to visualize in closer detail
minDuration = 0.5; % minimum duration of artifact in ms

type = 'dictionary';

useFixedEnd = 0;
%fixedDistance = 2;
fixedDistance = 4; % in ms
plotIt = 0;

%pre = 0.4096; % in ms
%post = 0.4096; % in ms

pre = 0.8; % started with 1
post = 0.5; % started with 0.2
% 2.8, 1, 0.5 was 3/19/2018

% these are the metrics used if the dictionary method is selected. The
% options are 'eucl', 'cosine', 'corr', for either euclidean distance,
% cosine similarity, or correlation for clustering and template matching.

distanceMetricDbscan = 'cosine';
distanceMetricSigMatch = 'eucl';
amntPreAverage = 3;
normalize = 'preAverage';
%normalize = 'firstSamp';

recoverExp = 0;

[processedSigTact,templateDictCellTact,templateTrialTact,startIndsTact,endIndsTact] = analyFunc.template_subtract(epochedECoGTact,'type',type,...
    'fs',fsData,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
    'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
    'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,'minDuration',minDuration);
%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(processedSigTact,epochedECoGTact,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')
