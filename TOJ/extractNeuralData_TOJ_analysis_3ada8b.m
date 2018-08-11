% this is from my z_constants
close all;clear all;clc
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\3ada8b\data\d10\MATLAB_conversions\3ada8b_TOJ';
sid = SIDS{6};
% load in data
folder_data = strcat(DATA_DIR);

load(fullfile(folder_data,'TOJ-1.mat'))
block = '1';
ECoG = 4*cat(2,ECO1.data,ECO2.data,ECO3.data);
%stim = Stim.data;
%tact = Tact.data;


ecoFs = ECO1.info.SamplingRateHz;
% fsStim = Stim.info.SamplingRateHz;
% fsTact = Tact.info.SamplingRateHz;
clearvars ECO1 ECO2 ECO3 Stim Tact

%  % get rid of bad channels
bads = [];
% bads = [79:98];
bads = [65:128];
goodVec = logical(ones(size(ECoG,2),1));
goodVec(bads) = 0;
ECoG = ECoG(:,goodVec);

% convert sampling rate

% fac = fsTact/ecoFs;
fac = 2;

stimChans = [4 3];


load(fullfile('3ada8b_TOJ_matlab.mat'));

%% define what to epoch around for centering on stim trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trainTimesConverted = round(trainTimes/fac) + round(ecoFs*0.8) ; % convert from the tactor sampling rate to the eco sampling rate, it is centered around the stim train

preTime = 1000; % ms
postTime = 2000; % ms
preSamps = round(preTime*ecoFs/1e3); % convert time to samps
postSamps = round(postTime*ecoFs/1e3); % convert time to samps
tEpoch = [-preSamps:postSamps-1]/ecoFs;

% get signal epochs
epochedECoG = getEpochSignal(ECoG,trainTimesConverted -preSamps,trainTimesConverted +postSamps); % break up the ECoG into chunks

%% behavioral data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delays manually used
delays = tactorStimDiff;

% tactor = 0
% stim = 1
firstFeel = whichPerceived;


%% process artifacts centered on stim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trainDuration = [0 500]; % this is how long the stimulation train was
xlims = [-200 2000]; % these are the x limits to visualize in plots
chanIntList = [2 10 11 12 18 19 20 21 13 5 6 14 22];
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
post = 1; % started with 0.2
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

[processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(epochedECoG,'type',type,...
    'fs',ecoFs,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
    'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
    'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,'minDuration',minDuration);
%%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(processedSig,epochedECoG,'fs',ecoFs,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')

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

tMorlet = linspace(-preTime,postTime,length(tMorlet))/1e3;

% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedData] = normalize_spectrogram(dataRef,powerout);
%%
individual = 0;
average = 1;
chanIntList = [1 2 3 4 5 12 13 30 33];
%chanIntList = 3;
trainDuration = [];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-300 300];
vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)

%%
% chanIntList = chanInt;
for chanInt = chanIntList
    visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,epochedECoG,chanInt,stimTime,response,individual,average)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HGPowerWavelet = squeeze(mean(squeeze(normalizedData(fMorlet < 150 & fMorlet > 70,:,:,:)),1));

%%
vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
%% hilb amp HG
processedSigHG = zeros(size(processedSig));
for trial = 1:size(processedSig,3)
    [amp] = log(hilbAmp(squeeze(processedSig(:,:,trial)), [70 150], fsData).^2);
    processedSigHG(:,:,trial) = amp;
end


%% process artifacts focused on tactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circular shift already processed signal

sigShifted = nan(size(processedSig));
for i = 1:length(tactorStimDiff)
    if ~isnan(tactorStimDiff(i))
        sigShifted(:,:,i) = circshift(processedSig(:,:,i),-round(ecoFs*tactorStimDiff(i)),1);
    end
end

avgResponseShift = nanmean(sigShifted,3);
avgResponse = nanmean(processedSig,3);
%%
chanInt = 30
figure
subplot(2,1,1)
plot(1e3*tEpoch,1e6*avgResponseShift(:,chanInt))
xlabel('time (ms)')
ylabel('voltage (\muV)')
hold on
%plot(tEpoch,avgResponseShift(:,chanInt+3))
title('aligned on digital touch probe onset')
set(gca,'fontsize',14)
subplot(2,1,2)
plot(1e3*tEpoch,1e6*avgResponse(:,chanInt))
hold on
%plot(tEpoch,avgResponse(:,chanInt+3))
title('aligned on stimulation train onset')
xlabel('time (ms)')
ylabel('voltage (\muV)')
set(gca,'fontsize',14)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(sigShifted,epochedECoG,'fs',ecoFs,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')

%% rereference
rerefMode = 'mean';
badChannels = stimChans;
processedSigReref = rereference_CAR_median(processedSig,rerefMode,badChannels);

%%
sigShiftedReref = nan(size(processedSigReref));
for i = 1:length(tactorStimDiff)
    if ~isnan(tactorStimDiff(i))
        sigShiftedReref(:,:,i) = circshift(processedSigReref(:,:,i),-round(ecoFs*tactorStimDiff(i)),1);
    end
end

avgResponseShift = nanmean(sigShiftedReref,3);
avgResponse = nanmean(processedSig,3);
%%
individual = 0;
average = 1;
chanIntList = [1 2 3 4 5 12 13 30 33];
%chanIntList = 3;
trainDuration = [];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-300 300];

vizFunc.small_multiples_time_series(sigShiftedReref,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)


%%
trainTimesTactor = trainTimes + round(tactorStimDiff*ecoFs)';

preTime = 1000; % ms
postTime = 2000; % ms
preSamps = round(preTime*ecoFs/1e3); % convert time to samps
postSamps = round(postTime*ecoFs/1e3); % convert time to samps
tEpoch = [-preSamps:postSamps-1]/ecoFs;

% get signal epochs
epochedECoGTact = getEpochSignal(ECoG,trainTimesTactor-preSamps,trainTimesTactor+postSamps); % break up the ECoG into chunks



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

% Below here is tactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    'fs',ecoFs,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
    'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
    'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,'minDuration',minDuration);
%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(processedSigTact,epochedECoGTact,'fs',ecoFs,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')
