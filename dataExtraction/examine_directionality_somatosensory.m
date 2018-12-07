%% script to calculate metrics of directionality from processed somatosensory evoked responses
%
% David.J.Caldwell 12/4/2018
loadIt = 1;
if loadIt
    load('3ada8b_tactor_only_processed_data.mat')
end

%%
%%%%%%% wavelet
timeRes = 0.01; % 10 ms bins

% [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
[poweroutTactorTotal,fMorlet,tMorlet,~] = waveletWrapper(processedSigTactorTotal,fsData,timeRes,stimChans);
%%
poweroutTactorTotal = poweroutTactorTotal(:,:,:,:);

tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
% normalize data
dataRefTactorTotal = poweroutTactorTotal(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedDataTactorTotal] = normalize_spectrogram(dataRefTactorTotal,poweroutTactorTotal);
%%
individual = 0;
average = 1;
chanIntList = [1 2 3 4 5 12 13 30 33];
%chanIntList = [1 2 3 4 5 12 13 30 33 73 74 75 76 77 78 79 80];

%chanIntList = 3;
trainDuration = [];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-160 160];
vizFunc.small_multiples_time_series(processedSigTactorTotal,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
%%
startInd = 65;
for index = 2:length(Montage.Montage)
    vizFunc.small_multiples_time_series(processedSigTactorTotal(:,startInd:(startInd+Montage.Montage(index)),:),...
        tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
    startInd = startInd + Montage.Montage(index);
end
%%
% chanIntList = chanInt;
for chanInt = chanIntList
    visualize_wavelet_channel(normalizedDataTactorTotal,tMorlet,fMorlet,processedSigTactorTotal,...
        tEpoch,epochedCortEcoTactorTotal,chanInt,stimTime,response,individual,average)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HGPowerWaveletTactorTotal = squeeze(mean(squeeze(normalizedDataTactorTotal(fMorlet < 150 & fMorlet > 70,:,:,:)),1));

%%
vizFunc.small_multiples_spectrogram(normalizedDataTactorTotal(:,:,1:64,:),tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
%%
startInd = 65;
for index = 2:length(Montage.Montage)
    vizFunc.small_multiples_spectrogram(normalizedDataTactorTotal(:,:,startInd:(startInd+Montage.Montage(index)-1),:),tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
    startInd = startInd + Montage.Montage(index);
end
%% hilb amp HG
processedSigHGTotal = zeros(size(processedSigTactorTotal));
for trial = 1:size(processedSigTactorTotal,3)
    [amp] = log(hilbAmp(squeeze(processedSigTactorTotal(:,:,trial)), [70 150], fsData).^2);
    processedSigHGTotal(:,:,trial) = amp;
end
%%
processedSigHGTotal = zeros(size(processedSigTactorTotal));
for trial = 1:size(processedSigTactorTotal,3)
    [amp] = (hilbAmp(squeeze(processedSigTactorTotal(:,:,trial)), [70 300], fsData).^2);
    processedSigHGnoLogTotal(:,:,trial) = amp;
end
%%
chanInt = 2;
figure
subplot(2,1,1)
plot(1e3*tEpoch,squeeze(mean(squeeze(processedSigHGTotal(:,chanInt,:)),2)))
xlabel('time (ms)')
ylabel('power (log(HG amplitude squared)')
xlim([-500 500])
vline(0)
title(['hilbert HG amplitude - channel ' num2str(chanInt)])
subplot(2,1,2)
plot(1e3*tMorlet,mean(squeeze(HGPowerWaveletTactorTotal(:,chanInt,:)),2))
xlim([-500 500])
vline(0)
xlabel('time (ms)')
ylabel('power normalized to baseline')
title(['average wavelet amplitude - channel ' num2str(chanInt)])
%%
figure
trials = 1:size(HGPowerWaveletTactorTotal,3);
time = tMorlet;
tLow = -0.2;
tHigh = 0.5;
imagesc(1e3*tMorlet(tMorlet>tLow & tMorlet < tHigh),trials,squeeze(HGPowerWaveletTactorTotal((tMorlet>tLow & tMorlet < tHigh),chanInt,:))')
colormap(flipud(bone))
axis('normal')
ylabel('trial')
xlabel('time (ms)')
colorbar()
title('average wavelet HG amplitude')
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
trainDuration = [];
modePlot = 'avg';
xlims = [-500 1000];
ylims = [-1 5];
HGPowerMeanTotal = mean(HGPowerWaveletTactorTotal(:,[1:64,65:72,85:92],:),3);
vizFunc.small_multiples_time_series(HGPowerMeanTotal/1e6,tMorlet,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot','ind','highlightRange',trainDuration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now for detecting interactions

processedSigTactorTotalSub = processedSigTactorTotal(:,[1:64,65:72,85:92],:);
%% phase slope index

freqbins = {[4:8],[8:12],[12:30],[70:150]};
%freqbins = {[70:150]};
%freqbins = {[8:12]};
% parameters for PSI-calculation
segleng=500;epleng=1000;
tBeginSub = -0.5; % seconds
tEndSub = 1; % seconds;

counter = 1;
clearvars psiTotal
for freq = freqbins
    freq = freq{:};
    
    for trial = 1:size(processedSigTactorTotalSub,3)
        [psi, stdpsi, psisum, stdpsisum]=data2psi(processedSigTactorTotalSub(tEpoch>=tBeginSub & tEpoch<=tEndSub,:,trial),segleng,epleng,freq);
        psi./(stdpsi+eps);
        psiTotal(:,:,trial,counter) = psi;
    end
    counter = counter + 1;
end

%%
counter = 1;
figure
for freq = freqbins
    freq = freq{:};
        subplot(2,2,counter)

    imagesc(mean(psiTotal(:,:,:,counter),3));
    xlabel('Channel')
    ylabel('Channel')
    set(gca,'fontsize',14)
    colorbar()
    title(['Phase slope index for somatosensory data for ' num2str(freq(1)) '-' num2str(freq(end)) ' Hz']);
    counter = counter + 1;
        colormap(gca,cmap)

end

%% phase slope index on average
freqbins = {[4:8],[8:12],[12:30],[70:150]};
%freqbins = {[70:150]};
%freqbins = {[8:12]};
% parameters for PSI-calculation
segleng=500;epleng=1000;
counter = 1;
clearvars psiTotalAverage
for freq = freqbins
    freq = freq{:};
    
    [psi, stdpsi, psisum, stdpsisum]=data2psi(mean(processedSigTactorTotalSub(tEpoch>=tBeginSub & tEpoch<=tEndSub,:,:),3),segleng,epleng,freq);
    psi./(stdpsi+eps);
    
    psiTotalAverage(:,:,counter) = psi;
    
    counter = counter + 1;
end

%%
counter = 1;
    figure

for freq = freqbins
    freq = freq{:};
    subplot(2,2,counter)
    imagesc(psiTotal(:,:,counter));
    xlabel('Channel')
    ylabel('Channel')
    set(gca,'fontsize',16)
    colorbar()
    title(['Phase slope index for somatosensory data for ' num2str(freq(1)) '-' num2str(freq(end)) ' Hz']);
    counter = counter + 1;
    colormap(gca,cmap)

end

%% correlation
subsetHGsignal = HGPowerMeanTotal(tMorlet>=tBeginSub & tMorlet <= tEndSub,:);
[c,lags] = xcorr(subsetHGsignal,'biased');
lags = 1e3*10*lags/fsData ;
cMat = reshape(c,length(lags),size(subsetHGsignal,2),[]); % now the second dimension has the cov between the channel in the 3rd dimension and all others

%% get normal differences
load('america');
cmap = cm;
CT = cm;


[cMax,ind] = max(abs(cMat),[],1);
lagsMax = lags(ind);

cMax = squeeze(cMax);
lagsMax = squeeze(lagsMax);

figure
subplot(2,1,1)
imagesc(cMax)
set(gca,'fontsize',16)
title('Maximum cross-correlation coefficient')
xlabel('Channel')
ylabel('Channel')
colormap(gca,pink)
colorbar()

subplot(2,1,2)
imagesc(lagsMax)
xlabel('Channel')
ylabel('Channel')
title('Lag (ms) at maximum cross-correlation coefficient')
set(gca,'fontsize',16)
colormap(gca,cmap)
colorbar()

%%

chanInt1 = 3;
chanInt2 = 22;
[cMaxTrial,indTrial] = max(abs(cMat(:,chanInt1,chanInt2)));
lagsMaxTrial = lags(indTrial);
figure
subplot(3,1,1)
plot(lags,cMat(:,chanInt1,chanInt2),'linewidth',3);
set(gca,'fontsize',16)
vline(lagsMaxTrial,'black')
title(['Cross-correlation between grid channels ' num2str(chanInt1) ' and ' num2str(chanInt2)])

chanInt1 = 3;
chanInt2 = 36;
[cMaxTrial,indTrial] = max(abs(cMat(:,chanInt1,chanInt2)));
lagsMaxTrial = lags(indTrial);
subplot(3,1,2)
plot(lags,cMat(:,chanInt1,chanInt2),'linewidth',3);
title(['Cross-correlation between grid channels ' num2str(chanInt1) ' and ' num2str(chanInt2)])
set(gca,'fontsize',16)
vline(lagsMaxTrial,'black')

chanInt1 = 3;
chanInt2 = 79;
[cMaxTrial,indTrial] = max(abs(cMat(:,chanInt1,chanInt2)));
lagsMaxTrial = lags(indTrial);
subplot(3,1,3)
plot(lags,cMat(:,chanInt1,chanInt2),'linewidth',3);
set(gca,'fontsize',16)
vline(lagsMaxTrial,'black')
title(['Cross-correlation between grid channel ' num2str(chanInt1) ' and RPT 6'])


set(gca,'fontsize',16)
vline(lagsMaxTrial,'black')
xlabel('time lag (ms)')
ylabel('correlation value')


return
%% significance with permutation test

nperms = 1000;
sizeSignal = size(subsetHGsignal);

cMaxShuffTotal = zeros([size(cMax),nperms]);
lagsMaxShuffTotal = zeros([size(cMax),nperms]);

for iterate = 1:nperms
    
    [~, out] = sort(rand(sizeSignal(1),sizeSignal(2)),1); % make random number matrix
    for chanInd = 1:size(subsetHGsignal,2)
        out(:,chanInd) = 1:size(subsetHGsignal,1);
        subsetHGsignalShuff = subsetHGsignal(out);
        [c,lags] = xcov(subsetHGsignalShuff,'biased');
        cMatShuff = reshape(c,length(lags),64,[]); % now the second dimension has the cov between the channel in the 3rd dimension and all others
        [cMaxShuff,ind] = max(abs(cMatShuff),[],1);
        lagsMaxShuff = lags(ind);
        cMaxShuff = squeeze(cMaxShuff);
        lagsMaxShuff = squeeze(lagsMaxShuff);
        
        cMaxShuffTotal(chanInd,:,iterate) = cMaxShuff(chanInd,:);
        lagsMaxShuffTotal(chanInd,:,iterate) = lagsMaxShuff(chanInd,:);
    end
    % cMatShuffTotal(:,:,:,ind) = cMatShuff;
end
%%



%%
