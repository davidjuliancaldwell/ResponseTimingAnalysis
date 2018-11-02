
loadFile = 1;
%%
if loadFile
    % load('G:\My Drive\GRIDLabDavidShared\ResponseTiming\3ada8b_TOJ_8_17_2018_brainData_DCSonly_tactorOnly.mat')
    load('G:\My Drive\GRIDLabDavidShared\ResponseTiming\3ada8b_priming_neural_block_1_cond4_processed.mat')
    %
    tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
    dataRefCond4 = poweroutCond4(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
    [normalizedDataCond4] = normalize_spectrogram(dataRefCond4,powerout);
    
    
    load('G:\My Drive\GRIDLabDavidShared\ResponseTiming\3ada8b_TOJ_HDBSCAN_tactorOnlyAsWell_11_1_2018.mat') %  TOJ with new processing
    
end
%%
chanIntList = [1,2,22,30,35];
for chanInt = chanIntList
    figure
    order = 3;
    framelen = 501;
    plot(1e3*tEpoch,1e6*sgolayfilt_complete(mean(squeeze(processedSigTactor(:,chanInt,:)),2),order,framelen),'linewidth',2)
    hold on
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Processed Channel ' num2str(chanInt)]);
    
    plot(1e3*tEpoch,1e6*sgolayfilt_complete(avgResponseShift(:,chanInt),order,framelen),'linewidth',2)
    plot(1e3*tEpoch,1e6*sgolayfilt_complete(avgResponse(:,chanInt),order,framelen),'linewidth',2)
    plot(1e3*tEpoch,1e6*sgolayfilt_complete(mean(squeeze(processedSigCond4(:,chanInt,:)),2),order,framelen),'linewidth',2)
    
    %ylims = [-(max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100) (max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100)];
    ylims = [-60 60];
    ylim(ylims);
    
    xlim([-200 700]);
    vline(mean(stimTime),'r','stim')
    
    set(gca,'fontsize',20)
    ylabel('signal (\mu V)')
    
    xlabel('time (ms)')
    legend('Haptic only','TOJ aligned on tactor','TOJ aligned on stimulation train','DCS only')
    title('Haptic only compared to Simultaneous Stimulation')
end


%%
% chanIntList = chanInt;
stimTime = 0;
response = 0;

preStim = 1000;
postSTim = 2000;
tMorletTactor = linspace(-preStim,postStim,size(normalizedDataTactor,2))/1e3;

for chanInt = chanIntList
    visualize_wavelet_channel_onlyProcessed(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,chanInt,stimTime,response,individual,average)
    
    visualize_wavelet_channel_onlyProcessed(normalizedDataShift,tMorlet,fMorlet,sigShifted,...
        tEpoch,chanInt,stimTime,response,individual,average)
    
        visualize_wavelet_channel_onlyProcessed(normalizedDataCond4,tMorlet,fMorlet,sigShifted,...
        tEpoch,chanInt,stimTime,response,individual,average)
    
    visualize_wavelet_channel_onlyProcessed(normalizedDataTactor,tMorletTactor,fMorlet,processedSigTactor,...
        tEpoch,chanInt,stimTime,response,individual,average)
end

return
%%
% if loadFile
%     load('3ada8b_priming_neural_block_1_cond4_processed.mat')
% end


%%

%ylims = [-(max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100) (max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100)];
ylims = [-60 60]
ylim(ylims);

xlim([-200 700]);
vline(mean(stimTime),'r','stim')

set(gca,'fontsize',20)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
dataRefTactor = poweroutTactor(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
[normalizedData] = normalize_spectrogram(dataRef,powerout);
[normalizeDataTactor] = normalize_spectrogram(dataRefTactor,poweroutTactor);

%%
% chanIntList = chanInt;
chanIntList = 22;
for chanInt = chanIntList
    visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,epochedECoG,chanInt,stimTime,response,individual,average)
end
%%
for chanInt = chanIntList
    visualize_wavelet_channel(normalizedDataTactor,tMorlet,fMorlet,processedSig,...
        tEpoch,epochedECoG,chanInt,stimTime,response,individual,average)
end


%%
vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);

%%
vizFunc.small_multiples_spectrogram(normalizedDataTactor,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
