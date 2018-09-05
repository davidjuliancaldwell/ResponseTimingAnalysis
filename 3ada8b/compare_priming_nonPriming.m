
loadFile = 1;
%%
if loadFile
    load('3ada8b_priming_neural_block_1_cond7_processed.mat')
end


processedSig7 = processedSig;

%%
if loadFile
    load('3ada8b_priming_neural_block_1_cond4_processed.mat')
end

processedSig4 = processedSig;


%%
%chanIntList = 3;
trainDuration = [];
rerefMode =  'selectedChannels';
xlims = [-100 500];
ylims = [-400 400];
badChannels = stimChans;
channelsToUse = [5:8,13:23,25:31,33:64]';
processedSigReref4 = rereference_CAR_median(processedSig4,rerefMode,badChannels,[],[],channelsToUse);
processedSigReref7 = rereference_CAR_median(processedSig7,rerefMode,badChannels,[],[],channelsToUse);
%%
individual = 0;
average = 1;
modePlot = 'avg';
vizFunc.small_multiples_time_series(processedSig7,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
vizFunc.small_multiples_time_series(processedSigReref4,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)


