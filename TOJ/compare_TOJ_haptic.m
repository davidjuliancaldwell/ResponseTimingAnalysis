
loadFile = 1;
%%
if loadFile
    load('G:\My Drive\GRIDLabDavidShared\ResponseTiming\3ada8b_TOJ_8_17_2018_brainData_DCSonly_tactorOnly.mat')
end

figure
chanInt = 11;

plot(1e3*tEpoch,1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))
hold on
xlabel('time (ms)');
ylabel('microvolts')
title(['Processed Channel ' ]);


plot(1e3*tEpoch,1e6*avgResponseShift(:,chanInt))
plot(1e3*tEpoch,1e6*avgResponse(:,chanInt))
%%
if loadFile
    load('3ada8b_priming_neural_block_1_cond4_processed.mat')
end

%%
plot(1e3*tEpoch,1e6*mean(squeeze(processedSigReref(:,chanInt,:)),2))

ylabel('signal (\mu V)')

xlabel('time (ms)')
legend('Haptic only','TOJ aligned on tactor','TOJ aligned on stimulation train','DCS only')
title('Haptic only compared to Simultaneous Stimulation')
%%
ylims = [-(max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100) (max(abs(1e6*mean(squeeze(processedSigTactor(:,chanInt,:)),2))) + 100)];
ylim(ylims);

xlim([-500 2000]);
vline(mean(stimTime),'r','stim')

set(gca,'fontsize',14')