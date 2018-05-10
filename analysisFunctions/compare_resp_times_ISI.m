function [mdl,mdlNoNuOt] = compare_resp_times_ISI(uniqueCond,buttonLocs,ISICellSecondsNoNuOt,ISICellSeconds)

% if want to change reaction times, do it here!
respLo = 0.100;
respHi = 1;


buttonLocsAllStim = [buttonLocs{3:7}];
ISICellSecondsNoNuOtAllStim = [ISICellSecondsNoNuOt{3:7}];
ISICellSecondsAllStim = [ISICellSeconds{3:7}];

%%
% individual scatters of RT
% make cell array for legends
uniqueCondText = cellstr(num2str(uniqueCond));
uniqueCondText{1} = 'tactor';
uniqueCondText{2} = 'no stimulation';
uniqueCondText{3} = 'off target stimulation';
uniqueCondText{4} = '100 ms train';
uniqueCondText{5} = '200 ms train';
uniqueCondText{6} = '400 ms train';
uniqueCondText{7} = '800 ms train';
%individual  histogram of each condition type
figure

for i = 1:length(uniqueCond)
    subplot(length(uniqueCond),1,i)
    a = scatter(1e3*ISICellSeconds{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi));
    title(['Scatter of RT vs. ISI ' uniqueCondText{i} ])
    xlabel('ISI (ms)')
    ylabel('RT')
    ylim([0 1000])
end

%% overall STIM rt

figure
scatter(1e3*ISICellSecondsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),...
    1e3*buttonLocsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi));
ylim([0 1000])

figure
mdl = fitlm(1e3*ISICellSecondsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),...
    1e3*buttonLocsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),'intercept',true);
hold on
plot(mdl)
ylim([0 1000])
title(['Scatter of RT vs. ISI ' ])
xlabel('ISI (ms)')
ylabel('RT (ms)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% individual scatters of RT
% make cell array for legends
uniqueCondText = cellstr(num2str(uniqueCond));
uniqueCondText{1} = 'tactor';
uniqueCondText{2} = 'no stimulation';
uniqueCondText{3} = 'off target stimulation';
uniqueCondText{4} = '100 ms train';
uniqueCondText{5} = '200 ms train';
uniqueCondText{6} = '400 ms train';
uniqueCondText{7} = '800 ms train';
%individual  histogram of each condition type
figure

for i = 1:length(uniqueCond)
    subplot(length(uniqueCond),1,i)
    a = scatter(1e3*ISICellSecondsNoNuOt{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi));
    title(['Scatter of RT vs. ISI ' uniqueCondText{i} ])
    xlabel('ISI (ms)')
    ylabel('RT')
    ylim([0 1000])
end

%% overall STIM rt

figure
scatter(1e3*ISICellSecondsNoNuOtAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),...
    1e3*buttonLocsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi));
ylim([0 1000])

figure
mdlNoNuOt = fitlm(1e3*ISICellSecondsNoNuOtAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),...
    1e3*buttonLocsAllStim(buttonLocsAllStim>respLo & buttonLocsAllStim<respHi),'intercept',true);
hold on
plot(mdlNoNuOt)
ylim([0 1000])

title(['Scatter of RT vs. ISI ' ])
xlabel('ISI (ms)')
ylabel('RT (ms)')

end