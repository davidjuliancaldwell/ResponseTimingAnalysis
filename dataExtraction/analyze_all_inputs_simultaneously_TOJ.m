function [] = analyze_all_inputs_simultaneously_TOJ(tactorData,buttonData,stim,stimFromFile,fsTact)

tTact = (0:length(tactorData)-1)/fsTact;
figure
plot(tTact,tactorData);

title('tactor data')

% look at button press

tButton = (0:length(buttonData)-1)/fsTact;
figure
plot(tButton,buttonData);

title('button data')

% look at stim from file saved

tStimFile = (0:length(stim)-1)/fsTact;
figure
plot(tStimFile,stimFromFile);
title('stim from file')

% look at all 3 + stim waveform

figure
ax1 = subplot(5,1,1);
plot(tTact,tactorData)
title('tactor data')

ax2 = subplot(5,1,2);
plot(tButton,buttonData);
title('button data')

ax3 = subplot(5,1,3);
plot(tStimFile,stimFromFile);
title('stim from file')

% assuming stim1 here is the channel where stim was being delivered
ax4 = subplot(5,1,4);
plot(tStimFile,stim(:,1));
title('S1 Stim Channel')

%
ax5 = subplot(5,1,5);
plot(tStimFile,stim(:,2));
title('Off Target Stim Channel')

%link axis
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

end