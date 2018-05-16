function [bursts,delay] = extract_stimulus_delivery(stim,sing,condType,trainTimes,fsStim,fsSing,plotIt)

% plot stim

if plotIt
    figure
    hold on
    for i = 1:size(stim,2)
        
        t = (0:length(stim)-1)/fsStim;
        subplot(3,2,i)
        plot(t*1e3,stim(:,i))
        title(sprintf('Channel %d',i))
        
    end
    
    xlabel('Time (ms)')
    ylabel('Amplitude (V)')
    subtitle('Stimulation Channels')
    
end

% 1st stim channel
Sing1 = sing(:,1);
% 2nd stim channel
Sing2 = sing(:,4);

samplesOfPulse = round(2*fsStim/1e3);

% build a burst table with the timing of stimuli
bursts = [];
bursts(1,:) = condType;
bursts(2,:) = trainTimes;
bursts(3,:) = trainTimes + samplesOfPulse;

stims1 = squeeze(getEpochSignal(Sing1,(bursts(2,condType>1)-1),(bursts(3,condType>1))+1));
t = (0:size(stims1,1)-1)/fsSing;
t = t*1e3;

stims2 = squeeze(getEpochSignal(Sing2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stims1,1)-1)/fsSing;
t = t*1e3;

if plotIt
    
    figure
    plot(t,stims1,'b','linewidth',2)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (mA)')
    ylim([(min(stims1(:))-100) (max(stims1(:))+100)])
    title('Current to be delivered for all trials  on 1st channel')
    
    figure
    plot(t,stims2,'b','linewidth',2)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (mA)')
    ylim([(min(stims2(:))-100) (max(stims2(:))+100)])
    title('Current to be delivered for all trials  on 2nd channel')
    
end

%% Plot stims with info from above

% 1st stimulation channel
stim1 = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,condType>1)-1),(bursts(3,condType>1))+1));
t = (0:size(stim1Epoched,1)-1)/fsStim;
t = t*1e3;

% 2nd stimulation channel
stim2= stim(:,4);
stim2Epoched = squeeze(getEpochSignal(stim2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stim2Epoched,1)-1)/fsStim;
t = t*1e3;

if plotIt
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery - 1st stim channel')
    
    % hold on
    plot(t,stims1)
    
    figure
    plot(t,stim2Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery - 2nd stim channel')
    
end

% get the delay in stim times

delay = round(0.2867*fsStim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType>1)-1+delay;
stimTimesEnd = bursts(3,condType>1)-1+delay;
stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
t = (0:size(stim1Epoched,1)-1)/fsStim;
t = t*1e3;

stim2Epoched = squeeze(getEpochSignal(stim2,stimTimesBegin,stimTimesEnd));


if plotIt
    
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in - 1st stim channel ')
    
    
    % hold on
    figure
    plot(t,stim2Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in - 2nd stim channel')
    
end


end