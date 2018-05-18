function [stimTimes,trainTimes] = extract_stimulation_times_TOJ_readIn(tact,fsStim,bads)

samplesOfPulse = round(2*fsStim/1e3);
% build a burst table with the timing of stimuli
preStim = round(fsStim*0.8);

% use stimDelay, which is 200 ms delayed 
[trainTimes] = find(tact(:,7)==1)- preStim;
logicalMask = logical(ones(size(trainTimes)));
logicalMask(bads) = 0;
trainTimes = trainTimes(logicalMask);
starts = trainTimes;
ends = trainTimes + samplesOfPulse;
delay = round(0.2867*fsStim/1e3);
% extract data
% try and account for delay for the stim times
stimTimes = starts+delay;
trainTimes=stimTimes;

end