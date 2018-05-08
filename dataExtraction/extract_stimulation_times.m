function [trainTimesTotal,stimFromFile,trainTimes,condType,uniqueCond] = extract_stimulation_times(tact,condType)

% stim cue from file
stimFromFile = tact(:,3);

% button press, start being onset of stimulation marker

trainTimes = find(stimFromFile~=0);

% pick condition type where stimulation was delivered

uniqueCond = unique(condType);
trainTimesTotal = {};

for i = 1:length(uniqueCond)
    trainTimesTotal{i} = trainTimes(condType==uniqueCond(i));
end


end