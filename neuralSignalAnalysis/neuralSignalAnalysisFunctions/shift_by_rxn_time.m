function [sigShifted,tShift] = shift_by_rxn_time(sig,buttonLocs,plotIt)

rxnTimes = buttonLocs{condInt};
[sorted,indexes] = sort(rxnTimes);
sortedSig = sig(:,:,indexes);

% calculate first response time, shift others based from there
sortedBasedOffFirst = sorted - sorted(1);
sigShifted = sortedSig;

numNonNan = sum(~isnan(rxnTimes));

for i = 2:numNonNan
    sigShifted(:,:,i) = circshift(sortedSig(:,:,i),sortedBasedOffFirst(i),1);
end

avgResponse = mean(sigShifted,3);
tShift = t - sorted(1)/fs_stim;

end