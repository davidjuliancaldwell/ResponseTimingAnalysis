function [evaluated_value,huberInside_mat,huberOutside_mat] = objective_func_ica_huber(raw_sig,processedSig,artifactSig,varargin)

numTrials = size(raw_sig,3);
numChans = size(raw_sig,2);

huberInside_mat = zeros(numChans,numTrials);
huberOutside_mat = zeros(numChans,numTrials);

% default params
fs = 1.2207e04;
plotIt = 0;

for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'plotit'
            plotIt = varargin{i+1};
        case 'fs'
            fs = varargin{i+1};
    end
end


for ind  = 1:numTrials
    
    % take diff of sig to find onset of stim train
    diff_sig = [zeros(size(raw_sig,2),1) diff(raw_sig(:,:,ind))']';
    
    % find channel that has the max signal
    [~,chanMax] = max(max(diff_sig));
    
    % find onset
    startInd = find(abs(zscore(diff_sig(:,chanMax)))>3,1,'first');
    endInd = find(abs(zscore(diff_sig(:,chanMax)))>3,1,'last');
    if plotIt
        figure
        hold on
        plot(diff_sig(:,chanMax))
        vline(startInd)
        vline(endInd)
    end
    
    % figure out when artifact is
    art = zeros(size(diff_sig,1),1);
    art(startInd:endInd) = 1;
    art = logical(art);
    
    % huber loss for region outside of artifact
    sigma = 1;
    huberOutsideArtifact = huber_loss(raw_sig(~art,:,ind),processedSig(~art,:,ind),sigma);
    
    % huber loss for region inside artifact
    sigma = 1;

    huberInsideArtifact = huber_loss(raw_sig(art,:,ind),artifactSig(art,:,ind),sigma);
    % build up matrix trial by trial
    huberInside_mat(:,ind) = huberInsideArtifact;
    huberOutside_mat(:,ind) = huberOutsideArtifact;
    
end

if plotIt
    figure
    hold on
    trial_inds = [1:numTrials];
    for ind = trial_inds
        scatter(repmat(ind,[size(huberOutside_mat,1),1]),huberOutside_mat(:,ind))
    end
    
    figure
    hold on
    trial_inds = [1:numTrials];
    for ind = trial_inds
        scatter(repmat(ind,[size(huberInside_mat,1),1]),huberInside_mat(:,ind))
    end
end
%
weight_inside = 1;
weight_outside = 1*weight_inside;

evaluated_value = weight_outside*mean(mean(huberOutsideArtifact))+weight_inside*mean(mean(huberInsideArtifact));

% power_vals = zeros;
% [f,p1] = spectralAnalysisComp(fs,processedSig);
% figure
% plot(f,log(p1))
% 
% [f,p1] = spectralAnalysisComp(fs,raw_sig);
% figure
% plot(f,log(p1))
% 
% 
% [f,p1] = spectralAnalysisComp(fs,raw_sig(~art,:,ind));
% figure
% plot(f,log(p1))
% 
% [f,p1] = spectralAnalysisComp(fs,processedSig(~art,:,ind));
% figure
% plot(f,log(p1))
% 
% [f,p1] = spectralAnalysisComp(fs,raw_sig(~art,:,ind));
% figure
% plot(f,log(p1))
% 
% [f,p1] = spectralAnalysisComp(fs,processedSig(~art,:,ind));
% figure
% plot(f,log(p1))

end