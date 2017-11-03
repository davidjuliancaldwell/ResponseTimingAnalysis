function evaluated_value = objective_func_ica(raw_sig,processedSig,varargin)

numTrials = size(raw_sig,3);
numChans = size(raw_sig,2);

meanSqError_mat = zeros(numChans,numTrials);
artifactSize_mat = zeros(numChans,numTrials);

% default params
fs = 1.2207e04;
plotIt = 0;

for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'plotIt'
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
    startInd = find(abs(zscore(diff_sig(:,chanMax)))>0.1,1,'first');
    endInd = find(abs(zscore(diff_sig(:,chanMax)))>0.1,1,'last');
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
    
    % mean square error for region outside of artifact
    meanSqError = mean((raw_sig(~art,:,ind) - processedSig(~art,:,ind)).^2);   % mean sq error for region outside of artifact
    
    
    % mean square error for artifact size
    artifactSize = std(processedSig(art,:,ind),[],1); % what is the standard deviation of the artifact - if too big, then there's lots of stim still there!
    
    % build up matrix trial by trial
    meanSqError_mat(:,ind) = meanSqError;
    artifactSize_mat(:,ind) = artifactSize;
    
    if plotIt
        figure
        hold on
        trial_inds = [1:numTrials];
        for ind = trial_inds
            scatter(repmat(ind,[size(meanSqError_mat,1),1]),meanSqError_mat(:,ind))
        end
        
        figure
        hold on
        trial_inds = [1:numTrials];
        for ind = trial_inds
            scatter(repmat(ind,[size(meanSqError_mat,1),1]),artifactSize_mat(:,ind))
        end
    end
    
end

% 
evaluated_value = 1e6*mean(mean(meanSqError_mat))+1e6*mean(mean(artifactSize_mat));


end