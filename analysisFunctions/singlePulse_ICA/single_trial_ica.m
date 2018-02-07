function [processedSig,reconArtifact] = single_trial_ica(raw_sig,varargin)
%USAGE:
% This function will perform an interpolation scheme for artifacts on a
% trial by trial, channel by channel basis, implementing either a linear
% interpolation scheme, or a pchip interpolation scheme
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% pre_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% post_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

% default parameters
fs = 1.2207e04;
plotIt = 0;
type = 'linear';
pre = 1;
post = 1;

nonlinear = 'gauss';
meanSub = 1;
orderPoly = 1;
tTotal = 0:size(raw_sig,1)-1;
freqFilter = false;
numComponentsSearch = 64;
scale_factor = 100;
numComponentsSearch = round(numComponentsSearch);
stimChans = [];
% define matrix of zeros
processedSig = zeros(size(raw_sig));
reconArtifact = zeros(size(raw_sig));

channelInt = 10;

for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'plotit'
            plotIt = varargin{i+1};
        case 'fs'
            fs = varargin{i+1};
        case 'type'
            type = varargin{i+1};
        case 'pre'
            pre = varargin{i+1};
        case 'post'
            post = varargin{i+1};
        case 'stimchans'
            stimChans = varargin{i+1};
            
    end
end

%% preprocess eco
presamps = round(pre/1e3 * fs); % pre time in sec

postsamps = round(post/1e3 * fs); %

% take diff of signal to find onset of stimulation train
diff_sig = permute(cat(3,zeros(size(raw_sig,2), size(raw_sig,3)),permute(diff(raw_sig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diff_sig)));
chanMax = chanMax(1);

for trial = 1:size(raw_sig,3)
    
    % find onset
    raw_sig_temp = raw_sig(:,:,trial);
    artifact_temp = reconArtifact(:,:,trial);
    inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>0.5);
    diff_bt_inds = [diff(inds)' 0];
    [~,inds_onset] = find(abs(zscore(diff_bt_inds))>0.5);
    start_inds = [inds(1)-presamps; inds(inds_onset+1)-presamps];
    end_inds = [ inds(inds_onset)+postsamps; inds(end)+postsamps];
    
    for sts = 1:length(start_inds)
        win = start_inds(sts):end_inds(sts);
        raw_sig_temp_toProcess = raw_sig(win,:,trial);
        [evaluated_value,huberInside_mat,huberOutside_mat,raw_sig_temp_clean,recon_artifact_matrix] =  ica_train_optimize_subtract_ResponseTiming_singleTrial(raw_sig_temp_toProcess,stimChans,fs,numComponentsSearch,scale_factor,meanSub,orderPoly,plotIt,channelInt,nonlinear);
        raw_sig_temp(win,:) = raw_sig_temp_clean;
        artifact_temp(win,:) = recon_artifact_matrix;
    end
    
    if plotIt
        figure
        plot(raw_sig_temp,'linewidth',2)
        hold on
        plot(raw_sig(:,channelInt,trial),'linewidth',2)
        
        vline(start_inds)
        vline(end_inds,'g')
    end
    processedSig(:,:,trial) = raw_sig_temp;
    reconArtifact = artifact_temp;
    fprintf(['iteration ' num2str(trial) ' complete \n'])

end

end


