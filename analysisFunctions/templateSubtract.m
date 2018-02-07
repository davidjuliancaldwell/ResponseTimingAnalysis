function [processedSig,template] = templateSubtract(raw_sig,varargin)
%USAGE:
% This function will perform a template subtraction scheme for artifacts on
% a trial by trial, channel by channel basis
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

% default parameters
fs = 1.2207e04;
plotIt = 0;
type = 'linear';
pre = 0.4096;
post = 0.4096;

% define matrix of zeros
processedSig = zeros(size(raw_sig));

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
    
    inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>0.5);
    diff_bt_inds = [diff(inds)' 0];
    [~,inds_onset] = find(abs(zscore(diff_bt_inds))>0.5);
    start_inds = [inds(1)-presamps; inds(inds_onset+1)-presamps];
    end_inds = [ inds(inds_onset)+postsamps; inds(end)+postsamps];
    if plotIt
        figure
        plot(diff_sig(:,chanMax,trial))
        vline(start_inds)
        vline(end_inds,'g')
    end
    
    for chan = 1:size(raw_sig,2)
        % get single trial
        raw_sig_temp = raw_sig(:,chan,trial);
        % default time to average
        default_time_average = end_inds(1)-start_inds(1)+1;
        
        avg_signal = zeros((end_inds(1)-start_inds(1)+1),length(start_inds));
        for sts = 1:length(start_inds)
            win = start_inds(sts):start_inds(sts)+default_time_average-1;
           % avg_signal(:,sts) = raw_sig_temp(win);
            avg_signal(:,sts) = raw_sig_temp(win) - mean(raw_sig_temp(start_inds(sts):start_inds(sts)+presamps-5));% take off pre period

        end
        
        % find average stimulus
        avg_signal_mean = mean(avg_signal,2);
        
        for sts = 1:length(start_inds)
            win = start_inds(sts):start_inds(sts)+default_time_average-1;
            raw_sig_temp(win) = raw_sig_temp(win) - avg_signal_mean;
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            hold on
            plot(raw_sig(:,chan,trial),'linewidth',2)
            
            vline(start_inds)
            vline(end_inds,'g')
        end
        
        processedSig(:,chan,trial) = raw_sig_temp;
        template{chan}{trial} = avg_signal;
    end
    
end


end