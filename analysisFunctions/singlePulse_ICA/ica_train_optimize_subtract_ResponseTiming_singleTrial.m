function [evaluated_value,huberInside_mat,huberOutside_mat,processedSig,recon_artifact_matrix] =  ica_train_optimize_subtract_ResponseTiming_singleTrial(data,stimChans,fs_data,numComponentsSearch,scale_factor,meanSub,orderPoly,plotIt,channelInt,nonlinear)
%USAGE: function [subtracted_sig_matrixS_I, subtracted_sig_cellS_I] =  ica_artifact_remove(t,data,stimChans,pre,post,fs_data,scale_factor,numComponentsSearch,plotIt,channelInt)
%This function will perform the fast_ica algorithm upon a data set in the
%format of m x n x p, where m is samples, n is channels, and p is the
%individual trial. This is for trains of stimuli
%
% data = samples x channels x trials
% tTotal =  time vector
% stimChans = stimulation channels, or any channels to ignore
% pre = the time point at which to begin extracting the signal
% post = the time point at which to stop extracting the signal
% fs_data = sampling rate (Hz)
% scale_factor = scaling factor tp ensure the ICA algorithm functions
%       correctly
%numComponentsSearch = the number of ICA components to search through for
%       artifacts that meet a certain profile
% plotIt = plot it or not
% channelInt = plot a channel if interested
% REQUIRES FastICA algorithm in path

% set scale factor


% plot intermediate steps

if (~exist('nonlinear','var'))
    nonlinear = 'gauss';
end

if (~exist('meanSub','var'))
    meanSub = 0;
end

if (~exist('orderPoly','var'))
    orderPoly = 2;
end

% make a time vector if one doesn't exist
if (~exist('tTotal','var'))
    tTotal = 0:size(data,1)-1;
end

% make a pre time condition to start from
% if this is not input, matching of artifact will fail
if (~exist('pre','var'))
    pre = tTotal(1);
end

% make a post time condition to start from
% if this is not input, matching of artifact will fail
if (~exist('post','var'))
    post = tTotal(end);
end

freqFilter = false;

numComponentsSearch = round(numComponentsSearch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stim channels, as we don't want to perform ICA on them

bads = [];
badTotal = [stimChans; bads];

% total channels
numChans = size(data,2);
% make logical good channels matrix to index
goods = zeros(numChans,1);

channelsOfInt = 1:numChans;

goods(channelsOfInt) = 1;
% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;
% make it logical
goods = logical(goods);

% make storage matrices
i_icasigS = {};
i_mixing_matS = {};
i_sep_matS = {};

% extract the data of interest
dataInt = data(:,goods,:);

% NOTE THIS IS DIFFERENT THAN BEFORE, WE WANT TO KEEP STIMULATION IN THERE
dataIntTime = dataInt((tTotal>=pre & tTotal<=post),:,:);

t = tTotal(tTotal>=pre & tTotal<=post); % get new subselected t vector

if meanSub == 1
    for i = 1:size(dataIntTime,2)
        for j = 1:size(dataIntTime,3)
            data_int_temp = squeeze(dataIntTime(:,i,j));
            [p,s,mu] = polyfit((1:numel(data_int_temp))',data_int_temp,orderPoly);
            f_y = polyval(p,(1:numel(data_int_temp))',[],mu);
            
            % subtract poly fit
            dataIntTime(:,i,j) = data_int_temp - f_y;
            
            %dataIntTime = dataIntTime - repmat(mean(data,1),size(data,1),1);
            
        end
        %         figure;
        %         plot(f_y)
    end
end

numTrials = size(dataIntTime,3);

length_signal = size(dataIntTime,1);

ratio_mat = zeros(sum(goods),numTrials);

rescale = false;

for i = 1:numTrials
    trial_dataIntTime = squeeze(dataIntTime(:,:,i));
    
    if ~rescale % not rescale
        
        sig_epoch = scale_factor.*trial_dataIntTime;
    else
        % rescale each to max - bing's suggestion
        sig_max = max(abs(dataIntTime(:,:,i)));
        ratio = 1./(sig_max/max(sig_max));
        ratio_mat(:,i) = ratio;
        sig_epoch = scale_factor.*trial_dataIntTime.*repmat(ratio,length_signal,1);
    end
    
    %[icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch','g',nonlinear,'approach',...
    %    'defl','stabilization','on','numOfIC',numComponentsSearch,'verbose','off');
    %[icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch','g','gauss','approach','symm','numOfIC',numComponentsSearch);
    [icasig_temp,mixing_mat_temp,sep_mat_temp] = fastica(sig_epoch','g',nonlinear,'approach',...
        'defl','stabilization','on','numOfIC',numComponentsSearch,'verbose','off','epsilon',0.00001,'maxNumIterations',5000,'maxFinetune',300);    
    
    i_icasigS{i} = icasig_temp;
    i_mixing_matS{i} = mixing_mat_temp;
    i_sep_matS{i} = sep_mat_temp;
    %fprintf(['iteration ' num2str(i) '\n'])
end

% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'maxFinetune'         (integer) Maximum number of iterations in 
%                       fine-tuning. Default 100.
%
% 'sampleSize'          (number) [0 - 1] Percentage of samples used in
%                       one iteration. Samples are chosen in random.
%                       Default is 1 (all samples).
%
    
%     numInt = min(size(icasig_temp,1),10);
% 
%     for j = 1:size(dataIntTime,3)
%         figure
%         for i = 1:numInt
%             subplot(numInt,1,i)
%             plot(i_icasigS{j}(i,:))
%             title(['ICA component # , ', num2str(i)])
%         end
%         subtitle(['Trial # ', num2str(j)])
% 
%     end

%% set ICA components that are like the artifact to zero (they occur near a certain time and have prominence)

% need to adjust this for case where it's close to zero but not quite
% equal?
numTrials = size(dataIntTime,3);


i_ica_kept = {};
i_ica_mix_kept = {};

for i = 1:numTrials
    numICs = size(i_icasigS{i},1);
    start_index = 1;
    for j = 1:numICs
        
        % have to tune this
        [pk_temp_pos,locs_temp_pos] = findpeaks(zscore(i_icasigS{i}(j,:)),fs_data,'MinPeakProminence',5);
        [pk_temp_neg,locs_temp_neg] = findpeaks(-1*zscore(i_icasigS{i}(j,:)),fs_data,'MinPeakProminence',5);
        %
        
        
        % DJC 4-17-2017 - add in 200 Hz stim peak
        if ((~isempty(pk_temp_pos) || ~isempty(pk_temp_neg)))
            
            i_ica_kept{i}(start_index,:) = i_icasigS{i}(j,:);
            i_ica_mix_kept{i}(:,start_index) = i_mixing_matS{i}(:,j);
            
            start_index = start_index+1;
        end
        
    end
    
    if isempty(i_ica_kept) || isempty(i_ica_kept{i})
        i_ica_kept{i} = zeros(size(i_icasigS{i}(j,:)));
        i_ica_mix_kept{i} = zeros(size(i_mixing_matS{i}(:,j)));
    end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct stim artifact across channels

% make matrix of reconstruction artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_art = zeros(size(dataIntTime,1),size(data,2));

for i = 1:numTrials
    
    recon_artifact_temp = (i_ica_mix_kept{i}*i_ica_kept{i})'./scale_factor;
    
    total_art(:,goods) = recon_artifact_temp;
    total_art(:,badTotal) = zeros(size(recon_artifact_temp,1),size(badTotal,2));
    
    recon_artifact{i} = total_art;
    recon_artifact_matrix(:,:,i) = total_art;
    num_modes_kept = size(i_ica_kept{i},1);
    
    
end

%% subtract each one of these components

processedSig_cell = {};

processedSig = zeros(size(dataIntTime,1),size(data,2),size(numTrials,1));

total_sig = zeros(size(dataIntTime,1),size(data,2));

for i = 1:numTrials
    
    
    combined_ica_recon = (i_ica_mix_kept{i}*i_ica_kept{i})';
    
    num_modes_kept = size(i_ica_kept{i},1);
    
    % unscaled
    if ~rescale
        subtracted_sig_ICA_temp = dataIntTime(:,:,i) - combined_ica_recon./scale_factor;
    else
        % rescaled
        
        subtracted_sig_ICA_temp = dataIntTime(:,:,i) - ((combined_ica_recon./scale_factor)./repmat(ratio_mat(:,i)',length_signal,1));
    end
    % add in bad channels back
    total_sig(:,goods) = subtracted_sig_ICA_temp;
    total_sig(:,badTotal) = zeros(size(subtracted_sig_ICA_temp,1),size(badTotal,2));
    
    
    processedSig_cell{i} = total_sig;
    processedSig(:,:,i) = total_sig;
    
    
    if plotIt
        figure
        plot(t,1e6*total_sig(:,channelInt),'LineWidth',2)
        hold on
        plot(t,1e6*data((tTotal>=pre & tTotal<=post),channelInt,i),'LineWidth',2)
        title(['Channel ', num2str(channelInt), ' Trial ', num2str(i), ' Number of ICA modes subtracted = ', num2str(num_modes_kept)])
        legend({'subtracted signal','original signal'})
        ylabel(['Signal \muV'])
        xlabel(['Time (ms)'])
        set(gca,'Fontsize',[14]),
        
        figure
        plot(t,1e6*total_sig(:,channelInt),'LineWidth',2)
        title(['Subtracted Signal for ', num2str(num_modes_kept), ' ICA modes, Channel ', num2str(channelInt), ' Trial ', num2str(i)])
        ylabel(['Signal \muV'])
        xlabel(['Time (ms)'])
        set(gca,'Fontsize',[14])
    end
    
    %
end
[evaluated_value,huberInside_mat,huberOutside_mat] = objective_func_ica_huber(data(:,goods,:),processedSig(:,goods,:),recon_artifact_matrix(:,goods,:),'fs',fs_data,'plotIt',0); % try using huber loss function
evaluated_value = evaluated_value*1e6;
end
