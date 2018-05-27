%% 5.12.2018 - David J. Caldwell
% analyze different response timing subjects one at a time

close all; clearvars ; clc
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks

% sid
% 1 - acabb1
% 2 - c19968
% 3 - 693ffd
% 4 - 2fd831
% 5 - a1355e
SIDSint = {'c19968','693ffd','2fd831'};
SIDSint = {'693ffd'};
%%
for i = SIDSint
    %%
    sid = i{:};
    DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';
    load(fullfile([sid 'pooledData_tactorSub.mat']));
    fsData = fs_data;
    %% combine the pooled data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [buttonLocsSamps,buttonLocs,data,tEpoch,uniqueCond] = combine_pooled_data(sid,epochedCortEco_cell,t_epoch);
    
    % additional parameters
    postStim = 2000;
    sampsPostStim = round(postStim/1e3*fs_data);
    
    preStim = 1000;
    sampsPreStim = round(preStim/1e3*fs_data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    switch sid
        case 'c19968'
            
            stimChans = [9 17 50 58];
            chanIntList = [1 10 51 42];
            chanInt = 10;
            
        case '693ffd'
            chanInt = 17;
            stimChans = [20 29];
            chanIntList = [21 28 19 36 44 43 30];
        case '2fd831'
            stimChans = [1 9 24 42];
            chanInt = 10;
            chanIntList = [2 10 51 42];
            
    end
    
    trainDuration = [];
    xlims = [-50 500];
    minDuration = 0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % params for DBSscan optimization
    
    type = 'dictionary';
    
    useFixedEnd = 0;
    %fixedDistance = 2;
    fixedDistance = 4; % in ms
    plotIt = 0;
    
    %pre = 0.4096; % in ms
    %post = 0.4096; % in ms
    
    pre = 0.8; % started with 1
    post = 0.6; % started with 0.2
    % 2.8, 1, 0.5 was 3/19/2018
    
    % these are the metrics used if the dictionary method is selected. The
    % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
    % cosine similarity, or correlation for clustering and template matching.
    
    distanceMetricDbscan = 'cosine';
    distanceMetricSigMatch = 'eucl';
    amntPreAverage = 3;
    normalize = 'preAverage';
    %normalize = 'firstSamp';
    
    recoverExp = 0;
    %%
    % condIntAns
    % -1 = tactor
    %  0 = null
    %  1 = off-target
    %  2-5 = conditions 1->4
    
    condInt = 1;
    condIntAns = uniqueCond(condInt);
    dataInt = data{condInt};
    response = buttonLocs{condInt};
    
    % get additional fake channels on 693ffd so it plots ok
    if strcmp(sid,'693ffd')
        dataInt = cat(2,dataInt,zeros(size(dataInt,1),64-size(dataInt,2),size(dataInt,3)));
    end
    
    buttonLocsInt = buttonLocs{condInt};
    %%%%%%%%%%%%%%%%%%
    if (condIntAns == 2 || condIntAns == 3 || condIntAns == 4 || condIntAns == 5)
        meanSub = 1;
        %orderPoly = 6;
        orderPoly = 3; %10-12-2017 - djc change
        if meanSub == 1
            for i = 1:size(dataInt,2)
                dataInt(:,i,:) = polyfit_subtract(squeeze(dataInt(:,i,:)),orderPoly);
            end
        end
        
        
        [processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(dataInt,'type',type,...
            'fs',fsData,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
            'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
            'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,'minDuration',minDuration);
        stimTime = zeros(size(processedSig,3),1); % it is centered around zero now
        
    elseif (condIntAns == -1)
        meanSub = 1;
        %orderPoly = 6;
        orderPoly = 3; %10-12-2017 - djc change
        if meanSub == 1
            for i = 1:size(dataInt,2)
                processedSig(:,i,:) = polyfit_subtract(squeeze(dataInt(:,i,:)),orderPoly);
            end
        else
            processedSig = dataInt;
        end
        
        %stimTime = 1e3*tactorLocsVec; %
        stimTime = zeros(size(processedSig,3),1); % it is centered around zero now
        t = t_epoch;
    end
    
    %%
    % notch 
%     for trial = 1:size(processedSig,3)
%     processedSig(:,:,trial) = notch(squeeze(processedSig(:,:,trial)),[60 120 180 240],fsData); 
%     end
    %%
    % visualization
    % of note - more visualizations are created here, including what the
    % templates look like on each channel, and what the discovered templates are
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trainDuration = [];
    vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
        tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
        'chanIntList',chanIntList,'modePlot','confInt')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% wavelet and plv
    
    %%
    %%%%%% PLV
    freqRange = [8 12];
    [plv] = plvWrapper(processedSig,fsData,freqRange,stimChans);
    %%
    %%%%%%% wavelet
    timeRes = 0.01; % 25 ms bins
    
    % [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
    [powerout,fMorlet,tMorlet,~] = waveletWrapper(processedSig,fsData,timeRes,stimChans);
    %
    tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
    % normalize data
    dataRef = powerout(:,tMorlet<0,:,:);
    %
    [normalizedData] = normalize_spectrogram(dataRef,powerout);
    %%
    chanIntList = [19 28];
    individual = 1;
    average = 1;
    for chanInt = chanIntList
        visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
            tEpoch,dataInt,chanInt,stimTime,response,individual,average)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
    %% hilb amp HG
    return
    processedSigHG = zeros(size(processedSig));
    for trial = 1:size(processedSig,3)
        [amp] = log(hilbAmp(squeeze(processedSig(:,:,trial)), [70 110], fsData).^2);
        processedSigHG(:,:,trial) = amp;
    end
    
vizFunc.small_multiples_time_series(processedSigHG,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',[-40 -20],'modePlot','avg','highlightRange',trainDuration)

    
    % sort by rxn time
    %%
    sortByRxnTime = 0;
    if sortByRxnTime
        [sigShifted,tShift] = shift_by_rxn_time(processedSig,buttonLocs)
        smallMultiples_responseTiming(avgResponse,tShift,'type1',stimChans,'type2',0,'average',1)
    end
end