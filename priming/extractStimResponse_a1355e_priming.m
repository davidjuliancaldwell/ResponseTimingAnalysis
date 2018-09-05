%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming';
sid = SIDS{5};

plotIt = 1;

% % ui box for input
% list_str = {'1st block','2nd block'};
%
% [s,v] = listdlg('PromptString','Pick experiment',...
%     'SelectionMode','single',...
%     'ListString',list_str);

for s = 2:2
    
    % load in data
    if (strcmp(sid, 'a1355e'))
        folder_data = strcat(DATA_DIR,'\2fd831');
        
        if s == 1
            load(fullfile(DATA_DIR,'ResponseTiming-3.mat'))
            block = '3';
        elseif s == 2
            load(fullfile(DATA_DIR,'ResponseTiming-4.mat'))
            block = '4';
        end
        
    end
    
    %% load in data of interest
    
    [stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);
    
    clear Stim Tact Sing
    
    
    %% figure out stim times
    
    condType = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_condition_primingPilot.txt');
    primedOption = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_primedOption_primingPilot.txt');
    train = dlmread('C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming\rxnTime_stimTrainDelivery_primingPilot.txt');
    
    [trainTimesTotal,stimFromFile,trainTimes,primedOption,uniqueCond] = extract_stimulation_times(tact,primedOption);
    
    %% extract stimulus data, find delay, and get timing of stimuli
    
    [bursts,delay] = extract_stimulus_delivery_primed(stim,sing,condType,primedOption,trainTimes,fsStim,fsSing,plotIt);
    
    %% extract data
    % try and account for delay for the stim times
    stimTimes = bursts(2,:)+delay;
    trainTimes=stimTimes;
    
    %% look at all simultaneously
    
    tactorData = tact(:,1);
    buttonData = tact(:,2);
    
    analyze_all_inputs_simultaneously(tactorData,buttonData,stim,stimFromFile,fsTact)
    %%
    [buttonLocs,buttonLocsSamps,~,~,tEpoch,epochedButton] = get_response_timing_segs(tactorData,uniqueCond,stim,buttonData,stimFromFile,fsStim,fsTact,trainTimesTotal,plotIt);
    
    %%
    % 9-13-2016 - script to compare response times once they've been calculated
    
    respLo = 0.100;
    respHi = 1;
    
    for i = 1:length(uniqueCond)
        
        trim = buttonLocs{i};
        trim = trim(trim>respLo & trim<respHi);
        zTrim = zscore(trim);
        buttonLocsThresh{i} = 1e3.*trim(abs(zTrim)<3);
        
    end
    
    %% Histogram
    
    % set number of bins
    nbins = 15;
    
    % make cell array for legends
    uniqueCondText = cellstr(num2str(uniqueCond));
    uniqueCondText{1} = 'no Prime';
    uniqueCondText{2} = 'Prime';
    
    %individual  histogram of each condition type
    figure
    
    for i = 1:length(uniqueCond)
        subplot(length(uniqueCond),1,i)
        a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
        title(['Histogram of reaction times for condition ' uniqueCondText{i} ])
        xlabel('Time (ms)')
        ylabel('Count')
        xlim([0 1000])
        a.BinWidth = 10;
        a.Normalization = 'probability';
    end
    subtitle([' block ', block])
    
    % overall histogram
    
    colormap lines;
    cmap = colormap ;
    
    figure
    hold on
    leg = {};
    for i = 1:length(uniqueCond)
        a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
        leg{end+1} = uniqueCondText{i};
        a.FaceColor = cmap(i,:);
        a.BinWidth = 10;
        a.Normalization = 'probability';
        xlim([0 800])
        
    end
    xlabel('Time (ms)')
    ylabel('Count')
    title(['Histogram of response Times for block ', block])
    
    legend(leg)
    
    % histogram just for tactor and 100,200,400,800
    
    
    %BOX PLOT
    
    % change colormap to matlab default lines
    
    combinedInfo = [];
    groups = [];
    colors = [];
    leg = {};
    temp = cmap(1,:);
    cmap(1,:) = cmap(2,:);
    cmap(2,:) = temp;
    
    keeps = [1 2];
    
    j = length(keeps);
    k = 1;
    
    for i = keeps
        combinedInfo = [buttonLocsThresh{i} combinedInfo ];
        groups = [(j).*ones(length(buttonLocsThresh{i}),1)' groups];
        colors(i,:) = cmap(i,:);
        leg{end+1} = uniqueCondText{i};
        j = j - 1;
        k = k + 1;
    end
    
    figure
    prettybox(combinedInfo,groups,colors,1,true)
    fig1 = gca;
    
    ylim([0 1000])
    fig1.XTick = [];
    legend(findobj(gca,'Tag','Box'),leg)
    ylabel('Response times (ms)')
    title(['Reaction Times for block ',block])
    
    %% statistics! kruskal wallis to start
    
    groupsKW = [];
    for i = keeps
        groupsKW = [cellstr(repmat(uniqueCondText{i},[length(buttonLocsThresh{i}),1])); groupsKW(:)];
    end
    
    [p,table,stats] = kruskalwallis(combinedInfo,groupsKW);
    [c,m,h,nms] = multcompare(stats);
    
    [nms num2cell(m)]
    
    c((c(:,6)<0.05),[1 2 6])
    
    
    %% save it
    
    current_direc = pwd;
    saveIt = 1;
    if saveIt
        save(fullfile(current_direc, [sid '_priming_behavior_block_' num2str(s) '.mat']),'buttonLocsSamps',...
            's','block','sid','primedOption','buttonLocs','tEpoch','stimTimes','fsStim','epochedButton',...
            'uniqueCond', 'respLo','respHi');
    end
    clearvars buttonLocsSamps primedOption buttonLocs tEpoch stimTimes fsStim epochedButton uniqueCond respLo respHi

end
