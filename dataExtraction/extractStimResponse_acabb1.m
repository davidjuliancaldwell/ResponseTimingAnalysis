%% starting with subject a
% load in subject

% this is from my z_constants
Z_ConstantsStimResponse;

sid = SIDS{1};
subjdir = getenv('SUBJECT_DIR');
DATA_DIR = fullfile(subjdir,'\ConvertedTDTfiles');

% ui box for input
%list_str = {'sensory stimulation','tactor stimulation','off target stimulation'};

%[s,v] = listdlg('PromptString','Pick experiment',...
%    'SelectionMode','single',...
%    'ListString',list_str);

for s = 1:2
    % load in data
    if (strcmp(sid, 'acabb1'))
        folder_data = strcat(DATA_DIR,'\',sid);
        
        if s == 1
            load(fullfile(folder_data,'ReactionTime-1.mat'))
        elseif s == 2
            load(fullfile(folder_data,'ReactionTime-3.mat'))
        elseif s == 3
            load(fullfile(folder_data,'ReactionTime-4.mat'))
            
        end
        
    end
    
    plotIt = 1;
    
    %% load in data of interest
    
    [stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,Wave,Tact);
    
    clear Stim Tact Sing
    
    if s~=2
        %% plot stim
        
        figure
        hold on
        for i = 1:size(stim,2)
            
            t = (0:length(stim)-1)/fsStim;
            subplot(2,2,i)
            plot(t*1e3,stim(:,i))
            title(sprintf('Channel %d',i))
            
            
        end
        
        
        xlabel('Time (ms)')
        ylabel('Amplitude (V)')
        
        subtitle('Stimulation Channels')
        %%Sing looks like the wave to be delivered, with amplitude in uA
        % Try working from this - do this if not tactor stim
        
        % build a burst table with the timing of stimuli
        bursts = [];
        
        Sing1 = sing(:,1);
        
        samplesOfPulse = round(2*fsStim/1e3);
        
        
        
        % trying something like A_BuildStimTables from BetaStim
        
        
        Sing1Mask = Sing1~=0;
        dmode = diff([0 Sing1Mask' 0 ]);
        
        
        dmode(end-1) = dmode(end);
        
        
        bursts(2,:) = find(dmode==1);
        bursts(3,:) = find(dmode==-1);
        
        stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
        t = (0:size(stims,1)-1)/fsStim;
        t = t*1e3;
        figure
        plot(t,stims,'b','linewidth',2)
        xlabel('Time (ms)');
        ylabel('Current to be delivered (mA)')
        ylim([(min(stims(:))-100) (max(stims(:))+100)])
        title('Current to be delivered for all trials')
        
        % delay loks to be 0.2867 ms from below.
        
        %% Plot stims with info from above
        
        stim1 = stim(:,1);
        stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,:)-1),(bursts(3,:))+1));
        t = (0:size(stim1Epoched,1)-1)/fsStim;
        t = t*1e3;
        figure
        plot(t,stim1Epoched)
        xlabel('Time (ms)');
        ylabel('Voltage (V)');
        title('Finding the delay between current output and stim delivery')
        
        % hold on
        
        plot(t,stims)
        
        % get the delay in stim times
        
        delay = round(0.2867*fsStim/1e3);
        
        % plot the appropriately delayed signal
        stimTimesBegin = bursts(2,:)-1+delay;
        stimTimesEnd = bursts(3,:)-1+delay;
        stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd+5));
        t = (0:size(stim1Epoched,1)-1)/fsStim;
        t = t*1e3;
        figure
        plot(t,stim1Epoched)
        xlabel('Time (ms)');
        ylabel('Voltage (V)');
        title('Stim voltage monitoring with delay added in')
        
    elseif s==2
        stim1 = stim(:,1);
    end
    
    
    % %% extract data
    % if s~=2
    %     % try and account for delay for the stim times
    %     stimTimes = bursts(2,:)-1+delay;
    %
    %     % DJC 7-7-2016, changed presamps and post samps to 1 second
    %     presamps = round(1 * fs_data); % pre time in sec
    %     postsamps = round(1 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after
    %
    %
    %     % sampling rate conversion between stim and data
    %     fac = fsStim/fs_data;
    %
    %     % find times where stims start in terms of data sampling rate
    %     sts = round(stimTimes / fac);
    % end
    
    %% look at all simultaneously
    
    tactorData = tact(:,1);
    buttonData = tact(:,2);
    
    %analyze_all_inputs_simultaneously(tactorData,buttonData,stim,stimFromFile,fsTact)
    tactorData = tact(:,1);
    t_tact = (0:length(tactorData)-1)/fsTact;
    figure
    plot(t_tact,tactorData);
    
    title('tactor data')
    
    % look at button press
    
    buttonData = tact(:,2);
    t_button = (0:length(buttonData)-1)/fsTact;
    figure
    plot(t_button,buttonData);
    
    title('button data')
    
    % look at stim from file saved
    
    stimFromFile = tact(:,3);
    t_stimFile = (0:length(stim)-1)/fsTact;
    figure
    plot(t_stimFile,stimFromFile);
    title('stim from file')
    
    % look at all 3 + stim waveform
    
    figure
    ax1 = subplot(4,1,1);
    plot(t_tact,tactorData)
    title('tactor data')
    
    ax2 = subplot(4,1,2);
    plot(t_button,buttonData);
    title('button data')
    
    ax3 = subplot(4,1,3);
    plot(t_stimFile,stimFromFile);
    title('stim from file')
    
    
    % assuming stim1 here is the channel where stim was being delivered
    ax4 = subplot(4,1,4);
    t_stimFile = t_stimFile(1:length(stim1)); % DJC 4-15-2018 for the off target condition
    plot(t_stimFile,stim1);
    
    %link axis
    linkaxes([ax1,ax2,ax3,ax4],'x')
    
    %% for 1st subject - only look at parts where t > 50 for sensory stim, t > 12 for tactor  (t_begin)
    
    if s==1
        tBegin = 40;
    elseif s ==2
        tBegin = 12;
    else
        tBegin = 0;
    end;
    
    tButton = (0:length(buttonData)-1)*fsStim;
    
    buttonData = buttonData(tButton>tBegin);
    tactorData = tactorData(tButton>tBegin);
    stimFromFile = stimFromFile(tButton>tBegin);
    stim1 = stim1(tButton>tBegin);
    
    tButton = tButton(tButton>tBegin);
    
    % set respLo and respHi, which are the values which for less or greater
    % than rxn times aren't analyzed
    
    respLo = 0.150;
    respHi = 1;
    
    %% 8-12-2016 - start quantifying data
    
    % vector of condition type - for first subject, looks like condition type
    % is what was used , rather than test_condition
    
    condType= dlmread('condition.txt');
    condTestType = dlmread('test_condition.txt');
    train = dlmread('testTrain.txt');
    
    % find button data peaks
    
    % set above certain threshold to 0.009
    buttonDataClip = buttonData;
    buttonDataClip(buttonData >= 0.009) = 0.009;
    [buttonPks,buttonLocs] = findpeaks(buttonDataClip,tButton,'MinpeakDistance',2,'Minpeakheight',0.008);
    
    figure
    findpeaks(buttonDataClip,tButton(tButton>tBegin),'MinpeakDistance',2,'Minpeakheight',0.008);
    hold on
    plot(tButton,stimFromFile,'g');
    plot(tButton,stim1,'r');
    
    % raw button
    %[buttonPks,buttonLocs] = findpeaks(buttonData,tButton,'MinpeakDistance',1.5)
    %findpeaks(buttonData,tButton,'MinpeakDistance',2,'Minpeakheight',10e-3)
    %% QUANTIFY RXN TIME TO CORTICAL STIM
    % get epochs for button press, with start being onset of stimulation marker
    if s == 1 || s == 3
        %
        trainTimes = find(stimFromFile~=0);
        
        % shrink condition type to be 120
        condType = condType(1:120);
        
        % pick condition type where stimulation was delivered
        if s == 1
            trainTimesCond1 = trainTimes(condType==0);
        elseif s == 2
            trainTimesCond1 = trainTimes(condType==0 | condType==1);
        elseif s == 3
            condType = condType(1:26);
            
            trainTimesCond1 = trainTimes(condType==0);
            trainTimes = trainTimes(1:26);
            
        end
        
        sampsEnd = round(2*fsStim);
        
        % epoched button press
        epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEnd)));
        
        figure
        tEpoch = [0:size(epochedButton,1)-1]/fsStim;
        plot(tEpoch,epochedButton);
        
        % vector of pks of button press
        
        buttonPksVecCort = zeros(size(epochedButton,2),1);
        buttonLocsVecCort = zeros(size(epochedButton,2),1);
        
        clear buttonLocsTemp tactorLocsTemp;
        
             buttonStart = 0.1;
        buttonEnd = 2.5;
        
        for i = 1:size(epochedButton,2)
            
            
            [ipt,residual] = findchangepts((epochedButton(tEpoch>buttonStart & tEpoch<buttonEnd,i)),'maxnumchanges',2);
            if isempty(ipt) || max(epochedButton(tEpoch>buttonStart & tEpoch<buttonEnd,i)) < 8e-3
                ipt = NaN;
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            
            buttonLocsTempSamps = ipt(1)+round(buttonStart*fsTact);
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            % [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
            if isempty(buttonLocsTemp)
                %    buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            % buttonPksVecCort(i) = buttonPksTemp;
            buttonLocsVecCort(i) = buttonLocsTemp;
        end
        
        % histogram of rxn times, assume 200 ms or greater  & less than 1 s
        figure
        % set number of bins
        nbins = 15;
        histogram(buttonLocsVecCort(buttonLocsVecCort>respLo & buttonLocsVecCort<respHi ),nbins)
        title('Histogram of reaction times')
        xlabel('Time (seconds')
        ylabel('Count')
        %%%%%%%%%%
        % DJC 4-15-2018
        trainTimesNull = trainTimes(condType == 2);
        if s == 3
            sampsEnd = 24414;
        end
        epochedButtonNull = squeeze(getEpochSignal(buttonDataClip,trainTimesNull,(trainTimesNull + sampsEnd)));
        
        figure
        tEpoch = [0:size(epochedButtonNull,1)-1]/fsStim;
        plot(tEpoch,epochedButtonNull);
        
        % vector of pks of button press
        
        buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
        buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
        
        clear buttonLocsTemp tactorLocsTemp;
        
   
        for i = 1:size(epochedButtonNull,2)
            
            % [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButtonNull(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
            
            [ipt,residual] = findchangepts((epochedButtonNull(tEpoch>buttonStart&tEpoch<buttonEnd,i)),'maxnumchanges',2);
            
            if isempty(ipt) || max((epochedButtonNull(tEpoch>buttonStart&tEpoch<buttonEnd,i))) < 8e-3
                ipt = NaN;
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            %  buttonLocsTempSamps = ipt(1)+round(1*fsTact);
            buttonLocsTempSamps = ipt(1)+round(buttonStart*fsTact);
            
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            
            if isempty(buttonLocsTemp)
                %     buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            
            %  buttonPksVecCortNull(i) = buttonPksTemp;
            buttonLocsVecCortNull(i) = buttonLocsTemp;
        end
        
        % histogram of rxn times, assume 200 ms or greater  & less than 1 s
        figure
        % set number of bins
        nbins = 15;
        histogram(buttonLocsVecCortNull(buttonLocsVecCortNull>respLo & buttonLocsVecCortNull<respHi ),nbins)
        title('Histogram of reaction times')
        xlabel('Time (seconds')
        ylabel('Count')
        
    end
    %% RXN TIME FOR TACTOR STIM
    % get epochs for button press, with start being onset of stimulation marker
    
    % include tactor delay
    useTactorDelay = 0;
    
    if useTactorDelay
        tactorDelaySamps = (1.04/1e3)*fsStim; % ms
        tactorDelaySecs = 1.04/1e3;
    else
        tactorDelaySamps = 0;
        tactorDelaySecs = 0;
    end
    
    %
    if s == 2
        trainTimes = find(stimFromFile~=0);
        
        % shrink condition type to be 120
        condType = condType(1:120);
        
        % pick condition type where stimulation was delivered
        if s == 1
            trainTimesCond1 = trainTimes(condType==0);
        elseif s == 2
            trainTimesCond1 = trainTimes(condType==0 | condType==1);
        end
        
        % different sample end for tactor and button to account for double
        % delay
        
        sampsEndButton = round(3.5*fsStim);
        sampsEndTactor = round(3.5*fsStim);
        
        % epoched button press
        epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEndButton)));
        
        % if tactor stim, epoch that too - FINISH THIS
        epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesCond1,(trainTimesCond1 + sampsEndTactor)));
        
        figure
        tEpochButton = [0:size(epochedButton,1)-1]/fsStim;
        tEpochTact = [0:size(epochedTactor,1)-1]/fsStim;
        
        plot(tEpochButton,epochedButton);
        
        % vector of pks of button press
        
        buttonPksVecTact = zeros(size(epochedButton,2),1);
        buttonLocsVecTact = zeros(size(epochedButton,2),1);
        
        % vector of pks of tactor press
        
        % tactorPksVecTact = zeros(size(epochedTactor,2),1);
        tactorLocsVecTact = zeros(size(epochedTactor,2),1);
        
        clear buttonLocsTemp tactorLocsTemp;
        
        buttonStart = 0;
        for i = 1:size(epochedButton,2)
            
            [ipt,residual] = findchangepts((epochedButton(tEpochButton>buttonStart,i)),'maxnumchanges',2);
            if isempty(ipt) || max(epochedButton(tEpochButton>buttonStart,i)) < 8e-3
                ipt = NaN;
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            
%                              figure
%                           findchangepts((epochedButton(tEpochButton>buttonStart,i)),'maxnumchanges',2)
%             
            buttonLocsTempSamps = ipt(1)+round(buttonStart*fsTact);
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            % [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton(:,i),t_epoch_button,'NPeaks',1,'Minpeakheight',0.008);
            [tactorPksTemp,tactorLocsTemp] = findpeaks(epochedTactor(:,i),tEpochTact,'NPeaks',1,'Minpeakheight',2);
            
            if isempty(buttonLocsTemp)
                %buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            
            if isempty(tactorLocsTemp)
                %    tactorPksTemp = NaN;
                tactorLocsTemp = NaN;
            end
            
            %  buttonPksVecTact(i) = buttonPksTemp;
            buttonLocsVecTact(i) = buttonLocsTemp;
            %  tactorPksVecTact(i) = tactorPksTemp;
            tactorLocsVecTact(i) = tactorLocsTemp - tactorDelaySecs;
        end
        %%
        % calculate differences
        
        buttonTactDiff = buttonLocsVecTact - tactorLocsVecTact;
        
        % histogram of rxn times for tacxtor , assume 200 ms or greater, and
        % less than 1 s
        figure
        %set number of bins for histogram
        nbins = 15;
        histogram(tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi),nbins);
        title('Histogram of reaction times for tactor')
        xlabel('Time (seconds')
        ylabel('Count')
        
        % histogram of rxn times for button press - tactor at each
        % corresponding epoch
        figure
        histogram(buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi),nbins)
        title('Histogram of reaction times for button relative to tactor onset')
        xlabel('Time (seconds')
        ylabel('Count')
        
        % save train delivery times for brain data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% null trial check
        
        % DJC 4-15-2018
        trainTimesNull = trainTimes(condType == 2);
        epochedButtonNull = squeeze(getEpochSignal(buttonDataClip,trainTimesNull,(trainTimesNull + sampsEndButton)));
        
        figure
        tEpoch = [0:size(epochedButtonNull,1)-1]/fsStim;
        plot(tEpoch,epochedButtonNull);
        
        % vector of pks of button press
        
        buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
        buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
        
        buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
        buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
        
        clear buttonLocsTemp tactorLocsTemp;
        
        for i = 1:size(epochedButtonNull,2)
            % [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButtonNull(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
            
            [ipt,residual] = findchangepts((epochedButtonNull(:,i)),'maxnumchanges',2);
            
            if isempty(ipt) || max((epochedButtonNull(:,i))) < 8e-3
                ipt = NaN;
                
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            %  buttonLocsTempSamps = ipt(1)+round(1*fsTact);
            buttonLocsTempSamps = ipt(1);
            
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            if isempty(buttonLocsTemp)
                %buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            %buttonPksVecCortNull(i) = %buttonPksTemp;
            buttonLocsVecCortNull(i) = buttonLocsTemp;
        end
        
        % histogram of rxn times, assume 200 ms or greater  & less than 1 s
        figure
        % set number of bins
        nbins = 15;
        histogram(buttonLocsVecCortNull(buttonLocsVecCortNull>respLo & buttonLocsVecCortNull<respHi ),nbins)
        title('Histogram of reaction times')
        xlabel('Time (seconds')
        ylabel('Count')
        
        
        
    end
    
    % % % % % % % %     % clear all variables except the ones that are useful for further
    % iterations
    
    clearvars -except buttonTactDiff buttonLocsVectTact tactorLocsVecTact buttonLocsVecCort respLo respHi SIDS DATA_DIR sid
    
end

respLo = 0.1; % djc 4-15-2018
respHi = 1; % djc 4-15-2018
tactorLocsVecTactTrim = tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi);
buttonTactDiffTrim = buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi);
buttonLocsVecCortTrim = buttonLocsVecCort(buttonLocsVecCort>respLo & buttonLocsVecCort<respHi );

zTact = zscore(tactorLocsVecTactTrim);
zDiff = zscore(buttonTactDiffTrim);
zCort = zscore(buttonLocsVecCortTrim);

%tactor = 1e3.*tactorLocsVecTactTrim(abs(zTact)<3);
%difference = 1e3.*buttonTactDiffTrim(abs(zDiff)<3);
%cort = 1e3.*buttonLocsVecCortTrim(abs(zCort)<3);

% DJC 12-9-2016 - do thresholding later if wanted, to match the other two
% subjects
tactor = 1e3.*tactorLocsVecTact;
difference = 1e3.*buttonTactDiff;
cort = 1e3.*buttonLocsVecCort;

current_direc = pwd;

%save(fullfile(current_direc, [sid '_compareResponse_tactorSub.mat']), 'tactor', 'difference', 'cort','tactorLocsVecTactTrim','buttonTactDiffTrim','buttonLocsVecCortTrim','buttonLocsVecCort','tactorLocsVecTact','buttonTactDiff','respLo','respHi');

save(fullfile(current_direc, [sid '_compareResponse_changePts_noDelay.mat']), 'tactor', 'difference', 'cort','tactorLocsVecTactTrim','buttonTactDiffTrim','buttonLocsVecCortTrim','buttonLocsVecCort','tactorLocsVecTact','buttonTactDiff','respLo','respHi');