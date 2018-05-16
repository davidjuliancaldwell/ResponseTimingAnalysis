%% starting with subject a

% % DJC Desktop
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';

%% load in subject

% this is from my z_constants
Z_ConstantsStimResponse;

sid = SIDS{1};
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';

% ui box for input
list_str = {'sensory stimulation','tactor stimulation','off target stimulation'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

% load in data
if (strcmp(sid, 'acabb1'))
    folder_data = strcat(DATA_DIR,'\acabb1');
    
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

[stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact);

clear Stim Tact Sing


if s~=2 
%% plot stim

figure
hold on
for i = 1:size(stim,2)
    
    t = (0:length(stim)-1)/fs_stim;
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
    
    samplesOfPulse = round(2*fs_stim/1e3);
    
    
    
    % trying something like A_BuildStimTables from BetaStim
    
    
    Sing1Mask = Sing1~=0;
    dmode = diff([0 Sing1Mask' 0 ]);
    
    
    dmode(end-1) = dmode(end);
    
    
    bursts(2,:) = find(dmode==1);
    bursts(3,:) = find(dmode==-1);
    
    stims = squeeze(getEpochSignal(Sing1,(bursts(2,:)-1),(bursts(3,:))+1));
    t = (0:size(stims,1)-1)/fs_sing;
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
    t = (0:size(stim1Epoched,1)-1)/fs_stim;
    t = t*1e3;
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery')
    
    % hold on
    
    plot(t,stims)
    
    % get the delay in stim times
    
    delay = round(0.2867*fs_stim/1e3);
    
    % plot the appropriately delayed signal
    stimTimesBegin = bursts(2,:)-1+delay;
    stimTimesEnd = bursts(3,:)-1+delay;
    stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd+5));
    t = (0:size(stim1Epoched,1)-1)/fs_stim;
    t = t*1e3;
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in')
    
elseif s==2
    stim1 = stim(:,1);
end


%% extract data
if s~=2
    % try and account for delay for the stim times
    stimTimes = bursts(2,:)-1+delay;
    
    % DJC 7-7-2016, changed presamps and post samps to 1 second
    presamps = round(1 * fs_data); % pre time in sec
    postsamps = round(1 * fs_data); % post time in sec, % modified DJC to look at up to 300 ms after
    
    
    % sampling rate conversion between stim and data
    fac = fs_stim/fs_data;
    
    % find times where stims start in terms of data sampling rate
    sts = round(stimTimes / fac);
end

%% look at all simultaneously

tactorData = tact(:,1);
buttonData = tact(:,2);

analyze_all_inputs_simultaneously(tactorData,buttonData,stim,stimFromFile,fsTact)

%% for 1st subject - only look at parts where t > 50 for sensory stim, t > 12 for tactor  (t_begin)

if s==1
    t_begin = 40;
elseif s ==2
    t_begin = 12;
else 
    t_begin = 0;
end;
buttonData = buttonData(t_button>t_begin);
tactorData = tactorData(t_tact>t_begin);
stimFromFile = stimFromFile(t_stimFile>t_begin);
stim1 = stim1(t_stimFile>t_begin);

t_button = t_button(t_button>t_begin);
t_tact = t_tact(t_tact>t_begin);
t_stimFile = t_stimFile(t_stimFile>t_begin);

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
[buttonPks,buttonLocs] = findpeaks(buttonDataClip,t_button,'MinpeakDistance',2,'Minpeakheight',0.008);

figure
findpeaks(buttonDataClip,t_button(t_button>t_begin),'MinpeakDistance',2,'Minpeakheight',0.008);
hold on
plot(t_stimFile,stimFromFile,'g');
plot(t_stimFile,stim1,'r');

% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)
%% QUANTIFY RXN TIME TO CORTICAL STIM
% get epochs for button press, with start being onset of stimulation marker
if s == 1 || s == 3
    %
    trainTimes = find(stimFromFile~=0);
    
    % shrink condition type to be 120
    condType = condType(1:120);
    
    % pick condition type where stimulation was delivered
    if s == 1 || s== 3
        trainTimesCond1 = trainTimes(condType==0);
    elseif s == 2
        trainTimesCond1 = trainTimes(condType==0 | condType==1);
    end
    
    sampsEnd = round(2*fs_stim);
    
    % epoched button press
    epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEnd)));
    
    figure
    t_epoch = [0:size(epochedButton,1)-1]/fs_stim;
    plot(t_epoch,epochedButton);
    
    % vector of pks of button press
    
    buttonPksVecCort = zeros(size(epochedButton,2),1);
    buttonLocsVecCort = zeros(size(epochedButton,2),1);
    
    clear buttonPksTemp buttonLocsTemp tactorPksTemp tactorLocsTemp;

    for i = 1:size(epochedButton,2)
        [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
        if isempty(buttonPksTemp)
            buttonPksTemp = NaN;
            buttonLocsTemp = NaN;
        end
        buttonPksVecCort(i) = buttonPksTemp;
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
    epochedButtonNull = squeeze(getEpochSignal(buttonDataClip,trainTimesNull,(trainTimesNull + sampsEnd)));
    
        figure
    t_epoch = [0:size(epochedButtonNull,1)-1]/fs_stim;
    plot(t_epoch,epochedButtonNull);
    
    % vector of pks of button press
    
    buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
    buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
    
    clear buttonPksTemp buttonLocsTemp tactorPksTemp tactorLocsTemp;

    for i = 1:size(epochedButtonNull,2)
        [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButtonNull(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
        if isempty(buttonPksTemp)
            buttonPksTemp = NaN;
            buttonLocsTemp = NaN;
        end
        buttonPksVecCortNull(i) = buttonPksTemp;
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
    
    sampsEndButton = round(3*fs_stim);
    sampsEndTactor = round(2*fs_stim);
    
    % epoched button press
    epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEndButton)));
    
    % if tactor stim, epoch that too - FINISH THIS
    epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesCond1,(trainTimesCond1 + sampsEndTactor)));
    
    figure
    t_epoch_button = [0:size(epochedButton,1)-1]/fs_stim;
    t_epoch_tact = [0:size(epochedTactor,1)-1]/fs_stim;
    
    plot(t_epoch_button,epochedButton);
    
    % vector of pks of button press
    
    buttonPksVecTact = zeros(size(epochedButton,2),1);
    buttonLocsVecTact = zeros(size(epochedButton,2),1);
    
    % vector of pks of tactor press
    
    tactorPksVecTact = zeros(size(epochedTactor,2),1);
    tactorLocsVecTact = zeros(size(epochedTactor,2),1);
    
    clear buttonPksTemp buttonLocsTemp tactorPksTemp tactorLocsTemp;

    
    for i = 1:size(epochedButton,2)
        [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton(:,i),t_epoch_button,'NPeaks',1,'Minpeakheight',0.008);
        [tactorPksTemp,tactorLocsTemp] = findpeaks(epochedTactor(:,i),t_epoch_tact,'NPeaks',1,'Minpeakheight',2);
        
        if isempty(buttonPksTemp)
            buttonPksTemp = NaN;
            buttonLocsTemp = NaN;
        end
        
        if isempty(tactorPksTemp)
            tactorPksTemp = NaN;
            tactorLocsTemp = NaN;
        end
        
        buttonPksVecTact(i) = buttonPksTemp;
        buttonLocsVecTact(i) = buttonLocsTemp;
        tactorPksVecTact(i) = tactorPksTemp;
        tactorLocsVecTact(i) = tactorLocsTemp;
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
    t_epoch = [0:size(epochedButtonNull,1)-1]/fs_stim;
    plot(t_epoch,epochedButtonNull);
    
    % vector of pks of button press
    
    buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
    buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
       
    buttonPksVecCortNull = zeros(size(epochedButtonNull,2),1);
    buttonLocsVecCortNull = zeros(size(epochedButtonNull,2),1);
    
    clear buttonPksTemp buttonLocsTemp tactorPksTemp tactorLocsTemp;

    for i = 1:size(epochedButtonNull,2)
        [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButtonNull(:,i),t_epoch,'NPeaks',1,'Minpeakheight',0.008);
        if isempty(buttonPksTemp)
            buttonPksTemp = NaN;
            buttonLocsTemp = NaN;
        end
        buttonPksVecCortNull(i) = buttonPksTemp;
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

% clear all variables except the ones that are useful for further
% iterations 

%clearvars -except buttonTactDiff buttonLocsVectTact tactorLocsVecTact buttonLocsVecCort respLo respHi SIDS DATA_DIR sid


