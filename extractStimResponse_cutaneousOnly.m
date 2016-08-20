%% only cutaneous reacton times - healthy humans - 8-18-2016

%% initialize output and meta dir - THIS IS ALL IN THE MAJRO
% % clear workspace
% close all; clear all; clc
%
% % set input output working directories - for David's PC right now
% % Z_ConstantsStimResponse;
%
% % add path for scripts to work with data tanks
% addpath('./scripts')
%
% % subject directory, change as needed
% % SUB_DIR = fullfile(myGetenv('subject_dir')); - for David's PC right now
%
% % data directory
%
% %PUT PATH TO DATA DIRECTORY WITH CONVERTED DATA FILES
%
% % DJC Desktop
% DATA_DIR = 'C:\Users\djcald\Data\ConvertedTDTfiles';
%
% % DJC Laptop
% %DATA_DIR = 'C:\Users\David\GoogleDriveUW\GRIDLabDavidShared\ResponseTiming';
%
% SIDS = {'acabb1'};

%% initialiaze names cell - only do this once!
% 

% close all;clear all;clc
% names = {};
% differenceTotal = {};

%% load in subject


% load in data interactively 
DATA_DIR = 'C:\Users\djcald\Google Drive\GRIDLabDavidShared\StimulationSpacing\Tactile Reaction Time';
cd (DATA_DIR)
dataTotal = uiimport('-file');

% get subject name
prompt = {'Enter subject name - and which hand for button'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'brandon_L'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
name = answer{1};

% make vector of total names 


names = {names{:},name};


%% load in data of interest
tact = dataTotal.Tact.data;
fs_tact = dataTotal.Tact.info.SamplingRateHz;

%% look at tactor

tactorData = tact(:,1);
t_tact = (0:length(tactorData)-1)/fs_tact;
% figure
% plot(t_tact,tactorData);
% 
% title('tactor data')

% look at button press

buttonData = tact(:,2);
t_button = (0:length(buttonData)-1)/fs_tact;
% figure
% plot(t_button,buttonData);
% 
% title('button data')

% look at stim from file saved

stimFromFile = tact(:,3);
t_stimFile = (0:length(buttonData)-1)/fs_tact;
% figure
% plot(t_stimFile,stimFromFile);
% title('stim from file')

% look at tone out

tone = tact(:,4);
t_tone = (0:length(buttonData)-1)/fs_tact;

% look at all 3 + tone 

% figure
% ax1 = subplot(4,1,1);
% plot(t_tact,tactorData)
% title('tactor data')
% 
% ax2 = subplot(4,1,2);
% plot(t_button,buttonData);
% title('button data')
% 
% ax3 = subplot(4,1,3);
% plot(t_stimFile,stimFromFile);
% title('stim from file')
% 
% ax4 = subplot(4,1,4);
% plot(t_tact,tone);
% 
% 
% 
% %link axis
% linkaxes([ax1,ax2,ax3,ax4],'x')

%% for 1st subject - 
% change t_begin if need be to only look at later trials 

t_begin = 0; 

buttonData = buttonData(t_button>t_begin);
tactorData = tactorData(t_tact>t_begin);
stimFromFile = stimFromFile(t_stimFile>t_begin);

t_button = t_button(t_button>t_begin);
t_tact = t_tact(t_tact>t_begin);
t_stimFile = t_stimFile(t_stimFile>t_begin);

% set respLo and respHi, which are the values which for less or greater
% than rxn times aren't analyzed

respLo = 0.100;
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

% figure
% findpeaks(buttonDataClip,t_button(t_button>t_begin),'MinpeakDistance',2,'Minpeakheight',0.008);
% hold on
% plot(t_stimFile,stimFromFile,'g');



% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)

%% RXN TIME FOR TACTOR STIMULATION

% last one seems wonky? 

trainTimes = find(stimFromFile~=0);
trainTimes = trainTimes(1:end-1);

% shrink condition type to be 120

% why not able to go to the end? 
condType = condType(1:119);

% pick condition type where stimulation was delivered

trainTimesCond1 = trainTimes(condType==0);


% different sample end for tactor and button to account for double
% delay

sampsEndButton = round(3*fs_tact);
sampsEndTactor = round(2*fs_tact);

% epoched button press
epochedButton = squeeze(getEpochSignal(buttonDataClip,trainTimesCond1,(trainTimesCond1 + sampsEndButton)));

% if tactor stim, epoch that too - FINISH THIS
epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesCond1,(trainTimesCond1 + sampsEndTactor)));

% figure
t_epoch_button = [0:size(epochedButton,1)-1]/fs_tact;
t_epoch_tact = [0:size(epochedTactor,1)-1]/fs_tact;

% plot(t_epoch_button,epochedButton);

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
% figure
% %set number of bins for histogram
% nbins = 15;
% histogram(tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi),nbins);
% title('Histogram of reaction times for tactor')
% xlabel('Time (seconds')
% ylabel('Count')
% 
% % histogram of rxn times for button press - tactor at each
% % corresponding epoch
% figure
% histogram(buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi),nbins)
% title('Histogram of reaction times for button relative to tactor onset')
% xlabel('Time (seconds')
% ylabel('Count')

% save train delivery times for brain data


% clear all variables except the ones that are useful for further
% iterations


tactorLocsVecTactTrim = tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi);
buttonTactDiffTrim = buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi);

zTact = zscore(tactorLocsVecTactTrim);
zDiff = zscore(buttonTactDiffTrim);

tactor = 1e3.*tactorLocsVecTactTrim(abs(zTact)<3);
difference = 1e3.*buttonTactDiffTrim(abs(zDiff)<3);

differenceTotal{end+1} = difference;

clearvars -except names differenceTotal iteration



