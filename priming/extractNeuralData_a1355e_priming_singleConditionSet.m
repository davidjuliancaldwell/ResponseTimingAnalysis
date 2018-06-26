%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;
% add path for scripts to work with data tanks
addpath('./scripts')

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\a1355e\data\d7\Converted_Matlab\ResponseTiming';
sid = SIDS{5};

plotIt = 1;

% ui box for input
list_str = {'1st block','2nd block'};

[s,v] = listdlg('PromptString','Pick experiment',...
    'SelectionMode','single',...
    'ListString',list_str);

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


%% neural data

% the eco data is crashing it right now
clearvars -except ECO1 ECO2 ECO3 Tact sid block s
eco1 = ECO1.data;
fsData = ECO1.info.SamplingRateHz;
ecoFs = fsData;
clear ECO1
eco2 = ECO2.data;
clear ECO2

eco3 = ECO3.data;
clear ECO3


data = [eco1 eco2 eco3];
clearvars eco1 eco2 eco3

% only 64 channels grid

data = data(:,1:64);

%
load([sid,'_Priming_block_',block,'.mat'])

% get train times

% look at stim from file saved (this is the sample where things were
% delivered
tact = Tact.data;
fs_tact = Tact.info.SamplingRateHz;
stimFromFile = tact(:,3);

% get stimulation times of delivery
trainTimes = find(stimFromFile~=0);

% convert sample times for eco

convertSamps = fs_tact/fsData;

trainTimesConvert = round(stimTimes/convertSamps);


%% cortical brain data
% ui box for input
prompt = {'Channel of interest?','condition'};
dlg_title = 'Channel of Interest';
num_lines = 1;
defaultans = {'8','0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

% options for this subject are 0 1
chanInt = str2num(answer{1});
condIntAns = str2num(answer{2});
condInt = find(uniqueCond==condIntAns);
% where to begin plotting with artifact
%artifact_end = round(0.05*eco_fs);
artifact_end = 0;

%where to end plotting
sampsEnd = round(2*ecoFs);

%presamps - where to begin looking for "rest" period (500 ms before?)
presamps = round(0.5*ecoFs);

trainTimesCell = {};
trainTimesCellThresh = {};

for i = 1:length(uniqueCond)
    
    trainTimesCell{i} = trainTimesConvert(primedOption==uniqueCond(i));
    trim = buttonLocs{i};
    trim = trim(trim>respLo & trim<respHi);
    zTrim = zscore(trim);
    if ~isempty(trainTimesCell{i}) % check to make sure not indexing empty cell
        %trainTimesCellThresh{i} = trainTimesCell{i}(abs(zTrim)<3); % z score
        % buttonLocsThresh = buttLocs{i}(abs(zTrim)<3);
        
        trainTimesCellThresh{i} = trainTimesCell{i};% no zscore
        buttonLocsThresh{i} = buttonLocs{i};% no zscore
        
    end
end

% ARTIFACT

post_stim = 2000;
samps_post_stim = round(post_stim/1e3*ecoFs);

pre_stim = 1000;
samps_pre_stim = round(pre_stim/1e3*ecoFs);

epochedCortEco = squeeze(getEpochSignal(data,(trainTimesCellThresh{condInt})-samps_pre_stim,(trainTimesCellThresh{condInt}+ samps_post_stim)));

response = buttonLocsThresh{condInt};
%t_epoch = [1:size(epochedCortEco,1)]/eco_fs;
t_epoch = (-samps_pre_stim:samps_post_stim-1)/ecoFs;

exampChan = mean(squeeze(epochedCortEco(:,chanInt,:)),2);

figure
plot(1e3*t_epoch,exampChan);
xlim([-1000 2000])
ylim([-10e-5 10e-5])
title(['Subject ' sid ' Channel ' num2str(chanInt) ' Condition ' num2str(condIntAns)])
%clear exampChan

dataInt = epochedCortEco;
fsData = ecoFs; 
stimChans = [16 24];
save('a1355e_examplePriming_Prime_low.mat','dataInt','fs_data','t_epoch','stimChans')
