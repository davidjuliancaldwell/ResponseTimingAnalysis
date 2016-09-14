%% script to extract response timing data - 8-2-2016


%%
% clear workspace
close all; clear all; clc

%%  initialize output and meta dir
% set input output working directories - for David's PC right now
Z_ConstantsStimResponse;

% add path for scripts to work with data tanks
addpath('./scripts')

% subject directory, change as needed
% SUB_DIR = fullfile(myGetenv('subject_dir')); - for David's PC right now

% data directory

%PUT PATH TO DATA DIRECTORY WITH CONVERTED DATA FILES

% DJC Desktop
DATA_DIR = 'C:\Users\djcald\Data\ConvertedTDTfiles';

% DJC Laptop
%DATA_DIR = 'C:\Users\David\GoogleDriveUW\GRIDLabDavidShared\ResponseTiming';

sid = SIDS{2};
%%
if strcmp(sid,'acabb1') % first subject 
    %% add the right path
    
    addpath(genpath('acabb1'))
    %% extract the response times
    
    extractStimResponse_acabb1
    %% compare the response times
    
    compareRespTimes_acabb1
    
    %% extract the neural data
    
    extractNeuralData_acabb1
end
%%
if strcmp(sid,'c19968') % second subject
    %% add the right path
    
    addpath(genpath('c19968'))
    %% extract the response times
    
    extractStimResponse_c19968
    %% compare the response times
    
    compareRespTimes_c19968
    
    %% extract the neural data
    
    %extractNeuralData_c19968
    
    
    
end




