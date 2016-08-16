%% script to extract response timing data - 8-2-2016


%%
% clear workspace
close all; clear all; clc

%%  initialize output and meta dir
% set input output working directories - for David's PC right now
% Z_ConstantsStimResponse;

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

SIDS = {'acabb1'};

%% extract the response times
extractStimResponse

% once you have run this and gotten your vectors of interest can do compare
% 

%% compare the response times 

compareRespTimes

%%

n