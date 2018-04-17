%% script to extract response timing data - 8-2-2016
% clear workspace
close all; clear all; clc

%  initialize output and meta dir
% set input output working directories - for David's PC right now
Z_ConstantsStimResponse;

% DJC Desktop
DATA_DIR = 'C:\Users\djcald.CSENETID\Data\ConvertedTDTfiles';

sid = SIDS{1};
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

%%
if strcmp(sid,'693ffd') % third subject
    %% add the right path
    
    addpath(genpath('693ffd'))
    %% extract the response times
    
    extractStimResponse_693ffd
    %% compare the response times
    
    compareRespTimes_693ffd
    
    %% extract the neural data
    
    %extractNeuralData_c19968
    
    
    
end

if strcmp(sid,'2fd831') % third subject
    %% add the right path
    
    addpath(genpath('2fd831'))
    %% extract the response times
    
    extractStimResponse_2fd831
    %% compare the response times
    
    compareRespTimes_2fd831
    
    %% extract the neural data
    
    %extractNeuralData_2fd831
end

if strcmp(sid,'a1355e')
    
    addpath(genpath('a1355e'))
    %% extract the response times
    
    extractStimResponse_a1355e
    %% compare the response times
    
    compareRespTimes_a1355e
end






