% first load in the MATLAB file

load('resultTrim.mat')

% sid of interest, experiment, block
% change these to match what you are interested in! 
sidInt = 'acabb1';
experimentInt = '200 ms ';
blockInt = 1;

% this makes a vector of indices which is then passed to responsetimems to
% get the response times of interest for a subject. 

rows = find(strcmp(sidInt,sid) & strcmp(experimentInt,experiment) & blockInt==block);

responseTimeSubj = responsetimems(rows)