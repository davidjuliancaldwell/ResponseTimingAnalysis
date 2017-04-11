function [powerout,f,t,phase_angle] = neural_analysis_calcs(signal,fs,time_res)
% This is a function to run on the response timing data collected by David
% Caldwell and Jeneva Cronin while in the GRIDLab. This uses a
% morletprocess script as implemented by James Wu
% signal:
%   time x channels x trials
% fs:
%   sampling rate in Hz
% time_res:
%   time resolution for morlet process
%
% OUTPUT 
% powerout - freq x time x channel x trial

%
% DJC 3-30-2017

% Morlet Process
%%%%%%%%%%%%%%%%%%%
num_trials = size(signal,3);
powerout = [];



for i = 1:num_trials
    
    data_temp = signal(:,:,i);
    % compute f and t once
    if i == 1
        [powerout_temp, f, t] = morletprocess(data_temp, fs, time_res);
        powerout(:,:,:,i) = powerout_temp;
    end
    [powerout_temp] = morletprocess( data_temp, fs, time_res);
    powerout(:,:,:,i) = powerout_temp;
    
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set phase angle output to be zero 

if ~exist('phase_angle','var')
   phase_angle = []; 
end

end