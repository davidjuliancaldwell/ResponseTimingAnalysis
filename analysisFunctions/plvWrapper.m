function [plv] = plvWrapper(signal,fs,freq_range)
% This is a function to run on the response timing data collected by David
% Caldwell and Jeneva Cronin while in the GRIDLab. This uses a
%a PLV function developed by Praneeth Namburi . 
% Required inputs are a:
% signal:
%   time x channels x trials
% fs:
%   sampling rate in Hz

%
% OUTPUT 
% PLV - time x channel x channel x numConditions 
%
% DJC 3-30-2017

desired_f = round(freq_range(1)+freq_range(2))/2;
period = 1/desired_f;
time_4oscil = period*4; % time total in seconds
order = round(time_4oscil*fs);

filtSpec.order = order;
filtSpec.range = freq_range;

% need channels x time x trials
signal_perm = permute(signal,[2,1,3]);
[plv] = pn_eegPLV(signal_perm,fs,filtSpec);


end