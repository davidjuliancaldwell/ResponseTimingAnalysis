function [powerout,f,coi] = wavelet_wrapper(signal,fs,badChans)
% This is a function to run on the response timing data collected by David
% Caldwell and Jeneva Cronin while in the GRIDLab. This uses a
% morletprocess script as implemented by James Wu
% signal:
%   time x channels x trials
% fs:
%   sampling rate in Hz
% time_res:
%   time resolution for morlet process
% badChans:
%   channels to ignore
%
% OUTPUT
% powerout - freq x time x channel x trial

%
% DJC 3-30-2017

% Morlet Process
%%%%%%%%%%%%%%%%%%%
numTrials = size(signal,3);
numChannels = size(signal,4);
timeRes = 0.050;

if ~exist('badChans','var')
    badChans = [];
end

badChansMask = ones(size(signal,2),1);
badChansMask(badChans) = 0;
badChansMask = logical(badChansMask);
signalTemp = signal(:,badChansMask,:);
freqLimits = [1 300];

totalChans = size(signal,2);

for i = 1:numTrials
    
    for j = 1:numChannels
        dataTempChan = squeeze(signalTemp(:,j,i));
        [powerOut_temp,f,coi] = cwt(dataTempChan,fs,'amor','frequencyLimits',freqLimits);
        [minf,maxf] = cwtfreqbounds(numel( dataTempChan ),fs,'wavelet','amor');
        numfreq = 10;
        freq = logspace(log10(minf/1e3),log10(maxf/1e3),numfreq);
        AX = gca;
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
        
        [powerout_temp2, f, t] = morletprocess(dataTempChan, fs, timeRes);
        poweroutTemp(:,:,:,i) = powerout_temp;
    end
    % [powerout_temp] = morletprocess( data_temp, fs, time_res);
    poweroutTemp(:,:,:,i) = powerout_temp;
    
    
end

powerout(:,:,badChansMask,:) = poweroutTemp;
powerout(:,:,~badChansMask,:) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set phase angle output to be zero

if ~exist('phase_angle','var')
    phase_angle = [];
end

end