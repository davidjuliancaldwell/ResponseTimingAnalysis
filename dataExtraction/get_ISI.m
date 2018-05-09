function [buttonLocs,buttonLocsSamps,tactorLocsVec,tactorLocsVecSamps,tEpoch] = get_ISI(tactorData,uniqueCond,stim,buttonData,stimFromFile,fsStim,fsTact,trainTimesTotal,plotIt)

tButton = (0:length(buttonData)-1)/fsTact;

% set above certain threshold to 0.009
buttonDataClip = buttonData;
buttonDataClip(buttonData >= 0.009) = 0.009;
[buttonPks,buttonLocs] = findpeaks(buttonDataClip,tButton,'MinpeakDistance',2,'Minpeakheight',0.008);

if plotIt
    figure
    findpeaks(buttonDataClip,tButton,'MinpeakDistance',2,'Minpeakheight',0.008);
    hold on
    plot(tButton,stimFromFile,'g');
    plot(tButton,stim(:,1),'r');
    plot(tButton,stim(:,2),'m');
    
    legend({'Button Data','Button Press Onset Peaks','Stimulation Times From File','S1 stim output','Off Target Stim Output'})
end

% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)
%% QUANTIFY RXN TIME TO CORTICAL STIM

sampsEnd = round(3.5*fsStim);

% epoched button press

epochedButton = {};

% the last trial of the epoched button press for the null condition is
% clipped - so omit that one

trainTimesTotal{2}(end) = [];

for i = 1:length(uniqueCond)
    epochedButton{i} = squeeze(getEpochSignal(buttonDataClip,trainTimesTotal{i},(trainTimesTotal{i} + sampsEnd)));
end

epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesTotal{1},trainTimesTotal{1} + sampsEnd));

tEpoch = [0:size(epochedTactor,1)-1]/fsStim;
tEpochSamps = [0:size(epochedTactor,1)-1];


buttonPks = {};
buttonLocs = {};

buttonPksSamps = {};
buttonLocsSamps = {};

buttonPksTempVec = [];
buttonLocsTempVec = [];

buttonPksTempVecSamps = [];
buttonLocsTempVecSamps = [];
%%
for i = 1:length(uniqueCond)
    
    % for stimulation condititions
    if uniqueCond(i)~=-1
        for j = 1:length(trainTimesTotal{i})
            [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton{i}(:,j),tEpoch,'NPeaks',1,'Minpeakheight',0.008);
            [buttonPksTempSamps,buttonLocsTempSamps] = findpeaks(epochedButton{i}(:,j),tEpochSamps,'NPeaks',1,'Minpeakheight',0.008); % get sample number DJC 10-12-2017
            
            if isempty(buttonPksTemp)
                buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
                buttonPksTempSamps = NaN;
                buttonLocsTempSamps = NaN;
            end
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
            
            % do samples too
            buttonPksTempVecSamps(j) = buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
        end
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        buttonPksSamps{i} = buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
        % for tactor target condition
    elseif uniqueCond(i)==-1
        for j = 1:length(trainTimesTotal{i})
            
            [buttonPksTemp,buttonLocsTemp] = findpeaks((epochedButton{i}(tEpoch>1,j)),tEpoch(tEpoch>1),'NPeaks',1,'Minpeakheight',0.008);
            [buttonPksTempSamps,buttonLocsTempSamps] = findpeaks((epochedButton{i}(tEpochSamps>24415,j)),tEpochSamps(tEpochSamps>24415),'NPeaks',1,'Minpeakheight',0.008);
            
            sprintf(['button ' num2str(buttonLocsTemp)])
            sprintf(['button ' num2str(buttonLocsTempSamps)])
            
            [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),tEpoch,'NPeaks',1,'Minpeakheight',1);
            [tactorPksTempSamps,tactorLocsTempSamps] = findpeaks((epochedTactor(:,j)),tEpochSamps,'NPeaks',1,'Minpeakheight',1);
            
            sprintf(['tactor ' num2str(tactorLocsTemp)])
            sprintf(['tactor ' num2str(tactorLocsTempSamps)])
            
            if isempty(buttonPksTemp)
                buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
                buttonPksTempSamps = NaN;
                buttonLocsTempSamps = NaN;
            end
            
            if isempty(tactorPksTemp)
                tactorPksTemp = NaN;
                tactorLocsTemp = NaN;
                tactorPksTempSamps = NaN;
                tactorLocsTempSamps = NaN;
                
            end
            
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
            
            % do samples too
            buttonPksTempVecSamps(j) = buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
            
            tactorPksVec(j) = tactorPksTemp;
            tactorLocsVec(j) = tactorLocsTemp;
            
            tactorPksVecSamps(j) = tactorPksTempSamps;
            tactorLocsVecSamps(j) = tactorLocsTempSamps;
        end
        
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        buttonPksSamps{i} = buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
    end
end

%%
% calculate differences - MAKE SURE YOU ONLY DO THIS ONCE

buttonTactDiff = buttonLocs{uniqueCond==-1} - tactorLocsVec;
buttonTactDiffSamps = buttonLocsSamps{uniqueCond==-1} - tactorLocsVecSamps;

buttonLocs{uniqueCond==-1} = buttonTactDiff;
buttonLocsSamps{uniqueCond==-1} = buttonTactDiffSamps;

end