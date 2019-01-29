function [buttonLocs,buttonLocsSamps,tactorLocsVec,tactorLocsVecSamps,tEpoch,epochedButton,epochedTactor,buttonTactDiffSamps] = get_response_timing_segs_newTactor(tactorData,uniqueCond,stim,buttonData,stimFromFile,fsStim,fsTact,trainTimesTotal,plotIt)

tButton = (0:length(buttonData)-1)/fsTact;
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

tactThresh = 0.97;
tactorData(tactorData > tactThresh) = tactThresh;
% raw button
%[buttonPks,buttonLocs] = findpeaks(buttonData,t_button,'MinpeakDistance',1.5)
%findpeaks(buttonData,t_button,'MinpeakDistance',2,'Minpeakheight',10e-3)
%% QUANTIFY RXN TIME TO CORTICAL STIM

sampsEnd = round(3.5*fsStim);

% epoched button press

epochedButton = {};

for i = 1:length(uniqueCond)
    epochedButton{i} = squeeze(getEpochSignal(buttonDataClip,trainTimesTotal{i},(trainTimesTotal{i} + sampsEnd)));
end

epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesTotal{1},trainTimesTotal{1} + sampsEnd));

tEpoch = [0:size(epochedTactor,1)-1]/fsStim;
tEpochSamps = [0:size(epochedTactor,1)-1];

% include tactor delay
useTactorDelay = 0;

if useTactorDelay
    tactorDelaySamps = (1.04/1e3)*fsStim; % ms
    tactorDelaySecs = 1.04/1e3;
else
    tactorDelaySamps = 0;
    tactorDelaySecs = 0;
end


%buttonPks = {};
buttonLocs = {};

%buttonPksSamps = {};
buttonLocsSamps = {};

%buttonPksTempVec = [];
buttonLocsTempVec = [];

%buttonPksTempVecSamps = [];
buttonLocsTempVecSamps = [];


%%
for i = 1:length(uniqueCond)
    
    % for stimulation condititions
    if uniqueCond(i)~=-1
        buttonStart = 0.1;
        buttonEnd = 2.5;
        
        for j = 1:length(trainTimesTotal{i})
            
            
            [ipt,residual] = findchangepts((epochedButton{i}((tEpoch>buttonStart & tEpoch < buttonEnd),j)),'maxnumchanges',2);
            
            if isempty(ipt) || max((epochedButton{i}((tEpoch>buttonStart & tEpoch < buttonEnd),j))) < 8e-3
                ipt = NaN;
                
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            %  buttonLocsTempSamps = ipt(1)+round(1*fsTact);
            buttonLocsTempSamps = ipt(1)+round(buttonStart*fsTact);
            
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            %            if j == 1
            %                figure
            %            end
            %                         if i == 5
            %                                         figure
            %                                          findchangepts((epochedButton{i}(:,j)),'maxnumchanges',2)
            %                             figure
            %                             plot(tEpoch((tEpoch>0.075 & tEpoch<1.5)),(epochedButton{i}((tEpoch>0.075 & tEpoch<1.5),j)))
            %                             figure
            %                             findchangepts((epochedButton{i}((tEpoch>0.075 & tEpoch<2),j)),'maxnumchanges',2);
            %                         end
            %
            %[%buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton{i}(:,j),tEpoch,'NPeaks',1,'Minpeakheight',0.008);
            %[%buttonPksTempSamps,buttonLocsTempSamps] = findpeaks(epochedButton{i}(:,j),tEpochSamps,'NPeaks',1,'Minpeakheight',0.008); % get sample number DJC 10-12-2017
            sprintf(['button ' num2str(buttonLocsTemp)])
            sprintf(['button ' num2str(buttonLocsTempSamps)])
            
            if isempty(buttonLocsTemp)
                %buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
                %buttonPksTempSamps = NaN;
                buttonLocsTempSamps = NaN;
            end
            %buttonPksTempVec(j) = %buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
            
            % do samples too
            %buttonPksTempVecSamps(j) = %buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
        end
        %buttonPks{i} = %buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        %buttonPksSamps{i} = %buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
        % for tactor target condition
    elseif uniqueCond(i)==-1
        buttonStart = 1;
        for j = 1:length(trainTimesTotal{i})
            % set above certain threshold to 0.009
            
            [ipt,residual] = findchangepts((epochedButton{i}(tEpoch>buttonStart,j)),'maxnumchanges',2);
            if isempty(ipt) || max((epochedButton{i}(tEpoch>1,j))) < 8e-3
                ipt = NaN;
                
            elseif length(ipt) >= 2
                if ipt(2)-ipt(1) < round(fsTact*5/1e3)
                    ipt = NaN;
                end
            end
            
            buttonLocsTempSamps = ipt(1)+round(buttonStart*fsTact);
            buttonLocsTemp = buttonLocsTempSamps/fsTact;
            
            %
            %             figure
            %
            %             findchangepts((epochedButton{i}(tEpoch>buttonStart,j)),'maxnumchanges',2);
            %             figure
            %             findpeaks((epochedTactor(:,j)),tEpochSamps,'NPeaks',1,'Minpeakheight',1);
            %             xlim([24414 length(epochedTactor(:,j))])
            
            
            %[%buttonPksTemp,buttonLocsTemp] = findpeaks((epochedButton{i}(tEpoch>1,j)),tEpoch(tEpoch>1),'NPeaks',1,'Minpeakheight',0.008);
            % [%buttonPksTempSamps,buttonLocsTempSamps] = findpeaks((epochedButton{i}(tEpochSamps>24415,j)),tEpochSamps(tEpochSamps>24415),'NPeaks',1,'Minpeakheight',0.008);
            
            sprintf(['button ' num2str(buttonLocsTemp)])
            sprintf(['button ' num2str(buttonLocsTempSamps)])
            
            [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),tEpoch,'NPeaks',1,'Minpeakheight',tactThresh-0.05);
            [tactorPksTempSamps,tactorLocsTempSamps] = findpeaks((epochedTactor(:,j)),tEpochSamps,'NPeaks',1,'Minpeakheight',tactThresh-0.05);
            
            sprintf(['tactor ' num2str(tactorLocsTemp)])
            sprintf(['tactor ' num2str(tactorLocsTempSamps)])
            
            if isempty(buttonLocsTemp)
                %buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
                %buttonPksTempSamps = NaN;
                buttonLocsTempSamps = NaN;
            end
            
            if isempty(tactorLocsTemp)
                % %tactorPksTemp = NaN;
                tactorLocsTemp = NaN;
                %%tactorPksTempSamps = NaN;
                tactorLocsTempSamps = NaN;
                
            end
            
            %buttonPksTempVec(j) = %buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
            
            % do samples too
            %buttonPksTempVecSamps(j) = %buttonPksTempSamps;
            buttonLocsTempVecSamps(j) = buttonLocsTempSamps;
            
            % account for tactor delay
            if ~useTactorDelay
                tactorLocsVec(j) = tactorLocsTemp;
                tactorLocsVecSamps(j) = tactorLocsTempSamps;
            else
                tactorLocsVec(j) = tactorLocsTemp - tactorDelaySecs;
                tactorLocsVecSamps(j) = tactorLocsTempSamps - tactorDelaySamps;
            end
            
            %tactorPksVec(j) = %tactorPksTemp;
            %tactorPksVecSamps(j) = %tactorPksTempSamps;
            
            
        end
        
        %buttonPks{i} = %buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        %buttonPksSamps{i} = %buttonPksTempVecSamps;
        buttonLocsSamps{i} = buttonLocsTempVecSamps;
        
        
    end
end

%%
% calculate differences - MAKE SURE YOU ONLY DO THIS ONCE

% in case of no tactor trials
if exist('tactorLocsVec','var')
    buttonTactDiff = buttonLocs{uniqueCond==-1} - tactorLocsVec;
    buttonTactDiffSamps = buttonLocsSamps{uniqueCond==-1} - tactorLocsVecSamps;
    
    buttonLocs{uniqueCond==-1} = buttonTactDiff;
    buttonLocsSamps{uniqueCond==-1} = buttonTactDiffSamps;
else
    tactorLocsVec = [];
    tactorLocsVecSamps = [];
    buttonTactDiff = [];
    buttonTactDiffSamps = [];
end

end