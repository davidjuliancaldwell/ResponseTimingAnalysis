function [ISICellSamps,ISICellSeconds,ISICondBefore,ISICellSampsNoNuOt,ISICellSecondsNoNuOt,ISICondBeforeNoNuOt] = get_ISI(condType,uniqueCond,tactorLocsVecSamps,stimFromFile,fsStim,trainTimesTotal,trainTimes)

% the last trial of the epoched button press for the null condition is
% clipped - so omit that one

% CHECK THIS DJC 5-8-2018
%trainTimesTotal{2}(end) = [];


ISICellSamps = {};
ISICellSeconds = {};
ISICondBefore = {};

ISICellSampsNoNuOt = {};
ISICellSecondsNoNuOt = {};
ISICondBeforeNoNuOt = {};

% convert train times to actual delivery of tactor

trainTimes(condType==-1) = trainTimes(condType==-1) + tactorLocsVecSamps;


%%
for i = 1:length(uniqueCond)
    
    if uniqueCond(i)~=-1
        trainTimeTemp = trainTimesTotal{i};
    elseif uniqueCond(i)==-1
        trainTimeTemp = trainTimesTotal{i} + tactorLocsVecSamps';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get values including null and off target
    trainTimesDiff = trainTimeTemp - trainTimes;
    trainTimesDiffNaN = trainTimesDiff;
    trainTimesDiffNaN(trainTimesDiffNaN < 0) = NaN;
    
    [~,I] = min(trainTimesDiffNaN,[],2);
    
    isi = trainTimeTemp - (trainTimes(I))';
    isi(isi<=0) = NaN;
    ISICellSamps{i} = isi';
    ISICellSeconds{i} = isi/fsStim';
    ISICondBefore{i} = condType(I);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get values without null and off target
    trainTimesNoNuOt = trainTimes(condType ~= 0 & condType ~= 1);
    condTypeNoNuOt = condType(condType ~= 0 & condType ~= 1);
    
    trainTimesDiff = trainTimeTemp -  trainTimesNoNuOt;
    trainTimesDiffNaN = trainTimesDiff;
    trainTimesDiffNaN(trainTimesDiffNaN < 0) = NaN;
    
    [~,I] = min(trainTimesDiffNaN,[],2);
    
    isi = trainTimeTemp - (trainTimesNoNuOt(I))';
    isi(isi<=0) = NaN;
    
    ISICellSampsNoNuOt{i} = isi';
    ISICellSecondsNoNuOt{i} = isi/fsStim';
    ISICondBeforeNoNuOt{i} = condTypeNoNuOt(I);
    
end

end