function [buttonLocsSampsCell,buttonLocsCell,data,tEpoch,uniqueCond] = combine_pooled_data_singleBlock(sid,epochedCortEco_cell,tEpoch,fileName)

buttonLocsSampsCellInd = {};
buttonlocsSamps_cell = {};
data = {};
tEpochGood = tEpoch;

    load(fileName)
    
    buttonLocsSampsCellInd = buttonLocsSamps; % samples
    buttonLocsCellInd = buttonLocs; % seconds

numConditions = length(buttonLocs);

for i = 1:numConditions
    buttonLocsSampsCell{i} = [buttonLocsSampsCellInd{i}] ;
    buttonLocsCell{i} =  [buttonLocsCellInd{i}];
    data{i} = cat(3,[epochedCortEco_cell{i}]);
    
end

tEpoch = tEpochGood;

end