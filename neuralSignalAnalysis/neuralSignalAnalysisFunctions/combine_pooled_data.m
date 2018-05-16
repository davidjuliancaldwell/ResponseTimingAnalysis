function [buttonLocsSampsCell,buttonLocsCell,data,tEpoch,uniqueCond] = combine_pooled_data(sid,epochedCortEco_cell,tEpoch)

block = [1,2];
buttonLocsSampsCellInd = {};
buttonlocsSamps_cell = {};
data = {};
tEpochGood = tEpoch;

for i = block
    load([sid,'_compareResponse_block_',num2str(i),'.mat'])
    
    buttonLocsSampsCellInd{i} = buttonLocsSamps; % samples
    buttonLocsCellInd{i} = buttonLocs; % seconds
end

numConditions = length(buttonLocs);

for i = 1:numConditions
    buttonLocsSampsCell{i} = [[buttonLocsSampsCellInd{1}{i}] [buttonLocsSampsCellInd{2}{i}]];
    buttonLocsCell{i} =  [[buttonLocsCellInd{1}{i}] [buttonLocsCellInd{2}{i}]];
    data{i} = cat(3,[epochedCortEco_cell{1}{i}], [epochedCortEco_cell{2}{i}]);
    
end

tEpoch = tEpochGood;

end