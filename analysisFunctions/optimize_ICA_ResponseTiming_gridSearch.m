function [best_numComponents,best_nonLinear,outputSig,recon_artifact] = optimize_ICA_ResponseTiming_gridSearch(data,varargin)

scale_factor = 1000;
stimChans = [];
meanSub = 0;
orderPoly = 2;
plotIt = 0;
channelInt = 45;

% lower bound for num components factor
lb = [3];
ub = [size(data,2)];
% starting scale factor
test_points = [lb:1:ub];
nonlinears_cell = {'pow3','tanh','gauss','skew'};


for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'orderpoly'
            orderPoly = varargin{i+1};
        case 'stimchans'
            stimChans = varargin{i+1};
        case 'meansub'
            meanSub = varargin{i+1};
        case 'examptrial'
            exampTrial = varargin{i+1};
        case 'fs'
            fs = varargin{i+1};
        case 'plotit'
            plotIt = varargin{i+1};
        case 'channelint'
            channelInt = varargin{i+1};
            
    end
    
end


matrix_vals = zeros(length(test_points),length(nonlinears_cell));

for i = 1:length(test_points)
    for j = 1:length(nonlinears_cell)
        nonlinear = nonlinears_cell{j};
        numComponentsSearch = test_points(i);
        [evaluated_value,huberInside_mat,huberOutside_mat,processedSig,recon_artifact] =  ica_train_optimize_subtract_ResponseTiming(data,stimChans,fs,numComponentsSearch,scale_factor,meanSub,orderPoly,plotIt,channelInt,nonlinear);
        matrix_vals(i,j) = evaluated_value;
        fprintf(['number of components  ' num2str(i) ' nonlinear ' num2str(j) ' complete \n'])

    end
end

[minMatrix,I] = min(matrix_vals(:));
[row, col] = ind2sub(size(matrix_vals),I);

best_numComponents = test_points(row);
best_nonLinear = nonlinears_cell{col};

[~,~,~,outputSig,recon_artifact] = ica_train_optimize_subtract_ResponseTiming(data,stimChans,fs,best_numComponents,scale_factor,meanSub,orderPoly,plotIt,channelInt,best_nonLinear);

end