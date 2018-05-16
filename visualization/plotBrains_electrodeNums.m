function plotBrains_electrodeNums(subjid,stims)

load(fullfile(getSubjDir(subjid), 'trodes.mat'));

clims = [-1 1];
% w = zeros(size(Grid, 1), 1);

%  w(stims) = -1;
% w(beta) = 1;

figure

%     % original plot
%  PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), clims, 20, 'recon_colormap', 1:size(Grid, 1), true);

%     % to just plot white labels
%map = [1 1 0; 0 0 0; 1 0 1];
% PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), [-1 1], 10, map,[],[]);
% colormap('flag')

% DJC - 1-17-2017 - plot brain with white numbers
%map = [1 1 1];
%map = [1 1 1; 1 1 1; 1 1 1];
map = [.2 1 0; 1 1 1; 1 0 1];

w = zeros(size(Grid, 1), 1);
w(stims) = 1;
clims = [-1 1];
PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), clims, 15, map, 1:size(Grid, 1), true);

end