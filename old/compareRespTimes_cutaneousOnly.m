%% 8-15-2016 - script to compare response times once they've been calculated 

% make sure you have vectors buttonLocsVecCort, buttonTactDiff,
% tactorLocsVecTact 

% set bounds on data, assume rxn time has to be greater than 200 ms and
% less than 1 s

current_direc = pwd;

%save(fullfile(current_direc, [sid '_compareResponse.mat']), 'tactor', 'difference', 'cort');


%% BOX PLOT

% change colormap to matlab default lines

colormap lines;
cmap = colormap ;


combinedInfo = [];
groups = [];
colors = [];


for i = 1:length(names)
    combinedInfo = [combinedInfo differenceTotal{i}];
    groups = [(i-1).*ones(length(differenceTotal{i}),1) groups];
    colors(i,:) = cmap(i);
end


figure
prettybox(combinedInfo,groups,colors,1,false)
fig1 = gca;

fig1.XTick = [];
legend(findobj(gca,'Tag','Box'),'Cortical Stimulation Response Times','Tactor Response Time','Experimenter Response Time')
ylabel('Response times (ms)')
title('Reaction Times')

%% HISTOGRAM

figure
nbins = 15;

for j = 1:length(names)
    a = histogram(differenceTotal{j,nbins);
    a.FaceColor = cmap(j,:);
    hold on 
end



title('Histogram of response Times')
legend(names)
xlabel('Response time (ms)')
ylabel('Count')
