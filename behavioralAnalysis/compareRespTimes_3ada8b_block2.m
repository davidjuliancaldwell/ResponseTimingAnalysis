%% 9-13-2016 - script to compare response times once they've been calculated

%% load in subject
close all;clear all;clc
% this is from my z_constants
Z_ConstantsStimResponse;

DATA_DIR = 'C:\Users\djcald.CSENETID\Data\Subjects\3ada8b\data\d9\MATLAB_conversions\3ada8b_ResponseTiming';
sid = SIDS{6};

% load subject data, need sid still
block = '2';
load([sid,'_compareResponse_block_',block,'.mat'])
% set bounds on data, assume rxn time has to be greater than 0.150
% if want to change reaction times, do it here! 
 respLo = 0.150;
 resphi = 1;

for i = 1:length(uniqueCond)
    
    trim = buttonLocs{i};
    trim = trim(trim>respLo & trim<respHi);
    zTrim = zscore(trim);
   % buttonLocsThresh{i} = 1e3.*trim(abs(zTrim)<3);
    buttonLocsThresh{i} = buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi);
end


%% Histogram

% set number of bins
nbins = 15;

% make cell array for legends
uniqueCondText = cellstr(num2str(uniqueCond));
uniqueCondText{1} = 'no stim';
uniqueCondText{2} = '2 pulses @ 3000 uA, 38 @ 1250 uA';
uniqueCondText{3} = '2 pulses @ 3000 uA';

%individual  histogram of each condition type
    figure

for i = 1:length(uniqueCond)
    subplot(length(uniqueCond),1,i)
    a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
    title(['Histogram of reaction times for condition ' uniqueCondText{i} ])
    xlabel('Time (ms)')
    ylabel('Count')
    xlim([0 1000])
    a.BinWidth = 10;
    a.Normalization = 'probability';
end
subtitle([' block ', block])

% overall histogram

colormap lines;
cmap = colormap ;

figure
hold on
leg = {};
for i = 1:length(uniqueCond)
    a = histogram(1e3*buttonLocs{i}(buttonLocs{i}>respLo & buttonLocs{i}<respHi),nbins);
    leg{end+1} = uniqueCondText{i};
    a.FaceColor = cmap(i,:);
    a.BinWidth = 10;
        a.Normalization = 'probability';
    xlim([0 800])

end
xlabel('Time (ms)')
ylabel('Count')
title(['Histogram of response Times for block ', block])

legend(leg)

% histogram just for tactor and 100,200,400,800

% keep
keeps = [2 3];

figure
hold on
leg = {};
for i = keeps
    a = histogram(1e3*buttonLocsThresh{i},nbins);
    leg{end+1} = uniqueCondText{i};
    a.FaceColor = cmap(i,:);
    a.BinWidth = 10;
        a.Normalization = 'probability';
    xlim([0 500])

end
xlabel('Time (ms)')
ylabel('Count')
title(['Histogram of response Times for block ', block])

legend(leg)

%% BOX PLOT

% change colormap to matlab default lines

colormap lines;
cmap = colormap ;


combinedInfo = [];
groups = [];
colors = [];
leg = {};

keeps = [2 3];


j = length(keeps);
k = 1;

for i = keeps
    combinedInfo = [1e3*buttonLocsThresh{i} combinedInfo ];
   groups = [(j).*ones(length(buttonLocsThresh{i}),1)' groups];
    colors(k,:) = cmap(k,:);
    leg{end+1} = uniqueCondText{i};
    j = j - 1;
    k = k + 1;
end

figure
prettybox(combinedInfo,groups,colors,1,true)
fig1 = gca;

ylim([0 1000])
fig1.XTick = [];
legend(findobj(gca,'Tag','Box'),leg)
ylabel('Response times (ms)')
title(['Reaction Times for block ',block])

%% statistics! kruskal wallis to start 

groupsKW = [];
for i = keeps
    groupsKW = [cellstr(repmat(uniqueCondText{i},[length(buttonLocsThresh{i}),1])); groupsKW(:)];
end

[p,table,stats] = kruskalwallis(combinedInfo,groupsKW);
[c,m,h,nms] = multcompare(stats);

[nms num2cell(m)]

c((c(:,6)<0.05),[1 2 6])

