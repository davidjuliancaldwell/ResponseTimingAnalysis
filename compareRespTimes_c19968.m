%% 9-13-2016 - script to compare response times once they've been calculated

% load subject data, need sid still
 block = '1';
load([sid,'_compareResponse_block_',block,'.mat'])
% set bounds on data, assume rxn time has to be greater than 0.150
% less than 1 s
%
% tactorLocsVecTactTrim = tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi);
% buttonTactDiffTrim = buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi);
% buttonLocsVecCortTrim = buttonLocsVecCort(buttonLocsVecCort>respLo & buttonLocsVecCort<respHi );
%
% zTact = zscore(tactorLocsVecTactTrim);
% zDiff = zscore(buttonTactDiffTrim);
% zCort = zscore(buttonLocsVecCortTrim);
%
% tactor = 1e3.*tactorLocsVecTactTrim(abs(zTact)<3);
% difference = 1e3.*buttonTactDiffTrim(abs(zDiff)<3);
% cort = 1e3.*buttonLocsVecCortTrim(abs(zCort)<3);

%

for i = 1:length(uniqueCond)
    % 12-10-2016
    respLo = 0.150;
    respHi = 1;
    
    
    trim = buttonLocs{i};
    trim = trim(trim>respLo & trim<respHi);
    zTrim = zscore(trim);
    buttonLocsThresh{i} = 1e3.*trim(abs(zTrim)<3);
    %buttonLocsThresh{i} = 1e3.*trim;

end


%% Histogram

% set number of bins
nbins = 15;

% make cell array for legends
uniqueCondText = cellstr(num2str(uniqueCond));
uniqueCondText{1} = 'tactor';
uniqueCondText{2} = 'no stimulation';
uniqueCondText{3} = 'off target stimulation';


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
keeps = [1 4 5 6 7];

figure
hold on
leg = {};
for i = keeps
    a = histogram(buttonLocsThresh{i},nbins);
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

% DJC - 12/10/2016
% Exclude 100 ms case 
keeps = [1 4 5 6 7];


colormap lines;
cmap = colormap ;


combinedInfo = [];
groups = [];
colors = [];
leg = {};


j = length(keeps);
k = 1;

for i = keeps
    combinedInfo = [buttonLocsThresh{i} combinedInfo ];
   groups = [(j).*ones(length(buttonLocsThresh{i}),1)' groups];
    colors(k,:) = cmap(k,:);
    leg{end+1} = uniqueCondText{i};
    j = j - 1;
    k = k + 1;
end

figure
prettybox(combinedInfo,groups,colors,1,true)
fig1 = gca;

ylim([140 600])
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

