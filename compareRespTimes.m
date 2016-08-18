%% 8-15-2016 - script to compare response times once they've been calculated 

% make sure you have vectors buttonLocsVecCort, buttonTactDiff,
% tactorLocsVecTact 

% set bounds on data, assume rxn time has to be greater than 200 ms and
% less than 1 s

tactorLocsVecTactTrim = tactorLocsVecTact(tactorLocsVecTact>respLo & tactorLocsVecTact<respHi);
buttonTactDiffTrim = buttonTactDiff(buttonTactDiff>respLo & buttonTactDiff<respHi);
buttonLocsVecCortTrim = buttonLocsVecCort(buttonLocsVecCort>respLo & buttonLocsVecCort<respHi );

zTact = zscore(tactorLocsVecTactTrim);
zDiff = zscore(buttonTactDiffTrim);
zCort = zscore(buttonLocsVecCortTrim);

tactor = 1e3.*tactorLocsVecTactTrim(abs(zTact)<3);
difference = 1e3.*buttonTactDiffTrim(abs(zDiff)<3);
cort = 1e3.*buttonLocsVecCortTrim(abs(zCort)<3);

current_direc = pwd;

save(fullfile(current_direc, [sid '_compareResponse.mat']), 'tactor', 'difference', 'cort');


%% BOX PLOT

combinedInfo = cat(1,cort,difference,tactor);
groups = cat(1,2.*ones(length(cort),1),ones(length(difference),1),0.*ones(length(tactor),1));
figure
prettybox(combinedInfo,groups,[1 0 0; 0 0 1; 0 1 0],1,false)
fig1 = gca;

fig1.XTick = [];
legend(findobj(gca,'Tag','Box'),'Cortical Stimulation Response Times','Tactor Response Time','Experimenter Response Time')
ylabel('Response times (ms)')
title('Reaction Times')

%% HISTOGRAM

figure
nbins = 15;
a = histogram(cort,nbins);
a.FaceColor = [0 1 0];
hold on
b = histogram(tactor,nbins);
b.FaceColor = [0 0 1];

hold on
c = histogram(difference,nbins);
c.FaceColor = [1 0 0];

title('Histogram of response Times')
legend({'Cortical Stimulation Response Times','Tactor Response Time','Experimenter Response Time'})
xlabel('Response time (ms)')
ylabel('Count')
