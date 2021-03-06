%% - 3-30-2016 - DJC - response timing file generator - updated 8-19-2016 
% this script will generate two text files. One of these has the sample
% number at which each stimulus train will be delivered. the second has the
% condition which should be read in 

% this is a dummy condition 

prompt = {'Enter subject name','What is the range of ITI?', 'What is the sample rate of the TDT?','Which file number is this?'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'dummy','[1.5,2.5]','24414','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
sid = answer{1};
ITI = str2num(answer{2});
fs = str2num(answer{3});
fileNum = answer{4};


%% make the timing file
% add 1 to ITI - changed from 0.8 which we originally thought. This is to
% account for 200 ms of pulse, then the 800 ms of downtime, then the ITI 

ITIlo = ITI(1)+1;
ITIhi = ITI(2)+1;

% number of trials

numTrials = 140;
randTimes = unifrnd(ITIlo,ITIhi,numTrials,1);


% here the vector is converted to the sample number where the stimulus
% train should start to be delivered 
sample = 1; % start with sample 1 
pts = [];
for i = 1:length(randTimes)
    sample = floor(sample + randTimes(i)*fs);
    pts = [pts; sample];
end

pts1 = pts;
clear pts

%% do timing again
randTimes = unifrnd(ITIlo,ITIhi,numTrials,1);

% here the vector is converted to the sample number where the stimulus
% train should start to be delivered 
sample = 1; % start with sample 1 
pts = [];
for i = 1:length(randTimes)
    sample = floor(sample + randTimes(i)*fs);
    pts = [pts; sample];
end

pts2 = pts;

%% make the conditions file

% tactor 
tact = repmat(-1,20,1);

% no stim
noStim = repmat(0,10,1);

% stim conditions

%numbers of stimuli to delivery 
numEachStim = 20;

stim2 = repmat(200,numEachStim,1);



vectorCond = [stim2];
vectorTact = [tact;noStim];

vectorCondRand1 = vectorCond(randperm(length(vectorCond)));
vectorTactRand1 = vectorTact(randperm(length(vectorTact)));

%% do that again, and stack them

vectorCondRand2 = vectorCond(randperm(length(vectorCond)));
vectorTactRand2 = vectorTact(randperm(length(vectorTact)));

%% put them together

ptsTotal = [pts1; pts2];
vectorCondRandTotal = [vectorCondRand1; vectorTactRand1 ;vectorCondRand2 ;vectorTactRand2];
vectorCondRandTotal = vectorCondRandTotal(randperm(length(vectorCondRandTotal)));

%% write these times to file for stim train delivery

filename = sprintf('%s_stimTrainDelivery_%s.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',pts);
fclose(fileID);

%% write these times to file for condition 

filename = sprintf('%s_condition_%s.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',vectorCondRandTotal);
fclose(fileID);