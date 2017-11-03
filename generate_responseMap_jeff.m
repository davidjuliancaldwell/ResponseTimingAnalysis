%% script 7-11-2017 to generate small multiple style plot that jeff was interested in

load('c19968_block_2_postICAprocessData.mat')
sig = subtracted_sig_matrixS_I;
avgResponse = mean(sig,3);

stimChans = [17 9];

%%
chansInt = [1 2 3 9 10 11 17 18 19 25 26 27];

smallMultiples_responseTiming(avgResponse,t,'type1',stimChans,'type2',0,'average',1)


%% 9-28-2017 - DJC 
load('c19968_block_1_postICAprocessData.mat')
%%
sig = processedSig;
avgResponse = mean(sig,3);

stimChans = [17 9];

%%
chansInt = [1 2 3 9 10 11 17 18 19 25 26 27];

smallMultiples_responseTiming(avgResponse,t,'type1',stimChans,'type2',0,'average',1)