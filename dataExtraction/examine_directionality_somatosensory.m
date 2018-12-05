%% script to calculate metrics of directionality from processed somatosensory evoked responses
%
% David.J.Caldwell 12/4/2018
loadIt = 0;
if loadIt
    load('3ada8b_tactor_only_processed_data.mat') 
end
%% phase slope index

freqbins = {[4:8],[8:12],[12:30],[70:170]};
freqbins = {[70:170]};
% parameters for PSI-calculation
segleng=500;epleng=1000;
for trial = 1:size(processedSigTactor,3)
    for freq = freqbins
        freq = freq{:};
        [psi, stdpsi, psisum, stdpsisum]=data2psi(processedSigTactor(:,:,trial),segleng,epleng,freq);
        psi./(stdpsi+eps);
    end
end

%%
figure
(imagesc(psi))