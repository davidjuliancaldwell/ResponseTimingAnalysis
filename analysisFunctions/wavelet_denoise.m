function outputSig = wavelet_denoise(inputSig,varargin)

numComponents = 15;
beginRecon = 0;
endRecon = numComponents;
plotIt = 0;
exampChan = 10;
exampTrial = 10;
fs = 1.2207e4;

for i=1:2:(length(varargin)-1)

    switch lower(varargin{i})
        case 'numcomponents'
            numComponents = varargin{i+1};
        case 'beginrecon'
            beginRecon = varargin{i+1};
        case 'endrecon'
            endRecon = varargin{i+1};
        case 'plotit'
            plotIt = varargin{i+1};
        case 'exampchan'
            exampChan = varargin{i+1};
        case 'examptrial'
            exampTrial = varargin{i+1};
                    case 'fs'
            fs = varargin{i+1};
    end
    
end

outputSig = zeros(size(inputSig));
for ind_chan = 1:size(inputSig,2)
    string = ['Doing analysis for channel ' num2str(ind_chan) ' '];
    textprogressbar(string)
    num_iter = size(inputSig,3);
    
    for ind_trial = 1:size(inputSig,3)
        
        textprogressbar((ind_trial/num_iter)*100);

        tempSig = squeeze(inputSig(:,ind_chan,ind_trial));
        wt = modwpt(tempSig,numComponents);
        wtrec = zeros(size(wt));
        wtrec(beginRecon:endRecon,:) = wt(beginRecon:endRecon,:);
        outputSig(:,ind_chan,ind_trial) = imodwpt(wtrec,'sym4');
    end
    
        textprogressbar(' Finished channel');

end

if plotIt
    figure
    orig = inputSig(:,exampChan,exampTrial);
    plot(orig,'linewidth',2)
    hold on
    new = outputSig(:,exampChan,exampTrial);
    plot(new,'linewidth',2)
    legend({'Original Signal','Denoised Signal'});
    ylim([-3e-4 3e-4])
    
    figure
    [f,P1] = spectralAnalysisComp(fs,squeeze(outputSig(:,exampChan,exampTrial)));
    plot(f,P1)
    figure
   plot(f,P1);

    set(gca,'xscale','log');
       set(gca,'yscale','log');

end


end

