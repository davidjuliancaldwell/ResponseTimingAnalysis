function [] = visualize_wavelet_channel(powerout,tMorlet,fMorlet,processedSig,tEpoch,dataInt,chanInt,stimTime,response,individual,average)
% set colormap using cbrewer
CT = cbrewer('div','RdBu',11);
% flip it so red is increase, blue is down
CT = flipud(CT);

if individual
    
    for i = 1:size(powerout,4)
        totalFig = figure;
        totalFig.Units = 'inches';
        totalFig.Position = [12.1806 3.4931 6.0833 7.8056];
        subplot(3,1,1);
        surf(1e3*tMorlet,fMorlet,powerout(:,:,chanInt,i),'edgecolor','none');
        view(0,90);
        axis tight;
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        title(['Wavelet decomposition Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        xlim([-200 1000]);
        set(gca,'fontsize',14)
        colormap(CT);
        set_colormap_threshold(gcf, [-0.5 0.5], [-6 6], [.5 .5 .5])
        hold on
        plot3([stimTime(i),stimTime(i)],[0 300],[1000,1000],'r','linewidth',2)
        plot3([1e3*response(i),1e3*response(i)],[0 300],[1000,1000],'g','linewidth',2)
        
        ylim([1 200])
        
        h1 = subplot(3,1,2);
        plot(1e3*tEpoch,1e6*processedSig(:,chanInt,i))
        vline(stimTime(i),'r','stim')
        xlabel('time (ms)');
        ylabel('microvolts')
        title(['Processed Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        vline(1e3*response(i),'g','response')
        ylims = [-(max(abs(1e6*processedSig(:,chanInt,i))) + 10) (max(abs(1e6*processedSig(:,chanInt,i))) + 10)];
        ylim(ylims);
        ylim_h1 = ylims;
        xlim([-200 1000]);
        set(gca,'fontsize',14)
        
        h2 = subplot(3,1,3);
        plot(1e3*tEpoch,1e6*dataInt(:,chanInt,i))
        vline(stimTime(i),'r','stim')
        xlabel('time (ms)');
        ylabel('microvolts')
        title(['Raw Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        vline(1e3*response(i),'g','response')
        ylim(ylim_h1);
        xlim([-200 1000]);
        set(gca,'fontsize',14);
        
        linkaxes([h1,h2],'xy');
        
    end
    
end
% now average
if average
    
    poweroutAvg = mean(squeeze(powerout(:,:,chanInt,:)),3);
    
    totalFig2 = figure;
    totalFig2.Units = 'inches';
    totalFig2.Position = [12.1806 3.4931 6.0833 7.8056];
    subplot(3,1,1);
    surf(1e3*tMorlet,fMorlet,poweroutAvg,'edgecolor','none');
    view(0,90);
    axis tight;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    title(['Wavelet decomposition Channel ' num2str(chanInt)]);
    xlim([-200 1000]);
    set(gca,'fontsize',14)
    colormap(CT);
    set_colormap_threshold(gcf, [-0.5 0.5], [-6 6], [.5 .5 .5])
    hold on
    plot3([mean(stimTime),mean(stimTime)],[0 300],[1000,1000],'r','linewidth',2)
    plot3([1e3*mean(response),1e3*mean(response)],[0 300],[1000,1000],'g','linewidth',2)
    
    ylim([1 200])
    %   colorbar;
    % vline(stimTime(i),'r','stim')
    %   vline(1e3*response(i),'g','response')
    
    % figure
    % helperCWTTimeFreqPlot(squeeze(powerout(:,:,chanInt,i)),1e3*tMorlet,fMorlet,'surf')
    
    %figure;
    h3 = subplot(3,1,2);
    plot(1e3*tEpoch,1e6*mean(squeeze(processedSig(:,chanInt,:)),2))
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Processed Channel ' ]);
    ylims = [-(max(abs(1e6*mean(squeeze(processedSig(:,chanInt,:)),2))) + 10) (max(abs(1e6*mean(squeeze(processedSig(:,chanInt,:)),2))) + 10)];
    ylim(ylims);
    ylim_h1 = ylims;
    xlim([-200 1000]);
    vline(mean(stimTime),'r','stim')
    vline(1e3*mean(response),'g','response')
    
    set(gca,'fontsize',14)
    
    h4 = subplot(3,1,3);
    plot(1e3*tEpoch,1e6*mean(squeeze(dataInt(:,chanInt,:)),2))
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Raw Channel ' num2str(chanInt)]);
    ylim(ylim_h1);
    xlim([-200 1000]);
    set(gca,'fontsize',14);
    vline(mean(stimTime),'r','stim')
    vline(1e3*mean(response),'g','response')
    
    linkaxes([h3,h4],'xy');
end

end