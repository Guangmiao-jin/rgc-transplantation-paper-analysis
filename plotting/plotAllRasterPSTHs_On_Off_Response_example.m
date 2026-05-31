function ax = plotAllRasterPSTHs_On_Off_Response_example(plotMetrics)

titleText = [{'Scotopic'}, {'Mesopic'}, {'Photopic'}];

figH = figure('units','normalized','outerposition',[0 0 1 1], 'Color','white', 'MenuBar','none');
%% plot raster per trial
nBlocks  = length(plotMetrics.trialSpikes);
for stimblk = 1:nBlocks
    if nBlocks == 1
        ax(stimblk) = subplot(2,1,stimblk); hold on
    else
        ax(stimblk) = subplot(2,3,stimblk); hold on
    end
    trialSpikesCnd = plotMetrics.trialSpikes{stimblk};
    xSpikePos = [];
    ySpikePos = [];
    % for each trial
    for tr = 1:length(trialSpikesCnd)
        trSpikes =trialSpikesCnd{tr};


        % build plottings
        xSpikePosTemp = repmat(trSpikes',2,1);
        xSpikePos = [xSpikePos xSpikePosTemp];

        ySpikePosTemp(1,:) = tr-1;                % Y-offset for raster plot
        ySpikePosTemp(2,:) = tr;
        ySpikePosTemp = repmat(ySpikePosTemp, 1, size(xSpikePosTemp,2));

        ySpikePos = [ySpikePos ySpikePosTemp];

        % wipe variables for next trial
        trSpikes = [];
        ySpikePosTemp = [];
        xSpikePosTemp = [];
        % disp(['Len X: ' num2str(length(xSpikePos)) ' Len Y: ' num2str(length(ySpikePos))])
    end

    plot(xSpikePos, ySpikePos, 'Color', 'k');

    % ---- light-on bar above raster ----
    lightWindow = [-2 0];           % 
    nTrials = length(trialSpikesCnd);

    barH   = 2;                  
    gap    = 1;                  
    y0     = nTrials + gap;        

    currAx = ax(stimblk);
    currAx.YLim = [0 nTrials + gap + barH + 0.2];

     patch(currAx, ...
    [lightWindow(1) lightWindow(2) lightWindow(2) lightWindow(1)], ...
    [y0 y0 y0+barH y0+barH], ...
    [1 0.9 0.1], 'EdgeColor','none', 'FaceAlpha',1, ...
    'HandleVisibility','on');

    text(mean(lightWindow), y0 + barH/2, 'light on', ...
    'Parent', currAx, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontWeight','bold', 'Color','k');

    darkWindow = [0 2];           % 

     patch(currAx, ...
    [darkWindow(1) darkWindow(2) darkWindow(2) darkWindow(1)], ...
    [y0 y0 y0+barH y0+barH], ...
    [0.5 0.5 0.5], 'EdgeColor','none', 'FaceAlpha',1, ...
    'HandleVisibility','on');

    text(mean(darkWindow), y0 + barH/2, 'light off', ...
    'Parent', currAx, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontWeight','bold', 'Color','k');
    set(currAx,'Box','on','LineWidth',1.2);


    currAx = ax(stimblk);
    currAx.XLim             = [-2 plotMetrics.binEdges(end)];
    %currAx.YLim             = [0 length(trialSpikesCnd)];

    currAx.XLabel.String  	= 'Time(s)';
    currAx.YLabel.String  	= 'Trials';
    xline(0, 'Color', 'r', 'LineWidth',2);
    title(titleText{stimblk});



    %% PSTH
    if nBlocks == 1
        ax(stimblk+1)   = subplot(2,1,stimblk+1);
    else
        ax(stimblk+3)   = subplot(2,3,stimblk+3);
    end
    binsize = 25; %ms
    % nbins               = (range(responseMetrics.binEdges)*1000)/binsize;                        % Bin duration in [ms]
    nobins              = 1000/binsize;                            % No of bins/sec

    meanPSTH = mean(plotMetrics.trialPSTHs{stimblk});
    countAverageSec     = (meanPSTH) * nobins;


    h                   = histogram('BinCounts', countAverageSec, 'BinEdges', plotMetrics.binEdges);
    h.FaceColor         = 'k';

    hold on
    xline(0, 'Color', 'r', 'LineWidth',2)

    mVal                = max(h.Values)+round(max(h.Values)*.1);
    if nBlocks == 1
        currAx = ax(stimblk+1);
    else
        currAx = ax(stimblk+3);
    end
    currAx.XLim             = [-2 plotMetrics.binEdges(end)];

    % fix for empty histogram
    if mVal == 0
        mVal = 1;
    end

    currAx.YLim             = [0 mVal];
    currAx.XLabel.String  	= 'Time(s)';
    currAx.YLabel.String  	= 'Average spikes per second';

end
if nBlocks ~= 1
    subplotEvenAxes(ax, [0 1 0] , [4 5 6])
end
end