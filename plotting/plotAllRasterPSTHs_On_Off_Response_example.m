function plotAllRasterPSTHs_combined_Response_example_excel(psthMatFile, neuronIndex, outputExcelPath)
% Usage:
%   plotAndExportPSTH_example('path/to/file_psth.mat', 42, 'output.xlsx')
%
% Reads pltcurve from .mat, extracts neuron n, 
% plots PSTH and exports raster + PSTH data to Excel.

    % =====================================================
    % 1. Read .mat and build plotMetrics
    % =====================================================
    data = load(psthMatFile);
    pltcurve = data.pltcurve;

    plotMetrics.trialPSTHs = pltcurve.trialPSTHs{neuronIndex};
    plotMetrics.binEdges   = pltcurve.PSTH_binEdges;

    % Unwrap binEdges if needed
    if isstruct(plotMetrics.binEdges)
        fns = fieldnames(plotMetrics.binEdges);
        plotMetrics.binEdges = plotMetrics.binEdges.(fns{1});
    end
    if iscell(plotMetrics.binEdges)
        plotMetrics.binEdges = plotMetrics.binEdges{1};
    end
    plotMetrics.binEdges = plotMetrics.binEdges(:)';

    % =====================================================
    % 2. Plot (your existing function)
    % =====================================================
    plotAllRasterPSTHs_combined_Response_example1(plotMetrics);

    % =====================================================
    % 3. Export to Excel
    % =====================================================
    binEdges = plotMetrics.binEdges(:)';
    binCtrs  = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    binW     = median(diff(binEdges));
    nBins    = numel(binCtrs);

    condLabels = {'Scotopic', 'Mesopic', 'Photopic'};
    nConds     = min(3, length(plotMetrics.trialPSTHs));

    % --- Sheet 1: PSTH summary (mean firing rate per condition) ---
    T_psth = table();
    T_psth.BinCenter_s = binCtrs(:);

    for k = 1:nConds
        psthMat = plotMetrics.trialPSTHs{1, k};
        meanFR  = mean(psthMat, 1, 'omitnan') / binW;
        T_psth.(sprintf('%s_MeanFR', condLabels{k})) = meanFR(:);
    end

    writetable(T_psth, outputExcelPath, 'Sheet', 'PSTH_Summary');

    % --- Sheet 2: Metadata ---
    T_meta = table();
    T_meta.Field = {
        'Source file'; 
        'Neuron index'; 
        'N conditions'; 
        'N bins'; 
        'Bin width (s)';
        'Time range (s)';
        'N trials (Scotopic)';
        'N trials (Mesopic)';
        'N trials (Photopic)'
    };

    nTrials_per_cond = cell(3, 1);
    for k = 1:nConds
        nTrials_per_cond{k} = num2str(size(plotMetrics.trialPSTHs{1,k}, 1));
    end
    for k = (nConds+1):3
        nTrials_per_cond{k} = 'N/A';
    end

    T_meta.Value = {
        psthMatFile;
        num2str(neuronIndex);
        num2str(nConds);
        num2str(nBins);
        num2str(binW);
        sprintf('%.2f to %.2f', binEdges(1), binEdges(end));
        nTrials_per_cond{1};
        nTrials_per_cond{2};
        nTrials_per_cond{3}
    };

    writetable(T_meta, outputExcelPath, 'Sheet', 'Metadata');

    fprintf('\nExcel exported to: %s\n', outputExcelPath);
    fprintf('Sheets:\n');
    fprintf('  1. PSTH_Summary  (bin centers + mean FR per condition)\n');
    fprintf('  2. Metadata\n');
end

function plotAllRasterPSTHs_combined_Response_example1(plotMetrics)
% Plot 3 conditions PSTH lines on ONE axes,
% and add a separate bottom stimulus pattern axis using patch.
% Input: plotMetrics only.
% -------- settings --------
condLabels = {'Scotopic','Mesopic','Photopic'};   % change if needed
lightWindow = [-2 0];   % light ON
darkWindow  = [0 2];    % light OFF
YFIX = [-2 60];        %y-axis
% smoothing (optional)
smooth_ms = 10;        % set [] or 0 to disable
% -------------------------
nBlocks = length(plotMetrics.trialPSTHs);
nShow   = min(3, nBlocks);
binEdges = plotMetrics.binEdges(:)';
binCtrs  = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binW     = median(diff(binEdges)); % seconds
if isempty(binW) || binW <= 0
    error('plotMetrics.binEdges invalid.');
end
% ---- compute rate traces ----
rateAll = cell(1, nShow);
% smoothing window in bins
if ~isempty(smooth_ms) && smooth_ms > 0
    winBins = max(3, round(smooth_ms / (binW*1000)));
    if mod(winBins,2)==0, winBins = winBins+1; end
else
    winBins = 0;
end
for k = 1:nShow
    psthMat = plotMetrics.trialPSTHs{1,k};     % [nTrials x nBins]
    meanCountsPerBin = mean(psthMat, 1);
    y = meanCountsPerBin / binW;             % spikes/s
    % optional smoothing (gaussian)
    if winBins > 0
        if exist('smoothdata','file') == 2
            y = smoothdata(y, 'gaussian', winBins);
        else
            % fallback gaussian conv
            sigma = winBins/6;
            xx = (-floor(winBins/2):floor(winBins/2));
            g = exp(-(xx.^2)/(2*sigma^2)); g = g/sum(g);
            y = conv(y, g, 'same');
        end
    end
    rateAll{k} = y;
end
% ---- figure layout: main axis + bottom stim axis ----
figure('Color','white','MenuBar','none',...
    'Units','normalized','OuterPosition',[0.05 0.1 0.65 0.75]);
left   = 0.12;
width  = 0.82;
bottom = 0.12;
top    = 0.92;
gap    = 0.03;
hBar   = 0.08;
hMain  = (top - bottom - hBar - gap);
% Main axes (all 3 lines)
axMain = axes('Position',[left, bottom + hBar + gap, width, hMain]); 
hold(axMain,'on');
% plot lines (same color black but different line styles to distinguish)
lineStyles = {':','--','-'};  
for k = 1:nShow
    plot(axMain, binCtrs, rateAll{k}, 'k', 'LineWidth', 2, 'LineStyle', lineStyles{k});
end
xline(axMain, 0, 'Color',[0.45 0.45 0.45], 'LineWidth', 1.2);
axMain.XLim = [lightWindow(1) darkWindow(2)];
axMain.YLim = YFIX;
hLen = 0.25;           % horizontal scale bar length (s)
vLen = 5;            % vertical scale bar length (Hz) 
xL = axMain.XLim;
yL = axMain.YLim;
% margins (fraction of range)
mx = 0.06 * range(xL);
my = 0.08 * range(yL);
% anchor (bottom-left of the scale bar)
x0 = xL(2) - mx - hLen;
y0 = yL(2) - my - vLen;
hold(axMain,'on');
% vertical bar
plot(axMain, [x0 x0], [y0 y0+vLen], 'k', 'LineWidth', 2);
% horizontal bar
plot(axMain, [x0 x0+hLen], [y0 y0], 'k', 'LineWidth', 2);
% labels
text(axMain, x0 + hLen/2, y0 - 0.02*range(yL), sprintf('%.2f s', hLen), ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontWeight','bold', 'Color','k');
text(axMain, x0 - 0.02*range(xL), y0 + vLen/2, sprintf('%g Hz', vLen), ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
    'FontWeight','bold', 'Color','k');
axMain.Box = 'off';
axMain.LineWidth = 1.2;
axMain.TickDir = 'out';
ylabel(axMain, 'Firing rate (spikes/s)');
% legend
legend(axMain, condLabels(1:nShow), 'Location','northeast', 'Box','off');
% hide x tick labels on main axes (since stim bar has its own axis)
axMain.XTickLabel = [];
% Stimulus axis (bottom)
axStim = axes('Position',[left, bottom, width, hBar]);
hold(axStim,'on');
axStim.XLim = [lightWindow(1) darkWindow(2)];
axStim.YLim = [0 1];
axStim.Box = 'off';
axStim.LineWidth = 1.2;
axStim.TickDir = 'out';
axStim.YTick = [];
% If you don't want any x-axis here either, hide it:
% axStim.XTick = []; axStim.XColor = 'none'; xlabel(axStim,'');
% If you want x-label only once, keep:
xlabel(axStim,'Time (s)');
% stimulus pattern: white then black
y0 = 0.25; h = 0.5;
patch(axStim, ...
    [lightWindow(1) lightWindow(2) lightWindow(2) lightWindow(1)], ...
    [y0 y0 y0+h y0+h], ...
    [1 0.9 0.1], 'EdgeColor','k', 'LineWidth', 1);
text(mean(lightWindow), y0 + h/2, 'light on', ...
    'Parent', axStim, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontWeight','bold', 'Color','k');
patch(axStim, ...
    [darkWindow(1) darkWindow(2) darkWindow(2) darkWindow(1)], ...
    [y0 y0 y0+h y0+h], ...
    [0.5 0.5 0.5], 'EdgeColor','k', 'LineWidth', 1);
text(mean(darkWindow), y0 + h/2, 'light off', ...
    'Parent', axStim, 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontWeight','bold', 'Color','k');
xline(axStim, 0, 'Color',[0.45 0.45 0.45], 'LineWidth', 1.2);
% Link x axes so zoom/pan stays aligned
linkaxes([axMain, axStim], 'x');
ax = [axMain, axStim];
end