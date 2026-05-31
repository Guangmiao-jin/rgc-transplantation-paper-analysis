function PlotDeltaRatePerCondition(T, outputFolder)

T_plot = T;

% exclude unresponsive cells
T_plot = T_plot(string(T_plot.CellType) ~= "unresponsive", :);

stageOrder = ["early degeneration", "late degeneration"];
groupOrder = ["SHAM", "NRL"];

% =====================================================
% all cell type combinations
% =====================================================
cellTypes = ["ON", "OFF"];  % original individual types

% All cell types combined（not including unresponsive）
% ON + OFF + ONOFF + unconventional
allCombinedTypes = ["ON", "OFF", "ONOFF", "unconventional"];

plotGroupOrder = [
    "early degeneration_SHAM"
    "early degeneration_NRL"
    "late degeneration_SHAM"
    "late degeneration_NRL"
];

if iscategorical(T_plot.LightCondition)
    conditions = categories(T_plot.LightCondition);
else
    conditions = unique(T_plot.LightCondition);
end

% ======================================================
% PART 1: original per condition + per cell type plots
% （ON cells / OFF cells，separated by light condition）
% ======================================================
for c = 1:numel(conditions)
    cond = conditions(c);
    
    for ct = 1:numel(cellTypes)
        cellType = cellTypes(ct);
        
        if iscategorical(T_plot.LightCondition)
            idxCond = T_plot.LightCondition == categorical(cond);
            condLabel = string(cond);
        else
            idxCond = T_plot.LightCondition == cond;
            condLabel = string(cond);
        end
        
        idxCell = string(T_plot.CellType) == cellType;
        S = T_plot(idxCond & idxCell, :);
        S = addPlotGroup(S, plotGroupOrder, T_plot, cellType, cond);
        
        fig = makePlot(S, plotGroupOrder, ...
            sprintf('%s cells, light condition %s', cellType, condLabel));
        
        saveas(fig, fullfile(outputFolder, ...
            sprintf('%s_cells_delta_rate_condition_%s.png', ...
            cellType, condLabel)));
        close(fig);
    end
end

% ======================================================
% PART 2: original three conditions combined，per cell type
% （ON cells / OFF cells，three light conditions combined）
% ======================================================
for ct = 1:numel(cellTypes)
    cellType = cellTypes(ct);
    
    idxCell = string(T_plot.CellType) == cellType;
    S1 = T_plot(idxCell, :);  % no filter condition, three combined
    S1 = addPlotGroup(S1, plotGroupOrder, T_plot, cellType, "");
    
    fig = makePlot(S1, plotGroupOrder, ...
        sprintf('%s cells, three light conditions combined', cellType));
    
    saveas(fig, fullfile(outputFolder, ...
        sprintf('%s_cells_delta_rate_three_conditions_combined.png', ...
        cellType)));
    close(fig);
end

% ======================================================
% PART 3: All cell types combined
% per condition
% ======================================================
for c = 1:numel(conditions)
    cond = conditions(c);
    
    if iscategorical(T_plot.LightCondition)
        idxCond = T_plot.LightCondition == categorical(cond);
        condLabel = string(cond);
    else
        idxCond = T_plot.LightCondition == cond;
        condLabel = string(cond);
    end
    
    % contain all non-unresponsive types
    idxCell = ismember(string(T_plot.CellType), allCombinedTypes);
    S_all = T_plot(idxCond & idxCell, :);
    S_all = addPlotGroup(S_all, plotGroupOrder, T_plot, "ALL", cond);
    
    fig = makePlot(S_all, plotGroupOrder, ...
        sprintf('All responsive cells, light condition %s', condLabel));
    
    saveas(fig, fullfile(outputFolder, ...
        sprintf('ALL_cells_delta_rate_condition_%s.png', condLabel)));
    close(fig);
end

% ======================================================
% PART 4: All cell types combined
% three conditions combined
% ======================================================
idxCell = ismember(string(T_plot.CellType), allCombinedTypes);
S_all_combined = T_plot(idxCell, :);
S_all_combined = addPlotGroup(S_all_combined, plotGroupOrder, ...
                              T_plot, "ALL", "");

fig = makePlot(S_all_combined, plotGroupOrder, ...
    'All responsive cells, three light conditions combined');

saveas(fig, fullfile(outputFolder, ...
    'ALL_cells_delta_rate_three_conditions_combined.png'));
close(fig);

fprintf('\nDone! All figures saved to: %s\n', outputFolder);

end

% ======================================================
% Helper function 1: addPlotGroup
% PlotGroup + NaN placeholder
% ======================================================
function S = addPlotGroup(S, plotGroupOrder, T_ref, cellType, cond)

    if height(S) == 0
        % if no data available
        warning('No data for cellType=%s, cond=%s', cellType, cond);
        return;
    end
    
    S.PlotGroup = strcat(string(S.Stage), "_", string(S.Group));
    S.PlotGroup = categorical(S.PlotGroup, plotGroupOrder, 'Ordinal', true);
    S.Xpos = double(S.PlotGroup);
    
    % fill in NaN placeholder for each group
    for pg = 1:numel(plotGroupOrder)
        thisGroup = plotGroupOrder(pg);
        if ~any(string(S.PlotGroup) == thisGroup)
            newRow = S(1, :);
            newRow.DeltaFR = NaN;
            
            if startsWith(thisGroup, "early")
                newRow.Stage = categorical("early degeneration");
            else
                newRow.Stage = categorical("late degeneration");
            end
            
            if endsWith(thisGroup, "SHAM")
                newRow.Group = categorical("SHAM");
            else
                newRow.Group = categorical("NRL");
            end
            
            newRow.PlotGroup = categorical(thisGroup, ...
                plotGroupOrder, 'Ordinal', true);
            newRow.Xpos = pg;
            S = [S; newRow];
        end
    end
end

% ======================================================
% Helper function 2: makePlot
% plotting code
% ======================================================
function fig = makePlot(S, plotGroupOrder, titleStr)

    fig = figure;
    hold on;
    
    % Boxchart + swarmchart
    boxchart(S.Xpos, S.DeltaFR);
    swarmchart(S.Xpos, S.DeltaFR, 35, 'filled');
    
    ylabel('\Delta firing rate (spikes/s)');
    title(titleStr);
    xticks(1:numel(plotGroupOrder));
    xticklabels({'early SHAM', 'early NRL', 'late SHAM', 'late NRL'});
    xtickangle(45);
    
    % Y range
    yVals = S.DeltaFR(~isnan(S.DeltaFR));
    if isempty(yVals)
        yTop = 1; yBottom = 0;
    else
        yTop = max(yVals); yBottom = min(yVals);
    end
    yRange = yTop - yBottom;
    if yRange == 0; yRange = 1; end
    
    % N labels
    for pg = 1:numel(plotGroupOrder)
        thisGroup = plotGroupOrder(pg);
        nThis = sum(string(S.PlotGroup) == thisGroup & ~isnan(S.DeltaFR));
        text(pg, yTop + 0.08*yRange, sprintf('n=%d', nThis), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 10);
    end
    
    % Stats
    comparisons = [
        "early degeneration_SHAM", "early degeneration_NRL"
        "late degeneration_SHAM",  "late degeneration_NRL"
        "early degeneration_SHAM", "late degeneration_NRL"
    ];
    
    yStar = yTop + 0.15 * yRange;
    
    for comp = 1:size(comparisons, 1)
        g1 = comparisons(comp, 1);
        g2 = comparisons(comp, 2);
        x1 = find(plotGroupOrder == g1);
        x2 = find(plotGroupOrder == g2);
        
        data1 = S.DeltaFR(string(S.PlotGroup) == g1);
        data2 = S.DeltaFR(string(S.PlotGroup) == g2);
        data1 = data1(~isnan(data1));
        data2 = data2(~isnan(data2));
        
        if numel(data1) >= 2 && numel(data2) >= 2
            p = ranksum(data1, data2);
            starText = pToStars(p);
            
            plot([x1 x1 x2 x2], ...
                [yStar yStar+0.04*yRange yStar+0.04*yRange yStar], ...
                'k-', 'LineWidth', 1);
            text(mean([x1 x2]), yStar+0.06*yRange, starText, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 12);
            
            fprintf('%s | %s vs %s: p = %.4g\n', titleStr, g1, g2, p);
            yStar = yStar + 0.18 * yRange;
        else
            fprintf('%s | %s vs %s: not tested (n1=%d, n2=%d)\n', ...
                titleStr, g1, g2, numel(data1), numel(data2));
        end
    end
    
    hold off;
end

% ======================================================
% pToStars: place star
% ======================================================
function stars = pToStars(p)
    if p < 0.001
        stars = "***";
    elseif p < 0.01
        stars = "**";
    elseif p < 0.05
        stars = "*";
    else
        stars = "ns";
    end
end