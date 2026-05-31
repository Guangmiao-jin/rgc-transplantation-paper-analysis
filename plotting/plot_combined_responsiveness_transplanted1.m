function outTbl = plot_combined_responsiveness_transplanted1(baseFolder)

stages   = {'early degeneration', 'late degeneration'};
behavPos1 = 'behavioural positive and ephys positive';
behavPos2 = 'behavioural positive and ephys negative';
behavNeg  = 'behavioural negative';

% show in order（9  positions，and 2 gaps）
plotGroupNames = { ...
    'NRL behav+ ephys+', ...
    'NRL behav+ ephys-', ...
    'NRL behav-', ...
    'Spacer1', ...
    'CRX behav+ ephys+', ...
    'CRX behav+ ephys-', ...
    'CRX behav-', ...
    'Spacer2', ...
    'SHAM' ...
};

% correspond to plotGroupNames（spacer: {}）
scanSpec = { ...
    {behavPos1,'NRL'}, ...
    {behavPos2,'NRL'}, ...
    {behavNeg ,'NRL'}, ...
    {}, ...
    {behavPos1,'CRX'}, ...
    {behavPos2,'CRX'}, ...
    {behavNeg ,'CRX'}, ...
    {}, ...
    {'SHAM'} ...
};

outRows = [];

for s = 1:numel(stages)
    stageName = stages{s};
    stagePath = fullfile(baseFolder, stageName);

    nG = numel(plotGroupNames);
    groups = string(plotGroupNames(:));
    isSpacer = startsWith(groups, "Spacer");

    % save results (for plotting + tables)
    res = struct();
    for i = 1:nG
        res(i).Stage = stageName;
        res(i).Group = plotGroupNames{i};
        res(i).nFiles = 0;
        res(i).bin_centers = [];
        res(i).mean_resp = [];
        res(i).smoothed_mean_resp = [];
        res(i).std_resp = [];
    end

    for i = 1:nG
        if isempty(scanSpec{i})
            continue; % spacer
        end
        spec = scanSpec{i};

        % path：
        % SHAM:  path/SHAM
        % others: path/{NRL|CRX}/{behav folder}
        if isscalar(spec)
            groupPath = fullfile(stagePath, spec{1});
        else
            behav = spec{1};
            treat = spec{2};
            groupPath = fullfile(stagePath, treat, behav);
        end

        if ~isfolder(groupPath)
            continue;
        end

        % searching _distance.mat in an ascending order
        matFiles = dir(fullfile(groupPath, '**', '*_distance.mat'));
        if isempty(matFiles)
            continue;
        end

        Data = [];  % struct array of analysisResults

        for j = 1:length(matFiles)
            filePath = fullfile(matFiles(j).folder, matFiles(j).name);
            try
                N1 = load(filePath);
                if ~isfield(N1, 'analysisResults')
                    continue;
                end
                N = N1.analysisResults;

                loaded.analysisResults.centre = N.centre;
                loaded.analysisResults.distances = N.distances;
                loaded.analysisResults.bin_edges = N.bin_edges;
                loaded.analysisResults.bin_centers = N.bin_centers;
                loaded.analysisResults.response_percentage = N.response_percentage;
                loaded.analysisResults.bin_stats = N.bin_stats;
                loaded.analysisResults.neuronData = N.neuronData;

                Data = [Data; loaded.analysisResults]; 
            catch ME
                warning('unable to load file %s: %s', matFiles(j).name, ME.message);
            end
        end

        if isempty(Data)
            continue;
        end

        % stats：single file or multiple files（align bins）
        if isscalar(Data)
            mean_resp   = Data(1).response_percentage;
            smoothed_mean_resp = smooth(mean_resp);
            std_resp    = nan(size(mean_resp));
            bin_centers = Data(1).bin_centers;
        else
            all_bin_edges = arrayfun(@(x) x.bin_edges, Data, 'UniformOutput', false);
            common_bins = get_common_bins(all_bin_edges);
            [mean_resp, std_resp, bin_centers] = calculate_day_stats(Data, common_bins);
            smoothed_mean_resp = smooth(mean_resp);
        end

        res(i).nFiles = numel(Data);
        res(i).bin_centers = bin_centers;
        res(i).mean_resp = mean_resp;
        res(i).smoothed_mean_resp = smoothed_mean_resp;
        res(i).std_resp = std_resp;
    end

    % ---- output: stage table（bin/mean/std in cell） ----
    StageCol = repmat(string(stageName), nG, 1);
    GroupCol = string({res.Group})';
    nFilesCol = [res.nFiles]';

    binCell  = arrayfun(@(x) {x.bin_centers}, res)';
    meanCell = arrayfun(@(x) {x.mean_resp},   res)';
    smoothed_meanCell = arrayfun(@(x) {x.smoothed_mean_resp},   res)';
    stdCell  = arrayfun(@(x) {x.std_resp},    res)';

    stageTbl = table(StageCol, GroupCol, nFilesCol, binCell, meanCell, smoothed_meanCell, stdCell, ...
        'VariableNames', {'Stage','Group','nFiles','bin_centers','mean_resp','smoothed_mean_response','std_resp'});

    outRows = [outRows; stageTbl];

    % ---- plotting：this stage one plot，plot all non-spacer only when group's nFiles>0  ----
    figure('Position', [100 100 1100 650], 'Color','w'); hold on;

    idxPlot = find(~isSpacer & nFilesCol>0);

    % color：align with how many lines to plot
    colors = lines(max(1,numel(idxPlot)));

    leg = {};
    cc = 0;

    for k = 1:numel(idxPlot)
        gi = idxPlot(k);
        if isempty(res(gi).bin_centers) || isempty(res(gi).smoothed_mean_resp)
            continue;
        end

        cc = cc + 1;
        bc = res(gi).bin_centers(:);
        mr = res(gi).smoothed_mean_resp(:);
        sr = res(gi).std_resp(:);

        valid = ~isnan(mr) & ~isnan(bc);
        bc = bc(valid); mr = mr(valid);

        % mean curve
        plot(bc, mr, 'o-', 'LineWidth', 1.8, 'MarkerSize', 5, ...
            'Color', colors(cc,:), 'MarkerFaceColor', colors(cc,:));

        % std shading
        if ~all(isnan(sr))
            sr = sr(:);
            sr = sr(valid);
            lo = mr - sr;
            hi = mr + sr;
            fill([bc; flipud(bc)], [lo; flipud(hi)], colors(cc,:), ...
                'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility','off');
        end

        leg{end+1} = sprintf('%s (n=%d)', res(gi).Group, res(gi).nFiles); 
    end

    xlabel('Distance from center (\mum)');
    ylabel('Responsive Percentage (%)');
    title(['Responsiveness vs Distance - ' stageName], 'Interpreter','none');
    grid on; box on;

    if ~isempty(leg)
        legend(leg, 'Location','bestoutside', 'Interpreter','none');
    end

    % x/y limits（adapt to each curve）
    allX = [];
    allY = [];
    for k = 1:numel(idxPlot)
        gi = idxPlot(k);
        allX = [allX; res(gi).bin_centers(:)]; 
        allY = [allY; res(gi).smoothed_mean_resp(:)];  
    end
    allX = allX(~isnan(allX));
    allY = allY(~isnan(allY));
    if ~isempty(allX), xlim([min(allX) max(allX)]); end
    if ~isempty(allY), ylim([0 max(allY)*1.05]); end

    % save
    outPng = fullfile(baseFolder, ['Responsiveness_vs_Distance_' strrep(stageName,' ','_') '.png']);
    saveas(gcf, outPng);
    close(gcf);
end

outTbl = outRows;
writetable(outTbl, fullfile(baseFolder, 'summary_responsiveness_over_distance_early_late.xlsx'));

end


%% ========== helper: find common bins ==========
function common_bins = get_common_bins(all_bin_edges)
min_edge = min(cellfun(@min, all_bin_edges));
max_edge = max(cellfun(@max, all_bin_edges));
bin_widths = cellfun(@(x) x(2)-x(1), all_bin_edges);
finest_width = min(bin_widths);

common_bins = min_edge:finest_width:max_edge;
if common_bins(end) < max_edge
    common_bins = [common_bins, common_bins(end)+finest_width];
end
end

%% ========== helper: align + mean/std ==========
function [mean_resp, std_resp, bin_centers] = calculate_day_stats(Data, common_bins)
bin_centers = common_bins(1:end-1) + diff(common_bins)/2;
nBins  = length(bin_centers);
nFiles = length(Data);
all_resp = nan(nFiles, nBins);

for i = 1:nFiles
    [~, bin_idx] = histc(Data(i).bin_centers, common_bins);
    valid_bins = bin_idx > 0 & bin_idx <= nBins;
    all_resp(i, bin_idx(valid_bins)) = Data(i).response_percentage(valid_bins);
end

mean_resp = nanmean(all_resp, 1);
std_resp  = nanstd(all_resp, 0, 1);
end
