function outTbl = plot_combined_responsiveness_transplanted(baseFolder)
% plot_combined_responsiveness_contribution
%
% Plots the spatial distribution of responsive RGCs as a function of
% distance from the optic nerve head (ONH), expressed as each distance
% bin's CONTRIBUTION to the retina's total responsive percentage.
%
% For each recording:
%      contribution(bin) = N_responsive_in_bin
%                          / N_total_neurons_in_whole_retina
%                          * 100
%
% By construction:
%      sum(contribution across bins) = overall responsive percentage
%
% This makes the area under each curve directly interpretable as the
% retina's total responsive percentage.
%
% LABEL MAPPING (analysisResults.neuronData.labels):
%      0 = noise           (excluded from total)
%      1 = ON              (responsive)
%      2 = OFF             (responsive)
%      3 = ON-OFF          (responsive)
%      4 = unconventional  (responsive)
%      5 = unresponsive    (counted in total, not in responsive)
%
% NOTE on the denominator:
%   The total neuron count excludes label 0 (noise) but includes label 5
%   (unresponsive). This matches the conventional definition of
%   "responsive % out of all detected, non-noise units".
%
% OUTPUT:
%   outTbl : table with one row per (Stage, Group), containing mean
%            contribution profile, SD, and number of retinas in that
%            group.

%% Define experimental stages and behavioural categories
stages    = {'early degeneration', 'late degeneration'};
behavPos1 = 'behavioural positive and ephys positive';
behavPos2 = 'behavioural positive and ephys negative';
behavNeg  = 'behavioural negative';

%% Define plotting order
% Spacers create visual gaps between treatment groups in the legend.
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

%% Define folder lookup for each plot group
% SHAM is searched directly; treated groups are searched under their
% behavioural-classification subfolder.
scanSpec = { ...
    {behavPos1, 'NRL'}, ...
    {behavPos2, 'NRL'}, ...
    {behavNeg , 'NRL'}, ...
    {}, ...
    {behavPos1, 'CRX'}, ...
    {behavPos2, 'CRX'}, ...
    {behavNeg , 'CRX'}, ...
    {}, ...
    {'SHAM'} ...
};

% Labels considered responsive.
respLabels = [1 2 3 4];

outRows = [];

%% Loop through degeneration stages
for s = 1:numel(stages)
    stageName = stages{s};
    stagePath = fullfile(baseFolder, stageName);

    nG       = numel(plotGroupNames);
    groups   = string(plotGroupNames(:));
    isSpacer = startsWith(groups, "Spacer");

    %% Initialise per-group results container
    res = struct();
    for i = 1:nG
        res(i).Stage = stageName;
        res(i).Group = plotGroupNames{i};
        res(i).nFiles = 0;
        res(i).bin_centers = [];
        res(i).mean_contrib = [];
        res(i).smoothed_mean_contrib = [];
        res(i).std_contrib = [];
        res(i).mean_total_pct = NaN;   % verification: should match sum(mean_contrib)
    end

    %% Loop through plotting groups within this stage
    for i = 1:nG
        if isempty(scanSpec{i})
            continue;   % spacer
        end
        spec = scanSpec{i};

        %% Resolve folder path
        % SHAM:  <stage>/SHAM
        % NRL/CRX: <stage>/<treatment>/<behavioural folder>
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

        %% Find all distance files recursively under this folder
        matFiles = dir(fullfile(groupPath, '**', '*_distance.mat'));
        if isempty(matFiles)
            continue;
        end

        %% Compute per-file contribution profiles
        Data = [];

        for j = 1:numel(matFiles)
            filePath = fullfile(matFiles(j).folder, matFiles(j).name);
            try
                N1 = load(filePath);
                if ~isfield(N1, 'analysisResults')
                    continue;
                end
                N = N1.analysisResults;

                %% Extract required fields
                edges     = N.bin_edges(:)';
                distances = N.distances(:);
                labels    = N.neuronData.labels(:);

                %% Denominator: all non-noise neurons in this retina
                isNonNoise   = labels ~= 0;
                nTotalRetina = sum(isNonNoise);
                if nTotalRetina == 0
                    continue;
                end

                %% Numerator: responsive neurons binned by distance
                isResp = ismember(labels, respLabels);
                count_per_bin = histcounts(distances(isResp), edges);

                %% Contribution per bin (%)
                % Each bin's contribution sums to the retina's overall
                % responsive percentage.
                contrib_per_bin = count_per_bin / nTotalRetina * 100;

                %% Store per-file profile
                entry.bin_edges    = edges;
                entry.bin_centers  = edges(1:end-1) + diff(edges)/2;
                entry.contrib      = contrib_per_bin;
                entry.total_pct    = 100 * sum(isResp) / nTotalRetina;
                Data = [Data; entry]; %#ok<AGROW>

            catch ME
                warning('Unable to process %s: %s', matFiles(j).name, ME.message);
            end
        end

        if isempty(Data)
            continue;
        end

        %% Aggregate across retinas in this group
        % Single retina: profile is taken as-is; SD is NaN.
        % Multiple retinas: align to a common bin grid, then average
        % each retina equally.
        if isscalar(Data)
            mean_contrib = Data(1).contrib;
            std_contrib  = nan(size(mean_contrib));
            bin_centers  = Data(1).bin_centers;
            mean_total_pct = Data(1).total_pct;
        else
            all_edges    = arrayfun(@(x) x.bin_edges, Data, 'UniformOutput', false);
            common_edges = get_common_bins(all_edges);
            [mean_contrib, std_contrib, bin_centers] = ...
                aggregate_contrib(Data, common_edges);
            mean_total_pct = mean([Data.total_pct], 'omitnan');
        end

        smoothed_mean_contrib = movmean(mean_contrib, 5, 'omitnan');

        %% Save into results struct
        res(i).nFiles                = numel(Data);
        res(i).bin_centers           = bin_centers;
        res(i).mean_contrib          = mean_contrib;
        res(i).smoothed_mean_contrib = smoothed_mean_contrib;
        res(i).std_contrib           = std_contrib;
        res(i).mean_total_pct        = mean_total_pct;
    end

    %% Build per-stage table
    StageCol     = repmat(string(stageName), nG, 1);
    GroupCol     = string({res.Group})';
    nFilesCol    = [res.nFiles]';
    totalPctCol  = [res.mean_total_pct]';

    binCell    = arrayfun(@(x) {x.bin_centers},           res)';
    contribCell= arrayfun(@(x) {x.mean_contrib},          res)';
    smoothCell = arrayfun(@(x) {x.smoothed_mean_contrib}, res)';
    sdCell     = arrayfun(@(x) {x.std_contrib},           res)';

    stageTbl = table(StageCol, GroupCol, nFilesCol, totalPctCol, ...
        binCell, contribCell, smoothCell, sdCell, ...
        'VariableNames', {'Stage','Group','nFiles','mean_total_responsive_pct', ...
            'bin_centers','mean_contribution_pct', ...
            'smoothed_mean_contribution_pct','std_contribution_pct'});

    outRows = [outRows; stageTbl]; %#ok<AGROW>

    %% Plot all non-empty groups for this stage
    figure('Position', [100 100 1100 650], 'Color','w'); hold on;
    idxPlot = find(~isSpacer & nFilesCol > 0);
    colors  = lines(max(1, numel(idxPlot)));

    leg = {};
    cc  = 0;
    for k = 1:numel(idxPlot)
        gi = idxPlot(k);
        if isempty(res(gi).bin_centers) || isempty(res(gi).smoothed_mean_contrib)
            continue;
        end

        cc = cc + 1;
        bc = res(gi).bin_centers(:);
        mc = res(gi).smoothed_mean_contrib(:);
        sd = res(gi).std_contrib(:);

        valid = ~isnan(mc) & ~isnan(bc);
        bc = bc(valid); mc = mc(valid);

        %% Mean curve
        plot(bc, mc, 'o-', 'LineWidth', 1.8, 'MarkerSize', 5, ...
            'Color', colors(cc,:), 'MarkerFaceColor', colors(cc,:));

        %% SD shading (only when SD is available, i.e. >1 retina)
        if ~all(isnan(sd))
            sd = sd(valid);
            lo = mc - sd;
            hi = mc + sd;
            fill([bc; flipud(bc)], [lo; flipud(hi)], colors(cc,:), ...
                'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility','off');
        end

        leg{end+1} = sprintf('%s (n=%d, total=%.1f%%)', ...
            res(gi).Group, res(gi).nFiles, res(gi).mean_total_pct); %#ok<AGROW>
    end

    xlabel('Distance from ONH (\mum)');
    ylabel('Contribution to total responsive % (per bin)');
    title(['Spatial distribution of responsiveness - ' stageName], ...
        'Interpreter','none');
    grid on; box on;

    if ~isempty(leg)
        legend(leg, 'Location','bestoutside', 'Interpreter','none');
    end

    %% Adapt axis limits to the actual data range
    allX = []; allY = [];
    for k = 1:numel(idxPlot)
        gi = idxPlot(k);
        allX = [allX; res(gi).bin_centers(:)];          %#ok<AGROW>
        allY = [allY; res(gi).smoothed_mean_contrib(:)];%#ok<AGROW>
    end
    allX = allX(~isnan(allX));
    allY = allY(~isnan(allY));
    if ~isempty(allX), xlim([min(allX) max(allX)]); end
    if ~isempty(allY), ylim([0 max(allY)*1.05]); end

    %% Save figure
    outPng = fullfile(baseFolder, ...
        ['ResponsivenessContribution_vs_Distance_' strrep(stageName,' ','_') '.png']);
    saveas(gcf, outPng);
    close(gcf);
end

%% Write summary spreadsheet
outTbl = outRows;
writetable(outTbl, fullfile(baseFolder, ...
    'summary_responsiveness_contribution_early_late.xlsx'));
end


%% ========== helper: find common bin edges across files ==========
function common_edges = get_common_bins(all_bin_edges)
% Returns a common edge vector spanning the union of all files' bins,
% using the finest bin width found across files.
min_edge     = min(cellfun(@min, all_bin_edges));
max_edge     = max(cellfun(@max, all_bin_edges));
bin_widths   = cellfun(@(x) x(2)-x(1), all_bin_edges);
finest_width = min(bin_widths);

common_edges = min_edge:finest_width:max_edge;
if common_edges(end) < max_edge
    common_edges = [common_edges, common_edges(end)+finest_width];
end
end


%% ========== helper: align + mean/std for contribution profiles ==========
function [mean_contrib, std_contrib, bin_centers] = ...
    aggregate_contrib(Data, common_edges)
% Aligns each file's contribution profile onto common bin edges and
% computes file-level mean and SD. Each retina contributes equally.
bin_centers = common_edges(1:end-1) + diff(common_edges)/2;
nBins  = length(bin_centers);
nFiles = length(Data);

all_contrib = nan(nFiles, nBins);

for i = 1:nFiles
    bin_idx    = discretize(Data(i).bin_centers, common_edges);
    valid_bins = ~isnan(bin_idx);
    all_contrib(i, bin_idx(valid_bins)) = Data(i).contrib(valid_bins);
end

mean_contrib = mean(all_contrib, 1, 'omitnan');
std_contrib  = std(all_contrib, 0, 1, 'omitnan');
end