function ax = plotAllRasterPSTHs_On_Off_Response_example(plotMetrics)
% plotAllRasterPSTHs_On_Off_Response_example
%
% Plots raster plots and PSTHs for one example neuron across light
% stimulation conditions.
%
% For each condition, the function generates:
%   1. A raster plot showing spike timing across trials
%   2. A PSTH showing the trial-averaged firing rate over time
%
% The expected layout is:
%   - If there is one condition:
%         2 rows × 1 column
%         row 1 = raster
%         row 2 = PSTH
%
%   - If there are multiple conditions:
%         2 rows × 3 columns
%         top row    = raster plots
%         bottom row = PSTHs
%
% INPUT:
%   plotMetrics : structure containing spike and PSTH data for one neuron.
%
%                 Required fields:
%
%                 plotMetrics.trialSpikes
%                     1 × nConditions cell array.
%                     Each cell contains trial-wise spike times.
%
%                     Example:
%                         plotMetrics.trialSpikes{condition}{trial}
%                         = spike times for one trial
%
%                 plotMetrics.trialPSTHs
%                     1 × nConditions cell array.
%                     Each cell contains a trials × time-bins PSTH matrix.
%
%                 plotMetrics.binEdges
%                     1 × nBins+1 vector of PSTH bin edges in seconds.
%
% OUTPUT:
%   ax : array of axes handles.
%
%        For three conditions:
%           ax(1:3) = raster axes
%           ax(4:6) = PSTH axes
%
%        For one condition:
%           ax(1) = raster axis
%           ax(2) = PSTH axis
%
% STIMULUS TIMING:
%   The time window is assumed to be:
%       -2 to 0 s : light ON period
%        0 to 2 s : light OFF period
%
%   A yellow bar marks the light ON period.
%   A grey bar marks the light OFF period.
%
% NOTES:
%   - x = 0 is marked with a red vertical line, corresponding to light
%     offset.
%   - For ON responses, light onset occurs at x = -2.
%   - For OFF responses, light offset occurs at x = 0.
%   - The helper function subplotEvenAxes() is assumed to be available
%     elsewhere in the MATLAB path.

%% Define condition titles
% These correspond to the three light conditions used in the experiment.

titleText = {'Scotopic', 'Mesopic', 'Photopic'};

%% Create full-screen figure
figH = figure( ...
    'units', 'normalized', ...
    'outerposition', [0 0 1 1], ...
    'Color', 'white', ...
    'MenuBar', 'none');

%% Determine the number of stimulus blocks / light conditions
nBlocks = length(plotMetrics.trialSpikes);

%% Loop through each light condition
for stimblk = 1:nBlocks

    %% Create raster subplot
    % For a single condition, use a 2 × 1 layout.
    % For multiple conditions, use a 2 × 3 layout:
    %   top row    = raster plots
    %   bottom row = PSTHs.

    if nBlocks == 1
        ax(stimblk) = subplot(2, 1, stimblk);
        hold on;
    else
        ax(stimblk) = subplot(2, 3, stimblk);
        hold on;
    end

    %% Extract trial-wise spike times for this condition
    trialSpikesCnd = plotMetrics.trialSpikes{stimblk};

    % xSpikePos and ySpikePos store line coordinates for all raster ticks.
    xSpikePos = [];
    ySpikePos = [];

    %% Build raster line coordinates trial by trial
    for tr = 1:length(trialSpikesCnd)

        trSpikes = trialSpikesCnd{tr};

        % Skip empty trials with no spikes.
        if isempty(trSpikes)
            continue;
        end

        % Each spike is drawn as a short vertical line.
        % xSpikePosTemp has two rows so each spike has a start and end point.
        xSpikePosTemp = repmat(trSpikes', 2, 1);

        % ySpikePosTemp defines the vertical extent of each raster tick.
        % The tick spans from trial-1 to trial.
        ySpikePosTemp(1, :) = tr - 1;
        ySpikePosTemp(2, :) = tr;
        ySpikePosTemp = repmat(ySpikePosTemp, 1, size(xSpikePosTemp, 2));

        % Append this trial's spike coordinates to the full raster arrays.
        xSpikePos = [xSpikePos xSpikePosTemp];
        ySpikePos = [ySpikePos ySpikePosTemp];

        % Clear temporary variables before the next trial.
        trSpikes = [];
        ySpikePosTemp = [];
        xSpikePosTemp = [];
    end

    %% Plot raster ticks
    plot(xSpikePos, ySpikePos, 'Color', 'k');

    %% Add light ON and light OFF bars above the raster
    % The stimulus is assumed to be:
    %   -2 to 0 s = light ON
    %    0 to 2 s = light OFF

    lightWindow = [-2 0];
    darkWindow  = [0 2];

    nTrials = length(trialSpikesCnd);

    % Bar height and vertical position above the raster trials.
    barH = 2;
    gap  = 1;
    y0   = nTrials + gap;

    currAx = ax(stimblk);

    % Extend y-axis to include stimulus bars above the raster.
    currAx.YLim = [0 nTrials + gap + barH + 0.2];

    %% Light ON bar
    patch(currAx, ...
        [lightWindow(1) lightWindow(2) lightWindow(2) lightWindow(1)], ...
        [y0 y0 y0 + barH y0 + barH], ...
        [1 0.9 0.1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1, ...
        'HandleVisibility', 'on');

    text(mean(lightWindow), y0 + barH / 2, 'light on', ...
        'Parent', currAx, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'Color', 'k');

    %% Light OFF bar
    patch(currAx, ...
        [darkWindow(1) darkWindow(2) darkWindow(2) darkWindow(1)], ...
        [y0 y0 y0 + barH y0 + barH], ...
        [0.5 0.5 0.5], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1, ...
        'HandleVisibility', 'on');

    text(mean(darkWindow), y0 + barH / 2, 'light off', ...
        'Parent', currAx, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'Color', 'k');

    %% Format raster axis
    set(currAx, 'Box', 'on', 'LineWidth', 1.2);

    currAx.XLim = [-2 plotMetrics.binEdges(end)];
    currAx.XLabel.String = 'Time (s)';
    currAx.YLabel.String = 'Trials';

    % Mark light offset at 0 s.
    xline(0, 'Color', 'r', 'LineWidth', 2);

    % Add condition title.
    if stimblk <= numel(titleText)
        title(titleText{stimblk});
    else
        title(sprintf('Condition %d', stimblk));
    end

    %% Create PSTH subplot
    if nBlocks == 1
        ax(stimblk + 1) = subplot(2, 1, stimblk + 1);
    else
        ax(stimblk + 3) = subplot(2, 3, stimblk + 3);
    end

    %% Calculate trial-averaged firing rate
    % plotMetrics.trialPSTHs{stimblk} is expected to be:
    %   trials × time bins
    %
    % The PSTH is converted to spikes/s by dividing by the bin width.

    binWidth = mean(diff(plotMetrics.binEdges));  % seconds

    meanPSTH = mean(plotMetrics.trialPSTHs{stimblk}, 1, 'omitnan');
    countAverageSec = meanPSTH ./ binWidth;

    %% Plot PSTH as a histogram
    h = histogram( ...
        'BinCounts', countAverageSec, ...
        'BinEdges', plotMetrics.binEdges);

    h.FaceColor = 'k';

    hold on;

    % Mark light offset at 0 s.
    xline(0, 'Color', 'r', 'LineWidth', 2);

    %% Format PSTH axis
    if nBlocks == 1
        currAx = ax(stimblk + 1);
    else
        currAx = ax(stimblk + 3);
    end

    currAx.XLim = [-2 plotMetrics.binEdges(end)];

    % Set y-axis limit slightly above the maximum PSTH value.
    mVal = max(h.Values) + round(max(h.Values) * 0.1);

    % Fix empty or silent PSTH.
    if isempty(mVal) || isnan(mVal) || mVal == 0
        mVal = 1;
    end

    currAx.YLim = [0 mVal];
    currAx.XLabel.String = 'Time (s)';
    currAx.YLabel.String = 'Average spikes per second';

end

%% Match subplot axes if multiple conditions are plotted
% subplotEvenAxes is a custom helper function and must be available on the
% MATLAB path. The argument [4 5 6] likely refers to the PSTH axes.
if nBlocks ~= 1
    subplotEvenAxes(ax, [0 1 0], [4 5 6]);
end

end