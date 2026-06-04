function plot_locations(clusterFilepath, neuronIDs, varargin)
% plot_locations
%
% Plots the spatial distribution of classified retinal ganglion cells (RGCs)
% and calculates distance-binned responsiveness from a user-selected centre.
%
% This function combines:
%   1. Cluster location information from an HS2/FLAME cluster file
%   2. Neuron response classification from a totalneurons file
%
% It then:
%   - Assigns each neuron a response label
%   - Allows the user to manually select a centre point
%   - Calculates the distance of each neuron from this centre
%   - Bins neurons by distance
%   - Calculates responsiveness across distance bins
%   - Saves the analysis into a *_distance.mat file
%   - Generates spatial and distance-response plots
%
% INPUTS:
%   clusterFilepath : path to the HS2/FLAME cluster file.
%
%                     If the cluster file exists, the function reads it
%                     using readHS2_FLAME().
%
%                     If the cluster file does not exist, the function tries
%                     to load a .mat file with the same basename.
%
%   neuronIDs       : path to the corresponding totalneurons file.
%
%                     This file should contain:
%                       totalneurons.num_OnNeurons.trans
%                       totalneurons.num_OnNeurons.sus
%                       totalneurons.num_OffNeurons.trans
%                       totalneurons.num_OffNeurons.sus
%                       totalneurons.num_OnOffNeurons
%                       totalneurons.num_notRespond
%                       totalneurons.num_unconventional, optional
%
% OPTIONAL PARAMETERS:
%   'BinWidth'      : width of distance bins in µm.
%                     Default = 50
%
%   'SmoothWindow'  : moving-average smoothing window for the responsiveness
%                     curve.
%                     Default = 3
%
% OUTPUTS:
%   This function does not return variables directly.
%
%   It saves:
%     1. *_distance.mat
%        containing analysisResults:
%           analysisResults.centre
%           analysisResults.distances
%           analysisResults.bin_edges
%           analysisResults.bin_centers
%           analysisResults.response_percentage
%           analysisResults.bin_stats
%           analysisResults.neuronData
%
%     2. *_responsiveness_vs_distance.png
%     3. *_responsiveness_vs_distance.tif
%     4. *_neuron_distribution.png
%     5. *_neuron_distribution.tif
%     6. *_neuron_distribution.pdf
%
% RESPONSE LABELS:
%   neuronData.labels:
%       0 = unlabelled / not assigned
%       1 = ON
%       2 = OFF
%       3 = ON-OFF
%       4 = unconventional
%       5 = non-responsive
%
% NOTE:
%   The selected centre is manually defined by the user and is usually the
%   optic nerve head or another anatomical reference point.

    %% Parse optional inputs
    p = inputParser;

    addParameter(p, 'BinWidth', 50, @isnumeric);
    addParameter(p, 'SmoothWindow', 3, @isnumeric);

    parse(p, varargin{:});

    %% Load cluster location data
    % If clusterFilepath does not exist, try loading a .mat file with the
    % same basename.
    %
    % Otherwise, read the HS2/FLAME cluster file using readHS2_FLAME().
    %
    % Expected fields in data:
    %   data.centres      : 2 × N matrix of neuron spatial coordinates
    %   data.channelNames : contains channel/cluster IDs

    if ~exist(clusterFilepath)
        matSavePath = extractBefore(clusterFilepath, '.');
        matSavePath = [matSavePath '.mat'];
        load(matSavePath);
    else
        data = readHS2_FLAME(clusterFilepath);
    end

    %% Define save-folder prefix from the cluster filepath
    % cellRasterFolder is used as the prefix for output files.
    % Example:
    %   clusterFilepath = 'recording_cluster.hdf5'
    %   cellRasterFolder = 'recording_cluster'

    cellRasterFolder = extractBefore(clusterFilepath, '.');

    %% Load neuron response classification file
    % neuronIDs should point to a totalneurons file.
    % The loaded structure is expected to contain totalneurons.

    s1 = load(neuronIDs);
    s = s1.totalneurons;

    %% Combine ON and OFF subtype IDs
    % ON cells include:
    %   - ON transient
    %   - ON sustained
    %
    % OFF cells include:
    %   - OFF transient
    %   - OFF sustained
    %
    % The try/catch is used because these arrays may be stored as either
    % column vectors or row vectors depending on the upstream processing.

    try
        ON_all = sort(cat(1, ...
            s.num_OnNeurons.trans, ...
            s.num_OnNeurons.sus));

        OFF_all = sort(cat(1, ...
            s.num_OffNeurons.trans, ...
            s.num_OffNeurons.sus));
    catch
        ON_all = sort(cat(2, ...
            s.num_OnNeurons.trans, ...
            s.num_OnNeurons.sus));

        OFF_all = sort(cat(2, ...
            s.num_OffNeurons.trans, ...
            s.num_OffNeurons.sus));
    end

    %% Extract unconventional neuron IDs if available
    % Some older totalneurons files may not contain num_unconventional.

    unconventional_all = [];

    if isfield(s, 'num_unconventional')
        unconventional_all = sort(s.num_unconventional(:));
    end

    %% Match classified neuron IDs to cluster/channel IDs
    % data.channelNames stores the channel/cluster IDs from the cluster file.
    % The response-classified IDs from totalneurons are matched against these
    % channel IDs using ismember().
    %
    % IDX_* contains the position of each classified neuron in channelIDs.
    % If no match is found, ismember returns 0.

    channelIDs = cell2mat(data.channelNames(4, :));

    [~, IDX_on]             = ismember(ON_all, channelIDs);
    [~, IDX_off]            = ismember(OFF_all, channelIDs);
    [~, IDX_on_off]         = ismember(s.num_OnOffNeurons, channelIDs);
    [~, IDX_not_respond]    = ismember(s.num_notRespond, channelIDs);
    [~, IDX_unconventional] = ismember(unconventional_all, channelIDs);

    %% Build neuronData structure
    % neuronData stores spatial coordinates, original channel IDs, and
    % response labels for all detected neurons/clusters.
    %
    % Label definitions:
    %   0 = unlabelled
    %   1 = ON
    %   2 = OFF
    %   3 = ON-OFF
    %   4 = unconventional
    %   5 = non-responsive

    neuronData = struct();

    neuronData.centres = data.centres;
    neuronData.channelIDs = channelIDs;
    neuronData.labels = zeros(size(channelIDs));

    neuronData.labels(IDX_on(IDX_on > 0)) = 1;
    neuronData.labels(IDX_off(IDX_off > 0)) = 2;
    neuronData.labels(IDX_on_off(IDX_on_off > 0)) = 3;
    neuronData.labels(IDX_unconventional(IDX_unconventional > 0)) = 4;
    neuronData.labels(IDX_not_respond(IDX_not_respond > 0)) = 5;

    %% Define distance-analysis output file
    % If the distance file already exists, the saved analysisResults are
    % loaded directly.
    %
    % If not, the user manually selects a centre point and the distance
    % analysis is performed.

    distanceFile = strcat(cellRasterFolder, "_distance.mat");

    if ~exist(distanceFile, 'file')

        %% Manually select the centre point
        % The user clicks one centre point on the neuron spatial map.
        % This is typically the optic nerve head or another reference centre.

        figure('Name', 'Select Center Point');

        scatter(data.centres(1, :), data.centres(2, :), ...
            30, [0.5 0.5 0.5], 'filled');

        title('Select a central point and press ENTER');

        [centre_x, centre_y] = ginput(1);

        close(gcf);

        %% Calculate distance from each neuron to the selected centre
        distances = sqrt((data.centres(1, :) - centre_x).^2 + ...
                         (data.centres(2, :) - centre_y).^2);

        %% Define distance bins
        % Bins start from 0 and extend to the maximum neuron distance.
        % The bin width is defined by the optional BinWidth parameter.

        bin_edges = 0:p.Results.BinWidth: ...
            ceil(max(distances) / p.Results.BinWidth) * p.Results.BinWidth;

        bin_centers = bin_edges(1:end-1) + p.Results.BinWidth / 2;

        %% Define responsive and included neurons
        % Responsive neurons include:
        %   1 = ON
        %   2 = OFF
        %   3 = ON-OFF
        %   4 = unconventional
        %
        % all_neurons includes responsive and non-responsive classified cells:
        %   1, 2, 3, 4, 5
        %
        % Label 0 neurons are excluded from the denominator.

        responsive_neurons = ismember(neuronData.labels, [1 2 3 4]);
        all_neurons = ismember(neuronData.labels, [1 2 3 4 5]);

        %% Calculate binned responsiveness
        [response_percentage, ~, bin_stats] = calculate_binned_response( ...
            distances, ...
            responsive_neurons, ...
            all_neurons, ...
            bin_edges);

        %% Save analysis results
        analysisResults = struct();

        analysisResults.centre = [centre_x, centre_y];
        analysisResults.distances = distances;
        analysisResults.bin_edges = bin_edges;
        analysisResults.bin_centers = bin_centers;
        analysisResults.response_percentage = response_percentage;
        analysisResults.bin_stats = bin_stats;
        analysisResults.neuronData = neuronData;

        save(distanceFile, 'analysisResults');

    else

        %% Load existing distance-analysis results
        load(distanceFile, 'analysisResults');

    end

    %% Plot responsiveness as a function of distance from centre
    plot_response_curve( ...
        analysisResults, ...
        p.Results.SmoothWindow, ...
        cellRasterFolder);

    %% Plot spatial distribution of classified neurons
    plot_neuron_distribution( ...
        analysisResults, ...
        cellRasterFolder);

end


function [response_percentage, valid_bins, bin_stats] = calculate_binned_response( ...
    distances, isresponsive, all, bin_edges)
% calculate_binned_response
%
% Calculates responsiveness across distance bins.
%
% INPUTS:
%   distances    : 1 × N vector of neuron distances from the selected centre
%
%   isresponsive : 1 × N logical vector
%                  true for responsive neurons:
%                       ON, OFF, ON-OFF, unconventional
%
%   all          : 1 × N logical vector
%                  true for all classified neurons included in the analysis:
%                       ON, OFF, ON-OFF, unconventional, non-responsive
%
%   bin_edges    : vector of distance bin edges
%
% OUTPUTS:
%   response_percentage : responsiveness value for each distance bin
%
%   valid_bins          : logical vector indicating whether each bin contains
%                         at least one included neuron
%
%   bin_stats           : structure containing:
%                           total_count
%                           responsive_count
%                           mean_distance
%
% IMPORTANT NOTE:
%   In the current implementation:
%
%       response_percentage(i) =
%           responsive_count_in_this_bin / total_number_of_all_classified_neurons * 100
%
%   This measures each bin's contribution to total responsiveness across the
%   whole retina.
%
%   If you want true within-bin responsiveness, use:
%
%       responsive_count_in_this_bin / total_count_in_this_bin * 100
%
%   instead.

    bin_stats = struct();

    response_percentage = zeros(1, length(bin_edges) - 1);
    valid_bins = false(1, length(bin_edges) - 1);

    %% Loop through distance bins
    for i = 1:length(bin_edges) - 1

        %% Identify neurons located within the current distance bin
        a = distances >= bin_edges(i) & distances < bin_edges(i + 1);

        %% Keep only classified neurons within this bin
        in_bin = all(a);

        %% Store bin-level counts
        bin_stats(i).total_count = sum(in_bin);
        bin_stats(i).responsive_count = sum(isresponsive(a));
        bin_stats(i).mean_distance = mean(distances(a));

        %% Calculate responsiveness for this bin
        if bin_stats(i).total_count > 0

            % Current version:
            % responsive neurons in this bin divided by all classified neurons
            % across the full retina.
            response_percentage(i) = ...
                bin_stats(i).responsive_count / sum(all) * 100;

            valid_bins(i) = true;

        else

            response_percentage(i) = NaN;

        end
    end
end


function plot_response_curve(results, smooth_window, saveFolder)
% plot_response_curve
%
% Plots binned responsiveness as a function of distance from the selected
% centre.
%
% INPUTS:
%   results       : analysisResults structure generated by plot_locations
%   smooth_window : moving-average smoothing window
%   saveFolder    : output file prefix
%
% OUTPUTS:
%   Saves:
%       *_responsiveness_vs_distance.png
%       *_responsiveness_vs_distance.tif

    figure('Position', [100 100 800 400]);

    %% Identify bins containing valid data
    valid_bins = ~isnan(results.response_percentage);

    %% Smooth and plot response curve
    % If smooth_window > 1, plot a moving-average-smoothed curve.
    % Currently, only the smoothed curve is plotted.

    if smooth_window > 1

        smoothed = movmean( ...
            results.response_percentage(valid_bins), ...
            smooth_window);

        hold on;

        plot(results.bin_centers(valid_bins), ...
             smoothed, ...
             'k-', ...
             'LineWidth', 2);
    end

    %% Format plot
    xlabel('Distance from center (μm)');
    ylabel('Responsiveness (%)');
    title('Cell Responsiveness vs. Distance from Center');
    grid on;

    %% Save figure
    % try/catch is used because saveFolder may sometimes be a char array
    % and sometimes a cell array depending on upstream handling.

    try
        saveas(gcf, [saveFolder, '_responsiveness_vs_distance.png']);
        saveas(gcf, [saveFolder, '_responsiveness_vs_distance.tif']);
    catch
        saveas(gcf, [saveFolder{:}, '_responsiveness_vs_distance.png']);
    end

    close(gcf);

end


function plot_neuron_distribution(results, saveFolder)
% plot_neuron_distribution
%
% Plots the spatial distribution of classified neurons.
%
% INPUTS:
%   results    : analysisResults structure generated by plot_locations
%   saveFolder : output file prefix
%
% OUTPUTS:
%   Saves:
%       *_neuron_distribution.png
%       *_neuron_distribution.tif
%       *_neuron_distribution.pdf
%
% COLOUR CODE:
%   red    = ON
%   blue   = OFF
%   green  = ON-OFF
%   orange = unconventional
%   grey   = non-responsive

    figure('Position', [100 100 600 600]);

    %% Define colours and labels
    colors = [
        1   0   0;      % ON, red
        0   0   1;      % OFF, blue
        0   1   0;      % ON-OFF, green
        1   0.6 0;      % unconventional, orange
        0.5 0.5 0.5     % non-responsive, grey
    ];

    labels = {'ON', 'OFF', 'ON-OFF', 'Unconventional', 'Non-responsive'};

    %% Plot each response class
    for i = 1:5

        idx = results.neuronData.labels == i;

        if i == 5

            % Non-responsive cells are plotted with partial transparency.
            scatter(results.neuronData.centres(1, idx), ...
                    results.neuronData.centres(2, idx), ...
                    50, colors(i, :), ...
                    'filled', ...
                    'MarkerFaceAlpha', 0.25);

        else

            % Responsive cells are plotted fully opaque.
            scatter(results.neuronData.centres(1, idx), ...
                    results.neuronData.centres(2, idx), ...
                    50, colors(i, :), ...
                    'filled', ...
                    'MarkerFaceAlpha', 1);
        end

        hold on;
    end

    %% Mark the selected centre
    plot(results.centre(1), results.centre(2), ...
         'kx', ...
         'MarkerSize', 15, ...
         'LineWidth', 2);

    %% Format plot
    title('RGC Type Distribution');
    legend([labels, {'Center'}], 'Location', 'best');

    grid on;
    axis equal;

    ax = gca;
    box(ax, 'on');

    xlim([0 4000]);
    ylim([0 4000]);

    % Hide tick labels to emphasise spatial pattern rather than exact axes.
    set(gca, 'XTickLabel', [], 'YTickLabel', []);

    %% Save figure
    saveas(gcf, [saveFolder, '_neuron_distribution.png']);
    saveas(gcf, [saveFolder, '_neuron_distribution.tif']);

    try
        saveas(gcf, [saveFolder, '_neuron_distribution.png']);

        exportgraphics(gcf, ...
            [saveFolder, '_neuron_distribution.pdf'], ...
            'ContentType', 'vector');

    catch
        saveas(gcf, [saveFolder{:}, '_neuron_distribution.png']);

        exportgraphics(gcf, ...
            [saveFolder{:}, '_neuron_distribution.pdf'], ...
            'ContentType', 'vector');
    end

    close(gcf);

end