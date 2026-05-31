function plot_locations(clusterFilepath, neuronIDs, varargin)
% Plots the locations of neurons in relation to the center of the retina on
% MEA. Also plots the cell type clusters on the MEA grid so we can align
% things. Requires user input to choose the center of the retina
%
% Inputs: clusterFilepath- fullfile path to *_cluster.hdf5
%
%         neuronsIDs - fullfile to '*_totalneuronsV2.mat'which contains the
%                      cluster ID classifications

%% defaults
p = inputParser;
addParameter(p, 'BinWidth', 50, @isnumeric); %
addParameter(p, 'SmoothWindow', 3, @isnumeric); %
parse(p, varargin{:});

%%  load data
if ~exist(clusterFilepath)
    matSavePath = extractBefore(clusterFilepath, '.');
    matSavePath = [matSavePath '.mat'];
    load(matSavePath);
else
    data = readHS2_FLAME(clusterFilepath);
end

cellRasterFolder = extractBefore(clusterFilepath, '.');
s1 = load(neuronIDs);
s = s1.totalneurons;

%% filter cell IDs into ON vs OFF etc
try
    ON_all = sort(cat(1, s.num_OnNeurons.trans, s.num_OnNeurons.sus));
    OFF_all = sort(cat(1, s.num_OffNeurons.trans, s.num_OffNeurons.sus));
catch
    ON_all = sort(cat(2, s.num_OnNeurons.trans, s.num_OnNeurons.sus));
    OFF_all = sort(cat(2, s.num_OffNeurons.trans, s.num_OffNeurons.sus));
end

unconventional_all = [];
if isfield(s, 'num_unconventional')
    unconventional_all = sort(s.num_unconventional(:));
end

% get all the channel indexs for the various types
channelIDs = cell2mat(data.channelNames(4,:));
[~, IDX_on] = ismember(ON_all, channelIDs);
[~, IDX_off] = ismember(OFF_all, channelIDs);
[~, IDX_on_off] = ismember(s.num_OnOffNeurons, channelIDs);
[~, IDX_not_respond] = ismember(s.num_notRespond, channelIDs);
[~, IDX_unconventional] = ismember(unconventional_all, channelIDs);

% put data into structure
neuronData = struct();
neuronData.centres = data.centres;
neuronData.channelIDs = channelIDs;
neuronData.labels = zeros(size(channelIDs)); % 1=ON, 2=OFF, etc.
neuronData.labels(IDX_on(IDX_on>0)) = 1;
neuronData.labels(IDX_off(IDX_off>0)) = 2;
neuronData.labels(IDX_on_off(IDX_on_off>0)) = 3;
neuronData.labels(IDX_unconventional(IDX_unconventional>0)) = 4;
neuronData.labels(IDX_not_respond(IDX_not_respond>0)) = 5;

distanceFile = strcat(cellRasterFolder, "_distance.mat");
% if distance file is not present, do the center selection
if ~exist(distanceFile, 'file')

    figure('Name', 'Select Center Point');
    scatter(data.centres(1,:), data.centres(2,:), 30, [0.5 0.5 0.5], 'filled');
    title('Select a central point and press ENTER');
    [centre_x, centre_y] = ginput(1);
    close(gcf);

    % get distances from the centre picked location
    distances = sqrt((data.centres(1,:) - centre_x).^2 + (data.centres(2,:) - centre_y).^2);

    % get the bin edges/centres for plotting
    bin_edges = 0:p.Results.BinWidth:ceil(max(distances)/p.Results.BinWidth)*p.Results.BinWidth;
    bin_centers = bin_edges(1:end-1) + p.Results.BinWidth/2;

    % calcuate responsivity 
    responsive_neurons = ismember(neuronData.labels, [1 2 3 4]); % ON, OFF, ON-OFF, unconventional
    [response_percentage, ~, bin_stats] = calculate_binned_response(...
        distances, responsive_neurons, bin_edges);

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
    load(distanceFile, 'analysisResults');
end

% plot distance from center for cluster responses
plot_response_curve(analysisResults, p.Results.SmoothWindow, cellRasterFolder);

% plot cell clusters on MEA chip classified by cell type
plot_neuron_distribution(analysisResults, cellRasterFolder);
end


