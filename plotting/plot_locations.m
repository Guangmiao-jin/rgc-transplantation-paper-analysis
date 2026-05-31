function plot_locations(clusterFilepath, neuronIDs, varargin)
    
    p = inputParser;
    addParameter(p, 'BinWidth', 50, @isnumeric); % 
    addParameter(p, 'SmoothWindow', 3, @isnumeric); % 
    parse(p, varargin{:});
    
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
    
 
    channelIDs = cell2mat(data.channelNames(4,:));
    [~, IDX_on] = ismember(ON_all, channelIDs);
    [~, IDX_off] = ismember(OFF_all, channelIDs);
    [~, IDX_on_off] = ismember(s.num_OnOffNeurons, channelIDs);
    [~, IDX_not_respond] = ismember(s.num_notRespond, channelIDs);
    [~, IDX_unconventional] = ismember(unconventional_all, channelIDs);
    
  
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
    if ~exist(distanceFile, 'file')
       
        figure('Name', 'Select Center Point');
        scatter(data.centres(1,:), data.centres(2,:), 30, [0.5 0.5 0.5], 'filled');
        title('Select a central point and press ENTER');
        [centre_x, centre_y] = ginput(1);
        close(gcf);
        
        
        distances = sqrt((data.centres(1,:) - centre_x).^2 + (data.centres(2,:) - centre_y).^2);
        
     
        bin_edges = 0:p.Results.BinWidth:ceil(max(distances)/p.Results.BinWidth)*p.Results.BinWidth;
        bin_centers = bin_edges(1:end-1) + p.Results.BinWidth/2;
        
        responsive_neurons = ismember(neuronData.labels, [1 2 3 4]); % ON, OFF, ON-OFF, unconventional
        all_neurons = ismember(neuronData.labels, [1 2 3 4 5]);
        [response_percentage, ~, bin_stats] = calculate_binned_response(...
            distances, responsive_neurons, all_neurons, bin_edges);
        
     
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
    
  
    plot_response_curve(analysisResults, p.Results.SmoothWindow, cellRasterFolder);
    
 
    plot_neuron_distribution(analysisResults, cellRasterFolder);
end

%% 
function [response_percentage, valid_bins, bin_stats] = calculate_binned_response(distances, isresponsive, all, bin_edges)
    bin_stats = struct();
    response_percentage = zeros(1, length(bin_edges)-1);
    valid_bins = false(1, length(bin_edges)-1);
    
    for i = 1:length(bin_edges)-1
        a = (distances >= bin_edges(i)) & (distances < bin_edges(i+1));
        in_bin = all(a);
        bin_stats(i).total_count = sum(in_bin);
        bin_stats(i).responsive_count = sum(isresponsive(a));
        bin_stats(i).mean_distance = mean(distances(a));
        
        if bin_stats(i).total_count > 0
            response_percentage(i) = bin_stats(i).responsive_count / sum(all) * 100;
            valid_bins(i) = true;
        else
            response_percentage(i) = NaN;
        end
    end
end

%% 
function plot_response_curve(results, smooth_window, saveFolder)
    figure('Position', [100 100 800 400]);
    
    valid_bins = ~isnan(results.response_percentage);
    
    % 
    if smooth_window > 1
        smoothed = movmean(results.response_percentage(valid_bins), smooth_window);
        hold on;
        plot(results.bin_centers(valid_bins), smoothed, 'k-', 'LineWidth', 2);
    end
    %ylim([0 12]);
    xlabel('Distance from center (μm)');
    ylabel('Responsiveness (%)');
    title('Cell Responsiveness vs. Distance from Center');
    grid on;

    try
        saveas(gcf, [saveFolder, '_responsiveness_vs_distance.png']);
        %saveas(gcf, [saveFolder, '_responsiveness_vs_distance.tif']);
    catch
        saveas(gcf, [saveFolder{:}, '_responsiveness_vs_distance.png']);
    end
    close(gcf);
end

%% 
function plot_neuron_distribution(results, saveFolder)
    %figure('Position', [100 100 600 600]);
    
    % 
    colors = [1 0 0; 0 0 1; 0 1 0; 1 0.6 0; 0.5 0.5 0.5];
    labels = {'ON', 'OFF', 'ON-OFF', 'Unconventional', 'Non-responsive'};
    
    % 
    for i = 1:5
        if i == 5
            idx = (results.neuronData.labels == i);
            scatter(results.neuronData.centres(1,idx), results.neuronData.centres(2,idx), ...
            50, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.25);
            hold on;
        else
            idx = (results.neuronData.labels == i);
            scatter(results.neuronData.centres(1,idx), results.neuronData.centres(2,idx), ...
            50, colors(i,:), 'filled', 'MarkerFaceAlpha', 1);
            hold on;
        end
    end
    
    % 
    plot(results.centre(1), results.centre(2), 'kx', 'MarkerSize', 15, 'LineWidth', 2);
   
    %title('RGC Type Distribution');
    %legend([labels, {'Center'}], 'Location', 'best');
    grid on;
    axis equal;
    ax = gca; % Get current axes
    box(ax, 'on');
    xlim([0 4000]);
    ylim([0 4000]);
    set(gca, 'XTickLabel',[], 'YTickLabel',[]);
    %saveas(gcf, [saveFolder, '_neuron_distribution.png']);
    %saveas(gcf, [saveFolder, '_neuron_distribution.tif']);
    try
        %saveas(gcf, [saveFolder, '_neuron_distribution.png']);
        exportgraphics(gcf, [saveFolder, '_neuron_distribution.pdf'],  'ContentType','vector');
    catch
        %saveas(gcf, [saveFolder{:}, '_neuron_distribution.png']);
        exportgraphics(gcf, [saveFolder{:}, '_neuron_distribution.pdf'],  'ContentType','vector');
    end
    close(gcf);
end
