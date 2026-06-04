function chipGridAlign(imgPath, matPath, varargin)
% chipGridAlign
%
% Manually aligns a 64 × 64 MEA chip grid to a retinal TIFF image, overlays
% responsive clusters on the aligned image, and recalculates each neuron's
% distance from a user-defined reference centre.
%
% INPUTS:
%   imgPath : path to the retinal TIFF image
%             e.g. '/path/to/retina_image.tif'
%
%   matPath : path to the corresponding *_distance.mat file containing
%             analysisResults.neuronData
%
% OPTIONAL PARAMETERS:
%   'BinWidth'     : bin width for distance-binned responsiveness analysis
%                    default = 50 µm
%
%   'SmoothWindow' : smoothing window size
%                    default = 3
%                    Note: currently parsed but not used in this function.
%
% OUTPUTS:
%   1. Saves an aligned TIFF image with responsive clusters labelled:
%        *_aligned.tif
%
%   2. Updates matPath with a new analysisResults structure containing:
%        analysisResults.centre
%        analysisResults.distances
%        analysisResults.bin_edges
%        analysisResults.bin_centers
%        analysisResults.response_percentage
%        analysisResults.bin_stats
%        analysisResults.neuronData
%        analysisResults.neuronVirtualCoords
%
%   3. Exports channelCoords to the MATLAB base workspace:
%        channelCoords = [chip_col_index, chip_row_index, image_x, image_y]
%
% WORKFLOW:
%   Step 1: User manually selects four corner channels on the TIFF image:
%           top-left, top-right, bottom-left, bottom-right.
%
%   Step 2: A 64 × 64 virtual grid is generated from these corner points.
%
%   Step 3: Responsive clusters are overlaid on the TIFF image.
%
%   Step 4: User manually selects a new reference centre, usually the optic
%           nerve head or another biologically meaningful centre.
%
%   Step 5: Distances from all neurons to the selected centre are recalculated,
%           and responsiveness is binned by distance.

%% Parse optional inputs
p = inputParser;
addParameter(p, 'BinWidth', 50, @isnumeric);
addParameter(p, 'SmoothWindow', 3, @isnumeric);
parse(p, varargin{:});

%% Read and display the TIFF image
I = imread(imgPath);

figure('Name', 'Chip Grid Calibration', 'NumberTitle', 'off');
hAx = axes;
imshow(I, [], 'Parent', hAx);

title(['Step 1: Zoom & pick the FOUR corner channels, ' ...
       'top left -> top right -> bottom left -> bottom right']);

%% Step 1: Manually select the four corner channel positions
% The user clicks four points on the TIFF image.
% Required order:
%   1 = top-left corner channel
%   2 = top-right corner channel
%   3 = bottom-left corner channel
%   4 = bottom-right corner channel
%
% Only the first three points are used to define the grid geometry below.
% The fourth point is useful for visual/manual confirmation.

cornerPt = [];   % 4 × 2 matrix storing [x, y] pixel coordinates

for k = 1:4
    h = drawpoint('Color', 'r', ...
                  'MarkerSize', 10, ...
                  'Label', sprintf('%d', k));
    cornerPt(k, :) = h.Position;
end

%% Step 2: Calculate the 64 × 64 grid coordinates in image pixel space
% Define the three main corner coordinates.
% c00 = top-left
% c0N = top-right
% cM0 = bottom-left

c00 = cornerPt(1, :);
c0N = cornerPt(2, :);
cM0 = cornerPt(3, :);

% Calculate the pixel displacement per channel step.
% Since the grid contains 64 channels per axis, there are 63 intervals.
%
% ux: one-channel step vector along the horizontal direction
% uy: one-channel step vector along the vertical direction

ux = (c0N - c00) / (64 - 1);
uy = (cM0 - c00) / (64 - 1);

% Generate chip indices from 0 to 63 in both x and y directions.
% Xidx and Yidx are 64 × 64 matrices describing the virtual grid position.

[Xidx, Yidx] = meshgrid(0:63, 0:63);

% Convert chip indices into image pixel coordinates.
% gridX and gridY store the pixel location of each virtual channel.

gridX = c00(1) + ux(1) * Xidx + uy(1) * Yidx;
gridY = c00(2) + ux(2) * Xidx + uy(2) * Yidx;

%% Step 3: Visualise the aligned grid on top of the TIFF image
hold on;

% Plot all 64 × 64 virtual channel positions.
scatter(gridX(:), gridY(:), 15, 'y', 'filled');

% Draw the four borders of the aligned grid.
plot([gridX(1, 1)   gridX(1, end)], ...
     [gridY(1, 1)   gridY(1, end)], ...
     'y-', 'LineWidth', 1.5);

plot([gridX(end, 1) gridX(end, end)], ...
     [gridY(end, 1) gridY(end, end)], ...
     'y-', 'LineWidth', 1.5);

plot([gridX(1, 1)   gridX(end, 1)], ...
     [gridY(1, 1)   gridY(end, 1)], ...
     'y-', 'LineWidth', 1.5);

plot([gridX(1, end) gridX(end, end)], ...
     [gridY(1, end) gridY(end, end)], ...
     'y-', 'LineWidth', 1.5);

title('Completed! Yellow dots indicate 64 × 64 channels');

% Pause here so the user can visually check whether the grid alignment is good.
pause;
close(gcf);

%% Step 4: Export grid coordinates to the MATLAB base workspace
% coords contains:
%   column 1: chip x-index, 0–63
%   column 2: chip y-index, 0–63
%   column 3: image x-coordinate in pixels
%   column 4: image y-coordinate in pixels

coords = [reshape(Xidx, [], 1), ...
          reshape(Yidx, [], 1), ...
          reshape(gridX, [], 1), ...
          reshape(gridY, [], 1)];

assignin('base', 'channelCoords', coords);

disp('Variable channelCoords = [chip_x chip_y image_x image_y] has been saved to the workspace');

%% Step 5: Load neuron classification and coordinate data
% This section reads the existing analysisResults structure from matPath.
% Expected fields:
%   S.neuronData.labels     : neuron response classification labels
%   S.neuronData.centres    : neuron centre coordinates in chip space, µm
%   S.neuronData.channelIDs : cluster/channel IDs

basename = extractBefore(imgPath, '.tif');
outTiff = [basename '_aligned.tif'];

S = load(matPath).analysisResults;

index     = S.neuronData.labels;      % 1 × N response labels
centre    = S.neuronData.centres;     % 2 × N neuron coordinates in µm
clusterID = S.neuronData.channelIDs;  % cluster/channel IDs

%% Convert neuron coordinates from chip space to approximate channel row/column
% centre is assumed to be in µm.
% Dividing by 60 converts µm coordinates into approximate channel indices,
% assuming 60 µm pitch.
%
% row and col are expected to range from 0 to 63.

row = round(centre(2, :) / 60);
col = round(centre(1, :) / 60);

% Keep only neurons whose mapped channel position falls inside the 64 × 64 grid.
valid = row >= 0 & row < 64 & col >= 0 & col < 64;

row = row(valid);
col = col(valid);
idx = index(valid);

%% Display the image again before drawing cluster overlays
figure('Name', 'Chip Grid Calibration', 'NumberTitle', 'off');
hAx = axes;
imshow(I, [], 'Parent', hAx);
title('Clusters overlayed');
hold on;

%% Step 6: Overlay responsive clusters on the aligned TIFF image
% Only labels 1–4 are drawn:
%   1 = ON
%   2 = OFF
%   3 = ON-OFF
%   4 = unconventional
%
% Each responsive cluster is drawn as a coloured rectangular outline.
% The cluster ID is also printed above the rectangle.

N = numel(idx);

% Estimate one channel width/height in image pixel units.
channelW = norm(ux);
channelH = norm(uy);

% Convert grayscale TIFF image to RGB so coloured annotations can be written.
Irgb = repmat(mat2gray(I), 1, 1, 3);

% Width of the rectangle border in pixels.
lineW = 2;

% Store virtual coordinates of labelled responsive clusters.
% Each row:
%   [clusterID, responseLabel, image_x, image_y]
clusterCoords = [];

for k = 1:N

    % Skip unresponsive neurons or labels outside the defined responsive classes.
    if idx(k) < 1 || idx(k) > 4
        continue;
    end

    %% Define colour according to response type
    switch idx(k)
        case 1
            color = [1 0 0];      % red: ON
        case 2
            color = [0 0 1];      % blue: OFF
        case 3
            color = [0 1 0];      % green: ON-OFF
        case 4
            color = [1 0.5 0];    % orange: unconventional
    end

    %% Find the top-left pixel coordinate of this neuron's corresponding channel
    % row/col are chip-space coordinates.
    % gridX/gridY are image pixel coordinates.
    %
    % The expression 64 - r0 flips the vertical axis to match the orientation
    % between chip coordinates and image coordinates.

    r0 = max(row(k), 0);
    c0 = max(col(k), 0);

    TLx = gridX(64 - r0, c0 + 1);
    TLy = gridY(64 - r0, c0 + 1);

    % Save cluster ID, response label, and image coordinate.
    clusterCoords = [clusterCoords; [clusterID(k), idx(k), TLx, TLy]];

    %% Draw a rectangular border around the cluster location
    % The rectangle size is set to approximately 2 channels wide/high.
    x1 = TLx;
    y1 = TLy;
    x2 = x1 + round(2 * channelW);
    y2 = y1 + round(2 * channelH);

    % Restrict the rectangle coordinates to the image boundaries.
    x1 = max(x1, 1);
    y1 = max(y1, 1);
    x2 = min(x2, size(I, 2));
    y2 = min(y2, size(I, 1));

    % Create a logical mask for the rectangle border only.
    % The inside of the rectangle is not filled.
    mask = false(size(I));

    % Top border
    mask(y1:y1 + lineW - 1, x1:x2) = true;

    % Bottom border
    mask(y2 - lineW + 1:y2, x1:x2) = true;

    % Left border
    mask(y1:y2, x1:x1 + lineW - 1) = true;

    % Right border
    mask(y1:y2, x2 - lineW + 1:x2) = true;

    %% Apply the colour to the RGB image
    for c = 1:3
        tmp = Irgb(:, :, c);
        tmp(mask) = color(c);
        Irgb(:, :, c) = tmp;
    end

    %% Add cluster ID text above the rectangle
    pos = [x1 + (x2 - x1) / 2, y1 - 12];

    Irgb = insertText(Irgb, pos, sprintf('%d', clusterID(k)), ...
                      'FontSize', 18, ...
                      'BoxOpacity', 0, ...
                      'AnchorPoint', 'Center', ...
                      'TextColor', color);
end

%% Save the annotated aligned image as a 16-bit TIFF
Irgb16 = uint16(Irgb * 65535);

imwrite(Irgb16, outTiff, ...
        'tiff', ...
        'Compression', 'none');

close(gcf);

%% Step 7: User selects a new reference centre
% This is usually the optic nerve head or another manually defined centre.
% The centre is selected in image pixel coordinates first, then converted
% back to chip-space coordinates in µm.

figure('Name', 'Chip Grid Calibration', 'NumberTitle', 'off');
hAx = axes;
imshow(I, [], 'Parent', hAx);

uiwait(msgbox({'Step 2', ...
               'Zoom as needed, then click the new reference center (yellow dot)', ...
               'Window will close automatically.'}, ...
               'Choose Center'));

hC = drawpoint('Color', 'c', ...
               'MarkerSize', 12, ...
               'Label', 'C');

centerPixel = hC.Position;

%% Convert selected image-pixel centre into chip-space µm coordinates
% Find the closest virtual grid channel to the clicked centre point.

Dg = hypot(gridX - centerPixel(1), ...
           gridY - centerPixel(2));

[~, idxMin] = min(Dg(:));

[rowC, colC] = ind2sub(size(Dg), idxMin);

% Convert nearest channel index into real chip coordinate.
% Multiplication by 60 assumes 60 µm inter-channel spacing.
%
% realCenter is stored as:
%   [x; y] in µm

realCenter = [colC - 1; ...
              63 - (rowC - 1)] * 60;

%% Step 8: Recalculate distances from all neurons to the selected centre
distances = sqrt((S.neuronData.centres(1, :) - realCenter(1)).^2 + ...
                 (S.neuronData.centres(2, :) - realCenter(2)).^2);

%% Step 9: Bin neurons by distance from the selected centre
% bin_edges defines distance bins from 0 to the maximum observed distance.
% The bin width is controlled by p.Results.BinWidth.

bin_edges = 0:p.Results.BinWidth: ...
            ceil(max(distances) / p.Results.BinWidth) * p.Results.BinWidth;

bin_centers = bin_edges(1:end-1) + p.Results.BinWidth / 2;

% Define responsive neurons.
% Labels:
%   1 = ON
%   2 = OFF
%   3 = ON-OFF
%   4 = unconventional
responsive_neurons = ismember(S.neuronData.labels, [1 2 3 4]);

% Calculate responsive percentage within each distance bin.
[response_percentage, ~, bin_stats] = calculate_binned_response( ...
    distances, ...
    responsive_neurons, ...
    bin_edges);

%% Step 10: Save updated analysis results
analysisResults = struct();

analysisResults.centre              = [realCenter(1), realCenter(2)];
analysisResults.distances           = distances;
analysisResults.bin_edges           = bin_edges;
analysisResults.bin_centers         = bin_centers;
analysisResults.response_percentage = response_percentage;
analysisResults.bin_stats           = bin_stats;
analysisResults.neuronData          = S.neuronData;
analysisResults.neuronVirtualCoords = clusterCoords;

save(matPath, 'analysisResults');

fprintf('Updated %s with new distances and binned responsiveness results.\n', matPath);

close(gcf);

end


function [response_percentage, valid_bins, bin_stats] = calculate_binned_response(distances, isresponsive, bin_edges)
% calculate_binned_response
%
% Calculates the percentage of responsive neurons in each distance bin.
%
% INPUTS:
%   distances    : 1 × N vector of neuron distances from the selected centre
%   isresponsive : 1 × N logical vector indicating whether each neuron is responsive
%   bin_edges    : vector defining distance bin edges
%
% OUTPUTS:
%   response_percentage : percentage of responsive neurons in each bin
%   valid_bins          : logical vector indicating bins containing at least one neuron
%   bin_stats           : structure containing:
%                           total_count
%                           responsive_count
%                           mean_distance

    bin_stats = struct();

    response_percentage = zeros(1, length(bin_edges) - 1);
    valid_bins = false(1, length(bin_edges) - 1);

    for i = 1:length(bin_edges) - 1

        % Find neurons whose distance falls within the current bin.
        in_bin = distances >= bin_edges(i) & distances < bin_edges(i + 1);

        % Store basic bin-level statistics.
        bin_stats(i).total_count = sum(in_bin);
        bin_stats(i).responsive_count = sum(isresponsive(in_bin));
        bin_stats(i).mean_distance = mean(distances(in_bin));

        % Calculate responsive percentage only for non-empty bins.
        if bin_stats(i).total_count > 0
            response_percentage(i) = ...
                bin_stats(i).responsive_count / bin_stats(i).total_count * 100;

            valid_bins(i) = true;
        else
            response_percentage(i) = NaN;
        end
    end
end