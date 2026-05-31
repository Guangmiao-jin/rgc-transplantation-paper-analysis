function chipGridAlign(imgPath, matPath,varargin)

p = inputParser;
addParameter(p, 'BinWidth', 50, @isnumeric); % 
addParameter(p, 'SmoothWindow', 3, @isnumeric); % 
parse(p, varargin{:});

%% read the tiff file
I = imread(imgPath);
%imSz = size(I);
figure('Name','Chip Grid Calibration','NumberTitle','off');
hAx = axes; imshow(I,[],'Parent',hAx); title('Step 1: Zoom & pick the FOUR corner channels, top left -> top right -> bottom left -> bottom right');

%% 1. extract the centres of four corners
cornerPt = [];                 % 4×2
for k = 1:4
    h = drawpoint('Color','r','MarkerSize',10,'Label',sprintf('%d',k));
    cornerPt(k,:) = h.Position;
end
%% 2. calculate the grids' distances
% assume cornerPt goes in the following order：c00, c0N, cM0, cMN
c00 = cornerPt(1,:);         % top left
c0N = cornerPt(2,:);         % top right
cM0 = cornerPt(3,:);         % bottom left

% calculate the distance between each channels horizontally and vertically
ux = (c0N - c00)/(64-1);     % horizontal distance
uy = (cM0 - c00)/(64-1);     % vertical distance

[Xidx,Yidx] = meshgrid(0:63,0:63); 
gridX = c00(1) + ux(1)*Xidx + uy(1)*Yidx;
gridY = c00(2) + ux(2)*Xidx + uy(2)*Yidx;

%% 3. visualization
hold on;
scatter(gridX(:),gridY(:),15,'y','filled');   % grids' points
% draw the borders
plot([gridX(1,1) gridX(1,end)],[gridY(1,1) gridY(1,end)],'y-','LineWidth',1.5);
plot([gridX(end,1) gridX(end,end)],[gridY(end,1) gridY(end,end)],'y-','LineWidth',1.5);
plot([gridX(1,1) gridX(end,1)],[gridY(1,1) gridY(end,1)],'y-','LineWidth',1.5);
plot([gridX(1,end) gridX(end,end)],[gridY(1,end) gridY(end,end)],'y-','LineWidth',1.5);
title('Completed！Yellow dosts indicate 64×64 channels'); pause;
close(gcf);

%% 4. coordinates output
coords = [reshape(Xidx,[],1) reshape(Yidx,[],1) reshape(gridX,[],1) reshape(gridY,[],1)];
assignin('base','channelCoords',coords);  % saved to workspace
disp('Variable channelCoords = [row col x y] has been saved to the workspace');

%% 5. read _distance.mat file and annotate responsive neurons on tiff image
basename = extractBefore(imgPath,'.tif');
outTiff = [basename '_aligned.tif'];
%outTiff2 = [basename '_aligned2.tif'];

S = load(matPath).analysisResults;          
index  = S.neuronData.labels;       % 1×N
centre = S.neuronData.centres;      % 2×N
clusterID = S.neuronData.channelIDs;
N      = numel(index);

% centre → channel (row,col)
row = round( centre(2,:)/60 );    % 0–63
col = round( centre(1,:)/60 );    % 0–63
valid = row>=0 & row<64 & col>=0 & col<64;  
row = row(valid);  col = col(valid);  idx = index(valid);

figure('Name','Chip Grid Calibration','NumberTitle','off');
hAx = axes; imshow(I,[],'Parent',hAx)
title('Clusters overlayed');
hold on;

%%%%%%%% test %%%%%%%%

N = numel(idx);

channelW = norm(ux);   % single channel width
channelH = norm(uy);   % single channel height

Irgb = repmat(mat2gray(I),1,1,3);   % double 0 1 RGB
lineW = 2;                          % pixel width
clusterCoords = [];
for k = 1:N
    if idx(k) < 1 || idx(k) > 4,  continue;  end

    % —— colour table ——————————————————————————
    switch idx(k)
        case 1,  color = [1 0 0];        % red
        case 2,  color = [0 0 1];        % blue
        case 3,  color = [0 1 0];        % green
        case 4,  color = [1 0.5 0];      % orange
    end

    % —— count the pixel coords from the top left ——————————————
    r0 = max(row(k),0); c0 = max(col(k),0);
    TLx = gridX(64-r0,c0+1); TLy = gridY(64-r0,c0+1);
    clusterCoords = [clusterCoords; [clusterID(k) idx(k) TLx TLy]]; 
    
    % —— border mask（only border，no filling） ——————————
    x1 = TLx; y1 = TLy;
    x2 = x1 + round(2*channelW); y2 = y1 + round(2*channelH);

    % set the border axies
    x1 = max(x1,1);  y1 = max(y1,1);
    x2 = min(x2,size(I,2)); 
    y2 = min(y2,size(I,1));

    mask = false(size(I));
    % top bottom left right
    mask(y1:y1+lineW-1 , x1:x2) = true;           % top
    mask(y2-lineW+1:y2, x1:x2) = true;            % bottom
    mask(y1:y2         , x1:x1+lineW-1) = true;   % left
    mask(y1:y2         , x2-lineW+1:x2) = true;   % right

    % α blend：for loop
    for c = 1:3
        tmp              = Irgb(:,:,c);
        tmp(mask)        = color(c);
        Irgb(:,:,c)      = tmp;
    end

    % —— text —— top of the rectangle ——————
    pos = [x1 + (x2-x1)/2 , y1 - 12];   % 12 px shift top，can adjust freely
    Irgb = insertText(Irgb, pos, sprintf('%d',clusterID(k)), ...
                      'FontSize', 18,         ...
                      'BoxOpacity', 0,        ...  % transparent box
                      'AnchorPoint','Center', ...
                      'TextColor', color);    % text according to the classification
end

% —— 16 bit TIFF ————————————————————
Irgb16 = uint16(Irgb * 65535);
imwrite(Irgb16, outTiff, 'tiff', 'Compression', 'none');
close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% user select the centre
figure('Name','Chip Grid Calibration','NumberTitle','off');
hAx = axes; imshow(I,[],'Parent',hAx)
uiwait(msgbox({'Step 2','Zoom as needed, then click the new reference center (yellow dot)', ...
               'Window will close automatically.'},'Choose Center'));
hC = drawpoint('Color','c','MarkerSize',12,'Label','C');
centerPixel = hC.Position;

% pixel coordiantes -> real um in chip
% find nearest channel (rowC,colC) → multiply by 60
Dg = hypot(gridX-centerPixel(1), gridY-centerPixel(2));
[~,idxMin] = min(Dg(:));
[rowC,colC] = ind2sub(size(Dg),idxMin);
realCenter = [colC-1; 63-(rowC-1)] * 60;        % (µm)

% new distances and new binned results
distances = sqrt((S.neuronData.centres(1,:) - realCenter(1)).^2 + (S.neuronData.centres(2,:) - realCenter(2)).^2);
     
bin_edges = 0:p.Results.BinWidth:ceil(max(distances)/p.Results.BinWidth)*p.Results.BinWidth;
bin_centers = bin_edges(1:end-1) + p.Results.BinWidth/2;
        
responsive_neurons = ismember(S.neuronData.labels, [1 2 3 4]); % ON, OFF, ON-OFF, unconventional
[response_percentage, ~, bin_stats] = calculate_binned_response(distances, responsive_neurons, bin_edges);
        
     
analysisResults = struct();
analysisResults.centre = [realCenter(1), realCenter(2)];
analysisResults.distances = distances;
analysisResults.bin_edges = bin_edges;
analysisResults.bin_centers = bin_centers;
analysisResults.response_percentage = response_percentage;
analysisResults.bin_stats = bin_stats; 
analysisResults.neuronData = S.neuronData; 
analysisResults.neuronVirtualCoords = clusterCoords;

save(matPath,'analysisResults');
fprintf('Updated %s with new dist (& binID if available).\n',matPath);

close(gcf);   % close all;

end

function [response_percentage, valid_bins, bin_stats] = calculate_binned_response(distances, isresponsive, bin_edges)
    bin_stats = struct();
    response_percentage = zeros(1, length(bin_edges)-1);
    valid_bins = false(1, length(bin_edges)-1);
    
    for i = 1:length(bin_edges)-1
        in_bin = (distances >= bin_edges(i)) & (distances < bin_edges(i+1));
        bin_stats(i).total_count = sum(in_bin);
        bin_stats(i).responsive_count = sum(isresponsive(in_bin));
        bin_stats(i).mean_distance = mean(distances(in_bin));
        
        if bin_stats(i).total_count > 0
            response_percentage(i) = bin_stats(i).responsive_count / bin_stats(i).total_count * 100;
            valid_bins(i) = true;
        else
            response_percentage(i) = NaN;
        end
    end
end