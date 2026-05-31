function chipGridAlign1(imgPath, matPath,varargin)

p = inputParser;
addParameter(p, 'BinWidth', 50, @isnumeric); % 
addParameter(p, 'SmoothWindow', 3, @isnumeric); % 
parse(p, varargin{:});

%% read the tiff file
I = imread(imgPath);
[H,W] = size(I);
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
S = load(matPath).analysisResults;          
index  = S.neuronData.labels;       % 1×N
centre = S.neuronData.centres;      % 2×N
clusterID = S.neuronData.channelIDs;

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

I8     = im2uint8(rescale(I,0,1));        % H×W  → uint8
Irgb8  = repmat(I8,1,1,3);                  % H×W×3
colors = uint8([255 0 0; 0 0 255; 0 255 0; 255 128 0]);  % r、b、g、orange
lineW  = 10;
fontSz = 50;                                

clusterCoords = [];
for k = 1:N
    if idx(k) < 1 || idx(k) > 4,  continue;  end
    % —— count the pixel coords from the top left ——————————————
        % ------ 4.1 calculate matrix pixel coordinates  ------
    TLx = gridX(64-row(k),col(k)+1);
    TLy = gridY(64-row(k),col(k)+1);
    x1  = max(1, round(TLx));
    y1  = max(1, round(TLy));
    x2  = min(W, round(TLx+2*channelW));
    y2  = min(H, round(TLy+2*channelH));


    % ------ 4.2 border colours  ------
        colRGB = colors(idx(k),:);
    for c = 1:3
        Irgb8(y1:y1+lineW-1 , x1:x2 , c) = colRGB(c);  % top
        Irgb8(y2-lineW+1:y2 , x1:x2 , c) = colRGB(c);  % bottom
        Irgb8(y1:y2 , x1:x1+lineW-1 , c) = colRGB(c);  % left
        Irgb8(y1:y2 , x2-lineW+1:x2 , c) = colRGB(c);  % right
    end

   txtStr  = sprintf('%d',clusterID(k));
    pad     = 4;                       
    % estimate the font, font≈ width pixel/1.5
    patchW  = round(fontSz*length(txtStr)*0.75) + 2*pad;
    patchH  = fontSz + 2*pad;
    patch   = zeros(patchH,patchW,'uint8');
    patchRGB= repmat(patch,1,1,3);     % patchH×patchW×3

    %    3.2 write on the patch by insertText
    patchRGB = insertText(patchRGB,[patchW/2 patchH/2],txtStr, ...
                'FontSize',fontSz,'BoxOpacity',0,'AnchorPoint','Center', ...
                'TextColor',colRGB);

    %    3.3 use α-blend patch in the bigger picture
    yTxt = max(1, y1 - patchH - 2);    % patch above the rectangle
    xTxt = round((x1+x2)/2 - patchW/2);
    xTxt = max(1, min(W-patchW+1, xTxt));      % strip the margin 

    subBlock = Irgb8( yTxt:yTxt+patchH-1 , xTxt:xTxt+patchW-1 , : );

% Non black regions（RGB!=0): allow colour
maskP = any(patchRGB,3);                       % patchH×patchW logical
mask3 = repmat(maskP,1,1,3);                   % expand to 3 channels

% use patchRGB to cover subBlock
subBlock(mask3) = patchRGB(mask3);

% big picture
Irgb8( yTxt:yTxt+patchH-1 , xTxt:xTxt+patchW-1 , : ) = subBlock;


end

% —— 16 bit TIFF ————————————————————
outTiff = [extractBefore(imgPath,'.tif') '_aligned.tif'];
imwrite(Irgb8, outTiff, 'tiff', ...
        'Compression','lzw', ...
        'RowsPerStrip', 64);              % >4 GB, use BigTIFF

fprintf('8-bit RGB overlay saved to %s (≈ %0.1f MB)\n', ...
        outTiff, dir(outTiff).bytes/1e6);


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

