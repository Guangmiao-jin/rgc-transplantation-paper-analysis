function [T, T1] = collectDeltaFR_transplant_sham(rootDir)

% rootDir example:
% '/Users/xxx/MEA_analysis/'

allRows = {};

allLatency = {};

stages = {'early degeneration', 'late degeneration'};
groups = {'NRL', 'SHAM'};

for s = 1:numel(stages)
    stageName = stages{s};

    for g = 1:numel(groups)
        groupName = groups{g};

        % For NRL, only use behavioural positive + ephys positive folder
        if strcmp(groupName, 'NRL')
            searchDir = fullfile(rootDir, stageName, groupName, ...
                'behavioural positive and ephys positive');
        else
            searchDir = fullfile(rootDir, stageName, groupName);
        end

        psthFiles = dir(fullfile(searchDir, '**', '*_psth.mat'));

        for f = 1:numel(psthFiles)

            psthPath = fullfile(psthFiles(f).folder, psthFiles(f).name);

            % Find matching totalneuronsV2 file
            totalFiles = dir(fullfile(psthFiles(f).folder, '*_totalneuronsV2.mat'));

            if isempty(totalFiles)
                warning('No totalneuronsV2 file found for %s', psthPath);
                continue;
            end

            totalPath = fullfile(totalFiles(f).folder, totalFiles(f).name);

            P1 = load(psthPath);
            P = P1.pltcurve;
            N1 = load(totalPath);
            N = N1.totalneurons;

            % Adjust these variable names if needed
            PSTH_binEdges = P.PSTH_binEdges;
            trialPSTHs = P.trialPSTHs;      
            num_OnNeurons = N.num_OnNeurons;
            num_OffNeurons = N.num_OffNeurons;
            num_OnOffNeurons = N.num_OnOffNeurons;
            num_unconventional = N.num_unconventional;

            % 0-based indexing in totalneuronsV2; +1 for MATLAB indexing
            onIDs = unique([
                num_OnNeurons.trans(:);
                num_OnNeurons.sus(:)
            ]) + 1;

            offIDs = unique([
                num_OffNeurons.trans(:);
                num_OffNeurons.sus(:)
            ]) + 1;

            onoffIDs = num_OnOffNeurons+ 1;
            unconvIDs = num_unconventional+ 1;

            % Remove invalid IDs just in case
            nUnits = numel(trialPSTHs);
            onIDs = onIDs(onIDs >= 1 & onIDs <= nUnits);
            offIDs = offIDs(offIDs >= 1 & offIDs <= nUnits);
            onoffIDs = onoffIDs(onoffIDs >= 1 & onoffIDs <= nUnits);
            unconvIDs = unconvIDs(unconvIDs >= 1 & unconvIDs <= nUnits);

            binWidth = mean(diff(PSTH_binEdges));

            binStarts = PSTH_binEdges(1:end-1);
            binEnds = PSTH_binEdges(2:end);

            onBins = binStarts >= -2 & binEnds <= 0;
            offBins = binStarts >= 0 & binEnds <= 2;

            % Process ON cells
            for i = 1:numel(onIDs)
                unitID = onIDs(i);

                for cond = 1:3
                    psth = trialPSTHs{unitID}{cond};  % 91 x 160

                    meanPSTH = mean(psth, 1, 'omitnan');
                    fr = meanPSTH ./ binWidth;

                    FR_on = mean(fr(onBins), 'omitnan');
                    FR_off = mean(fr(offBins), 'omitnan');

                    deltaFR = FR_on - FR_off;

                    allRows(end+1, :) = {
                        stageName, groupName, 'ON', cond, unitID, ...
                        deltaFR, FR_on, FR_off, psthFiles(f).name
                    };
                    [onsetLatency, peakLatency, peakFR, threshold] = computeResponseLatency( psth, PSTH_binEdges, [-2 0], [1 2], -2);
                    allLatency(end+1,: )= {stageName, groupName, 'ON', cond, unitID, ...
                       onsetLatency, peakLatency, peakFR, threshold, psthFiles(f).name};
                end
            end

            % Process OFF cells
            for i = 1:numel(offIDs)
                unitID = offIDs(i);

                for cond = 1:3
                    psth = trialPSTHs{unitID}{cond};  % 91 x 160

                    meanPSTH = mean(psth, 1, 'omitnan');
                    fr = meanPSTH ./ binWidth;

                    FR_on = mean(fr(onBins), 'omitnan');
                    FR_off = mean(fr(offBins), 'omitnan');

                    deltaFR = FR_off - FR_on;

                    allRows(end+1, :) = {
                        stageName, groupName, 'OFF', cond, unitID, ...
                        deltaFR, FR_on, FR_off, psthFiles(f).name
                    };
                    [onsetLatency, peakLatency, peakFR, threshold] = computeResponseLatency( psth, PSTH_binEdges, [0 2], [-1 0], 0);
                       allLatency(end+1,: )= {stageName, groupName, 'OFF', cond, unitID, ...
                       onsetLatency, peakLatency, peakFR, threshold, psthFiles(f).name};

                end
            end

            % ON-OFF cells
            for i = 1:numel(onoffIDs)
                unitID = onoffIDs(i);

                for cond = 1:3
                    psth = trialPSTHs{unitID}{cond};  % 91 x 160

                    meanPSTH = mean(psth, 1, 'omitnan');
                    fr = meanPSTH ./ binWidth;

                    FR_on = mean(fr(onBins), 'omitnan');
                    FR_off = mean(fr(offBins), 'omitnan');

                    deltaFR = abs(FR_on - FR_off);

                    allRows(end+1, :) = {
                        stageName, groupName, 'ONOFF', cond, unitID, ...
                        deltaFR, FR_on, FR_off, psthFiles(f).name
                    };
                    % =====================================================
                    % Latency：use whole stimulus window
                    % because ON-OFF type has peak at ON/OFF timing
                    % use the widest window to capture any response
                    % =====================================================
                    [onsetLatency_on, peakLatency_on, peakFR_on, threshold_on] = computeResponseLatency( psth, PSTH_binEdges, [-2 -2.5], [1 2], -2);
                    [onsetLatency_off, peakLatency_off, peakFR_off, threshold_off] = computeResponseLatency( psth, PSTH_binEdges, [0 0.5], [-1 0], 0);
                    if ~isnan(onsetLatency_on) && ~isnan(onsetLatency_off)
                        onsetLatency = (onsetLatency_on + onsetLatency_off)/2;
                        peakLatency = (peakLatency_on + peakLatency_off)/2;
                        peakFR = (peakFR_on + peakFR_off)/2;
                        threshold = (threshold_on + threshold_off)/2;
                    elseif ~isnan(onsetLatency_on) && isnan(onsetLatency_off)
                        onsetLatency = onsetLatency_on;
                        peakLatency = peakLatency_on;
                        peakFR = peakFR_on;
                        threshold = threshold_on;
                     elseif isnan(onsetLatency_on) && ~isnan(onsetLatency_off)
                        onsetLatency = onsetLatency_off;
                        peakLatency = peakLatency_off;
                        peakFR = peakFR_off;
                        threshold = threshold_off;
                    else
                        onsetLatency = NaN;
                        peakLatency = NaN;
                        peakFR = NaN;
                        threshold = NaN;
                    end
                    allLatency(end+1, :) = {stageName, groupName, 'ONOFF', cond, unitID, ...
                     onsetLatency, peakLatency, peakFR, threshold, ...
                     psthFiles(f).name};
                end
            end

             % unconventional cells
            for i = 1:numel(unconvIDs)
                unitID = unconvIDs(i);
                for cond = 1:3
                    psth = trialPSTHs{unitID}{cond};
                    meanPSTH = mean(psth, 1, 'omitnan');
                    fr = meanPSTH ./ binWidth;
                    FR_on  = mean(fr(onBins),  'omitnan');
                    FR_off = mean(fr(offBins), 'omitnan');
                    deltaFR = abs(FR_on - FR_off);
                    allRows(end+1, :) = {stageName, groupName, 'unconventional', cond, unitID, ...
                    deltaFR, FR_on, FR_off, psthFiles(f).name};

                    % =====================================================
                    % Latency：use whole stimulus window
                    % because unconventional response does not have certain ON/OFF timing
                    % use the widest window to capture any response
                    % =====================================================
                    [onsetLatency_on, peakLatency_on, peakFR_on, threshold_on] = computeResponseLatency( psth, PSTH_binEdges, [-2 -2.5], [1 2], -2);
                    [onsetLatency_off, peakLatency_off, peakFR_off, threshold_off] = computeResponseLatency( psth, PSTH_binEdges, [0 0.5], [-1 0], 0);
                    if ~isnan(onsetLatency_on) && ~isnan(onsetLatency_off)
                        onsetLatency = (onsetLatency_on + onsetLatency_off)/2;
                        peakLatency = (peakLatency_on + peakLatency_off)/2;
                        peakFR = (peakFR_on + peakFR_off)/2;
                        threshold = (threshold_on + threshold_off)/2;
                    elseif ~isnan(onsetLatency_on) && isnan(onsetLatency_off)
                        onsetLatency = onsetLatency_on;
                        peakLatency = peakLatency_on;
                        peakFR = peakFR_on;
                        threshold = threshold_on;
                     elseif isnan(onsetLatency_on) && ~isnan(onsetLatency_off)
                        onsetLatency = onsetLatency_off;
                        peakLatency = peakLatency_off;
                        peakFR = peakFR_off;
                        threshold = threshold_off;
                    else
                        onsetLatency = NaN;
                        peakLatency = NaN;
                        peakFR = NaN;
                        threshold = NaN;
                    end
                    allLatency(end+1, :) = {stageName, groupName, 'unconventional', cond, unitID, ...
                     onsetLatency, peakLatency, peakFR, threshold, ...
                     psthFiles(f).name};
               end
           end
        end
    end
end

T = cell2table(allRows, ...
    'VariableNames', {'Stage', 'Group', 'CellType', 'LightCondition', ...
    'UnitID', 'DeltaFR', 'FR_ON', 'FR_OFF', 'FileID'});

T.Stage = categorical(T.Stage);
T.Group = categorical(T.Group);
T.CellType = categorical(T.CellType);
T.LightCondition = categorical(T.LightCondition);

T1 = cell2table(allLatency, ...
    'VariableNames', {'Stage', 'Group', 'CellType', 'LightCondition', ...
    'UnitID','onsetLatency', 'peakLatency', 'peakFR', 'threshold','FileID'} );
T1.Stage = categorical(T1.Stage);
T1.Group = categorical(T1.Group);
T1.CellType = categorical(T1.CellType);
T1.LightCondition = categorical(T1.LightCondition);

end

function [onsetLatency, peakLatency, peakFR, threshold] = computeResponseLatency( ...
    psth, binEdges, responseWindow, baselineWindow, eventTime)
% psth: trials x time bins, e.g. 91 x 160
% binEdges: 1 x 161
% responseWindow: e.g. [-2 0] for ON, [0 2] for OFF
% baselineWindow: e.g. [1 2] for ON, [-1 0] for OFF
% eventTime: -2 for ON onset, 0 for OFF onset

    binWidth = mean(diff(binEdges));
    binCenters = binEdges(1:end-1) + binWidth/2;

    % Trial-averaged firing rate
    meanPSTH = mean(psth, 1, 'omitnan');
    fr = meanPSTH ./ binWidth;

    % Optional smoothing
    %smoothBins = 3;
    %frSmooth = movmean(fr, smoothBins, 'omitnan');

    responseBins = binCenters >= responseWindow(1) & binCenters < responseWindow(2);
    baselineBins = binCenters >= baselineWindow(1) & binCenters < baselineWindow(2);

    baselineFR = fr(baselineBins);

    if isempty(baselineFR) || all(isnan(baselineFR))
        onsetLatency = NaN;
        peakLatency = NaN;
        peakFR = NaN;
        threshold = NaN;
        return;
    end

    baselineMean = mean(baselineFR, 'omitnan');
    baselineSD = std(baselineFR, 'omitnan');

    % Threshold can be adjusted: 2 SD is sensitive, 3 SD is stricter
    threshold = baselineMean + 2 * baselineSD;

    respFR = fr(responseBins);
    respTimes = binCenters(responseBins);

    if isempty(respFR) || all(isnan(respFR))
        onsetLatency = NaN;
        peakLatency = NaN;
        peakFR = NaN;
        return;
    end

    % Require consecutive bins above threshold to avoid single-bin noise
    aboveThreshold = respFR > threshold;
    minConsecutiveBins = 2;

    onsetIdx = NaN;

    for k = 1:(numel(aboveThreshold) - minConsecutiveBins + 1)
        if all(aboveThreshold(k:k+minConsecutiveBins-1))
            onsetIdx = k;
            break;
        end
    end

    if isnan(onsetIdx)
        onsetLatency = NaN;
    else
        onsetLatency = respTimes(onsetIdx) - eventTime;
    end

    % Peak latency
    [peakFR, peakIdx] = max(respFR, [], 'omitnan');
    peakLatency = respTimes(peakIdx) - eventTime;
end