function [T, T1] = collectDeltaFR_transplant_sham(rootDir)
% collectDeltaFR_transplant_sham
%
% Collects firing-rate and latency metrics from transplanted and sham
% retinal MEA recordings.
%
% This function searches through the transplantation analysis folder,
% loads each recording's PSTH file and corresponding totalneuronsV2 file,
% then extracts response metrics for classified responsive RGCs.
%
% For each responsive neuron and each light condition, the function
% calculates:
%   1. Delta firing rate:
%          DeltaFR = abs(peak firing rate during response window ...
%                        - baseline mean firing rate)
%
%   2. Peak firing rate during the response window
%
%   3. Baseline mean firing rate
%
%   4. Response latency metrics:
%          onsetLatency
%          peakLatency
%          peakFR
%          threshold
%
% OUTPUTS:
%   T  : table containing DeltaFR, PeakStimFR, and BaselineMeanFR
%
%   T1 : table containing onset latency, peak latency, peak firing rate,
%        and threshold
%
% EXPECTED FOLDER STRUCTURE:
%   rootDir/
%       early degeneration/
%           NRL/
%               behavioural positive and ephys positive/
%                   *_psth.mat
%                   *_totalneuronsV2.mat
%           SHAM/
%               *_psth.mat
%               *_totalneuronsV2.mat
%
%       late degeneration/
%           NRL/
%               behavioural positive and ephys positive/
%                   *_psth.mat
%                   *_totalneuronsV2.mat
%           SHAM/
%               *_psth.mat
%               *_totalneuronsV2.mat
%
% NOTE:
%   For the NRL transplanted group, only the
%   'behavioural positive and ephys positive' folder is analysed.
%
%   For the SHAM group, all PSTH files under the SHAM folder are analysed.
%
% rootDir example:
%   rootDir = '/Users/xxx/MEA_analysis/';

%% Initialise output containers
% allRows stores firing-rate measurements.
% allLatency stores latency measurements.

allRows = {};
allLatency = {};

%% Define experimental stages and groups to analyse
stages = {'early degeneration', 'late degeneration'};
groups = {'NRL', 'SHAM'};

%% Loop through degeneration stages
for s = 1:numel(stages)

    stageName = stages{s};

    %% Loop through treatment groups
    for g = 1:numel(groups)

        groupName = groups{g};

        %% Define folder to search
        % For NRL transplanted retinas, only analyse recordings from the
        % behavioural-positive and ephys-positive subgroup.
        %
        % For SHAM retinas, search directly within the SHAM folder.

        if strcmp(groupName, 'NRL')
            searchDir = fullfile(rootDir, stageName, groupName, ...
                'behavioural positive and ephys positive');
        else
            searchDir = fullfile(rootDir, stageName, groupName);
        end

        %% Find all PSTH files recursively under the selected folder
        psthFiles = dir(fullfile(searchDir, '**', '*_psth.mat'));

        %% Loop through all PSTH files
        for f = 1:numel(psthFiles)

            psthPath = fullfile(psthFiles(f).folder, psthFiles(f).name);

            %% Find the matching totalneuronsV2 file in the same folder
            % The totalneuronsV2 file contains classified neuron indices:
            %   ON transient/sustained
            %   OFF transient/sustained
            %   ON-OFF
            %   unconventional

            [~, psthBase, ~] = fileparts(psthFiles(f).name);
            psthBase = erase(psthBase, '_psth');   % strip the _psth suffix
            totalName = [psthBase '_totalneuronsV2.mat'];
            totalPath = fullfile(psthFiles(f).folder, totalName);
            if ~isfile(totalPath)
                warning('No matching totalneuronsV2 file for %s', psthPath);
                continue;
            end
            %% Load PSTH and neuron classification data
            P1 = load(psthPath);
            P = P1.pltcurve;

            N1 = load(totalPath);
            N = N1.totalneurons;

            %% Extract required variables
            PSTH_binEdges = P.PSTH_binEdges;
            trialPSTHs = P.trialPSTHs;

            num_OnNeurons       = N.num_OnNeurons;
            num_OffNeurons      = N.num_OffNeurons;
            num_OnOffNeurons    = N.num_OnOffNeurons;
            num_unconventional  = N.num_unconventional;

            %% Convert classified neuron IDs to MATLAB indexing
            % totalneuronsV2 stores neuron IDs using 0-based indexing.
            % MATLAB uses 1-based indexing, so add 1 before using them
            % to access trialPSTHs.

            onIDs = unique([
                num_OnNeurons.trans(:);
                num_OnNeurons.sus(:)
            ]) + 1;

            offIDs = unique([
                num_OffNeurons.trans(:);
                num_OffNeurons.sus(:)
            ]) + 1;

            onoffIDs = num_OnOffNeurons + 1;
            unconvIDs = num_unconventional + 1;

            %% Remove invalid neuron IDs
            % This prevents indexing errors if any classified ID is outside
            % the available trialPSTHs range.

            nUnits = numel(trialPSTHs);

            onIDs     = onIDs(onIDs >= 1 & onIDs <= nUnits);
            offIDs    = offIDs(offIDs >= 1 & offIDs <= nUnits);
            onoffIDs  = onoffIDs(onoffIDs >= 1 & onoffIDs <= nUnits);
            unconvIDs = unconvIDs(unconvIDs >= 1 & unconvIDs <= nUnits);

            %% Define PSTH timing information
            % PSTH_binEdges define the time bin edges.
            % trialPSTHs{unitID}{cond} is expected to be trials × time bins,
            % for example 91 × 160.
            %
            % In this dataset:
            %   ON response window  = -2 to 0 s
            %   OFF response window =  0 to 2 s

            binWidth = mean(diff(PSTH_binEdges));

            binStarts = PSTH_binEdges(1:end-1);
            binEnds   = PSTH_binEdges(2:end);

            onBins  = binStarts >= -2 & binEnds <= 0;
            offBins = binStarts >= 0  & binEnds <= 2;

            %% Process ON cells
            % For ON cells:
            %   response window = [-2, 0]
            %   baseline window = [1, 2]
            %   event time      = -2
            %
            % The event time is used so latency is calculated relative to
            % light onset.

            for i = 1:numel(onIDs)

                unitID = onIDs(i);

                for cond = 1:3

                    psth = trialPSTHs{unitID}{cond};

                    [deltaFR, peakStimFR, baselineMeanFR] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [-2 0], [1 2]);

                    allRows(end+1, :) = { ...
                        stageName, groupName, 'ON', cond, unitID, ...
                        deltaFR, peakStimFR, baselineMeanFR, ...
                        psthFiles(f).name};

                    [onsetLatency, peakLatency, peakFR, threshold] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [-2 0], [1 2], -2);

                    allLatency(end+1, :) = { ...
                        stageName, groupName, 'ON', cond, unitID, ...
                        onsetLatency, peakLatency, peakFR, threshold, ...
                        psthFiles(f).name};
                end
            end

            %% Process OFF cells
            % For OFF cells:
            %   response window = [0, 2]
            %   baseline window = [-1, 0]
            %   event time      = 0
            %
            % The event time is used so latency is calculated relative to
            % light offset.

            for i = 1:numel(offIDs)

                unitID = offIDs(i);

                for cond = 1:3

                    psth = trialPSTHs{unitID}{cond};

                    [deltaFR, peakStimFR, baselineMeanFR] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [0 2], [-1 0]);

                    allRows(end+1, :) = { ...
                        stageName, groupName, 'OFF', cond, unitID, ...
                        deltaFR, peakStimFR, baselineMeanFR, ...
                        psthFiles(f).name};

                    [onsetLatency, peakLatency, peakFR, threshold] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [0 2], [-1 0], 0);

                    allLatency(end+1, :) = { ...
                        stageName, groupName, 'OFF', cond, unitID, ...
                        onsetLatency, peakLatency, peakFR, threshold, ...
                        psthFiles(f).name};
                end
            end

            %% Process ON-OFF cells
            % ON-OFF cells can respond at both light onset and light offset.
            %
            % For DeltaFR:
            %   The function calculates ON-window DeltaFR and OFF-window
            %   DeltaFR separately, then keeps the larger response.
            %
            % For latency:
            %   ON and OFF latency are calculated separately. If both are
            %   detected, the mean value is used. If only one is detected,
            %   that value is kept.

            for i = 1:numel(onoffIDs)

                unitID = onoffIDs(i);

                for cond = 1:3

                    psth = trialPSTHs{unitID}{cond};

                    %% Calculate ON-window DeltaFR
                    [deltaFR_on, peakStimFR_on, baselineMeanFR_on] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [-2 0], [1 2]);

                    %% Calculate OFF-window DeltaFR
                    [deltaFR_off, peakStimFR_off, baselineMeanFR_off] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [0 2], [-1 0]);

                    %% Keep the larger ON or OFF response
                    if deltaFR_on >= deltaFR_off
                        deltaFR = deltaFR_on;
                        peakStimFR = peakStimFR_on;
                        baselineMeanFR = baselineMeanFR_on;
                    else
                        deltaFR = deltaFR_off;
                        peakStimFR = peakStimFR_off;
                        baselineMeanFR = baselineMeanFR_off;
                    end

                    allRows(end+1, :) = { ...
                        stageName, groupName, 'ONOFF', cond, unitID, ...
                        deltaFR, peakStimFR, baselineMeanFR, ...
                        psthFiles(f).name};

                    %% Calculate ON and OFF latency separately
                    [onsetLatency_on, peakLatency_on, peakFR_on, threshold_on] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [-2 0], [1 2], -2);

                    [onsetLatency_off, peakLatency_off, peakFR_off, threshold_off] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [0 2], [-1 0], 0);

                    %% Combine ON and OFF latency estimates
                    if ~isnan(onsetLatency_on) && ~isnan(onsetLatency_off)

                        onsetLatency = (onsetLatency_on + onsetLatency_off) / 2;
                        peakLatency  = (peakLatency_on  + peakLatency_off)  / 2;
                        peakFR       = (peakFR_on       + peakFR_off)       / 2;
                        threshold    = (threshold_on    + threshold_off)    / 2;

                    elseif ~isnan(onsetLatency_on) && isnan(onsetLatency_off)

                        onsetLatency = onsetLatency_on;
                        peakLatency  = peakLatency_on;
                        peakFR       = peakFR_on;
                        threshold    = threshold_on;

                    elseif isnan(onsetLatency_on) && ~isnan(onsetLatency_off)

                        onsetLatency = onsetLatency_off;
                        peakLatency  = peakLatency_off;
                        peakFR       = peakFR_off;
                        threshold    = threshold_off;

                    else

                        onsetLatency = NaN;
                        peakLatency  = NaN;
                        peakFR       = NaN;
                        threshold    = NaN;
                    end

                    allLatency(end+1, :) = { ...
                        stageName, groupName, 'ONOFF', cond, unitID, ...
                        onsetLatency, peakLatency, peakFR, threshold, ...
                        psthFiles(f).name};
                end
            end

            %% Process unconventional cells
            % Unconventional cells do not fit canonical ON, OFF, or ON-OFF
            % response categories.
            %
            % For DeltaFR:
            %   Both ON and OFF windows are tested, and the larger response
            %   is retained.
            %
            % For latency:
            %   ON and OFF latency are estimated separately. This is a
            %   pragmatic approach because unconventional cells may show
            %   their strongest response at either timing.

            for i = 1:numel(unconvIDs)

                unitID = unconvIDs(i);

                for cond = 1:3

                    psth = trialPSTHs{unitID}{cond};

                    %% Calculate ON-window DeltaFR
                    [deltaFR_on, peakStimFR_on, baselineMeanFR_on] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [-2 0], [1 2]);

                    %% Calculate OFF-window DeltaFR
                    [deltaFR_off, peakStimFR_off, baselineMeanFR_off] = ...
                        computeAbsPeakDeltaFR( ...
                            psth, PSTH_binEdges, [0 2], [-1 0]);

                    %% Keep the larger ON or OFF response
                    if deltaFR_on >= deltaFR_off
                        deltaFR = deltaFR_on;
                        peakStimFR = peakStimFR_on;
                        baselineMeanFR = baselineMeanFR_on;
                    else
                        deltaFR = deltaFR_off;
                        peakStimFR = peakStimFR_off;
                        baselineMeanFR = baselineMeanFR_off;
                    end

                    allRows(end+1, :) = { ...
                        stageName, groupName, 'unconventional', cond, unitID, ...
                        deltaFR, peakStimFR, baselineMeanFR, ...
                        psthFiles(f).name};

                    %% Calculate ON and OFF latency separately
                    [onsetLatency_on, peakLatency_on, peakFR_on, threshold_on] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [-2 0], [1 2], -2);

                    [onsetLatency_off, peakLatency_off, peakFR_off, threshold_off] = ...
                        computeResponseLatency( ...
                            psth, PSTH_binEdges, [0 2], [-1 0], 0);

                    %% Combine ON and OFF latency estimates
                    if ~isnan(onsetLatency_on) && ~isnan(onsetLatency_off)

                        onsetLatency = (onsetLatency_on + onsetLatency_off) / 2;
                        peakLatency  = (peakLatency_on  + peakLatency_off)  / 2;
                        peakFR       = (peakFR_on       + peakFR_off)       / 2;
                        threshold    = (threshold_on    + threshold_off)    / 2;

                    elseif ~isnan(onsetLatency_on) && isnan(onsetLatency_off)

                        onsetLatency = onsetLatency_on;
                        peakLatency  = peakLatency_on;
                        peakFR       = peakFR_on;
                        threshold    = threshold_on;

                    elseif isnan(onsetLatency_on) && ~isnan(onsetLatency_off)

                        onsetLatency = onsetLatency_off;
                        peakLatency  = peakLatency_off;
                        peakFR       = peakFR_off;
                        threshold    = threshold_off;

                    else

                        onsetLatency = NaN;
                        peakLatency  = NaN;
                        peakFR       = NaN;
                        threshold    = NaN;
                    end

                    allLatency(end+1, :) = { ...
                        stageName, groupName, 'unconventional', cond, unitID, ...
                        onsetLatency, peakLatency, peakFR, threshold, ...
                        psthFiles(f).name};
                end
            end
        end
    end
end

%% Convert firing-rate results into a table
T = cell2table(allRows, ...
    'VariableNames', {'Stage', 'Group', 'CellType', 'LightCondition', ...
    'UnitID', 'DeltaFR', 'PeakStimFR', 'BaselineMeanFR', 'FileID'});

% Convert grouping variables to categorical for easier plotting/statistics.
T.Stage = categorical(T.Stage);
T.Group = categorical(T.Group);
T.CellType = categorical(T.CellType);
T.LightCondition = categorical(T.LightCondition);

%% Convert latency results into a table
T1 = cell2table(allLatency, ...
    'VariableNames', {'Stage', 'Group', 'CellType', 'LightCondition', ...
    'UnitID', 'onsetLatency', 'peakLatency', 'peakFR', 'threshold', 'FileID'});

% Convert grouping variables to categorical for easier plotting/statistics.
T1.Stage = categorical(T1.Stage);
T1.Group = categorical(T1.Group);
T1.CellType = categorical(T1.CellType);
T1.LightCondition = categorical(T1.LightCondition);

end


function [deltaFR, peakStimFR, baselineMeanFR] = computeAbsPeakDeltaFR( ...
    psth, binEdges, responseWindow, baselineWindow)
% computeAbsPeakDeltaFR
%
% Calculates the absolute change in firing rate between the response window
% and the baseline window.
%
% INPUTS:
%   psth           : trials × time bins PSTH matrix
%   binEdges       : 1 × N+1 vector of PSTH bin edges
%   responseWindow : [start end] time window used to detect response peak
%   baselineWindow : [start end] time window used to calculate baseline FR
%
% OUTPUTS:
%   deltaFR        : abs(peakStimFR - baselineMeanFR)
%   peakStimFR     : maximum firing rate within the response window
%   baselineMeanFR : mean firing rate within the baseline window
%
% NOTE:
%   The PSTH values are converted to firing rate by dividing by bin width.

    %% Calculate bin width and bin centres
    binWidth = mean(diff(binEdges));
    binCenters = binEdges(1:end-1) + binWidth / 2;

    %% Average PSTH across trials and convert to firing rate
    meanPSTH = mean(psth, 1, 'omitnan');
    fr = meanPSTH ./ binWidth;

    %% Identify response and baseline bins
    responseBins = binCenters >= responseWindow(1) & binCenters < responseWindow(2);
    baselineBins = binCenters >= baselineWindow(1) & binCenters < baselineWindow(2);

    respFR = fr(responseBins);
    baselineFR = fr(baselineBins);

    %% Return NaN if either response or baseline window is empty
    if isempty(respFR) || all(isnan(respFR)) || ...
       isempty(baselineFR) || all(isnan(baselineFR))

        deltaFR = NaN;
        peakStimFR = NaN;
        baselineMeanFR = NaN;
        return;
    end

    %% Calculate peak response, baseline mean, and absolute DeltaFR
    peakStimFR = max(respFR, [], 'omitnan');
    baselineMeanFR = mean(baselineFR, 'omitnan');

    deltaFR = abs(peakStimFR - baselineMeanFR);

end


function [onsetLatency, peakLatency, peakFR, threshold] = computeResponseLatency( ...
    psth, binEdges, responseWindow, baselineWindow, eventTime)
% computeResponseLatency
%
% Estimates response onset latency and peak latency from a trial-averaged
% firing-rate trace.
%
% INPUTS:
%   psth           : trials × time bins PSTH matrix
%   binEdges       : 1 × N+1 vector of PSTH bin edges
%   responseWindow : [start end] response detection window
%                    Example:
%                       [-2 0] for ON response
%                       [0 2]  for OFF response
%
%   baselineWindow : [start end] baseline window
%                    Example:
%                       [1 2]  for ON response
%                       [-1 0] for OFF response
%
%   eventTime      : stimulus transition time used as latency zero
%                    Example:
%                       -2 for ON response
%                        0 for OFF response
%
% OUTPUTS:
%   onsetLatency : time from eventTime to the first sustained threshold
%                  crossing
%
%   peakLatency  : time from eventTime to the maximum firing rate within
%                  the response window
%
%   peakFR       : maximum firing rate within the response window
%
%   threshold    : baselineMean + 2 × baselineSD
%
% METHOD:
%   1. Average PSTH across trials.
%   2. Convert PSTH to firing rate using bin width.
%   3. Calculate baseline mean and SD.
%   4. Define threshold as baselineMean + 2 × baselineSD.
%   5. Detect onset as the first time point where at least two consecutive
%      bins exceed the threshold.
%   6. Detect peak latency as the time of maximum firing rate within the
%      response window.

    %% Calculate bin width and bin centres
    binWidth = mean(diff(binEdges));
    binCenters = binEdges(1:end-1) + binWidth / 2;

    %% Average PSTH across trials and convert to firing rate
    meanPSTH = mean(psth, 1, 'omitnan');
    fr = meanPSTH ./ binWidth;

    %% Identify response and baseline bins
    responseBins = binCenters >= responseWindow(1) & binCenters < responseWindow(2);
    baselineBins = binCenters >= baselineWindow(1) & binCenters < baselineWindow(2);

    baselineFR = fr(baselineBins);

    %% Return NaN if the baseline window is empty or invalid
    if isempty(baselineFR) || all(isnan(baselineFR))

        onsetLatency = NaN;
        peakLatency = NaN;
        peakFR = NaN;
        threshold = NaN;
        return;
    end

    %% Calculate response threshold from baseline activity
    baselineMean = mean(baselineFR, 'omitnan');
    baselineSD = std(baselineFR, 'omitnan');

    % Threshold can be adjusted.
    % 2 SD is more sensitive; 3 SD would be stricter.
    threshold = baselineMean + 2 * baselineSD;

    %% Extract firing rate during response window
    respFR = fr(responseBins);
    respTimes = binCenters(responseBins);

    %% Return NaN if response window is empty or invalid
    if isempty(respFR) || all(isnan(respFR))

        onsetLatency = NaN;
        peakLatency = NaN;
        peakFR = NaN;
        return;
    end

    %% Detect response onset
    % To reduce false detection from one noisy bin, onset must be defined
    % by at least two consecutive bins above threshold.

    aboveThreshold = respFR > threshold;
    minConsecutiveBins = 2;

    onsetIdx = NaN;

    for k = 1:(numel(aboveThreshold) - minConsecutiveBins + 1)
        if all(aboveThreshold(k:k + minConsecutiveBins - 1))
            onsetIdx = k;
            break;
        end
    end

    %% Convert onset index to latency relative to event time
    if isnan(onsetIdx)
        onsetLatency = NaN;
    else
        onsetLatency = respTimes(onsetIdx) - eventTime;
    end

    %% Calculate peak latency
    [peakFR, peakIdx] = max(respFR, [], 'omitnan');
    peakLatency = respTimes(peakIdx) - eventTime;

end