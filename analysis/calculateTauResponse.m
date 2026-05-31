function  [responseTau, meanPostA2] = calculateTauResponse(meanStimAlignedTr, binsize)
% Create time to response decay for each light condition. It is eseentially
% the time taken to go from peak response to expontential decay of that
% firing rate by finiding the time for peak response/e. 
% See 
% Classification of Retinal Ganglion Cells: A Statistical Approach
% Stephen M. Carcieri, Adam L. Jacobs, and Sheila Nirenberg
% Journal of Neurophysiology 2003 90:3, 1704-1713
%
% Input: meanStimAlignedTr- cell array of trial averaged histogram
%                           response for each light level used
%
%        binsize - size of histogram bins in ms
%
% Output: responseTau - vector of response decay times

%% defaults

EPS_ABS        = 1e-9;    % gets around divide by zero errors                     
MIN_POST_BINS  = 0;        

%% response decay tau

% for each stim condition
for i =1:length(meanStimAlignedTr)
    temp = meanStimAlignedTr{i};
    [peakResponse, peakInxd] = max(temp);

    tauDecayLevel = peakResponse / exp(1);
    postPeakTrace = temp(peakInxd:end);

    thresholdBin  = find(postPeakTrace < tauDecayLevel, 1, 'first');

    if ~isempty(thresholdBin)
        tauBin_on = peakInxd + thresholdBin - 1; % time bin at which A2 level is achieved
    else
        tauBin_on = numel(temp);
    end

    responseTau(i) = (tauBin_on - peakInxd) * (binsize/1000);  % s % time at which A2 level is achieved


%% CHECK THIS WITH GJ!!!!!!!!!!!!!!!!!!!!!!!!
    % Post-τ average firing rate
    [meanPostA2_on(i), meanPostA2(i)] = safe_post_stats( ...
        temp, tauBin_on, peakResponse, EPS_ABS, MIN_POST_BINS);

end
end