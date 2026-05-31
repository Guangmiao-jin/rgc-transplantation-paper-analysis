function [postMean, postOverTau] = safe_post_stats(series, tauBin, peakVal, epsAbs, minPostBins)
% finds the decay rate of 
% postmean is the mean rate of firing after the Ti until stim ends. 
% postOveTau - postmean rate / peak for reasons?????? check 
% ratio of post-tau mean firing rate/ tau time firing ( high ratio is fast
% decay, low ratio is slow decay. 

series = series(:)';
N = numel(series);
if N==0 || ~isfinite(tauBin) || tauBin<1 || tauBin>N
    postMean    = NaN;
    postOverTau = NaN;
    return;
end

% post-τ region
startIdx = min(tauBin + 1, N);
if minPostBins > 0
    stopIdx = min(startIdx + minPostBins - 1, N);
else
    stopIdx = N;
end

seg = series(startIdx:stopIdx);
seg = seg(isfinite(seg));
if isempty(seg)
    postMean = 0;             %
else
    postMean = mean(seg, 'omitnan');
end

%
denomTau = series(tauBin);
if isfinite(denomTau) && denomTau > epsAbs
    denom = denomTau;
elseif isfinite(peakVal) && peakVal > epsAbs
    denom = peakVal;
else
    denom = epsAbs;
end

if all(seg==0) && (~isfinite(denomTau) || denomTau<=epsAbs)
    postOverTau = 0;
else
    postOverTau = postMean / denom;
end
end
