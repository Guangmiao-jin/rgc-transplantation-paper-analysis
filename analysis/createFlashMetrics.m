function [pltcurve, responseMetrics] = createFlashMetrics(data, stimON_Events, stimOFF_Events, ~, ~, prestimTimePSTH, postStimTimePSTH)
% function [pltcurve, responseMetrics, responseSum] = createFlashMetrics(data, stimON_Events, stimOFF_Events, ~, ~, prestimTimePSTH, postStimTimePSTH)


% preAlignTime = 1; %s
% postAlignTime = 1; %s
binsize = 25; % ms
fs = data.Sampling;
data.recLen = max(data.times)/fs;

%% Run through all the clusters
% for each cluster
spikes = data.spiketimestamps;
 % parfor i = 1:length(data.spiketimestamps)
 for i = 1:length(data.spiketimestamps)

    spikeFrames = spikes{i};
    % out = random_shift(spikeFrames,stimON_Events,preAlignTime, posAlignTime,binsize,fs, 1000);
    %% ISI violation value
    ISI_vio(i)= calculateISI_violations(spikeFrames, fs, data.recLen);
    %% create trial based PSTHs
    warning('off','MATLAB:colon:operandsNotRealScalar'); % stops Warning: Colon operands must be real scalars. This warning will become an error in a future release.
    [trialPSTHs{i}, trialSpikeStruct{i}, binEdges{i}] = createTrialPSTHs(spikeFrames, fs, binsize, stimON_Events, stimOFF_Events, prestimTimePSTH, postStimTimePSTH);
    [trialPSTHs50] = createTrialPSTHs(spikeFrames, fs, 50, stimON_Events, stimOFF_Events, prestimTimePSTH, postStimTimePSTH); % need 50ms bin PSTHs for quality index
    
    % create ON and OFF aligned PSTHs/histcounts
    [PSTHsON_Aligned, trSpikesStructON_Aligned, meanStimAlignedTrON, edgesON_Aligned] = createEventAlignedPSTHs(spikeFrames, fs, binsize, stimON_Events,  2);
    [PSTHsOFF_Aligned, trSpikesStructOFF_Aligned, meanStimAlignedTrOFF, edgesON_Aligned] = createEventAlignedPSTHs(spikeFrames, fs, binsize, stimOFF_Events,  2);
    
    warning('on','MATLAB:colon:operandsNotRealScalar');

    %% get response quality
    QI(i,:) = retinaResponseQuality(trialPSTHs50);

    %% ISI coefficient of variation for ON and OFF period
    ISIon_cv{i} = calculateISI_CV(trSpikesStructON_Aligned);
    ISIoff_cv{i} = calculateISI_CV(trSpikesStructOFF_Aligned);

    %% max responses for ON/ OFF
    peakHeight_on{i} = cellfun(@max, meanStimAlignedTrON);
    peakHeight_off{i} = cellfun(@max, meanStimAlignedTrOFF);

    %% ratio of mean response for ON vs OFF
    ratio{i} = cellfun(@mean, meanStimAlignedTrON)./cellfun(@mean, meanStimAlignedTrOFF);

    %% calculate bias index
    BI{i} = calculateBiasIndex(peakHeight_on{i}, peakHeight_off{i});

    %% tau response metrics (peak response / e ) time
    [Ti_on{i},ratioPost_on{i} ] = calculateTauResponse(meanStimAlignedTrON, binsize);
    [Ti_off{i}, ratioPost_off{i}] = calculateTauResponse(meanStimAlignedTrOFF, binsize);

 end

%% Put everything into struct
responseMetrics.ISI_vio = ISI_vio;
responseMetrics.ISIon_cv = ISIon_cv;
responseMetrics.ISIoff_cv = ISIoff_cv;
responseMetrics.trialSpikes = trialSpikeStruct;
responseMetrics.trialPSTHs = trialPSTHs;
responseMetrics.PSTH_binEdges = binEdges{1};
responseMetrics.responseQuality = QI;
responseMetrics.BI = BI;
responseMetrics.TI_on = Ti_on;
responseMetrics.TI_off = Ti_off;
responseMetrics.ratio = ratio;
responseMetrics.onslope = ratioPost_on;
responseMetrics.offslope =  ratioPost_off;

pltcurve.trialPSTHs = trialPSTHs;
pltcurve.PSTH_binEdges = binEdges{1};
end