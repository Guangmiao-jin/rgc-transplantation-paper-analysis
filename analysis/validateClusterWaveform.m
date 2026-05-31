function [validIds] = validateClusterWaveform(data, ISI_vio)
% Check if waveforms are valid by checking whether any postive peak is
% larger than the negative deflection associated with spikes. Remove any
% clusters with ISI violation larger that 5%
%
% Inputs: data- data from herdingspikes
%         ISI_vio- inter spike interval violation rates for each cluster
% Outputs: validIds- Index numbers of valid clusters

validIds = [];
for idx = 1:length(data.centres)
    meanwaveforms = data.waveformClusterMeans(idx,:);
    ISI_vio_temp = ISI_vio(idx);

    % detect negative peaks and positive peaks
    [invertedPeaks, t0] = findpeaks(-meanwaveforms);
    numNegativeVally = sum(invertedPeaks>0);
    [~,t1] = max(meanwaveforms);

    % save those neuron ids which fall into given conditions
    if data.channelNames{6,idx} >= 250 % more than 250 spikes total

        if numNegativeVally == 1 && abs(invertedPeaks(1)) > 1000 % make sure there is only one negative peak and absolute value is greater than 1000
            peaktimediff = t1-t0(1);

            if peaktimediff > 0 && ISI_vio_temp < 0.05
                validIds = [validIds; idx];
            end
        end
    end
end

end