function ISI_vio = calculateISI_violations(spikeFrames, fs, reclen)
% Calculates ISI violation rate based on 
%
% Quality Metrics to Accompany Spike Sorting of Extracellular Signals
% Daniel N. Hill, Samar B. Mehta and David Kleinfeld
% https://doi.org/10.1523/JNEUROSCI.0971-11.2011
%
% See https://spikeinterface.readthedocs.io/en/stable/modules/qualitymetrics/isi_violations.html
%
% Inputs: spikeFrames- vector of spike frames
%         fs - sampling rate
%         reclen - recording length in seconds
%
% Outputs: ISI_vio - ISI violation metric

%%
ISI_threshold = 0.001; % 1 ms ISI threshold

ISI = diff(spikeFrames/fs);

num_vio = sum(ISI < ISI_threshold);
ISI_min = 1/fs;
ISI_vio = abs((num_vio * reclen) / (2*length(spikeFrames)^2*(ISI_threshold-ISI_min)));
end