function [trialHists, trSpikesStructON,  edges] = createTrialPSTHs(spikeFrames, fs, binsize, stimON_Events, stimOFF_Events, preAlignTime, postAlignTime)
% Creates trial by trial PSTHs and aligned spike structure
%
% Inputs: spikeFrames- vector of spike frames for cluster in question
%         fs - recording sampling rate
%         binsize - bin size for histogram binning
%         stimON_Events - cell array with ON events split into blocks,
%                         usually 4 x 1
%         stimOFF_Events - cell array with OFF events split into blocks,
%                         usually 4 x 1
%         preAlignTime - prestimulus aligment time for trials (usually
%                        2 sec)
%         postAlignTime - poststimulus aligment time for trials (usually
%                        2 sec)
%
% Outputs: trialHists -  trial by trial histogram with bin size based on
%                        'binsize' input
%          trSpikesStructON - trial by trial spike structure zeroed to
%                             stimON with prestim and poststim time
%          edges - bin edges for plotting the histograms in trialHists

%%
if isscalar(stimON_Events) % if only one condition (LEGACY)
    stimblk = 1;
    preAlignInFrames = preAlignTime * fs;
    postAlignInFrames = postAlignTime * fs;

    % for each trial
    for tr = 1:length(stimON_Events{stimblk})

        % trial times in frames
        trStart = stimON_Events{stimblk}(tr)- preAlignInFrames;
        trEnd = stimON_Events{stimblk}(tr)+ postAlignInFrames;

        % find inclusions for trial start and end
        trStartIndx = find(spikeFrames >trStart,1, 'first');
        trEndIndx = find(spikeFrames <trEnd,1, 'last');

        % get the spikes
        trSpikes = spikeFrames(trStartIndx:  trEndIndx);
        trSpikes = trSpikes / fs; % convert to sec
        trSpikes = trSpikes- (stimON_Events{stimblk}(tr)/ fs); % rezero to alignment event

        % put into cell array for psth
        trSpikesStructON{stimblk}{tr} = trSpikes;
        trialLenSec(tr) = (trEnd-trStart)/fs;

        % wipe variables for next trial
        trSpikes = [];
    end

    %% PSTH
    % trialLen            = mean(trialLenSec) * 1000;                % trial length ms
    binActual =[0 :(binsize/1000): trialLenSec]- preAlignTime;
    % nbins               = round(trialLen/binsize);                        % Bin duration in [ms]
    % nobins              = 1000/binsize;                            % No of bins/sec

    for iTrial = 1:length(trSpikesStructON{stimblk})
        [trialHistTemp, edges] = histcounts(trSpikesStructON{stimblk}{iTrial},binActual);
        trialHists{stimblk}(iTrial,:) = trialHistTemp;
    end

else
    % for stimblk = 1:length(stimON_Events)-1
    for stimblk = 1:3 % HACK 20250623 MS

        % alignment times in sec
        preAlignInFrames = preAlignTime * fs;
        postAlignInFrames = postAlignTime * fs;

        % for each trial
        for tr = 1:length(stimON_Events{stimblk})

            % trial times in frames
            trStart = stimON_Events{stimblk}(tr)- preAlignInFrames;
            trEnd = stimON_Events{stimblk}(tr)+ postAlignInFrames;

            % find inclusions for trial start and end
            trStartIndx = find(spikeFrames >trStart,1, 'first');
            trEndIndx = find(spikeFrames <trEnd,1, 'last');

            % get the spikes
            trSpikes = spikeFrames(trStartIndx:  trEndIndx);
            trSpikes = trSpikes / fs; % convert to sec
            trSpikes = trSpikes- (stimON_Events{stimblk}(tr)/ fs); % rezero to alignment event

            % put into cell array for psth
            trSpikesStructON{stimblk}{tr} = trSpikes;
            trialLenSec(tr) = (trEnd-trStart)/fs;

            % wipe variables for next trial
            trSpikes = [];
        end

        %% PSTH
        trialLen            = mean(trialLenSec) * 1000;                % trial length ms
        binActual = [0 :(binsize/1000): trialLenSec]- preAlignTime;
        nbins               = round(trialLen/binsize);                        % Bin duration in [ms]
        nobins              = 1000/binsize;                            % No of bins/sec

        for iTrial = 1:length(trSpikesStructON{stimblk})
            [trialHistTemp, edges] = histcounts(trSpikesStructON{stimblk}{iTrial},binActual);
            trialHists{stimblk}(iTrial,:) = trialHistTemp;
        end
    end
end
end