function stimPerBlock = splitStimEvents2Blocks(stimFrames,data)
% Splits the stimFrames up into blocks based on the time limit between
% blocks
%
% Inputs: stimFrames - vector of frames of the event in question
%
%         data - data structure containing the extracted data from HS2
%
% Output: outBlocks - cell array with the events split into blocks, usually
%                     4 x 1


%% split on frames into blocks
blockLimit = 10 * data.Sampling; % 10s

diffOn = diff(stimFrames);

% first stim on
blockStarts = 1;

% block stim on starts
stimBreaks = [blockStarts; find(diffOn > blockLimit)+1];

% block stim on ends
stimStopBreaks = [ find(diffOn > blockLimit) ;length(diffOn)+1];


for i =1:length(stimBreaks)
    stimPerBlock{i,:} = stimFrames(stimBreaks(i):stimStopBreaks(i));
end

% % % plotting for testing
% blockStarts = [stimFrames(stimBreaks)];
% blockEnds = [stimFrames(stimStopBreaks)];
% scatter(stimFrames, ones(length(stimFrames))*100);
% hold on
% scatter(blockStarts, repmat(100.5, 1, length(blockStarts)));
% scatter(blockEnds, repmat(100.5, 1, length(blockEnds)));

end