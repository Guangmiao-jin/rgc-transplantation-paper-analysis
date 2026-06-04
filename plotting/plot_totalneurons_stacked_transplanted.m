function outTbl = plot_totalneurons_stacked_transplanted(baseFolder)
%
% Summarises and plots the distribution of classified RGC response types
% across transplanted and sham rd10 retinal MEA recordings.
%
% This function recursively searches for *_totalneuronsV2.mat files under
% each experimental subgroup, extracts the number of neurons classified as:
%   - ON
%   - OFF
%   - ON-OFF
%   - unconventional
%   - non-responsive
%
% It then pools counts across recordings within each subgroup, calculates
% within-group percentages, generates one stacked bar plot per degeneration
% stage, and exports a summary table to Excel.
%
% INPUT:
%   baseFolder : root directory of the rd10 transplantation analysis folder.
%
%                Expected folder structure:
%
%                baseFolder/
%                    early degeneration/
%                        NRL/
%                            behavioural positive and ephys positive/
%                            behavioural positive and ephys negative/
%                            behavioural negative/
%                        CRX/
%                            behavioural positive and ephys positive/
%                            behavioural positive and ephys negative/
%                            behavioural negative/
%                        SHAM/
%
%                    late degeneration/
%                        NRL/
%                            behavioural positive and ephys positive/
%                            behavioural positive and ephys negative/
%                            behavioural negative/
%                        CRX/
%                            behavioural positive and ephys positive/
%                            behavioural positive and ephys negative/
%                            behavioural negative/
%                        SHAM/
%
% OUTPUT:
%   outTbl : summary table containing, for each stage and group:
%              - number of files analysed
%              - total number of classified cells
%              - raw counts of each response category
%              - percentage of each response category
%
%   The function also saves:
%       1. One stacked bar plot per stage:
%              Neuron_types_early_degeneration.png
%              Neuron_types_late_degeneration.png
%
%       2. One Excel summary table:
%              summary_rd10_transplanted_early_late.xlsx
%
% RESPONSE CATEGORIES:
%   OnNeurons      = ON transient + ON sustained
%   OffNeurons     = OFF transient + OFF sustained
%   OnOffNeurons   = ON-OFF responsive cells
%   unconventional = non-canonical responsive cells
%   notRespond     = non-responsive cells
%
% NOTE:
%   Percentages are calculated within each subgroup:
%
%       category percentage =
%           category count / total classified cells in that subgroup * 100
%
%   Spacer groups are used only to create visual gaps in the plot and are
%   excluded from plotting/statistical interpretation.
% -----------------------------
stages = {'early degeneration', 'late degeneration'};
behavPos1 = 'behavioural positive and ephys positive';
behavPos2 = 'behavioural positive and ephys negative';
behavNeg = 'behavioural negative';

mainCategories = {'OnNeurons', 'OffNeurons', 'OnOffNeurons', 'unconventional', 'notRespond'};

% Desired plotting order: two NRL bars, gap, two CRX bars, gap, one SHAM bar
plotGroupNames = {'NRL behav+ ephys+','NRL behav+ ephys-',....
    'NRL behav-','Spacer1',...
    'CRX behav+ ephys+','CRX behav+ ephys-',....
    'CRX behav-','Spacer2','SHAM'};

% Corresponding folders to scan under each stage
scanSpec = { ...
    {behavPos1,'NRL'}, ...
    {behavPos2,'NRL'}, ...
    {behavNeg,'NRL'}, ...
    {}, ...
    {behavPos1,'CRX'}, ...
    {behavPos2,'CRX'}, ...
    {behavNeg,'CRX'}, ...
    {}, ...
    {'SHAM'} ...
};

% -----------------------------
% Initialise output table container
% -----------------------------
outRows = [];

% -----------------------------
% Process each stage separately and generate one plot per stage
% -----------------------------
for s = 1:numel(stages)
    stageName = stages{s};
    stagePath = fullfile(baseFolder, stageName);

    nG = numel(plotGroupNames);
    data = struct();

    for i = 1:nG
        data(i).group = plotGroupNames{i};
        data(i).mainValues = zeros(1, numel(mainCategories));
        data(i).nFiles = 0;
        data(i).stage = stageName;
    end

    for i = 1:nG
        if isempty(scanSpec{i})
            continue; % spacer
        end

        spec = scanSpec{i};

        % 1) SHAM: spec = {'SHAM'}
        % 2) Other groups: spec = {behavioural status, treatment}
        if isscalar(spec)
            groupPath = fullfile(stagePath, spec{1});
        else
            behav = spec{1};
            treat = spec{2};
            groupPath = fullfile(stagePath, treat, behav);
        end

        if ~isfolder(groupPath)
            continue; % leave this group blank if the folder is missing
        end

        % Use recursive search ('**') in case files are stored in subfolders
        matFiles = dir(fullfile(groupPath, '**', '*_totalneuronsV2.mat'));
        if isempty(matFiles)
            continue;
        end

        for j = 1:length(matFiles)
            matPath = fullfile(matFiles(j).folder, matFiles(j).name);
            matData = load(matPath);

            if ~isfield(matData, 'totalneurons')
                continue;
            end
            matData1 = matData.totalneurons;

%% Extract response-category counts from this totalneuronsV2 file
% ON and OFF categories are stored as structures containing transient and
% sustained subtypes, so these are combined explicitly.
%
% Other categories are stored as vectors of neuron IDs.
%
% unique() is used for ON/OFF to avoid double-counting if any neuron appears
% in both transient and sustained fields.

mainData = zeros(1, numel(mainCategories));

for k = 1:numel(mainCategories)

    switch mainCategories{k}

        case 'OnNeurons'
            if isfield(matData1, 'num_OnNeurons')
                ids = [
                    matData1.num_OnNeurons.trans(:);
                    matData1.num_OnNeurons.sus(:)
                ];
                mainData(k) = numel(unique(ids));
            end

        case 'OffNeurons'
            if isfield(matData1, 'num_OffNeurons')
                ids = [
                    matData1.num_OffNeurons.trans(:);
                    matData1.num_OffNeurons.sus(:)
                ];
                mainData(k) = numel(unique(ids));
            end

        case 'OnOffNeurons'
            if isfield(matData1, 'num_OnOffNeurons')
                mainData(k) = numel(matData1.num_OnOffNeurons(:));
            end

        case 'unconventional'
            if isfield(matData1, 'num_unconventional')
                mainData(k) = numel(matData1.num_unconventional(:));
            end

        case 'notRespond'
            if isfield(matData1, 'num_notRespond')
                mainData(k) = numel(matData1.num_notRespond(:));
            end
    end
end

            data(i).mainValues = data(i).mainValues + mainData;
            data(i).nFiles = data(i).nFiles + 1;
        end
    end

    % -----------------------------
    % Generate output table rows for this stage
    % -----------------------------
    groups = string({data.group})';
    isSpacer = startsWith(groups, "Spacer");

    rawMain = vertcat(data.mainValues);
    nFiles  = [data.nFiles]';

    totalCells = sum(rawMain, 2);

    % Calculate within-group percentages
    mainPct = nan(size(rawMain));
    for i = 1:nG
        if isSpacer(i), continue; end
        if totalCells(i) > 0
            mainPct(i,:) = (rawMain(i,:) ./ totalCells(i)) * 100;
        end
    end

    stageCol = repmat(string(stageName), nG, 1);

    stageTbl = table(stageCol, groups, nFiles, totalCells, ...
        rawMain(:,1), rawMain(:,2), rawMain(:,3), rawMain(:,4), rawMain(:,5), ...
        mainPct(:,1), mainPct(:,2), mainPct(:,3), mainPct(:,4), mainPct(:,5), ...
        'VariableNames', { ...
            'Stage','Group','nFiles','TotalCells', ...
            'OnNeurons','OffNeurons','OnOffNeurons','unconventional','notRespond', ...
            'Pct_On','Pct_Off','Pct_OnOff','Pct_unconventional','Pct_notRespond' ...
        });

    outRows = [outRows; stageTbl]; 

    % -----------------------------
    % Plot stacked bars for non-spacer groups while keeping spacer positions as gaps
    % -----------------------------
    % Use numeric x positions to ensure stable spacing without relying on categorical axes
    xPos = [1 2 3 4 5 6 7 8 9]; % Spacer1 = 4, Spacer2 = 8
    idxPlot = ~isSpacer;    % groups to plot as bars

    % Plot bars only for the non-spacer groups
    xBar = xPos(idxPlot);
    YBar = mainPct(idxPlot, :)';  % transpose to match the stacked bar input format

    figure('Color','w','Position',[100 100 950 600]); hold on;

    h = bar(xBar, YBar', 'stacked', 'BarWidth', 0.7);

    title(['Neuron Type Distribution - ' stageName], 'Interpreter','none');
    ylabel('Percentage of Neurons (%)');
    set(gca,'XTick', xPos);
    set(gca,'XTickLabel', {'NRL-behav+ ephys+','NRL-behav+ ephys-','NRL-behav-', ' '....
        'CRX-behav+ ephys+','CRX-behav+ ephys-','CRX-behav-',' ','SHAM'});
    xlim([0.3 9.7]);
    grid on; box off;

    legend(mainCategories, 'Location','bestoutside', 'Interpreter','none');

    colors1 = lines(numel(mainCategories));
    for k = 1:numel(h)
        h(k).FaceColor = colors1(k,:);
    end

    % Add nFiles / TotalCells labels above each bar.
    % This is important to avoid over-interpreting percentages from small sample sizes.
    barGroups = find(idxPlot);
    for bi = 1:numel(barGroups)
        gi = barGroups(bi); % original group index
        if totalCells(gi) <= 0, continue; end
        txt = sprintf('nFiles=%d\\nCells=%d', nFiles(gi), totalCells(gi));
        text(xPos(gi), 102, txt, 'HorizontalAlignment','center', 'FontSize', 9);
    end
    ylim([0 115]);

    % Save the plot
    outPng = fullfile(baseFolder, ['Neuron_types_' strrep(stageName,' ','_') '.png']);
    saveas(gcf, outPng);
    close(gcf);
end

% Combine all stage-specific output tables
outTbl = outRows;

% Write the summary table to Excel
writetable(outTbl, fullfile(baseFolder, 'summary_rd10_transplanted_early_late.xlsx'));

end