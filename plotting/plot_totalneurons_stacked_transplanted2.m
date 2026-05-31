function outTbl = plot_totalneurons_stacked_transplanted2(baseFolder)
% baseFolder: root directory of rd10_transplanted
%
% Directory structure:
% rd10_transplanted/
% early degeneration/
% NRL/behav pos or behav neg
% CRX/behav pos or behav neg
% SHAM/
% late degeneration/
% NRL/behav pos or behav neg
% CRX/behav pos or behav neg
% SHAM/
%
% Output: outTbl (MATLAB table)
% -----------------------------
% Define 6  groups + 1 spacer to create a visual gap in the plot
% -----------------------------
stages = {'early degeneration', 'late degeneration'};
behavPos1 = 'behavioural positive and ephys positive';
behavPos2 = 'behavioural positive and ephys negative';
behavNeg = 'behavioural negative';

mainCategories = {'OnNeurons', 'OffNeurons', 'OnOffNeurons', 'unconventional', 'notRespond'};

% Desired plotting order: two NRL bars, gap, two CRX bars, gap, one SHAM bar
plotGroupNames = {'NRL behav+ ephys+','NRL behav+ ephys+',....
    'NRL behav-','Spacer1',...
    'CRX behav+ ephys+','CRX behav+ ephys+',....
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

            mainData = zeros(1, numel(mainCategories));
            for k = 1:numel(mainCategories)
                fieldName = ['num_' mainCategories{k}];
                if isfield(matData1, fieldName)
                    if isstruct(matData1.(fieldName)) && ismember(mainCategories{k}, {'OnNeurons','OffNeurons'})
                        mainData(k) = mainData(k) + length(struct2array(matData1.(fieldName)));
                    else
                        mainData(k) = mainData(k) + length(matData1.(fieldName));
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