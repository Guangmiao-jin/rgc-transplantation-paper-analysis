function data = processRetinaFlashStimData(clusterFilepath)
% Processing the sorted cluster file output from Herdingspikes. Creates
% flash stimulus aligned PTSHs and raster plots. Classifies neurons into
% ON, OFF, ON_OFF, and unconvential
%
% Inputs: clusterFilepath- fullfile path to *_cluster.hdf5

%% defaults
prestimTimeZScore = 1; % sec
stimTimeZScore = 1.5; % sec

prestimTimePSTH = 2; % sec
postStimTimePSTH = 2; % sec

responseQualityThreshold = 0.15;

cellRasterFolder = extractBefore(clusterFilepath, '.');
binsize = 25; % ms

%% get the appropriate paths for stim on and off triggers
filepathPrefix = extractBefore(clusterFilepath, '_cluster');
stimOnFile = [filepathPrefix{:} '_triggerON.npy'];
stimOffFile = [filepathPrefix{:} '_triggerOFF.npy'];

%% load data
matSavePath = extractBefore(clusterFilepath, '.');
matSavePath = [matSavePath{:} '.mat'];

if ~exist(matSavePath)
    data = readHS2_FLAME(clusterFilepath);
    save(matSavePath, "data", "-v7.3");
else
    load(matSavePath);
end

stimOnFrames = double(readNPY(stimOnFile));
stimOffFrames = double(readNPY(stimOffFile));

%% split on frames into blocks
stimOnPerBlock = splitStimEvents2Blocks(stimOnFrames,data);

%% split off frames into blocks
stimOffPerBlock = splitStimEvents2Blocks(stimOffFrames,data);

%% make all the metrics we use to seperate out the cells
[pltcurve, responseMetrics] = createFlashMetrics(data, stimOnPerBlock, stimOffPerBlock, prestimTimeZScore, stimTimeZScore, prestimTimePSTH, postStimTimePSTH);

%% check cluster waveforms to make sure they are valid
[validIds] = validateClusterWaveform(data, responseMetrics.ISI_vio);

%% start the plotting and filter cluster IDs into cell types

% set up blank struct
totalneurons.num_OnNeurons.trans = [];
totalneurons.num_OnNeurons.sus = [];

totalneurons.num_OffNeurons.trans = [];
totalneurons.num_OffNeurons.sus = [];

totalneurons.num_OnOffNeurons = [];
totalneurons.num_unconventional = [];

totalneurons.num_notRespond = [];

count = 0;
for i = 1:length(validIds) % go through all valid neurons
    count = count +1;

    curCl = validIds(i);
    clusterID = curCl-1; % for file naming, HS2 is 0 indexe

    disp(['On ' num2str(count) ' of ' num2str(length(validIds))]);

    % put all cluster responses into structure
    clusterResponses.responseQuality = responseMetrics.responseQuality(curCl,:);
    clusterResponses.trialSpikes = responseMetrics.trialSpikes{curCl};
    clusterResponses.trialPSTHs = responseMetrics.trialPSTHs{curCl};
    clusterResponses.binEdges = responseMetrics.PSTH_binEdges;
    clusterResponses.BI = responseMetrics.BI{curCl};
    clusterResponses.TI_on =responseMetrics.TI_on{curCl};
    clusterResponses.TI_off = responseMetrics.TI_off{curCl};
    clusterResponses.ratio =responseMetrics.ratio{curCl};
    clusterResponses.onslope = responseMetrics.onslope{curCl};
    clusterResponses.offslope= responseMetrics.offslope{curCl};

%% giant filter to classify the neurons
    if max(clusterResponses.responseQuality)>responseQualityThreshold
        if mean(clusterResponses.BI,'omitnan') >= 1/3 % ON
            if mean(clusterResponses.TI_on,'omitnan') <= 0.1 && mean(clusterResponses.onslope,'omitnan') <= 0.45
                totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                cls = 'OnNeurons/trans';
            elseif mean(clusterResponses.TI_on,'omitnan') <= 0.1 && mean(clusterResponses.onslope,'omitnan') > 0.45
                totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                cls = 'OnNeurons/sus';
            elseif mean(clusterResponses.TI_on,'omitnan') > 0.1 && mean(clusterResponses.onslope,'omitnan') <= 0.45
                totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                cls = 'OnNeurons/trans';
            elseif mean(clusterResponses.TI_on,'omitnan') > 0.1 && mean(clusterResponses.onslope,'omitnan') > 0.45
                totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                cls = 'OnNeurons/sus';
            end
        elseif mean(clusterResponses.BI,'omitnan') <= -1/3 % OFF
            if mean(clusterResponses.TI_off,'omitnan') <= 0.1 && mean(clusterResponses.offslope,'omitnan') <= 0.45
                totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                cls = 'OffNeurons/trans';
            elseif mean(clusterResponses.TI_off,'omitnan') <= 0.1 && mean(clusterResponses.offslope,'omitnan') > 0.45
                totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                cls = 'OffNeurons/sus';
            elseif mean(clusterResponses.TI_off,'omitnan') > 0.1 && mean(clusterResponses.offslope,'omitnan') <= 0.45
                totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                cls = 'OffNeurons/trans';
            elseif mean(clusterResponses.TI_off,'omitnan') > 0.1 && mean(clusterResponses.offslope,'omitnan') > 0.45
                totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                cls = 'OffNeurons/sus';
            end
        else
            if mean(clusterResponses.ratio,'omitnan') >= 5.5 % ON based on ratio
                if mean(clusterResponses.TI_on,'omitnan') <= 0.1 && mean(clusterResponses.onslope,'omitnan') <= 0.45
                    totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/trans';
                elseif mean(clusterResponses.TI_on,'omitnan') <= 0.1 && mean(clusterResponses.onslope,'omitnan') > 0.45
                    totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/sus';
                elseif mean(clusterResponses.TI_on,'omitnan') > 0.1 && mean(clusterResponses.onslope,'omitnan') <= 0.45
                    totalneurons.num_OnNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/trans';
                elseif mean(clusterResponses.TI_on,'omitnan') > 0.1 && mean(clusterResponses.onslope,'omitnan') > 0.45
                    totalneurons.num_OnNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OnNeurons/sus';
                end
            elseif mean(clusterResponses.ratio,'omitnan') <= 0.15 % Off based on ratio
                if mean(clusterResponses.TI_off,'omitnan') <= 0.1 && mean(clusterResponses.offslope,'omitnan') <= 0.45
                    totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/trans';
                elseif mean(clusterResponses.TI_off,'omitnan') <= 0.1 && mean(clusterResponses.offslope,'omitnan') > 0.45
                    totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/sus';
                elseif mean(clusterResponses.TI_off,'omitnan') > 0.1 && mean(clusterResponses.offslope,'omitnan') <= 0.45
                    totalneurons.num_OffNeurons.trans(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/trans';
                elseif mean(clusterResponses.TI_off,'omitnan') > 0.1 && mean(clusterResponses.offslope,'omitnan') > 0.45
                    totalneurons.num_OffNeurons.sus(end+1) = data.channelNames{4,curCl};
                    cls = 'OffNeurons/sus';
                end
            else % anything inbetween BI -1/3 and + 1/3, ratio 0.15 and 5.5 is ON_OFF (ratio is obervaation based)
                totalneurons.num_OnOffNeurons(end+1) = data.channelNames{4,curCl};
                cls = 'OnOffNeurons';
            end
        end
    elseif max(clusterResponses.responseQuality)>0.05 && max(clusterResponses.responseQuality) <= responseQualityThreshold % catch for unconventional
        totalneurons.num_unconventional(end+1) = data.channelNames{4,curCl};
        cls = 'unconventional';
    else % catch crap
        totalneurons.num_notRespond(end+1) = data.channelNames{4,curCl};
        cls = 'notRespond';
    end

    totalneurons.allneurons = [totalneurons.num_OnNeurons.trans, totalneurons.num_OnNeurons.sus, ...
        totalneurons.num_OffNeurons.trans, totalneurons.num_OffNeurons.sus,...
        totalneurons.num_OnOffNeurons, totalneurons.num_unconventional, totalneurons.num_notRespond];
    %% start actually plotting 

    outputFolder = fullfile([cellRasterFolder{:} '_PSTHPlotsV2'], cls);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    [ax, tlo] = plotAllRasterPSTHs_On_Off_Response(clusterResponses,binsize);
    sgtitle(tlo, ['Cluster ID: ' num2str(clusterID) ' Spks: ' num2str(data.channelNames{6,curCl})], 'FontWeight','bold');

    saveName = sprintf('%s\\cluster%04d.png',outputFolder,clusterID);

    warning('off','all')
    winHandle = gethwnd(gcf);
    cmndstr = sprintf('%s','MiniCap.exe -save ','"',saveName,'"',...
        ' -compress 9', ' -capturehwnd ', num2str(winHandle),' -exit');
    system(cmndstr);
    close
    warning('on','all')
end
%% save the structures
save(strcat(cellRasterFolder,"_responseMetrics.mat"),'responseMetrics','-v7.3');
save(strcat(cellRasterFolder,"_totalneuronsV2.mat"), 'totalneurons');  %
save(strcat(cellRasterFolder,"_psth.mat"), 'pltcurve','-v7.3');  %
end
