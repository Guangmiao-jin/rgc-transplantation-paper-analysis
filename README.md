# rgc-transplantation-paper-analysis

This repository contains HerdingSpikes2-based Python pre-processing scripts and MATLAB analysis scripts for RGC light-response classification and transplantation-related electrophysiology analyses.

## Contributors

This analysis pipeline was developed by Guangmiao Jin, with contributions from Dr Michael Savage.

- **Guangmiao Jin** developed the overall pipeline structure, selected and adjusted key analysis and classification parameters, implemented the RGC response classification workflow, and organised the plotting comparison analyses between transplanted group and SHAM group.
- **Dr Michael Savage** contributed to the Python stimulation pulse extraction script, MATLAB code optimisation, code annotation and commenting, and technical review of the analysis pipeline.

## Project background

In this project, human induced pluripotent stem cell (hiPSC) organoid-derived cone-like (CRX+) or rod-like (NRL+) photoreceptors were transplanted into rd10 mice, a slow retinal degeneration model. Following behavioural assessment at either 12 or 26 weeks after transplantation, retinas were isolated from the eyes and recorded using multi-electrode array (MEA) electrophysiology.

## Repository structure

- **pre-processing**: adapts the HerdingSpikes2 Lightning workflow ([HS2](https://github.com/mhhennig/HS2)), which is based on the SpikeInterface sorting framework. This step performs spike detection, spike extraction, spike sorting and stimulation pulse extraction from raw MEA recordings.
- **analysis**: aligns stimulation time points with spike times, separates noise from electrophysiologically valid sorted units, extracts response metrics, and classifies selected sorted units into ON, OFF, ON–OFF, unconventional and unresponsive response types.
- **plotting**: generates single-recording plots, including sorted unit locations, response subtype percentages, pie charts of classified cell types, and transient/sustained response subtype summaries.

## Pre-processing

Pre-processing adapts the HerdingSpikes2 Lightning workflow (HS2), which is based on the SpikeInterface sorting framework. This step performs spike detection, spike extraction, spike sorting and stimulation pulse extraction from raw MEA recordings.

Raw MEA recordings were collected using the BioCAM DupleX system (3Brain, Lanquart, Switzerland) and saved in .brw format. The final pre-processing output, including spike times and stimulation pulse trains, is saved in _cluster.hdf5 format.

## MATLAB analysis
For each recording session, place the corresponding _cluster.hdf5 file in the single-recording session folder and run: `processRetinaFlashStimData.m`

This script aligns stimulation time points with spike times, separates noise from electrophysiologically valid sorted units, extracts response metrics, and classifies selected sorted units into ON, OFF, ON–OFF, unconventional and unresponsive response types.

The main analysis outputs are:

- `_responseMetrics.mat`
Contains extracted response parameters, including bias index, inter-spike interval (ISI) coefficient of variation, ISI violation rate, tau value, peri-stimulus time histogram (PSTH), post-tau firing rate, and the ON/OFF spike ratio.
- `_totalneuronsV2.mat`
Contains the indices of classified cell clusters, including ON transient, ON sustained, OFF transient, OFF sustained, ON–OFF, unconventional and unresponsive units.
- `_psth.mat`
Contains PSTH values, bin edges and trial-wise PSTH data used for plotting individual cell response traces.

## Indexing note

Individual unit plot filenames are based on Python/HS2 indexing, which starts from 0. MATLAB indexing starts from 1. Therefore, unit numbers in individual plot filenames and MATLAB cell-array indices may differ by 1. When using indices from _totalneuronsV2.mat to access trialPSTHs, confirm whether the stored IDs are 0-based or MATLAB-indexed for the specific pipeline version.

## Input and output data format
### Required input files for downstream plotting

After running processRetinaFlashStimData.m, the following files are required for most downstream plotting functions:

- `_responseMetrics.mat`
- `_totalneuronsV2.mat`
- `_psth.mat`

For spatial plotting and retinal image alignment, the corresponding _cluster.hdf5 file and, where available, the retinal TIFF image are also required.

## Single-recording plotting functions
### Spatial plots

Run:

`plot_locations(clusterFilepath, neuronIDs)`

using the corresponding _cluster.hdf5 file and _totalneuronsV2.mat file.

This function generates:

- a spatial plot of classified RGCs projected onto the MEA layout;
- the percentage of responsive units as a function of distance from the retinal centre, defined here as the optic nerve head;
- `a _distance.mat` file containing spatial coordinates, distance values, binned responsiveness and classified neuron labels for downstream plotting.

## Optional MEA–retina image alignment

If a live image of the retina on the MEA chip was acquired during recording, the image can be aligned to the MEA grid using:

- `chipGridAlign.m`

or, for large TIFF images:

- `chipGridAlign1.m`

These are semi-manual alignment scripts. The user selects four reference coordinates on the chip image: top left, top right, bottom left and bottom right. The script reconstructs the 64 × 64 MEA grid and aligns the retinal image to the electrode layout. The user can then manually select a new reference centre, usually the optic nerve head, for recalculating neuron distances.

## Single-cell raster and PSTH examples

Run:

`plotAllRasterPSTHs_On_Off_Response_example(plotMetrics)`

This function plots raster plots and PSTHs for one example cell across light conditions. The expected plotMetrics structure contains:

| Field | Type | Description |
|-------|------|-------------|
| `trialSpikes` | 1 × nConditions cell array | Trial-wise spike times for each light condition |
| `trialPSTHs` | 1 × nConditions cell array | 1 × nConditions cell array	PSTH matrices with dimensions nTrials × nTimeBins |
| `binEdges` | vector | PSTH bin edges in seconds |

The function generates one raster plot and one PSTH plot per light condition. In the standard flash-stimulation protocol, the light ON window is from -2 to 0 s and the light OFF window is from 0 to 2 s.

## Transplantation group-level folder organisation

For transplantation analysis, organise files by degeneration stage, treatment group and behavioural/electrophysiological subgroup.

Example structure:

rd10_transplanted/
├── early degeneration/
│   ├── NRL/
│   │   ├── behavioural positive and ephys positive/
│   │   ├── behavioural positive and ephys negative/
│   │   └── behavioural negative/
│   ├── CRX/
│   │   ├── behavioural positive and ephys positive/
│   │   ├── behavioural positive and ephys negative/
│   │   └── behavioural negative/
│   └── SHAM/
│
└── late degeneration/
    ├── NRL/
    │   ├── behavioural positive and ephys positive/
    │   ├── behavioural positive and ephys negative/
    │   └── behavioural negative/
    ├── CRX/
    │   ├── behavioural positive and ephys positive/
    │   ├── behavioural positive and ephys negative/
    │   └── behavioural negative/
    └── SHAM/
## Group-level plotting functions
### Cell-type composition across transplantation groups

Run:

`plot_totalneurons_stacked_transplanted2(baseFolder)`

This function recursively searches for _totalneuronsV2.mat files under each experimental subgroup, extracts the number of classified RGCs in each response category, and generates stacked bar plots for early and late degeneration stages.

Response categories include:

- ON neurons
- OFF neurons
- ON–OFF neurons
- unconventional neurons
- non-responsive neurons

The function also exports a summary Excel table containing raw counts and within-group percentages.

# Transplantation Plotting Pipeline

## Delta Firing Rate Plots

First, collect firing-rate metrics using:

```matlab
[T, T1] = collectDeltaFR_transplant_sham(rootDir)
```

The output table `T` contains the following columns:

- `Stage`
- `Group`
- `CellType`
- `LightCondition`
- `UnitID`
- `DeltaFR`
- `PeakStimFR`
- `BaselineMeanFR`
- `FileID`

Then run:

```matlab
PlotDeltaRatePerCondition(T, outputFolder)
```

This function generates boxplot and swarmplot figures for delta firing rate. It creates four sets of plots:

1. ON and OFF cells separately for each light condition.
2. ON and OFF cells separately with all three light conditions combined.
3. All responsive cell types combined for each light condition.
4. All responsive cell types combined across all three light conditions.

The plotted groups are:

- early SHAM
- early NRL
- late SHAM
- late NRL

Pairwise comparisons are performed using **Wilcoxon rank-sum tests**.

## Latency Plots

The same collector function also returns `T1`, which contains latency-related metrics:

- `onsetLatency`
- `peakLatency`
- `peakFR`
- `threshold`

Run:

```matlab
PlotLatencyPerCondition(T1, outputFolder)
```

This function generates boxplot and swarmplot figures for peak latency. The plotting structure mirrors `PlotDeltaRatePerCondition`, including per-condition, combined-condition, ON/OFF-specific and all-responsive-cell plots.

> **Note:** The current plotting function uses `peakLatency` as the latency metric. If `onsetLatency` is required, the plotting variable should be modified accordingly.

## Relationship to the Baseline rd10 Pipeline

The transplantation plotting workflow follows the same logic as the baseline rd10 pipeline:

- Single-recording classification is performed using `processRetinaFlashStimData.m`.
- Classified cell indices are stored in `_totalneuronsV2.mat`.
- PSTH and raster information is stored in `_psth.mat`.
- Spatial information is stored in `_distance.mat`.
- Downstream plotting functions summarise classified response types, spatial responsiveness, firing-rate changes and latency metrics.

The main difference is that transplantation experiments are grouped by **degeneration stage**, **treatment group** and **behavioural/electrophysiological outcome**, rather than by postnatal age alone.
