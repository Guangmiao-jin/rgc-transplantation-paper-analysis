# rgc-transplantation-paper-analysis

This repository contains HerdingSpikes2-based Python pre-processing scripts and MATLAB analysis scripts for RGC light-response classification and transplantation-related electrophysiology analyses.

## Contributors

This analysis pipeline was developed by Guangmiao Jin, with contributions from Dr Michael Savage.

- **Guangmiao Jin** developed the overall pipeline structure, selected and adjusted key analysis and classification parameters, implemented the RGC response classification workflow, and organised the plotting comparison analyses between transplanted group and SHAM group.
- **Dr Michael Savage** contributed to the Python stimulation pulse extraction script, MATLAB code optimisation, code annotation and commenting, and technical review of the analysis pipeline.

## Project background

In this project, human induced pluripotent stem cells (hiPSC) organoids-derived cones (CRX+) or rods (NRL+) were transplanted into rd10 mice which represent a slow degeneration model. After behavioural test after either 12 weeks or 26 weeks, multi-electrode array (MEA) recording was conducted on the transplanted mice retinas which were isolated from the eyes.

## Repository structure

- **pre-processing**: adapts the HerdingSpikes2 Lightning workflow ([HS2](https://github.com/mhhennig/HS2)), which is based on the SpikeInterface sorting framework. This step performs spike detection, spike extraction, spike sorting and stimulation pulse extraction from raw MEA recordings.
- **analysis**: aligns stimulation time points with spike times, separates noise from electrophysiologically valid sorted units, extracts response metrics, and classifies selected sorted units into ON, OFF, ON–OFF, unconventional and unresponsive response types.
- **plotting**: generates single-recording plots, including sorted unit locations, response subtype percentages, pie charts of classified cell types, and transient/sustained response subtype summaries.
