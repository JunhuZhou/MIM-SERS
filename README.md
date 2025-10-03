# Raman Spectra Analysis Pipeline (MATLAB)

This repository provides a MATLAB implementation of a digital Raman spectra analysis pipeline, including baseline correction, classification, and visualization. The baseline correction step is based on the adaptive iteratively reweighted penalized least squares (airPLS) algorithm. An example of Raman spectra dataset is provided.

## Overview
The pipeline performs the following steps on Raman spectral data:
1. Despiking (cosmic ray spike removal)
2. Savitzky–Golay smoothing (noise reduction)
3. Baseline correction (airPLS algorithm)
4. Classification based on peak detection relative to reference background
5. Outlier filtering (removal of overexposed spectra and statistical filtering)
6. Visualization (heatmap, spectral plots, mean ± standard deviation)  

## Requirements
- MATLAB R2019b or later (older versions may also work)
- Signal Processing Toolbox
- Input Raman spectra file in Excel format (.xlsx)

## Input Format
The input file must be an Excel sheet with the following structure:
- Column 1: Raman shift values (in cm⁻¹)
- Column 2 → N: Intensity values for each spectrum

## Output
### The pipeline generates the following files in the same directory as the input file:
- *_Despiked_Spectra.xlsx — spectra after spike removal
- *_Smoothed_Spectra.xlsx — spectra after Savitzky–Golay smoothing
- *_Corrected_Spectra.xlsx — spectra after baseline correction (airPLS)
- *_Classification.xlsx — classification results (0 = negative, 1 = positive)
- *_Mean_STD_Spectra.xlsx — mean spectrum with standard deviation
### Plots generated:
- Heatmap of classification results
- Original vs corrected spectra
- Mean spectrum ± 1 standard deviation

## How to Use:
1. Run `MIM_SERS_Analysis.m`.
2. When prompted, select the input Excel file.
3. The processed spectra and classification results will be saved automatically in the same directory.

## Citation
If you use this code, please cite:
- Zhou, Junhu, Xin Qi, and John XJ Zhang. "Controlled Synthesis of Metal-Insulator-Metal Nanoparticles for Enhanced Raman Spectroscopy." Nanoscale (2025). DOI: https://doi.org/10.1039/D5NR03536H


