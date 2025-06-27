# Raman Spectra Analysis Pipeline

This MATLAB script automates the processing of Raman spectra:

1. Despiking  
2. Savitzky-Golay smoothing  
3. airPLS baseline correction  
4. Peak-based classification  
5. Heatmap of results  
6. Statistical analysis (mean ± standard deviation)  
7. High-quality plots  

### Files Generated:
- `*_Despiked_Spectra.xlsx`  
- `*_Smoothed_Spectra.xlsx`  
- `*_Corrected_Spectra.xlsx`  
- `*_Classification.xlsx`  
- `*_Mean_STD_Spectra.xlsx`

### Requirements:
- MATLAB R2020b+  
- Signal Processing Toolbox  

### How to Use:
1. Run `raman_analysis_pipeline.m`.
2. Select an `.xlsx` file with Raman shift in column 1 and spectra in columns 2–101.
3. Outputs saved in the same folder.

### Reference:
Zhang, Z.-M., Chen, S., & Liang, Y.-Z. (2010). *Baseline correction using adaptive iteratively reweighted penalized least squares*. Analyst, 135(5), 1138–1146.

### License:
MIT
