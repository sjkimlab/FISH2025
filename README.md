## FISH2025

This repo contains analysis code and test dataset with analysis results for smFISH experiment. 

Please add the folder and its subfolder to MATLAB path before running.

1. Sort image files after exporting nd2 files into TIFF files. Images from each channel (phase contrast, fluorescence 1, fluorescence 2) should be separated into different folders. This step is done by *FISHfilePrep.m*. This script also prepared your fluorescence channel to be ready for u-track analysis in step 3.
2. Analyze cell phase contrast images to get cell outlines using Oufti
3. Identify spots in fluorescence images using u-track. 
4. Combine Oufti and u-track analysis results and perform further analysis. This step is done by *FISHdataAnalysis.m*
5. Make plots. Use scripts in the 'Plotting' folder to visualize analysis results (2D histogram of spot localization, spot number histogram, normalized spot localization along cell axes)

