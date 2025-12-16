## About

This repo contains analysis scripts and analysis results for smFISH experiment. I also provided a test dataset that you could run the analysis on your own.



#### System Requirement

1. MATLAB (R2024a or later)
2. Oufti (compiled version)
3. u-track (v2.3) 



#### Analysis Code

- *FISHfilePrep.m*: this script sorts raw TIFF files into dedicated folders, preparing the data for subsequent oufti & u-track analysis.
- *FISHdataAnalysis.m*: this script combines the output from oufti & u-track analyses and performs further analysis, such as assigning foci to cells and calculating normalized cellular position.

#### Plotting Code

find the scripts in the *Plotting* subfolder 

- *plotSpotNorm_FISH.m*: plot normalized spot localization along both cell axes.
- *plot2DNorm_FISH.m*: plot 2D histogram of normalized spot localizations.
- *plotSpotN_FISH.m*: plot histogram of spot numbers per cell.

#### DataSet

- test: folder containing example dataset for 2 strains (SK519, SK591). 
  - Each subfolder contains raw .nd2 files and TIFF files exported using NIS-Element AR software. Note that each .nd2 image consists of 3 channels: Cy5, Cy3, and phase contrast, so there're 3 TIFF images corresponding to each .nd2 image.
- test_analyzed: the same test folder but after running the analysis
  - TIFF files are sorted into dedicated folders: 1ph, 2Cy5, 3Cy3. The Cy3_uTrack folder contains the uTrack analysis results. The oufti results are saved as *SK519-mesh.mat*. Cy3_spotsMesh contains the combined data from oufti (cell mesh) and uTrack (spot position) analyses.
- output: folder contains the final analysis results and plots.
  - plots: subfolder that contains the example plots of the analysis.
  - Final analysis results: *FISH SK519 250311.mat, FISH SK591 250311.mat*
    - These result files are used by plotting scripts for making plots, such as 2D histogram of spot localization, spot number histogram, normalized spot localization along cell axes.



## Running Workflow

1. Download MATLAB and add the **FISH2025** folder and its subfolder to MATLAB path.
2. Change the current working directory to the folder that contains the raw dataset. For example, enter `cd( 'yourPath\FISH2025\Dataset\test\SK519')` in the MATLAB command window
3. run *FISHfilePrep.m* to sort the TIFF files and prepare them for subsequent oufti & u-track analysis.
4. run Oufti for the phase contrast images and save the results in the same folder as *SK519-mesh.mat*
5. run u-track for the **Cy3_uTrack** folder to detect the foci positions and save the results inside the **Cy3_uTrack** folder
6. run *FISHdataAnalysis.m* to analyze the data, make sure the working directory is still the folder containing the raw dataset as in Step 2. The script will save the analysis results in the folder `'..\output'`. The saved file is named as `['FISH ' strain ' ' Date]`.
7. Make plots. Change the current working directory to the parent folder that contains **output** folder. This folder was created in Step 6. Use scripts in the **Plotting** folder to visualize analysis results (2D histogram of spot localization, spot number histogram, normalized spot localization along cell axes)
