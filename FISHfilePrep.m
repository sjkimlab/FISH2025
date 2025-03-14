%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/28/2023
    Last updated date: 3/10/2025

Description: this script organized the FISH tiff files output into
corresponding folder for phase, Cy5, Cy3. Then it creates folders for Cy3
specifically with renamed files for uTrack analysis

The .nd files are converted into tiff files using Nikon NIS Elements software.

The tiff files are named based on their channel numbers:
    c1 - Cy5 | c2 - Cy3 | c3 - phase
This script moves the corresponding tiff files into the corresponding folders
    2Cy5 | 3Cy3 | 1ph
---------------------------------------------------------------------------
%}


% naming rules: epi001.nd2, epi010.nd2, etc.

% when exported as mono TIFF images for each channel, the output will be
% named as epi001c1, epi001c2, epi001c3.tif

% ----- please go to the directory that stores nd2 images -----


%% 1. Organize tif files in different channels

clear, clc

% folderName corresponding to channel 1, 2, 3
folderName = { '2Cy5' '3Cy3' '1ph'}; 

% loop through all 3 channels
for c = 1: 3
    
    % split images from different channels into folders
    listTif = dir( sprintf( 'epi*c%d*', c)); % find images of this channel
    mkdir( folderName{c}) % create a folder for this channel
    
    for k = 1: length( listTif)
        
        name = listTif(k).name; % original name from FISH imaging (e.g. epi001c1.tif)
        movefile( name, [ folderName{c} '\' name]) % move into 2Cy5 Folder
    end
    fprintf( '~~~  %s tiff files Moved  ~~~\n\n', folderName{c})
end


%% 2. Rename the epi Cy3 images for uTrack analysis

    path = [ pwd '\'];
    folderName = 'Cy3_uTrack';
    
    % copy the epi Cy3 channel to another folder named Cy3_time
    copyfile( [path '\3Cy3'], [path folderName]);  cd( [path folderName])
    
    % rename epi files so that uTrack can treat them as a time-stack (I
    % only use this to run uTrack spotDetection more efficiently, but not
    % link them to form tracks)
    %       epi001c2.tif --> epiCy3t01.tif
    %       epi002c2.tif --> epiCy3t02.tif    
    %
    listEpi = dir( 'epi*');
    for k = 1: length( listEpi)        
        name = listEpi(k).name; % original name of FISH file (e.g. epi001c2.tif)
        newName = sprintf( 'epiCy3t%02.f.tif', k); 
        movefile( name, newName)
        fprintf( '~~~  %s -> %s  ~~~\n', name, newName)
    end
    
    fprintf( '~~~~~~    File Conversion Complete    ~~~~~~\n\n')

% change directory to parent folder
cd ..\


%% 3. Run uTrack on the FISH Cy3 data for spot detection

fprintf( '~~~ Run uTrack, select input channel folder ''Cy3_uTrack'', set the output folder the same ~~~\n\n')
% uTrack output will be saved to the same folder containing Cy3 tiff files

% for FISH data, I'm currently using following parameter fow spot detection
%   Gaussian std = 1.28 pixel (Hamamatsu camera)
%   Alpha (Detection) = 0.015
%   Alpha (fitting) = 0.05 (default)

% movieSelectorGUI % launch u-track