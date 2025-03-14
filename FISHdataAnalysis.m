%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 3/11/2025

Description: this script is for FISH data analysis, please run it section
by section.

====== Output =======
spotCell: cell number for each spots if it's inside any cell, nan for none
spotNorm: [xNorm, lNorm] within the cell for each spots, [nan nan] for none
spotAmp: signal amplitude of all spots
cellSpots: the number of spots in each cells, 0 if no signal
cellSpotsStat: unique cell ID vs number of spots
cellArea & cellLength & cellWid: properties of each cells (unit: um)

---------------------------------------------------------------------------
----------- please go to the directory that stores nd2 images -------------
---------------------------------------------------------------------------
%}


%% 1. Combine uTrack & oufti output and save into files 

clear, clc, close all

% specify the strain name
Date = input( ' Please input the Date of the experiment (like: 250310):  ', 's'); fprintf( '\n')
strain = input( ' Please input the strain number (like: SK519):  ', 's'); fprintf( '\n')
% Date = '250311';    strain = 'SK591';

% set the output folder to store analsis result
% fishPath = ''; % change to your folder
fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\2025\test\output\';

if ~exist( fishPath, 'dir'), mkdir( fishPath), end % create save folder

% find uTrack output folders
uTrackFolder = dir( 'Cy3_uTrack*');
    
    % load uTrack output - movieInfo (containes spot coordinates)
    load( [ uTrackFolder.name '\TrackingPackage\GaussianMixtureModels\detections_for_channel_1\Channel_1_detection_result'],...
        'movieInfo')
    
    tmp = split( uTrackFolder.name, '_');    % e.g. 'Cy3_uTrack'
    
    % load oufti mesh output (cell meshes, e.g. 'test-mesha.mat')
    meshFile = dir( '*-mesh.mat'); load( meshFile.name, 'cellList', 'cellListN')
    
    % set up cell mesh and remove empty cells
    cellMeshAll = setCellMesh( cellList, cellListN);
    cellListN = cellfun( @length, cellMeshAll); % number of cell in each image
    
    totalSpots = sum( cellfun( @length, {movieInfo.amp}'));
    totalCells = sum( cellListN);
    fprintf( '~~~  %8s |  spot number:%5d,  cell number:%5d', uTrackFolder.name, totalSpots, totalCells)
    
    % save uTrack spot output 'movieInfo' & other info into Cy3_spotsMesh_time
    save( 'Cy3_spotsMesh', '-regexp', '^(?!(tmp|cellList)$).') % save variables except 'tmp' & 'cellList'
    
    fprintf( '   ~~~   spotsMesh Saved   ~~~\n')
    
fprintf( '\n~~~~~~ uTrack & oufti output extracted and saved to Cy3_spotsMesh ~~~~~~\n\n')


%% 2. Calculate the spotNorm & other quantities, Save to files

clear

pixelSize = 64.5; % unit: nm, change to the pixel size of your microscope

% load spotsMesh file
list = dir( 'Cy3_spotsMesh.mat');   load( list.name)
    
spotCellAll = {};   spotNormAll = {};   spotAmpAll = {};
cellArea = [];      cellLength = [];    cellSpots = [];   cellSpotStat = [];

cellNcum = cumsum( [ 0; cellListN]); % this is for assigning unique cell number

for i = 1: length( movieInfo) % image number
    
    cellMesh = cellMeshAll{ i};
    % skip this image if no cells are detected
    if isempty( fieldnames( cellMesh)), continue, end
    
    spotPos = [movieInfo(i).xCoord(:,1) movieInfo(i).yCoord(:,1)]; % [x, y] coordinate
    spotAmp = movieInfo(i).amp(:,1); % spot amplitude
            
    nSpots = length( spotAmp);  % total number of spots detected
    spotCell = nan( nSpots, 1); % which cell the spots are in
    spotNorm = nan( nSpots, 2); % normalized position of spots [lNorm, xNorm]
    
    for Cell = 1: cellListN( i) % number of cells in this image            
        
        % a. find which cell the spots belong to
        meshOut = double( cellMesh( Cell).meshOut);
            % compare spot positions & mesh to find which cells they belong to
            inCell = inpolygon( spotPos(:,1), spotPos(:,2), meshOut(:,1), meshOut(:,2)); % return 0 or 1 in a vector
            spotCell( inCell) = Cell; % assign the spots to cell
        
        % b. find normalized position of the spots in cells
        badCell = false;
        plotFlag = false; % flag for normalized positions plots
        for spotNum = find( inCell)'
            pt = spotPos( spotNum, :); % spot coordinate
            [ spotNorm( spotNum, :), badCell] = findNormPos( pt, cellMesh( Cell), plotFlag);
        end
        
        if badCell                
            warning( '   ~~~ Image #%d Cell #%d mesh curvature has problem ~~~\n', i, cellMesh( Cell).cellId)
        end
        if cellMesh( Cell).area/ cellMesh( Cell).length * pixelSize/1000 < 0.3
            warning( '   ~~~ Image #%d Cell #%d width is abnormal (too thin) ~~~\n', i, cellMesh( Cell).cellId)
        end
    end
    
    % count the # of spots in each cell (integer bins)
    [N, ~] = histcounts( spotCell, 'binLimit', [0.5 length( cellMesh)+0.5],...
        'binMethod', 'integers');
    % histogram( N)

    % store the information from this image and stack them together,
    % cell number adds up over images
    spotCellAll = [ spotCellAll; spotCell + cellNcum(i)]; % cell number assigned to spots (accumulated)
    spotNormAll = [ spotNormAll; spotNorm]; % spot normalized position
    spotAmpAll = [ spotAmpAll; spotAmp]; % spot intensity
    cellArea = [ cellArea; [cellMesh.area]'];
    cellLength = [cellLength; [cellMesh.length]'];
    cellSpots = [ cellSpots; N'];        % # of spots in all cells
    cellSpotStat{i,1} = [ [cellMesh.cellId]' N']; % count each cell's spot number for all movies
end

cellArea = cellArea* pixelSize^2/ 10^6;     % unit: um^2
cellLength = cellLength * pixelSize/ 1000;  % unit: um
cellWid = cellArea./ cellLength;            % unit: um

spotCell = cell2mat( spotCellAll); % cell number assigned to spots (accumulated)
spotNorm = cell2mat( spotNormAll); % spot normalized position
spotAmp  = cell2mat( spotAmpAll);  % spot intensity

% save the analysis result
fishName = [ 'FISH ' strain ' ' Date];
save( [fishPath fishName],  'pixelSize', 'totalSpots', 'totalCells', ...
    'cellArea', 'cellLength', 'cellWid', 'cellSpots', 'cellSpotStat', ...
    'spotCell', 'spotNorm', 'spotAmp', 'Date', 'strain', 'fishPath', 'fishName')

fprintf( '~~~  FISH spotNorm calculation completed & saved  ~~~\n')

fprintf( '\n~~~~~~ Analysis Done ~~~~~~\n\n')


