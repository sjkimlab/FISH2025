%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated date: 3/11/2025

Description: this script make plots for FISH data. It plots the 2D Norm
histogram.

FISHdataAnalysis output:
    spotCell: cell number for each spots if it's inside any cell
    spotNorm: [xNorm, lNorm] within the cell for each spots
    cellSpots: the number of spots in each cells, 0 if no signal
    cellArea & cellLength & cellWid: properties of each cells (unit: um)

---------------------------------------------------------------------------
%}    
    
clear, clc, close all

% fishPath = ''; % change to your FISH output folder
fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\2025\test\output\';

% find which files to plot 
list = dir( [fishPath 'FISH*']);
plotNum = getPlotNum( list);

% plotting parameters
nBinY = 5;  nBinX = 3* nBinY; % number of bins for 2D Norm heatmap
maxN( 1, length( plotNum)) = 0; % stores the max value for colorbar plot
norm2DScale = 0.07; % colormap limit

c = 0;

for i =  plotNum % different files
    
    c = c + 1;
    % set up the figures
    figure( 'Position', [ 370*c 500 360 150]) 

    % load FISH analysis file
    load( [list(i).folder '\' list(i).name])

    % ~~~ condition for cell variables ~~~~
    % goodCells = 1: totalCells;  extra = ' (all cells)';
    goodCells = find( cellSpots == 1);   extra = ' (1 spot)';
    goodCells = find( cellSpots == 2);   extra = ' (2 spot)';

    condSpots = ismember( spotCell, goodCells); % flag for spots in selected cells
    xNorm = spotNorm( condSpots, 1);
    lNorm = spotNorm( condSpots, 2);

    % display # of cells & spots used in analysis
    fprintf( ' ~~~ %s: %5d spots, %5d cells has signal, %5d good cells ~~~\n',...
        fishName, sum( cellSpots), sum( cellSpots>0), length( goodCells))
    
    
    % hist for a quarter of cell, flip twice to form whole cell
    [N, xEdges, yEdges] = histcounts2( abs( lNorm-0.5)*2, abs( xNorm)/2, [ nBinX, nBinY],...
        'XBinLimits', [0 1], 'YBinLimits', [0 0.5], 'Normalization', 'probability');
    % N: 15 rows * 5 columns, it should be flipped (transpose') when plotted by pcolor
    % N' is the bottom right part of the matrix, should flip to fill
    
    % flip: upside down,    flip(N,2): left-right
    allN = [ flip( flip(N',2))   flip(N');   
                   flip(N',2)         N' ];
    
    allN( end+1, end+1) = 0; % the color of each grid was determined by the first [X,Y] 

    [l, r] = size( N);
    X = -l: 1: l;   Y = -r: 1: r;        
    
    s = pcolor( X, Y, allN); hold on 

    s.LineWidth = 1; % edge line width
    set( gca, 'xtick',[], 'ytick',[]) % no tick or labels
    caxis( [0 norm2DScale]) % set colormap limit
    axis image, colormap jet

    maxN( c) = max( N(:));
    
    % plot cell mesh
    rMesh = r*0.95;     th = (0: pi/100 : pi) + pi/2;
    x = [ rMesh*cos( th)+ r-l,  rMesh*cos( th+pi)+ l-r r-l];
    y = [ rMesh* sin( th),  rMesh* sin( th+pi) rMesh];
    plot( x, y, 'w', 'lineWidth', 3)       
        
    
    % Plot Setting        
    % title( sprintf( '%s %s', strain, Date), 'FontSize', 16)
    title( sprintf( '%s %s', strain, extra), 'FontSize', 16)
    fprintf( '\n~~~  maxN for colorbar should be %.3f ~~~\n', max( maxN(:)))
    
    fprintf( '\n')
end

% colorbar
% set( gca, 'FontName', 'Arial')
