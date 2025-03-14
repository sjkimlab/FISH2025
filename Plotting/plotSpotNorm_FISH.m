%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 12/14/2024
    Last updated date: 3/11/2025

Description: this script make spotNorm plots for FISH data
1. absolute normalized x position (along cell short axis)
2. normalized L position (along cell long axis)

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
xNormBin = 0.12; lNormBin = 0.0625; % for spotNorm plotting for all timePoints

% set up the figures
figure( 'Position', [400 350 400 350]) % for xNorm
figure( 'Position', [810 350 400 350]) % for LNorm

for i = plotNum % different files
    
    % load FISH analysis file
    load( [list(i).folder '\' list(i).name])

    % ~~~ condition for cell variables ~~~~
    goodCells = 1: totalCells;  extra = ' (all cells)';
    % goodCells = find( cellSpots == 1);   extra = ' (1 spot)';
    % goodCells = find( cellSpots == 2);   extra = ' (2 spot)';

    condSpots = ismember( spotCell, goodCells); % flag for spots in selected cells
    xNorm = spotNorm( condSpots, 1);
    lNorm = spotNorm( condSpots, 2);

    % display # of cells & spots used in analysis
    fprintf( ' ~~~ %s: %5d spots, %5d cells has signal ~~~\n',...
        fishName, sum( cellSpots), sum( cellSpots>0))

    % bootstrap to get errorbar
    getxNorm = @( xNorm) histcounts( abs( xNorm), 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    mX = bootstrp( 1000, getxNorm, xNorm);

    % plot abs( xNorm)
    figure(1)
    [~, edges] = histcounts( abs( xNorm), 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    tmp = movmean( edges, 2);   centers = tmp( 2:end);
    errorbar( centers, mean( mX), std( mX), 'LineWidth', 2.5, 'DisplayName',...
        sprintf( '%s, %d spots', strain, sum( ~isnan( xNorm))))
    hold on

    % plot lNorm
    figure(2)
    [N, edges] = histcounts( lNorm, 'BinWidth', lNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'DisplayName',...
        sprintf( '%s, %d cells', strain, sum( cellSpots>0)))
    hold on
end

fprintf( '\n')

%% Plot Setting

% xNorm
figure(1)
set( gca, 'FontSize', 15)
xlabel( sprintf( 'abs(xNorm), bin=%g', xNormBin))
ylabel( 'Probability')
legend( 'Location', 'northeast', 'FontSize', 12)
ylim( [0 0.3])

% LNorm
figure(2)
set( gca, 'FontSize', 15)
xlabel( sprintf( 'LNorm, bin=%g', lNormBin))
ylabel( 'Probability')
legend( 'Location', 'south', 'FontSize', 12)
ylim( [0 0.13])


