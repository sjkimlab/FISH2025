%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 11/28/2023
    Last updated date: 3/11/2025

Description: this script is for plotting bar graph of FISH data spot number
per cell
---------------------------------------------------------------------------
%}

clear, clc, close all

% fishPath = ''; % change to your FISH output folder
fishPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\2025\test\output\';

% find which files to plot
list = dir( [fishPath 'FISH*']);
plotNum = getPlotNum( list);

% plotting parameters
signalFlag = false; % true: cells with signal, false: include cells with no signal

% set up the figures & color
figure( 'Position', [400 400 280 320])
colorList = get( gca,'colororder');

Nstat = cell( length( plotNum),1);  d = 0;

for i = plotNum

    % load FISH analysis file
    load( [list(i).folder '\' list(i).name])

    % ~~~ condition for cell variables ~~~~
    goodCells = 1: totalCells;  extra = ' (all cells)';
    if logical( signalFlag)
        goodCells = find( cellSpots > 0);   extra = ' (1+ spot)';
    end
    cSpots = cellSpots( goodCells); % number of spots (>0) in cells

    % bootstrap to get errorbar
    nReps = 1000;    spotN = 0: 5;  N = nan( nReps, length( spotN));
    for k = 1: nReps
        p = unidrnd( length( cSpots), length( cSpots), 1); % sample with replacement
        sample = cSpots( p);
        [N(k,:), ~] = histcounts( sample, 'binMethod', 'integers', 'binLimit', [-0.5 5.5], 'Normalization', 'probability');
    end

    % plot bar graph
    d = d + 1;
    bar( spotN, mean(N), 1, 'FaceColor', colorList( d,:), 'FaceAlpha', 0.4, 'EdgeColor', 'w', 'LineWidth', 1.5, ...
        'DisplayName', strain), hold on
    errorbar( spotN, mean( N), std( N), 'Color', [0 0 0], 'CapSize', 10,...
        'LineWidth', 1, 'LineStyle', 'none', 'HandleVisibility', 'off')

    % statistics of the bar percentage
    Nstat{ d} = [ mean(N); std(N)];
end

%% Plot Setting
figure( gcf)
set( gca, 'FontSize', 14)
xlabel( 'Number of spots per cell')
ylabel( 'Probability')
legend( 'Location', 'northeast', 'FontSize', 13)

if ~signalFlag
    xlim( [-0.5 4.8])
else
    xlim( [0.5 4.8])
end
ylim( [0 0.6])


