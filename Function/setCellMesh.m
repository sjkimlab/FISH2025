
function cellMeshAll = setCellMesh( cellList, cellListN)
%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 3/11/2025

Description: this function set up cell mesh and remove empty cells
---------------------------------------------------------------------------
%}

    cellMeshAll = {};
    for round = 1: size( cellListN, 2) % number of images

        clear cellMesh

        nCells = size( cellList.meshData{ round}, 2);
        if nCells == 0  % in case no cell are detected
            cellMesh( 1, 1) = struct;
            cellMeshAll = [cellMeshAll; cellMesh];
            continue
        end
        cellMesh( nCells, 1) = struct;  badCells = [];
        
        for Cell = 1: nCells
            
            mesh = cellList.meshData{ round}{ Cell}.mesh; % cell outline: [n,4] - (x1, y1, x2, y2), (right, left)

            if length( mesh) > 10 % some cellMesh only have 6 points (maybe error?)
                meshMid = [ mean( mesh( :, [1 3]), 2), mean( mesh( :, [2 4]), 2)]; % midline along the long axis
                % save the cell mesh info for later use
                cellMesh( Cell).mesh = mesh;
                cellMesh( Cell).meshOut = [ mesh(:, 1:2); flipud( mesh(:, 3:4))]; % reshape the mesh matrix to form a circle [2n, 2]
                cellMesh( Cell).meshMid = meshMid;
                cellMesh( Cell).gridLen = vecnorm( diff( meshMid), 2, 2); % length of each grid (L direction)
                cellMesh( Cell).gridLenCum = cumsum( vecnorm( diff( meshMid), 2, 2)); % cumulative length of each grid (L direction)
                cellMesh( Cell).area = double( polyarea( cellMesh( Cell).meshOut(:,1), cellMesh( Cell).meshOut(:,2)));
                cellMesh( Cell).length = double( sum( cellMesh( Cell).gridLen));
                if Cell > length( cellList.cellId{ round}) % sometimes oufti has this error
                    cellMesh( Cell).cellId = nan;
                    continue
                end
                cellMesh( Cell).cellId = double( cellList.cellId{ round}(Cell));
                
            else
                badCells = [badCells Cell];
                % Sometimes the Oufti went wrong and some cells has no information
                % fprintf( '~~~~~~ Image %d Cell #%d has problem with its mesh! ~~~~~~\n',round, Cell)
            end
        end
        cellMesh( badCells) = [];
        cellMeshAll = [cellMeshAll; cellMesh];
    end
end
