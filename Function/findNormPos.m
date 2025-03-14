
function [spotNorm, badCell] = findNormPos( pt, cellMesh, plotFlag)
%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 8/9/2023
    Last updated at 3/11/2025

Description: this function calculates the subcellular position of spots
---------------------------------------------------------------------------
%}

    badCell = false;
    mesh = cellMesh.mesh;
    meshMid = cellMesh.meshMid;
    lenCum = [0; cellMesh.gridLenCum]; % add 0, so that the index matches with the mesh
    
    % calculate lNorm, match idx with mesh
    [~, dist] = findPerpFoot( pt, mesh(:, 3:4), mesh(:, 1:2)); % -: below, +: above, should be - to +
    bra = find( abs( diff( sign( dist))) == 2); % find the grid idx where distance changes sign, first & end = nan
    % dist(bra) & dist(ket) are the distance of the point to its neighboring minor axis

    % Parallel Methods
    if isempty( bra) % not sandwiched between 2 segment edges, must be in two caps then 
        if sum( dist< 0) == 0       % all positive, it's in the 1st segment
            bra = 1;
        elseif sum( dist> 0) == 0   % all negative, it's in the last segment                    
            bra = length( dist)- 1;
        else                        % the point lands on the segment!
            spotNorm = [nan nan];
            return
        end
    elseif length( bra) > 1 % bad cell, outline curls back
        badCell = true;
        spotNorm = [nan nan];
        return
    end
    ket = bra + 1;
    
    % use parallel line of segment edges to find intersection with midline
    % (LNorm) & outline (xNorm), modified @7/23/2023 by YHW 
    a = mesh( [bra ket],:);
    
    % find the mean slope of the neighboring two segments
    slopes = (a(:,4)-a(:,2))./ (a(:,3)-a(:,1));
    if sum( sign( slopes)) == 0 % two slopes have opposite sign 
        k = abs( diff( slopes))* sign( sum( slopes));
    else
        k = mean( slopes, 'omitnan');
    end
    
    [D, xPos] = findIntersect( pt, k, meshMid(ket,:), meshMid(bra,:));
    % find the intersection point of pt-D & the cell outline (cell width at lPos for this pt)
    if xPos > 0 % on the right side of the mid line
        [IntPt, ~] = findIntersect( pt, k, mesh( bra, 3:4), mesh( ket, 3:4));
    else        % on the left side of the mid line
        [IntPt, ~] = findIntersect( pt, k, mesh( bra, 1:2), mesh( ket, 1:2));
    end            
    xNorm = xPos/ norm( IntPt-D); % normalized by the width at that point                
    lPos = lenCum( bra) + norm( D- meshMid( bra,:)); % real L value by portion
    lNorm = lPos/ lenCum( end); % normalized by the total length of the major axis
    
    if abs(xNorm) > 1
        warning( 'xNorm outside the [-1 1] region, Cell %d', cellMesh.cellId)         
        xNorm = nan;
        plotFlag = true;
    end
    if abs( lNorm-0.5) > 0.5
        warning( 'lNorm outside the [0 1] region, Cell %d', cellMesh.cellId)
        lNorm = nan;
        plotFlag = true;
    end
    
    spotNorm = [xNorm, lNorm];
    
    if plotFlag
        close, figure()
        meshOut = cellMesh.meshOut;        
        plot( meshOut(:,1), meshOut(:,2), 'b', 'LineWidth', 1), hold on
        plot( meshMid(:,1), meshMid(:,2), 'b', 'LineWidth', 1)
        for n = 1: length( mesh)-1
            plot( mesh( n, [1 3]), mesh( n, [2 4]), 'b', 'LineWidth', 1)
        end

        % position along the long axis
        plot( meshMid(:,1), meshMid(:,2), 'k', 'LineWidth', 7) % total length
        plot( [ meshMid( 1:bra, 1); D(1)], [ meshMid( 1:bra, 2); D(2)],...
            'r', 'LineWidth', 3)                    
        plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 3)

        % position along the short axis
        plot( [IntPt(1) D(1)], [IntPt(2) D(2)], 'm', 'LineWidth', 10)
        plot( [D(1) pt(1)], [D(2) pt(2)], 'w', 'LineWidth', 3)                    
        scatter( IntPt(1), IntPt(2), 100, 'c', 'filled') % signal point
        scatter( D(1), D(2), 100, 'c', 'filled') % perpendicular foot
        scatter( pt(1), pt(2), 100, 'w', 'filled') % signal point
        axis equal
        figure( gcf)
        hold off
        input( '~~~~ Enter to view next cell: ');
    end
end

function [D, dist] = findPerpFoot( pt, B, C)
% this function returns the coordinate of the perpendicular foot D so that 
% AD perpendicular to BC (everything in 2D) and the distance of pt to line BC
% dist > 0 if pt is on the right side of line BC (pt, B, C: counterclockwise)
% dist < 0 if pt is on the  left side of line BC (pt, B, C: clockwise)

    AB = B - pt; % vector
    BC = C - B;  % vector

    area = AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1); % cross product
    side = sign( area); % -1: left side, +1: right side

    normVec = [ BC(:,2) -BC(:,1)]; % normal vector of BC in 2D 
    unitNormVec = normVec./ vecnorm( normVec, 2, 2); % unit normal vector of BC
    AD = dot(unitNormVec, AB, 2).* unitNormVec; % AD is perpendicular to BC, dot product
    D = pt + AD; % D point of intersection, Perpendicular Foot
    dist = side.* vecnorm( AD, 2, 2);
end    


function [IntPt, dist] = findIntersect( pt, k1, B, C)
% this function find the intersection point between two lines, one line is
% determined by a point pt and a slope k1, the other line is determined by
% two points B & C

    IntPt = nan( 1, 2);
    b1 = pt(2) - k1*pt(1); % pass through pt

    % kx1 + b = y1
    % kx2 + b = y2    
    % k = (y2-y1)/ (x2-x1);  b = y1 - kx1;
    k2 = ( C(2)- B(2))/ (C(1)- B(1));
    b2 = C(2) - k2* C(1);

    % k1x + b1 = y
    % k2x + b2 = y
    % x = - (b1-b2)/ (k1-k2)
    % y = k1x + b1 = k2x + b2
    IntPt(1) = -( b1-b2)/ (k1-k2);
    IntPt(2) = k2* IntPt(1) + b2;

    if sum( isnan(IntPt)) || sum( isinf(IntPt))
        % fprintf('    ~~~ Intersection has problem! ~~~\n')
        IntPt = [nan nan];
    end

    % cross product, right hand rule opposite, +1: clockwise, -1: counter        
    AB = B - pt; % vector AB
    BC = C - B;  % vector BC
    side = sign( AB(:,1).*BC(:,2) - AB(:,2).*BC(:,1)); 
    dist = side.* vecnorm( IntPt-pt, 2, 2);
end