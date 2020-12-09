function [ grainLogical ] = getGrainLogical( grainInfo, x, y )
%GETGRAINLOGICAL creates a logical variable describing which pixels should be changed
%   called by nucGrainSim_vX

%get the centers
x_cent = grainInfo.posX;
y_cent = grainInfo.posY;

%get the new corners
y1 = y_cent - round(grainInfo.size/2);
y2 = y_cent + round(grainInfo.size/2);
x1 = x_cent - round(grainInfo.size/2);
x2 = x_cent + round(grainInfo.size/2);

rotAngle = grainInfo.rot;

if rem(rotAngle,90)==0 %if there is no rotation
    grainLogical = y >= y1 & y <= y2 & x >= x1 & x <= x2;

else %most cases there will be some rotation
    vertices = [x1, x1, x2, x2; y1, y2, y2, y1];
    center = repmat([x_cent y_cent],4,1)' ;

    %general error handler
    try
        [unique_slopes,yint1,yint2] = prepSquareRot(vertices,rotAngle,center);
    catch ME
        assignin('base','grainError',grainInfo)
        warning('there was an error in nucGrowthSim. See assignin base grainError')
    end

    %easiest way I found to locate the indexes of a rotated box
    grainLogical = y >= unique_slopes(1).*x + min(yint1) & y <= unique_slopes(1).*x + max(yint1) & y >= unique_slopes(2).*x + min(yint2) & y <= unique_slopes(2).*x + max(yint2);

end


end

