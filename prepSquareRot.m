function [ unique_slopes, yint1, yint2] = prepSquareRot( vertices, angle, center )
%PREPSQUAREROT Prepares the mathematics needed to draw a rotated square on
%a pixel image.

roundingthreshold = 0.0000001;
rotationArray = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];

% rotate the square
rotatedArray = rotationArray * (vertices - center) + center;

% find slopes and intercepts of each of the edges
lr1 = [ones(2,1) rotatedArray(1,1:2)']\rotatedArray(2,1:2)'; % [yint; slope]
lr2 = [ones(2,1) rotatedArray(1,2:3)']\rotatedArray(2,2:3)'; % [yint; slope]
lr3 = [ones(2,1) rotatedArray(1,3:4)']\rotatedArray(2,3:4)'; % [yint; slope]
lr4 = [ones(2,1) [rotatedArray(1,4) rotatedArray(1,1)]']\[rotatedArray(2,4) rotatedArray(2,1)]'; % [yint; slope]

lr = horzcat(lr1,lr2,lr3,lr4);

slopes = lr(2,:);
unique_slopes = uniquetol(slopes,roundingthreshold); %is also ordered, so negative slope is slope1
yints = lr(1,:);

yint1 = yints(abs(slopes - unique_slopes(1)) < roundingthreshold); % a pair of yintercepts for slope 1
yint2 = yints(abs(slopes - unique_slopes(2)) < roundingthreshold); % a pair of yintercepts for slope 2


end

