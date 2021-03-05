function [ doesIntersect ] = csimmath_linesIntersect( line1m, line1b, line1x, line2m, line2b, line2x )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% Takes two lines (slope, intercept, and x limits)
%%% Determines if those two lines intersect

% %% a test environment
% line1m = 0;
% line1b = 12;
% line1x = [7, 10];
% 
% line2m = -1;
% line2b = 15;
% line2x = [2, 12];
% 
% plotX = [line1x(1):0.001:line1x(2)];
% figure; hold on
% plot(plotX,line1m .* plotX + line1b)
% plot(plotX,line2m .* plotX + line2b)


%% perform the function

%if you set y to be equal, what is the x value?
commonX = (line2b - line1b) / (line1m - line2m);

%is the x value in the region of both x limits?
doesIntersect = csimmath_isBetween(commonX, line1x) && csimmath_isBetween(commonX, line2x);


end

