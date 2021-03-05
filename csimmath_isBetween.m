function isBetween = csimmath_isBetween( num, limits )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% A simple function that takes a number and determines if it is between
%%% two numbers
%%% Outputs a boolean

% num should be a number
% limits should be a list [num num]

isBetween = num >= min(limits) && num <= max(limits);

end

