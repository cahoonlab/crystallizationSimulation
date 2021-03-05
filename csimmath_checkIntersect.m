function [ collidedID ] = csimmath_checkIntersect( crystals )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% Takes a struct of crystals uses a subfunction to determine if any
%%% of their edges intersect

collidedID = []; %initialize list of collided crystals

for dd = 1:size([crystals.id],2)
    
    collidedID = horzcat(collidedID, ...
        csimmath_checkIntersectOne(crystals, crystals(dd))...
        );

end; clear dd


collidedID = unique(collidedID); %output

end

