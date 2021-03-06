function collidedID = csimmath_checkIntersectOne( crystals, suspect )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% Takes a struct of crystals and one suspect crystal
%%% Determines if the suspect intersects with any others

collidedID = []; %initialize list of all colliding crystals

suspectLine = {suspect.m1, suspect.b11, suspect.x11;...
    suspect.m1, suspect.b12, suspect.x12;...
    suspect.m2, suspect.b21, suspect.x21;...
    suspect.m2, suspect.b22, suspect.x22};

%avoid unecessary work by first determining proximity
%first part is distance between centers
%second part is length of both half-diagonals
spaceBetween = sqrt(...
    ([crystals.Cx] - suspect.Cx).^2 + ([crystals.Cy] - suspect.Cy).^2 ...
    )...
    - sqrt(2).*[crystals.L]./2 - sqrt(2)*suspect.L/2;


for cc = 1:size([crystals.id],2) %iterate through crystals

    %first check proximity (and don't compare identical crystals)
    if spaceBetween(cc) <= 0 && ~eq(crystals(cc).id, suspect.id)
        %four edges, do this four times
        for tt = 1:4
            if csimmath_linesIntersect(...
                    suspectLine{tt,1}, suspectLine{tt,2}, suspectLine{tt,3},...
                    crystals(cc).m1, crystals(cc).b11, crystals(cc).x11)...
                    ||...
                    csimmath_linesIntersect(...
                    suspectLine{tt,1}, suspectLine{tt,2}, suspectLine{tt,3},...
                    crystals(cc).m1, crystals(cc).b12, crystals(cc).x12)...
                    ||...
                    csimmath_linesIntersect(...
                    suspectLine{tt,1}, suspectLine{tt,2}, suspectLine{tt,3},...
                    crystals(cc).m2, crystals(cc).b21, crystals(cc).x21)...
                    ||...
                    csimmath_linesIntersect(...
                    suspectLine{tt,1}, suspectLine{tt,2}, suspectLine{tt,3},...
                    crystals(cc).m2, crystals(cc).b22, crystals(cc).x22)

                collidedID(end+1) = crystals(cc).id; %append
            end    
        end
    end
end


collidedID = unique(collidedID); %output list of ids


end

