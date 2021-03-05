function embedded = csimmath_checkEmbed( crystals,suspect )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% Takes a list of crystals and a suspect
%%% Determines if the suspect is embedded in any other crystals
%%% Outputs boolean

% %% a test environment
% crystals = placeCrystal(...
%     simSummary, 'Cx', 0, 'Cy', 0, 'rot', 30);
% crystals(2) = placeCrystal(...
%     simSummary, 'Cx', 3e-5, 'Cy', 0, 'rot', 45); crystals(2).id = 2;


%% do function

%suspect center position
embryoCx = suspect.Cx;
embryoCy = suspect.Cy;

embedded = false; %default


%iterate through the list of crystals
for ii = 1:size([crystals.id],2)
    
    %skip the entry for the identical crystal
    if ~eq(crystals(ii).id, suspect.id)
        
        %is the center between the two y lines two x lines?
        if csimmath_isBetween(embryoCy, sort(...
                crystals(ii).m1 * embryoCx +...
                [crystals(ii).b11, crystals(ii).b12])...
                ) &&...
                csimmath_isBetween(embryoCx, sort(...
                [embryoCy - crystals(ii).b21,...
                embryoCy - crystals(ii).b22] ./ crystals(ii).m2)...
                )
        
            embedded = true;
            break %don't bother testing the rest
        end
    end
end; clear ii


end

