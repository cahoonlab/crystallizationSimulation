function [ soothesayerSays ] = soothesayer( grains, growingGrains, soothesayer_distance, soothesayer_grainSize, x, y, soothesayer_image, color_grainActive, color_grainFinished )


for gg = growingGrains %go through all the active grains in the struct
    %determine which part of the image to change.
    soothesayer_grainInfo = grains(gg); 
    soothesayer_grainInfo.size = soothesayer_grainSize; %don't actually want to change the grain struct
    grainLogical = getGrainLogical(soothesayer_grainInfo, x, y);

    %Use the logical to change the color of the newly grown grain.
    %Not sure yet if there has been a collision
    soothesayer_image(grainLogical) = color_grainActive; 
end %done going through all active grains


%now determine if there were grain collisions
allGrains = soothesayer_image == color_grainActive | soothesayer_image == color_grainFinished;
L = bwlabeln(allGrains, 8); %label spots (all grains)
S = regionprops(L, 'Area'); %measure area of spots
listofspotnumbersbygrain = zeros(length(grains),1); %initialize list of all the spot numbers

for ss = 1:length(S) %go through all spots
    
    for gg = 1:length(grains) %need to go through all grains
        listofspotnumbersbygrain(gg) = L(grains(gg).posY,grains(gg).posX); %get the spot number for each grain (based on center point)
    end; clear gg
    
end; clear ss

for gg = growingGrains
    if length(find(listofspotnumbersbygrain == listofspotnumbersbygrain(gg))) > 1 %if the future involves collision
%         soothesayer_image(L == listofspotnumbersbygrain(gg)) = color_grainFinished; %color the connected grains dormant colored.
        continue %why bother with the above? Keep for debugging.
    else %no collision in the future for this grain.
        grains(gg).waitUntil = soothesayer_distance;
    end
end; clear gg

soothesayerSays = [grains.waitUntil];

end

