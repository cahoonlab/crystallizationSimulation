function [simParams, grains, summary] = nucGrowthSim_v12 ( nuc_density_mm2, varargin )
%new in v3: rotation
%new in v4: measure area to determine if collisions have taken place
%new in v5: predict the future
%new in v6: predict the future better
%new in v7: include particles in the outer region that affect final growth size but aren't included in the sample statistics (pseudo periodic boundary conditions)
%new in v8: different optimizer for size based on improving resolution but
    %staying within numGrains and pixel limits
%new in v9: a test to see if a large initial soothesayer sweep would speed
    %things up. It didn't.
%new in v10: a completely different system. Went grain by grain, trying
    %timesteps between the active and dormant timestep until it found the
    %collision. I was hoping this would be faster, but because every grain has
    %to draw the entire landscape from scratch (including all dormant grains),
    %it was actually rather slow.
%new in v11: a copy of v8, with some implementation of v10 inside of it.
    %We're back to iterating chronalogically, but the soothesayer attempts to be
    %smarter by testing at midpoints.
%new in v12: realized the initial soothesayer wasn't behaving as expected.
    %Now it just gives a good initial guess of an end point and adapts as
    %needed.
simVersion = 12;


%%%%% NOTE: nucleation density should be in #/mm^2

%% user defined stuff

%default values (actually defined in varargin)
saveIm = false; %do you want to save image frames of the growth? (this saves ALL frames for purpose of animation)
pxLimit = [5000, 10000]; %how large would you like your final image? This is the dimension for one edge. The arrays get pretty large!
growthRate_px = 1; %in pixels (smaller value will give more exact result but take longer). I don't do anything besides 1.
startSize_px = 1; %in pixels (smaller value will get more exact result but might be more correct. Maybe?)
permitSoothesaying = false; %do you want to speed up the simulation? No if saving frames for animation.
soothesaying_limit = 10; %keep soothesayer going until deltaTimeStep of this value. Somewhere between 5-30 seems good.
bufferWidthMultiplier = 0.5; %how thick of buffer around the ROI? ex. 0.5 would put 1/2 of the ROI width on each side of the ROI
grainOptimizerLimits = [50 100]; %How many grains do you allow in your simulation? (Inside the ROI, specifically)



%% get the user defined values
for vv = 1:nargin-1 %-1 because of the 1 required field
    if strcmpi(varargin{vv},'saveIm')
        saveIm = varargin{vv+1};
    elseif strcmpi(varargin{vv},'pxLimit')
        pxLimit = varargin{vv+1};
    elseif strcmpi(varargin{vv},'growthrate_px')
        growthRate_px = varargin{vv+1};
    elseif strcmpi(varargin{vv},'startSize_px')
        startSize_px = varargin{vv+1};
    elseif strcmpi(varargin{vv},'soothesayer')
        if strcmpi(varargin{vv+1},'on')
            permitSoothesaying = true;
        elseif strcmpi(varargin{vv+1},'off')
            permitSoothesaying = false;
        else
            warning('unexpected parameter for soothesayer property')
        end
    elseif strcmpi(varargin{vv}, 'soothesayerLimit')
        soothesaying_limit = varargin{vv+1};
    elseif strcmpi(varargin{vv},'bufferWidthMultiplier')
        bufferWidthMultiplier = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'grainOptimizerLimits')
        grainOptimizerLimits = varargin{vv+1};
    end
end; clear vv varargin

%save these in output log
simParams.simVersion = simVersion;
simParams.saveIm = saveIm;
simParams.pxLimit = pxLimit;
simParams.growthRate_px = growthRate_px;
simParams.startSize_px = startSize_px;
simParams.soothesayer = permitSoothesaying;
simParams.soothesayerLimit = soothesaying_limit;
simParams.bufferWidthMultiplier = bufferWidthMultiplier;
simParams.grainOptimizerLimits = grainOptimizerLimits;

clear simVersion


%% Initialize some colors
color_grainActive = 0.8;
color_bufferRegion = 0.6;
color_fieldOfView = 0.2;
color_grainFinished = 0.5;

%store them in output log
simParams.colors.color_grainActive = color_grainActive;
simParams.colors.color_bufferRegion = color_bufferRegion;
simParams.colors.color_fieldOfView = color_fieldOfView;
simParams.colors.color_grainFinished = color_grainFinished;


%% create a matrix of the size of the simulation area. 
%Can't go too large or you'll run out of memory. don't want to be too small or the results get thrown off.
%To work around this, we scale each pixel to some value depending on the
%nucleation density and number of pixels desired in the simulation.

numGrainsROI = grainOptimizerLimits(1):5:grainOptimizerLimits(2); %helper for determining how many grains to include. doesn't affect final nucleation density
scaling_exp = 0:6; %0 is mm pixels, 6 is nm pixels
pxSize = {'1mm','100um','10um','1um','100nm','10nm','1nm'}; %pixel size (10.^scaling_exp) is in units of px/mm

%optimize the scaling factor, the number of grains, and total number of pixels.
area_optimizer_numPixels = zeros(length(numGrainsROI), length(pxSize)); %init
for ii = 1:length(numGrainsROI) %get the size to be manageable
    area_fieldOfView_mm2 = numGrainsROI(ii)/nuc_density_mm2; %in mm^2, determined by nucleation density and number of grains to put in field of view
    area_fieldOfView_px2 = area_fieldOfView_mm2.*((10.^scaling_exp).^2); %in px^2, for processing as pixels
    width_fieldOfView_px1 = round(sqrt(area_fieldOfView_px2)); %in px (length of one side)
    width_total_px1 = round(width_fieldOfView_px1.*(bufferWidthMultiplier*2 + 1)); %give buffer region (pseudo periodic boundary)
    area_optimizer_numPixels(ii,:) = width_total_px1;
end; clear ii

%choose the best values
foundValue = false;
for bb = length(pxSize):-1:1 %give preference for higher resolution pixels
    for aa = 1:length(numGrainsROI) %preference for smaller number of grains
        if area_optimizer_numPixels(aa, bb) >= pxLimit(1) && area_optimizer_numPixels(aa, bb) <= pxLimit(2)
            foundValue = true;
            numGrainsROI = numGrainsROI(aa);
            scaling_exp = scaling_exp(bb);
           break 
        end
    end; clear aa
    if foundValue
        break
    end
end; clear bb

%if your pixel limits are too narrow
if ~foundValue
    for bb = length(pxSize):-1:1 %give preference for higher resolution pixels
        for aa = length(numGrainsROI):-1:1 %in this case, prefer larger number of grains
            if area_optimizer_numPixels(aa, bb) <= pxLimit(2)*1.1
                foundValue = true;
                numGrainsROI = numGrainsROI(aa);
                scaling_exp = scaling_exp(bb);
               break 
            end
        end; clear aa
        if foundValue
            break
        end
    end; clear bb
end; clear foundValue
simParams.scaling_exp = scaling_exp; %store in log

%calculate the values to use later
area_fieldOfView_mm2 = numGrainsROI/nuc_density_mm2; %in mm^2, determined by nucleation density and number of grains to put in field of view
area_fieldOfView_px2 = area_fieldOfView_mm2*((10^scaling_exp)^2); %in px^2, for processing as pixels
width_fieldOfView_px1 = round(sqrt(area_fieldOfView_px2)); %in px (length of one side)
width_total_px1 = round(width_fieldOfView_px1.*(bufferWidthMultiplier*2 + 1)); %give buffer region (pseudo periodic boundary)
area_total_px2 = round(width_total_px1^2);
area_total_mm2 = area_total_px2 / (10^scaling_exp)^2; %px^2 / (px^2/mm^2)

%store in log
simParams.widthROI_px1 = width_fieldOfView_px1;
simParams.widthTotal_px1 = width_total_px1;


%calculate number of grains in the entire region
num_grains_total = round(area_total_mm2 * nuc_density_mm2); %in numGrains

clear area_fieldOfView_mm2 area_fieldOfView_px2 area_optimizer_numPixels area_total_mm2 area_total_px2 bufferWidthMultiplier

%Create an array for the grain coordinates
[x, y] = meshgrid(1:width_total_px1, 1:width_total_px1);
edge_LB = width_total_px1/2-width_fieldOfView_px1/2;
edge_RT = width_total_px1/2+width_fieldOfView_px1/2;
    
fprintf('%i nucleation points chosen with pixel size %s.\n',numGrainsROI,pxSize{scaling_exp})

%% Output an image of the starting positions
imtitlebase = sprintf('NucDens%0.4f_Grains%0.0f_GrowthRate%0.0fx%s',nuc_density_mm2,numGrainsROI,growthRate_px,pxSize{scaling_exp});
imtitle = sprintf('%s_0_init.png',imtitlebase);
doesfileexist = eq(exist(imtitle,'file'),2);
filenum = 0; tempimtitle = imtitle;
while doesfileexist
    filenum = filenum+1;
    tempimtitle = sprintf('%s_%i_init.png',imtitlebase,filenum);
    doesfileexist = eq(exist(tempimtitle,'file'),2);
end
imtitlebase = sprintf('%s_%i',imtitlebase,filenum);

clear tempimtitle filenum doesfileexist imtitle 


%% initialize struct with beginning data
clear grains
grains(num_grains_total) = struct(); %initialize struct
for gg = 1:length(grains)
    grains(gg).grow = true; %is it currently growing-enabled?
    %we'll get their positions later
%     grains(gg).posX = rand()*system_size_um1; %in um
%     grains(gg).posY = rand()*system_size_um1; %in um
%     grains(gg).inROI = False;
    grains(gg).rot = rand()*360; %in degrees
    grains(gg).size = startSize_px; % in px (side length)
    grains(gg).waitUntil = 1; %iteration number for acceleration (soothesayer) feature
end; clear gg


%% initialize beginning nucleation positions

%got some help here: https://www.mathworks.com/matlabcentral/answers/329166-plotting-random-points-within-boundary

tic %counter for whole process

%initialize an image with color of the buffer region
ourImage = color_bufferRegion * ones(width_total_px1, width_total_px1);
%color the field of view it's proper color
ourImage(x >= (width_total_px1/2-width_fieldOfView_px1/2) & x <= (width_total_px1/2+width_fieldOfView_px1/2) & y >= (width_total_px1/2-width_fieldOfView_px1/2) & y <= (width_total_px1/2+width_fieldOfView_px1/2)) = color_fieldOfView;


% Place the small squares inside one by one.
numG = 1; %grain identifier
numBytes = 0; %for printing progress
failsafe = 0; %to prevent eternal loops
while numG <= length(grains) && failsafe < length(grains) * 100
    failsafe = failsafe + 1;
  
    % Get a random row, column location for nucleation
    x_cent = randi(width_total_px1, 1);
    y_cent = randi(width_total_px1, 1);
    
	% See if any of the centers are outside the system. If so, mark them as
	% being outside the region of interest
    inROI = true;
    if x_cent <= edge_LB || x_cent >= edge_RT
        inROI = false;
    end
    
    if y_cent <= edge_LB || y_cent >= edge_RT
        inROI = false;
    end
  

    %calculate their corner coordinates
    y1 = y_cent - round(grains(numG).size/2);
    y2 = y_cent + round(grains(numG).size/2);
    x1 = x_cent - round(grains(numG).size/2);
    x2 = x_cent + round(grains(numG).size/2);
 
  
    if rem(grains(numG).rot,90)==0 %if there is no rotation
        grainLogical = y >= y1 & y <= y2 & x >= x1 & x <= x2;
        
    else %most cases there will be some rotation
        vertices = [x1, x1, x2, x2; y1, y2, y2, y1];
        center = repmat([x_cent y_cent],4,1)';
        angle = grains(numG).rot;

        [unique_slopes,yint1,yint2] = prepSquareRot(vertices,angle,center);
        
        temp1 = unique_slopes(1).*x;
        temp2 = unique_slopes(2).*x;
        grainLogical = y >= temp1 + min(yint1) & y <= temp1 + max(yint1) & y >= temp2 + min(yint2) & y <= temp2 + max(yint2);
        clear temp1 temp2
    end
    
    %now make sure none of this grain intersects with another that has already been placed
    %if it does, stop and try a new place.
    if ~isempty(find(ourImage(grainLogical) == color_grainActive, 1))
        continue %goes to the next loop iteration
    end
  
  
    %if we get here, we're free to place the square.
    ourImage(grainLogical) = color_grainActive;

    %accept the position in our grain struct
    grains(numG).posX = x_cent;
    grains(numG).posY = y_cent;
    grains(numG).inROI = inROI;

    if rem(numG,5)==0 %so we don't update as often. Can get kinda crazy.
        fprintf(repmat('\b',1,numBytes)) %delete last line from command window
        numBytes = fprintf('initializing grain %i/%i.\n',numG,length(grains));
    end
    

    
    numG = numG + 1; %increment

end


clear grainLogical angle center inROI unique_slopes vertices x1 x2 y1 y2 y_cent x_cent yint2 yint1


imwrite(ourImage,sprintf('%s_init.png',imtitlebase)) %write an initial image
if saveIm
    imwrite(ourImage,sprintf('%s_step%04d.png',imtitlebase,0)) %save first image
end



%% grow the grains


%initialize some while loop stuff
timestep = 0;
nextWaitUntil = 1;
growingGrains = find([grains.grow] == 1); %list of all the grains that are still active
numBytes = 0;
numBytesSoothesayer = 1;
numBytesExpand = 2;

%init soothesayer end point. Most are done before ~2000. But what if they
%don't? Soothesayer will adjust below if it makes it beyond this point.
%This is just to give an initial end goal to give good half-points.
%Example: Starts at 2000, then tests 1000, then 500, then 250... I could
%find the definite end point by going grain by grain and finding the max
%number of iterations it takes for every grain to collide with its
%neighbors in their embryo form, but I think that would be slower than this
%method and yield the same result.
soothesayerEnd = 2000; 


while ~isempty(growingGrains) %while grains are still growing
    timestep = timestep + 1;
    
    new_grainSize = startSize_px + growthRate_px * timestep; %new proposed grain size
    
    
    %acceleration feature. See if we can not redraw certain grains for a
    %while. Basically it looks X steps in the future and if a grain
    %has not collided with other grains, it tells the system not to pay
    %attention to that grain for the next several steps.
    if timestep >= nextWaitUntil && permitSoothesaying %if this timestep will contain a soothesayer event
        numBytes = fprintf('Soothesayer active on timestep %i\n', timestep); %display notification
        
        if timestep > soothesayerEnd %because of the note above. Sometimes the end point is not always accurate.
            numBytes = fprintf('Increasing soothesayer end point from %i to %i\n', soothesayerEnd, soothesayerEnd + 500);
            soothesayerEnd = soothesayerEnd + 500;
        end
        
        nextWaitUntil = min([min([grains([grains.waitUntil] > timestep).waitUntil]) soothesayerEnd]); %next time to visit the soothesayer.
        deltaStep = nextWaitUntil - timestep; %If the timestep is close to the next soothesayer event, just go ahead and grow the grains. Likely there are only a few grains to grow anyway.
        while deltaStep > soothesaying_limit %this value can be set by the user
            soothesayerStep = ceil(mean([nextWaitUntil timestep])); %find the new timestep to explore
            soothesayerSize = startSize_px + growthRate_px * soothesayerStep; %what is the grain size at this timestep?
            candidateGrains = find([grains.grow] == 1 & [grains.waitUntil] <= soothesayerStep); %all active grains that may collide before this step
            
            if ~isempty(candidateGrains)
                if eq(numBytesSoothesayer, numBytes)
                    fprintf(repmat('\b',1,numBytes)) %delete last line from command window
                end
                numBytes = fprintf('\tenvisioning %i\n', soothesayerStep);
                numBytesSoothesayer = numBytes;
                soothesayerSays = soothesayer(grains, candidateGrains, soothesayerStep, soothesayerSize, x, y, ourImage, color_grainActive, color_grainFinished); %look into the future. reports identity of everyone who has not yet collided
            end
            
            newWaitUntil = num2cell(max([grains.waitUntil], soothesayerSays)); %update waitUntil array
            [grains.waitUntil] = newWaitUntil{:};            
            
            %sometimes the soothesayer lands on a timestep where there are
            %no survivors (but locally, not globally). Prevent eternal
            %loops by making sure the soothesayer knows we have tried this
            %timestep before
            if isempty(find([grains.waitUntil] == soothesayerStep, 1))
                nextWaitUntil = soothesayerStep;
            else
                nextWaitUntil = min([min([grains([grains.waitUntil] > timestep).waitUntil]) soothesayerEnd]);
            end


            deltaStep = nextWaitUntil - timestep; %calculate new deltaStep

        end %while loop for this grain number
        
        
%         fprintf('\tfinished soothesaying\n')
    end %soothesayer
    
    
        
    
    %Cool. Now we have a much shorter list of grains that are actually
    %growing on this particular timestep. Most of them won't collide until
    %later on, so we don't need to draw them. Some of them have already
    %collided, and we already have them burned into our image.
    
    growingGrainsNOW = find([grains.waitUntil] == timestep & [grains.grow] == 1);
    

    for gg = growingGrainsNOW %go through all the activeNOW grains in the struct
        
        grains(gg).size = new_grainSize; %in px, grow
        
        %determine which part of the image to change.
        grainLogical = getGrainLogical(grains(gg), x, y);
        
        %Use the logical to change the color of the newly grown grain.
        %Not sure yet if there has been a collision
        ourImage(grainLogical) = color_grainActive; 
      
    end %done going through all active grains
    
    if ~isempty(growingGrainsNOW)
        if eq(numBytes,numBytesExpand)
            fprintf(repmat('\b',1,numBytes)) %delete last line from command window
        end
        numBytes = fprintf('Expanded %i activeNOW. frame %i. ',length(growingGrainsNOW),timestep);
        if permitSoothesaying
            numBytes = numBytes + fprintf('next soothesayer %i. ',nextWaitUntil);
        end
        numBytes = numBytes + fprintf('total remaining: %i\n', length(growingGrains));
        numBytesExpand = numBytes;
    end

    if ~isempty(growingGrainsNOW)
        %now determine if there were grain collisions
        allGrains = ourImage == color_grainActive | ourImage == color_grainFinished;
        L = bwlabeln(allGrains, 8); %label spots (all grains)
        S = regionprops(L, 'Area'); %measure area of spots
        listofspotnumbersbygrain = zeros(length(grains),1);
        for ss = 1:length(S) %go through all spots
            for gg = 1:length(grains) %need to go through all grains
                listofspotnumbersbygrain(gg) = L(grains(gg).posY,grains(gg).posX);
            end; clear gg
        end; clear ss
    end
   
    for gg = growingGrainsNOW
        if length(find(listofspotnumbersbygrain == listofspotnumbersbygrain(gg))) > 1
            ourImage(L == listofspotnumbersbygrain(gg)) = color_grainFinished; %color the connected grains dormant colored.
            grains(gg).grow = false;
        else
            grains(gg).waitUntil = timestep + 1;
        end
    end; clear gg    
    
    
    if saveIm && ~isempty(growingGrainsNOW)
        imwrite(ourImage,sprintf('%s_step%04d.png',imtitlebase,timestep))
    end
    
    growingGrains = find([grains.grow] == 1); %update growingGrains
    
end %after there are no more actively growing grains

timetook = toc; 
imwrite(ourImage,sprintf('%s_final.png',imtitlebase)) %write a final image

%save in log
summary.timetook = timetook; 
summary.iterations = timestep;
simParams.imTitleBase = imtitlebase;


%summarize some of the data for convenience
gsizeROI_px = [grains([grains.inROI]==1).size]; %grain sizes of grains in ROI (length of one side)
gsizeROI_um = gsizeROI_px./10^scaling_exp.*1000; %put back in mm, then to um
fprintf('The average grain size is %0.3f um (after %0.1f min)\n',mean(gsizeROI_um),timetook/60)

summary.size_um = gsizeROI_um; %grain sizes of just those in the ROI (um)
summary.areaROI_mm2 = (simParams.widthROI_px1/10^scaling_exp)^2; %mm^2
summary.areaTotal_mm2 = (simParams.widthTotal_px1/10^scaling_exp)^2; %mm^2
summary.nucDensROI = length(gsizeROI_um) / summary.areaROI_mm2; %#/mm^2
summary.nucDensTotal = length(grains) / summary.areaTotal_mm2; %#/mm^2



end