function [ simSummary ] = csimmath_main( simTemp_invK, varargin )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% The main function for crystalization simulator ("mathematician")
%%% Only input it needs is the temperature in inverse Kelvin. However, you
%%% can also input other optional variable pairs including:

% 'grainNum', 100  (About how many grains you want in the simulation)
% 'calcZone', 0.75 (Ratio of the calculation region to simulation region)
% 'NRslopeChange', 0 (Explore other nucleation rate dependencies)
% 'timeMult', 1    (Determines frame rate. 0.5 would be 2 per s)
% 'debugPlot', false (If true, will generate figure per time step)
% 'recordHistory', false (If true, save figures & data every 10 steps)
% 'simNum', 0      (for identifying purposes)
% 'embryoSize', 5e-6 (How large should embryo be in cm)
% 'placeTryLimit', 1500 (How many times to try placing embryo)


%%% ---Version Notes---
%%%v2 - deactivated crystals still "grow" so that little grains don't grow
%%%in the space that they would have expanded to in a physical experiment.
%%%However, the final grain size is calculated based on the size that the
%%%crystals were when they were deactivated.

%%%v3 - instituted a fixed embryo size for any temperature

%%%v4 - removed ability for the simulation to stop based on nucleation
%%%density. Now it only stops if embryos can't find a nucleation location.


%%% ---Dependency Notes---
% csimmath_checkEmbed
% csimmath_checkIntersect
% csimmath_checkIntersectOne
% csimmath_growCrystals
% csimmath_isBetween
% csimmath_linesIntersect
% csimmath_placeCrystal
% savePlots


% The key to calculating growth rates, nucleation rates/density
% These are from fits of experimental data.
% In parenthesis is the error (+-)
key = struct();
key.growthRateSlope = -4.883454786150889e+03; %um/s (4.116065627696693e+02)
key.growthRateIntercept = 16.906828507291355; %um/s (1.152522538849132)
key.nucRateSlope = 1.723812577524499e+04; %#/cm^2 (2.040625805787124e+03)
key.nucRateIntercept = -44.811317244828640; %#/cm^2 (5.698342801669105)
key.nucDensSlope = 1.397499834192973e+04; %#/cm^2 (1.811070665687012e+03)
key.nucDensIntercept = -32.223663979782350; %#/cm^2 (5.016732676739276)


%calculate the nucleation/growth parameters
growthRate = exp(...
    key.growthRateSlope * simTemp_invK + key.growthRateIntercept...
    ) * 10^-4; % cm/s
nucleationRate = exp(...
    key.nucRateSlope * simTemp_invK + key.nucRateIntercept...
    ); % #/cm^2/s
nucleationDensity = exp(...
    key.nucDensSlope * simTemp_invK + key.nucDensIntercept...
    ); % #/cm^2


grainNumTarget = 100; %default (just used to establish size of zones)
fracZoneCalc = 0.75; %default
for vv = 1:length(varargin)
    if strcmpi(varargin{vv}, 'grainNum')
        grainNumTarget = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'calcZone')
        fracZoneCalc = varargin{vv+1};
        if fracZoneCalc > 1
            fracZoneCalc = fracZoneCalc/100;
        end
    elseif strcmpi(varargin{vv}, 'NRslopeChange')
        key.nucRateSlope = key.nucRateSlope + varargin{vv+1};
    end
end; clear vv


%determine size of simulation
simZoneSize_area = grainNumTarget / nucleationDensity; %cm^2
simZoneSize_width = sqrt(simZoneSize_area); %cm
calcZoneSize_area = fracZoneCalc * simZoneSize_area; %cm^2
calcZoneSize_width = sqrt(calcZoneSize_area); %cm


%simulation timing
nucleationFrequency = 1 / (nucleationRate * simZoneSize_area); % s/1 nucleation
timeMultiplier = 1; %s

%fixed embryo size
embryoSize = 5e-6; %50 nm (in cm)


debugPlot = false; %default
recordHistory = false; %default
simNum = 0; %default
placeTryNumLimit = 1500; %default
for vv = 1:length(varargin)
    if strcmpi(varargin{vv}, 'timeMult')
        timeMultiplier = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'debugPlot')
        debugPlot = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'recordHistory')
        recordHistory = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'simNum')
        simNum = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'embryoSize')
        embryoSize = varargin{vv+1};
    elseif strcmpi(varargin{vv}, 'placeTryLimit')
        placeTryNumLimit = varargin{vv+1};
    end
end; clear vv


%organize data structure to make it easier to pass data in and through
simSummary = struct();
simSummary.TinvK = simTemp_invK;
simSummary.key = key;
simSummary.GR = growthRate;
simSummary.NR = nucleationRate;
simSummary.ND = nucleationDensity;
simSummary.grainNumTarget = grainNumTarget;
simSummary.simZoneSize_area = simZoneSize_area;
simSummary.simZoneSize_width = simZoneSize_width;
simSummary.calcZoneSize_area = calcZoneSize_area;
simSummary.calcZoneSize_width = calcZoneSize_width;
simSummary.nucFreq = nucleationFrequency;
simSummary.nextNuc = nucleationFrequency;
simSummary.timeMultiplier = timeMultiplier;
simSummary.embryoSize = embryoSize;

clear simTemp_invK growthRate nucleationRate nucleationDensity grainNumTarget simZoneSize_area simZoneSize_width calcZoneSize_area calcZoneSize_width nucleationFrequency timeMultiplier key embryoSize



%% Initialize the first embryo
tic; %start counting time

simSummary.timeStep = 1; %there is no time 0
simSummary.crystalCount = 0;

%generate embryo
embryo = csimmath_placeCrystal(simSummary);
embryo.active = true;

%This crystals struct used to be a history of ALL time steps. It can
%probably be removed for future versions.
clear crystals %initialize struct
crystals = embryo; %add the embryo to the list of crystals
crystals(simSummary.grainNumTarget * 100).id = [];
%(over-initializing to improve performance)



%% Start the growth loop

growing = true; %establish loop controller
while growing
    
    %make a subset of the crystals struct that is only the crystals present
    %on the most recent timestep.
    currentCrystals = crystals(...
        eq([crystals.timeStep], max([crystals.timeStep]))...
        );
    
    %advance time
    simSummary.timeStep = simSummary.timeStep + 1;
    simSummary.crystalCount = size(currentCrystals,2);
    %apply that new time to the current list of crystals
    [currentCrystals.timeStep] = deal(simSummary.timeStep);
    %display that time in the command window (occasionally)
    if rem(simSummary.timeStep,100)==0
        fprintf('%i crystals present (%i active).\n',...
            simSummary.crystalCount, sum([currentCrystals.active]))
    end
    
    %grow the existing crystals
    currentCrystals = csimmath_growCrystals(simSummary, currentCrystals);
    
    %check to see which ones have collided
    collided = csimmath_checkIntersect(currentCrystals);
    if ~isempty(collided) %if some collided, mark them as deactivated
        [~,collidedIdx] = ismember(collided,[currentCrystals.id]);
        for mm = collidedIdx
            if currentCrystals(mm).active
                currentCrystals(mm).active = false;
                %(new in v2) record final size and vertices
                currentCrystals(mm).L_final = currentCrystals(mm).L;
                currentCrystals(mm).xcoords_final = currentCrystals(mm).xcoords;
                currentCrystals(mm).ycoords_final = currentCrystals(mm).ycoords;
            end
        end; clear mm
    end
    
    
    
    %if it's time for a nucleation event, place a new crystal
    if simSummary.timeStep * simSummary.timeMultiplier >= simSummary.nextNuc
        
        %update next nucleation flag
        simSummary.nextNuc = simSummary.nextNuc + simSummary.nucFreq; 
        
        %Start a loop. Need to make sure the new crystal isn't already
        %colliding with or embedded within another crystal
        crystalPlaced = false; %initialize
        placeTryNum = 0;
        while ~crystalPlaced
            %generate a candidate embryo
            embryo = csimmath_placeCrystal(simSummary);
            %is it embedded within an existing crystal?
            embedded = csimmath_checkEmbed(currentCrystals, embryo);
            %is it intersecting with other crystals?
            intersected = csimmath_checkIntersectOne(currentCrystals, embryo);

            if ~embedded && isempty(intersected)
                %success! The candidate embryo is born. Add to list.
                embryo.placeAttempts = placeTryNum;
                embryo.active = true;
                currentCrystals(size(currentCrystals,2)+1) = embryo;
                crystalPlaced = true;
            end
            
            placeTryNum = placeTryNum + 1; %increment
            
            %determine if we should stop trying
            %in other words, is the simulation zone full?
            if placeTryNum > placeTryNumLimit
                crystalPlaced = true;
                growing = false;
                simSummary.endCause = 'placeTryNumLimit';
            end

        end %waiting for nucleus to be placed correctly

    end %time to add nucleus

    
    %are all the crystals deactivated?
    %and are we far from the next nucleation event?
    if sum([currentCrystals.active]) == 0 &&...
            simSummary.nextNuc / simSummary.timeMultiplier - simSummary.timeStep > 10
        
        %If so, skip to next nucleation event (time saver)
        simSummary.timeStep = ...
            round(simSummary.nextNuc / simSummary.timeMultiplier)-1;
    end


    %Just informational
    simSummary.NDnow = size(currentCrystals,2) / simSummary.simZoneSize_area;
    
    
    %I'm not sure I need to write this every time, but it works
    for ii = find([currentCrystals.active])
        currentCrystals(ii).L_final = currentCrystals(ii).L;
        currentCrystals(ii).xcoords_final = currentCrystals(ii).xcoords;
        currentCrystals(ii).ycoords_final = currentCrystals(ii).ycoords;
    end; clear ii
    
    
    if recordHistory
        %Add the currentCrystals to the crystals history (old from v3)
%         crystals(...
%             length([crystals.id])+1 :...
%             length([crystals.id]) + size(currentCrystals,2)...
%             ) = currentCrystals;
        
        %new in v3 (because I realized file sizes were getting HUGE!)
        crystals = currentCrystals;
        if rem(simSummary.timeStep,10) == 0
            
            %Plot the grains layout
            figure; hold on
            plot(polyshape([-simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2], [simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2]))
            for pp = 1:size([crystals.id],2)
                plot(polyshape(crystals(pp).xcoords_final, crystals(pp).ycoords_final))
            end; clear pp
            axis equal
            xlim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
            ylim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
            title(sprintf('sim%i (step%i): %0.1f *1000/K (%0.2f K)',simNum, simSummary.timeStep, simSummary.TinvK*1000, 1/simSummary.TinvK))
            savePlots(sprintf('sim%istep%i',simNum,simSummary.timeStep))
            close
            
            %Plot the grain shadows layout
            figure; hold on
            plot(polyshape([-simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2], [simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2]))
            for pp = 1:size([crystals.id],2)
                plot(polyshape(crystals(pp).xcoords, crystals(pp).ycoords))
            end; clear pp
            axis equal
            xlim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
            ylim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
            title(sprintf('sim%i (step%i, ext): %0.1f *1000/K (%0.2f K)',simNum, simSummary.timeStep, simSummary.TinvK*1000, 1/simSummary.TinvK))
            savePlots(sprintf('sim%istep%i(ext)',simNum,simSummary.timeStep))
            close
            
            %save the current data
            save(sprintf('sim%istep%i',simNum,simSummary.timeStep),'crystals','simSummary')
        end
        
    else
        crystals = currentCrystals;
    end

    
    if debugPlot
        disp(simSummary.timeStep +1)
        plotCrystals = crystals( eq([crystals.timeStep], max([crystals.timeStep])) );
        
        figure; hold on
        plot(polyshape([-simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2], [simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2]))
        for ii = 1:size([plotCrystals.id],2)
            plot(polyshape(plotCrystals(ii).xcoords, plotCrystals(ii).ycoords))
        end; clear ii
        axis equal
        xlim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
        ylim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
    end
        
end; clear growing


%update the simulation summary
simSummary.timeTook = toc;
simSummary.crystals = crystals([crystals.id]>0); %get rid of empties

%which crystals are in the calculation zone? Create boolean
simSummary.calcZoneBool =...
    eq([crystals.timeStep], max([crystals.timeStep])) &...
    abs([crystals.Cx]) <= simSummary.calcZoneSize_width/2 &...
    abs([crystals.Cy]) <= simSummary.calcZoneSize_width/2;

%summarize the calculation region statistics
calcCrystals = crystals(simSummary.calcZoneBool);
simSummary.calcZoneND = size(calcCrystals,2)/simSummary.calcZoneSize_area;
simSummary.calcZoneL = [calcCrystals.L_final];
simSummary.calcZoneLshad = [calcCrystals.L];

%perform a correction for the void space (only in calculation zone):
%(A_shadow - A_crystal) / sumOfAll(A_shadow - A_crystal) * A_void
AcalcCryst = simSummary.calcZoneL.^2;
AcalcVoid = simSummary.calcZoneSize_area - sum(AcalcCryst);
Ashad_min_Acrys = simSummary.calcZoneLshad.^2 - AcalcCryst;
Acor = Ashad_min_Acrys ./ sum(Ashad_min_Acrys) .* AcalcVoid;
simSummary.calcZoneLcor = simSummary.calcZoneL + sqrt(Acor);




end

