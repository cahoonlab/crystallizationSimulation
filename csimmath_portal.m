%%%%%%  Portal for using the crystallization simulator "Mathematician"
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)

numSimIterations = 500;

%temperature region of interest
Trange_invK = 2.4E-3:1E-5:3.6E-3; %inverse Kelvin

clear cSim %clear any old struct
cSim(numSimIterations) = struct(); %initialize struct
for ii = 1:numSimIterations
    %get a random temperature in our region of interest
    simTemp_invK = rand()*(max(Trange_invK) - min(Trange_invK)) + min(Trange_invK);
    
    
    simSummary = csimmath_main(simTemp_invK,...
        'timeMult', 0.05,...
        'debugPlot', false,...
        'recordHistory', false,...
        'simNum', ii,...
        'grainNum', 500,...
        'calcZone', .50,...
        'embryoSize', 5e-6,...
        'placeTryNumLimit', 1500);

    
    %for sim#1, initialize combined struct
    if ii==1 && numSimIterations>1
        cSim = simSummary;
        cSim(numSimIterations).TinvK = [];
    elseif numSimIterations==1
        cSim = simSummary;
    else
        cSim(ii) = simSummary;
    end
    
    %plot and save
    crystals = simSummary.crystals;
    plotCrystals = crystals( eq([crystals.timeStep], max([crystals.timeStep])) );
    figure; hold on
    plot(polyshape([-simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2], [simSummary.calcZoneSize_width/2, simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2, -simSummary.calcZoneSize_width/2]))
    for pp = 1:size([plotCrystals.id],2)
        plot(polyshape(plotCrystals(pp).xcoords_final, plotCrystals(pp).ycoords_final))
    end; clear pp
    axis equal
    xlim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
    ylim([-simSummary.simZoneSize_width/2 simSummary.simZoneSize_width/2])
    title(sprintf('sim%i: %0.1f *1000/K (%0.2f K)',ii, simSummary.TinvK*1000, 1/simSummary.TinvK))
    savePlots(sprintf('sim%i',ii))
    
    fprintf('--------------------------\n-------------------------\n')
    fprintf('sim #%i took %0.0f s\n',ii, simSummary.timeTook)
    fprintf('--------------------------\n-------------------------\n')
    save(sprintf('cSimProg%i',ii),'cSim')
    close()
end; clear ii



