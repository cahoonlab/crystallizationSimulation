%% The key to temperature/nucleation density conversion

%% first load fit_nucDens.mat
func_NucDens2TinvK = @(ND) (log(ND) - fit_nucDens.intercept) ./ fit_nucDens.slope;
func_TinvK2NucDens = @(TinvK) exp(fit_nucDens.slope .* TinvK + fit_nucDens.intercept); %#/mm2


%% calculate a grain size with the simulator

numSimIterations = 1000;

%temperature region of interest
TinvK = 2.4E-3:1E-5:3.6E-3; %inverse Kelvin


clear gsSim %clear any old struct
gsSim(numSimIterations) = struct(); %initialize struct
for ii = 1:numSimIterations
    %get a random temperature in our region of interest
    TinvK_current = rand()*(max(TinvK) - min(TinvK)) + min(TinvK);
    %then get the associated nucleation density
    nucDens = func_TinvK2NucDens(TinvK_current);
    
    %do the simulation
    [simParams, simGrains, simSummary] = ...
        nucGrowthSim_v12(nucDens,...
        'saveIm',false,... %be careful with this. Saves every frame as image.
        'growthRate_px',1,... %1 is the most accurate, but takes longest
        'pxLimit',[5000, 10000],...
        'startSize_px',1,... %1 is the most accurate (maybe?), but takes longer
        'soothesayer','on',... %speeds things up a lot!
        'soothesayerLimit',15,...
        'bufferWidthMultiplier',.4,...
        'grainOptimizerLimits',[50, 100]); 
    
    %add results to the struct
    gsSim(ii).nuc_dens = simSummary.nucDensROI; %just that in ROI
    gsSim(ii).size_um = simSummary.size_um; %just those in ROI
    gsSim(ii).simParams = simParams;
    gsSim(ii).simGrains = simGrains;
    gsSim(ii).simSummary = simSummary;
    fprintf('--------------------------\n-------------------------\n')
    fprintf('the last simulation took %0.1 minutes\n',simSummary.timetook/60)
    fprintf('--------------------------\n-------------------------\n')
    save(sprintf('gsSimProg%i',ii),'gsSim')
end; clear rr

