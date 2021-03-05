# crystallizationSimulation
A MATLAB simulation of crystal nucleation and growth

Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)

Created in MATLAB R2016a (though the more convenient plotting portion of it needs 2020b)


csimmath_portal.m is provided as a convenience for performing the simulation with the helper functions. The results are saved after each new simulation to avoid sudden computer-restart catastrophe.

The main function, csimmath_main.m, performs a simulation for a single inputted temperature. The temperature is then used with hard-coded exponential relationships (contained in the structure variable called 'key') that calculate the simulation nucleation rate and growth velocity based on experimental values, as well as predict approximately how large the simulation should be to provide a manageable number of crystals for computation. The crystals in the simulation zone are nucleated at time steps dictated by the nucleation rate with the assistance of csimmath_placeCrystal.m. Nuclei cannot be placed inside any existing crystal as checked by csimmath_checkEmbed.m or at a location intersecting any existing crystal as checked by csimmath_checkIntersect.m, csimmath_checkIntersectOne.m, csimmath_linesIntersect.m, and csimmath_isBetween.m. The crystals are then grown at a constant given velocity with csimmath_growCrystals.m and continually checked for intersection. Once no suitable locations can be found for nucleation, the simulation is terminated and a structure variable is output with all simulation parameters and final crystal descriptions. An image is generated and saved with a shortcut function savePlots.m.
