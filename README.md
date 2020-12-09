# crystallizationSimulation
A MATLAB simulation of crystal nucleation and growth

Created in MATLAB R2016a


crystallizationSimulationPortal_pub.m is provided as a convenience for performing the simulation with the helper functions. First load the MATLAB variable file 'fit_nucDens.mat' that gives you the key to converting temperature to nucleation density (unless you want to simulate something with a different nucleation density).

After choosing a nucleation density, an optimal simulation size is determined in order to keep memory use reasonable. The simulation zone is divided into two areas: a calculation region and a boundary region. The boundary region exists to provide a psuedo-periodic boundary condition.

Crystals are then initialized randomly in the simulation zone assuming instantaneous nucleation and then grow one pixel at a time at uniform growth rates. Once two crystals collide, they change color, become deactivated, and no longer grow. Once all crystals are deactivated, the simulation ends. Only the crystals in the calculation region should be measured.

To improve performance, a soothesayer function was instated. The soothesayer looks far ahead into future iterations to see if any active crystals have collided. It looks gradually closer and closer and gains an idea of approximately when each crystal will collide. Then the simulator moves frame by frame, growing each crystal that is known to collide before the next frame known by the soothesayer. In this way, many frames can be skipped.

Two images are included herein to illustrate the initialized and finalized simulations.

This simulator is conveniently simplified. However, now that I have experimental data on nucleation rates and growth velocity based on temperature, I think I would remove the assumption of instantaneous nucleation and create a more physically-realistic model. However, this would also be much more time intensive. Though perhaps I could do it mathematically instead of pixel-based. Hmm...
