Here we want to compare results for the population based and the individual based simulators.
I added "simulation.cpp" and "two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R" because you need to run these, before playing around with plottings.
Population based simulator looks to have an upward bias compared to the true individual simulator.
This difference in the plottings is very clear if you start with the two populations in equilibrium (p and q equal to c(0.25, 0.25, 0.25, 0.25) ).
That's simply becuase the range of the y-axis is much narrower if you start with r^2 = 0.
