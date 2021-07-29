Here we want to compare results for the population based and the individual based simulators.

I added "simulation.cpp", "two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R" and "two_populations_fixed_matrix_variable_N_pop_based.R" because you need to run these, before playing around with plottings.

I generally expected population based simulators to have a downward bias compared to the true individual simulator (as it happened with 1-population simulators).

This was the case with "two_populations_fixed_matrix_variable_N_pop_based.R" (reproduction first, migration second), as the algorithm followed the same structure of the 1-pop population based simulator (i.e. (1) non-stochastic recombination  (2) multinomial sampling  (3) calculating r^2 for generation).

The very interesting result was with "two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R". Here the algorithm had a slightly more complex structure

This difference in the plottings is very clear if you start with the two populations in equilibrium (p and q equal to c(0.25, 0.25, 0.25, 0.25) ).
That's simply becuase the range of the y-axis is much narrower if you start with r^2 = 0.
