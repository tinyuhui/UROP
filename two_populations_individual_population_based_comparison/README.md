Here we want to compare results for the population based and the individual based simulators.

I added "simulation.cpp", "two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R" and "two_populations_fixed_matrix_variable_N_pop_based.R" because you need to run these, before playing around with plottings.

I generally expected population based simulators to have a downward bias compared to the true individual simulator (as it happened with 1-population simulators).

This was the case with "two_populations_fixed_matrix_variable_N_pop_based.R" (reproduction first, migration second), as the algorithm followed the same structure of the 1-pop population based simulator:

(1) non-stochastic recombination (using fixed proportion c)  (2) multinomial sampling for reproduction  (3) calculation of r^2 for generation 

The very interesting result was with "two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R". Here the algorithm had a slightly more complex structure:

(1) multinomial sampling for migration  (2) non-stochastic recombination (using fixed proportion c)  (3) multinomial sampling for reproduction  (4) calculation of r^2 for generation 

 In this case, we use multinomial sampling twice for each iteration, increasing the variance (and hence the r^2), making it an upward biased estimator (for the individual based r^2).
 
 The plottings for this second algorithm, make a solid argument in favour of this hypotesis. 
 
I also analysed the case, where population 2 is 10 times bigger than population 1. The plots behaved as expected, showing that smaller population sizes have bigger fluctuations in the frequencies around the mean (think SLLN), and that migration partially mitigates this variance.

The difference in the r^2 equilibrium values is very clear if you start with the two populations in equilibrium (p and q equal to c(0.25, 0.25, 0.25, 0.25) ).
That's simply becuase the range of the y-axis is much narrower if you start with r^2 = 0. Therefore all the plots I made, had populations starting in equilibrium.
