Population based simulators for 2 migrating populations.

"two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R" is the one I have used in later analyses.
It has the right sequence of migration (first) and reproduction (second) at each generation. Moreover, we use multinomial sampling twice for each iteration, increasing the variance (and hence the r^2), making in an UNBIASED estimator

The other two scripts have reproduction first and then migration. When migration is set to zero, they are they same algorithm as for 1-population population based simulator.
