In this section, I try to see how well the "Multinomial/Normal Test" responds to migration signals. And how much difference there is for tests run on each of the 8 frequencies.

To this purpose, I test a bunch of cases with different migration rates, while all of the other parameters are fixed.

Other params: N = 1000, M = 1000, generations = 100, c = 0.1, initial frequencies in equilibrium c(0.25,0.25,0.25,0.25)

Migration rates are c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

For each of these rates I do 1000 simulations. For each simulation I simulate the 100 generations using the C++ script, then I consider the frequencies p and q at generation 100, and make a LR test on these final frequencies. I store all the 1000 r^2 values (they are the p-values of the LR tests), and then I plot the empirical probability distribution function for each migration rate.

Then I compare this distribution to the case where the is no migration (m = 0.001) to see how well-separated the two are.

These are done in "single_plots_started_in_equilibrium", "comparison_migration_no_migration".

Even though there is a clear trend, I think the other approach (chi-sq contingency table) gives much better separabiity.



"disequilbrium, no migration, 8 tests plots" shows the tests run on the 8 different frequencies, when run on two populations with no migration, and very different initial frequencies.

"equilibrium, no migration, 8 test plots" does the same thing with two pops started in equilibrium with no migration

"equilibrium, full migration, 8 test plots" does the same thing with two pops started in equilibrium with 0.5 migration

For the last 2 symmetric scenarios, the 8 tests give essentially the same results. For the asymmetric case, the differences are notable.

