In this section, I try to see how well the Likelihood-Ratio Test Statistic developed in the subfolder 'theory and proposals', picks up migration signals.

To this purpose, I test a bunch of cases with different migration rates, while all of the other parameters are fixed.

Other params: N = 1000, M = 1000, generations = 100, c = 0.1, initial frequencies in equilibrium c(0.25,0.25,0.25,0.25)

Migration rates are c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

For each of these rates I do 1000 simulations. For each simulation I simulate the 100 generations using the C++ script, then I consider the frequencies p and q at generation 100, 
and make a LR test on these final frequencies. I store all the 1000 r^2 values (they are the p-values of the LR tests), and then I plot the empirical probability distribution
function for each migration rate.

Then I compare this distribution to the case where the is no migration (m = 0.001) to see how well-separated the two are.

High values of my 'r^2' definitiely confirm that there is non-existent/extremely-low migration between the two pop.
Low values of 'r^2' give more uncertain answers: Why? Think about two independent populations started at equilibrium c(0.25,0.25,0.25,0.25). The final frequencies are generally 0.25 +- 0.3
If p_{00} = 0.28 and q_{00} = 0.23, you are sure that there is no migration.
However it can be the case that p_{00} = q_{00} = 0.28, even if the two pop are disjoint.
