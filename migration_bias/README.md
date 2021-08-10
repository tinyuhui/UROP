The main function takes fixed population sizes N_val, M_val, simulates 100 generations of them for "simulation" times, takes the average of the last 50 generations (when we are at equilibrium already) 
Now we have equilibrium values of r^2 for each (c, migration) value specified.

We use Park (2) to estimate "N" from the value of the "r^2" and "c"
Then we see how much it is biased upwards by computing the ratio of the estimated N, with the real N

Plots are given for N = 1000, M = 1000; and for N = 1000, M = 10000 (those with the suffix "_10x")

For example:
we can see that as migration grows to say 0.3, for "normal values of c", when tend to overestimate the value of "N" by 30% approximately

The heatmap is definitely easier to read than the 3d plots
