
equivalence_checker = function(p1, p2, p3, recombination_rate, iterations, n, simulations){

  r_sq_simulator1 = numeric(iterations)
  r_sq_simulator2 = numeric(iterations)
  r_sq_simulator3 = numeric(iterations)
  p4 = 1 - p1 - p2 - p3
  
  for ( i in 1:simulations){
    out1 = (unlist(sim.ld(p0 = c(p1, p2, p3, p4), N = rep(n, iterations), c = recombination_rate)[3]))[1:iterations + 1]
    r_sq_simulator1 = r_sq_simulator1 + out1
    
    out2 = two_loci_population_simulator(p1, p2, p3, recombination_rate, iterations, n)
    r_sq_simulator2 = r_sq_simulator2 + out2
    
    m = constrained_matrix(p1, p2, p3, n)
    out3 = simulator(m, recombination_rate, iterations)
    r_sq_simulator3 = r_sq_simulator3 + out3
  }
  
  return( matrix(c(r_sq_simulator1, r_sq_simulator2, r_sq_simulator3)/simulations, iterations, 3))
}


m = equivalence_checker(0.25, 0.25, 0.25, recombination_rate = 0.1, iterations = 25, n = 500, simulations = 500)

plot(m[, 1], xlab = 'Generations', ylab='r^2', pch = 19)
points(m[,2], pch = 19, col = 'red')
points(m[,3], pch = 19, col = 'blue')

#We can see that the two population based estimators are consistent with the individual based simulator (column three)


#       [,1]         [,2]         [,3]
# [1,]  0.0009947231 0.0009021366 0.0009290699
# [2,]  0.0017667911 0.0017265993 0.0017321151
# [3,]  0.0022018534 0.0023226320 0.0022356280
# [4,]  0.0027649115 0.0029814429 0.0029679434
# [5,]  0.0031465718 0.0033210725 0.0037726829
# [6,]  0.0037094977 0.0036507369 0.0039989675
# [7,]  0.0040067995 0.0041341858 0.0046609976
# [8,]  0.0040525378 0.0042650001 0.0049770850
# [9,]  0.0042444846 0.0045186646 0.0050417765
# [10,] 0.0043152849 0.0042403982 0.0050063445
# [11,] 0.0046640770 0.0041383387 0.0051133968
# [12,] 0.0046949440 0.0041751373 0.0054101931
# [13,] 0.0046504303 0.0044287186 0.0054863638
# [14,] 0.0049953930 0.0042188969 0.0052963367
# [15,] 0.0052347686 0.0043244153 0.0051648280
# [16,] 0.0053256039 0.0045607063 0.0049697518
# [17,] 0.0052935846 0.0046150118 0.0051203274
# [18,] 0.0053621809 0.0045451172 0.0051162966
# [19,] 0.0052693253 0.0045744306 0.0052470544
# [20,] 0.0054146613 0.0048614588 0.0051883937
# [21,] 0.0051438040 0.0048196572 0.0051671263
# [22,] 0.0051642333 0.0046674185 0.0050886400
# [23,] 0.0052387252 0.0045654868 0.0049388270
# [24,] 0.0051488264 0.0047106257 0.0048762714
# [25,] 0.0049936027 0.0050024842 0.0049333426


m_2 = equivalence_checker(0.25, 0.1, 0.25, recombination_rate = 0.1, iterations = 25, n = 500, simulations = 500)

plot(m_2[, 1], xlab = 'Generations', ylab='r^2', pch = 19)
points(m_2[,2], pch = 19, col = 'red')
points(m_2[,3], pch = 19, col = 'blue')
legend("topright", c( "Population sim. Fran", "Population sim. Tin-Yu", "Individual sim. Fran") , col = c('black', 'red', 'blue') , pch = c(19,19,19))
title("Average r^2 for 500 simulations")