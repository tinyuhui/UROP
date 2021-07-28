

average_r_sq_population_based = function(p = c(0.25, 0.25, 0.25, 0.25), q = c(0.25, 0.25, 0.25, 0.25), 
                                         generations = 50, N = rep(1000, generations), M = rep(1000, generations),
                                         c = 0.5, migration = matrix(c(1,0,0,1), nrow = 2), simulations = 5000){
  
  #inputs: p = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population 1 to have
  #inputs: q = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population 2 to have
  #generations: length(N)
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i for pop 1
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i for pop 2
  #inputs: c = recombination rate
  #inputs: migration: 2x2 migration matrix; migration[1,2] is in the interval [0,1]
  #inputs: simulations = number of simulations you want to average
  #output: average r_squared values for each generation over the (5000) simulations using POPULATION BASED SIMULATOR sim2.ld
  
  r_sq1 = numeric(generations)  #vector containing the r^2 for t = 1, ..., generations for pop 1
  r_sq2 = numeric(generations)  #vector containing the r^2 for t = 1, ..., generations for pop 2
  for ( i in 1:simulations){
    out = sim2.ld(p, q, N, M, c, migration)
    out1 = unlist(out[5])
    out2 = unlist(out[6])
    r_sq1 = r_sq1 + out1
    r_sq2 = r_sq2 + out2
  }
  return(list(r_sq1/simulations, r_sq2/simulations))  #dividing by the number of simulations to get the average
  
}



comparison = function(p = c(0.25, 0.25, 0.25, 0.25), q = c(0.25, 0.25, 0.25, 0.25), 
                      generations = 50, N = rep(1000, generations), M = rep(1000, generations),
                      c = 0.5, migration = matrix(c(1,0,0,1), nrow = 2, ncol = 2) ,simulations = 5000){
  
  #the function simply groups together the 4 r^2 vectors and returns a (generations)x4 matrix
  #1st columns is the average r^2 for pop based sim, for pop 1
  #2nd column  is the average r^2 for ind based sim, for pop 1
  out_a = average_r_sq_population_based(p, q, generations, N, M, c, migration, simulations)
  out_b = average_r_sq_individual_based(p, q, N, M, c, migration, simulations)
  out1 = unlist(out_a[1])
  out2 = unlist(out_b[1])
  out3 = unlist(out_a[2])
  out4 = unlist(out_b[2])
  print(typeof(out1))
  return( matrix(c(out1, out2, out3, out4), generations, 4))
}





#MAIN FUNCTION YOU NEED TO CALL IF YOU WANT TO PLAY AROUND WITH PLOTTINGS
plottings = function(p = c(0.25, 0.25, 0.25, 0.25), q = c(0.25, 0.25, 0.25, 0.25), 
                     generations = 50, N = rep(1000, generations), M = rep(1000, generations),
                     vec = c(0.1, 0.2, 0.3, 0.4, 0.5), migration = matrix(c(1,0,0,1), nrow = 2, ncol = 2) ,simulations = 5000){
  
  #inputs: p = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population to have
  #        q, same for pop 2
  #inputs: vec = vector of recombination rates you want to plot
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i
  #        M same for pop 2
  #inputs: migration is the 2x2 migration matrix
  #inputs: simulations = number of simulations you want to average
  #output: plots of the 4 r^2 for each value of recombination rate c provided
  
  for (c in vec)
  {
    m = comparison(p ,q, generations, N, M, c, migration, simulations)
    
    plot(m[, 1], xlab = 'Generations', ylab='r^2', ylim = c(0, 1.2*max(m[,1])), pch = 19)
    points(m[,2], col = 'red', pch = 19)
    points(m[,3], col = 'blue', pch = 19)
    points(m[,4], col = 'green', pch = 19)
    legend("topright", c( "pop1 pop based", "pop1 indiv based", 
                          "pop2 pop based", "pop2 indiv based") , 
           col = c('black', 'red', 'blue', 'green') , pch = c(19,19,19,19))
    title(paste("N =", as.character(N[1]) ,",M =", as.character(M[1]), ",c =", as.character(c), 
                ",migration =", as.character(migration[1,2]), ",simulations =", as.character(simulations), sep = " "))
  }
}







## EXAMPLES ##################################################################

Rcpp::sourceCpp("./simulation.cpp")   #add relevant path here
#ALSO RUN two_populations_fixed_matrix_variable_N_pop_based_RIGHT_ORDER.R



p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25,0.25,0.25,0.25)
N = rep(1000, 100)
M = rep(1000, 100)
c = 0.2
migration =  matrix(data = c(0.9,0.1, 0.1, 0.9), nrow = 2)
simulations = 1000

#Checking C++ simulator works
simulation(p0, q0, N, M, c, migration)

#Checking C++ average calculates correctly
average_r_sq_individual_based(p0, q0, N, M, c, migration, simulations)

#Checking R average calculates correctly
average_r_sq_population_based(p0, q0, length(N), N, M, c, migration, simulations)

#Checking comparisons puts results in a matrix with 4 columns
comparison(p0, q0, 100, N, M, c, migration, simulations)

#Making the actual plots
plottings(p0, q0, 100, N, M, c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), migration, simulations)





