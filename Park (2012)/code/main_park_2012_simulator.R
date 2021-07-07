#You need constrained_matrix.R and LD_pop_TY.R for the script to run, plus modify the path of sim_ld.dll in line 4

# Importing the compiled C code
dyn.load('C:\\Users\\angus\\OneDrive\\Desktop\\urop\\Park (2012)\\sim_ld.dll')

average_r_sq_population_based = function(p = c(0.25, 0.25, 0.25, 0.25), c , generations, N = rep(1000, generations), simulations = 5000){
  
  #inputs: p = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population to have
  #inputs: c = recombination rate
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i
  #inputs: simulations = number of simulations you want to average
  #output: average r_squared values for each generation over the (5000) simulations
  
  r_sq_simulator1 = numeric(generations + 1)  #vector containing the r^2 for t = 0, ..., generations
  
  for ( i in 1:simulations){
    out = (unlist(sim.ld(p0 = p, N = N, c = c)[3])) #the output of sim.ld is a list, with 3rd element equal to the vector of r^2
    r_sq_simulator1 = r_sq_simulator1 + out
    
  }
  return(r_sq_simulator1/simulations)  #dividing by the number of simulations to get the average

}


average_r_sq_individual_based = function(p = c(0.25, 0.25, 0.25, 0.25), c , generations, N = rep(1000, generations), simulations = 5000){
  
  #inputs: p = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population to have
  #inputs: c = recombination rate
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i
  #inputs: simulations = number of simulations you want to average
  #output: average r_squared values for each generation over the (5000) simulations
  
  
  r_sq_simulator2 = numeric(generations + 1)
  
  for (i in 1:simulations){
    # wlog I am setting up an initial population matrix with size equal to N[1]
    # the proportion of the 4 alleles are p
    initial = constrained_matrix(p, N[1]) 
    
    
    #calling C simulator
    mode(N) = 'integer'
    mode(c) = 'numeric'
    mode(initial) = 'raw'
    result = .Call('sim_ld', initial, N, c )
    
    #result output is a list of size generations + 2 (becuase you have generations 0, ..., n and the r^2 at the end)
    r_sq_simulator2 = r_sq_simulator2 + as.numeric(unlist(result[generations + 2]))  
  }
  
  return(r_sq_simulator2/simulations)
}  


recurrence_1_Park = function(p = c(0.25, 0.25, 0.25, 0.25), c , generations, N = rep(1000, generations)){
  
  q = numeric(generations + 1)
  D = p[1]*p[4] - p[2]*p[3]
  q[1] = D^2/((p[1]+p[2])*(p[3]+p[4])*(p[1] + p[3])*(p[2] + p[4]) )
  for (i in 1:generations){
    
    #remark N starts from generation 1 not zero, small subtlety
    q[i+1] = 1/(2*N[i]) + (1 - 1/(2*N[i]))* ((1-c)^2)* q[i]
  } 
  
  return (q)
}


recurrence_2_Park = function(p = c(0.25, 0.25, 0.25, 0.25), c , generations, N = rep(1000, generations)){
  
  q = numeric(generations + 1)
  D = p[1]*p[4] - p[2]*p[3]
  q[1] = D^2/((p[1]+p[2])*(p[3]+p[4])*(p[1] + p[3])*(p[2] + p[4]) )
  N = c(N[1], N) #just adding population size at time zero, which we assume to be equal to N[1]
  
  for (i in 1:generations){
    
    #remark N starts from generation 1 not zero, small subtlety
    q[i+1] = q[i]*(1 - c/N[i])*(1 - 1/(2*N[i+1]))*(1-c)^2 + c^2 / (2*N[i]) + 1/(2*N[i+1])  
    } 
  
  return (q)
} 





comparison = function(p = c(0.25, 0.25, 0.25, 0.25), c , generations = 50, N = rep(1000, generations), simulations = 5000){
  
  #the function simply groups together the 4 r^2 vectors and returns a (generations + 1)x4 matrix
  out1 = average_r_sq_population_based(p, c, generations, N, simulations)
  out2 = average_r_sq_individual_based(p, c, generations, N, simulations)
  out3 = recurrence_1_Park(p, c, generations, N)
  out4 = recurrence_2_Park(p, c, generations, N)
  
  return( matrix(c(out1, out2, out3, out4), generations +1, 4))
}



#MAIN FUNCTION YOU NEED TO CALL IF YOU WANT TO PLAY AROUND WITH PARK'S EQUATIONS
plottings = function(p = c(0.25, 0.25, 0.25, 0.25), vec , generations = 50, N = rep(1000, generations), simulations = 5000){
  
  #inputs: p = vector of size 4,  represent the (0,0), (0,1), (1,0), (1,1) frequencies you want the initial population to have
  #inputs: c = recombination rate
  #inputs: N = vector of size = #generations, N[i] is the number of individuals in generation i
  #inputs: simulations = number of simulations you want to average
  #output: plots of the 4 r^2 for each value of recombination rate c provided
  
  for (c in vec)
  {
    m = comparison(p , c , generations, N, simulations)
    
    plot(m[, 1], xlab = 'Generations', ylab='r^2', ylim = c(0, 1.2*max(m[,1])), pch = 19)
    points(m[,2], col = 'red', pch = 19)
    points(m[,3], col = 'blue', pch = 19)
    points(m[,4], col = 'green', pch = 19)
    legend("topright", c( "Population based sim.", "Individual based sim.", 
                             "1st recursion Park", "2nd recursion Park") , 
                             col = c('black', 'red', 'blue', 'green') , pch = c(19,19,19,19))
    title(paste("N jumps from 1000 to 5000 at generation 25; c =" , as.character(c), sep = " "))
  }
}


#sample plots
#Here N is fixed at 1000
vec = c(0.1, 0.2, 0.3, 0.4, 0.5) #vector of recombination rates we want to plot
plottings(vec = vec)

#Here N jumps to 5000 at generation 25
gen = c(rep(1000, 25), rep(5000, 25))
plottings(vec = vec, generations = 50, N = gen)
