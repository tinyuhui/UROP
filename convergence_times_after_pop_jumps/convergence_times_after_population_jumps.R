
convergence_plotter = function(M, N, precision = 0.1, vec_length = 100, C = seq(from = 0.01, to = 0.5, by = 0.01)){
  #inputs: M = population before the jump, we assume r^2 is in equilibrium (value obtained by solving
  #            equation (1) in Park 2012)
  #inputs: N = population after the jump
  #inputs: precision in [0,1] = you consider r^2 to have reached convergence, when the r^2value is within
  #                             a factor of "precision" to its equilibrium value
  #                             ex. if equilibrium value is 0.025 and precision is 0.1
  #                                 we stop when r^2 enters [0.0225, 0.0275]
  #inputs: vec_length = how many generations you want to simulate, simulation is very fast, 
  #                     as I am using Park (2012) recurrence. If c is low, it's better to use a high number 
  #                     of vec_length to ensure convergence even for low "precision"
  #inputs: C = vector of recombination rates you want to simulate
  
  #outputs: a plot of the #generations needed for convergence, against recombination rates c
  #outputs: a vector f the #generations needed for convergence
  
  out = numeric(length(C))
  for (j in 1:length(C)){
    
    c = C[j]
    q = numeric(vec_length)
    q[1] = (c^2/(2*M) + 1/(2*M))/(1 - (1- c/M)*(1 - 1/(2*M))*(1 - c)^2)
    
    for (i in 1: (vec_length - 1)){
      q[i+1] = q[i]*(1 - c/N)*(1 - 1/(2*N))*(1-c)^2 + c^2 / (2*N) + 1/(2*N)  
    }
    
    equilibrium = q[vec_length]
    
    for (i in 2: (vec_length)){
      if (abs(q[i] - equilibrium) < precision*equilibrium){
        out[j] = i - 1
        break
      }  
        
    
    }
  }
  
  plot(C,out, xlab = 'recombination rate c', ylab='generations for convergence', 
       pch = 19)
  title(paste("N jumps from " , as.character(M), "to ", as.character(N), ", precision = ", 
              as.character(precision*100), "%",
              sep = " "))
  
  return(out)
}