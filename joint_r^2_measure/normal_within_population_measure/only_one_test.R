## BETWEEN-POP R^2 

#importing the C++ simulator
Rcpp::sourceCpp("C:/Users/angus/OneDrive/Desktop/urop/C++/simulation_2.cpp")




likelihood_test = function(p, q, size1, size2){
  # Given the 2 vector of alleles frequencies and the final size2 
  # of the two population, and the degrees of freedom of the chi_sq distribution
  # It returns the p-value for the likelihood ratio test I developed (see pdfs),
  # which is indeed our joint r^2 measure
  
  
  #initializing the vectors that contain the predicted frequencies under H0
  # i.e. they are the product of the 3 marginals for each cell
  exp = numeric(4)
  
  
  exp = (size1*p[1] + size2*q[1])/(size1 + size2)
  
  test = numeric(8)
  
  
  test = (-1 + 2* pnorm(sqrt((size1 + size2) / (exp*(1-exp)))*abs(exp - p[1])))
  
  return(test)
}





signal = function(p0, q0, N, M, c, migration){
  #Given the initial parameters, and the degrees of freedom of the chi_sq,
  #it simulates the 2 populations, then considers the final generation (where 
  #we assume we have reached equilibrium), and performs a LR test on the final
  #frequencies, and returns the p-value for the test
  
  #running the individual based simulator
  out = simulation(p0, q0, N, M, c, migration)
  
  #find the population sizes at the final generation
  final_generation_number = length(N)
  N_size = N[final_generation_number]
  M_size = M[final_generation_number]
  #print(c('initial p', p0))
  #print(c('initial q', q0))
  
  #get the final generation allele frequencies for the 2 pop
  p_values = matrix(unlist(out[1]), nrow = 4)[1:4,final_generation_number]
  q_values = matrix(unlist(out[2]), nrow = 4)[1:4,final_generation_number]
  #print(c('p', p_values))
  #print(c('q', q_values))
  #print(c('N_size', N_size))
  #print(c('M_size', M_size))
  
  #carry out the LR test
  return (likelihood_test(p_values, q_values, N_size, M_size))
}


create_pdf = function(p0, q0, N, M, c, migration, simulations, no_breaks, title = 'no title'){
  
  # Given initial inputs for the individual based simulator, the degrees of freedom
  # of the chi_sq distr, how many simulations to carry out, and how many breaks (i.e. bars)
  # the histogram should have,
  # It plots an approximate probability distribution function for the test statistics, given these inputs
  # And It returns a vector containing all the simulated p-values
  
  #doing the iterations
  values_of_tests = numeric(simulations)
  for (i in 1:simulations){
    values_of_tests[i] = signal(p0, q0, N, M, c, migration)
  }
  
  #df = approxfun(density(values_of_tests))
  #plot(density(values_of_tests))
  
  #plotting the histogram
  if (title != 'no title'){
    tit = paste("m =", as.character(migration[1,2]), title, sep = " ")
  }else{
    tit = 'Distribution of the simulated r^2 values'
  }
  
  hist(values_of_tests, breaks = no_breaks, col = 'lightblue', xlab = 'r^2 simulated values', 
       main = tit)
  
  return(values_of_tests)
  
}


compare_pdfs = function(simulated1, simulated2, number_breaks, migration_value = 'unknown'){
  #We want to see how well the test manages to catch migration signals
  #We plot two pdfs for different migration rates and see how well-separated they are
  #In this case we plot the two histograms (with separate width of the bars), and then put it on the same
  #axis
  
  #inputs simulated1 and 2 are vectors containing all the simulated p_values using "create_pdf"
  #number_breaks tells us how many bars should each histogram have
  c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
  c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
  hist1 = hist(simulated1, breaks = number_breaks, plot = FALSE)
  hist2 = hist(simulated2, breaks = number_breaks, plot = FALSE)
  range_vec = range(c(hist1$breaks, hist2$breaks))
  maximum = max(c(hist1$count, hist2$count))
  plot(hist1, col = c1, xlim = range_vec, ylim = c(0,maximum), xlab = 'r^2 simulated values',
       main = paste("m =", as.character(migration_value), sep = " "))
  plot(hist2, add = TRUE, col = c2)
  
}

compare_pdfs2 = function(simulated1, simulated2, number_breaks, migration_value = 'unknown'){
  #We want to see how well the test manages to catch migration signals
  #We plot two pdfs for different migration rates and see how well-separated they are
  #In this case we plot the two histograms (with a common width of the bars), and then put it on the same
  #axis
  
  #inputs simulated1 and 2 are vectors containing all the simulated p_values using "create_pdf"
  #number_breaks tells us how many bars should each histogram have
  c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
  c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
  b <- min(c(simulated1, simulated2)) - 0.01  # Set the minimum for the breakpoints
  #print(b)
  e <- max(c(simulated1, simulated2)) + 0.1 # Set the maximum for the breakpoints
  #print(e)
  ax <- seq(b,e, length.out = number_breaks)
  #print(ax)
  hist1 = hist(simulated1, breaks = ax, plot = FALSE)
  hist2 = hist(simulated2, breaks = ax, plot = FALSE)
  maximum = max(c(hist1$count, hist2$count))
  plot(hist1, col = c1, xlim = c(b,e), ylim = c(0,maximum), xlab = 'r^2 simulated values',
       main = paste("m =", as.character(migration_value), sep = " "))
  plot(hist2, add = TRUE, col = c2)
}