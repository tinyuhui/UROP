

two_loci_population_simulator = function(p1, p2, p3, recombination_rate, iterations, n){
  
  #inputs: p_i are the 4 frequencies for the two biallelic loci (p4 ommitted to avoid overparametrization)
  #inputs: recombination_rate (also called 'c')
  #inputs: n is the size of the diploid population
  
  calculate_r_sq_population_based = function(vec){
    #input: vector of 4 bi-allelic frequencies
    #output: r^2
    
    D = prob[1]*prob[4] - prob[2]*prob[3]
    r_sq = D^2 / ((prob[1] + prob[3])*(prob[2] + prob[4])*(prob[1] + prob[2])*(prob[4] + prob[3]))
    return(r_sq)
  }
  
  
  #initial vector of frequencies
  prob = c(p1, p2, p3, 1 - p1 - p2 - p3)
  
  #initializing empty vector (to be filled with r_sq values)
  r_sq = numeric(iterations)
  c = recombination_rate
  vec = c(-c, c, c, -c)
  
  for (i in 1:iterations){
    
    #calculating raw LD measure
    D = prob[1]*prob[4] - prob[2]*prob[3]
    
    #updating probabilities due to recombination
    prob = prob + D*vec
    
    #simulating new generation frequencies with multinomial distr
    prob = rmultinom(n = 1, size = 2*n, prob = prob)/(2*n)
    
    #calculating r^2 for generation i
    r_sq[i] = calculate_r_sq_population_based(prob)
    
  }
  return(r_sq)
}