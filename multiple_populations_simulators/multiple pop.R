sim3.ld = function(populations = 3, generations = 20, p0 = matrix(c(0.25), ncol = 4, nrow  = populations),
                   N = matrix(1000, nrow = populations, ncol = generations),
                   c = 0.5, migration = diag(x = 1, nrow = populations, ncol = populations))
  #inputs: populations: number of different population of the system
  #inputs: generations: how many generations you want to simulate 
  #                     (we start from gen 0, so will get "generations + 1" output)
  
  #inputs: p0: populations x 4 matrix; each row is the initial haplotypes 
  #            frequencies at generation 0 for a population
  
  #inputs: N: populations x generations matrix, N[j,i] represents the size of 
  #           population j at generation i, BEFORE MIGRATION OCCURS
  
  #inputs: c is the recombination rate for the two loci

  #inputs: migration: FIXED migration matrix of size populations x populations;
  #        M[i,j] represents the proportion of population 1 that migrates
  #        to population 2 at each generation
  
  #output: list of size 3 containing: simulated haplotypes frequencies for all the 
  #                                   populations for all generations
  #                                   simulated D and r^2 values for all the populations for all generations

{
  
  #haplotye frequency 3d array. p[j,i,2] is the freq. of haplotype 2 at a time point i, for population j
  #dimension are populations x generations +1 x 4
  p = array(dim = c(populations, generations + 1, 4)) #3d array
  
  #initializing p at time zero, using initial conditions
  p[,1,]<-p0
  
  
  # matrix of raw LD measures. D[j,i] is the raw LD meas. for population j at generation i
  D = array(dim = c(populations,  generations +1))
  D[,1]<-p[,1,1]*p[,1,4]-p[,1,2]*p[,1,3]
  
  
  # PROPOGATION
  for (i in 1:t)
  { 
    #matrix containing probabilities used in multinomial distribution later
    exp_p = matrix(nrow = populations, ncol = 4)
    for (j in 1:populations){
      #entering new probabilities for generation i, using the population-based approximation
      exp_p[j, ] = p[j,i,] +c*c(-1, 1, 1, -1)*D[j,i]
    }
    
    #size[j,i] indicates how many gametes will move from pop j to pop i
    size = matrix(nrow = populations, ncol = populations)
    for (j in 1:populations){
      size[j,] = round(2*N[j,i]*migration[j,])
    }
    
    #Multinomial sampling for randoom genetic drift, generation i is reproducing
    
    #3d array, pop_flows[k,j,] is a 4x1 vector with the number of gametes of each type,
    #moving from pop k to pop j
    pop_flows = array(dim = c(populations,populations,4))
    for (k in 1:populations){
      for (j in 1:populations){
        pop_flows[k,j,] = rmultinom(1, size=size[k,j], prob=exp_p[k,]) 
                                                                        
      }
    }
    
    #matrix with entry [k,j] indicating the number of alles of type j at population k
    pop_after_migration = matrix(nrow = populations, ncol = 4)
    pop_after_migration = colSums(pop_flows, dim = 1)
    
    
    #getting haplotypes frequencies from pop_after_migration, by simply normalising each row
    for (j in 1:populations){
    # Calculating vectors of frequencies, after migration occurs
      p[j,i+1,] = pop_after_migration[j,]/(pop_after_migration[j,1] + pop_after_migration[j,2] + 
                                             pop_after_migration[j,3] + pop_after_migration[j,4])
      
      
      # Calculate raw LD for the generation
      D[j,i+1]<-p[j,i+1,1]*p[j,i+1,4]-p[j,i+1,2]*p[j,i+1,3]
      
    }
  }
  
  # Calculate r2 for both population
  r2 = matrix(nrow = populations, ncol = generations + 1)
  for (j in 1:populations){
    r2[j,] = D[j,]^2/((p[j,,1]+p[j,,2])*(p[j,,3]+p[j,,4])*(p[j,,1]+p[j,,3])*(p[j,,2]+p[j,,4]))
  }
  
  return(list(p=p, D=D, r2=r2))
}