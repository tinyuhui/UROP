sim2.ld = function(p0 = c(0.25, 0.25, 0.25, 0.25), q0 = c(0.25, 0.25, 0.25, 0.25), N = rep(1000, 20),
                   M =rep(1000, 20), c = 0.5, migration = matrix(c(1,0,0,1), nc = 2, byrow = T))
  
  #inputs: p0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 1
  #inputs: q0: 4x1 vector; initial haplotypes frequencies at generation 0 for population 2
  #inputs: N: vector representing the size of population 1 at generation i, BEFORE MIGRATION OCCURS
  #inputs: M: vector representing the size of population 2 at generation i, BEFORE MIGRATION OCCURS
  # M and N should have the same size, which is the number of generations to simulate
  #inputs: c is the recombination rate for the two loci
  #inputs: migration: FIXED migration matrix; M[1,2] represents the proportion of population 1 that migrates
  #        to population 2 at each generation
  #inputs: initial_pop_sizes: population sizes of generation 0. I didn't put it together with N and M,
  #                           so that there is agreement with the earlier simulator (where recombination
  #                           and migration are inverted)
  
  #output: list of size 6 containing: simulated haplotypes frequencies for the two populations for all generations
  #                                   simulated D and r^2 values for the two populations for all generations

{
  t<-length(N)
  
  #haplotye frequency matrix for population 1. Each columns is the haplotype freq. at a time point
  #4 rows, t columns
  p<-matrix(nr=4, nc=t)
  p[,1]<-p0
  
  #same for population 2
  q<-matrix(nr=4, nc=t)
  q[,1]<-q0
  
  #Raw LD vector for population1, size = t
  D1<-rep(NA, t)
  D1[1]<-p[1,1]*p[4,1]-p[2,1]*p[3,1]
  
  D2<-rep(NA, t)
  D2[1]<-q[1,1]*q[4,1]-q[2,1]*q[3,1]
  
  #vector s.t size[i] indicates: 
  size1 = floor(2*N*migration[1,1]) #how many gametes originated in pop1 stay there at generation i
  size2 = 2*N - size1               #how many gametes originated in pop1 migrate at generation i
  
  size3 = floor(2*M*migration[2,2]) #how many gametes originated in pop2 stay there at generation i
  size4 = 2*M - size3               #how many gametes originated in pop2 migrate at generation i
  
  # PROPOGATION
  for (i in 1:(t - 1))
  {
    #Expected gametic freq in the two population, before recombination, and before migration
    #No recombination has yet occurred
    exp_p<-p[,i]
    exp_q<-q[,i]
    
    #Multinomial sampling for individuals who are migrating, migration happening
    
    population_1_stays <- rmultinom(1, size=size1[i], prob=exp_p)
    population_1_migrates <- rmultinom(1, size=size2[i], prob = exp_p)
    
    population_2_stays <- rmultinom(1, size=size3[i], prob=exp_q)
    population_2_migrates <- rmultinom(1, size=size4[i], prob = exp_q)
    
    
    # Calculating vectors of frequencies, after migration occurs, but before reproduction
    # These are temporary vectors, because we only store p, q and D1, D2 after reproduction
    p_temp = (population_1_stays + population_2_migrates)/(size1[i] + size4[i])
    q_temp = (population_2_stays + population_1_migrates)/(size2[i] + size3[i])
    
    # Calculate raw LD before reproduction
    D1_temp<-p_temp[1]*p_temp[4]-p_temp[2]*p_temp[3]
    D2_temp<-q_temp[1]*q_temp[4]-q_temp[2]*q_temp[3]
    
    #Updating frequencies due to recombination / genetic drift
    exp_p = p_temp + c(-1, 1, 1, -1)*c*D1_temp
    exp_q = q_temp + c(-1, 1, 1, -1)*c*D2_temp
    
    #Allele frequencies after reproduction
    p[,i+1] = rmultinom(1, size= N[i+1], prob=exp_p) / N[i+1]
    q[,i+1] = rmultinom(1, size= M[i+1], prob=exp_q) / M[i+1]
    
  }
  
  # Calculate raw LD for both populations
  D1 <- p[1,]*p[4,]-p[2,]*p[3,]
  D2 <- q[1,]*q[4,]-q[2,]*q[3,]
  
  # Calculate r2 for both populations
  r2_1<-D1^2/((p[1,]+p[2,])*(p[3,]+p[4,])*(p[1,]+p[3,])*(p[2,]+p[4,]))
  r2_2<-D2^2/((q[1,]+q[2,])*(q[3,]+q[4,])*(q[1,]+q[3,])*(q[2,]+q[4,]))
  
  return(list(p=p,q=q, D1=D1, D2=D2, r2_1=r2_1, r2_2=r2_2))
}