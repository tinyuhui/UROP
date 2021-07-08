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
  
  #output: list of size 6 containing: simulated haplotypes frequencies for the two populations for all generations
  #                                   simulated D and r^2 values for the two populations for all generations
  
{
  t<-length(N)

  #haplotye frequency matrix for population 1. Each columns is the haplotype freq. at a time point
  #4 rows, t+1 columns
  p<-matrix(nr=4, nc=t+1)
  p[,1]<-p0
  
  #same for population 2
  q<-matrix(nr=4, nc=t+1)
  q[,1]<-q0
  
  #Raw LD vector for population1, size = t+1
  D1<-rep(NA, t+1)
  D1[1]<-p[1,1]*p[4,1]-p[2,1]*p[3,1]
  
  D2<-rep(NA, t+1)
  D2[1]<-q[1,1]*q[4,1]-q[2,1]*q[3,1]
  
  # PROPOGATION
  for (i in 1:t)
  {
    #Calculate the expected gametic freq in the two population, after recombination, but before migration
    exp_p<-p[,i]+c(-1, 1, 1, -1)*c*D1[i]
    exp_q<-q[,i]+c(-1, 1, 1, -1)*c*D2[i]
    
    #Multinomial sampling for randoom genetic drift, generation i is reproducing
    
    size1 = floor(2*N[i]*migration[1,1]) #how many gametes originated in pop1 stay there
    size2 = 2*N[i] - size1               #how many gametes originated in pop1 migrate
    population_1_stays <- rmultinom(1, size=size1, prob=exp_p)
    population_1_migrates <- rmultinom(1, size=size2, prob = exp_p)
    
    size3 = floor(2*M[i]*migration[2,2]) #how many gametes originated in pop2 stay there
    size4 = 2*M[i] - size3               #how many gametes originated in pop2 migrate
    population_2_stays <- rmultinom(1, size=size3, prob=exp_q)
    population_2_migrates <- rmultinom(1, size=size4, prob = exp_q)

        
    # Calculating vectors of frequencies, after migration occurs
    p[,i+1] = (population_1_stays + population_2_migrates)/(size1 + size4)
    q[,i+1] = (population_2_stays + population_1_migrates)/(size2 + size3)
  
    # Calculate raw LD for the generation
    D1[i+1]<-p[1,i+1]*p[4,i+1]-p[2,i+1]*p[3,i+1]
    D2[i+1]<-q[1,i+1]*q[4,i+1]-q[2,i+1]*q[3,i+1]
  }
  
  # Calculate r2 for both population
  r2_1<-D1^2/((p[1,]+p[2,])*(p[3,]+p[4,])*(p[1,]+p[3,])*(p[2,]+p[4,]))
  r2_2<-D2^2/((q[1,]+q[2,])*(q[3,]+q[4,])*(q[1,]+q[3,])*(q[2,]+q[4,]))
  
  return(list(p=p,q=q, D1=D1, D2=D2, r2_1=r2_1, r2_2=r2_2))
}