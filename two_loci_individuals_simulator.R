

produce_gamete = function(old_matrix, recombination_rate, parent){

  # inputs: a matrix with 2n rows, representing the genome of n diploid organisms
  # inputs : parent in {1, ..., n} (which corresponds to rows 2m-1, 2m for some m <= n), 
  # inputs: recombination rate in [0, 0.5], corresponds to the probability of an odd number of crossing overs happening between the two loci
  # output: 2x1 matrix corresponding to a gamete of the parent
  
    x = runif(1)
  if (x < recombination_rate/2){
    gamete_1 = c(old_matrix[2*parent - 1, 1], old_matrix[2*parent, 2])   #the first allele comes from row 2n -1, the second from row 2n
  } else if (x < recombination_rate){
    gamete_1 = c(old_matrix[2*parent, 1], old_matrix[2*parent - 1, 2])  #same as above but switched
    
  } else if (x < (1 + recombination_rate)/2){
    gamete_1 = c(old_matrix[2*parent - 1, 1], old_matrix[2*parent - 1, 2])  #no crossing over, both alleles from row 2n -1
  } else {
    gamete_1 = c(old_matrix[2*parent, 1], old_matrix[2*parent, 2])
  }
  return(gamete_1)
}



calculate_r_sq = function(matr){  
  
  #input: m x 2 binary matrix (only 2 loci for now)
  #output: r^2 for the 2 loci
  
  m = nrow(matr)
  p1 = 0
  p2 = 0
  p3 = 0
  p4 = 0
  for (i in 1:m){
    if( matr[i, 1] == 0 & matr[i, 2] == 0){
      p1 = p1 + 1
    }else if (matr[i, 1] == 0 & matr[i, 2] == 1){
      p2 = p2 + 1
    }else if (matr[i, 1] == 1 & matr[i, 2] == 0){
      p3 = p3 + 1
    } else {
      p4 = p4 + 1
    }
  }
  p1 = p1/m
  p2 = p2/m
  p3 = p3/m
  p4 = p4/m
  
  raw_LD = p1*p4 - p2*p3
  r_sq = raw_LD^2/((p1 + p2)*(p3 + p4)*(p1 + p3)*(p2 + p4))
  return(r_sq)
}






simulator = function(matrix, recombination_rate, iterations){
  
  #inputs: 2n x 2 binary matrix, represents genomic data at time zero
  #inputs: recombination rate in [0, 0.5], corresponds to the probability of an odd number of crossing overs happening between the two loci
  #inputs: number of iterations/generations to simulate
  #output: vector of size equal to 'iterations', representing the evolution of the r^2
  
  r_squared = numeric(iterations)
  old_matrix = matrix
  n = nrow(matrix)/2
  new_matrix = matrix(0,2*n, 2)
  
  for (j in 1:iterations){
    for (i in 1:n){
      parents = sample(1:n, 2, replace = FALSE)
      print(c(parents, 'parents'))
      haplotype_1 = produce_gamete(old_matrix, recombination_rate, parents[1])
      haplotype_2 = produce_gamete(old_matrix, recombination_rate, parents[2])
      print(c(haplotype_1, haplotype_2))
      new_matrix[2*i -1, ] = haplotype_1
      new_matrix[2*i, ] = haplotype_2
    }
    
    print(new_matrix)
    r_squared[j] = calculate_r_sq(new_matrix)
    old_matrix = new_matrix
    new_matrix = matrix(0,2*n, 2)
  }
  return(r_squared)
}