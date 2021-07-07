constrained_matrix = function(p, n){
  
  #inputs: p = c(p1, p2, p3, p4) in [0,1], represent the (0,0), (0,1), (1,0) frequencies you want the simulated matrix to have
  #inputs: n is the number of diploid individuals
  #output: a simulated matrix with the required p1, p2, p3 frequencies. 
  
  #Formally let Row_i be uniform on {(0,0), ... (1,1)}. Then conditioned on the number of (0,0)s, ... , (1,1)s, 
  #the random matrix sampled below has the correct distribution
  
  p1 = p[1]
  p2 = p[2]
  p3 = p[3]
  
  matr = matrix(numeric(1), 2*n, 2)
  q1 = floor(2*n*p1)
  q2 = floor(2*n*p2)
  q3 = floor(2*n*p3)
  q4 = 2*n - q1 - q2 -q3
  
  if(abs(q4 - floor(2*n*(1 - p1 - p2 -p3)))/ (2*n) > 0.1){
    print("imprecise approximation due to discretization error")
    
    
  }
  
  s = sample(1:(2*n))
  
  for (i in 1:q1){
    matr[s[i],] = c(0,0)
  }
  
  for (i in (q1 +1):(q1 + q2)){
    matr[s[i],] = c(0,1)
  }
  
  for (i in (q1 + q2 +1):(q1 + q2 +q3)){
    matr[s[i],] = c(1,0)
  }
  
  if ((q1 + q2 + q3 +1) <= (2*n)){
    for (i in (q1 + q2 + q3 +1):(2*n)){
      matr[s[i],] = c(1,1)
    }
  }
  
  return(matr)
}