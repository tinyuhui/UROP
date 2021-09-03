simple_gene_drive_simulator = function(drive = 0.95, pop_size = 1000, initial_vector_size = 100, iterations = 20){
  
  male = c(pop_size/2, 0, initial_vector_size) 
  
  female = c(pop_size/2, 0, 0)   
  
  p = male / (pop_size/2 + initial_vector_size)
  
  q = female / (pop_size/2)
  
  f = male / (pop_size + initial_vector_size)
  
  g = female / (pop_size + initial_vector_size)
  
  out = matrix(nrow = iterations, ncol = 5)
  
  out[1,] = c(f, g[1:2])
  
  for (i in 2:iterations){
    
    male = c(p[1]*q[1] + (1-drive)*(p[2]*q[1] + p[1]*q[2])
             +((1-drive)^2)*p[2]*q[2],
             drive*(p[1]*q[2] + p[2]*q[1]) + 
               2*drive*(1-drive)*(p[2]*q[2]) +
               (1-drive)*(p[2]*q[3] + p[3]*q[2])+
               p[1]*q[3] + p[3]*q[1],
             (drive^2)*p[2]*q[2] + 
               drive*(p[2]*q[3] + p[3]*q[2]) + 
               p[3]*q[3])/2
    female = c(male[1], male[2], 0)
    
    p = male/sum(male)
    q = female/sum(female)
    
    f = male/(sum(male) + sum(female))
    g = female/(sum(male) + sum(female))
    
    out[i,] = c(f, g[1:2])
    
    
  }
  plots(out)
  return (out)
}



plots = function(m){
  
    
  plot(m[, 1], xlab = 'Generations', ylab='frequecny', ylim = c(0, 1), pch = 19)
  points(m[,2], col = 'red', pch = 19)
  points(m[,3], col = 'blue', pch = 19)
  points(m[,4], col = 'green', pch = 19)
  points(m[,5], col = 'orange', pch = 19)
  legend("topright", c( "M AA", "M AB", 
                           "M BB", "F AA", "F AB") , 
         col = c('black', 'red', 'blue', 'green', 'orange') , pch = c(19,19,19,19,19))

}  

