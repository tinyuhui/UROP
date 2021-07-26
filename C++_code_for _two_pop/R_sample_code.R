library("Rcpp")

# Here you can evaluate cpp expressions directly in R
evalCpp("2 + 2")
evalCpp("Rcpp::seq(1, 4)")  
evalCpp("seq(1, 4)")        

#Compiling the code, you can now call all the functions in the script from your R console
#If you use Rstudio the functions will appear in the global environment

Rcpp::sourceCpp("./simulation.cpp")   #add relevant path here


#defining the inputs of simulation function
p0 = c(0.5,0.3,0,0.2)
q0 = c(0.1,0.3,0,0.6)
N = rep(5000, 100)
M = rep(5000, 100)
c = 0.2
migration =  matrix(data = c(0.8,0.2, 0.2, 0.8), nrow = 2)

simulation(p0, q0, N, M, c, migration)


#more extreme case with 1 million gametes, takes 10-15 secs on my computerp0 = c(0.5,0.3,0,0.2)
q0 = c(0.1,0.3,0,0.6)
N = rep(500000, 100)
M = rep(500000, 100)
c = 0.2
migration =  matrix(data = c(0.8,0.2, 0.2, 0.8), nrow = 2)

simulation(p0, q0, N, M, c, migration)