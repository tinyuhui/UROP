#need this for 3d plot
library(car)

#importing "average_r_sq_at_equilibrium" function "from simulation2.cpp"
library("Rcpp")
Rcpp::sourceCpp("C:/Users/angus/OneDrive/Desktop/urop/C++/simulation_2.cpp")




#Recombination rates we want
vals_of_c = c(0.01,0.02,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

#Migration rates we want
vals_of_migration = c(0.01,0.02,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

#meshgrid of the 13x10 values of (c, migration) we want to analyse
x = rep(c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),10)
y = c(rep(0.01,13),rep(0.03,13),rep(0.05,13),rep(0.1,13),rep(0.15,13),rep(0.2,13),
      rep(0.25,13),rep(0.3,13), rep(0.4,13),rep(0.5,13))

n = length(x)

#population sizes
N_val = 1000
M_val = 1000

simulations = 100




#main function:
# It takes fixed population sizes N_val, M_val, simulates 100 generations of them
# for "simulation" times, takes the average of the last 50 generations (when we are at equilibrium already)
# Now we have equilibrium values of r^2 for each (c, migration) value specified
# We use Park (2) to estimate "N" from the value N_val and "c"
# Then we see how much it is biased upwards by computing the ratio of the estimated N, with the real N

bias = function(N_val, M_val, x, y, simulations){
  #inputs: N_val = size of population 1
  #        M_val = size of population 2
  #        x = meshgrid of the c values we want (as above)
  #        y = meshgrid of the migration values we want (as above)
  #        simulations
  #outputs: list of 4:
  #         1. vector of computed r_sq values
  #         2. vector of computed N_hat/N ratios
  #         3., 4. same things but in matrix form nrow = #migration_rates ncol = #recombination_rates
  
  p0 = c(0.25,0.25,0.25,0.22)
  q0 = c(0.25,0.25,0.25,0.25)
  N = rep(N_val, 100)
  M = rep(M_val, 100)
  n = length(x)
  
  r = numeric(n)
  
  z = numeric(n)
  
  for (i in 1:n){
    #initializing migraton matrix
    migration = matrix(data = c(1 - y[i],y[i], y[i], 1-y[i]), nrow = 2)
    
    #calculating r_sq at equilibrium
    r[i] = average_r_sq_at_equilibrium(p0, q0, N, M, x[i], migration, simulations, 50)
    
    
    #solving the quadratic equation Park (2) to estimate N
    a = 2*(r[i]) - 2*(r[i]) * ((1-x[i])^2)
    b = (r[i]) * ((1-x[i])^2) * (2*x[i] + 1) - (1+x[i]^2)
    d = r[i] * (1-x[i])^2 * x[i]
    
    z[i] = quadratic(a,b,d)[1] / N_val
    
  }
  
  #turning vectors into matrices for readability
  r2 = matrix(data = r, nrow =  10, byrow = TRUE)
  z2 = matrix(data = z, nrow =  10, byrow = TRUE)
  return(list(r, z, r2, z2))
}

#3d plot, not very easy to visualize
scatter3d(x, y, z, df.smooth = 15, fit = "smooth", xlab = "recombination rate", ylab = "migration rate")

out = bias(N_val, M_val, x, y, simulations)
z2 = unlist(out[2])




#turning matrix into df
df = data.frame(z2, row.names = c("m = 0.01","m = 0.03","m = 0.05","m = 0.1",
                                  "m = 0.15","m = 0.2","m = 0.25","m = 0.3",
                                  "m = 0.4","m = 0.5"))
vec = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
for (i in 1:length(vec)){
  vec[i] = paste("c =", as.character(vec[i]), sep = " ")
}

colnames(df) = vec


write.csv(df,"C:\\Users\\angus\\OneDrive\\Desktop\\urop\\migration_bias.csv")



#function that predicts r^2 given N and c using Park (2)
pred_r2 = function(c, N){
  return((1+c^2)/(2*N)/(1-(1-c/N)*(1-1/(2*N))*(1-c)^2))
}

































quadratic <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
    return(result)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
    return(x)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}