# We start with three examples, for each we plot the approximate distribution
# of the test statistic
# 1. two populations with 1% migration, started with very different initial freq
# the values of r^2 after 100 generations are around 0.4
# 2. two populations, both started in equilibrium, but with 0.1% migration.
#     This will be our base case when comparing different migration rates
# 3. This case is exactly the same as above, but we have migration of 0.5, so 
#   essentially one population

# We then compare the distribution of case 2 with different migration rates, 
# to see how well they are separable. Results are promising

#################################################
p0 = c(0,0.99,0,0.1)
q0 = c(0.1,0,0.99,0)
N = rep(5000, 100)
M = rep(5000, 100)
c = 0
migration =  matrix(data = c(0.99,0.01, 0.01, 0.99), nrow = 2)

#out = simulation(p0, q0, N, M, c, migration)
#out2 = signal(p0, q0, N, M, c, migration, 2)
p_values = create_pdf(p0, q0, N, M, c, migration, 2,simulations = 1000, no_breaks = 100,
                      title = 'start in disequilibrium')

######################################################################
#independent-ish populations 
p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25, 0.25, 0.25, 0.25)
N = rep(1000, 100)
M = rep(1000, 100)
c = 0.1
migration =  matrix(data = c(0.999,0.001, 0.001, 0.999), nrow = 2)

#out = simulation(p0, q0, N, M, c, migration)
#out2 = signal(p0, q0, N, M, c, migration, 2)
p_values2 = create_pdf(p0, q0, N, M, c, migration, dof = 2, simulations = 1000, no_breaks = 100,
                       title = 'start in equilibrium')

########################################################################
#full-migration 0.003, 0.0004, 0.007, 0.0005
p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25, 0.25, 0.25, 0.25)
N = rep(1000, 50)
M = rep(1000, 50)
c = 0.1
migration =  matrix(data = c(0.5,0.5, 0.5, 0.5), nrow = 2)

p_values3 = create_pdf(p0, q0, N, M, c, migration, dof = 2, simulations = 1000, no_breaks = 100,
                       title = 'start in equilibrium')

compare_pdfs(p_values2, p_values3, 50, migration_value = 0.5)
compare_pdfs2(p_values2, p_values3, 200, migration_value = 0.5)

############################################################################

migration_rates = c(0.05,0.1,0.2,0.3,0.4,0.5)
simulations = 1000
p_values = matrix(numeric(1), nrow = 6, ncol = simulations)
for (i in 1:6){
  m = migration_rates[i]
  p0 = c(0.25,0.25,0.25,0.25)
  q0 = c(0.25, 0.25, 0.25, 0.25)
  N = rep(1000, 100)
  M = rep(1000, 100)
  c = 0.1
  migration =  matrix(data = c(1-m,m,m,1-m), nrow = 2)
  
  p_values[i,] = create_pdf(p0, q0, N, M, c, migration, dof = 2, simulations = 1000, no_breaks = 100,
                         title = 'start in equilibrium')
  compare_pdfs(p_values2, p_values[i,], 50, migration_value = m)
  compare_pdfs2(p_values2, p_values[i,], 200, migration_value = m)
  
}

