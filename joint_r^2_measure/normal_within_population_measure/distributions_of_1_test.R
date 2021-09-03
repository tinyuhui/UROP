p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25, 0.25, 0.25, 0.25)
N = rep(1000, 100)
M = rep(1000, 100)
c = 0.1
migration =  matrix(data = c(0.999,0.001, 0.001, 0.999), nrow = 2)

#out = simulation(p0, q0, N, M, c, migration)
#out2 = signal(p0, q0, N, M, c, migration, 2)
p_values2 = create_pdf(p0, q0, N, M, c, migration, simulations = 1000, no_breaks = 100,
                       title = 'start in equilibrium')

#####################################################################


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
  
  p_values[i,] = create_pdf(p0, q0, N, M, c, migration, simulations = 1000, no_breaks = 100,
                            title = 'start in equilibrium')
  compare_pdfs(p_values2, p_values[i,], 50, migration_value = m)
  compare_pdfs2(p_values2, p_values[i,], 200, migration_value = m)
  
}
