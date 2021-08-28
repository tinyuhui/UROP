#################################################
p0 = c(0,0.99,0,0.1)
q0 = c(0.1,0,0.99,0)
N = rep(5000, 100)
M = rep(5000, 100)
c = 0
migration =  matrix(data = c(0.99,0.01, 0.01, 0.99), nrow = 2)

#out = simulation(p0, q0, N, M, c, migration)
#out2 = signal(p0, q0, N, M, c, migration, 2)
p_values = create_pdf(p0, q0, N, M, c, migration,simulations = 1000, no_breaks = 100,
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
p_values2 = create_pdf(p0, q0, N, M, c, migration, simulations = 1000, no_breaks = 100,
                       title = 'start in equilibrium')


########################################################################
#full-migration 0.003, 0.0004, 0.007, 0.0005
p0 = c(0.25,0.25,0.25,0.25)
q0 = c(0.25, 0.25, 0.25, 0.25)
N = rep(1000, 50)
M = rep(1000, 50)
c = 0.1
migration =  matrix(data = c(0.5,0.5, 0.5, 0.5), nrow = 2)

p_values3 = create_pdf(p0, q0, N, M, c, migration, simulations = 1000, no_breaks = 100,
                       title = 'start in equilibrium')

