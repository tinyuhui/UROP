try1 = sim2.ld(p0 = c(0.5, 0, 0, 0.5), c = 0.01, migration = matrix(c(0.7,0.3,0,1), nc = 2, byrow = T))

plot(unlist(try1[5]), pch = 19 , xlab = 'Generations', ylab='r^2', ylim = c(0, 1))
points(unlist(try1[6]), col = 'red', pch = 19)

legend("topright", c( "Population 1", "Population 2") , 
       col = c('black', 'red') , pch = c(19,19))
title(paste("Population 1 in complete disequilibrium, c = 0.01, migration only from 1 to 2"))

try2 = sim2.ld(p0 = c(0.5, 0, 0, 0.5), c = 0.25, migration = matrix(c(0.8,0.2,0.2,0.8), nc = 2, byrow = T))

plot(unlist(try2[5]), pch = 19 , xlab = 'Generations', ylab='r^2', ylim = c(0, 1))
points(unlist(try2[6]), col = 'red', pch = 19)

legend("topright", c( "Population 1", "Population 2") , 
       col = c('black', 'red') , pch = c(19,19))
title(paste("P1 in complete disequilibrium, P2 in eq., c = 0.25, symm migration"))

try3 = sim2.ld(p0 = c(0.5, 0, 0, 0.5), c = 0.01, migration = matrix(c(0.8,0.2,0.2,0.8), nc = 2, byrow = T))

plot(unlist(try3[5]), pch = 19 , xlab = 'Generations', ylab='r^2', ylim = c(0, 1))
points(unlist(try3[6]), col = 'red', pch = 19)

legend("topright", c( "Population 1", "Population 2") , 
       col = c('black', 'red') , pch = c(19,19))
title(paste("P1 in complete disequilibrium, P2 in eq., c = 0.01, symm migration"))
