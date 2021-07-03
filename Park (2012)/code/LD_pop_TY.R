# POPULATION BASED LD SIMULATOR
# ONLY TRACK FOUR HAPLOTYPE FREQUENCIES FORWARD IN TIME
# NUMBER OF GENERATIONS DEPENDS ON length(N)
# c IS RECOMBINATION RATE
# N[1] IS POPULATION SIZE
# N[2] OS NUMBER OF ITERATIONS

sim.ld<-function(p0=c(0.25, 0.25, 0.25, 0.25), N=rep(1000, 20), c=0.5)
{
  t<-length(N)
  # HAPLOTYPE FREQ MATRIX. EACH COLUMN IS THE HAPLOTYPE FREQ AT A TIME POINT. 4 ROWS, (t+1) COLUMNS. 
  p<-matrix(nr=4, nc=t+1)
  p[,1]<-p0
  # RAW LD MEASURE
  D<-rep(NA, t+1)
  D[1]<-p[1,1]*p[4,1]-p[2,1]*p[3,1]
  # PROPOGATION
  for (i in 1:t)
  {
    # CALCULATE THE EXPECTED GAMETIC FREQ IN THE GAMETE POOL, AFTER RECOMBINATION
    exp_p<-p[,i]+c(-1, 1, 1, -1)*c*D[i]
    # MULTINOMIAL SAMPLING FOR RANDOM GENETIC DRIFT
    p[,i+1]<-rmultinom(1, size=2*N[i], prob=exp_p)/(2*N[i])
    # CALCULATE THE NEW D
    D[i+1]<-p[1,i+1]*p[4,i+1]-p[2,i+1]*p[3,i+1]
  }
  # CALCULATE r2
  r2<-D^2/((p[1,]+p[2,])*(p[3,]+p[4,])*(p[1,]+p[3,])*(p[2,]+p[4,]))
  return(list(p=p, D=D, r2=r2))
}