X <- lhs::maximinLHS(n=50,k=2)
Sigma <- diag(1,50)
for(i in 1:49) {
  for(j in (i+1):50) {
    cor.ij <- exp(-sum((X[i,]-X[j,])^2)/.1)
    Sigma[i,j] <- cor.ij
    Sigma[j,i] <- cor.ij
  }
}
Z <- MASS::mvrnorm(n=1,mu=rep(0,50),Sigma=Sigma)
mod <- laGP::newGPsep(X=X,Z=Z,d=2,g=1e-3)

XX <- matrix(runif(100),ncol=2)
system.time(replicate(1000,laGP::predGPsep(gpsepi = mod,XX = XX),))
XX2 <- matrix(runif(2),ncol=2)
system.time(replicate(50000,laGP::predGPsep(gpsepi = mod,XX = XX2),))
laGP::deleteGPsep(mod)

Sample.sizes <- c(2,4,8,16,32,64,128,256,512,1024)
Sample.reps <- 1024*100/Sample.sizes
Sample.times <- rep(-1,10)
for(i in 1:length(Sample.sizes)) {
  XX <- matrix(runif(Sample.sizes[i]),ncol=2)
  S.time <- system.time(replicate(Sample.reps[i],laGP::predGPsep(gpsepi = mod,XX = XX),)) -> temp
  Sample.times[i] <- S.time[3]
}
rbind(Sample.sizes,Sample.times)
laGP::deleteGPsep(mod)
