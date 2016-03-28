#f <- function(xx){apply(xx,1,function(yy){exp(-.005*yy[1]-.5*(yy[2]+.03*yy[1]^2-3)^2)})}
f <- function(xx){apply(xx,1,function(yy){exp(-.005*(yy[1]*40-20)^2-.5*(yy[2]*15-10+.03*(yy[1]*40-20)^2-3)^2)})}

# Get contour plot
nc <- 100 # number contour each dimension
fx <- fy <- seq(0,1,length.out = nc)
fz <- matrix(0,nc,nc)
for (xi in 1:nc) for(yi in 1:nc) fz[xi,yi] <- f(matrix(c(fx[xi],fy[yi]),1,2))
contour(fx,fy,fz,nlevels = 5)
filled.contour(fx,fy,fz,nlevels = 5)


qq <- function(xx){f(xx)^-.5}
f_min <- Vectorize(function(xnew,xall,k=4) {
  #min(
  qq(xnew)^k*sum(sapply(xall,function(xx){(qq(xx)/abs(xx-xnew))^k}))
  #,1e100) # Thought Inf was a problem, still might be
}, vectorize.args='xnew')
curve(f_min(x,c(.75001)))

# Initialize with mode
curve(f)
#X <- c(.75)
xavail <- runif(1000)
xmodeind <- which.max(f(xavail))
X <- xavail[xmodeind]
text(X,f(X),1,col='red')
xavail <- xavail[-xmodeind]
for(i in 2:30) {
  #f_opt <- optimize(function(xx){f_min(xx,X)},interval = c(0,1),temp='SANN')
  #xnew <- f_opt$min
  f_opt <- which.min(f_min(xavail,X))
  xnew <- xavail[f_opt]
  X <- c(X,xnew)
  xavail <- xavail[-f_opt]
  text(xnew,f(xnew),i,col='red')
  points(xnew,0,col=2)
}
#curve(f_min(x,X),ylim=c(1e6,1e11))
hist(X,freq=F)
curve(f,add=T,col='red')