banana <- function(xx){
  if(is.matrix(xx)) {
    return(apply(xx,1,function(yy){exp(-.005*(yy[1]*40-20)^2-.5*(yy[2]*15-10+.03*(yy[1]*40-20)^2-3)^2)}))
  } else {
    yy <- xx
    exp(-.005*(yy[1]*40-20)^2-.5*(yy[2]*15-10+.03*(yy[1]*40-20)^2-3)^2)
  }
}