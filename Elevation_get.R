if (F) {
  par(mfrow=c(2,1))
  plot(1:4)
  my.filled.contour.func(banana)
  
  
  require(XML)
  require(RCurl)
  data_parse <- xmlParse("http://forecast.weather.gov/MapClick.php?lat=29.803&lon=-82.411&FcstType=digitalDWML")
  data_parse <- xmlParse("https://maps.googleapis.com/maps/api/elevation/xml?locations=39.7391536,-104.9847034|36.455556,-116.866667&key=AIzaSyDpHsTEEAWJqE1MbWnMjwJI_Yo54CUviw0")
  xml_data <- xmlToList(data)
  
  elev_url <- "https://maps.googleapis.com/maps/api/elevation/xml?locations=39.7391536,-104.9847034|36.455556,-116.866667&key=AIzaSyDpHsTEEAWJqE1MbWnMjwJI_Yo54CUviw0"
  elev_url <- "https://maps.googleapis.com/maps/api/elevation/xml?locations=39.7391536,-104.9847034|36.455556,-116.866667&key=AIzaSyDpHsTEEAWJqE1MbWnMjwJI_Yo54CUviw0"
  elev_html <- getURL(elev_url)
  elev_xml <- xmlParse(elev_html)
  elevs <- unname(xmlToDataFrame(nodes = getNodeSet((elev_xml),'//result/elevation')))
  resols <- unname(xmlToDataFrame(nodes = getNodeSet((elev_xml),'//result/resolution')))
  lats <- unname(xmlToDataFrame(nodes = getNodeSet((elev_xml),'//result/location/lat')))
  longs <- unname(xmlToDataFrame(nodes = getNodeSet((elev_xml),'//result/location/lng')))
  data.frame(lat=lats,lng=longs,elev=elevs,res=resols)
  #elev_xml[[1]]
  
  elev_xml <- xmlTreeParse(elev_html)
  lat <-  elev_xml[[1]][[1]][[2]][[1]][[1]][[1]]
  long <- elev_xml[[1]][[1]][[2]][[1]][[2]][[1]]
  elev <- elev_xml[[1]][[1]][[2]][[2]][[1]]
  resol <- elev_xml[[1]][[1]][[2]][[3]][[1]]
  #lat <-  XML::xmlTreeParse(elev_html)[[1]][[1]][[2]][[1]][[1]][[1]]
  #long <- XML::xmlTreeParse(elev_html)[[1]][[1]][[2]][[1]][[2]][[1]]
  #elev <- XML::xmlTreeParse(elev_html)[[1]][[1]][[2]][[2]][[1]]
  num_pts <- 2
  for(inum in 1:num_pts) {
    lat <-  as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[1]][[1]][[1]])[6])
    long <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[1]][[2]][[1]])[6])
    elev <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[2]][[1]])[6])
    resol <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[3]][[1]])[6])
    data.frame(lat,long,elev,resol)
  }
}

get.elevation <- function(X) {
  # Limits on requests
  #  2,500 free requests per day
  #  512 locations per request
  #  10 requests per second
  #  Lower resolution for larger requests
  # X: n by 2 matrix
  X <- matrix(round(X,6),ncol=2)
  #X <- X[,c(2,1)]
  #browser()
  # build elevation string
  #elev_url <- "https://maps.googleapis.com/maps/api/elevation/xml?locations=39.7391536,-104.9847034|36.455556,-116.866667&key=AIzaSyDpHsTEEAWJqE1MbWnMjwJI_Yo54CUviw0"
  #X <- matrix(c(39.7391536,-104.9847034,36.455556,-116.866667),ncol=2,byrow=T)
  elev_url <- "https://maps.googleapis.com/maps/api/elevation/xml?locations="
  if(length(X)>2) elev_url <- paste0(elev_url,paste(apply(X,1,function(Xrow)paste(Xrow[2],Xrow[1],sep = ',')),sep=',',collapse = '|'))#"39.7391536,-104.9847034|36.455556,-116.866667")
  else {elev_url <- paste0(elev_url,X[2],',',X[1])}
  elev_url <- paste0(elev_url,"&key=AIzaSyDpHsTEEAWJqE1MbWnMjwJI_Yo54CUviw0")
  
  elev_html <- getURL(elev_url)
  #elev_xml <- xmlTreeParse(elev_html)
  
  #for(inum in 1:(length(X)/2)) {
    #browser()
    #print(elev_xml)
    #browser()
  #  lat <-  as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[1]][[1]][[1]])[6])
  #  long <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[1]][[2]][[1]])[6])
  #  elev <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[2]][[1]])[6])
  #  resol <- as.numeric(as.character(elev_xml[[1]][[1]][[1+inum]][[3]][[1]])[6])
  #  newelevdf <- data.frame(lat,long,elev,resol)
  #  if(inum==1) elev_df <- newelevdf
  #  else elev_df <- rbind(elev_df,newelevdf)
    #newelevdf <- data.frame(lat,long,elev,resol)
  #}
  #browser()
  if (grepl(pattern = 'illegal request',x = elev_html)) {
    stop('Illegal request (probably too long a request)')
  }
  elev_parse <- xmlParse(elev_html)
  status <- xmlToDataFrame(nodes = getNodeSet((elev_parse),'//status'),stringsAsFactors = F)
  if (xmlToDataFrame(nodes = getNodeSet((elev_parse),'//status'))$text[1]=='OK') {
    elevs <- unname(xmlToDataFrame(nodes = getNodeSet((elev_parse),'//result/elevation'),colClasses = 'numeric'))
    resols <- unname(xmlToDataFrame(nodes = getNodeSet((elev_parse),'//result/resolution'),colClasses = 'numeric'))
    lats <- unname(xmlToDataFrame(nodes = getNodeSet((elev_parse),'//result/location/lat'),colClasses = 'numeric'))
    longs <- unname(xmlToDataFrame(nodes = getNodeSet((elev_parse),'//result/location/lng'),colClasses = 'numeric'))
    elev_df <- data.frame(lat=lats,lng=longs,elev=elevs,res=resols)
  } else {
    #browser()
    stop('API request was not successful (status was not OK)')
  }
  return(elev_df)
}
rescale1D <- function(x,xlim,oldxlim=NULL) {
  if (is.null(oldxlim)) {oldxlim <- c(min(x),max(x))}
  (x-oldxlim[1])/(oldxlim[2]-oldxlim[1]) * (xlim[2] - xlim[1]) + xlim[1]
}
rescale2D <- function(x,xlim,ylim,oldxlim=NULL,oldylim=NULL) {
  if (dim(x)[2] != 2) {stop('x must have two columns')}
  x[,1] <- rescale1D(x[,1],xlim,oldxlim)
  x[,2] <- rescale1D(x[,2],ylim,oldylim)
  return(x)
}

if (F) {
  setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('C:/Users/cbe117/School/DOE/SMED/SMED-Code/myfilledcontour.R')
  source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  require(RCurl)
  require(XML)
  StLouisCntyMN <- c(47.882,-92.476)
  SWSanJose <- c(36.504,-123.106)
  NWCol <- c(40.999198, -109.054831)
  SECol <- c(36.996650, -102.056540)
  Colorado.xlim <- c(-109.054831,-102.056540)
  Colorado.ylim <- c(36.996650,40.999198)
  X.LHS <- lhs::maximinLHS(n=40,k=2)
  X.Colorado <- t(apply(X.LHS,1,function(xrow) return(c(xrow[1]*(Colorado.xlim[2]-Colorado.xlim[1])+Colorado.xlim[1],xrow[2]*(Colorado.ylim[2]-Colorado.ylim[1])+Colorado.ylim[1]))))
  X.Colorado2 <- rescale2D(X.LHS,Colorado.xlim,Colorado.ylim,0:1,0:1)
  cbind(X.Colorado,X.Colorado2)
  get.elevation(c(40,-100))
  # Google elevation only takes 92 at a time, 93+ will return an error
  #my.filled.contour.func(get.elevation,n=21,xcontlim=Colorado.xlim,ycontlim=Colorado.ylim,batchmax = 90,out.col.name = 'elev')
  contourfilled.func(get.elevation,n=21,xcontlim=Colorado.xlim,ycontlim=Colorado.ylim,batchmax = 90,out.col.name = 'elev')
  contourfilled.func(get.elevation,n=41,xcontlim=Colorado.xlim,ycontlim=Colorado.ylim,batchmax = 90,out.col.name = 'elev')
}