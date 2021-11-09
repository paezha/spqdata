scan.bernoulli <- function(xf = NULL, nv = NULL, coor = NULL){
  N <- length(xf)
  cx <- coor[,1]
  cy <- coor[,2]
  nv <- trunc(N/2)
  nn <- cbind(1:N,spdep::knearneigh(cbind(cx,cy), k = nv)$nn)
  XF <- matrix(xf[nn], ncol = (nv+1), nrow = N)
  oz <- t(apply(XF==1, 1 , cumsum))
  nz <- t(matrix(rep(1:(nv+1),N),nrow = (nv+1)))
  O <- sum(xf==1)
  a <- log(oz/nz)
  a[a==-Inf]<-0
  b <- log(1-oz/nz)
  b[b==-Inf]<-0
  c <- log((O-oz)/(N-nz))
  c[c==-Inf]<-0
  d <- log(1-(O-oz)/(N-nz))
  d[d==-Inf]<-0
  lnlz <- exp(oz*a + (nz-oz)*b + (O-oz)*c + (N-O-nz+oz)*d)
  # lnlz <- oz*log(oz/nz) + (nz-oz)*log(1-oz/nz) + (O-oz)*log((O-oz)/(N-nz)) + (N-O-nz+oz)*log(1-(O-oz)/(N-nz))
  lnlz0 <- exp(O*log(O) + (N-O)*log(N-O) - N*log(N))
  lnlz <- log(lnlz/lnlz0)
  a <- which(lnlz == max(lnlz), arr.ind = TRUE)[1,]
  Zmlc <- nn[a[1],1:a[2]]
  loglik <- max(lnlz)
  cases <- sum(xf[Zmlc]==1)
  expec.cases <- a[2]*(O/N)

}
