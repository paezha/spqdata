#'
#' @title Compute the global spatial runs test.
#'
#' @description This function compute the global spatial runs test for spatial independence of a
#' categorical spatial data set.
#' @param data an (optional) data frame or a sf object containing the variable to testing for.
#' @param listw una lista de vecinos (tipo knn o nb) o una matrix W que indique el orden de cada $m_i-entorno$
#' (por ejemplo de inversa distancia). Para calcular el número de rachas en cada m_i-entorno debe establecerse
#' un orden, por ejemplo del vecino más próximo al más lejano.
#' @param fx a factor (optional).
#' @param formula a symbolic description of the factor (optional).
#' @param alternative a character string specifying the alternative hypothesis, must be one
#' of "two.sided" (default), "greater" or "less".
#' @param distr A string. Distribution of the test "asymptotic" (default) or "bootstrap"
#' @param nsim Number of permutations to obtain confidence intervals (CI).
#' Default value is NULL to don`t get CI of number of runs.
#' @param control List of additional control arguments.
#' @usage sp.runs.test(formula = NULL, data = NULL, na.action, fx = NULL,
#' listw = listw, alternative = "two.sided" , regular = FALSE,
#' distr = "asymptotic", nsim = NULL,control = list())
#' @keywords spatial association, qualitative variable, runs test
#' @details El orden de las vecindades ($m_i-entornos$) es crítico. \cr
#' Para obetener el número de rachas observadas en cada $m_i-entorno$,
#' cada elemento debe asociarse a un conjunto de vecinos ordenados por proximidad.
#' Tres clases de listas pueden incluirse para identificar $m_i-entornos$:
#'
#' \tabular{ll}{
#'     \code{knn} \tab Matrices de la clase knn que consideran los vecinos por orden de proximidad.
#'     Ver \code{\link{knn2knn_order}}\cr
#'     \code{nb} \tab Si los vecinos se obtienen a partir de un objeto sf, el código internamente llamará
#'     a la función \code{\link{nb2nb_order}} los ordenará en orden de proximidad de los centroides.  \cr
#'     \code{matrix} \tab Si se introduce simplemente una matriz basasa en la inversa de la distancia,
#'     también se llamará internamente a la función \code{\link{nb2nb_order}} para transformar la matriz
#'     de la clase matrix a una matriz de la clase nb con vecinos ordenados. \cr
#'     }
#'
#' Two alternative sets of arguments can be included in this function to compute the spatial runs test:
#'
#'   \tabular{ll}{
#'     \code{Option 1} \tab A factor (fx) and a list of neighborhood (\code{listw}) of the class knn. \cr
#'     \code{Option 2} \tab A sf object (data) and formula to specify the factor. A list of neighborhood (listw) \cr
#'     }
#'
#' @return A object of the \emph{htest} and \emph{sprunstest} class
#'   \tabular{ll}{
#'     \code{data.name} \tab a character string giving the names of the data.\cr
#'     \code{method} \tab the type of test applied ().\cr
#'     \code{SR} \tab total number of runs \cr
#'     \code{dnr} \tab empirical distribution of the number of runs \cr
#'     \code{statistic} \tab Value of the homogeneity runs statistic. Negative sign indicates global homogeneity \cr
#'     \code{alternative} \tab a character string describing the alternative hypothesis. \cr
#'     \code{p.value} \tab p-value of the SRQ \cr
#'     \code{pseudo.value} \tab the pseudo p-value of the SRQ test if nsim is not NULL\cr
#'     \code{MeanNeig} \tab Mean of the Maximum number of neighborhood \cr
#'     \code{MaxNeig} \tab Maximum number of neighborhood \cr
#'     \code{listw} \tab The list of neighborhood \cr
#'     \code{nsim} \tab number of boots (only for boots version) \cr
#'     \code{SRGP} \tab nsim simulated values of statistic. \cr
#'     \code{SRLP} \tab matrix with the number of runs for eacl localization. \cr
#'     }
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{seedinit} \tab Numerical value for the seed (only for boot version). Default value seedinit=123 \cr
#'       }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @seealso
#' \code{\link{local.sp.runs.test}}, \code{\link{dgp.spq}}, \code{\link{Q.test}},
#' @export
#'
#' @examples
#'
#' # Case 1: SRQ test based on factor and knn
#' rm(list = ls())
#' n <- 100
#' cx <- runif(n)
#' cy <- runif(n)
#' x <- cbind(cx,cy)
#' listw <- knearneigh(cbind(cx,cy), k=3)
#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' fx <- dgp.spq(listw = listw, p = p, rho = rho)
#' srq <- sp.runs.test(fx = fx, listw = listw)
#' print(srq)
#' plot(srq)
#' # Version boots
#' control <- list(seedinit = 1255)
#' srq <- sp.runs.test(fx = fx, listw = listw, distr = "bootstrap" , nsim = 299, control = control)
#' print(srq)
#' plot(srq)
#'
#' # Case 2: SRQ test with formula, a sf object (points) and knn
#' rm(list = ls())
#' data("FastFood")
#' x <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' listw <- spdep::knearneigh(x, k=4)
#' formula <- ~ Type
#' srq <- sp.runs.test(formula = formula, data = FastFood.sf, listw = listw)
#' print(srq)
#' plot(srq)
#' # Version boots
#' srq <- sp.runs.test(formula = formula, data = FastFood.sf, listw = listw, distr = "bootstrap", nsim = 199)
#' print(srq)
#' plot(srq)
#'
#' # Case 3: SRQ test (permutation) using formula with a sf object (polygons) and nb
#' rm(list = ls())
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' listw <- spdep::poly2nb(as(nc,"Spatial"), queen = FALSE)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' co <- sf::st_coordinates(sf::st_centroid(nc))
#' nc$fx <- dgp.spq(listw = listw, p = p, rho = rho)
#' plot(nc["fx"])
#' formula <- ~ fx
#' srq <- sp.runs.test(formula = formula, data = nc, listw = listw, distr = "bootstrap", nsim = 399)
#' print(srq)
#' plot(srq)
#'
#' # Case 4: SRQ test (Asymptotic) using formula with a sf object (polygons) and nb
#' rm(list = ls())
#' # PARA PROBAR ELEMENTOS SIN VECINOS
#' data(Spain)
#' listw <- spdep::poly2nb(spain.sf, queen = FALSE)
#' plot(spain.sf["Coast"])
#' formula <- ~ Coast
#' srq <- sp.runs.test(formula = formula, data = spain.sf, listw = listw)
#' print(srq)
#' plot(srq)
#' # Version boots
#' srq <- sp.runs.test(formula = formula, data = spain.sf, listw = listw, distr = "bootstrap", nsim = 299)
#' print(srq)
#' plot(srq)
#'
#' # Case 5: SRQ test based on a distance matrix (inverse distance)
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' data <- as.data.frame(cbind(cx,cy))
#' data <- st_as_sf(data,coords = c("cx","cy"))
#' n = dim(data)[1]
#' dis <- 1/matrix(as.numeric(st_distance(data,data)),ncol=n,nrow=n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.10))*dis
#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' fx <- dgp.spq(listw = dis , p = p, rho = rho)
#' srq <- sp.runs.test(fx = fx, listw = dis)
#' print(srq)
#' plot(srq)
#'
#' srq <- sp.runs.test(fx = fx, listw = dis, data = data)
#' print(srq)
#' plot(srq)
#' # Version boots
#' srq <- sp.runs.test(fx = fx, listw = dis, data = data, distr = "bootstrap", nsim = 299)
#' print(srq)
#' plot(srq)
#'
#' # Case 6: SRQ test based on a distance matrix (inverse distance)
#' rm(list = ls())
#' data("FastFood")
#' n = dim(FastFood.sf)[1]
#' dis <- 1000000/matrix(as.numeric(st_distance(FastFood.sf,FastFood.sf)), ncol = n, nrow = n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.005))*dis
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' co <- sf::st_coordinates(sf::st_centroid(FastFood.sf))
#' FastFood.sf$fx <- dgp.spq(p = p, listw = dis, rho = rho)
#' plot(FastFood.sf["fx"])
#' formula <- ~ fx
#' # Version boots
#' srq <- sp.runs.test(formula = formula, data = FastFood.sf, listw = dis, distr = "bootstrap", nsim = 299)
#' print(srq)
#' plot(srq)


sp.runs.test <-  function(formula = NULL, data = NULL, na.action, fx = NULL,
                       listw = listw, alternative = "two.sided" , regular = FALSE, distr = "asymptotic", nsim = NULL,
                      control = list()){

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  distr <- match.arg(distr, c("asymptotic", "bootstrap"))

  # Solo admite matrices knn, nb o matrix
  if (class(listw)[1] != "knn"){
    if (class(listw)[1] != "nb"){
      if (class(listw)[1] != "matrix"){
      stop ("listw must be is an object of the class: knn, nb or matrix")
      }
    }
  }

# Controls
  con <- list(seedinit = 123)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  seedinit <- con$seedinit
################################################
## Tratamiento de la matrix
################################################

# Si se trata de un objeto sf y la matrix es tipo 'nb' / "matrix" hay que ordenar los m_i-entornos
if (sum(class(data)[1]=="sf")==1){
if (class(listw)[1]=='nb'){ # hay que ordenar los elementos
  listw <- nb2nb_order(listw=listw, sf = data)
}
if (class(listw)[1]=='matrix'){ # hay que ordenar los elementos
      listw <- mat2listw(listw)$neighbours
      class(listw) <- "nb"
      listw <- nb2nb_order(listw=listw, sf = data)
}
}

if (sum(class(data)[1]=="sf")==0){
  if (class(listw)[1]=='matrix'){ # hay que ordenar los elementos
    listw <- nb2nb_order(listw=listw)
  }
}
################################################
## Tratamiento del input de los datos
################################################
# Selecciona los argumentos. Bien con (formula + data) o bien incluye la variable (fx)
  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mxf <- get_all_vars(formula, data)
  } else if (!is.null(fx)) {
    mxf <- fx
    # if (!is.matrix(mxf)) mxf <- as.matrix(mxf, ncol = 1)
    mxf <- as.data.frame(mxf)
    for (i in 1:ncol(mxf)) {
      if (!is.factor(mxf[,i])) mxf[,i] <- as.factor(mxf[,i])
    }
  } else stop("data wrong")

  # fx debe ser un factor. Lo transformo en var numerica para calcular
  if (is.factor(mxf[,1])){
    levels(mxf[,1]) <- as.character(1:length(levels(mxf[,1])))
    y <- as.numeric(mxf[,1])
  }
  if (is.character(mxf[,1])){
    y <- as.factor(mxf[,1])
    levels(fx[,1]) <- as.character(1:length(levels(mxf[,1])))
    y <- as.numeric(mxf[,1])
  }
  if (is.numeric(mxf[,1])){
    stop("Only factors are admitted")
  }

################################################
## Empezamos las cuentas del test
################################################

# Calculo valores previos para obtener media y varianza estadístico
nv <- creation_nvar_SR(listw = listw)
q <- max(y)
n <- length(y)
# Cont is a binary variable that takes on the value of 1 if data are
# continuous and 0 if data are categorical.

m <- numeric() # matrix(0,nrow=q,ncol=1)
pprod <- numeric() # matrix(0,nrow=q,ncol=1)
for (i in 1:q){
  m[i] <- sum(y==i)
  pprod[i]<- m[i]*(n-m[i])
}
if (class(listw)[1]=="knn"){
lnnb <- matrix(dim(listw$nn)[2],ncol = 1,nrow = dim(listw$nn)[1])}
if (class(listw)[1]=="nb"){
lnnb <- rowSums(nb2mat(listw,style='B',zero.policy = TRUE))
}

MaxNeig <- max(lnnb)+1 # El elemento en cuestion + sus vecinos

# here we categorize the original data set y into the q categories
# compute the m_k needed for the computation of mean and variance
# pprod is needed for the computation of p
p <- sum(pprod)/(n*(n-1))


##### COMPUTING THE VARIANCE #####
## case 1 ##
aux1 <- numeric()
aux31 <- numeric()
aux3 <- numeric()

t1=0;
for (k in 1:q){
  for (c in 1:q){
    t1=t1+1
    aux1[t1]=m[k]*m[c]*(n-m[c]-1)
    aux31[t1]=m[k]*m[c]*((m[k]-1)*(n-m[k]-1)+(m[c]-1)*(n-m[c]-1))
    if(k==c){
      aux1[t1]=0
      aux31[t1]=0
    }
  }
}

t3=0
aux3<-numeric()
for (k in 1:q){
  for (c in 1:q){
    for (d in 1:q){
      t3=t3+1
      aux3[t3]=m[k]*m[c]*m[d]*(n-m[d]-2)
      if (c==k){aux3[t3]=0}
      if (d==k){aux3[t3]=0}
      if (d==c){aux3[t3]=0}
    }
  }
}

var1=1/(n*(n-1)*(n-2)*(n-3))*(sum(aux3)+sum(aux31));
var2=1/(n*(n-1)*(n-2))*sum(aux1);
var3=p;

varSR=p*(1-p)*sum(lnnb)+nv[1]*var1+nv[2]*var2+nv[3]*var3-(nv[1]+nv[2]+nv[3])*p^2

# Here we compute the runs starting at each location and it sum is the total number of runs
nruns <- matrix(0,ncol = 1,nrow = n)
for (i in 1:n){
  if (lnnb[i]!= 0){ # Solo calcula los test locales si el elemento tiene vecinos
  if (class(listw)[1]=="knn"){
    runs <- y[c(i,listw$nn[i,])]}
  if (class(listw)[1]=="nb"){
    runs <- y[c(i,listw[[i]])]}
nruns[i] <- 1 + sum(abs(diff(runs))>0)
  }
}
# The distribution of the number of runs
dnr <- table(factor(nruns, levels = c(1:MaxNeig))) # dnr <- table(SRQlocal[,1],exclude = 0) # Excluimos localizaciones con 0 vecinos


SR <- sum(nruns)

# The mean of the statistic
meanSR <- n+p*sum(lnnb)
vec <- c(SR,meanSR,varSR)
names(vec) <- c("Total runs","Mean total runs","Variance total runs")
# The SRQ global test statistic which is N(0,1) distributed
SRQ <- (SR-meanSR)/sqrt(varSR)
if (alternative =="two.sided"){
p.value <- 2*(1-pnorm(abs(SRQ), mean = 0, sd = 1))
} else if (alternative =="less"){
p.value <- pnorm(SRQ, mean = 0, sd = 1)
} else if (alternative =="greater"){
p.value <- 1-pnorm(SRQ, mean = 0, sd = 1)
}

############################################################################
# Para la obtención de los intervalos de confianza por boots permutacional
############################################################################
if ((is.null(nsim) == FALSE) && (distr != "asymptotic")){
  if (!is.null(seedinit)) set.seed(seedinit)
SRGP <- matrix(0,ncol = 1,nrow = nsim)
SRLP <- matrix(0,ncol = nsim, nrow = n)
    for (i in 1:nsim){
    yp <- y[sample(1:n)]
    srqp <- SR_test_boots(fx = yp, listw = listw, nv = nv)
    SRGP[i] <- srqp$SRglobal
    SRLP[,i] <- srqp$nruns
    }

vec <- c(SR,mean(colSums(SRLP)),sd(colSums(SRLP))^2)
names(vec) <- c("Observed Total runs","Mean total runs boots","Variance total runs boots")

if (alternative =="greater"){
pseudo.value <- sum(SRGP > SRQ)/(nsim+1)
}
else if (alternative =="less"){
  pseudo.value <- sum(SRGP < SRQ)/(nsim+1)
}
else if (alternative =="two.sided"){
  pseudo.value <- (sum(SRGP < -abs(SRQ)) + sum(SRGP > abs(SRQ)))/(nsim+1)
  # p.valueSRQB <- (sum(SRGP < SRQ) + sum(SRGP > SRQ))/(nsim+1)
} else stop("alternative wrong")
}

# El número total de rachas es la suma de las rachas en cada m_i-entorno
# colSums(SRLP) da el número de rachas de las nsim repeticiones
# SRGP es un vector con los valores de SRGlobal bajo aleatoriedad
############################################################
# Salida en función si se piden intervalos IC o no
names(SRQ) <- "sp.runs test"
data.name <- names(mxf)
if ((is.null(nsim) == FALSE) && (distr != "asymptotic")){
srq <- list(method = "Runs test of spatial dependence for qualitative data. Boots version",
               SR = SR, dnr = dnr, statistic  = SRQ, alternative = alternative,
               SRGP = SRGP, p.value = pseudo.value,
               data.name = data.name, estimate = vec,
               nsim = nsim, SRLP = SRLP, MeanNeig=sum(lnnb)/n,
               listw = listw, MaxNeig = MaxNeig)

}
else
{
srq <- list(method = "Runs test of spatial dependence for qualitative data. Asymptotic version",
               SR = SR, dnr = dnr, statistic = SRQ, p.value = p.value, alternative = alternative,
               data.name = data.name, estimate = vec,
               MeanNeig=sum(lnnb)/n,weights= listw, MaxNeig = MaxNeig)


}
class(srq) <- c("htest","sprunstest")
srq
}
