#' @title Generation of qualitative process with spatial structure
#' @description The purpose of the function \code{dgp.spq} is to generate a random dataset
#' with the dimensions and spatial structure decided by the user. This function may be useful in pure
#' simulation experiments or with the aim of showing specific properties and characteristics of a
#' spatial qualitative dataset ...
#' @usage dgp.spq(listw = listw, p = p,  rho = rho, control = list())
#' @param listw A \code{listw} object of the class nb, knn, listw o matrix created for example by
#'   \code{\link[spdep]{nb2listw}} from \pkg{spatialreg} package; if
#'   \code{\link[spdep]{nb2listw}} not given, set to
#'   the same spatial weights as the \code{listw} argument. It can
#'   also be a spatial weighting matrix of order \emph{(NxN)} instead of
#'   a \code{listw} object. Default = \code{NULL}.
#' @param rho the level of spatial dependence (values between -1 y 1)
#' @param p a vector with the percentage of elements of each categories. The lengths must be the number of categories.
#' The sum of the elements of vector must be 1.
#' @param control List of additional control arguments. See control argument section.
#' @return a factor of length N with levels the fist natural numbers.
#' @details Aquí Antonio escribe una linda historia ...
#'
#' La forma de generar datos es la descrita en Páez et al. 2010 (pag 291)
#'
#' $$ Y = (I- rho W)^{-1} epsilon $$
#'
#' where $epsilon$ = N(0,1) and where I is the N x N identity matrix, q is a parameter of spatial dependence,
#' and W is a connectivity matrix that determines the set of spatial relationships among points.
#' the continuous spatially autocorrelated variable Y is used to define a discrete spatial process
#' as follows. Let $b_{ij}$ be defined by ...
#' @section Control arguments:
#' \describe{
#' \item{seedinit}{seed to generate the data sets}
#' }
#' @seealso
#' \code{\link{qtest}}, \code{\link{sp.runs.test}}, \code{\link{m_surr_no}}
#' @keywords m_surround, q-test, spatial run test
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Ruiz M, López FA, A Páez. (2010). \emph{Testing for spatial association of qualitative
#'     data using symbolic dynamics}. Journal of Geographical Systems. 12 (3) 281-309
#'   }
#' @export
#' @examples
#' #
#' rm(list = ls())
#' N <- 10
#' cx <- runif(N)
#' cy <- runif(N)
#' coor <- cbind(cx,cy)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' listw <- spdep::nb2listw(knn2nb(knearneigh(cbind(cx,cy), k = 4)))
#' xf <- dgp.spq(list = listw, p = p, rho = rho)
#'
#' rm(list = ls())
#' data(Spain)
#' listw <- spdep::poly2nb(spain.sf, queen = FALSE)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.9
#' xf <- dgp.spq(p = p, listw = listw, rho = rho)
#' spain.sf$xf <- xf
#' plot(spain.sf["xf"])

dgp.spq <- function(listw = listw, p = p,  rho = rho, control = list()) {

  ################################################
  # Controls
  con <- list(seedinit = NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  seedinit <- con$seedinit

  if (sum(p)!=1)  stop("The sum of p must be equal to 1")

if (class(listw)[1]=="knn"){listw <- spdep::nb2listw(knn2nb(listw))}
if (class(listw)[1]=="listw"){listw <- spdep::listw2mat(listw)}
if (class(listw)[1]=="nb"){listw <- spdep::nb2mat(listw, style = "W",zero.policy = TRUE)}
if (class(listw)[1]=="matrix") {
  listw <- listw/matrix(rowSums(listw),ncol = dim(listw)[1],nrow =dim(listw)[1])
  listw[is.na(listw)]<-0
  }

set.seed(seedinit)

n <- dim(listw)[1]
listw <- as(listw,"dgCMatrix")
k = 1
y <- Matrix::solve(Matrix::Diagonal(n)-rho*listw)%*%matrix(rnorm(n*k,mean = 0, sd = 1), n, k)  # y <- Matrix::solve(diag(n)-rho*listw)%*%rnorm(n,1)
y <- as.matrix(y)
Y <- cut(y,quantile(y,c(0,cumsum(p))),include.lowest=TRUE)
levels(Y) <- as.character(1:length(p))
levels(Y) <- LETTERS[1:length(p)]

return(Y)
}
