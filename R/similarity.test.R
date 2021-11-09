#'
#' @title Compute the similarity test.
#'
#' @description This function compute the similarity test for categorical spatial process.
#'
#' @param data an (optional) data frame or a sf object containing the variable to testing for.
#' @param formula a symbolic description of the factor (optional).
#' @param fx a factor (optional).
#' @param nsim Number of permutations.
#' @param alternative a character string specifying the type of cluster, must be one
#' of "High" (default), "Both" or "Low".
#' @param distr A string. Distribution of the test "asymptotic" (default) or "bootstrap"
#' @param control List of additional control arguments.
#' @usage similarity.test(formula = NULL, data = NULL, na.action, fx = NULL, listw = listw,
#' alternative = "two.sided", distr = "asymptotic", nsim = NULL, control = list())
#' @keywords spatial association, qualitative variable, scan test, clusters
#' @details
#' Escribir aquí cosas
#'
#' @return A object of the \emph{htest}
#'   \tabular{ll}{
#'     \code{data.name} \tab a character string giving the names of the data.\cr
#'     \code{statistic} \tab Value of the similarity test \cr
#'     \code{N} \tab total number of observations.\cr
#'     \code{Zmlc} \tab Elements in the Most Likelihood Cluster.\cr
#'     \code{alternative} \tab a character string describing the alternative hypothesis. \cr
#'     \code{p.value} \tab p-value of the similarity test \cr
#'     \code{similiarity.mc} \tab values of the similarity test in each permutation. \cr
#'     }
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{seedinit} \tab Numerical value for the seed (only for boot version). Default value seedinit=123 \cr
#'       }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Farber, S., Marin, M. R., & Paez, A. (2015).
#'     Testing for spatial independence using similarity relations.
#'       \emph{Geographical Analysis}. 47(2), 97-120.
#'   }
#' @seealso
#' \code{\link{sp.runs.test}}, \code{\link{dgp.spq}}, \code{\link{Q.test}}, , \code{\link{scan.test}}
#' @export
#'
#' @examples
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' listw <- knearneigh(cbind(cx,cy), k = 3)
#' p <- c(1/4,1/4,1/4,1/4)
#' rho <- 0.5
#' fx <- dgp.spq(p = p, listw = listw, rho = rho)
#' W <- (nb2mat(knn2nb(listw)) >0)*1
#' similarity <- similarity.test(fx = fx, data = FastFood.sf, listw = listw)
#' print(similarity)
#'
#' # Case 2: SRQ test with formula, a sf object (points) and knn
#' rm(list = ls())
#' data("FastFood")
#' coor <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' listw <- spdep::knearneigh(coor, k = 4)
#' formula <- ~ Type
#' similarity <- similarity.test(formula = formula, data = FastFood.sf, listw = listw)
#' print(similarity)
#'
#' # Case 3:
#' rm(list = ls())
#' data("Spain")
#' listw <- spdep::poly2nb(as(spain.sf,"Spatial"), queen = FALSE)
#' formula <- ~ MenWoman
#' similarity <- similarity.test(formula = formula, data = spain.sf, listw = listw)
#' print(similarity)


similarity.test <-  function(formula = NULL, data = NULL, na.action, fx = NULL, listw = listw,
                             alternative = "two.sided", distr = "asymptotic", nsim = NULL,
                          control = list()){

# chequeo
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  distr <- match.arg(distr, c("asymptotic", "bootstrap"))

  # Solo admite matrices knn, nb o matrix
  if (class(listw)[1] != "nb"){
    if (class(listw)[1] != "matrix"){
      if (class(listw)[1] != "knn"){
        stop ("listw must be is an object of the class: knn, nb or matrix")
      }
    }
  }
  # Solo admite matrices knn, nb o matrix
  if (class(listw)[1] == "knn"){W <- (nb2mat(knn2nb(listw), zero.policy = TRUE) >0)*1}
  if (class(listw)[1] == "nb"){W <- (nb2mat(listw, zero.policy = TRUE) >0)*1}


# Selecciona los argumentos. Bien con (formula + data) o bien incluye la variable (fx)
if (!is.null(formula) && !is.null(data)) {
  if (inherits(data, "Spatial")) data <- as(data, "sf")
  mfx <- get_all_vars(formula, data)[,1]
  data.name <- names(get_all_vars(formula,data))
} else if (!is.null(fx)) {
  mfx <- fx
  if (is.null(names(fx))) data.name <- "fx"
} else stop("data wrong")

if (!is.factor(mfx))
  stop(paste(deparse(substitute(fx)), "is not a factor"))

  # Controls
  con <- list(seedinit = 123)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  seedinit <- con$seedinit


## Calculo test
N <- length(mfx)
q <- length(levels(mfx))
m <- rep(0,q)
m <- table(mfx)
S0=sum(W)
S1=2*S0;

sw <- rowSums(W)^2
S2 <- 2*sum(sw)
SS <- (matrix(mfx,ncol = N,nrow = N)==t(matrix(mfx,ncol = N,nrow = N)))*W
est <- sum(SS)/N
p0 <- t(m)%*%(m-1)/(N*(N-1))
varij <- p0*(1-p0)
covijk <- (m*(m-1))%*%(m-2)/(N*(N-1)*(N-2))-p0^2
num <- 0
for (c in 1:q){
  m2 <- m[-c]
  num <- num + m[c]*(m[c]-1)* ( (m[c]-2)*(m[c]-3) + t(m2)%*%(m2-1))
}
covijks <-  num/(N*(N-1)*(N-2)*(N-3))-p0^2

# MEAN OF THE ESTATISTIC
meanest <- S0/N*p0

# VARIANCE OF THE ESTATISTIC
varest <- 1/N^2*(S1*varij + (S2-2*S1)*covijk + (S0^2-S2+S1)*covijks)

#####################################
## Similarity mc
#####################################
if (distr == "bootstrap"){
  if (!is.null(seedinit)) set.seed(seedinit)
  similarity.mc <- rep(0,nsim)
  for (f in 1:nsim){
    fxp <- mfx[sample(N)]
    SS <- (matrix(fxp, ncol = N, nrow = N)==t(matrix(fxp, ncol = N, nrow = N)))*W
    estp <- sum(SS)/N
    similarity.mc[f] <- (estp-meanest)/varest^0.5
  }

  SSnormal <- (est-meanest)/varest^0.5
  names(SSnormal) <- "Similarity-test"

  if (alternative == "two.sided")
    p.value <- sum(abs(similarity.mc) > est)/(nsim+1)
  else if (alternative == "greater")
    p.value <- sum(similarity.mc > est)/(nsim+1)
  else p.value <- sum(similarity.mc < est)/(nsim+1)
  if (!is.finite(SSnormal) || p.value < 0 || p.value > 1)
    warning("Out-of-range p-value: reconsider test arguments")

}
if (distr == "asymptotic"){
# NORMAL APPROXIMATION OF THE ESTATISTIC
SSnormal <- (est-meanest)/varest^0.5
names(SSnormal) <- "Similarity-test"
similarity.mc <- NULL
if (alternative == "two.sided")
  p.value <- 2 * pnorm(abs(SSnormal), lower.tail = FALSE)
else if (alternative == "greater")
  p.value <- pnorm(SSnormal, lower.tail = FALSE)
else p.value <- pnorm(SSnormal )
if (!is.finite(SSnormal) || p.value < 0 || p.value > 1)
  warning("Out-of-range p-value: reconsider test arguments")
}

srq <- list(method = paste0("Similarity test of spatial dependence for qualitative data. Distribution: ", distr),
            statistic = SSnormal, p.value = p.value, N = N, alternative = alternative, distr = distr,
            similarity.mc = similarity.mc,
                    data.name = data.name#, estimate = vec
            )
class(srq) <- c("htest")
srq
}
