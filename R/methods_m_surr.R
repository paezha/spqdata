#' @name methods_m_surr
#' @title Method for class m_surr
#' @description A function to plots the m-surround give an object of the
#' class \emph{m_surr} obtain with the code \code{m.surround}.\cr
#'  The \code{plot()} function allows the user view the configuration of the m-surrounding.
#'  The argument \code{type} select the type o visualization. \cr
#'  The \code{print()} print the matrix of the m-surrounding. \cr.
#'  The \code{summary} give information about the characteristics of th m-surrounding \cr.
#'
#' @param x object of class \emph{m_surr}
#' @param type numeric. 1 (default) to get the plot with igraph.
#' @param object object of class \emph{m_surr}.
#' 2 plot W matrix with network
#' @param ... further arguments passed to or from other methods.
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @examples
#'
#' # Obtain m-surroundings with degree of overlapping r
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' m = 4
#' r = 2
#' msurr_points <- m.surround(x = x, m = m, r = r,control = list(dtmaxabs = 0.5))
#' plot(msurr_points, type = 1)
#' plot(msurr_points, type = 2)
#' print(msurr_points)
#'
#' rm(list = ls())
#' data("FastFood")
#' m = 6
#' r = 1
#' msurr_points <-  m.surround(x = FastFood.sf, m = m, r = r, distance = "Euclidean", control = list(dtmaxpc = .01))
#' plot(msurr_points, type = 1)
#' plot(msurr_points, type = 2)
#' print(msurr_points)

#'

NULL

#' @name summary.m_surr
#' @rdname methods_m_surr
#'
#' @export
summary.m_surr <- function(object, ...) {
  z <- object
  N <- z$N
  R <- z$R
  stopifnot(inherits(z, "m_surr"))
  no_symbolized <- (1:N)[!(1:N) %in% unique(matrix(z$ms,ncol = 1))]
  mh <- z$ms
  overlaping <- matrix(0,ncol = R, nrow = R)
  for ( i in 1:R){
    for (j in 1:R){
      overlaping[i,j] <- sum(mh[i,] %in% mh[j,])
    }
  }
  diag(overlaping) <-0

  # Print output
  cat("\nCharacteristics of m-surrounding:\n\n")
  cat(paste0("Number of m-surrounding (R): ",dim(z$ms)[1],"\n"))
  cat(paste0("Length of m-surrounding (m): ",dim(z$ms)[2],"\n"))
  cat(paste0("Number no-symbolized observations: ",length(no_symbolized),"\n"))
  cat("\nList of no-symbolized observations:\n")
  cat(no_symbolized)
  cat("\n\nList of the degree overlaping:\n")
  a <- table(rowSums(overlaping))
  mean.overlaping <- sum(as.numeric(names(a))*a)/sum(a)
  for (i in 1:length(a)){
    cat(paste0("    There are ",a[i], " m-surrounding that have intersection with ",names(a)[i]," m-surrounding","\n"))
  }
  cat(paste0("Mean degree of overlaping: ",round(mean.overlaping,4),"\n"))
}
#'

#' @name plot.m_surr
#' @rdname methods_m_surr
#'
#' @export
plot.m_surr <- function(x, type =1 ,...){

  m_surr <- x
  W <- matrix(0, ncol =  m_surr$N, nrow =  m_surr$N)
  for (i in 1:dim(m_surr$ms)[1]){W[m_surr$ms[i,1],m_surr$ms[i,2:m]] <- 1}
  if (type ==1){
    g1 = igraph::graph.adjacency(W)
    list <- as.list(as.data.frame(t(m_surr$ms)))
    lo <- igraph::layout.norm(st_coordinates(st_centroid(m_surr$x[1])))
    mycolor <- as.matrix(0,ncol = 1, nrow = m_surr$N)
    mycolor <- "red"
    mycolor[m_surr$ms[,1]] = "black"
    plot(g1, margin = 0, edge.width=1,
         vertex.label.font = 1,
         vertex.label.cex = .5,
         vertex.label.dist = 1,
         mark.groups=list,
         edge.arrow.mode=0,layout=lo,vertex.size=3,
         vertex.color = mycolor ,edge.color='black') # ,sub="black points are origin of m-surround"
    title(main=paste("m-surrounding;", "m = ", m_surr$m, "and r = ", m_surr$r,"\n black points are origin of m-surrounding"),cex.main=1,col.main="black")
  }
  if (type ==2){
    W <- mat2listw(W)
    plot(st_geometry(m_surr$x))
    par(new = TRUE)
    plot(W,st_coordinates(st_centroid(m_surr$x)), col = "red")
    title(main=paste("m-surrounding;", "m = ", m_surr$m, "and r = ", m_surr$r))
  }
}

#' @name print.m_surr
#' @rdname methods_m_surr
#' @export
#'
print.m_surr <- function(x,...){
print(x$ms)
}



