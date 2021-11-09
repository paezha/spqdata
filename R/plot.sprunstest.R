#' @name plot.sprunstest
#' @rdname plot.sprunstest
#'
#' @title Plot the empirical distribution of runs
#'
#' @details Plot the histogram with the empirical distribution
#' of the runs
#' @usage
#' ## S3 method for class 'sprunstest'
#' plot.sprunstest(srq = srq)
#' @param srq  A object of class \emph{sprunstest}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paez@@gmail.com} \cr
#'   Manuel Ruiz  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'
#' @seealso
#' \code{\link{sp.runs.test}}.
#'
#'
#' @examples
#' # Fastfood example. sf (points)
#' data("FastFood")
#' x <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' listw <- spdep::knearneigh(x, k = 2)
#' formula <- ~ Type
#' srq <- sp.runs.test(formula = formula, data = FastFood.sf, listw = listw, nsim = 299)
#' plot(srq)
#'
#' # Spain example (poligons with 0 neinghbourhood)
#' data("Spain")
#' listw <- spdep::poly2nb(as(spain.sf,"Spatial"), queen = FALSE)
#' formula <- ~ Older65
#' srq <- sp.runs.test(formula = formula, data = spain.sf, listw = listw, nsim = 299)
#' plot(srq)
#' formula <- ~ MenWoman
#' srq <- sp.runs.test(formula = formula, data = spain.sf, listw = listw, nsim = 299)
#' plot(srq)

#' @export
#'

plot.sprunstest <- function(srq = srq,...){

#####################
### Plot Q Global
#####################

if (is.null(srq$nsim)){
  runs <- srq$MaxNeig # max(as.numeric(levels(as.data.frame(srq$dnr)$Var1)))
  xxx <- as.data.frame(srq$dnr)
  # xxx$Freqr <- as.data.frame(srq$dnr)$Freq/sum(as.data.frame(srq$dnr)$Freq)
  ggplot(data=xxx, aes(x=Var1, y = Freq)) +
    geom_bar(stat="identity",color = "black", fill = "steelblue") +
    labs(x = "Number of runs", y = "Frequency") +
    theme_bw()

}
  else {
runs <- srq$MaxNeig # max(as.numeric(levels(as.data.frame(srq$dnr)$Var1)))
xxx <- as.data.frame(srq$dnr)
aa <- matrix(0,ncol = runs, nrow = srq$nsim)
for (i in 1:srq$nsim){
aa[i,] <- as.data.frame(table(factor(srq$SRLP[,i], levels = c(1:runs))))$Freq
}
hh <- matrix(0, ncol = 2, nrow = runs)
for (i in 1:runs){
hh[i,] <- quantile(aa[,i],c(0.05,.95))
}
xxx$min <- hh[,1]
xxx$max <- hh[,2]
xxx$mean <- colMeans(aa)

ggplot(data=xxx, aes(x=Var1, y = Freq)) +
  geom_bar(stat="identity",color = "black", fill = "steelblue") +
   labs(x = "Number of runs", y = "Frequency") +
   geom_errorbar(data = xxx, aes(x =  Var1, ymin = min,ymax = max), width=0.3, colour="red", alpha=0.9, size=1.1) +
   geom_point(data=xxx, aes(x = Var1, y = mean), size = 3, shape = 21, fill = "white") +
   theme_bw()
  }
}





