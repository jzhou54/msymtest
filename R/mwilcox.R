#' @title
#' Modified Wilcoxon test when the center of symmetry is unknown.
#'
#' @description
#' Testing whether a data set comes from a symmetric distribution when the center of symmetry is unknown.
#'
#' @param x data set need to be tested
#'
#' @return
#' @export
#'
#' @examples
#' x = rchisq(100,1)
#' mod.wilcox.test(x)
#'
#'
#'
## 1. modified wilxon test
mod.wilcox.test<- function(x, conf.level = 0.95, ...) {
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  METHOD <- "Modified Wilcoxon signed rank test"
  DNAME <- deparse(substitute(x))
  alternative <- "not symmetric"

  n <- as.double(length(x))
  m <- mean(x)
  sigma2 <- var(x)
  xc <- x-m
  ## Test statistic
  STATISTIC <- setNames(as.numeric(wilcox.test(x,mu=m)$statistic), "W")

  ## Tn selection
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1])/1.34
  Tn<-log(n)/( 3 * 1.06 * min(sqrt(var(x)), h))

  ## Estimation of theta
  S<-function(u) sum( sin(2*pi*(xc[xc!=u]-u)*Tn)/(2*pi*(xc[xc!=u]-u)))+
    sum(sin(2*pi*(xc+u)*Tn)/(2*pi*(xc+u)))
  SV<-Vectorize(S)
  hat_theta <-2*sum(SV(xc))/n^2 + 2*Tn/n

  ## Estimation of tau
  xs <- sort(xc)  # V(i)
  S1 <- seq(from=1, to=n, by=1) # i
  hat_tau <- sum(xs*S1)/n^2

  ## Asymptotic mean and variance
  E <- n*(n+1)/4
  V <- n*(n+1)*(2*n+1)/24 - n*(n-1)*(n-3) * hat_theta * hat_tau +
    (n-1)*(n-2)*(n-3)*(n-4)*sigma2/(4*n)*(hat_theta)^2

  ## The resulting p-value
  pval <- 2*(1-pnorm(abs(E-STATISTIC)/sqrt(V)))

  RVAL <- list(statistic = STATISTIC,
               var = V,
               p.value = as.numeric(pval),
               method = METHOD,
               data.name = DNAME,
               alternative = alternative)
  class(RVAL) <- "htest"
  RVAL
  #return(list("STATISTIC"=STATISTIC, "var"=V, "p.value"=pval))
}
