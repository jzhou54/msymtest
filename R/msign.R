#' @title
#' Modified sign test when the center of symmetry is unknown
#'
#' @description
#' Testing whether a data set comes from a symmetric distribution when the center of symmetry is unknown.
#'
#' @param x
#' data set need to be tested
#'
#' @return
#' @export
#'
#' @examples
#' x = rchisq(100,1)
#' mod.sign.test(x)
#'
#'

## modified sign test
mod.sign.test <- function(x) {
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  METHOD <- "Modified sign test"
  DNAME <- deparse(substitute(x))
  alternative <- "not symmetric"

  n <- as.double(length(x))
  m <- mean(x)

  ## Test statistic
  STATISTIC <- setNames(sum(1*(x<m)),"S")

  ## Estimation of w
  a <- 1 * (x>m-n^(-1/5)) * (x < m+n^(-1/5))
  D <- sum(a)
  A <- max(c(1,D))
  VHW <- (n^(3/10)/A)^2
  hat_w <- sqrt(1/(4*n*VHW))

  ## probability weight moment for x
  CE <- mean((x-m)*(x<m))

  ## Asymptotic variance
  V <- 1/4 + var(x)*(hat_w)^2 + 2*hat_w*CE

  # The resulting p-value
  pval <- 2*(1-pnorm(abs(STATISTIC-n/2)/sqrt(n*V)))

  RVAL <- list(statistic = STATISTIC,
               var = V,
               p.value = as.numeric(pval),
               method = METHOD,
               data.name = DNAME,
               alternative = alternative)
  class(RVAL) <- "htest"
  RVAL
}
