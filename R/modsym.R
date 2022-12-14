#' @title
#' Modified symmetric test
#'
#' @description
#' Provides one symmetric test that incorporated 2 methods, "W" and "S". Also, paired data can be tested.
#'
#' @param x data set to be tested
#' @param y another data set
#' @param paired logical. whether data set x and y are paired
#' @param alternative
#' a character string specifying the alternative hypothesis, must be one of "two sided" (default), "greater", or "less".
#' greater: test whether positively skewed (right-skewed),
#' less : test whether negatively skewed (left-skewed)
#' @param method a character string specifying which symmetric test to be used, "W" is modified wilcoxon sign rank test,
#' and "S" is modified sign test.
#'
#' @return
#' @export
#'
#' @examples
#' x = rchisq(100,1)
#' y = rchisq(100,1)
#' mod.sym.test(x, method="wilcox")
#' mod.sym.test(x, method="sign", alternative="greater")
#'
#' mod.sym.test(x, y, method="wilcox")
#'
mod.sym.test <- function(x, y=NULL,
           alternative = c("two.sided", "less", "greater"),
           method = "W") {
    alternative <- match.arg(alternative)
    alt.text <- switch (alternative,
                        "two.sided" = "not symmetric",
                        "less" = "left-skewed",
                        "greater" = "right-skewed" )
    if(!is.numeric(x)) stop("'x' must be numeric")
    if(!is.null(y)) {
      if(!is.numeric(y)) stop("'y' must be numeric")
      DNAME <- paste(deparse(substitute(x)), "and",
                     deparse(substitute(y)))
      
      if(length(x) != length(y))
         stop("'x' and 'y' must have the same length")
       OK <- complete.cases(x, y)
       x <- x[OK] - y[OK]
       y <- NULL
     
   }else {
      DNAME <- deparse(substitute(x))
      x <- x[is.finite(x)]
   }

    if(length(x) < 1L)
      stop("not enough (finite) 'x' observations")
    n <- as.double(length(x))
    m <- mean(x)
    sigma2 <- var(x)
    xc <- x-m


    if (method == "wilcox"){
      METHOD <- "Modified wilcoxon signed rank test"

      ## Test statistic
      STATISTIC <- setNames(as.numeric(wilcox.test(x,mu=m)$statistic), "W")

      ## Tn selection
      r <- quantile(x, c(0.25, 0.75))
      h <- (r[2] - r[1])/1.34
      Tn<-log(n)/( 3 * 1.06 * min(sqrt(var(x)), h))

      ## Estimation of theta
      S <- function(u) sum( sin(2*pi*(xc[xc!=u]-u)*Tn)/(2*pi*(xc[xc!=u]-u)))+
        sum(sin(2*pi*(xc+u)*Tn)/(2*pi*(xc+u)))
      SV <- Vectorize(S)
      hat_theta <- 2*sum(SV(xc))/n^2 + 2*Tn/n

      ## Estimation of tau
      xs <- sort(xc)  # V(i)
      S1 <- seq(from=1, to=n, by=1) # i
      hat_tau <- sum(xs*S1)/n^2

      ## Asymptotic mean and variance
      E <- n*(n+1)/4
      V <- n*(n+1)*(2*n+1)/24 - n*(n-1)*(n-3) * hat_theta * hat_tau +
        (n-1)*(n-2)*(n-3)*(n-4)*sigma2/(4*n)*(hat_theta)^2

      ## The resulting p-value
      pval <- switch (alternative,
                      "two.sided" =  2 * ( 1 - pnorm(abs(E-STATISTIC)/sqrt(V))),
                      "greater" = 1 - pnorm( (E -STATISTIC)/sqrt(V) ),
                      "less" = 1 - pnorm(-(E -STATISTIC)/sqrt(V))
                      )

    } else if (method == "sign") {
      METHOD <- "Modified sign test"

      ## Test statistic
      STATISTIC <- setNames(sum(1*(x<m)), "S")

      ## Estimation of w
      a <- 1 * (x>m-n^(-1/5)) * (x < m+n^(-1/5))
      D <- sum(a)
      A <- max(c(1,D))
      VHW <- (n^(3/10)/A)^2
      hat_w <- sqrt(1/(4*n*VHW))
      E <- n/2

      ## probability weight moment for x
      CE <- mean((x-m)*(x<m))

      ## Asymptotic variance
      V <- 1/4 + var(x)*(hat_w)^2 + 2 *hat_w*CE

      pval <- switch (alternative,
                       "two.sided" =  2 * ( 1 - pnorm(abs(STATISTIC-E)/sqrt(n*V))),
                       "greater" = 1 - pnorm( (STATISTIC-E)/sqrt(n*V) ),
                       "less" = 1 - pnorm(-(STATISTIC-E)/sqrt(n*V))
                      )

    }


    pval <- setNames(pval, "p.value")
    # pval <- 2*(1-pnorm(abs(STATISTIC-n/2)/sqrt(n*V)))
    RVAL <- list(statistic = STATISTIC,
                 var = V,
                 p.value = as.numeric(pval),
                 method = METHOD,
                 data.name = DNAME,
                 alternative = alt.text)
    class(RVAL) <- "htest"
    RVAL


  }

