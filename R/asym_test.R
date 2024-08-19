#' Large Sample Approximation (Asymptotic) Test of Correlation Coefficients
#'
#' This function performs a large sample approximation test of correlation coefficients, ensuring control over
#' type I error under general scenarios when the sample size exceeds 200. It is suitable for cases where two
#' variables are dependent but uncorrelated.
#' @details
#' #' The test supports the following correlation coefficients: Pearson correlation coefficient, Weighted
#' Pearson correlation coefficient, Spearman correlation coefficient, and Lin's concordance correlation coefficient (CCC)
#'
#' For Pearson, weighted Pearson, and Spearman correlation coefficients, the test supports a zero null
#' hypothesis. The alternative hypothesis can be either one-sided or two-sided.
#'
#' For Lin's concordance correlation coefficient (CCC), the test accommodates a more general null hypothesis.
#' Currently, the test only supports a one-sided alternative hypothesis (greater).
#'
#' @param x a \code{numeric} vector.
#' @param y a \code{numeric} vector.
#' @param r0 a \code{numeric} denoting the CCC under the null hypothesis. It should be in the range between
#' -1 and 1. This parameter will be ignored for the tests of Pearson, weighted Pearson, or Spearman's
#' correlation coefficient.
#' @param w \code{numeric} vector denoting the weights of the elements in vectors \code{x} and \code{y}.
#' @param method the correlation coefficient to be tested, options include Pearson's correlation coefficient
#' (\code{Pearson}), weighted Pearson correlation coefficient (\code{wtdPearson}),
#' Spearman's correlation coefficient (\code{Spearman}), Lin's concordance correlation coefficient (\code{CCC}).
#' @param alternative the alternative hypothesis, can be \code{two.sided}, \code{less}, or \code{greater}.
#' @return
#' \describe{
#' \item{\code{estimate}}{the estimated correlation coefficient.}
#' \item{\code{p.value}}{the p-value from the studentized test.}
#' \item{\code{method}}{the method for measuring correlation coefficient.}
#' \item{\code{alternative}}{the alternative hypothesis.}
#' }
#'
#' @author Mengyu Fang, Han Yu, Alan Hutson
#' @references
#'
#' Lawrence, I., & Lin, K. (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics, 255-268.
#'
#' Serfling, R. J. (2009). Approximation theorems of mathematical statistics. John Wiley & Sons.
#'
#'@importFrom stats cov cor var cov.wt na.omit complete.cases pnorm
#'
#'@examples
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnorm(250)
#' asym_test(x, y, method = "Pearson", alternative = "greater")
#'
#' asym_test(x, y, method = "Spearman", alternative = "two.sided")
#'
#' asym_test(x, y, w = rep(0.004,250), method = "wtdPearson", alternative = "less")
#'
#' asym_test(x, y, r0 = -0.5, method = "CCC", alternative = "greater")
#' @export

asym_test <- function(x, y, r0 = 0, w = NULL,
                      method =  c("Pearson", "wtdPearson", "Spearman", "CCC"),
                      alternative = c("two.sided", "less", "greater")) {
  # missing value
  if (any(is.na(x)) || any(is.na(y))) {
    warning("Missing values detected; paired data with missing values will be omitted.")
    valid_indices <- complete.cases(x, y)
    x <- x[valid_indices]
    y <- y[valid_indices]
    if (!is.null(w)) {
      w <- w[valid_indices]
    }
  }


  # "method" input
  method <- match.arg(method)
  if (!method %in% c("Pearson", "wtdPearson", "Spearman", "CCC")){
    warning("Invalid method specified. Using 'Pearson' as default." )
    method <- "Pearson"
  }

  # "alternative" input
  alternative <- match.arg(alternative)
  if (!alternative %in% c("two.sided", "less", "greater")){
    warning("Invalid alternative hypothesis specified. Using 'two.sided' as default.")
    alternative <- "two.sided"
  }

  if (method == "Pearson" || method == "Spearman"){
    if (r0 != 0){
      warning("Non-zero null hypothesis is not supported for Pearson/Spearman. Using H0: r0 = 0")
    }
    if (!is.null(w)){
      warning("Argument 'w' is unused for Pearson/Spearman/CCC methods and will be ignored.")
    }
    res <- test_rho_serfling(x, y, method, alternative)
  } else if (method =="wtdPearson") {
    if (r0 != 0){
      warning("Non-zero null hypothesis is not supported for wighted Pearson. Using default H0: r0 = 0")
    }
    if (is.null(w)){
      warning("No weight assigned. Using default Pearson correlation")
      res <- test_rho_serfling(x, y,"Pearson", alternative)
    } else {
      res <- test_rho_serfling_w(x, y, w, method, alternative)
    }
  } else if (method == "CCC") {
    if (!is.null(w)){
      warning("Argument 'w' is unused for Pearson/Spearman/CCC methods and will be ignored.")
    }
    if (r0 == 0){
      res <- test_rho_c(x, y, rc0 = 0)
    } else if (r0 >= -1 && r0 <= 1) {
      res <- test_rho_c(x, y, rc0 = r0)
    } else {
      stop("r0 must be in the range between -1 and 1.")
    }
    if (alternative != "greater") {
      warning("Only one-sided (greater) hypothesis is supported for CCC. Setting alternative to 'greater'.")
    }
  }
  return(res)
}









