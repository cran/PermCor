#' Robust Permutation Tests of Correlation Coefficients
#'
#' This function performs robust permutation tests for various correlation coefficients, providing
#' reliable type I error control under general scenarios, especially when the sample size is small and
#' two variables are dependent but uncorrelated.
#'
#' @details
#' #' The test supports the following correlation coefficients: Pearson correlation coefficient, Weighted
#' Pearson correlation coefficient, Spearman correlation coefficient, and Lin's concordance correlation coefficient (CCC)
#'
#' For Pearson, weighted Pearson, and Spearman correlation coefficients, the test supports a zero null
#' hypothesis. The alternative hypothesis can be either one-sided or two-sided.
#'
#' For Lin's concordance correlation coefficient (CCC), the test accommodates a more general null
#' hypothesis. Currently, the test only supports a one-sided alternative hypothesis (greater).
#'
#' @param x a \code{numeric} vector.
#' @param y a \code{numeric} vector.
#' @param B an \code{integer} number of permutations.
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
#' Lawrence, I., & Lin, K. (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics, 255-268.
#'
#' DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for correlation and regression coefficients. Journal of the American Statistical Association, 112(519), 1211-1220.
#'
#' Hutson, A. D., & Yu, H. (2021). A robust permutation test for the concordance correlation coefficient. Pharmaceutical Statistics, 20(4), 696-709.
#'
#' Yu, H., & Hutson, A. D. (2024). A robust Spearman correlation coefficient permutation test. Communications in Statistics-Theory and Methods, 53(6), 2141-2153.
#'
#' Yu, H., & Hutson, A. D. (2024). Inferential procedures based on the weighted Pearson correlation coefficient test statistic. Journal of Applied Statistics, 51(3), 481-496.
#'
#'
#' @importFrom stats cov cor var cov.wt na.omit complete.cases
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(20)
#' y <- rnorm(20)
#' perm_test(x, y, B = 500, method = "Pearson", alternative = "greater")
#'
#' perm_test(x, y, B = 500, method = "Spearman", alternative = "two.sided")
#'
#' perm_test(x, y, B = 500, w = rep(0.05,20), method = "wtdPearson", alternative = "less")
#'
#' perm_test(x, y, B = 500, r0 = -0.5, method = "CCC", alternative = "greater")
#'
#' @export

perm_test <- function(x, y, B = 1000, r0 = 0, w = NULL,
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
    res <- perm_cor(x, y, B, method, alternative)
  } else if (method =="wtdPearson") {
    if (r0 != 0){
      warning("Non-zero null hypothesis is not supported for wighted Pearson. Using default H0: r0 = 0")
    }
    if (is.null(w)){
      warning("No weight assigned. Using default Pearson correlation")
      res <- perm_cor(x, y, B, "Pearson", alternative)
    } else {
      res <- perm_cor_w(x, y, w, B, method, alternative)
    }
  } else if (method == "CCC") {
    if (!is.null(w)){
      warning("Argument 'w' is unused for Pearson/Spearman/CCC methods and will be ignored.")
    }
    if (r0 == 0){
      res <- perm_ccc(x, y, B)
    } else if (r0 >= -1 && r0 <= 1) {
      res <- perm_ccc_2(x, y, rc0 = r0, B)
    } else {
      stop("r0 must be in the range between -1 and 1.")
    }
    if (alternative != "greater") {
      warning("Only one-sided (greater) hypothesis is supported for CCC. Setting alternative to 'greater'.")
    }
  }
  return(res)
}








