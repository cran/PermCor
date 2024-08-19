# permutation test for pearson or spearman correlation coefficients
perm_cor <- function(x, y, B = 500, method = c("Pearson", "Spearman"), alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  if (method == "Spearman") {
    x <- rank(x)
    y <- rank(y)
  }
  rho_c <- cor(x, y, method = "pearson")
  rho_s <- rho_stu_func(x, y) # studentized correlation
  rho_s_star <- replicate(B, rho_stu_func(sample(x), y), simplify = TRUE)
  if (alternative == "two.sided") {
    p_value <- mean(rho_s_star > abs(rho_s) | rho_s_star < -abs(rho_s))
  } else if (alternative == "less") {
    p_value <- mean(rho_s_star < rho_s)
  } else if (alternative == "greater") {
    p_value <- mean(rho_s_star > rho_s)
  }
  res <- list(estimate = rho_c, p.value = p_value, method = method, alternative = alternative)
  return(res)
}


# permutation test for weighted pearson correlation coefficients
perm_cor_w <- function(x, y, w, B = 500, method, alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  df <- data.frame(x, y)
  se <- sigma_serfling(x, y)*sqrt(sum(w^2)*n)
  wcor0 <- cov.wt(df, wt = w, cor = TRUE)$cor[1,2]
  wcor <- (wcor0/se)[1,1]
  wcor_star <- c()
  for (i in 1:B) {
    perm_1 <- sample(n)
    perm_2 <- sample(n)
    df_star <- df
    df_star[,1] <- df[perm_1,1]
    df_star[,2] <- df[perm_2,2]
    wcor_star[i] <- cov.wt(df_star, wt = w, cor = TRUE)$cor[1,2]
    se_star <- sigma_serfling(df_star[,1], df_star[,2])*sqrt(sum(w^2)*n)
    wcor_star[i] <- wcor_star[i]/se_star
  }
  if (alternative == "two.sided"){
    p_value <- mean(wcor_star > abs(wcor) | wcor_star < -abs(wcor))
  } else if (alternative == "less") {
    p_value <- mean(wcor_star < wcor)
  } else if (alternative == "greater"){
    p_value <- mean(wcor_star > wcor)
  }
  res <- list(estimate = wcor0, p.value = p_value, method = "wtdPearson", alternative = alternative)
  return(res)
}

# permutation test for ccc
perm_ccc <- function(x, y, B = 1000) {
  n <- length(x)
  rho_c <- rho_ccc_func(x, y)
  rho_s <- rho_ccc_s_func(x, y)
  rho_s_star <- replicate(B, rho_ccc_s_func(sample(x), y), simplify = TRUE)
  p_value <- mean(rho_s_star > rho_s)
  res <- list(estimate = rho_c, p.value = p_value, method = "CCC", alternative = "greater")
  return(res)
}


perm_ccc_2 <- function(x, y, rc0, B = 1000) {

  n <- length(x)
  u <- (x - mean(x)) / sqrt(var(x))
  v <- (y - mean(y)) / sqrt(var(y))
  stat_sample <- stat_func(x, y, rc0)

  ###########
  rho_sample <- cor(u, v)
  C_0 <- 2*sqrt(var(x)*var(y))/(var(x) + var(y) + (mean(x) - mean(y))^2)
  rho_hat_0 <- rc0/C_0
  if (rho_hat_0 >= 1) return(res <- list(estimate = rho_ccc_func(x, y), p.value = 1, method = "CCC", alternative = "greater"))

  # de-correlate
  u_1 <- u
  v_1 <- (u - v/rho_sample)/sqrt(1 + 1/rho_sample^2)
  # re-correlate
  u_2 <- u_1
  v_2 <- rho_hat_0*u_1 + sqrt(1 - rho_hat_0^2)*v_1
  # re-scale and shift
  x_2 <- u_2*sqrt(var(x)) + mean(x)
  y_2 <- v_2*sqrt(var(y)) + mean(y)
  #######

  rho_star_vec <- c()
  for (b in 1:n) {
    ind_x <- setdiff(1:n, b)
    x_star <- x_2[ind_x]
    y_star <- y_2[ind_x]
    rho_star_vec[b] <- stat_func(x_star, y_star, rc0)
  }

  tau_hat <- sqrt(var(rho_star_vec))
  stat_star_vec <- c()

  for (b in 1:B) {
    ind_x <- sample(1:n, n, replace = TRUE)
    ind_y <- sample(1:n, n, replace = TRUE)
    u_star <- u[ind_x]
    v_star <- v[ind_y]
    stat_star_vec[b] <- cor(u_star, v_star)
  }

  p_value <- mean(stat_star_vec*tau_hat*sqrt(n*(n-1)) > stat_sample)
  res <- list(estimate = rho_ccc_func(x, y), p.value = p_value, method = "CCC", alternative = "greater")
  return(res)
}


# Pearson
test_rho_serfling <- function(x, y, method = c("Pearson", "Spearman"), alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  if (method == "Spearman") {
    x <- rank(x)
    y <- rank(y)
  }
  r <- cor(x, y)
  se <- sigma_serfling(x, y)
  if (alternative == "two.sided"){
    p_value <- 2 * (1 - pnorm(abs(r), 0, se))
  } else if (alternative == "less"){
    p_value <- pnorm(r, 0, se)
  } else if (alternative == "greater"){
    p_value <- 1 - pnorm(r, 0, se)
  }
  res <- list(estimate = r, p.value = p_value, method = method, alternative = alternative)
  return(res)
}

# Weighted Pearson
test_rho_serfling_w <- function(x, y, w, method, alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  rho <- cov.wt(data.frame(x, y), wt = w, cor = TRUE)$cor[1,2]
  se <- sigma_serfling(x, y)*sqrt(sum(w^2)*n)
  if (alternative == "two.sided"){
    p_value <- 2 * (1 - pnorm(abs(rho), 0, se))
  } else if (alternative == "less"){
    p_value <- pnorm(rho, 0, se)
  } else if (alternative == "greater"){
    p_value <- 1 - pnorm(rho, 0, se)
  }
  res <- list(estimate = rho, p.value = p_value, method = method, alternative = alternative)
  return(res)
}

# CCC
rho_c_func <- function(x, y) {
  x2=x^2;
  y2=y^2;
  xy=x*y;

  xbar=mean(x)
  ybar=mean(y)

  W <- matrix(NA, nrow = 2, ncol = 2)
  fct <- 1

  W[1,1]=var(x) * fct;
  W[1,2]=cov(x, y) * fct;
  W[2,1]=W[1,2];
  W[2,2]=var(y) * fct;

  rho_c <- 2*W[1,2]/(W[1,1] + W[2,2] + (xbar - ybar)^2)
  return(rho_c)
}



sigma_c_func <- function(x, y) {

  n <- length(x)
  x2=x^2;
  y2=y^2;
  xy=x*y;

  xbar=mean(x)
  ybar=mean(y)
  x2bar=mean(x2)
  y2bar=mean(y2)
  xybar=mean(xy)

  W <- matrix(NA, nrow = 5, ncol = 5)

  fct <- 1

  W[1,1]=var(x) * fct;
  W[1,2]=cov(x, y) * fct;
  W[1,3]=cov(x, x2) * fct;
  W[1,4]=cov(x, y2) * fct;
  W[1,5]=cov(x, xy) * fct;

  W[2,1]=W[1,2];
  W[2,2]=var(y) * fct;
  W[2,3]=cov(y, x2) * fct;
  W[2,4]=cov(y, y2) * fct;
  W[2,5]=cov(y, xy) * fct;

  W[3,1]=W[1,3];
  W[3,2]=W[2,3];
  W[3,3]=var(x2) * fct;
  W[3,4]=cov(x2, y2) * fct;
  W[3,5]=cov(x2, xy) * fct;

  W[4,1]=W[1,4];
  W[4,2]=W[2,4];
  W[4,3]=W[3,4];
  W[4,4]=var(y2) * fct;
  W[4,5]=cov(y2, xy) * fct;

  W[5,1]=W[1,5];
  W[5,2]=W[2,5];
  W[5,3]=W[3,5];
  W[5,4]=W[4,5];
  W[5,5]=var(xy) * fct;

  d <- c()
  denom <- x2bar + y2bar -2*xbar*ybar

  d[1] = 2*ybar*(2*xybar - x2bar - y2bar)/denom^2
  d[2] = 2*xbar*(2*xybar - x2bar - y2bar)/denom^2
  d[3] = -2*(xybar - xbar*ybar)/denom^2
  d[4] = d[3]
  d[5] = 2/denom

  sigma_c <- sqrt(t(d) %*% W %*% d/n)
  return(sigma_c)
}



test_rho_c <- function(x, y, rc0) {
  rho_c <- rho_c_func(x, y)
  sd_c <- sigma_c_func(x, y)
  p_value <- 1 - pnorm(rho_c-rc0, 0, sd_c)
  res <- list(estimate = rho_c_func(x, y), p.value = p_value, method = "CCC", alternative = "greater")
  return(res)
}


