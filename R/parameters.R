# Constant for studentization
tau <- function(x, y) {
  xbar <- mean(x)
  ybar <- mean(y)
  mu_20_vec <- (x - xbar)^2
  mu_02_vec <- (y - ybar)^2
  mu_22_vec <- mu_20_vec * mu_02_vec
  mu_20 <- mean(mu_20_vec)
  mu_02 <- mean(mu_02_vec)
  mu_22 <- mean(mu_22_vec)
  tau_hat <- sqrt(mu_22/(mu_20*mu_02))
  return(tau_hat)
}

# Studentized Pearson correlation coefficient
rho_stu_func <- function(x,y) {
  rho_ori <- cor(x, y, method = "pearson")
  tau_hat <- tau(x, y)
  rho_stu <- rho_ori/tau_hat
  return(rho_stu)
}

# Studentized Spearman correlation coefficient
rho_stu_func_s <- function(x, y) {
  return(rho_stu_func(rank(x), rank(y)))
}

# CCC
stat_func <- function(x, y, rc0) {
  C_hat <- 2*sqrt(var(x)*var(y))/(var(x) + var(y) + (mean(x) - mean(y))^2)
  rho_hat_0 <- rc0/C_hat
  u <- (x - mean(x)) / sqrt(var(x))
  v <- (y - mean(y)) / sqrt(var(y))
  v2 <- (v/rho_hat_0-u)/sqrt(1 + 1/rho_hat_0^2)
  rho_hat <- cor(u, v2)
}


rho_ccc_func <- function(x, y) {
  xbar=mean(x)
  ybar=mean(y)
  Sigma <- cov(data.frame(x, y), method = "pearson")
  rho_c <- 2*Sigma[1,2]/(Sigma[1,1] + Sigma[2,2] + (xbar - ybar)^2)
  return(rho_c)
}

# studentized CCC
rho_ccc_s_func <- function(x, y) {
  rho_c <- rho_ccc_func(x, y)
  tau_hat <- tau(x, y)
  rho_s <- rho_c/tau_hat
  return(rho_s)
}

# Unweighted asymptotic normal approximation #
sigma_serfling <- function(x, y) {
  n <- length(x)

  x_bar <- mean(x)
  y_bar <- mean(y)
  x2 <- x^2
  y2 <- y^2
  xy <- x*y
  x2_bar <- mean(x2)
  y2_bar <- mean(y2)
  xy_bar <- mean(xy)

  df <- data.frame(x, y, x2, y2, xy)
  w <- cov(df)

  sigma_x <- sqrt(var(x))
  sigma_y <- sqrt(var(y))
  rho <- cor(x, y)

  d <- c()

  d[1] <- (rho*x_bar)/sigma_x^2 - y_bar/(sigma_x*sigma_y)
  d[2] <- (rho*y_bar)/sigma_y^2 - x_bar/(sigma_x*sigma_y)
  d[3] <- -rho/(2*sigma_x^2)
  d[4] <- -rho/(2*sigma_y^2)
  d[5] <- 1/(sigma_x*sigma_y)

  se <- sqrt(t(d) %*% w %*% d/n)
  return(se)
}

