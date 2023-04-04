#' @noRd
rnorm_bound <- function(n, mean, sd = NULL, lo = -Inf, up = Inf, cv = NULL, ...) {
  if(is.null(cv) & is.null(mean)) stop("At least need sd or cv")
  if(!is.null(cv)) sd <- abs(cv * mean)
  x <- rnorm(n, mean, sd, ...)
  x[x < lo] <- NA
  x[x > up] <- NA
  x <- na.omit(x)
  len_diff <- n - length(x)
  while(len_diff > 0) {
    x_new <- rnorm(len_diff, mean, sd, ...)
    x <- c(x, x_new)
    x[x < lo] <- NA
    x[x > up] <- NA
    x <- na.omit(x)
    len_diff <- n - length(x)
  }
  x
}

#' @noRd
vec_to_mat <- function(x, n, by = 1) {

  if(by == 1) { # for row
    result <- t(t(rep(1, n))) %*% x
    colnames(result) <- names(x)
  } else if(by == 2) {
    result <- t(t(x)) %*% rep(1, n)
    rownames(result) <- names(x)
  }

  return(result)

}


#' @noRd
def_frac <- function(name, value) {
  v <- numeric(length = length(phytodive_iop$name_phyto))
  names(v) <- names(phytodive_iop$name_phyto)
  v[name] <- value
  v
}

#### IOP to Rrs ####

#' @noRd
IOP2AOP_Albert_Mobley_03 <- Albert_Mobley_03 <- function(a, bb, theta_s = 30, theta_v = 0, windspd = 7) {

  x = bb / (bb + a)

  p1 = 0.0512
  p2 = 4.6659
  p3 = -7.8387
  p4 = 5.4571
  p5 = 0.1098
  p6 = -0.0044
  p7 = 0.4021

  rrs  <-  p1 * x * (1 + p2 * x + p3 * x^2 + p4 * x^3) *
    (1 + p5 / cos(theta_s/180*pi)) *
    (1 + p7 / cos(theta_v/180*pi)) *
    (1 + p6 * windspd)

  Rrs <- 0.526 * rrs / (1 - 2.164 * rrs)

  Rrs

}

#' @noRd
IOP2AOP_Lee_11 <- function(bbw,
                           bbp,
                           at,
                           bbt,
                           gw0 = 0.05881474,
                           gw1 = 0.05062697,
                           gp0 = 0.03997009,
                           gp1 = 0.1398902) {

  # g coefficients are for SZA == 30, VAA == 0, and VZA == 0

  k = at + bbt

  (gw0 + gw1 * bbw / k) * (bbw/k) + (gp0 + gp1 * bbp/k) * (bbp/k)

}



