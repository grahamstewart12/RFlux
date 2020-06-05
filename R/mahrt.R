
mahrt <- function(x) {
  x <- na.omit(x, na.rm = TRUE)
  ifelse(nrow(x) > 20000, wl <- 5995, wl <- 2995)
  cov <- cov(x[, 1], x[, 2], use = "complete.obs")
  covs <- zoo::rollapply(
    x, FUN = function(x) cov(x[, 1], x[, 2], use = "complete.obs"), 
    width = wl, by = wl, by.column = FALSE
  )
  cov_w <- zoo::rollapply(
    x, FUN = function(x) cov(x[, 1], x[, 2], use = "complete.obs"), 
    width = wl / 6, by = wl / 6, by.column = FALSE
  )
  sigma_b <- sqrt(sum((covs - cov)^2) / (length(covs) - 1))
  cov_w_len <- length(cov_w)
  ifelse(
    cov_w_len >= 6, 
    sigma_w1 <- sqrt(1 / 5 * sum((cov_w[1:6] - covs[1])^2)), 
    sigma_w1 <- NA
  )
  ifelse(
    cov_w_len >= 12, 
    sigma_w2 <- sqrt(1 / 5 * sum((cov_w[7:12] - covs[2])^2)), 
    sigma_w2 <- NA
  )
  ifelse(
    cov_w_len >= 18, 
    sigma_w3 <- sqrt(1 / 5 * sum((cov_w[13:18] - covs[3])^2)), 
    sigma_w3 <- NA
  )
  ifelse(
    cov_w_len >= 24, 
    sigma_w4 <- sqrt(1 / 5 * sum((cov_w[19:24] - covs[4])^2)), 
    sigma_w4 <- NA
  )
  ifelse(
    cov_w_len >= 30, 
    sigma_w5 <- sqrt(1 / 5 * sum((cov_w[25:30] - covs[5])^2)), 
    sigma_w5 <- NA
  )
  ifelse(
    cov_w_len > 30 & length(cov_w) <= 36, 
    sigma_w6 <- sqrt(1 / cov_w_len * sum((cov_w[31:cov_w_len] - covs[6])^2)), 
    sigma_w6 <- NA
  )
  sigma_wi <- c(sigma_w1, sigma_w2, sigma_w3, sigma_w4, sigma_w5, sigma_w6)
  sigma_w <- mean(sigma_wi, na.rm = TRUE)
  stat <- sigma_b / (sigma_w / sqrt(length(na.omit(sigma_wi))))
  
  list("M98" = stat)
}
