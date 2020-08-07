
mahrt <- function(x) {
  
  x <- na.omit(x, na.rm = TRUE)
  
  wl <- if (nrow(x) > 20000) 5995 else 2995
  
  cov <- cov(x[, 1], x[, 2], use = "complete.obs")
  
  covs <- wapply(
    x, FUN = function(x) cov(x[, 1], x[, 2], use = "complete.obs"), 
    width = wl, by = wl
  )
  cov_w <- wapply(
    x, FUN = function(x) cov(x[, 1], x[, 2], use = "complete.obs"), 
    width = wl / 6, by = wl / 6
  )
  
  sigma_b <- sqrt(sum((covs - cov)^2) / (length(covs) - 1))
  cov_w_len <- length(cov_w)
  
  sigma_wi <- rep(NA, 6)
  
  if (cov_w_len >= 6) {
    sigma_wi[1] <- sqrt(1/5 * sum((cov_w[1:6] - covs[1])^2))
  }
  
  if (cov_w_len >= 12) {
    sigma_wi[2] <- sqrt(1/5 * sum((cov_w[7:12] - covs[2])^2))
  }
  
  if (cov_w_len >= 18) {
    sigma_wi[3] <- sqrt(1/5 * sum((cov_w[13:18] - covs[3])^2))
  }
  
  if (cov_w_len >= 24) {
    sigma_wi[4] <- sqrt(1/5 * sum((cov_w[19:24] - covs[4])^2))
  }
  
  if (cov_w_len >= 30) {
    sigma_wi[5] <- sqrt(1/5 * sum((cov_w[25:30] - covs[5])^2))
  }
  
  if (cov_w_len > 30 & length(cov_w) <= 36) {
    sigma_wi[6] <- sqrt(1 / cov_w_len * sum((cov_w[31:cov_w_len] - covs[6])^2))
  }
  
  sigma_w <- mean(sigma_wi, na.rm = TRUE)
  stat <- sigma_b / (sigma_w / sqrt(length(na.omit(sigma_wi))))
  
  stat
}

wapply <- function(x, width, by = NULL, FUN = NULL, ...) {
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width
  
  len <- nrow(x)
  seq1 <- seq(1, len - width + 1, by = by)
  seq2 <- lapply(seq1, function(x) x:(x + width - 1))
  
  out <- lapply(seq2, function(a) FUN(x[a, ], ...))
  base:::simplify2array(out, higher = TRUE)
}

