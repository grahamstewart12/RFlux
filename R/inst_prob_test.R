
inst_prob_test <- function(x) {
  
  # This script seems to be set up for mole fractions, not densities
  # (i.e. closed-path systems, not open-path)
  # Also doesn't consider CH4, which has a much smaller variance
  
  # CHANGE: sigma_min from 1e-2 to 1e-4 to better represent densities & ch4
  sigma_min <- 0.0001
  
  # CHANGE: diff_min from 1e-3 to 1e-5 to better represent densities & ch4
  diff_min <- 0.00001
  
  flucts <- x - mean(x, na.rm = TRUE)
  sigma_f <- max(sigma_min, robustbase::Qn(na.omit(flucts)))
  n_spike1 <- length(which(abs(flucts) > 5 * sigma_f))
  n_spike2 <- length(which(abs(flucts) > 10 * sigma_f))

  x_n <- na.omit(x)
  d0 <- diff(x)
  ind <- which(abs(d0) < diff_min)
  x_r <- replace(x, ind + 1, NA)
  d1 <- na.omit(diff(x_r))

  K_VM97 <- as.numeric(
    3 + timeDate::kurtosis(egcm::detrend(x_n), na.rm = TRUE)
  )
  S_VM97 <- as.numeric(
    timeDate::skewness(egcm::detrend(x_n), na.rm = TRUE)
  )
  KID0 <- as.numeric(3 + timeDate::kurtosis(d0, na.rm = TRUE))
  KID1 <- as.numeric(3 + timeDate::kurtosis(diff(na.omit(x_r)), na.rm = TRUE))
  
  if (length(d1) > 1000) {
    sigma_d <- max(sigma_min, robustbase::Qn(d1))
    n_spike3 <- length(which(abs(d0) > 5 * sigma_d))
    n_spike4 <- length(which(abs(d0) > 10 * sigma_d))
  } else {
    n_spike3 <- NA
    n_spike4 <- NA
  }

  list(
    "Skew" = S_VM97, "Kurt" = K_VM97, "KID0" = KID0, "KID1" = KID1, 
    "HF4" = n_spike1, "HF1" = n_spike2, "HD4" = n_spike3, "HD1" = n_spike4
  )
}
