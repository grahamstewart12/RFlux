
qc_stat <- function(path_rawdata, ext_tstamp = c("START", "END"), 
                    path_output = NULL, file_name = NULL) {

  ##############################################################################
  # 	Raw Data Processing
  ##############################################################################

  # Read data if passed as a file connection
  if (inherits(path_rawdata, "data.frame")) {
    raw_data <- as.data.frame(path_rawdata)
    tstamp <- as.character(lubridate::ymd_hms(stringr::str_sub(
      paste(raw_data$DATE[1], raw_data$TIME[1]), 1, -5
    ), tz = "GMT"))
  } else {
    if (ext_tstamp == "START") {
      tstamp <- as.character(lubridate::ymd_hms(
        stringr::str_sub(path_rawdata, -30, -14), tz = "GMT"
      ))
    } 
    if (ext_tstamp == "END") {
      tstamp <- as.character(lubridate::ymd_hms(
        stringr::str_sub(path_rawdata, -30, -14), tz = "GMT"
      ) - 1800)
    }
    raw_data <- data.table::fread(
      path_rawdata, sep = ",", header = TRUE, data.table = FALSE, 
      na.strings = c("-9999")
    )
  }
  
  n <- nrow(raw_data)
  HZ <- if (n > 19000) 20 else 10
  N <- max(HZ * 60 * 30, n)
  delta <- if (n < N) N - n else 0
  
  fluxes <- c(H = "T_SONIC", FC = "CO2", LE = "H2O", Fch4 = "CH4")
  nf <- length(fluxes)
  fmr <- vctrs::vec_init(numeric(), nf)
  lgd <- vctrs::vec_init(numeric(), nf)
  for (i in 1:nf) {
    if (length(which(is.na(raw_data$W + raw_data[[fluxes[i]]]))) == n) {
      fmr[i] <- 100
      lgd[i] <- 1800
    } else {
      stat_temp <- imputeTS::statsNA(
        raw_data$W + raw_data[[fluxes[i]]], printOnly = FALSE
      )
      fmr[i] <- (stat_temp$numberNAs + delta) / N * 100
      lgd[i] <- max(delta, stat_temp$naGapLongest, na.rm = TRUE) / HZ
    }
  }
  
  
  L <- 25 # lag max for LSR test
  ccf2 <- purrr::partial(ccf, lag.max = L, plot = FALSE, na.action = na.pass)
  
  if (is.null(raw_data$V) | length(which(is.na(raw_data$V))) > N * 0.95) {
    ipt_v <- rep(NA, 8)
  } else {
    ipt_v <- inst_prob_test(raw_data$V)
  }
  
  if (fmr[1] > 15 | lgd[1] > 180) {
    ipt_w <- rep(NA, 8)
    SADiag <- NA
  } else {
    ipt_w <- inst_prob_test(raw_data$W)
    SADiag <- length(which(raw_data$SA_DIAG == 0))
  }
  
  D0 <- vctrs::vec_init(NA, nf)
  lrt <- vctrs::vec_init(NA, nf)
  ipt <- purrr::map(vctrs::vec_init(list(), nf), ~ vctrs::vec_init(NA, 8))
  ipt <- vctrs::vec_init(list(rep(NA, 8)), nf)
  M98 <- vctrs::vec_init(NA, nf)
  cov <- vctrs::vec_init(NA, nf)
  for (i in 1:nf) {
    if (fmr[i] <= 15 | lgd[i] <= 180) {
      ind_w <- which(diff(raw_data$W) == 0) + 1
      ind_f <- which(diff(raw_data[[fluxes[i]]]) == 0) + 1
      D0[i] <- max(length(ind_w), length(ind_f))
      cor_ori <- ccf2(raw_data$W, raw_data[[fluxes[i]]])
      if (D0[i] < N * 0.9 & D0[i] > 0) {
        cor_sub <- ccf2(
          replace(raw_data$W, ind_w, NA), 
          replace(raw_data[[fluxes[i]]], ind_f, NA)
        )
        lrt[i] <- summary(lm(cor_ori$acf ~ cor_sub$acf - 1))$r.squared
      } else if (D0[i] < N * 0.9) {
        lrt[i] <- 1
      } else {
        lrt[i] <- -1
      }
      ipt[[i]] <- inst_prob_test(raw_data[[fluxes[i]]])
      M98[i] <- mahrt(cbind(raw_data$W, raw_data[[fluxes[i]]]))
      cov[i] <- cov(raw_data$W, raw_data[[fluxes[i]]], use = "complete.obs")
    }
  }
  
  
  
  results <- unlist(c(
    tstamp, SADiag, 
    fmr[1], lgd[1], fmr[2], lgd[2], fmr[3], lgd[3], fmr[4], lgd[4],
    ipt_v, ipt_w,
    ipt[1], cov[1], D0[1], lrt[1], M98[1],
    ipt[2], cov[2], D0[2], lrt[2], M98[2],
    ipt[3], cov[3], D0[3], lrt[3], M98[3],
    ipt[4], cov[4], D0[4], lrt[4], M98[4]
  ), use.names = FALSE)
  
  ipt_names <- c("Skew", "Kurt", "KID0", "KID1", "HF5", "HF10", "HD5", "HD10")

  names(results) <- c(
    "TSTAMP", "SADiag", 
    "FMR_H", "LGD_H", "FMR_Fc", "LGD_Fc", "FMR_LE", "LGD_LE", "FMR_Fch4", "LGD_Fch4",
    paste0(ipt_names, "_v"),
    paste0(ipt_names, "_w"),
    paste0(ipt_names, "_ts"), "COV_wts", "N0_H", "LSR_H", "M98_H",
    paste0(ipt_names, "_co2"), "COV_wco2", "N0_Fc", "LSR_Fc", "M98_Fc",
    paste0(ipt_names, "_h2o"), "COV_wh2o", "N0_LE", "LSR_LE", "M98_LE",
    paste0(ipt_names, "_ch4"), "COV_wch4", "N0_Fch4", "LSR_Fch4", "M98_Fch4"
  )

  if (!is.null(path_output) & !is.null(file_name)) {
    write.table(
      t(results), paste0(path_output, "/", file_name, ".csv"), 
      quote = FALSE, sep = ",", row.names = FALSE
    )
  } 

  results
}
