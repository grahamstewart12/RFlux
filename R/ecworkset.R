
ecworkset <- function(path_ep, path_qc, path_md, path_stat, 
                      path_output = NULL, FileName = NULL) {
  ep <- data.table::fread(
    path_ep, header = FALSE, skip = 3, na.strings = c("-9999", "-9999.0"), 
    data.table = FALSE
  )
  colnames(ep) <- as.vector(
    scan(path_ep, nlines = 1, skip = 1, sep = ",", what = "character")
  )
  timestamp_ep <- as.POSIXct(
    paste(ep[, 2], ep[, 3]), format = "%Y-%m-%d %H:%M", tz = "GMT"
  ) - 1800

  qc <- data.table::fread(
    path_qc, header = FALSE, skip = 3, na.strings = c("-9999", "-9999.0"), 
    data.table = FALSE
  )
  colnames(qc) <- as.vector(
    scan(path_qc, nlines = 1, skip = 1, sep = ",", what = "character")
  )
  timestamp_qc <- as.POSIXct(
    paste(qc[, 2], qc[, 3]), format = "%Y-%m-%d %H:%M", tz = "GMT"
  ) - 1800

  md <- data.table::fread(
    path_md, header = FALSE, skip = 1, na.strings = c("-9999", "-9999.0"), 
    data.table = FALSE
  )
  colnames(md) <- as.vector(
    scan(path_md, nlines = 1, skip = 0, sep = ",", what = "character")
  )
  timestamp_md <- as.POSIXct(
    paste(md[, 2], md[, 3]), format = "%Y-%m-%d %H:%M", tz = "GMT"
  ) - 1800

  stat <- data.table::fread(
    path_stat, header = TRUE, data.table = FALSE, na.strings = c(NA, "-9999")
  )
  timestamp_stat <- as.POSIXct(
    as.character(stat[, 1]), format = "%Y%m%d%H%M", tz = "GMT"
  )
  stat.xts <- xts(stat[, -1], order.by = timestamp_stat)

  VarSelOut <- c(
    which(colnames(ep) == "daytime"),
    which(colnames(ep) == "file_records"),
    which(colnames(ep) == "used_records"),
    which(colnames(ep) == "(z-d)/L"),
    which(colnames(ep) == "H"),
    which(colnames(ep) == "qc_H"),
    which(colnames(ep) == "rand_err_H"),
    which(colnames(ep) == "LE"),
    which(colnames(ep) == "qc_LE"),
    which(colnames(ep) == "rand_err_LE"),
    which(colnames(ep) == "co2_flux"),
    which(colnames(ep) == "qc_co2_flux"),
    which(colnames(ep) == "rand_err_co2_flux"),
    which(colnames(ep) == "co2_strg"),
    which(colnames(ep) == "u*"),
    which(colnames(ep) == "u_rot"),
    which(colnames(ep) == "w_rot"),
    which(colnames(ep) == "wind_dir"),
    which(colnames(ep) == "wind_speed"),
    which(colnames(ep) == "u_var"),
    which(colnames(ep) == "v_var"),
    which(colnames(ep) == "w_var"),
    which(colnames(ep) == "L"),
    which(colnames(ep) == "sonic_temperature"),
    which(colnames(ep) == "air_temperature"),
    which(colnames(ep) == "air_pressure"),
    which(colnames(ep) == "air_density"),
    which(colnames(ep) == "air_heat_capacity"),
    which(colnames(ep) == "ts_var"),
    which(colnames(ep) == "co2_mole_fraction"),
    which(colnames(ep) == "co2_var"),
    which(colnames(ep) == "h2o_mole_fraction"),
    which(colnames(ep) == "h2o_var"),
    which(colnames(ep) == "co2_scf"),
    which(colnames(ep) == "h2o_scf"),
    which(colnames(ep) == "H_scf")
  )


  VarSelQC <- c(
    which(colnames(qc) == "dev(w)")[1], ## Instationary Test
    which(colnames(qc) == "dev(ts)")[1],
    which(colnames(qc) == "dev(co2)"),
    which(colnames(qc) == "dev(h2o)"),
    which(colnames(qc) == "dev(w/ts)"),
    which(colnames(qc) == "dev(w/co2)"),
    which(colnames(qc) == "dev(w/h2o)"),
    which(colnames(qc) == "dev(w)")[2]
  ) ## ITC

  VarSelMD <- c(
    which(colnames(md) == "acquisition_frequency"),
    which(colnames(md) == "canopy_height"),
    which(colnames(md) == "master_sonic_height"),
    which(colnames(md) == "master_sonic_north_offset")
  )

  ##############################################################################
  #
  ### IRGA Diagnostic Management
  #
  ##############################################################################
  IRGA <- substr(md[1, "co2_irga_model"], 1, 6)

  if (IRGA == "li7200") {
    GA_TCellDiag <- apply(
      cbind(ep[, "t_out_LI-7200"], ep[, "t_in_LI-7200"]), 1, 
      function(x) min(x, na.rm = TRUE)
    )
    GA_DiagVar <- cbind(
      ep[, "head_detect_LI-7200"], GA_TCellDiag, ep[, "aux_in_LI-7200"], 
      ep[, "delta_p_LI-7200"], ep[, "chopper_LI-7200"], 
      ep[, "detector_LI-7200"], ep[, "pll_LI-7200"], ep[, "sync_LI-7200"]
    )
    GA_Diag <- apply(GA_DiagVar, MARGIN = 1, function(x) sum(x, na.rm = TRUE))
  }

  if (IRGA == "li7500") {
    GA_DiagVar <- cbind(
      ep[, "chopper_LI-7500"], ep[, "detector_LI-7500"], ep[, "pll_LI-7500"], 
      ep[, "sync_LI-7500"]
    )
    GA_Diag <- apply(GA_DiagVar, MARGIN = 1, function(x) sum(x, na.rm = TRUE))
  }

  GA_Diag.xts <- xts(GA_Diag, order.by = timestamp_ep)
  if (IRGA != "li7200" & IRGA != "li7500") {
    GA_Diag.xts <- xts(rep(0, length(timestamp_ep)), order.by = timestamp_ep)
  } 

  ##############################################################################

  ep.xts <- xts(ep[, VarSelOut], order.by = timestamp_ep)
  qc.xts <- xts(qc[, VarSelQC], order.by = timestamp_qc)
  md.xts <- xts(md[, VarSelMD], order.by = timestamp_md)
  stat.xts <- xts(stat[, -1], order.by = timestamp_stat)


  ecworkset.xts <- merge(ep.xts, GA_Diag.xts, qc.xts, md.xts, stat.xts)
  set2exp <- data.frame(
    format(time(ecworkset.xts), "%Y%m%d%H%M", tz = "GMT"), 
    zoo::coredata(ecworkset.xts)
  )
  colnames(set2exp) <- c(
    "TIMESTAMP", "DNtime", "file_records", "used_records", "Stability", "H", 
    "qcH", "ruH", "LE", "qcLE", "ruLE", "CO2flux", "qcCO2flux", "ruCO2", 
    "CO2str", "ustar", "U", "W", "WDir", "WSpeed", "u_var", "v_var", "w_var", 
    "MOL", "TSonic", "AT", "AP", "rho", "cp", "ts_var", "CO2mf", "co2_var", 
    "H2Omf", "h2o_var", "CO2scf", "H2Oscf", "Hscf", "GADiag", "dw", "dts", 
    "dco2", "dh2o", "dwts", "dwco2", "dwh2o", "itc_w", "acquisition_frequency", 
    "canopy_height", "SA_height", "SA_NorthOffset", colnames(stat)[-1]
  )


  if (!is.null(path_output) & !is.null(FileName)) {
    data.table::fwrite(
      set2exp, paste0(path_output, "/", FileName, ".csv"), sep = ",", 
      row.names = FALSE, col.names = TRUE, quote = FALSE, na = "-9999"
    )
  } 

  set2exp
}
