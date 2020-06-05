
eddypro_run <- function(site_id, path_eddypro_bin, path_eddypro_projfiles, 
                        show_log = TRUE) {
  workdir <- getwd()
  setwd(path_eddypro_projfiles)
  setwd("..")
  eddypro_path <- getwd()
  dir.create(paste0(getwd(), "/tmp"))

  if (Sys.info()[["sysname"]] == "Darwin") op_sys <- " -s mac "
  if (Sys.info()[["sysname"]] == "Windows") op_sys <- " -s win "
  if (Sys.info()[["sysname"]] == "Linux") op_sys <- " -s linux "

  file.copy(
    from = path_eddypro_bin, to = eddypro_path, 
    recursive = TRUE, overwrite = TRUE
  )

  setwd(paste0(eddypro_path, "/bin"))
  getwd()

  rp_command <- paste0(
    "./eddypro_rp", op_sys, path_eddypro_projfiles, "/", site_id, ".eddypro"
  )
  rp <- system(rp_command, intern = !show_log)

  fcc_command <- paste0(
    "./eddypro_fcc", op_sys, path_eddypro_projfiles, "/", site_id, ".eddypro"
  )
  fcc <- system(fcc_command, intern = !show_log)

  unlink(paste0(eddypro_path, "/bin"), recursive = TRUE)
  unlink(paste0(eddypro_path, "/tmp"), recursive = TRUE)
  setwd(workdir)
  cat(
    "\n The results of EddyPro run are stored in ", path_eddypro_projfiles, 
    "\n", sep = ""
  )
}
