{
  # Open this saved file in RStudio, place the cursor inside this block,
  # and click Run once without selecting individual lines. Source also works.
  .ecat_root <- local({
    is_release <- function(path) {
      length(path) == 1L && !is.na(path) && nzchar(path) &&
        all(file.exists(file.path(path, c("E_CAT.R", "utils/ecat_core.R",
          "example_data/simulation_01.csv", "example_data/simulation_02.csv"))))
    }
    from_file <- function(path) {
      if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) return(NULL)
      folder <- dirname(normalizePath(path, winslash = "/", mustWork = FALSE))
      if (is_release(folder)) folder else NULL
    }
    # source() provides the path in an evaluation frame.
    for (frame in rev(sys.frames())) {
      folder <- from_file(frame$ofile)
      if (!is.null(folder)) return(folder)
    }
    # RStudio Run sends code to the console without a source-file frame.
    if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") &&
        !requireNamespace("rstudioapi", quietly = TRUE)) {
      tryCatch(install.packages("rstudioapi", repos = "https://cloud.r-project.org"),
               error = function(e) message("Editor path lookup unavailable: ", conditionMessage(e)))
    }
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      editor <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                         error = function(e) NULL)
      folder <- from_file(editor)
      if (!is.null(folder)) return(folder)
    }
    if (!interactive()) {
      arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
      if (length(arg)) {
        folder <- from_file(sub("^--file=", "", arg[1]))
        if (!is.null(folder)) return(folder)
      }
    }
    # Also support pasted code and Run with the repository as working directory.
    if (is_release(getwd())) return(normalizePath(getwd(), winslash = "/"))
    stop("Cannot locate E-CAT files. Open the saved example.R in its release folder, ",
         "or set the working directory to that folder. Keep both subfolders together.")
  })
  if (!requireNamespace("energy", quietly = TRUE)) {
    message("Installing the required R package: energy")
    tryCatch(install.packages("energy", repos = "https://cloud.r-project.org"),
             error = function(e) message(conditionMessage(e)))
    if (!requireNamespace("energy", quietly = TRUE))
      stop("Package installation failed. Check your internet connection and run install.packages('energy').")
  }
  source(file.path(.ecat_root, "E_CAT.R"), local = TRUE)
  examples <- data.frame(filename = c("simulation_01.csv", "simulation_02.csv"),
                         n = c(200L, 500L), strong = c(2L, 4L),
                         weak = c(2L, 3L), invalid = c(2L, 7L), K = c(2L, 3L))
  ecat_results <- setNames(vector("list", nrow(examples)), examples$filename)
  for (i in seq_len(nrow(examples))) {
    filename <- examples$filename[i]
    data <- read.csv(file.path(.ecat_root, "example_data", filename))
    result <- ecat(data, K = examples$K[i])
    ecat_results[[filename]] <- result
  }
  ecat_summary <- data.frame(
    filename = examples$filename,
    n = examples$n,
    K = examples$K,
    beta_hat = vapply(ecat_results, function(x) x$beta_hat, numeric(1)),
    true_hat = 1,
    row.names = NULL
  )
  # One fixed sample per design: MSE equals the single squared error.
  ecat_summary$mse <- (ecat_summary$beta_hat - ecat_summary$true_hat)^2
  ecat_summary$mean_runtime_sec <- vapply(ecat_results, function(x) x$runtime_sec, numeric(1))
  print(ecat_summary, row.names = FALSE)
}
