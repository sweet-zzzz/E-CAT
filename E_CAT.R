# Public E-CAT interface and command-line entry point.
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
.ecat_engine <- new.env(parent = environment())
sys.source(file.path(.ecat_root, "utils", "ecat_core.R"), envir = .ecat_engine)

# Fit E-CAT to a numeric data frame. Column names may be arbitrary.
ecat <- function(data, instruments = grep("^IV", names(data), value = TRUE),
                 treatment = "Treatment", outcome = "Outcome",
                 covariates = character(0), K = 2L, lambda_moment = 1,
                 tuning.1st = 1, alpha = 0.05, max_combs = 2000L,
                 verbose = FALSE) {
  if (!requireNamespace("energy", quietly = TRUE))
    stop("Install the required package with install.packages('energy').")
  start_time <- proc.time()
  if (!is.data.frame(data) || anyDuplicated(names(data)))
    stop("data must be a data frame with unique column names.")
  cols <- c(treatment, outcome, instruments, covariates)
  if (length(treatment) != 1L || length(outcome) != 1L ||
      length(instruments) < 2L || anyDuplicated(cols) ||
      !all(cols %in% names(data))) stop("Specify distinct, existing columns and at least two instruments.")
  if (!all(vapply(data[cols], is.numeric, logical(1))) ||
      !all(is.finite(as.matrix(data[cols]))))
    stop("All selected columns must be numeric and finite; handle missing values first.")
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(K) || K != as.integer(K) || K < 2 || K > length(instruments))
    stop("K must be an integer between 2 and the number of instruments.")
  if (!scalar(lambda_moment) || lambda_moment < 0 ||
      !scalar(tuning.1st) || tuning.1st <= 0 ||
      !scalar(alpha) || alpha <= 0 || alpha >= 1 ||
      !scalar(max_combs) || max_combs < 1 || max_combs != floor(max_combs))
    stop("Invalid scoring, screening, confidence-level, or search-limit parameter.")
  Z <- data[instruments]
  names(Z) <- paste0("IV", seq_along(Z))
  W <- if (length(covariates)) data[covariates] else NULL
  if (!is.null(W)) names(W) <- paste0("W", seq_along(W))
  if (nrow(data) <= ncol(Z) + length(covariates) + 1L)
    stop("Too few observations for the supplied instruments and covariates.")
  if (!is.null(W) && qr(cbind(1, as.matrix(W)))$rank < ncol(W) + 1L)
    stop("Covariates and intercept must be linearly independent.")
  prep <- .ecat_engine$residualize_by_covariates_ecat(data[[treatment]], data[[outcome]], Z, W)
  if (qr(as.matrix(prep$IV_set))$rank < ncol(Z))
    stop("Instrument matrix must have full column rank after covariate adjustment.")
  fit <- .ecat_engine$ECAT_JointRank_TSHT(
    data[[treatment]], data[[outcome]], Covariates = W, IV_set = Z,
    K = K, tuning.1st = tuning.1st, alpha = alpha, verbose = verbose,
    distance = "dcorr", lambda_moment = lambda_moment, max_combs = max_combs)
  if (!is.finite(fit$beta_hat) ||
      (!is.null(fit$joint_info) && !is.finite(fit$joint_info$best_score)))
    stop("E-CAT could not obtain a finite estimate or candidate score.")
  fit$selected_instruments <- instruments[fit$Valid_IV_subset]
  fit$screened_instruments <- instruments[fit$Strong_IV_set]
  fit$confidence_interval <- as.numeric(fit$estimate_detail$ci)
  fit$n <- nrow(data)
  if (!is.null(fit$joint_info)) {
    fit$joint_info$candidates$instrument_names <- vapply(
      fit$joint_info$candidates$subset,
      function(idx) paste(fit$screened_instruments[idx], collapse = ","), character(1))
  }
  elapsed <- proc.time() - start_time
  fit$runtime_sec <- unname(as.numeric(elapsed["elapsed"]))
  fit$cpu_time_sec <- unname(as.numeric(elapsed["user.self"] + elapsed["sys.self"]))
  fit
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input <- if (length(args)) args[1] else file.path(.ecat_root, "example_data", "simulation_01.csv")
  k <- if (length(args) >= 2L) as.numeric(args[2]) else 2L
  fit <- ecat(read.csv(input, check.names = FALSE), K = k)
  cat("Selected instruments:", paste(fit$selected_instruments, collapse = ", "), "\n")
  cat("Effect estimate:", format(fit$beta_hat, digits = 8), "\n")
  cat("Runtime (seconds):", fit$runtime_sec, "\n")
  cat("CPU time (seconds):", fit$cpu_time_sec, "\n")
  cat("Nominal 95% interval:", paste(format(fit$confidence_interval, digits = 8), collapse = " to "), "\n")
}
