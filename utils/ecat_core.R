# E-CAT numerical routines.

normalize_distance_matrix_ecat <- function(Dmat) {
  vals <- Dmat[upper.tri(Dmat)]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(Dmat)

  min_val <- min(vals)
  max_val <- max(vals)
  if (!is.finite(min_val) || !is.finite(max_val) || max_val <= min_val) {
    diag(Dmat) <- 0
    return(Dmat)
  }

  Dnorm <- Dmat
  off_diag <- upper.tri(Dnorm) | lower.tri(Dnorm)
  Dnorm[off_diag] <- (Dnorm[off_diag] - min_val) / (max_val - min_val)
  diag(Dnorm) <- 0
  Dnorm
}

residualize_by_covariates_ecat <- function(Treatment, Outcome, IV_set, Covariates = NULL) {
  if (is.null(Covariates)) {
    return(list(
      Treatment = as.numeric(Treatment),
      Outcome = as.numeric(Outcome),
      IV_set = as.data.frame(IV_set)
    ))
  }

  data <- data.frame(Treatment = Treatment, Outcome = Outcome, Covariates, IV_set)
  covariate_names <- names(Covariates)
  target_names <- setdiff(names(data), covariate_names)

  residuals_df <- data.frame(matrix(ncol = length(target_names), nrow = nrow(data)))
  colnames(residuals_df) <- target_names

  for (col_name in target_names) {
    formula <- as.formula(paste(col_name, "~", paste(covariate_names, collapse = "+")))
    fit <- lm(formula, data = data)
    residuals_df[[col_name]] <- residuals(fit)
  }

  list(
    Treatment = as.numeric(residuals_df$Treatment),
    Outcome = as.numeric(residuals_df$Outcome),
    IV_set = residuals_df[, grep("^IV", names(residuals_df)), drop = FALSE]
  )
}

strong_iv_screen_ecat <- function(Treatment, IV_set, tuning.1st = 2.01) {
  Z <- as.matrix(IV_set)
  D <- as.numeric(Treatment)
  n <- length(D)
  pz <- ncol(Z)

  if (pz == 0) stop("IV_set has no columns.")

  W <- Z
  covW <- crossprod(W) / n
  U <- solve(covW)
  WUMat <- W %*% U

  qrW <- qr(W)
  ITT_D <- qr.coef(qrW, D)
  resid_D <- qr.resid(qrW, D)
  SigmaSqD <- sum(resid_D^2) / (n - ncol(W))
  se_vec <- sqrt(SigmaSqD * colSums(WUMat^2) / n)

  tuning_val <- if (is.null(tuning.1st)) 2.01 else tuning.1st
  threshold_vec <- se_vec * sqrt(tuning_val * log(pz) / n)
  SHat <- which(abs(ITT_D) >= threshold_vec)

  if (length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. Defaulting to all IVs.")
    SHat <- seq_len(pz)
  }

  list(
    SHat = SHat,
    ITT_D = ITT_D,
    se_vec = se_vec,
    threshold = threshold_vec,
    tuning = tuning_val,
    rule = 'TSHT_default_first_stage'
  )
}

estimate_TSHT_ecat <- function(Y, D, Z, VHat, X = NULL, intercept = FALSE, alpha = 0.05) {
  if (!is.null(X)) {
    W <- cbind(Z, X)
  } else {
    W <- Z
  }

  if (intercept) {
    W <- cbind(W, 1)
  }

  n <- length(Y)
  pz <- ncol(Z)

  covW <- t(W) %*% W / n
  WUMat <- W %*% solve(covW)[, 1:pz, drop = FALSE]
  qrW <- qr(W)

  ITT_Y <- qr.coef(qrW, Y)[1:pz]
  ITT_D <- qr.coef(qrW, D)[1:pz]
  SigmaSqY <- sum(qr.resid(qrW, Y)^2) / (n - ncol(W))
  SigmaSqD <- sum(qr.resid(qrW, D)^2) / (n - ncol(W))
  SigmaYD <- sum(qr.resid(qrW, Y) * qr.resid(qrW, D)) / (n - ncol(W))

  A <- t(WUMat) %*% WUMat / n
  AVHat <- solve(A[VHat, VHat, drop = FALSE])

  betaHat <- (t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) /
    (t(ITT_D[VHat]) %*% AVHat %*% ITT_D[VHat])

  SigmaSq <- SigmaSqY + betaHat^2 * SigmaSqD - 2 * betaHat * SigmaYD
  betaVarHat <- SigmaSq *
    (t(ITT_D[VHat]) %*% AVHat %*% (t(WUMat) %*% WUMat / n)[VHat, VHat, drop = FALSE] %*% AVHat %*% ITT_D[VHat]) /
    (t(ITT_D[VHat]) %*% AVHat %*% ITT_D[VHat])^2

  ci <- c(betaHat - qnorm(1 - alpha / 2) * sqrt(betaVarHat / n),
          betaHat + qnorm(1 - alpha / 2) * sqrt(betaVarHat / n))

  list(betaHat = betaHat, betaVarHat = betaVarHat, ci = ci)
}

pairwise_cat_distance_ecat <- function(IV1, IV2, Treatment, Outcome, distance = c('dcorr', 'corr', 'effect')) {
  distance <- match.arg(distance)

  beta1 <- cov(Outcome, IV1) / cov(Treatment, IV1)
  beta2 <- cov(Outcome, IV2) / cov(Treatment, IV2)

  if (distance == 'effect') {
    if (!is.finite(beta1) || !is.finite(beta2)) return(Inf)
    return((beta1 - beta2)^2)
  }

  A1 <- Outcome - beta1 * Treatment
  A2 <- Outcome - beta2 * Treatment

  if (distance == 'dcorr') {
    return(energy::dcor(A1, IV2) + energy::dcor(A2, IV1))
  }

  corr_1 <- suppressWarnings(stats::cor(A1, IV2, use = 'complete.obs'))
  corr_2 <- suppressWarnings(stats::cor(A2, IV1, use = 'complete.obs'))
  corr_1 <- ifelse(is.finite(corr_1), abs(corr_1), 1)
  corr_2 <- ifelse(is.finite(corr_2), abs(corr_2), 1)
  corr_1 + corr_2
}

build_cat_distance_matrix_ecat <- function(IV_set, Treatment, Outcome, distance = c('dcorr', 'corr', 'effect'), normalize = TRUE) {
  distance <- match.arg(distance)
  Z <- as.data.frame(IV_set)
  p <- ncol(Z)

  if (p == 1) {
    mat <- matrix(0, 1, 1)
    colnames(mat) <- colnames(Z)
    rownames(mat) <- colnames(Z)
    return(mat)
  }

  Dmat <- matrix(0, nrow = p, ncol = p)
  colnames(Dmat) <- colnames(Z)
  rownames(Dmat) <- colnames(Z)

  for (i in seq_len(p)) {
    for (j in seq(i, p)) {
      if (i == j) {
        Dmat[i, j] <- 0
      } else {
        dij <- pairwise_cat_distance_ecat(Z[[i]], Z[[j]], Treatment, Outcome, distance = distance)
        Dmat[i, j] <- dij
        Dmat[j, i] <- dij
      }
    }
  }

  if (normalize) {
    Dmat <- normalize_distance_matrix_ecat(Dmat)
  }

  Dmat
}

subset_score_from_dmat_ecat <- function(subset_idx, Dmat) {
  subset_idx <- as.integer(subset_idx)
  subset_idx <- subset_idx[!is.na(subset_idx)]
  if (length(subset_idx) <= 1) return(0)
  submat <- Dmat[subset_idx, subset_idx, drop = FALSE]
  mean(submat[upper.tri(submat)])
}

select_best_subset_within_group_ecat <- function(group_idx, Dmat, K) {
  score_one_subset <- function(subset_global) {
    subset_score_from_dmat_ecat(subset_global, Dmat)
  }

  if (length(group_idx) <= K) {
    score <- score_one_subset(group_idx)
    return(list(best_subset_global = group_idx, best_score = score))
  }

  local_combs <- combn(length(group_idx), K)
  scores <- numeric(ncol(local_combs))
  for (c in seq_len(ncol(local_combs))) {
    local_idx <- local_combs[, c]
    subset_global <- group_idx[local_idx]
    scores[c] <- score_one_subset(subset_global)
  }
  best_col <- which.min(scores)
  list(best_subset_global = group_idx[local_combs[, best_col]], best_score = scores[best_col])
}

tsht_subset_moment_objective_ecat <- function(IV_set, subset_i, Outcome, Treatment) {
  subset_i <- as.integer(subset_i)
  subset_i <- subset_i[!is.na(subset_i)]
  if (length(subset_i) == 0) return(Inf)

  Z_all <- as.matrix(IV_set)
  Z <- as.matrix(IV_set[, subset_i, drop = FALSE])
  Y <- as.numeric(Outcome)
  D <- as.numeric(Treatment)

  beta_hat <- tryCatch(
    as.numeric(estimate_TSHT_ecat(as.matrix(Y), as.matrix(D), Z_all, subset_i)$betaHat),
    error = function(e) NA_real_
  )
  if (!is.finite(beta_hat)) return(Inf)

  resid <- Y - beta_hat * D
  moments <- colMeans(Z * resid)
  z_scale <- apply(Z, 2, stats::sd)
  z_scale[!is.finite(z_scale) | z_scale <= 0] <- 1
  resid_scale <- stats::sd(resid)
  if (!is.finite(resid_scale) || resid_scale <= 0) resid_scale <- 1
  scaled_moments <- moments / (z_scale * resid_scale)

  as.numeric(mean(scaled_moments^2))
}

rank01_ecat <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  if (sum(ok) == 1) {
    out[ok] <- 0
    return(out)
  }
  out[ok] <- (rank(x[ok], ties.method = 'average') - 1) / (sum(ok) - 1)
  out
}

collect_all_subsets_ecat <- function(Dmat, K, max_combs = 2000) {
  p <- nrow(Dmat)
  if (p < K) return(data.frame())
  cmb <- combn(p, K)
  if (ncol(cmb) > max_combs) {
    scores <- apply(cmb, 2, subset_score_from_dmat_ecat, Dmat = Dmat)
    cmb <- cmb[, order(scores)[seq_len(max_combs)], drop = FALSE]
  }
  rows <- vector('list', ncol(cmb))
  for (i in seq_len(ncol(cmb))) {
    subset_i <- cmb[, i]
    rows[[i]] <- data.frame(
      subset_key = paste(subset_i, collapse = ','),
      subset = I(list(subset_i)),
      distance_score = subset_score_from_dmat_ecat(subset_i, Dmat),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

choose_joint_rank_tsht_subset_ecat <- function(candidate_subsets,
                                               IV_set,
                                               Outcome,
                                               Treatment,
                                               lambda_moment = 1) {
  if (nrow(candidate_subsets) == 0) {
    return(list(best_subset_global = integer(0), best_score = Inf,
                candidates = candidate_subsets, lambda_moment = lambda_moment))
  }

  moment_q <- rep(NA_real_, nrow(candidate_subsets))
  for (i in seq_len(nrow(candidate_subsets))) {
    moment_q[i] <- tryCatch(
      tsht_subset_moment_objective_ecat(IV_set, candidate_subsets$subset[[i]], Outcome, Treatment),
      error = function(e) NA_real_
    )
  }

  candidate_subsets$tsht_moment_objective <- moment_q
  candidate_subsets$distance_rank01 <- rank01_ecat(candidate_subsets$distance_score)
  candidate_subsets$moment_rank01 <- rank01_ecat(candidate_subsets$tsht_moment_objective)
  candidate_subsets$joint_rank_score <- candidate_subsets$distance_rank01 +
    lambda_moment * candidate_subsets$moment_rank01

  ok <- which(is.finite(candidate_subsets$joint_rank_score))
  if (length(ok) == 0) {
    return(list(best_subset_global = candidate_subsets$subset[[1]], best_score = Inf,
                candidates = candidate_subsets, lambda_moment = lambda_moment))
  }
  best_row <- ok[which.min(candidate_subsets$joint_rank_score[ok])]

  list(
    best_subset_global = candidate_subsets$subset[[best_row]],
    best_score = candidate_subsets$joint_rank_score[best_row],
    candidates = candidate_subsets,
    best_candidate_row = best_row,
    lambda_moment = lambda_moment
  )
}

ECAT_JointRank_TSHT <- function(Treatment,
                                Outcome,
                                Covariates = NULL,
                                IV_set,
                                K,
                                tuning.1st = NULL,
                                alpha = 0.05,
                                verbose = TRUE,
                                distance = c('dcorr', 'corr', 'effect'),
                                lambda_moment = 1,
                                max_combs = 2000) {
  distance <- match.arg(distance)
  prep <- residualize_by_covariates_ecat(Treatment, Outcome, IV_set, Covariates)
  Treatment_r <- prep$Treatment
  Outcome_r <- prep$Outcome
  IV_set_r <- prep$IV_set

  screen <- strong_iv_screen_ecat(Treatment_r, IV_set_r, tuning.1st = tuning.1st)
  SHat <- screen$SHat
  joint_info <- NULL

  if (length(SHat) <= K) {
    if (verbose) cat('|S_hat| <= K, directly using S_hat.\n')
    chosen_group_global <- SHat
    chosen_global <- SHat
  } else {
    strong_IV_set <- IV_set_r[, SHat, drop = FALSE]
    distance_Dmat <- build_cat_distance_matrix_ecat(
      strong_IV_set, Treatment_r, Outcome_r, distance = distance, normalize = TRUE
    )
    candidate_subsets <- collect_all_subsets_ecat(distance_Dmat, K = K, max_combs = max_combs)
    joint_info <- choose_joint_rank_tsht_subset_ecat(
      candidate_subsets = candidate_subsets,
      IV_set = strong_IV_set,
      Outcome = Outcome_r,
      Treatment = Treatment_r,
      lambda_moment = lambda_moment
    )

    if (length(joint_info$best_subset_global) < K || all(is.na(joint_info$best_subset_global))) {
      fallback_info <- select_best_subset_within_group_ecat(seq_along(SHat), Dmat = distance_Dmat, K = K)
      chosen_global <- SHat[fallback_info$best_subset_global]
    } else {
      chosen_global <- SHat[joint_info$best_subset_global]
    }
    chosen_group_global <- SHat
  }

  est_detail <- estimate_TSHT_ecat(as.matrix(Outcome_r), as.matrix(Treatment_r), as.matrix(IV_set_r), chosen_global, alpha = alpha)
  beta_hat <- est_detail$betaHat

  if (verbose) {
    cat('Strong IV set (S_hat):', SHat, '\n')
    cat('Candidate pool:', chosen_group_global, '\n')
    cat('Final E-CAT-JointRank-TSHT subset:', chosen_global, '\n')
    cat('Estimated Beta (E-CAT-JointRank-TSHT):', beta_hat, '\n')
  }

  list(
    beta_hat = as.numeric(beta_hat),
    Valid_IV_subset = chosen_global,
    Strong_IV_set = SHat,
    selected_cluster = chosen_group_global,
    joint_info = joint_info,
    estimate_detail = est_detail,
    K = K,
    distance = distance,
    method = 'E-CAT-JointRank-TSHT',
    lambda_moment = lambda_moment
  )
}
