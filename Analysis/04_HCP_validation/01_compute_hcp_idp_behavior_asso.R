#!/usr/bin/env Rscript

# =========================================
# This script computes HCP IDP-behavior association vectors.
#
# For each selected HCP behavioral phenotype, every structural IDP is tested
# against the phenotype using residualized partial correlation. The default
# method is Spearman implemented by rank-transforming the IDP and phenotype and
# then residualizing both on covariates before Pearson correlation of residuals.
#
# The key outputs are results_<phenotype>.csv files. Each file is an IDP-wise
# association vector containing r, p, regression coefficients, and FDR columns;
# these vectors are then aligned with UKB CCA brain loadings by
# 02_compute_bwas_pattern_similarity.py.
# =========================================

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

parse_args <- function(args) {
  cfg <- list(
    idp_csv = "data/CCA结构指标IDP/HCP_idp_sMRI.csv",
    behavior = "data/HCP_YA_behavior/HCP_behavior.csv",
    restricted = "data/HCP_YA_behavior/RESTRICTED_s1200.csv",
    outdir = "results/hcp_idp_cognition_validation",
    method = "spearman",
    min_n = 100,
    covars = "Age_in_Yrs,Gender",
    cognition_fields = "CogFluidComp_Unadj,WM_Task_Acc,VSPLOT_CRTE,ListSort_Unadj,PMAT24_A_CR",
    require_complete_selected_idps = "TRUE",
    idp_subset = "",
    max_idp = 0
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    keyval <- sub("^--", "", a)
    pos <- regexpr("=", keyval, fixed = TRUE)[1]
    if (pos < 0) next
    key <- substr(keyval, 1, pos - 1)
    val <- substr(keyval, pos + 1, nchar(keyval))
    if (key %in% names(cfg)) cfg[[key]] <- val
  }
  cfg$min_n <- as.integer(cfg$min_n)
  cfg$max_idp <- as.integer(cfg$max_idp)
  cfg$require_complete_selected_idps <- toupper(cfg$require_complete_selected_idps) %in% c("TRUE", "T", "1", "YES", "Y")
  cfg
}

split_list <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

drop_invariant_cols <- function(d) {
  if (ncol(d) == 0) return(d)
  keep <- rep(FALSE, ncol(d))
  for (j in seq_len(ncol(d))) {
    x <- d[[j]]
    if (is.factor(x) || is.character(x)) {
      keep[j] <- length(unique(as.character(x[!is.na(x)]))) >= 2
    } else {
      v <- suppressWarnings(var(as.numeric(x), na.rm = TRUE))
      keep[j] <- is.finite(v) && v > 0
    }
  }
  d[, keep, drop = FALSE]
}

complete_cases_cov <- function(d, n) {
  if (is.null(d) || ncol(d) == 0) return(rep(TRUE, n))
  complete.cases(d)
}

residualize <- function(y, d) {
  out <- rep(NA_real_, length(y))
  keep <- is.finite(y) & complete_cases_cov(d, length(y))
  if (sum(keep) < 3) return(out)
  d2 <- drop_invariant_cols(d[keep, , drop = FALSE])
  if (ncol(d2) == 0) {
    out[keep] <- y[keep] - mean(y[keep], na.rm = TRUE)
    return(out)
  }
  fit <- lm(y[keep] ~ ., data = d2)
  out[keep] <- resid(fit)
  out
}

partial_corr_resid <- function(x, y, d, method = "spearman", min_n = 100) {
  ok <- is.finite(x) & is.finite(y) & complete_cases_cov(d, length(x))
  if (sum(ok) < min_n) return(list(n = sum(ok), r = NA_real_, p = NA_real_))
  x2 <- as.numeric(x[ok]); y2 <- as.numeric(y[ok]); d2 <- d[ok, , drop = FALSE]
  if (method == "spearman") {
    x2 <- rank(x2, ties.method = "average")
    y2 <- rank(y2, ties.method = "average")
  }
  rx <- residualize(x2, d2)
  ry <- residualize(y2, d2)
  use <- is.finite(rx) & is.finite(ry)
  if (sum(use) < min_n || length(unique(rx[use])) < 2 || length(unique(ry[use])) < 2) {
    return(list(n = sum(use), r = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(rx[use], ry[use], method = "pearson"))
  list(n = sum(use), r = unname(ct$estimate), p = ct$p.value)
}

linear_regression_assoc <- function(x, y, d, min_n = 100) {
  ok <- is.finite(x) & is.finite(y) & complete_cases_cov(d, length(x))
  n <- sum(ok)
  if (n < min_n) {
    return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  d2 <- drop_invariant_cols(d[ok, , drop = FALSE])
  dat2 <- data.frame(y = as.numeric(y[ok]), score = as.numeric(x[ok]), check.names = FALSE)
  if (ncol(d2) > 0) dat2 <- cbind(dat2, d2)
  if (length(unique(dat2$score)) < 2 || length(unique(dat2$y)) < 2) {
    return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  fit <- lm(y ~ score + ., data = dat2)
  co <- summary(fit)$coefficients
  if (!("score" %in% rownames(co))) {
    return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  list(
    n = n,
    beta = unname(co["score", "Estimate"]),
    se = unname(co["score", "Std. Error"]),
    t = unname(co["score", "t value"]),
    p = unname(co["score", "Pr(>|t|)"])
  )
}

bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (sum(ok) > 0) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

cfg <- parse_args(args)
if (!cfg$method %in% c("spearman", "pearson")) stop("--method must be spearman/pearson")

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
if (is.na(script_path) || !nzchar(script_path)) script_path <- "hcp_idp_cognition_validation.R"
writeLines(
  c(
    paste0("run_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("working_dir: ", getwd()),
    paste0("command: Rscript ", shQuote(script_path), " ", paste(args, collapse = " ")),
    "",
    "[parsed_parameters]",
    paste(names(cfg), unlist(cfg), sep = " = ")
  ),
  con = file.path(cfg$outdir, "run_cmd.txt"),
  useBytes = TRUE
)

cat("Reading files...\n")
idp <- read.csv(cfg$idp_csv, stringsAsFactors = FALSE, check.names = FALSE)
beh <- read.csv(cfg$behavior, stringsAsFactors = FALSE, check.names = FALSE)
res <- read.csv(cfg$restricted, stringsAsFactors = FALSE, check.names = FALSE)

subject_col <- if ("Subject" %in% names(idp)) "Subject" else if ("subject" %in% names(idp)) "subject" else ""
if (!nzchar(subject_col)) stop("IDP csv must include Subject or subject column")
idp$Subject <- as.numeric(idp[[subject_col]])
if (subject_col != "Subject") idp[[subject_col]] <- NULL
idp <- idp[!duplicated(idp$Subject), , drop = FALSE]

idp_names_orig <- setdiff(names(idp), "Subject")
if (nzchar(cfg$idp_subset)) {
  keep_idp <- split_list(cfg$idp_subset)
  idp_names_orig <- intersect(keep_idp, idp_names_orig)
}
if (cfg$max_idp > 0) idp_names_orig <- head(idp_names_orig, cfg$max_idp)
if (length(idp_names_orig) == 0) stop("No IDP columns selected")

behm <- merge(beh, res, by = "Subject", all = FALSE, suffixes = c("_beh", "_res"))
overlap <- intersect(idp_names_orig, names(behm))
idp_name_map <- setNames(idp_names_orig, idp_names_orig)
if (length(overlap) > 0) {
  mapped_overlap <- paste0(overlap, "__idp")
  names(idp)[match(overlap, names(idp))] <- mapped_overlap
  idp_name_map[overlap] <- mapped_overlap
}

idp <- idp[, c("Subject", unname(idp_name_map[idp_names_orig])), drop = FALSE]

if (cfg$require_complete_selected_idps) {
  idp_complete <- complete.cases(idp[, setdiff(names(idp), "Subject"), drop = FALSE])
  cat("IDP-complete subjects retained:", sum(idp_complete), "of", nrow(idp), "\n")
  idp <- idp[idp_complete, , drop = FALSE]
}

dat <- merge(idp, behm, by = "Subject", all = FALSE)

cat("Merged rows:", nrow(dat), "\n")
cat("Unique subjects:", length(unique(dat$Subject)), "\n")

covars <- split_list(cfg$covars)
covars <- covars[covars %in% names(dat)]
for (cv in covars) if (is.character(dat[[cv]])) dat[[cv]] <- as.factor(dat[[cv]])
if (length(covars) > 0) {
  cat("Covariates used:", paste(covars, collapse = ", "), "\n")
} else {
  cat("Covariates used: none\n")
}

phenos <- split_list(cfg$cognition_fields)
phenos <- phenos[phenos %in% names(dat)]
if (length(phenos) == 0) stop("No cognition phenotypes found in merged data")
cat("Phenotypes tested:", paste(phenos, collapse = ", "), "\n")
cat("IDP predictors tested:", length(idp_names_orig), "\n")
if (length(overlap) > 0) {
  cat("Renamed overlapping IDP columns:", length(overlap), "\n")
}

d_cov <- dat[, covars, drop = FALSE]
rows <- vector("list", length(idp_names_orig) * length(phenos))
idx <- 0L
cat("Running association tests...\n")
for (tr in idp_names_orig) {
  tr_col <- unname(idp_name_map[[tr]])
  x <- suppressWarnings(as.numeric(dat[[tr_col]]))
  for (ph in phenos) {
    y <- suppressWarnings(as.numeric(dat[[ph]]))
    st <- partial_corr_resid(x, y, d_cov, method = cfg$method, min_n = cfg$min_n)
    lm_st <- linear_regression_assoc(x, y, d_cov, min_n = cfg$min_n)
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      score_trait = tr,
      phenotype = ph,
      n = st$n,
      r = st$r,
      p = st$p,
      lm_beta = lm_st$beta,
      lm_se = lm_st$se,
      lm_t = lm_st$t,
      lm_p = lm_st$p,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- do.call(rbind, rows)
res_df$p_fdr_by_pheno <- NA_real_
res_df$p_fdr_global <- bh(res_df$p)
res_df$lm_p_fdr_by_pheno <- NA_real_
res_df$lm_p_fdr_global <- bh(res_df$lm_p)
for (ph in unique(res_df$phenotype)) {
  ii <- which(res_df$phenotype == ph)
  res_df$p_fdr_by_pheno[ii] <- bh(res_df$p[ii])
  res_df$lm_p_fdr_by_pheno[ii] <- bh(res_df$lm_p[ii])
}

sig_df <- res_df[is.finite(res_df$p_fdr_by_pheno) & res_df$p_fdr_by_pheno < 0.05, , drop = FALSE]
lm_sig_df <- res_df[is.finite(res_df$lm_p_fdr_by_pheno) & res_df$lm_p_fdr_by_pheno < 0.05, , drop = FALSE]
top_df <- res_df[is.finite(res_df$r), , drop = FALSE]
top_df <- top_df[order(top_df$phenotype, -abs(top_df$r), top_df$p), , drop = FALSE]
top_df <- do.call(rbind, lapply(split(top_df, top_df$phenotype), function(d) head(d, 100)))
lm_top_df <- res_df[is.finite(res_df$lm_beta), , drop = FALSE]
lm_top_df <- lm_top_df[order(lm_top_df$phenotype, -abs(lm_top_df$lm_t), lm_top_df$lm_p), , drop = FALSE]
lm_top_df <- do.call(rbind, lapply(split(lm_top_df, lm_top_df$phenotype), function(d) head(d, 100)))

write.csv(res_df, file.path(cfg$outdir, "idp_cognition_partialcorr_all.csv"), row.names = FALSE)
write.csv(sig_df, file.path(cfg$outdir, "idp_cognition_partialcorr_fdr005_by_pheno.csv"), row.names = FALSE)
write.csv(top_df, file.path(cfg$outdir, "idp_cognition_partialcorr_top100_absr_by_pheno.csv"), row.names = FALSE)
write.csv(res_df, file.path(cfg$outdir, "idp_cognition_regression_all.csv"), row.names = FALSE)
write.csv(lm_sig_df, file.path(cfg$outdir, "idp_cognition_regression_fdr005_by_pheno.csv"), row.names = FALSE)
write.csv(lm_top_df, file.path(cfg$outdir, "idp_cognition_regression_top100_abst_by_pheno.csv"), row.names = FALSE)

for (ph in unique(res_df$phenotype)) {
  sub <- res_df[res_df$phenotype == ph, , drop = FALSE]
  sub <- sub[order(sub$p, -abs(sub$r)), , drop = FALSE]
  write.csv(sub, file.path(cfg$outdir, paste0("results_", ph, ".csv")), row.names = FALSE)
}

cat("Done.\n")
cat("Output:", cfg$outdir, "\n")
cat("All tests:", nrow(res_df), "\n")
cat("FDR<0.05 by phenotype:", nrow(sig_df), "\n")
cat("Regression FDR<0.05 by phenotype:", nrow(lm_sig_df), "\n")
