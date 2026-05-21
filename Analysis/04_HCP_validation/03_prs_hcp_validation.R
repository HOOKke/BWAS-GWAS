#!/usr/bin/env Rscript

# =========================================
# This script performs the HCP PRS external-validation analysis for Fig.8.
#
# Workflow:
#   1. Read mode-specific PRS score files across p-value thresholds.
#   2. Merge PRS scores with HCP behavior, restricted phenotype, and covariate
#      tables.
#   3. Residualize PRS and behavior values for age, sex/gender, handedness,
#      race, and genetic PCs 1-10.
#   4. Estimate adjusted partial Spearman associations for each PRS threshold.
#   5. Select the best threshold per PRS-phenotype pair.
#   6. Group participants into quantiles using residualized PRS scores and plot
#      mean residualized behavior values in each quantile.
#
# The compact GitHub release smoke test uses code/validation and the released
# Source Data workbook. This script is retained as the upstream full-data
# analysis reference and requires the original local HCP/PRS input files.
# =========================================

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

parse_args <- function(args) {
  # Default file paths mirror the full working repository. They can be
  # overridden with --key=value command-line options.
  cfg <- list(
    prs1 = "data/hcp_prs/IDP_PRS/IDP1_PRS_all.csv",
    prs2 = "data/hcp_prs/IDP_PRS/IDP2_PRS_all.csv",
    behavior = "data/HCP_YA_behavior/HCP_behavior.csv",
    restricted = "data/HCP_YA_behavior/RESTRICTED_s1200.csv",
    covar_csv = "data/HCP_YA_behavior/hcp_covariates_with_genePC.csv",
    mode1_list_txt = "results/BWAS层面验证/IDP_Bhv/reusable_hcp_phenotype_lists/mode1_hcp_bwas_plot.txt",
    mode2_list_txt = "results/BWAS层面验证/IDP_Bhv/reusable_hcp_phenotype_lists/mode2_hcp_bwas_plot.txt",
    outdir = "results/hcp_prs_response_letter_validation_20260501",
    method = "spearman",
    min_n = 100,
    covars = "Age_in_Yrs,Gender,Handedness,Race,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
    n_groups = 4,
    dpi = 320,
    point_size = 3.3,
    err_width = 0.14,
    panel_width = 3.9,
    panel_height = 3.9,
    panel_facet_ncol = 4,
    font_family = "Hiragino Sans GB"
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) next
    k <- kv[1]
    v <- kv[2]
    if (k %in% names(cfg)) cfg[[k]] <- v
  }
  cfg$min_n <- as.integer(cfg$min_n)
  cfg$n_groups <- as.integer(cfg$n_groups)
  cfg$dpi <- as.integer(cfg$dpi)
  cfg$point_size <- as.numeric(cfg$point_size)
  cfg$err_width <- as.numeric(cfg$err_width)
  cfg$panel_width <- as.numeric(cfg$panel_width)
  cfg$panel_height <- as.numeric(cfg$panel_height)
  cfg$panel_facet_ncol <- as.integer(cfg$panel_facet_ncol)
  cfg
}

drop_invariant_cols <- function(d) {
  # Remove constant covariates before fitting residualization models.
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

residualize <- function(y, d) {
  # Return residuals from y ~ covariates while preserving original row order.
  out <- rep(NA_real_, length(y))
  keep <- is.finite(y) & complete.cases(d)
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
  # Adjust x and y for covariates, then correlate residuals. For Spearman
  # associations, rank-transform x and y before residualization.
  ok <- is.finite(x) & is.finite(y) & complete.cases(d)
  if (sum(ok) < min_n) return(list(n = sum(ok), r = NA_real_, p = NA_real_))
  x2 <- as.numeric(x[ok])
  y2 <- as.numeric(y[ok])
  d2 <- d[ok, , drop = FALSE]
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

safe_quantile_groups <- function(x, n_groups) {
  # Bin residualized PRS scores into quantiles; duplicated breaks are collapsed.
  q <- unique(quantile(x, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE, type = 7))
  if (length(q) < 3) return(rep(NA_integer_, length(x)))
  as.integer(cut(x, breaks = q, include.lowest = TRUE, labels = FALSE))
}

bh <- function(p) {
  # Benjamini-Hochberg helper used for exploratory FDR summaries.
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (sum(ok) > 0) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

sanitize_filename <- function(x) {
  # Convert phenotype names into filesystem-safe plot filenames.
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "plot"
  x
}

read_list_txt <- function(path) {
  # Read a whitespace-delimited phenotype list file.
  vals <- scan(path, what = "", quiet = TRUE)
  vals <- trimws(vals)
  unique(vals[nzchar(vals)])
}

display_label <- function(x) {
  # Map HCP variable names to final figure-facing labels.
  label_map <- c(
    DSM_Anxi_Raw = "DSM anxiety",
    DSM_Avoid_Raw = "DSM avoidance",
    DSM_Depr_Raw = "DSM depression",
    DSM_Somp_Raw = "DSM somatic problems",
    FearAffect_Unadj = "Fear affect",
    FearSomat_Unadj = "Fear somatic",
    Loneliness_Unadj = "Loneliness",
    PercHostil_Unadj = "Perceived hostility",
    PercReject_Unadj = "Perceived rejection",
    PercStress_Unadj = "Perceived stress",
    SSAGA_Agoraphobia = "SSAGA agoraphobia",
    SSAGA_Depressive_Ep = "SSAGA depressive episode",
    SSAGA_Depressive_Sx = "SSAGA depressive symptoms",
    SSAGA_PanicDisorder = "SSAGA panic disorder",
    ASR_Attn_Raw = "ASR attention problems",
    ASR_Extn_Raw = "ASR externalizing",
    ASR_Thot_Raw = "ASR thought problems",
    DSM_Adh_Raw = "DSM ADHD",
    IWRD_TOT = "Verbal episodic memory",
    PMAT24_A_CR = "Fluid intelligence",
    ProcSpeed_Unadj = "Processing speed",
    VSPLOT_TC = "Spatial orientation accuracy",
    AngAffect_Unadj = "Anger affect",
    AngAggr_Unadj = "Anger / aggression",
    DSM_Hype_Raw = "DSM hyperactivity",
    DDisc_AUC_200 = "Delay discounting AUC 200"
  )
  out <- unname(label_map[x])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- x[miss]
  out
}

compact_covar_note <- function(covars) {
  # Create a concise covariate label for figure captions and source-data tables.
  if (identical(
    covars,
    c("Age_in_Yrs", "Gender", "Handedness", "Race", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
  )) {
    return("Age, Gender, Handedness, Race, PC1-10")
  }
  paste(covars, collapse = ", ")
}

save_plot_triplet <- function(plot_obj, out_png, out_pdf, out_svg, width, height, cfg, dpi) {
  # Save each figure in PNG, PDF, and SVG formats for manuscript assembly.
  ggsave(
    out_png,
    plot_obj,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    out_pdf,
    plot_obj,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white",
    family = cfg$font_family
  )
  ggsave(
    out_svg,
    plot_obj,
    device = svglite::svglite,
    width = width,
    height = height,
    bg = "white",
    system_fonts = list(sans = cfg$font_family)
  )
}

plot_one_pair <- function(subdf, out_png, out_pdf, out_svg, cfg) {
  # Render one phenotype-specific PRS quantile panel. The y-axis is the mean
  # residualized behavioral phenotype value within each residualized-PRS bin.
  suppressWarnings(suppressMessages(library(ggplot2)))
  first <- subdf[1, , drop = FALSE]
  mode_name <- first$mode
  col_fill <- if (mode_name == "mode1") "#E74C3C" else "#4A90E2"
  col_line <- if (mode_name == "mode1") "#C0392B" else "#2E6FD8"
  pt_shape <- if (mode_name == "mode1") 21 else 25
  x_lab <- if (mode_name == "mode1") "HCP PRS_mode1 quantile" else "HCP PRS_mode2 quantile"

  p <- ggplot(subdf, aes(x = quantile_group, y = mean)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.35, color = "#8A8F98") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = cfg$err_width, linewidth = 0.55, color = "black") +
    geom_point(shape = pt_shape, size = cfg$point_size, stroke = 0.55, fill = col_fill, color = col_line) +
    scale_x_continuous(breaks = 1:cfg$n_groups, limits = c(0.8, cfg$n_groups + 0.2)) +
    labs(
      title = first$display_label,
      subtitle = paste0(
        if (mode_name == "mode1") "PRSmode1" else "PRSmode2",
        " | best p-threshold = ", first$pthres,
        " | partial r = ", sprintf("%.3f", first$r),
        " | p = ", sprintf("%.2g", first$p),
        "\nBest-threshold FDR = ", sprintf("%.2g", first$p_fdr_best_targeted),
        " | n = ", first$n_plot
      ),
      x = x_lab,
      y = "Behavioral phenotype (residualized)",
      caption = paste0("Covariates: ", first$covar_note)
    ) +
    theme_classic(base_size = 11.5, base_family = cfg$font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 14.8),
      plot.subtitle = element_text(size = 9.8, lineheight = 0.96),
      plot.caption = element_text(size = 8.9, hjust = 0),
      axis.title = element_text(size = 12.4),
      axis.title.y = element_text(size = 12.4, margin = margin(r = 16)),
      axis.text = element_text(size = 11.0),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 30)
    )

  save_plot_triplet(p, out_png, out_pdf, out_svg, cfg$panel_width, cfg$panel_height, cfg, cfg$dpi)
}

plot_contact_sheet <- function(plot_df, mode_name, out_png, out_pdf, out_svg, cfg) {
  # Render the four-panel mode-level figure used in the final Fig.8 display.
  suppressWarnings(suppressMessages(library(ggplot2)))
  sub <- plot_df[plot_df$mode == mode_name, , drop = FALSE]
  if (nrow(sub) == 0) return(invisible(NULL))
  sub$pair_label <- factor(sub$pair_label, levels = unique(sub$pair_label))
  fill_map <- c(mode1 = "#E74C3C", mode2 = "#4A90E2")
  color_map <- c(mode1 = "#C0392B", mode2 = "#2E6FD8")
  shape_map <- c(mode1 = 21, mode2 = 25)

  p <- ggplot(sub, aes(x = quantile_group, y = mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#8A8F98") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = cfg$err_width, linewidth = 0.38, color = "black") +
    geom_point(aes(fill = mode, color = mode, shape = mode), size = 2.4, stroke = 0.45) +
    facet_wrap(~pair_label, scales = "free_y", ncol = cfg$panel_facet_ncol) +
    scale_fill_manual(values = fill_map) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_x_continuous(breaks = 1:cfg$n_groups, limits = c(0.8, cfg$n_groups + 0.2)) +
    labs(
      x = if (mode_name == "mode1") "HCP PRS_mode1 quantile" else "HCP PRS_mode2 quantile",
      y = "Behavioral phenotype (residualized)"
    ) +
    theme_classic(base_size = 10.5, base_family = cfg$font_family) +
    theme(
      strip.text = element_text(size = 11.4, face = "bold"),
      legend.position = "none",
      axis.title = element_text(size = 12.4),
      axis.title.y = element_text(size = 12.4, margin = margin(r = 16)),
      axis.text = element_text(size = 10.7),
      panel.spacing.x = grid::unit(0.34, "in"),
      panel.spacing.y = grid::unit(0.26, "in"),
      plot.margin = margin(t = 12, r = 14, b = 12, l = 26)
    )

  n_panels <- length(unique(sub$pair_label))
  n_rows <- ceiling(n_panels / cfg$panel_facet_ncol)
  sheet_width <- cfg$panel_width * cfg$panel_facet_ncol + 0.8
  sheet_height <- cfg$panel_height * n_rows + 0.55
  save_plot_triplet(p, out_png, out_pdf, out_svg, sheet_width, sheet_height, cfg, cfg$dpi)
}

cfg <- parse_args(args)
if (!cfg$method %in% c("spearman", "pearson")) stop("--method must be spearman/pearson")
if (cfg$n_groups < 2) stop("--n_groups must be >= 2")

# Create output folders and record the exact command-line configuration.
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots", "mode1"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots", "mode2"), recursive = TRUE, showWarnings = FALSE)

script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
if (is.na(script_path) || !nzchar(script_path)) script_path <- "hcp_prs_response_letter_validation.R"
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

mode1_list <- read_list_txt(cfg$mode1_list_txt)
mode2_list <- read_list_txt(cfg$mode2_list_txt)
writeLines(mode1_list, file.path(cfg$outdir, "mode1_validation_list.txt"))
writeLines(mode2_list, file.path(cfg$outdir, "mode2_validation_list.txt"))

list_df <- data.frame(
  mode = c(rep("mode1", length(mode1_list)), rep("mode2", length(mode2_list))),
  phenotype = c(mode1_list, mode2_list),
  stringsAsFactors = FALSE
)
list_df$phenotype_order <- seq_len(nrow(list_df))
write.csv(list_df, file.path(cfg$outdir, "validation_phenotype_lists.csv"), row.names = FALSE)

# Load PRS scores, HCP behavior variables, restricted variables, and covariates.
prs1 <- read.csv(cfg$prs1, stringsAsFactors = FALSE, check.names = FALSE)
prs2 <- read.csv(cfg$prs2, stringsAsFactors = FALSE, check.names = FALSE)
beh <- read.csv(cfg$behavior, stringsAsFactors = FALSE, check.names = FALSE)
res <- read.csv(cfg$restricted, stringsAsFactors = FALSE, check.names = FALSE)
cov <- read.csv(cfg$covar_csv, stringsAsFactors = FALSE, check.names = FALSE)

prs1s <- prs1[, c("FID", "pthres", "SCORE"), drop = FALSE]
prs2s <- prs2[, c("FID", "pthres", "SCORE"), drop = FALSE]
names(prs1s) <- c("Subject", "pthres", "PRS_IDP1")
names(prs2s) <- c("Subject", "pthres", "PRS_IDP2")
prs <- merge(prs1s, prs2s, by = c("Subject", "pthres"), all = FALSE)
prs$Subject <- as.numeric(prs$Subject)

behm <- merge(beh, res, by = "Subject", all = FALSE, suffixes = c("_beh", "_res"))
dat <- merge(prs, behm, by = "Subject", all = FALSE)
dat <- merge(dat, cov, by = "Subject", all = FALSE, suffixes = c("", "_cov"))

# Restrict to covariates present in the merged analysis table and encode
# categorical covariates as factors for lm().
covars <- trimws(unlist(strsplit(cfg$covars, ",", fixed = TRUE)))
covars <- covars[covars != ""]
covars <- covars[covars %in% names(dat)]
if (length(covars) == 0) stop("No valid covariates found.")
for (cv in covars) if (is.character(dat[[cv]])) dat[[cv]] <- as.factor(dat[[cv]])

all_needed <- unique(c(mode1_list, mode2_list))
missing_vars <- setdiff(all_needed, names(dat))
writeLines(
  c(
    if (length(missing_vars) == 0) "All requested phenotypes are present in merged HCP data."
    else paste("Missing requested phenotypes:", paste(missing_vars, collapse = ", ")),
    paste("Covariates:", paste(covars, collapse = ", "))
  ),
  file.path(cfg$outdir, "method_notes.txt")
)
if (length(missing_vars) > 0) stop("Missing requested phenotypes in merged data: ", paste(missing_vars, collapse = ", "))

# Define which PRS score is tested against which targeted phenotype list.
mode_map <- list(
  PRS_IDP1 = list(mode = "mode1", phenos = mode1_list),
  PRS_IDP2 = list(mode = "mode2", phenos = mode2_list)
)

rows <- list()
idx <- 0L
pthres_vals <- sort(unique(dat$pthres))

# Main all-threshold association loop: within each PRS threshold, compute
# covariate-adjusted partial Spearman correlations for all targeted phenotypes.
for (pt in pthres_vals) {
  sub <- dat[dat$pthres == pt, , drop = FALSE]
  d_cov <- sub[, covars, drop = FALSE]
  for (prs_name in names(mode_map)) {
    x <- as.numeric(sub[[prs_name]])
    phenos <- mode_map[[prs_name]]$phenos
    mode_name <- mode_map[[prs_name]]$mode
    for (ph in phenos) {
      y <- suppressWarnings(as.numeric(sub[[ph]]))
      st <- partial_corr_resid(x, y, d_cov, method = cfg$method, min_n = cfg$min_n)
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        mode = mode_name,
        prs_trait = prs_name,
        phenotype = ph,
        pthres = pt,
        n = st$n,
        r = st$r,
        p = st$p,
        stringsAsFactors = FALSE
      )
    }
  }
}
res_df <- do.call(rbind, rows)
res_df$display_label <- display_label(res_df$phenotype)

# FDR columns are retained for exploratory comparison; the final release also
# includes separate Bonferroni summaries created by script 02.
grp <- interaction(res_df$pthres, res_df$prs_trait, drop = TRUE)
res_df$p_fdr_within_targeted <- NA_real_
for (g in unique(grp)) {
  ii <- which(grp == g)
  res_df$p_fdr_within_targeted[ii] <- bh(res_df$p[ii])
}
res_df$p_fdr_global_by_prs <- NA_real_
for (prs_name in unique(res_df$prs_trait)) {
  ii <- which(res_df$prs_trait == prs_name)
  res_df$p_fdr_global_by_prs[ii] <- bh(res_df$p[ii])
}

ord <- order(res_df$prs_trait, res_df$phenotype, res_df$p, -abs(res_df$r), res_df$pthres)

# Select the best threshold per PRS-phenotype pair for quantile visualization.
best_df <- res_df[ord, , drop = FALSE]
best_df <- best_df[!duplicated(best_df[, c("prs_trait", "phenotype")]), , drop = FALSE]
best_df$p_fdr_best_targeted <- NA_real_
for (prs_name in unique(best_df$prs_trait)) {
  ii <- which(best_df$prs_trait == prs_name)
  best_df$p_fdr_best_targeted[ii] <- bh(best_df$p[ii])
}
best_df$phenotype_order <- list_df$phenotype_order[match(paste(best_df$mode, best_df$phenotype), paste(list_df$mode, list_df$phenotype))]
best_df <- best_df[order(best_df$mode, best_df$phenotype_order, best_df$p, -abs(best_df$r)), , drop = FALSE]

write.csv(res_df, file.path(cfg$outdir, "prs_targeted_partialcorr_all_thresholds.csv"), row.names = FALSE)
write.csv(best_df, file.path(cfg$outdir, "prs_targeted_partialcorr_best_threshold.csv"), row.names = FALSE)
write.csv(best_df[best_df$mode == "mode1", , drop = FALSE], file.path(cfg$outdir, "mode1_prs_best_threshold.csv"), row.names = FALSE)
write.csv(best_df[best_df$mode == "mode2", , drop = FALSE], file.path(cfg$outdir, "mode2_prs_best_threshold.csv"), row.names = FALSE)

plot_rows <- list()
summary_rows <- list()
plot_idx <- 0L
summary_idx <- 0L

# Build quantile plotting data from the selected threshold for each pair:
# participants are grouped by residualized PRS, and y-values summarize the
# residualized behavioral phenotype within each quantile.
for (i in seq_len(nrow(best_df))) {
  rowi <- best_df[i, , drop = FALSE]
  sub <- dat[dat$pthres == rowi$pthres, c("Subject", rowi$prs_trait, rowi$phenotype, covars), drop = FALSE]
  score_vec <- suppressWarnings(as.numeric(sub[[rowi$prs_trait]]))
  pheno_vec <- suppressWarnings(as.numeric(sub[[rowi$phenotype]]))
  d_cov <- sub[, covars, drop = FALSE]
  score_res <- residualize(score_vec, d_cov)
  pheno_res <- residualize(pheno_vec, d_cov)
  qgrp <- safe_quantile_groups(score_res, cfg$n_groups)

  tmp <- data.frame(
    Subject = sub$Subject,
    quantile_group = qgrp,
    score_res = score_res,
    pheno_res = pheno_res,
    stringsAsFactors = FALSE
  )
  tmp <- tmp[!is.na(tmp$quantile_group) & is.finite(tmp$pheno_res), , drop = FALSE]
  if (nrow(tmp) < cfg$min_n) next

  ag <- aggregate(pheno_res ~ quantile_group, tmp, function(x) c(n = length(x), mean = mean(x), se = sd(x) / sqrt(length(x))))
  ag2 <- data.frame(
    quantile_group = ag$quantile_group,
    n_group = ag$pheno_res[, "n"],
    mean = ag$pheno_res[, "mean"],
    se = ag$pheno_res[, "se"],
    stringsAsFactors = FALSE
  )
  ag2$lwr <- ag2$mean - 1.96 * ag2$se
  ag2$upr <- ag2$mean + 1.96 * ag2$se
  ag2$mode <- rowi$mode
  ag2$prs_trait <- rowi$prs_trait
  ag2$phenotype <- rowi$phenotype
  ag2$display_label <- rowi$display_label
  ag2$pthres <- rowi$pthres
  ag2$r <- rowi$r
  ag2$p <- rowi$p
  ag2$p_fdr_best_targeted <- rowi$p_fdr_best_targeted
  ag2$n_assoc <- rowi$n
  ag2$n_plot <- nrow(tmp)
  ag2$covar_note <- compact_covar_note(covars)
  ag2$pair_id <- paste(rowi$mode, rowi$phenotype, sep = "__")
  ag2$pair_label <- rowi$display_label
  plot_idx <- plot_idx + 1L
  plot_rows[[plot_idx]] <- ag2

  summary_idx <- summary_idx + 1L
  summary_rows[[summary_idx]] <- data.frame(
    mode = rowi$mode,
    prs_trait = rowi$prs_trait,
    phenotype = rowi$phenotype,
    display_label = rowi$display_label,
    pthres = rowi$pthres,
    n_assoc = rowi$n,
    r = rowi$r,
    p = rowi$p,
    p_fdr_best_targeted = rowi$p_fdr_best_targeted,
    n_plot = nrow(tmp),
    pair_id = paste(rowi$mode, rowi$phenotype, sep = "__"),
    plot_png = file.path("plots", rowi$mode, paste0(sanitize_filename(paste(rowi$mode, rowi$phenotype, sep = "__")), ".png")),
    plot_pdf = file.path("plots", rowi$mode, paste0(sanitize_filename(paste(rowi$mode, rowi$phenotype, sep = "__")), ".pdf")),
    plot_svg = file.path("plots", rowi$mode, paste0(sanitize_filename(paste(rowi$mode, rowi$phenotype, sep = "__")), ".svg")),
    stringsAsFactors = FALSE
  )
}

plot_df <- do.call(rbind, plot_rows)
summary_df <- do.call(rbind, summary_rows)
plot_df$phenotype_order <- list_df$phenotype_order[match(paste(plot_df$mode, plot_df$phenotype), paste(list_df$mode, list_df$phenotype))]
summary_df$phenotype_order <- list_df$phenotype_order[match(paste(summary_df$mode, summary_df$phenotype), paste(list_df$mode, list_df$phenotype))]
plot_df <- plot_df[order(plot_df$mode, plot_df$phenotype_order, plot_df$quantile_group), , drop = FALSE]
summary_df <- summary_df[order(summary_df$mode, summary_df$phenotype_order, summary_df$p, -abs(summary_df$r)), , drop = FALSE]

write.csv(plot_df, file.path(cfg$outdir, "best_threshold_quantile_plot_data.csv"), row.names = FALSE)
write.csv(summary_df, file.path(cfg$outdir, "best_threshold_plot_summary.csv"), row.names = FALSE)

# Export individual quantile panels for each selected PRS-phenotype pair.
for (i in seq_len(nrow(summary_df))) {
  subdf <- plot_df[plot_df$pair_id == summary_df$pair_id[i], , drop = FALSE]
  stem <- sanitize_filename(paste(summary_df$mode[i], summary_df$phenotype[i], sep = "__"))
  out_png <- file.path(cfg$outdir, "plots", summary_df$mode[i], paste0(stem, ".png"))
  out_pdf <- file.path(cfg$outdir, "plots", summary_df$mode[i], paste0(stem, ".pdf"))
  out_svg <- file.path(cfg$outdir, "plots", summary_df$mode[i], paste0(stem, ".svg"))
  plot_one_pair(subdf, out_png, out_pdf, out_svg, cfg)
}

# Export the final mode-level quantile contact sheets.
plot_contact_sheet(
  plot_df,
  "mode1",
  file.path(cfg$outdir, "mode1_best_threshold_quantile_panels.png"),
  file.path(cfg$outdir, "mode1_best_threshold_quantile_panels.pdf"),
  file.path(cfg$outdir, "mode1_best_threshold_quantile_panels.svg"),
  cfg
)
plot_contact_sheet(
  plot_df,
  "mode2",
  file.path(cfg$outdir, "mode2_best_threshold_quantile_panels.png"),
  file.path(cfg$outdir, "mode2_best_threshold_quantile_panels.pdf"),
  file.path(cfg$outdir, "mode2_best_threshold_quantile_panels.svg"),
  cfg
)

# Write a compact run summary for provenance.
writeLines(
  c(
    "HCP targeted PRS validation for response-letter use",
    paste0("Mode1 phenotypes: ", length(mode1_list)),
    paste0("Mode2 phenotypes: ", length(mode2_list)),
    paste0("All-threshold tests: ", nrow(res_df)),
    paste0("Best-threshold phenotype pairs: ", nrow(best_df)),
    paste0("Quantile plots generated: ", nrow(summary_df)),
    paste0("Covariates: ", paste(covars, collapse = ", ")),
    paste0("Quantile groups: ", cfg$n_groups),
    "Threshold-selection rule: smallest partial-correlation p-value within each PRS-trait x phenotype pair; ties broken by larger |r|."
  ),
  file.path(cfg$outdir, "analysis_summary.txt")
)

cat("Saved outputs to:", cfg$outdir, "\n")
cat("All-threshold tests:", nrow(res_df), "\n")
cat("Best-threshold pairs:", nrow(best_df), "\n")
cat("Plots generated:", nrow(summary_df), "\n")
