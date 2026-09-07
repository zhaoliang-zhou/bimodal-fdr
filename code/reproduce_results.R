#!/usr/bin/env Rscript

# Reproducibility script for the empirical-Bayes auxiliary-information analyses.
#
# Usage from this directory:
#   Rscript reproduce_results.R
#   Rscript reproduce_results.R --data=/path/to/LLD_data.RData --output=results
#
# The input RData file must contain fMRI.est.all and DTI.est.all. The script
# writes all generated figures, estimates, and numerical checks to the output
# directory. It does not change the input data.

required_packages <- c("ggplot2", "mgcv")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

get_script_dir <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  script_argument <- command_line[grepl("^--file=", command_line)]
  if (length(script_argument) > 0L) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_argument[[1L]]))
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  source_files <- vapply(
    sys.frames(),
    function(frame) if (!is.null(frame$ofile)) frame$ofile else "",
    character(1)
  )
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files) > 0L) {
    return(dirname(normalizePath(tail(source_files, 1L), mustWork = TRUE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_document <- rstudioapi::getActiveDocumentContext()
    if (!is.null(active_document$path) && nzchar(active_document$path)) {
      return(dirname(normalizePath(active_document$path, mustWork = TRUE)))
    }
  }

  getwd()
}

script_dir <- get_script_dir()

arguments <- commandArgs(trailingOnly = TRUE)
argument_value <- function(name) {
  prefix <- paste0("--", name, "=")
  hit <- arguments[startsWith(arguments, prefix)]
  if (length(hit) == 0L) return(NULL)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

data_file <- argument_value("data")
if (is.null(data_file)) {
  data_candidates <- c(
    file.path(script_dir, "LLD_data.RData"),
    file.path(script_dir, "data", "LLD_data.RData")
  )
  existing_candidates <- data_candidates[file.exists(data_candidates)]
  if (length(existing_candidates) == 0L) {
    stop(
      "No input RData file was found. Supply one with --data=/path/to/LLD_data.RData.\n",
      "Searched:\n  ", paste(data_candidates, collapse = "\n  ")
    )
  }
  data_file <- existing_candidates[[1L]]
}
data_file <- normalizePath(data_file, mustWork = TRUE)

output_dir <- argument_value("output")
if (is.null(output_dir)) {
  output_dir <- file.path(script_dir, "results")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
figures_dir <- file.path(output_dir, "Figures")
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

source(file.path(script_dir, "utilities.R"), local = TRUE)

analysis_environment <- new.env(parent = emptyenv())
loaded_objects <- load(data_file, envir = analysis_environment)
required_objects <- c("fMRI.est.all", "DTI.est.all")
if (!all(required_objects %in% loaded_objects)) {
  stop(
    "The input RData file must contain: ",
    paste(required_objects, collapse = " and "),
    ". Found: ",
    paste(loaded_objects, collapse = ", ")
  )
}
fmri.df <- get("fMRI.est.all", envir = analysis_environment)
dti.df <- get("DTI.est.all", envir = analysis_environment)

analytical_dat <- create_analytical_dataset(
  fmri.df = fmri.df,
  dti.df = dti.df,
  p_cutoff = 0.05
)
analytical_dat <- analytical_dat[
  is.finite(analytical_dat$fmri_t_value) &
    is.finite(analytical_dat$dti_t_value) &
    analytical_dat$fmri_t_value > 0 &
    analytical_dat$dti_t_value > 0,
  ,
  drop = FALSE
]
analytical_dat$is_non_null <- as.integer(analytical_dat$group == "non-null")

# Estimate the Gamma models by matching the first two moments, as in the
# analysis document. The shifted location is the smallest positive non-null
# absolute FC statistic.
x1 <- analytical_dat$fmri_t_value[analytical_dat$group == "non-null"]
if (length(x1) < 3L) {
  stop("At least three non-null FC statistics are required for estimation.")
}
mu1 <- min(x1)
est <- estimate_mom_from_analytical_dat(
  analytical_dat = analytical_dat,
  mu1 = mu1,
  dti_group = "non-null"
)

alpha0_mom <- est$alpha0_mom
beta0_mom <- est$beta0_mom
alpha1_mom <- est$alpha1_mom
beta1_mom <- est$beta1_mom
alpha_s_mom <- est$alpha_s_mom
beta_s_mom <- est$beta_s_mom
alpha_bimodal_mom <- est$alpha_bimodal_mom
pi0 <- 0.93
pi0_t2 <- 0.9083
alpha_level <- 0.20

estimates <- data.frame(
  parameter = c("pi0", "pi0_t2", "mu1", "alpha0", "beta0", "alpha1", "beta1", "alpha_s", "beta_s", "alpha_bimodal"),
  estimate = c(pi0, pi0_t2, mu1, alpha0_mom, beta0_mom, alpha1_mom, beta1_mom, alpha_s_mom, beta_s_mom, alpha_bimodal_mom)
)
utils::write.csv(estimates, file.path(output_dir, "estimates.csv"), row.names = FALSE)
utils::write.csv(analytical_dat, file.path(output_dir, "analytical_data.csv"), row.names = FALSE)

save_plot <- function(plot_object, filename, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

# Likelihood-ratio curve.
t_grid <- seq(
  min(analytical_dat$fmri_t_value),
  max(analytical_dat$fmri_t_value),
  length.out = 2000
)
lr_data <- data.frame(
  t = t_grid,
  likelihood_ratio = lr_gamma(
    t_grid,
    alpha0 = alpha0_mom,
    beta0 = beta0_mom,
    mu1 = mu1,
    alpha_alt = alpha1_mom,
    beta_alt = beta1_mom
  )
)
lr_plot <- ggplot2::ggplot(lr_data, ggplot2::aes(t, likelihood_ratio)) +
  ggplot2::geom_line(linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = mu1, linetype = "dashed", color = "gray35") +
  ggplot2::labs(
    title = "Likelihood-ratio function",
    x = "Absolute FC statistic",
    y = expression(f[1](t) / f[0](t))
  ) +
  ggplot2::theme_bw(base_size = 13)
save_plot(lr_plot, "likelihood_ratio.png")

# LFDR step-up thresholds under the FC-only and enriched alternative models.
c_star_fc <- find_c_star_from_t1(
  t1 = analytical_dat$fmri_t_value,
  alpha0 = alpha0_mom,
  beta0 = beta0_mom,
  mu1 = mu1,
  alpha_alt = alpha1_mom,
  beta_alt = beta1_mom,
  pi0 = pi0,
  alpha = alpha_level
)
c_star_fc_sc <- find_c_star_from_t1(
  t1 = analytical_dat$fmri_t_value,
  alpha0 = alpha0_mom,
  beta0 = beta0_mom,
  mu1 = mu1,
  alpha_alt = alpha_bimodal_mom,
  beta_alt = beta1_mom,
  pi0 = pi0,
  alpha = alpha_level
)
thresholds <- data.frame(
  model = c("FC-only", "FC-plus-SC"),
  threshold = c(c_star_fc$c_star, c_star_fc_sc$c_star),
  number_rejected = c(c_star_fc$n_reject, c_star_fc_sc$n_reject),
  average_lfdr = c(c_star_fc$avg_lfdr, c_star_fc_sc$avg_lfdr)
)
utils::write.csv(thresholds, file.path(output_dir, "thresholds.csv"), row.names = FALSE)

# Enrichment power comparison using the same FC-only rejection region.
enrichment_plot <- plot_enrichment_power_inequality(
  c_values = seq(0.01, 0.99, by = 0.01),
  alpha0 = alpha0_mom,
  beta0 = beta0_mom,
  mu1 = mu1,
  alpha1 = alpha1_mom,
  beta1 = beta1_mom,
  alpha_s = alpha_s_mom,
  beta_s = beta_s_mom,
  pi0 = pi0,
  alpha_level = alpha_level
)
save_plot(enrichment_plot$power_plot, "enrichment_power.png")
save_plot(enrichment_plot$difference_plot, "enrichment_power_difference.png")

# Stochastic-dominance density and survival curves.
t_max <- max(
  analytical_dat$fmri_t_value,
  mu1 + qgamma(0.999, shape = alpha_bimodal_mom, rate = beta1_mom)
)
dominance_grid <- seq(mu1, t_max, length.out = 2000)
dominance_data <- rbind(
  data.frame(
    t = dominance_grid,
    value = dgamma(dominance_grid - mu1, shape = alpha1_mom, rate = beta1_mom),
    model = "FC"
  ),
  data.frame(
    t = dominance_grid,
    value = dgamma(dominance_grid - mu1, shape = alpha_bimodal_mom, rate = beta1_mom),
    model = "FC-plus-SC"
  )
)
dominance_density_plot <- ggplot2::ggplot(
  dominance_data,
  ggplot2::aes(t, value, linetype = model)
) +
  ggplot2::geom_line(color = "black", linewidth = 0.5) +
  ggplot2::labs(title = "Non-null density comparison", x = "Primary statistic", y = "Density", linetype = NULL) +
  ggplot2::theme_bw(base_size = 13)
save_plot(dominance_density_plot, "stochastic_dominance_density.png")

survival_data <- rbind(
  data.frame(
    t = dominance_grid,
    value = pgamma(dominance_grid - mu1, shape = alpha1_mom, rate = beta1_mom, lower.tail = FALSE),
    model = "FC"
  ),
  data.frame(
    t = dominance_grid,
    value = pgamma(dominance_grid - mu1, shape = alpha_bimodal_mom, rate = beta1_mom, lower.tail = FALSE),
    model = "FC-plus-SC"
  )
)
survival_plot <- ggplot2::ggplot(
  survival_data,
  ggplot2::aes(t, value, linetype = model)
) +
  ggplot2::geom_line(color = "black", linewidth = 0.5) +
  ggplot2::labs(title = "Non-null survival-function comparison", x = "Primary statistic", y = "Survival probability", linetype = NULL) +
  ggplot2::theme_minimal(base_size = 12)
save_plot(survival_plot, "stochastic_dominance_survival.png")

dominance_check <- data.frame(
  survival_ordering_holds = all(
    survival_data$value[survival_data$model == "FC-plus-SC"] >=
      survival_data$value[survival_data$model == "FC"] - 1e-12
  )
)
utils::write.csv(dominance_check, file.path(output_dir, "stochastic_dominance_check.csv"), row.names = FALSE)

# Conditional Storey estimate using a smooth model for P(p > lambda | SC).
lambda <- 0.5
storey_dat <- data.frame(
  fmri_p_value = fmri.df$p.value,
  dti_t_value = abs(dti.df$t)
)
storey_dat <- storey_dat[
  is.finite(storey_dat$fmri_p_value) &
    is.finite(storey_dat$dti_t_value) &
    storey_dat$fmri_p_value >= 0 & storey_dat$fmri_p_value <= 1 &
    storey_dat$dti_t_value > 0,
  ,
  drop = FALSE
]
storey_dat$p_above_lambda <- as.integer(storey_dat$fmri_p_value > lambda)
pi0_storey_marginal <- min(max(mean(storey_dat$p_above_lambda) / (1 - lambda), 0), 1)
conditional_storey_fit <- mgcv::gam(
  p_above_lambda ~ s(dti_t_value, k = 6, bs = "cr"),
  data = storey_dat,
  family = stats::binomial(link = "logit"),
  method = "REML"
)
storey_dat$pi0_t2_storey <- pmin(
  pmax(
    stats::predict(conditional_storey_fit, newdata = storey_dat, type = "response") / (1 - lambda),
    0
  ),
  1
)
pi0_sc_storey <- mean(storey_dat$pi0_t2_storey)
storey_grid <- data.frame(
  dti_t_value = seq(min(storey_dat$dti_t_value), max(storey_dat$dti_t_value), length.out = 500)
)
storey_link <- stats::predict(conditional_storey_fit, newdata = storey_grid, type = "link", se.fit = TRUE)
storey_grid$pi0 <- pmin(pmax(stats::plogis(storey_link$fit) / (1 - lambda), 0), 1)
storey_grid$lower <- pmin(pmax(stats::plogis(storey_link$fit - 1.96 * storey_link$se.fit) / (1 - lambda), 0), 1)
storey_grid$upper <- pmin(pmax(stats::plogis(storey_link$fit + 1.96 * storey_link$se.fit) / (1 - lambda), 0), 1)
storey_plot <- ggplot2::ggplot(storey_grid, ggplot2::aes(dti_t_value, pi0)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "gray80") +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_hline(yintercept = pi0, linetype = "dashed", color = "gray35") +
  ggplot2::labs(title = "Conditional Storey null-proportion estimate", x = "Absolute SC statistic", y = "Estimated null proportion") +
  ggplot2::theme_classic(base_size = 13)
save_plot(storey_plot, "conditional_storey_pi0.png")
utils::write.csv(
  data.frame(marginal_storey = pi0_storey_marginal, sc_informed_storey = pi0_sc_storey),
  file.path(output_dir, "storey_estimates.csv"),
  row.names = FALSE
)

# Figure 1: reproduce the four Gaussian power curves from
# theorem3_power_curves.R. These curves use the numerical illustration's
# N(0, 1) null, N(2, 1) baseline alternative, and N(2.4, 1) enriched
# alternative, with prior null proportions 0.90 and 0.88.
figure1_pi0 <- 0.90
figure1_pi0_t2 <- 0.88
figure1_mu0 <- 2.0
figure1_mu_enriched <- 2.4
figure1_c_star <- 0.3917
figure1_c_grid <- seq(0.05, 0.50, length.out = 1000)

gaussian_power <- function(c, pi, mu) {
  likelihood_ratio_threshold <- (pi / (1 - pi)) * ((1 - c) / c)
  rejection_boundary <-
    (log(likelihood_ratio_threshold) + mu^2 / 2) / mu
  stats::pnorm(rejection_boundary - mu, lower.tail = FALSE)
}

figure1_power <- data.frame(
  c = rep(figure1_c_grid, times = 4),
  power = c(
    gaussian_power(figure1_c_grid, figure1_pi0_t2, figure1_mu_enriched),
    gaussian_power(figure1_c_grid, figure1_pi0, figure1_mu_enriched),
    gaussian_power(figure1_c_grid, figure1_pi0_t2, figure1_mu0),
    gaussian_power(figure1_c_grid, figure1_pi0, figure1_mu0)
  ),
  curve = factor(
    rep(
      c("joint", "nonnull_enriched", "prior_refined", "baseline"),
      each = length(figure1_c_grid)
    ),
    levels = c("joint", "nonnull_enriched", "prior_refined", "baseline")
  )
)

figure1_cutoff <- function(c, pi, mu) {
  likelihood_ratio_threshold <- (pi / (1 - pi)) * ((1 - c) / c)
  (log(likelihood_ratio_threshold) + mu^2 / 2) / mu
}

figure1_power_at_star <- c(
  gaussian_power(figure1_c_star, figure1_pi0_t2, figure1_mu_enriched),
  gaussian_power(figure1_c_star, figure1_pi0, figure1_mu_enriched),
  gaussian_power(figure1_c_star, figure1_pi0_t2, figure1_mu0),
  gaussian_power(figure1_c_star, figure1_pi0, figure1_mu0)
)
figure1_cutoffs <- c(
  joint = figure1_cutoff(figure1_c_star, figure1_pi0_t2, figure1_mu_enriched),
  nonnull_enriched = figure1_cutoff(figure1_c_star, figure1_pi0, figure1_mu_enriched),
  prior_refined = figure1_cutoff(figure1_c_star, figure1_pi0_t2, figure1_mu0),
  baseline = figure1_cutoff(figure1_c_star, figure1_pi0, figure1_mu0)
)
figure1_star <- data.frame(
  c = rep(figure1_c_star, 4),
  power = figure1_power_at_star,
  label = sprintf("%.3f", figure1_power_at_star),
  label_x = figure1_c_star + 0.011,
  label_y = figure1_power_at_star + c(0.05, -0.006, 0.05, -0.006)
)

figure1_legend <- as.expression(list(
  bquote("Top (joint): " * f[1](t[1], t[2]) * ", " * pi[0](t[2]) * "  (" * C[t[2]] == "[" * .(sprintf("%.3f", figure1_cutoffs[["joint"]])) * ", " * infinity * ")"),
  bquote("Upper-middle (non-null enriched): " * f[1](t[1], t[2]) * ", " * pi[0] * "  (" * A[t[2]] == "[" * .(sprintf("%.3f", figure1_cutoffs[["nonnull_enriched"]])) * ", " * infinity * ")"),
  bquote("Lower-middle (prior refined): " * f[1](t[1]) * ", " * pi[0](t[2]) * "  (" * B[0] == "[" * .(sprintf("%.3f", figure1_cutoffs[["prior_refined"]])) * ", " * infinity * ")"),
  bquote("Bottom (baseline): " * f[1](t[1]) * ", " * pi[0] * "  (" * A[0] == "[" * .(sprintf("%.3f", figure1_cutoffs[["baseline"]])) * ", " * infinity * ")")
))

figure1_plot <- ggplot2::ggplot(
  figure1_power,
  ggplot2::aes(x = c, y = power, linetype = curve)
) +
  ggplot2::geom_line(color = "black", linewidth = 0.5) +
  ggplot2::scale_linetype_manual(
    values = c(joint = "solid", nonnull_enriched = "dotdash", prior_refined = "longdash", baseline = "dotted"),
    breaks = levels(figure1_power$curve),
    labels = figure1_legend
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0.05, 0.50, by = 0.05),
    labels = function(x) sprintf("%.2f", x)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 0.6, by = 0.2),
    labels = function(y) sprintf("%.1f", y)
  ) +
  ggplot2::coord_cartesian(xlim = c(0.05, 0.50), ylim = c(0, 0.65), clip = "off") +
  ggplot2::labs(
    title = "Power curves for the four procedures",
    x = "LFDR threshold c",
    y = "Power",
    linetype = NULL
  ) +
  ggplot2::guides(
    linetype = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(color = "black", linewidth = 0.5, linetype = c("solid", "dotdash", "longdash", "dotted"))
    )
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 6.5),
    legend.key.width = grid::unit(1.5, "cm"),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::geom_vline(xintercept = figure1_c_star, color = "gray35", linetype = "dashed", linewidth = 0.55) +
  ggplot2::geom_point(
    data = figure1_star,
    ggplot2::aes(x = c, y = power),
    color = "black",
    size = 1.4,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = figure1_star,
    ggplot2::aes(x = label_x, y = label_y, label = label),
    color = "gray35",
    size = 3,
    hjust = 0,
    inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = figure1_c_star + 0.03,
    y = min(figure1_power$power) + 0.02,
    label = paste0("c*=", formatC(figure1_c_star, format = "f", digits = 3)),
    size = 3,
    color = "gray25",
    vjust = -0.6
  )

save_plot(figure1_plot, "Theorem3_Power_Curves.png", dpi = 600)
utils::write.csv(
  data.frame(
    procedure = c("joint", "nonnull_enriched", "prior_refined", "baseline"),
    power_at_c_star = figure1_power_at_star,
    rejection_cutoff = unname(figure1_cutoffs)
  ),
  file.path(output_dir, "figure1_power_values.csv"),
  row.names = FALSE
)

# Figure 2(a): generate the TP2 power-ratio comparison across prior-proportion
# levels. The resulting file is tp2_prior_ratio.png.
tp2_prior_plot <- plot_tp2_prior_ratio_inequality(
  c_values = seq(0.01, 0.99, by = 0.01),
  alpha0 = alpha0_mom,
  beta0 = beta0_mom,
  mu1 = mu1,
  alpha1 = alpha1_mom,
  beta1 = beta1_mom,
  alpha_s = alpha_s_mom,
  beta_s = beta_s_mom,
  pi02 = 0.98,
  pi01 = 0.95,
  t21 = 0,
  t22 = 2,
  c_star = 0.35,
  x_limits = c(0.15, 0.75),
  print_c_star = TRUE
)
save_plot(tp2_prior_plot$ratio_plot, "tp2_prior_ratio.png", dpi = 600)

# Figure 2(b): generate the TP2 power-ratio comparison across auxiliary-
# information levels. The resulting file is tp2_auxiliary_ratio.png.
tp2_auxiliary_plot <- plot_tp2_auxiliary_ratio_inequality(
  c_values = seq(0.01, 0.99, by = 0.01),
  alpha0 = alpha0_mom,
  beta0 = beta0_mom,
  mu1 = mu1,
  alpha1 = alpha1_mom,
  beta1 = beta1_mom,
  pi02 = 0.98,
  pi01 = 0.95,
  alpha_s = alpha_s_mom,
  beta_s = beta_s_mom,
  t21 = 0,
  t22 = 2,
  c_star = 0.35,
  x_limits = c(0, 0.75),
  print_c_star = TRUE
)
save_plot(tp2_auxiliary_plot$ratio_plot, "tp2_auxiliary_ratio.png", dpi = 600)

summary_lines <- c(
  paste0("Input data: ", data_file),
  paste0("Rows after filtering: ", nrow(analytical_dat)),
  paste0("Null rows: ", sum(analytical_dat$group == "null")),
  paste0("Non-null rows: ", sum(analytical_dat$group == "non-null")),
  paste0("FC-only threshold: ", signif(c_star_fc$c_star, 8)),
  paste0("FC-only rejections: ", c_star_fc$n_reject),
  paste0("FC-plus-SC threshold: ", signif(c_star_fc_sc$c_star, 8)),
  paste0("FC-plus-SC rejections: ", c_star_fc_sc$n_reject),
  paste0("Marginal Storey estimate: ", signif(pi0_storey_marginal, 8)),
  paste0("SC-informed Storey estimate: ", signif(pi0_sc_storey, 8)),
  paste0("Stochastic-ordering check: ", dominance_check$survival_ordering_holds),
  paste0("TP2 prior-ratio ordering check: ", tp2_prior_plot$ordering_holds),
  paste0("TP2 auxiliary-ratio ordering check: ", tp2_auxiliary_plot$ordering_holds)
)
writeLines(summary_lines, con = file.path(output_dir, "results_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
