# Power helpers for the empirical-Bayes auxiliary-information analyses.
#
# Source this file before calling the analysis or plotting functions below.
# The default model uses:
#   f0(t1)       = Gamma(alpha0, beta0)
#   f1(t1)       = shifted Gamma(mu1, alpha1, beta1)
#   f1(t1, t2)   = shifted Gamma(mu1, alpha1 + alpha_s / beta_s, beta1)
#   k(c, pi0)    = (pi0 / (1 - pi0)) * ((1 - c) / c)

pi0_default <- 0.93
pi0_t2_default <- 0.9083

create_analytical_dataset <- function(fmri.df, dti.df, p_cutoff = 0.05) {
  required_cols <- c("t", "p.value")
  missing_fmri <- setdiff(required_cols, names(fmri.df))
  missing_dti <- setdiff("t", names(dti.df))

  if (length(missing_fmri) > 0) {
    stop("fmri.df is missing required column(s): ", paste(missing_fmri, collapse = ", "))
  }
  if (length(missing_dti) > 0) {
    stop("dti.df is missing required column(s): ", paste(missing_dti, collapse = ", "))
  }
  if (length(fmri.df$t) != length(dti.df$t)) {
    stop("fmri.df$t and dti.df$t must have the same length.")
  }

  analytical_dat <- data.frame(
    group = ifelse(fmri.df$p.value <= p_cutoff, "non-null", "null"),
    fmri_t_value = abs(fmri.df$t),
    dti_t_value = abs(dti.df$t)
  )

  analytical_dat$group <- factor(
    analytical_dat$group,
    levels = c("null", "non-null")
  )

  analytical_dat
}

estimate_mom_from_analytical_dat <- function(analytical_dat,
                                             mu1 = 2,
                                             dti_group = "non-null",
                                             tol = 1e-8) {
  required_cols <- c("group", "fmri_t_value", "dti_t_value")
  missing_cols <- setdiff(required_cols, names(analytical_dat))

  if (length(missing_cols) > 0) {
    stop("analytical_dat is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  clean_dat <- analytical_dat[
    is.finite(analytical_dat$fmri_t_value) &
      is.finite(analytical_dat$dti_t_value) &
      analytical_dat$fmri_t_value > 0 &
      analytical_dat$dti_t_value > 0,
    ,
    drop = FALSE
  ]

  clean_dat$group <- factor(clean_dat$group, levels = c("null", "non-null"))

  x0 <- clean_dat$fmri_t_value[clean_dat$group == "null"]
  x1 <- clean_dat$fmri_t_value[clean_dat$group == "non-null"]

  y1 <- x1 - mu1
  y1 <- y1[is.finite(y1) & y1 > 0]

  t_s <- clean_dat$dti_t_value[clean_dat$group == dti_group]
  t_s <- t_s[is.finite(t_s) & t_s > tol]

  if (length(x0) < 2) {
    stop("Need at least two null fMRI t-values to estimate alpha0 and beta0.")
  }
  if (length(y1) < 2) {
    stop("Need at least two shifted non-null fMRI t-values greater than mu1.")
  }
  if (length(t_s) < 2) {
    stop("Need at least two DTI t-values in dti_group to estimate alpha_s and beta_s.")
  }

  m0 <- mean(x0)
  v0 <- var(x0)
  m1 <- mean(y1)
  v1 <- var(y1)
  mean_s <- mean(t_s)
  var_s <- var(t_s)

  alpha0_mom <- m0^2 / v0
  beta0_mom <- m0 / v0
  alpha1_mom <- m1^2 / v1
  beta1_mom <- m1 / v1
  alpha_s_mom <- mean_s^2 / var_s
  beta_s_mom <- mean_s / var_s
  alpha_bimodal_mom <- alpha1_mom + alpha_s_mom / beta_s_mom

  list(
    analytical_dat = clean_dat,
    mu1 = mu1,
    alpha0_mom = alpha0_mom,
    beta0_mom = beta0_mom,
    alpha1_mom = alpha1_mom,
    beta1_mom = beta1_mom,
    alpha_s_mom = alpha_s_mom,
    beta_s_mom = beta_s_mom,
    alpha_bimodal_mom = alpha_bimodal_mom,
    n_null = length(x0),
    n_nonnull = length(x1),
    n_shifted_nonnull = length(y1),
    n_dti = length(t_s)
  )
}

k_fun <- function(c, pi0 = pi0_default) {
  if (any(c <= 0 | c >= 1, na.rm = TRUE)) {
    stop("c must be in (0, 1).")
  }
  if (pi0 <= 0 || pi0 >= 1) {
    stop("pi0 must be in (0, 1).")
  }

  (pi0 / (1 - pi0)) * ((1 - c) / c)
}

c_from_k <- function(k, pi0 = pi0_default) {
  if (any(k < 0, na.rm = TRUE)) {
    stop("k must be non-negative.")
  }
  if (pi0 <= 0 || pi0 >= 1) {
    stop("pi0 must be in (0, 1).")
  }

  a <- pi0 / (1 - pi0)
  a / (a + k)
}

f0_gamma <- function(t1, alpha0, beta0) {
  ifelse(t1 <= 0, 0, dgamma(t1, shape = alpha0, rate = beta0))
}

f1_shifted_gamma <- function(t1, mu1, alpha1, beta1) {
  ifelse(t1 <= mu1, 0, dgamma(t1 - mu1, shape = alpha1, rate = beta1))
}

f1_bimodal_gamma <- function(t1, mu1, alpha1, beta1, alpha_s, beta_s) {
  alpha_bimodal <- alpha1 + alpha_s / beta_s
  f1_shifted_gamma(t1, mu1 = mu1, alpha1 = alpha_bimodal, beta1 = beta1)
}

lr_gamma <- function(t1, alpha0, beta0, mu1, alpha_alt, beta_alt) {
  f0 <- f0_gamma(t1, alpha0 = alpha0, beta0 = beta0)
  f1 <- f1_shifted_gamma(t1, mu1 = mu1, alpha1 = alpha_alt, beta1 = beta_alt)

  ifelse(f0 <= 0, Inf, f1 / f0)
}

lfdr_from_lr <- function(lr, pi0 = pi0_default) {
  1 / (1 + ((1 - pi0) / pi0) * lr)
}

find_t_star_gamma <- function(c,
                              alpha0,
                              beta0,
                              mu1,
                              alpha_alt,
                              beta_alt,
                              pi0 = pi0_default,
                              x_max = NULL,
                              n_grid = 20000) {
  k <- k_fun(c, pi0 = pi0)

  if (is.null(x_max)) {
    x_max <- max(
      qgamma(1 - 1e-8, shape = alpha0, rate = beta0),
      mu1 + qgamma(1 - 1e-8, shape = alpha_alt, rate = beta_alt)
    )
  }

  x_min <- max(.Machine$double.eps, mu1 + 1e-8)

  for (attempt in seq_len(8)) {
    x_grid <- seq(x_min, x_max, length.out = n_grid)
    lr_vals <- lr_gamma(
      x_grid,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha_alt,
      beta_alt = beta_alt
    )

    idx <- which(lr_vals >= k)[1]
    if (!is.na(idx)) {
      if (idx == 1) {
        return(x_grid[idx])
      }

      root_fun <- function(x) {
        lr_gamma(
          x,
          alpha0 = alpha0,
          beta0 = beta0,
          mu1 = mu1,
          alpha_alt = alpha_alt,
          beta_alt = beta_alt
        ) - k
      }

      root <- tryCatch(
        uniroot(
          root_fun,
          lower = x_grid[idx - 1],
          upper = x_grid[idx],
          tol = sqrt(.Machine$double.eps)
        )$root,
        error = function(e) x_grid[idx]
      )
      return(root)
    }

    x_max <- x_max * 2
  }

  Inf
}

reject_prob_null <- function(t_star, alpha0, beta0) {
  if (is.infinite(t_star)) {
    return(0)
  }
  if (t_star <= 0) {
    return(1)
  }

  pgamma(t_star, shape = alpha0, rate = beta0, lower.tail = FALSE)
}

power_shifted_gamma <- function(t_star, mu1, alpha_alt, beta_alt) {
  if (is.infinite(t_star)) {
    return(0)
  }
  if (t_star <= mu1) {
    return(1)
  }

  pgamma(t_star - mu1, shape = alpha_alt, rate = beta_alt, lower.tail = FALSE)
}

compute_power_three <- function(c_values,
                                alpha0,
                                beta0,
                                mu1,
                                alpha1,
                                beta1,
                                alpha_s,
                                beta_s,
                                pi0 = pi0_default) {
  alpha_bimodal <- alpha1 + alpha_s / beta_s

  rows <- lapply(c_values, function(c_val) {
    t_star_fc <- find_t_star_gamma(
      c = c_val,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1,
      pi0 = pi0
    )

    t_star_bimodal <- find_t_star_gamma(
      c = c_val,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha_bimodal,
      beta_alt = beta1,
      pi0 = pi0
    )

    data.frame(
      c = c_val,
      k = k_fun(c_val, pi0 = pi0),
      model = c(
        "f0_null_rejection_under_fc_rule",
        "f1_fc_power",
        "f1_fc_sc_power"
      ),
      alpha_alt = c(NA_real_, alpha1, alpha_bimodal),
      beta_alt = c(NA_real_, beta1, beta1),
      t_star = c(t_star_fc, t_star_fc, t_star_bimodal),
      probability = c(
        reject_prob_null(t_star_fc, alpha0 = alpha0, beta0 = beta0),
        power_shifted_gamma(t_star_fc, mu1 = mu1, alpha_alt = alpha1, beta_alt = beta1),
        power_shifted_gamma(
          t_star_bimodal,
          mu1 = mu1,
          alpha_alt = alpha_bimodal,
          beta_alt = beta1
        )
      )
    )
  })

  do.call(rbind, rows)
}

find_c_star_from_lfdr <- function(lfdr, alpha = 0.05, pi0 = pi0_default) {
  lfdr <- lfdr[is.finite(lfdr)]
  if (length(lfdr) == 0) {
    return(list(
      c_star = 0,
      k_star = 0,
      k_star_c_pi0 = Inf,
      n_reject = 0,
      avg_lfdr = NA_real_
    ))
  }

  lfdr_ordered <- sort(lfdr)
  avg_lfdr <- cumsum(lfdr_ordered) / seq_along(lfdr_ordered)
  ok <- which(avg_lfdr <= alpha)

  if (length(ok) == 0) {
    return(list(
      c_star = 0,
      k_star = 0,
      k_star_c_pi0 = Inf,
      n_reject = 0,
      avg_lfdr = NA_real_
    ))
  }

  k_star_index <- max(ok)
  c_star <- lfdr_ordered[k_star_index]

  list(
    c_star = c_star,
    k_star = k_star_index,
    k_star_c_pi0 = k_fun(c_star, pi0 = pi0),
    n_reject = k_star_index,
    avg_lfdr = avg_lfdr[k_star_index]
  )
}

find_c_star_from_t1 <- function(t1,
                                alpha0,
                                beta0,
                                mu1,
                                alpha_alt,
                                beta_alt,
                                pi0 = pi0_default,
                                alpha = 0.05) {
  lr <- lr_gamma(
    t1,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha_alt = alpha_alt,
    beta_alt = beta_alt
  )
  lfdr <- lfdr_from_lr(lr, pi0 = pi0)
  out <- find_c_star_from_lfdr(lfdr, alpha = alpha, pi0 = pi0)
  out$k_star_c_pi0 <- if (out$c_star > 0) k_fun(out$c_star, pi0 = pi0) else Inf
  out$pi0 <- pi0
  out$alpha <- alpha
  out
}

plot_power_curves_three <- function(c_values = seq(0.01, 0.99, by = 0.01),
                                    alpha0,
                                    beta0,
                                    mu1,
                                    alpha1,
                                    beta1,
                                    alpha_s,
                                    beta_s,
                                    pi0 = pi0_default,
                                    alpha_level = 0.05,
                                    t1 = NULL,
                                    mark_c_star = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_power_curves_three().")
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("alpha_level must be in (0, 1).")
  }
  if (any(c_values <= 0 | c_values >= 1, na.rm = TRUE)) {
    stop("All c_values must be in (0, 1).")
  }

  alpha_bimodal <- alpha1 + alpha_s / beta_s

  power_df <- compute_power_three(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi0 = pi0
  )

  power_df$curve <- factor(
    power_df$model,
    levels = c(
      "f0_null_rejection_under_fc_rule",
      "f1_fc_power",
      "f1_fc_sc_power"
    ),
    labels = c(
      "Null: f0(t1)",
      "FC only: f1(t1)",
      "FC + SC: f1(t1, t2)"
    )
  )

  p <- ggplot2::ggplot(
    power_df,
    ggplot2::aes(x = c, y = probability, color = curve, linetype = curve)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Power Curves Across LFDR Thresholds",
      subtitle = sprintf("pi0 = %.2f; FDR alpha level = %.2f", pi0, alpha_level),
      x = "LFDR threshold c",
      y = "Probability",
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  c_star_df <- data.frame()

  if (mark_c_star && !is.null(t1)) {
    c_star_fc <- find_c_star_from_t1(
      t1 = t1,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1,
      pi0 = pi0,
      alpha = alpha_level
    )

    c_star_fc_sc <- find_c_star_from_t1(
      t1 = t1,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha_bimodal,
      beta_alt = beta1,
      pi0 = pi0,
      alpha = alpha_level
    )

    c_star_df <- data.frame(
      model = c("FC only", "FC + SC"),
      c_star = c(c_star_fc$c_star, c_star_fc_sc$c_star),
      k_star = c(c_star_fc$k_star, c_star_fc_sc$k_star),
      k_star_c_pi0 = c(c_star_fc$k_star_c_pi0, c_star_fc_sc$k_star_c_pi0),
      n_reject = c(c_star_fc$n_reject, c_star_fc_sc$n_reject),
      avg_lfdr = c(c_star_fc$avg_lfdr, c_star_fc_sc$avg_lfdr)
    )

    c_star_plot_df <- c_star_df[
      is.finite(c_star_df$c_star) & c_star_df$c_star > 0,
      ,
      drop = FALSE
    ]

    if (nrow(c_star_plot_df) > 0) {
      p <- p +
        ggplot2::geom_vline(
          data = c_star_plot_df,
          ggplot2::aes(xintercept = c_star),
          color = "gray35",
          linetype = "dashed",
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = c_star_plot_df,
          ggplot2::aes(
            x = c_star,
            y = 0.95,
            label = paste0(model, " c* = ", round(c_star, 3))
          ),
          angle = 90,
          vjust = -0.4,
          hjust = 1,
          size = 3,
          color = "gray25",
          inherit.aes = FALSE
        )
    }
  }

  list(
    plot = p,
    power_data = power_df,
    c_star_data = c_star_df
  )
}

compute_enrichment_power_inequality <- function(c_values,
                                        alpha0,
                                        beta0,
                                        mu1,
                                        alpha1,
                                        beta1,
                                        alpha_s,
                                        beta_s,
                                        pi0 = pi0_default) {
  alpha_bimodal <- alpha1 + alpha_s / beta_s

  rows <- lapply(c_values, function(c_val) {
    t_star_fc <- find_t_star_gamma(
      c = c_val,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1,
      pi0 = pi0
    )

    power_fc <- power_shifted_gamma(
      t_star_fc,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1
    )

    power_fc_sc_same_region <- power_shifted_gamma(
      t_star_fc,
      mu1 = mu1,
      alpha_alt = alpha_bimodal,
      beta_alt = beta1
    )

    data.frame(
      c = c_val,
      k = k_fun(c_val, pi0 = pi0),
      t_star = t_star_fc,
      power_fc = power_fc,
      power_fc_sc_same_region = power_fc_sc_same_region,
      power_difference = power_fc_sc_same_region - power_fc,
      ordering_holds = power_fc_sc_same_region >= power_fc - 1e-12
    )
  })

  do.call(rbind, rows)
}

plot_enrichment_power_inequality <- function(c_values = seq(0.01, 0.99, by = 0.01),
                                     alpha0,
                                     beta0,
                                     mu1,
                                     alpha1,
                                     beta1,
                                     alpha_s,
                                     beta_s,
                                     pi0 = pi0_default,
                                     alpha_level = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_enrichment_power_inequality().")
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("alpha_level must be in (0, 1).")
  }

  enrichment_df <- compute_enrichment_power_inequality(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi0 = pi0
  )

  power_long <- rbind(
    data.frame(
      c = enrichment_df$c,
      curve = "FC",
      power = enrichment_df$power_fc
    ),
    data.frame(
      c = enrichment_df$c,
      curve = "FC+SC",
      power = enrichment_df$power_fc_sc_same_region
    )
  )

  p_power <- ggplot2::ggplot(
    power_long,
    ggplot2::aes(x = c, y = power, linetype = curve)
  ) +
    ggplot2::geom_line(linewidth = 0.5, color = "black") +
    ggplot2::scale_linetype_manual(
      values = c("FC+SC" = "solid", "FC" = "dotted")
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::coord_cartesian(xlim = c(0.1, 0.4)) +
    ggplot2::labs(
      title = "Comparison of power curves of FC and FC+SC",
      subtitle = bquote(pi[0] == .(round(pi0, 2))),
      #subtitle = sprintf("Both curves use the same FC-only rejection region A_k; pi0 = %.2f", pi0),
      x = "LFDR threshold c",
      y = "Power",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom"
    )

  p_difference <- ggplot2::ggplot(
    enrichment_df,
    ggplot2::aes(x = c, y = power_difference)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "gray45", linetype = "dashed") +
    ggplot2::geom_line(color = "#2C7FB8", linewidth = 1) +
    ggplot2::labs(
      title = "Theorem 1 Difference",
      subtitle = sprintf("Power difference: FC+SC minus FC; FDR alpha level = %.2f", alpha_level),
      x = "LFDR threshold c",
      y = "Power difference"
    ) +
    ggplot2::coord_cartesian(xlim = c(0.1, 0.4)) +
    ggplot2::theme_classic(base_size = 12)

  list(
    power_plot = p_power,
    difference_plot = p_difference,
    enrichment_data = enrichment_df
  )
}

# Three power curves for the baseline, prior-refined, and jointly enriched models.
#
# Bottom:
#   E_f1(t1)[I_Ak(c, pi0)(t1)]
# Middle:
#   E_f1(t1)[I_Bk(c, pi0(t2))(t1)]
# Top:
#   E_f1(t1,t2)[I_Bk(c, pi0(t2))(t1)]
#
# Importantly, the top and middle curves use the same B_k rejection
# threshold. The secondary modality changes the non-null distribution
# in the top curve, but it does not change the B_k threshold again.
compute_three_model_power_curves <- function(c_values,
                                      alpha0,
                                      beta0,
                                      mu1,
                                      alpha1,
                                      beta1,
                                      alpha_s,
                                      beta_s,
                                      pi0 = pi0_default,
                                      pi0_t2 = pi0_t2_default,
                                      tolerance = 1e-12) {
  if (length(c_values) == 0 || any(!is.finite(c_values))) {
    stop("c_values must contain at least one finite value.")
  }
  if (any(c_values <= 0 | c_values >= 1)) {
    stop("All c_values must be in (0, 1).")
  }
  if (pi0 <= 0 || pi0 >= 1 || pi0_t2 <= 0 || pi0_t2 >= 1) {
    stop("pi0 and pi0_t2 must both be in (0, 1).")
  }
  if (pi0_t2 > pi0) {
    stop(
      "The prior-refined model requires pi0_t2 <= pi0; received pi0_t2 = ",
      pi0_t2,
      " and pi0 = ",
      pi0,
      "."
    )
  }

  density_parameters <- c(
    alpha0 = alpha0,
    beta0 = beta0,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s
  )
  if (any(!is.finite(density_parameters)) || any(density_parameters <= 0)) {
    stop("All Gamma shape and rate parameters must be finite and positive.")
  }
  if (!is.finite(mu1) || mu1 < 0) {
    stop("mu1 must be finite and non-negative.")
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be finite and non-negative.")
  }

  alpha_bimodal <- alpha1 + alpha_s / beta_s

  rows <- lapply(c_values, function(c_val) {
    # A_k uses the marginal null proportion and the FC-only likelihood ratio.
    t_star_A <- find_t_star_gamma(
      c = c_val,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1,
      pi0 = pi0
    )

    # B_k uses pi0(t2), but its likelihood ratio is still f1(t1) / f0(t1).
    t_star_B <- find_t_star_gamma(
      c = c_val,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1,
      pi0 = pi0_t2
    )

    power_bottom <- power_shifted_gamma(
      t_star = t_star_A,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1
    )

    power_middle <- power_shifted_gamma(
      t_star = t_star_B,
      mu1 = mu1,
      alpha_alt = alpha1,
      beta_alt = beta1
    )

    power_top <- power_shifted_gamma(
      t_star = t_star_B,
      mu1 = mu1,
      alpha_alt = alpha_bimodal,
      beta_alt = beta1
    )

    data.frame(
      c = c_val,
      pi0 = pi0,
      pi0_t2 = pi0_t2,
      k_A = k_fun(c_val, pi0 = pi0),
      k_B = k_fun(c_val, pi0 = pi0_t2),
      t_star_A = t_star_A,
      t_star_B = t_star_B,
      power_bottom = power_bottom,
      power_middle = power_middle,
      power_top = power_top,
      top_minus_middle = power_top - power_middle,
      middle_minus_bottom = power_middle - power_bottom,
      region_nested = t_star_B <= t_star_A + tolerance,
      first_inequality_holds = power_top >= power_middle - tolerance,
      second_inequality_holds = power_middle >= power_bottom - tolerance,
      ordering_holds =
        power_top >= power_middle - tolerance &
        power_middle >= power_bottom - tolerance
    )
  })

  do.call(rbind, rows)
}

plot_three_model_power_curves <- function(c_values = seq(0.01, 0.99, by = 0.01),
                                   alpha0,
                                   beta0,
                                   mu1,
                                   alpha1,
                                   beta1,
                                   alpha_s,
                                   beta_s,
                                   pi0 = pi0_default,
                                   pi0_t2 = pi0_t2_default,
                                   c_star = NULL,
                                   print_c_star = TRUE,
                                   x_limits = NULL,
                                   tolerance = 1e-12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_three_model_power_curves().")
  }
  if (is.list(c_star) && !is.null(c_star$c_star)) {
    c_star <- c_star$c_star
  }
  if (!is.null(c_star)) {
    if (length(c_star) != 1 || !is.finite(c_star) || c_star <= 0 || c_star >= 1) {
      stop("c_star must be NULL or one finite value in (0, 1).")
    }
  }
  if (!is.null(x_limits)) {
    if (length(x_limits) != 2 || any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be NULL or a finite increasing vector of length two.")
    }
  }

  power_comparison_df <- compute_three_model_power_curves(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi0 = pi0,
    pi0_t2 = pi0_t2,
    tolerance = tolerance
  )

  curve_levels <- c(
    "top_bimodal_and_pi0_t2",
    "middle_pi0_t2_only",
    "bottom_fc_only"
  )

  curve_labels <- c(
    "Top: f1(t1, t2), Bk with pi0(t2)",
    "Middle: f1(t1), Bk with pi0(t2)",
    "Bottom: f1(t1), Ak with pi0"
  )

  power_long <- rbind(
    data.frame(
      c = power_comparison_df$c,
      curve_id = "top_bimodal_and_pi0_t2",
      power = power_comparison_df$power_top
    ),
    data.frame(
      c = power_comparison_df$c,
      curve_id = "middle_pi0_t2_only",
      power = power_comparison_df$power_middle
    ),
    data.frame(
      c = power_comparison_df$c,
      curve_id = "bottom_fc_only",
      power = power_comparison_df$power_bottom
    )
  )

  power_long$curve <- factor(
    power_long$curve_id,
    levels = curve_levels,
    labels = curve_labels
  )

  curve_plotmath_labels <- expression(
    "Top: " * f[1](t[1], t[2]) * ", " * B[k] * " with " * pi[0](t[2]),
    "Middle: " * f[1](t[1]) * ", " * B[k] * " with " * pi[0](t[2]),
    "Bottom: " * f[1](t[1]) * ", " * A[k] * " with " * pi[0]
  )

  difference_long <- rbind(
    data.frame(
      c = power_comparison_df$c,
      difference_id = "top_minus_middle",
      difference = power_comparison_df$top_minus_middle
    ),
    data.frame(
      c = power_comparison_df$c,
      difference_id = "middle_minus_bottom",
      difference = power_comparison_df$middle_minus_bottom
    )
  )

  difference_long$comparison <- factor(
    difference_long$difference_id,
    levels = c("top_minus_middle", "middle_minus_bottom"),
    labels = c(
      "Top minus middle",
      "Middle minus bottom"
    )
  )

  c_star_power_comparison <- NULL
  c_star_power_data <- data.frame()

  if (!is.null(c_star)) {
    c_star_power_comparison <- compute_three_model_power_curves(
      c_values = c_star,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s,
      pi0 = pi0,
      pi0_t2 = pi0_t2,
      tolerance = tolerance
    )

    c_star_power_data <- data.frame(
      c_star = rep(c_star, 3),
      curve_id = curve_levels,
      curve = factor(curve_labels, levels = curve_labels),
      power = c(
        c_star_power_comparison$power_top,
        c_star_power_comparison$power_middle,
        c_star_power_comparison$power_bottom
      )
    )

    if (print_c_star) {
      cat("\nModel powers at c_star = ", signif(c_star, 6), "\n", sep = "")
      print(
        data.frame(
          curve = curve_labels,
          power_at_c_star = signif(c_star_power_data$power, 6)
        ),
        row.names = FALSE
      )
      cat("\n")
    }

    if (!is.null(x_limits) && (c_star < x_limits[1] || c_star > x_limits[2])) {
      warning("c_star is outside x_limits, so the vertical line will not be visible.")
    }
  }

  power_plot <- ggplot2::ggplot(
    power_long,
    ggplot2::aes(x = c, y = power, linetype = curve)
  ) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Top: f1(t1, t2), Bk with pi0(t2)" = "solid",
        "Middle: f1(t1), Bk with pi0(t2)" = "longdash",
        "Bottom: f1(t1), Ak with pi0" = "dotted"
      ),
      breaks = curve_labels,
      labels = curve_plotmath_labels
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggplot2::labs(
      title = expression("Power curve comparison using SC in" ~ pi[0] ~ "and non-null density"),
      subtitle = bquote(
        pi[0] == .(pi0) * ";" ~~
          pi[0](t[2]) == .(round(pi0_t2, 2))
      ),
      x = "LFDR threshold c",
      y = "Power",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0)
  )
  if (!is.null(c_star)) {
    power_plot <- power_plot +
      ggplot2::geom_vline(
        xintercept = c_star,
        color = "gray35",
        linetype = "dotdash",
        linewidth = 0.55
      ) +
      ggplot2::annotate(
        "text",
        x = c_star+0.02,
        y = 0.035,
        label = paste0("c*=", formatC(c_star, format = "f", digits = 2)),
        size = 4,
        color = "gray25",
        vjust = 0
      )
  }

  difference_plot <- ggplot2::ggplot(
    difference_long,
    ggplot2::aes(x = c, y = difference, linetype = comparison)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "gray55",
      linetype = "dashed"
    ) +
    ggplot2::geom_line(color = "black", linewidth = 0.75) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Top minus middle" = "solid",
        "Middle minus bottom" = "longdash"
      )
    ) +
    ggplot2::labs(
      title = "Pairwise Power Differences",
      subtitle = "Nonnegative differences indicate the expected ordering",
      x = "LFDR threshold c",
      y = "Power difference",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(x_limits)) {
    power_plot <- power_plot + ggplot2::coord_cartesian(xlim = x_limits)
    difference_plot <- difference_plot + ggplot2::coord_cartesian(xlim = x_limits)
  }

  list(
    power_plot = power_plot,
    difference_plot = difference_plot,
    power_data = power_long,
    difference_data = difference_long,
    c_star_power_data = c_star_power_data,
    c_star_power_comparison = c_star_power_comparison,
    power_comparison_data = power_comparison_df,
    region_nested = all(power_comparison_df$region_nested),
    first_ordering_holds = all(power_comparison_df$first_inequality_holds),
    second_ordering_holds = all(power_comparison_df$second_inequality_holds),
    ordering_holds = all(power_comparison_df$ordering_holds)
  )
}

# TP2 power-ratio inequality across prior-proportion levels.
#
#   G(c, pi02, t22)   G(c, pi02, t21)
#   --------------- >= ---------------
#   G(c, pi01, t22)   G(c, pi01, t21)
#
# with pi02 >= pi01 and t22 >= t21.
#
# In the numerical model:
#   t21 = 0  means SC absent, so alpha(t21) = alpha1.
#   t22 = 1  means SC present at the average modeled effect, so
#             alpha(t22) = alpha1 + alpha_s / beta_s.
#
auxiliary_shape <- function(t2, alpha1, alpha_s, beta_s) {
  if (any(!is.finite(t2)) || any(t2 < 0)) {
    stop("t2 must be finite and non-negative.")
  }
  if (!is.finite(alpha1) || !is.finite(alpha_s) || !is.finite(beta_s) ||
      alpha1 <= 0 || alpha_s <= 0 || beta_s <= 0) {
    stop("alpha1, alpha_s, and beta_s must be finite and positive.")
  }

  alpha1 + alpha_s * t2 / beta_s
}

power_at_threshold <- function(c,
                               pi0,
                               t2,
                               alpha0,
                               beta0,
                               mu1,
                               alpha1,
                               beta1,
                               alpha_s,
                               beta_s) {
  if (length(c) != 1 || !is.finite(c) || c <= 0 || c >= 1) {
    stop("c must be one finite value in (0, 1).")
  }
  if (length(pi0) != 1 || !is.finite(pi0) || pi0 <= 0 || pi0 >= 1) {
    stop("pi0 must be one finite value in (0, 1).")
  }
  if (length(t2) != 1 || !is.finite(t2) || t2 < 0) {
    stop("t2 must be one finite non-negative value.")
  }

  t_star <- find_t_star_gamma(
    c = c,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha_alt = alpha1,
    beta_alt = beta1,
    pi0 = pi0
  )

  alpha_t2 <- auxiliary_shape(
    t2 = t2,
    alpha1 = alpha1,
    alpha_s = alpha_s,
    beta_s = beta_s
  )

  power_shifted_gamma(
    t_star = t_star,
    mu1 = mu1,
    alpha_alt = alpha_t2,
    beta_alt = beta1
  )
}

compute_tp2_prior_ratio_inequality <- function(c_values,
                                             alpha0,
                                             beta0,
                                             mu1,
                                             alpha1,
                                             beta1,
                                             alpha_s,
                                             beta_s,
                                             pi02 = pi0_default,
                                             pi01 = pi0_t2_default,
                                             t21 = 0,
                                             t22 = 1,
                                             tolerance = 1e-12) {
  if (length(c_values) == 0 || any(!is.finite(c_values))) {
    stop("c_values must contain at least one finite value.")
  }
  if (any(c_values <= 0 | c_values >= 1)) {
    stop("All c_values must be in (0, 1).")
  }
  if (pi01 <= 0 || pi01 >= 1 || pi02 <= 0 || pi02 >= 1) {
    stop("pi01 and pi02 must both be in (0, 1).")
  }
  if (pi02 < pi01) {
    stop("pi02 must be greater than or equal to pi01.")
  }
  if (t22 < t21) {
    stop("t22 must be greater than or equal to t21.")
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be finite and non-negative.")
  }

  safe_ratio <- function(numerator, denominator) {
    ifelse(
      is.finite(numerator) & is.finite(denominator) & denominator > 0,
      numerator / denominator,
      NA_real_
    )
  }

  rows <- lapply(c_values, function(c_val) {
    G_pi02_t22 <- power_at_threshold(
      c = c_val,
      pi0 = pi02,
      t2 = t22,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s
    )

    G_pi01_t22 <- power_at_threshold(
      c = c_val,
      pi0 = pi01,
      t2 = t22,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s
    )

    G_pi02_t21 <- power_at_threshold(
      c = c_val,
      pi0 = pi02,
      t2 = t21,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s
    )

    G_pi01_t21 <- power_at_threshold(
      c = c_val,
      pi0 = pi01,
      t2 = t21,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s
    )

    ratio_sc_present <- safe_ratio(G_pi02_t22, G_pi01_t22)
    ratio_sc_absent <- safe_ratio(G_pi02_t21, G_pi01_t21)
    ratio_difference <- ratio_sc_present - ratio_sc_absent

    data.frame(
      c = c_val,
      pi02 = pi02,
      pi01 = pi01,
      t22 = t22,
      t21 = t21,
      alpha_t22 = auxiliary_shape(t22, alpha1, alpha_s, beta_s),
      alpha_t21 = auxiliary_shape(t21, alpha1, alpha_s, beta_s),
      G_pi02_t22 = G_pi02_t22,
      G_pi01_t22 = G_pi01_t22,
      G_pi02_t21 = G_pi02_t21,
      G_pi01_t21 = G_pi01_t21,
      ratio_sc_present = ratio_sc_present,
      ratio_sc_absent = ratio_sc_absent,
      ratio_difference = ratio_difference,
      ordering_holds = is.finite(ratio_difference) &
        ratio_difference >= -tolerance
    )
  })

  do.call(rbind, rows)
}

plot_tp2_prior_ratio_inequality <- function(c_values = seq(0.01, 0.99, by = 0.01),
                                          alpha0,
                                          beta0,
                                          mu1,
                                          alpha1,
                                          beta1,
                                          alpha_s,
                                          beta_s,
                                          pi02 = pi0_default,
                                          pi01 = pi0_t2_default,
                                          t21 = 0,
                                          t22 = 1,
                                          c_star = NULL,
                                          x_limits = NULL,
                                          print_c_star = TRUE,
                                          print_summary = TRUE,
                                          tolerance = 1e-12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_tp2_prior_ratio_inequality().")
  }
  if (is.list(c_star) && !is.null(c_star$c_star)) {
    c_star <- c_star$c_star
  }
  if (!is.null(c_star)) {
    if (length(c_star) != 1 || !is.finite(c_star) || c_star <= 0 || c_star >= 1) {
      stop("c_star must be NULL or one finite value in (0, 1).")
    }
  }
  if (!is.null(x_limits)) {
    if (length(x_limits) != 2 || any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be NULL or a finite increasing vector of length two.")
    }
  }

  tp2_prior_df <- compute_tp2_prior_ratio_inequality(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi02 = pi02,
    pi01 = pi01,
    t21 = t21,
    t22 = t22,
    tolerance = tolerance
  )

  if (print_summary) {
    cat("\nTP2 prior-ratio numerical check\n")
    cat("pi02 = ", pi02, ", pi01 = ", pi01, ", t22 = ", t22, ", t21 = ", t21, "\n", sep = "")
    cat("All finite checks hold: ", all(tp2_prior_df$ordering_holds, na.rm = TRUE), "\n", sep = "")
    cat("Minimum ratio difference = ", signif(min(tp2_prior_df$ratio_difference, na.rm = TRUE), 6), "\n\n", sep = "")
  }

  ratio_curve_levels <- c(
    "SC present: G(c, pi02, t22) / G(c, pi01, t22)",
    "SC absent: G(c, pi02, t21) / G(c, pi01, t21)"
  )

  ratio_curve_labels <- expression(
    "SC present: " * G(c, pi["02"], t["22"]) * " / " * G(c, pi["01"], t["22"]),
    "SC absent: " * G(c, pi["02"], t["21"]) * " / " * G(c, pi["01"], t["21"])
  )

  ratio_long <- rbind(
    data.frame(
      c = tp2_prior_df$c,
      ratio = tp2_prior_df$ratio_sc_present,
      curve = "SC present: G(c, pi02, t22) / G(c, pi01, t22)"
    ),
    data.frame(
      c = tp2_prior_df$c,
      ratio = tp2_prior_df$ratio_sc_absent,
      curve = "SC absent: G(c, pi02, t21) / G(c, pi01, t21)"
    )
  )

  ratio_long$curve <- factor(ratio_long$curve, levels = ratio_curve_levels)

  c_star_tp2_prior_data <- NULL
  c_star_ratio_data <- data.frame()

  if (!is.null(c_star)) {
    c_star_tp2_prior_data <- compute_tp2_prior_ratio_inequality(
      c_values = c_star,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s,
      pi02 = pi02,
      pi01 = pi01,
      t21 = t21,
      t22 = t22,
      tolerance = tolerance
    )

    c_star_ratio_data <- data.frame(
      c_star = rep(c_star, 2),
      curve = factor(ratio_curve_levels, levels = ratio_curve_levels),
      ratio = c(
        c_star_tp2_prior_data$ratio_sc_present,
        c_star_tp2_prior_data$ratio_sc_absent
      )
    )

    if (print_c_star) {
      cat("\nTP2 prior-ratio power ratios at c_star = ", signif(c_star, 6), "\n", sep = "")
      print(
        data.frame(
          curve = ratio_curve_levels,
          power_ratio_at_c_star = signif(c_star_ratio_data$ratio, 6)
        ),
        row.names = FALSE
      )
      cat("\n")
    }

    if (!is.null(x_limits) && (c_star < x_limits[1] || c_star > x_limits[2])) {
      warning("c_star is outside x_limits, so the vertical line will not be visible.")
    }
  }

  ratio_plot <- ggplot2::ggplot(
    ratio_long,
    ggplot2::aes(x = c, y = ratio, linetype = curve)
  ) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::scale_linetype_manual(
      values = c(
        "SC present: G(c, pi02, t22) / G(c, pi01, t22)" = "solid",
        "SC absent: G(c, pi02, t21) / G(c, pi01, t21)" = "longdash"
      ),
      breaks = ratio_curve_levels,
      labels = ratio_curve_labels
    ) +
    ggplot2::labs(
      title = "Comparison of power ratio curves with and without SC",
      subtitle = bquote(
        pi["02"] == .(round(pi02, 2)) * ";" ~~
          pi["01"] == .(round(pi01, 2)) * ";" ~~
          t[21] == .(t21) * "," ~~
          t[22] == .(t22)
      ),
      x = "LFDR threshold c",
      y = "Power ratio",
      linetype = NULL
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(
          color = "black",
          linewidth = 0.75,
          linetype = c("solid", "longdash")
        )
      )
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.width = grid::unit(1.5, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(c_star)) {
    ratio_plot <- ratio_plot +
      ggplot2::geom_vline(
        xintercept = c_star,
        color = "gray35",
        linetype = "dotdash",
        linewidth = 0.55
      ) +
      ggplot2::annotate(
        "text",
        x = c_star+0.045,
        y = min(ratio_long$ratio, na.rm = TRUE)+0.2,
        label = paste0("c*=", formatC(c_star, format = "f", digits = 2)),
        size = 4,
        color = "gray25",
        vjust = -0.6
      )
  }

  difference_plot <- ggplot2::ggplot(
    tp2_prior_df,
    ggplot2::aes(x = c, y = ratio_difference)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "gray55",
      linetype = "dashed"
    ) +
    ggplot2::geom_line(color = "black", linewidth = 0.75) +
    ggplot2::labs(
      title = "TP2 Prior-Ratio Difference",
      subtitle = "SC-present ratio minus SC-absent ratio",
      x = "LFDR threshold c",
      y = "Ratio difference"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(x_limits)) {
    ratio_plot <- ratio_plot + ggplot2::coord_cartesian(xlim = x_limits)
    difference_plot <- difference_plot + ggplot2::coord_cartesian(xlim = x_limits)
  }

  list(
    ratio_plot = ratio_plot,
    difference_plot = difference_plot,
    c_star_ratio_data = c_star_ratio_data,
    c_star_tp2_prior_data = c_star_tp2_prior_data,
    tp2_prior_data = tp2_prior_df,
    ordering_holds = all(tp2_prior_df$ordering_holds, na.rm = TRUE)
  )
}

# TP2 power-ratio inequality across auxiliary-information levels.
#
#   G(c, pi02, t22)   G(c, pi01, t22)
#   --------------- >= ---------------
#   G(c, pi02, t21)   G(c, pi01, t21)
#
# with pi02 >= pi01 and t22 >= t21.
#
# This is the complementary ratio view: instead of comparing
# pi02/pi01 power ratios across SC levels; the complementary comparison
# present/absent power ratio across null-proportion levels.
compute_tp2_auxiliary_ratio_inequality <- function(c_values,
                                             alpha0,
                                             beta0,
                                             mu1,
                                             alpha1,
                                             beta1,
                                             alpha_s,
                                             beta_s,
                                             pi02 = pi0_default,
                                             pi01 = pi0_t2_default,
                                             t21 = 0,
                                             t22 = 1,
                                             tolerance = 1e-12) {
  if (length(c_values) == 0 || any(!is.finite(c_values))) {
    stop("c_values must contain at least one finite value.")
  }
  if (any(c_values <= 0 | c_values >= 1)) {
    stop("All c_values must be in (0, 1).")
  }
  if (pi01 <= 0 || pi01 >= 1 || pi02 <= 0 || pi02 >= 1) {
    stop("pi01 and pi02 must both be in (0, 1).")
  }
  if (pi02 < pi01) {
    stop("pi02 must be greater than or equal to pi01.")
  }
  if (t22 < t21) {
    stop("t22 must be greater than or equal to t21.")
  }
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be finite and non-negative.")
  }

  tp2_prior_df <- compute_tp2_prior_ratio_inequality(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi02 = pi02,
    pi01 = pi01,
    t21 = t21,
    t22 = t22,
    tolerance = tolerance
  )

  safe_ratio <- function(numerator, denominator) {
    ifelse(
      is.finite(numerator) & is.finite(denominator) & denominator > 0,
      numerator / denominator,
      NA_real_
    )
  }

  ratio_high_pi0 <- safe_ratio(tp2_prior_df$G_pi02_t22, tp2_prior_df$G_pi02_t21)
  ratio_low_pi0 <- safe_ratio(tp2_prior_df$G_pi01_t22, tp2_prior_df$G_pi01_t21)
  ratio_difference <- ratio_high_pi0 - ratio_low_pi0

  data.frame(
    c = tp2_prior_df$c,
    pi02 = tp2_prior_df$pi02,
    pi01 = tp2_prior_df$pi01,
    t22 = tp2_prior_df$t22,
    t21 = tp2_prior_df$t21,
    alpha_t22 = tp2_prior_df$alpha_t22,
    alpha_t21 = tp2_prior_df$alpha_t21,
    G_pi02_t22 = tp2_prior_df$G_pi02_t22,
    G_pi01_t22 = tp2_prior_df$G_pi01_t22,
    G_pi02_t21 = tp2_prior_df$G_pi02_t21,
    G_pi01_t21 = tp2_prior_df$G_pi01_t21,
    ratio_high_pi0 = ratio_high_pi0,
    ratio_low_pi0 = ratio_low_pi0,
    ratio_difference = ratio_difference,
      ordering_holds = is.finite(ratio_difference) &
      ratio_difference >= -tolerance
  )
}

plot_tp2_auxiliary_ratio_inequality <- function(c_values = seq(0.01, 0.99, by = 0.01),
                                          alpha0,
                                          beta0,
                                          mu1,
                                          alpha1,
                                          beta1,
                                          alpha_s,
                                          beta_s,
                                          pi02 = pi0_default,
                                          pi01 = pi0_t2_default,
                                          t21 = 0,
                                          t22 = 1,
                                          c_star = NULL,
                                          x_limits = NULL,
                                          print_c_star = TRUE,
                                          print_summary = TRUE,
                                          tolerance = 1e-12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_tp2_auxiliary_ratio_inequality().")
  }
  if (is.list(c_star) && !is.null(c_star$c_star)) {
    c_star <- c_star$c_star
  }
  if (!is.null(c_star)) {
    if (length(c_star) != 1 || !is.finite(c_star) || c_star <= 0 || c_star >= 1) {
      stop("c_star must be NULL or one finite value in (0, 1).")
    }
  }
  if (!is.null(x_limits)) {
    if (length(x_limits) != 2 || any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be NULL or a finite increasing vector of length two.")
    }
  }

  tp2_auxiliary_df <- compute_tp2_auxiliary_ratio_inequality(
    c_values = c_values,
    alpha0 = alpha0,
    beta0 = beta0,
    mu1 = mu1,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha_s = alpha_s,
    beta_s = beta_s,
    pi02 = pi02,
    pi01 = pi01,
    t21 = t21,
    t22 = t22,
    tolerance = tolerance
  )

  if (print_summary) {
    cat("\nTP2 auxiliary-ratio numerical check\n")
    cat("pi02 = ", pi02, ", pi01 = ", pi01, ", t22 = ", t22, ", t21 = ", t21, "\n", sep = "")
    cat("All finite checks hold: ", all(tp2_auxiliary_df$ordering_holds, na.rm = TRUE), "\n", sep = "")
    cat("Minimum ratio difference = ", signif(min(tp2_auxiliary_df$ratio_difference, na.rm = TRUE), 6), "\n\n", sep = "")
  }

  ratio_curve_levels <- c(
    "High pi0: G(c, pi02, t22) / G(c, pi02, t21)",
    "Low pi0: G(c, pi01, t22) / G(c, pi01, t21)"
  )

  ratio_curve_labels <- expression(
    "High " * pi[0] * ": " * G(c, pi["02"], t["22"]) * " / " * G(c, pi["02"], t["21"]),
    "Low " * pi[0] * ": " * G(c, pi["01"], t["22"]) * " / " * G(c, pi["01"], t["21"])
  )

  ratio_long <- rbind(
    data.frame(
      c = tp2_auxiliary_df$c,
      ratio = tp2_auxiliary_df$ratio_high_pi0,
      curve = "High pi0: G(c, pi02, t22) / G(c, pi02, t21)"
    ),
    data.frame(
      c = tp2_auxiliary_df$c,
      ratio = tp2_auxiliary_df$ratio_low_pi0,
      curve = "Low pi0: G(c, pi01, t22) / G(c, pi01, t21)"
    )
  )

  ratio_long$curve <- factor(ratio_long$curve, levels = ratio_curve_levels)

  c_star_tp2_auxiliary_data <- NULL
  c_star_ratio_data <- data.frame()

  if (!is.null(c_star)) {
    c_star_tp2_auxiliary_data <- compute_tp2_auxiliary_ratio_inequality(
      c_values = c_star,
      alpha0 = alpha0,
      beta0 = beta0,
      mu1 = mu1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha_s = alpha_s,
      beta_s = beta_s,
      pi02 = pi02,
      pi01 = pi01,
      t21 = t21,
      t22 = t22,
      tolerance = tolerance
    )

    c_star_ratio_data <- data.frame(
      c_star = rep(c_star, 2),
      curve = factor(ratio_curve_levels, levels = ratio_curve_levels),
      ratio = c(
        c_star_tp2_auxiliary_data$ratio_high_pi0,
        c_star_tp2_auxiliary_data$ratio_low_pi0
      )
    )

    if (print_c_star) {
      cat("\nTP2 auxiliary-ratio power ratios at c_star = ", signif(c_star, 6), "\n", sep = "")
      print(
        data.frame(
          curve = ratio_curve_levels,
          power_ratio_at_c_star = signif(c_star_ratio_data$ratio, 6)
        ),
        row.names = FALSE
      )
      cat("\n")
    }

    if (!is.null(x_limits) && (c_star < x_limits[1] || c_star > x_limits[2])) {
      warning("c_star is outside x_limits, so the vertical line will not be visible.")
    }
  }

  ratio_plot <- ggplot2::ggplot(
    ratio_long,
    ggplot2::aes(x = c, y = ratio, linetype = curve)
  ) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::scale_linetype_manual(
      values = c(
        "High pi0: G(c, pi02, t22) / G(c, pi02, t21)" = "solid",
        "Low pi0: G(c, pi01, t22) / G(c, pi01, t21)" = "longdash"
      ),
      breaks = ratio_curve_levels,
      labels = ratio_curve_labels
    ) +
    ggplot2::labs(
      title = expression("Comparison of power ratio curves with and without using SC in" ~ pi[0]),
      subtitle = bquote(
        pi["02"] == .(round(pi02, 2)) * ";" ~~
          pi["01"] == .(round(pi01, 2)) * ";" ~~
          t[21] == .(t21) * "," ~~
          t[22] == .(t22)
      ),
      x = "LFDR threshold c",
      y = "Power ratio",
      linetype = NULL
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(
          color = "black",
          linewidth = 0.75,
          linetype = c("solid", "longdash")
        )
      )
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.width = grid::unit(1.5, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(c_star)) {
    ratio_plot <- ratio_plot +
      ggplot2::geom_vline(
        xintercept = c_star,
        color = "gray35",
        linetype = "dotdash",
        linewidth = 0.55
      ) +
      ggplot2::annotate(
        "text",
        x = c_star + 0.065,
        y = min(ratio_long$ratio, na.rm = TRUE) + 50,
        label = paste0("c*=", formatC(c_star, format = "f", digits = 2)),
        size = 4,
        color = "gray25",
        vjust = -0.6
      )
  }

  difference_plot <- ggplot2::ggplot(
    tp2_auxiliary_df,
    ggplot2::aes(x = c, y = ratio_difference)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "gray55",
      linetype = "dashed"
    ) +
    ggplot2::geom_line(color = "black", linewidth = 0.75) +
    ggplot2::labs(
      title = "TP2 Auxiliary-Ratio Difference",
      subtitle = "High-pi0 ratio minus low-pi0 ratio",
      x = "LFDR threshold c",
      y = "Ratio difference"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(x_limits)) {
    ratio_plot <- ratio_plot + ggplot2::coord_cartesian(xlim = x_limits)
    difference_plot <- difference_plot + ggplot2::coord_cartesian(xlim = x_limits)
  }

  list(
    ratio_plot = ratio_plot,
    difference_plot = difference_plot,
    c_star_ratio_data = c_star_ratio_data,
    c_star_tp2_auxiliary_data = c_star_tp2_auxiliary_data,
    tp2_auxiliary_data = tp2_auxiliary_df,
    ordering_holds = all(tp2_auxiliary_df$ordering_holds, na.rm = TRUE)
  )
}

# Example using the current RData-based analytical dataset:
#
# load("path/to/analysis.RData")
# fmri.df <- fMRI.est.all
# dti.df <- DTI.est.all
#
# analytical_dat <- create_analytical_dataset(fmri.df, dti.df, p_cutoff = 0.05)
# est <- estimate_mom_from_analytical_dat(analytical_dat, mu1 = 2)
#
# alpha0_mom <- est$alpha0_mom
# beta0_mom <- est$beta0_mom
# alpha1_mom <- est$alpha1_mom
# beta1_mom <- est$beta1_mom
# alpha_s_mom <- est$alpha_s_mom
# beta_s_mom <- est$beta_s_mom
# alpha_bimodal_mom <- est$alpha_bimodal_mom
# mu1 <- est$mu1
#
# power_df <- compute_power_three(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi0 = 0.93
# )
#
# c_star_fc <- find_c_star_from_t1(
#   t1 = analytical_dat$fmri_t_value,
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha_alt = alpha1_mom,
#   beta_alt = beta1_mom,
#   pi0 = 0.93,
#   alpha = 0.05
# )
#
# power_plot <- plot_power_curves_three(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi0 = 0.93,
#   alpha_level = 0.05,
#   t1 = analytical_dat$fmri_t_value
# )
# power_plot$plot
# power_plot$c_star_data
#
# enrichment_plot <- plot_enrichment_power_inequality(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi0 = 0.93,
#   alpha_level = 0.05
# )
# enrichment_plot$power_plot
# enrichment_plot$difference_plot
# all(enrichment_plot$enrichment_data$ordering_holds)
#
# three_model_plot <- plot_three_model_power_curves(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi0 = 0.93,
#   pi0_t2 = 0.9083,
#   c_star = c_star_fc,
#   x_limits = c(0.10, 0.40)
# )
# three_model_plot$power_plot
# three_model_plot$c_star_power_data
# three_model_plot$difference_plot
# three_model_plot$ordering_holds
#
# tp2_prior_plot <- plot_tp2_prior_ratio_inequality(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi02 = 0.93,
#   pi01 = 0.9083,
#   t21 = 0,
#   t22 = 1,
#   c_star = c_star_fc,
#   x_limits = c(0.10, 0.40)
# )
# tp2_prior_plot$ratio_plot
# tp2_prior_plot$c_star_ratio_data
# tp2_prior_plot$difference_plot
# tp2_prior_plot$ordering_holds
# summary(tp2_prior_plot$tp2_prior_data$ratio_difference)
#
# tp2_auxiliary_plot <- plot_tp2_auxiliary_ratio_inequality(
#   c_values = seq(0.01, 0.99, by = 0.01),
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha1 = alpha1_mom,
#   beta1 = beta1_mom,
#   alpha_s = alpha_s_mom,
#   beta_s = beta_s_mom,
#   pi02 = 0.93,
#   pi01 = 0.9083,
#   t21 = 0,
#   t22 = 1,
#   c_star = c_star_fc,
#   x_limits = c(0.10, 0.40)
# )
# tp2_auxiliary_plot$ratio_plot
# tp2_auxiliary_plot$c_star_ratio_data
# tp2_auxiliary_plot$difference_plot
# tp2_auxiliary_plot$ordering_holds
# summary(tp2_auxiliary_plot$tp2_auxiliary_data$ratio_difference)
#
# c_star_fc_sc <- find_c_star_from_t1(
#   t1 = analytical_dat$fmri_t_value,
#   alpha0 = alpha0_mom,
#   beta0 = beta0_mom,
#   mu1 = mu1,
#   alpha_alt = alpha_bimodal_mom,
#   beta_alt = beta1_mom,
#   pi0 = 0.93,
#   alpha = 0.05
# )
