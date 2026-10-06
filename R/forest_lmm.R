library(dplyr)
library(ggplot2)
library(patchwork)

# ---- 1. Collect results from all models -------------------------------------
# Each run_lmm() result has $results (term, estimate, lo, hi, p)
# and $model_info (n_ids, n_obs).
collect_lmm <- function(models, term = "klot_elisa:time") {
  bind_rows(lapply(names(models), function(nm) {
    m <- models[[nm]]
    m$results %>%
      filter(.data$term == !!term) %>%
      mutate(model = nm,
             n_ids = m$model_info$n_ids,
             n_obs = m$model_info$n_obs)
  }))
}

# ---- 2. Forest plot ---------------------------------------------------------
# rows: tibble with columns model (name in `models`), label, group
plot_forest_lmm <- function(models, rows, term = "klot_elisa:time",
                            xlab = "Difference in yearly change per SD Klotho (points)",
                            est_lab = "β (95% CI)",
                            digits = 2, xlim = NULL) {

  res <- collect_lmm(models, term)

  # one row per model + one header row per group, in plotting order
  plot_df <- rows %>%
    left_join(res, by = "model") %>%
    group_by(group = factor(group, levels = unique(rows$group))) %>%
    group_modify(~ bind_rows(tibble(label = as.character(.y$group), header = TRUE),
                             mutate(.x, header = FALSE))) %>%
    ungroup() %>%
    mutate(
      row    = rev(seq_len(n())),                          # top to bottom
      crude  = !header & grepl("crude", label, ignore.case = TRUE),
      est_txt = if_else(header, "",
                        sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"),
                                estimate, lo, hi)),
      n_txt  = if_else(header, "",
                       paste0(prettyNum(n_ids, big.mark = ","), " / ", prettyNum(n_obs, big.mark = ",")))
    )

  pts <- filter(plot_df, !header)
  ylim <- c(0.5, max(plot_df$row) + 0.9)   # headroom for column titles
  col_main <- "#2C5F8A"; col_crude <- "#8C949D"; txt <- "#2B3440"

  base_theme <- theme_void(base_size = 13) +
    theme(plot.margin = margin(5, 5, 5, 5))

  # left: labels
  p_lab <- ggplot(plot_df, aes(y = row)) +
    geom_text(aes(x = if_else(header, 0, 0.08), label = label,
                  fontface = if_else(header, "bold", "plain")),
              hjust = 0, size = 4.2, colour = txt) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    base_theme

  # middle: estimates
  rng <- if (is.null(xlim)) range(c(pts$lo, pts$hi, 0), na.rm = TRUE) else xlim
  p_fp <- ggplot(pts, aes(y = row, x = estimate)) +
    geom_vline(xintercept = 0, colour = "#8C949D", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, yend = row,
                     colour = crude, linewidth = crude)) +
    geom_point(aes(colour = crude), size = 3.4) +
    scale_colour_manual(values = c(`FALSE` = col_main, `TRUE` = col_crude), guide = "none") +
    scale_linewidth_manual(values = c(`FALSE` = 1.1, `TRUE` = 0.7), guide = "none") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(xlim = rng) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          axis.line.x = element_line(colour = txt, linewidth = 0.4),
          axis.ticks.x = element_line(colour = txt, linewidth = 0.4),
          axis.title.x = element_text(colour = "#5E6976", size = 11),
          plot.margin = margin(5, 5, 5, 5))

  # right: estimate and N text
  p_txt <- ggplot(plot_df, aes(y = row)) +
    geom_text(aes(x = 0, label = est_txt), hjust = 0, size = 4.2, colour = txt) +
    geom_text(aes(x = 0.58, label = n_txt), hjust = 0, size = 4.2, colour = txt) +
    annotate("text", x = 0, y = max(plot_df$row) + 0.7, label = est_lab,
             hjust = 0, size = 4, colour = "#5E6976") +
    annotate("text", x = 0.58, y = max(plot_df$row) + 0.7, label = "Persons / obs.",
             hjust = 0, size = 4, colour = "#5E6976") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    base_theme

  p_lab + p_fp + p_txt + plot_layout(widths = c(1.3, 1.6, 1.6))
}

# ---- 3. Use it --------------------------------------------------------------
# rows_lmm <- tibble::tribble(
#   ~model,                ~label,          ~group,
#   "res_3mse_m0",         "Crude",         "Overall",
#   "res_3mse_m1",         "Model 1",       "Overall",
#   "res_3mse_m2",         "Model 2",       "Overall",
#   "res_3mse_m3",         "Model 3",       "Overall",
#   "res_3mse_m1_apoe4T",  "ε4 carriers",  "Stratified by APOE ε4 (Model 1)",
#   "res_3mse_m1_apoe4F",  "Non-carriers",  "Stratified by APOE ε4 (Model 1)"
# )
#
# models_3mse <- mget(rows_lmm$model)          # grabs the objects by name
#
# # Longitudinal: difference in yearly change
# p_slope <- plot_forest_lmm(models_3mse, rows_lmm, term = "klot_elisa:time",
#                            xlab = "Difference in yearly 3MSE change per SD Klotho (points)")
#
# # Cross-sectional: difference at year 5
# p_level <- plot_forest_lmm(models_3mse, rows_lmm, term = "klot_elisa",
#                            xlab = "Difference in 3MSE at year 5 per SD Klotho (points)")
#
# ggsave(file.path(data_path, "Analysis", "ELISA", "02_LMM", "forest_3mse_slope.png"),
#        p_slope, width = 10, height = 3.8, dpi = 300, bg = "white")
