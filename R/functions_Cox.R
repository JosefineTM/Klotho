########################################
#   Cox Regression Utility Functions   #
########################################
#
# - runcox():          one Cox model per exposure, all terms + PH check + metadata
# - run_models():      several covariate sets (crude, M1, M2, ...) on the same sample
# - covar_change():    change-in-estimate: what each covariate does to the exposure HR
# - cox_gt():          gt table for one model
# - cox_gt_models():   side-by-side gt table for several models
# - plot_models_forest(), plot_covariates(), plot_covar_change(): figures
# - volcano_plot():    volcano plot (high-throughput use)
# - runcox_int():      interaction models
#
# Requirements: data has columns ID, tto<outcome>, <outcome>
########################################

library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(gt)


# ---- helper: tidy a coxph fit with robust SEs --------------------------------
.tidy_cox <- function(fit) {
  s  <- summary(fit)
  co <- s$coefficients
  ci <- s$conf.int
  se_col <- if ("robust se" %in% colnames(co)) "robust se" else "se(coef)"
  tibble(
    term    = rownames(co),
    coef    = co[, "coef"],
    HR      = co[, "exp(coef)"],
    SE      = co[, se_col],
    CI_low  = ci[, "lower .95"],
    CI_high = ci[, "upper .95"],
    pval    = co[, "Pr(>|z|)"]
  ) %>%
    mutate(HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_low, CI_high))
}


# ---- Run Cox regression for each exposure ------------------------------------
# covar:         character vector of covariates (NULL = crude)
# complete_vars: if given, rows with NA in any of these are dropped first
#                (use to keep the same sample across models)
runcox <- function(covar = NULL, data, vars, outcome, complete_vars = NULL, ...) {
  
  if (!"ID" %in% names(data)) stop("data needs an ID column (e.g. rename(ID = idno))")
  
  if (!is.null(complete_vars)) {
    data <- data %>% drop_na(all_of(complete_vars))
  }
  
  results <- list()
  terms   <- list()
  meta    <- list()
  
  for (j in vars) {
    
    modelForm <- as.formula(
      paste0("Surv(tto", outcome, ", ", outcome, ") ~ ",
             paste(c(j, covar), collapse = " + "))
    )
    
    fit_overall <- tryCatch(
      coxph(modelForm, data = data, id = ID, robust = TRUE),
      error = function(e) { message("Model failed for ", j, ": ", e$message); NULL }
    )
    if (is.null(fit_overall)) next
    
    td <- .tidy_cox(fit_overall)
    
    # all terms (exposure + covariates)
    terms[[j]] <- td %>%
      mutate(Name = j,
             is_exposure = term %in% c(j, paste0("`", j, "`")),
             .before = 1)
    
    # exposure row only
    results[[j]] <- td %>%
      filter(term %in% c(j, paste0("`", j, "`"))) %>%
      mutate(Name = j, n = fit_overall$n, events = fit_overall$nevent, .before = 1) %>%
      select(-term)
    
    # proportional hazards check
    zph <- tryCatch(cox.zph(fit_overall)$table, error = function(e) NULL)
    ph_exp  <- if (!is.null(zph) && j %in% rownames(zph)) zph[j, "p"] else NA
    ph_glob <- if (!is.null(zph)) zph["GLOBAL", "p"] else NA
    
    meta[[j]] <- tibble(
      variable     = j,
      formula      = paste(deparse(modelForm), collapse = ""),
      n_used       = fit_overall$n,
      n_events     = fit_overall$nevent,
      n_missing    = nrow(data) - fit_overall$n,
      method       = fit_overall$method,
      ph_p_exposure = ph_exp,
      ph_p_global   = ph_glob
    )
  }
  
  results <- bind_rows(results) %>%
    mutate(logHR = log(HR),
           padj  = p.adjust(pval, method = "fdr")) %>%
    arrange(pval)
  
  list(
    results = results,
    terms   = bind_rows(terms),
    meta    = meta,              # kept for backward compatibility
    meta_df = bind_rows(meta)
  )
}


# ---- Run several covariate sets ----------------------------------------------
# models: named list, e.g. list("Crude" = NULL, "Model 1" = cov_m1, "Model 2" = cov_m2)
# same_sample = TRUE: every model uses people complete on ALL covariates in all models
run_models <- function(models, data, vars, outcome, same_sample = TRUE, ...) {
  
  complete_vars <- if (same_sample) {
    unique(c(unlist(models), paste0("tto", outcome), outcome))
  } else NULL
  
  fits <- imap(models, ~ runcox(covar = .x, data = data, vars = vars,
                                outcome = outcome, complete_vars = complete_vars, ...))
  
  lvl <- names(models)
  list(
    results = imap_dfr(fits, ~ mutate(.x$results, model = .y, .before = 1)) %>%
      mutate(model = factor(model, levels = lvl)),
    terms   = imap_dfr(fits, ~ mutate(.x$terms,   model = .y, .before = 1)) %>%
      mutate(model = factor(model, levels = lvl)),
    meta    = imap_dfr(fits, ~ mutate(.x$meta_df, model = .y, .before = 1)) %>%
      mutate(model = factor(model, levels = lvl)),
    fits    = fits
  )
}


# ---- Change-in-estimate: what does each covariate do? ------------------------
# direction = "add":  base_covar + one covariate at a time
# direction = "drop": full (base_covar + add_covar) minus one covariate at a time
# pct_change_HR = % change in the exposure HR relative to the reference model
covar_change <- function(data, vars, outcome, base_covar, add_covar,
                         direction = c("add", "drop"), same_sample = TRUE) {
  
  direction <- match.arg(direction)
  full <- unique(c(base_covar, add_covar))
  complete_vars <- if (same_sample) unique(c(full, paste0("tto", outcome), outcome)) else NULL
  
  ref_covar <- if (direction == "add") base_covar else full
  ref <- runcox(ref_covar, data, vars, outcome, complete_vars)$results %>%
    select(Name, HR_ref = HR)
  
  map_dfr(add_covar, function(cv) {
    cov_now <- if (direction == "add") c(base_covar, cv) else setdiff(full, cv)
    runcox(cov_now, data, vars, outcome, complete_vars)$results %>%
      transmute(Name, covariate = cv, HR, CI_low, CI_high, pval)
  }) %>%
    left_join(ref, by = "Name") %>%
    mutate(pct_change_HR = 100 * (HR / HR_ref - 1),
           direction = direction)
}


# ---- gt table: one model -----------------------------------------------------
cox_gt <- function(model, model_name, covariates, p_cutoff = 0.05) {
  
  model %>%
    select(Name, HR, HR_CI, pval, padj, any_of(c("n", "events"))) %>%
    filter(pval <= p_cutoff) %>%
    separate(Name, into = c("Protein", "Platform"), sep = "_", extra = "merge") %>%
    gt() %>%
    tab_header(title = paste0("Cox regression results - ", model_name)) %>%
    fmt_number(columns = HR, decimals = 2) %>%
    fmt(columns = c(pval, padj), fns = function(x) format.pval(x, digits = 2, eps = 0.001)) %>%
    cols_label(HR = "Hazard ratio", HR_CI = "HR (95% CI)",
               pval = "p", padj = "FDR-adjusted p") %>%
    tab_options(table.font.names = "Arial",
                heading.title.font.size = 14, table.font.size = 12) %>%
    tab_footnote(if (length(covariates) == 0) "Unadjusted" else
      paste0("Adjusted for ", paste(covariates, collapse = ", ")))
}


# ---- gt table: several models side by side -----------------------------------
cox_gt_models <- function(results, models, title = "Circulating Klotho and incident dementia") {
  
  mods <- levels(results$model)
  
  wide <- results %>%
    select(model, Name, HR_CI, pval) %>%
    pivot_wider(names_from = model, values_from = c(HR_CI, pval),
                names_glue = "{model}__{.value}") %>%
    separate(Name, into = c("Protein", "Platform"), sep = "_", extra = "merge")
  
  col_order <- c("Protein", "Platform",
                 as.vector(rbind(paste0(mods, "__HR_CI"), paste0(mods, "__pval"))))
  
  n_info <- results %>% group_by(model) %>%
    summarise(txt = paste0("N = ", min(n), "–", max(n), ", events = ",
                           min(events), "–", max(events)), .groups = "drop")
  
  footnotes <- map_chr(mods, function(m) {
    cv <- models[[m]]
    paste0(m, ": ", if (length(cv) == 0) "unadjusted" else paste(cv, collapse = ", "),
           " (", n_info$txt[n_info$model == m], ")")
  })
  
  tbl <- wide %>%
    select(all_of(col_order)) %>%
    gt() %>%
    tab_spanner_delim(delim = "__") %>%
    fmt(columns = ends_with("__pval"),
        fns = function(x) format.pval(x, digits = 2, eps = 0.001)) %>%
    cols_label(.list = setNames(
      as.list(rep(c("HR (95% CI)", "p"), length(mods))),
      as.vector(rbind(paste0(mods, "__HR_CI"), paste0(mods, "__pval")))
    )) %>%
    tab_header(title = title, subtitle = "Hazard ratio per 1 SD higher protein") %>%
    tab_options(table.font.names = "Arial", table.font.size = 12)
  
  for (f in footnotes) tbl <- tbl %>% tab_source_note(f)
  tbl
}


# ---- Forest plot: exposure HR across models ----------------------------------
plot_models_forest <- function(results, x_lab = "Hazard ratio per 1 SD (95% CI)") {
  ggplot(results, aes(x = HR, y = Name, xmin = CI_low, xmax = CI_high, colour = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.35) +
    scale_x_log10() +
    scale_colour_brewer(palette = "Set2") +
    labs(x = x_lab, y = NULL, colour = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


# ---- Forest plot: all covariates in one model --------------------------------
plot_covariates <- function(terms, protein, model_name) {
  terms %>%
    filter(Name == protein, model == model_name) %>%
    mutate(term = factor(term, levels = rev(term))) %>%
    ggplot(aes(x = HR, y = term, xmin = CI_low, xmax = CI_high, colour = is_exposure)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(size = 0.3) +
    scale_x_log10() +
    scale_colour_manual(values = c(`TRUE` = "#990033", `FALSE` = "#4C78A8"), guide = "none") +
    labs(x = "Hazard ratio (95% CI)", y = NULL,
         title = paste0(protein, " - ", model_name)) +
    theme_minimal()
}


# ---- Heatmap: change-in-estimate ----------------------------------------------
plot_covar_change <- function(chg) {
  lim <- max(abs(chg$pct_change_HR), na.rm = TRUE)
  ggplot(chg, aes(x = covariate, y = Name, fill = pct_change_HR)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.1f", pct_change_HR)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "#990033",
                         midpoint = 0, limits = c(-lim, lim)) +
    labs(x = NULL, y = NULL, fill = "% change\nin HR",
         title = paste0("Change in exposure HR when covariate is ",
                        ifelse(unique(chg$direction) == "add", "added", "dropped"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ---- Volcano plot (fixed) -----------------------------------------------------
volcano_plot <- function(data, HR_col = "HR", covar, x_lab, log = FALSE,
                         adj_p_values = TRUE, max_overlaps = 25) {
  
  null_val <- if (log) 0 else 1
  p_col <- if (adj_p_values) "padj" else "pval"
  
  plot_data <- data %>%
    mutate(
      x = .data[[HR_col]],
      p = .data[[p_col]],
      label = sub("_.*", "", Name),
      colors = case_when(
        p < 0.05 & x > null_val ~ "Significantly higher risk",
        p < 0.05 & x < null_val ~ "Significantly lower risk",
        TRUE                    ~ "Not significant"
      )
    )
  
  ggplot(plot_data, aes(x = x, y = -log10(p), colour = colors)) +
    geom_point(shape = 20, size = 4) +
    scale_colour_manual(values = c("Significantly higher risk" = "#990033",
                                   "Significantly lower risk"  = "#69b3a2",
                                   "Not significant"           = "darkgray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey20") +
    geom_vline(xintercept = null_val, linetype = "dotted", colour = "grey50") +
    geom_text_repel(data = filter(plot_data, p < 0.05), aes(label = label),
                    colour = "black", size = 3, max.overlaps = max_overlaps) +
    labs(caption = paste0("Covariates: ", paste(covar, collapse = ", ")),
         y = paste0("-log10(", if (adj_p_values) "FDR-adjusted p" else "p", ")"),
         x = paste("Hazard ratio\n", x_lab)) +
    theme_minimal() +
    theme(legend.position = "none")
}


# ---- Interaction models (weights optional) -----------------------------------
runcox_int <- function(covar, data, vars, outcome, int_var, intvar_other = "",
                       weights_col = NULL, ...) {
  
  covar  <- covar[!covar %in% int_var]
  result <- list()
  meta   <- list()
  
  if (!is.null(weights_col)) data$.wt <- data[[weights_col]]
  
  get_coef_by_idx <- function(fit, idx) {
    if (length(idx) == 0) return(list(HR = NA, se = NA, CI_lower = NA, CI_upper = NA, pval = NA))
    i <- idx[1]
    co <- fit$coefficients
    se_col <- if ("robust se" %in% colnames(co)) "robust se" else "se(coef)"
    list(HR = co[i, "exp(coef)"], se = co[i, se_col],
         CI_lower = fit$conf.int[i, "lower .95"], CI_upper = fit$conf.int[i, "upper .95"],
         pval = co[i, "Pr(>|z|)"])
  }
  
  for (j in vars) {
    modelForm <- as.formula(
      paste0("Surv(tto", outcome, ", ", outcome, ") ~ ",
             paste(c(j, int_var, covar, paste0(j, ":", int_var)), collapse = " + "))
    )
    
    fit_overall <- tryCatch({
      if (is.null(weights_col)) coxph(modelForm, data = data, id = ID, robust = TRUE)
      else coxph(modelForm, data = data, weights = .wt, id = ID, robust = TRUE)
    }, error = function(e) { message("Model failed for ", j, ": ", e$message); NULL })
    if (is.null(fit_overall)) next
    
    fit <- summary(fit_overall)
    all_names <- rownames(fit$coefficients)
    
    target_int <- paste0(j, ":", int_var, intvar_other)
    target_var <- paste0(int_var, intvar_other)
    
    int_idx <- which(all_names == target_int)
    if (length(int_idx) == 0) int_idx <- which(all_names == paste0(target_var, ":", j))
    var_idx  <- which(all_names == target_var)
    prot_idx <- which(all_names == j)
    
    interact <- get_coef_by_idx(fit, int_idx)
    intvar   <- get_coef_by_idx(fit, var_idx)
    prot     <- get_coef_by_idx(fit, prot_idx)
    
    result[[j]] <- data.frame(
      Name = j,
      HR_int = interact$HR, se_int = interact$se, CI_lower_int = interact$CI_lower,
      CI_upper_int = interact$CI_upper, pval_int = interact$pval,
      HR_prot = prot$HR, se_prot = prot$se, CI_lower_prot = prot$CI_lower,
      CI_upper_prot = prot$CI_upper, pval_prot = prot$pval,
      HR_intvar = intvar$HR, se_intvar = intvar$se, CI_lower_intvar = intvar$CI_lower,
      CI_upper_intvar = intvar$CI_upper, pval_intvar = intvar$pval
    )
    
    meta[[j]] <- tibble(variable = j, n_used = fit_overall$n, n_events = fit_overall$nevent,
                        formula = paste(deparse(modelForm), collapse = ""))
  }
  
  result_df <- bind_rows(result) %>%
    mutate(padj_int  = p.adjust(pval_int, method = "fdr"),
           HR_CI_int = sprintf("%.2f (%.2f, %.2f)", HR_int, CI_lower_int, CI_upper_int))
  
  list(results = result_df, meta = meta, meta_df = bind_rows(meta))
}