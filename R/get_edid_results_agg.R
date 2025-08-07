#' Get aggregated efficient DiD event study or calendar time results
#'
#' Aggregates the efficient ATT(g,t)s into event study or calendar time results.
#'
#' @param dat The data frame returned by [prep_edid_data()].
#' @param mod The EDiD model returned by [estimate_edid_model()].
#' @param anticip The number of treatment anticipation periods. This must align
#'   with the number specified in [prep_edid_data()]. The treatment group IDs
#'   will be accordingly re-aligned with the original dataset and aggregated
#'   treatment effects computed with [get_edid_results_agg()] will also be
#'   computed accordingly. Unlike in [prep_edid_data()], no default is
#'   provided as a safety measure against the scenario where one forgets to
#'   re-specify non-zero anticipation periods.
#' @param type A string indicating the type of aggregation, either `es` for
#'   event study (default) or `cal` for calendar time.
#' @param res Optionally supply the results data frame returned by
#'   [edid_get_results_attgt()]. The default `NULL` will run
#'   [edid_get_results_attgt()] under the hood to obtain the data frame.
#' @param cluster Boolean specifying whether to cluster standard errors at a
#'   level higher than the unit level. The default is set to `FALSE`, which
#'   results in standard errors clustered at the unit level. A cluster variable
#'   must have been specified in [prep_edid_data()] or the function will throw
#'   an error.
#' @param biters An integer specifying number of bootstrap iterations to
#'   perform.
#' @param seed An integer used to set the random seed for bootstrapping.
#'
#' @returns
#' A data frame with the number of rows equal to the number of time
#' periods for which treatment effects are aggregated. The data frame will
#' contain the following 14 columns:
#'
#' * `e` or `t`: The event (`e`) or calendar (`t`) time period.
#' * `att`: The efficient ATT(g,e) or ATT(g,t) aggregated across all treatment
#' groups.
#' * `se_boot` & `se_analytic`: The bootstrapped and analytic standard errors.
#' * `crit_val_95` & `crit_val_90`: The 95% and 90% simultaneous confidence
#' band critical values obtained during the bootstrapping procedure.
#' * `ci_*_*_*`: The lower (upper) 95% (90%) bootstrapped (anlaytic)
#' simultaneous confidence bands.
#'
#' @export
#'
#' @examples
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   id_var = "id",
#'   treat_time_var = "treat_adopt_time",
#'   time_var = "time",
#'   num_t_pre = 3
#' )
#' ytilde_list <- get_ytilde(edid_data)
#' inf_func_list <- get_influence_func(edid_data, ytilde_list)
#' edid_mod <- estimate_edid_model(
#'   edid_data, ytilde_list, inf_func_list
#' )
#'
#' # Obtain event study results directly after model estimation
#' res_es <- get_edid_results_agg(
#'   edid_data, edid_mod, anticip = 0
#' )
#'
#' # Use individual ATT(g,t) results
#' res_attgt <- get_edid_results_attgt(
#'   edid_data, edid_mod, anticip = 0
#' )
#' res_es <- get_edid_results_agg(
#'   edid_data, edid_mod, anticip = 0, res = res_attgt
#' )
#'
#' # Get calendar time aggregated results
#' res_cal <- get_edid_results_agg(
#'   edid_data, edid_mod, anticip = 0, type = 'cal'
#' )
get_edid_results_agg <- function(dat,
                                 mod,
                                 anticip,
                                 type = 'es',  # 'cal' other option
                                 res = NULL,
                                 cluster = FALSE,
                                 biters = 1000,
                                 seed = 6) {
  if (cluster & !('clust' %in% names(dat))) stop(
    'data must include cluster ID col named "clust"; specify in prep_edid_data()'
  )
  if (!(type %in% c('es', 'cal'))) stop('"type =" must be "es" or "cal"')
  if (is.null(res)) res <- get_edid_results_attgt(dat, mod, anticip)
  eif <- purrr::map(mod, `[[`, 3)
  if (type == 'cal') {  # Drop anticipation times
    not_anticip <- which(res$t >= res$g)
    res <- res[not_anticip, ]
    eif <- eif[not_anticip]
  }
  if (type == 'es') res <- dplyr::mutate(res, t = e, q_gt = q_ge)  # e -> t

  out <- dplyr::summarise(res, att = sum(att * q_gt), .by = t)  # Get ATT(t)'s

  eif <- purrr::map(min(res$t):max(res$t), \(t) {  # Get ATT(t) EIFs
    eif_gt <- eif[res$t == t]  # EIFs for all g in t
    att_gt <- res$att[res$t == t]  # ATT for all g in t
    g_in_t <- res$g[res$t == t]  # All g in t
    pi_g <- res$pi_g[res$t == t]  # Sample pcts for all g in t
    pi_gt_tot <- sum(pi_g)  # Sum of sample pcts for g in t
    q_gt <- pi_g / pi_gt_tot  # Time aggregation weights for g in t

    wt_infunc_1 <- sapply(g_in_t, \(g) {  # Inf func for time agg weights
      pi_g <- unique(res$pi_g[res$g == g])
      (BMisc::TorF(dat$g == g) - pi_g) / pi_gt_tot
    })
    wt_infunc_2 <- rowSums(sapply(g_in_t, \(g) {
      pi_g <- unique(res$pi_g[res$g == g])
      BMisc::TorF(dat$g == g) - pi_g
    })) %*% t(res$pi_g[res$t == t] / pi_gt_tot^2)
    wt_infunc <- wt_infunc_1 - wt_infunc_2

    q_gt_mat <- as.matrix(q_gt)  # Compute final weighted influence func for t
    eif_gt_mat <- matrix(unlist(eif_gt), ncol = length(eif_gt))
    eif_t_att <- eif_gt_mat %*% q_gt_mat
    eif_t_att + wt_infunc %*% as.matrix(att_gt)
  })

  eif <- matrix(unlist(eif), ncol = length(eif))  # Bootstrap SEs
  n <- nrow(eif)
  if (cluster) {
    n <- length(unique(dat$clust))
    cluster <- unique(dat[, c('id', 'clust')])[ ,2]
    n_per_cluster <- aggregate(cluster, by = cluster, length)[ ,2]
    eif <- rowsum(eif, cluster$clust, reorder = T) / n_per_cluster
  }
  boots <- sqrt(n) * withr::with_seed(
    seed, BMisc::multiplier_bootstrap(eif, biters)
  )
  ndg.dim <- (!is.na(colSums(boots))) &
    (colSums(boots^2) > sqrt(.Machine$double.eps)*10)
  bSigma <- apply(boots, 2, \(b) {
    (quantile(b, .75, type = 1, na.rm = T)
     - quantile(b, .25, type = 1, na.rm = T)) /
      (qnorm(.75) - qnorm(.25))
  })
  bT <- suppressWarnings(apply(boots, 1, \(b) {
    max(abs(b / bSigma), na.rm = T)
  }))
  bT <- bT[is.finite(bT)]
  crit_val_95 <- quantile(bT, .95, type = 1, na.rm = T)[[1]]
  crit_val_90 <- quantile(bT, .90, type = 1, na.rm = T)[[1]]
  se_boot <- rep(NA, length(ndg.dim))
  se_boot[ndg.dim] <- as.numeric(bSigma) / sqrt(n)
  se_analytic <- c()  # Analytic SEs
  for (i in 1:ncol(eif)) se_analytic[i] <- sqrt(mean(eif[ ,i]^2) / n)

  out <- out |> dplyr::bind_cols(  # Output results data frame
    se_boot = se_boot, se_analytic = se_analytic,
    crit_val_95 = crit_val_95, crit_val_90 = crit_val_90
  ) |>
    dplyr::mutate(
      ci_low_95_boot = att - se_boot * crit_val_95,
      ci_up_95_boot = att + se_boot * crit_val_95,
      ci_low_90_boot = att - se_boot * crit_val_90,
      ci_up_90_boot = att + se_boot * crit_val_90,
      ci_low_95_analytic = att - se_analytic * crit_val_95,
      ci_up_95_analytic = att + se_analytic * crit_val_95,
      ci_low_90_analytic = att - se_analytic * crit_val_90,
      ci_up_90_analytic = att + se_analytic * crit_val_90
    )
  if (type == 'es') out <- dplyr::rename(out, e = t)
  out
}
