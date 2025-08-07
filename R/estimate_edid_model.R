#' Estimate the efficient DiD model
#'
#' Computes each ATT(g, t, g', t_pre), their efficiency weights, and combines
#' them into the final ATT(g,t)s as described in
#' \href{https://arxiv.org/abs/2506.17729}{Chen, Sant'Anna and Xie (2025)}.
#' Also computes standard errors, both analytic and bootstrapped, and 95% and
#' 90% simultaneous confidence bands.
#'
#' @details
#' The bootstrapped standard errors and simultaneous confidence bands are
#' computed following the procedure described in
#' \href{https://www.sciencedirect.com/science/article/abs/pii/S0304407620303948}{Callaway and Sant'Anna (2021)}
#' with code adapted directly from the authors'
#' \href{https://bcallaway11.github.io/did/}{did package} that implements the
#' procedure outlined in the paper.
#'
#' @param dat The data frame returned by [prep_edid_data()].
#' @param ytilde The list of Ytilde values returned by [get_ytilde()].
#' @param inf_func The list of influence function values returned by
#'   [get_influence_func()].
#' @param cluster Optional boolean to override the default standard error
#'   clustering behavior. The default `NULL` will cluster standard errors
#'   if the data frame has a column named "`clust`", i.e. if a `cluster_var`
#'   was specified in [prep_edid_data()]. This can be overridden by setting
#'   this argument to `FALSE`, producing the default standard errors
#'   (which are still clustered at the unit level). Setting the argument to
#'   `TRUE` without having specified a cluster variable in [prep_edid_data()]
#'   will throw an error.
#' @param biters An integer specifying number of bootstrap iterations to
#'   perform.
#' @param seed An integer used to set the random seed for bootstrapping.
#'
#' @returns
#' A list with one element per treatment group and post-period (g, t). Each
#' element's name is a sequence of two time periods, the first denoting the
#' treatment group (g) and the second denoting the post-period (t). Note that
#' the treatment group numbers are equal to the treatment adoption times in the
#' original dataset minus anticipation periods. Each list element contains the
#' following sixteen items:
#'
#' * 1: The set of ATTs for the ATT(g,t), one for each combination of
#' not-yet-treated comparison groups (g') and their period 2's (t_pre).
#' * 2: The efficiency weights for each ATT(g, t, g', t_pre) used to combine
#' the ATTs into the efficient ATT(g,t).
#' * 3: The efficient influence function values for each unit.
#' * 4: The efficient ATT(g,t).
#' * 5-6: The bootstrapped and analytic standard errors respectively.
#' * 7-8: The 95% and 90% simultaneous confidence band critical values
#' respectively.
#' * 9-16: The 95% and 90% simultaneous confidence bands for bootstrapped and
#' analytic standard errors respectively.
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
#'
#' edid_mod <- estimate_edid_model(
#'   edid_data, ytilde_list, inf_func_list
#' )
#'
#' # Adjust bootstrap iterations & random seed generation
#' edid_mod <- estimate_edid_model(
#'   edid_data,
#'   ytilde_list,
#'   inf_func_list,
#'   biters = 5000,
#'   seed = 66
#' )
#'
#' # With clustered standard errors
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   id_var = "id",
#'   treat_time_var = "treat_adopt_time",
#'   time_var = "time",
#'   num_t_pre = 3,
#'   cluster_var = "cluster"
#' )
#' ytilde_list <- get_ytilde(edid_data)
#' inf_func_list <- get_influence_func(edid_data, ytilde_list)
#'
#' edid_mod <- estimate_edid_model(
#'   edid_data,
#'   ytilde_list,
#'   inf_func_list,
#'   cluster = TRUE
#' )
estimate_edid_model <- function(dat,
                                ytilde,
                                inf_func,
                                cluster = NULL,
                                biters = 1000,
                                seed = 6) {
  if (is.null(cluster) & 'clust' %in% names(dat)) {
    cluster <- TRUE
  } else {
    cluster <- FALSE
  }
  if (cluster & !('clust' %in% names(dat))) stop(
    'data must include cluster ID col named "clust"; specify in prep_edid_data()'
  )
  tx_groups <- unique(stringr::str_extract(names(ytilde), '\\d+'))
  out <- purrr::map(tx_groups, \(g) {
    tpost <- unique(stringr::str_extract(
      names(ytilde), stringr::str_c('^', g, '_(\\d+)'), group = 1
    ))
    tpost <- tpost[!is.na(tpost)]
    out <- purrr::map(tpost, \(tpost) {
      out <- list()  # Subset the proper Ytilde and influence functions
      ytilde <- ytilde[grepl(stringr::str_c('^', g, '_', tpost), names(ytilde))]
      inf_func <- inf_func[grepl(stringr::str_c('^', g, '_', tpost), names(inf_func))]
      inf_func <- matrix(unlist(inf_func), nrow = length(inf_func), byrow = T)
      ytilde <- purrr::map(ytilde, mean)  # Get ATTs for each yr & comparison

      out[['att']] <- matrix(unlist(ytilde))
      V <- cov(t(inf_func))  # Compute efficiency weights
      one <- matrix(rep(1, nrow(V)), nrow = nrow(V))
      out[['wt']] <- (t(one) %*% solve(V)) / as.numeric(
        t(one) %*% solve(V) %*% one
      )
      out[['eif']] <- out[['wt']] %*% inf_func  # Efficient influence function
      out[['eatt']] <- out[['wt']] %*% out[['att']]  # Efficient ATT
      out
    })
    names(out) <- tpost
    out
  })
  names(out) <- tx_groups
  out <- purrr::list_flatten(out)  # one element per ATTgt

  eif <- purrr::map(out, `[[`, 3)  # Bootstrap SEs & get simult. conf. band
  eif <- matrix(unlist(eif), ncol = length(eif))
  n <- nrow(eif)
  if (cluster) {
    n <- length(unique(dat$clust))
    cluster <- unique(dat[ ,c('id', 'clust')])[ ,2]
    n_per_cluster <- aggregate(cluster, by = cluster, length)[ ,2]
    eif <- rowsum(eif, cluster$clust, reorder = T) / n_per_cluster
  }
  boots <- sqrt(n) * withr::with_seed(
    seed, BMisc::multiplier_bootstrap(eif, biters)
  )
  ndg.dim <- (!is.na(colSums(boots))) &  # safety measure
    (colSums(boots^2) > sqrt(.Machine$double.eps)*10)
  bSigma <- apply(boots, 2, \(b) {  # estimate main diag of sqrt(Sigma)
    (quantile(b, .75, type = 1, na.rm = T)
     - quantile(b, .25, type = 1, na.rm = T)) /
      (qnorm(.75) - qnorm(.25))
  })
  bT <- suppressWarnings(apply(boots, 1, \(b) {  # simult. conf. band crit val
    max(abs(b / bSigma), na.rm = T)
  }))
  bT <- bT[is.finite(bT)]
  crit_val_95 <- quantile(bT, .95, type = 1, na.rm = T)[[1]]
  crit_val_90 <- quantile(bT, .90, type = 1, na.rm = T)[[1]]
  se <- rep(NA, length(ndg.dim))
  se[ndg.dim] <- as.numeric(bSigma) / sqrt(n)  # get final SEs

  for (i in 1:length(out)) {
    out[[i]][['se_boot']] <- se[i]
    out[[i]][['se_analytic']] <- sd(eif[ ,i]) / sqrt(n)  # Analytic SE
    out[[i]][['crit_val_95']] <- crit_val_95  # Compute CIs
    out[[i]][['crit_val_90']] <- crit_val_90
    out[[i]][['ci_low_95_boot']] <- out[[i]][['eatt']] - crit_val_95 * se[i]
    out[[i]][['ci_up_95_boot']] <- out[[i]][['eatt']] + crit_val_95 * se[i]
    out[[i]][['ci_low_90_boot']] <- out[[i]][['eatt']] - crit_val_90 * se[i]
    out[[i]][['ci_up_90_boot']] <- out[[i]][['eatt']] + crit_val_90 * se[i]
    out[[i]][['ci_low_95_analytic']] <- out[[i]][['eatt']] - crit_val_95 *
      out[[i]][['se_analytic']]
    out[[i]][['ci_up_95_analytic']] <- out[[i]][['eatt']] + crit_val_95 *
      out[[i]][['se_analytic']]
    out[[i]][['ci_low_90_analytic']] <- out[[i]][['eatt']] - crit_val_90 *
      out[[i]][['se_analytic']]
    out[[i]][['ci_up_90_analytic']] <- out[[i]][['eatt']] + crit_val_90 *
      out[[i]][['se_analytic']]
    out[[i]]['boot'] <- NULL
  }
  out
}
