#' Extract the efficient ATT(g,t) results from an EDiD model
#'
#' Extracts the efficient ATT(g,t)s, standard errors, and simultaneous
#' confidence bands from an efficient DiD model and puts them in a data frame.
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
#'
#' @returns
#' A data frame with the number of rows equal to the number of efficient
#' ATT(g,t)s estimated. The data frame with contain the following 19
#' columns:
#' * `g`: The treatment group, defined by the treatment adoption time in the
#' original dataset.
#' * `t`: The time period.
#' * `e`: The "event time," i.e. time relative to `g`.
#' * `att`: The efficient ATT(g,t) estimate.
#' * `se_boot` & `se_analytic`: The bootstrapped and analytic standard errors.
#' * `crit_val_95` & `crit_val_90`: The 95% and 90% simultaneous confidence
#' band critical values obtained during the bootstrapping procedure.
#' * `ci_*_*_*`: The lower (upper) 95% (90%) bootstrapped (anlaytic)
#' simultaneous confidence bands.
#' * `pi_g`: The proportion of the sample comprised of members of treatment
#' group `g`.
#' * `q_ge`: The proportion of all units treated in event time `e` that belong
#' to treatment group `g`.
#' * `q_gt`: The proportion of all units treated in time `t` that belong to
#' treatment group `g`. Units with anticipation periods in time `t` are
#' omitted from this calculation and have `NA` values for this variable.
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
#' res_attgt <- get_edid_results_attgt(
#'   edid_data, edid_mod, anticip = 0
#' )
#'
#' # With an anticipation period
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   id_var = "id",
#'   treat_time_var = "treat_adopt_time",
#'   time_var = "time",
#'   num_t_pre = 2,
#'   anticip = 1
#' )
#' ytilde_list <- get_ytilde(edid_data)
#' inf_func_list <- get_influence_func(edid_data, ytilde_list)
#' edid_mod <- estimate_edid_model(
#'   edid_data, ytilde_list, inf_func_list
#' )
#'
#' res_attgt <- get_edid_results_attgt(
#'   edid_data, edid_mod, anticip = 1
#' )
get_edid_results_attgt <- function(dat, mod, anticip) {
  out <- purrr::map(mod, `[[`, 4) |> unlist()  # Put model list into data frame
  out <- data.frame(att = out) |> tibble::rownames_to_column('g_t')
  other_cols <- purrr::map(mod, \(gt) {
    gt <- gt[5:16]  # remove individual attgts, wts, & inf func
    purrr::map(gt, as.data.frame) |> purrr::list_cbind()
  }) |>
    purrr::list_rbind() |>
    tidyr::unpack(everything(), names_sep = ".") |>
    dplyr::rename_with(~ stringr::str_remove(.x, '[.].+'))
  out <- dplyr::bind_cols(out, other_cols)

  pi_G <- dplyr::summarise(dat, pi = dplyr::n() / nrow(dat), .by = g)
  out |>
    dplyr::mutate(
      g = as.numeric(stringr::str_extract(g_t, '^\\d+')),
      t = as.numeric(stringr::str_extract(g_t, '_(\\d+)', group = 1)),
      e = t - g,
      pi_g = pi_G$pi[match(g, pi_G$g)]  # tx group sample proportions
    ) |>
    dplyr::mutate(q_ge = pi_g / sum(pi_g), .by = e) |>  # g event time shares
    dplyr::mutate(
      g = g + anticip,
      e = e - anticip,
      q_gt = ifelse(t >= g, pi_g / sum(pi_g[t >= g]), NA),  # g cal time shares
      .by = t
    ) |>
    dplyr::relocate(e) |> dplyr::relocate(t) |> dplyr::relocate(g) |>
    dplyr::select(-g_t)
}
