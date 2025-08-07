#' Prepare a dataset for EDiD modeling
#'
#' This takes a panel dataset with the data you want to model and returns a
#' wide dataset with the variables needed for modeling, renamed to work with
#' the modeling functions.
#'
#' @param dat A data frame with the data you want to model.
#' @param y_var A character string with the name of the outcome variable name
#'   in the data frame.
#' @param treat_time_var A string with the name of a numeric variable in the
#'   data frame that indicates the first period of treatment. Never treated
#'   units should have a value of 0 on this variable.
#' @param id_var A string with the unit-level ID variable name in the
#'   data frame. The variable can be numeric or character.
#' @param time_var A string with name of the numeric time variable in the
#'   data frame.
#' @param num_t_pre An integer specifying the number of pre-treatment
#'   periods to use for computing treatment effects.
#' @param cluster_var An optional string with the name of a variable in the
#'   data frame to cluster standard errors on. `NULL` defaults to clustering at
#'   the unit level. The variable can be numeric or character.
#' @param anticip An integer specifying the number of pre-treatment
#'   anticipation periods. Treatment adoption will be re-defined as the first
#'   anticipation period for modeling purposes. The treatment adoption period
#'   will be re-aligned with the initial dataset when extracting results and
#'   computing aggregated effects. The default is set to `0`.
#'
#' @returns
#' A wide-format, balanced-panel data frame that can be used with the
#' efficient DiD estimator functions. The row length will be the number of
#' units with data observed in all time periods. The columns will include:
#'
#' * `id`: each unit's ID
#' * `g`: the unit's treatment group
#' * `y[time]`: one column per time period containing outcome data for that
#' period
#' * `g[treatment adoption time]`: one column per treatment group containing a
#' binary indicator of treatment group membership.
#' @export
#'
#' @examples
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   treat_time_var = "treat_adopt_time",
#'   id_var = "id",
#'   time_var = "time",
#'   num_t_pre = 3
#' )
#'
#' # Clustered standard errors
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   treat_time_var = "treat_adopt_time",
#'   id_var = "id",
#'   time_var = "time",
#'   num_t_pre = 3,
#'   cluster_var = "cluster"
#' )
#'
#' # Specify anticipation
#' edid_data <- prep_edid_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   treat_time_var = "treat_adopt_time",
#'   id_var = "id",
#'   time_var = "time",
#'   num_t_pre = 2,
#'   anticip = 1
#' )
prep_edid_data <- function(dat,
                           y_var,
                           treat_time_var,
                           id_var,
                           time_var,
                           num_t_pre,
                           cluster_var = NULL,
                           anticip = 0) {
  dat[[treat_time_var]] <- ifelse(  # Set treat var to first anticipation period
    dat[[treat_time_var]] > 0,
    dat[[treat_time_var]] - anticip, dat[[treat_time_var]]
  )  # Set t1 & t_last
  min_g <- min(dat[[treat_time_var]][dat[[treat_time_var]] > 0], na.rm = T)
  min_data_t <- min(dat[[time_var]], na.rm = T)
  first_t <- min_g - num_t_pre
  if (first_t < min_data_t) stop(
    'specified number of pre-treatment periods extends beyond first period of panel'
  )
  last_t <- max(dat[[time_var]], na.rm = T)
  last_t <- find_last_t(last_t, dat, y_var, time_var)  # so last t != all NA Y
  dat <- dat |>  # Prep data
    dplyr::filter(get(time_var) %in% (first_t):(last_t)) |>
    dplyr::select(
      id = dplyr::all_of(id_var), clust = dplyr::all_of(cluster_var),
      t = dplyr::all_of(time_var), y = dplyr::all_of(y_var),
      g = dplyr::all_of(treat_time_var)
    ) |>
    dplyr::mutate(y = ifelse(is.infinite(y), NA, y)) |>
    dplyr::filter(!is.na(y)) |>  # drop missing Y
    BMisc::makeBalancedPanel('id', 't') |>
    dplyr::mutate(g = ifelse(g > last_t, 0, g)) |>  # treat after last t = never
    tidyr::pivot_wider(names_from = t, values_from = y) |>  # wide format
    dplyr::rename_with(~ stringr::str_replace(.x, '^(\\d)', 'y\\1'))  # prefix y cols
  tx_group_dummies <- purrr::map(sort(unique(dat$g[!is.na(dat$g)])), \(g) {
    stats::setNames(
      data.frame(as.integer(ifelse(dat$g == g, 1, 0))), stringr::str_c('g', g)
    )
  }) |>
    purrr::list_cbind()
  dplyr::bind_cols(dat, tx_group_dummies)
}


#' Find the last time period with observed outcome data
#'
#' This recursive function is used in `prep_edid_data()` to ensure that the
#' last period retained does not have all `NA` data for the outcome.
#'
#' @param t A single number that is initialized as the last time period in
#'   the dataset.
#' @param dat The dataset supplied to `prep_edid_data()`.
#' @param y_var The outcome variable supplied to `prep_edid_data()`.
#' @param time_var The time variable supplied to `prep_edid_data()`.
#'
#' @returns A single number, the last time period in the dataset with
#'   the outcome observed.
#'
#' @noRd
find_last_t <- function(t, dat, y_var, time_var) {
  if (all(is.na(dat[[y_var]][dat[[time_var]] == t]))) {
    find_last_t(t-1, dat, y_var, time_var)
  }  else {
    t
  }
}
