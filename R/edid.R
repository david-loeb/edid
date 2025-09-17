#' Run the full efficient DiD modeling process
#'
#' This runs the entire efficient DiD process from the data setup to the
#' extraction of final results. It runs each individual package function under
#' the hood and returns the data frames produced by the
#' [get_results_attgt()] and [get_results_agg()] functions, the full
#' set of model elements produced by [estimate_model()], or both.
#'
#' @details
#' Note that the list of model elements returned by setting
#' `get_full_mod = TRUE` does not include the Ytilde or influence function
#' values for each unit. To obtain those, run [get_ytilde()] or
#' [get_influence_func()] on a data frame prepared with [prep_data()].
#'
#' @param dat A data frame with the data you want to model.
#' @param y A character string with the name of the outcome variable name
#'   in the data frame.
#' @param treat_time A string with the name of a numeric variable in the
#'   data frame that indicates the first period of treatment. Never treated
#'   units should have a value of 0 on this variable.
#' @param id A string with the unit-level ID variable name in the
#'   data frame. The variable can be numeric or character.
#' @param time A string with name of the numeric time variable in the
#'   data frame.
#' @param num_t_pre An integer specifying the number of pre-treatment
#'   periods to use for computing treatment effects.
#' @param cluster An optional string with the name of a variable in the
#'   data frame to cluster standard errors on. `NULL` defaults to clustering at
#'   the unit level. The variable can be numeric or character.
#' @param anticip An integer specifying the number of pre-treatment
#'   anticipation periods. Treatment adoption will be re-defined as the first
#'   anticipation period for modeling purposes. The treatment adoption period
#'   will be re-aligned with the initial dataset when extracting results and
#'   computing aggregated effects. The default is set to `0`.
#' @param get_attgt Boolean specifying whether to return the individual
#'   efficient ATT(g,t) results. Default is `TRUE`.
#' @param get_es Boolean specifying whether to return event study results.
#'   Default is `TRUE`.
#' @param get_cal Boolean specifying whether to return calendar time
#'   aggregation results. Default is `TRUE`.
#' @param get_full_mod Boolean specifying whether to return the list with the
#'   full set of model elements, including the individual ATT(g,t)s and
#'   weights used to compute the efficient ATT(g,t)s. Default is `FALSE`.
#' @param biters An integer specifying number of bootstrap iterations to
#'   perform.
#' @param seed An integer used to set the random seed for bootstrapping.
#'
#' @returns
#' The function will either return a data frame of the efficient DiD ATT
#' results, a list of full model elements, or a list containing both of these
#' pieces of output, with the data frame of results as element one and the
#' model elements list as element two, depending on the combination of
#' requested results.
#'
#' The results data frame will be returned if any of `get_attgt`, `get_es`, or
#' `get_cal` are set to `TRUE`. The data frame will have a number of rows equal
#' to the total number of ATT estimates among the requested results sets. The
#' data frame is the result of binding the rows of the data frames returned by
#' the [get_results_attgt()] and [get_results_agg()] functions; see
#' their documentation for descriptions of the columns. One additional column,
#' `type`, will be added as the first column of the data frame indicating the
#' type of estimate (`attgt`, `es` or `cal`).
#'
#' The list of full model elements is the output of the [estimate_model()]
#' function; see that function's documentation for a description of the output.
#'
#' @export
#'
#' @examples
#' # Return the efficient ATT(g,t)s, event study, and calendar
#' # time aggregated results.
#' edid_res <- edid(
#'   edid_ex_data,
#'   y = "outcome",
#'   treat_time = "treat_adopt_time",
#'   id = "id",
#'   time = "time",
#'   num_t_pre = 3
#' )
#'
#' # Return the event study results only
#' edid_res_all <- edid(
#'   edid_ex_data,
#'   y = "outcome",
#'   treat_time = "treat_adopt_time",
#'   id = "id",
#'   time = "time",
#'   num_t_pre = 3,
#'   get_attgt = FALSE,
#'   get_cal = FALSE
#' )
#'
#' # Return the ATT(g,t) data frame and the full model results
#' edid_res_all <- edid(
#'   edid_ex_data,
#'   y = "outcome",
#'   treat_time = "treat_adopt_time",
#'   id = "id",
#'   time = "time",
#'   num_t_pre = 3,
#'   get_es = FALSE,
#'   get_cal = FALSE,
#'   get_full_mod = TRUE
#' )
#'
#' # Get results with clustered standard errors
#' edid_res_all <- edid(
#'   edid_ex_data,
#'   y = "outcome",
#'   treat_time = "treat_adopt_time",
#'   id = "id",
#'   time = "time",
#'   num_t_pre = 3,
#'   cluster = "cluster"
#' )
#'
#' # Specify an anticipation period
#' edid_res_all <- edid(
#'   edid_ex_data,
#'   y = "outcome",
#'   treat_time = "treat_adopt_time",
#'   id = "id",
#'   time = "time",
#'   num_t_pre = 2,
#'   anticip = 1
#' )
edid <- function(dat,
                 y,
                 treat_time,
                 id,
                 time,
                 num_t_pre,
                 cluster = NULL,
                 anticip = 0,
                 get_attgt = TRUE,
                 get_es = TRUE,
                 get_cal = TRUE,
                 get_full_mod = FALSE,
                 biters = 1000,
                 seed = 6) {
  if (!any(get_attgt, get_es, get_cal, get_full_mod)) stop(
    'must specify at least one set of results to return'
  )
  dat <- prep_data(
    dat, y, treat_time, id, time, num_t_pre, cluster, anticip
  )
  ytilde <- get_ytilde(dat)
  inf_func <- get_influence_func(dat, ytilde)
  if (!is.null(cluster)) cluster <- TRUE else cluster <- FALSE
  mod <- estimate_model(dat, ytilde, inf_func, cluster, biters, seed)
  res <- get_results_attgt(dat, mod, anticip)
  out <- list()
  if (get_attgt) out[['attgt']] <- res |> dplyr::mutate(type = 'attgt')
  if (get_es) {
    out[['es']] <- get_results_agg(
      dat, mod, anticip, 'es', res, cluster, biters, seed
    ) |>
      dplyr::mutate(type = 'es')
  }
  if (get_cal) {
    out[['cal']] <- get_results_agg(
      dat, mod, anticip, 'cal', res, cluster, biters, seed
    ) |>
      dplyr::mutate(type = 'cal')
  }
  if (any(get_attgt, get_es, get_cal)) {
    res_df <- purrr::list_rbind(out) |> dplyr::relocate(type)
    if (!get_attgt & get_es & get_cal) {
      res_df <- dplyr::relocate(res_df, t, .before = e)
    }
  }

  if (get_full_mod) {
    if (any(get_attgt, get_es, get_cal)) {
      out <- list()
      out[['res']] <- res_df
      out[['mod']] <- mod
    } else {
      out <- mod
    }
  } else {
    out <- res_df
  }
  out
}
