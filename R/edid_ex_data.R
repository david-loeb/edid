#' EDiD Example Data
#'
#' A simulated dataset to demonstrate the edid functions.
#'
#' @format ## `edid_ex_data`
#' A long-format data frame with 4,000 rows and 6 columns:
#' \describe{
#'   \item{id}{Unit ID, with 500 units.}
#'   \item{cluster}{Cluster ID, with 10 units per cluster.}
#'   \item{time}{Time variable, with 8 time periods.}
#'   \item{treated}{Binary treatment indicator that = 1 when a unit is treated. This variable is not actually used, but is there to demonstrate its relationship to the treatment group identifier variable that the function uses.}
#'   \item{treat_adopt_time}{Time period when a unit first became treated. This can be thought of as the treatment group identifier variable.}
#'   \item{outcome}{Outcome variable for which treatment effects are estimated.}
#' }
#' @source Simulated by package author
'edid_ex_data'
