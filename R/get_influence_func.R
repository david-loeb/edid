#' Compute the influence function values for group, time and comparison
#'
#' Computes the influence function as described in
#' \href{https://arxiv.org/abs/2506.17729}{Chen, Sant'Anna and Xie (2025)} for
#' each combination of treatment group (g), post-treatment period (t) and
#' control group comparison (g', t_pre) set.
#'
#' @param dat The data frame returned by [prep_edid_data()].
#' @param ytilde The list of Ytilde values returned by [get_ytilde()].
#'
#' @returns
#' A list with one element per (g, t, g', t_pre). Each element contains a
#' vector the length of the number of units containing each unit's influence
#' function value. The name of each list element contains a sequence of four
#' time period numbers, separated by "`_`". The numbers correspond to the
#' following:
#'
#' 1. Treatment group (g). Note that these are equal to the treatment adoption
#' times in the original dataset minus anticipation periods.
#' 2. Post-treatment period (t)
#' 3. Not-yet-treated comparison group (g')
#' 4. Period 2 for the not-yet-treated comparison group (t_pre in Chen et al.).
#' Note that this value is referred to as `tmid` in the code to avoid confusion
#' with the pre-treatment periods for the treatment group.
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
#'
#' inf_func_list <- get_influence_func(edid_data, ytilde_list)
get_influence_func <- function(dat, ytilde) {
  ids <- names(ytilde)
  min_g <- min(as.numeric(stringr::str_extract(ids, '^\\d+')))
  min_t1 <- min(
    as.numeric(stringr::str_extract(
      ids[as.numeric(stringr::str_extract(ids, '^\\d+')) == min_g], '\\d+$'
    ))
  )
  num_t_pre <- min_g - min_t1
  purrr::imap(ytilde, \(ytilde, id) {
    g <- as.numeric(stringr::str_extract(id, '\\d+'))
    tpost <- as.numeric(stringr::str_extract(id, '_(\\d+)', group = 1))
    gc <- as.numeric(stringr::str_extract(id, '_\\d+_(\\d+)', group = 1))
    tmid <- as.numeric(stringr::str_extract(id, '\\d+$'))
    t1 <- g - num_t_pre
    pi_g <- mean(dat$g == g)

    # Influence function formula
    ytilde - dat[[stringr::str_c('g', g)]] / pi_g * (  # Ytilde - Gg / sample %
      mean(  # x ATT: E[y_tpost - y_t1 | G = g] -
        dat[[stringr::str_c('y', tpost)]][dat$g == g]
        - dat[[stringr::str_c('y', t1)]][dat$g == g]
      ) - (
        mean(  # (E[y_tpost - y_tmid | G = never] +
          dat[[stringr::str_c('y', tpost)]][dat$g == 0]
          - dat[[stringr::str_c('y', tmid)]][dat$g == 0]
        ) + mean(  # E[y_tmid - y_t1 | G = g'])
          dat[[stringr::str_c('y', tmid)]][dat$g == gc]
          - dat[[stringr::str_c('y', t1)]][dat$g == gc]
        )
      )
    )
  })
}
