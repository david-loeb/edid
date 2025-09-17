#' Compute the Ytilde values for group, time and comparison
#'
#' Computes the "Ytilde" function as described in
#' \href{https://arxiv.org/abs/2506.17729}{Chen, Sant'Anna and Xie (2025)} for
#' each combination of treatment group (g), post-treatment period (t) and
#' control group comparison (g', t_pre) set.
#'
#' @param dat The data frame returned by [prep_data()].
#'
#' @returns
#' A list with one element per (g, t, g', t_pre). Each element contains a
#' vector the length of the number of units containing each unit's Ytilde
#' value. The name of each list element contains a sequence of four time
#' period numbers, separated by "`_`". The numbers correspond to the
#' following:
#'
#' 1. Treatment group (g). Note that these are equal to the treatment adoption
#' times in the original dataset minus anticipation periods.
#' 2. Post-treatment period (t).
#' 3. Not-yet-treated comparison group (g').
#' 4. Period 2 for the not-yet-treated comparison group (t_pre in Chen et al.).
#' Note that this value is referred to as `tmid` in the code to avoid confusion
#' with the pre-treatment periods for the treatment group.
#'
#' @export
#'
#' @examples
#' edid_data <- prep_data(
#'   edid_ex_data,
#'   y_var = "outcome",
#'   id_var = "id",
#'   treat_time_var = "treat_adopt_time",
#'   time_var = "time",
#'   num_t_pre = 3
#' )
#'
#' ytilde_list <- get_ytilde(edid_data)
get_ytilde <- function(dat) {
  num_t_pre <- min(dat$g[dat$g > 0]) -
    min(as.numeric(
      stringr::str_sub(names(dat)[grepl('^y', names(dat))], 2, -1)
    ))
  last_t <- max(as.numeric(
    stringr::str_sub(names(dat)[grepl('^y', names(dat))], 2, -1)
  ))
  tx_groups <- sort(unique(dat$g)[unique(dat$g) > 0])
  pi <- list()  # Each treat group's sample proportion
  for (g in tx_groups) pi[[g]] <- mean(dat$g == g)
  pi_0 <- mean(dat$g == 0)
  out <- purrr::map(tx_groups, \(g) {  # Loop thru each tx group
    t1 <- g - num_t_pre  # Time periods & nyt ctrl group setup
    tpost <- g:last_t
    gc <- tx_groups[tx_groups - 1 > t1]
    out <- purrr::map(tpost, \(tpost) {  # Loop thru t_post > nyt ctrls > t_mid
      out <- purrr::map(gc, \(gc) {
        if (gc == g) tmid <- t1:(gc-1) else tmid <- (t1+1):(gc-1)
        out <- purrr::map(tmid, \(tmid) {

          # Ytilde formula; Treated group block
          dat[[stringr::str_c('g', g)]] / pi[[g]] * (  # wt * y_tpost_i - y_t1_i
            dat[[stringr::str_c('y', tpost)]] - dat[[stringr::str_c('y', t1)]]
            - mean(  # - E[y_tpost - y_tmid | G = 0]
              dat[[stringr::str_c('y', tpost)]][dat$g == 0]
              - dat[[stringr::str_c('y', tmid)]][dat$g == 0]
            ) - mean(  # - E[y_tmid - y_t1 | G = g']
              dat[[stringr::str_c('y', tmid)]][dat$g == gc]
              - dat[[stringr::str_c('y', t1)]][dat$g == gc]
            )  # Never treated block
          ) - (pi[[g]] / pi_0) * (dat$g0 / pi[[g]]) * (  # wt * y_post_i - y_mid_i
            dat[[stringr::str_c('y', tpost)]] - dat[[stringr::str_c('y', tmid)]]
            - mean(  # - E[y_tpost - y_tmid | G = 0]
              dat[[stringr::str_c('y', tpost)]][dat$g == 0]
              - dat[[stringr::str_c('y', tmid)]][dat$g == 0]
            )  # Not yet treated block; wt * y_post_i - y_mid_i
          ) - (pi[[g]] / pi[[gc]]) * (dat[[stringr::str_c('g', gc)]] / pi[[g]]) * (
            dat[[stringr::str_c('y', tmid)]] - dat[[stringr::str_c('y', t1)]]
            - mean(  #  - E[y_tmid - y_t1 | G = g']
              dat[[stringr::str_c('y', tmid)]][dat$g == gc]
              - dat[[stringr::str_c('y', t1)]][dat$g == gc]
            )
          )
        })
        names(out) <- tmid
        out
      })
      names(out) <- gc
      out
    })
    names(out) <- tpost
    out
  })
  names(out) <- tx_groups
  purrr::list_flatten(out) |> purrr::list_flatten() |> purrr::list_flatten()
}
