edid_ex_data <- data.frame(
  id = rep(1:500, each = 8),
  cluster = rep(1:50, each = 80),
  time = rep(1:8, 500),
  treated = c(
    rep(0,1200),  # g0 (never treated)
    rep(c(rep(0,3), rep(1,5)), 120),  # g4
    rep(c(rep(0,4), rep(1,4)), 80),   # g5
    rep(c(rep(0,5), rep(1,3)), 100),  # g6
    rep(c(rep(0,6), rep(1,2)), 50)    # g7
  ),
  treat_adopt_time = c(
    rep(0,1000), rep(9,200),  # g0, some treated after last period
    rep(4,960), rep(5,640), rep(6,800), rep(7,400)  # treated groups
  )
) |>
  dplyr::mutate(  # General outcome time trends & dynamic treatment effects
    outcome = 10 + time + treated * time + rnorm(4000)
  )

usethis::use_data(edid_ex_data, overwrite = TRUE)
