
# edid

The edid package implements the efficient difference-in-differences
estimator developed by [Chen, Sant’Anna, and Xie
(2025)](https://arxiv.org/abs/2506.17729) in “Efficient
Difference-in-Differences and Event Study Estimators.” I wrote the code
for a personal research project, as I could not find a publicly
available implementation, so I thought I’d share it in case others would
like to use the estimator as well. The functionality is quite limited; I
describe the limitations below. I imagine the paper authors will release
a much more robust implementation at some point. However, there does not
currently appear to be another publicly available implementation, so I
hope this is better than nothing!

## Installation

Install the package from this GitHub page by running either the {pak} or
{devtools} code below in R.

``` r
# install.packages("pak")
pak::pak("david-loeb/edid")

# install.packages("devtools")
devtools::install_github("david-loeb/edid")
```

## Overview

The package contains a series of functions designed to be run
sequentially, beginning with data setup and ending with a data frame of
results. The function `edid()` executes the entire process by running
each of these functions under the hood. This is the most straightforward
way to use the package. The code below demonstrates how to use `edid()`.

``` r
library(edid)

edid_results <- edid(
  edid_ex_data,
  y_var = "outcome",
  treat_time_var = "treat_adopt_time",
  id_var = "id",
  time_var = "time",
  num_t_pre = 3
)
```

The result will be a data frame with the efficient $ATT(g,t)$, event
study, and calendar time aggregation results. You can also request a
more complete set of model components, including the individual
$ATT(g,t)$ and weights used to compute the efficient $ATT(g,t)$ by
setting `get_full_mod = TRUE` in the `edid()` function.

## Full set of package functions

A list of all the package functions with a brief description of each is
below. See the package documentation for a detailed explanation of each
function. The functions are run in the order listed to execute the
modeling process.

- `prep_edid_data()`: Takes a panel dataset and converts it to the
  proper format for the modeling functions.
- `get_ytilde()`: Returns the $\tilde{Y}_{g',t\_{pre}}^{att(g,t)}$
  values for each unit.
- `get_influence_func()`: Returns the influence function values for each
  unit.
- `estimate_edid_model()`: Estimates the efficient $ATT(g,t)$, standard
  errors, and confidence intervals.
- `get_edid_results_attgt()`: Extracts the efficient $ATT(g,t)$ and
  related results into a data frame.
- `get_edid_results_agg()`: Aggregates the results into event study
  and/or calendar time effects and returns the results in a data frame.
- `edid()`: Runs all of the above functions and returns any or all of
  the results of the three previous functions.

## Limitations

edid’s biggest limitation is that it does not handle covariates. It also
requires that all treatment groups use the same number of pre-treatment
periods to estimate treatment effects, and it requires balanced panels.
It is probably limited in other ways that aren’t coming to mind. I am
definitely open to expanding the functionality so don’t hesitate to
reach out with requests / issues / feedback of any kind.

## Acknowledgements

First I would like to thank the authors Drs. Xiaohong Chen, Pedro H. C.
Sant’Anna, and Haitian Xie for their fantastic work developing this
estimator. It has already proved to be quite useful in my own work. I
highly recommend reading [their
paper](https://arxiv.org/abs/2506.17729). I also want to thank
Dr. Brantly Callaway (and Pedro Sant’Anna again) for the multiplier
bootstrap procedure developed in their renowned [2021
paper](https://www.sciencedirect.com/science/article/abs/pii/S0304407620303948)
“Difference-in-Differences With Multiple Time Periods.” I used not only
the procedure but much of the code directly from their implementation in
the [{did} package](https://bcallaway11.github.io/did/).
