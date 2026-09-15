# Changelog

## SteadyStateBVAR 0.2.0

CRAN release: 2026-09-15

- [`IRF()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/IRF.md)
  argument `shock` renamed to `impulse`, for better clarity. If you have
  `shock =` in existing calls to
  [`IRF()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/IRF.md),
  rename it to `impulse =`.

- [`forecast()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/forecast.md),
  [`conditional_forecast()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/conditional_forecast.md),
  and
  [`IRF()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/IRF.md)
  now transform each posterior draw before summarizing, rather than
  summarizing and then transforming. This fixes incorrect
  prediction/credible intervals and median point estimates/predictions
  when using `growth_rate_idx` argument. Mean point
  estimates/predictions are unaffected.

- Now `d_pred` is automatically created in
  [`fit()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
  to make the user-experience smoother. However, it is still left as an
  argument.

- Fixed
  [`IRF()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/IRF.md)
  plotting: specifying only one of `response`/`impulse` (leaving the
  other `NULL`) previously ignored the specified index and plotted all
  response-impulse combinations. Now the specified index is correctly
  held fixed while the other dimension varies.

- New
  [`steady_state_priors_plot()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/steady_state_priors_plot.md)
  function for visualizing steady-state priors.

- [`forecast()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/forecast.md)
  now overlays posterior steady-state estimates on forecast plots.

- Stan code has been optimized. This includes assuming that all prior
  covariance matrices are diagonal.

- Added a vignette
  [`vignette("SteadyStateBVAR-intro")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/SteadyStateBVAR-intro.md)

- Function documentations, vignettes, and the README have been revamped
  quite extensively.

## SteadyStateBVAR 0.1.1

CRAN release: 2026-07-28

- Stan code now uses the new array syntax
  ([@andrjohns](https://github.com/andrjohns),
  [\#4](https://github.com/markjwbecker/SteadyStateBVAR/issues/4)).

- Small rework of plotting in
  [`forecast()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/forecast.md)
  and corrected documentation.

## SteadyStateBVAR 0.1.0

CRAN release: 2026-07-24

- Initial CRAN release.

## SteadyStateBVAR 0.0.0.9000

- Initial development version.
- Preparing first CRAN submission.
