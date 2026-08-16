
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SteadyStateBVAR

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/SteadyStateBVAR)](https://CRAN.R-project.org/package=SteadyStateBVAR)
[![CRAN
checks](https://badges.cranchecks.info/worst/SteadyStateBVAR.svg)](https://cran.r-project.org/web/checks/check_results_SteadyStateBVAR.html)
[![R-CMD-check](https://github.com/markjwbecker/SteadyStateBVAR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markjwbecker/SteadyStateBVAR/actions/workflows/R-CMD-check.yaml)
[![CRAN
downloads/day](https://cranlogs.r-pkg.org/badges/last-day/SteadyStateBVAR)](https://cranlogs.r-pkg.org/badges/last-day/SteadyStateBVAR)
[![CRAN
downloads/week](https://cranlogs.r-pkg.org/badges/last-week/SteadyStateBVAR)](https://cranlogs.r-pkg.org/badges/last-week/SteadyStateBVAR)
[![CRAN
downloads/month](https://cranlogs.r-pkg.org/badges/last-month/SteadyStateBVAR)](https://cranlogs.r-pkg.org/badges/last-month/SteadyStateBVAR)
[![CRAN downloads
total](https://cranlogs.r-pkg.org/badges/grand-total/SteadyStateBVAR)](https://cranlogs.r-pkg.org/badges/grand-total/SteadyStateBVAR)
<!-- badges: end -->

With this package, the user can estimate the steady-state BVAR model of
Villani (2009). The steady-state BVAR is simply a BVAR rewritten in
mean-adjusted form. The benefit of the mean-adjusted parametrization is
that it allows the user to specify prior beliefs about the unconditional
mean, or *steady state* of the VAR system. The model has proven very
useful for forecasting of macroeconomic variables, and is routinely used
in many central banks and other finanicial institutions (Gustafsson and
Villani, 2025).

After estimation, the user can produce forecasts (unconditional and
conditional) and impulse response functions (orthogonalized and
generalized). The goal of `SteadyStateBVAR` is to use modern Bayesian
tools (Stan) to: i) estimate the model as specified in the original
paper, and ii) extend the model in different ways. Previously,
extensions of the model seemed to be limited by what Mattias Villani had
time to derive.

See for example Clark (2011), which extends the steady-state BVAR model
to include Random Walk stochastic volatility: “*In a methodological
sense, this paper extends the estimator of Villani (2009) to include
stochastic volatility.*” Later on, he writes: “*(Special thanks are due
to Mattias Villani for providing the formulas for posterior means and
variances of* $\Pi$ *and* $\Psi$*, which generalize the
constant-variance formulas of Villani 2009.)*”

Another example is Dieppe, Legrand, and van Roye (2016): “*Villani
(2009) only provides derivation in the case of the normal-diffuse prior
distribution, so the incoming analysis will be restricted to this
case.*”

Times are different now, and with the help of Stan, we are essentially
limited only by our imagination. At the time of writing, the package
provides three versions of the steady-state BVAR model: i) the
homoscedastic model, i.e. the original model in Villani (2009), ii) the
Random Walk stochastic volatility model, i.e. the one in Clark (2011),
and iii) an AR(1) stochastic volatility model.

## Installation

You can install SteadyStateBVAR (version 0.1.1) from CRAN with:

``` r
install.packages("SteadyStateBVAR")
```

You can also install the development version of SteadyStateBVAR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("markjwbecker/SteadyStateBVAR")
```

## Vignettes

Feel free to have a look at the package vignettes by running:

``` r
vignette("Homoscedastic-steady-state-BVAR")
vignette("RW-stochastic-volatility-steady-state-BVAR")
vignette("AR1-stochastic-volatility-steady-state-BVAR")
```

## Example

To estimate the model in Section 4.1 of Villani (2009), produce
(un)conditional forecasts and perform impulse response analysis, simply
run the following code:

``` r
library(SteadyStateBVAR)
data("Villani2009")
yt <- Villani2009

#hold out last two observations to facilitate comparisons
#to Figures 1-3 in Villani (2009)
yt <- ts(yt[1:102, ], start = start(yt), frequency = frequency(yt))

bvar_obj <- bvar(data = yt)

#Use a dummy to model Sweden’s change in monetary policy in the 1990s
#(move to inflation targeting and flexible exchange rate)
bp <- which(time(yt) == 1992.75) #breakpoint
dummy_variable <- c(rep(1,bp), rep(0,nrow(yt)-bp))

bvar_obj <- setup(bvar_obj,
                  p=4,
                  deterministic = "constant_and_dummy",
                  dummy = dummy_variable)

lambda_1 <- 0.2 #overall tightness
lambda_2 <- 0.5 #cross-equation tightness
lambda_3 <- 1.0 #lag decay rate

#fol_pm = first own lag prior means
fol_pm=c(0,   #delta y_f
         0,   #pi_f
         0.9, #i_f
         0,   #delta y
         0,   #pi
         0.9, #i
         0.9  #q
         )

#95% prior probability intervals (normal distribution)
#See Table I in Villani (2009)
#These are the "steady-state priors"
#psi_1 = Psi col 1
#psi_2 = Psi col 2

theta_Psi <- 
  c(
  ppi( 2.00,  3.00,  annualized_growthrate=TRUE)$mean,   #psi_1: delta y_f
  ppi( 1.50,  2.50,  annualized_growthrate=TRUE)$mean,   #psi_1: pi_f
  ppi( 4.50,  5.50,  annualized_growthrate=FALSE)$mean,  #psi_1: i_f
  ppi( 2.00,  2.50,  annualized_growthrate=TRUE)$mean,   #psi_1: delta y
  ppi( 1.70,  2.30,  annualized_growthrate=TRUE)$mean,   #psi_1: pi
  ppi( 4.00,  4.50,  annualized_growthrate=FALSE)$mean,  #psi_1: i
  ppi( 3.85,  4.00,  annualized_growthrate=FALSE)$mean,  #psi_1: q
  ppi(-1.00,  1.00,  annualized_growthrate=TRUE)$mean,   #psi_2: delta y_f
  ppi( 1.50,  2.50,  annualized_growthrate=TRUE)$mean,   #psi_2: pi_f
  ppi( 1.50,  2.50,  annualized_growthrate=FALSE)$mean,  #psi_2: i_f
  ppi(-1.00,  1.00,  annualized_growthrate=TRUE)$mean,   #psi_2: delta y
  ppi( 4.30,  5.70,  annualized_growthrate=TRUE)$mean,   #psi_2: pi
  ppi( 3.00,  5.50,  annualized_growthrate=FALSE)$mean,  #psi_2: i
  ppi(-0.50,  0.50,  annualized_growthrate=FALSE)$mean   #psi_2: q
  )

Omega_Psi <- 
  diag(
  c(
  ppi( 2.00,  3.00,  annualized_growthrate=TRUE)$var,    #psi_1: delta y_f
  ppi( 1.50,  2.50,  annualized_growthrate=TRUE)$var,    #psi_1: pi_f
  ppi( 4.50,  5.50,  annualized_growthrate=FALSE)$var,   #psi_1: i_f
  ppi( 2.00,  2.50,  annualized_growthrate=TRUE)$var,    #psi_1: delta y
  ppi( 1.70,  2.30,  annualized_growthrate=TRUE)$var,    #psi_1: pi
  ppi( 4.00,  4.50,  annualized_growthrate=FALSE)$var,   #psi_1: i
  ppi( 3.85,  4.00,  annualized_growthrate=FALSE)$var,   #psi_1: q
  ppi(-1.00,  1.00,  annualized_growthrate=TRUE)$var,    #psi_2: delta y_f
  ppi( 1.50,  2.50,  annualized_growthrate=TRUE)$var,    #psi_2: pi_f
  ppi( 1.50,  2.50,  annualized_growthrate=FALSE)$var,   #psi_2: i_f
  ppi(-1.00,  1.00,  annualized_growthrate=TRUE)$var,    #psi_2: delta y
  ppi( 4.30,  5.70,  annualized_growthrate=TRUE)$var,    #psi_2: pi
  ppi( 3.00,  5.50,  annualized_growthrate=FALSE)$var,   #psi_2: i
  ppi(-0.50,  0.50,  annualized_growthrate=FALSE)$var    #psi_2: q
  )
  )

bvar_obj <- priors(bvar_obj,
                   lambda_1,
                   lambda_2,
                   lambda_3,
                   fol_pm,
                   theta_Psi,
                   Omega_Psi,
                   Jeffreys=TRUE) #FALSE for uninformative inverse-Wishart

p <- bvar_obj$setup$p
k <- bvar_obj$setup$k
kf <- 3 #first three variables in yt are foreign

restriction_matrix <- matrix(1, k*p, k)

for(i in 1:p){
  rows <- ((i-1)*k + kf + 1) : (i*k)
  cols <- 1:kf
  restriction_matrix[rows, cols] <- 0
}
print(restriction_matrix)
#block exogeneity for foreign variables
bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)

H <- 12 #forecast horizon
(d_pred <- cbind(rep(1, 12), 0)) #future d_t values


#fit the model
bvar_obj <- fit(bvar_obj,
                H = H,
                d_pred = d_pred,
                iter = 10000,
                warmup = 2500,
                chains = 2,
                cores = 2)

#posterior summaries
summary(bvar_obj , stat = "mean")

#you can look at the stanfit object directly
stan_fit <- bvar_obj$fit$stan
print(stan_fit)

#unconditional forecasts
#see last forecasts in Figures 1-3 in Villani (2009)
fcst <- forecast(bvar_obj,
                 pi = 0.68, #pi = prediction interval
                 fcst_type = "mean",
                 growth_rate_idx = c(4,5), #convert QoQ forecasts to YoY
                 plot_idx = c(4,5,6))

#conditional forecasts
#Toy scenario: inflation gets really high
#What will happen to domestic interest rate?
conditions <- data.frame(
              var        = rep(5,12),
              horizon    = rep(1:12),
              value      = c(1.0,1.5,2.0,1.8,
                             1.5,1.2,1.0,1.0,
                             rep(0.5,4)) #QoQ scale for inflation here
              )
              
cond_fcst <- conditional_forecast(bvar_obj,
                    conditions,
                    pi=0.68,
                    fcst_type = "mean",
                    plot_idx = c(5,6),
                    growth_rate_idx = c(5)) #convert QoQ forecasts to YoY

#impulse response analysis
irf <- IRF(bvar_obj,
           H=20,
           response=5,#inflation
           shock=6, #interest rate
           type="median",
           method="OIRF",
           ci=0.68,
           growth_rate_idx=5) #YoY inflation instead of QoQ
```

## References

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327–341.

Dieppe, A., Legrand, R., and van Roye, B. (2016). The BEAR toolbox.
Working Paper Series, No. 1934. European Central Bank.

Gustafsson, O., and Villani, M. (2025). Variational inference for
steady-state BVARs. arXiv preprint arXiv:2506.09271.

Karlsson, S. (2013). Forecasting with Bayesian vector autoregression.
In: Elliott, G. and Timmermann, A. (eds), *Handbook of Economic
Forecasting*. Elsevier B.V., Vol. 2, Part B, pp. 791–897.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
