# 6MWD and muscle fat fraction in Duchenne muscular dystrophy (Kim 2023)

``` r

library(nlmixr2lib)
library(rxode2)
library(dplyr)
library(tidyr)
library(ggplot2)

rxode2::rxSetSeed(20230803)
set.seed(20230803)
```

Kim et al. (2023) modelled the natural history of Duchenne muscular
dystrophy (DMD) in ambulatory boys from the ImagingDMD study. Two
endpoints were modelled jointly: the 6-minute walk distance (6MWD), the
usual primary endpoint of ambulatory DMD trials, and the fat fraction of
a lower-limb muscle measured by MR spectroscopy. The paper fitted two
multivariate models, one pairing 6MWD with the vastus lateralis fat
fraction (FFVL) and one pairing it with the soleus fat fraction (FFSOL).
The package ships them as two model files that share this vignette:

``` r

mod_vl <- rxode2::rxode2(nlmixr2lib::readModelDb("Kim_2023_DMD_6MWD_FFVL"))
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_sol <- rxode2::rxode2(nlmixr2lib::readModelDb("Kim_2023_DMD_6MWD_FFSOL"))
#> ℹ parameter labels from comments will be replaced by 'label()'
```

Reference: Kim S, Willcocks RJ, Daniels MJ, Morales JF, Yoon DY,
Triplett WT, Barnard AM, Senesac CR, Lott DJ, Forbes SC, Walter GA,
Sweeney HL, Vandenborne K (2023). Multivariate modeling of magnetic
resonance biomarkers and clinical outcome measures for Duchenne muscular
dystrophy clinical trials. *CPT Pharmacometrics Syst Pharmacol*
12(10):1437-1449.
[doi:10.1002/psp4.13021](https://doi.org/10.1002/psp4.13021)

There is no drug in these models. The model time is the boy’s **age in
years**, so a simulation data set carries age in its `time` column, has
no dosing rows, and supplies four time-fixed baseline covariates:

| Column | Meaning | Units |
|----|----|----|
| `WALK_DIST_BL` | the boy’s own 6MWD at the baseline visit | m |
| `FF_VL_BL` / `FF_SOL_BL` | the boy’s own vastus lateralis / soleus fat fraction at baseline | fraction 0-1 |
| `AGE_BASE` | age at the baseline visit; below 7 years selects the “young” covariate exponents | years |
| `CONMED_STEROID` | corticosteroid use at baseline (1 = yes) | binary |

## Population

The analysis data set held 118 boys with DMD aged 5-20 years at
assessment, with 567 observations of each of 6MWD, FFVL and FFSOL (Kim
2023 Results “Data summary”). Boys were ambulatory at the derived
baseline visit, the first visit at which all three measures were
available; baseline ages were 5-13 years. Boys taking experimental
dystrophin-restoring therapies, boys likely to express dystrophin
(milder, Becker-like mutations) and boys still walking after age 15 were
excluded. By the 7-year baseline-age cut-off and steroid use, the cohort
split into 10 young and 8 old boys without steroids and 27 young and 73
old boys on steroids (Figure 2f). 6MWD values missing because a boy
could no longer walk were treated as censored (M4 method). The data came
from three US sites (University of Florida, Oregon Health & Science
University, Children’s Hospital of Philadelphia) with annual assessments
over 1-9 years.

``` r

str(mod_vl$population)
#> List of 10
#>  $ species       : chr "human"
#>  $ n_subjects    : int 118
#>  $ n_studies     : int 1
#>  $ n_observations: chr "567 observations of each of 6MWD, FFVL and FFSOL (Results 'Data summary'); missing 6MWD from inability to walk "| __truncated__
#>  $ age_range     : chr "5-20 years at assessment; baseline 5-13 years"
#>  $ sex_female_pct: num 0
#>  $ disease_state : chr "Duchenne muscular dystrophy, ambulatory at the derived baseline visit (completed the 10-m walk/run in 45 s or l"| __truncated__
#>  $ dose_range    : chr "Not applicable -- natural-history model, no drug."
#>  $ regions       : chr "United States (ImagingDMD: University of Florida, Oregon Health & Science University, Children's Hospital of Philadelphia)"
#>  $ notes         : chr "Figure 2f subgroups: no steroid + young 10, no steroid + old 8, steroid + young 27, steroid + old 73 (young = b"| __truncated__
```

## Model structure

For boy *i* at age *Age* (Kim 2023 Equations 4 and 5):

- 6MWD = S0_6MWD x (1 - DPmax_6MWD x Age^g / (DPT50^g + Age^g)), with
  DPmax_6MWD fixed to 1, so the walk distance falls from its
  early-childhood level S0_6MWD to zero along a very steep sigmoid (Hill
  coefficient about 53).
- FF = S0_FF + DPmax_FF x Age^g / (DPT50^g + Age^g): the fat fraction
  rises from S0_FF by at most DPmax_FF. DPmax_FF is logit-normal, on
  (0, 1) for FFVL and on (0, 2) for FFSOL (Supplement 2 `max=2`).

All other parameters are log-normal (Equation 2). Covariates enter as
power terms of the baseline scores centred on group medians and as
exp(beta) for steroid use (Equation 3). Every baseline-score effect has
a separate exponent and centring value for boys younger and older than 7
years at baseline; the centring values are printed only in Supplement 2.
The two endpoints are linked by correlated random effects on the ages at
half-maximal change of 6MWD and of the fat fraction and on the 6MWD
intercept. The residual error is proportional for 6MWD and Monolix
`combined1` for the fat fraction.

## Source trace

| Element | Source location | Model file |
|----|----|----|
| 6MWD sigmoid Imax, DPmax fixed 1 | Equation 4; Table 1; Supplement 2 `Imax_pop` FIXED | both |
| Fat-fraction sigmoid Emax | Equation 5; Table 1 | both |
| Logit-normal DPmax_FF on (0, 1) / (0, 2) | Equation 1; Supplement 2 `Emax` definitions | FFVL / FFSOL |
| Power and exponential covariate model | Equation 3; Table 1 betas | both |
| Young / old split at a baseline age of 7 years | Methods “Time metric, dependent variables, and covariates” | both |
| Covariate centring values (363, 367, 364 m; 0.04, 0.17, 0.06, 0.09) | Supplement 2 `[COVARIATE]` EQUATION blocks | both |
| Covariate-to-parameter mapping | Supplement 2 `[INDIVIDUAL]` DEFINITION blocks | both |
| Typical values, betas, omegas, correlations, residual SDs | Table 1, one column per model | both |
| Proportional 6MWD and combined1 fat-fraction residual error | Results “Univariate models”; Supplement 2 `[LONGITUDINAL]` | both |
| Simulation covariate profiles | Supplement 3 Table S5 | vignette |

The parameter-level trace is in the trailing comment of every
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) line of
the two model files.

## Structural verification

The typical-value solution must equal Equations 4 and 5 evaluated
directly. This compares the solver against the paper’s closed form with
the same parameters, so the bound is tight.

``` r

ages <- c(0.5, 5, 8, 10, 11, 12, 13, 15, 18)
cov_vl <- data.frame(
  WALK_DIST_BL = 330, FF_VL_BL = 0.12, AGE_BASE = 8, CONMED_STEROID = 1
)
ev_vl <- data.frame(id = 1, time = ages, cov_vl)
typ_vl <- rxode2::rxSolve(rxode2::zeroRe(mod_vl), ev_vl) |> as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalhill_walk_dist', 'etalogitemax_ff_vl', 'etalhill_ff_vl', 'etale0_ff_vl', 'etalt50max_walk_dist', 'etale0_walk_dist', 'etalt50max_ff_vl'

p <- mod_vl$theta
s0w <- exp(p[["le0_walk_dist"]]) *
  (0.12 / 0.17)^p[["e_ff_vl_bl_old_e0_walk_dist"]] *
  (330 / 367)^p[["e_walk_dist_bl_old_e0_walk_dist"]]
t50w <- exp(p[["lt50max_walk_dist"]]) *
  exp(p[["e_conmed_steroid_t50max_walk_dist"]]) *
  (0.12 / 0.17)^p[["e_ff_vl_bl_old_t50max_walk_dist"]]
gw <- exp(p[["lhill_walk_dist"]]) *
  (0.12 / 0.17)^p[["e_ff_vl_bl_old_hill_walk_dist"]]
t50f <- exp(p[["lt50max_ff_vl"]]) *
  (0.12 / 0.17)^p[["e_ff_vl_bl_old_t50max_ff_vl"]]
gf <- exp(p[["lhill_ff_vl"]])
walk_cf <- s0w * (1 - ages^gw / (t50w^gw + ages^gw))
ff_cf <- exp(p[["le0_ff_vl"]]) +
  plogis(p[["logitemax_ff_vl"]]) * ages^gf / (t50f^gf + ages^gf)

stopifnot(
  max(abs(typ_vl$walk_dist - walk_cf)) < 1e-6,
  max(abs(typ_vl$ff_vl - ff_cf)) < 1e-9
)
data.frame(age = ages, walk_dist = typ_vl$walk_dist, ff_vl = typ_vl$ff_vl)
#>    age    walk_dist      ff_vl
#> 1  0.5 3.357397e+02 0.03200000
#> 2  5.0 3.357397e+02 0.03286363
#> 3  8.0 3.357397e+02 0.06177826
#> 4 10.0 3.357389e+02 0.16959180
#> 5 11.0 3.356613e+02 0.26825968
#> 6 12.0 3.303925e+02 0.38235610
#> 7 13.0 1.866617e+02 0.48963503
#> 8 15.0 3.943674e-01 0.63568630
#> 9 18.0 5.488860e-05 0.71900716
```

For the soleus model the young-group 6MWD effect on DPmax_FF acts on the
logit of DPmax_FF / 2, and DPmax_FF stays below its bound of 2 for any
baseline:

``` r

ev_sol_young <- data.frame(
  id = 1:4, time = 20,
  WALK_DIST_BL = c(150, 300, 363, 500), FF_SOL_BL = 0.06,
  AGE_BASE = 6, CONMED_STEROID = 0
)
# zeroRe() removes every omega, so rxSolve() warns that a multi-subject
# solve has none; each row here is a typical subject by design.
typ_sol_young <- suppressWarnings(
  rxode2::rxSolve(rxode2::zeroRe(mod_sol), ev_sol_young)
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalhill_walk_dist', 'etalogitemax_ff_sol', 'etalhill_ff_sol', 'etale0_ff_sol', 'etalt50max_walk_dist', 'etale0_walk_dist', 'etalt50max_ff_sol'
p_sol <- mod_sol$theta
emax_cf <- 2 * plogis(
  p_sol[["logitemax_ff_sol"]] +
    p_sol[["e_walk_dist_bl_young_emax_ff_sol"]] * log(c(150, 300, 363, 500) / 363)
)
stopifnot(
  all(emax_cf > 0 & emax_cf < 2),
  abs(emax_cf[3] - 0.49) < 1e-9
)
data.frame(WALK_DIST_BL = ev_sol_young$WALK_DIST_BL, DPmax_FF = emax_cf)
#>   WALK_DIST_BL  DPmax_FF
#> 1          150 1.0205539
#> 2          300 0.5889115
#> 3          363 0.4900000
#> 4          500 0.3507129
```

## Published derived quantities

The paper reports that baseline steroid use delays the age at which 6MWD
has fallen by half by about 10%, and delays the corresponding age for
FFSOL by a similar amount (Results “Introducing combined covariates into
the model”).

``` r

steroid <- data.frame(
  model = c("6MWD + FFVL", "6MWD + FFSOL", "6MWD + FFSOL"),
  parameter = c("DPT50, 6MWD", "DPT50, 6MWD", "DPT50, FFSOL"),
  delay_factor = c(
    exp(mod_vl$theta[["e_conmed_steroid_t50max_walk_dist"]]),
    exp(mod_sol$theta[["e_conmed_steroid_t50max_walk_dist"]]),
    exp(mod_sol$theta[["e_conmed_steroid_t50max_ff_sol"]])
  )
)
stopifnot(all(steroid$delay_factor > 1.05 & steroid$delay_factor < 1.2))
knitr::kable(steroid, digits = 3)
```

| model        | parameter    | delay_factor |
|:-------------|:-------------|-------------:|
| 6MWD + FFVL  | DPT50, 6MWD  |        1.094 |
| 6MWD + FFSOL | DPT50, 6MWD  |        1.139 |
| 6MWD + FFSOL | DPT50, FFSOL |        1.150 |

## Simulated trajectories for the Table S5 profiles (Figure 5)

Figure 5 shows the median and the 40th-60th percentile band of 500
simulated profiles: 50 replicates of 10 boys whose covariate profiles
are listed in Supplement 3 Table S5. All ten boys are in the old
baseline-age group. The profiles are reproduced here with 20 replicates
(100 boys per steroid group), simulating between-subject variability.

``` r

s5_vl <- data.frame(
  profile = 1:10,
  WALK_DIST_BL = rep(c(295, 349), 5),
  FF_VL_BL = rep(c(0.17, 0.08), 5),
  CONMED_STEROID = rep(c(0, 1), 5)
)
s5_sol <- data.frame(
  profile = 1:10,
  WALK_DIST_BL = c(295, 349, 216, 295, 349, 216, 295, 349, 216, 295),
  FF_SOL_BL = c(0.11, 0.04, 0.17, 0.11, 0.04, 0.17, 0.11, 0.04, 0.17, 0.11),
  CONMED_STEROID = c(0, 1, 1, 0, 1, 1, 0, 1, 1, 0)
)
# Table S5 gives the baseline-age group only ('Old'); any baseline age of
# 7 years or more selects the same exponents.
s5_vl$AGE_BASE <- 9
s5_sol$AGE_BASE <- 9

n_rep <- 20
age_grid <- seq(5, 15, by = 0.25)

make_events <- function(profiles, n_rep, ages) {
  cohort <- profiles[rep(seq_len(nrow(profiles)), n_rep), ]
  cohort$id <- seq_len(nrow(cohort))
  tidyr::crossing(cohort, time = ages) |>
    dplyr::arrange(id, time) |>
    as.data.frame()
}

sim_vl <- rxode2::rxSolve(
  mod_vl, make_events(s5_vl, n_rep, age_grid),
  keep = "CONMED_STEROID"
) |> as.data.frame()
sim_sol <- rxode2::rxSolve(
  mod_sol, make_events(s5_sol, n_rep, age_grid),
  keep = "CONMED_STEROID"
) |> as.data.frame()

band <- function(sim, var) {
  sim |>
    dplyr::mutate(steroid = ifelse(CONMED_STEROID == 1, "Yes", "No")) |>
    dplyr::group_by(steroid, time) |>
    dplyr::summarise(
      lo = quantile(.data[[var]], 0.4),
      med = median(.data[[var]]),
      hi = quantile(.data[[var]], 0.6),
      .groups = "drop"
    )
}

bands <- dplyr::bind_rows(
  band(sim_vl, "walk_dist") |> dplyr::mutate(panel = "(a) 6MWD (m), FFVL model"),
  band(sim_vl, "ff_vl") |> dplyr::mutate(panel = "(b) FFVL"),
  band(sim_sol, "walk_dist") |> dplyr::mutate(panel = "(c) 6MWD (m), FFSOL model"),
  band(sim_sol, "ff_sol") |> dplyr::mutate(panel = "(d) FFSOL")
)
```

``` r

ggplot(bands, aes(time, med, colour = steroid, fill = steroid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, colour = NA) +
  geom_line() +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = c(No = "#1b2a4a", Yes = "#6ab023")) +
  scale_fill_manual(values = c(No = "#1b2a4a", Yes = "#6ab023")) +
  labs(x = "Age at assessment (years)", y = NULL, colour = "Steroid use", fill = "Steroid use") +
  theme_bw()
```

![Replicates Figure 5 of Kim 2023: median and 40th-60th percentile band
of simulated 6MWD and fat fraction for the Table S5 profiles, by steroid
use.](Kim_2023_DMD_6MWD_FF_files/figure-html/figure5-1.png)

Replicates Figure 5 of Kim 2023: median and 40th-60th percentile band of
simulated 6MWD and fat fraction for the Table S5 profiles, by steroid
use.

The medians are compared with values read off Figure 5 by the
maintainers. The early-childhood 6MWD plateau is set almost entirely by
S0_6MWD, which has little variability (omega 0.06), so it is the
sharpest test of the intercepts and the baseline-score effects; the fat
fraction at 15 years tests DPmax_FF and the age sigmoid.

``` r

fig5 <- tibble::tribble(
  ~panel, ~steroid, ~time, ~published,
  "(a) 6MWD (m), FFVL model", "No", 5, 312,
  "(a) 6MWD (m), FFVL model", "Yes", 5, 348,
  "(b) FFVL", "No", 10, 0.20,
  "(b) FFVL", "Yes", 10, 0.15,
  "(b) FFVL", "No", 15, 0.62,
  "(b) FFVL", "Yes", 15, 0.55,
  "(c) 6MWD (m), FFSOL model", "No", 5, 305,
  "(c) 6MWD (m), FFSOL model", "Yes", 5, 290,
  "(d) FFSOL", "No", 10, 0.165,
  "(d) FFSOL", "Yes", 10, 0.11,
  "(d) FFSOL", "No", 15, 0.51,
  "(d) FFSOL", "Yes", 15, 0.37
)
cmp <- fig5 |>
  dplyr::left_join(bands, by = c("panel", "steroid", "time")) |>
  dplyr::mutate(pct_diff = 100 * (med - published) / published)

stopifnot(
  # The plateau is set by the intercept parameters and baseline scores.
  all(abs(cmp$pct_diff[cmp$time == 5]) < 6),
  # Fat-fraction medians at 10 and 15 years: robust to cohort draw.
  abs(median(cmp$pct_diff)) < 10,
  quantile(abs(cmp$pct_diff), 0.9) < 20
)
cmp |>
  dplyr::select(panel, steroid, time, published, med, pct_diff) |>
  dplyr::rename(
    Panel = panel, `Steroid use` = steroid, `Age (years)` = time,
    `Figure 5 median` = published, `Simulated median` = med,
    `Difference (%)` = pct_diff
  ) |>
  knitr::kable(digits = 3)
```

| Panel | Steroid use | Age (years) | Figure 5 median | Simulated median | Difference (%) |
|:---|:---|---:|---:|---:|---:|
| \(a\) 6MWD (m), FFVL model | No | 5 | 312.000 | 308.437 | -1.142 |
| \(a\) 6MWD (m), FFVL model | Yes | 5 | 348.000 | 350.417 | 0.694 |
| \(b\) FFVL | No | 10 | 0.200 | 0.235 | 17.514 |
| \(b\) FFVL | Yes | 10 | 0.150 | 0.158 | 5.049 |
| \(b\) FFVL | No | 15 | 0.620 | 0.655 | 5.670 |
| \(b\) FFVL | Yes | 15 | 0.550 | 0.584 | 6.182 |
| \(c\) 6MWD (m), FFSOL model | No | 5 | 305.000 | 307.727 | 0.894 |
| \(c\) 6MWD (m), FFSOL model | Yes | 5 | 290.000 | 297.406 | 2.554 |
| \(d\) FFSOL | No | 10 | 0.165 | 0.166 | 0.618 |
| \(d\) FFSOL | Yes | 10 | 0.110 | 0.115 | 4.107 |
| \(d\) FFSOL | No | 15 | 0.510 | 0.458 | -10.155 |
| \(d\) FFSOL | Yes | 15 | 0.370 | 0.350 | -5.337 |

The largest difference, FFVL at 10 years without steroids, falls on the
steepest part of the curve, where a small error in reading the figure
moves the value most.

The 6MWD panels of Figure 5 level off at about 10 m rather than reaching
zero. The model itself has no floor (DPmax_6MWD is 1): values below the
walking limit were censored in the fit, and the published simulation
appears to have been plotted at that limit. The steep 6MWD decline is
therefore compared through the age at which the median falls below 150
m. In the simulation, steroid use delays that age by about 1.5 years in
both models; in Figure 5, whose curves are drawn through whole-year
ages, the delay is about one year.

``` r

loa <- bands |>
  dplyr::filter(grepl("6MWD", panel)) |>
  dplyr::group_by(panel, steroid) |>
  dplyr::summarise(age_below_150m = min(time[med < 150]), .groups = "drop")
loa_delay <- loa |>
  tidyr::pivot_wider(names_from = steroid, values_from = age_below_150m) |>
  dplyr::mutate(delay = Yes - No)
stopifnot(all(loa_delay$delay > 0.5 & loa_delay$delay < 2.5))
loa_delay |>
  dplyr::rename(
    Panel = panel, `No steroid (years)` = No, `Steroid (years)` = Yes,
    `Delay (years)` = delay
  ) |>
  knitr::kable()
```

| Panel | No steroid (years) | Steroid (years) | Delay (years) |
|:---|---:|---:|---:|
| \(a\) 6MWD (m), FFVL model | 11.50 | 13.00 | 1.5 |
| \(c\) 6MWD (m), FFSOL model | 10.75 | 12.25 | 1.5 |

## Baseline-age group

Each baseline score acts through different exponents in boys younger and
older than 7 years at baseline. For the same baseline 6MWD of 300 m, the
6MWD intercept is lower in the young group because its exponent is
larger (0.89 versus 0.63 in the FFVL model) and its centring value is
363 m rather than 367 m:

``` r

ev_age <- data.frame(
  id = 1:2, time = 5, WALK_DIST_BL = 300, FF_VL_BL = 0.17,
  AGE_BASE = c(6, 8), CONMED_STEROID = 1
)
# Two typical subjects; see the note on the zeroRe() warning above.
suppressWarnings(rxode2::rxSolve(rxode2::zeroRe(mod_vl), ev_age)) |>
  as.data.frame() |>
  dplyr::mutate(group = ifelse(AGE_BASE < 7, "young", "old")) |>
  dplyr::select(group, walk_dist, ff_vl)
#> ℹ omega/sigma items treated as zero: 'etalhill_walk_dist', 'etalogitemax_ff_vl', 'etalhill_ff_vl', 'etale0_ff_vl', 'etalt50max_walk_dist', 'etale0_walk_dist', 'etalt50max_ff_vl'
#>   group walk_dist      ff_vl
#> 1 young  298.7781 0.05956019
#> 2   old  311.7989 0.03318700
```

## PKNCA

This is a natural-history disease-progression model with no drug input,
so there is no concentration-time profile for non-compartmental
analysis. The validation above uses the closed-form solution, the
published steroid effects and the paper’s own simulation (Figure 5)
instead.

## Assumptions and deviations

- **Boundary of the baseline-age groups.** The paper splits boys at a
  baseline age of 7 years without saying which group a boy of exactly
  7.0 years joins. The models use young = `AGE_BASE < 7`. Supply a
  baseline age just below 7 to force the young group if a source defines
  the cut-off the other way.
- **Group-specific covariate columns.** Monolix carried each baseline
  score in two columns (for example `BSSMWTYOUNG` and `BSSMWTOLD`), so
  that only one group’s effect applies to any boy. The models take a
  single column per score (`WALK_DIST_BL`, `FF_VL_BL`, `FF_SOL_BL`) and
  select the group’s exponent and centring value from `AGE_BASE`, which
  gives the same individual parameters.
- **Values from Table 1; transforms from Supplement 2.** The
  `<PARAMETER>` values in the Supplement 2 Monolix code are the initial
  values of the run, not the final estimates, and are not used. The
  supplement does supply the covariate centring values, the logit-normal
  bounds of DPmax_FF and the covariate-to-parameter mapping, none of
  which are printed in the paper.
- **Logit-scale covariate on DPmax_FF (FFSOL).** Monolix applies
  covariate effects on the transformed scale, so the young-group
  baseline-6MWD effect (-1.32) is additive on logit(DPmax_FF / 2), not a
  power term on DPmax_FF.
- **Censoring.** 6MWD values missing because a boy could no longer walk
  were handled as censored (M4) in the fit. The models predict the
  uncensored walk distance, which reaches 0 m; the censoring limit is
  not part of the model.
- **Simulation variability.** Figure 5 includes between-subject
  variability, residual variability and parameter uncertainty. The
  replication here includes between-subject variability only, with 20
  rather than 50 replicates of the Table S5 profiles, and Table S5 gives
  the baseline-age group but not the age, so a baseline age of 9 years
  was used for every profile.
- **Figure 5 values** in the comparison table were read off the
  published figure by the maintainers and carry reading error of a few
  percent.
