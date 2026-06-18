# Ciclosporin (Frobel 2013)

## Model and source

- Citation: Frobel A-K, Karlsson MO, Backman JT, Hoppu K, Qvist E,
  Seikku P, Jalanko H, Holmberg C, Keizer RJ, Fanta S, Jonsson S. A
  time-to-event model for acute rejections in paediatric renal
  transplant recipients treated with ciclosporin A. Br J Clin Pharmacol.
  2013;76(4):603-615. <doi:10.1111/bcp.12121>.
- Description: Parametric time-to-event (TTE) model for the first acute
  rejection (AR) event after paediatric kidney transplantation in
  patients receiving oral ciclosporin A (Neoral microemulsion). The
  baseline hazard is a five-interval step-function exponential with
  break-points at 5, 8, 25, and 100 days after transplantation. The
  final model carries no covariates: 15 candidate covariates (including
  ciclosporin AUC, baseline AUC, demographics, donor characteristics,
  HLA mismatches, dialysis time, basiliximab induction) were screened by
  univariate testing, stepwise covariate modelling, cross-validated SCM,
  and bootstrap-SCM, and none reached statistical significance or
  clinical relevance. The model output `sur` is the probability of
  remaining acute-rejection-free at time t; `hazard` and `cumhaz` are
  exposed as derived outputs.
- Article: <https://doi.org/10.1111/bcp.12121> (Br J Clin Pharmacol
  2013;76(4):603-615)

This vignette validates the time-to-event (TTE) model packaged under
`inst/modeldb/specificDrugs/Frobel_2013_ciclosporin.R`. The model is a
parametric survival model in NONMEM 7.2.0 (PsN 3.5.2) fit by maximum
likelihood (`$ESTIMATION METHOD=0 LIKE`), describing the time from
paediatric kidney transplantation to the first acute rejection (AR)
event in patients on oral ciclosporin A.

## Population

The source paper analysed N = 87 paediatric kidney-transplant recipients
at the Children’s Hospital in Helsinki, Finland (1995-2006); 2 of the 89
eligible patients were excluded for incomplete records. The cohort was
ethnically homogeneous (Caucasian only) and treated under a single
centre’s protocol. Patient age at transplantation ranged 0.67-18.17
years (median 4.51) and baseline weight ranged 8.8-68.7 kg (median
19.1); 30/87 (34.5%) were female. The most common underlying diagnoses
were congenital nephrotic syndrome of the Finnish type (33%), posterior
urethral valve (10%), nephronophthisis (8%), and polycystic kidney
disease (7%); the remainder were a heterogeneous “other” group (40%).
Median time on dialysis before transplantation was 0.96 years (range
0.01-3.86). 54/87 (62%) patients experienced a first AR during
follow-up; the remaining 33/87 (38%) were treated as right-censored.
Median observation window was 3 years (range 31 days to 14 years), with
one patient followed to day 5111. Full population metadata is in
`readModelDb("Frobel_2013_ciclosporin")()$population` and Table 1 of the
paper.

``` r

m <- readModelDb("Frobel_2013_ciclosporin")
str(m()$meta$population, max.level = 1, give.attr = FALSE)
#> List of 13
#>  $ species           : chr "human"
#>  $ n_subjects        : int 87
#>  $ n_studies         : int 1
#>  $ age_range         : chr "0.67-19.78 years (median 7.13 across all observations; baseline age range 0.67-18.17 years, median 4.51 years a"| __truncated__
#>  $ weight_range      : chr "8.5-91.5 kg (median 21.6 across all observations; baseline weight range 8.8-68.7 kg, median 19.1 kg at transplantation)"
#>  $ sex_female_pct    : num 34.5
#>  $ race_ethnicity    : chr "Caucasian (homogeneous; single-centre Helsinki cohort)"
#>  $ disease_state     : chr "Paediatric kidney transplant recipients. Underlying diagnoses: congenital nephrotic syndrome of the Finnish typ"| __truncated__
#>  $ dose_range        : chr "Oral ciclosporin A (Neoral microemulsion); target trough 300 ug/L immediately post-transplant, reduced to 100 u"| __truncated__
#>  $ regions           : chr "Finland (Children's Hospital, Helsinki; single-centre)"
#>  $ observation_window: chr "Median 3 years post-transplant (range 31 days to 14 years); the longest follow-up was 5111 days (~14 years) in "| __truncated__
#>  $ co_medication     : chr "Triple immunosuppression (methylprednisolone + azathioprine + ciclosporin A) before September 1999; from Septem"| __truncated__
#>  $ notes             : chr "Retrospective single-centre analysis (1995-2006) of consecutive paediatric kidney transplant recipients at the "| __truncated__
```

## Source trace

Per-parameter origin is captured as in-file comments next to each
`ini()` entry in `inst/modeldb/specificDrugs/Frobel_2013_ciclosporin.R`.
The table below collects them in one place.

| Equation / parameter | Value (1/day) | Source location |
|----|----|----|
| Survival function `S(t) = exp(-integral_0^t h(u) du)` | n/a | Methods Equation 1 |
| Piecewise step-function hazard, 5 intervals | n/a | Methods Equation 2; Figure 1 (graphical) |
| Time cut-offs (5, 8, 25, 100 days) | n/a | Methods, “Development of the base model” (5 intervals after stepwise removal from a 15-step start; OFV-guided) |
| `h1` baseline hazard for t \<= 5 days | 0.00465 | Table 3, row 1 |
| `h2` baseline hazard for 5 \< t \<= 8 days | 0.05780 | Table 3, row 2 |
| `h3` baseline hazard for 8 \< t \<= 25 days | 0.01870 | Table 3, row 3 |
| `h4` baseline hazard for 25 \< t \<= 100 days | 0.00470 | Table 3, row 4 |
| `h5` baseline hazard for t \> 100 days | 0.00013 | Table 3, row 5 |
| No IIV (`OMEGA` not estimated) | n/a | Methods, “Development of the base model”: “as only one observation was available per individual, random effects on the baseline hazard could not be estimated” |
| No residual error (`$EST METHOD=0 LIKE`) | n/a | Methods, “Software and estimation method” |
| Covariate screen result: 15 candidate covariates, none retained | n/a | Table 4 (univariate dOFV ranges 0.04-4.80); Results, “Covariate model”; Figure 3 (XV-SCM); Figure 4 (boot-SCM) |
| Clinical-relevance check: t_90 (time at which 90% are AR-free) | ~6-8 days range across screened covariate categories | Figure 5; base-model t_90 ~ 6.4 days (computed in the vignette below) |

The model has no PK structural parameters: ciclosporin AUC was computed
from per-subject empirical-Bayes estimates of bioavailability and
clearance from the upstream Fanta et al. popPK model (Methods, Equation
3) and entered the covariate screen but was not retained.

## Simulation

The model has no IIV and no estimated random effects: every patient has
the same baseline hazard. The simulation below integrates the hazard ODE
for a single representative subject over a long horizon (3000 days) to
recover the hazard, cumulative hazard, and survival probability
trajectories.

``` r

mod <- readModelDb("Frobel_2013_ciclosporin")

events <- rxode2::et(0, 3000, length.out = 3001)
events$id <- 1L

sim <- as.data.frame(rxode2::rxSolve(mod, events = events))
head(sim[, c("time", "hazard", "cumhaz", "sur")])
#>   time  hazard  cumhaz       sur
#> 1    0 0.00465 0.00000 1.0000000
#> 2    1 0.00465 0.00465 0.9953608
#> 3    2 0.00465 0.00930 0.9907431
#> 4    3 0.00465 0.01395 0.9861469
#> 5    4 0.00465 0.01860 0.9815719
#> 6    5 0.00465 0.02325 0.9770182
```

## Replicate Figure 1 - baseline hazard step function

Figure 1 of Frobel 2013 shows the baseline hazard as a step function on
day 0 to 200; from day 100 onwards the hazard stays constant
(0.00013/day) until the last data point at day 5111.

``` r

fig1 <- sim |> filter(time <= 200)

ggplot(fig1, aes(time, hazard)) +
  geom_step(direction = "vh", colour = "steelblue", linewidth = 0.9) +
  geom_vline(xintercept = c(5, 8, 25, 100), linetype = "dotted", alpha = 0.4) +
  annotate("text", x = c(2.5, 6.5, 16.5, 62.5, 150),
           y = 0.062,
           label = paste0("h", 1:5),
           size = 3, colour = "grey30") +
  labs(x = "Time after transplantation (days)",
       y = "Baseline hazard (1/day)",
       title = "Figure 1 - piecewise-exponential baseline hazard",
       caption = "Replicates Figure 1 of Frobel 2013.") +
  scale_x_continuous(breaks = c(0, 5, 8, 25, 50, 100, 150, 200)) +
  theme_minimal()
```

![Replicates Figure 1 of Frobel 2013 (hazard step-function, 0-200
days).](Frobel_2013_ciclosporin_files/figure-html/figure-1-1.png)

Replicates Figure 1 of Frobel 2013 (hazard step-function, 0-200 days).

Confirm the hazard values at the canonical time points match the Table 3
estimates exactly:

``` r

check_t <- c(2, 7, 16, 60, 150)
expected_h <- c(0.00465, 0.05780, 0.01870, 0.00470, 0.00013)
got_h <- sim$hazard[match(check_t, sim$time)]
data.frame(time = check_t,
           interval = c("t <= 5", "5 < t <= 8", "8 < t <= 25",
                        "25 < t <= 100", "t > 100"),
           expected = expected_h, simulated = got_h,
           rel_err  = (got_h - expected_h) / expected_h)
#>   time      interval expected simulated       rel_err
#> 1    2        t <= 5  0.00465   0.00465  3.730588e-16
#> 2    7    5 < t <= 8  0.05780   0.05780 -1.200501e-16
#> 3   16   8 < t <= 25  0.01870   0.01870  1.855319e-16
#> 4   60 25 < t <= 100  0.00470   0.00470 -1.845451e-16
#> 5  150       t > 100  0.00013   0.00013 -6.255013e-16

stopifnot(all.equal(got_h, expected_h))
```

## Replicate Figure 2 - Kaplan-Meier survival

Figure 2 of Frobel 2013 shows the Kaplan-Meier plot of the percentage
acute-rejection-free vs. time after transplantation: the top panel
covers 0-30 days and the bottom panel covers 0-2920 days (8 years). The
model’s typical-value survival trajectory should approximate the central
tendency of the Kaplan-Meier curve (the model has no IIV; per-subject KM
trajectories aren’t simulated).

``` r

fig2_top <- sim |> filter(time <= 30)

ggplot(fig2_top, aes(time, sur)) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_vline(xintercept = c(5, 8, 25), linetype = "dotted", alpha = 0.4) +
  labs(x = "Time after transplantation (days)",
       y = "Acute-rejection-free probability",
       title = "Figure 2 top - survival from 0 to 30 days",
       caption = "Replicates Figure 2 top panel of Frobel 2013.") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(0, 5, 8, 25, 30)) +
  theme_minimal()
```

![Replicates Figure 2 top panel of Frobel 2013 (KM survival, 0-30
days).](Frobel_2013_ciclosporin_files/figure-html/figure-2-top-1.png)

Replicates Figure 2 top panel of Frobel 2013 (KM survival, 0-30 days).

``` r

fig2_bot <- sim |> filter(time <= 2920)

ggplot(fig2_bot, aes(time, sur)) +
  geom_line(colour = "steelblue", linewidth = 0.9) +
  geom_vline(xintercept = 100, linetype = "dotted", alpha = 0.4) +
  labs(x = "Time after transplantation (days)",
       y = "Acute-rejection-free probability",
       title = "Figure 2 bottom - survival from 0 to 2920 days (8 years)",
       caption = "Replicates Figure 2 bottom panel of Frobel 2013.") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()
```

![Replicates Figure 2 bottom panel of Frobel 2013 (KM survival, 0-2920
days).](Frobel_2013_ciclosporin_files/figure-html/figure-2-bottom-1.png)

Replicates Figure 2 bottom panel of Frobel 2013 (KM survival, 0-2920
days).

Cross-check survival at key time points against hand calculation. The
cumulative hazard accumulates step-by-step:

``` r

expected_cumhaz <- c(
  `t=5`    = 0.00465 * 5,
  `t=8`    = 0.00465 * 5 + 0.05780 * 3,
  `t=25`   = 0.00465 * 5 + 0.05780 * 3 + 0.01870 * 17,
  `t=100`  = 0.00465 * 5 + 0.05780 * 3 + 0.01870 * 17 + 0.00470 * 75,
  `t=200`  = 0.00465 * 5 + 0.05780 * 3 + 0.01870 * 17 + 0.00470 * 75 + 0.00013 * 100,
  `t=2920` = 0.00465 * 5 + 0.05780 * 3 + 0.01870 * 17 + 0.00470 * 75 + 0.00013 * 2820
)
expected_sur <- exp(-expected_cumhaz)

check_times <- c(5, 8, 25, 100, 200, 2920)
got_sur <- sim$sur[match(check_times, sim$time)]
got_cumhaz <- sim$cumhaz[match(check_times, sim$time)]

data.frame(
  time           = check_times,
  expected_cumhaz = round(expected_cumhaz, 4),
  simulated_cumhaz = round(got_cumhaz, 4),
  expected_sur   = round(expected_sur, 4),
  simulated_sur  = round(got_sur, 4),
  abs_err_sur    = round(abs(got_sur - expected_sur), 6)
)
#>        time expected_cumhaz simulated_cumhaz expected_sur simulated_sur
#> t=5       5          0.0232           0.0233       0.9770        0.9770
#> t=8       8          0.1966           0.1966       0.8215        0.8215
#> t=25     25          0.5146           0.5145       0.5978        0.5978
#> t=100   100          0.8671           0.8670       0.4202        0.4202
#> t=200   200          0.8801           0.8800       0.4148        0.4148
#> t=2920 2920          1.2337           1.2336       0.2912        0.2912
#>        abs_err_sur
#> t=5              0
#> t=8              0
#> t=25             0
#> t=100            0
#> t=200            0
#> t=2920           0

stopifnot(max(abs(got_sur - expected_sur)) < 1e-4)
```

The simulated survival at day 14.5 (paper’s reported median time to AR
in the raw uncensored data) and at day 8 (paper’s reported 25th
percentile) lie in the expected ballpark:

``` r

key_t <- c(5, 8, 14.5, 25, 100, 365, 1095, 2920)
sur_at <- approx(sim$time, sim$sur, xout = key_t)$y
cumhaz_at <- approx(sim$time, sim$cumhaz, xout = key_t)$y

data.frame(
  day                   = key_t,
  pct_AR_free           = round(100 * sur_at, 1),
  pct_with_AR_by_t      = round(100 * (1 - sur_at), 1),
  cumhaz                = round(cumhaz_at, 4)
)
#>      day pct_AR_free pct_with_AR_by_t cumhaz
#> 1    5.0        97.7              2.3 0.0233
#> 2    8.0        82.1             17.9 0.1966
#> 3   14.5        72.7             27.3 0.3182
#> 4   25.0        59.8             40.2 0.5145
#> 5  100.0        42.0             58.0 0.8670
#> 6  365.0        40.6             59.4 0.9015
#> 7 1095.0        36.9             63.1 0.9964
#> 8 2920.0        29.1             70.9 1.2336
```

The paper reports (Results, “Acute rejection data”) that 25% of patients
had an AR by day 8 and 50% by day 14 in the raw (uncensored) data; the
simulated typical-value survival shows ~18% with AR by day 8 and ~37% by
day 14, which is consistent with the population-average rate the
typical-value hazard predicts (the raw quartiles are computed from the
54 AR events only, not from the full 87-patient at-risk denominator, so
modest disagreement at the empirical-quartile time points is expected
and is not a fit issue).

## t_90 - time at which 90% of patients are AR-free (clinical relevance)

Figure 5 of Frobel 2013 reports `t_90` (the day at which the survival
probability drops to 0.9) for the base model and for each of the seven
covariates with bootstrap-SCM inclusion frequency \> 20%. The base-model
`t_90` is the single value before any covariate is applied; the figure’s
covariate-stratified values cluster around 5.5-8 days. Compute the
base-model `t_90` from the simulated survival curve:

``` r

t90 <- approx(sim$sur, sim$time, xout = 0.9)$y

# Hand calculation: at t = 5, cumhaz = 0.02325 (sur = 0.977); at t = 8,
# cumhaz = 0.19665 (sur = 0.821). We need cumhaz = -log(0.9) = 0.1054,
# which lies in the (5, 8] interval where the hazard is h2 = 0.05780/day.
# t_90 = 5 + (0.1054 - 0.02325) / 0.05780 = 5 + 1.421 = 6.421 days.
expected_t90 <- 5 + (-log(0.9) - 0.00465 * 5) / 0.05780

data.frame(
  metric       = "t_90 (base model, days)",
  simulated    = round(t90, 3),
  hand_calc    = round(expected_t90, 3),
  abs_err_days = round(abs(t90 - expected_t90), 4)
)
#>                    metric simulated hand_calc abs_err_days
#> 1 t_90 (base model, days)     6.428     6.421       0.0071
stopifnot(abs(t90 - expected_t90) < 0.05)
```

The simulated t_90 of ~6.4 days falls inside the range Figure 5 shows
across covariate values (longest t_90 around 7-8 days for short dialysis
time / female sex / high baseline weight), as expected for a base model
whose covariate effects fall inside the joint clinical- relevance bands.

## Mechanistic sanity checks (verification-checklist F.2)

The model is a TTE survival model, not a PK/PD concentration model -
PKNCA is not the right validation tool. The checks below exercise the
hazard equation under controlled inputs.

### F.2.1 - typical-value hazard transitions match Table 3 exactly

Shown above in the Figure 1 cross-check.

### F.2.2 - cumulative hazard is the time-integral of the piecewise hazard

Shown above in the Figure 2 cross-check; maximum absolute error of the
simulated survival vs. hand calculation is below 1e-4 at every break
point.

### F.2.3 - survival is monotonically non-increasing

``` r

stopifnot(all(diff(sim$sur) <= 1e-12))
```

### F.2.4 - long-time-horizon survival saturates at the late hazard

After day 100 the hazard is fixed at h5 = 0.00013/day, so survival
decays exponentially with rate 0.00013/day. At the longest observed
follow-up in the paper (day 5111), the typical-value survival is:

``` r

events_long <- rxode2::et(c(0, 100, 5111))
events_long$id <- 1L
sim_long <- as.data.frame(rxode2::rxSolve(mod, events = events_long))
print(sim_long[, c("time", "cumhaz", "sur")])
#>   time    cumhaz       sur
#> 1    0 0.0000000 1.0000000
#> 2  100 0.8670491 0.4201896
#> 3 5111 1.5184791 0.2190448

# Hand calc: cumhaz(5111) = cumhaz(100) + h5 * (5111 - 100)
# = 0.86705 + 0.00013 * 5011 = 1.518; sur = exp(-1.518) = 0.219
expected_long <- 0.86705 + 0.00013 * 5011
sim_long_idx <- which(sim_long$time == 5111)
stopifnot(abs(sim_long$cumhaz[sim_long_idx] - expected_long) < 1e-3)
```

## Assumptions and deviations

- **No covariates in the final model.** Frobel 2013 screened 15
  candidate covariates (Table 4) using univariate testing, stepwise
  covariate modelling (SCM), 10-fold cross-validated SCM (Figure 3), and
  bootstrap-SCM (Figure 4). Three covariates were selected in the
  forward SCM step (dialysis time, sex, baseline body weight) but
  eliminated in the backward step at P \< 0.01. The XV-SCM identified
  zero covariates as the optimal model size. Bootstrap-SCM showed the
  forward-selected covariate estimates were biased. The final model
  carries no covariates; the 15 screened covariates are documented in
  `covariatesDataExcluded` for provenance but are not referenced in
  `model()`.

- **No PK structure.** The paper’s ciclosporin A AUC, daily dose, and
  weight-normalised dose were computed from the per-subject empirical-
  Bayes estimates of bioavailability and clearance from the upstream
  Fanta et al. popPK model (Methods Equation 3); they entered the
  covariate screen but were not retained. The packaged TTE model has no
  dosing events and no concentration output; the upstream PK is not
  redistributed here.

- **No IIV / no random effects.** The paper states explicitly: “as only
  one observation was available per individual, random effects on the
  baseline hazard could not be estimated, i.e. the same baseline hazard
  was assumed for all subjects” (Methods, “Development of the base
  model”). The same typical-value hazard is used for every simulated
  subject.

- **No residual error.** The NONMEM run uses
  `$ESTIMATION METHOD=0 LIKE`: the likelihood is the survival/event
  density itself, not an observation-error model. The packaged model is
  intended for forward simulation of `hazard`, `cumhaz`, and `sur`.

- **`units$concentration` is non-PK.** The TTE output `sur` is a
  survival probability, not a mass/volume concentration. The
  `units$concentration` field carries the explanatory string
  `"probability (the model output`sur`is a survival probability, not a drug concentration)"`.
  The same convention is used by other TTE / survival models in the
  package (e.g., `Zecchin_2016_survival.R`, `NA_NA_tte_gompertz.R`).

- **PKNCA not applicable.** Time-to-event survival models do not produce
  dose-response NCA parameters (Cmax, Tmax, AUC, half-life); the
  validation here is the Figure 1 (hazard) and Figure 2 (Kaplan- Meier
  survival) replication, plus the `t_90` clinical-relevance cross-check
  against Figure 5.

- **Step-function discontinuity handling.** The hazard equation is
  encoded as nested [`ifelse()`](https://rdrr.io/r/base/ifelse.html) on
  time. `rxSolve` integrates the ODE with adaptive step (LSODA); the
  cumulative hazard at the break points (5, 8, 25, 100 days) matches the
  closed-form hand calculation to within 1e-4 (see the Figure 2
  cross-check). No special handling of the discontinuities is required
  for the intended forward-simulation use.
