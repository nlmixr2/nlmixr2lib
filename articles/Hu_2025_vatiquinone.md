# Vatiquinone (Hu 2025)

## Model and source

- Citation: Hu Y, Gao L, Lee L, Cherry JJ, Kong R. Characterizing
  Population Pharmacokinetics of Vatiquinone in Healthy Volunteers and
  Patients with Friedreich’s Ataxia. Pharmaceuticals. 2025;18(9):1339.
  <doi:10.3390/ph18091339>
- Description: Two-compartment population pharmacokinetic model for
  vatiquinone (PTC743), a first-in-class 15-lipoxygenase inhibitor
  developed for Friedreich’s ataxia and other mitochondrial diseases,
  with parallel zero-order and first-order oral absorption and linear
  elimination (Hu 2025; 343 participants and 4,608 quantifiable plasma
  samples pooled from eight phase I/II/III studies in adult healthy
  volunteers and adult and pediatric patients with Friedreich’s ataxia
  or other mitochondrial diseases). 74.4% of the absorbed dose enters
  through a first-order arm delayed by a 2.79 h lag time; the remaining
  25.6% enters the central compartment as a 6.03 h zero-order input, a
  dual pathway the authors attribute to the compound’s extreme
  lipophilicity (cLogP 7.8) and partial lymphatic uptake via
  chylomicrons. Vatiquinone exposure is dominated by prandial state:
  relative to the reference medium-fat meal, a liquid PediaSure
  supplement gives 6.9% and the fasted state 3.6% of the reference
  exposure, i.e. a medium-fat meal raises exposure roughly 14-fold and
  28-fold respectively. Strong CYP3A4 modulation moves apparent
  clearance in both directions (itraconazole to 23.5%, rifampicin to
  202% of the monotherapy value). Patients with Friedreich’s ataxia
  carry a 50.1% lower relative bioavailability and a 40.6% lower
  apparent clearance, which combine to a net 19% lower steady-state AUC.
  Apparent clearance also scales with body weight (power exponent 0.915,
  reference 65 kg) and inversely with body mass index (exponent -0.975,
  reference 21.6 kg/m^2), and the central volume scales linearly with
  body weight.
- Article: <https://doi.org/10.3390/ph18091339>
- Supplement (Tables S1-S7, Figures S1-S5): downloadable from the
  article page at <https://www.mdpi.com/article/10.3390/ph18091339/s1>.

Vatiquinone (PTC743) is a 15-lipoxygenase inhibitor developed for
Friedreich’s ataxia and other mitochondrial diseases. What makes its
population PK unusual is the size of the prandial effect: the compound
has a cLogP of 7.8, is insoluble in water, and is thought to reach the
systemic circulation partly through chylomicron-mediated lymphatic
transport, so exposure moves roughly 28-fold between the fasted state
and a medium-fat meal. That single covariate dominates the model.

## Population

343 participants contributed 4608 quantifiable plasma samples across 8
phase I/II/III studies (Hu 2025 Table 1). The cohort is 116 (33.8%)
adult healthy volunteers, 173 (50.4%) patients with Friedreich’s ataxia
(adult and pediatric), and 54 (15.7%) pediatric patients with other
mitochondrial diseases (epilepsy indication).

Age ranged over 1-67 years (median 20); 144/343 (42.0%) pediatric (\<18
years), including 29 participants under 7 years of age, body weight over
6.30-119 kg (median 58.6) and body mass index over 11.9-41.4 kg/m^2
(median 21.6); 51.9% were female. Race and ethnicity: White 279 (81.3%),
Black or African American 35 (10.2%), Asian 12 (3.5%), not reported 17
(5.0%); Hispanic or Latino ethnicity 76 (22.2%). Dosing: 120-1,400 mg
(median 400 mg), given as single doses or three times daily. Capsule 284
(82.8%) and oral solution 59 (17.2%). Pediatric studies dosed 15 mg/kg
below 13 kg body weight and 200 mg at or above 13 kg.

The pooled design has one feature that matters for every prandial
conclusion below. Of the eight studies, only EPI743-12-001 varied the
meal: it was an 18-subject three-way crossover of fasted, liquid
PediaSure and medium-fat meal at 300 mg. Every other study dosed
exclusively with a medium-fat meal (Supplementary Table S1: 331 of 343
participants medium-fat, 6 fasted, 6 liquid). The two prandial covariate
coefficients are therefore identified by that one small crossover, even
though they are the largest effects in the model.

## Source trace

Every `ini()` value and every non-obvious `model()` equation, with the
place in Hu 2025 it came from.

| Quantity | Value | Source |
|:---|:---|:---|
| FK0 (first-order fraction) | 0.744 | Table 2 theta 1 |
| Ka | 0.200 1/h | Table 2 theta 2 |
| TK0 (zero-order duration D2) | 6.034 h | Table 2 theta 3 |
| TLAG1 (first-order lag) | 2.787 h | Table 2 theta 4 |
| V/F | 180.748 L | Table 2 theta 5 |
| CL/F | 162.721 L/h | Table 2 theta 6 |
| V2/F | 4852.69 L | Table 2 theta 7 |
| Q/F | 67.896 L/h | Table 2 theta 8 |
| Residual (log-scale SD) | 1.062 | Table 2 row 9 ‘Additive residual’ |
| Itraconazole on CL/F | -1.446 | Table 2 theta 10 |
| Rifampicin on CL/F | 0.704 | Table 2 theta 11 |
| Liquid PediaSure on FK0 | -2.671 | Table 2 theta 12 |
| Fasted on FK0 | -3.324 | Table 2 theta 13 |
| Body weight on CL/F | 0.915 | Table 2 theta 14 |
| Friedreich ataxia on CL/F | -0.406 | Table 2 theta 15 |
| Friedreich ataxia on FK0 | -0.501 | Table 2 theta 16 |
| BMI on CL/F | -0.975 | Table 2 theta 17 |
| Body weight on V/F | 1 (fixed) | Section 2.3 V equation (printed literal; no Table 2 row) |
| IIV on CL/F | 0.191 (variance) | Table 2 IIV block; printed %CV 45.880 = sqrt(exp(0.191) - 1) |
| Reference body weight | 65 kg | Section 2.3 text |
| Reference BMI | 21.6 kg/m^2 | Section 2.3 text (population median; Supplementary Table S2 confirms) |
| Covariate equation forms | see below | Section 2.3 displayed equations |
| Concentration unit | ng/mL | Figure 2 VPC axis; Supplementary Figures S1 and S4 axes |
| Two-compartment, dual absorption | structure | Figure 5 schematic; Section 4.6 |

Source trace for Hu_2025_vatiquinone. {.table}

The four covariate equations are printed in Hu 2025 Section 2.3 and are
reproduced in `model()` exactly as given:

    Ka_i   = TVKa * exp(eta_Ka)
    FK0_i  = TVFK0 * e^(LQD_i * t12) * e^(FST_i * t13) * (1 + FA_i * t16)
    CL_i   = TVCL * (1 + FA_i * t15) * e^(ITR_i * t10) * e^(RFM_i * t11)
                    * (BWT_i / 65)^t14 * (BMI_i / 21.6)^t17 * exp(eta_CL)
    V_i    = TVV * (BWT_i / 65)^1 * exp(eta_V)

Note the deliberate mixture of functional forms: the two disease effects
are linear deviations `(1 + FA * theta)` while the prandial and
comedication effects are log-linear `exp(x * theta)`. That is what lets
the Discussion quote theta 15 and theta 16 straight back as “a 50%
reduction in relative bioavailability (FK0) and a 40% decrease in
clearance”.

## Model structure and the dosing idiom

The absorption model is two parallel arms feeding one central
compartment (Hu 2025 Figure 5): a first-order arm carrying 74.4% of the
dose, delayed by a 2.79 h lag, and a zero-order arm carrying the
remaining 25.6% straight into central over 6.03 h.

Simulating this requires **two dose records per administration**: a
bolus to `depot` and a modelled-duration record (`rate = -2`) of the
same amount to `central`. The model’s `f()` statements split the dose
between them, so nothing is double counted.

``` r

# Reference covariate condition: healthy volunteer, 65 kg, BMI 21.6, medium-fat
# meal, no comedication. This is the reference for every ratio in Hu 2025 Table 3.
ref_cov <- list(
  WT = 65, BMI = 21.6, DIS_FRDA = 0, FED = 1, FED_LIQUIDSUPP = 0,
  CONMED_ITRACONAZOLE = 0, CONMED_RIFAMPICIN = 0
)

# Build an event table for one subject: paired depot / central dose records
# plus observation rows on the `central` ODE state.
vatiq_events <- function(id, dose, obs_times, ii = 0, addl = 0, covariates = ref_cov) {
  dosing <- data.frame(
    id = id, time = 0, amt = dose, evid = 1L,
    cmt = c("depot", "central"), rate = c(0, -2), ii = ii, addl = addl
  )
  obs <- data.frame(
    id = id, time = obs_times, amt = NA_real_, evid = 0L,
    cmt = "central", rate = 0, ii = 0, addl = 0
  )
  cbind(rbind(dosing, obs), as.data.frame(covariates))
}
```

## Structural check: steady-state mass balance

For a linear model at true steady state the AUC over one day of
three-times-daily dosing is exactly `frel * 3 * Dose / CL`, whatever the
absorption model does. Checking the solved AUC against that closed form
gates the whole disposition block, the dose-splitting idiom and the
ng/mL unit conversion at once. Because both sides use the same parameter
values, a tight bound is the right gate here.

The terminal half-life of this model is about 70 h at the reference
condition, so “steady state” needs a long run: 300 three-times-daily
doses (100 days).

``` r

mod_typ <- rxode2::zeroRe(readModelDb("Hu_2025_vatiquinone"))
n_dose_ss <- 300
tau <- 8

ev_ss <- vatiq_events(
  id = 1, dose = 400,
  obs_times = seq(tau * n_dose_ss - 24, tau * n_dose_ss, by = 0.1),
  ii = tau, addl = n_dose_ss - 1
)
sol_ss <- rxode2::rxSolve(mod_typ, ev_ss, returnType = "data.frame", addDosing = FALSE)
#> ℹ omega/sigma items treated as zero: 'etalcl'
sol_ss <- sol_ss[!is.na(sol_ss$Cc) & sol_ss$time >= tau * n_dose_ss - 24, ]

trap_auc <- function(time, conc) sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)

auc_sim <- trap_auc(sol_ss$time, sol_ss$Cc)
auc_closed <- 1 * 3 * 400 / 162.721 * 1000   # frel * 3 doses * Dose / CL, mg/L -> ng/mL

c(simulated = auc_sim, closed_form = auc_closed,
  rel_diff = abs(auc_sim / auc_closed - 1))
#>    simulated  closed_form     rel_diff 
#> 7.374018e+03 7.374586e+03 7.702809e-05

stopifnot(abs(auc_sim / auc_closed - 1) < 1e-3)
```

The model is a true ODE system, not a `linCmt()` analytic solution that
would silently discard the explicit `d/dt` statements:

``` r

stopifnot(length(ui$linCmt) == 0L || !isTRUE(ui$linCmt))
stopifnot(identical(ui$state, c("depot", "central", "peripheral1")))
```

## Reproducing Hu 2025 Table 3: covariate effects on steady-state exposure

Table 3 is the paper’s own model-based simulation of Cmax,ss, Cmin,ss
and AUC0-24h,ss at 400 mg three times daily, expressed as ratios to the
reference condition. It is the single most informative published result
for validating this encoding, because it constrains eleven covariate
coefficients at once.

``` r

scenarios <- list(
  "Reference (healthy volunteer)"     = list(),
  "Friedreich ataxia (both effects)"  = list(DIS_FRDA = 1),
  "Itraconazole"                      = list(CONMED_ITRACONAZOLE = 1),
  "Rifampicin"                        = list(CONMED_RIFAMPICIN = 1),
  "Liquid PediaSure meal"             = list(FED_LIQUIDSUPP = 1),
  "Fasted"                            = list(FED = 0),
  "Body weight 17.14 kg (5th pctl)"   = list(WT = 17.14),
  "Body weight 96.66 kg (95th pctl)"  = list(WT = 96.66),
  "BMI 13.82 kg/m2 (5th pctl)"        = list(BMI = 13.82),
  "BMI 30.67 kg/m2 (95th pctl)"       = list(BMI = 30.67)
)

ss_metrics <- function(covariates, n_dose = n_dose_ss) {
  ev <- vatiq_events(
    id = 1, dose = 400,
    obs_times = seq(tau * n_dose - 24, tau * n_dose, by = 0.1),
    ii = tau, addl = n_dose - 1, covariates = covariates
  )
  s <- rxode2::rxSolve(mod_typ, ev, returnType = "data.frame", addDosing = FALSE)
  s <- s[!is.na(s$Cc) & s$time >= tau * n_dose - 24, ]
  c(auc = trap_auc(s$time, s$Cc), cmax = max(s$Cc), cmin = min(s$Cc))
}

ss_raw <- vapply(
  scenarios,
  function(x) ss_metrics(utils::modifyList(ref_cov, x)),
  numeric(3)
)
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
ss_ratio <- sweep(ss_raw, 1, ss_raw[, 1], "/")
```

Hu 2025 Table 3 reports mean ratios with 90% confidence intervals from a
smoothed parametric bootstrap. The simulated ratios below are single
typical-value predictions, so they are compared against the published
mean and checked for containment in the published interval.

| Scenario | Cmax,ss (sim) | Cmax,ss (Hu 2025) | Cmin,ss (sim) | Cmin,ss (Hu 2025) | AUC0-24,ss (sim) | AUC0-24,ss (Hu 2025) | AUC 90% CI low | AUC 90% CI high | In published CI |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| Reference (healthy volunteer) | 1.000 | 1.00 | 1.000 | 1.00 | 1.000 | 1.00 | 1.00 | 1.00 | yes |
| Friedreich ataxia (both effects) | 0.781 | 0.76 | 0.922 | 0.89 | 0.840 | 0.81 | 0.73 | 0.91 | yes |
| Itraconazole | 3.572 | 3.04 | 5.102 | 4.15 | 4.246 | 3.52 | 2.87 | 4.27 | yes |
| Rifampicin | 0.546 | 0.55 | 0.430 | 0.44 | 0.495 | 0.50 | 0.43 | 0.58 | yes |
| Liquid PediaSure meal | 0.069 | 0.07 | 0.069 | 0.07 | 0.069 | 0.07 | 0.06 | 0.09 | yes |
| Fasted | 0.036 | 0.04 | 0.036 | 0.04 | 0.036 | 0.04 | 0.03 | 0.05 | yes |
| Body weight 17.14 kg (5th pctl) | 3.177 | 2.88 | 3.707 | 3.16 | 3.386 | 2.99 | 2.37 | 3.76 | yes |
| Body weight 96.66 kg (95th pctl) | 0.701 | 0.71 | 0.690 | 0.70 | 0.696 | 0.70 | 0.65 | 0.76 | yes |
| BMI 13.82 kg/m2 (5th pctl) | 0.689 | 0.70 | 0.591 | 0.60 | 0.647 | 0.65 | 0.55 | 0.78 | yes |
| BMI 30.67 kg/m2 (95th pctl) | 1.341 | 1.33 | 1.509 | 1.47 | 1.408 | 1.38 | 1.21 | 1.58 | yes |

Simulated vs published steady-state exposure ratios (Hu 2025 Table 3).
{.table}

``` r

# Every simulated AUC ratio must lie inside the published 90% confidence interval.
stopifnot(all(tbl3$sim_auc >= tbl3$auc_lo), all(tbl3$sim_auc <= tbl3$auc_hi))

# The two purely-bioavailability effects and the induction effect are pure
# scalar multipliers on the whole profile, so they are reproduced essentially
# exactly rather than merely within the interval.
exact_rows <- tbl3$Scenario %in%
  c("Liquid PediaSure meal", "Fasted", "Rifampicin", "Body weight 96.66 kg (95th pctl)",
    "BMI 13.82 kg/m2 (5th pctl)")
stopifnot(max(abs(tbl3$sim_auc[exact_rows] - tbl3$pub_auc[exact_rows])) < 0.02)
```

### Why the identical Cmax / Cmin / AUC ratios settle the FK0 encoding

Hu 2025 writes the prandial and disease-on-bioavailability effects on a
single parameter FK0, which the paper describes both as “absorption
fraction via the first-order absorption process” (Table 2 footnote) and
as “relative bioavailability” (Section 2.3). Those are different
quantities, and only one reading reproduces the paper’s own numbers.

If the covariates moved only the *split* between the two absorption
arms, total exposure could not change at all: both arms deliver to the
same central compartment, so their fractions would still sum to one and
AUC would be identical. Table 3 instead reports the FK0-mediated ratios
as 0.50 / 0.50 / 0.50 (Friedreich ataxia), 0.07 / 0.07 / 0.07 (liquid
PediaSure) and 0.04 / 0.04 / 0.04 (fasted). Three different exposure
metrics moving by an identical factor is the signature of a pure scalar
multiplier applied to the whole concentration-time curve, which happens
only when the multiplier scales **both** arms and leaves their ratio
untouched.

The model file therefore separates the two jobs: `logitffo` carries the
estimated 74.4% / 25.6% split, and `lfdepot` carries the relative
bioavailability that the covariates act on. The check below confirms
that the split is genuinely untouched by the prandial covariates while
total exposure scales by exactly the published factor.

``` r

shape_ratio <- ss_ratio["cmax", ] / ss_ratio["auc", ]
prandial <- c("Liquid PediaSure meal", "Fasted", "Friedreich ataxia (both effects)")

# For the two pure-bioavailability scenarios the profile SHAPE is unchanged,
# so Cmax, Cmin and AUC all scale by the same factor (shape ratio == 1).
stopifnot(all(abs(shape_ratio[c("Liquid PediaSure meal", "Fasted")] - 1) < 1e-6))

# Friedreich ataxia also moves clearance, so its shape ratio is NOT 1.
stopifnot(abs(shape_ratio[["Friedreich ataxia (both effects)"]] - 1) > 0.01)

round(shape_ratio[prandial], 6)
#>            Liquid PediaSure meal                           Fasted 
#>                          1.00000                          1.00000 
#> Friedreich ataxia (both effects) 
#>                          0.92976
```

### The two scenarios that sit above the published mean

Itraconazole (simulated AUC ratio 4.25 against a published 3.52) and the
17.14 kg body weight (3.39 against 2.99) are the two rows furthest from
the published mean, though both remain inside the published 90%
interval. They are also precisely the two scenarios that *lower*
clearance and therefore *lengthen* the terminal half-life: from about 70
h at the reference to roughly 139 h with itraconazole and 120 h at 17.14
kg.

That suggests the published simulations were run over a study-like
dosing period rather than to true mathematical steady state, which would
systematically understate the exposure of the slowest-clearing arms.
Re-running the same two ratios over shorter horizons tests that
directly.

``` r

horizons <- c(21, 42, 90, 300)
horizon_tbl <- lapply(horizons, function(nd) {
  base <- ss_metrics(ref_cov, n_dose = nd)[["auc"]]
  data.frame(
    Doses = nd,
    Days = round(nd * tau / 24),
    Itraconazole = round(ss_metrics(utils::modifyList(ref_cov, list(CONMED_ITRACONAZOLE = 1)), nd)[["auc"]] / base, 2),
    `WT 17.14 kg` = round(ss_metrics(utils::modifyList(ref_cov, list(WT = 17.14)), nd)[["auc"]] / base, 2),
    check.names = FALSE
  )
}) %>% bind_rows()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'

knitr::kable(horizon_tbl, caption = "Simulated AUC0-24 ratio versus length of the dosing run.")
```

| Doses | Days | Itraconazole | WT 17.14 kg |
|------:|-----:|-------------:|------------:|
|    21 |    7 |         3.17 |        2.75 |
|    42 |   14 |         3.73 |        3.11 |
|    90 |   30 |         4.16 |        3.35 |
|   300 |  100 |         4.25 |        3.39 |

Simulated AUC0-24 ratio versus length of the dosing run. {.table}

The published means of 3.52 and 2.99 fall between the 7-day and 14-day
rows, which is where a simulation replicating an actual trial’s dosing
duration would land. The ratios rise monotonically towards the true
steady-state values as the run lengthens, so the difference is an
artefact of the simulation horizon rather than a discrepancy in the
encoded coefficients. No parameter was adjusted.

``` r

stopifnot(all(diff(horizon_tbl$Itraconazole) > 0))
stopifnot(all(diff(horizon_tbl$`WT 17.14 kg`) > 0))
# The published means are bracketed by the 7-day and 14-day horizons.
stopifnot(horizon_tbl$Itraconazole[1] < 3.52, horizon_tbl$Itraconazole[2] > 3.52)
stopifnot(horizon_tbl$`WT 17.14 kg`[1] < 2.99, horizon_tbl$`WT 17.14 kg`[2] > 2.99)
```

## Virtual cohort and PKNCA validation

Supplementary Table S3 reports a published dose-linearity assessment:
power coefficients from regressing log(Cmax) and log(AUC0-24h) on
log(dose) for 98 participants dosed once over 200 to 1400 mg. That is
the paper’s own non-compartmental result and the natural PKNCA target.

The virtual cohort matches the analysis-population demographics from
Supplementary Tables S1 and S2: body weight and BMI drawn to reproduce
the reported medians and ranges, roughly half female (sex is not a model
covariate), and healthy-volunteer disease status with a medium-fat meal,
which is the condition under which the single-dose linearity data were
collected.

``` r

rxode2::rxSetSeed(20250906)
set.seed(20250906)

doses <- c(200, 300, 400, 1400)
n_per_arm <- 50          # 50 per arm, well under the 200-per-arm cap
obs_grid <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24)

# Adult healthy-volunteer demographics (Supplementary Table S2, EPI743/CNS/NEU
# healthy-volunteer studies): body weight median ~73 kg, BMI median ~26.
cohort <- expand.grid(arm = seq_len(n_per_arm), dose = doses) %>%
  mutate(
    id = row_number(),
    WT = round(pmin(pmax(rlnorm(n(), log(73), 0.18), 50), 115), 1),
    BMI = round(pmin(pmax(rnorm(n(), 26, 3.2), 18), 34), 1)
  )

ev_cohort <- lapply(seq_len(nrow(cohort)), function(i) {
  vatiq_events(
    id = cohort$id[i], dose = cohort$dose[i], obs_times = obs_grid,
    covariates = utils::modifyList(ref_cov, list(WT = cohort$WT[i], BMI = cohort$BMI[i]))
  )
}) %>% bind_rows()

sim_cohort <- rxode2::rxSolve(
  readModelDb("Hu_2025_vatiquinone"), ev_cohort,
  returnType = "data.frame", addDosing = FALSE
) %>%
  left_join(distinct(cohort, id, dose), by = "id")

stopifnot(!anyNA(sim_cohort$sim), all(sim_cohort$sim > 0))
nrow(sim_cohort)
#> [1] 3400
```

`sim` carries the residual error and is the analogue of an observed
concentration; `Cc` is the individual prediction without residual error.
The dose-linearity regression is run on `sim` because Supplementary
Table S3 was computed from observed concentrations.

``` r

# `conc` here holds `sim`, i.e. the individual prediction WITH residual error,
# because Supplementary Table S3 was computed from observed concentrations.
conc_df <- sim_cohort %>%
  filter(!is.na(sim)) %>%
  transmute(id, dose_level = factor(dose), time, conc = sim)

# Time-zero anchor: the grid already contains time 0, but add it defensively so
# PKNCA never reports "AUC range starting before the first measurement".
conc_df <- conc_df %>%
  bind_rows(conc_df %>% distinct(id, dose_level) %>% mutate(time = 0, conc = 0)) %>%
  distinct(id, dose_level, time, .keep_all = TRUE) %>%
  arrange(id, time)

dose_df <- cohort %>% transmute(id, dose_level = factor(dose), amt = dose, time = 0)

conc_obj <- PKNCA::PKNCAconc(conc_df, conc ~ time | dose_level + id,
                             concu = "ng/mL", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | dose_level + id, doseu = "mg")

intervals <- data.frame(
  start = 0, end = 24, cmax = TRUE, tmax = TRUE, auclast = TRUE
)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
nca_wide <- as.data.frame(nca_res) %>%
  select(dose_level, id, PPTESTCD, PPORRES) %>%
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES) %>%
  mutate(dose = as.numeric(as.character(dose_level)))

stopifnot(nrow(nca_wide) == nrow(cohort), !anyNA(nca_wide$cmax), !anyNA(nca_wide$auclast))
summary(nca_wide[, c("cmax", "tmax", "auclast")])
#>       cmax               tmax           auclast       
#>  Min.   :   66.45   Min.   : 0.500   Min.   :  442.5  
#>  1st Qu.:  457.08   1st Qu.: 4.000   1st Qu.: 1476.2  
#>  Median :  921.50   Median : 5.500   Median : 2591.8  
#>  Mean   : 1450.08   Mean   : 5.525   Mean   : 4387.7  
#>  3rd Qu.: 1972.44   3rd Qu.: 7.000   3rd Qu.: 5534.1  
#>  Max.   :13254.86   Max.   :14.000   Max.   :34485.2
```

``` r

fit_cmax <- lm(log(cmax) ~ log(dose), data = nca_wide)
fit_auc <- lm(log(auclast) ~ log(dose), data = nca_wide)

linearity <- tibble::tibble(
  Parameter = c("Cmax", "AUC0-24h"),
  `Power coefficient (sim)` = round(c(coef(fit_cmax)[2], coef(fit_auc)[2]), 3),
  `Standard error (sim)` = round(c(summary(fit_cmax)$coefficients[2, 2],
                                   summary(fit_auc)$coefficients[2, 2]), 3),
  `Power coefficient (Hu 2025 Table S3)` = c(1.09, 1.04),
  `Standard error (Hu 2025)` = c(0.15, 0.097),
  `95% CI (Hu 2025)` = c("0.80, 1.39", "0.85, 1.23")
)
knitr::kable(linearity, caption = "Dose proportionality: simulated vs Hu 2025 Supplementary Table S3.")
```

| Parameter | Power coefficient (sim) | Standard error (sim) | Power coefficient (Hu 2025 Table S3) | Standard error (Hu 2025) | 95% CI (Hu 2025) |
|:---|---:|---:|---:|---:|:---|
| Cmax | 0.863 | 0.073 | 1.09 | 0.150 | 0.80, 1.39 |
| AUC0-24h | 0.922 | 0.050 | 1.04 | 0.097 | 0.85, 1.23 |

Dose proportionality: simulated vs Hu 2025 Supplementary Table S3.
{.table}

``` r

sim_power <- c(coef(fit_cmax)[2], coef(fit_auc)[2])

# The model is linear, so the true power coefficient is exactly 1. Gate on the
# centre of the estimate, not on any single subject's extreme, and require the
# published confidence intervals to contain the simulated value.
stopifnot(all(abs(sim_power - 1) < 0.35))
stopifnot(sim_power[1] > 0.80, sim_power[1] < 1.39)
stopifnot(sim_power[2] > 0.85, sim_power[2] < 1.23)

# Underlying model exposure is exactly dose proportional. This must be checked
# on TYPICAL VALUES at a fixed covariate condition: taking the maximum over a
# random cohort would compare a different extreme subject in each dose arm and
# says nothing about linearity.
dn <- lapply(doses, function(d) {
  s <- rxode2::rxSolve(
    mod_typ, vatiq_events(id = 1, dose = d, obs_times = obs_grid),
    returnType = "data.frame", addDosing = FALSE
  )
  data.frame(dose = d, peak_dn = max(s$Cc[!is.na(s$Cc)]) / d)
}) %>% bind_rows()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> ℹ omega/sigma items treated as zero: 'etalcl'
stopifnot(diff(range(dn$peak_dn)) / mean(dn$peak_dn) < 1e-8)
```

## Concentration scale against the published VPC

Hu 2025 Figure 2 shows the visual predictive check on a log scale over
the whole pooled dataset. The observed and simulated medians sit between
roughly 200 and 500 ng/mL, with the 5th percentile near 10 to 50 and the
95th near 2,000 to 5,000. Since 400 mg three times daily is the modal
regimen (212 of 343 participants), the typical-value steady-state
profile at that regimen should fall inside the median band.

``` r

ss_range <- range(sol_ss$Cc)
c(min = ss_range[1], max = ss_range[2])
#>      min      max 
#> 231.5623 402.4004

stopifnot(ss_range[1] > 100, ss_range[2] < 900)
```

That wide observed spread is also the reason the residual error is so
large. The model carries a random effect on clearance only; with no
random effect on absorption in a compound whose exposure moves 28-fold
with the meal, and with one 300 mg study peaking near 150 ng/mL (Figure
1C) while another 400 mg study peaks nearer 1,500 ng/mL (Figure 1B), all
of that unexplained between-study variability lands in the residual. A
log-scale residual SD of 1.062 corresponds to roughly a 2.9-fold spread
per standard deviation, which is what Figure 2 shows.

    #> ℹ omega/sigma items treated as zero: 'etalcl'
    #> ℹ omega/sigma items treated as zero: 'etalcl'
    #> ℹ omega/sigma items treated as zero: 'etalcl'
    #> Warning in scale_y_log10(): log-10 transformation introduced infinite values.

![Typical-value single-dose profiles at 300 mg under the three prandial
conditions of study EPI743-12-001, the crossover that identifies both
prandial coefficients (compare Hu 2025 Figure
1C).](Hu_2025_vatiquinone_files/figure-html/profile-plot-1.png)

Typical-value single-dose profiles at 300 mg under the three prandial
conditions of study EPI743-12-001, the crossover that identifies both
prandial coefficients (compare Hu 2025 Figure 1C).

``` r

peaks <- prandial_prof %>% group_by(Condition) %>% summarise(cmax = max(Cc), .groups = "drop")
peak_ref <- peaks$cmax[peaks$Condition == "Medium-fat meal"]

# Hu 2025 abstract: a medium-fat meal raises exposure ~14-fold over liquid
# PediaSure and ~25-fold over the fasted state.
ratio_liquid <- peak_ref / peaks$cmax[peaks$Condition == "Liquid PediaSure"]
ratio_fasted <- peak_ref / peaks$cmax[peaks$Condition == "Fasted"]
round(c(vs_liquid = ratio_liquid, vs_fasted = ratio_fasted), 1)
#> vs_liquid vs_fasted 
#>      14.5      27.8

stopifnot(abs(ratio_liquid - exp(2.671)) < 0.01, abs(ratio_fasted - exp(3.324)) < 0.01)
```

The medium-fat meal raises the peak 14.5-fold over the liquid supplement
and 27.8-fold over the fasted state, against the “14-fold” and “25-fold”
quoted in Hu 2025 Section 2.4. The 27.8 against 25 gap is a rounding of
the published summary statement, not a coefficient difference:
`exp(3.324)` is 27.8 exactly.

## Assumptions and deviations

1.  **FK0 re-parameterised into a split plus a relative
    bioavailability.** Hu 2025 gives the single parameter FK0 two jobs –
    the 74.4% / 25.6% absorption-arm split and the carrier of the
    prandial and disease bioavailability effects. The model file
    separates them into `logitffo` and `lfdepot`. This is not a
    modelling choice but the only reading consistent with the paper’s
    own Table 3, as shown above: covariates that moved only the split
    could not change exposure at all. All eleven published ratios are
    reproduced within their 90% confidence intervals.

2.  **TLAG2 is set to zero.** Figure 5 labels a lag time for the
    zero-order absorption process, drawn earlier than TLAG1, but Table 2
    contains no TLAG2 row and the theta numbering 1 to 17 is gapless, so
    no value is reported anywhere in the article or its supplement. The
    zero-order input therefore starts at the dose time. No value was
    invented. The practical effect is that the model has no absolute
    absorption lag: concentrations begin rising immediately through the
    zero-order arm, whereas Hu 2025 Section 2.2 describes “a visible lag
    of 2-3 h”. Note that Figure 1D, cited as the evidence for that lag,
    comes from a study whose first post-dose sample is at 2 h, so the
    figure cannot in fact resolve an absolute lag.

3.  **Residual error read as an additive error on log-transformed
    concentration.** The Methods never write out an `$ERROR` block.
    Table 2’s “Additive residual” of 1.062 is numbered 9 inside the
    theta sequence and reported with a %RSE and a 90% CI in theta style,
    which is the NONMEM idiom `W = THETA(9); Y = LOG(IPRED) + W*EPS(1)`.
    A linear-scale additive residual of 1.062 ng/mL would be about 0.05%
    of the median observed concentration and cannot be reconciled with
    the width of the Figure 2 VPC, whose 5th-to-95th spread implies a
    total log-scale SD near 1.4. Encoded as `Cc ~ lnorm(expSd)`. Whether
    1.062 is the SD or the variance changes little numerically (1.062
    against 1.031); the log-versus-linear scale is the decision that
    matters.

4.  **Inter-individual variability on Ka and V is omitted.** The Section
    2.3 equations carry eta terms on Ka and on V, but Table 2 reports
    estimates for neither and the text never returns to them. Their
    magnitudes are simply unreported. They are omitted rather than
    written as `~ fixed(0)` because a zero-variance diagonal makes OMEGA
    singular and breaks the Cholesky sampler used by `rxSolve`. Only the
    reported IIV on CL/F (variance 0.191, printed %CV 45.88) is encoded.
    No IIV is reported on FK0, TK0, TLAG1, Q/F or V2/F.

5.  **The fasted indicator is carried on the complement of `FED`.** The
    canonical column is oriented fed = 1, whereas Hu 2025’s FST
    indicator is oriented fasted = 1, so `model()` applies the
    coefficient to `(1 - FED)`. The model’s prandial reference is the
    medium-fat meal (`FED = 1`, `FED_LIQUIDSUPP = 0`), not the fasted
    state.

6.  **Patients with other mitochondrial diseases share the
    healthy-volunteer reference level.** `DIS_FRDA = 0` pools 116
    healthy volunteers with 54 patients with other mitochondrial
    diseases. Hu 2025’s Discussion reports that the latter group’s
    exposures were not distinguishable from healthy volunteers, and
    Section 2.3 defines the indicator against healthy volunteers only,
    so no separate level is modelled.

7.  **Two dose records per administration are required.** Any simulation
    must supply a bolus to `depot` and a `rate = -2` record of the same
    amount to `central`; supplying only one halves the dose. This is a
    property of encoding parallel absorption in `rxode2`, not of the
    source model.

8.  **Body weight and BMI act on clearance simultaneously and in
    opposite directions** (exponents 0.915 and -0.975) and are
    correlated in the cohort, so neither exponent is interpretable as a
    marginal effect. Simulations that vary one while holding the other
    fixed – as Hu 2025 Table 3 and this vignette both do – are reporting
    a conditional, not a marginal, effect.

9.  **No absolute bioavailability data exist.** Every volume and
    clearance is an apparent (`/F`) value, and relative bioavailability
    is anchored at 1 at the reference prandial and disease condition.

10. **The prandial coefficients rest on 18 subjects.** Only study
    EPI743-12-001 varied the meal (6 fasted, 6 liquid, 6 medium-fat
    records among 343 participants). The two largest effects in the
    model are therefore the most weakly supported, a limitation the
    paper does not discuss.
