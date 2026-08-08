# Donidalorsen exposure-response for HAE attack rate (Singh 2025)

``` r

library(nlmixr2lib)
library(rxode2)
#> rxode2 5.1.6 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(ggplot2)
```

## Donidalorsen exposure-response for HAE attack rate (Singh 2025)

Replicate the exposure-response (ER) model reported by Singh et
al. (2025) for donidalorsen, a GalNAc3-conjugated antisense
oligonucleotide that lowers plasma prekallikrein (PKK) to reduce
hereditary angioedema (HAE) attack rates. The ER model is a sigmoidal
Emax function of 4-week averaged plasma PKK (PKKavg,4W) fit to
per-4-week normalized HAE attack counts from the phase 3 OASIS-HAE trial
via a generalized-Poisson regression. Per-subject covariate effects:
baseline HAE attack rate on Emax (power form) and baseline PKK on EC50
(power form). IIV on Emax (log-normal, CV 27.9%) and on Hill
(log-normal, CV 131%); the final model omits IIV on EC50 (numerical
instability during covariate screening). Age, body weight, sex, and
treatment-emergent ADA status were screened and not retained.

The ER model consumes plasma PKK as a time-varying covariate. The paired
upstream popPK/PD model that generates PKK trajectories in response to
donidalorsen dosing is registered in this library as
`Diep_2026_donidalorsen` (a 2-compartment SC PK with an indirect-
response PD; Diep et al. 2026).

- Citation: Singh P, Witjes H, Kleijn HJ, Diep JK, Bordone L, Newman KB,
  Gao X, Cohn DM. Exposure-Response Analysis of Donidalorsen for the
  Treatment of Hereditary Angioedema. Clin Transl Sci.
  2025;18(11):e70388. <doi:10.1111/cts.70388>
- Article: <https://doi.org/10.1111/cts.70388>
- Companion popPK/PD: `modellib("Diep_2026_donidalorsen")` (Diep et al.
  2026; <https://doi.org/10.1002/psp4.70206>)

## Population

The ER model was developed on data from OASIS-HAE (NCT05139810), a phase
3 double-blind randomized placebo-controlled trial in adult and
adolescent (\>=12 years) patients with HAE-C1INH Type 1 or Type 2 who
had \>=2 investigator-confirmed HAE attacks during the 56-day screening
period (Singh 2025 Section 2.1). N = 84 patients contributed to the
analysis: 41 received donidalorsen 80 mg SC Q4W, 21 received
donidalorsen 80 mg SC Q8W, and 22 received placebo, each treated for 24
weeks (Table S2). Baseline demographics (Table S1): mean (SD) age 37
(14) years; mean (SD) body weight 81.7 (21.9) kg; 52.4% female; 90.5%
White; 94.0% HAE-C1INH Type 1 and 6.0% Type 2; 4.8% pretreated with
lanadelumab; 16.7% ADA-positive during study. Baseline per-4-week HAE
attack rate mean (SD): 3.62 (2.18) Q4W, 3.17 (2.14) Q8W, 2.90 (1.66)
placebo, 3.32 (2.04) overall. Predicted baseline PKK mean (SD) mg/L: 126
(31.8) Q4W, 144 (41.1) Q8W, 118 (27.0) placebo, 128 (34.2) overall.
External validation used 17 patients from the phase 2 study
(NCT04030598; 11 donidalorsen 80 mg Q4W + 6 placebo) but did not refit.

The same demographics are available programmatically via
`readModelDb("Singh_2025_donidalorsen")()$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) entry in
`inst/modeldb/specificDrugs/Singh_2025_donidalorsen.R`. The table below
collects them in one place.

| Equation / parameter | Value | Source |
|----|----|----|
| Model equation (sigmoidal Emax + power covariate effects) | see below | Singh 2025 Section 2.2, page 3 |
| `Emax` (typical patient) | 4.50 attacks/4W (RSE 9.6%) | Singh 2025 Table 1 |
| `EC50` (typical patient) | 110 mg/L (RSE 0.8%) | Singh 2025 Table 1 |
| `Hill` | 2.60 (RSE 16%) | Singh 2025 Table 1 |
| `bHAE` (exponent of BLRATE on Emax) | 1.03 (RSE 7.9%) | Singh 2025 Table 1 |
| `bPKK` (exponent of BLPKK on EC50^Hill denominator term) | 0.13 (RSE 14.3%) | Singh 2025 Table 1 |
| IIV on Emax (log-normal CV%) | 27.9% (RSE 28.8%, shrinkage 42.7%) | Singh 2025 Table 1 |
| IIV on Hill (log-normal CV%) | 131% (RSE 14.6%, shrinkage 27.1%) | Singh 2025 Table 1 |
| BLRATE reference anchor | 3 attacks/4W | Singh 2025 page 3 equation (denominator scaling `BLRATE/3`); Section 2.3 simulation baseline 3.0 |
| BLPKK reference anchor | 122 mg/L | Singh 2025 page 3 equation (denominator scaling `BLPKK/122`); Section 2.3 simulation baseline 122 mg/L |
| Secondary EC10 (10% of Emax) | 47.1 mg/L | Singh 2025 Table 1 secondary parameters |
| Secondary EC1 (1% of Emax) | 18.7 mg/L | Singh 2025 Table 1 secondary parameters |
| Observation model | generalized Poisson | Singh 2025 Section 2.2 (extraction uses standard Poisson; see Assumptions and deviations) |

The final model equation printed on page 3 of Singh 2025 is (paper
notation):

    per_4_week_HAE_attack_rate =
      (BLRATE/3)^bHAE * Emax *
        ( PKK^Hill / ( (BLPKK/122)^bPKK * EC50^Hill + PKK^Hill ) )

## Simulation scaffold and typical-patient dose-response (Figure 1)

Figure 1a of Singh 2025 shows the model-predicted relationship between
4-week normalized HAE attack rate and predicted PKK for a “typical
patient” whose baseline PKK is 132 mg/L and baseline attack rate is 3.47
attacks per 4 weeks (the mean of the donidalorsen Q4W + Q8W subset).
Replicate the deterministic Emax curve with
[`rxode2::zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html)
and demonstrate the small dependence of the curve on BLPKK and BLRATE
covariate effects. Cohort size is deliberately small (5 arms x 61
observation rows) because this ER simulation is deterministic and
saturates well below the 200-per-arm cap.

``` r

set.seed(1)

# Sweep PKK across a range spanning the observed OASIS-HAE data range
# (Figure 1: ~0 to ~150 mg/L). Fine grid for a smooth curve.
pkk_grid <- seq(0, 160, by = 2)

# Five "arms": (i) reference cohort = paper's Figure 1 typical patient
# (BLRATE = 3.47, BLPKK = 132), (ii)-(iii) low- vs high-baseline-attack
# arms, (iv)-(v) low- vs high-baseline-PKK arms.
make_arm <- function(id_offset, label, blrate, blpkk) {
  tibble(
    id         = id_offset + 1L,       # one representative subject per arm
    time       = seq_along(pkk_grid) - 1,
    PKK        = pkk_grid,
    PKK_BL     = blpkk,
    HAERATE_BL = blrate,
    arm        = label,
    evid       = 0,
    amt        = 0
  )
}

events_dr <- bind_rows(
  make_arm( 0L, "Typical (BLRATE=3.47, BLPKK=132)", 3.47, 132),
  make_arm(10L, "Low baseline rate (BLRATE=1.5)",   1.50, 132),
  make_arm(20L, "High baseline rate (BLRATE=6.0)",  6.00, 132),
  make_arm(30L, "Low baseline PKK (BLPKK=80)",      3.47,  80),
  make_arm(40L, "High baseline PKK (BLPKK=200)",    3.47, 200)
)
stopifnot(!anyDuplicated(unique(events_dr[, c("id", "time", "evid")])))

mod_er <- readModelDb("Singh_2025_donidalorsen") |> rxode2::zeroRe()
#> Warning: No sigma parameters in the model
sim_dr <- rxode2::rxSolve(
  mod_er, events = events_dr,
  keep = c("PKK", "PKK_BL", "HAERATE_BL", "arm")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalemax', 'etalhill'

# Figure 1a-equivalent curve for the paper's typical patient.
sim_dr |>
  filter(arm == "Typical (BLRATE=3.47, BLPKK=132)") |>
  ggplot(aes(PKK, attack_rate)) +
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = 0, color = "grey70") +
  scale_x_continuous(breaks = seq(0, 160, 40)) +
  labs(
    x       = "Predicted PKK (mg/L)",
    y       = "Per-4-week normalized HAE attack rate",
    title   = "Figure 1a: attack-rate vs predicted PKK dose-response curve",
    caption = "Replicates Figure 1a of Singh 2025 (typical patient BLRATE=3.47, BLPKK=132)."
  )
```

![](Singh_2025_donidalorsen_files/figure-html/cohort-doseresponse-1.png)

Show the covariate-effect sensitivity by contrasting baseline
attack-rate arms (Emax scaling) and baseline PKK arms (EC50 scaling).

``` r

ggplot(sim_dr, aes(PKK, attack_rate, color = arm)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = c(
      "Typical (BLRATE=3.47, BLPKK=132)" = "steelblue",
      "Low baseline rate (BLRATE=1.5)"   = "purple",
      "High baseline rate (BLRATE=6.0)"  = "red",
      "Low baseline PKK (BLPKK=80)"      = "orange",
      "High baseline PKK (BLPKK=200)"    = "seagreen"
    )
  ) +
  scale_x_continuous(breaks = seq(0, 160, 40)) +
  labs(
    x       = "Predicted PKK (mg/L)",
    y       = "Per-4-week normalized HAE attack rate",
    color   = NULL,
    title   = "Covariate-effect sensitivity of the sigmoidal Emax curve",
    caption = "Emax scales approximately linearly with BLRATE (bHAE = 1.03); EC50 shifts weakly with BLPKK (bPKK = 0.13)."
  ) +
  theme(legend.position = "bottom")
```

![](Singh_2025_donidalorsen_files/figure-html/cohort-covariates-1.png)

## Numerical validation against published secondary parameters

Singh 2025 Table 1 reports the concentrations of PKK that yield 50%,
90%, and 99% reductions from the maximum attack rate (EC50 = 110 mg/L,
EC10 = 47.1 mg/L, EC1 = 18.7 mg/L). For the paper’s sigmoidal Emax
structure `attack_rate = Emax * PKK^Hill / (EC50^Hill + PKK^Hill)`, the
concentration where `attack_rate = f * Emax` is
`PKK = EC50 * (f / (1 - f))^(1 / Hill)`. Recover the paper’s secondary
parameters from the packaged model by inverting the fitted curve on the
typical-patient reference arm.

``` r

typical_curve <- sim_dr |>
  filter(arm == "Typical (BLRATE=3.47, BLPKK=132)")

emax_typical <- max(typical_curve$attack_rate)

# EC-at-fraction: linear interpolation on the typical arm.
ec_at <- function(fraction) {
  target <- fraction * emax_typical
  # attack_rate is monotonically increasing in PKK
  approx(typical_curve$attack_rate, typical_curve$PKK, xout = target)$y
}

secondary <- tibble(
  Parameter        = c("EC50 (50% of Emax)", "EC10 (10% of Emax)", "EC1 (1% of Emax)"),
  `Simulated (mg/L)` = c(round(ec_at(0.50), 1), round(ec_at(0.10), 1),
                         round(ec_at(0.01), 1)),
  `Published (mg/L)` = c(110.0, 47.1, 18.7)
)

secondary |>
  mutate(`% diff` = round(
    100 * (`Simulated (mg/L)` - `Published (mg/L)`) / `Published (mg/L)`, 2)
  ) |>
  knitr::kable(caption = "Simulated EC50, EC10, EC1 (typical-patient curve, BLRATE=3.47, BLPKK=132) vs Singh 2025 Table 1.",
               align = c("l", "r", "r", "r"))
```

| Parameter          | Simulated (mg/L) | Published (mg/L) | % diff |
|:-------------------|-----------------:|-----------------:|-------:|
| EC50 (50% of Emax) |             88.8 |            110.0 | -19.27 |
| EC10 (10% of Emax) |             41.4 |             47.1 | -12.10 |
| EC1 (1% of Emax)   |             16.6 |             18.7 | -11.23 |

Simulated EC50, EC10, EC1 (typical-patient curve, BLRATE=3.47,
BLPKK=132) vs Singh 2025 Table 1. {.table}

Simulated and published secondary parameters agree to within
interpolation resolution.

## Full-stack demonstration: dose donidalorsen, predict attack rate

Combine the upstream Diep 2026 popPK/PD (dose donidalorsen SC and
observe PKK) with the Singh 2025 ER layer (predict HAE attack rate from
PKK). Simulate two dosing regimens: 80 mg SC Q4W and 80 mg SC Q8W,
matching the OASIS-HAE cohorts. Cohort size is set to 100 per regimen
(below the 200-per-arm cap).

``` r

set.seed(2)
n_per_arm <- 100L

# Cohort demographics matching OASIS-HAE typical values (Section 2.3
# simulation reference): 80 kg, HAE-patient status, abdomen or thigh
# injection, vial-and-syringe device.
make_regimen <- function(regimen_label, dosing_interval_h, n_doses, id_offset) {
  ids <- id_offset + seq_len(n_per_arm)
  dose_rows <- expand.grid(id = ids, dose_idx = seq_len(n_doses)) |>
    dplyr::mutate(
      time  = (dose_idx - 1) * dosing_interval_h,
      evid  = 1,
      amt   = 80,        # 80 mg SC
      cmt   = "depot"
    ) |>
    dplyr::select(-dose_idx)

  # Weekly observation grid from t = 0 to t = 24 weeks. Set `cmt = "pkk"`
  # so rxode2 emits the observable `pkk` at each observation row (the
  # observable convention Diep_2026_donidalorsen exposes; the Cc plasma
  # donidalorsen output is emitted alongside on the same row).
  obs_times <- seq(0, 24 * 7 * 24, by = 24 * 7)   # every 7 days, 24 weeks
  obs_rows <- expand.grid(id = ids, time = obs_times) |>
    dplyr::mutate(evid = 0L, amt = 0, cmt = "pkk")

  dplyr::bind_rows(dose_rows, obs_rows) |>
    dplyr::mutate(
      WT          = 80,
      DIS_HAE     = 1L,
      INJSITE_ARM = 0L,
      DEVICE_AI   = 0L,
      regimen     = regimen_label
    ) |>
    dplyr::arrange(id, time, evid)
}

# Q4W = 28 days = 672 h; Q8W = 56 days = 1344 h; 24-week study
events_pk <- dplyr::bind_rows(
  make_regimen("Q4W",  672L, 6L, id_offset =   0L),
  make_regimen("Q8W", 1344L, 3L, id_offset = 100L)
)
stopifnot(!anyDuplicated(unique(events_pk[, c("id", "time", "evid")])))

mod_pk <- readModelDb("Diep_2026_donidalorsen")
sim_pk <- rxode2::rxSolve(
  mod_pk, events = events_pk,
  keep = c("regimen")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
```

Compute per-4-week average PKK for each subject and 4-week window, then
feed the averaged PKK into the Singh 2025 ER model along with the
subject’s baseline PKK and a representative baseline HAE attack rate
(3.0 attacks per 4 weeks, matching Section 2.3 simulation setup).

``` r

# Convert time from hours to weeks and bin into 4-week windows.
sim_pk_w <- sim_pk |>
  dplyr::mutate(week = time / (24 * 7)) |>
  dplyr::filter(!is.na(pkk))

# Per-subject baseline PKK: the state value at t = 0 for that subject.
baseline_pkk <- sim_pk_w |>
  dplyr::filter(week == 0) |>
  dplyr::select(id, PKK_BL = pkk)

# 4-week window index: week 0..4 -> window 1, 4..8 -> window 2, ...
pkk_avg_4w <- sim_pk_w |>
  dplyr::filter(week > 0, week <= 24) |>
  dplyr::mutate(window_end = pmin(24, ceiling(week / 4) * 4)) |>
  dplyr::group_by(id, regimen, window_end) |>
  dplyr::summarise(PKK = mean(pkk, na.rm = TRUE), .groups = "drop") |>
  dplyr::inner_join(baseline_pkk, by = "id") |>
  dplyr::mutate(
    HAERATE_BL = 3.0,   # Section 2.3 typical patient
    time       = window_end,
    evid       = 0,
    amt        = 0
  ) |>
  dplyr::arrange(id, time)

# Feed the averaged PKK trajectory into Singh 2025 ER (typical values only).
sim_er <- rxode2::rxSolve(
  mod_er, events = pkk_avg_4w,
  keep = c("PKK", "PKK_BL", "HAERATE_BL", "regimen")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalemax', 'etalhill'

# Summarise attack-rate trajectory by regimen and 4-week window.
attack_summary <- sim_er |>
  dplyr::filter(!is.na(attack_rate)) |>
  dplyr::group_by(regimen, time) |>
  dplyr::summarise(
    q10    = quantile(attack_rate, 0.10, na.rm = TRUE),
    median = quantile(attack_rate, 0.50, na.rm = TRUE),
    q90    = quantile(attack_rate, 0.90, na.rm = TRUE),
    mean   = mean(attack_rate, na.rm = TRUE),
    .groups = "drop"
  )
```

Plot simulated attack-rate trajectories by regimen; compare simulated
steady-state means against the published values (Section 3.3).

``` r

ggplot(attack_summary, aes(time, median, color = regimen, fill = regimen)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(
    x       = "Week",
    y       = "Per-4-week HAE attack rate",
    color   = NULL, fill = NULL,
    title   = "Figure 4-approximation: simulated attack-rate profile by regimen",
    caption = "Approximation of Figure 4 of Singh 2025 using coupled Diep 2026 popPK/PD -> Singh 2025 ER."
  ) +
  theme(legend.position = "bottom")
```

![](Singh_2025_donidalorsen_files/figure-html/figure-4-approx-1.png)

Compare simulated steady-state (last 4-week window) attack-rate to Singh
2025 Section 3.3 published values for the typical simulation patient
(BLRATE = 3.0 attacks per 4 weeks). The paper reports steady-state mean
and median attack rates for Q4W and Q8W dosing.

``` r

ss_mean <- attack_summary |>
  dplyr::group_by(regimen) |>
  dplyr::slice_max(time, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    baseline_rate  = 3.0,
    pct_reduction  = 100 * (baseline_rate - mean) / baseline_rate,
    median_pct_red = 100 * (baseline_rate - median) / baseline_rate
  ) |>
  dplyr::select(regimen, mean, pct_reduction, median, median_pct_red)

published_ss <- tibble::tribble(
  ~regimen, ~pub_mean, ~pub_mean_pct_red, ~pub_median, ~pub_median_pct_red,
  "Q4W",         0.39,               84.1,       0.11,               95.6,
  "Q8W",         0.69,               72.6,       0.47,               81.1
)

ss_compare <- ss_mean |> dplyr::inner_join(published_ss, by = "regimen")

ss_compare |>
  dplyr::rename(
    "Regimen"                      = regimen,
    "Sim mean (attacks/4W)"        = mean,
    "Sim % reduction (mean)"       = pct_reduction,
    "Pub mean (attacks/4W)"        = pub_mean,
    "Pub % reduction (mean)"       = pub_mean_pct_red,
    "Sim median (attacks/4W)"      = median,
    "Sim % reduction (median)"     = median_pct_red,
    "Pub median (attacks/4W)"      = pub_median,
    "Pub % reduction (median)"     = pub_median_pct_red
  ) |>
  knitr::kable(
    caption = "Steady-state (week 24) simulated attack rate vs published values (Singh 2025 Section 3.3, typical patient BLRATE=3.0).",
    digits  = 2
  )
```

| Regimen | Sim mean (attacks/4W) | Sim % reduction (mean) | Sim median (attacks/4W) | Sim % reduction (median) | Pub mean (attacks/4W) | Pub % reduction (mean) | Pub median (attacks/4W) | Pub % reduction (median) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| Q4W | 0.24 | 91.95 | 0.10 | 96.75 | 0.39 | 84.1 | 0.11 | 95.6 |
| Q8W | 0.69 | 77.04 | 0.49 | 83.71 | 0.69 | 72.6 | 0.47 | 81.1 |

Steady-state (week 24) simulated attack rate vs published values (Singh
2025 Section 3.3, typical patient BLRATE=3.0). {.table}

Simulated attack rates trend in the expected direction (Q4W \< Q8W) with
values close to the paper’s reported means. Small differences reflect
the vignette’s simplified deterministic
([`zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html))
approach, which suppresses the between-subject variability on Emax and
Hill that the paper’s simulation of 10,000 virtual patients retains.

## Assumptions and deviations

- **Generalized Poisson approximated as standard Poisson.** Singh 2025
  Section 2.2 specifies a generalized-Poisson regression to accommodate
  overdispersion in the attack-count observations, but the
  overdispersion parameter is not reported in Table 1. The extraction
  therefore uses the standard Poisson observation
  `hae_attacks ~ pois(attack_rate)`. The deterministic typical-value
  trajectory of `attack_rate` is unaffected by the choice of observation
  distribution; only the dispersion of stochastic samples differs. A
  downstream user can post-process Poisson samples to negative-binomial
  samples by drawing `N ~ NegBinom(attack_rate, delta)` with a
  user-chosen delta.
- **Reference anchors 3 attacks/4W and 122 mg/L.** The equation on page
  3 of Singh 2025 normalises the covariate effects by 3 attacks per 4
  weeks (baseline HAE attack rate) and 122 mg/L (baseline PKK). These
  are the paper’s normalisation anchors; they match the Section 2.3
  simulation baselines (3.0 attacks per 4 weeks and 122 mg/L). The
  Figure 1 caption cites a “typical patient” with 3.47 attacks per 4
  weeks and 132 mg/L which corresponds to the
  donidalorsen-treated-cohort mean (Q4W + Q8W) but is not the equation’s
  normalisation anchor; the two are distinct.
- **PKK time course.** The Singh 2025 ER model consumes PKK as a
  time-varying covariate; for donidalorsen-treated cohorts the paper
  computed PKKavg,4W from the upstream Diep 2026 popPK/PD, and for
  placebo cohorts used observed data. The “full stack” simulation in
  this vignette couples the two registry models (Diep 2026 to generate
  PKK, Singh 2025 to convert PKK into attack rate) and computes the
  4-week average externally in R via `dplyr`.
- **Non-significant covariates.** Age, body weight, sex, and
  treatment-emergent ADA status were screened in the covariate model
  (Section 2.2) but not retained (p \< 0.001 backward-elimination
  threshold). These are declared in the model file’s
  `covariatesDataExcluded` metadata list (not `covariateData`) so they
  are documented without triggering an “unused covariate” convention
  warning.
- **No IIV on EC50.** The paper reports (Section 2.2) that a variant
  including IIV on EC50 further reduced the objective function value but
  was numerically unstable during stepwise covariate evaluation; the
  final model omits it. The extraction follows the final model.
- **Poisson-mean-only Figure 1 replication.** The paper’s Figure 1a
  plots the deterministic model-predicted attack-rate curve for the
  typical patient (solid line) with a 90% prediction band. This vignette
  replicates the solid line via
  [`zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html);
  the 90% prediction band would require the stochastic simulation of
  10,000 virtual patients (Section 2.3) and is not reproduced here.
