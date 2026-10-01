# Cimetidine (Hanke 2020)

## Model and source

Hanke 2020 is mainly a whole-body PBPK paper. It builds PK-Sim / MoBi
models of metformin and cimetidine and uses them to describe the
metformin *SLC22A2* 808G\>T drug-gene interaction, the
cimetidine-metformin drug-drug interaction, and metformin exposure in
chronic kidney disease. Alongside that work the authors ran a **NONMEM
population-PK analysis of cimetidine**. It is reported in full in
Electronic Supplementary Material (ESM) section 3, with the final
control stream in section 3.5. Its purpose was to describe the double
peaks that cimetidine plasma profiles show after oral dosing in the
fasted state, and so to derive the split-dose input used in the PBPK
simulations.

That population-PK analysis is what this package provides. The PK-Sim
whole-body PBPK layer is not reproduced. Its cell-to-plasma partition
coefficients are given only as a calculation method (“Diverse”,
“Calculated”, “R&R” in ESM Table S4.3.1; “PK-Sim” for metformin in Table
S2.3.1), and its organ volumes, blood flows and transporter expression
are PK-Sim database outputs that appear nowhere in the paper.

``` r

mod <- readModelDb("Hanke_2020_cimetidine")
```

- Citation: Hanke N, Turk D, Selzer D, Ishiguro N, Ebner T, Wiebe S,
  Muller F, Stopfer P, Nock V, Lehr T. A Comprehensive Whole-Body
  Physiologically Based Pharmacokinetic Drug-Drug-Gene Interaction Model
  of Metformin and Cimetidine in Healthy Adults and Renally Impaired
  Individuals. Clin Pharmacokinet. 2020;59(11):1419-1431.
  <doi:10.1007/s40262-020-00896-w>. Population-PK parameters from
  Electronic Supplementary Material Table S3.4.1; structure and
  residual-error model from the NONMEM control stream in ESM section
  3.5.
- Article: <https://doi.org/10.1007/s40262-020-00896-w>
- ESM 1 (model documentation, Tables S3.3.1-S3.4.2 and the NONMEM code):
  distributed with the article and available through Europe PMC
  (`PMC7658088`, supplementary file `40262_2020_896_MOESM1_ESM.pdf`).

Two-compartment population PK model for intravenous and fasted oral
cimetidine in adults (Hanke 2020, Electronic Supplementary Material
section 3), fitted in NONMEM to mean concentration-time profiles from 25
published studies (100-800 mg). The fasted-state double peak is
described by splitting each oral dose into two portions that share one
first-order absorption rate constant: the first portion (fraction 1 -
VF2, typical 71.2 percent) is absorbed from depot without delay and the
second (fraction VF2, typical 28.8 percent) from depot2 after a lag time
of 1.54 h. Total oral bioavailability is 90.2 percent and elimination is
first order from the central compartment. Random effects are
between-study variability on VF2 (logit scale), the second-portion lag
time, CL and Vc; there are no covariates. This is the population-PK
analysis the paper used to derive the split-dose input for its PK-Sim
whole-body PBPK model; the PBPK layer itself is not reproduced here.

## Population

ESM Table S3.3.1 lists the 25 published studies the model was fitted to.
Nine gave cimetidine intravenously (100-400 mg as a bolus or 2-30 min
infusion) and 16 gave it orally in the fasted state (200-800 mg as a
solution, capsule or tablet). Subjects were 19-80 years old. ESM Table
S4.2.1 marks each cohort as healthy volunteers or peptic ulcer patients.
The table has 215 subject-arm entries, and subjects overlap between the
intravenous and oral arms of the crossover studies.

The analysis was fitted to **study-mean** concentration-time profiles,
plus three published individual profiles (ESM section 3.3.1). The random
effects are therefore interstudy variability (ISV), not between-subject
variability (section 3.3.2). A simulated “subject” from this model
stands for one study’s mean profile. The model has no covariates.

## Source trace

| Element | Value | Source |
|----|----|----|
| Two-compartment disposition, first-order elimination from central | – | ESM section 3.4; `$DES` in section 3.5 |
| Oral dose split into two portions, shared `KA1`, lag on the second only | – | ESM section 3.4, Figure S3.4.1; `$PK`/`$DES` |
| `lka` | 0.753 1/h | Table S3.4.1 `KA` |
| `lfdepot` (FTOT) | 0.902 | Table S3.4.1 `FTOT` = 90.2% |
| `logitfrac` (VF2) | logit(0.288) | Table S3.4.1 `VF2`; `$PK` `PHI_2 = LOG(VF2/(1-VF2))` |
| `ltlag2` (ALAG2) | 1.54 h | Table S3.4.1 `ALAG` |
| `lcl` | 41.2 L/h | Table S3.4.1 `CL` |
| `lvc` | 32.6 L | Table S3.4.1 `V3` |
| `lq` | 45.4 L/h | Table S3.4.1 `Q` |
| `lvp` | 46 L | Table S3.4.1 `V4` |
| `etalogitfrac` | log(1 + 0.848^2) | Table S3.4.1 `ISV VF2` = 84.8 %CV (logit scale, `$PK` `VF2_2`) |
| `etaltlag2` | log(1 + 0.205^2) | Table S3.4.1 `ISV ALAG` = 20.5 %CV |
| `etalcl` | log(1 + 0.216^2) | Table S3.4.1 `ISV CL` = 21.6 %CV |
| `etalvc` | log(1 + 0.387^2) | Table S3.4.1 `ISV V3` = 38.7 %CV |
| `propSd` | 0.122 | Table S3.4.1 `Prop RE` = 12.2% |
| `addSd` | 0.017972 mg/L, fixed | `$SIGMA` `0.000323 FIX`, square-rooted; Table S3.4.1 `Add RE` 1.8 (fixed) |
| `Cc = central / vc` | mg/L | `$PK` `S3 = V3` |
| `Cc ~ add + prop` | – | `$ERROR` `Y = IPRED + IPRED*EPS(1) + EPS(2)` |

## Dosing convention

The NONMEM model splits an oral dose by bioavailability fractions:
`F1 = (1 - VF2) * FTOT` for `DEPOT1` and `F2 = VF2 * FTOT` for `DEPOT2`,
with `ALAG2` on `DEPOT2` only. An oral dose is therefore given as **two
dose records of the full dose at the same time**, one to `depot` and one
to `depot2`. An intravenous dose is a single record to `central`.

``` r

obs_times <- sort(unique(c(seq(0, 6, by = 0.05), seq(6.25, 12, by = 0.25), 14, 16, 20, 24)))

oral_events <- function(id, dose, treatment) {
  dplyr::bind_rows(
    data.frame(id = id, time = 0, amt = dose, evid = 1, cmt = c("depot", "depot2")),
    data.frame(id = id, time = obs_times, amt = 0, evid = 0, cmt = "central")
  ) |>
    dplyr::mutate(treatment = treatment)
}

iv_events <- function(id, dose, treatment) {
  dplyr::bind_rows(
    data.frame(id = id, time = 0, amt = dose, evid = 1, cmt = "central"),
    data.frame(id = id, time = obs_times, amt = 0, evid = 0, cmt = "central")
  ) |>
    dplyr::mutate(treatment = treatment)
}
```

## Typical-value simulation

``` r

mod_typ <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

typ_design <- tibble::tribble(
  ~treatment,       ~route, ~dose,
  "200 mg IV bolus", "iv",   200,
  "200 mg oral",     "oral", 200,
  "400 mg oral",     "oral", 400,
  "800 mg oral",     "oral", 800
)
typ_design$id <- seq_len(nrow(typ_design))

ev_typ <- dplyr::bind_rows(lapply(seq_len(nrow(typ_design)), function(i) {
  d <- typ_design[i, ]
  if (d$route == "iv") {
    iv_events(d$id, d$dose, d$treatment)
  } else {
    oral_events(d$id, d$dose, d$treatment)
  }
}))

sim_typ <- rxode2::rxSolve(mod_typ, ev_typ, keep = "treatment", returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalogitfrac', 'etaltlag2', 'etalcl', 'etalvc'
#> Warning: multi-subject simulation without without 'omega'
```

``` r

ggplot(sim_typ, aes(time, Cc * 1000, colour = treatment)) +
  geom_line() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 12)) +
  labs(x = "Time (h)", y = "Cimetidine (ng/mL)", colour = NULL)
#> Warning in scale_y_log10(): log-10 transformation introduced infinite values.
```

![Typical-value cimetidine profiles. The fasted oral profiles show the
second peak from the delayed dose
portion.](Hanke_2020_cimetidine_files/figure-html/typical-plot-1.png)

Typical-value cimetidine profiles. The fasted oral profiles show the
second peak from the delayed dose portion.

## Replicate the per-study fits

ESM Table S3.4.2 gives the study-specific (empirical Bayes) estimates of
the second-portion fraction VF2 and lag time ALAG for every fasted oral
study. The ESM does not print the study-specific CL or V3, so those stay
at their typical values here. Each study is simulated by writing the
matching `etalogitfrac` and `etaltlag2` values into its event rows and
solving the zero-random-effects model, so the data columns supply the
etas.

The observed Cmax and AUC come from ESM Table S4.5.2 (“Obs” columns).
Those are the observed study means the model was fitted to. The table’s
“Pred” columns are the PBPK model’s predictions and are not used here.
The two Burland 1975 single-subject arms and the Tiseo 1998
multiple-dose arm are left out: Burland reports one subject per arm, and
Tiseo’s AUC covers a multiple-dose interval. The three Bodemar 1979
tablet cohorts are also left out, because Table S4.5.2 does not say
which of its three rows matches which Table S3.4.2 row.

``` r

# ESM Table S3.4.2 (vf2, alag) joined to ESM Table S4.5.2 observed Cmax and
# AUClast (ng/mL, h*ng/mL), matched on dose, formulation and reference.
studies <- tibble::tribble(
  ~study,                         ~dose, ~vf2,  ~alag, ~cmax_obs, ~auclast_obs,
  "Bodemar 1981, 200 mg",         200,   0.426, 1.68,  868.72,    3706.71,
  "Kanto 1981, 200 mg",           200,   0.276, 2.46,  918.16,    3408.48,
  "Mihaly 1984, 200 mg",          200,   0.496, 1.39,  601.30,    2954.89,
  "Walkenstein 1978, 300 mg sol", 300,   0.148, 1.64,  1062.50,   4933.56,
  "D'Angio 1986, 300 mg",         300,   0.233, 1.93,  1723.10,   6172.41,
  "Walkenstein 1978, 300 mg tab", 300,   0.282, 1.66,  1430.00,   4855.28,
  "Bodemar 1981, 400 mg",         400,   0.481, 1.41,  1934.73,   7169.56,
  "Grahnen 1979, 400 mg",         400,   0.133, 1.53,  2220.10,   7603.51,
  "Bodemar 1981, 800 mg",         800,   0.393, 1.64,  3682.78,   14177.73
)
studies$id <- seq_len(nrow(studies))

ev_study <- dplyr::bind_rows(lapply(seq_len(nrow(studies)), function(i) {
  oral_events(studies$id[i], studies$dose[i], studies$study[i])
})) |>
  dplyr::mutate(
    etalogitfrac = qlogis(studies$vf2[id]) - qlogis(0.288),
    etaltlag2 = log(studies$alag[id] / 1.54),
    etalcl = 0,
    etalvc = 0
  )

sim_study <- rxode2::rxSolve(mod_typ, ev_study, keep = "treatment", returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalogitfrac', 'etaltlag2', 'etalcl', 'etalvc'
#> Warning: multi-subject simulation without without 'omega'
```

``` r

ggplot(sim_study, aes(time, Cc * 1000)) +
  geom_line(colour = "steelblue") +
  geom_hline(
    data = studies |> dplyr::rename(treatment = study),
    aes(yintercept = cmax_obs), linetype = "dashed", colour = "grey40"
  ) +
  facet_wrap(~treatment, scales = "free_y") +
  coord_cartesian(xlim = c(0, 8)) +
  labs(x = "Time (h)", y = "Cimetidine (ng/mL)")
```

![Study-specific profiles from the Table S3.4.2 estimates, with the
observed Cmax (points) from Table S4.5.2. Compare with ESM Figures
S3.4.3-S3.4.16.](Hanke_2020_cimetidine_files/figure-html/per-study-plot-1.png)

Study-specific profiles from the Table S3.4.2 estimates, with the
observed Cmax (points) from Table S4.5.2. Compare with ESM Figures
S3.4.3-S3.4.16.

The dashed line is the observed Cmax of each study.

## Stochastic simulation across studies

Here 200 virtual studies of a 400 mg fasted oral dose are drawn from the
interstudy variability. The band shows how widely study-mean profiles
are expected to vary, which is the variability the model describes.

``` r

n_studies <- 200
ev_vpc <- dplyr::bind_rows(lapply(seq_len(n_studies), function(i) {
  oral_events(i, 400, "400 mg oral")
}))
rxode2::rxSetSeed(20200525)
sim_vpc <- rxode2::rxSolve(mod, ev_vpc, keep = "treatment", returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'

vpc_band <- sim_vpc |>
  dplyr::group_by(time) |>
  dplyr::summarise(
    q05 = quantile(Cc, 0.05) * 1000,
    q50 = quantile(Cc, 0.50) * 1000,
    q95 = quantile(Cc, 0.95) * 1000,
    .groups = "drop"
  )
```

``` r

ggplot(vpc_band, aes(time)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = q50), colour = "steelblue") +
  coord_cartesian(xlim = c(0, 12)) +
  labs(x = "Time (h)", y = "Cimetidine (ng/mL)")
```

![Median and 90% interval of 200 simulated study-mean profiles after 400
mg cimetidine orally in the fasted state (IPRED, without residual
error).](Hanke_2020_cimetidine_files/figure-html/vpc-plot-1.png)

Median and 90% interval of 200 simulated study-mean profiles after 400
mg cimetidine orally in the fasted state (IPRED, without residual
error).

## PKNCA validation

``` r

nca_for <- function(sim, events) {
  conc <- sim |>
    dplyr::filter(!is.na(Cc)) |>
    dplyr::mutate(Cc = Cc * 1000) |>
    dplyr::select(id, time, Cc, treatment)
  # One dose row per id: for oral dosing the depot and depot2 records carry
  # the same full dose.
  dose <- events |>
    dplyr::filter(evid == 1) |>
    dplyr::distinct(id, time, amt, treatment)
  PKNCA::pk.nca(PKNCA::PKNCAdata(
    PKNCA::PKNCAconc(conc, Cc ~ time | treatment + id),
    PKNCA::PKNCAdose(dose, amt ~ time | treatment + id),
    intervals = data.frame(
      start = 0, end = Inf,
      cmax = TRUE, tmax = TRUE, aucinf.obs = TRUE, half.life = TRUE
    )
  ))
}

nca_typ <- as.data.frame(nca_for(sim_typ, ev_typ)) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "aucinf.obs", "half.life")) |>
  dplyr::select(treatment, PPTESTCD, PPORRES)

nca_study <- as.data.frame(nca_for(sim_study, ev_study)) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "aucinf.obs")) |>
  dplyr::select(treatment, PPTESTCD, PPORRES)
```

### Structural checks

The typical-value AUC0-inf has a closed form, `Dose / CL` intravenously
and `FTOT * Dose / CL` orally, whatever the split fraction and lag. The
terminal half-life is `ln 2 / beta`, where `beta` is the smaller
eigenvalue of the two-compartment system. Both sides use the same
parameters, so the only difference is numerical error in the solve and
the NCA, and a tight bound is appropriate.

``` r

p <- list(cl = 41.2, vc = 32.6, q = 45.4, vp = 46, f = 0.902)
k10 <- p$cl / p$vc
k12 <- p$q / p$vc
k21 <- p$q / p$vp
beta <- ((k10 + k12 + k21) - sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21)) / 2

expected <- typ_design |>
  dplyr::mutate(
    aucinf_expected = ifelse(route == "iv", 1, p$f) * dose / p$cl * 1000,
    half_life_expected = log(2) / beta
  )

check_typ <- nca_typ |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
  dplyr::left_join(expected, by = "treatment")

stopifnot(
  all(abs(check_typ$aucinf.obs / check_typ$aucinf_expected - 1) < 0.01),
  all(abs(check_typ$half.life / check_typ$half_life_expected - 1) < 0.02)
)

# Every fasted oral study profile has two peaks: the immediate portion and
# the delayed portion released after ALAG. These are deterministic
# (typical-value) solves, so the count is exact.
n_peaks <- sim_study |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(n_peaks = sum(diff(sign(diff(Cc))) < 0), .groups = "drop")
stopifnot(all(n_peaks$n_peaks == 2))

knitr::kable(
  check_typ |>
    dplyr::select(treatment, cmax, tmax, aucinf.obs, aucinf_expected, half.life, half_life_expected) |>
    dplyr::rename(
      "Treatment" = treatment,
      "Cmax (ng/mL)" = cmax,
      "Tmax (h)" = tmax,
      "AUC0-inf (ng*h/mL)" = aucinf.obs,
      "Closed-form AUC0-inf" = aucinf_expected,
      "t1/2 (h)" = half.life,
      "Closed-form t1/2" = half_life_expected
    ),
  digits = 2,
  caption = "Typical-value NCA against the closed forms."
)
```

| Treatment | Cmax (ng/mL) | Tmax (h) | AUC0-inf (ng\*h/mL) | Closed-form AUC0-inf | t1/2 (h) | Closed-form t1/2 |
|:---|---:|---:|---:|---:|---:|---:|
| 200 mg IV bolus | 6134.97 | 0.00 | 4855.37 | 4854.37 | 1.81 | 1.81 |
| 200 mg oral | 889.12 | 2.05 | 4377.93 | 4378.64 | 1.83 | 1.81 |
| 400 mg oral | 1778.24 | 2.05 | 8755.87 | 8757.28 | 1.83 | 1.81 |
| 800 mg oral | 3556.47 | 2.05 | 17511.73 | 17514.56 | 1.83 | 1.81 |

Typical-value NCA against the closed forms. {.table}

### Comparison against published NCA

The comparison is against the observed study means in ESM Table S4.5.2.
The table’s AUC is AUC from dosing to the **last measured
concentration** (AUClast, ESM section 1.3), and it does not give the
last sampling times, so the simulated AUC0-inf is expected to run
somewhat above it.

``` r

published <- studies |>
  dplyr::transmute(treatment = study, cmax = cmax_obs, aucinf.obs = auclast_obs)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_study,
  reference = published,
  by = "treatment",
  units = c(cmax = "ng/mL", tmax = "h", aucinf.obs = "ng*h/mL"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = "Study-specific simulated NCA vs. the observed means in Hanke 2020 ESM Table S4.5.2 (reference AUC is AUClast). * differs from the reference by more than 20%.",
  align = c("l", "l", "r", "r", "r")
)
```

| NCA parameter | treatment | Reference | Simulated | % diff |
|:---|:---|---:|---:|---:|
| Cmax (ng/mL) | Bodemar 1981, 200 mg | 869 | 895 | +3.0% |
| Cmax (ng/mL) | Kanto 1981, 200 mg | 918 | 776 | -15.5% |
| Cmax (ng/mL) | Mihaly 1984, 200 mg | 601 | 947 | +57.6%\* |
| Cmax (ng/mL) | Walkenstein 1978, 300 mg sol | 1060 | 1370 | +28.9%\* |
| Cmax (ng/mL) | D’Angio 1986, 300 mg | 1720 | 1230 | -28.5%\* |
| Cmax (ng/mL) | Walkenstein 1978, 300 mg tab | 1430 | 1300 | -9.0% |
| Cmax (ng/mL) | Bodemar 1981, 400 mg | 1930 | 1880 | -2.6% |
| Cmax (ng/mL) | Grahnen 1979, 400 mg | 2220 | 1860 | -16.3% |
| Cmax (ng/mL) | Bodemar 1981, 800 mg | 3680 | 3570 | -3.0% |
| AUC0-∞ (obs) (ng\*h/mL) | Bodemar 1981, 200 mg | 3710 | 4380 | +18.1% |
| AUC0-∞ (obs) (ng\*h/mL) | Kanto 1981, 200 mg | 3410 | 4380 | +28.4%\* |
| AUC0-∞ (obs) (ng\*h/mL) | Mihaly 1984, 200 mg | 2950 | 4380 | +48.2%\* |
| AUC0-∞ (obs) (ng\*h/mL) | Walkenstein 1978, 300 mg sol | 4930 | 6570 | +33.1%\* |
| AUC0-∞ (obs) (ng\*h/mL) | D’Angio 1986, 300 mg | 6170 | 6570 | +6.4% |
| AUC0-∞ (obs) (ng\*h/mL) | Walkenstein 1978, 300 mg tab | 4860 | 6570 | +35.3%\* |
| AUC0-∞ (obs) (ng\*h/mL) | Bodemar 1981, 400 mg | 7170 | 8760 | +22.1%\* |
| AUC0-∞ (obs) (ng\*h/mL) | Grahnen 1979, 400 mg | 7600 | 8760 | +15.2% |
| AUC0-∞ (obs) (ng\*h/mL) | Bodemar 1981, 800 mg | 14200 | 17500 | +23.5%\* |

Study-specific simulated NCA vs. the observed means in Hanke 2020 ESM
Table S4.5.2 (reference AUC is AUClast). \* differs from the reference
by more than 20%. {.table style="width:100%;"}

``` r

attr(cmp, "footnote")
#> [1] "* differs from reference by more than ±20%."
```

``` r

ratios <- nca_study |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
  dplyr::left_join(published, by = "treatment", suffix = c("_sim", "_obs")) |>
  dplyr::mutate(
    r_cmax = cmax_sim / cmax_obs,
    r_auc = aucinf.obs_sim / aucinf.obs_obs
  )

# Intravenous studies: typical AUC0-inf = Dose / CL against the observed
# AUClast of the seven plasma-sampled intravenous arms in Table S4.5.2
# (Walkenstein 1978 sampled whole blood and is left out).
iv_obs <- tibble::tribble(
  ~study,             ~dose, ~auclast_obs,
  "Grahnen 1979",     100,   1691.08,
  "Bodemar 1981",     200,   5179.28,
  "Mihaly 1984",      200,   3781.58,
  "Larsson 1982",     200,   4383.99,
  "Morgan 1983 5 min", 200,  3902.27,
  "Morgan 1983 30 min", 200, 4222.39,
  "Lebert 1981",      300,   6049.42
) |>
  dplyr::mutate(r_auc = (dose / 41.2 * 1000) / auclast_obs)

# Bioavailability implied by the observed data alone: dose-normalised
# observed AUClast, oral over intravenous (medians over studies). It needs no
# model output, so it checks the transcribed FTOT against the data it was
# fitted to; the AUClast truncation largely cancels in the ratio.
f_obs <- median(studies$auclast_obs / studies$dose) /
  median(iv_obs$auclast_obs / iv_obs$dose)

gate <- tibble::tibble(
  Check = c(
    "Oral Cmax, median simulated/observed",
    "Oral AUC0-inf / observed AUClast, median",
    "IV AUC0-inf / observed AUClast, median",
    "Observed oral/IV dose-normalised AUC (vs FTOT = 0.902)"
  ),
  Value = c(median(ratios$r_cmax), median(ratios$r_auc), median(iv_obs$r_auc), f_obs)
)
knitr::kable(gate, digits = 2, caption = "Summary of the per-study comparison.")
```

| Check                                                  | Value |
|:-------------------------------------------------------|------:|
| Oral Cmax, median simulated/observed                   |  0.97 |
| Oral AUC0-inf / observed AUClast, median               |  1.24 |
| IV AUC0-inf / observed AUClast, median                 |  1.20 |
| Observed oral/IV dose-normalised AUC (vs FTOT = 0.902) |  0.88 |

Summary of the per-study comparison. {.table}

``` r


stopifnot(
  # Cmax depends on ka, the split, the lag, Vc and the peripheral exchange:
  # a mis-transcribed volume or rate constant moves it by tens of percent.
  abs(median(ratios$r_cmax) - 1) < 0.15,
  # AUC0-inf must not fall below the observed AUClast in the median, and a
  # clearance or bioavailability error of a factor 1.5 would breach the cap.
  median(ratios$r_auc) > 1,
  median(ratios$r_auc) < 1.5,
  median(iv_obs$r_auc) > 1,
  median(iv_obs$r_auc) < 1.5,
  # The data-implied bioavailability agrees with the fitted FTOT.
  abs(f_obs / 0.902 - 1) < 0.15
)
```

Study-specific **Cmax** agrees closely with the observed means (median
ratio about 0.97), even though the study-specific CL and V3 are held at
their typical values. Individual studies differ by up to about 60%.
Mihaly 1984 differs most: its observed peak of 601 ng/mL is low for 200
mg. Part of that spread belongs to the study-level V3 variability (38.7
%CV), which the ESM does not print per study.

**AUC** runs above the observed AUClast by a similar amount on both
routes: a median of about 24% orally and 20% intravenously. Because the
gap is the same size with and without absorption, it lies in the
disposition and not in the split-dose input. Part of it is the
AUClast-versus-AUC0-inf difference, since the reference AUC stops at the
last sample. The rest is the between-study clearance variability: CL has
21.6 %CV, while every simulated study here uses the typical 41.2 L/h.
The observed data alone imply a bioavailability of about 0.88
(dose-normalised oral over intravenous AUClast), which agrees with the
fitted FTOT of 0.902. The maintainers did not adjust any parameter to
improve these comparisons.

## Assumptions and deviations

- **Only the population-PK layer of the paper is provided.** The PK-Sim
  / MoBi whole-body PBPK models of metformin and cimetidine are not
  reproduced. Their partition coefficients, organ physiology and
  transporter expression are platform database outputs that the paper
  does not tabulate. The cimetidine-metformin DDI and the
  renal-impairment scaling live entirely in the PBPK layer.
- **%CV to variance.** Table S3.4.1 reports each ISV as a %CV. The
  variances here use `omega^2 = log(1 + CV^2)`, the log-normal relation
  that the ESM itself uses elsewhere (Table S9.0.1 footnote: “35 % CV …
  (= 1.40 GSD)”). The same column and conversion are applied to the
  logit-scale VF2 eta. If the authors instead reported
  `100 * sqrt(omega^2)` for that row, the variance would be 0.719 rather
  than 0.542. The spread of the study-specific VF2 estimates in Table
  S3.4.2 (SD of the logit-scale deviations about 0.80) cannot separate
  the two readings. The `$OMEGA` values printed in the control stream
  are initial estimates and were not used.
- **Additive residual error.** The control stream fixes the additive
  `$SIGMA` at 0.000323, i.e. an SD of 0.018 in the model’s concentration
  unit. Table S3.4.1 prints it as “1.8”, scaled by 100 like its
  percent-valued proportional row.
- **Concentration unit.** The control stream scales the central amount
  by `S3 = V3` with volumes in litres. With doses in mg this gives mg/L,
  which is also the only unit in which a 0.018 additive SD is plausible.
  The vignette multiplies by 1000 to compare with the ng/mL values of
  Table S4.5.2.
- **ALAG uncertainty.** Table S3.4.1 prints an RSE of 0% for ALAG. The
  control stream does not fix it (`$THETA (0, 1.35)`, no FIX flag) and
  the text says all parameters were estimated, so it is treated as
  estimated.
- **Dose records.** An oral dose must be given as two records of the
  full dose, to `depot` and `depot2`. A single record to `depot` alone
  would deliver only `(1 - VF2) * FTOT` of the dose.
- **Matrix.** Some source studies measured whole blood rather than
  plasma (Table S4.5.2). The population-PK analysis pooled them without
  a matrix term, and the model output is labelled plasma.
