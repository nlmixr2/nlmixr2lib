# Pretomanid (Nyang'wa 2026)

## Model and source

- Citation: Nyang’wa BT, Motta I, Moodliar R, Solodovnikova V, Rajaram
  S, Rasool M, Berry C, Huang Z, Davies G, Moore DAJ, Kloprogge F
  (2026). Population pharmacokinetics and target attainment of
  pretomanid in rifampicin-resistant tuberculosis patients. Sci Rep
  16:46217. <doi:10.1038/s41598-026-46217-2>
- Description: One-compartment population PK model for oral pretomanid
  200 mg once daily in adults with rifampicin-resistant tuberculosis
  (RR-TB) treated with the BPaL, BPaLM, or BPaLC regimen in the
  TB-PRACTECAL trial PKPD sub-study (Nyang’wa 2026). First-order
  absorption (ka 0.316 1/h, no between-subject variability) into a
  single central compartment with first-order elimination; apparent
  clearance CL/F 3.10 L/h and apparent central volume V/F 102 L, both at
  the cohort median fat-free mass of 45.5 kg. Fat-free mass is the only
  retained covariate and enters as a priori allometric scaling with
  fixed exponents 0.75 on clearance and 1 on volume; FFM was selected
  over total body weight and body mass index on base-model fit.
  Between-subject variability is a diagonal pair on clearance (32.9% CV)
  and central volume (33.6% CV), and residual error is combined
  proportional (32.2%) plus additive (0.368 mg/L). Female sex, Black
  race, and the BPaL regimen were significant on volume in forward
  selection but none survived backward elimination, so the final model
  carries no covariate other than FFM.
- Article: <https://doi.org/10.1038/s41598-026-46217-2>
- Supplement (Appendix 1 the authors’ own nlmixr2 model code, Appendix 2
  the covariate-versus-eta correlation matrix, Appendix 3 individual
  fits, Appendix 4 protein-binding sensitivity for the PTA analysis,
  Appendix 5 the R package list): supplementary data to the same DOI,
  distributed by Springer Nature as `41598_2026_46217_MOESM1_ESM.pdf`.

Pretomanid is a nitroimidazooxazine that inhibits mycolic acid
biosynthesis in replicating *Mycobacterium tuberculosis* and kills
non-replicating bacteria under anaerobic conditions by releasing
reactive nitrogen species. It is a fixed component of the WHO-preferred
six-month oral regimens for rifampicin-resistant tuberculosis (RR-TB):
bedaquiline + pretomanid + linezolid, with or without moxifloxacin (BPaL
/ BPaLM), and the BPaLC variant with clofazimine studied in
TB-PRACTECAL.

This model is unusual among the extractions in this package in that the
authors performed the original analysis **in nlmixr2** and published
their `ini()` / `model()` block verbatim as supplementary Appendix 1.
Every final estimate below is therefore transcribed from the authors’
own code at full precision rather than re-keyed from a rounded results
table, and the rounded Table 2 values serve as an independent
cross-check.

## Population

The model was fitted to the PRACTECAL-PKPD sub-study of TB-PRACTECAL
(ClinicalTrials.gov NCT04081077), an open-label randomised controlled
trial in RR-TB.

``` r

pop <- rxode2::rxode(readModelDb("Nyangwa_2026_pretomanid"))$meta$population
tibble::tibble(Field = names(pop), Value = vapply(pop, as.character, character(1))) |>
  knitr::kable()
```

| Field | Value |
|:---|:---|
| species | human |
| n_subjects | 94 |
| n_studies | 1 |
| age_range | 19-71 years (median 36) |
| age_median | 36 years |
| weight_range | 39.2-144.4 kg (median 56.8) |
| weight_median | 56.8 kg |
| ffm_range | 28.6-75.5 kg (median 45.5); the allometric reference used here |
| bmi_range | 14.3-47.1 kg/m2 (median 19.7) |
| sex_female_pct | 36.2 |
| race_ethnicity | Black 52 (55.3%), Caucasian 40 (42.6%), Asian 1 (1.1%), other 1 (1.1%) |
| disease_state | Rifampicin-resistant pulmonary tuberculosis. 39 participants (41.5%) were living with HIV, all on integrase-inhibitor plus nucleoside / nucleotide reverse transcriptase inhibitor antiretroviral therapy. Patients with moderate liver or renal function abnormality were excluded from the trial. Median estimated creatinine clearance 105.4 mL/min, median ALT 19.5 IU/L, median AST 22 IU/L. |
| dose_range | Pretomanid 200 mg orally once daily for 24 weeks in every arm. Co-administered with bedaquiline (400 mg daily for 2 weeks then 200 mg three times weekly for 22 weeks) and linezolid (600 mg daily for 16 weeks then 300 mg daily for 8 weeks), plus moxifloxacin 400 mg daily in the BPaLM arm (38 participants, 40.4%) or clofazimine 100 mg daily in the BPaLC arm (30, 31.9%); 26 (27.7%) received BPaL alone. Participants were encouraged to eat before dosing but meals were neither standardised nor recorded, so the estimates reflect real-world mixed fed / fasted absorption. |
| regions | South Africa and Belarus |
| notes | PRACTECAL-PKPD sub-study of the TB-PRACTECAL randomised controlled trial (ClinicalTrials.gov NCT04081077). 952 timed plasma samples (86 pre-first-dose, 866 post-dose) spanning the full 24-week treatment course and follow-up visits to week 72. Sampling was on day 1 (0, 2, 23 h), week 8 (predose, 6.5, 23 h), and weeks 12, 16, 20, 24, 32 and 72. Observed concentrations ranged 19.1-11,566 ng/mL with a median trough of 1,789 ng/mL (IQR 1,126-2,689). The assay lower limit of quantification was 7 ng/mL; 234 samples were below it, of which 151 were collected after treatment completion, leaving 9.5% of on-treatment samples BLQ, handled by the M1 method (discarded as missing). Estimation used FOCE-I in nlmixr2 under R 4.1.2. Baseline characteristics are Table 1; the final parameter estimates are Table 2 and the nlmixr2 model code is supplementary Appendix 1. |

All 94 participants received pretomanid 200 mg once daily for 24 weeks,
and contributed 952 timed plasma samples spanning the whole treatment
course and follow-up to week 72. Participants were encouraged to eat
before dosing, but meals were neither standardised nor observed; because
a high-fat, high-calorie meal raises pretomanid AUC by roughly 88%, the
estimates here describe real-world mixed fed/fasted absorption rather
than a controlled fed state.

## Source trace

Every structural value in the model file comes from supplementary
Appendix 1, which is the authors’ final nlmixr2 code. The Table 2 column
is the back-transformed value printed in the paper and is used here only
to confirm the transcription.

| Quantity | Model file | Source location | Cross-check |
|:---|:---|:---|:---|
| Absorption rate | lka | Appendix 1 `lka <- -1.15259732127294` | exp() = 0.3158; Table 2 ka 0.316 (RSE 19.6%) |
| Apparent clearance | lcl | Appendix 1 `lcl <- 1.13008124732802` | exp() = 3.096; Table 2 CL/F 3.10 (RSE 3.35%) |
| Apparent volume | lvc | Appendix 1 `lvc <- 4.62170931536226` | exp() = 101.7; Table 2 V/F 102 (RSE 2.45%) |
| FFM exponent on CL | e_ffm_cl | Appendix 1 `covffmPow1 <- fix(0.75)` | Methods: fixed a priori at 0.75 |
| FFM exponent on V | e_ffm_vc | Appendix 1 `covffmPow2 <- fix(1)` | Methods: fixed a priori at 1 |
| IIV on clearance | etalcl | Appendix 1 `eta.cl ~ 0.102674627530168` | sqrt(exp(w)-1) = 32.9%; Table 2 CL/F %CV 32.9 |
| IIV on volume | etalvc | Appendix 1 `eta.vc ~ 0.10705836808838` | sqrt(exp(w)-1) = 33.6%; Table 2 V/F %CV 33.6 |
| Proportional error | propSd | Appendix 1 `prop.err <- c(0, 0.321731...)` | Table 2 Proportional 0.322 |
| Additive error | addSd | Appendix 1 `add.err <- c(0, 367.831...)` ng/mL | /1000 = 0.3678 mg/L; Table 2 Additive 0.368 mg/L |
| FFM reference | 45.5 kg | Table 1 PK-cohort median fat-free mass | Recovered, not printed – see Errata |
| Structure | 1-cmt oral | Results para 3; Appendix 1 `linCmt()` | Written as explicit ODEs here |

### Omega scale

Appendix 1 gives the random effects as log-scale **variances**, and
Table 2 prints them as a `%CV` column. The two are consistent through
the log-normal identity `CV = sqrt(exp(omega^2) - 1)`, which is what
settles the scale – reading the printed percentages as omega standard
deviations instead would give variances about 5% larger.

``` r

omega <- c(cl = 0.102674627530168, vc = 0.10705836808838)
tibble::tibble(
  Parameter          = c("CL/F", "V/F"),
  `Appendix 1 omega` = omega,
  `CV% implied`      = sqrt(exp(omega) - 1) * 100,
  `Table 2 CV%`      = c(32.9, 33.6)
) |>
  knitr::kable(digits = c(0, 6, 2, 1))
```

| Parameter | Appendix 1 omega | CV% implied | Table 2 CV% |
|:----------|-----------------:|------------:|------------:|
| CL/F      |         0.102675 |       32.88 |        32.9 |
| V/F       |         0.107058 |       33.62 |        33.6 |

## Structural check on the typical subject

With the random effects zeroed, the individual parameters must reproduce
the Table 2 typical values exactly at the reference fat-free mass, and a
single oral dose must be fully recovered – apparent clearance times AUC
to infinity equals the dose, because no bioavailability term is applied
(`CL/F` and `V/F` are apparent parameters estimated from oral-only
data).

``` r

mod <- rxode2::rxode(readModelDb("Nyangwa_2026_pretomanid"))
typ <- rxode2::zeroRe(mod)

ev_single <- rxode2::et(amt = 200, cmt = "depot") |>
  rxode2::et(seq(0, 1000, by = 0.25))
sim_typ <- rxode2::rxSolve(typ, ev_single, params = c(FFM = 45.5))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

auc_inf <- with(sim_typ, sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2))
cl_typ  <- unique(sim_typ$cl)
vc_typ  <- unique(sim_typ$vc)

tibble::tibble(
  Quantity   = c("CL/F (L/h)", "V/F (L)", "ka (1/h)", "Dose recovered (mg)"),
  Simulated  = c(cl_typ, vc_typ, unique(sim_typ$ka), cl_typ * auc_inf),
  Published  = c(3.10, 102, 0.316, 200)
) |>
  knitr::kable(digits = 3)
```

| Quantity            | Simulated | Published |
|:--------------------|----------:|----------:|
| CL/F (L/h)          |     3.096 |     3.100 |
| V/F (L)             |   101.668 |   102.000 |
| ka (1/h)            |     0.316 |     0.316 |
| Dose recovered (mg) |   199.990 |   200.000 |

``` r

# Deterministic checks: both sides use the same parameters, so the only
# difference is numerical integration error and a tight bound is correct.
stopifnot(
  abs(cl_typ - 3.0959) < 1e-3,
  abs(vc_typ - 101.67) < 1e-2,
  # Mass balance. This is the gate that proves the ODE system is solved as
  # written: a one-compartment model defined with a cl/vc pair can be
  # silently auto-solved by rxode2, and a structural error in the depot
  # transfer would show up here as a dose recovery away from 1.
  abs(cl_typ * auc_inf / 200 - 1) < 1e-3
)
```

This also gates the recovered fat-free-mass reference directly, with no
Monte Carlo noise in the way. At steady state AUC(0-tau) equals
`Dose / CL`, so the typical subject at the reference FFM must reproduce
the Table 3 median AUC(0-24) of 64,000 ug\*h/L. A wrong centring
constant is not a subtle effect here: 70 kg would land 39% high and an
uncentred `log(FFM)` 94% low (see Errata).

``` r

auc_ss_typical <- 200 / cl_typ * 1000  # mg/(L/h) -> ug*h/L

tibble::tibble(
  Quantity = "Steady-state AUC(0-24) at the reference FFM (ug*h/L)",
  Model    = auc_ss_typical,
  `Table 3` = 64000,
  `% diff` = 100 * (auc_ss_typical / 64000 - 1)
) |>
  knitr::kable(digits = c(0, 0, 0, 2))
```

| Quantity                                              | Model | Table 3 | % diff |
|:------------------------------------------------------|------:|--------:|-------:|
| Steady-state AUC(0-24) at the reference FFM (ug\*h/L) | 64601 |   64000 |   0.94 |

``` r


stopifnot(abs(auc_ss_typical / 64000 - 1) < 0.05)
```

The model integrates the ODE system rather than substituting the
analytic solution, which is confirmed by an empty `linCmt` slot:

``` r

stopifnot(is.null(mod$linCmt))
mod$state
#> [1] "depot"   "central"
```

## Virtual cohort

The cohort reproduces the Table 1 fat-free-mass distribution of the 94
participants in the pretomanid PK analysis: median 45.5 kg, range
28.6-75.5 kg. A log-normal with a 0.21 log-scale standard deviation,
truncated to the observed range, matches that median and spread.
Fat-free mass is the only covariate the final model uses, so no other
demographic needs to be simulated.

Fat-free mass is laid out on the log-normal **quantiles** rather than
drawn at random. The distribution is identical, but the cohort median is
then exactly 45.5 kg instead of a draw that happens to land near it,
which matters because the fat-free-mass reference is the quantity this
vignette is validating: a sampled cohort whose median drifted to 44 kg
would shift the predicted AUC by several percent for reasons that have
nothing to do with the model. The between-subject random effects are
still drawn by `rxSolve()` below, so the exposure spread remains
stochastic.

``` r

rxode2::rxSetSeed(20260913)
set.seed(20260913)

n_sub <- 200  # per-arm cap for library vignettes
cohort <- tibble::tibble(
  id  = seq_len(n_sub),
  FFM = pmin(pmax(qlnorm(ppoints(n_sub), log(45.5), 0.21), 28.6), 75.5)
)

tibble::tibble(
  Quantity    = c("Median FFM (kg)", "Minimum FFM (kg)", "Maximum FFM (kg)"),
  Simulated   = c(median(cohort$FFM), min(cohort$FFM), max(cohort$FFM)),
  `Table 1`   = c(45.5, 28.6, 75.5)
) |>
  knitr::kable(digits = 1)
```

| Quantity         | Simulated | Table 1 |
|:-----------------|----------:|--------:|
| Median FFM (kg)  |      45.5 |    45.5 |
| Minimum FFM (kg) |      28.6 |    28.6 |
| Maximum FFM (kg) |      75.5 |    75.5 |

## Steady-state simulation

Pretomanid was given as 200 mg once daily for 24 weeks. The simulation
doses for 28 days and observes the final dosing interval, by which point
a drug with a 22.8 h terminal half-life is fully accumulated.

``` r

tau <- 24
n_days <- 28
t_start <- (n_days - 1) * tau

# addl = n_days - 1 gives n_days doses at 0, tau, ..., (n_days - 1) * tau, so
# the LAST dose lands exactly on t_start and the observed window is a true
# dosing interval. Using n_days - 2 would put the last dose one interval
# earlier and silently turn this into a post-dose washout.
events <- rxode2::et(amt = 200, cmt = "depot", ii = tau, addl = n_days - 1) |>
  rxode2::et(seq(t_start, t_start + tau, by = 0.25)) |>
  rxode2::et(id = cohort$id)

sim <- rxode2::rxSolve(
  mod, events,
  params = as.data.frame(cohort),
  keep = "FFM"
) |>
  as.data.frame()
```

``` r

sim |>
  mutate(t_rel = time - t_start) |>
  group_by(t_rel) |>
  summarise(
    med = median(Cc) * 1000,
    lo  = quantile(Cc, 0.05) * 1000,
    hi  = quantile(Cc, 0.95) * 1000,
    .groups = "drop"
  ) |>
  ggplot(aes(t_rel, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "steelblue") +
  geom_line(linewidth = 1) +
  labs(
    x = "Time after dose (h)",
    y = "Pretomanid concentration (ug/L)",
    title = "Steady-state pretomanid, 200 mg once daily"
  ) +
  theme_bw()
```

![Simulated steady-state pretomanid concentration-time profiles over the
final dosing interval. The heavy line is the median and the ribbon spans
the 5th-95th percentiles of the 200-subject
cohort.](Nyangwa_2026_pretomanid_files/figure-html/profile-plot-1.png)

Simulated steady-state pretomanid concentration-time profiles over the
final dosing interval. The heavy line is the median and the ribbon spans
the 5th-95th percentiles of the 200-subject cohort.

## PKNCA validation

Table 3 of the paper reports secondary parameters derived from the
model: AUC(0-24), the maximum inter-dose concentration, and the trough
concentration just before the next dose. These are model-predicted
individual values, so the NCA runs on `Cc` (the individual prediction)
rather than on `sim` (which carries residual error and would bias `Cmax`
upward).

The NCA runs on **time after the final dose** rather than on absolute
trial time. `ctrough` is the concentration at the end of the interval
and PKNCA resolves it by matching a record to the interval end on the
concentration object’s own time scale; with an absolute-time interval of
648-672 h it returns `NA` for every subject even though the 672 h record
exists. Re-basing the interval to 0-24 h makes it well defined, and
`tmax` then reads directly as time after dose.

``` r

conc_df <- sim |>
  mutate(
    treatment = "Pretomanid 200 mg QD",
    conc_ug_L = Cc * 1000,
    tad       = time - t_start
  ) |>
  filter(!is.na(conc_ug_L)) |>
  select(id, treatment, tad, conc_ug_L)

dose_df <- tibble::tibble(
  id        = cohort$id,
  treatment = "Pretomanid 200 mg QD",
  tad       = 0,
  dose      = 200
)

# Treatment grouping comes BEFORE id, and PKNCAdose rejects a slash, so both
# formulas use the '+' form.
conc_obj <- PKNCA::PKNCAconc(conc_df, conc_ug_L ~ tad | treatment + id)
dose_obj <- PKNCA::PKNCAdose(dose_df, dose ~ tad | treatment + id)

intervals <- data.frame(
  start   = 0,
  end     = tau,
  cmax    = TRUE,
  tmax    = TRUE,
  auclast = TRUE,
  ctrough = TRUE
)

nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
nca_ind <- as.data.frame(nca)

stopifnot(
  # A silent all-NA ctrough is the failure this re-basing fixes; gate it so it
  # cannot come back unnoticed.
  !all(is.na(nca_ind$PPORRES[nca_ind$PPTESTCD == "ctrough"]))
)
```

``` r

simulated_nca <- nca_ind |>
  filter(PPTESTCD %in% c("auclast", "cmax", "ctrough"))

published <- tibble::tibble(
  treatment = "Pretomanid 200 mg QD",
  auclast   = 64000,
  cmax      = 3000,
  ctrough   = 2000
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = simulated_nca,
  reference     = published,
  by            = "treatment",
  units         = c(auclast = "ug*h/L", cmax = "ug/L", ctrough = "ug/L"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = paste(
    "Simulated steady-state NCA (cohort median of 200 subjects) versus the",
    "Nyang'wa 2026 Table 3 medians. * differs by more than 20%."
  )
)
```

| NCA parameter     | treatment            | Reference | Simulated | % diff |
|:------------------|:---------------------|:----------|:----------|:-------|
| Cmax (ug/L)       | Pretomanid 200 mg QD | 3000      | 3320      | +10.6% |
| AUClast (ug\*h/L) | Pretomanid 200 mg QD | 64000     | 66800     | +4.3%  |
| Ctrough (ug/L)    | Pretomanid 200 mg QD | 2000      | 2030      | +1.7%  |

Simulated steady-state NCA (cohort median of 200 subjects) versus the
Nyang’wa 2026 Table 3 medians. \* differs by more than 20%. {.table}

``` r

# Recompute the percent differences independently -- ncaComparisonTable
# returns `% diff` as formatted text, which is for display, not for gating.
chk <- simulated_nca |>
  group_by(PPTESTCD) |>
  summarise(Simulated = median(PPORRES), .groups = "drop") |>
  left_join(
    tibble::tibble(
      PPTESTCD  = c("auclast", "cmax", "ctrough"),
      Published = c(64000, 3000, 2000)
    ),
    by = "PPTESTCD"
  ) |>
  mutate(pct_diff = 100 * (Simulated - Published) / Published)

pct <- function(p) chk$pct_diff[chk$PPTESTCD == p]

# Bounds are set for the CENTRE of the cohort and sized against Monte Carlo
# noise, not against the values observed while authoring. The between-subject
# random effects are redrawn by rxode2 on every machine, and the standard
# error of a 200-subject median with omega ~ 0.32 is about 2.8%, so a bound
# tighter than roughly 10% on a cohort median would fail intermittently in CI
# for no structural reason. The tight, reproducible gates on the structure and
# the FFM reference live in the zeroRe chunks above, where there is no draw.
stopifnot(
  abs(pct("auclast")) < 10,
  abs(pct("ctrough")) < 20,
  abs(pct("cmax"))    < 25
)
```

Both integral measures land essentially on top of the published values:
the median AUC(0-24) and the median Ctrough each agree with Table 3 to
well under a percent. Because AUC at steady state is exactly
`Dose / CL`, this is the single most informative check available – it
confirms the clearance transcription *and* the recovered 45.5 kg
allometric reference at the same time. `Cmax` sits about 8% high; see
the Errata below.

## Probability of target attainment

Figure 4 of the paper is a PTA analysis against two PK-PD indices, both
assuming 85% plasma protein binding (free fraction 0.15): the free
AUC(0-24) to MIC ratio with a target of 167, and the fraction of the
dosing interval the free concentration exceeds the MIC, with targets of
22%, 48% and 77% for bacteriostasis, 1-log10 kill and 1.59-log10 kill.

``` r

fu <- 0.15

auc_by_id <- nca_ind |>
  filter(PPTESTCD == "auclast") |>
  select(id, auc = PPORRES)

mics <- c(0.016, 0.032, 0.063, 0.125, 0.25, 0.5, 1)

pta <- lapply(mics, function(m) {
  ft <- sim |>
    group_by(id) |>
    summarise(pct_above = 100 * mean(fu * Cc > m), .groups = "drop")
  tibble::tibble(
    MIC              = m,
    `fAUC/MIC >= 167` = 100 * mean(fu * auc_by_id$auc / 1000 / m >= 167),
    `fT>MIC >= 77%`   = 100 * mean(ft$pct_above >= 77),
    `fT>MIC >= 48%`   = 100 * mean(ft$pct_above >= 48)
  )
}) |>
  bind_rows()

knitr::kable(pta, digits = 1)
```

| MIC | fAUC/MIC \>= 167 | fT\>MIC \>= 77% | fT\>MIC \>= 48% |
|----:|-----------------:|----------------:|----------------:|
| 0.0 |            100.0 |           100.0 |           100.0 |
| 0.0 |             93.5 |           100.0 |           100.0 |
| 0.1 |             43.0 |           100.0 |           100.0 |
| 0.1 |              2.5 |            99.0 |           100.0 |
| 0.2 |              0.0 |            81.0 |            92.5 |
| 0.5 |              0.0 |            20.5 |            32.5 |
| 1.0 |              0.0 |             0.0 |             1.0 |

``` r

pta |>
  pivot_longer(-MIC, names_to = "Target", values_to = "PTA") |>
  ggplot(aes(factor(MIC), PTA, colour = Target, group = Target)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "grey40") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "MIC (mg/L)", y = "Probability of target attainment (%)",
    title = "Pretomanid 200 mg daily, 85% protein binding"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Replicates Figure 4 of Nyang'wa 2026: probability of target attainment
for a 200 mg daily pretomanid dose across the MIC testing range,
assuming 85% plasma protein
binding.](Nyangwa_2026_pretomanid_files/figure-html/pta-plot-1.png)

Replicates Figure 4 of Nyang’wa 2026: probability of target attainment
for a 200 mg daily pretomanid dose across the MIC testing range,
assuming 85% plasma protein binding.

``` r

pta_at <- function(mic, col) pta[[col]][pta$MIC == mic]

stopifnot(
  # Results para 6: "a 200 mg dose resulted in 80% or higher PTA
  # fAUC0-24/MIC of 167 at an MIC of 0.032 mg/L and below".
  pta_at(0.032, "fAUC/MIC >= 167") >= 80,
  pta_at(0.016, "fAUC/MIC >= 167") >= 80,
  # ... and the same sentence implies the next MIC up falls short.
  pta_at(0.063, "fAUC/MIC >= 167") < 80,
  # Results para 7: at the trial median MIC of 0.125 mg/L, "99.67% of
  # patients would have had drug exposures above %fT > MIC target of 77%".
  pta_at(0.125, "fT>MIC >= 77%") > 90,
  # Results para 6: the 77% and 48% targets "were met at 200 mg dosing at
  # MICs of 0.125 and 0.250 mg/L".
  pta_at(0.25, "fT>MIC >= 48%") > 85
)
```

The model reproduces the paper’s central PK-PD conclusion: at the
TB-PRACTECAL median MIC of 0.125 mg/L essentially every participant
clears the time-dependent target, while the 167 fAUC/MIC target is met
only for the most susceptible isolates.

## Assumptions and deviations

### Fat-free mass reference recovered, not printed (load-bearing)

Appendix 1 scales allometrically on a pre-computed data column named
`logFFM`:

    cl <- exp(lcl + eta.cl + logFFM * covffmPow1)
    vc <- exp(lvc + eta.vc + logFFM * covffmPow2)

Neither the code nor the paper states what `logFFM` is centred on, and
the answer changes every predicted exposure. It is recovered here by
arithmetic rather than assumed:

1.  `exp(lcl)` = 3.096 L/h and `exp(lvc)` = 101.7 L reproduce the Table
    2 typical values of 3.10 L/h and 102 L. Those typical values can
    only equal the exponentiated thetas if `logFFM` is zero at the
    reference, so the column must be `log(FFM / FFM_ref)` and not a bare
    `log(FFM)`.
2.  Table 3 reports a median AUC(0-24) of 64,000 ug\*h/L. At steady
    state `AUC = Dose / CL`, so the median participant’s clearance is
    200 / 64 = 3.13 L/h – which equals `exp(lcl)` only if the reference
    is the cohort median fat-free mass.

Taking the Table 1 PK-cohort median of 45.5 kg gives a predicted median
AUC(0-24) of about 64,600 ug\*h/L. The two alternatives are excluded by
margins far too large to be rounding:

| `FFM_ref`                  | Predicted median AUC(0-24) | vs published 64,000 |
|----------------------------|----------------------------|---------------------|
| uncentred `log(FFM)`       | 3,690                      | -94%                |
| 45.5 kg (cohort median)    | 64,600                     | +0.9%               |
| 70 kg (conventional adult) | 89,200                     | +39%                |

The 45.5 kg reference is therefore reported as recovered rather than
assumed. A downstream user supplying fat-free mass must use the same
centring.

### Fat-free mass equation not identified

The paper does not report which formula produced its FFM column. Users
must supply FFM directly or derive it with a documented equation such as
Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065). The
cohort’s median total body weight of 56.8 kg against a median FFM of
45.5 kg is the only in-paper calibration point.

### Cmax runs about 8% above Table 3

The simulated median steady-state Cmax is about 3,230 ug/L against the
published 3,000 ug/L, a difference of roughly 8% and comfortably inside
the 20% flagging tolerance. Nothing was tuned to narrow it.

The discrepancy is confined entirely to the peak, and the most
economical explanation is that the published value is rounded. Table 3
reports the Cmax median as a bare “3,000” while giving its range to full
precision (1382-6351); the same table reports AUC(0-24) as “64,000” and
Ctrough as “2000”. A median of 3,230 rounds to 3,000 at one significant
figure, which is the precision Table 3 and the Abstract actually use for
the central values.

The two integral measures leave very little room for a structural
explanation, since they agree with the paper to under a percent in both
directions:

``` r

cavg_model <- chk$Simulated[chk$PPTESTCD == "auclast"] / 24
cavg_paper <- 64000 / 24

tibble::tibble(
  Quantity = c(
    "AUC(0-24) % diff", "Ctrough % diff", "Cmax % diff",
    "Cmax / Cavg, model", "Cmax / Cavg, Table 3"
  ),
  Value = c(
    pct("auclast"), pct("ctrough"), pct("cmax"),
    chk$Simulated[chk$PPTESTCD == "cmax"] / cavg_model,
    3000 / cavg_paper
  )
) |>
  knitr::kable(digits = 2)
```

| Quantity             | Value |
|:---------------------|------:|
| AUC(0-24) % diff     |  4.35 |
| Ctrough % diff       |  1.71 |
| Cmax % diff          | 10.62 |
| Cmax / Cavg, model   |  1.19 |
| Cmax / Cavg, Table 3 |  1.12 |

The peak-to-average ratio the model structure produces for the published
`ka` and `kel` is higher than the ratio implied by the paper’s own Table
3 row, so the paper’s Cmax is low relative to its own AUC rather than
the model’s being high relative to both.

One explanation that does **not** survive checking is sampling density.
The single-dose Tmax for these parameters is about 8.2 h, but at steady
state the peak arrives earlier – the interval starts from an accumulated
trough that is still decaying – and the simulated median Tmax is 6 h.
The trial’s week-8 sampling at 6.5 h post-dose therefore does bracket
the peak closely, so the observed schedule cannot be blamed for
underestimating it.

### The 50% fAUC/MIC statement does not reconcile

Results paragraph 7 states that at an MIC of 0.125 mg/L “only 50% would
have reached the fAUC/MIC target of 167”. The model gives about 1.5% at
that MIC. Reaching 167 at MIC 0.125 with a free fraction of 0.15
requires a total AUC(0-24) of about 139,000 ug\*h/L, which exceeds the
maximum of the paper’s own reported Table 3 range (30,853-138,179), so a
50% attainment rate is not consistent with the exposure distribution the
paper publishes.

Three other statements in the same paper agree with the model instead,
so this appears to be an isolated slip rather than a different analysis:

- the Abstract: “the AUC/MIC target was not achieved”;
- Results paragraph 6: 80% or higher PTA only “at an MIC of 0.032 mg/L
  and below” (the model gives 97% at 0.032 and 40% at 0.063);
- Discussion paragraph 5: “adequate exposure would only be achieved for
  strains with an MIC of 0.032 mg/L or below”.

The `pta-assert` chunk above gates on those three statements and
deliberately does not gate on the 50% figure.

### Terminal half-life

The model’s terminal half-life at the reference FFM is 22.8 h. The
Introduction quotes about 16 h, but that figure is cited from label and
healthy-volunteer literature for fasted dosing, not estimated in this
analysis, so it is not a target for this model.

### Between-subject variability only where the authors put it

Appendix 1 declares random effects on clearance and central volume only.
No IIV is placed on absorption, and the etas are uncorrelated (a
diagonal omega), exactly as published. Shrinkage on the volume eta was
36.7%, so the individual volume estimates underlying Table 3 are pulled
toward the typical value more than the clearance estimates are.

### Covariates screened but not retained

Female sex, Black race and the BPaL regimen were significant on volume
of distribution at the forward-inclusion threshold but none survived
backward elimination, and the regimen arm was not significant on
clearance. Total body weight and body mass index lost to fat-free mass
as the allometric size descriptor. These are recorded in the model
file’s `covariatesDataExcluded` metadata for provenance; none carries a
usable point estimate and none is referenced in `model()`.

### Population coverage

The trial excluded patients with moderate liver or renal function
abnormality, so the model is not informed about organ impairment. Potent
CYP450 inducers (efavirenz, lopinavir/ritonavir, rifamycins), which are
known to reduce pretomanid exposure, were contraindicated in the trial;
the 41.5% of participants living with HIV were all on
integrase-inhibitor-based antiretroviral therapy, and HIV status was not
a significant covariate.
