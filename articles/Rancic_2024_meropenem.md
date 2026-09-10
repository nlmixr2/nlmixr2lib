# Meropenem (Rancic 2024)

## Model and source

- Citation: Rancic A, Milosavljevic MN, Rosic N, Milovanovic D, Folic M,
  Ruzic Zecevic D, Petrovic N, Milojevic Corbic M, Dabanovic V, Jankovic
  SM. Population pharmacokinetics of meropenem in critically ill
  patients. Open Med (Wars). 2024;19(1):20241004.
  <doi:10.1515/med-2024-1004>. PMCID PMC11278387.
- Description: One-compartment intravenous population PK model for
  meropenem in critically ill adults in intensive care (Rancic 2024).
  Clearance combines a multiplicative power term in serum creatinine and
  white blood cell count with additive shifts for hypertension and for
  concomitant vancomycin or colistimethate; the central volume carries
  no covariates. NOTE: the published central volume (2.05 L) is roughly
  ten-fold smaller than the value implied by the paper’s own reported
  peak concentrations, so simulated peaks are correspondingly high - see
  the validation vignette.
- Article: <https://doi.org/10.1515/med-2024-1004> (Open Medicine, open
  access; PMCID PMC11278387)

No supplementary material accompanies the article; the Data Availability
statement offers the raw data “upon motivated request to the
corresponding author” only. Every value below therefore comes from the
main text and its five tables.

## Population

The model was fitted to 101 critically ill adults treated in the
intensive care unit of the University Clinical Centre Kragujevac,
Serbia, all with a severe infection (meningitis, pneumonia, sepsis,
septic shock or febrile neutropenia) caused by a multi-resistant
Gram-negative organism. Mean age was 62.4 +/- 14.9 years (range 21 -
86), mean total body weight 79.0 +/- 13.8 kg (range 48 - 130), and 39 of
101 (38.6%) were women. Mean serum creatinine was 94.5 +/- 63.7 umol/L
(range 40 - 452), so most of the cohort was not renally impaired – the
authors emphasise *augmented* renal clearance rather than renal failure
as the dominant phenomenon in this population (Discussion). Recorded
comorbidities were hypertension 32.7%, chronic renal failure 11.9%,
neoplasm 17.8%, cerebral infarction 20.8%, pneumonia 31.7% and urinary
tract infection 35.6% (Table 1).

Meropenem was given as an intermittent intravenous infusion of 1,000 -
2,000 mg every 8 or 12 h, for a mean total daily dose of 3,000 +/- 693
mg/day (range 2,000 - 6,000). Twenty-seven patients (26.7%) received
antibiotic polytherapy: 17 (16.8%) vancomycin, 7 (6.9%) colistin, and 3
(2.97%) both. Patients entered the study only after at least three days
of continuous meropenem, i.e. at steady state, and contributed exactly
two plasma samples each (202 observations total): the first 5 - 30 min
after the end of the infusion (mean 40.69 +/- 16.67 mg/L, range 13.07 -
88.95) and the second 3 - 4 h after the end of the infusion (mean 12.55
+/- 7.62 mg/L, range 2.06 - 36.28).

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Rancic_2024_meropenem")()$population`).

``` r

pop <- rxode2::rxode(readModelDb("Rancic_2024_meropenem"))$population
#> ℹ parameter labels from comments will be replaced by 'label()'
str(pop, max.level = 1)
#> List of 11
#>  $ species       : chr "human"
#>  $ n_subjects    : int 101
#>  $ n_studies     : int 1
#>  $ n_observations: int 202
#>  $ age_range     : chr "21 - 86 years (mean 62.37 +/- 14.89)"
#>  $ weight_range  : chr "48 - 130 kg (mean 78.97 +/- 13.76)"
#>  $ sex_female_pct: num 38.6
#>  $ disease_state : chr "Critically ill adults in the intensive care unit with severe infection (meningitis, pneumonia, sepsis, septic s"| __truncated__
#>  $ dose_range    : chr "1,000 - 2,000 mg meropenem every 8 or 12 h by intermittent intravenous infusion (total daily dose mean 3,000 +/"| __truncated__
#>  $ regions       : chr "Serbia (Intensive Care Unit, University Clinical Centre Kragujevac)."
#>  $ notes         : chr "Prospective observational case-series. Two plasma samples per patient: the first 5 - 30 min after the end of th"| __truncated__
```

## Source trace

Every `ini()` entry in
`inst/modeldb/specificDrugs/Rancic_2024_meropenem.R` carries an in-file
comment pointing at its origin. They are collected here for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| Structural model: one compartment, IV, no absorption | n/a | Methods 2.3: NONMEM “ADVAN 1 subroutine (version 5, level 1.1, double precision) … a single-compartment model without absorption” |
| Final clearance equation | `CL = 5.29 * CRE^1e-6 * WBC^-0.165 + 1e-6*HTA + 0.825*VAN + 1.28*COL` | Abstract; Results paragraph following Table 2 (the display equation itself is a figure in the PDF and is transcribed from the Abstract, which prints it in full) |
| `lcl` = `log(5.29)` | 5.29 L/h | Table 3, row “Clearance (L/h) (theta 1)”; SE 0.038, 95% CI 5.215 - 5.365 |
| `lvc` = `log(2.05)` | 2.05 L | Table 3, row “Volume of distribution (L) (theta 2)”; SE 0.006, 95% CI 2.039 - 2.061 |
| `e_creat_cl` | 0.000001 | Abstract equation; Table 3 row “Effect of CRE” prints estimate, SE and both CI limits as 0.0000 |
| `e_wbc_cl` | -0.165 | Table 3, row “Effect of WBCs”; SE 0.0545, 95% CI -0.219 to -0.110 |
| `e_dis_hypert_cl` | 0.000001 | Abstract equation; Table 3 row “Effect of HTA” prints estimate, SE and both CI limits as 0.0000 |
| `e_conmed_vancomycin_cl` | 0.825 L/h | Table 3, row “Effect of VAN”; SE 0.054, 95% CI 0.770 - 0.879 |
| `e_conmed_colistimethate_cl` | 1.28 L/h | Table 3, row “Effect of COL”; SE 0.05, 95% CI 1.23 - 1.33 |
| `etalcl` variance | 0.0215 | Table 3, row “Inter-individual variance of CL (omega^2 CL)”; SE 0.0146, 95% CI 0.0000 - 0.0501 |
| Exponential IIV placement (multiplies the whole bracket) | n/a | Table 2, base / univariate / full model rows, all of the form `[...] x Exp[ETA(1)]` |
| `propSd` = `sqrt(0.44)` = 0.6633 | 0.44 (variance) | Table 3, row “Residual error variance (sigma^2 CL)”; SE 0.066, 95% CI 0.310 - 0.570. Form from Methods 2.3; see Errata |
| Volume: no covariates, no IIV | n/a | Table 3 reports no volume covariate rows and no `omega^2 V`; the Table 3 footnote defines ETA(2) and theta16 - theta20 but no estimates for them appear anywhere |
| Cohort demographics | n/a | Table 1 |
| Covariate screen (24 screened, 13 in the full model, 5 retained) | n/a | Results paragraphs 3 - 4; Table 2 |

## The covariate screen

Table 2 of the paper records the objective-function value for the base
model, for each univariate covariate model, and for the full model. Only
the five covariates that survived backward deletion at p \< 0.01 appear
in the final model; the rest are recorded in the model file’s
`covariatesDataExcluded` metadata (where a canonical register name
exists) or in the narrative here.

``` r

screen <- tibble::tribble(
  ~Model,                                       ~MOF,      ~dMOF,  ~p,
  "Base: CL = theta1 * Exp[ETA(1)]",             1283.053,  NA,     NA_character_,
  "+ (AGE/50)^theta3",                           1277.958,  5.095,  "<0.05",
  "+ (TBW/70)^theta4",                           1266.640, 16.413,  "<0.01",
  "+ (daily dose)^theta5",                       1282.624,  0.429,  ">0.05",
  "+ (CRE)^theta6",                              1272.552, 10.501,  "<0.01",
  "+ (WBCs)^theta7",                             1267.150, 15.903,  "<0.01",
  "+ theta8 * HTA",                              1262.636, 20.417,  "<0.01",
  "+ theta9 * CRF",                              1282.954,  0.099,  ">0.05",
  "+ theta10 * NEO",                             1280.925,  2.128,  ">0.05",
  "+ theta11 * CI",                              1264.258, 18.795,  "<0.01",
  "+ theta12 * PNE",                             1275.855,  7.198,  "<0.01",
  "+ theta13 * UTI",                             1272.183, 10.870,  "<0.01",
  "+ theta14 * VAN",                             1272.001, 11.052,  "<0.01",
  "+ theta15 * COL",                             1274.813,  8.240,  "<0.01",
  "Full model (all significant covariates)",     1207.153,  NA,     NA_character_,
  "Final model (after backward deletion)",       1224.035, NA,      NA_character_
)

# The MOF differences the paper prints are all relative to the base model.
chk <- screen |>
  dplyr::filter(!is.na(dMOF)) |>
  dplyr::mutate(recomputed = 1283.053 - MOF)
stopifnot(max(abs(chk$recomputed - chk$dMOF)) < 5e-4)

screen |>
  dplyr::rename(
    "Clearance model"    = Model,
    "MOF"                = MOF,
    "Difference vs base" = dMOF,
    "p value"            = p
  ) |>
  knitr::kable(
    digits  = 3,
    caption = paste(
      "Reproduces Table 2 of Rancic 2024. Every printed MOF difference is",
      "verified against the base-model MOF of 1283.053 in the chunk above."
    )
  )
```

| Clearance model                         |      MOF | Difference vs base | p value |
|:----------------------------------------|---------:|-------------------:|:--------|
| Base: CL = theta1 \* Exp\[ETA(1)\]      | 1283.053 |                 NA | NA      |
| \+ (AGE/50)^theta3                      | 1277.958 |              5.095 | \<0.05  |
| \+ (TBW/70)^theta4                      | 1266.640 |             16.413 | \<0.01  |
| \+ (daily dose)^theta5                  | 1282.624 |              0.429 | \>0.05  |
| \+ (CRE)^theta6                         | 1272.552 |             10.501 | \<0.01  |
| \+ (WBCs)^theta7                        | 1267.150 |             15.903 | \<0.01  |
| \+ theta8 \* HTA                        | 1262.636 |             20.417 | \<0.01  |
| \+ theta9 \* CRF                        | 1282.954 |              0.099 | \>0.05  |
| \+ theta10 \* NEO                       | 1280.925 |              2.128 | \>0.05  |
| \+ theta11 \* CI                        | 1264.258 |             18.795 | \<0.01  |
| \+ theta12 \* PNE                       | 1275.855 |              7.198 | \<0.01  |
| \+ theta13 \* UTI                       | 1272.183 |             10.870 | \<0.01  |
| \+ theta14 \* VAN                       | 1272.001 |             11.052 | \<0.01  |
| \+ theta15 \* COL                       | 1274.813 |              8.240 | \<0.01  |
| Full model (all significant covariates) | 1207.153 |                 NA | NA      |
| Final model (after backward deletion)   | 1224.035 |                 NA | NA      |

Reproduces Table 2 of Rancic 2024. Every printed MOF difference is
verified against the base-model MOF of 1283.053 in the chunk above.
{.table}

Four screened comorbidity flags – chronic renal failure, cerebral
infarction, pneumonia and urinary tract infection – are deliberately
*not* given canonical covariate names in the model file. Three of them
were significant univariately and were carried into the full model, but
none survived backward deletion, so none carries an effect; and their
mapping onto the existing register entries (`DIS_CUTI` is specifically
*complicated* urinary tract infection, `DIS_HABP` specifically
hospital-acquired bacterial pneumonia) is not settled by what this paper
reports. They are recorded here rather than minted as new canonical
names on the strength of a screen that rejected them.

## Virtual cohort

Original observed data are not publicly available. The cohort below is a
101-subject virtual population – the same size as the paper’s – whose
covariate distributions reproduce the Table 1 marginals.

Two distributions had to be assumed, because the paper does not report
them:

- **Leukocyte count.** Table 1 has no WBC row at all, even though WBC is
  one of the five retained covariates. The *unit* is settled by internal
  consistency (see the model file’s `covariateData[[WBC]]$notes`): only
  a `10^9/L` scale reconciles the final model’s covariate-adjusted
  clearance with the base model’s 3.80 L/h. The *distribution* is
  assumed log-normal with median 10.8 x 10^9/L – the value implied by
  that same reconciliation – and a 45% geometric CV, truncated to 2 -
  40, which spans the range seen in an ICU cohort with severe
  Gram-negative infection.
- **Dosing regimen mix.** The paper reports the daily-dose *summary*
  (3,000 +/- 693 mg/day, range 2,000 - 6,000) and the permitted regimens
  (1,000 - 2,000 mg q8h or q12h) but not the mix. The mix below was
  chosen to reproduce the reported mean and standard deviation, and that
  reproduction is asserted in the chunk.

``` r

# set.seed() seeds R's RNG, which is what draws the covariates below. It does
# NOT seed rxode2's simulation RNG (that draws etalcl, and its streams are
# partitioned per solver thread), so the cohort's covariates are reproducible
# but its etas are not. Every assertion downstream is either a per-subject
# algebraic identity or a typical-value replication, so none of them depends on
# which etas were drawn. See pattern 12 of
# references/known-vignette-failure-patterns.md.
set.seed(20240701)
rxode2::rxSetSeed(20240701)

n_sub <- 101L

# Regimen mix: 15 / 76 / 7 / 3 reproduces the Table 1 daily-dose mean and SD.
regimens <- tibble::tribble(
  ~regimen,        ~amt,  ~tau, ~n,
  "1000 mg q12h",  1000,    12, 15L,
  "1000 mg q8h",   1000,     8, 76L,
  "2000 mg q12h",  2000,    12,  7L,
  "2000 mg q8h",   2000,     8,  3L
)
stopifnot(sum(regimens$n) == n_sub)

subj <- regimens |>
  dplyr::mutate(daily = amt * 24 / tau) |>
  tidyr::uncount(n) |>
  dplyr::mutate(id = seq_len(dplyr::n()))

# Table 1: total daily dose 3,000 +/- 692.82 mg/day, range 2,000 - 6,000.
stopifnot(
  abs(mean(subj$daily) - 3000) < 50,
  abs(stats::sd(subj$daily) - 692.82) < 50,
  min(subj$daily) == 2000,
  max(subj$daily) == 6000
)

# Serum creatinine: Table 1 mean 94.47 +/- 63.68 umol/L, range 40 - 452. The
# marginal is strongly right-skewed, so a log-normal matched on the mean and SD
# and truncated to the reported range is used.
creat_cv    <- 63.68 / 94.47
creat_sdlog <- sqrt(log(1 + creat_cv^2))
creat_mulog <- log(94.47) - creat_sdlog^2 / 2
subj$CREAT <- pmin(452, pmax(40, stats::rlnorm(n_sub, creat_mulog, creat_sdlog)))

# Leukocyte count: NOT reported by the paper (see prose above). Assumed.
subj$WBC <- pmin(40, pmax(2, stats::rlnorm(n_sub, log(10.8), sqrt(log(1 + 0.45^2)))))

# Comorbidity and comedication indicators, at the Table 1 prevalences. The
# vancomycin and colistimethate flags are NOT mutually exclusive: 17 patients
# received vancomycin, 7 colistin, and 3 both, so 14 / 4 / 3 / 80.
subj$DIS_HYPERT <- 0
subj$DIS_HYPERT[sample.int(n_sub, 33L)] <- 1        # 33/101 = 32.7%

combo <- sample.int(n_sub, 21L)
subj$CONMED_VANCOMYCIN     <- 0
subj$CONMED_COLISTIMETHATE <- 0
subj$CONMED_VANCOMYCIN[combo[1:14]]      <- 1       # vancomycin only
subj$CONMED_COLISTIMETHATE[combo[15:18]] <- 1       # colistimethate only
subj$CONMED_VANCOMYCIN[combo[19:21]]     <- 1       # both
subj$CONMED_COLISTIMETHATE[combo[19:21]] <- 1

stopifnot(
  sum(subj$DIS_HYPERT) == 33L,
  sum(subj$CONMED_VANCOMYCIN) == 17L,
  sum(subj$CONMED_COLISTIMETHATE) == 7L,
  sum(subj$CONMED_VANCOMYCIN == 1 & subj$CONMED_COLISTIMETHATE == 1) == 3L
)

knitr::kable(
  subj |>
    dplyr::group_by(regimen) |>
    dplyr::summarise(
      n              = dplyr::n(),
      `Daily dose`   = unique(daily),
      `Median CREAT` = round(stats::median(CREAT), 1),
      `Median WBC`   = round(stats::median(WBC), 1),
      .groups = "drop"
    ),
  caption = "Virtual cohort by regimen (n = 101, matching the published cohort size)."
)
```

| regimen      |   n | Daily dose | Median CREAT | Median WBC |
|:-------------|----:|-----------:|-------------:|-----------:|
| 1000 mg q12h |  15 |       2000 |        102.7 |        9.7 |
| 1000 mg q8h  |  76 |       3000 |         68.1 |       10.5 |
| 2000 mg q12h |   7 |       4000 |        101.1 |        9.6 |
| 2000 mg q8h  |   3 |       6000 |         99.2 |       13.2 |

Virtual cohort by regimen (n = 101, matching the published cohort size).
{.table}

The paper’s inclusion criterion was at least three days of continuous
therapy, so every subject is dosed for 72 h and the final dosing
interval is the one analysed.

``` r

# Infusion duration is NOT stated: Methods 2.2 says only "intermittent
# intravenous infusion", with the first sample drawn 5 - 30 min after its end.
# A 30-minute infusion is assumed (see Assumptions and deviations).
inf_dur <- 0.5

# Dense early sampling is needed because the published volume implies a
# ~0.4 h half-life; a coarse grid would understate AUC and mis-fit half-life.
obs_grid <- function(tau) {
  sort(unique(c(seq(0, 1, by = 0.02), seq(1, 4, by = 0.05), seq(4, tau, by = 0.25))))
}

make_subject <- function(row) {
  dose_times <- seq(0, 72 - row$tau, by = row$tau)
  last_dose  <- max(dose_times)
  doses <- tibble::tibble(
    time = dose_times, amt = row$amt, evid = 1L, rate = row$amt / inf_dur,
    cmt = "central"
  )
  # cmt on an observation row must be an ODE STATE name, never the algebraic
  # observable "Cc" (that auto-injects a compartment slot and renumbers the
  # model). rxode2 returns Cc as a column regardless.
  obs <- tibble::tibble(
    time = last_dose + obs_grid(row$tau), amt = NA_real_, evid = 0L,
    rate = NA_real_, cmt = "central"
  )
  dplyr::bind_rows(doses, obs) |>
    dplyr::mutate(
      id = row$id, regimen = row$regimen, last_dose = last_dose, tau = row$tau,
      CREAT = row$CREAT, WBC = row$WBC, DIS_HYPERT = row$DIS_HYPERT,
      CONMED_VANCOMYCIN = row$CONMED_VANCOMYCIN,
      CONMED_COLISTIMETHATE = row$CONMED_COLISTIMETHATE
    )
}

events <- do.call(
  dplyr::bind_rows,
  lapply(seq_len(nrow(subj)), function(i) make_subject(subj[i, ]))
) |>
  dplyr::arrange(id, time, dplyr::desc(evid))

# Disjoint (id, time, evid) keys guard against the silent multi-cohort ID
# collision that merges subjects and doubles the dose.
stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod <- readModelDb("Rancic_2024_meropenem")

# Tight solver tolerances: the published volume gives a ~0.4 h half-life, so
# the trough at the end of an 8 - 12 h interval is ~1e-9 of the peak. At
# default tolerances that tail is solver noise and can go negative, which makes
# PKNCA's log-down trapezoid return NaN.
sim <- rxode2::rxSolve(
  mod, events = events,
  keep   = c("regimen", "last_dose", "tau"),
  atol   = 1e-14, rtol = 1e-12,
  returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'

# `Cc` is the IPRED (no residual error); `sim`/`ipredSim` columns carry it.
stopifnot(all(sim$Cc >= 0, na.rm = TRUE))
```

### Replicating the published clearance equation

The strictest available check on the transcription is the clearance
equation itself: it is printed in full in the Abstract and every
coefficient is tabulated in Table 3. The model’s typical-value clearance
is compared against a direct evaluation of the published equation over a
grid of covariate scenarios. This is a pure algebraic identity – the two
sides use the same numbers – so the tolerance is at machine precision,
not a cohort-derived envelope.

``` r

scen <- tidyr::expand_grid(
  CREAT                 = c(40, 94.47, 452),
  WBC                   = c(2, 7.4, 10.8, 25),
  DIS_HYPERT            = c(0, 1),
  CONMED_VANCOMYCIN     = c(0, 1),
  CONMED_COLISTIMETHATE = c(0, 1)
) |>
  dplyr::mutate(id = seq_len(dplyr::n()))

# The published equation, typed out from the Abstract of Rancic 2024.
published_cl <- function(CRE, WBC, HTA, VAN, COL) {
  5.29 * CRE^0.000001 * WBC^(-0.165) +
    0.000001 * HTA + 0.825 * VAN + 1.28 * COL
}
scen$CL_published <- published_cl(
  scen$CREAT, scen$WBC, scen$DIS_HYPERT,
  scen$CONMED_VANCOMYCIN, scen$CONMED_COLISTIMETHATE
)

mod_typ <- rxode2::zeroRe(rxode2::rxode(mod), "omega")
#> ℹ parameter labels from comments will be replaced by 'label()'
ev_typ <- scen |>
  dplyr::mutate(time = 0, amt = 1000, evid = 1L, rate = 2000, cmt = "central")
sim_typ <- rxode2::rxSolve(mod_typ, ev_typ, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: column 'CREAT' has only 'NA' values for id '1'
#> Warning: column 'WBC' has only 'NA' values for id '1'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '1'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '1'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '1'
#> Warning: column 'CREAT' has only 'NA' values for id '2'
#> Warning: column 'WBC' has only 'NA' values for id '2'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '2'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '2'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '2'
#> Warning: column 'CREAT' has only 'NA' values for id '3'
#> Warning: column 'WBC' has only 'NA' values for id '3'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '3'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '3'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '3'
#> Warning: column 'CREAT' has only 'NA' values for id '4'
#> Warning: column 'WBC' has only 'NA' values for id '4'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '4'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '4'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '4'
#> Warning: column 'CREAT' has only 'NA' values for id '5'
#> Warning: column 'WBC' has only 'NA' values for id '5'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '5'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '5'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '5'
#> Warning: column 'CREAT' has only 'NA' values for id '6'
#> Warning: column 'WBC' has only 'NA' values for id '6'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '6'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '6'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '6'
#> Warning: column 'CREAT' has only 'NA' values for id '7'
#> Warning: column 'WBC' has only 'NA' values for id '7'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '7'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '7'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '7'
#> Warning: column 'CREAT' has only 'NA' values for id '8'
#> Warning: column 'WBC' has only 'NA' values for id '8'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '8'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '8'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '8'
#> Warning: column 'CREAT' has only 'NA' values for id '9'
#> Warning: column 'WBC' has only 'NA' values for id '9'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '9'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '9'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '9'
#> Warning: column 'CREAT' has only 'NA' values for id '10'
#> Warning: column 'WBC' has only 'NA' values for id '10'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '10'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '10'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '10'
#> Warning: column 'CREAT' has only 'NA' values for id '11'
#> Warning: column 'WBC' has only 'NA' values for id '11'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '11'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '11'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '11'
#> Warning: column 'CREAT' has only 'NA' values for id '12'
#> Warning: column 'WBC' has only 'NA' values for id '12'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '12'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '12'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '12'
#> Warning: column 'CREAT' has only 'NA' values for id '13'
#> Warning: column 'WBC' has only 'NA' values for id '13'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '13'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '13'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '13'
#> Warning: column 'CREAT' has only 'NA' values for id '14'
#> Warning: column 'WBC' has only 'NA' values for id '14'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '14'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '14'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '14'
#> Warning: column 'CREAT' has only 'NA' values for id '15'
#> Warning: column 'WBC' has only 'NA' values for id '15'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '15'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '15'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '15'
#> Warning: column 'CREAT' has only 'NA' values for id '16'
#> Warning: column 'WBC' has only 'NA' values for id '16'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '16'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '16'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '16'
#> Warning: column 'CREAT' has only 'NA' values for id '17'
#> Warning: column 'WBC' has only 'NA' values for id '17'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '17'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '17'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '17'
#> Warning: column 'CREAT' has only 'NA' values for id '18'
#> Warning: column 'WBC' has only 'NA' values for id '18'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '18'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '18'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '18'
#> Warning: column 'CREAT' has only 'NA' values for id '19'
#> Warning: column 'WBC' has only 'NA' values for id '19'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '19'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '19'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '19'
#> Warning: column 'CREAT' has only 'NA' values for id '20'
#> Warning: column 'WBC' has only 'NA' values for id '20'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '20'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '20'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '20'
#> Warning: column 'CREAT' has only 'NA' values for id '21'
#> Warning: column 'WBC' has only 'NA' values for id '21'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '21'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '21'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '21'
#> Warning: column 'CREAT' has only 'NA' values for id '22'
#> Warning: column 'WBC' has only 'NA' values for id '22'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '22'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '22'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '22'
#> Warning: column 'CREAT' has only 'NA' values for id '23'
#> Warning: column 'WBC' has only 'NA' values for id '23'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '23'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '23'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '23'
#> Warning: column 'CREAT' has only 'NA' values for id '24'
#> Warning: column 'WBC' has only 'NA' values for id '24'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '24'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '24'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '24'
#> Warning: column 'CREAT' has only 'NA' values for id '25'
#> Warning: column 'WBC' has only 'NA' values for id '25'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '25'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '25'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '25'
#> Warning: column 'CREAT' has only 'NA' values for id '26'
#> Warning: column 'WBC' has only 'NA' values for id '26'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '26'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '26'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '26'
#> Warning: column 'CREAT' has only 'NA' values for id '27'
#> Warning: column 'WBC' has only 'NA' values for id '27'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '27'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '27'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '27'
#> Warning: column 'CREAT' has only 'NA' values for id '28'
#> Warning: column 'WBC' has only 'NA' values for id '28'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '28'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '28'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '28'
#> Warning: column 'CREAT' has only 'NA' values for id '29'
#> Warning: column 'WBC' has only 'NA' values for id '29'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '29'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '29'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '29'
#> Warning: column 'CREAT' has only 'NA' values for id '30'
#> Warning: column 'WBC' has only 'NA' values for id '30'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '30'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '30'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '30'
#> Warning: column 'CREAT' has only 'NA' values for id '31'
#> Warning: column 'WBC' has only 'NA' values for id '31'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '31'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '31'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '31'
#> Warning: column 'CREAT' has only 'NA' values for id '32'
#> Warning: column 'WBC' has only 'NA' values for id '32'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '32'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '32'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '32'
#> Warning: column 'CREAT' has only 'NA' values for id '33'
#> Warning: column 'WBC' has only 'NA' values for id '33'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '33'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '33'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '33'
#> Warning: column 'CREAT' has only 'NA' values for id '34'
#> Warning: column 'WBC' has only 'NA' values for id '34'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '34'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '34'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '34'
#> Warning: column 'CREAT' has only 'NA' values for id '35'
#> Warning: column 'WBC' has only 'NA' values for id '35'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '35'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '35'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '35'
#> Warning: column 'CREAT' has only 'NA' values for id '36'
#> Warning: column 'WBC' has only 'NA' values for id '36'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '36'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '36'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '36'
#> Warning: column 'CREAT' has only 'NA' values for id '37'
#> Warning: column 'WBC' has only 'NA' values for id '37'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '37'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '37'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '37'
#> Warning: column 'CREAT' has only 'NA' values for id '38'
#> Warning: column 'WBC' has only 'NA' values for id '38'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '38'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '38'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '38'
#> Warning: column 'CREAT' has only 'NA' values for id '39'
#> Warning: column 'WBC' has only 'NA' values for id '39'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '39'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '39'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '39'
#> Warning: column 'CREAT' has only 'NA' values for id '40'
#> Warning: column 'WBC' has only 'NA' values for id '40'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '40'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '40'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '40'
#> Warning: column 'CREAT' has only 'NA' values for id '41'
#> Warning: column 'WBC' has only 'NA' values for id '41'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '41'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '41'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '41'
#> Warning: column 'CREAT' has only 'NA' values for id '42'
#> Warning: column 'WBC' has only 'NA' values for id '42'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '42'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '42'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '42'
#> Warning: column 'CREAT' has only 'NA' values for id '43'
#> Warning: column 'WBC' has only 'NA' values for id '43'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '43'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '43'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '43'
#> Warning: column 'CREAT' has only 'NA' values for id '44'
#> Warning: column 'WBC' has only 'NA' values for id '44'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '44'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '44'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '44'
#> Warning: column 'CREAT' has only 'NA' values for id '45'
#> Warning: column 'WBC' has only 'NA' values for id '45'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '45'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '45'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '45'
#> Warning: column 'CREAT' has only 'NA' values for id '46'
#> Warning: column 'WBC' has only 'NA' values for id '46'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '46'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '46'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '46'
#> Warning: column 'CREAT' has only 'NA' values for id '47'
#> Warning: column 'WBC' has only 'NA' values for id '47'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '47'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '47'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '47'
#> Warning: column 'CREAT' has only 'NA' values for id '48'
#> Warning: column 'WBC' has only 'NA' values for id '48'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '48'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '48'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '48'
#> Warning: column 'CREAT' has only 'NA' values for id '49'
#> Warning: column 'WBC' has only 'NA' values for id '49'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '49'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '49'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '49'
#> Warning: column 'CREAT' has only 'NA' values for id '50'
#> Warning: column 'WBC' has only 'NA' values for id '50'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '50'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '50'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '50'
#> Warning: column 'CREAT' has only 'NA' values for id '51'
#> Warning: column 'WBC' has only 'NA' values for id '51'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '51'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '51'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '51'
#> Warning: column 'CREAT' has only 'NA' values for id '52'
#> Warning: column 'WBC' has only 'NA' values for id '52'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '52'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '52'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '52'
#> Warning: column 'CREAT' has only 'NA' values for id '53'
#> Warning: column 'WBC' has only 'NA' values for id '53'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '53'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '53'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '53'
#> Warning: column 'CREAT' has only 'NA' values for id '54'
#> Warning: column 'WBC' has only 'NA' values for id '54'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '54'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '54'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '54'
#> Warning: column 'CREAT' has only 'NA' values for id '55'
#> Warning: column 'WBC' has only 'NA' values for id '55'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '55'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '55'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '55'
#> Warning: column 'CREAT' has only 'NA' values for id '56'
#> Warning: column 'WBC' has only 'NA' values for id '56'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '56'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '56'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '56'
#> Warning: column 'CREAT' has only 'NA' values for id '57'
#> Warning: column 'WBC' has only 'NA' values for id '57'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '57'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '57'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '57'
#> Warning: column 'CREAT' has only 'NA' values for id '58'
#> Warning: column 'WBC' has only 'NA' values for id '58'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '58'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '58'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '58'
#> Warning: column 'CREAT' has only 'NA' values for id '59'
#> Warning: column 'WBC' has only 'NA' values for id '59'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '59'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '59'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '59'
#> Warning: column 'CREAT' has only 'NA' values for id '60'
#> Warning: column 'WBC' has only 'NA' values for id '60'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '60'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '60'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '60'
#> Warning: column 'CREAT' has only 'NA' values for id '61'
#> Warning: column 'WBC' has only 'NA' values for id '61'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '61'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '61'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '61'
#> Warning: column 'CREAT' has only 'NA' values for id '62'
#> Warning: column 'WBC' has only 'NA' values for id '62'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '62'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '62'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '62'
#> Warning: column 'CREAT' has only 'NA' values for id '63'
#> Warning: column 'WBC' has only 'NA' values for id '63'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '63'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '63'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '63'
#> Warning: column 'CREAT' has only 'NA' values for id '64'
#> Warning: column 'WBC' has only 'NA' values for id '64'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '64'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '64'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '64'
#> Warning: column 'CREAT' has only 'NA' values for id '65'
#> Warning: column 'WBC' has only 'NA' values for id '65'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '65'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '65'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '65'
#> Warning: column 'CREAT' has only 'NA' values for id '66'
#> Warning: column 'WBC' has only 'NA' values for id '66'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '66'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '66'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '66'
#> Warning: column 'CREAT' has only 'NA' values for id '67'
#> Warning: column 'WBC' has only 'NA' values for id '67'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '67'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '67'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '67'
#> Warning: column 'CREAT' has only 'NA' values for id '68'
#> Warning: column 'WBC' has only 'NA' values for id '68'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '68'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '68'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '68'
#> Warning: column 'CREAT' has only 'NA' values for id '69'
#> Warning: column 'WBC' has only 'NA' values for id '69'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '69'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '69'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '69'
#> Warning: column 'CREAT' has only 'NA' values for id '70'
#> Warning: column 'WBC' has only 'NA' values for id '70'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '70'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '70'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '70'
#> Warning: column 'CREAT' has only 'NA' values for id '71'
#> Warning: column 'WBC' has only 'NA' values for id '71'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '71'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '71'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '71'
#> Warning: column 'CREAT' has only 'NA' values for id '72'
#> Warning: column 'WBC' has only 'NA' values for id '72'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '72'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '72'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '72'
#> Warning: column 'CREAT' has only 'NA' values for id '73'
#> Warning: column 'WBC' has only 'NA' values for id '73'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '73'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '73'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '73'
#> Warning: column 'CREAT' has only 'NA' values for id '74'
#> Warning: column 'WBC' has only 'NA' values for id '74'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '74'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '74'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '74'
#> Warning: column 'CREAT' has only 'NA' values for id '75'
#> Warning: column 'WBC' has only 'NA' values for id '75'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '75'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '75'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '75'
#> Warning: column 'CREAT' has only 'NA' values for id '76'
#> Warning: column 'WBC' has only 'NA' values for id '76'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '76'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '76'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '76'
#> Warning: column 'CREAT' has only 'NA' values for id '77'
#> Warning: column 'WBC' has only 'NA' values for id '77'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '77'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '77'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '77'
#> Warning: column 'CREAT' has only 'NA' values for id '78'
#> Warning: column 'WBC' has only 'NA' values for id '78'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '78'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '78'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '78'
#> Warning: column 'CREAT' has only 'NA' values for id '79'
#> Warning: column 'WBC' has only 'NA' values for id '79'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '79'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '79'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '79'
#> Warning: column 'CREAT' has only 'NA' values for id '80'
#> Warning: column 'WBC' has only 'NA' values for id '80'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '80'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '80'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '80'
#> Warning: column 'CREAT' has only 'NA' values for id '81'
#> Warning: column 'WBC' has only 'NA' values for id '81'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '81'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '81'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '81'
#> Warning: column 'CREAT' has only 'NA' values for id '82'
#> Warning: column 'WBC' has only 'NA' values for id '82'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '82'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '82'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '82'
#> Warning: column 'CREAT' has only 'NA' values for id '83'
#> Warning: column 'WBC' has only 'NA' values for id '83'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '83'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '83'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '83'
#> Warning: column 'CREAT' has only 'NA' values for id '84'
#> Warning: column 'WBC' has only 'NA' values for id '84'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '84'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '84'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '84'
#> Warning: column 'CREAT' has only 'NA' values for id '85'
#> Warning: column 'WBC' has only 'NA' values for id '85'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '85'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '85'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '85'
#> Warning: column 'CREAT' has only 'NA' values for id '86'
#> Warning: column 'WBC' has only 'NA' values for id '86'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '86'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '86'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '86'
#> Warning: column 'CREAT' has only 'NA' values for id '87'
#> Warning: column 'WBC' has only 'NA' values for id '87'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '87'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '87'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '87'
#> Warning: column 'CREAT' has only 'NA' values for id '88'
#> Warning: column 'WBC' has only 'NA' values for id '88'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '88'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '88'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '88'
#> Warning: column 'CREAT' has only 'NA' values for id '89'
#> Warning: column 'WBC' has only 'NA' values for id '89'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '89'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '89'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '89'
#> Warning: column 'CREAT' has only 'NA' values for id '90'
#> Warning: column 'WBC' has only 'NA' values for id '90'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '90'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '90'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '90'
#> Warning: column 'CREAT' has only 'NA' values for id '91'
#> Warning: column 'WBC' has only 'NA' values for id '91'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '91'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '91'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '91'
#> Warning: column 'CREAT' has only 'NA' values for id '92'
#> Warning: column 'WBC' has only 'NA' values for id '92'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '92'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '92'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '92'
#> Warning: column 'CREAT' has only 'NA' values for id '93'
#> Warning: column 'WBC' has only 'NA' values for id '93'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '93'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '93'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '93'
#> Warning: column 'CREAT' has only 'NA' values for id '94'
#> Warning: column 'WBC' has only 'NA' values for id '94'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '94'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '94'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '94'
#> Warning: column 'CREAT' has only 'NA' values for id '95'
#> Warning: column 'WBC' has only 'NA' values for id '95'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '95'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '95'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '95'
#> Warning: column 'CREAT' has only 'NA' values for id '96'
#> Warning: column 'WBC' has only 'NA' values for id '96'
#> Warning: column 'DIS_HYPERT' has only 'NA' values for id '96'
#> Warning: column 'CONMED_VANCOMYCIN' has only 'NA' values for id '96'
#> Warning: column 'CONMED_COLISTIMETHATE' has only 'NA' values for id '96'

cl_model <- sim_typ |>
  dplyr::group_by(id) |>
  dplyr::summarise(CL_model = dplyr::first(cl), V_model = dplyr::first(vc),
                   .groups = "drop")

cmp_cl <- dplyr::left_join(scen, cl_model, by = "id")
rel_err <- abs(cmp_cl$CL_model - cmp_cl$CL_published) / cmp_cl$CL_published

cat(sprintf(
  "Clearance-equation replication over %d covariate scenarios: max relative error %.3g\n",
  nrow(cmp_cl), max(rel_err)
))
#> Clearance-equation replication over 96 covariate scenarios: max relative error 2.14e-15
# Algebraic identity, not a simulated statistic: machine precision is correct.
stopifnot(max(rel_err) < 1e-10)

# The volume carries no covariates, so it must be exactly 2.05 L everywhere.
stopifnot(max(abs(cmp_cl$V_model - 2.05)) < 1e-10)
```

``` r

cmp_cl |>
  dplyr::filter(CREAT == 94.47, WBC %in% c(7.4, 10.8)) |>
  dplyr::transmute(
    `CREAT (umol/L)`      = CREAT,
    `WBC (10^9/L)`        = WBC,
    HTA                   = DIS_HYPERT,
    VAN                   = CONMED_VANCOMYCIN,
    COL                   = CONMED_COLISTIMETHATE,
    `CL published (L/h)`  = round(CL_published, 4),
    `CL model (L/h)`      = round(CL_model, 4)
  ) |>
  knitr::kable(
    caption = paste(
      "Typical-value clearance from the packaged model against a direct",
      "evaluation of the Rancic 2024 Abstract equation, at the cohort-mean",
      "creatinine. The CREAT exponent of 1e-6 makes the creatinine term",
      "1.0000 to four decimal places at every plausible value, which is why",
      "the covariate is numerically inert."
    )
  )
```

| CREAT (umol/L) | WBC (10^9/L) | HTA | VAN | COL | CL published (L/h) | CL model (L/h) |
|---------------:|-------------:|----:|----:|----:|-------------------:|---------------:|
|          94.47 |          7.4 |   0 |   0 |   0 |             3.8022 |         3.8022 |
|          94.47 |          7.4 |   0 |   0 |   1 |             5.0822 |         5.0822 |
|          94.47 |          7.4 |   0 |   1 |   0 |             4.6272 |         4.6272 |
|          94.47 |          7.4 |   0 |   1 |   1 |             5.9072 |         5.9072 |
|          94.47 |          7.4 |   1 |   0 |   0 |             3.8022 |         3.8022 |
|          94.47 |          7.4 |   1 |   0 |   1 |             5.0822 |         5.0822 |
|          94.47 |          7.4 |   1 |   1 |   0 |             4.6272 |         4.6272 |
|          94.47 |          7.4 |   1 |   1 |   1 |             5.9072 |         5.9072 |
|          94.47 |         10.8 |   0 |   0 |   0 |             3.5723 |         3.5723 |
|          94.47 |         10.8 |   0 |   0 |   1 |             4.8523 |         4.8523 |
|          94.47 |         10.8 |   0 |   1 |   0 |             4.3973 |         4.3973 |
|          94.47 |         10.8 |   0 |   1 |   1 |             5.6773 |         5.6773 |
|          94.47 |         10.8 |   1 |   0 |   0 |             3.5723 |         3.5723 |
|          94.47 |         10.8 |   1 |   0 |   1 |             4.8523 |         4.8523 |
|          94.47 |         10.8 |   1 |   1 |   0 |             4.3973 |         4.3973 |
|          94.47 |         10.8 |   1 |   1 |   1 |             5.6773 |         5.6773 |

Typical-value clearance from the packaged model against a direct
evaluation of the Rancic 2024 Abstract equation, at the cohort-mean
creatinine. The CREAT exponent of 1e-6 makes the creatinine term 1.0000
to four decimal places at every plausible value, which is why the
covariate is numerically inert. {.table}

Two facts about the covariate model are visible in that table. First,
the creatinine and hypertension effects are estimated at (essentially)
zero and change clearance by less than one part in 10^5, even though
both survived backward deletion at p \< 0.01. Second, the vancomycin and
colistimethate effects are *additive* shifts, so at a typical clearance
near 3.6 L/h they are worth +23% and +36% respectively, and +59%
together.

## Steady-state profiles

``` r

# rxSolve returns observation records only (dose rows are not echoed back), so
# every row here is already an observation.
sim_int <- sim |>
  dplyr::mutate(tad = time - last_dose)

sim_int |>
  dplyr::group_by(regimen, tad) |>
  dplyr::summarise(
    Q05 = stats::quantile(Cc, 0.05), Q50 = stats::median(Cc),
    Q95 = stats::quantile(Cc, 0.95), .groups = "drop"
  ) |>
  ggplot(aes(tad, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = c(40.69, 12.55), linetype = "dashed", colour = "firebrick") +
  facet_wrap(~regimen) +
  scale_y_log10() +
  labs(
    x = "Time after the last dose (h)", y = "Meropenem Cc (mg/L, log scale)",
    caption = "Dashed lines: mean observed C1 and C2 from Rancic 2024 Table 1."
  ) +
  theme_bw()
```

![Simulated steady-state meropenem profiles over the final dosing
interval, by regimen. Ribbons are the 5th - 95th percentiles across the
101-subject cohort; the dashed lines mark the mean concentrations the
paper actually observed (40.69 mg/L 5 - 30 min after the end of the
infusion, 12.55 mg/L 3 - 4 h after
it).](Rancic_2024_meropenem_files/figure-html/profiles-1.png)

Simulated steady-state meropenem profiles over the final dosing
interval, by regimen. Ribbons are the 5th - 95th percentiles across the
101-subject cohort; the dashed lines mark the mean concentrations the
paper actually observed (40.69 mg/L 5 - 30 min after the end of the
infusion, 12.55 mg/L 3 - 4 h after it).

The predicted profile is far more peaked than the observed
concentrations: it overshoots the observed peak and has already fallen
below the observed 3 - 4 h concentration by about an hour. That is the
signature of an under-estimated volume with a correctly estimated
clearance, and it is quantified in the Errata below.

## PKNCA validation

``` r

# Only `!is.na(Cc)` -- adding `time > 0` or `Cc > 0` would drop the time-zero
# row that anchors AUC.
sim_nca <- sim_int |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, regimen, tau, time = tad, Cc)

# Time-zero row is produced by the observation grid itself (obs_grid starts at
# 0); assert it rather than assume it.
stopifnot(
  sim_nca |> dplyr::group_by(id) |> dplyr::summarise(has0 = any(time == 0)) |>
    dplyr::pull(has0) |> all()
)

conc_obj <- PKNCA::PKNCAconc(as.data.frame(sim_nca), Cc ~ time | regimen + id)

dose_df <- sim_int |>
  dplyr::distinct(id, regimen) |>
  dplyr::left_join(subj |> dplyr::select(id, amt), by = "id") |>
  dplyr::mutate(time = 0)
dose_obj <- PKNCA::PKNCAdose(as.data.frame(dose_df), amt ~ time | regimen + id)

# One interval per regimen, spanning that regimen's dosing interval.
intervals <- regimens |>
  dplyr::transmute(
    regimen, start = 0, end = tau,
    cmax = TRUE, tmax = TRUE, auclast = TRUE, ctrough = TRUE, half.life = TRUE
  ) |>
  as.data.frame()

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))

nca_wide <- as.data.frame(nca_res) |>
  dplyr::select(regimen, id, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)
stopifnot(nrow(nca_wide) == n_sub, !anyNA(nca_wide$auclast), !anyNA(nca_wide$half.life))
```

### Internal identities

These two checks compare the NCA output against closed forms built from
each subject’s *own* simulated `cl` and `vc`. Both sides use the same
drawn parameters, so the only difference is numerical (trapezoidal and
log-linear regression error) and a tight bound is the correct assertion
– these are not cohort statistics whose spread depends on which etas
were drawn.

``` r

per_subj <- sim_int |>
  dplyr::group_by(id, regimen, tau) |>
  dplyr::summarise(cl = dplyr::first(cl), vc = dplyr::first(vc), .groups = "drop") |>
  dplyr::left_join(subj |> dplyr::select(id, amt), by = "id") |>
  dplyr::left_join(nca_wide, by = c("id", "regimen")) |>
  dplyr::mutate(
    auc_closed  = amt / cl,                 # steady-state AUC over one interval
    thalf_closed = log(2) * vc / cl,
    auc_pct     = 100 * (auclast - auc_closed) / auc_closed,
    thalf_pct   = 100 * (half.life - thalf_closed) / thalf_closed,
    accum       = 1 / (1 - exp(-(cl / vc) * tau))
  )

cat(sprintf("AUC over the interval vs Dose/CL: median %.4f%%, max |diff| %.4f%%\n",
            stats::median(per_subj$auc_pct), max(abs(per_subj$auc_pct))))
#> AUC over the interval vs Dose/CL: median -0.0074%, max |diff| 0.0190%
cat(sprintf("NCA half-life vs log(2)*V/CL:     median %.4f%%, max |diff| %.4f%%\n",
            stats::median(per_subj$thalf_pct), max(abs(per_subj$thalf_pct))))
#> NCA half-life vs log(2)*V/CL:     median 0.0000%, max |diff| 0.0000%
cat(sprintf("Steady-state accumulation factor: max %.8f (1.0 = no accumulation)\n",
            max(per_subj$accum)))
#> Steady-state accumulation factor: max 1.00001773 (1.0 = no accumulation)

stopifnot(
  # Trapezoidal / log-down error on a dense grid over a mono-exponential decay.
  # Realised max 0.017%; the error scales with (k*h)^3 and k varies with the
  # drawn etas (14.7% CV), so 0.2% leaves an order of magnitude of headroom.
  max(abs(per_subj$auc_pct)) < 0.2,
  # Log-linear terminal fit against the analytic half-life. A log-linear
  # regression on a mono-exponential decay recovers log(2)/k exactly whatever k
  # is, so this realises 0.0000% and does not widen with the eta draw.
  max(abs(per_subj$thalf_pct)) < 0.05,
  # With a ~0.4 h half-life and an 8 - 12 h interval, the model accumulates not
  # at all, so a single interval IS the steady-state interval. Realised
  # ~1 + 1e-5 across the cohort; 1.001 leaves three orders of magnitude of
  # headroom over the eta spread while still going red for any volume large
  # enough to matter (V = 21 L would give an accumulation factor of ~1.35).
  max(per_subj$accum) < 1.001
)
```

``` r

per_subj |>
  dplyr::group_by(regimen) |>
  dplyr::summarise(
    n              = dplyr::n(),
    `Cmax (mg/L)`  = round(stats::median(cmax), 1),
    `Tmax (h)`     = round(stats::median(tmax), 2),
    `AUCtau (mg*h/L)` = round(stats::median(auclast), 1),
    `Ctrough (mg/L)`  = signif(stats::median(ctrough), 3),
    `t1/2 (h)`     = round(stats::median(half.life), 3),
    `CL (L/h)`     = round(stats::median(cl), 2),
    .groups = "drop"
  ) |>
  knitr::kable(caption = "Simulated steady-state NCA by regimen (median across the cohort).")
```

| regimen | n | Cmax (mg/L) | Tmax (h) | AUCtau (mg\*h/L) | Ctrough (mg/L) | t1/2 (h) | CL (L/h) |
|:---|---:|---:|---:|---:|---:|---:|---:|
| 1000 mg q12h | 15 | 306.9 | 0.5 | 241.2 | 0.00e+00 | 0.343 | 4.15 |
| 1000 mg q8h | 76 | 318.2 | 0.5 | 263.5 | 2.97e-04 | 0.374 | 3.79 |
| 2000 mg q12h | 7 | 657.8 | 0.5 | 575.3 | 2.20e-06 | 0.409 | 3.48 |
| 2000 mg q8h | 3 | 667.3 | 0.5 | 598.8 | 3.30e-03 | 0.425 | 3.34 |

Simulated steady-state NCA by regimen (median across the cohort).
{.table style="width:100%;"}

### Comparison against the published concentrations

Rancic 2024 reports no NCA parameters, but Table 1 reports the two
observed concentrations per patient, and the first of them – drawn 5 -
30 min after the end of the infusion – is a direct observation of Cmax
for an intravenous infusion. It is compared against the simulated Cmax
below.

``` r

# Table 1: C1 = 40.69 +/- 16.67 mg/L, drawn 5 - 30 min after the end of the
# infusion, pooled over all regimens.
published <- tibble::tibble(regimen = regimens$regimen, cmax = 40.69)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_res,
  reference     = published,
  by            = "regimen",
  params        = "cmax",
  units         = c(cmax = "mg/L"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = paste(
    "Simulated steady-state Cmax against the mean concentration Rancic 2024",
    "observed 5 - 30 min after the end of the infusion (Table 1).",
    "* marks a difference above 20%. Every row is starred; this is a known,",
    "reproducible deviation of the published model from the paper's own data",
    "and is analysed in the Errata below. No parameter was tuned."
  )
)
```

| NCA parameter | regimen      | Reference | Simulated | % diff     |
|:--------------|:-------------|:----------|:----------|:-----------|
| Cmax (mg/L)   | 1000 mg q12h | 40.7      | 307       | +654.3%\*  |
| Cmax (mg/L)   | 1000 mg q8h  | 40.7      | 318       | +681.9%\*  |
| Cmax (mg/L)   | 2000 mg q12h | 40.7      | 658       | +1516.6%\* |
| Cmax (mg/L)   | 2000 mg q8h  | 40.7      | 667       | +1540.1%\* |

Simulated steady-state Cmax against the mean concentration Rancic 2024
observed 5 - 30 min after the end of the infusion (Table 1). \* marks a
difference above 20%. Every row is starred; this is a known,
reproducible deviation of the published model from the paper’s own data
and is analysed in the Errata below. No parameter was tuned. {.table}

``` r

# NOTE ON THE DIRECTION OF THE ARGUMENT. For a one-compartment model the peak
# is Dose/V, so a SMALL volume makes concentrations too HIGH. The observed
# peaks therefore bound V from ABOVE, not below: C <= Dose/V rearranges to
# V <= Dose/C, i.e. the largest observed peak gives V <= 1000/88.95 = 11.2 L,
# which 2.05 L satisfies. A ceiling argument does NOT falsify the published
# volume, and it is worth being explicit about that because the inequality is
# easy to state backwards.
#
# What does falsify 2.05 L is that the model then MISSES the paper's own data
# in two independent, direction-safe ways.

# Typical clearance at the cohort-average covariates (this reproduces the
# paper's own base-model value of 3.80 L/h).
cl_typ <- 5.29 * 10.8^(-0.165) + 0.825 * 0.168 + 1.28 * 0.069

# (a) DECAY RATE -- the strongest check, because the observed side never
#     references V. Two samples per patient pin the elimination rate constant
#     directly: kel = log(C1/C2) / dt. Table 1 gives C1 = 40.69 mg/L and
#     C2 = 12.55 mg/L; Methods 2.2 draws the first 5 - 30 min after the end of
#     the infusion and the second 3 - 4 h after it, so dt spans 2.5 - 3.92 h.
#     The SHORTEST dt gives the fastest decay and hence the SMALLEST implied
#     volume, so it is the conservative end of the range.
dt_obs       <- c(short = 3 - 0.5, long = 4 - 1 / 12)
kel_obs      <- log(40.69 / 12.55) / dt_obs
v_from_decay <- cl_typ / kel_obs

# (b) PEAK MAGNITUDE -- to REPRODUCE (not merely exceed) the mean observed peak
#     at the smallest dose used, the model needs V = Dose / C.
v_from_peak <- 1000 / 40.69

cat(sprintf("Published V                                 : %.2f L\n", 2.05))
#> Published V                                 : 2.05 L
cat(sprintf("V implied by the observed C1/C2 decay       : %.1f - %.1f L\n",
            min(v_from_decay), max(v_from_decay)))
#> V implied by the observed C1/C2 decay       : 8.1 - 12.7 L
cat(sprintf("V implied by the mean observed peak         : %.1f L\n", v_from_peak))
#> V implied by the mean observed peak         : 24.6 L
cat(sprintf("Literature value cited in the Introduction  : %d L\n", 21L))
#> Literature value cited in the Introduction  : 21 L
cat(sprintf("Model half-life, log(2)*V/CL                : %.2f h\n",
            log(2) * 2.05 / cl_typ))
#> Model half-life, log(2)*V/CL                : 0.37 h
cat(sprintf("Half-life implied by the observed decay     : %.2f - %.2f h\n",
            min(log(2) / kel_obs), max(log(2) / kel_obs)))
#> Half-life implied by the observed decay     : 1.47 - 2.31 h
cat(sprintf("Simulated median Cmax, 1000 mg q8h          : %.1f mg/L\n",
            stats::median(per_subj$cmax[per_subj$regimen == "1000 mg q8h"])))
#> Simulated median Cmax, 1000 mg q8h          : 318.2 mg/L

# Deterministic: every number above comes from the paper's own printed values,
# none from the simulated cohort. Both independent routes put the volume at
# least 3.5-fold above the published 2.05 L, at the conservative end of the
# sampling-time range.
stopifnot(
  min(v_from_decay) > 3.5 * 2.05,
  v_from_peak       > 3.5 * 2.05
)
```

## Assumptions and deviations

- **Leukocyte-count units and distribution.** The paper never states the
  units of `WBC` and Table 1 has no WBC row, yet `WBC^(-0.165)` is one
  of only two multiplicative covariate terms in the final model. The
  unit was settled by internal consistency rather than assumed: solving
  `5.29 * WBC^(-0.165) + 0.23 = 3.80` (the base model’s typical
  clearance, with 0.23 the cohort-average contribution of the additive
  vancomycin and colistimethate terms at their Table 1 prevalences)
  gives `WBC = 10.8 x 10^9/L`, a typical leukocytosis. On a cells/uL
  scale the same arithmetic would need `WBC ~ 10,800` and would give a
  typical clearance of 1.13 L/h – less than a third of the base-model
  value – so cells/uL is excluded. The simulated distribution
  (log-normal, median 10.8, 45% CV, truncated to 2 - 40) is an
  assumption of this vignette, not a published one.
- **Infusion duration.** Methods 2.2 says only “intermittent intravenous
  infusion”. A 30-minute infusion is assumed here, consistent with the
  first sample being drawn 5 - 30 min after its end. The assumption
  affects the simulated Cmax by only a few percent at the published
  half-life of 0.4 h; it does not affect AUC, clearance, or any
  assertion in this vignette.
- **Dosing-regimen mix.** The paper reports the daily-dose mean, SD and
  range but not the mix of 1,000/2,000 mg and q8h/q12h. The 15/76/7/3
  mix used here was chosen to reproduce the published mean and SD, which
  the cohort chunk asserts.
- **Creatinine distribution.** A log-normal matched on the Table 1 mean
  and SD and truncated to the reported 40 - 452 umol/L range. Because
  the fitted creatinine exponent is 1e-6, this choice has no effect on
  any simulated quantity.
- **Screened-but-dropped comorbidities.** Chronic renal failure,
  cerebral infarction, pneumonia and urinary tract infection are
  documented in the covariate-screen table above rather than given
  canonical covariate names, for the reason set out in that section.

## Errata and internal inconsistencies in the source

The paper is internally inconsistent in four places. All four are
recorded faithfully in the model file rather than corrected, and none
was tuned.

**1. The published volume of distribution is falsified by the paper’s
own data.** Table 3 reports `V = 2.05 L` (95% CI 2.039 - 2.061) and the
bootstrap in Table 5 reports 1.809 L, so the value is reproducible
within the paper and is not a single-cell typo. The base model’s 3.52 L
is of the same order, so the discrepancy is systemic rather than a Table
3 misprint.

It is worth being precise about *why* 2.05 L is wrong, because the
obvious argument runs the wrong way. For a one-compartment model the
peak is `Dose/V`, so a small volume makes concentrations too **high**;
`C <= Dose/V` rearranges to `V <= Dose/C`, meaning the observed peaks
bound the volume from **above** (the largest reported peak, 88.95 mg/L
at the smallest dose, gives `V <= 11.2 L`). A ceiling argument therefore
does *not* rule out 2.05 L. What rules it out is that the model, at that
volume, misses the paper’s own Table 1 in two independent directions:

- **Decay rate.** The two samples per patient pin the elimination rate
  constant without reference to `V` at all: `kel = log(C1/C2)/dt`. With
  `C1 = 40.69 mg/L`, `C2 = 12.55 mg/L` and the 2.5 - 3.9 h sampling gap
  of Methods 2.2, the observed half-life is 1.7 - 2.7 h and the implied
  volume `CL/kel` is 8 - 13 L. The published model gives
  `log(2)*V/CL = 0.4 h`, four to five-fold too fast.
- **Peak magnitude.** Reproducing (not merely exceeding) the mean
  observed peak of 40.69 mg/L at the smallest dose used requires
  `V = 24.6 L`.

Both routes put the volume several-fold above 2.05 L, and the literature
value the paper’s own Introduction cites is 21 L. The consequence is
visible in the figure and in the NCA comparison: simulated Cmax is
roughly eight-fold above the observed mean peak, the simulated half-life
is 0.4 h instead of the 1.7 - 2.7 h the observed C1/C2 pair implies, and
the simulated 4-hour concentration is far below the observed one.
**Clearance is unaffected** – `Dose/(CL*tau)` at the typical clearance
gives an average steady-state concentration of about 35 mg/L against an
observed mean of the two samples of 26.6 mg/L – so the model remains
usable for clearance-driven quantities (AUC, total exposure, average
concentration) and should not be used for peak, trough, or
time-above-MIC without substituting a defensible volume.

**2. The sign of the leukocyte effect contradicts the paper’s own
prose.** Table 3 reports the WBC exponent as -0.165 with a 95% CI of
-0.219 to -0.110, i.e. clearance *decreases* as the leukocyte count
rises, and the Abstract’s equation prints `WBCs^(-0.165)`. The
Discussion and Conclusion state the opposite: WBC is listed among “the
factors that are significantly associated with *increased* clearance”.
The signed estimate and its confidence interval are encoded, per the
standing rule that a printed equation beats prose.

**3. The variability terms are reported as variances but restated as
percentages.** Table 3 rows read “Inter-individual variance of CL
(omega^2 CL) = 0.0215” and “Residual error variance (sigma^2 CL) =
0.44”, both with standard errors and confidence limits, while the
Results prose says “Inter- and intra-individual variability were 2.15
and 44%, respectively” – i.e. it multiplies each variance by 100 and
calls the product a percentage. The table’s explicit `omega^2` /
`sigma^2` notation is taken as primary, so the model carries
`etalcl ~ 0.0215` (a 14.7% CV) and a residual SD of
`sqrt(0.44) = 0.6633`.

The *form* of the residual needed a further decision. Methods 2.3
describes an “additive” residual on concentration, but an additive SD of
0.66 mg/L on the linear scale is 1.6% of the mean observed peak and
cannot generate the paper’s own final-model RMSPE of 15.86 (Table 4) or
its observed between-patient CVs of 41% on C1 and 61% on C2 (Table 1);
it is roughly 25-fold too small. In NONMEM, an “additive” residual
applied to log-transformed observations is exactly a proportional
residual once back-transformed to the linear concentration scale, and
that is how it is encoded here (`Cc ~ prop(0.6633)`). A user who wants
the literal linear-scale additive reading can substitute
`Cc ~ add(0.6633)` in a copy of the model file.

**4. The bootstrap does not reproduce the final model.** Table 5 reports
the bootstrap means as “comparable to” the final model, but they are
not: the clearance intercept is 3.864 against 5.29, the volume 1.809
against 2.05, and the WBC effect `+0.001` against `-0.165` – a sign
change on the covariate that carries most of the covariate model. The
bootstrap also reports five effects on the volume of distribution
(`theta16` - `theta20`) for which the final model reports no estimates
at all. Only the Table 3 final-model estimates are encoded; the
bootstrap column is not a second model and is not used.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3         tidyr_1.3.2           dplyr_1.2.1          
#> [4] rxode2_5.1.6          PKNCA_0.12.1          nlmixr2lib_0.3.2.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.60           bslib_0.12.0       
#>  [4] lattice_0.22-9      vctrs_0.7.3         tools_4.6.1        
#>  [7] generics_0.1.4      parallel_4.6.1      tibble_3.3.1       
#> [10] symengine_0.2.13    pkgconfig_2.0.3     data.table_1.18.6.1
#> [13] checkmate_2.3.4     RColorBrewer_1.1-3  S7_0.2.2           
#> [16] desc_1.4.3          RcppParallel_6.2.1  lifecycle_1.0.5    
#> [19] compiler_4.6.1      farver_2.1.2        textshaping_1.0.5  
#> [22] fontawesome_0.5.3   htmltools_0.5.9     sys_3.4.3          
#> [25] sass_0.4.10         yaml_2.3.12         pillar_1.11.1      
#> [28] pkgdown_2.2.1       crayon_1.5.3        jquerylib_0.1.4    
#> [31] whisker_0.4.1       openssl_2.4.2       cachem_1.1.0       
#> [34] nlme_3.1-169        tidyselect_1.2.1    digest_0.6.39      
#> [37] lotri_1.0.4         purrr_1.2.2         labeling_0.4.3     
#> [40] rxode2ll_2.0.16     fastmap_1.2.0       grid_4.6.1         
#> [43] cli_3.6.6           dparser_1.3.1-13    magrittr_2.0.5     
#> [46] withr_3.0.3         scales_1.4.0        backports_1.5.1    
#> [49] rmarkdown_2.32      otel_0.2.0          askpass_1.2.1      
#> [52] ragg_1.5.2          memoise_2.0.1       evaluate_1.0.5     
#> [55] knitr_1.52          rex_1.2.2           PreciseSums_0.7    
#> [58] rlang_1.3.0         downlit_0.4.5       Rcpp_1.1.2         
#> [61] glue_1.8.1          xml2_1.6.0          jsonlite_2.0.0     
#> [64] R6_2.6.1            systemfonts_1.3.2   fs_2.1.0
```
