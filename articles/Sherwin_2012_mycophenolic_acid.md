# Mycophenolic acid (Sherwin 2012)

``` r

library(nlmixr2lib)
library(rxode2)
#> rxode2 5.1.2 using 2 threads (see ?getRxThreads)
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
library(PKNCA)
#> 
#> Attaching package: 'PKNCA'
#> The following object is masked from 'package:stats':
#> 
#>     filter
```

## Mycophenolic acid (MPA) with enterohepatic recirculation in pediatric SLE

Mycophenolic acid (MPA) is the active immunosuppressant moiety of
mycophenolate mofetil (MMF), administered orally to suppress lymphocyte
proliferation in autoimmune disease and after solid-organ
transplantation. MPA is converted into the inactive metabolite
7-O-MPA-glucuronide (MPAG) by UDP-glucuronosyltransferases, and MPAG
undergoes enterohepatic recirculation: it is excreted into the gut lumen
via biliary transport (MRP2), de-glucuronidated back to MPA by
intestinal bacteria, and reabsorbed. The resulting secondary absorption
peak between 4 and 9 h post-dose can contribute as much as 40% (range
10-60%) to the systemic MPA AUC.

Sherwin et al. (2012) developed the first population PK model with
enterohepatic recirculation specifically in pediatric and adolescent
patients with childhood-onset systemic lupus erythematosus (cSLE). The
cohort comprised 19 outpatients on a stable oral MMF regimen (typical
dose 600 mg/m^2 BID), with PK profiles drawn on a single fasting visit
and standardized meals at +1 h and +4 h post-dose. The final model has
the following structure:

- A series of `NN = 8.2` transit compartments (Savic 2007
  parameterisation) delivers the MMF / MPA-equivalent dose into a gut
  compartment with mean transit time `MTT = 1.1 h`;
- The gut compartment empties into the MPA central compartment at a
  fixed first-order rate `Ka = 1.5 1/h`;
- MPA disposition is two-compartment (`CL1/F = 25.3 L/h`,
  `V3/F = 20.9 L`, `CL2/F = 19.8 L/h`, `V4/F = 234 L`);
- Of the total MPA elimination, a fraction `FM = 0.85` (fixed) is
  converted to MPAG and enters a one-compartment MPAG central pool
  (`V_M = V3 = 20.9 L`; apparent renal `CLM/F = 2.5 L/h`); the remaining
  0.15 represents AcMPAG formation and is not carried in the model;
- Of the total MPAG elimination, a fraction `FMPAG = 0.65` (fixed) is
  biliary and accumulates in a gallbladder compartment; the complement
  is renal;
- The gallbladder empties to the gut at a fast rate during two fixed
  windows post-dose (1-2 h and 4-6 h, matching the study meals); only
  `EHC = 0.35` of the emptied content (fixed) re-enters the gut for
  re-absorption as MPA, generating the secondary peak.

The packaged model produces two outputs: `Cc` (MPA plasma concentration
in mg/L) and `Cc_mpag` (MPAG plasma concentration in mg/L). No
covariates were retained in the final model; bodyweight, age, sex, race,
ethnicity, and disease duration were all screened and rejected (see
`population` and `covariatesDataExcluded` metadata on the model and the
Assumptions and deviations section below).

- Citation: Sherwin CMT, Sagcal-Gironella ACP, Fukuda T, Brunner HI,
  Vinks AA. Development of population PK model with enterohepatic
  circulation for mycophenolic acid in patients with childhood-onset
  systemic lupus erythematosus. British Journal of Clinical
  Pharmacology. 2012;73(5):727-740.
  <doi:10.1111/j.1365-2125.2011.04140.x>.
- Article: <https://doi.org/10.1111/j.1365-2125.2011.04140.x>

## Population

The model was fit to 186 MPA + 186 MPAG plasma concentrations (372
total) from 19 cSLE outpatients (Sherwin 2012 Table 1):

| Characteristic | Value |
|----|----|
| Age (years, mean +/- SD) | 16.9 +/- 4 (range 10.6-28.2) |
| Children (2-12 y) / Adolescents (12-21 y) / Adults (\>21 y) | 2 / 14 / 3 |
| Weight (kg, mean +/- SD) | 66.6 +/- 15 (range 43.4-103) |
| SLE disease duration (y) | 3.3 +/- 3 (range 0.2-12.8) |
| MMF treatment duration (y) | 1.5 +/- 1.3 (range 0.14-6.4) |
| MMF daily dose (mg, mean +/- SD) | 1973 +/- 634 (range 1000-3000) |
| Female / Male | 18 / 1 (95% / 5%) |
| African American / Caucasian | 11 / 8 (58% / 42%) |
| Hispanic / Non-Hispanic | 4 / 15 (21% / 79%) |
| Subjects on prednisone | 18 / 19 (95%) |

PK samples were drawn at pre-dose, 20 min, 40 min, 1, 1.5, 2, 3, 4, 6,
and 9 h post-dose on a single fasting visit. Standardized meals were
given at +1 h and +4 h post-dose to time the bile-release windows. The
same metadata is available programmatically through
`readModelDb("Sherwin_2012_mycophenolic_acid")$population`.

## Source trace

Per-parameter origin is recorded as in-file comments next to each
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) entry in
`inst/modeldb/specificDrugs/Sherwin_2012_mycophenolic_acid.R`. The table
collects them for review.

| Equation / parameter | Value (paper) | Value (file) | Source |
|----|----|----|----|
| `lka` (Ka, fixed) | 1.5 1/h | `fixed(log(1.5))` | Table 3 (fixed) |
| `lcl` (CL1 MPA / F) | 25.3 L/h | `log(25.3)` | Table 3 |
| `lvc` (V3 MPA / F = V_MPAG) | 20.9 L | `log(20.9)` | Table 3 |
| `lq` (CL2 MPA / F) | 19.8 L/h | `log(19.8)` | Table 3 |
| `lvp` (V4 MPA / F) | 234 L | `log(234)` | Table 3 |
| `lcl_mpag` (CLM MPAG / F, apparent renal) | 2.5 L/h | `log(2.5)` | Table 3 |
| `lmtt` (Mean transit time MTT) | 1.1 h | `log(1.1)` | Table 3 |
| `lnn` (Number of transit compartments NN) | 8.2 | `log(8.2)` | Table 3 |
| `e_fm` (FM, fixed) | 0.85 | `fixed(0.85)` | Table 3 (fixed at 85%) |
| `e_fmpag` (FMPAG, fixed) | 0.65 | `fixed(0.65)` | Table 3 (fixed at 65%) |
| `e_ehc` (EHC, fixed) | 0.35 | `fixed(0.35)` | Table 3 (fixed at 35%) |
| `etalcl` IIV CL1 MPA | 48.6% CV | `log(1+0.486^2)=0.21201` | Table 3 |
| `etalvc` IIV V3 MPA | 59.2% CV | `log(1+0.592^2)=0.30048` | Table 3 |
| `etalq` IIV CL2 MPA | 42.9% CV | `log(1+0.429^2)=0.16895` | Table 3 |
| `etalvp` IIV V4 MPA | 60.0% CV | `log(1+0.600^2)=0.30748` | Table 3 |
| `etalcl_mpag` IIV CLM MPAG | 55.9% CV | `log(1+0.559^2)=0.27193` | Table 3 |
| `propSd` Residual MPA | 41.2% CV | `0.412` | Table 3 |
| `propSd_mpag` Residual MPAG | 45.4% CV | `0.454` | Table 3 |
| Transit-compartment input (Savic) | n/a | `transit(nn, mtt, 1.0)` | Eq. for Ktr; Figure 3 |
| MPA -\> MPAG metabolic conversion | n/a | `fm * kel * central` | Figure 3 caption |
| MPAG biliary uptake to gallbladder | n/a | `fmpag * kel_mpag * central_mpag` | Figure 3 caption |
| Gallbladder emptying windows (1-2 h, 4-6 h) | n/a | `(tpost >= 1) * (tpost <= 2) + (tpost >= 4) * (tpost <= 6)` | Results “Gallbladder emptying was modelled to simulate two release times, 1 and 4 h post dose … ‘turned off’ at 2 h and 6 h post dose” |
| EHC reabsorption fraction | 35% | `ehc * empty_rate * gallbladder` into depot | Figure 3 caption “EHC = biliary recirculation of MPAG into gut was fixed at 35%” |

## Virtual cohort

Original observed data are not publicly available. The cohort below
simulates a typical-value steady-state series at the cohort-mean MMF
dose (1000 mg BID, corresponding to 739 mg MPA-equivalent BID by the MW
ratio MPA/MMF = 320.3 / 433.5), plus a small stochastic cohort (n = 20)
for variability bands. Dosing is q12h to match the hardcoded
`tau = 12 h` gallbladder-emptying periodicity in the model.

``` r

set.seed(20120501L) # paper Accepted Article online date 2011-11-07; rounded forward

mmf_mg       <- 1000              # 1 g MMF BID (typical adult-equivalent of 600 mg/m^2 in this cohort)
mw_ratio     <- 320.3 / 433.5     # MPA MW / MMF MW
dose_mpa_eq  <- mmf_mg * mw_ratio # = 738.9 mg MPA-equivalent per dose
n_doses      <- 24L               # 12 days BID dosing -> steady state
tau_h        <- 12

obs_times <- sort(unique(c(
  seq(0, n_doses * tau_h - tau_h, by = 4),       # coarse pre-steady-state grid
  seq((n_doses - 1) * tau_h, n_doses * tau_h, by = 0.25)  # dense last interval
)))

n_typical   <- 1L
n_stoch     <- 20L

# Single typical-value subject (id 1)
typical_df <- data.frame(id = 1L)

# n_stoch stochastic subjects (id 2..21)
stoch_df <- data.frame(id = 1L + seq_len(n_stoch))

build_events <- function(id_df) {
  id_df |>
    dplyr::rowwise() |>
    do({
      row <- .
      et_obj <- rxode2::et() |>
        rxode2::et(amt = dose_mpa_eq, time = 0, ii = tau_h,
                   addl = n_doses - 1L, cmt = "depot") |>
        rxode2::et(obs_times, cmt = "Cc")
      et_obj$id <- row$id
      et_obj
    }) |>
    dplyr::ungroup() |>
    as.data.frame()
}

events_typ   <- build_events(typical_df)
events_stoch <- build_events(stoch_df)

stopifnot(!anyDuplicated(unique(events_typ[, c("id", "time", "evid")])))
stopifnot(!anyDuplicated(unique(events_stoch[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod         <- readModelDb("Sherwin_2012_mycophenolic_acid")
mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_typ <- rxode2::rxSolve(mod_typical, events = events_typ,
                           returnType = "data.frame", addDosing = FALSE,
                           atol = 1e-8, rtol = 1e-6)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalcl_mpag'
# rxSolve drops the id column for single-subject runs; restore it so
# downstream PKNCA group keys (id, cohort) remain valid.
if (!"id" %in% names(sim_typ)) sim_typ$id <- 1L

sim_stoch <- rxode2::rxSolve(mod, events = events_stoch,
                             returnType = "data.frame", addDosing = FALSE,
                             atol = 1e-8, rtol = 1e-6)
#> ℹ parameter labels from comments will be replaced by 'label()'
if (!"id" %in% names(sim_stoch)) sim_stoch$id <- 1L
```

## Replicate published Figure 1 and Figure 2 (concentration profile shape)

Sherwin 2012 Figure 1 (MPA) and Figure 2 (MPAG) show observed mean +/-
SD plasma concentrations vs. time across 0-9 h post-dose. The
characteristic features that the model needed to reproduce are:

- An MPA sharp initial peak around 0.5-1 h post-dose (absorption);
- A smaller MPA secondary peak at 4-9 h post-dose driven by
  enterohepatic recirculation (the meal-time bile release);
- An MPAG peak delayed by ~1.5-2 h relative to MPA (formation-rate
  limited);
- MPAG concentrations 5-10x higher than MPA across the dosing interval.

The figure below shows the simulated typical-value MPA and MPAG profiles
over the final (steady-state) dosing interval. The vertical dashed bands
mark the two gallbladder-emptying windows. Compare the qualitative shape
against Sherwin 2012 Figures 1 and 2.

``` r

last_typ <- sim_typ |>
  dplyr::filter(time >= (n_doses - 1L) * tau_h & time <= n_doses * tau_h) |>
  dplyr::mutate(tpost = time - (n_doses - 1L) * tau_h) |>
  tidyr::pivot_longer(c(Cc, Cc_mpag),
                      names_to = "analyte", values_to = "conc_mg_L") |>
  dplyr::mutate(analyte = factor(analyte,
    levels = c("Cc", "Cc_mpag"),
    labels = c("MPA (mg/L)", "MPAG (mg/L)")))

meal_windows <- data.frame(
  xmin = c(1, 4), xmax = c(2, 6),
  label = c("Bile window 1", "Bile window 2"))

ggplot(last_typ, aes(tpost, conc_mg_L)) +
  geom_rect(data = meal_windows,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.15, fill = "tan") +
  geom_line(linewidth = 1, colour = "steelblue") +
  facet_wrap(~analyte, scales = "free_y", ncol = 1) +
  labs(x = "Time post-dose (h, final BID interval)",
       y = "Concentration (mg/L)",
       title = "Steady-state typical-value profile: MPA and MPAG",
       caption = "Replicates qualitative shape of Sherwin 2012 Figures 1 (MPA) and 2 (MPAG). Tan bands: bile-release windows (1-2 h, 4-6 h post-dose).")
```

![Typical-value steady-state plasma MPA (top) and MPAG (bottom)
concentration profile over the final BID dosing interval (1000 mg MMF
BID = 739 mg MPA-equivalent BID). Shaded bands mark the bile-release
windows (1-2 h and 4-6 h post-dose). The model reproduces the
qualitative features of Sherwin 2012 Figures 1 and 2: a sharp MPA
initial peak, an MPA secondary peak driven by EHC, an MPAG peak delayed
by ~1.5-2 h, and MPAG concentrations multiple-fold higher than
MPA.](Sherwin_2012_mycophenolic_acid_files/figure-html/figure-1-2-profile-1.png)

Typical-value steady-state plasma MPA (top) and MPAG (bottom)
concentration profile over the final BID dosing interval (1000 mg MMF
BID = 739 mg MPA-equivalent BID). Shaded bands mark the bile-release
windows (1-2 h and 4-6 h post-dose). The model reproduces the
qualitative features of Sherwin 2012 Figures 1 and 2: a sharp MPA
initial peak, an MPA secondary peak driven by EHC, an MPAG peak delayed
by ~1.5-2 h, and MPAG concentrations multiple-fold higher than MPA.

## Visual predictive band (stochastic cohort)

A small n=20 stochastic cohort with the published IIV gives a sense of
the inter-individual spread (median and 10th-90th percentile bands).

``` r

sim_stoch |>
  dplyr::filter(time >= (n_doses - 1L) * tau_h & time <= n_doses * tau_h) |>
  dplyr::mutate(tpost = time - (n_doses - 1L) * tau_h) |>
  tidyr::pivot_longer(c(Cc, Cc_mpag),
                      names_to = "analyte", values_to = "conc_mg_L") |>
  dplyr::mutate(analyte = factor(analyte,
    levels = c("Cc", "Cc_mpag"),
    labels = c("MPA (mg/L)", "MPAG (mg/L)"))) |>
  dplyr::group_by(analyte, tpost) |>
  dplyr::summarise(
    Q10 = quantile(conc_mg_L, 0.10, na.rm = TRUE),
    Q50 = quantile(conc_mg_L, 0.50, na.rm = TRUE),
    Q90 = quantile(conc_mg_L, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(tpost, Q50)) +
  geom_rect(data = meal_windows,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, alpha = 0.15, fill = "tan") +
  geom_ribbon(aes(ymin = Q10, ymax = Q90), alpha = 0.3, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) +
  facet_wrap(~analyte, scales = "free_y", ncol = 1) +
  labs(x = "Time post-dose (h, final BID interval)",
       y = "Concentration (mg/L)",
       title = "Stochastic VPC band - final BID interval",
       caption = "Median (line) + 10th-90th percentile (ribbon), n = 20.")
```

![Visual predictive band (median, 10th-90th percentile) of MPA and MPAG
plasma concentrations over the final BID dosing interval at 1000 mg MMF
BID. Stochastic cohort n = 20 with the published IIV %CV per
parameter.](Sherwin_2012_mycophenolic_acid_files/figure-html/vpc-band-1.png)

Visual predictive band (median, 10th-90th percentile) of MPA and MPAG
plasma concentrations over the final BID dosing interval at 1000 mg MMF
BID. Stochastic cohort n = 20 with the published IIV %CV per parameter.

## PKNCA validation

The source paper does not tabulate NCA-derived `Cmax` / `Tmax` / `AUC` /
`half-life` values. The model itself reports CL1/F = 25.3 L/h, which
constrains the expected steady-state AUC by mass-balance: for 1000 mg
MMF (739 mg MPA-equivalent) BID,
`AUC_ss = Dose / (CL/F * tau_dose_interval / tau)`. Per the source
paper, dose-interval AUC = Dose / CL = 739 / 25.3 = 29.2 mg\*h/L.

``` r

sim_nca <- sim_typ |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::filter(time >= (n_doses - 1L) * tau_h & time <= n_doses * tau_h) |>
  dplyr::mutate(tpost = time - (n_doses - 1L) * tau_h,
                cohort = "Sherwin_2012_MPA")

# Guarantee a time=0 row per (id, cohort) so PKNCA can anchor AUC0-*.
sim_nca <- dplyr::bind_rows(
  sim_nca,
  sim_nca |> dplyr::distinct(id, cohort) |>
    dplyr::mutate(tpost = 0, Cc = 0, Cc_mpag = 0)
) |>
  dplyr::distinct(id, cohort, tpost, .keep_all = TRUE) |>
  dplyr::arrange(id, cohort, tpost)

# MPA NCA on the last steady-state interval
nca_one <- function(analyte_col, conc_label) {
  conc_df <- sim_nca |>
    dplyr::select(id, tpost, cohort,
                  conc = !!rlang::sym(analyte_col)) |>
    dplyr::filter(!is.na(conc))
  conc_obj <- PKNCA::PKNCAconc(conc_df, conc ~ tpost | cohort + id)

  dose_df <- sim_nca |>
    dplyr::distinct(id, cohort) |>
    dplyr::mutate(tpost = 0, amt = dose_mpa_eq)
  dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ tpost | cohort + id)

  intervals <- data.frame(start = 0, end = tau_h,
                          cmax = TRUE, tmax = TRUE, auclast = TRUE)
  data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
  res  <- PKNCA::pk.nca(data)
  out  <- as.data.frame(res$result)
  out$analyte <- conc_label
  out
}

nca_all <- bind_rows(
  nca_one("Cc",      "MPA"),
  nca_one("Cc_mpag", "MPAG")
)
```

``` r

nca_summary <- nca_all |>
  dplyr::group_by(analyte, PPTESTCD) |>
  dplyr::summarise(value = mean(PPORRES, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = value) |>
  dplyr::select(analyte, dplyr::any_of(c("cmax", "tmax", "auclast")))

knitr::kable(
  nca_summary,
  caption = paste(
    "Simulated typical-value steady-state NCA over the final BID dosing",
    "interval (0-12 h post-dose, 739 mg MPA-equivalent dose). Cmax in",
    "mg/L; tmax in h; AUClast in mg*h/L."
  ),
  digits = 3
)
```

| analyte |   cmax | tmax | auclast |
|:--------|-------:|-----:|--------:|
| MPA     | 11.052 | 1.75 |  36.172 |
| MPAG    | 15.147 | 3.00 | 108.996 |

Simulated typical-value steady-state NCA over the final BID dosing
interval (0-12 h post-dose, 739 mg MPA-equivalent dose). Cmax in mg/L;
tmax in h; AUClast in mg\*h/L. {.table}

### EHC contribution to MPA AUC

For an EHC model, `AUC_ss(0-tau)` is NOT equal to `Dose / CL1` because
the recycled MPA from the gallbladder returns mass to the central
compartment and inflates exposure. The source paper’s introduction
states “enterohepatic recirculation … can contribute to an increase in
exposure to MPA of 40% (range 10-60%).” The simulated EHC contribution
is therefore computed as `(AUC_sim - Dose / CL1) / (Dose / CL1)`:

``` r

mpa_auc_sim <- nca_summary |>
  dplyr::filter(analyte == "MPA") |>
  dplyr::pull(auclast)

paper_cl1 <- 25.3  # L/h, Sherwin 2012 Table 3
no_ehc_auc <- dose_mpa_eq / paper_cl1
ehc_pct    <- 100 * (mpa_auc_sim - no_ehc_auc) / no_ehc_auc

cmp <- data.frame(
  Quantity = c("MPA AUC0-12 simulated (mg*h/L)",
               "MPA AUC0-12 if no EHC = Dose/CL1 (mg*h/L)",
               "EHC contribution to MPA AUC (%)"),
  Value    = c(round(mpa_auc_sim, 2),
               round(no_ehc_auc, 2),
               round(ehc_pct, 1)),
  check.names = FALSE
)

knitr::kable(
  cmp,
  caption = paste(
    "Simulated MPA AUC0-12 vs. the AUC that would arise in the absence",
    "of EHC (= Dose / CL1). The difference quantifies the model's EHC",
    "contribution to MPA exposure; the source paper's introduction",
    "reports a typical 40% EHC contribution to MPA AUC (range 10-60%)."
  ),
  align = c("l", "r")
)
```

| Quantity                                   | Value |
|:-------------------------------------------|------:|
| MPA AUC0-12 simulated (mg\*h/L)            | 36.17 |
| MPA AUC0-12 if no EHC = Dose/CL1 (mg\*h/L) | 29.20 |
| EHC contribution to MPA AUC (%)            | 23.90 |

Simulated MPA AUC0-12 vs. the AUC that would arise in the absence of EHC
(= Dose / CL1). The difference quantifies the model’s EHC contribution
to MPA exposure; the source paper’s introduction reports a typical 40%
EHC contribution to MPA AUC (range 10-60%). {.table}

The simulated EHC contribution falls within the source paper’s reported
10-60% range for the typical-value parameter set. The full set of NCA
values (`Cmax`, `Tmax`, `AUClast`) for both analytes is shown in the
table above, alongside the visual replication of Figures 1 and 2.

## Assumptions and deviations

- **Residual error: combined vs proportional.** The source paper’s
  Methods section describes residual variability as a combined
  additive + proportional error model (Eq. 2), but Table 3 reports only
  a single `error CV%` per analyte (MPA 41.2%, MPAG 45.4%) with no
  separate additive component. The packaged model encodes the values as
  proportional residual SDs (`propSd = 0.412`, `propSd_mpag = 0.454`)
  and omits the additive component. If the source paper’s NONMEM control
  stream had `EPS(1)` proportional and `EPS(2)` additive but `EPS(2)`
  was fixed to zero or numerically negligible, the proportional-only
  encoding is exact; otherwise the packaged model carries slightly less
  residual SD at low concentrations than the original fit.

- **No additional MPA renal elimination (k_30 = 0).** Sherwin 2012
  Figure 3 caption distinguishes `k_30 = renal eliminated MPA` from
  `(FM * k_35) = fraction of MPA metabolized to MPAG`, but Table 3
  reports only a single apparent MPA clearance `CL1/F = 25.3 L/h` with
  no separate `k_30 * V3` value. The packaged model interprets `CL1/F`
  as the total apparent MPA clearance and partitions 0.85 of the
  eliminated flux to MPAG (via `fm * kel`) with the remaining 0.15 lost
  (interpretable as the AcMPAG + small renal MPA bucket). MPA renal
  clearance in adult humans is

  ~ 1% of total CL, so the `k_30 = 0` simplification is physiologically
  near-exact for the modeled flux of MPA into MPAG.

- **FMPAG = 0.65 is the biliary fraction.** Figure 3 caption describes
  `(FMPAG * k_56) = ... 65% excretion of transferred MPAG from central compartment to gall bladder, with 35% accumulation in the gallbladder`
  and `(1 - FMPAG) * k_50 = fraction of renal eliminated MPAG`, which
  together imply FMPAG = 0.65 is the biliary fraction and CLM/F = 2.5
  L/h is the apparent renal CL of MPAG (=
  `(1 - FMPAG) * total MPAG elimination`). The body of the paper also
  contains the sentence “It was assumed that approximately 65% of MPAG
  is excreted unchanged by renal elimination” which appears to
  contradict the figure caption; the packaged model follows the figure
  caption (and Table 3’s `FMPAG = 65%` label) because this is the
  parameterisation consistent with the named partitioning in Figure 3
  and with prominent EHC.

- **EHC = 0.35 applied at the gallbladder-to-gut step.** The paper text
  states “EHC = biliary recirculation of MPAG into gut was fixed at 35%”
  with the explicit destination “into gut”, which the packaged model
  encodes as `ehc * empty_rate * gallbladder` being added to the gut
  (depot) ODE during the meal-time windows; the complement
  `(1 - ehc) * empty_rate * gallbladder` is excreted to feces and not
  tracked. An alternative reading would apply EHC = 0.35 at the
  bile-to-gallbladder uptake step instead (with 100% of the accumulated
  gallbladder content reaching the gut on emptying), which is
  mathematically equivalent for the steady-state ratio of EHC
  contribution to systemic MPA exposure but redistributes the mass in
  transit. Either reading reproduces the qualitative secondary-peak
  behaviour seen in the source figures.

- **Gallbladder emptying rate constant `k_empty`.** The source paper
  specifies only the bile-release time windows (1-2 h and 4-6 h
  post-dose) without an explicit emptying rate constant. The packaged
  model uses `k_empty = 5 1/h` during the open windows, giving
  near-complete first-order emptying within each window. Slower rates
  would push more bile through to the second window without changing the
  steady-state EHC contribution.

- **Hardcoded BID dosing interval `tau = 12 h`.** The
  gallbladder-emptying windows are computed from
  `tpost = t - floor(t / 12) * 12`, baking the BID periodicity into the
  model. Users wanting a different dose interval should edit the `tau`
  constant directly in
  [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html).
  Single-dose simulation is still valid because the modulo arithmetic
  evaluates correctly for the first interval.

- **MMF -\> MPA mass conversion is external to the model.** The dose
  unit declared in `units` is MPA-mass-equivalent mg. Users passing an
  MMF dose in mg should multiply by 0.739 (= 320.3 / 433.5, the MW ratio
  MPA/MMF) before placing the dose on the `depot` compartment. The
  vignette demonstrates this conversion explicitly
  (`dose_mpa_eq <- mmf_mg * mw_ratio`).

- **`V_MPAG = V3 MPA` is a structural assumption.** Sherwin 2012 Table 3
  lists `V_M = V3 MPA` as a fixed structural relationship (the MPAG
  central volume of distribution equals the MPA central volume). The
  packaged model implements this by writing `v_mpag <- vc` inside
  [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html),
  which means the IIV `etalvc` applied to V3 MPA also applies
  (identically) to V_MPAG.

- **No covariates retained.** Bodyweight, age, sex, race, ethnicity, and
  disease duration were all screened in the covariate analysis (Sherwin
  2012 Results, Covariate analysis). None were retained: bodyweight
  produced significant improvement in a simpler two-compartment
  first-order-absorption submodel but failed when added to the full EHC
  model (attributed to over-parameterisation), and the others were not
  significant. The packaged model therefore has an empty `covariateData`
  list and documents the screened-but-rejected covariates in
  `covariatesDataExcluded` so the provenance of the covariate screen is
  preserved.

- **Baseline-only chemistries (serum albumin, AST, ALT, serum
  creatinine, urine protein:creatinine) not included.** Sherwin 2012
  acknowledges in the Discussion that other MPA popPK studies (Sam & Joy
  2008 in adult glomerulonephritis) found serum albumin and estimated
  creatinine clearance to affect MPA clearance, but these chemistries
  were collected at the screening visit only – NOT at the PK visit – so
  the authors judged it inappropriate to include them as time-of-PK
  covariates. The packaged model preserves this restriction.

- **No co-medication effects.** 18 of 19 subjects were on prednisone
  (mean 17.2 mg/day), 3 on high-dose i.v. methylprednisolone, 17 on
  hydroxychloroquine, 11 on NSAIDs. The paper notes that corticosteroids
  may induce UGT-glucuronidation of MPA, but did not have sufficient
  variation in the cohort to estimate the effect; corticosteroid
  co-medication is therefore treated as a uniform background condition.

- **Standard errors of parameter estimates.** Sherwin 2012 reports
  parameter estimates with 95% confidence intervals derived from 200
  non-parametric bootstrap runs (Table 3, 96.5% bootstrap success rate).
  The packaged model carries only the point estimates from the final
  model column of Table 3.
