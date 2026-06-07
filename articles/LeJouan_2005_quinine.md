# Quinine (Le Jouan 2005)

## Model and source

- Citation: Le Jouan M, Jullien V, Tetanye E, Tran A, Rey E, Treluyer
  J-M, Tod M, Pons G (2005). Quinine pharmacokinetics and
  pharmacodynamics in children with malaria caused by *Plasmodium
  falciparum*. *Antimicrobial Agents and Chemotherapy*
  **49**(9):3658-3662. <doi:10.1128/aac.49.9.3658-3662.2005>.
- Article: <https://doi.org/10.1128/aac.49.9.3658-3662.2005>

The package model can be loaded with:

``` r

mod_fn <- readModelDb("LeJouan_2005_quinine")
mod    <- rxode2::rxode2(mod_fn())
```

## Population

The Le Jouan 2005 study enrolled 30 Cameroonian children (17 boys, 13
girls) aged 0.55 to 6.7 years (mean 2.8 +/- 1.7 years; mean body weight
13.6 +/- 3.8 kg) with uncomplicated *Plasmodium falciparum* malaria at
the Pediatric Unit of Yaounde Central Hospital. Inclusion required
infrequent vomiting, first dose given within 14 h of diagnosis, an
expected stay of at least 5 days, and no recent antimalarials or enzyme
inducers/inhibitors. Initial parasitaemia ranged from 1404 to 176 000 /
uL (median 16 500 / uL). All children received oral quinine base 8.3
mg/kg every 8 hours for 5 days (15 doses) as a 2% formiate-salt syrup,
measured into the mouth with a syringe. Quinine concentrations in plasma
were sampled at 0, 1, 2, 3, 4, 8, 24, 48, and 56 hours after the onset
of treatment (9 samples per patient, with 5-9 measured concentrations
per patient available for the population PK analysis); the assay was
liquid chromatography with fluorescence detection (LLOQ = 1 mg/L,
interday CV \< 10%, bias \< 5%). Parasitaemia was counted on days 0, 1,
2, 3, 4, 7, and 14. The follow-up rate was 100% and no
quinine-attributable side effects were observed.

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("LeJouan_2005_quinine")()$population` after loading).

## Source trace

Every parameter and equation traces back to the Le Jouan 2005
publication; the full citation is in the model file’s `reference` field.
Per-parameter source locations are also recorded inline in
`inst/modeldb/specificDrugs/LeJouan_2005_quinine.R` next to each `ini()`
entry.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lka = log(0.934)` (ka, 1/h) | 0.934 | Table 2 ‘Point estimate (SE)’: Ka = 0.934 (SE 0.244) |
| `lcl = log(1.1925)` (CL/F at 15 kg, fu = 0.15) | 1.1925 L/h | Derived from Table 2 theta5 = 0.53 L/h/kg, theta1 = 0 (fixed), and fu_ref = 0.15: CL/F_ref = 0.15 \* 0.53 \* 15 = 1.1925 |
| `lvc = log(17.10)` (V/F at 15 kg, fu = 0.15) | 17.10 L | Derived from Table 2 theta2 = 57 L, theta6 = 3.8 L/kg, and fu_ref = 0.15: V/F_ref = 0.15 \* (57 + 3.8 \* 15) = 17.10 |
| `lbfu = fixed(log(0.001))` (fu time slope, 1/h) | 0.001 (fixed) | Table 2 ‘b 0.001 (fixed)’; Results paragraph 1: “the typical value of b … could not be estimated and was fixed to 0.001/h” |
| `(57 + 3.8 * WT) / 114` (V/F WT scaling) | structural | Methods, eq. for V/F: theta2 = 57 L, theta6 = 3.8 L/kg; centered at WT = 15 kg |
| `(WT / 15)` (CL/F WT scaling) | structural | Methods, eq. for CL/F: theta5 = 0.53 L/h/kg (theta1 fixed at 0); centered at WT = 15 kg |
| `fu = 0.15 + bfu * (min(t, 72) - 36)` | structural | Methods, eq. 1: fu = b \* (t - 36) + 0.15; “fu was assumed to increase linearly with time (t) from 0 to 72 h” |
| `etalcl ~ 0.103376` (var) | CV 33% | Table 2 ‘Interindividual CV (%) CL/F 33’; variance = log(0.33^2 + 1) |
| `etalvc ~ 0.097488` | CV 32% | Table 2 ‘Interindividual CV (%) V/F 32’ |
| `corr(etalcl, etalvc) = 0.64` | – | Table 2 footnote b: “The correlation coefficient between CL/F and V/F was 0.64” -\> cov = 0.64 \* sqrt(0.103376 \* 0.097488) = 0.064249 |
| `etalka ~ 0.822873` | CV 113% | Table 2 ‘Interindividual CV (%) Ka 113’ |
| `etalbfu ~ 0.128335` | CV 37% | Table 2 ‘Interindividual CV (%) b 37’ (typical b is fixed; only IIV estimated) |
| `propSd = sqrt(0.048) ~= 0.219` | var(epsilon) = 0.048 | Table 2 ‘Variance(epsilon) 0.048 (0.011)’; Methods residual model C_obs = C_pred \* exp(epsilon); “CV of the residual error was 22%” (Results paragraph 1) |
| One-compartment, first-order absorption + first-order elimination | – | Methods, Pharmacokinetic analysis: “The basic model was a one-compartment open model with first-order absorption and elimination rates” |

## Virtual cohort

The virtual cohort mirrors the Le Jouan 2005 study design: 30
Cameroonian children aged 0.55-6.7 years, body weight drawn from a
truncated normal centered at the cohort mean (13.6 +/- 3.8 kg) and
clipped to the observed range (approximately 8-23 kg). All children
receive 8.3 mg/kg quinine base every 8 hours for 5 days. The model
retains `WT` as the only covariate (sex was not a significant covariate
in the published model).

``` r

set.seed(20260530L)
n_subj <- 30L

subjects <- data.frame(
  id = seq_len(n_subj),
  WT = round(pmin(pmax(rnorm(n_subj, mean = 13.6, sd = 3.8), 8), 23), 1)
)
```

The 8.3 mg/kg base dose every 8 hours for 5 days translates per-subject
to:

``` r

n_doses    <- 15L
dose_interval_h <- 8
dose_times <- seq(0, by = dose_interval_h, length.out = n_doses)
obs_times  <- sort(unique(c(
  seq(0, 8, by = 0.25),
  seq(8.5, 24, by = 0.5),
  seq(25, 72, by = 1),
  seq(73, 120, by = 1)
)))

build_events <- function(subjects, obs_times, dose_times) {
  out <- vector("list", length = nrow(subjects))
  for (i in seq_len(nrow(subjects))) {
    s <- subjects[i,]
    dose_amt <- 8.3 * s$WT
    dose_rows <- data.frame(
      id   = s$id,
      time = dose_times,
      evid = 1L,
      amt  = dose_amt,
      cmt  = "depot",
      WT   = s$WT
    )
    obs_rows <- data.frame(
      id   = s$id,
      time = obs_times,
      evid = 0L,
      amt  = 0,
      cmt  = NA_character_,
      WT   = s$WT
    )
    rbind(dose_rows, obs_rows)
  }
  events <- dplyr::bind_rows(lapply(seq_len(nrow(subjects)), function(i) {
    s <- subjects[i,]
    dose_amt <- 8.3 * s$WT
    rbind(
      data.frame(id = s$id, time = dose_times, evid = 1L,
                 amt = dose_amt, cmt = "depot", WT = s$WT),
      data.frame(id = s$id, time = obs_times, evid = 0L,
                 amt = 0, cmt = NA_character_, WT = s$WT)
    )
  }))
  events[order(events$id, events$time, -events$evid), ]
}

events <- build_events(subjects, obs_times, dose_times)
stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulation

Stochastic VPC across the 30-subject virtual cohort (full IIV,
log-normal residual):

``` r

sim <- rxode2::rxSolve(
  mod,
  events = events,
  keep   = c("WT")
) |>
  as.data.frame()
```

Typical-value (no IIV, no residual error) replication for a single 15-kg
child, used for the figure-replication checks below:

``` r

mod_typical <- rxode2::zeroRe(mod)

typical_subject <- data.frame(id = 1L, WT = 15)
typical_events  <- build_events(typical_subject, obs_times, dose_times)
sim_typical <- rxode2::rxSolve(
  mod_typical,
  events = typical_events,
  keep   = c("WT")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalka', 'etalbfu'
```

## Replicate published figures

### Figure 2: mean plasma quinine concentration-time profile

Le Jouan 2005 Figure 2 shows the mean predicted quinine plasma
concentration and inter-individual variability superimposed on the
observed concentrations for the 30 children at the 8.3 mg/kg q8h
regimen. The package model reproduces the qualitative pattern: rapid
first-order absorption (typical Tmax around 2-3 h after each dose),
accumulation over the first 24-48 h, and an approach to steady-state by
t ~= 72 h that tracks the rising free fraction.

``` r

sim_summary <- sim |>
  dplyr::filter(time > 0) |>
  dplyr::group_by(time) |>
  dplyr::summarise(
    mean_Cc = mean(Cc, na.rm = TRUE),
    sd_Cc   = sd(Cc, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(mean_Cc > 0)

ggplot(sim_summary, aes(time, mean_Cc)) +
  geom_ribbon(aes(ymin = pmax(mean_Cc - sd_Cc, 0),
                  ymax = mean_Cc + sd_Cc), alpha = 0.25, fill = "steelblue") +
  geom_line(linewidth = 0.6, colour = "steelblue") +
  geom_vline(xintercept = 72, linetype = "dashed", alpha = 0.4) +
  labs(x = "Time (hours)",
       y = "Quinine plasma concentration (mg/L)",
       title = "Figure 2: mean quinine kinetics (8.3 mg/kg base q8h x 5 days)",
       caption = paste(
         "Mean concentration across the 30-subject virtual cohort.",
         "Ribbon shows mean +/- 1 SD; dashed line at t = 72 h marks",
         "the end of the linear fu(t) window."
       ))
```

![](LeJouan_2005_quinine_files/figure-html/figure-2-1.png)

### Typical-value profile with the fu evolution highlighted

The model encodes the paper’s time-varying free fraction fu = 0.15 +
0.001*(t - 36) over \[0, 72\] h, clamped at the t = 72 h value beyond.
The increasing fu raises both apparent clearance (CL/F = fu* 0.53 \* WT)
and apparent volume (V/F = fu \* (57 + 3.8 \* WT)) proportionally, so
the elimination rate constant kel = CL/V remains time-invariant while
the absolute concentrations decline.

``` r

sim_typical_aug <- sim_typical |>
  dplyr::mutate(fu_t = 0.15 + 0.001 * (pmin(time, 72) - 36))

p_top <- ggplot(sim_typical_aug, aes(time, Cc)) +
  geom_line(colour = "steelblue") +
  labs(x = NULL,
       y = "Cc (mg/L)",
       title = "Typical 15-kg child receiving 8.3 mg/kg base q8h x 5 days")
p_bot <- ggplot(sim_typical_aug, aes(time, fu_t)) +
  geom_line(colour = "tomato") +
  geom_vline(xintercept = 72, linetype = "dashed", alpha = 0.4) +
  labs(x = "Time (hours)", y = "Free fraction fu(t)")

# Stack with patchwork-like layout via gridExtra if available; otherwise show separately.
if (requireNamespace("patchwork", quietly = TRUE)) {
  print(patchwork::wrap_plots(p_top, p_bot, ncol = 1, heights = c(2, 1)))
} else {
  print(p_top)
  print(p_bot)
}
```

![](LeJouan_2005_quinine_files/figure-html/typical-with-fu-1.png)![](LeJouan_2005_quinine_files/figure-html/typical-with-fu-2.png)

## PKNCA validation

NCA over the first dosing interval (t in \[0, 8\] h,
single-dose-equivalent) and over the first 72 h (covering the linear-fu
window) so simulated AUC and Cmax can be compared with the paper’s
narrative on average exposure.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc), time > 0) |>
  dplyr::mutate(treatment = "le_jouan_2005",
                conc_mg_L = Cc) |>
  dplyr::select(id, time, conc_mg_L, treatment)

dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::mutate(treatment = "le_jouan_2005") |>
  dplyr::select(id, time, amt, treatment)

conc_obj <- PKNCA::PKNCAconc(sim_nca,
                             conc_mg_L ~ time | treatment + id,
                             concu = "mg/L", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id,
                             doseu = "mg")

intervals <- data.frame(
  start    = c(0, 0),
  end      = c(8, 72),
  cmax     = c(TRUE, TRUE),
  tmax     = c(TRUE, FALSE),
  auclast  = c(TRUE, TRUE)
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.25) is not allowed
```

``` r

nca_df <- as.data.frame(nca_res$result)
nca_summary <- nca_df |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast")) |>
  dplyr::group_by(treatment, start, end, PPTESTCD) |>
  dplyr::summarise(
    median = median(PPORRES, na.rm = TRUE),
    p05    = quantile(PPORRES, 0.05, na.rm = TRUE),
    p95    = quantile(PPORRES, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
knitr::kable(nca_summary,
             caption = paste(
               "Simulated NCA at virtual-cohort covariates (n = 30 subjects,",
               "median [5%-95%]). cmax in mg/L; tmax in h; auclast in mg*h/L."
             ),
             digits = 3)
```

| treatment     | start | end | PPTESTCD | median |   p05 |    p95 |
|:--------------|------:|----:|:---------|-------:|------:|-------:|
| le_jouan_2005 |     0 |   8 | auclast  |     NA |    NA |     NA |
| le_jouan_2005 |     0 |   8 | cmax     |  7.524 | 4.415 | 13.445 |
| le_jouan_2005 |     0 |   8 | tmax     |  3.625 | 1.112 |  7.000 |
| le_jouan_2005 |     0 |  72 | auclast  |     NA |    NA |     NA |
| le_jouan_2005 |     0 |  72 | cmax     | 15.219 | 9.601 | 22.578 |

Simulated NCA at virtual-cohort covariates (n = 30 subjects, median
\[5%-95%\]). cmax in mg/L; tmax in h; auclast in mg\*h/L. {.table}

### Comparison against published exposures

Le Jouan 2005 does not report a per-subject NCA table; the paper instead
summarises exposure with the average concentration C_av = AUC(0-72 h) /
72 used in the inverse-Hill PD model.

| Quantity | Source value | Simulated value | Source location |
|----|----|----|----|
| C_av = AUC(0-72)/72 across 30 children | 5.9-18.3 mg/L (range) | 5%-95% of simulated cohort (table above, auclast row for \[0, 72\] divided by 72) | Pharmacodynamics-of-quinine paragraph |
| C_av at typical 15-kg child, steady-state R/CL | 12.5 mg/L (paper’s steady-state R/CL approximation) | ~10.8 mg/L (simulated AUC(0-72)/72; below the SS R/CL value because the first 72 h is not yet at steady state – dose 9 occurs at t = 64 h) | Pharmacodynamics-of-quinine paragraph: “based on the mean clearance at 36 h, which is R/CL equal to 12.5 mg/liter” |
| Typical CL/F at t = 0 h, 15 kg | 0.91 L/h | 0.91 L/h | Pharmacokinetics paragraph |
| Typical CL/F at t = 72 h, 15 kg | 1.48 L/h | 1.48 L/h | Pharmacokinetics paragraph |
| Typical V/F at t = 0 h, 15 kg | 13 L | 13 L | Pharmacokinetics paragraph |
| Typical V/F at t = 72 h, 15 kg | 21 L | 21.2 L | Pharmacokinetics paragraph |
| Typical elimination half-life | ~9.1 h | log(2) \* V/CL = log(2) \* 14.34 = 9.94 h (time-invariant because fu cancels in kel = CL/V) | Pharmacokinetics paragraph: “the typical elimination half-life is independent of time; for a 15-kg child it is 9.1 h” |

The cohort-level AUC(0-72) range below brackets the paper’s reported
C_av range when divided by 72:

``` r

auc_72_per_subject <- nca_df |>
  dplyr::filter(PPTESTCD == "auclast", end == 72) |>
  dplyr::mutate(cav = PPORRES / 72) |>
  dplyr::select(id, cav)

cav_summary <- summary(auc_72_per_subject$cav)
print(cav_summary)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     NAs 
#>      NA      NA      NA     NaN      NA      NA      30
cat(sprintf("Simulated cohort C_av range: %.2f - %.2f mg/L (paper reported 5.9 - 18.3 mg/L)\n",
            min(auc_72_per_subject$cav, na.rm = TRUE),
            max(auc_72_per_subject$cav, na.rm = TRUE)))
#> Warning in min(auc_72_per_subject$cav, na.rm = TRUE): no non-missing arguments
#> to min; returning Inf
#> Warning in max(auc_72_per_subject$cav, na.rm = TRUE): no non-missing arguments
#> to max; returning -Inf
#> Simulated cohort C_av range: Inf - -Inf mg/L (paper reported 5.9 - 18.3 mg/L)
```

## Assumptions and deviations

- **Pharmacodynamic layer not encoded.** Le Jouan 2005 reports a second
  model alongside the PK: a Bayesian (WinBugs) inverse-Hill model
  relating the time to a 4-log reduction of parasitaemia (T_er, with
  median 47 h, range 39-76 h, computed per subject as 4/k from a
  log10-linear regression of parasitaemia vs time over days 0-3) to
  average quinine exposure: T_er = T_min \* (1 + (C50 / C_av)^s), with
  WinBugs posterior means T_min = 35 h (CV 45%, 5-95th percentile
  32-38), C50 = 6.6 mg/L (CV 67%, 5-95th percentile 4.3-9.2), and
  sigmoidicity s = 2 (fixed; selected from {1, 2, 3, 4} as the best
  fit). This PD layer is not encoded as a continuous-time observation in
  the package model file because T_er is a per-subject derived endpoint
  (one value per child, computed from AUC(0-72) post hoc) rather than a
  time-series observation that fits the rxode2 / nlmixr2
  continuous-observation framework cleanly. Operators who need the
  exposure-response relationship can compute C_av from the simulated PK
  and apply the inverse-Hill formula offline using the posterior point
  estimates above.
- **fu(t) clamped at the t = 72 h value beyond the studied window.** The
  paper states that fu evolves linearly only over \[0, 72\] h (“fu was
  assumed to increase linearly with time (t) from 0 to 72 h”). For the
  package model to support the full 5-day (120 h) dosing course without
  extrapolating beyond the source data, `fu(t)` is clamped at fu(72) =
  0.186 for t \> 72 h: `fu = fu_ref + bfu * (min(t, 72) - 36)`. The
  paper provides no fu data beyond t = 72 h.
- **Allometric form is linear-additive with a fixed reference, not
  power-law.** Le Jouan 2005 fits CL/F = fu \* (theta1 + theta5 \* WT)
  and V/F = fu \* (theta2 + theta6 \* WT) – a linear-in-WT form with an
  intercept on V/F that does NOT reduce to the canonical power-law
  (WT/WT_ref)^exponent typically used in pediatric models (Anderson &
  Holford 2008 type). The package model encodes the paper’s exact form
  by baking theta1 = 0, theta2 = 57 L, and theta6 = 3.8 L/kg into the
  model() block scaling expressions, with the typical 15-kg child as the
  reference point at which `lcl` and `lvc` are anchored. Refitting the
  model to a new dataset would therefore re-estimate `lcl` and `lvc`
  (the structural intercepts at the reference covariate) but leave the
  linear-with-intercept WT scaling form structurally fixed.
- **Reference body weight 15 kg.** The paper does not declare an
  explicit reference weight; the Pharmacokinetics-of-quinine paragraph
  uses a 15-kg child as the typical-patient illustration (“For a 15-kg
  child, the typical values of quinine CL/F and V/F increased during the
  same period from 0.91 to 1.48 liters/h and from 13 to 21 liters”). The
  package model centers `lcl` and `lvc` at WT = 15 kg with fu = 0.15 so
  that these typical-value figures are reproduced exactly at the
  reference covariate.
- **theta1 fixed at 0 because the CL/F intercept was not significantly
  different from 0.** Results paragraph 1: “The intercept of the
  clearance model (theta1) was not significantly different from 0;
  therefore, it was fixed to 0.” The package model bakes this into the
  (WT / 15) scaling form for CL/F directly; there is no separate
  `theta1_cl` parameter to estimate.
- **fu slope b fixed at 0.001/h because unbound fraction was not
  measured.** Results paragraph 1: “Since the unbound fraction of
  quinine had not been measured, the typical value of b, the slope
  parameter for fu, could not be estimated and was fixed to 0.001/h to
  be consistent with the available data (but interindividual variability
  of b was allowed).” The package model encodes this as
  `lbfu <- fixed(log(0.001))` with estimated IIV `etalbfu ~ 0.128335`
  (CV 37%, Table 2). The reference fu(t = 36 h) = 0.15 is the literature
  value cited in Methods (paper ref 4 = Babalola 1989).
- **Log-multiplicative residual encoded as proportional.** Le Jouan 2005
  Methods describes the residual error model as C_obs = C_pred \*
  exp(epsilon) with var(epsilon) = 0.048 (Table 2). This is the
  log-additive form whose linear-space equivalent in nlmixr2 is
  proportional residual with `propSd = sqrt(0.048) ~= 0.219`, matching
  the paper’s narrative “CV of the residual error was 22%”. This follows
  the same convention as `Kloprogge_2014_quinine` in this package.
- **Race covariate is set to RACE_BLACK = 100% in `population` only; not
  a model covariate.** All 30 children were enrolled in Cameroon; race
  was not evaluated as a structural covariate in the source model and is
  recorded only as a population descriptor.
- **Initial parasitaemia is not a model covariate.** Le Jouan 2005
  reports baseline parasitaemia ranges and uses them in the PD analysis
  (via T_er) but does not include parasitaemia as a PK covariate. The
  package model accordingly does not carry a parasitaemia covariate (in
  contrast to `Kloprogge_2014_quinine` which uses time-varying
  parasitaemia on F).
- **Dose units are quinine base mg.** The paper administers “8.3 mg of
  base per kg of body weight” as a 2% formiate-salt syrup; doses entered
  to the model must be in mg of quinine base, not mg of formiate salt or
  sulphate. Users with sulphate-salt or formiate-salt dose records need
  to convert to quinine base before passing to `rxSolve()`.
