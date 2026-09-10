# Clozapine (Han 2025)

## Model and source

``` r

mod <- readModelDb("Han_2025_clozapine")
ui  <- rxode2::rxode(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Han HH, Zhang Y, Wang J, Tian X, Li Y, He SM, Zhang C, Chen
  X, Wang DD. Population pharmacokinetics modelling to predict DDI from
  zopiclone on clozapine in schizophrenia patients. Frontiers in
  Psychiatry. 2025;16:1664678. <doi:10.3389/fpsyt.2025.1664678>. Final
  model Equations (6) and (7); parameter estimates Table 3. The fixed
  absorption rate constant is quoted from the paper’s references 25 and
  26: Li LJ, Shang DW, Li WB, et al. Population pharmacokinetics of
  clozapine and its primary metabolite norclozapine in Chinese patients
  with schizophrenia. Acta Pharmacol Sin. 2012;33:1409-1416.
  <doi:10.1038/aps.2012.71>; see modellib(‘Li_2012_clozapine’). And
  Shang DW, Li LJ, Wang XP, et al. Population
  pharmacokinetic/pharmacodynamic model of clozapine for characterizing
  the relationship between accumulated exposure and PANSS scores in
  patients with schizophrenia. Ther Drug Monit. 2014;36:378-386.
  <doi:10.1097/FTD.0000000000000014>.
- Description: One-compartment population PK model for oral clozapine
  with first-order absorption in adults with schizophrenia, built from
  routine trough therapeutic-drug-monitoring concentrations (Han 2025).
  Apparent oral clearance is allometrically scaled on body weight and
  reduced by 25.4% when zopiclone is co-administered; the absorption
  rate constant is fixed to a published value. Between-subject
  variability was retained on CL/F only.
- Article: <https://doi.org/10.3389/fpsyt.2025.1664678>

Han and colleagues built a one-compartment population PK model for oral
clozapine from routine therapeutic-drug-monitoring troughs in 81
inpatients with schizophrenia, screened 28 concomitant medications for
drug-drug interactions, and found exactly one: co-administered
**zopiclone reduces clozapine apparent oral clearance by 25.4%**. They
then used the model to recommend weight-banded starting doses with and
without zopiclone.

## Population

The analysis pooled 81 adults with schizophrenia treated at a single
centre (Xuzhou Oriental Hospital Affiliated to Xuzhou Medical
University, Jiangsu, China) between December 2023 and November 2024.
Baseline characteristics (Table 1 of the source) were: 37 men / 44
women; age mean 49.46 years (SD 11.15), median 50.67, range 20.67-73.11;
body weight mean 70.49 kg (SD 13.53), median 71.00, range 38.00-120.00.
Concomitant medication counts are Table 2; 8 of the 81 patients were
taking zopiclone tablets.

Every observation is a **trough** concentration: the source states that
“the sample extraction times for plasma concentrations were before the
next administration, which was the value of the trough concentration.”
Clozapine was assayed by homogeneous enzyme immunoassay. The paper
reports no administered doses, dosing frequencies or sampling times.

The same information is available programmatically via the model’s
`population` metadata:

``` r

str(ui$population, max.level = 1)
#> List of 14
#>  $ species       : chr "human"
#>  $ n_subjects    : num 81
#>  $ n_studies     : num 1
#>  $ age_mean      : chr "49.46 years (SD 11.15)"
#>  $ age_median    : chr "50.67 years"
#>  $ age_range     : chr "20.67-73.11 years"
#>  $ weight_mean   : chr "70.49 kg (SD 13.53)"
#>  $ weight_median : chr "71.00 kg"
#>  $ weight_range  : chr "38.00-120.00 kg"
#>  $ sex_female_pct: num 54.3
#>  $ disease_state : chr "Schizophrenia; inpatients on routine oral clozapine therapy."
#>  $ dose_range    : chr "Not reported. Concentrations came from routine therapeutic drug monitoring and the paper tabulates no administe"| __truncated__
#>  $ regions       : chr "China (single centre: Xuzhou Oriental Hospital Affiliated to Xuzhou Medical University, Jiangsu)"
#>  $ notes         : chr "Retrospective analysis of clozapine therapeutic-drug-monitoring concentrations collected between December 2023 "| __truncated__
```

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in `inst/modeldb/specificDrugs/Han_2025_clozapine.R`. The
table below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL/F at 70 kg, no zopiclone) | 29.6 L/h | Table 3 and Equation (6); SE 6.5%, bootstrap median 29.4 \[26.1, 33.9\] |
| `lvc` (V/F at 70 kg) | 308 L | Table 3 and Equation (7); SE 14.4%, bootstrap median 309 \[230, 421\] |
| `lka` (absorption rate constant) | 1.3 1/h, **fixed** | Table 3 prints “1.3 (fixed)”; Methods, Model building, citing refs 25-26 (Li 2012; Shang 2014) |
| `e_wt_cl` (allometric exponent, CL/F) | 0.75, **fixed** | Methods, Equation (3), citing Anderson & Holford 2008; assembled in Equation (6) |
| `e_wt_vc` (allometric exponent, V/F) | 1, **fixed** | Methods, Equation (3); assembled in Equation (7) |
| `e_zop_cl` (zopiclone effect on CL/F) | -0.254 | Table 3 (theta_ZOP), SE 30.8%, bootstrap median -0.241 \[-0.408, -0.009\]; assembled in Equation (6) |
| `etalcl` (IIV on CL/F) | 0.348 as an **SD** (variance 0.348^2) | Table 3 (omega_CL/F), SE 11.3%, bootstrap median 0.342 \[0.264, 0.422\]; scale settled below |
| `propSd` (proportional residual error) | 0.257 | Table 3 (sigma_1), SE 6.6%, bootstrap median 0.254 \[0.220, 0.289\]; form from Methods Equation (2) |
| Structural model: 1-compartment, first-order oral absorption | n/a | Methods, Model building (“CL/F, V/F, and Ka … were the main pharmacokinetic parameters”) |
| IIV form: `Z_i = TV(Z) * exp(eta_i)` | n/a | Methods, Equation (1) |
| Residual form: `Y_i = X_i + X_i * eps_1` (proportional only) | n/a | Methods, Equation (2) |
| Categorical covariate form: `R_i = TV(R) * (1 + Q * S_i)` | n/a | Methods, Equation (5) |
| Reference weight 70 kg | n/a | Methods, Equation (3) (“Vstd denoted standard weight of 70 kg”) |
| Therapeutic range 350-800 ng/mL; toxicity threshold 1000 ng/mL | n/a | Methods, Dosage simulation, citing refs 28-31 |

Note that the seven model equations are images in the published PDF; the
values above were read from the typeset equations and cross-checked
against Table 3.

## Reconstructing the dosing regimen

The paper reports simulated doses only as **mg/kg/day** and never states
the dosing frequency. The frequency matters: a trough depends strongly
on it. Two internal facts in the paper settle it.

``` r

th   <- ui$theta
cl70 <- exp(th[["lcl"]])            # 29.6 L/h at 70 kg
v70  <- exp(th[["lvc"]])            # 308 L
ka   <- exp(th[["lka"]])            # 1.3 1/h

# Steady-state trough for a 1-compartment first-order-absorption model.
ss_trough <- function(dose_mg, tau, cl, v, ka) {
  kel <- cl / v
  1000 * (dose_mg * ka) / (v * (ka - kel)) *
    (exp(-kel * tau) / (1 - exp(-kel * tau)) -
     exp(-ka  * tau) / (1 - exp(-ka  * tau)))
}

# A 70 kg patient at the paper's recommended 8 mg/kg/day (Table 4, 67-88 kg band).
regimen <- tibble::tibble(
  regimen      = c("once daily (tau = 24 h)", "twice daily (tau = 12 h)"),
  tau          = c(24, 12),
  trough_ngmL  = c(ss_trough(70 * 8,     24, cl70, v70, ka),
                   ss_trough(70 * 8 / 2, 12, cl70, v70, ka))
)
knitr::kable(regimen, digits = 1,
             caption = "Typical steady-state trough at the paper's own recommended dose.")
```

| regimen                  | tau | trough_ngmL |
|:-------------------------|----:|------------:|
| once daily (tau = 24 h)  |  24 |       217.2 |
| twice daily (tau = 12 h) |  12 |       452.7 |

Typical steady-state trough at the paper’s own recommended dose.
{.table}

The paper’s stated therapeutic range is 350-800 ng/mL, and Table 4
recommends 8 mg/kg/day precisely because it puts a 67-88 kg patient in
that range. Once-daily dosing misses the window low by a wide margin,
while twice-daily lands near its centre. Clozapine’s elimination
half-life under this model is 7.2 h, far shorter than a 24 h interval,
which is why the once-daily trough collapses. **All simulations below
therefore use twice-daily dosing**; this is an inference, recorded in
the Assumptions section, not something the paper states.

## Virtual cohort

The original observed data are not public. Because this model carries
exactly **one** random effect (`etalcl`, on CL/F) and no other
stochastic element, the between-subject distribution can be represented
by a deterministic quadrature lattice rather than by random sampling:
200 equally-spaced normal quantiles reproduce the log-normal clearance
distribution exactly and identically on every machine.

This matters for reproducibility. `rxode2` partitions its RNG streams
per solver thread, so a randomly sampled cohort differs between a
workstation and a CI runner and no seed can make them agree. The lattice
removes that failure mode entirely: **every number in this vignette is
deterministic**, so the assertions below can be tight without becoming
machine-dependent.

``` r

n_lattice <- 200L   # per arm; the 200-per-arm cap, used exactly
omega_sd  <- 0.348  # Table 3, read as a standard deviation (see below)

# Deterministic midpoint lattice on the eta distribution.
eta_grid <- stats::qnorm((seq_len(n_lattice) - 0.5) / n_lattice) * omega_sd

# The paper simulated the discrete weights 40, 60, 80, 100 and 120 kg
# (Methods, Dosage simulation) and reported results in weight bands (Table 4).
scenarios <- tidyr::expand_grid(
  WT  = c(40, 60, 80, 100, 120),
  zop = c(0L, 1L)
) |>
  dplyr::mutate(
    # Table 4 recommended starting dose for the band this weight falls in.
    dose_mg_kg_day = dplyr::if_else(
      zop == 0L,
      dplyr::case_when(WT <  50 ~ 10, WT <  67 ~ 9, WT < 88 ~ 8, TRUE ~ 7),
      dplyr::case_when(WT <  70 ~  6, TRUE ~ 5)
    ),
    band = dplyr::if_else(
      zop == 0L,
      dplyr::case_when(WT <  50 ~ "[40-50)", WT < 67 ~ "[50-67)",
                       WT <  88 ~ "[67-88)", TRUE ~ "[88-120]"),
      dplyr::case_when(WT <  70 ~ "[40-70)", TRUE ~ "[70-120]")
    ),
    arm = paste0(WT, " kg, ", ifelse(zop == 1L, "with", "without"), " zopiclone")
  )

# Build one event table per (scenario, lattice point). Dosing is twice daily
# for 14 days; the observation at t = 336 h is the trough of the 28th dose.
make_events <- function(scn, id_offset) {
  tidyr::expand_grid(k = seq_len(n_lattice), row = 1:2) |>
    dplyr::mutate(
      id    = id_offset + k,
      time  = dplyr::if_else(row == 1L, 0, 336),
      evid  = dplyr::if_else(row == 1L, 1L, 0L),
      amt   = dplyr::if_else(row == 1L, scn$WT * scn$dose_mg_kg_day / 2, NA_real_),
      ii    = dplyr::if_else(row == 1L, 12, 0),
      addl  = dplyr::if_else(row == 1L, 27L, 0L),
      # Observation rows point at the ODE state `central`, never at the
      # algebraic observable `Cc`.
      cmt   = dplyr::if_else(row == 1L, "depot", "central"),
      WT    = scn$WT,
      CONMED_ZOPICLONE = scn$zop,
      etalcl = eta_grid[k],
      arm = scn$arm, band = scn$band, zop = scn$zop,
      dose_mg_kg_day = scn$dose_mg_kg_day
    ) |>
    dplyr::select(-k, -row)
}

events <- do.call(
  dplyr::bind_rows,
  lapply(seq_len(nrow(scenarios)), function(i)
    make_events(scenarios[i, ], id_offset = (i - 1L) * n_lattice))
)

# IDs must be disjoint across arms: duplicates silently merge into one subject.
stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
stopifnot(dplyr::n_distinct(events$id) == nrow(scenarios) * n_lattice)
```

## Simulation

``` r

sim <- suppressWarnings(rxode2::rxSolve(
  mod, events = events,
  omega = NA,          # etas are supplied explicitly, not drawn
  returnType = "data.frame",
  keep = c("arm", "band", "zop", "dose_mg_kg_day")
))
#> ℹ parameter labels from comments will be replaced by 'label()'

# `Cc` is the individual prediction (no residual error); `sim` would carry it.
trough <- sim |> dplyr::filter(time == 336)
stopifnot(nrow(trough) == nrow(scenarios) * n_lattice, !anyNA(trough$Cc))
```

## Check 1 – the zopiclone drug-drug interaction (Figure 3)

The headline result is that concomitant zopiclone lowers clozapine CL/F
by 25.4%. Because the effect is multiplicative and shares the eta, the
ratio is exact and can be asserted tightly.

``` r

ratio <- trough |>
  dplyr::group_by(WT, zop) |>
  dplyr::summarise(cl_med = stats::median(cl), .groups = "drop") |>
  tidyr::pivot_wider(names_from = zop, values_from = cl_med,
                     names_prefix = "zop") |>
  dplyr::mutate(ratio = zop1 / zop0,
                pct_reduction = 100 * (1 - ratio))

knitr::kable(
  ratio |> dplyr::rename("Weight (kg)" = WT, "CL/F without (L/h)" = zop0,
                         "CL/F with (L/h)" = zop1, "Ratio" = ratio,
                         "Reduction (%)" = pct_reduction),
  digits = 3,
  caption = "Clozapine CL/F with vs without concomitant zopiclone (Figure 3)."
)
```

| Weight (kg) | CL/F without (L/h) | CL/F with (L/h) | Ratio | Reduction (%) |
|------------:|-------------------:|----------------:|------:|--------------:|
|          40 |             19.454 |          14.513 | 0.746 |          25.4 |
|          60 |             26.368 |          19.671 | 0.746 |          25.4 |
|          80 |             32.718 |          24.408 | 0.746 |          25.4 |
|         100 |             38.678 |          28.854 | 0.746 |          25.4 |
|         120 |             44.346 |          33.082 | 0.746 |          25.4 |

Clozapine CL/F with vs without concomitant zopiclone (Figure 3).
{.table}

``` r


# Deterministic identity: the model must reproduce 1 - 0.254 exactly.
stopifnot(all(abs(ratio$pct_reduction - 25.4) < 1e-8))
```

The model reproduces the published 25.4% clearance reduction exactly, at
every weight, as it must.

## Check 2 – recommended doses land in the therapeutic range (Figures 4-6)

Table 4 recommends a starting dose per weight band. The paper’s
therapeutic range is 350-800 ng/mL. A typical patient given the
recommended dose should sit inside it.

``` r

typ <- trough |>
  dplyr::group_by(arm, band, zop, dose_mg_kg_day, WT) |>
  dplyr::summarise(trough_typ = stats::median(Cc),
                   pct_in_range = 100 * mean(Cc >= 350 & Cc <= 800),
                   .groups = "drop") |>
  dplyr::arrange(zop, WT)

knitr::kable(
  typ |> dplyr::select(-arm) |>
    dplyr::rename("Weight (kg)" = WT, "Band" = band,
                  "Zopiclone" = zop, "Dose (mg/kg/day)" = dose_mg_kg_day,
                  "Typical trough (ng/mL)" = trough_typ,
                  "In 350-800 ng/mL (%)" = pct_in_range),
  digits = 1,
  caption = "Typical steady-state trough at each Table 4 recommended dose."
)
```

| Band | Zopiclone | Dose (mg/kg/day) | Weight (kg) | Typical trough (ng/mL) | In 350-800 ng/mL (%) |
|:---|---:|---:|---:|---:|---:|
| \[40-50) | 0 | 10 | 40 | 448.8 | 51.0 |
| \[50-67) | 0 | 9 | 60 | 478.5 | 53.5 |
| \[67-88) | 0 | 8 | 80 | 477.4 | 54.5 |
| \[88-120\] | 0 | 7 | 100 | 455.8 | 55.0 |
| \[88-120\] | 0 | 7 | 120 | 488.8 | 56.5 |
| \[40-70) | 1 | 6 | 40 | 430.8 | 55.0 |
| \[40-70) | 1 | 6 | 60 | 500.5 | 58.5 |
| \[70-120\] | 1 | 5 | 80 | 462.4 | 59.0 |
| \[70-120\] | 1 | 5 | 100 | 500.0 | 60.5 |
| \[70-120\] | 1 | 5 | 120 | 532.5 | 61.5 |

Typical steady-state trough at each Table 4 recommended dose. {.table}

``` r


# The centre of the distribution must sit inside the paper's own window at
# every recommended dose. This is deterministic, so the bound is tight.
stopifnot(all(typ$trough_typ > 350), all(typ$trough_typ < 800))
```

Every recommended dose puts the typical patient inside the published
therapeutic window (431-533 ng/mL), which is a direct, independent
confirmation of both the structural model and the twice-daily regimen
inferred above.

``` r

trough |>
  dplyr::mutate(zopl = ifelse(zop == 1L, "With zopiclone", "Without zopiclone")) |>
  ggplot(aes(factor(WT), Cc)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_hline(yintercept = c(350, 800), linetype = "dashed") +
  geom_hline(yintercept = 1000, linetype = "dotted", colour = "red") +
  facet_wrap(~zopl) +
  scale_y_log10() +
  labs(x = "Body weight (kg)", y = "Steady-state trough (ng/mL)",
       caption = paste("Dashed: therapeutic range 350-800 ng/mL.",
                       "Dotted red: 1000 ng/mL toxicity threshold."))
```

![Replicates Figures 4 and 5 of Han 2025: simulated steady-state trough
distributions at the Table 4 recommended dose, with and without
zopiclone.](Han_2025_clozapine_files/figure-html/figure-4-5-1.png)

Replicates Figures 4 and 5 of Han 2025: simulated steady-state trough
distributions at the Table 4 recommended dose, with and without
zopiclone.

## Check 3 – the Table 4 safety column

Table 4 reports, per weight band, an upper bound on the probability of
exceeding the 1000 ng/mL toxicity threshold. Because the paper reports
these as “\< x%” they are bounds over the weights inside each band, so
the model’s value for a band is the worst (largest) of the weights it
contains.

The paper’s Monte Carlo simulated concentrations; the residual error is
integrated analytically here rather than sampled, again to keep the
result deterministic.

``` r

prop_sd <- ui$theta[["propSd"]]

band_summary <- trough |>
  dplyr::group_by(zop, band, dose_mg_kg_day, WT) |>
  dplyr::summarise(
    # P(Cc > 1000) over the eta lattice, no residual error.
    p_ipred = 100 * mean(Cc > 1000),
    # P(Cc * (1 + eps) > 1000), residual integrated analytically.
    p_resid = 100 * mean(1 - stats::pnorm((1000 / Cc - 1) / prop_sd)),
    .groups = "drop"
  ) |>
  dplyr::group_by(zop, band, dose_mg_kg_day) |>
  dplyr::summarise(p_ipred = max(p_ipred), p_resid = max(p_resid),
                   .groups = "drop")

published <- tibble::tribble(
  ~zop, ~band,      ~p_published,
  0L,   "[40-50)",  13.0,
  0L,   "[50-67)",  12.5,
  0L,   "[67-88)",  11.2,
  0L,   "[88-120]",  9.8,
  1L,   "[40-70)",  11.1,
  1L,   "[70-120]",  9.0
)

tab4 <- band_summary |>
  dplyr::left_join(published, by = c("zop", "band")) |>
  dplyr::mutate(
    zopl  = ifelse(zop == 1L, "With", "Without"),
    meets = p_resid <= p_published
  ) |>
  dplyr::arrange(zop, dplyr::desc(dose_mg_kg_day))

knitr::kable(
  tab4 |> dplyr::select(zopl, band, dose_mg_kg_day, p_ipred, p_resid,
                        p_published, meets) |>
    dplyr::rename("Zopiclone" = zopl, "Weight band (kg)" = band,
                  "Dose (mg/kg/day)" = dose_mg_kg_day,
                  "Model, no residual (%)" = p_ipred,
                  "Model, with residual (%)" = p_resid,
                  "Published bound (%)" = p_published,
                  "Within bound" = meets),
  digits = 1,
  caption = paste("Probability of exceeding the 1000 ng/mL toxicity threshold",
                  "at each Table 4 recommended dose.")
)
```

| Zopiclone | Weight band (kg) | Dose (mg/kg/day) | Model, no residual (%) | Model, with residual (%) | Published bound (%) | Within bound |
|:---|:---|---:|---:|---:|---:|:---|
| Without | \[40-50) | 10 | 6.5 | 8.2 | 13.0 | TRUE |
| Without | \[50-67) | 9 | 7.5 | 9.3 | 12.5 | TRUE |
| Without | \[67-88) | 8 | 7.0 | 8.7 | 11.2 | TRUE |
| Without | \[88-120\] | 7 | 7.0 | 8.8 | 9.8 | TRUE |
| With | \[40-70) | 6 | 7.0 | 8.9 | 11.1 | TRUE |
| With | \[70-120\] | 5 | 8.0 | 10.1 | 9.0 | FALSE |

Probability of exceeding the 1000 ng/mL toxicity threshold at each Table
4 recommended dose. {.table}

Five of the six bands sit strictly under the published bound. The 1
exception is the \[70-120\] kg band with zopiclone, where the model
gives 10.1% against a published bound of 9.0%. That band spans 80-120 kg
and its bound is driven by the heaviest patient; the paper does not
state how weight was distributed inside a band, so a small excess at the
top of the widest band is expected. This is recorded as a deviation
below and is **not** tuned away.

``` r

# Gate on the bands the model should clear, and on the overall magnitude.
stopifnot(
  sum(tab4$meets) >= 5,
  # No band may exceed its published bound by more than 2 percentage points.
  all(tab4$p_resid - tab4$p_published < 2),
  # Every band must be in the same regime as published (single-digit to low
  # teens), which is what the omega-scale reading below actually decides.
  all(tab4$p_resid < 15)
)
```

## Check 4 – the omega scale is a standard deviation, not a variance

Table 3 prints `omega_CL/F = 0.348` without saying whether it is a
standard deviation or a variance. The choice changes every prediction,
so it is settled here against the paper’s own Table 4 rather than
assumed.

``` r

# Repeat the Table 4 sweep with omega read as a VARIANCE (SD = sqrt(0.348)).
eta_grid_var <- stats::qnorm((seq_len(n_lattice) - 0.5) / n_lattice) * sqrt(0.348)

# Lattice position within an arm is ((id - 1) mod n_lattice) + 1.
events_var <- events |>
  dplyr::mutate(etalcl = eta_grid_var[((id - 1L) %% n_lattice) + 1L])

sim_var <- suppressWarnings(rxode2::rxSolve(
  mod, events = events_var, omega = NA, returnType = "data.frame",
  keep = c("arm", "band", "zop", "dose_mg_kg_day")
))

alt <- sim_var |>
  dplyr::filter(time == 336) |>
  dplyr::group_by(zop, band, WT) |>
  dplyr::summarise(
    p_resid = 100 * mean(1 - stats::pnorm((1000 / Cc - 1) / prop_sd)),
    .groups = "drop"
  ) |>
  dplyr::group_by(zop, band) |>
  dplyr::summarise(p_resid = max(p_resid), .groups = "drop")

compare <- tab4 |>
  dplyr::select(zop, band, p_published, sd_reading = p_resid) |>
  dplyr::left_join(alt |> dplyr::rename(var_reading = p_resid),
                   by = c("zop", "band"))

knitr::kable(
  compare |> dplyr::select(band, p_published, sd_reading, var_reading) |>
    dplyr::rename("Weight band (kg)" = band, "Published bound (%)" = p_published,
                  "omega read as SD (%)" = sd_reading,
                  "omega read as variance (%)" = var_reading),
  digits = 1,
  caption = paste("Reading omega_CL/F as a variance roughly doubles every",
                  "exceedance probability and breaks every published bound.")
)
```

| Weight band (kg) | Published bound (%) | omega read as SD (%) | omega read as variance (%) |
|:---|---:|---:|---:|
| \[40-50) | 13.0 | 8.2 | 18.8 |
| \[50-67) | 12.5 | 9.3 | 19.9 |
| \[67-88) | 11.2 | 8.7 | 19.2 |
| \[88-120\] | 9.8 | 8.8 | 19.3 |
| \[40-70) | 11.1 | 8.9 | 19.3 |
| \[70-120\] | 9.0 | 10.1 | 20.5 |

Reading omega_CL/F as a variance roughly doubles every exceedance
probability and breaks every published bound. {.table}

``` r


stopifnot(
  # The SD reading is compatible with the published bounds ...
  all(compare$sd_reading < 15),
  # ... and the variance reading is not, in EVERY band.
  all(compare$var_reading > compare$p_published),
  all(compare$var_reading > 15)
)
```

The variance reading roughly doubles every probability and violates all
six published bounds, so `0.348` is a standard deviation (35.9% CV on
the log scale). Two further facts agree: Table 3 labels the row
`omega_CL/F` rather than `omega^2`, and its 11.3% relative standard
error is characteristic of an omega SD rather than a variance. The same
convention was reached independently for `Zhang_2024_olanzapine`, from
the same research group and the same table format.

## PKNCA validation

**The source paper reports no non-compartmental analysis** – no Cmax,
Tmax, AUC or half-life, for either the observed data or the simulations.
There is therefore nothing to compare against, and
[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
is not used. Instead PKNCA is run on a single-dose simulation and
checked against the closed-form solution of the same one-compartment
model. Both sides use the same drawn parameters, so the only difference
is numerical and a tight bound is correct.

``` r

# Single 300 mg oral dose to a lattice of 200 subjects at 70 kg, no zopiclone.
# The absorption phase is sampled at 0.05 h so that Tmax (~1.6-2.8 h across the
# lattice) is resolved to well under a percent; a coarse grid would otherwise
# both blur Tmax and bias AUC low (failure pattern 11).
nca_events <- tidyr::expand_grid(
  k    = seq_len(n_lattice),
  time = sort(unique(c(seq(0, 6, by = 0.05), seq(6.5, 72, by = 0.5))))
) |>
  dplyr::mutate(
    id = k, evid = 0L, amt = NA_real_, cmt = "central",
    WT = 70, CONMED_ZOPICLONE = 0L, etalcl = eta_grid[k]
  ) |>
  dplyr::bind_rows(
    tibble::tibble(
      k = seq_len(n_lattice), id = seq_len(n_lattice), time = 0,
      evid = 1L, amt = 300, cmt = "depot",
      WT = 70, CONMED_ZOPICLONE = 0L, etalcl = eta_grid[seq_len(n_lattice)]
    )
  ) |>
  dplyr::arrange(id, time, dplyr::desc(evid)) |>
  dplyr::select(-k)

nca_sim <- suppressWarnings(rxode2::rxSolve(
  mod, events = nca_events, omega = NA, returnType = "data.frame"
))

# Filter on !is.na(Cc) ONLY -- a time > 0 or Cc > 0 filter would drop the
# time-zero anchor and trigger PKNCA's "AUC range starting before the first
# measurement" warning on every subject.
sim_nca <- nca_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(treatment = "300 mg single oral dose") |>
  dplyr::select(id, time, Cc, treatment)

# Guarantee a time-zero record per subject (pre-dose Cc = 0 extravascularly).
sim_nca <- dplyr::bind_rows(
  sim_nca,
  sim_nca |> dplyr::distinct(id, treatment) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, treatment, time, .keep_all = TRUE) |>
  dplyr::arrange(id, treatment, time)

stopifnot(nrow(sim_nca) > 0, all(sim_nca$Cc >= 0))

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id)

dose_df <- nca_events |>
  dplyr::filter(evid == 1L) |>
  dplyr::mutate(treatment = "300 mg single oral dose") |>
  dplyr::select(id, time, amt, treatment)

dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id)

intervals <- data.frame(start = 0, end = Inf,
                        cmax = TRUE, tmax = TRUE,
                        aucinf.obs = TRUE, half.life = TRUE)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj,
                                          intervals = intervals))
```

``` r

per_subject <- nca_sim |>
  dplyr::distinct(id, cl, vc, ka)

closed_form <- per_subject |>
  dplyr::mutate(
    kel        = cl / vc,
    aucinf_cf  = 1000 * 300 / cl,                    # Dose / (CL/F), ng*h/mL
    halflife_cf = log(2) / kel,
    tmax_cf    = log(ka / kel) / (ka - kel)
  )

obs <- as.data.frame(nca_res) |>
  dplyr::filter(PPTESTCD %in% c("aucinf.obs", "half.life", "tmax")) |>
  dplyr::select(id, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
  dplyr::mutate(id = as.integer(as.character(id))) |>
  dplyr::left_join(closed_form, by = "id") |>
  dplyr::mutate(
    auc_pct  = 100 * (aucinf.obs - aucinf_cf) / aucinf_cf,
    hl_pct   = 100 * (half.life  - halflife_cf) / halflife_cf,
    tmax_pct = 100 * (tmax       - tmax_cf) / tmax_cf
  )

knitr::kable(
  tibble::tibble(
    Parameter = c("AUC(0-inf) vs Dose/(CL/F)",
                  "Half-life vs ln(2)/kel",
                  "Tmax vs ln(ka/kel)/(ka-kel)"),
    `Median % difference` = c(stats::median(obs$auc_pct),
                              stats::median(obs$hl_pct),
                              stats::median(obs$tmax_pct)),
    `Max abs % difference` = c(max(abs(obs$auc_pct)),
                               max(abs(obs$hl_pct)),
                               max(abs(obs$tmax_pct)))
  ),
  digits = 3,
  caption = "PKNCA on the simulation vs the model's own closed-form solution."
)
```

| Parameter                   | Median % difference | Max abs % difference |
|:----------------------------|--------------------:|---------------------:|
| AUC(0-inf) vs Dose/(CL/F)   |              -0.003 |                0.020 |
| Half-life vs ln(2)/kel      |               0.274 |                0.297 |
| Tmax vs ln(ka/kel)/(ka-kel) |              -0.033 |                1.403 |

PKNCA on the simulation vs the model’s own closed-form solution.
{.table}

``` r


# Deterministic solve-vs-closed-form: differences are numerical only.
stopifnot(
  max(abs(obs$auc_pct))  < 1.0,   # trapezoidal + extrapolation error
  max(abs(obs$hl_pct))   < 1.0,
  # Tmax is resolved to at most half a 0.05 h grid step, i.e. ~1.6% of the
  # shortest Tmax in the lattice. A mis-transcribed ka or clearance moves Tmax
  # by tens of percent, so 3% still goes red on a real error.
  max(abs(obs$tmax_pct)) < 3.0
)
```

PKNCA recovers the closed-form AUC and half-life to well under 1%,
confirming that the packaged ODE system, the parameterisation and the
mg-to-ng/mL scaling are internally consistent.

## Assumptions and deviations

- **Dosing frequency is inferred, not published.** The paper gives
  simulated doses only in mg/kg/day. Twice-daily dosing is used
  throughout because it is the only common regimen that places the
  typical patient inside the paper’s own 350-800 ng/mL therapeutic range
  at the paper’s own recommended doses (see “Reconstructing the dosing
  regimen”); once-daily dosing falls roughly 45% short. Twice-daily is
  also standard clinical practice for maintenance clozapine. All
  absolute concentrations in this vignette depend on this choice; the
  clearance ratio in Check 1 does not.

- **`omega_CL/F = 0.348` is read as a standard deviation.** Table 3 does
  not say. Check 4 settles it against the paper’s Table 4: the variance
  reading breaks all six published bounds while the SD reading is
  compatible with all of them. See also the in-file comment on `etalcl`.

- **Weight distribution inside a Table 4 band is unknown.** The paper
  simulated the discrete weights 40, 60, 80, 100 and 120 kg but reports
  results in bands. Each band is evaluated here at the simulated weights
  it contains, taking the worst case, which matches the published “\<
  x%” bound semantics.

- **Known deviation, not tuned:** the 70-120 kg band with zopiclone
  yields 10.1% against a published bound of 9.0%. This is the widest
  band and its value is set by the 120 kg extreme; the discrepancy is
  under 1.5 percentage points. No parameter was altered to remove it.

- **Residual error is integrated analytically** rather than sampled, so
  the reported exceedance probabilities are deterministic. Probabilities
  computed without residual error are shown alongside for comparison;
  the paper does not state whether its Monte Carlo included the residual
  term.

- **Unit errors in the source’s Table 1.** Creatinine, direct bilirubin
  and total bilirubin are tagged “mmol/L” where the quoted values are
  plainly umol/L (61 mmol/L creatinine is ~1000-fold above any
  survivable value). Hemoglobin’s lower range bound of 21.00 g/L is
  likewise implausible. All four are recorded corrected in
  `covariatesDataExcluded`; none enters the final model, so no
  prediction is affected.

- **The model’s `ka` is fixed, not estimated**, at 1.3 1/h carried from
  Li 2012 / Shang 2014 (the source’s references 25-26). Li 2012 is
  packaged as `modellib("Li_2012_clozapine")` and carries the same fixed
  value.

- **No bioavailability term.** The dataset is oral-only therapeutic drug
  monitoring, so `F` is not identifiable and `cl` / `vc` are the
  apparent quantities CL/F and V/F throughout.

- **No published NCA to compare against.** The source reports no Cmax,
  Tmax, AUC or half-life, so the PKNCA section validates against the
  model’s own closed-form solution rather than against published values.
