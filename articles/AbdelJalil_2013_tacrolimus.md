# Tacrolimus (Abdel Jalil 2013)

## Model and source

- Citation: Abdel Jalil MH, Hawwa AF, McKiernan PJ, Shields MD, McElnay
  JC. Population pharmacokinetic and pharmacogenetic analysis of
  tacrolimus in paediatric liver transplant patients. Br J Clin
  Pharmacol. 2014;77(1):130-140. <doi:10.1111/bcp.12174>
- Description: One-compartment population PK model for oral tacrolimus
  in paediatric liver transplant recipients (Abdel Jalil 2013), with
  first-order absorption (ka fixed at the literature value 4.5 1/h), an
  apparent volume of distribution fixed at the literature value 30 L/kg,
  allometric (WT/13.2 kg)^0.75 scaling on apparent clearance with the
  theory-based exponent fixed, multiplicative exponential effects of
  time post-transplantation (days) and CYP3A5\*1 carrier status on CL/F,
  exponential (log-normal) inter-individual variability on CL/F, and
  proportional residual error.
- Article: <https://doi.org/10.1111/bcp.12174>

## Population

The model was developed from 628 tacrolimus pre-dose whole-blood trough
concentrations collected retrospectively from therapeutic drug
monitoring records of 43 paediatric patients during their first year
after liver transplantation at Birmingham Children’s Hospital, UK (Abdel
Jalil 2013 Table 2). Mean (range) age was 5 (0.65-17.56) years and mean
(range) body weight 21.6 (6.06-69.95) kg; 44.2% (19/43) of subjects were
female. The race distribution was 31 Caucasian (72.1%), 9 Asian (20.9%),
and 3 African Caribbean (7.0%). Tacrolimus was administered orally as
capsules or oral suspension, predominantly twice daily, in combination
with low-dose steroids; the dose range was 0.4-30 mg/day (0.03-1.23
mg/kg/day) and the median time since transplantation 153 days (range
15-364). Trough concentrations were measured by enzyme immunoassay
(Abbott IMx; LLQ 2 ng/mL, linear to 30 ng/mL). The cohort was genotyped
for CYP3A5*3 and three ABCB1 SNPs; the CYP3A5 genotype distribution was*
1/*1 n = 2 (4.7%),* 1/*3 n = 9 (20.9%),* 3/\*3 n = 32 (74.4%), in
Hardy-Weinberg equilibrium.

The same information is available programmatically via
`readModelDb("AbdelJalil_2013_tacrolimus")$population`.

## Source trace

Every parameter in the model file carries an inline source-location
comment. The table below collects the entries in one place.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lka` (ka, fixed) | 4.5 1/h | Methods ‘Selecting the base model’ (literature reference 9) |
| `lvc` (V/F per kg, fixed) | 30 L/kg | Methods ‘Selecting the base model’ (literature references 7, 8) |
| `lcl` (CL/F at WT 13.2 kg, *3/*3, POD 0) | 12.92 L/h | Table 5, theta_1 row (RSE 4.8%) |
| `e_wt_cl` (allometric exponent, fixed) | 0.75 | Methods ‘Selecting the base model’ (theory-based allometric exponent) |
| `e_pod_cl` (theta_2) | -0.00158 /d | Table 5, theta_2 row (RSE 13.8%) |
| `e_cyp3a5_expr_cl` (theta_3) | 0.4282 | Table 5, theta_3 row (RSE 25.5%) |
| IIV CL/F (omega^2) | 0.16 | Table 5, omega^2 row (RSE 19.8%; ~40% CV) |
| Proportional residual error (sigma^2) | 0.125 | Table 5, sigma^2 row (RSE 11.9%; ~35.4% CV) |
| Compartment structure (1-cmt oral, first-order absorption + elimination) | – | Results ‘Population pharmacokinetic analysis’ paragraph |
| Final CL/F covariate equation | `12.92 * (WT/13.2)^0.75 * exp(-0.00158 * TPT) * exp(0.4282 * CYP3A5)` | Results ‘Covariate selection and final pharmacokinetic model’ equation; Abstract Results equation |
| Reference weight (study median) | 13.2 kg | Methods ‘Selecting the base model’ allometric paragraph |

## Virtual cohort

The published dataset is not openly available, so the virtual cohort
below mirrors the demographics in Abdel Jalil 2013 Table 2 and the
CYP3A5 genotype distribution stated in Results. Two genotype strata are
built so the time-related clearance simulation reproduces Figure 3 of
the paper.

``` r

set.seed(20130506)

n_per_geno <- 100L

make_cohort <- function(n, cyp3a5_expr, label, id_offset = 0L) {
  tibble(
    id          = id_offset + seq_len(n),
    # Weight distribution centred on the study median 13.2 kg; the geometric
    # spread is chosen so the 5th-95th percentile spans roughly the observed
    # 6-70 kg range. A log-normal draw keeps weights strictly positive.
    WT          = pmin(pmax(exp(rnorm(n, mean = log(13.2), sd = 0.55)),
                            5.0), 75.0),
    POD         = 100,   # POD set per-scenario; placeholder for cohort construction
    CYP3A5_EXPR = cyp3a5_expr,
    cohort      = label
  )
}

demo <- bind_rows(
  make_cohort(n_per_geno, cyp3a5_expr = 0L,
              label = "*3/*3 (nonexpresser)",
              id_offset = 0L),
  make_cohort(n_per_geno, cyp3a5_expr = 1L,
              label = "*1 carrier",
              id_offset = n_per_geno)
)
stopifnot(!anyDuplicated(demo$id))
```

## Simulation

Two scenarios are simulated:

1.  A typical-value clearance sweep across post-transplant day for two
    genotype groups (Figure 3 of the paper).
2.  A stochastic steady-state simulation of twice-daily oral tacrolimus
    over five days, used for the trough comparison and the PKNCA tables.

``` r

mod <- rxode2::rxode2(readModelDb("AbdelJalil_2013_tacrolimus"))
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_typical <- mod |> rxode2::zeroRe()
```

### Figure 3 – time-related changes in apparent clearance

Figure 3 shows the time-related changes in tacrolimus apparent clearance
for two typical individuals (one CYP3A5*3/*3, one CYP3A5\*1 carrier)
over the first post-transplant year, at a fixed dose of 2.3 mg/day (0.17
mg/kg every 12 h, paper’s caption). Apparent clearance is a
deterministic function of weight, POD, and CYP3A5 genotype in this
model, so a single typical individual per genotype recovers the figure.

``` r

pod_grid <- seq(0, 365, by = 5)

clearance_by_pod <- function(WT, CYP3A5_EXPR, pod_grid,
                              theta_cl = 12.92, theta_pod = -0.00158,
                              theta_cyp = 0.4282, ref_wt = 13.2,
                              allo = 0.75, omega_sq = 0.16) {
  tv <- theta_cl * (WT / ref_wt)^allo *
        exp(theta_pod * pod_grid) *
        exp(theta_cyp * CYP3A5_EXPR)
  # 95% IIV interval from a log-normal eta with variance omega_sq
  sd_log <- sqrt(omega_sq)
  data.frame(
    POD     = pod_grid,
    CL_typ  = tv,
    CL_low  = tv * exp(qnorm(0.025) * sd_log),
    CL_high = tv * exp(qnorm(0.975) * sd_log)
  )
}

fig3 <- bind_rows(
  clearance_by_pod(WT = 13.2, CYP3A5_EXPR = 0, pod_grid = pod_grid) |>
    mutate(cohort = "*3/*3 (nonexpresser)"),
  clearance_by_pod(WT = 13.2, CYP3A5_EXPR = 1, pod_grid = pod_grid) |>
    mutate(cohort = "*1 carrier")
)

ggplot(fig3, aes(POD, CL_typ, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = CL_low, ymax = CL_high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365)) +
  labs(x = "Time post-transplantation (days)",
       y = "Apparent oral clearance CL/F (L/h)",
       color = "CYP3A5 genotype", fill = "CYP3A5 genotype",
       title = "Time-related changes in tacrolimus apparent clearance",
       caption = "Replicates Figure 3 of Abdel Jalil 2013 (typical paediatric subject at WT = 13.2 kg).")
```

![Replicates Figure 3 of Abdel Jalil 2013: time-related changes in
tacrolimus apparent clearance for a typical paediatric patient by CYP3A5
genotype. Solid lines are typical-value clearance; the shaded ribbon is
the 95% confidence range derived from the model IIV (omega^2 =
0.16).](AbdelJalil_2013_tacrolimus_files/figure-html/figure-3-1.png)

Replicates Figure 3 of Abdel Jalil 2013: time-related changes in
tacrolimus apparent clearance for a typical paediatric patient by CYP3A5
genotype. Solid lines are typical-value clearance; the shaded ribbon is
the 95% confidence range derived from the model IIV (omega^2 = 0.16).

### Worked-example dose calculations from the paper

The paper closes its Results section with worked dose calculations for a
13.2 kg CYP3A5 nonexpresser at 14 and 40 weeks post-transplant,
targeting a 5 ng/mL trough. The table below replicates those
calculations using the model file’s parameter values to confirm the
model reproduces the paper-quoted clearances exactly.

``` r

worked_cl <- function(POD_days, CYP3A5_EXPR, WT = 13.2,
                      theta_cl = 12.92, theta_pod = -0.00158,
                      theta_cyp = 0.4282, ref_wt = 13.2, allo = 0.75) {
  theta_cl * (WT / ref_wt)^allo *
    exp(theta_pod * POD_days) *
    exp(theta_cyp * CYP3A5_EXPR)
}

scenarios <- tibble::tibble(
  scenario = c("13.2 kg, *3/*3, 14 weeks post-transplant",
               "13.2 kg, *3/*3, 40 weeks post-transplant",
               "13.2 kg, *1 carrier, 14 weeks post-transplant",
               "13.2 kg, *1 carrier, 40 weeks post-transplant"),
  POD      = c(98, 280, 98, 280),
  CYP3A5_EXPR = c(0, 0, 1, 1)
) |>
  rowwise() |>
  mutate(
    CL_predicted_Lph = round(worked_cl(POD, CYP3A5_EXPR), 2),
    Dose_mg_per_day  = round(CL_predicted_Lph * 5 / 1000 * 24, 2)
  ) |>
  ungroup()

paper_quoted <- tibble::tibble(
  scenario = c("13.2 kg, *3/*3, 14 weeks post-transplant",
               "13.2 kg, *3/*3, 40 weeks post-transplant",
               "13.2 kg, *1 carrier, 14 weeks post-transplant",
               "13.2 kg, *1 carrier, 40 weeks post-transplant"),
  CL_paper_Lph = c(11.05, 8.29, NA_real_, NA_real_),
  Dose_paper_mg_per_day = c(1.3, 1.0, 2.0, 1.5)
)

knitr::kable(
  dplyr::left_join(scenarios, paper_quoted, by = "scenario"),
  caption = "Worked-example apparent clearance and daily dose to maintain a 5 ng/mL trough. Paper-quoted values are from Abdel Jalil 2013 Results, paragraph after Table 4."
)
```

| scenario | POD | CYP3A5_EXPR | CL_predicted_Lph | Dose_mg_per_day | CL_paper_Lph | Dose_paper_mg_per_day |
|:---|---:|---:|---:|---:|---:|---:|
| 13.2 kg, *3/*3, 14 weeks post-transplant | 98 | 0 | 11.07 | 1.33 | 11.05 | 1.3 |
| 13.2 kg, *3/*3, 40 weeks post-transplant | 280 | 0 | 8.30 | 1.00 | 8.29 | 1.0 |
| 13.2 kg, \*1 carrier, 14 weeks post-transplant | 98 | 1 | 16.98 | 2.04 | NA | 2.0 |
| 13.2 kg, \*1 carrier, 40 weeks post-transplant | 280 | 1 | 12.74 | 1.53 | NA | 1.5 |

Worked-example apparent clearance and daily dose to maintain a 5 ng/mL
trough. Paper-quoted values are from Abdel Jalil 2013 Results, paragraph
after Table 4. {.table}

### Stochastic steady-state simulation

A second simulation gives a 12-hour steady-state dosing interval for
each genotype stratum at POD = 100 days with twice-daily dosing of 0.1
mg/kg. This is used for the PKNCA tables and the trough comparison.

``` r

build_events <- function(demo, dose_mg_per_kg = 0.1, sim_hours = 120,
                          POD_value = 100) {
  demo_pod <- demo |> mutate(POD = POD_value,
                              amt_dose = dose_mg_per_kg * WT)
  doses <- demo_pod |>
    mutate(amt = amt_dose, evid = 1L, cmt = "depot",
           ii = 12, addl = 9L, time = 0) |>
    select(id, time, amt, evid, cmt, ii, addl, cohort,
           WT, POD, CYP3A5_EXPR)
  # Observation grid: dense over the last (steady-state) interval, 96-108 h.
  obs_times <- sort(unique(c(0,
                              seq(96, 96 + 12, by = 0.5))))
  obs <- demo_pod |>
    select(id, cohort, WT, POD, CYP3A5_EXPR) |>
    tidyr::crossing(time = obs_times) |>
    mutate(amt = NA_real_, evid = 0L, cmt = NA_character_,
           ii = NA_real_, addl = NA_integer_)
  bind_rows(doses, obs) |>
    arrange(id, time, desc(evid))
}

events_ss <- build_events(demo, dose_mg_per_kg = 0.1, POD_value = 100)
```

``` r

sim_iiv <- rxode2::rxSolve(
  mod, events = events_ss,
  keep   = c("cohort"),
  nStud  = 1
) |> as.data.frame()

sim_typ <- rxode2::rxSolve(
  mod_typical, events = events_ss,
  keep = c("cohort")
) |> as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> Warning: multi-subject simulation without without 'omega'
```

## PKNCA validation

A steady-state NCA over the last (10th) dosing interval gives Cmax,
Tmax, Cmin, AUC0-12, and Cavg by CYP3A5 genotype. Abdel Jalil 2013 only
reports pre-dose trough concentrations (no full profile), so the NCA
outputs here are reference values for downstream users who want a
single-snapshot characterisation of the simulated steady-state interval.
Cmin is the value that maps to the paper’s reported trough measurements.

``` r

last_dose_time <- 96  # 10th dose at t=96h; steady-state window 96-108 h

nca_window <- sim_iiv |>
  filter(time >= last_dose_time, time <= last_dose_time + 12) |>
  mutate(time_after_dose = time - last_dose_time) |>
  filter(!is.na(Cc)) |>
  select(id, time = time_after_dose, Cc, cohort)

dose_df <- demo |>
  mutate(time = 0, amt = 0.1 * WT) |>
  select(id, time, amt, cohort)

conc_obj <- PKNCA::PKNCAconc(nca_window, Cc ~ time | cohort + id,
                             concu = "ng/mL", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(dose_df,    amt ~ time | cohort + id,
                             doseu = "mg")

intervals <- data.frame(start = 0, end = 12,
                        cmax = TRUE, tmax = TRUE, cmin = TRUE,
                        auclast = TRUE, cav = TRUE, ctrough = TRUE)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- suppressMessages(suppressWarnings(PKNCA::pk.nca(nca_data)))

nca_summary <- summary(nca_res)
knitr::kable(
  nca_summary,
  caption = "Steady-state NCA on the simulated cohort (12 h dosing interval at POD = 100 days, 0.1 mg/kg twice daily, last interval after 10 doses), stratified by CYP3A5 genotype."
)
```

| Interval Start | Interval End | cohort | N | AUClast (h\*ng/mL) | Cmax (ng/mL) | Cmin (ng/mL) | Tmax (h) | Cav (ng/mL) | Ctrough (ng/mL) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|
| 0 | 12 | \*1 carrier | 100 | 76.1 \[40.3\] | 7.95 \[30.8\] | 4.75 \[53.2\] | 1.00 \[1.00, 1.00\] | 6.34 \[40.3\] | 4.80 \[54.3\] |
| 0 | 12 | *3/*3 (nonexpresser) | 100 | 107 \[37.9\] | 10.5 \[30.7\] | 7.27 \[46.0\] | 1.00 \[1.00, 1.00\] | 8.96 \[37.9\] | 7.44 \[47.8\] |

Steady-state NCA on the simulated cohort (12 h dosing interval at POD =
100 days, 0.1 mg/kg twice daily, last interval after 10 doses),
stratified by CYP3A5 genotype. {.table}

### Comparison against published trough

Abdel Jalil 2013 Table 2 reports a study-population mean tacrolimus
trough concentration of 8.93 ng/mL (range 1.2-26.4 ng/mL) collected
approximately 12 hours after a dose, pooled across all post-transplant
days and genotypes. The table below compares that observed mean and
range with the simulated steady-state Cmin (the model output that
matches the paper’s trough definition), stratified by CYP3A5 genotype.
Because the published distribution is pooled across all genotypes and
POD values, the comparison is qualitative; the typical dose chosen for
the simulation (0.1 mg/kg twice daily) is in the middle of the paper’s
observed dose range (0.03-1.23 mg/kg/day).

``` r

trough_sim <- sim_iiv |>
  filter(time == last_dose_time) |>
  group_by(cohort) |>
  summarise(Q10 = quantile(Cc, 0.10),
            median = quantile(Cc, 0.50),
            Q90 = quantile(Cc, 0.90),
            .groups = "drop")

trough_typ <- sim_typ |>
  filter(time == last_dose_time) |>
  group_by(cohort) |>
  summarise(typical = mean(Cc), .groups = "drop")

tbl <- tibble::tibble(
  metric = c("Abdel Jalil 2013 Table 2 observed trough (mean, range)",
             "Simulated typical-value trough, *3/*3 (ng/mL)",
             "Simulated typical-value trough, *1 carrier (ng/mL)",
             "Simulated cohort median, *3/*3 (10-90 pct)",
             "Simulated cohort median, *1 carrier (10-90 pct)"),
  value  = c(sprintf("%.2f (range %.1f-%.1f) ng/mL", 8.93, 1.2, 26.4),
             sprintf("%.2f",
                     trough_typ$typical[trough_typ$cohort == "*3/*3 (nonexpresser)"]),
             sprintf("%.2f",
                     trough_typ$typical[trough_typ$cohort == "*1 carrier"]),
             sprintf("%.2f (%.2f-%.2f)",
                     trough_sim$median[trough_sim$cohort == "*3/*3 (nonexpresser)"],
                     trough_sim$Q10[trough_sim$cohort    == "*3/*3 (nonexpresser)"],
                     trough_sim$Q90[trough_sim$cohort    == "*3/*3 (nonexpresser)"]),
             sprintf("%.2f (%.2f-%.2f)",
                     trough_sim$median[trough_sim$cohort == "*1 carrier"],
                     trough_sim$Q10[trough_sim$cohort    == "*1 carrier"],
                     trough_sim$Q90[trough_sim$cohort    == "*1 carrier"]))
)
knitr::kable(tbl, caption = "Simulated steady-state trough vs Abdel Jalil 2013 Table 2 observed trough range.")
```

| metric | value |
|:---|:---|
| Abdel Jalil 2013 Table 2 observed trough (mean, range) | 8.93 (range 1.2-26.4) ng/mL |
| Simulated typical-value trough, *3/*3 (ng/mL) | 7.64 |
| Simulated typical-value trough, \*1 carrier (ng/mL) | 4.95 |
| Simulated cohort median, *3/*3 (10-90 pct) | 7.51 (4.05-12.16) |
| Simulated cohort median, \*1 carrier (10-90 pct) | 4.94 (2.60-8.86) |

Simulated steady-state trough vs Abdel Jalil 2013 Table 2 observed
trough range. {.table}

The simulated cohort medians sit well within the paper’s reported range
(1.2-26.4 ng/mL), and the *1-carrier cohort median is approximately
exp(-0.4282) ~ 65% of the* 3/*3 cohort median, matching the paper’s
CYP3A5*1-carrier multiplier on CL/F.

## Assumptions and deviations

- **Volume of distribution and absorption rate are fixed to literature
  values, not estimated.** Per Abdel Jalil 2013 Methods, V/F was fixed
  at 30 L/kg (literature references 7 and 8 in the paper) and ka at 4.5
  1/h (literature reference 9). The paper attributes this to the data
  consisting only of pre-dose troughs, which do not contain information
  about the distribution phase or the absorption rate. The model file
  uses `fixed()` on both `lvc` and `lka` so the provenance is explicit.
- **Allometric exponent fixed at the theory-based 0.75.** The paper
  attempted to estimate the exponent as an additional theta but found it
  did not improve the fit (Abdel Jalil 2013 Results “Population
  pharmacokinetic analysis” paragraph), so the theory-based value was
  retained. The model file uses `fixed(0.75)` on `e_wt_cl`.
- **Reference weight is the study median 13.2 kg, NOT 70 kg.** Abdel
  Jalil 2013 standardised allometric clearance on the median weight of
  the paediatric study population. Users porting this model into pooled
  adult+paediatric simulations should keep the 13.2 kg reference unless
  they also re-fit the typical-value clearance to a 70 kg adult.
- **CYP3A5 expresser status pools *1/*1 and *1/*3.** The paper coded a
  single CYP3A5 indicator (= 1 if at least one *1 allele, 0 if* 3/*3),
  consistent with the canonical `CYP3A5_EXPR` covariate. The cohort had
  only 2* 1/\*1 homozygotes (4.7%), so the homozygote vs heterozygote
  distinction was not separately estimable.
- **Inter-occasion variability (IOV) is omitted from the model.** Abdel
  Jalil 2013 Results state that IOV on clearance was below 15% and did
  not result in any significant change in the parameter estimates, so it
  was dropped from the final model. The model file therefore has no IOV
  slots.
- **No ABCB1 covariate effect on PK.** The paper screened C1236T,
  G2677T, and C3435T (and the T-T-T haplotype) as potential covariates
  on CL/F but none was retained in the final model. They are documented
  in the model file’s `covariatesDataExcluded` list for provenance.
- **Time post-transplantation effect is monotonic.** The paper’s data
  ranged from 15-364 days post-transplant; the model’s exponential decay
  on CL/F is identifiable only within that window. Users simulating POD
  = 0 (day of transplant) should be aware the model has no support for
  the early-period transient changes that Fukudo et al. (the paper’s
  reference 11) reported within the first 50 days post-transplant.
- **Concentrations measured by enzyme immunoassay, not LC-MS/MS.**
  Tacrolimus immunoassays cross-react with metabolites (estimated \<
  minimum detectable sensitivity in adult kidney transplant patients,
  per Abdel Jalil 2013 reference 33). Users predicting concentrations to
  compare against an LC-MS/MS-measured dataset should expect a small
  positive bias from the immunoassay-trained typical-value estimate.
- **Vignette uses 100 subjects per CYP3A5 stratum.** This is small
  enough to render the vignette well under 5 minutes (the pkgdown gate)
  but large enough to give stable percentiles for the steady-state
  trough comparison.
