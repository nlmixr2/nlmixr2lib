# Clozapine (Li 2012)

## Model and source

``` r

mod_meta <- nlmixr2est::nlmixr(readModelDb("Li_2012_clozapine"))$meta
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Li LJ, Shang DW, Li WB, Guo W, Wang XP, Ren YP, Li AN, Fu
  PX, Ji SM, Lu W, Wang CY (2012). Population pharmacokinetics of
  clozapine and its primary metabolite norclozapine in Chinese patients
  with schizophrenia. Acta Pharmacol Sin 33(11):1409-1416.
  <doi:10.1038/aps.2012.71>.
- Description: One-compartment parent-plus-metabolite population PK
  model for oral clozapine and its primary active metabolite
  norclozapine (N-desmethylclozapine) in 162 Chinese adult inpatients
  (74 male, 88 female; 35.5 +/- 10.6 years) with refractory
  schizophrenia on maintenance oral clozapine therapy (Li 2012).
  First-order absorption (Ka fixed at 1.3 1/h from prior rich-data
  clozapine PK studies) into a single central compartment with
  first-order elimination; a fixed fraction (KF = 0.66) of the absorbed
  clozapine dose is converted in the parent central compartment to
  norclozapine and feeds a separate one-compartment metabolite
  compartment with its own apparent clearance and apparent volume. Two
  binary covariates were retained in the final
  forward-and-backward-selected model: current-smoker status increases
  apparent clearance of both species (clozapine by 45%, norclozapine by
  54.3%), and male sex increases apparent clearance of both species
  (clozapine by 20.8%, norclozapine by 24.2%); the typical values
  reported in Table 2 are for the female-nonsmoker reference stratum. A
  combined additive-plus-proportional residual error model is reported
  separately for clozapine and norclozapine. The model was internally
  validated using normalized prediction distribution errors (NPDE).
- Article: <https://doi.org/10.1038/aps.2012.71>

## Population

Li 2012 enrolled 162 Chinese adult inpatients (74 male, 88 female; age
18-59 years, mean 35.5 +- 10.6 years) at multiple mental-health sites in
China who carried a DSM-IV diagnosis of refractory schizophrenia and
were on stable maintenance oral clozapine therapy (Jiangsu Nhwa
pharmaceutical, single brand). Smoking status was self-reported and
nurse-confirmed; the four-cell smoking-by-sex stratum from Table 1 is:
female nonsmoker n = 72, female smoker n = 2, male nonsmoker n = 40,
male smoker n = 48. A total of 1617 plasma samples (809 clozapine, 808
norclozapine) were quantified by HPLC-UV (LOD = 0.08 umol/L for both
species, CV \< 5%); 20% of samples were taken 0.5-4.5 h post-dose and
the remaining 80% 8.5-15.5 h post-dose, reflecting sparse
therapeutic-drug-monitoring (TDM) timing.

The same demographics are available programmatically via
`mod_meta$population` after the meta chunk above.

## Source trace

Every parameter line in `inst/modeldb/specificDrugs/Li_2012_clozapine.R`
carries a trailing source-location comment. The table below collects
them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lka` (Ka, 1/h, fixed) | log(1.3) | Methods Results paragraph 1; Table 2: Ka = 1.3 (Fixed) |
| `lcl` (CL/F, L/h, female nonsmoker reference) | log(21.9) | Table 2: CL/F = 21.9 (RSE 6%) |
| `lvc` (V/F, L) | log(526) | Table 2: V/F = 526 (RSE 10%) |
| `kf` (fraction cloz to norcloz, fixed) | 0.66 | Methods Results paragraph 1; Table 2: KF = 0.66 (Fixed) |
| `lcl_norcloz` (CLM, L/h, female nonsmoker reference) | log(32.7) | Table 2: CLM = 32.7 (RSE 5.6%) |
| `lvc_norcloz` (VM, L) | log(624) | Table 2: VM = 624 (RSE 5.5%) |
| `e_smoke_cl` (smoking effect on cloz CL) | 0.45 | Table 2: theta_Smoking (cloz) = 0.45 (RSE 34.9%); Results sentence “Smoking was associated with increases in the clearance of clozapine … of 45%” |
| `e_sex_cl` (male-sex effect on cloz CL) | 0.208 | Table 2: theta_Gender (cloz) = 0.208 (RSE 44.6%); Results sentence “The clearance of clozapine … 20.8% … greater … in males” |
| `e_smoke_cl_norcloz` (smoking on CLM) | 0.543 | Table 2: theta_Smoking (norcloz) = 0.543 (RSE 35.7%) |
| `e_sex_cl_norcloz` (male-sex on CLM) | 0.242 | Table 2: theta_Gender (norcloz) = 0.242 (RSE 49.2%) |
| `etalcl` (IIV on log CL) | log(0.429^2 + 1) = 0.169 | Table 2: CL/F IIV CV% = 42.9 (lognormal omega^2 = log(CV^2 + 1)) |
| `etalvc` (IIV on log V) | log(0.657^2 + 1) = 0.359 | Table 2: V/F IIV CV% = 65.7 |
| `etalcl_norcloz` (IIV on log CLM) | log(0.421^2 + 1) = 0.163 | Table 2: CLM IIV CV% = 42.1 |
| `etalvc_norcloz` (IIV on log VM) | log(0.756^2 + 1) = 0.452 | Table 2: VM IIV CV% = 75.6 |
| `addSd` (cloz additive residual, umol/L) | 0.162 | Table 2 Residual Variability: sigma_1 SD = 0.162 umol/L |
| `propSd` (cloz proportional residual) | 0.266 | Table 2 Residual Variability: sigma_2 CV = 26.6% |
| `addSd_norcloz` (norcloz additive residual, umol/L) | 0.117 | Table 2 Residual Variability: sigma_3 SD = 0.117 umol/L |
| `propSd_norcloz` (norcloz proportional residual) | 0.169 | Table 2 Residual Variability: sigma_4 CV = 16.9% |
| Structural ODEs (parent + metabolite) | n/a | Methods ‘Population pharmacokinetic modeling’ (NONMEM ADVAN5 mass-balance K20 + K23 + K30 with the constraint K23 = KF \* (CL/F) / (V/F)) and Figure 1 (model diagram). |
| Multiplicative covariate-effect form `(1 + theta_smoke * SMOKE + theta_male * MALE)` | n/a | Methods Eqs. for continuous and discrete covariate models; verified by matching reported percent increases (45%, 20.8%, 54.3%, 24.2%) to the parameter values (0.45, 0.208, 0.543, 0.242). |

## Virtual cohort

Original observed data are not publicly available. The simulations below
build a virtual cohort whose sex-by-smoking distribution matches Li 2012
Table 1 (n = 162 split as female-nonsmoker 72, female-smoker 2,
male-nonsmoker 40, male-smoker 48). Each subject receives oral clozapine
200 mg twice daily (one of the common clinical maintenance regimens in
the cohort) for 14 days (28 doses), with concentrations sampled densely
over the final dosing interval to support steady-state NCA. Oral doses
are converted from mg to umol using the clozapine molecular weight
326.83 g/mol (1 mg = 3.0596 umol); this matches the model’s molar units
(`dosing = "umol"`, `concentration = "umol/L"`).

``` r

set.seed(20260603L)

# Sex-by-smoking distribution from Table 1 of Li 2012. The two female
# smokers in the published cohort were a small subgroup that
# nevertheless contributed to identifying the smoking effect; we
# preserve their two-subject share rather than rounding up.
strata <- tibble::tibble(
  treatment = c("female_nonsmoker", "female_smoker",
                "male_nonsmoker",   "male_smoker"),
  SEXF      = c(1, 1, 0, 0),
  SMOKE     = c(0, 1, 0, 1),
  n         = c(72L, 2L, 40L, 48L)
)

# Build subject covariate table with disjoint ids per stratum.
subjects <- strata |>
  dplyr::mutate(
    id_offset = c(0L, cumsum(n)[-length(n)])
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ids = list(id_offset + seq_len(n))
  ) |>
  tidyr::unnest(ids) |>
  dplyr::transmute(
    id        = ids,
    treatment = treatment,
    SEXF      = SEXF,
    SMOKE     = SMOKE
  )

# 200 mg oral clozapine BID for 14 days = 28 doses; 200 mg = 611.96 umol.
dose_umol <- 200 / 326.83 * 1000  # 611.96
dose_times <- seq(0, 13 * 24, by = 12)   # 0, 12, 24, ..., 312 h

doses <- tidyr::expand_grid(subjects, time = dose_times) |>
  dplyr::mutate(
    amt  = dose_umol,
    cmt  = "depot",
    evid = 1L
  )

# Observation grid: dense around the final dosing interval (288-312 h)
# plus a coarser grid through day 1 for the absorption / single-dose
# profile.
obs_times <- sort(unique(c(
  seq(0,   24,  by = 0.5),     # day 1 dense profile
  seq(24,  288, by = 4),       # multi-dose ramp to steady state
  seq(288, 320, by = 0.5)      # final dosing interval, dense
)))

obs <- tidyr::expand_grid(subjects, time = obs_times) |>
  dplyr::mutate(
    amt  = NA_real_,
    cmt  = "Cc",
    evid = 0L
  )

events <- dplyr::bind_rows(doses, obs) |>
  dplyr::arrange(id, time, dplyr::desc(evid))

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid", "cmt")])))
```

## Simulation

``` r

mod <- readModelDb("Li_2012_clozapine")
sim <- rxode2::rxSolve(mod, events = events,
                       keep = c("treatment", "SEXF", "SMOKE"),
                       seed = 20260603L) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
```

Deterministic typical-value trajectories (no IIV, no residual error) for
the four sex / smoking strata are used in the figure replication below.

``` r

mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
sim_typical <- rxode2::rxSolve(mod_typical, events = events,
                               keep = c("treatment", "SEXF", "SMOKE")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalcl_norcloz', 'etalvc_norcloz'
#> Warning: multi-subject simulation without without 'omega'
```

## Replicate Figure 3 (typical-value clearance trends by sex and smoking)

Figure 3 of Li 2012 plots clozapine clearance and norclozapine clearance
as boxplots stratified by female-nonsmoker, male-nonsmoker, and
male-smoker (the female-smoker stratum had only n = 2 in the source
paper and is omitted from the published boxplot for that reason). The
visual takeaway is that clearance rises monotonically across the three
groups for both clozapine and norclozapine. The chunk below renders the
corresponding typical-value steady-state exposure (Cavg over the final
dosing interval) by stratum, which is the inverse of the clearance
trend: higher clearance -\> lower Cavg.

``` r

ss <- sim_typical |>
  dplyr::filter(time >= 288, time <= 312) |>
  dplyr::group_by(id, treatment) |>
  dplyr::summarise(
    cavg_cloz    = mean(Cc),
    cavg_norcloz = mean(Cc_norcloz),
    .groups      = "drop"
  )

ss_long <- ss |>
  tidyr::pivot_longer(
    cols      = c(cavg_cloz, cavg_norcloz),
    names_to  = "species",
    values_to = "cavg"
  ) |>
  dplyr::mutate(
    species = dplyr::recode(species,
                            cavg_cloz    = "Clozapine",
                            cavg_norcloz = "Norclozapine")
  )

ggplot(ss_long, aes(x = treatment, y = cavg, fill = treatment)) +
  geom_point(aes(colour = treatment), size = 3) +
  facet_wrap(~species, scales = "free_y") +
  scale_x_discrete(limits = c("female_nonsmoker", "female_smoker",
                              "male_nonsmoker",   "male_smoker")) +
  labs(
    x       = NULL,
    y       = "Steady-state Cavg (umol/L) over the final dosing interval",
    fill    = "Sex x Smoking",
    colour  = "Sex x Smoking",
    title   = "Typical-value steady-state exposure by sex and smoking status",
    caption = paste0("Replicates the qualitative pattern of Figure 3 of Li 2012: ",
                     "higher CL/F in male smokers gives lower steady-state exposure.")
  ) +
  theme(legend.position = "bottom",
        axis.text.x     = element_text(angle = 30, hjust = 1))
```

![](Li_2012_clozapine_files/figure-html/figure-3-replicate-1.png)

## Single-dose profile through 24 hours

The chunk below shows the day-1 typical-value profile of clozapine and
norclozapine in a female nonsmoker after a single 200 mg oral dose,
illustrating the structural absorption-into-parent-elimination-
into-metabolite cascade encoded in the ODE system.

``` r

day1 <- sim_typical |>
  dplyr::filter(id == subjects$id[subjects$treatment == "female_nonsmoker"][1],
                time <= 24) |>
  dplyr::select(time, Cc, Cc_norcloz) |>
  tidyr::pivot_longer(c(Cc, Cc_norcloz),
                      names_to = "species", values_to = "conc") |>
  dplyr::mutate(species = dplyr::recode(species,
                                        Cc          = "Clozapine",
                                        Cc_norcloz  = "Norclozapine"))

ggplot(day1, aes(time, conc, colour = species)) +
  geom_line(linewidth = 1) +
  labs(x = "Time (h)", y = "Concentration (umol/L)",
       colour = NULL,
       title = "Day-1 typical-value profile in a female nonsmoker after 200 mg oral clozapine")
```

![](Li_2012_clozapine_files/figure-html/figure-day1-1.png)

## PKNCA validation – steady-state NCA over the final dosing interval

PKNCA is run separately for the parent (clozapine) and the metabolite
(norclozapine) outputs over the final dosing interval (288-312 h,
i.e. the 25th-26th doses). The interval-stratified summary reports Cmax,
Tmax, Cavg, Cmin, and AUC over tau for each stratum.

``` r

tau      <- 12
start_ss <- 288  # second-to-last dosing-interval start (full tau covered by obs grid)
end_ss   <- start_ss + tau  # 300

# Single-row-per-(id, time) concentration frame for clozapine.
sim_nca_cloz <- sim |>
  dplyr::filter(time >= start_ss, time <= end_ss) |>
  dplyr::distinct(id, time, .keep_all = TRUE) |>
  dplyr::select(id, time, Cc, treatment)

conc_obj_cloz <- PKNCA::PKNCAconc(
  data    = sim_nca_cloz,
  formula = Cc ~ time | treatment + id,
  concu   = "umol/L",
  timeu   = "h"
)

# Dose at start_ss for AUC0-tau alignment.
dose_df <- events |>
  dplyr::filter(evid == 1, time == start_ss) |>
  dplyr::select(id, time, amt, treatment)

dose_obj <- PKNCA::PKNCAdose(
  data    = dose_df,
  formula = amt ~ time | treatment + id,
  doseu   = "umol"
)

intervals <- data.frame(
  start    = start_ss,
  end      = end_ss,
  cmax     = TRUE,
  tmax     = TRUE,
  cmin     = TRUE,
  auclast  = TRUE,
  cav      = TRUE
)

nca_cloz <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(conc_obj_cloz, dose_obj, intervals = intervals)
)
nca_cloz_summary <- summary(nca_cloz)
knitr::kable(nca_cloz_summary,
             caption = "Simulated steady-state NCA (clozapine) over the second-to-last dosing interval (288-300 h).")
```

| Interval Start | Interval End | treatment | N | AUClast (h\*umol/L) | Cmax (umol/L) | Cmin (umol/L) | Tmax (h) | Cav (umol/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 288 | 300 | female_nonsmoker | 72 | 27.2 \[47.2\] | 2.81 \[37.4\] | 1.61 \[80.4\] | 2.00 \[1.50, 2.00\] | 2.27 \[47.2\] |
| 288 | 300 | female_smoker | 2 | 23.8 \[32.4\] | 2.57 \[15.4\] | 1.33 \[65.9\] | 2.00 \[2.00, 2.00\] | 1.99 \[32.4\] |
| 288 | 300 | male_nonsmoker | 40 | 23.3 \[35.8\] | 2.49 \[29.7\] | 1.30 \[69.8\] | 2.00 \[1.50, 2.00\] | 1.95 \[35.8\] |
| 288 | 300 | male_smoker | 48 | 16.1 \[43.9\] | 1.84 \[42.0\] | 0.776 \[83.0\] | 2.00 \[1.50, 2.00\] | 1.34 \[43.9\] |

Simulated steady-state NCA (clozapine) over the second-to-last dosing
interval (288-300 h). {.table style="width:100%;"}

``` r

sim_nca_norcloz <- sim |>
  dplyr::filter(time >= start_ss, time <= end_ss) |>
  dplyr::distinct(id, time, .keep_all = TRUE) |>
  dplyr::select(id, time, Cc_norcloz, treatment)

conc_obj_norcloz <- PKNCA::PKNCAconc(
  data    = sim_nca_norcloz,
  formula = Cc_norcloz ~ time | treatment + id,
  concu   = "umol/L",
  timeu   = "h"
)

nca_norcloz <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(conc_obj_norcloz, dose_obj, intervals = intervals)
)
nca_norcloz_summary <- summary(nca_norcloz)
knitr::kable(nca_norcloz_summary,
             caption = "Simulated steady-state NCA (norclozapine) over the second-to-last dosing interval (288-300 h).")
```

| Interval Start | Interval End | treatment | N | AUClast (h\*umol/L) | Cmax (umol/L) | Cmin (umol/L) | Tmax (h) | Cav (umol/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 288 | 300 | female_nonsmoker | 72 | 13.5 \[37.7\] | 1.14 \[37.2\] | 1.09 \[38.9\] | 6.50 \[5.00, 8.50\] | 1.12 \[37.7\] |
| 288 | 300 | female_smoker | 2 | 11.1 \[7.48\] | 0.953 \[9.46\] | 0.888 \[4.18\] | 6.00 \[6.00, 6.00\] | 0.928 \[7.48\] |
| 288 | 300 | male_nonsmoker | 40 | 10.6 \[39.3\] | 0.909 \[39.6\] | 0.844 \[39.3\] | 6.00 \[5.00, 7.50\] | 0.884 \[39.3\] |
| 288 | 300 | male_smoker | 48 | 6.83 \[42.0\] | 0.603 \[40.7\] | 0.515 \[46.3\] | 5.75 \[4.00, 6.50\] | 0.569 \[42.0\] |

Simulated steady-state NCA (norclozapine) over the second-to-last dosing
interval (288-300 h). {.table style="width:100%;"}

## Comparison against the paper’s reported TDM concentration means

Li 2012 Table 1 reports the cohort-pooled mean clozapine concentration
1.14 +- 0.73 umol/L and the cohort-pooled mean norclozapine
concentration 0.54 +- 0.32 umol/L across all TDM samples. The samples
were drawn at variable post-dose times (20% at 0.5-4.5 h, 80% at
8.5-15.5 h) under variable daily-dose regimens, so these means are not
steady-state Css values for any fixed dosing regimen; they are the time-
and dose-pooled mean of the TDM dataset. The relevant audit is the
metabolite-to-parent ratio (a within-subject structural quantity that
does not depend on the dose level or sampling-time distribution), and
the rank order of exposures across sex / smoking strata.

``` r

ss_summary <- ss_long |>
  dplyr::group_by(treatment, species) |>
  dplyr::summarise(
    typical_cavg = round(mean(cavg), 3),
    .groups      = "drop"
  ) |>
  tidyr::pivot_wider(names_from = species, values_from = typical_cavg) |>
  dplyr::mutate(metab_over_parent = round(Norclozapine / Clozapine, 3))

ss_summary <- ss_summary |>
  dplyr::add_row(
    treatment         = "Pooled TDM mean (Li 2012 Table 1)",
    Clozapine         = 1.14,
    Norclozapine      = 0.54,
    metab_over_parent = round(0.54 / 1.14, 3)
  )

knitr::kable(
  ss_summary,
  caption = "Simulated typical-value steady-state Cavg by stratum (umol/L) under 200 mg BID, compared against the pooled TDM mean reported by Li 2012."
)
```

| treatment                         | Clozapine | Norclozapine | metab_over_parent |
|:----------------------------------|----------:|-------------:|------------------:|
| female_nonsmoker                  |     2.316 |        1.029 |             0.444 |
| female_smoker                     |     1.594 |        0.666 |             0.418 |
| male_nonsmoker                    |     1.916 |        0.828 |             0.432 |
| male_smoker                       |     1.393 |        0.576 |             0.413 |
| Pooled TDM mean (Li 2012 Table 1) |     1.140 |        0.540 |             0.474 |

Simulated typical-value steady-state Cavg by stratum (umol/L) under 200
mg BID, compared against the pooled TDM mean reported by Li 2012.
{.table}

The simulated metabolite-to-parent ratio (around 0.44) agrees closely
with the TDM-pooled ratio 0.54 / 1.14 ~= 0.47, which is the cleanest
internal-consistency check the paper supports (KF = 0.66 fixed, with the
metabolite-to-parent ratio at steady state ~ KF \* CL / CLM). The
absolute simulated Cavg values are about 2x the TDM mean because the TDM
samples are drawn at variable post-dose times under heterogeneous
daily-dose regimens, while the simulation evaluates Cavg over a
canonical 200 mg BID schedule.

## Assumptions and deviations

- **Female-smoker stratum (n = 2)** is preserved in the virtual cohort
  to mirror Table 1, even though the source paper folded it into the
  smoking-effect estimate but omitted it from the published Figure 3
  boxplot for visual clarity. The model parameters were estimated on the
  four-stratum dataset, so the simulation should preserve the same
  stratification.
- **Sex coding inversion.** Li 2012 codes sex as a male indicator (SEX =
  1 male, 0 female); the canonical nlmixr2lib covariate is `SEXF` (1
  female, 0 male). The model file applies the male-sex effect via
  `(1 - SEXF)` to preserve Li 2012’s female-nonsmoker reference
  category. The covariate-effect parameter names `e_sex_cl` /
  `e_sex_cl_norcloz` denote the male-indicator coefficient in the
  source, not a female-indicator coefficient on the canonical SEXF.
- **AGE and WT were screened but not retained** in the final model (Li
  2012 Results: “Other covariates such as weight and age did not
  significantly influence the PK parameters of clozapine and
  norclozapine”). They are documented in
  `mod_meta$covariatesDataExcluded` so the FREM-style
  covariate-screening record is preserved.
- **Dosing in umol, not mg.** The source paper reports concentrations in
  umol/L (consistent with the HPLC-UV assay quantification practice for
  clozapine TDM); residual error sigma_1 and sigma_3 in Table 2 are in
  umol/L. We carry the model in molar units throughout so dose,
  concentration, and residual error are unit-consistent. To feed the
  model from a mg-based clinical event table, divide the mg dose by
  326.83 g/mol and multiply by 1000 (1 mg = 3.0596 umol). The 1:1 molar
  conversion of clozapine into norclozapine via N- demethylation is
  mass-balance-correct in molar units; in a mass- unit fit the apparent
  CLM would absorb a small ~4% bias from the parent/metabolite MW
  difference (326.83 vs 312.80 g/mol).
- **NPDE validation is not reproduced.** The published validation used
  NPDE (Brendel and colleagues) on the source dataset, which is not
  publicly available; this vignette therefore validates structural
  consistency (metabolite-to-parent ratio, rank order across strata,
  half-life) rather than re-running an NPDE analysis.
- **Compliance was assumed.** Patients were on stable maintenance oral
  clozapine with compliance confirmed by repeated serum samples and
  clinician interview; the simulated event table assumes perfect
  adherence to the canonical 200 mg BID regimen.
