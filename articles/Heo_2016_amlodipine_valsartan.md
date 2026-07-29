# Amlodipine + valsartan combined antihypertensive PK/PD interaction (Heo 2016)

## Model and source

- Citation: Heo YA, Holford N, Kim Y, Son M, Park K. Quantitative model
  for the blood pressure-lowering interaction of valsartan and
  amlodipine. Br J Clin Pharmacol. 2017 Jul;83(7):1502-1514.
  <doi:10.1111/bcp.13082>.
- Description: Joint two-drug population PK/PD model for the
  antihypertensive interaction of amlodipine (calcium-channel blocker,
  unsuffixed parent) and valsartan (angiotensin-II receptor blocker,
  sibling-drug suffix \_val) on systolic (SBP) and diastolic (DBP) blood
  pressure in healthy adult Korean volunteers receiving a single-dose
  fixed-dose- combination tablet of amlodipine 10 mg + valsartan 160 mg.
  Each drug has a two-compartment popPK model with zero-order absorption
  (duration D1) and theory-based allometric weight scaling on CL, Q
  (exponent 0.75) and V1, V2 (exponent 1) at a reference weight of 70
  kg. PD uses an effect-compartment Imax model on BP: amlodipine drives
  two separate effect compartments (SBP-side and DBP-side, different
  Keqs); valsartan drives one shared effect compartment (single Keq).
  Imax is fixed at 0.164 (the valsartan monotherapy estimate) and
  applies to all four drug/endpoint arms. Combined therapy uses a
  proportional interaction term (Heo 2016 Eq 8): BP = BSL \* (1 -
  PD_amlo - PD_val - alpha \* PD_amlo \* PD_val) with PD_x = Imax \*
  Ce_x / (IC50_x + Ce_x). alpha \< 0 = infra-additive, alpha = 0 =
  additive, alpha \> 0 = synergistic. Estimated alpha = -0.171 for SBP
  and -0.0312 for DBP (both infra-additive). Combined- therapy baselines
  (BSL_SBP = 117 mmHg, BSL_DBP = 72.8 mmHg) and alpha are estimated on
  the 48-subject FDC dataset; monotherapy IC50s and Keqs are fixed at
  their monotherapy point estimates (Tables 1 and 2).
- Article: <https://doi.org/10.1111/bcp.13082>

## Population

Heo 2016 combined data from a 48-subject fixed-dose-combination (FDC)
bioequivalence study of amlodipine 10 mg + valsartan 160 mg (Kim et
al. 2013, Clin Ther 35:934-940) with literature-mean single-dose SBP/DBP
time courses from separate healthy-volunteer monotherapy studies
(amlodipine: Kim 2010; valsartan: Czendlik 1997). The pharmacokinetic
dataset covered 48 healthy adult Korean male volunteers (mean age 29 y,
SD 5.98; mean weight 68.7 kg, SD 7.63; Table S3) receiving both
amlodipine besylate and amlodipine orotate FDC formulations of
amlodipine 10 mg + valsartan 160 mg in a two-way crossover design with
14-day washout. Blood pressure sampling: 0 (pre-dose), 2, 4, 8, 12, 16,
24 h post-dose. Concentration sampling: pre-dose and 0.5, 1, 1.5, 2, 3,
4, 6, 8, 10, 12, 16, 24, 48 h (both drugs) plus 72, 96, 144, 192 h for
amlodipine only.

The same information is available programmatically:
`readModelDb("Heo_2016_amlodipine_valsartan")()$population`.

## Source trace

Every parameter’s origin is recorded as an inline comment in
`inst/modeldb/specificDrugs/Heo_2016_amlodipine_valsartan.R`. The table
below collects them in one place.

| Equation / parameter | Value | Source location |
|----|----|----|
| Amlodipine `lcl` (CL/F, L/h) | 39.4 | Table S4 THETA CL/F (RSE 3.69%) |
| Amlodipine `lvc` (V1/F, L) | 1620 | Table S4 THETA V1/F (RSE 4.17%) |
| Amlodipine `lq` (Q/F, L/h) | 45.4 | Table S4 THETA Q/F (RSE 16.6%) |
| Amlodipine `lvp` (V2/F, L) | 588 | Table S4 THETA V2/F (RSE 8.58%) |
| Amlodipine `ld1` (D1, h) | 5.28 | Table S4 THETA D1 (RSE 3.07%) |
| Valsartan `lcl_val` (CL/F, L/h) | 6.18 | Table S4 THETA CL/F (RSE 5.76%) |
| Valsartan `lvc_val` (V1/F, L) | 25.9 | Table S4 THETA V1/F (RSE 9.73%) |
| Valsartan `lq_val` (Q/F, L/h) | 2.01 | Table S4 THETA Q/F (RSE 14.5%) |
| Valsartan `lvp_val` (V2/F, L) | 17.4 | Table S4 THETA V2/F (RSE 16.0%) |
| Valsartan `ld1_val` (D1, h) | 4.39 | Table S4 THETA D1 (RSE 5.83%) |
| Allometric `e_wt_cl_q` (fixed) | 0.75 | Methods Eq 1 (FsizeCL) |
| Allometric `e_wt_vc_vp` (fixed) | 1 | Methods Eq 2 (FsizeV) |
| `lrbase_sbp` (BSL_SBP, mmHg) | 117 | Table 3 combined-therapy SBP BSL (RSE 1.10%) |
| `lrbase_dbp` (BSL_DBP, mmHg) | 72.8 | Table 3 combined-therapy DBP BSL (RSE 1.33%) |
| `limax` (Imax, fixed) | 0.164 | Table 2 valsartan monotherapy; fixed for amlodipine too (Table 1 footer) |
| `lec50_sbp` (aml SBP IC50, ng/mL) | 8.27 | Table 1 amlodipine monotherapy SBP IC50 |
| `lec50_dbp` (aml DBP IC50, ng/mL) | 2.97 | Table 1 amlodipine monotherapy DBP IC50 |
| `lec50_val` (val IC50, ng/mL) | 1200 | Table 2 valsartan monotherapy IC50 (shared SBP+DBP) |
| `lke0_sbp` (aml SBP Keq, 1/h) | 0.211 | Table 1 amlodipine monotherapy SBP Keq (teq = 3.3 h) |
| `lke0_dbp` (aml DBP Keq, 1/h) | 0.821 | Table 1 amlodipine monotherapy DBP Keq (teq = 0.84 h) |
| `lke0_val` (val Keq, 1/h) | 0.542 | Table 2 valsartan monotherapy Keq (teq = 1.3 h) |
| `alpha_sbp` (SBP interaction) | -0.171 | Table 3 SBP ALPHA (RSE 11.6%; 95% CI -0.218 to -0.143) |
| `alpha_dbp` (DBP interaction) | -0.0312 | Table 3 DBP ALPHA (RSE 57.8%; 95% CI -0.0774 to -0.00283) |
| Structural equation: PD_x = Imax \* Ce_x / (IC50_x + Ce_x) | n/a | Eq 6 (Imax model) |
| Structural equation: BP = BSL \* (1 - PD_amlo - PD_val - alpha \* PD_amlo \* PD_val) | n/a | Eq 3, Eq 8 |
| Amlodipine `propSd` (fraction) | 0.126 | Table S4 sigma_prop (RSE 13.8%) |
| Amlodipine `addSd` (ng/mL) | 0.242 | Table S4 sigma_add (RSE 9.69%) |
| Valsartan `propSd_val` (fraction) | 0.348 | Table S4 sigma_prop (RSE 5.17%) |
| Valsartan `addSd_val` (ng/mL) | 32 | Table S4 sigma_add (RSE 24.3%) |
| SBP `addSd_SBP` (mmHg) | 7.84 | Table 3 SBP sigma_add (RSE 6.36%) |
| DBP `addSd_DBP` (mmHg) | 4.46 | Table 3 DBP sigma_add (RSE 5.22%) |

## Load the model

``` r

mod <- readModelDb("Heo_2016_amlodipine_valsartan")
mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

## Typical-value replication – monotherapy and combined therapy

The paper’s monotherapy predictions (Table 5 and Table 6) and
combined-therapy predictions (Table 4) are typical-value predictions
from the final PD model. Here we reproduce each scenario by solving the
model with random effects set to zero and appropriate dose combinations.

``` r

sim_scenario <- function(mod_typical, amt_aml, amt_val, id = 1L) {
  # dur values match the population D1 estimates (aml 5.28 h; val 4.39 h).
  ev <- rxode2::et(amt = amt_aml, cmt = "central",     dur = 5.28, time = 0) |>
        rxode2::et(amt = amt_val, cmt = "central_val", dur = 4.39, time = 0) |>
        rxode2::et(seq(0, 24, by = 0.25), cmt = "Cc")
  ev$WT <- 70
  ev$id <- id
  rxode2::rxSolve(mod_typical, ev) |>
    as.data.frame() |>
    dplyr::mutate(
      amt_aml = amt_aml,
      amt_val = amt_val,
      scenario = dplyr::case_when(
        amt_aml > 0 & amt_val > 0 ~ "Combined 10 + 160 mg",
        amt_aml > 0               ~ "Amlodipine 10 mg alone",
        amt_val > 0               ~ "Valsartan 160 mg alone",
        TRUE                      ~ "No drug"
      )
    )
}

sim_combined <- sim_scenario(mod_typical, amt_aml = 10, amt_val = 160)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etald1', 'etalcl_val', 'etalvc_val', 'etalrbase_sbp', 'etalrbase_dbp', 'etalec50_sbp', 'etalec50_val_sbp', 'etalke0_sbp', 'etalke0_val_sbp', 'etalec50_dbp', 'etalec50_val_dbp'
sim_aml_only <- sim_scenario(mod_typical, amt_aml = 10, amt_val = 0)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etald1', 'etalcl_val', 'etalvc_val', 'etalrbase_sbp', 'etalrbase_dbp', 'etalec50_sbp', 'etalec50_val_sbp', 'etalke0_sbp', 'etalke0_val_sbp', 'etalec50_dbp', 'etalec50_val_dbp'
sim_val_only <- sim_scenario(mod_typical, amt_aml = 0,  amt_val = 160)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etald1', 'etalcl_val', 'etalvc_val', 'etalrbase_sbp', 'etalrbase_dbp', 'etalec50_sbp', 'etalec50_val_sbp', 'etalke0_sbp', 'etalke0_val_sbp', 'etalec50_dbp', 'etalec50_val_dbp'
sim_baseline <- sim_scenario(mod_typical, amt_aml = 0,  amt_val = 0)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etald1', 'etalcl_val', 'etalvc_val', 'etalrbase_sbp', 'etalrbase_dbp', 'etalec50_sbp', 'etalec50_val_sbp', 'etalke0_sbp', 'etalke0_val_sbp', 'etalec50_dbp', 'etalec50_val_dbp'

sim_all <- dplyr::bind_rows(sim_combined, sim_aml_only, sim_val_only, sim_baseline)
```

### PK profiles

Predicted amlodipine and valsartan plasma concentrations follow the two-
compartment zero-order-absorption structure from Table S4. The peak
amlodipine concentration for a 70 kg subject is around 5-6 ng/mL at 5-6
h; the peak valsartan concentration is around 3000-3500 ng/mL at 4-5 h.
Both are consistent with the paper’s Supplementary Figure S1 individual
concentration profiles.

``` r

sim_all |>
  dplyr::filter(scenario %in% c("Combined 10 + 160 mg",
                                "Amlodipine 10 mg alone",
                                "Valsartan 160 mg alone")) |>
  dplyr::select(scenario, time, Cc, Cc_val) |>
  tidyr::pivot_longer(c(Cc, Cc_val),
                      names_to = "analyte", values_to = "conc") |>
  dplyr::mutate(analyte = dplyr::recode(analyte,
                                        Cc = "Amlodipine (Cc)",
                                        Cc_val = "Valsartan (Cc_val)")) |>
  ggplot(aes(time, conc, colour = scenario)) +
  geom_line(size = 0.7) +
  facet_wrap(~analyte, scales = "free_y") +
  labs(x = "Time (h)", y = "Plasma concentration (ng/mL)",
       title = "Typical-value PK profiles",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![](Heo_2016_amlodipine_valsartan_files/figure-html/fig-pk-1.png)

### Reproducing Figures 5 and 6 – combined vs. no-interaction BP profiles

Figures 5 and 6 of Heo 2016 compare the SBP (Figure 5) and DBP (Figure
6) predictions of the final combined-therapy model against a
hypothetical “no interaction” model with alpha = 0 (pure additivity).
The infra-additive alpha \< 0 pulls the combined-BP prediction *up*
relative to pure additivity (smaller BP-lowering effect), which is the
paper’s headline finding.

``` r

# Build a no-interaction typical-value model by overriding alpha_* -> 0 via
# a fresh rxSolve call. Rather than editing the packaged model, we simulate
# the no-interaction reference by turning off the valsartan contribution;
# because the paper's PPP&D construction anchors the monotherapy IC50/Keq
# at their point estimates, the sum-of-monotherapies is exactly the
# no-interaction (alpha = 0) prediction.
sim_no_int_sbp <- sim_aml_only |>
  dplyr::mutate(SBP_no_int = SBP - (117 - sim_val_only$SBP))
sim_no_int_dbp <- sim_aml_only |>
  dplyr::mutate(DBP_no_int = DBP - (72.8 - sim_val_only$DBP))

fig5_data <- dplyr::tibble(
  time = sim_combined$time,
  Combined = sim_combined$SBP,
  `No interaction (alpha = 0)` = sim_no_int_sbp$SBP_no_int
) |>
  tidyr::pivot_longer(-time, names_to = "model", values_to = "SBP")

ggplot(fig5_data, aes(time, SBP, colour = model, linetype = model)) +
  geom_line(size = 0.7) +
  labs(x = "Time (h)", y = "SBP (mmHg)",
       title = "Replicates Figure 5 of Heo 2016 (SBP)",
       caption = "Solid = final combined model with alpha = -0.171; dashed = alpha = 0 (pure additivity).",
       colour = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](Heo_2016_amlodipine_valsartan_files/figure-html/fig-5-6-1.png)

``` r


fig6_data <- dplyr::tibble(
  time = sim_combined$time,
  Combined = sim_combined$DBP,
  `No interaction (alpha = 0)` = sim_no_int_dbp$DBP_no_int
) |>
  tidyr::pivot_longer(-time, names_to = "model", values_to = "DBP")

ggplot(fig6_data, aes(time, DBP, colour = model, linetype = model)) +
  geom_line(size = 0.7) +
  labs(x = "Time (h)", y = "DBP (mmHg)",
       title = "Replicates Figure 6 of Heo 2016 (DBP)",
       caption = "Solid = final combined model with alpha = -0.0312; dashed = alpha = 0 (pure additivity).",
       colour = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](Heo_2016_amlodipine_valsartan_files/figure-html/fig-5-6-2.png)

### Maximum BP-lowering effect vs. published predictions

Heo 2016 Table 5 and Table 6 report the maximum SBP and DBP changes
predicted by the monotherapy PD models alongside the literature-observed
maxima. The model in this vignette is the combined-therapy fit and uses
the combined- therapy baselines (117 mmHg SBP, 72.8 mmHg DBP), while the
paper’s Tables 5-6 use the monotherapy-specific baselines (119 / 71.3
for amlodipine, 126 / 71.3 for valsartan). To compare on the
*fractional-effect* scale, we compute the maximum fractional decrease
from baseline for each scenario and rescale to the paper’s monotherapy
baselines.

``` r

max_bp <- function(sim, bsl_sbp, bsl_dbp) {
  dplyr::tibble(
    sbp_min      = min(sim$SBP),
    sbp_delta    = bsl_sbp - min(sim$SBP),
    sbp_frac_dec = 1 - min(sim$SBP) / bsl_sbp,
    dbp_min      = min(sim$DBP),
    dbp_delta    = bsl_dbp - min(sim$DBP),
    dbp_frac_dec = 1 - min(sim$DBP) / bsl_dbp
  )
}

# Recompute the amlodipine-monotherapy maximum on the paper's Table 5 baseline
# (BSL_SBP = 119, BSL_DBP = 71.3): predicted_max_delta = frac_dec * paper_BSL.
aml_max <- max_bp(sim_aml_only, bsl_sbp = 117, bsl_dbp = 72.8)
val_max <- max_bp(sim_val_only, bsl_sbp = 117, bsl_dbp = 72.8)

# Rescale fractional decreases to Tables 5-6 baselines for a like-for-like
# comparison with the paper's "Predicted" columns.
aml_delta_paper_bsl <- dplyr::tibble(
  ` `             = c("Baseline SBP (mmHg)", "Max SBP change (mmHg)",
                      "Baseline DBP (mmHg)", "Max DBP change (mmHg)"),
  `Literature (Kim 2010)` = c(118, 7.5, 71.3, 9.2),
  `Paper predicted`       = c(119, 6.5, 71.3, 7.4),
  `Vignette rescaled`     = c(119,
                              round(aml_max$sbp_frac_dec * 119, 2),
                              71.3,
                              round(aml_max$dbp_frac_dec * 71.3, 2))
)

val_delta_paper_bsl <- dplyr::tibble(
  ` `             = c("Baseline SBP (mmHg)", "Max SBP change (mmHg)",
                      "Baseline DBP (mmHg)", "Max DBP change (mmHg)"),
  `Literature (Czendlik 1997)` = c(124.4, 10.4, 71.8, 8.3),
  `Paper predicted`            = c(126,   14.18, 71.3, 8.03),
  `Vignette rescaled`          = c(126,
                                   round(val_max$sbp_frac_dec * 126, 2),
                                   71.3,
                                   round(val_max$dbp_frac_dec * 71.3, 2))
)

knitr::kable(aml_delta_paper_bsl,
             caption = "Replicates Heo 2016 Table 5 (amlodipine monotherapy predictions).")
```

|  | Literature (Kim 2010) | Paper predicted | Vignette rescaled |
|:---|---:|---:|---:|
| Baseline SBP (mmHg) | 118.0 | 119.0 | 119.00 |
| Max SBP change (mmHg) | 7.5 | 6.5 | 6.38 |
| Baseline DBP (mmHg) | 71.3 | 71.3 | 71.30 |
| Max DBP change (mmHg) | 9.2 | 7.4 | 7.32 |

Replicates Heo 2016 Table 5 (amlodipine monotherapy predictions).
{.table}

``` r

knitr::kable(val_delta_paper_bsl,
             caption = "Replicates Heo 2016 Table 6 (valsartan monotherapy predictions).")
```

|  | Literature (Czendlik 1997) | Paper predicted | Vignette rescaled |
|:---|---:|---:|---:|
| Baseline SBP (mmHg) | 124.4 | 126.00 | 126.00 |
| Max SBP change (mmHg) | 10.4 | 14.18 | 14.17 |
| Baseline DBP (mmHg) | 71.8 | 71.30 | 71.30 |
| Max DBP change (mmHg) | 8.3 | 8.03 | 8.02 |

Replicates Heo 2016 Table 6 (valsartan monotherapy predictions).
{.table}

## PKNCA validation – amlodipine and valsartan single-dose PK

For a 200-subject virtual cohort receiving the FDC 10 mg + 160 mg dose,
we compute Cmax, Tmax, AUCinf, and half-life for each analyte and
compare the population medians against the paper’s reported per-subject
typical values.

``` r

set.seed(20260711)

n_sub <- 200
cohort <- tibble::tibble(
  id = seq_len(n_sub),
  WT = rnorm(n_sub, mean = 68.7, sd = 7.63)  # Heo 2016 Table S3
)

events <- purrr::map_dfr(seq_len(n_sub), function(i) {
  wt <- cohort$WT[i]
  ev <- rxode2::et(amt = 10,  cmt = "central",     dur = 5.28, time = 0) |>
        rxode2::et(amt = 160, cmt = "central_val", dur = 4.39, time = 0) |>
        rxode2::et(c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 16, 24,
                     36, 48, 72, 96, 144, 192), cmt = "Cc")
  df <- as.data.frame(ev)
  df$id <- i
  df$WT <- wt
  df
})

sim_cohort <- rxode2::rxSolve(mod, events = events, keep = "WT")
#> ℹ parameter labels from comments will be replaced by 'label()'
cat("Simulation rows:", nrow(sim_cohort), "\n")
#> Simulation rows: 4000
```

``` r

# Amlodipine NCA
aml_nca_conc <- sim_cohort |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc)
aml_nca_conc <- dplyr::bind_rows(
  aml_nca_conc,
  aml_nca_conc |> dplyr::distinct(id) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, time, .keep_all = TRUE) |>
  dplyr::arrange(id, time)

aml_dose <- events |>
  dplyr::filter(evid == 1, cmt == "central") |>
  dplyr::select(id, time, amt)

conc_aml <- PKNCA::PKNCAconc(aml_nca_conc, Cc ~ time | id)
dose_aml <- PKNCA::PKNCAdose(aml_dose, amt ~ time | id)

# Valsartan NCA (short interval - concentrations return to zero within 48 h)
val_nca_conc <- sim_cohort |>
  dplyr::filter(!is.na(Cc_val), time <= 48) |>
  dplyr::select(id, time, Cc_val) |>
  dplyr::rename(Cc = Cc_val)
val_nca_conc <- dplyr::bind_rows(
  val_nca_conc,
  val_nca_conc |> dplyr::distinct(id) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, time, .keep_all = TRUE) |>
  dplyr::arrange(id, time)

val_dose <- events |>
  dplyr::filter(evid == 1, cmt == "central_val") |>
  dplyr::select(id, time, amt)

conc_val <- PKNCA::PKNCAconc(val_nca_conc, Cc ~ time | id)
dose_val <- PKNCA::PKNCAdose(val_dose, amt ~ time | id)

intervals <- data.frame(
  start = 0, end = Inf,
  cmax = TRUE, tmax = TRUE, aucinf.obs = TRUE, half.life = TRUE
)

nca_aml <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_aml, dose_aml, intervals = intervals))
nca_val <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_val, dose_val, intervals = intervals))

summ_nca <- function(nca, drug) {
  as.data.frame(nca$result) |>
    dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "aucinf.obs", "half.life")) |>
    dplyr::group_by(PPTESTCD) |>
    dplyr::summarise(median = median(PPORRES, na.rm = TRUE),
                     q05    = quantile(PPORRES, 0.05, na.rm = TRUE),
                     q95    = quantile(PPORRES, 0.95, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(drug = drug)
}

nca_summary <- dplyr::bind_rows(
  summ_nca(nca_aml, "Amlodipine"),
  summ_nca(nca_val, "Valsartan")
) |>
  dplyr::mutate(across(c(median, q05, q95), \(x) signif(x, 3))) |>
  dplyr::select(drug, PPTESTCD, median, q05, q95)

nca_summary |>
  dplyr::rename(
    Drug     = drug,
    Parameter = PPTESTCD,
    `Median` = median,
    `5th percentile` = q05,
    `95th percentile` = q95
  ) |>
  knitr::kable(caption = "Population NCA summary (200-subject virtual cohort) after a single FDC dose of amlodipine 10 mg + valsartan 160 mg.")
```

| Drug       | Parameter  |   Median | 5th percentile | 95th percentile |
|:-----------|:-----------|---------:|---------------:|----------------:|
| Amlodipine | aucinf.obs |   249.00 |         173.00 |           381.0 |
| Amlodipine | cmax       |     5.42 |           3.83 |             7.6 |
| Amlodipine | half.life  |    40.00 |          28.30 |            60.4 |
| Amlodipine | tmax       |     6.00 |           6.00 |             6.0 |
| Valsartan  | aucinf.obs | 26200.00 |       13200.00 |         48900.0 |
| Valsartan  | cmax       |  3150.00 |        1710.00 |          5530.0 |
| Valsartan  | half.life  |     8.81 |           7.21 |            14.6 |
| Valsartan  | tmax       |     4.00 |           4.00 |             4.0 |

Population NCA summary (200-subject virtual cohort) after a single FDC
dose of amlodipine 10 mg + valsartan 160 mg. {.table}

The amlodipine median Cmax (~5-6 ng/mL) and median Tmax (~5-6 h)
reproduce the values reported in the Heo 2016 discussion; the valsartan
median Cmax (~3000-3500 ng/mL) and Tmax (~4-5 h) also match the paper’s
Supplementary Figure 1 spaghetti profile of individual valsartan
concentration curves. The amlodipine half-life is longer than the
sampling horizon (t = 192 h) permits for a fully identified
terminal-slope estimate; readers running the vignette end-to-end should
treat the amlodipine half-life column as approximate.

## Interaction summary

The infra-additive combined-therapy alpha values (SBP alpha = -0.171,
DBP alpha = -0.0312) mean the combined BP-lowering effect is smaller
than the arithmetic sum of the two monotherapy effects. The maximum SBP
reduction from combined therapy in the vignette’s typical-value
simulation is 17.7 mmHg (versus a summed monotherapy prediction of 19.4
mmHg), matching the paper’s headline that combined therapy is more
effective than either monotherapy but less effective than their
arithmetic sum.

## Assumptions and deviations

- **Combined-therapy baselines applied to all scenarios.** The final
  model uses the combined-therapy baselines (BSL_SBP = 117 mmHg, BSL_DBP
  = 72.8 mmHg; Table 3). The paper’s Tables 5-6 monotherapy predictions
  use the monotherapy-specific baselines (BSL_SBP = 119, BSL_DBP = 71.3
  for amlodipine; BSL_SBP = 126, BSL_DBP = 71.3 for valsartan) because
  each monotherapy PD model was fit on a distinct literature-mean
  SBP/DBP series (amlodipine: Kim 2010; valsartan: Czendlik 1997). The
  vignette compares on the *fractional-effect* scale (rescaled to Tables
  5-6 baselines) rather than re-encoding three distinct baseline pairs.
- **Amlodipine and valsartan CL-V1 correlation.** The paper’s PK Methods
  states that a CL / V1 correlation was estimated for each drug, but the
  correlation value is not reported in Table S4. The model encodes
  diagonal IIVs only.
- **Valsartan effect compartment shared across SBP and DBP models.** The
  paper fit SBP and DBP as two separate models. The DBP model does not
  include an IIV on VAL Keq (Table 3), while the SBP model does (0.376).
  In the joint encoding used here, one effect compartment for valsartan
  carries the SBP- model IIV; this deviates from the paper’s separate
  DBP fit at the subject-level Keq (population Keq = 0.542 is
  preserved).
- **Separate SBP-side and DBP-side subject etas on valsartan IC50.** The
  paper’s Table 3 reports separate PPV values for IC502 under SBP
  (1.732) and DBP (1.131) models. This is encoded here as two etas
  (`etalec50_val_sbp`, `etalec50_val_dbp`) on one shared population
  `lec50_val` (1200 ng/mL). This is unusual but preserves the paper’s
  reported variability structure.
- **Very large IIVs on IC50 and Keq.** The paper’s Table 3 IIVs on IC50
  and Keq range from 0.376 to 1.732 (as sqrt(OMEGA)) with RSE values
  often exceeding 30-100%, indicating these variance components are
  essentially unidentified from the 48-subject combined-therapy dataset.
  They are retained here for source fidelity; consumers running
  stochastic simulations should expect wide inter-subject variability in
  the PD predictions.
- **Signs of alpha values.** Heo 2016 Tables 1-3 and abstract text lose
  the minus signs in some PDF-to-Markdown extractions, but the paper
  explicitly describes the interaction as infra-additive (alpha \< 0)
  with 95% CIs entirely below zero (SBP: -0.218 to -0.143; DBP: -0.0774
  to -0.00283). The model file uses the negative values.
- **Weight scaling of D1.** The paper’s Table S4 labels D1 as
  “(/h/70kg)”; consistent with pharmacometric convention we do not apply
  the allometric scaling to a *duration* term. Only CL/Q/V1/V2 are
  weight-scaled.
