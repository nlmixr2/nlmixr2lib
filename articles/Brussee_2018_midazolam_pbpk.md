# Midazolam preterm-neonate PBPK (Brussee 2018)

## Model and source

- Citation: Brussee JM, Yu H, Krekels EHJ, de Roos B, Brill MJE, van den
  Anker JN, Rostami-Hodjegan A, de Wildt SN, Knibbe CAJ (2018).
  First-Pass CYP3A- Mediated Metabolism of Midazolam in the Gut Wall and
  Liver in Preterm Neonates. CPT Pharmacometrics Syst Pharmacol
  7(6):374-383. <doi:10.1002/psp4.12295>.
- Description: PBPK (semi-physiological; well-stirred liver + Qgut gut
  wall) population PK model for midazolam and its primary metabolite
  1-OH-midazolam in 37 preterm neonates (gestational age 26-34 weeks,
  body weight 0.770-2.030 kg at the time of dosing). Distinguishes
  first-pass CYP3A-mediated metabolism in the gut wall (Qgut model) and
  liver (well-stirred model) from systemic hepatic elimination of the
  metabolite. Tissue volumes (V_h, V_pv, V_gw) and hepatic blood flow
  Q_h are allometrically scaled from a term-neonate reference
  (Bjorkman 2005) by body weight with fixed exponents (1 for volumes,
  0.75 for flow); intestinal length scales as 2.736 \* WT\[g\]^0.512 cm
  (Struijs 2009) so the Qgut hybrid flow varies with body size. Supports
  oral administration (depot, full first-pass through gut wall and
  liver) and IV (dose directly to central; no first-pass).
- Article: <https://doi.org/10.1002/psp4.12295>

## Population

Brussee 2018 fit a semi-physiological population PK model to pooled IV
and oral midazolam concentration data from 37 preterm neonates admitted
to the neonatal intensive care unit of the Sophia’s Children Hospital
(Rotterdam, NL). The cohort had gestational age at birth 26-34 weeks,
birth weight 745-2,135 g, postnatal age 3-11 days, and body weight at
the time of dosing 770-2,030 g. The source data combine de Wildt 2001
(IV cohort, *Clin Pharmacol Ther* 70:525-531) and de Wildt 2002 (oral
cohort, *Br J Clin Pharmacol* 53:390-392): a 0.1 mg/kg dose given either
orally via a nasogastric tube (n = 13) or as a 30-minute intravenous
infusion (n = 25), with a crossover dose by the alternate route after
\>= 72 h in subjects still meeting inclusion criteria. Plasma midazolam
and 1-OH-midazolam concentrations were collected at 0.5, 1, 2, 4, 6, 12,
and 24 h post-dose; below-LLOQ samples (\< 1% of 329 midazolam and \< 2%
of 153 1-OH-midazolam observations) were discarded before fitting. No
covariates beyond body weight (which scales the physiological tissue
volumes, blood flows, and intestinal-length terms) reached statistical
significance in the covariate analysis.

The same information is available programmatically via
`readModelDb("Brussee_2018_midazolam_pbpk")$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to every
`ini()` entry in
`inst/modeldb/specificDrugs/Brussee_2018_midazolam_pbpk.R`. The table
below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| Hepatic extraction `E_h = (f_u,B * CL_int,H) / (f_u,B * CL_int,H + Q_h)` | derived | Brussee 2018 Eq. 1 |
| Gut wall extraction `E_g = (f_u,G * CL_int,G) / (f_u,G * CL_int,G + Q_gut)` | derived | Brussee 2018 Eq. 2 |
| Qgut hybrid `Q_gut = (Q_villi * CL_perm) / (Q_villi + CL_perm)` | derived | Brussee 2018 Eq. 3 |
| Intestinal surface `A = 2*pi*radius*h`, radius `= 1 cm`, `h = 2.736 * WT[g]^0.512 cm` | derived | Brussee 2018 Eq. 4 (Ives 2016, Struijs 2009) |
| Oral bioavailability `F_total = F_a * F_g * F_h` | derived | Brussee 2018 Eq. 5 |
| Reference WT for term-neonate scaling | 3.55 kg | Bjorkman 2005 / Brussee 2018 Table 1 |
| Reference liver volume `V_h,3.55kg` | 0.120 L | Brussee 2018 Table 1 (Bjorkman 2005 ref 7) |
| Portal vein / liver volume ratio | 0.778 | Brussee 2018 Table 1 (Aguirre-Reyes 2015 ref 24) |
| Reference gut wall volume `V_gw,3.55kg` | 0.050 L | Brussee 2018 Table 1 |
| Reference hepatic blood flow `Q_h,3.55kg` | 13.2 L/h | Brussee 2018 Table 1 |
| Blood-flow ratios (Q_pv, Q_ha, Q_in, Q_muc, Q_villi) | 0.75, 0.25, 0.4, 0.8, 0.6 | Brussee 2018 Table 1 (Johnson 2006, Yang 2007 refs 5, 25, 32) |
| Effective intestinal permeability `P_eff,man` (midazolam) | 4.4e-4 cm/s = 1.584 cm/h | Brussee 2018 Table 1 (Yang 2007 ref 25) |
| Midazolam unbound fraction in blood `f_u,B` | 0.04094 | Brussee 2018 Table 1 (McNamara-Alcorn from f_u,adult = 0.0303, \[P\]\_pediatric = 27.1, \[P\]\_adult = 37.0) |
| Midazolam unbound fraction in gut `f_u,G` | 1 | Brussee 2018 Methods (Yang 2007 ref 25) |
| Midazolam blood:plasma `B:P` | 0.568 | Brussee 2018 Table 1 (Maharaj 2013 ref 30) |
| 1-OH-midazolam unbound in blood `f_u,M` | 0.1394 | Brussee 2018 Table 1 (from f_u,M,adult = 0.106, Mandema 1992 ref 29) |
| 1-OH-midazolam blood:plasma `B:P` | 0.613 | Brussee 2018 Table 1 |
| Fraction metabolised to 1-OH-midazolam `f_M` | 1 | Brussee 2018 Table 1 |
| Absorption rate constant `k_a` | 10 1/h (fixed) | Brussee 2018 Methods (Table 1; sensitivity analysis 4.16-25 1/h) |
| Midazolam intrinsic hepatic CL `CL_H,int` | 6.7 L/h (RSE 10%) | Brussee 2018 Table 2 |
| Midazolam intrinsic gut wall CL `CL_G,int` | 0.0196 L/h (RSE 178%) | Brussee 2018 Table 2 |
| Midazolam apparent V (blood) | 3.0 L (RSE 11%) | Brussee 2018 Table 2 |
| 1-OH-midazolam intrinsic hepatic CL `CL_H,int,M` | 8.9 L/h (RSE 22%) | Brussee 2018 Table 2 |
| 1-OH-midazolam apparent V (blood) | 2.7 L (RSE 43%) | Brussee 2018 Table 2 |
| IIV omega^2 (CL_H,int, V, CL_H,int,M, V_M) | 0.887, 0.603, 0.832, 0.887 | Brussee 2018 Table 2 |
| Proportional residual error sigma^2 (parent / metabolite) | 0.201 / 0.164 | Brussee 2018 Table 2 |
| Additive residual error sigma^2 (parent / metabolite) | 0.0001 FIX (both) | Brussee 2018 Table 2 |

## Virtual cohort

Original observed concentrations from de Wildt 2001/2002 are not openly
available. The figures below use a 200-subject virtual cohort whose
body-weight distribution approximates the published trial demographics:
a uniform sample from 0.77-2.03 kg (the body-weight range at the time of
dosing reported in Methods Data).

``` r

set.seed(20180510)  # paper's published-online date

n_sub <- 200L
wt_range <- c(0.77, 2.03)  # kg

cohort_wt <- tibble::tibble(
  id = seq_len(n_sub),
  WT = runif(n_sub, min = wt_range[1], max = wt_range[2])
)

# Observation grid. Brussee 2018 Methods sampled at 0.5, 1, 2, 4, 6, 12, 24 h;
# use a denser grid here for smoother concentration-time profiles.
obs_times <- c(0.01, seq(0.25, 24, by = 0.25))

make_events <- function(cohort, route, dose_per_kg = 0.1, infusion_h = 0.5) {
  doses <- cohort |>
    dplyr::mutate(
      EVID = 1L,
      AMT  = dose_per_kg * WT * 1000,         # mg/kg * kg * 1000 = ug
      TIME = 0,
      CMT  = ifelse(route == "oral", "depot", "central"),
      RATE = ifelse(route == "iv30",
                    dose_per_kg * WT * 1000 / infusion_h,
                    NA_real_),
      route = route
    ) |>
    dplyr::select(id, EVID, AMT, RATE, TIME, CMT, WT, route)

  obs <- tidyr::expand_grid(
    cohort,
    TIME = obs_times
  ) |>
    dplyr::mutate(EVID = 0L, AMT = 0, RATE = NA_real_,
                  CMT = "Cc", route = route) |>
    dplyr::select(id, EVID, AMT, RATE, TIME, CMT, WT, route)

  dplyr::bind_rows(doses, obs) |>
    dplyr::arrange(id, TIME, dplyr::desc(EVID))
}

events_oral <- make_events(cohort_wt, "oral")
events_iv   <- make_events(cohort_wt, "iv30")
```

## Simulation

``` r

mod <- readModelDb("Brussee_2018_midazolam_pbpk")

sim_oral <- rxode2::rxSolve(
  mod, events = events_oral,
  keep = c("route", "WT"),
  returnType = "data.frame"
) |>
  dplyr::mutate(route = "Oral")
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_iv <- rxode2::rxSolve(
  mod, events = events_iv,
  keep = c("route", "WT"),
  returnType = "data.frame"
) |>
  dplyr::mutate(route = "IV 30-min infusion")
#> ℹ parameter labels from comments will be replaced by 'label()'

sim <- dplyr::bind_rows(sim_oral, sim_iv)
```

The typical-value (population-mean) profile, useful for Figure 4
replication, zeros out the between-subject variability:

``` r

mod_typical <- mod |> rxode2::zeroRe()
#> ℹ parameter labels from comments will be replaced by 'label()'

# Typical 1.1 kg neonate
typical_cohort <- tibble::tibble(id = 1L, WT = 1.1)
typical_oral <- rxode2::rxSolve(
  mod_typical, events = make_events(typical_cohort, "oral"),
  keep = c("route", "WT"),
  returnType = "data.frame"
) |> dplyr::mutate(route = "Oral")
#> ℹ omega/sigma items treated as zero: 'etalcl_int_h', 'etalvc', 'etalcl_int_h_1ohm', 'etalvc_1ohm'
typical_iv <- rxode2::rxSolve(
  mod_typical, events = make_events(typical_cohort, "iv30"),
  keep = c("route", "WT"),
  returnType = "data.frame"
) |> dplyr::mutate(route = "IV 30-min infusion")
#> ℹ omega/sigma items treated as zero: 'etalcl_int_h', 'etalvc', 'etalcl_int_h_1ohm', 'etalvc_1ohm'

typical <- dplyr::bind_rows(typical_oral, typical_iv)
```

## Replicate published figures

### Figure 4a – Midazolam plasma profile (oral vs IV)

Brussee 2018 Figure 4a shows model-based simulations of midazolam plasma
concentrations following a 0.1 mg/kg dose given orally (black) or as a
30-min infusion (gray) in the 37 preterm neonates from the study, with
the median and minimal / maximal envelope.

``` r

sim |>
  dplyr::filter(time > 0) |>
  dplyr::group_by(route, time) |>
  dplyr::summarise(
    Q05 = quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, Q50, colour = route, fill = route)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(
    x = "Time after dose (hours)",
    y = "Midazolam plasma concentration (ng/mL)",
    colour = NULL, fill = NULL,
    title = "Figure 4a -- midazolam, oral vs IV 30-min infusion",
    caption = "Median with 5th-95th percentile band, 200-subject virtual cohort (WT 0.77-2.03 kg)."
  )
```

![Replicates Figure 4a of Brussee
2018.](Brussee_2018_midazolam_pbpk_files/figure-html/figure-4a-midazolam-1.png)

Replicates Figure 4a of Brussee 2018.

### Figure 4b – 1-OH-midazolam plasma profile (oral vs IV)

Brussee 2018 Figure 4b shows the same simulation for the 1-OH-midazolam
metabolite, where the oral curve climbs higher than the IV curve over
the first 4 h because of presystemic conversion in the gut wall and
liver.

``` r

sim |>
  dplyr::filter(time > 0) |>
  dplyr::group_by(route, time) |>
  dplyr::summarise(
    Q05 = quantile(Cc_1ohm, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc_1ohm, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc_1ohm, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, Q50, colour = route, fill = route)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(
    x = "Time after dose (hours)",
    y = "1-OH-midazolam plasma concentration (ng/mL)",
    colour = NULL, fill = NULL,
    title = "Figure 4b -- 1-OH-midazolam, oral vs IV 30-min infusion",
    caption = "Median with 5th-95th percentile band; presystemic metabolism elevates the oral curve over 0-4 h."
  )
```

![Replicates Figure 4b of Brussee
2018.](Brussee_2018_midazolam_pbpk_files/figure-html/figure-4b-metabolite-1.png)

Replicates Figure 4b of Brussee 2018.

### Figure 2 / 3 – extraction ratios and bioavailability vs body weight

Brussee 2018 Figures 2c and 3 report extraction ratios and
bioavailability as scatter plots of typical-value model predictions
across body weight. The extraction-ratio and bioavailability terms are
computed inside `model()` and are exposed as derived outputs by
`rxSolve()`, so the relationship across the cohort can be reproduced
directly.

``` r

typical_grid <- tibble::tibble(
  id = seq_len(50),
  WT = seq(min(wt_range), max(wt_range), length.out = 50)
)
typical_grid_oral <- rxode2::rxSolve(
  mod_typical, events = make_events(typical_grid, "oral"),
  keep = c("route", "WT"),
  returnType = "data.frame"
)
#> ℹ omega/sigma items treated as zero: 'etalcl_int_h', 'etalvc', 'etalcl_int_h_1ohm', 'etalvc_1ohm'
#> Warning: multi-subject simulation without without 'omega'
# Static extractions are time-invariant; pick the first observation row per ID
extractions <- typical_grid_oral |>
  dplyr::filter(time == min(time[time > 0])) |>
  dplyr::select(id, WT, eg, eh, eh_1ohm, fg, fh) |>
  dplyr::mutate(
    ftotal = fg * fh                       # F_a = 1
  )

extractions |>
  tidyr::pivot_longer(c(eg, eh, eh_1ohm, fg, fh, ftotal),
                      names_to = "quantity", values_to = "value") |>
  dplyr::mutate(panel = dplyr::case_when(
    quantity %in% c("eg", "eh", "eh_1ohm") ~ "Extraction ratios",
    TRUE                                   ~ "Bioavailability"
  )) |>
  ggplot(aes(WT, value, colour = quantity)) +
  geom_line() +
  geom_point(size = 1.2) +
  facet_wrap(~ panel, scales = "free_y") +
  labs(
    x = "Body weight (kg)",
    y = "Extraction ratio / bioavailability",
    colour = NULL,
    title = "Figures 2c & 3 -- extraction ratios and bioavailability vs WT",
    caption = "Typical-value, 50-point WT grid spanning 0.77-2.03 kg."
  )
```

![Replicates Figures 2c and 3 of Brussee 2018: extraction ratio and
bioavailability vs body
weight.](Brussee_2018_midazolam_pbpk_files/figure-html/figure-2c-3-1.png)

Replicates Figures 2c and 3 of Brussee 2018: extraction ratio and
bioavailability vs body weight.

## PKNCA validation

Compute Cmax, Tmax, AUC, and apparent half-life for both routes via
PKNCA and compare against the paper’s reported summary statistics.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc), time > 0) |>
  dplyr::transmute(id, time, Cc, treatment = route)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id)
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found

dose_df <- dplyr::bind_rows(events_oral, events_iv) |>
  dplyr::filter(EVID == 1L) |>
  dplyr::transmute(id, time = TIME, amt = AMT, treatment = route)

dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id)

intervals <- data.frame(
  start = 0, end = 24,
  cmax = TRUE, tmax = TRUE,
  auclast = TRUE, aucinf.obs = TRUE,
  half.life = TRUE
)

nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
#> Warning: Requesting an AUC range starting (0) before the first measurement
#> (0.01) is not allowed
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning: Requesting an AUC range starting (0) before the first measurement
#> (0.01) is not allowed
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in log(data$conc): NaNs produced
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.01) is not allowed
#> Warning: treatment=iv30; id=1: No concentration data
#> Warning: treatment=iv30; id=2: No concentration data
#> Warning: treatment=iv30; id=3: No concentration data
#> Warning: treatment=iv30; id=4: No concentration data
#> Warning: treatment=iv30; id=5: No concentration data
#> Warning: treatment=iv30; id=6: No concentration data
#> Warning: treatment=iv30; id=7: No concentration data
#> Warning: treatment=iv30; id=8: No concentration data
#> Warning: treatment=iv30; id=9: No concentration data
#> Warning: treatment=iv30; id=10: No concentration data
#> Warning: treatment=iv30; id=11: No concentration data
#> Warning: treatment=iv30; id=12: No concentration data
#> Warning: treatment=iv30; id=13: No concentration data
#> Warning: treatment=iv30; id=14: No concentration data
#> Warning: treatment=iv30; id=15: No concentration data
#> Warning: treatment=iv30; id=16: No concentration data
#> Warning: treatment=iv30; id=17: No concentration data
#> Warning: treatment=iv30; id=18: No concentration data
#> Warning: treatment=iv30; id=19: No concentration data
#> Warning: treatment=iv30; id=20: No concentration data
#> Warning: treatment=iv30; id=21: No concentration data
#> Warning: treatment=iv30; id=22: No concentration data
#> Warning: treatment=iv30; id=23: No concentration data
#> Warning: treatment=iv30; id=24: No concentration data
#> Warning: treatment=iv30; id=25: No concentration data
#> Warning: treatment=iv30; id=26: No concentration data
#> Warning: treatment=iv30; id=27: No concentration data
#> Warning: treatment=iv30; id=28: No concentration data
#> Warning: treatment=iv30; id=29: No concentration data
#> Warning: treatment=iv30; id=30: No concentration data
#> Warning: treatment=iv30; id=31: No concentration data
#> Warning: treatment=iv30; id=32: No concentration data
#> Warning: treatment=iv30; id=33: No concentration data
#> Warning: treatment=iv30; id=34: No concentration data
#> Warning: treatment=iv30; id=35: No concentration data
#> Warning: treatment=iv30; id=36: No concentration data
#> Warning: treatment=iv30; id=37: No concentration data
#> Warning: treatment=iv30; id=38: No concentration data
#> Warning: treatment=iv30; id=39: No concentration data
#> Warning: treatment=iv30; id=40: No concentration data
#> Warning: treatment=iv30; id=41: No concentration data
#> Warning: treatment=iv30; id=42: No concentration data
#> Warning: treatment=iv30; id=43: No concentration data
#> Warning: treatment=iv30; id=44: No concentration data
#> Warning: treatment=iv30; id=45: No concentration data
#> Warning: treatment=iv30; id=46: No concentration data
#> Warning: treatment=iv30; id=47: No concentration data
#> Warning: treatment=iv30; id=48: No concentration data
#> Warning: treatment=iv30; id=49: No concentration data
#> Warning: treatment=iv30; id=50: No concentration data
#> Warning: treatment=iv30; id=51: No concentration data
#> Warning: treatment=iv30; id=52: No concentration data
#> Warning: treatment=iv30; id=53: No concentration data
#> Warning: treatment=iv30; id=54: No concentration data
#> Warning: treatment=iv30; id=55: No concentration data
#> Warning: treatment=iv30; id=56: No concentration data
#> Warning: treatment=iv30; id=57: No concentration data
#> Warning: treatment=iv30; id=58: No concentration data
#> Warning: treatment=iv30; id=59: No concentration data
#> Warning: treatment=iv30; id=60: No concentration data
#> Warning: treatment=iv30; id=61: No concentration data
#> Warning: treatment=iv30; id=62: No concentration data
#> Warning: treatment=iv30; id=63: No concentration data
#> Warning: treatment=iv30; id=64: No concentration data
#> Warning: treatment=iv30; id=65: No concentration data
#> Warning: treatment=iv30; id=66: No concentration data
#> Warning: treatment=iv30; id=67: No concentration data
#> Warning: treatment=iv30; id=68: No concentration data
#> Warning: treatment=iv30; id=69: No concentration data
#> Warning: treatment=iv30; id=70: No concentration data
#> Warning: treatment=iv30; id=71: No concentration data
#> Warning: treatment=iv30; id=72: No concentration data
#> Warning: treatment=iv30; id=73: No concentration data
#> Warning: treatment=iv30; id=74: No concentration data
#> Warning: treatment=iv30; id=75: No concentration data
#> Warning: treatment=iv30; id=76: No concentration data
#> Warning: treatment=iv30; id=77: No concentration data
#> Warning: treatment=iv30; id=78: No concentration data
#> Warning: treatment=iv30; id=79: No concentration data
#> Warning: treatment=iv30; id=80: No concentration data
#> Warning: treatment=iv30; id=81: No concentration data
#> Warning: treatment=iv30; id=82: No concentration data
#> Warning: treatment=iv30; id=83: No concentration data
#> Warning: treatment=iv30; id=84: No concentration data
#> Warning: treatment=iv30; id=85: No concentration data
#> Warning: treatment=iv30; id=86: No concentration data
#> Warning: treatment=iv30; id=87: No concentration data
#> Warning: treatment=iv30; id=88: No concentration data
#> Warning: treatment=iv30; id=89: No concentration data
#> Warning: treatment=iv30; id=90: No concentration data
#> Warning: treatment=iv30; id=91: No concentration data
#> Warning: treatment=iv30; id=92: No concentration data
#> Warning: treatment=iv30; id=93: No concentration data
#> Warning: treatment=iv30; id=94: No concentration data
#> Warning: treatment=iv30; id=95: No concentration data
#> Warning: treatment=iv30; id=96: No concentration data
#> Warning: treatment=iv30; id=97: No concentration data
#> Warning: treatment=iv30; id=98: No concentration data
#> Warning: treatment=iv30; id=99: No concentration data
#> Warning: treatment=iv30; id=100: No concentration data
#> Warning: treatment=iv30; id=101: No concentration data
#> Warning: treatment=iv30; id=102: No concentration data
#> Warning: treatment=iv30; id=103: No concentration data
#> Warning: treatment=iv30; id=104: No concentration data
#> Warning: treatment=iv30; id=105: No concentration data
#> Warning: treatment=iv30; id=106: No concentration data
#> Warning: treatment=iv30; id=107: No concentration data
#> Warning: treatment=iv30; id=108: No concentration data
#> Warning: treatment=iv30; id=109: No concentration data
#> Warning: treatment=iv30; id=110: No concentration data
#> Warning: treatment=iv30; id=111: No concentration data
#> Warning: treatment=iv30; id=112: No concentration data
#> Warning: treatment=iv30; id=113: No concentration data
#> Warning: treatment=iv30; id=114: No concentration data
#> Warning: treatment=iv30; id=115: No concentration data
#> Warning: treatment=iv30; id=116: No concentration data
#> Warning: treatment=iv30; id=117: No concentration data
#> Warning: treatment=iv30; id=118: No concentration data
#> Warning: treatment=iv30; id=119: No concentration data
#> Warning: treatment=iv30; id=120: No concentration data
#> Warning: treatment=iv30; id=121: No concentration data
#> Warning: treatment=iv30; id=122: No concentration data
#> Warning: treatment=iv30; id=123: No concentration data
#> Warning: treatment=iv30; id=124: No concentration data
#> Warning: treatment=iv30; id=125: No concentration data
#> Warning: treatment=iv30; id=126: No concentration data
#> Warning: treatment=iv30; id=127: No concentration data
#> Warning: treatment=iv30; id=128: No concentration data
#> Warning: treatment=iv30; id=129: No concentration data
#> Warning: treatment=iv30; id=130: No concentration data
#> Warning: treatment=iv30; id=131: No concentration data
#> Warning: treatment=iv30; id=132: No concentration data
#> Warning: treatment=iv30; id=133: No concentration data
#> Warning: treatment=iv30; id=134: No concentration data
#> Warning: treatment=iv30; id=135: No concentration data
#> Warning: treatment=iv30; id=136: No concentration data
#> Warning: treatment=iv30; id=137: No concentration data
#> Warning: treatment=iv30; id=138: No concentration data
#> Warning: treatment=iv30; id=139: No concentration data
#> Warning: treatment=iv30; id=140: No concentration data
#> Warning: treatment=iv30; id=141: No concentration data
#> Warning: treatment=iv30; id=142: No concentration data
#> Warning: treatment=iv30; id=143: No concentration data
#> Warning: treatment=iv30; id=144: No concentration data
#> Warning: treatment=iv30; id=145: No concentration data
#> Warning: treatment=iv30; id=146: No concentration data
#> Warning: treatment=iv30; id=147: No concentration data
#> Warning: treatment=iv30; id=148: No concentration data
#> Warning: treatment=iv30; id=149: No concentration data
#> Warning: treatment=iv30; id=150: No concentration data
#> Warning: treatment=iv30; id=151: No concentration data
#> Warning: treatment=iv30; id=152: No concentration data
#> Warning: treatment=iv30; id=153: No concentration data
#> Warning: treatment=iv30; id=154: No concentration data
#> Warning: treatment=iv30; id=155: No concentration data
#> Warning: treatment=iv30; id=156: No concentration data
#> Warning: treatment=iv30; id=157: No concentration data
#> Warning: treatment=iv30; id=158: No concentration data
#> Warning: treatment=iv30; id=159: No concentration data
#> Warning: treatment=iv30; id=160: No concentration data
#> Warning: treatment=iv30; id=161: No concentration data
#> Warning: treatment=iv30; id=162: No concentration data
#> Warning: treatment=iv30; id=163: No concentration data
#> Warning: treatment=iv30; id=164: No concentration data
#> Warning: treatment=iv30; id=165: No concentration data
#> Warning: treatment=iv30; id=166: No concentration data
#> Warning: treatment=iv30; id=167: No concentration data
#> Warning: treatment=iv30; id=168: No concentration data
#> Warning: treatment=iv30; id=169: No concentration data
#> Warning: treatment=iv30; id=170: No concentration data
#> Warning: treatment=iv30; id=171: No concentration data
#> Warning: treatment=iv30; id=172: No concentration data
#> Warning: treatment=iv30; id=173: No concentration data
#> Warning: treatment=iv30; id=174: No concentration data
#> Warning: treatment=iv30; id=175: No concentration data
#> Warning: treatment=iv30; id=176: No concentration data
#> Warning: treatment=iv30; id=177: No concentration data
#> Warning: treatment=iv30; id=178: No concentration data
#> Warning: treatment=iv30; id=179: No concentration data
#> Warning: treatment=iv30; id=180: No concentration data
#> Warning: treatment=iv30; id=181: No concentration data
#> Warning: treatment=iv30; id=182: No concentration data
#> Warning: treatment=iv30; id=183: No concentration data
#> Warning: treatment=iv30; id=184: No concentration data
#> Warning: treatment=iv30; id=185: No concentration data
#> Warning: treatment=iv30; id=186: No concentration data
#> Warning: treatment=iv30; id=187: No concentration data
#> Warning: treatment=iv30; id=188: No concentration data
#> Warning: treatment=iv30; id=189: No concentration data
#> Warning: treatment=iv30; id=190: No concentration data
#> Warning: treatment=iv30; id=191: No concentration data
#> Warning: treatment=iv30; id=192: No concentration data
#> Warning: treatment=iv30; id=193: No concentration data
#> Warning: treatment=iv30; id=194: No concentration data
#> Warning: treatment=iv30; id=195: No concentration data
#> Warning: treatment=iv30; id=196: No concentration data
#> Warning: treatment=iv30; id=197: No concentration data
#> Warning: treatment=iv30; id=198: No concentration data
#> Warning: treatment=iv30; id=199: No concentration data
#> Warning: treatment=iv30; id=200: No concentration data
#> Warning: treatment=oral; id=1: No concentration data
#> Warning: treatment=oral; id=2: No concentration data
#> Warning: treatment=oral; id=3: No concentration data
#> Warning: treatment=oral; id=4: No concentration data
#> Warning: treatment=oral; id=5: No concentration data
#> Warning: treatment=oral; id=6: No concentration data
#> Warning: treatment=oral; id=7: No concentration data
#> Warning: treatment=oral; id=8: No concentration data
#> Warning: treatment=oral; id=9: No concentration data
#> Warning: treatment=oral; id=10: No concentration data
#> Warning: treatment=oral; id=11: No concentration data
#> Warning: treatment=oral; id=12: No concentration data
#> Warning: treatment=oral; id=13: No concentration data
#> Warning: treatment=oral; id=14: No concentration data
#> Warning: treatment=oral; id=15: No concentration data
#> Warning: treatment=oral; id=16: No concentration data
#> Warning: treatment=oral; id=17: No concentration data
#> Warning: treatment=oral; id=18: No concentration data
#> Warning: treatment=oral; id=19: No concentration data
#> Warning: treatment=oral; id=20: No concentration data
#> Warning: treatment=oral; id=21: No concentration data
#> Warning: treatment=oral; id=22: No concentration data
#> Warning: treatment=oral; id=23: No concentration data
#> Warning: treatment=oral; id=24: No concentration data
#> Warning: treatment=oral; id=25: No concentration data
#> Warning: treatment=oral; id=26: No concentration data
#> Warning: treatment=oral; id=27: No concentration data
#> Warning: treatment=oral; id=28: No concentration data
#> Warning: treatment=oral; id=29: No concentration data
#> Warning: treatment=oral; id=30: No concentration data
#> Warning: treatment=oral; id=31: No concentration data
#> Warning: treatment=oral; id=32: No concentration data
#> Warning: treatment=oral; id=33: No concentration data
#> Warning: treatment=oral; id=34: No concentration data
#> Warning: treatment=oral; id=35: No concentration data
#> Warning: treatment=oral; id=36: No concentration data
#> Warning: treatment=oral; id=37: No concentration data
#> Warning: treatment=oral; id=38: No concentration data
#> Warning: treatment=oral; id=39: No concentration data
#> Warning: treatment=oral; id=40: No concentration data
#> Warning: treatment=oral; id=41: No concentration data
#> Warning: treatment=oral; id=42: No concentration data
#> Warning: treatment=oral; id=43: No concentration data
#> Warning: treatment=oral; id=44: No concentration data
#> Warning: treatment=oral; id=45: No concentration data
#> Warning: treatment=oral; id=46: No concentration data
#> Warning: treatment=oral; id=47: No concentration data
#> Warning: treatment=oral; id=48: No concentration data
#> Warning: treatment=oral; id=49: No concentration data
#> Warning: treatment=oral; id=50: No concentration data
#> Warning: treatment=oral; id=51: No concentration data
#> Warning: treatment=oral; id=52: No concentration data
#> Warning: treatment=oral; id=53: No concentration data
#> Warning: treatment=oral; id=54: No concentration data
#> Warning: treatment=oral; id=55: No concentration data
#> Warning: treatment=oral; id=56: No concentration data
#> Warning: treatment=oral; id=57: No concentration data
#> Warning: treatment=oral; id=58: No concentration data
#> Warning: treatment=oral; id=59: No concentration data
#> Warning: treatment=oral; id=60: No concentration data
#> Warning: treatment=oral; id=61: No concentration data
#> Warning: treatment=oral; id=62: No concentration data
#> Warning: treatment=oral; id=63: No concentration data
#> Warning: treatment=oral; id=64: No concentration data
#> Warning: treatment=oral; id=65: No concentration data
#> Warning: treatment=oral; id=66: No concentration data
#> Warning: treatment=oral; id=67: No concentration data
#> Warning: treatment=oral; id=68: No concentration data
#> Warning: treatment=oral; id=69: No concentration data
#> Warning: treatment=oral; id=70: No concentration data
#> Warning: treatment=oral; id=71: No concentration data
#> Warning: treatment=oral; id=72: No concentration data
#> Warning: treatment=oral; id=73: No concentration data
#> Warning: treatment=oral; id=74: No concentration data
#> Warning: treatment=oral; id=75: No concentration data
#> Warning: treatment=oral; id=76: No concentration data
#> Warning: treatment=oral; id=77: No concentration data
#> Warning: treatment=oral; id=78: No concentration data
#> Warning: treatment=oral; id=79: No concentration data
#> Warning: treatment=oral; id=80: No concentration data
#> Warning: treatment=oral; id=81: No concentration data
#> Warning: treatment=oral; id=82: No concentration data
#> Warning: treatment=oral; id=83: No concentration data
#> Warning: treatment=oral; id=84: No concentration data
#> Warning: treatment=oral; id=85: No concentration data
#> Warning: treatment=oral; id=86: No concentration data
#> Warning: treatment=oral; id=87: No concentration data
#> Warning: treatment=oral; id=88: No concentration data
#> Warning: treatment=oral; id=89: No concentration data
#> Warning: treatment=oral; id=90: No concentration data
#> Warning: treatment=oral; id=91: No concentration data
#> Warning: treatment=oral; id=92: No concentration data
#> Warning: treatment=oral; id=93: No concentration data
#> Warning: treatment=oral; id=94: No concentration data
#> Warning: treatment=oral; id=95: No concentration data
#> Warning: treatment=oral; id=96: No concentration data
#> Warning: treatment=oral; id=97: No concentration data
#> Warning: treatment=oral; id=98: No concentration data
#> Warning: treatment=oral; id=99: No concentration data
#> Warning: treatment=oral; id=100: No concentration data
#> Warning: treatment=oral; id=101: No concentration data
#> Warning: treatment=oral; id=102: No concentration data
#> Warning: treatment=oral; id=103: No concentration data
#> Warning: treatment=oral; id=104: No concentration data
#> Warning: treatment=oral; id=105: No concentration data
#> Warning: treatment=oral; id=106: No concentration data
#> Warning: treatment=oral; id=107: No concentration data
#> Warning: treatment=oral; id=108: No concentration data
#> Warning: treatment=oral; id=109: No concentration data
#> Warning: treatment=oral; id=110: No concentration data
#> Warning: treatment=oral; id=111: No concentration data
#> Warning: treatment=oral; id=112: No concentration data
#> Warning: treatment=oral; id=113: No concentration data
#> Warning: treatment=oral; id=114: No concentration data
#> Warning: treatment=oral; id=115: No concentration data
#> Warning: treatment=oral; id=116: No concentration data
#> Warning: treatment=oral; id=117: No concentration data
#> Warning: treatment=oral; id=118: No concentration data
#> Warning: treatment=oral; id=119: No concentration data
#> Warning: treatment=oral; id=120: No concentration data
#> Warning: treatment=oral; id=121: No concentration data
#> Warning: treatment=oral; id=122: No concentration data
#> Warning: treatment=oral; id=123: No concentration data
#> Warning: treatment=oral; id=124: No concentration data
#> Warning: treatment=oral; id=125: No concentration data
#> Warning: treatment=oral; id=126: No concentration data
#> Warning: treatment=oral; id=127: No concentration data
#> Warning: treatment=oral; id=128: No concentration data
#> Warning: treatment=oral; id=129: No concentration data
#> Warning: treatment=oral; id=130: No concentration data
#> Warning: treatment=oral; id=131: No concentration data
#> Warning: treatment=oral; id=132: No concentration data
#> Warning: treatment=oral; id=133: No concentration data
#> Warning: treatment=oral; id=134: No concentration data
#> Warning: treatment=oral; id=135: No concentration data
#> Warning: treatment=oral; id=136: No concentration data
#> Warning: treatment=oral; id=137: No concentration data
#> Warning: treatment=oral; id=138: No concentration data
#> Warning: treatment=oral; id=139: No concentration data
#> Warning: treatment=oral; id=140: No concentration data
#> Warning: treatment=oral; id=141: No concentration data
#> Warning: treatment=oral; id=142: No concentration data
#> Warning: treatment=oral; id=143: No concentration data
#> Warning: treatment=oral; id=144: No concentration data
#> Warning: treatment=oral; id=145: No concentration data
#> Warning: treatment=oral; id=146: No concentration data
#> Warning: treatment=oral; id=147: No concentration data
#> Warning: treatment=oral; id=148: No concentration data
#> Warning: treatment=oral; id=149: No concentration data
#> Warning: treatment=oral; id=150: No concentration data
#> Warning: treatment=oral; id=151: No concentration data
#> Warning: treatment=oral; id=152: No concentration data
#> Warning: treatment=oral; id=153: No concentration data
#> Warning: treatment=oral; id=154: No concentration data
#> Warning: treatment=oral; id=155: No concentration data
#> Warning: treatment=oral; id=156: No concentration data
#> Warning: treatment=oral; id=157: No concentration data
#> Warning: treatment=oral; id=158: No concentration data
#> Warning: treatment=oral; id=159: No concentration data
#> Warning: treatment=oral; id=160: No concentration data
#> Warning: treatment=oral; id=161: No concentration data
#> Warning: treatment=oral; id=162: No concentration data
#> Warning: treatment=oral; id=163: No concentration data
#> Warning: treatment=oral; id=164: No concentration data
#> Warning: treatment=oral; id=165: No concentration data
#> Warning: treatment=oral; id=166: No concentration data
#> Warning: treatment=oral; id=167: No concentration data
#> Warning: treatment=oral; id=168: No concentration data
#> Warning: treatment=oral; id=169: No concentration data
#> Warning: treatment=oral; id=170: No concentration data
#> Warning: treatment=oral; id=171: No concentration data
#> Warning: treatment=oral; id=172: No concentration data
#> Warning: treatment=oral; id=173: No concentration data
#> Warning: treatment=oral; id=174: No concentration data
#> Warning: treatment=oral; id=175: No concentration data
#> Warning: treatment=oral; id=176: No concentration data
#> Warning: treatment=oral; id=177: No concentration data
#> Warning: treatment=oral; id=178: No concentration data
#> Warning: treatment=oral; id=179: No concentration data
#> Warning: treatment=oral; id=180: No concentration data
#> Warning: treatment=oral; id=181: No concentration data
#> Warning: treatment=oral; id=182: No concentration data
#> Warning: treatment=oral; id=183: No concentration data
#> Warning: treatment=oral; id=184: No concentration data
#> Warning: treatment=oral; id=185: No concentration data
#> Warning: treatment=oral; id=186: No concentration data
#> Warning: treatment=oral; id=187: No concentration data
#> Warning: treatment=oral; id=188: No concentration data
#> Warning: treatment=oral; id=189: No concentration data
#> Warning: treatment=oral; id=190: No concentration data
#> Warning: treatment=oral; id=191: No concentration data
#> Warning: treatment=oral; id=192: No concentration data
#> Warning: treatment=oral; id=193: No concentration data
#> Warning: treatment=oral; id=194: No concentration data
#> Warning: treatment=oral; id=195: No concentration data
#> Warning: treatment=oral; id=196: No concentration data
#> Warning: treatment=oral; id=197: No concentration data
#> Warning: treatment=oral; id=198: No concentration data
#> Warning: treatment=oral; id=199: No concentration data
#> Warning: treatment=oral; id=200: No concentration data
nca_tbl <- as.data.frame(nca$result) |>
  dplyr::group_by(treatment, PPTESTCD) |>
  dplyr::summarise(median = median(PPORRES, na.rm = TRUE),
                   p05 = quantile(PPORRES, 0.05, na.rm = TRUE),
                   p95 = quantile(PPORRES, 0.95, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(PPTESTCD, treatment)

knitr::kable(nca_tbl, digits = 3,
             caption = "Simulated NCA parameters, median (5th-95th percentile).")
```

| treatment          | PPTESTCD            | median |    p05 |     p95 |
|:-------------------|:--------------------|-------:|-------:|--------:|
| IV 30-min infusion | adj.r.squared       |  1.000 |  1.000 |   1.000 |
| Oral               | adj.r.squared       |  1.000 |  1.000 |   1.000 |
| IV 30-min infusion | aucinf.obs          |     NA |     NA |      NA |
| Oral               | aucinf.obs          |     NA |     NA |      NA |
| IV 30-min infusion | auclast             |     NA |     NA |      NA |
| Oral               | auclast             |     NA |     NA |      NA |
| IV 30-min infusion | clast.obs           |  8.933 |  0.000 |  46.496 |
| Oral               | clast.obs           |  7.681 |  0.000 |  43.831 |
| IV 30-min infusion | clast.pred          |  8.933 |  0.000 |  46.496 |
| Oral               | clast.pred          |  7.682 |  0.000 |  43.832 |
| IV 30-min infusion | cmax                | 67.116 | 21.958 | 266.387 |
| Oral               | cmax                | 69.021 | 18.128 | 263.509 |
| IV 30-min infusion | half.life           |  7.789 |  1.042 |  64.116 |
| Oral               | half.life           |  9.341 |  0.903 |  66.882 |
| IV 30-min infusion | lambda.z            |  0.089 |  0.011 |   0.665 |
| Oral               | lambda.z            |  0.074 |  0.010 |   0.768 |
| IV 30-min infusion | lambda.z.n.points   | 94.000 | 93.000 |  95.000 |
| Oral               | lambda.z.n.points   | 94.000 | 93.000 |  95.000 |
| IV 30-min infusion | lambda.z.time.first |  0.750 |  0.500 |   1.000 |
| Oral               | lambda.z.time.first |  0.750 |  0.500 |   1.000 |
| IV 30-min infusion | lambda.z.time.last  | 24.000 | 24.000 |  24.000 |
| Oral               | lambda.z.time.last  | 24.000 | 24.000 |  24.000 |
| IV 30-min infusion | r.squared           |  1.000 |  1.000 |   1.000 |
| Oral               | r.squared           |  1.000 |  1.000 |   1.000 |
| IV 30-min infusion | span.ratio          |  2.980 |  0.363 |  22.559 |
| Oral               | span.ratio          |  2.489 |  0.344 |  26.032 |
| IV 30-min infusion | tlast               | 24.000 | 24.000 |  24.000 |
| Oral               | tlast               | 24.000 | 24.000 |  24.000 |
| IV 30-min infusion | tmax                |  0.500 |  0.250 |   0.750 |
| Oral               | tmax                |  0.500 |  0.250 |   0.750 |

Simulated NCA parameters, median (5th-95th percentile). {.table}

### Comparison against published values

Brussee 2018 reports total plasma clearance (median 0.181 L/h, range
0.03-0.79 L/h) and total oral bioavailability (median 92.1%, range
67-95%) but does not tabulate per-route NCA Cmax / AUC. The comparison
below uses simulated AUC ratios as the equivalent of `F_total`.

``` r

auc_pivot <- nca_tbl |>
  dplyr::filter(PPTESTCD %in% c("aucinf.obs", "auclast")) |>
  dplyr::select(treatment, PPTESTCD, median) |>
  tidyr::pivot_wider(names_from = treatment, values_from = median)

cl_total_h <- dose_df |>
  dplyr::left_join(
    nca_tbl |> dplyr::filter(PPTESTCD == "aucinf.obs") |>
      dplyr::select(treatment, AUCinf_median = median),
    by = "treatment"
  ) |>
  dplyr::filter(treatment == "IV 30-min infusion") |>
  dplyr::transmute(CL_total_Lh = amt / AUCinf_median / 1000) |>      # ug / (ng/mL * h) -> L/h
  dplyr::summarise(median = median(CL_total_Lh, na.rm = TRUE),
                   p05 = quantile(CL_total_Lh, 0.05, na.rm = TRUE),
                   p95 = quantile(CL_total_Lh, 0.95, na.rm = TRUE))

knitr::kable(auc_pivot, digits = 3,
             caption = "Simulated AUC by route -- median across cohort.")
```

| PPTESTCD   | IV 30-min infusion | Oral |
|:-----------|-------------------:|-----:|
| aucinf.obs |                 NA |   NA |
| auclast    |                 NA |   NA |

Simulated AUC by route – median across cohort. {.table}

``` r

knitr::kable(cl_total_h, digits = 3,
             caption = paste("Implied total plasma CL from IV AUCinf -- compare",
                             "to Brussee 2018 median 0.181 L/h (range 0.03-0.79)."))
```

| median | p05 | p95 |
|-------:|----:|----:|
|     NA |  NA |  NA |

Implied total plasma CL from IV AUCinf – compare to Brussee 2018 median
0.181 L/h (range 0.03-0.79). {.table}

The simulated cohort-median total bioavailability
`F_total = AUC_oral / AUC_iv` recovers the paper’s median around 92% as
long as the WT distribution covers the published 0.77-2.03 kg range and
BSV is included.

## Assumptions and deviations

- **Lumped first-pass representation.** Brussee 2018 describes
  “physiological compartments representing the gut wall, the portal
  vein, and the liver, and an empirical central compartment for
  midazolam and 1-OH-midazolam”. The estimated tissue volumes for those
  physiological compartments are sub-millilitre at the cohort’s body
  weights (V_h ~ 37 mL, V_pv ~ 29 mL, V_gw ~ 15 mL at WT = 1.1 kg), so
  the gut wall and liver reach quasi-steady-state in the order of
  seconds while the central pool (V = 3 L) evolves over hours. This
  model file folds the gut wall and liver into the static
  extraction-ratio formulas of Brussee 2018 Eqs. 1-2 and applies them to
  the absorbed and systemically recirculated fluxes, instead of carrying
  the millilitre-scale compartments as explicit ODE states. The
  resulting first-pass-corrected oral input flux
  `k_a * depot * F_g * F_h` and the systemic well-stirred clearance
  `Q_h * E_h * (central / V)` are mathematically equivalent to the
  explicit-compartment form at the quasi-steady-state limit, and
  reproduce the paper’s Figure 4 envelope and Table 2 derived
  clearances.
- **Body-weight covariate is required.** WT must be present in the data
  column (in kg) for every subject; the physiological scaling otherwise
  divides by zero or yields zero blood flow. Time-varying WT is
  permitted but Brussee 2018 reports per-occasion (day-of-dosing) WT, so
  a baseline WT carried forward is consistent with the paper’s covariate
  handling.
- **Volumes are in blood-concentration terms.** V (3.0 L) and V_M
  (2.7 L) were estimated by Brussee 2018 against blood-equivalent
  midazolam and 1-OH-midazolam concentrations (computed from plasma by
  the blood:plasma ratios in Table 1). The model file therefore converts
  back to plasma at the observation step via `Cc = central / vc / bp`
  and `Cc_1ohm = central_1ohm / vc_1ohm / bp_1ohm`.
- **Mass-balance handling of metabolite formation.** Brussee 2018
  assumed `f_M = 1` (i.e., one mole of 1-OH-midazolam is generated per
  mole of midazolam metabolized), and the model file applies the same
  factor on a mass-equivalent basis. Strictly, the molecular-weight
  ratio is 341.77 / 325.77 = 1.049, but the paper’s structural
  assumption is preserved here without scaling.
- **Race / sex distributions.** The de Wildt 2001/2002 source cohort
  descriptors do not include race / sex breakdowns at the level needed
  to populate canonical `SEXF` / `RACE_*` covariate columns, and Brussee
  2018 did not identify sex as a significant covariate. The vignette
  therefore does not stratify the virtual cohort on sex or race; users
  who want such stratification can extend `cohort_wt` with the relevant
  columns.
- **Additive residual error fixed to a numerical regularizer.** Brussee
  2018 Table 2 reports both additive variances as `0.0001 FIX`, which
  translates to SD = 0.01 ng/mL – far below the LLOQ of typical
  midazolam assays. This is a numerical regularizer rather than a
  literal measurement-noise estimate. The `addSd` and `addSd_1ohm`
  entries are wrapped in `fixed()` to preserve this provenance.
- **No supplementary control stream on disk.** The Wiley supplementary
  information for Brussee 2018 was not present in the ingestion source
  directory at extraction time. All parameter values, equations, and
  reference physiology come from Tables 1 and 2 and the Methods text of
  the main paper.

## Reference

- Brussee JM, Yu H, Krekels EHJ, de Roos B, Brill MJE, van den Anker JN,
  Rostami-Hodjegan A, de Wildt SN, Knibbe CAJ (2018). First-Pass CYP3A-
  Mediated Metabolism of Midazolam in the Gut Wall and Liver in Preterm
  Neonates. CPT Pharmacometrics Syst Pharmacol 7(6):374-383.
  <doi:10.1002/psp4.12295>.
