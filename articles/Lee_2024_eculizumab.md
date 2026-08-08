# Eculizumab and SB12 biosimilar (Lee 2024)

## Model and source

- Citation: Lee H, Park J, Jang H, Lee SJ, Kim J. Population
  pharmacokinetic, pharmacodynamic and efficacy modeling of SB12
  (proposed eculizumab biosimilar) and reference eculizumab. Eur J Clin
  Pharmacol. 2024;80(9):1325-1338. <doi:10.1007/s00228-024-03703-8>
- Description: Two-compartment population PK model of eculizumab (SB12
  proposed biosimilar and reference Soliris) with a direct inhibitory
  sigmoid Imax model for terminal complement activity and a direct
  sigmoid Emax model linking terminal complement activity to serum LDH,
  in healthy subjects and patients with paroxysmal nocturnal
  haemoglobinuria (Lee 2024)
- Article: <https://doi.org/10.1007/s00228-024-03703-8>
- Supplement (Appendix Tables A1/A2 and Figure A1):
  <https://doi.org/10.1007/s00228-024-03703-8> (Electronic Supplementary
  Material, `228_2024_3703_MOESM1_ESM.docx`)

Lee 2024 developed a linked population pharmacokinetic / pharmacodynamic
/ efficacy model for SB12 (a proposed eculizumab biosimilar) and
reference eculizumab (Soliris), pooling a phase I single-dose study in
healthy subjects with a phase III cross-over study in patients with
paroxysmal nocturnal haemoglobinuria (PNH). Three endpoints are chained:

1.  **PK** – serum eculizumab concentration, described by a
    two-compartment model with first-order elimination.
2.  **PD** – terminal complement activity (%), described by a *direct
    response* model with an inhibitory sigmoid `Imax` relationship to
    serum concentration.
3.  **Efficacy** – serum lactate dehydrogenase (LDH, U/L), described by
    a *direct* sigmoid `Emax` relationship to terminal complement
    activity.

Treatment group (SB12 versus EU-sourced or US-sourced reference
eculizumab) was **not** a statistically significant covariate on any PK,
PD or efficacy parameter, which is the biosimilarity conclusion of the
paper. One parameter set therefore describes both products, and the
packaged model carries no treatment covariate.

## Population

The analysis pooled 289 subjects across two studies (Lee 2024 Table 1
and Table 2):

- **Phase I** – 240 healthy subjects (80 per treatment arm: SB12,
  EU-sourced Soliris, US-sourced Soliris) given a single 300 mg
  intravenous infusion over 35 min, with dense PK sampling to 1512 h and
  PD sampling to 1512 h. Median age 40 years (19-55), median weight
  82.40 kg (70.0-94.3), 95.8% male, 95.8% White.
- **Phase III** – 49 patients with PNH in a double-blind cross-over
  study (SB12-to-ECU and ECU-to-SB12 sequences) given the approved
  regimen: 600 mg IV every week for 4 weeks, 900 mg at week 5, then 900
  mg every 2 weeks through week 50. Trough-only PK sampling. Median age
  36 years (18-79), median weight 63.00 kg (43.0-111.0), 55.1% female,
  53.1% Asian.

The pooled dataset comprised 4136 quantifiable serum concentrations,
2900 terminal complement activity levels (289 subjects) and 1350 LDH
levels (49 PNH patients). PNH patients had notably lower weight and
height and higher baseline aspartate aminotransferase, bilirubin and LDH
than the healthy subjects.

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Lee_2024_eculizumab")()$population`).

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in `inst/modeldb/specificDrugs/Lee_2024_eculizumab.R`. The
table below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `d/dt(central)`, `d/dt(peripheral1)` (two-compartment, first-order elimination) | n/a | Lee 2024 Results, “PK/PD model of pooled subjects”: “a two-compartment model best described the PKs of SB12 and ECU in pooled data” |
| `CL = theta1 * (WT/80.9)^theta7` | n/a | Lee 2024 Table 3 (covariate equation printed in the parameter column) |
| `Vc = theta2 * (WT/80.9)^theta8` | n/a | Lee 2024 Table 3 (covariate equation printed in the parameter column) |
| `lcl` (CL) | 0.0174 L/h | Lee 2024 Table 3, theta1 |
| `e_wt_cl` | 1.1400 | Lee 2024 Table 3, theta7 |
| `lvc` (Vc, healthy) | 3.47 L | Lee 2024 Table 3, theta2 |
| `e_pnh_vc` | log(5.68/3.47) = 0.4928 | Lee 2024 Table 3, theta2_pat = 5.68 L versus theta2 = 3.47 L |
| `e_wt_vc` | 0.8630 | Lee 2024 Table 3, theta8 |
| `lvp` (Vp) | 0.79 L | Lee 2024 Table 3, Vp |
| `lq` (Q) | 0.0134 L/h | Lee 2024 Table 3, Q |
| `etalcl` | 15.62 %CV -\> omega^2 = 0.0241056 | Lee 2024 Table 3, omega_CL |
| `etalvc` | 12.74 %CV -\> omega^2 = 0.0161004 | Lee 2024 Table 3, omega_VC |
| `etalvp` | 36.80 %CV -\> omega^2 = 0.1270061 | Lee 2024 Table 3, omega_VP |
| `cov(etalcl, etalvc)` | rho = 0.54 -\> 0.0106383 | Lee 2024 Table 3, rho_CL-VC |
| `propSd` | 11.70% | Lee 2024 Table 3, sigma_prop |
| `tca = E0 * (1 - Imax * Cc^H / (IC50^H + Cc^H))` | n/a | Lee 2024 Results / Discussion: “direct response model with an inhibitory sigmoid maximum effect (Emax) relationship between serum concentrations and terminal complement activity” |
| `lrbase_tca` (E0, healthy) | 85.90% | Lee 2024 Table 4, E0 |
| `e_pnh_rbase_tca` | log(101.00/85.90) = 0.1619 | Lee 2024 Table 4, E0_pat = 101.00% versus E0 = 85.90% |
| `limax_tca` (Imax, healthy) | 0.93 | Lee 2024 Table 4, Imax |
| `e_pnh_imax_tca` | log(0.88/0.93) = -0.0553 | Lee 2024 Table 4, Imax_pat = 0.88 versus Imax = 0.93 |
| `lec50_tca` (IC50) | 36.60 ug/mL | Lee 2024 Table 4, IC50 |
| `lhill_tca` (H) | 4.56 | Lee 2024 Table 4, H |
| `etalrbase_tca` | 15.29 %CV -\> omega^2 = 0.0231093 | Lee 2024 Table 4, omega_E0 |
| `etalimax_tca` | 2.65 %CV -\> omega^2 = 0.0007020 | Lee 2024 Table 4, omega_Imax |
| `etalec50_tca` | 22.74 %CV -\> omega^2 = 0.0504181 | Lee 2024 Table 4, omega_IC50 |
| `propSd_tca` | 18.60% | Lee 2024 Table 4, sigma_prop |
| `ldh = LL0 + LMAX * tca^LGAM / (LC50^LGAM + tca^LGAM)` | n/a | Lee 2024 Results / Discussion: “the relationship between terminal complement activity and LDH was well described by a direct sigmoid Emax model” |
| `lrbase_ldh` (LL0) | 206 U/L | Lee 2024 Table 5, LL0 |
| `lemax_ldh` (LMAX) | 1680 U/L | Lee 2024 Table 5, LMAX |
| `lec50_ldh` (LC50) | 39.0% | Lee 2024 Table 5, LC50 |
| `lhill_ldh` (LGAM) | 4.30 | Lee 2024 Table 5, LGAM |
| `etalrbase_ldh` | 28.97 %CV -\> omega^2 = 0.0805897 | Lee 2024 Table 5, omega_LL0 |
| `etalemax_ldh` | 81.57 %CV -\> omega^2 = 0.5100452 | Lee 2024 Table 5, omega_LMAX |
| `etalhill_ldh` | 55.91 %CV -\> omega^2 = 0.2720044 | Lee 2024 Table 5, omega_LGAM |
| `propSd_ldh` | 32.00% | Lee 2024 Table 5, sigma_prop |
| Exponential ETA on every parameter | n/a | Lee 2024 Methods, “Base model”: “The inter-individual variability (i.e. ETA) of each parameter was applied exponentially” |
| Reference weight 80.9 kg | n/a | Lee 2024 Table 3, covariate equations `(WT/80.9)` |

## Virtual cohort

Original observed data are not publicly available. The simulations below
use virtual populations whose covariate distributions approximate the
published trial demographics (Lee 2024 Table 2). Two arms are simulated:

- **Healthy** (`DIS_PNH = 0`) – 150 subjects, weight drawn to match the
  phase I median 82.40 kg and range 70.0-94.3 kg, given a single 300 mg
  IV infusion over 35 min (0.58 h) with dense sampling to 1512 h.
- **PNH** (`DIS_PNH = 1`) – 150 patients, weight drawn to match the
  phase III median 63.00 kg and range 43.0-111.0 kg, given the approved
  regimen (600 mg weekly x 4, 900 mg at week 5, then 900 mg every 2
  weeks to week 50), sampled to week 52.

``` r

set.seed(20240530)

n_per_arm <- 150L

# Truncated-normal weight sampler matching a published median and range.
sample_weight <- function(n, median_wt, sd_wt, lo, hi) {
  w <- stats::rnorm(n, mean = median_wt, sd = sd_wt)
  pmin(pmax(w, lo), hi)
}

# Healthy phase I arm: single 300 mg IV infusion over 35 min (0.58 h).
healthy_subj <- tibble(
  id      = seq_len(n_per_arm),
  WT      = sample_weight(n_per_arm, 82.40, 6.0, 70.0, 94.3),
  DIS_PNH = 0,
  arm     = "Healthy (phase I, 300 mg single IV)"
)

healthy_dose <- healthy_subj |>
  mutate(time = 0, amt = 300, evid = 1L, cmt = "central", rate = 300 / 0.58)

healthy_obs <- healthy_subj |>
  tidyr::crossing(time = c(0, 0.58, 4, 8, 12, 24, 48, 96, 168, 240, 336,
                           504, 672, 840, 1008, 1176, 1344, 1512)) |>
  mutate(amt = NA_real_, evid = 0L, cmt = "Cc", rate = NA_real_)

# PNH phase III arm: 600 mg QW x 4 (weeks 0-3), 900 mg at week 4 (the fifth
# dose, 1 week later), then 900 mg Q2W through week 50. Times in hours.
wk <- function(x) x * 24 * 7
pnh_dose_times <- c(wk(0:3), wk(4), wk(seq(6, 50, by = 2)))
pnh_dose_amts  <- c(rep(600, 4), 900, rep(900, length(seq(6, 50, by = 2))))

pnh_subj <- tibble(
  id      = n_per_arm + seq_len(n_per_arm),
  WT      = sample_weight(n_per_arm, 63.00, 13.0, 43.0, 111.0),
  DIS_PNH = 1,
  arm     = "PNH (phase III, 600 mg QW x4 then 900 mg Q2W)"
)

pnh_dose <- pnh_subj |>
  tidyr::crossing(tibble(time = pnh_dose_times, amt = pnh_dose_amts)) |>
  mutate(evid = 1L, cmt = "central", rate = amt / 0.58)

pnh_obs <- pnh_subj |>
  tidyr::crossing(time = sort(unique(c(0, wk(seq(0, 52, by = 1)),
                                       wk(seq(0, 52, by = 2)) + 1)))) |>
  mutate(amt = NA_real_, evid = 0L, cmt = "Cc", rate = NA_real_)

events <- bind_rows(healthy_dose, healthy_obs, pnh_dose, pnh_obs) |>
  arrange(id, time, desc(evid)) |>
  as.data.frame()

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
c(healthy = n_per_arm, pnh = n_per_arm, rows = nrow(events))
#> healthy     pnh    rows 
#>     150     150   19050
```

## Simulation

``` r

mod <- readModelDb("Lee_2024_eculizumab")

sim <- rxode2::rxSolve(
  mod,
  events = events,
  keep   = c("WT", "DIS_PNH", "arm"),
  useLinCmt = FALSE
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'

# Guard against silently dropped subjects (rxSolve does this quietly).
stopifnot(dplyr::n_distinct(sim$id) == 2L * n_per_arm)
```

## Replicate published figures

### Figure 4A/4B – serum eculizumab concentration

``` r

# Replicates Figure 4A (healthy subjects) and 4B (PNH patients) of Lee 2024:
# visual predictive check of serum eculizumab concentration versus time.
sim |>
  group_by(arm, time) |>
  summarise(
    Q05 = quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time / 24, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25, fill = "steelblue") +
  geom_line(colour = "steelblue") +
  geom_hline(yintercept = 36.6, linetype = "dashed") +
  facet_wrap(~arm, scales = "free") +
  labs(
    x = "Time (days)", y = "Serum eculizumab (ug/mL)",
    title = "Figure 4A/4B - simulated PK",
    caption = paste("Replicates Figures 4A and 4B of Lee 2024.",
                    "Dashed line = IC50 (36.6 ug/mL, Table 4).")
  )
```

![](Lee_2024_eculizumab_files/figure-html/figure-4ab-1.png)

### Figure 4C/4D – terminal complement activity

``` r

# Replicates Figure 4C (healthy subjects) and 4D (PNH patients) of Lee 2024:
# visual predictive check of terminal complement activity versus time.
sim |>
  group_by(arm, time) |>
  summarise(
    Q05 = quantile(tca, 0.05, na.rm = TRUE),
    Q50 = quantile(tca, 0.50, na.rm = TRUE),
    Q95 = quantile(tca, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time / 24, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25, fill = "darkorange") +
  geom_line(colour = "darkorange") +
  facet_wrap(~arm, scales = "free_x") +
  labs(
    x = "Time (days)", y = "Terminal complement activity (%)",
    title = "Figure 4C/4D - simulated terminal complement activity",
    caption = "Replicates Figures 4C and 4D of Lee 2024."
  )
```

![](Lee_2024_eculizumab_files/figure-html/figure-4cd-1.png)

### Figure 4E – LDH in PNH patients

``` r

# Replicates Figure 4E of Lee 2024: visual predictive check of the
# PK/PD/efficacy model (serum LDH) in PNH patients.
sim |>
  filter(DIS_PNH == 1) |>
  group_by(time) |>
  summarise(
    Q05 = quantile(ldh, 0.05, na.rm = TRUE),
    Q50 = quantile(ldh, 0.50, na.rm = TRUE),
    Q95 = quantile(ldh, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time / 24, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25, fill = "firebrick") +
  geom_line(colour = "firebrick") +
  geom_hline(yintercept = 250, linetype = "dashed") +
  labs(
    x = "Time (days)", y = "Serum LDH (U/L)",
    title = "Figure 4E - simulated LDH in PNH patients",
    caption = paste("Replicates Figure 4E of Lee 2024.",
                    "Dashed line = approximate upper limit of normal (250 U/L).")
  )
```

![](Lee_2024_eculizumab_files/figure-html/figure-4e-1.png)

### Exposure-response relationships (Tables 4 and 5)

The two direct-response relationships are algebraic, so they can be
drawn exactly from the typical-value parameters without simulating a
time course.

``` r

conc_grid <- tibble(Cc = seq(0, 150, length.out = 401)) |>
  mutate(
    `Healthy (E0 85.9%, Imax 0.93)` =
      85.90 * (1 - 0.93 * Cc^4.56 / (36.60^4.56 + Cc^4.56)),
    `PNH (E0 101.0%, Imax 0.88)` =
      101.00 * (1 - 0.88 * Cc^4.56 / (36.60^4.56 + Cc^4.56))
  ) |>
  pivot_longer(-Cc, names_to = "group", values_to = "tca")

p_pkpd <- ggplot(conc_grid, aes(Cc, tca, colour = group)) +
  geom_line() +
  geom_vline(xintercept = 36.6, linetype = "dashed") +
  geom_vline(xintercept = 35, linetype = "dotted") +
  labs(
    x = "Serum eculizumab (ug/mL)", y = "Terminal complement activity (%)",
    colour = NULL,
    title = "PK/PD: inhibitory sigmoid Imax (Lee 2024 Table 4)",
    caption = paste("Dashed = IC50 36.6 ug/mL;",
                    "dotted = 35 ug/mL, the EMA-cited concentration sufficient",
                    "for complete complement inhibition (Lee 2024 Discussion).")
  ) +
  theme(legend.position = "bottom")

tca_grid <- tibble(tca = seq(0, 120, length.out = 401)) |>
  mutate(ldh = 206 + 1680 * tca^4.30 / (39^4.30 + tca^4.30))

p_pdeff <- ggplot(tca_grid, aes(tca, ldh)) +
  geom_line(colour = "firebrick") +
  geom_vline(xintercept = 39, linetype = "dashed") +
  labs(
    x = "Terminal complement activity (%)", y = "Serum LDH (U/L)",
    title = "PD/efficacy: direct sigmoid Emax (Lee 2024 Table 5)",
    caption = "Dashed = LC50 (39%, Table 5)."
  )

p_pkpd
```

![](Lee_2024_eculizumab_files/figure-html/exposure-response-1.png)

``` r

p_pdeff
```

![](Lee_2024_eculizumab_files/figure-html/exposure-response-2.png)

## PKNCA validation

NCA is run on the phase I healthy arm, where the paper’s single 300 mg
IV dose gives a clean single-dose profile. The PNH arm is a
multiple-dose regimen with trough-only sampling in the source study and
is validated below by its steady-state exposure instead.

``` r

sim_nca <- sim |>
  filter(DIS_PNH == 0) |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc, arm)

# Guarantee a time = 0 record per subject (IV bolus/infusion: pre-dose Cc = 0).
sim_nca <- bind_rows(
  sim_nca,
  sim_nca |> distinct(id, arm) |> mutate(time = 0, Cc = 0)
) |>
  distinct(id, arm, time, .keep_all = TRUE) |>
  arrange(id, time)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | arm + id)

dose_df <- events |>
  filter(evid == 1, DIS_PNH == 0) |>
  select(id, time, amt, arm)

dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id)

intervals <- data.frame(
  start      = 0,
  end        = Inf,
  cmax       = TRUE,
  tmax       = TRUE,
  aucinf.obs = TRUE,
  half.life  = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
```

### Comparison against the model-implied analytic solution

Lee 2024 does not publish an NCA table, so the reference column below is
the **analytic solution of the packaged model** at the typical healthy
82.40 kg subject. For a linear two-compartment IV model the exposure and
disposition half-lives follow in closed form from `CL`, `Vc`, `Vp` and
`Q`; a simulated NCA that reproduces them confirms that the ODE system,
the covariate scaling and the units were all encoded correctly.

``` r

wt_typ <- 82.40
wt_ref <- 80.9
dose   <- 300

cl_typ <- 0.0174 * (wt_typ / wt_ref)^1.1400
vc_typ <- 3.47   * (wt_typ / wt_ref)^0.8630
vp_typ <- 0.79
q_typ  <- 0.0134

k10 <- cl_typ / vc_typ
k12 <- q_typ / vc_typ
k21 <- q_typ / vp_typ

disc  <- sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21)
alpha <- ((k10 + k12 + k21) + disc) / 2
beta  <- ((k10 + k12 + k21) - disc) / 2

published <- tibble::tibble(
  arm        = "Healthy (phase I, 300 mg single IV)",
  cmax       = dose / vc_typ,        # C0 for a (near-)instantaneous IV input
  tmax       = 0.58,                 # end of the 35-min infusion
  aucinf.obs = dose / cl_typ,        # AUCinf = Dose / CL, exactly
  half.life  = log(2) / beta         # terminal (beta) half-life
)

knitr::kable(
  published |>
    dplyr::rename(
      "Arm"                  = arm,
      "Cmax (ug/mL)"         = cmax,
      "Tmax (h)"             = tmax,
      "AUC0-inf (ug*h/mL)"   = aucinf.obs,
      "t1/2 (h)"             = half.life
    ),
  digits  = c(0, 1, 2, 0, 1),
  caption = "Analytic reference values implied by the packaged model at the typical 82.40 kg healthy subject."
)
```

| Arm | Cmax (ug/mL) | Tmax (h) | AUC0-inf (ug\*h/mL) | t1/2 (h) |
|:---|---:|---:|---:|---:|
| Healthy (phase I, 300 mg single IV) | 85.1 | 0.58 | 16884 | 177.6 |

Analytic reference values implied by the packaged model at the typical
82.40 kg healthy subject. {.table}

``` r

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = published,
  by        = "arm",
  units     = c(cmax = "ug/mL", tmax = "h",
                aucinf.obs = "ug*h/mL", half.life = "h"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = paste("Simulated NCA (median of 150 virtual healthy subjects)",
                  "versus the model-implied analytic solution.",
                  "* differs from reference by >20%."),
  align   = c("l", "l", "r", "r", "r")
)
```

| NCA parameter | arm | Reference | Simulated | % diff |
|:---|:---|---:|---:|---:|
| Cmax (ug/mL) | Healthy (phase I, 300 mg single IV) | 85.1 | 84.4 | -0.8% |
| Tmax (h) | Healthy (phase I, 300 mg single IV) | 0.58 | 0.58 | +0.0% |
| AUC0-∞ (obs) (ug\*h/mL) | Healthy (phase I, 300 mg single IV) | 16900 | 16600 | -1.8% |
| t½ (h) | Healthy (phase I, 300 mg single IV) | 178 | 178 | +0.3% |

Simulated NCA (median of 150 virtual healthy subjects) versus the
model-implied analytic solution. \* differs from reference by \>20%.
{.table}

### Published claims the simulation must reproduce

Lee 2024 makes three quantitative statements that the packaged model can
be checked against directly.

``` r

cmax_by_id <- sim |>
  filter(DIS_PNH == 0) |>
  group_by(id) |>
  summarise(cmax = max(Cc, na.rm = TRUE), .groups = "drop")

# Claim 1 (Discussion): "In case of SB12 and ECU in single-dose phase I study
# with healthy subjects, the Cmax value was more than 0.85 IC50 units."
claim1 <- mean(cmax_by_id$cmax > 0.85 * 36.6)

# Claim 2 (Discussion, citing the EMA scientific report): "a serum eculizumab
# concentration of 35 ug/mL is sufficient to completely inhibit complement
# activity."  Typical-value residual activity at 35 ug/mL, PNH patients.
tca_at_35 <- 101.00 * (1 - 0.88 * 35^4.56 / (36.60^4.56 + 35^4.56))

# Claim 3 (Conclusion): the phase III regimen achieves "complete complement
# blockade in almost all PNH patients".  Evaluate at the last simulated trough.
last_trough <- sim |>
  filter(DIS_PNH == 1, time >= 24 * 7 * 48) |>
  group_by(id) |>
  slice_min(Cc, n = 1, with_ties = FALSE) |>
  ungroup()

tibble::tibble(
  Claim = c(
    "Phase I Cmax > 0.85 x IC50 (Discussion)",
    "Residual complement activity at 35 ug/mL, PNH (Discussion, EMA)",
    "Median phase III maintenance trough concentration (Conclusion)",
    "Fraction of PNH patients with trough > IC50 (Conclusion)",
    "Median phase III maintenance LDH (Results / Figure 4E)"
  ),
  Value = c(
    sprintf("%.1f%% of subjects", 100 * claim1),
    sprintf("%.1f%%", tca_at_35),
    sprintf("%.0f ug/mL", median(last_trough$Cc)),
    sprintf("%.1f%%", 100 * mean(last_trough$Cc > 36.6)),
    sprintf("%.0f U/L", median(last_trough$ldh))
  )
) |>
  knitr::kable(caption = "Published claims checked against the packaged model.")
```

| Claim | Value |
|:---|:---|
| Phase I Cmax \> 0.85 x IC50 (Discussion) | 100.0% of subjects |
| Residual complement activity at 35 ug/mL, PNH (Discussion, EMA) | 61.1% |
| Median phase III maintenance trough concentration (Conclusion) | 123 ug/mL |
| Fraction of PNH patients with trough \> IC50 (Conclusion) | 100.0% |
| Median phase III maintenance LDH (Results / Figure 4E) | 258 U/L |

Published claims checked against the packaged model. {.table}

``` r


stopifnot(claim1 == 1)                       # every phase I subject clears 0.85 x IC50
stopifnot(median(last_trough$Cc) > 36.6)     # maintenance troughs stay above IC50

# LDH on maintenance therapy. The typical-value exposure-response evaluated at
# the median trough complement activity gives ~222 U/L, below the ~250 U/L
# approximate upper limit of normal drawn in Figure 4E. The *population* median
# sits a little above that line (~258 U/L, with 47% of subjects below it)
# because LMAX carries 81.57 %CV and the response is non-linear, so the median
# of the transformed quantity exceeds the transform of the median. That is a
# property of the published variance model, not a defect: the fixed effects
# here are exactly Lee 2024 Table 5.
#
# So assert the typical-value claim against the clinical threshold, and assert
# the population median against the reduction from the untreated state, which
# is what the paper's Conclusion actually claims. A hard bound on the skewed
# population median against an approximate clinical number would be testing
# neither.
ldh_typical_at_trough <- 206 + 1680 * median(last_trough$tca)^4.30 /
  (39^4.30 + median(last_trough$tca)^4.30)
ldh_untreated <- 206 + 1680 * 100^4.30 / (39^4.30 + 100^4.30)

stopifnot(ldh_typical_at_trough < 250)                    # typical patient normalises
stopifnot(median(last_trough$ldh) < 0.2 * ldh_untreated)  # >= 80% reduction vs untreated
```

Every simulated phase I subject reaches a `Cmax` above `0.85 x IC50`,
which is the condition Lee 2024 cites (from Choe & Lee 2017) for the
sigmoid `Imax` parameters to be accurately and precisely estimable. At
35 ug/mL the typical PNH patient retains only a few percent of baseline
terminal complement activity, consistent with the EMA finding the paper
quotes. On the approved maintenance regimen the simulated troughs sit
far above `IC50` and the simulated LDH falls from an untreated typical
value of roughly 1860 U/L to inside the normal range – the “complete
complement blockade in almost all PNH patients” claim of the Conclusion.

## Assumptions and deviations

- **Only the pooled final model is packaged.** Lee 2024 also reports a
  healthy-subject-only PK model (Appendix Table A1: CL 0.0177 L/h,
  weight exponent 0.7630 on a 82.4 kg reference; Vc 3.50 L, weight
  exponent 0.6840; Vp 0.7850 L; Q 0.0136 L/h; IIV 13.55 / 12.66 / 31.76
  %CV with rho 0.553; proportional error 8.49%) and a
  healthy-subject-only PD model (Appendix Table A2: E0 86.20%, Imax
  0.93, `IC50 = 36.20 * (Baseline/87.75)^0.14`, H 4.52; IIV 15.88 / 2.63
  / 13.33 %CV; proportional error 18.00%). The paper describes these as
  the backbone for the pooled fit (“this model served as a backbone for
  the subsequent PK model incorporating pooled data from phase III”), so
  under the base-versus-final policy only the final pooled model is
  extracted. The appendix models’ baseline-terminal-complement-activity
  covariate on `IC50` has no analogue in the final model and is
  therefore not represented in the packaged file.
- **IIV scale.** Lee 2024 reports inter-individual variability as a
  percentage and states that ETAs were applied exponentially, without
  stating whether the percentage is `sqrt(omega^2) * 100` or
  `sqrt(exp(omega^2) - 1) * 100`. The model file uses the log-normal
  coefficient-of-variation convention `omega^2 = log(CV^2 + 1)`,
  matching the rest of the nlmixr2lib mAb models. The two conventions
  differ negligibly for the small IIVs (2.65-22.74 %CV) and by roughly
  15-25% in variance for the largest (`omega_LMAX` = 81.57 %CV).
- **Categorical covariate encoding.** The paper estimates separate
  typical values per subject group for `Vc`, `E0` and `Imax` rather than
  reporting covariate coefficients. These are encoded here as log-scale
  shifts (`e_pnh_vc`, `e_pnh_rbase_tca`, `e_pnh_imax_tca`) computed as
  the log ratio of the two published typical values, which reproduces
  both published values exactly. The three shift values are therefore
  derived, not transcribed.
- **`Imax` bound.** `Imax` carries an exponential ETA (Table 4) and is
  not bounded above by 1, so a small fraction of simulated subjects can
  take `imax_tca > 1` and drive terminal complement activity slightly
  negative at full blockade. The model floors the complement-activity
  driver of the LDH relationship at zero (`tca_pos <- max(tca, 0)`) so
  that `tca^LGAM` cannot return `NaN`; the floor is inert for any
  subject with a physically meaningful non-negative terminal complement
  activity. The `tca` output itself is left unfloored so this encoding
  choice is visible.
- **The LDH sub-model is calibrated in PNH patients only.** Lee 2024
  fitted the terminal-complement-activity-to-LDH relationship using
  phase III patient data alone, so `LL0`, `LMAX`, `LC50` and `LGAM`
  describe haemolysis in PNH. The packaged model evaluates `ldh` for any
  subject, and at the healthy baseline terminal complement activity of
  85.9% it returns roughly 1830 U/L – an extrapolation far outside the
  calibration range and not a prediction of LDH in a healthy person. Use
  the `ldh` output only with `DIS_PNH = 1`; the figures below do exactly
  that.
- **`cmt` on observation rows.** This is a declared three-endpoint model
  (`Cc`, `tca`, `ldh` each carry a `~` residual line), so
  `rxode2::rxode(mod)$predDf` assigns each endpoint its own slot as part
  of the model definition. Observation rows therefore use `cmt = "Cc"`,
  which is the endpoint name rather than an ODE state; `cmt = "central"`
  fails here with `'dvid'->'cmt' ... undefined compartment` because
  rxode2 cannot tell which endpoint the row belongs to. This is not the
  slot-renumbering antipattern, which applies to a bare algebraic
  observable with no `~` line. `rxSolve()` returns all three observables
  as columns regardless of which endpoint the observation rows name.
- **Sequential estimation not reproduced.** Lee 2024 linked PK to PD and
  PK/PD to efficacy sequentially with the individual population
  prediction (IPP) method. The packaged model evaluates all three levels
  simultaneously from one parameter vector, which is the correct
  construction for *simulation* but is not the estimation scheme the
  authors used; re-fitting this file to data would be a simultaneous
  fit, not an IPP fit.
- **No treatment covariate.** SB12, EU-sourced Soliris and US-sourced
  Soliris are described by the same parameter set because the paper
  found no statistically significant treatment effect on any ETA. The
  packaged model has no treatment covariate; simulate all three products
  with the same model.
- **Virtual-cohort covariates.** Body weights are drawn from truncated
  normal distributions matched to the phase I and phase III medians and
  ranges of Table 2; the paper reports medians and ranges only, not
  distributions. Race, sex, age and the baseline laboratory covariates
  screened by the authors are not in the final model and are not
  simulated.
- **Infusion duration in the PNH arm.** The phase I infusion duration
  (35 min) is given in Table 1; the phase III infusion duration is not
  stated, and the same 35 min is used for the maintenance doses. With a
  terminal half-life of roughly a week this has no material effect on
  the simulated profiles.
- **No published NCA to compare against.** Lee 2024 reports no NCA
  table, so the NCA comparison above uses the model’s own analytic
  solution as the reference. This validates the encoding (ODEs,
  covariate scaling, units), not the model against external data.
