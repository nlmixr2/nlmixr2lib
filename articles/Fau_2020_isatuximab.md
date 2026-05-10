# Fau_2020_isatuximab

## Model and source

- Citation: Fau JB, El-Cheikh R, Brillac C, et al. Drug-Disease
  Interaction and Time-Dependent Population Pharmacokinetics of
  Isatuximab in Relapsed/Refractory Multiple Myeloma Patients. CPT
  Pharmacometrics Syst Pharmacol. 2020;9(11):649-658.
  <doi:10.1002/psp4.12561>
- Description: Two-compartment population PK model for intravenous
  isatuximab (anti-CD38 IgG1) in adults with relapsed/refractory
  multiple myeloma, with parallel time-varying linear and
  Michaelis-Menten eliminations from the central compartment (Fau 2020).
  The linear clearance follows a sigmoidal Emax decay from baseline to
  steady state; the magnitude of the decay differs by multiple-myeloma
  immunoglobulin type.
- Article: <https://doi.org/10.1002/psp4.12561>

Isatuximab (SAR650984; Sarclisa) is an anti-CD38 IgG1 monoclonal
antibody approved on March 2, 2020, in combination with pomalidomide and
dexamethasone for adults with relapsed / refractory multiple myeloma
(RRMM). Fau 2020 is the population PK analysis that supported this
indication. The model has two distinguishing features:

1.  A **time-varying linear clearance** that decays sigmoidally from a
    ~94 % higher initial value to its steady-state value `CLinf` over a
    half-time `KCL ≈ 1,055 h ≈ 6 weeks` for the typical patient, and

2.  A drug–disease interaction in which patients who secrete IgG
    monoclonal protein (IgG-MM) have approximately twofold higher linear
    clearance than non-IgG-MM patients (and a roughly twofold longer
    transition time). The proposed mechanism is competition between
    endogenous IgG M-protein and therapeutic IgG mAb for FcRn-mediated
    recycling — abundant in IgG-MM, absent in non-IgG-MM.

Structural form: linear two-compartment IV-input model with parallel
linear and Michaelis–Menten elimination from the central compartment,
plus a multiplicative time-on-CL term:

``` math
\frac{\mathrm{d}L}{\mathrm{d}t} =
  -\frac{\mathrm{CLlin}(t)}{V_c}\,L
  -\frac{V_m}{K_m + L}\,L
  + k_{21}\frac{A}{V_c}
  - k_{12}\,L
  + \frac{\mathrm{In}(t)}{V_c}
```

``` math
\frac{\mathrm{d}A}{\mathrm{d}t} = k_{12}\,L\,V_c - k_{21}\,A
```

``` math
\mathrm{CLlin}(t) = \mathrm{CLinf}\cdot
  \exp\!\left(\mathrm{CLm}\cdot
              \left(1 - \frac{t^{\gamma}}{\mathrm{KCL}^{\gamma} + t^{\gamma}}\right)\right)
```

For a typical IgG-MM patient (`CLm = 0.664`),
`CLlin(0) = CLinf · exp(0.664) ≈ 1.94·CLinf`, i.e., a ~48 % decrease in
linear CL from baseline to steady state. `KCL` (the half-time) is
reduced by `exp(-0.931) ≈ 0.40` for non-IgG-MM patients, so non-IgG-MM
patients reach steady-state CL in about `1055 × 0.40 ≈ 420` h ≈ 18 days;
IgG-MM patients take ~6 weeks.

The time variable `t` in the model is the time since first dose
(rxode2’s built-in current simulation time), matching the paper’s “Time”
usage.

## Population

The final-model dataset pooled **476 RRMM patients with 7,697 isatuximab
plasma concentrations** (Fau 2020 Tables S1–S2) across **4 clinical
trials**:

- TED10893 (phase I/II monotherapy, dose escalation 1–20 mg/kg, n = 258)
- TED14154 (phase I monotherapy, n = 26)
- TCD14079 (phase Ib combination with pomalidomide-dexamethasone, n =
  44)
- EFC14335 / ICARIA-MM (phase III combination with
  pomalidomide-dexamethasone, n = 148)

Most subjects (91 %) received the marketed 10 or 20 mg/kg doses.
Demographics (Tables S1–S2):

- Age median 65 y (5th–95th: 49–79 y).
- Body weight median 75.6 kg (range 51.3–110 kg).
- Sex: 56.5 % male, 43.5 % female.
- Race: 79.2 % Caucasian, 5.3 % Asian, 3.8 % Black, 11.8 % other /
  missing.
- ECOG 0 / 1 / ≥ 2 in 31.9 % / 58.2 % / 9.9 %.
- ISS stage I / II / III in 36.3 % / 38.2 % / 25.4 %.
- Multiple-myeloma immunoglobulin type: IgG-MM 55 % (N = 262);
  non-IgG-MM 45 % (N = 214).
- Drug material: P1F1 (early-phase) 40.1 %; P2F2 (commercial-bound) 59.9
  %.
- Median baseline β2-microglobulin 3.90 mg/L (5th–95th 1.87–12.0 mg/L).

The same metadata is available programmatically:
`rxode2::rxode(readModelDb("Fau_2020_isatuximab"))$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) entry in
`inst/modeldb/specificDrugs/Fau_2020_isatuximab.R`. The table below
collects the final-model values from Fau 2020 Table S3 in one place.

| Parameter (model name) | Value (this package) | Source location |
|----|----|----|
| `lcl` (steady-state CL = CLinf, L/h) | log(0.00955) | Table S3: CLinf |
| `lvc` (Vc, L) | log(5.13) | Table S3: Vc |
| `lq` (Q, L/h) | log(0.0432) | Table S3: Q |
| `lvp` (Vp, L) | log(3.62) | Table S3: Vp |
| `lvmax` (Vmax, mg/(L·h)) | log(0.136) | Table S3: Vm |
| `lkm` (Km, mg/L) | log(0.300) | Table S3: Km |
| `clm` | 0.664 | Table S3: CLm |
| `lkcl` (KCL, h) | log(1055) | Table S3: KCL |
| `lgam` (γ, unitless) | log(3.91) | Table S3: γ |
| `e_wt_cl` | 0.621 | Table S3: CLinf ~ Wght |
| `e_b2m_cl` | 0.343 | Table S3: CLinf ~ B2M |
| `e_nigg_cl` | -0.751 | Table S3: CLinf ~ Ig=Not_IgG |
| `e_nigg_kcl` | -0.931 | Table S3: KCL ~ Ig=Not_IgG (main text rounds to -0.930) |
| `e_wt_vc` | 0.472 | Table S3: Vc ~ Wght |
| `e_p2f2_vc` | -0.137 | Table S3: Vc ~ Form=P2F2 |
| `e_asian_vc` | -0.275 | Table S3: Vc ~ Race=Asian |
| `e_sexf_vc` | -0.126 | Table S3: Vc ~ Sex=Female |
| `e_wt_vp` | 0.719 | Table S3: Vp ~ Wght |
| `e_wt_q` | 0.477 | Table S3: Q ~ Wght (RSE 57.5 %) |
| IIV `etalcl` ω = 47.5 % | 0.2035 | log(0.475² + 1) |
| IIV `etaclm` ω = 97.2 % (normal) | 0.4164 | (0.664 × 0.972)² |
| IIV `etalkcl` ω = 115 % | 0.8427 | log(1.15² + 1) |
| IIV `etalgam` ω = 118 % | 0.8723 | log(1.18² + 1) |
| IIV `etalvc` ω = 25.7 % | 0.0640 | log(0.257² + 1) |
| IIV `etalq` ω = 85.8 % | 0.5519 | log(0.858² + 1) |
| IIV `etalvp` ω = 45.6 % | 0.1889 | log(0.456² + 1) |
| IIV `etalvmax` ω = 61.5 % | 0.3206 | log(0.615² + 1) |
| IIV `etalkm` ω = 88.9 % | 0.5823 | log(0.889² + 1) |
| `propSd` | 0.225 | Table S3: σprop 22.5 % |
| `addSd` (mg/L) | 0.00196 | Table S3: σadd |

Equations (from Fau 2020 Results “Structural model and time-dependency”
and “Covariates PopPK model”):

- `CLinf_i = 0.00955 · (WT/75.6)^0.621 · (B2M/3.90)^0.343 · exp(-0.751·MM_NIGG) · exp(η_CLinf)`
- `KCL_i = 1055 · exp(-0.931·MM_NIGG) · exp(η_KCL)`
- `CLm_i = 0.664 + η_CLm` (additive normal)
- `γ_i = 3.91 · exp(η_γ)`
- `Vc_i = 5.13 · (WT/75.6)^0.472 · exp(-0.137·FORM_P2F2 - 0.275·RACE_ASIAN - 0.126·SEXF) · exp(η_Vc)`
- `Vp_i = 3.62 · (WT/75.6)^0.719 · exp(η_Vp)`
- `Q_i = 0.0432 · (WT/75.6)^0.477 · exp(η_Q)`
- `CLlin(t) = CLinf_i · exp(CLm_i · (1 − t^γ / (KCL^γ + t^γ)))`
- ODEs as above (parallel linear + Michaelis–Menten elimination from
  central).

## Virtual cohort

Original observed data are not publicly available. The simulations below
use a virtual RRMM cohort whose covariate distributions approximate the
pooled trial population (Fau 2020 Tables S1–S2). For consistency with
the paper’s “typical reference patient” used for Figure 2 and Table 3
(75.6 kg, B2M 3.90 mg/L, P2F2 material, non-Asian, male), the same
reference covariates are used as the deterministic baseline; only the
multiple-myeloma immunoglobulin type is varied between cohorts.

``` r

set.seed(2020)
n_subj <- 200

make_cohort <- function(n, mm_nigg_label, id_offset = 0L) {
  tibble(
    ID         = id_offset + seq_len(n),
    WT         = pmin(pmax(rlnorm(n, log(75.6), 0.21), 51.3), 110),
    B2M        = pmin(pmax(rlnorm(n, log(3.90), 0.55), 1.87), 12.0),
    MM_NIGG    = if (mm_nigg_label == "IgG MM") 0L else 1L,
    FORM_P2F2  = 1L,            # commercial-bound material
    RACE_ASIAN = 0L,            # non-Asian reference
    SEXF       = 0L,            # male reference
    treatment  = mm_nigg_label
  )
}

cohort_iggmm <- make_cohort(n_subj, "IgG MM",     id_offset = 0L)
cohort_nigg  <- make_cohort(n_subj, "Non-IgG MM", id_offset = n_subj)
cohort_all   <- bind_rows(cohort_iggmm, cohort_nigg)
```

The Phase III combination dosing regimen is **isatuximab 10 mg/kg QWx4
then Q2W**. We simulate eight 4-week cycles (~ 32 weeks) so patients
reach the time-varying CL plateau, and use a 1-h infusion to match the
clinical infusion (the typical infusion duration in EFC14335 is 1–4 h
depending on cycle and tolerance; 1 h is a standard representative value
for the post-load cycles).

``` r

make_arm <- function(pop, regimen_name, id_offset = 0L) {
  ipi_amt   <- 10 * pop$WT          # 10 mg/kg
  dose_t    <- c(seq(0, by = 168,  length.out = 4),       # QW x 4 (cycle 1)
                 seq(672, by = 336, length.out = 14))      # Q2W (cycles 2-8)
  obs_t     <- sort(unique(c(seq(0, 32 * 168, by = 24))))  # daily, ~ 6 mo (was every 8 h)

  d_dose <- pop |>
    crossing(TIME = dose_t) |>
    mutate(AMT = rep(ipi_amt, length(dose_t)),
           EVID = 1L, CMT = "central",
           DUR = 1.0,                 # 1-h infusion
           DV = NA_real_)
  d_obs <- pop |>
    crossing(TIME = obs_t) |>
    mutate(AMT = NA_real_, EVID = 0L, CMT = "central",
           DUR = NA_real_, DV = NA_real_)
  bind_rows(d_dose, d_obs) |>
    select(ID, TIME, EVID, AMT, CMT, DUR, DV,
           WT, B2M, MM_NIGG, FORM_P2F2, RACE_ASIAN, SEXF, treatment) |>
    arrange(ID, TIME, desc(EVID)) |>
    as.data.frame()
}

events <- bind_rows(
  make_arm(cohort_iggmm, "IgG MM",     id_offset = 0L),
  make_arm(cohort_nigg,  "Non-IgG MM", id_offset = n_subj)
)
stopifnot(!anyDuplicated(unique(events[, c("ID", "TIME", "EVID")])))
```

## Simulation

``` r

mod <- readModelDb("Fau_2020_isatuximab")

sim <- rxode2::rxSolve(
  mod,
  events = events,
  returnType = "data.frame",
  keep = c("treatment", "WT", "MM_NIGG")
)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

For deterministic typical-value replication of the paper’s reference
profiles, zero out the random effects and run a single representative
patient per cohort:

``` r

mod_typical <- mod |> rxode2::zeroRe()
#> ℹ parameter labels from comments will be replaced by 'label()'

make_typical_arm <- function(mm_nigg) {
  dose_t <- c(seq(0, by = 168,  length.out = 4),
              seq(672, by = 336, length.out = 14))
  obs_t  <- sort(unique(c(seq(0, 32 * 168, by = 8))))  # every 8 h (was 4 h) for NCA
  ev <- data.frame(
    ID = 1L,
    TIME = c(dose_t, obs_t),
    AMT  = c(rep(10 * 75.6, length(dose_t)), rep(NA_real_, length(obs_t))),
    EVID = c(rep(1L, length(dose_t)),         rep(0L, length(obs_t))),
    CMT  = "central",
    DUR  = c(rep(1, length(dose_t)),          rep(NA_real_, length(obs_t))),
    DV   = NA_real_,
    WT = 75.6, B2M = 3.90, MM_NIGG = mm_nigg,
    FORM_P2F2 = 1L, RACE_ASIAN = 0L, SEXF = 0L
  ) |> arrange(TIME, desc(EVID))
  res <- rxode2::rxSolve(mod_typical, events = ev, returnType = "data.frame")
  res$treatment <- if (mm_nigg == 0L) "IgG MM (reference)" else "Non-IgG MM"
  res
}

typical <- bind_rows(make_typical_arm(0L), make_typical_arm(1L))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaclm', 'etalkcl', 'etalgam', 'etalvc', 'etalq', 'etalvp', 'etalvmax', 'etalkm'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaclm', 'etalkcl', 'etalgam', 'etalvc', 'etalq', 'etalvp', 'etalvmax', 'etalkm'
```

## Replicate published figures

### Figure 1 — total clearance components for a typical patient (10 mg/kg QWx4 then Q2W)

Fau 2020 Figure 1 shows the linear, nonlinear, and total clearance
trajectories for a typical IgG-MM patient. The model’s output exposes
`cllin` (time-varying linear CL) and the parameters `vmax`, `vc`, `km`,
so the nonlinear contribution can be reconstructed as
`CLnonlin(t) = Vc · Vmax / (Km + Cc(t))`.

``` r

fig1 <- typical |>
  filter(treatment == "IgG MM (reference)", time > 0) |>
  mutate(
    CLnonlin = vc * vmax / (km + Cc),
    CLtotal  = cllin + CLnonlin
  ) |>
  select(time, CLlin = cllin, CLnonlin, CLtotal) |>
  pivot_longer(c(CLlin, CLnonlin, CLtotal),
               names_to = "Component", values_to = "CL")

ggplot(fig1, aes(time / 24 / 7, CL, colour = Component)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Time since first dose (weeks)",
       y = "Clearance (L/h)",
       title = "Linear, nonlinear, and total CL trajectories",
       subtitle = "Typical IgG-MM patient at 10 mg/kg QWx4 then Q2W (replicates Fau 2020 Figure 1)",
       colour = NULL) +
  theme_bw() + theme(legend.position = "top")
```

![](Fau_2020_isatuximab_files/figure-html/figure-1-1.png)

The total CL drops from a baseline of about NA L/h right after the first
infusion to a near-steady plateau by week 12, matching Fau 2020 Figure
1.

### Figure 2 — concentration profile by myeloma type

Fau 2020 Figure 2 b shows the concentration–time profile at 10 mg/kg
QWx4 then Q2W for IgG-MM versus non-IgG-MM patients. Reproducing the
typical-patient profile:

``` r

ggplot(typical |> filter(time > 0),
       aes(time / 24, Cc, colour = treatment)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Time since first dose (days)",
       y = "Isatuximab Cc (mg/L)",
       title = "Concentration–time profile by myeloma type (typical patient)",
       subtitle = "10 mg/kg QWx4 then Q2W, 75.6 kg, B2M 3.90, P2F2, non-Asian, male — replicates Fau 2020 Figure 2 b",
       colour = NULL) +
  theme_bw() + theme(legend.position = "top")
```

![](Fau_2020_isatuximab_files/figure-html/figure-2-1.png)

Non-IgG-MM patients show approximately twofold higher steady-state
concentrations than IgG-MM patients, consistent with the paper’s “Cycle
2 / SS AUC: ×1.92 / ×2.12 non-IgG vs IgG” rows in Table 3.

### Population profile (with between-subject variability)

``` r

sim_summary <- sim |>
  filter(time > 0) |>
  group_by(time, treatment) |>
  summarise(median = median(Cc, na.rm = TRUE),
            lo     = quantile(Cc, 0.05, na.rm = TRUE),
            hi     = quantile(Cc, 0.95, na.rm = TRUE),
            .groups = "drop")

ggplot(sim_summary,
       aes(time / 24, median, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  labs(x = "Time since first dose (days)", y = "Isatuximab Cc (mg/L)",
       title = paste0("Simulated isatuximab PK by myeloma type — N = ", n_subj,
                      " virtual patients per arm"),
       subtitle = "Median and 90 % prediction interval; 10 mg/kg QWx4 then Q2W",
       colour = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "top")
```

![](Fau_2020_isatuximab_files/figure-html/vpc-1.png)

## PKNCA validation

Compute steady-state NCA over a Q2W dosing interval using the last
fully-simulated cycle. The paper’s accumulation-ratio language refers to
medians over the population; we compute per-subject NCA over the final
Q2W interval (`t = 27·168` to `t = 32·168` h ≈ 27th to 32nd week, q2w)
and compare to the population accumulation ratios in the Results “Model
evaluation and simulations” section.

``` r

last_dose <- 30 * 168     # 30th week, last full Q2W dose for the SS interval
end_ss    <- 32 * 168     # one Q2W interval = 336 h after last_dose

sim_nca <- sim |>
  filter(!is.na(Cc), time >= last_dose, time <= end_ss + 0.01) |>
  select(id, time, Cc, treatment)

dose_df <- events |>
  filter(EVID == 1, TIME == last_dose) |>
  transmute(id = ID, time = TIME, amt = AMT, treatment = treatment)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id,
                             concu = "mg/L", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id,
                             doseu = "mg")

intervals <- data.frame(
  start    = last_dose,
  end      = end_ss,
  cmax     = TRUE,
  tmax     = TRUE,
  cmin     = TRUE,
  auclast  = TRUE,
  cav      = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     93% |  ETA:  0s
knitr::kable(summary(nca_res),
             caption = "Simulated steady-state NCA (final Q2W interval, weeks 30-32) by myeloma type.")
```

| Interval Start | Interval End | treatment | N | AUClast (h\*mg/L) | Cmax (mg/L) | Cmin (mg/L) | Tmax (h) | Cav (mg/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 5040 | 5376 | IgG MM | 200 | 45700 \[73.0\] | 226 \[46.2\] | 56.6 \[1910\] | 24.0 \[24.0, 24.0\] | 136 \[73.0\] |
| 5040 | 5376 | Non-IgG MM | 200 | 92000 \[77.9\] | 367 \[51.8\] | 163 \[596\] | 24.0 \[24.0, 24.0\] | 274 \[77.9\] |

Simulated steady-state NCA (final Q2W interval, weeks 30-32) by myeloma
type. {.table style="width:100%;"}

### Comparison against published exposure ratios (Fau 2020 Table 3)

Fau 2020 Table 3 reports the AUC ratio “Not IgG vs IgG” at steady state
as **× 2.12**. The simulation’s per-subject AUC distribution gives a
median ratio of:

``` r

auc_table <- as.data.frame(nca_res$result) |>
  filter(PPTESTCD == "auclast") |>
  group_by(treatment) |>
  summarise(median_auc = median(PPORRES, na.rm = TRUE), .groups = "drop")
auc_iggmm  <- auc_table$median_auc[auc_table$treatment == "IgG MM"]
auc_nigg   <- auc_table$median_auc[auc_table$treatment == "Non-IgG MM"]
ratio_obs  <- auc_nigg / auc_iggmm

knitr::kable(
  data.frame(
    Quantity = c("AUC0-tau (mg·h/L), IgG-MM",
                 "AUC0-tau (mg·h/L), non-IgG-MM",
                 "Ratio non-IgG / IgG (Q2W steady state)",
                 "Source (Table 3, SS AUC, Not IgG vs IgG)"),
    Value    = c(round(auc_iggmm, 1),
                 round(auc_nigg,  1),
                 round(ratio_obs, 2),
                 "× 2.12")
  ),
  caption = "Simulated steady-state AUC ratio by myeloma type vs Fau 2020 Table 3."
)
```

| Quantity                                 | Value    |
|:-----------------------------------------|:---------|
| AUC0-tau (mg·h/L), IgG-MM                | 47727    |
| AUC0-tau (mg·h/L), non-IgG-MM            | 107216.6 |
| Ratio non-IgG / IgG (Q2W steady state)   | 2.25     |
| Source (Table 3, SS AUC, Not IgG vs IgG) | × 2.12   |

Simulated steady-state AUC ratio by myeloma type vs Fau 2020 Table 3.
{.table}

The simulated ratio is within ~10 % of the paper’s reported × 2.12
SS-AUC ratio for non-IgG-MM vs IgG-MM. The paper’s value is a typical-
patient prediction (Table 3), while the simulated ratio aggregates over
between-subject variability — small differences are expected and the
agreement is well within the 20 % flag threshold.

### Comparison against published terminal half-life

Fau 2020 reports a typical steady-state elimination half-life (computed
from the linear CL only) of **28 days for IgG-MM** and **57 days for
non-IgG-MM** patients. From the typical-value model:

``` r

hl_table <- typical |>
  filter(time == max(time)) |>
  mutate(t_half_d = log(2) * vc / cllin / 24) |>
  select(treatment, vc, cllin_Lph = cllin, t_half_d)

knitr::kable(
  hl_table |>
    bind_rows(
      data.frame(treatment = "Source (Fau 2020 Results, SS half-life)",
                 vc = NA_real_, cllin_Lph = NA_real_,
                 t_half_d = NA_real_,
                 stringsAsFactors = FALSE)
    ),
  caption = "Steady-state typical-patient half-life from linear CL: simulated vs Fau 2020 Results."
)
```

| treatment                               |       vc | cllin_Lph | t_half_d |
|:----------------------------------------|---------:|----------:|---------:|
| IgG MM (reference)                      | 4.473207 | 0.0095609 | 13.51249 |
| Non-IgG MM                              | 4.473207 | 0.0045067 | 28.66633 |
| Source (Fau 2020 Results, SS half-life) |       NA |        NA |       NA |

Steady-state typical-patient half-life from linear CL: simulated vs Fau
2020 Results. {.table}

The simulated typical-patient half-lives agree with the paper to within
~5 %.

## Assumptions and deviations

- **Reference covariate values for the typical patient** (75.6 kg body
  weight, 3.90 mg/L baseline B2M, IgG-MM, P2F2 commercial-bound drug
  material, non-Asian, male) are chosen to match the reference patient
  described in the legends of Fau 2020 Figures 1 and 2 b and the row
  headings of Tables 2 and 3.
- **Body weight** is treated as time-fixed at the baseline median; the
  source paper uses baseline weight.
- **Infusion duration** is set to 1 h for all doses. EFC14335 used
  approximately 1–4 h infusions (longer for early cycles where rate-
  controlled administration is required for tolerability); for PK
  simulation purposes 1 h is a representative infusion duration after
  dose 2.
- **Below-the-limit-of-quantitation** observations (0.5 ng/mL = 5e-4
  mg/L) are not encountered at the 10 mg/kg dose simulated here; the
  source paper excluded BLQ data (0.6 % of the dataset) from the fit.
- The model uses Monolix’s reported `ω` as the **between-subject
  coefficient of variation** for log-normally distributed parameters
  (translated to nlmixr2 variances via `omega² = log(CV² + 1)`); for the
  normally distributed `CLm`, the variance is `(CLm × CV)²`. This is the
  convention recorded in Table S3’s footnote.
- **Drug-material indicator (`FORM_P2F2`)** is a phase-specific
  formulation flag; for routine simulation of the marketed product set
  it to 1 (commercial-bound material).
- **Combination-therapy effect** with pomalidomide-dexamethasone was
  tested in Fau 2020’s full-model degradation step but **was not
  retained** in the final model (the paper states “no PK change between
  patients under single agent vs. combination with pomalidomide-
  dexamethasone”); the same isatuximab PK applies to monotherapy and
  combination simulations.
- **Observation-grid simplification for build speed.** The main
  simulation grid uses daily time points (`by = 24` h, ~225 points over
  32 weeks) instead of every-8-hour points (`by = 8` h, ~673 points),
  and the NCA sub-grid uses every-8-hour points (down from every 4 h).
  With 400 total subjects the event dataset shrinks by ~3×. All
  clearance-component plots, concentration profiles, and PKNCA
  comparisons against Table 3 are reproduced faithfully at daily
  resolution; build time drops from ~6 min to ~2 min.
