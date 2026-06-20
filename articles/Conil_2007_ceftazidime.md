# Ceftazidime (Conil 2007)

## Model and source

``` r

mod_meta <- nlmixr2est::nlmixr(readModelDb("Conil_2007_ceftazidime"))$meta
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Conil JM, Georges B, Lavit M, Laguerre J, Samii K, Houin G,
  Saivin S. A population pharmacokinetic approach to ceftazidime use in
  burn patients: influence of glomerular filtration, gender and
  mechanical ventilation. Br J Clin Pharmacol. 2007;64(1):27-35.
  <doi:10.1111/j.1365-2125.2007.02857.x>
- Description: Two-compartment IV population PK model for ceftazidime in
  adult burn-ICU patients, with creatinine clearance on CL and sex /
  mechanical ventilation / creatinine clearance on the peripheral volume
  V2 (Conil 2007)
- Article (DOI): <https://doi.org/10.1111/j.1365-2125.2007.02857.x>

This vignette validates the packaged `Conil_2007_ceftazidime` model – a
two-compartment IV population PK model for ceftazidime in adult burn-ICU
patients with creatinine clearance on CL and sex, mechanical
ventilation, and creatinine clearance on the peripheral volume V2 –
against the source publication’s Table 1 (baseline demographics), Table
2 (basic-model parameter estimates), and Table 3 (mean PK parameters in
the four sex / ventilation strata).

## Population

Conil 2007 studied 50 adult burn patients (38 male, 12 female; mean age
52.3 +/- 20.7 years; mean body weight 71.3 +/- 13.5 kg) at the
University Hospital of Toulouse-Rangueil Burns Unit over a 4-year period
(Conil 2007 Methods, “Subjects and sampling”; Table 1). Patients were
treated during the secondary phase of their burn injury for local
infection or sepsis. The cohort had a mean burned surface area of 23 +/-
13.5% of total body surface, UBS index 64.8 +/- 50.0, and Baux index
75.6 +/- 22.4. Mean serum creatinine was in the normal range (76.5 +/-
21.8 umol/L), mean creatinine clearance was 105.3 +/- 39.3 mL/min (range
33-191), and 16 patients (32%) required mechanical ventilation during
ceftazidime treatment (Table 1).

Each patient received 6 g/24 h of ceftazidime administered as either
three 2 g doses every 8 h or six 1 g doses every 4 h; each dose was
given as a 20-minute IV infusion via electric syringe. Plasma
ceftazidime was assayed by reversed-phase HPLC with UV detection at 260
nm (LLOQ 1 mg/L, calibration linear 1-100 mg/L); 237 concentrations were
available from the 50 patients (mean 4.7 samples per patient).

The same information is available programmatically via the model’s
`population` metadata:

``` r

str(mod_meta$population)
#> List of 14
#>  $ species       : chr "human"
#>  $ n_subjects    : int 50
#>  $ n_studies     : int 1
#>  $ age_range     : chr "15-90 years"
#>  $ age_median    : chr "52.3 years (mean +/- 20.7 SD)"
#>  $ weight_range  : chr "50-108 kg"
#>  $ weight_median : chr "71.3 kg (mean +/- 13.5 SD)"
#>  $ sex_female_pct: num 24
#>  $ race_ethnicity: chr "Not reported (French university-hospital burn-ICU population, Toulouse)"
#>  $ disease_state : chr "Burn injury during the secondary phase; treated for local infection or sepsis with ceftazidime; mean burned sur"| __truncated__
#>  $ dose_range    : chr "Ceftazidime 6 g/24 h administered as either three 2 g doses every 8 h or six 1 g doses every 4 h; each dose giv"| __truncated__
#>  $ regions       : chr "France (Burns Unit, University Hospital of Toulouse-Rangueil)"
#>  $ mech_vent_pct : num 32
#>  $ notes         : chr "Single-center observational study over 4 years, 50 patients with 237 serum ceftazidime concentrations (mean 4.7"| __truncated__
```

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in `inst/modeldb/specificDrugs/Conil_2007_ceftazidime.R`.
The table below collects them in one place. Final-model structural
equations and IIV / residual values come from Conil 2007 Results (p. 31,
“Final model” subsection); basic-model carry-forward values for Q IIV,
V2 IIV, and the proportional residual come from Table 2 and the “Basic
model” subsection of Results.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL intercept; non-renal CL) | log(1.08) | Conil 2007 Results p. 31 final-model eqn `CL = 1.08 + 0.0536 * CLCR` |
| `e_crcl_cl` (CL slope on CRCL) | 0.0536 L/h per mL/min | Conil 2007 Results p. 31 final-model eqn |
| `lvc` (V1) | log(18.81) | Conil 2007 Results p. 31 final-model eqn `V1 = 18.81 L`; Table 3 reports 18.8 in all four strata |
| `lvp` (V2 baseline scale) | log(2.69) | Conil 2007 Results p. 31 final-model eqn `TVV2 = 2.69 * (1 + 1.43*SEX) * (1 + 2.44*VENT) * (1 + 0.00414*CLCR)` |
| `e_sexf_vp` (SEXF fractional effect on V2) | 1.43 | Conil 2007 Results p. 31 final-model eqn |
| `e_mech_vent_vp` (MECH_VENT fractional effect on V2) | 2.44 | Conil 2007 Results p. 31 final-model eqn |
| `e_crcl_vp` (CRCL fractional effect on V2) | 0.00414 per mL/min | Conil 2007 Results p. 31 final-model eqn |
| `lq` (Q) | log(6.881) | Conil 2007 Results p. 31 final-model eqn `Q = 6.881 L/h`; Table 3 reports 6.9 in all four strata |
| `etalcl ~ 0.02528` | log(0.16^2 + 1) | Conil 2007 Results p. 31: “Interindividual variability in the clearance was decreased to 16%” |
| `etalvc ~ 0.01676` | log(0.13^2 + 1) | Conil 2007 Results p. 31: “and that of the central volume of distribution to 13%” |
| `etalq ~ 2.6520` | log(3.63^2 + 1) | Conil 2007 Table 2 basic-model row “Inter-compartmental clearance” CV 363% (final model silent; carried forward per operator standing instruction) |
| `etalvp ~ 1.5280` | log(1.90^2 + 1) | Conil 2007 Table 2 basic-model row “Distribution volume of the peripheral compartment” CV 190% (final model silent; carried forward) |
| `propSd <- 0.38` | 0.38 | Conil 2007 Results p. 30 basic-model proportional CV 38% (final model silent; carried forward) |
| `d/dt(central)` / `d/dt(peripheral1)` | n/a | Conil 2007 Results p. 30 (“A two-compartment model… ADVAN3 TRANS4”); standard 2-cmt linear ODE form |
| `Cc ~ prop(propSd)` | n/a | Conil 2007 Results p. 30 (“a proportional error model best fitted the data”) |

## Virtual cohort

The original observed ceftazidime concentrations are not publicly
available. The virtual cohort below reproduces the four sex /
mechanical- ventilation strata of Conil 2007 Table 3, using the
stratum-specific median CRCL the paper reports (males not ventilated 114
mL/min, females not ventilated 103, males ventilated 93, females
ventilated 87). Each subject receives a single 2 g IV infusion over 20
minutes (rate = 6000 mg/h), sampled densely enough to characterize Cmax
during the distribution phase and the terminal half-life out to 24 h.

``` r

set.seed(20260618)

# Per-stratum cohort definition matching Conil 2007 Table 3.
strata <- tibble::tribble(
  ~stratum,                 ~n,  ~SEXF, ~MECH_VENT, ~CRCL_median,
  "Males, not ventilated",  100L,    0,         0,          114,
  "Females, not ventilated", 100L,   1,         0,          103,
  "Males, ventilated",      100L,    0,         1,           93,
  "Females, ventilated",    100L,    1,         1,           87
)

# Slight CRCL variability around the per-stratum median so the cohort is
# heterogeneous enough for stable NCA percentile bands; SD 10 mL/min is
# tight relative to the cohort-wide range (33-191) but realistic for an
# already-homogenized stratum.
make_stratum <- function(s, idx_offset) {
  n <- s$n
  crcl <- pmax(20, rnorm(n, mean = s$CRCL_median, sd = 10))
  tibble::tibble(
    id        = idx_offset + seq_len(n),
    stratum   = s$stratum,
    SEXF      = s$SEXF,
    MECH_VENT = s$MECH_VENT,
    CRCL      = crcl
  )
}

cov_tab <- bind_rows(
  make_stratum(strata[1, ], 0L),
  make_stratum(strata[2, ], 100L),
  make_stratum(strata[3, ], 200L),
  make_stratum(strata[4, ], 300L)
)

# Dosing: 2 g IV over 20 min; sampling 0 (pre-dose) through 24 h spanning
# the distribution and terminal phases (cohort half-life ~ 2-8 h per
# stratum, so 24 h covers 3-12 terminal half-lives).
dose_amt    <- 2000          # mg
infusion_h  <- 20 / 60       # 20 minutes
dose_rate   <- dose_amt / infusion_h
sample_times <- c(0, 1/3, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 24)

make_subject <- function(row) {
  doses <- tibble::tibble(
    id   = row$id,        time = 0,
    evid = 1L,            amt  = dose_amt,
    rate = dose_rate,     dv   = NA_real_
  )
  obs <- tibble::tibble(
    id   = row$id,        time = sample_times,
    evid = 0L,            amt  = NA_real_,
    rate = NA_real_,      dv   = NA_real_
  )
  bind_rows(doses, obs) |>
    mutate(
      stratum   = row$stratum,
      SEXF      = row$SEXF,
      MECH_VENT = row$MECH_VENT,
      CRCL      = row$CRCL
    ) |>
    arrange(time, desc(evid))
}

events <- bind_rows(lapply(seq_len(nrow(cov_tab)), function(i) {
  make_subject(cov_tab[i, ])
}))

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod         <- readModelDb("Conil_2007_ceftazidime")
mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_stoch <- rxode2::rxSolve(
  object = mod, events = events,
  keep   = c("stratum", "SEXF", "MECH_VENT", "CRCL")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_typical <- rxode2::rxSolve(
  object = mod_typical, events = events,
  keep   = c("stratum", "SEXF", "MECH_VENT", "CRCL")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp'
#> Warning: multi-subject simulation without without 'omega'
```

## Replicate Table 3 – typical-value PK parameters by stratum

Conil 2007 Table 3 reports the typical-value CL, V1, Q, V2, total V,
elimination rate `kel = CL / V1`, and apparent half-life
`t1/2 = ln(2) / kel` for each of the four sex / ventilation strata
evaluated at the stratum-specific median CRCL. The packaged model
reproduces those values directly from the final-model equations:

``` r

stratum_typicals <- strata |>
  mutate(
    TVCL = 1.08 + 0.0536 * CRCL_median,
    TVV1 = 18.81,
    TVQ  = 6.881,
    TVV2 = 2.69 *
      (1 + 1.43    * SEXF)      *
      (1 + 2.44    * MECH_VENT) *
      (1 + 0.00414 * CRCL_median),
    TVtot   = TVV1 + TVV2,
    kel     = TVCL / TVV1,
    half_life_apparent = log(2) / kel
  ) |>
  transmute(
    Stratum            = stratum,
    n_paper            = c(27, 7, 11, 5),
    CRCL_mL_min        = CRCL_median,
    `CL (L/h)`         = round(TVCL, 2),
    `V1 (L)`           = round(TVV1, 2),
    `Q (L/h)`          = round(TVQ, 2),
    `V2 (L)`           = round(TVV2, 2),
    `V (L)`            = round(TVtot, 2),
    `kel (1/h)`        = round(kel, 3),
    `t1/2 apparent (h)` = round(half_life_apparent, 2)
  )

knitr::kable(
  stratum_typicals,
  caption = paste(
    "Conil 2007 Table 3 reproduction. n_paper is the stratum size in the",
    "published cohort; the simulated values come from the packaged",
    "model's final-model equations evaluated at the per-stratum",
    "median CRCL."
  ),
  align = c("l", "r", "r", "r", "r", "r", "r", "r", "r", "r")
)
```

| Stratum | n_paper | CRCL_mL_min | CL (L/h) | V1 (L) | Q (L/h) | V2 (L) | V (L) | kel (1/h) | t1/2 apparent (h) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Males, not ventilated | 27 | 114 | 7.19 | 18.81 | 6.88 | 3.96 | 22.77 | 0.382 | 1.81 |
| Females, not ventilated | 7 | 103 | 6.60 | 18.81 | 6.88 | 9.32 | 28.13 | 0.351 | 1.98 |
| Males, ventilated | 11 | 93 | 6.06 | 18.81 | 6.88 | 12.82 | 31.63 | 0.322 | 2.15 |
| Females, ventilated | 5 | 87 | 5.74 | 18.81 | 6.88 | 30.59 | 49.40 | 0.305 | 2.27 |

Conil 2007 Table 3 reproduction. n_paper is the stratum size in the
published cohort; the simulated values come from the packaged model’s
final-model equations evaluated at the per-stratum median CRCL. {.table}

Compared against the published Table 3 the model reproduces every value
to within rounding tolerance (males not vent: CL 7.19 vs 7.2 published,
V2 3.96 vs 4.0; females not vent: CL 6.60 vs 6.6, V2 9.32 vs 9.3; males
vent: CL 6.07 vs 6.1, V2 12.82 vs 12.8; females vent: CL 5.74 vs 5.7, V2
30.59 vs 30.6).

## Typical-value concentration-time profiles by stratum

``` r

sim_typical |>
  filter(time > 0) |>
  group_by(stratum, time) |>
  summarise(Cc = median(Cc, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(time, Cc, colour = stratum)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  labs(
    x = "Time (h)",
    y = "Ceftazidime in serum (mg/L, log scale)",
    title    = "Typical-value concentration-time profiles by stratum",
    subtitle = paste(
      "Single 2 g IV infusion over 20 min;",
      "per-stratum median CRCL from Conil 2007 Table 3"
    ),
    colour   = "Stratum"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](Conil_2007_ceftazidime_files/figure-html/typical-curves-1.png)

## Stochastic VPC by stratum

``` r

sim_stoch |>
  filter(time > 0) |>
  group_by(stratum, time) |>
  summarise(
    Q05  = quantile(Cc, 0.05, na.rm = TRUE),
    Q50  = quantile(Cc, 0.50, na.rm = TRUE),
    Q95  = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), fill = "gray70", alpha = 0.5) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ stratum, nrow = 2) +
  scale_y_log10() +
  labs(
    x = "Time (h)",
    y = "Ceftazidime in serum (mg/L, log scale)",
    title    = "Stochastic VPC by stratum",
    subtitle = "Median (line) and 5th-95th percentile band (ribbon)"
  ) +
  theme_minimal()
```

![](Conil_2007_ceftazidime_files/figure-html/vpc-1.png)

## PKNCA validation

PKNCA computes Cmax, Tmax, AUC0-Inf, and the terminal half-life on the
stochastic simulation. The Conil 2007 publication does not tabulate NCA
values (its Table 3 is structural-model-based with an apparent half-life
`ln(2) / (CL / V1)` that excludes distributional decay), so the PKNCA
values are presented here against Table 3 as a coarse cross-reference
and against the dose-divided-by-CL expectation as an exact sanity check.
The NCA terminal half-life is the true terminal half-life of the
two-compartment model and is expected to be longer than Table 3’s
apparent half-life, particularly for strata with a large peripheral
volume V2.

``` r

sim_nca <- sim_stoch |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc, stratum)

# Ensure a time=0 row per (id, stratum); for IV bolus / infusion the
# pre-dose concentration is 0.
sim_nca <- bind_rows(
  sim_nca,
  sim_nca |> distinct(id, stratum) |> mutate(time = 0, Cc = 0)
) |>
  distinct(id, stratum, time, .keep_all = TRUE) |>
  arrange(id, stratum, time)

dose_df <- events |>
  filter(evid == 1L) |>
  select(id, time, amt, stratum)

conc_obj <- PKNCA::PKNCAconc(
  data    = as.data.frame(sim_nca),
  formula = Cc ~ time | stratum + id,
  concu   = "mg/L",
  timeu   = "h"
)
dose_obj <- PKNCA::PKNCAdose(
  data    = as.data.frame(dose_df),
  formula = amt ~ time | stratum + id,
  doseu   = "mg"
)

intervals <- data.frame(
  start       = 0,
  end         = Inf,
  cmax        = TRUE,
  tmax        = TRUE,
  aucinf.obs  = TRUE,
  half.life   = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- suppressWarnings(PKNCA::pk.nca(nca_data))
```

### Comparison against Table 3

The reference values are derived from Conil 2007 Table 3: the apparent
half-life is `ln(2) / (CL / V1)` (Table 3 row “Half-life (h)”), and
`AUCinf` is `Dose / CL` at the per-stratum typical CL. Cmax is not
tabulated in the paper – Table 3 reports structural parameters but no
peak concentration – so the simulated Cmax appears here without a
reference value for comparison.

``` r

published <- tibble::tribble(
  ~stratum,                 ~aucinf.obs, ~half.life,
  "Males, not ventilated",  2000 / 7.2,  1.81,
  "Females, not ventilated",2000 / 6.6,  1.98,
  "Males, ventilated",      2000 / 6.1,  2.15,
  "Females, ventilated",    2000 / 5.7,  2.27
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_res,
  reference     = published,
  by            = "stratum",
  params        = c("aucinf.obs", "half.life"),
  units         = c(aucinf.obs = "mg*h/L", half.life = "h"),
  tolerance_pct = 30
)

knitr::kable(
  cmp,
  caption = paste(
    "Simulated NCA vs Conil 2007 Table 3 reference values.",
    "AUC0-Inf reference = Dose / CL with the per-stratum typical CL;",
    "t1/2 reference = Table 3 apparent half-life (= ln(2) * V1 / CL,",
    "the central-compartment elimination half-life, NOT the true",
    "terminal half-life). The two-compartment model's true terminal",
    "half-life is longer than the apparent half-life, particularly for",
    "strata with a large peripheral V2, so the NCA half-life is",
    "expected to exceed Table 3 by a stratum-dependent factor."
  ),
  align = c("l", "l", "r", "r", "r")
)
```

| NCA parameter          | stratum                 | Reference | Simulated |    % diff |
|:-----------------------|:------------------------|----------:|----------:|----------:|
| AUC0-∞ (obs) (mg\*h/L) | Males, not ventilated   |       278 |       279 |     +0.6% |
| AUC0-∞ (obs) (mg\*h/L) | Females, not ventilated |       303 |       299 |     -1.4% |
| AUC0-∞ (obs) (mg\*h/L) | Males, ventilated       |       328 |       333 |     +1.4% |
| AUC0-∞ (obs) (mg\*h/L) | Females, ventilated     |       351 |       342 |     -2.5% |
| t½ (h)                 | Males, not ventilated   |      1.81 |       2.5 |  +38.0%\* |
| t½ (h)                 | Females, not ventilated |      1.98 |       3.5 |  +76.6%\* |
| t½ (h)                 | Males, ventilated       |      2.15 |      4.29 |  +99.7%\* |
| t½ (h)                 | Females, ventilated     |      2.27 |      6.99 | +208.0%\* |

Simulated NCA vs Conil 2007 Table 3 reference values. AUC0-Inf reference
= Dose / CL with the per-stratum typical CL; t1/2 reference = Table 3
apparent half-life (= ln(2) \* V1 / CL, the central-compartment
elimination half-life, NOT the true terminal half-life). The
two-compartment model’s true terminal half-life is longer than the
apparent half-life, particularly for strata with a large peripheral V2,
so the NCA half-life is expected to exceed Table 3 by a
stratum-dependent factor. {.table style="width:100%;"}

``` r

fn <- attr(cmp, "footnote")
if (!is.null(fn)) cat(fn)
#> * differs from reference by more than ±30%.
```

### Simulated Cmax by stratum

``` r

cmax_tbl <- as.data.frame(nca_res$result) |>
  filter(PPTESTCD == "cmax") |>
  group_by(stratum) |>
  summarise(
    `Median Cmax (mg/L)` = round(median(PPORRES, na.rm = TRUE), 1),
    `Q05  Cmax (mg/L)`   = round(quantile(PPORRES, 0.05, na.rm = TRUE), 1),
    `Q95  Cmax (mg/L)`   = round(quantile(PPORRES, 0.95, na.rm = TRUE), 1),
    .groups = "drop"
  )

knitr::kable(
  cmax_tbl,
  caption = paste(
    "Simulated Cmax after a single 2 g 20-minute IV infusion.",
    "Conil 2007 does not tabulate Cmax in Table 3; the values are",
    "presented here as a sanity check that simulated peak",
    "concentrations are in the expected clinical range",
    "(~80-130 mg/L for a 2 g IV dose into V1 ~ 18.8 L)."
  ),
  align = c("l", "r", "r", "r")
)
```

| stratum                 | Median Cmax (mg/L) | Q05 Cmax (mg/L) | Q95 Cmax (mg/L) |
|:------------------------|-------------------:|----------------:|----------------:|
| Females, not ventilated |               94.9 |            67.6 |           116.4 |
| Females, ventilated     |               91.4 |            60.2 |           113.4 |
| Males, not ventilated   |               92.7 |            77.1 |           111.7 |
| Males, ventilated       |               93.8 |            65.4 |           114.5 |

Simulated Cmax after a single 2 g 20-minute IV infusion. Conil 2007 does
not tabulate Cmax in Table 3; the values are presented here as a sanity
check that simulated peak concentrations are in the expected clinical
range (~80-130 mg/L for a 2 g IV dose into V1 ~ 18.8 L). {.table
style="width:100%;"}

## Assumptions and deviations

- **IIV on Q and V2 and the residual error are carried forward from the
  basic model.** Conil 2007 Results p. 31 explicitly reports only the
  final-model IIV on CL (16%) and V1 (13%); the section is silent on IIV
  for Q and V2 in the final model and on the final-model residual error.
  Per the operator’s standing instruction (sidecar response,
  2026-06-17), silence is interpreted as “unchanged from the basic
  model” rather than “dropped to zero.” The packaged model therefore
  carries the basic-model values forward: Q IIV 363% (Table 2 row
  “Inter-compartmental clearance”), V2 IIV 190% (Table 2 row
  “Distribution volume of the peripheral compartment”), and proportional
  residual CV 38% (Results p. 30, “The residual variability… was 38%”).
  The 363% Q IIV is unusually large but matches the reported basic-model
  value verbatim; downstream users who simulate VPCs will see a wide
  spread on Q-driven distribution.

- **MECH_VENT canonical newly ratified.** Conil 2007 uses a binary
  indicator `VENT` (0 = not ventilated, 1 = mechanically ventilated) on
  the peripheral volume V2. No existing entry in
  `inst/references/covariate-columns.md` covered the concept (the
  `OXYSUP_HIGH` entry pools mechanical ventilation with high-flow oxygen
  in the Lin 2024 casirivimab analysis but is structurally a different
  stratification). This PR ratifies the new canonical `MECH_VENT` – a
  binary mechanical-ventilation treatment-status indicator with general
  scope, following the precedent of `HEMODIAL` and `CRRT_STATUS`. The
  orientation matches the source paper directly (0 = no MV, 1 = MV).

- **CL covariate equation is additive linear on the linear scale.**
  Conil 2007 final-model equation `CL = 1.08 + 0.0536 * CLCR` is an
  additive intercept-plus-slope model with no centering / scaling /
  normalization of CRCL. The packaged model encodes this as
  `cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)` with
  `exp(lcl) = 1.08` L/h and `e_crcl_cl = 0.0536` L/h per mL/min. At zero
  CRCL the model gives CL = 1.08 L/h (non-renal baseline), which is
  biologically meaningful (a small fraction of ceftazidime is cleared by
  non-renal routes); at the cohort-wide range (33-191 mL/min) CL spans
  2.85-11.32 L/h, consistent with the paper’s reported four-stratum
  range of 5.7-7.2 L/h at the per-stratum median CRCL.

- **V2 covariates are chained multiplicative fractional changes, not
  multiplicative percent or exponential terms.** Conil 2007 final-model
  equation
  `V2 = 2.69 * (1 + 1.43*SEX) * (1 + 2.44*VENT) * (1 + 0.00414*CLCR)` is
  a chained product of (1 + e \* COV) factors. This is the form the
  paper uses literally (Results p. 31, equation block); the packaged
  model preserves the chain rather than collapsing it to a single
  exponential or power form. The four-strata reproduction in “Replicate
  Table 3” confirms the form is correct.

- **CRCL stored under the canonical `CRCL` column despite NOT being
  BSA-normalized.** The canonical `CRCL` column accepts raw
  Cockcroft-Gault mL/min when the source paper does not BSA-normalize
  (precedent: `Delattre_2010_amikacin.R` for raw Cockcroft-Gault in a
  Belgian ICU cohort). Conil 2007 Methods records the Cockcroft-Gault
  equation as their CLCR derivation, with the standard 0.85 correction
  factor applied to women; no BSA normalization is mentioned. The
  per-model `covariateData[[CRCL]]$units` and `notes` document this.
  Reference values listed in the canonical entry (80, 90, 100
  mL/min/1.73 m^2) are for BSA-normalized estimates and should not be
  used as a comparison anchor for Conil 2007’s raw mL/min values.

- **Independent (diagonal) IIV between CL, V1, V2, and Q.** Conil 2007
  reports the IIV model as `theta_i = theta_pop * exp(eta_i)` (Results
  p. 30, basic-model section) with eta drawn from a normal distribution.
  The paper does not state whether OMEGA was diagonal or
  block-structured and reports a single CV per parameter with no
  off-diagonal covariance estimates, consistent with diagonal OMEGA. The
  packaged model uses diagonal IIV; this is consistent with the reported
  information but cannot be cross-checked against the original NONMEM
  control stream (not on disk).

- **`omega^2 = log(CV^2 + 1)`.** Conil 2007 reports interindividual
  variability as CV%; the corresponding log-normal variance was computed
  via `omega^2 = log(CV^2 + 1)` – the standard NONMEM/PsN
  back-transformation – and entered as the `eta...` initial value.

- **Race / ethnicity not modeled.** Conil 2007 does not report race
  composition; the cohort is a French university-hospital burn-ICU
  population. The analysis did not test race as a covariate and so no
  race effect is included.

- **NCA terminal half-life differs from Table 3 apparent half-life.**
  Conil 2007 Table 3 reports `t1/2 = ln(2) * V1 / CL` (computed from
  `kel = CL / V1`), which is the central-compartment elimination half-
  life – NOT the true terminal half-life of the two-compartment model.
  The simulated NCA half-life captures the true terminal slope and is
  longer than Table 3’s apparent half-life by a stratum- dependent
  factor (largest in the ventilated-females stratum, where V2 is largest
  relative to V1). The 30% tolerance flag used in the comparison table
  accommodates this methodological gap rather than flagging a model
  defect.

- **Dose units.** The model uses `mg` for dose and `mg/L` for
  concentration (paper convention for ceftazidime). With dose in mg and
  volumes in L, the ratio `central / vc` directly gives mg/L; no scale
  factor is applied.

- **Single-dose simulation in this vignette.** Conil 2007’s data
  collection captured trough and peak samples 24 h and 48 h after
  treatment start across q4h or q8h dosing regimens. The vignette
  simulates a single 2 g dose for NCA simplicity; multi-dose simulation
  (3 doses q8h or 6 doses q4h) is mechanically straightforward (add
  additional dose rows in the event table).
