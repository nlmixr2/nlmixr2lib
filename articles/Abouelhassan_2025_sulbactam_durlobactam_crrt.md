# Sulbactam-durlobactam during CRRT (Abouelhassan 2025)

## Model and source

- Citation: Abouelhassan Y, Shen Y, Chen A, Ye X, Nicolau DP, Kuti JL.
  Ex vivo assessment of sulbactam-durlobactam clearance during
  continuous renal replacement therapy to guide dosing recommendations.
  Antimicrob Agents Chemother. 2025;69(1):e01674-23.
  <doi:10.1128/aac.01674-23>. CRRT clearance regression from the Results
  (‘CRRT clearance studies’); PK parameters and dispersion from
  supplement Tables S11 and S12.
- Description (sulbactam): Two-compartment population PK simulation
  model for sulbactam (given as sulbactam-durlobactam) in critically ill
  adults with ventilator-associated bacterial pneumonia, with and
  without continuous renal replacement therapy (CRRT) (Abouelhassan
  2025). Without CRRT, total clearance is a non-renal arm plus a
  creatinine-clearance-driven renal arm, scaled by body weight and an
  infection-type shift (supplement Table S11). On CRRT the renal arm is
  replaced by the transmembrane clearance measured in an ex vivo CRRT
  circuit, a linear function of the effluent flow rate (CLTM = 0.114 +
  0.973 x effluent, L/h), added to the non-renal clearance (supplement
  Table S12). Unbound concentration uses a fixed 38% protein binding.
- Description (durlobactam): Two-compartment population PK simulation
  model for durlobactam (given as sulbactam-durlobactam) in critically
  ill adults with hospital-acquired or ventilator-associated bacterial
  pneumonia, with and without continuous renal replacement therapy
  (CRRT) (Abouelhassan 2025). Without CRRT, total clearance is a
  non-renal arm plus a creatinine-clearance-driven renal arm, scaled by
  body weight (supplement Table S11). On CRRT the renal arm is replaced
  by the transmembrane clearance measured in an ex vivo CRRT circuit, a
  linear function of the effluent flow rate (CLTM = 0.134 + 0.760 x
  effluent, L/h), added to the non-renal clearance (supplement Table
  S12). Unbound concentration uses a fixed 10% protein binding.
- Article: <https://doi.org/10.1128/aac.01674-23>

Abouelhassan et al. measured the transmembrane clearance (CL_(TM)) of
sulbactam and durlobactam in an ex vivo continuous renal replacement
therapy (CRRT) circuit and regressed it on the effluent flow rate. They
then replaced the renal clearance of an established two-compartment
population PK model with CL_(TM) and ran Monte Carlo simulations to
choose CRRT doses. The paper’s two drugs were modelled and simulated
independently, so the library carries one model per drug
(`Abouelhassan_2025_sulbactam`, `Abouelhassan_2025_durlobactam`). Each
model holds both clearance models of the paper and switches between them
on `RRT_CRRT_STATUS`:

- `RRT_CRRT_STATUS = 0` – the reference population of
  ventilator-associated bacterial pneumonia (VABP) patients not on CRRT
  (supplement Table S11).
- `RRT_CRRT_STATUS = 1` – CRRT patients, with clearance equal to the
  non-renal clearance plus CL_(TM) at the prescribed effluent flow
  (supplement Table S12).

## Population

No patient was studied in this paper. The ex vivo experiments circulated
one litre of heparinized bovine blood (haematocrit 38.7%, albumin 3.6
g/dL) through a Prismaflex circuit at a blood flow of 200 mL/min, in
CVVH and CVVHD modes, with M100 and HF1400 hemofilters and effluent
flows of 1, 2 and 3 L/h, each in triplicate (Methods, ‘Ex vivo CRRT
model’). The base population PK model was fitted to 162 patients from
the Phase 2 and Phase 3 sulbactam-durlobactam trials; its demographics
are not reported in this paper. The simulated populations were VABP
patients with a mean creatinine clearance of 80 mL/min (range 30-129)
and weight 76.4 +/- 4.7 kg, and CRRT patients of the same weight with
effluent flows of 1-5 L/h (supplement Tables S11 and S12).

The same information is available programmatically via
`readModelDb("Abouelhassan_2025_sulbactam")()$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in
`inst/modeldb/specificDrugs/Abouelhassan_2025_sulbactam.R` and
`inst/modeldb/specificDrugs/Abouelhassan_2025_durlobactam.R`.

| Parameter | Sulbactam | Durlobactam | Source location |
|----|----|----|----|
| `lcl_nonren` (non-renal CL, L/h) | 4.79 | 3.18 | Table S11 ‘Total Body CL’; Table S12 ‘CL_NR’; supplement Methods |
| `lcl` (population PK CL, L/h) | 13.6 | 9.36 | Table S11 ‘Total Body CL’ |
| `f_renal` (renal fraction) | 0.648 | 0.66 | Table S11 ‘Total Body CL’ |
| `e_crcl_cl_renal` (exponent on CLcr/90) | 0.932 | 0.694 | Table S11 ‘Total Body CL’ |
| `e_wt_cl` (exponent on WT/75, non-CRRT CL) | 1.03 | 0.664 | Table S11 ‘Total Body CL’ |
| `e_vabp_cl` (VABP shift on CL, non-CRRT) | -0.294 | – | Table S11 footnote a |
| `lvc` (Vc at 75 kg, L) | 12.1 | 12.5 | Table S11 ‘Vc’ |
| `e_wt_vc` (exponent on WT/75) | 0.828 | 0.518 | Table S11 ‘Vc’ |
| `e_vabp_vc` / `e_habp_vabp_vc` (infection shift on Vc) | 1.38 | 1.52 | Table S11 footnote b; Table S12 ‘Vc’ |
| `lq` (Q, L/h) | 7.82 | 4.42 | Tables S11 and S12 ‘Q’ |
| `lvp` (Vp, L) | 6.98 | 5.83 | Tables S11 and S12 ‘Vp’ |
| `lcl_crrt` (CL_(TM) intercept, L/h) | 0.114 | 0.134 | Results ‘CRRT clearance studies’; Fig. S1 |
| `e_effl_cl_crrt` (CL_(TM) slope on effluent, L/h per L/h) | 0.973 | 0.760 | Results ‘CRRT clearance studies’; Table 2; Fig. S1 |
| `etalvp` (variance) | 0.647 (CV 95.4%) | 0.597 (CV 90.4%) | Tables S11/S12 ‘Vp’ mean +/- SD |
| `etalcl_nonren` (variance) | 0.0862 (CV 30%) | 0.0862 (CV 30%) | supplement Methods; Table S12 ‘CL_NR’ |
| `etalcl_crrt` (variance) | 0.0143 (CV 12%) | 0.0143 (CV 12%) | supplement Methods |
| `fu` (unbound fraction) | 0.62 | 0.90 | supplement Methods (protein binding 38% / 10%) |
| `propSd` | 0 (fixed) | 0 (fixed) | not reported (simulation study) |
| Non-CRRT CL equation | – | – | Table S11 ‘Total Body CL’ |
| CRRT CL = CL_NR + CL_TM | – | – | Table S12 ‘Total Body CL’; supplement Methods |
| Two-compartment ODEs | – | – | supplement Methods (‘two compartment population pharmacokinetic analysis’) |

## Transmembrane clearance against effluent flow

The CRRT clearance is a deterministic function of effluent flow. Solving
the typical-value models at the effluent flows of supplement Tables S3
and S4 reproduces the ‘ER only’ column of those tables.

``` r

mod_sul <- readModelDb("Abouelhassan_2025_sulbactam")
mod_dur <- readModelDb("Abouelhassan_2025_durlobactam")

effluent_grid <- c(1, 2, 2.5, 3, 4, 5)

typical_cltm <- function(mod, drug) {
  ev <- data.frame(
    id = seq_along(effluent_grid),
    time = 0,
    evid = 0,
    amt = 0,
    cmt = "central",
    WT = 75,
    CRCL = 90,
    DIS_HABP = 0,
    DIS_VABP = 1,
    RRT_CRRT_STATUS = 1,
    RRT_CRRT_EFFLUENT_FLOW = effluent_grid * 1000
  )
  sim <- rxode2::rxSolve(rxode2::zeroRe(mod), ev, returnType = "data.frame")
  data.frame(drug = drug, effluent = effluent_grid, cl_crrt = sim$cl_crrt)
}

cltm <- bind_rows(
  typical_cltm(mod_sul, "Sulbactam"),
  typical_cltm(mod_dur, "Durlobactam")
) |>
  mutate(
    published = c(
      1.09, 2.06, 2.55, 3.03, 4.01, 4.98, # Table S3 'ER only'
      0.89, 1.65, 2.03, 2.41, 3.17, 3.93 # Table S4 'ER only'
    )
  )
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalvp', 'etalcl_nonren', 'etalcl_crrt'
#> Warning: multi-subject simulation without without 'omega'
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalvp', 'etalcl_nonren', 'etalcl_crrt'
#> Warning: multi-subject simulation without without 'omega'

# Same parameters on both sides: the only difference is the table's rounding.
stopifnot(all(abs(cltm$cl_crrt - cltm$published) < 0.006))

cltm |>
  mutate(cl_crrt = signif(cl_crrt, 4)) |>
  rename(
    "Drug" = drug,
    "Effluent flow (L/h)" = effluent,
    "Model CL_TM (L/h)" = cl_crrt,
    "Published CL_TM (L/h)" = published
  ) |>
  knitr::kable(caption = "Typical CL_TM vs supplement Tables S3 and S4 ('ER only').")
```

| Drug        | Effluent flow (L/h) | Model CL_TM (L/h) | Published CL_TM (L/h) |
|:------------|--------------------:|------------------:|----------------------:|
| Sulbactam   |                 1.0 |             1.087 |                  1.09 |
| Sulbactam   |                 2.0 |             2.060 |                  2.06 |
| Sulbactam   |                 2.5 |             2.547 |                  2.55 |
| Sulbactam   |                 3.0 |             3.033 |                  3.03 |
| Sulbactam   |                 4.0 |             4.006 |                  4.01 |
| Sulbactam   |                 5.0 |             4.979 |                  4.98 |
| Durlobactam |                 1.0 |             0.894 |                  0.89 |
| Durlobactam |                 2.0 |             1.654 |                  1.65 |
| Durlobactam |                 2.5 |             2.034 |                  2.03 |
| Durlobactam |                 3.0 |             2.414 |                  2.41 |
| Durlobactam |                 4.0 |             3.174 |                  3.17 |
| Durlobactam |                 5.0 |             3.934 |                  3.93 |

Typical CL_TM vs supplement Tables S3 and S4 (‘ER only’). {.table}

``` r

ggplot(cltm, aes(effluent, cl_crrt, colour = drug)) +
  geom_line() +
  geom_point(aes(y = published), shape = 1, size = 2.5) +
  labs(x = "Effluent flow rate (L/h)", y = "Transmembrane clearance (L/h)", colour = NULL) +
  theme_bw()
```

![Replicates Fig. S1 of Abouelhassan 2025: CL_TM against effluent flow
rate (lines: model; points: Tables S3 and
S4).](Abouelhassan_2025_sulbactam_durlobactam_crrt_files/figure-html/cltm-plot-1.png)

Replicates Fig. S1 of Abouelhassan 2025: CL_TM against effluent flow
rate (lines: model; points: Tables S3 and S4).

## Virtual cohort

The paper simulated 1,000 patients per scenario; the cohorts below use
100 per scenario. Weight is drawn from N(76.4, 4.7) kg for every
scenario, and the non-CRRT VABP creatinine clearance from a uniform
distribution over 30-129 mL/min (the paper reports a mean of 80 and the
range only). Each scenario is a 3-h infusion of 1 g of each drug at the
stated interval.

``` r

rxode2::rxSetSeed(20250101)
set.seed(20250101)

n_per <- 100L

scenarios <- tibble::tribble(
  ~scenario,          ~crrt, ~effluent, ~tau,
  "No CRRT, 1g q6h",  0,     0,         6,
  "1 L/h, 1g q8h",    1,     1,         8,
  "1 L/h, 1g q12h",   1,     1,         12,
  "2 L/h, 1g q6h",    1,     2,         6,
  "2 L/h, 1g q8h",    1,     2,         8,
  "2.5 L/h, 1g q6h",  1,     2.5,       6,
  "2.5 L/h, 1g q8h",  1,     2.5,       8,
  "3 L/h, 1g q6h",    1,     3,         6,
  "3 L/h, 1g q8h",    1,     3,         8,
  "4 L/h, 1g q6h",    1,     4,         6,
  "4 L/h, 1g q8h",    1,     4,         8,
  "5 L/h, 1g q6h",    1,     5,         6,
  "5 L/h, 1g q4h",    1,     5,         4
) |>
  mutate(
    scenario = factor(scenario, levels = scenario),
    # Table S5/S6 footnote b: q12h was evaluated over 48-72 h, all others 24-48 h
    win_start = ifelse(tau == 12, 48, 24),
    win_end = win_start + 24
  )

subjects <- scenarios |>
  slice(rep(seq_len(n()), each = n_per)) |>
  mutate(
    id = row_number(),
    WT = rnorm(n(), 76.4, 4.7),
    CRCL = ifelse(crrt == 0, runif(n(), 30, 129), 90),
    DIS_HABP = 0,
    DIS_VABP = 1,
    RRT_CRRT_STATUS = crrt,
    RRT_CRRT_EFFLUENT_FLOW = effluent * 1000
  )

dose_rows <- subjects |>
  group_by(id) |>
  reframe(
    time = seq(0, win_end - tau, by = tau),
    across(c(scenario, WT, CRCL, DIS_HABP, DIS_VABP, RRT_CRRT_STATUS, RRT_CRRT_EFFLUENT_FLOW))
  ) |>
  mutate(evid = 1, amt = 1000, rate = 1000 / 3, cmt = "central")

obs_rows <- subjects |>
  group_by(id) |>
  reframe(
    time = sort(unique(c(0, seq(win_start, win_end, by = 0.1)))),
    across(c(scenario, WT, CRCL, DIS_HABP, DIS_VABP, RRT_CRRT_STATUS, RRT_CRRT_EFFLUENT_FLOW))
  ) |>
  mutate(evid = 0, amt = 0, rate = 0, cmt = "central")

events <- bind_rows(dose_rows, obs_rows) |>
  arrange(id, time, desc(evid)) |>
  as.data.frame()

stopifnot(!anyDuplicated(events[events$evid == 1, c("id", "time")]))
```

## Simulation

``` r

simulate_drug <- function(mod, drug) {
  rxode2::rxSolve(mod, events, keep = "scenario", returnType = "data.frame") |>
    mutate(drug = drug)
}
sim <- bind_rows(
  simulate_drug(mod_sul, "Sulbactam"),
  simulate_drug(mod_dur, "Durlobactam")
) |>
  left_join(select(scenarios, scenario, tau, win_start, win_end), by = "scenario")
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ parameter labels from comments will be replaced by 'label()'
```

``` r

sim |>
  filter(
    scenario %in% c("No CRRT, 1g q6h", "2 L/h, 1g q8h", "3 L/h, 1g q6h"),
    time >= 24, time <= 48
  ) |>
  group_by(drug, scenario, time) |>
  summarise(
    p05 = quantile(Cc, 0.05), p50 = median(Cc), p95 = quantile(Cc, 0.95),
    .groups = "drop"
  ) |>
  ggplot(aes(time, p50, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.15, colour = NA) +
  geom_line() +
  facet_wrap(~drug) +
  labs(x = "Time (h)", y = "Total plasma concentration (mg/L)", colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
```

![Median and 5th-95th percentile total plasma concentrations over the
24-48 h window for the regimens the paper recommends (1g q8h for
effluent \< 3 L/h, 1g q6h for 3-5 L/h) and the non-CRRT
reference.](Abouelhassan_2025_sulbactam_durlobactam_crrt_files/figure-html/profiles-1.png)

Median and 5th-95th percentile total plasma concentrations over the
24-48 h window for the regimens the paper recommends (1g q8h for
effluent \< 3 L/h, 1g q6h for 3-5 L/h) and the non-CRRT reference.

## PKNCA: total-drug AUC over the evaluation window

``` r

# A 0.25-h subsample of the simulation grid keeps PKNCA fast; the linear
# trapezoidal AUC moves by well under 1% against the full grid.
conc_df <- sim |>
  filter(!is.na(Cc), abs(time * 4 - round(time * 4)) < 1e-8) |>
  select(id, time, Cc, drug, scenario)
dose_df <- events |>
  filter(evid == 1) |>
  select(id, time, amt, scenario) |>
  tidyr::expand_grid(drug = c("Sulbactam", "Durlobactam"))

intervals <- scenarios |>
  distinct(scenario, win_start, win_end) |>
  transmute(scenario, start = win_start, end = win_end, aucint.last = TRUE)

conc_obj <- PKNCA::PKNCAconc(conc_df, Cc ~ time | drug + scenario / id)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | drug + scenario + id)
nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))

auc_ind <- as.data.frame(nca_res) |>
  filter(PPTESTCD == "aucint.last") |>
  select(drug, scenario, id, auc = PPORRES)
```

The paper reports the mean, SD and 25th/75th percentiles of the
total-drug AUC₂₄₋₄₈ (Tables S5 and S6), not the median that
[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
would compute from individual results, so the simulated results are
first aggregated to the per-scenario mean.

``` r

published <- tibble::tribble(
  ~scenario,          ~Sulbactam, ~Durlobactam,
  "No CRRT, 1g q6h",  466, 452,
  "1 L/h, 1g q8h",    533, 735,
  "1 L/h, 1g q12h",   362, 505,
  "2 L/h, 1g q6h",    607, 838,
  "2 L/h, 1g q8h",    456, 630,
  "2.5 L/h, 1g q6h",  563, 777,
  "2.5 L/h, 1g q8h",  420, 581,
  "3 L/h, 1g q6h",    526, 726,
  "3 L/h, 1g q8h",    399, 544,
  "4 L/h, 1g q6h",    469, 636,
  "4 L/h, 1g q8h",    348, 478,
  "5 L/h, 1g q6h",    420, 566,
  "5 L/h, 1g q4h",    622, 844
) |>
  pivot_longer(-scenario, names_to = "drug", values_to = "aucint.last") |>
  mutate(group = paste(drug, scenario, sep = ": ")) |>
  select(group, aucint.last)

sim_mean <- auc_ind |>
  group_by(drug, scenario) |>
  summarise(aucint.last = mean(auc), .groups = "drop") |>
  mutate(group = paste(drug, scenario, sep = ": ")) |>
  select(group, aucint.last)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = sim_mean,
  reference = published,
  by = "group",
  units = c(aucint.last = "mg*h/L"),
  tolerance_pct = 20
)
#> Warning: ncaParamLabel(): unknown PKNCA code(s) returned as-is: 'aucint.last'
knitr::kable(
  cmp,
  caption = "Mean total-drug AUC24-48 (AUC48-72 for q12h): simulated vs supplement Tables S5 and S6. * differs by more than 20%."
)
```

| NCA parameter | group | Reference | Simulated | % diff |
|:---|:---|:---|:---|:---|
| aucint.last (mg\*h/L) | Sulbactam: No CRRT, 1g q6h | 466 | 466 | -0.0% |
| aucint.last (mg\*h/L) | Durlobactam: No CRRT, 1g q6h | 452 | 467 | +3.3% |
| aucint.last (mg\*h/L) | Sulbactam: 1 L/h, 1g q8h | 533 | 532 | -0.2% |
| aucint.last (mg\*h/L) | Durlobactam: 1 L/h, 1g q8h | 735 | 705 | -4.1% |
| aucint.last (mg\*h/L) | Sulbactam: 1 L/h, 1g q12h | 362 | 357 | -1.5% |
| aucint.last (mg\*h/L) | Durlobactam: 1 L/h, 1g q12h | 505 | 494 | -2.1% |
| aucint.last (mg\*h/L) | Sulbactam: 2 L/h, 1g q6h | 607 | 577 | -5.0% |
| aucint.last (mg\*h/L) | Durlobactam: 2 L/h, 1g q6h | 838 | 833 | -0.6% |
| aucint.last (mg\*h/L) | Sulbactam: 2 L/h, 1g q8h | 456 | 431 | -5.4% |
| aucint.last (mg\*h/L) | Durlobactam: 2 L/h, 1g q8h | 630 | 612 | -2.9% |
| aucint.last (mg\*h/L) | Sulbactam: 2.5 L/h, 1g q6h | 563 | 526 | -6.6% |
| aucint.last (mg\*h/L) | Durlobactam: 2.5 L/h, 1g q6h | 777 | 757 | -2.6% |
| aucint.last (mg\*h/L) | Sulbactam: 2.5 L/h, 1g q8h | 420 | 411 | -2.1% |
| aucint.last (mg\*h/L) | Durlobactam: 2.5 L/h, 1g q8h | 581 | 578 | -0.6% |
| aucint.last (mg\*h/L) | Sulbactam: 3 L/h, 1g q6h | 526 | 532 | +1.1% |
| aucint.last (mg\*h/L) | Durlobactam: 3 L/h, 1g q6h | 726 | 708 | -2.5% |
| aucint.last (mg\*h/L) | Sulbactam: 3 L/h, 1g q8h | 399 | 388 | -2.7% |
| aucint.last (mg\*h/L) | Durlobactam: 3 L/h, 1g q8h | 544 | 535 | -1.6% |
| aucint.last (mg\*h/L) | Sulbactam: 4 L/h, 1g q6h | 469 | 457 | -2.6% |
| aucint.last (mg\*h/L) | Durlobactam: 4 L/h, 1g q6h | 636 | 623 | -2.0% |
| aucint.last (mg\*h/L) | Sulbactam: 4 L/h, 1g q8h | 348 | 343 | -1.3% |
| aucint.last (mg\*h/L) | Durlobactam: 4 L/h, 1g q8h | 478 | 470 | -1.7% |
| aucint.last (mg\*h/L) | Sulbactam: 5 L/h, 1g q6h | 420 | 410 | -2.3% |
| aucint.last (mg\*h/L) | Durlobactam: 5 L/h, 1g q6h | 566 | 557 | -1.5% |
| aucint.last (mg\*h/L) | Sulbactam: 5 L/h, 1g q4h | 622 | 635 | +2.1% |
| aucint.last (mg\*h/L) | Durlobactam: 5 L/h, 1g q4h | 844 | 838 | -0.7% |

Mean total-drug AUC24-48 (AUC48-72 for q12h): simulated vs supplement
Tables S5 and S6. \* differs by more than 20%. {.table}

``` r


pct <- merge(sim_mean, published, by = "group", suffixes = c("_sim", "_pub")) |>
  mutate(pct_diff = 100 * (aucint.last_sim / aucint.last_pub - 1))
stopifnot(
  nrow(pct) == 26,
  # A mis-transcribed clearance, slope or dose shifts every scenario of a drug
  # by tens of percent.
  abs(median(pct$pct_diff)) < 10,
  quantile(abs(pct$pct_diff), 0.9) < 20
)
```

## Probability of target attainment

The targets are 50% fT\>MIC for sulbactam and fAUC/MIC \>= 10 for
durlobactam, evaluated on the unbound concentration `Cu` over the same
24-h window.

``` r

mics <- c(2, 4, 8, 16, 32)

sim_win <- sim |>
  filter(time >= win_start, time <= win_end) |>
  arrange(drug, id, time)

fauc_ind <- sim_win |>
  group_by(drug, scenario, id) |>
  summarise(
    fauc = sum(diff(time) * (head(Cu, -1) + tail(Cu, -1)) / 2),
    .groups = "drop"
  )

# On the uniform grid, the share of grid intervals whose left end is above the
# MIC approximates fT>MIC.
ft_ind <- sim_win |>
  filter(time < win_end) |>
  select(drug, scenario, id, Cu) |>
  tidyr::expand_grid(mic = mics) |>
  group_by(drug, scenario, id, mic) |>
  summarise(ft_above = mean(Cu > mic), .groups = "drop")

pta_ind <- ft_ind |>
  left_join(fauc_ind, by = c("drug", "scenario", "id")) |>
  mutate(
    attained = ifelse(drug == "Sulbactam", ft_above >= 0.5, fauc / mic >= 10)
  )

pta <- pta_ind |>
  group_by(drug, scenario, mic) |>
  summarise(pta = mean(attained), .groups = "drop")

pta |>
  mutate(mic = paste0("MIC ", mic)) |>
  pivot_wider(names_from = mic, values_from = pta) |>
  rename("Drug" = drug, "Scenario" = scenario) |>
  knitr::kable(digits = 2, caption = "Simulated PTA (compare supplement Tables S5 and S6).")
```

| Drug        | Scenario        | MIC 2 | MIC 4 | MIC 8 | MIC 16 | MIC 32 |
|:------------|:----------------|------:|------:|------:|-------:|-------:|
| Durlobactam | No CRRT, 1g q6h |     1 |     1 |  1.00 |   1.00 |   0.98 |
| Durlobactam | 1 L/h, 1g q8h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 1 L/h, 1g q12h  |     1 |     1 |  1.00 |   1.00 |   0.92 |
| Durlobactam | 2 L/h, 1g q6h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 2 L/h, 1g q8h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 2.5 L/h, 1g q6h |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 2.5 L/h, 1g q8h |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 3 L/h, 1g q6h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 3 L/h, 1g q8h   |     1 |     1 |  1.00 |   1.00 |   0.98 |
| Durlobactam | 4 L/h, 1g q6h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 4 L/h, 1g q8h   |     1 |     1 |  1.00 |   1.00 |   0.98 |
| Durlobactam | 5 L/h, 1g q6h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Durlobactam | 5 L/h, 1g q4h   |     1 |     1 |  1.00 |   1.00 |   1.00 |
| Sulbactam   | No CRRT, 1g q6h |     1 |     1 |  0.99 |   0.10 |   0.00 |
| Sulbactam   | 1 L/h, 1g q8h   |     1 |     1 |  0.98 |   0.24 |   0.00 |
| Sulbactam   | 1 L/h, 1g q12h  |     1 |     1 |  0.58 |   0.01 |   0.00 |
| Sulbactam   | 2 L/h, 1g q6h   |     1 |     1 |  1.00 |   0.31 |   0.00 |
| Sulbactam   | 2 L/h, 1g q8h   |     1 |     1 |  0.91 |   0.01 |   0.00 |
| Sulbactam   | 2.5 L/h, 1g q6h |     1 |     1 |  1.00 |   0.22 |   0.00 |
| Sulbactam   | 2.5 L/h, 1g q8h |     1 |     1 |  0.85 |   0.01 |   0.00 |
| Sulbactam   | 3 L/h, 1g q6h   |     1 |     1 |  0.97 |   0.17 |   0.00 |
| Sulbactam   | 3 L/h, 1g q8h   |     1 |     1 |  0.79 |   0.00 |   0.00 |
| Sulbactam   | 4 L/h, 1g q6h   |     1 |     1 |  0.99 |   0.02 |   0.00 |
| Sulbactam   | 4 L/h, 1g q8h   |     1 |     1 |  0.63 |   0.00 |   0.00 |
| Sulbactam   | 5 L/h, 1g q6h   |     1 |     1 |  0.94 |   0.00 |   0.00 |
| Sulbactam   | 5 L/h, 1g q4h   |     1 |     1 |  1.00 |   0.62 |   0.00 |

Simulated PTA (compare supplement Tables S5 and S6). {.table}

``` r


# The paper's headline claim: every simulated regimen achieves PTA > 0.9 at the
# 4/4 mg/L susceptibility breakpoint.
stopifnot(all(pta$pta[pta$mic == 4] > 0.9))
```

Durlobactam attains its target in nearly every subject at every MIC, as
in Table S6. For sulbactam the discriminating MICs are 8 and 16 mg/L:

``` r

pta_sul_pub <- tibble::tribble(
  ~scenario,          ~mic8, ~mic16,
  "No CRRT, 1g q6h",  0.14,  0,
  "1 L/h, 1g q8h",    1,     0.21,
  "1 L/h, 1g q12h",   0.53,  0,
  "2 L/h, 1g q6h",    1,     0.41,
  "2 L/h, 1g q8h",    1,     0.04,
  "2.5 L/h, 1g q6h",  1,     0.271,
  "2.5 L/h, 1g q8h",  0.98,  0.01,
  "3 L/h, 1g q6h",    1,     0.16,
  "3 L/h, 1g q8h",    0.93,  0,
  "4 L/h, 1g q6h",    1,     0.02,
  "4 L/h, 1g q8h",    0.71,  0,
  "5 L/h, 1g q6h",    1,     0,
  "5 L/h, 1g q4h",    1,     0.54
)

pta |>
  filter(drug == "Sulbactam", mic %in% c(8, 16)) |>
  mutate(scenario = as.character(scenario), mic = paste0("sim", mic)) |>
  pivot_wider(names_from = mic, values_from = pta) |>
  select(-drug) |>
  left_join(pta_sul_pub, by = "scenario") |>
  select(scenario, mic8, sim8, mic16, sim16) |>
  rename(
    "Scenario" = scenario,
    "Published PTA, MIC 8" = mic8,
    "Simulated PTA, MIC 8" = sim8,
    "Published PTA, MIC 16" = mic16,
    "Simulated PTA, MIC 16" = sim16
  ) |>
  knitr::kable(digits = 2, caption = "Sulbactam PTA (50% fT>MIC) at MIC 8 and 16 mg/L: simulated vs supplement Table S5.")
```

| Scenario | Published PTA, MIC 8 | Simulated PTA, MIC 8 | Published PTA, MIC 16 | Simulated PTA, MIC 16 |
|:---|---:|---:|---:|---:|
| No CRRT, 1g q6h | 0.14 | 0.99 | 0.00 | 0.10 |
| 1 L/h, 1g q8h | 1.00 | 0.98 | 0.21 | 0.24 |
| 1 L/h, 1g q12h | 0.53 | 0.58 | 0.00 | 0.01 |
| 2 L/h, 1g q6h | 1.00 | 1.00 | 0.41 | 0.31 |
| 2 L/h, 1g q8h | 1.00 | 0.91 | 0.04 | 0.01 |
| 2.5 L/h, 1g q6h | 1.00 | 1.00 | 0.27 | 0.22 |
| 2.5 L/h, 1g q8h | 0.98 | 0.85 | 0.01 | 0.01 |
| 3 L/h, 1g q6h | 1.00 | 0.97 | 0.16 | 0.17 |
| 3 L/h, 1g q8h | 0.93 | 0.79 | 0.00 | 0.00 |
| 4 L/h, 1g q6h | 1.00 | 0.99 | 0.02 | 0.02 |
| 4 L/h, 1g q8h | 0.71 | 0.63 | 0.00 | 0.00 |
| 5 L/h, 1g q6h | 1.00 | 0.94 | 0.00 | 0.00 |
| 5 L/h, 1g q4h | 1.00 | 1.00 | 0.54 | 0.62 |

Sulbactam PTA (50% fT\>MIC) at MIC 8 and 16 mg/L: simulated vs
supplement Table S5. {.table}

The CRRT rows follow Table S5: high attainment at 8 mg/L that falls with
effluent flow under q8h dosing, and low attainment at 16 mg/L. The
simulated attainment at 8 mg/L under q8h is somewhat lower than
published at the higher effluent flows (for example about 0.8 vs 0.93 at
3 L/h), consistent with the untruncated log-normal dispersion used here
(see Assumptions and deviations).

The non-CRRT row is **not reproduced**: Table S5 gives a PTA of 0.14 at
8 mg/L, the simulation gives close to 1. The published row is hard to
reconcile with the same table: its mean AUC₂₄₋₄₈ of 466 mg*h/L implies
an average unbound concentration of 0.62 x 466 / 24 = 12 mg/L over the
window, and the 2 L/h q8h CRRT row, with a lower mean AUC (456 mg*h/L),
reports a PTA of 1 at the same MIC. The mean AUC of this row is
reproduced to within a few percent, so the clearance is right; the PTA
difference must come from the shape of the profile or the PTA
computation for this row, which the paper does not describe. The model
was not adjusted to match it.

## Assumptions and deviations

- **Dispersion as log-normal.** The authors simulated with the Pmetrics
  non-parametric simulator from the means and SDs of Tables S11 and S12,
  with parameter limits (for example Vp 1-30 L, CL_NR 2-7 L/h for
  sulbactam). The library carries the SDs as log-normal between-subject
  variability, `omega^2 = log(CV^2 + 1)`, and does not truncate. The
  very wide Vp dispersion (CV 95% and 90%) therefore yields some large
  Vp values the authors’ limits would have clipped; this affects the
  shape of the profile, not the AUC.
- **No variability on Vc and Q.** Tables S11 and S12 print ranges for Vc
  and Q (for example Vc 20-45 L for sulbactam) but no SD, so no
  between-subject variability is carried on them; Vc varies only through
  body weight.
- **Variability on the non-renal clearance applies on CRRT only.** Table
  S11 gives the non-CRRT clearance no dispersion, so the 30% CV of the
  supplement Methods enters only the CRRT clearance; in the non-CRRT arm
  the spread in clearance comes from creatinine clearance and weight
  alone.
- **Creatinine clearance distribution.** The paper gives only the mean
  (80 mL/min) and range (30-129 mL/min) of the simulated VABP
  population; the cohort above draws it uniformly over the range, whose
  mean is 79.5 mL/min.
- **Infection type.** The base model’s infection-type effects are
  printed only for the simulated population: a VABP shift of -0.294 on
  sulbactam clearance and HABP/VABP shifts of 1.38 (sulbactam) and 1.52
  (durlobactam) on Vc. They are encoded on `DIS_VABP` (sulbactam) and
  `DIS_HABP + DIS_VABP` (durlobactam). The sulbactam clearance shift for
  HABP is not printed, so HABP is not encoded for sulbactam. Set
  `DIS_VABP = 1` to reproduce the paper; the reference category
  (`DIS_VABP = 0`) of the unpublished base model is not identified in
  the source.
- **Base population PK model.** The base model (162 Phase 2/3 patients)
  is not published separately in this paper and is cited as data on
  file. It is close to, but not the same as, the later published model
  of Cammarata 2024 (`Cammarata_2024_sulbactam_durlobactam`, 373
  subjects, creatinine clearance BSA-normalized with reference 100
  mL/min/1.73 m^2). The non-renal clearances of 4.79 and 3.18 L/h are
  numerically equal to (1 - f_renal) x CL of the base model.
- **Final CL_(TM) model only.** The multiple regressions on filter type,
  CRRT mode and replacement-fluid dilution (Table 2 and Table S2) were
  used by the authors only for sensitivity analyses (Tables S7 and S8);
  the final model that drives the dose recommendations uses effluent
  flow alone and is the one encoded.
- **No residual error.** The source is a simulation study with no
  residual error model; `propSd` is fixed to 0.
- **PTA definitions.** fT\>MIC is computed as the fraction of the 0.1-h
  grid intervals with an unbound concentration above the MIC, and fAUC
  by the linear trapezoidal rule, over the 24-h evaluation window. The
  paper does not say over which window PTA was assessed; the AUC window
  of Tables S5 and S6 is used.
- **Sulbactam PTA without CRRT not reproduced.** Supplement Table S5
  reports a PTA of 0.14 at 8 mg/L for the non-CRRT 1g q6h reference; the
  model gives close to 1 while reproducing that row’s mean AUC. See the
  discussion under “Probability of target attainment”.
- No erratum or correction notice was found for this article (Europe
  PMC, checked 2026-11-02).
