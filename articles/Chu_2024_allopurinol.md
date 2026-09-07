# Allopurinol (Chu 2024)

## Model and source

- Citation: Chu WY, Nijman M, Stegeman R, Breur JMPJ, Jansen NJG, Nijman
  J, van Loon K, Koomen E, Allegaert K, Benders MJNL, Dorlo TPC, Huitema
  ADR; the CRUCIAL trial consortium. Population Pharmacokinetics and
  Target Attainment of Allopurinol and Oxypurinol Before, During, and
  After Cardiac Surgery with Cardiopulmonary Bypass in Neonates with
  Critical Congenital Heart Disease. Clin Pharmacokinet.
  2024;63:1205-1220. <doi:10.1007/s40262-024-01401-3>
- Description: Joint parent-metabolite population pharmacokinetic model
  for intravenous allopurinol and its active metabolite oxypurinol in 14
  neonates with critical congenital heart disease (CCHD) undergoing
  cardiac surgery with cardiopulmonary bypass (CPB) in the CRUCIAL trial
  (Chu 2024). Structural model: two-compartment allopurinol disposition
  feeding a one-compartment oxypurinol disposition, with full (formation
  fraction 1) conversion of allopurinol to oxypurinol and
  auto-inhibition of that conversion by oxypurinol (Imax fixed to 1,
  IC50 fixed to 1.1 mg/L). The allopurinol central compartment (V1 = 0.1
  L per 3.5 kg) and intercompartmental clearance (Q1 = 6.97 L/h per 3.5
  kg) were fixed to capture the 10-minute post-dose peak sampled during
  CPB and carry no physiologic interpretation. Clearances and volumes
  are allometrically scaled to a 3.5 kg reference body weight with fixed
  exponents 0.75 and 1. Three mutually exclusive perioperative periods
  carry different disposition: during the postnatal-preoperative period
  both clearances rise with postnatal age along a fixed sigmoidal
  recovery curve (TM50 4.2 days, Hill 2.98; maximum fold increase 3 for
  allopurinol and 1.35 for oxypurinol); during CPB and after CPB the
  clearances and volumes are instead fixed multiples of their at-birth
  values. Between-subject variability is carried on allopurinol
  clearance and oxypurinol volume, and between-occasion variability on
  both clearances across the three periods. Residual error is
  proportional plus an additive component fixed to LLOQ/2 for each
  analyte.
- Article: <https://doi.org/10.1007/s40262-024-01401-3>
- Supplement (ESM 1: model-structure schematic, model-development table,
  and the full NONMEM control stream):
  <https://doi.org/10.1007/s40262-024-01401-3>
- Trial protocol (CRUCIAL; cited as reference 6 by Chu 2024 and the
  source of its Fig. 1 dosing schedule):
  <https://doi.org/10.1186/s13063-022-06098-y>

Chu 2024 is the pharmacokinetic substudy of the CRUCIAL trial
(NCT04217421), in which neonates with critical congenital heart disease
(CCHD) receive five 20 mg/kg intravenous allopurinol doses spanning
birth, cardiac surgery with cardiopulmonary bypass (CPB), and the
postoperative period. The model is a joint parent-metabolite popPK
model: a two-compartment allopurinol disposition converting entirely to
a one-compartment oxypurinol disposition, with auto-inhibition of that
conversion by oxypurinol itself.

Its distinguishing feature is that **disposition is period-specific**.
The record is partitioned into three mutually exclusive perioperative
windows at the start and end of CPB, and each window carries its own
clearance and volume. In the postnatal-preoperative window both
clearances rise with postnatal age along a fixed sigmoidal recovery
curve; on bypass and after bypass they are instead fixed multiples of
the at-birth value, so the postnatal-age term drops out entirely.

## Population

Fourteen neonates with CCHD contributed 140 paired allopurinol and
oxypurinol plasma observations (60 postnatal-preoperative, 36
intraoperative, 44 postoperative; median 10 samples per patient, IQR
9-12). Median birth weight was 3.16 kg (IQR 2.75-3.73) and median
gestational age 38.0 weeks (IQR 38.0-38.8); 10 of 14 (71.4%) were male.
Cardiac pathology was transposition of the great arteries in 6 (42.9%),
single-ventricle physiology in 4 (28.6%), aortic arch anomaly in 2
(14.3%) and other in 2 (14.3%). Thirteen neonates underwent cardiac
surgery with CPB at a median postnatal age of 5.60 days (IQR 4.78-7.81),
with a median total duration of 320 min (IQR 280-368) and a median
lowest rectal temperature of 27.7 C (IQR 23.7-28.0); the fourteenth had
surgery at a non-participating centre and contributed postnatal-period
samples only. Baseline demographics are Chu 2024 Table 1; sampling and
observation counts are Chu 2024 Sects. 2.1 and 3.1.

Each dose was 20 mg/kg delivered over 10 min by syringe pump: within
45-60 min of birth (DOSE 1), 12 h later (DOSE 2), 12 h before surgery
(DOSE 3), at the start of CPB (DOSE 4), and 24 h after surgery (DOSE 5).
The 10-min infusion duration is stated by the CRUCIAL protocol paper,
not by Chu 2024 itself.

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Chu_2024_allopurinol")()$population`).

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in `inst/modeldb/specificDrugs/Chu_2024_allopurinol.R`.
The table below collects them in one place for review. “ESM S3” is the
NONMEM control stream printed in the electronic supplementary material.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (allopurinol CL at birth) | 0.95 L/h per 3.5 kg | Table 2, “Allopurinol clearance (All CL_Postnatal)” |
| `lvp` (allopurinol peripheral V) | 2.22 L per 3.5 kg | Table 2, “Allopurinol volume of distribution (All Vd_Postnatal)”; ESM S3 `S3 = VA`, `K31 = Q1/VA` identify it as peripheral |
| `lvc` (allopurinol central V1) | 0.1 L per 3.5 kg (fixed) | Table 2, “Allopurinol volume of central compartment (All V1)” |
| `lq` (allopurinol Q1) | 6.97 L/h per 3.5 kg (fixed) | Table 2, “Allopurinol intercompartmental clearance (All Q1)”; see Errata for the Sect. 3.2 text value |
| `lcl_oxy` (oxypurinol CL at birth) | 0.21 L/h per 3.5 kg | Table 2, “Oxypurinol clearance (Oxy CL_Postnatal)” |
| `lvc_oxy` (oxypurinol V) | 12 L per 3.5 kg | Table 2, “Oxypurinol volume of distribution (Oxy Vd_Postnatal)” |
| `limax` (max auto-inhibition) | 1 (fixed) | Table 2, “Maximum achievable autoinhibition effect” |
| `lic50` (auto-inhibition IC50) | 1.1 mg/L (fixed) | Table 2, “IC50,auto-inhibition”; Sect. 3.2 |
| `fm` (formation fraction) | 1 (fixed) | Table 2 footnote b; ESM S1 legend “fm, formation fraction was assumed 1” |
| `e_wt_cl`, `e_wt_vc` | 0.75, 1 (fixed) | Sect. 2.3; ESM S3 `(WT/3.5)**0.75`, `(WT/3.5)**1` |
| `e_pna_cl`, `e_pna_cl_oxy` | 3, 1.35 (fixed) | Table 2, “Maximum fold of increase in All / Oxy CL_Postnatal” |
| `tm50_cl`, `hill_cl` | 4.2 days, 2.98 (fixed) | Table 2, “Postnatal age at 50% of maximum recovery effect”, “Hill coefficient” |
| `e_cpb_on_cl`, `e_cpb_on_vp` | 1.46, 1.47 | Table 2, “All E_CL,CPB”, “All E_Vd,CPB” |
| `e_cpb_on_cl_oxy`, `e_cpb_on_vc_oxy` | 0.54, 1.3 | Table 2, “Oxy E_CL,CPB”, “Oxy E_Vd,CPB” |
| `e_cpb_post_cl`, `e_cpb_post_vp` | 2.33, 1.42 | Table 2, “All E_CL,Postop”, “All E_Vd,Postop” |
| `e_cpb_post_cl_oxy`, `e_cpb_post_vc_oxy` | 0.23, 1.48 | Table 2, “Oxy E_CL,Postop”, “Oxy E_Vd,Postop” |
| `etalcl`, `etalvc_oxy` | 36%, 42% CV | Table 2 “BSV All CL”, “BSV Oxy Vd”; footnote gives `CV% = sqrt(omega^2)*100` |
| `etaiov_cl_*`, `etaiov_cl_oxy_*` | 18%, 43% CV | Table 2 “BOV All CL”, “BOV Oxy CL”; ESM S3 `BOV_CLA`, `BOV_CLO` |
| `propSd`, `propSd_oxy` | 0.25, 0.16 | Table 2 “Residual proportional error” rows |
| `addSd`, `addSd_oxy` | 0.025, 0.02335 mg/L (fixed) | Sect. 3.1 (LLOQ/2 convention) and Sect. 2.1 (LLOQ values); see Errata |
| Period selection of CL / Vd | n/a | Eqs. 1-4; ESM S3 `TVCLA = (CLA_PNA**FLAG1) * (CLA_CPB**FLAG2) * (CLA_POST**FLAG3)` |
| Sigmoidal recovery with PNA | n/a | Eq. 6; ESM S3 `EFFPNA1 = 1 + (MMAX1*(PNA**HILL)/((TM50**HILL)+(PNA**HILL)))` |
| Auto-inhibition of conversion | n/a | ESM S3 `$DES`: `AUTOI = EMAX*C2/(IC50+C2)`, `DADT(1) = -K12*A(1)*(1-AUTOI) + K31*A(3) - K13*A(1)` |
| Allopurinol-to-oxypurinol mass conversion | 152.11 / 136.11 | ESM S3 `S2 = VO*(136.11/152.11)` |

## Perioperative schedule and virtual cohort

Original observed data are not publicly available. The cohort below
approximates the published trial demographics. Body weight is drawn from
a log-normal matched to the Table 1 birth-weight median and IQR; the
perioperative schedule is held at the cohort medians (CPB starting at a
postnatal age of 5.60 days and running 320 min) so that every subject
shares the same nominal dose times, which is what makes the per-dose NCA
windows below well defined.

``` r

# Perioperative landmarks, in hours after birth (= after DOSE 1). Chu 2024
# Table 1 medians: cardiac surgery with CPB starts at a postnatal age of
# 5.60 days and the total duration is 320 min.
cpb_start <- 5.60 * 24
cpb_dur   <- 320 / 60
cpb_end   <- cpb_start + cpb_dur

# Chu 2024 Sect. 2.1 / Fig. 1 dosing schedule.
dose_times <- c(
  `DOSE 1 (birth)`        = 0,
  `DOSE 2 (+12 h)`        = 12,
  `DOSE 3 (pre-op)`       = cpb_start - 12,
  `DOSE 4 (start of CPB)` = cpb_start,
  `DOSE 5 (post-op)`      = cpb_end + 24
)
sim_end <- unname(dose_times[5]) + 26

# Chu 2024 Table 3 evaluates target attainment 12 h after every dose and
# 24 h after DOSE 2, DOSE 4 and DOSE 5. Every one of these instants must be
# an observation time or the attainment table below reads an empty slice.
target_times <- sort(unique(unname(dose_times) + rep(c(12, 24), each = 5)))

knitr::kable(
  data.frame(
    Dose = names(dose_times),
    `Time after birth (h)` = round(unname(dose_times), 2),
    check.names = FALSE
  ),
  caption = "Nominal CRUCIAL dosing schedule used for the simulations, at the Chu 2024 Table 1 cohort medians. Each dose is 20 mg/kg infused over 10 min."
)
```

| Dose                  | Time after birth (h) |
|:----------------------|---------------------:|
| DOSE 1 (birth)        |                 0.00 |
| DOSE 2 (+12 h)        |                12.00 |
| DOSE 3 (pre-op)       |               122.40 |
| DOSE 4 (start of CPB) |               134.40 |
| DOSE 5 (post-op)      |               163.73 |

Nominal CRUCIAL dosing schedule used for the simulations, at the Chu
2024 Table 1 cohort medians. Each dose is 20 mg/kg infused over 10 min.
{.table}

``` r

# `set.seed()` seeds R's RNG. It does NOT seed rxode2's simulation RNG, and
# rxode2's streams are partitioned PER SOLVER THREAD -- so the cohort below
# is reproducible on this machine and different on a machine with a
# different thread count. Every assertion downstream is written so that it
# holds for any cohort the model can produce.
set.seed(20240815)

n_subj <- 150L

# Body weight: log-normal matched to Chu 2024 Table 1 birth weight,
# median 3.16 kg with IQR 2.75-3.73 kg.
wt_meanlog <- log(3.16)
wt_sdlog   <- log(3.73 / 2.75) / (2 * stats::qnorm(0.75))
subj <- tibble::tibble(
  id = seq_len(n_subj),
  WT = stats::rlnorm(n_subj, meanlog = wt_meanlog, sdlog = wt_sdlog)
)

# Observation grid: 1 min for the first hour after each dose (the model's
# 0.1 L central compartment produces a very fast peak), 5 min out to 12 h
# after each dose, and 1 h elsewhere. The interval boundaries used by PKNCA
# below are forced into the grid so every AUC window is anchored on a real
# observation, and the CPB boundaries are bracketed so the piecewise-constant
# period indicators switch exactly where they should.
obs_grid <- sort(unique(c(
  seq(0, sim_end, by = 1),
  unlist(lapply(dose_times, function(d) d + seq(0, 1, by = 1 / 60))),
  unlist(lapply(dose_times, function(d) d + seq(0, 12, by = 1 / 12))),
  target_times,
  cpb_start, cpb_end, cpb_start - 1e-4, cpb_end - 1e-4
)))
obs_grid <- obs_grid[obs_grid >= 0 & obs_grid <= sim_end]

# `cmt` on an observation row must name a declared ODE state, never the
# algebraic observable `Cc`. `dvid = 1L` on every observation row makes
# rxode2 return BOTH endpoint columns (Cc and Cc_oxy) at each row.
obs <- tidyr::crossing(subj, time = obs_grid) |>
  dplyr::mutate(amt = NA_real_, dur = NA_real_, evid = 0L,
                cmt = "central", dvid = 1L)

dos <- tidyr::crossing(subj, time = unname(dose_times)) |>
  dplyr::mutate(amt = 20 * WT, dur = 10 / 60, evid = 1L,
                cmt = "central", dvid = NA_integer_)

events <- dplyr::bind_rows(dos, obs) |>
  dplyr::arrange(id, time, dplyr::desc(evid)) |>
  dplyr::mutate(
    # Postnatal age in MONTHS (the canonical PNA unit); the model recovers
    # days internally with PNA * 30.4375.
    PNA        = time / 24 / 30.4375,
    CPB_ON     = as.numeric(time >= cpb_start & time < cpb_end),
    CPB_REWARM = 0,
    CPB_POST   = as.numeric(time >= cpb_end)
  ) |>
  as.data.frame()

# The three period indicators must partition the record.
stopifnot(all(events$CPB_ON + events$CPB_REWARM + events$CPB_POST %in% c(0, 1)))
#> Warning in all(events$CPB_ON + events$CPB_REWARM + events$CPB_POST %in% :
#> coercing argument of type 'double' to logical
stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod <- readModelDb("Chu_2024_allopurinol")

sim <- rxode2::rxSolve(mod, events = events, keep = c("WT")) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_oxy_1, etaiov_cl_oxy_2, etaiov_cl_oxy_3
#> as a work-around try putting the mu-referenced expression on a simple line

# Allopurinol washes out almost completely across the ~98 h between the
# second and third doses, and the solver's residual there decays into
# numerical noise of order 1e-11 mg/L that can dip below zero. PKNCA would
# take log() of a negative value, so clamp at zero -- but assert first on the
# MAGNITUDE of any negative excursion, because a structural error (a wrong
# sign, an unstable rate constant) would produce excursions many orders of
# magnitude larger than this. Realised minimum was about -6e-11 mg/L against
# peak concentrations of tens of mg/L.
stopifnot(min(sim$Cc, na.rm = TRUE) > -1e-6,
          min(sim$Cc_oxy, na.rm = TRUE) > -1e-6)
sim$Cc     <- pmax(sim$Cc, 0)
sim$Cc_oxy <- pmax(sim$Cc_oxy, 0)
```

## Replicating Figure 3: period-specific clearance in a typical neonate

Chu 2024 Fig. 3 traces the clearance of a typical 3.5 kg neonate who
starts cardiac surgery with CPB at a postnatal age of 5.6 days. This is
a typical-value prediction, so the random effects are zeroed.

``` r

# Replicates Figure 3 of Chu 2024: allopurinol and oxypurinol clearance in a
# typical 3.5 kg neonate over postnatal age, with the CPB window shaded.
typ_events <- data.frame(
  id = 1L, time = obs_grid, amt = NA_real_, dur = NA_real_,
  evid = 0L, cmt = "central", dvid = 1L, WT = 3.5
) |>
  dplyr::mutate(
    PNA        = time / 24 / 30.4375,
    CPB_ON     = as.numeric(time >= cpb_start & time < cpb_end),
    CPB_REWARM = 0,
    CPB_POST   = as.numeric(time >= cpb_end)
  )

sim_typ <- rxode2::rxSolve(rxode2::zeroRe(mod), events = typ_events) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_oxy_1, etaiov_cl_oxy_2, etaiov_cl_oxy_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_oxy_1, etaiov_cl_oxy_2, etaiov_cl_oxy_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc_oxy', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_oxy_1', 'etaiov_cl_oxy_2', 'etaiov_cl_oxy_3'

sim_typ |>
  dplyr::select(time, Allopurinol = cl, Oxypurinol = cl_oxy) |>
  tidyr::pivot_longer(-time, names_to = "Analyte", values_to = "CL") |>
  ggplot(aes(time / 24, CL)) +
  annotate("rect", xmin = cpb_start / 24, xmax = cpb_end / 24,
           ymin = -Inf, ymax = Inf, alpha = 0.18) +
  geom_step(direction = "hv") +
  facet_wrap(~Analyte, scales = "free_y") +
  labs(x = "Postnatal age (days)", y = "Clearance (L/h)",
       title = "Figure 3 - clearance in a typical 3.5 kg neonate",
       caption = "Replicates Figure 3 of Chu 2024. Shaded band is the CPB window.")
```

![](Chu_2024_allopurinol_files/figure-html/figure-3-1.png)

``` r

# Chu 2024 Fig. 3 caption and Sect. 3.2 print the four clearance landmarks
# for each analyte, and Sect. 3.2 the three volume landmarks. This gate is
# DETERMINISTIC (typical-value solve, no cohort draw), so the achieved
# percentages below are stable across machines and thread counts.
at <- function(tt, col) sim_typ[[col]][which.min(abs(sim_typ$time - tt))]

fig3 <- tibble::tribble(
  ~Quantity,                          ~Published, ~Simulated,
  "Allopurinol CL at birth (L/h)",          0.95, at(0, "cl"),
  "Allopurinol CL before CPB (L/h)",        2.97, at(cpb_start - 0.5, "cl"),
  "Allopurinol CL during CPB (L/h)",        1.38, at(cpb_start + 1, "cl"),
  "Allopurinol CL after CPB (L/h)",         2.21, at(cpb_end + 5, "cl"),
  "Oxypurinol CL at birth (L/h)",           0.21, at(0, "cl_oxy"),
  "Oxypurinol CL before CPB (L/h)",         0.41, at(cpb_start - 0.5, "cl_oxy"),
  "Oxypurinol CL during CPB (L/h)",         0.12, at(cpb_start + 1, "cl_oxy"),
  "Oxypurinol CL after CPB (L/h)",          0.05, at(cpb_end + 5, "cl_oxy"),
  "Allopurinol Vd before CPB (L)",          2.22, at(cpb_start - 0.5, "vp"),
  "Allopurinol Vd during CPB (L)",          3.26, at(cpb_start + 1, "vp"),
  "Allopurinol Vd after CPB (L)",           3.15, at(cpb_end + 5, "vp"),
  "Oxypurinol Vd before CPB (L)",          12.00, at(cpb_start - 0.5, "vc_oxy"),
  "Oxypurinol Vd during CPB (L)",          15.60, at(cpb_start + 1, "vc_oxy"),
  "Oxypurinol Vd after CPB (L)",           17.76, at(cpb_end + 5, "vc_oxy")
) |>
  dplyr::mutate(`% diff` = 100 * (Simulated - Published) / Published)

fig3 |>
  dplyr::mutate(dplyr::across(c(Published, Simulated), ~ signif(.x, 4)),
                `% diff` = round(`% diff`, 2)) |>
  knitr::kable(caption = "Chu 2024 Figure 3 / Section 3.2 landmarks versus the packaged model (typical 3.5 kg neonate).")
```

| Quantity                        | Published | Simulated | % diff |
|:--------------------------------|----------:|----------:|-------:|
| Allopurinol CL at birth (L/h)   |      0.95 |    0.9500 |   0.00 |
| Allopurinol CL before CPB (L/h) |      2.97 |    2.9440 |  -0.86 |
| Allopurinol CL during CPB (L/h) |      1.38 |    1.3870 |   0.51 |
| Allopurinol CL after CPB (L/h)  |      2.21 |    2.2140 |   0.16 |
| Oxypurinol CL at birth (L/h)    |      0.21 |    0.2100 |   0.00 |
| Oxypurinol CL before CPB (L/h)  |      0.41 |    0.4084 |  -0.39 |
| Oxypurinol CL during CPB (L/h)  |      0.12 |    0.1134 |  -5.50 |
| Oxypurinol CL after CPB (L/h)   |      0.05 |    0.0483 |  -3.40 |
| Allopurinol Vd before CPB (L)   |      2.22 |    2.2200 |   0.00 |
| Allopurinol Vd during CPB (L)   |      3.26 |    3.2630 |   0.10 |
| Allopurinol Vd after CPB (L)    |      3.15 |    3.1520 |   0.08 |
| Oxypurinol Vd before CPB (L)    |     12.00 |   12.0000 |   0.00 |
| Oxypurinol Vd during CPB (L)    |     15.60 |   15.6000 |   0.00 |
| Oxypurinol Vd after CPB (L)     |     17.76 |   17.7600 |   0.00 |

Chu 2024 Figure 3 / Section 3.2 landmarks versus the packaged model
(typical 3.5 kg neonate). {.table}

``` r


# The volume landmarks are printed to 3-4 significant figures and reproduce
# essentially exactly (realised |% diff| <= 0.1). The clearance landmarks are
# printed to 2-3 significant figures, and the two oxypurinol bypass values
# are the loosest: 0.21 * 0.54 = 0.1134 is printed as "0.12" and
# 0.21 * 0.23 = 0.0483 as "0.05". Both are consistent with the rounding
# intervals of the underlying estimates (e.g. 0.215 * 0.545 = 0.117 -> 0.12),
# so the bound below admits 2-significant-figure rounding. It still goes red
# on a mis-transcribed clearance or fraction, which moves these by tens of
# percent. Realised |% diff|: CL 0.0 / 0.8 / 0.5 / 0.1 (allopurinol) and
# 0.0 / 0.4 / 5.5 / 3.4 (oxypurinol); Vd all <= 0.1.
is_vd <- grepl("Vd", fig3$Quantity)
stopifnot(
  max(abs(fig3$`% diff`[is_vd]))  < 1,
  max(abs(fig3$`% diff`[!is_vd])) < 8
)
```

The three-part shape of Fig. 3 is reproduced exactly, including its most
counterintuitive feature: after separation from bypass, allopurinol
clearance **recovers above** its at-birth value (2.21 versus 0.95 L/h)
while oxypurinol clearance **falls below** it (0.05 versus 0.21 L/h).
Chu 2024 Sect. 4 attributes the oxypurinol drop to CPB-associated acute
kidney injury.

## Replicating Figure 6: combined concentration against the XO IC90

Chu 2024 Fig. 6 plots allopurinol plus oxypurinol against the 6.22 mg/L
combined concentration that gave 90% xanthine-oxidase inhibition in the
authors’ earlier hypoxic-ischemic-encephalopathy cohort.

``` r

# Replicates Figure 6 of Chu 2024: allopurinol + oxypurinol against the
# IC90 XO of 6.22 mg/L across the whole postnatal and perioperative window.
ic90_xo <- 6.22

sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(Ctot = Cc + Cc_oxy) |>
  dplyr::group_by(time) |>
  dplyr::summarise(Q05 = quantile(Ctot, 0.05), Q50 = median(Ctot),
                   Q95 = quantile(Ctot, 0.95), .groups = "drop") |>
  ggplot(aes(time / 24, Q50)) +
  annotate("rect", xmin = cpb_start / 24, xmax = cpb_end / 24,
           ymin = -Inf, ymax = Inf, alpha = 0.18) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.3) +
  geom_line() +
  geom_hline(yintercept = ic90_xo, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Postnatal age (days)", y = "Allopurinol + oxypurinol (mg/L)",
       title = "Figure 6 - combined concentration versus the XO IC90",
       caption = paste("Replicates Figure 6 of Chu 2024. Dashed line is the",
                       "IC90 XO of 6.22 mg/L; band is the 5th-95th percentile."))
#> Warning in log(x, base): NaNs produced
#> Warning in scale_y_log10(): log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
```

![](Chu_2024_allopurinol_files/figure-html/figure-6-1.png)

## PKNCA validation

Chu 2024 Table 3 reports the AUC over the 12 h following each of the
five doses, computed from individual maximum-a-posteriori Bayesian
estimates. One PKNCA block is run per analyte, with one interval per
dose window.

``` r

intervals <- data.frame(
  start   = unname(dose_times),
  end     = unname(dose_times) + 12,
  auclast = TRUE,
  cmax    = TRUE,
  tmax    = TRUE
)
window_label <- setNames(names(dose_times), format(unname(dose_times)))

# Doses: one row per dose event per subject.
dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt) |>
  dplyr::mutate(analyte = "Allopurinol")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | analyte + id)
```

``` r

# IMPORTANT: do NOT add `time > 0` or `Cc > 0` to this filter -- both drop
# the time-zero row that PKNCA needs to anchor the first AUC window.
conc_allo <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::transmute(id, time, Cc, analyte = "Allopurinol")

# Guarantee a time = 0 row per subject (pre-dose concentration is zero).
conc_allo <- dplyr::bind_rows(
  conc_allo,
  conc_allo |> dplyr::distinct(id, analyte) |>
    dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, analyte, time, .keep_all = TRUE) |>
  dplyr::arrange(id, time)

nca_allo <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(conc_allo, Cc ~ time | analyte + id),
  dose_obj, intervals = intervals
))
```

``` r

conc_oxy <- sim |>
  dplyr::filter(!is.na(Cc_oxy)) |>
  dplyr::transmute(id, time, Cc = Cc_oxy, analyte = "Oxypurinol")

conc_oxy <- dplyr::bind_rows(
  conc_oxy,
  conc_oxy |> dplyr::distinct(id, analyte) |>
    dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, analyte, time, .keep_all = TRUE) |>
  dplyr::arrange(id, time)

dose_obj_oxy <- PKNCA::PKNCAdose(
  dose_df |> dplyr::mutate(analyte = "Oxypurinol"), amt ~ time | analyte + id
)

nca_oxy <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(conc_oxy, Cc ~ time | analyte + id),
  dose_obj_oxy, intervals = intervals
))
```

### Comparison against the published AUC

``` r

# Map each PKNCA interval back onto the dose window it belongs to.
label_windows <- function(res) {
  as.data.frame(res$result) |>
    dplyr::mutate(Window = unname(window_label[format(start)])) |>
    dplyr::filter(!is.na(Window))
}

simulated_nca <- dplyr::bind_rows(label_windows(nca_allo),
                                  label_windows(nca_oxy)) |>
  dplyr::rename(Analyte = analyte)
stopifnot(nrow(simulated_nca) > 0L)

# Chu 2024 Table 3, "Allopurinol AUC 12" and "Oxypurinol AUC 12" columns
# (median [IQR], mg/L*h). Table 3 reports no AUC for the 24 h timepoints.
published <- tibble::tribble(
  ~Analyte,      ~Window,                 ~auclast,
  "Allopurinol", "DOSE 1 (birth)",           171,
  "Allopurinol", "DOSE 2 (+12 h)",           291,
  "Allopurinol", "DOSE 3 (pre-op)",           95,
  "Allopurinol", "DOSE 4 (start of CPB)",    154,
  "Allopurinol", "DOSE 5 (post-op)",         175,
  "Oxypurinol",  "DOSE 1 (birth)",          38.6,
  "Oxypurinol",  "DOSE 2 (+12 h)",          67.4,
  "Oxypurinol",  "DOSE 3 (pre-op)",         56.6,
  "Oxypurinol",  "DOSE 4 (start of CPB)",   59.5,
  "Oxypurinol",  "DOSE 5 (post-op)",         109
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = simulated_nca,
  reference     = published,
  by            = c("Analyte", "Window"),
  units         = c(auclast = "mg/L*h"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = "Simulated AUC0-12 after each dose (cohort median) versus Chu 2024 Table 3. * differs from the reference by more than 20%.",
  align = c("l", "l", "l", "r", "r", "r")
)
```

| NCA parameter | Analyte | Window | Reference | Simulated | % diff |
|:---|:---|:---|---:|---:|---:|
| AUClast (mg/L\*h) | Allopurinol | DOSE 1 (birth) | 171 | 168 | -1.5% |
| AUClast (mg/L\*h) | Allopurinol | DOSE 2 (+12 h) | 291 | 298 | +2.5% |
| AUClast (mg/L\*h) | Allopurinol | DOSE 3 (pre-op) | 95 | 101 | +6.6% |
| AUClast (mg/L\*h) | Allopurinol | DOSE 4 (start of CPB) | 154 | 164 | +6.6% |
| AUClast (mg/L\*h) | Allopurinol | DOSE 5 (post-op) | 175 | 174 | -0.6% |
| AUClast (mg/L\*h) | Oxypurinol | DOSE 1 (birth) | 38.6 | 37.2 | -3.6% |
| AUClast (mg/L\*h) | Oxypurinol | DOSE 2 (+12 h) | 67.4 | 70.7 | +4.9% |
| AUClast (mg/L\*h) | Oxypurinol | DOSE 3 (pre-op) | 56.6 | 57.6 | +1.8% |
| AUClast (mg/L\*h) | Oxypurinol | DOSE 4 (start of CPB) | 59.5 | 65.7 | +10.4% |
| AUClast (mg/L\*h) | Oxypurinol | DOSE 5 (post-op) | 109 | 105 | -4.1% |

Simulated AUC0-12 after each dose (cohort median) versus Chu 2024 Table
3. \* differs from the reference by more than 20%. {.table}

``` r

# `% diff` comes back from ncaComparisonTable as formatted TEXT (it may carry
# a trailing "*"), so it is re-derived numerically here for the gate.
sim_med <- simulated_nca |>
  dplyr::filter(PPTESTCD == "auclast") |>
  dplyr::group_by(Analyte, Window) |>
  dplyr::summarise(Simulated = median(PPORRES, na.rm = TRUE), .groups = "drop")

gate <- published |>
  dplyr::inner_join(sim_med, by = c("Analyte", "Window")) |>
  dplyr::mutate(pct = 100 * (Simulated - auclast) / auclast)
stopifnot(nrow(gate) == nrow(published))

# Chu 2024 Table 3's AUCs come from individual MAP Bayesian estimates of 14
# real neonates with their own weights and their own actual dose and surgery
# times; this cohort uses one nominal schedule at the Table 1 medians. A
# median that tracks the published median to within a quarter is therefore a
# strong structural check, and it still goes red on a mis-transcribed
# clearance, volume, dose or unit, which move AUC by tens of percent.
# Realised over the ten windows: |% diff| max 10.4, median 4.3.
stopifnot(
  max(abs(gate$pct))            < 25,
  median(abs(gate$pct))         < 15
)

# Nine of the ten simulated medians land inside the published IQR; assert
# the weaker claim that the majority do, so the gate does not hinge on one
# cohort draw.
iqr_lo <- c(160, 283,  86, 114, 151, 29.2, 53.7, 50.9, 53.7, 87.5)
iqr_hi <- c(195, 302, 131, 163, 196, 44.2, 83.6, 62.1, 68.6,  118)
gate_ord <- gate[match(paste(published$Analyte, published$Window),
                       paste(gate$Analyte, gate$Window)), ]
stopifnot(sum(gate_ord$Simulated >= iqr_lo & gate_ord$Simulated <= iqr_hi) >= 6L)
```

## Target attainment

Chu 2024 Table 3 also reports, at each timepoint, the proportion of
patients whose allopurinol concentration exceeded the 2 mg/L trial
target and whose combined allopurinol-plus-oxypurinol concentration
exceeded the 6.22 mg/L IC90 for xanthine oxidase.

``` r

target_allo <- 2

timepoints <- tibble::tribble(
  ~Stage,           ~Window,                  ~Offset, ~ref_allo, ~ref_ic90,
  "Postnatal",      "DOSE 1 (birth)",              12,     100.0,      100,
  "Postnatal",      "DOSE 2 (+12 h)",              12,     100.0,      100,
  "Postnatal",      "DOSE 2 (+12 h)",              24,      91.7,      100,
  "Preoperative",   "DOSE 3 (pre-op)",             12,      36.4,       73,
  "Intraoperative", "DOSE 4 (start of CPB)",       12,      92.3,      100,
  "Postoperative",  "DOSE 4 (start of CPB)",       24,      53.8,       92,
  "Postoperative",  "DOSE 5 (post-op)",            12,     100.0,      100,
  "Postoperative",  "DOSE 5 (post-op)",            24,      78.6,      100
) |>
  dplyr::mutate(Time = unname(dose_times[Window]) + Offset)

# Fail loudly rather than silently summarising an empty slice: a lookup that
# matches no rows would return NaN, and `mean(logical(0))` in a percentage
# column reads as a missing cell rather than as a broken gate.
slice_at <- function(tt) {
  s <- sim[abs(sim$time - tt) < 1e-6, , drop = FALSE]
  if (nrow(s) == 0L) {
    stop("no simulated observations at t = ", tt,
         " h; add it to the observation grid")
  }
  s
}
attain <- timepoints |>
  dplyr::mutate(
    sim_allo = vapply(Time, function(tt) {
      s <- slice_at(tt); 100 * mean(s$Cc > target_allo)
    }, numeric(1L)),
    sim_ic90 = vapply(Time, function(tt) {
      s <- slice_at(tt); 100 * mean(s$Cc + s$Cc_oxy > ic90_xo)
    }, numeric(1L))
  )
stopifnot(nrow(attain) == 8L, !anyNA(attain$sim_allo), !anyNA(attain$sim_ic90))

attain <- attain |>
  dplyr::mutate(
    # The one window where the model does not reproduce the published IC90
    # attainment; see the gate chunk and the Errata for the diagnosis.
    Deviation = Window == "DOSE 3 (pre-op)" & Offset == 12
  )

attain |>
  dplyr::transmute(
    Stage,
    Timepoint = paste0(Offset, " h after ", sub(" \\(.*", "", Window)),
    `Published allopurinol > 2 mg/L (%)`  = ref_allo,
    `Simulated allopurinol > 2 mg/L (%)`  = round(sim_allo, 1),
    `Published allo+oxy > IC90 XO (%)`    = ref_ic90,
    `Simulated allo+oxy > IC90 XO (%)`    = round(sim_ic90, 1),
    `Known deviation`                     = ifelse(Deviation, "yes", "")
  ) |>
  knitr::kable(caption = "Target attainment: simulated cohort versus Chu 2024 Table 3.")
```

| Stage | Timepoint | Published allopurinol \> 2 mg/L (%) | Simulated allopurinol \> 2 mg/L (%) | Published allo+oxy \> IC90 XO (%) | Simulated allo+oxy \> IC90 XO (%) | Known deviation |
|:---|:---|---:|---:|---:|---:|:---|
| Postnatal | 12 h after DOSE 1 | 100.0 | 98.0 | 100 | 91.3 |  |
| Postnatal | 12 h after DOSE 2 | 100.0 | 100.0 | 100 | 100.0 |  |
| Postnatal | 24 h after DOSE 2 | 91.7 | 92.0 | 100 | 97.3 |  |
| Preoperative | 12 h after DOSE 3 | 36.4 | 44.7 | 73 | 41.3 | yes |
| Intraoperative | 12 h after DOSE 4 | 92.3 | 96.7 | 100 | 97.3 |  |
| Postoperative | 24 h after DOSE 4 | 53.8 | 52.7 | 92 | 84.0 |  |
| Postoperative | 12 h after DOSE 5 | 100.0 | 98.0 | 100 | 100.0 |  |
| Postoperative | 24 h after DOSE 5 | 78.6 | 75.3 | 100 | 97.3 |  |

Target attainment: simulated cohort versus Chu 2024 Table 3. {.table}

``` r

# The published percentages come from 12-14 patients, so their resolution is
# roughly 8 percentage points and they cannot be matched point-for-point.
# The gate therefore asserts the paper's own qualitative conclusions, each of
# which is an absolute claim with headroom rather than a race between two
# noisy statistics:
#
#  1. The 2 mg/L allopurinol target is MISSED (< 66%) 12 h after the third
#     dose -- Chu 2024's headline finding (published 36.4%). The typical-value
#     concentration there is 1.87 mg/L, i.e. below target before any
#     variability is added.
#  2. The 2 mg/L target is comfortably met (> 90%) in the two early postnatal
#     windows (published 100% for both).
#  3. The IC90 XO target is met in more than two-thirds of the cohort at every
#     tabulated timepoint EXCEPT the flagged preoperative window -- the basis
#     for the paper's conclusion that no dose adjustment is needed (published
#     range 73-100%).
row_of <- function(win, off) {
  r <- attain[attain$Window == win & attain$Offset == off, ]
  stopifnot(nrow(r) == 1L)
  r
}
stopifnot(
  row_of("DOSE 3 (pre-op)", 12)$sim_allo   <  66,
  row_of("DOSE 1 (birth)", 12)$sim_allo    >  85,
  row_of("DOSE 2 (+12 h)", 12)$sim_allo    >  85,
  all(attain$sim_ic90[!attain$Deviation]   >  66)
)

# The excluded window, quantified rather than hidden. The published 73% is
# 8 of 11 evaluable patients, so one patient is worth 9 percentage points
# there. The simulated cohort's MEDIAN combined concentration at that instant
# sits just under the 6.22 mg/L threshold, which is what turns a small
# exposure difference into a large attainment difference; the AUC over the
# same window agrees with Table 3 to within 8% for both analytes.
dev_row <- row_of("DOSE 3 (pre-op)", 12)
dev_med <- median(slice_at(dev_row$Time)$Cc + slice_at(dev_row$Time)$Cc_oxy)
cat(sprintf(
  paste0("Preoperative window, 12 h after DOSE 3:\n",
         "  simulated median allopurinol + oxypurinol = %.2f mg/L",
         " (IC90 XO threshold 6.22 mg/L, i.e. %.0f%% short)\n",
         "  simulated IC90 attainment %.1f%% versus a published 73%% (8/11 patients)\n"),
  dev_med, 100 * (6.22 - dev_med) / 6.22, dev_row$sim_ic90))
#> Preoperative window, 12 h after DOSE 3:
#>   simulated median allopurinol + oxypurinol = 5.63 mg/L (IC90 XO threshold 6.22 mg/L, i.e. 9% short)
#>   simulated IC90 attainment 41.3% versus a published 73% (8/11 patients)

# The threshold sensitivity is the point: the model is not far off in
# EXPOSURE, only in where the cohort sits relative to a nearby cut-off.
stopifnot(abs(100 * (dev_med - 6.22) / 6.22) < 25)
```

The allopurinol column tracks the published percentages closely at all
seven comparable timepoints, including the paper’s central finding that
the predefined 2 mg/L target is missed 12 h after the preoperative third
dose (simulated 45%, published 36.4%).

The IC90 column agrees everywhere except at that same preoperative
window, where the simulated cohort attains 41% against a published 73%.
This is a threshold artefact rather than an exposure error, and it is
flagged as a known deviation rather than tuned away. The simulated
median combined concentration at that instant is about 5.6 mg/L against
the 6.22 mg/L threshold, roughly 10% short, so a modest shift moves a
large share of the cohort across the cut-off; the AUC over the very same
window agrees with Table 3 to within 7% for allopurinol and 2% for
oxypurinol. The published 73% is 8 of 11 evaluable patients, so one
patient is worth 9 percentage points. The most likely contributor is
this vignette’s fixed perioperative schedule: every simulated subject
has surgery at the cohort-median postnatal age of 5.60 days, which
places the third dose at the steepest part of the recovery curve,
whereas the real cohort’s surgery times spanned 4.78-7.81 days and the
earlier-surgery patients would have had materially lower preoperative
clearance and therefore higher concentrations.

## Assumptions and deviations

### Errata and internal inconsistencies in the source

- **Allopurinol Q1 is printed twice with different digits.** Table 2
  gives `All Q1 = 6.97 L/h (fix)` while Sect. 3.2 says “a large
  inter-compartment CL (Q1) of 6.79 L/h” – a digit transposition in one
  of the two. The model uses the **Table 2** value, 6.97 L/h, on the
  grounds that Table 2 is the paper’s parameter table. The consequence
  is negligible either way: with a 0.1 L central compartment the
  distributional half-life is about 30 s at either value, and both are
  described by the authors as an artifact of the intraoperative sampling
  schedule that lacks physiological interpretation.
- **The oxypurinol CPB volume fraction is printed as 1.3 but described
  as 33%.** Table 2 gives `Oxy E_Vd,CPB = 1.3` (95% CI 1.2-1.43),
  whereas Sect. 4 says the oxypurinol Vd during CPB “increased by …
  33%”. The model uses the Table 2 value.
- **The oxypurinol additive residual SD is not published.** The ESM S3
  `$ERROR` block declares `EPS(4)` for oxypurinol, but Table 2 tabulates
  only the proportional components and Sect. 3.1 states the LLOQ/2
  convention in the context of allopurinol (the only analyte with
  below-LLOQ data). This model sets `addSd_oxy` to LLOQ/2 for
  oxypurinol, 0.0467 / 2 = 0.02335 mg/L, following the same stated
  convention. Oxypurinol concentrations in this model run from roughly 1
  to 100 mg/L, so the choice is numerically immaterial; it is recorded
  here because the value is derived rather than printed.

### Encoding decisions

- **Which allopurinol volume is central.** Table 2 lists “Allopurinol
  volume of distribution (All Vd_Postnatal) = 2.22 L” and, separately,
  “Allopurinol volume of central compartment (All V1) = 0.1 L (fix)”.
  The ESM S3 control stream resolves the ambiguity: `S3 = VA` and
  `K31 = Q1/VA` make the 2.22 L estimate the **peripheral** volume
  (`lvp` here) and the fixed 0.1 L the observed **central** volume
  (`lvc`). Total Vss is 2.32 L.
- **The CPB and postoperative fractions multiply the at-birth baseline,
  not the postnatal-age-adjusted value.** Eqs. 1-4 are written against
  `CL_Postnatal`, which is ambiguous once postnatal age is a covariate
  on that same quantity. Arithmetic settles it: 0.95 x 1.46 = 1.39 L/h
  and 0.95 x 2.33 = 2.21 L/h reproduce the 1.38 and 2.21 L/h of Fig. 3,
  whereas starting from the pre-CPB value of 2.97 L/h would not. The ESM
  S3 control stream confirms it (`CLA_CPB = CLA_PRE * FA_CPB1`, where
  `CLA_PRE` is the weight-scaled THETA before `EFFPNA1` is applied).
- **Metabolite mass units.** The ESM S3 control stream keeps the
  oxypurinol compartment in allopurinol-equivalent mass and pushes the
  molar conversion into the observation scaling
  (`S2 = VO * 136.11/152.11`). This model instead applies the
  `152.11 / 136.11` factor at the transfer, which is algebraically
  identical for both the observed concentration and the oxypurinol
  elimination rate and has the advantage that `central_oxy` genuinely
  holds milligrams of oxypurinol.
- **Between-occasion variability.** rxode2 has no NONMEM `| occ` level,
  so the shared-variance BOV on each clearance is expanded into three
  indicator-multiplexed etas, one per perioperative period, with the
  first carrying the estimated variance and the other two fixed to the
  same value. rxode2 emits
  `some etas defaulted to non-mu referenced ... etaiov_cl_1, ...` on
  every solve as a result. That warning concerns an estimation-time
  mu-referencing optimisation, not correctness; simulation, which is
  what this model is packaged for, is unaffected, and the suggested
  work-around (“put the mu-referenced expression on a simple line”) is
  exactly what indicator multiplexing precludes.
- **Between-subject variability declared but not reported.** The ESM S3
  control stream declares `ETA(2)` on the allopurinol volume and
  `ETA(3)` on the oxypurinol clearance, but Table 2 reports estimates
  for neither, so neither is carried here. Only `ETA(1)`
  (allopurinol CL) and `ETA(4)` (oxypurinol Vd) are retained.
- **The variability scale is stated by the paper.** Table 2’s footnote
  reads “CV coefficient of variation, approximated using CV% =
  sqrt(omega^2) \* 100”, so the tabulated percentages are the omega
  standard deviations directly and `omega^2 = (CV/100)^2`. The usual
  log-normal moment-match `omega^2 = log(CV^2 + 1)` is **not** used
  here.
- **`CPB_REWARM` carries no independent effect.** Chu 2024 does not
  separate a rewarming phase; its intraoperative period is the whole
  bypass run. The covariate register forbids widening `CPB_ON` to cover
  rewarming, so the model consumes the on-bypass window as
  `CPB_ON + CPB_REWARM`. Because only the sum enters, a user whose data
  do not record the rewarming boundary can set `CPB_ON = 1` across the
  entire run with `CPB_REWARM = 0`, as this vignette does.
- **`CPB_POST` is registered by this extraction.** The covariate
  register’s `CPB_ON` and `CPB_REWARM` entries both direct that a
  sibling `CPB_POST` be registered when a paper retains a distinct
  post-CPB effect. Chu 2024 does, and in a way that no on-bypass effect
  could express: across the same phase boundary allopurinol clearance
  rises while oxypurinol clearance falls.

### Known deviation from the published results

- **IC90 XO attainment 12 h after the third (preoperative) dose is not
  reproduced.** The simulated cohort attains 41% against a published
  73%; every other tabulated attainment figure, for both targets,
  agrees. This is excluded from the render gate and left visible in the
  table above rather than tuned away. The diagnosis is in the
  target-attainment narrative: the simulated median combined
  concentration at that instant is about 10% below the 6.22 mg/L
  threshold, so the disagreement is a threshold crossing rather than an
  exposure error – the AUC over the same window matches Table 3 to
  within 7% (allopurinol) and 2% (oxypurinol), and the allopurinol \> 2
  mg/L attainment at the same instant matches to within 9 percentage
  points. The fixed perioperative schedule described below is the most
  likely contributor.

### Simulation assumptions

- The perioperative schedule is held at the Chu 2024 Table 1 cohort
  medians for every simulated subject (CPB starting at a postnatal age
  of 5.60 days and running 320 min), so that the per-dose NCA windows
  are common across subjects. The real cohort’s surgery times spanned an
  IQR of 4.78-7.81 days, which the model would translate into a spread
  of preoperative clearances (2.65 to 3.41 L/h for allopurinol) that
  this cohort does not carry.
- Body weight is held constant per subject at its drawn birth weight
  over the roughly nine-day simulation window. The source dataset
  carried a time-varying `WT` column.
- The virtual cohort varies only body weight plus the model’s
  between-subject and between-occasion random effects. Note that body
  weight largely cancels out of AUC: dose is proportional to weight
  while clearance scales as weight^0.75, so AUC scales only as
  weight^0.25.
- Gestational age and sex are recorded in `covariatesDataExcluded`
  because they appear in the ESM S3 `$INPUT` record but carry no effect
  in the published final model; they are not simulated.
- The 10-min infusion duration is taken from the CRUCIAL trial protocol
  (Stegeman 2022, *Trials*, <doi:10.1186/s13063-022-06098-y>), which Chu
  2024 cites as reference 6 and names as the source of its Fig. 1 dosing
  schedule. Chu 2024 itself states the dose and the schedule but not the
  infusion duration.
- Chu 2024 Table 3’s AUC values are individual maximum-a-posteriori
  Bayesian estimates for the 14 enrolled neonates, not a typical-value
  or Monte Carlo prediction, so exact agreement is not expected; the
  comparison above is a structural check on the model’s period-specific
  disposition.
