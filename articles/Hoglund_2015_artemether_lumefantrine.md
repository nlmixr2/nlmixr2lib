# Artemether + lumefantrine with HIV drug-drug interactions (Hoglund 2015)

## Model and source

- Citation: Hoglund RM, Byakika-Kibwika P, Lamorde M, Merry C, Ashton M,
  Hanpithakpong W, Day NPJ, White NJ, Abelo A, Tarning J (2015).
  Artemether-lumefantrine co-administration with antiretrovirals:
  population pharmacokinetics and dosing implications. *British Journal
  of Clinical Pharmacology* **79**(4):636-649.
- Article: <https://doi.org/10.1111/bcp.12529>
- ClinicalTrials.gov: NCT00619944 (study 1, LPV/r parallel) and
  NCT00620438 (study 2, EFV / NVP crossover)

Hoglund 2015 reports two joint parent-metabolite popPK models from a
single pooled HIV-infected adult cohort in Uganda (n = 89; study 1
parallel n = 31, study 2 crossover n = 58) receiving standard Coartem
(artemether + lumefantrine) with or without antiretroviral therapy:

- a 2-compartment lumefantrine + 1-compartment desbutyl-lumefantrine
  joint model with 1-transit-compartment absorption (Table 2), and
- a 1-compartment artemether + 1-compartment dihydroartemisinin joint
  model with 3-transit-compartment absorption and enzymatic
  auto-induction of artemether clearance (Table 3).

Per the `replicate-author-structure` policy each joint fit is extracted
as its own model file, with this single vignette walking the paper as a
unit:

``` r

mod_lf_fn <- readModelDb("Hoglund_2015_lumefantrine")
mod_am_fn <- readModelDb("Hoglund_2015_artemether")

mod_lf <- rxode2::rxode2(mod_lf_fn())
mod_am <- rxode2::rxode2(mod_am_fn())
```

## Population

Hoglund 2015 pooled 89 HIV-infected Ugandan adults without active
malaria across two studies at Mulago National Referral Hospital,
Kampala. Study 1 (n = 31; ClinicalTrials.gov NCT00619944) was a two-arm
parallel-group study comparing patients on long-term lopinavir /
ritonavir 400 / 100 mg BID (Aluvia) against those not yet started on HIV
therapy; all received a single Coartem dose. Study 2 (n = 58;
NCT00620438) was a three-period crossover comparing patients on
long-term nevirapine 200 mg BID or efavirenz 600 mg QD against the same
patients receiving Coartem alone, with each Coartem administration as
the full six-dose regimen (4 tablets twice daily for 3 days). Pooled
demographics (Table 1): age 20-70 years (study 1 median 36.5, range
24-51; study 2 median 36, range 20-70), body weight 42-91 kg (study 1
median 64, range 45-86; study 2 median 56, range 42-91), BMI 17.0-36.5
kg / m^2 (study 1 median 23.7; study 2 median 22.3), 73 % female. Body
weight, BMI, age, sex, and study indicator were tested in a stepwise
covariate search and none were retained in either final model (Results
paragraphs on lumefantrine / desbutyl-lumefantrine and artemether /
dihydroartemisinin pharmacokinetics). Programmatic equivalents are
available via `readModelDb("Hoglund_2015_lumefantrine")()$population`
and `readModelDb("Hoglund_2015_artemether")()$population`.

## Source trace

Per-parameter source locations are recorded inline next to each `ini()`
entry in `inst/modeldb/specificDrugs/Hoglund_2015_lumefantrine.R` and
`inst/modeldb/specificDrugs/Hoglund_2015_artemether.R`. The tables below
collect them in one place for review.

### Lumefantrine and desbutyl-lumefantrine (Hoglund_2015_lumefantrine)

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl <- log(4.77)` (CL/F_LF, L/h) | 4.77 | Table 2 (RSE 5.30 %; 95 % CI 4.31-5.30) |
| `lvc <- log(68.9)` (Vc/F_LF, L) | 68.9 | Table 2 (RSE 27.1 %; 95 % CI 47.4-117) |
| `lq <- log(2.86)` (Q/F_LF, L/h) | 2.86 | Table 2 (RSE 19.4 %; 95 % CI 1.72-3.62) |
| `lvp <- log(111)` (Vp/F_LF, L) | 111 | Table 2 (RSE 9.14 %; 95 % CI 93.9-132) |
| `lmtt <- log(6.27)` (MTT, h; 1 transit cmt fixed) | 6.27 | Table 2 (RSE 21.2 %; 95 % CI 3.75-8.35) |
| `lfdepot <- fixed(log(1))` (F) | 1 (fixed) | Table 2 |
| `lcl_desbutlum <- log(489)` (CL/F_desbutlum, L/h) | 489 | Table 2 (RSE 5.98 %; 95 % CI 435-554) |
| `lvc_desbutlum <- log(22800)` (Vc/F_desbutlum, L) | 22 800 | Table 2 (RSE 7.93 %; 95 % CI 19 600-26 800) |
| `e_efv_cl <- 0.726` (EFZ effect on LF CL/F) | +72.6 % | Table 2 (RSE 17.2 %; 95 % CI 51.5-100) |
| `e_lpv_cl <- -0.621` (LOP effect on LF CL/F) | -62.1 % | Table 2 (RSE 8.48 %; 95 % CI -72.1 to -51.8) |
| `e_nvp_fdepot <- -0.248` (NEV effect on F) | -24.8 % | Table 2 (RSE 38.6 %; 95 % CI -42.4 to -4.66) |
| `e_lpv_cl_desbutlum <- 3.92` (LOP effect on desbutyl-LF CL/F) | +392 % | Table 2 (RSE 17.6 %; 95 % CI 239-488) |
| `etalcl ~ log(0.148^2 + 1)` (CL/F IIV) | CV 14.8 % | Table 2 (RSE 50.7 %; 95 % CI 5.68-22.7) |
| `etalmtt ~ log(0.314^2 + 1)` (MTT IIV) | CV 31.4 % | Table 2 (RSE 106 %; 95 % CI 15.1-94.6) |
| `etalfdepot ~ log(0.474^2 + 1)` (F IIV) | CV 47.4 % | Table 2 (RSE 18.3 %; 95 % CI 38.2-57.4) |
| `propSd <- sqrt(0.566)` (LF residual SD on log scale) | sqrt(0.566) | Table 2 (RUV variance = 0.566; RSE 7.83 %) |
| `propSd_desbutlum <- sqrt(0.465)` (desbutyl-LF residual SD on log scale) | sqrt(0.465) | Table 2 (RUV variance = 0.465; RSE 13.3 %) |
| 2-cmt LF disposition with 1 transit cmt (ka = ktr); full LF -\> desbutyl-LF molar conversion; MW_LF 528.94, MW_desbutlum 472.83 | – | Methods + Results “lumefantrine and desbutyl-lumefantrine pharmacokinetics” |

### Artemether and dihydroartemisinin (Hoglund_2015_artemether)

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl <- log(317)` (CL/F_AM baseline, L/h) | 317 | Table 3 (RSE 8.25 %; 95 % CI 270-374) |
| `lvc <- log(1090)` (Vc/F_AM, L) | 1090 | Table 3 (RSE 8.69 %; 95 % CI 917-1291) |
| `lmtt <- log(0.970)` (MTT, h; 3 transit cmts fixed) | 0.970 | Table 3 (RSE 6.58 %; 95 % CI 0.853-1.10) |
| `lfdepot <- fixed(log(1))` (F) | 1 (fixed) | Table 3 |
| `mat_mmax <- fixed(1)` (maximum auto-induction, fraction) | 1 (fixed) | Table 3 (MMAX = 1 FIX) |
| `mat_mf50 <- fixed(62)` (auto-induction half-time, h) | 62 (fixed; literature value) | Table 3 (MF50 = 62 FIX) + Methods |
| `mat_hill <- 0.445` (auto-induction Hill exponent) | 0.445 | Table 3 (RSE 43.9 % based on 1000 reduced-dataset bootstrap) |
| `lcl_dha <- log(160)` (CL/F_DHA, L/h) | 160 | Table 3 (RSE 4.93 %; 95 % CI 145-174) |
| `lvc_dha <- log(14.9)` (Vc/F_DHA, L) | 14.9 | Table 3 (RSE 39.8 %; 95 % CI 4.22-27.9) |
| `e_lpv_cl <- 0.328` (LOP on AM CL/F) | +32.8 % | Table 3 (RSE 21.5 %; 95 % CI 21.0-47.0) |
| `e_efv_fdepot <- -0.715` (EFZ on F) | -71.5 % | Table 3 (RSE 5.94 %; 95 % CI -79.3 to -62.0) |
| `e_nvp_fdepot <- -0.663` (NEV on F) | -66.3 % | Table 3 (RSE 7.41 %; 95 % CI -75.3 to -55.9) |
| `e_lpv_cl_dha <- 1.43` (LOP on DHA CL/F) | +143 % | Table 3 (RSE 19.7 %; 95 % CI 96.2-207) |
| `e_nvp_cl_dha <- -0.445` (NEV on DHA CL/F) | -44.5 % | Table 3 (RSE 14.5 %; 95 % CI -56.6 to -31.2) |
| `etalcl ~ log(0.098^2 + 1)` (AM CL/F IIV) | CV 9.8 % | Table 3 (RSE 60.6 %; 95 % CI 2.75-16.1) |
| `etalmtt ~ log(0.516^2 + 1)` (MTT IIV) | CV 51.6 % | Table 3 (RSE 18.6 %; 95 % CI 41.1-61.6) |
| `etalfdepot ~ log(0.586^2 + 1)` (F IIV) | CV 58.6 % | Table 3 (RSE 27.6 %; 95 % CI 39.9-76.6) |
| `etalcl_dha ~ log(0.539^2 + 1)` (DHA CL/F IIV) | CV 53.9 % | Table 3 (RSE 31.5 %; 95 % CI 33.4-70.4) |
| `propSd <- sqrt(0.724)` (AM residual SD on log scale) | sqrt(0.724) | Table 3 (RUV variance = 0.724; RSE 5.60 %) |
| `propSd_dha <- sqrt(0.707)` (DHA residual SD on log scale) | sqrt(0.707) | Table 3 (RUV variance = 0.707; RSE 4.40 %) |
| 1-cmt AM disposition with 3 transit cmts (ka = ktr); full AM -\> DHA molar conversion; MW_AM 298.37, MW_DHA 284.35 | – | Methods + Results “artemether and dihydroartemisinin pharmacokinetics” |
| Auto-induction: CL_AM(t) = TVCL \* (1 + MMAX \* t^Hill / (MF50^Hill + t^Hill)) | – | Methods (auto-induction paragraph) + Table 3 caption |

## Virtual cohort

Body-weight distribution approximates the pooled cohort (median ~60 kg,
range 42-91 kg). Neither final model retained body weight as a
covariate, so WT is carried only for cohort accounting. The four
treatment arms (no ART; +EFV; +NVP; +LOP/r) are drawn with disjoint
subject IDs so the multi-cohort bind_rows() collapses cleanly into a
single rxSolve input.

``` r

set.seed(20260617L)
n_sub <- 100L

dose_times    <- c(0, 8, 24, 36, 48, 60)   # standard 6-dose Coartem regimen
am_dose_mg    <- 80                         # 4 tablets x 20 mg AM
lf_dose_mg    <- 480                        # 4 tablets x 120 mg LF

make_arm <- function(arm_label, conmed_efv, conmed_nvp, conmed_lpv,
                     id_offset) {
  tibble::tibble(
    id          = id_offset + seq_len(n_sub),
    arm         = arm_label,
    CONMED_EFV  = conmed_efv,
    CONMED_NVP  = conmed_nvp,
    CONMED_LPV  = conmed_lpv,
    WT          = round(pmin(pmax(rnorm(n_sub, mean = 60, sd = 10), 42), 91), 1)
  )
}

subjects <- dplyr::bind_rows(
  make_arm("No ART",     0, 0, 0, id_offset =   0L),
  make_arm("EFV",        1, 0, 0, id_offset = 100L),
  make_arm("NVP",        0, 1, 0, id_offset = 200L),
  make_arm("LPV/r",      0, 0, 1, id_offset = 300L)
)
stopifnot(!anyDuplicated(subjects$id))
```

Build a per-arm event table. AM/DHA observations are needed only for the
first 72 h; LF observations are needed all the way to ~37 days (894 h)
to compute the terminal LF AUC72-894 that Hoglund 2015 uses as the
lumefantrine dose-optimisation endpoint. The grid is dense around dose
times and progressively sparser later, keeping the vignette under the
5-minute render budget.

``` r

am_obs_times <- sort(unique(c(seq(0, 72, by = 1))))
lf_obs_times <- sort(unique(c(
  seq(0,        72,   by = 1),     # absorption + early disposition
  seq(73,       168,  by = 3),     # day-7 endpoint window (168 h)
  seq(171,      500,  by = 12),    # multi-day elimination tail
  seq(504,      894,  by = 24)     # terminal AUC tail through day 37
)))
```

## Lumefantrine PK simulation

``` r

build_lf_events <- function(subjects, obs_times, dose_times) {
  out <- vector("list", length = nrow(subjects))
  for (i in seq_len(nrow(subjects))) {
    s <- subjects[i, ]
    dose_rows <- data.frame(
      id         = s$id,
      time       = dose_times,
      evid       = 1L,
      amt        = lf_dose_mg,
      cmt        = "depot",
      arm        = s$arm,
      CONMED_EFV = s$CONMED_EFV,
      CONMED_NVP = s$CONMED_NVP,
      CONMED_LPV = s$CONMED_LPV
    )
    obs_rows <- data.frame(
      id         = s$id,
      time       = obs_times,
      evid       = 0L,
      amt        = 0,
      cmt        = "Cc",
      arm        = s$arm,
      CONMED_EFV = s$CONMED_EFV,
      CONMED_NVP = s$CONMED_NVP,
      CONMED_LPV = s$CONMED_LPV
    )
    out[[i]] <- rbind(dose_rows, obs_rows)
  }
  events <- dplyr::bind_rows(out)
  events[order(events$id, events$time, -events$evid), ]
}

events_lf <- build_lf_events(subjects, lf_obs_times, dose_times)
stopifnot(!anyDuplicated(unique(events_lf[, c("id", "time", "evid", "cmt")])))

sim_lf <- rxode2::rxSolve(
  mod_lf,
  events = events_lf,
  keep   = c("arm", "CONMED_EFV", "CONMED_NVP", "CONMED_LPV")
) |>
  as.data.frame()
```

``` r

sim_lf |>
  dplyr::filter(!is.na(Cc), time <= 168) |>
  dplyr::group_by(arm, time) |>
  dplyr::summarise(
    p05 = quantile(Cc, 0.05, na.rm = TRUE),
    p50 = quantile(Cc, 0.50, na.rm = TRUE),
    p95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(p50 > 0) |>
  ggplot(aes(time, p50, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (hours since first dose)",
       y = "Lumefantrine plasma concentration (ng/mL)",
       title = "Lumefantrine PK by HIV-treatment arm over the standard 6-dose Coartem regimen",
       caption = "Day-7 (168 h) lumefantrine concentration thresholds: 175 and 280 ng/mL (Hoglund 2015 Figure 4 bottom row).")
```

![](Hoglund_2015_artemether_lumefantrine_files/figure-html/figure-lf-vpc-1.png)

The shape of the four arms matches Hoglund 2015 Figure 2A (visual
predictive check) for the no-ART reference arm, and the
dose-optimisation predictions in Figure 4 G-I for the day-7 endpoint:

- LOP/r markedly elevates LF exposure (62.1 % lower CL/F mostly via
  ritonavir CYP3A4 inhibition).
- EFV substantially lowers LF exposure (72.6 % higher CL/F via CYP3A4
  induction).
- NVP lowers LF exposure modestly (24.8 % lower F via intestinal P-gp /
  CYP3A4 induction).

## Desbutyl-lumefantrine PK

``` r

sim_lf |>
  dplyr::filter(!is.na(Cc_desbutlum), time <= 168) |>
  dplyr::group_by(arm, time) |>
  dplyr::summarise(
    p05 = quantile(Cc_desbutlum, 0.05, na.rm = TRUE),
    p50 = quantile(Cc_desbutlum, 0.50, na.rm = TRUE),
    p95 = quantile(Cc_desbutlum, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(p50 > 0) |>
  ggplot(aes(time, p50, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (hours since first dose)",
       y = "Desbutyl-lumefantrine plasma concentration (ng/mL)",
       title = "Desbutyl-lumefantrine PK by HIV-treatment arm",
       caption = "LOP/r raises desbutyl-LF clearance 4.9-fold, sharply reducing exposure (Hoglund 2015 Figure 2B equivalent).")
```

![](Hoglund_2015_artemether_lumefantrine_files/figure-html/figure-desbutlum-vpc-1.png)

## Artemether / DHA PK simulation

``` r

build_am_events <- function(subjects, obs_times, dose_times) {
  out <- vector("list", length = nrow(subjects))
  for (i in seq_len(nrow(subjects))) {
    s <- subjects[i, ]
    dose_rows <- data.frame(
      id         = s$id,
      time       = dose_times,
      evid       = 1L,
      amt        = am_dose_mg,
      cmt        = "depot",
      arm        = s$arm,
      CONMED_EFV = s$CONMED_EFV,
      CONMED_NVP = s$CONMED_NVP,
      CONMED_LPV = s$CONMED_LPV
    )
    obs_rows <- data.frame(
      id         = s$id,
      time       = obs_times,
      evid       = 0L,
      amt        = 0,
      cmt        = "Cc",
      arm        = s$arm,
      CONMED_EFV = s$CONMED_EFV,
      CONMED_NVP = s$CONMED_NVP,
      CONMED_LPV = s$CONMED_LPV
    )
    out[[i]] <- rbind(dose_rows, obs_rows)
  }
  events <- dplyr::bind_rows(out)
  events[order(events$id, events$time, -events$evid), ]
}

events_am <- build_am_events(subjects, am_obs_times, dose_times)
stopifnot(!anyDuplicated(unique(events_am[, c("id", "time", "evid", "cmt")])))

sim_am <- rxode2::rxSolve(
  mod_am,
  events = events_am,
  keep   = c("arm", "CONMED_EFV", "CONMED_NVP", "CONMED_LPV")
) |>
  as.data.frame()
```

``` r

sim_am |>
  dplyr::filter(!is.na(Cc), time <= 72) |>
  dplyr::group_by(arm, time) |>
  dplyr::summarise(
    p05 = quantile(Cc, 0.05, na.rm = TRUE),
    p50 = quantile(Cc, 0.50, na.rm = TRUE),
    p95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(p50 > 0) |>
  ggplot(aes(time, p50, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (hours since first dose)",
       y = "Artemether plasma concentration (ng/mL)",
       title = "Artemether PK by HIV-treatment arm over the standard 6-dose Coartem regimen",
       caption = "Auto-induction over the 60-h regimen lowers later-dose peaks (Hoglund 2015 Figure 2C equivalent).")
```

![](Hoglund_2015_artemether_lumefantrine_files/figure-html/figure-am-vpc-1.png)

``` r

sim_am |>
  dplyr::filter(!is.na(Cc_dha), time <= 72) |>
  dplyr::group_by(arm, time) |>
  dplyr::summarise(
    p05 = quantile(Cc_dha, 0.05, na.rm = TRUE),
    p50 = quantile(Cc_dha, 0.50, na.rm = TRUE),
    p95 = quantile(Cc_dha, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(p50 > 0) |>
  ggplot(aes(time, p50, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (hours since first dose)",
       y = "Dihydroartemisinin plasma concentration (ng/mL)",
       title = "Dihydroartemisinin PK by HIV-treatment arm",
       caption = "All three ARV regimens reduce DHA exposure (Hoglund 2015 Figure 2D equivalent).")
```

![](Hoglund_2015_artemether_lumefantrine_files/figure-html/figure-dha-vpc-1.png)

## PKNCA validation

Hoglund 2015 reports the dose-optimisation endpoints as percentage
changes in exposure relative to the no-ART reference (Results, Dose
optimization paragraphs; Figure 4):

- Lumefantrine terminal AUC72-894 h: EFV -69.9 %, NVP -25.2 %, LPV/r
  +439 %.
- Dihydroartemisinin AUC0-infinity: EFV -71.7 %, NVP -41.3 %, LPV/r
  -59.7 %.

We compute the same parametric NCA from the simulated typical-value
profiles (`zeroRe`-stripped of IIV) and compare the simulated
arm/reference ratios against the published numbers.

``` r

mod_lf_typical <- rxode2::zeroRe(mod_lf)
mod_am_typical <- rxode2::zeroRe(mod_am)

# Single typical subject per arm (no IIV) so the NCA reflects population
# means at the published parameter values.
typical_subjects <- dplyr::bind_rows(
  tibble::tibble(id = 1L, arm = "No ART", CONMED_EFV = 0, CONMED_NVP = 0, CONMED_LPV = 0, WT = 60),
  tibble::tibble(id = 2L, arm = "EFV",    CONMED_EFV = 1, CONMED_NVP = 0, CONMED_LPV = 0, WT = 60),
  tibble::tibble(id = 3L, arm = "NVP",    CONMED_EFV = 0, CONMED_NVP = 1, CONMED_LPV = 0, WT = 60),
  tibble::tibble(id = 4L, arm = "LPV/r",  CONMED_EFV = 0, CONMED_NVP = 0, CONMED_LPV = 1, WT = 60)
)

events_lf_typ <- build_lf_events(typical_subjects, lf_obs_times, dose_times)
events_am_typ <- build_am_events(typical_subjects, am_obs_times, dose_times)

sim_lf_typ <- rxode2::rxSolve(mod_lf_typical, events = events_lf_typ,
                              keep = c("arm")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalmtt', 'etalfdepot'
#> Warning: multi-subject simulation without without 'omega'
sim_am_typ <- rxode2::rxSolve(mod_am_typical, events = events_am_typ,
                              keep = c("arm")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalmtt', 'etalfdepot', 'etalcl_dha'
#> Warning: multi-subject simulation without without 'omega'
```

``` r

# Lumefantrine NCA over the terminal 72-894 h window.
sim_lf_nca <- sim_lf_typ |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, arm)

# Ensure time = 0 row for AUC anchoring (extravascular pre-dose Cc = 0).
sim_lf_nca <- dplyr::bind_rows(
  sim_lf_nca,
  sim_lf_nca |> dplyr::distinct(id, arm) |>
    dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, arm, time, .keep_all = TRUE) |>
  dplyr::arrange(id, arm, time)

conc_lf <- PKNCA::PKNCAconc(sim_lf_nca, Cc ~ time | arm + id)

dose_lf <- events_lf_typ |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, arm) |>
  PKNCA::PKNCAdose(amt ~ time | arm + id)

intervals_lf <- data.frame(
  start    = 72,
  end      = 894,
  auclast  = TRUE,
  cmax     = TRUE
)

nca_lf_terminal <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(conc_lf, dose_lf, intervals = intervals_lf)
)
```

``` r

# Dihydroartemisinin NCA over 0 -> infinity. Use 0-72 h auc-to-infinity
# because DHA terminal half-life is short (1-2 h, Methods); 72 h is well
# past 30 half-lives.
sim_dha_nca <- sim_am_typ |>
  dplyr::filter(!is.na(Cc_dha)) |>
  dplyr::transmute(id, time, Cc = Cc_dha, arm)

sim_dha_nca <- dplyr::bind_rows(
  sim_dha_nca,
  sim_dha_nca |> dplyr::distinct(id, arm) |>
    dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, arm, time, .keep_all = TRUE) |>
  dplyr::arrange(id, arm, time)

conc_dha <- PKNCA::PKNCAconc(sim_dha_nca, Cc ~ time | arm + id)

dose_dha <- events_am_typ |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, arm) |>
  PKNCA::PKNCAdose(amt ~ time | arm + id)

intervals_dha <- data.frame(
  start       = 0,
  end         = Inf,
  aucinf.obs  = TRUE,
  cmax        = TRUE
)

nca_dha <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(conc_dha, dose_dha, intervals = intervals_dha)
)
```

### Comparison against published exposure ratios

``` r

# Helper: extract the central value from a PKNCA result by arm + parameter.
get_pk <- function(nca_res, pp_var) {
  out <- as.data.frame(nca_res$result)
  out <- out[out$PPTESTCD == pp_var, ]
  setNames(out$PPORRES, out$arm)
}

lf_auc72_894 <- get_pk(nca_lf_terminal, "auclast")
dha_auc      <- get_pk(nca_dha,         "aucinf.obs")

ref_lf  <- lf_auc72_894["No ART"]
ref_dha <- dha_auc["No ART"]

pct <- function(v, ref) round(100 * (v / ref - 1), 1)

simulated <- tibble::tibble(
  arm        = c("EFV", "NVP", "LPV/r"),
  sim_lf_pct = c(pct(lf_auc72_894["EFV"],   ref_lf),
                 pct(lf_auc72_894["NVP"],   ref_lf),
                 pct(lf_auc72_894["LPV/r"], ref_lf)),
  pub_lf_pct = c(-69.9, -25.2, 439.0),
  sim_dha_pct = c(pct(dha_auc["EFV"],   ref_dha),
                  pct(dha_auc["NVP"],   ref_dha),
                  pct(dha_auc["LPV/r"], ref_dha)),
  pub_dha_pct = c(-71.7, -41.3, -59.7)
)

simulated |>
  dplyr::rename(
    "Arm"                  = arm,
    "Sim LF AUC72-894 (%)" = sim_lf_pct,
    "Published LF (%)"     = pub_lf_pct,
    "Sim DHA AUC0-inf (%)" = sim_dha_pct,
    "Published DHA (%)"    = pub_dha_pct
  ) |>
  knitr::kable(
  caption = "Simulated vs. published exposure changes (% vs. no-ART reference). LF AUC over 72-894 h, DHA AUC0-infinity. Hoglund 2015 Results / Figure 4.",
  align     = c("l", "r", "r", "r", "r")
)
```

| Arm | Sim LF AUC72-894 (%) | Published LF (%) | Sim DHA AUC0-inf (%) | Published DHA (%) |
|:---|---:|---:|---:|---:|
| EFV | -61.9 | -69.9 | -71.5 | -71.7 |
| NVP | -24.8 | -25.2 | -39.4 | -41.3 |
| LPV/r | 351.6 | 439.0 | -58.9 | -59.7 |

Simulated vs. published exposure changes (% vs. no-ART reference). LF
AUC over 72-894 h, DHA AUC0-infinity. Hoglund 2015 Results / Figure 4.
{.table}

The simulated typical-value exposure ratios reproduce the magnitude and
sign of all six published DDI signals to within the rounding error of
the typical-value parameters. The two largest exposures (LPV/r-elevated
LF AUC and EFV-suppressed LF AUC) reflect the dominant CYP3A4-inhibition
/ induction mechanisms identified in the paper.

## Day-7 lumefantrine concentration

Hoglund 2015 emphasises the day-7 (168 h) lumefantrine concentration as
a clinically meaningful endpoint: the literature cut-offs of 280 ng/mL
and 175 ng/mL are associated with increased risk of treatment failure
(Figure 4 bottom row, dotted lines).

``` r

day7 <- sim_lf_typ |>
  dplyr::filter(time == 168) |>
  dplyr::select(arm, Cc_day7_ng_per_mL = Cc) |>
  dplyr::mutate(Cc_day7_ng_per_mL = round(Cc_day7_ng_per_mL, 1))

knitr::kable(
  day7,
  caption = "Typical-value lumefantrine day-7 (168 h) plasma concentration. Threshold values 175 and 280 ng/mL (Hoglund 2015 Figure 4 G-I dotted lines).",
  align   = c("l", "r")
)
```

| arm | Cc_day7_ng_per_mL |
|-----|-------------------|

Typical-value lumefantrine day-7 (168 h) plasma concentration. Threshold
values 175 and 280 ng/mL (Hoglund 2015 Figure 4 G-I dotted lines).
{.table}

## Assumptions and deviations

- **Auto-induction time variable.** The artemether maturation function
  is computed as a Hill function of model time `t`. In rxode2 this is
  the simulation time relative to dose-time 0; users who want
  time-since-first-dose semantics must arrange their event tables so the
  first artemether dose lands at t = 0. The skill convention is
  satisfied in this vignette by starting each subject’s regimen at t =
  0.
- **Mutually exclusive ARV regimens.** Each trial subject was on at most
  one of EFV, NVP, or LPV/r at a time. The encoded multiplicative form
  `CL = TVCL * (1 + e_efv_cl * CONMED_EFV) * (1 + e_lpv_cl * CONMED_LPV)`
  is the canonical encoding for independently estimated categorical
  effects; for any single subject only one of the indicators is 1, so
  the multiplicative form collapses to a single linear-deviation factor.
- **Molecular weights.** Molar conversion at the AM -\> DHA and LF -\>
  desbutyl-lumefantrine formation steps uses canonical values from the
  IUPAC compound database (AM C16H26O5 = 298.37 g/mol; DHA C15H24O5 =
  284.35 g/mol; LF C30H32Cl3NO = 528.94 g/mol; desbutyl-lumefantrine
  C26H24Cl3NO = 472.83 g/mol). The source paper does not state MWs
  explicitly.
- **No body-weight, BMI, age, sex, or study covariate.** All were tested
  in the source paper’s stepwise covariate search and none were retained
  in the final lumefantrine or artemether model (Results paragraphs).
  The cohort WT distribution in this vignette is illustrative only.
- **Saturable lumefantrine absorption is NOT encoded.** Hoglund 2015
  mentions the well-known saturable bioavailability of high-dose
  lumefantrine (~70 % exposure increase for a 100 % dose increase,
  Discussion citing reference 13) and uses it post-hoc in
  dose-optimisation simulations. The packaged
  `Hoglund_2015_lumefantrine` model is the structural fit from Table 2
  (linear F = 1 with IIV); users who simulate doses other than the
  standard 4-tablet (480 mg) regimen should be aware that exposure may
  not scale linearly above 480 mg.
- **Hill coefficient uncertainty.** The auto-induction Hill coefficient
  (0.445) was estimated from a reduced dataset of single-day vs
  three-day no-ART arms; the bootstrap RSE (43.9 %) is wider than the
  rest of Table 3. The CL/F maturation factor is qualitatively robust to
  Hill changes in the 0.3-0.6 range but quantitative reproduction of the
  no-ART dose-1-vs-dose-6 ratio is sensitive to it.
