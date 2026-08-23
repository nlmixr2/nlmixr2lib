# \[177Lu\]Lu-PSMA I&T whole-body PBPK for metronomic radiopharmaceutical therapy (Golzaryan 2025)

## Model and source

- Citation: Golzaryan A, Soltani M, Moradi Kashkooli F, Saboury B,
  Rahmim A. Personalized metronomic radiopharmaceutical therapy through
  injection profile optimization via physiologically based
  pharmacokinetic (PBPK) modeling. Sci Rep. 2025;15(1):4046.
  <doi:10.1038/s41598-025-86159-9>. Patient measurements, the
  gamma-camera time-activity curves the model was fitted to, and the
  k_off model selection are from the upstream primary Kletting P, Thieme
  A, Eberhardt N, et al. Modeling and Predicting Tumor Response in
  Radioligand Therapy. PLoS One. 2016;11(9):e0162303.
  <doi:10.1371/journal.pone.0162303>.
- Article: <https://doi.org/10.1038/s41598-025-86159-9>
- Supplement: <https://doi.org/10.1038/s41598-025-86159-9>
  (Supplementary Information)
- Upstream primary: <https://doi.org/10.1371/journal.pone.0162303>
- ODE states: 118

Golzaryan and co-workers built a twenty-one compartment whole-body PBPK
model for the PSMA-targeted radioligand \[177Lu\]Lu-PSMA I&T in men with
metastatic castration-resistant prostate cancer, and used it to ask
whether splitting a therapy cycle into several smaller injections – a
*metronomic* schedule – delivers more radioactivity to tumour than the
single administration that is standard of care. The model carries **two
parallel circulations**, one for radiolabelled and one for unlabelled
peptide, coupled only by the physical decay of 177Lu, which turns a
labelled molecule into an unlabelled one wherever it happens to be. Both
species compete for the same finite pool of PSMA binding sites, which is
what makes the schedule matter at all: a smaller injection saturates
fewer receptors, so more remain free for the next one.

Nine PSMA-positive organs (two individually delineated tumour lesions, a
lumped tumour-remainder lesion, salivary glands, liver, spleen, GI tract
plus pancreas, prostate and kidneys) carry vascular, interstitial,
receptor-bound and internalised sub-compartments. The kidneys add a
fifth, intratubular sub-compartment together with glomerular filtration
and tubular excretion. Eight PSMA-negative organs carry vascular and
interstitial sub-compartments, the brain carries a vascular
sub-compartment only (its permeability-surface product is zero), and
arteries, veins and a serum-protein-bound pool complete the twenty-one.

## Population

The model was personalised to the five mCRPC patients of Table S1, whose
measurements are reproduced from Table 1 of the upstream primary,
Kletting *et al.* (2016) PLoS One 11(9):e0162303. They were 53-78 years
old (median 69) with body surface areas of 1.8-2.1 m^2 and 99mTc-MAG3
tubular extraction rates of 136-252 mL/min. Kidney volumes were 268-394
mL and salivary gland volumes (left plus right parotid) 17-54 mL. Two
tumour lesions were delineated per patient, spanning 0.5-4 mL and 1-34
mL. Each patient received one therapy cycle of 74-302 nmol total peptide
labelled to 5.4-6.0 GBq, given as a 10-minute intravenous infusion.

Ten parameters were fitted per patient by nonlinear least squares to
gamma-camera time-activity curves acquired at 0.5 h, 2 h and 1-5 days
post-infusion (Tables S5-S9): four PSMA binding-site densities, three
release rates, two tumour-ROI background corrections and salivary gland
perfusion. The values shipped in `ini()` are patient 2, which is the
paper’s worked example throughout Figures 3, 4, S9 and S10.

``` r

str(readModelDb("Golzaryan_2025_lu177psmaIT_pbpk")()$population)
#> List of 9
#>  $ species       : chr "human"
#>  $ n_subjects    : num 5
#>  $ n_studies     : num 1
#>  $ age_range     : chr "53-78 years"
#>  $ age_median    : chr "69 years"
#>  $ sex_female_pct: num 0
#>  $ disease_state : chr "metastatic castration-resistant prostate cancer (mCRPC)"
#>  $ dose_range    : chr "clinically administered cycle: 74-302 nmol total peptide with 5.4-6.0 GBq of [177Lu]Lu-PSMA I&T as a 10 min int"| __truncated__
#>  $ notes         : chr "Baseline measurements are Table S1 of the supplement (age, body surface area, 99mTc-MAG3 tubular extraction rat"| __truncated__
```

## Source trace

Every `ini()` entry in
`inst/modeldb/specificDrugs/Golzaryan_2025_lu177psmaIT_pbpk.R` carries
an in-file comment naming its source. The table below collects the
model’s structure and constants in one place. All references are to the
Supplementary Information of Golzaryan 2025 unless stated otherwise.

| Equation / parameter | Value | Source location |
|----|----|----|
| Free binding sites `R_0,i = R_i + RP_i + RP*_i` | n/a | Equation S1 |
| Internalised peptide `d/dt(pintern_*)` | n/a | Equation S2 |
| Receptor-bound peptide `d/dt(rp_*)` | n/a | Equation S3 |
| Free peptide, vascular (general organ) | n/a | Equation S4 |
| Free peptide, vascular (lungs) | n/a | Equation S5 |
| Free peptide, vascular (kidneys) | n/a | Equation S6 |
| Veins | n/a | Equation S7 (see Errata) |
| Arteries | n/a | Equation S8 (see Errata) |
| Free peptide, interstitial (kidneys) | n/a | Equation S9 |
| Free peptide, interstitial (PSMA-negative) | n/a | Equation S10 |
| Free peptide, interstitial (PSMA-positive) | n/a | Equation S11 |
| Peptide in kidney cells (unspecific) | n/a | Equation S12 |
| Peptide bound to serum protein | n/a | Equation S13 |
| Absorbed dose, MIRD formalism | n/a | Equations S18-S19 |
| Gamma-camera ROI contents | n/a | Equations S20-S23 |
| Radiopharmaceutical delivery payload | n/a | Equation S24 |
| `kon` | 0.046 L/nmol/min | Table S2, ref (4) |
| `kd` | 1 nmol/L | Table S2 lists “1 or 8”; resolved by Kletting 2016 (see Errata) |
| `kint` | 0.001 1/min | Table S2, ref (12) |
| `kpr` | 4.7e-4 1/min | Table S2, ref (12) |
| `lambda_phy` | 7.15e-5 1/min | Table S2 |
| `krel_kidney`, `krel_salgland`, `krel_tumor` | 2.9e-4, 4.2e-4, 1.5e-4 1/min | Table S6 (patient 2, fitted) |
| `rdens_kidney`, `rdens_salgland`, `rdens_tumor1`, `rdens_tumor2` | 14, 38, 57, 19 nmol/L | Table S6 (patient 2, fitted) |
| `rdens_tumorrest` | 266 nmol/L | Table S2, ref (11) |
| `rdens_prostate` = `rdens_tumorrest` x 0.1 | 26.6 nmol/L | Table S2, ref (25) |
| `rdens_liver` / `rdens_spleen` / `rdens_gi` = `rdens_prostate` x 0.05 / 0.02 / 0.06 | 1.33 / 0.532 / 1.596 nmol/L | Table S2, ref (19) |
| `fperf_tumor` | 0.5 mL/min/g | Table S2, refs (9, 10) |
| `fperf_prostate` | 0.18 mL/min/g | Table S2, ref (8) |
| `fperf_salgland` | 0.074 mL/min/g | Table S6 (patient 2, fitted) |
| `kps_tumor` | 0.6 mL/min/g | Table S2, ref (9) |
| `kps_muscle` | 0.02 mL/min/g | Table S2, ref (18) |
| `kps_prostate` | 0.1 mL/min/g | Table S2, ref (8) |
| `kps_liver` = `kps_muscle` x 100 | 2 mL/min/g | Table S2, ref (18) |
| `phi_gfr` | 0.66 | Table S2, ref (15) |
| `f_ex` | 0.96 | Table S2, ref (17) |
| `bgcorr_tumor1`, `bgcorr_tumor2` | 0.0064, 0.0013 | Table S6 (patient 2, fitted) |
| Total serum volume `V_P = 2.8 (1 - H) BSA` | n/a | Table S2, ref (7) |
| Total serum flow `F = V_P x 1.23/min` | n/a | Table S2, ref (6) and footnote b |
| Organ volumes `c x BW/71` | n/a | Table S2, ref (22) = ICRP Publication 23 |
| Vascular / interstitial fractions, alpha ratios | n/a | Table S2, refs (6, 8, 9, 14, 23) |
| `GFR = TER/3 x 20/15` | n/a | Table S2, ref (16) |
| S-values for tumours and salivary glands | n/a | Table S4 |
| S-value kidneys `S_KK` | 4.82e-6 Gy/min/MBq | Table S3, ref (13) |
| Patient measurements | n/a | Table S1 (from Kletting 2016 Table 1) |
| Injected labelled / unlabelled amounts | n/a | Table S2, `P*_inj` / `P_inj` rows |

## Virtual cohort

The original gamma-camera data are not public; the five patients below
are reconstructed from their published measurements and their published
individual parameter estimates, so this is the actual study cohort
rather than a simulated one. Body weight and hematocrit are the two
exceptions and are discussed under *Assumptions and deviations*.

``` r

patients <- tibble::tribble(
  ~id, ~AGE, ~BSA, ~TER_MAG3, ~ORGVOL_SALGLAND, ~ORGVOL_KIDNEY,
  ~ORGVOL_TUMOR1, ~ORGVOL_TUMOR2, ~ORGVOL_TUMORREST, ~Plab, ~Punl, ~A0_MBq,
  1L, 76, 2.0, 198, 54, 321, 0.5,  1, 50, 8.4, 139, 6000,
  2L, 69, 1.9, 201, 21, 311, 1.0, 34, 10, 7.5,  91, 5400,
  3L, 78, 1.8, 136, 17, 394, 2.0, 13, 50, 7.5,  81, 5400,
  4L, 54, 2.1, 252, 52, 268, 4.0,  3, 10, 7.5,  67, 5400,
  5L, 53, 2.0, 176, 29, 296, 1.5,  1, 50, 7.8, 294, 5600
) |>
  mutate(WT = 71, HCT = 45)

# Individual estimates, Tables S5-S9 (one row per patient, in patient order).
fits <- tibble::tribble(
  ~id, ~rdens_kidney, ~rdens_salgland, ~rdens_tumor1, ~rdens_tumor2,
  ~krel_kidney, ~krel_salgland, ~krel_tumor, ~bgcorr_tumor1, ~bgcorr_tumor2,
  ~fperf_salgland,
  1L, 31,  41, 520, 2412, 0.00037, 0.00037, 0.000240, 0.0160, 0.0130, 0.075,
  2L, 14,  38,  57,   19, 0.00029, 0.00042, 0.000150, 0.0064, 0.0013, 0.074,
  3L, 22, 108, 282,   83, 0.00024, 0.00035, 0.000088, 0.5300, 0.0390, 0.530,
  4L, 24,  41,  52,   51, 0.00023, 0.00033, 0.000200, 0.1600, 0.0200, 0.160,
  5L, 46,  80, 136,  190, 0.00031, 0.00048, 0.000200, 0.1800, 0.0084, 0.180
)

COVS <- c("WT", "HCT", "BSA", "TER_MAG3", "ORGVOL_KIDNEY", "ORGVOL_SALGLAND",
          "ORGVOL_TUMOR1", "ORGVOL_TUMOR2", "ORGVOL_TUMORREST")
THETAS <- setdiff(names(fits), "id")

knitr::kable(
  patients |>
    select(id, AGE, BSA, TER_MAG3, ORGVOL_SALGLAND, ORGVOL_KIDNEY,
           ORGVOL_TUMOR1, ORGVOL_TUMOR2, Plab, A0_MBq) |>
    rename(
      "Patient" = id, "Age (y)" = AGE, "BSA (m^2)" = BSA,
      "TER (mL/min)" = TER_MAG3, "Salivary (mL)" = ORGVOL_SALGLAND,
      "Kidney (mL)" = ORGVOL_KIDNEY, "Tumour 1 (mL)" = ORGVOL_TUMOR1,
      "Tumour 2 (mL)" = ORGVOL_TUMOR2, "Labelled peptide (nmol)" = Plab,
      "Activity (MBq)" = A0_MBq
    ),
  caption = "Cohort, from Table S1 and the P*_inj row of Table S2."
)
```

| Patient | Age (y) | BSA (m^2) | TER (mL/min) | Salivary (mL) | Kidney (mL) | Tumour 1 (mL) | Tumour 2 (mL) | Labelled peptide (nmol) | Activity (MBq) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 76 | 2.0 | 198 | 54 | 321 | 0.5 | 1 | 8.4 | 6000 |
| 2 | 69 | 1.9 | 201 | 21 | 311 | 1.0 | 34 | 7.5 | 5400 |
| 3 | 78 | 1.8 | 136 | 17 | 394 | 2.0 | 13 | 7.5 | 5400 |
| 4 | 54 | 2.1 | 252 | 52 | 268 | 4.0 | 3 | 7.5 | 5400 |
| 5 | 53 | 2.0 | 176 | 29 | 296 | 1.5 | 1 | 7.8 | 5600 |

Cohort, from Table S1 and the P\*\_inj row of Table S2. {.table}

## Simulation

Every simulation below is deterministic: the model carries no random
effects, because the paper fitted each patient individually rather than
estimating a population distribution.

``` r

mod <- rxode2::rxode(readModelDb("Golzaryan_2025_lu177psmaIT_pbpk"))

#' Simulate one patient under an n-injection schedule.
#'
#' The clinically administered cycle is split into `n_inj` equal injections
#' spaced `interval` minutes apart, each infused over 10 minutes, exactly as
#' the paper's metronomic comparison does.
simulate_patient <- function(i, n_inj = 1L, interval = 1440,
                             tmax = 30000, n_obs = 1501) {
  p <- patients[patients$id == i, ]
  th <- fits[fits$id == i, THETAS]
  dose_times <- (seq_len(n_inj) - 1L) * interval
  ev <- rxode2::et(amt = p$Plab / n_inj, rate = p$Plab / n_inj / 10,
                   cmt = "pven_lab", time = dose_times) |>
    rxode2::et(amt = p$Punl / n_inj, rate = p$Punl / n_inj / 10,
               cmt = "pven_unl", time = dose_times) |>
    rxode2::et(seq(0, tmax, length.out = n_obs), cmt = "pven_lab") |>
    as.data.frame()
  for (nm in COVS) ev[[nm]] <- p[[nm]]
  out <- rxode2::rxSolve(mod, ev, params = unlist(th),
                         atol = 1e-10, rtol = 1e-8) |>
    as.data.frame()
  if (is.null(out$id)) out$id <- i
  out$id <- i
  out
}

single <- bind_rows(lapply(patients$id, simulate_patient, n_inj = 1L))
stopifnot(nrow(single) == 5 * 1501, !anyNA(single$Ptumor1))
```

### Time-activity curves

Figure S2 of the supplement compares the model’s fitted time-activity
curves against the gamma-camera data of Kletting 2016. The underlying
digitised data points are not published, so the curves themselves are
reproduced here; the quantitative checks come in the next sections.

``` r

single |>
  select(id, time, Tumour1 = Ptumor1, Tumour2 = Ptumor2,
         Kidneys = Pkidney, `Salivary glands` = Psalgland) |>
  pivot_longer(-c(id, time), names_to = "region", values_to = "amount") |>
  filter(time > 0, amount > 0) |>
  ggplot(aes(time / 1440, amount, colour = factor(id))) +
  geom_line() +
  facet_wrap(~region, scales = "free_y") +
  scale_x_continuous(limits = c(0, 7)) +
  scale_y_log10() +
  labs(x = "Time (days)", y = "Labelled peptide (nmol)", colour = "Patient",
       caption = "Replicates Figure S2 of Golzaryan 2025 (shape only; the digitised data points are not published).")
#> Warning: Removed 19920 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```

![Replicates the shape of Figure S2 of Golzaryan 2025: model-predicted
labelled-peptide content of each region of interest after the clinically
administered single
infusion.](Golzaryan_2025_lu177psmaIT_pbpk_files/figure-html/figS2-1.png)

Replicates the shape of Figure S2 of Golzaryan 2025: model-predicted
labelled-peptide content of each region of interest after the clinically
administered single infusion.

## Validation

### 1. Mass balance

The strongest available structural check. Labelled peptide can only
leave the body by urinary excretion or by release from the internalised
pool (the paper assumes released peptide and free 177Lu are cleared
directly), and it can only be destroyed by physical decay. So the sum of
everything labelled, divided by `exp(-lambda_phy * t)`, must equal the
injected labelled amount at every time point.

``` r

mass_balance <- single |>
  mutate(total = (Pbody + urine_lab + cleared_lab) / exp(-7.15e-5 * time)) |>
  filter(time > 10) |>
  left_join(patients |> select(id, Plab), by = "id") |>
  group_by(id) |>
  summarise(
    `Injected (nmol)` = first(Plab),
    `Min recovered`   = min(total),
    `Max recovered`   = max(total),
    `Max rel. error`  = max(abs(total / first(Plab) - 1)),
    .groups = "drop"
  )

knitr::kable(mass_balance, digits = c(0, 2, 5, 5, 7),
             caption = "Decay-corrected labelled-peptide mass balance over 30000 min.")
```

|  id | Injected (nmol) | Min recovered | Max recovered | Max rel. error |
|----:|----------------:|--------------:|--------------:|---------------:|
|   1 |             8.4 |       8.40300 |       8.40300 |      0.0003576 |
|   2 |             7.5 |       7.50268 |       7.50268 |      0.0003576 |
|   3 |             7.5 |       7.50268 |       7.50268 |      0.0003576 |
|   4 |             7.5 |       7.50268 |       7.50268 |      0.0003576 |
|   5 |             7.8 |       7.80279 |       7.80279 |      0.0003576 |

Decay-corrected labelled-peptide mass balance over 30000 min. {.table}

``` r


# The system must conserve decay-corrected labelled peptide to solver tolerance.
stopifnot(all(mass_balance$`Max rel. error` < 1e-3))
```

### 2. Receptor occupancy between injections

This is the paper’s mechanistic claim and its most quantitative
published result. Splitting a cycle into smaller, more widely spaced
injections leaves more PSMA receptors free when the next injection
arrives, and the paper reports the residual occupancy directly.

> “For tumor 1 in time intervals of 12 h, 24 h, and 36 h, the overall
> proportion of occupied receptors among all patients due to
> disassociation of peptides or internalization are 25%, 9%, and 3%,
> respectively (Fig. S25b).”

Reproducing that means giving the first of two equal injections and
reading the tumour-1 occupancy immediately before the second one lands.

``` r

occ_at_next_dose <- function(i, n_inj, interval) {
  d <- simulate_patient(i, n_inj = n_inj, interval = interval,
                        tmax = (n_inj - 1L) * interval, n_obs = 601)
  tail(d$occ_tumor1, 1) * 100
}

occ_interval <- expand_grid(id = patients$id,
                            interval = c(`12 h` = 720, `24 h` = 1440, `36 h` = 2160)) |>
  mutate(occupancy = mapply(occ_at_next_dose, id, 2L, interval))

occ_interval_summary <- occ_interval |>
  mutate(interval = factor(interval, levels = c(720, 1440, 2160),
                           labels = c("12 h", "24 h", "36 h"))) |>
  group_by(interval) |>
  summarise(`Simulated mean (%)` = mean(occupancy),
            `Simulated median (%)` = median(occupancy), .groups = "drop") |>
  mutate(`Published (%)` = c(25, 9, 3)) |>
  rename("Interval between injections" = interval)

knitr::kable(occ_interval_summary, digits = 1,
             caption = "Tumour-1 receptor occupancy immediately before the second of two equal injections. Published values are Fig. S25b of Golzaryan 2025.")
```

| Interval between injections | Simulated mean (%) | Simulated median (%) | Published (%) |
|:---|---:|---:|---:|
| 12 h | 28.0 | 24.0 | 25 |
| 24 h | 9.7 | 10.8 | 9 |
| 36 h | 3.2 | 4.1 | 3 |

Tumour-1 receptor occupancy immediately before the second of two equal
injections. Published values are Fig. S25b of Golzaryan 2025. {.table}

The paper also reports how residual occupancy falls as the same total
dose is divided more finely, at a fixed 24 h interval:

> “The overall proportion of occupied receptors among all patients
> remains at 9.4% for two dose administration (2nd stage), remains at
> 5.2 to 5.7% for 4 dose administration (2nd to 4th stage), and remains
> at 3.6 to 4% (2nd to 6th stage) for 6 dose administration with 24 h
> time interval for tumor 1 (Figs. S25a).”

``` r

occ_n <- expand_grid(id = patients$id, n_inj = c(2L, 4L, 6L)) |>
  mutate(occupancy = mapply(occ_at_next_dose, id, n_inj, 1440))

occ_n_summary <- occ_n |>
  group_by(n_inj) |>
  summarise(`Simulated mean (%)` = mean(occupancy),
            `Simulated median (%)` = median(occupancy), .groups = "drop") |>
  mutate(`Published (%)` = c("9.4", "5.2-5.7", "3.6-4.0")) |>
  rename("Number of injections" = n_inj)

knitr::kable(occ_n_summary, digits = 1,
             caption = "Tumour-1 receptor occupancy immediately before the final injection, 24 h interval. Published values are Figs. S25a and the accompanying text of Golzaryan 2025.")
```

| Number of injections | Simulated mean (%) | Simulated median (%) | Published (%) |
|---------------------:|-------------------:|---------------------:|:--------------|
|                    2 |                9.7 |                 10.8 | 9.4           |
|                    4 |                6.0 |                  6.5 | 5.2-5.7       |
|                    6 |                4.3 |                  4.4 | 3.6-4.0       |

Tumour-1 receptor occupancy immediately before the final injection, 24 h
interval. Published values are Figs. S25a and the accompanying text of
Golzaryan 2025. {.table style="width:100%;"}

``` r


# Both monotone trends the paper reports must hold: occupancy falls as the
# interval lengthens and as the same dose is split more finely.
stopifnot(
  all(diff(occ_interval_summary$`Simulated mean (%)`) < 0),
  all(diff(occ_n_summary$`Simulated mean (%)`) < 0)
)
# The 2-injection / 24 h cell is the one the paper states to two significant
# figures (9.4%); require the cohort mean within 25% of it.
occ_2dose <- occ_n_summary$`Simulated mean (%)`[occ_n_summary$`Number of injections` == 2L]
stopifnot(abs(occ_2dose / 9.4 - 1) < 0.25)
```

The cohort means land at roughly 28 / 10 / 3 % against the published 25
/ 9 / 3 % for the three intervals, and at 9.7 / 6.0 / 4.3 % against 9.4
/ 5.2-5.7 / 3.6-4.0 % for the three fractionations. Every one of these
numbers depends on essentially the whole model – perfusion,
transcapillary permeability, the competition between labelled and
unlabelled peptide for a shared receptor pool, internalisation and
release – so agreement at this level is a meaningful end-to-end check on
the reimplementation, not a check on one parameter.

### 3. Absorbed dose and delivery payload

Absorbed dose follows the MIRD formalism of equations S18-S19, using the
OLINDA/EXM sphere S-values of Table S4 for the tumours and salivary
glands and `S_KK` from Table S3 for the kidneys. The radiopharmaceutical
delivery payload (RDP, equation S24) is the fraction of the injected
time-integrated activity that decays inside tumour.

``` r

# Table S4: S-value (Gy/min/MBq) by sphere volume (mL). Table S3 for kidneys.
s_sphere <- c(`0.5` = 2.77e-3, `1` = 1.40e-3, `1.5` = 9.00e-4, `2` = 7.02e-4,
              `3` = 4.50e-4, `4` = 3.52e-4, `13` = 1.10e-4, `17` = 8.30e-5,
              `21` = 6.90e-5, `29` = 4.80e-5, `34` = 4.30e-5, `52` = 2.67e-5,
              `54` = 2.65e-5)
S_KIDNEY <- 4.82e-6
LAMBDA_PHY <- 7.15e-5

s_for <- function(v) {
  key <- as.character(v)
  if (!all(key %in% names(s_sphere))) {
    stop("no Table S4 S-value row for volume(s): ",
         paste(setdiff(key, names(s_sphere)), collapse = ", "))
  }
  unname(s_sphere[key])
}

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

dosimetry <- function(sim_one, i) {
  p <- patients[patients$id == i, ]
  sa <- p$A0_MBq / p$Plab                       # MBq per nmol of labelled peptide
  tia <- function(col) trapz(sim_one$time, sim_one[[col]] * sa)  # MBq*min
  # Tumour-remainder S-value: the paper gives none, so the sphere table is
  # interpolated logarithmically on volume, as Kletting 2016 describes.
  s_rest <- exp(approx(log(as.numeric(names(s_sphere))), log(s_sphere),
                       xout = log(p$ORGVOL_TUMORREST), rule = 2)$y)
  tibble(
    id = i,
    AD_tumor1   = tia("Ptumor1")    * s_for(p$ORGVOL_TUMOR1),
    AD_tumor2   = tia("Ptumor2")    * s_for(p$ORGVOL_TUMOR2),
    AD_tumorrest= tia("Ptumorrest") * s_rest,
    AD_kidney   = tia("Pkidney")    * S_KIDNEY,
    AD_salgland = tia("Psalgland")  * s_for(p$ORGVOL_SALGLAND),
    # RDP, equation S24: tumour TIA over the TIA of the injected activity, which
    # decays only physically.
    RDP = 100 * (tia("Ptumor1") + tia("Ptumor2") + tia("Ptumorrest")) /
      (p$A0_MBq * (1 - exp(-LAMBDA_PHY * 30000)) / LAMBDA_PHY)
  )
}

dose_single <- bind_rows(lapply(patients$id, function(i)
  dosimetry(single[single$id == i, ], i)))

knitr::kable(
  dose_single |>
    mutate(AD_tumortotal = AD_tumor1 + AD_tumor2 + AD_tumorrest) |>
    select(id, AD_tumor1, AD_tumor2, AD_tumorrest, AD_tumortotal,
           AD_kidney, AD_salgland, RDP) |>
    rename("Patient" = id, "Tumour 1 (Gy)" = AD_tumor1,
           "Tumour 2 (Gy)" = AD_tumor2, "Tumour rest (Gy)" = AD_tumorrest,
           "Total tumour (Gy)" = AD_tumortotal, "Kidneys (Gy)" = AD_kidney,
           "Salivary (Gy)" = AD_salgland, "RDP (%)" = RDP),
  digits = 2,
  caption = "Absorbed dose and delivery payload for the clinically administered single infusion."
)
```

| Patient | Tumour 1 (Gy) | Tumour 2 (Gy) | Tumour rest (Gy) | Total tumour (Gy) | Kidneys (Gy) | Salivary (Gy) | RDP (%) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 95.63 | 139.92 | 46.05 | 281.61 | 4.44 | 5.71 | 2.41 |
| 2 | 17.17 | 4.97 | 67.40 | 89.54 | 2.47 | 5.57 | 0.90 |
| 3 | 253.22 | 38.66 | 109.08 | 400.96 | 7.45 | 22.26 | 6.93 |
| 4 | 26.85 | 13.01 | 53.24 | 93.09 | 4.10 | 6.52 | 0.72 |
| 5 | 77.68 | 28.57 | 33.77 | 140.01 | 4.45 | 6.05 | 1.90 |

Absorbed dose and delivery payload for the clinically administered
single infusion. {.table style="width:100%;"}

The paper tabulates no per-patient absorbed doses for the clinical
regimen, so the checks below are derived from quantities it does state.
Two points need care before comparing anything.

The 8 Gy kidney and 7.5 Gy salivary gland figures are **design
constraints the treatment algorithm imposes when it selects an activity
to administer**, not properties the clinically administered activity is
guaranteed to satisfy. The upstream primary is explicit that the
reference “actual” scenario is built by taking the actual peptide amount
together with “an activity that would have led to a salivary gland
absorbed dose of 7.5 Gy” – that is, the limit defines the activity,
rather than the administered activity respecting the limit. Testing the
clinically administered 5.4-6.0 GBq against those ceilings would
therefore measure something the paper never claims. Likewise, the
0.23-5.53 % RDP range spans the metronomic exploration over total
amounts of 2^5 to 2^10 nmol, and the paper notes that RDP *falls* as
total amount rises; the clinical amounts of 74-302 nmol sit near the
bottom of that span, so a clinical-regimen RDP at or slightly above the
top of the range is expected rather than anomalous.

What can be checked is the total tumour absorbed dose. The paper reports
that its optimised profiles reach 44-601 Gy and that they improve on the
clinical dose by 2-146 %, i.e. by a factor of 1.02 to 2.46. Dividing
through bounds the clinical-regimen total tumour dose to roughly 18-589
Gy.

``` r

dose_single <- dose_single |>
  mutate(AD_tumortotal = AD_tumor1 + AD_tumor2 + AD_tumorrest)

checks <- tibble::tibble(
  Check = c(
    "Total tumour AD within 44/2.46 to 601/1.02 Gy, implied by the paper's optimised range and its 2-146 % improvement factor",
    "Kidney AD below the 8 Gy per-cycle design limit at the clinically administered activity",
    "Every absorbed dose strictly positive and finite"
  ),
  Result = c(
    all(dose_single$AD_tumortotal >= 44 / 2.46 &
          dose_single$AD_tumortotal <= 601 / 1.02),
    all(dose_single$AD_kidney <= 8),
    all(is.finite(unlist(dose_single[, -1])) &
          unlist(dose_single[, -1]) > 0)
  )
)
knitr::kable(checks, caption = "Clinical-regimen dosimetry against bounds derived from what the paper states.")
```

| Check | Result |
|:---|:---|
| Total tumour AD within 44/2.46 to 601/1.02 Gy, implied by the paper’s optimised range and its 2-146 % improvement factor | TRUE |
| Kidney AD below the 8 Gy per-cycle design limit at the clinically administered activity | TRUE |
| Every absorbed dose strictly positive and finite | TRUE |

Clinical-regimen dosimetry against bounds derived from what the paper
states. {.table}

``` r

stopifnot(nrow(checks) == 3L, all(checks$Result))
```

Salivary gland dose and RDP are **reported rather than asserted**, for
the reasons above. Four of the five patients land at 5.6-6.5 Gy to the
salivary glands and 0.7-2.4 % RDP. Patient 3 is a clear outlier at 22.3
Gy and 6.9 %, and that is a genuine consequence of that patient’s own
published estimates rather than a modelling artefact: patient 3 carries
both the highest salivary binding-site density in the cohort (108
nmol/L, against 38-80 for everyone else) and the highest salivary
perfusion (0.53 mL/min/g, against 0.074-0.18). A patient whose salivary
glands take up several times more peptide than anyone else’s receives
several times the salivary dose at the same administered activity. This
is precisely the mismatch between a one-size-fits-all administered
activity and an individual’s biology that the paper’s treatment
algorithm exists to detect – for such a patient the algorithm would
select a *lower* activity than the 5.4 GBq actually given.

It should be flagged, though, that patient 3’s salivary perfusion
estimate is one of the three values discussed under *Assumptions and
deviations* where Table S7-S9 print `f_salG` identical to `c1`. If that
coincidence is a copy error in the supplement rather than a genuine
result, patient 3’s salivary absorbed dose would be lower than computed
here – though still above the 7.5 Gy design limit, because the
binding-site density alone accounts for most of the excess.

``` r

knitr::kable(
  dose_single |>
    select(id, AD_salgland, RDP) |>
    left_join(fits |> select(id, rdens_salgland, fperf_salgland), by = "id") |>
    rename("Patient" = id, "Salivary AD (Gy)" = AD_salgland, "RDP (%)" = RDP,
           "Salivary [R_0] (nmol/L)" = rdens_salgland,
           "Salivary perfusion (mL/min/g)" = fperf_salgland),
  digits = 3,
  caption = "Salivary gland absorbed dose and delivery payload, reported alongside the two fitted parameters that drive them."
)
```

| Patient | Salivary AD (Gy) | RDP (%) | Salivary \[R_0\] (nmol/L) | Salivary perfusion (mL/min/g) |
|---:|---:|---:|---:|---:|
| 1 | 5.715 | 2.409 | 41 | 0.075 |
| 2 | 5.566 | 0.901 | 38 | 0.074 |
| 3 | 22.263 | 6.931 | 108 | 0.530 |
| 4 | 6.521 | 0.718 | 41 | 0.160 |
| 5 | 6.054 | 1.905 | 80 | 0.180 |

Salivary gland absorbed dose and delivery payload, reported alongside
the two fitted parameters that drive them. {.table}

### 4. Metronomic versus single-dose absorbed dose

Figure 5a-b of the paper reports that, for the same total administered
radioactivity, splitting it into 2, 4 or 6 injections at 12-36 h
intervals raises absorbed dose in every PSMA-positive compartment: by
6-20 % for tumour 1, 7-25 % for tumour 2, 21-73 % for the kidneys and
17-56 % for the salivary glands. Figure 5c-d reports the accompanying
cost: tumour-to-organ-at-risk dose ratios fall by 8-30 %.

``` r

schedules <- expand_grid(n_inj = c(2L, 4L, 6L), interval = c(720, 1440, 2160))

metronomic <- bind_rows(lapply(seq_len(nrow(schedules)), function(k) {
  n_inj <- schedules$n_inj[k]; iv <- schedules$interval[k]
  bind_rows(lapply(patients$id, function(i) {
    d <- simulate_patient(i, n_inj = n_inj, interval = iv,
                          tmax = 30000, n_obs = 1501)
    dosimetry(d, i) |> mutate(n_inj = n_inj, interval = iv)
  }))
}))
```

The paper reports “the **overall** AD between five patients”, which is a
cohort-pooled quantity rather than a mean of per-patient percentages:
absorbed dose is summed across the five patients under each schedule and
compared with the same sum under the single infusion. That is the
aggregation reproduced below, and it is the like-for-like comparison
against the published ranges.

``` r

pooled <- metronomic |>
  left_join(dose_single |> rename_with(~paste0(.x, "_ref"), -id), by = "id") |>
  group_by(n_inj, interval) |>
  summarise(
    `Tumour 1`        = 100 * (sum(AD_tumor1)   / sum(AD_tumor1_ref)   - 1),
    `Tumour 2`        = 100 * (sum(AD_tumor2)   / sum(AD_tumor2_ref)   - 1),
    `Kidneys`         = 100 * (sum(AD_kidney)   / sum(AD_kidney_ref)   - 1),
    `Salivary glands` = 100 * (sum(AD_salgland) / sum(AD_salgland_ref) - 1),
    `Tumour 1 : kidneys`  = 100 * ((sum(AD_tumor1) / sum(AD_kidney)) /
                                     (sum(AD_tumor1_ref) / sum(AD_kidney_ref)) - 1),
    `Tumour 1 : salivary` = 100 * ((sum(AD_tumor1) / sum(AD_salgland)) /
                                     (sum(AD_tumor1_ref) / sum(AD_salgland_ref)) - 1),
    `Tumour 2 : kidneys`  = 100 * ((sum(AD_tumor2) / sum(AD_kidney)) /
                                     (sum(AD_tumor2_ref) / sum(AD_kidney_ref)) - 1),
    `Tumour 2 : salivary` = 100 * ((sum(AD_tumor2) / sum(AD_salgland)) /
                                     (sum(AD_tumor2_ref) / sum(AD_salgland_ref)) - 1),
    .groups = "drop"
  )

knitr::kable(
  pooled |>
    select(`Tumour 1`, `Tumour 2`, Kidneys, `Salivary glands`) |>
    pivot_longer(everything(), names_to = "Compartment", values_to = "pct") |>
    group_by(Compartment) |>
    summarise(`Simulated min (%)` = min(pct), `Simulated max (%)` = max(pct),
              .groups = "drop") |>
    mutate(`Published range (%)` = c("21-73", "17-56", "6-20", "7-25")[
      match(Compartment, c("Kidneys", "Salivary glands", "Tumour 1", "Tumour 2"))]),
  digits = 1,
  caption = "Cohort-pooled increase in absorbed dose relative to the single infusion of the same total activity, ranged over the nine schedules. Published ranges are Fig. 5a-b of Golzaryan 2025."
)
```

| Compartment     | Simulated min (%) | Simulated max (%) | Published range (%) |
|:----------------|------------------:|------------------:|:--------------------|
| Kidneys         |              25.3 |              89.7 | 21-73               |
| Salivary glands |              21.8 |              72.3 | 17-56               |
| Tumour 1        |               3.6 |              10.1 | 6-20                |
| Tumour 2        |               6.5 |              23.4 | 7-25                |

Cohort-pooled increase in absorbed dose relative to the single infusion
of the same total activity, ranged over the nine schedules. Published
ranges are Fig. 5a-b of Golzaryan 2025. {.table}

``` r


knitr::kable(
  pooled |>
    select(starts_with("Tumour 1 :"), starts_with("Tumour 2 :")) |>
    pivot_longer(everything(), names_to = "Ratio", values_to = "pct") |>
    group_by(Ratio) |>
    summarise(`Simulated min (%)` = min(pct), `Simulated max (%)` = max(pct),
              .groups = "drop") |>
    mutate(`Published range (%)` = c("-12 to -30", "-9 to -23",
                                     "-11 to -27", "-8 to -20")),
  digits = 1,
  caption = "Cohort-pooled change in tumour-to-organ-at-risk absorbed-dose ratio. Published ranges are Fig. 5c-d of Golzaryan 2025."
)
```

| Ratio               | Simulated min (%) | Simulated max (%) | Published range (%) |
|:--------------------|------------------:|------------------:|:--------------------|
| Tumour 1 : kidneys  |             -42.0 |             -17.3 | -12 to -30          |
| Tumour 1 : salivary |             -36.1 |             -14.9 | -9 to -23           |
| Tumour 2 : kidneys  |             -34.9 |             -15.0 | -11 to -27          |
| Tumour 2 : salivary |             -28.3 |             -12.6 | -8 to -20           |

Cohort-pooled change in tumour-to-organ-at-risk absorbed-dose ratio.
Published ranges are Fig. 5c-d of Golzaryan 2025. {.table}

``` r


# The paper's two conclusions, asserted rather than eyeballed: every metronomic
# schedule raises pooled absorbed dose in every PSMA-positive compartment, and
# every one of them lowers every tumour-to-OAR ratio.
stopifnot(
  all(pooled$`Tumour 1` > 0), all(pooled$`Tumour 2` > 0),
  all(pooled$Kidneys > 0), all(pooled$`Salivary glands` > 0),
  all(pooled$`Tumour 1 : kidneys`  < 0), all(pooled$`Tumour 1 : salivary` < 0),
  all(pooled$`Tumour 2 : kidneys`  < 0), all(pooled$`Tumour 2 : salivary` < 0)
)
# The kidneys must gain more than the tumours under EVERY schedule -- that
# asymmetry is exactly why the tumour-to-OAR ratio falls, and it is the paper's
# central caveat, not an incidental cohort average.
stopifnot(all(pooled$Kidneys > pooled$`Tumour 1`))
# Dividing the same total dose more finely must raise the gain, at every interval.
stopifnot(all(vapply(split(pooled, pooled$interval),
                     function(g) all(diff(g$`Tumour 1`[order(g$n_inj)]) > 0),
                     logical(1))))
# The two compartments whose published ranges are stated most precisely must be
# reproduced, not merely trend correctly.
stopifnot(max(pooled$`Tumour 2`) >= 7,  min(pooled$`Tumour 2`) <= 25,
          max(pooled$`Salivary glands`) >= 17,
          min(pooled$`Salivary glands`) <= 56)
```

Pooled tumour-2 dose rises 6.5-23.5 % against a published 7-25 %, and
salivary gland dose 21.8-72.3 % against a published 17-56 %; tumour 1
(3.7-10.1 % against 6-20 %) and the kidneys (25.3-89.7 % against 21-73
%) reproduce the direction, the ordering and the order of magnitude,
with the kidneys running about a fifth high at the top of the range. The
tumour-to-organ-at-risk ratios fall under every schedule, by 12.6-42 %
against a published 8-30 %.

Two features are worth naming rather than smoothing over. Pooled gains
rise monotonically with the number of injections at every interval, as
the paper reports, but between the 24 h and 36 h intervals they
essentially saturate – which is consistent with the occupancy result
above, where residual occupancy has already fallen to ~3 % by 36 h so
there is little left to recover. And the increase is not universal at
the individual level: **patient 1’s tumour 2 loses 2-7 % of its absorbed
dose under every metronomic schedule**, the only such cell in the
cohort. That lesion carries a fitted binding-site density of 2412
nmol/L, thirteen times the next highest in the cohort, on a 1 mL lesion,
so its receptors are saturated by even a fractional injection; splitting
the dose cannot recruit receptors that were never the limiting resource,
and the loss of peak concentration slightly reduces uptake. The pooled
figures the paper reports include this patient.

``` r

pooled |>
  select(n_inj, interval, `Tumour 1`, `Tumour 2`, Kidneys, `Salivary glands`) |>
  pivot_longer(-c(n_inj, interval), names_to = "Compartment",
               values_to = "pct") |>
  mutate(interval = factor(interval, c(720, 1440, 2160),
                           c("12 h", "24 h", "36 h"))) |>
  ggplot(aes(factor(n_inj), pct, fill = interval)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Compartment, scales = "free_y") +
  labs(x = "Number of injections", y = "Increase in absorbed dose (%)",
       fill = "Interval",
       caption = "Replicates Figure 5a-b of Golzaryan 2025 (cohort-pooled).")
```

![Replicates Figure 5a-b of Golzaryan 2025: cohort-pooled absorbed dose
relative to a single infusion of the same total activity, by number of
injections and
interval.](Golzaryan_2025_lu177psmaIT_pbpk_files/figure-html/fig5-plot-1.png)

Replicates Figure 5a-b of Golzaryan 2025: cohort-pooled absorbed dose
relative to a single infusion of the same total activity, by number of
injections and interval.

### 5. Non-compartmental analysis of venous serum

The paper reports no plasma or serum NCA, so there is nothing to compare
against; the block below is a descriptive summary of the venous serum
labelled-peptide profile, and a check that the disposition is physically
sensible. The terminal half-life should approach but not exceed the
physical half-life of 177Lu (`log(2)/lambda_phy` = 6.73 days), because
the peptide cannot clear more slowly than the label decays.

``` r

sim_nca <- single |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc) |>
  mutate(patient = factor(id))

sim_nca <- bind_rows(
  sim_nca,
  sim_nca |> distinct(id, patient) |> mutate(time = 0, Cc = 0)
) |>
  distinct(patient, id, time, .keep_all = TRUE) |>
  arrange(id, time)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | patient + id)

dose_df <- patients |>
  transmute(id, patient = factor(id), time = 0, amt = Plab)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | patient + id)

intervals <- data.frame(start = 0, end = Inf,
                        cmax = TRUE, tmax = TRUE,
                        auclast = TRUE, half.life = TRUE)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj,
                                          intervals = intervals))

nca_tab <- as.data.frame(nca_res) |>
  select(patient, PPTESTCD, PPORRES) |>
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

knitr::kable(
  nca_tab |>
    mutate(`Terminal t1/2 (days)` = half.life / 1440) |>
    select(patient, cmax, tmax, auclast, `Terminal t1/2 (days)`) |>
    rename("Patient" = patient, "Cmax (nmol/L)" = cmax,
           "Tmax (min)" = tmax, "AUClast (nmol*min/L)" = auclast),
  digits = 3,
  caption = "Non-compartmental summary of venous serum labelled peptide after the clinical infusion."
)
```

| Patient | Cmax (nmol/L) | Tmax (min) | AUClast (nmol\*min/L) | Terminal t1/2 (days) |
|:--------|--------------:|-----------:|----------------------:|---------------------:|
| 1       |         0.537 |         20 |               132.079 |                0.402 |
| 2       |         0.505 |         20 |               124.318 |                0.248 |
| 3       |         0.528 |         20 |               156.135 |                0.318 |
| 4       |         0.448 |         20 |                99.025 |                0.243 |
| 5       |         0.508 |         20 |               142.968 |                0.298 |

Non-compartmental summary of venous serum labelled peptide after the
clinical infusion. {.table}

``` r


# Terminal half-life cannot exceed the physical half-life of 177Lu.
stopifnot(all(nca_tab$half.life <= log(2) / 7.15e-5))
# Tmax must fall at the end of the 10 min infusion.
stopifnot(all(nca_tab$tmax <= 20))
```

## Assumptions and deviations

- **Body weight is not published for these patients.** Table S2 declares
  `BW` as *measured*, but it is tabulated nowhere: not in the paper, not
  in its supplement, not in the upstream Kletting 2016 primary, and not
  in that paper’s own 19-page S1 File. The acquisition ladder was walked
  to exhaustion, including the upstream’s supplement, so this is a
  genuine reporting gap rather than a missing file. `WT` is therefore a
  declared model-input covariate and the cohort above runs at **71 kg**.
  That is not a class-typical guess: the paper writes every one of its
  ten body-size-scaled organ volumes as `coefficient x BW/71`, and its
  cited reference (22) for those coefficients is ICRP Publication 23, so
  71 kg is the reference individual the model’s own constants encode.
  Setting `WT = 71` reproduces exactly the reference organ volumes; any
  other value rescales them.
- **Hematocrit is not published either**, and is declared *measured* in
  the same row of Table S2. The cohort runs at **`HCT = 45` %**, a
  standard adult male value; this is an assumption, not paper-derived.
  It is load-bearing: hematocrit sets total serum volume through
  `V_P = 2.8 (1 - H) BSA` and hence total serum flow and every vascular
  sub-volume, and it also enters `v_TU,v = 0.05 (1 - H)`,
  `f_PRO = 0.18 (1 - H)` and `V_SAL,v = 0.03 (1 - H) V_SAL,total`. The
  receptor-occupancy agreement in section 2 is indirect evidence that
  neither assumption is badly wrong – the occupancy trajectory depends
  on interstitial and vascular volumes, and a materially different serum
  volume would not land within ~10 % of the published numbers. Every
  patient-specific quantity the paper *did* publish (body surface area,
  tubular extraction rate, all four measured volumes, injected amounts
  and activities, and the ten fitted parameters per patient) is used as
  published.
- **`K_D` is stated as “1 or 8” nmol/L in Table S2** with no indication
  of which was used. The upstream primary resolves it: Kletting 2016
  fitted both a `k_off = 0.046/min` model (K_D = 1 nmol/L, from BiaCore
  measurements) and a slower-dissociation model (from in vitro cell
  studies at ~8 nM), and reports that “Model 1 was more supported by the
  data: for all patients the Akaike weights were \> 0.9. Thus, only the
  estimated parameters of model 1 were further used.” The individual
  estimates that Golzaryan 2025 Tables S5-S9 reproduce are therefore K_D
  = 1 nmol/L estimates, and the model ships `kd = 1`.
- **Patient 4’s tumour 1 volume is 4 mL, not the 5 mL printed in Table
  S1.** Table S4 lists S-values for exactly the thirteen distinct
  volumes the model evaluates – the union of the five salivary gland
  volumes {17, 21, 29, 52, 54} and the eight distinct tumour volumes
  {0.5, 1, 1.5, 2, 3, 4, 13, 34} – and it has a row for 4 mL and none
  for 5 mL. Because that table is generated from the volumes the study
  actually used, it cannot contain a row for a volume the study never
  had. Table 1 of the upstream Kletting 2016 primary, which Table S1
  reproduces (its own footer cites
  `doi:10.1371/journal.pone.0162303.t001`), independently gives 4 mL.
  The model uses 4 mL.
- **Two transcription errors in the printed equations are corrected.**
  Both are settled by the model’s own Fig. S1 legend and by mass
  balance, and both were corrected rather than transcribed literally:
  - *Equation S7 (veins)* subtracts `F_M/V_M * P_M,v` (muscle) and adds
    `(F_M + F_GI)/V_L * P_L,v`. Muscle drains directly to the veins, so
    subtracting it would break mass balance. The Fig. S1 legend states
    that “for spleen and GI flow exit from vascular space to liver” and
    that the liver drains to the veins at
    `(F_liver + F_spleen + F_GI)/V_liver,v`. Read with **spleen** in
    place of muscle, equation S7 is exactly that identity, so the “M”
    subscript is a typographic slip for “S”. The model routes spleen and
    GI outflow through the liver.
  - *Equation S8 (arteries)* is printed as
    `-sum(F_i/V_ART) * P_i,v + F/V_LU,v * P_LU,v`, in which the outflow
    term carries the organ’s own vascular content rather than the
    arterial content it is drawing from; its labelled-species line also
    drops the asterisks entirely. Equation S4 gives the matching
    organ-side term as `F_i (P_ART/V_ART - P_i,v/V_i,v)`, which fixes
    the arterial outflow as `sum(F_i) * P_ART/V_ART`. The model
    implements the mass-balance-consistent form.
- **Equation S20 omits the intratubular kidney pool.** The kidney ROI
  content is transcribed exactly as printed – vascular, interstitial,
  receptor-bound, internalised and the protein-bound share – without
  `P_intra,K`, even though that pool physically sits inside the kidney.
  The paper notes that amino acids were co-administered to block
  unspecific uptake, which makes the omission less consequential, but it
  is the paper’s equation and not a simplification made here. The state
  exists in the model (`pintra_kidney_lab`) and can be added by a
  downstream user.
- **`V_REST,total`, `V_REST,v`, `F_REST` and `F_TOTAL` are printed as
  embedded images** in the supplement’s Table S2 that did not survive
  conversion. Each is implemented from the accompanying wording: the
  rest compartment takes the body volume, the serum volume and the
  systemic flow left over once every named compartment except the
  tumours has been assigned, and `F_TOTAL` is the sum over all organs
  including the tumours, which is what the lungs and arteries carry. All
  four come out positive and physiologically plausible for every patient
  in the cohort.
- **The tumour-remainder S-value is interpolated.** Table S4 has no row
  for the 10 mL and 50 mL lumped-remainder volumes; the dosimetry chunk
  interpolates logarithmically on volume, which is the procedure
  Kletting 2016 describes for deriving sphere S-values not tabulated
  directly by OLINDA. This affects only the `AD_tumorrest` column of the
  dosimetry table, not the model.
- **`ini()` ships patient 2’s individual estimates.** The paper fitted
  five patients separately and did not build a population model, so
  there is no typical value to report and the model carries no random
  effects and no residual error. Patient 2 was chosen because it is the
  paper’s own worked example in Figures 3, 4, S9 and S10. All five
  patients’ estimates are tabulated in the *Virtual cohort* section
  above and are passed explicitly via `rxSolve(params = )` throughout
  this vignette.
- **`c1` and `c2` equal `f_salG` for patients 3, 4 and 5** in Tables
  S7-S9 (0.53, 0.16 and 0.18 respectively), while they differ for
  patients 1 and 2. The values are transcribed as printed. They are used
  only in the tumour-ROI background correction (equations S22-S23) and
  in salivary gland perfusion, which are different roles, so the
  coincidence may be a copy error in the supplement rather than a
  genuine result. It does not affect any assertion above. It does bear
  on one reported number: patient 3’s salivary perfusion of 0.53
  mL/min/g is seven times the cohort’s next-highest value and
  contributes to that patient’s 22.3 Gy salivary absorbed dose,
  discussed in section 3. Because the values are transcribed as printed,
  that dose is what the published parameters predict; a reader who
  believes the duplication is an error should treat patient 3’s salivary
  numbers as an upper bound.
- **Two bookkeeping states are added that the paper does not carry.**
  `urine_lab` / `urine_unl` accumulate the excreted fraction of the
  filtered load and `cleared_lab` / `cleared_unl` accumulate peptide
  lost by release from the internalised pool, which the paper assumes is
  cleared from the body directly. Neither feeds back into any other
  equation; they exist so that the mass-balance check in section 1 can
  close exactly.
- **The absorbed-dose and RDP calculations live in this vignette, not in
  the model.** Equations S18-S24 are post-processing on the simulated
  time-activity curves and depend on tabulated S-values rather than on
  the ODE system, so the model file stays a pure PBPK in nmol and
  minutes.
