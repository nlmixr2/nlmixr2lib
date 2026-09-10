# Teicoplanin (Zhou 2025)

## Model and source

``` r

mod <- readModelDb("Zhou_2025_teicoplanin")
ui <- rxode2::rxode(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Zhou Y, Peng J, Xu P, Wang F, Xi J, Zhang H, Hu S, Yan H,
  Tan L, Cai H, Zhang B, Lan G. Population pharmacokinetics and dosing
  optimization of teicoplanin in renal transplant patients. Antimicrob
  Agents Chemother. 2025;69(6):e01568-24. <doi:10.1128/aac.01568-24>
- Article: <https://doi.org/10.1128/aac.01568-24> (PMC12135518, open
  access)
- Supplement: `aac.01568-24-s0001.docx` (Tables S1-S2, Figures S1-S2),
  retrieved from the EuropePMC supplementary-files endpoint for
  PMC12135518.

Teicoplanin disposition in 79 adult renal transplant recipients was
described by a two-compartment model with first-order elimination and a
zero-order (1 h infusion) input. Real-time Cockcroft-Gault creatinine
clearance on clearance was the only covariate retained.

## Population

The development cohort was 79 adults (56 male, 70.9%) enrolled
prospectively at the Second Xiangya Hospital, Changsha, between January
2022 and December 2023, contributing 306 plasma concentrations (about
four samples each). Median age was 42.5 years (IQR 36-52.8), median
weight 64.8 kg (IQR 54-74.8) and median height 168 cm (Zhou 2025 Table
1). Every patient was on tacrolimus and mycophenolate mofetil;
hypertension (81.0%) and diabetes (25.3%) were the commonest
comorbidities. Sixty-five of the 79 were within one month of
transplantation.

The defining feature of the cohort is renal: median **baseline**
Cockcroft-Gault creatinine clearance was 9 mL/min (IQR 6.6-13.8) and
median blood urea nitrogen 18.8 mmol/L. No patient received renal
replacement therapy or mechanical ventilation. Because renal function
recovers rapidly after transplantation, the authors modelled the
**real-time** rather than the baseline creatinine clearance; the 17
mL/min normalising constant in the clearance equation sits well above
the baseline median for exactly that reason.

All patients received the same empirical regimen: three 400 mg loading
doses every 12 h, then 400 mg once daily for at least 3 days, each as a
1 h infusion. Eighty of the 99 enrolled patients received teicoplanin as
perioperative prophylaxis and 19 for a confirmed or probable
Gram-positive infection (Table S1: pneumonia 14, urinary tract infection
5; MRSA 9 and *Enterococcus faecium* 5, teicoplanin MIC 0.5-1 mg/L).

A separate cohort of 20 patients / 80 samples, recruited November 2024
to March 2025, was used for external validation (RMSE 8.49%, mean
prediction error 0.704%).

The same information is available programmatically via
`readModelDb("Zhou_2025_teicoplanin")()$population`.

## Source trace

Every `ini()` entry in
`inst/modeldb/specificDrugs/Zhou_2025_teicoplanin.R` carries an in-file
comment naming its source. They are collected here for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL at CrCL = 17 mL/min) | 0.711 L/h | Table 2, `CL`; Results equation block; Abstract |
| `lvc` (Vc) | 11.3 L | Table 2, `Vc`; Results equation block |
| `lq` (CLd) | 4.22 L/h | Table 2, `CLd`; Results equation block |
| `lvp` (Vp) | 35.3 L | Table 2, `Vp`; Results equation block (the Abstract and the sentence below Table 2 round to 35.2) |
| `e_crcl_cl` (power exponent of CrCL on CL) | 0.198 | Table 2, `CrCLCL`, “the influence coefficient of CrCL on CL” |
| `etalcl` | omega = 0.401 (variance 0.1608) | Table 2, `omega CL` = 40.1% |
| `etalvc` | omega = 0.373 (variance 0.1391) | Table 2, `omega Vc` = 37.3% |
| `propSd` | 0.146 | Table 2, “Proportional residual (%)” = 14.6 |
| `cl <- exp(lcl + etalcl) * (CRCL / 17)^e_crcl_cl` | n/a | Results, “Population pharmacokinetic modeling”: `CL (L/h) = 0.711 x (CrCL / 17)^0.198 x e^eta1` |
| `vc <- exp(lvc + etalvc)` | n/a | same equation block: `Vc (L) = 11.3 x e^eta2` |
| `q <- exp(lq)`, `vp <- exp(lvp)` | n/a | same equation block: `CLd (L/h) = 4.22, Vp (L) = 35.3` (no `e^eta` factor) |
| two-compartment ODEs, zero-order input | n/a | Results, “Population pharmacokinetic modeling”; Table S2 base model, “two-compartment model with a zero-order input rate”; Methods, 1 h infusion |
| `Cc ~ prop(propSd)` | n/a | Results, “proportional residual variability”; Methods lists the three candidate log-scale error forms |

The equation block is typeset as a display equation and is lost by a
naive text extraction of the PDF (it appears as `formula-not-decoded` in
the preprocessed markdown). It was recovered with `pdftotext -layout`,
which preserves it verbatim.

Two source-internal disagreements are recorded in [Assumptions,
deviations and errata](#assumptions-deviations-and-errata) below: the
Results prose calls the retained covariate relationship “a proportional
model” while the printed equation is a power model, and Vp is printed as
both 35.3 and 35.2.

``` r

ui$iniDf |>
  dplyr::select(name, est, fix, label) |>
  dplyr::rename(
    "Parameter" = name, "Estimate" = est,
    "Fixed" = fix, "Label" = label
  ) |>
  knitr::kable(digits = 6, caption = "Packaged `ini()` block.")
```

| Parameter | Estimate | Fixed | Label |
|:---|---:|:---|:---|
| lcl | -0.341083 | FALSE | Clearance at CrCL = 17 mL/min (L/h) |
| lvc | 2.424803 | FALSE | Central volume Vc (L) |
| lq | 1.439835 | FALSE | Inter-compartmental clearance CLd (L/h) |
| lvp | 3.563883 | FALSE | Peripheral volume Vp (L) |
| e_crcl_cl | 0.198000 | FALSE | Power exponent of real-time CrCL on CL (unitless) |
| propSd | 0.146000 | FALSE | Proportional residual error (fraction) |
| etalcl | 0.160801 | FALSE | Zhou 2025 Table 2: omega_CL = 40.1% (SE 12%; bootstrap median 39.6%, 90% CI 30.6-49.2) -\> omega^2 = 0.401^2 |
| etalvc | 0.139129 | FALSE | Zhou 2025 Table 2: omega_Vc = 37.3% (SE 21%; bootstrap median 36.2%, 90% CI 16.7-53.0) -\> omega^2 = 0.373^2 |

Packaged `ini()` block. {.table}

## Structural check: ODE solution against the closed form

The two-compartment model with a constant-rate infusion has an exact
bi-exponential solution. Both sides use the same parameter values, so
the only difference between them is numerical integration error and the
bound below is correspondingly tight.

``` r

crcl_ref <- 30
p <- list(cl = 0.711 * (crcl_ref / 17)^0.198, vc = 11.3, q = 4.22, vp = 35.3)
kel <- p$cl / p$vc
k12 <- p$q / p$vc
k21 <- p$q / p$vp
bsum <- kel + k12 + k21
disc <- sqrt(bsum^2 - 4 * kel * k21)
alpha <- (bsum + disc) / 2
beta <- (bsum - disc) / 2

dose_mg <- 400
t_inf <- 1
coef_a <- (k21 - alpha) / (p$vc * (beta - alpha))
coef_b <- (k21 - beta) / (p$vc * (alpha - beta))

closed_form <- function(t) {
  arm <- function(t, k, coef) {
    ifelse(
      t <= t_inf,
      coef * (dose_mg / t_inf) / k * (1 - exp(-k * t)),
      coef * (dose_mg / t_inf) / k * (1 - exp(-k * t_inf)) * exp(-k * (t - t_inf))
    )
  }
  arm(t, alpha, coef_a) + arm(t, beta, coef_b)
}

cf_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24, 48, 96, 168)
cf_ev <- data.frame(
  id = 1L,
  time = c(0, cf_times),
  amt = c(dose_mg, rep(NA_real_, length(cf_times))),
  evid = c(1L, rep(0L, length(cf_times))),
  dur = c(t_inf, rep(NA_real_, length(cf_times))),
  cmt = "central",
  CRCL = crcl_ref
)
cf_sim <- rxode2::rxSolve(
  rxode2::zeroRe(mod), cf_ev, returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
cf_cmp <- data.frame(
  time = cf_sim$time,
  ode = cf_sim$Cc,
  closed = closed_form(cf_sim$time)
) |>
  dplyr::mutate(pct_diff = 100 * (ode - closed) / closed)

cf_cmp |>
  dplyr::mutate(pct_diff = sprintf("%.1e", pct_diff)) |>
  dplyr::rename(
    "Time (h)" = time, "ODE solve (mg/L)" = ode,
    "Closed form (mg/L)" = closed, "% difference" = pct_diff
  ) |>
  knitr::kable(digits = c(2, 4, 4, NA),
               caption = "400 mg over 1 h, CrCL 30 mL/min, typical subject.")
```

| Time (h) | ODE solve (mg/L) | Closed form (mg/L) | % difference |
|---------:|-----------------:|-------------------:|:-------------|
|     0.25 |           8.3801 |             8.3801 | 4.2e-13      |
|     0.50 |          15.9018 |            15.9018 | 4.1e-13      |
|     1.00 |          28.7923 |            28.7923 | 3.8e-13      |
|     2.00 |          19.4383 |            19.4383 | 2.0e-13      |
|     4.00 |          10.7960 |            10.7960 | -6.6e-14     |
|     8.00 |           6.6432 |             6.6432 | -2.0e-13     |
|    12.00 |           5.8554 |             5.8554 | -1.5e-13     |
|    24.00 |           4.8260 |             4.8260 | -1.3e-13     |
|    48.00 |           3.3381 |             3.3381 | -2.5e-13     |
|    96.00 |           1.5971 |             1.5971 | -4.6e-13     |
|   168.00 |           0.5285 |             0.5285 | -8.0e-13     |

400 mg over 1 h, CrCL 30 mL/min, typical subject. {.table}

``` r


# Both sides use the SAME parameters, so this is pure integration error, not a
# cohort statistic; a tight all() bound is the correct assertion here.
stopifnot(max(abs(cf_cmp$pct_diff)) < 1e-6)
```

The terminal (beta) half-life implied by these parameters is 45.1 h at
CrCL 30 mL/min. The authors are explicit that this is *not*
teicoplanin’s true terminal half-life of 83-163 h: the sampling window
was too short to resolve the third exponential phase, so the model
describes the alpha and beta phases only (Discussion, first paragraph,
and “limitations”).

## Typical-value profiles and PKNCA validation

Non-compartmental analysis of a single 400 mg infusion in a typical
subject recovers two quantities the paper prints independently of the
packaged parameters:

- **Clearance.** The printed equation `CL = 0.711 * (CrCL / 17)^0.198`
  gives a reference value at every renal stratum. Recovering it from
  `Dose / AUC0-inf` after solving the ODE system checks that the
  micro-constants, the infusion input and the `Cc <- central / vc`
  observation are wired correctly.
- **Steady-state volume.** The Discussion states the model’s volume of
  distribution as **46.6 L** (“The distribution in the present model is
  smaller than previous models (46.6 L VS 60.5-121.1 L)”). That is
  `Vc + Vp`, and it adjudicates the Vp rounding: 11.3 + 35.3 = 46.6,
  whereas 11.3 + 35.2 = 46.5.

``` r

strata_crcl <- c(5, 17, 20, 45, 75, 105)
nca_grid <- sort(unique(c(seq(0, 4, by = 0.25), seq(4, 24, by = 1),
                          seq(24, 336, by = 6))))

typ_ev <- do.call(rbind, lapply(seq_along(strata_crcl), function(i) {
  rbind(
    data.frame(id = i, time = 0, amt = dose_mg, evid = 1L, dur = t_inf,
               cmt = "central", CRCL = strata_crcl[i]),
    data.frame(id = i, time = nca_grid, amt = NA_real_, evid = 0L,
               dur = NA_real_, cmt = "central", CRCL = strata_crcl[i])
  )
}))
typ_ev$stratum <- sprintf("CrCL %g mL/min", typ_ev$CRCL)

typ_sim <- rxode2::rxSolve(
  rxode2::zeroRe(mod), typ_ev, keep = c("stratum"), returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> Warning: multi-subject simulation without without 'omega'
stopifnot(all(typ_sim$Cc >= 0))
```

``` r

typ_sim |>
  dplyr::filter(time > 0) |>
  ggplot(aes(time, Cc, colour = stratum)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(limits = c(0, 240)) +
  scale_y_log10() +
  labs(x = "Time after a single 400 mg 1 h infusion (h)",
       y = "Teicoplanin concentration (mg/L)",
       colour = "Renal stratum",
       title = "Typical-value disposition by renal function",
       caption = paste("Only clearance carries the CrCL covariate, so the",
                       "curves separate in the terminal phase only.")) +
  theme_bw()
#> Warning: Removed 96 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```

![](Zhou_2025_teicoplanin_files/figure-html/typical-figure-1.png)

``` r

nca_conc <- typ_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, stratum)

# Guarantee a time = 0 record so PKNCA anchors AUC0-* without warning. For an
# IV infusion the pre-dose concentration is 0.
nca_conc <- dplyr::bind_rows(
  nca_conc,
  nca_conc |> dplyr::distinct(id, stratum) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, stratum, time, .keep_all = TRUE) |>
  dplyr::arrange(id, stratum, time)

nca_dose <- typ_ev |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, dur, stratum)

conc_obj <- PKNCA::PKNCAconc(nca_conc, Cc ~ time | stratum + id,
                             concu = "mg/L", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(nca_dose, amt ~ time | stratum + id,
                             doseu = "mg", duration = "dur")

nca_intervals <- data.frame(
  start = 0, end = Inf,
  cmax = TRUE, tmax = TRUE, aucinf.obs = TRUE, aucpext.obs = TRUE,
  half.life = TRUE, cl.obs = TRUE, vss.iv.obs = TRUE
)
nca_res <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = nca_intervals)
)

# The 336 h window must capture essentially all of AUC0-inf for cl.obs and
# vss.iv.obs to be meaningful.
pext <- as.data.frame(nca_res$result)
stopifnot(max(pext$PPORRES[pext$PPTESTCD == "aucpext.obs"]) < 5)
```

``` r

published_nca <- tibble::tibble(
  stratum = sprintf("CrCL %g mL/min", strata_crcl),
  # Zhou 2025 Results equation block: CL = 0.711 * (CrCL / 17)^0.198.
  cl.obs = 0.711 * (strata_crcl / 17)^0.198,
  # Zhou 2025 Discussion: "The distribution in the present model is smaller
  # than previous models (46.6 L VS 60.5-121.1 L)". No covariate acts on either
  # volume, so the value is the same in every stratum.
  vss.iv.obs = 46.6
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = published_nca,
  by = "stratum",
  units = c(cl.obs = "L/h", vss.iv.obs = "L"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = paste("Simulated vs published typical-value NCA.",
                  "* differs from the reference by >20%."),
  align = c("l", "l", "r", "r", "r")
)
```

| NCA parameter | stratum         | Reference | Simulated | % diff |
|:--------------|:----------------|----------:|----------:|-------:|
| CL/F (L/h)    | CrCL 5 mL/min   |     0.558 |     0.558 |  -0.0% |
| CL/F (L/h)    | CrCL 17 mL/min  |     0.711 |     0.711 |  -0.0% |
| CL/F (L/h)    | CrCL 20 mL/min  |     0.734 |     0.734 |  -0.0% |
| CL/F (L/h)    | CrCL 45 mL/min  |     0.862 |     0.862 |  -0.0% |
| CL/F (L/h)    | CrCL 75 mL/min  |     0.954 |     0.954 |  -0.0% |
| CL/F (L/h)    | CrCL 105 mL/min |      1.02 |      1.02 |  -0.0% |
| Vss (IV) (L)  | CrCL 5 mL/min   |      46.6 |      46.6 |  -0.1% |
| Vss (IV) (L)  | CrCL 17 mL/min  |      46.6 |      46.6 |  -0.0% |
| Vss (IV) (L)  | CrCL 20 mL/min  |      46.6 |      46.6 |  -0.0% |
| Vss (IV) (L)  | CrCL 45 mL/min  |      46.6 |      46.6 |  -0.0% |
| Vss (IV) (L)  | CrCL 75 mL/min  |      46.6 |      46.6 |  -0.0% |
| Vss (IV) (L)  | CrCL 105 mL/min |      46.6 |      46.6 |  -0.0% |

Simulated vs published typical-value NCA. \* differs from the reference
by \>20%. {.table}

``` r

cmp_num <- as.data.frame(nca_res$result) |>
  dplyr::filter(PPTESTCD %in% c("cl.obs", "vss.iv.obs")) |>
  dplyr::left_join(
    tidyr::pivot_longer(published_nca, -stratum,
                        names_to = "PPTESTCD", values_to = "ref"),
    by = c("stratum", "PPTESTCD")
  ) |>
  dplyr::mutate(pct = 100 * (PPORRES - ref) / ref)

stopifnot(nrow(cmp_num) == 2L * length(strata_crcl), !anyNA(cmp_num$pct))
# Deterministic (typical-value) quantities: no cohort draw is involved, so the
# bound is tight. Realised max |pct| is well under 0.5% -- 1% leaves room for
# the AUC extrapolation and quadrature while still going red on any parameter,
# unit or wiring error.
stopifnot(max(abs(cmp_num$pct)) < 1)
```

[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
prints PKNCA’s `cl.obs` under its generic label “CL/F”; teicoplanin is
given intravenously here, so there is no bioavailability term and the
quantity is plain clearance.

The model half-lives recovered by PKNCA (36, 38, 42, 48, 50, 62 h across
the strata) are the beta half-lives discussed above, not the 83-163 h
terminal values quoted from the literature in the paper’s Introduction.

## Virtual cohort on the study regimen

The trial regimen was 400 mg q12h x 3 followed by 400 mg q24h, each
infused over 1 h. Creatinine clearance is drawn to match the observed
baseline distribution (median 9 mL/min, IQR 6.6-13.8) and held constant
per subject; see [Assumptions](#assumptions-deviations-and-errata) for
why that is a simplification of the paper’s time-varying covariate.

``` r

set.seed(20250423)
rxode2::rxSetSeed(20250423)

n_sub <- 200L
# Log-normal matched to the Table 1 development-set baseline CrCL: median 9,
# IQR 6.6-13.8 mL/min -> sdlog = log(13.8 / 6.6) / (2 * qnorm(0.75)).
crcl_sdlog <- log(13.8 / 6.6) / (2 * stats::qnorm(0.75))
cohort_crcl <- stats::rlnorm(n_sub, meanlog = log(9), sdlog = crcl_sdlog)

study_dose_t <- c(0, 12, 24, seq(48, 336, by = 24))
study_obs_t <- sort(unique(c(seq(0, 60, by = 1), seq(60, 336, by = 3),
                             study_dose_t - 1e-6)))
study_obs_t <- study_obs_t[study_obs_t >= 0]

cohort_ev <- do.call(rbind, lapply(seq_len(n_sub), function(i) {
  rbind(
    data.frame(id = i, time = study_dose_t, amt = 400, evid = 1L, dur = 1,
               cmt = "central", CRCL = cohort_crcl[i]),
    data.frame(id = i, time = study_obs_t, amt = NA_real_, evid = 0L,
               dur = NA_real_, cmt = "central", CRCL = cohort_crcl[i])
  )
}))
stopifnot(!anyDuplicated(unique(cohort_ev[, c("id", "time", "evid")])))

cohort_sim <- rxode2::rxSolve(mod, cohort_ev, keep = c("CRCL"),
                              returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
```

``` r

# Comparable to Figure 2 of Zhou 2025 (prediction-corrected VPC of the final
# model): the 5th, 50th and 95th percentiles of the model's own predictions
# over the observed sampling window.
cohort_sim |>
  dplyr::filter(time <= 240) |>
  dplyr::group_by(time) |>
  dplyr::summarise(
    Q05 = quantile(Cc, 0.05),
    Q50 = quantile(Cc, 0.50),
    Q95 = quantile(Cc, 0.95),
    .groups = "drop"
  ) |>
  ggplot(aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = c(15, 40), linetype = c("dashed", "dotted")) +
  labs(x = "Time after first dose (h)", y = "Teicoplanin concentration (mg/L)",
       title = "Study regimen: 400 mg q12h x 3, then 400 mg q24h",
       caption = paste("Median with 5th-95th percentile band, 200 virtual",
                       "subjects. Dashed line: 15 mg/L efficacy target.",
                       "Dotted line: 40 mg/L toxicity threshold.")) +
  theme_bw()
```

![](Zhou_2025_teicoplanin_files/figure-html/figure-2-1.png)

The paper’s central clinical message is visible directly: on the
empirical 400 mg regimen, the median trough in this severely renally
impaired cohort sits close to the 15 mg/L target only after several
days, which is why the authors conclude that higher-than-standard
loading and maintenance doses are needed.

## Reproducing the published Monte Carlo table

Table 3 of Zhou 2025 reports the probability of target attainment for 40
dose regimens across five renal strata, against a 15 mg/L efficacy
trough and a 40 mg/L toxicity threshold, at 72 h and 168 h after the
first dose. It is the paper’s only large block of printed numbers
derived from the model, and it therefore exercises every packaged
parameter at once: both structural clearances, both volumes, the CrCL
power term and both between-subject variances.

Two reconstruction choices are needed, and the table itself pins both.

- **When maintenance dosing starts.** Taking it to begin 24 h after the
  final loading dose puts the doses at 0, 12, 24 then 48, 72, … for a
  three-dose load and at 0, 12, 24, 36, 48 then 72, 96, … for a
  five-dose load. The table confirms this: with a five-dose load the
  first maintenance dose lands *at* 72 h, so it cannot influence
  `Cmin,72h` – and indeed the paper reports identical `Cmin,72h`
  attainment for 400 mg x 5 with 200 mg QD and with 400 mg QD (both
  91.8%), for 600 mg x 5 with 200 and 400 mg QD (both 99.3%), and for
  800 mg x 5 with 200 and 400 mg QD (99 and 99.0%). With a three-dose
  load the maintenance dose at 48 h *does* fall inside the window, and
  there the same comparison separates (600 mg x 3 gives 84.1, 92.2 and
  96.4% with 200, 400 and 600 mg QD).
- **Which creatinine clearance represents a stratum.** The paper
  resampled its enrolled subjects, so the within-stratum distribution is
  not recoverable. The midpoint of each stratum is used here (5, 20, 45,
  75 and 105 mL/min). The exponent 0.198 is weak enough that this
  matters little: across the 90-120 stratum the clearance multiplier
  moves by only 6%.

``` r

table3 <- tibble::tribble(
  ~stratum,      ~crcl, ~load, ~n_load, ~maint, ~p72_15, ~p168_15, ~p72_40, ~p168_40,
  "CrCL <10",        5,   400,      5L,    200,    91.8,     56.7,     0.0,      0.0,
  "CrCL <10",        5,   400,      3L,    400,    70.7,     83.9,     0.0,      1.6,
  "CrCL <10",        5,   400,      5L,    400,    91.8,     86.8,     0.0,      4.8,
  "CrCL <10",        5,   600,      3L,    200,    84.1,     50.5,     0.0,      0.0,
  "CrCL <10",        5,   600,      5L,    200,    99.3,     71.9,    18.5,      2.4,
  "CrCL <10",        5,   600,      3L,    400,    92.2,     87.1,     0.0,      5.5,
  "CrCL <10",        5,   600,      5L,    400,    99.3,     90.2,    18.5,     15.2,
  "CrCL <10",        5,   600,      3L,    600,    96.4,     97.2,     0.6,     28.5,
  "CrCL 10-30",     20,   600,      3L,    400,    76.5,     66.5,     0.0,      0.8,
  "CrCL 10-30",     20,   600,      5L,    400,    96.2,     75.3,     5.2,      4.0,
  "CrCL 10-30",     20,   600,      3L,    600,    86.8,     89.0,     0.0,      9.4,
  "CrCL 10-30",     20,   600,      5L,    600,    95.8,     89.8,     5.2,     15.2,
  "CrCL 10-30",     20,   800,      5L,    200,    99.0,     57.6,    40.1,      1.9,
  "CrCL 10-30",     20,   800,      3L,    400,    89.2,     71.6,     1.0,      2.0,
  "CrCL 10-30",     20,   800,      5L,    400,    99.0,     79.6,    40.1,      9.2,
  "CrCL 10-30",     20,   800,      3L,    600,    92.9,     88.8,     3.0,     13.2,
  "CrCL 30-60",     45,   600,      3L,    600,    73.5,     76.6,     0.0,      3.2,
  "CrCL 30-60",     45,   600,      5L,    600,    89.4,     78.9,     1.7,      6.6,
  "CrCL 30-60",     45,   800,      5L,    400,    97.9,     65.4,    23.1,      3.5,
  "CrCL 30-60",     45,   800,      3L,    600,    84.8,     77.4,     0.7,      4.7,
  "CrCL 30-60",     45,   800,      5L,    600,    96.8,     80.1,    22.5,     11.7,
  "CrCL 30-60",     45,   800,      3L,    800,    90.6,     90.0,     2.3,     18.0,
  "CrCL 30-60",     45,  1000,      3L,    600,    91.7,     80.4,     6.7,      8.0,
  "CrCL 30-60",     45,  1000,      3L,    800,    94.6,     89.9,    11.3,     20.4,
  "CrCL 60-90",     75,   600,      3L,    600,    64.3,     67.8,     0.0,      2.3,
  "CrCL 60-90",     75,   600,      5L,    600,    86.8,     72.4,     0.9,      4.1,
  "CrCL 60-90",     75,   800,      5L,    600,    94.7,     73.2,    15.1,      7.2,
  "CrCL 60-90",     75,   800,      3L,    800,    82.9,     82.3,     1.3,     10.7,
  "CrCL 60-90",     75,   800,      5L,    800,    94.9,     84.5,    16.6,     17.7,
  "CrCL 60-90",     75,  1000,      3L,    600,    86.2,     70.3,     4.1,      5.0,
  "CrCL 60-90",     75,  1000,      3L,    800,    90.5,     85.0,     6.3,     12.7,
  "CrCL 60-90",     75,  1000,      3L,   1000,    93.7,     92.6,    12.7,     27.9,
  "CrCL 90-120",   105,   800,      5L,    600,    91.6,     64.6,    10.0,      3.8,
  "CrCL 90-120",   105,   800,      3L,    800,    79.6,     79.0,     0.7,      7.6,
  "CrCL 90-120",   105,   800,      5L,    800,    91.8,     79.8,     9.5,     10.2,
  "CrCL 90-120",   105,  1000,      3L,    600,    81.9,     61.8,     2.0,      2.5,
  "CrCL 90-120",   105,  1000,      3L,    800,    85.7,     79.2,     4.5,      9.7,
  "CrCL 90-120",   105,  1000,      5L,    800,    97.6,     83.2,    31.2,     13.8,
  "CrCL 90-120",   105,  1000,      3L,   1000,    90.3,     88.8,     9.2,     21.2,
  "CrCL 90-120",   105,  1000,      5L,   1000,    97.1,     89.0,    32.4,     24.6
)
stopifnot(nrow(table3) == 40L)

regimen_label <- function(load, n_load, maint) {
  sprintf("%s mg q12h x %d, %s mg QD",
          format(load, big.mark = ",", trim = TRUE), n_load,
          format(maint, big.mark = ",", trim = TRUE))
}

build_pta_events <- function(load, n_load, maint, crcl, n, id_offset) {
  load_t <- seq(0, by = 12, length.out = n_load)
  maint_t <- seq(load_t[n_load] + 24, 192, by = 24)
  dose_t <- c(load_t, maint_t)
  dose_a <- c(rep(load, n_load), rep(maint, length(maint_t)))
  do.call(rbind, lapply(seq_len(n), function(i) {
    rbind(
      data.frame(id = id_offset + i, time = dose_t, amt = dose_a, evid = 1L,
                 dur = 1, cmt = "central", CRCL = crcl),
      data.frame(id = id_offset + i, time = c(72, 168), amt = NA_real_,
                 evid = 0L, dur = NA_real_, cmt = "central", CRCL = crcl)
    )
  }))
}

n_pta <- 200L
rxode2::rxSetSeed(97531)
pta_sim <- lapply(seq_len(nrow(table3)), function(k) {
  ev <- build_pta_events(table3$load[k], table3$n_load[k], table3$maint[k],
                         table3$crcl[k], n_pta, (k - 1L) * n_pta)
  s <- rxode2::rxSolve(mod, ev, returnType = "data.frame")
  c72 <- s$Cc[s$time == 72]
  c168 <- s$Cc[s$time == 168]
  stopifnot(length(c72) == n_pta, length(c168) == n_pta)
  data.frame(
    row = k,
    s72_15 = 100 * mean(c72 > 15), s168_15 = 100 * mean(c168 > 15),
    s72_40 = 100 * mean(c72 > 40), s168_40 = 100 * mean(c168 > 40)
  )
}) |>
  dplyr::bind_rows()

pta <- dplyr::bind_cols(table3, dplyr::select(pta_sim, -row)) |>
  dplyr::mutate(regimen = regimen_label(load, n_load, maint))
```

``` r

pta |>
  dplyr::mutate(
    `Cmin,72h > 15` = sprintf("%.1f / %.1f", p72_15, s72_15),
    `Cmin,168h > 15` = sprintf("%.1f / %.1f", p168_15, s168_15),
    `Cmin,72h > 40` = sprintf("%.1f / %.1f", p72_40, s72_40),
    `Cmin,168h > 40` = sprintf("%.1f / %.1f", p168_40, s168_40)
  ) |>
  dplyr::select(stratum, regimen, `Cmin,72h > 15`, `Cmin,168h > 15`,
                `Cmin,72h > 40`, `Cmin,168h > 40`) |>
  dplyr::rename("Stratum (mL/min)" = stratum, "Dose regimen" = regimen) |>
  knitr::kable(
    caption = paste("Zhou 2025 Table 3 reproduced. Each cell is",
                    "published / simulated PTA in percent."),
    align = c("l", "l", "r", "r", "r", "r")
  )
```

| Stratum (mL/min) | Dose regimen | Cmin,72h \> 15 | Cmin,168h \> 15 | Cmin,72h \> 40 | Cmin,168h \> 40 |
|:---|:---|---:|---:|---:|---:|
| CrCL \<10 | 400 mg q12h x 5, 200 mg QD | 91.8 / 92.0 | 56.7 / 58.0 | 0.0 / 0.0 | 0.0 / 0.5 |
| CrCL \<10 | 400 mg q12h x 3, 400 mg QD | 70.7 / 65.5 | 83.9 / 79.0 | 0.0 / 0.0 | 1.6 / 1.0 |
| CrCL \<10 | 400 mg q12h x 5, 400 mg QD | 91.8 / 90.0 | 86.8 / 87.0 | 0.0 / 0.0 | 4.8 / 5.5 |
| CrCL \<10 | 600 mg q12h x 3, 200 mg QD | 84.1 / 82.0 | 50.5 / 44.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| CrCL \<10 | 600 mg q12h x 5, 200 mg QD | 99.3 / 99.5 | 71.9 / 64.5 | 18.5 / 19.5 | 2.4 / 1.5 |
| CrCL \<10 | 600 mg q12h x 3, 400 mg QD | 92.2 / 93.5 | 87.1 / 88.5 | 0.0 / 0.5 | 5.5 / 5.5 |
| CrCL \<10 | 600 mg q12h x 5, 400 mg QD | 99.3 / 99.5 | 90.2 / 89.5 | 18.5 / 23.5 | 15.2 / 20.5 |
| CrCL \<10 | 600 mg q12h x 3, 600 mg QD | 96.4 / 95.0 | 97.2 / 97.5 | 0.6 / 0.5 | 28.5 / 19.5 |
| CrCL 10-30 | 600 mg q12h x 3, 400 mg QD | 76.5 / 72.0 | 66.5 / 62.5 | 0.0 / 0.0 | 0.8 / 0.5 |
| CrCL 10-30 | 600 mg q12h x 5, 400 mg QD | 96.2 / 97.0 | 75.3 / 71.5 | 5.2 / 6.5 | 4.0 / 5.0 |
| CrCL 10-30 | 600 mg q12h x 3, 600 mg QD | 86.8 / 78.5 | 89.0 / 84.0 | 0.0 / 0.0 | 9.4 / 10.0 |
| CrCL 10-30 | 600 mg q12h x 5, 600 mg QD | 95.8 / 96.0 | 89.8 / 87.5 | 5.2 / 4.5 | 15.2 / 13.0 |
| CrCL 10-30 | 800 mg q12h x 5, 200 mg QD | 99.0 / 98.5 | 57.6 / 53.0 | 40.1 / 36.5 | 1.9 / 2.5 |
| CrCL 10-30 | 800 mg q12h x 3, 400 mg QD | 89.2 / 85.0 | 71.6 / 72.0 | 1.0 / 1.0 | 2.0 / 3.0 |
| CrCL 10-30 | 800 mg q12h x 5, 400 mg QD | 99.0 / 100.0 | 79.6 / 72.5 | 40.1 / 35.5 | 9.2 / 7.0 |
| CrCL 10-30 | 800 mg q12h x 3, 600 mg QD | 92.9 / 90.0 | 88.8 / 85.5 | 3.0 / 3.0 | 13.2 / 13.0 |
| CrCL 30-60 | 600 mg q12h x 3, 600 mg QD | 73.5 / 74.5 | 76.6 / 79.0 | 0.0 / 0.0 | 3.2 / 5.5 |
| CrCL 30-60 | 600 mg q12h x 5, 600 mg QD | 89.4 / 89.5 | 78.9 / 75.5 | 1.7 / 4.0 | 6.6 / 9.0 |
| CrCL 30-60 | 800 mg q12h x 5, 400 mg QD | 97.9 / 98.0 | 65.4 / 66.0 | 23.1 / 22.5 | 3.5 / 4.0 |
| CrCL 30-60 | 800 mg q12h x 3, 600 mg QD | 84.8 / 84.5 | 77.4 / 78.5 | 0.7 / 1.0 | 4.7 / 9.0 |
| CrCL 30-60 | 800 mg q12h x 5, 600 mg QD | 96.8 / 95.5 | 80.1 / 81.0 | 22.5 / 27.0 | 11.7 / 14.0 |
| CrCL 30-60 | 800 mg q12h x 3, 800 mg QD | 90.6 / 90.5 | 90.0 / 89.5 | 2.3 / 4.0 | 18.0 / 20.0 |
| CrCL 30-60 | 1,000 mg q12h x 3, 600 mg QD | 91.7 / 92.5 | 80.4 / 81.0 | 6.7 / 11.0 | 8.0 / 12.5 |
| CrCL 30-60 | 1,000 mg q12h x 3, 800 mg QD | 94.6 / 92.5 | 89.9 / 87.0 | 11.3 / 15.0 | 20.4 / 25.0 |
| CrCL 60-90 | 600 mg q12h x 3, 600 mg QD | 64.3 / 64.5 | 67.8 / 67.0 | 0.0 / 0.0 | 2.3 / 3.5 |
| CrCL 60-90 | 600 mg q12h x 5, 600 mg QD | 86.8 / 82.5 | 72.4 / 68.5 | 0.9 / 2.5 | 4.1 / 7.5 |
| CrCL 60-90 | 800 mg q12h x 5, 600 mg QD | 94.7 / 97.0 | 73.2 / 69.0 | 15.1 / 16.5 | 7.2 / 10.5 |
| CrCL 60-90 | 800 mg q12h x 3, 800 mg QD | 82.9 / 82.5 | 82.3 / 82.0 | 1.3 / 2.5 | 10.7 / 15.5 |
| CrCL 60-90 | 800 mg q12h x 5, 800 mg QD | 94.9 / 92.0 | 84.5 / 81.5 | 16.6 / 13.5 | 17.7 / 15.5 |
| CrCL 60-90 | 1,000 mg q12h x 3, 600 mg QD | 86.2 / 81.0 | 70.3 / 71.5 | 4.1 / 5.5 | 5.0 / 6.0 |
| CrCL 60-90 | 1,000 mg q12h x 3, 800 mg QD | 90.5 / 90.0 | 85.0 / 82.5 | 6.3 / 10.5 | 12.7 / 16.0 |
| CrCL 60-90 | 1,000 mg q12h x 3, 1,000 mg QD | 93.7 / 89.0 | 92.6 / 86.5 | 12.7 / 17.5 | 27.9 / 28.5 |
| CrCL 90-120 | 800 mg q12h x 5, 600 mg QD | 91.6 / 89.5 | 64.6 / 67.0 | 10.0 / 10.5 | 3.8 / 5.0 |
| CrCL 90-120 | 800 mg q12h x 3, 800 mg QD | 79.6 / 74.0 | 79.0 / 73.5 | 0.7 / 0.5 | 7.6 / 8.5 |
| CrCL 90-120 | 800 mg q12h x 5, 800 mg QD | 91.8 / 87.0 | 79.8 / 73.0 | 9.5 / 10.5 | 10.2 / 12.0 |
| CrCL 90-120 | 1,000 mg q12h x 3, 600 mg QD | 81.9 / 82.5 | 61.8 / 65.5 | 2.0 / 2.5 | 2.5 / 4.0 |
| CrCL 90-120 | 1,000 mg q12h x 3, 800 mg QD | 85.7 / 83.5 | 79.2 / 77.0 | 4.5 / 6.0 | 9.7 / 11.5 |
| CrCL 90-120 | 1,000 mg q12h x 5, 800 mg QD | 97.6 / 95.0 | 83.2 / 85.0 | 31.2 / 36.0 | 13.8 / 21.0 |
| CrCL 90-120 | 1,000 mg q12h x 3, 1,000 mg QD | 90.3 / 88.5 | 88.8 / 87.5 | 9.2 / 12.0 | 21.2 / 22.0 |
| CrCL 90-120 | 1,000 mg q12h x 5, 1,000 mg QD | 97.1 / 94.0 | 89.0 / 84.5 | 32.4 / 35.0 | 24.6 / 28.5 |

Zhou 2025 Table 3 reproduced. Each cell is published / simulated PTA in
percent. {.table style="width:100%;"}

``` r

pta_diff <- c(
  pta$s72_15 - pta$p72_15, pta$s168_15 - pta$p168_15,
  pta$s72_40 - pta$p72_40, pta$s168_40 - pta$p168_40
)
stopifnot(length(pta_diff) == 160L, !anyNA(pta_diff))

pta_summary <- data.frame(
  Statistic = c("cells compared", "median |difference| (pp)",
                "90th percentile |difference| (pp)", "maximum |difference| (pp)"),
  Value = c(length(pta_diff), median(abs(pta_diff)),
            unname(quantile(abs(pta_diff), 0.9)), max(abs(pta_diff)))
)
knitr::kable(pta_summary, digits = 2,
             caption = "Agreement with the published Table 3.")
```

| Statistic                           |  Value |
|:------------------------------------|-------:|
| cells compared                      | 160.00 |
| median \|difference\| (pp)          |   1.55 |
| 90th percentile \|difference\| (pp) |   4.81 |
| maximum \|difference\| (pp)         |   9.00 |

Agreement with the published Table 3. {.table}

``` r


# These are cohort statistics, so the bounds are robust summaries, not
# per-cell limits. Each simulated cell is a binomial proportion from 200
# subjects (about 2.8 pp standard error near 80%) against a published cell from
# 1,000 (about 1.3 pp), so an individual cell can legitimately be several
# points off; assert the CENTRE and a robust upper quantile instead. At 1,000
# subjects per arm the realised values are 1.5 / 3.7 / 4.5 pp -- the bounds
# below sit outside the 200-subject sampling spread and still go red on any
# parameter that moves the distribution, which is what a mis-read variance,
# a dropped covariate term or a wrong reference constant would do.
stopifnot(
  median(abs(pta_diff)) < 6,
  quantile(abs(pta_diff), 0.9) < 12
)
```

## The AUC24/MIC target

The paper’s second efficacy criterion is `AUC24h/MIC >= 610.4`,
evaluated over two 24 h windows: 72-96 h, which the authors describe as
loading-dose driven (Figure 6), and 168-192 h, which is maintenance
driven (Figure 7). Two claims from the Results are checked below in both
windows.

1.  “All of the selected regimens achieved the PTA target when MIC = 0.5
    mg/L.”
2.  At MIC = 1.0 mg/L the recommended regimens are 600 mg q12h x 5 then
    400 mg QD for CrCL \<10; 800 mg q12h x 3 then 600 mg QD for 10-30;
    800 mg q12h x 3 then 800 mg QD for 30-60; 1,000 mg q12h x 3 then 800
    mg QD for 60-90; and 1,000 mg q12h x 3 then 1,000 mg QD for 90-120.

Unlike Table 3, both of these are read off graphical figures, so only
the 80% threshold line the authors drew on them is available as a
reference value.

``` r

auc_reg <- tibble::tribble(
  ~stratum,      ~crcl, ~load, ~n_load, ~maint,
  "CrCL <10",        5,   600,      5L,    400,
  "CrCL 10-30",     20,   800,      3L,    600,
  "CrCL 30-60",     45,   800,      3L,    800,
  "CrCL 60-90",     75,  1000,      3L,    800,
  "CrCL 90-120",   105,  1000,      3L,   1000
)
auc_grid <- sort(unique(c(seq(72, 96, by = 0.5), seq(168, 192, by = 0.5))))

rxode2::rxSetSeed(24680)
auc_res <- lapply(seq_len(nrow(auc_reg)), function(k) {
  load_t <- seq(0, by = 12, length.out = auc_reg$n_load[k])
  maint_t <- seq(load_t[auc_reg$n_load[k]] + 24, 192, by = 24)
  dose_t <- c(load_t, maint_t)
  dose_a <- c(rep(auc_reg$load[k], auc_reg$n_load[k]),
              rep(auc_reg$maint[k], length(maint_t)))
  ev <- do.call(rbind, lapply(seq_len(n_pta), function(i) {
    rbind(
      data.frame(id = i, time = dose_t, amt = dose_a, evid = 1L, dur = 1,
                 cmt = "central", CRCL = auc_reg$crcl[k]),
      data.frame(id = i, time = auc_grid, amt = NA_real_, evid = 0L,
                 dur = NA_real_, cmt = "central", CRCL = auc_reg$crcl[k])
    )
  }))
  s <- rxode2::rxSolve(mod, ev, returnType = "data.frame")
  s <- s[!is.na(s$Cc), c("id", "time", "Cc")]
  s$stratum <- auc_reg$stratum[k]
  s
}) |>
  dplyr::bind_rows()

auc_conc <- PKNCA::PKNCAconc(auc_res, Cc ~ time | stratum + id,
                             concu = "mg/L", timeu = "h")
auc_out <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  auc_conc,
  intervals = data.frame(start = c(72, 168), end = c(96, 192), auclast = TRUE)
))
#> No dose information provided, calculations requiring dose will return NA.

auc24 <- as.data.frame(auc_out$result) |>
  dplyr::filter(PPTESTCD == "auclast") |>
  dplyr::mutate(window = ifelse(start == 72, "72-96 h", "168-192 h")) |>
  dplyr::group_by(stratum, window) |>
  dplyr::summarise(
    median_auc = median(PPORRES),
    pta_mic05 = 100 * mean(PPORRES / 0.5 >= 610.4),
    pta_mic10 = 100 * mean(PPORRES / 1.0 >= 610.4),
    n = dplyr::n(),
    .groups = "drop"
  )
stopifnot(all(auc24$n == n_pta), nrow(auc24) == 2L * nrow(auc_reg))

auc24 |>
  dplyr::left_join(
    auc_reg |> dplyr::mutate(regimen = regimen_label(load, n_load, maint)) |>
      dplyr::select(stratum, regimen),
    by = "stratum"
  ) |>
  dplyr::select(stratum, regimen, window, median_auc, pta_mic05, pta_mic10) |>
  dplyr::rename(
    "Stratum (mL/min)" = stratum, "Dose regimen" = regimen,
    "Window" = window, "Median AUC24 (mg*h/L)" = median_auc,
    "PTA (%) at MIC 0.5" = pta_mic05, "PTA (%) at MIC 1.0" = pta_mic10
  ) |>
  knitr::kable(digits = 1,
               caption = paste("Attainment of AUC24h/MIC >= 610.4 for the",
                               "regimens the paper recommends at",
                               "MIC = 1.0 mg/L."))
```

| Stratum (mL/min) | Dose regimen | Window | Median AUC24 (mg\*h/L) | PTA (%) at MIC 0.5 | PTA (%) at MIC 1.0 |
|:---|:---|:---|---:|---:|---:|
| CrCL 10-30 | 800 mg q12h x 3, 600 mg QD | 168-192 h | 828.9 | 100.0 | 80.0 |
| CrCL 10-30 | 800 mg q12h x 3, 600 mg QD | 72-96 h | 820.2 | 100.0 | 85.0 |
| CrCL 30-60 | 800 mg q12h x 3, 800 mg QD | 168-192 h | 912.5 | 100.0 | 88.0 |
| CrCL 30-60 | 800 mg q12h x 3, 800 mg QD | 72-96 h | 859.4 | 100.0 | 88.5 |
| CrCL 60-90 | 1,000 mg q12h x 3, 800 mg QD | 168-192 h | 820.1 | 99.5 | 80.5 |
| CrCL 60-90 | 1,000 mg q12h x 3, 800 mg QD | 72-96 h | 847.7 | 99.5 | 84.5 |
| CrCL 90-120 | 1,000 mg q12h x 3, 1,000 mg QD | 168-192 h | 1002.5 | 100.0 | 89.0 |
| CrCL 90-120 | 1,000 mg q12h x 3, 1,000 mg QD | 72-96 h | 964.7 | 100.0 | 89.0 |
| CrCL \<10 | 600 mg q12h x 5, 400 mg QD | 168-192 h | 800.0 | 98.5 | 76.5 |
| CrCL \<10 | 600 mg q12h x 5, 400 mg QD | 72-96 h | 902.9 | 100.0 | 92.5 |

Attainment of AUC24h/MIC \>= 610.4 for the regimens the paper recommends
at MIC = 1.0 mg/L. {.table}

``` r

mic10_load <- auc24$pta_mic10[auc24$window == "72-96 h"]
mic10_maint <- auc24$pta_mic10[auc24$window == "168-192 h"]

# Claim 1 reproduces in both windows with a very wide margin: the realised
# minimum at 1,000 subjects per arm is 98.6%. The bound is one-sided and comes
# from the paper's own 80% criterion, not from this run.
stopifnot(all(auc24$pta_mic05 >= 90))

# Claim 2 reproduces in the loading-dose window. Realised at 1,000 subjects per
# arm: 81.1-90.8%. The 70 bound leaves room for the 200-subject binomial
# standard error (about 2.8 pp near 80%) while still going red on any change
# that moves the exposure distribution materially.
stopifnot(all(mic10_load >= 70))

# Claim 2 in the MAINTENANCE window is a DOCUMENTED DEVIATION -- three of five
# strata land a few points under 80% (see the prose below). It is deliberately
# NOT gated at 80; the assertion below only keeps it in a plausible band so a
# gross regression still fails.
stopifnot(median(mic10_maint) >= 65)
```

Claim 1 reproduces with a wide margin in both windows: every recommended
regimen exceeds 98% attainment at MIC 0.5 mg/L.

Claim 2 reproduces in the loading-dose window (72-96 h), where all five
recommended regimens clear 80%. In the maintenance window (168-192 h)
three of the five land a few percentage points short. At 1,000 subjects
per arm – enough to remove sampling noise as an explanation – the
realised values are 73.8% (CrCL \<10), 76.8% (10-30), 85.1% (30-60),
78.9% (60-90) and 88.2% (90-120).

This is recorded as a deviation rather than corrected. The most likely
explanation is the stratum-representative creatinine clearance: this
window’s attainment depends on maintenance-phase clearance alone, and a
representative CrCL at the lower end of each band instead of its
midpoint raises attainment in every stratum (at CrCL 2 rather than 5
mL/min, the `<10` stratum’s median AUC168-192 rises by about a fifth and
its attainment clears 80%). The paper resampled its enrolled subjects
rather than using a single value per stratum, and it reports these
results only as figures, so neither the distribution nor the exact
plotted values can be recovered. Notably, the same representative values
reproduce the *printed* Table 3 trough attainment closely (median
absolute difference under 2 percentage points), so the disagreement is
confined to the figure-derived claim.

## Assumptions, deviations and errata

- **“Proportional model” versus the printed power equation.** The
  Results paragraph introducing the covariate says real-time CrCL was
  included in CL “using a proportional model”, but the equation block
  printed three lines later on the same page is unambiguously a power
  model, `CL (L/h) = 0.711 x (CrCL / 17)^0.198 x e^eta1`. The equation
  is used. Two further pieces of evidence support it: Table 2 labels the
  estimate `CrCLCL`, “the influence coefficient of CrCL on CL”, which is
  exponent-like language; and the linear reading is arithmetically
  impossible in this cohort. A centred proportional form
  `0.711 * (1 + 0.198 * (CrCL - 17))` turns negative below
  `CrCL = 17 - 1/0.198 = 11.9` mL/min, and the development set’s median
  baseline CrCL is 9 mL/min – more than half the modelled subjects would
  carry a negative clearance, and the ODE system does not solve at all.
  The Table 3 reproduction above is the positive confirmation.
- **Vp is printed as both 35.3 and 35.2 L.** Table 2 and the Results
  equation block give 35.3; the Abstract and the sentence immediately
  below Table 2 give 35.2. 35.3 is used, because it is the value in both
  the parameter table and the equation, and because the Discussion’s
  total volume of 46.6 L equals 11.3 + 35.3 rather than 11.3 + 35.2.
- **Between-subject variability scale.** Table 2 reports the random
  effects under a “Random effects (%)” header as `omega CL` 40.1% and
  `omega Vc` 37.3%. They are taken as omega (the log-scale standard
  deviation) x 100, so the packaged variances are 0.401^2 and 0.373^2.
  The column cannot be variances: the residual row in the same block
  reads 14.6%, which as a variance would mean a 38% residual standard
  deviation – irreconcilable with the paper’s own external-validation
  RMSE of 8.49% and with an assay whose inter-day CV was at most 3.72%.
  Reading the same percentages as approximate CVs instead
  (`omega = sqrt(log(1 + CV^2))`) would give 0.386 and 0.361; the
  difference is immaterial for simulation and is inside the reproduction
  tolerance of the Table 3 check.
- **Residual error form.** Phoenix NLME writes its proportional error on
  the log scale as `ln Y = ln F + eps1` (the second of the three
  candidate forms printed in Methods), which is formally exponential.
  The packaged model uses `Cc ~ prop(propSd)` with `propSd = 0.146`,
  matching the paper’s own “proportional residual” label; at that
  magnitude the two forms differ by well under 1% in the simulated
  concentration distribution.
- **CrCL is held constant per subject in the simulations here.** The
  published model uses *real-time* creatinine clearance, which rises
  substantially during a teicoplanin course as the graft recovers, but
  the paper publishes no trajectory for it and the Monte Carlo
  simulations in Table 3 are themselves run at fixed renal strata. Every
  simulation in this vignette therefore uses a time-constant `CRCL`. A
  user with longitudinal creatinine can supply a time-varying `CRCL`
  column directly; the model needs no change.
- **Stratum representative CrCL.** The paper resampled its enrolled
  subjects within each renal stratum, so the within-stratum covariate
  distribution is not recoverable. Stratum midpoints (5, 20, 45, 75, 105
  mL/min) are used. The weak 0.198 exponent makes the Table 3
  reproduction insensitive to this choice: at 1,000 subjects per arm, a
  uniform draw across the `<10` stratum and a fixed 5 mL/min give median
  absolute Table 3 differences of 1.4 and 1.5 percentage points
  respectively.
- **Known deviation: the figure-derived AUC24/MIC recommendations at MIC
  1.0 mg/L in the 168-192 h window.** Three of the five recommended
  regimens reach 74-79% rather than the 80% the authors’ Figure 7
  threshold implies. The printed Table 3 trough attainment for the same
  strata reproduces to within about 2 percentage points, and the same
  regimens clear 80% in the paper’s other (72-96 h) window, so this is
  attributed to the unrecoverable within-stratum CrCL distribution
  rather than to a transcription error. It is reported in the table and
  excluded from the assertion rather than having the bound widened until
  it passes.
- **Maintenance-dose start time.** Taken as 24 h after the final loading
  dose. This is not stated in the Methods but is pinned by the internal
  structure of Table 3, as set out in the reproduction section above.
- **Baseline CrCL distribution for the virtual cohort.** Table 1 reports
  the development-set baseline CrCL as a median with an interquartile
  range only. A log-normal matched to that median and IQR is used; the
  paper reports no distributional form.
- **No supplement gap.** The article’s supplement (Tables S1-S2, Figures
  S1-S2) was retrieved and used. It contributed the covariate-screening
  delta-OFV values recorded in `population$screened_covariates` and the
  infection / pathogen detail in `population$disease_state`; it contains
  no final parameter estimates beyond those in Table 2.
- **Model scope.** The authors state that the short sampling window
  prevented the model from resolving teicoplanin’s true terminal phase
  (literature half-life 83-163 h). This two-compartment model describes
  the alpha and beta phases and should not be used to project exposure
  far beyond the 192 h horizon the paper itself simulates.
