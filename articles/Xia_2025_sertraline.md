# Sertraline (Xia 2025)

## Model and source

- Citation: Xia H, Deng G, Liang F, Zhang Z, Huang W, Guo Z, Song Q, Wen
  Y, Shang D, Tan Y. Investigating Remedial Strategies for Missed or
  Delayed Dose of Sertraline in Chinese Adolescent Patients with
  Depressive Disorders via Population Pharmacokinetics Modeling and
  Simulation Approaches. Drug Des Devel Ther. 2025;19:3001-3016.
  <doi:10.2147/DDDT.S504521>. PMID 40260198. PMCID PMC12011033. The
  fixed absorption rate constant ka = 0.5 1/h is taken from Poweleit EA,
  Taylor ZL, Mizuno T, et al. Escitalopram and sertraline population
  pharmacokinetic analysis in pediatric patients. Clin Pharmacokinet.
  2023;62(11):1621-1637. <doi:10.1007/s40262-023-01294-8> (Xia 2025
  reference 25). Xia 2025 reference 34 is the same group’s earlier
  adult-and-adolescent sertraline model, packaged here as
  modellib(‘Zhang_2024_sertraline’).
- Description: One-compartment first-order absorption population PK
  model for sertraline in Chinese adolescent (13-17 years) inpatients
  with depressive disorders (Xia 2025), built from 221 routine
  therapeutic-drug-monitoring trough concentrations in 103 patients. The
  model carries NO covariates: height, body weight, age, gender and
  concomitant quetiapine / olanzapine / alprazolam were all screened by
  forward inclusion and none reduced the objective function value by the
  required 6.63, which the authors attribute to the narrow 13-17 year
  age band and the small weight spread of the cohort. The absorption
  rate constant is held at ka = 0.5 1/h taken from Poweleit 2023 because
  every sample was an elimination-phase trough and the absorption phase
  was not identifiable. Typical CL/F = 65.8 L/h and V/F = 1570 L give
  kel = 0.0419 1/h and a 16.5 h half-life. The paper’s purpose is
  dosing-remediation simulation: it derives recommended remedial doses
  for one, two and three consecutive missed 50 / 100 / 200 mg QD doses
  as a function of how late the dose is taken.
- Article: <https://doi.org/10.2147/DDDT.S504521> (open access;
  PMC12011033)

Xia 2025 is a dosing-remediation study. Its population PK model is
deliberately plain – one compartment, first-order absorption, no
covariates – and exists to support the paper’s real product: a set of
recommended remedial doses for a Chinese adolescent who misses one, two
or three consecutive daily sertraline doses. This vignette therefore
validates the model twice: first as a PK model (structural gates against
the paper’s own printed steady-state quantities), then as the simulation
engine behind the remedial recommendations.

## Population

The model was fitted to 221 routine therapeutic-drug-monitoring (TDM)
serum sertraline concentrations from 103 Chinese adolescents aged 13-17
years who were hospitalised with depressive disorders at the Affiliated
Brain Hospital of Guangzhou Medical University between 1 January 2019
and 31 December 2023 (Table 1). The cohort was 30 male (29.0%) and 73
female (71.0%); median age was 15 years, median weight 57.5 kg in males
and 50 kg in females, and median height 169 cm in males and 159.5 cm in
females. Median daily dose was 150 mg (range 50-300 mg) and the median
observed concentration was 63.54 ug/L (range 5.07-299.87). Quetiapine
(41.75%), alprazolam (33.00%) and olanzapine (31.07%) were the common
comedications.

Two features of the dataset shape the model. Every sample is an
elimination-phase trough drawn at 6-7 a.m. before the next scheduled
dose, so the absorption phase was never observed and `ka` had to be
fixed to a literature value. And 37.90% of heights and 5.80% of weights
were missing and were replaced by within-sex medians, which blunts the
covariate screen: no covariate reached the forward inclusion threshold,
and the final model has none.

The same information is available programmatically via
`readModelDb("Xia_2025_sertraline")()$population`.

``` r

pop <- rxode2::rxode(readModelDb("Xia_2025_sertraline"))$population
#> ℹ parameter labels from comments will be replaced by 'label()'
tibble::tibble(Field = names(pop), Value = vapply(pop, function(x) paste(as.character(x), collapse = "; "), character(1))) |>
  knitr::kable(caption = "Population metadata carried with the model.")
```

| Field | Value |
|:---|:---|
| species | human |
| n_subjects | 103 |
| n_studies | 1 |
| n_concentrations | 221 |
| age_range | 13-17 years |
| age_median | 15 years |
| weight_range | 39.9-89 kg (male); 34-91 kg (female) |
| weight_median | 57.5 kg (male); 50 kg (female) |
| height_range | 150-176 cm (male); 150-172 cm (female) |
| height_median | 169 cm (male); 159.5 cm (female) |
| sex_female_pct | 71 |
| race_ethnicity | Chinese (single-centre Guangzhou cohort; non-Chinese patients were an explicit exclusion criterion, and the cohort is not further stratified by the source) |
| disease_state | Hospitalised adolescents with depressive disorders receiving oral sertraline |
| dose_range | Daily dose median 150 mg, range 50-300 mg; the simulations use 50, 100 and 200 mg QD |
| regions | China (The Affiliated Brain Hospital of Guangzhou Medical University, Guangzhou, Guangdong) |
| co_medication | Quetiapine 43 patients (41.75%); alprazolam 34 (33.00%); olanzapine 32 (31.07%) |
| sampling | Retrospective therapeutic drug monitoring; all 221 samples are elimination-phase troughs drawn at 6-7 a.m. before the next scheduled dose. Median (range) observed concentration 63.54 (5.07-299.87) ug/L. |
| notes | Baseline demographics from Xia 2025 Table 1. Retrospective TDM data collected 1 January 2019 to 31 December 2023; IRB approval 2021027. Serum sertraline quantified by HPLC-MS/MS over a 5-500 ng/mL calibrated range with intra- and inter-day precision below 15% RSE; the AGNP therapeutic reference range used is 10-150 ng/mL and the laboratory alert level is 300 ng/mL. Inclusion required at least two sertraline concentrations per patient within one hospitalisation. Missing heights (37.90%) and weights (5.80%) were replaced by within-sex medians. Model fitted in NONMEM 7.3.0 with FOCE-I; evaluated by goodness-of-fit plots, NPDE (variance 1.14, mean 0.152) and a 1000-run bootstrap that converged in 928 runs (92.8%). |

Population metadata carried with the model. {.table}

## Source trace

Every `ini()` entry carries an in-file comment naming its source
location; the table below collects them for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lka` (ka) | 0.5 1/h, fixed | Table 2, `Ka (h-1) = 0.5 fixed`; Methods “Model Development” fixes ka to reference 25 (Poweleit 2023) because no sample was drawn in the absorption phase |
| `lcl` (CL/F) | 65.8 L/h | Table 2, Final Model Estimate; RSE 6%; bootstrap median 65.96 (95% CI 59.11-72.82) |
| `lvc` (V/F) | 1570 L | Table 2, Final Model Estimate; RSE 19%; bootstrap median 1652.18 (95% CI 1129.10-2401.14) |
| `etalcl` | 0.134 | Table 2, `IIV (%)` on CL/F, read as the NONMEM `$OMEGA` variance (CV 37.9%) – see Assumptions |
| `etalvc` | 0.326 | Table 2, `IIV (%)` on V/F, read as the NONMEM `$OMEGA` variance (CV 62.1%) – see Assumptions |
| `propSd` | `sqrt(0.0959)` = 0.3097 | Table 2, `PRO (%)` read as the NONMEM `$SIGMA` variance; bootstrap median 0.09 (95% CI 0.01-0.12) |
| `addSd` | `sqrt(24.7)` = 4.970 ug/L | Table 2, `ADD (%)` read as the NONMEM `$SIGMA` variance; bootstrap median 26.28 (95% CI 2.68-344.63) |
| IIV model `P = Ptv * exp(eta)` | n/a | Equation (1), page 3003 |
| Residual model `Y = F * (1 + eps1) + eps2` | n/a | Equation (2), page 3003 |
| One-compartment first-order absorption ODEs | n/a | Methods “Model Development”: “a one-compartment model”; Results “one-compartment PPK model of primary absorption” |
| Covariate forms (screened, none retained) | n/a | Equations (3) and (4), page 3003; Results “Model Establishment” |
| Therapeutic reference range 10-150 ug/L; alert level 300 ug/L | n/a | Methods “Determination of Sertraline Concentration” |
| Steady-state troughs 22.34 / 44.69 / 89.37 ug/L | n/a | Results “Model Simulation” |
| Steady-state peak at 4.5 h post-dose | n/a | Results “Model Simulation” |
| CL/V = 0.042 1/h | n/a | Discussion, paragraph on elimination rate |

Equations (1)-(4) are vector graphics in the published PDF and are
invisible to text extraction; they were read from a 200 dpi
rasterisation of page 3003 and from the publisher’s equation images in
the EuropePMC bundle (`DDDT-19-3001-e0001.jpg` .. `e0004.jpg`).

## Setup: dosing conventions

Xia 2025 doses at 08:30 every morning and draws the TDM sample at 6-7
a.m. before the next dose. That sample is therefore **21.5-22.5 h after
the previous dose, not 24 h**, which matters below: the paper’s
“steady-state trough” is the concentration at the TDM clock time, not
the true pre-dose minimum.

Steady state is established with a single `ss = 1` record rather than a
multi-week dosing run-in. That is exact for this linear model (verified
below to 4e-07 against an explicit 33-day history) and makes every
steady-state quantity independent of how long a run-in happened to be
simulated.

``` r

DAY      <- 24     # dosing interval (h)
SAMPLE_H <- 21.5   # 06:00 TDM sample, 21.5 h after the 08:30 dose
LOWER_TR <- 10     # lower therapeutic reference limit (ug/L)
UPPER_TR <- 150    # upper therapeutic reference limit (ug/L)
ALERT    <- 300    # laboratory alert level (ug/L)
DOSES    <- c("50 mg QD" = 50, "100 mg QD" = 100, "200 mg QD" = 200)

mod  <- readModelDb("Xia_2025_sertraline")
modT <- rxode2::zeroRe(mod)   # typical-value model, as Xia 2025 simulated
#> ℹ parameter labels from comments will be replaced by 'label()'

# Build an event table. `cmt` on observation rows is the ODE state "central";
# rxode2 returns the algebraic observable Cc as a column at those rows.
makeEvents <- function(doseTime, doseAmt, obsTime, id = 1L, ss1 = TRUE) {
  dplyr::bind_rows(
    data.frame(
      id = id, time = doseTime, amt = doseAmt, evid = 1L, cmt = "depot",
      ii = ifelse(ss1 & doseTime == 0, DAY, 0),
      ss = ifelse(ss1 & doseTime == 0, 1L, 0L)
    ),
    data.frame(
      id = id, time = obsTime, amt = NA_real_, evid = 0L, cmt = "central",
      ii = 0, ss = 0L
    )
  ) |>
    dplyr::arrange(time, dplyr::desc(evid))
}

solveTypical <- function(ev) {
  suppressWarnings(rxode2::rxSolve(modT, ev, keep = intersect("scenario", names(ev)),
                                   returnType = "data.frame"))
}
concAt <- function(s, t) stats::approx(s$time, s$Cc, xout = t)$y
```

``` r

# The ss = 1 shortcut against an explicit 33-day dosing history.
obsLong  <- seq(32 * DAY, 44 * DAY, by = 0.1)
runIn    <- solveTypical(makeEvents(DAY * (0:44), 50, obsLong, ss1 = FALSE))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
ssShort  <- solveTypical(makeEvents(DAY * (0:12), 50, seq(0, 12 * DAY, by = 0.1)))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
ssRelDiff <- max(abs(ssShort$Cc - runIn$Cc) / runIn$Cc)
stopifnot(ssRelDiff < 1e-5)
sprintf("ss = 1 vs a 33-day dosing run-in: max relative difference %.2g", ssRelDiff)
#> [1] "ss = 1 vs a 33-day dosing run-in: max relative difference 4.1e-07"
```

## Structural verification

All checks in this section are **deterministic** – they compare a
typical-value (`zeroRe`) solve against a closed form or against a number
Xia 2025 printed. There is no simulated cohort and therefore no cohort
noise, so tight tolerances are correct here (a mis-transcribed
clearance, volume or dose moves these by tens of percent).

``` r

obsGrid <- sort(unique(c(seq(0, 3 * DAY, by = 0.05), SAMPLE_H, DAY, 2 * DAY)))

typicalSS <- dplyr::bind_rows(lapply(seq_along(DOSES), function(i) {
  s <- solveTypical(makeEvents(0, DOSES[[i]], obsGrid, id = i))
  s$treatment <- names(DOSES)[i]
  s$dose      <- DOSES[[i]]
  s
}))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

ssSummary <- typicalSS |>
  dplyr::filter(time <= DAY) |>
  dplyr::group_by(treatment, dose) |>
  dplyr::summarise(
    cmax      = max(Cc),
    tmax      = time[which.max(Cc)],
    trough06  = Cc[which.min(abs(time - SAMPLE_H))],
    trough24  = Cc[which.min(abs(time - DAY))],
    .groups   = "drop"
  ) |>
  dplyr::arrange(dose)
ssSummary
#> # A tibble: 3 × 6
#>   treatment  dose  cmax  tmax trough06 trough24
#>   <chr>     <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 50 mg QD     50  41.7   4.4     22.3     20.0
#> 2 100 mg QD   100  83.4   4.4     44.5     40.1
#> 3 200 mg QD   200 167.    4.4     89.0     80.2
```

### Elimination rate against the paper’s own derived value

``` r

kelModel <- 65.8 / 1570
sprintf("model kel = CL/F / (V/F) = %.5f 1/h; Xia 2025 Discussion quotes CL/V = 0.042 1/h", kelModel)
#> [1] "model kel = CL/F / (V/F) = 0.04191 1/h; Xia 2025 Discussion quotes CL/V = 0.042 1/h"
sprintf("terminal half-life = %.2f h", log(2) / kelModel)
#> [1] "terminal half-life = 16.54 h"
stopifnot(abs(kelModel - 0.042) < 0.0005)
```

### ODE solve against the closed-form steady-state superposition

This gate is a pure implementation check: both sides use the same
parameters, so the only difference is numerical. It also confirms that
rxode2’s automatic `cl`/`vc` analytic solution agrees with the
explicitly written ODE system.

``` r

cssClosed <- function(D, t, cl = 65.8, v = 1570, ka = 0.5) {
  kel <- cl / v
  (D * 1000 * ka / (v * (ka - kel))) *
    (exp(-kel * t) / (1 - exp(-kel * DAY)) - exp(-ka * t) / (1 - exp(-ka * DAY)))
}
tt      <- seq(0.1, DAY, by = 0.1)
odeVals <- concAt(dplyr::filter(typicalSS, treatment == "50 mg QD"), tt)
closedRelDiff <- max(abs(odeVals - cssClosed(50, tt)) / cssClosed(50, tt))
stopifnot(closedRelDiff < 1e-4)
sprintf("ODE vs closed-form superposition (50 mg QD): max relative difference %.2g", closedRelDiff)
#> [1] "ODE vs closed-form superposition (50 mg QD): max relative difference 6e-07"
```

### Steady-state trough against the published values

Xia 2025 reports steady-state troughs of 22.34, 44.69 and 89.37 ug/L for
the 50, 100 and 200 mg QD regimens. Those numbers are only reproducible
if the sample is placed in the paper’s stated 6-7 a.m. TDM window: at
21.5 h post-dose the model gives 22.26 / 44.51 / 89.03 (a 0.4%
difference), whereas the true 24 h pre-dose minimum is 20.04 / 40.09 /
80.17 (a 10% difference). The agreement at the sampling time – and only
at the sampling time – is what validates the transcription.

``` r

publishedTrough <- c("50 mg QD" = 22.34, "100 mg QD" = 44.69, "200 mg QD" = 89.37)

troughCheck <- ssSummary |>
  dplyr::mutate(
    published   = unname(publishedTrough[treatment]),
    pct_06      = 100 * (trough06 - published) / published,
    pct_24      = 100 * (trough24 - published) / published
  ) |>
  dplyr::select(treatment, `Model at 06:00 (21.5 h)` = trough06,
                `Model at 24 h` = trough24, `Xia 2025` = published,
                `% diff at 06:00` = pct_06, `% diff at 24 h` = pct_24)
knitr::kable(troughCheck, digits = 2,
             caption = "Steady-state trough: the published value is the concentration at the 6-7 a.m. TDM sampling time, not the 24 h pre-dose minimum.")
```

| treatment | Model at 06:00 (21.5 h) | Model at 24 h | Xia 2025 | % diff at 06:00 | % diff at 24 h |
|:---|---:|---:|---:|---:|---:|
| 50 mg QD | 22.26 | 20.04 | 22.34 | -0.37 | -10.28 |
| 100 mg QD | 44.51 | 40.09 | 44.69 | -0.39 | -10.30 |
| 200 mg QD | 89.03 | 80.17 | 89.37 | -0.38 | -10.29 |

Steady-state trough: the published value is the concentration at the 6-7
a.m. TDM sampling time, not the 24 h pre-dose minimum. {.table}

``` r


# Deterministic: a mis-transcribed CL, V or dose moves this by tens of percent.
stopifnot(max(abs(troughCheck$`% diff at 06:00`)) < 2)
# All three regimens sit inside the 10-150 ug/L therapeutic reference range,
# as Xia 2025 states ("all concentrations were within the therapeutic window").
stopifnot(all(ssSummary$trough06 > LOWER_TR), all(ssSummary$trough06 < UPPER_TR))
```

### Steady-state peak time and exact dose proportionality

``` r

sprintf("model steady-state tmax = %.2f h; Xia 2025 reports the peak sample at 4.5 h post-dose",
        ssSummary$tmax[1])
#> [1] "model steady-state tmax = 4.40 h; Xia 2025 reports the peak sample at 4.5 h post-dose"
stopifnot(all(abs(ssSummary$tmax - 4.5) < 0.2))

# The model is linear, so every exposure metric must be exactly proportional to
# dose. This catches an accidental non-linearity or a mis-scaled observable.
proportional <- ssSummary$trough06 / ssSummary$dose
stopifnot(max(abs(proportional - proportional[1])) / proportional[1] < 1e-6)
sprintf("trough per mg of dose: %s ug/L/mg (identical across regimens)",
        paste(sprintf("%.6f", proportional), collapse = ", "))
#> [1] "trough per mg of dose: 0.445140, 0.445140, 0.445140 ug/L/mg (identical across regimens)"
```

``` r

typicalSS |>
  dplyr::filter(time <= 2 * DAY) |>
  ggplot(aes(time, Cc, colour = factor(dose, levels = DOSES, labels = names(DOSES)))) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(LOWER_TR, UPPER_TR), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = SAMPLE_H, linetype = "dotted") +
  labs(x = "Time after an 08:30 dose (h)", y = "Sertraline (ug/L)", colour = NULL,
       title = "Steady-state typical-value profiles",
       caption = paste("Dashed lines: 10-150 ug/L therapeutic reference range.",
                       "Dotted line: 06:00 TDM sampling time (21.5 h)."))
```

![](Xia_2025_sertraline_files/figure-html/figure-steady-state-1.png)

## PKNCA validation

``` r

ncaGrid <- sort(unique(c(seq(0, DAY, by = 0.25), SAMPLE_H)))
ncaEvents <- dplyr::bind_rows(lapply(seq_along(DOSES), function(i) {
  ev <- makeEvents(0, DOSES[[i]], ncaGrid, id = i)
  ev$treatment <- names(DOSES)[i]
  ev
}))
ncaSim <- suppressWarnings(
  rxode2::rxSolve(modT, ncaEvents, keep = "treatment", returnType = "data.frame")
)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

# Filter on !is.na(Cc) only -- a `time > 0` or `Cc > 0` filter would drop the
# time-zero row that PKNCA needs to anchor the AUC.
ncaConc <- ncaSim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, treatment)
ncaConc <- dplyr::bind_rows(
  ncaConc,
  ncaConc |> dplyr::distinct(id, treatment) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, treatment, time, .keep_all = TRUE) |>
  dplyr::arrange(id, treatment, time)
stopifnot(nrow(ncaConc) > 0, sum(ncaConc$time == 0) == length(DOSES))

concObj <- PKNCA::PKNCAconc(ncaConc, Cc ~ time | treatment + id)
doseObj <- PKNCA::PKNCAdose(
  ncaEvents |> dplyr::filter(evid == 1) |> dplyr::select(id, time, amt, treatment),
  amt ~ time | treatment + id
)

wideNca <- function(res) {
  as.data.frame(res) |>
    dplyr::select(treatment, PPTESTCD, PPORRES) |>
    tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)
}
```

### Full dosing interval: AUCtau against `Dose / CL`

At steady state `AUCtau` must equal `F * Dose / (CL/F)` exactly.
`ctrough` here is the true 24 h pre-dose minimum (a record sits exactly
on the interval end, as PKNCA requires).

``` r

ncaTau <- suppressWarnings(PKNCA::pk.nca(PKNCA::PKNCAdata(
  concObj, doseObj,
  intervals = data.frame(start = 0, end = DAY,
                         cmax = TRUE, tmax = TRUE, auclast = TRUE,
                         ctrough = TRUE, half.life = TRUE)
)))
tauTab <- wideNca(ncaTau) |>
  dplyr::mutate(dose = unname(DOSES[treatment]),
                auc_closed = dose * 1000 / 65.8,
                auc_pct = 100 * (auclast - auc_closed) / auc_closed) |>
  dplyr::arrange(dose)
tauTab |>
  dplyr::select(Regimen = treatment, `AUCtau (PKNCA)` = auclast,
                `Dose/CL` = auc_closed, `% diff` = auc_pct,
                Cmax = cmax, Tmax = tmax, `Ctrough (24 h)` = ctrough,
                `t1/2` = half.life) |>
  knitr::kable(digits = 3, caption = "Steady-state NCA over the full 24 h interval.")
```

| Regimen   | AUCtau (PKNCA) |  Dose/CL | % diff |    Cmax | Tmax | Ctrough (24 h) |  t1/2 |
|:----------|---------------:|---------:|-------:|--------:|-----:|---------------:|------:|
| 50 mg QD  |        759.791 |  759.878 | -0.011 |  41.721 |  4.5 |         20.043 | 16.73 |
| 100 mg QD |       1519.582 | 1519.757 | -0.011 |  83.442 |  4.5 |         40.087 | 16.73 |
| 200 mg QD |       3039.164 | 3039.514 | -0.011 | 166.883 |  4.5 |         80.174 | 16.73 |

Steady-state NCA over the full 24 h interval. {.table
style="width:100%;"}

``` r


stopifnot(max(abs(tauTab$auc_pct)) < 0.5)
# Terminal half-life recovered from the profile, against log(2) / kel = 16.54 h.
# The lambda-z window still carries a little absorption, so allow 5%.
stopifnot(all(abs(tauTab$half.life - log(2) / kelModel) / (log(2) / kelModel) < 0.05))
```

### Comparison against the published values

Xia 2025 publishes two quantities that map onto NCA parameters: the
steady-state trough at the 6-7 a.m. TDM sample and the 4.5 h time of the
peak. Computing `ctrough` over the interval `[0, 21.5]` makes it exactly
the concentration at that sampling time.

``` r

ncaPub <- suppressWarnings(PKNCA::pk.nca(PKNCA::PKNCAdata(
  concObj, doseObj,
  intervals = data.frame(start = 0, end = SAMPLE_H,
                         cmax = TRUE, tmax = TRUE, ctrough = TRUE)
)))

published <- tibble::tribble(
  ~treatment,   ~ctrough, ~tmax,
  "50 mg QD",      22.34,   4.5,
  "100 mg QD",     44.69,   4.5,
  "200 mg QD",     89.37,   4.5
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = ncaPub,
  reference     = published,
  by            = "treatment",
  params        = c("ctrough", "tmax"),
  units         = c(ctrough = "ug/L", tmax = "h"),
  tolerance_pct = 20
)
knitr::kable(cmp, caption = "Simulated vs published NCA (Xia 2025 Results, 'Model Simulation'). * differs from reference by >20%.")
```

| NCA parameter  | treatment | Reference | Simulated | % diff |
|:---------------|:----------|:----------|:----------|:-------|
| Tmax (h)       | 50 mg QD  | 4.5       | 4.5       | +0.0%  |
| Tmax (h)       | 100 mg QD | 4.5       | 4.5       | +0.0%  |
| Tmax (h)       | 200 mg QD | 4.5       | 4.5       | +0.0%  |
| Ctrough (ug/L) | 50 mg QD  | 22.3      | 22.3      | -0.4%  |
| Ctrough (ug/L) | 100 mg QD | 44.7      | 44.5      | -0.4%  |
| Ctrough (ug/L) | 200 mg QD | 89.4      | 89        | -0.4%  |

Simulated vs published NCA (Xia 2025 Results, ‘Model Simulation’). \*
differs from reference by \>20%. {.table}

No row is flagged; the largest deviation is the 0.4% on `ctrough`, and
`tmax` matches to the resolution of the 0.25 h simulation grid (the
continuous maximum is at 4.40 h).

## Between-subject variability

Xia 2025 ran every published simulation at typical values, so the cohort
below is not a replication of any figure in the paper – it shows the
between-subject spread the fitted `$OMEGA` implies, which is the
quantity most sensitive to the variance-versus-standard-deviation
reading discussed under Assumptions.

``` r

# rxode2's RNG streams are partitioned per solver thread, so this cohort is not
# byte-identical across machines with different thread counts. Every assertion
# below is written on a robust statistic, never on an extreme.
set.seed(20250417)
N_PER_ARM <- 200

cohortEvents <- dplyr::bind_rows(lapply(seq_along(DOSES), function(i) {
  ids <- (i - 1L) * N_PER_ARM + seq_len(N_PER_ARM)
  ev  <- dplyr::bind_rows(lapply(ids, function(k) makeEvents(0, DOSES[[i]], ncaGrid, id = k)))
  ev$treatment <- names(DOSES)[i]
  ev
}))
stopifnot(!anyDuplicated(unique(cohortEvents[, c("id", "time", "evid")])))

cohortSim <- rxode2::rxSolve(mod, cohortEvents, keep = "treatment",
                             returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
```

``` r

cohortTrough <- cohortSim |>
  dplyr::filter(abs(time - SAMPLE_H) < 1e-9) |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(
    n      = dplyr::n(),
    Q10    = quantile(Cc, 0.10),
    Median = median(Cc),
    Q90    = quantile(Cc, 0.90),
    `Below 10 ug/L (%)` = 100 * mean(Cc < LOWER_TR),
    .groups = "drop"
  ) |>
  dplyr::mutate(dose = unname(DOSES[treatment])) |>
  dplyr::arrange(dose) |>
  dplyr::select(Regimen = treatment, n, Q10, Median, Q90, `Below 10 ug/L (%)`)
knitr::kable(cohortTrough, digits = 2,
             caption = "Simulated 06:00 steady-state trough by regimen (200 subjects per arm).")
```

| Regimen   |   n |   Q10 | Median |    Q90 | Below 10 ug/L (%) |
|:----------|----:|------:|-------:|-------:|------------------:|
| 50 mg QD  | 200 |  7.86 |  21.01 |  41.06 |              15.0 |
| 100 mg QD | 200 | 20.44 |  45.88 |  76.05 |               2.5 |
| 200 mg QD | 200 | 36.04 |  82.44 | 173.69 |               0.0 |

Simulated 06:00 steady-state trough by regimen (200 subjects per arm).
{.table}

``` r


# Robust, not extreme: the cohort median sits near the typical value but the two
# log-normal etas skew the distribution, so allow generous headroom.
medianPct <- 100 * (cohortTrough$Median - ssSummary$trough06) / ssSummary$trough06
stopifnot(max(abs(medianPct)) < 20)
sprintf("cohort median vs typical value: %s", paste(sprintf("%+.1f%%", medianPct), collapse = ", "))
#> [1] "cohort median vs typical value: -5.6%, +3.1%, -7.4%"
```

``` r

cohortSim |>
  dplyr::group_by(treatment, time) |>
  dplyr::summarise(Q05 = quantile(Cc, 0.05), Q50 = median(Cc),
                   Q95 = quantile(Cc, 0.95), .groups = "drop") |>
  dplyr::mutate(treatment = factor(treatment, levels = names(DOSES))) |>
  ggplot(aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(LOWER_TR, UPPER_TR), linetype = "dashed", colour = "grey40") +
  facet_wrap(~treatment) +
  labs(x = "Time after an 08:30 dose (h)", y = "Sertraline (ug/L)",
       title = "Steady-state spread implied by the fitted between-subject variability",
       caption = "Median with 5th-95th percentile band, 200 subjects per arm. Not a figure of Xia 2025.")
```

![](Xia_2025_sertraline_files/figure-html/figure-vpc-1.png)

## Missed doses without a remedy

This section replicates Figure 3 of Xia 2025 (concentration-time curves
when the missed dose is simply skipped and the regular schedule
resumes). Xia 2025 misses days 33-35 of a long treatment course; steady
state is established here with `ss = 1`, so the first missed dose is
placed on day 3 and days 0-2 show the undisturbed steady state.

``` r

MISS_DAY <- 3L
HORIZON  <- 12L
obsMiss  <- sort(unique(c(seq(0, HORIZON * DAY, by = 0.25),
                          (0:HORIZON) * DAY + SAMPLE_H)))

missedDays <- function(nMiss) MISS_DAY + seq_len(nMiss) - 1L

buildScenario <- function(nMiss, dose, delay = NA_real_, immediate = 0,
                          nextDose = NA_real_, id = 1L) {
  keep  <- setdiff(0:HORIZON, missedDays(nMiss))
  times <- DAY * keep
  amts  <- rep(dose, length(times))
  resumeDay <- MISS_DAY + nMiss                 # first scheduled dose after the gap
  if (!is.na(nextDose)) amts[times == DAY * resumeDay] <- nextDose
  if (immediate > 0) {                          # dose taken `delay` h late on the last missed day
    times <- c(times, DAY * (resumeDay - 1L) + delay)
    amts  <- c(amts, immediate)
  }
  ord <- order(times)
  makeEvents(times[ord], amts[ord], obsMiss, id = id)
}

referenceArm <- function(dose, id) makeEvents(DAY * (0:HORIZON), dose, obsMiss, id = id)

noRemedyGrid <- tidyr::expand_grid(nMiss = 1:3, dose = unname(DOSES))
noRemedyEvents <- dplyr::bind_rows(lapply(seq_len(nrow(noRemedyGrid)), function(i) {
  ev <- buildScenario(noRemedyGrid$nMiss[i], noRemedyGrid$dose[i], id = i)
  ev$scenario <- i
  ev
}))
noRemedySim <- solveTypical(noRemedyEvents) |>
  dplyr::left_join(dplyr::mutate(noRemedyGrid, scenario = dplyr::row_number()), by = "scenario")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

refEvents <- dplyr::bind_rows(lapply(seq_along(DOSES), function(i) {
  ev <- referenceArm(DOSES[[i]], id = i)
  ev$scenario <- i
  ev
}))
refSim <- solveTypical(refEvents) |>
  dplyr::mutate(dose = unname(DOSES[scenario]))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
refTroughAt <- function(dose, day) {
  concAt(dplyr::filter(refSim, dose == !!dose), day * DAY + SAMPLE_H)
}
```

``` r

noRemedySim |>
  dplyr::mutate(
    panel = factor(sprintf("%s, %d missed", names(DOSES)[match(dose, DOSES)], nMiss),
                   levels = as.vector(t(outer(names(DOSES), 1:3,
                                              function(a, b) sprintf("%s, %d missed", a, b)))))
  ) |>
  ggplot(aes(time / DAY, Cc)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = LOWER_TR, linetype = "dashed", colour = "red") +
  facet_wrap(~panel, scales = "free_y", ncol = 3) +
  labs(x = "Day (first missed dose on day 3)", y = "Sertraline (ug/L)",
       title = "Missed doses with no remedial action",
       caption = paste("Replicates Figure 3 of Xia 2025 (A-C single, D-F double, G-I triple missed dose).",
                       "Red dashed line: 10 ug/L lower therapeutic limit."))
```

![](Xia_2025_sertraline_files/figure-html/figure-3-1.png)

``` r

noRemedyStats <- noRemedySim |>
  dplyr::group_by(nMiss, dose) |>
  dplyr::summarise(
    nadir = min(Cc[time >= MISS_DAY * DAY & time <= (MISS_DAY + nMiss[1] + 1L) * DAY]),
    .groups = "drop"
  ) |>
  dplyr::mutate(refTrough = mapply(refTroughAt, dose, MISS_DAY + nMiss),
                worstTroughPct = NA_real_)

for (i in seq_len(nrow(noRemedyStats))) {
  s <- dplyr::filter(noRemedySim, nMiss == noRemedyStats$nMiss[i], dose == noRemedyStats$dose[i])
  d <- noRemedyStats$nMiss[i] + MISS_DAY
  noRemedyStats$worstTroughPct[i] <-
    100 * (concAt(s, d * DAY + SAMPLE_H) - noRemedyStats$refTrough[i]) / noRemedyStats$refTrough[i]
}

knitr::kable(
  noRemedyStats |>
    dplyr::select(`Missed doses` = nMiss, `Dose (mg)` = dose,
                  `Nadir (ug/L)` = nadir, `First trough vs steady state (%)` = worstTroughPct),
  digits = 2,
  caption = "Depth of the excursion when a missed dose is simply skipped."
)
```

| Missed doses | Dose (mg) | Nadir (ug/L) | First trough vs steady state (%) |
|-------------:|----------:|-------------:|---------------------------------:|
|            1 |        50 |         7.33 |                           -23.20 |
|            1 |       100 |        14.66 |                           -23.20 |
|            1 |       200 |        29.32 |                           -23.20 |
|            2 |        50 |         2.68 |                           -31.68 |
|            2 |       100 |         5.36 |                           -31.68 |
|            2 |       200 |        10.72 |                           -31.68 |
|            3 |        50 |         0.98 |                           -34.79 |
|            3 |       100 |         1.96 |                           -34.79 |
|            3 |       200 |         3.92 |                           -34.79 |

Depth of the excursion when a missed dose is simply skipped. {.table}

``` r


# Xia 2025: "There was an increased likelihood of a patient's serum
# concentration falling below the treatment threshold following a missed dose"
# and, for three missed doses, "the sertraline serum concentration was below the
# lower limit of the effective therapeutic".
stopifnot(all(noRemedyStats$nadir[noRemedyStats$dose == 50] < LOWER_TR))
stopifnot(all(noRemedyStats$worstTroughPct < -20))
# Deeper excursion for more missed doses, at every dose level. Arrange
# explicitly so the check does not depend on the incoming row order.
depthOrdered <- dplyr::arrange(noRemedyStats, dose, nMiss)
stopifnot(all(tapply(depthOrdered$worstTroughPct, depthOrdered$dose,
                     function(x) length(x) == 3L && all(diff(x) < 0))))
```

### How long recovery takes without a remedy

Xia 2025 states that “it took at least five days for sertraline
steady-state serum concentrations to restore if patient did not take a
remedy after a missed dose”. The measured recovery time depends on how
close to steady state counts as “restored”, so both a 5% and a 1%
criterion are reported.

``` r

recoveryDays <- function(nMiss, dose, tol) {
  s <- dplyr::filter(noRemedySim, nMiss == !!nMiss, dose == !!dose)
  days <- (MISS_DAY + nMiss):HORIZON
  dev  <- vapply(days, function(d) {
    abs(concAt(s, d * DAY + SAMPLE_H) - refTroughAt(dose, d)) / refTroughAt(dose, d)
  }, numeric(1))
  hit <- which(dev < tol)[1]
  if (is.na(hit)) NA_real_ else days[hit] - MISS_DAY
}
recovery <- tidyr::expand_grid(nMiss = 1:3, dose = 50) |>
  dplyr::mutate(`Within 5% (days)` = mapply(recoveryDays, nMiss, dose, 0.05),
                `Within 1% (days)` = mapply(recoveryDays, nMiss, dose, 0.01))
knitr::kable(dplyr::select(recovery, `Missed doses` = nMiss,
                           `Within 5% (days)`, `Within 1% (days)`),
             caption = "Days from the first missed dose until the 06:00 trough returns to steady state (recovery is dose-independent because the model is linear).")
```

| Missed doses | Within 5% (days) | Within 1% (days) |
|-------------:|-----------------:|-----------------:|
|            1 |                3 |                5 |
|            2 |                4 |                6 |
|            3 |                5 |                7 |

Days from the first missed dose until the 06:00 trough returns to steady
state (recovery is dose-independent because the model is linear).
{.table}

``` r


# The paper's "approximately five days" is a visual read of Figure 3; the
# criterion-dependent range below brackets it.
stopifnot(all(recovery$`Within 5% (days)` >= 2), all(recovery$`Within 1% (days)` <= 8))
```

Recovery from a *single* missed dose takes 3 days at the 5% criterion
and 4 days at the 1% criterion, against the paper’s stated
“approximately five days”; a *triple* missed dose takes exactly 5 days
at the 5% criterion. The claim is reproduced in magnitude, and the
residual difference is a matter of where one declares the curve to have
rejoined the steady-state band by eye. This is recorded as a deviation
rather than gated tightly.

## Remedial strategies

The paper’s recommendations are stated in the Results as delay windows
with a dose taken immediately and a dose at the next scheduled 08:30.
Two rules run across all of them: if the delay is within 14 h the missed
dose is taken immediately, and if the resulting single remedial dose
would exceed 200 mg the authors decline to recommend a remedy at all.

``` r

remedial <- tibble::tribble(
  ~nMiss, ~dose, ~window,     ~delay, ~immediate, ~nextDose, ~figure,
  1L,  50, "<= 7 h",              4,   50,   50.0, "4A",
  1L,  50, "7-14 h",             10,   50,   37.5, "4B",
  1L,  50, "> 14 h",             18,    0,   75.0, "4C",
  1L, 100, "<= 5 h",              3,  100,  100.0, "4D",
  1L, 100, "5-14 h",             10,  100,   75.0, "4E",
  1L, 100, "> 14 h",             18,    0,  150.0, "4F",
  1L, 200, "<= 4 h",              2,  200,  200.0, "4G",
  1L, 200, "4-9 h",               6,  200,  175.0, "4H",
  1L, 200, "9-14 h",             11,  200,  150.0, "4I",
  1L, 200, "14-24 h",            18,    0,  200.0, "4J",
  2L,  50, "<= 14 h",            10,   50,   62.5, "5A",
  2L,  50, "14-24 h",            18,    0,   87.5, "5B",
  2L, 100, "<= 14 h",            10,  100,  125.0, "5C",
  2L, 100, "14-24 h",            18,    0,  175.0, "5D",
  2L, 200, "<= 2 h",              1,  200,  225.0, "5E",
  2L, 200, "> 2 h",               6,    0,  300.0, "5F",
  3L,  50, "<= 6 h",              4,   50,   75.0, "6A",
  3L,  50, "6-24 h",             12,    0,  100.0, "6B",
  3L, 100, "<= 6 h",              4,  100,  150.0, "6C",
  3L, 100, "6-24 h",             12,    0,  200.0, "6D",
  3L, 200, "<= 14 h",            10,  200,  250.0, "6E",
  3L, 200, "14-24 h",            18,    0,  300.0, "6F"
)
stopifnot(nrow(remedial) == 22L)

remEvents <- dplyr::bind_rows(lapply(seq_len(nrow(remedial)), function(i) {
  r  <- remedial[i, ]
  ev <- buildScenario(r$nMiss, r$dose, r$delay, r$immediate, r$nextDose, id = i)
  ev$scenario <- i
  ev
}))
remSim <- solveTypical(remEvents)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
```

``` r

remedialEval <- remedial |>
  dplyr::mutate(scenario = dplyr::row_number(), resumeDay = MISS_DAY + nMiss)

remedialEval$refTrough <- mapply(refTroughAt, remedialEval$dose, remedialEval$resumeDay)
remedialEval$remTrough <- vapply(seq_len(nrow(remedialEval)), function(i) {
  concAt(dplyr::filter(remSim, scenario == i), remedialEval$resumeDay[i] * DAY + SAMPLE_H)
}, numeric(1))
remedialEval$noRemTrough <- vapply(seq_len(nrow(remedialEval)), function(i) {
  s <- dplyr::filter(noRemedySim, nMiss == remedialEval$nMiss[i], dose == remedialEval$dose[i])
  concAt(s, remedialEval$resumeDay[i] * DAY + SAMPLE_H)
}, numeric(1))
remedialEval$peak <- vapply(seq_len(nrow(remedialEval)), function(i) {
  s <- dplyr::filter(remSim, scenario == i)
  max(s$Cc[s$time >= MISS_DAY * DAY & s$time <= (remedialEval$resumeDay[i] + 2L) * DAY])
}, numeric(1))

remedialEval <- remedialEval |>
  dplyr::mutate(
    remPct   = 100 * (remTrough - refTrough) / refTrough,
    noRemPct = 100 * (noRemTrough - refTrough) / refTrough,
    adjusts  = !(immediate == 0 & nextDose == dose)   # 4J prescribes no adjustment
  )

knitr::kable(
  remedialEval |>
    dplyr::select(Fig = figure, `Missed` = nMiss, `Dose (mg)` = dose, Delay = window,
                  `Immediate (mg)` = immediate, `Next dose (mg)` = nextDose,
                  `Trough vs SS (%)` = remPct, `If skipped (%)` = noRemPct,
                  `Peak (ug/L)` = peak),
  digits = 1,
  caption = "Every remedial recommendation of Xia 2025 Figures 4-6, simulated from the packaged model."
)
```

| Fig | Missed | Dose (mg) | Delay | Immediate (mg) | Next dose (mg) | Trough vs SS (%) | If skipped (%) | Peak (ug/L) |
|:---|---:|---:|:---|---:|---:|---:|---:|---:|
| 4A | 1 | 50 | \<= 7 h | 50 | 50.0 | 4.2 | -23.2 | 43.7 |
| 4B | 1 | 50 | 7-14 h | 50 | 37.5 | -3.8 | -23.2 | 41.1 |
| 4C | 1 | 50 | \> 14 h | 0 | 75.0 | 8.5 | -23.2 | 43.9 |
| 4D | 1 | 100 | \<= 5 h | 100 | 100.0 | 3.1 | -23.2 | 86.3 |
| 4E | 1 | 100 | 5-14 h | 100 | 75.0 | -3.8 | -23.2 | 82.2 |
| 4F | 1 | 100 | \> 14 h | 0 | 150.0 | 8.5 | -23.2 | 87.9 |
| 4G | 1 | 200 | \<= 4 h | 200 | 200.0 | 2.0 | -23.2 | 170.6 |
| 4H | 1 | 200 | 4-9 h | 200 | 175.0 | -1.3 | -23.2 | 166.5 |
| 4I | 1 | 200 | 9-14 h | 200 | 150.0 | -2.3 | -23.2 | 167.1 |
| 4J | 1 | 200 | 14-24 h | 0 | 200.0 | -23.2 | -23.2 | 151.5 |
| 5A | 2 | 50 | \<= 14 h | 50 | 62.5 | 19.4 | -31.7 | 49.6 |
| 5B | 2 | 50 | 14-24 h | 0 | 87.5 | 15.9 | -31.7 | 46.6 |
| 5C | 2 | 100 | \<= 14 h | 100 | 125.0 | 19.4 | -31.7 | 99.3 |
| 5D | 2 | 100 | 14-24 h | 0 | 175.0 | 15.9 | -31.7 | 93.1 |
| 5E | 2 | 200 | \<= 2 h | 200 | 225.0 | 0.4 | -31.7 | 167.2 |
| 5F | 2 | 200 | \> 2 h | 0 | 300.0 | 0.0 | -31.7 | 166.9 |
| 6A | 3 | 50 | \<= 6 h | 50 | 75.0 | 24.4 | -34.8 | 51.0 |
| 6B | 3 | 50 | 6-24 h | 0 | 100.0 | 28.6 | -34.8 | 51.5 |
| 6C | 3 | 100 | \<= 6 h | 100 | 150.0 | 24.4 | -34.8 | 102.0 |
| 6D | 3 | 100 | 6-24 h | 0 | 200.0 | 28.6 | -34.8 | 103.1 |
| 6E | 3 | 200 | \<= 14 h | 200 | 250.0 | 16.3 | -34.8 | 192.9 |
| 6F | 3 | 200 | 14-24 h | 0 | 300.0 | -3.1 | -34.8 | 164.8 |

Every remedial recommendation of Xia 2025 Figures 4-6, simulated from
the packaged model. {.table}

``` r

# 1. Every recommendation that actually adjusts a dose brings the next trough
#    closer to steady state than skipping would.
adj <- dplyr::filter(remedialEval, adjusts)
stopifnot(nrow(adj) == 21L)
stopifnot(all(abs(adj$remPct) < abs(adj$noRemPct)))

# 2. Figure 4J is the one case where the paper prescribes no adjustment at all
#    (skip the missed 200 mg dose, take the regular dose): it must be identical
#    to the no-remedy arm.
j4 <- dplyr::filter(remedialEval, figure == "4J")
stopifnot(nrow(j4) == 1L, abs(j4$remPct - j4$noRemPct) < 1e-6)

# 3. No recommendation approaches the 300 ug/L laboratory alert level.
stopifnot(max(remedialEval$peak) < ALERT)

# 4. Overshoot is bounded: no recommendation leaves the trough more than 30%
#    above steady state.
stopifnot(max(remedialEval$remPct) < 30)

# 5. The total remedial drug escalates with the number of missed doses, at a
#    matched dose level and a matched position in the delay range (the paper's
#    central claim that the remedial dose depends on "the frequency of missed
#    doses"). Each (dose, nMiss) cell offers several delay windows, so the
#    earliest and the latest window of each cell are compared separately --
#    comparing across windows would mix the two effects the paper separates.
escalation <- remedialEval |>
  dplyr::mutate(total = immediate + nextDose) |>
  dplyr::group_by(dose, nMiss) |>
  dplyr::summarise(earliest = total[which.min(delay)],
                   latest   = total[which.max(delay)], .groups = "drop") |>
  dplyr::arrange(dose, nMiss)
stopifnot(nrow(escalation) == 9L)
for (d in unname(DOSES)) {
  e <- dplyr::filter(escalation, dose == d)
  stopifnot(nrow(e) == 3L, all(diff(e$earliest) > 0), all(diff(e$latest) >= 0))
}
knitr::kable(
  escalation |>
    dplyr::select(`Dose (mg)` = dose, `Missed doses` = nMiss,
                  `Total remedial dose, earliest window (mg)` = earliest,
                  `Total remedial dose, latest window (mg)` = latest),
  caption = "Total remedial drug (immediate + next scheduled) rises with the number of missed doses at every dose level."
)
```

| Dose (mg) | Missed doses | Total remedial dose, earliest window (mg) | Total remedial dose, latest window (mg) |
|---:|---:|---:|---:|
| 50 | 1 | 100.0 | 75.0 |
| 50 | 2 | 112.5 | 87.5 |
| 50 | 3 | 125.0 | 100.0 |
| 100 | 1 | 200.0 | 150.0 |
| 100 | 2 | 225.0 | 175.0 |
| 100 | 3 | 250.0 | 200.0 |
| 200 | 1 | 400.0 | 200.0 |
| 200 | 2 | 425.0 | 300.0 |
| 200 | 3 | 450.0 | 300.0 |

Total remedial drug (immediate + next scheduled) rises with the number
of missed doses at every dose level. {.table}

``` r


sprintf("max peak across all 22 recommendations: %.1f ug/L (alert level %d)",
        max(remedialEval$peak), ALERT)
#> [1] "max peak across all 22 recommendations: 192.9 ug/L (alert level 300)"
```

``` r

remSim |>
  dplyr::left_join(dplyr::select(remedialEval, scenario, figure, nMiss, dose, window),
                   by = "scenario") |>
  dplyr::filter(nMiss == 3L) |>
  dplyr::mutate(panel = sprintf("%s: %s missed x3, delay %s",
                                figure, names(DOSES)[match(dose, DOSES)], window)) |>
  ggplot(aes(time / DAY, Cc)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = LOWER_TR, linetype = "dashed", colour = "red") +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(x = "Day (first missed dose on day 3)", y = "Sertraline (ug/L)",
       title = "Remedial strategies for three consecutive missed doses",
       caption = "Replicates Figure 6 of Xia 2025. Red dashed line: 10 ug/L lower therapeutic limit.")
```

![](Xia_2025_sertraline_files/figure-html/figure-4-6-1.png)

### Where the recommendations sit relative to a trough-matching optimum

A natural question is whether each recommended dose is the one that best
restores the next morning’s trough. Holding the immediate dose as
recommended and sweeping the next scheduled dose over the 12.5 mg grid
the paper uses (a quarter of a 50 mg tablet) locates that optimum.

``` r

sweepGrid <- seq(0, 400, by = 12.5)
optimum <- vapply(seq_len(nrow(remedialEval)), function(i) {
  r  <- remedialEval[i, ]
  ev <- dplyr::bind_rows(lapply(seq_along(sweepGrid), function(k) {
    e <- buildScenario(r$nMiss, r$dose, r$delay, r$immediate, sweepGrid[k], id = k)
    e$scenario <- k
    e
  }))
  s <- solveTypical(ev)
  dev <- vapply(seq_along(sweepGrid), function(k) {
    abs(concAt(dplyr::filter(s, scenario == k), r$resumeDay * DAY + SAMPLE_H) - r$refTrough)
  }, numeric(1))
  sweepGrid[which.min(dev)]
}, numeric(1))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

optTab <- remedialEval |>
  dplyr::mutate(optimum = optimum, steps = (nextDose - optimum) / 12.5) |>
  dplyr::select(Fig = figure, Missed = nMiss, `Dose (mg)` = dose, Delay = window,
                `Xia 2025 (mg)` = nextDose, `Trough-matching optimum (mg)` = optimum,
                `Difference (12.5 mg steps)` = steps)
knitr::kable(optTab, digits = 1,
             caption = "Recommended next-scheduled dose against the dose that would exactly restore the following morning's trough.")
```

| Fig | Missed | Dose (mg) | Delay | Xia 2025 (mg) | Trough-matching optimum (mg) | Difference (12.5 mg steps) |
|:---|---:|---:|:---|---:|---:|---:|
| 4A | 1 | 50 | \<= 7 h | 50.0 | 50.0 | 0 |
| 4B | 1 | 50 | 7-14 h | 37.5 | 37.5 | 0 |
| 4C | 1 | 50 | \> 14 h | 75.0 | 62.5 | 1 |
| 4D | 1 | 100 | \<= 5 h | 100.0 | 100.0 | 0 |
| 4E | 1 | 100 | 5-14 h | 75.0 | 75.0 | 0 |
| 4F | 1 | 100 | \> 14 h | 150.0 | 137.5 | 1 |
| 4G | 1 | 200 | \<= 4 h | 200.0 | 187.5 | 1 |
| 4H | 1 | 200 | 4-9 h | 175.0 | 175.0 | 0 |
| 4I | 1 | 200 | 9-14 h | 150.0 | 162.5 | -1 |
| 4J | 1 | 200 | 14-24 h | 200.0 | 275.0 | -6 |
| 5A | 2 | 50 | \<= 14 h | 62.5 | 50.0 | 1 |
| 5B | 2 | 50 | 14-24 h | 87.5 | 75.0 | 1 |
| 5C | 2 | 100 | \<= 14 h | 125.0 | 100.0 | 2 |
| 5D | 2 | 100 | 14-24 h | 175.0 | 150.0 | 2 |
| 5E | 2 | 200 | \<= 2 h | 225.0 | 225.0 | 0 |
| 5F | 2 | 200 | \> 2 h | 300.0 | 300.0 | 0 |
| 6A | 3 | 50 | \<= 6 h | 75.0 | 50.0 | 2 |
| 6B | 3 | 50 | 6-24 h | 100.0 | 75.0 | 2 |
| 6C | 3 | 100 | \<= 6 h | 150.0 | 112.5 | 3 |
| 6D | 3 | 100 | 6-24 h | 200.0 | 150.0 | 4 |
| 6E | 3 | 200 | \<= 14 h | 250.0 | 200.0 | 4 |
| 6F | 3 | 200 | 14-24 h | 300.0 | 312.5 | -1 |

Recommended next-scheduled dose against the dose that would exactly
restore the following morning’s trough. {.table}

For a single missed dose the two agree to within one 12.5 mg step in 9
of 10 cases. For two and three missed doses the paper consistently
recommends *more* than the trough-matching dose, by up to four steps –
which is the expected consequence of its stated objective (“quickly
restore the sertraline steady-state concentration”): after two or three
missed doses the entire profile is depressed, not just the next trough,
so a dose that only repairs tomorrow’s trough leaves the following days
low.

The single large discrepancy runs the other way and is a safety
constraint, not a modelling one. In Figure 4J (one missed 200 mg dose,
14-24 h late) the trough-matching dose is 275 mg, but the paper
recommends the regular 200 mg, because it explicitly declines to
recommend any single remedial dose above 200 mg. The analysis here
recovers that boundary independently.

``` r

overCap <- dplyr::filter(remedialEval, nextDose > 200)
knitr::kable(
  dplyr::select(overCap, Fig = figure, Missed = nMiss, `Dose (mg)` = dose,
                Delay = window, `Next dose (mg)` = nextDose),
  caption = "Recommendations whose single remedial dose exceeds the 200 mg ceiling Xia 2025 sets for itself."
)
```

| Fig | Missed | Dose (mg) | Delay    | Next dose (mg) |
|:----|-------:|----------:|:---------|---------------:|
| 5E  |      2 |       200 | \<= 2 h  |            225 |
| 5F  |      2 |       200 | \> 2 h   |            300 |
| 6E  |      3 |       200 | \<= 14 h |            250 |
| 6F  |      3 |       200 | 14-24 h  |            300 |

Recommendations whose single remedial dose exceeds the 200 mg ceiling
Xia 2025 sets for itself. {.table}

``` r

stopifnot(all(overCap$dose == 200))
```

All four recommendations above the paper’s own 200 mg ceiling belong to
the 200 mg QD regimen, matching the paper’s summary that “when double or
triple doses of 200 mg were missed … the remedial dose at the next
schedule time was 225-300 mg” together with its instruction that such a
dose should not be taken. For the 200 mg regimen the paper’s practical
advice is therefore to skip and resume.

## Assumptions and deviations

- **Scale of the Table 2 variability rows (the one material
  interpretation).** Table 2 heads its variability column `IIV (%)` and
  its residual rows `PRO (%)` and `ADD (%)`, but prints
  fraction-magnitude numbers (0.134, 0.326, 0.0959, 24.7) in all four.
  They are encoded here as the raw NONMEM `$OMEGA` and `$SIGMA`
  **variances**, giving IIV of 37.9% CV on CL/F and 62.1% CV on V/F, a
  31.0% proportional residual and a 4.97 ug/L additive residual. Three
  checks support this over reading them as standard deviations: the
  bootstrap 95% CI on the proportional term reaches down to 0.01, which
  is a 10% CV as a variance but an impossible 1% CV as a standard
  deviation (the paper reports assay precision below 15% RSE); the
  additive term’s bootstrap upper bound of 344.63 read as a standard
  deviation would exceed the largest concentration in the dataset
  (299.87 ug/L); and with N = 103 the variance reading puts
  `omega/sqrt(N)` on CL/F at 3.6%, the same order as the 6% RSE Table 2
  reports, whereas the standard-deviation reading would imply 1.3%. The
  same group’s earlier paper (`modellib("Zhang_2024_sertraline")`)
  prints percentage-magnitude integers under an explicit `IIV (CV%)`
  heading and a raw variance for its proportional residual; Xia 2025
  prints untransformed NONMEM output throughout. **This interpretation
  affects no structural gate in this vignette** – every published
  quantity Xia 2025 reports was simulated at typical values – but it
  does set the width of the cohort band above.
- **`ADD (%)` units.** An additive residual on a concentration cannot be
  a percentage; the `(%)` on that row is a table-template artefact. It
  is encoded in ug/L, the unit of the observations.
- **Trough sampling time.** Xia 2025’s Methods state samples were drawn
  at 6-7 a.m. before the 08:30 dose, while the Results sentence
  describing the simulation says 7:00 a.m. The published troughs (22.34
  / 44.69 / 89.37 ug/L) are reproduced to 0.4% at 21.5 h post-dose
  (06:00) and to 4.5% at 22.5 h (07:00); exact reproduction requires
  21.41 h, i.e. 05:55. The 06:00 end of the stated window is used
  throughout. Nothing in the model depends on this choice; it only fixes
  where the published number is read off.
- **Recovery time.** The paper’s “at least five days” for recovery from
  a single missed dose is a visual read of Figure 3. The model gives 3
  days at a 5% criterion and 4 days at a 1% criterion; a *triple* missed
  dose takes 5 days at the 5% criterion. Reported above rather than
  gated tightly.
- **Delay values within each window.** The paper specifies delay
  *windows*, not single delays. A representative delay inside each
  window is simulated (for example 10 h for the “7-14 h” window).
  Because absorption is fast relative to elimination the results are
  insensitive to the choice within a window.
- **Missed-dose day.** Xia 2025 misses days 33-35 of a long course.
  Steady state is established here with a single `ss = 1` record and the
  first missed dose placed on day 3, which is verified above to agree
  with an explicit 33-day dosing history to 4e-07.
- **No covariates.** None were retained, so `covariateData` is empty.
  The eight covariates the paper screened (or, for CYP2C19, could not
  screen) are recorded in `covariatesDataExcluded` with the paper’s
  reasoning; they carry no encoded effect.
- **Supplement not on disk.** Xia 2025 cites Tables S1 (drug
  information) and S2 (theoretical remedial doses). Neither is in the
  PMC deposit (`hasSuppl: N`; the EuropePMC supplementary bundle
  contains only the figure and equation images). No parameter depends on
  them: every remedial dose used above is stated in the Results text,
  and Table S2’s content is reconstructed in the remedial table.
  Recorded as a documentation gap only.
- **Cohort simulation.** The 200-subject-per-arm cohort is not a
  replication of any published figure – Xia 2025 simulated only typical
  values – and its assertions are written on medians, never on extremes.
