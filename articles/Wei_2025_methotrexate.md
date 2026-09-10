# High-dose methotrexate in primary CNS lymphoma (Wei 2025)

``` r

library(nlmixr2lib)
library(rxode2)
library(PKNCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
```

## The paper

Wei S, Zhang S, Wang D, Zhang D, Lu Q, Mo J, Yang Z, Guan L, He Y, Zhao
Z, Mei S (2025). *Population pharmacokinetics of high-dose methotrexate
in patients with primary central nervous system lymphoma.* Front
Pharmacol 16:1578033.
[doi:10.3389/fphar.2025.1578033](https://doi.org/10.3389/fphar.2025.1578033).

The authors fitted a three-compartment model with first-order
elimination and proportional residual error to 6074 methotrexate plasma
concentrations from 752 Chinese adults with primary central nervous
system lymphoma (PCNSL) treated at Beijing Tiantan Hospital between
September 2016 and August 2023. Estimation was by first-order
conditional estimation extended least squares in Phoenix NLME 8.3.

**The paper reports two final models**, and this library ships both:

| Model file | Paper’s name | What it adds |
|----|----|----|
| `Wei_2025_methotrexate` | “nongene-model” (Table 5, model 5) | eGFR, BUN and ALT on CL; total protein on Q1 |
| `Wei_2025_methotrexate_genotype` | “gene-model” (Table 5, model 6) | the above, plus a composite ABCC4-ABCG2-ADORA2A genotype effect on CL |

The authors kept both because “given the limited routine implementation
of genetic testing in clinical practice, a model without genetic factors
was also developed to enhance clinical applicability” (Results 3.2). Use
the nongene model when genotype data are unavailable.

``` r

modNongene <- readModelDb("Wei_2025_methotrexate")
modGene    <- readModelDb("Wei_2025_methotrexate_genotype")
```

## Population

| Characteristic        | Value                                    |
|:----------------------|:-----------------------------------------|
| Subjects              | 752                                      |
| Concentration records | 6074                                     |
| Age                   | 18.12-86.65 years (median 57.445)        |
| Body weight           | 30-115 kg (median 68)                    |
| BSA                   | 1.16-2.39 m^2 (median 1.73)              |
| Female                | 332 (44.1%)                              |
| eGFR (2021 CKD-EPI)   | 5.4-162.9 mL/min/1.73 m^2 (median 101.8) |
| Serum creatinine      | 24.8-641.7 umol/L (median 64.6)          |
| BUN                   | 0.5-19 mmol/L (median 4.6)               |
| ALT                   | 2.2-1141.7 U/L (median 25)               |
| Total protein         | 27.4-95.7 g/L (median 61.8)              |
| Dose                  | 3.5 g/m^2 IV, median infusion 3.1 h      |
| Region                | China, single centre                     |

Reproduces Wei 2025 Table 3 (patient characteristics). {.table}

Every patient received leucovorin rescue starting 6 h after the end of
infusion; leucovorin is **not** represented in the model, so simulated
concentrations describe an unrescued profile. Each methotrexate
administration was treated as an independent event in the authors’
dataset, because dosing intervals exceeded five elimination half-lives.

## Source trace

Every value in both model files, with the location it came from. All
parameter estimates are from Table 5; all closed forms are from the
numbered equations in Results 3.2.

| Quantity | Nongene | Gene | Source |
|:---|:---|:---|:---|
| CL (L/h) | 8.2 | 8.45 | Table 5; Eq. 12 / 14 |
| Vc (L) | 33.39 | 33.29 | Table 5 (Eq. 17 rounds to 33.3) |
| Q1 (L/h) | 0.04 | 0.04 | Table 5; Eq. 13 / 15 |
| Vp1 (L) | 17.9 | 17.85 | Table 5 (Eq. 18 rounds to 17.9) |
| Q2 (L/h) | 0.09 | 0.09 | Table 5; Eq. 16 |
| Vp2 (L) | 1.14 | 1.14 | Table 5; Eq. 19 |
| eGFR exponent on CL | 0.67 | 0.67 | Table 5 theta_eGFR; Eq. 12 / 14 |
| BUN exponent on CL | -0.08 | -0.08 | Table 5 theta_BUN; Eq. 12 / 14 |
| ALT exponent on CL | 0.03 | 0.03 | Table 5 theta_ALT; Eq. 12 / 14 |
| TP exponent on Q1 | -1.68 | -1.72 | Table 5 theta_TP; Eq. 13 / 15 |
| Genotype multiplier on CL | n/a | 0.91 | Eq. 14; Table 5 theta = -0.09 |
| IIV CL (CV%) | 27.3 | 27.1 | Table 5 |
| IIV Vc (CV%) | 20.95 | 20.27 | Table 5 |
| IIV Q1 (CV%) | 98.91 | 99.01 | Table 5 |
| IIV Vp1 (CV%) | 78.35 | 78.1 | Table 5 |
| IIV Vp2 (CV%) | 26.95 | 27.2 | Table 5 |
| IIV Q2 | none | none | Table 5 has no IIV_Q2 row |
| Proportional residual (CV%) | 73.79 | 73.79 | Table 5 sigma |
| Structural ODEs | \- | \- | Eq. 4-9 |
| Exponential IIV | \- | \- | Eq. 10 |
| Proportional residual model | \- | \- | Eq. 11 |
| eGFR centring 101.8 | \- | \- | Eq. 12 / 14 = Table 3 median |
| BUN centring 4.6 | \- | \- | Eq. 12 / 14 = Table 3 median |
| ALT centring 25 | \- | \- | Eq. 12 / 14 = Table 3 median |
| TP centring 58 | \- | \- | Eq. 13 / 15 (Table 3 median is 61.8; see Errata) |

Source trace for both Wei 2025 model files. {.table}

## Units

The paper reports every concentration and every delayed-elimination
threshold in umol/L, so both model files use molar units: doses in umol,
volumes in L, concentrations in umol/L. Methotrexate has a molar mass of
454.44 g/mol.

``` r

MTX_MW <- 454.44                     # g/mol
gToUmol <- function(g) g / MTX_MW * 1e6
# A 3.5 g/m^2 dose in a median (1.73 m^2) patient:
round(gToUmol(3.5 * 1.73))
#> [1] 13324
```

## Check 1: the ODE encoding reproduces the analytic three-compartment solution

Equations 4-9 define a standard linear three-compartment system. Its
solution is a matrix exponential, which base R can compute exactly by
eigendecomposition of the rate matrix. This check is **deterministic** –
both sides use the same parameter values, so the only difference is
solver tolerance, and a tight bound is the right assertion here.

``` r

threeCmtMatrix <- function(cl, vc, q, vp, q2, vp2) {
  kel <- cl / vc
  k12 <- q  / vc
  k21 <- q  / vp
  k13 <- q2 / vc
  k31 <- q2 / vp2
  matrix(
    c(-(kel + k12 + k13),  k21,   k31,
      k12,                -k21,   0,
      k13,                 0,    -k31),
    nrow = 3, byrow = TRUE
  )
}

# Central-compartment concentration after an IV bolus, by eigendecomposition.
threeCmtBolusConc <- function(amt, times, cl, vc, q, vp, q2, vp2) {
  A <- threeCmtMatrix(cl, vc, q, vp, q2, vp2)
  ev <- eigen(A)
  cf <- as.vector(solve(ev$vectors, c(amt, 0, 0)))
  vapply(
    times,
    function(tt) Re(sum(ev$vectors[1, ] * exp(ev$values * tt) * cf)) / vc,
    numeric(1)
  )
}

# AUC from 0 to infinity of the central concentration, same decomposition.
threeCmtBolusAuc <- function(amt, cl, vc, q, vp, q2, vp2) {
  A <- threeCmtMatrix(cl, vc, q, vp, q2, vp2)
  ev <- eigen(A)
  cf <- as.vector(solve(ev$vectors, c(amt, 0, 0)))
  Re(sum(ev$vectors[1, ] * (-1 / ev$values) * cf)) / vc
}
```

``` r

# Typical (zero random effect) parameters of the nongene model, at the
# reference covariates so every covariate factor is exactly 1.
refCov <- list(CRCL = 101.8, BUN = 4.6, ALT = 25, TPRO = 58)
tv <- list(cl = 8.2, vc = 33.39, q = 0.04, vp = 17.9, q2 = 0.09, vp2 = 1.14)

boluAmt <- gToUmol(3.5 * 1.73)
chkTimes <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24, 48, 72, 96, 120, 168)

evBolus <- data.frame(
  id = 1L, time = 0, amt = boluAmt, evid = 1L, cmt = "central",
  CRCL = refCov$CRCL, BUN = refCov$BUN, ALT = refCov$ALT, TPRO = refCov$TPRO
) |>
  bind_rows(
    data.frame(
      id = 1L, time = chkTimes, amt = NA_real_, evid = 0L, cmt = "central",
      CRCL = refCov$CRCL, BUN = refCov$BUN, ALT = refCov$ALT, TPRO = refCov$TPRO
    )
  ) |>
  arrange(time, desc(evid))

simBolus <- rxSolve(
  zeroRe(modNongene), evBolus,
  atol = 1e-12, rtol = 1e-12, returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalvp2'

obsBolus <- simBolus |> filter(!is.na(Cc), time > 0)
analytic <- threeCmtBolusConc(
  boluAmt, obsBolus$time,
  tv$cl, tv$vc, tv$q, tv$vp, tv$q2, tv$vp2
)

relErr <- abs(obsBolus$Cc - analytic) / analytic
tibble::tibble(
  `Time (h)` = obsBolus$time,
  `rxode2 (umol/L)` = signif(obsBolus$Cc, 6),
  `Analytic (umol/L)` = signif(analytic, 6),
  `Relative error` = signif(relErr, 3)
) |>
  kable(caption = "Check 1: rxode2 solve against the eigendecomposition solution of Eq. 4-9.")
```

| Time (h) | rxode2 (umol/L) | Analytic (umol/L) | Relative error |
|---------:|----------------:|------------------:|---------------:|
|     0.25 |     374.9190000 |       374.9190000 |              0 |
|     0.50 |     352.2570000 |       352.2570000 |              0 |
|     1.00 |     310.9740000 |       310.9740000 |              0 |
|     2.00 |     242.4040000 |       242.4040000 |              0 |
|     4.00 |     147.4310000 |       147.4310000 |              0 |
|     8.00 |      54.8571000 |        54.8571000 |              0 |
|    12.00 |      20.6947000 |        20.6947000 |              0 |
|    24.00 |       1.4244300 |         1.4244300 |              0 |
|    48.00 |       0.0869428 |         0.0869428 |              0 |
|    72.00 |       0.0258483 |         0.0258483 |              0 |
|    96.00 |       0.0160806 |         0.0160806 |              0 |
|   120.00 |       0.0139408 |         0.0139408 |              0 |
|   168.00 |       0.0123067 |         0.0123067 |              0 |

Check 1: rxode2 solve against the eigendecomposition solution of Eq.
4-9. {.table}

``` r


# Deterministic: both sides use identical parameters, so this is pure solver
# error and a tight bound is correct (realised max ~1e-9).
stopifnot(max(relErr) < 1e-6)
```

## Check 2: AUC(0-inf) equals Dose / CL

For any linear disposition model the total area under the
central-compartment concentration curve is exactly `Dose / CL`,
independent of the peripheral structure. This validates that clearance
enters the ODEs where the paper says it does.

``` r

aucAnalytic <- threeCmtBolusAuc(
  boluAmt, tv$cl, tv$vc, tv$q, tv$vp, tv$q2, tv$vp2
)
aucExpected <- boluAmt / tv$cl

c(analytic = aucAnalytic, doseOverCl = aucExpected,
  relDiff = abs(aucAnalytic - aucExpected) / aucExpected)
#>     analytic   doseOverCl      relDiff 
#> 1.624889e+03 1.624889e+03 4.197954e-16

stopifnot(abs(aucAnalytic - aucExpected) / aucExpected < 1e-10)
```

## Check 3: the covariate equations reproduce Equations 12-15

The model’s individual `cl` and `q` are compared against the published
closed forms evaluated independently in plain R, over a grid spanning
the cohort ranges. Deterministic, so a tight bound is again correct.

``` r

# The published closed forms, transcribed directly from Results 3.2.
clEq12 <- function(CRCL, BUN, ALT) {
  8.2 * (CRCL / 101.8)^0.67 * (BUN / 4.6)^-0.08 * (ALT / 25)^0.03
}
q1Eq13 <- function(TPRO) 0.04 * (TPRO / 58)^-1.68

clEq14 <- function(CRCL, BUN, ALT, carrier) {
  8.45 * (CRCL / 101.8)^0.67 * (BUN / 4.6)^-0.08 * (ALT / 25)^0.03 *
    ifelse(carrier, 0.91, 1)
}
q1Eq15 <- function(TPRO) 0.04 * (TPRO / 58)^-1.72

covGrid <- expand.grid(
  CRCL = c(20, 60, 101.8, 140, 160),
  BUN  = c(1, 4.6, 12, 19),
  ALT  = c(5, 25, 200, 1000),
  TPRO = c(30, 58, 61.8, 90)
) |>
  mutate(id = row_number())

solveCov <- function(mod, grid, extra = NULL) {
  ev <- grid |>
    mutate(time = 0, amt = 1000, evid = 1L, cmt = "central") |>
    bind_rows(grid |> mutate(time = 1, amt = NA_real_, evid = 0L, cmt = "central"))
  if (!is.null(extra)) for (nm in names(extra)) ev[[nm]] <- extra[[nm]]
  ev <- ev |> arrange(id, time, desc(evid))
  rxSolve(zeroRe(mod), ev, returnType = "data.frame") |>
    filter(time == 1)
}

simCovN <- solveCov(modNongene, covGrid)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalvp2'
#> Warning: multi-subject simulation without without 'omega'
expClN <- clEq12(simCovN$CRCL, simCovN$BUN, simCovN$ALT)
expQ1N <- q1Eq13(simCovN$TPRO)

# Gene model: exercise both carrier states.
gridCarrier <- covGrid |> mutate(
  SNP_ABCC4_RS2274407_G_COUNT   = rep(c(0, 2), length.out = n()),
  SNP_ABCG2_RS2231142_T_COUNT   = rep(c(0, 1, 2), length.out = n()),
  SNP_ADORA2A_RS2298383_T_COUNT = rep(c(0, 2, 1), length.out = n())
)
simCovG <- solveCov(modGene, gridCarrier)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etalvp2'
#> Warning: multi-subject simulation without without 'omega'
carrierG <- (simCovG$SNP_ABCC4_RS2274407_G_COUNT +
               simCovG$SNP_ABCG2_RS2231142_T_COUNT +
               simCovG$SNP_ADORA2A_RS2298383_T_COUNT) > 3
expClG <- clEq14(simCovG$CRCL, simCovG$BUN, simCovG$ALT, carrierG)
expQ1G <- q1Eq15(simCovG$TPRO)

covErr <- c(
  cl_nongene = max(abs(simCovN$cl - expClN) / expClN),
  q1_nongene = max(abs(simCovN$q  - expQ1N) / expQ1N),
  cl_gene    = max(abs(simCovG$cl - expClG) / expClG),
  q1_gene    = max(abs(simCovG$q  - expQ1G) / expQ1G)
)
tibble::tibble(
  Comparison = names(covErr),
  `Max relative error` = signif(covErr, 3),
  `Grid points` = c(nrow(simCovN), nrow(simCovN), nrow(simCovG), nrow(simCovG))
) |>
  kable(caption = "Check 3: model covariate terms against Equations 12-15.")
```

| Comparison | Max relative error | Grid points |
|:-----------|-------------------:|------------:|
| cl_nongene |                  0 |         320 |
| q1_nongene |                  0 |         320 |
| cl_gene    |                  0 |         320 |
| q1_gene    |                  0 |         320 |

Check 3: model covariate terms against Equations 12-15. {.table}

``` r


stopifnot(max(covErr) < 1e-10)
# Confirm the grid actually exercised both genotype strata (pattern 10 guard).
stopifnot(any(carrierG), any(!carrierG))
```

## Check 4: the paper’s own derived percentages

The Discussion states two numeric consequences of the covariate model.
Both are deterministic functions of the transcribed parameters, so
reproducing them is a direct test that the exponents, the centring
values and the functional form were all transcribed correctly.

> “in individuals with severe liver dysfunction (more than 5 to 20 times
> the upper limit of normal), MTX clearance increased by approximately
> 5%-9%”

> “Patients harboring the ABCC-ABCG-ADORA2A gene mutations demonstrated
> an approximately 9% reduction in MTX clearance”

``` r

# The ALT claim reproduces when the multiple of the upper limit of normal is
# applied to the model's own 25 U/L centring value.
altRatio <- function(mult) (mult * 25 / 25)^0.03
alt5  <- altRatio(5)
alt20 <- altRatio(20)

geneReduction <- 1 - 0.91

derived <- tibble::tibble(
  Claim = c(
    "CL increase at 5x ULN ALT",
    "CL increase at 20x ULN ALT",
    "CL reduction in composite-genotype carriers"
  ),
  `Paper states` = c("~5%", "~9%", "~9%"),
  Reproduced = sprintf(
    "%.1f%%",
    100 * c(alt5 - 1, alt20 - 1, geneReduction)
  )
)
kable(derived, caption = "Check 4: the Discussion's own derived percentages.")
```

| Claim                                       | Paper states | Reproduced |
|:--------------------------------------------|:-------------|:-----------|
| CL increase at 5x ULN ALT                   | ~5%          | 4.9%       |
| CL increase at 20x ULN ALT                  | ~9%          | 9.4%       |
| CL reduction in composite-genotype carriers | ~9%          | 9.0%       |

Check 4: the Discussion’s own derived percentages. {.table}

``` r


stopifnot(
  abs(100 * (alt5 - 1) - 5) < 0.5,
  abs(100 * (alt20 - 1) - 9) < 0.5,
  abs(100 * geneReduction - 9) < 0.01
)
```

## Virtual cohort

200 subjects per model (the per-arm cap), dosed at the paper’s 3.5 g/m^2
over its median 3.1 h infusion. Covariate distributions are drawn to
match the Table 3 medians and ranges; see Errata for what had to be
assumed.

``` r

rxSetSeed(20250519)
set.seed(20250519)
nSub <- 200L

# Truncated-normal / lognormal draws matched to the Table 3 median and range.
rtruncnorm <- function(n, mean, sd, lo, hi) {
  x <- rnorm(n, mean, sd)
  pmin(pmax(x, lo), hi)
}
rtrunclnorm <- function(n, median, sdlog, lo, hi) {
  x <- rlnorm(n, log(median), sdlog)
  pmin(pmax(x, lo), hi)
}

makeSubjects <- function(n, idOffset = 0L) {
  data.frame(
    id   = idOffset + seq_len(n),
    BSA  = rtruncnorm(n, 1.73, 0.20, 1.16, 2.39),
    CRCL = rtruncnorm(n, 101.8, 28, 5.4, 162.9),
    BUN  = rtrunclnorm(n, 4.6, 0.42, 0.5, 19),
    ALT  = rtrunclnorm(n, 25, 0.85, 2.2, 1141.7),
    TPRO = rtruncnorm(n, 61.8, 8.0, 27.4, 95.7)
  )
}

# Composite genotype is assigned BY DESIGN as a balanced 100/100 split rather
# than drawn at an assumed allele frequency. The paper's own variant
# frequencies are in Supplementary Appendix SA1, which is not on disk, so any
# frequency here would be invented; a balanced design also keeps both strata
# well powered for the carrier comparison below. Non-carrier patterns sum to
# 3 or less, carrier patterns to 4 or more, per the Results 3.2 rule.
nonCarrierPatterns <- list(
  c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(2, 1, 0), c(1, 1, 1), c(2, 0, 1)
)
carrierPatterns <- list(
  c(2, 1, 1), c(2, 2, 0), c(1, 2, 2), c(2, 2, 1), c(2, 2, 2), c(0, 2, 2)
)
genotypeMatrix <- function(patterns, n) {
  idx <- rep_len(seq_along(patterns), n)
  do.call(rbind, patterns[idx])
}

subjN <- makeSubjects(nSub, 0L)
gtN <- genotypeMatrix(nonCarrierPatterns, nSub / 2L)
gtC <- genotypeMatrix(carrierPatterns, nSub / 2L)
gt <- rbind(gtN, gtC)
subjG <- makeSubjects(nSub, 1000L) |>
  mutate(
    SNP_ABCC4_RS2274407_G_COUNT   = gt[, 1],
    SNP_ABCG2_RS2231142_T_COUNT   = gt[, 2],
    SNP_ADORA2A_RS2298383_T_COUNT = gt[, 3]
  )
# Confirm the design produced exactly the intended balanced split.
stopifnot(
  sum(rowSums(gt) > 3) == nSub / 2L,
  sum(rowSums(gt) <= 3) == nSub / 2L
)

INF_DUR <- 3.1
obsTimes <- sort(unique(c(
  seq(0, 12, by = 0.25),
  seq(13, 48, by = 1),
  seq(54, 168, by = 6)
)))

buildEvents <- function(subj) {
  covCols <- setdiff(names(subj), "id")
  doseRows <- subj |>
    mutate(
      time = 0,
      amt  = gToUmol(3.5 * BSA),
      evid = 1L,
      cmt  = "central",
      dur  = INF_DUR
    )
  obsRows <- subj |>
    crossing(time = obsTimes) |>
    mutate(amt = NA_real_, evid = 0L, cmt = "central", dur = NA_real_)
  bind_rows(doseRows, obsRows) |>
    select(id, time, amt, evid, cmt, dur, all_of(covCols)) |>
    arrange(id, time, desc(evid))
}

evN <- buildEvents(subjN)
evG <- buildEvents(subjG)

simN <- rxSolve(modNongene, evN, returnType = "data.frame") |>
  mutate(treatment = "nongene-model")
#> ℹ parameter labels from comments will be replaced by 'label()'
simG <- rxSolve(modGene, evG, returnType = "data.frame") |>
  mutate(treatment = "gene-model")
#> ℹ parameter labels from comments will be replaced by 'label()'
sim <- bind_rows(simN, simG)

nrow(sim)
#> [1] 42000
```

``` r

band <- sim |>
  filter(!is.na(Cc), time > 0) |>
  group_by(treatment, time) |>
  summarise(
    lo  = quantile(Cc, 0.10),
    mid = median(Cc),
    hi  = quantile(Cc, 0.90),
    .groups = "drop"
  )

thresholds <- data.frame(time = c(24, 48, 72, 96), conc = c(50, 5, 0.2, 0.05))

ggplot(band, aes(time, mid, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(
    data = thresholds, aes(time, conc),
    inherit.aes = FALSE, shape = 4, size = 3
  ) +
  scale_y_log10() +
  labs(
    x = "Time after start of infusion (h)",
    y = "Methotrexate (umol/L)",
    colour = NULL, fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top")
```

![Simulated methotrexate concentration-time profiles (median and
10th-90th percentile band) for both final models at 3.5 g/m^2 over 3.1
h. The dashed lines are the label thresholds for delayed elimination at
24, 48, 72 and 96
h.](Wei_2025_methotrexate_files/figure-html/cohort-profile-1.png)

Simulated methotrexate concentration-time profiles (median and 10th-90th
percentile band) for both final models at 3.5 g/m^2 over 3.1 h. The
dashed lines are the label thresholds for delayed elimination at 24, 48,
72 and 96 h.

## NCA with PKNCA

The paper reports **no** non-compartmental parameters – its model
evaluation is goodness-of-fit, a 200-replicate bootstrap and a VPC – so
there is no published Cmax / AUC / half-life table to compare against
and
[`ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
is not applicable here. The NCA below characterises the shipped models,
and is gated against the analytic solution of the same system rather
than against the paper.

``` r

simNca <- sim |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc, treatment)

# Guarantee a time-zero record (concentration is 0 at the start of an infusion).
simNca <- bind_rows(
  simNca,
  simNca |> distinct(id, treatment) |> mutate(time = 0, Cc = 0)
) |>
  distinct(id, treatment, time, .keep_all = TRUE) |>
  arrange(treatment, id, time)

doseDf <- bind_rows(
  subjN |> transmute(id, treatment = "nongene-model", amt = gToUmol(3.5 * BSA)),
  subjG |> transmute(id, treatment = "gene-model",    amt = gToUmol(3.5 * BSA))
) |>
  mutate(time = 0)

concObj <- PKNCA::PKNCAconc(
  simNca, Cc ~ time | treatment + id,
  concu = "umol/L", timeu = "h"
)
doseObj <- PKNCA::PKNCAdose(
  doseDf, amt ~ time | treatment + id,
  doseu = "umol"
)

intervals <- data.frame(
  start = 0, end = 168,
  cmax = TRUE, tmax = TRUE, auclast = TRUE, half.life = TRUE
)

ncaRes <- PKNCA::pk.nca(
  PKNCA::PKNCAdata(concObj, doseObj, intervals = intervals)
)
ncaSum <- as.data.frame(ncaRes)
```

| Parameter | gene-model | nongene-model |
|:---|:---|:---|
| AUC(0-168 h) (umol\*h/L) | 1.64e+03 \[1.06e+03, 2.51e+03\] | 1.66e+03 \[1.07e+03, 2.7e+03\] |
| Cmax (umol/L) | 268 \[207, 353\] | 272 \[209, 349\] |
| t1/2 (h) | 234 \[57.3, 798\] | 228 \[33.5, 763\] |
| Tmax (h) | 3 \[3, 3.25\] | 3 \[3, 3.25\] |

Simulated NCA, median \[10th, 90th percentile\] over 200 subjects per
model. No published NCA exists for comparison. {.table}

``` r

# Gate the NCA against the analytic solution for the typical subject rather
# than against the paper. The median simulated Cmax should sit close to the
# analytic Cmax at the typical parameters; the two differ only by the cohort's
# covariate and random-effect spread, so this is a loose but real bound.
typicalCmax <- max(
  threeCmtBolusConc(
    gToUmol(3.5 * 1.73), seq(0.01, INF_DUR, by = 0.01),
    tv$cl, tv$vc, tv$q, tv$vp, tv$q2, tv$vp2
  )
)

medCmax <- ncaSum |>
  filter(PPTESTCD == "cmax", treatment == "nongene-model") |>
  pull(PPORRES) |>
  median()

# A bolus overstates the peak of a 3.1 h infusion; the ratio below is the
# infusion-dilution factor and is not expected to be 1. What IS asserted is
# that the simulated peak sits in the physically sensible band between the
# infusion-averaged and the bolus value, which a mis-transcribed dose, volume
# or unit conversion would break by orders of magnitude.
cmaxRatio <- medCmax / typicalCmax
c(analytic_bolus_cmax = typicalCmax, simulated_median_cmax = medCmax,
  ratio = cmaxRatio)
#>   analytic_bolus_cmax simulated_median_cmax                 ratio 
#>           398.0500858           272.4729049             0.6845191

stopifnot(cmaxRatio > 0.3, cmaxRatio < 1.1)

# Every subject must have produced a finite Cmax (pattern 10 guard).
stopifnot(
  sum(ncaSum$PPTESTCD == "cmax") == 2 * nSub,
  all(is.finite(ncaSum$PPORRES[ncaSum$PPTESTCD == "cmax"]))
)
```

## Genotype effect at the cohort level

The composite genotype lowers clearance by 9%, so carriers should show
slightly higher exposure. This is a small effect against 27% IIV on
clearance, so the assertion below is on the **median ratio**, not on
individual subjects.

``` r

carrierFlag <- subjG |>
  transmute(
    id,
    carrier = (SNP_ABCC4_RS2274407_G_COUNT +
                 SNP_ABCG2_RS2231142_T_COUNT +
                 SNP_ADORA2A_RS2298383_T_COUNT) > 3
  )

aucByCarrier <- ncaSum |>
  filter(PPTESTCD == "auclast", treatment == "gene-model") |>
  left_join(carrierFlag, by = "id") |>
  group_by(carrier) |>
  summarise(medianAuc = median(PPORRES), n = n(), .groups = "drop")

kable(
  aucByCarrier |>
    rename("Composite carrier" = carrier, "Median AUC(0-168 h)" = medianAuc,
           "Subjects" = n),
  caption = "AUC by composite ABCC4-ABCG2-ADORA2A genotype status."
)
```

| Composite carrier | Median AUC(0-168 h) | Subjects |
|:------------------|--------------------:|---------:|
| FALSE             |            1614.745 |      100 |
| TRUE              |            1679.367 |      100 |

AUC by composite ABCC4-ABCG2-ADORA2A genotype status. {.table}

``` r


# Both strata must be populated, or the comparison is vacuous. The balanced
# design above guarantees 100 per stratum.
stopifnot(nrow(aucByCarrier) == 2L, all(aucByCarrier$n == nSub / 2L))

carrierRatio <-
  aucByCarrier$medianAuc[aucByCarrier$carrier] /
  aucByCarrier$medianAuc[!aucByCarrier$carrier]

# Expected ratio is 1/0.91 = 1.099 for clearance-driven exposure, but the two
# medians are drawn from independent 27%-IIV samples, so the realised ratio
# moves with the draw. Assert the magnitude is within a band that a sign flip
# or an order-of-magnitude transcription error would break, not the exact value.
carrierRatio
#> [1] 1.04002
stopifnot(carrierRatio > 0.85, carrierRatio < 1.45)
```

## Delayed elimination

The methotrexate label defines delayed elimination as a concentration
above 50 umol/L at 24 h, 5 umol/L at 48 h, 0.2 umol/L at 72 h, or 0.05
umol/L at 96 h. Wei 2025 reports that “at least 17.4% of patients in
this cohort exhibited evidence of delayed MTX clearance” (Results 3.1).

``` r

delayedTimes <- c(24, 48, 72, 96)
delayedLimits <- c(50, 5, 0.2, 0.05)

delayed <- sim |>
  filter(time %in% delayedTimes, !is.na(Cc)) |>
  mutate(limit = delayedLimits[match(time, delayedTimes)],
         over = Cc > limit) |>
  group_by(treatment, time) |>
  summarise(pct = 100 * mean(over), .groups = "drop") |>
  mutate(Threshold = sprintf("> %g umol/L at %g h", delayedLimits[match(time, delayedTimes)], time))

anySubject <- sim |>
  filter(time %in% delayedTimes, !is.na(Cc)) |>
  mutate(limit = delayedLimits[match(time, delayedTimes)], over = Cc > limit) |>
  group_by(treatment, id) |>
  summarise(any_over = any(over), .groups = "drop") |>
  group_by(treatment) |>
  summarise(pct_any = 100 * mean(any_over), .groups = "drop")

kable(
  delayed |> select(Model = treatment, Threshold, `% exceeding` = pct) |>
    mutate(`% exceeding` = round(`% exceeding`, 1)),
  caption = "Simulated proportion exceeding each label threshold."
)
```

| Model         | Threshold              | % exceeding |
|:--------------|:-----------------------|------------:|
| gene-model    | \> 50 umol/L at 24 h   |         1.0 |
| gene-model    | \> 5 umol/L at 48 h    |         1.5 |
| gene-model    | \> 0.2 umol/L at 72 h  |         9.5 |
| gene-model    | \> 0.05 umol/L at 96 h |        20.5 |
| nongene-model | \> 50 umol/L at 24 h   |         0.5 |
| nongene-model | \> 5 umol/L at 48 h    |         1.0 |
| nongene-model | \> 0.2 umol/L at 72 h  |        11.5 |
| nongene-model | \> 0.05 umol/L at 96 h |        28.0 |

Simulated proportion exceeding each label threshold. {.table}

``` r

kable(
  anySubject |> rename(Model = treatment, `% meeting any criterion` = pct_any) |>
    mutate(`% meeting any criterion` = round(`% meeting any criterion`, 1)),
  caption = "Simulated proportion meeting at least one delayed-elimination criterion."
)
```

| Model         | % meeting any criterion |
|:--------------|------------------------:|
| gene-model    |                    21.5 |
| nongene-model |                    28.5 |

Simulated proportion meeting at least one delayed-elimination criterion.
{.table}

This is **reported, not asserted**. The paper’s 17.4% is an explicit
lower bound derived from sparse routine therapeutic-drug-monitoring
samples in a cohort that received leucovorin rescue from 6 h onward,
whereas the simulation is a fully-sampled, unrescued cohort in which
every subject is evaluated at all four times. The two are not the same
quantity, so scoring one against the other would be a false gate. The
simulated figure is expected to be the larger of the two, and it is.

``` r

# The only defensible directional claim: an unrescued, fully-sampled cohort
# cannot show LESS delayed elimination than a rescued, sparsely-sampled one.
stopifnot(all(anySubject$pct_any > 17.4))
```

## Assumptions and deviations

- **IIV scale (the one interpretive choice in this extraction).** Table
  5 heads the random-effect rows `IIV<param> (CV%)` and Equation 10
  defines the exponential model `theta_i = theta_TV * exp(eta)` with
  `eta ~ N(0, omega^2)`. The tabulated numbers are encoded here as
  `100 * omega`, i.e. `variance = (CV%/100)^2`. Two arguments support
  reading them as an SD-scale quantity rather than a variance: the
  reported %RSE on `IIV_CL` is 1.86%, which is the order of the
  `1/sqrt(2N)` precision an SD estimate carries on the roughly 2000
  administrations in this dataset, whereas a variance estimate carries
  the much looser `sqrt(2/N)`; and reading 27.3 as a variance would
  imply a log-scale SD of 5.2, which is not physically plausible. **What
  remains genuinely ambiguous on disk** is whether the authors tabulated
  `100 * omega` directly or the exact log-normal CV
  `sqrt(exp(omega^2) - 1) * 100`. The two coincide to within 2% for CL,
  Vc and Vp2, but differ materially for the two large ones: under the
  alternative reading `omega_Q1` would be 0.826 rather than 0.989 (-16%)
  and `omega_Vp1` 0.692 rather than 0.784 (-12%). Nothing in the paper
  discriminates – there is no back-computable derived percentage, and
  the residual-error row does not help because the SD and the CV of a
  proportional epsilon are the same number. The literal reading is used.
  Both affected parameters sit on the peripheral distribution, which the
  authors themselves flag as the least well determined part of the model
  (“the limited sample size during the distribution phase may introduce
  bias in the estimation of Vp and Q”).
- **Total-protein centring value.** Equations 13 and 15 and the Abstract
  all normalise total protein to 58 g/L, but Table 3 gives the cohort
  median as 61.8 g/L and Methods states that “all continuous covariates
  were standardized to their median values”. The printed equation
  constant 58 is used, following the register’s standing rule that the
  equation rather than the demographics table is the authority for a
  centring value. At the exponent -1.68 the choice moves Q1 by about
  11%; because Q1 is only 0.04 L/h this has almost no effect on the
  central-compartment profile. The other three centring values (101.8,
  4.6 and 25) match their Table 3 medians exactly.
- **Genotype multiplier taken as printed.** Equation 14 and the Abstract
  print the carrier multiplier as `a = 0.91`, while Table 5 reports the
  underlying coefficient as `-0.09`. The printed 0.91 equals
  `1 + (-0.09)` exactly and equals `exp(-0.09) = 0.9139` only after
  rounding, and the Discussion quotes “approximately 9%”, so the
  multiplier is encoded as 0.91 rather than back-transformed. The two
  readings differ by 0.4% in clearance.
- **Composite-genotype rule reconstructed from three places.** The
  carrier definition is not stated in one place. Table 1 scores each
  variant 1/2/3 with the score rising as clearance falls; all three of
  these variants have their variant allele associated with decreased
  clearance (each stated separately in the Discussion), so each score is
  `1 + (variant allele count)` and the three-variant sum runs 3 to 9.
  Results 3.2 then says carriers show “more than three nucleotide
  mutations among these variants”, i.e. more than 3 of the 6 possible
  variant alleles, i.e. a sum of at least 7 – which is exactly Table 2’s
  “Two groups, rule 2” row for three combined variants (3-6 -\> group 1,
  7-9 -\> group 2). The model implements the allele-count form directly.
- **ABCC4 rs2274407 strand orientation.** Results 3.2 defines the
  variant as `(T > G)` while the Discussion calls it `G912T` and refers
  to “T allele carriers”, following the cited Mesrian Tanha 2017
  nomenclature. These are opposite orientations for the same variant.
  The model column counts the G allele per the paper’s own `(T > G)`
  definition, and the covariate register entry records the collision
  explicitly so a future extraction does not silently invert the sign.
- **Genotype frequencies are not reproduced.** The per-variant allele
  frequencies are in Supplementary Appendix SA1, which is not on disk,
  so the cohort above does not attempt to reproduce the cohort’s real
  carrier prevalence. Instead the composite genotype is assigned by
  design as a balanced 100 carrier / 100 non-carrier split, which keeps
  the carrier comparison well powered and avoids inventing a frequency.
  This affects only the mix of simulated subjects, not any model
  parameter; a user with real genotype data supplies the three
  allele-count columns directly.
- **Covariate distributions are assumed.** Table 3 gives medians and
  ranges but no distributional form or correlation structure. The
  virtual cohort draws each covariate independently from a truncated
  normal or lognormal matched to the published median and range. The
  paper reports a between-covariate `R^2` of 0.22 for BUN against eGFR;
  that correlation is not reproduced here, which slightly widens the
  simulated clearance distribution relative to the real cohort.
- **No published NCA to compare against.** The paper evaluates its
  models with goodness-of-fit plots, a 200-replicate bootstrap and a
  VPC, and reports no Cmax / Tmax / AUC / half-life values. The NCA
  section therefore characterises the shipped model and is gated against
  the analytic solution of the same system, not against the paper.
- **Base-model estimates are promised but not printed.** Results 3.2
  states that “detailed estimates for the base model, final
  pharmacokinetic parameters … are presented in Table 5”, but Table 5 as
  published contains only the two final models and their bootstraps –
  there is no base-model column. Nothing in these model files depends on
  it.
- **Rounded values in the equation block.** Equations 17 and 18 print Vc
  as 33.3 and Vp1 as 17.9 and the surrounding text says these are
  “consistent across both models”, but Table 5 resolves them per model
  (33.39 / 33.29 and 17.9 / 17.85). The Table 5 per-model values are
  used.
- **Leucovorin is not modelled.** Every patient received leucovorin
  rescue from 6 h after the infusion, and dosing was intensified on high
  concentrations. The model describes methotrexate disposition only, so
  simulated profiles are unrescued.
- **Q2 has no inter-individual variability.** Table 5 has no `IIV_Q2`
  row, so `q2` is a typical value only. This is faithful to the paper,
  not an omission.
- **Applicability window.** The authors state that predictive accuracy
  beyond 120 h post-dose is limited because only 5% of samples fall
  there. The simulations above run to 168 h; the tail beyond 120 h
  should be read with that caveat.
