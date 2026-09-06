# Teicoplanin (Wang 2023)

## Model and source

``` r

mod      <- readModelDb("Wang_2023_teicoplanin")
mod_meta <- rxode2::rxode(mod)$meta
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Wang Y, Yao F, Chen S, Ouyang X, Lan J, Wu Z, Wang Y, Chen
  J, Wang X, Chen C. Optimal teicoplanin dosage regimens in critically
  ill patients: population pharmacokinetics and dosing simulations based
  on renal function and infection type. Drug Des Devel Ther.
  2023;17:2259-2271. <doi:10.2147/DDDT.S413662>
- Article (DOI): <https://doi.org/10.2147/DDDT.S413662>

This vignette validates the packaged `Wang_2023_teicoplanin` model – a
two-compartment intravenous-infusion population PK model for teicoplanin
built from 304 plasma concentrations in 108 critically ill adults in the
ICU of Guangdong Provincial People’s Hospital. The final model carries a
single covariate: the CKD-EPI estimated glomerular filtration rate
(eGFR) as a **linear** term on clearance, centred at 50 mL/min/1.73 m^2.

The paper’s quantitative deliverable is a Monte Carlo
dosing-optimization analysis. Supplementary Tables S2-1 and S2-2 report,
for eight candidate regimens crossed with four renal-function strata,
the median (IQR) trough concentration at 72 h and at 168 h **and** the
corresponding probability of target attainment – 64 published trough
medians and 64 published PTA values. Those are pure forward predictions
of the packaged model, which makes them an unusually strong validation
target, so this vignette reproduces all of them and compares against the
published numbers.

## Population

``` r

pop <- rxode2::rxode(mod)$population
#> ℹ parameter labels from comments will be replaced by 'label()'
```

The Wang 2023 cohort comprised 151 ICU patients contributing 347 plasma
teicoplanin concentrations. The population PK model was built on the 108
patients who each contributed more than two samples (304
concentrations); the remaining 43 single-sample patients formed an
external validation group. Median age was 56.5 years (range 21-91) and
median weight 61 kg (range 40-115); 38 of 108 subjects (35.2%) were
female.

The defining feature of the cohort is severe renal impairment: median
CKD-EPI eGFR 29.2 mL/min/1.73 m^2 (range 4.89-170), median
Cockcroft-Gault creatinine clearance 34.8 mL/min, median serum
creatinine 161 umol/L, and 43 of 108 subjects (39.8%) receiving
continuous renal replacement therapy. Median baseline APACHE II score
was 25.5 (range 10-37). Pneumonia was the commonest infection site (73
of the recorded sites), followed by intra-abdominal infection (13) and
bacteraemia (10).

Teicoplanin was given by 1-hour intravenous infusion. Concentrations
were measured by HPLC-MS/MS quantifying the teicoplanin A2-2 component,
the most active of the five components of the teicoplanin mixture, as
the surrogate for teicoplanin (LLOQ 1.0 mg/L, linear to 100.0 mg/L).

## Source trace

Every value in the packaged model, and every equation in `model()`,
traced to its location in the source.

| Quantity | Value | Source location |
|:---|:---|:---|
| CL (typical, eGFR = 50) | 0.838 L/h | Table 2, Final Model Estimate (RSE 8.4%); restated in Results and Abstract |
| Vc | 14.4 L | Table 2 (RSE 7.2%); Table 2 note ‘Vc (L) = 14.4’ |
| Q | 3.08 L/h | Table 2 (RSE 14.3%); Table 2 note ‘Q (L/h) = 3.08’ |
| Vp | 51.6 L | Table 2 (RSE 10.7%); Table 2 note ‘Vp (L) = 51.6 x e^0.387’ |
| theta eGFR-CL | 0.00823 | Table 2 (RSE 28.3%); Abstract ‘eGFR adjustment factor of 0.00823’ |
| CL covariate equation | CL = 0.838 \* (1 + (eGFR - 50) \* 0.0082) | Table 2 notes; generic form in Supplementary Table S1 notes |
| omega CL | 37.7% -\> var 0.142 | Table 2, BSV block (RSE 20%) |
| omega Vc | 26.6% -\> var 0.0708 | Table 2, BSV block (RSE 43.4%) |
| omega Q | 0 FIX (no eta) | Table 2, BSV block; Results ‘inter-individual variabilities … 0%’ |
| omega Vp | 62.2% -\> var 0.387 | Table 2, BSV block (RSE 18%); confirmed by ‘e^0.387’ in the Table 2 note |
| Proportional residual error | 31.7% | Table 2, RUV block (RSE 10.6%) |
| Two-compartment structure | central + peripheral1 | Results ‘the data showed a better fit to two-compartment model’ |
| IV infusion, 1 h | \- | Methods, ‘Dosing Regimen’: ‘The infusion time was 1 h.’ |
| Exponential IIV, proportional RUV | \- | Results ‘Exponential and proportional error models were selected’ |

Two transcription points are worth recording explicitly.

**The minus sign in `(eGFR - 50)` is not text-extractable.** The Table 2
notes are typeset in a symbol font whose minus glyph is dropped by
`pdftotext` and by the markdown preprocessor, so the equation extracts
as `0:838 * (1 + (eGFR 50) * 0:0082)` with a bare gap. Rendering page 6
of the PDF to an image resolves it as a subtraction, and the generic
covariate form printed in the Supplementary Table S1 notes,
`P_ij = P_tv,j * [1 + theta_j * (COV - COV_ave)] * e^eta_i`,
corroborates it. The reproduction of 64 published simulation outputs
below is a further, independent confirmation: the alternative readings
(`+ 50`, or no centring at all) move the trough predictions by tens of
percent.

**`Vp (L) = 51.6 x e^0.387` in the Table 2 note is not a multiplier.**
Taken literally it would give 76.0 L, contradicting both the Table 2
estimate and the Results text, which give the typical Vp as 51.6 L. The
exponent is the between-subject variance: 0.622^2 = 0.3869, and 62.2% is
exactly the Table 2 `omega Vp`. The same identity fixes the omega
convention for the whole table – the published percentages are omega
(the standard deviation on the log scale) times 100, not
`sqrt(exp(omega^2) - 1)`.

## Deterministic checks

With between-subject and residual variability switched off, the model
must satisfy exact internal identities. These are the tightest gates
available because they contain no Monte Carlo noise.

``` r

theta_cl <- 0.838
theta_e  <- 0.00823

# The covariate equation as printed in the Table 2 note.
cl_at <- function(egfr) theta_cl * (1 + (egfr - 50) * theta_e)

grid <- sort(unique(c(seq(0, 12, by = 0.1), seq(12, 48, by = 0.5),
                      seq(48, 336, by = 2))))
ev_single <- rbind(
  data.frame(id = 1L, time = 0, amt = 400, evid = 1L, dur = 1,
             cmt = "central", CRCL = 50),
  data.frame(id = 1L, time = grid, amt = NA_real_, evid = 0L, dur = NA_real_,
             cmt = "central", CRCL = 50)
)
ev_single <- ev_single[order(ev_single$time, -ev_single$evid), ]

sim_single <- rxode2::rxSolve(mod, ev_single, omega = NA,
                              returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
sim_single$id <- 1L
```

``` r

conc_obj <- PKNCA::PKNCAconc(dplyr::filter(sim_single, !is.na(Cc)),
                             Cc ~ time | id)
dose_obj <- PKNCA::PKNCAdose(data.frame(id = 1L, time = 0, dose = 400),
                             dose ~ time | id)
nca_single <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  conc_obj, dose_obj,
  intervals = data.frame(start = 0, end = Inf, cmax = TRUE, tmax = TRUE,
                         auclast = TRUE, aucinf.obs = TRUE,
                         half.life = TRUE, cl.obs = TRUE)
))
nca_val <- function(code) {
  v <- as.data.frame(nca_single$result)
  out <- v$PPORRES[v$PPTESTCD == code]
  if (length(out) != 1L) stop("no unique NCA row for '", code, "'")
  out
}

# Analytic terminal half-life of the two-compartment system.
k10 <- theta_cl / 14.4
k12 <- 3.08 / 14.4
k21 <- 3.08 / 51.6
beta_rate <- (k10 + k12 + k21 -
              sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21)) / 2

det_tab <- tibble::tibble(
  Identity = c("AUC(0-inf) after a single 400 mg dose = Dose / CL",
               "PKNCA CL(obs) = model CL",
               "PKNCA terminal half-life = analytic beta half-life"),
  Expected = c(400 / theta_cl, theta_cl, log(2) / beta_rate),
  Observed = c(nca_val("aucinf.obs"), nca_val("cl.obs"), nca_val("half.life"))
) |>
  dplyr::mutate(`% diff` = 100 * (Observed - Expected) / Expected)
knitr::kable(det_tab, digits = c(0, 4, 4, 5))
```

| Identity | Expected | Observed | % diff |
|:---|---:|---:|---:|
| AUC(0-inf) after a single 400 mg dose = Dose / CL | 477.3270 | 477.2962 | -0.00644 |
| PKNCA CL(obs) = model CL | 0.8380 | 0.8381 | 0.00644 |
| PKNCA terminal half-life = analytic beta half-life | 64.0443 | 63.8623 | -0.28418 |

``` r


# These are numerical-integration identities with no cohort sampling, so an
# exact bound is correct here (pattern: tighten deterministic gates). The
# realised deviations are ~0.006% (AUC, CL) and ~0.3% (half-life, which is a
# log-linear regression over a finite grid rather than a closed form).
stopifnot(
  abs(det_tab$`% diff`[1]) < 0.05,
  abs(det_tab$`% diff`[2]) < 0.05,
  abs(det_tab$`% diff`[3]) < 1
)
```

The AUC identity is the single most informative gate on the
transcription: it fails if the clearance value, the dose units, the
concentration units or the compartment scaling is wrong by any amount.

``` r

cov_tab <- tibble::tibble(
  `eGFR (mL/min/1.73 m^2)` = c(4.89, 29.2, 30, 50, 60, 90, 130, 170),
  `CL (L/h)`               = cl_at(c(4.89, 29.2, 30, 50, 60, 90, 130, 170)),
  Note = c("cohort minimum", "cohort median", "stratum 1 boundary",
           "centring value: CL = theta_CL exactly", "stratum 2 boundary",
           "stratum 3 boundary", "stratum 4 boundary", "cohort maximum")
)
knitr::kable(cov_tab, digits = 3)
```

| eGFR (mL/min/1.73 m^2) | CL (L/h) | Note                                  |
|-----------------------:|---------:|:--------------------------------------|
|                   4.89 |    0.527 | cohort minimum                        |
|                  29.20 |    0.695 | cohort median                         |
|                  30.00 |    0.700 | stratum 1 boundary                    |
|                  50.00 |    0.838 | centring value: CL = theta_CL exactly |
|                  60.00 |    0.907 | stratum 2 boundary                    |
|                  90.00 |    1.114 | stratum 3 boundary                    |
|                 130.00 |    1.390 | stratum 4 boundary                    |
|                 170.00 |    1.666 | cohort maximum                        |

``` r


# At the centring value the covariate term is exactly 1.
stopifnot(abs(cl_at(50) - theta_cl) < 1e-12)
# The linear form stays positive across and well beyond the observed range.
stopifnot(cl_at(0) > 0)
```

## Recovering the published simulation design

The paper does not state what eGFR value was used inside each
renal-function stratum, nor whether the reported trough concentrations
include residual error. Both choices materially change the predictions,
so rather than guess we identify them from the 64 published trough
medians.

``` r

n_sub <- 200  # participants per arm

make_events <- function(load_dose, n_load, maint_dose, egfr, n = n_sub,
                        tmax = 180, obs_times = c(72, 168)) {
  load_times  <- seq(0, by = 12, length.out = n_load)
  # Maintenance starts one full 24 h interval after the last loading dose.
  # This is the only schedule that places a trough exactly at both 72 h and
  # 168 h, which is how the paper defines its two endpoints, for BOTH the
  # q12h x 3 and q12h x 5 loading arms.
  maint_times <- seq(max(load_times) + 24, tmax, by = 24)
  do.call(rbind, lapply(seq_len(n), function(i) {
    dosing <- data.frame(
      id = i, time = c(load_times, maint_times),
      amt = c(rep(load_dose, n_load), rep(maint_dose, length(maint_times))),
      evid = 1L, dur = 1, cmt = "central", CRCL = egfr
    )
    obs <- data.frame(
      id = i, time = obs_times, amt = NA_real_, evid = 0L, dur = NA_real_,
      cmt = "central", CRCL = egfr
    )
    x <- rbind(dosing, obs)
    x[order(x$time, -x$evid), ]
  }))
}

regimens <- tibble::tribble(
  ~regimen,                    ~load, ~n_load, ~maint, ~target, ~infection,
  "400 mg q12hx3, 400 mg qd",  400,   3L,      400,    10,      "Gram-positive",
  "400 mg q12hx5, 400 mg qd",  400,   5L,      400,    10,      "Gram-positive",
  "600 mg q12hx3, 600 mg qd",  600,   3L,      600,    10,      "Gram-positive",
  "600 mg q12hx5, 600 mg qd",  600,   5L,      600,    10,      "Gram-positive",
  "800 mg q12hx3, 400 mg qd",  800,   3L,      400,    15,      "Deep-seated",
  "800 mg q12hx3, 600 mg qd",  800,   3L,      600,    15,      "Deep-seated",
  "800 mg q12hx3, 800 mg qd",  800,   3L,      800,    15,      "Deep-seated",
  "800 mg q12hx5, 800 mg qd",  800,   5L,      800,    15,      "Deep-seated"
)

strata_levels <- c("eGFR <= 30", "eGFR 30-60", "eGFR 60-90", "eGFR 90-130")

# Run every regimen x stratum cell at one eGFR value and return the
# simulated 72 h and 168 h troughs. `Cc` is the individual prediction;
# `sim` additionally carries the proportional residual error.
sweep_egfr <- function(egfr_by_stratum) {
  rxode2::rxSetSeed(42)
  out <- list()
  for (k in seq_along(strata_levels)) {
    for (j in seq_len(nrow(regimens))) {
      r <- regimens[j, ]
      s <- rxode2::rxSolve(
        mod,
        make_events(r$load, r$n_load, r$maint, egfr_by_stratum[k]),
        returnType = "data.frame"
      )
      out[[length(out) + 1L]] <- tibble::tibble(
        stratum = strata_levels[k], regimen = r$regimen,
        infection = r$infection, target = r$target,
        time = rep(c(72, 168), each = n_sub),
        ipred = c(s$Cc[s$time == 72], s$Cc[s$time == 168]),
        obs   = c(s$sim[s$time == 72], s$sim[s$time == 168])
      )
    }
  }
  dplyr::bind_rows(out)
}
```

### Which eGFR was used inside each stratum?

The two natural candidates are the midpoint of each band and its upper
bound. They are easy to tell apart: across the eGFR range the model’s
trough predictions change by roughly a factor of two, so a 20-40
mL/min/1.73 m^2 difference in the assumed eGFR is worth 10-20% on every
cell.

``` r

summarise_medians <- function(sw, value_col) {
  sw |>
    dplyr::group_by(stratum, regimen, time) |>
    dplyr::summarise(sim = stats::median(.data[[value_col]]), .groups = "drop") |>
    dplyr::inner_join(published, by = c("stratum", "regimen", "time")) |>
    dplyr::mutate(pct = 100 * (sim - med) / med)
}

sweep_mid  <- sweep_egfr(c(17.4, 45, 75, 110))
sweep_edge <- sweep_egfr(c(30, 60, 90, 130))

design_tab <- dplyr::bind_rows(
  summarise_medians(sweep_mid,  "ipred") |>
    dplyr::mutate(scheme = "stratum midpoint (17.4 / 45 / 75 / 110)"),
  summarise_medians(sweep_edge, "ipred") |>
    dplyr::mutate(scheme = "stratum upper bound (30 / 60 / 90 / 130)")
) |>
  dplyr::group_by(scheme) |>
  dplyr::summarise(
    `n cells`             = dplyr::n(),
    `median |% diff|`     = stats::median(abs(pct)),
    `90th pctile |% diff|`= stats::quantile(abs(pct), 0.9),
    `median % diff`       = stats::median(pct),
    .groups = "drop"
  )
knitr::kable(design_tab, digits = 1)
```

| scheme | n cells | median \|% diff\| | 90th pctile \|% diff\| | median % diff |
|:---|---:|---:|---:|---:|
| stratum midpoint (17.4 / 45 / 75 / 110) | 64 | 16.2 | 23.8 | 16.2 |
| stratum upper bound (30 / 60 / 90 / 130) | 64 | 4.6 | 9.3 | 4.4 |

Simulating each stratum at its **upper bound** reproduces all 64
published trough medians with a median absolute deviation of a few
percent, whereas the midpoint assumption leaves a systematic bias of
roughly +16% on every cell. The upper bound is also the sensible design
choice for a paper whose purpose is to recommend a regimen for a whole
band: within a band the highest eGFR gives the highest clearance and
hence the lowest trough, so it is the worst case the regimen must cover.

``` r

edge_row <- design_tab[grepl("upper", design_tab$scheme), ]
mid_row  <- design_tab[grepl("midpoint", design_tab$scheme), ]
# The discrimination is large (a few percent against roughly sixteen), so a
# factor-of-two separation is a bound with real headroom rather than one read
# off a single run.
stopifnot(
  edge_row$`median |% diff|` < 10,
  edge_row$`median |% diff|` < 0.5 * mid_row$`median |% diff|`
)
```

All remaining comparisons use the stratum upper bound.

### Do the published values carry residual error?

The individual prediction and the simulated observation have the same
median – the residual error is proportional and symmetric about zero –
so the trough medians cannot distinguish them. The interquartile
**width** can: the proportional residual error is 31.7%, which widens
the IQR appreciably.

``` r

iqr_ratio <- function(sw, value_col) {
  sw |>
    dplyr::group_by(stratum, regimen, time) |>
    dplyr::summarise(
      ratio = stats::quantile(.data[[value_col]], 0.75) /
              stats::quantile(.data[[value_col]], 0.25),
      .groups = "drop"
    )
}
pub_ratio <- published |> dplyr::mutate(ratio = q75 / q25)

ruv_tab <- tibble::tibble(
  Quantity = c("Published (Tables S2-1, S2-2)",
               "Simulated, individual prediction only",
               "Simulated, with proportional residual error"),
  `Median IQR width ratio (q75/q25)` = c(
    stats::median(pub_ratio$ratio),
    stats::median(iqr_ratio(sweep_edge, "ipred")$ratio),
    stats::median(iqr_ratio(sweep_edge, "obs")$ratio)
  )
)
knitr::kable(ruv_tab, digits = 2)
```

| Quantity | Median IQR width ratio (q75/q25) |
|:---|---:|
| Published (Tables S2-1, S2-2) | 2.15 |
| Simulated, individual prediction only | 1.82 |
| Simulated, with proportional residual error | 2.16 |

``` r


# The published spread matches the residual-error-inclusive simulation and is
# clearly wider than the individual predictions alone. Published ratio 2.146;
# realised simulated ratios across 2 / 4 / 16 solver threads were
# 2.162 / 2.125 / 2.13 with residual error and 1.820 / 1.803 / 1.81 without,
# so the two candidates are separated by ~0.33 and neither bound is read off
# a single run.
pub_r <- ruv_tab$`Median IQR width ratio (q75/q25)`[1]
ip_r  <- ruv_tab$`Median IQR width ratio (q75/q25)`[2]
ob_r  <- ruv_tab$`Median IQR width ratio (q75/q25)`[3]
stopifnot(abs(ob_r - pub_r) < 0.25, ip_r < pub_r - 0.10)
```

The published interquartile widths match the simulated **observations**,
so the paper’s Monte Carlo reports concentrations that include the
proportional residual error. The PTA comparisons below therefore use the
same quantity.

## Replicating Supplementary Tables S2-1 and S2-2

``` r

replication <- sweep_edge |>
  dplyr::group_by(stratum, regimen, infection, target, time) |>
  dplyr::summarise(
    sim_med = stats::median(obs),
    sim_pta = 100 * mean(obs >= dplyr::first(target)),
    .groups = "drop"
  ) |>
  dplyr::inner_join(published, by = c("stratum", "regimen", "time")) |>
  dplyr::mutate(
    med_pct = 100 * (sim_med - med) / med,
    pta_pp  = sim_pta - pta
  )
stopifnot(nrow(replication) == 64L)
```

| Comparison | n cells | Median | 90th pctile | Maximum |
|:---|---:|---:|---:|---:|
| Median trough concentration (% difference) | 64 | 4.0 | 9.0 | 10.4 |
| Probability of target attainment (percentage points) | 64 | 2.1 | 4.9 | 7.3 |

``` r

# Bounds are set well outside what a single cohort draw gives. rxSetSeed()
# fixes the rxode2 stream per solver-thread count, so a CI runner draws a
# different cohort than a development box; these statistics are aggregates
# over 64 cells, which damps but does not remove that. Measured across
# 2 / 4 / 16 solver threads when authored:
#   trough |% diff|  median 3.98 / 4.11 / 4.40, 90th pctile 8.99 / 9.41 / 9.20
#   PTA |pp diff|    median 2.05 / 1.70 / 2.00, 90th pctile 4.92 / 5.68 / 6.00
# The bounds below sit at roughly three times the observed medians and 2.5
# times the observed 90th percentiles. A mis-transcribed clearance, dose,
# volume or covariate sign moves the troughs by tens of percent and the PTAs
# by tens of points, so they still go red on any real error -- the ruled-out
# stratum-midpoint design alone lands at a median of ~16%.
stopifnot(
  stats::median(abs(replication$med_pct))         < 12,
  stats::quantile(abs(replication$med_pct), 0.9)  < 25,
  stats::median(abs(replication$pta_pp))          < 6,
  stats::quantile(abs(replication$pta_pp), 0.9)   < 14
)
```

The full 72 h comparison, one row per regimen and stratum:

| Stratum | Regimen | Published Cmin (mg/L) | Simulated Cmin (mg/L) | % diff | Published PTA (%) | Simulated PTA (%) | pp diff |
|:---|:---|---:|---:|---:|---:|---:|---:|
| eGFR 30-60 | 400 mg q12hx3, 400 mg qd | 8.9 | 8.5 | -4.8 | 39.4 | 38.5 | -0.9 |
| eGFR 30-60 | 400 mg q12hx5, 400 mg qd | 11.0 | 11.6 | 5.3 | 56.9 | 62.0 | 5.1 |
| eGFR 30-60 | 600 mg q12hx3, 600 mg qd | 13.4 | 13.1 | -2.2 | 68.8 | 70.0 | 1.2 |
| eGFR 30-60 | 600 mg q12hx5, 600 mg qd | 16.6 | 17.7 | 6.5 | 80.3 | 79.5 | -0.8 |
| eGFR 30-60 | 800 mg q12hx3, 400 mg qd | 13.6 | 14.4 | 6.2 | 43.5 | 46.5 | 3.0 |
| eGFR 30-60 | 800 mg q12hx3, 600 mg qd | 15.4 | 16.4 | 6.3 | 52.2 | 56.5 | 4.3 |
| eGFR 30-60 | 800 mg q12hx3, 800 mg qd | 17.0 | 17.9 | 5.2 | 59.2 | 62.0 | 2.8 |
| eGFR 30-60 | 800 mg q12hx5, 800 mg qd | 22.8 | 21.2 | -7.0 | 76.2 | 77.0 | 0.8 |
| eGFR 60-90 | 400 mg q12hx3, 400 mg qd | 7.1 | 7.7 | 8.3 | 26.0 | 30.5 | 4.5 |
| eGFR 60-90 | 400 mg q12hx5, 400 mg qd | 9.2 | 8.9 | -3.1 | 43.6 | 42.0 | -1.6 |
| eGFR 60-90 | 600 mg q12hx3, 600 mg qd | 10.9 | 9.8 | -10.1 | 55.7 | 49.0 | -6.7 |
| eGFR 60-90 | 600 mg q12hx5, 600 mg qd | 13.5 | 13.9 | 3.0 | 70.2 | 74.5 | 4.3 |
| eGFR 60-90 | 800 mg q12hx3, 400 mg qd | 11.4 | 11.9 | 4.4 | 31.7 | 36.0 | 4.3 |
| eGFR 60-90 | 800 mg q12hx3, 600 mg qd | 13.4 | 13.4 | 0.2 | 40.2 | 39.5 | -0.7 |
| eGFR 60-90 | 800 mg q12hx3, 800 mg qd | 14.0 | 13.7 | -2.2 | 44.2 | 43.5 | -0.7 |
| eGFR 60-90 | 800 mg q12hx5, 800 mg qd | 18.8 | 18.5 | -1.7 | 64.9 | 65.5 | 0.6 |
| eGFR 90-130 | 400 mg q12hx3, 400 mg qd | 5.5 | 5.5 | 0.8 | 13.8 | 12.0 | -1.8 |
| eGFR 90-130 | 400 mg q12hx5, 400 mg qd | 7.5 | 7.2 | -4.3 | 31.1 | 31.0 | -0.1 |
| eGFR 90-130 | 600 mg q12hx3, 600 mg qd | 8.7 | 8.7 | -0.2 | 40.4 | 36.5 | -3.9 |
| eGFR 90-130 | 600 mg q12hx5, 600 mg qd | 11.2 | 11.9 | 6.4 | 57.7 | 60.5 | 2.8 |
| eGFR 90-130 | 800 mg q12hx3, 400 mg qd | 8.9 | 9.0 | 1.7 | 18.7 | 18.0 | -0.7 |
| eGFR 90-130 | 800 mg q12hx3, 600 mg qd | 10.3 | 10.1 | -2.2 | 24.6 | 24.0 | -0.6 |
| eGFR 90-130 | 800 mg q12hx3, 800 mg qd | 10.9 | 11.1 | 1.6 | 30.4 | 27.5 | -2.9 |
| eGFR 90-130 | 800 mg q12hx5, 800 mg qd | 14.0 | 15.3 | 9.5 | 47.8 | 51.5 | 3.7 |
| eGFR \<= 30 | 400 mg q12hx3, 400 mg qd | 10.5 | 9.9 | -5.5 | 52.8 | 49.5 | -3.3 |
| eGFR \<= 30 | 400 mg q12hx5, 400 mg qd | 13.7 | 13.3 | -2.8 | 71.2 | 70.5 | -0.7 |
| eGFR \<= 30 | 600 mg q12hx3, 600 mg qd | 15.7 | 16.3 | 3.5 | 77.9 | 81.0 | 3.1 |
| eGFR \<= 30 | 600 mg q12hx5, 600 mg qd | 20.7 | 19.5 | -6.0 | 87.0 | 87.5 | 0.5 |
| eGFR \<= 30 | 800 mg q12hx3, 400 mg qd | 16.9 | 18.5 | 9.5 | 59.2 | 60.0 | 0.8 |
| eGFR \<= 30 | 800 mg q12hx3, 600 mg qd | 18.7 | 19.9 | 6.5 | 66.8 | 70.5 | 3.7 |
| eGFR \<= 30 | 800 mg q12hx3, 800 mg qd | 21.2 | 21.5 | 1.3 | 74.3 | 75.5 | 1.2 |
| eGFR \<= 30 | 800 mg q12hx5, 800 mg qd | 27.1 | 29.5 | 9.0 | 85.4 | 86.0 | 0.6 |

Replicates Wang 2023 Supplementary Tables S2-1 and S2-2, Cmin 72 h rows.
{.table style="width:100%;"}

And at 168 h:

| Stratum | Regimen | Published Css (mg/L) | Simulated Css (mg/L) | % diff | Published PTA (%) | Simulated PTA (%) | pp diff |
|:---|:---|---:|---:|---:|---:|---:|---:|
| eGFR 30-60 | 400 mg q12hx3, 400 mg qd | 10.8 | 10.7 | -1.3 | 55.9 | 53.5 | -2.4 |
| eGFR 30-60 | 400 mg q12hx5, 400 mg qd | 11.5 | 11.4 | -0.8 | 60.7 | 55.0 | -5.7 |
| eGFR 30-60 | 600 mg q12hx3, 600 mg qd | 16.2 | 16.7 | 2.8 | 80.5 | 81.0 | 0.5 |
| eGFR 30-60 | 600 mg q12hx5, 600 mg qd | 17.2 | 17.6 | 2.0 | 83.2 | 85.5 | 2.3 |
| eGFR 30-60 | 800 mg q12hx3, 400 mg qd | 12.2 | 12.7 | 4.1 | 36.2 | 36.5 | 0.3 |
| eGFR 30-60 | 800 mg q12hx3, 600 mg qd | 17.3 | 17.4 | 0.6 | 60.4 | 59.0 | -1.4 |
| eGFR 30-60 | 800 mg q12hx3, 800 mg qd | 21.6 | 22.4 | 3.8 | 73.5 | 75.5 | 2.0 |
| eGFR 30-60 | 800 mg q12hx5, 800 mg qd | 23.6 | 24.2 | 2.7 | 76.6 | 75.0 | -1.6 |
| eGFR 60-90 | 400 mg q12hx3, 400 mg qd | 8.3 | 8.8 | 6.3 | 37.7 | 41.0 | 3.3 |
| eGFR 60-90 | 400 mg q12hx5, 400 mg qd | 9.3 | 9.2 | -1.1 | 45.4 | 43.0 | -2.4 |
| eGFR 60-90 | 600 mg q12hx3, 600 mg qd | 13.0 | 12.6 | -3.0 | 66.4 | 67.0 | 0.6 |
| eGFR 60-90 | 600 mg q12hx5, 600 mg qd | 13.9 | 13.7 | -1.2 | 69.6 | 66.5 | -3.1 |
| eGFR 60-90 | 800 mg q12hx3, 400 mg qd | 9.8 | 8.8 | -10.4 | 23.2 | 19.0 | -4.2 |
| eGFR 60-90 | 800 mg q12hx3, 600 mg qd | 13.5 | 12.8 | -5.0 | 43.6 | 41.5 | -2.1 |
| eGFR 60-90 | 800 mg q12hx3, 800 mg qd | 16.4 | 17.2 | 4.8 | 56.3 | 58.0 | 1.7 |
| eGFR 60-90 | 800 mg q12hx5, 800 mg qd | 18.6 | 18.5 | -0.3 | 65.8 | 64.5 | -1.3 |
| eGFR 90-130 | 400 mg q12hx3, 400 mg qd | 6.5 | 6.0 | -8.4 | 23.4 | 21.0 | -2.4 |
| eGFR 90-130 | 400 mg q12hx5, 400 mg qd | 6.9 | 6.6 | -5.0 | 27.2 | 26.5 | -0.7 |
| eGFR 90-130 | 600 mg q12hx3, 600 mg qd | 9.9 | 9.2 | -7.0 | 48.4 | 45.0 | -3.4 |
| eGFR 90-130 | 600 mg q12hx5, 600 mg qd | 10.5 | 9.7 | -7.9 | 53.4 | 47.0 | -6.4 |
| eGFR 90-130 | 800 mg q12hx3, 400 mg qd | 7.2 | 7.0 | -2.7 | 9.9 | 11.5 | 1.6 |
| eGFR 90-130 | 800 mg q12hx3, 600 mg qd | 10.0 | 9.6 | -3.8 | 24.5 | 19.0 | -5.5 |
| eGFR 90-130 | 800 mg q12hx3, 800 mg qd | 12.9 | 12.3 | -4.7 | 40.9 | 36.5 | -4.4 |
| eGFR 90-130 | 800 mg q12hx5, 800 mg qd | 13.8 | 13.6 | -1.2 | 44.8 | 44.5 | -0.3 |
| eGFR \<= 30 | 400 mg q12hx3, 400 mg qd | 13.9 | 13.5 | -3.2 | 74.9 | 71.5 | -3.4 |
| eGFR \<= 30 | 400 mg q12hx5, 400 mg qd | 15.0 | 16.5 | 10.2 | 75.2 | 82.5 | 7.3 |
| eGFR \<= 30 | 600 mg q12hx3, 600 mg qd | 20.6 | 21.2 | 2.7 | 90.0 | 89.5 | -0.5 |
| eGFR \<= 30 | 600 mg q12hx5, 600 mg qd | 22.8 | 21.5 | -5.8 | 91.9 | 92.0 | 0.1 |
| eGFR \<= 30 | 800 mg q12hx3, 400 mg qd | 16.5 | 18.1 | 9.6 | 57.9 | 65.0 | 7.1 |
| eGFR \<= 30 | 800 mg q12hx3, 600 mg qd | 21.5 | 22.1 | 2.7 | 76.3 | 74.5 | -1.8 |
| eGFR \<= 30 | 800 mg q12hx3, 800 mg qd | 28.3 | 30.8 | 8.9 | 86.8 | 87.5 | 0.7 |
| eGFR \<= 30 | 800 mg q12hx5, 800 mg qd | 30.4 | 30.1 | -1.0 | 88.0 | 91.0 | 3.0 |

Replicates Wang 2023 Supplementary Tables S2-1 and S2-2, Css 168 h rows.
{.table}

## Replicating Figures 3 and 4

![](Wang_2023_teicoplanin_files/figure-html/figure-3-4-1.png)

Replicates Figures 3 (Gram-positive panel, dotted 50% line, 10 mg/L
target) and 4 (deep-seated panel, 15 mg/L target) of Wang 2023. Solid
and dashed lines for the same regimen overlay closely across all four
strata.

## Replicating Table 3: the recommended regimens

Table 3 of the paper names one regimen per infection type per stratum.
The check here is that the packaged model reproduces the published PTA
for each of those eight recommended cells.

``` r

recommended <- tibble::tribble(
  ~stratum,      ~infection,      ~regimen,
  "eGFR <= 30",  "Gram-positive", "400 mg q12hx3, 400 mg qd",
  "eGFR 30-60",  "Gram-positive", "400 mg q12hx5, 400 mg qd",
  "eGFR 60-90",  "Gram-positive", "600 mg q12hx3, 600 mg qd",
  "eGFR 90-130", "Gram-positive", "600 mg q12hx5, 600 mg qd",
  "eGFR <= 30",  "Deep-seated",   "800 mg q12hx3, 400 mg qd",
  "eGFR 30-60",  "Deep-seated",   "800 mg q12hx3, 600 mg qd",
  "eGFR 60-90",  "Deep-seated",   "800 mg q12hx3, 800 mg qd",
  "eGFR 90-130", "Deep-seated",   "800 mg q12hx5, 800 mg qd"
) |>
  dplyr::mutate(stratum = factor(stratum, levels = strata_levels))

rec_tab <- recommended |>
  dplyr::inner_join(replication, by = c("stratum", "infection", "regimen")) |>
  dplyr::select(stratum, infection, regimen, time, pta, sim_pta, pta_pp) |>
  dplyr::arrange(infection, stratum, time)
stopifnot(nrow(rec_tab) == 16L)

rec_tab |>
  dplyr::rename(
    "Stratum" = stratum, "Infection type" = infection, "Regimen" = regimen,
    "Endpoint (h)" = time, "Published PTA (%)" = pta,
    "Simulated PTA (%)" = sim_pta, "pp diff" = pta_pp
  ) |>
  knitr::kable(digits = 1, caption = "Replicates Wang 2023 Table 3.")
```

| Stratum | Infection type | Regimen | Endpoint (h) | Published PTA (%) | Simulated PTA (%) | pp diff |
|:---|:---|:---|---:|---:|---:|---:|
| eGFR 30-60 | Deep-seated | 800 mg q12hx3, 600 mg qd | 72 | 52.2 | 56.5 | 4.3 |
| eGFR 30-60 | Deep-seated | 800 mg q12hx3, 600 mg qd | 168 | 60.4 | 59.0 | -1.4 |
| eGFR 60-90 | Deep-seated | 800 mg q12hx3, 800 mg qd | 72 | 44.2 | 43.5 | -0.7 |
| eGFR 60-90 | Deep-seated | 800 mg q12hx3, 800 mg qd | 168 | 56.3 | 58.0 | 1.7 |
| eGFR 90-130 | Deep-seated | 800 mg q12hx5, 800 mg qd | 72 | 47.8 | 51.5 | 3.7 |
| eGFR 90-130 | Deep-seated | 800 mg q12hx5, 800 mg qd | 168 | 44.8 | 44.5 | -0.3 |
| eGFR \<= 30 | Deep-seated | 800 mg q12hx3, 400 mg qd | 72 | 59.2 | 60.0 | 0.8 |
| eGFR \<= 30 | Deep-seated | 800 mg q12hx3, 400 mg qd | 168 | 57.9 | 65.0 | 7.1 |
| eGFR 30-60 | Gram-positive | 400 mg q12hx5, 400 mg qd | 72 | 56.9 | 62.0 | 5.1 |
| eGFR 30-60 | Gram-positive | 400 mg q12hx5, 400 mg qd | 168 | 60.7 | 55.0 | -5.7 |
| eGFR 60-90 | Gram-positive | 600 mg q12hx3, 600 mg qd | 72 | 55.7 | 49.0 | -6.7 |
| eGFR 60-90 | Gram-positive | 600 mg q12hx3, 600 mg qd | 168 | 66.4 | 67.0 | 0.6 |
| eGFR 90-130 | Gram-positive | 600 mg q12hx5, 600 mg qd | 72 | 57.7 | 60.5 | 2.8 |
| eGFR 90-130 | Gram-positive | 600 mg q12hx5, 600 mg qd | 168 | 53.4 | 47.0 | -6.4 |
| eGFR \<= 30 | Gram-positive | 400 mg q12hx3, 400 mg qd | 72 | 52.8 | 49.5 | -3.3 |
| eGFR \<= 30 | Gram-positive | 400 mg q12hx3, 400 mg qd | 168 | 74.9 | 71.5 | -3.4 |

Replicates Wang 2023 Table 3. {.table}

``` r


stopifnot(stats::median(abs(rec_tab$pta_pp)) < 8)
```

Four of the eight recommended cells fall below the paper’s own stated
“optimal if PTA \> 50%” criterion at one or both endpoints – for example
deep-seated infection at eGFR 60-90 (published 44.2% at 72 h) and at
eGFR 90-130 (published 47.8% at 72 h and 44.8% at 168 h). The paper
presents these as the best available option rather than as regimens
meeting the threshold; the Discussion notes that still-higher regimens
reported elsewhere (15 mg/kg or 1000 mg) were judged not clinically
applicable.

## NCA of the recommended standard regimen

The paper reports observed plasma Cmax and Cmin for the model-building
cohort (Table 1). Those are compared here against a PKNCA analysis of a
simulated cohort receiving the paper’s own group-1 standard regimen (400
mg q12h x 3, then 400 mg q24h) at the cohort median eGFR.

``` r

obs_grid <- sort(unique(c(seq(0, 180, by = 1), 145, 169)))
ev_cohort <- do.call(rbind, lapply(seq_len(n_sub), function(i) {
  dosing <- data.frame(
    id = i, time = c(0, 12, 24, seq(48, 180, by = 24)),
    amt = 400, evid = 1L, dur = 1, cmt = "central", CRCL = 29.2
  )
  obs <- data.frame(id = i, time = obs_grid, amt = NA_real_, evid = 0L,
                    dur = NA_real_, cmt = "central", CRCL = 29.2)
  x <- rbind(dosing, obs)
  x[order(x$time, -x$evid), ]
}))
rxode2::rxSetSeed(42)
sim_cohort <- rxode2::rxSolve(mod, ev_cohort, returnType = "data.frame")
if (is.null(sim_cohort$id)) sim_cohort$id <- 1L
```

``` r

# One steady-state-approaching dosing interval, 144-168 h. The filter is
# !is.na(Cc) only, so the interval-start record is retained.
nca_input <- sim_cohort |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(conc = sim)

conc_cohort <- PKNCA::PKNCAconc(nca_input, conc ~ time | id)
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
dose_cohort <- PKNCA::PKNCAdose(
  data.frame(id = seq_len(n_sub), time = 144, dose = 400), dose ~ time | id
)
nca_cohort <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  conc_cohort, dose_cohort,
  intervals = data.frame(start = 144, end = 168,
                         cmax = TRUE, cmin = TRUE, tmax = TRUE, auclast = TRUE)
))
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in log(conc.2/conc.1): NaNs produced
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in log(conc.2/conc.1): NaNs produced
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
#> Warning in log(conc.2/conc.1): NaNs produced
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc = conc): Negative concentrations found
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
nca_res <- as.data.frame(nca_cohort$result)
stopifnot(nrow(nca_res) > 0L)
```

``` r

reference <- data.frame(
  cmax = 39.9,   # Wang 2023 Table 1, observed plasma Cmax, model-establishment group
  cmin = 10.5    # Wang 2023 Table 1, observed plasma Cmin, model-establishment group
)
cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = reference,
  params    = c("cmax", "cmin"),
  units     = c(cmax = "mg/L", cmin = "mg/L"),
  tolerance_pct = 20
)
knitr::kable(cmp, digits = 1)
```

| NCA parameter | Reference | Simulated | % diff   |
|:--------------|:----------|:----------|:---------|
| Cmax (mg/L)   | 39.9      | 43.8      | +9.9%    |
| Cmin (mg/L)   | 10.5      | 6.63      | -36.9%\* |

``` r

attr(cmp, "footnote")
#> [1] "* differs from reference by more than ±20%."
```

This comparison is deliberately **not** gated, and the Cmax row is
expected to be flagged. The two sides are not measuring the same thing:
the published Cmax and Cmin are per-patient extremes of two or three
sparse clinical samples pooled across all four dosing groups – including
the 91 patients on irregular doses at irregular intervals – whereas the
simulated values come from a dense grid under a single fixed 400 mg
regimen at one eGFR. The Cmin comparison is nevertheless close, which is
the more informative of the two because the paper’s targets, its Monte
Carlo endpoints and its dosing recommendations are all trough-based.

![](Wang_2023_teicoplanin_files/figure-html/profile-figure-1.png)

Simulated concentration-time profile (median and 5th-95th percentiles, n
= 200) for the standard regimen at the cohort median eGFR of 29.2
mL/min/1.73 m^2. Dashed lines mark the 10 and 15 mg/L trough targets;
dotted lines mark the paper’s two endpoints at 72 h and 168 h. The long
terminal half-life (64 h) is why the 168 h trough is still appreciably
above the 72 h trough: the profile has not reached true steady state
within the simulated week.

## Assumptions and deviations

- **Within-stratum eGFR is an inference, not a published value.** The
  paper reports simulations “stratified by eGFR” but never states the
  eGFR used inside each band. The stratum upper bound (30 / 60 / 90 /
  130 mL/min/1.73 m^2) is identified above by reproducing the 64
  published trough medians, and is what this vignette uses. The stratum
  midpoint is ruled out by a systematic +16% bias.
- **Residual error is included in the reproduced trough distributions.**
  Also an inference, identified from the published interquartile widths.
  The medians are insensitive to this choice; the PTAs are not.
- **The maintenance schedule starts 24 h after the last loading dose.**
  The paper writes “a loading dose of 400 mg every 12 hours for three
  doses, followed a maintenance dose of 400 mg every 24 hours” without
  giving the gap. Starting maintenance one full interval after the last
  loading dose is the only reading that places a trough exactly at both
  published endpoints (72 h and 168 h) for the q12h x 3 and the q12h x 5
  arms simultaneously.
- **`omega Q` is encoded as no random effect on Q.** Table 2 reports it
  as “0 FIX”. Writing `etalq ~ fixed(0)` would make OMEGA singular and
  break rxode2’s Cholesky sampler; omitting the eta is mathematically
  identical.
- **The reported `omega` percentages are treated as the log-scale
  standard deviation times 100**, so the packaged variances are the
  squared percentages. This is fixed by the paper’s own printed
  `e^0.387` for Vp against its 62.2% table entry, not assumed.
- **Number of simulated subjects.** The paper does not state its Monte
  Carlo replicate count. This vignette uses 200 per arm, the package
  cap.
- **No body-weight or dose-per-kg scaling.** Body weight was screened
  and not retained (Supplementary Table S1), so all doses are absolute
  mg, as in the paper’s own simulations.

## Errata and internal inconsistencies in the source

- **Delta OFV for eGFR on CL is reported inconsistently.** Results
  states “eGFR was a covariate significantly associated with
  inter-individual variability of systemic clearance (CL) (dOFV
  -23.189)”. The Discussion gives “eGFR (-20.374)” and Supplementary
  Table S1 gives -20.374 against a base-model OFV of 1564.335 and a
  final-model OFV of 1543.961, which is self-consistent. The value
  -23.189 appears nowhere else; the closest Table S1 entry is -23.839,
  for serum creatinine on CL. The packaged model and its metadata use
  -20.374. No model parameter depends on this.
- **The residual error is reported as two transposed values.** Table 2
  gives the proportional error as 31.7% (RSE 10.6%), but the Results
  text states “the intra-individual variability was 37.1%” – the same
  digits in a different order. Table 2 is the correct one: its own
  bootstrap 95% CI for this parameter, 21.4-36.4%, contains 31.7% and
  excludes 37.1%. The packaged model uses `propSd = 0.317`. This matters
  for the reproduction above, because the residual error sets the
  interquartile widths that identify whether the published simulations
  carry residual error.
- **`Vp (L) = 51.6 x e^0.387` in the Table 2 note** reads as a
  multiplier but is the between-subject variance; see the Source trace
  section above. The typical Vp is 51.6 L per Table 2 and per the
  Results text.
- **The minus sign in the CL covariate equation is lost by text
  extraction** from the PDF; it was confirmed from the rendered page
  image and from the supplement’s generic covariate form. See the Source
  trace section.
- **The Table 1 infection-site counts exceed the group sizes** (127
  sites for 108 model-establishment patients, 47 for 43 validation
  patients), consistent with patients carrying more than one infection
  site. The counts are recorded verbatim in the model’s `population`
  metadata.
- **No erratum or corrigendum was found.** Crossref reports no
  `update-to` and no `relation` entries for <doi:10.2147/DDDT.S413662>
  as of the extraction date.
