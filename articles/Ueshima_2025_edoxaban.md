# Edoxaban (Ueshima 2025)

## Model and source

``` r

mod <- readModelDb("Ueshima_2025_edoxaban")
ui  <- rxode2::rxode(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Ueshima S, Hira D, Matsuda S, Michihata R, Tabuchi Y, Ozawa
  T, Itoh H, Iguchi M, Akao M, Aizawa T, Kashiwa A, Shizuta S, Makiyama
  T, Nakagawa Y, Horie M, Terada T, Katsura T. Population
  pharmacokinetics and pharmacogenomics of edoxaban in Japanese adults
  with atrial fibrillation. J Pharm Health Care Sci. 2025;11:46.
  <doi:10.1186/s40780-025-00453-2>
- Description: One-compartment population PK model for oral edoxaban in
  Japanese adults with non-valvular atrial fibrillation, fitted to
  sparse real-world therapeutic-drug-monitoring data (one steady-state
  sample per patient, 9-47 h after the last dose). Absorption is not
  modelled: with no data on the absorption phase the authors used NONMEM
  ADVAN1 TRANS2, so a dose enters the central compartment directly.
  Apparent oral clearance carries a power effect of raw Cockcroft-Gault
  creatinine clearance; the apparent volume of distribution is fixed to
  the ENGAGE AF-TIMI 48 Asian value with a body-weight exponent fixed
  at 1. CYP3A5\*3 and the ABCB1 1236C\>T / 2677G\>T,A / 3435C\>T
  polymorphisms (and the 2677/3435 haplotype) were screened and none
  affected edoxaban pharmacokinetics (Ueshima 2025).
- Article: <https://doi.org/10.1186/s40780-025-00453-2> (open access,
  PMC12131839)

## Population

Ueshima 2025 is a retrospective real-world therapeutic-drug-monitoring
study of 131 Japanese adults with atrial fibrillation treated with
once-daily oral edoxaban (Lixiana) at three Japanese centres between
January 2017 and August 2019. Each patient contributed exactly **one**
plasma sample, drawn at steady state 9 to 47 h after the last dose, so
the analysis dataset holds 131 observations from 131 patients and
contains no absorption-phase information at all. That single fact drives
every structural choice in the model: NONMEM ADVAN1 TRANS2 (one
compartment, no absorption, dose straight into the central compartment),
an apparent volume of distribution fixed to an external value, and an
apparent clearance as the only estimated structural parameter.

Baseline characteristics (Ueshima 2025 Table 1, median with range): age
72.2 years (35.2-92.6), body weight 61.4 kg (37.3-97.3), sex 83 male /
48 female, serum creatinine 0.89 mg/dL (0.56-2.00), Cockcroft-Gault
creatinine clearance 61.8 mL/min (20.0-135.5), AST 25 IU/L (11-92), ALT
19 IU/L (5-77). Daily doses were 15 mg (n = 3), 30 mg (n = 70) and 60 mg
(n = 58). Sixteen patients took a CYP3A4 and/or P-glycoprotein inhibitor
(amiodarone 5, diltiazem 6, verapamil 5). Observed edoxaban
concentrations were 19.5 ng/mL (1.7-152.0), with a lower limit of
quantification of 1 ng/mL.

The paper’s headline results are **negative** pharmacogenomic findings:
`CYP3A5*3` (rs776746) and the `ABCB1` 1236C\>T / 2677G\>T,A / 3435C\>T
polymorphisms, and the 2677/3435 haplotype, were all screened and none
affected edoxaban pharmacokinetics. Those screened-but-not-retained
covariates are recorded in the model file’s `covariatesDataExcluded`
list so the provenance of the screen survives with the model.

``` r

str(ui$population)
#> List of 17
#>  $ species         : chr "human"
#>  $ n_subjects      : int 131
#>  $ n_studies       : int 1
#>  $ n_observations  : int 131
#>  $ age_range       : chr "35.2-92.6 years"
#>  $ age_median      : chr "72.2 years"
#>  $ weight_range    : chr "37.3-97.3 kg"
#>  $ weight_median   : chr "61.4 kg"
#>  $ sex_female_pct  : num 36.6
#>  $ race_ethnicity  : chr "Japanese"
#>  $ disease_state   : chr "Adults with atrial fibrillation on chronic once-daily oral edoxaban (Lixiana) for prevention of cardioembolic stroke"
#>  $ dose_range      : chr "15-60 mg once daily (15 mg n = 3, 30 mg n = 70, 60 mg n = 58)"
#>  $ renal_function  : chr "Cockcroft-Gault creatinine clearance median 61.8 mL/min, range 20.0-135.5; serum creatinine median 0.89 mg/dL, range 0.56-2.00"
#>  $ hepatic_function: chr "AST median 25 IU/L (range 11-92); ALT median 19 IU/L (range 5-77)"
#>  $ co_medication   : chr "CYP3A4 and/or P-glycoprotein inhibitors in 16 of 131 patients: amiodarone 5, diltiazem 6, verapamil 5"
#>  $ regions         : chr "Three centres in Japan: Shiga University of Medical Science Hospital, National Hospital Organization Kyoto Medi"| __truncated__
#>  $ notes           : chr "Retrospective real-world therapeutic-drug-monitoring cohort. ONE blood sample per patient, drawn at steady stat"| __truncated__
```

## Source trace

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL/F at CRCL 61.8 mL/min) | 28.2 L/h | Table 2, `theta1` (95% CI 26.9-29.5; bootstrap median 28.1); Eq. 10 |
| `e_crcl_cl` (power exponent on CRCL/61.8) | 0.692 | Table 2, `theta3` (95% CI 0.582-0.801; bootstrap median 0.693); Eq. 10 |
| `lvc` (Vd/F at WT 70 kg) | 336.4 L, fixed | Table 2, `theta2` (“336.4 fixed”); Methods Eq. 3 = 193 L central volume x 1.743 Asian fold-change from reference \[22\] (Krekels 2016, ENGAGE AF-TIMI 48) |
| `e_wt_vc` (exponent on WT/70) | 1, fixed | Methods, PPK modeling: “with the exponent fixed to 1 due to a lack of individual data”; Eq. 3 |
| `etalcl` (IIV on CL/F) | 26.4 CV%, variance 0.0697 | Table 2, `omega1` (95% CI 18.8-34.0; bootstrap median 25.9). Footnote c: the parenthesised 27.4 is the eta SHRINKAGE, not an RSE |
| `expSd` (residual) | 58.7 CV%, fixed | Table 2, `sigma` (“58.7 fixed”); Methods Eq. 1 and text: 14.6% (healthy subjects) + 44.1% (AF-patient increment), both from reference \[22\] |
| `d/dt(central)` (1 cmt, no absorption) | n/a | Methods, PPK modeling: “a 1-compartment model without first-order absorption was employed … (ADVAN1 TRANS2)” |
| `Cc = central / vc * 1000` | n/a | Unit conversion only: mg / L = ug/mL; the paper reports ng/mL throughout (Table 1, Figs. 2-3) |
| Reference CRCL 61.8 mL/min | n/a | Table 1 cohort median, and printed inside Eq. 10 and the Table 2 parameter header |
| Reference WT 70 kg | n/a | Eq. 3 (inherited from the ENGAGE AF-TIMI 48 parameterisation, not a cohort median) |

Note on the variability scale: Table 2 labels the two rows
`omega1 (CV%)` and `sigma (CV%)`, and footnote b states they are “the
coefficient of variations of the inter-individual variability for CL/F
and residual variability”. The reported number is therefore the
standard-deviation parameter itself expressed as a percentage – the
usual NONMEM `%CV = 100 * sqrt(omega^2)` convention – so
`omega1 = 0.264` (variance 0.0697) and `sigma = 0.587` on the log scale.
The “two variability terms” gate below tests that reading term by term.

## Structural check against the paper’s own printed clearances

The Discussion prints CL/F for three renal-function scenarios, computed
by the authors from Eq. 10. Those three numbers are an exact answer key
for the transcription of `lcl` and `e_crcl_cl`, and they are
deterministic – no simulation, no cohort – so the assertion can be
exact.

``` r

cl_paper <- function(crcl) 28.2 * (crcl / 61.8)^0.692   # Ueshima 2025 Eq. 10
v_paper  <- function(wt)   336.4 * (wt / 70)            # Ueshima 2025 Eq. 3

answer_key <- tibble::tibble(
  CRCL          = c(100, 60, 30),
  scenario      = c("Normal renal function", "Mild renal impairment",
                    "Severe renal impairment"),
  published_cl  = c(39.3, 27.6, 17.1),   # Ueshima 2025 Discussion, paragraph 3
  model_cl      = round(cl_paper(c(100, 60, 30)), 1)
)
knitr::kable(
  answer_key |>
    dplyr::rename(
      "CLcr (mL/min)"        = CRCL,
      "Scenario"             = scenario,
      "Published CL/F (L/h)" = published_cl,
      "Encoded CL/F (L/h)"   = model_cl
    ),
  caption = "Ueshima 2025 Discussion CL/F answer key vs. the encoded equation."
)
```

| CLcr (mL/min) | Scenario | Published CL/F (L/h) | Encoded CL/F (L/h) |
|---:|:---|---:|---:|
| 100 | Normal renal function | 39.3 | 39.3 |
| 60 | Mild renal impairment | 27.6 | 27.6 |
| 30 | Severe renal impairment | 17.1 | 17.1 |

Ueshima 2025 Discussion CL/F answer key vs. the encoded equation.
{.table}

``` r


stopifnot(identical(answer_key$model_cl, answer_key$published_cl))
```

The paper additionally compares those to the ENGAGE AF-TIMI 48 values
(34.6 / 24.8 / 18.5 L/h); we do not assert against those, since they
come from a different model.

## Virtual cohort

Original observed data are not publicly available. Two virtual cohorts
are used below.

The first replicates the paper’s own Monte Carlo design (Ueshima 2025,
Model-based simulation): patients at fixed creatinine clearances of 30,
60, 90 and 120 mL/min, dosed 30 mg once daily at CLcr 30 and 60 mg once
daily otherwise, per the Japanese package insert. The paper does not
state the body weight used for that simulation; 70 kg (the model’s own
reference weight, so Vd/F = 336.4 L exactly) reproduces the Figure 3
medians and is used here – see Assumptions and deviations.

Steady state is imposed directly with `ss = 1`, `ii = 24` rather than by
simulating a burn-in, so the profile is the exact steady-state solution.

``` r

# set.seed() seeds R's RNG, NOT rxode2's simulation RNG, and rxode2's streams
# are partitioned per solver thread -- so this cohort differs between a
# 16-thread workstation and a 2-core CI runner. Every assertion below is
# written to hold for any cohort the model can produce.
set.seed(20250908)

n_per_arm <- 200L
tau       <- 24

arms <- tibble::tibble(
  arm  = factor(
    c("CLcr 30, 30 mg", "CLcr 60, 60 mg", "CLcr 90, 60 mg", "CLcr 120, 60 mg"),
    levels = c("CLcr 30, 30 mg", "CLcr 60, 60 mg", "CLcr 90, 60 mg", "CLcr 120, 60 mg")
  ),
  CRCL = c(30, 60, 90, 120),
  dose = c(30, 60, 60, 60)
)

make_arm <- function(i, id_offset) {
  ids <- id_offset + seq_len(n_per_arm)
  covs <- list(WT = 70, CRCL = arms$CRCL[i], arm = arms$arm[i])
  dosing <- tibble::tibble(
    id = ids, time = 0, amt = arms$dose[i], evid = 1L, cmt = "central",
    ii = tau, ss = 1L
  )
  obs <- tidyr::crossing(id = ids, time = seq(0, tau, by = 0.25)) |>
    dplyr::mutate(amt = NA_real_, evid = 0L, cmt = "central", ii = 0, ss = 0L)
  dplyr::bind_rows(dosing, obs) |>
    dplyr::mutate(WT = covs$WT, CRCL = covs$CRCL, arm = covs$arm) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

events <- dplyr::bind_rows(
  lapply(seq_len(nrow(arms)), function(i) make_arm(i, (i - 1L) * n_per_arm))
)
stopifnot(
  !anyDuplicated(unique(events[, c("id", "time", "evid")])),
  length(unique(events$id)) == n_per_arm * nrow(arms)
)
```

## Simulation

``` r

sim <- rxode2::rxSolve(mod, events = as.data.frame(events),
                       keep = c("arm", "WT", "CRCL")) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'

# `Cc` is the individual prediction (structural model + eta, no residual
# error); `sim` additionally carries the log-normal residual and is therefore
# the column that corresponds to a simulated OBSERVED concentration, which is
# what Ueshima 2025 Figure 3 plots.
stopifnot(all(c("Cc", "sim", "cl", "vc") %in% names(sim)), all(sim$Cc > 0))
```

### Gate: the packaged model really is the paper’s ODE

The model writes an explicit `d/dt(central)` but names its parameters
`cl` and `vc`, which lets rxode2 recognise a one-compartment linear
system and solve it analytically. That substitution is only safe because
the ODE’s rate constant IS `cl / vc`; the two checks below confirm the
packaged object behaves as the paper’s equations say it should, rather
than silently solving something else.

``` r

# 1. Typical-value profile must equal the analytical steady-state solution of
#    Eqs. 3 and 10. Both sides use the same parameters, so the only difference
#    is numerical -- a tight bound is correct here.
mod0 <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
# One representative subject per arm (the first ID of each arm's block).
typical_ids <- (seq_len(nrow(arms)) - 1L) * n_per_arm + 1L
sim0 <- rxode2::rxSolve(
  mod0,
  events = as.data.frame(dplyr::filter(events, id %in% typical_ids)),
  keep = c("arm", "WT", "CRCL")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> Warning: multi-subject simulation without without 'omega'

closed_form <- function(t, dose, crcl, wt = 70) {
  k <- cl_paper(crcl) / v_paper(wt)
  (dose / v_paper(wt)) * 1000 * exp(-k * t) / (1 - exp(-k * tau))
}

sim0 <- sim0 |>
  dplyr::mutate(
    dose     = arms$dose[match(CRCL, arms$CRCL)],
    analytic = closed_form(time, dose, CRCL),
    pct_diff = 100 * (Cc - analytic) / analytic
  )
cat(sprintf("max |%% diff| typical-value solve vs analytical solution: %.5f%%\n",
            max(abs(sim0$pct_diff))))
#> max |% diff| typical-value solve vs analytical solution: 0.00010%
stopifnot(max(abs(sim0$pct_diff)) < 0.01)

# 2. Perturbing the estimated clearance must move the profile. A model whose
#    ODE had been discarded in favour of an unrelated solved form would not
#    respond correctly to this.
mod_slow <- rxode2::zeroRe(mod) |> rxode2::ini(lcl = log(28.2 / 2))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ change initial estimate of `lcl` to `2.64617479738412`
sim_slow <- rxode2::rxSolve(
  mod_slow,
  events = as.data.frame(dplyr::filter(events, id == n_per_arm + 1L))
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl'
trough_ref  <- sim0$Cc[sim0$CRCL == 60 & abs(sim0$time - tau) < 1e-8]
trough_slow <- sim_slow$Cc[abs(sim_slow$time - tau) < 1e-8]
# Halving CL/F raises the 24 h trough of a mono-exponential; the analytical
# ratio for this arm is deterministic, so assert it exactly.
ratio_expected <- closed_form(tau, 60, 60) / {
  k2 <- (cl_paper(60) / 2) / v_paper(70)
  (60 / v_paper(70)) * 1000 * exp(-k2 * tau) / (1 - exp(-k2 * tau))
}
cat(sprintf("trough ratio (published CL / halved CL): observed %.4f, analytic %.4f\n",
            trough_ref / trough_slow, ratio_expected))
#> trough ratio (published CL / halved CL): observed 0.2718, analytic 0.2718
stopifnot(abs(trough_ref / trough_slow / ratio_expected - 1) < 1e-3)
```

## Replicate Figure 3 – simulated steady-state trough concentrations

Ueshima 2025 Figure 3 plots box-and-whisker summaries of the predicted
trough concentration 24 h after the last dose at four creatinine
clearances, over the grey band and dashed line that mark the
interquartile range (13.9-47.4 ng/mL) and median (25.0 ng/mL) of
observed trough concentrations in Asians reported by Chao et
al. (reference \[10\]).

``` r

chao <- list(median = 25.0, q25 = 13.9, q75 = 47.4)   # Ueshima 2025 Results / Fig. 3 caption
lloq <- 1.0                                            # Ueshima 2025 Edoxaban assay

trough <- sim |>
  dplyr::filter(abs(time - tau) < 1e-8) |>
  dplyr::mutate(dose = arms$dose[match(CRCL, arms$CRCL)])

ggplot(trough, aes(x = arm, y = sim)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = chao$q25, ymax = chao$q75,
           fill = "grey70", alpha = 0.45) +
  geom_hline(yintercept = chao$median, linetype = "dashed") +
  geom_hline(yintercept = lloq, linetype = "dotdash") +
  geom_boxplot(width = 0.5, outlier.size = 0.6) +
  scale_y_log10() +
  labs(x = NULL, y = "Predicted trough concentration (ng/mL)",
       title = "Figure 3 - steady-state trough by creatinine clearance",
       caption = paste("Replicates Figure 3 of Ueshima 2025. Grey band and dashed line:",
                       "observed IQR and median in Asians (Chao et al.).",
                       "Dot-dash line: assay LLOQ 1 ng/mL.")) +
  theme_bw()
```

![Replicates Figure 3 of Ueshima
2025.](Ueshima_2025_edoxaban_files/figure-html/figure-3-1.png)

Replicates Figure 3 of Ueshima 2025.

``` r

fig3 <- trough |>
  dplyr::group_by(arm, CRCL, dose) |>
  dplyr::summarise(
    median_sim = median(sim),
    q25_sim    = quantile(sim, 0.25),
    q75_sim    = quantile(sim, 0.75),
    .groups    = "drop"
  ) |>
  dplyr::mutate(
    typical  = closed_form(tau, dose, CRCL),
    pct_diff = 100 * (median_sim - typical) / typical,
    iqr_over_median = (q75_sim - q25_sim) / median_sim
  )
knitr::kable(
  fig3 |>
    dplyr::mutate(dplyr::across(median_sim:iqr_over_median, \(x) round(x, 3))) |>
    dplyr::rename(
      "Arm"                        = arm,
      "CLcr (mL/min)"              = CRCL,
      "Dose (mg)"                  = dose,
      "Simulated median (ng/mL)"   = median_sim,
      "Simulated Q1 (ng/mL)"       = q25_sim,
      "Simulated Q3 (ng/mL)"       = q75_sim,
      "Typical value (ng/mL)"      = typical,
      "% diff (median vs typical)" = pct_diff,
      "IQR / median"               = iqr_over_median
    ),
  caption = "Simulated steady-state trough concentrations and the analytical typical value."
)
```

| Arm | CLcr (mL/min) | Dose (mg) | Simulated median (ng/mL) | Simulated Q1 (ng/mL) | Simulated Q3 (ng/mL) | Typical value (ng/mL) | % diff (median vs typical) | IQR / median |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| CLcr 30, 30 mg | 30 | 30 | 36.877 | 21.113 | 61.606 | 37.350 | -1.267 | 1.098 |
| CLcr 60, 60 mg | 60 | 60 | 30.571 | 19.777 | 50.745 | 28.866 | 5.909 | 1.013 |
| CLcr 90, 60 mg | 90 | 60 | 15.585 | 6.734 | 28.058 | 14.162 | 10.046 | 1.368 |
| CLcr 120, 60 mg | 120 | 60 | 7.470 | 3.614 | 17.758 | 7.703 | -3.027 | 1.893 |

Simulated steady-state trough concentrations and the analytical typical
value. {.table}

``` r


# The cohort median estimates the typical value; assert on the CENTRE, with
# headroom for which subjects a given rxode2 build / thread count happens to
# draw. 200 subjects with a total log-scale SD of ~0.86 give a median standard
# error of about 8%.
stopifnot(max(abs(fig3$pct_diff)) < 25)
```

### Gate: the two variability terms reach the simulation at the encoded scale

Reading `omega1 (CV%) = 26.4` and `sigma (CV%) = 58.7` as standard
deviations rather than as lognormal coefficients of variation (or,
worse, as variances) is the one transcription decision in Table 2 that a
structural check cannot catch. Three checks pin it down. Note that the
trough is a strongly non-linear function of clearance –
`d log C(24) / d log CL = -24k / (1 - exp(-24k))`, about -2.3 for the 60
mL/min arm – so the trough spread is roughly 2.4x the clearance spread
and cannot be compared to `omega1` directly. Test each term where it
enters instead.

``` r

# A. Every individual profile is the closed-form solution evaluated at that
#    individual's own cl and vc. This is deterministic per subject, so the
#    bound is numerical.
ind <- sim |>
  dplyr::filter(abs(time - tau) < 1e-8) |>
  dplyr::mutate(
    dose     = arms$dose[match(CRCL, arms$CRCL)],
    k_i      = cl / vc,
    analytic = (dose / vc) * 1000 * exp(-k_i * tau) / (1 - exp(-k_i * tau)),
    rel      = Cc / analytic - 1
  )
cat(sprintf("max |Cc / individual closed form - 1| over %d subjects: %.3g\n",
            nrow(ind), max(abs(ind$rel))))
#> max |Cc / individual closed form - 1| over 800 subjects: 1.11e-06
stopifnot(max(abs(ind$rel)) < 1e-4)

# B. IIV enters on cl only, so sd(log(cl)) recovers omega1 = 0.264 directly.
sd_log_cl <- sd(log(ind$cl / (28.2 * (ind$CRCL / 61.8)^0.692)))
cat(sprintf("sd(log eta multiplier on CL/F): %.4f (encoded omega1 = 0.264)\n", sd_log_cl))
#> sd(log eta multiplier on CL/F): 0.2659 (encoded omega1 = 0.264)
stopifnot(abs(sd_log_cl / 0.264 - 1) < 0.15)

# C. `sim` is `Cc` times the log-normal residual, so their log ratio recovers
#    sigma = 0.587 directly.
sd_log_eps <- sd(log(ind$sim / ind$Cc))
cat(sprintf("sd(log(sim / Cc)): %.4f (encoded sigma = 0.587)\n", sd_log_eps))
#> sd(log(sim / Cc)): 0.5702 (encoded sigma = 0.587)
stopifnot(abs(sd_log_eps / 0.587 - 1) < 0.15)
```

Ueshima 2025 (Discussion, final paragraph) reports that the predicted
troughs at CLcr 90 and 120 mL/min “tended to be lower than the median of
previously observed concentrations”, while the lower-clearance arms sat
within the observed range. Both halves of that claim are reproduced:

``` r

claim <- fig3 |>
  dplyr::mutate(
    within_chao_iqr = median_sim > chao$q25 & median_sim < chao$q75,
    below_chao_med  = median_sim < chao$median
  )
knitr::kable(
  claim |>
    dplyr::select(arm, median_sim, within_chao_iqr, below_chao_med) |>
    dplyr::mutate(median_sim = round(median_sim, 1)) |>
    dplyr::rename(
      "Arm"                       = arm,
      "Simulated median (ng/mL)"  = median_sim,
      "Within Chao IQR 13.9-47.4" = within_chao_iqr,
      "Below Chao median 25.0"    = below_chao_med
    ),
  caption = "Ueshima 2025 Discussion claim about Figure 3, reproduced."
)
```

| Arm | Simulated median (ng/mL) | Within Chao IQR 13.9-47.4 | Below Chao median 25.0 |
|:---|---:|:---|:---|
| CLcr 30, 30 mg | 36.9 | TRUE | FALSE |
| CLcr 60, 60 mg | 30.6 | TRUE | FALSE |
| CLcr 90, 60 mg | 15.6 | TRUE | TRUE |
| CLcr 120, 60 mg | 7.5 | FALSE | TRUE |

Ueshima 2025 Discussion claim about Figure 3, reproduced. {.table}

``` r

stopifnot(
  claim$within_chao_iqr[claim$CRCL %in% c(30, 60)],
  claim$below_chao_med[claim$CRCL %in% c(90, 120)]
)
```

## PKNCA validation

Steady-state non-compartmental analysis over the 24 h dosing interval.
The model is dosed directly into the central compartment, so the
time-zero record is the steady-state peak, not a pre-dose zero – no
synthetic `Cc = 0` row is added, and its presence is asserted rather
than assumed.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, arm)

# Every subject must already carry a record at the interval start; PKNCA
# otherwise silently anchors AUC at the first available time.
stopifnot(
  all(tapply(sim_nca$time, sim_nca$id, \(x) any(abs(x) < 1e-8))),
  all(tapply(sim_nca$time, sim_nca$id, \(x) any(abs(x - tau) < 1e-8)))
)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | arm + id)

dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, arm)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id)

intervals <- data.frame(
  start     = 0,
  end       = tau,
  cmax      = TRUE,
  tmax      = TRUE,
  cmin      = TRUE,
  auclast   = TRUE,
  cav       = TRUE,
  ctrough      = TRUE,
  half.life = TRUE
)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
```

### Comparison against the paper’s printed equations

Ueshima 2025 publishes no NCA table – it is a sparse-sampling
therapeutic-drug-monitoring study with one concentration per patient.
The reference column below is therefore the **analytical steady-state
solution of the paper’s own printed equations** (Eq. 3 and Eq. 10),
which is the thing this vignette exists to check: does the packaged
encoding reproduce what the paper printed? For a one-compartment model
dosed into the central compartment,

- `AUC0-24,ss` = 1000 x Dose / (CL/F),
- `Cmax,ss` = 1000 x (Dose / (Vd/F)) / (1 - exp(-k x 24)) at `Tmax = 0`,
- `Ctrough` = `Cmin` = `Cmax,ss` x exp(-k x 24),
- `Cav` = `AUC0-24,ss` / 24, and
- `t1/2` = log(2) / k, with `k` = (CL/F) / (Vd/F).

The comparison uses the typical-value (zero random effects) simulation,
so both sides are deterministic and the agreement should be
numerical-precision tight.

``` r

nca_typ_conc <- sim0 |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, arm)
nca_typ <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(nca_typ_conc, Cc ~ time | arm + id),
  PKNCA::PKNCAdose(
    events |>
      dplyr::filter(evid == 1, id %in% nca_typ_conc$id) |>
      dplyr::select(id, time, amt, arm),
    amt ~ time | arm + id
  ),
  intervals = intervals
))

published_eq <- arms |>
  dplyr::mutate(
    k         = cl_paper(CRCL) / v_paper(70),
    auclast   = 1000 * dose / cl_paper(CRCL),
    cmax      = 1000 * (dose / v_paper(70)) / (1 - exp(-k * tau)),
    tmax      = 0,
    cmin      = cmax * exp(-k * tau),
    ctrough      = cmin,
    cav       = auclast / tau,
    half.life = log(2) / k
  ) |>
  dplyr::select(arm, cmax, tmax, cmin, ctrough, cav, auclast, half.life)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_typ,
  reference     = published_eq,
  by            = "arm",
  units         = c(cmax = "ng/mL", cmin = "ng/mL", ctrough = "ng/mL",
                    cav = "ng/mL", auclast = "ng*h/mL",
                    tmax = "h", half.life = "h"),
  tolerance_pct = 20
)
knitr::kable(
  cmp,
  caption = paste("Typical-value steady-state NCA vs. the analytical solution of",
                  "Ueshima 2025 Eqs. 3 and 10. * differs by >20%."),
  align   = c("l", "l", "r", "r", "r")
)
```

| NCA parameter      | arm             | Reference | Simulated | % diff |
|:-------------------|:----------------|----------:|----------:|-------:|
| Cmax (ng/mL)       | CLcr 30, 30 mg  |       127 |       127 |  +0.0% |
| Cmax (ng/mL)       | CLcr 60, 60 mg  |       207 |       207 |  +0.0% |
| Cmax (ng/mL)       | CLcr 90, 60 mg  |       193 |       193 |  +0.0% |
| Cmax (ng/mL)       | CLcr 120, 60 mg |       186 |       186 |  +0.0% |
| Cmin (ng/mL)       | CLcr 30, 30 mg  |      37.4 |      37.4 |  -0.0% |
| Cmin (ng/mL)       | CLcr 60, 60 mg  |      28.9 |      28.9 |  -0.0% |
| Cmin (ng/mL)       | CLcr 90, 60 mg  |      14.2 |      14.2 |  -0.0% |
| Cmin (ng/mL)       | CLcr 120, 60 mg |       7.7 |       7.7 |  -0.0% |
| Tmax (h)           | CLcr 30, 30 mg  |         0 |         0 |      — |
| Tmax (h)           | CLcr 60, 60 mg  |         0 |         0 |      — |
| Tmax (h)           | CLcr 90, 60 mg  |         0 |         0 |      — |
| Tmax (h)           | CLcr 120, 60 mg |         0 |         0 |      — |
| AUClast (ng\*h/mL) | CLcr 30, 30 mg  |      1750 |      1750 |  -0.0% |
| AUClast (ng\*h/mL) | CLcr 60, 60 mg  |      2170 |      2170 |  -0.0% |
| AUClast (ng\*h/mL) | CLcr 90, 60 mg  |      1640 |      1640 |  -0.0% |
| AUClast (ng\*h/mL) | CLcr 120, 60 mg |      1340 |      1340 |  -0.0% |
| t½ (h)             | CLcr 30, 30 mg  |      13.6 |      13.6 |  -0.0% |
| t½ (h)             | CLcr 60, 60 mg  |      8.44 |      8.44 |  -0.0% |
| t½ (h)             | CLcr 90, 60 mg  |      6.37 |      6.37 |  +0.0% |
| t½ (h)             | CLcr 120, 60 mg |      5.22 |      5.22 |  +0.0% |
| Cavg (ng/mL)       | CLcr 30, 30 mg  |      73.1 |      73.1 |  -0.0% |
| Cavg (ng/mL)       | CLcr 60, 60 mg  |      90.5 |      90.5 |  -0.0% |
| Cavg (ng/mL)       | CLcr 90, 60 mg  |      68.3 |      68.3 |  -0.0% |
| Cavg (ng/mL)       | CLcr 120, 60 mg |        56 |        56 |  -0.0% |
| Ctrough (ng/mL)    | CLcr 30, 30 mg  |      37.4 |      37.4 |  -0.0% |
| Ctrough (ng/mL)    | CLcr 60, 60 mg  |      28.9 |      28.9 |  -0.0% |
| Ctrough (ng/mL)    | CLcr 90, 60 mg  |      14.2 |      14.2 |  -0.0% |
| Ctrough (ng/mL)    | CLcr 120, 60 mg |       7.7 |       7.7 |  -0.0% |

Typical-value steady-state NCA vs. the analytical solution of Ueshima
2025 Eqs. 3 and 10. \* differs by \>20%. {.table}

``` r


# `% diff` from ncaComparisonTable() is FORMATTED CHARACTER, so assert on
# numbers recomputed here rather than parsing the display column.
typ_wide <- as.data.frame(nca_typ) |>
  dplyr::select(arm, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)
chk <- published_eq |>
  dplyr::left_join(typ_wide, by = "arm", suffix = c("_ref", "_sim"))

# Drive the comparison off an explicit parameter list and demand that both
# columns exist for every one of them. A silently-missing column would
# otherwise drop that parameter from `worst` and leave a gate that cannot go
# red.
compare_params <- c("cmax", "cmin", "ctrough", "cav", "auclast", "half.life")
pct_by_param <- vapply(compare_params, function(p) {
  ref <- chk[[paste0(p, "_ref")]]
  sm  <- chk[[paste0(p, "_sim")]]
  if (is.null(ref) || is.null(sm) || anyNA(ref) || anyNA(sm)) {
    stop("missing or NA reference/simulated column for '", p, "'")
  }
  max(abs(100 * (sm - ref) / ref))
}, numeric(1))
print(round(pct_by_param, 4))
#>      cmax      cmin   ctrough       cav   auclast half.life 
#>     0e+00     1e-04     1e-04     1e-04     1e-04     0e+00
worst <- max(pct_by_param)
cat(sprintf("worst |%% diff| across %s: %.4f%%\n",
            paste(compare_params, collapse = ", "), worst))
#> worst |% diff| across cmax, cmin, ctrough, cav, auclast, half.life: 0.0001%
stopifnot(worst < 0.5, all(chk$tmax_sim == 0))
```

### Cohort-level NCA

The same NCA over the full 200-subject-per-arm cohort. Between-subject
variability is on clearance only, so exposure metrics vary while `Tmax`
and `Cmax` per unit volume do not.

``` r

nca_summary <- as.data.frame(nca_res) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "ctrough", "auclast", "half.life")) |>
  dplyr::group_by(arm, PPTESTCD) |>
  dplyr::summarise(
    median = median(PPORRES),
    q05    = quantile(PPORRES, 0.05),
    q95    = quantile(PPORRES, 0.95),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Parameter = nlmixr2lib::ncaParamLabel(
      PPTESTCD,
      units = c(cmax = "ng/mL", ctrough = "ng/mL",
                auclast = "ng*h/mL", half.life = "h")
    )
  ) |>
  dplyr::select(arm, Parameter, median, q05, q95)

knitr::kable(
  nca_summary |>
    dplyr::mutate(dplyr::across(median:q95, \(x) signif(x, 4))) |>
    dplyr::rename("Arm" = arm, "Median" = median,
                  "5th percentile" = q05, "95th percentile" = q95),
  caption = "Cohort steady-state NCA, median and 90% interval across subjects."
)
```

| Arm             | Parameter          |   Median | 5th percentile | 95th percentile |
|:----------------|:-------------------|---------:|---------------:|----------------:|
| CLcr 30, 30 mg  | AUClast (ng\*h/mL) | 1822.000 |       1120.000 |        2602.000 |
| CLcr 30, 30 mg  | Cmax (ng/mL)       |  129.000 |        104.700 |         159.100 |
| CLcr 30, 30 mg  | Ctrough (ng/mL)    |   39.850 |         15.470 |          69.890 |
| CLcr 30, 30 mg  | t½ (h)             |   14.160 |          8.702 |          20.230 |
| CLcr 60, 60 mg  | AUClast (ng\*h/mL) | 2191.000 |       1425.000 |        3184.000 |
| CLcr 60, 60 mg  | Cmax (ng/mL)       |  207.800 |        187.700 |         241.300 |
| CLcr 60, 60 mg  | Ctrough (ng/mL)    |   29.460 |          9.303 |          62.900 |
| CLcr 60, 60 mg  | t½ (h)             |    8.515 |          5.537 |          12.370 |
| CLcr 90, 60 mg  | AUClast (ng\*h/mL) | 1681.000 |       1029.000 |        2519.000 |
| CLcr 90, 60 mg  | Cmax (ng/mL)       |  193.500 |        181.200 |         218.200 |
| CLcr 90, 60 mg  | Ctrough (ng/mL)    |   15.170 |          2.827 |          39.890 |
| CLcr 90, 60 mg  | t½ (h)             |    6.534 |          3.998 |           9.788 |
| CLcr 120, 60 mg | AUClast (ng\*h/mL) | 1346.000 |        879.900 |        2231.000 |
| CLcr 120, 60 mg | Cmax (ng/mL)       |  186.100 |        179.700 |         209.100 |
| CLcr 120, 60 mg | Ctrough (ng/mL)    |    7.735 |          1.386 |          30.700 |
| CLcr 120, 60 mg | t½ (h)             |    5.230 |          3.419 |           8.672 |

Cohort steady-state NCA, median and 90% interval across subjects.
{.table}

``` r


# The cohort median AUC estimates 1000 * Dose / CL/F. Assert on the centre.
auc_med <- nca_summary |>
  dplyr::filter(grepl("^AUC", Parameter)) |>
  dplyr::left_join(published_eq |> dplyr::select(arm, auclast), by = "arm") |>
  dplyr::mutate(pct = 100 * (median - auclast) / auclast)
cat(sprintf("AUC0-24 cohort median vs analytical, %% diff by arm: %s\n",
            paste(round(auc_med$pct, 1), collapse = ", ")))
#> AUC0-24 cohort median vs analytical, % diff by arm: 3.8, 0.9, 2.5, 0.1
stopifnot(max(abs(auc_med$pct)) < 25)
```

## Replicate Figure 2 – concentration vs. time after last dose

Ueshima 2025 Figure 2 is a prediction-corrected VPC of the observed
concentrations against time after the last dose, with the 5th, 50th and
95th percentiles of 1,000 simulated datasets. The panel below simulates
the study cohort itself – 131 patients at the observed dose split, with
body weight and creatinine clearance drawn to match the Table 1 medians
and ranges – and plots the same three percentiles.

``` r

n_study <- 131L
# Log-normal draws centred on the Table 1 medians, truncated to the Table 1
# ranges. The paper reports only median and range, so the dispersion is an
# assumption; see Assumptions and deviations.
rtrunc_lnorm <- function(n, med, lo, hi, sdlog) {
  x <- rep(NA_real_, n)
  todo <- seq_len(n)
  while (length(todo)) {
    cand <- med * exp(stats::rnorm(length(todo), 0, sdlog))
    ok <- cand >= lo & cand <= hi
    x[todo[ok]] <- cand[ok]
    todo <- todo[!ok]
  }
  x
}

study <- tibble::tibble(
  id   = seq_len(n_study),
  WT   = rtrunc_lnorm(n_study, 61.4, 37.3, 97.3, 0.17),
  CRCL = rtrunc_lnorm(n_study, 61.8, 20.0, 135.5, 0.38),
  amt  = rep(c(15, 30, 60), times = c(3L, 70L, 58L))   # Ueshima 2025 Table 1
)

study_ev <- dplyr::bind_rows(
  study |> dplyr::mutate(time = 0, evid = 1L, cmt = "central", ii = tau, ss = 1L),
  tidyr::crossing(id = study$id, time = seq(0, 48, by = 0.5)) |>
    dplyr::left_join(study |> dplyr::select(id, WT, CRCL), by = "id") |>
    dplyr::mutate(amt = NA_real_, evid = 0L, cmt = "central", ii = 0, ss = 0L)
) |>
  dplyr::arrange(id, time, dplyr::desc(evid))

# Observations beyond 24 h describe the profile after the LAST dose, i.e. a
# missed / delayed dose -- which is exactly the situation Figure 2's 24-47 h
# samples represent -- so dosing is not continued past time 0.
study_sim <- rxode2::rxSolve(mod, events = as.data.frame(study_ev),
                             keep = c("WT", "CRCL")) |>
  as.data.frame()

study_pct <- study_sim |>
  dplyr::group_by(time) |>
  dplyr::summarise(
    Q05 = quantile(sim, 0.05), Q50 = quantile(sim, 0.50),
    Q95 = quantile(sim, 0.95), .groups = "drop"
  )

ggplot(study_pct, aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = lloq, linetype = "dotdash") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 48, by = 12)) +
  labs(x = "Time after last dose (h)", y = "Edoxaban concentration (ng/mL)",
       title = "Figure 2 - simulated 5th / 50th / 95th percentiles",
       caption = paste("Replicates Figure 2 of Ueshima 2025.",
                       "Dot-dash line: assay LLOQ 1 ng/mL.")) +
  theme_bw()
```

![Replicates Figure 2 of Ueshima
2025.](Ueshima_2025_edoxaban_files/figure-html/figure-2-1.png)

Replicates Figure 2 of Ueshima 2025.

``` r

# Ueshima 2025 Table 1: the 131 observed concentrations, sampled 9-47 h after
# the last dose, had a median of 19.5 ng/mL and a range of 1.7-152.0.
#
# This CANNOT be turned into a point comparison of medians. The observed median
# is taken over the sampling times that actually occurred, and the paper does
# not report their distribution; Figure 2 shows them clustered near 24 h with a
# tail out to 48 h and a group near 12 h. The grid below weights 9-47 h
# uniformly, which puts more mass on the late, low-concentration end than the
# real sampling did, so the simulated median sits below the observed one by
# construction. What the simulation can legitimately be asked is whether the
# published summary statistics fall inside the distribution it produces.
window <- study_sim |> dplyr::filter(time >= 9, time <= 47) |> dplyr::pull(sim)
window_q <- quantile(window, c(0.05, 0.25, 0.50, 0.75, 0.95))
print(round(window_q, 1))
#>    5%   25%   50%   75%   95% 
#>   0.4   3.7  12.8  34.1 107.1
cat(sprintf("simulated median over the 9-47 h window: %.1f ng/mL (observed 19.5; %.1f-%.1f is the simulated IQR)\n",
            window_q[["50%"]], window_q[["25%"]], window_q[["75%"]]))
#> simulated median over the 9-47 h window: 12.8 ng/mL (observed 19.5; 3.7-34.1 is the simulated IQR)

# Containment checks, all robust to which cohort a given machine draws:
stopifnot(
  # The observed median lies inside the simulated central 90%.
  19.5 > window_q[["5%"]], 19.5 < window_q[["95%"]],
  # The simulated distribution spans the observed range rather than sitting
  # entirely inside or entirely outside it.
  min(window) < 1.7, max(window) > 152.0
)

# Trough-time comparison, which removes the sampling-time weighting problem:
# the median 24 h concentration should be the same order as the observed
# median, since 24 h is where the paper's sampling density peaks.
med24 <- median(study_sim$sim[abs(study_sim$time - tau) < 1e-8])
cat(sprintf("simulated median at 24 h: %.1f ng/mL\n", med24))
#> simulated median at 24 h: 15.5 ng/mL
stopifnot(med24 > 19.5 / 3, med24 < 19.5 * 3)
```

## Assumptions and deviations

- **Body weight for the Figure 3 replication.** Ueshima 2025 does not
  state the body weight used in its Monte Carlo simulation. 70 kg is
  used here, which is the model’s own reference weight and therefore
  gives Vd/F = 336.4 L exactly; it reproduces the Figure 3 medians. The
  cohort median weight is 61.4 kg, which would lower Vd/F to 295 L and
  raise the trough by about 20%.
- **Covariate distributions for the Figure 2 replication.** Table 1
  reports only medians and ranges for body weight and creatinine
  clearance, so the log-normal dispersions (`sdlog` 0.17 and 0.38) are
  chosen to fill the published ranges and are an assumption, not a
  published quantity. Weight and creatinine clearance are drawn
  independently even though the Cockcroft-Gault equation makes them
  correlated in reality; this widens the simulated spread slightly.
- **No non-renal clearance arm.** The paper fitted CL/F as the sum of an
  apparent renal and an apparent non-renal component (Eq. 4), but the
  non-renal population mean converged to 0.01 L/h and was judged
  negligible, so the published final model (Eq. 10) is the single power
  term encoded here. The paper itself notes this is likely an
  underestimate – an intravenous study puts renal clearance at 49.1% of
  total – because no urine or metabolite samples were collected. A user
  applying this model outside the fitted creatinine-clearance range of
  20.0-135.5 mL/min will therefore extrapolate a term that carries the
  whole of clearance.
- **Vd/F and the residual error are not this paper’s estimates.** Both
  were fixed to values from Krekels 2016 (ENGAGE AF-TIMI 48, reference
  \[22\]): Vd/F = 336.4 L (193 L x 1.743 Asian fold-change) and sigma =
  58.7 CV% (14.6% healthy-subject residual + 44.1% AF-patient
  increment). They are encoded with `fixed()` and are recorded as such
  in the model file. The weight exponent on Vd/F is likewise `fixed(1)`,
  not an allometric 0.75 or 1 chosen by this analysis – the paper fixed
  it “due to a lack of individual data”.
- **No absorption, and therefore no `Tmax`.** The model dose enters
  `central` directly (NONMEM ADVAN1 TRANS2), so the simulated `Tmax` is
  0 and `Cmax` is the immediate post-dose concentration. This is a
  property of the published model, not of the encoding: with only one
  sample per patient on the elimination phase, no absorption parameter
  was identifiable. Simulated peak concentrations should not be compared
  to observed edoxaban `Cmax`, which occurs 1-2 h post-dose in richly
  sampled studies.
- **Pharmacogenomic covariates are documented, not implemented.**
  `CYP3A5*3` and the three `ABCB1` polymorphisms, the 2677/3435
  haplotype, sex, age, AST, ALT and concomitant CYP3A4 / P-glycoprotein
  inhibitor use were all screened and none was retained. They appear in
  the model file’s `covariatesDataExcluded` metadata and are absent from
  `model()`. The paper flags the small number of exposed patients (2
  with `CYP3A5*1/*1`, 16 on an inhibitor) as a limitation of the
  negative co-medication and genotype findings.
- **Figure 3 arms are simulated at steady state via `ss = 1`** rather
  than by a multi-day burn-in, which gives the exact steady-state
  solution rather than an approach to it.
