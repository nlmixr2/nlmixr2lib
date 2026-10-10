# Wuzhi-capsule schizandrol lignans as CYP3A4/CYP3A5 inhibitors (He 2022)

## Model and source

- Citation: He Q, Bu F, Wang Q, Li M, Lin J, Tang Z, Mak WY, Zhuang X,
  Zhu X, Lin HS, Xiang X. Examination of the Impact of CYP3A4/5 on
  Drug-Drug Interaction between Schizandrol A/Schizandrol B and
  Tacrolimus (FK-506): A Physiologically Based Pharmacokinetic Modeling
  Approach. Int J Mol Sci. 2022;23(9):4485. <doi:10.3390/ijms23094485>.
  PMCID: PMC9103789. Reversible inhibition constants (Ki): Results 2.2
  and Figure 2C,D (SZB on CYP3A4 Ki = 2.18 uM, on CYP3A5 Ki = 2.03 uM);
  pooled CYP3A Ki = 5.82 uM (Figure 2B). NOTE: the Discussion text swaps
  the CYP3A4 and CYP3A5 Ki values (it prints CYP3A4 Ki = 2.03, CYP3A5 Ki
  = 2.18); Results 2.2 and the Figure 2C,D annotations are used here
  (CYP3A4 = 2.18, CYP3A5 = 2.03). Time-dependent inactivation constants:
  Results 2.3, Discussion, and the figure annotations (CYP3A4 kinact =
  0.37 /min, KI = 0.69 uM, Figure 5B; CYP3A5 kinact = 0.009 /min, KI =
  0.5 uM, Figure 5D; pooled CYP3A kinact = 0.044 /min, KI = 0.43 uM,
  Figure 3D). IC50 shift 21.39 (IC50 11.98 uM no preincubation, 0.56 uM
  after NADPH preincubation): Results 2.1. Michaelis-Menten inactivation
  form kobs = kinact\*I/(KI+I): Materials and Methods 4.4 and Equation
  (3).
- Article: <https://doi.org/10.3390/ijms23094485>
- PubMed Central open-access copy:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9103789/>

He 2022 co-administer the Wuzhi capsule (a standardised ethanolic
extract of *Schisandra sphenanthera*) with tacrolimus (FK-506), a
combination used in Chinese transplant practice because the extract’s
lignans inhibit CYP3A and raise tacrolimus exposure, cutting the
required tacrolimus dose. An earlier paper by the same group (He 2021,
PMC7997453) measured the two most abundant lignans, schisantherin A and
schisandrin A. This paper measures the other two constituents,
**schizandrol A (SZA)** and **schizandrol B (SZB)**: their reversible
(RI) and time-dependent (TDI) inhibition of CYP3A4 and CYP3A5
separately, using CYP3A5-genotyped human liver microsomes with CYP3cide
to isolate CYP3A5.

This paper contributes two model files, one per lignan:

``` r

mod_sza <- rxode2::rxode2(readModelDb("He_2022_schizandrolA_cyp3a"))
mod_szb <- rxode2::rxode2(readModelDb("He_2022_schizandrolB_cyp3a"))

# Fixed and estimated ini() values are not returned as columns by rxSolve(), so
# read them straight off the packaged model when a check needs one.
ini_val <- function(mod, nm) {
  stopifnot(nm %in% mod$iniDf$name)
  mod$iniDf$est[mod$iniDf$name == nm]
}
```

## Scope: what this paper does and does not contribute

He 2022 has two halves, and only one of them is reproducible.

The **first half is the authors’ own bench work**: reversible-inhibition
(RI) and time-dependent-inhibition (TDI) assays that measure the
inhibitory potency of SZA and SZB on CYP3A4 and CYP3A5, each isoform in
its own microsomal system, plus a pooled-HLM CYP3A screen. Every fitted
constant is printed in the Results and Discussion text and annotated on
Figures 2-5. That is the layer this vignette validates and that the two
model files carry.

The **second half is a Simcyp (version 16) whole-body physiologically
based model** that integrates those constants (together with the
physicochemical and disposition parameters of Table 3) to predict the
fold-increase in tacrolimus AUC in CYP3A5 expressers and non-expressers
(Table 2). That layer is a vendor platform model and is **not**
extracted: the whole-body ODE system is Simcyp’s rather than the
authors’, the tissue-partition (Kp) values come from a built-in method,
the organ volumes/flows and enzyme abundances are Simcyp database
outputs, and no project file is deposited. Reconstructing it would
require substituting physiology from outside the source. The section
“Published results not reproduced here” records what that half produced.

A load-bearing consequence for the in vitro layer: the two lignans and
the two isoforms are inhibited by different combinations of the two
mechanisms, and each model file records exactly the combination the
paper fitted.

| Lignan | Isoform | Reversible (Ki, uM) | Time-dependent (kinact /min, KI uM) |
|:---|:---|:---|:---|
| Schizandrol A (SZA) | CYP3A4 | little (not fitted) | 0.024, 15.38 |
| Schizandrol A (SZA) | CYP3A5 | little (not fitted) | none (Figure 4C) |
| Schizandrol B (SZB) | CYP3A4 | 2.18 | 0.37, 0.69 |
| Schizandrol B (SZB) | CYP3A5 | 2.03 | 0.009, 0.5 |

Isoform-resolved inhibition mechanisms fitted by He 2022 (Results 2.2,
2.3; Figures 2, 4, 5). A pooled-HLM CYP3A screen (Figures 2B, 3) is also
carried in each file. {.table}

## In vitro systems

| Assay | CYP3A4 system | CYP3A5 system | Pooled CYP3A | Probe |
|:---|:---|:---|:---|:---|
| Reversible inhibition (RI) | CYP3A5*3/*3 HLM, 0.2 mg/mL | CYP3A5*1/*3 HLM + 1.2 uM CYP3cide, 0.2 mg/mL | Pooled HLM (22 donors), 0.2 mg/mL | Tacrolimus disappearance (Dixon plot) |
| Time-dependent inhibition (TDI) | CYP3A5*3/*3 HLM, 0.5 mg/mL | CYP3A5*1/*3 HLM + CYP3cide, 0.5 mg/mL | Pooled HLM (22 donors), 0.5 mg/mL | Testosterone 6-beta-hydroxylation (double-reciprocal plot) |

In vitro systems of He 2022 (Methods 4.3, 4.4). {.table}

CYP3cide (PF-04981517) is a selective mechanism-based inactivator of
CYP3A4; it is used to knock out CYP3A4 so that the residual testosterone
or tacrolimus turnover in the CYP3A5*1/*3 microsomes is attributable to
CYP3A5.

## Source trace

Every value in the two model files, with the location it came from.

| Model | Quantity | Value | Source location |
|:---|:---|:---|:---|
| He_2022_schizandrolA_cyp3a | ki_inact_3a4 | 15.38 uM | Figure 4B annotation; Results 2.3 |
| He_2022_schizandrolA_cyp3a | kinact_3a4 | 0.024 /min | Figure 4B annotation; Results 2.3 |
| He_2022_schizandrolA_cyp3a | ki_inact_3a (pooled) | 15.625 uM | Figure 3B annotation; Results 2.3 |
| He_2022_schizandrolA_cyp3a | kinact_3a (pooled) | 0.029 /min | Figure 3B annotation; Results 2.3 |
| He_2022_schizandrolA_cyp3a | (CYP3A5 TDI) | none | Figure 4C (no inactivation observed) |
| He_2022_schizandrolA_cyp3a | (any RI) | little | Figure 2A (no significant reversible inhibition) |
| He_2022_schizandrolB_cyp3a | ki_3a4 | 2.18 uM | Figure 2C annotation; Results 2.2 |
| He_2022_schizandrolB_cyp3a | ki_3a5 | 2.03 uM | Figure 2D annotation; Results 2.2 |
| He_2022_schizandrolB_cyp3a | ki_3a (pooled) | 5.82 uM | Figure 2B annotation; Results 2.2 |
| He_2022_schizandrolB_cyp3a | ki_inact_3a4 | 0.69 uM | Figure 5B annotation; Results 2.3 |
| He_2022_schizandrolB_cyp3a | kinact_3a4 | 0.37 /min | Figure 5B annotation; Results 2.3 |
| He_2022_schizandrolB_cyp3a | ki_inact_3a5 | 0.5 uM | Figure 5D annotation; Results 2.3 |
| He_2022_schizandrolB_cyp3a | kinact_3a5 | 0.009 /min | Figure 5D annotation; Results 2.3 |
| He_2022_schizandrolB_cyp3a | ki_inact_3a (pooled) | 0.43 uM | Figure 3D annotation; Results 2.3 |
| He_2022_schizandrolB_cyp3a | kinact_3a (pooled) | 0.044 /min | Figure 3D annotation; Results 2.3 |
| both | addSd (fixed 0) | 0 | Not reported anywhere; see Assumptions |
| Inactivation-rate form | kobs = kinact\*I/(KI+I) | \- | Methods 4.4.3, Equation (3) |
| Inactivation ODE | d/dt(E) = -kobs\*E | \- | log-linear decay of Figures 3A, 4A, 5A, 5C |
| Reversible form | 1/(1 + I/Ki) | \- | Competitive inhibition, Dixon plots (Figure 2) |

Source trace for every ini() value and every model equation. {.table}

The fitted constants live in the annotations of Figures 2-5 and are
repeated in the Results and Discussion text; they are the authors’ own
reported estimates, not digitised from the plotted points. **One
transcription caveat:** the Discussion paragraph swaps the reversible
CYP3A4 and CYP3A5 Ki values for SZB (it prints CYP3A4 Ki = 2.03, CYP3A5
Ki = 2.18), whereas Results 2.2 and the Figure 2C/2D annotations give
CYP3A4 Ki = 2.18 and CYP3A5 Ki = 2.03. The figure/Results order is used
here.

## Structural identities of the inactivation model

The time-dependent inactivation is a first-order decay whose rate
constant `kobs` is a hyperbolic (Michaelis-Menten) function of the
inhibitor concentration. These identities are exact by construction and
pin the transcription of `kinact` and `KI`.

``` r

kobs <- function(mod, suffix, I) {
  kinact <- ini_val(mod, paste0("kinact_", suffix))
  KI     <- ini_val(mod, paste0("ki_inact_", suffix))
  kinact * I / (KI + I)
}

inactivation_systems <- list(
  list(mod_sza, "3a4"), list(mod_sza, "3a"),
  list(mod_szb, "3a4"), list(mod_szb, "3a5"), list(mod_szb, "3a")
)

# 1. kobs(0) = 0 exactly (no inhibitor => no inactivation).
for (s in inactivation_systems) stopifnot(kobs(s[[1]], s[[2]], 0) == 0)

# 2. kobs(KI) = kinact/2 exactly (half-maximal at I = KI).
for (s in inactivation_systems) {
  m <- s[[1]]; suf <- s[[2]]
  KI <- ini_val(m, paste0("ki_inact_", suf))
  kinact <- ini_val(m, paste0("kinact_", suf))
  stopifnot(abs(kobs(m, suf, KI) - kinact / 2) < 1e-12)
}

# 3. Plateau: kobs approaches kinact at saturating I.
for (s in inactivation_systems) {
  stopifnot(abs(kobs(s[[1]], s[[2]], 1e7) - ini_val(s[[1]], paste0("kinact_", s[[2]]))) < 1e-5)
}
```

## Known-answer check: the published inhibitory-efficiency ratios

The paper prints the inactivation efficiency `kinact / KI` for every
enzyme system (Discussion), in mL/min/umol. With `kinact` in /min and
`KI` in uM (= umol/L), `kinact / KI` is in L/min/umol,
i.e. `1000 * kinact / KI` is the printed mL/min/umol value. Reproducing
all five is a direct check that every `kinact` and `KI` pair was
transcribed correctly, because an error in either moves the ratio.

``` r

efficiency <- function(mod, suffix) {
  1000 * ini_val(mod, paste0("kinact_", suffix)) / ini_val(mod, paste0("ki_inact_", suffix))
}

eff_tbl <- tibble::tibble(
  System = c("SZA pooled CYP3A", "SZA CYP3A4",
             "SZB pooled CYP3A", "SZB CYP3A4", "SZB CYP3A5"),
  `Model kinact/KI (mL/min/umol)` = c(
    efficiency(mod_sza, "3a"), efficiency(mod_sza, "3a4"),
    efficiency(mod_szb, "3a"), efficiency(mod_szb, "3a4"), efficiency(mod_szb, "3a5")
  ),
  `Published (mL/min/umol)` = c(1.856, 1.56, 102.33, 536.2, 18)
)

# The published ratios are rounded; the model recomputes them to full precision.
stopifnot(with(eff_tbl, max(abs(`Model kinact/KI (mL/min/umol)` - `Published (mL/min/umol)`) /
  `Published (mL/min/umol)`) < 0.01))

eff_tbl |>
  dplyr::mutate(`Model kinact/KI (mL/min/umol)` = round(`Model kinact/KI (mL/min/umol)`, 3)) |>
  knitr::kable(caption = "Reproduction of the five published kinact/KI efficiency ratios (He 2022 Discussion).")
```

| System           | Model kinact/KI (mL/min/umol) | Published (mL/min/umol) |
|:-----------------|------------------------------:|------------------------:|
| SZA pooled CYP3A |                         1.856 |                   1.856 |
| SZA CYP3A4       |                         1.560 |                   1.560 |
| SZB pooled CYP3A |                       102.326 |                 102.330 |
| SZB CYP3A4       |                       536.232 |                 536.200 |
| SZB CYP3A5       |                        18.000 |                  18.000 |

Reproduction of the five published kinact/KI efficiency ratios (He 2022
Discussion). {.table}

SZB is one to two orders of magnitude more efficient a CYP3A inactivator
than SZA, and SZB is far more efficient on CYP3A4 (536) than on CYP3A5
(18) - the asymmetry the paper attributes to CYP3A5’s higher affinity
but lower turnover.

## IC50 shift (Results 2.1)

The IC50-shift screen (IC50 without preincubation divided by IC50 after
30-minute NADPH preincubation) is the paper’s first-pass flag for
time-dependent inhibition; a value above 1.5 is taken as suggestive.
These are reproduced from the printed IC50 values, not carried as model
parameters.

``` r

ic50_shift <- tibble::tibble(
  Lignan = c("SZA", "SZB"),
  `IC50 no preinc (uM)` = c(63.46, 11.98),
  `IC50 after preinc (uM)` = c(40.45, 0.56),
  `Shift` = c(63.46 / 40.45, 11.98 / 0.56),
  `Published shift` = c(1.57, 21.39)
)
stopifnot(with(ic50_shift, max(abs(Shift - `Published shift`)) < 0.05))
ic50_shift |>
  dplyr::mutate(Shift = round(Shift, 2)) |>
  knitr::kable(caption = "IC50 shift reproduced from the printed IC50 values (He 2022 Results 2.1).")
```

| Lignan | IC50 no preinc (uM) | IC50 after preinc (uM) | Shift | Published shift |
|:-------|--------------------:|-----------------------:|------:|----------------:|
| SZA    |               63.46 |                  40.45 |  1.57 |            1.57 |
| SZB    |               11.98 |                   0.56 | 21.39 |           21.39 |

IC50 shift reproduced from the printed IC50 values (He 2022 Results
2.1). {.table}

## Replicating the inactivation time courses (Figures 3A, 4A, 5A, 5C)

The ODE `d/dt(E) = -kobs * E` with `E(0) = 1` has the closed-form
solution `E(t) = exp(-kobs * t)`, which is the straight line the authors
plot on a log axis. The model is solved over the assay’s preincubation
grid at each measured inhibitor concentration.

``` r

solve_decay <- function(mod, isoform, conc, covname, tmax = 30) {
  ev <- do.call(rbind, lapply(seq_along(conc), function(i) {
    d <- data.frame(id = i, time = seq(0, tmax, by = 1), evid = 0L, cmt = paste0("enzyme_", isoform))
    d[[covname]] <- conc[i]
    d
  }))
  state <- paste0("enzyme_", isoform)
  out <- rxode2::rxSolve(mod, ev, returnType = "data.frame")
  out$conc <- conc[out$id]
  out$activity <- out[[state]]
  out
}

# SZA on CYP3A4 (Figure 4A): I = 0, 10, 16, 20, 25, 32 uM.
sza_3a4 <- solve_decay(mod_sza, "3a4", c(0, 10, 16, 20, 25, 32), "CP_SZA_UM")
#> Warning: multi-subject simulation without without 'omega'

# SZB on CYP3A4 (Figure 5A): I = 0, 0.5, 1, 2, 4, 8 uM.
szb_3a4 <- solve_decay(mod_szb, "3a4", c(0, 0.5, 1, 2, 4, 8), "CP_SZB_UM")
#> Warning: multi-subject simulation without without 'omega'
# SZB on CYP3A5 (Figure 5C): I = 0, 2, 4, 8, 16 uM.
szb_3a5 <- solve_decay(mod_szb, "3a5", c(0, 2, 4, 8, 16), "CP_SZB_UM")
#> Warning: multi-subject simulation without without 'omega'
```

``` r

# The ODE solution must equal exp(-kobs*t) to solver tolerance.
check_closed_form <- function(df, mod, suffix) {
  df$kobs <- kobs(mod, suffix, df$conc)
  df$analytic <- exp(-df$kobs * df$time)
  max(abs(df$activity - df$analytic))
}
stopifnot(check_closed_form(sza_3a4, mod_sza, "3a4") < 1e-5)
stopifnot(check_closed_form(szb_3a4, mod_szb, "3a4") < 1e-5)
stopifnot(check_closed_form(szb_3a5, mod_szb, "3a5") < 1e-5)
```

![Modelled log-percent CYP activity remaining during preincubation.
Replicates the linear-regression lines of He 2022 Figures 4A (SZA on
CYP3A4), 5A (SZB on CYP3A4) and 5C (SZB on CYP3A5). SZB's steeper CYP3A4
decay reflects its far larger
kinact/KI.](He_2022_tacrolimus_wuzhi_schizandrol_cyp3a_files/figure-html/decay_plot-1.png)

Modelled log-percent CYP activity remaining during preincubation.
Replicates the linear-regression lines of He 2022 Figures 4A (SZA on
CYP3A4), 5A (SZB on CYP3A4) and 5C (SZB on CYP3A5). SZB’s steeper CYP3A4
decay reflects its far larger kinact/KI.

SZA has no time-dependent inactivation of CYP3A5 (Figure 4C), so
`enzyme_3a5` in the SZA model holds at its baseline of 1 for all time;
that is confirmed directly:

``` r

sza_3a5 <- solve_decay(mod_sza, "3a5", c(0, 10, 20, 30, 40, 50), "CP_SZA_UM")
#> Warning: multi-subject simulation without without 'omega'
stopifnot(all(abs(sza_3a5$activity - 1) < 1e-8))
```

## Reversible-inhibition factors

Only SZB is a reversible inhibitor (SZA shows little reversible
inhibition, Figure 2A). The Dixon-plot inhibition constants enter a
victim-drug clearance as the competitive factor `1/(1 + I/Ki)`, exposed
by the SZB model as `riFactor_*`. At `I = Ki` the factor is exactly one
half. The two isoform Ki values are close (CYP3A4 2.18 uM, CYP3A5 2.03
uM), with the slightly smaller CYP3A5 Ki making it the marginally
stronger reversible target.

``` r

ri_factor <- function(I, Ki) 1 / (1 + I / Ki)

# Exact identity at I = Ki for both isoforms and the pooled constant.
for (nm in c("ki_3a4", "ki_3a5", "ki_3a")) {
  ki <- ini_val(mod_szb, nm)
  stopifnot(abs(ri_factor(ki, ki) - 0.5) < 1e-12)
}

# SZB reversibly inhibits CYP3A5 marginally more strongly than CYP3A4: the
# CYP3A5 Ki (2.03 uM) is smaller than the CYP3A4 Ki (2.18 uM), so at a common
# concentration the CYP3A5 activity loss is slightly larger.
I_common <- 2
cyp3a4_loss <- 1 - ri_factor(I_common, ini_val(mod_szb, "ki_3a4"))
cyp3a5_loss <- 1 - ri_factor(I_common, ini_val(mod_szb, "ki_3a5"))
stopifnot(cyp3a5_loss > cyp3a4_loss)
c(SZB_CYP3A4_activity_loss = cyp3a4_loss, SZB_CYP3A5_activity_loss = cyp3a5_loss)
#> SZB_CYP3A4_activity_loss SZB_CYP3A5_activity_loss 
#>                0.4784689                0.4962779
```

## Published results not reproduced here

The Simcyp whole-body half predicted the fold-increase (AUC ratio, AUCR)
in tacrolimus exposure when co-dosed with each lignan, split by CYP3A5
genotype and by which inhibition mechanism was switched on (Table 2). In
CYP3A5 non-expressers, multiple doses of SZB raised tacrolimus AUC by
57% (AUCR 1.57) under the combined RI+TDI case, while SZA raised it only
slightly (AUCR 1.16); the effect was smaller in expressers. These
numbers depend on the physicochemical and disposition inputs of Table 3,
the tissue-partition method, the organ volumes/flows and enzyme
abundances of the Simcyp physiology library, none of which is an nlmixr2
ODE model, so they are recorded here rather than reproduced.

## Assumptions and deviations

- **No residual-error model.** The source fits the RI and TDI assays by
  linear and nonlinear regression and reports only the resulting
  constants; there is no residual-error model and no assay CV. Per the
  library convention the additive residual SD is fixed at zero, so each
  model returns the deterministic published curve.
- **No inter-individual variability.** The reported `+/-` terms are
  standard errors of the regression, not between-donor variability, so
  they are recorded in `population` and not encoded as an omega. Pooled
  and genotyped microsomes are single systems, not populations.
- **Reversible Ki values for SZB.** The Discussion swaps the CYP3A4 and
  CYP3A5 reversible Ki values relative to Results 2.2 and the Figure
  2C/2D annotations; the figure/Results order is used here (CYP3A4 Ki =
  2.18 uM, CYP3A5 Ki = 2.03 uM).
- **Constants read from figure annotations.** The fitted constants are
  printed in the annotations of Figures 2-5 and repeated in the Results
  and Discussion text. They are the authors’ own estimates.
- **kinact time unit.** `kinact` is carried in per-minute units,
  matching the per-minute preincubation-time axis of the assays; the
  model time unit is minutes.
- **No enzyme turnover.** The CYP degradation rate constant `kdeg`,
  which governs the in-vivo magnitude of TDI, is a Simcyp default and is
  not printed. No turnover term is carried, so the models describe the
  inactivation phase of the in-vitro assay only; resynthesis is
  negligible over the 30-minute incubation.
- **Mechanisms carried per isoform as fitted.** SZA has no reversible
  inhibition (Figure 2A) and no CYP3A5 inactivation (Figure 4C), so that
  model carries only the CYP3A4 (and pooled) inactivation terms and
  `enzyme_3a5` holds at baseline. SZB inhibits both isoforms by both
  mechanisms. Each file also carries the pooled-HLM CYP3A screen as
  `ini()` constants exposed through derived outputs.
- **In-vivo use requires an unbound driving concentration.** Applying
  either model to a plasma profile requires supplying the unbound
  inhibitor concentration at the site of enzyme interaction as
  `CP_SZA_UM` / `CP_SZB_UM`; the in-vitro incubation concentration is
  not a plasma concentration. \`\`\`
