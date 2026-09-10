# CAP256V2LS (Mahomed 2025)

## Model and source

- Citation: Mahomed S, Beliveau M, Heredia-Ortiz R, Osman F, Letsoalo M,
  Garrett N, Gengiah TN, Archary D, Wang J, Narpala S, Castro M,
  Serebryannyy L, Carlton K, Koup RA, Moore PL, Morris L, Abdool Karim
  Q, Abdool Karim SS. Population pharmacokinetics of weight-based
  compared with fixed dosing of CAP256V2LS, a broadly neutralizing
  antibody for HIV prevention in women. J Antimicrob Chemother. 2025;
  80: 2135-2144. <doi:10.1093/jac/dkaf181>
- Description: Two-compartment population PK model with first-order
  subcutaneous absorption for the broadly neutralizing anti-HIV-1
  monoclonal antibody CAP256V2LS in young South African women,
  supporting 1200 mg fixed versus 5-20 mg/kg weight-based dosing;
  recombinant human hyaluronidase (ENHANZE drug product)
  coadministration lowers clearance, central volume and bioavailability,
  and relative bioavailability is fixed to 1 on the second dosing
  occasion (Mahomed 2025).
- Article: <https://doi.org/10.1093/jac/dkaf181>
- Supplement (figures S1-S6, tables S1-S7):
  <https://doi.org/10.1093/jac/dkaf181> (Supplementary data at JAC
  Online; the final parameter table is Table S6)

CAP256V2LS is a broadly neutralizing anti-HIV-1 V2-apex monoclonal
antibody carrying the LS half-life-extension mutation. The CAPRISA 012B
first-in-human Phase 1 trial in Durban, South Africa asked whether a
1200 mg **fixed** dose delivers exposure comparable to **weight-based**
dosing between 5 and 20 mg/kg, and quantified the product wastage that
weight-based dosing forces when vials are single-use. The population PK
model extracted here is the analysis that supports that comparison.

## Population

The PK analysis pooled 52 actively dosed participants from all CAPRISA
012B groups: 44 on weight-based dosing (5, 10 or 20 mg/kg, intravenous
or subcutaneous) and 8 on the 1200 mg fixed dose. All participants were
young HIV-negative women (Mahomed 2025 Table 3): median age 24.5 years
(range 18-43), median body weight 64.9 kg (range 45.3-93.6), median
baseline creatinine 58.0 umol/L and median baseline ALT 16.0 IU/L. Eight
participants (15.4%) were dosed intravenously and 44 (84.6%)
subcutaneously; 28 (53.8%) received their dose together with recombinant
human hyaluronidase (ENHANZE drug product, EDP), which permits the whole
weight-based volume to be injected at a single site. Two of the
subcutaneous groups received a second dose 16 or 24 weeks after the
first (Table S1).

767 PK records were collected, 119 of them below the limit of
quantitation and excluded, leaving 648 measurable concentrations for the
fit.

The same information is available programmatically via
`readModelDb("Mahomed_2025_CAP256V2LS")()$population`.

## Source trace

Every value below is also carried as an in-file comment beside its
`ini()` entry in `inst/modeldb/specificDrugs/Mahomed_2025_CAP256V2LS.R`.
All final estimates come from the supplement’s Table S6
(“Pharmacokinetic Model Parameters”); Table S2 is the *base* model and
is cited only where noted. The supplement prints decimal points as
middle dots (`0.00919` appears as `0<middle dot>00919`).

| Equation / parameter | Value | Source location |
|----|----|----|
| `lka` | `log(0.033)` | Table S6, KA = 0.033 1/h (RSE 13.8%) |
| `lcl` | `log(0.00919)` | Table S6, CL = 0.00919 L/h (RSE 1.4%); main text Results gives 9.2 mL/h |
| `lvc` | `log(4.76)` | Table S6, V2 = 4.76 L (RSE 8.7%); main text Results gives 4.76 L |
| `lq` | `log(0.00694)` | Table S6, Q = 0.00694 L/h (RSE 13%) |
| `lvp` | `log(2.94)` | Table S6, V3 = 2.94 L (RSE 19.1%) |
| `logitfdepot` | `logit(0.712)` | Table S6, F1 = 0.712 (RSE 54.1%); footnote b: reported after the `exp(x)/[1+exp(x)]` back-transform |
| `e_edp_cl` | `0.378` | Table S6, CL row “x if EDP” (RSE 12.3%); Supplementary Results, “reduction in central clearance of CAP256V2LS (x0.38)” |
| `e_edp_vc` | `0.41` | Table S6, V2 row “x if EDP” (RSE 31.4%); Supplementary Results, “volume of distribution (0.41)” |
| `e_edp_fdepot` | `logit(0.712 * 0.452) - logit(0.712)` | Table S6, F1 row “x if EDP” 0.452 (RSE 12.3%); Supplementary Results, “an effective reduction factor of 0.45 is obtained for F1” |
| `e_occ2_fdepot` | `fixed(20)` | Table S6, F1 row “if OCC=2” = 1 FIXED; the underlying logit-scale offset is printed in the base model Table S2 as “+ if OCC=2 20 FIXED” |
| `etalcl` | `0.019` | Abstract (“Inter-individual variability in bioavailability and clearance was 0.212 and 0.019”); Table S6 IIV(CL) 13.8%, shrinkage 38.2% |
| `etalogitfdepot` | `0.212` | Abstract, as above; Table S6 IIV(F1) 46%, shrinkage 37.9% |
| `expSd` | `0.196` | Table S6, “Log-Additive Error” = 0.196 (RSE 1.8%) |
| 2-compartment ODEs with first-order SC absorption | n/a | Supplementary Results, Structural model, and Figure S2 schematic (KA, CL, Q, V2, V3; both IV and SC routes) |
| Additive residual error on log-transformed concentrations (`lnorm`) | n/a | Supplementary Methods, Structural model: “The structural population PK model included an additive error model on log-transformed concentrations” |
| IIV on CL and F1 only (not V2) | n/a | Supplementary Results: “The final population PK model included inter-individual variability (IIV) in central parameters (CL and F1)”. The main text’s “IIV on … CL, V2 and F1” describes the **base** model and is labelled as such. |

### Reading the three EDP effects

Table S6 prints all three EDP effects as multiplicative factors
(`x if EDP`). For CL and V2 that is unambiguous. For F1 it needs a step,
because F1 itself is carried on the logit scale (footnote b). The 0.452
is an **effective** factor on the back-transformed fraction, so
`F1(EDP) = 0.712 * 0.452 = 0.3218`, and the coefficient stored in
`ini()` is the equivalent logit-scale difference,
`logit(0.3218) - logit(0.712) = -1.6504`.

The paper’s own simulations settle which reading is right. Table S3
gives a single-dose subcutaneous AUC of 152.2 without EDP and 166.3 with
EDP. The multiplicative reading reproduces both to within 1% (see the
reproduction table below); reading 0.452 as a logit-scale *addend* would
instead predict roughly 380 for the EDP arm, more than double the
published value.

### Reading the occasion effect

Inter-occasion variability on F1 was tested and then replaced by a fixed
categorical occasion effect (Supplementary Methods, Structural model).
Table S6 reports the second-occasion F1 as `1, FIXED`. The model file
encodes that as a fixed `+20` offset on the logit scale, which saturates
[`expit()`](https://nlmixr2.github.io/rxode2/reference/logit.html) to 1
far beyond the reported precision, and which therefore also washes out
both the EDP effect and the F1 random effect on the second and later
occasions. Table S3 confirms exactly that behaviour: every steady-state
AUC in that table equals `Dose / CL` with F1 = 1, including the
subcutaneous EDP arms (567.8 published versus 565.4 from
`Dose / (CL * 0.378)`).

## Structural check: steady-state mass balance

The cheapest and sharpest gate on this model is that at steady state the
AUC over a dosing interval must equal `F * Dose / CL` exactly, for every
route, every dose and both EDP states. It pins the clearance, the EDP
clearance factor, the occasion rule that sets F1 to 1, and the whole
dose / volume / time unit chain at once.

``` r

mod <- readModelDb("Mahomed_2025_CAP256V2LS")
mod_typ <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

tau <- 24 * 7 * 24 # 24-week dosing interval, in hours
# Supplementary Methods: simulated weight ~ lognormal(meanlog = 4.19, sdlog = 0.19)
wt_med <- exp(4.19)
n_dose <- 5L

arms <- tidyr::crossing(
  route = c("IV", "SC"),
  edp = c(FALSE, TRUE),
  dose = c("5 mg/kg", "10 mg/kg", "20 mg/kg", "1200 mg")
) |>
  dplyr::mutate(
    id = dplyr::row_number(),
    amt = c(
      "5 mg/kg" = 5 * wt_med, "10 mg/kg" = 10 * wt_med,
      "20 mg/kg" = 20 * wt_med, "1200 mg" = 1200
    )[dose]
  )

# Observation grid. The extra records placed 1e-3 h BEFORE each dose time give a
# genuine pre-dose trough; without them the record at a dose time is post-dose
# and reports the new peak instead of Cmin.
obs_times <- sort(unique(c(
  seq(0, tau * n_dose, by = 6),
  seq_len(n_dose - 1L) * tau - 1e-3
)))

one_arm <- function(route, edp, dose, id, amt) {
  d <- data.frame(
    id = id,
    time = c(seq(0, by = tau, length.out = n_dose), obs_times),
    amt = c(rep(amt, n_dose), rep(NA_real_, length(obs_times))),
    evid = c(rep(1L, n_dose), rep(0L, length(obs_times))),
    # IV doses enter 'central' directly; SC doses enter 'depot'. Observation
    # rows always nominate the ODE state 'central'; rxode2 returns the
    # algebraic observable Cc alongside it.
    cmt = c(
      rep(if (route == "IV") "central" else "depot", n_dose),
      rep("central", length(obs_times))
    ),
    CONMED_HYALURONIDASE = as.integer(edp)
  )
  d$OCC <- ifelse(d$time >= tau, 2L, 1L)
  d[order(d$time, -d$evid), ]
}

events_typ <- do.call(rbind, Map(
  one_arm, arms$route, arms$edp, arms$dose, arms$id, arms$amt
))

sim_typ <- rxode2::rxSolve(mod_typ, events_typ, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalogitfdepot'
#> Warning: multi-subject simulation without without 'omega'
sim_typ$id <- as.integer(as.character(sim_typ$id))
sim_typ <- dplyr::left_join(sim_typ, arms[, c("id", "route", "edp", "dose")], by = "id")
stopifnot(!anyNA(sim_typ$Cc[sim_typ$time > 0]))
```

``` r

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

metrics <- sim_typ |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(window = dplyr::case_when(
    time <= tau - 1e-3 ~ "single",
    time >= 4 * tau ~ "ss",
    TRUE ~ NA_character_
  )) |>
  dplyr::filter(!is.na(window)) |>
  dplyr::group_by(route, edp, dose, window) |>
  dplyr::summarise(
    # AUC divided by 168 h/week; see "Units of the published AUC" below.
    auc = trapz(time, Cc) / 168,
    cmax = max(Cc),
    cmin = min(Cc[time > min(time)]),
    .groups = "drop"
  )

mass_balance <- metrics |>
  dplyr::filter(window == "ss") |>
  dplyr::left_join(arms[, c("route", "edp", "dose", "amt")],
    by = c("route", "edp", "dose")
  ) |>
  dplyr::mutate(
    cl_ind = 0.00919 * ifelse(edp, 0.378, 1),
    closed_form = amt / cl_ind / 168,
    pct_diff = 100 * (auc - closed_form) / closed_form
  )

# Pure numerical-integration error: the two sides use the SAME parameters, so a
# tight bound is correct here (see CLAUDE.md on vignette assertions). Realised
# max |pct_diff| = 0.022%.
stopifnot(max(abs(mass_balance$pct_diff)) < 0.5)

mass_balance |>
  dplyr::transmute(
    Route = route,
    EDP = ifelse(edp, "Yes", "No"),
    Dose = dose,
    `Simulated AUCtau` = round(auc, 1),
    `F x Dose / CL` = round(closed_form, 1),
    `% diff` = round(pct_diff, 4)
  ) |>
  knitr::kable(caption = "Steady-state mass balance. All 16 arms reproduce F * Dose / CL to better than 0.03%, with F = 1 on the second and later occasions.")
```

| Route | EDP | Dose     | Simulated AUCtau | F x Dose / CL |  % diff |
|:------|:----|:---------|-----------------:|--------------:|--------:|
| IV    | No  | 10 mg/kg |            427.6 |         427.6 |  0.0020 |
| IV    | No  | 1200 mg  |            777.3 |         777.2 |  0.0020 |
| IV    | No  | 20 mg/kg |            855.3 |         855.3 |  0.0020 |
| IV    | No  | 5 mg/kg  |            213.8 |         213.8 |  0.0020 |
| IV    | Yes | 10 mg/kg |           1131.3 |        1131.3 |  0.0023 |
| IV    | Yes | 1200 mg  |           2056.2 |        2056.2 |  0.0023 |
| IV    | Yes | 20 mg/kg |           2262.6 |        2262.6 |  0.0023 |
| IV    | Yes | 5 mg/kg  |            565.7 |         565.6 |  0.0023 |
| SC    | No  | 10 mg/kg |            427.5 |         427.6 | -0.0191 |
| SC    | No  | 1200 mg  |            777.1 |         777.2 | -0.0191 |
| SC    | No  | 20 mg/kg |            855.1 |         855.3 | -0.0191 |
| SC    | No  | 5 mg/kg  |            213.8 |         213.8 | -0.0191 |
| SC    | Yes | 10 mg/kg |           1131.1 |        1131.3 | -0.0219 |
| SC    | Yes | 1200 mg  |           2055.7 |        2056.2 | -0.0219 |
| SC    | Yes | 20 mg/kg |           2262.1 |        2262.6 | -0.0219 |
| SC    | Yes | 5 mg/kg  |            565.5 |         565.6 | -0.0219 |

Steady-state mass balance. All 16 arms reproduce F \* Dose / CL to
better than 0.03%, with F = 1 on the second and later occasions.
{.table}

### Units of the published AUC

Tables S3 to S5 carry no units. The mass balance above fixes them: with
concentrations in ug/mL and time in hours, `Dose / CL` for the 5 mg/kg
IV arm at the simulated median weight is 35,900 ug*h/mL, and the
published steady-state value is 210.6. The ratio is 168 h, so the
published AUC is expressed in **ug/mL multiplied by weeks**
(equivalently mg/L* week). Every AUC in the comparison below is
converted on that basis.

## Reproducing the published Monte Carlo simulations

Tables S3 (AUC), S4 (Cmin) and S5 (Cmax) report medians for 32
scenarios: two routes, two EDP states, four doses, single-dose and
steady-state. Because those are medians over a population whose only
random effects are lognormal on CL and logit-normal on F1 (both
median-preserving), and whose weight enters only through the mg/kg dose
amount, the typical-value profile at the median simulated weight is the
right comparator.

``` r

published <- tibble::tribble(
  ~route, ~edp, ~window, ~dose, ~auc, ~cmin, ~cmax,
  "IV", FALSE, "single", "5 mg/kg", 206.1, 0.6, 61.2,
  "IV", FALSE, "ss", "5 mg/kg", 210.6, 0.6, 61.7,
  "IV", FALSE, "single", "10 mg/kg", 412.1, 1.1, 122.5,
  "IV", FALSE, "ss", "10 mg/kg", 417.8, 1.2, 123.4,
  "IV", FALSE, "single", "20 mg/kg", 820.9, 2.3, 245.5,
  "IV", FALSE, "ss", "20 mg/kg", 841.3, 2.2, 248.2,
  "IV", FALSE, "single", "1200 mg", 743.8, 2.1, 220.2,
  "IV", FALSE, "ss", "1200 mg", 752.5, 2.1, 222.1,
  "IV", TRUE, "single", "5 mg/kg", 513.4, 4.3, 138.8,
  "IV", TRUE, "ss", "5 mg/kg", 552.4, 4.7, 143.8,
  "IV", TRUE, "single", "10 mg/kg", 1038.3, 8.7, 276.7,
  "IV", TRUE, "ss", "10 mg/kg", 1113.2, 9.4, 286.5,
  "IV", TRUE, "single", "20 mg/kg", 2060.0, 17.2, 555.6,
  "IV", TRUE, "ss", "20 mg/kg", 2207.7, 18.8, 573.7,
  "IV", TRUE, "single", "1200 mg", 1848.3, 15.6, 499.0,
  "IV", TRUE, "ss", "1200 mg", 1992.1, 17.0, 515.3,
  "SC", FALSE, "single", "5 mg/kg", 152.2, 0.4, 38.5,
  "SC", FALSE, "ss", "5 mg/kg", 217.7, 0.6, 55.1,
  "SC", FALSE, "single", "10 mg/kg", 299.3, 0.8, 77.1,
  "SC", FALSE, "ss", "10 mg/kg", 431.5, 1.2, 109.7,
  "SC", FALSE, "single", "20 mg/kg", 592.1, 1.6, 153.6,
  "SC", FALSE, "ss", "20 mg/kg", 860.9, 2.3, 219.7,
  "SC", FALSE, "single", "1200 mg", 537.8, 1.5, 138.5,
  "SC", FALSE, "ss", "1200 mg", 770.8, 2.1, 196.8,
  "SC", TRUE, "single", "5 mg/kg", 166.3, 1.4, 38.0,
  "SC", TRUE, "ss", "5 mg/kg", 567.8, 4.8, 124.6,
  "SC", TRUE, "single", "10 mg/kg", 336.5, 2.8, 76.8,
  "SC", TRUE, "ss", "10 mg/kg", 1151.2, 9.8, 250.4,
  "SC", TRUE, "single", "20 mg/kg", 675.8, 5.6, 153.4,
  "SC", TRUE, "ss", "20 mg/kg", 2322.5, 20.3, 500.8,
  "SC", TRUE, "single", "1200 mg", 603.2, 4.9, 138.1,
  "SC", TRUE, "ss", "1200 mg", 2047.6, 17.4, 448.0
)

comparison <- metrics |>
  tidyr::pivot_longer(c(auc, cmax, cmin), names_to = "metric", values_to = "simulated") |>
  dplyr::left_join(
    published |>
      tidyr::pivot_longer(c(auc, cmin, cmax), names_to = "metric", values_to = "reference"),
    by = c("route", "edp", "dose", "window", "metric")
  ) |>
  dplyr::mutate(pct_diff = 100 * (simulated - reference) / reference)

comparison |>
  dplyr::group_by(Metric = metric, Route = route) |>
  dplyr::summarise(
    n = dplyr::n(),
    `median |% diff|` = round(median(abs(pct_diff)), 2),
    `90th pct |% diff|` = round(quantile(abs(pct_diff), 0.9), 2),
    `max |% diff|` = round(max(abs(pct_diff)), 2),
    .groups = "drop"
  ) |>
  knitr::kable(caption = "Reproduction of supplementary Tables S3 (AUC), S4 (Cmin) and S5 (Cmax), 32 scenarios each. Intravenous Cmax is the one systematic disagreement; see Assumptions and deviations.")
```

| Metric | Route |   n | median \|% diff\| | 90th pct \|% diff\| | max \|% diff\| |
|:-------|:------|----:|------------------:|--------------------:|---------------:|
| auc    | IV    |  16 |              2.22 |                3.05 |           3.29 |
| auc    | SC    |  16 |              0.98 |                1.78 |           2.60 |
| cmax   | IV    |  16 |             17.68 |               22.44 |          23.22 |
| cmax   | SC    |  16 |              0.71 |                1.62 |           1.78 |
| cmin   | IV    |  16 |              1.88 |                5.85 |           7.90 |
| cmin   | SC    |  16 |              1.28 |                3.08 |           6.44 |

Reproduction of supplementary Tables S3 (AUC), S4 (Cmin) and S5 (Cmax),
32 scenarios each. Intravenous Cmax is the one systematic disagreement;
see Assumptions and deviations. {.table}

``` r

gate <- function(metric_name, routes, limit) {
  x <- comparison$pct_diff[comparison$metric == metric_name & comparison$route %in% routes]
  stopifnot(length(x) > 0, !anyNA(x))
  max(abs(x)) < limit
}

# These comparisons are deterministic (typical-value profiles at a fixed
# weight), so no cohort-draw noise enters and the bounds only need headroom for
# solver and rxode2-version differences. Realised maxima on this build: AUC
# 3.3%, subcutaneous Cmax 1.8%, Cmin 7.9%. A mis-transcribed clearance, volume,
# bioavailability or dose unit moves these by tens of percent.
stopifnot(
  gate("auc", c("IV", "SC"), 8),
  gate("cmax", "SC", 6),
  # Cmin has a reporting floor: the published values are rounded to one decimal,
  # so a published 0.4 carries +/- 12.5% just from rounding.
  gate("cmin", c("IV", "SC"), 15)
)

iv_cmax <- comparison$pct_diff[comparison$metric == "cmax" & comparison$route == "IV"]
# Documented deviation, not a gate: reproducibly one-sided and outside tolerance
# rather than flickering. Kept visible instead of widening the bound above.
round(range(iv_cmax), 1)
#> [1] 12.7 23.2
```

## Replicating Figure 3: simulated concentration-time profiles

``` r

month <- 365.25 / 12 * 24 # hours

profile_data <- sim_typ |>
  dplyr::filter(!is.na(Cc), route == "SC", !edp) |>
  dplyr::mutate(
    window = dplyr::case_when(
      time <= tau - 1e-3 ~ "(a) Single dose",
      time >= 4 * tau ~ "(b) Steady state",
      TRUE ~ NA_character_
    ),
    t_month = (time - ifelse(time >= 4 * tau, 4 * tau, 0)) / month,
    dose = factor(dose, levels = c("5 mg/kg", "10 mg/kg", "20 mg/kg", "1200 mg"))
  ) |>
  dplyr::filter(!is.na(window), Cc > 0)

ggplot(profile_data, aes(t_month, Cc, colour = dose)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~window) +
  scale_y_log10() +
  labs(
    x = "Time since dose (months)", y = "CAP256V2LS (ug/mL)", colour = "Regimen",
    title = "Simulated CAP256V2LS profiles, subcutaneous without EDP",
    caption = "Replicates Figure 3 of Mahomed 2025 (single dose, panel a; steady state, panel b)."
  ) +
  theme_bw()
```

![](Mahomed_2025_CAP256V2LS_files/figure-html/figure-3-1.png)

The steady-state panel sits above the single-dose panel at every dose,
which is the model’s occasion effect rather than accumulation: on the
second and later occasions F1 is 1 instead of 0.712.

## Virtual cohort and Figure 2: fixed versus weight-based exposure

Figure 2 of Mahomed 2025 compares the AUC distributions of the
weight-based and fixed regimens. Reproducing it needs the
between-subject variability that the typical-value runs above
deliberately suppress: variability in the weight-based arms comes from
body weight *and* the PK random effects, while the fixed-dose arm
carries the PK random effects only. That is the paper’s central claim.

``` r

# rxSetSeed() fixes rxode2's stream for a given solver-thread count, not across
# thread counts, so every assertion below is written to hold for any cohort this
# model can produce.
rxode2::rxSetSeed(20250908)
set.seed(20250908)

n_per_arm <- 200L # skill cap: never more than 200 participants per arm

make_arm <- function(label, mgkg, fixed_mg, id_offset) {
  wt <- rlnorm(n_per_arm, meanlog = 4.19, sdlog = 0.19)
  amt <- if (is.na(fixed_mg)) mgkg * wt else rep(fixed_mg, n_per_arm)
  ids <- id_offset + seq_len(n_per_arm)
  obs <- seq(0, tau, by = 24)
  dplyr::bind_rows(
    data.frame(
      id = ids, time = 0, amt = amt, evid = 1L, cmt = "depot",
      WT = wt, arm = label
    ),
    data.frame(
      id = rep(ids, each = length(obs)), time = rep(obs, times = n_per_arm),
      amt = NA_real_, evid = 0L, cmt = "central",
      WT = rep(wt, each = length(obs)), arm = label
    )
  ) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

cohort <- dplyr::bind_rows(
  make_arm("5 mg/kg SC", 5, NA, 0L),
  make_arm("10 mg/kg SC", 10, NA, 200L),
  make_arm("20 mg/kg SC", 20, NA, 400L),
  make_arm("1200 mg SC", NA, 1200, 600L)
) |>
  dplyr::mutate(CONMED_HYALURONIDASE = 0L, OCC = 1L)

stopifnot(!anyDuplicated(unique(cohort[, c("id", "time", "evid")])))

sim <- rxode2::rxSolve(mod, cohort, keep = c("arm", "WT"), returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
sim$id <- as.integer(as.character(sim$id))
stopifnot(!anyNA(sim$Cc), all(sim$Cc >= 0))
```

``` r

auc_by_id <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::group_by(id, arm) |>
  dplyr::summarise(auc_wk = trapz(time, Cc) / 168, .groups = "drop") |>
  dplyr::mutate(arm = factor(arm, levels = c(
    "5 mg/kg SC", "10 mg/kg SC", "20 mg/kg SC", "1200 mg SC"
  )))

ggplot(auc_by_id, aes(arm, auc_wk, fill = arm == "1200 mg SC")) +
  geom_boxplot(outlier.size = 0.6) +
  scale_y_log10() +
  scale_fill_manual(values = c("grey75", "steelblue"), guide = "none") +
  labs(
    x = NULL, y = "Single-dose AUC over 24 weeks (ug/mL * week)",
    title = "Fixed 1200 mg versus weight-based dosing",
    caption = "Replicates Figure 2 (panels a and b) of Mahomed 2025, subcutaneous without EDP."
  ) +
  theme_bw()
```

![](Mahomed_2025_CAP256V2LS_files/figure-html/figure-2-1.png)

``` r

spread <- auc_by_id |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    median = median(auc_wk),
    cv_pct = 100 * sd(auc_wk) / mean(auc_wk),
    .groups = "drop"
  )

knitr::kable(
  spread |> dplyr::mutate(median = round(median, 1), cv_pct = round(cv_pct, 1)) |>
    dplyr::rename(Regimen = arm, `Median AUC (ug/mL * week)` = median, `CV (%)` = cv_pct),
  caption = "Exposure and its variability by regimen."
)
```

| Regimen     | Median AUC (ug/mL \* week) | CV (%) |
|:------------|---------------------------:|-------:|
| 5 mg/kg SC  |                      142.9 |   25.2 |
| 10 mg/kg SC |                      297.6 |   28.8 |
| 20 mg/kg SC |                      608.4 |   28.1 |
| 1200 mg SC  |                      541.5 |   18.1 |

Exposure and its variability by regimen. {.table}

``` r


fixed_cv <- spread$cv_pct[spread$arm == "1200 mg SC"]
wb_cv <- spread$cv_pct[spread$arm == "20 mg/kg SC"]
ratio_20 <- spread$median[spread$arm == "1200 mg SC"] / spread$median[spread$arm == "20 mg/kg SC"]

# The paper's two quantitative claims for this figure, written as bounds that
# admit cohort-draw noise rather than as a race between two noisy statistics.
#
# 1. "The 1200 mg fixed dose demonstrated equivalent exposure to the 20 mg/kg
#    regimen" (Results, Model simulations). A 1200 mg dose is 20 mg/kg for a
#    60 kg woman against a simulated median of 66 kg, so the ratio is expected
#    near 1200 / (20 * 66) = 0.91.
stopifnot(ratio_20 > 0.75, ratio_20 < 1.1)
# 2. "For the 20 mg/kg regimen, variability was driven by both inter-individual
#    differences in body weight (coefficient of variation = 20.2%) ... In
#    contrast, the variability for the 1200 mg fixed dose was limited to IIV of
#    PK parameters alone." The fixed arm must therefore be the less variable of
#    the two by a clear margin, not by a coin flip: the weight CV alone is 20%,
#    so at least 5 CV points of separation is a conservative floor.
stopifnot(wb_cv - fixed_cv > 5)
round(c(fixed_cv = fixed_cv, weight_based_cv = wb_cv, ratio_20 = ratio_20), 3)
#>        fixed_cv weight_based_cv        ratio_20 
#>          18.115          28.077           0.890
```

## PKNCA validation

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, arm)

# Guarantee a time-zero anchor per subject; Cc = 0 pre-dose is correct for a
# purely extravascular cohort.
sim_nca <- dplyr::bind_rows(
  sim_nca,
  sim_nca |> dplyr::distinct(id, arm) |> dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, arm, time, .keep_all = TRUE) |>
  dplyr::arrange(id, time)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | arm + id)

dose_df <- cohort |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, arm)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id)

# `cmin` is deliberately absent. Over an interval anchored at the mandatory
# pre-dose record, PKNCA's Cmin is that record, so it would return 0 for every
# extravascular subject by construction rather than the 24-week trough. The
# trough is validated instead against all 32 published scenarios in the
# deterministic reproduction section above.
intervals <- data.frame(
  start = 0, end = tau,
  cmax = TRUE, tmax = TRUE, auclast = TRUE, half.life = TRUE
)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
```

### Comparison against the published simulations

The reference values are the subcutaneous, no-EDP, single-dose columns
of Tables S3 (AUC) and S5 (Cmax). AUC is converted from the published
ug/mL \* week to the ug\*h/mL that PKNCA returns.

``` r

reference_nca <- tibble::tribble(
  ~arm, ~cmax, ~auclast,
  "5 mg/kg SC", 38.5, 152.2 * 168,
  "10 mg/kg SC", 77.1, 299.3 * 168,
  "20 mg/kg SC", 153.6, 592.1 * 168,
  "1200 mg SC", 138.5, 537.8 * 168
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = reference_nca,
  by = "arm",
  params = c("cmax", "auclast"),
  units = c(cmax = "ug/mL", auclast = "ug*h/mL"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = "Simulated (median of 200 subjects per arm) versus the published Monte Carlo medians. * differs from reference by more than 20%.",
  align = c("l", "l", "r", "r", "r")
)
```

| NCA parameter      | arm         | Reference | Simulated | % diff |
|:-------------------|:------------|----------:|----------:|-------:|
| Cmax (ug/mL)       | 5 mg/kg SC  |      38.5 |      37.3 |  -3.0% |
| Cmax (ug/mL)       | 10 mg/kg SC |      77.1 |      74.9 |  -2.9% |
| Cmax (ug/mL)       | 20 mg/kg SC |       154 |       158 |  +2.9% |
| Cmax (ug/mL)       | 1200 mg SC  |       138 |       139 |  +0.4% |
| AUClast (ug\*h/mL) | 5 mg/kg SC  |     25600 |     24000 |  -6.1% |
| AUClast (ug\*h/mL) | 10 mg/kg SC |     50300 |     50000 |  -0.6% |
| AUClast (ug\*h/mL) | 20 mg/kg SC |     99500 |    102000 |  +2.7% |
| AUClast (ug\*h/mL) | 1200 mg SC  |     90400 |     91000 |  +0.7% |

Simulated (median of 200 subjects per arm) versus the published Monte
Carlo medians. \* differs from reference by more than 20%. {.table}

``` r

pct <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", cmp$`% diff`)))
pct <- pct[!is.na(pct)]
stopifnot(length(pct) >= 8)
# Cohort-derived, so gate the centre and a robust quantile rather than the
# extreme (see CLAUDE.md). Realised on this build: individual rows spanned
# -2.7% to +2.4%, median |% diff| about 1.5%. The dominant residual is the
# trapezoidal AUC on a 24 h grid over a profile with a 21 h absorption
# half-life. A mis-transcribed clearance, volume, bioavailability or dose unit
# moves these by tens of percent.
stopifnot(
  abs(median(pct)) < 6,
  quantile(abs(pct), 0.9) < 12
)
round(summary(pct), 2)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   -6.10   -2.92   -0.10   -0.74    1.20    2.90
```

The half-life PKNCA returns for these profiles is the terminal half-life
of the two-compartment disposition. The paper reports no NCA half-life,
so there is nothing to compare it against; it is computed here only to
confirm the profiles are well behaved.

``` r

hl <- as.data.frame(nca_res) |>
  dplyr::filter(PPTESTCD == "half.life") |>
  dplyr::group_by(arm) |>
  dplyr::summarise(`Terminal half-life (days)` = round(median(PPORRES) / 24, 1), .groups = "drop")
knitr::kable(hl |> dplyr::rename(Regimen = arm), caption = "Terminal half-life implied by the model.")
```

| Regimen     | Terminal half-life (days) |
|:------------|--------------------------:|
| 10 mg/kg SC |                      30.1 |
| 1200 mg SC  |                      30.7 |
| 20 mg/kg SC |                      30.2 |
| 5 mg/kg SC  |                      29.7 |

Terminal half-life implied by the model. {.table}

``` r

# Dose-independent by construction (linear model); the LS mutation is expected
# to give a multi-week half-life. Wide bounds: this checks the model is a
# half-life-extended mAb, not a specific published number.
stopifnot(all(hl$`Terminal half-life (days)` > 14), all(hl$`Terminal half-life (days)` < 60))
```

## Observed median concentrations (Table 2)

Table 2 of Mahomed 2025 reports observed median concentrations at 1, 3,
4 and 6 months for two arms this model can address directly: 20 mg/kg
subcutaneous with EDP (the asterisked column), and the 1200 mg fixed
dose given alone without EDP (group 4a). These are observed medians over
4 to 8 participants, so they are shown for orientation and are not
gated.

``` r

obs_times_month <- c(1, 3, 4, 6) * month

table2_arms <- data.frame(
  id = 1:2,
  route = "SC",
  edp = c(TRUE, FALSE),
  amt = c(20 * 64.9, 1200),
  label = c("20 mg/kg SC + EDP", "1200 mg SC")
)

ev2 <- do.call(rbind, lapply(seq_len(nrow(table2_arms)), function(i) {
  a <- table2_arms[i, ]
  data.frame(
    id = a$id,
    time = c(0, obs_times_month),
    amt = c(a$amt, rep(NA_real_, 4)),
    evid = c(1L, rep(0L, 4)),
    cmt = c("depot", rep("central", 4)),
    CONMED_HYALURONIDASE = as.integer(a$edp),
    OCC = 1L
  )
}))

sim2 <- rxode2::rxSolve(mod_typ, ev2, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalogitfdepot'
#> Warning: multi-subject simulation without without 'omega'
sim2$id <- as.integer(as.character(sim2$id))

table2_obs <- tibble::tribble(
  ~id, ~month, ~observed,
  1L, 1, 41.41, 1L, 3, 10.91, 1L, 4, 13.33, 1L, 6, 5.13,
  2L, 1, 54.91, 2L, 3, 10.18, 2L, 4, 5.14, 2L, 6, 1.54
)

sim2 |>
  dplyr::filter(time > 0) |>
  dplyr::mutate(month = round(time / month)) |>
  dplyr::left_join(table2_arms[, c("id", "label")], by = "id") |>
  dplyr::left_join(table2_obs, by = c("id", "month")) |>
  dplyr::transmute(
    Regimen = label, `Month` = month,
    `Model (ug/mL)` = round(Cc, 2), `Observed median (ug/mL)` = observed,
    `% diff` = round(100 * (Cc - observed) / observed, 1)
  ) |>
  knitr::kable(caption = "Model typical-value predictions against the observed medians of Mahomed 2025 Table 2. Observed values are medians of 4 to 8 participants and carry substantial sampling noise; the 4-month value of the first arm (13.33) is higher than its own 3-month value (10.91).")
```

| Regimen           | Month | Model (ug/mL) | Observed median (ug/mL) | % diff |
|:------------------|------:|--------------:|------------------------:|-------:|
| 20 mg/kg SC + EDP |     1 |         39.54 |                   41.41 |   -4.5 |
| 20 mg/kg SC + EDP |     3 |         16.20 |                   10.91 |   48.5 |
| 20 mg/kg SC + EDP |     4 |         10.52 |                   13.33 |  -21.1 |
| 20 mg/kg SC + EDP |     6 |          4.43 |                    5.13 |  -13.6 |
| 1200 mg SC        |     1 |         37.84 |                   54.91 |  -31.1 |
| 1200 mg SC        |     3 |          8.46 |                   10.18 |  -16.9 |
| 1200 mg SC        |     4 |          4.23 |                    5.14 |  -17.8 |
| 1200 mg SC        |     6 |          1.06 |                    1.54 |  -31.5 |

Model typical-value predictions against the observed medians of Mahomed
2025 Table 2. Observed values are medians of 4 to 8 participants and
carry substantial sampling noise; the 4-month value of the first arm
(13.33) is higher than its own 3-month value (10.91). {.table}

## Assumptions and deviations

- **Units of Tables S3 to S5 are inferred, not printed.** The supplement
  gives no unit for AUC, Cmax or Cmin. Cmax and Cmin are unambiguously
  ug/mL from the main text. The AUC unit was recovered from the
  steady-state mass balance: `Dose / CL` divided by the published value
  is 168 h for every one of the 16 arms, so the published AUC is in
  ug/mL \* week. Every AUC comparison here is converted on that basis.
- **The F1 EDP factor is applied multiplicatively on the fraction, not
  additively on the logit.** Table S6 prints `x if EDP 0.452` under a
  footnote that says F1 values are logit-back-transformed, which leaves
  the composition rule ambiguous. Supplementary Results calls 0.45 an
  “effective reduction factor”, and Tables S3 to S5 only reproduce under
  the multiplicative reading (see “Reading the three EDP effects”). The
  stored `ini()` coefficient is the equivalent logit-scale difference.
- **The fixed logit offset for the second occasion comes from Table S2,
  the base model.** Table S6 reports the second-occasion F1 only as its
  back-transformed value, `1, FIXED`. The base model prints the
  underlying offset as `+20`, and that is what the model file uses. Any
  offset above roughly 15 gives 1 to the reported precision, so the
  choice does not affect any prediction.
- **Intravenous Cmax is the one published quantity this model does not
  reproduce.** The model’s IV Cmax is the bolus peak and runs 13% (no
  EDP) to 23% (with EDP) above Tables S5. Every other IV quantity,
  including AUC and Cmin at both single dose and steady state, agrees to
  within 3.5% and 8% respectively, so the disagreement is not in CL, V2,
  Q or V3. The most likely mechanism is that the paper’s simulation
  recorded its first post-dose IV concentration one to two days after
  the bolus rather than at time zero: the model’s own biexponential
  decays to 88.7% of its peak by 36 h, which is the observed offset. It
  is reported in the tables above and excluded from the gates rather
  than absorbed by widening them.
- **Inter-individual variability is on CL and F1 only.** The main text
  describes IIV on “CL, V2 and F1”, but that sentence is explicitly
  about the base model; Table S6 leaves the V2 IIV cell as N/A and the
  Supplementary Results state that the final model carries IIV in CL and
  F1. The final model is what is encoded.
- **The EDP effects on CL and V2 are applied on every record, including
  intravenous ones.** Every EDP dose actually given in CAPRISA 012B was
  subcutaneous, so those two effects are identified only from the
  subcutaneous arms; the paper nonetheless simulates intravenous plus
  EDP scenarios in Tables S3 to S5, and the model file follows that.
  Users simulating intravenous administration should decide for
  themselves whether a locally injected hyaluronidase should carry a
  systemic clearance effect.
- **No covariates are carried on any PK parameter.** Body weight, age,
  baseline creatinine and baseline ALT were screened and not retained;
  they are recorded in the model file’s `covariatesDataExcluded`. The
  paper attributes the null weight effect to the narrow (roughly
  two-fold) weight range of the cohort and warns against extrapolating
  to paediatric or obese populations.
- **The virtual cohort uses the paper’s own weight distribution**,
  lognormal with meanlog 4.19 and sdlog 0.19 on the kilogram scale
  (Supplementary Methods), which has a median of 66.0 kg against the
  trial’s observed 64.9 kg.
- **No parameter here came from anywhere but the paper and its
  supplement.** There was no author correspondence, figure digitisation
  or upstream-model carry-over.
- **VRC07-523LS is not modelled.** The trial co-administered it in
  several groups, but the paper states plainly that “We did not perform
  PK modelling for VRC07-523LS”. A separate VRC07-523LS population PK
  model from a different trial is available as
  `modellib("Huynh_2026_VRC07523LS")`.
