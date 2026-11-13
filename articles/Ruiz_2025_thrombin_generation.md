# Thrombin generation in normal and trauma plasma (Ruiz 2025)

## Model and source

- Citation: Ruiz L, Jaramillo S, Calvo A, Torrente MA, Tassies D,
  Reverter JC, Blasi A, Troconiz I. Dynamics of thrombin generation:
  Filling the gap between the system pharmacology theory and clinical
  practice in clinical pharmacology and therapeutics. Pharmacol Res
  Perspect. 2025;13(1):e70014. <doi:10.1002/prp2.70014>
- Description: In vitro (calibrated automated thrombogram, platelet-poor
  plasma, 5 pM recombinant tissue factor). Semi-mechanistic model of
  thrombin generation in normal subjects and trauma patients: tissue
  factor activates prothrombin directly (second order, with thrombin
  positive feedback) and through a seven-compartment transit chain
  (second order); thrombin is degraded first order. All rate constants
  and their IIV are estimated separately for normal and trauma subjects.
- Article: <https://doi.org/10.1002/prp2.70014> (open access;
  PMC11687557)

Ruiz et al. reduced the large systems-biology coagulation models of
Wajima (2009) and Nayak (2015) to a semi-mechanistic model of the
thrombin-generation curve measured by the calibrated automated
thrombogram (CAT). Tissue factor (TF) activates prothrombin (FII) to
thrombin (FIIa) through two pathways:

- a direct pathway, second order in TF and FII, amplified by the
  thrombin already formed (`ksynIIa * TF * FII * (1 + FIIa)`), which
  produces the early burst; and
- a delayed pathway in which TF passes through a chain of seven transit
  compartments (`ktran`) whose last member activates FII
  (`ksynIIa2 * T7 * FII`), which produces most of the thrombin.

Thrombin is removed first order (`kdegIIa`). TF degradation (`kdegTF`)
was negligible and fixed to zero. Every rate constant, and its
inter-individual variance, was estimated separately for normal subjects
and trauma patients.

## Population

The data are the thrombin-generation curves of Menezes et al.: 20 normal
(control) subjects and 40 trauma patients, with 1978 and 3736 thrombin
observations respectively over 0-25 min. Each assay used platelet-poor
plasma triggered with 5 pM recombinant tissue factor and 4 uM
procoagulant phospholipids. Baseline factor II, V, VII, VIII, IX, X and
antithrombin were measured in nM (Ruiz 2025 Table 1); only factor II
enters the model, as the initial condition of the prothrombin state.
Normal subjects had a median FII of 1108.5 nM (725.1-1338.06) and trauma
patients 1052.8 nM (474.09-1784.83). No age, weight or sex data are
reported. The same information is available programmatically via
`readModelDb("Ruiz_2025_thrombin_generation")()$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in
`inst/modeldb/endogenous/Ruiz_2025_thrombin_generation.R`.

| Equation / parameter | Value (normal / trauma) | Source location |
|----|----|----|
| `d/dt(tissue_factor)` | `-kdeg_tf * tissue_factor` | Equation 1; Figure 2 |
| `d/dt(effect1..effect7)` | `ktr * (input - effect<n>)`, input to `effect1` = `tissue_factor` | Figure 2 (T1-T7) |
| `d/dt(prothrombin)` | `-ksyn_tf * TF * FII * (1 + FIIa) - ksyn_effect * T7 * FII` | Figure 2 |
| `d/dt(thrombin)` | `+ksyn_tf * TF * FII * (1 + FIIa) + ksyn_effect * T7 * FII - kdeg_thrombin * FIIa` | Figure 2 |
| `tf0` | 0.005 nM, fixed | Methods, ‘final concentration of 5 pM of recombinant tissue factor’ |
| `kdeg_tf` | 0, fixed | Results and Discussion, ‘fixed to 0’ |
| `prothrombin(0)` | `PROTHROMBIN_BASE` | Discussion, model ‘requires the values of FII and TF’; Table 1 |
| `lksyn_thrombin_tf` / `_trauma` | log(4.7e-3) / log(0.15) nM^-1 min^-1 | Table 2, ksynIIa |
| `lkdeg_thrombin` / `_trauma` | log(1.04) / log(0.8) min^-1 | Table 2, kdegIIa |
| `lktr` / `lktr_trauma` | log(1.4) / log(1.8) min^-1 | Table 2, ktran |
| `lksyn_thrombin_effect` / `_trauma` | log(47) / log(45) nM^-1 min^-1 | Table 2, ksynII2a |
| `etalksyn_thrombin_tf_*` | 2.40^2 / 0.51^2 (CV 240% / 51%) | Table 2 IIV, footnote CV = sqrt(omega^2) x 100 |
| `etalkdeg_thrombin_*` | 0.20^2 / 0.25^2 (CV 20% / 25%) | Table 2 IIV |
| `etalktr_*` | 0.23^2 / 0.42^2 (CV 23% / 42%) | Table 2 IIV |
| `etalksyn_thrombin_effect_*` | 0.14^2 / 0.34^2 (CV 14% / 34%) | Table 2 IIV |
| `propSd` | 0.18 (before 15 min) | Table 2, ‘Proportinal error’ |
| `addSd` | 2.4 nM (before 15 min) | Table 2, Additive error Time \< 15 min |
| `addSd_late` | 2.4 nM (from 15 min) | Table 2, Additive error Time \> 15 min |
| Residual structure | combined before 15 min, additive after | Results, ‘Mechanism-based model’ |

### Units

Time is in minutes and every state is a molar concentration (nM) in the
assay plasma. The transit states carry the TF signal, so they share its
nM units.

| Term | Units | Product |
|----|----|----|
| `kdeg_tf * tissue_factor` | min^-1 x nM | nM/min |
| `ktr * effect<n>` | min^-1 x nM | nM/min |
| `ksyn_thrombin_effect * effect7 * prothrombin` | nM^-1 min^-1 x nM x nM | nM/min |
| `ksyn_thrombin_tf * tissue_factor * prothrombin` | nM^-1 min^-1 x nM x nM | nM/min |
| `(1 + thrombin)` | 1 + nM | see note |
| `kdeg_thrombin * thrombin` | min^-1 x nM | nM/min |

The thrombin feedback factor `(1 + FIIa)` adds a dimensionless 1 to a
concentration. It is reproduced exactly as printed in Figure 2, with
thrombin entering as its value in nM; the fitted `ksynIIa` is only
meaningful in that form, so it must not be rescaled or rewritten.

## Model structure checks

These checks use typical values (`zeroRe()`) and a single subject per
condition at the Table 1 median factor II.

``` r

mod <- readModelDb("Ruiz_2025_thrombin_generation")
mod_typ <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalksyn_thrombin_tf_trauma, etalkdeg_thrombin_trauma, etalktr_trauma, etalksyn_thrombin_effect_trauma
#> as a work-around try putting the mu-referenced expression on a simple line
#> Warning: No sigma parameters in the model
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalksyn_thrombin_tf_trauma, etalkdeg_thrombin_trauma, etalktr_trauma, etalksyn_thrombin_effect_trauma
#> as a work-around try putting the mu-referenced expression on a simple line

typ_cov <- tibble(
  id = 1:2,
  condition = c("Normal", "Trauma"),
  DIS_TRAUMA = c(0, 1),
  PROTHROMBIN_BASE = c(1108.5, 1052.8)
)
typ_ev <- tidyr::crossing(typ_cov, time = seq(0, 25, by = 0.01)) |>
  mutate(evid = 0, cmt = "thrombin") |>
  arrange(id, time)

typ <- rxode2::rxSolve(mod_typ, typ_ev, keep = "condition",
                       returnType = "data.frame", atol = 1e-10, rtol = 1e-10)
#> ℹ omega/sigma items treated as zero: 'etalksyn_thrombin_tf_normal', 'etalksyn_thrombin_tf_trauma', 'etalkdeg_thrombin_normal', 'etalkdeg_thrombin_trauma', 'etalktr_normal', 'etalktr_trauma', 'etalksyn_thrombin_effect_normal', 'etalksyn_thrombin_effect_trauma'
#> Warning: multi-subject simulation without without 'omega'
```

### Transit chain against its closed form

With `kdeg_tf = 0` the tissue factor is constant, so the seventh transit
state is the input times a gamma(7, `ktr`) cumulative distribution:
`effect7(t) = tf0 * pgamma(t, 7, ktr)`. Both sides use the same
parameters, so the difference is pure integration error.

``` r

ktr_typ <- c(Normal = 1.4, Trauma = 1.8)
transit_chk <- typ |>
  mutate(closed = 0.005 * pgamma(time, shape = 7, rate = ktr_typ[condition]))
stopifnot(
  nrow(transit_chk) == 2 * 2501,
  max(abs(typ$tissue_factor - 0.005)) < 1e-12,
  max(abs(transit_chk$effect7 - transit_chk$closed)) < 1e-7
)
transit_chk |>
  group_by(condition) |>
  summarise(max_abs_err_nM = max(abs(effect7 - closed)), .groups = "drop") |>
  knitr::kable(digits = 12)
```

| condition | max_abs_err_nM |
|:----------|---------------:|
| Normal    |          2e-12 |
| Trauma    |          3e-12 |

### Mass balance

Prothrombin is converted one-for-one into thrombin, and thrombin is
removed only by `kdeg_thrombin`, so
`prothrombin(t) + thrombin(t) + kdeg_thrombin * integral(thrombin) = PROTHROMBIN_BASE`.
The integral is a trapezoid on the 0.01-min grid.

``` r

kdeg_typ <- c(Normal = 1.04, Trauma = 0.8)
mb <- typ |>
  group_by(condition) |>
  arrange(time, .by_group = TRUE) |>
  mutate(
    removed = kdeg_typ[condition] *
      cumsum(c(0, diff(time) * (head(thrombin, -1) + tail(thrombin, -1)) / 2)),
    total = prothrombin + thrombin + removed,
    rel_err = total / PROTHROMBIN_BASE - 1
  ) |>
  summarise(max_abs_rel_err = max(abs(rel_err)),
            fii_remaining_25min = last(prothrombin), .groups = "drop")
stopifnot(nrow(mb) == 2, all(mb$max_abs_rel_err < 1e-4))
knitr::kable(mb, digits = 6)
```

| condition | max_abs_rel_err | fii_remaining_25min |
|:----------|----------------:|--------------------:|
| Normal    |           0e+00 |            9.829762 |
| Trauma    |           1e-06 |            3.346256 |

### Typical thrombin-generation curves

``` r

ggplot(typ, aes(time, thrombin, colour = condition)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(Normal = "#d7301f", Trauma = "#1b9e9e")) +
  labs(x = "Time (min)", y = "Thrombin, FIIa (nM)", colour = NULL) +
  theme_bw()
```

![Typical-value thrombin-generation curves at the Table 1 median factor
II. Compare with the median lines of the Figure S2 VPCs of Ruiz
2025.](Ruiz_2025_thrombin_generation_files/figure-html/typical-plot-1.png)

Typical-value thrombin-generation curves at the Table 1 median factor
II. Compare with the median lines of the Figure S2 VPCs of Ruiz 2025.

``` r


typ_pk <- typ |>
  group_by(condition) |>
  summarise(peak_nM = max(thrombin), time_to_peak_min = time[which.max(thrombin)],
            .groups = "drop")
knitr::kable(typ_pk, digits = 2)
```

| condition | peak_nM | time_to_peak_min |
|:----------|--------:|-----------------:|
| Normal    |  125.74 |             7.25 |
| Trauma    |  217.86 |             5.79 |

The typical trauma curve peaks higher and earlier than the normal one,
as in the paper’s Figures 1, 3 and S2: trauma raises the direct
activation rate about 30-fold and speeds up the transit chain, while
thrombin is removed more slowly.

``` r

pk_n <- typ_pk$peak_nM[typ_pk$condition == "Normal"]
pk_t <- typ_pk$peak_nM[typ_pk$condition == "Trauma"]
tp_n <- typ_pk$time_to_peak_min[typ_pk$condition == "Normal"]
tp_t <- typ_pk$time_to_peak_min[typ_pk$condition == "Trauma"]
stopifnot(
  length(c(pk_n, pk_t, tp_n, tp_t)) == 4,
  pk_t > 1.4 * pk_n,
  tp_t < tp_n - 1
)
```

## Virtual cohort

The observed curves are not public. The cohort below has 200 subjects
per condition, with factor II drawn log-normally around each group’s
Table 1 median and clamped to the group’s reported range. The spread (a
log-SD of one fifth of the log range) is the maintainers’ choice, as the
paper reports only medians and ranges.

``` r

set.seed(2025)
rxode2::rxSetSeed(2025)
n_per_arm <- 200

draw_fii <- function(n, med, lo, hi) {
  pmin(pmax(med * exp(rnorm(n, 0, log(hi / lo) / 5)), lo), hi)
}

cohort <- bind_rows(
  tibble(id = seq_len(n_per_arm), condition = "Normal", DIS_TRAUMA = 0,
         PROTHROMBIN_BASE = draw_fii(n_per_arm, 1108.5, 725.1, 1338.06)),
  tibble(id = n_per_arm + seq_len(n_per_arm), condition = "Trauma", DIS_TRAUMA = 1,
         PROTHROMBIN_BASE = draw_fii(n_per_arm, 1052.8, 474.09, 1784.83))
)

ev <- tidyr::crossing(cohort, time = seq(0, 25, by = 0.25)) |>
  mutate(evid = 0, cmt = "thrombin") |>
  arrange(id, time)

sim <- rxode2::rxSolve(mod, ev, keep = "condition", returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalksyn_thrombin_tf_trauma, etalkdeg_thrombin_trauma, etalktr_trauma, etalksyn_thrombin_effect_trauma
#> as a work-around try putting the mu-referenced expression on a simple line
stopifnot(nrow(sim) == nrow(ev), !anyNA(sim$thrombin))
```

## Visual predictive check (Figure S2)

The lines show the 5th, 50th and 95th percentiles of the simulated
observations (with residual error) at each time. The points mark the
observed percentiles of the Figure S2 VPC at their peaks, digitised by
the maintainers from the supplement image.

``` r

vpc <- sim |>
  group_by(condition, time) |>
  summarise(p05 = quantile(sim, 0.05), p50 = median(sim), p95 = quantile(sim, 0.95),
            .groups = "drop")

fig_s2 <- tribble(
  ~condition, ~stat, ~time, ~value,
  "Normal", "p50", 8.5, 105,
  "Normal", "p95", 5.5, 195,
  "Normal", "p05", 8.5, 60,
  "Trauma", "p50", 5.0, 197,
  "Trauma", "p95", 5.0, 360,
  "Trauma", "p05", 7.5, 50
)
```

``` r

vpc_long <- vpc |>
  pivot_longer(c(p05, p50, p95), names_to = "stat", values_to = "value")
ggplot(vpc_long, aes(time, value, linetype = stat)) +
  geom_line(aes(colour = condition)) +
  geom_point(data = fig_s2, aes(shape = stat), size = 2.5) +
  facet_wrap(~condition, scales = "free_y") +
  scale_colour_manual(values = c(Normal = "#d7301f", Trauma = "#1b9e9e"), guide = "none") +
  scale_linetype_manual(values = c(p05 = "dashed", p50 = "solid", p95 = "dashed")) +
  labs(x = "Time (min)", y = "FIIa (nM)", linetype = "Percentile", shape = "Figure S2") +
  theme_bw()
```

![Replicates the VPC panels of Figure S2 of Ruiz 2025 (lower half).
Lines: simulated percentiles; points: observed percentiles digitised
from Figure
S2.](Ruiz_2025_thrombin_generation_files/figure-html/vpc-plot-1.png)

Replicates the VPC panels of Figure S2 of Ruiz 2025 (lower half). Lines:
simulated percentiles; points: observed percentiles digitised from
Figure S2.

``` r

vpc_cmp <- vpc_long |>
  group_by(condition, stat) |>
  summarise(sim_peak = max(value), sim_time = time[which.max(value)], .groups = "drop") |>
  inner_join(fig_s2 |> rename(fig_value = value, fig_time = time),
             by = c("condition", "stat")) |>
  mutate(pct_diff = 100 * (sim_peak / fig_value - 1))
stopifnot(nrow(vpc_cmp) == 6)
vpc_cmp |>
  select(condition, stat, sim_peak, fig_value, pct_diff, sim_time, fig_time) |>
  rename(
    "Condition" = condition, "Percentile" = stat,
    "Simulated peak (nM)" = sim_peak, "Figure S2 peak (nM)" = fig_value,
    "Difference (%)" = pct_diff, "Simulated time (min)" = sim_time,
    "Figure S2 time (min)" = fig_time
  ) |>
  knitr::kable(digits = 1)
```

| Condition | Percentile | Simulated peak (nM) | Figure S2 peak (nM) | Difference (%) | Simulated time (min) | Figure S2 time (min) |
|:---|:---|---:|---:|---:|---:|---:|
| Normal | p05 | 41.7 | 60 | -30.5 | 6.2 | 8.5 |
| Normal | p50 | 116.9 | 105 | 11.3 | 7.5 | 8.5 |
| Normal | p95 | 212.2 | 195 | 8.8 | 7.0 | 5.5 |
| Trauma | p05 | 55.3 | 50 | 10.6 | 6.0 | 7.5 |
| Trauma | p50 | 184.4 | 197 | -6.4 | 5.5 | 5.0 |
| Trauma | p95 | 447.9 | 360 | 24.4 | 4.2 | 5.0 |

The simulated medians sit within about 15% of the observed medians in
both groups. The extreme percentiles depend on the assumed spread of
factor II and on the very large IIV of `ksynIIa` in normal subjects (CV
240%), and the Figure S2 percentiles are computed in coarse time bins
that blunt the peaks, so they are reported but checked only loosely.

``` r

med_cmp <- vpc_cmp |> filter(stat == "p50")
stopifnot(
  # Median curve: a mis-transcribed rate constant moves this by far more.
  nrow(med_cmp) == 2,
  all(abs(med_cmp$pct_diff) < 30),
  # The trauma median peaks earlier than the normal median.
  med_cmp$sim_time[med_cmp$condition == "Trauma"] <
    med_cmp$sim_time[med_cmp$condition == "Normal"]
)
```

## Individual profiles (Figure 3)

``` r

show_ids <- c(seq_len(12), n_per_arm + seq_len(12))
sim |>
  filter(id %in% show_ids) |>
  ggplot(aes(time)) +
  geom_point(aes(y = sim, colour = condition), size = 0.4) +
  geom_line(aes(y = thrombin), colour = "grey30") +
  facet_wrap(~id, ncol = 6) +
  scale_colour_manual(values = c(Normal = "#d7301f", Trauma = "#1b9e9e")) +
  labs(x = "Time (min)", y = "FIIa (nM)", colour = NULL) +
  theme_bw() +
  theme(strip.text = element_blank())
```

![Twelve simulated individual curves per condition, in the layout of
Figure 3 of Ruiz 2025 (individual predictions as lines, simulated
observations as
points).](Ruiz_2025_thrombin_generation_files/figure-html/indiv-plot-1.png)

Twelve simulated individual curves per condition, in the layout of
Figure 3 of Ruiz 2025 (individual predictions as lines, simulated
observations as points).

## Thrombin-generation parameters (PKNCA)

The paper reports no non-compartmental summaries, so this table is
descriptive. Peak height (`cmax`), time to peak (`tmax`) and the area
under the curve over the 25-min analysis window (`auclast`, the
truncated endogenous thrombin potential) are computed from the
individual predictions with PKNCA. The individual peaks in Figure 3 of
the paper span roughly 70-200 nM in normal subjects and 150-480 nM in
trauma patients.

``` r

conc <- sim |>
  filter(!is.na(thrombin)) |>
  select(id, condition, time, thrombin)
conc_obj <- PKNCA::PKNCAconc(conc, thrombin ~ time | condition + id)
intervals <- data.frame(start = 0, end = 25, cmax = TRUE, tmax = TRUE, auclast = TRUE)
nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, intervals = intervals))
#> No dose information provided, calculations requiring dose will return NA.
nca_tab <- as.data.frame(nca$result) |>
  group_by(condition, PPTESTCD) |>
  summarise(median = median(PPORRES), p10 = quantile(PPORRES, 0.1),
            p90 = quantile(PPORRES, 0.9), .groups = "drop")
stopifnot(nrow(nca_tab) == 6)
nca_tab |>
  mutate(PPTESTCD = recode(PPTESTCD, cmax = "Peak height (nM)",
                           tmax = "Time to peak (min)",
                           auclast = "AUC 0-25 min (nM x min)")) |>
  rename("Condition" = condition, "Parameter" = PPTESTCD, "Median" = median,
         "10th percentile" = p10, "90th percentile" = p90) |>
  knitr::kable(digits = 1)
```

| Condition | Parameter               | Median | 10th percentile | 90th percentile |
|:----------|:------------------------|-------:|----------------:|----------------:|
| Normal    | AUC 0-25 min (nM x min) | 1074.9 |           769.2 |          1419.0 |
| Normal    | Peak height (nM)        |  130.5 |            90.2 |           203.0 |
| Normal    | Time to peak (min)      |    7.2 |             5.5 |             8.8 |
| Trauma    | AUC 0-25 min (nM x min) | 1277.0 |           832.3 |          2156.4 |
| Trauma    | Peak height (nM)        |  226.9 |           114.0 |           419.0 |
| Trauma    | Time to peak (min)      |    5.5 |             4.0 |             7.8 |

## Assumptions and deviations

- **Tissue factor input.** The Methods give the recombinant TF as 5 pM;
  it is encoded as `tf0 = 0.005` nM to match the model’s nM units. A
  typical-value solve with this input reproduces the Figure S2 median
  curves, whereas reading 5 pM as 5 nM exhausts all factor II within
  half a minute.
- **Sign of the TF equation.** Figure 2 prints `dTF/dt = kdegTF * TF`,
  while Equation 1 has `-kdegTF * TF`. The Equation 1 sign is used; it
  makes no difference to any simulation because `kdegTF` is fixed to
  zero.
- **Proportional residual error.** Table 2 lists ‘Proportinal error (%)
  0.18’. The value is read as a fraction (18%), not 0.18%: the
  observation-versus-individual-prediction scatter of Figure S2 spreads
  by roughly +/-50 nM around predictions of 200 nM, which an 18% term
  explains and a 0.18% term (about 0.4 nM) does not.
- **Residual error switch.** The Results describe a combined error model
  ‘up to 15 min’ and an additive one ‘beyond 15 min’; observations at
  exactly 15 min use the additive model. The two additive SDs are
  separate parameters that both round to 2.4 nM in Table 2.
- **Inter-individual variability.** Table 2 reports IIV as CV% with the
  footnote CV = sqrt(omega^2) x 100, so each variance is (CV/100)^2.
  Each condition has its own eta for every rate constant; no covariances
  were estimated. The `ksynIIa` IIV in normal subjects (240%, shrinkage
  63%) is very large and dominates the spread of the normal-subject
  curves.
- **Thrombin feedback term.** The `(1 + FIIa)` factor of Figure 2 is
  kept exactly as printed, with FIIa in nM (see Units).
- **Covariate naming.** Subject condition is carried by `DIS_TRAUMA` (1
  = trauma patient, 0 = normal subject) and baseline factor II by
  `PROTHROMBIN_BASE` (nM). Table 1 reports the other coagulation
  factors, but none entered the final model.
- **Virtual cohort.** Factor II is drawn log-normally around each
  group’s Table 1 median and clamped to its range; the spread is assumed
  because only medians and ranges are reported. The observed Figure S2
  percentiles were digitised from the supplement image and are
  approximate.
- **Errata.** No erratum or correction was found for this article
  (checked 2026-11-09).
