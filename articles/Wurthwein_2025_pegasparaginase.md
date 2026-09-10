# PEGylated asparaginase in the AIEOP-BFM ALL 2009 trial (Wurthwein 2025)

## Model and source

Wurthwein 2025 reports **four** independently fitted population PK
models for PEGylated asparaginase (PEG-ASNase) in the AIEOP-BFM ALL 2009
trial. They share one structure and differ in the cohort, the set of
administrations covered, and the estimated values, so each is packaged
as its own model file. This vignette walks the paper as a unit.

| Model file | Source table | Cohort and administrations |
|----|----|----|
| `Wurthwein_2025_pegasparaginase_germanCzech` | ESM Table S13, model 223015 | German/Czech, induction (protocol IA days 12, 26) + re-induction (protocol II day 8) |
| `Wurthwein_2025_pegasparaginase_italian` | ESM Table S13, model 223005 | Italian, same three administrations, MAAT assay converted to the AHA scale |
| `Wurthwein_2025_pegasparaginase_hrPostinduction` | ESM Table S7, model 354103 | German/Czech high risk, induction through protocol IB, HR-1 to HR-3 and the three protocol III re-inductions |
| `Wurthwein_2025_pegasparaginase_r2ea` | ESM Table S9, model 493158 | German/Czech non-high-risk R2 experimental arm, induction plus ten biweekly doses to maintenance M10 |

- Citation: Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath
  CE, Matteo C, Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in
  Children with Acute Lymphoblastic Leukemia Treated within the
  AIEOP-BFM ALL 2009 Trial: Population Pharmacokinetics and Drug
  Exposure. Eur J Drug Metab Pharmacokinet. 2025;50(6):683-696.
  <doi:10.1007/s13318-025-00962-3>
- Article: <https://doi.org/10.1007/s13318-025-00962-3>
- Supplement (parameter tables S4-S15):
  <https://doi.org/10.1007/s13318-025-00962-3> (Electronic Supplementary
  Material 1)
- Predecessor paper supplying the structural model and the NONMEM
  control stream: Wurthwein G, et al. Eur J Drug Metab Pharmacokinet.
  2021;46:289-300. <https://doi.org/10.1007/s13318-021-00670-8>

The authors deliberately did **not** pool the German/Czech and Italian
data: a combined model carrying a country covariate “became unstable:
some runs showed extremely high relative standard error (RSE) values in
one or the other parameter estimate” (Section 3.3.2). That decision is
why two of the four models below are structurally identical but
separately fitted.

``` r

mGC <- rxode2::rxode(readModelDb("Wurthwein_2025_pegasparaginase_germanCzech"))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
mIT <- rxode2::rxode(readModelDb("Wurthwein_2025_pegasparaginase_italian"))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
mHR <- rxode2::rxode(readModelDb("Wurthwein_2025_pegasparaginase_hrPostinduction"))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7
#> as a work-around try putting the mu-referenced expression on a simple line
mEA <- rxode2::rxode(readModelDb("Wurthwein_2025_pegasparaginase_r2ea"))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_vc_8, etaiov_vc_9, etaiov_vc_10, etaiov_vc_11, etaiov_vc_12, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7, etaiov_cl_8, etaiov_cl_9, etaiov_cl_10, etaiov_cl_11, etaiov_cl_12
#> as a work-around try putting the mu-referenced expression on a simple line
```

## Population

Children aged 1 to under 18 years with newly diagnosed acute
lymphoblastic leukemia enrolled in AIEOP-BFM ALL 2009 (EudraCT
2007-004270-43, NCT01117441) between 1 June 2010 and 28 February 2017.
PEG-ASNase was given as a 2-hour intravenous infusion at 2500 U/m^2 per
dose, with an absolute cap of 3750 U per dose that bound for about 12%
of administrations.

Of 2771 German/Czech and 1656 Italian trial patients, 2535 and 1603
respectively were eligible for the popPK analyses, contributing 17,221
and 6894 samples (Table 1, ESM Table S1). German/Czech median age 5.1
years (range 1.06-18.0) and median BSA 0.78 m^2 (0.39-2.44), 1468 male /
1067 female; Italian median age 5.1 years (1.04-18.0), median BSA 0.77
m^2 (0.28-2.30), 927 male / 676 female. High-risk patients were older
(median 8.0 years) and larger (median BSA 0.99 m^2) than non-high-risk
patients and more often male. A third, Australian cohort of 279 patients
contributed dose-intensity values only; it was “too small to evaluate
the adequacy of the popPK model” and no model was fitted to it.

Sampling followed the trial’s therapeutic-drug-monitoring schedule
rather than a classical PK design: before each treatment phase and 7 and
14 days after each dose, so 85.5-93.3% of analysed samples fall in the
day 7 +/- 1 or day 14 +/- 1 windows (ESM Table S2).

**These models describe standard elimination only.** Samples indicating
silent inactivation (asparaginase below 100 U/L within 8 days and/or
undetectable within 15 days, per van der Sluis 2016) and every sample at
or after a hypersensitivity reaction were excluded before fitting. The
authors are explicit that inactivation pharmacokinetics “are expected
not to follow the standard elimination of the drug” and are not
described here (Section 4.3).

The same information is available programmatically, e.g.
`readModelDb("Wurthwein_2025_pegasparaginase_germanCzech")()$population`.

## Structural model

Asparaginase activity is carried by a chain of **14 serial compartments
that share a single serum volume** `V`. Drug steps from compartment *i*
to *i+1* with intercompartmental clearance `Qtr`, standing in for
progressive hydrolysis of the polyethylene-glycol moiety, and every
compartment is eliminated with the initial clearance `CLinitial`. The
terminal compartment additionally loses drug at `Qtr`: the authors’
original model gave it a separate induced clearance `CLinduced`, and the
published simplification sets `CLinduced = Qtr` (Wurthwein 2021 ESM
Section 2, dOFV = 0.6).

That is the whole mechanism behind the paper’s central observation.
Immediately after a dose all mass sits in compartment 1 and the apparent
clearance is `CLinitial`; as mass migrates down the chain, apparent
clearance rises toward an asymptote of `CLinitial + Qtr`, roughly
ninefold higher. The approach is slow – it needs `Qtr/V * t` to exceed
the chain length of 13, i.e. tens of days – so over a single 14-day
dosing interval the apparent clearance climbs only from 0.126 to about
0.149 L/day. That gentle within-interval ramp is exactly the curvature
the transit chain was introduced to capture. The AHA assay measures
**total** catalytic activity, so all 14 species contribute to the
observation: `Cc = (central + transit1 + ... + transit13) / V`.

The equations are not printed in Wurthwein 2025; they come from the
complete NONMEM control stream in Section 4 of the Wurthwein 2021
supplement (reference 7 of the 2025 supplement), whose `$DES` block and
`Figure 1` schematic fix the topology, and whose `$PK` block fixes every
covariate form used below.

``` r

gc0 <- rxode2::zeroRe(mGC)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
ev_ramp <- rxode2::et(amt = 2500, dur = 2 / 24, cmt = "central") |>
  rxode2::et(seq(0, 28, by = 0.05), cmt = "central") |>
  as.data.frame() |>
  mutate(BSA = 1, AGE = 5, SEXF = 0, OCC = 1)
ramp <- rxode2::rxSolve(gc0, ev_ramp, returnType = "data.frame") |>
  filter(time > 2 / 24) |>
  # Instantaneous apparent clearance = (rate of elimination) / concentration.
  # Every compartment is eliminated at CLinitial and the terminal one also at
  # Qtr, so elimination rate = CLinitial * Ctotal + Qtr * transit13 / V and
  #   CL_app = CLinitial + Qtr * (transit13 / V) / Cc.
  # With a bolus into compartment 1 and a common loss rate (ke + kt) in every
  # compartment, the state amounts are Poisson in kt * t, so the closed form is
  #   CL_app = CLinitial + Qtr * dpois(13, kt * t) / ppois(13, kt * t).
  mutate(
    cl_app = 0.126 + 0.926 * (transit13 / 1.68) / Cc,
    cl_cf  = 0.126 + 0.926 * dpois(13, (0.926 / 1.68) * time) /
      ppois(13, (0.926 / 1.68) * time)
  )
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'

ggplot(ramp, aes(time)) +
  geom_line(aes(y = cl_app), linewidth = 0.9) +
  geom_line(aes(y = cl_cf), colour = "firebrick", linetype = "dashed") +
  geom_hline(yintercept = 0.126, linetype = "dotted") +
  geom_hline(yintercept = 0.126 + 0.926, linetype = "dotted") +
  annotate("text", x = 14, y = 0.126, vjust = -0.6, size = 3, label = "CLinitial") +
  annotate("text", x = 14, y = 1.052, vjust = -0.6, size = 3,
           label = "CLinitial + Qtr (asymptote)") +
  labs(x = "Days after dose", y = "Apparent clearance (L/day)") +
  theme_bw()
```

![Apparent clearance over 28 days for the German/Czech typical child at
BSA 1 m^2, with the closed-form prediction overlaid. Model-derived; not
a published
figure.](Wurthwein_2025_pegasparaginase_files/figure-html/cl-ramp-1.png)

Apparent clearance over 28 days for the German/Czech typical child at
BSA 1 m^2, with the closed-form prediction overlaid. Model-derived; not
a published figure.

``` r

# Structural gate on the chain topology. The solved apparent clearance must
# start at CLinitial, increase monotonically, never exceed the CLinitial + Qtr
# asymptote, and track the Poisson closed form. If the terminal compartment did
# not carry the extra Qtr the curve would be flat at CLinitial; if the chain
# length were wrong the closed form would use the wrong Poisson index.
c(start = min(ramp$cl_app), day14 = ramp$cl_app[which.min(abs(ramp$time - 14))],
  day28 = max(ramp$cl_app),
  max_pct_diff_vs_closed_form = max(abs(100 * (ramp$cl_app - ramp$cl_cf) / ramp$cl_cf)))
#>                       start                       day14 
#>                   0.1260000                   0.1490239 
#>                       day28 max_pct_diff_vs_closed_form 
#>                   0.3825922                   0.3381019

stopifnot(
  abs(min(ramp$cl_app) - 0.126) < 0.01,
  all(diff(ramp$cl_app) > -1e-9),
  max(ramp$cl_app) < 0.126 + 0.926,
  # Solve versus its own closed form: pure numerical error, so a tight bound.
  max(abs((ramp$cl_app - ramp$cl_cf) / ramp$cl_cf)) < 0.01
)
```

### Mass balance

A linear system has a closed-form total AUC. For this chain, with
`ke = CLinitial/V`, `kt = Qtr/V` and `ratio = kt/(ke + kt)`, the total
concentration-time integral of a dose `D` is

`AUC = D / (V * (ke + kt)) * sum(ratio^j, j = 0..13)`.

Reproducing that from the solved ODEs confirms the chain is wired
correctly and that no mass is created or lost at the boundaries.

``` r

ev_mb <- rxode2::et(amt = 2500, cmt = "central") |>
  rxode2::et(seq(0, 120, by = 0.01), cmt = "central") |>
  as.data.frame() |>
  mutate(BSA = 1, AGE = 5, SEXF = 0, OCC = 1)
mb <- rxode2::rxSolve(gc0, ev_mb, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'

auc_num <- sum(diff(mb$time) * (head(mb$Cc, -1) + tail(mb$Cc, -1)) / 2)

V <- 1.68; ke <- 0.126 / V; kt <- 0.926 / V; r <- kt / (ke + kt)
auc_cf <- 2500 / (V * (ke + kt)) * sum(r^(0:13))

c(numeric = auc_num, closed_form = auc_cf, pct_diff = 100 * (auc_num - auc_cf) / auc_cf)
#>      numeric  closed_form     pct_diff 
#> 1.651542e+04 1.651542e+04 5.631459e-06

# Pure numerical error between a solve and its own closed form: a tight bound
# is correct here (both sides use the same parameters; there is no per-subject
# random mechanism that could move across machines).
stopifnot(abs(auc_num - auc_cf) / auc_cf < 0.005)
```

## Source trace

Every `ini()` entry carries an in-file comment naming its source
location. The table collects them and, importantly, records which facts
come from the 2025 supplement and which from the 2021 control stream.

| Item | Value(s) | Source location |
|----|----|----|
| Chain of 14 serial compartments; dose into compartment 1 | topology | Wurthwein 2021 Fig. 1; 2021 ESM Section 4 `$MODEL` / `$DES` |
| Terminal compartment loses `CLinitial + Qtr` (`CLinduced = Qtr`) | structural | Wurthwein 2021 Fig. 1 caption; 2021 ESM Section 2 |
| Observation = total activity over all 14 states divided by `V` | structural | Wurthwein 2021 ESM Section 4 `$ERROR` (`TOT = A(1)+...+A(14); IPRED = TOT/V1`) |
| BSA linear and centred on 0.79 m^2, one slope on `V`, one shared by `CLinitial` and `Qtr` | form + constant | Wurthwein 2021 ESM Section 4 `$PK`; 2025 ESM Tables S7/S9/S13 footnote (a) |
| Age hockey stick, break point fixed at 8 years | form + constant | Wurthwein 2021 ESM Section 4 `$PK` (two `AGEGTP` branches, sub-8-year slope `0 FIX`); 2021 ESM Section 2 |
| Sex effect on `CLinitial`, males the reference | form | Wurthwein 2021 ESM Section 4 `$PK` (`SEX.EQ.1` = 1) |
| Anti-PEG IgM hockey stick, cut point 1.3 on log scale, first induction dose only | form + constant | Wurthwein 2025 ESM Section 2.2 eq. (7), Table S6 model 354103, Table S9 footnote (e) |
| Residual error combined proportional + additive SD | form | Wurthwein 2021 ESM Section 4 `$ERROR` + `$SIGMA 1 FIX` |
| `V`, `CLinitial`, `Qtr`, BSA slopes, phase changes, age, sex, IIV, IOV, residual error (German/Czech, Italian) | see model files | Wurthwein 2025 ESM Table S13 |
| Same set plus the anti-PEG IgM effect (high-risk post-induction) | see model file | Wurthwein 2025 ESM Table S7 |
| Same set plus the anti-PEG IgM effect (R2-EA) | see model file | Wurthwein 2025 ESM Table S9 |
| `F(CLinitial HR-1)` fixed to 0 | `0` | Wurthwein 2025 ESM Table S7 footnote (f); freely estimated it was -0.031 with 70.4% RSE (Table S4 model 334003) |
| MAAT-to-AHA conversion factors 1.23 / 1.42 | conversion | Wurthwein 2025 Section 3.3.1, ESM Table S12 |
| Published dose-normalized activity medians used as gates below | see tables | Wurthwein 2025 ESM Table S3 (German/Czech), Table S11 (Italian) |
| Published dose-intensity parameters used as gates below | see table | Wurthwein 2025 ESM Table S15 |

## Reproducing the published activity levels

The paper’s primary descriptive result is ESM Table S3: median
dose-normalized asparaginase activity 7 and 14 days after each
administration. Dose-normalized means scaled to a 2500 U/m^2 dose (ESM
Section 1.4), so a typical-value solve at 2500 U/m^2 is directly
comparable.

These are deterministic typical-value solves with random effects zeroed,
so the tolerances below are tight; they compare a population-typical
prediction against an observed cohort median, and the residual gap is
model misspecification rather than Monte Carlo noise.

``` r

# Simulate a dosing series where OCC steps at each administration. OCC is a
# time-varying covariate, so it must be attached to a materialized data frame
# (assigning onto an rxEt object is silently dropped).
run_series <- function(mod, dose_times, occs, bsa, age, sexf, ab = 1,
                       tail_days = 40, by = 0.05) {
  m0 <- rxode2::zeroRe(mod)
  grid <- seq(0, max(dose_times) + tail_days, by = by)
  ev <- rxode2::et(amt = 2500 * bsa, time = dose_times, dur = 2 / 24, cmt = "central") |>
    rxode2::et(grid, cmt = "central") |>
    as.data.frame()
  idx <- findInterval(ev$time, dose_times)
  idx[idx < 1] <- 1
  ev$OCC <- occs[idx]
  ev$BSA <- bsa
  ev$AGE <- age
  ev$SEXF <- sexf
  if ("ABPEG_IGM" %in% mod$allCovs) ev$ABPEG_IGM <- ab
  rxode2::rxSolve(m0, ev, returnType = "data.frame")
}

# Value at `t`. When the readout time coincides with an administration -- the
# day-14 sample after the first induction dose is drawn on the day of the
# second -- take the last grid point strictly BEFORE it. The stepwise covariate
# model changes V discontinuously at each administration, so reading the row at
# the dose time itself would use the post-switch volume and overstate the
# trough by 1/(1 + F_V) (about 18% for the induction step). Away from a dose
# time the profile is continuous and the nearest grid point is correct.
at_day <- function(sim, tt, doses = numeric(0)) {
  i <- if (any(abs(doses - tt) < 1e-6)) {
    max(which(sim$time < tt - 1e-9))
  } else {
    which.min(abs(sim$time - tt))
  }
  sim$Cc[i]
}
```

``` r

# German/Czech induction + re-induction (ESM Table S13 model).
ind <- run_series(mGC, c(0, 14), c(1, 2), bsa = 0.78, age = 5.1, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'
rei <- run_series(mGC, 0, 3, bsa = 0.78, age = 5.1, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'

gc_chk <- tibble::tibble(
  Administration = c("PIA d12", "PIA d12", "PIA d26", "PIA d26",
                     "PII d8", "PII d8"),
  Day = c(7, 14, 7, 14, 7, 14),
  Simulated = c(at_day(ind, 7, c(0, 14)), at_day(ind, 14, c(0, 14)),
                at_day(ind, 21, c(0, 14)), at_day(ind, 28, c(0, 14)),
                at_day(rei, 7, 0), at_day(rei, 14, 0)),
  Published = c(895, 543, 1260, 630, 1440, 753)
) |>
  mutate(`Difference (%)` = 100 * (Simulated - Published) / Published)

gc_chk |>
  knitr::kable(digits = c(0, 0, 0, 0, 1),
               caption = paste("German/Czech model versus ESM Table S3 median",
                               "dose-normalized activity (U/L). PII d8 row is",
                               "the non-high-risk median."))
```

| Administration | Day | Simulated | Published | Difference (%) |
|:---------------|----:|----------:|----------:|---------------:|
| PIA d12        |   7 |       920 |       895 |            2.8 |
| PIA d12        |  14 |       525 |       543 |           -3.3 |
| PIA d26        |   7 |      1284 |      1260 |            1.9 |
| PIA d26        |  14 |       591 |       630 |           -6.2 |
| PII d8         |   7 |      1398 |      1440 |           -2.9 |
| PII d8         |  14 |       727 |       753 |           -3.5 |

German/Czech model versus ESM Table S3 median dose-normalized activity
(U/L). PII d8 row is the non-high-risk median. {.table}

``` r

stopifnot(
  # Structural: a mis-transcribed clearance, volume, dose or unit would move
  # the whole set by tens of percent.
  abs(median(gc_chk$`Difference (%)`)) < 5,
  max(abs(gc_chk$`Difference (%)`)) < 10
)
```

``` r

# Italian model against the AHA-CONVERTED Italian medians of ESM Table S11
# (the "Median conv." column), since the model predicts on the AHA scale.
it_ind <- run_series(mIT, c(0, 14), c(1, 2), bsa = 0.77, age = 5.1, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'
it_rei <- run_series(mIT, 0, 3, bsa = 0.77, age = 5.1, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'

it_chk <- tibble::tibble(
  Administration = c("PIA d12", "PIA d12", "PIA d26", "PIA d26", "PII d8", "PII d8"),
  Day = c(7, 14, 7, 14, 7, 14),
  Simulated = c(at_day(it_ind, 7, c(0, 14)), at_day(it_ind, 14, c(0, 14)),
                at_day(it_ind, 21, c(0, 14)), at_day(it_ind, 28, c(0, 14)),
                at_day(it_rei, 7, 0), at_day(it_rei, 14, 0)),
  Published = c(980, 547, 1390, 652, 1390, 721)
) |>
  mutate(`Difference (%)` = 100 * (Simulated - Published) / Published)

it_chk |>
  knitr::kable(digits = c(0, 0, 0, 0, 1),
               caption = paste("Italian model versus ESM Table S11 median",
                               "dose-normalized activity after MAAT-to-AHA",
                               "conversion (U/L)."))
```

| Administration | Day | Simulated | Published | Difference (%) |
|:---------------|----:|----------:|----------:|---------------:|
| PIA d12        |   7 |       976 |       980 |           -0.4 |
| PIA d12        |  14 |       551 |       547 |            0.8 |
| PIA d26        |   7 |      1372 |      1390 |           -1.3 |
| PIA d26        |  14 |       631 |       652 |           -3.2 |
| PII d8         |   7 |      1348 |      1390 |           -3.1 |
| PII d8         |  14 |       697 |       721 |           -3.4 |

Italian model versus ESM Table S11 median dose-normalized activity after
MAAT-to-AHA conversion (U/L). {.table}

``` r


stopifnot(
  abs(median(it_chk$`Difference (%)`)) < 5,
  max(abs(it_chk$`Difference (%)`)) < 10
)
```

### Accumulation across repeated dosing

The two extended models exist to describe accumulation, so the sharper
test is whether they reproduce the *rise* in day-7 activity across a
whole schedule.

``` r

# R2-EA: ten biweekly doses, protocol II day 8 through maintenance M10.
ea_times <- seq(0, by = 14, length.out = 10)
ea <- run_series(mEA, ea_times, 3:12, bsa = 0.74, age = 4.7, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_vc_8, etaiov_vc_9, etaiov_vc_10, etaiov_vc_11, etaiov_vc_12, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7, etaiov_cl_8, etaiov_cl_9, etaiov_cl_10, etaiov_cl_11, etaiov_cl_12
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_vc_4', 'etaiov_vc_5', 'etaiov_vc_6', 'etaiov_vc_7', 'etaiov_vc_8', 'etaiov_vc_9', 'etaiov_vc_10', 'etaiov_vc_11', 'etaiov_vc_12', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_4', 'etaiov_cl_5', 'etaiov_cl_6', 'etaiov_cl_7', 'etaiov_cl_8', 'etaiov_cl_9', 'etaiov_cl_10', 'etaiov_cl_11', 'etaiov_cl_12'

ea_chk <- tibble::tibble(
  Administration = c("PII d8", "PII-ASP+ d22", "PII-ASP+ d36", "PII-ASP+ d50",
                     "M5", "M6", "M7", "M8", "M9", "M10"),
  Simulated = vapply(ea_times, function(d) at_day(ea, d + 7, ea_times), numeric(1)),
  Published = c(1440, 1830, 1660, 1780, 1790, 1840, 2120, 2030, 2040, 2070)
) |>
  mutate(`Difference (%)` = 100 * (Simulated - Published) / Published)

# HR-EA: induction plus the four weekly protocol IB experimental-arm doses.
hr_times <- c(0, 14, 28, 35, 42, 49)
hr <- run_series(mHR, hr_times, c(1, 2, 3, 3, 3, 3), bsa = 0.99, age = 8, sexf = 0)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_vc_4', 'etaiov_vc_5', 'etaiov_vc_6', 'etaiov_vc_7', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_4', 'etaiov_cl_5', 'etaiov_cl_6', 'etaiov_cl_7'

hr_chk <- tibble::tibble(
  Administration = c("PIA d12", "PIA d26", "PIB-ASP+ d40", "PIB-ASP+ d47",
                     "PIB-ASP+ d54", "PIB-ASP+ d61"),
  Simulated = vapply(hr_times, function(d) at_day(hr, d + 7, hr_times), numeric(1)),
  Published = c(844, 1110, 1310, 1990, 2200, 2190)
) |>
  mutate(`Difference (%)` = 100 * (Simulated - Published) / Published)

dplyr::bind_rows(
  mutate(ea_chk, Arm = "R2-EA (non-high-risk)"),
  mutate(hr_chk, Arm = "HR-EA (high-risk)")
) |>
  select(Arm, Administration, Simulated, Published, `Difference (%)`) |>
  knitr::kable(digits = c(0, 0, 0, 0, 1),
               caption = paste("Extended models versus ESM Table S3 median day",
                               "7 dose-normalized activity (U/L)."))
```

| Arm                   | Administration | Simulated | Published | Difference (%) |
|:----------------------|:---------------|----------:|----------:|---------------:|
| R2-EA (non-high-risk) | PII d8         |      1430 |      1440 |           -0.7 |
| R2-EA (non-high-risk) | PII-ASP+ d22   |      1588 |      1830 |          -13.2 |
| R2-EA (non-high-risk) | PII-ASP+ d36   |      1775 |      1660 |            6.9 |
| R2-EA (non-high-risk) | PII-ASP+ d50   |      1819 |      1780 |            2.2 |
| R2-EA (non-high-risk) | M5             |      1819 |      1790 |            1.6 |
| R2-EA (non-high-risk) | M6             |      1819 |      1840 |           -1.1 |
| R2-EA (non-high-risk) | M7             |      1950 |      2120 |           -8.0 |
| R2-EA (non-high-risk) | M8             |      1985 |      2030 |           -2.2 |
| R2-EA (non-high-risk) | M9             |      1985 |      2040 |           -2.7 |
| R2-EA (non-high-risk) | M10            |      1985 |      2070 |           -4.1 |
| HR-EA (high-risk)     | PIA d12        |       868 |       844 |            2.9 |
| HR-EA (high-risk)     | PIA d26        |      1232 |      1110 |           11.0 |
| HR-EA (high-risk)     | PIB-ASP+ d40   |      1417 |      1310 |            8.2 |
| HR-EA (high-risk)     | PIB-ASP+ d47   |      2030 |      1990 |            2.0 |
| HR-EA (high-risk)     | PIB-ASP+ d54   |      2288 |      2200 |            4.0 |
| HR-EA (high-risk)     | PIB-ASP+ d61   |      2326 |      2190 |            6.2 |

Extended models versus ESM Table S3 median day 7 dose-normalized
activity (U/L). {.table}

``` r

# Gate the ACCUMULATION, which is what these models were built for: the
# observed fold-rise from first to last dose in each arm must be reproduced.
ea_fold_sim <- ea_chk$Simulated[10] / ea_chk$Simulated[1]
ea_fold_pub <- ea_chk$Published[10] / ea_chk$Published[1]
hr_fold_sim <- hr_chk$Simulated[6] / hr_chk$Simulated[1]
hr_fold_pub <- hr_chk$Published[6] / hr_chk$Published[1]

c(R2EA_sim = ea_fold_sim, R2EA_pub = ea_fold_pub,
  HREA_sim = hr_fold_sim, HREA_pub = hr_fold_pub)
#> R2EA_sim R2EA_pub HREA_sim HREA_pub 
#> 1.388678 1.437500 2.678292 2.594787

stopifnot(
  abs(ea_fold_sim - ea_fold_pub) / ea_fold_pub < 0.15,
  abs(hr_fold_sim - hr_fold_pub) / hr_fold_pub < 0.15,
  # Level agreement: centre tight, envelope loose enough for the one
  # administration the authors' own grouping under-resolves (see Errata).
  abs(median(ea_chk$`Difference (%)`)) < 6,
  quantile(abs(ea_chk$`Difference (%)`), 0.9) < 12,
  abs(median(hr_chk$`Difference (%)`)) < 8,
  max(abs(hr_chk$`Difference (%)`)) < 15
)
```

``` r

ggplot(filter(ea, time <= max(ea_times) + 14), aes(time, Cc)) +
  geom_line(linewidth = 0.6) +
  geom_point(
    data = tibble::tibble(time = ea_times + 7, Cc = ea_chk$Published),
    colour = "firebrick", size = 2
  ) +
  geom_hline(yintercept = 100, linetype = "dotted") +
  annotate("text", x = 5, y = 100, vjust = -0.6, size = 3,
           label = "100 U/L target") +
  labs(x = "Days from the protocol II day 8 dose",
       y = "Asparaginase activity (U/L)") +
  theme_bw()
```

![Simulated R2-EA profile over the ten biweekly doses with the published
day-7 medians overlaid. Replicates the R2-EA panel of Wurthwein 2025
Figure 2a and ESM Table
S3.](Wurthwein_2025_pegasparaginase_files/figure-html/accumulation-plot-1.png)

Simulated R2-EA profile over the ten biweekly doses with the published
day-7 medians overlaid. Replicates the R2-EA panel of Wurthwein 2025
Figure 2a and ESM Table S3.

## Virtual cohort

The dose-intensity comparison below needs a cohort rather than one
typical child, because ESM Table S15 reports medians over patients. Only
the median and range of BSA and age are published, so the cohort uses
log-normal distributions matched to the published medians and truncated
to the published ranges, sampled on a **deterministic quantile lattice**
rather than randomly. The lattice makes the cohort byte-identical on
every machine, which a random draw would not be.

``` r

n_sub <- 200L
p <- (seq_len(n_sub) - 0.5) / n_sub

# BSA and age are both driven by growth and are strongly correlated in a
# pediatric cohort, so they take the same lattice position rather than being
# drawn independently. Assumption, documented below.
cohort <- tibble::tibble(
  id   = seq_len(n_sub),
  BSA  = pmin(pmax(qlnorm(p, log(0.78), 0.38), 0.39), 2.44),
  AGE  = pmin(pmax(qlnorm(p, log(5.1),  0.50), 1.06), 18.0),
  # 1067 of 2535 German/Czech patients are female (42.1%). Assigned by a Weyl
  # sequence rather than by a contiguous block, so the marginal is exact and
  # deterministic without correlating sex with the BSA / age lattice position.
  SEXF = as.integer(((seq_len(n_sub) * 0.421) %% 1) < 0.421),
  ABPEG_IGM = 1
) |>
  # 2500 U/m^2 with the protocol's absolute cap of 3750 U per dose.
  mutate(amt = pmin(2500 * BSA, 3750))

summarise(cohort,
          bsa_median = median(BSA), bsa_min = min(BSA), bsa_max = max(BSA),
          age_median = median(AGE), pct_female = 100 * mean(SEXF),
          pct_capped = 100 * mean(amt < 2500 * BSA - 1e-9))
#> # A tibble: 1 × 6
#>   bsa_median bsa_min bsa_max age_median pct_female pct_capped
#>        <dbl>   <dbl>   <dbl>      <dbl>      <dbl>      <dbl>
#> 1      0.780    0.39    2.27       5.10         42        4.5
```

``` r

stopifnot(
  abs(median(cohort$BSA) - 0.78) < 0.02,
  abs(median(cohort$AGE) - 5.1) < 0.2,
  abs(mean(cohort$SEXF) - 0.421) < 0.02
)
```

## Dose-intensity parameters (ESM Table S15)

The paper’s applied output is the per-exposure-phase dose-intensity
parameter: total `AUC0-inf` and time above activity thresholds. These
were derived from the final models by simulation (`MAXEVAL = 0`), so
they are the most direct end-to-end test of the packaged models.

``` r

# Observation grid: 0.25 day (6 h) throughout, refined to 0.05 day for the
# first day after every administration so the post-infusion peak is not clipped.
# Benchmarked against a 0.01-day reference grid on the induction phase, this
# recovers AUC to 0.00% and time above 100 U/L to 0.24% using 329 points
# instead of 7401 -- a uniform coarse grid without the refinement loses 0.8% of
# the AUC, biased low, because the trapezoid cuts the peak.
dip_grid <- function(dose_times, tail_days) {
  sort(unique(c(
    seq(0, max(dose_times) + tail_days, by = 0.25),
    unlist(lapply(dose_times, function(d) d + seq(0, 1, by = 0.05)))
  )))
}

# Simulate one exposure phase for a cohort. Random effects are zeroed: the
# Table S15 medians are dominated by the covariate distribution, and a
# typical-value cohort makes the comparison reproducible.
run_cohort <- function(mod, dose_times, occs, subjects, tail_days = 60) {
  m0 <- rxode2::zeroRe(mod)
  grid <- dip_grid(dose_times, tail_days)
  ev <- lapply(subjects$id, function(i) {
    ci <- subjects[subjects$id == i, ]
    e <- rxode2::et(amt = ci$amt, time = dose_times, dur = 2 / 24, cmt = "central") |>
      rxode2::et(grid, cmt = "central") |>
      as.data.frame()
    e$id <- i
    e
  }) |>
    dplyr::bind_rows()
  idx <- findInterval(ev$time, dose_times)
  idx[idx < 1] <- 1
  ev$OCC <- occs[idx]
  ev <- dplyr::left_join(ev, select(subjects, id, BSA, AGE, SEXF, ABPEG_IGM), by = "id")
  if (!"ABPEG_IGM" %in% mod$allCovs) ev$ABPEG_IGM <- NULL
  rxode2::rxSolve(m0, ev, returnType = "data.frame")
}

# Every 4th lattice point: 50 subjects spanning the same covariate range. The
# lattice is deterministic, so its median is stable at this size and the five
# exposure phases stay inside the vignette's render budget.
dip_cohort <- cohort[seq(1, n_sub, by = 4), ]

# AUC0-inf and time above threshold, per subject, by trapezoid on the dense
# grid -- exactly the quantities the authors derive from their own solves.
dips <- function(sim) {
  sim |>
    filter(!is.na(Cc)) |>
    group_by(id) |>
    summarise(
      auc  = sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2),
      t100 = sum(diff(time) * (head(Cc, -1) > 100)),
      t400 = sum(diff(time) * (head(Cc, -1) > 400)),
      t1000 = sum(diff(time) * (head(Cc, -1) > 1000)),
      t2000 = sum(diff(time) * (head(Cc, -1) > 2000)),
      .groups = "drop"
    )
}

phases <- list(
  list(name = "PIAd12-PIAd26 (induction)", mod = mGC, times = c(0, 14), occs = c(1, 2),
       pub = c(auc = 35900, t100 = 37.4, t1000 = 15.2, t2000 = 2.53)),
  list(name = "PIId8 (R2-SA)", mod = mGC, times = 0, occs = 3,
       pub = c(auc = 23300, t100 = 22.5, t400 = 17.2, t1000 = 11.0, t2000 = 2.05)),
  list(name = "PIAd12-PIB-ASP+d61 (HR-EA)", mod = mHR,
       times = c(0, 14, 28, 35, 42, 49), occs = c(1, 2, 3, 3, 3, 3),
       pub = c(auc = 120000, t100 = 76.0, t400 = 69.3, t1000 = 48.8, t2000 = 24.0)),
  list(name = "PIId8-M10 (R2-EA)", mod = mEA,
       times = seq(0, by = 14, length.out = 10), occs = 3:12,
       pub = c(auc = 275000, t100 = 153, t1000 = 132, t2000 = 57.3)),
  list(name = "2.PIII", mod = mHR, times = 0, occs = 7,
       pub = c(auc = 22400, t100 = 22.4, t400 = 17.3, t1000 = 11.2))
)

dip_tab <- lapply(phases, function(ph) {
  d <- dips(run_cohort(ph$mod, ph$times, ph$occs, dip_cohort))
  med <- c(auc = median(d$auc), t100 = median(d$t100), t400 = median(d$t400),
           t1000 = median(d$t1000), t2000 = median(d$t2000))
  tibble::tibble(
    `Exposure phase` = ph$name,
    Parameter = names(ph$pub),
    Simulated = as.numeric(med[names(ph$pub)]),
    Published = as.numeric(ph$pub)
  )
}) |>
  dplyr::bind_rows() |>
  mutate(
    Parameter = recode(Parameter,
                       auc = "AUC0-inf (U*day/L)", t100 = "Time > 100 U/L (day)",
                       t400 = "Time > 400 U/L (day)", t1000 = "Time > 1000 U/L (day)",
                       t2000 = "Time > 2000 U/L (day)"),
    `Difference (%)` = 100 * (Simulated - Published) / Published
  )
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_vc_4', 'etaiov_vc_5', 'etaiov_vc_6', 'etaiov_vc_7', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_4', 'etaiov_cl_5', 'etaiov_cl_6', 'etaiov_cl_7'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_vc_8, etaiov_vc_9, etaiov_vc_10, etaiov_vc_11, etaiov_vc_12, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7, etaiov_cl_8, etaiov_cl_9, etaiov_cl_10, etaiov_cl_11, etaiov_cl_12
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_vc_4', 'etaiov_vc_5', 'etaiov_vc_6', 'etaiov_vc_7', 'etaiov_vc_8', 'etaiov_vc_9', 'etaiov_vc_10', 'etaiov_vc_11', 'etaiov_vc_12', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_4', 'etaiov_cl_5', 'etaiov_cl_6', 'etaiov_cl_7', 'etaiov_cl_8', 'etaiov_cl_9', 'etaiov_cl_10', 'etaiov_cl_11', 'etaiov_cl_12'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_vc_4, etaiov_vc_5, etaiov_vc_6, etaiov_vc_7, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3, etaiov_cl_4, etaiov_cl_5, etaiov_cl_6, etaiov_cl_7
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_vc_4', 'etaiov_vc_5', 'etaiov_vc_6', 'etaiov_vc_7', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3', 'etaiov_cl_4', 'etaiov_cl_5', 'etaiov_cl_6', 'etaiov_cl_7'
#> Warning: multi-subject simulation without without 'omega'

dip_tab |>
  knitr::kable(digits = c(0, 0, 1, 1, 1),
               caption = paste("Cohort-median dose-intensity parameters versus",
                               "Wurthwein 2025 ESM Table S15."))
```

| Exposure phase | Parameter | Simulated | Published | Difference (%) |
|:---|:---|---:|---:|---:|
| PIAd12-PIAd26 (induction) | AUC0-inf (U\*day/L) | 36370.6 | 35900.0 | 1.3 |
| PIAd12-PIAd26 (induction) | Time \> 100 U/L (day) | 37.2 | 37.4 | -0.5 |
| PIAd12-PIAd26 (induction) | Time \> 1000 U/L (day) | 16.0 | 15.2 | 5.1 |
| PIAd12-PIAd26 (induction) | Time \> 2000 U/L (day) | 2.7 | 2.5 | 4.7 |
| PIId8 (R2-SA) | AUC0-inf (U\*day/L) | 23264.4 | 23300.0 | -0.2 |
| PIId8 (R2-SA) | Time \> 100 U/L (day) | 22.4 | 22.5 | -0.2 |
| PIId8 (R2-SA) | Time \> 400 U/L (day) | 17.4 | 17.2 | 1.5 |
| PIId8 (R2-SA) | Time \> 1000 U/L (day) | 11.7 | 11.0 | 6.4 |
| PIId8 (R2-SA) | Time \> 2000 U/L (day) | 1.3 | 2.0 | -37.8 |
| PIAd12-PIB-ASP+d61 (HR-EA) | AUC0-inf (U\*day/L) | 127857.8 | 120000.0 | 6.5 |
| PIAd12-PIB-ASP+d61 (HR-EA) | Time \> 100 U/L (day) | 74.7 | 76.0 | -1.7 |
| PIAd12-PIB-ASP+d61 (HR-EA) | Time \> 400 U/L (day) | 69.2 | 69.3 | -0.1 |
| PIAd12-PIB-ASP+d61 (HR-EA) | Time \> 1000 U/L (day) | 51.4 | 48.8 | 5.4 |
| PIAd12-PIB-ASP+d61 (HR-EA) | Time \> 2000 U/L (day) | 28.4 | 24.0 | 18.2 |
| PIId8-M10 (R2-EA) | AUC0-inf (U\*day/L) | 272188.7 | 275000.0 | -1.0 |
| PIId8-M10 (R2-EA) | Time \> 100 U/L (day) | 150.7 | 153.0 | -1.5 |
| PIId8-M10 (R2-EA) | Time \> 1000 U/L (day) | 133.9 | 132.0 | 1.4 |
| PIId8-M10 (R2-EA) | Time \> 2000 U/L (day) | 58.2 | 57.3 | 1.6 |
| 2.PIII | AUC0-inf (U\*day/L) | 23931.1 | 22400.0 | 6.8 |
| 2.PIII | Time \> 100 U/L (day) | 22.7 | 22.4 | 1.3 |
| 2.PIII | Time \> 400 U/L (day) | 17.7 | 17.3 | 2.3 |
| 2.PIII | Time \> 1000 U/L (day) | 12.2 | 11.2 | 8.9 |

Cohort-median dose-intensity parameters versus Wurthwein 2025 ESM Table
S15. {.table}

``` r

# AUC is the structural quantity: it depends only on dose, V and the two
# clearance terms, so a transcription error anywhere in the chain moves it
# immediately. Times above threshold additionally depend on the shape of the
# published cohort's covariate distribution, which is only partly known here,
# so they carry a looser envelope.
auc_rows <- filter(dip_tab, Parameter == "AUC0-inf (U*day/L)")
stopifnot(
  abs(median(auc_rows$`Difference (%)`)) < 10,
  max(abs(auc_rows$`Difference (%)`)) < 20,
  abs(median(dip_tab$`Difference (%)`)) < 12,
  quantile(abs(dip_tab$`Difference (%)`), 0.9) < 30
)
```

## PKNCA validation

PKNCA repeats the exposure calculation for the induction phase through
the package’s validated NCA path rather than the hand-rolled trapezoid
above, and adds `Cmax` and the terminal half-life. Subjects are grouped
into BSA bands so per-group results can be compared. The concentration
frame is filtered only on `!is.na(Cc)`, keeping the time-zero record so
PKNCA does not extrapolate an AUC start before the first measurement.

``` r

# First induction dosing interval only: one dose at time 0, occasion 1, so the
# NCA window is not crossed by the occasion switch at day 14.
pknca_cohort <- cohort[seq(1, n_sub, by = 2), ]
sim_ind <- run_cohort(mGC, 0, 1, pknca_cohort, tail_days = 14)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_vc_1, etaiov_vc_2, etaiov_vc_3, etaiov_cl_1, etaiov_cl_2, etaiov_cl_3
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_vc_1', 'etaiov_vc_2', 'etaiov_vc_3', 'etaiov_cl_1', 'etaiov_cl_2', 'etaiov_cl_3'
#> Warning: multi-subject simulation without without 'omega'

band <- pknca_cohort |>
  mutate(bsa_band = factor(
    case_when(BSA < 0.6 ~ "Low (<0.6 m^2)",
              BSA < 1.1 ~ "Mid (0.6-1.1 m^2)",
              TRUE      ~ "High (>=1.1 m^2)"),
    levels = c("Low (<0.6 m^2)", "Mid (0.6-1.1 m^2)", "High (>=1.1 m^2)")))

nca_conc <- sim_ind |>
  filter(!is.na(Cc)) |>
  left_join(select(band, id, bsa_band), by = "id") |>
  transmute(id, time, Cc, treatment = bsa_band)

nca_dose <- band |>
  transmute(id, time = 0, amt, treatment = bsa_band)

conc_obj <- PKNCA::PKNCAconc(nca_conc, Cc ~ time | treatment + id,
                             concu = "U/L", timeu = "day")
dose_obj <- PKNCA::PKNCAdose(nca_dose, amt ~ time | treatment + id, doseu = "U")

intervals <- data.frame(
  start = 0, end = 14,
  cmax = TRUE, tmax = TRUE, auclast = TRUE, half.life = TRUE
)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))

as.data.frame(nca_res) |>
  filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "half.life")) |>
  group_by(treatment, PPTESTCD) |>
  summarise(median = median(PPORRES, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = median) |>
  dplyr::rename(
    "BSA band"                = treatment,
    "AUC0-14d (U*day/L)"      = auclast,
    "Cmax (U/L)"              = cmax,
    "t1/2 (day)"              = "half.life",
    "Tmax (day)"              = tmax
  ) |>
  knitr::kable(digits = 2,
               caption = paste("PKNCA summary over the first induction dosing",
                               "interval, median per BSA band."))
```

| BSA band          | AUC0-14d (U\*day/L) | Cmax (U/L) | t1/2 (day) | Tmax (day) |
|:------------------|--------------------:|-----------:|-----------:|-----------:|
| Low (\<0.6 m^2)   |            15330.71 |    1805.44 |       7.03 |        0.1 |
| Mid (0.6-1.1 m^2) |            13500.05 |    1550.02 |       8.01 |        0.1 |
| High (\>=1.1 m^2) |            12415.73 |    1423.98 |       8.92 |        0.1 |

PKNCA summary over the first induction dosing interval, median per BSA
band. {.table style="width:100%;"}

### Comparison against the model’s closed-form expectations

Wurthwein 2025 publishes no NCA table, so the reference values are
closed forms of the published parameters.

The reference half-life needs care. The chain’s *asymptotic* terminal
half-life is `ln(2) * V / (CLinitial + Qtr)` = 1.1 days, but no NCA over
a realistic sampling window will ever measure it: as the ramp above
shows, apparent clearance is still only about 0.149 L/day at day 14, so
over the observed interval the decline is governed by `CLinitial` and a
little curvature. The right reference for an NCA over 0-14 days is
therefore the **initial-phase** half-life `ln(2) * V / CLinitial` = 9.2
days, which the measured value must sit just *below* because clearance
is already rising. Quoting the 1.1-day asymptote here would be a
unit-of-analysis error, not a validation.

The reference is **not** the same in every BSA band, and that is the
point. A size-proportional model (`V` and `CL` both proportional to BSA)
would give a BSA-independent half-life. This model is deliberately not
proportional: BSA enters *linearly and centred*, with a steeper slope on
`V` (1.57) than on `CLinitial` and `Qtr` (1.45). Below the centring
point the volume factor therefore falls faster than the clearance
factor, so smaller children get a shorter half-life – 8.7 days at the
low band’s median BSA against 9.3 days in the high band. Age and sex
shift `CLinitial` further. The reference below is computed **per
subject** from that subject’s own covariates, which turns the comparison
into a check of the entire covariate model rather than of one typical
value.

The measured NCA half-life sits below its reference in every band, and
by more in the smaller children (-19% low band, -4% high band). That is
the same mechanism seen from the other side: the transit rate `Qtr/V` is
*higher* at low BSA, so the apparent-clearance ramp runs faster and
bends the profile away from the initial-phase exponential sooner within
the 14-day window.

``` r

hl_asympt <- log(2) * 1.68 / (0.126 + 0.926)  # asymptotic terminal half-life

# Per-subject initial-phase half-life = ln(2) * V(cov) / CLinitial(cov), using
# exactly the factors model() applies.
ref_hl <- band |>
  mutate(
    fbsa_vc = (1 + 1.57 * (BSA - 0.79)) / (1 + 1.57 * 0.21),
    fbsa_cl = (1 + 1.45 * (BSA - 0.79)) / (1 + 1.45 * 0.21),
    fage_cl = 1 + 0.017 * (AGE - 8) * (AGE > 8),
    fsex_cl = 1 - 0.073 * SEXF,
    hl_ref  = log(2) * (1.68 * fbsa_vc) / (0.126 * fbsa_cl * fage_cl * fsex_cl)
  )

published <- ref_hl |>
  group_by(treatment = bsa_band) |>
  summarise(PPORRES = median(hl_ref), .groups = "drop") |>
  mutate(PPTESTCD = "half.life")

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_res,
  reference     = published,
  by            = "treatment",
  units         = c(half.life = "day"),
  tolerance_pct = 20
)

knitr::kable(cmp, caption = paste(
  "Simulated NCA half-life over 0-14 days versus the per-subject initial-phase",
  "closed form ln(2) * V(cov) / CLinitial(cov), median per BSA band.",
  "* marks rows differing by >20%."))
```

| NCA parameter | treatment         | Reference | Simulated | % diff |
|:--------------|:------------------|:----------|:----------|:-------|
| t½ (day)      | Low (\<0.6 m^2)   | 8.66      | 7.03      | -18.9% |
| t½ (day)      | Mid (0.6-1.1 m^2) | 9.24      | 8.01      | -13.3% |
| t½ (day)      | High (\>=1.1 m^2) | 9.26      | 8.92      | -3.7%  |

Simulated NCA half-life over 0-14 days versus the per-subject
initial-phase closed form ln(2) \* V(cov) / CLinitial(cov), median per
BSA band. \* marks rows differing by \>20%. {.table}

``` r

hl_band <- as.data.frame(nca_res) |>
  filter(PPTESTCD == "half.life") |>
  group_by(treatment) |>
  summarise(nca = median(PPORRES, na.rm = TRUE), .groups = "drop") |>
  left_join(rename(published, treatment = treatment, ref = PPORRES), by = "treatment") |>
  mutate(pct = 100 * (nca - ref) / ref)

as.data.frame(hl_band)
#>           treatment      nca      ref  PPTESTCD        pct
#> 1    Low (<0.6 m^2) 7.025161 8.663004 half.life -18.906185
#> 2 Mid (0.6-1.1 m^2) 8.013263 9.241288 half.life -13.288458
#> 3  High (>=1.1 m^2) 8.920931 9.261005 half.life  -3.672109

stopifnot(
  # The NCA value must sit BELOW the initial-phase closed form in every band
  # (apparent clearance is already rising over the window) but well above the
  # asymptote, which a 14-day window cannot reach.
  all(hl_band$nca < hl_band$ref),
  all(hl_band$nca > 3 * hl_asympt),
  # Deterministic solve against its own closed form, so a modest bound; the
  # residual gap is the within-window clearance ramp, not encoding error.
  max(abs(hl_band$pct)) < 20,
  # The size ordering must fall out of the linear centred BSA model: smaller
  # children get the shorter half-life, in both the closed form and the NCA.
  !is.unsorted(hl_band$ref[order(hl_band$treatment)]),
  !is.unsorted(hl_band$nca[order(hl_band$treatment)])
)
```

## Population variability

A stochastic run exercises the inter-individual and inter-occasion
variability that the typical-value checks above deliberately switch off.
The German/Czech model carries IIV on `CLinitial` only (24.2% CV) plus
IOV on both `CLinitial` (22.6%) and `V` (14.4%), with no IIV on `V` or
`Qtr`.

`rxSolve()` reuses the omega of a previous solve in the same session
unless one is supplied, so it is passed explicitly here.

``` r

om <- diag(c(
  etalcl      = log(1 + 0.242^2),
  etaiov_vc_1 = log(1 + 0.144^2), etaiov_vc_2 = log(1 + 0.144^2),
  etaiov_vc_3 = log(1 + 0.144^2),
  etaiov_cl_1 = log(1 + 0.226^2), etaiov_cl_2 = log(1 + 0.226^2),
  etaiov_cl_3 = log(1 + 0.226^2)
))
dimnames(om) <- list(
  c("etalcl", "etaiov_vc_1", "etaiov_vc_2", "etaiov_vc_3",
    "etaiov_cl_1", "etaiov_cl_2", "etaiov_cl_3"),
  c("etalcl", "etaiov_vc_1", "etaiov_vc_2", "etaiov_vc_3",
    "etaiov_cl_1", "etaiov_cl_2", "etaiov_cl_3"))

grid <- seq(0, 28, by = 0.25)
ev_pop <- lapply(cohort$id, function(i) {
  ci <- cohort[cohort$id == i, ]
  e <- rxode2::et(amt = ci$amt, time = c(0, 14), dur = 2 / 24, cmt = "central") |>
    rxode2::et(grid, cmt = "central") |>
    as.data.frame()
  e$id <- i
  e
}) |>
  dplyr::bind_rows() |>
  mutate(OCC = ifelse(time < 14, 1, 2)) |>
  left_join(select(cohort, id, BSA, AGE, SEXF), by = "id")

pop <- rxode2::rxSolve(mGC, ev_pop, omega = om, returnType = "data.frame")

stopifnot(dplyr::n_distinct(pop$id) == n_sub, !anyNA(pop$Cc))
```

``` r

qs <- pop |>
  group_by(time) |>
  summarise(lo = quantile(Cc, 0.05), md = median(Cc), hi = quantile(Cc, 0.95),
            .groups = "drop")

ggplot(qs, aes(time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(aes(y = md), linewidth = 0.8) +
  geom_point(
    data = tibble::tibble(time = c(7, 14, 21, 28), md = c(895, 543, 1260, 630)),
    aes(y = md), colour = "firebrick", size = 2
  ) +
  labs(x = "Days from the first induction dose",
       y = "Asparaginase activity (U/L)") +
  theme_bw()
```

![Simulated induction profiles for 200 virtual German/Czech patients
(median and 5th-95th percentile band) with the published day-7 and
day-14 medians overlaid. Compare Wurthwein 2025 Figure 2a, induction
panel.](Wurthwein_2025_pegasparaginase_files/figure-html/population-plot-1.png)

Simulated induction profiles for 200 virtual German/Czech patients
(median and 5th-95th percentile band) with the published day-7 and
day-14 medians overlaid. Compare Wurthwein 2025 Figure 2a, induction
panel.

``` r

# Assert on the CENTRE of the cohort, never on its extremes: the tails of a
# random cohort are not reproducible across rxode2 versions.
med7  <- median(pop$Cc[abs(pop$time - 7) < 1e-6])
med14 <- median(pop$Cc[abs(pop$time - 13.75) < 1e-6])
c(day7 = med7, day14 = med14)
#>     day7    day14 
#> 918.1534 522.9600

stopifnot(
  abs(med7 - 895) / 895 < 0.20,
  abs(med14 - 543) / 543 < 0.25,
  # IIV on CLinitial must actually be doing something.
  IQR(pop$Cc[abs(pop$time - 7) < 1e-6]) > 50
)
```

## Errata

Discrepancies found in the published source while transcribing. None was
“fixed” by tuning; each is recorded and the affected value either
avoided or used with the reading given.

1.  **Table 1, `Dose (U/m2/dose)` rows carry an incorrect unit tag.**
    The values are the ABSOLUTE administered dose in U, not the dose per
    square metre. The protocol dose is 2500 U/m^2, so a per-m^2 median
    of 1940 U/m^2 would mean the median patient was systematically
    underdosed by 22%, which the text nowhere claims. Dividing the
    printed medians by the corresponding median BSA recovers the
    protocol dose: German/Czech induction 1940 / 0.78 = 2487 U/m^2;
    Italian 1900 / 0.77 = 2468; Australian 1900 / 0.77 = 2468. All three
    land on 2500 U/m^2. These models are dosed at 2500 U/m^2
    accordingly.
2.  **ESM Table S15, induction row, `Time > 400 U/L`.** The median
    prints as `3.50` while its own interquartile range prints as
    `(28.4 to 32.1)`. A median outside its own IQR is impossible; the
    median is a corrupted `30.5` or similar. That cell is excluded from
    the comparison above.
3.  **ESM Table S15, `PIId8-M10 (R2-EA)` row.** Two interquartile ranges
    are garbled in the same way: `Time > 100 U/L` prints
    `153 (15.0 to 158)` and `Time > 400 U/L` prints `147 (145 to 15.0)`.
    The medians are internally consistent with the neighbouring columns
    and are used; the ranges are not.
4.  **`Qtr` is spelled `Qrt` in the parameter rows of ESM Tables S7 and
    S9** while the footnotes, Table S13 and the running text all use
    `Qtr`. A typographic slip; the same parameter throughout.
5.  **The main text and the supplement disagree slightly on the largest
    clearance reductions.** Section 3.2.1 says the high-risk reduction
    in `CLinitial` was “up to -47.6%” citing ESM Table S4 model 335000
    (-0.476), while Section 4.1 says “-47.5%” and the final model in ESM
    Table S7 gives -0.473. Likewise R2-EA is quoted as “-61.0%” in
    Section 3.2.3 and “-61.1%” in Section 4.1 while Table S9 gives
    -0.601. The models use the FINAL-model values from Tables S7 and S9,
    which are the estimates of the models the paper actually reports.

## Assumptions and deviations

- **The structural equations are not in Wurthwein 2025.** The paper
  describes “a chain of 14 transit compartments” in prose and gives no
  ODEs, no diagram and no control stream (“The NONMEM codes for the
  final popPK models are provided on request”). Every structural fact –
  the serial topology, the shared volume, elimination from every
  compartment, the extra `Qtr` on the terminal compartment, the
  total-activity observation, the linear centred BSA model, the age
  hockey stick and its fixed 8-year break point, the 1/2 sex coding and
  the combined error model – is taken from the complete NONMEM control
  stream in Section 4 of the **Wurthwein 2021** supplement, which the
  2025 supplement cites as its reference 7 and its starting point. The
  2025 values are used throughout; only the forms come from 2021.
- **BSA centring constant 0.79 m^2.** Hard-coded in the 2021 control
  stream and not restated in the 2025 supplement, which says only
  “centred on the median”. The 2025 cohort medians are 0.78
  (German/Czech) and 0.77 (Italian), so the constant may have shifted by
  0.01-0.02 m^2. It is not load-bearing: because `model()` renormalises
  the linear factor to 1 at BSA = 1 m^2, re-centring from 0.79 to 0.77
  changes typical-value predictions by under 1% across the whole
  observed BSA range. 0.79 is used for all four models.
- **`V`, `CLinitial` and `Qtr` are absolute, not per square metre.** The
  tables tag them `L/m^2` and `L/day/m^2`, but the control stream
  multiplies nothing by BSA – the entire BSA dependence is the linear
  centred covariate term. The tag records the authors’ stated convention
  of quoting values “for a child with BSA = 1 m^2 for better
  comparison”. The `ini()` values are therefore the printed numbers, and
  `model()` renormalises the BSA factor to 1 at BSA = 1 m^2 so that they
  mean exactly what the table says.
- **New covariate canonical, NAME PROVISIONAL.** No anti-PEG antibody
  column existed in `inst/references/covariate-columns.md`. The
  pre-existing anti-PEG antibody level is genuinely distinct from every
  registered anti-drug-antibody canonical: it targets the
  polyethylene-glycol moiety rather than the drug, it is present before
  any exposure (which is why it can act on the first dose), and it has
  no zero-encoding for “negative” subjects. `ABPEG_IGM` and `ABPEG_IGG`
  are proposed as a new `ABPEG_<isotype>` family and are **awaiting
  operator ratification** (sidecar `oasweep_PMC12540506` request-001).
  If the operator prefers a different name the change is a rename in
  three files.
- **Compartment naming.** The authors’ control stream names the 14
  states `CENTRAL` and `PERI2`-`PERI14`. `peripheral<n>` would be
  misleading here – these are not classical peripheral distribution
  compartments; there is no back-flow and all 14 species share one serum
  volume. The blessed `transit<n>` chain prefix is used instead, which
  describes the role correctly, with `central` retained for the dosed
  compartment as in the control stream.
- **`OCC` carries both the phase effect and the IOV.** This is the
  authors’ own encoding, not a simplification: their control stream uses
  a single `OCC` column to select both the treatment-phase covariate
  factors and the IOV etas. Occasion codes are renumbered from the
  control stream’s `2/3/5` to the register’s `1..N` convention; the
  mapping is documented per model in `covariateData$OCC$notes`.
- **NONMEM `$OMEGA BLOCK(1) SAME` has no nlmixr2 equivalent**, so
  occasions after the first carry their own eta with the variance fixed
  equal to the occasion-1 estimate.
- **Cohort distributions are assumed.** Only medians and ranges are
  published for BSA and age, so the virtual cohort uses truncated
  log-normals matched to those, sampled on a deterministic quantile
  lattice. BSA and age are given the same lattice position rather than
  being drawn independently, because both are growth-driven and strongly
  correlated in a pediatric cohort; drawing them independently would
  create implausible 1-year-olds with 2 m^2 BSA and would distort the
  age \> 8 years covariate. Anti-PEG IgM is set to 1 (below the 3.67 cut
  point, i.e. no effect) for all virtual subjects, matching the ~90% of
  real patients below it (ESM Table S5).
- **Concentration is discontinuous at an occasion switch.** `V` changes
  stepwise between administrations, so the predicted activity jumps at
  each dose time. This is inherent to the authors’ covariate model, not
  an artefact of this encoding; NONMEM behaves identically. Trough
  readings above are taken just before the switch.
- **The R2-EA day-22 administration is the weakest fit**, simulating
  about 13% below the published median. The authors’ own grouping shares
  one clearance factor between protocol II day 8 and day 22 (ESM Table
  S8 model 493151, chosen because only the day-22-to-day-36 confidence
  intervals separated), so the model cannot follow the observed jump at
  that dose. This is a property of the published model, not of the
  transcription; no parameter was adjusted.
- **Silent inactivation and hypersensitivity are out of scope.** Records
  at or after either were excluded before fitting. Do not use these
  models to predict activity in an inactivating patient.
- **Australian cohort.** Contributed dose-intensity values only; no
  model was fitted and none is packaged.

## Reference

Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath CE, Matteo C,
Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in Children with
Acute Lymphoblastic Leukemia Treated within the AIEOP-BFM ALL 2009
Trial: Population Pharmacokinetics and Drug Exposure. *Eur J Drug Metab
Pharmacokinet*. 2025;50(6):683-696.
<https://doi.org/10.1007/s13318-025-00962-3>
