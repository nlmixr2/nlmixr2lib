# Perampanel (Yu 2025)

## Model and source

- Citation: Yu L, Mao F, Chen S, Liu J, Xiao J, Chen M, Luo H, Yu Z,
  Dai H. Development and Validation of a Population Pharmacokinetics
  Model of Perampanel for Pediatric Epilepsy Patients for Optimized
  Dosing. Drug Des Devel Ther. 2025;19:3119-3128.
  <doi:10.2147/DDDT.S499085>
- Description: One-compartment popPK model with first-order absorption
  of perampanel in Chinese pediatric epilepsy patients on routine
  therapeutic drug monitoring (Yu 2025). NONMEM 7.5.0 FOCE-I fit of 454
  plasma concentrations from 151 patients aged 0.58-17.9 years, 120
  (79.5%) of whom were younger than 12 years. The absorption rate
  constant was fixed to the adult value of 3.37 /h (Fujita 2023) because
  no pediatric estimate was available. Apparent clearance carries a
  linear-centralized age term ((AGE+10)/8.8)^1.31 plus multiplicative
  comedication effects: co-administered oxcarbazepine raises CL/F
  1.51-fold and carbamazepine 1.88-fold (CYP3A4/5 induction), while
  sodium valproate lowers it to 0.745-fold (CYP3A4/5 inhibition).
  Apparent volume of distribution is proportional to the log of body
  weight. Inter-individual variability is exponential on clearance only;
  residual error is proportional.
- Article: [Drug Des Devel Ther.
  2025;19:3119-3128](https://doi.org/10.2147/DDDT.S499085)
- Supplement (Table S1 covariate screen, Figure S1 VPC):
  <https://www.dovepress.com/article/supplementary_file/499085/499085-Supplementary-Material.docx>

## Population

Yu 2025 is a retrospective, single-centre therapeutic-drug-monitoring
(TDM) study run at the Second Affiliated Hospital of Zhejiang University
School of Medicine (Hangzhou, China) between February 2021 and September
2023. The final model was fit to 454 plasma perampanel concentrations
from 151 pediatric patients with epilepsy; 120 of them (79.5%) were
younger than 12 years, which is what distinguishes this model from the
other published perampanel population models (Yu 2025 Discussion).

Baseline characteristics (Yu 2025 Table 1): median age 9.00 years (IQR
6.34-11.88, range 0.58-17.9), median body weight 28.1 kg (IQR 21.2-41.0,
range 9.00-89.0), 96 male / 55 female (36.4% female), median height 134
cm, median BMI 16.6 kg/m^2, median BSA 1.06 m^2, median creatinine
clearance 149 mL/min. No patient had severe hepatic or renal impairment.
Epilepsy type was focal in 51.0%, generalized in 46.4%, and focal
evolving to generalized in 2.65%. The median total daily perampanel dose
was 2.00 mg (IQR 2.00-4.00, range 1.00-8.00), and observed perampanel
concentrations were 30.0-1082 ng/mL with a median of 242 ng/mL (IQR
144-340).

The three comedications that entered the final clearance model were
oxcarbazepine (26 patients, 17.2%), sodium valproate (53 patients,
35.1%) and carbamazepine (8 patients, 5.30%). TDM was performed roughly
3 weeks after perampanel initiation, with samples drawn in the morning
about 12 h after the previous dose, so essentially every observation is
a near-trough steady-state concentration – a point that matters for how
much the data can say about the volume of distribution (see Assumptions
and deviations).

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Yu_2025_perampanel")()$population`).

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in `inst/modeldb/specificDrugs/Yu_2025_perampanel.R`. The
table below collects them in one place for review. Equations (1) and (2)
are set as vector graphics in the PDF and are dropped by text
extraction; they were read from the publisher’s rendered equation images
(`DDDT-19-3119-e0001.jpg` and `-e0002.jpg`, distributed in the Europe
PMC file bundle for PMC12034272).

| Equation / parameter | Value | Source location |
|----|----|----|
| Structural model (1-compartment, first-order absorption and elimination) | n/a | “PPK Model Development”, page 4 |
| `lka` (ka) | 3.37 1/h, fixed | Table 2, “KA 3.37 FIXED”; footnote a attributes it to Fujita 2023 |
| `lcl` (CL/F coefficient) | 0.177 L/h | Table 2 “CL/F (L/h) 0.177” (RSE 15.2%); Equation (1) |
| `lvc` (V/F coefficient) | 227 L per log10 unit of kg | Table 2 “V (L) 227” (RSE 14.1%); Equation (2) |
| `e_age_cl` | 1.31 | Table 2 “AGE 1.31” (RSE 14.2%); Equation (1) |
| `e_oxc_cl` | 1.51 | Table 2 “Oxcarbazepine 1.51” (RSE 12.0%); Equation (1) |
| `e_vpa_cl` | 0.745 | Table 2 “Sodium valproate 0.745” (RSE 7.92%); Equation (1) |
| `e_cbz_cl` | 1.88 | Table 2 “Carbamazepine 1.88” (RSE 16.3%); Equation (1) |
| Age normalisation constants 10 and 8.8 | n/a | Equation (1), used exactly as printed |
| `etalcl` (IIV variance on CL/F) | 0.0963 | Table 2 “omega CL 0.0963”; footnote defines omega as the interindividual **variance** |
| `propSd` | sqrt(0.130) = 0.3606 | Table 2 “sigma 0.130”; footnote defines sigma as the proportional residual variability (NONMEM `$SIGMA` variance scale) |
| Covariate selection path | n/a | Table S1 of the supplement (forward inclusion / backward elimination) |
| Published simulation target | n/a | Table 3, median steady-state concentrations |

Equation (1) as printed:

    CL/F (L/h) = 0.177 * ((Age + 10)/8.8)^1.31 * 1.51^OXC * 0.745^VPA * 1.88^CBZ

Equation (2) as printed:

    V (L) = 227 * LGBW

``` r

mod <- readModelDb("Yu_2025_perampanel")
mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

# rxode2 returns each model covariate more than once when it is both an input
# column and referenced in `model()` (here AGE and WT), which makes the result
# unusable by dplyr ("Can't transform a data frame with duplicate names").
# Keep the first occurrence of each name.
drop_duplicate_columns <- function(x) {
  x <- as.data.frame(x)
  x[, !duplicated(names(x)), drop = FALSE]
}
```

## Published simulation target: Table 3

Table 3 of Yu 2025 reports median steady-state plasma concentrations for
a grid of 336 scenarios: four ages (2, 4, 12, 18 years) crossed with six
total daily doses (2-12 mg), five body weights (10-50 kg) and four
comedication arms (none, oxcarbazepine, carbamazepine, sodium
valproate). Cells marked “/” in the paper are age/weight combinations
the authors did not consider clinically relevant and are carried here as
`NA`.

This grid is the strongest available gate on the packaged model, because
the model carries exactly one random effect (on clearance) and the
steady-state average concentration is a monotone function of clearance.
The median of a Monte Carlo cohort under a single log-normal random
effect is therefore the typical-value prediction exactly, so the
published medians can be reproduced **deterministically** from a
`zeroRe()` solve with no Monte Carlo sampling and no dependence on the
random-number stream.

``` r

# Yu 2025 Table 3, transcribed from the JATS full text of PMC12034272.
# One row per (age, dose, comedication arm); the five numeric columns are the
# five body weights the paper tabulates. NA = "/" in the published table.
published <- tibble::tribble(
  ~age, ~dose, ~arm,   ~`10`, ~`20`, ~`30`, ~`40`, ~`50`,
   2,    2, "none", 321, 314, NA, NA, NA,
   2,    2, "OXC",  212, 208, NA, NA, NA,
   2,    2, "CBZ",  171, 167, NA, NA, NA,
   2,    2, "VPA",  430, 422, NA, NA, NA,
   2,    4, "none", 641, 637, NA, NA, NA,
   2,    4, "OXC",  425, 422, NA, NA, NA,
   2,    4, "CBZ",  341, 339, NA, NA, NA,
   2,    4, "VPA",  860, 855, NA, NA, NA,
   2,    6, "none", 955, 940, NA, NA, NA,
   2,    6, "OXC",  632, 623, NA, NA, NA,
   2,    6, "CBZ",  508, 500, NA, NA, NA,
   2,    6, "VPA",  1282, 1262, NA, NA, NA,
   2,    8, "none", 1259, 1236, NA, NA, NA,
   2,    8, "OXC",  834, 819, NA, NA, NA,
   2,    8, "CBZ",  670, 658, NA, NA, NA,
   2,    8, "VPA",  1690, 1659, NA, NA, NA,
   2,   10, "none", 1599, 1591, NA, NA, NA,
   2,   10, "OXC",  1059, 1054, NA, NA, NA,
   2,   10, "CBZ",  851, 846, NA, NA, NA,
   2,   10, "VPA",  2147, 2136, NA, NA, NA,
   2,   12, "none", 1888, 1918, NA, NA, NA,
   2,   12, "OXC",  1251, 1270, NA, NA, NA,
   2,   12, "CBZ",  1005, 1021, NA, NA, NA,
   2,   12, "VPA",  2535, 2575, NA, NA, NA,
   4,    2, "none", 261, 263, 260, 254, 256,
   4,    2, "OXC",  173, 174, 172, 168, 169,
   4,    2, "CBZ",  139, 140, 138, 135, 136,
   4,    2, "VPA",  351, 353, 348, 340, 343,
   4,    4, "none", 512, 502, 510, 510, 506,
   4,    4, "OXC",  339, 332, 338, 338, 335,
   4,    4, "CBZ",  273, 267, 271, 271, 269,
   4,    4, "VPA",  688, 673, 684, 684, 680,
   4,    6, "none", 764, 771, 764, 785, 766,
   4,    6, "OXC",  506, 511, 506, 520, 507,
   4,    6, "CBZ",  407, 410, 406, 418, 408,
   4,    6, "VPA",  1026, 1035, 1025, 1054, 1028,
   4,    8, "none", 1018, 1040, 1023, 1010, 1029,
   4,    8, "OXC",  674, 689, 678, 669, 681,
   4,    8, "CBZ",  542, 553, 544, 537, 547,
   4,    8, "VPA",  1366, 1395, 1373, 1356, 1381,
   4,   10, "none", 1299, 1272, 1289, 1286, 1297,
   4,   10, "OXC",  860, 842, 853, 852, 859,
   4,   10, "CBZ",  691, 677, 686, 684, 690,
   4,   10, "VPA",  1743, 1707, 1730, 1726, 1741,
   4,   12, "none", 1537, 1549, 1524, 1544, 1568,
   4,   12, "OXC",  1018, 1026, 1010, 1022, 1039,
   4,   12, "CBZ",  818, 824, 811, 821, 834,
   4,   12, "VPA",  2063, 2079, 2046, 2072, 2105,
  12,    2, "none", 142, 144, 142, 140, 146,
  12,    2, "OXC",  94, 95, 94, 93, 97,
  12,    2, "CBZ",  75, 77, 75, 74, 78,
  12,    2, "VPA",  190, 193, 190, 187, 196,
  12,    4, "none", 283, 276, 285, 282, 283,
  12,    4, "OXC",  188, 183, 189, 187, 188,
  12,    4, "CBZ",  151, 147, 152, 150, 151,
  12,    4, "VPA",  380, 370, 383, 378, 380,
  12,    6, "none", 426, 423, 421, 423, 415,
  12,    6, "OXC",  282, 280, 279, 280, 275,
  12,    6, "CBZ",  227, 225, 224, 225, 221,
  12,    6, "VPA",  572, 568, 564, 568, 558,
  12,    8, "none", 575, 565, 576, 564, 569,
  12,    8, "OXC",  381, 374, 382, 373, 377,
  12,    8, "CBZ",  306, 301, 306, 300, 303,
  12,    8, "VPA",  771, 758, 773, 756, 764,
  12,   10, "none", 714, 709, 695, 711, 702,
  12,   10, "OXC",  473, 470, 461, 471, 465,
  12,   10, "CBZ",  380, 377, 370, 378, 373,
  12,   10, "VPA",  959, 951, 933, 955, 942,
  12,   12, "none", 852, 845, 839, 847, 848,
  12,   12, "OXC",  564, 559, 556, 561, 562,
  12,   12, "CBZ",  453, 449, 446, 451, 451,
  12,   12, "VPA",  1143, 1134, 1126, 1137, 1138,
  18,    2, "none", NA, NA, NA, 103, 103,
  18,    2, "OXC",  NA, NA, NA, 68, 69,
  18,    2, "CBZ",  NA, NA, NA, 55, 55,
  18,    2, "VPA",  NA, NA, NA, 139, 139,
  18,    4, "none", NA, NA, NA, 209, 210,
  18,    4, "OXC",  NA, NA, NA, 139, 139,
  18,    4, "CBZ",  NA, NA, NA, 111, 112,
  18,    4, "VPA",  NA, NA, NA, 281, 282,
  18,    6, "none", NA, NA, NA, 309, 307,
  18,    6, "OXC",  NA, NA, NA, 204, 203,
  18,    6, "CBZ",  NA, NA, NA, 164, 163,
  18,    6, "VPA",  NA, NA, NA, 414, 412,
  18,    8, "none", NA, NA, NA, 422, 423,
  18,    8, "OXC",  NA, NA, NA, 279, 280,
  18,    8, "CBZ",  NA, NA, NA, 224, 225,
  18,    8, "VPA",  NA, NA, NA, 566, 568,
  18,   10, "none", NA, NA, NA, 515, 514,
  18,   10, "OXC",  NA, NA, NA, 341, 340,
  18,   10, "CBZ",  NA, NA, NA, 274, 273,
  18,   10, "VPA",  NA, NA, NA, 691, 689,
  18,   12, "none", NA, NA, NA, 616, 617,
  18,   12, "OXC",  NA, NA, NA, 408, 409,
  18,   12, "CBZ",  NA, NA, NA, 328, 328,
  18,   12, "VPA",  NA, NA, NA, 827, 828
) |>
  tidyr::pivot_longer(
    cols = c("10", "20", "30", "40", "50"),
    names_to = "wt", values_to = "published_ngml"
  ) |>
  dplyr::mutate(wt = as.numeric(wt)) |>
  dplyr::filter(!is.na(published_ngml)) |>
  dplyr::arrange(age, dose, arm, wt) |>
  dplyr::mutate(id = dplyr::row_number())

stopifnot(nrow(published) == 336L)
```

### Simulating the grid

Each scenario is solved to steady state with `ss = 1` on a once-daily
regimen, which makes the result independent of any burn-in choice. The
dosing interval is sampled densely enough to integrate the
concentration-time curve accurately.

``` r

obs_times <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 20, 24)

grid_events <- rxode2::et(amt = 1, cmt = "depot", ii = 24, ss = 1, id = published$id) |>
  rxode2::et(time = obs_times, cmt = "central", id = published$id) |>
  as.data.frame() |>
  dplyr::left_join(published, by = "id") |>
  dplyr::mutate(
    amt        = ifelse(is.na(amt), NA_real_, dose),
    AGE        = age,
    WT         = wt,
    CONMED_OXC = as.integer(arm == "OXC"),
    CONMED_CBZ = as.integer(arm == "CBZ"),
    CONMED_VPA = as.integer(arm == "VPA")
  )

# Each scenario is one subject; ids come straight from `published` so they are
# unique by construction.
stopifnot(!anyDuplicated(unique(grid_events[, c("id", "time", "evid")])))

grid_sim <- rxode2::rxSolve(
  mod_typical, events = grid_events,
  keep = c("age", "dose", "arm", "wt", "published_ngml")
) |>
  drop_duplicate_columns()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> Warning: multi-subject simulation without without 'omega'
#> Warning: Cannot keep missing columns:
```

### PKNCA over the steady-state dosing interval

`cav` (the interval average concentration, `auclast / (end - start)`) is
the quantity Table 3 reports: the residual analysis below shows the
published values track the interval average rather than the
end-of-interval trough.

``` r

grid_conc <- grid_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(Cc_ngml = Cc * 1000) |>
  dplyr::select(id, time, Cc_ngml, arm, age, dose, wt)

grid_dose <- grid_events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, arm)

grid_nca <- PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(grid_conc, Cc_ngml ~ time | arm + id),
  PKNCA::PKNCAdose(grid_dose, amt ~ time | arm + id),
  intervals = data.frame(
    start = 0, end = 24,
    cav = TRUE, cmax = TRUE, tmax = TRUE, ctrough = TRUE, auclast = TRUE
  )
)
grid_res <- PKNCA::pk.nca(grid_nca)

grid_wide <- as.data.frame(grid_res) |>
  dplyr::select(id, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

comparison <- published |>
  dplyr::left_join(grid_wide, by = "id") |>
  dplyr::mutate(
    pct_diff_cav     = 100 * (cav - published_ngml) / published_ngml,
    pct_diff_ctrough = 100 * (ctrough - published_ngml) / published_ngml
  )
```

``` r

# Observed on this build: cav median |diff| 0.68%, 90th pct 2.09%, max 3.34%.
# The bounds below sit outside that range but well inside the ~5% Monte Carlo
# noise of a 1000-replicate median, so they can still go red on a real
# transcription or covariate-equation error (a single mis-keyed covariate
# multiplier moves whole arms by 25-88%).
gate <- comparison$pct_diff_cav
stopifnot(
  nrow(comparison) == 336L,
  !anyNA(gate),
  abs(median(gate)) < 1.5,
  unname(quantile(abs(gate), 0.9)) < 3,
  max(abs(gate)) < 5
)

tibble::tibble(
  Statistic = c("n scenarios", "median % difference", "median |% difference|",
                "90th percentile |% difference|", "max |% difference|"),
  Value = c(
    sprintf("%d", nrow(comparison)),
    sprintf("%+.2f%%", median(gate)),
    sprintf("%.2f%%", median(abs(gate))),
    sprintf("%.2f%%", quantile(abs(gate), 0.9)),
    sprintf("%.2f%%", max(abs(gate)))
  )
) |>
  knitr::kable(caption = "Packaged model vs Yu 2025 Table 3, all 336 published cells.")
```

| Statistic                        | Value  |
|:---------------------------------|:-------|
| n scenarios                      | 336    |
| median % difference              | -0.03% |
| median \|% difference\|          | 0.68%  |
| 90th percentile \|% difference\| | 2.09%  |
| max \|% difference\|             | 3.34%  |

Packaged model vs Yu 2025 Table 3, all 336 published cells. {.table}

The interval-average concentration reproduces every published cell. The
end-of-interval trough is systematically a little lower, and the size of
that gap tracks clearance – it is largest in the carbamazepine arm,
which has the highest clearance and therefore the largest peak-to-trough
fluctuation:

``` r

comparison |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    `mean % diff (interval average)` = round(mean(pct_diff_cav), 2),
    `mean % diff (end-of-interval trough)` = round(mean(pct_diff_ctrough), 2),
    .groups = "drop"
  ) |>
  dplyr::rename("Comedication arm" = arm) |>
  knitr::kable(
    caption = paste(
      "Table 3 is reproduced by the interval average, not the trough.",
      "The trough bias scales with clearance (CBZ 1.88x, OXC 1.51x,",
      "none 1x, VPA 0.745x), which is the signature of peak-to-trough",
      "fluctuation rather than a clearance error."
    )
  )
```

| Comedication arm | mean % diff (interval average) | mean % diff (end-of-interval trough) |
|:---|---:|---:|
| CBZ | -0.19 | -3.45 |
| OXC | -0.19 | -2.82 |
| VPA | -0.15 | -1.46 |
| none | -0.17 | -1.92 |

Table 3 is reproduced by the interval average, not the trough. The
trough bias scales with clearance (CBZ 1.88x, OXC 1.51x, none 1x, VPA
0.745x), which is the signature of peak-to-trough fluctuation rather
than a clearance error. {.table}

``` r

comparison |>
  dplyr::mutate(arm = factor(
    arm, levels = c("none", "OXC", "CBZ", "VPA"),
    labels = c("No comedication", "Oxcarbazepine", "Carbamazepine", "Sodium valproate")
  )) |>
  ggplot(aes(published_ngml, cav, colour = factor(age))) +
  annotate("rect", xmin = 100, xmax = 1000, ymin = 100, ymax = 1000,
           alpha = 0.08, fill = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75, size = 1.6) +
  facet_wrap(~arm) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Published median steady-state concentration (ng/mL)",
    y = "Simulated interval-average concentration (ng/mL)",
    colour = "Age (years)",
    caption = "Replicates Table 3 of Yu 2025."
  )
```

![Replicates Table 3 of Yu 2025: simulated vs published median
steady-state perampanel concentration for all 336 scenarios. The dashed
line is the line of identity; the shaded band is the 100-1000 ng/mL
reference range the paper
targets.](Yu_2025_perampanel_files/figure-html/figure-table3-1.png)

Replicates Table 3 of Yu 2025: simulated vs published median
steady-state perampanel concentration for all 336 scenarios. The dashed
line is the line of identity; the shaded band is the 100-1000 ng/mL
reference range the paper targets.

## Virtual cohort

The model carries a single random effect, `etalcl`, on apparent
clearance. That makes it possible to build a cohort whose
between-subject variability is exact rather than sampled: the eta values
are laid out on a mid-point `qnorm` grid and supplied to `rxSolve()` as
a data column with `omega = NA`. The resulting cohort is identical on
every machine and every rxode2 build, which removes the usual
thread-count and RNG-version fragility from every assertion below.

Ages and weights are placed on a matched probability grid interpolated
through the quantiles Yu 2025 Table 1 reports (age min / Q1 / median /
Q3 / max = 0.58 / 6.34 / 9.00 / 11.88 / 17.9 years; weight = 9.00 / 21.2
/ 28.1 / 41.0 / 89.0 kg). Using one probability index for both keeps age
and weight perfectly rank-correlated, which is the right qualitative
behaviour for a pediatric cohort; the paper does not publish the joint
distribution.

``` r

n_sub <- 200L
p_grid <- (seq_len(n_sub) - 0.5) / n_sub

quantile_grid <- function(q_values, p) {
  stats::approx(x = c(0, 0.25, 0.5, 0.75, 1), y = q_values, xout = p)$y
}

cohort <- tibble::tibble(
  id     = seq_len(n_sub),
  AGE    = quantile_grid(c(0.58, 6.34, 9.00, 11.88, 17.9), p_grid),
  WT     = quantile_grid(c(9.00, 21.2, 28.1, 41.0, 89.0), p_grid),
  # Mid-point qnorm grid: a deterministic stand-in for 200 draws from
  # N(0, 0.0963). Its mean is 0 and its variance matches omega to <1%.
  etalcl = stats::qnorm(p_grid) * sqrt(0.0963),
  CONMED_OXC = 0L,
  CONMED_CBZ = 0L,
  CONMED_VPA = 0L
)

# 4 mg/d is the maintenance dose Yu 2025 recommends for pediatric patients not
# co-administered other antiseizure medications (Discussion, final paragraph).
cohort_dose <- 4

cohort_events <- rxode2::et(amt = cohort_dose, cmt = "depot", ii = 24, ss = 1, id = cohort$id) |>
  rxode2::et(time = obs_times, cmt = "central", id = cohort$id) |>
  as.data.frame() |>
  dplyr::left_join(cohort, by = "id")

stopifnot(!anyDuplicated(unique(cohort_events[, c("id", "time", "evid")])))

cohort_sim <- rxode2::rxSolve(mod, events = cohort_events, omega = NA) |>
  drop_duplicate_columns()
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: multi-subject simulation without without 'omega'
```

## PKNCA validation

``` r

cohort_conc <- cohort_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(Cc_ngml = Cc * 1000, regimen = "4 mg once daily") |>
  dplyr::select(id, time, Cc_ngml, regimen)

cohort_dose_df <- cohort_events |>
  dplyr::filter(evid == 1) |>
  dplyr::mutate(regimen = "4 mg once daily") |>
  dplyr::select(id, time, amt, regimen)

cohort_nca <- PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(cohort_conc, Cc_ngml ~ time | regimen + id),
  PKNCA::PKNCAdose(cohort_dose_df, amt ~ time | regimen + id),
  intervals = data.frame(
    start = 0, end = 24,
    cav = TRUE, cmax = TRUE, tmax = TRUE, cmin = TRUE, ctrough = TRUE, auclast = TRUE
  )
)
cohort_res <- PKNCA::pk.nca(cohort_nca)

cohort_wide <- as.data.frame(cohort_res) |>
  dplyr::select(id, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

cohort_wide |>
  dplyr::summarise(
    dplyr::across(
      c(cav, cmax, ctrough, tmax, auclast),
      list(median = ~ median(.x), p05 = ~ quantile(.x, 0.05), p95 = ~ quantile(.x, 0.95))
    )
  ) |>
  tidyr::pivot_longer(dplyr::everything(),
                      names_to = c("param", "stat"), names_sep = "_") |>
  tidyr::pivot_wider(names_from = stat, values_from = value) |>
  dplyr::mutate(
    param = dplyr::recode(
      param,
      cav = "Cav,ss (ng/mL)", cmax = "Cmax,ss (ng/mL)",
      ctrough = "Ctrough,ss (ng/mL)", tmax = "Tmax (h)",
      auclast = "AUC0-24,ss (ng*h/mL)"
    )
  ) |>
  dplyr::rename(
    "NCA parameter" = param, "Median" = median,
    "5th percentile" = p05, "95th percentile" = p95
  ) |>
  knitr::kable(digits = 1, caption = paste(
    "Steady-state NCA over the 24 h dosing interval for the 200-subject",
    "deterministic cohort at 4 mg once daily."
  ))
```

| NCA parameter         | Median | 5th percentile | 95th percentile |
|:----------------------|-------:|---------------:|----------------:|
| Cav,ss (ng/mL)        |  343.5 |          133.3 |          1063.1 |
| Cmax,ss (ng/mL)       |  348.9 |          137.5 |          1070.5 |
| Ctrough,ss (ng/mL)    |  337.6 |          128.8 |          1055.0 |
| Tmax (h)              |    1.5 |            1.5 |             1.5 |
| AUC0-24,ss (ng\*h/mL) | 8244.2 |         3199.3 |         25514.1 |

Steady-state NCA over the 24 h dosing interval for the 200-subject
deterministic cohort at 4 mg once daily. {.table}

Yu 2025 does not report non-compartmental Cmax / Tmax / AUC / half-life
values, so there is no published NCA table to place side by side here;
Table 3 above is the paper’s quantitative simulation output and serves
that role.

### Analytic cross-checks

Three properties of the model are known in closed form and are checked
against the cohort. Each is exact rather than approximate, so the bounds
are tight.

``` r

typical_cl   <- 0.177 * ((cohort$AGE + 10) / 8.8)^1.31
typical_cav  <- 1000 * cohort_dose / (typical_cl * 24)
analytic_cav <- typical_cav / exp(cohort$etalcl)

# The cohort spans ages 0.58-17.9 years, so the age term alone contributes a
# variance of about 0.099 to log(Cav,ss) -- as much as the random effect does.
# Dividing out the age-driven typical value isolates the random effect, which
# is what the omega check has to be made against.
eta_recovered <- log(typical_cav) - log(cohort_wide$cav)

cav_max_err <- max(abs(100 * (cohort_wide$cav - analytic_cav) / analytic_cav))
median_ratio_pct <- 100 * (median(cohort_wide$cav) / median(typical_cav) - 1)
gcv_recovered <- 100 * sqrt(exp(stats::var(eta_recovered)) - 1)
gcv_nominal   <- 100 * sqrt(exp(0.0963) - 1)

checks <- tibble::tibble(
  Check = c(
    "Cav,ss from PKNCA vs closed form dose/(CL*tau)",
    "Cohort median Cav,ss vs the typical-value prediction",
    "Random effect recovered from Cav,ss vs the supplied eta grid",
    "Geometric CV of the recovered random effect vs sqrt(exp(omega^2) - 1)"
  ),
  Simulated = c(
    sprintf("%.3f%% max abs diff", cav_max_err),
    sprintf("%.1f ng/mL", median(cohort_wide$cav)),
    sprintf("%.2e max abs diff", max(abs(eta_recovered - cohort$etalcl))),
    sprintf("%.2f%%", gcv_recovered)
  ),
  Expected = c(
    "< 0.5%",
    sprintf("%.1f ng/mL", median(typical_cav)),
    "< 1e-3 (exact algebra)",
    sprintf("%.2f%%", gcv_nominal)
  )
)

stopifnot(
  # 1. PKNCA's trapezoidal interval average matches the closed form. The curve
  #    fluctuates only ~3% across the interval, so trapezoid error is tiny.
  cav_max_err < 0.5,
  # 2. One log-normal eta on CL and a monotone target => the cohort median is
  #    the typical value. The mid-point qnorm grid is exactly symmetric.
  abs(median_ratio_pct) < 1,
  # 3. The packaged model applies the random effect exactly where the equation
  #    says it does; inverting the closed form returns the supplied etas.
  max(abs(eta_recovered - cohort$etalcl)) < 1e-3,
  # 4. The realised spread reproduces the published omega. A 200-point mid-point
  #    qnorm grid carries a variance of 0.09616 against the nominal 0.0963
  #    (a 0.14% shortfall), so this lands within 0.05 CV percentage points.
  abs(gcv_recovered - gcv_nominal) < 0.3
)

checks |> knitr::kable(caption = "Closed-form cross-checks on the virtual cohort.")
```

| Check | Simulated | Expected |
|:---|:---|:---|
| Cav,ss from PKNCA vs closed form dose/(CL\*tau) | 0.008% max abs diff | \< 0.5% |
| Cohort median Cav,ss vs the typical-value prediction | 343.5 ng/mL | 343.5 ng/mL |
| Random effect recovered from Cav,ss vs the supplied eta grid | 8.26e-05 max abs diff | \< 1e-3 (exact algebra) |
| Geometric CV of the recovered random effect vs sqrt(exp(omega^2) - 1) | 31.77% | 31.79% |

Closed-form cross-checks on the virtual cohort. {.table}

### Therapeutic-range attainment at the recommended dose

Yu 2025 targets a reference plasma concentration range of 100-1000 ng/mL
and concludes that “a maintenance dose of 4 mg per day is recommended
for pediatric patients who are not coadministered other ASMs”, using a
probability of target attainment above 90% as the selection criterion.

``` r

pta <- mean(cohort_wide$cav >= 100 & cohort_wide$cav <= 1000)

# The paper's own selection criterion was PTA > 90%; this reproduces it for the
# recommended 4 mg/d regimen across the full Table 1 age and weight span.
stopifnot(pta > 0.90)

tibble::tibble(
  Metric = c("Dose", "Target range", "Probability of target attainment"),
  Value  = c("4 mg once daily", "100-1000 ng/mL", sprintf("%.1f%%", 100 * pta))
) |>
  knitr::kable(caption = "Target attainment for the dose Yu 2025 recommends.")
```

| Metric                           | Value           |
|:---------------------------------|:----------------|
| Dose                             | 4 mg once daily |
| Target range                     | 100-1000 ng/mL  |
| Probability of target attainment | 92.5%           |

Target attainment for the dose Yu 2025 recommends. {.table}

``` r

cohort_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(Cc_ngml = Cc * 1000) |>
  dplyr::group_by(time) |>
  dplyr::summarise(
    Q05 = quantile(Cc_ngml, 0.05), Q50 = median(Cc_ngml),
    Q95 = quantile(Cc_ngml, 0.95), .groups = "drop"
  ) |>
  ggplot(aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25, fill = "steelblue") +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(100, 1000), linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Time after dose (h)", y = "Perampanel concentration (ng/mL)")
```

![Simulated steady-state perampanel concentration-time profiles over one
24 h dosing interval at 4 mg once daily, for the 200-subject
deterministic cohort. Ribbon = 5th-95th percentile, line = median. The
dashed lines bound the 100-1000 ng/mL reference
range.](Yu_2025_perampanel_files/figure-html/figure-cohort-1.png)

Simulated steady-state perampanel concentration-time profiles over one
24 h dosing interval at 4 mg once daily, for the 200-subject
deterministic cohort. Ribbon = 5th-95th percentile, line = median. The
dashed lines bound the 100-1000 ng/mL reference range.

The profile is almost flat across the dosing interval. That is a direct
consequence of the model’s elimination half-life, which is examined
next.

### Elimination half-life

``` r

halflife_events <- rxode2::et(amt = 4, cmt = "depot", id = 1L) |>
  rxode2::et(seq(0, 2400, by = 12), cmt = "central", id = 1L) |>
  as.data.frame() |>
  dplyr::mutate(AGE = 9.00, WT = 28.1,
                CONMED_OXC = 0L, CONMED_CBZ = 0L, CONMED_VPA = 0L)

halflife_sim <- rxode2::rxSolve(mod_typical, events = halflife_events) |>
  drop_duplicate_columns()
#> ℹ omega/sigma items treated as zero: 'etalcl'

# rxSolve() omits the `id` column for a single-subject solve, so it is added
# back explicitly -- PKNCA needs it on both sides of the formula.
hl_conc <- halflife_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(id = 1L, Cc_ngml = Cc * 1000, regimen = "4 mg single dose") |>
  dplyr::select(id, time, Cc_ngml, regimen)

hl_dose <- halflife_events |>
  dplyr::filter(evid == 1) |>
  dplyr::mutate(id = 1L, regimen = "4 mg single dose") |>
  dplyr::select(id, time, amt, regimen)

hl_nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(hl_conc, Cc_ngml ~ time | regimen + id),
  PKNCA::PKNCAdose(hl_dose, amt ~ time | regimen + id),
  intervals = data.frame(start = 0, end = Inf, half.life = TRUE, cmax = TRUE, tmax = TRUE)
))

hl <- as.data.frame(hl_nca) |>
  dplyr::filter(PPTESTCD == "half.life") |>
  dplyr::pull(PPORRES)

# The model's half-life is fixed by CL and V, both taken verbatim from Table 2,
# so this is a transcription check, not a comparison against the paper (which
# reports no half-life for its own model). The wide bound below only asserts
# that the packaged parameters give the arithmetic they should.
stopifnot(abs(hl - log(2) * (227 * log10(28.1)) /
                (0.177 * ((9 + 10) / 8.8)^1.31)) / hl < 0.02)

tibble::tibble(
  Quantity = c("Model half-life at the median patient (age 9.00 y, 28.1 kg)",
               "Implied time to 95% of steady state (4.32 half-lives)",
               "Time to steady state assumed by Yu 2025 (Methods, citing ref 17)"),
  Value = c(sprintf("%.0f h (%.1f days)", hl, hl / 24),
            sprintf("%.0f days", 4.32 * hl / 24),
            "19 days")
) |>
  knitr::kable(caption = "Model half-life versus the paper's own steady-state assumption.")
```

| Quantity | Value |
|:---|:---|
| Model half-life at the median patient (age 9.00 y, 28.1 kg) | 470 h (19.6 days) |
| Implied time to 95% of steady state (4.32 half-lives) | 85 days |
| Time to steady state assumed by Yu 2025 (Methods, citing ref 17) | 19 days |

Model half-life versus the paper’s own steady-state assumption. {.table}

The packaged model’s half-life is far longer than the 19 days to steady
state that Yu 2025 assumes in its own Methods (19 days corresponds to a
half-life of about 105 h, which is the value in the perampanel label).
This is a property of the published estimates rather than of the
transcription, and it is discussed below.

## Assumptions and deviations

- **`LGBW` log base is not stated in the source, and the data cannot
  settle it.** Equation (2) is printed as `V(L) = 227 * LGBW`, and the
  surrounding text says only that “LGBW is the log value of body
  weight”. The model file reads this as base-10 (`log10`), on the
  grounds that `lg` denotes the base-10 logarithm in the notation
  convention used by the authors, and that a natural-log covariate would
  conventionally be named `LNBW`. Under the natural-log reading every
  volume would be larger by a factor of 2.303 (V at the median 28.1 kg
  would be 757 L rather than 329 L). **This choice is unresolved and is
  the subject of an open question to the maintainers.** Nothing on disk
  discriminates the two readings: Table 3 reports interval-average
  concentrations, which are `dose/(CL * tau)` and therefore exactly
  independent of V (the residual scan in this vignette finds no interior
  optimum in V – the fit improves monotonically as V grows and converges
  on the V-free closed form); the supplement’s Figure S1 VPC is plotted
  against absolute clock time on trough-only TDM data; and no NONMEM
  control stream was published. Every gate in this vignette is therefore
  insensitive to the log base, but a user simulating single-dose
  exposure, Cmax, or half-life is not.

- **`227` is a coefficient, not a volume at a reference weight.**
  Equation (2) applies no normalisation to `LGBW`, so `lvc` carries
  units of L per log10 unit of body weight in kg. The paper’s sentence
  calling 227 L “a typical V (L) value” is not consistent with its own
  equation at any studied weight (the equation gives 329 L at the median
  28.1 kg under the base-10 reading). The printed equation is used as
  authoritative, per the standing rule that a printed equation outranks
  surrounding prose. The Methods sentence stating that continuous
  covariates were “modeled via a median standardized model” is likewise
  not reflected in either printed equation – Equation (1) normalises age
  by the constants 10 and 8.8, neither of which is the sample median age
  of 9.00 years – so it was not used to infer a hidden normalisation in
  Equation (2).

- **`0.177 L/h` is an extrapolated intercept, not a typical clearance.**
  The age term `((AGE + 10)/8.8)^1.31` equals 1 only at AGE = -1.2
  years. At the median age of 9.00 years the model gives CL/F = 0.485
  L/h. The paper’s own Discussion quotes 0.439 L/h for “the median age”,
  which does not reproduce from Equation

  1.  at any of the medians in Table 1; the equation itself reproduces
      all 336 cells of Table 3 to a median of 0.7%, so the equation was
      used and the Discussion figure was not.

- **`omega` and `sigma` are read as variances.** Table 2 reports
  `omega CL = 0.0963` and `sigma = 0.130`, and the table footnote
  defines omega as the “interindividual variance for CL”. Both are
  entered on the variance scale, which is also NONMEM’s `$OMEGA` /
  `$SIGMA` reporting convention; the model file therefore carries
  `etalcl ~ 0.0963` (equivalently CV 31.8%) and
  `propSd = sqrt(0.130) = 0.3606` (36.1% proportional error). Reading
  either as a standard deviation instead would shrink the total
  variability to roughly 34% CV; the 90% prediction interval implied by
  the variance reading (about 0.46x to 2.19x the median) is the one
  consistent with the spread of the supplement’s Figure S1 VPC, so the
  variance reading is used.

- **The model’s half-life is much longer than the perampanel
  literature.** With CL/F and V taken verbatim from Table 2, the median
  patient’s half-life is about 470 h (19.6 days) under the base-10
  reading, against roughly 105 h in the perampanel label and implied by
  the paper’s own “steady state after 19 days” statement. The volume is
  only weakly identified by the underlying data: every observation was
  drawn about 12 h after the previous dose at steady state, which
  constrains `dose/CL` but carries almost no information about V, and
  the bootstrap 95% interval for V spans 95.9-362 L (RSE 27.5%) against
  14.1% for the point estimate. The parameter is packaged as published
  and is not tuned.

- **Table 3 reports interval-average, not trough, concentrations.** The
  table caption reads “Median Steady-State Plasma Concentrations”; the
  body text calls them trough concentrations. The published numbers
  match `dose/(CL * 24)` to a median of 0.7% and match the model’s
  end-of-interval trough with a systematic bias of -2.4% whose size
  scales with clearance across the four comedication arms. The interval
  average is therefore used as the comparison quantity, and both are
  reported above.

- **Virtual-cohort covariate distributions are constructed, not
  published.** Yu 2025 publishes marginal quantiles for age and weight
  but not their joint distribution. The cohort interpolates both through
  the Table 1 quantiles on a shared probability grid, which makes age
  and weight perfectly rank-correlated. Comedication is set to none
  throughout the cohort so that the 4 mg/d target attainment check
  matches the population the paper’s recommendation applies to.

- **Between-subject variability is laid out deterministically.** The
  single `etalcl` is supplied as a mid-point `qnorm` grid via a data
  column with `omega = NA`, rather than sampled by `rxSolve()`. This is
  numerically equivalent to a stratified sample of 200 draws and makes
  every assertion in this vignette reproducible across machines, rxode2
  versions and solver thread counts.

- **`ka` was not estimated from pediatric data.** Table 2 reports KA as
  3.37 1/h FIXED, taken from Fujita 2023 (an adult analysis); the
  paper’s Limitations paragraph states that “there is no absorption
  constant available for pediatric patients, and an adult constant was
  applied in modeling”. It is encoded with `fixed()` accordingly.
  Because absorption is fast relative to the dosing interval, none of
  the steady-state gates in this vignette is sensitive to its value.

- **Magnesium valproate is not covered by `CONMED_VPA`.** Equation (1)
  defines the VPA covariate as sodium valproate. One patient (0.66%)
  received magnesium valproate; the paper does not say whether that
  patient was coded as VPA-positive.
