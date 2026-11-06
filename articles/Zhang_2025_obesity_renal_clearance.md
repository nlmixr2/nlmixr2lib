# Renal clearance of hypothetical drugs in obesity (Zhang 2025)

## Model and source

- Citation: Zhang T, Calvier EAM, Krekels EHJ, Knibbe CAJ. The Impact of
  Drug Properties and Severity of Obesity on Renal Drug Clearance
  Through Glomerular Filtration and Active Tubular Secretion: A
  Systematic Analysis Using PBPK Modeling. Pharm Res.
  2025;42(7):1079-1088. <doi:10.1007/s11095-025-03885-5>. PMID 40571872.
  PMCID PMC12304071. Eqs. 1-3 of the main text; the subject-specific
  equations and the drug- and system-specific Eqs. 1-3 are from the
  Supplementary Material (sections 1 and 2). The blood-weight and
  heart-rate functions inside the printed cardiac-output equation are
  from Berton M, Bettonte S, Stader F, Battegay M, Marzolini C. Clin
  Pharmacokinet. 2022;61(9):1251-1270
  (<doi:10.1007/s40262-022-01132-3>), Table 1.
- Description: PBPK (renal clearance, hypothetical drug). Renal
  clearance (CLr) of a hypothetical small-molecule drug eliminated by
  glomerular filtration (GF) and active tubular secretion (ATS) in
  adults with normal weight, overweight or obesity (BMI 20 to 60 kg/m^2;
  Zhang 2025). CLr = CLGF + CLATS, with CLGF = fu \* GFR and CLATS =
  (Qr - GFR) \* fu \* CLint_sec / (Qr + fu \* CLint_sec / BP). Renal
  blood flow, GFR, kidney weight, albumin, alpha-1 acid glycoprotein and
  hematocrit follow the obese population regressions of Berton 2022 in
  BMI, height, age, sex and body surface area. fu scales with the
  binding protein (albumin or AAG) from its normal-weight value, BP from
  fu, hematocrit and the red-cell partition coefficient Kp, and
  CLint_sec = CLint_ATS \* PTCPGK \* kidney weight \* relative
  transporter activity (rTA). The drug properties (normal-weight fu, Kp,
  CLint_ATS, binding protein) and rTA are fixed placeholders that the
  user replaces. The model has no compartments, ODEs, dosing,
  between-subject variability or residual error; each record gives the
  clearances implied by its covariates and the time column is not used.
- Article: <https://doi.org/10.1007/s11095-025-03885-5>

Zhang 2025 asks how obesity changes the renal clearance (CLr) of drugs
that are filtered at the glomerulus and secreted by tubular
transporters, and whether the two common dose-scaling rules, allometric
scaling with exponent 0.75 (AS0.75) and flat dosing (exponent 0),
predict that change. It uses a static physiologically based renal
clearance model:

- CLr = CLGF + CLATS (Eq. 1), with CLGF = fu x GFR (Eq. 2) and CLATS =
  (Qr - GFR) x fu x CLint_sec / (Qr + fu x CLint_sec / BP) (Eq. 3).
- The subject quantities (renal blood flow Qr, GFR, kidney weight,
  albumin, alpha-1 acid glycoprotein, hematocrit) follow the
  obese-population regressions of Berton 2022 in BMI, height, age and
  sex (Supplementary section 1).
- The drug quantities follow Supplementary Eqs. 1-3: fu scales with the
  concentration of the binding protein relative to the normal-weight
  subject, BP = 1 + hematocrit x (fu x Kp - 1), and CLint_sec =
  CLint_ATS x PTCPGK x kidney weight x rTA, where rTA is the transporter
  activity relative to normal weight.

The model has no time course: every record returns the clearances
implied by its covariates and drug properties. The drug properties
(`fu_nw`, `kp_rbc`, `clint_ats`, `bind_aag`) and `rta` are fixed `ini()`
placeholders. They can be supplied per record as data columns of the
same name, which is how this vignette sweeps the paper’s hypothetical
drugs.

## Population

No patients were studied. The paper defines six typical adults with BMI
20 (the normal-weight reference), 25, 30, 40, 50 and 60 kg/m^2, all 172
cm tall (the NHANES adult median) and 50 years old, with every
sex-dependent quantity averaged over men and women (Methods;
Supplementary section 1). Their body weights are 59.2 to 177.5 kg. The
hypothetical drugs combine a normal-weight fu of 5, 25, 50, 75, 95 or
100%, binding to albumin (HSA) or to AAG, a Kp of 0.35, 1, 2, 3, 4, 5,
10, 20, 30 or 40, and 32 equidistant CLint_ATS values from 2 to 500
uL/min/mg protein. Each obese subject was also given a relative
transporter activity (rTA) of 20, 50, 100, 150, 200 or 250%.

## Source trace

| Quantity | Value / equation | Source |
|----|----|----|
| CLr, CLGF, CLATS | Eqs. 1-3 as in the model description | Main text Eqs. 1-3 |
| Contribution of GF and ATS (%) | 100 x CLGF / CLr, 100 x CLATS / CLr | Main text Eqs. 4-5 |
| Body weight | BMI x (HT/100)^2 | Supplement section 1 |
| Fat mass | (BMI - FFMI) x (HT/100)^2, FFMI = 30.86 - 3.09 Sex - 0.13 HT + 0.14 BW | Supplement section 1 (sign, see Assumptions) |
| Blood weight | 2.58 + 0.057 BW - 0.061 x FM/BW | Berton 2022 Table 1 (fraction, see Assumptions) |
| Cardiac output | HR x (43.13 + 12.96 x blood weight) x (0.65 - 0.0018 BMI) / 1000 | Supplement section 1 |
| `hr` | 80 beats/min | Berton 2022 Table 1 (mean of the truncated normal) |
| Renal blood flow | (20.57 - 1.76 Sex) x CO / 100 | Supplement section 1 |
| HSA, AAG | 47.58 - 0.16 BMI; 0.48 + 0.012 BMI | Supplement section 1 |
| Hematocrit | 0.45 - 0.04 Sex | Supplement section 1 |
| Kidney weight | -0.11 - 0.048 Sex + 0.003 Age + 0.18 BSA | Supplement section 1 |
| GFR | 13.02 - 0.51 Age + 0.93 x FM/BW + 93.56 x Qr | Supplement section 1 |
| fu | 1 / (1 + (1 - fu_nw) x PP / (fu_nw x PP_nw)) | Supplement Eq. 1 |
| CLint_sec | CLint_ATS x PTCPGK x KW x rTA | Supplement Eq. 2 |
| BP | 1 + Hct x (fu x Kp - 1) | Supplement Eq. 3 |
| `ptcpgk` | 99.4 | Supplement section 2 (Simcyp V22) |
| `bmi_nw` | 20 kg/m^2 | Methods (normal-weight subject) |
| `fu_nw`, `kp_rbc`, `clint_ats`, `bind_aag`, `rta` | placeholders 0.5, 1, 50, 0 (HSA), 1 | Methods (swept values) |

## Typical subjects

The kidney-weight regression needs body surface area, but neither paper
prints the BSA formula. Berton 2022 states that it uses the
Ashby-Thompson equations. This vignette uses their sex-specific forms
(0.01624 x W^0.4725 x H^0.5231 for men and 0.01522 x W^0.4921 x H^0.5231
for women, W in kg, H in cm) and averages them, like every other
sex-dependent quantity.

``` r

bsa_ashby <- function(wt, ht) {
  (0.01624 * wt^0.4725 * ht^0.5231 + 0.01522 * wt^0.4921 * ht^0.5231) / 2
}
typical <- function(bmi) {
  wt <- bmi * 1.72^2
  data.frame(BMI = bmi, HT = 172, AGE = 50, SEXF = 0.5, BSA = bsa_ashby(wt, 172))
}
solve_rows <- function(d) {
  d$id <- seq_len(nrow(d))
  d$time <- 0
  d$evid <- 0
  as.data.frame(rxode2::rxSolve(ui, d, returnType = "data.frame"))
}
bmis <- c(20, 25, 30, 40, 50, 60)
kps <- c(0.35, 1, 2, 3, 4, 5, 10, 20, 30, 40)
subjects <- typical(bmis)
```

## Checks against the paper

### Subject-specific parameters (Supplementary section 1)

``` r

phys <- solve_rows(subjects)
printed <- data.frame(
  BMI = bmis,
  bw = c(59.2, 74, 88.8, 118, 148, 178),
  fm = c(14.1, 22.7, 31.4, 48.8, 66.1, 83.4),
  co = c(5.86, 6.29, 6.81, 7.59, 8.31, 8.96),
  qr = c(1.15, 1.24, 1.34, 1.49, 1.64, 1.76),
  hsa = c(44.4, 43.6, 42.8, 41.2, 40.0, 38.0),
  aag = c(0.72, 0.78, 0.84, 0.96, 1.08, 1.20),
  hct = 0.43,
  kw = c(0.31, 0.35, 0.38, 0.44, 0.48, 0.53),
  gfr = c(95.7, 104, 113, 128, 141, 153)
)
phys_cmp <- printed |>
  pivot_longer(-BMI, names_to = "quantity", values_to = "printed") |>
  left_join(
    phys |>
      select(BMI, all_of(setdiff(names(printed), "BMI"))) |>
      pivot_longer(-BMI, names_to = "quantity", values_to = "model"),
    by = c("BMI", "quantity")
  ) |>
  mutate(pct_diff = 100 * (model - printed) / printed)
phys_cmp |>
  select(quantity, BMI, model) |>
  pivot_wider(names_from = BMI, values_from = model) |>
  knitr::kable(digits = 3, caption = "Model values for the six typical subjects")
```

| quantity |     20 |      25 |      30 |      40 |      50 |      60 |
|:---------|-------:|--------:|--------:|--------:|--------:|--------:|
| bw       | 59.168 |  73.960 |  88.752 | 118.336 | 147.920 | 177.504 |
| fm       | 14.086 |  22.752 |  31.417 |  48.748 |  66.079 |  83.410 |
| co       |  5.899 |   6.338 |   6.763 |   7.568 |   8.310 |   8.989 |
| qr       |  1.161 |   1.248 |   1.332 |   1.490 |   1.636 |   1.770 |
| hsa      | 44.380 |  43.580 |  42.780 |  41.180 |  39.580 |  37.980 |
| aag      |  0.720 |   0.780 |   0.840 |   0.960 |   1.080 |   1.200 |
| hct      |  0.430 |   0.430 |   0.430 |   0.430 |   0.430 |   0.430 |
| kw       |  0.315 |   0.349 |   0.380 |   0.434 |   0.481 |   0.524 |
| gfr      | 96.407 | 104.572 | 112.444 | 127.312 | 141.017 | 153.560 |

Model values for the six typical subjects {.table}

``` r

phys_cmp |>
  group_by(quantity) |>
  summarise(max_abs_pct_diff = max(abs(pct_diff))) |>
  rename("Quantity" = quantity, "Largest |difference| vs Supplement (%)" = max_abs_pct_diff) |>
  knitr::kable(digits = 2)
```

| Quantity | Largest \|difference\| vs Supplement (%) |
|:---------|-----------------------------------------:|
| aag      |                                     0.00 |
| bw       |                                     0.28 |
| co       |                                     0.77 |
| fm       |                                     0.23 |
| gfr      |                                     0.74 |
| hct      |                                     0.00 |
| hsa      |                                     1.05 |
| kw       |                                     1.66 |
| qr       |                                     1.00 |

Every quantity is within 2% of the printed table, and most of the
differences are the two-digit rounding of the table. The cardiac output
alone scatters by up to 0.8% in both directions, which no single
constant explains (see Assumptions).

``` r

stopifnot(max(abs(phys_cmp$pct_diff)) < 2)
```

### Clearance by glomerular filtration (Table S1)

``` r

s1 <- read.csv(text = "
protein,BMI,f5,f50,f95
HSA,25,5.27,52.31,98.58
HSA,30,5.87,57.69,107.83
HSA,40,6.86,66.25,121.77
HSA,50,7.86,74.55,134.71
HSA,60,8.87,82.46,146.42
AAG,25,4.8,49.77,98.09
AAG,30,4.89,52.29,106.75
AAG,40,4.85,54.74,119.34
AAG,50,4.78,56.41,130.71
AAG,60,4.68,57.38,140.68
") |>
  pivot_longer(c(f5, f50, f95), names_to = "fu", values_to = "printed") |>
  mutate(fu_nw = as.numeric(sub("f", "", fu)) / 100)
s1_grid <- s1 |>
  distinct(protein, BMI, fu_nw) |>
  left_join(subjects, by = "BMI") |>
  mutate(bind_aag = as.numeric(protein == "AAG"))
s1_sim <- solve_rows(s1_grid)
s1_cmp <- s1 |>
  mutate(model = s1_sim$cl_gf[match(
    paste(protein, BMI, fu_nw),
    paste(s1_grid$protein, s1_grid$BMI, s1_grid$fu_nw)
  )]) |>
  mutate(pct_diff = 100 * (model - printed) / printed)
s1_cmp |>
  select(protein, BMI, fu_nw, printed, model, pct_diff) |>
  rename(
    "Binding protein" = protein, "fu normal weight" = fu_nw,
    "Table S1 CLGF (mL/min)" = printed, "Model CLGF (mL/min)" = model,
    "Difference (%)" = pct_diff
  ) |>
  knitr::kable(digits = 2)
```

| Binding protein | BMI | fu normal weight | Table S1 CLGF (mL/min) | Model CLGF (mL/min) | Difference (%) |
|:---|---:|---:|---:|---:|---:|
| HSA | 25 | 0.05 | 5.27 | 5.32 | 0.94 |
| HSA | 25 | 0.50 | 52.31 | 52.76 | 0.86 |
| HSA | 25 | 0.95 | 98.58 | 99.43 | 0.87 |
| HSA | 30 | 0.05 | 5.87 | 5.82 | -0.82 |
| HSA | 30 | 0.50 | 57.69 | 57.25 | -0.76 |
| HSA | 30 | 0.95 | 107.83 | 107.02 | -0.76 |
| HSA | 40 | 0.05 | 6.86 | 6.83 | -0.38 |
| HSA | 40 | 0.50 | 66.25 | 66.04 | -0.32 |
| HSA | 40 | 0.95 | 121.77 | 121.38 | -0.32 |
| HSA | 50 | 0.05 | 7.86 | 7.86 | -0.02 |
| HSA | 50 | 0.50 | 74.55 | 74.54 | -0.01 |
| HSA | 50 | 0.95 | 134.71 | 134.69 | -0.01 |
| HSA | 60 | 0.05 | 8.87 | 8.90 | 0.30 |
| HSA | 60 | 0.50 | 82.46 | 82.75 | 0.35 |
| HSA | 60 | 0.95 | 146.42 | 146.94 | 0.36 |
| AAG | 25 | 0.05 | 4.80 | 4.85 | 0.94 |
| AAG | 25 | 0.50 | 49.77 | 50.19 | 0.85 |
| AAG | 25 | 0.95 | 98.09 | 98.93 | 0.86 |
| AAG | 30 | 0.05 | 4.89 | 4.85 | -0.74 |
| AAG | 30 | 0.50 | 52.29 | 51.90 | -0.75 |
| AAG | 30 | 0.95 | 106.75 | 105.94 | -0.76 |
| AAG | 40 | 0.05 | 4.85 | 4.83 | -0.32 |
| AAG | 40 | 0.50 | 54.74 | 54.56 | -0.32 |
| AAG | 40 | 0.95 | 119.34 | 118.96 | -0.32 |
| AAG | 50 | 0.05 | 4.78 | 4.78 | 0.00 |
| AAG | 50 | 0.50 | 56.41 | 56.41 | -0.01 |
| AAG | 50 | 0.95 | 130.71 | 130.70 | -0.01 |
| AAG | 60 | 0.05 | 4.68 | 4.70 | 0.44 |
| AAG | 60 | 0.50 | 57.38 | 57.59 | 0.36 |
| AAG | 60 | 0.95 | 140.68 | 141.18 | 0.35 |

``` r

stopifnot(max(abs(s1_cmp$pct_diff)) < 2)
```

### Clearance by active tubular secretion (Table S2)

Table S2 prints each cell as a central value with a range in brackets.
Its caption says the drugs have Kp values of 0.35, 1 and 4, but those
three Kp values cannot produce the printed upper limits: at fu 50% and
CLint_ATS 500 the bracket reaches 4545.7 mL/min at BMI 25, against 1364
mL/min at Kp 4. The cells are reproduced instead by the **median,
minimum and maximum over all ten Kp values** of the Methods, which is
what is compared here.

``` r

s2 <- read.csv(text = "
protein,BMI,clint,fu,med,lo,hi
HSA,25,5,5,8,7.99,8.05
HSA,25,5,50,76.67,72.28,79.58
HSA,25,5,95,143.16,127.34,149.94
HSA,25,50,5,73.06,71.97,77.05
HSA,25,50,50,549.8,382.95,744.97
HSA,25,50,95,971.74,527.03,1400.99
HSA,25,500,5,391.35,362.02,541.51
HSA,25,500,50,1438.67,671.63,4545.7
HSA,25,500,95,2313.58,768.11,8457.57
HSA,30,5,5,8.84,8.83,8.9
HSA,30,5,50,84.06,79.17,87.28
HSA,30,5,95,155.68,138.4,163.09
HSA,30,50,5,80.66,79.41,85.17
HSA,30,50,50,601.45,417.22,816.73
HSA,30,50,95,1055.18,571.27,1523.38
HSA,30,500,5,428.46,395.55,596.64
HSA,30,500,50,1567.53,728.08,4974.21
HSA,30,500,95,2506.29,831.28,9180.75
HSA,40,5,5,10.5,10.48,10.57
HSA,40,5,50,97.96,91.99,101.85
HSA,40,5,95,178.42,158.1,187.16
HSA,40,50,5,95.24,93.64,100.97
HSA,40,50,50,694.1,475.47,951.24
HSA,40,50,95,1198.53,642.9,1744.9
HSA,40,500,5,493.26,453.09,698.63
HSA,40,500,50,1776.96,815.39,5730.58
HSA,40,500,95,2806.61,927.24,10404.15
HSA,50,5,5,12.06,12.03,12.14
HSA,50,5,50,110.37,103.43,114.84
HSA,50,5,95,197.67,174.9,207.47
HSA,50,50,5,108.92,106.95,115.85
HSA,50,50,50,778.09,528.12,1071.71
HSA,50,50,95,1322.93,706.67,1932.72
HSA,50,500,5,553.64,506.14,795.56
HSA,50,500,50,1973.79,896.04,6426.37
HSA,50,500,95,3080.1,1015.4,11474.47
HSA,60,5,5,13.69,13.66,13.79
HSA,60,5,50,122.82,114.83,127.9
HSA,60,5,95,216.2,190.97,227.08
HSA,60,50,5,123.1,120.69,131.42
HSA,60,50,50,860.79,578.44,1192.44
HSA,60,50,95,1440.52,765.71,2113.41
HSA,60,500,5,613,557.49,894.86
HSA,60,500,50,2161.01,970.12,7110.4
HSA,60,500,95,3330.65,1095.38,12481.3
AAG,25,5,5,7.29,7.28,7.33
AAG,25,5,50,73,69.06,75.71
AAG,25,5,95,142.44,126.78,149.19
AAG,25,50,5,67.04,66.19,70.3
AAG,25,50,50,526.12,372.73,708.88
AAG,25,50,95,967.26,525.75,1393.93
AAG,25,500,5,371.19,346.59,499.33
AAG,25,500,50,1389.8,665.27,4330.43
AAG,25,500,95,2304.24,767.18,8415.52
AAG,30,5,5,7.39,7.38,7.42
AAG,30,5,50,76.32,72.37,79.11
AAG,30,5,95,154.13,137.18,161.45
AAG,30,50,5,68.22,67.45,71.26
AAG,30,50,50,551.58,395.52,740.64
AAG,30,50,95,1045.47,568.52,1508.11
AAG,30,500,5,386.41,363.07,509.69
AAG,30,500,50,1464.8,714.63,4521.12
AAG,30,500,95,2486.14,829.26,9089.84
AAG,40,5,5,7.45,7.44,7.48
AAG,40,5,50,81.2,77.26,84.16
AAG,40,5,95,174.9,155.35,183.43
AAG,40,50,5,69.24,68.6,71.92
AAG,40,50,50,587.03,428.78,786.67
AAG,40,50,95,1176.7,636.8,1710.19
AAG,40,500,5,405.22,384.39,518.46
AAG,40,500,50,1559.57,786.76,4761.07
AAG,40,500,95,2761.9,922.79,10199.82
AAG,50,5,5,7.37,7.36,7.39
AAG,50,5,50,83.93,80.15,86.92
AAG,50,5,95,191.88,170.37,201.32
AAG,50,50,5,68.85,68.34,71.17
AAG,50,50,50,609.53,453.96,812.12
AAG,50,50,95,1287.1,696.7,1875.53
AAG,50,500,5,416.78,398.41,518.92
AAG,50,500,50,1633.57,850.77,4903.67
AAG,50,500,95,3007.09,1008.14,11139.21
AAG,60,5,5,7.27,7.27,7.3
AAG,60,5,50,86.02,82.41,89.03
AAG,60,5,95,207.83,184.44,218.18
AAG,60,50,5,68.32,67.88,70.34
AAG,60,50,50,627.07,475.02,831.39
AAG,60,50,95,1388.94,751.46,2030.7
AAG,60,500,5,425.24,408.95,517.82
AAG,60,500,50,1692.53,907.27,5003.7
AAG,60,500,95,3226.26,1085.02,11999.04
")

# The representative drugs at all ten Kp values, in every typical subject
# (BMI 20 is the normal-weight reference for Table S3 below).
rep_grid <- expand.grid(
  BMI = bmis, fu_nw = c(0.05, 0.5, 0.95), bind_aag = c(0, 1),
  clint_ats = c(5, 50, 500), kp_rbc = kps
) |>
  left_join(subjects, by = "BMI")
rep_sim <- solve_rows(rep_grid)

s2_model <- rep_sim |>
  filter(BMI > 20) |>
  mutate(protein = ifelse(bind_aag == 1, "AAG", "HSA"), fu = round(100 * fu_nw)) |>
  group_by(protein, BMI, clint = clint_ats, fu) |>
  summarise(
    med_m = median(cl_ats), lo_m = min(cl_ats), hi_m = max(cl_ats),
    .groups = "drop"
  )
s2_cmp <- s2 |>
  left_join(s2_model, by = c("protein", "BMI", "clint", "fu")) |>
  mutate(
    d_med = 100 * (med_m - med) / med,
    d_lo = 100 * (lo_m - lo) / lo,
    d_hi = 100 * (hi_m - hi) / hi
  )
s2_cmp |>
  summarise(across(c(d_med, d_lo, d_hi), ~ max(abs(.x)))) |>
  rename(
    "Median, largest |difference| (%)" = d_med,
    "Minimum, largest |difference| (%)" = d_lo,
    "Maximum, largest |difference| (%)" = d_hi
  ) |>
  knitr::kable(digits = 2, caption = "Model versus Table S2, all 90 cells (rTA 100%)")
```

| Median, largest \|difference\| (%) | Minimum, largest \|difference\| (%) | Maximum, largest \|difference\| (%) |
|---:|---:|---:|
| 0.65 | 0.72 | 0.63 |

Model versus Table S2, all 90 cells (rTA 100%) {.table}

``` r

stopifnot(
  max(abs(s2_cmp$d_med)) < 2,
  max(abs(s2_cmp$d_lo)) < 2,
  max(abs(s2_cmp$d_hi)) < 2
)
```

### Allometric exponent for exact scaling (Table S3)

Eq. 8 gives the exponent that would scale CLr from the normal-weight
subject exactly: ln(CLr obese / CLr normal weight) / ln(BW obese / BW
normal weight). As for Table S2, each cell is the median and range over
the ten Kp values, with rTA 100%.

``` r

s3 <- read.csv(text = "
protein,BMI,clint,fu,med,lo,hi
HSA,25,5,5,0.49,0.49,0.49
HSA,25,5,50,0.45,0.45,0.46
HSA,25,5,95,0.42,0.41,0.42
HSA,25,50,5,0.51,0.51,0.52
HSA,25,50,50,0.44,0.40,0.48
HSA,25,50,95,0.41,0.37,0.45
HSA,25,500,5,0.43,0.42,0.48
HSA,25,500,50,0.37,0.34,0.44
HSA,25,500,95,0.34,0.33,0.40
HSA,30,5,5,0.53,0.52,0.53
HSA,30,5,50,0.48,0.48,0.49
HSA,30,5,95,0.44,0.44,0.45
HSA,30,50,5,0.53,0.52,0.54
HSA,30,50,50,0.47,0.44,0.49
HSA,30,50,95,0.43,0.40,0.45
HSA,30,500,5,0.46,0.45,0.50
HSA,30,500,50,0.42,0.39,0.46
HSA,30,500,95,0.39,0.38,0.42
HSA,40,5,5,0.55,0.55,0.55
HSA,40,5,50,0.49,0.49,0.50
HSA,40,5,95,0.45,0.44,0.45
HSA,40,50,5,0.55,0.54,0.56
HSA,40,50,50,0.48,0.45,0.51
HSA,40,50,95,0.43,0.41,0.46
HSA,40,500,5,0.47,0.46,0.52
HSA,40,500,50,0.42,0.39,0.47
HSA,40,500,95,0.39,0.38,0.43
HSA,50,5,5,0.56,0.56,0.56
HSA,50,5,50,0.50,0.50,0.51
HSA,50,5,95,0.45,0.44,0.45
HSA,50,50,5,0.56,0.56,0.57
HSA,50,50,50,0.49,0.45,0.51
HSA,50,50,95,0.44,0.41,0.46
HSA,50,500,5,0.48,0.47,0.54
HSA,50,500,50,0.44,0.40,0.48
HSA,50,500,95,0.40,0.39,0.43
HSA,60,5,5,0.58,0.58,0.58
HSA,60,5,50,0.52,0.51,0.52
HSA,60,5,95,0.45,0.45,0.46
HSA,60,50,5,0.58,0.57,0.59
HSA,60,50,50,0.50,0.46,0.53
HSA,60,50,95,0.44,0.42,0.46
HSA,60,500,5,0.50,0.48,0.55
HSA,60,500,50,0.45,0.41,0.50
HSA,60,500,95,0.40,0.39,0.44
AAG,25,5,5,0.07,0.07,0.07
AAG,25,5,50,0.23,0.23,0.24
AAG,25,5,95,0.39,0.38,0.40
AAG,25,50,5,0.12,0.11,0.13
AAG,25,50,50,0.25,0.25,0.27
AAG,25,50,95,0.39,0.35,0.42
AAG,25,500,5,0.19,0.11,0.22
AAG,25,500,50,0.22,0.20,0.28
AAG,25,500,95,0.32,0.32,0.38
AAG,30,5,5,0.08,0.08,0.08
AAG,30,5,50,0.24,0.24,0.25
AAG,30,5,95,0.42,0.41,0.42
AAG,30,50,5,0.11,0.09,0.12
AAG,30,50,50,0.25,0.25,0.29
AAG,30,50,95,0.41,0.39,0.43
AAG,30,500,5,0.20,0.11,0.23
AAG,30,500,50,0.25,0.22,0.33
AAG,30,500,95,0.37,0.36,0.40
AAG,40,5,5,0.05,0.05,0.05
AAG,40,5,50,0.22,0.22,0.23
AAG,40,5,95,0.42,0.41,0.42
AAG,40,50,5,0.08,0.07,0.09
AAG,40,50,50,0.24,0.23,0.28
AAG,40,50,95,0.41,0.39,0.43
AAG,40,500,5,0.19,0.09,0.22
AAG,40,500,50,0.23,0.20,0.33
AAG,40,500,95,0.37,0.36,0.40
AAG,50,5,5,0.02,0.02,0.02
AAG,50,5,50,0.20,0.20,0.21
AAG,50,5,95,0.42,0.41,0.42
AAG,50,50,5,0.06,0.04,0.06
AAG,50,50,50,0.22,0.21,0.27
AAG,50,50,95,0.41,0.39,0.43
AAG,50,500,5,0.17,0.07,0.20
AAG,50,500,50,0.23,0.19,0.33
AAG,50,500,95,0.37,0.37,0.40
AAG,60,5,5,0.00,0.0053,0.01
AAG,60,5,50,0.19,0.19,0.20
AAG,60,5,95,0.42,0.41,0.42
AAG,60,50,5,0.04,0.02,0.05
AAG,60,50,50,0.21,0.20,0.27
AAG,60,50,95,0.41,0.40,0.43
AAG,60,500,5,0.16,0.06,0.19
AAG,60,500,50,0.22,0.17,0.33
AAG,60,500,95,0.38,0.37,0.40
")

ref <- rep_sim |>
  filter(BMI == 20) |>
  select(fu_nw, bind_aag, clint_ats, kp_rbc, clr0 = clr, bw0 = bw)
rep_rel <- rep_sim |>
  filter(BMI > 20) |>
  left_join(ref, by = c("fu_nw", "bind_aag", "clint_ats", "kp_rbc")) |>
  mutate(expo = log(clr / clr0) / log(bw / bw0))
s3_cmp <- rep_rel |>
  mutate(protein = ifelse(bind_aag == 1, "AAG", "HSA"), fu = round(100 * fu_nw)) |>
  group_by(protein, BMI, clint = clint_ats, fu) |>
  summarise(
    med_m = median(expo), lo_m = min(expo), hi_m = max(expo),
    .groups = "drop"
  ) |>
  right_join(s3, by = c("protein", "BMI", "clint", "fu")) |>
  mutate(d_med = med_m - med, d_lo = lo_m - lo, d_hi = hi_m - hi)
s3_cmp |>
  summarise(across(c(d_med, d_lo, d_hi), ~ max(abs(.x)))) |>
  rename(
    "Median, largest |difference|" = d_med,
    "Minimum, largest |difference|" = d_lo,
    "Maximum, largest |difference|" = d_hi
  ) |>
  knitr::kable(digits = 3, caption = "Model versus Table S3, all 90 cells (exponent units)")
```

| Median, largest \|difference\| | Minimum, largest \|difference\| | Maximum, largest \|difference\| |
|---:|---:|---:|
| 0.033 | 0.034 | 0.032 |

Model versus Table S3, all 90 cells (exponent units) {.table}

``` r

s3_cmp |>
  group_by(BMI) |>
  summarise(mean_d_med = mean(d_med)) |>
  rename("Mean difference of the median exponent" = mean_d_med) |>
  knitr::kable(digits = 3)
```

| BMI | Mean difference of the median exponent |
|----:|---------------------------------------:|
|  25 |                                  0.002 |
|  30 |                                 -0.018 |
|  40 |                                 -0.010 |
|  50 |                                 -0.006 |
|  60 |                                 -0.005 |

``` r

stopifnot(
  max(abs(s3_cmp$d_med)) < 0.035,
  max(abs(s3_cmp$d_lo)) < 0.035,
  max(abs(s3_cmp$d_hi)) < 0.035
)
```

The exponents agree within 0.035 but are systematically low at BMI 30.
The exponent divides a log clearance ratio by a small log weight ratio
(ln 1.5 at BMI 30), so it magnifies small clearance differences: a 1%
difference in CLr moves it by about 0.025. The source is the cardiac
output, which the model puts 0.7% above the printed value at BMI 20 and
0.7% below it at BMI 30 (see the first check). Feeding the printed
cardiac outputs into the model, through the heart rate each would
require, closes the gap to the two-decimal rounding of the table:

``` r

hr_eff <- 80 * printed$co / phys$co[match(printed$BMI, phys$BMI)]
names(hr_eff) <- printed$BMI
round(hr_eff, 2)
#>    20    25    30    40    50    60 
#> 79.48 79.39 80.55 80.24 80.00 79.74
rep_sim_co <- solve_rows(rep_grid |> mutate(hr = hr_eff[as.character(BMI)]))
ref_co <- rep_sim_co |>
  filter(BMI == 20) |>
  select(fu_nw, bind_aag, clint_ats, kp_rbc, clr0 = clr, bw0 = bw)
s3_co <- rep_sim_co |>
  filter(BMI > 20) |>
  left_join(ref_co, by = c("fu_nw", "bind_aag", "clint_ats", "kp_rbc")) |>
  mutate(
    expo = log(clr / clr0) / log(bw / bw0),
    protein = ifelse(bind_aag == 1, "AAG", "HSA"), fu = round(100 * fu_nw)
  ) |>
  group_by(protein, BMI, clint = clint_ats, fu) |>
  summarise(
    med_m = median(expo), lo_m = min(expo), hi_m = max(expo),
    s2_med_m = median(cl_ats), .groups = "drop"
  ) |>
  right_join(s3, by = c("protein", "BMI", "clint", "fu")) |>
  left_join(
    s2 |> select(protein, BMI, clint, fu, s2_med = med),
    by = c("protein", "BMI", "clint", "fu")
  )
s3_co_max <- c(
  s3_median = max(abs(s3_co$med_m - s3_co$med)),
  s3_minimum = max(abs(s3_co$lo_m - s3_co$lo)),
  s3_maximum = max(abs(s3_co$hi_m - s3_co$hi)),
  s2_median_pct = max(abs(100 * (s3_co$s2_med_m - s3_co$s2_med) / s3_co$s2_med))
)
round(s3_co_max, 3)
#>     s3_median    s3_minimum    s3_maximum s2_median_pct 
#>         0.009         0.010         0.010         0.551
stopifnot(
  all(s3_co_max[1:3] < 0.011),
  s3_co_max[4] < 1
)
```

### Relative renal clearance (Figure 3)

``` r

rtas <- c(0.2, 0.5, 1, 1.5, 2, 2.5)
fig_grid <- expand.grid(
  BMI = bmis[-1], fu_nw = c(0.05, 0.5, 0.95), bind_aag = c(0, 1),
  clint_ats = c(5, 50, 500), kp_rbc = 1, rta = rtas
) |>
  left_join(subjects, by = "BMI")
fig_ref <- rep_sim |>
  filter(BMI == 20, kp_rbc == 1) |>
  select(fu_nw, bind_aag, clint_ats, clr0 = clr, bw0 = bw)
fig_sim <- solve_rows(fig_grid) |>
  left_join(fig_ref, by = c("fu_nw", "bind_aag", "clint_ats")) |>
  mutate(
    rel = 100 * clr / clr0,
    pe_as075 = 100 * (clr0 * (bw / bw0)^0.75 - clr) / clr,
    pe_flat = 100 * (clr0 - clr) / clr,
    protein = ifelse(bind_aag == 1, "AAG-bound", "HSA-bound"),
    fu_lab = paste0("fu nw ", 100 * fu_nw, "%"),
    clint_lab = factor(paste0("CLint_ATS ", clint_ats), levels = paste0("CLint_ATS ", c(5, 50, 500))),
    rTA = factor(paste0(100 * rta, "%"), levels = paste0(100 * rtas, "%"))
  )
ggplot(fig_sim, aes(BMI, rel, colour = rTA)) +
  geom_point() +
  geom_hline(yintercept = c(70, 130), linetype = "dashed") +
  facet_grid(clint_lab ~ protein + fu_lab, scales = "free_y") +
  scale_colour_viridis_d() +
  labs(x = "Body mass index (kg/m^2)", y = "Relative obese CLr (% of normal-weight CLr)") +
  theme_bw()
```

![Replicates Figure 3 of Zhang 2025: CLr in the obese subject as a
percentage of CLr in the normal-weight subject, Kp = 1, for each rTA.
Dashed lines mark 70% and
130%.](Zhang_2025_obesity_renal_clearance_files/figure-html/fig3-1.png)

Replicates Figure 3 of Zhang 2025: CLr in the obese subject as a
percentage of CLr in the normal-weight subject, Kp = 1, for each rTA.
Dashed lines mark 70% and 130%.

Points read from Figure 3 (rTA 250% unless stated) against the model:

``` r

fig3_read <- data.frame(
  protein = c("HSA-bound", "HSA-bound", "HSA-bound", "AAG-bound", "AAG-bound", "AAG-bound"),
  fu_nw = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.95),
  clint_ats = c(5, 5, 50, 5, 500, 5),
  rta = c(2.5, 0.2, 2.5, 2.5, 0.2, 2.5),
  BMI = c(60, 25, 60, 60, 60, 60),
  read = c(355, 60, 370, 190, 35, 260)
)
fig3_cmp <- fig3_read |>
  left_join(
    fig_sim |> select(protein, fu_nw, clint_ats, rta, BMI, rel),
    by = c("protein", "fu_nw", "clint_ats", "rta", "BMI")
  )
fig3_cmp |>
  mutate(fu_nw = 100 * fu_nw, rta = 100 * rta) |>
  rename(
    "Binding" = protein, "fu nw (%)" = fu_nw, "CLint_ATS" = clint_ats, "rTA (%)" = rta,
    "Read from Figure 3 (%)" = read, "Model (%)" = rel
  ) |>
  knitr::kable(digits = 0)
```

| Binding   | fu nw (%) | CLint_ATS | rTA (%) | BMI | Read from Figure 3 (%) | Model (%) |
|:----------|----------:|----------:|--------:|----:|-----------------------:|----------:|
| HSA-bound |         5 |         5 |     250 |  60 |                    355 |       354 |
| HSA-bound |         5 |         5 |      20 |  25 |                     60 |        58 |
| HSA-bound |         5 |        50 |     250 |  60 |                    370 |       379 |
| AAG-bound |         5 |         5 |     250 |  60 |                    190 |       189 |
| AAG-bound |         5 |       500 |      20 |  60 |                     35 |        38 |
| AAG-bound |        95 |         5 |     250 |  60 |                    260 |       257 |

``` r

stopifnot(max(abs(fig3_cmp$rel - fig3_cmp$read)) < 15)
```

### Prediction error of AS0.75 and flat dosing (Figure 4)

The prediction error (Eq. 7) compares the CLr scaled from the
normal-weight subject (Eq. 6) with the model CLr in the obese subject.

``` r

fig4 <- fig_sim |>
  select(protein, fu_lab, clint_lab, rTA, BMI, pe_as075, pe_flat) |>
  pivot_longer(c(pe_as075, pe_flat), names_to = "method", values_to = "pe") |>
  mutate(method = ifelse(method == "pe_as075", "AS0.75", "Flat dosing")) |>
  filter(abs(pe) <= 100)
ggplot(fig4, aes(BMI, pe, colour = rTA)) +
  geom_point() +
  geom_hline(yintercept = c(-30, 30), linetype = "dashed") +
  facet_grid(method + clint_lab ~ protein + fu_lab) +
  scale_colour_viridis_d() +
  labs(x = "Body mass index (kg/m^2)", y = "Percentage prediction error (%)") +
  theme_bw()
```

![Replicates Figure 4 of Zhang 2025: percentage prediction error of CLr
scaled by AS0.75 (top) and by flat dosing (bottom), Kp = 1. Dashed lines
mark +/-30%; values outside +/-100% are not shown, as in the
paper.](Zhang_2025_obesity_renal_clearance_files/figure-html/fig4-1.png)

Replicates Figure 4 of Zhang 2025: percentage prediction error of CLr
scaled by AS0.75 (top) and by flat dosing (bottom), Kp = 1. Dashed lines
mark +/-30%; values outside +/-100% are not shown, as in the paper.

### Headline results over the full drug grid

The full grid (6 normal-weight fu x 2 binding proteins x 10 Kp x 32
CLint_ATS values, in five obese subjects under six rTA values) tests the
Results statements. CLint_ATS is taken as 32 linearly equidistant values
from 2 to 500 uL/min/mg protein.

``` r

full <- expand.grid(
  BMI = bmis, fu_nw = c(0.05, 0.25, 0.5, 0.75, 0.95, 1), bind_aag = c(0, 1),
  kp_rbc = kps, clint_ats = seq(2, 500, length.out = 32), rta = rtas
) |>
  filter(BMI > 20 | rta == 1) |>
  left_join(subjects, by = "BMI")
full_sim <- solve_rows(full)
full_ref <- full_sim |>
  filter(BMI == 20) |>
  select(fu_nw, bind_aag, kp_rbc, clint_ats, clr0 = clr, bw0 = bw)
full_rel <- full_sim |>
  filter(BMI > 20) |>
  left_join(full_ref, by = c("fu_nw", "bind_aag", "kp_rbc", "clint_ats")) |>
  mutate(
    rel = clr / clr0,
    expo = log(clr / clr0) / log(bw / bw0),
    pe_as075 = 100 * (clr0 * (bw / bw0)^0.75 - clr) / clr,
    pe_flat = 100 * (clr0 - clr) / clr
  )

# rTA 100%: largest change in CLr and largest |PE| per BMI
full_rel |>
  filter(rta == 1) |>
  group_by(BMI) |>
  summarise(
    change = 100 * max(abs(rel - 1)),
    as075 = max(abs(pe_as075)),
    flat = max(abs(pe_flat))
  ) |>
  rename(
    "Largest |change in CLr| (%)" = change,
    "Largest |PE| AS0.75 (%)" = as075, "Largest |PE| flat (%)" = flat
  ) |>
  knitr::kable(digits = 1)
```

| BMI | Largest \|change in CLr\| (%) | Largest \|PE\| AS0.75 (%) | Largest \|PE\| flat (%) |
|---:|---:|---:|---:|
| 25 | 12.3 | 16.7 | 10.9 |
| 30 | 24.0 | 33.0 | 19.4 |
| 40 | 46.5 | 65.2 | 31.7 |
| 50 | 68.4 | 97.5 | 40.6 |
| 60 | 90.4 | 130.3 | 47.5 |

``` r


expo_rta1 <- full_rel |>
  filter(rta == 1) |>
  group_by(bind_aag) |>
  summarise(min = min(expo), max = max(expo))
expo_rta1
#> # A tibble: 2 × 3
#>   bind_aag      min   max
#>      <dbl>    <dbl> <dbl>
#> 1        0  0.330   0.586
#> 2        1 -0.00927 0.452
range(full_rel$expo)
#> [1] -6.270954  4.171533

bmi25 <- full_rel |> filter(BMI == 25, rta == 1)
stopifnot(
  # 'When BMI is below 30 and rTA remains unchanged ... differed less than 30%'
  max(abs(bmi25$rel - 1)) < 0.3,
  # '... both scaling methods were systematically accurate'
  max(abs(bmi25$pe_as075)) < 30,
  max(abs(bmi25$pe_flat)) < 30,
  # 'exponent ... between 0.33 and 0.59 for drugs binding to HSA'
  abs(expo_rta1$min[expo_rta1$bind_aag == 0] - 0.33) < 0.01,
  abs(expo_rta1$max[expo_rta1$bind_aag == 0] - 0.59) < 0.01,
  # 'range of allometric exponents ... falls within -6.26 to 4.17'
  abs(min(full_rel$expo) - (-6.26)) < 0.05,
  abs(max(full_rel$expo) - 4.17) < 0.05
)
```

The model gives a change in CLr of at most 12% and prediction errors
within 17% for both scaling rules at BMI 25 with unchanged transporter
activity, as the Results state. At BMI 30 the AS0.75 error already
exceeds 30% for some drugs. With rTA 100% the exact-scaling exponent
spans 0.33 to 0.59 for albumin-bound drugs, as printed. For AAG-bound
drugs the model spans -0.009 to 0.45, against 0.0053 to 0.43 in the
text. The text range matches the Table S3 cells (fu 5, 50 and 95%;
CLint_ATS 5, 50 and 500). The wider model range comes from grid drugs
outside those cells: the minimum from CLint_ATS 2 (fu 5%, Kp 40, BMI 60)
and the maximum from a normal-weight fu of 100% (Kp 40, BMI 60). Over
every scenario the exponent spans -6.27 to 4.17, against -6.26 to 4.17
in the Results.

## PKNCA

Not applicable. The model returns renal clearances, not
concentration-time profiles, so there is nothing to run NCA on. The
comparisons above with the supplementary table of subject parameters,
Tables S1 to S3, Figure 3 and the Results statements replace it.

## Assumptions and deviations

- **Sex-averaged typical subject.** Berton 2022 codes Sex as 0 for men
  and 1 for women, the same as `SEXF`. The paper averaged every
  sex-dependent value over the two sexes. Sex enters each regression
  linearly, so `SEXF = 0.5` reproduces that average exactly; computing
  men and women separately and averaging the outputs differs by less
  than 0.01 mL/min in GFR.
- **Fat-mass sign.** The Supplement prints fat mass as
  `[(30.86 - 3.09 Sex - 0.13 HT + 0.14 BW) - BMI] x (HT/100)^2`, which
  is negative for every typical subject. The printed fat masses (14.1 to
  83.4 kg) are `(BMI - FFMI) x (HT/100)^2`, which the model uses.
- **Blood weight and the adipose fraction.** The printed cardiac-output
  equation contains a blood weight that the Supplement does not define.
  It is taken from Berton 2022 Table 1,
  `2.58 + 0.057 BW - 0.061 x adipose weight percentage`, with the
  adipose share entered as a fraction (FM/BW), as in the printed GFR
  equation. As a percentage it would give a cardiac output of 4.98 L/min
  at BMI 20 instead of the printed 5.86 L/min.
- **Heart rate.** Berton 2022 samples heart rate from a truncated normal
  distribution with mean 80 beats/min; the model uses 80. The model
  cardiac outputs scatter around the printed ones by -0.8% to +0.7%,
  with no systematic trend in BMI, which suggests the printed values
  carry some sampling or rounding of their own. Renal blood flow and GFR
  inherit this difference (at most 1%). The heart rates that would
  reproduce the printed cardiac outputs are 79.4 to 80.6 beats/min. The
  difference is small for clearances but visible in the Table S3
  exponents at BMI 30 (up to 0.034); with the printed cardiac outputs
  the exponents match within 0.01.
- **Kidney-weight sex coefficient.** The Supplement prints
  `-0.048 x Sex` and Berton 2022 Table 1 prints `-0.48 x Sex`. The
  printed kidney weights need -0.048, which the model uses.
- **Body surface area.** Neither paper prints the BSA formula for the
  typical subjects. BSA is a covariate of the model. In this vignette it
  is the average of the male and female Ashby-Thompson equations
  (Obesity 2020;28:706-717, <doi:10.1002/oby.22743>), which Berton 2022
  names as its BSA method; with it the model kidney weights match the
  printed ones within the two-digit rounding. The Ashby-Thompson
  abstract describes height in metres, but its coefficients give a
  realistic BSA only with height in centimetres, which is what is used.
- **GFR units.** The Supplement writes the renal-blood-flow term as
  `93.56 x Qr / 1000` next to a renal blood flow tabulated in L/min. The
  printed GFRs need Qr in L/min inside the regression (that is, Qr in
  mL/min divided by 1000), and the adipose share as a fraction.
- **PTCPGK units.** PTCPGK = 99.4 is printed without units, and
  CLint_ATS is in uL/min/mg protein. With kidney weight in kg, CLint_ATS
  x PTCPGK x kidney weight is in mL/min, and this reproduces all 90
  cells of Table S2 within 2%.
- **Table S2 and S3 ranges.** The captions say the drugs have Kp values
  of 0.35, 1 and 4, but the printed medians and ranges are those over
  all ten Kp values of the Methods (see above). One Table S3 cell prints
  a median of 0.00 below its range of 0.0053 to 0.01 (AAG, BMI 60, fu
  5%, CLint_ATS 5); this is two-decimal rounding of a median near 0.005.
- **CLint_ATS grid.** The Methods give 32 equidistant values from 2 to
  500 uL/min/mg protein. Linear spacing is assumed; the figures and
  tables use the representative values 5, 50 and 500, which are not on
  that grid.
- **Clearance above renal blood flow.** For drugs with a large Kp (BP
  well above 1) and high CLint_ATS, CLATS, a plasma clearance, can
  exceed renal blood flow. This follows from Eq. 3 and appears in the
  printed Table S2 ranges (up to 12,481 mL/min).
- **Figure 2 and Supplementary Figures S1 and S2** show the same
  quantities as Tables S1 to S3 and Figures 3 and 4, which are
  reproduced here; they are not replotted.
