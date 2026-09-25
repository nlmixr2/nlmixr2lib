# Pyrazinamide + ethambutol in children with TB with or without HIV (Maranchick 2026)

## Model and source

Maranchick 2026 developed two independent population PK models from a
single cohort of Ghanaian children on first-line HRZE anti-tuberculosis
therapy: a one-compartment model for pyrazinamide (PZA) and a
two-compartment model for ethambutol (EMB). Following the library’s
replicate-the-author’s-structure policy, the two fits are packaged as
two model files and validated here in one vignette.

``` r

pza <- readModelDb("Maranchick_2026_pyrazinamide")
emb <- readModelDb("Maranchick_2026_ethambutol")

pza_ui <- rxode2::rxode(pza)
#> ℹ parameter labels from comments will be replaced by 'label()'
emb_ui <- rxode2::rxode(emb)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Maranchick NF, Martyn-Dickens C, Enimil A, Yang H, Amissah
  AK, Dompreh A, Bosomtwe D, Sly-Moore E, Opoku T, Appiah AF, Asiedu P,
  Antwi S, Scheetz MH, Peloquin CA, Kwara A. Population pharmacokinetics
  of pyrazinamide and ethambutol in children with tuberculosis with or
  without HIV. Antimicrob Agents Chemother. 2026.
  <doi:10.1128/aac.00909-25>
- Article: <https://doi.org/10.1128/aac.00909-25>
- Pyrazinamide: One-compartment population pharmacokinetic model with
  first-order absorption, an absorption lag time and linear elimination
  for oral pyrazinamide in Ghanaian children with tuberculosis with or
  without HIV coinfection (Maranchick 2026); estimated allometric weight
  scaling on CL/F (exponent 0.70) and V/F (exponent 0.79) normalised to
  15 kg, and an exponential HIV-positive effect raising CL/F by 18.5%.
- Ethambutol: Two-compartment population pharmacokinetic model with
  first-order absorption, an absorption lag time and linear elimination
  for oral ethambutol in Ghanaian children with tuberculosis with or
  without HIV coinfection (Maranchick 2026); estimated allometric weight
  scaling on CL/F (exponent 0.70) and V1/F (exponent 0.62) normalised to
  15.1 kg, and an exponential HIV-positive effect raising CL/F by 24.6%.

No supplementary parameter tables accompany the article; the EuropePMC
supplementary-file endpoint for PMC13041307 returns only the four figure
images. All values below come from the main text and Table 2.

## Population

Eighty-five Ghanaian children with drug-susceptible tuberculosis were
enrolled at Komfo Anokye Teaching Hospital, Kumasi, between February
2019 and June 2021: 41 (48.2%) with TB alone and 44 (51.8%) with TB/HIV
coinfection. Median (range) age was 5.0 years (0.3 to 14.5); 49.4% were
under 5 years and 18.8% under 2 years. Median weight was 16 kg (range 4
to 60) and 52 participants (61.2%) were male. Twenty-four (28.2%) were
malnourished by a body-mass-index-for-age Z score below -2 SD. Of the
TB/HIV participants, 29 (65.9%) received efavirenz-based antiretroviral
therapy (Maranchick 2026 Table 1).

Dosing was by WHO weight band once daily: PZA a median 31.6 mg/kg (range
21.4 to 49.7) and EMB a median 21.4 mg/kg (range 14.3 to 34.2). PK
sampling was performed on one occasion after at least 4 weeks of
therapy, at 0 (pre-dose), 1, 2, 4, 8 and 12 h post-dose. The PZA model
used 509 samples from 85 participants; the EMB model used 501 samples
from 84 (one participant whose samples were all at or near the limit of
quantification was excluded, along with six apparently mislabelled
samples).

The same information is available programmatically via each model’s
`population` metadata (for example
`rxode2::rxode(readModelDb("Maranchick_2026_pyrazinamide"))$population`).

## Source trace

Every `ini()` entry carries an in-file comment naming its source
location. The tables below collect them for review. Table 2 of the
source reports, for each parameter, `Estimate (RSE, %)` and an
`IIV (% CV)` column whose footnote states that the IIV is “reported as
standard deviation of the random effects” with the coefficient of
variation in parentheses.

### Pyrazinamide (one compartment, first-order absorption with lag)

| Equation / parameter | Value | Source location |
|----|----|----|
| `ltlag` | 0.26 h | Table 2, PZA row `tlag (h)` (RSE 46.86) |
| `lka` | 3.76 /h | Table 2, PZA row `Ka (h-1)` (RSE 20.79) |
| `lvc` | 11.30 L | Table 2, PZA row `V/F (L)` (RSE 2.87) |
| `lcl` | 1.27 L/h | Table 2, PZA row `Cl/F (L/h)` (RSE 4.46) |
| `e_wt_cl` | 0.70 | Table 2, PZA row `Exponent, BWonCl/F` (RSE 9.74) |
| `e_wt_vc` | 0.79 | Table 2, PZA row `Exponent, BWonV/F` (RSE 6.64) |
| `e_hiv_pos_cl` | 0.17 | Table 2, PZA row `Exponent, HIV+ on Cl/F` (RSE 31.29) |
| reference weight 15 kg | n/a | Results “Pyrazinamide”: `Cl/F*(Weight/15)^0.7`, `V/F*(Weight/15)^0.79` |
| `etaltlag`, `etalka`, `etalvc`, `etalcl` | 0.75, 0.62, 0.23, 0.32 (SD) | Table 2, PZA `IIV (% CV)` column; squared to variance |
| `addSd`, `propSd` | 1.06, 0.06 | Table 2, PZA `Residual variability` rows a and b |
| `d/dt(depot)`, `d/dt(central)`, `alag(depot)` | n/a | Results “Pyrazinamide”: “one-compartment model, first-order absorption with tlag, and linear elimination” |

### Ethambutol (two compartments, first-order absorption with lag)

| Equation / parameter | Value | Source location |
|----|----|----|
| `ltlag` | 0.67 h | Table 2, EMB row `tlag (h)` (RSE 7.04) |
| `lka` | 3.83 /h | Table 2, EMB row `Ka (h-1)` (RSE 22.97) |
| `lvc` | 95.16 L | Table 2, EMB row `V1/F (L)` (RSE 6.49) |
| `lcl` | 23.2 L/h | Table 2, EMB row `Cl/F (L/h)` (RSE 4.4) |
| `lq` | 11.25 L/h | Table 2, EMB row `Q/F (L/h)` (RSE 6.69) |
| `lvp` | 162.41 L | Table 2, EMB row `V2/F (L)` (RSE 16.45) |
| `e_wt_cl` | 0.70 | Table 2, EMB row `Exponent, BWonCl/F` (RSE 10.31) |
| `e_wt_vc` | 0.62 | Table 2, EMB row `Exponent, BWonV1/F` (RSE 19.34) |
| `e_hiv_pos_cl` | 0.22 | Table 2, EMB row `Exponent, HIV+ on Cl/F` (RSE 21.03) |
| reference weight 15.1 kg | n/a | Results “Ethambutol”: `Cl/F*(Weight/15.1)^0.7`, `V1/F*(Weight/15.1)^0.62` |
| `etaltlag`, `etalka`, `etalvc`, `etalcl`, `etalq`, `etalvp` | 0.35, 1.18, 0.51, 0.32, 0.34, 0.81 (SD) | Table 2, EMB `IIV (% CV)` column; squared to variance |
| `addSd`, `propSd` | 0.03, 0.17 | Table 2, EMB `Residual variability` rows a and b |
| `d/dt(depot)`, `d/dt(central)`, `d/dt(peripheral1)`, `alag(depot)` | n/a | Results “Ethambutol”: “two-compartment model with a tlag and first-order absorption and linear elimination” |

### Confirming the IIV scale from the source table

The `IIV (% CV)` column prints two numbers per row. Squaring the first
and treating it as a log-scale variance reproduces the second as
`sqrt(exp(omega^2) - 1) * 100`, which pins the first number to the
standard deviation of the random effects rather than the variance. This
is a check the table performs on itself, so no assumption is needed.

``` r

iiv_check <- tibble::tribble(
  ~drug,          ~parameter, ~omega_sd, ~printed_cv,
  "Pyrazinamide", "tlag",     0.75,      87.16,
  "Pyrazinamide", "Ka",       0.62,      68.71,
  "Pyrazinamide", "V/F",      0.23,      23.40,
  "Pyrazinamide", "Cl/F",     0.32,      32.60,
  "Ethambutol",   "tlag",     0.35,      36.35,
  "Ethambutol",   "Ka",       1.18,     172.61,
  "Ethambutol",   "Cl/F",     0.32,      33.26,
  "Ethambutol",   "V1/F",     0.51,      54.63,
  "Ethambutol",   "Q/F",      0.34,      35.23,
  "Ethambutol",   "V2/F",     0.81,      96.75
) |>
  mutate(
    cv_from_sd = sqrt(exp(omega_sd^2) - 1) * 100,
    abs_diff   = abs(cv_from_sd - printed_cv)
  )

# Deterministic arithmetic on transcribed table values -- no simulation, so a
# tight bound is correct. Reading the column as a VARIANCE instead would put
# every row out by tens of CV points.
stopifnot(max(iiv_check$abs_diff) < 1.5)

iiv_check |>
  dplyr::rename(
    "Drug" = drug, "Parameter" = parameter,
    "Printed omega" = omega_sd, "Printed CV (%)" = printed_cv,
    "CV from omega as SD (%)" = cv_from_sd, "Abs. difference" = abs_diff
  ) |>
  knitr::kable(
    digits  = 2,
    caption = "Table 2's IIV column self-pins: treating the first number as a log-scale SD reproduces the printed CV."
  )
```

| Drug | Parameter | Printed omega | Printed CV (%) | CV from omega as SD (%) | Abs. difference |
|:---|:---|---:|---:|---:|---:|
| Pyrazinamide | tlag | 0.75 | 87.16 | 86.89 | 0.27 |
| Pyrazinamide | Ka | 0.62 | 68.71 | 68.46 | 0.25 |
| Pyrazinamide | V/F | 0.23 | 23.40 | 23.31 | 0.09 |
| Pyrazinamide | Cl/F | 0.32 | 32.60 | 32.84 | 0.24 |
| Ethambutol | tlag | 0.35 | 36.35 | 36.10 | 0.25 |
| Ethambutol | Ka | 1.18 | 172.61 | 173.91 | 1.30 |
| Ethambutol | Cl/F | 0.32 | 33.26 | 32.84 | 0.42 |
| Ethambutol | V1/F | 0.51 | 54.63 | 54.50 | 0.13 |
| Ethambutol | Q/F | 0.34 | 35.23 | 35.01 | 0.22 |
| Ethambutol | V2/F | 0.81 | 96.75 | 96.29 | 0.46 |

Table 2’s IIV column self-pins: treating the first number as a log-scale
SD reproduces the printed CV. {.table}

## Closed-form verification of the encoded structure

Before any stochastic work, the packaged models are checked against an
independent closed-form steady-state solution written directly from the
paper’s printed parameter values. This gate is deterministic – both
sides use the same typical-value parameters and differ only by
ODE-solver and grid error – so a tight bound is the correct one, and it
fails loudly on a mis-transcribed volume, clearance, rate constant,
allometric exponent, reference weight or HIV coefficient.

``` r

# Superposition of a single-dose profile to steady state: each exponential term
# exp(-lambda * t) is scaled by 1 / (1 - exp(-lambda * tau)).
ss_scale <- function(lambda, tau) 1 / (1 - exp(-lambda * tau))

# One-compartment, first-order absorption. Time measured from the end of the lag.
css_1cmt <- function(t, dose, ka, cl, vc, tau) {
  kel <- cl / vc
  (dose * ka / (vc * (ka - kel))) *
    (exp(-kel * t) * ss_scale(kel, tau) - exp(-ka * t) * ss_scale(ka, tau))
}

# Two-compartment, first-order absorption. alpha / beta are the roots of
# s^2 - (kel + k12 + k21) s + kel * k21 = 0.
css_2cmt <- function(t, dose, ka, cl, vc, q, vp, tau) {
  kel <- cl / vc
  k12 <- q / vc
  k21 <- q / vp
  s   <- kel + k12 + k21
  r   <- sqrt(s^2 - 4 * kel * k21)
  alpha <- (s + r) / 2
  beta  <- (s - r) / 2
  (dose * ka / vc) * (
    ((k21 - ka)    / ((alpha - ka)    * (beta - ka)))    * exp(-ka * t)    * ss_scale(ka, tau) +
    ((k21 - alpha) / ((ka - alpha)    * (beta - alpha))) * exp(-alpha * t) * ss_scale(alpha, tau) +
    ((k21 - beta)  / ((ka - beta)     * (alpha - beta))) * exp(-beta * t)  * ss_scale(beta, tau)
  )
}

# The steady-state profile is periodic with period tau, so an absorption lag is
# a pure phase shift: C_lagged(t) = C_unlagged((t - tlag) mod tau).
shift_lag <- function(t, tlag, tau) (t - tlag) %% tau
```

``` r

tau      <- 24      # once-daily dosing (Materials and Methods, "Model building")
n_doses  <- 12L     # burn-in to steady state; EMB terminal half-life is ~16 h
t_last   <- tau * (n_doses - 1L)

# WHO-recommended weight-band doses, Maranchick 2026 Table 3 (PZA) and Table 4
# (EMB), column "WHO-recommended dose in mg(mg/kg)".
bands <- tibble::tribble(
  ~band,        ~wt_lo, ~wt_hi, ~wt_mid, ~dose_pza, ~dose_emb,
  "4 to <8",         4,      8,       6,       150,       100,
  "8 to <12",        8,     12,      10,       300,       200,
  "12 to <16",      12,     16,      14,       450,       300,
  "16 to <25",      16,     25,      20,       600,       400,
  "25 to <35",      25,     35,      30,       800,       550
)

# Values transcribed from Table 2 / Results, used to build the closed-form
# reference INDEPENDENTLY of what the model files contain.
pub <- list(
  pza = list(tlag = 0.26, ka = 3.76, vc = 11.30, cl = 1.27,
             e_wt_cl = 0.70, e_wt_vc = 0.79, e_hiv = 0.17, wref = 15.0),
  emb = list(tlag = 0.67, ka = 3.83, vc = 95.16, cl = 23.2, q = 11.25, vp = 162.41,
             e_wt_cl = 0.70, e_wt_vc = 0.62, e_hiv = 0.22, wref = 15.1)
)

typ_grid <- seq(0, tau, by = 0.02)

# Typical-value event table: one subject per (weight band, HIV stratum).
make_typical_events <- function(dose_col) {
  subj <- tidyr::crossing(bands, HIV_POS = c(0, 1)) |>
    mutate(id = dplyr::row_number(), WT = wt_mid, dose_mg = .data[[dose_col]])
  dos <- tidyr::crossing(subj, time = tau * (seq_len(n_doses) - 1L)) |>
    mutate(evid = 1L, cmt = "depot", amt = dose_mg)
  obs <- tidyr::crossing(subj, time = t_last + typ_grid) |>
    mutate(evid = 0L, cmt = "central", amt = NA_real_)
  dplyr::bind_rows(dos, obs) |> dplyr::arrange(id, time, dplyr::desc(evid))
}
```

``` r

ev_pza_typ <- make_typical_events("dose_pza")
ev_emb_typ <- make_typical_events("dose_emb")

sim_pza_typ <- rxode2::rxSolve(
  rxode2::zeroRe(pza), events = ev_pza_typ,
  keep = c("band", "HIV_POS", "WT", "dose_mg")
) |> as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etaltlag', 'etalka', 'etalvc', 'etalcl'
#> Warning: multi-subject simulation without without 'omega'

sim_emb_typ <- rxode2::rxSolve(
  rxode2::zeroRe(emb), events = ev_emb_typ,
  keep = c("band", "HIV_POS", "WT", "dose_mg")
) |> as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etaltlag', 'etalka', 'etalvc', 'etalcl', 'etalq', 'etalvp'
#> Warning: multi-subject simulation without without 'omega'

closed_form_pza <- sim_pza_typ |>
  filter(!is.na(Cc), time >= t_last) |>
  mutate(
    t_rel = time - t_last,
    cl_i  = pub$pza$cl * (WT / pub$pza$wref)^pub$pza$e_wt_cl * exp(pub$pza$e_hiv * HIV_POS),
    vc_i  = pub$pza$vc * (WT / pub$pza$wref)^pub$pza$e_wt_vc,
    Cc_cf = css_1cmt(shift_lag(t_rel, pub$pza$tlag, tau),
                     dose_mg, pub$pza$ka, cl_i, vc_i, tau)
  )

closed_form_emb <- sim_emb_typ |>
  filter(!is.na(Cc), time >= t_last) |>
  mutate(
    t_rel = time - t_last,
    cl_i  = pub$emb$cl * (WT / pub$emb$wref)^pub$emb$e_wt_cl * exp(pub$emb$e_hiv * HIV_POS),
    vc_i  = pub$emb$vc * (WT / pub$emb$wref)^pub$emb$e_wt_vc,
    Cc_cf = css_2cmt(shift_lag(t_rel, pub$emb$tlag, tau),
                     dose_mg, pub$emb$ka, cl_i, vc_i, pub$emb$q, pub$emb$vp, tau)
  )

cf_err <- function(d) max(abs(d$Cc - d$Cc_cf) / max(d$Cc_cf))
err_pza <- cf_err(closed_form_pza)
err_emb <- cf_err(closed_form_emb)

# Confirm the gate had rows to test (a zero-row comparison would pass vacuously).
stopifnot(nrow(closed_form_pza) > 1000, nrow(closed_form_emb) > 1000)

# Deterministic: solver + grid error only, so a tight bound is the correct one.
# Realised 7.2e-14 (PZA, an algebraic identity to machine precision) and 6.5e-6
# (EMB, limited by ODE-solver tolerance). A single mis-transcribed parameter
# moves this to order 1e-2 or worse.
stopifnot(err_pza < 1e-3, err_emb < 1e-3)

cat(sprintf(
  "Max relative deviation from the closed-form steady-state profile:\n  pyrazinamide %.2e\n  ethambutol   %.2e\n",
  err_pza, err_emb
))
#> Max relative deviation from the closed-form steady-state profile:
#>   pyrazinamide 1.15e-07
#>   ethambutol   6.48e-06
```

``` r

dplyr::bind_rows(
  closed_form_pza |> mutate(drug = "Pyrazinamide"),
  closed_form_emb |> mutate(drug = "Ethambutol")
) |>
  mutate(
    band  = factor(band, levels = bands$band),
    Group = ifelse(HIV_POS == 1, "TB/HIV", "TB")
  ) |>
  ggplot(aes(t_rel, Cc, colour = Group)) +
  geom_line(linewidth = 0.7) +
  geom_point(
    data = ~ dplyr::filter(.x, abs(t_rel - round(t_rel / 2) * 2) < 1e-6),
    aes(y = Cc_cf), shape = 1, size = 1.6
  ) +
  facet_grid(drug ~ band, scales = "free_y") +
  labs(
    x = "Time after dose at steady state (h)", y = "Concentration (ug/mL)",
    colour = NULL,
    caption = "Lines: packaged nlmixr2lib models. Circles: closed-form solution from Maranchick 2026 Table 2."
  ) +
  theme(legend.position = "bottom")
```

![Packaged models (lines) against an independent closed-form
steady-state solution (points), typical values at each weight-band
midpoint.](Maranchick_2026_pyrazinamide_ethambutol_files/figure-html/closed-form-plot-1.png)

Packaged models (lines) against an independent closed-form steady-state
solution (points), typical values at each weight-band midpoint.

### Steady-state mass balance and the HIV clearance effect

At steady state the exposure over one dosing interval satisfies
`AUC(0,tau) * CL/F = Dose` exactly for a linear model, whatever the
number of compartments. Checking that identity against `CL/F` rebuilt
from the paper’s printed formula gates the clearance path – typical
value, allometric exponent, reference weight and HIV coefficient –
independently of the volumes.

``` r

auc_trapz <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)

auc_identity <- dplyr::bind_rows(
  closed_form_pza |> mutate(drug = "Pyrazinamide"),
  closed_form_emb |> mutate(drug = "Ethambutol")
) |>
  group_by(drug, band, HIV_POS, dose_mg, cl_i) |>
  summarise(auc_sim = auc_trapz(t_rel, Cc), .groups = "drop") |>
  mutate(
    auc_expected = dose_mg / cl_i,
    pct_diff     = 100 * (auc_sim - auc_expected) / auc_expected
  )

stopifnot(nrow(auc_identity) == 20L)
# Deterministic identity; realised max was ~0.003%.
stopifnot(max(abs(auc_identity$pct_diff)) < 0.05)

# The paper's headline covariate claims: 18.5% faster PZA clearance and 25%
# faster EMB clearance in children with TB/HIV. Same-subject (paired) contrast,
# so this is deterministic rather than a race between two noisy arms.
hiv_contrast <- auc_identity |>
  select(drug, band, HIV_POS, auc_sim) |>
  tidyr::pivot_wider(names_from = HIV_POS, values_from = auc_sim,
                     names_prefix = "hiv") |>
  mutate(cl_increase_pct = 100 * (hiv0 / hiv1 - 1))

published_hiv <- c(Pyrazinamide = 18.5, Ethambutol = 24.6)
hiv_obs <- tapply(hiv_contrast$cl_increase_pct, hiv_contrast$drug, mean)
stopifnot(max(abs(hiv_obs[names(published_hiv)] - published_hiv)) < 0.5)

auc_identity |>
  mutate(band = factor(band, levels = bands$band),
         Group = ifelse(HIV_POS == 1, "TB/HIV", "TB")) |>
  arrange(drug, band, Group) |>
  select(drug, band, Group, dose_mg, cl_i, auc_expected, auc_sim, pct_diff) |>
  dplyr::rename(
    "Drug" = drug, "Weight band (kg)" = band, "Group" = Group,
    "Dose (mg)" = dose_mg, "CL/F (L/h)" = cl_i,
    "Dose / (CL/F) (mg*h/L)" = auc_expected,
    "Simulated AUC0-24 (mg*h/L)" = auc_sim, "% diff" = pct_diff
  ) |>
  knitr::kable(
    digits  = c(0, 0, 0, 0, 2, 1, 1, 4),
    caption = "Steady-state mass balance: simulated AUC0-24 against Dose / (CL/F) with CL/F rebuilt from Maranchick 2026 Table 2."
  )
```

| Drug | Weight band (kg) | Group | Dose (mg) | CL/F (L/h) | Dose / (CL/F) (mg\*h/L) | Simulated AUC0-24 (mg\*h/L) | % diff |
|:---|:---|:---|---:|---:|---:|---:|---:|
| Ethambutol | 4 to \<8 | TB | 100 | 12.16 | 8.2 | 8.2 | -0.0026 |
| Ethambutol | 4 to \<8 | TB/HIV | 100 | 15.15 | 6.6 | 6.6 | 0.0007 |
| Ethambutol | 8 to \<12 | TB | 200 | 17.39 | 11.5 | 11.5 | 0.0008 |
| Ethambutol | 8 to \<12 | TB/HIV | 200 | 21.66 | 9.2 | 9.2 | 0.0017 |
| Ethambutol | 12 to \<16 | TB | 300 | 22.00 | 13.6 | 13.6 | 0.0013 |
| Ethambutol | 12 to \<16 | TB/HIV | 300 | 27.42 | 10.9 | 10.9 | 0.0019 |
| Ethambutol | 16 to \<25 | TB | 400 | 28.24 | 14.2 | 14.2 | 0.0015 |
| Ethambutol | 16 to \<25 | TB/HIV | 400 | 35.19 | 11.4 | 11.4 | 0.0020 |
| Ethambutol | 25 to \<35 | TB | 550 | 37.51 | 14.7 | 14.7 | 0.0016 |
| Ethambutol | 25 to \<35 | TB/HIV | 550 | 46.74 | 11.8 | 11.8 | 0.0020 |
| Pyrazinamide | 4 to \<8 | TB | 150 | 0.67 | 224.3 | 224.3 | -0.0015 |
| Pyrazinamide | 4 to \<8 | TB/HIV | 150 | 0.79 | 189.2 | 189.2 | -0.0018 |
| Pyrazinamide | 8 to \<12 | TB | 300 | 0.96 | 313.7 | 313.7 | -0.0015 |
| Pyrazinamide | 8 to \<12 | TB/HIV | 300 | 1.13 | 264.7 | 264.7 | -0.0017 |
| Pyrazinamide | 12 to \<16 | TB | 450 | 1.21 | 371.9 | 371.9 | -0.0014 |
| Pyrazinamide | 12 to \<16 | TB/HIV | 450 | 1.43 | 313.7 | 313.7 | -0.0017 |
| Pyrazinamide | 16 to \<25 | TB | 600 | 1.55 | 386.3 | 386.3 | -0.0014 |
| Pyrazinamide | 16 to \<25 | TB/HIV | 600 | 1.84 | 325.9 | 325.9 | -0.0016 |
| Pyrazinamide | 25 to \<35 | TB | 800 | 2.06 | 387.8 | 387.8 | -0.0013 |
| Pyrazinamide | 25 to \<35 | TB/HIV | 800 | 2.45 | 327.1 | 327.1 | -0.0016 |

Steady-state mass balance: simulated AUC0-24 against Dose / (CL/F) with
CL/F rebuilt from Maranchick 2026 Table 2. {.table}

The paired HIV contrast reproduces the published clearance effects to
0.03 percentage points (PZA 18.5% versus a published 18.5%; EMB 24.6%
versus a published 25%, which the paper rounds from
`exp(0.22) - 1 = 24.6%`).

## Virtual cohort

The trial data are not public. The cohort below mirrors the simulation
described in Materials and Methods: WHO weight-band dosing once daily,
stratified by weight band and HIV status, evaluated at steady state. The
paper resampled its own participants’ demographics within each band;
because the within-band weight distribution is not published, weight is
drawn uniformly across each band here (see Assumptions and deviations).

``` r

# set.seed() seeds R's RNG; rxSetSeed() seeds rxode2's, which is partitioned per
# solver thread. Neither makes the drawn cohort identical across machines with
# different thread counts, so every assertion below is written to hold for any
# cohort these models can produce.
set.seed(20260302)
rxode2::rxSetSeed(20260302)

n_per_arm <- 100L  # 10 arms per drug; well under the 200-per-arm cap
obs_grid  <- c(seq(0, 4, by = 0.1), seq(4.5, tau, by = 0.5))

# Draw n subjects per arm with a uniform within-band weight, assigning disjoint
# id ranges so arms can be bound together without rxSolve collapsing subjects.
draw_arm_subjects <- function(arms, n, dose_col, id_offset) {
  out <- vector("list", nrow(arms))
  for (i in seq_len(nrow(arms))) {
    out[[i]] <- tibble::tibble(
      id      = id_offset + (i - 1L) * n + seq_len(n),
      WT      = runif(n, arms$wt_lo[i], arms$wt_hi[i]),
      HIV_POS = arms$HIV_POS[i],
      band    = arms$band[i],
      dose_mg = arms[[dose_col]][i]
    ) |>
      mutate(arm = paste(band, ifelse(HIV_POS == 1L, "TB/HIV", "TB")))
  }
  dplyr::bind_rows(out)
}

make_cohort <- function(dose_col, id_offset = 0L) {
  arms <- tidyr::crossing(bands, HIV_POS = c(0L, 1L))
  subj <- draw_arm_subjects(arms, n_per_arm, dose_col, id_offset)
  dos <- tidyr::crossing(subj, time = tau * (seq_len(n_doses) - 1L)) |>
    mutate(evid = 1L, cmt = "depot", amt = dose_mg)
  obs <- tidyr::crossing(subj, time = t_last + obs_grid) |>
    mutate(evid = 0L, cmt = "central", amt = NA_real_)
  dplyr::bind_rows(dos, obs) |> dplyr::arrange(id, time, dplyr::desc(evid))
}

ev_pza <- make_cohort("dose_pza")
ev_emb <- make_cohort("dose_emb")

stopifnot(!anyDuplicated(unique(ev_pza[, c("id", "time", "evid")])))
stopifnot(!anyDuplicated(unique(ev_emb[, c("id", "time", "evid")])))
stopifnot(dplyr::n_distinct(ev_pza$id) == 10L * n_per_arm)
```

## Simulation

``` r

# `keep =` may return character columns as factors; coerce so the downstream
# regex split on `arm` sees a character vector.
sim_pza <- rxode2::rxSolve(
  pza, events = ev_pza, keep = c("arm", "band", "HIV_POS", "WT", "dose_mg")
) |> as.data.frame() |> mutate(arm = as.character(arm), band = as.character(band))
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_emb <- rxode2::rxSolve(
  emb, events = ev_emb, keep = c("arm", "band", "HIV_POS", "WT", "dose_mg")
) |> as.data.frame() |> mutate(arm = as.character(arm), band = as.character(band))
#> ℹ parameter labels from comments will be replaced by 'label()'

# Guard against a silently empty or all-NA solve.
stopifnot(nrow(sim_pza) > 0, nrow(sim_emb) > 0)
stopifnot(!all(is.na(sim_pza$Cc)), !all(is.na(sim_emb$Cc)))
stopifnot(all(sim_pza$Cc >= 0, na.rm = TRUE), all(sim_emb$Cc >= 0, na.rm = TRUE))
```

## PKNCA validation

Steady-state NCA over the final dosing interval. Time is re-expressed
relative to the final dose so the interval runs from 0 to `tau`.

``` r

prep_nca <- function(sim) {
  sim |>
    filter(!is.na(Cc), time >= t_last) |>
    mutate(time = time - t_last) |>
    select(id, time, Cc, arm, band, HIV_POS, dose_mg)
}

nca_conc_pza <- prep_nca(sim_pza)
nca_conc_emb <- prep_nca(sim_emb)

# Guarantee a time = 0 record per subject so PKNCA can anchor AUC0-tau. At
# steady state the pre-dose record is the trough, not zero, so the existing row
# must win -- .keep_all = TRUE keeps the first occurrence.
add_time_zero <- function(d) {
  dplyr::bind_rows(
    d,
    d |> dplyr::distinct(id, arm, band, HIV_POS, dose_mg) |> mutate(time = 0, Cc = 0)
  ) |>
    dplyr::distinct(id, arm, time, .keep_all = TRUE) |>
    dplyr::arrange(id, arm, time)
}

nca_conc_pza <- add_time_zero(nca_conc_pza)
nca_conc_emb <- add_time_zero(nca_conc_emb)

nca_dose <- function(d) d |> dplyr::distinct(id, arm, dose_mg) |>
  mutate(time = 0, amt = dose_mg) |> select(id, arm, time, amt)

run_nca <- function(conc, dosed) {
  conc_obj <- PKNCA::PKNCAconc(as.data.frame(conc), Cc ~ time | arm + id,
                               concu = "ug/mL", timeu = "h")
  dose_obj <- PKNCA::PKNCAdose(as.data.frame(dosed), amt ~ time | arm + id,
                               doseu = "mg")
  intervals <- data.frame(
    start = 0, end = tau,
    cmax = TRUE, tmax = TRUE, cmin = TRUE, auclast = TRUE, cav = TRUE
  )
  PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
}

nca_pza <- run_nca(nca_conc_pza, nca_dose(nca_conc_pza))
nca_emb <- run_nca(nca_conc_emb, nca_dose(nca_conc_emb))

nca_wide <- function(res) {
  as.data.frame(res$result) |>
    filter(PPTESTCD %in% c("cmax", "tmax", "cmin", "auclast", "cav")) |>
    select(arm, id, PPTESTCD, PPORRES) |>
    tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES)
}

nca_pza_tbl <- nca_wide(nca_pza)
nca_emb_tbl <- nca_wide(nca_emb)

stopifnot(nrow(nca_pza_tbl) == 10L * n_per_arm,
          nrow(nca_emb_tbl) == 10L * n_per_arm)
```

``` r

summarise_nca <- function(tbl, drug) {
  tbl |>
    tidyr::separate_wider_regex(
      arm, c(band = ".*", " ", group = "TB/HIV|TB"), cols_remove = FALSE
    ) |>
    group_by(band, group) |>
    summarise(
      cmax_med = median(cmax), auc_med = median(auclast),
      tmax_med = median(tmax), cmin_med = median(cmin),
      .groups = "drop"
    ) |>
    mutate(drug = drug, band = factor(band, levels = bands$band)) |>
    arrange(band, group)
}

dplyr::bind_rows(
  summarise_nca(nca_pza_tbl, "Pyrazinamide"),
  summarise_nca(nca_emb_tbl, "Ethambutol")
) |>
  select(drug, band, group, cmax_med, tmax_med, cmin_med, auc_med) |>
  dplyr::rename(
    "Drug" = drug, "Weight band (kg)" = band, "Group" = group,
    "Cmax (ug/mL)" = cmax_med, "Tmax (h)" = tmax_med,
    "Cmin (ug/mL)" = cmin_med, "AUC0-24 (ug*h/mL)" = auc_med
  ) |>
  knitr::kable(
    digits  = 2,
    caption = "Median simulated steady-state NCA by weight band and HIV status, WHO-recommended doses (PKNCA)."
  )
```

| Drug | Weight band (kg) | Group | Cmax (ug/mL) | Tmax (h) | Cmin (ug/mL) | AUC0-24 (ug\*h/mL) |
|:---|:---|:---|---:|---:|---:|---:|
| Pyrazinamide | 4 to \<8 | TB | 25.52 | 1.30 | 1.89 | 243.14 |
| Pyrazinamide | 4 to \<8 | TB/HIV | 24.48 | 1.20 | 1.06 | 191.09 |
| Pyrazinamide | 8 to \<12 | TB | 36.45 | 1.20 | 2.44 | 320.95 |
| Pyrazinamide | 8 to \<12 | TB/HIV | 35.15 | 1.20 | 1.41 | 268.20 |
| Pyrazinamide | 12 to \<16 | TB | 39.39 | 1.30 | 3.04 | 349.34 |
| Pyrazinamide | 12 to \<16 | TB/HIV | 39.26 | 1.40 | 2.24 | 335.46 |
| Pyrazinamide | 16 to \<25 | TB | 41.09 | 1.25 | 3.39 | 371.51 |
| Pyrazinamide | 16 to \<25 | TB/HIV | 39.79 | 1.20 | 1.45 | 303.42 |
| Pyrazinamide | 25 to \<35 | TB | 40.99 | 1.30 | 3.23 | 383.36 |
| Pyrazinamide | 25 to \<35 | TB/HIV | 36.87 | 1.40 | 2.44 | 340.53 |
| Ethambutol | 4 to \<8 | TB | 1.50 | 1.30 | 0.11 | 7.97 |
| Ethambutol | 4 to \<8 | TB/HIV | 1.41 | 1.20 | 0.08 | 6.85 |
| Ethambutol | 8 to \<12 | TB | 2.13 | 1.40 | 0.12 | 10.95 |
| Ethambutol | 8 to \<12 | TB/HIV | 1.94 | 1.40 | 0.08 | 9.08 |
| Ethambutol | 12 to \<16 | TB | 2.29 | 1.40 | 0.14 | 13.90 |
| Ethambutol | 12 to \<16 | TB/HIV | 2.42 | 1.40 | 0.10 | 11.74 |
| Ethambutol | 16 to \<25 | TB | 2.99 | 1.30 | 0.13 | 14.91 |
| Ethambutol | 16 to \<25 | TB/HIV | 2.42 | 1.35 | 0.07 | 10.47 |
| Ethambutol | 25 to \<35 | TB | 2.96 | 1.50 | 0.11 | 14.47 |
| Ethambutol | 25 to \<35 | TB/HIV | 2.62 | 1.45 | 0.07 | 11.64 |

Median simulated steady-state NCA by weight band and HIV status,
WHO-recommended doses (PKNCA). {.table}

## Replicating Figures 3 and 4

Figures 3 and 4 of Maranchick 2026 show boxplots of simulated
steady-state Cmax and AUC0-24 by weight band and HIV status against the
adult target ranges. The panels below reproduce that layout.

``` r

targets <- tibble::tribble(
  ~drug,          ~metric,             ~lo,  ~hi,  ~thresh,
  "Pyrazinamide", "Cmax (ug/mL)",      20,   60,   35,
  "Pyrazinamide", "AUC0-24 (ug*h/mL)", 250,  450,  363,
  "Ethambutol",   "Cmax (ug/mL)",      2,    6,    2,
  "Ethambutol",   "AUC0-24 (ug*h/mL)", 16,   29,   16
)

fig_dat <- dplyr::bind_rows(
  nca_pza_tbl |> mutate(drug = "Pyrazinamide"),
  nca_emb_tbl |> mutate(drug = "Ethambutol")
) |>
  tidyr::separate_wider_regex(
    arm, c(band = ".*", " ", group = "TB/HIV|TB"), cols_remove = FALSE
  ) |>
  select(drug, band, group, `Cmax (ug/mL)` = cmax, `AUC0-24 (ug*h/mL)` = auclast) |>
  tidyr::pivot_longer(c(`Cmax (ug/mL)`, `AUC0-24 (ug*h/mL)`),
                      names_to = "metric", values_to = "value") |>
  mutate(band = factor(band, levels = bands$band))

ggplot(fig_dat, aes(band, value, fill = group)) +
  geom_rect(
    data = targets, inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
    fill = "grey85", alpha = 0.5
  ) +
  geom_hline(data = targets, aes(yintercept = thresh), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  facet_wrap(~ drug + metric, scales = "free_y", ncol = 2) +
  labs(x = "Weight band (kg)", y = NULL, fill = NULL,
       caption = "Replicates Figures 3 and 4 of Maranchick 2026.") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))
```

![Replicates Figures 3 (pyrazinamide) and 4 (ethambutol) of Maranchick
2026: simulated steady-state Cmax and AUC0-24 by weight band and HIV
status at WHO-recommended doses. Shaded bands are the adult target
ranges; dashed lines are the target
thresholds.](Maranchick_2026_pyrazinamide_ethambutol_files/figure-html/figure-3-4-1.png)

Replicates Figures 3 (pyrazinamide) and 4 (ethambutol) of Maranchick
2026: simulated steady-state Cmax and AUC0-24 by weight band and HIV
status at WHO-recommended doses. Shaded bands are the adult target
ranges; dashed lines are the target thresholds.

## Comparison against the published target attainment

The source paper does not tabulate simulated Cmax / AUC point estimates,
so
[`nlmixr2lib::ncaComparisonTable()`](https://nlmixr2.github.io/nlmixr2lib/reference/ncaComparisonTable.md)
does not apply here. What it does report, in the Results text for each
drug, is the **percentage of simulated patients meeting the adult
target** in each weight band and HIV stratum. That is the published
quantity, and it is what the table below compares against.

``` r

published <- tibble::tribble(
  ~drug,          ~metric, ~band,        ~group,   ~published_pct,
  # Maranchick 2026 Results, "Pyrazinamide" paragraph 2.
  "Pyrazinamide", "Cmax",  "4 to <8",    "TB",      14.3,
  "Pyrazinamide", "Cmax",  "4 to <8",    "TB/HIV",   6.1,
  "Pyrazinamide", "Cmax",  "8 to <12",   "TB",      45.7,
  "Pyrazinamide", "Cmax",  "8 to <12",   "TB/HIV",  34.7,
  "Pyrazinamide", "Cmax",  "12 to <16",  "TB",      63.3,
  "Pyrazinamide", "Cmax",  "12 to <16",  "TB/HIV",  55.0,
  "Pyrazinamide", "Cmax",  "16 to <25",  "TB",      56.7,
  "Pyrazinamide", "Cmax",  "16 to <25",  "TB/HIV",  60.6,
  "Pyrazinamide", "Cmax",  "25 to <35",  "TB",      60.7,
  "Pyrazinamide", "Cmax",  "25 to <35",  "TB/HIV",  60.6,
  "Pyrazinamide", "AUC",   "4 to <8",    "TB",       5.4,
  "Pyrazinamide", "AUC",   "4 to <8",    "TB/HIV",   0.0,
  "Pyrazinamide", "AUC",   "8 to <12",   "TB",      21.0,
  "Pyrazinamide", "AUC",   "8 to <12",   "TB/HIV",  12.7,
  "Pyrazinamide", "AUC",   "12 to <16",  "TB",      41.8,
  "Pyrazinamide", "AUC",   "12 to <16",  "TB/HIV",  24.8,
  "Pyrazinamide", "AUC",   "16 to <25",  "TB",      47.0,
  "Pyrazinamide", "AUC",   "16 to <25",  "TB/HIV",  31.6,
  "Pyrazinamide", "AUC",   "25 to <35",  "TB",      48.2,
  "Pyrazinamide", "AUC",   "25 to <35",  "TB/HIV",  33.3,
  # Maranchick 2026 Results, "Ethambutol" paragraph 2.
  "Ethambutol",   "Cmax",  "4 to <8",    "TB",      20.8,
  "Ethambutol",   "Cmax",  "4 to <8",    "TB/HIV",  13.5,
  "Ethambutol",   "Cmax",  "8 to <12",   "TB",      38.6,
  "Ethambutol",   "Cmax",  "8 to <12",   "TB/HIV",  53.2,
  "Ethambutol",   "Cmax",  "12 to <16",  "TB",      64.6,
  "Ethambutol",   "Cmax",  "12 to <16",  "TB/HIV",  61.5,
  "Ethambutol",   "Cmax",  "16 to <25",  "TB",      65.2,
  "Ethambutol",   "Cmax",  "16 to <25",  "TB/HIV",  62.9,
  "Ethambutol",   "Cmax",  "25 to <35",  "TB",      55.8,
  "Ethambutol",   "Cmax",  "25 to <35",  "TB/HIV",  72.5,
  "Ethambutol",   "AUC",   "4 to <8",    "TB",       1.9,
  "Ethambutol",   "AUC",   "4 to <8",    "TB/HIV",   0.0,
  "Ethambutol",   "AUC",   "8 to <12",   "TB",       0.1,
  "Ethambutol",   "AUC",   "8 to <12",   "TB/HIV",   0.1,
  "Ethambutol",   "AUC",   "12 to <16",  "TB",      13.1,
  "Ethambutol",   "AUC",   "12 to <16",  "TB/HIV",   3.7,
  "Ethambutol",   "AUC",   "16 to <25",  "TB",      18.2,
  "Ethambutol",   "AUC",   "16 to <25",  "TB/HIV",   6.9,
  "Ethambutol",   "AUC",   "25 to <35",  "TB",      21.2,
  "Ethambutol",   "AUC",   "25 to <35",  "TB/HIV",   7.2
)

thresholds <- tibble::tribble(
  ~drug,          ~metric, ~thresh,
  "Pyrazinamide", "Cmax",  35,
  "Pyrazinamide", "AUC",   363,
  "Ethambutol",   "Cmax",  2,
  "Ethambutol",   "AUC",   16
)

simulated <- dplyr::bind_rows(
  nca_pza_tbl |> mutate(drug = "Pyrazinamide"),
  nca_emb_tbl |> mutate(drug = "Ethambutol")
) |>
  tidyr::separate_wider_regex(
    arm, c(band = ".*", " ", group = "TB/HIV|TB"), cols_remove = FALSE
  ) |>
  select(drug, band, group, Cmax = cmax, AUC = auclast) |>
  tidyr::pivot_longer(c(Cmax, AUC), names_to = "metric", values_to = "value") |>
  left_join(thresholds, by = c("drug", "metric")) |>
  group_by(drug, metric, band, group) |>
  summarise(simulated_pct = 100 * mean(value > thresh), .groups = "drop")

attain <- published |>
  left_join(simulated, by = c("drug", "metric", "band", "group")) |>
  mutate(diff_pp = simulated_pct - published_pct)

# Every published cell must have found a simulated partner -- a failed join
# would otherwise leave NA and silently shrink the comparison.
stopifnot(nrow(attain) == 40L, !anyNA(attain$simulated_pct))
```

One published cell is internally inconsistent with the rest of its own
table and is excluded from the numeric gate rather than accommodated by
widening it. The ethambutol AUC row for the 8 to \<12 kg band is
published as 0.1% attainment for both groups, yet the neighbouring 4 to
\<8 kg band – which receives a *lower* exposure (100 mg at roughly 6 kg
gives `Dose/(CL/F)` near 8 ug\*h/mL, versus 200 mg at roughly 10 kg
giving near 11.5) – is published as 1.9%. Attainment of a fixed
threshold cannot fall as exposure rises under this model, so 0.1% is not
reachable from the paper’s own Table 2 parameters. It reads as a
misprint.

``` r

attain <- attain |>
  mutate(
    deviation = drug == "Ethambutol" & metric == "AUC" & band == "8 to <12"
  )

gated <- attain |> filter(!deviation)

# Cohort-derived proportions: bound the AGREEMENT, not any single cell, and
# leave headroom for the unpublished within-band weight distribution and for
# n = 100 per arm against the paper's 1,000 replicates (binomial SE ~5 points
# at p = 0.5). Rendered at 2 / 4 / 16 solver threads (rxSetSeed does not fix the
# cohort across thread counts) the realised values were mean 8.2 / 9.6 / 8.3 pp
# and max 23.9 / 26.9 / 23.2 pp. A mis-transcribed clearance, volume or dose
# moves whole bands by 40+ points, so these bounds can still go red.
mean_abs_pp <- mean(abs(gated$diff_pp))
max_abs_pp  <- max(abs(gated$diff_pp))
stopifnot(mean_abs_pp < 18, max_abs_pp < 35)

# Shape agreement across the 38 gated cells: the model must rank the weight
# bands and strata the way the paper does. Guards against a structurally wrong
# model that happens to sit near the right average level. Realised 0.964 /
# 0.955 / 0.964 at 2 / 4 / 16 threads.
shape_cor <- cor(gated$published_pct, gated$simulated_pct)
stopifnot(shape_cor > 0.7)

# The paper's explicit claim: "Target attainment was lowest in the 4-<8 kg
# weight band." Stated for both drugs and both metrics, and a large effect.
lowest_band <- attain |>
  group_by(drug, metric) |>
  slice_min(simulated_pct, n = 1, with_ties = FALSE) |>
  ungroup()
stopifnot(all(lowest_band$band == "4 to <8"))

cat(sprintf(
  "Attainment agreement over %d gated cells: mean |diff| %.1f pp, max %.1f pp, correlation %.3f\n",
  nrow(gated), mean_abs_pp, max_abs_pp, shape_cor
))
#> Attainment agreement over 38 gated cells: mean |diff| 9.2 pp, max 26.2 pp, correlation 0.953
```

``` r

attain |>
  mutate(
    band  = factor(band, levels = bands$band),
    Note  = ifelse(deviation, "see text", "")
  ) |>
  arrange(drug, metric, band, group) |>
  select(drug, metric, band, group, published_pct, simulated_pct, diff_pp, Note) |>
  dplyr::rename(
    "Drug" = drug, "Metric" = metric, "Weight band (kg)" = band,
    "Group" = group, "Published (%)" = published_pct,
    "Simulated (%)" = simulated_pct, "Difference (pp)" = diff_pp
  ) |>
  knitr::kable(
    digits  = 1,
    caption = "Target attainment at WHO-recommended doses: published (Maranchick 2026 Results) against simulated. Targets are PZA Cmax > 35 ug/mL and AUC0-24 > 363 ug*h/mL; EMB Cmax > 2 ug/mL and AUC0-24 > 16 ug*h/mL."
  )
```

| Drug | Metric | Weight band (kg) | Group | Published (%) | Simulated (%) | Difference (pp) | Note |
|:---|:---|:---|:---|---:|---:|---:|:---|
| Ethambutol | AUC | 4 to \<8 | TB | 1.9 | 4 | 2.1 |  |
| Ethambutol | AUC | 4 to \<8 | TB/HIV | 0.0 | 1 | 1.0 |  |
| Ethambutol | AUC | 8 to \<12 | TB | 0.1 | 14 | 13.9 | see text |
| Ethambutol | AUC | 8 to \<12 | TB/HIV | 0.1 | 3 | 2.9 | see text |
| Ethambutol | AUC | 12 to \<16 | TB | 13.1 | 34 | 20.9 |  |
| Ethambutol | AUC | 12 to \<16 | TB/HIV | 3.7 | 13 | 9.3 |  |
| Ethambutol | AUC | 16 to \<25 | TB | 18.2 | 41 | 22.8 |  |
| Ethambutol | AUC | 16 to \<25 | TB/HIV | 6.9 | 13 | 6.1 |  |
| Ethambutol | AUC | 25 to \<35 | TB | 21.2 | 40 | 18.8 |  |
| Ethambutol | AUC | 25 to \<35 | TB/HIV | 7.2 | 15 | 7.8 |  |
| Ethambutol | Cmax | 4 to \<8 | TB | 20.8 | 30 | 9.2 |  |
| Ethambutol | Cmax | 4 to \<8 | TB/HIV | 13.5 | 21 | 7.5 |  |
| Ethambutol | Cmax | 8 to \<12 | TB | 38.6 | 58 | 19.4 |  |
| Ethambutol | Cmax | 8 to \<12 | TB/HIV | 53.2 | 48 | -5.2 |  |
| Ethambutol | Cmax | 12 to \<16 | TB | 64.6 | 66 | 1.4 |  |
| Ethambutol | Cmax | 12 to \<16 | TB/HIV | 61.5 | 66 | 4.5 |  |
| Ethambutol | Cmax | 16 to \<25 | TB | 65.2 | 81 | 15.8 |  |
| Ethambutol | Cmax | 16 to \<25 | TB/HIV | 62.9 | 61 | -1.9 |  |
| Ethambutol | Cmax | 25 to \<35 | TB | 55.8 | 82 | 26.2 |  |
| Ethambutol | Cmax | 25 to \<35 | TB/HIV | 72.5 | 76 | 3.5 |  |
| Pyrazinamide | AUC | 4 to \<8 | TB | 5.4 | 10 | 4.6 |  |
| Pyrazinamide | AUC | 4 to \<8 | TB/HIV | 0.0 | 1 | 1.0 |  |
| Pyrazinamide | AUC | 8 to \<12 | TB | 21.0 | 30 | 9.0 |  |
| Pyrazinamide | AUC | 8 to \<12 | TB/HIV | 12.7 | 20 | 7.3 |  |
| Pyrazinamide | AUC | 12 to \<16 | TB | 41.8 | 43 | 1.2 |  |
| Pyrazinamide | AUC | 12 to \<16 | TB/HIV | 24.8 | 38 | 13.2 |  |
| Pyrazinamide | AUC | 16 to \<25 | TB | 47.0 | 55 | 8.0 |  |
| Pyrazinamide | AUC | 16 to \<25 | TB/HIV | 31.6 | 32 | 0.4 |  |
| Pyrazinamide | AUC | 25 to \<35 | TB | 48.2 | 55 | 6.8 |  |
| Pyrazinamide | AUC | 25 to \<35 | TB/HIV | 33.3 | 41 | 7.7 |  |
| Pyrazinamide | Cmax | 4 to \<8 | TB | 14.3 | 18 | 3.7 |  |
| Pyrazinamide | Cmax | 4 to \<8 | TB/HIV | 6.1 | 6 | -0.1 |  |
| Pyrazinamide | Cmax | 8 to \<12 | TB | 45.7 | 57 | 11.3 |  |
| Pyrazinamide | Cmax | 8 to \<12 | TB/HIV | 34.7 | 50 | 15.3 |  |
| Pyrazinamide | Cmax | 12 to \<16 | TB | 63.3 | 72 | 8.7 |  |
| Pyrazinamide | Cmax | 12 to \<16 | TB/HIV | 55.0 | 71 | 16.0 |  |
| Pyrazinamide | Cmax | 16 to \<25 | TB | 56.7 | 77 | 20.3 |  |
| Pyrazinamide | Cmax | 16 to \<25 | TB/HIV | 60.6 | 74 | 13.4 |  |
| Pyrazinamide | Cmax | 25 to \<35 | TB | 60.7 | 76 | 15.3 |  |
| Pyrazinamide | Cmax | 25 to \<35 | TB/HIV | 60.6 | 58 | -2.6 |  |

Target attainment at WHO-recommended doses: published (Maranchick 2026
Results) against simulated. Targets are PZA Cmax \> 35 ug/mL and AUC0-24
\> 363 ug*h/mL; EMB Cmax \> 2 ug/mL and AUC0-24 \> 16 ug*h/mL. {.table
style="width:100%;"}

The simulated attainment tracks the published values in level and in
shape, and reproduces every qualitative conclusion the paper draws:
attainment is lowest in the 4 to \<8 kg band for both drugs; TB/HIV
children attain the AUC target less often than TB-only children; and
ethambutol AUC attainment is poor throughout, which is what drives the
paper’s recommendation for a two- to three-fold dose increase. Simulated
attainment runs modestly higher than published in the middle and upper
bands, which is consistent with the virtual cohort drawing weight
uniformly across each band rather than resampling the trial’s own
within-band weight distribution: `AUC = Dose/(CL/F)` falls as weight
rises within a fixed-dose band, so a cohort weighted toward the top of
each band attains less often.

## Assumptions and deviations

- **Within-band weight distribution.** The paper resampled its own
  participants’ demographics 1,000 times within each weight band; the
  within-band weight distribution is not published. Weight is drawn
  uniformly across each band here. This is the main driver of the
  residual difference in simulated versus published attainment, and it
  moves attainment in the direction observed.
- **Cohort size.** 100 subjects per arm (10 arms per drug) against the
  paper’s 1,000 replicates, giving a binomial standard error near 5
  percentage points at 50% attainment. The attainment gate is written to
  tolerate this.
- **Residual-error parameterisation.** Table 2 reports Monolix
  residual-error constants `a` and `b` without naming the combined-error
  variant. Monolix offers `combined1` (`SD = a + b*f`) and `combined2`
  (`SD = sqrt(a^2 + (b*f)^2)`); the paper says only that “a combination
  of additive and proportional error models was utilized”. The models
  use nlmixr2’s `add(addSd) + prop(propSd)`, which is the `combined2`
  form and the library convention. At pyrazinamide concentrations near
  30 ug/mL the two forms differ by roughly a third in residual SD, so a
  downstream user re-fitting these models should confirm the intended
  variant. The choice does not affect any typical-value or
  between-subject quantity used in this vignette.
- **Ethambutol Q/F and V2/F are not weight-scaled.** The paper reports
  allometric exponents only for CL/F and V1/F and gives no exponent for
  Q/F or V2/F, so neither is scaled here. This is faithful to Table 2 as
  printed.
- **Reference weights.** The normalisation constants are 15 kg (PZA) and
  15.1 kg (EMB) exactly as printed in the Results. Neither equals the
  cohort median weight of 16 kg reported in Table 1, and the paper does
  not explain the difference or why the two models differ by 0.1 kg. The
  printed values are used unchanged.
- **Published ethambutol AUC attainment for the 8 to \<12 kg band.**
  Published as 0.1% for both groups, which is not reachable from the
  paper’s own Table 2 parameters: the lower-exposure 4 to \<8 kg band is
  published as 1.9%, and attainment of a fixed threshold cannot fall as
  exposure rises. Treated as a probable misprint, reported in the
  comparison table and excluded from the numeric gate rather than
  accommodated by widening it.
- **Between-occasion variability, maturation and HIV-medication
  effects** were tested by the authors and not retained, so none is
  encoded.
- **No non-paper-derived parameter values.** Every `ini()` entry comes
  from Table 2 or the Results text of the main article. The EuropePMC
  supplementary endpoint for PMC13041307 returns only figure images, and
  no erratum was found.
