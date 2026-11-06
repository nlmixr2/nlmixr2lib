# Anifrolumab, sifalimumab, daxdilimab and litifilimab type I IFN QSP in SLE (Volkova 2025)

## Model and source

``` r

mod <- readModelDb("Volkova_2025_SLE_ifn_qsp")
ui <- rxode2::rxode(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_typ <- rxode2::zeroRe(ui)
#> Warning: No sigma parameters in the model
```

- Citation: Volkova A, Sokolov V, Tettamanti F, Verma M, Ugolkov Y,
  Peskov K, Tang W, Kimko H (2025). An Integrative Mechanistic Model of
  Type 1 IFN-Mediated Inflammation in Systemic Lupus Erythematosus. CPT
  Pharmacometrics Syst Pharmacol 14(7):1225-1236.
  <doi:10.1002/psp4.70043>. PMID: 40364448; PMCID: PMC12256560.
  Equations from Data S3 (Supplementary Files 1-2, Equations 1-25);
  parameter values from Table S2, the Monolix model code in Data S2 and
  the typical baseline IFN-alpha from Table S3.
- Article: <https://doi.org/10.1002/psp4.70043> (open access,
  PMC12256560)

QSP. Volkova 2025 quantitative systems pharmacology model of type I
interferon (IFN1)-mediated inflammation in adults with systemic lupus
erythematosus (SLE), coupled to the PK of four IFN1-pathway-targeting
monoclonal antibodies. Twenty ODEs. The biological core (7 ODEs) carries
IFN1-producing immune cells (IC, proportional to plasmacytoid dendritic
cells), free IFN-alpha and IFN-beta, free IFNAR1 receptor, the
IFNAR1-IFN-alpha and IFNAR1-IFN-beta complexes and the 21-gene IFN gene
signature (IFNGS), which is stimulated log-linearly by the summed
IFN1-IFNAR1 complexes above a residual IFN1-independent floor. Two
feedbacks: IFNGS relative to baseline amplifies IFN1 production per IC
and accelerates IC loss from blood. Drugs: anifrolumab (anti-IFNAR1;
two-compartment PK with linear plus saturable elimination, IFNAR1
binding with complex internalisation, IV or SC), sifalimumab
(anti-IFN-alpha; two-compartment PK with IFN-alpha binding, IV or SC),
daxdilimab (anti-ILT7, pDC-depleting; one-compartment SC PK, Emax
stimulation of IC loss) and litifilimab (anti-BDCA2; two-compartment SC
PK with lag and Michaelis-Menten elimination, Hill inhibition of IFN1
production). The system starts at its steady-state baseline. IIV on the
IFNGS synthesis rate constant and on anifrolumab bioavailability, linear
elimination and Vmax; residual error magnitudes are not reported, so
none is encoded. Doses are in mg for all four drugs.

The model is a quantitative systems pharmacology (QSP) model. Its
biological core is an endogenous steady-state system, so the validation
below combines the endogenous pattern (steady-state hold, closed-form
limits) with replication of the paper’s reported drug effects and a
PKNCA summary of the simulated anifrolumab exposure.

## Population

The systemic IFN1 / IFN gene signature (IFNGS) parameters and the
anifrolumab PK sub-model were calibrated on patient-level data from the
Phase IIb MUSE trial (NCT01438489; 254 adults with moderate-to-severe
SLE; anifrolumab 300 or 1000 mg IV every 4 weeks, or placebo).
Daxdilimab PK and its effect on IFN1-producing immune cells were fitted
to aggregated Phase I data (NCT03817424, 5-150 mg SC). Sifalimumab and
litifilimab PK and drug-effect parameters were taken from the published
population PK(/PD) models of Zheng 2016 and Hartmann 2020. External
validation used aggregated data from five anifrolumab, four sifalimumab,
one daxdilimab and one litifilimab trial (Table S1). The paper reports
no demographics.

``` r

str(ui$population)
#> List of 6
#>  $ species      : chr "human"
#>  $ n_subjects   : int 254
#>  $ n_studies    : int 11
#>  $ disease_state: chr "Adults with moderate-to-severe systemic lupus erythematosus (the litifilimab and daxdilimab validation trials a"| __truncated__
#>  $ dose_range   : chr "Anifrolumab 100-1000 mg IV Q4W and 150-300 mg SC Q2W; sifalimumab up to 30 mg/kg IV and SC regimens; daxdilimab"| __truncated__
#>  $ notes        : chr "Patient-level calibration data: 254 patients of the Phase IIb MUSE anifrolumab trial (NCT01438489; anifrolumab "| __truncated__
```

The only covariate is the baseline 21-gene IFNGS status, `BGENE21_HIGH`
(1 = baseline IFNGS at or above 2.75, chosen from the bimodal MUSE
distribution). IFNGS-low patients have an IFNGS synthesis rate constant
exp(-1.57) = 4.8-fold lower than IFNGS-high patients.

## Source trace

### Equations

| Model element | Source |
|----|----|
| Immune cells `ic`, IFNGS feedback and daxdilimab stimulation of loss | Data S3 Equation 1; Data S2 code `ddt_pDC` |
| Free IFN-alpha `ifna` with IFNGS feedback, litifilimab inhibition and sifalimumab binding | Equation 2; code `ddt_IFNa` |
| Free IFN-beta `ifnb` | Equation 3; code `ddt_IFNb` |
| Free IFNAR1 `ifnar1` with anifrolumab binding | Equation 4; code `ddt_IFNAR1` |
| IFNAR1-IFN-alpha / IFN-beta complexes | Equations 5-6 |
| IFNGS, log-linear in the summed complexes above the floor `ifngs_ll` | Equation 7; code `ddt_IFNGS` (amount scaled by 10^8) |
| Anifrolumab depot / central / peripheral with linear + saturable elimination | Equations 8-10 |
| Anifrolumab-IFNAR1 complex | Equation 11 |
| Sifalimumab depot / central / peripheral and sifalimumab-IFN-alpha complex | Equations 12-15 |
| Daxdilimab depot / central | Equations 16-17 |
| Litifilimab depot / central / peripheral with saturable elimination and lag | Equations 18-20; code `Tlag = Lag_Lit` |
| Steady-state synthesis rates and baseline complexes / IFNGS | Equations 21-25; code “SS for IFN in plasma” block |
| IFN-beta baseline = IFN-alpha baseline x (1 - 0.99) / 0.99 | Supplementary File 2; code `IFNb_BL_pgmL` |
| IFNGS-low effect `exp(beta * GRNUM)` on kin_IFNGS | code `kin_IFNGS_t`; Supplementary File 3 |
| Dose conversions (mg to nmol / ug) and bioavailabilities | code `iv()` statements |

### Parameters

``` r

ini_df <- as.data.frame(ui$iniDf)
ini_df |>
  dplyr::filter(is.na(neta1)) |>
  dplyr::transmute(
    Parameter = name,
    Value = signif(ifelse(grepl("^l(?!ogit)", name, perl = TRUE), exp(est), est), 4),
    Label = label,
    Fixed = fix
  ) |>
  knitr::kable(caption = "ini() values (log-scale entries shown back-transformed). Every value is traced line by line to Table S2, the Data S2 Monolix code or Table S3 in the model file comments.")
```

| Parameter | Value | Label | Fixed |
|:---|---:|:---|:---|
| bl_ifna | 6.350e+00 | Baseline plasma IFN-alpha (pg/mL) | TRUE |
| ifna_fraction | 9.900e-01 | Fraction of baseline IFN1 that is IFN-alpha (unitless) | TRUE |
| mw_ifna | 1.800e+04 | IFN-alpha molecular weight (g/mol) | TRUE |
| mw_ifnb | 2.200e+04 | IFN-beta molecular weight (g/mol) | TRUE |
| bl_ifnar1 | 2.600e-01 | Baseline free IFNAR1 concentration (nmol/L) | TRUE |
| kout_ifnar1 | 2.420e+01 | IFNAR1 elimination rate constant (1/day) | TRUE |
| kout_ifna_ifnar1 | 1.660e+01 | IFNAR1-IFN-alpha complex elimination rate constant (1/day) | TRUE |
| kout_ifnb_ifnar1 | 1.660e+01 | IFNAR1-IFN-beta complex elimination rate constant (1/day) | TRUE |
| kon_ifna | 4.320e+01 | IFN-alpha/IFNAR1 association rate constant (1/(nmol/L\*day)) | TRUE |
| kon_ifnb | 4.320e+01 | IFN-beta/IFNAR1 association rate constant (1/(nmol/L\*day)) | TRUE |
| kd_ifna | 5.000e+02 | IFN-alpha/IFNAR1 equilibrium dissociation constant (nmol/L) | TRUE |
| kd_ifnb | 1.000e+02 | IFN-beta/IFNAR1 equilibrium dissociation constant (nmol/L) | TRUE |
| kout_ifna | 8.300e+00 | Free IFN-alpha elimination rate constant (1/day) | TRUE |
| kout_ifnb | 4.160e+00 | Free IFN-beta elimination rate constant (1/day) | TRUE |
| lkin_ifngs | 2.700e+00 | IFNGS synthesis rate constant, IFNGS-high patients (1/day) | FALSE |
| e_bgene21_low_kin_ifngs | -1.570e+00 | IFNGS-low status effect on log kin_IFNGS (unitless) | FALSE |
| kout_ifngs | 1.920e+00 | IFNGS elimination rate constant (1/day) | TRUE |
| kn | 1.000e+00 | Slope of IFNGS stimulation by IFN1-IFNAR1 complexes (1/(nmol\*1e8)) | TRUE |
| ifngs_ll | 7.640e-01 | Residual IFNGS in the absence of IFN1 signalling (unitless) | TRUE |
| lvc | 2.750e+00 | Plasma volume shared by all four drugs and the IFN1 species (L) | TRUE |
| bl_ic | 1.000e+07 | Baseline IFN1-producing immune cells (cell) | TRUE |
| kout_ic | 1.000e-01 | Immune cell elimination rate constant (1/day) | TRUE |
| lka_anifrolumab | 2.000e-01 | Anifrolumab SC absorption rate constant (1/day) | TRUE |
| lk12_anifrolumab | 4.690e+00 | Anifrolumab central-to-peripheral rate constant (1/day) | FALSE |
| lk21_anifrolumab | 3.980e+00 | Anifrolumab peripheral-to-central rate constant (1/day) | FALSE |
| logitfbio_anifrolumab | 1.325e+00 | Logit of anifrolumab IV bioavailability (unitless) | FALSE |
| fsc_anifrolumab | 8.700e-01 | Anifrolumab SC bioavailability relative to IV (unitless) | TRUE |
| lkel_anifrolumab | 4.860e-02 | Anifrolumab linear elimination rate constant (1/day) | FALSE |
| lvm_anifrolumab | 2.140e+01 | Anifrolumab maximum saturable elimination rate (nmol/day) | FALSE |
| lkm_anifrolumab | 3.890e+01 | Anifrolumab amount at half-maximal saturable elimination (nmol) | FALSE |
| mw_anifrolumab | 1.451e+05 | Anifrolumab molecular weight (g/mol) | TRUE |
| kd_anifrolumab | 6.700e-02 | Anifrolumab/IFNAR1 equilibrium dissociation constant (nmol/L) | TRUE |
| kon_anifrolumab | 2.330e+02 | Anifrolumab/IFNAR1 association rate constant (1/(nmol/L\*day)) | TRUE |
| kint_anifrolumab | 7.710e+01 | Anifrolumab-IFNAR1 complex internalisation rate constant (1/day) | TRUE |
| lka_sifalimumab | 2.000e-01 | Sifalimumab SC absorption rate constant (1/day) | TRUE |
| lfdepot_sifalimumab | 7.000e-01 | Sifalimumab SC bioavailability (unitless) | TRUE |
| lcl_sifalimumab | 1.900e-01 | Sifalimumab linear clearance (L/day) | TRUE |
| lvp_sifalimumab | 1.890e+00 | Sifalimumab peripheral volume (L) | TRUE |
| lq_sifalimumab | 3.200e-01 | Sifalimumab intercompartmental clearance (L/day) | TRUE |
| kd_sifalimumab | 1.800e-01 | Sifalimumab/IFN-alpha equilibrium dissociation constant (nmol/L) | TRUE |
| koff_sifalimumab | 1.123e+01 | Sifalimumab/IFN-alpha dissociation rate constant (1/day) | TRUE |
| mw_sifalimumab | 1.470e+05 | Sifalimumab molecular weight (g/mol) | TRUE |
| kout_sifalimumab_ifna | 3.300e-02 | Sifalimumab-IFN-alpha complex elimination rate constant (1/day) | TRUE |
| lka_daxdilimab | 5.730e-01 | Daxdilimab SC absorption rate constant (1/day) | FALSE |
| lfdepot_daxdilimab | 1.470e-01 | Daxdilimab SC bioavailability (unitless) | FALSE |
| lkel_daxdilimab | 4.850e-02 | Daxdilimab elimination rate constant (1/day) | FALSE |
| mw_daxdilimab | 1.450e+05 | Daxdilimab molecular weight (g/mol) | TRUE |
| lec50_daxdilimab | 5.200e-01 | Daxdilimab concentration at half-maximal stimulation of IC loss (nmol/L) | FALSE |
| emax_daxdilimab | 4.000e+00 | Daxdilimab maximal stimulation of IC loss (unitless) | TRUE |
| lka_litifilimab | 2.340e-01 | Litifilimab SC absorption rate constant (1/day) | TRUE |
| lfdepot_litifilimab | 7.130e-01 | Litifilimab SC bioavailability (unitless) | TRUE |
| tlag_litifilimab | 9.000e-02 | Litifilimab absorption lag time (day) | TRUE |
| lkel_litifilimab | 7.100e-02 | Litifilimab linear elimination rate constant (1/day) | TRUE |
| lvp_litifilimab | 2.310e+00 | Litifilimab peripheral volume (L) | TRUE |
| lq_litifilimab | 4.896e-01 | Litifilimab intercompartmental clearance (L/day) | TRUE |
| lkm_litifilimab | 5.470e+02 | Litifilimab concentration at half-maximal saturable elimination (ug/L) | TRUE |
| lvm_litifilimab | 2.568e+02 | Litifilimab maximal saturable elimination rate (ug/day) | TRUE |
| mw_litifilimab | 1.460e+05 | Litifilimab molecular weight (g/mol) | TRUE |
| ic50_litifilimab | 3.500e+02 | Litifilimab concentration at half-maximal inhibition of IFN1 production (ug/L) | TRUE |
| imax_litifilimab | 8.750e-01 | Litifilimab maximal inhibition of IFN1 production (fraction) | TRUE |
| hill_litifilimab | 3.630e+00 | Litifilimab Hill coefficient (unitless) | TRUE |

ini() values (log-scale entries shown back-transformed). Every value is
traced line by line to Table S2, the Data S2 Monolix code or Table S3 in
the model file comments. {.table}

IIV (Table S2 CV, converted to log-normal variance as
`omega^2 = log(1 + CV^2)`; anifrolumab bioavailability is logit-normal
with SD 1.05):

``` r

ini_df |>
  dplyr::filter(!is.na(neta1), neta1 == neta2) |>
  dplyr::transmute(Eta = name, Variance = signif(est, 4), `SD (log or logit scale)` = signif(sqrt(est), 3)) |>
  knitr::kable()
```

| Eta                      | Variance | SD (log or logit scale) |
|:-------------------------|---------:|------------------------:|
| etalkin_ifngs            |  0.04196 |                   0.205 |
| etalogitfbio_anifrolumab |  1.10200 |                   1.050 |
| etalkel_anifrolumab      |  0.07548 |                   0.275 |
| etalvm_anifrolumab       |  0.09749 |                   0.312 |

## Steady-state hold

With no drug the system must sit exactly at the baseline that Equations
21-25 construct. A deterministic solve over 400 days is a
numerical-identity check.

``` r

ev_ss <- data.frame(id = 1:2, time = 0, evid = 0, amt = 0, cmt = "ifngs", BGENE21_HIGH = c(1, 0)) |>
  dplyr::bind_rows(data.frame(
    id = rep(1:2, each = 8), time = rep(seq(50, 400, by = 50), 2), evid = 0, amt = 0,
    cmt = "ifngs", BGENE21_HIGH = rep(c(1, 0), each = 8)
  ))
ss <- rxode2::rxSolve(mod_typ, ev_ss, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalkin_ifngs', 'etalogitfbio_anifrolumab', 'etalkel_anifrolumab', 'etalvm_anifrolumab'
ss |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    `IFNGS status` = ifelse(BGENE21_HIGH[1] == 1, "high", "low"),
    `Baseline IFNGS` = signif(ifngs_bl[1], 4),
    `Max abs IFNGS change (%)` = signif(max(abs(ifngs_cfb)), 2),
    `Max abs IFN-alpha change (%)` = signif(max(abs(ifna_cfb)), 2),
    `Max abs IC change (%)` = signif(max(abs(ic_cfb)), 2)
  ) |>
  knitr::kable()
```

| id | IFNGS status | Baseline IFNGS | Max abs IFNGS change (%) | Max abs IFN-alpha change (%) | Max abs IC change (%) |
|---:|:---|---:|---:|---:|---:|
| 1 | high | 6.360 | 0 | 0 | 0 |
| 2 | low | 1.928 | 0 | 0 | 0 |

``` r


stopifnot(
  max(abs(ss$ifngs_cfb)) < 1e-6,
  max(abs(ss$ifna_cfb)) < 1e-6,
  max(abs(ss$ic_cfb)) < 1e-6
)
```

The typical IFNGS-high baseline is 6.36 (with baseline IFN-alpha 6.35
pg/mL) and the IFNGS-low baseline is 1.93.

## Maximum IFNGS suppression: closed form

At steady state the two feedbacks cancel for IFN1 production: immune
cells settle at `IC_BL / (rho (1 + E_DAX))`, where
`rho = IFNGS / IFNGS_BL`, and IFN1 production per cell is multiplied by
`rho (1 - I_LIT)`, so total IFN1 production is
`IC_BL (1 - I_LIT) / (1 + E_DAX)`. IFN-IFNAR1 binding is far below
saturation (Kd 100-500 nmol/L against a picomolar IFN pool), so each
IFN1-IFNAR1 complex scales linearly with its free ligand and with free
IFNAR1. The maximal effect of each drug then follows in closed form from
Equation 7 at steady state:

- anifrolumab at saturation blocks all IFNAR1, so IFNGS falls to
  `ifngs_ll`;
- sifalimumab at saturation removes the IFN-alpha complexes only,
  leaving the IFN-beta term;
- daxdilimab at saturation divides IFN1 production by
  `1 + Emax_DAX = 5`;
- litifilimab at saturation multiplies IFN1 production by
  `1 - Imax_LIT = 0.125`.

``` r

p <- setNames(ini_df$est, ini_df$name)
vc <- exp(p[["lvc"]])
ifna_bl <- p[["bl_ifna"]] * vc / p[["mw_ifna"]]
ifnb_bl <- p[["bl_ifna"]] / p[["ifna_fraction"]] * (1 - p[["ifna_fraction"]]) * vc / p[["mw_ifnb"]]
ifnar1_bl <- p[["bl_ifnar1"]] * vc
cx_a <- p[["kon_ifna"]] * ifna_bl * ifnar1_bl / vc / (p[["kd_ifna"]] * p[["kon_ifna"]] + p[["kout_ifna_ifnar1"]])
cx_b <- p[["kon_ifnb"]] * ifnb_bl * ifnar1_bl / vc / (p[["kd_ifnb"]] * p[["kon_ifnb"]] + p[["kout_ifnb_ifnar1"]])
kin <- exp(p[["lkin_ifngs"]])
ifngs_ss <- function(cx) (kin * log(1 + p[["kn"]] * cx * 1e8) + p[["kout_ifngs"]] * p[["ifngs_ll"]]) / p[["kout_ifngs"]]
bl <- ifngs_ss(cx_a + cx_b)
max_closed <- c(
  anifrolumab = 100 * (p[["ifngs_ll"]] - bl) / bl,
  sifalimumab = 100 * (ifngs_ss(cx_b) - bl) / bl,
  daxdilimab = 100 * (ifngs_ss((cx_a + cx_b) / (1 + p[["emax_daxdilimab"]])) - bl) / bl,
  litifilimab = 100 * (ifngs_ss((cx_a + cx_b) * (1 - p[["imax_litifilimab"]])) - bl) / bl
)
max_closed
#> anifrolumab sifalimumab  daxdilimab litifilimab 
#>   -87.98749   -63.13878   -33.99115   -43.25848
```

## Exposure-IFNGS relationship (Figure 5, typical values)

Figure 5 benchmarks the four drugs at typical parameter values. The
model is solved here at a series of constant drug concentrations (a
2000-day zero-order infusion into the central compartment, read at day
600, after every state has reached steady state), for an IFNGS-high
patient.

``` r

drugs <- c(
  anifrolumab = "central_anifrolumab", sifalimumab = "central_sifalimumab",
  daxdilimab = "central_daxdilimab", litifilimab = "central_litifilimab"
)
rates <- 10^seq(-4, 5, by = 0.125)
grid <- expand.grid(rate = rates, drug = names(drugs), stringsAsFactors = FALSE)
grid$id <- seq_len(nrow(grid))
ev_er <- dplyr::bind_rows(
  data.frame(id = grid$id, time = 0, amt = grid$rate * 2000, rate = grid$rate, cmt = drugs[grid$drug], evid = 1),
  data.frame(id = grid$id, time = 600, amt = 0, rate = 0, cmt = drugs[grid$drug], evid = 0)
) |>
  dplyr::mutate(BGENE21_HIGH = 1) |>
  dplyr::arrange(id, time)
er <- rxode2::rxSolve(mod_typ, ev_er, returnType = "data.frame") |>
  dplyr::left_join(grid, by = "id") |>
  dplyr::mutate(conc = dplyr::case_when(
    drug == "anifrolumab" ~ Cc_anifrolumab,
    drug == "sifalimumab" ~ Cc_sifalimumab,
    drug == "daxdilimab" ~ Cc_daxdilimab,
    drug == "litifilimab" ~ Cc_litifilimab
  ))
#> ℹ omega/sigma items treated as zero: 'etalkin_ifngs', 'etalogitfbio_anifrolumab', 'etalkel_anifrolumab', 'etalvm_anifrolumab'

er_summary <- er |>
  dplyr::group_by(drug) |>
  dplyr::summarise(
    max_sim = min(ifngs_cfb),
    ic50_sim = exp(approx(ifngs_cfb, log(conc), xout = min(ifngs_cfb) / 2, ties = mean)$y)
  ) |>
  dplyr::mutate(
    max_closed = max_closed[drug],
    max_pub = c(anifrolumab = -88, sifalimumab = -63.1, daxdilimab = -34, litifilimab = -43.3)[drug],
    ic50_pub = c(anifrolumab = 132, sifalimumab = 7187, daxdilimab = 8.4, litifilimab = 83.6)[drug]
  )

er_summary |>
  dplyr::transmute(
    Drug = drug,
    `Max IFNGS change, closed form (%)` = round(max_closed, 1),
    `Max IFNGS change, simulated (%)` = round(max_sim, 1),
    `Max IFNGS change, published (%)` = max_pub,
    `IC50 simulated (ng/mL)` = signif(ic50_sim, 4),
    `IC50 published (ng/mL)` = ic50_pub,
    `IC50 ratio sim/pub` = round(ic50_sim / ic50_pub, 2)
  ) |>
  knitr::kable(caption = "Maximum IFNGS suppression and half-maximal concentration at typical parameter values. Published values: Results, 'Evaluation of Treatment Benefit' section.")
```

| Drug | Max IFNGS change, closed form (%) | Max IFNGS change, simulated (%) | Max IFNGS change, published (%) | IC50 simulated (ng/mL) | IC50 published (ng/mL) | IC50 ratio sim/pub |
|:---|---:|---:|---:|---:|---:|---:|
| anifrolumab | -88.0 | -88.0 | -88.0 | 132.70 | 132.0 | 1.01 |
| daxdilimab | -34.0 | -34.0 | -34.0 | 32.53 | 8.4 | 3.87 |
| litifilimab | -43.3 | -43.3 | -43.3 | 458.10 | 83.6 | 5.48 |
| sifalimumab | -63.1 | -63.1 | -63.1 | 27830.00 | 7187.0 | 3.87 |

Maximum IFNGS suppression and half-maximal concentration at typical
parameter values. Published values: Results, ‘Evaluation of Treatment
Benefit’ section. {.table}

``` r


stopifnot(
  # Closed form vs published: same parameters, pure arithmetic.
  all(abs(er_summary$max_closed - er_summary$max_pub) < 0.5),
  # Simulated plateau vs closed form (sifalimumab approaches its plateau slowly
  # because the complex half-life is 21 days and binding is not saturated at
  # the top of the grid).
  all(abs(er_summary$max_sim - er_summary$max_closed) < 1),
  # Anifrolumab IC50 matches the published 132 ng/mL.
  abs(er_summary$ic50_sim[er_summary$drug == "anifrolumab"] / 132 - 1) < 0.02
)
```

``` r

er |>
  dplyr::left_join(dplyr::select(er_summary, drug, ic50_sim), by = "drug") |>
  dplyr::filter(conc / ic50_sim >= 1e-4, conc / ic50_sim <= 1e3) |>
  ggplot(aes(conc / ic50_sim, ifngs_cfb, colour = drug)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  labs(x = "Relative concentration (C / IC50)", y = "Change in IFNGS, %", colour = NULL) +
  theme_bw()
```

![Replicates the typical-value panel of Figure 5 of Volkova 2025: IFNGS
change versus concentration relative to each drug's simulated
IC50.](Volkova_2025_SLE_ifn_qsp_files/figure-html/fig5-1.png)

Replicates the typical-value panel of Figure 5 of Volkova 2025: IFNGS
change versus concentration relative to each drug’s simulated IC50.

All four maximal effects reproduce the published values, which pins the
system parameters and each drug’s mechanism. The anifrolumab IC50 also
reproduces (132.7 vs 132 ng/mL). The sifalimumab, daxdilimab and
litifilimab IC50s do not: under a constant-concentration definition they
are 3.8- to 5.7-fold higher than the published values. The paper does
not say how the exposure in Figure 5 was defined for these three drugs;
see Assumptions and deviations.

## Clinical regimens at typical values

Typical-value (IIV zeroed) steady-state responses for the regimens
discussed in the Results and the sensitivity analysis. Doses run for 48
weeks and the trough is read at day 336, one dosing interval after the
last dose.

``` r

mk_reg <- function(id, label, amt, ii, n, cmt, high = 1, extra = NULL) {
  dose_t <- c(seq(0, by = ii, length.out = n), extra)
  dplyr::bind_rows(
    data.frame(id = id, time = dose_t, amt = amt, cmt = cmt, evid = 1),
    data.frame(id = id, time = seq(0, 364, by = 1), amt = 0, cmt = cmt, evid = 0)
  ) |>
    dplyr::mutate(regimen = label, BGENE21_HIGH = high)
}
ev_reg <- dplyr::bind_rows(
  mk_reg(1, "Anifrolumab 300 mg IV Q4W, IFNGS-high", 300, 28, 12, "central_anifrolumab"),
  mk_reg(2, "Anifrolumab 1000 mg IV Q4W, IFNGS-high", 1000, 28, 12, "central_anifrolumab"),
  mk_reg(3, "Anifrolumab 300 mg IV Q4W, IFNGS-low", 300, 28, 12, "central_anifrolumab", high = 0),
  mk_reg(4, "Anifrolumab 100 mg IV Q4W, IFNGS-high", 100, 28, 12, "central_anifrolumab"),
  mk_reg(5, "Anifrolumab 150 mg IV Q4W, IFNGS-high", 150, 28, 12, "central_anifrolumab"),
  mk_reg(6, "Sifalimumab 600 mg IV Q2W, IFNGS-high", 600, 14, 24, "central_sifalimumab"),
  mk_reg(7, "Daxdilimab 50 mg SC Q4W, IFNGS-high", 50, 28, 12, "depot_daxdilimab"),
  mk_reg(8, "Daxdilimab 150 mg SC Q4W, IFNGS-high", 150, 28, 12, "depot_daxdilimab"),
  mk_reg(9, "Litifilimab 450 mg SC Q4W + week 2, IFNGS-high", 450, 28, 12, "depot_litifilimab", extra = 14)
) |>
  dplyr::arrange(id, time, dplyr::desc(evid))
reg <- rxode2::rxSolve(mod_typ, ev_reg, returnType = "data.frame", keep = "regimen")
#> ℹ omega/sigma items treated as zero: 'etalkin_ifngs', 'etalogitfbio_anifrolumab', 'etalkel_anifrolumab', 'etalvm_anifrolumab'

reg_tab <- reg |>
  dplyr::filter(time == 336) |>
  dplyr::distinct(regimen, .keep_all = TRUE) |>
  dplyr::select(regimen, ifngs_cfb, ifna_cfb, ic_cfb)
reg_tab |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(v) round(v, 1))) |>
  dplyr::rename(
    Regimen = regimen,
    `IFNGS change at trough (%)` = ifngs_cfb,
    `Free IFN-alpha change at trough (%)` = ifna_cfb,
    `Immune-cell change at trough (%)` = ic_cfb
  ) |>
  knitr::kable()
```

| Regimen | IFNGS change at trough (%) | Free IFN-alpha change at trough (%) | Immune-cell change at trough (%) |
|:---|---:|---:|---:|
| Anifrolumab 300 mg IV Q4W, IFNGS-high | -87.1 | 1.4 | 686.3 |
| Anifrolumab 1000 mg IV Q4W, IFNGS-high | -87.8 | -1.0 | 712.3 |
| Anifrolumab 300 mg IV Q4W, IFNGS-low | -59.8 | 0.6 | 149.8 |
| Anifrolumab 100 mg IV Q4W, IFNGS-high | -60.9 | 102.2 | 423.2 |
| Anifrolumab 150 mg IV Q4W, IFNGS-high | -82.5 | 22.0 | 597.7 |
| Sifalimumab 600 mg IV Q2W, IFNGS-high | -51.2 | -95.8 | 117.3 |
| Daxdilimab 50 mg SC Q4W, IFNGS-high | -33.0 | -78.9 | -68.6 |
| Daxdilimab 150 mg SC Q4W, IFNGS-high | -33.6 | -79.6 | -69.3 |
| Litifilimab 450 mg SC Q4W + week 2, IFNGS-high | -43.3 | -87.5 | 76.2 |

Comparison with the values the paper reports:

- Anifrolumab 300 mg Q4W, IFNGS-high: published individual-prediction
  trough -86.1% (95% CI -87 to -85.1); typical value -87.1%. The 1000 mg
  regimen adds almost nothing, as the paper notes (“no discernible
  difference between 300 mg and 1000 mg”).
- Anifrolumab, IFNGS-low: published -30.1% is a mean of individual post
  hoc predictions; individual baseline IFN-alpha (lower in IFNGS-low
  patients) and per-subject IFNGS floors are not available, so the
  typical value is not comparable.
- Daxdilimab: free IFN-alpha falls by about 80%, matching the published
  “IFNa decrease close to 80% for both 50 mg and 150 mg”. The published
  trough IFNGS response (-27%, virtual-population mean) is a little
  smaller than the typical value.
- Litifilimab: published IFNGS -44.1% and IFN-alpha -88.5%
  (virtual-population means); typical values -43.3% and -87.5%.

``` r

stopifnot(
  # Deterministic solves at the same parameters as the closed form.
  abs(reg_tab$ifngs_cfb[1] - (-86.1)) < 2,
  abs(reg_tab$ifna_cfb[7] - (-80)) < 3,
  abs(reg_tab$ifna_cfb[8] - (-80)) < 3,
  abs(reg_tab$ifngs_cfb[9] - (-44.1)) < 2,
  abs(reg_tab$ifna_cfb[9] - (-88.5)) < 2
)
```

``` r

reg |>
  dplyr::filter(id %in% c(1, 3, 6, 7, 9)) |>
  dplyr::select(time, regimen, `IFNGS` = ifngs_cfb, `Free IFN-alpha` = ifna_cfb) |>
  tidyr::pivot_longer(c(`IFNGS`, `Free IFN-alpha`)) |>
  ggplot(aes(time / 7, value, colour = regimen)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  labs(x = "Time (weeks)", y = "Change from baseline, %", colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical")
```

![Typical-value analogue of Figure 2C (anifrolumab) and Figure 3 (other
drugs): IFNGS and free IFN-alpha change from baseline. Anifrolumab
produces a step fall in IFN-alpha at treatment start followed by
recovery over about 16 weeks, as described in the
Results.](Volkova_2025_SLE_ifn_qsp_files/figure-html/fig2c-1.png)

Typical-value analogue of Figure 2C (anifrolumab) and Figure 3 (other
drugs): IFNGS and free IFN-alpha change from baseline. Anifrolumab
produces a step fall in IFN-alpha at treatment start followed by
recovery over about 16 weeks, as described in the Results.

``` r

ifna_ani <- reg |> dplyr::filter(id == 1)
stopifnot(
  # Step decline at treatment start ...
  ifna_ani$ifna_cfb[ifna_ani$time == 7] < -50,
  # ... and recovery to within 25% of baseline by week 16.
  abs(ifna_ani$ifna_cfb[ifna_ani$time == 112]) < 25
)
```

## Anifrolumab population simulation

Stochastic simulation with the four published random effects (IFNGS
synthesis rate, anifrolumab IV bioavailability, linear elimination and
Vmax) for IFNGS-high patients at the four IV Q4W doses evaluated in the
paper, 200 virtual patients per arm. As in the paper’s virtual
populations, every patient is IFNGS-high. Baseline IFN-alpha is held at
its typical value, because the variance of its log-normal distribution
is not reported.

``` r

n_per_arm <- 200
doses <- c(100, 150, 300, 1000)
# Dense sampling in the first and the steady-state dosing interval (with a
# point just after each of those doses so Cmax is captured), weekly between.
obs_t <- sort(unique(c(seq(0, 28, by = 0.5), 0.01, seq(35, 329, by = 7), seq(336, 364, by = 0.5), 336.01)))
ev_pop <- dplyr::bind_rows(lapply(seq_along(doses), function(j) {
  ids <- (j - 1) * n_per_arm + seq_len(n_per_arm)
  dplyr::bind_rows(
    data.frame(id = rep(ids, each = 14), time = rep(seq(0, 364, by = 28), n_per_arm), amt = doses[j], evid = 1),
    data.frame(id = rep(ids, each = length(obs_t)), time = rep(obs_t, n_per_arm), amt = 0, evid = 0)
  ) |>
    dplyr::mutate(treatment = paste0(doses[j], " mg IV Q4W"))
})) |>
  dplyr::filter(!(evid == 1 & time == 364)) |>
  dplyr::mutate(cmt = "central_anifrolumab", BGENE21_HIGH = 1) |>
  # Observation rows sort before a dose at the same time, so a sample at a
  # dosing time is the pre-dose trough.
  dplyr::arrange(id, time, evid)

rxode2::rxSetSeed(20250701)
pop <- rxode2::rxSolve(ui, ev_pop, returnType = "data.frame", keep = "treatment")
```

``` r

pop |>
  dplyr::filter(time >= 336) |>
  dplyr::group_by(treatment, time) |>
  dplyr::summarise(
    med = median(Cc_anifrolumab), lo = quantile(Cc_anifrolumab, 0.05),
    hi = quantile(Cc_anifrolumab, 0.95), .groups = "drop"
  ) |>
  ggplot(aes(time - 336, med, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line() +
  scale_y_log10() +
  labs(x = "Time after dose at steady state (days)", y = "Anifrolumab (ng/mL)", colour = NULL, fill = NULL) +
  theme_bw()
```

![Analogue of Figure 2B: simulated anifrolumab concentrations (median
and 90% interval) over the last dosing interval at steady
state.](Volkova_2025_SLE_ifn_qsp_files/figure-html/fig2b-1.png)

Analogue of Figure 2B: simulated anifrolumab concentrations (median and
90% interval) over the last dosing interval at steady state.

Steady-state trough IFNGS change (day 336, before the next dose) by
dose:

``` r

pop_trough <- pop |>
  dplyr::filter(time == 336) |>
  dplyr::distinct(id, .keep_all = TRUE) |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(
    mean_cfb = mean(ifngs_cfb), median_cfb = median(ifngs_cfb),
    q10 = quantile(ifngs_cfb, 0.1), q90 = quantile(ifngs_cfb, 0.9)
  )
pop_trough |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(v) round(v, 1))) |>
  dplyr::rename(
    Regimen = treatment, `Mean (%)` = mean_cfb, `Median (%)` = median_cfb,
    `10th percentile (%)` = q10, `90th percentile (%)` = q90
  ) |>
  knitr::kable()
```

| Regimen        | Mean (%) | Median (%) | 10th percentile (%) | 90th percentile (%) |
|:---------------|---------:|-----------:|--------------------:|--------------------:|
| 100 mg IV Q4W  |    -48.1 |      -55.7 |               -82.3 |                -2.7 |
| 1000 mg IV Q4W |    -87.8 |      -88.0 |               -90.5 |               -84.7 |
| 150 mg IV Q4W  |    -67.2 |      -78.5 |               -86.7 |               -25.3 |
| 300 mg IV Q4W  |    -85.0 |      -86.8 |               -89.4 |               -81.7 |

``` r


stopifnot(
  # Exposure-driven ordering of the centre of each arm: lower doses lose
  # IFNGS suppression at trough because saturable clearance takes over.
  median(pop$ifngs_cfb[pop$time == 336 & pop$treatment == "300 mg IV Q4W"]) < -80,
  median(pop$ifngs_cfb[pop$time == 336 & pop$treatment == "100 mg IV Q4W"]) >
    median(pop$ifngs_cfb[pop$time == 336 & pop$treatment == "300 mg IV Q4W"]) + 10
)
```

The paper reports, for its 100 mg Q4W virtual population, steady-state
“trough and peak” IFNGS changes of -66.5% and -41.3%; the order of the
two labels appears transposed, because suppression is necessarily
weakest at the trough. The simulated 100 mg trough mean (-48.1%) is
close to -41.3%. For 150 mg Q4W the paper reports predose values no
deeper than -54.9%; the simulated mean is deeper. The published virtual
populations also sampled baseline IFN-alpha from a log-normal
distribution whose spread is not reported, which these simulations
cannot reproduce.

## PKNCA of simulated anifrolumab

The paper reports no NCA for any drug, so there is no published table to
compare against. NCA over the first and the steady-state dosing interval
documents the simulated anifrolumab exposure and its nonlinearity (the
saturable elimination term has Km = 38.9 nmol, about 5.6 mg in the
central compartment).

``` r

conc_nca <- pop |>
  dplyr::filter(!is.na(Cc_anifrolumab)) |>
  dplyr::select(id, time, Cc = Cc_anifrolumab, treatment)
dose_nca <- ev_pop |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, treatment)

conc_obj <- PKNCA::PKNCAconc(conc_nca, Cc ~ time | treatment + id, concu = "ng/mL", timeu = "day")
dose_obj <- PKNCA::PKNCAdose(dose_nca, amt ~ time | treatment + id, doseu = "mg")
intervals <- data.frame(
  start = c(0, 336), end = c(28, 364),
  cmax = TRUE, cmin = TRUE, auclast = TRUE, cav = TRUE
)
nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
nca_tab <- as.data.frame(nca$result) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "cmin", "cav")) |>
  dplyr::group_by(treatment, start, PPTESTCD) |>
  dplyr::summarise(median = signif(median(PPORRES), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = median) |>
  dplyr::mutate(interval = ifelse(start == 0, "Dose 1 (0-28 d)", "Steady state (336-364 d)"))
nca_tab |>
  dplyr::select(treatment, interval, cmax, cmin, cav) |>
  dplyr::rename(
    Regimen = treatment, Interval = interval,
    `Cmax (ng/mL)` = cmax, `Cmin (ng/mL)` = cmin, `Cavg (ng/mL)` = cav
  ) |>
  knitr::kable(caption = "Median simulated NCA parameters for anifrolumab (no published NCA to compare against).")
```

| Regimen | Interval | Cmax (ng/mL) | Cmin (ng/mL) | Cavg (ng/mL) |
|:---|:---|---:|---:|---:|
| 100 mg IV Q4W | Dose 1 (0-28 d) | 29100 | 473 | 5340 |
| 100 mg IV Q4W | Steady state (336-364 d) | 29200 | 508 | 5310 |
| 1000 mg IV Q4W | Dose 1 (0-28 d) | 338000 | 55400 | 92400 |
| 1000 mg IV Q4W | Steady state (336-364 d) | 401000 | 120000 | 175000 |
| 150 mg IV Q4W | Dose 1 (0-28 d) | 45200 | 1790 | 9410 |
| 150 mg IV Q4W | Steady state (336-364 d) | 46400 | 2310 | 10500 |
| 300 mg IV Q4W | Dose 1 (0-28 d) | 98600 | 11000 | 24700 |
| 300 mg IV Q4W | Steady state (336-364 d) | 110000 | 22600 | 38600 |

Median simulated NCA parameters for anifrolumab (no published NCA to
compare against). {.table}

``` r


cav_ss <- nca_tab |> dplyr::filter(start == 336)
stopifnot(
  # Saturable elimination makes exposure more than dose-proportional between
  # 100 and 300 mg.
  cav_ss$cav[cav_ss$treatment == "300 mg IV Q4W"] / cav_ss$cav[cav_ss$treatment == "100 mg IV Q4W"] > 3
)
```

## Assumptions and deviations

- **Two Table S2 values differ from the deposited Monolix code (Data
  S2); the code values are used.** Anifrolumab molecular weight is
  145,119.2 g/mol in the code but 148 kDa in Table S2 and the Data S3
  text; the anifrolumab-IFNAR1 internalisation rate is 77.1 /day in the
  code and 77.4 /day in Table S2. The molecular weight matters: with
  145.1 kDa the model reproduces the published anifrolumab IC50 of 132
  ng/mL (132.7 ng/mL above), whereas 148 kDa gives 135.3 ng/mL. The
  internalisation-rate difference changes the IC50 by less than 0.1%.
- **Sifalimumab, daxdilimab and litifilimab IC50s (Figure 5) are not
  reproduced.** With constant concentrations at steady state they are
  3.8- to 5.7-fold above the published 7187, 8.4 and 83.6 ng/mL, while
  their maximal effects reproduce exactly. The drug-specific potency
  parameters (sifalimumab Kd / koff / complex elimination, daxdilimab
  EC50, litifilimab IC50 / Hill) agree between Table S2 and the code.
  The paper does not define the exposure used on the Figure 5 axis for
  these drugs, so the discrepancy is recorded rather than resolved; no
  parameter has been adjusted.
- **Baseline IFN-alpha is a fixed parameter, `bl_ifna = 6.35` pg/mL**
  (the Table S3 typical value). The paper treats it as an individual
  per-subject value: Table S2 calls it a covariate, and the Monolix code
  gives it a log-normal individual distribution whose variance is not
  reported. Set `bl_ifna` to a patient’s measured baseline IFN-alpha to
  individualise.
- **The residual IFNGS floor is fixed at 0.764.** Table S2 defines it
  per subject as min(lowest observed IFNGS, 0.764); 0.764 is the value
  for any patient whose observed IFNGS never falls below it.
- **No residual error.** The paper counts four residual-error parameters
  but does not report their values, so the model carries none.
- **IIV on baseline IFNAR1 is not encoded.** Data S3 Supplementary File
  3 mentions random effects on baseline IFNAR1 and kin_IFNGS in the
  first calibration step. Table S2 (IFNAR1_BL = 0.26 nmol/L, no CV,
  literature source), the abstract (“one random effect” among the
  systemic parameters) and the code (no IFNAR1 random effect) agree that
  only kin_IFNGS carries one, and that is what is implemented.
- **Transcription notes on Table S2.** The description column swaps the
  directions of anifrolumab k12 and k21 (Equation 9 and the code use k12
  for central to peripheral); the descriptions of KD_SIF and koff_SIF
  are copied from the Vm / Km rows; the litifilimab Imax / IC50 rows say
  “effect on ICs turnover”, whereas Equations 2-3 and the code place the
  litifilimab effect on IFN1 production, which is what is implemented.
  The main text twice writes “axatilimab” where the anifrolumab PK
  sub-model is meant.
- **Shared plasma volume.** All four drugs and the IFN1 species use the
  same Vc = 2.75 L, as in the code, including litifilimab and
  sifalimumab whose source PK models had their own central volumes.
- **IFNGS-low covariate.** The source indicator is GRNUM = 1 for
  IFNGS-low; the model uses the canonical `BGENE21_HIGH` and applies the
  published coefficient to `1 - BGENE21_HIGH`.
- **Doses are in mg for all drugs.** Sifalimumab doses reported in mg/kg
  must be converted by the user; the regimen table above uses a flat 600
  mg.
- No erratum was found for this article (Europe PMC search, 2026-11-05).
