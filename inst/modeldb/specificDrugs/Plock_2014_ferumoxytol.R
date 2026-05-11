Plock_2014_ferumoxytol <- function() {
  description <- "Two-compartment population PK model with Michaelis-Menten elimination for IV ferumoxytol in healthy adults and adults with chronic kidney disease (Plock 2014). Encodes the typical non-dialysing-patient form; the haemodialysis-driven time-varying central volume (VSLOPE) and the within-session weight-loss effect on V1 (WLO) are described in the vignette but not enabled in this model file."
  reference   <- paste(
    "Plock N, Facius A, Lahu G, Wood N, Frigo T, Deveney A, Aceves P.",
    "Population Pharmacokinetic Meta-Analysis to Bridge Ferumoxytol Plasma",
    "Pharmacokinetics Across Populations. Clin Pharmacokinet. 2015;54(4):385-395.",
    "doi:10.1007/s40262-014-0203-9.",
    sep = " "
  )
  vignette <- "Plock_2014_ferumoxytol"
  units    <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Plock 2014 Eq 1 linear centered relation on V1:",
        "V1 = V1_pop * (1 + (WT - 80) * 0.614 / 100).",
        "Reference value 80 kg is the cohort-pooled median (Table 3:",
        "median 79 kg in HV studies A and B, 83.9 kg in CKD study C;",
        "pooled n-weighted central tendency ~80 kg)."
      ),
      source_name        = "WGT"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = paste(
        "Plock 2014 Eq 2 fractional-change relation on V1 with female as the",
        "non-reference category: V1 = V1_pop * (1 + SEXF * -18.3 / 100).",
        "Source paper encoded SEX as 1 = male / 2 = female; here mapped to the",
        "canonical SEXF (1 = female, 0 = male)."
      ),
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 111,
    n_studies      = 3,
    age_range      = "18-77 years (pooled across the three studies; study C [CKD] median 64 years)",
    age_median     = "studies A and B (HV) medians 31 and 30 years; study C (CKD) median 64 years",
    weight_range   = "46.5-115 kg (pooled across studies A, B, C; Table 3)",
    weight_median  = "studies A and B (HV) medians both 79 kg; study C (CKD) median 83.9 kg",
    sex_female_pct = 45,
    race_ethnicity = c(White = 16, Black = 73, Hispanic = 9, Asian_Other = 2),
    disease_state  = "Healthy volunteers (studies A and B) and adults with chronic kidney disease stage 5D on haemodialysis (study C). Used to bridge ferumoxytol PK from healthy and CKD-on-HD populations to the broader iron-deficiency-anaemia population.",
    dose_range     = "Study A: ascending 1, 2, 4 mg/kg single IV (rates 30 mg/20 s, 30 mg/10 s, 30 mg/s; one cohort 60 mg/min). Study B: 2 x 510 mg IV (17 mL over 17 s) administered 24 h apart. Study C: 125 or 250 mg single IV over 5 min within 30 min after dialysis start.",
    regions        = "United States.",
    notes          = paste(
      "Pooled meta-analysis of three Phase I studies (Table 1). 91 HV PK",
      "subjects from studies A (n=33) and B (n=58) and 20 CKD-stage-5D-on-HD",
      "subjects from study C. 1,686 observations after pre-dose exclusion,",
      "29 below LLOQ. LLOQ 5.83, 6.0, and 11.16 ug/mL for studies A, B,",
      "and C respectively. Complete covariate information for all 111",
      "subjects. The haemodialysis-specific effects (3-h linear V1 decline",
      "during dialysis with VSLOPE = -0.198 L/h, 3-h symmetric recovery,",
      "and the WLO covariate on initial V1 with effect 7.22 % per kg of",
      "intra-dialysis weight loss) are not enabled in this model file;",
      "see the vignette for parameter values and an example simulation."
    )
  )

  ini({
    # Structural PK parameters (Plock 2014 Table 4 final model)
    lvmax <- log(16.5);    label("Michaelis-Menten Vmax (mg/h)")                                    # Table 4: Vmax = 16.5 mg/h (RSE 4.8%)
    lkm   <- log(96.7);    label("Michaelis-Menten Km (mg/L)")                                       # Table 4: Km = 96.7 mg/L (RSE 7.2%)
    lvc   <- log(2.78);    label("Central volume V1 at WT = 80 kg, male (L)")                        # Table 4: V1 = 2.78 L (RSE 3.3%)
    lq    <- log(0.0289);  label("Intercompartmental clearance Q (L/h)")                             # Table 4: Q  = 0.0289 L/h (RSE 8.6%)
    lvp   <- log(0.348);   label("Peripheral volume V2 (L)")                                         # Table 4: V2 = 0.348 L (RSE 11.5%)

    # Covariate effects on V1 (linear / fractional-change forms; Plock 2014 Methods Eqs 1-2)
    e_wt_vc   <- 0.614;    label("WT effect on V1 (% per kg, linear centered at 80 kg)")             # Table 4: WGT on V1 = 0.614 %/kg (RSE 24.9%)
    e_sexf_vc <- -18.3;    label("SEXF effect on V1 (% reduction for females)")                      # Table 4: V1 change for females = -18.3% (RSE 16.7%)

    # IIV (Plock 2014 Table 4 final model; CV% on the variance scale per table footnote)
    # Internal log-normal variance: omega^2 = log(1 + CV^2)
    etalvmax ~ 0.03016    # Table 4: BSV Vmax = 17.5 %CV  -> log(1 + 0.175^2)
    etalvc   ~ 0.02783    # Table 4: BSV V1   = 16.8 %CV  -> log(1 + 0.168^2)
    etalvp   ~ 0.24509    # Table 4: BSV V2   = 52.7 %CV  -> log(1 + 0.527^2)

    # Residual error (combined; Plock 2014 Table 4)
    propSd <- 0.0785;      label("Proportional residual error (fraction)")                           # Table 4: 7.85 %CV
    addSd  <- 1.40;        label("Additive residual error (ug/mL)")                                  # Table 4: 1.40 ug/mL
  })

  model({
    # Reference body weight for the WT centering (matches e_wt_vc above; see covariateData$WT$notes)
    wt_ref <- 80

    # Individual PK parameters
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)
    vc   <- exp(lvc + etalvc) *
              (1 + e_wt_vc   / 100 * (WT - wt_ref)) *
              (1 + e_sexf_vc / 100 * SEXF)
    q    <- exp(lq)
    vp   <- exp(lvp + etalvp)

    # Two-compartment IV with saturable Michaelis-Menten elimination
    # (Plock 2014 Eqs 3-4; the time-varying-V1 mechanism is omitted -- see vignette).
    Cc <- central / vc
    d/dt(central)     <- -q / vc * central + q / vp * peripheral1 -
                           vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
