Rascher_2025_empagliflozin <- function() {
  description <- paste(
    "Population pharmacokinetic model for oral empagliflozin in paediatric",
    "patients aged 10-17 years with type 2 diabetes mellitus (DINAMO trial).",
    "Two-compartment model with sequential zero-order (duration D1) then",
    "first-order (ka) absorption and linear elimination, re-estimated in a",
    "Bayesian (NUTS) framework from the adult Baron 2016 structure.",
    "Fixed allometric weight scaling (exponent 0.75 on CL/F and Q/F, 1 on",
    "V2/F and V3/F, reference 70 kg); power eGFR effect on CL/F (eGFR",
    "truncated at 120 mL/min/1.73 m^2, reference 90); multiplicative Black,",
    "Asian and female effects on CL/F. IIV on CL/F only."
  )
  reference <- paste(
    "Rascher J, Cheng S, Johnston C, Haertter S, Marquard J, Tartakovsky I,",
    "Laffel L, Wojtyniak J-G.",
    "Pharmacokinetics and pharmacodynamics of empagliflozin in paediatric",
    "patients aged 10-17 years with type 2 diabetes mellitus.",
    "Br J Clin Pharmacol. 2025;91(8):2390-2400. doi:10.1002/bcp.70096."
  )
  vignette <- "Rascher_2025_empagliflozin"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed allometric scaling with reference 70 kg: exponent 0.75 on CL/F and Q/F, 1 on V2/F and V3/F (Methods 'PopPK for adult and paediatric patients'; Supplement CL/F equation WT_CL = (WT/70)^theta).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalised)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "eGFR computed with the Zappitelli equation in DINAMO (Table S2 footnote;",
        "the adult predecessor analysis used MDRD). Values above 120 mL/min/1.73 m^2",
        "are truncated to 120 inside the model before entering the power effect",
        "(EGFR/90)^theta on CL/F (Supplement CL/F equation and definition of EGFR_i)."
      ),
      source_name = "EGFR"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative effect cl_female^SEXF on CL/F (Table 3 FEMALE CL/F = exp(theta) = 1.19). The Supplement equation applies the effect when FEMALE = 1, so the reference is male; see the vignette for the conflicting Results sentence.",
      source_name = "FEMALE"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Multiplicative effect cl_black^RACE_BLACK on CL/F (Table 3 BLACK CL/F = 0.897).",
      source_name = "BLACK"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Multiplicative effect cl_asian^RACE_ASIAN on CL/F (Table 3 ASIAN CL/F = 0.934).",
      source_name = "ASIAN"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "The adult predecessor model carried AGE_CL = theta * ln(AGE/55) on CL/F (Supplement CL/F equation); the paediatric final model removed it (Discussion 'PopPK model'; Table S3 shows no DINAMO estimate)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 74L,
    n_studies = 1L,
    n_observations = 223L,
    studies = "DINAMO (NCT03429543), phase 3 randomised, double-blind trial in young people with type 2 diabetes.",
    age_range = "10-17 years (median 14.5; Table S2)",
    weight_range = "42.5-169 kg (median 90.0; Table S2)",
    egfr_range = "89.7-241 mL/min/1.73 m^2 (median 125, Zappitelli equation; Table S2)",
    sex_female_pct = 62.2,
    race_ethnicity = c(
      White = 50.0,
      Black = 32.4,
      Asian = 5.4,
      `American Indian or Alaska Native` = 4.1,
      `Native Hawaiian or Other Pacific Islander` = 1.4,
      Multiple = 5.4,
      Unknown = 1.4
    ),
    disease_state = "Type 2 diabetes mellitus (paediatric, 10-17 years)",
    dose_range = "Empagliflozin 10 mg (n = 46), 25 mg (n = 16) or 10 mg re-randomised to 25 mg (n = 12) orally once daily (Table 1).",
    regions = "Multinational (DINAMO).",
    notes = paste(
      "278 PK samples, 55 (19.8%) below the limit of quantification and excluded",
      "(Table 1; Supplement). Reference subject for the typical values: male",
      "(per the Supplement equation), non-Black, non-Asian, eGFR 90",
      "mL/min/1.73 m^2, 70 kg."
    )
  )

  ini({
    # Table 3 (DINAMO model posterior medians; parameters estimated on the
    # log scale and printed back-transformed).
    lcl <- log(6.74); label("Apparent clearance CL/F for the reference subject (L/h)") # Table 3 CL/F = 6.74 L/h (95% CDI 5.71, 7.90)
    lvc <- log(4.12); label("Apparent central volume V2/F for the reference subject (L)") # Table 3 V2/F = 4.12 L (95% CDI 1.30, 6.63)
    lka <- log(0.239); label("First-order absorption rate constant ka (1/h)") # Table 3 KA = 0.239 1/h (95% CDI 0.232, 0.246)
    lq <- log(5.51); label("Apparent intercompartmental clearance Q/F for the reference subject (L/h)") # Table 3 Q/F = 5.51 L/h (95% CDI 5.22, 5.83)
    lvp <- log(71.7); label("Apparent peripheral volume V3/F for the reference subject (L)") # Table 3 V3/F = 71.7 L (95% CDI 67.6, 76.2)
    ld1 <- log(0.326); label("Duration of the zero-order input into the depot D1 (h)") # Table 3 D1 = 0.326 h (95% CDI 0.199, 0.542)

    # Fixed allometric exponents (Methods: clearances 0.75, volumes 1)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL/F (unitless)") # Methods 'PopPK' paragraph 1: clearances exponent 0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent of (WT/70) on Q/F (unitless)") # Methods 'PopPK' paragraph 1: clearances exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V2/F (unitless)") # Methods 'PopPK' paragraph 1: volumes exponent 1
    e_wt_vp <- fixed(1); label("Allometric exponent of (WT/70) on V3/F (unitless)") # Methods 'PopPK' paragraph 1: volumes exponent 1

    # Covariate effects on CL/F
    e_crcl_cl <- 0.407; label("Power exponent of (min(eGFR, 120)/90) on CL/F (unitless)") # Table 3 EGFR CL/F = 0.407 (95% CDI 0.363, 0.454)
    cl_black <- 0.897; label("Multiplicative effect of Black race on CL/F (ratio)") # Table 3 BLACK CL/F = exp(theta 8) = 0.897
    cl_asian <- 0.934; label("Multiplicative effect of Asian race on CL/F (ratio)") # Table 3 ASIAN CL/F = exp(theta 9) = 0.934
    cl_female <- 1.19; label("Multiplicative effect of female sex on CL/F (ratio)") # Table 3 FEMALE CL/F = exp(theta 10) = 1.19

    # IIV: Table 3 IIV-CL/F 33.3 %CV with CV = sqrt(exp(omega^2) - 1)
    etalcl ~ 0.1052 # Table 3 omega CL/F 33.3 %CV -> omega^2 = log(1 + 0.333^2) = 0.1052

    # Residual error (Table 3; sigma reported as CV = sqrt(sigma^2) or SD = sqrt(sigma^2))
    propSd <- 0.474; label("Proportional residual error, non-outlier observations (fraction)") # Table 3 Sigma 1,1 = 47.4 %CV
    addSd <- 1.67; label("Additive residual error, non-outlier observations (nmol/L)") # Table 3 Sigma 2,2 SD = 1.67 nmol/L
  })
  model({
    # Empagliflozin C23H27ClO7, 450.91 g/mol: converts mg/L to nmol/L, the
    # concentration unit of Table 3 and of AUC50 in the companion HbA1c model.
    mw_empa <- 450.91
    nmol_per_mg <- 1e6 / mw_empa

    # Supplement: EGFR_i is truncated at 120 mL/min/1.73 m^2, reference 90.
    egfr_cl <- CRCL
    if (CRCL > 120) egfr_cl <- 120

    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (egfr_cl / 90)^e_crcl_cl *
      cl_black^RACE_BLACK *
      cl_asian^RACE_ASIAN *
      cl_female^SEXF
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    ka <- exp(lka)
    d1 <- exp(ld1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Sequential zero-order then first-order absorption: each dose is infused
    # into the depot over D1 (dose records need rate = -2) and the depot
    # empties into central at first-order ka.
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    dur(depot) <- d1

    Cc <- central / vc * nmol_per_mg
    Cc ~ prop(propSd) + add(addSd)
  })
}
