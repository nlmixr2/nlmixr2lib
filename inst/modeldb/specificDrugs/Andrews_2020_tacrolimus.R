Andrews_2020_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and a fixed absorption lag time for twice-daily oral immediate-release tacrolimus (Prograft capsules or Modigraf granules) in paediatric kidney transplant recipients during the first 6 weeks post-transplantation (Andrews 2020 FINAL model; improved re-fit of Andrews 2017 on an expanded cohort of 95 children). Apparent oral clearance CL/F scales with body weight at an estimated allometric exponent of 0.62 referenced to 70 kg; apparent central and peripheral volumes V1/F and V2/F scale with a fixed exponent of 1; Q/F and ka are not scaled. CL/F additionally varies with CYP3A5 expresser status (1.4-fold for *1/*1 or *1/*3), haematocrit (power -0.60 centred at 0.29 L/L) and serum creatinine (power -0.1 centred at 84 umol/L). Diagonal IIV on ka, CL/F, V1/F and V2/F. Combined additive + proportional residual error with separate magnitudes for immunoassay and LC-MS/MS samples. Inter-occasion variability on CL/F (20.1% CV) is reported but not encoded (no occasion definition in the source). The companion reduced-covariate starting-dose model is Andrews_2020_tacrolimus_startingdose."
  reference <- "Andrews LM, de Winter BCM, Cornelissen EAM, de Jong H, Hesselink DA, Schreuder MF, Bruggemann RJM, van Gelder T, Cransberg K. A Population Pharmacokinetic Model Does Not Predict the Optimal Starting Dose of Tacrolimus in Pediatric Renal Transplant Recipients in a Prospective Study: Lessons Learned and Model Improvement. Clin Pharmacokinet. 2020;59(5):591-603. doi:10.1007/s40262-019-00831-8 (published online 26 October 2019)."
  vignette <- "Andrews_2020_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying within subject over the first 6 weeks post-transplant. Allometric power scaling referenced to 70 kg: estimated exponent 0.62 on CL/F (Table 4 final model, Equation 2), fixed exponent 1 on V1/F and V2/F (Results 3.2.1). Table 4 labels Q/F in 'L/h' (not 'L/h/70 kg'), so Q/F is not scaled; ka is not scaled. Model-building cohort median 32.0 kg, range 10.4-87.5 kg (Table 3).",
      source_name = "weight"
    ),
    HCT = list(
      description = "Haematocrit, packed red blood cell volume fraction",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Andrews 2020 reports haematocrit as a volume fraction in L/L (Table 3 median 0.29, range 0.16-0.52) and centres Equation (2) at 0.29 L/L. The canonical HCT column is on the percent scale, so the model centres at 29 (%) -- numerically identical to (HCT_fraction / 0.29). Multiply an L/L column by 100 before use. Enters CL/F as (HCT/29)^-0.60.",
      source_name = "Hematocrit"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying ('serum creatinine during follow-up', Results 3.2.1). Enters CL/F as (CREAT/84)^-0.1 per Equation (2); 84 umol/L is the model-building cohort median (Table 3, range 12-1454).",
      source_name = "Creatinine"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if a non-expresser (*3/*3 or *3/*7).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 non-expresser, *3/*3)",
      notes = "Time-fixed germline genotype. Equation (2) writes the CL/F factor as 1.0 for *3/*3 and 1.4 for *1/*3 or *1/*1. Table 3 lists *3/*7 (n = 2) and Unknown (n = 27, 28%) strata; the paper does not state how Unknown-genotype children were coded. Following the Andrews 2017 predecessor (which pooled Unknown with *3/*3), a user without genotype should set CYP3A5_EXPR = 0.",
      source_name = "CYP3A5"
    ),
    IMMUNOASSAY = list(
      description = "Per-sample bioanalytical assay indicator: 1 if the tacrolimus whole-blood concentration was measured by immunoassay; 0 if measured by LC-MS/MS.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (LC-MS/MS reference method)",
      notes = "Per-sample (per-row) indicator. Table 3: 64 immunoassay samples (4.8%) and 1274 LC-MS/MS samples (95.2%). Switches the additive and proportional residual-error magnitudes (Table 4). Set to 0 for simulation of LC-MS/MS-measured concentrations.",
      source_name = "assay"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 95L,
    n_studies = 2L,
    n_observations = 1338L,
    age_range = "1.6-17.9 years",
    age_median = "11.4 years",
    weight_range = "10.4-87.5 kg",
    weight_median = "32.0 kg",
    sex_female_pct = 39,
    race_ethnicity = c(Caucasian = 74, African = 9, Asian = 2, Other = 14),
    disease_state = "Paediatric kidney transplant recipients during the first 6 weeks post-transplantation, treated with basiliximab, twice-daily immediate-release tacrolimus, mycophenolic acid and a 5-day course of glucocorticoids (TWIST protocol). 78% living-donor, 22% deceased-donor grafts.",
    dose_range = "Twice-daily (q12h) oral tacrolimus; starting dose 0.3 mg/kg/day (historic patients) or per the Andrews 2017 algorithm (0.27-1.33 mg/kg/day, prospective-trial patients), then adjusted by therapeutic drug monitoring to a C0 target of 10-15 ng/mL.",
    regions = "Netherlands (Erasmus MC-Sophia Children's Hospital, Rotterdam; Radboudumc Amalia Children's Hospital, Nijmegen).",
    cyp3a5_distribution = "*1/*1 3 (3%), *1/*3 11 (12%), *3/*3 52 (55%), *3/*7 2 (2%), Unknown 27 (28%) (Table 3).",
    laboratory = "Haematocrit 0.29 L/L (0.16-0.52); creatinine 84 umol/L (12-1454); eGFR 63 mL/min (2.9-274) (Table 3, median and range).",
    assay = "Whole-blood tacrolimus by LC-MS/MS (95.2% of samples) or immunoassay (4.8%) (Table 3).",
    sampling_window = "1338 samples, 0-42 days post-transplant (Table 3); at least one PK profile in 90 patients at a median 12.2 days post-transplant (Results 3.2).",
    iov_structure = "Inter-occasion variability on CL/F of 20.1% CV (Table 4 final model) is not encoded structurally: the paper defines no occasion (it refers to the Andrews 2017 methods), matching the handling of the Andrews 2017 and Andrews 2019 tacrolimus models in this library.",
    notes = "Model-building cohort of 95 children: 45 from the Andrews 2017 original cohort (with additional retrospectively retrieved PK data), 16 from the prospective trial reported in the same paper, and 34 further children transplanted March 2012 - October 2017 (Methods 2.2)."
  )

  ini({
    # Final-model estimates, Andrews 2020 Table 4 column 'Final model (RSE %)
    # [shrinkage]'. All clearances and volumes are apparent (CL/F, V/F). The
    # typical values refer to a 70 kg CYP3A5 non-expresser with haematocrit
    # 0.29 L/L and serum creatinine 84 umol/L.
    ltlag <- fixed(log(0.41)); label("Absorption lag time tlag (h)") # Table 4 't lag (h) FIX' = 0.41
    lka <- log(1.7); label("Absorption rate constant ka (1/h)") # Table 4 final 'k a (L/h)' = 1.7 (RSE 11%); the 'L/h' unit label is a typo for 1/h
    lcl <- log(36.6); label("Apparent oral clearance CL/F at 70 kg, CYP3A5 non-expresser, HCT 0.29 L/L, creatinine 84 umol/L (L/h)") # Table 4 final 'CL/ F (L/h/70 kg)' = 36.6 (RSE 12%); also the leading coefficient of Equation (2)
    lvc <- log(496); label("Apparent central volume V1/F at 70 kg (L)") # Table 4 final 'V 1 / F (L/70 kg)' = 496 (RSE 22%)
    lq <- log(31.7); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 4 final 'Q / F (L/h)' = 31.7 (RSE 19%)
    lvp <- log(1270); label("Apparent peripheral volume V2/F at 70 kg (L)") # Table 4 final 'V 2 / F (L/70 kg)' = 1270 (RSE 13%)

    # Allometric exponents. Results 3.2.1: 'Allometric scaling with an
    # estimated exponent on CL/F and fixed exponent on volume of distribution
    # of the central compartment (1) and volume of distribution of the
    # peripheral compartment (1)'.
    e_wt_cl <- 0.62; label("Allometric exponent of (WT/70) on CL/F (unitless)") # Table 4 final 'Allometric scaling on CL' = 0.62 (RSE 20%); Equation (2)
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V1/F (unitless)") # Results 3.2.1, fixed exponent (1)
    e_wt_vp <- fixed(1); label("Allometric exponent of (WT/70) on V2/F (unitless)") # Results 3.2.1, fixed exponent (1)

    # Covariate effects on CL/F, Equation (2):
    #   CL/F = 36.6 * (weight/70)^0.62
    #          * [(1.0, if CYP3A5*3/*3) or (1.4, if CYP3A5*1/*3 or *1/*1)]
    #          * (Hematocrit/0.29)^-0.6 * (Creatinine/84)^-0.1
    e_cyp3a5_expr_cl <- 1.4; label("CYP3A5 expresser (*1/*1 or *1/*3) multiplier on CL/F") # Table 4 final 'CYP3A5*1/*1 or *1/*3' = 1.4; Equation (2)
    e_hct_cl <- -0.60; label("Haematocrit power exponent on CL/F, centred at 0.29 L/L (unitless)") # Table 4 final 'Hematocrit (L/L)' = -0.60; Equation (2)
    e_creat_cl <- -0.1; label("Serum creatinine power exponent on CL/F, centred at 84 umol/L (unitless)") # Table 4 final 'Creatinine (umol/L)' = -0.1; Equation (2)

    # Diagonal IIV. Table 4 reports IIV as %CV; variances on the log scale are
    # omega^2 = log(1 + CV^2) (same convention as the Andrews 2017 and 2019
    # tacrolimus models in this library):
    #   ka   183%  -> log(1 + 1.83^2)  = 1.469923
    #   CL/F 42.1% -> log(1 + 0.421^2) = 0.163174
    #   V1/F 99.6% -> log(1 + 0.996^2) = 0.689147
    #   V2/F 85.2% -> log(1 + 0.852^2) = 0.545751
    etalka ~ 1.469923 # Table 4 final IIV k a = 183% (RSE 11%) [shrinkage 20%]
    etalcl ~ 0.163174 # Table 4 final IIV CL/F = 42.1% (RSE 10%) [shrinkage 5%]
    etalvc ~ 0.689147 # Table 4 final IIV V1/F = 99.6% (RSE 12%) [shrinkage 10%]
    etalvp ~ 0.545751 # Table 4 final IIV V2/F = 85.2% (RSE 15%) [shrinkage 22%]

    # Combined additive + proportional residual error with separate
    # magnitudes per assay (Results 3.2.1). Table 4 row 'Residual variability
    # Additional' is the additive SD in ng/mL; 'Proportional' is a fraction.
    addSd_immuno <- 1.27; label("Additive residual SD for immunoassay samples (ng/mL)") # Table 4 final additive Immunoassay = 1.27
    addSd_lcms <- 0.87; label("Additive residual SD for LC-MS/MS samples (ng/mL)") # Table 4 final additive LC-MS/MS = 0.87
    propSd_immuno <- 0.11; label("Proportional residual SD for immunoassay samples (fraction)") # Table 4 final proportional Immunoassay = 0.11
    propSd_lcms <- 0.23; label("Proportional residual SD for LC-MS/MS samples (fraction)") # Table 4 final proportional LC-MS/MS = 0.23
  })

  model({
    # CYP3A5 expresser multiplier on CL/F; non-expresser reference = 1.0
    # (Equation 2).
    f_cyp3a5 <- 1 + (e_cyp3a5_expr_cl - 1) * CYP3A5_EXPR

    # Haematocrit (canonical percent scale; 29% = 0.29 L/L) and serum
    # creatinine power effects on CL/F (Equation 2).
    f_hct <- (HCT / 29)^e_hct_cl
    f_creat <- (CREAT / 84)^e_creat_cl

    wt70 <- WT / 70

    tlag <- exp(ltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * wt70^e_wt_cl * f_cyp3a5 * f_hct * f_creat
    vc <- exp(lvc + etalvc) * wt70^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * wt70^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Whole-blood tacrolimus in ng/mL: dose in mg and vc in L give mg/L, so
    # multiply by 1000.
    Cc <- central / vc * 1000

    # Per-sample assay-conditional residual error (IMMUNOASSAY: 1 =
    # immunoassay, 0 = LC-MS/MS).
    addSd <- addSd_immuno * IMMUNOASSAY + addSd_lcms * (1 - IMMUNOASSAY)
    propSd <- propSd_immuno * IMMUNOASSAY + propSd_lcms * (1 - IMMUNOASSAY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
