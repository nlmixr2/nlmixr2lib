Le_2025_meropenem <- function() {
  description <- "Two-compartment IV population PK model for meropenem in critically ill children with and without sepsis (4 weeks to <21 years) given 20 mg/kg every 8 h as 30-min infusions (Le 2025). Central and peripheral volumes scale linearly with body weight and clearance and intercompartmental clearance allometrically (exponent 0.75), all with reference weight 70 kg. The central volume additionally increases linearly with serum cystatin C (normalised to 0.45 mg/L) and clearance follows a power function of modified-Schwartz eGFR (reference 140 mL/min/1.73 m^2). A single random effect is shared by the central and peripheral volumes; residual error is proportional."
  reference <- "Le J, Huynh J, Vo B, Mai A, Mak RH, Momper JD, Capparelli EV, Harvey H, Avedissian S, Bradley E, Sitapati A, Singh K, Bradley JS. Variability in Meropenem Distribution and Clearance in Children with Sepsis: Population-Based Pharmacokinetics with Assessment of Renal Biomarkers. Clin Pharmacokinet. 2025;64(5):769-777. doi:10.1007/s40262-025-01495-3"
  vignette <- "Le_2025_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Meropenem is given as an IV infusion, so the dose enters `central`
  # directly. The paper calls the measured concentrations "serum meropenem
  # concentrations" throughout but describes the sample handling as "plasma
  # frozen at minus 70 C" with calibrators prepared in human plasma, so the
  # matrix is not settled by the source and the entry is left unverified.
  # `peripheral1` is a mathematical distribution compartment.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Le 2025 Table 1: median 36.3 kg (range 7.2-98.0). Reference weight 70 kg in all four structural equations of the Table 2 footnote: (Wt/70) on V1 and V2, (Wt/70)^0.75 on CL and Q.",
      source_name = "Wt"
    ),
    CYSC = list(
      description = "Serum cystatin C (research ELISA assay)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source analysis (nearest value carried forward or backward). The paper reports cystatin C in ng/mL (Online Resource 3: median 451.1, range 178.3-1824.1 ng/mL), so CYSC = SCys (ng/mL) / 1000. The Table 2 footnote divides SCys by 0.45, which is the cohort median expressed in mg/L (0.451 mg/L); read literally in ng/mL the term would be about 1000 and V1 several hundred times too large, so the normaliser is taken as 0.45 mg/L. Measured by LS Bio ELISA (catalog LS-F436), which reads lower than the clinical immunoturbidimetric assays behind the adult reference range of about 0.6-1.1 mg/L.",
      source_name = "SCys"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the modified (bedside) Schwartz equation, eGFR (mL/min/1.73 m^2) = 0.413 * height (cm) / serum creatinine (mg/dL)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source analysis. Online Resource 1 gives the modified Schwartz equation for children <18 years and CKD-EPI (creatinine) for the two participants aged >=18 years. Online Resource 3: median 138.92 (range 22.85-364.82) mL/min/1.73 m^2. Reference value 140 mL/min/1.73 m^2 from the Table 2 footnote clearance equation.",
      source_name = "eGFR-MS"
    )
  )

  # Screened in the univariable step (Online Resource 4) but not retained in
  # the final model; documentation only, not referenced in model().
  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Online Resource 4: dOFV -5.57 on CL and -1.23 on V1. Enters the final model only through eGFR-MS.",
      source_name = "SCr"
    ),
    ALB = list(
      description = "Serum albumin (lowest value)",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Online Resource 4: dOFV -5.28 on V1 for the low value; described in Results as a significant univariable covariate for V1 but not retained after stepwise selection.",
      source_name = "Albumin - Low Value"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Online Resource 4: dOFV -0.05 on CL and -0.49 on V1 (continuous); an age-group categorisation gave -9.04 on CL. Not retained.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Online Resource 4: dOFV -0.19 on CL and -1.56 on V1. Not retained.",
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    n_concentrations = 309L,
    age_range = "0.6-19.6 years (eligibility 4 weeks to <21 years)",
    age_median = "11.8 years",
    weight_range = "7.2-98.0 kg",
    weight_median = "36.3 kg",
    sex_female_pct = 33.3,
    race_ethnicity = c(`Hispanic (ethnicity)` = 74, Caucasian = 70, `Other/unknown` = 30),
    disease_state = "Children admitted to the PICU for >=72 h: 19 cases with sepsis or septic shock requiring vasopressors and 8 matched controls receiving meropenem for other indications (88% appendicitis) without sepsis or hypotension. PELOD score median 2 (range 0-11).",
    dose_range = "Meropenem 20 mg/kg IV every 8 h over 30 min (maximum 2000 mg per dose), standard of care",
    regions = "United States (Rady Children's Hospital, San Diego, CA)",
    renal_function = "Modified-Schwartz eGFR median 138.9 (range 22.9-364.8) mL/min/1.73 m^2, spanning acute kidney injury to augmented renal clearance; serum creatinine median 0.33 (0.09-2.57) mg/dL; serum cystatin C median 451.1 (178.3-1824.1) ng/mL",
    notes = "Prospective study, March 2019 to January 2023 (Le 2025 Methods, Table 1, Online Resource 3). Intensive sampling over 8 h after a dose on days 1-3 for cases (pre-dose, 0.5, 1, 2, 4 and 8 h after the start of infusion on day 1; pre-dose, 1, 2, 4 and 8 h on days 2 and 3) and on day 1 only for controls; 271 concentrations from cases and 38 from controls. NONMEM 7.3, FOCE with interaction; 1000-sample bootstrap. Serum neutrophil gelatinase-associated lipocalin (SNgal, median 180.7 ng/mL) was also screened (Online Resource 4: dOFV -6.89 on CL, -3.57 on V1) and not retained."
  )

  ini({
    # Structural parameters (Le 2025 Table 2), for a 70-kg subject with
    # serum cystatin C 0.45 mg/L (V1 multiplier 1 + 0.449 = 1.449) and
    # eGFR 140 mL/min/1.73 m^2.
    lvc <- log(9.74); label("Central volume at WT = 70 kg before the cystatin C term (L)") # Table 2 thetaV1 = 9.74 (SE 3.90; bootstrap 95% CI 5.07-15.22)
    lvp <- log(4.53); label("Peripheral volume at WT = 70 kg (L)") # Table 2 thetaV2 = 4.53 (SE 1.23)
    lcl <- log(9.72); label("Clearance at WT = 70 kg and eGFR = 140 mL/min/1.73 m^2 (L/h)") # Table 2 thetaCL = 9.72 (SE 2.37)
    lq <- log(1.66); label("Intercompartmental clearance at WT = 70 kg (L/h)") # Table 2 thetaQ = 1.66 (SE 0.400)

    # Size scaling. Methods: weight was applied in the base model "with
    # allometric scaling of 0.75 for CLMERO and isometric for VMERO"; the
    # exponents carry no estimate in Table 2, so they are fixed.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless)") # Table 2 footnote: (Wt/70)^0.75 on CL and Q
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/70) for V1 and V2 (unitless)") # Table 2 footnote: (Wt/70) on V1 and V2

    # Covariate effects (estimated).
    e_cysc_vc <- 0.449; label("Linear slope of (CYSC/0.45) on the central volume (unitless)") # Table 2 thetaV1-SCys = 0.449 (SE 0.428); footnote V1 = thetaV1 * (Wt/70) * (1 + thetaV1-SCYS * SCys/0.45)
    e_crcl_cl <- 0.169; label("Power exponent of (CRCL/140) on clearance (unitless)") # Table 2 thetaCL-eGFR-MS = 0.169 (SE 0.0967); footnote CL = thetaCL * (Wt/70)^0.75 * (eGFR/140)^thetaCL-eGFR-MS

    # Between-subject variability. Table 2 prints the raw NONMEM OMEGA
    # (variances): the bootstrap 95% CI relative width of the eta V1/V2 row,
    # (0.464 - 0.153) / 0.287 = 1.08, matches the variance prediction
    # 3.92 * sqrt(2/27) = 1.07 and not the SD prediction 0.53. One eta is
    # shared by V1 and V2 ('The inter-subject variability for V1-MERO and
    # V2-MERO was assumed to be equal'); Online Resource 3 reports identical
    # quartile-wise case-vs-control p-values for V1 and V2, which a single
    # shared eta produces and two independent etas would not.
    etalcl ~ 0.487 # Table 2 row 'etaCL' = 0.487 (SE 0.103)
    etalvc ~ 0.306 # Table 2 row 'etaV1, etaV2' = 0.306 (SE 0.0739); shared by vc and vp

    # Residual error. Table 2 prints a single epsilon. FOCE with ETA-EPS
    # interaction was used (Methods) and the Figure 1a observed-vs-IPRED
    # scatter widens with concentration, so the error is proportional; the
    # value is read on the same raw NONMEM (variance) scale as the OMEGA rows.
    propSd <- 0.3647; label("Proportional residual error (fraction)") # Table 2 epsilon = 0.133 (SE 0.0205), variance; sqrt(0.133) = 0.3647
  })
  model({
    # Individual parameters (Le 2025 Table 2 footnote):
    #   V1 = thetaV1 * (Wt/70) * (1 + thetaV1-SCys * SCys/0.45)
    #   V2 = thetaV2 * (Wt/70)
    #   CL = thetaCL * (Wt/70)^0.75 * (eGFR/140)^thetaCL-eGFR-MS
    #   Q  = thetaQ  * (Wt/70)^0.75
    cysc_vc <- 1 + e_cysc_vc * CYSC / 0.45
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * cysc_vc
    vp <- exp(lvp + etalvc) * (WT / 70)^e_wt_vc_vp
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (CRCL / 140)^e_crcl_cl
    q <- exp(lq) * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L, so central / vc is mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
