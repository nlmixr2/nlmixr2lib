Yoo_2025_benzylpenicillin <- function() {
  description <- "Two-compartment population PK model for intravenous benzylpenicillin (penicillin G) in healthy adult volunteers of the CHIPS human-challenge trial (Yoo 2025), fitted jointly to a 600 mg IV bolus and to low-dose continuous IV infusions targeting 3-20 ng/mL; clearance is proportional to cystatin-C-based eGFR in absolute mL/min, the central volume is 1.96-fold larger during continuous infusion, and clearance rises by an additive 12.1 L/h during infusion at the 3 ng/mL target"
  reference <- "Yoo O, Salman S, Hla TK, Osowicki J, Page-Sharp M, Marsh JA, Barr R, Azzopardi K, Morici M, Batty KT, Enkel SL, Kado J, Hatchuel L, Fulurija A, McCarthy JS, Snelling T, Steer AC, Carapetis J, Manning L. Population pharmacokinetics of penicillin G: insights into increased clearance at low concentrations to guide development of improved long-acting formulations for syphilis and prevention of rheumatic fever. Antimicrob Agents Chemother. 2025;69(7):e00269-25. doi:10.1128/aac.00269-25"
  vignette <- "Yoo_2025_benzylpenicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Benzylpenicillin was assayed in plasma by a validated LC-MS assay
  # (Methods 'Study design').
  compartmentData <- list(
    central = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling applied a priori with reference 70 kg: exponent 1 on Vc and Vp, 0.75 on Q, none on CL (Methods 'Population pharmacokinetic analysis'; Table 4 footnote). Cohort median 69 kg (range 50-107; Table 1).",
      source_name = "Weight"
    ),
    CRCL = list(
      description = "Cystatin-C-based CKD-EPI 2012 eGFR de-normalized to absolute mL/min",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "NOT BSA-normalized. The source computes eGFRcys with the CKD-EPI cystatin C age-sex equation (2012) in mL/min/1.73 m^2 and converts it to mL/min with the participant's Mosteller BSA, CRCL = eGFRcys * BSA / 1.73, BSA = sqrt(height_cm * weight_kg / 3600) (Methods 'Benzylpenicillin clearance predictor', equation 2). CL enters through the origin as theta * CRCL / 60 (Table 4). Selected over Cockcroft-Gault CLcr, creatinine CKD-EPI 2021 and combined creatinine-cystatin C CKD-EPI 2021 (Table 3). Cohort cystatin C 0.60 mg/L (range 0.39-0.87; Table 1), i.e. eGFRcys of roughly 130-140 mL/min/1.73 m^2 and about 140 mL/min absolute.",
      source_name = "eGFRcys"
    ),
    REGIMEN_CONTINF = list(
      description = "Continuous IV infusion record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IV bolus / short injection: the 600 mg dose given over 3 min)",
      notes = "1 on the low-dose continuous-infusion phase of the CHIPS trial (all five target-Css arms), 0 on the preceding single-dose bolus phase. Gates the central-volume ratio Vc,Linf = Vc * 1.96 (Table 4 'Ratio,Vc'). Time-varying within a subject: each participant contributed bolus records 7-35 days before randomization and infusion records afterwards (Methods 'Study design').",
      source_name = "bolus vs low-dose continuous infusion"
    ),
    DOSE_LOW = list(
      description = "3 ng/mL target-Css infusion arm indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the bolus phase, and the 6, 9, 12 and 20 ng/mL infusion arms)",
      notes = "1 only on the continuous-infusion records of the arm randomized to a target steady-state concentration of 3 ng/mL (infusion rates 54.0-136.7 ug/h; Table 2); 0 on the bolus records of the same participants. Gates the additive clearance step CL_Css3 = CL + 12.1 L/h (Table 4 'Difference,CL,Css3'). The source defines the step by the categorical target Css, having rejected an infusion-rate threshold at 100 ug/h because its estimate was negative (Results 'Population pharmacokinetic analysis'); no target between 3 and 6 ng/mL was studied.",
      source_name = "target Css = 3 ng/mL"
    )
  )

  # Screened in the covariate analysis but not retained (Results 'Population
  # pharmacokinetic analysis': 'additional covariates such as body surface
  # area (BSA), age, and sex showed no significant relationship').
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area (Mosteller)",
      units = "m^2",
      type = "continuous",
      notes = "Median 1.79 m^2 (range 1.49-2.35; Table 1). Not retained as a direct covariate; it enters only through the de-normalization of eGFRcys into CRCL."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Median 25 years (range 18-38; Table 1). Not retained as a direct covariate; it enters only through the CKD-EPI cystatin C equation."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "30 of 72 participants female (42%; Table 1). Not retained as a direct covariate; it enters only through the CKD-EPI cystatin C equation."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 72,
    n_studies = 1,
    n_observations = 893,
    age_range = "18-38 years",
    age_median = "25 years",
    weight_range = "50-107 kg",
    weight_median = "69 kg",
    sex_female_pct = 42,
    race_ethnicity = "Not reported in the source.",
    disease_state = "Healthy adult volunteers (males and non-pregnant, non-lactating females aged 18-40 years without risk factors for severe Streptococcus pyogenes disease) of the CHIPS controlled human pharyngitis challenge trial.",
    renal_function = "Normal. Serum creatinine median 70 umol/L (range 44-101); cystatin C median 0.59 mg/L (range 0.39-0.87; Table 1).",
    dose_range = "Single 600 mg IV benzylpenicillin bolus over 3 min in 10 mL (all 72 participants), then, 7-35 days later, a loading dose plus continuous IV infusion at 46.9-1133.6 ug/h (20 mL/h) individualized to a target steady-state concentration of 3, 6, 9, 12 or 20 ng/mL (47 participants).",
    regions = "Perth, Western Australia, Australia",
    notes = "504 bolus samples (0.25, 0.5, 1, 2, 3, 4 and 6 h post-dose) and 389 continuous-infusion samples (every 12 h). Table 2 reports 10, 9, 10, 9 and 9 infusion participants for the 3, 6, 9, 12 and 20 ng/mL arms (47 total); the Results text lists 9, 9, 8, 8 and 9. NONMEM 7.5.1, FOCE-I, on log-transformed concentrations. Trial conducted August 2022 to July 2023."
  )

  ini({
    # Clearance through the origin in absolute eGFRcys: CL = theta * CRCL / 60
    # (Table 4 row 'CL, l.h-1 = theta x eGFRcys/60'). No weight effect on CL.
    lcl <- log(13.3); label("Clearance at eGFRcys = 60 mL/min (L/h)") # Table 4 theta = 13.3 (RSE 2%; bootstrap 13.31 [12.77-13.92])
    lvc <- log(14.6); label("Central volume of distribution during the bolus phase, 70 kg (L)") # Table 4 VC = 14.6 L/70 kg (RSE 3%)
    lq <- log(5.34); label("Intercompartmental clearance, 70 kg (L/h)") # Table 4 Q = 5.34 L/h/70 kg (RSE 5%)
    lvp <- log(5.77); label("Peripheral volume of distribution, 70 kg (L)") # Table 4 Vp = 5.77 L/70 kg (RSE 4%)

    # Allometric exponents fixed a priori (Methods: 'Allometric scaling based
    # on body size was applied a priori, with exponents of 0.75 for Q and 1
    # for volume parameters').
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc (unitless)") # Methods 'Population pharmacokinetic analysis'
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on Vp (unitless)") # Methods 'Population pharmacokinetic analysis'
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)") # Methods 'Population pharmacokinetic analysis'

    e_regimen_continf_vc <- 1.96; label("Ratio of Vc during low-dose continuous infusion to Vc during the bolus (unitless)") # Table 4 Ratio,Vc = 1.96 (RSE 4%); Vc,Linf = Vc x Ratio,Vc
    e_dose_low_cl <- 12.1; label("Additive increase in CL during infusion at the 3 ng/mL target Css (L/h)") # Table 4 Difference,CL,Css3 = 12.1 L/h (RSE 19%); CL_Css3 = CL + Difference

    # IIV. Table 4 footnote: 'inter-individual variability and residual
    # error are presented as 100% x variability estimate'; read as the SD of
    # the log-normal eta (CV%), so variance = (value/100)^2. See the vignette
    # Assumptions section for the reading and its check against Table 2.
    etalcl ~ 0.144^2 # Table 4 IIV CL = 14.4% [shrinkage 6%]
    etalvc ~ 0.169^2 # Table 4 IIV VC = 16.9% [shrinkage 22%]
    etalvp ~ 0.137^2 # Table 4 IIV Vp = 13.7% [shrinkage 15%]
    etae_dose_low_cl ~ 0.429^2 # Table 4 IIV Difference,CL,Css3 = 42.9% [shrinkage 69%]

    # Log-transformed data with an additive error model (Methods; Results)
    # is proportional error on the linear scale; Table 4 labels the row
    # 'Proportional error (%)'.
    propSd <- 0.165; label("Proportional residual error (fraction)") # Table 4 proportional error = 16.5% [shrinkage 10%] (RSE 1%)
  })
  model({
    # Absolute-mL/min eGFRcys through the origin, reference 60 mL/min.
    cl_base <- exp(lcl + etalcl) * (CRCL / 60)
    # Step increase in clearance on the 3 ng/mL infusion records only;
    # log-normal IIV on the size of the step.
    cl_step <- e_dose_low_cl * exp(etae_dose_low_cl)
    cl <- cl_base + cl_step * DOSE_LOW
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * e_regimen_continf_vc^REGIMEN_CONTINF
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
