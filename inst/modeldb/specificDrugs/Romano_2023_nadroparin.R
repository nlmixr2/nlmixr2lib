Romano_2023_nadroparin <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and first-order elimination for nadroparin (a low-molecular-weight heparin) given as thromboprophylaxis to critically ill COVID-19 patients in the intensive care unit (Romano 2023). The model is fitted to plasma anti-Xa activity (IU/mL), so clearance and volume are apparent anti-Xa quantities (CL/F, Vd/F). Apparent clearance (2230 mL/h at the reference patient) carries inter-individual variability and five covariate effects: power functions of C-reactive protein (reference 100 mg/L), D-dimer (reference 10 mg/L) and CKD-EPI 2009 eGFR (reference 80 mL/min/1.73 m^2, computed inside the model from serum creatinine, age and sex), and multiplicative factors of 0.749 for vasopressor use and 0.775 for corticosteroid use. Residual error is combined additive (0.0859 IU/mL) plus proportional (20%)."
  reference <- paste(
    "Romano LGR, Hunfeld NGM, Kruip MJHA, Endeman H, Preijers T.",
    "Population pharmacokinetics of nadroparin for thromboprophylaxis in COVID-19",
    "intensive care unit patients.",
    "Br J Clin Pharmacol. 2023;89(5):1617-1628.",
    "doi:10.1111/bcp.15634",
    sep = " "
  )
  vignette <- "Romano_2023_nadroparin"

  units <- list(
    time = "h",
    dosing = "IU",
    concentration = "IU/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "nadroparin", units = "IU", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nadroparin", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRP = list(
      description = "C-reactive protein, standard assay",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source data (Romano 2023 Table 1 footnote: continuous characteristics",
        "other than age, height and SOFA score were measured repeatedly within each patient).",
        "Enters CL/F as (CRP / 100)^0.182; the reference 100 mg/L is the value the Table 2",
        "footnote names for the typical patient, not the cohort median (Table 1 median 61.0 mg/L,",
        "IQR 32.5-101; modelling cohort median 76.0 mg/L). The control stream (Appendix SA)",
        "substitutes the reference value 100 mg/L when CRP is recorded as 0, i.e. missing;",
        "that imputation is a data-handling rule and is not reproduced here.",
        sep = " "
      ),
      source_name = "CRP"
    ),
    DDIMER = list(
      description = "Plasma D-dimer concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source records D-dimer in mg/L (Romano 2023 Table 1, median 1.24 mg/L, IQR",
        "0.85-2.18; Figure 5 x-axis 'D-Dimer (mg/L)' spanning 0-10). The control stream",
        "column DDMR enters CL/F as (DDMR / 10)^0.117 with DDMR in mg/L. The canonical DDIMER",
        "unit is ng/mL, so the model divides by 10000 ng/mL (= 10 mg/L); 1 mg/L = 1000 ng/mL.",
        "The reference is the Table 2 footnote value, which lies at the top of the observed",
        "range rather than at the cohort median. The assay reporting basis (fibrinogen-",
        "equivalent units vs D-dimer units) is not stated. Time-varying in the source data.",
        sep = " "
      ),
      source_name = "DDMR (control stream); 'D-dimer' (Table 1, Table 2)"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only to compute the CKD-EPI 2009 eGFR inside model(), exactly as Appendix SA does:",
        "eGFR = 141 * min(SCr/kappa, 1)^alpha * max(SCr/kappa, 1)^-1.209 * 0.993^AGE,",
        "times 1.018 for women, with kappa = 61.9 umol/L (women) or 79.6 umol/L (men) and",
        "alpha = -0.329 (women) or -0.411 (men). The kappa values are 0.7 and 0.9 mg/dL",
        "expressed in umol/L, so SCr must be supplied in umol/L. Time-varying in the source data.",
        sep = " "
      ),
      source_name = "CREAT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only in the CKD-EPI 2009 eGFR term 0.993^AGE (Appendix SA). Table 1 median 63.0 years (IQR 53.0-69.7).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Used only in the CKD-EPI 2009 eGFR term (Appendix SA selects kappa and alpha by sex and",
        "multiplies the female value by 1.018). The control stream column GEN is coded 1 = male",
        "('GEN=1 is male'), so SEXF = 1 - GEN. Romano 2023 Table 1 does not tabulate sex.",
        sep = " "
      ),
      source_name = "GEN (control stream; 1 = male)"
    ),
    CONMED_INOTROPE = list(
      description = "Vasopressor or inodilator use; 1 = in use, 0 = not in use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no vasopressor)",
      notes = paste(
        "Romano 2023 Methods 2.1: vasopressor (including inodilator) use comprised",
        "norepinephrine, epinephrine or enoximone; the Discussion notes doses usually did not",
        "exceed about 0.1 ug/kg/min. Enters CL/F as 0.749^CONMED_INOTROPE (25.1% lower CL/F",
        "with use). 89% of the modelling cohort used vasopressors at the start (Table 1).",
        "Time-varying in the source data (control stream column VASO).",
        sep = " "
      ),
      source_name = "VASO (control stream); 'Vasopressors' (Table 2 equation)"
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid use; 1 = in use, 0 = not in use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no corticosteroid)",
      notes = paste(
        "Romano 2023 Methods 2.1: (methyl)prednisolone or dexamethasone. Enters CL/F as",
        "0.775^CONMED_STEROID (22.5% lower CL/F with use). 65% of the modelling cohort used",
        "corticosteroids at the start (Table 1). Time-varying in the source data (control",
        "stream column CORT).",
        sep = " "
      ),
      source_name = "CORT (control stream); 'Corticosteroids' (Table 2 equation)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a renal-function-related size descriptor (Romano 2023 Methods 2.5) and not retained; no estimate is reported. Table 1 median 89.0 kg (IQR 78.5-97.3).",
      source_name = "Bodyweight (Table 1)"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5, 'ideal/lean body weight') and not retained; no estimate is reported. Table 1 median 60.7 kg (IQR 50.0-67.3).",
      source_name = "Lean body mass calculated (Table 1)"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5, 'length') and not retained; no estimate is reported. Table 1 median 1.75 m (IQR 1.65-1.80).",
      source_name = "Height (Table 1)"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5) and not retained; no estimate is reported. Table 1 median 29.6 kg/m^2 (IQR 27.5-32.3).",
      source_name = "Body mass index (Table 1)"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a hepatic-function marker (Methods 2.5) and not retained; no estimate is reported. Table 1 median 20.0 g/L (IQR 18.5-23.0), modelling cohort only.",
      source_name = "Albumin (Table 1)"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a hepatic-function marker (Methods 2.5) and not retained; no estimate is reported. Table 1 median 49.5 U/L (IQR 32.0-81.0), modelling cohort only.",
      source_name = "Alanine transaminase (Table 1)"
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a coagulopathy marker (Methods 2.5, 'thrombocytes') and not retained; no estimate is reported. Table 1 median 345 x 10^9/L (IQR 265-417), modelling cohort only.",
      source_name = "Thrombocytes (Table 1)"
    ),
    LACT = list(
      description = "Arterial lactate",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods 2.5) and not retained; no estimate is reported. Table 1 median 1.15 mmol/L (IQR 1.00-1.40), modelling cohort only.",
      source_name = "Arterial lactate (Table 1)"
    ),
    BLOOD_GROUP_O = list(
      description = "Blood group O; 1 = group O, 0 = other",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-O)",
      notes = "Screened as a coagulopathy marker (Methods 2.5, 'blood group') and not retained; no estimate is reported. Table 1: 21 of 90 patients (23%) were group O; blood group was unknown for two patients.",
      source_name = "Blood group O (Table 1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 1L,
    n_observations = "280 anti-Xa samples in the modelling cohort (21 below the 0.1 IU/mL limit of quantification, handled with the M3 method); a separate validation cohort of 25 patients contributed 147 samples (Romano 2023 Table 1)",
    age_range = "modelling cohort median 62.0 years (IQR 53.0-70.0)",
    age_median = "62.0 years (modelling cohort)",
    weight_range = "modelling cohort IQR 78.0-100 kg",
    weight_median = "89.0 kg (modelling cohort)",
    race_ethnicity = "Not reported (single centre in Rotterdam, the Netherlands)",
    disease_state = "Critically ill adults with PCR-confirmed SARS-CoV-2 infection admitted to the intensive care unit; 87% mechanically ventilated and 9% on renal replacement therapy at first anti-Xa sampling (modelling cohort, Table 1)",
    dose_range = "Nadroparin 5700 IU subcutaneously twice daily per the local ICU protocol (dosing frequency lowered when CKD-EPI eGFR < 30 mL/min)",
    regions = "The Netherlands (Erasmus University Medical Center, Rotterdam; single centre)",
    renal_function = "CKD-EPI eGFR modelling cohort median 94.4 mL/min/1.73 m^2 (IQR 78.6-111)",
    inflammation = "CRP modelling cohort median 76.0 mg/L (IQR 35.5-124); D-dimer median 1.29 mg/L (IQR 0.87-2.19)",
    co_medication = "Vasopressors at the start 89%, corticosteroids at the start 65% (modelling cohort)",
    notes = paste(
      "Retrospective single-centre observational cohort of ICU admissions between 1 March 2020",
      "and 30 January 2021 (Romano 2023 Methods 2.1). Anti-Xa activity was targeted 4 h after",
      "a dose, twice weekly (Sysmex CS-5100; LLOQ 0.1 IU/mL, ULOQ 4.0 IU/mL). The model was",
      "fitted in NONMEM 7.4 with Laplacian conditional estimation and the M3 method for",
      "below-quantification samples.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Romano 2023 Table 2 'Final model' column; the
    # same values are the final $THETA initials in Appendix SA. Dose (IU) /
    # volume (mL) gives anti-Xa activity directly in IU/mL. Typical values are
    # for a patient with CRP 100 mg/L, D-dimer 10 mg/L and eGFR 80 mL/min/1.73
    # m^2 on neither corticosteroids nor vasopressors (Table 2 footnote).
    lka <- log(0.276); label("Absorption rate constant ka (1/h)") # Romano 2023 Table 2 final model ka = 0.276 1/h (RSE 28%); Appendix SA THETA(3)
    lcl <- log(2230); label("Apparent anti-Xa clearance CL/F at the reference patient (mL/h)") # Romano 2023 Table 2 final model CL = 2230 mL/h (RSE 12%); Appendix SA THETA(4)
    lvc <- log(11000); label("Apparent anti-Xa volume of distribution Vd/F (mL)") # Romano 2023 Table 2 final model Vd = 11 000 mL (RSE 22%); Appendix SA THETA(5)

    # Covariate effects on CL/F (Table 2 equation for CL_i; Appendix SA
    # COV_CL, COV_DDMR, COV_VASO, COV_CL2, COV_CORT).
    e_crp_cl <- 0.182; label("Power exponent of (CRP / 100 mg/L) on CL/F (unitless)") # Romano 2023 Table 2 'CL-CRP' 0.182 (RSE 21%); Appendix SA THETA(6)
    e_ddimer_cl <- 0.117; label("Power exponent of (D-dimer / 10 mg/L) on CL/F (unitless)") # Romano 2023 Table 2 'CL-D-dimer' 0.117 (RSE 56%); Appendix SA THETA(7)
    e_inotrope_cl <- 0.749; label("Multiplicative factor on CL/F with vasopressor use (theta^CONMED_INOTROPE)") # Romano 2023 Table 2 'CL-use of vasopressors' 0.749 (RSE 8%); Appendix SA THETA(8)
    e_crcl_cl <- 0.368; label("Power exponent of (CKD-EPI eGFR / 80 mL/min/1.73 m^2) on CL/F (unitless)") # Romano 2023 Table 2 'CL-GFR CKD-EPI' 0.368 (RSE 33%); Appendix SA THETA(9)
    e_steroid_cl <- 0.775; label("Multiplicative factor on CL/F with corticosteroid use (theta^CONMED_STEROID)") # Romano 2023 Table 2 'CL-use of corticosteroids' 0.775 (RSE 11%); Appendix SA THETA(10)

    # IIV on CL/F only (Results 3.2). Appendix SA $OMEGA 0.115 is the variance;
    # it matches Table 2's 34.9% CV as log(1 + 0.349^2) = 0.1149.
    etalcl ~ 0.115 # Romano 2023 Appendix SA $OMEGA 0.115; Table 2 IIV on CL 34.9% CV (RSE 21%, shrinkage 14%)

    # Residual error. Appendix SA: W = SQRT(THETA(1)^2 + THETA(2)^2 * IPRED^2),
    # Y = IPRED + W * EPS(1) with $SIGMA 1 FIX, so both THETAs are standard
    # deviations combined in quadrature (nlmixr2's default combined2 form).
    addSd <- 0.0859; label("Additive residual error SD (IU/mL)") # Romano 2023 Appendix SA THETA(1) = 0.0859; Table 2 prints 0.086 IU/mL (RSE 13%)
    propSd <- 0.2; label("Proportional residual error SD (fraction)") # Romano 2023 Appendix SA THETA(2) = 0.2; Table 2 'Proportional residual error' 0.2 (RSE 18%)
  })

  model({
    # CKD-EPI 2009 eGFR (mL/min/1.73 m^2) from serum creatinine in umol/L,
    # transcribed from Appendix SA: kappa 61.9 / 79.6 umol/L (= 0.7 / 0.9
    # mg/dL) and alpha -0.329 / -0.411 for women / men, female factor 1.018.
    kappa_scr <- 79.6 * (1 - SEXF) + 61.9 * SEXF
    alpha_scr <- -0.411 * (1 - SEXF) - 0.329 * SEXF
    egfr <- 141 * min(CREAT / kappa_scr, 1)^alpha_scr * max(CREAT / kappa_scr, 1)^(-1.209) * 0.993^AGE * (1 + 0.018 * SEXF)

    # Romano 2023 Table 2 equation:
    # CL_i = 2230 * (CRP/100)^0.182 * (D-dimer/10)^0.117 * (GFR/80)^0.368
    #        * 0.749^Vasopressors * 0.775^Corticosteroids * exp(eta_i),
    # with D-dimer in mg/L; DDIMER here is in ng/mL, so 10 mg/L = 10000 ng/mL.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (CRP / 100)^e_crp_cl *
      (DDIMER / 10000)^e_ddimer_cl *
      (egfr / 80)^e_crcl_cl *
      e_inotrope_cl^CONMED_INOTROPE *
      e_steroid_cl^CONMED_STEROID
    vc <- exp(lvc)

    kel <- cl / vc

    # Appendix SA uses ADVAN5 with DEPOT -> CENTRAL (K12 = KA) and elimination
    # from CENTRAL (K20 = CL/V2); its declared PERIPH compartment has no rate
    # constants, so the model is one-compartment. The IPRED = F + BL baseline
    # term references TVBL, which the control stream never defines, and Results
    # 3.2 states a baseline parameter was not evaluated, so no baseline is added.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in IU and volume in mL, so central / vc is anti-Xa activity in IU/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
