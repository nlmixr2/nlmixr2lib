Kim_2024_teicoplanin <- function() {
  description <- "Three-compartment population PK model for intravenous teicoplanin in healthy Korean adults with normal renal function, with a power effect of BSA-adjusted CKD-EPI creatinine eGFR on total clearance and a power effect of body weight on the second peripheral volume (Kim 2024)"
  reference <- "Kim YK, Jo KM, Lee JH, Jang JH, Choe EJ, Kang G, Zang DY, Lee DH. Beyond One-Size-Fits-All: Tailoring Teicoplanin Regimens for Normal Renal Function Patients Using Population Pharmacokinetics and Monte Carlo Simulation. Pharmaceutics. 2024;16(4):499. doi:10.3390/pharmaceutics16040499"
  vignette <- "Kim_2024_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kim 2024 Section 3.2, which names
  # V1 (central), V2 (first peripheral) and V3 (second peripheral) as the
  # three distribution volumes of the plasma teicoplanin model.
  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the 2021 CKD-EPI creatinine equation, de-normalized to the individual's body surface area (raw mL/min, NOT per 1.73 m^2)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source symbol CE (Kim 2024 Table 2 footnote: 'estimated glomerular filtration rate calculated using the CKD-EPI equation based on creatinine levels, adjusted for body surface area'). Built in two steps per Kim 2024 Table 1 footnotes c and e: first the creatinine-only CKD-EPI equation gives eGFR in mL/min/1.73 m^2 (female: 142 * min(CR/0.7,1)^-0.241 * max(CR/0.7,1)^-1.200 * 0.9938^Age * 1.012; male: 142 * min(CR/0.9,1)^-0.302 * max(CR/0.9,1)^-1.200 * 0.9938^Age, with CR in mg/dL), then it is multiplied by BSA/1.73 m^2 to give an absolute filtration rate in mL/min. Stored under the canonical CRCL column, which pools creatinine-based and tracer-measured filtration estimates; this model uses the RAW, de-normalized mL/min variant following the precedent of Delattre_2010_amikacin.R, Georges_2009_ceftazidime.R and Chen_2023_nemonoxacin.R. Supplying a BSA-normalized value (mL/min/1.73 m^2) would misstate the covariate, because the 105.27 mL/min reference is on the de-normalized scale. Population values (Kim 2024 Table 1, 'Adjusted eGFR by CKD-EPICR for BSA'): mean 103 mL/min (CV 15.5%), median 105 (IQR 96.7-113). The paper explicitly tested and rejected Cockcroft-Gault CLCR, MDRD eGFR, and the combined creatinine-cystatin C CKD-EPI eGFR in favour of this covariate (Kim 2024 Section 4). Effect form: power scaling (CRCL / 105.27)^e_crcl_cl on CL.",
      source_name        = "CE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed); single-dose study in healthy volunteers, so no within-subject weight change. Population values (Kim 2024 Table 1): mean 64.8 kg (CV 19.9%), median 67.9 (IQR 51.3-73.4). Effect form: power scaling (WT / 67.85)^e_wt_vp2 on the SECOND peripheral volume V3 only. Kim 2024 did NOT carry an allometric term on CL, Vc, Vp, Q or Q2 - the stepwise covariate search retained weight on V3 alone (dOFV 16.6 on removal; Section 3.2). The 67.85 kg reference is the cohort median weight (Table 1 reports 67.9 kg to three significant figures).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    n_observations = 96L,
    age_range      = "19-55 years (protocol inclusion criterion)",
    age_median     = "32.0 years (IQR 30.0-40.3)",
    weight_median  = "67.9 kg (IQR 51.3-73.4)",
    height_median  = "164 cm (IQR 158-169)",
    bsa_median     = "1.77 m^2 (IQR 1.52-1.85)",
    sex_female_pct = 50,
    race_ethnicity = "Korean (single-centre study in the Republic of Korea); not tabulated by the source",
    disease_state  = "Healthy adults with normal renal function; no congenital or chronic disease, screened by medical history, vital signs, physical examination, haematology, blood chemistry, urinalysis and infectious serology",
    dose_range     = "200 mg teicoplanin in 100 mL normal saline as a single 30-minute intravenous infusion",
    regions        = "Republic of Korea (Hallym University Sacred Heart Hospital, Anyang)",
    renal_function = "Normal: BSA-adjusted CKD-EPI creatinine eGFR median 105 mL/min (IQR 96.7-113); serum creatinine median 0.875 mg/dL; cystatin C median 0.760 mg/dL",
    notes          = "Prospective single-dose study conducted July-September 2023 (IRB 2023-05-013). Eight samples per subject at 33, 36, 45 and 90 min and 4, 8, 48-120 and 168-240 h after the start of infusion; 96 plasma concentrations in total. Assay: HPLC-MS/MS (Shimadzu LC-40, Phenomenex Gemini C18, SCIEX 4500 QTRAP) with vancomycin internal standard and 1/x^2-weighted calibration. Estimation: NONMEM 7.5 FOCE-I; model evaluation by VPC (1000 replicates) and 2000-sample nonparametric bootstrap via PsN 5.3.1. Final model OFV -209.055 (one-compartment 305.997, two-compartment -65.658, three-compartment base -186.520). See Kim 2024 Table 1 for baseline demographics."
  )

  ini({
    # Structural parameters - Kim 2024 Table 2, "Structural model" block.
    # NONMEM parameterisation: CL = theta1 * (CE/105.27)^theta2, V1 = theta3,
    # Q2 = theta4, V2 = theta5, Q3 = theta6, V3 = theta7 * (WT/67.85)^theta8.
    # NONMEM V1/V2/V3 map to canonical vc/vp/vp2, and Q2/Q3 to q/q2.
    lcl  <- log(0.693); label("Total clearance at CRCL = 105.27 mL/min (CL, L/h)")                    # Kim 2024 Table 2 theta1 = 0.693 (RSE 2.97%; bootstrap 0.693, 95% CI 0.653-0.74)
    lvc  <- log(3.96);  label("Central volume of distribution (V1, L)")                               # Kim 2024 Table 2 theta3 = 3.96 (RSE 8.41%; bootstrap 3.97, 95% CI 3.15-4.62)
    lq   <- log(4.45);  label("Intercompartmental clearance central <-> peripheral1 (Q2, L/h)")       # Kim 2024 Table 2 theta4 = 4.45 (RSE 11.6%; bootstrap 4.45, 95% CI 3.63-5.86)
    lvp  <- log(8.24);  label("First peripheral volume of distribution (V2, L)")                      # Kim 2024 Table 2 theta5 = 8.24 (RSE 8.32%; bootstrap 8.33, 95% CI 7.07-9.85)
    lq2  <- log(1.76);  label("Intercompartmental clearance central <-> peripheral2 (Q3, L/h)")       # Kim 2024 Table 2 theta6 = 1.76 (RSE 9.7%; bootstrap 1.75, 95% CI 1.44-2.13)
    lvp2 <- log(69.8);  label("Second peripheral volume of distribution at WT = 67.85 kg (V3, L)")    # Kim 2024 Table 2 theta7 = 69.8 (RSE 8.74%; bootstrap 69.7, 95% CI 55.8-82.6)

    # Covariate effects - both are estimated power exponents (not fixed
    # allometric constants); each carries an RSE and a bootstrap CI in Table 2,
    # so neither is wrapped in fixed().
    e_crcl_cl <- 0.785; label("Power exponent on (CRCL / 105.27) for CL (unitless)")  # Kim 2024 Table 2 theta2 = 0.785 (RSE 16.2%; bootstrap 0.789, 95% CI 0.422-1.16)
    e_wt_vp2  <- 1.73;  label("Power exponent on (WT / 67.85) for V3 (unitless)")     # Kim 2024 Table 2 theta8 = 1.73 (RSE 22.6%; bootstrap 1.73, 95% CI 0.67-2.44)

    # Inter-individual variability. Kim 2024 Table 2 reports each IIV as a
    # percentage under the heading "Interindividual variability"; NONMEM
    # log-normal omegas are converted to the internal variance scale with
    # omega^2 = log(1 + CV^2).
    #
    # The BSV on CL, V1, V2 and V3 was FIXED by the authors because the
    # corresponding relative standard errors exceeded 25% (Kim 2024 Section 3.2
    # and Section 4 limitation 2; the "f" superscript in Table 2). Those four
    # therefore carry fixed(); Q2 and Q3 were estimated (RSE 20.2% and 18.9%,
    # with bootstrap CIs) and do not.
    etalcl  ~ fixed(0.0077667) # Kim 2024 Table 2 CL IIV 8.83%; log(1 + 0.0883^2)
    etalvc  ~ fixed(0.0550978) # Kim 2024 Table 2 V1 IIV 23.8%; log(1 + 0.238^2)
    etalq   ~ 0.1015895        # Kim 2024 Table 2 Q2 IIV 32.7% (RSE 20.2%; bootstrap 30.7, 95% CI 14.2-42.2); log(1 + 0.327^2)
    etalvp  ~ fixed(0.0555492) # Kim 2024 Table 2 V2 IIV 23.9%; log(1 + 0.239^2)
    etalq2  ~ 0.0917584        # Kim 2024 Table 2 Q3 IIV 31.0% (RSE 18.9%; bootstrap 29.6, 95% CI 16.7-41.4); log(1 + 0.310^2)
    etalvp2 ~ fixed(0.0056691) # Kim 2024 Table 2 V3 IIV 7.54%; log(1 + 0.0754^2)

    # Residual error. Table 2 reports a single "Proportional error (%)" of
    # 6.33% with no additive component, so the error model is purely
    # proportional on the linear concentration scale.
    propSd <- 0.0633; label("Proportional residual error (fraction)") # Kim 2024 Table 2 proportional error 6.33% (RSE 13.1%; bootstrap 6.22, 95% CI 4.63-7.94)
  })

  model({
    # Individual PK parameters. Both covariate effects are power functions
    # centred on the cohort medians reported in Kim 2024 Table 1
    # (BSA-adjusted CKD-EPI creatinine eGFR 105 mL/min, weight 67.9 kg);
    # the paper's own centring constants 105.27 and 67.85 are used verbatim.
    cl  <- exp(lcl + etalcl) * (CRCL / 105.27)^e_crcl_cl
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq)
    vp  <- exp(lvp + etalvp)
    q2  <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2) * (WT / 67.85)^e_wt_vp2

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment intravenous model. Teicoplanin was given only as a
    # 30-minute IV infusion into the central compartment (Kim 2024 Section 2.2),
    # so there is no depot and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Dose in mg, vc in L -> mg/L, the unit used throughout Kim 2024.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
