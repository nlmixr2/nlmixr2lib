Alihodzic_2022_ciprofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous ciprofloxacin in",
    "critically ill adults on veno-arterial or veno-venous extracorporeal",
    "membrane oxygenation (ECMO) (Alihodzic 2022), with zero-order infusion",
    "into the central compartment, first-order elimination, and linear",
    "centred effects of plasma bicarbonate and of CKD-EPI estimated",
    "glomerular filtration rate (actual body surface area, mL/min) on",
    "clearance. Interindividual variability on CL, V1 and Q,",
    "interoccasion variability on CL with each administration an occasion,",
    "and a proportional residual error."
  )
  reference <- paste(
    "Alihodzic D, Wicha SG, Frey OR, Konig C, Baehr M, Jarczak D, Kluge S,",
    "Langebrake C. Ciprofloxacin in Patients Undergoing Extracorporeal",
    "Membrane Oxygenation (ECMO): A Population Pharmacokinetic Study.",
    "Pharmaceutics. 2022;14(5):965.",
    "doi:10.3390/pharmaceutics14050965.",
    sep = " "
  )
  vignette <- "Alihodzic_2022_ciprofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "ciprofloxacin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ciprofloxacin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    BICARBONATE = list(
      description = "Plasma bicarbonate concentration",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear centred effect on CL, (1 + 0.0571 * (BICARB - 25.7))",
        "(Alihodzic 2022 Equation 1; Table 2 row 'COV_CL_BICARB (L/mmol)",
        "0.0571'). 25.7 mmol/L is the centring value printed in Equation 1;",
        "the Table 1 per-patient median is 26.6 mmol/L (IQR 24.3-29.7,",
        "range 18.0-39.0). The first covariate included (dOFV -22.5),",
        "lowering IIV on CL from 46.6% to 29.3% (Results 3.2). The model",
        "only describes the alkalotic side: the Discussion notes that the",
        "effect of acidosis on clearance was not characterised, and the",
        "linear form reaches zero clearance at 8.2 mmol/L, far below the",
        "observed range. Collected from the electronic patient record",
        "(Methods 2.3); the paper does not state whether it was treated",
        "as time-varying."
      ),
      source_name = "BICARB"
    ),
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate by the creatinine CKD-EPI",
        "equation computed with the patient's ACTUAL body surface area,",
        "in absolute mL/min (NOT indexed to 1.73 m^2)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear centred effect on CL, (1 + 0.0038 * (eGFR - 40.83))",
        "(Alihodzic 2022 Equation 1; Table 2 row 'COV_CL_eGFR (min/mL)",
        "0.0038'). Methods 2.3 computed MDRD, Cockcroft-Gault and CKD-EPI",
        "each with the actual and with a 1.73 m^2-normalised body surface",
        "area; Results 3.2 retained 'the CKD-EPI equation using the actual",
        "body surface'. 40.83 mL/min is the centring value printed in",
        "Equation 1; the Table 1 per-patient median is 46.8 mL/min (IQR",
        "25.4-79.6, range 10.7-143). Supply the de-indexed value",
        "(eGFR x BSA / 1.73), not the usual mL/min/1.73 m^2 report."
      ),
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for interoccasion variability on clearance (one occasion per administration)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Alihodzic 2022 Results 3.2: 'each administration was defined as a",
        "new occasion'. OCC = k labels the dosing interval that starts at",
        "the k-th administration (records before the first dose take",
        "OCC = 1). The paper reports one IOV magnitude and no occasion",
        "count; 21 occasions are encoded here, covering the 7-day sampling",
        "window of Methods 2.4 at the 8-hourly regimen. Occasions 2-21 are",
        "fixed equal to occasion 1 (NONMEM $OMEGA BLOCK(1) SAME). Records",
        "with OCC outside 1..21 carry no IOV."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Allometric scaling on total body weight (exponents fixed at 0.75 on CL and Q, 1 on V1 and V2; Methods 2.5) was tested but 'did not improve the model substantially' and was not retained (Results 3.2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods 2.5) but not retained; Results 3.2 states none of the other 55 tested covariates were significant."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 2.5) but not retained; only 1 of 17 patients was female (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods 2.5) but not retained; enters only through the CKD-EPI eGFR. Cohort median 1.3 mg/dL (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened as a marker of liver function (Methods 2.5) but not retained; cohort median 21.3 g/L (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Collected (Methods 2.3) and screened but not retained; cohort median 128 mg/L (Table 1)."
    ),
    RRT = list(
      description = "Renal replacement therapy indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 2.5) but not retained; eight patients received continuous veno-venous hemofiltration (Discussion)."
    ),
    ECMO_PUMP_SPEED = list(
      description = "ECMO centrifugal pump speed",
      units = "RPM (revolutions per minute)",
      type = "continuous",
      notes = "Screened with the other ECMO parameters (cannulation, blood and gas flow, duration) but not retained (Results 3.2); cohort median 3460 RPM (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 17L,
    n_studies = 1L,
    age_range = "25-73 years",
    age_median = "57 years",
    weight_range = "66-123 kg",
    weight_median = "90 kg",
    sex_female_pct = 5.9,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Critically ill adults on ECMO (8 veno-arterial, 8 veno-venous,",
      "1 veno-veno-arterial) receiving ciprofloxacin for nosocomial",
      "infection. ECMO indications: cardiogenic shock (7), ARDS (6),",
      "global respiratory insufficiency, pulmonary fibrosis, cardiac",
      "decompensation and bridge to transplant (1 each). Median SOFA 13,",
      "SAPS II 48, APACHE II 32. Median CKD-EPI eGFR 46.8 mL/min",
      "(10.7-143); median plasma bicarbonate 26.6 mmol/L (18.0-39.0)."
    ),
    dose_range = paste(
      "400 mg as a 30-min intravenous infusion every 8 h (12 patients) or",
      "every 12 h (5 patients)."
    ),
    regions = "Germany (single centre: University Medical Center Hamburg-Eppendorf).",
    notes = paste(
      "Prospective observational study (DRKS00010765), July 2016 to",
      "February 2018. Peak (1 h post-infusion), mid-dose (4 h) and trough",
      "samples on each of three days within 7 days of ECMO therapy; 139",
      "plasma samples analysed (median 8 per patient, range 3-17). Total",
      "ciprofloxacin by HPLC, LLOQ 0.5 mg/L. Fitted in NONMEM 7.4.1 with",
      "FOCE-I (Table 1, Methods 2.1-2.5)."
    )
  )

  ini({
    # Structural parameters: Alihodzic 2022 Table 2, final model. CL is the
    # typical value at plasma bicarbonate 25.7 mmol/L and eGFR 40.83 mL/min
    # (Results 3.2).
    lcl <- log(12.9)
    label("Clearance at bicarbonate 25.7 mmol/L and eGFR 40.83 mL/min (L/h)")
    # Table 2 row 'CL (L/h)' final 12.9 (7.6)
    lvc <- log(73.4)
    label("Central volume of distribution (L)")
    # Table 2 row 'V 1 (L)' final 73.4 (15)
    lq <- log(38.8)
    label("Intercompartmental clearance (L/h)")
    # Table 2 row 'Q(L/h)' final 38.8 (30)
    lvp <- log(87.6)
    label("Peripheral volume of distribution (L)")
    # Table 2 row 'V 2 (L)' final 87.6 (16)

    # Linear covariate slopes on CL (Equation 1).
    e_bicarbonate_cl <- 0.0571
    label("Fractional change in CL per mmol/L of plasma bicarbonate (L/mmol)")
    # Table 2 row 'COV_CL_BICARB (L/mmol)' 0.0571 (12)
    e_crcl_cl <- 0.0038
    label("Fractional change in CL per mL/min of CKD-EPI eGFR (min/mL)")
    # Table 2 row 'COV_CL_eGFR (min/mL)' 0.0038 (29)

    # Interindividual variability. Table 2 reports CV% = sqrt(exp(OMEGA) - 1)
    # x 100 (footnote), so omega^2 = log(CV^2 + 1). No IIV on V2 in the final
    # model (it tended to zero, Results 3.2).
    etalcl ~ 0.0644423
    # Table 2 row 'IIV CL (CV%)' final 25.8; log(0.258^2 + 1)
    etalvc ~ 0.255046
    # Table 2 row 'IIV V 1 (CV%)' final 53.9; log(0.539^2 + 1)
    etalq ~ 0.110003
    # Table 2 row 'IIV Q(CV%)' final 34.1; log(0.341^2 + 1)

    # Interoccasion variability on CL, one occasion per administration
    # (Results 3.2). One shared variance; occasions 2-21 repeat it.
    etaiov_cl_1 ~ 0.0291548
    # Table 2 row 'IOV CL (CV%)' final 17.2; log(0.172^2 + 1)
    etaiov_cl_2 ~ fixed(0.0291548)
    etaiov_cl_3 ~ fixed(0.0291548)
    etaiov_cl_4 ~ fixed(0.0291548)
    etaiov_cl_5 ~ fixed(0.0291548)
    etaiov_cl_6 ~ fixed(0.0291548)
    etaiov_cl_7 ~ fixed(0.0291548)
    etaiov_cl_8 ~ fixed(0.0291548)
    etaiov_cl_9 ~ fixed(0.0291548)
    etaiov_cl_10 ~ fixed(0.0291548)
    etaiov_cl_11 ~ fixed(0.0291548)
    etaiov_cl_12 ~ fixed(0.0291548)
    etaiov_cl_13 ~ fixed(0.0291548)
    etaiov_cl_14 ~ fixed(0.0291548)
    etaiov_cl_15 ~ fixed(0.0291548)
    etaiov_cl_16 ~ fixed(0.0291548)
    etaiov_cl_17 ~ fixed(0.0291548)
    etaiov_cl_18 ~ fixed(0.0291548)
    etaiov_cl_19 ~ fixed(0.0291548)
    etaiov_cl_20 ~ fixed(0.0291548)
    etaiov_cl_21 ~ fixed(0.0291548)

    # Residual error: proportional (Results 3.2).
    propSd <- 0.121
    label("Proportional residual error (fraction)")
    # Table 2 row 'Prop. sigma (CV%)' final 12.1 (21.9)
  })

  model({
    # Occasion indicators; a record with OCC outside 1..21 carries no IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)
    oc11 <- (OCC == 11)
    oc12 <- (OCC == 12)
    oc13 <- (OCC == 13)
    oc14 <- (OCC == 14)
    oc15 <- (OCC == 15)
    oc16 <- (OCC == 16)
    oc17 <- (OCC == 17)
    oc18 <- (OCC == 18)
    oc19 <- (OCC == 19)
    oc20 <- (OCC == 20)
    oc21 <- (OCC == 21)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10 + oc11 * etaiov_cl_11 + oc12 * etaiov_cl_12 +
      oc13 * etaiov_cl_13 + oc14 * etaiov_cl_14 + oc15 * etaiov_cl_15 +
      oc16 * etaiov_cl_16 + oc17 * etaiov_cl_17 + oc18 * etaiov_cl_18 +
      oc19 * etaiov_cl_19 + oc20 * etaiov_cl_20 + oc21 * etaiov_cl_21

    # Equation 1: CL = CLtyp * (1 + 0.0571 * (BICARB - 25.7)) *
    # (1 + 0.0038 * (eGFR - 40.83)), with lognormal IIV and IOV on CL.
    cl <- exp(lcl + etalcl + iov_cl) *
      (1 + e_bicarbonate_cl * (BICARBONATE - 25.7)) *
      (1 + e_crcl_cl * (CRCL - 40.83))
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 30-min intravenous infusions enter central (Methods 2.3).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
