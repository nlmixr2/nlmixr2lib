Gebreyesus_2024_cefazolin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for total plasma cefazolin in",
    "South African children under 25 kg undergoing cardiac surgery with",
    "cardiopulmonary bypass (CPB), with a third compartment for the CPB",
    "circuit. The circuit compartment exchanges with the central compartment",
    "only while the patient is on bypass; its volume is the circuit priming",
    "volume and its flow a fixed adult cardiac output (336 L/h, allometric",
    "to 70 kg), both multiplied by an estimated scaling factor (0.679).",
    "Clearances and volumes scale allometrically on body weight (reference",
    "10 kg) and clearance is proportional to bedside-Schwartz eGFR",
    "(reference 120 mL/min/1.73 m^2). A dose added to the CPB priming fluid",
    "is a bolus into the circuit compartment."
  )
  reference <- paste(
    "Gebreyesus MS, Dresner A, Wiesner L, Coetzee E, Verschuuren T,",
    "Wasmann R, Denti P. Dose optimization of cefazolin in South African",
    "children undergoing cardiac surgery with cardiopulmonary bypass.",
    "CPT Pharmacometrics Syst Pharmacol 2024;13(9):1595-1605.",
    "doi:10.1002/psp4.13196."
  )
  vignette <- "Gebreyesus_2024_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Supplementary Material S3 control stream, $MODEL: COMP = (CENTRAL),
  # (PERIPHERAL), (CPB). Concentrations were measured in plasma by LC-MS/MS
  # (Supplementary Material S1 'Drug assay'). The CPB compartment holds the
  # priming fluid plus blood in the circuit; its volume and flow are scaled
  # by a factor the authors read as the plasma proportion of whole blood
  # (Results), so it is the lumped plasma-equivalent circuit state.
  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_plasma = list(analyte = "cefazolin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Total body weight. Allometric scaling of CL and Q (exponent 0.75)",
        "and of Vc and Vp (exponent 1) to a 10 kg reference (the rounded",
        "cohort median), and of the CPB circuit flow (exponent 0.75) to a",
        "70 kg adult reference (Supplementary Material S3 $PK). Fat-free mass",
        "was tested as the size descriptor and not retained."
      ),
      source_name = "WTKG"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the bedside (modified) Schwartz equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "eGFR = 0.413 x height (cm) / serum creatinine (mg/dL)",
        "(Supplementary Material S1 Eq. 1; computed inside the control stream",
        "as CRCL_SCHW = 0.413*HTCM/CREAT). Enters clearance as the linear",
        "ratio CRCL / 120, with 120 mL/min/1.73 m^2 an adult reference rather",
        "than the cohort median (116; range 48-159, Table 2). Pre-operative",
        "value; postoperative creatinine was not available, so it is",
        "time-fixed within subject."
      ),
      source_name = "CRCL_SCHW (from HTCM, CREAT)"
    ),
    CPB_ON = list(
      description = "Cardiopulmonary bypass phase indicator (before rewarming begins)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (off the bypass circuit)",
      notes = paste(
        "Time-varying. The control stream's single SWITCH data item turns",
        "the CPB compartment on (SWITCH = 1) at the start of bypass and off",
        "(SWITCH = 2) at its end (Supplementary Material S1 'Population",
        "pharmacokinetic modelling'); the paper does not split out rewarming.",
        "The circuit is connected whenever CPB_ON + CPB_REWARM == 1, which is",
        "the switch the model uses. A user with no rewarming record sets",
        "CPB_ON = 1 for the whole time on bypass and CPB_REWARM = 0."
      ),
      source_name = "SWITCH (CPB_ON in $PK)"
    ),
    CPB_REWARM = list(
      description = "Cardiopulmonary bypass rewarming phase indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not in the rewarming phase)",
      notes = paste(
        "Time-varying. The patient remains on the bypass circuit during",
        "rewarming, so the model's circuit switch is CPB_ON + CPB_REWARM.",
        "Mutually exclusive with CPB_ON. Set to 0 throughout when the",
        "rewarming phase is not recorded separately."
      ),
      source_name = "SWITCH (CPB_ON in $PK)"
    ),
    CIRCUIT_PRIME_VOL = list(
      description = "Cardiopulmonary bypass circuit priming volume",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-subject priming fluid volume of the CPB circuit; the circuit",
        "compartment volume is CIRCUIT_PRIME_VOL x 0.679 (control stream",
        "V3 = VCPB0*THETA(8)). Median 0.40 (range 0.25-0.55) L in the study",
        "cohort (Table 2). The paper's simulations used the hospital's",
        "standard: 0.4 L for infants < 10 kg, 0.6 L for 10-20 kg and 0.8 L",
        "for 20-30 kg (Methods 'Monte Carlo simulations'). Read only while",
        "the circuit is connected."
      ),
      source_name = "VCPB0"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    n_observations = 249L,
    age_range = "1-94 months",
    age_median = "17 months (Table 2; 19.5 months in the Abstract)",
    weight_range = "2.0-22 kg",
    weight_median = "8.5 kg (Table 2; 8.7 kg in the Abstract)",
    height_range = "45-117 cm",
    sex_female_pct = 41,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Children under 25 kg undergoing cardiac surgery with cardiopulmonary",
      "bypass for congenital heart disease (e.g. atrioventricular septal",
      "defect, pulmonary atresia, tetralogy of Fallot). Expected CPB time",
      "< 30 min, beta-lactam allergy and preexisting hepatic or renal",
      "dysfunction were exclusions. Surgical duration 2.38 (0.567-5.63) h."
    ),
    renal_function = "Bedside Schwartz eGFR median 116 (range 48-159) mL/min/1.73 m^2; serum creatinine 0.294 (0.215-0.407) mg/dL (Table 2).",
    dose_range = paste(
      "50 mg/kg IV cefazolin over 45 s after induction of anaesthesia, a",
      "second dose 4-6 h later and at most three doses; 6 of 22 children also",
      "received 50 mg/kg into the CPB priming fluid (Methods 'Study design',",
      "Results)."
    ),
    regions = "South Africa (Red Cross War Memorial Children's Hospital, Cape Town)",
    notes = paste(
      "Observational study; total plasma cefazolin by LC-MS/MS, LLOQ",
      "0.586 mg/L. Fit in NONMEM 7.5.0 with FOCE-I. Baseline characteristics",
      "in Table 2."
    )
  )

  ini({
    # Final estimates: Table 3, identical to $THETA / $OMEGA of the control
    # stream in Supplementary Material S3. Typical values are for a 10 kg
    # child with eGFR 120 mL/min/1.73 m^2.
    lcl <- log(0.856); label("Clearance for 10 kg and eGFR 120 mL/min/1.73 m^2 (L/h)") # Table 3 'Clearance' 0.856 (RSE 10%); S3 THETA(1)
    lvc <- log(1.07); label("Central volume of distribution for 10 kg (L)") # Table 3 'Central volume of distribution' 1.07 (RSE 9%); S3 THETA(2)
    lq <- log(1.41); label("Intercompartmental clearance for 10 kg (L/h)") # Table 3 'Intercompartmental clearance' 1.41 (RSE 13%); S3 THETA(6)
    lvp <- log(1.61); label("Peripheral volume of distribution for 10 kg (L)") # Table 3 'Peripheral volume of distribution' 1.61 (RSE 18%); S3 THETA(5)

    # CPB circuit: flow fixed to an adult cardiac output, volume taken from
    # the per-subject priming volume, both scaled by one estimated factor.
    qbfr <- fixed(336); label("CPB circuit flow for a 70 kg adult, set to adult cardiac output (L/h)") # Table 3 'Fixed to adult cardiac output of 336 L/h'; S3 THETA(7) (336) FIX
    fcircuit <- 0.679; label("Scaling factor on the CPB circuit volume and flow (unitless)") # Table 3 'Scaling factor for the CPB compartment' 0.679 (RSE 15%); S3 THETA(8)

    # Allometric exponents fixed a priori (Supplementary Material S1).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL, Q and the CPB flow (unitless)") # S3 ALLMCL_WT and ALLMQ2 exponent 0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # S3 ALLMV_WT exponent 1

    # BSV variances from S3 $OMEGA; Table 3 reports sqrt(exp(omega^2) - 1).
    etalcl ~ 0.0523 # S3 $OMEGA 'BSV CL' 0.0523; Table 3 BSV 23.2% CV
    etalvc ~ 0.033 # S3 $OMEGA 'BSV V' 0.033; Table 3 BSV 18.3% CV
    etalvp ~ 0.0779 # S3 $OMEGA 'BSVV2' 0.0779; Table 3 BSV 28.5% CV

    # Residual error: W = sqrt(ADD^2 + PROP^2) with $SIGMA 1 FIX.
    propSd <- 0.155; label("Proportional residual error (fraction)") # Table 3 'Proportional error' 15.5% (RSE 7%); S3 THETA(3)
    addSd <- fixed(0.1172); label("Additive residual error (mg/L)") # Table 3 footnote f: 20% of LLOQ 0.586 mg/L; S3 ADD = THETA(4) + LLOQ*0.2 with THETA(4) = 0 FIX
  })
  model({
    # Circuit connected (S3: SWITCH == 1 sets CPB_ON = 1, SWITCH == 2 resets it)
    cpb_switch <- CPB_ON + CPB_REWARM

    cl <- exp(lcl + etalcl) * (WT / 10)^e_wt_cl_q * (CRCL / 120)
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc_vp
    q <- exp(lq) * (WT / 10)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 10)^e_wt_vc_vp

    # S3: Q2 = THETA(7) * (WTKG/70)**0.75 * THETA(8); V3 = VCPB0 * THETA(8)
    qcircuit <- qbfr * (WT / 70)^e_wt_cl_q * fcircuit
    vcircuit <- CIRCUIT_PRIME_VOL * fcircuit

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_central_circuit <- qcircuit / vc
    k_circuit_central <- qcircuit / vcircuit

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 +
      cpb_switch * (-k_central_circuit * central + k_circuit_central * circuit_plasma)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(circuit_plasma) <- cpb_switch * (k_central_circuit * central - k_circuit_central * circuit_plasma)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
