Duong_2017_AT9283 <- function() {
  description <- "Two-compartment IV population PK model for AT9283 (aurora kinase inhibitor) in adults and children with leukaemia or solid tumours (Duong 2017): allometric body-weight scaling on all four disposition parameters (CL, Vc, Q, Vp) with a power effect of estimated GFR on CL. Population residual error switches between adults (combined additive + proportional) and children (additive only) via the CHILD binary indicator."
  reference <- paste(
    "Duong JK, Griffin MJ, Hargrave D, Vormoor J, Edwards D, Boddy AV.",
    "A population pharmacokinetic model of AT9283 in adults and children",
    "to predict the maximum tolerated dose in children with leukaemia.",
    "Br J Clin Pharmacol. 2017;83(8):1713-1722.",
    "doi:10.1111/bcp.13260.",
    sep = " "
  )
  vignette <- "Duong_2017_AT9283"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used for allometric scaling of all four PK parameters (CL, Vc, Q,",
        "Vp) standardised to a reference body weight of 70 kg (Duong 2017",
        "Methods Equations 2 and 3). Allometric exponents fixed a priori",
        "at 0.75 on clearances (CL, Q) and 1.0 on volumes (Vc, Vp), citing",
        "Anderson and Holford 2008. Cohort range from Table 2: 8.9-120.5",
        "kg across the pooled adult and paediatric trials (medians:",
        "adults solid tumour 73.6 kg, adults leukaemia 67.4 kg, children",
        "solid tumour 29.2 kg, children leukaemia 16.1 kg)."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalised mL/min/1.73 m^2; MDRD formula for adults, bedside Schwartz formula for children under 18 years)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Estimated GFR by the Modification of Diet in Renal Disease",
        "(MDRD) formula for adults (Duong 2017 Equation 4) and the bedside",
        "Schwartz formula for children aged under 18 years (Equation 5).",
        "Both formulas produce BSA-normalised values in mL/min/1.73 m^2.",
        "Power effect on clearance with reference 100 mL/min, equivalent",
        "to 6 L/h in the source paper's Methods text (Equation 6); the",
        "GFR exponent is 0.453 (Table 3). Stored under canonical CRCL",
        "because the register lists MDRD eGFR and bedside-Schwartz GFR as",
        "accepted source assays. Cohort GFR range from Table 2: 31.9-299.4",
        "mL/min/1.73 m^2 (adults median 77.1, children median 132.9)."
      ),
      source_name        = "GFR"
    ),
    CHILD = list(
      description        = "Child age-cohort indicator (1 = paediatric subject under 18 years of age, 0 = adult)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = paste(
        "Selects the residual error structure: adults (CHILD = 0) use a",
        "combined additive + proportional residual error, children",
        "(CHILD = 1) use an additive-only residual error (Duong 2017",
        "Results section and Table 3). The split reflects assay and",
        "sample-collection differences between the adult Astex",
        "Pharmaceuticals trials and the paediatric Cancer Research UK",
        "trials. The age cutoff (under 18 years) also drives the choice",
        "of renal-function estimating equation (MDRD for adults, bedside",
        "Schwartz for children); CRCL is computed at data-preparation",
        "time outside the model."
      ),
      source_name        = "Study group (adult vs paediatric trial)"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 92L,
    n_adults            = 53L,
    n_children          = 39L,
    n_studies           = 4L,
    n_concentrations    = 1770L,
    age_range           = "1-86 years (pooled across four Phase I trials; Table 2)",
    age_median          = "Adults solid tumour: 63 years; adults leukaemia: 54 years; children solid tumour: 9 years; children leukaemia: 3 years (Table 2)",
    weight_range        = "8.9-120.5 kg (Table 2)",
    weight_median       = "Adults solid tumour: 73.6 kg; adults leukaemia: 67.4 kg; children solid tumour: 29.2 kg; children leukaemia: 16.1 kg (Table 2)",
    bsa_range           = "0.44-2.50 m^2 (Table 2)",
    bmi_range           = "13.2-43.6 kg/m^2 (Table 2)",
    crcl_range          = "31.9-299.4 mL/min/1.73 m^2 (Table 2); adults median 77.1 mL/min/1.73 m^2 (most with mild-to-moderately reduced kidney function), children median 132.9 mL/min/1.73 m^2 (predominantly normal or elevated kidney function)",
    sex_female_pct      = "Not reported per group in Duong 2017",
    race_ethnicity      = "Not reported",
    disease_state       = "Patients with relapsed or refractory solid tumours (n = 61) or relapsed or refractory leukaemias (n = 31) across four Phase I dose-escalation trials. Adults had a mix of normal and mild-to-moderately reduced kidney function; children had predominantly normal or elevated kidney function.",
    dose_range          = "4.5-486 mg/m^2 per 72 h as a continuous intravenous infusion every 21 days. Doses adjusted to body surface area by the Mosteller formula. Identified MTDs (mg/m^2 per 72 h): 27 in adults with solid tumours, 324 in adults with leukaemia, 55.5 in children with solid tumours; the paediatric leukaemia trial was terminated before the MTD was reached, with simulated MTD estimated at 30 mg/kg per 72 h.",
    regions             = "USA and UK (Astex Pharmaceuticals sponsored the adult studies; Cancer Research UK sponsored the paediatric studies)",
    trials              = "Arkenau et al. (NCT00443976, adults solid tumour); Foran et al. (NCT00522990, adults leukaemia); Moreno et al. (NCT00985868, children solid tumour); Cancer Research UK (NCT01431664, children leukaemia)",
    sampling            = "Adult trials: 0.5, 1, 4, 8, 12, 22, 32, 46, 56, 70, 72, 72.05, 72.25, 72.30, 72.45, 73, 74, 75, 76, 78, 80, 84, 96 h and day 8 in cycles 1 and 2. Paediatric trials: 0, 4, 24, 48, 70, 73, 76 and 96 h after the start of the cycle-1 infusion (Table 1).",
    assay               = "Validated LC-MS/MS assay; calibration range 0.1-500 ng/mL; lower limit of quantification 0.1 ng/mL. Approximately 3 percent of samples were below the LLOQ and were excluded from the analysis.",
    notes               = "Demographics from Duong 2017 Tables 1 and 2. 1770 plasma AT9283 concentrations from 92 patients used in the final population PK analysis. NONMEM v7.3 with first-order conditional estimation with interaction (FOCE-I); observations log-transformed prior to fitting. Final model selection driven by drop in objective function value with backward elimination (DELTA OFV > 6.63, P < 0.01)."
  )

  ini({
    # ===== Structural PK parameters (Duong 2017 Table 3, final-model column) =====
    # Reference subject: WT = 70 kg, CRCL = 100 mL/min/1.73 m^2 (= 6 L/h).
    lcl <- log(32.3); label("Typical clearance CL at WT=70 kg, CRCL=100 mL/min/1.73 m^2 (L/h)")  # Duong 2017 Table 3 final model: CL = 32.3 L/h/70 kg (RSE 5%; bootstrap median 32.2, 95% CI 30.0-34.9)
    lvc <- log(58.6); label("Typical central volume Vc at WT=70 kg (L)")                          # Duong 2017 Table 3 final model: Vc = 58.6 L/70 kg (RSE 7%; bootstrap median 58.5, 95% CI 50.1-64.9)
    lq  <- log(38.5); label("Typical intercompartmental clearance Q at WT=70 kg (L/h)")           # Duong 2017 Table 3 final model: Q  = 38.5 L/h/70 kg (RSE 12%; bootstrap median 39.2, 95% CI 32.5-49.2)
    lvp <- log(162);  label("Typical peripheral volume Vp at WT=70 kg (L)")                       # Duong 2017 Table 3 final model: Vp = 162 L/70 kg (RSE 6%; bootstrap median 162, 95% CI 148.3-179.5)

    # ===== Allometric exponents (Duong 2017 Methods Equations 2 and 3) =====
    # Standard Anderson-Holford 2008 values held fixed: 0.75 on CL and Q,
    # 1.0 on Vc and Vp, all standardised to 70 kg. The paper does not report
    # an RSE for these exponents; the prose "An allometric weight model was
    # used to standardize all pharmacokinetic parameters to a body weight of
    # 70 kg" combined with Equations 2-3 indicates the canonical 0.75 / 1.0
    # were held fixed during estimation, matching the convention in
    # AbdulAziz_2016_doripenem.R and other allometric-fixed registry models.
    e_wt_cl_q   <- fixed(0.75); label("Allometric exponent on (WT/70) shared by CL and Q (unitless, fixed)")  # Duong 2017 Methods, Equation 2 (clearances)
    e_wt_vc_vp  <- fixed(1.0);  label("Allometric exponent on (WT/70) shared by Vc and Vp (unitless, fixed)") # Duong 2017 Methods, Equation 3 (volumes)

    # ===== Renal-function effect on CL (Duong 2017 Equation 6, Table 3) =====
    # Power effect of estimated GFR on CL, normalised to 100 mL/min (= 6 L/h):
    #   CL = CL_pop * (WT/70)^0.75 * (CRCL/100)^theta_EC
    e_crcl_cl   <- 0.453; label("Power exponent on (CRCL / 100) for CL (unitless)") # Duong 2017 Table 3 final model: GFR exponent = 0.453 (RSE 23%; bootstrap median 0.452, 95% CI 0.206-0.606)
    crcl_ref_cl <- 100;   label("Reference CRCL for the GFR power effect on CL (mL/min/1.73 m^2)") # Duong 2017 Methods: "normalized to a standard of 6 l h-1 (100 ml min-1)"

    # ===== Inter-individual variability (Duong 2017 Table 3 final-model estimates) =====
    # CV%; omega^2 = log(1 + CV^2). Paper states the IIV was estimated with a
    # full covariance matrix (correlations between all four etas), but only
    # the diagonal CVs are reported in Table 3. Encoded here as diagonal-only
    # IIV; the missing off-diagonals are documented in the vignette's
    # Assumptions and deviations section.
    etalcl ~ 0.16882   # Duong 2017 Table 3 final model: IIV CL = 42.9% (RSE 9%; bootstrap median 43.1, 95% CI 36.2-49.6); log(1 + 0.429^2) = 0.16882
    etalvc ~ 0.08510   # Duong 2017 Table 3 final model: IIV Vc = 29.8% (RSE 17%; bootstrap median 30.5, 95% CI 22.0-40.8); log(1 + 0.298^2) = 0.08510
    etalq  ~ 0.46038   # Duong 2017 Table 3 final model: IIV Q  = 77.0% (RSE 18%; bootstrap median 74.1, 95% CI 44.5-98.7); log(1 + 0.770^2) = 0.46038
    etalvp ~ 0.14148   # Duong 2017 Table 3 final model: IIV Vp = 38.9% (RSE 13%; bootstrap median 38.6, 95% CI 30.4-47.7); log(1 + 0.389^2) = 0.14148

    # ===== Residual error (Duong 2017 Table 3 final-model estimates) =====
    # Population-specific residual error structure: adults use combined
    # additive + proportional, children use additive only. Switched in
    # model() by the CHILD binary covariate.
    addSdAdult  <- 0.166; label("Adult additive residual error (ng/mL)")            # Duong 2017 Table 3 final model: Adults additive = 0.166 ng/mL (RSE 6%; bootstrap median 0.163, 95% CI 0.145-0.181)
    propSdAdult <- 0.499; label("Adult proportional residual error (fraction)")     # Duong 2017 Table 3 final model: Adults proportional = 49.9% (RSE 24%; bootstrap median 50.0, 95% CI 28.1-75.8)
    addSdChild  <- 0.359; label("Paediatric additive residual error (ng/mL)")       # Duong 2017 Table 3 final model: Children additive = 0.359 ng/mL (RSE 8%; bootstrap median 0.359, 95% CI 0.308-0.409)
  })

  model({
    # ===== Individual PK parameters =====
    # Allometric weight scaling on all four parameters (Equations 2 and 3)
    # plus a power effect of estimated GFR on CL (Equation 6).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * (CRCL / crcl_ref_cl)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # ===== Micro-constants =====
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ===== ODE system =====
    # AT9283 is administered as a continuous 72-hour intravenous infusion
    # every 21 days (Duong 2017 Methods); dose is delivered directly into
    # the central compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ===== Observation =====
    # central is in mg, vc is in L, so central/vc is in mg/L = ug/mL.
    # Multiply by 1000 to report Cc in ng/mL (matching the calibration range
    # 0.1-500 ng/mL in Duong 2017 Methods).
    Cc <- (central / vc) * 1000

    # ===== Population-specific residual error =====
    # Adults (CHILD = 0): combined additive + proportional residual error.
    # Children (CHILD = 1): additive-only residual error; proportional term
    # collapses to zero.
    addSd_pop  <- addSdAdult  * (1 - CHILD) + addSdChild * CHILD
    propSd_pop <- propSdAdult * (1 - CHILD)
    Cc ~ add(addSd_pop) + prop(propSd_pop)
  })
}
