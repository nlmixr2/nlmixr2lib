Jaber_2023_cyclophosphamide <- function() {
  description <- "Population PK model for phosphoramide mustard, the cytotoxic end metabolite of cyclophosphamide, after a single 2-h intravenous infusion of cyclophosphamide 50 mg/kg in adults undergoing allogeneic hematopoietic cell transplantation (Jaber 2023). Formation is delayed through two transit compartments sharing one rate constant (cyclophosphamide, then 4-hydroxycyclophosphamide) into a one-compartment phosphoramide mustard disposition with parallel linear renal clearance, proportional to weight-normalised Cockcroft-Gault creatinine clearance, and Michaelis-Menten metabolic elimination; theory-based allometry on total body weight. Only phosphoramide mustard was measured, so its clearance and volume are apparent (scaled by the fraction of the cyclophosphamide dose converted)."
  reference <- "Jaber MM, Takahashi T, Kirstein MN, Al-Kofahi M, Jacobson PA, Brundage RC. Influence of Renal Function on Phosphoramide Mustard Exposure: A Nonlinear Mixed-Effects Analysis. J Clin Pharmacol. 2023;63(1):135-141. doi:10.1002/jcph.2144"
  vignette <- "Jaber_2023_cyclophosphamide"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg: exponent 0.75 on renal clearance and Vmax, 1 on volume and 0.25 on the transit rate constant (Jaber 2023 Methods Equations 2-4; Table 2). Also enters the Cockcroft-Gault creatinine clearance and the weight normalisation of that clearance on the renal-clearance term. Cohort median 84 kg (range 47-177; Table 1). The cyclophosphamide dose used total body weight in 72% of patients and adjusted body weight (IBW + 0.5 x (TBW - IBW)) in the 28% whose weight exceeded 150% of ideal body weight; the model itself uses WT throughout.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only inside the Cockcroft-Gault creatinine clearance, (140 - AGE) * WT / (72 * CREAT), as written in the supplementary NONMEM $PK block. Cohort median 63 years (range 21-75; Table 1).",
      source_name = "AGE"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only inside the Cockcroft-Gault creatinine clearance (supplementary NONMEM $PK block, data column SRC). The 72 constant of the Cockcroft-Gault equation requires mg/dL. Eligibility required serum creatinine <= 2.0 mg/dL (Methods 'Study Setting').",
      source_name = "SRC"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Applies the Cockcroft-Gault 0.85 female factor (supplementary NONMEM $PK: 'IF(SEX==1) CRCL = 0.85*CRCL'), so the source SEX column is already coded 1 = female and SEXF = SEX with no transformation. 45% female (Table 1). Sex has no other effect in the model.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot_iv = list(analyte = "cyclophosphamide", units = "mg", specimen = "not applicable", verified = TRUE),
    transit1 = list(analyte = "4-hydroxycyclophosphamide", units = "mg", specimen = "not applicable", verified = TRUE),
    central = list(analyte = "phosphoramide mustard", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 85,
    n_studies = 1,
    n_observations = "506 plasma phosphoramide mustard concentrations, 11.4% below the 50 ng/mL limit of quantification (handled with the M3 method)",
    age_range = "21-75 years",
    age_median = "63 years",
    weight_range = "47-177 kg",
    weight_median = "84 kg",
    sex_female_pct = 45,
    renal_function = "Cockcroft-Gault creatinine clearance median 77 mL/min (range 32-150; one patient at 32, the next lowest 42); capped at 150 mL/min",
    disease_state = "Adults with high-risk hematologic malignancy (leukemia 49%, myelodysplastic syndrome 19%, lymphoma 11%, other 21%) undergoing allogeneic hematopoietic cell transplantation; antithymocyte globulin in the regimen for 33%",
    dose_range = "Single cyclophosphamide 50 mg/kg IV infused over 2 h on day -7 (dosed on total body weight in 72%, adjusted body weight in 28%), followed by fludarabine 30-40 mg/m^2/day on days -6 to -2 and 200 cGy total body irradiation on day -1",
    regions = "United States (University of Minnesota, single center)",
    notes = "Baseline characteristics from Jaber 2023 Table 1. 110 patients enrolled; 25 excluded (complications precluding PK 19, regimen changes 2, no PM samples 4). Sampling at 0, 4, 6, 8, 23, 26 and 47 h after the start of the infusion."
  )

  ini({
    # ---- Structural parameters (Jaber 2023 Table 2; supplementary NONMEM $THETA final estimates) ----
    # Typical values for a 70-kg patient; renal clearance additionally at
    # creatinine clearance 85 mL/min.
    lcl_renal <- log(14.9); label("Renal (linear) clearance of phosphoramide mustard at 70 kg and CrCL 85 mL/min (L/h)") # Table 2 ClR = 14.9 L/h (RSE 9.3%); $THETA(1) 14.913
    lvc <- log(525); label("Apparent volume of distribution of phosphoramide mustard at 70 kg (L)") # Table 2 Vc = 525 L (RSE 3.8%); $THETA(2) 524.657
    lktr <- log(1.3); label("Formation (transit) rate constant shared by both transit steps at 70 kg (1/h)") # Table 2 KTR = 1.3 1/h (RSE 6.7%); $THETA(3) 1.30505
    lvmax <- log(81.2); label("Maximum metabolic elimination rate of phosphoramide mustard at 70 kg (mg/h)") # Table 2 VMax = 81.2 mg/h (RSE 6.3%); $THETA(7) 81.2164
    lkm <- log(0.51); label("Phosphoramide mustard concentration at half-maximal metabolic rate (mg/L)") # Table 2 KM = 0.51 mg/L (RSE 15.4%); $THETA(8) 0.514354

    # ---- Covariate effects (Jaber 2023 Equations 2-4; Table 2 footnote b, held constant) ----
    e_crcl_cl_renal <- fixed(1); label("Power exponent on weight-normalised CrCL, (CrCL/85 * 70/WT), for renal clearance (unitless)") # Methods: estimated 0.97 then set to 1 (proportionality); Table 2 'ClR ~ CrCL power' = 1; $THETA(4) 1 FIX
    e_wt_cl_renal <- fixed(0.75); label("Allometric exponent of body weight on renal clearance (unitless)") # Table 2 'ClR ~ W power' = 0.75; Equation 3
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on volume (unitless)") # Table 2 'Vc ~ W power' = 1
    e_wt_ktr <- fixed(0.25); label("Allometric exponent of body weight on the transit rate constant (unitless)") # Table 2 'KTR ~ W power' = 0.25; Equation 2
    e_wt_vmax <- fixed(0.75); label("Allometric exponent of body weight on Vmax (unitless)") # Table 2 'VMax ~ W power' = 0.75; Equation 4

    # ---- Between-subject variability (supplementary NONMEM $OMEGA, variances of exponential etas) ----
    # Table 2 reports BSV% = 100 * sqrt(omega^2): 39.8, 32.9, 27.3 and 29.0%.
    # Vmax carries no BSV.
    etalcl_renal ~ 0.15861 # $OMEGA 'IIV CL0' 0.15861; Table 2 BSV ClR 39.8%
    etalvc ~ 0.108285 # $OMEGA 'IIV V' 0.108285; Table 2 BSV Vc 32.9%
    etalktr ~ 0.0744403 # $OMEGA 'IIV KTR' 0.0744403; Table 2 BSV KTR 27.3%
    etalkm ~ 0.0842706 # $OMEGA 'IIV KM' 0.0842706; Table 2 BSV KM 29.0%

    # ---- Residual error (combined, SD = sqrt(prop^2 + add^2); $SIGMA 1 held constant) ----
    propSd <- 0.097; label("Proportional residual error (fraction)") # Table 2 RUVProp = 9.7% (RSE 11.5%); $THETA(5) 0.0970147
    addSd <- 54.9; label("Additive residual error (ng/mL)") # Table 2 RUVAdd = 54.9 ng/mL (RSE 10.9%); $THETA(6) 54.8995
  })

  model({
    # 1. Cockcroft-Gault creatinine clearance (mL/min), capped at 150 mL/min,
    #    exactly as the supplementary NONMEM $PK block computes it:
    #    CRCL = (140 - AGE) * WT / (SRC * 72); IF(SEX==1) CRCL = 0.85 * CRCL;
    #    IF(CRCL >= 150) CRCL = 150.
    crcl_cg <- (140 - AGE) * WT / (72 * CREAT) * 0.85^SEXF
    crcl <- min(crcl_cg, 150)

    # 2. Individual parameters (Jaber 2023 Equations 2-5). Creatinine clearance
    #    is normalised by weight, (CrCL/85) * (70/WT), before the allometric
    #    weight term, so the net weight exponent on renal clearance is
    #    0.75 - 1 = -0.25 at the fixed CrCL exponent of 1.
    cl_renal <- exp(lcl_renal + etalcl_renal) *
      ((crcl / 85) * (70 / WT))^e_crcl_cl_renal *
      (WT / 70)^e_wt_cl_renal
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ktr <- exp(lktr + etalktr) * (WT / 70)^e_wt_ktr
    vmax <- exp(lvmax) * (WT / 70)^e_wt_vmax
    km <- exp(lkm + etalkm)

    # 3. Micro-constant of the renal (linear) pathway
    kel <- cl_renal / vc

    # 4. ODEs (Jaber 2023 Equations 7-9; NONMEM $DES). The cyclophosphamide
    #    dose is infused into depot_iv (NONMEM A(1)); transit1 is the
    #    4-hydroxycyclophosphamide pool (A(3)); central is phosphoramide
    #    mustard (A(2)). Amounts in mg throughout; km and the concentration
    #    driving Michaelis-Menten elimination are in mg/L.
    cpm_mgl <- central / vc
    d/dt(depot_iv) <- -ktr * depot_iv
    d/dt(transit1) <- ktr * depot_iv - ktr * transit1
    d/dt(central) <- ktr * transit1 - kel * central - vmax * cpm_mgl / (km + cpm_mgl)

    # 5. Observation: phosphoramide mustard plasma concentration in ng/mL
    #    (NONMEM S2 = V/1000). Below-quantification data (LLOQ 50 ng/mL) were
    #    fit with the M3 method in the source analysis.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
