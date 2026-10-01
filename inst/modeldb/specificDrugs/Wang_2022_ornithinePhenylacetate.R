Wang_2022_ornithinePhenylacetate <- function() {
  description <- "Population PK model for intravenous L-ornithine phenylacetate (L-OPA) in adults with cirrhosis or hepatic encephalopathy (Wang 2022). Phenylacetic acid (PAA): one compartment with Michaelis-Menten conversion to phenylacetylglutamine (PAGN), the only PAA elimination route. PAGN: one compartment with first-order clearance. L-ornithine (ORN): one compartment with linear clearance plus an additive endogenous baseline. Covariates: body weight and Child-Pugh class on PAA Vmax and volume; creatinine clearance on PAGN clearance; sex, weight, creatinine clearance and Child-Pugh class on ORN clearance; Child-Pugh class on the ORN baseline."
  reference <- "Wang X, Vilchez RA. Population Pharmacokinetic Analysis to Assist Dose Selection of the L-Ornithine Salt of Phenylacetic Acid. Clin Pharmacokinet. 2022;61(4):515-526. doi:10.1007/s40262-021-01075-1"
  vignette <- "Wang_2022_ornithinePhenylacetate"
  units <- list(time = "h", dosing = "mmol", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "phenylacetic acid", units = "mmol", specimen = "plasma", verified = TRUE),
    central_pagn = list(analyte = "phenylacetylglutamine", units = "mmol", specimen = "plasma", verified = TRUE),
    central_ornithine = list(
      analyte = "L-ornithine (exogenous, dose-derived)",
      units = "mmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on PAA Vmax and PAA volume normalized to 83 kg (ESM PAA control stream, (WT/83)**THETA), and on ORN clearance normalized to 75 kg (ESM ORN control stream and Table S2, (WT/75)**THETA(11)). Patient weight range 45-153 kg (Section 3.2).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Enters ORN clearance only. The source column SEX is 1 = male (ESM ORN control stream: THETA(1)*(1-SEX) + SEX*THETA(7), with THETA(1) labelled CL_female and THETA(7) CL_male; Table S2 CL female 18.2, male 25.4 L/h), so SEXF = 1 - SEX.",
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min (not stated to be BSA-normalised)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect normalized to 90 mL/min on PAGN clearance (ESM PAGN code: (CLCR/90)**THETA(3), uncapped) and on ORN clearance (ESM ORN code: capped at 90 mL/min before the power term, so renal function >= 90 mL/min carries no effect). The estimating equation is not stated in the paper. Lowest value in the patient data 26 mL/min (Section 3.1).",
      source_name = "CLCR"
    ),
    HEPIMP_MOD = list(
      description = "Child-Pugh class B indicator (1 = Child-Pugh B, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Child-Pugh A when HEPIMP_SEV is also 0)",
      notes = "Child-Pugh classification (Section 3.1). Mutually exclusive with HEPIMP_SEV. All modelled subjects had cirrhosis, so the all-zero reference is Child-Pugh A, not normal hepatic function. The source codes class as a single column (PAA code CPA0B1C2 = 0/1/2; ORN code CP). Child-Pugh B alone reduces PAA Vmax and ORN baseline; B and C pooled (HEPIMP_MOD + HEPIMP_SEV) raise PAA volume and lower ORN clearance.",
      source_name = "CPA0B1C2 / CP"
    ),
    HEPIMP_SEV = list(
      description = "Child-Pugh class C indicator (1 = Child-Pugh C, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Child-Pugh A when HEPIMP_MOD is also 0)",
      notes = "Child-Pugh classification (Section 3.1). Mutually exclusive with HEPIMP_MOD. See HEPIMP_MOD notes for the pooled B/C effects.",
      source_name = "CPA0B1C2 / CP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 152L,
    n_studies = 2L,
    age_median = "59.0 years (female), 57.0 years (male)",
    weight_range = "45-153 kg",
    weight_mean = "73.9 kg (SD 18.0) female, 86.3 kg (SD 20.1) male",
    sex_female_pct = 38.8,
    race_ethnicity = "Caucasian or unknown ethnicity (all 152 patients)",
    disease_state = "Adults with stable cirrhosis (study OCR002-HE201) or hospitalized with cirrhosis and an acute episode of hepatic encephalopathy (study OCR002-HE209).",
    dose_range = "L-OPA 1-40 g IV over 4 h or 10-40 g IV over 24 h (HE201); 10, 15 or 20 g/24 h continuous IV infusion for 5 days by Child-Pugh score (HE209).",
    hepatic_function = "Child-Pugh A 32 (21%), B 47 (31%), C 73 (48%)",
    renal_function = "Normal 78 (51%), mild impairment 41 (27%), moderate 30 (20%), severe 3 (2%)",
    notes = "Final patient models (Table 2 for PAA; ESM Table S2 for ORN; Section 3.4 and ESM PAGN code for PAGN). Demographics from Table 1. The paper also fitted healthy-subject models (46 subjects from studies OCR002-HV201 and MNK61051112, used to test Caucasian versus Chinese/Japanese ethnicity) but did not tabulate their parameter estimates, so only the patient models are encoded here."
  )

  ini({
    # ---- Phenylacetic acid (PAA): Table 2, final model in patients ----
    # Reference weight 83 kg and Child-Pugh A, from the ESM 'PAA Model in
    # Patients' control stream ((WT/83)**THETA; CPA0B1C2 = 0).
    lvmax <- log(12.4); label("Maximum PAA-to-PAGN conversion rate Vmax at 83 kg, Child-Pugh A (mmol/h)") # Table 2: Vmax 12.4 mmol/h
    lkm <- log(1.33); label("Michaelis-Menten constant Km for PAA conversion (mmol/L)") # Table 2: Km 1.33 mM (179 ug/mL)
    lvc <- log(24.4); label("PAA volume of distribution at 83 kg, Child-Pugh A (L)") # Table 2: VPAA 24.4 L
    e_wt_vmax <- 0.97; label("Power exponent of body weight on Vmax (unitless)") # Table 2: Weight power Vmax 0.97
    e_wt_vc <- 0.79; label("Power exponent of body weight on PAA volume (unitless)") # Table 2: Weight power VPAA 0.79
    e_hepimp_mod_vmax <- 0.63; label("Vmax ratio, Child-Pugh B versus A (unitless)") # Table 2: Vmax ratio C-P B/C-P A 0.63
    e_hepimp_sev_vmax <- 0.39; label("Vmax ratio, Child-Pugh C versus A (unitless)") # Table 2: C-P C/C-P A 0.39
    e_hepimp_modsev_vc <- 1.89; label("PAA volume ratio, Child-Pugh B or C versus A (unitless)") # Table 2: VPAA ratio C-P BC/C-P A 1.89

    # ---- Phenylacetylglutamine (PAGN): Section 3.4 and ESM 'NONMEM Code for PAGN' ----
    # Fitted sequentially on the individual PAA parameters. Reference
    # creatinine clearance 90 mL/min ((CLCR/90)**THETA(3)).
    lcl_pagn <- log(14.9); label("PAGN clearance at CRCL 90 mL/min (L/h)") # Section 3.4: 'apparent clearance for PAGN ... 14.9 L/h'
    lvc_pagn <- log(33.2); label("PAGN volume of distribution (L)") # Section 3.4: 'volume of distribution was 33.2 L'
    e_crcl_cl_pagn <- 0.8; label("Power exponent of creatinine clearance on PAGN clearance (unitless)") # Section 3.4: 'CL PAGN proportional to CLcr^0.8'

    # ---- L-ornithine (ORN): ESM Table S2, final model in patients ----
    # Reference: male, 75 kg, CRCL >= 90 mL/min, Child-Pugh A.
    lcl_ornithine <- log(25.4); label("Exogenous ORN clearance, male, 75 kg, CRCL >= 90 mL/min, Child-Pugh A (L/h)") # Table S2: theta7 = 25.4, 'CL, L/h (Male)'
    lvc_ornithine <- log(64.8); label("ORN volume of distribution (L)") # Table S2: 'V, L' 64.8
    lc0_ornithine <- log(13.6); label("Endogenous ORN baseline plasma concentration, Child-Pugh A (ug/mL)") # Table S2: Baseline ORN, Child-Pugh A 13.6 ug/mL
    e_sexf_cl_ornithine <- 18.2 / 25.4; label("ORN clearance ratio, female versus male (unitless)") # Table S2: theta1 = 18.2 (female) over theta7 = 25.4 (male)
    e_wt_cl_ornithine <- 0.824; label("Power exponent of body weight on ORN clearance (unitless)") # Table S2: theta11 = 0.824
    e_crcl_cl_ornithine <- 0.614; label("Power exponent of creatinine clearance (capped at 90 mL/min) on ORN clearance (unitless)") # Table S2: theta6 = 0.614
    e_hepimp_modsev_cl_ornithine <- 0.719; label("ORN clearance ratio, Child-Pugh B or C versus A (unitless)") # Table S2: theta8 = 0.719, 'Child-Pugh B/C, adjust by a coefficient'
    e_hepimp_mod_c0_ornithine <- 11.4 / 13.6; label("ORN baseline ratio, Child-Pugh B versus A (unitless)") # Table S2: Baseline ORN Child-Pugh B 11.4 over A 13.6 ug/mL
    e_hepimp_sev_c0_ornithine <- 9.56 / 13.6; label("ORN baseline ratio, Child-Pugh C versus A (unitless)") # Table S2: Baseline ORN Child-Pugh C 9.56 over A 13.6 ug/mL

    # ---- Inter-individual variability ----
    # Table 2 and Table S2 give BSV as a percentage only. Converted as
    # omega^2 = log(1 + CV^2). The control streams estimate OMEGA BLOCK(3)
    # for PAA and BLOCK(2) for ORN CL-V, but no covariance is reported, so
    # the blocks are encoded diagonally.
    etalvmax ~ 0.2074 # Table 2: BSV Vmax 48%; log(1 + 0.48^2)
    etalkm ~ 0.6133 # Table 2: BSV Km 92%; log(1 + 0.92^2)
    etalvc ~ 0.1553 # Table 2: BSV VPAA 41%; log(1 + 0.41^2)
    # PAGN: the control stream carries ETA(1) on VPAGN and ETA(2) on CLPAGN,
    # but no variance is reported anywhere in the paper or ESM.
    etalcl_pagn ~ fixed(0) # not reported; see vignette
    etalvc_pagn ~ fixed(0) # not reported; see vignette
    etalcl_ornithine ~ 0.3616 # Table S2: BSV CL 66%; log(1 + 0.66^2)
    etalvc_ornithine ~ 0.6931 # Table S2: BSV V 100%; log(1 + 1.00^2)
    etalc0_ornithine ~ 0.2476 # Table S2: BSV Baseline 53%; log(1 + 0.53^2)

    # ---- Residual error ----
    # W = SQRT(IPRED**2*THETA(4)**2 + THETA(5)**2) with SIGMA 1 FIX in both
    # control streams, so THETA(4) and THETA(5) are the proportional and
    # additive SDs.
    propSd <- 0.35; label("PAA proportional residual SD (fraction)") # Table 2: PAA proportional error 35%
    addSd <- 2.23; label("PAA additive residual SD (ug/mL)") # Table 2: PAA additive error 2.23 ug/mL
    propSd_pagn <- fixed(0); label("PAGN proportional residual SD (fraction; not reported)") # not reported; see vignette
    addSd_pagn <- fixed(0); label("PAGN additive residual SD (ug/mL; not reported)") # not reported; see vignette
    propSd_ornithine <- 0.4; label("ORN proportional residual SD (fraction)") # Table S2: error-model row 0.4 (SE 0.0029); additive THETA(5) is 0 FIX in the ORN code
  })

  model({
    # Molecular weights used by the source control streams to convert the
    # model's mmol/L concentrations to ug/mL.
    mw_paa <- 135.142 # ESM PAA code comment: 'MW of PAA = 135.142 g/mol'
    mw_pagn <- 264.281 # ESM PAA code comment: 'MW of PAGN = 264.281 g/mol'
    mw_ornithine <- 132.163 # ESM ORN code: IPRED = (A(1)/V)*132.163 + BASE

    # Child-Pugh B or C pooled (the two indicators are mutually exclusive)
    hepimp_modsev <- HEPIMP_MOD + HEPIMP_SEV

    # Creatinine clearance capped at 90 mL/min for ORN (ESM ORN code:
    # IF (CLCR.GE.90) CLCRi = 90)
    crcl_cap_ornithine <- min(CRCL, 90)

    # PAA
    vmax <- exp(lvmax + etalvmax) * (WT / 83)^e_wt_vmax *
      e_hepimp_mod_vmax^HEPIMP_MOD * e_hepimp_sev_vmax^HEPIMP_SEV
    km <- exp(lkm + etalkm)
    vc <- exp(lvc + etalvc) * (WT / 83)^e_wt_vc * e_hepimp_modsev_vc^hepimp_modsev

    # PAGN
    cl_pagn <- exp(lcl_pagn + etalcl_pagn) * (CRCL / 90)^e_crcl_cl_pagn
    vc_pagn <- exp(lvc_pagn + etalvc_pagn)

    # ORN (exogenous, dose-derived part only; the endogenous baseline is added
    # to the prediction)
    cl_ornithine <- exp(lcl_ornithine + etalcl_ornithine) *
      e_sexf_cl_ornithine^SEXF *
      (WT / 75)^e_wt_cl_ornithine *
      (crcl_cap_ornithine / 90)^e_crcl_cl_ornithine *
      e_hepimp_modsev_cl_ornithine^hepimp_modsev
    vc_ornithine <- exp(lvc_ornithine + etalvc_ornithine)
    c0_ornithine <- exp(lc0_ornithine + etalc0_ornithine) *
      e_hepimp_mod_c0_ornithine^HEPIMP_MOD * e_hepimp_sev_c0_ornithine^HEPIMP_SEV

    # Michaelis-Menten PAA -> PAGN flux (mmol/h); concentration in mmol/L
    conv_paa <- vmax * (central / vc) / (km + central / vc)

    # L-OPA dissociates 1:1 into ORN and PAA: every infusion is given as
    # equimolar doses into central (PAA) and central_ornithine (ORN).
    d/dt(central) <- -conv_paa
    d/dt(central_pagn) <- conv_paa - cl_pagn * central_pagn / vc_pagn
    d/dt(central_ornithine) <- -cl_ornithine * central_ornithine / vc_ornithine

    Cc <- central / vc * mw_paa
    Cc_pagn <- central_pagn / vc_pagn * mw_pagn
    Cc_ornithine <- central_ornithine / vc_ornithine * mw_ornithine + c0_ornithine

    Cc ~ add(addSd) + prop(propSd)
    Cc_pagn ~ add(addSd_pagn) + prop(propSd_pagn)
    Cc_ornithine ~ prop(propSd_ornithine)
  })
}
