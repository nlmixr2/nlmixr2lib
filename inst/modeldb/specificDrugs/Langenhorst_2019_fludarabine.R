Langenhorst_2019_fludarabine <- function() {
  description <- "Three-compartment IV population PK model for the circulating fludarabine metabolite F-ara-A in children and adults (n = 258, age 0.3-74 years) receiving fludarabine phosphate as part of myeloablative conditioning prior to allogeneic hematopoietic cell transplantation (Langenhorst 2019). Total clearance is decomposed into a non-renal component (3.2 L/h at 70 kg) and an eGFR-driven renal component with slope 0.78 (unitless), both allometrically scaled to actual body weight at a fixed 0.75 exponent referenced to 70 kg: CL = (3.2 + eGFR_L_per_h * 0.78) * (BW/70)^0.75 L/h, where eGFR is expressed in L/h (converted from Cockcroft-Gault / Schwartz eGFR in mL/min/1.73 m^2 by * 60/1000). Central volume V1 (39 L), first peripheral V2 (20 L), and second peripheral V3 (50 L) all scale allometrically at a fixed 1.0 exponent to 70 kg. Inter-compartmental clearances Q2 (8.6 L/h, V1-V2) and Q3 (3.8 L/h, V1-V3) scale at a fixed 0.75 exponent to 70 kg. Because the random effects on the volume triplet (V1, V2, V3) and on the clearance triplet (CL, Q2, Q3) were highly correlated in the original NONMEM fit, the paper encodes a SINGLE shared eta on {V1, V2, V3} (48% CV) and a SINGLE shared eta on {CL, Q2, Q3} (23% CV) -- implemented here as identical eta terms across the three parameters in each group (perfect correlation, matching the paper's coding). Inter-occasion variability (12% CV on CL/Q2/Q3, 31% CV on V1/V2/V3, one occasion per dose) reported by Langenhorst 2019 Table 2 is NOT encoded structurally: the source paper defines each dose as an occasion but nlmixr2lib omits IOV when no operational occasion column is defined for downstream use (Brooks 2021 / Andrews 2017 precedent); downstream users can add an OCC indicator and per-occasion etas in rxode2. Residual error is proportional (6.3% CV) on the linear concentration scale."
  reference <- "Langenhorst JB, Dorlo TPC, van Maarseveen EM, Nierkens S, Kuball J, Boelens JJ, van Kesteren C, Huitema ADR. Population Pharmacokinetics of Fludarabine in Children and Adults during Conditioning Prior to Allogeneic Hematopoietic Cell Transplantation. Clin Pharmacokinet. 2019;58(5):627-637. doi:10.1007/s40262-018-0715-9"
  vignette <- "Langenhorst_2019_fludarabine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed within the pharmacokinetic sampling window (Langenhorst 2019 Table 1: median 60 kg across 258 patients, range 4.3-130 kg; pediatric subgroup 4.3-96 kg, adult subgroup 47-130 kg). Used as the size descriptor for allometric scaling of all six structural parameters (CL, Q2, Q3, V1, V2, V3) with fixed exponents 0.75 for clearances and 1.0 for volumes, referenced to 70 kg (Langenhorst 2019 Section 3.2 'Structural and Stochastic Model' and Table 2). Alternative body-size descriptors (fat-free mass, body surface area) were tested but did not improve model fit (Section 3.3, dOFV +60 and +68 respectively). Actual body weight (not ideal / adjusted / fat-free) is the descriptor to use.",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalized (mL/min/1.73 m^2), capped at a maturation-adjusted maximum",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts BSA-normalized creatinine-based renal-function estimates, following the eGFR precedent of Cirincione 2017 MDRD eGFR, Bajaj 2017 CKD-EPI eGFR, and Andrews 2017 adapted-Schwartz eGFR). Source paper's column is eGFR. Calculated from the mean serum creatinine value available between day -7 and day 0 prior to fludarabine infusion (Langenhorst 2019 Section 2.5). The Cockroft-Gault formula was used for adults (women >= 17 years, men >= 14 years; result BSA-normalized as 1.73 * CrCl / BSA to give mL/min/1.73 m^2) and the Schwartz formula for pediatric patients (women < 17 years, men < 14 years). Values were capped at a maturation-adjusted maximum: eGFR_max = 140 mL/min/1.73 m^2, ramping UP from 35 (25% of 140) at birth to 140 at 1.5 years of age, then held constant to age 30 years, then declining by 8 mL/min/1.73 m^2 per decade thereafter. Population range 25-140 mL/min/1.73 m^2, median 120 (pediatric median 140; adult median 110). In model() the CRCL value is converted from mL/min/1.73 m^2 to L/h via * 60/1000 (= * 0.06) so that the paper's unit-less Slope_pop = 0.78 (parameter e_crcl_cl) multiplies directly to give the renal-clearance contribution at the 70-kg reference. Used only in the CL equation; V1, V2, V3, Q2, Q3 depend only on WT.",
      source_name        = "eGFR"
    )
  )

  population <- list(
    species                 = "human",
    n_subjects              = 258L,
    n_studies               = 1L,
    n_pk_samples            = 2605L,
    n_doses                 = 596L,
    age_range               = "0.3-74 years",
    age_median              = "18 years",
    weight_range            = "4.3-130 kg",
    weight_median           = "60 kg",
    sex_female_pct          = 38.0,
    race_ethnicity          = "Not reported in source",
    disease_state           = "Allogeneic hematopoietic cell transplantation recipients across a full range of HCT indications (benign disorders 27%, acute leukemia 45%, lymphoma 6.6%, myelodysplastic syndrome 12%, plasma cell disorder 8.9%). Pediatric (age <= 20 years) 52%, adult (age > 20 years) 48%.",
    renal_function          = "eGFR range 25-140 mL/min/1.73 m^2 (median 120); Cockroft-Gault for adults, Schwartz for pediatric, both BSA-normalized to mL/min/1.73 m^2 and capped at a maturation-adjusted maximum described in Section 2.5",
    dose_range              = "160 mg/m^2 cumulative fludarabine phosphate (n = 197) OR 40 mg/m^2 cumulative fludarabine phosphate combined with 120 mg/m^2 clofarabine (n = 61); administered as a 1-h IV infusion of fludarabine phosphate (F-ara-AMP prodrug, rapidly converted to circulating F-ara-A) once daily for 4 days (day -5 to day -2 pre-transplant)",
    concomitant_medications = "Busulfan (myeloablative, targeted to 90 mg*h/L cumulative AUC_T0-inf; 30 mg*h/L for Fanconi anemia) as a 3-h IV infusion directly following each fludarabine infusion. Rabbit ATG in the unrelated-donor HCT setting (10 mg/kg for children < 30 kg, 7.5 mg/kg for children > 30 kg over 4 daily doses from day -9 for children; 6 mg/kg over four 12-h infusions from day -12 for adults). Clofarabine 30 mg/m^2/day x 4 doses preceding fludarabine in the low-dose pediatric-malignancy subgroup. Clemastine, paracetamol, and 2 mg/kg prednisolone (max 100 mg) were administered IV prior to ATG. Neither clofarabine co-administration nor ATG was retained as a covariate in the final PK model (Langenhorst 2019 Section 3.3).",
    regions                 = "Single centre: University Medical Centre Utrecht, The Netherlands",
    sampling_window         = "Samples drawn on days 1 (42%), 2 (17%), 3 (4%), and 4 (37%) of conditioning. Routine TDM sampling at 4, 5, 6, and 7 h after end of fludarabine infusion; additional samples 7-24 h post-infusion for a subset; from January 2016 onwards, additional samples 15-45 min after end of fludarabine infusion (pre-busulfan). Median 10 samples per patient (range 3-19); 116/596 doses (19%) had peak samples (< 3 h post infusion end), 117/596 (20%) had trough samples (> 8 h).",
    assay                   = "F-ara-A (the circulating metabolite of fludarabine) quantified in plasma by validated liquid chromatography mass spectrometry (LC-MS), LLOQ 0.001 mg/L; per FDA / EMA bioanalytical validation guidelines. None of the 2605 samples were below the LLOQ.",
    iov_structure           = "Each dose + subsequent sampling was defined as a separate occasion in the paper. IOV on CL/Q2/Q3 (12% CV) and V1/V2/V3 (31% CV) was reported (Langenhorst 2019 Table 2). NOT encoded structurally in this file - the source paper does not define an operational occasion column suitable for downstream model-library use (Brooks 2021 / Andrews 2017 precedent); downstream users who want IOV can add an OCC indicator column and per-occasion etas in rxode2.",
    notes                   = "Retrospective analysis of PK samples acquired May 2010 - January 2017 during routine busulfan TDM (protocol UMCU 11/063). Non-linear mixed-effects modelling in NONMEM 7.3.0 with FOCE-I. Pirana 2.9.5 + R 3.3.3 for workflow / visualization. Non-parametric bootstrap n = 1000 (95% success); prediction-corrected VPCs and NPDE from 1000 simulations."
  )

  ini({
    # Structural PK parameters - all population estimates from
    # Langenhorst 2019 Table 2 (final model for a 70-kg reference
    # subject). Clearance is decomposed into a non-renal component
    # and an eGFR-driven renal component; the paper's Eq. 7 /
    # Table 2 formula is:
    #   CL = (CL_nonrenal_70 + eGFR_L_per_h * Slope_pop) * (BW/70)^0.75
    # where eGFR_L_per_h = CRCL * 60/1000 is the eGFR (mL/min/1.73 m^2)
    # converted to L/h so that Slope_pop is unitless. Reference
    # weight is 70 kg; reference eGFR for the equation is 0 (the
    # non-renal component is what remains when eGFR = 0).
    lcl_nonren <- log(3.2) ; label("Non-renal clearance at 70 kg reference (L/h)")               # Langenhorst 2019 Table 2 Cl_70kg,non-renal = 3.2 L/h (95% CI 1.6-4.9, RSE 20%; bootstrap median 3.2)
    lvc        <- log(39)  ; label("Central volume V1 at 70 kg reference (L)")                    # Langenhorst 2019 Table 2 V1,70kg = 39 L (95% CI 33-45; bootstrap median 39)
    lvp        <- log(20)  ; label("First peripheral volume V2 at 70 kg reference (L)")           # Langenhorst 2019 Table 2 V2,70kg = 20 L (95% CI 17-24; bootstrap median 21)
    lvp2       <- log(50)  ; label("Second peripheral volume V3 at 70 kg reference (L)")          # Langenhorst 2019 Table 2 V3,70kg = 50 L (95% CI 41-58; bootstrap median 51)
    lq         <- log(8.6) ; label("Inter-compartmental clearance V1<->V2 (Q2) at 70 kg reference (L/h)") # Langenhorst 2019 Table 2 Q2,70kg = 8.6 L/h (95% CI 6.8-10; bootstrap median 8.8)
    lq2        <- log(3.8) ; label("Inter-compartmental clearance V1<->V3 (Q3) at 70 kg reference (L/h)") # Langenhorst 2019 Table 2 Q3,70kg = 3.8 L/h (95% CI 2.9-4.6; bootstrap median 3.7)

    # Covariate effect: eGFR slope on CL. Represents the unit-less
    # "fraction of eGFR accounting for renal clearance of fludarabine"
    # (Langenhorst 2019 Section 3.3, paragraph after Eq. 7). eGFR is
    # converted from CRCL (mL/min/1.73 m^2) to L/h inside model() by
    # * 60/1000 before multiplication so the coefficient stays as
    # reported.
    e_crcl_cl <- 0.78 ; label("CRCL slope on CL (unitless; multiplied by CRCL converted to L/h)") # Langenhorst 2019 Table 2 Slope_pop = 0.78 (95% CI 0.57-1.0, RSE 11%; bootstrap median 0.79)

    # Inter-individual variability. Langenhorst 2019 Section 3.2:
    # "the random effects on volume (V1, V2, V3) and clearances (CL,
    # Q2, Q3) were highly correlated. Therefore, single random
    # effects (both IIV and IOV) were estimated for V1, V2, and V3,
    # and for CL, Q2, Q3, respectively." Encoded here as identical
    # eta terms shared across the three parameters in each triplet
    # (perfect correlation, matching the paper's coding). Log-normal
    # IIV: omega^2 = log(CV^2 + 1).
    #   CL/Q2/Q3 group: 23% CV -> omega^2 = log(0.23^2 + 1) = 0.05155
    #   V1/V2/V3 group: 48% CV -> omega^2 = log(0.48^2 + 1) = 0.20735
    etalcl ~ 0.05155  # Langenhorst 2019 Table 2 IIV CL/Q2/Q3 = 23% CV (95% CI 15-31%, shrinkage 7%); shared across CL, Q2, Q3
    etalvc ~ 0.20735  # Langenhorst 2019 Table 2 IIV V1/V2/V3 = 48% CV (95% CI 36-60%, shrinkage 15%); shared across V1, V2, V3

    # Residual error. Proportional model on the linear concentration
    # scale (mg/L). Langenhorst 2019 Table 2: 6.3% CV.
    propSd <- 0.063 ; label("Proportional residual error (fraction)")  # Langenhorst 2019 Table 2 proportional residual error = 6.3% (95% CI 4.3-8.3%, shrinkage 20%)
  })

  model({
    # Convert CRCL from mL/min/1.73 m^2 to L/h so that the unit-less
    # Slope_pop = e_crcl_cl (0.78) multiplies directly to give the
    # renal-clearance contribution at the 70-kg reference:
    #   * 60 min/h * 0.001 L/mL = * 0.06.
    # Langenhorst 2019 Eq. 7 / Table 2 formula for CL.
    egfr_lh <- CRCL * 60 / 1000

    # Individual PK parameters. A single shared eta drives all three
    # clearance parameters and a separate single shared eta drives
    # all three volume parameters (Langenhorst 2019 Section 3.2 --
    # perfect correlation across each triplet). Allometric exponents
    # 0.75 (clearances) and 1.0 (volumes) are hardcoded per the
    # paper's explicit choice: the exponents were TESTED for
    # estimation and the estimated values were "very close to 0.75
    # and 1.0, respectively, and did not result in a relevant
    # improvement in model fit" so the fixed exponents were kept
    # (Section 3.3, paragraph after Eq. 7).
    cl  <- (exp(lcl_nonren) + egfr_lh * e_crcl_cl) * (WT / 70)^0.75 * exp(etalcl)
    q   <- exp(lq)  * (WT / 70)^0.75 * exp(etalcl)
    q2  <- exp(lq2) * (WT / 70)^0.75 * exp(etalcl)

    vc  <- exp(lvc)  * (WT / 70) * exp(etalvc)
    vp  <- exp(lvp)  * (WT / 70) * exp(etalvc)
    vp2 <- exp(lvp2) * (WT / 70) * exp(etalvc)

    # Micro-constants for the three-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment IV linear PK model. Fludarabine phosphate is
    # administered as a 1-h IV infusion into the central compartment
    # via the rate/dur column in the event table; the model does not
    # carry a separate depot compartment. F-ara-AMP is very rapidly
    # and fully converted to the circulating F-ara-A metabolite that
    # is quantified in plasma (Langenhorst 2019 Introduction /
    # Section 2.3), so the model treats the infused dose as directly
    # entering the central F-ara-A compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Observation. Plasma concentration of the circulating F-ara-A
    # metabolite in central. Dose in mg divided by vc in L gives
    # mg/L (= ug/mL), the units used throughout the paper (LLOQ
    # 0.001 mg/L; AUC reported in mg*h/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
