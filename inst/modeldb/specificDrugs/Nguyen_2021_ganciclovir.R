Nguyen_2021_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for IV ganciclovir and oral",
    "valganciclovir in a pediatric population (Nguyen 2021), with allometric",
    "body-weight scaling on all disposition parameters, a power effect of eGFR on",
    "clearance, and a multiplicative critical-illness factor on clearance.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Nguyen 2021 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Nguyen T, Oualha M, Briand C, Bendavid M, Beranger A, Benaboud S,",
    "Treluyer JM, Zheng Y, Foissac F, Winter S, et al. Population pharmacokinetics",
    "of intravenous ganciclovir and oral valganciclovir in a pediatric population",
    "to optimize dosing regimens. Antimicrob Agents Chemother. 2021;65(3):e02254-20.",
    "doi:10.1128/AAC.02254-20.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on every disposition parameter, normalized to the",
        "cohort median of 11.7 kg (Yang 2023 Table 2 reports median weight",
        "11.7 kg, range 2.6-80 kg): exponent 0.75 on CL and Q, exponent 1 on Vc",
        "and Vp. Weight was retained on CL, Vc, Q and Vp (Yang 2023 Table 4).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate (eGFR), BSA-normalized. Yang 2023",
        "Table 3 footnote defines eGFR as 'the estimated glomerular filtration",
        "rate (mL/min/1.73 m^2)'; the review does not state which estimating",
        "equation Nguyen 2021 used, so the assay form is to be confirmed against",
        "the primary publication.",
        sep = " "
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on clearance normalized to a reference of",
        "167 mL/min/1.73 m^2 with exponent 0.763: (eGFR/167)^0.763. The reference",
        "value is high for an adult scale but is consistent with a cohort whose",
        "median age is 2.5 years, where BSA-normalized eGFR routinely exceeds",
        "150 mL/min/1.73 m^2. eGFR was the only laboratory covariate tested and it",
        "was retained on CL (Yang 2023 Table 4). Stored under the canonical CRCL",
        "column, which accepts BSA-normalized creatinine-based estimates.",
        sep = " "
      ),
      source_name        = "eGFR"
    ),
    DIS_CRITILL = list(
      description        = "Critical illness / ICU-admission indicator (1 = critically ill, 0 = not critically ill)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not critically ill)",
      notes              = paste(
        "Yang 2023 Table 3 footnote: 'critically ill: 1 for critically ill",
        "patients and 0 for others'. Applied as a multiplicative power factor",
        "0.806^DIS_CRITILL on clearance, so a critically ill child has clearance",
        "19.4% lower than a non-critically-ill child of the same weight and eGFR.",
        "Critical illness was tested as a demographic covariate and retained on CL",
        "(Yang 2023 Table 4). The review's Section 3.3.3 notes that of the two",
        "studies testing transplant status, CMV presentation or critical illness,",
        "only transplant and CMV shedding reached significance -- but Table 3",
        "nonetheless prints the critical-illness factor as part of this model's",
        "final CL equation, and Table 4 lists it among the significant covariates",
        "on CL, so it is encoded here. Founding example for the DIS_CRITILL",
        "canonical.",
        sep = " "
      ),
      source_name        = "critically ill"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 105L,
    n_studies      = 1L,
    n_observations = 374L,
    age_median     = "2.5 years (range 0.01-17.3 years)",
    weight_median  = "11.7 kg (range 2.6-80 kg)",
    sex_female_pct = 43.8,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Pediatric patients receiving ganciclovir or valganciclovir; the cohort",
      "includes critically ill children, whose contrast against the",
      "non-critically-ill children in the same dataset is carried by the",
      "DIS_CRITILL covariate on clearance.",
      sep = " "
    ),
    dose_range     = paste(
      "IV ganciclovir a median 10 mg/kg/day (range 1.2-15.4) and oral",
      "valganciclovir a median 36 mg/kg/day (range 14.6-83.8), both administered",
      "twice daily.",
      sep = " "
    ),
    regions        = "France (prospective).",
    bioassay       = "LC-MS, LLOQ 0.05 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2 (105 subjects, 59 male /",
      "46 female, 374 observations); the sampling schedule is recorded as not",
      "reported. Estimated in Monolix with the SAEM algorithm. Covariates tested:",
      "weight, age, height, BSA, sex, critical illness and eGFR, with forward",
      "inclusion at p < 0.05 and backward elimination at p < 0.01; weight, eGFR",
      "and critical illness were retained on CL and weight on Vc, Q and Vp",
      "(Yang 2023 Table 4). Internal validation by GOF, NPDE and pcVPC; externally",
      "validated in 35 subjects. This study did NOT convert valganciclovir doses",
      "to ganciclovir equivalents, so oral doses are supplied as valganciclovir",
      "milligrams and the estimated bioavailability F = 0.438 absorbs the",
      "molecular-weight ratio. The primary study's stated purpose was to design",
      "dosing regimens targeting a preventive AUC0-24h of 40-80 mg*h/L and a",
      "curative AUC0-24h of 80-120 mg*h/L.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Nguyen et al. (2021) row. Reference
    # subject: WT = 11.7 kg (the cohort median), eGFR = 167 mL/min/1.73 m^2,
    # not critically ill. Clearances in L/h, volumes in L, ka in 1/h. All values
    # are typical (population) values of the final model; the review reports no
    # standard errors. No absorption lag time is reported for this model.
    lcl <- log(2.55) ; label("Clearance at WT = 11.7 kg, eGFR = 167 mL/min/1.73 m^2, not critically ill (CL, L/h)") # Yang 2023 Table 3 (Nguyen 2021): CL = 2.55 * (BW/11.7)^0.75 * (eGFR/167)^0.763 * 0.806^critically ill
    lvc <- log(5.96) ; label("Central volume of distribution at WT = 11.7 kg (Vc, L)")                              # Yang 2023 Table 3 (Nguyen 2021): Vc = 5.96 * (BW/11.7)
    lq  <- log(0.222); label("Inter-compartmental clearance at WT = 11.7 kg (Q, L/h)")                              # Yang 2023 Table 3 (Nguyen 2021): Q = 0.222 * (BW/11.7)^0.75
    lvp <- log(1.29) ; label("Peripheral volume of distribution at WT = 11.7 kg (Vp, L)")                           # Yang 2023 Table 3 (Nguyen 2021): Vp = 1.29 * (BW/11.7)
    lka <- log(0.506); label("First-order oral absorption rate constant (ka, 1/h)")                                 # Yang 2023 Table 3 (Nguyen 2021): Ka = 0.506

    # Oral bioavailability of ganciclovir from a valganciclovir milligram dose
    # (no molecular-weight conversion was applied by this study).
    lfdepot <- log(0.438); label("Oral bioavailability of ganciclovir from valganciclovir (F, fraction)")            # Yang 2023 Table 3 (Nguyen 2021): F = 0.438

    # Covariate effects. The eGFR exponent 0.763 is a non-canonical estimated
    # value; the body-weight exponents 0.75 (CL, Q) and 1 (Vc, Vp) are the
    # canonical allometric values imposed by the authors and are encoded as fixed.
    e_crcl_cl        <- 0.763    ; label("Power exponent of eGFR on CL (unitless; reference 167 mL/min/1.73 m^2)")   # Yang 2023 Table 3 (Nguyen 2021): (eGFR/167)^0.763
    e_wt_cl_q        <- fixed(0.75); label("Allometric exponent of WT on CL and Q (unitless; reference 11.7 kg)")    # Yang 2023 Table 3 (Nguyen 2021): (BW/11.7)^0.75 on CL and Q
    e_wt_vc_vp       <- fixed(1)   ; label("Allometric exponent of WT on Vc and Vp (unitless; linear, reference 11.7 kg)") # Yang 2023 Table 3 (Nguyen 2021): (BW/11.7) on Vc and Vp
    e_dis_critill_cl <- 0.806    ; label("Multiplicative critical-illness factor on CL (unitless; non-critically-ill reference)") # Yang 2023 Table 3 (Nguyen 2021): 0.806^critically ill

    # Between-subject variability. Yang 2023 Methods define the reported
    # percentage as %CV = sqrt(omega^2) * 100%, i.e. omega equals the tabulated
    # percentage divided by 100, so the variance is (BSV% / 100)^2. Q, Vp, ka and
    # F carry no BSV in the source table.
    etalcl ~ 0.236196  # Yang 2023 Table 3 (Nguyen 2021): BSV CL = 48.6% -> 0.486^2
    etalvc ~ 0.219961  # Yang 2023 Table 3 (Nguyen 2021): BSV Vc = 46.9% -> 0.469^2

    # Residual unexplained variability: proportional.
    propSd <- 0.477; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Nguyen 2021): 47.7% proportional error
  })

  model({
    # The critical-illness effect is a multiplicative power factor on the base
    # (indicator-free) clearance, so DIS_CRITILL = 0 recovers the reference value.
    cl <- exp(lcl + etalcl) * (WT / 11.7)^e_wt_cl_q * (CRCL / 167)^e_crcl_cl *
      e_dis_critill_cl^DIS_CRITILL
    vc <- exp(lvc + etalvc) * (WT / 11.7)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 11.7)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 11.7)^e_wt_vc_vp
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral valganciclovir doses enter `depot` as valganciclovir milligrams;
    # IV ganciclovir doses go directly to `central`.
    f(depot) <- fdepot

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
