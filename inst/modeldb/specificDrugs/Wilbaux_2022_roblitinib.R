Wilbaux_2022_roblitinib <- function() {
  description <- "Two-compartment population PK model for oral roblitinib (FGF401), a selective FGFR4 inhibitor, in adults with hepatocellular carcinoma or other FGF19-FGFR4-expressing solid tumors (Wilbaux 2022). The paper describes a delayed zero-order absorption directly into the central compartment (lag time Tlag before absorption starts, then duration Tk0 of the zero-order input rate) with linear elimination. Categorical covariate effects on CL/F and V1/F for female sex (SEXF) and Asian race (RACE_ASIAN) with non-Asian male as the reference; on Tk0 for fed vs fasted (FED). Continuous power-form covariate effects: body weight on V1/F (exponent 0.332), BMI on Tk0 (exponent -1.66), and administered dose on Tk0 (exponent 0.983). The source paper reports these covariates as identified but not clinically relevant based on simulated exposure metrics. Diagonal IIV on Tlag, Tk0, CL/F, V1/F, and V2/F (Q/F fixed at 0 IIV per Table 1). Combined additive + proportional residual error."
  reference <- "Wilbaux M, Yang S, Jullion A, Demanse D, Graus Porta D, Myers A, Meille C, Gu Y. Integration of Pharmacokinetics, Pharmacodynamics, Safety, and Efficacy into Model-Informed Dose Selection in Oncology First-in-Human Study: A Case of Roblitinib (FGF401). Clin Pharmacol Ther. 2022;112(6):1330-1339. doi:10.1002/cpt.2752. PMID 36131557."
  vignette <- "Wilbaux_2022_roblitinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form scaling on V1/F, exponent 0.332 (Wilbaux 2022 Table 1). Reference weight 70 kg (rounded standard: the source paper does not report a centering value). Time-fixed per subject at baseline.",
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form scaling on Tk0 (duration of zero-order absorption), exponent -1.66 (Wilbaux 2022 Table 1). Reference BMI 25 kg/m^2 (rounded standard: the source paper does not report a centering value). Time-fixed per subject at baseline.",
      source_name        = "BMI"
    ),
    DOSE = list(
      description        = "Administered oral dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record covariate carrying the administered oral dose in mg (per canonical DOSE use case (a): per-subject assigned dose level). Power-form scaling on Tk0, exponent 0.983 (Wilbaux 2022 Table 1). Reference dose 100 mg (rounded mid-range of the tested doses 50/80/120/150 mg: the source paper does not report a centering value). Distinct from the event-table amt column.",
      source_name        = "DOSE"
    ),
    FED = list(
      description        = "Fed vs fasted dosing indicator (1 = fed with light meal within 30 min of dosing, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Wilbaux 2022 Methods, 'Data' section: 'food effect cohorts, where patients received FGF401 doses within 30 minutes following a light meal, were allowed to be opened in parallel in the dose escalation part.' Table 1 reports Tk0 = 0.811 h fasted and Tk0 = 1.58 h fed; encoded as fractional multiplier on Tk0: fed_d1 = 1 + e_fed_d1 * FED, with e_fed_d1 = 1.58/0.811 - 1 = 0.948.",
      source_name        = "FED"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Wilbaux 2022 Table 1 reports separate typical values for female subgroups: CL/F female = 15.1 L/h vs non-Asian male 19.7 L/h; V1/F female = 84.5 L vs non-Asian male 110 L. Encoded as fractional multipliers: e_sexf_cl = 15.1/19.7 - 1 = -0.234 (23.4% lower CL/F in females); e_sexf_v1 = 84.5/110 - 1 = -0.232 (23.2% lower V1/F in females).",
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator (1 = Asian, 0 = non-Asian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Wilbaux 2022 Table 1 reports separate typical values for Asian subgroups: CL/F Asian = 16.3 L/h vs non-Asian male 19.7 L/h; V1/F Asian = 84.5 L vs non-Asian male 110 L. Encoded as fractional multipliers: e_asian_cl = 16.3/19.7 - 1 = -0.173 (17.3% lower CL/F in Asians); e_asian_v1 = 84.5/110 - 1 = -0.232 (23.2% lower V1/F in Asians).",
      source_name        = "RACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 160L,
    n_studies      = 1L,
    disease_state  = "Adults with hepatocellular carcinoma (n = 127) or other FGF19-FGFR4-expressing solid malignancies (n = 33); Phase I/II first-in-human study NCT02325739.",
    dose_range     = "Dose escalation: 50, 80, 120, and 150 mg once daily fasted, and 80 and 120 mg once daily fed (light meal). Expansion part: 120 mg once daily fasted.",
    regions        = "Multicenter international phase I/II study (ClinicalTrials.gov identifier NCT02325739).",
    notes          = "PopPK model was developed on the totality of longitudinal FGF401 plasma-concentration data from all 160 patients across the dose-escalation and dose-expansion cohorts. Detailed baseline demographics (median weight, BMI, age, race distribution) are not tabulated in the extracted paper; the companion Wilbaux 2022 CPT-PSP 11:1122-1134 (doi:10.1002/psp4.12842) reports the detailed popPK development and the Chan 2022 J Exp Clin Cancer Res 41:189 reports the clinical FIH characteristics. Reference values used in this implementation for the continuous covariates are rounded standards (WT 70 kg, BMI 25 kg/m^2, DOSE 100 mg) because the source paper does not report the centering values used at fitting time; see vignette Assumptions and deviations."
  )

  ini({
    # Structural PK parameters at the reference subject: non-Asian male,
    # fasted, WT = 70 kg, BMI = 25 kg/m^2, DOSE = 100 mg. Wilbaux 2022 Table 1
    # reports each fixed effect as a subgroup typical value; the reference
    # typical values below correspond to the "non-Asian male" (CL/F, V1/F)
    # and "fasted condition" (Tk0) rows.
    ltlag <- log(0.268) ; label("Absorption lag time Tlag (h)")                                       # Wilbaux 2022 Table 1: Tlag = 0.268 h, RSE 11%
    ld1   <- log(0.811) ; label("Zero-order absorption duration Tk0 at the reference subject (h)")   # Wilbaux 2022 Table 1: Tk0 fasted = 0.811 h, RSE 10%
    lcl   <- log(19.7)  ; label("Apparent clearance CL/F at the non-Asian male reference (L/h)")     # Wilbaux 2022 Table 1: CL/F non-Asian male = 19.7 L/h, RSE 4%
    lvc   <- log(110)   ; label("Apparent central volume V1/F at 70 kg non-Asian male reference (L)") # Wilbaux 2022 Table 1: V1/F non-Asian male = 110 L, RSE 3%
    lq    <- log(5.59)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                  # Wilbaux 2022 Table 1: Q/F = 5.59 L/h, RSE 6%
    lvp   <- log(49.2)  ; label("Apparent peripheral volume V2/F (L)")                               # Wilbaux 2022 Table 1: V2/F = 49.2 L, RSE 7%

    # Categorical covariate effects (fractional-change form; the reference
    # subgroup has multiplier 1 and the alternate subgroup has multiplier
    # 1 + effect). Coefficients derived from the Wilbaux 2022 Table 1
    # subgroup typical values (non-Asian male / fasted reference).
    e_fed_d1   <-  0.948 ; label("Fractional change in Tk0 for fed vs fasted (unitless)")             # Wilbaux 2022 Table 1: Tk0 fed = 1.58 h vs fasted 0.811 h -> (1.58/0.811 - 1) = 0.948
    e_sexf_cl  <- -0.234 ; label("Fractional change in CL/F for female vs male (unitless)")           # Wilbaux 2022 Table 1: CL/F female = 15.1 vs non-Asian male 19.7 -> (15.1/19.7 - 1) = -0.234
    e_asian_cl <- -0.173 ; label("Fractional change in CL/F for Asian vs non-Asian (unitless)")       # Wilbaux 2022 Table 1: CL/F Asian = 16.3 vs non-Asian male 19.7 -> (16.3/19.7 - 1) = -0.173
    e_sexf_vc  <- -0.232 ; label("Fractional change in V1/F for female vs male (unitless)")           # Wilbaux 2022 Table 1: V1/F female = 84.5 vs non-Asian male 110 -> (84.5/110 - 1) = -0.232
    e_asian_vc <- -0.232 ; label("Fractional change in V1/F for Asian vs non-Asian (unitless)")       # Wilbaux 2022 Table 1: V1/F Asian = 84.5 vs non-Asian male 110 -> (84.5/110 - 1) = -0.232

    # Continuous covariate power exponents. The source paper reports the
    # exponents in Table 1 but does not report the centering (reference)
    # values used during fitting; standard rounded reference values are
    # used here (WT 70 kg, BMI 25 kg/m^2, DOSE 100 mg). See vignette
    # Assumptions and deviations.
    e_bmi_d1  <- -1.66  ; label("Power exponent of BMI on Tk0 (unitless; reference 25 kg/m^2)")      # Wilbaux 2022 Table 1: BMI effect on Tk0 = -1.66, RSE 25%
    e_dose_d1 <-  0.983 ; label("Power exponent of DOSE on Tk0 (unitless; reference 100 mg)")       # Wilbaux 2022 Table 1: Dose effect on Tk0 = 0.983, RSE 28%
    e_wt_vc   <-  0.332 ; label("Power exponent of WT on V1/F (unitless; reference 70 kg)")          # Wilbaux 2022 Table 1: Weight effect on V1/F = 0.332, RSE 33%

    # Inter-individual variability. Wilbaux 2022 Table 1 column header
    # is "IIV: SD of random effect (% RSE)", so the tabulated values are
    # omega (SD on the log scale). Variances = omega^2. No inter-eta
    # correlations reported. Q/F IIV is fixed at 0 in Table 1 ("0 FIX")
    # and is omitted here.
    etaltlag ~ 0.3969  # Wilbaux 2022 Table 1 IIV Tlag omega = 0.63 -> omega^2 = 0.3969
    etald1   ~ 0.4096  # Wilbaux 2022 Table 1 IIV Tk0  omega = 0.64 -> omega^2 = 0.4096
    etalcl   ~ 0.0841  # Wilbaux 2022 Table 1 IIV CL/F omega = 0.29 -> omega^2 = 0.0841
    etalvc   ~ 0.0225  # Wilbaux 2022 Table 1 IIV V1/F omega = 0.15 -> omega^2 = 0.0225
    etalvp   ~ 0.4624  # Wilbaux 2022 Table 1 IIV V2/F omega = 0.68 -> omega^2 = 0.4624

    # Residual error: combined additive + proportional as reported in
    # Wilbaux 2022 Table 1 "Combined residual error: additive (ng.mL^-1) /
    # proportional (%): 4.1 (7) / 33 (2)". Additive is on the ng/mL
    # concentration scale (see Cc unit conversion in model()); the
    # proportional value 33% is the fractional CV.
    addSd  <- 4.1  ; label("Additive residual error (ng/mL)")                                        # Wilbaux 2022 Table 1: additive residual error = 4.1 ng/mL, RSE 7%
    propSd <- 0.33 ; label("Proportional residual error (fraction)")                                 # Wilbaux 2022 Table 1: proportional residual error = 33%, RSE 2%
  })

  model({
    # Categorical covariate multipliers (fractional-change form).
    fed_d1   <- 1 + e_fed_d1   * FED
    sexf_cl  <- 1 + e_sexf_cl  * SEXF
    asian_cl <- 1 + e_asian_cl * RACE_ASIAN
    sexf_vc  <- 1 + e_sexf_vc  * SEXF
    asian_vc <- 1 + e_asian_vc * RACE_ASIAN

    # Continuous covariate power multipliers (rounded-standard centering).
    bmi_d1  <- (BMI  / 25)  ^ e_bmi_d1
    dose_d1 <- (DOSE / 100) ^ e_dose_d1
    wt_vc   <- (WT   / 70)  ^ e_wt_vc

    # Individual PK parameters.
    tlag_i <- exp(ltlag + etaltlag)
    d1_i   <- exp(ld1   + etald1)  * fed_d1  * bmi_d1 * dose_d1
    cl     <- exp(lcl   + etalcl)  * sexf_cl * asian_cl
    vc     <- exp(lvc   + etalvc)  * sexf_vc * asian_vc * wt_vc
    q      <- exp(lq)
    vp     <- exp(lvp   + etalvp)

    # Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Delayed zero-order oral absorption directly into central:
    # the dose is delivered at a constant rate over duration d1_i starting
    # after a lag of tlag_i. Event table must use rate = -2 and
    # cmt = "central" (or the equivalent slot integer) for the oral doses.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(central) <- tlag_i
    dur(central)  <- d1_i

    # Concentration: dose in mg, vc in L -> mg/L; multiply by 1000 to get
    # ng/mL matching the Wilbaux 2022 residual-error units (Table 1
    # additive = 4.1 ng/mL).
    Cc <- (central / vc) * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
