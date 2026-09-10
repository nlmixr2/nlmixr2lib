Huang_2025_colistinSulfate <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous colistin sulfate in",
    "critically ill adults receiving continuous veno-venous hemodiafiltration",
    "(CVVHDF) for acute kidney injury (Huang 2025; n = 20 Chinese ICU",
    "patients, 86 plasma concentrations spanning 0.09-2.32 mg/L). Linear",
    "elimination from the central compartment with a 0.5-h or 2-h",
    "intravenous infusion input. Serum cystatin C and body weight both enter",
    "clearance as power functions centred on the cohort medians 2.31 mg/L",
    "and 65 kg (exponents -0.626 and 1.03); clearance falls as cystatin C",
    "rises (worse residual renal function) and rises with body weight.",
    "Cockcroft-Gault creatinine clearance was screened and NOT retained --",
    "the authors attribute this to creatinine (113 Da) being efficiently",
    "removed by CVVHDF while cystatin C (13.25 kDa) is not, so only",
    "cystatin C still tracks residual renal function on dialysis.",
    "Inter-individual variability was estimated on CL and V2 only (IIV on V1",
    "and Q collapsed toward zero with shrinkage > 30%); residual error is",
    "exponential. Fixed allometric scaling was tested and rejected as not",
    "supported by the data. Colistin sulfate is administered as the active",
    "drug and must not be confused with colistimethate sodium (CMS), the",
    "inactive prodrug modelled in Plachouras 2009, Mohamed 2012,",
    "Jacobs 2016 and Karaiskos 2015. Dose unit conversion: 1 million units",
    "(MU) = 44 mg."
  )
  reference <- paste(
    "Huang T, Luo Y, Wu Y, Niu L, Xiao Y, Wu T, Chen X, Liu Y, Lu J,",
    "Zhu D, Liu T (2025).",
    "Population pharmacokinetics of colistin sulfate in patients on",
    "continuous veno-venous hemodiafiltration.",
    "Science Progress 108(1):1-20.",
    "doi:10.1177/00368504251325334.",
    sep = " "
  )
  vignette <- "Huang_2025_colistinSulfate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CYSC = list(
      description        = paste(
        "Serum cystatin C. Retained as the sole renal-function covariate on",
        "clearance. Cockcroft-Gault creatinine clearance (CrCL), serum",
        "creatinine, blood urea nitrogen and uric acid were all screened and",
        "none reached the forward-inclusion criterion on CL. The authors'",
        "mechanistic explanation (Discussion) is that CVVHDF combines",
        "convection and diffusion and therefore clears small solutes",
        "efficiently: creatinine (113 Da) is removed by the filter, so serum",
        "creatinine and any CrCL derived from it no longer reflect residual",
        "renal function, whereas cystatin C (13.25 kDa) is not appreciably",
        "cleared (cited clearance 17 mL/min, under 30% of its production",
        "rate) and so remains a usable marker of residual renal function",
        "during CRRT."
      ),
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL centred on 2.31 mg/L, the cohort MEDIAN, stated",
        "in the Results text immediately below the final-model equation",
        "('The values 2.31 and 65 were the medians for CYSC and WT'). This",
        "is distinct from the Table 1 cohort MEAN of 2.87 +/- 1.91 mg/L --",
        "do not substitute the mean. Observed range 1.05-5.11 mg/L (Figure 5",
        "note, which lists the 5th/25th/50th/75th/95th percentiles used for",
        "the dosing simulations). The exponent is negative (-0.626), so",
        "clearance FALLS as cystatin C rises; the paper states this",
        "explicitly. This is the first registered model in which cystatin C",
        "is the renal covariate in a continuous-renal-replacement-therapy",
        "population."
      ),
      source_name        = "CysC"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL centred on 65 kg, the cohort MEDIAN (Results",
        "text below the final-model equation), NOT the Table 1 cohort mean",
        "of 68.6 +/- 14.1 kg and NOT the conventional 70 kg. The Methods",
        "describe a FIXED allometric model 'P = theta_p * (WT/70)^K' with K",
        "fixed at 0.75 on CL/Q and 1 on V1/V2, but Results state that model",
        "'was not suitable for our data, so it was not considered finally'",
        "(Table 2 confirms: dOFV 958.59 -> 954.37, p > 0.05, Reserve = NO).",
        "The exponent 1.03 in the final model is therefore ESTIMATED (RSE",
        "21.7%), applies to CL only, and is centred on 65 -- not the",
        "textbook fixed 0.75 on 70 kg. Simulation range 50-80 kg (Figure 5",
        "note). Weight does not enter V1, V2 or Q in the final model."
      ),
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    # Methods "Quantification of colistin sulfate concentrations": samples
    # were centrifuged and "the PLASMA concentrations of colistin sulfate
    # were measured" by UHPLC-MS/MS. Colistin A and colistin B were
    # quantified separately (LOQ 0.027 and 0.053 mg/L) and, because they
    # share structure, molecular weight, activity and PK, "the plasma
    # concentration of colistin was derived by summing the concentrations of
    # colistin A and B". The assayed analyte is therefore active colistin
    # itself, not a colistimethate-derived metabolite.
    central     = list(
      analyte  = "colistin sulfate (sum of colistin A and colistin B)",
      units    = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "colistin sulfate (sum of colistin A and colistin B)",
      units    = "mg", specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    n_observations = 86,
    age_mean       = "50.5 +/- 14.1 years",
    weight_mean    = "68.6 +/- 14.1 kg",
    weight_median  = "65 kg",
    weight_range   = "50-80 kg (simulation range, Figure 5 note)",
    sex_female_pct = 25,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill adults with confirmed or suspected carbapenem-resistant",
      "organism (CRO) infection, all receiving CVVHDF for acute kidney injury",
      "for at least 48 h. Severe illness: SOFA 9.27 +/- 4.08, APACHE II",
      "26.6 +/- 10.1; 90% on vasoactive agents, 90% mechanically ventilated,",
      "90% with lung infection and 75% with multi-site infection. Pathogens",
      "were all CRO (meropenem MIC >= 8 mg/L): Klebsiella pneumoniae 35%,",
      "Acinetobacter baumannii 30%, Pseudomonas aeruginosa 30%,",
      "Escherichia coli 10%."
    ),
    renal_function = paste(
      "All subjects on CVVHDF (Prismaflex, 1.5 m2 AN69-ST150 polyacrylonitrile",
      "filter) with blood flow 139.3 +/- 29.1 mL/min, dialysate flow",
      "1019.8 +/- 262.9 mL/h, replacement fluid 1002.2 +/- 332.7 mL/h,",
      "ultrafiltration 109.9 +/- 62.1 mL/h and 25-35 mL/kg/h dose intensity.",
      "Residual diuresis 0.17 +/- 0.25 mL/h/kg. Cystatin C 2.87 +/- 1.91 mg/L",
      "(median 2.31), serum creatinine 140.3 +/- 66.9 umol/L, Cockcroft-Gault",
      "CrCL 59.7 +/- 25.9 mL/min. No CVVHDF flow parameter (BFR, DFR, RFR,",
      "UFR) reached significance on any PK parameter; the authors attribute",
      "this to a plateau effect at the high flow rates used."
    ),
    dose_range     = paste(
      "1.0-2.0 MU daily intravenously (1 MU = 44 mg, i.e. 44-88 mg/day),",
      "divided q8h (85%) or q12h (15%); 1.5 MU daily in 85% of subjects.",
      "Infusion duration 0.5 h (90%) or 2 h (10%). A loading dose (twice the",
      "maintenance dose) was given to 50%. Three subjects with lung infection",
      "additionally received nebulized colistin sulfate."
    ),
    regions        = "China (single centre, Nanning, Guangxi)",
    notes          = paste(
      "Prospective single-centre observational study, May 2023 - January 2024",
      "(ChiCTR2300072191). Baseline demographics in Table 1. Sampling was",
      "rich within a maintenance-phase dosing interval: a pre-dose trough",
      "plus four to six post-dose samples at nominally 0.5, 1, 2, 4 and 6 h.",
      "Observed Css,min 0.30 +/- 0.22 mg/L and trapezoidal AUCss,0-24h",
      "12.51 +/- 6.41 mg*h/L (obtainable for only 16 of 20 subjects).",
      "Exclusions: no blood samples available, death within 72 h of starting",
      "colistin sulfate, or ECMO."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Huang 2025 Table 3 "Final model / Estimate",
    # cross-checked against the final-model equation block printed in the
    # Results text (publisher equation image `...-eq2.jpg`, which renders
    # the same four lines: CL, V1, V2, Q).
    #
    # A two-compartment model with first-order linear elimination was chosen
    # over one-compartment (Results "Population PK model"; supplement
    # Table S1 "Structure model" compares the two).
    # ----------------------------------------------------------------------
    lcl <- log(3.69);  label("Clearance at the reference covariates, CYSC = 2.31 mg/L and WT = 65 kg (L/h)")  # Table 3 CL = 3.69 L/h (RSE 22.2%; bootstrap median 3.18, 95% CI 1.14-5.34)
    lvc <- log(20.5);  label("Central volume of distribution, V1 (L)")                                        # Table 3 V1 = 20.50 L (RSE 10.1%; bootstrap median 20.58, 95% CI 14.44-27.73)
    lvp <- log(33.2);  label("Peripheral volume of distribution, V2 (L)")                                      # Table 3 V2 = 33.20 L (RSE 17.3%; bootstrap median 34.71, 95% CI 24.44-50.87)
    lq  <- log(25.3);  label("Intercompartmental clearance, Q (L/h)")                                          # Table 3 Q = 25.30 L/h (RSE 26.5%; bootstrap median 27.15, 95% CI 11.78-40.66)

    # ----------------------------------------------------------------------
    # Covariate effects on CL -- BOTH are power functions.
    #
    # ERRATUM / wording conflict, resolved in favour of the printed equation:
    # the Results narrative states "CYSC and WT are included on CL in power
    # form and exponential form, respectively", but the paper's own
    # final-model equation is
    #     CL(L/h) = 3.69 * (CYSC/2.31)^-0.626 * (WT/65)^1.03
    # in which WT is unambiguously a POWER term, not an exponential one.
    # The equation is confirmed character-for-character from the publisher's
    # native equation rendering (EuropePMC supplementaryFiles
    # `10.1177_00368504251325334-eq2.jpg`), so this is a printed value and
    # not a digitisation. An exponential reading is also arithmetically
    # impossible: exp(1.03 * WT/65) would multiply CL by 2.80 at the
    # reference weight instead of 1, contradicting the stated typical CL of
    # 3.69 L/h, and exp(1.03 * (WT - 65)) would inflate CL by 5e6-fold
    # across the 50-80 kg simulation range. The word "exponential" is a
    # narrative slip for "exponent". Power form is used here.
    # ----------------------------------------------------------------------
    e_cysc_cl <- -0.626; label("Power exponent on (CYSC / 2.31 mg/L) for CL (unitless)")  # Table 3 theta_CysC = -0.626 (RSE 16.8%; bootstrap median -0.537, 95% CI -0.824 to -0.222)
    e_wt_cl   <-  1.03;  label("Power exponent on (WT / 65 kg) for CL (unitless)")        # Table 3 theta_wt = 1.030 (RSE 21.7%; bootstrap median 1.189, 95% CI 0.761-2.30). ESTIMATED, not the fixed 0.75 allometric value -- the fixed-allometry model was rejected (Table 2, p > 0.05).

    # ----------------------------------------------------------------------
    # Inter-individual variability -- on CL and V2 ONLY.
    #
    # Results "Population PK model": IIV on V1 and Q "yielded values nearing
    # zero, with shrinkage values surpassing 30%" and did not significantly
    # change the OFV, so "interindividual variability was incorporated solely
    # into CL and V2". Supplement Table S1 records the same sequence
    # (IIV(CL+V1) -> "IIV on V1 was close to zero and shrinkage 99.8%";
    # IIV(CL+Q) -> "shrinkage of IIV-Q 37.9%").
    #
    # SCALE: Table 3 reports these as percentages in the "Final model /
    # Estimate" column and as bare numbers in the "Bootstrap / Median"
    # column -- 31% vs 0.311, 54.5% vs 0.557, 20.3% vs 0.193. The bare
    # numbers are omega STANDARD DEVIATIONS, not variances: 0.311 reproduces
    # "31%" directly, whereas reading 0.311 as a variance would imply a
    # 60.4% CV, which is double the printed figure and would leave the
    # bootstrap median inconsistent with a final estimate the paper
    # describes as "close". nlmixr2's ini() takes VARIANCES, so each SD is
    # squared here.
    # ----------------------------------------------------------------------
    etalcl ~ 0.31^2   # Table 3 IIV_CL  = 31%   as an omega SD (RSE 18.0%, eta-shrinkage 1.6%;  bootstrap median 0.311, 95% CI 0.180-0.426) -> variance 0.0961
    etalvp ~ 0.545^2  # Table 3 IIV_V2  = 54.5% as an omega SD (RSE 20.4%, eta-shrinkage 28.3%; bootstrap median 0.557, 95% CI 0.214-0.893) -> variance 0.297025

    # ----------------------------------------------------------------------
    # Residual error -- EXPONENTIAL.
    #
    # Methods: "Residual variability was described by additive,
    # proportional, exponential, or mixed models"; Results: covariate-model
    # IIV and "residual variability modeled using an exponential model".
    # Supplement Table S1 confirms exponential RSV was carried from the
    # structural-model step onward, with the additive form failing
    # ("NONMEM Running Error") and proportional/mixed forms not retained.
    # An exponential residual C_obs = C_pred * exp(eps) is lnorm() in
    # nlmixr2, i.e. an additive SD on the log-concentration scale.
    # ----------------------------------------------------------------------
    expSd <- 0.203; label("Residual SD on the log-transformed concentration scale (exponential error)")  # Table 3 RSV = 20.3% as a log-scale SD (RSE 15.6%, epsilon-shrinkage 16.1%; bootstrap median 0.193, 95% CI 0.127-0.263)
  })

  model({
    # 1. Individual parameters. Covariates enter CL only; V1, V2 and Q carry
    #    no covariate in the final model.
    cl <- exp(lcl + etalcl) * (CYSC / 2.31)^e_cysc_cl * (WT / 65)^e_wt_cl
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Intravenous infusion is dosed directly into `central`
    #    (0.5 h in 90% of subjects, 2 h in the remainder).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observation and error
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
