Yuan_2021_busulfan <- function() {
  description <- paste(
    "One-compartment IV population PK model with first-order elimination for",
    "intravenous busulfan in Chinese children (age 0.5-15.2 years) undergoing",
    "allogeneic haematopoietic stem cell transplantation. Clearance carries a",
    "power effect of body surface area (exponent 0.83, referenced to the median",
    "0.67 m^2), a power effect of aspartate aminotransferase (exponent -0.21,",
    "referenced to the median 29.10 U/L), and a log-linear GSTA1 diplotype",
    "effect (the *A/*B group clears busulfan 17.3% more slowly than the *A/*A",
    "reference); volume carries a power effect of body surface area (exponent",
    "0.92). Exponential between-subject variability on CL and V and a combined",
    "additive-plus-proportional residual error. This is the first paediatric",
    "busulfan popPK model to incorporate GSTA1 genotypes in a Chinese",
    "population (Yuan 2021)."
  )
  reference <- paste(
    "Yuan J, Sun N, Feng X, He H, Mei D, Zhu G, Zhao L. Optimization of",
    "Busulfan Dosing Regimen in Pediatric Patients Using a Population",
    "Pharmacokinetic Model Incorporating GST Mutations. Pharmgenomics Pers",
    "Med. 2021;14:253-268. doi:10.2147/PGPM.S289834.",
    sep = " "
  )
  vignette <- "Yuan_2021_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the paper: "Plasma concentrations of Bu
  # were determined using a high performance liquid chromatography-tandem mass
  # spectrometry (HPLC-MS/MS)" (Methods, "Bu Determination and Genotyping"),
  # and busulfan was given as a 2-hour constant-rate IV infusion into the
  # systemic circulation (Methods, "Patients and Treatment Regimens").
  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on both CL (exponent 0.83) and V (exponent 0.92),",
        "referenced to the cohort median 0.67 m^2 (Table 1). Equations 8 and 9",
        "(p. 259): CL = CLpop * (BSA/0.67)^covBSA(CL) * (AST/29.10)^covAST(CL) *",
        "exp(covGSTA1(CL)*GSTA1) * exp(etaCL) and V = Vpop * (BSA/0.67)^covBSA(V)",
        "* exp(etaV). BSA was the most predictive covariate for CL and V,",
        "explaining 25.50% and 24.17% of the observed IIV respectively",
        "(Discussion). Cohort BSA range 0.28-1.50 m^2 (Table 1); the paper",
        "computes BSA from the Stevenson formula BSA = 0.0061*height(cm) +",
        "0.0128*weight(kg) - 0.1529 (Discussion)."
      ),
      source_name = "BSA"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (exponent -0.21), referenced to the cohort median",
        "29.10 U/L (Table 1); Equation 8 (p. 259). AST was negatively correlated",
        "with busulfan CL: 'CL declined 38.34% when AST increased from 12.7 to",
        "127.4 U/L' (Discussion), consistent with (127.4/12.7)^-0.21 = 0.617.",
        "The paper reports the AST covariate in IU/L (Table 1: 29.10 IU/L,",
        "range 12.70-127.40); recorded here under the SI canonical U/L, which",
        "the register documents as interchangeable with IU/L for this analyte.",
        "AST marks hepatocyte lesion, and busulfan is eliminated mainly by the",
        "liver (Discussion)."
      ),
      source_name = "AST"
    ),
    GSTA1_PM = list(
      description = "GSTA1 poor-metabolizer (*A/*B) diplotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GSTA1 *A/*A, the higher-clearance reference group)",
      notes = paste(
        "1 = subject carries the GSTA1 *A/*B diplotype (the poor-metabolizer",
        "group, 17.3% lower busulfan clearance); 0 = GSTA1 *A/*A reference.",
        "The paper genotyped GSTA1 rs3957356 (-52 G>A) and rs3957357 (-69 C>T),",
        "which together define haplotype *A and *B (Methods 'Bu Determination",
        "and Genotyping'); GSTA1-52C and -69T define *B. In Equation 8, GSTA1 =",
        "1 for the *A/*B group and GSTA1 = 0 for the *A/*A group (text below",
        "Equation 9, p. 259). The single *B/*B homozygote and 4 subjects with",
        "missing GSTA1 were excluded, leaving 54 *A/*A and 15 *A/*B in the",
        "n = 69 model-building cohort (Table 2). Because the source resolves",
        "only the *A/*B-vs-*A/*A dichotomy (not the rapid stratum), GSTA1_PM is",
        "used alone -- see the register entry, which prescribes GSTA1_PM alone",
        "for a source that genotypes only the *B haplotype."
      ),
      source_name = "GSTA1 (0 = *A/*A, 1 = *A/*B)"
    )
  )

  # Screened during covariate analysis but not retained in the final model.
  # Documented here so the paper's covariate screen is preserved without
  # carrying "declared but not referenced" convention warnings. The full
  # screen (Methods 'PPK Analysis') tested sex, age, body weight, BSA, ALT,
  # AST, ALP, TBIL, creatinine, creatinine clearance, GSTA1, GSTM1, primary
  # disease and fludarabine co-administration on CL and V; only BSA, AST and
  # GSTA1 survived backward elimination.
  covariatesDataExcluded <- list(
    DIS_MALIGNANT = list(
      description = "Malignant primary-disease indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Types of primary disease (malignant vs non-malignant) entered the",
        "forward step on V but were eliminated in backward elimination",
        "(delta-2LL = 6.333, p > 0.01; Results 'PPK Model'). 26 of 69 subjects",
        "(37.68%) had a malignant diagnosis (Table 1). Name is descriptive only",
        "-- no canonical register entry was created, because the covariate never",
        "enters model()."
      )
    ),
    GSTM1_VARIANT = list(
      description = "GSTM1 A>C (rs3754446) variant indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "GSTM1 genotypes (rs3754446) were tested and had no significant effect",
        "on PK parameters (Results 'PPK Model'), 'likely because the function of",
        "the GSTM1 enzyme involved in Bu metabolism was less than the GSTA1",
        "enzyme' (Discussion). Name is descriptive only -- no canonical register",
        "entry was created, because the covariate never enters model()."
      )
    ),
    CONMED_FLUDARABINE = list(
      description = "Concomitant fludarabine indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Fludarabine co-administration was tested as a candidate covariate and",
        "failed to significantly reduce -2LL during forward selection",
        "(Discussion). Name is descriptive only -- no canonical register entry",
        "was created, because the covariate never enters model()."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 69L,
    n_studies = 1L,
    age_range = "0.50-15.18 years",
    age_median = "4.90 years",
    weight_range = "5.00-48.00 kg",
    weight_median = "16.50 kg",
    sex_female_pct = 46.38,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "Chinese children receiving intravenous busulfan before allogeneic",
      "haematopoietic stem cell transplantation. 37.68% malignant, 62.32%",
      "non-malignant primary diseases (Table 1). Conditioning regimens: Bu+CTX",
      "1.45%, Bu+CTX+FLU 37.68%, Bu+CTX+FLU+VP16 26.09%, Bu+CTX+Ara-C 2.90%,",
      "Bu+CTX+Ara-C+Me-CCNU 30.43%, Bu+CTX+FLU+Dac 1.45%. GSTA1 diplotypes in",
      "the model-building cohort: *A/*A 78.26% (n=54), *A/*B 21.74% (n=15)."
    ),
    dose_range = paste(
      "Intravenous busulfan (Busulfex) as a 2-hour infusion every 6 hours for",
      "three or four days (12 or 16 doses total), first dose 7-9 days before",
      "HSCT. Weight-banded dosing per the EMA regimen: 1.0 mg/kg for < 9 kg,",
      "1.2 mg/kg for 9-16 kg, 1.1 mg/kg for 16-23 kg, 0.95 mg/kg for 23-34 kg,",
      "0.8 mg/kg for > 34 kg. Target AUC0-6h 1125 uM.min (window",
      "900-1350 uM.min)."
    ),
    regions = "China (Beijing Children's Hospital)",
    notes = paste(
      "Prospective cohort collected March 2019 - April 2020; 76 patients",
      "enrolled, 69 included in the final PPK analysis after excluding 4 with",
      "missing GSTA1, 1 GSTA1 *B/*B homozygote and 2 further exclusions.",
      "398 plasma busulfan concentrations. Sampling after the first infusion:",
      "pre-dose and 0.5, 1, 2, 2.5, 4 and 6 h. Assay: HPLC-MS/MS, LLOQ",
      "10 ng/mL, range 10-10000 ng/mL. Estimation: Phoenix NLME 8.0 (Certara),",
      "FOCE ELS. External validation on 81 concentrations from 14 children.",
      "Parameters standardised to median BSA 0.67 m^2, median AST 29.10 U/L",
      "and GSTA1 *A/*A (Results 'PPK Model'). Reported CL/V (per BSA 0.67 m^2)",
      "renormalise to 11.08 L/h per 70 kg for GSTA1 *A/*A (Discussion)."
    )
  )

  ini({
    # ---- Structural parameters --------------------------------------------
    # Table 3 final-model estimates, standardised to the median BSA (0.67 m^2),
    # median AST (29.10 U/L) and the GSTA1 *A/*A reference (Results 'PPK Model';
    # text below Equation 9, p. 259). CLpop and Vpop are therefore the typical
    # values at the covariate reference point, exactly as they enter Equations
    # 8 and 9.
    lcl <- log(4.79); label("Typical busulfan clearance at reference covariates (L/h)")  # Yuan 2021 Table 3 final model CL 4.79 L/h (CV 4.00%; bootstrap 5.09, 95% CI 4.21-6.19)
    lvc <- log(14.80); label("Typical busulfan volume of distribution at reference covariates (L)")  # Yuan 2021 Table 3 final model V 14.80 L (CV 4.02%; bootstrap 15.85, 95% CI 13.28-19.62)

    # ---- Covariate effects on clearance -----------------------------------
    # Equation 8 (p. 259), rendered from the source figure:
    #   CL = CLpop * (BSA/0.67)^covBSA(CL) * (AST/29.10)^covAST(CL)
    #        * exp(covGSTA1(CL) * GSTA1) * exp(etaCL)
    # covBSA(CL) and covAST(CL) are dimensionless power exponents referenced to
    # the cohort medians; covGSTA1(CL) is a log-scale coefficient applied to the
    # GSTA1 *A/*B indicator (1 for *A/*B, 0 for *A/*A). exp(-0.19) = 0.827,
    # i.e. the 17.3% lower CL the paper reports for *A/*B (Results 'PPK Model').
    covbsa_cl <- 0.83; label("Power exponent of body surface area on CL (unitless)")  # Yuan 2021 Table 3 'Cov BSA (CL)' 0.83 (CV 8.72%; bootstrap 0.83, 95% CI 0.67-0.97)
    covast_cl <- -0.21; label("Power exponent of aspartate aminotransferase on CL (unitless)")  # Yuan 2021 Table 3 'Cov AST (CL)' -0.21 (CV -31.55%; bootstrap -0.21, 95% CI -0.34 to -0.09)
    e_gsta1_pm_cl <- -0.19; label("Log-scale CL coefficient for GSTA1 *A/*B vs *A/*A")  # Yuan 2021 Table 3 'Cov GSTA1 (CL)' -0.19 (CV -34.11%; bootstrap -0.19, 95% CI -0.33 to -0.06); exp(-0.19) = 0.827 = 17.3% lower CL

    # ---- Covariate effect on volume ---------------------------------------
    # Equation 9 (p. 259): V = Vpop * (BSA/0.67)^covBSA(V) * exp(etaV).
    covbsa_vc <- 0.92; label("Power exponent of body surface area on V (unitless)")  # Yuan 2021 Table 3 'Cov BSA (V)' 0.92 (CV 9.14%; bootstrap 0.91, 95% CI 0.75-1.08)

    # ---- Between-subject variability --------------------------------------
    # Table 3 reports IIV as %CV under an exponential (log-normal) random-effect
    # model, P_i = P * exp(eta_i) (Methods 'PPK Analysis'). Converted to the
    # internal log-scale variance with omega^2 = log(CV^2 + 1):
    #   IIV CL 18.65% -> log(0.1865^2 + 1) = 0.0341910
    #   IIV V  23.63% -> log(0.2363^2 + 1) = 0.0543345
    etalcl ~ 0.0341910  # Yuan 2021 Table 3 final model 'omega CL' 18.65% (CV 28.08%; shrinkage 0.207)
    etalvc ~ 0.0543345  # Yuan 2021 Table 3 final model 'omega V' 23.63% (CV 25.42%; shrinkage 0.130)

    # ---- Residual error ----------------------------------------------------
    # Combined additive + proportional residual error, Equation 3 (p. 258):
    # OBS = PRED * (1 + eps2) + eps1. Table 3 final model: eps1 (additive) =
    # 0.043 ug/mL, eps2 (proportional) = 7.8%. Concentrations are in ug/mL
    # (= mg/L), so the additive term is 0.043 mg/L. Neither the proportional
    # nor the additive error model alone fit well; only the combined model gave
    # an adequate fit (Results 'PPK Model').
    addSd <- 0.043; label("Additive residual error (mg/L)")  # Yuan 2021 Table 3 final model 'eps1' 0.043 ug/mL
    propSd <- 0.078; label("Proportional residual error (fraction)")  # Yuan 2021 Table 3 final model 'eps2' 7.8%
  })

  model({
    # ---- 1. Individual parameters -----------------------------------------
    # Equations 8 and 9 (p. 259). BSA and AST enter as dimensionless power
    # terms referenced to the cohort medians (0.67 m^2 and 29.10 U/L); the
    # GSTA1 *A/*B indicator enters CL log-linearly. The reference subject
    # (BSA 0.67, AST 29.10, GSTA1 *A/*A) reduces every covariate factor to 1,
    # recovering CLpop and Vpop.
    cl <- exp(lcl + etalcl) *
      (BSA / 0.67)^covbsa_cl *
      (AST / 29.10)^covast_cl *
      exp(e_gsta1_pm_cl * GSTA1_PM)
    vc <- exp(lvc + etalvc) * (BSA / 0.67)^covbsa_vc

    # ---- 2. Micro-constants -----------------------------------------------
    kel <- cl / vc

    # ---- 3. ODE system ----------------------------------------------------
    # One compartment with first-order elimination (Results 'PPK Model': "A
    # one-compartment model with first-order elimination best described the
    # data"). Busulfan is infused directly into the systemic circulation, so
    # there is no depot and no bioavailability term.
    d/dt(central) <- -kel * central

    # ---- 4. Observation and error -----------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
