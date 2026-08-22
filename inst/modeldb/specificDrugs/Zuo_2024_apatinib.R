Zuo_2024_apatinib <- function() {
  description <- paste0(
    "One-compartment population PK model for oral apatinib in Chinese adult patients with solid tumours ",
    "(Zuo 2024), developed from steady-state trough therapeutic-drug-monitoring samples. First-order ",
    "absorption with Ka fixed at 0.08 1/h (the absorption phase was not sampled), a power effect of ",
    "aspartate aminotransferase on apparent clearance ((AST/26.6 U/L)^-0.298), and a four-level ",
    "concomitant-antineoplastic-regimen multiplier on CL/F (reference: apatinib plus an ",
    "immune-checkpoint-inhibitor monoclonal antibody; paclitaxel 0.58, other cytotoxic agents 1.60, ",
    "apatinib monotherapy 1.38). Inter-individual variability is estimated on CL/F only; residual error ",
    "is combined additive plus proportional."
  )
  reference <- paste(
    "Zuo L, Ling J, Hu N, Chen R. (2024).",
    "Establishment and validation of a population pharmacokinetic model for apatinib in patients with tumors.",
    "BMC Cancer 24:1338. doi:10.1186/s12885-024-13118-4"
  )
  vignette <- "Zuo_2024_apatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zuo 2024 "Method of analyzing patients'
  # blood samples": apatinib was quantified in SERUM by LC-MS/MS, and the only
  # dosed route is oral (250 mg tablets), so `depot` is the gut lumen.
  compartmentData <- list(
    depot   = list(analyte = "apatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apatinib", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Serum aspartate aminotransferase, a marker of hepatocellular injury / liver function",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at the biochemistry panel drawn within 1 week before the apatinib dose (Zuo 2024",
        "inclusion criteria). Enters CL/F as a median-normalised power term (AST/26.6)^-0.298; 26.6 U/L",
        "is the study-population median (Table 1, range 10.5-187 U/L). Apatinib is cleared primarily by",
        "hepatic CYP3A4/5, so higher AST (worse hepatocellular function) maps to lower CL/F.",
        "AST was the only biochemical covariate retained after backward elimination."
      ),
      source_name        = "AST"
    ),
    CONMED_PACLITAXEL = list(
      description        = "Concomitant paclitaxel indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (apatinib co-administered with an immune-checkpoint-inhibitor monoclonal antibody)",
      notes              = paste(
        "1 = the apatinib sample was drawn while the patient was co-treated with paclitaxel (Zuo 2024",
        "co-administration group B; 29 of 189 classified samples). Multiplicative effect on CL/F:",
        "CL/F is multiplied by 0.58 (i.e. 42% lower) relative to the mAb reference group.",
        "One of three mutually exclusive indicators; all three zero identifies the reference group."
      ),
      source_name        = "CM (category B)"
    ),
    CONMED_ANTINEO_OTHER = list(
      description        = "Concomitant non-taxane cytotoxic antineoplastic indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (apatinib co-administered with an immune-checkpoint-inhibitor monoclonal antibody)",
      notes              = paste(
        "1 = the apatinib sample was drawn while the patient was co-treated with a cytotoxic",
        "antineoplastic other than paclitaxel - platinum agents, capecitabine, or the",
        "tegafur/gimeracil/oteracil potassium combination (S-1) (Zuo 2024 co-administration group C;",
        "43 of 189 classified samples). Multiplicative effect on CL/F: CL/F is multiplied by 1.60",
        "(i.e. 60% higher) relative to the mAb reference group.",
        "One of three mutually exclusive indicators; all three zero identifies the reference group."
      ),
      source_name        = "CM (category C)"
    ),
    CONMED_ANTINEO_NONE = list(
      description        = "Apatinib monotherapy indicator (no concomitant antineoplastic)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (apatinib co-administered with an immune-checkpoint-inhibitor monoclonal antibody)",
      notes              = paste(
        "1 = the apatinib sample was drawn while the patient was on apatinib monotherapy, with no",
        "concomitant antineoplastic agent (Zuo 2024 co-administration group D; 6 of 189 classified",
        "samples - the smallest group, and its 1.38 multiplier carries the largest RSE in Table 2 at",
        "38%). Multiplicative effect on CL/F: CL/F is multiplied by 1.38 relative to the mAb reference",
        "group. One of three mutually exclusive indicators; all three zero identifies the reference group."
      ),
      source_name        = "CM (category D)"
    )
  )

  # Covariates Zuo 2024 screened through stepwise covariate modelling but did
  # NOT retain in the final model. TBILI reached the full covariate model in
  # forward inclusion but was dropped in backward exclusion (the OFV increase on
  # its removal was < 10.83); every other name below failed forward inclusion
  # (delta OFV < 3.84). Documentation only - none of these are referenced in
  # model(). Screened covariates with no canonical register column (serum total
  # bile acids, globulin, adenosine deaminase, cholinesterase, uric acid, anion
  # gap, prealbumin, tumour type) are recorded in population$notes instead of
  # being minted as unregistered column names.
  covariatesDataExcluded <- list(
    AGE   = list(description = "Subject age", units = "years", type = "continuous",
                 notes = "Screened; not retained. Median 64 years (range 27-86), Zuo 2024 Table 1."),
    WT    = list(description = "Body weight", units = "kg", type = "continuous",
                 notes = "Screened; not retained. Median 55 kg (range 39-85), Zuo 2024 Table 1."),
    SEXF  = list(description = "Female sex indicator", units = "(binary)", type = "binary",
                 notes = "Screened; not retained. 31 of 91 patients female, Zuo 2024 Results."),
    ALT   = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                 notes = "Screened; not retained. Median 17.4 U/L (range 5.1-124.3), Zuo 2024 Table 1."),
    GGT   = list(description = "Gamma-glutamyltransferase", units = "U/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported."),
    ALP   = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported."),
    LDH   = list(description = "Lactate dehydrogenase", units = "U/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported."),
    TBILI = list(description = "Total bilirubin", units = "umol/L", type = "continuous",
                 notes = paste("Screened; entered the full covariate model in forward inclusion together with",
                               "AST and CM, then dropped in backward exclusion (OFV increase < 10.83).",
                               "No point estimate is published, so no effect can be encoded."),
                 source_name = "TBIL"),
    TPRO  = list(description = "Total serum protein", units = "g/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported."),
    ALB   = list(description = "Serum albumin", units = "g/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported."),
    CYSC  = list(description = "Serum cystatin C", units = "mg/L", type = "continuous",
                 notes = "Screened; not retained. Median 0.99 mg/L (range 0.43-2.98), Zuo 2024 Table 1."),
    CREAT = list(description = "Serum creatinine", units = "umol/L", type = "continuous",
                 notes = "Screened; not retained. Median 63 umol/L (range 35-207), Zuo 2024 Table 1.",
                 source_name = "Scr"),
    BUN   = list(description = "Blood urea", units = "mmol/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported.",
                 source_name = "UREA"),
    CRCL  = list(description = "Estimated glomerular filtration rate (BSA-normalised)",
                 units = "mL/min/1.73 m^2", type = "continuous",
                 notes = "Screened; not retained. Median 95.32 (range 14.08-144.48), Zuo 2024 Table 1.",
                 source_name = "eGFR"),
    GLU   = list(description = "Blood glucose", units = "mmol/L", type = "continuous",
                 notes = "Screened; not retained. No summary statistics reported.")
  )

  population <- list(
    species        = "human",
    n_subjects     = 91L,
    n_studies      = 1L,
    n_observations = 199L,
    age_range      = "27-86 years",
    age_median     = "64 years",
    weight_range   = "39-85 kg",
    weight_median  = "55 kg",
    sex_female_pct = 34.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Adult inpatients with solid tumours receiving oral apatinib as part of routine care;",
      "hepatic function spans normal to markedly abnormal (AST 10.5-187 U/L, ALT 5.1-124.3 U/L)",
      "and renal function spans normal to severely impaired (eGFR 14.08-144.48 mL/min/1.73 m^2,",
      "serum creatinine 35-207 umol/L)."
    ),
    dose_range     = "Apatinib 250 mg orally once daily, as monotherapy or in combination.",
    regions        = "China (single centre: the Third Affiliated Hospital of Soochow University, Changzhou).",
    sampling       = paste(
      "Sparse therapeutic drug monitoring. All patients were at steady state and every sample was a",
      "pre-dose trough drawn immediately before the next administration; on average two samples per",
      "patient. No absorption-phase or post-peak data were collected, which is why Ka was fixed and",
      "no inter-individual variability could be estimated on V/F or Ka."
    ),
    co_medication  = paste(
      "Four mutually exclusive co-administration groups classified by the drug given with apatinib",
      "(Zuo 2024 Table 1 footnote, 111:29:43:6 samples in groups A:B:C:D): A = immune-checkpoint-inhibitor",
      "monoclonal antibodies (camrelizumab, sintilimab, etc.); B = paclitaxel; C = other agents",
      "(platinum, capecitabine, or tegafur/gimeracil/oteracil potassium); D = apatinib monotherapy."
    ),
    bioanalytical  = paste(
      "LC-MS/MS (AB SCIEX Triple Quad 4500MD with a Jasper HPLC), duloxetine internal standard,",
      "apatinib transition m/z 398.1 -> 211.9. Calibrated range 1.00-500 ug/L; intra-day CV 2.79-5.02%",
      "and inter-day CV 6.20-9.94%."
    ),
    notes          = paste(
      "Baseline demographics from Zuo 2024 Table 1. Data were collected between October 2021 and",
      "December 2023. Covariates screened by stepwise covariate modelling but not retained were WT,",
      "AGE, SEX, tumour type (CANCER), ALT, GGT, ALP, LDH, serum total bile acids (TBA), total protein",
      "(TP), albumin (ALB), globulin (GLB), total bilirubin (TBIL), adenosine deaminase (ADA),",
      "cholinesterase (CHE), urea (UREA), cystatin C (CYSC), serum creatinine (SCR), uric acid (UA),",
      "anion gap (AG), blood glucose (GLU), prealbumin (PA), and eGFR. Only AST and the",
      "co-administration category survived. The 111:29:43:6 co-administration split sums to 189",
      "samples, 10 fewer than the 199 total, consistent with the 10 records the authors held out for",
      "external validation (MPE -2.80%, MAE 18.75% for the final model)."
    )
  )

  ini({
    # Structural PK -- Zuo 2024 Table 2 final-model estimates and the boxed
    # final-model equation printed immediately after Table 2. Apparent clearance
    # CL/F in L/h and apparent volume V/F in L (bioavailability is not separately
    # identifiable from oral-only data, so F is folded into both). Concentrations
    # are reported in ng/mL after the 1000x conversion from internal mg/L
    # (dose mg, vc L) inside model().
    #
    # Ka was NOT estimated: "The data from the absorption phase were not
    # available; thus, Ka was fixed according to the initiate information [10]
    # available in order to avoid the effect of Ka fluctuation on the stability
    # of model, and its interindividual variability was not estimated."
    lka <- fixed(log(0.08)) ; label("First-order absorption rate constant Ka (1/h)")    # Zuo 2024 Abstract and final-model equation box: "Ka = 0.08" (fixed)
    lcl <- log(56.7)        ; label("Apparent clearance CL/F at AST 26.6 U/L, mAb co-administration (L/h)")  # Zuo 2024 Table 2 typical CL/F = 56.7 L/h (RSE 33%)
    lvc <- log(674)         ; label("Apparent volume of distribution V/F (L)")          # Zuo 2024 Table 2 typical V/F = 674 L (RSE 34%)

    # Covariate effects on CL/F -- Zuo 2024 final-model equation box:
    #   CL/F = 56.7 * (AST/26.6)^-0.298 * theta
    # The AST exponent is taken from the boxed equation (-0.298); Table 2 prints
    # the same estimate rounded to -0.30 (RSE 17%, bootstrap 95% CI -0.38 to -0.21).
    e_ast_cl <- -0.298 ; label("AST power exponent on CL/F (AST/26.6 U/L; unitless)")   # Zuo 2024 final-model equation box, exponent -0.298 (Table 2 rounds to -0.30)

    # Co-administration category multipliers (theta) on CL/F. The reference
    # category is "combined with monoclonal antibodies", for which Table 2 fixes
    # theta = 1 (printed as "1" with no RSE); it therefore has no free parameter
    # here and is recovered when all three indicators are 0.
    e_conmed_paclitaxel_cl    <- 0.58 ; label("CL/F multiplier for concomitant paclitaxel (vs mAb reference; unitless)")            # Zuo 2024 Table 2 "Combined with paclitaxel" = 0.58 (RSE 9%)
    e_conmed_antineo_other_cl <- 1.60 ; label("CL/F multiplier for other concomitant cytotoxics (vs mAb reference; unitless)")      # Zuo 2024 Table 2 "Combined with other drugs" = 1.60 (RSE 27%)
    e_conmed_antineo_none_cl  <- 1.38 ; label("CL/F multiplier for apatinib monotherapy (vs mAb reference; unitless)")              # Zuo 2024 Table 2 "Monotherapy" = 1.38 (RSE 38%)

    # Inter-individual variability -- exponential IIV model (Zuo 2024 Eq. 1,
    # P_i = P_pop * exp(eta_i)); Table 2 reports it as "Interindividual
    # variation, %" = 19.57% for CL/F. Converting the lognormal CV to a
    # log-scale variance: omega^2 = log(0.1957^2 + 1) = 0.037583.
    # (Reading 19.57% instead as omega itself would give omega^2 = 0.038299,
    # a 0.9% difference on the SD scale; the choice is immaterial here.)
    # IIV on V/F and Ka was NOT estimated -- only steady-state troughs were
    # collected, so neither is identifiable (Zuo 2024 Modeling, Discussion).
    etalcl ~ 0.037583                                                                    # Zuo 2024 Table 2 IIV on CL/F = 19.57% CV (RSE 47%)

    # Residual unexplained variability -- Zuo 2024 selected the combined
    # (additive + proportional) error model (Eq. 2) over additive, proportional
    # and exponential alternatives. Table 2 reports the proportional term as a
    # percentage and the additive term in ng/mL, i.e. both on the SD scale.
    propSd <- 0.3067 ; label("Proportional residual error (fraction)")                    # Zuo 2024 Table 2 residual proportional = 30.67% (RSE 33%)
    addSd  <- 16.4   ; label("Additive residual error (ng/mL)")                           # Zuo 2024 Table 2 residual additive = 16.4 ng/mL (RSE 44%)
  })

  model({
    # 1. Derived covariate terms.
    # AST effect on CL/F: median-normalised power model (Zuo 2024 Eq. 4 applied
    # with COVmedian = 26.6 U/L, the Table 1 study median).
    ast_cl <- (AST / 26.6)^e_ast_cl

    # Co-administration effect on CL/F. The three indicators are mutually
    # exclusive, so at most one exponent is 1 and the product collapses to the
    # single applicable theta; all three zero gives 1, the mAb reference group.
    conmed_cl <-
      e_conmed_paclitaxel_cl^CONMED_PACLITAXEL *
      e_conmed_antineo_other_cl^CONMED_ANTINEO_OTHER *
      e_conmed_antineo_none_cl^CONMED_ANTINEO_NONE

    # 2. Individual PK parameters. IIV is on CL/F only.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * ast_cl * conmed_cl
    vc <- exp(lvc)

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. One-compartment oral PK with first-order absorption and first-order
    # elimination. F is not separately identifiable from oral-only data and is
    # absorbed into CL/F and V/F, so no f(depot) term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation. Apatinib serum concentration: dose mg, vc L -> mg/L;
    # multiply by 1000 to report ng/mL, the unit of the Table 2 additive
    # residual error and of the LC-MS/MS calibration range (1.00-500 ug/L).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
