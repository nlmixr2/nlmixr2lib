Tuey_2024_cholecalciferol <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral cholecalciferol",
    "(vitamin D3) and its three major metabolites in 29 adults with chronic",
    "kidney disease and vitamin D insufficiency or deficiency (Tuey 2024).",
    "Cholecalciferol follows two-compartment disposition with first-order",
    "absorption plus a constant zero-order endogenous production rate;",
    "25-hydroxyvitamin D3 (25D3), 1,25-dihydroxyvitamin D3 (1,25D3) and",
    "24,25-dihydroxyvitamin D3 (24,25D3) are each one-compartment models with",
    "first-order formation and first-order elimination, chained in sequence",
    "VitD3 -> 25D3 -> {1,25D3, 24,25D3}. Every species carries an estimated",
    "endogenous pre-dose baseline concentration used as the compartment",
    "initial condition. The formed fractions are fixed for identifiability",
    "(all VitD3 clearance forms 25D3; 1.7 percent of 25D3 clearance forms",
    "1,25D3 and the remainder forms 24,25D3). Concentrations are on a molar",
    "scale (nmol/L) so the four analytes could be fitted in one data set.",
    "Proportional residual error on each analyte. No covariate reached",
    "statistical significance."
  )
  reference <- paste(
    "Tuey SM, Ghimire A, Guzy S, Prebehalla L, Roque AA, Roda G,",
    "West RE 3rd, Chonchol MB, Shah N, Nolin TD, Joy MS.",
    "Population Pharmacokinetic Model of Vitamin D3 and Metabolites in",
    "Chronic Kidney Disease Patients with Vitamin D Insufficiency and",
    "Deficiency. Int J Mol Sci. 2024;25(22):12279.",
    "doi:10.3390/ijms252212279.",
    sep = " "
  )
  vignette <- "Tuey_2024_cholecalciferol"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All four analytes were measured in plasma by
  # UHPLC-MS/MS (Tuey 2024 Methods 4.1 and 4.4); amounts are in nmol
  # because Methods 4.5 converted every observation from ng/mL to nmol/L
  # before fitting.
  compartmentData <- list(
    depot           = list(analyte = "cholecalciferol", units = "nmol", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "cholecalciferol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "cholecalciferol", units = "nmol", specimen = "plasma", verified = TRUE),
    central_25d3    = list(analyte = "25-hydroxyvitamin D3", units = "nmol", specimen = "plasma", verified = TRUE),
    central_125d3   = list(analyte = "1,25-dihydroxyvitamin D3", units = "nmol", specimen = "plasma", verified = TRUE),
    central_2425d3  = list(analyte = "24,25-dihydroxyvitamin D3", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Tuey 2024 Methods 4.7 screened a broad covariate panel by visual
  # inspection of IIV-versus-covariate plots followed by stepwise forward
  # addition (p < 0.05) and backward elimination (p < 0.01). Results 2.3:
  # "No covariates evaluated led to significant influences on parent and
  # metabolite pharmacokinetics". None is therefore referenced in model(),
  # and no point estimate is published for any of them; they are recorded
  # here so the paper's covariate screen is not lost.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all parameters and specifically flagged on the 1,25D3",
        "baseline concentration by covariate-plot inspection (Results 2.3);",
        "not retained. Median 92.0 kg, range 70.7-135.3 kg (Table 1)."
      ),
      source_name        = "Weight (kg)"
    ),
    BMI = list(
      description        = "Body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Flagged on the 1,25D3 volume of distribution by covariate-plot",
        "inspection (Results 2.3); not retained. Median 32.6 kg/m^2, range",
        "25.6-43.4 (Table 1); every subject had BMI > 25 kg/m^2, which the",
        "Discussion cites as a reason the obesity effect known from the",
        "literature could not be resolved here."
      ),
      source_name        = "BMI (kg/m2)"
    ),
    AGE = list(
      description        = "Age at enrolment.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (Methods 4.7); not retained. Median 61 years, range 29-73 (Table 1).",
      source_name        = "Age (years)"
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened as 'gender' (Methods 4.7); not retained. 17 of 29 (59",
        "percent) female (Table 1)."
      ),
      source_name        = "Gender"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White)",
      notes              = paste(
        "Screened (Methods 4.7); not retained. Table 1 reports only two race",
        "categories: White 19 of 29 (66 percent), Black 10 of 29 (34 percent)."
      ),
      source_name        = "Race"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic ethnicity indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Hispanic)",
      notes              = paste(
        "Screened as 'ethnicity' (Methods 4.7); not retained. 4 of 29 (14",
        "percent) Hispanic, 24 non-Hispanic, 1 declined to disclose",
        "(Table 1 footnote a)."
      ),
      source_name        = "Ethnicity"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Flagged on the 1,25D3 and 24,25D3 baseline concentrations by",
        "covariate-plot inspection (Results 2.3) and discussed as",
        "biologically plausible for 24,25D3, but not retained by the",
        "stepwise procedure. Median 37 mL/min/1.73 m^2, range 11-97",
        "(Table 1); CKD stages 1-5 represented."
      ),
      source_name        = "eGFR (mL/min/1.73m2)"
    ),
    PTH = list(
      description        = "Serum parathyroid hormone concentration.",
      units              = "not reported",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a serum protein covariate (Methods 4.7); not retained.",
        "No summary statistic is tabulated for it in Tuey 2024. No canonical",
        "entry exists in inst/references/covariate-columns.md because no",
        "model in the library has yet retained it; the name here is",
        "documentation of the paper's screen only, not a canonical claim."
      ),
      source_name        = "PTH"
    ),
    FGF23 = list(
      description        = "Serum fibroblast growth factor 23 concentration.",
      units              = "not reported",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a serum protein covariate (Methods 4.7); not retained.",
        "No summary statistic is tabulated for it in Tuey 2024. No canonical",
        "entry exists in inst/references/covariate-columns.md because no",
        "model in the library has yet retained it; the name here is",
        "documentation of the paper's screen only, not a canonical claim."
      ),
      source_name        = "FGF-23"
    ),
    SNP_CYP2R1_RS12794714 = list(
      description        = "CYP2R1 rs12794714 genotype (25-hydroxylase).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "G/G",
      notes              = paste(
        "Screened (Methods 4.7); not retained. G/G 15 (51 percent), G/A 13",
        "(45 percent), A/A 0, not determined 1 (Table 1)."
      ),
      source_name        = "CYP2R1 rs12794714"
    ),
    SNP_CYP27B1_RS10877012 = list(
      description        = "CYP27B1 rs10877012 genotype (1-alpha-hydroxylase).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "C/C",
      notes              = paste(
        "Screened (Methods 4.7); not retained. C/C 15 (51 percent), C/A 11",
        "(38 percent), A/A 2 (7 percent), not determined 1 (Table 1)."
      ),
      source_name        = "CYP27B1 rs10877012"
    ),
    SNP_CYP24A1_RS6013897 = list(
      description        = "CYP24A1 rs6013897 genotype (24-hydroxylase).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "A/A",
      notes              = paste(
        "Screened (Methods 4.7); not retained. A/A 18 (62 percent), A/S 8",
        "(28 percent), S/S 2 (7 percent), not determined 1 (Table 1)."
      ),
      source_name        = "CYP24A1 rs6013897"
    ),
    SNP_GC_RS7041 = list(
      description        = "GC rs7041 genotype (vitamin D binding protein).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "G/G",
      notes              = paste(
        "Screened (Methods 4.7); not retained. G/G 17 (59 percent), G/A 9",
        "(31 percent), A/A 2 (7 percent), not determined 1 (Table 1)."
      ),
      source_name        = "GC_VDBP rs7041"
    ),
    SNP_VDR_RS2228570 = list(
      description        = "VDR rs2228570 genotype (vitamin D receptor, FokI).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "G/G",
      notes              = paste(
        "Screened (Methods 4.7); not retained. G/G 17 (59 percent), A/G 6",
        "(21 percent), A/A 5, not determined 1 (Table 1)."
      ),
      source_name        = "VDR rs2228570"
    ),
    SNP_VDR_RS7968585 = list(
      description        = "VDR rs7968585 genotype (vitamin D receptor).",
      units              = "(genotype)",
      type               = "categorical",
      reference_category = "G/G",
      notes              = paste(
        "Screened (Methods 4.7); not retained. G/G 7 (24 percent), G/A 14",
        "(49 percent), A/A 7 (24 percent), not determined 1 (Table 1)."
      ),
      source_name        = "VDR rs7968585"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29,
    n_studies      = 1,
    age_range      = "29-73 years",
    age_median     = "61 years",
    weight_range   = "70.7-135.3 kg",
    weight_median  = "92.0 kg",
    sex_female_pct = 59,
    race_ethnicity = c(White = 66, Black = 34),
    disease_state  = paste(
      "Chronic kidney disease (stages 1-5, non-dialysis) with vitamin D",
      "insufficiency or deficiency, defined as total 25D3 below 30 ng/mL.",
      "Median baseline 25D3 18 ng/mL (range 7-29); median eGFR 37",
      "mL/min/1.73 m^2 (range 11-97); median BMI 32.6 kg/m^2 (range",
      "25.6-43.4). CKD stage distribution: 1 (3 percent), 5 (17 percent),",
      "14 (48 percent), 8 (28 percent), 1 (3 percent) for stages 1-5."
    ),
    dose_range     = "Single 5000 I.U. oral dose of cholecalciferol (equivalent to 325 nmol)",
    regions        = "United States (University of Colorado and University of Pittsburgh)",
    notes          = paste(
      "Baseline characteristics from Tuey 2024 Table 1; ClinicalTrials.gov",
      "NCT02360644. Serial plasma sampling at baseline and 0.5, 1, 2, 4, 8,",
      "12, 24, 48, 168 and 336 h (Methods 4.1). 310 plasma concentrations",
      "were used for model development, of which 212 parent VitD3",
      "observations (over 72 percent) were below the limit of",
      "quantification and handled with the Beal M3 method. Fitted in",
      "Phoenix NLME v8.3 with the QRPEM engine (Methods 4.5, Results 2.2)."
    )
  )

  ini({
    # Cholecalciferol (parent) -- Tuey 2024 Table 2. The CV% printed beside
    # each estimate in Table 2 is the relative standard error of the
    # estimate, not an IIV magnitude: Results 2.2 reads those same numbers
    # as precision ("estimated with adequate precision with the exception
    # of Vm2 and Vm3", whose Table 2 CV% are 206.8 and 140.5).
    # NOTE on units: every baseline below is taken from Table 2's nmol/L
    # column, which is the scale the model was fitted on (Methods 4.5
    # converted all observations to nmol/L). The ng/mL values given
    # parenthetically in the Results 2.2 prose do NOT reconcile with the
    # nmol/L values at the paper's own stated molecular weights and are
    # not used here -- see the vignette Errata section.
    lrbase <- log(0.98); label("Endogenous baseline cholecalciferol concentration (nmol/L)")  # Table 2, C0 = 0.98 nmol/L, RSE 41.7%
    lka    <- fixed(log(0.054)); label("Absorption rate constant (1/h)")  # Table 2, ka fixed to 0.054/h; Discussion confirms 0.054 (inadequate absorption-phase data)
    ksyn   <- fixed(0.55); label("Zero-order endogenous cholecalciferol production rate (nmol/h)")  # Table 2, kendog fixed to 0.55 nmol/h; Discussion "fixed endogenous rate of 0.55 nmol/h"
    lvc    <- log(21.3); label("Apparent central volume of distribution, Vc/F (L)")  # Table 2, Vc/FVitD3 = 21.3 L, RSE 22.2%
    lcl    <- log(1.4); label("Apparent clearance, CL/F (L/h)")  # Table 2, CL/FVitD3 = 1.4 L/h, RSE 42.4%
    lvp    <- fixed(log(50)); label("Apparent peripheral volume of distribution, Vp/F (L)")  # Table 2, Vp/FVitD3 fixed to 50 L
    lq     <- fixed(log(0.44)); label("Apparent intercompartmental clearance, Q/F (L/h)")  # Table 2, Q/FVitD3 fixed to 0.44 L/h; Discussion restates 0.44 L/h

    # 25-hydroxyvitamin D3 -- Tuey 2024 Table 2
    lrbase_25d3 <- log(43.5); label("Endogenous baseline 25D3 concentration (nmol/L)")  # Table 2, C0m1 = 43.5 nmol/L, RSE 4.1% (= 17.4 ng/mL, consistent with the Table 1 median of 18 ng/mL)
    lvc_25d3    <- log(58.3); label("25D3 volume of distribution, Vm1 (L)")  # Table 2, Vm1 = 58.3 L, RSE 14.8%
    lcl_25d3    <- log(0.02); label("25D3 clearance, CLm1 (L/h)")  # Table 2, CLm1 = 0.02 L/h, RSE 52.2%

    # 1,25-dihydroxyvitamin D3 -- Tuey 2024 Table 2
    lrbase_125d3 <- log(0.20); label("Endogenous baseline 1,25D3 concentration (nmol/L)")  # Table 2, C0m2 = 0.20 nmol/L, RSE 6.9%
    lvc_125d3    <- log(71.5); label("1,25D3 volume of distribution, Vm2 (L)")  # Table 2, Vm2 = 71.5 L, RSE 206.8%
    lcl_125d3    <- log(0.08); label("1,25D3 clearance, CLm2 (L/h)")  # Table 2, CLm2 = 0.08 L/h, RSE 47.7%

    # 24,25-dihydroxyvitamin D3 -- Tuey 2024 Table 2
    lrbase_2425d3 <- log(2.2); label("Endogenous baseline 24,25D3 concentration (nmol/L)")  # Table 2, C0m3 = 2.2 nmol/L, RSE 9.4%
    lvc_2425d3    <- log(105.2); label("24,25D3 volume of distribution, Vm3 (L)")  # Table 2, Vm3 = 105.2 L, RSE 140.5%
    lcl_2425d3    <- log(0.40); label("24,25D3 clearance, CLm3 (L/h)")  # Table 2, CLm3 = 0.40 L/h, RSE 53.4%

    # Formed fractions, both fixed for identifiability because no
    # metabolite was administered (Results 2.2, Methods 4.6). Both values
    # were taken from the Sawyer 2022 vitamin D PBPK model built on the
    # healthy-control arm of the same clinical study.
    fm_25d3  <- fixed(1); label("Fraction of cholecalciferol clearance forming 25D3 (unitless)")  # Results 2.2 / Methods 4.6, fm1 fixed to 1; no alternative VitD3 elimination pathway assumed
    fm_125d3 <- fixed(0.017); label("Fraction of 25D3 clearance forming 1,25D3 (unitless)")  # Results 2.2 / Methods 4.6, fm2 fixed to 0.017 from the Sawyer 2022 PBPK model; remainder 1 - fm2 forms 24,25D3

    # IIV. Results 2.2: "IIV terms for C0, Vc/FVitD3, and CL/FVitD3 were
    # included." Tuey 2024 publishes no OMEGA magnitude for any of the
    # three (Table 2's CV% column is the estimate's RSE, not an IIV), so
    # each is carried at fixed(0) to preserve the model's random-effect
    # structure without inventing a variance. Replace with the published
    # values if the authors release them.
    etalrbase ~ fixed(0)  # Results 2.2 declares IIV on C0; magnitude not published
    etalvc    ~ fixed(0)  # Results 2.2 declares IIV on Vc/FVitD3; magnitude not published
    etalcl    ~ fixed(0)  # Results 2.2 declares IIV on CL/FVitD3; magnitude not published

    # Residual error -- proportional on every analyte (Results 2.2,
    # Table 2). Table 2 prints each sigma as a percentage.
    propSd         <- 0.125; label("Proportional residual error, cholecalciferol (fraction)")  # Table 2, sigma1 = 12.5%
    propSd_25d3    <- 0.657; label("Proportional residual error, 25D3 (fraction)")  # Table 2, sigma2 = 65.7%
    propSd_125d3   <- 0.172; label("Proportional residual error, 1,25D3 (fraction)")  # Table 2, sigma3 = 17.2%
    propSd_2425d3  <- 0.166; label("Proportional residual error, 24,25D3 (fraction)")  # Table 2, sigma4 = 16.6%
  })

  model({
    # 1. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    cl_25d3   <- exp(lcl_25d3)
    vc_25d3   <- exp(lvc_25d3)
    cl_125d3  <- exp(lcl_125d3)
    vc_125d3  <- exp(lvc_125d3)
    cl_2425d3 <- exp(lcl_2425d3)
    vc_2425d3 <- exp(lvc_2425d3)

    rbase         <- exp(lrbase + etalrbase)
    rbase_25d3    <- exp(lrbase_25d3)
    rbase_125d3   <- exp(lrbase_125d3)
    rbase_2425d3  <- exp(lrbase_2425d3)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Endogenous pre-dose baselines as compartment initial conditions.
    # Each state is an amount (nmol), so the baseline concentration is
    # multiplied by that species' volume. peripheral1 is started at
    # distribution equilibrium with central (equal concentrations, since
    # k12 * central = k21 * peripheral1 implies central/vc = peripheral1/vp):
    # the subject has been at their endogenous baseline for life, so a
    # zero peripheral store would create a spurious distributional
    # transient at t = 0. Tuey 2024 reports only the central baseline C0.
    central(0)         <- rbase * vc
    peripheral1(0)     <- rbase * vp
    central_25d3(0)    <- rbase_25d3 * vc_25d3
    central_125d3(0)   <- rbase_125d3 * vc_125d3
    central_2425d3(0)  <- rbase_2425d3 * vc_2425d3

    # 4. Observations (concentrations drive the sequential formation terms)
    Cc         <- central / vc
    Cc_25d3    <- central_25d3 / vc_25d3
    Cc_125d3   <- central_125d3 / vc_125d3
    Cc_2425d3  <- central_2425d3 / vc_2425d3

    # 5. ODE system (Tuey 2024 Figure 1)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot + ksyn - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # All cholecalciferol clearance forms 25D3 (fm_25d3 = 1); 25D3
    # clearance splits into 1,25D3 (fm_125d3) and 24,25D3 (1 - fm_125d3).
    d/dt(central_25d3)    <- fm_25d3 * cl * Cc - cl_25d3 * Cc_25d3
    d/dt(central_125d3)   <- fm_125d3 * cl_25d3 * Cc_25d3 - cl_125d3 * Cc_125d3
    d/dt(central_2425d3)  <- (1 - fm_125d3) * cl_25d3 * Cc_25d3 - cl_2425d3 * Cc_2425d3

    # 6. Residual error
    Cc        ~ prop(propSd)
    Cc_25d3   ~ prop(propSd_25d3)
    Cc_125d3  ~ prop(propSd_125d3)
    Cc_2425d3 ~ prop(propSd_2425d3)
  })
}
