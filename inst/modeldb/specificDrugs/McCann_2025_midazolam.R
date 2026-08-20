McCann_2025_midazolam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for midazolam in 93",
    "children and young adults (2.1-20.1 years), 47.3% of whom had obesity,",
    "who received midazolam per standard of care as IV bolus and/or IV",
    "infusion in the Pediatric Trials Network real-world opportunistic",
    "sampling study (McCann 2025). All disposition parameters are",
    "allometrically scaled by total body weight standardised to 70 kg, with",
    "the exponents fixed a priori at 0.75 for the clearances (CL, Q) and 1",
    "for the volumes (V1, V2). Interindividual variability was supported",
    "only on clearance and is very large (185.4% CV); IIV on the remaining",
    "disposition parameters shrank to 100% and was removed. Residual",
    "variability is proportional (50.4%). Stepwise covariate modeling",
    "retained no covariate beyond total body weight - obesity status,",
    "BMI, BMI percentile for age and sex, extended BMI, postnatal and",
    "postmenstrual age (sigmoidal Emax maturation), hepatic and renal",
    "laboratory values, race, ethnicity, sex, and CYP3A4 inducer /",
    "inhibitor use were all screened and rejected. The paper's second",
    "analysis, a PK-Sim whole-body PBPK evaluation, applies the",
    "third-party Open Systems Pharmacology midazolam library model",
    "unchanged and publishes none of its equations or parameters, so it",
    "is not represented here (see vignette Assumptions and deviations)."
  )
  reference <- paste(
    "McCann S, Helfer VE, Balevic SJ, Muller WJ, van den Anker JN,",
    "Al-Uzri A, Meyer ML, Anderson SG, Turdalieva S, Edginton AN,",
    "Gonzalez D; On Behalf of the Best Pharmaceuticals for Children Act",
    "Pediatric Trials Network Steering Committee.",
    "Physiologically Based and Population Pharmacokinetic Modeling of",
    "Midazolam in Children With Obesity Using Real-World Data.",
    "Clin Transl Sci. 2025;18(5):e70247.",
    "doi:10.1111/cts.70247.",
    sep = " "
  )
  vignette <- "McCann_2025_midazolam"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states hold midazolam itself; the study
  # quantified midazolam plasma concentrations only (no metabolite was
  # modelled), and doses are recorded in mg of midazolam.
  compartmentData <- list(
    central     = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Allometric",
        "scaling per Supplementary Equation S2,",
        "P_i = theta_pop * (W_i / W_standard)^beta, with W_standard =",
        "70 kg for total body weight and beta fixed a priori at 0.75 for",
        "the clearances and 1 for the volumes (Methods 'Population PK",
        "Analysis'; Supplementary Table S4 writes the base model out",
        "explicitly as CL = theta_CL * (WT/70)^0.75, V1 = theta_V1 *",
        "(WT/70), Q = theta_Q * (WT/70)^0.75, V2 = theta_V2 * (WT/70)).",
        "Total body weight beat lean body weight and fat-free mass on OFV",
        "(from the unscaled base model at OFV 1834.599: delta-OFV -47.827",
        "for WT, -41.949 for LBW, -39.901 for FFM; Supplementary Table",
        "S4) and gave the largest reduction in clearance IIV (185.4% vs",
        "198.8% for LBW and 202.6% for FFM), so WT was carried forward.",
        "Note that Results 3.2 of the main text transposes the LBW and",
        "FFM delta-OFVs relative to Supplementary Table S4; the",
        "supplement is self-consistent (1834.599 - 1792.650 = 41.949 for",
        "LBW, 1834.599 - 1794.698 = 39.901 for FFM) and agrees with the",
        "IIV ordering, so the supplement values are used here.",
        "Estimating the exponents rather than fixing them did not improve",
        "the fit significantly and produced implausible values (1.26 for",
        "the clearances, 0.95 for the volumes), so the fixed exponents",
        "were retained (Discussion). Cohort weight range 7.2-193 kg,",
        "median 48.1 kg (Table 1). Time-varying weight is not indicated",
        "by the paper; weight is recorded per participant."
      ),
      source_name        = "WT"
    )
  )

  # Covariates screened by stepwise covariate modeling (Methods 'Population
  # PK Analysis'; full run-by-run results in Supplementary Table S4) but NOT
  # retained in the final model: "None of the evaluated covariates were
  # significant" (Results 3.2). Documented here so the paper's covariate
  # screen is preserved without declaring covariates that model() never
  # references.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 23.4 kg/m^2 (range 12.6-96.9; Table 1)."
    ),
    BMI_PCTL = list(
      description = "BMI percentile for age and sex, and the derived extended BMI (participant BMI divided by the 95th-percentile BMI for age and sex, expressed as a percentage; extended BMI >= 100% denotes obesity)",
      units       = "percentile / percent",
      type        = "continuous",
      notes       = "Screened on CL as both the raw percentile and the extended BMI; not retained. Non-canonical key - this concept has no entry in inst/references/covariate-columns.md and is documentation-only here, so no new canonical name is being proposed."
    ),
    OBESE = list(
      description = "Obesity status, defined as BMI >= the 95th percentile for age and sex",
      units       = "binary (1 = with obesity, 0 = without obesity or missing)",
      type        = "categorical",
      notes       = "Screened on CL both alone and jointly with a severe-obesity indicator (BMI >= 99th percentile); neither was retained (Supplementary Table S4 footnotes b and c). 44/93 participants (47.3%) had obesity. Non-canonical key; documentation-only."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "years in this paper (canonical PNA is months)",
      type        = "continuous",
      notes       = "Screened on CL through the sigmoidal Emax maturation function of Supplementary Equation S8, both with the Hill coefficient estimated (delta-OFV -1.273) and fixed to 1 (delta-OFV +11.122); neither reached the delta-OFV -3.84 forward-inclusion threshold."
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units       = "weeks in this paper (canonical PAGE is months)",
      type        = "continuous",
      notes       = "Screened on CL through the same sigmoidal Emax maturation function (delta-OFV -1.159 estimated Hill, +89.416 with Hill fixed to 1); not retained. Cohort median 740 weeks (range 149-1090; Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "reported in g/dL by the source (canonical ALB is g/L)",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Values in Supplementary Table S3."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "reported in mg/dL by the source (canonical CREAT accepts umol/L or mg/dL)",
      type        = "continuous",
      notes       = "Screened on CL; not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "reported in mg/dL by the source (canonical TBILI is umol/L)",
      type        = "continuous",
      notes       = "Screened on CL together with direct bilirubin; neither was retained."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "reported in mg/dL by the source (canonical DBIL is umol/L)",
      type        = "continuous",
      notes       = "Screened on CL; not retained."
    ),
    SEXF = list(
      description = "Sex",
      units       = "binary (1 = female)",
      type        = "categorical",
      notes       = "Screened on CL; not retained. 45/93 participants (48.4%) were female (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL in several groupings (a multiplicative shift for all non-White categories relative to White/missing, delta-OFV -12.132 on 4 df; a White vs non-White dichotomy, -0.6 on 1 df; and a three-level Black / other-non-White / White grouping, -8.876 on 2 df - Supplementary Table S4 footnotes e, f, g). The first and third pass forward inclusion but neither survives backward elimination, whose p < 0.01 thresholds are 13.28 on 4 df and 9.21 on 2 df (Methods 'Population PK Analysis'), so none was retained - 'None of the evaluated covariates were significant' (Results 3.2). 19/93 participants (20.4%) were Black or African American (Table 1)."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL as ETHN (Supplementary Table S4 footnote d); not retained. 15/93 participants (16.1%) were Hispanic or Latino (Table 1)."
    ),
    CYP3A4_INDUCER = list(
      description = "Concomitant use of a drug with potential to induce CYP3A4 (phenobarbital, phenytoin, dexamethasone)",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL; not retained (Supplementary Table S4 footnote h). Non-canonical key; documentation-only."
    ),
    CYP3A4_INHIBITOR = list(
      description = "Concomitant use of a drug with potential to inhibit CYP3A4 (amiodarone, amlodipine, azithromycin, dexmedetomidine, diazepam, hydromorphone, fentanyl and others per Supplementary Table S4 footnote i)",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL both as the pooled inhibitor flag and, separately, as dexmedetomidine use alone (delta-OFV -0.338); neither was retained. Non-canonical key; documentation-only."
    ),
    VASOPRESSOR = list(
      description = "Concomitant vasopressor use",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL (delta-OFV -0.002; Supplementary Table S4); not retained. Non-canonical key; documentation-only."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 93L,
    n_studies        = 1L,
    n_concentrations = 164L,
    age_range        = "2.1-20.1 years postnatal age",
    age_median       = "13.4 years postnatal age",
    weight_range     = "7.2-193 kg",
    weight_median    = "48.1 kg",
    sex_female_pct   = 48.4,
    race_ethnicity   = c(
      White                                  = 71.0,
      `Black or African American`            = 20.4,
      `Non-White / non-Black`                = 5.4,
      `Unknown or not reported`              = 3.2,
      `Hispanic or Latino (ethnicity)`       = 16.1
    ),
    disease_state    = paste(
      "Hospitalised children and young adults (< 21 years at enrolment)",
      "receiving intravenous midazolam per standard of care for sedation,",
      "enrolled without a protocol-specified indication. 44 of 93",
      "participants (47.3%) had obesity, defined as a BMI at or above the",
      "95th percentile for age and sex. Cohort BMI median 23.4 kg/m^2",
      "(range 12.6-96.9); one participant had an implausibly high recorded",
      "BMI of 96.9 kg/m^2 (Table 1 footnote c)."
    ),
    dose_range       = paste(
      "Standard-of-care intravenous midazolam, site specific. Bolus doses",
      "median 0.07 mg/kg (range 0.0007-0.21) in participants without",
      "obesity and 0.05 mg/kg (0.001-0.20) in those with obesity;",
      "infusion rates median 0.1 mg/h/kg in both groups (ranges 0.002-3.0",
      "and 0.002-2.0 mg/h/kg respectively). 38 participants received",
      "bolus doses only, 11 infusions only, and 44 both",
      "(Supplementary Table S2)."
    ),
    regions          = "United States, Canada, and the United Kingdom (18 sites; Supplementary Table S1)",
    obesity_status   = "44 with obesity (47.3%), 46 without obesity, 3 unknown / unevaluable",
    notes            = paste(
      "Best Pharmaceuticals for Children Act Pediatric Trials Network",
      "opportunistic real-world study, ClinicalTrials.gov NCT01431326.",
      "96 participants were enrolled and 169 PK samples collected; three",
      "participants were excluded (one who received an oral dose, two",
      "missing key dosing information), and two further samples were",
      "excluded (one an unusually high concentration 163.8 h after the",
      "dose, the other higher than the preceding sample with no dose in",
      "between), leaving 93 participants with 164 plasma concentrations,",
      "all above the lower limit of quantification (Supplementary Results",
      "'Participant Information'). Observed concentrations spanned",
      "0.13-17825.89 ng/mL, most drawn within the first hour after a dose.",
      "Samples per participant median 2 (range 1-4). Fitted in NONMEM 7.5",
      "with FOCE-I; run management in Pirana 2.9.9, bootstrap (1000",
      "replicates, 99.3% successful) in PsN 4.9.0. Shrinkage was < 6% for",
      "the clearance random effect and 15% for the residual error; all",
      "final estimates had RSE < 40%. Demographics from Table 1;",
      "laboratory values in Supplementary Table S3; empirical Bayes",
      "estimates of clearance by age and obesity subgroup in",
      "Supplementary Table S5."
    )
  )

  ini({
    # ================================================================
    # Structural disposition - McCann 2025 Table 2 'Final population PK
    # parameter estimates'. All four values are reported standardised to
    # a 70 kg individual (the table's own units carry the /70 kg), so
    # these are the typical values at WT = 70 kg.
    # ================================================================
    lcl <- log(14.9)
    label("Typical clearance CL at WT = 70 kg (L/h)")                    # Table 2 structural model: CL = 14.9 L/h/70 kg (RSE 14%; bootstrap mean 14.68, 95% CI 10.43-19.59)
    lvc <- log(6.09)
    label("Typical central volume V1 at WT = 70 kg (L)")                 # Table 2 structural model: V1 = 6.09 L/70 kg (RSE 37%; bootstrap mean 7.09, 95% CI 0.44-22.79)
    lq  <- log(45.7)
    label("Typical intercompartmental clearance Q at WT = 70 kg (L/h)")  # Table 2 structural model: Q = 45.7 L/h/70 kg (RSE 23%; bootstrap mean 44.71, 95% CI 7.14-81.43)
    lvp <- log(34.5)
    label("Typical peripheral volume V2 at WT = 70 kg (L)")              # Table 2 structural model: V2 = 34.5 L/70 kg (RSE 17%; bootstrap mean 34.07, 95% CI 15.76-57.89)

    # ================================================================
    # Allometric exponents, fixed a priori (Methods 'Population PK
    # Analysis': "The exponent beta was fixed to 0.75 for clearance
    # parameters and 1 for distribution volumes"). Supplementary
    # Equation S2 gives the scaling form and W_standard = 70 kg for
    # total body weight. Estimating them instead was rejected
    # (Discussion): the estimates (1.26 for clearances, 0.95 for
    # volumes) did not improve the fit significantly and disagreed with
    # reported values, so the fixed exponents were retained.
    # ================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT / 70) for CL and Q (unitless)")    # Methods 'Population PK Analysis'; Supplementary Equation S2 and Table S4
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on (WT / 70) for V1 and V2 (unitless)")   # Methods 'Population PK Analysis'; Supplementary Equation S2 and Table S4

    # ================================================================
    # Interindividual variability. Exponential IIV (Supplementary
    # Equation S1: P_ij = theta_pop,j * exp(eta_ij)) on clearance only -
    # "IIV was added only on CL once estimates of variability in other
    # parameters resulted in shrinkage of 100%" (Results 3.2). Table 2
    # reports IIV as %CV; the table's own footnote gives the conversion
    # CV(%) = sqrt(exp(omega^2) - 1) * 100, so
    # omega^2 = log(1 + CV^2) = log(1 + 1.854^2) = 1.49005.
    # ================================================================
    etalcl ~ 1.49005
    # Table 2 interindividual variability: CL = 185.4% CV (RSE 9%, shrinkage 6%; bootstrap mean 187.8%, 95% CI 127.0-291.3%)

    # ================================================================
    # Residual variability - proportional only ("A two-compartment PK
    # model with proportional residual error described the midazolam
    # concentration data well", Results 3.2). Table 2 reports it as a
    # percentage, i.e. a 0.504 fractional SD.
    # ================================================================
    propSd <- 0.504
    label("Proportional residual error (fraction)")                      # Table 2 residual error: proportional error = 50.4% (RSE 9%, shrinkage 15%; bootstrap mean 48.5%, 95% CI 36.1-59.1%)
  })

  model({
    # ----------------------------------------------------------------
    # Individual parameters. Allometric scaling on total body weight
    # standardised to 70 kg, exponents fixed (Supplementary Equation S2
    # and the base model written out in Supplementary Table S4):
    #   CL = theta_CL * (WT/70)^0.75
    #   V1 = theta_V1 * (WT/70)
    #   Q  = theta_Q  * (WT/70)^0.75
    #   V2 = theta_V2 * (WT/70)
    # IIV is carried on CL only.
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_cl
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc

    # ----------------------------------------------------------------
    # Micro-constants.
    # ----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----------------------------------------------------------------
    # ODE system. Midazolam was given intravenously only (bolus and/or
    # infusion); dose enters the central compartment directly.
    # ----------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------
    # Observation. Doses are in mg and volumes in L, so central / vc is
    # mg/L; the study reports plasma midazolam in ng/mL (concentrations
    # 0.13-17825.89 ng/mL, Results 3.1), and 1 mg/L = 1000 ng/mL.
    # ----------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
