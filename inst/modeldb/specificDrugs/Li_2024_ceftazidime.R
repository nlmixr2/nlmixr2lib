Li_2024_ceftazidime <- function() {
  description <- "One-compartment IV population PK model for ceftazidime in critically ill children (0.03-15 years) admitted to a paediatric intensive care unit (Li 2024). Clearance scales as a power function of body weight (exponent 0.90, reference 70 kg) and of modified-Schwartz estimated GFR (exponent 0.38, reference 116.93 mL/min/1.73 m^2); the volume of distribution scales linearly with body weight (exponent fixed to 1, reference 70 kg). Residual variability is additive."
  reference <- "Li M, Gao L, Wang Z, Zeng L, Chen C, Wang J, Li S, Liu M, Wang Y. Population pharmacokinetics and dose optimization of ceftazidime in critically ill children. Front Pharmacol. 2024;15:1470350. doi:10.3389/fphar.2024.1470350"
  vignette <- "Li_2024_ceftazidime"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Ceftazidime was given intravenously (Li 2024 Sect. 2.2),
  # so the dose enters `central` directly and there is no depot state.
  # `central` is verified: Li 2024 Sect. 2.3 states that ceftazidime was
  # extracted from SERUM and assayed by HPLC-UV at 230 nm, and Table 1 reports
  # the observations as "Ceftazidime serum concentration".
  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Li 2024 Table 1 mean 21.07 kg (SD 14.70), median 18.35, range 2.80-95. Enters BOTH structural parameters as a power of (WT/70): exponent 0.90 on CL and exponent 1 (fixed) on Vd, per Li 2024 Eqs. 5 and 6. The 70 kg reference is stated in the Results text under Table 3 ('normalized by a weight of 70 kg') and is a nominal adult reference, NOT the cohort median -- every subject in this paediatric cohort except one sits below it, so the two structural estimates (27.83 L and 7.76 L/h) are extrapolations to a 70 kg subject rather than values observed in the data. A consequence of the exponent pair worth carrying forward: because Vd scales with WT^1 while CL scales with WT^0.90, the elimination rate constant kel = CL/Vd scales with WT^(-0.10), so SMALLER children eliminate ceftazidime faster per unit volume. Li 2024 Sect. 3.4 reports exactly that pattern from its own simulations ('ceftazidime concentrations were lower for those weighing less than 10 kg').",
      source_name        = "weight"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the modified Schwartz equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = 0.413 * height (cm) / serum creatinine (mg/dL)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Li 2024 Sect. 2.2 gives the modified Schwartz formula (citing Schwartz 2009); Table 1 reports mean 116.41 (SD 31.35), median 116.93, range 40.42-197.15. Reference value 116.93 is the cohort MEDIAN and is the normalizing constant in Li 2024 Eq. 5; the Results text under Table 3 states the estimates were 'normalized by ... a median eGFR of 116.93 mL/min/1.73 m^2'. Power effect (CRCL/116.93)^0.38 on CL. Note the unit trap in the source: Table 1 reports serum creatinine in umol/L (median 33.9) whereas the Schwartz formula consumes mg/dL; the median subject reproduces the reported median eGFR only under the mg/dL reading (0.413 * 111 cm / (33.9 / 88.4) mg/dL = 119.6, against the reported median of 116.93), which settles the formula's input unit. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates. Renal strata (Li 2024 Sect. 3.1, n = 88): moderate insufficiency 30-60 in 5 children, mild 60-90 in 9, normal 90-120 in 31, augmented 120-200 in 43 -- so nearly half the cohort had augmented renal clearance and the model carries far more information about supranormal than about impaired filtration.",
      source_name        = "eGFR"
    )
  )

  # Screened in the Li 2024 covariate analysis (Sect. 2.5 lists the full
  # candidate set) but NOT retained in the final model, so these are
  # documentation only and are not referenced in model(). Li 2024 Sect. 3.2
  # states that age, weight, height, BSA, BMI, Cys-C and eGFR passed the
  # preliminary hypothesis test (p < 0.05) and entered forward selection, of
  # which only weight and eGFR survived; the remaining candidates below were
  # screened but did not reach the forward-selection step. The per-covariate
  # hypothesis-test statistics are in Supplementary Table S1, which is not
  # part of the PMC open-access deposit -- see the vignette Errata.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 5.43 years (SD 4.10), median 5.17, range 0.03-15. Reached forward selection (Sect. 3.2) but was not retained. Li 2024 Table 2 steps 5-9 additionally tested an age-dependent exponent maturation model (Model V, OFV 551.16) against the retained simple-exponent model (Model II, OFV 551.56); Model V was rejected because 'model structure was unstable' despite the marginally lower OFV.",
      source_name = "age"
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 108.33 cm (SD 32.45), median 111, range 50-172. Reached forward selection (Sect. 3.2) but was not retained as a standalone covariate. Height is nonetheless an INPUT to the retained CRCL covariate through the modified Schwartz equation.",
      source_name = "height"
    ),
    BSA = list(
      description = "Body surface area, Mosteller formula: BSA (m^2) = sqrt(height (cm) * weight (kg) / 3600)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 0.78 m^2 (SD 0.37), median 0.76, range 0.20-2.09; the Mosteller formula is given in the Table 1 preamble. Reached forward selection (Sect. 3.2) but was not retained; collinear with body weight, which was.",
      source_name = "BSA"
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Named among the covariates that entered forward selection in Li 2024 Sect. 3.2 but not retained. BMI is not tabulated in Li 2024 Table 1, so no cohort summary is available.",
      source_name = "BMI"
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Listed as a candidate renal-function marker in Li 2024 Sect. 2.5 and named among the covariates entering forward selection in Sect. 3.2, but not retained -- the modified-Schwartz creatinine-based eGFR won the renal-function slot. Cystatin C is not tabulated in Li 2024 Table 1, so no cohort summary or unit is stated in the source; the unit above is the conventional one and is not a paper-sourced value.",
      source_name = "Cys-C"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      reference_category = "male",
      notes       = "Li 2024 Table 1 reports 56 male / 32 female (36.4% female). Screened per Sect. 2.5; did not reach forward selection.",
      source_name = "gender"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 34.99 umol/L (SD 12.10), median 33.9, range 15.9-66.2. Screened per Sect. 2.5; did not reach forward selection as a standalone covariate. It is nonetheless an INPUT to the retained CRCL covariate through the modified Schwartz equation, which consumes it in mg/dL rather than the umol/L of Table 1.",
      source_name = "SCR"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 4.69 mmol/L (SD 2.51), median 4.19, range 1-14.1. Screened per Sect. 2.5; did not reach forward selection.",
      source_name = "BUN"
    ),
    UA = list(
      description = "Serum uric acid",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 317.11 umol/L (SD 160.70), median 265.3, range 38-789. Screened per Sect. 2.5; did not reach forward selection. Recorded under the source paper's own column name because serum uric acid has no entry in inst/references/covariate-columns.md; no canonical is proposed here since no shipped model retains it as a covariate (same treatment as Nie_2023_nalbuphine.R).",
      source_name = "UA"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 62.45 U/L (SD 220.96), median 20.5, range 4-1672. Screened as a hepatic-function candidate per Sect. 2.5; did not reach forward selection. The mean sits far above the median because of a small number of extreme values, which is the expected shape in a PICU cohort.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 141.11 U/L (SD 785.87), median 29, range 9-7758. Screened as a hepatic-function candidate per Sect. 2.5; did not reach forward selection.",
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Li 2024 Table 1 mean 22.04 umol/L (SD 52.55), median 10.1, range 2.9-408.7. Screened as a hepatic-function candidate per Sect. 2.5; did not reach forward selection.",
      source_name = "TBIL"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 88L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 100L,
    age_range        = "0.03-15 years",
    age_median       = "5.17 years (mean 5.43, SD 4.10)",
    weight_range     = "2.80-95 kg",
    weight_median    = "18.35 kg (mean 21.07, SD 14.70)",
    sex_female_pct   = 36.4,
    disease_state    = "Critically ill children in a paediatric intensive care unit, diagnosed with or suspected of having a bacterial infection and treated with intravenous ceftazidime for at least 3 consecutive days. Excluded: enrolment in another clinical trial, incomplete dose information, intolerance to ceftazidime. No patient discontinued ceftazidime for an adverse event.",
    dose_range       = "Intravenous ceftazidime 25-100 mg/kg, tailored to clinical status; median loading dose 48.38 mg/kg (range 21.05-76)",
    regions          = "China (Wuhan Children's Hospital / Wuhan Maternal and Child Healthcare Hospital, Tongji Medical College, Huazhong University of Science and Technology)",
    renal_function   = "Modified-Schwartz eGFR mean 116.41 mL/min/1.73 m^2 (SD 31.35), median 116.93, range 40.42-197.15. Strata: moderate insufficiency (30-60) 5 children, mild insufficiency (60-90) 9, normal (90-120) 31, augmented (120-200) 43.",
    notes            = "Prospective, open-label population PK study run October 2019 to January 2020 (Li 2024 Sect. 2.1). Opportunistic sampling gave 1-3 samples per patient, 100 serum concentrations from 88 patients; observed concentrations ranged 0.046-74.84 ug/mL with a median of 3.35. Serum ceftazidime was assayed by validated HPLC-UV (230 nm) over a linear range of 0.025-100 ug/mL with intra- and inter-day precision below 10%. The model was fit in Phoenix NLME 8.2. Covariate selection used forward selection (dOFV < -3.84, p < 0.05) then backward elimination (dOFV < 6.635, p < 0.01); no covariate was removed in backward elimination. Five allometric / maturation parameterizations were then compared (Li 2024 Table 2 steps 5-9) and the simple-exponent model was retained. The final model was evaluated by a 1000-replicate nonparametric bootstrap and by NPDE (mean 0.026, variance 0.9788; global test p = 0.721). Registered as MR-42-22-000220; ethics approval 2021R153-E03. A two-compartment model reduced the OFV by 4.63% but was rejected because the opportunistic sampling design could not support its parameter estimates (Li 2024 Sect. 3.2)."
  )

  ini({
    # Structural parameters (Li 2024 Table 3). The reference subject weighs
    # 70 kg and has a modified-Schwartz eGFR of 116.93 mL/min/1.73 m^2; the
    # Results text under Table 3 states both normalizing constants explicitly.
    lcl <- log(7.76);  label("Clearance at WT = 70 kg and CRCL = 116.93 mL/min/1.73 m^2 (L/h)")  # Li 2024 Table 3: 7.76 L/h (RSE 10.40%, bootstrap median 7.79, 95% CI 6.14-10.03)
    lvc <- log(27.83); label("Volume of distribution at WT = 70 kg (L)")                          # Li 2024 Table 3: 27.83 L (RSE 7.85%, bootstrap median 27.96, 95% CI 21.04-37.70)

    # Covariate exponents (Li 2024 Table 3 theta1-theta3, whose footnote maps
    # each theta to its covariate/parameter pair). theta1 is reported as
    # "1 (fixed)" in both the final-model and bootstrap columns, so it is
    # wrapped in fixed(); theta2 and theta3 were estimated and are not.
    e_wt_vc   <- fixed(1); label("Power exponent on (WT/70) for Vc (unitless)")                 # Li 2024 Table 3 theta1: 1 (fixed); appears in Eq. 6 as an unwritten exponent of 1 on (weight/70)
    e_wt_cl   <- 0.90;     label("Power exponent on (WT/70) for CL (unitless)")                 # Li 2024 Table 3 theta2: 0.90 (RSE 6.41%); Eq. 5
    e_crcl_cl <- 0.38;     label("Power exponent on (CRCL/116.93) for CL (unitless)")           # Li 2024 Table 3 theta3: 0.38 (RSE 33.59%); Eq. 5

    # Inter-individual variability. Li 2024 Eq. 1 defines the exponential IIV
    # model Pi = theta * exp(eta_i) with eta_i ~ N(0, omega^2), and the Table 3
    # rows are labelled omega^2 (the footnote defines the unsquared omega as
    # the "square root of inter-individual variance"), so these values are
    # VARIANCES on the log scale and are used as printed -- no CV%-to-variance
    # back-transformation is required.
    etalcl ~ 0.06  # Li 2024 Table 3: omega^2_CL = 0.06 (RSE 21.40%, 95% CI 0.03-0.09); ~24.9% CV
    etalvc ~ 0.04  # Li 2024 Table 3: omega^2_Vd = 0.04 (RSE 39.77%, 95% CI 0.002-0.07); ~20.2% CV

    # Residual variability. Li 2024 Sect. 3.2 states that the additive model
    # (Eq. 2) outperformed the proportional (Eq. 3) and combined (Eq. 4)
    # models, so the final model carries an additive error only.
    addSd <- 1.20; label("Additive residual error (ug/mL)")  # Li 2024 Table 3: sigma = 1.20 mg/L (RSE 14.60%, bootstrap median 1.17, 95% CI 0.67-1.52)
  })
  model({
    # Li 2024 Eqs. 5 and 6:
    #   CL (L/h) = 7.76  * (weight/70)^0.90 * (eGFR/116.93)^0.38 * exp(eta_CL)
    #   Vd (L)   = 27.83 * (weight/70)^1                         * exp(eta_Vd)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 116.93)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L, so central/vc is mg/L == ug/mL, matching the
    # units of the Li 2024 Table 3 additive residual error (mg/L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
