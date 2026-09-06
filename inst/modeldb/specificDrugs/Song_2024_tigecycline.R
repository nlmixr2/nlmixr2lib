Song_2024_tigecycline <- function() {
  description <- "Two-compartment IV population PK model for tigecycline in critically ill adult patients undergoing continuous renal replacement therapy, with no covariates retained in the final model (Song 2024)"
  reference <- "Song S, Liu J, Su W, Yu H, Feng B, Wu Y, Guo F, Yu Z. Population pharmacokinetics of tigecycline for critically ill patients undergoing continuous renal replacement therapy. Drug Des Devel Ther. 2024;18:4459-4469. doi:10.2147/DDDT.S473080"
  vignette <- "Song_2024_tigecycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Song 2024 Methods (plasma tigecycline
  # concentration collected per the protocol of the companion study, ref 14)
  # and the two-compartment structure of Table 2.
  compartmentData <- list(
    central     = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # The final model contains NO covariates. Song 2024 Results: "no covariates
  # were found to adequately explain the viability [sic: variability] in the
  # pharmacokinetic parameters of tigecycline". AST reached forward-inclusion
  # significance on V2 but was removed; see covariatesDataExcluded below.
  covariateData <- list()

  # Covariates the paper collected and screened but did NOT retain. Documented
  # here (rather than in covariateData) so the provenance of the covariate
  # screen survives without triggering a "declared but not referenced"
  # convention warning. Canonical names per inst/references/covariate-columns.md.
  covariatesDataExcluded <- list(
    AST = list(
      description        = "Serum aspartate aminotransferase activity. The ONLY covariate that reached forward-inclusion significance: Song 2024 Results, 'the inclusion of AST in the peripheral volume distribution led to a decrease in the OFV of 13.719'. It was then removed from the final model because 'its inclusion had no impact on the PTA or CFR outcomes' (Results) and the change was judged 'statistically significant but not clinically significant' (Discussion).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column AST, Table 1 (observed range 8-106 U/L). The effect form and coefficient are not printed anywhere in the paper - only the delta OFV is given - so the AST-on-V2 relationship cannot be reconstructed even if a user wanted it. Screening detail is in Table S1, which is not distributed through PubMed Central."
    ),
    WT = list(
      description        = "Body weight. Collected per Methods ('The patients' physio-pathological data, such as age, sex, weight (WT), body mass index (BMI), ...') and screened in the stepwise covariate search; not retained.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column WT, Table 1 (range 42-100 kg, median 70 kg). Note this differs from the companion non-CRRT model Su_2024_tigecycline, where WT was retained as a power term on both volumes."
    ),
    AGE = list(
      description        = "Age. Collected per Methods and screened; not retained.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column AGE, Table 1 (range 22-87 years, median 71 years)."
    ),
    SEXF = list(
      description        = "Sex, 1 = female. Collected per Methods and screened; not retained.",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "male (SEXF = 0)",
      notes              = "Source column SEX, coded M/F in Table 1; SEXF = 1 for the 6 of 21 patients recorded as F. The canonical SEXF orientation (1 = female) matches the source directly, so no value inversion is needed."
    ),
    BMI = list(
      description        = "Body mass index. Collected per Methods and screened; not retained.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BMI, Table 1 (range 16.4-37.5 kg/m^2)."
    ),
    BUN = list(
      description        = "Blood urea nitrogen. Reported per patient in Table 1 and screened; not retained.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BUN, Table 1 in SI units mmol/L (range 3.64-37.3)."
    ),
    CREAT = list(
      description        = "Serum creatinine. Reported per patient in Table 1 and screened; not retained.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column SCr, Table 1 in SI units umol/L (range 47-469). In a CRRT-dependent cohort serum creatinine reflects the dialysis prescription more than native renal function, which is consistent with its failure to explain any tigecycline PK variability here."
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity. Reported per patient in Table 1 and screened; not retained.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column ALT, Table 1 (range 3-150 U/L)."
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase activity. Reported per patient in Table 1 and screened; not retained.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column ALP, Table 1 (range 59-184 U/L)."
    ),
    TPRO = list(
      description        = "Total serum protein. Reported per patient in Table 1 and screened; not retained.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column TP, Table 1 in SI units g/L (range 42.7-66.6)."
    ),
    GGT = list(
      description        = "Serum gamma-glutamyl transferase activity. Reported per patient in Table 1 and screened; not retained.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column GGT, Table 1 (range 16-181 U/L). Retained on Q in the companion non-CRRT model Su_2024_tigecycline but not here."
    ),
    TBILI = list(
      description        = "Total serum bilirubin. Reported per patient in Table 1 and screened; not retained.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column TBIL, Table 1 in SI units umol/L (range 11.4-113.9). Bilirubin was the covariate Broeker et al retained on tigecycline CL in their CRRT cohort (Song 2024 Discussion, ref 12); Song 2024 explicitly did not reproduce that finding."
    ),
    ALB = list(
      description        = "Serum albumin. Reported per patient in Table 1 and screened; not retained.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column ALB, Table 1 in SI units g/L (range 23.0-35.7). Retained on V2 in the companion non-CRRT model Su_2024_tigecycline but not here."
    ),
    URINE_VOL_24H = list(
      description        = "24-hour urine volume, the residual-native-renal-function marker in a CRRT cohort. Named explicitly in Methods as a collected covariate; not retained.",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column UV, Table 1 (range 0-6950 mL; one anuric patient, several oliguric). Screened but not retained, consistent with the paper's overall negative covariate result."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = 167L,
    age_range      = "22-87 years",
    age_median     = "71 years",
    weight_range   = "42-100 kg",
    weight_median  = "70 kg",
    bmi_range      = "16.4-37.5 kg/m^2",
    sex_female_pct = 28.6,
    race_ethnicity = "Not reported (single-center Chinese ICU cohort)",
    disease_state  = "Critically ill adults receiving intermittent intravenous tigecycline while on continuous renal replacement therapy. Patients who did not receive CRRT continuously throughout tigecycline treatment were excluded.",
    dose_range     = "Intermittent intravenous tigecycline every 12 h; the dose was chosen by the treating physician. Dose regimens evaluated by Monte Carlo simulation were 50, 100 and 150 mg every 12 h.",
    regions        = "China (Sir Run Run Shaw Hospital, School of Medicine, Zhejiang University, Hangzhou)",
    renal_function = "All patients CRRT-dependent. Serum creatinine 47-469 umol/L, BUN 3.64-37.3 mmol/L, 24-hour urine volume 0-6950 mL (Table 1).",
    hepatic_function = "Albumin 23.0-35.7 g/L; total bilirubin 11.4-113.9 umol/L; AST 8-106 U/L; ALT 3-150 U/L; GGT 16-181 U/L; ALP 59-184 U/L; total protein 42.7-66.6 g/L (Table 1).",
    crrt           = "20 of 21 patients on CVVH and 1 on CVVHDF; CRRT intensity 21.0-51.2 mL/kg/h; anticoagulation heparin in 12, regional sodium citrate in 5, none in 4 (Table 1). CRRT modality, intensity and anticoagulant were all screened as covariates and none was significant; the authors attribute this partly to only one CVVHDF patient being available.",
    notes          = "Baseline demographics per Song 2024 Table 1, which lists all 21 patients individually. Single-center prospective study with intensive sampling, July 2019 - July 2023. The dosing and sampling protocol is not restated in this paper; it is carried from the companion study (ref 14, Su et al Front Pharmacol 2024;15:1342947), where tigecycline was given as a 30-minute intravenous infusion every 12 h and samples were drawn immediately before the seventh dose and at 0.5, 1, 2, 3, 4, 6 and 12 h post-dose. Model built in NONMEM 7.5.0 with PDx-Pop 5.3.1; two-compartment structure selected over one-compartment on AIC (2188.14 vs 2425.331). Stepwise covariate search used forward inclusion at delta OFV > 3.84 and backward elimination at delta OFV > 10.83; no covariate survived into the final model."
  )

  ini({
    # Structural parameters, Song 2024 Table 2, "CRRT / Final Model /
    # Estimate [RSE (%)]" column. The bootstrap medians in the adjacent column
    # (4.27, 31.6, 35.0, 98.5) agree with every estimate to within 2%.
    #
    # NOTE on lcl: Table 2, the Discussion ("4.22 L/h, 34.8 L/h, 30.9 L and
    # 98.7 L") and the bootstrap 95% CI (3.47-5.07) all give CL = 4.22 L/h.
    # The Abstract prints 4.42 L/h for the same quantity while agreeing with
    # Table 2 on the other three parameters. The Abstract is a transcription
    # slip; 4.22 governs. See the vignette "Assumptions and deviations".
    lcl <- log(4.22); label("Clearance (L/h)")                       # Song 2024 Table 2: CL = 4.22 L/h [RSE 10.4%]
    lvc <- log(30.9); label("Central volume of distribution (L)")    # Song 2024 Table 2: V1 = 30.9 L [RSE 15.9%]
    lq  <- log(34.8); label("Intercompartmental clearance (L/h)")    # Song 2024 Table 2: Q  = 34.8 L/h [RSE 8.94%]
    lvp <- log(98.7); label("Peripheral volume of distribution (L)") # Song 2024 Table 2: V2 = 98.7 L [RSE 11.8%]

    # Inter-individual variability, Song 2024 Table 2 "Inter-individual
    # variability" block, printed as "omega CL (%) 22.4" and
    # "omega V1 (%) 55.2". Exponential (log-normal) IIV per Methods. Only CL
    # and V1 carry an eta; Table 2 lists no omega for Q or V2.
    #
    # SCALE: these are the OMEGA VARIANCES expressed in percent, i.e.
    # omega^2_CL = 0.224 and omega^2_V1 = 0.552 (log-scale SDs 0.473 and
    # 0.743), NOT coefficients of variation. The paper's own Monte Carlo
    # target-attainment results settle this with no free parameters: the six
    # PTA values printed in Results, inverted through the standard normal,
    # imply a median steady-state AUC0-24 of 23.1 mg*h/L under the variance
    # reading (against 100 mg/day / 4.22 L/h = 23.7, a 2.7% miss) but only
    # 15.1 mg*h/L under the CV reading (a 36% miss), and reproduce the printed
    # PTAs to a mean absolute error of 3.5 percentage points versus 10.6.
    # The full arithmetic is in the vignette "Assumptions and deviations".
    etalcl ~ 0.224 # Song 2024 Table 2: omega CL (%) 22.4, read as omega^2 = 0.224
    etalvc ~ 0.552 # Song 2024 Table 2: omega V1 (%) 55.2, read as omega^2 = 0.552

    # Residual error: proportional (Methods, "a proportional model [was] used
    # to describe ... residual variability"; Table 2 footnote, "sigma,
    # residual variability for proportional error"). Table 2 prints
    # "Residual variability sigma (%) 1.58", read on the same variance scale
    # as the omegas above: sigma^2 = 0.0158, so propSd = sqrt(0.0158) = 0.126.
    # Read literally as a 1.58% proportional CV the value is below the
    # precision of the LC-MS/MS assay itself and is falsified by the scatter
    # in the paper's own Figure 1 DV-vs-IPRED panel.
    propSd <- 0.1257; label("Proportional residual error (fraction)") # Song 2024 Table 2: sigma = 1.58 (%), read as sigma^2 = 0.0158
  })
  model({
    # 1. Individual PK parameters. No covariates are retained in the final
    # model, so each parameter is the typical value times its exponential eta.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Tigecycline is given as a 30-minute intravenous infusion
    # into the central compartment (protocol carried from the companion study,
    # Su et al Front Pharmacol 2024;15:1342947).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observation and error. Dose in mg, volumes in L -> central/vc has
    # units mg/L, which is also the unit in which the paper's PK/PD targets
    # are expressed (AUC/MIC with MIC in mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
