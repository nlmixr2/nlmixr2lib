Tamura_2023_remdesivir <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for intravenous remdesivir",
    "and its major circulating nucleoside metabolite GS-441524 in Japanese",
    "adults hospitalised with moderate-to-severe COVID-19 (Tamura 2023).",
    "One compartment for each compound, coupled in series, with first-order",
    "elimination throughout (NONMEM ADVAN13). Only the GS-441524 data were",
    "fitted: remdesivir fell below the assay limit of quantification within",
    "5 h of the infusion, so its renal clearance, metabolic clearance and",
    "volume could not be estimated and were all fixed to the values reported",
    "for cohort 5 (150 mg) of a previous single-dose study in healthy",
    "subjects. The parent-to-metabolite flux carries an explicit",
    "molecular-weight conversion (GS-441524 291.26 g/mol divided by",
    "remdesivir 602.58 g/mol) because the model is written in mass rather",
    "than molar units. The fraction of metabolised remdesivir that appears",
    "in plasma as measurable GS-441524 (Fm) is not identifiable from plasma",
    "data alone, so the metabolite clearance and volume are apparent",
    "parameters (source CLm/Fm and Vm/Fm) and the metabolite state carries",
    "an apparent amount; metabolite concentrations are nevertheless",
    "predicted correctly because the same unknown factor scales the apparent",
    "amount and the apparent volume. The only retained covariate is eGFR on",
    "GS-441524 apparent clearance (power, exponent 0.745, reference",
    "68 mL/min/1.73 m^2). Correlated between-subject variability was",
    "estimated on GS-441524 apparent clearance and apparent volume; residual",
    "error is proportional. The paper found no relationship between",
    "GS-441524 exposure and either recovery rate or transaminase elevation.")
  reference <- "Tamura R, Irie K, Nakagawa A, Muroi H, Eto M, Ikesue H, et al. Population pharmacokinetics and exposure-clinical outcome relationship of remdesivir major metabolite GS-441524 in patients with moderate and severe COVID-19. CPT Pharmacometrics Syst Pharmacol. 2023;12(4):513-521. doi:10.1002/psp4.12936"
  vignette <- "Tamura_2023_remdesivir"

  # The source NONMEM control stream (Supporting Information, "NONMEM model
  # code") sets S1 = VP and S2 = VM with volumes in L, and the accompanying
  # "Dataset sample" gives AMT = 200000 for the 200 mg day-1 dose. Amounts
  # are therefore in ug, so amount / volume returns ug/L = ng/mL directly,
  # which is the scale the paper reports concentrations on.
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # `central` holds a true remdesivir mass amount. `central_gs441524` holds
  # an APPARENT amount: the true GS-441524 mass divided by Fm, the
  # unidentifiable fraction of metabolised remdesivir reaching plasma as
  # GS-441524, which the source folds into its CLm/Fm and Vm/Fm parameters.
  # Dividing that apparent amount by the apparent volume `vc_gs441524`
  # returns the true plasma concentration, because the same unknown factor
  # scales numerator and denominator. Verified against the NONMEM control
  # stream in the Supporting Information (COMP=(CENTRAL), COMP=(METCENTRAL);
  # S1 = VP, S2 = VM).
  compartmentData <- list(
    central          = list(analyte = "remdesivir", units = "ug", specimen = "plasma", verified = TRUE),
    central_gs441524 = list(analyte = "GS-441524", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalised.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tamura 2023 Methods 'Patients' lists eGFR among the patient characteristics collected; the paper does not name the estimating equation, and the cohort is Japanese, so the Japanese Society of Nephrology eGFR equation is the likely but unstated source. Enters as a power effect on GS-441524 apparent clearance only, normalised to the study-population median of 68 mL/min/1.73 m^2, per the final-model equation printed in Results 'Covariate analysis': CLm/Fm (L/h) = 11.0 * (eGFR / 68)^0.745. In the source NONMEM control stream the column is named GFR and is time-fixed (one value per subject in the dataset sample). Observed range 33-113 mL/min/1.73 m^2 (Table 1), so the model is uninformed outside that interval; the paper's own Monte Carlo simulations exercise it only at 30, 68 and 113. Height and age also passed forward inclusion on CLm/Fm but were dropped in backward elimination.",
      source_name        = "GFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Tamura 2023 Methods 'Covariate analysis') but not retained in the final model; no point estimate is reported."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Passed forward inclusion (p < 0.05) on CLm/Fm but was removed in backward elimination (p < 0.01); no point estimate is reported (Tamura 2023 Results 'Covariate analysis')."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Reported as 'BH' in Results 'Covariate analysis'. Passed forward inclusion on CLm/Fm but was removed in backward elimination; no point estimate is reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened but not retained (Tamura 2023 Methods 'Covariate analysis')."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 2.8 g/dL (range 1-4.2), Table 1."
    ),
    BILI = list(
      description = "Total serum bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 0.5 mg/dL (range 0.2-2.9), Table 1."
    ),
    AST = list(
      description = "Serum aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 45 IU/L (range 14-276), Table 1."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 36 IU/L (range 9-130), Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained; renal function entered the final model through eGFR (CRCL) instead."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Tamura 2023 Methods 'Covariate analysis')."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 1.74 m^2 (range 1.36-2.03), Abstract and Results."
    ),
    WHO_ORDINAL = list(
      description = "WHO clinical-status ordinal score (1-7)",
      units       = "(ordinal)",
      type        = "categorical",
      notes       = "Screened as 'clinical status' but not retained. Cohort distribution: score 5 in 17 patients (43.6%), score 6 in 7 (17.9%), score 7 in 15 (38.5%), Table 1. Not a registered canonical covariate column; listed here as documentation of the covariate screen only."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    n_observations = "102 serum samples across the two analytes in 39 patients (Abstract; Results 'Base model'). Only the GS-441524 concentrations entered the population PK fit -- every remdesivir sample drawn more than 5 h after the infusion was below the 10 ng/mL limit of quantification, so the remdesivir parameters were fixed rather than estimated (Methods 'Pharmacokinetic analysis').",
    age_range      = "42-85 years",
    age_median     = "70 years",
    weight_range   = "41.8-84 kg",
    weight_median  = "65.2 kg",
    sex_female_pct = 25.6,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "Moderate-to-severe COVID-19 confirmed by real-time PCR, all hospitalised and all requiring oxygen: WHO ordinal score 5 (requiring oxygen) in 17 patients (43.6%), score 6 (high-flow oxygen or non-invasive ventilation) in 7 (17.9%) and score 7 (invasive ventilation and/or ECMO) in 15 (38.5%). Recovery ratio at day 28 was 56.1% and mortality 7.7%.",
    renal_function = "Median eGFR 68 mL/min/1.73 m^2 (range 33-113): 4 patients (10%) at 30-44, 10 (25%) at 45-59, 19 (49%) at 60-89 and 6 (15%) at 90 or above (Table 1). No patient received renal replacement therapy in the analysed cohort.",
    hepatic_function = "Median AST 45 IU/L (range 14-276), ALT 36 IU/L (9-130), total bilirubin 0.5 mg/dL (0.2-2.9), albumin 2.8 g/dL (1-4.2) (Table 1).",
    dose_range     = "Licensed regimen only: remdesivir 200 mg intravenously on day 1 followed by 100 mg once daily on days 2-5, each infused over 60 min. One patient discontinued on day 4.",
    regions        = "Single centre, Kobe City Medical Center General Hospital, Kobe, Japan; 16 May 2020 to 31 March 2021.",
    notes          = "Retrospective observational study using residual serum from routine arterial blood-gas testing, so sampling was opportunistic rather than protocol-scheduled. Concentrations were measured by LC-MS/MS over a 10-2000 ng/mL calibration range (LOQ 10 ng/mL). Baseline demographics are in Tamura 2023 Table 1. Estimation was FOCE-I in NONMEM 7.4.1 with ADVAN13; the final model was checked by a 1000-resample nonparametric bootstrap (99.2% success) and a 1000-replicate prediction-corrected VPC using PsN 4.9.0. Note that samples were ARTERIAL, which the authors flag as a possible cause of the positive bias in the remdesivir goodness-of-fit plot immediately after dosing (Discussion, limitations)."
  )

  ini({
    # ------------------------------------------------------------------
    # REMDESIVIR (PARENT) STRUCTURAL PARAMETERS -- ALL FIXED
    # ------------------------------------------------------------------
    # Tamura 2023 Methods 'Pharmacokinetic analysis': "Renal clearance
    # (CLp) and distribution volume (Vp) of remdesivir and metabolic
    # clearance of remdesivir to GS-441524 conversion (CLpm) could not be
    # estimated because the PK data for remdesivir are limited. Therefore,
    # the PK data of remdesivir measured in this study was not included in
    # the PopPK analysis, and the PK parameters of remdesivir (CLp, CLpm,
    # and Vp) were set according to a previous single dose study (cohort 5,
    # 150 mg)". All three carry an explicit FIX flag in Table 2 and in the
    # $THETA block of the control stream in the Supporting Information.
    #
    # The two clearance arms are an exact metabolic-formation / non-
    # formation decomposition of total parent clearance: CLpm is the flux
    # that forms GS-441524 and CLp is renal excretion of unchanged
    # remdesivir, and the control stream uses precisely that split
    # (K10 = CLP/VP for loss out of the system, K12 = CLPM/VP into the
    # metabolite state).

    lcl_met <- fixed(log(51.8))
    label("Remdesivir metabolic (GS-441524-forming) clearance CLpm (L/h)")     # Tamura 2023 Table 2: CLpm = 51.8 L/h FIX (= 863 mL/min); Methods 'Pharmacokinetic analysis' and Table 2 footnote, from cohort 5 (150 mg) of the previous single-dose study

    lcl_nonmet <- fixed(log(4.69))
    label("Remdesivir renal (non-forming) clearance CLp (L/h)")                # Tamura 2023 Table 2: CLp = 4.69 L/h FIX (= 78.1 mL/min); Methods 'Pharmacokinetic analysis' and Table 2 footnote

    lvc <- fixed(log(73.4))
    label("Remdesivir distribution volume Vp (L)")                             # Tamura 2023 Table 2: Vp = 73.4 L FIX; Methods 'Pharmacokinetic analysis' and Table 2 footnote

    # ------------------------------------------------------------------
    # GS-441524 (METABOLITE) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Apparent parameters: the source reports CLm/Fm and Vm/Fm because the
    # metabolic ratio Fm is not identifiable from plasma data alone. Both
    # are carried here exactly as reported; predicted metabolite
    # concentrations are unaffected by the unknown scaling (see the
    # compartmentData note above). The clearance is the typical value at
    # the covariate reference point, eGFR = 68 mL/min/1.73 m^2.

    lcl_gs441524 <- log(11.0)
    label("GS-441524 apparent clearance CLm/Fm (L/h) at eGFR 68 mL/min/1.73 m^2")  # Tamura 2023 Table 2: CLm/Fm = 11.0 L/h (RSE 6.6%; bootstrap median 11.1, 2.5th-97.5th 9.5-12.5). Reference eGFR from the Results 'Covariate analysis' final-model equation.

    lvc_gs441524 <- log(271)
    label("GS-441524 apparent distribution volume Vm/Fm (L)")                  # Tamura 2023 Table 2: Vm/Fm = 271 L (RSE 10.8%; bootstrap median 268, 2.5th-97.5th 210-336)

    # ------------------------------------------------------------------
    # COVARIATE EFFECT
    # ------------------------------------------------------------------
    # The single retained covariate, as a power model on the metabolite's
    # apparent clearance. Tamura 2023 Results 'Covariate analysis' prints
    # the final-model equation as
    #   CLm/Fm (L/h) = 11.0 * (eGFR / 68)^0.745
    # and the control stream codes it as
    #   CLM = THETA(3)*EXP(ETA(1))*(GFR/68)**THETA(6).
    # Age and height passed forward inclusion but were dropped in backward
    # elimination, so no covariate acts on Vm/Fm or on any remdesivir
    # parameter.

    e_crcl_cl_gs441524 <- 0.745
    label("Power exponent for eGFR on GS-441524 apparent clearance (unitless)")  # Tamura 2023 Table 2 "CLm/Fm on eGFR" = 0.745 (RSE 32.9%; bootstrap median 0.820, 2.5th-97.5th 0.317-1.29); same exponent in the Abstract and the Results equation

    # ------------------------------------------------------------------
    # BETWEEN-SUBJECT VARIABILITY
    # ------------------------------------------------------------------
    # Tamura 2023 Methods 'Pharmacokinetic analysis': "A log-normally
    # distributed intersubject variability (ISV) model was used ... and was
    # included in CLm/Fm and Vm/Fm", theta_i = theta_Tv * exp(eta_i) with
    # "variance of omega^2". Results 'Base model' adds that "the
    # interindividual random effects were considered to be CLm/Fm, Vm/Fm,
    # and their correlation", and the control stream declares
    # $OMEGA BLOCK(2). Table 2 therefore reports variances on the log
    # scale, not %CV: the Abstract's ISV% of 43.0% and 58.1% are recovered
    # as sqrt(0.185) = 0.4301 and sqrt(0.338) = 0.5814, which is the
    # sqrt(omega^2) convention rather than sqrt(exp(omega^2) - 1).
    #
    # The off-diagonal is read as the OMEGA(2,1) COVARIANCE element, which
    # is what a NONMEM $OMEGA BLOCK(2) reports and what the row label
    # "omega_CLm x omega_Vm" denotes. It implies a correlation of
    # 0.240 / sqrt(0.185 * 0.338) = 0.960 -- high, but mechanistically
    # expected here, because CLm/Fm and Vm/Fm share the same unidentifiable
    # Fm in their denominators and therefore share its variability. The
    # resulting matrix is positive definite (determinant 0.00493).

    etalcl_gs441524 + etalvc_gs441524 ~ c(0.185, 0.240, 0.338)                 # Tamura 2023 Table 2: omega^2 CLm = 0.185 (RSE 39.7%, shrinkage 3.9%), omega_CLm x omega_Vm = 0.240 (RSE 8.2%), omega^2 Vm = 0.338 (RSE 30.8%, shrinkage 6.1%); NONMEM $OMEGA BLOCK(2) ordering

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Tamura 2023 Results 'Base model': "the proportional residual error
    # model was selected". The control stream codes it as
    # Y = F*(1+EPS(1)), so Table 2's sigma^2 is the variance of a
    # proportional error and the nlmixr2 SD is its square root:
    # sqrt(0.0378) = 0.19442. Applies to GS-441524 only, the sole fitted
    # analyte.

    propSd_gs441524 <- 0.19442
    label("GS-441524 proportional residual SD (fraction)")                     # Tamura 2023 Table 2: sigma^2 proportional error = 0.0378 (RSE 19.6%, shrinkage 19.9%; bootstrap median 0.036, 2.5th-97.5th 0.023-0.051) -> SD = sqrt(0.0378)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters
    # ------------------------------------------------------------------
    # No covariate and no between-subject variability acts on any
    # remdesivir parameter -- all three were fixed from external data.
    cl_met    <- exp(lcl_met)
    cl_nonmet <- exp(lcl_nonmet)
    vc        <- exp(lvc)

    cl_gs441524 <- exp(lcl_gs441524 + etalcl_gs441524) *
      (CRCL / 68)^e_crcl_cl_gs441524
    vc_gs441524 <- exp(lvc_gs441524 + etalvc_gs441524)

    # ------------------------------------------------------------------
    # Micro-constants
    # ------------------------------------------------------------------
    # Transcribed one-for-one from the $PK block of the control stream in
    # the Supporting Information:
    #   K10 = CLP/VP    (renal loss of unchanged remdesivir)
    #   K12 = CLPM/VP   (conversion to GS-441524)
    #   K20 = CLM/VM    (GS-441524 elimination)
    kel          <- cl_nonmet / vc
    kform        <- cl_met / vc
    kel_gs441524 <- cl_gs441524 / vc_gs441524

    # Mass-to-mass stoichiometric conversion for the formation flux. The
    # model is written in mass units (ug), so one mass unit of remdesivir
    # metabolised yields MW(GS-441524) / MW(remdesivir) mass units of
    # GS-441524. Hardcoded in the $DES block of the control stream as the
    # literal factor 291.26/602.58; kept in that form here so the
    # transcription is checkable against the source line.
    mwRatio <- 291.26 / 602.58

    # ------------------------------------------------------------------
    # ODE system
    # ------------------------------------------------------------------
    # Control stream $DES, transcribed one-for-one:
    #   DADT(1) = -K12*A(1) - K10*A(1)
    #   DADT(2) =  K12*A(1)*291.26/602.58 - K20*A(2)
    # Remdesivir is dosed into `central` as a 60-minute intravenous
    # infusion (the dataset sample sets AMT = RATE, giving a 1 h duration).
    d/dt(central)          <- -(kel + kform) * central
    d/dt(central_gs441524) <-  mwRatio * kform * central -
      kel_gs441524 * central_gs441524

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    # Amounts in ug divided by volumes in L give ug/L = ng/mL, the scale
    # the paper reports on. Only GS-441524 carries a residual error model:
    # the remdesivir observations were excluded from the fit, so `Cc` is a
    # prediction from wholly fixed parameters rather than a fitted output.
    Cc          <- central          / vc
    Cc_gs441524 <- central_gs441524 / vc_gs441524

    Cc_gs441524 ~ prop(propSd_gs441524)
  })
}
