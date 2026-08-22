Ozdin_2025_dexamethasone <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "dexamethasone released from dexamethasone sodium phosphate (DSP)",
    "encapsulated in autologous erythrocytes (eDSP, the EryDex System) given",
    "as a monthly intravenous infusion, in 133 pooled subjects: 18 healthy",
    "adults from a phase 1 single-dose study (NCT01925859) and 115 patients",
    "with ataxia telangiectasia (109 pediatric, 6 adult) from the phase 3",
    "ATTeST study (NCT02770807). Dexamethasone is released from the infused",
    "red blood cells over 20-30 days, so eDSP behaves as a sustained-release",
    "depot; the authors evaluated depot, dual-absorption, Michaelis-Menten and",
    "semi-physiological release structures and retained the simpler",
    "two-compartment disposition model, which estimates systemic exposure",
    "adequately without reproducing the underlying triphasic release. Body",
    "weight scales clearances (exponent fixed to 0.75 on CL and Q) and volumes",
    "(exponent fixed to 1 on V1 and V2) about a 70 kg reference. Clearance is",
    "about 10 percent lower in patients with ataxia telangiectasia than in",
    "healthy adults. Between-subject variability is carried on CL and V1 only.",
    "Combined proportional plus additive residual error is stratified between",
    "the phase 1 healthy cohort and the phase 3 patient cohort."
  )
  reference <- paste(
    "Ozdin D, Kheibarshekan L, Mambrini G, Tremblay PO. Population",
    "Pharmacokinetic Modeling and Pediatric Exposure of Dexamethasone Sodium",
    "Phosphate Encapsulated in Erythrocytes (eDSP) Administered Monthly for",
    "Treatment of Neurological Symptoms of Patients With Ataxia",
    "Telangiectasia. CPT Pharmacometrics Syst Pharmacol. 2025;14(11):1882-1893.",
    "doi:10.1002/psp4.70103"
  )
  vignette <- "Ozdin_2025_dexamethasone"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The canonical residual-error matcher recognises only the bare names
  # propSd / addSd / expSd, so the four cohort-stratified residual SDs are
  # declared here (same pattern as Shoji_2011_pregabalin.R).
  paper_specific_residual_sds <- c(
    "propSdHealthy", "addSdHealthy",
    "propSdPatient", "addSdPatient"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  #
  # The dose administered is the mass of DSP (the phosphate prodrug)
  # encapsulated in the final infusion bag, but the measured analyte is
  # dexamethasone. The final model carries no molecular-weight conversion:
  # Ozdin 2025 supplement Figure S4 introduces FD = 0.76 (516.4 / 392.464
  # g/mol) only in the REJECTED semi-physiological model, and the retained
  # two-compartment model (supplement "Final PK Model - NONMEM File") doses
  # the central compartment directly with the DSP amount. The states
  # therefore hold dexamethasone amount expressed on the DSP mass basis.
  compartmentData <- list(
    central     = list(analyte = "dexamethasone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dexamethasone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline (time-fixed) body weight, source column WTBL. Allometric",
        "power scaling about a 70 kg reference, with all four exponents held",
        "FIXED at their canonical values rather than estimated: 0.75 shared by",
        "CL and Q, 1 shared by V1 and V2 (Ozdin 2025 Table 2 rows 'Weight",
        "effect on CL/V1/Q/V2', each marked FIXED; Methods section 2.4). The",
        "supplement control stream writes the exponents as literal constants,",
        "(WTBL/70)**0.75 and (WTBL/70)**1. Cohort range 15.1-95.7 kg (overall",
        "median 25.0 kg, Table 1); the single 2-to-under-6-years subject",
        "weighed 15.1 kg. Weight was the only covariate retained in the final",
        "model besides disease status, and it was fixed during structural base",
        "model development rather than tested in the stepwise covariate search."
      ),
      source_name        = "WTBL"
    ),
    DIS_HEALTHY = list(
      description        = paste(
        "Healthy-participant cohort indicator; 0 = patient with ataxia",
        "telangiectasia"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (patient with ataxia telangiectasia enrolled in the phase 3 ATTeST",
        "study, NCT02770807)"
      ),
      notes              = paste(
        "Time-fixed per subject. Two distinct roles in this model.",
        "(1) Multiplicative power-form effect on clearance:",
        "e_dis_healthy_cl^DIS_HEALTHY with e_dis_healthy_cl = 1 / 0.899 =",
        "1.112, i.e. healthy adults clear dexamethasone about 11 percent",
        "faster than patients with ataxia telangiectasia at the same body",
        "weight (equivalently, about 10 percent lower CL in patients, the",
        "wording used in Ozdin 2025 Results section 3.3). Disease status was",
        "the only covariate that reached significance in the stepwise search",
        "(forward p < 0.01, backward p < 0.001) and the only one retained.",
        "(2) Binary stratifier selecting which of the two combined",
        "proportional-plus-additive residual error magnitudes applies",
        "(propSdHealthy / addSdHealthy versus propSdPatient / addSdPatient;",
        "Ozdin 2025 Table 2 'Residual error' rows, labelled Phase 1 and",
        "Phase 3).",
        "Orientation: the source uses two reverse-coded columns, PTNT",
        "(0 = healthy, 1 = patient) for the clearance effect and PHAS",
        "(1 = phase 1, 3 = phase 3) for the residual-error branch. Every",
        "phase 1 subject is a healthy adult volunteer and every phase 3",
        "subject is an ataxia-telangiectasia patient (Methods section 2.1;",
        "Table 1 tabulates the two studies as 'Phase 1-Healthy volunteers'",
        "and 'Phase 3 patients with ataxia telangiectasia'), so PTNT and PHAS",
        "partition the pooled dataset identically and both are re-expressed",
        "on the single canonical column as DIS_HEALTHY = 1 - PTNT. The",
        "structural typical lcl is shifted to the patient state",
        "log(20.8 * 0.899) = log(18.70) so that DIS_HEALTHY = 1 restores the",
        "printed healthy typical of 20.8 L/h, following the",
        "Chen_2023_nemonoxacin.R and Yoneyama_2017_emicizumab.R precedents.",
        "Cohort split: 18 healthy (13.5 percent), 115 patients (86.5",
        "percent).",
        "The authors attribute the slightly lower patient clearance to known",
        "fluctuation of CYP3A4 activity with age and body composition",
        "(Results section 3.5)."
      ),
      source_name        = "PTNT"
    )
  )

  # Covariates screened during the stepwise covariate analysis (Methods
  # section 2.4) and NOT retained in the final model. Ozdin 2025 Results
  # section 3.3 reports that the ETA-versus-covariate plots showed no
  # statistically significant relationship for any of these, and only
  # disease status survived forward selection and backward elimination.
  # Documentation only -- these names are deliberately not referenced in
  # model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and V1; not retained (Results section 3.3, 'no",
        "statistically significant relationships were observed between PK",
        "parameters and continuous covariates such as BMI, age, and dose').",
        "Cohort range 5.00-55.0 years, overall median 10.0 years (Table 1).",
        "Age group (2 to under 6, 6 to under 10, 10 to under 17, adult) was",
        "screened as a separate categorical stratification of this column and",
        "was also not retained. No maturation function was incorporated;",
        "Results section 3.5 cautions that the 2-to-under-6-years group rests",
        "on a single subject."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and V1; not retained (Results section 3.3). Cohort",
        "range 9.60-29.5 kg/m^2, overall median 15.2 kg/m^2 (Table 1).",
        "Correlated with WT, and the Methods state that highly correlated",
        "covariates were not entered on the same parameter."
      )
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Tabulated in the baseline covariate summary (Table 1; cohort range",
        "66.8-184 cm, overall median 129 cm) and used to derive BMI, but not",
        "tested as a covariate in its own right and not in the final model."
      )
    ),
    SEXF = list(
      description = "Sex",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL and V1; not retained (Results section 3.3, 'the ETA",
        "versus categorical covariate plots did not indicate any",
        "statistically significant relationships'). Cohort 62 of 133 female",
        "(47 percent), 71 male (53 percent) (Table 1)."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL and V1; not retained (Results section 3.3). Cohort 76",
        "of 133 White (57 percent) (Table 1)."
      )
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL and V1; not retained (Results section 3.3). Cohort 57",
        "of 133 African American and/or Black (43 percent) (Table 1). Race is",
        "fully dichotomous in this pooled dataset, so RACE_WHITE and",
        "RACE_BLACK are complementary."
      )
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL and V1 as the paper's 'ethnicity' covariate; not",
        "retained (Results section 3.3). Cohort 6 of 133 Hispanic or Latino",
        "(5 percent), 127 non-Hispanic or Latino (95 percent) (Table 1)."
      )
    ),
    DOSE_DSP_MG = list(
      description = "Mass of DSP encapsulated in red blood cells per infusion",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "Screened as the paper's 'dose level' covariate on CL and V1; not",
        "retained (Results section 3.3), consistent with the",
        "dose-proportional AUC reported in the Discussion. Overall mean 11.6",
        "mg (CV 56.3 percent) (Table 1); nominal levels are half-low 4.2 mg,",
        "low 8.3 mg and high 17.4 mg (Methods section 2.6). This quantity is",
        "the dose amount itself, so it is carried on the dosing records rather",
        "than as a covariate column."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 133,
    n_studies      = 2,
    age_range      = "5.00-55.0 years",
    age_median     = "10.0 years",
    weight_range   = "15.1-95.7 kg",
    weight_median  = "25.0 kg",
    sex_female_pct = 47,
    race_ethnicity = c(White = 57, Black = 43),
    disease_state  = paste(
      "ataxia telangiectasia (115 subjects, phase 3 ATTeST NCT02770807) and",
      "healthy adults (18 subjects, phase 1 NCT01925859)"
    ),
    dose_range     = paste(
      "4.2, 8.3 and 17.4 mg DSP encapsulated in red blood cells per",
      "intravenous infusion; single dose over approximately 10 min in phase",
      "1, monthly infusions over approximately 40 min in phase 3"
    ),
    regions        = "international, multi-center (phase 3); single-center (phase 1)",
    notes          = paste(
      "Baseline demographics from Ozdin 2025 Table 1. Age groups: adults",
      "n = 24 (18 healthy from phase 1 plus 6 patients from phase 3),",
      "10 to under 17 years n = 44, 6 to under 10 years n = 64, and 2 to",
      "under 6 years n = 1. The single subject in the youngest group is the",
      "reason the paper repeatedly cautions against extrapolating to children",
      "aged 2 to under 6 years, and no maturation function was included.",
      "1719 dexamethasone concentrations were available, of which 68 percent",
      "were above the limit of quantification (Table S1); pre-dose",
      "concentrations were set to zero and one subject was excluded for",
      "missing dosing information (Methods section 2.3). Estimated in NONMEM",
      "7.4 with FOCE-I; final objective function value 8584.636."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural disposition parameters, Ozdin 2025 Table 2 "Typical values"
    # (final estimates). The supplement control stream $THETA block lists
    # INITIAL estimates only (CL 19.3324, V1 121.007, Q 0.323498, V2 9.94667,
    # PTNT 0.5), so Table 2 is the authority for every final value here.
    #
    # lcl is shifted from the paper's healthy-cohort typical of 20.8 L/h to
    # the ataxia-telangiectasia patient state (20.8 * 0.899 = 18.70 L/h) so
    # that the canonical DIS_HEALTHY = 1 term below restores 20.8 L/h. See
    # the DIS_HEALTHY covariateData notes for the orientation argument.
    lcl <- log(20.8 * 0.899); label("Clearance in patients with ataxia telangiectasia at 70 kg (L/h)")  # Table 2: CL 20.8 L/h (RSE 8%) x theta_PTNT 0.899
    lvc <- log(122);          label("Central volume of distribution at 70 kg (L)")                      # Table 2: V1 122 L (RSE 9%)
    lq  <- log(0.358);        label("Inter-compartmental clearance at 70 kg (L/h)")                     # Table 2: Q 0.358 L/h (RSE 13%)
    lvp <- log(16.8);         label("Peripheral volume of distribution at 70 kg (L)")                   # Table 2: V2 16.8 L (RSE 9%)

    # ---------------------------------------------------------------
    # Allometric exponents on body weight about a 70 kg reference. All four
    # were held FIXED at canonical values, not estimated (Table 2 marks each
    # "FIXED" with no RSE and no bootstrap CI; Methods section 2.4 "Body
    # weight was incorporated as a fixed allometric scaling factor"). CL and
    # Q share the 0.75 exponent; V1 and V2 share the exponent of 1.
    e_wt_cl_q   <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)")   # Table 2: "Weight effect on CL/Q -> * (Weight/70)^theta", 0.75 FIXED
    e_wt_vc_vp  <- fixed(1);    label("Allometric exponent of body weight on V1 and V2 (unitless)")  # Table 2: "Weight effect on V1/V2 -> * (Weight/70)^theta", 1 FIXED

    # ---------------------------------------------------------------
    # Disease-status effect on clearance. The source form is
    # TVCL = THETA(1) * (WTBL/70)**0.75 * THETA(9)**PTNT with THETA(9) =
    # 0.899 and PTNT = 1 for patients (supplement control stream $PK;
    # Table 2 row "Patient effect on CL -> * theta_PTNT"). Re-expressed on
    # the canonical DIS_HEALTHY = 1 - PTNT orientation, the multiplier
    # becomes 1 / 0.899 = 1.112 raised to DIS_HEALTHY.
    e_dis_healthy_cl <- 1 / 0.899; label("Multiplicative effect of healthy status on CL (ratio)")  # Table 2: theta_PTNT 0.899 (RSE 10%), reciprocal for the DIS_HEALTHY orientation

    # ---------------------------------------------------------------
    # Between-subject variability, Ozdin 2025 Table 2 "Between Subject
    # Variability (CV%)". The table does not state whether these are exact
    # log-normal CVs or the approximation 100 * sqrt(omega^2); the exact-CV
    # reading is used, so the internal variance is omega^2 = log(1 + CV^2).
    # The decisive evidence is that log(1 + 0.707^2) = 0.4054 reproduces the
    # supplement control stream's $OMEGA(2,2) starting value of 0.404706 to
    # within 0.2 percent, where the alternative reading (0.4998) would be 23
    # percent away. That agreement is meaningful rather than coincidental
    # because the control stream's starting values are evidently a previous
    # run's near-converged estimates: the phase 3 residual-error starting
    # values (PROP2 0.296021, ADD2 24.2317) are likewise within 0.5 percent of
    # the Table 2 finals. The vignette additionally simulates both readings
    # against the Table 3 exposure spread; that check favours this reading in
    # all three age bands but only by one to two Monte Carlo standard errors,
    # so it corroborates rather than decides. IIV was carried on CL and V1 only -- the
    # control stream fixes $OMEGA to 0 for Q and V2, and Results section 3.3
    # states that adding IIV on the second-compartment parameters did not
    # improve the fit. The $OMEGA block is diagonal, so CL and V1 are
    # uncorrelated.
    etalcl ~ log(1 + 0.414^2); label("IIV on clearance")                    # Table 2: ETA_CL 41.4% CV (RSE 51%, shrinkage 26%) -> omega^2 = 0.1582
    etalvc ~ log(1 + 0.707^2); label("IIV on central volume of distribution") # Table 2: ETA_V1 70.7% CV (RSE 21%, shrinkage 7%) -> omega^2 = 0.4054

    # ---------------------------------------------------------------
    # Residual variability. Combined proportional plus additive error, with
    # the magnitudes stratified by study cohort. The supplement control
    # stream $ERROR block is
    #   IF (PHAS.EQ.1) W = (PROP1*PROP1*IPRED*IPRED + ADD1*ADD1)**0.5
    #   IF (PHAS.EQ.3) W = (PROP2*PROP2*IPRED*IPRED + ADD2*ADD2)**0.5
    #   Y = IPRED + W*EPS(1)          with $SIGMA 1 FIX
    # so PROP1/PROP2 are proportional SDs on the fraction scale and
    # ADD1/ADD2 are additive SDs in ng/mL. Phase 1 maps onto the healthy
    # cohort and phase 3 onto the patient cohort (see DIS_HEALTHY notes).
    # The very small phase 1 additive SD reflects the dense sampling out to
    # day 42 in healthy adults, where the terminal concentrations are orders
    # of magnitude below the phase 3 trough range.
    propSdHealthy <- 0.181;  label("Proportional residual SD, healthy phase 1 cohort (fraction)")  # Table 2: Proportional Error (%)-Phase 1 = 18.1% (RSE 10%)
    addSdHealthy  <- 0.0013; label("Additive residual SD, healthy phase 1 cohort (ng/mL)")         # Table 2: Additive error (ng/mL)-Phase 1 = 0.0013 (RSE 18%)
    propSdPatient <- 0.295;  label("Proportional residual SD, patient phase 3 cohort (fraction)")  # Table 2: Proportional Error (%)-Phase 3 = 29.5% (RSE 23%)
    addSdPatient  <- 24.3;   label("Additive residual SD, patient phase 3 cohort (ng/mL)")         # Table 2: Additive error (ng/mL)-Phase 3 = 24.3 (RSE 58%)
  })

  model({
    # ---------------------------------------------------------------
    # Individual PK parameters. Supplement control stream $PK block:
    #   TVCL = THETA(1) * (WTBL/70)**0.75 * THETA(9)**PTNT ; CL = TVCL*EXP(ETA(1))
    #   TVV1 = THETA(2) * (WTBL/70)**1                     ; V1 = TVV1*EXP(ETA(2))
    #   TVQ  = THETA(3) * (WTBL/70)**0.75                  ; Q  = TVQ *EXP(ETA(3))
    #   TVV2 = THETA(4) * (WTBL/70)**1                     ; V2 = TVV2*EXP(ETA(4))
    # ETA(3) and ETA(4) sit on $OMEGA blocks fixed to 0, so Q and V2 carry
    # no between-subject variability here.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * e_dis_healthy_cl^DIS_HEALTHY
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    # Micro-constants for the two-compartment central-disposition model
    # (NONMEM ADVAN3 TRANS4: central (V1) <-> peripheral (V2), elimination
    # from central; supplement Figure S2).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. eDSP is infused intravenously, so the dose enters the
    # central compartment directly; the retained final model has no depot or
    # absorption process (Results section 3.2 rejects the depot,
    # dual-absorption and semi-physiological alternatives).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma dexamethasone concentration. Dose in mg and vc in L give
    # central / vc in mg/L = ug/mL, so the factor of 1000 converts to the
    # ng/mL reported by the bioanalytical assay (Methods section 2.1: assay
    # range 0.5-250 ng/mL) and used throughout Table 2, Table 3 and Tables
    # S1-S4. The source control stream instead sets S1 = V1 and carries the
    # dose amount in ug; the two encodings are identical.
    Cc <- 1000 * central / vc

    # ---------------------------------------------------------------
    # Cohort-stratified combined proportional plus additive residual error.
    # Exactly one of the two indicator products is 1 for any given subject.
    propSd <- propSdHealthy * DIS_HEALTHY + propSdPatient * (1 - DIS_HEALTHY)
    addSd  <- addSdHealthy  * DIS_HEALTHY + addSdPatient  * (1 - DIS_HEALTHY)

    Cc ~ prop(propSd) + add(addSd)
  })
}
