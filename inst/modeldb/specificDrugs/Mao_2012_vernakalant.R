Mao_2012_vernakalant <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous vernakalant",
    "hydrochloride (RSD1235) in 605 pooled subjects: 597 patients with",
    "atrial fibrillation or atrial flutter from five phase 2/3 trials",
    "(Scene 2, ACT I-IV) and 8 healthy volunteers from a phase 1 study.",
    "First-order elimination from the central compartment. All structural",
    "parameters are estimated per kilogram of body weight (CL and Q in",
    "L/h/kg, Vc and Vp in L/kg), so body weight enters linearly rather",
    "than allometrically. Clearance carries CYP2D6 phenotype, subject",
    "status (patient vs healthy), age, and serum creatinine effects;",
    "intercompartmental clearance carries a subject-status effect; the",
    "volumes carry no covariates. Interindividual variability is a full",
    "4x4 correlated block on CL, Vc, Vp, and Q, and the proportional",
    "residual error is stratified into healthy-volunteer and patient",
    "strata. Parameters are for vernakalant FREE BASE, so dose amounts",
    "must be supplied as free base (see units notes)."
  )
  reference <- paste(
    "Mao ZL, Townsend RW, Gao Y, Wheeler JJ, Kastrissios H, Keirns J.",
    "Population pharmacokinetics of vernakalant hydrochloride injection",
    "(RSD1235) in patients with atrial fibrillation or atrial flutter.",
    "J Clin Pharmacol. 2012;52(7):1042-1053.",
    "doi:10.1177/0091270011408425",
    sep = " "
  )
  vignette <- "Mao_2012_vernakalant"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are mg of vernakalant FREE BASE -- a dose
  # labelled as vernakalant hydrochloride must be divided by 1.104 (the
  # 385.93 / 349.50 g/mol salt-to-base molecular weight ratio given in
  # Mao 2012 Results, "Model Development") before it is passed as `amt`.
  compartmentData <- list(
    central     = list(analyte = "vernakalant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vernakalant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mao 2012 estimated every structural parameter on a per-kilogram",
        "basis (CL and Q in L/h/kg, Vc and Vp in L/kg), so weight enters",
        "as a LINEAR multiplier rather than an allometric power. This",
        "follows from the dosing design: vernakalant was administered as",
        "mg/kg in all five AF/AFL trials, and the Discussion states that",
        "'body weight did not have an effect on the pharmacokinetics",
        "because the dose of vernakalant hydrochloride was administered",
        "on a per kilogram basis'. Body weight was screened as an",
        "additional covariate on Vp during univariate screening but was",
        "dropped in backward elimination. Cohort weight means by study",
        "range 75.8 to 87.4 kg (Mao 2012 Table I). Time-fixed at",
        "baseline. Because dose, CL, Q, Vc, and Vp all scale linearly",
        "with weight, simulated concentrations after a mg/kg dose are",
        "weight-invariant in this model."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL, entered on the log scale as",
        "theta8 * log(AGE / 60) so that CL scales as (AGE / 60)^-0.284.",
        "The reference of 60 years is the base case used throughout Mao",
        "2012's diagnostic and sensitivity analyses (Figure 3 caption:",
        "'normalized for patients, aged 60 years'; Sensitivity Analysis:",
        "'the base case was a 60-year-old patient'). The negative",
        "exponent means older subjects have LOWER clearance and",
        "therefore higher exposure. Sensitivity range explored by the",
        "paper: 38 years (5th percentile) to 81 years (95th percentile).",
        "Cohort age means: 63.4 years for the AF/AFL patients and 31.2",
        "years for the healthy volunteers (Mao 2012 Table I). Time-fixed",
        "at baseline."
      ),
      source_name        = "Age"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL, entered on the log scale as",
        "theta9 * log(CREAT / 100) so that CL scales as",
        "(CREAT / 100)^-0.345. The reference of 100 umol/L is the base",
        "case used throughout Mao 2012's diagnostic and sensitivity",
        "analyses (Figure 3 caption: 'serum creatinine of 100 umol/L';",
        "Sensitivity Analysis base case). The negative exponent means",
        "subjects with higher serum creatinine (poorer renal function)",
        "have LOWER clearance and therefore higher exposure. Sensitivity",
        "range explored by the paper: 62 umol/L (5th percentile) to 133",
        "umol/L (95th percentile). Cohort means by study range 83.5 to",
        "105.3 umol/L (Mao 2012 Table I). NOTE the unit: this model uses",
        "umol/L, not mg/dL; divide mg/dL by 0.0113 to convert.",
        "Time-fixed at baseline."
      ),
      source_name        = "serum creatinine"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with atrial fibrillation or atrial flutter)",
      notes              = paste(
        "Mao 2012 calls this covariate 'subject status' and codes it as",
        "(Subject = Volunteer), which is already the canonical",
        "DIS_HEALTHY orientation -- no re-expression was needed. It has",
        "three distinct roles in this model. (1) It gates the",
        "intercompartmental clearance directly:",
        "Q = exp(theta4 + theta10 * DIS_HEALTHY) * WT, giving 1.34",
        "L/h/kg in patients and 2.63 L/h/kg in healthy volunteers. (2)",
        "It selects which CYP2D6 extensive-metabolizer coefficient",
        "applies to CL -- theta5 for healthy EM subjects and theta6 for",
        "patient EM/UM subjects -- so the CYP2D6 effect on clearance is",
        "an INTERACTION with subject status, not a main effect. (3) It",
        "stratifies the proportional residual error: healthy volunteers",
        "use propSd_hv = sqrt(0.056) = 0.237 and patients use",
        "propSd_pt = sqrt(0.079) = 0.281 (Mao 2012 Table II, sigma^2",
        "Volunteer and sigma^2 Patient rows, reported as 24% and 28%).",
        "Cohort split: 8 healthy volunteers (1.3%) from the phase 1",
        "study and 597 AF/AFL patients (98.7%) from Scene 2 and ACT",
        "I-IV. The paper explicitly cautions that 'this comparison is",
        "limited by the small number of volunteers included in the study",
        "cohort'. The healthy stratum is also fully confounded with age",
        "(volunteer mean 31.2 vs patient mean 63.4 years), sex (all 8",
        "volunteers were men) and dosing (a flat 240 mg rather than",
        "mg/kg)."
      ),
      source_name        = "Subject = Volunteer"
    ),
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (poor, ultrarapid, or ungenotyped subject)",
      notes              = paste(
        "Part of the three-indicator CYP2D6 encoding used here",
        "(CYP2D6_EM, CYP2D6_PM, CYP2D6_UM). Mao 2012 fits FOUR clearance",
        "strata and the IMPLICIT REFERENCE -- all three indicators 0 --",
        "is the 'CYP2D6 genotype NOT COLLECTED' group, not a phenotype:",
        "genotype was obtained for only 311 of 605 subjects (51.4%), and",
        "the paper estimates a separate intercept (theta1) for the 294",
        "ungenotyped subjects, concluding from CL = 0.35 L/h/kg that",
        "'these patients were primarily EMs/URMs'. The extensive- and",
        "ultrarapid-metabolizer strata are POOLED for patients (theta6",
        "is labelled 'Patient and EM on CL' in Table II but written as",
        "theta6 * (Patient & EM/URM) in the covariate equation, and",
        "Table III reports the resulting 0.41 L/h/kg for 'EM/URM,",
        "Patient'), which is why CYP2D6_UM must be carried separately",
        "rather than folded into CYP2D6_EM: the canonical CYP2D6_EM is",
        "defined as 0 for ultrarapid metabolizers. Mao 2012 fits no",
        "intermediate-metabolizer stratum, so CYP2D6_IM is deliberately",
        "not used; Table I footnote b notes that one subject classified",
        "as EM was subsequently identified as an intermediate",
        "metabolizer. Genotyped cohort: 291 EM (93.6%), 16 PM (5.1%), 4",
        "URM (1.3%) of 311 (Mao 2012 Results, 'Study Population').",
        "Time-fixed per subject (germline genotype-derived phenotype)."
      ),
      source_name        = "EM"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, ultrarapid, or ungenotyped subject)",
      notes              = paste(
        "Single main effect on CL (theta7 = -0.535), applied to poor",
        "metabolizers regardless of subject status -- unlike the",
        "extensive-metabolizer effect, which interacts with",
        "DIS_HEALTHY. Table III accordingly reports one PM clearance,",
        "0.20 L/h/kg, for 'Either' patient status. Poor metabolizers",
        "therefore clear vernakalant about 50% more slowly than patient",
        "extensive metabolizers (0.20 vs 0.41 L/h/kg), which the paper",
        "translates into an approximately 15% higher AUC(0-90 min) and",
        "an approximately 8% higher Cmax -- differences it judges",
        "clinically unimportant. 16 of the 311 genotyped subjects (5.1%)",
        "were PMs, consistent with the expected prevalence in a",
        "predominantly white cohort. Time-fixed per subject (germline",
        "genotype-derived phenotype)."
      ),
      source_name        = "PM"
    ),
    CYP2D6_UM = list(
      description        = "CYP2D6 ultrarapid-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, poor, or ungenotyped subject)",
      notes              = paste(
        "Mao 2012 POOLS ultrarapid with extensive metabolizers for the",
        "patient clearance stratum (theta6 applies to",
        "'Patient & EM/URM'), so in this model CYP2D6_UM enters CL only",
        "through the pooled patient term and carries no separate",
        "coefficient of its own. It must nonetheless be a distinct",
        "column because the canonical CYP2D6_EM is defined as 0 for",
        "ultrarapid metabolizers -- setting CYP2D6_EM = 1 for a URM",
        "subject to obtain the pooled effect would contradict the",
        "register definition. All 4 URM subjects (1.3% of the 311",
        "genotyped) were patients in ACT IV; no healthy volunteer was a",
        "URM, so the paper's theta5 healthy-metabolizer term is defined",
        "for EM only and a healthy URM subject is outside the range of",
        "the data. Time-fixed per subject (germline genotype-derived",
        "phenotype)."
      ),
      source_name        = "URM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 605L,
    n_studies      = 6L,
    n_observations = 2524L,
    age_range      = "study means 60.9-68.3 years in patients; 31.2 years in healthy volunteers",
    age_median     = "60 years (model reference / base case)",
    weight_range   = "study means 75.8-87.4 kg",
    weight_median  = "80 kg (the paper's reference body weight for reporting L/h and L)",
    sex_female_pct = 29.8,
    race_ethnicity = c(White = 96.1, Black = 1.0, Asian = 1.3, Other = 1.5),
    disease_state  = paste(
      "597 patients with atrial fibrillation (90.1%) or atrial flutter",
      "(9.9%) lasting >3 hours to <=45 days, including post-cardiac-surgery",
      "AF/AFL (ACT II); 15.6% had a history of congestive heart failure,",
      "27.5% coronary artery disease, and 50.9% hypertension. Plus 8",
      "healthy male volunteers from the phase 1 study."
    ),
    dose_range     = paste(
      "Patients: vernakalant hydrochloride 3 mg/kg as a 10-minute IV",
      "infusion; if AF/AFL persisted after a 15-minute observation period,",
      "a second 10-minute infusion of 2 mg/kg (i.e. starting 25 minutes",
      "after the start of the first). 201 subjects (33%) received one",
      "infusion and 404 (67%) received two. Healthy volunteers: a single",
      "flat 240 mg 10-minute IV infusion. Doses are labelled as the",
      "HYDROCHLORIDE salt; divide by 1.104 to obtain the free-base amount",
      "this model expects."
    ),
    renal_function = "serum creatinine study means 83.5-105.3 umol/L; 5th-95th percentile 62-133 umol/L",
    cyp2d6_status  = paste(
      "Genotyped in 311 of 605 subjects (51.4%): 291 extensive (93.6%),",
      "16 poor (5.1%), 4 ultrarapid (1.3%). The remaining 294 subjects",
      "(48.6%) were not genotyped and form the model's reference stratum."
    ),
    cohort_split   = "597 patients (98.7%) + 8 healthy volunteers (1.3%) = 605 total",
    regions        = "Not stated; multicenter phase 2/3 program (Scene 2, ACT I-IV) plus a phase 1 study",
    notes          = paste(
      "Baseline demographics from Mao 2012 Table I (six per-study",
      "columns). Of 4518 plasma concentrations from 818 enrolled",
      "subjects, 904 placebo samples, 160 oral-dose samples, 902",
      "otherwise-excluded samples (774 below the limit of quantitation)",
      "and 28 outliers were removed, leaving 2524 observations from 605",
      "subjects. Race percentages are for the 597-subject AF/AFL PATIENT",
      "cohort, summed across the five patient columns of Table I (574",
      "white, 6 black, 8 Asian, 9 other); they exclude the 8 phase 1",
      "volunteers, and the 96.1% white reproduces the figure the paper",
      "states in Results, 'Study Population'. Across all 605 subjects",
      "the split is instead 95.9 / 1.2 / 1.5 / 1.5. Sex percentage is",
      "180 of 605 women, all of them patients; all 8 volunteers were",
      "men, and the patient-only figure the paper quotes is 30.2%.",
      "Estimation was by FOCE in NONMEM VI."
    )
  )

  ini({
    # Structural parameters. Mao 2012 Table II reports the thetas already
    # on the LOG scale (the covariate equations are written as
    # P = exp(theta + ... + eta)), so the values below are transcribed
    # verbatim rather than wrapped in log(). Every parameter is PER
    # KILOGRAM of body weight and is for vernakalant FREE BASE; the
    # multiplication by WT happens in model().
    #
    # Back-transformed reference values (patient, CYP2D6 genotype not
    # collected, age 60 y, serum creatinine 100 umol/L) reproduce Mao 2012
    # Table III exactly: exp(-1.050) = 0.350 L/h/kg, exp(-0.670) = 0.512
    # L/kg, exp(0.038) = 1.039 L/kg, exp(0.291) = 1.338 L/h/kg.
    lcl <- -1.050; label("Clearance intercept, log(L/h/kg); patient, CYP2D6 not genotyped, age 60 y, serum creatinine 100 umol/L")  # Mao 2012 Table II theta1 (CL intercept) = -1.050 +/- 0.036; Table III CL = 0.35 L/h/kg
    lvc <- -0.670; label("Central volume Vc, log(L/kg)")                                                                            # Mao 2012 Table II theta2 (Vc) = -0.670 +/- 0.044; Table III Vc = 0.51 L/kg
    lvp <-  0.038; label("Peripheral volume Vp, log(L/kg)")                                                                         # Mao 2012 Table II theta3 (Vp) = 0.038 +/- 0.027; Table III Vp = 1.04 L/kg
    lq  <-  0.291; label("Intercompartmental clearance Q intercept, log(L/h/kg); patient reference")                                # Mao 2012 Table II theta4 (Q intercept) = 0.291 +/- 0.098; Table III Q = 1.34 L/h/kg in patients

    # Covariate effects on clearance. The CYP2D6 extensive-metabolizer
    # effect is an INTERACTION with subject status: theta5 applies to
    # healthy EM subjects and theta6 to patient EM/UM subjects. The
    # poor-metabolizer effect (theta7) is a main effect applied
    # regardless of subject status. All three are additive shifts on the
    # log scale, relative to the ungenotyped reference stratum.
    e_cyp2d6_em_healthy_cl <- 0.472; label("Effect of healthy volunteer AND CYP2D6 extensive metabolizer on log-CL (unitless)")   # Mao 2012 Table II theta5 (Volunteer and EM on CL) = 0.472 +/- 0.103; Table III CL = 0.56 L/h/kg
    e_cyp2d6_em_patient_cl <- 0.156; label("Effect of patient AND CYP2D6 extensive/ultrarapid metabolizer on log-CL (unitless)")  # Mao 2012 Table II theta6 (Patient and EM on CL) = 0.156 +/- 0.038; Table III CL = 0.41 L/h/kg for EM/URM patients
    e_cyp2d6_pm_cl         <- -0.535; label("Effect of CYP2D6 poor metabolizer on log-CL, either subject status (unitless)")      # Mao 2012 Table II theta7 (PM on CL) = -0.535 +/- 0.061; Table III CL = 0.20 L/h/kg for PMs
    e_age_cl               <- -0.284; label("Power exponent on (AGE / 60) for CL (unitless)")                                     # Mao 2012 Table II theta8 (Age on CL) = -0.284 +/- 0.080; covariate equation theta8 * log(Age / 60)
    e_creat_cl             <- -0.345; label("Power exponent on (CREAT / 100 umol/L) for CL (unitless)")                           # Mao 2012 Table II theta9 (Serum creatinine on CL) = -0.345 +/- 0.083; covariate equation theta9 * log(serum creatinine / 100)

    # Covariate effect on intercompartmental clearance.
    e_dis_healthy_q <- 0.676; label("Effect of healthy volunteer on log-Q (unitless)")  # Mao 2012 Table II theta10 (Volunteer on Q) = 0.676 +/- 0.185; Table III Q = 2.63 L/h/kg in volunteers

    # Interindividual variability: a full 4x4 correlated block on CL, Vc,
    # Vp, and Q. Mao 2012 Table II lists the block COLUMN-WISE (variance
    # then that column's remaining covariances): omega^2 Patient [= CL],
    # CL-Vc, CL-Vp, CL-Q, then Vc, Vc-Vp, Vc-Q, then Vp, Vp-Q ["omega^2
    # p,Q"], then Q. nlmixr2 wants the lower triangle ROW-WISE, which is
    # the reordering below. The diagonal reproduces the Table II
    # "Variability, %" column exactly: sqrt(0.159) = 40%, sqrt(0.446) =
    # 67%, sqrt(0.121) = 35%, sqrt(0.480) = 69%. The implied correlations
    # (0.69 for CL-Vc, 0.84 for Vp-Q) are close to the 0.60 and 0.83
    # reported in Model Assessment, which are correlations of the
    # shrunken empirical-Bayes estimates rather than of the omega block
    # itself. The matrix is positive definite (eigenvalues 0.685, 0.424,
    # 0.086, 0.011).
    etalcl + etalvc + etalvp + etalq ~ c(
      0.159,
      0.184, 0.446,
      0.042, 0.016, 0.121,
      0.048, 0.118, 0.203, 0.480
    )  # Mao 2012 Table II omega^2 rows (Patient/CL, CL-Vc, CL-Vp, Cl-Q, Vc, Vc-Vp, Vc-Q, Vp, p-Q, Q)

    # Residual error: proportional, stratified by subject status. Mao 2012
    # reports sigma^2 as a variance with the corresponding CV in the
    # "Variability, %" column (24% and 28%), so the proportional SD is
    # the square root: sqrt(0.056) = 0.237 and sqrt(0.079) = 0.281.
    propSd_hv <- sqrt(0.056); label("Proportional residual SD for healthy volunteers (fraction)")  # Mao 2012 Table II sigma^2 Volunteer = 0.056 +/- 0.009 (24%)
    propSd_pt <- sqrt(0.079); label("Proportional residual SD for patients (fraction)")            # Mao 2012 Table II sigma^2 Patient = 0.079 +/- 0.009 (28%)
  })

  model({
    # 1. Derived covariate terms.
    #
    # The CYP2D6 effect on clearance is stratified by subject status, and
    # the four strata Mao 2012 fits are mutually exclusive:
    #   * all indicators 0        -> genotype not collected (the reference)
    #   * healthy AND EM          -> theta5
    #   * patient AND (EM or UM)  -> theta6
    #   * PM, either status       -> theta7
    # CYP2D6_EM and CYP2D6_UM are mutually exclusive by construction, so
    # their sum is itself a valid 0/1 indicator for the pooled EM/UM
    # stratum the paper's theta6 applies to.
    cyp2d6_em_healthy <- DIS_HEALTHY * CYP2D6_EM
    cyp2d6_em_patient <- (1 - DIS_HEALTHY) * (CYP2D6_EM + CYP2D6_UM)

    # 2. Individual parameters. Mao 2012's covariate equations, with the
    # per-kilogram estimates scaled to absolute units by body weight:
    #   CL = exp(theta1 + theta5*(Volunteer & EM) + theta6*(Patient & EM/URM)
    #            + theta7*(PM) + theta8*log(Age/60)
    #            + theta9*log(serum creatinine/100) + eta_CL)
    #   Vc = exp(theta2 + eta_Vc)
    #   Vp = exp(theta3 + eta_Vp)
    #   Q  = exp(theta4 + theta10*(Subject = Volunteer) + eta_Q)
    cl <- exp(lcl + etalcl +
                e_cyp2d6_em_healthy_cl * cyp2d6_em_healthy +
                e_cyp2d6_em_patient_cl * cyp2d6_em_patient +
                e_cyp2d6_pm_cl * CYP2D6_PM +
                e_age_cl * log(AGE / 60) +
                e_creat_cl * log(CREAT / 100)) * WT
    vc <- exp(lvc + etalvc) * WT
    vp <- exp(lvp + etalvp) * WT
    q  <- exp(lq + etalq + e_dis_healthy_q * DIS_HEALTHY) * WT

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system: open two-compartment mammillary disposition with
    # first-order elimination from the central compartment (Mao 2012
    # Figure 1).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation. central is mg of vernakalant free base and vc is L,
    # so central/vc is mg/L = ug/mL; the factor of 1000 converts to the
    # ng/mL used throughout Mao 2012 (Table IV Cmax 3371-4558 ng/mL).
    Cc <- 1000 * central / vc

    # 6. Subject-status-stratified proportional residual error.
    propSdEff <- propSd_hv * DIS_HEALTHY + propSd_pt * (1 - DIS_HEALTHY)
    Cc ~ prop(propSdEff)
  })
}
