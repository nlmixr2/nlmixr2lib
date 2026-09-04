Bamgboye_2026_topiramate <- function() {
  description <- "Three-compartment population pharmacokinetic model with linear elimination for a novel intravenous topiramate (TPM) formulation, developed from 246 plasma concentrations in 20 adult patients with epilepsy or migraine who received a single 25 mg dose of stable-isotope-labeled IV TPM infused over 10 minutes on top of their maintenance oral TPM therapy (Bamgboye 2026). Clearances (CL, Q2, Q3) scale allometrically with body weight with a fixed exponent of 0.75 and volumes (V1, V2, V3) with a fixed exponent of 1, both referenced to 70 kg. Concomitant enzyme-inducing antiepileptic comedication (carbamazepine, phenytoin, or oxcarbazepine) is the only retained covariate and raises central clearance by 63% in power form (1.63^CONMED_EIAED); age, height, sex, and creatinine clearance were screened and not retained. Inter-individual variability is estimated on CL, V1, and V2 as a correlated block, and residual unexplained variability is proportional (12.8%CV). This is the first population PK model of intravenous topiramate in patients, and the analysis supported the conclusion that IV TPM loading doses do not require adjustment for enzyme-inducing comedications because those comedications affect clearance but not volume of distribution."
  reference <- paste(
    "Bamgboye AO, Coles LD, Suriyapakorn B, Mishra U, Kriel RL, Leppik IE,",
    "White JR, Cloyd JC.",
    "Population pharmacokinetic modeling of intravenous topiramate in patients",
    "with epilepsy or migraine.",
    "J Clin Pharmacol. 2026;66(4):e70191. doi:10.1002/jcph.70191.",
    sep = " "
  )
  vignette <- "Bamgboye_2026_topiramate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The model is intravenous-only (no depot); doses land
  # directly in `central` as a 10-minute infusion.
  compartmentData <- list(
    central     = list(analyte = "topiramate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "topiramate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "topiramate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling referenced to 70 kg, applied to every disposition parameter: exponent fixed at 0.75 for the clearances (CL, Q2, Q3) and at 1 for the volumes (V1, V2, V3). Bamgboye 2026 Methods states the exponents were fixed a priori 'to balance model complexity with physiological plausibility'; Table 2 prints the scaling inline in each parameter row (e.g. 'CL (L/h) 1.31 x (WT/70)^0.75'). Weight is time-fixed in the source analysis (single-dose study).",
      source_name        = "WT"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiepileptic drug indicator; 1 = the patient was receiving an enzyme-inducing comedication, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no enzyme-inducing comedication)",
      notes              = "The 7 of 20 patients (35%) coded 1 were taking carbamazepine (n = 4), phenytoin (n = 1), oxcarbazepine (n = 1), or both carbamazepine and phenytoin (n = 1) (Bamgboye 2026 Methods, Study Population). The source deliberately collapses all inducing comedications into a single binary covariate rather than per-drug terms because of the small sample size (Bamgboye 2026 Discussion, limitations). Effect is multiplicative on CL in power form: CL is multiplied by 1.63^CONMED_EIAED, i.e. 63% higher clearance on an inducer. Note this cohort's inducer list includes oxcarbazepine, which is a weaker inducer than carbamazepine or phenytoin; the source acknowledges the uniform-effect assumption as a limitation.",
      source_name        = "Inducer"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years; screened as a covariate on CL and V1 in the forward-inclusion step and not retained (dOFV = -0.95, below the 3.84 threshold).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not retained (Bamgboye 2026 Table S3 model 3; Discussion). Cohort age 26-74 years, mean 39.8."
    ),
    HT = list(
      description = "Body height at baseline; screened as a covariate on CL and V1 in the forward-inclusion step and not retained (dOFV = -1.53, below the 3.84 threshold).",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not retained (Bamgboye 2026 Table S3 model 4). Cohort height 155-182 cm, mean 169."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation; screened on CL in addition to the retained inducer effect and not retained (dOFV = -1.84, below the 3.84 threshold).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened; not retained (Bamgboye 2026 Table S3 model 6; Discussion). The source attributes the absent effect to the near-normal renal function of the cohort (46-206 mL/min, median 102.5) and states the model may not generalize to renal impairment.",
      source_name = "CrCL"
    ),
    SEXF = list(
      description = "Female-sex indicator; listed among the tested covariates in Methods and not retained in the final model.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained (Bamgboye 2026 Methods, covariate analysis; Discussion). Cohort 13 female (65%) / 7 male (35%). Sex does not appear as a numbered row in Table S3, so no dOFV is reported for it."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = 246L,
    age_range      = "26-74 years",
    age_median     = "41.5 years",
    weight_range   = "54.5-150.3 kg",
    weight_median  = "85.2 kg",
    height_range   = "155-182 cm",
    height_median  = "170 cm",
    sex_female_pct = 65,
    race_ethnicity = c(White = 100),
    disease_state  = "Adults on maintenance oral topiramate therapy for epilepsy or migraine. Pregnant or breastfeeding individuals, and those with a history of intolerance to intravenous administration or known topiramate hypersensitivity, were excluded.",
    dose_range     = "Single 25 mg dose of stable-isotope-labeled intravenous topiramate (10 mg/mL in 10% sulfobutyl ether beta-cyclodextrin), infused over 10 minutes, given alongside the patient's usual morning oral topiramate dose.",
    renal_function = "Creatinine clearance (Cockcroft-Gault) 46-206 mL/min, median 102.5 (mean 106.7, SD 34.69). No patients with clinically important renal impairment.",
    co_medication  = "7 of 20 patients (35%) on enzyme-inducing comedications: carbamazepine (n = 4), phenytoin (n = 1), oxcarbazepine (n = 1), carbamazepine + phenytoin (n = 1). All patients were on maintenance oral topiramate, which the stable-isotope label allows to be quantified separately from the IV dose (distinct m/z 338 vs 344 transitions), so no oral washout was required.",
    regions        = "Single centre, United States (University of Minnesota; IND #78993).",
    notes          = "Demographics from Bamgboye 2026 Table 1. Rich sampling: predose and 5, 15, 30 min and 1, 2, 4, 6, 12, 24, 48, 72, 96 h post-dose. The source attributes its three-compartment structure (versus the one- or two-compartment models typical of oral topiramate analyses) to this early-time-point-dense sampling; see modellib('Ahmed_2015_topiramate') for the two-compartment IV + oral model in healthy volunteers that this paper cites as reference 21."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural disposition parameters - final model estimates from
    # Bamgboye 2026 Table 2, "Fixed effects estimate" column, reported
    # as typical values for a 70 kg person. Source-to-canonical name
    # map: V1 -> `lvc` (central volume); V2 -> `lvp` (first peripheral
    # volume); V3 -> `lvp2` (second peripheral volume); Q2 -> `lq`
    # (central <-> peripheral1 intercompartmental clearance);
    # Q3 -> `lq2` (central <-> peripheral2). Table 2 reports linear-scale
    # point estimates, so they are log-transformed here.
    # ------------------------------------------------------------------
    lcl  <- log(1.31)  ; label("Central clearance CL at 70 kg (L/h)")                          # Bamgboye 2026 Table 2: CL = 1.31 x (WT/70)^0.75, RSE 9.51%, bootstrap 95% CI 1.01-1.53
    lvc  <- log(9.84)  ; label("Central volume of distribution V1 at 70 kg (L)")               # Bamgboye 2026 Table 2: V1 = 9.84 x (WT/70), RSE 8.7%, bootstrap 95% CI 8.49-11.0 (see vignette Errata: the Results-section equation misprints this as 15.6)
    lq   <- log(197)   ; label("Intercompartmental clearance Q2 (central <-> peripheral1) at 70 kg (L/h)") # Bamgboye 2026 Table 2: Q2 = 197 x (WT/70)^0.75, RSE 6.73%, bootstrap 95% CI 181-223
    lvp  <- log(39.1)  ; label("First peripheral volume of distribution V2 at 70 kg (L)")      # Bamgboye 2026 Table 2: V2 = 39.1 x (WT/70), RSE 5.1%, bootstrap 95% CI 36.5-41.8
    lq2  <- log(0.6)   ; label("Intercompartmental clearance Q3 (central <-> peripheral2) at 70 kg (L/h)") # Bamgboye 2026 Table 2: Q3 = 0.6 x (WT/70)^0.75, RSE 41.7%, bootstrap 95% CI 0.4-1.2
    lvp2 <- log(9.01)  ; label("Second peripheral volume of distribution V3 at 70 kg (L)")     # Bamgboye 2026 Table 2: V3 = 9.01 x (WT/70), RSE 18.8%, bootstrap 95% CI 6.41-44.3

    # ------------------------------------------------------------------
    # Allometric exponents. Both were FIXED a priori, not estimated:
    # Bamgboye 2026 Methods states "a fixed allometric exponent of 0.75
    # for clearances and 1 for volumes". A single shared exponent is used
    # across all three clearances (CL, Q2, Q3) and all three volumes
    # (V1, V2, V3), so the shared-exponent naming form
    # `e_<cov>_<param1>_<param2>` applies (parameter-names.md).
    # ------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent of body weight on all clearances CL, Q2, Q3 (unitless)") # Bamgboye 2026 Methods, PopPK model development; printed inline in every clearance row of Table 2
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent of body weight on all volumes V1, V2, V3 (unitless)")    # Bamgboye 2026 Methods, PopPK model development; printed inline in every volume row of Table 2

    # ------------------------------------------------------------------
    # Covariate effect: enzyme-inducing comedication on CL, in power
    # form on the binary indicator. Bamgboye 2026 Results prints the
    # relationship explicitly:
    #   CL (L/h) = 1.31 x (WT/70)^0.75 x 1.63^Inducer
    # so the coefficient is a multiplicative ratio (63% higher CL on an
    # inducer), NOT a log-scale or fractional-change coefficient. Kept in
    # the source's power form per parameter-names.md (power-style effect).
    # ------------------------------------------------------------------
    e_eiaed_cl <- 1.63 ; label("Fold-change in CL with enzyme-inducing comedication (ratio, power form)") # Bamgboye 2026 Table 2 "Inducer ~ CL" = 1.63, RSE 16.3%, bootstrap 95% CI 1.12-2.30; Results equation CL = 1.31 x (WT/70)^0.75 x 1.63^Inducer

    # ------------------------------------------------------------------
    # Inter-individual variability. Exponential (log-normal) IIV on CL,
    # V1, and V2 only - the source states "an exponential IIV model for
    # central clearance and volume" and Table 2 lists IIV rows for CL,
    # V1, and V2 but none for Q2, Q3, or V3.
    #
    # Table 2 reports IIV as %CV with the footnote %CV = sqrt(e^omega^2 - 1),
    # so the internal log-scale variance is omega^2 = log(1 + CV^2):
    #   CL: 33.5%CV -> log(1 + 0.335^2) = 0.1063625
    #   V1: 18.1%CV -> log(1 + 0.181^2) = 0.0322358
    #   V2: 19.3%CV -> log(1 + 0.193^2) = 0.0365720
    #
    # Table 2's three "Correlation (Omega X,Y)" rows are correlation
    # COEFFICIENTS, not covariances; the off-diagonal covariance is
    # therefore r * omega_X * omega_Y. Reading them as covariances is
    # arithmetically impossible: 0.06 as cov(CL,V1) implies
    # r = 0.06 / sqrt(0.1063625 * 0.0322358) = 1.025 > 1, and the
    # resulting matrix has a negative eigenvalue (-0.0056). The
    # correlation reading (used here) is positive definite.
    #   cov(CL,V1) = 0.06 * sqrt(0.1063625 * 0.0322358) = 0.00351330
    #   cov(CL,V2) = 0.03 * sqrt(0.1063625 * 0.0365720) = 0.00187107
    #   cov(V1,V2) = 0.03 * sqrt(0.0322358 * 0.0365720) = 0.00103007
    # ------------------------------------------------------------------
    # Element-by-element source trace for the block below (comments cannot be
    # placed inside the c(...) itself: rxode2's ini() comment-to-label rewriter
    # turns a trailing comment on an element line into a label() call and the
    # block then fails to parse):
    #   [1,1] 0.1063625  Table 2 "IIV on CL (CV%)"  33.5 (RSE 40%,   bootstrap 20.2-52.6, shrinkage 1.48%)
    #   [2,1] 0.00351330 Table 2 "Correlation (Omega CL,V1)" 0.06 (RSE 51.7%)
    #   [2,2] 0.0322358  Table 2 "IIV on V1 (CV%)"  18.1 (RSE 74.2%, bootstrap 5.5-32.0,  shrinkage 1.50%)
    #   [3,1] 0.00187107 Table 2 "Correlation (Omega CL,V2)" 0.03 (RSE 60.4%)
    #   [3,2] 0.00103007 Table 2 "Correlation (Omega V1,V2)" 0.03 (RSE 52.8%)
    #   [3,3] 0.0365720  Table 2 "IIV on V2 (CV%)"  19.3 (RSE 37.6%, bootstrap 11.9-24.1, shrinkage 4.21%)
    etalcl + etalvc + etalvp ~ c(
      0.1063625,
      0.00351330, 0.0322358,
      0.00187107, 0.00103007, 0.0365720
    )

    # ------------------------------------------------------------------
    # Residual unexplained variability: proportional only (the source
    # tested additive, proportional, and combined, and selected
    # proportional by likelihood-ratio test). Table 2 prints the row as
    # "Proportional (%CV)  0.02(12.8)": the variance rounded to two
    # decimals, with the %CV in parentheses. The %CV carries three
    # significant figures and is the more precise of the two printed
    # numbers (0.128^2 = 0.0164, which rounds to the printed 0.02), so
    # the SD passed to prop() is taken from it.
    # ------------------------------------------------------------------
    propSd <- 0.128 ; label("Proportional residual error SD (fraction)") # Bamgboye 2026 Table 2: Proportional 0.02 (12.8 %CV), RSE 13.7%, bootstrap 95% CI 0.01-0.02, epsilon-shrinkage 6.28% (Table S2)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual disposition parameters. Allometric weight scaling on
    # every parameter, referenced to 70 kg; the enzyme-inducer effect
    # applies to CL only. IIV is a correlated log-normal block on CL,
    # V1, and V2; Q2, Q3, and V3 have no IIV in the source model.
    #   CL (L/h) = 1.31 x (WT/70)^0.75 x 1.63^Inducer   [Results equation]
    # ------------------------------------------------------------------
    cl  <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl_q * e_eiaed_cl ^ CONMED_EIAED
    vc  <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc_vp
    q   <- exp(lq)           * (WT / 70) ^ e_wt_cl_q
    vp  <- exp(lvp + etalvp) * (WT / 70) ^ e_wt_vc_vp
    q2  <- exp(lq2)          * (WT / 70) ^ e_wt_cl_q
    vp2 <- exp(lvp2)         * (WT / 70) ^ e_wt_vc_vp

    # ------------------------------------------------------------------
    # 2. Three-compartment ODE system with linear elimination from the
    # central compartment (NONMEM ADVAN11/TRANS4 equivalent). The study
    # drug is given intravenously only, so there is no depot state and
    # no bioavailability term; doses enter `central` directly (as a
    # 10-minute infusion in the source study).
    # ------------------------------------------------------------------
    d/dt(central)     <- -(cl / vc) * central -
                          (q  / vc) * central + (q  / vp ) * peripheral1 -
                          (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <-  (q  / vc) * central - (q  / vp ) * peripheral1
    d/dt(peripheral2) <-  (q2 / vc) * central - (q2 / vp2) * peripheral2

    # ------------------------------------------------------------------
    # 3. Observation. Doses are in mg and volumes in L, so central/vc is
    # in mg/L; multiply by 1000 to report ng/mL, the units used for every
    # concentration axis in the source (Figures 1-4) and by the LC-MS
    # assay (Bamgboye 2026 Methods).
    # ------------------------------------------------------------------
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
