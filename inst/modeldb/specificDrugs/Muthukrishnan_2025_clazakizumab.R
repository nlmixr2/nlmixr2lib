Muthukrishnan_2025_clazakizumab <- function() {
  description <- "Population PK and PK-PD model for the anti-interleukin-6 monoclonal antibody clazakizumab given as a 3-minute IV bolus to adults with end-stage kidney disease undergoing maintenance dialysis (POSIBIL 6 ESKD phase 2b, NCT05485961; Muthukrishnan 2025). Two-compartment linear PK with allometric weight scaling on CL and V1 and a baseline free IL-6 covariate on CL; sequential-fit indirect-response inhibition of hs-CRP zero-order production (kin) with Imax fixed to 1, an estimated Hill coefficient, and a Manly-transformed IIV on IC50. Baseline hs-CRP enters kin as a power covariate so the pre-dose steady state E0 = kin/kout tracks the observed baseline distribution."
  reference <- "Muthukrishnan VY, Kerbusch T, Strong LE, Kleijn HJ, Pfister M, Chang AM, Acharya M, Nandy P, McCune JS. Population Pharmacokinetic and Pharmacokinetic-Pharmacodynamic Analysis for Clazakizumab in Patients With End-Stage Kidney Disease Undergoing Dialysis. Clin Transl Sci. 2025. doi:10.1111/cts.70381"
  vignette <- "Muthukrishnan_2025_clazakizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL", hs_CRP = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline weight used for power scaling on CL and V1 with reference 95 kg (PK-dataset overall median; Table S1, N=95). Muthukrishnan 2025 states the weight allometric exponents were estimated freely (0.979 on CL, 0.784 on V1) rather than fixed to 0.75/1; the paper notes 'Any effects by sex are likely already described through the weight-based allometry model' (Results 3.1). The same median 95 kg holds in the larger PD dataset (Table S2, N=126).",
      source_name        = "WT"
    ),
    IL6 = list(
      description        = "Baseline free serum interleukin-6 concentration",
      units              = "ng/L (equivalently pg/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline free IL-6 used for a power effect on CL with reference 6.80 ng/L (PK-dataset overall median; Table S1). Retained after producing a 20-point OFV drop (p < 0.001; Results 3.1). Missing in 3/95 subjects (3.2%). Consistent with clazakizumab's mechanism (anti-IL-6 mAb): higher baseline free IL-6 modestly increases CL.",
      source_name        = "IL6"
    ),
    CRP = list(
      description        = "Baseline high-sensitivity C-reactive protein (predose Day 1)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline hs-CRP (Roche Cardiac CRP Latex High Sensitive assay on Cobas c502; particle-enhanced immunoturbidimetry). Value = predose Day 1 sample per Muthukrishnan 2025 Methods 2.4; the paper's 'Baseline hs-CRP: Predose Day 1 (mg/L)' median = 8.15 mg/L (PD-dataset overall, Table S2, N=126) is used as the reference for the power scaling on kin. Study inclusion required baseline hs-CRP >= 2 mg/L (assay linearity 0.15-20.0 mg/L extended to 0.2-300 mg/L with 1:15 dilution). Because the exponent 0.809 < 1 and E0 = kin/kout at steady state, the model-implied E0 partially compensates but does not exactly reproduce each subject's observed baseline; the PD residual absorbs the offset.",
      source_name        = "CRP"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a baseline covariate on CL and V1 (Methods 2.3) but not retained in the final popPK model. Also screened for the PK-PD model (on kin or IC50; Methods 2.4) but not retained. Table S1 PK-dataset median 65 years [31-83]."
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened on CL, V1, and V2 but not retained (p > 0.001; Results 3.1). Paper attributes any residual sex signal to the weight-based allometry that is already in the model. PK dataset 66.3% male; PD dataset 66.7% male."
    ),
    RACE = list(
      description = "Self-reported race (White / Black or African American / Asian / American Indian or Alaska Native / Native Hawaiian or Other Pacific Islander / Other / Multiple)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on CL and V1 but not retained (Methods 2.3, Table S1). PK-dataset composition per Table S1: White 57.9%, Black or African American 35.8%, Asian 1.1%, American Indian/Alaska Native 1.1%, Native Hawaiian/Pacific Islander 1.1%, Other 1.1%, Multiple 2.1%."
    ),
    AST = list(
      description = "Aspartate aminotransferase (baseline)",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a baseline liver-enzyme covariate on CL (Methods 2.3) but not retained. Table S1 PK-dataset median 13.0 U/L [6.00-33.0]."
    ),
    ALT = list(
      description = "Alanine aminotransferase (baseline)",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a baseline liver-enzyme covariate on CL (Methods 2.3) but not retained. Table S1 PK-dataset median 12.0 U/L [4.00-51.0]; missing 2.1%."
    ),
    ALP = list(
      description = "Alkaline phosphatase (baseline)",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a baseline liver-enzyme covariate on CL (Methods 2.3) but not retained. Table S1 PK-dataset median 83.0 U/L [39.0-561]."
    ),
    BILI = list(
      description = "Total bilirubin (baseline)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as a baseline hepatic covariate on CL (Methods 2.3) but not retained. Table S1 PK-dataset median 0.450 mg/dL [0.240-0.959]."
    ),
    ALB = list(
      description = "Serum albumin (baseline)",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened on CL (PK) and on kin/IC50 (PK-PD) but not retained. Table S1/S2 median 38 g/L [33-47]. Note the population is on dialysis so albumin is skewed low; had it been retained it would likely have entered as a power effect analogous to Li 2018 PF-04236921."
    )
  )

  population <- list(
    species        = "human",
    n_subjects_pk  = 95L,
    n_subjects_pd  = 126L,
    n_subjects     = 126L,
    n_studies      = 1L,
    trial          = "POSIBIL 6 ESKD phase 2b (NCT05485961)",
    age_range      = "31-83 years (PK dataset); 31-86 years (PD dataset with placebo)",
    age_median     = "65 years (PK dataset); 66 years (PD dataset)",
    weight_range   = "53-163 kg",
    weight_median  = "95 kg",
    sex_female_pct = 33.7,
    race_ethnicity = c(
      White                                   = 57.9,
      `Black or African American`             = 35.8,
      Asian                                   = 1.1,
      `American Indian / Alaska Native`       = 1.1,
      `Native Hawaiian / Pacific Islander`    = 1.1,
      Other                                   = 1.1,
      Multiple                                = 2.1
    ),
    disease_state  = "End-stage kidney disease receiving maintenance dialysis; enrolment required cardiovascular disease and/or diabetes plus baseline hs-CRP >= 2 mg/L. Dosing occurred during dialysis; all patients had chronic inflammation (median baseline hs-CRP 8.15 mg/L, range 0.7-215 mg/L).",
    dose_range     = "Clazakizumab 2.5, 5, or 10 mg IV every 4 weeks (Q4W) as a 3-minute IV bolus via the return venous line of the hemodialysis circuit; placebo control; up to 6 doses (minimum treatment 12 weeks / 3 doses).",
    regions        = "Multi-regional phase 2b (POSIBIL 6 ESKD).",
    n_observations = c(pk = 891L, pd = 1401L),
    notes          = paste(
      "PK analysis dataset: 891 measurable clazakizumab concentrations from 95 subjects (Results 3.1).",
      "PD analysis dataset: 1401 measurable hs-CRP observations from 126 subjects (Results 3.2); 8.9% of post-first-dose hs-CRP observations were BLQ and handled by M3 imputation.",
      "Assay ranges: clazakizumab enzyme immunoassay 20.0-1250 ng/mL (measures free / unbound antibody); hs-CRP Roche Cardiac CRP (Latex) High Sensitive 0.15-20.0 mg/L extended to 0.2-300 mg/L with automatic 1:15 dilution.",
      "Sequential 2-stage PK-PD fit: individual post-hoc clazakizumab concentrations from the popPK model drove the PD fit (Methods 2.4).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # PK structural parameters. Muthukrishnan 2025 Table 1 final serum
    # clazakizumab popPK estimates. Reference covariate values are the
    # PK-dataset overall medians from Table S1: WT = 95 kg, IL6 = 6.80 ng/L.
    # Doses in mg; central / peripheral state amounts in mg; CL and Q in
    # L/day; V1 and V2 in L. Observed concentration Cc in ng/mL after unit
    # conversion in model() (1 mg/L = 1000 ng/mL).
    # =====================================================================
    lcl <- log(0.228); label("Clearance CL (L/day) at reference WT=95 kg, IL6=6.80 ng/L")  # Muthukrishnan 2025 Table 1
    lvc <- log(4.05);  label("Central volume of distribution V1 (L) at reference WT=95 kg")  # Muthukrishnan 2025 Table 1
    lq  <- log(0.443); label("Intercompartmental clearance Q (L/day)")                       # Muthukrishnan 2025 Table 1
    lvp <- log(3.51);  label("Peripheral volume of distribution V2 (L)")                     # Muthukrishnan 2025 Table 1

    # Body-weight allometric exponents. Estimated freely (not fixed to
    # 0.75/1); the Discussion notes weight-based allometry likely accounts
    # for the absence of a retained sex effect.
    e_wt_cl <- 0.979; label("Power exponent of WT/95 on CL (unitless)")   # Muthukrishnan 2025 Table 1 (WT~CL)
    e_wt_vc <- 0.784; label("Power exponent of WT/95 on V1 (unitless)")   # Muthukrishnan 2025 Table 1 (WT~V1)

    # Baseline free IL-6 covariate effect on CL.
    e_il6_cl <- 0.19; label("Power exponent of IL6/6.80 on CL (unitless)")  # Muthukrishnan 2025 Table 1 (Free IL-6~CL)

    # ---------------------------------------------------------------------
    # PK inter-individual variability. Muthukrishnan 2025 Table 1 reports
    # CV% back-transformed from the NONMEM variance-scale omega via
    # CV = sqrt(exp(variance) - 1) per the Table 1 footnote. Reversing that
    # gives the underlying variance = ln(1 + CV^2):
    #   CL  CV 48.7% -> omega^2 = ln(1 + 0.487^2) = 0.21290
    #   V1  CV 36.3% -> omega^2 = ln(1 + 0.363^2) = 0.12378
    #   V2  CV 72.4% -> omega^2 = ln(1 + 0.724^2) = 0.42156
    # Off-diagonal is reported directly on the covariance scale as
    # 'V1~CL 0.152' (%RSE 27.9%, 95% CI 0.0687-0.235; the CI arithmetic
    # 0.152 +/- 1.96*0.152*0.279 reproduces the reported CI, confirming
    # 0.152 is the covariance and not the correlation). Implied correlation
    # 0.152 / sqrt(0.21290 * 0.12378) = 0.936, high but expected for the
    # jointly weight-scaled CL and V1.
    # Q has no eta (Results 3.1: "IIV was ... not able to be reliably
    # estimated on intercompartmental clearance (Q)").
    # ---------------------------------------------------------------------
    etalcl + etalvc ~ c(0.21290,
                        0.152, 0.12378)  # Muthukrishnan 2025 Table 1 (CL, V1) IIV block
    etalvp ~ 0.42156                     # Muthukrishnan 2025 Table 1 CV(eta_V2)

    # PK residual error -- combined additive + proportional on linear ng/mL
    # (Results 3.1: "Proportional and combined additive + proportional
    # error models were evaluated, and the combined error model was
    # selected"; Table 1).
    addSd  <- 23.1;  label("Additive residual error (ng/mL)")             # Muthukrishnan 2025 Table 1 Additive
    propSd <- 0.266; label("Proportional residual error (fraction)")      # Muthukrishnan 2025 Table 1 Proportional 26.6%

    # =====================================================================
    # PK-PD structural parameters. Muthukrishnan 2025 Table 2 final PK-PD
    # model on ln-transformed hs-CRP. Indirect-response with drug
    # inhibiting the zero-order production rate (kin):
    #   d/dt(effect) = kin * (1 - Imax * Cc^gamma / (IC50^gamma + Cc^gamma))
    #                - kout * effect
    #   effect(0)   = kin / kout   (steady-state pre-dose, Methods 2.4)
    # Imax was fixed to 1 as part of the paper's model-refinement strategy
    # (Study Highlights, Discussion). Sequential 2-stage fit using post-hoc
    # PK from Table 1.
    # =====================================================================
    lkout <- log(0.381); label("First-order hs-CRP elimination rate kout (1/day)")                          # Muthukrishnan 2025 Table 2
    lkin  <- log(3.76);  label("Zero-order hs-CRP production rate kin ([mg/L]/day) at reference CRP=8.15")  # Muthukrishnan 2025 Table 2
    lic50 <- log(3.39);  label("Clazakizumab IC50 for kin inhibition (ng/mL)")                              # Muthukrishnan 2025 Table 2
    limax <- fixed(log(1)); label("Maximum fractional inhibition of kin (Imax)")                # Muthukrishnan 2025 Table 2 "1.00 Fixed"
    lhill <- log(0.523); label("Hill coefficient on the clazakizumab-vs-kin sigmoid (unitless)")            # Muthukrishnan 2025 Table 2 Hill factor

    # Baseline hs-CRP covariate effect on kin.
    e_crp_kin <- 0.809; label("Power exponent of CRP/8.15 on kin (unitless)")  # Muthukrishnan 2025 Table 2 (Baseline hs-CRP~kin)

    # ---------------------------------------------------------------------
    # PD inter-individual variability. Only IC50 has an eta and it uses a
    # Manly transformation (Muthukrishnan 2025 Study Highlights / Methods
    # 2.4 / Results 3.2) to normalise a heavy-tailed raw eta distribution:
    #   ETA_Manly = (exp(lambda * eta) - 1) / lambda   (with lambda != 0)
    #   IC50_i    = IC50 * exp(ETA_Manly)
    # The shape parameter lambda is estimated separately from the eta
    # variance itself.
    # CV(eta_IC50) = 1002% is a linear-space CV back-transformed the same
    # way as the PK IIVs: omega^2 = ln(1 + 10.02^2) = 4.61932.
    # No IIV on kout, kin, Imax, or the Hill coefficient.
    # ---------------------------------------------------------------------
    manly_ic50 <- -0.143; label("Manly-transformation shape parameter lambda on eta_IC50 (unitless)")  # Muthukrishnan 2025 Table 2 (Shape parameter on IC50 IIV)
    etalic50 ~ 4.61932  # Muthukrishnan 2025 Table 2 CV(eta_IC50) = 1002% -> ln(1 + 10.02^2)

    # PD residual: "Residual error was estimated as additive error on the
    # log scale" (Table 2 note); Table 2 reports it in the CV% column as
    # 72.6% under the same sqrt(exp(sigma^2)-1) presentation used for the
    # IIVs. Back-transforming: sigma^2 = ln(1 + 0.726^2) = 0.42355 so
    # sigma_log = 0.6508. That log-scale SD is what nlmixr2 lnorm() takes.
    expSd_hsCRP <- 0.6508; label("Log-normal (log-scale additive) residual SD on hs-CRP")  # Muthukrishnan 2025 Table 2 Proportional (%) 72.6 back-transformed via ln(1+CV^2)
  })

  model({
    # =====================================================================
    # Individual PK parameters. Muthukrishnan 2025 Methods 2.3 continuous-
    # covariate form P_ki = theta_k * (X_i / M(X))^theta_j. Reference
    # covariate values are the PK-dataset medians (Table S1): WT = 95 kg,
    # IL6 = 6.80 ng/L.
    # =====================================================================
    cl <- exp(lcl + etalcl) * (WT / 95)^e_wt_cl * (IL6 / 6.80)^e_il6_cl
    vc <- exp(lvc + etalvc) * (WT / 95)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # =====================================================================
    # 2-compartment PK with IV administration to the central compartment.
    # Doses are given as a 3-minute IV bolus into the hemodialysis return
    # venous line (Methods 2.1); modelled as instantaneous input to
    # `central` via the dose record (cmt = "central" in event tables).
    # =====================================================================
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Cc in ng/mL. Doses in mg; central amount in mg; vc in L so central/vc
    # is in mg/L. 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # =====================================================================
    # Individual PD parameters.
    #   kin scales with observed baseline hs-CRP via a power effect.
    #   IC50 uses a Manly-transformed eta (paper Methods 2.4 / Results 3.2).
    # =====================================================================
    kout <- exp(lkout)
    kin  <- exp(lkin) * (CRP / 8.15)^e_crp_kin
    imax <- exp(limax)
    hill <- exp(lhill)

    eta_manly_ic50 <- (exp(manly_ic50 * etalic50) - 1) / manly_ic50
    ic50 <- exp(lic50) * exp(eta_manly_ic50)

    # =====================================================================
    # Indirect-response hs-CRP ODE (paper Equation for E in Results 3.2)
    # with steady-state pre-dose initial condition E(0) = kin/kout.
    # =====================================================================
    inh <- imax * Cc^hill / (ic50^hill + Cc^hill)

    d/dt(effect) <- kin * (1 - inh) - kout * effect
    effect(0)    <- kin / kout

    hsCRP <- effect

    # PK residual: combined additive + proportional on linear ng/mL.
    Cc    ~ add(addSd) + prop(propSd)
    # PD residual: additive on log scale (equivalently, lognormal on the
    # linear-space observation).
    hsCRP ~ lnorm(expSd_hsCRP)
  })
}
