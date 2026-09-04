Yin_2026_ceftazidime <- function() {
  description <- "One-compartment IV population PK model for ceftazidime in Chinese neonates sampled by quantitative dried blood spot (qDBS) micro-sampling (Yin 2026). Concentrations are WHOLE BLOOD, not plasma. Body weight enters a priori as allometric scaling on CL (exponent fixed at 0.75) and V (exponent fixed at 1) referenced to 70 kg, and postmenstrual age drives a Rhodin-type sigmoidal glomerular-filtration maturation function on CL (TM50 fixed at 47.7 weeks, Hill coefficient fixed at 3.4). No covariate beyond WT and PMA was retained by stepwise selection."
  reference <- paste(
    "Yin Q, Chen Y, Lv J, Zheng Z, Zhao X, Chen H, Xie F, Yi H, Chen Q.",
    "Population pharmacokinetic modeling and simulation-informed",
    "ceftazidime dosing in Chinese neonates using quantitative dried blood",
    "spot micro-sampling.",
    "Antimicrob Agents Chemother. 2026;70(6):e01810-25.",
    "doi:10.1128/aac.01810-25.",
    sep = " "
  )
  vignette <- "Yin_2026_ceftazidime"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # The model was fitted directly to capillary WHOLE BLOOD concentrations
    # measured on Capitainer qDBS cards (Yin 2026 "Sample collection and
    # determination" and Discussion: "Whole-blood concentrations obtained
    # from qDBS samples were used directly to develop the ceftazidime PopPK
    # model"). CL and V are therefore blood-referenced, not plasma-referenced.
    central = list(analyte = "ceftazidime", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered a priori (not by stepwise selection) as allometric scaling",
        "on both CL and V referenced to 70 kg, per Yin 2026 equations 1 and",
        "2 and Methods 'PopPK modeling': 'Clearance (CL) and volume of",
        "distribution (V) were scaled allometrically based on WT, with",
        "exponents of 0.75 for CL and 1 for V'. Both exponents were fixed,",
        "not estimated. Cohort median 3.2 kg (range 1.8-4.2 kg; Table 1).",
        "Treated as time-fixed at study entry -- the study sampled within a",
        "few days of birth (median postnatal age 1 day) so body weight does",
        "not change materially over the sampling window."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WEEKS, not the register-default months: Yin 2026 writes equation 1",
        "in weeks and its fixed constants (TM50 = 47.7 weeks, Hill = 3.4)",
        "are only meaningful on the week scale. See the PAGE entry in",
        "inst/references/covariate-columns.md, which explicitly permits a",
        "model to declare weeks when its source equations are written in",
        "weeks. Entered a priori (not by stepwise selection) as a sigmoidal",
        "maturation function on CL only; PMA has no effect on V. The",
        "maturation form was adapted from the Rhodin et al. 2009 glomerular",
        "filtration rate model (Yin 2026 reference 25) because ceftazidime is",
        "primarily renally eliminated (Discussion). Cohort median 39.7 weeks",
        "(range 32.7-41.9 weeks; Table 1). Time-varying in principle, but",
        "time-fixed in practice over this study's few-day sampling window.",
        "The simulation subgroups were PMA 32-35, 35-38 and 38-42 weeks."
      ),
      source_name        = "PMA"
    )
  )

  # Covariates that Yin 2026 collected and screened but did NOT retain in the
  # final model. Documented here for provenance; they are deliberately absent
  # from covariateData because they are never referenced in model().
  # Yin 2026 Results 'PopPK model building': "Only the a priori covariates PMA
  # and WT were retained in the final model, as stepwise covariate selection
  # did not identify any additional significant predictors."
  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Collected (Yin 2026 Methods 'Study design and population') and",
        "entered the stepwise screen (forward P < 0.05, backward P < 0.01)",
        "but not retained. Cohort median 39.4 weeks (range 32.6-41.3;",
        "Table 1). Li et al. (Yin 2026 reference 7) did retain GA."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days in this paper (canonical PNA is months)",
      type        = "continuous",
      notes       = paste(
        "Collected and screened but not retained. Cohort median 1 day",
        "(range 1-4 days; Table 1) -- the cohort spans too narrow a",
        "postnatal-age range for a PNA effect to be identifiable. Li et al.",
        "(reference 7) and van der Veer et al. (reference 11) both retained",
        "PNA in wider-ranging cohorts."
      ),
      source_name = "PNA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Collected and screened but not retained. Cohort median 56.7 umol/L",
        "(range 23.7-88.0; Table 1). The only covariate with missing data",
        "(2 of 72 neonates, 2.8%); missing values were mean-imputed from the",
        "remaining 70 neonates (Discussion). Renal maturation is instead",
        "carried by the fixed PMA sigmoid, so no residual creatinine effect",
        "was detectable."
      ),
      source_name = "SCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 1L,
    n_observations = 140L,
    age_range      = "Postmenstrual age 32.7-41.9 weeks; postnatal age 1-4 days (inclusion required postnatal age <= 28 days)",
    age_median     = "Postmenstrual age 39.7 weeks; postnatal age 1 day",
    weight_range   = "1.8-4.2 kg",
    weight_median  = "3.2 kg",
    sex_female_pct = 40.3,
    race_ethnicity = c(Chinese = 100),
    disease_state  = paste(
      "Neonates with confirmed or suspected bacterial infection requiring",
      "intravenous ceftazidime, with available intravenous access.",
      "14 of 72 (19.4%) were preterm births. Exclusions: known",
      "hypersensitivity to ceftazidime or other cephalosporins,",
      "life-threatening infection, missing key baseline data, or absent",
      "parental consent (Methods 'Study design and population')."
    ),
    renal_function = paste(
      "Serum creatinine median 56.7 umol/L (range 23.7-88.0). Renal function",
      "is represented in the model only through the fixed Rhodin PMA",
      "maturation sigmoid; measured creatinine was screened and dropped."
    ),
    co_medication  = paste(
      "Vasoactive agents 13/72 (18.1%), ampicillin 13/72 (18.1%), diuretics",
      "3/72 (4.2%), prenatal steroid exposure 4/72 (5.6%), hydrocortisone",
      "2/72 (2.8%), meropenem 1/72 (1.4%) (Table 1)."
    ),
    dose_range     = paste(
      "Therapeutic ceftazidime as standard neonatal care. The paper does not",
      "tabulate the administered study regimens; the simulated regimens were",
      "25, 50, 75 and 100 mg/kg given every 6, 8 or 12 h as 30-minute",
      "intravenous infusions (Methods 'Model-based simulation')."
    ),
    regions        = "China (single centre: Xiamen Maternity and Child Health Care Hospital, Xiamen)",
    notes          = paste(
      "Prospective, single-centre, open-label study run 2022-2024. Baseline",
      "demographics per Yin 2026 Table 1 (medians with ranges). 140",
      "capillary whole-blood samples (median 2 per patient, range 1-4) were",
      "collected from the heel with the Capitainer quantitative DBS device",
      "at 0-0.5, 2-3 and 4-6 h after a dose; a metered 10 uL aliquot was",
      "applied per card. Ceftazidime was quantified by UPLC-PDA over",
      "0.25-50.0 ug/mL (accuracy 92.5-101.4%, precision < 10.4%). Estimation",
      "was FOCE-I in NONMEM 7.5 with PsN 5.3.0 and Pirana 2.9.9. A",
      "one-compartment model was selected over two-compartment; adding the a",
      "priori WT and PMA covariates gave OFV -153.8, 31.3 points below the",
      "base model. Internal evaluation: 1000-simulation VPC (94.3% of",
      "observations inside the 90% prediction interval) and a 1000-replicate",
      "non-parametric bootstrap (967 converged). No external validation was",
      "performed. NOTE: the Abstract and Methods say concentrations were",
      "quantified by UPLC-MS/MS while the Methods 'Sample collection and",
      "determination' section says UPLC with photodiode array detection",
      "(UPLC-PDA); the assay detector is reported inconsistently but this",
      "does not affect any model parameter."
    )
  )

  ini({
    # ===== Structural PK (Yin 2026 Table 2, "Final model" column) ==========
    # One-compartment IV disposition parameterised as CL and V, both
    # standardised to a 70 kg reference weight. Values are WHOLE BLOOD
    # referenced (the model was fitted to qDBS whole-blood concentrations).
    lcl <- log(19.6); label("Typical blood clearance standardised to 70 kg and full maturation (L/h)")  # Table 2: CL = 19.6 L/h/70 kg (RSE 4.2%; bootstrap median 19.5, 95% CI 17.9-21.2)
    lvc <- log(70.6); label("Typical blood volume of distribution standardised to 70 kg (L)")           # Table 2: V  = 70.6 L/70 kg   (RSE 7.5%; bootstrap median 71.3, 95% CI 61.3-82.8)

    # ===== Allometric exponents (fixed a priori, not estimated) ============
    # Yin 2026 Methods 'PopPK modeling': "Clearance (CL) and volume of
    # distribution (V) were scaled allometrically based on WT, with exponents
    # of 0.75 for CL and 1 for V". Neither exponent appears in Table 2, i.e.
    # neither carries an RSE or a bootstrap interval -- both were held fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")  # Methods 'PopPK modeling', equation 1: (WT / 70 kg)^0.75
    e_wt_vc <- fixed(1.00); label("Allometric exponent on V (unitless)")   # Methods 'PopPK modeling', equation 2: (WT / 70 kg)^1

    # ===== Renal maturation on CL (fixed to Rhodin 2009 literature values) =
    # Yin 2026 equation 1 carries the Hill sigmoid
    #   PMA^Hill / (PMA^Hill + TM50^Hill)
    # on postmenstrual age. Both constants were fixed, not estimated:
    # Methods 'PopPK modeling' -- "TM50, representing the PMA corresponding to
    # 50% maturation of clearance, was fixed at 47.7 weeks, while the Hill
    # coefficient was fixed at 3.4" -- and Discussion -- "the maturation model
    # was adapted from the glomerular filtration rate model proposed by
    # Rhodin et al. (25), with TM50 (47.7 weeks) and the Hill coefficient
    # (3.4) fixed to literature values, as our data set did not support
    # reliable estimation of these parameters." Neither value is in Table 2.
    tmat50   <- fixed(47.7); label("Postmenstrual age at 50% mature renal clearance (weeks)")  # Methods 'PopPK modeling' + Discussion; Rhodin et al. 2009 (Yin 2026 reference 25)
    hill_mat <- fixed(3.40); label("Hill coefficient of the renal maturation sigmoid (unitless)")  # Methods 'PopPK modeling' + Discussion; Rhodin et al. 2009 (Yin 2026 reference 25)

    # ===== IIV (exponential / log-normal on CL and V) ======================
    # Yin 2026 Methods 'PopPK modeling': "The interindividual variability
    # associated with PopPK parameters was characterized using an exponential
    # model, which is assumed to follow a log-normal distribution with a mean
    # of 0 and a variance of omega^2", and equations 1-2 carry exp(eta_CL) and
    # exp(eta_V).
    #
    # Table 2 reports IIV as a percent coefficient of variation, and its
    # footnote c gives the conversion explicitly:
    #   CV% = sqrt(exp(omega^2) - 1) * 100%
    # (the radical is dropped by pdftotext but the "2" superscript on omega
    # and the "- 1" term make the standard log-normal CV formula unambiguous;
    # no other formula reproduces a percent from a variance).
    # Inverting: omega^2 = log(1 + CV^2). The paper reports no CL-V
    # covariance, so the etas are encoded as independent (diagonal omega).
    etalcl ~ 0.0760  # Table 2: IIV on CL = 28.1% CV (RSE 21.0%, shrinkage 23.3%; bootstrap 27.4%, 95% CI 14.6-37.0) -> omega^2 = log(1 + 0.281^2) = 0.07600
    etalvc ~ 0.1653  # Table 2: IIV on V  = 42.4% CV (RSE 22.9%, shrinkage 33.9%; bootstrap 41.8%, 95% CI 15.0-60.8) -> omega^2 = log(1 + 0.424^2) = 0.16532

    # ===== Residual error =================================================
    # Yin 2026 Methods: "Residual variability was also modeled exponentially."
    # Table 2 nevertheless reports the estimate on the row "Proportional error
    # (%)" and its footnote a explains why: "An exponential error model
    # (additive error model in the log-transformed domain) was used to
    # characterize the residual unexplained variability, which approximates to
    # a proportional error in the normal domain." The value is therefore
    # encoded as nlmixr2's linear-space proportional error, matching both the
    # units the paper reports it in and the interpretation the paper's own
    # footnote assigns to it. See the vignette Errata for the (small)
    # numerical difference from a strict lnorm() encoding at this CV.
    propSd <- 0.223; label("Proportional residual error (fraction)")  # Table 2: proportional error = 22.3% (RSE 45%, shrinkage 31.6%; bootstrap 22.2%, 95% CI 13.1-32.1)
  })

  model({
    # ----- 1. Derived covariate terms -----------------------------------
    # Yin 2026 equation 1, sigmoidal (Hill) renal maturation on PMA:
    #   PMA^Hill / (PMA^Hill + TM50^Hill)
    # At the cohort median PMA of 39.7 weeks this evaluates to 0.349, i.e.
    # a typical study neonate expresses ~35% of the fully mature clearance.
    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    # ----- 2. Individual PK parameters ----------------------------------
    # Yin 2026 equation 1: CL  = TVCL * (WT / 70 kg)^0.75 * fmat * exp(etaCL)
    # Yin 2026 equation 2: V   = TVV  * (WT / 70 kg)^1            * exp(etaV)
    # PMA acts on CL only; V carries size scaling alone.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # ----- 3. Micro-constants -------------------------------------------
    kel <- cl / vc

    # ----- 4. ODE system -------------------------------------------------
    # Ceftazidime is given as an intermittent intravenous infusion (30 min in
    # the paper's simulations) directly into the central compartment; the
    # infusion duration or rate comes from the event table.
    d/dt(central) <- -kel * central

    # ----- 5. Observation and error --------------------------------------
    # WHOLE BLOOD ceftazidime concentration: dose in mg, vc in L -> mg/L.
    # To compare against a plasma-based model, Yin 2026 equations 3-5 convert
    # with the ratio BP = Cplasma / Cblood (reported as 0.72, reference 26)
    # and an unbound fraction fu:
    #   Cplasma          = BP * Cblood
    #   Cplasma,unbound  = fu * BP * Cblood
    # fu is never given a numeric value anywhere in the paper, so the unbound
    # conversion is not reproducible here and is deliberately NOT encoded in
    # the model. See the vignette Errata.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
