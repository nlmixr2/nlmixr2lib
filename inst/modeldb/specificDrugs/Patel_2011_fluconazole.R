Patel_2011_fluconazole <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous fluconazole in 10",
    "critically ill anuric adults receiving continuous venovenous",
    "hemodiafiltration (CVVHDF) (Patel 2011). Total fluconazole clearance",
    "from the central compartment is partitioned into a CVVHDF-route arm",
    "(CL_CVVHDF, 1.66 L/h typical, encoded as lcl_renal) and a non-CVVHDF",
    "arm (CL_NCVVHDF, 1.01 L/h typical, encoded as lcl_nonren) that are",
    "fitted simultaneously to the plasma concentration-time profile and the",
    "cumulative amount of fluconazole in the CVVHDF effluent. Dialysis-filter",
    "membranes in use for more than 48 hours reduce CVVHDF efficiency to 36.8",
    "percent of the fresh-filter baseline (FILT_AGE_HI indicator on the",
    "CL_CVVHDF arm, multiplicative effect 0.368, bootstrap 95 percent",
    "confidence interval 0.326 to 0.426; informed by 1 of 10 patients with",
    "a > 48 h filter). Input into the central compartment is zero-order over",
    "an estimated infusion duration D1 (typical 0.689 h, near the 60-min",
    "nominal infusion). No subject-level covariates (age, total body weight,",
    "sex, APACHE II score) reached the OFV-3.84 retention threshold on CL,",
    "Vc, Q, or Vp, so no demographic covariates are encoded in this file."
  )
  reference <- paste(
    "Patel K, Roberts JA, Lipman J, Tett SE, Deldot ME, Kirkpatrick CM.",
    "Population pharmacokinetics of fluconazole in critically ill patients",
    "receiving continuous venovenous hemodiafiltration: using Monte Carlo",
    "simulations to predict doses for specified pharmacodynamic targets.",
    "Antimicrob Agents Chemother. 2011;55(12):5868-5874.",
    "doi:10.1128/AAC.00424-11"
  )
  vignette <- "Patel_2011_fluconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FILT_AGE_HI = list(
      description        = "Indicator that the in-use CVVHDF hemofilter membrane has been operating for more than 48 hours; reduces CL_CVVHDF efficiency to 36.8 percent of the fresh-filter baseline (Patel 2011 Table 2: ffCL_CVVHDF = 0.368, bootstrap 95 percent CI 0.326-0.426; Delta-OBJ = -11.46 vs base model).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Patel 2011 Methods state that filter age was encoded from the recorded",
        "filter-swap timestamps. Of the 10 enrolled patients, 1 (patient 4)",
        "had the hemofilter in use for more than 48 hours at the start of the",
        "PK sampling profile; the other 9 patients started PK sampling with a",
        "fresh filter (< 48 h in use). The bootstrap 95 percent CI on the",
        "effect (0.326-0.426) is informed by that single subject combined with",
        "the dense plasma + effluent sampling over the 12 h profile. To",
        "simulate the > 48 h sub-cohort, set FILT_AGE_HI = 1 in the event",
        "table; for the typical fresh-filter case use FILT_AGE_HI = 0",
        "(reference category). Subject-level in the source data; future",
        "datasets that resolve mid-profile filter swaps would carry",
        "FILT_AGE_HI as a within-subject time-varying covariate."
      ),
      source_name        = "FILT_AGE_HI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "51-76 years (median 67)",
    weight_range   = "50-104 kg (median 80)",
    sex_female_pct = 50.0,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill anuric adults requiring continuous venovenous",
      "hemodiafiltration (CVVHDF) for renal failure of any cause and",
      "prescribed fluconazole for a suspected fungal infection. Diagnoses",
      "on admission: medical (5), emergency surgery (5); sites of infection:",
      "urinary tract (3), lung (2), intraabdominal (3), blood (1), indwelling",
      "vascular catheter (1); causative organism Candida albicans in 7,",
      "Candida parapsilosis in 1, nonfungal in 1, not listed in 1. APACHE II",
      "scores 17-44 (median 29). 8 of 10 had normal liver function; 2 had",
      "abnormally elevated ALT or AST (4x upper limit). Serum albumin",
      "11-30 g/L (low in all 10 patients)."
    ),
    dose_range     = paste(
      "200 mg intravenous fluconazole twice daily as a 60-min infusion. PK",
      "sampling: plasma at 0.5, 1, 2, 3, 4, 6, 8, and 12 h after the dose;",
      "CVVHDF effluent collected hourly over 12 h. Patients 3, 4, and 8 were",
      "sampled on the first day of treatment (initial profile); the other 7",
      "patients were sampled on day 3 or day 5 (steady-state profile)."
    ),
    regions        = paste(
      "Royal Brisbane and Women's Hospital, Queensland, Australia",
      "(28-month enrollment window). The University of Queensland REC and",
      "Royal Brisbane Hospital REC approved the protocol."
    ),
    notes          = paste(
      "Baseline demographics from Patel 2011 Table 1. Dialysis prescription",
      "(uniform across patients): predilution filtration solution 2 L/h plus",
      "dialysate 1 L/h giving 3 L/h CVVHDF effluent, 999 mL/h fluid input /",
      "effluent rates controlled by IMED PC4 volumetric pumps, blood flow",
      "200 mL/min via Gambro BMM-10 pump through a Hospal AN69HF",
      "hemofilter. With the exception of patient 4 (filter > 48 h at start),",
      "all patients had filters in use < 48 h. All patients except 1 and 5",
      "began dialysis before fluconazole dosing. NONMEM 6.1 FOCE-I with",
      "ADVAN6; nonparametric bootstrap n = 1000 (Patel 2011 Methods)."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL CLEARANCE PARAMETERS -- Patel 2011 Table 2 (final model;
    # estimate column, bootstrap median in parentheses for reference). The
    # paper partitions total fluconazole clearance from the central
    # compartment into a CVVHDF-route arm (CL_CVVHDF) and a non-CVVHDF arm
    # (CL_NCVVHDF). The CVVHDF arm is encoded as lcl_renal (the canonical
    # name for the dialysis-route extracorporeal "renal-mimetic" arm; see
    # the Multi-component clearance section of parameter-names.md), and the
    # non-CVVHDF arm is encoded as lcl_nonren. These map to the published
    # parameters as:
    #   lcl_renal  <- log(CL_CVVHDF)  = log(1.66)  (Table 2)
    #   lcl_nonren <- log(CL_NCVVHDF) = log(1.01)  (Table 2)
    # The typicals are anchored to the fresh-filter case (FILT_AGE_HI = 0);
    # the > 48 h filter case is applied as a multiplicative reduction on
    # the lcl_renal arm via e_filt_age_hi_cl_renal in model().
    # -----------------------------------------------------------------------
    lcl_renal  <- log(1.66); label("Typical CL_CVVHDF (dialysis-route fluconazole CL) at FILT_AGE_HI = 0 (L/h)")  # Patel 2011 Table 2: CL_CVVHDF = 1.66 L/h
    lcl_nonren <- log(1.01); label("Typical CL_NCVVHDF (non-CVVHDF fluconazole CL) (L/h)")                          # Patel 2011 Table 2: CL_NCVVHDF = 1.01 L/h
    lvc        <- log(31.7); label("Typical central volume of distribution Vc (L)")                                  # Patel 2011 Table 2: Vc = 31.7 L
    lvp        <- log(21.9); label("Typical peripheral volume of distribution Vp (L)")                               # Patel 2011 Table 2: Vp = 21.9 L
    lq         <- log(27.6); label("Typical inter-compartmental clearance Q (L/h)")                                  # Patel 2011 Table 2: Q = 27.6 L/h
    ldur       <- log(0.689); label("Typical zero-order IV infusion duration D1 (h)")                                # Patel 2011 Table 2: D1 = 0.689 h

    # -----------------------------------------------------------------------
    # COVARIATE EFFECTS
    #
    # ffCL_CVVHDF = 0.368 (Patel 2011 Table 2): multiplicative factor on the
    # CVVHDF clearance arm when the hemofilter has been in use for more than
    # 48 h. Encoded as e_filt_age_hi_cl_renal^FILT_AGE_HI in model() so that
    # FILT_AGE_HI = 0 (fresh filter) gives multiplier 1 and FILT_AGE_HI = 1
    # (> 48 h filter) gives multiplier 0.368 (Patel 2011 Discussion: filters
    # in use > 48 h reduced CVVHDF efficiency to 37 percent of total
    # fluconazole clearance). The effect is estimated (not fixed); the
    # bootstrap 95 percent CI is 0.326-0.426 (Patel 2011 Table 2).
    # -----------------------------------------------------------------------
    e_filt_age_hi_cl_renal <- 0.368; label("Multiplicative factor on CL_CVVHDF when FILT_AGE_HI = 1 (filter > 48 h)")  # Patel 2011 Table 2: ffCL_CVVHDF = 0.368

    # -----------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Patel 2011 Table 2, %CV column.
    #
    # Methods state "Between-subject variability (BSV) was calculated using
    # an exponential variability model and was assumed to follow a lognormal
    # distribution." Under the lognormal interpretation, %CV relates to the
    # log-scale variance as omega^2 = log(1 + CV^2) (see
    # references/verification-checklist.md and parameter-names.md).
    #
    #   BSV CL_CVVHDF  = 19.8 % CV -> omega^2 = log(1 + 0.198^2)
    #   BSV CL_NCVVHDF = 77.1 % CV -> omega^2 = log(1 + 0.771^2)
    #   BSV Vc         = 22.9 % CV -> omega^2 = log(1 + 0.229^2)
    #   BSV D1         = 23.0 % CV -> omega^2 = log(1 + 0.230^2)
    #
    # No IIV on Q or Vp (none reported in Patel 2011 Table 2). The paper
    # does not report correlations between etas; etas are treated as
    # independent. The reported eta-shrinkages on CL (0.03), Vc (0.10), and
    # D1 (0.08) are small enough that the IIVs are considered well
    # estimated (Patel 2011 Results: "Evaluation of eta shrinkage on
    # clearance, volume of distribution, and infusion duration suggested
    # that none of the parameters were poorly estimated").
    # -----------------------------------------------------------------------
    etalcl_renal  ~ log(1 + 0.198^2)  # Patel 2011 Table 2: BSV CL_CVVHDF = 19.8 % CV
    etalcl_nonren ~ log(1 + 0.771^2)  # Patel 2011 Table 2: BSV CL_NCVVHDF = 77.1 % CV
    etalvc        ~ log(1 + 0.229^2)  # Patel 2011 Table 2: BSV Vc = 22.9 % CV
    etaldur       ~ log(1 + 0.230^2)  # Patel 2011 Table 2: BSV D1 = 23.0 % CV

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Patel 2011 Table 2.
    #
    # Plasma: combined exponential + additive (Patel 2011 Methods:
    # "Residual unexplained variability (RUV) was modeled using a combined
    # exponential and additive random error"). NONMEM's "exponential
    # random error" is log-additive on the observation, which maps to a
    # proportional residual on the linear concentration scale in nlmixr2
    # (see references/parameter-names.md residual-error section).
    #   addSd  = 0.239 mg/L  (RUVSDP, Table 2)
    #   propSd = 0.0367      (RUVCVP = 3.67 % CV, Table 2; fraction units)
    #
    # CVVHDF effluent: additive only (Patel 2011 Results: "The residual
    # variability for the amount of fluconazole in the CVVHDF effluent
    # was best described by an additive error model"). The structural
    # observation in this file is the cumulative amount of fluconazole
    # in the effluent (urineAmt, mg), per the Krekels 2015 paracetamol
    # and Taubert 2018 finafloxacin urine-compartment precedent. Patel
    # 2011 Table 2 reports RUVSDC = 2.84 with units "mg/liters"; the
    # paper's Methods text describes the fit as being on "cumulative
    # amounts of fluconazole in the CVVHDF effluent", so the structural
    # observation is mass (mg), not concentration (mg/L). The "mg/liters"
    # unit label in Table 2 is a paper-side reporting convention (the
    # effluent assay standard curves were prepared at mg/L); the value
    # 2.84 is applied here as a mass SD on the cumulative effluent
    # state, consistent with the model the paper actually fit. See
    # vignette Errata for the unit-labelling discussion.
    # -----------------------------------------------------------------------
    addSd            <- 0.239 ; label("Additive plasma residual SD (mg/L)")                  # Patel 2011 Table 2: RUVSDP = 0.239 mg/L
    propSd           <- 0.0367; label("Proportional plasma residual SD (fraction)")          # Patel 2011 Table 2: RUVCVP = 3.67 % CV
    addSd_urineAmt   <- 2.84  ; label("Additive cumulative-effluent residual SD (mg)")       # Patel 2011 Table 2: RUVSDC = 2.84
  })

  model({
    # -------------------------------------------------------------------
    # Individual PK parameters.
    #
    # CL_CVVHDF carries lognormal IIV and a multiplicative reduction when
    # FILT_AGE_HI = 1 (e_filt_age_hi_cl_renal^FILT_AGE_HI evaluates to 1
    # when FILT_AGE_HI = 0 and to 0.368 when FILT_AGE_HI = 1; this is the
    # same power-encoding precedent used by Goti 2018 vancomycin for
    # 0.7^HEMODIAL and 0.5^HEMODIAL).
    # -------------------------------------------------------------------
    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * e_filt_age_hi_cl_renal^FILT_AGE_HI
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl        <- cl_renal + cl_nonren
    vc        <- exp(lvc + etalvc)
    vp        <- exp(lvp)
    q         <- exp(lq)
    dur_inf   <- exp(ldur + etaldur)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -------------------------------------------------------------------
    # ODE system.
    #
    # Two-compartment IV-infusion disposition with central + peripheral1.
    # Cumulative CVVHDF effluent amount (urine) is fed by the CL_CVVHDF
    # arm only -- the non-CVVHDF arm represents elimination by routes
    # other than the dialysis circuit (residual metabolism / biliary), so
    # only the CL_CVVHDF rate contributes to the cumulative-effluent
    # mass balance. Dose records target central with dur(central) <-
    # dur_inf for zero-order delivery over the (subject-individual) D1.
    # -------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(urine)       <-  (cl_renal / vc) * central

    dur(central) <- dur_inf

    # -------------------------------------------------------------------
    # Observations.
    # Cc       = plasma concentration in mg/L (dose mg, vc L -> mg/L).
    # urineAmt = cumulative fluconazole amount in the CVVHDF effluent in
    #            mg, paired with addSd_urineAmt residual error. The paper
    #            calls this "cumulative amount of fluconazole in the
    #            CVVHDF effluent" (Methods, CVVHDF model paragraph).
    # -------------------------------------------------------------------
    Cc       <- central / vc
    urineAmt <- urine

    Cc       ~ add(addSd) + prop(propSd)
    urineAmt ~ add(addSd_urineAmt)
  })
}
