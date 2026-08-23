Lam_2025_ondansetron <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for ondansetron in neonates with neonatal opioid withdrawal syndrome (NOWS). Bayesian (Metropolis-Hastings MCMC, NONMEM METHOD=BAYES) update of an intravenous infant reference model: informative normal priors on log CL, V, V2, Q were computed from the reference model at the cohort median birth weight (3.19 kg) and postmenstrual age (39 1/7 weeks) using allometric scaling plus a clearance-maturation function, and a weakly informative prior was used for the previously uncharacterised oral absorption rate constant KA. The selected final model applies NO individual-level covariate effects; the allometric and maturation terms enter only through the priors. Oral bioavailability is fixed at 0.62. Transplacental maternal transfer is carried by initialising the neonatal central compartment at t = 0 to the mother's observed plasma concentration at delivery (covariate CP0_MAT_NGML). Inter-individual variability is estimated on CL only, and residual error is additive."
  reference   <- "Lam K, Mondick JT, Peltz G, Wu M, Kraft WK. Bayesian Population Pharmacokinetic Modeling of Ondansetron for Neonatal Opioid Withdrawal Syndrome. Clin Transl Sci. 2025;18(2):e70147. doi:10.1111/cts.70147"
  vignette    <- "Lam_2025_ondansetron"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Appendix S1 NONMEM control
  # stream ($MODEL COMP=(DEPOT) / COMP=(BABYVC) / COMP=(PERI)) and the
  # Methods (oral ondansetron 0.07 mg/kg; plasma ondansetron assay).
  compartmentData <- list(
    depot       = list(analyte = "ondansetron", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ondansetron", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ondansetron", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # The one covariate the final model reads. Appendix S1 $PK ends with
  #   S2 = V / 1000
  #   A_0(2) = C0 * S2
  # i.e. the neonatal central compartment is initialised, at t = 0, to the
  # amount implied by the mother's observed plasma concentration. Results 3.1:
  # 'More than half of the neonates had pre-dose concentrations due to
  # maternal transfer, which were captured by initializing the neonatal
  # central compartment to the observed maternal concentrations.' This is a
  # structural feature of the published model (it lives in $PK, not in the
  # dataset), so it is encoded in model() and declared here.
  covariateData <- list(
    CP0_MAT_NGML = list(
      description        = paste(
        "Observed maternal plasma ondansetron concentration at delivery, used to",
        "initialise the neonate's central compartment at t = 0 (transplacental",
        "transfer).",
        sep = " "
      ),
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Maternal blood was drawn within 30 min of delivery (Methods 2.1);",
        "mothers received ondansetron 8 mg i.v. within 4 h of delivery. Cohort",
        "mean (SD) 28.2 (23.5) ng/mL over 29 maternal plasma samples (Table 1).",
        "More than half of the neonates had a measurable pre-dose concentration",
        "(Results 3.1). The umbilical cord was not modelled as an intermediate",
        "compartment between the maternal and neonatal central compartments",
        "(Discussion), so this single number carries the whole of in-utero",
        "exposure. Set it to 0 for a neonate with no maternal transfer.",
        sep = " "
      ),
      source_name        = "C0 (Appendix S1 $INPUT), ng/mL"
    )
  )

  # The selected final model carries no individual-level covariate effects
  # (Results 3.1: 'the model without covariate effects was selected'; the
  # Appendix S1 $PK block is CL=EXP(MU_1+ETA(1)), V=EXP(XMU2), V2=EXP(XMU3),
  # CLD=EXP(XMU4), KA=EXP(XMU5) with no WT / PMA / PNA terms). Every
  # demographic the paper screened is therefore recorded here for provenance
  # only -- none of these names is referenced in model().
  covariatesDataExcluded <- list(
    WT_BIRTH = list(
      description        = "Birth weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened, not retained at the individual level. Birth weight enters the",
        "ANALYSIS only through the informative priors: the infant reference model",
        "normalised all parameters to 10.4 kg with fixed allometric exponents",
        "(0.75 on CL and Q, 1 on V and V2), and the neonate priors were evaluated",
        "at the cohort median birth weight of 3.19 kg (Table 2 footnote; Methods",
        "2.2.1). The CL-prior sensitivity analysis swept birth weight over",
        "0.5 / 1.595 / 5 kg (Table 2); ELPD favoured larger birth weights but the",
        "diagnostic plots did not, so the reference-weight prior was retained.",
        "Cohort mean (SD) birth weight 3.1 (0.4) kg (Table 1).",
        sep = " "
      ),
      source_name        = "WT (Appendix S1 $INPUT)"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened, not retained at the individual level. Postmenstrual age enters",
        "the ANALYSIS only through the CL prior, via the infant reference model's",
        "maturation function 1 - beta_CL * exp(-(AGE - 1) * ln(2) / T_CL) with",
        "beta_CL = 0.76 and T_CL = 3.82 months (Methods 2.2.1; Table 3 rows",
        "'B CL' and 'T CL'), evaluated at the cohort median postmenstrual age of",
        "39 1/7 weeks (Table 2 footnote). PMA was preferred over PNA despite PNA",
        "giving higher ELPD, because the PNA-based models showed diagnostic",
        "misspecification (Results 3.1; Discussion). The Table 2 sensitivity",
        "sweep over PMA 37 0/7 and 42 4/7 weeks changed ELPD little.",
        "Source column is in weeks; the canonical PAGE unit is months.",
        sep = " "
      ),
      source_name        = "PMA (Appendix S1 $INPUT), weeks"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and rejected. CL priors computed with PNA (cohort median 2 days)",
        "gave higher ELPD than their PMA counterparts (Table 2, 'No covariate",
        "structures with PNA on CL prior' ELPD 676.9 and 'Full covariate structures",
        "with PNA on CL prior' ELPD 723.4 vs final model 591.5), but the",
        "corresponding diagnostic plots depicted misspecification, so PMA was used",
        "(Results 3.1). Source column is in days; the canonical PNA unit is months.",
        sep = " "
      ),
      source_name        = "CA (Appendix S1 $INPUT), days"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried in the analysis dataset ($INPUT GA) and used to stratify the",
        "exploratory exposure-response figures into early term (37 0/7-38 6/7",
        "weeks), term (39 0/7-40 6/7 weeks), late term (41 0/7-41 6/7 weeks) and",
        "post term (>= 42 0/7 weeks) per the joint working group definitions",
        "(Methods 2.3; Figures S16-S17). It is not a covariate on any PK",
        "parameter. Cohort mean (SD) 38.6 (1.1) weeks (Table 1).",
        sep = " "
      ),
      source_name        = "GA (Appendix S1 $INPUT)"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Carried in the analysis dataset ($INPUT SEX) and used only to stratify",
        "the exploratory exposure-response figures (Methods 2.3; Figures S18-S19).",
        "It is not a covariate on any PK parameter. Cohort is 18 male / 18 female",
        "(Table 1). The source SEX column's coding is not stated in the paper; the",
        "canonical SEXF orientation (1 = female) is recorded here for downstream",
        "users and no value transformation is applied by this model.",
        sep = " "
      ),
      source_name        = "SEX (Appendix S1 $INPUT)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "neonates; first blood sample generally within the first 12 h of life, subsequent samples approximately every 24 h. Median postnatal age 2 days",
    age_median     = "postmenstrual age 39 1/7 weeks (median); gestational age at birth 38.6 (1.1) weeks (mean (SD))",
    weight_range   = "birth weight 3.1 (0.4) kg (mean (SD)); cohort median birth weight 3.19 kg",
    weight_median  = "3.19 kg (birth weight)",
    sex_female_pct = 50,
    disease_state  = "Neonates with in-utero opioid exposure at risk of neonatal opioid withdrawal syndrome (NOWS), born to mothers with opioid use disorder and at least 3 weeks of daily opioid exposure before delivery. Substudy of the double-blind, placebo-controlled, multicenter trial NCT01965704 (98 mother/neonate dyads randomised 1:1).",
    dose_range     = "Neonates: ondansetron 0.07 mg/kg orally once every 24 h starting the day of birth, up to five doses (median 2, range 1-5; mean (SD) dose 0.21 (0.03) mg), first dose within 4-8 h of delivery. Mothers: ondansetron 8 mg intravenously within 4 h of delivery, repeated once if delivery had not occurred within 4 h (median 1, range 1-2 maternal doses).",
    regions        = "United States (multicenter)",
    bioanalysis    = "Plasma ondansetron quantified by an assay adapted from a similar study and described in the primary trial report; lower limit of quantification 1.0 ng/mL. Observations below the limit of quantification were excluded from the analysis (Methods 2.2).",
    maternal_transfer = "More than half of the neonates had measurable pre-dose ondansetron concentrations from transplacental maternal transfer. These were handled by initialising the neonatal central compartment to the observed maternal plasma concentration (Appendix S1 $PK: A_0(2) = C0 * S2). Maternal plasma ondansetron concentration at delivery: mean (SD) 28.2 (23.5) ng/mL, 29 maternal plasma samples (Table 1; Results 3, 3.1). The umbilical cord was not modelled as an intermediate compartment (Discussion).",
    notes          = paste(
      "36 neonates and 109 neonatal plasma samples (median 3 per neonate, range 1-5)",
      "plus 29 maternal plasma samples entered the PK analysis; 35 neonates entered",
      "the exploratory exposure-response analysis (one had no Finnegan scores).",
      "Exclusions: one dyad randomised to placebo whose mother received i.v.",
      "ondansetron; one home delivery with no study drug; one neonate with only BLQ",
      "measurements; two neonates who never received ondansetron because of QTc",
      "prolongation; four neonates who received at least one i.v. ondansetron dose;",
      "and one neonate born prematurely (Results 3). Height or length at birth",
      "49.7 (2.7) cm (Table 1).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters: posterior medians of the final updated
    # Bayesian neonate model (Lam 2025 Table 3, 'Posterior median'
    # column). The table note states that all parameters estimated in
    # the log domain were back transformed for clarity, so the values
    # below are the linear-scale medians and are re-logged here.
    #
    # These are ABSOLUTE (not weight-normalised) values for the neonate
    # cohort: the reference model's allometric scaling to 10.4 kg and
    # its clearance-maturation term were used to build the informative
    # priors at the cohort median birth weight (3.19 kg) and median
    # postmenstrual age (39 1/7 weeks), and the selected final model
    # applies no individual-level covariate effects (Results 3.1;
    # Appendix S1 $PK).
    lka     <- log(0.19); label("Absorption rate constant (1/h)")            # Table 3, KA row: 0.19 (95% CDI 0.15, 0.23)
    lcl     <- log(0.58); label("Clearance (L/h)")                           # Table 3, CL row: 0.58 (95% CDI 0.51, 0.67)
    lvc     <- log(0.29); label("Central volume of distribution (L)")        # Table 3, V row: 0.29 (95% CDI 0.26, 0.32)
    lvp     <- log(0.91); label("Peripheral volume of distribution (L)")     # Table 3, V2 row: 0.91 (95% CDI 0.75, 1.11)
    lq      <- log(6.15); label("Intercompartmental clearance (L/h)")        # Table 3, Q row: 6.15 (95% CDI 5.63, 6.74)

    # Oral bioavailability. Methods 2.2: 'Bioavailability (F) was fixed
    # to 0.62 for all model runs'; Table 3 reports F = 0.62 (FIXED).
    # NOTE: the Appendix S1 control stream declares THETA(6) = 0.49 FIX
    # with the comment ';F' -- 0.49 is logit(0.62) = log(0.62 / 0.38) --
    # but no F1 = ... assignment appears in its $PK block, so the printed
    # stream never applies it. The paper's prose and parameter table are
    # taken as authoritative here (see the vignette Errata section for the
    # exposure-based check that supports F = 0.62 over F = 1).
    lfdepot <- fixed(log(0.62)); label("Oral bioavailability (fraction)")    # Table 3, F row; Methods 2.2

    # ---------------------------------------------------------------
    # Inter-individual variability. Only IIV on CL was estimated
    # ('Given the available data, only the interindividual variability
    # (IIV) on CL was estimated', Methods 2.2), carried from the infant
    # reference model as an inverse-Wishart informative prior
    # ($OMEGAP 0.323, $OMEGAPD 14). The value below is the posterior
    # median of the log-scale variance. Table 3's parenthetical 84.2%
    # is the derived CV via the table note
    # %CV = sqrt(exp(posterior median) - 1) * 100; that formula applied
    # to the printed median 0.54 gives 84.6%, the 0.4-point difference
    # reflecting rounding of the median for display.
    # Trailing comment kept quote-free on purpose: this is the one ini() line
    # whose code part has no double quote of its own, so rxode2 converts its
    # trailing comment into a label() call, and a quote in the comment text
    # breaks that conversion.
    etalcl  ~ 0.54                                                            # Table 3, IIV (per cent CV) on CL: 0.54 (84.2) (95 per cent CDI 0.32, 0.92)

    # ---------------------------------------------------------------
    # Residual error: additive, retained from the infant reference model
    # as an inverse-Wishart informative prior ($SIGMAP 0.0257,
    # $SIGMAPD 5); Methods 2.2, Y_ij = Yhat_ij + eps_ij.
    # Table 3 reports the sigma posterior median 80.9 (a variance in
    # (ng/mL)^2) with a parenthetical SD of 9.2. The table note defines
    # SD of sigma = sqrt(posterior median); sqrt(80.9) = 8.99, so the
    # value below is derived from the posterior median via the paper's
    # own stated relationship. The 2% disagreement with the printed 9.2
    # is recorded in the vignette Errata.
    addSd   <- 8.99; label("Additive residual error SD (ng/mL)")             # Table 3, "Residual error: Additive error (SD)" row: 80.9 (9.2) (95% CDI 59.3, 114.6)
  })

  model({
    # 1. Individual parameters. IIV is on CL only.
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    q   <- exp(lq)

    # 2. Micro-constants, matching the Appendix S1 $DES block:
    #    K12 = KA, K20 = CL/V, K23 = CLD/V, K32 = CLD/V2.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Transplacental maternal transfer. Appendix S1 $PK: S2 = V / 1000 and
    #    A_0(2) = C0 * S2, so the neonatal central compartment starts at the
    #    amount corresponding to the mother's observed plasma concentration.
    #    With `central` in mg, `vc` in L and CP0_MAT_NGML in ng/mL, the amount
    #    is CP0_MAT_NGML * vc / 1000 mg, which makes Cc(0) == CP0_MAT_NGML.
    central(0) <- CP0_MAT_NGML * vc / 1000

    # 4. ODE system (Appendix S1 $DES; ADVAN13 with three compartments
    #    DEPOT / BABYVC / PERI).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability on the oral depot.
    f(depot) <- exp(lfdepot)

    # 6. Observation. `central` is in mg and `vc` in L, so central / vc is
    #    mg/L = ug/mL; the factor 1000 converts to the paper's ng/mL, and
    #    reproduces the control stream's scaling S2 = V / 1000.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
