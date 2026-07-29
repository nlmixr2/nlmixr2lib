Bertrand_2011_S33138 <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for the investigational",
    "antipsychotic S33138 (parent) and its active metabolite S35424 in adults",
    "with schizophrenia (Bertrand 2011). The final selected structural model",
    "is a two-compartment back-transformation form with a presystemic dose",
    "apportionment Fp into the parent depot vs (1 - Fp) into the metabolite",
    "depot, a shared first-order absorption rate Ka, and four linear",
    "elimination / interconversion clearances: parent elimination via other",
    "pathways (CLpo), parent-to-metabolite formation (CLpm), metabolite",
    "elimination via other pathways (CLmo), and metabolite back-transformation",
    "to parent (CLmp). The two volumes (parent Vp and metabolite Vm) are set",
    "equal to a single volume V for identifiability per the source paper.",
    "CYP2D6 poor metabolizers carry a 34% decrease in CLmo (the genetic",
    "covariate retained in the final model). Linear dose-level effects on",
    "the bioavailability f and the parent-fraction Fp are encoded with a 10 mg",
    "reference dose. Parameter values are from Table IV 'With the genetic",
    "covariate' column (SAEM in MONOLIX, closed-form coding, N = 99 patients",
    "with available CYP2D6 genotyping)."
  )
  reference <- paste(
    "Bertrand J, Laffont CM, Mentre F, Chenel M, Comets E. Development of a",
    "complex parent-metabolite joint population pharmacokinetic model.",
    "AAPS J. 2011 Sep;13(3):390-404. doi:10.1208/s12248-011-9282-9."
  )
  vignette <- "Bertrand_2011_S33138"
  units <- list(
    time = "hour",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Current administered dose level of S33138 in mg per dose record. Used as a continuous covariate driving the linear dose effects on the bioavailability f (Bertrand 2011 Eq. 1) and the parent presystemic fraction Fp (Eq. 2). Reference dose 10 mg; observed levels 5, 10, and 20 mg.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record (occasion) dose level supplied as a covariate column alongside the AMT event for the same dose. The current implementation expects DOSE = 5, 10, or 20 (mg) at every record so that the linear dose-effect terms (DOSE - 10) evaluate correctly. The DOSE covariate intentionally uses mg (paper-native unit for the dose-effect coefficient values reported in Table IV) even though the dosing AMT itself is in nmol (molar unit matching the molar concentrations); the user supplies DOSE in mg per-record and converts the mg dose to nmol for the AMT event externally (1 mg S33138 = 1e6 / 319.4 = 3130.87 nmol).",
      source_name        = "Dose"
    ),
    CYP2D6_PM = list(
      description        = "Indicator for the CYP2D6 poor-metabolizer phenotype (1 = CYP2D6 PM, defined by carriers of two non-functional alleles among CYP2D6 *3, *4, *6, *7, *8; 0 = CYP2D6 extensive metabolizer, intermediate metabolizer, or ultrarapid metabolizer). Time-fixed per subject (germline genotype-derived phenotype).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6 extensive, intermediate, or ultrarapid metabolizer)",
      notes              = "Bertrand 2011 'Concentration Measurements and Genetic Polymorphisms': '12 patients were classified as CYP2D6 PM' among 99 subjects with available CYP2D6 genotyping. CYP2D6 PM carries a 34% decrease in the metabolite elimination clearance CLmo (Eq. 3 of the source paper; Table IV with-genetic-covariate column: beta_CLmo,CYP2D6 = -0.42 on log L/h, p = 0.015 by permutation Wald test).",
      source_name        = "CYP2D6 phenotype (PM vs EM)"
    )
  )

  covariatesDataExcluded <- list(
    CYP2C19_PM = list(
      description        = "Indicator for the CYP2C19 poor-metabolizer phenotype (1 = CYP2C19 PM by *2 / *3 allele genotyping; 0 = non-PM). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 non-PM)",
      notes              = "Bertrand 2011 Results 'Covariate Model': 'No effect of the CYP2C19 polymorphisms was found, probably due to the small number of PM' (only 2 of 99 genotyped patients were CYP2C19 PM). Screened during forward selection on CLpo (the expected CYP2C19-related elimination route) but not retained in the final model. Recorded here for completeness so a downstream user knows the covariate was tested and rejected.",
      source_name        = "CYP2C19 phenotype (PM vs EM)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_studies      = 1L,
    n_observations = "626 S33138 + 626 S35424 plasma concentrations from W4 and W8 occasions in the genotyped cohort (Results 'Data'); 8 of 1252 measured were below the limit of quantification and discarded.",
    age_range      = "22-64 years",
    age_median     = "40 years (mean across the N = 101 phase II cohort; Results 'Data')",
    weight_range   = "43.6-120.0 kg",
    weight_median  = "69.0 kg (mean across the N = 101 phase II cohort)",
    sex_female_pct = 46.5,
    race_ethnicity = "Not separately reported in the source paper; the CYP2D6 PM frequency (12 / 99 ~ 12%) is consistent with a European-ancestry cohort (5-10% PM frequency cited in the paper Discussion).",
    disease_state  = "Adults with schizophrenia enrolled in a randomized double-blind multicentre Phase II trial comparing oral S33138 against risperidone (gold-standard comparator).",
    dose_range     = "Once-daily oral S33138 at 5 mg (n = 35), 10 mg (n = 31), or 20 mg (n = 35) for 8 weeks. Phase II patients only; the risperidone arm is excluded from the joint PK modeling cohort.",
    regions        = "International (multicentre Phase II trial; specific countries not separately reported).",
    sampling_window= "Two PK occasions per patient: W4 (week 4 of treatment) and W8 (week 8 of treatment). At each occasion, blood samples were collected pre-dose and at 1, 3, and 6 h post-dose, with exact administration and sampling times recorded. Assumes steady state at each occasion with a 24 h dosing interval.",
    cyp2d6_distribution = "12 / 99 (12.1%) CYP2D6 PM: 4 of the 5 mg cohort, 6 of the 10 mg cohort, 2 of the 20 mg cohort.",
    cyp2c19_distribution = "2 / 99 (2.0%) CYP2C19 PM: 1 in the 5 mg cohort, 1 in the 20 mg cohort.",
    external_validation = "Phase I randomized double-blind tolerance study in 23 healthy men receiving 10, 20, or 30 mg S33138 once daily for 14 days; 203 PK observations for parent and metabolite, sampled at pre-dose and 0.33, 0.66, 1, 1.5, 2, 3, 4, 8, and 12 h post-dose at steady state. Reported in the source paper as the external evaluation cohort; the final model adequately predicted central tendency but underestimated variability in healthy volunteers (paper Discussion).",
    notes          = "Plasma concentrations measured by LC/MS-MS with LOQ 0.56 nmol/L for parent and 0.44 nmol/L for metabolite. The microdose-study volume ratio Vm / Vp was 2.4 [0.5 - 6] (paper Discussion); the final model adopts Vm = Vp for identifiability. The within-subject variance on CLpo (gamma_CLpo = 0.82 on log scale) is reported in Table IV but not encoded in this typical-value library model; this is documented in the vignette's Assumptions and deviations section."
  )

  ini({
    # All parameter values are from Bertrand 2011 Table IV, column 'With the
    # genetic covariate' (N = 99 patients with available CYP2D6 genotyping;
    # SAEM algorithm in MONOLIX with the closed-form encoding of the four-
    # compartment back-transformation + first-pass + dose-apportionment
    # model). See the vignette for the full source trace.

    # Absorption rate constant. Fixed at its W4-only estimate of 8.06 1/h per
    # the paper Results 'Once included the data at W8, the absorption constant
    # rate had to be fixed to its estimate on data at W4 for stability
    # purposes.' IIV on log Ka remains estimated.
    lka          <- fixed(log(8.06));        label("Absorption rate constant Ka (1/h), fixed")                                  # Table IV K = 8.06 (no RSE)

    # Single volume of distribution shared by parent and metabolite (Vm set
    # equal to Vp for identifiability per Methods 'Joint Pharmacokinetic
    # Model'); estimates of V, Fp, and the four clearances are 'apparent
    # parameters tagged by an asterisk * as a reminder of their reliance on
    # the assumption made on volumes' (paper text).
    lvc          <- log(19.4);               label("Apparent shared volume V = Vp = Vm (L)")                                   # Table IV V = 19.4 (RSE 5%)

    # Paper's bioavailability `f` ('fraction of dose after absorption'). Pop
    # value is fixed at 1.00 per Methods 'For identifiability purposes, the
    # fraction of dose available after absorption (f) was set to 1.' The
    # encoded `lfdepot` carries the log of f. IIV on log f is retained
    # (omega_f = 0.27) so individuals deviate around 1 even though the
    # population typical value is fixed.
    lfdepot      <- fixed(log(1.00));        label("Bioavailability f (fraction, -), fixed")                                    # Table IV f = 1.00 (fixed, no RSE)

    # Linear dose effect on f (Eq. 1): f(Dose) = f_pop * exp(e_dose_fdepot *
    # (DOSE - 10)), with DOSE in mg and 10 mg the reference. Source-paper
    # Table IV labels the units of beta_f,D as 'nmol-1' but the Results
    # numerical predictions ('f* was 10% higher for a 5 mg dose and 19% lower
    # for a 20 mg dose') match beta_f,D = -0.02 mg-1 to within rounding; the
    # 'nmol-1' label is interpreted as a paper typo for 'mg-1'. See vignette
    # Errata for the unit-correction rationale.
    e_dose_fdepot <- -0.02;                  label("Linear dose effect on log(f) per mg (1/mg)")                                # Table IV beta_f,D = -0.02 (RSE 24%); unit corrected from 'nmol-1' to 'mg-1'

    # Logit-transformed parent presystemic fraction Fp ('fraction of parent
    # reaching systemic circulation after absorption'). Methods 'a logit
    # model to force 0 <= Fp <= 1'. No IIV on logit(Fp) in the final model
    # (paper Results: 'between-subject variance was found to improve BIC for
    # all parameters except the fraction of dose that escaped first-pass
    # effect (Fp*), CLmp*, and CLpm*').
    logitfp      <- logit(0.87);             label("Logit of parent presystemic fraction Fp (logit scale, -)")                  # Table IV Fp = 0.87 (RSE 2%)

    # Linear dose effect on logit(Fp) (Eq. 2): logit(Fp(Dose)) = logit(Fp_pop)
    # + e_dose_fp * (DOSE - 10). Same unit interpretation (mg-1, not nmol-1
    # as labelled in Table IV) as e_dose_fdepot. The paper Results predict
    # 'Fp* was 22% higher for a 5 mg dose and 33% lower for a 20 mg dose';
    # the logit-additive form encoded here produces a smaller change
    # (~3.4% at 5 mg) than the Results paragraph implies, suggesting the
    # paper computed the Results percentages from a multiplicative-
    # exponential form (Fp = Fp_pop * exp(beta * (Dose-10))) that violates
    # the Fp <= 1 constraint at the 5 mg dose. See vignette Errata.
    e_dose_fp    <- -0.06;                   label("Linear dose effect on logit(Fp) per mg (1/mg)")                            # Table IV beta_Fp,D = -0.06 (RSE 27%); unit corrected from 'nmol-1' to 'mg-1'

    # Parent S33138 clearance via other pathways (CLpo). Paper-mechanistic
    # clearance distinct from canonical `lcl` because the parent has two
    # parallel elimination flows (other pathways + transformation to
    # metabolite). Tagged with '*' in the paper as apparent under Vm = Vp.
    lclpo        <- log(0.67);               label("Parent clearance via other pathways CLpo (L/h)")                            # Table IV CLpo = 0.67 (RSE 11%)

    # Parent S33138 to metabolite S35424 transformation clearance (CLpm).
    # Paper-mechanistic; analogous to a metabolite-formation clearance.
    lclpm        <- log(2.09);               label("Parent-to-metabolite transformation clearance CLpm (L/h)")                  # Table IV CLpm = 2.09 (RSE 5%)

    # Metabolite S35424 clearance via other pathways (CLmo). Routed through
    # CYP2D6 (paper Discussion).
    lclmo        <- log(0.50);               label("Metabolite clearance via other pathways CLmo (L/h)")                        # Table IV CLmo = 0.50 (RSE 7%)

    # CYP2D6 PM effect on log(CLmo). Beta = -0.42 corresponds to a 34%
    # decrease (exp(-0.42) = 0.657, so CLmo_PM = 0.66 * CLmo_EM). Paper
    # Results 'The genetic effect analysis indicated that CLmo* was decreased
    # by 34% (p value = 0.015, Wald test by permutation) in CYP2D6 PM
    # patients.'
    e_cyp2d6_pm_clmo <- -0.42;               label("CYP2D6 PM effect on log(CLmo) (log L/h)")                                  # Table IV beta_CLmo,CYP2D6 = -0.42 (RSE 40%); -34% in CYP2D6 PM

    # Metabolite S35424 back-transformation to parent clearance (CLmp). The
    # back-transformation flux is the mechanism that drives the long
    # terminal half-life of the parent observed in this study (Discussion).
    lclmp        <- log(0.09);               label("Metabolite back-transformation clearance CLmp (L/h)")                       # Table IV CLmp = 0.09 (RSE 12%)

    # Between-subject variability. Source paper reports the standard
    # deviation of the log-normal (exponential / logit) random effect;
    # variance for nlmixr2 = SD^2. IIVs that the paper dropped from the
    # final model (Fp, CLpm, CLmp) are NOT encoded here per the source
    # paper's final structure.
    etalfdepot   ~ 0.0729                                                                                                       # Table IV omega_f    = 0.27  (RSE 12%) -> 0.27^2
    etalka       ~ 1.9881                                                                                                       # Table IV omega_Ka   = 1.41  (RSE 25%) -> 1.41^2 (Ka pop is fixed but IIV is estimated)
    etalvc       ~ 0.0676                                                                                                       # Table IV omega_V    = 0.26  (RSE 13%) -> 0.26^2
    etalclpo     ~ 0.2116                                                                                                       # Table IV omega_CLpo = 0.46  (RSE 38%) -> 0.46^2
    etalclmo     ~ 0.25                                                                                                         # Table IV omega_CLmo = 0.50  (RSE  8%) -> 0.50^2

    # Residual error. Parent S33138: proportional only (paper Results 'a
    # proportional and a combined error model were selected for the parent
    # drug and its metabolite, respectively'; a_p was fixed to 0 per Table
    # IV footnote). Metabolite S35424: combined additive + proportional.
    propSd          <- 0.30;                 label("Parent proportional residual SD (fraction)")                                # Table IV b_p = 0.30 (RSE 3%); a_p fixed to 0
    addSd_s35424    <- 66.5;                 label("Metabolite additive residual SD (nmol/L)")                                  # Table IV a_m = 66.5 nmol/L (RSE 14%)
    propSd_s35424   <- 0.06;                 label("Metabolite proportional residual SD (fraction)")                            # Table IV b_m = 0.06 (RSE 9%)
  })

  model({
    # ---------------------------------------------------------------
    # Per-record covariate-derived expressions. DOSE is the per-record
    # administered dose in mg (5 / 10 / 20 in the source cohort);
    # CYP2D6_PM is the binary phenotype indicator. Dosing AMT is
    # supplied in nmol (molar) so the depot / central / metabolite
    # states carry molar amounts and concentrations land in nmol/L
    # when divided by the shared volume V (L); 1 mg of S33138 = 1e6 /
    # 319.4 = 3130.87 nmol, which is the conversion the user applies
    # externally when building the dosing dataset (see vignette).
    # ---------------------------------------------------------------
    dose_eff <- DOSE - 10
    cyp_eff  <- e_cyp2d6_pm_clmo * CYP2D6_PM

    # ---------------------------------------------------------------
    # Individual parameters. Exponential / logit transforms match the
    # paper's stated IIV form. The dose effect on f is applied inside
    # exp() (multiplicative on the linear scale per Eq. 1); the dose
    # effect on Fp is applied inside expit() on the logit scale per
    # Eq. 2 and the paper's logit-IIV statement.
    # ---------------------------------------------------------------
    ka        <- exp(lka     + etalka)
    vc        <- exp(lvc     + etalvc)
    f_bio     <- exp(lfdepot + etalfdepot + e_dose_fdepot * dose_eff)
    fp        <- expit(logitfp           + e_dose_fp     * dose_eff)
    clpo      <- exp(lclpo   + etalclpo)
    clpm      <- exp(lclpm)
    clmo      <- exp(lclmo   + etalclmo  + cyp_eff)
    clmp      <- exp(lclmp)

    # ---------------------------------------------------------------
    # First-order rate constants. With Vm = Vp = V, the four
    # clearances divided by V give the four rate constants of the
    # Appendix Eq. 4 system (kpo, kpm, kmo, kmp).
    # ---------------------------------------------------------------
    kpo <- clpo / vc
    kpm <- clpm / vc
    kmo <- clmo / vc
    kmp <- clmp / vc

    # ---------------------------------------------------------------
    # ODE system. The final selected model (Fig. 2 right-most panel)
    # has Kap = Kam = Ka so the two absorption depots collapse to a
    # single lumped `depot`; the dose-apportionment fraction Fp
    # appears at the depot outflow into the parent central, with
    # (1 - Fp) routed to the metabolite central. This encoding is
    # mathematically equivalent to two separate depots with the dose
    # split between them when Kap = Kam.
    #
    # State units: depot in nmol (after MW conversion from mg via the
    # bioavailability hook); central and central_s35424 in nmol.
    # ---------------------------------------------------------------
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * fp       * depot - (kpo + kpm) * central       + kmp * central_s35424
    d/dt(central_s35424)  <-  ka * (1 - fp) * depot - (kmo + kmp) * central_s35424 + kpm * central

    # Bioavailability hook. f_bio is the paper's `f` (population typical
    # value fixed at 1, with exponential IIV and a linear dose effect).
    # AMT is in nmol so no molecular-weight conversion is applied here.
    f(depot) <- f_bio

    # ---------------------------------------------------------------
    # Observation outputs (concentrations in nmol/L). Parent S33138
    # uses the canonical bare `Cc`; metabolite S35424 carries the
    # `_s35424` metabolite suffix. Residual error matches the paper:
    # proportional for parent, combined additive + proportional for
    # metabolite (with a_p = 0 fixed per Table IV footnote).
    # ---------------------------------------------------------------
    Cc         <- central        / vc
    Cc_s35424  <- central_s35424 / vc

    Cc         ~ prop(propSd)
    Cc_s35424  ~ add(addSd_s35424) + prop(propSd_s35424)
  })
}
