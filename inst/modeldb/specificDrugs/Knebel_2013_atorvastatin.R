Knebel_2013_atorvastatin <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for atorvastatin (ATV) and",
    "its principal active metabolite o-hydroxyatorvastatin (o-ATV) in 6-17",
    "year-old children and adolescents with heterozygous familial",
    "hypercholesterolemia (Knebel 2013). Two compartments for each analyte:",
    "first-order oral absorption into an ATV depot, an ATV central/peripheral",
    "pair, and an o-ATV central/peripheral pair fed by the ATV metabolic",
    "flux. The fraction of atorvastatin clearance forming o-ATV (fm) was",
    "fixed to 1 for mathematical identifiability, so all metabolite",
    "parameters are apparent values proportional to the true ones. All",
    "clearances and volumes are apparent (X/F) and are allometrically scaled",
    "by body weight with exponents fixed to 0.75 (clearances) and 1",
    "(volumes) against a 70 kg reference. The first-order absorption rate",
    "constant is reparameterised as Ka = exp(lka + L2), where L2 is the",
    "smaller hybrid disposition rate constant of the parent two-compartment",
    "system, to avoid flip-flop kinetics. Relative bioavailability F1 is 1",
    "for Tanner Stage 2 or above (reference) and 0.752 for Tanner Stage 1.",
    "Dose and both concentrations are molar equivalents. The model was fitted",
    "to the paediatric data alone using NONMEM VI with the PRIOR subroutine;",
    "informative adult priors supported Ka, Vp/F, Vcm/fm, Qm/fm and Vpm/fm",
    "while CL/F, Vc/F and CLm/fm were driven by the paediatric data. Below-",
    "quantitation-limit observations were handled by the M3 censored-",
    "likelihood method, which nlmixr2 does not reproduce here; the residual",
    "error is encoded as the proportional model the paper reports."
  )
  reference <- paste(
    "Knebel W, Gastonguay MR, Malhotra B, El-Tahtawy A, Jen F, Gandelman K.",
    "Population pharmacokinetics of atorvastatin and its active metabolites",
    "in children and adolescents with heterozygous familial",
    "hypercholesterolemia: selective use of informative prior distributions",
    "from adults.",
    "J Clin Pharmacol. 2013;53(5):505-516.",
    "doi:10.1002/jcph.66.",
    sep = " "
  )
  vignette <- "Knebel_2013_atorvastatin"

  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = paste(
      "nmol/L (nM) for Cc (atorvastatin) and Cc_oatv",
      "(o-hydroxyatorvastatin); the analysis dataset converted dose and all",
      "concentrations to molar equivalents using MW 558.64 g/mol for ATV and",
      "574.64 g/mol for o-ATV (Knebel 2013 Methods, Data Assembly)"
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All amounts are nmol because the source analysis
  # dataset was converted to molar equivalents (Knebel 2013 Methods, Data
  # Assembly and Data Analysis).
  compartmentData <- list(
    depot = list(
      analyte = "atorvastatin", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "atorvastatin", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "atorvastatin", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_oatv = list(
      analyte = "o-hydroxyatorvastatin", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_oatv = list(
      analyte = "o-hydroxyatorvastatin", units = "nmol",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Fixed allometric scaling against a 70 kg reference on all eight",
        "disposition parameters of both analytes, with exponents fixed to",
        "0.75 for CL/F, Q/F, CLm/fm and Qm/fm and to 1 for Vc/F, Vp/F,",
        "Vcm/fm and Vpm/fm (Knebel 2013 Supplemental Equation S2; Tables 1",
        "and 2). The exponents were fixed rather than estimated because a",
        "fixed allometric model is the standard approach in paediatrics",
        "(Knebel 2013 Methods, Analysis Methods for Pediatric Data). Cohort",
        "weight ranged 25-99.4 kg (median 47 kg; Supplemental Table S1), so",
        "the 70 kg reference sits in the upper part of the observed range.",
        "Weight is the only covariate retained in the final model besides",
        "Tanner Stage."
      ),
      source_name        = "WT"
    ),
    TANNER_1 = list(
      description        = paste(
        "Tanner Stage 1 (prepubertal) indicator: 1 = Tanner Stage 1,",
        "0 = Tanner Stage 2 or above."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Tanner Stage 2 or above)",
      notes              = paste(
        "Modifies relative bioavailability F1 multiplicatively:",
        "F1 = 1 for Tanner Stage 2 or above (fixed reference) and 0.752 for",
        "Tanner Stage 1 (Knebel 2013 Table 2, F1-Tanner Stage 1 = theta12 =",
        "0.752, bootstrap 95% CI 0.577-1.01). The cohort had 15 Tanner Stage",
        "1 and 24 Tanner Stage 2-or-above patients. The paper is explicit",
        "that this covariate is a composite proxy: 'the Tanner Stage",
        "classification in this analysis represents a composite of potential",
        "factors such as age, different atorvastatin doses, formulation,",
        "and/or treatment compliance' (Discussion). Tanner Stage 1 patients",
        "received 5 mg daily of a chewable tablet and Tanner Stage 2-or-above",
        "patients 10 mg daily of the marketed tablet, so dose level and",
        "formulation are confounded with the stage indicator. The paper",
        "argues formulation is unlikely to be the driver because the two",
        "formulations are bioequivalent in adults. The 95% CI includes the",
        "null value of 1, so the effect is not firmly established."
      ),
      source_name        = "Tanner Stage"
    )
  )

  # Screened but NOT retained in the final model: age was one of the
  # pre-specified covariates of interest for CL/F under the full covariate
  # model approach, but weight entered as a fixed allometric relationship and
  # no separate age effect was carried forward (Knebel 2013 Results,
  # Population Pharmacokinetic Modeling Results - Atorvastatin: "Weight was
  # present in the base model as a fixed allometric relationship so Tanner
  # Stage was added to the base model as a predictor of F1 to create a full
  # covariate model"; Discussion: "physiologic or maturational changes in PK
  # due to age are unlikely within the age range studied here"). No age
  # coefficient is reported, so none can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Pre-specified as a candidate predictor of CL/F (Knebel 2013",
        "Methods, Data Assembly and Data Analysis) but not retained in the",
        "final model and no point estimate reported. Cohort age ranged 6-17",
        "years (median 12; Supplemental Table S1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    n_observations = 310L,
    age_range      = "6-17 years",
    age_median     = "12 years",
    weight_range   = "25-99.4 kg",
    weight_median  = "47 kg",
    disease_state  = paste(
      "Genetically verified heterozygous familial hypercholesterolemia with",
      "baseline low-density lipoprotein cholesterol > 4 mmol/L."
    ),
    dose_range     = paste(
      "Atorvastatin once daily by mouth for 8 weeks. Tanner Stage 1: 5 mg",
      "of a chewable tablet; Tanner Stage 2 or above: 10 mg of the marketed",
      "tablet. Doses were doubled at Week 4 in subjects who had not reached",
      "target LDL-C < 3.35 mmol/L."
    ),
    maturation     = "15 Tanner Stage 1 and 24 Tanner Stage 2-or-above patients",
    notes          = paste(
      "Single open-label 8-week paediatric trial (Pfizer study A2581172).",
      "Sparse sampling: 8 blood samples per subject - one between 4 and 12 h",
      "postdose at each of Weeks 2 and 6, and predose plus 1 h and 2 h",
      "postdose at each of Weeks 4 and 8. Assay linear range 0.250-100",
      "ng/mL (LLOQ 0.250 ng/mL); 51 of 310 atorvastatin observations (16%)",
      "and 37 (12%) of the o-hydroxyatorvastatin observations were below",
      "the limit of quantitation and were retained as censored data under",
      "the M3 method. p-hydroxyatorvastatin was BQL in more than 80% of",
      "samples and was therefore not modelled. Sex and race distributions",
      "are not reported: Supplemental Table S1 gives only weight and age.",
      "The three adult studies cited in the paper contributed prior",
      "distributions only, not observations, to the final model. Baseline",
      "demographics: Knebel 2013 Supplemental Table S1."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters - final combined parent-metabolite model
    # (Knebel 2013 Table 2 / Supplemental Table S6, Run 340), all at the
    # 70 kg reference weight. Supplemental Equation S2 holds every fixed
    # effect on the log scale
    # (e.g. CL/F_i = exp(theta1 + log(WT_i/70) * 0.75 + eta_CL/F_i)),
    # so Table 2's linear-scale point estimates are exp(theta) and the
    # log() wrappers below reproduce theta directly.
    # -----------------------------------------------------------------
    lcl <- log(699);  label("Atorvastatin apparent oral clearance CL/F at 70 kg (L/h)")               # Knebel 2013 Table 2: CL/F (theta1) = 699 L/h, bootstrap 95% CI (570, 881)
    lvc <- log(1020); label("Atorvastatin apparent central volume Vc/F at 70 kg (L)")                 # Knebel 2013 Table 2: Vc/F (theta2) = 1020 L, bootstrap 95% CI (209, 2210)
    lq  <- log(227);  label("Atorvastatin apparent inter-compartmental clearance Q/F at 70 kg (L/h)") # Knebel 2013 Supplemental Table S6 Run 340 (Q/F = 227) and Discussion p.514 ("227 L/hour (80.2, 470)"); Table 2 prints 277, a digit transposition - see vignette Errata
    lvp <- log(1960); label("Atorvastatin apparent peripheral volume Vp/F at 70 kg (L)")              # Knebel 2013 Table 2: Vp/F (theta5) = 1960 L, bootstrap 95% CI (1390, 2460); informative adult prior
    lka <- log(0.2);  label("Atorvastatin absorption rate offset, exp(lka) (1/h)")                    # Knebel 2013 Table 2: Ka (theta4) = 0.2 1/h, bootstrap 95% CI (0.139, 0.304); informative adult prior

    # -----------------------------------------------------------------
    # o-hydroxyatorvastatin (o-ATV) apparent parameters. Because fm was
    # fixed to 1 (below), each is proportional to the true value with the
    # true fraction metabolised as the proportionality constant (Knebel
    # 2013 Results, Combined Parent and Metabolite Model).
    # -----------------------------------------------------------------
    lcl_oatv <- log(616);  label("o-ATV apparent oral clearance CLm/fm at 70 kg (L/h)")                 # Knebel 2013 Table 2: CLm/fm (theta10) = 616 L/h, bootstrap 95% CI (518, 748); no prior
    lvc_oatv <- log(401);  label("o-ATV apparent central volume Vcm/fm at 70 kg (L)")                   # Knebel 2013 Table 2: Vcm/fm (theta6) = 401 L, bootstrap 95% CI (272, 715); informative prior from the ATV-only fit
    lq_oatv  <- log(368);  label("o-ATV apparent inter-compartmental clearance Qm/fm at 70 kg (L/h)")   # Knebel 2013 Table 2: Qm/fm (theta7) = 368 L/h, bootstrap 95% CI (248, 562); informative prior from the ATV-only fit
    lvp_oatv <- log(2040); label("o-ATV apparent peripheral volume Vpm/fm at 70 kg (L)")                # Knebel 2013 Table 2: Vpm/fm (theta8) = 2040 L, bootstrap 95% CI (1740, 2250); informative prior from the ATV-only fit

    # Fraction of atorvastatin clearance forming o-ATV. Supplemental
    # Equation S2 defines f = theta9 and fm = 1 - f; the reported fm is 1,
    # i.e. theta9 = 0, so the whole of CL/F is metabolic formation of o-ATV
    # and there is no competing parent-elimination route in the fitted
    # model. Fixed, not estimated, to make the metabolite sub-model
    # identifiable.
    fm <- fixed(1); label("Fraction of atorvastatin clearance forming o-ATV (unitless)")  # Knebel 2013 Table 2 and Supplemental Table S6 Run 340: fm = 1; Results: "the results of a number of interim models indicated fixing the value of fm to 1 was the best approach to ensure mathematical identifiability"

    # -----------------------------------------------------------------
    # Relative bioavailability. Tanner Stage 2-or-above is the reference
    # level and was fixed to 1; Tanner Stage 1 is estimated relative to it.
    # -----------------------------------------------------------------
    lfdepot           <- fixed(log(1)); label("Relative bioavailability F1 for Tanner Stage 2 or above (unitless)") # Knebel 2013 Table 2: F1-Tanner Stage 2 (theta11) = 1, precision NA (reference level)
    e_tanner_1_fdepot <- 0.752;         label("Relative bioavailability F1 for Tanner Stage 1, vs Tanner Stage 2 or above (unitless)") # Knebel 2013 Table 2: F1-Tanner Stage 1 (theta12) = 0.752, bootstrap 95% CI (0.577, 1.01)

    # -----------------------------------------------------------------
    # Fixed allometric exponents on body weight, reference 70 kg (Knebel
    # 2013 Supplemental Equation S2; Tables 1 and 2 row annotations
    # "CL/F ~ (WT/70)^0.75", "Vc/F ~ (WT/70)^1", etc.). Reported without
    # uncertainty because the paper fixed them: "a power model with the
    # exponents fixed to 0.75 for CL/F and Q/F and 1 for Vc/F and Vp/F"
    # (Results), and "CLm/fm, Vcm/fm, Qm/fm and Vpm/fm were allometrically
    # scaled by weight using a power model with fixed exponents".
    # -----------------------------------------------------------------
    e_wt_cl      <- fixed(0.75); label("Allometric exponent of (WT/70) on CL/F (unitless)")     # Knebel 2013 Supplemental Equation S2, CL/F line: log(WT/70) * 0.75
    e_wt_vc      <- fixed(1);    label("Allometric exponent of (WT/70) on Vc/F (unitless)")     # Knebel 2013 Supplemental Equation S2, Vc/F line: log(WT/70)
    e_wt_q       <- fixed(0.75); label("Allometric exponent of (WT/70) on Q/F (unitless)")      # Knebel 2013 Supplemental Equation S2, Q/F line: log(WT/70) * 0.75
    e_wt_vp      <- fixed(1);    label("Allometric exponent of (WT/70) on Vp/F (unitless)")     # Knebel 2013 Supplemental Equation S2, Vp/F line: log(WT/70)
    e_wt_cl_oatv <- fixed(0.75); label("Allometric exponent of (WT/70) on CLm/fm (unitless)")   # Knebel 2013 Supplemental Equation S2, CLm/fm line: log(WT/70) * 0.75
    e_wt_vc_oatv <- fixed(1);    label("Allometric exponent of (WT/70) on Vcm/fm (unitless)")   # Knebel 2013 Supplemental Equation S2, Vcm/fm line: log(WT/70)
    e_wt_q_oatv  <- fixed(0.75); label("Allometric exponent of (WT/70) on Qm/fm (unitless)")    # Knebel 2013 Supplemental Equation S2, Qm/fm line: log(WT/70) * 0.75
    e_wt_vp_oatv <- fixed(1);    label("Allometric exponent of (WT/70) on Vpm/fm (unitless)")   # Knebel 2013 Supplemental Equation S2, Vpm/fm line: log(WT/70)

    # -----------------------------------------------------------------
    # Inter-individual variability. Exponential (log-normal) IIV on CL/F,
    # Vc/F and CLm/fm, with a covariance between CL/F and Vc/F (Knebel
    # 2013 Results; Supplemental Equation S2 carries eta on exactly these
    # three parameters). ini() takes variances.
    #
    # The three diagonal entries below are the squared omega values from
    # Supplemental Table S7 (Run 340: ETA1 = 0.463 on CL/F, ETA2 = 1.058
    # on Vc/F, ETA3 = 0.434 on CLm/fm), which agree with Table 2's
    # reported CV% column (46.3%, 106%, 43.3% - the paper reports CV% as
    # 100 * sqrt(variance)).
    #
    # Table 2's printed Omega_1.1 point estimate of 0.124 is a digit
    # transposition of 0.214: sqrt(0.124) = 35.2%, not the 46.3% CV the
    # same row states, and 0.124 falls just below its own bootstrap 95% CI
    # of (0.125, 0.293), whereas 0.463^2 = 0.214 reproduces the CV exactly
    # and sits inside the CI. See vignette Errata.
    #
    # A trailing comment inside a multi-line c(...) breaks the rxode2
    # parser, so the per-entry traces are given here rather than inline.
    #   entry 1 = Omega_1.1 CL/F     = 0.463^2 = 0.214369
    #   entry 2 = Omega_1.2 COV      = 0.185 (Table 2; bootstrap CI (-0.139, 0.35)); implied correlation 0.185 / (0.463 * 1.058) = 0.378
    #   entry 3 = Omega_2.2 Vc/F     = 1.058^2 = 1.119364 (Table 2 prints 1.12)
    # -----------------------------------------------------------------
    etalcl + etalvc ~ c(
      0.214369,
      0.185, 1.119364
    )
    etalcl_oatv ~ 0.188356  # 0.434^2; Knebel 2013 Supplemental Table S7 Run 340 ETA3 = 0.434, Table 2 Omega_33 CLm/fm = 0.188 (CV 43.3%), bootstrap 95% CI (0.0868, 0.288)

    # -----------------------------------------------------------------
    # Residual error - proportional for both analytes (Knebel 2013
    # Supplemental Equation S2: C_ij = Chat_ij * (1 + eps_pij)). Table 2
    # reports the sigma variances; ini() takes SDs, so each value below is
    # sqrt(variance). The paper's own "Residual variability as CV" rows
    # (41.2% and 40.6%) confirm the same convention.
    # -----------------------------------------------------------------
    propSd         <- 0.412311; label("Atorvastatin proportional residual SD (fraction)")           # sqrt(0.17);  Knebel 2013 Table 2: sigma_1.1 pro-ATV = 0.17, bootstrap 95% CI (0.124, 0.207); reported CV 41.2% (35.2, 45.5)
    propSd_oatv <- 0.407431; label("o-hydroxyatorvastatin proportional residual SD (fraction)")   # sqrt(0.166); Knebel 2013 Table 2: sigma_2.2 pro-o-ATV = 0.166; reported CV 40.6% (34.2, 45.5); the printed variance CI lower bound of 0.177 is inconsistent with both the point estimate and the CV CI and is read as 0.117 - see vignette Errata
  })

  model({
    # -----------------------------------------------------------------
    # 1. Fixed allometric size factors, reference weight 70 kg (Knebel
    #    2013 Supplemental Equation S2; the log(WT/70) * exponent terms
    #    inside each exp()).
    # -----------------------------------------------------------------
    ref_wt <- 70

    wt_cl      <- (WT / ref_wt)^e_wt_cl
    wt_vc      <- (WT / ref_wt)^e_wt_vc
    wt_q       <- (WT / ref_wt)^e_wt_q
    wt_vp      <- (WT / ref_wt)^e_wt_vp
    wt_cl_oatv <- (WT / ref_wt)^e_wt_cl_oatv
    wt_vc_oatv <- (WT / ref_wt)^e_wt_vc_oatv
    wt_q_oatv  <- (WT / ref_wt)^e_wt_q_oatv
    wt_vp_oatv <- (WT / ref_wt)^e_wt_vp_oatv

    # -----------------------------------------------------------------
    # 2. Individual parameters. Etas sit on CL/F, Vc/F and CLm/fm only.
    # -----------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc + etalvc) * wt_vc
    q  <- exp(lq)           * wt_q
    vp <- exp(lvp)          * wt_vp

    cl_oatv <- exp(lcl_oatv + etalcl_oatv) * wt_cl_oatv
    vc_oatv <- exp(lvc_oatv)               * wt_vc_oatv
    q_oatv  <- exp(lq_oatv)                * wt_q_oatv
    vp_oatv <- exp(lvp_oatv)               * wt_vp_oatv

    # -----------------------------------------------------------------
    # 3. Absorption rate constant, reparameterised to avoid flip-flop
    #    kinetics. L2 is the smaller of the two hybrid rate constants
    #    (beta) of the PARENT two-compartment disposition system, built
    #    from the three transfer rates T1 = CL/Vc, T23 = Q/Vc and
    #    T32 = Q/Vp (Knebel 2013 Supplemental Equation S2; Results: "The
    #    L2 was calculated and incorporated into the calculation of Ka to
    #    avoid 'flip-flop' kinetics").
    # -----------------------------------------------------------------
    t1  <- cl / vc
    t23 <- q  / vc
    t32 <- q  / vp
    l2  <- ((t1 + t23 + t32) -
              sqrt((t1 + t23 + t32)^2 - 4 * t1 * t32)) / 2
    ka  <- exp(lka + l2)

    # -----------------------------------------------------------------
    # 4. Micro-constants, following the NONMEM compartment numbering of
    #    Supplemental Equation S2 (1 = ATV depot, 2 = ATV central,
    #    5 = ATV peripheral, 3 = o-ATV central, 4 = o-ATV peripheral,
    #    0 = eliminated).
    #
    #    K23 is the metabolic formation flux and K20 the competing
    #    parent-elimination flux; with fm fixed to 1, K20 is identically
    #    zero and the whole of CL/F forms o-ATV. Because the analysis
    #    dataset is in molar equivalents, the flux transfers 1 mol of ATV
    #    to 1 mol of o-ATV with no molecular-weight correction.
    #
    #    Equation S2 prints K30 = (CLm/fm) * fm_ATV / (Vcm/fm), reusing the
    #    fm_ATV factor from the K20 line immediately above. With fm fixed to
    #    1 that factor is (1 - fm) = 0, which would make K30 zero and leave
    #    o-ATV never eliminated -- impossible, and contradicted by Figure 1
    #    and by CLm/fm being estimated at 616 L/h with a tight CI. The
    #    factor is read as a copy-paste artifact; see vignette Errata.
    # -----------------------------------------------------------------
    k25 <- q / vc                 # ATV central -> ATV peripheral
    k52 <- q / vp                 # ATV peripheral -> ATV central
    k23 <- cl * fm       / vc     # ATV central -> o-ATV central (formation)
    k20 <- cl * (1 - fm) / vc     # ATV central -> eliminated (non-o-ATV routes)
    k34 <- q_oatv  / vc_oatv      # o-ATV central -> o-ATV peripheral
    k43 <- q_oatv  / vp_oatv      # o-ATV peripheral -> o-ATV central
    k30 <- cl_oatv / vc_oatv      # o-ATV central -> eliminated

    # -----------------------------------------------------------------
    # 5. ODE system (Knebel 2013 Figure 2 diagram and Supplemental
    #    Equation S2 rate constants).
    # -----------------------------------------------------------------
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot -
                               (k20 + k23 + k25) * central +
                               k52 * peripheral1
    d/dt(peripheral1)      <-  k25 * central - k52 * peripheral1
    d/dt(central_oatv)     <-  k23 * central -
                               (k30 + k34) * central_oatv +
                               k43 * peripheral1_oatv
    d/dt(peripheral1_oatv) <-  k34 * central_oatv - k43 * peripheral1_oatv

    # -----------------------------------------------------------------
    # 6. Relative bioavailability. The power form reproduces the paper's
    #    two-level categorical exactly: the factor is
    #    e_tanner_1_fdepot = 0.752 when TANNER_1 = 1 and 1 when
    #    TANNER_1 = 0 (Knebel 2013 Table 2, F1 rows).
    # -----------------------------------------------------------------
    f(depot) <- exp(lfdepot) * e_tanner_1_fdepot^TANNER_1

    # -----------------------------------------------------------------
    # 7. Observations. Both are plasma molar concentrations.
    # -----------------------------------------------------------------
    Cc      <- central      / vc
    Cc_oatv <- central_oatv / vc_oatv

    Cc      ~ prop(propSd)
    Cc_oatv ~ prop(propSd_oatv)
  })
}
