Pejcic_2024_clopidogrel <- function() {
  description <- paste(
    "Joint semi-physiological population PK model for oral clopidogrel",
    "(CLO) and its inactive carboxylic acid metabolite (CLO-CA) in healthy",
    "Caucasian adults, pooled from two 2-way crossover bioequivalence",
    "studies of 150 mg (2 x 75 mg) single oral doses (Pejcic 2024).",
    "Absorption is an Erlang-type transit chain (depot -> transit1 ->",
    "transit2 -> liver) with a single transit rate constant Ktr = 3 / MTT.",
    "The first-pass effect is represented by a hepatic compartment (Vh =",
    "1.5 L/70 kg, liver plasma flow Qh = 50 L/h, both fixed) that receives",
    "the whole absorbed dose and also receives the parent returning from",
    "the systemic circulation at CLP, so that all parent drug is",
    "ultimately metabolised there. Hepatic outflow Qh/Vh is partitioned",
    "into three branches whose fractions sum to one: FP to systemic",
    "clopidogrel, FiaM to CLO-CA, and FaM to the active thiol CLO-TH",
    "(a sink; CLO-TH concentrations were not measured). The fractions use",
    "the source's softmax reparameterisation FiaM = FR1 / (1 + FR1 + FR2),",
    "FaM = FR2 / (1 + FR1 + FR2), FP = 1 / (1 + FR1 + FR2) with FaM fixed",
    "at 12%, giving FiaM = 87.27% (Study 1) and 86.87% (Study 2). A",
    "molecular-weight factor of 0.9565 (307.80 / 321.82) scales the",
    "CLO-CA formation flux. Clopidogrel is one-compartment; CLO-CA is",
    "two-compartment. Body weight is applied allometrically (exponent",
    "0.75 on clearances, 1 on volumes) to a 70 kg reference. MTT, the",
    "generic-product relative bioavailability, FR1, the IIV and IOV",
    "magnitudes and the proportional residual error are all study-",
    "specific, selected by the STUDY_CLO_BE2 indicator; inter-occasion",
    "variability on bioavailability and MTT reflects the crossover",
    "design.")
  reference <- "Pejcic Z, Topic Vucenovic V, Miljkovic B, Vucicevic KM. Integrating Clopidogrel's First-Pass Effect in a Joint Semi-Physiological Population Pharmacokinetic Model of the Drug and Its Inactive Carboxylic Acid Metabolite. Pharmaceutics. 2024;16(5):685. doi:10.3390/pharmaceutics16050685"
  vignette <- "Pejcic_2024_clopidogrel"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL for clopidogrel (Cc); ug/mL for clopidogrel carboxylic acid (Cc_cloca)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling to a 70 kg reference: exponent 0.75 on CLP, CLiaM and QiaM; exponent 1 on Vc,P, Vc,iaM, Vp,iaM and the hepatic volume Vh (Pejcic 2024 Methods section 2.2 final paragraph, and the '/70 kg' unit labels in Table 2). Liver plasma flow Qh is the one flow term NOT scaled: Table 2 labels it 'Qh (L/h)' with no '/70 kg', unlike every other clearance and volume row. Cohort weights were 74.1 +/- 13.56 kg (range 47-100) per Table 1.",
      source_name        = "Body-weight"
    ),
    STUDY_CLO_BE2 = list(
      description        = "Bioequivalence-study indicator for the Pejcic 2024 pooled clopidogrel analysis: 1 = Study 2, 0 = Study 1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Study 1; n = 24, sampling to 48 h, 14 samples per subject per period).",
      notes              = "Selects every study-specific parameter in the joint model: MTT (0.470 vs 0.410 h), the generic relative bioavailability Fgen (1.08 vs 0.960), the fraction parameter FR1 (119 vs 76.8), the IIV on F and FR1, the IOV on F and MTT, and the proportional residual error (41.95% vs 29.39%). Pejcic 2024 Discussion paragraph 3: 'we opted to estimate certain parameters separately for each study ... as the data were derived from different study conditions'. Time-fixed within subject; each subject participated in exactly one study.",
      source_name        = "Study 1 / Study 2"
    ),
    FORM_CLO_GENERIC = list(
      description        = "Clopidogrel product indicator: 1 = the study's generic 75 mg film-coated test tablet, 0 = the Plavix 75 mg film-coated reference tablet.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Plavix 75 mg film-coated tablets, Sanofi Winthrop Industrie, Ambares, France -- the reference medicine in both studies, whose bioavailability F is the fixed 100% anchor).",
      notes              = "Switches bioavailability between the fixed reference value F = 1 and the estimated study-specific relative bioavailability Fgen (Fgen_st1 = 1.08, Fgen_st2 = 0.960; Pejcic 2024 Table 2). Both 95% CIs included one. Each subject received both products, one per period of the 2-way crossover, so this covariate is time-varying within subject and is paired with OCC.",
      source_name        = "formulation (generic / reference)"
    ),
    OCC = list(
      description        = "Crossover period indicator used for inter-occasion variability: 1 = first study period, 2 = second study period.",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- decomposed inside model() into mutually exclusive indicators oc1 = (OCC == 1) and oc2 = (OCC == 2).",
      notes              = "Both studies used the same 2-treatment, 2-period, 2-sequence crossover design, so there are exactly two occasions per subject. IOV was placed on the absorption parameters F and MTT (Pejcic 2024 Methods section 2.2: 'considering the cross-over design of bioequivalence studies, inter-occasional variability (IOV) was incorporated in the absorption parameters'); adding it reduced OFV by 284.202 and AIC by 264.562.",
      source_name        = "period"
    )
  )

  # Demographic and biochemical variables that Pejcic 2024 collected and
  # tabulated (Table 1) but deliberately did NOT carry into the model:
  # "Classical covariate model-building was not applied, as the
  # inclusion/exclusion criteria ensured a homogenous population of
  # healthy subjects, and there were no differences in the demographic
  # characteristics and biochemical parameters among the studies"
  # (Methods section 2.2). Body weight is the sole exception and is a
  # real covariate above.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Table 1: 31.94 +/- 8.51 years (range 19-54). Collected, not modelled."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "Table 1: 21 of 50 subjects female (42.00%). Collected, not modelled."
    ),
    HT = list(
      description = "Height", units = "cm", type = "continuous",
      notes = "Table 1: 177.26 +/- 9.06 cm (range 155-194). Collected, not modelled."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Table 1: 23.40 +/- 2.66 kg/m^2 (range 19.10-29.30). An inclusion criterion (Study 1 19-26, Study 2 19-29), not a model covariate."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = "Table 1: 9.1 +/- 4.43 umol/L (range 3-25). Collected as a hepatic-function screen, not modelled."
    ),
    CREAT = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      notes = "Table 1: 81.5 +/- 17.70 umol/L (range 53-114). Collected, not modelled."
    ),
    ALT = list(
      description = "Alanine transaminase", units = "U/L", type = "continuous",
      notes = "Table 1: 26.2 +/- 9.85 U/L (range 11-52). Collected as a hepatic-function screen, not modelled."
    ),
    AST = list(
      description = "Aspartate transaminase", units = "U/L", type = "continuous",
      notes = "Table 1: 23.7 +/- 4.80 U/L (range 16-35). Collected as a hepatic-function screen, not modelled."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "clopidogrel", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "clopidogrel", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "clopidogrel", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    liver = list(
      analyte = "clopidogrel", units = "mg",
      specimen = "tissue", verified = TRUE
    ),
    central = list(
      analyte = "clopidogrel", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_cloca = list(
      analyte = "clopidogrel carboxylic acid", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_cloca = list(
      analyte = "clopidogrel carboxylic acid", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 2L,
    n_observations = "841 non-zero clopidogrel and 1149 non-zero clopidogrel carboxylic acid plasma concentrations, from 1556 protocol-scheduled samples (54.04% and 73.84% respectively). Samples below the LLOQ were excluded; 87.48% of the excluded clopidogrel BLLOQs were at or after 6 h, and 78.50% of the excluded CLO-CA BLLOQs were at or after 12 h.",
    age_range      = "19-54 years (Table 1; mean 31.94 +/- 8.51). Per the Methods narrative, Study 1 enrolled subjects aged 21-53 and Study 2 aged 19-42.",
    weight_range   = "47-100 kg (Table 1; mean 74.1 +/- 13.56)",
    height_range   = "155-194 cm (Table 1; mean 177.26 +/- 9.06)",
    bmi_range      = "19.10-29.30 kg/m^2 (Table 1; mean 23.40 +/- 2.66). Inclusion criteria: BMI 19-26 kg/m^2 in Study 1 and 19-29 kg/m^2 in Study 2.",
    sex_female_pct = 42.0,
    race_ethnicity = c(White = 100),
    disease_state  = "Healthy adult volunteers meeting predefined inclusion / exclusion criteria for a bioequivalence study; no disease. Baseline hepatic and renal screening values were all within the normal range (Table 1: total bilirubin 9.1 +/- 4.43 umol/L, serum creatinine 81.5 +/- 17.70 umol/L, ALT 26.2 +/- 9.85 U/L, AST 23.7 +/- 4.80 U/L).",
    dose_range     = "Single oral 150 mg dose (2 x 75 mg film-coated tablets) under fasting conditions, given once as the generic test product and once as Plavix 75 mg reference, per the randomisation scheme of a 2-treatment, 2-period, 2-sequence crossover. The 150 mg (double) dose was chosen so that concentrations would be high enough for the assay to quantify.",
    regions        = "Serbia (Military Medical Academy, Belgrade); subjects of Caucasian origin.",
    notes          = "Study 1: n = 24, sampling to 48 h post-dose, 14 samples per subject per period (672 samples). Study 2: n = 26, sampling to 36 h post-dose, 17 samples per subject per period (884 samples). Plasma clopidogrel and CLO-CA were assayed by a validated HPLC-MS method with LLOQ 0.5 ng/mL (clopidogrel) and 0.1 ug/mL (CLO-CA). Observed mean maximum concentrations were 7.2 ng/mL for clopidogrel and 4.9 ug/mL for CLO-CA. No statistically significant differences in subject characteristics between the two studies (Table 1). Data were a secondary use of two sponsor-provided bioequivalence datasets; the underlying data are not publicly available."
  )

  ini({
    # ------------------------------------------------------------------
    # CLOPIDOGREL (PARENT) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Pejcic 2024 Table 2, "Dataset / Estimate" column. All clearances
    # and volumes are reported per 70 kg and are scaled allometrically
    # in model(); the sole exception is Qh, which Table 2 labels
    # "Qh (L/h)" with no per-70-kg normalisation.
    #
    # Structural confirmation of the parent parameterisation: the paper
    # reports a terminal clopidogrel half-life of 1.69 (1.46-1.92) h
    # (Discussion paragraph 1). log(2) * Vc,P / CLP = 0.693 * 218 /
    # 89.5 = 1.688 h, and propagating the Vc,P 95% CI of 188-248 L
    # gives 1.456-1.921 h -- an exact reproduction of both the point
    # estimate and the interval. This confirms that CLP is the
    # first-order loss of parent from the central compartment back into
    # the liver (the Figure 1 "CLP / Vc,P" return arrow), not a hepatic
    # extraction clearance acting on the liver compartment.

    lcl <- fixed(log(89.5))
    label("Clopidogrel clearance from central back into the liver, CLP (L/h/70 kg)")   # Pejcic 2024 Table 2: CL_P = 89.5 FIX. Methods section 2.2: "Initially, clearance of the parent drug (CLP) was estimated at 89.5 L/h, and subsequently, it was fixed at the estimated value in order to stabilize the model."

    lvc <- log(218)
    label("Clopidogrel volume of distribution, Vc,P (L/70 kg)")                        # Pejcic 2024 Table 2: V_c,P = 218 (95% CI 188-248; SIR median 217)

    # ------------------------------------------------------------------
    # HEPATIC (FIRST-PASS) COMPARTMENT
    # ------------------------------------------------------------------
    # Both values were fixed to physiological constants taken from
    # Jiang 2016 (source reference 24), not estimated.

    lqh <- fixed(log(50))
    label("Liver plasma flow, Qh (L/h)")                                               # Pejcic 2024 Table 2: Q_h = 50 FIX; Methods section 2.2: "The model employed a physiological liver volume of 1.5 L and a liver plasma flow (Qh) of 50 L/h [24]." Table 2 gives the unit as "L/h" (not "L/h/70 kg"), so this term is NOT allometrically scaled.

    lvh <- fixed(log(1.5))
    label("Hepatic (liver) compartment volume, Vh (L/70 kg)")                          # Pejcic 2024 Table 2: V_h = 1.5 FIX, unit "L/70 kg" -- allometrically scaled with exponent 1 like the other volumes

    # ------------------------------------------------------------------
    # ABSORPTION -- ERLANG TRANSIT CHAIN, STUDY-SPECIFIC MTT
    # ------------------------------------------------------------------
    # Two transit compartments sit between the depot and the liver, so
    # the dose makes three Ktr-governed transfers and the Savic 2007
    # relation is Ktr = (n_transit + 1) / MTT = 3 / MTT. The paper
    # confirms this arithmetic in Discussion paragraph 4: "Based on the
    # estimated MTTs and two transit compartments, the transit rate
    # constants (Ktr) equaled 6.38 h-1 and 7.32 h-1" -- and 3 / 0.470 =
    # 6.383, 3 / 0.410 = 7.317.

    lmttStudy1 <- log(0.470)
    label("Mean transit time, Study 1 (h)")                                            # Pejcic 2024 Table 2: MTT_st1 = 0.470 (95% CI 0.425-0.515)

    lmttStudy2 <- log(0.410)
    label("Mean transit time, Study 2 (h)")                                            # Pejcic 2024 Table 2: MTT_st2 = 0.410 (95% CI 0.381-0.439)

    # ------------------------------------------------------------------
    # BIOAVAILABILITY
    # ------------------------------------------------------------------
    # F for the Plavix reference product was fixed at 100%; only the
    # generic products' relative bioavailability was estimated, once
    # per study. Both 95% CIs included one, i.e. bioequivalence was
    # not rejected in either study.

    lfdepot <- fixed(log(1))
    label("Bioavailability of the Plavix reference product, F (unitless)")                        # Pejcic 2024 Methods section 2.2: "Bioavailability for the reference medicine (F) was fixed at 100%"

    lfgenStudy1 <- log(1.08)
    label("Relative bioavailability of the Study 1 generic product, Fgen_st1 (unitless)")         # Pejcic 2024 Table 2: F_gen_st1 = 1.08 (95% CI 0.993-1.17)

    lfgenStudy2 <- log(0.960)
    label("Relative bioavailability of the Study 2 generic product, Fgen_st2 (unitless)")         # Pejcic 2024 Table 2: F_gen_st2 = 0.960 (95% CI 0.818-1.10)

    # ------------------------------------------------------------------
    # FIRST-PASS BRANCH FRACTIONS
    # ------------------------------------------------------------------
    # Pejcic 2024 Equations 1-3 partition the hepatic outflow into
    # three branches using a softmax (multinomial-logit) form with the
    # intact-parent branch as the reference category:
    #
    #   FiaM = FR1 / (1 + FR1 + FR2)      (Eq. 1)  -> CLO-CA
    #   FaM  = FR2 / (1 + FR1 + FR2)      (Eq. 2)  -> CLO-TH (sink)
    #   FP   =   1 / (1 + FR1 + FR2)      (Eq. 3)  -> systemic parent
    #
    # FaM was fixed at 12% and FR1 estimated per study, so FR2 is a
    # derived quantity: inverting Eq. 2 gives FR2 = FaM * (1 + FR1) /
    # (1 - FaM). Substituting reproduces the paper's reported fractions
    # exactly -- FR1 = 119 gives FiaM = 87.27% and FR1 = 76.8 gives
    # FiaM = 86.87%, matching Results paragraph 4 to all four reported
    # digits.
    #
    # FR1 is stored on the log scale because the source gives it a
    # log-normal IIV (Table 2 reports IIV(FR1) as a CV%), and because
    # log(FR1) is exactly the logit of the CLO-CA share of the
    # non-CLO-TH flux: FR1 = FiaM / FP, so log(FR1) =
    # logit(FiaM / (FiaM + FP)). That identity is why the canonical
    # `logitfm` name (logit-transformed fraction metabolised) is the
    # right one here rather than a new parameter root; the softmax
    # parameterisation is the same one used by Xie_2019_agomelatine
    # (`ba3 / (1 + ba3 + ba4)`) for its two-metabolite liver split.

    fm_h4 <- fixed(0.12)
    label("Fraction of hepatic outflow metabolised to the active thiol CLO-TH, FaM")   # Pejcic 2024 Methods section 2.2: "FaM was fixed to the estimated value of 12%, which corresponded to the literature data [1-4]." A sensitivity analysis at FaM = 10% and 15% changed the OFV by fewer than two units (Discussion paragraph 3).

    logitfmStudy1 <- log(119)
    label("Logit of the CLO-CA share of non-CLO-TH hepatic outflow, Study 1 (log FR1_st1)")   # Pejcic 2024 Table 2: FR1_st1 = 119 (95% CI 84.3-154); FR1 is the odds FiaM/FP, so log(119) = 4.779 is the corresponding logit

    logitfmStudy2 <- log(76.8)
    label("Logit of the CLO-CA share of non-CLO-TH hepatic outflow, Study 2 (log FR1_st2)")   # Pejcic 2024 Table 2: FR1_st2 = 76.8 (95% CI 64.8-88.8); log(76.8) = 4.341

    # ------------------------------------------------------------------
    # CLOPIDOGREL CARBOXYLIC ACID (CLO-CA) -- TWO-COMPARTMENT
    # ------------------------------------------------------------------

    lcl_cloca <- log(8.70)
    label("CLO-CA clearance, CLiaM (L/h/70 kg)")                                       # Pejcic 2024 Table 2: CL_iaM = 8.70 (95% CI 7.38-10.0)

    lvc_cloca <- log(23.7)
    label("CLO-CA central volume, Vc,iaM (L/70 kg)")                                   # Pejcic 2024 Table 2: V_c,iaM = 23.7 (95% CI 19.7-27.7)

    lq_cloca <- log(10.8)
    label("CLO-CA intercompartmental clearance, QiaM (L/h/70 kg)")                     # Pejcic 2024 Table 2: Q_iaM = 10.8 (95% CI 8.02-13.6)

    lvp_cloca <- log(61.3)
    label("CLO-CA peripheral volume, Vp,iaM (L/70 kg)")                                # Pejcic 2024 Table 2: V_p,iaM = 61.3 (95% CI 50.3-72.3)

    # ------------------------------------------------------------------
    # ALLOMETRIC EXPONENTS
    # ------------------------------------------------------------------

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on clearances with body weight (unitless)")             # Pejcic 2024 Methods section 2.2 final paragraph: "the incorporation of body-weight via allometric scaling was applied to corresponding clearance(s) using an allometric exponent of 0.75"

    e_wt_vc <- fixed(1)
    label("Allometric exponent on volumes with body weight (unitless)")                # Pejcic 2024 Methods section 2.2 final paragraph: "... and volume(s) of distribution with an exponent of 1"

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Pejcic 2024 Table 2 reports IIV and IOV as "Estimate CV (%)" for
    # log-normally distributed parameters (Methods section 2.2:
    # "Inter-individual (IIV) random effects were described by a
    # log-normal distribution of the parameters"). The table carries no
    # footnote giving the CV formula, so the standard log-normal
    # relation omega^2 = log(1 + CV^2) is used throughout -- the same
    # convention documented explicitly in the sibling clopidogrel
    # extraction Danielak_2017_clopidogrel.R. See the vignette Errata
    # for the alternative omega^2 = CV^2 reading.
    #
    # IIV was retained only on Vc,P, Vc,iaM, F and FR1; the paper notes
    # that "When added to the other parameters, IIV resulted in
    # imprecise estimates and over-parameterization of the model."
    # F and FR1 carry separate IIV magnitudes per study.

    etalvc ~ 0.190577                                                                  # Pejcic 2024 Table 2: IIV(V_c,P) = 45.82% CV -> omega^2 = log(1 + 0.4582^2) = 0.190577
    etalvc_cloca ~ 0.060907                                                            # Pejcic 2024 Table 2: IIV(V_c,iaM) = 25.06% CV -> omega^2 = log(1 + 0.2506^2) = 0.060907
    etalfdepotStudy1 ~ 0.167197                                                        # Pejcic 2024 Table 2: IIV(F_st1) = 42.66% CV -> omega^2 = log(1 + 0.4266^2) = 0.167197
    etalfdepotStudy2 ~ 0.064830                                                        # Pejcic 2024 Table 2: IIV(F_st2) = 25.88% CV -> omega^2 = log(1 + 0.2588^2) = 0.064830
    etalogitfmStudy1 ~ 0.425257                                                        # Pejcic 2024 Table 2: IIV(FR1_st1) = 72.80% CV -> omega^2 = log(1 + 0.7280^2) = 0.425257
    etalogitfmStudy2 ~ 0.074753                                                        # Pejcic 2024 Table 2: IIV(FR1_st2) = 27.86% CV -> omega^2 = log(1 + 0.2786^2) = 0.074753

    # ------------------------------------------------------------------
    # INTER-OCCASION VARIABILITY (2-PERIOD CROSSOVER)
    # ------------------------------------------------------------------
    # IOV was placed on the two absorption parameters F and MTT
    # (Methods section 2.2). Each subject contributes two occasions
    # (the two crossover periods), which share a single variance in
    # NONMEM via $OMEGA BLOCK(1) SAME; the second occasion's slot is
    # therefore fixed to the first occasion's estimate.

    etaiov_fdepot_1_Study1 ~ 0.007767                                                  # Pejcic 2024 Table 2: IOV(F_st1) = 8.83% CV -> omega^2 = log(1 + 0.0883^2) = 0.007767
    etaiov_fdepot_2_Study1 ~ fixed(0.007767)                                           # NONMEM $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance
    etaiov_fdepot_1_Study2 ~ 0.052602                                                  # Pejcic 2024 Table 2: IOV(F_st2) = 23.24% CV -> omega^2 = log(1 + 0.2324^2) = 0.052602
    etaiov_fdepot_2_Study2 ~ fixed(0.052602)                                           # NONMEM $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance
    etaiov_mtt_1_Study1 ~ 0.062711                                                     # Pejcic 2024 Table 2: IOV(MTT_st1) = 25.44% CV -> omega^2 = log(1 + 0.2544^2) = 0.062711
    etaiov_mtt_2_Study1 ~ fixed(0.062711)                                              # NONMEM $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance
    etaiov_mtt_1_Study2 ~ 0.072800                                                     # Pejcic 2024 Table 2: IOV(MTT_st2) = 27.48% CV -> omega^2 = log(1 + 0.2748^2) = 0.072800
    etaiov_mtt_2_Study2 ~ fixed(0.072800)                                              # NONMEM $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Methods section 2.2: "a study-specific proportional error model
    # described the intra-individual variability". Table 2 lists
    # exactly two residual-error rows -- Wp(st1) and Wp(st2) -- with no
    # analyte dimension, so a single proportional error is shared
    # between clopidogrel and CLO-CA within each study. This is
    # consistent with Figure 2, whose goodness-of-fit panels pool the
    # two analytes on one set of axes, and with a proportional error
    # being unitless (the two analytes are reported in different units,
    # ng/mL and ug/mL).

    propSdStudy1 <- 0.4195
    label("Proportional residual SD, Study 1 (fraction)")                              # Pejcic 2024 Table 2: Wp(st1) = 41.95% (RSE 3.4%; SIR median 41.98%)

    propSdStudy2 <- 0.2939
    label("Proportional residual SD, Study 2 (fraction)")                              # Pejcic 2024 Table 2: Wp(st2) = 29.39% (RSE 5.7%; SIR median 29.51%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Constants and covariate decomposition
    # ------------------------------------------------------------------
    # Molecular-weight correction applied to the CLO-CA formation flux
    # so that a mass of clopidogrel converts to the corresponding mass
    # of the (lighter) carboxylic acid.
    mwRatio <- 0.9565                    # Pejcic 2024 Methods section 2.2: "To account for the differences in molecular weights between CLO and CLO-CA (321.82 and 307.80, respectively), a scaling factor of 0.9565 was applied to the formation rate of CLO-CA." 307.80 / 321.82 = 0.95644.

    st2 <- STUDY_CLO_BE2                 # 1 = Study 2, 0 = Study 1
    st1 <- 1 - STUDY_CLO_BE2
    oc1 <- (OCC == 1)                    # crossover period 1
    oc2 <- (OCC == 2)                    # crossover period 2

    # ------------------------------------------------------------------
    # 2. Study-specific typical values and random effects
    # ------------------------------------------------------------------
    lmttTv    <- st1 * lmttStudy1    + st2 * lmttStudy2
    lfgenTv   <- st1 * lfgenStudy1   + st2 * lfgenStudy2
    logitfmTv <- st1 * logitfmStudy1 + st2 * logitfmStudy2

    etaFdepot <- st1 * etalfdepotStudy1 + st2 * etalfdepotStudy2
    etaFm     <- st1 * etalogitfmStudy1 + st2 * etalogitfmStudy2

    iovFdepot <- st1 * (oc1 * etaiov_fdepot_1_Study1 + oc2 * etaiov_fdepot_2_Study1) +
                 st2 * (oc1 * etaiov_fdepot_1_Study2 + oc2 * etaiov_fdepot_2_Study2)
    iovMtt    <- st1 * (oc1 * etaiov_mtt_1_Study1    + oc2 * etaiov_mtt_2_Study1) +
                 st2 * (oc1 * etaiov_mtt_1_Study2    + oc2 * etaiov_mtt_2_Study2)

    # ------------------------------------------------------------------
    # 3. Individual parameters (allometric scaling to 70 kg)
    # ------------------------------------------------------------------
    wtCl <- (WT / 70)^e_wt_cl
    wtV  <- (WT / 70)^e_wt_vc

    cl       <- exp(lcl) * wtCl
    vc       <- exp(lvc + etalvc) * wtV
    qh       <- exp(lqh)                 # Table 2 unit is "L/h", not "L/h/70 kg" -- deliberately unscaled
    vh       <- exp(lvh) * wtV
    mtt      <- exp(lmttTv + iovMtt)
    ktr      <- 3 / mtt                  # Savic 2007 (source reference 35): Ktr = (n_transit + 1) / MTT with n_transit = 2

    cl_cloca <- exp(lcl_cloca) * wtCl
    vc_cloca <- exp(lvc_cloca + etalvc_cloca) * wtV
    q_cloca  <- exp(lq_cloca) * wtCl
    vp_cloca <- exp(lvp_cloca) * wtV

    # ------------------------------------------------------------------
    # 4. First-pass branch fractions (Pejcic 2024 Equations 1-3)
    # ------------------------------------------------------------------
    # FaM is fixed, so FR2 is back-solved from Eq. 2 rather than being a
    # free parameter; using the individual FR1 keeps FaM at exactly 12%
    # for every subject, which is what "FaM was fixed to ... 12%" says.
    fr1  <- exp(logitfmTv + etaFm)
    fr2  <- fm_h4 * (1 + fr1) / (1 - fm_h4)
    fiam <- fr1 / (1 + fr1 + fr2)        # Eq. 1 -- fraction to CLO-CA
    fp   <- 1   / (1 + fr1 + fr2)        # Eq. 3 -- fraction escaping as intact parent

    # ------------------------------------------------------------------
    # 5. ODE system (Pejcic 2024 Figure 1)
    # ------------------------------------------------------------------
    # depot -> transit1 -> transit2 -> liver, each transfer at Ktr. The
    # liver receives the whole absorbed dose plus the parent returning
    # from the systemic circulation at CLP / Vc,P, and empties at
    # Qh / Vh into three branches (FP, FiaM, FaM) that sum to one. The
    # FaM branch is a pure sink: CLO-TH was not measured in either
    # study and is not carried as a state.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(liver)    <-  ktr * transit2 + cl / vc * central - qh / vh * liver
    d/dt(central)  <-  fp * qh / vh * liver - cl / vc * central

    d/dt(central_cloca) <- mwRatio * fiam * qh / vh * liver -
                           cl_cloca / vc_cloca * central_cloca -
                           q_cloca  / vc_cloca * central_cloca +
                           q_cloca  / vp_cloca * peripheral1_cloca
    d/dt(peripheral1_cloca) <- q_cloca / vc_cloca * central_cloca -
                               q_cloca / vp_cloca * peripheral1_cloca

    # ------------------------------------------------------------------
    # 6. Bioavailability
    # ------------------------------------------------------------------
    # The reference product carries the fixed F = 1 anchor; the generic
    # carries the study's estimated relative bioavailability. IIV and
    # IOV act on whichever product was given.
    frel     <- exp(lfdepot) * (1 - FORM_CLO_GENERIC) +
                exp(lfgenTv) * FORM_CLO_GENERIC
    f(depot) <- frel * exp(etaFdepot + iovFdepot)

    # ------------------------------------------------------------------
    # 7. Observations
    # ------------------------------------------------------------------
    # Amounts are in mg and volumes in L, so amount / volume is mg/L =
    # ug/mL. CLO-CA was reported in ug/mL and needs no conversion;
    # clopidogrel was reported in ng/mL and is multiplied by 1000.
    Cc       <- central / vc * 1000                  # ng/mL
    Cc_cloca <- central_cloca / vc_cloca             # ug/mL

    # One proportional residual error per study, shared by both analytes
    # (Table 2 has no analyte dimension on Wp).
    propSd       <- st1 * propSdStudy1 + st2 * propSdStudy2
    propSd_cloca <- st1 * propSdStudy1 + st2 * propSdStudy2

    Cc       ~ prop(propSd)
    Cc_cloca ~ prop(propSd_cloca)
  })
}
