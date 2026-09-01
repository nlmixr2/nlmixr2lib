Assmus_2025_benznidazole_qpcr <- function() {
  description <- paste(
    "Exposure-response model linking cumulative benznidazole exposure to",
    "post-treatment Trypanosoma cruzi qPCR positivity in adults with",
    "chronic indeterminate Chagas disease (Assmus 2025; BENDITA trial",
    "NCT03378661, modified intention-to-treat population n = 201",
    "including placebo). The population PK layer is the transit-absorption",
    "one-compartment model of Assmus 2025 Table 2 (see",
    "modellib('Assmus_2025_benznidazole')), extended with a cumulative",
    "AUC state. The pharmacodynamic layer is the paper's univariate beta",
    "binomial regression on qPCR positivity, the per-patient proportion of",
    "qPCR-positive blood samples collected after the end of treatment:",
    "logit(p) = b0 + b1 * AUCinf + b2 * Ct, where Ct is the baseline",
    "qPCR cycle threshold (higher Ct means lower parasite density). The",
    "exposure odds ratio is 0.9995 per mg*h/L and the Ct odds ratio is",
    "0.9021 per cycle. Note the paper's own conclusion: this association",
    "is driven almost entirely by the placebo-versus-active contrast, and",
    "no significant pharmacokinetic driver of response remained after",
    "placebo was excluded. The beta binomial overdispersion parameter was",
    "not reported, so the pharmacodynamic layer is deterministic."
  )
  reference <- paste(
    "Assmus F, Cruz C, Watson JA, White NJ, Adehin A, Hoglund RM,",
    "Blum de Oliveira B, Barreira F, Scandale I, Tarning J.",
    "Population pharmacokinetic-pharmacodynamic analysis of",
    "benznidazole monotherapy and combination therapy with",
    "fosravuconazole in chronic Chagas disease (BENDITA).",
    "PLoS Neglected Tropical Diseases. 2025;19(9):e0013522.",
    "doi:10.1371/journal.pntd.0013522.",
    "PK layer shared with modellib('Assmus_2025_benznidazole').",
    sep = " "
  )
  vignette <- "Assmus_2025_benznidazole"

  # Bookkeeping state that integrates Cc so the exposure-response layer
  # can read AUCinf off the solve. Same idiom as auc_wk8_12 in
  # Kuchimanchi_2018_evolocumab_ldlc.R and auc_free in
  # Lallemand_2023_benzylpenicillin_horse.R; not a biological compartment.
  paper_specific_compartments <- c("auc_central")

  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling applied a priori on CL/F (exponent 0.75) and V/F (exponent 1.0), standardized to the 65 kg cohort median (Assmus 2025 Methods, Population pharmacokinetic analysis (i)). Enrolment was restricted to 50-80 kg.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) is the register-canonical reference, but Assmus 2025 reports typical values for a FEMALE subject; the male effect is applied via (1 - SEXF) to preserve the published structural estimates verbatim.",
      notes              = "Retained on relative oral bioavailability F only: men had 12.9% lower F than women (Assmus 2025 Table 2). PK-layer covariate; sex was not a predictor in the exposure-response layer.",
      source_name        = "SEX"
    ),
    CONMED_FOSRAVUCONAZOLE = list(
      description        = "Concomitant fosravuconazole (E1224) administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (benznidazole monotherapy)",
      notes              = paste(
        "PK-layer covariate: +17.7% on CL/F (Assmus 2025 Table 2).",
        "Assmus 2025 also adjusted the exposure-response subgroup",
        "analysis for combination therapy and found no significant",
        "effect on qPCR positivity (Results, Exposure-parasitological",
        "response model), so no PD-layer coefficient is reported and",
        "none is encoded here."
      ),
      source_name        = "E1224"
    ),
    CT_TCRUZI_BASE = list(
      description        = "Baseline Trypanosoma cruzi qPCR cycle threshold (Ct)",
      units              = "(cycles)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mean cycle threshold of nine measurements (three blood samples",
        "assayed in triplicate) at the screening visit, before dosing.",
        "The scale is INVERSE to parasite burden: a higher Ct means",
        "fewer amplification cycles' worth of template, i.e. a LOWER",
        "T. cruzi density. Samples with no amplification are reported at",
        "the assay ceiling of 40. mITT median 37.7 (range 30.2-40.0;",
        "Assmus 2025 Table 1). The negative logit coefficient therefore",
        "means that patients with a HIGHER baseline parasite load (lower",
        "Ct) had a higher probability of post-treatment qPCR positivity",
        "(Assmus 2025 Results, Exposure-parasitological response model)."
      ),
      source_name        = "Ct"
    )
  )

  # Alternative exposure metrics that Assmus 2025 screened as univariate
  # predictors of qPCR positivity (Table 3) but that are not encoded in
  # model(). AUCinf is carried as the model's exposure driver because it
  # is the paper's leading metric (Table 3 column 2, Fig 8a) and the only
  # one that follows directly from integrating the PK layer. The others
  # are either post-hoc NCA summaries or trial-design quantities.
  covariatesDataExcluded <- list(
    CMAX_BZN = list(
      description        = "Individual peak benznidazole concentration in dried blood spots",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Alternative univariate predictor (Assmus 2025 Table 3: OR 0.640 [0.581, 0.706], p < 0.001 including placebo; OR 0.966 [0.863, 1.073], p = 0.529 excluding placebo). Derivable from the PK layer by simulation but not encoded as a state; see the validation vignette.",
      source_name        = "CMAX"
    ),
    T_ABOVE_TARGET = list(
      description        = "Time with benznidazole above a putative target concentration",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Alternative univariate predictor evaluated at 3 mg/L and 6 mg/L in plasma (2.5 and 5 mg/L in DBS) and at the in vitro IC90,DBS of 7.61 mg/L (29.3 uM). Assmus 2025 Table 3 carries the two plasma thresholds: at 3 mg/L, OR 0.883 [0.851, 0.915], p < 0.001 including placebo and OR 0.983 [0.954, 1.008], p = 0.217 excluding placebo; at 6 mg/L, OR 0.946 [0.920, 0.973], p < 0.001 including placebo and OR 0.991 [0.967, 1.011], p = 0.408 excluding placebo. Per S5 Fig, parameter estimates and goodness of fit were highly sensitive to the threshold, and the in vivo target remains unknown.",
      source_name        = "T>target"
    ),
    DUR_BZN_WEEKS = list(
      description        = "Weeks of benznidazole treatment (a week counts if at least one dose was taken)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Trial-design duration metric, not a PK quantity. Assmus 2025 Table 3: OR 0.449 [0.366, 0.550], p < 0.001 including placebo; OR 0.854 [0.709, 1.009], p = 0.077 excluding placebo. Of all metrics screened, duration showed the strongest residual trend after placebo exclusion, prompting the paper's conclusion that time, not only exposure, may drive response.",
      source_name        = "weeks"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "benznidazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "benznidazole", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "benznidazole", units = "mg",
      specimen = "whole blood", verified = TRUE
    ),
    auc_central = list(
      analyte = "benznidazole", units = "mg*h/L",
      specimen = "not applicable", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 201L,
    n_studies      = 1L,
    age_range      = "18-50 years",
    age_median     = "33 years",
    weight_range   = "50-80 kg",
    weight_median  = "64.5 kg",
    sex_female_pct = 70.1,
    race_ethnicity = "All participants were Bolivian.",
    disease_state  = paste(
      "Adults (18-50 years, 50-80 kg) with chronic indeterminate Chagas",
      "disease, confirmed by serological testing and a positive",
      "qualitative PCR result. Subjects with signs or symptoms of the",
      "chronic cardiac or digestive form were excluded."
    ),
    dose_range     = paste(
      "Six active oral benznidazole arms plus placebo: 150 mg twice",
      "daily for 8, 4 or 2 weeks; 150 mg once daily for 4 weeks alone or",
      "with fosravuconazole; and 300 mg once weekly (split into two",
      "doses) for 8 weeks with fosravuconazole."
    ),
    regions        = "Bolivia (Cochabamba, Tarija and Sucre).",
    endpoint       = paste(
      "qPCR positivity: the per-patient proportion of T. cruzi",
      "qPCR-positive blood samples collected after the end of treatment",
      "(EOT, the assigned treatment duration plus a two-week grace",
      "period) over 12 months of follow-up. Post-EOT visits fell at",
      "weeks 2, 3, 4, 6, 10 and 12 and at 4, 6 and 12 months, with the",
      "number of eligible visits varying by arm. A visit counted as",
      "positive if at least one of three 5 mL blood samples, each",
      "assayed in triplicate, tested positive."
    ),
    notes          = paste(
      "The PD layer was fitted in R 4.2.2 with glmmTMB using a beta",
      "binomial distribution parameterised per Morris 1997; the",
      "estimated overdispersion coefficient was not reported. The PK",
      "layer (n = 175 actively treated subjects) is unchanged from",
      "modellib('Assmus_2025_benznidazole'). PD population is the mITT",
      "set (n = 201, including 30 placebo subjects), which excludes 5",
      "subjects with possible treatment misallocation and 4 with no",
      "post-EOT follow-up. Assmus 2025 also reports a subgroup fit",
      "excluding placebo and one subject who took only 4 doses (n = 170,",
      "ordinary binomial with dispersion fixed to 1) in which the",
      "exposure effect was no longer significant (AUCinf OR 0.9999",
      "[0.9998, 1.0000], p = 0.282), and a sensitivity analysis in a",
      "refined mITT set (n = 186; S3 Table). Demographics from Assmus",
      "2025 Table 1 (PK/PD mITT column)."
    )
  )

  ini({
    # ---- Pharmacokinetic layer -------------------------------------
    # Identical to modellib("Assmus_2025_benznidazole"); Assmus 2025
    # Table 2, reported for a female adult weighing 65 kg. All values
    # are apparent (CL/F, V/F) and reference dried blood spot (whole
    # blood) concentrations.
    lmtt <- log(0.75);  label("Mean absorption transit time MTT (h)")            # Assmus 2025 Table 2: MTT 0.75 h (RSE 9.1%; bootstrap 95% CI 0.62-0.89)
    lcl  <- log(1.30);  label("Apparent oral clearance CL/F (L/h)")              # Assmus 2025 Table 2: CL/F 1.30 L/h (RSE 2.7%; bootstrap 95% CI 1.24-1.38)
    lvc  <- log(31.6);  label("Apparent central volume of distribution V/F (L)") # Assmus 2025 Table 2: V/F 31.6 L (RSE 2.9%; bootstrap 95% CI 29.9-33.5)

    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (fraction)")  # Assmus 2025 Table 2: "Relative oral bioavailability, F: 1 fixed"

    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/65) for CL/F (unitless)")  # Assmus 2025 Methods: allometric function on clearance (exponent 0.75), standardized to 65 kg
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on (WT/65) for V/F (unitless)")   # Assmus 2025 Methods: "... and volume (exponent 1) parameters"

    e_conmed_fosravuconazole_cl <- 0.177;  label("Fractional increase in CL/F with concomitant fosravuconazole (unitless)")  # Assmus 2025 Table 2: "Co-administration of E1224 on CL/F (%) 17.7" (RSE 28.9%; bootstrap 95% CI 8.18-27.3)
    e_sexf_fdepot               <- -0.129; label("Fractional change in F in men (SEXF = 0) relative to the female reference (unitless)")  # Assmus 2025 Table 2: "Sex effect on F (reference: female) (%) -12.9" (RSE 23.0%; bootstrap 95% CI -18.6 to -6.9)

    etalmtt    ~ log(1 + 0.616^2)  # Assmus 2025 Table 2: IIV MTT 61.6% CV (RSE 15.6%; bootstrap 95% CI 35.0-75.3)
    etalcl     ~ log(1 + 0.189^2)  # Assmus 2025 Table 2: IIV CL/F 18.9% CV (RSE 16.8%; bootstrap 95% CI 12.6-25.3)
    etalfdepot ~ log(1 + 0.102^2)  # Assmus 2025 Table 2: IIV F 10.2% CV (RSE 31.9%; bootstrap 95% CI 2.97-15.5)

    expSd <- sqrt(0.076); label("Residual error SD on the log-transformed concentration scale (log-normal)")  # Assmus 2025 Table 2: "Variance of residual error, sigma 0.076" (RSE 9.3%; bootstrap 95% CI 0.051-0.106)

    # ---- Exposure-response layer -----------------------------------
    # Assmus 2025 Equation 2: logit(p) = b0 + b1 * x + b2 * Ct, fitted by
    # beta binomial regression in the mITT population (n = 201, placebo
    # included) with x = AUCinf. Table 3 reports the two slopes as odds
    # ratios per unit predictor, so they enter as their natural logs.
    e_auc_logit <- log(0.9995)
    label("Additive shift on the qPCR-positivity logit per mg*h/L of cumulative AUCinf (logit units / (mg*h/L))")
    # Assmus 2025 Table 3, mITT including placebo, AUCinf column: OR 0.9995 [0.9994, 0.9997], p < 0.001.

    e_ct_tcruzi_base_logit <- log(0.9021)
    label("Additive shift on the qPCR-positivity logit per baseline Ct cycle (logit units / cycle)")
    # Assmus 2025 Table 3, mITT including placebo, AUCinf column: OR for the Ct value 0.9021 [0.8024, 1.0142], p = 0.085.

    # FIGURE-DERIVED (not reported in the paper's text or tables). Assmus
    # 2025 reports the two slopes of Equation 2 as odds ratios but never
    # reports the intercept b0, without which the absolute probability
    # level is not reproducible. It was recovered by digitising the solid
    # "Including placebo" fitted curve of Assmus 2025 Fig 8a (the caption
    # calls it grey; it is the darker of the two fitted lines, distinct
    # from the black placebo-arm observed points): the plotted
    # linear predictor at AUCinf = 0 and the median baseline Ct is
    # -0.388, i.e. p = 0.404 for an untreated patient of median parasite
    # load. The digitisation is self-validating -- refitting the traced
    # curve freely returns an exposure odds ratio of 0.99955 against the
    # published 0.9995, agreeing to four significant figures, which
    # confirms the axis calibration. b0 is then recovered by removing the
    # Ct contribution at the mITT median Ct of 37.7 (Assmus 2025 Table 1).
    # Uncertainty is roughly +/- 0.05 logit units (p(0) between 0.40 and
    # 0.42 depending on the fitting window). Fig 8 does not state which Ct
    # value the plotted curves use; the median is assumed, consistent with
    # Assmus 2025 Methods, which states that predictions "utilis[e] the
    # median Ct value at baseline".
    b0_qpcr <- -0.388 - log(0.9021) * 37.7
    label("Logit intercept for qPCR positivity at zero benznidazole exposure and Ct = 0 (logit units)")
  })

  model({
    # ---- Pharmacokinetic layer -------------------------------------
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl *
      (1 + e_conmed_fosravuconazole_cl * CONMED_FOSRAVUCONAZOLE)
    vc  <- exp(lvc) * (WT / 65)^e_wt_vc

    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_sexf_fdepot * (1 - SEXF))

    # Assmus 2025 Fig 2: ktr = (1 + n) / MTT with n = 1 transit
    # compartment, so gut -> transit1 -> central are two sequential
    # transfers at rate ktr and the mean transit time is 2 / ktr = MTT.
    ktr <- (1 + 1) / mtt
    kel <- cl / vc

    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1
    d/dt(central)  <-  ktr * transit1 - kel * central

    f(depot) <- fdepot

    Cc <- central / vc

    # ---- Exposure accumulation -------------------------------------
    # Assmus 2025 Methods, Definition of exposure metrics: AUCinf is the
    # area under the benznidazole blood concentration-time curve to
    # infinity for the whole treatment course. Integrating Cc from time
    # zero and reading the state out once absorption and elimination are
    # complete reproduces it; Results confirm the metric is cumulative
    # over the course ("approximate 2-fold increase when comparing 2
    # weeks to 4 weeks, and 4 weeks to 8 weeks").
    d/dt(auc_central) <- Cc

    # ---- Exposure-response layer -----------------------------------
    # Assmus 2025 Equation 2. p is the probability that a single
    # post-treatment follow-up sample is qPCR positive. Reported as a
    # deterministic quantity: the beta binomial overdispersion parameter
    # that governs the spread of the per-patient proportions around p was
    # estimated but never reported, so no residual-variability term can
    # be encoded without inventing one. The paper's Equation 3,
    # p* = 1 - (1 - p)^N for the probability of at least one positive
    # result across N follow-up visits, is a pure post-processing step on
    # p and is applied in the validation vignette rather than here (N is
    # a per-patient visit count, not a model state).
    pqpcr <- expit(b0_qpcr + e_auc_logit * auc_central +
                     e_ct_tcruzi_base_logit * CT_TCRUZI_BASE)

    # Benznidazole concentration in dried blood spots (mg/L). This is the
    # only endpoint with a published residual-error model; pqpcr above is
    # returned as a derived output column.
    Cc ~ lnorm(expSd)
  })
}
