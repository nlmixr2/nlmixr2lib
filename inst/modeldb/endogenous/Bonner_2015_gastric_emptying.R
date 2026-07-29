Bonner_2015_gastric_emptying <- function() {
  description <- paste(
    "Population meta-analysis of the postprandial gastric emptying (GE)",
    "time course in humans, pooled across 49 published studies of",
    "1457 subjects spanning 28 weeks gestational age (very-low-birth-",
    "weight preterm neonates) through adults (Bonner 2015). Structural",
    "model is a double Weibull mixture on the fraction of test meal",
    "remaining in the stomach vs. time-after-ingestion",
    "(Eq. 1: GE = (100 - PR) * exp(-(t/gamma1)^beta1) +",
    "PR * exp(-(t/gamma2)^beta2)). The five test meal categories",
    "(aqueous, breast milk, formula, semi-solid, solid) are the only",
    "significant covariate; postnatal age and gestational age were",
    "tested and found NOT significant. The five reported theta_meal",
    "coefficients enter as multiplicative scalings on the fast-phase",
    "Weibull scale gamma1, with aqueous as the renormalised reference",
    "(theta_meal / theta_Aqueous). This attachment point, the natural-",
    "percent reading of PR, and the CV-percent reading of the Table 3",
    "omega^2 column are documented best-effort interpretations because",
    "the paper's Methods do not specify which parameter the meal-type",
    "theta modifies, the PR reporting scale, or the omega^2 unit",
    "convention, and neither the main paper nor the two supplementary",
    "material files contain the NONMEM control stream. See the",
    "validation vignette's Errata for the full ambiguity list and the",
    "operator-approved (2026-07-24) best-effort choices. Not a drug",
    "PK model -- no dose, no plasma compartment, no drug PK",
    "observation; the single state is the percentage of test meal",
    "remaining, initialised to 100 at t = 0 by construction. Reported",
    "meal-type mean gastric residence times (MGRT, Eq. 4) from the",
    "paper's 1000-individual simulations: 45 min (aqueous), 57 min",
    "(breast milk), 64 min (formula), 87 min (semi-solid), 98 min",
    "(solid)."
  )
  reference <- paste(
    "Bonner JJ, Vajjah P, Abduljalil K, Jamei M,",
    "Rostami-Hodjegan A, Tucker GT, Johnson TN.",
    "Does age affect gastric emptying time?",
    "A model-based meta-analysis of data from premature neonates",
    "through to adults. Biopharm Drug Dispos.",
    "2015 May;36(4):245-57. doi:10.1002/bdd.1937."
  )
  vignette <- "Bonner_2015_gastric_emptying"
  paper_specific_compartments <- c("gastric_remaining")

  units <- list(
    time          = "min",
    dosing        = "none (the % remaining is initialised to 100 at t = 0 by construction of Eq. 1; test meal is not represented as a dose)",
    concentration = "% of test meal remaining in the stomach (0 - 100)"
  )

  covariateData <- list(
    MEAL_AQUEOUS = list(
      description        = "Binary indicator for aqueous test meal (water, sugar solutions, fruit juice per Bonner 2015 Methods 'Covariate selection and evaluation'). 1 = aqueous test meal for the study record; 0 = otherwise. Selects theta_Aqueous in the meal-type mixture; because aqueous is the renormalised reference (theta_meal / theta_Aqueous), MEAL_AQUEOUS = 1 gives a meal-type ratio of exactly 1 on gamma1 for that record.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not aqueous)",
      notes              = "Bonner 2015 Table 3 theta_Aqueous = 0.697. Mutually exclusive with MEAL_BREASTMILK, MEAL_FORMULA, MEAL_SEMISOLID, MEAL_SOLID -- exactly one of the five MEAL_* indicators is 1 per record. See vignette Errata for the meal-type attachment-point best-effort choice (multiplies gamma1 with aqueous as the renormalised reference).",
      source_name        = "aqueous solution"
    ),
    MEAL_BREASTMILK = list(
      description        = "Binary indicator for breast-milk test meal. 1 = breast milk; 0 = otherwise. Multiplicative scaling factor on the fast-phase Weibull scale gamma1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not breast milk)",
      notes              = "Bonner 2015 Table 3 theta_Breast_milk = 0.959. Mutually exclusive with the other MEAL_* indicators.",
      source_name        = "breast milk"
    ),
    MEAL_FORMULA = list(
      description        = "Binary indicator for formula (any variety, including nutritional shakes per Methods).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not formula)",
      notes              = "Bonner 2015 Table 3 theta_Form = 1.15. Mutually exclusive with the other MEAL_* indicators.",
      source_name        = "formula"
    ),
    MEAL_SEMISOLID = list(
      description        = "Binary indicator for semi-solid meal (pudding, rice cereal, or oatmeal per Methods).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not semi-solid)",
      notes              = "Bonner 2015 Table 3 theta_Semi_solid = 1.61. Mutually exclusive with the other MEAL_* indicators.",
      source_name        = "semi-solid meals"
    ),
    MEAL_SOLID = list(
      description        = "Binary indicator for solid meal (e.g., pancakes, eggs, chicken liver, sandwich meals per Bonner 2015 supplementary Table 1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not solid)",
      notes              = "Bonner 2015 Table 3 theta_Solid = 1.99. Mutually exclusive with the other MEAL_* indicators.",
      source_name        = "solid meals"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description        = "Postnatal age in weeks. Tested as a covariate on GE (paper Methods 'Covariate selection and evaluation') after allowance for meal type and rejected because the objective function value did not change materially: BASE MODEL + Food types + postnatal Age OFV 1875.121 vs Food types alone OFV 1875.252 (Table 2), delta OFV = -0.13, not significant. Documented here so downstream users know age was screened.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Methods; not retained in the final model. Postmenstrual age (postnatal + gestational, in weeks) was used only for graphical purposes and not as a modelled covariate."
    ),
    GA = list(
      description        = "Gestational age at birth in weeks. Tested as a covariate on GE after allowance for meal type (paper Table 2: BASE MODEL + Food types + Gestational age OFV 1875.211 vs Food types alone 1875.252, delta OFV = -0.04). Rejected as non-significant.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Same rationale as PNA: screened per Methods, not retained. For preterm neonates gestational age at birth was added to postnatal age when computing postmenstrual age for graphical inspection."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1457L,
    n_studies      = 49L,
    age_range      = "28 weeks gestational age (VLBW preterm neonates) through adults up to ~84 years (1008 months). Median across paediatric studies is low (many neonatal / preterm cohorts); adult subjects (325 individuals) span ~17-84 years. See supplementary Table 1 for the per-study 'Age in months mean (range)' column.",
    weight_range   = "Not reported per-study in the pooled dataset. Population physically spans preterm neonates (<= ~1 kg) through adults.",
    sex_female_pct = "Not fully specified across the dataset. Of the 1457 modelling subjects, 637 were paediatrics with sex not specified; 495 paediatrics with sex specified (256 M / 228 F, plus 11 not-clearly-classified); 325 adults with sex specified (160 M / 165 F). The 165 adult women include both pre- and postmenopausal ages; most studies did not distinguish.",
    regions        = "Not reported per-study. 49 modelling studies drawn from PubMed / Embase (English-language, human) published 1975 - 2012.",
    dose_range     = "No drug dosing. Test meals administered orally, by nasogastric tube, or by orogastric tube; span 5 - 10% dextrose / sugar solutions, breast milk, term / preterm infant formula, rice cereal / pudding, pancakes, egg sandwiches, chicken liver, and other solids (see supplementary Table 1 for the per-study test-meal detail).",
    disease_state  = "Healthy preterm neonates through adults. Exclusions per Methods: obese subjects, subjects on GI-motility drugs (e.g., metoclopramide, cisapride), any disease other than apnoea, and subjects with confirmed gastro-oesophageal reflux (subjects referred for suspected GOR were retained).",
    notes          = "Model-development set: n = 1457 subjects across 49 studies. Independent validation set: n = 468 subjects across 17 studies (marked '*' in supplementary Table 1). Sampling times 20 - 300 min after meal ingestion; time-points-per-study ranged 1 - 19 (many single-time-point neonatal studies); subjects-per-study ranged 6 - 186. GE measurement methods pooled across scintigraphy (Tc-99m or In-111 -- residual-error weight T = 2), dilution (phenol red or polyethylene glycol), ultrasound, MRI, and applied potential tomography (T = 1). PR is estimated on a logit-transformed scale to keep the parameter in (0, 100). Fit in NONMEM 7.2 with FOCE-I."
  )

  ini({
    # ---------------------------------------------------------------------
    # Double Weibull structural parameters. Bonner 2015 Table 3 (p. 251).
    # gamma1 and gamma2 are Weibull scales in minutes; beta1 and beta2 are
    # unitless shapes. gamma1 < gamma2 is the paper's convention for the
    # "fast" and "slow" emptying phases, so gamma1 has the larger weight
    # (100 - PR) and gamma2 the smaller weight PR. PR is estimated on a
    # logit-transformed scale per Methods to constrain it to (0, 100);
    # llogit_pr = log(PR / (100 - PR)) with PR here read on the natural
    # percent scale (see Errata for the reporting-scale best-effort).
    # ---------------------------------------------------------------------
    llogit_pr <- log(0.26 / (100 - 0.26)); label("Logit of the fast-phase-completed remaining fraction PR/100 (unitless; logit(0.0026) = -5.951)")  # Table 3: PR (%) = 0.26 -- natural-percent scale (best-effort, response-003 2026-07-24; see Errata)
    lbeta1    <- log(0.816);              label("Log of the fast-phase Weibull shape beta1 (unitless)")                                              # Table 3: beta1 = 0.816
    lbeta2    <- log(2.48);               label("Log of the slow-phase Weibull shape beta2 (unitless)")                                              # Table 3: beta2 = 2.48
    lgamma1   <- log(37.6);               label("Log of the fast-phase Weibull scale gamma1 (min)")                                                  # Table 3: gamma1 = 37.6 min
    lgamma2   <- log(63.7);               label("Log of the slow-phase Weibull scale gamma2 (min)")                                                  # Table 3: gamma2 = 63.7 min

    # ---------------------------------------------------------------------
    # Meal-type coefficients. Bonner 2015 Table 3 reports five theta_meal
    # coefficient estimates. The paper does not state which Eq. 1
    # parameter the meal-type theta modifies. The operator-approved
    # (response-003, 2026-07-24) best-effort choice implemented here is
    # that theta_meal multiplies gamma1 with aqueous as the renormalised
    # reference (theta_meal_ratio = theta_meal / theta_Aqueous), so an
    # aqueous record leaves gamma1 unchanged and a solid record scales
    # gamma1 by 1.99 / 0.697 = 2.855. Alternative attachment points
    # (gamma2, both, or a shape parameter) are documented in the
    # validation vignette's Errata.
    # ---------------------------------------------------------------------
    ltheta_aq  <- log(0.697); label("Log of the aqueous test-meal coefficient theta_Aqueous (unitless)")            # Table 3: theta_Aqueous = 0.697
    ltheta_bm  <- log(0.959); label("Log of the breast-milk test-meal coefficient theta_Breast_milk (unitless)")    # Table 3: theta_Breast_milk = 0.959
    ltheta_fm  <- log(1.15);  label("Log of the formula test-meal coefficient theta_Form (unitless)")               # Table 3: theta_Form = 1.15
    ltheta_ss  <- log(1.61);  label("Log of the semi-solid test-meal coefficient theta_Semi_solid (unitless)")      # Table 3: theta_Semi_solid = 1.61
    ltheta_sol <- log(1.99);  label("Log of the solid test-meal coefficient theta_Solid (unitless)")                # Table 3: theta_Solid = 1.99

    # ---------------------------------------------------------------------
    # Between-study variability. Bonner 2015 Table 3 "Variability,
    # omega^2 (RSE)" column. The paper does not state the unit convention
    # for the omega^2 values (114, 38.6, 14.1, 58.7, 19.2), and the
    # NONMEM control stream is not on disk to disambiguate. The operator-
    # approved (response-003, 2026-07-24) best-effort interpretation
    # implemented here is that the reported values are CV% (percent
    # coefficient of variation of the corresponding structural
    # parameter), converted to log-scale variance via
    # omega^2_log = log(1 + (CV% / 100)^2). See vignette Errata for the
    # alternative readings (raw log-scale variance -- gives infeasibly
    # large BSV; or (CV%)^2 -- gives near-zero BSV).
    # ---------------------------------------------------------------------
    etallogit_pr ~ 0.833   # Table 3: omega^2 PR    = 114  -> CV = 114% -> log(1 + 1.14^2)   = 0.833
    etalbeta1    ~ 0.139   # Table 3: omega^2 beta1 = 38.6 -> CV = 38.6% -> log(1 + 0.386^2) = 0.139
    etalbeta2    ~ 0.0197  # Table 3: omega^2 beta2 = 14.1 -> CV = 14.1% -> log(1 + 0.141^2) = 0.0197
    etalgamma1   ~ 0.296   # Table 3: omega^2 gamma1 = 58.7 -> CV = 58.7% -> log(1 + 0.587^2) = 0.296
    etalgamma2   ~ 0.0362  # Table 3: omega^2 gamma2 = 19.2 -> CV = 19.2% -> log(1 + 0.192^2) = 0.0362

    # ---------------------------------------------------------------------
    # Residual error. Paper Eq. 3:
    #   GE_ij = GE_hat_ij + (theta_w / sqrt(N_ij * T_i)) * epsilon_ij
    # with Var(epsilon) fixed at 1, theta_w = 11.1 (in units of %
    # remaining), N_ij the sample size at time-point ij, and T_i the
    # test-type weight (2 for scintigraphy, 1 for other methods). The
    # per-record N_ij and T_i weighting is a fit-time construct that
    # nlmixr2's built-in error models do not natively encode; the
    # additive residual SD is therefore left at the unweighted theta_w
    # value in this implementation. See vignette Errata.
    # ---------------------------------------------------------------------
    addSd <- 11.1; label("Weighted-additive residual SD theta_w (% remaining); Eq. 3 N * T weighting NOT reproduced")  # Table 3: theta_w = 11.1
  })

  model({
    # ---------------------------------------------------------------------
    # Individual structural parameters. IIV enters on the log scale for
    # gamma / beta and on the logit scale for PR, matching paper Eq. 2
    # (log-normal exponential-error model for log-transformed parameters)
    # and the Methods statement that PR is estimated on a logit-
    # transformed scale.
    # ---------------------------------------------------------------------
    beta1_i    <- exp(lbeta1 + etalbeta1)
    beta2_i    <- exp(lbeta2 + etalbeta2)
    gamma1_i   <- exp(lgamma1 + etalgamma1)
    gamma2_i   <- exp(lgamma2 + etalgamma2)
    logit_pr_i <- llogit_pr + etallogit_pr
    pr_frac_i  <- 1 / (1 + exp(-logit_pr_i))
    pr_pct_i   <- 100 * pr_frac_i

    # ---------------------------------------------------------------------
    # Meal-type coefficient. Aqueous is the renormalised reference:
    # theta_meal_ratio = theta_meal / theta_Aqueous, so MEAL_AQUEOUS = 1
    # gives ratio = 1 by construction, and MEAL_SOLID = 1 gives ratio =
    # 1.99 / 0.697 = 2.855. Exactly one MEAL_* indicator is expected to
    # be 1 per record (the event table must enforce mutual exclusivity).
    # This attachment point is the operator-approved (response-003,
    # 2026-07-24) best-effort choice among the three plausible
    # placements (gamma1, gamma2, or both); see vignette Errata.
    # ---------------------------------------------------------------------
    theta_meal_raw <- exp(ltheta_aq)  * MEAL_AQUEOUS +
                      exp(ltheta_bm)  * MEAL_BREASTMILK +
                      exp(ltheta_fm)  * MEAL_FORMULA +
                      exp(ltheta_ss)  * MEAL_SEMISOLID +
                      exp(ltheta_sol) * MEAL_SOLID
    theta_meal_ratio <- theta_meal_raw / exp(ltheta_aq)
    gamma1_eff <- gamma1_i * theta_meal_ratio

    # ---------------------------------------------------------------------
    # Double Weibull observation (paper Eq. 1):
    #   GE(t) = (100 - PR) * exp(-(t / gamma1_eff)^beta1) +
    #                 PR   * exp(-(t / gamma2_i)  ^beta2)
    # A tiny regularisation offset on t handles the shape < 1 case
    # (paper beta1 = 0.816) at t = 0 in numerical simulation; at t = 0
    # the algebraic value is exactly 100 by construction.
    # ---------------------------------------------------------------------
    t_eff <- time + 1e-9
    gastric_remaining <- (100 - pr_pct_i) * exp(-(t_eff / gamma1_eff)^beta1_i) +
                pr_pct_i        * exp(-(t_eff / gamma2_i)  ^beta2_i)

    gastric_remaining ~ add(addSd)
  })
}
