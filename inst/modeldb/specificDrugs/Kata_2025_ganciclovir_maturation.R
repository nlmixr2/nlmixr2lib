Kata_2025_ganciclovir_maturation <- function() {
  description <- paste(
    "Postnatal-age maturation model for ganciclovir (GCV) clearance in an",
    "extremely low birth weight preterm neonate with congenital",
    "cytomegalovirus infection (Kata 2025 J Pharm Health Care Sci). This is",
    "the paper's only originally-estimated model: the four phase-wise",
    "Bayesian clearance estimates, obtained by re-using the Acosta 2007",
    "neonatal ganciclovir population PK model without re-estimation (shipped",
    "here as Acosta_2007_ganciclovir), were normalised to a body surface area",
    "of 1.73 m^2 and regressed on postnatal age with a sigmoid hyperbolic",
    "(Hill) maturation function of the Anderson and Holford form,",
    "CL = cl_pna0 + cl_matspan * PNA^hill / (pna50^hill + PNA^hill).",
    "There are no compartments and no dosing: the rxode2 time axis IS",
    "postnatal age in days, and the model predicts body-surface-area",
    "normalised GCV clearance in mL/min/1.73m^2 at each observation time.",
    "The fitted curve rises from 8.80 mL/min/1.73m^2 at birth to an",
    "asymptote of 230.8, crossing the ~100 mL/min/1.73m^2 normal adult",
    "glomerular filtration rate late in the treatment course, which is the",
    "'catch-up' renal maturation the paper reports.",
    sep = " "
  )
  reference <- paste(
    "Kata K, Inomata S, Nishikawa M, Ide H, Nakamura K, Yoshida T,",
    "Taguchi M.",
    "Renal maturation and catch-up clearance of ganciclovir in a preterm",
    "neonate: Bayesian pharmacokinetic analysis using a population model.",
    "J Pharm Health Care Sci. 2025;11:89.",
    "doi:10.1186/s40780-025-00496-5.",
    "The sigmoid hyperbolic maturation form is cited by Kata 2025 to",
    "Anderson BJ, Holford NHG. Mechanism-based concepts of size and",
    "maturity in pharmacokinetics. Annu Rev Pharmacol Toxicol.",
    "2008;48:303-332.",
    sep = " "
  )
  vignette <- "Kata_2025_ganciclovir"
  units <- list(
    time          = "day (postnatal age, PNA)",
    dosing        = "n/a (no exogenous dosing; clearance-maturation regression)",
    concentration = "mL/min/1.73m^2 (body-surface-area normalised GCV clearance)"
  )

  # No covariateData: postnatal age is the model's time axis rather than a
  # covariate column, following the shipped structural sibling
  # inst/modeldb/endogenous/Wu_2024_gfr_maturation.R. Body weight and body
  # surface area enter only upstream, in the companion PK model and in the
  # 1.73 m^2 normalisation the authors applied before this regression.

  population <- list(
    species        = "human",
    n_subjects     = 1L,
    n_studies      = 1L,
    age_range      = "postnatal age 9 to 112 days during GCV/VGCV therapy",
    gestational_age = "27 weeks 0 days at birth",
    weight_range   = "0.530 to 2.010 kg during therapy (birth weight 0.556 kg)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported (single Japanese centre)",
    disease_state  = paste(
      "Congenital cytomegalovirus infection, confirmed by CMV DNA in urine",
      "in the early neonatal period, in an extremely low birth weight",
      "preterm neonate delivered by caesarean section for non-reassuring",
      "fetal status."
    ),
    dose_range     = "n/a (this model regresses clearance on postnatal age)",
    regions        = "Japan (Toyama University Hospital)",
    notes          = paste(
      "The regression has only four clearance observations, one per",
      "treatment phase, at approximately PNA 30, 51, 72 and 93 days",
      "(Table 1), so five parameters are estimated from four points plus the",
      "two 24-h-urine creatinine clearance reference values at PNA 6 and",
      "PNA 52 shown as open squares in Fig. 2. Confidence intervals were",
      "computed by the authors as estimate +/- 1.96 * SE (Methods,",
      "Maturation analysis). NONMEM 7.5.1, FOCE-I."
    )
  )

  ini({
    # ---- Maturation function parameters (Kata 2025 Table 2) ----
    # Kata 2025 Eq. 5, Methods "Maturation analysis", p. 3:
    #   CL = theta1 + theta2 * PNA^theta3 / (theta4^theta3 + PNA^theta3)
    # in mL/min/1.73m^2, with PNA the postnatal age in days. Table 2 gives
    # the estimates and 95% CIs; the Table 2 footnote marks theta1 and
    # theta2 as mL/min/1.73m^2 (footnote a) and theta4 as days (footnote b).
    lcl_pna0 <- log(8.80)
    label("Clearance intercept at postnatal age 0, theta1 (mL/min/1.73m^2)")  # Kata 2025 Table 2: theta1 = 8.80 (95% CI 8.26-9.34); prose "theta1 is the intercept of this curve, which corresponds to the GCV clearance at PNA 0"
    lcl_matspan <- log(222)
    label("Absolute clearance gained between birth and full maturation, theta2 (mL/min/1.73m^2)")  # Kata 2025 Table 2: theta2 = 222 (95% CI 169-274); Table 2 footnote calls it the scaling factor
    lhill <- log(2.79)
    label("Hill coefficient on postnatal age in the maturation function, theta3 (unitless)")  # Kata 2025 Table 2: theta3 = 2.79 (95% CI 2.43-3.15); prose "theta3 (Hill coefficient) describes the slope of maturation"
    lpna50 <- log(95.8)
    label("Postnatal age at half of the clearance maturation span, theta4 (days)")  # Kata 2025 Table 2: theta4 = 95.8 (95% CI 79.5-112); prose "theta4 describes the maturation half-time"; Results "the age at 50% maturity of GCV CL, was 95.8 (days)"

    # ---- Residual error ----
    # Kata 2025 Methods, Maturation analysis: "The variability was modeled
    # using a proportional error model incorporating a random variable eta,
    # which follows a normal distribution with mean 0 and variance omega^2."
    # Table 2 reports omega^2 = 0.0126, so the SD is sqrt(0.0126) = 0.1122
    # (11.2 %). The Table 2 FOOTNOTE instead calls omega^2 "an interindividual
    # variance"; with a single subject the two levels are not identifiable
    # from the source. Encoded at the residual level per the Methods sentence,
    # which is the one that actually describes the error structure. See the
    # vignette Errata.
    propSd <- sqrt(0.0126)
    label("Proportional residual SD on body-surface-area normalised clearance (fraction)")  # Kata 2025 Table 2: omega^2 = 0.0126 (95% CI 0.00817-0.0170); SD = sqrt(0.0126) = 0.1122
  })

  model({
    # ---- Time axis ----
    # The rxode2 time variable IS postnatal age in days, matching units$time
    # and the shipped sibling Wu_2024_gfr_maturation.R. Kata 2025 Eq. 5 is
    # written directly in PNA, with no reference age and no normalisation.
    pna_d <- time

    cl_pna0    <- exp(lcl_pna0)
    cl_matspan <- exp(lcl_matspan)
    hill       <- exp(lhill)
    pna50      <- exp(lpna50)

    # ---- Sigmoid hyperbolic (Hill) maturation, Kata 2025 Eq. 5 ----
    # Additive-on-the-absolute-scale Anderson and Holford form: the curve
    # runs from cl_pna0 at PNA 0 to cl_pna0 + cl_matspan asymptotically
    # (8.80 -> 230.8 mL/min/1.73m^2), reaching cl_pna0 + cl_matspan/2 at
    # PNA = pna50. Reproduces the Results: 17.2, 41.4, 77.8 and
    # 115.2 mL/min/1.73m^2 at PNA 30, 51, 72 and 93, i.e. the late-phase
    # value "exceeded 100 mL/min/1.73m^2, which approximates the normal
    # glomerular filtration rate".
    CL_bsa <- cl_pna0 +
      cl_matspan * pna_d^hill / (pna50^hill + pna_d^hill)

    CL_bsa ~ prop(propSd)
  })
}
