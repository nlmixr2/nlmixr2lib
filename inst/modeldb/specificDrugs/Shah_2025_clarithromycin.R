Shah_2025_clarithromycin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for clarithromycin in",
    "critically ill adults, with a priori allometric body-weight scaling on",
    "all disposition parameters. No covariate effects were retained and no",
    "auto-inhibition of clearance was detectable in this population."
  )
  reference <- paste(
    "Shah RV, Kipper K, Baker EH, Barker CIS, Oldfield I, Davidson HC,",
    "Swire CC, Philips BJ, Johnston A, Rhodes A, Sharland M, Standing JF,",
    "Lonsdale DO (2025).",
    "Intravenous Clarithromycin in Critically Ill Adults: A Population",
    "Pharmacokinetic Study. Antibiotics 14(6):559.",
    "doi:10.3390/antibiotics14060559.",
    sep = " "
  )
  vignette <- "Shah_2025_clarithromycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: Shah 2025 Section 4 states that total
  # clarithromycin concentrations were measured in plasma separated from
  # arterial blood samples by UHPLC-MS/MS; Section 2.2 selects a
  # two-compartment disposition model.
  compartmentData <- list(
    central     = list(analyte = "clarithromycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "clarithromycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate in the final model. Added a priori, not by",
        "stepwise selection (Shah 2025 Section 4: 'Weight was added a priori",
        "with an allometric exponent of 0.75 for clearance parameters and 1",
        "for compartment volumes'). Every structural theta in Table 2 is",
        "reported per 70 kg, which fixes the reference weight at 70 kg.",
        "Cohort median 80 kg (IQR 65-95, range 53-120; Table 1)."
      ),
      source_name        = "Weight (kg)"
    )
  )

  # Covariates screened by the authors but NOT retained in the final model, plus
  # the renal-replacement-therapy state that was handled by data exclusion
  # rather than by a covariate effect. Documented here for provenance only; none
  # is referenced in model(). Shah 2025 Section 2.2: "The addition of albumin,
  # creatinine, the presence of liver disease, sex, height and age to parameters
  # did not provide any significant improvement in the model fit."
  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Tested; did not significantly improve model fit (Shah 2025",
        "Section 2.2). Cohort median 25 g/L (IQR 21-29, range 12-38;",
        "Table 1) -- markedly hypoalbuminaemic, as expected in critical",
        "illness. Relevant to interpretation because clarithromycin is",
        "roughly 80% protein bound (Section 1), and the paper's own",
        "target-attainment simulations apply a fixed 80% bound fraction",
        "rather than an albumin-dependent one (Section 4)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Tested on clearance; did not significantly improve model fit",
        "(Shah 2025 Section 2.2). Cohort median 89 umol/L (IQR 71-121,",
        "range 40-276; Table 1). Contrast with the same group's",
        "benzylpenicillin model, where a serum-creatinine power effect on",
        "CL WAS retained -- see modellib('Shah_2023_benzylpenicillin')."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested; did not significantly improve model fit (Shah 2025",
        "Section 2.2). Cohort was 12 male : 7 female (Table 1). The paper",
        "reports the covariate only as 'sex' and does not state a reference",
        "level; because the effect was rejected, no coefficient or reference",
        "category needs to be resolved."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Tested; did not significantly improve model fit (Shah 2025",
        "Section 2.2). Cohort median 173 cm (IQR 166-178, range 150-192;",
        "Table 1). One participant had no height recorded (Table 1",
        "footnote), so height was also not available for the full cohort."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested; did not significantly improve model fit (Shah 2025",
        "Section 2.2). Cohort median 66 years (IQR 56.7-72.7, range",
        "25-85.8; Table 1)."
      )
    ),
    LIVER_DISEASE = list(
      description = "Presence of liver disease",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested; did not significantly improve model fit (Shah 2025",
        "Section 2.2). The paper gives no definition of the flag and no",
        "count of affected participants; Table 1 reports ALT (median 34",
        "U/L, IQR 24-48, range 9-166) but no liver-disease indicator.",
        "Documentation-only label -- LIVER_DISEASE is deliberately NOT",
        "registered in inst/references/covariate-columns.md because the",
        "final model does not use it and the source supplies no operational",
        "definition to register. It is NOT the same concept as the",
        "registered LIVER_RESECT_MAJOR canonical (major hepatic resection)."
      )
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Renal replacement therapy active",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Never tested as a covariate. Handled by data exclusion instead:",
        "Shah 2025 Section 4, 'Periods during which time participants were",
        "receiving renal replacement therapy were excluded from the",
        "analysis'; Section 2.1 reports that this removed 18 samples from 5",
        "participants (139 samples from 22 participants collected, 121",
        "samples from 19 participants analysed). The packaged model is",
        "therefore valid only OFF renal replacement therapy and must not be",
        "used to predict exposure during RRT."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19,
    n_studies      = 1,
    n_samples      = 121,
    age_range      = "25-85.8 years",
    age_median     = "66 years",
    weight_range   = "53-120 kg",
    weight_median  = "80 kg",
    height_range   = "150-192 cm",
    height_median  = "173 cm",
    sex_female_pct = 36.8,
    race_ethnicity = c(
      White           = 63.2,
      `Black British` = 5.3,
      Asian           = 5.3,
      `Not stated`    = 26.3
    ),
    disease_state  = paste(
      "Critical illness requiring intensive care. Infection sources were",
      "chest (14), ear/nose/throat (2), skin (1), central nervous system (1),",
      "gastrointestinal (1) and unknown (1); participants could have more",
      "than one source. Median APACHE II 20 points (IQR 16-23, range 0-28)"
    ),
    dose_range     = "500 mg IV every 12 h (all participants)",
    regions        = "United Kingdom (single centre: St George's Hospital, London)",
    renal_function = paste(
      "Serum creatinine median 89 umol/L (IQR 71-121, range 40-276).",
      "Periods of renal replacement therapy were excluded a priori: 18 of",
      "139 samples, from 5 of 22 participants"
    ),
    hepatic_function = "ALT median 34 U/L (IQR 24-48, range 9-166)",
    severity       = paste(
      "APACHE II median 20 points (IQR 16-23, range 0-28); 12 participants",
      "received vasopressors; 12 had periods of invasive ventilation, 4 of",
      "non-invasive ventilation and 14 of self-ventilation (categories",
      "overlap over time); 12 had at least one plasma pH below 7.35.",
      "90-day outcome: 13 alive, 3 infection-attributable deaths, 3 deaths",
      "not attributable to infection"
    ),
    co_medication  = paste(
      "Concomitant drugs were screened against the British National",
      "Formulary; none was predicted to affect clarithromycin PK",
      "(Section 2.1)"
    ),
    notes          = paste(
      "Sub-study of the ABDose observational antibiotic PK/PD study",
      "(REC 14/LO/1999), the same study that supplied",
      "modellib('Shah_2023_benzylpenicillin'). Baseline demographics:",
      "Shah 2025 Table 1. Sampling was opportunistic around the indicative",
      "schedule in Table 3 (0.5/1, 2, 6, 11.75 h in dosing interval 1 and",
      "0.5/1, 8, 10, 11.75 h in dosing interval 2), maximum 8 samples per",
      "participant; not every participant had a peak concentration measured.",
      "Clarithromycin was given as an infusion over 1-2 h in local practice",
      "(Section 4). Model fit in NONMEM 7.5.1. This is the first published",
      "description of intravenous clarithromycin PK in humans."
    )
  )

  ini({
    # Structural parameters: Shah 2025 Table 2, 'Fixed effects', reported as
    # mean parameter estimates for a 70 kg adult. No covariate effects other
    # than the a priori allometric weight scaling were retained.
    lcl <- log(8.17); label("Clearance (L/h/70 kg)")                       # Table 2, theta_CL = 8.17 L/h/70 kg (17% RSE; bootstrap median 7.8, 95% CI 5.6-10.7)
    lvc <- log(25.7); label("Central volume of distribution (L/70 kg)")    # Table 2, theta_V1 = 25.7 L/70 kg (29% RSE; bootstrap median 25.3, 95% CI 5.4-45.0)
    lq  <- log(62.0); label("Intercompartmental clearance (L/h/70 kg)")    # Table 2, theta_Q = 62.0 L/h/70 kg (18% RSE; bootstrap median 62.4, 95% CI 43.4-112.5)
    # Table 2 prints the theta_V2 row unit as "L/h/70 kg", but V2 is a volume
    # and the unit is a typographical slip in the header. Two independent
    # closed-form checks on the reported value 60.6 confirm it is litres:
    #   Vss = V1 + V2 = 25.7 + 60.6 = 86.3 L/70 kg, exactly the "volume of
    #     distribution at steady state in this study (86.3 L/70 kg)" quoted in
    #     Section 3; and
    #   the terminal half-life of the two-compartment system built from
    #     CL 8.17, V1 25.7, Q 62.0, V2 60.6 is log(2)/beta = 7.81 h, exactly the
    #     "Half life: 7.8 h" printed beneath Table 2.
    # Neither identity holds if 60.6 is read as a clearance.
    lvp <- log(60.6); label("Peripheral volume of distribution (L/70 kg)")  # Table 2, theta_V2 = 60.6 L/70 kg (16% RSE; bootstrap median 62.0, 95% CI 44.7-108.0)

    # Allometric exponents: fixed a priori, not estimated. Shah 2025 Section 4:
    # "Weight was added a priori with an allometric exponent of 0.75 for
    # clearance parameters and 1 for compartment volumes, as previously
    # described". One exponent is shared across the two clearance parameters
    # (CL, Q) and one across the two volumes (V1, V2), hence the paired
    # e_<cov>_<param1>_<param2> canonical form.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")   # Section 4, Materials and Methods
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")  # Section 4, Materials and Methods

    # IIV: Shah 2025 Table 2, 'OMEGA'. The row labels are "eta^2 CL" and
    # "eta^2 V1", i.e. the values are already log-scale VARIANCES and need no
    # %CV back-transformation. Cross-check against the individual-estimate
    # ranges in the same table confirms the variance reading: with omega =
    # sqrt(0.53) = 0.728, the reported CL range 2.1-55.7 L/h/70 kg about a
    # typical 8.17 spans -1.9 to +2.6 SD over 19 subjects; with omega =
    # sqrt(1.55) = 1.245, the reported V1 range 3.1-766.4 L/70 kg about a
    # typical 25.7 spans -1.7 to +2.7 SD. Reading either number as a %CV would
    # make the published individual ranges unreachable.
    # No IIV on Q or V2 -- Table 2 reports OMEGA terms for CL and V1 only.
    etalcl ~ 0.53  # Table 2, eta^2 CL = 0.53 (32% RSE; bootstrap median 0.51, 95% CI 0.21-0.87) -> 84% CV
    etalvc ~ 1.55  # Table 2, eta^2 V1 = 1.55 (70% RSE; bootstrap median 1.58, 95% CI 0.15-7.95) -> 193% CV

    # Residual error: Shah 2025 Table 2, 'SIGMA', reported as the sigma^2
    # variance of a proportional error model; nlmixr2 parameterises residual
    # error by standard deviation, so this is entered as sqrt(sigma^2).
    propSd <- sqrt(0.034); label("Proportional residual error (fraction)")  # Table 2, sigma^2 proportional = 0.034 (30% RSE; bootstrap median 0.0317, 95% CI 0.0161-0.0553) -> SD 0.1844
  })

  model({
    # Individual parameters. Reference weight 70 kg -- every Table 2 theta is
    # reported "/70 kg".
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment intravenous disposition (Shah 2025 Section 2.2: "A
    # two-compartment model was found to provide the best fit"; one-, two- and
    # three-compartment models were tested per Section 4). Section 2.2 also
    # records that a Michaelis-Menten elimination model did not improve the fit
    # and that interoccasion variability showed no evidence of auto-inhibition,
    # so elimination is encoded as first-order and time-invariant.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
