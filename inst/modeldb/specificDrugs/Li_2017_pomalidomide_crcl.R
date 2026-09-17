Li_2017_pomalidomide_crcl <- function() {
  description <- "One-compartment population PK model with first-order absorption for orally administered pomalidomide in patients with relapsed or refractory multiple myeloma (rrMM) and various degrees of renal impairment (Li 2017, Table 5, CrCl column). Renal function enters as a CONTINUOUS covariate on apparent clearance through the paper's 'reverse-hockey-stick' relationship: CL/F rises linearly with Cockcroft-Gault creatinine clearance up to a breakpoint and is constant above it, so CL/F = cl_nonren + slope * min(CRCL, crcl_cap). The intercept of 3.71 L/h is the non-renal clearance, about 70% of the 5.13 L/h total body clearance, and the shallow slope is the paper's evidence that pomalidomide clearance is insensitive to renal function. Interindividual variability on Ka, V/F, the intercept and the slope; residual error is additive on the log scale (log-normal). Unlike the categorical parameterisation (Li_2017_pomalidomide_renalcat) this model has NO absorption lag time -- Table 5 does not report one. See Li_2017_pomalidomide_egfr for the same structure fitted with MDRD eGFR as the renal marker."
  reference <- "Li Y, Wang X, O'Mara E, Dimopoulos MA, Sonneveld P, Weisel KC, Matous J, Siegel DS, Shah JJ, Kueenburg E, Sternas L, Cavanaugh C, Zaki M, Palmisano M, Zhou S. Population pharmacokinetics of pomalidomide in patients with relapsed or refractory multiple myeloma with various degrees of impaired renal function. Clin Pharmacol Adv Appl. 2017;9:133-145. doi:10.2147/CPAA.S144606"
  vignette <- "Li_2017_pomalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RAW, NOT BSA-normalized (mL/min, not mL/min/1.73 m^2). Li 2017 Equation 3: CrCl = [(140 - age) * body weight] / (72 * serum creatinine) * (0.85 for females), capped at a physiological ceiling of 150 mL/min. Cohort median 28.3 mL/min (range 8.7-115.4), Table 2. Enters the reverse-hockey-stick clearance relationship of Equation 7 as an UNCENTERED linear term with an upper clamp at crcl_cap = 37.7 mL/min; supplying a BSA-normalized value instead would silently rescale the renal arm. The sibling model Li_2017_pomalidomide_egfr uses the same CRCL column to carry MDRD eGFR in mL/min/1.73 m^2 with its own breakpoint, so the two are not interchangeable.",
      source_name = "CrCl"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pomalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pomalidomide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 63,
    n_studies = 2,
    age_range = "46-86 years",
    age_median = "69 years",
    weight_range = "39.9-116.6 kg",
    weight_median = "78.2 kg",
    disease_state = "relapsed or refractory multiple myeloma (rrMM) with normal, moderately impaired, or severely impaired renal function (with and without hemodialysis)",
    dose_range = "2-4 mg once daily, oral",
    renal_function = "CrCl median 28.3 mL/min (range 8.7-115.4), Cockcroft-Gault, Table 2.",
    co_medication = "low-dose dexamethasone",
    regions = "North America (USA, Canada) and Europe (France, Germany, Greece, Italy, Netherlands, Spain, UK, Austria)",
    notes = "Pooled intensive and sparse plasma sampling from studies CC-4047-MM-008 (n = 20) and CC-4047-MM-013 (n = 43). Baseline demographics are Li 2017 Table 2. Sex distribution and race/ethnicity are not reported in the source."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Li 2017 Table 5 ('final model with renal
    # function treated as a continuous variable'), 'CrCl as renal function
    # marker' Estimate column.
    #
    # Table 5 reports NO lag time, so this model has none -- the Tlag of
    # 0.154 h belongs to the categorical parameterisation (Table 3) only.
    # ------------------------------------------------------------------
    lka <- log(0.68); label("Absorption rate constant (1/h)") # Table 5 CrCl column: Ka = 0.68 1/h (bootstrap CI 0.522-0.931)
    lvc <- log(60.7); label("Apparent central volume of distribution V/F (L)") # Table 5 CrCl column: V/F = 60.7 L (bootstrap CI 54.393-68.717)

    # ------------------------------------------------------------------
    # Reverse-hockey-stick renal-function relationship. Li 2017
    # Equation 7:
    #   CL/F = intercept + slope * CrCl   (when CrCl <= CrCl0)
    #   CL/F = constant                   (when CrCl >  CrCl0)
    # The plateau 'constant' is not separately tabulated, and no such
    # parameter appears anywhere in Table 5, so the relationship is the
    # continuous (hinged) one: the plateau is the value the linear arm
    # reaches at the breakpoint, 3.71 + 0.0469 * 37.7 = 5.478 L/h. That
    # agrees with the 5.13 L/h normal-renal-function CL/F of the
    # categorical model (Table 3) to within 7%, which is the arithmetic
    # confirmation of the continuity reading.
    #
    # The intercept is the paper's non-renal clearance, hence lcl_nonren:
    # Discussion, 'the intercept of the total pomalidomide clearance
    # (CL/F) vs renal function marker (3.71 and 3.96 L/h for CrCl and
    # eGFR, respectively) curve was approximately 70% of pomalidomide
    # total clearance (5.13 L/h), suggesting that non-renal clearance
    # contributed roughly 70% of the total body clearance'.
    # ------------------------------------------------------------------
    lcl_nonren <- log(3.71); label("Non-renal apparent clearance, the intercept of the CL/F vs CrCl relationship (L/h)") # Table 5 CrCl column: Intercept = 3.71 L/h (bootstrap CI 2.321-4.796)
    e_crcl_cl_renal <- 0.0469; label("Slope of the renal clearance arm with respect to CrCl (L/h per mL/min)") # Table 5 CrCl column: Slope = 0.0469 (bootstrap CI 0.013-0.095)
    lcrcl_cap <- log(37.7); label("Breakpoint above which the CL/F vs CrCl relationship is flat (mL/min)") # Table 5 CrCl column: CrCl0 = 37.7 mL/min (bootstrap CI 35.027-50.343)

    # ------------------------------------------------------------------
    # Interindividual variability. Li 2017 Equation 1 (P_i = P * exp(eta)),
    # applied to every parameter that Table 5 gives an omega^2 for --
    # including the intercept and the slope of the renal relationship, so
    # both are log-normally distributed across subjects. All values are
    # LOG-scale variances, transcribed without back-transformation.
    #
    # Table 5 reports these as four independent diagonal omegas; no
    # covariances are given for the continuous parameterisation (unlike
    # Table 3, which carries a V/F-CL/F block).
    # ------------------------------------------------------------------
    etalka ~ 0.986 # Table 5 CrCl column: omega^2 (Ka) = 0.986 (bootstrap CI 0.638-1.467)
    etalvc ~ 0.0293 # Table 5 CrCl column: omega^2 (V/F) = 0.0293 (bootstrap CI 0-0.081)
    etalcl_nonren ~ 0.163 # Table 5 CrCl column: omega^2 (intercept) = 0.163 (bootstrap CI 0.076-0.383)
    etae_crcl_cl_renal ~ 0.3 # Table 5 CrCl column: omega^2 (slope) = 0.3 (bootstrap CI 0-0.813, lower bound at zero)

    # ------------------------------------------------------------------
    # Residual error. Li 2017 Equation 2: ln(C_obs) = ln(C_pred) + eps,
    # additive on the natural-log scale = nlmixr2's lnorm() error model.
    # Table 5 reports the variance delta^2 = 0.351, so the log-scale SD is
    # sqrt(0.351) = 0.5925.
    # ------------------------------------------------------------------
    expSd <- 0.5925; label("Log-scale residual standard deviation (log ng/mL)") # Table 5 CrCl column: delta^2 = 0.351 (bootstrap CI 0.217-0.471); sqrt(0.351) = 0.5925
  })

  model({
    # 1. Individual parameters.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)

    # 2. Reverse-hockey-stick clearance, Li 2017 Equation 7. Clamping the
    #    covariate axis at crcl_cap gives the linear arm below the
    #    breakpoint and the flat arm above it in a single expression, with
    #    the two arms meeting continuously at CRCL = crcl_cap.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    crcl_cap <- exp(lcrcl_cap)
    cl_renal <- e_crcl_cl_renal * exp(etae_crcl_cl_renal) * min(CRCL, crcl_cap)
    cl <- cl_nonren + cl_renal

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. One-compartment model with first-order oral absorption. No
    #    absorption lag time is reported for this parameterisation.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Observation, in the ng/mL Li 2017 reports (mg/L * 1000).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
