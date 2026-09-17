Li_2017_pomalidomide_egfr <- function() {
  description <- "One-compartment population PK model with first-order absorption for orally administered pomalidomide in patients with relapsed or refractory multiple myeloma (rrMM) and various degrees of renal impairment (Li 2017, Table 5, eGFR column). Renal function enters as a CONTINUOUS covariate on apparent clearance through the paper's 'reverse-hockey-stick' relationship: CL/F rises linearly with MDRD-estimated glomerular filtration rate up to a breakpoint and is constant above it, so CL/F = cl_nonren + slope * min(CRCL, crcl_cap). The intercept of 3.96 L/h is the non-renal clearance and the shallow slope is the paper's evidence that pomalidomide clearance is insensitive to renal function. Interindividual variability on Ka, V/F, the intercept and the slope; residual error is additive on the log scale (log-normal). This model has NO absorption lag time -- Table 5 does not report one. Companion to Li_2017_pomalidomide_crcl, which is the same structure fitted with Cockcroft-Gault creatinine clearance as the renal marker."
  reference <- "Li Y, Wang X, O'Mara E, Dimopoulos MA, Sonneveld P, Weisel KC, Matous J, Siegel DS, Shah JJ, Kueenburg E, Sternas L, Cavanaugh C, Zaki M, Palmisano M, Zhou S. Population pharmacokinetics of pomalidomide in patients with relapsed or refractory multiple myeloma with various degrees of impaired renal function. Clin Pharmacol Adv Appl. 2017;9:133-145. doi:10.2147/CPAA.S144606"
  vignette <- "Li_2017_pomalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate from the Modification of Diet in Renal Disease (MDRD) equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "BSA-normalized (mL/min/1.73 m^2), the canonical CRCL normalisation. Li 2017 Equation 4: eGFR = 175 * (Scr,std)^-1.154 * (age)^-0.203 * (0.742 if female) * (1.212 if African American), with Scr,std the serum creatinine measured by a standardized assay. Cohort median 27 mL/min/1.73 m^2 (range 5-84), Table 2. Enters the reverse-hockey-stick clearance relationship of Equation 7 as an UNCENTERED linear term with an upper clamp at crcl_cap = 26.0 mL/min/1.73 m^2. NOTE: the female factor is printed as '0.724' in the Li 2017 Equation 4 typesetting, whereas the published MDRD-175 equation uses 0.742; the printed value appears to be a digit transposition in the source. Only the covariate DERIVATION is affected, never a model parameter, and downstream users supplying their own eGFR column are unaffected. The sibling model Li_2017_pomalidomide_crcl uses the same CRCL column to carry raw Cockcroft-Gault CrCl in mL/min, so the two are not interchangeable.",
      source_name = "eGFR"
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
    renal_function = "eGFR median 27 mL/min/1.73 m^2 (range 5-84), MDRD, Table 2.",
    co_medication = "low-dose dexamethasone",
    regions = "North America (USA, Canada) and Europe (France, Germany, Greece, Italy, Netherlands, Spain, UK, Austria)",
    notes = "Pooled intensive and sparse plasma sampling from studies CC-4047-MM-008 (n = 20) and CC-4047-MM-013 (n = 43). Baseline demographics are Li 2017 Table 2. Sex distribution and race/ethnicity are not reported in the source."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Li 2017 Table 5 ('final model with renal
    # function treated as a continuous variable'), 'eGFR as renal function
    # marker' Estimate column.
    #
    # Table 5 reports NO lag time, so this model has none -- the Tlag of
    # 0.154 h belongs to the categorical parameterisation (Table 3) only.
    # ------------------------------------------------------------------
    lka <- log(0.678); label("Absorption rate constant (1/h)") # Table 5 eGFR column: Ka = 0.678 1/h (90% bootstrap CI 0.526-0.948)
    lvc <- log(60.5); label("Apparent central volume of distribution V/F (L)") # Table 5 eGFR column: V/F = 60.5 L (90% bootstrap CI 54.502-69.225)

    # ------------------------------------------------------------------
    # Reverse-hockey-stick renal-function relationship. Li 2017
    # Equation 7:
    #   CL/F = intercept + slope * eGFR   (when eGFR <= eGFR0)
    #   CL/F = constant                   (when eGFR >  eGFR0)
    # The plateau 'constant' is not separately tabulated and no such
    # parameter appears in Table 5, so the relationship is the continuous
    # (hinged) one: the plateau is the value the linear arm reaches at the
    # breakpoint, 3.96 + 0.0483 * 26.0 = 5.216 L/h, which agrees with the
    # 5.13 L/h normal-renal-function CL/F of the categorical model
    # (Table 3) to within 2%.
    #
    # The intercept is the paper's non-renal clearance -- Discussion:
    # 'the intercept ... (3.71 and 3.96 L/h for CrCl and eGFR,
    # respectively) ... was approximately 70% of pomalidomide total
    # clearance (5.13 L/h)'.
    # ------------------------------------------------------------------
    lcl_nonren <- log(3.96); label("Non-renal apparent clearance, the intercept of the CL/F vs eGFR relationship (L/h)") # Table 5 eGFR column: Intercept = 3.96 L/h (90% bootstrap CI 2.577-5.086)
    e_crcl_cl_renal <- 0.0483; label("Slope of the renal clearance arm with respect to eGFR (L/h per mL/min/1.73 m^2)") # Table 5 eGFR column: Slope = 0.0483 (90% bootstrap CI 0.001-0.115)
    lcrcl_cap <- log(26.0); label("Breakpoint above which the CL/F vs eGFR relationship is flat (mL/min/1.73 m^2)") # Table 5 eGFR column: eGFR0 = 26.0 mL/min/1.73 m^2 (90% bootstrap CI 25.01-30.899)

    # ------------------------------------------------------------------
    # Interindividual variability. Li 2017 Equation 1 (P_i = P * exp(eta)),
    # applied to each parameter Table 5 gives an omega^2 for, including
    # the intercept and slope of the renal relationship. All values are
    # LOG-scale variances, transcribed without back-transformation. Table 5
    # reports four independent diagonal omegas and no covariances for the
    # continuous parameterisation.
    # ------------------------------------------------------------------
    etalka ~ 0.99 # Table 5 eGFR column: omega^2 (Ka) = 0.99 (90% bootstrap CI 0.64-1.46)
    etalvc ~ 0.0278 # Table 5 eGFR column: omega^2 (V/F) = 0.0278 (90% bootstrap CI 0-0.082)
    etalcl_nonren ~ 0.147 # Table 5 eGFR column: omega^2 (intercept) = 0.147 (90% bootstrap CI 0.062-0.297)
    etae_crcl_cl_renal ~ 0.516 # Table 5 eGFR column: omega^2 (slope) = 0.516 (90% bootstrap CI 0-3.837, very imprecise)

    # ------------------------------------------------------------------
    # Residual error. Li 2017 Equation 2: ln(C_obs) = ln(C_pred) + eps,
    # additive on the natural-log scale = nlmixr2's lnorm() error model.
    # Table 5 reports the variance delta^2 = 0.351, so the log-scale SD is
    # sqrt(0.351) = 0.5925.
    # ------------------------------------------------------------------
    expSd <- 0.5925; label("Log-scale residual standard deviation (log ng/mL)") # Table 5 eGFR column: delta^2 = 0.351 (90% bootstrap CI 0.221-0.492); sqrt(0.351) = 0.5925
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
