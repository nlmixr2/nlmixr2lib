Zhao_2024_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral tacrolimus in adult renal transplant recipients co-administered voriconazole, built from within-15-days-post-transplant therapeutic drug monitoring. The absorption rate constant is fixed and the residual error is additive because the data are almost entirely troughs. Apparent clearance and apparent volume of distribution both fall as the measured voriconazole concentration rises (CYP3A4 inhibition), and apparent clearance also falls as serum creatinine rises."
  reference <- paste(
    "Zhao Y-C, Sun Z-H, Li J-K, Liu H-Y, Zhang B-K, Xie X-B, Fang C-H,",
    "Sandaradura I, Peng F-H, Yan M. Individualized dosing parameters for",
    "tacrolimus in the presence of voriconazole: a real-world PopPK study.",
    "Front Pharmacol. 2024;15:1439232. doi:10.3389/fphar.2024.1439232.",
    "Parameter values are the final-model estimates in Table 3 ('Final model'",
    "block), cross-checked against the bootstrap summary in Table 4.",
    sep = " "
  )
  vignette <- "Zhao_2024_tacrolimus"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Zhao 2024 Methods section 2.2 states that tacrolimus was
  # measured in WHOLE-BLOOD samples by chemiluminescence microparticle
  # immunoassay, so V/F and CL/F are whole-blood apparent parameters. (The
  # Discussion's passing phrase "the population's typical V/F, derived from
  # plasma concentration data" contradicts the Methods and is not followed
  # here; a whole-blood assay is also what makes V/F = 2690 L comparable to the
  # other whole-blood tacrolimus models tabulated in Zhao 2024 Table 7.)
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CONC_VORI_NGML = list(
      description        = "Measured voriconazole concentration in the co-administered patient (the source's C_VRC); 0 for a patient not receiving voriconazole",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING within subject: voriconazole concentrations were measured alongside every tacrolimus concentration (Zhao 2024 inclusion criterion 5, 'at least three measurements of tacrolimus and voriconazole concentrations'). UNIT CONVERSION: Zhao 2024 reports C_VRC in ug/mL (Table 1 median 0.00 [0.00, 0.50]; Discussion gives the observed range as 0-3.38 ug/mL) whereas the canonical column is ng/mL, so the model divides by 1000 before applying the published coefficients. The covariate is EXPONENTIAL, not power-scaled: the cohort median is exactly 0.00 ug/mL (many patients contributed pre-voriconazole samples) so no ratio-to-median form is possible, and exp(theta * C_VRC) reduces to 1 at C_VRC = 0, which is the untreated reference. Zhao 2024 warns that the model should not be extrapolated far above the observed 3.38 ug/mL ceiling: the published dose-recommendation tables run to 7 ug/mL, and the exponential term amplifies steeply there. Distinct from the binary CONMED_VORICONAZOLE, which is the right column for a model that estimates an on/off coefficient instead of a concentration-driven one.",
      source_name        = "CVRC"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING within subject; collected with every tacrolimus concentration (Zhao 2024 section 2.2). Enters as the power scaling (CREAT / 237)^-0.40, so apparent clearance FALLS as serum creatinine rises. 237 umol/L is the cohort median (Zhao 2024 Table 1, IQR 162.9-648.0) and is the value the authors held fixed for the voriconazole dose-recommendation simulation in Table 5, which is what identifies it as the normalising constant; the paper never prints the covariate equation itself. The cohort is drawn from the first 15 post-operative days, when graft function is still recovering, which is why the creatinine range extends far above the normal adult range. See the vignette's Assumptions and deviations section for the paper's internal contradiction over whether this coefficient acts on CL/F or on V/F.",
      source_name        = "CREA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19L,
    n_studies      = 1L,
    age_range      = "median 44 years (IQR 37.5-52.5)",
    weight_range   = "median 63 kg (IQR 51-72)",
    sex_female_pct = 21.1,
    race_ethnicity = c(Chinese = 100.0),
    disease_state  = "Adults (>= 18 years) within 15 days of a first renal transplant, on an oral triple immunosuppressive regimen (tacrolimus + mycophenolate mofetil + glucocorticoid) and receiving concomitant voriconazole for invasive fungal infection",
    dose_range     = "Oral tacrolimus, clinician-titrated; observed dose per administration median 3.00 mg (IQR 1.50-3.50)",
    regions        = "Single centre, The Second Xiangya Hospital of Central South University, Changsha, Hunan, People's Republic of China",
    notes          = "Retrospective non-interventional cohort, January 2016 to March 2021. 167 tacrolimus whole-blood concentrations (8-9 per patient), the majority drawn within 30 min before a dose, i.e. troughs. Observed tacrolimus concentration median 7.90 ng/mL (IQR 5.55-10.78); time after operation median 8 days (IQR 4-11); voriconazole concentration median 0.00 ug/mL (IQR 0.00-0.50), observed range 0-3.38 ug/mL; serum creatinine median 237 umol/L (IQR 162.9-648.0); albumin median 34.0 g/L; haematocrit median 26.2%. Genotypes were an inclusion requirement but neither was retained: CYP3A5 *1/*3 10 (52.6%) and *3/*3 9 (47.4%) with no *1/*1 carriers, CYP2C19 *1/*1 6 and *1/*2 9. Renal source DBD 17, DCD 1, living 1. Estimated in Phoenix NLME 8.1; covariates selected by forward addition (p <= 0.01) and backward elimination (p < 0.001); validated by 1000-sample bootstrap and prediction-corrected VPC. Registered as ChiCTR2100048712."
  )

  ini({
    # -------------------------------------------------------------------------
    # Structural parameters. Zhao 2024 Table 3, 'Final model' block (repeated
    # with bootstrap summaries in Table 4).
    #
    # Almost every sample is a pre-dose trough, so the data carry no absorption
    # information and Ka was fixed: Table 3 marks it "Ka (Fixed)" in the base
    # model with Stderr 0.00 and a degenerate 8.39-8.39 confidence interval,
    # and the Discussion confirms "the Ka value was fixed during the model's
    # development, which restricted our ability to explore the impact of other
    # covariates on Ka".
    #
    # The abstract calls 8.39 /h the "elimination rate constant (Ka)" and
    # section 3.2.1 repeats that wording; both are slips. Table 7, which
    # compares this model against six published tacrolimus popPK models, lists
    # the same 8.39 /h under the heading whose footnote reads "Ka, absorption
    # rate constant", alongside Ka values of 0.419-13.1 /h from the comparator
    # models. An elimination rate constant of 8.39 /h would give tacrolimus a
    # 5-minute half-life, which is irreconcilable with the reported
    # CL/F = 42.87 L/h and V/F = 2690 L (kel = 0.0159 /h, t1/2 = 43.5 h).
    # -------------------------------------------------------------------------
    lka <- fixed(log(8.39)); label("Absorption rate constant Ka (1/h)")                          # Zhao 2024 Table 3 final model: Ka = 8.39 1/h, fixed (Stderr 0.00, CI 8.39-8.39)
    lvc <- log(2690);        label("Apparent volume of distribution V/F at zero voriconazole (L)")  # Zhao 2024 Table 3 final model: V/F = 2690 L (bootstrap mean 2655, 95% CI 1480-4060)
    lcl <- log(42.87);       label("Apparent clearance CL/F at zero voriconazole and CREAT = 237 umol/L (L/h)")  # Zhao 2024 Table 3 final model: CL/F = 42.87 L/h (bootstrap mean 42.00, 95% CI 30-50)

    # -------------------------------------------------------------------------
    # Covariate effects. Zhao 2024 never prints the covariate equations; the
    # only description of their functional form is the Table 3 and Table 4
    # footnotes, which call each coefficient an "exponent for <covariate>, as a
    # covariate for <parameter>". The forms below are the only ones consistent
    # with that wording, with Phoenix NLME's covariate parameterisations, and
    # with the covariate distributions:
    #
    #   * C_VRC has a cohort median of exactly 0.00 ug/mL, so no
    #     (C_VRC / median)^theta form exists. The coefficients are applied
    #     exponentially, exp(theta * C_VRC), which is 1 at C_VRC = 0.
    #   * CREAT is strictly positive with median 237 umol/L, and the authors
    #     held it at "the median of 237 umol/L" for the Table 5 simulation,
    #     which identifies 237 as the normalising constant of a power term.
    #
    # Both signs are negative, matching the Results and Conclusion: higher
    # voriconazole concentration lowers CL/F and V/F, and higher serum
    # creatinine lowers CL/F.
    # -------------------------------------------------------------------------
    e_conc_vori_ngml_cl <- -0.28; label("Exponential coefficient of voriconazole concentration (ug/mL) on CL/F (unitless)")  # Zhao 2024 Table 3/4: Theta VRC-CL = -0.28 (Stderr 0.03; bootstrap mean -0.34, 95% CI -0.99 to -0.21)
    e_conc_vori_ngml_vc <- -0.20; label("Exponential coefficient of voriconazole concentration (ug/mL) on V/F (unitless)")   # Zhao 2024 Table 3/4: Theta VRC-V  = -0.20 (Stderr 0.04; bootstrap mean -0.20, 95% CI -0.36 to -0.03)
    e_creat_cl          <- -0.40; label("Power exponent of (CREAT / 237 umol/L) on CL/F (unitless)")                          # Zhao 2024 Table 3/4: Theta CREA = -0.40 (Stderr 0.11; bootstrap mean -0.40, 95% CI -0.72 to -0.05); placed on CL/F, see vignette Assumptions

    # -------------------------------------------------------------------------
    # Inter-individual variability. Zhao 2024 Table 4 reports these as
    # "omega^2 V" and "omega^2 CL", i.e. already on the variance scale of the
    # log-normal random effects, so they are used directly with no CV%
    # conversion. (The Table 4 "CV%" column for these two rows, 6.67 and 6.32,
    # is the relative standard error of the estimate, not the CV of the
    # parameter distribution: the variances themselves imply CV of 14% and 42%.)
    # -------------------------------------------------------------------------
    etalvc ~ 0.02  # Zhao 2024 Table 4: omega^2 V  = 0.02 (RSE 6.67%)
    etalcl ~ 0.16  # Zhao 2024 Table 4: omega^2 CL = 0.16 (RSE 6.32%)

    # -------------------------------------------------------------------------
    # Residual error. Section 3.2.1: "a one-compartment model with first-order
    # absorption and elimination, along with an additive residual model, was
    # selected". Table 4 reports a single sigma = 3.50 for the final model;
    # Phoenix NLME reports the additive residual as a standard deviation, and
    # 3.50 ng/mL against a median observed concentration of 7.90 ng/mL is the
    # large residual the Discussion comments on ("the residual variability
    # indicated by the sigma value in Table 4 suggest that factors beyond
    # voriconazole concentration and CREA may contribute to the observed inter-
    # and intra-individual variability").
    # -------------------------------------------------------------------------
    addSd <- 3.50; label("Additive residual error for tacrolimus whole-blood concentration (ng/mL)")  # Zhao 2024 Table 4: sigma = 3.50 (bootstrap mean 3.41, 95% CI 3.06-3.74)
  })

  model({
    # -----------------------------------------------------------------------
    # The published coefficients are on the ug/mL scale the paper reports
    # voriconazole in; the canonical covariate column is ng/mL.
    # -----------------------------------------------------------------------
    cvrc <- CONC_VORI_NGML / 1000

    ka <- exp(lka)
    vc <- exp(lvc + etalvc) * exp(e_conc_vori_ngml_vc * cvrc)
    cl <- exp(lcl + etalcl) *
      exp(e_conc_vori_ngml_cl * cvrc) *
      (CREAT / 237)^e_creat_cl

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # -----------------------------------------------------------------------
    # Observation. Doses are in mg and V/F is in L, so central / vc is in
    # mg/L; tacrolimus whole-blood concentrations are reported in ng/mL and
    # 1 mg/L = 1000 ng/mL.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc

    Cc ~ add(addSd)
  })
}
