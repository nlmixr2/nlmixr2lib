AlZubaydi_2024_gabapentin <- function() {
  description <- paste(
    "One-compartment population PK model for gabapentin in hospitalized adults",
    "developed from therapeutic drug monitoring concentrations, with first-order",
    "absorption (ka fixed), saturable dose-dependent bioavailability, linear",
    "elimination, and serum creatinine on clearance"
  )
  reference <- paste(
    "Al-Zubaydi F, Wassef A, Kagan L, Brunetti L.",
    "Development of a Population Pharmacokinetic Gabapentin Model Leveraging",
    "Therapeutic Drug Monitoring Concentrations.",
    "Pharmaceutics. 2024;16(12):1514. doi:10.3390/pharmaceutics16121514"
  )
  vignette <- "AlZubaydi_2024_gabapentin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine on the day the gabapentin TDM sample was drawn",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power model on Cl, centered on the cohort MEAN serum creatinine of",
        "1.3 mg/dL (Al-Zubaydi 2024 Table 1: 1.3 +/- 1.0 mg/dL, range 0.4-3.8).",
        "Methods state that all continuous covariates were log-transformed and",
        "centered using mean values, which for a log-scale parameter is the",
        "power form cl = cl_pop * (CREAT / 1.3)^beta with beta = -0.89",
        "(Table 2, beta Cl_SCR). SCr gave the largest likelihood improvement of",
        "the four renal-function descriptors screened on Cl (Table S2:",
        "-2LL 452.85, delta -43.59 vs the 496.44 base model; eGFR -35.24,",
        "CrCl -20.70, AKI -13.83), and backward deletion of SCr worsened",
        "-2LL by 43.05 (Table S2 line 13), far above the 6.63 retention",
        "threshold. Table 1 quotes the laboratory",
        "reference interval as 0.66-1.25 mg/dL. Note the cohort mean 1.3 mg/dL",
        "is reported to only two significant figures, so the centering constant",
        "carries that rounding."
      ),
      source_name        = "SCr"
    )
  )

  # Covariates the paper screened but did NOT retain in the final model. These
  # are documentation only -- the negative findings for diabetes and obesity are
  # the paper's headline result (Abstract; Results Section 3.2; Discussion), so
  # the screen is worth preserving even though no point estimates were published
  # for any of these effects.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Actual body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested on Vd; not retained. Table S2 model 6: -2LL 456.83, i.e.",
        "delta +3.98 AGAINST the 452.85 covariate model, so adding WT did not",
        "improve the fit at all. Al-Zubaydi 2024 Results Section 3.2; Table 1",
        "mean 84.3 +/- 25.9 kg, range 44.2-195.0 kg."
      ),
      source_name = "WT"
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested on Vd; not retained (Table S2 model 5: -2LL 457.20,",
        "delta +4.35 against the covariate model). Table 1 mean",
        "59.9 +/- 10.3 kg, range 43.2-82.2 kg. Adjusted body weight was",
        "screened alongside it (Table S2 model 7: 456.41, delta +3.56, NS) but",
        "has no separate canonical column here since neither was retained."
      ),
      source_name = "IBW"
    ),
    LBM = list(
      description = "Lean body weight (Janmahasatian equation)",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested on Vd; not retained (Table S2 model 8: -2LL 457.44,",
        "delta +4.59 against the covariate model -- the worst of the five body",
        "size metrics screened). Table 1 mean 54.9 +/- 13.4 kg, range",
        "30.5-91.9 kg; Janmahasatian et al. 2005 formula (Al-Zubaydi 2024",
        "Methods Section 2.1, reference 29)."
      ),
      source_name = "LBW"
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Tested on Vd; not retained (Table S2 model 9: -2LL 457.24,",
        "delta +4.39 against the covariate model). Table 1 mean",
        "30.0 +/- 7.9 kg/m^2, range 18-60."
      ),
      source_name = "BMI"
    ),
    CRCL = list(
      description = paste(
        "Renal function estimate -- both Cockcroft-Gault creatinine clearance",
        "and MDRD eGFR were screened under this concept"
      ),
      units       = "mL/min (CrCl) or mL/min/1.73 m^2 (eGFR)",
      type        = "continuous",
      notes       = paste(
        "Both tested on Cl and both rejected in favour of raw serum creatinine,",
        "which gave the larger OFV reduction (Al-Zubaydi 2024 Results",
        "Section 3.2 and Discussion; Table S2 models 2 and 3: CrCl -2LL 475.74",
        "delta -20.70, eGFR -2LL 461.20 delta -35.24, versus SCr -2LL 452.85",
        "delta -43.59 -- all three significant at p<0.01, SCr best).",
        "Table 1: CrCl 91.67 +/- 67.6",
        "mL/min (range 13.3-192.7), computed by Cockcroft-Gault with the weight",
        "factor set to ideal body weight when actual weight was at least 30%",
        "above ideal; eGFR 71.0 +/- 42.4 mL/min/1.73 m^2 (range 16.0-127) by",
        "the MDRD equation."
      ),
      source_name = "CrCl / eGFR"
    ),
    FPG = list(
      description = "Fasting plasma glucose",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Tested on ka and Vd (main text Results Section 3.2); not retained",
        "(Table S2 model 11 tabulates FPG on Vd: -2LL 456.77, delta +3.92",
        "against the covariate model). Table 1 mean 127.6 +/- 78.2 mg/dL,",
        "range 76.0-251.0."
      ),
      source_name = "FPG"
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened among the continuous covariates (Methods Section 2.3);",
        "not retained. Table 1 mean 65.7 +/- 16.4 years, range 22-93."
      ),
      source_name = "Age"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened among the categorical covariates (Methods Section 2.3);",
        "not retained. Table 1: 51 of 82 (62.2%) female."
      ),
      source_name = "Sex"
    ),
    T2DM = list(
      description = "Diagnosis of type 2 diabetes (1 = yes)",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested on ka and Vd; not retained -- one of the paper's two primary",
        "negative findings (Abstract; Results Section 3.2; Conclusions).",
        "Table S2 model 10 tabulates DM on Vd: -2LL 456.83, delta +3.98",
        "against the covariate model, i.e. no improvement.",
        "Table 1: 26 of 82 (31.7%) with type 2 diabetes, determined from",
        "admission notes. Name is a documentation-only placeholder; it is not",
        "a ratified canonical in inst/references/covariate-columns.md because",
        "the paper publishes no point estimate for the effect."
      ),
      source_name = "Type 2 diabetes"
    ),
    OBESE = list(
      description = "Obesity status (1 = BMI at least 30 kg/m^2)",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested on Vd; not retained -- the paper's second primary negative",
        "finding. Table 1: 18 of 82 (22.0%) obese. Documentation-only",
        "placeholder name, as for T2DM."
      ),
      source_name = "Obesity"
    ),
    AKI = list(
      description = "Presence of acute kidney injury (1 = yes)",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested on Cl; not retained (serum creatinine was selected instead).",
        "Table S2 model 4: -2LL 482.61, delta -13.83, p<0.01 -- a real",
        "improvement over the base model, but the weakest of the four renal",
        "descriptors and superseded by SCr. Defined as a rise in serum",
        "creatinine of more than 0.3 mg/dL over the preceding 48 h.",
        "Table 1: 18 of 82 (21.9%). Documentation-only placeholder name."
      ),
      source_name = "AKI"
    ),
    OBESE_T2DM = list(
      description = paste(
        "Combined obesity/diabetes category (obese + diabetic, obese only,",
        "diabetic only, or neither)"
      ),
      units       = "(4-level factor)",
      type        = "categorical",
      notes       = paste(
        "Tested on ka and Vd; not retained (Table S2 model 12: -2LL 456.49,",
        "delta +3.64, flagged NS). Table 1 (cont.): obesity +",
        "diabetes 8 (9.8%), obesity only 22 (26.8%), 'healthy' 52 (63.4%).",
        "Documentation-only placeholder name."
      ),
      source_name = "Combined covariate of obesity and diabetes"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "gabapentin",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "gabapentin",
      units    = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 82,
    n_studies      = 1,
    age_range      = "22-93 years",
    age_mean       = "65.7 years",
    weight_range   = "44.2-195.0 kg",
    weight_mean    = "84.3 kg",
    sex_female_pct = 62.2,
    race_ethnicity = c(White = 89.0, Black = 7.3, Other = 2.4, Asian = 1.2),
    disease_state  = paste(
      "Hospitalized adults receiving oral gabapentin; the medical indication",
      "was not recorded in the dataset and the authors describe the cohort as",
      "a random general inpatient population"
    ),
    renal_function = paste(
      "Serum creatinine 1.3 +/- 1.0 mg/dL (range 0.4-3.8); MDRD eGFR",
      "71.0 +/- 42.4 mL/min/1.73 m^2 (range 16.0-127); Cockcroft-Gault CrCl",
      "91.67 +/- 67.6 mL/min (range 13.3-192.7); 18 of 82 (21.9%) with acute",
      "kidney injury within 48 h of sampling"
    ),
    dose_range     = paste(
      "Single oral dose 100-1200 mg (median 300 mg); total daily dose",
      "100-2700 mg (median 900 mg)"
    ),
    regions        = "United States (New Jersey)",
    notes          = paste(
      "Retrospective therapeutic drug monitoring cohort, Robert Wood Johnson",
      "University Hospital Somerset, 1 January 2009 to 7 December 2023",
      "(Methods Section 2.1). 123 gabapentin TDM concentrations from 108",
      "patients were available; 82 patients met the inclusion/exclusion",
      "criteria and contributed to the model (Results Section 3.1). The",
      "number of concentrations contributed by those 82 patients is not",
      "reported. Observed concentrations 7.97 +/- 7.8 ug/mL, range 0.6-56.2",
      "(Table 1); assay reference range 2.0-20.0 ug/mL with 25 ug/mL or above",
      "considered toxic (Methods Section 2.2). 31.7% type 2 diabetes, 22.0%",
      "obese (BMI at least 30 kg/m^2). Baseline demographics: Table 1.",
      "Estimation in MONOLIX 2020R1 (SAEM); final estimates and a",
      "1000-replicate bootstrap: Table 2. Structural and covariate model",
      "selection: Supplementary Tables S1 and S2. The selected structural",
      "model ('1st-order absorption with no Tlag and nonlinear F') had",
      "-2LL 496.44 and BICc 517.81, beating first-order absorption with a",
      "dose-independent F (507.42), the two Tlag variants (508.91, 499.57) and",
      "a transit-compartment model (505.39)."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters (Al-Zubaydi 2024 Table 2, "Fixed effects")
    # -----------------------------------------------------------------------
    # ka was NOT estimated: the TDM dataset carried too few absorption-phase
    # samples, so it was held at a literature value (Methods Section 2.3;
    # Table 2 reports "0.778 (fixed)").
    lka <- fixed(log(0.778)); label("Absorption rate constant ka (1/h)")            # Table 2: ka = 0.778 1/h (fixed); Methods Section 2.3 fixes it per references 16 and 32
    lcl <- log(5.73);         label("Clearance Cl (L/h)")                           # Table 2: Cl = 5.73 L/h (RSE 17.62%); bootstrap median 7.01, 95% CI 4.3-12.54
    lvc <- log(44.61);        label("Central volume of distribution Vd (L)")         # Table 2: Vd = 44.61 L (RSE 15.94%); bootstrap median 45.99, 95% CI 32.84-65.91

    # -----------------------------------------------------------------------
    # Saturable dose-dependent bioavailability, Equation 1:
    #     F = Dmax / (D50 + Dose)
    # Both constants were fixed to the previously published values of
    # Carlsson et al. 2009 (reference 20); neither was estimated here, so
    # neither appears in Table 2.
    # -----------------------------------------------------------------------
    ldmax_fdepot <- fixed(log(823));  label("Dmax: maximal absorbed gabapentin amount (mg)")                    # Equation 1 text: Dmax = 823, fixed to Carlsson 2009 (reference 20); Table S1 confirms "823 (fixed)" in both nonlinear-F columns; printed with a mg/day unit label -- see the vignette Errata note
    ld50_fdepot  <- fixed(log(1120)); label("D50: dose at which absorption is 50% saturated (mg)")              # Equation 1 text: D50 = 1120, fixed to Carlsson 2009 (reference 20); Table S1 confirms "1120 (fixed)"; same mg/day label caveat

    # -----------------------------------------------------------------------
    # Covariate effect
    # -----------------------------------------------------------------------
    e_creat_cl <- -0.89; label("Exponent for serum creatinine on Cl, centered at 1.3 mg/dL (unitless)")         # Table 2: beta Cl_SCR = -0.89 (RSE 16.94%); bootstrap median -1.12, 95% CI -1.66 to -0.58

    # -----------------------------------------------------------------------
    # Between-subject variability. MONOLIX reports lognormal random effects as
    # STANDARD DEVIATIONS on the log scale -- Table 2's footnote states
    # "omega: standard deviation" -- while nlmixr2 `ini()` takes VARIANCES, so
    # each value below is the published omega squared.
    # -----------------------------------------------------------------------
    etalcl ~ 0.0784; label("IIV on Cl (variance, log scale)")                       # Table 2: omega Cl = 0.28 SD (RSE 28.59%) -> 0.28^2 = 0.0784; approx 28.6% CV. Table S1 repeats the "omega: standard deviation" footnote
    etalvc ~ 0.5929; label("IIV on Vd (variance, log scale)")                       # Table 2: omega Vd = 0.77 SD (RSE 19.98%) -> 0.77^2 = 0.5929; approx 90% CV

    # -----------------------------------------------------------------------
    # Residual error: MONOLIX "constant" error model, y = f + a * eps, so the
    # reported a is an additive standard deviation on the concentration scale.
    # -----------------------------------------------------------------------
    addSd <- 2.03; label("Additive residual error (ug/mL)")                         # Table 2: a = 2.03 (RSE 17.5%); bootstrap median 1.93, 95% CI 1.1-2.65; "a: constant error model"
  })

  model({
    # ---- 1. Saturable dose-dependent bioavailability (Equation 1) ----
    # `podo(depot)` returns the amount of the dose currently being administered
    # BEFORE f(depot) scaling is applied, so F is evaluated self-consistently
    # at each dose event. Encoding the dose this way rather than as a separate
    # dose covariate column means the nonlinearity cannot be silently lost when
    # a simulation changes the dose. The paper's "Dose" is the last gabapentin
    # dose before the sample collection, i.e. the per-administration amount
    # (Equation 1 text; Table 1 footnote d).
    dmax_fdepot <- exp(ldmax_fdepot)
    d50_fdepot <- exp(ld50_fdepot)
    fdepot <- dmax_fdepot / (d50_fdepot + podo(depot))

    # ---- 2. Individual parameters ----
    ka <- exp(lka)
    # Log-transformed, mean-centered serum creatinine on Cl is the power form
    # (CREAT / 1.3)^e_creat_cl; 1.3 mg/dL is the cohort mean (Table 1).
    cl <- exp(lcl + etalcl) * (CREAT / 1.3)^e_creat_cl
    vc <- exp(lvc + etalvc)

    # ---- 3. Micro-constants ----
    kel <- cl / vc

    # ---- 4. ODE system: one compartment, first-order absorption, no lag ----
    # The final model has no lag time; Tlag was evaluated as a fixed 0.31 h but
    # the no-lag model was selected (Methods Section 2.3; Results Section 3.2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 5. Bioavailability ----
    f(depot) <- fdepot

    # ---- 6. Observation and error ----
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
