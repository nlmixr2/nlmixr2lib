Peng_2025_meropenem <- function() {
  description <- "One-compartment IV population PK model for prolonged-infusion meropenem in 21 Chinese critically ill adults on continuous venovenous hemofiltration (Peng 2025). Total clearance is the sum of an estimated endogenous (body) clearance of 2.89 L/h and the individually measured CRRT clearance supplied as the data column QEFF; Cockcroft-Gault creatinine clearance acts on the endogenous arm through an exponential term centered at the 13.6 mL/min cohort median. Central volume 26.0 L with no inter-individual variability (fixed to zero by the authors). Age, sex, body weight, APACHE II, CRP, PCT, sepsis, renal-failure category and anuria were screened but not retained."
  reference <- "Peng Y, Liu Y, Cheng Z, Zhang Q, Xie F, Zhu S, Li S. Population pharmacokinetics of prolonged infusion for meropenem: tailoring dosing recommendations for Chinese critically ill patients on continuous renal replacement therapy with consideration for renal function. Drug Des Devel Ther. 2025;19:1105-1117. doi:10.2147/DDDT.S489603"
  vignette <- "Peng_2025_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation; raw mL/min, NOT BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Peng 2025 Table 1: median 13.6 mL/min (IQR 6.9-22.2) across the 21-subject cohort;",
        "Methods list it as 'Cockcroft-Gault-based creatinine clearance (CLCR)' among the screened",
        "continuous covariates. The 13.6 mL/min cohort median is the centering constant CLCR_median",
        "of Equation 1. Retained as the SOLE covariate, on the endogenous clearance arm only, via the",
        "exponential form exp(theta_CLCR * (CLCR - 13.6)); age and sex were significant univariately",
        "but dropped in backward elimination once CLCR was in the model (Results, p. 1109).",
        "The cohort is uniformly renally impaired -- 13 chronic renal failure, 7 AKI, 1 normal --",
        "so the fitted range is roughly 7-22 mL/min, while the paper's own dosing simulations",
        "extrapolate the term out to 10-50 mL/min. All subjects were on CRRT, so this column is the",
        "RESIDUAL native renal function on top of the extracorporeal clearance carried by QEFF."
      ),
      source_name        = "CLCR"
    ),
    QEFF = list(
      description        = "Individually determined CRRT (continuous venovenous hemofiltration) clearance of meropenem",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Peng 2025 Methods 'Quantification of Meropenem and Its CRRT Clearance': CL_CRRT is the",
        "ultrafiltrate flow rate Q_uf (mL/h) multiplied by the sieving coefficient Sc, where Sc is the",
        "effluent-to-plasma meropenem concentration ratio measured at the same time point.",
        "Supplied to the model as a per-subject data column rather than estimated, and ADDED to the",
        "endogenous clearance in Equation 1. Table 1 cohort medians: ultrafiltrate 2477.5 mL/h",
        "(IQR 2406.6-2559.0), Sc 0.75 (IQR 0.72-0.88); the median CRRT dose was 25.76 mL/h/kg.",
        "The paper's Monte Carlo simulations standardise every virtual subject to a 25 mL/h/kg CRRT",
        "dose at the 65 kg cohort-median weight, i.e. Q_uf = 1625 mL/h, which at the median Sc of 0.75",
        "gives QEFF = 1.219 L/h. Note the mL/h -> L/h conversion: a weight-normalized CRRT",
        "prescription must be multiplied by body weight AND by Sc, then divided by 1000, before it",
        "enters this column."
      ),
      source_name        = "CL_CRRT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = "Peng 2025 Results p. 1109: significant on CLbody in the univariable screen, but 'age and sex did not remain significant covariates for CLbody after accounting for creatinine clearance'. Table 1 median 58.0 years (IQR 54.0-71.0)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Peng 2025 Results p. 1109: significant on CLbody univariately, dropped in backward elimination once CLCR was included. Table 1: 16 male (76.2%), 5 female (23.8%)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Peng 2025 Methods: screened as a continuous covariate; not retained. Table 1 median 65.0 kg (IQR 60.0-70.0). The 65 kg median is nevertheless load-bearing OUTSIDE the model: the paper's simulations convert the weight-normalized 25 mL/h/kg CRRT dose into the absolute ultrafiltrate flow that feeds QEFF."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation (APACHE) II score",
      units       = NA_character_,
      type        = "continuous",
      notes       = "Peng 2025 Methods: screened; not retained. Table 1 median 19 (IQR 15-26)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Peng 2025 Methods: screened; not retained. Cohort values not tabulated; the paper reports only that C-reactive protein was among the collected laboratory examinations."
    ),
    PCT = list(
      description = "Procalcitonin",
      units       = "ng/mL",
      type        = "continuous",
      notes       = "Peng 2025 Methods: screened; not retained. Cohort values not tabulated. Documentation only -- procalcitonin has no entry in inst/references/covariate-columns.md, and none is proposed here because the covariate was not retained and the paper reports no values or effect estimate for it; the units above are the conventional clinical ones, not a source-paper statement."
    ),
    DIS_SEPSIS = list(
      description = "Sepsis (versus other infection) indicator",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Peng 2025 Methods: 'sepsis or other infections' screened as a categorical covariate; not retained. Table 1: sepsis / septic shock 5 (23.8%), severe pneumonia 5 (23.8%), uremia 3 (14.3%), other 8 (38.1%)."
    ),
    DIS_ARF = list(
      description = "Acute kidney injury (versus chronic renal failure) indicator",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Peng 2025 Methods: 'AKI or chronic renal failure' screened as a categorical covariate; not retained. Table 1: chronic renal failure 13 (61.9%), AKI 7 (33.3%), normal 1 (4.8%)."
    ),
    ANURIA = list(
      description = "Anuric (versus non-anuric) urine-output indicator",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Peng 2025 Methods: 'urine output (anuric or non-anuric)' screened as a categorical covariate; not retained. Per-category counts not tabulated; the paper reports only that patients had 'minimal residual renal function' (Limitations, p. 1115)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 21L,
    n_studies        = 1L,
    age_median       = "58.0 years (IQR 54.0-71.0)",
    weight_median    = "65.0 kg (IQR 60.0-70.0)",
    sex_female_pct   = 23.8,
    race_ethnicity   = "Chinese (single-center cohort, Changsha, Hunan)",
    disease_state    = "Adults in a respiratory intensive care unit requiring continuous renal replacement therapy (CVVH) and receiving meropenem as standard antimicrobial therapy. Median APACHE II 19 (IQR 15-26). Main diagnoses: sepsis / septic shock 5 (23.8%), severe pneumonia 5 (23.8%), uremia 3 (14.3%), other including renal failure 8 (38.1%). Patients with major surgery in the preceding 4 weeks, or who were pregnant, were excluded.",
    renal_function   = "Chronic renal failure 13 (61.9%), acute kidney injury 7 (33.3%), normal 1 (4.8%). Median Cockcroft-Gault creatinine clearance 13.6 mL/min (IQR 6.9-22.2). All subjects on CVVH; median ultrafiltrate flow 2477.5 mL/h (IQR 2406.6-2559.0), median meropenem sieving coefficient 0.75 (IQR 0.72-0.88), median CRRT dose 25.76 mL/h/kg.",
    dose_range       = "Meropenem 1 g intravenously every 8-12 h as a 2-3 h prolonged infusion. Observed regimens: 1 g q8h 10 (47.6%), 1 g q12h 8 (38.1%), 1 g q6h 2 (9.5%), 0.5 g q8h 1 (4.8%).",
    regions          = "China (Third Xiangya Hospital, Central South University, Changsha). Enrolled May 2021 to April 2023.",
    n_concentrations = 94L,
    notes            = "Prospective single-center study. Plasma and CRRT-effluent meropenem measured by validated UPLC-PDA. NONMEM 7.5 FOCE-I with PsN 5.3 and Pirana 3.0. Model evaluated by prediction-corrected VPC (1000 replicates) and a 1000-sample non-parametric bootstrap (Table 2). Monte Carlo dosing simulations used mrgsolve 1.4.1 with 10,000 virtual subjects per scenario, all fixed at 65 kg and a 25 mL/h/kg CRRT dose."
  )

  ini({
    # Structural parameters: Peng 2025 Table 2, "Final Model Estimate (RSE%)"
    # column. Equation 1 (p. 1109) is the complete final structural model:
    #
    #   CL_total (L/h) = theta_CLbody
    #                    * exp(theta_CLCR * (CLCR - CLCR_median))
    #                    * exp(eta_CLbody)
    #                    + CL_CRRT
    #
    # CLCR_median = 13.6 mL/min, the Table 1 cohort median, so theta_CLbody is
    # the endogenous clearance of the MEDIAN-renal-function subject and not of
    # an anuric one. CL_CRRT is measured per subject (Q_uf * Sc) and enters as
    # the data column QEFF -- it is not estimated.
    # Table 2: theta_CLbody = 2.89 L/h (RSE 10.4%; bootstrap median 2.88, 95% CI 2.39-3.39)
    lcl <- log(2.89); label("Endogenous (body) clearance at CRCL = 13.6 mL/min (L/h)")
    # Table 2: theta_V = 26.0 L (RSE 12.7%; bootstrap median 26.78, 95% CI 19.51-32.41)
    lvc <- log(26.0); label("Central volume of distribution (L)")

    # Covariate effect. Exponential (log-linear) form, applied to the
    # endogenous clearance arm only and centered on the cohort median:
    #   cl <- exp(lcl) * exp(e_crcl_cl * (CRCL - 13.6))
    # Units of the coefficient are per mL/min. Positive sign: better residual
    # renal function raises endogenous clearance.
    # Table 2: theta_CLCR = 0.0183 (RSE 19.1%; bootstrap median 0.016, 95% CI 0.0065-0.030)
    e_crcl_cl <- 0.0183; label("Exponential coefficient of CRCL on endogenous clearance (per mL/min)")

    # Inter-individual variability. Peng 2025 Results p. 1109: "incorporating
    # log-normal inter-individual variability in endogenous clearance".
    # Table 2 reports IIV as %CV with the footnote formula
    #   CV(%) = sqrt(exp(omega^2) - 1) * 100
    # so omega^2 = log(CV^2 + 1). Final model eta_CLbody = 42.1 %CV
    # (RSE 37.4%, shrinkage 6.7%; bootstrap median 40.2, 95% CI 21.4-59.0).
    etalcl ~ 0.163214  # log(0.421^2 + 1)

    # No eta on the central volume. Peng 2025 Results p. 1109: "Due to
    # limitations in our dataset, we were unable to precisely estimate IIV on
    # the distribution volume (V). Therefore, to maintain the numerical
    # stability of the model, the variability of this population parameter was
    # fixed at zero, resulting in an increase of 9.94 points in the OFV."
    # Table 2's "eta_V (%CV)" row is a dash in every column. Encoded by
    # OMITTING the random effect rather than by a zero-variance eta, which is
    # numerically identical and avoids a singular OMEGA at solve time.

    # Combined proportional + additive residual error (Peng 2025 Results
    # p. 1109: "Residual variability was effectively captured by a combined
    # proportional and additive error model"; Table 2 final-model column).
    # The proportional row is a %CV on the residual SD scale, NOT the
    # log-normal transform used for the IIV rows: the bootstrap column prints
    # a median of 26.3(%) alongside a 95% CI of 0.018-0.103, which is the CI on
    # the VARIANCE, and sqrt(0.0692) = 0.263 recovers the median exactly.
    # Table 2: proportional error 28.6 %CV (RSE 18.0%, shrinkage 8.0%)
    propSd <- 0.286; label("Proportional residual error (fraction)")
    # Table 2: additive error 0.128 mg/L (RSE 9.8%, shrinkage 7.9%; bootstrap median 0.132)
    addSd <- 0.128; label("Additive residual error (mg/L)")
  })
  model({
    # Endogenous ("body") clearance: the estimated arm of Equation 1. The
    # exponential covariate term is centered on the 13.6 mL/min Table 1 cohort
    # median, so a median-renal-function subject has cl_body = exp(lcl) =
    # 2.89 L/h.
    cl_body <- exp(lcl + e_crcl_cl * (CRCL - 13.6) + etalcl)

    # Total clearance, Equation 1: the measured CRRT clearance is ADDED to the
    # endogenous arm. QEFF is a per-subject data column in L/h (Q_uf * Sc), so
    # no unit conversion happens here; convert mL/h to L/h on ingestion.
    #
    # The TOTAL clearance -- not the endogenous arm -- must be the variable
    # named `cl`. rxode2 recognises a `cl` / `vc` pair in a one-compartment
    # model and solves the system ANALYTICALLY from those two variables,
    # silently discarding the explicit `d/dt(central)` below. Writing the
    # endogenous arm as `cl` and the sum as, say, `cl_total` therefore
    # produces a model that reports a correct `cl_total` column while
    # simulating as though QEFF were zero: at CRCL = 25 and QEFF = 1.219 L/h
    # that understates the elimination rate constant by 26% (0.137 vs 0.184
    # per h) with no warning of any kind.
    cl <- cl_body + QEFF

    # Central volume carries no random effect (see ini()).
    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment model with linear elimination. Peng 2025 Discussion
    # p. 1113 notes the deviation from the commonly reported two-compartment
    # meropenem model: "the extended infusion method in our patient cohort
    # potentially obscured the distribution phase and hindered the capture of a
    # peripheral compartment."
    d/dt(central) <- -kel * central

    # Dose in mg and vc in L give mg/L, matching the UPLC-PDA plasma
    # concentrations the model was fitted to. The paper's PK/PD targets are
    # written on free concentrations; meropenem plasma protein binding is
    # negligible (about 2%) and the paper applies no unbound-fraction
    # correction, so Cc is used directly against the MIC targets.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
