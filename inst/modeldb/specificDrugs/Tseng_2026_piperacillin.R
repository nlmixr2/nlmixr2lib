Tseng_2026_piperacillin <- function() {
  description <- paste(
    "One-compartment population PK model for piperacillin in adults with low",
    "body weight (BMI <= 18.5 kg/m^2) receiving piperacillin-tazobactam",
    "(Tseng 2026). Clearance is scaled by a power function of the 2021 CKD-EPI",
    "creatinine-cystatin C estimated glomerular filtration rate; body weight was",
    "screened and not retained. Interindividual variability is carried on",
    "clearance only, the volume random effect having been removed for 83.2%",
    "shrinkage. A fixed unbound fraction of 0.7 converts the modelled total",
    "plasma concentration to the free concentration that drives the paper's",
    "100% fT > MIC target-attainment simulations."
  )
  reference <- paste(
    "Tseng YJ, Juan L, Tai CH, Wu CC. Optimal piperacillin/tazobactam dosing in",
    "adults with low body weight: a population pharmacokinetic and",
    "simulation-based study. Drug Des Devel Ther. 2026;20:1-12.",
    "doi:10.2147/DDDT.S602835. PMCID: PMC13157127."
  )
  vignette <- "Tseng_2026_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate from the 2021 CKD-EPI",
        "creatinine-cystatin C equation"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BSA-NORMALIZED eGFR in mL/min/1.73 m^2 from the 2021 CKD-EPI combined",
        "creatinine-cystatin C equation (Methods, Population PK Modeling,",
        "references 17-19). Enters the Table 2 formula as the power ratio",
        "(CRCL / 104.368)^1.12 on CL. Three renal descriptors were screened",
        "head to head and all three improved the fit, but the combined",
        "creatinine-cystatin C equation won: Cockcroft-Gault CLcr, 2012 CKD-EPI",
        "cystatin C eGFR and 2021 CKD-EPI creatinine-cystatin C eGFR gave OFV",
        "reductions of 6.48, 10.50 and 12.13 respectively (Results, Population",
        "PK Modeling). The paper's rationale is specific to this cohort and is",
        "worth carrying forward: serum creatinine is a product of skeletal",
        "muscle metabolism, so in low-body-weight and sarcopenic adults",
        "creatinine-based equations OVERESTIMATE renal function, whereas",
        "cystatin C is produced at a constant rate by all nucleated cells and",
        "is largely independent of muscle mass (Discussion). Cohort median",
        "(IQR) 109.5 (36) mL/min/1.73 m^2 by this equation, against 90.5",
        "(35.25) for the cystatin-C-only equation and 80 (48.4) mL/min for",
        "Cockcroft-Gault CLcr (Table 1). Note that the reference constant",
        "104.368 in the Table 2 formula does NOT equal the Table 1 median of",
        "109.5 even though the text describes it as 'normalized to the median",
        "eGFR value' -- see the note on e_crcl_cl and the vignette Errata.",
        "Patients who developed acute kidney injury or required renal",
        "replacement therapy were excluded, so the model carries no information",
        "about unstable renal function or dialysis. The paper's own dosing",
        "simulations exercised the range 20-130 mL/min/1.73 m^2."
      ),
      source_name        = "eGFRCre-Cys"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened and NOT retained -- the headline negative finding of the",
        "paper. Body weight was tested both as a free covariate and as",
        "allometric scaling with exponents FIXED at 0.75 on CL and 1.0 on Vd,",
        "and neither met the retention criteria (Methods, Population PK",
        "Modeling; Abstract: 'Body weight did not significantly improve model",
        "performance'). Cohort median (IQR) 42 (9.2) kg. The Discussion",
        "attributes the null result to the restricted weight range of a cohort",
        "selected entirely for low body weight rather than to an absence of a",
        "true weight effect, and explicitly calls for studies spanning the full",
        "weight spectrum. A downstream user must NOT read this model as",
        "evidence that piperacillin disposition is weight-independent."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained (Methods, Population PK Modeling). BMI is",
        "also the cohort's inclusion criterion (<= 18.5 kg/m^2), so its range",
        "is restricted by design; 6 of 29 participants (20.7%) had BMI <= 15",
        "(Table 1). No coefficient is reported anywhere in the paper."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained (Methods, Population PK Modeling). Cohort",
        "median (IQR) 64 (14) years; 14 of 29 (48.3%) were older than 65",
        "(Table 1). No coefficient is reported anywhere in the paper."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and not retained (Methods, Population PK Modeling). Table 1",
        "reports 15 of 29 male (51.72%), i.e. 14 female (48.28%). No",
        "coefficient is reported anywhere in the paper."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "piperacillin",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29,
    n_studies      = 1,
    n_observations = 55,
    age_median_iqr = "64 (14) years",
    height_median_iqr = "160 (17) cm",
    weight_median_iqr = "42 (9.2) kg",
    bmi_criterion  = "BMI <= 18.5 kg/m^2 by inclusion criterion; 6/29 (20.7%) had BMI <= 15",
    sex_female_pct = 48.28,
    disease_state  = paste(
      "Hospitalized adults with low body weight receiving",
      "piperacillin-tazobactam; 21 of 29 (72.4%) were admitted to intensive",
      "care. Median (IQR) Charlson Comorbidity Index 5 (4). Patients who",
      "developed acute kidney injury during antibiotic therapy or required",
      "renal replacement therapy were excluded."
    ),
    renal_function = paste(
      "Median (IQR) 2021 CKD-EPI creatinine-cystatin C eGFR 109.5 (36)",
      "mL/min/1.73 m^2; 2012 CKD-EPI cystatin C eGFR 90.5 (35.25)",
      "mL/min/1.73 m^2; Cockcroft-Gault CLcr 80 (48.4) mL/min; serum creatinine",
      "0.5 (0.2) mg/dL; cystatin C 0.84 (0.3) mg/L (available for 24 of 29",
      "participants). Renal function in this cohort skews preserved to",
      "augmented, which is the point of the paper: creatinine-based equations",
      "overestimate it when muscle mass is low."
    ),
    dose_range     = paste(
      "Observed data came from routine therapeutic drug monitoring under the",
      "institutional standard of piperacillin-tazobactam 4500 mg every 6 hours",
      "as a 1-hour infusion for patients with CLcr >= 40 mL/min, with",
      "renal-function-guided modification and non-protocolized low-body-weight",
      "adjustment at the treating physician's discretion. The model was then",
      "applied to simulated regimens of 2250, 3375 and 4500 mg every 6, 8 or 12",
      "hours given as 1-hour or 4-hour infusions."
    ),
    regions        = "Taiwan (single center: National Taiwan University Hospital, Taipei)",
    notes          = paste(
      "Prospective observational therapeutic-drug-monitoring study conducted",
      "January 2020 - December 2022; IRB 201907124RINC. Adults aged 20 years or",
      "older with BMI <= 18.5 kg/m^2. Baseline demographics per Table 1",
      "(continuous variables as median (interquartile range)). Two samples were",
      "drawn at random times within a dosing interval after at least four",
      "consecutive doses to ensure steady state; only one sample was obtained",
      "from three patients, giving 55 samples from 29 patients. TOTAL plasma",
      "piperacillin was quantified by validated UHPLC-ESI-MS/MS with a TZP-d5",
      "internal standard. Only piperacillin was modelled -- tazobactam was not",
      "measured. Estimation was performed in Monolix 2024R1, so the reported",
      "omega is the standard deviation of the normal random effect on the log",
      "scale."
    )
  )

  ini({
    # Structural parameters -- Tseng 2026 Table 2, "Original Dataset" estimates.
    # Bootstrap medians (n = 500) agree closely and are quoted alongside.
    lvc <- log(11.95)
    label("Central volume of distribution (L)")  # Table 2, tvV = 11.95 L, RSE 36.5%, bootstrap median 12.22 (95% CI 8.75-19.78)

    lcl <- log(4.89)
    label("Clearance at the reference eGFR of 104.368 mL/min/1.73 m^2 (L/h)")  # Table 2, tvCL = 4.89 L/hr, RSE 19.0%, bootstrap median 4.98 (95% CI 3.99-6.63)

    # Power exponent of renal function on CL. The reference constant 104.368 is
    # taken verbatim from the Table 2 formula footnote
    #   "CL = tvCL *((eGFR/104.368)^theta_eGFR)*exp^eta_CL"
    # even though the surrounding prose calls it "the median eGFR value" while
    # Table 1 reports a median of 109.5 mL/min/1.73 m^2. The printed equation is
    # authoritative per the standing text-vs-equation rule, and it is also the
    # value that reproduces the paper's own Figure 3 target-attainment table
    # (see the vignette). Not fixed: Table 2 reports an RSE of 29.5% for it.
    e_crcl_cl <- 1.12
    label("Power exponent of CKD-EPI Cr-CysC eGFR (CRCL / 104.368) on CL (unitless)")  # Table 2, theta_eGFR = 1.12, RSE 29.5%, bootstrap median 1.15 (95% CI 0.37-1.76)

    # Fraction of piperacillin unbound in plasma. Fixed, not estimated: "the
    # total concentrations were multiplied by the reported unbound fraction of
    # 0.7 to determine free drug concentrations" (Methods, Drug Administration
    # and Sampling). It is a literature value, not a measurement made in this
    # study -- the Discussion lists the fixed unbound fraction as a limitation
    # because free piperacillin was never assayed. Carried in the model because
    # it is load-bearing for the paper's 100% fT > MIC endpoint.
    fu <- fixed(0.7)
    label("Fraction of piperacillin unbound in plasma (unitless)")  # Methods, Drug Administration and Sampling

    # IIV on CL only. Monolix reports omega as the standard deviation of the
    # normal random effect on the log scale, so the variance below is 0.35^2.
    # The Table 2 CV% column confirms the convention exactly: an omega of
    # 0.3507 (of which 0.35 is the printed rounding) gives
    # sqrt(exp(omega^2) - 1) = 36.17%, the tabulated CV. Shrinkage 0.027.
    etalcl ~ 0.1225  # Table 2, omega_CL = 0.35 (RSE 21.0%, CV 36.17%, shrinkage 0.027)

    # No IIV on Vc: "omega_V was removed from the final model due to a high
    # shrinkage value of 0.832" (Table 2 footnote; Results, Population PK
    # Modeling). This is a deliberate structural choice in the source, not a
    # gap in the transcription -- do not add an eta on lvc.

    # Residual error. Additive, proportional and combined structures were all
    # evaluated; proportional gave the best fit (Methods, Population PK
    # Modeling; Results, Population PK Modeling).
    propSd <- 0.28
    label("Proportional residual error (fraction)")  # Table 2, sigma_prop = 0.28, RSE 13.1%, bootstrap median 0.28 (95% CI 0.22-0.35)
  })

  model({
    # 1. Individual parameters, Tseng 2026 Table 2 formula footnote:
    #      V  = tvV
    #      CL = tvCL * ((eGFR / 104.368)^theta_eGFR) * exp(eta_CL)
    #    Volume carries neither a covariate nor a random effect.
    cl <- exp(lcl + etalcl) * (CRCL / 104.368)^e_crcl_cl
    vc <- exp(lvc)

    # 2. One-compartment disposition with first-order elimination. One- and
    #    two-compartment models were compared and the one-compartment model was
    #    retained; the Discussion notes that the sparse two-samples-per-interval
    #    TDM design limited the ability to discriminate between them.
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # 3. Observation. The assay measured TOTAL plasma piperacillin and the model
    #    was fit on that scale, so Cc is the total concentration. Cfree is the
    #    unbound concentration the paper's PK/PD target is defined on: target
    #    attainment is "100% fT > MIC", evaluated as the steady-state trough of
    #    Cfree against the MIC. The neurotoxicity threshold of 361 mg/L is, by
    #    contrast, a TOTAL trough concentration -- compare it against Cc, not
    #    Cfree.
    Cc <- central / vc
    Cfree <- fu * Cc

    Cc ~ prop(propSd)
  })
}
