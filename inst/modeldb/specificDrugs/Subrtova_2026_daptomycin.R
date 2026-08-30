Subrtova_2026_daptomycin <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for",
    "intravenous daptomycin in adult patients with serious Gram-positive",
    "infections (infective endocarditis, bone and joint infection, sepsis,",
    "catheter-associated infection), developed from routine therapeutic drug",
    "monitoring data at a single Czech centre (143 serum concentrations from",
    "31 patients, 2022-2025; patients on renal replacement therapy excluded).",
    "Clearance is scaled by a power function of the time-varying CKD-EPI 2021",
    "estimated glomerular filtration rate normalised to 90 mL/min:",
    "CL = 0.69 * (eGFR/90)^0.40 L/h. eGFR was the only retained covariate --",
    "age, height, total / ideal / adjusted body weight, lean body mass, body",
    "surface area, sex, serum creatinine, Cockcroft-Gault creatinine",
    "clearance and serum urea were all tested and rejected, and the authors",
    "argue the absence of a body-size effect is genuine rather than a",
    "sample-size artefact. Inter-individual variability is log-normal on the",
    "volume of distribution, on clearance, and also on the eGFR power",
    "exponent itself: eGFR entered the Monolix model as a time-varying",
    "regressor rather than as a built-in covariate, which makes the exponent",
    "an ordinary individual parameter carrying its own random effect.",
    "Residual variability is proportional. The paper converts the",
    "clearance-eGFR relationship into a once-daily maintenance-dose nomogram",
    "targeting AUC24 = 832 mg*h/L at MIC 1 mg/L and 665.5 mg*h/L at MIC",
    "0.5 mg/L (the midpoints of the 666-998 and 333-998 mg*h/L therapeutic",
    "windows), with an initial loading dose of 1.36 times the maintenance",
    "dose to offset the day-1 accumulation shortfall."
  )
  reference <- paste(
    "Subrtova P, Rozsivalova P, Halvova P, Malakova J, Paterova P,",
    "Ryskova L, Maly J, Michalek P, Slanar O, Sima M. Population",
    "pharmacokinetic model and dosing nomogram for daptomycin in adult",
    "patients with serious Gram-positive infections: emphasizing the role of",
    "loading doses and renal function-based adjustment. Antimicrob Agents",
    "Chemother. 2026;70(3):e01532-25. doi:10.1128/aac.01532-25."
  )
  vignette <- "Subrtova_2026_daptomycin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate calculated with the CKD-EPI",
        "2021 (race-free) creatinine equation. Time-varying: entered the",
        "Monolix model as a regressor, so each concentration record carries",
        "the eGFR contemporaneous with that sample rather than a single",
        "baseline value."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference value 90 mL/min, described in the Results as 'normal",
        "renal function status'; the model applies (CRCL/90)^e_crcl_cl to",
        "clearance (Subrtova 2026 Table 2, Fixed effects block).",
        "SIZE NORMALISATION: the paper labels this column 'mL/min' in both",
        "Table 1 and the Table 2 footnote, but the CKD-EPI 2021 equation",
        "natively returns mL/min/1.73 m^2 and the paper never states that it",
        "de-indexed for body surface area. The cohort values reported in",
        "Table 1 (median 81.6, IQR 58.9-98.9, range 16.3-130.1) are the ones",
        "the model was fitted to and are what should be supplied, whichever",
        "normalisation they actually carry. Supplying a BSA-de-indexed value",
        "in a cohort whose median BSA is 1.98 m^2 would inflate eGFR by",
        "about 14% and, through the 0.40 exponent, clearance by about 5%.",
        "Two other renal-function measures were tested and rejected in",
        "favour of this one: serum creatinine and Cockcroft-Gault creatinine",
        "clearance (Table 1 median 54.2 mL/min), as was the CKD-EPI 2012",
        "eGFR (Table 1 median 73.5 mL/min). Do not substitute any of those",
        "for this column -- Cockcroft-Gault in particular runs roughly 30%",
        "lower in this cohort and would bias clearance downward."
      ),
      source_name        = "eGFR CKD-EPI 2021"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested as a continuous covariate on both Vd and CL and rejected;",
        "the Discussion argues explicitly that the cohort's 48-124 kg spread",
        "was wide enough to expose a clinically meaningful body-size effect",
        "had one existed. Retained here as documentation of the covariate",
        "screen only."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Tested (range 158-191 cm) and rejected (Subrtova 2026 Results, Population pharmacokinetic analysis)."
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested (median 68 kg) and rejected."
    ),
    ABW = list(
      description = "Adjusted body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested (median 73 kg) and rejected. No canonical entry exists in",
        "inst/references/covariate-columns.md for adjusted body weight; the",
        "key is used here as documentation of the covariate screen only and",
        "is deliberately NOT proposed as a new canonical, since the model",
        "never references it."
      )
    ),
    LBM = list(
      description = "Lean body mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested (median 59 kg) and rejected."
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Tested (median 1.98 m^2) and rejected."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested (median 63 years, range 29-93) and rejected."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "The only categorical covariate tested (7 of 31 patients female) and rejected."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Tested as a time-varying regressor (median 90 umol/L, range 36-324)",
        "and rejected in favour of CKD-EPI 2021 eGFR."
      )
    ),
    BUN = list(
      description = "Serum urea (the canonical BUN column accepts mmol/L urea; 1 mmol/L urea ~= 2.80 mg/dL BUN)",
      units       = "mmol/L",
      type        = "continuous",
      notes       = paste(
        "Tested as a time-varying regressor (Table 1 median 7.1 mmol/L, IQR",
        "4.6-13.7, range 2.2-29.7) and rejected in favour of CKD-EPI 2021",
        "eGFR. Recorded under the canonical name BUN rather than the paper's",
        "own label 'Urea'."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "daptomycin",
      units    = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 31L,
    n_studies        = 1L,
    n_observations   = 143L,
    age_range        = "29-93 years",
    age_median       = "63 years",
    weight_range     = "48-124 kg",
    weight_median    = "84 kg",
    height_range     = "158-191 cm",
    sex_female_pct   = 22.6,
    disease_state    = paste(
      "Serious Gram-positive bacterial infection: infective endocarditis",
      "(n = 14, 45.2%), bone and joint infection (n = 6, 19.4%), sepsis",
      "(n = 3, 9.7%), catheter-associated infection (n = 3, 9.7%).",
      "Causative organisms Staphylococcus spp. (n = 21, 67.7%),",
      "Enterococcus spp. (n = 4, 12.9%), Corynebacterium spp. (n = 2, 6.5%)."
    ),
    renal_function   = paste(
      "eGFR CKD-EPI 2021 median 81.6 mL/min (IQR 58.9-98.9, range",
      "16.3-130.1); eGFR CKD-EPI 2012 median 73.5 mL/min (range 14.9-132.1);",
      "Cockcroft-Gault creatinine clearance median 54.2 mL/min (range",
      "15.2-176.0). Patients receiving renal replacement support were",
      "EXCLUDED, so the model carries no information about dialysis or CRRT."
    ),
    dose_range       = paste(
      "350-1,000 mg (median 750 mg), i.e. 5-13 mg/kg (median 10 mg/kg),",
      "every 24 or 48 h as a 30- or 60-min intravenous infusion; the initial",
      "regimen was set by the attending physician and subsequently adjusted",
      "by therapeutic drug monitoring."
    ),
    regions          = "Czech Republic (single centre: University Hospital Hradec Kralove)",
    notes            = paste(
      "Retrospective, single-centre, cross-sectional TDM study, May 2022 to",
      "July 2025. Baseline demographics in Subrtova 2026 Table 1. Sampling",
      "was routine TDM: 80 trough samples (1 h before the next dose) and 63",
      "peak samples (2 h after the start of the infusion), 1-19 records per",
      "patient (mean 4.6). Assay UPLC-MS/MS, LLOQ 0.32 mg/L, calibrated to",
      "127.7 mg/L. Estimation in Monolix 2024R1 (SAEM); simulations in",
      "Simulx 2024R1. Daptomycin MIC was measured in 19 of 31 cases (median",
      "and mode 0.5 mg/L, range 0.016-1.5). Note that only 31 patients",
      "informed three fixed effects and three random effects, and that all",
      "concentrations are TOTAL (not unbound) daptomycin -- the paper states",
      "unbound concentrations were not measured, consistent with the PK/PD",
      "target being defined on total drug."
    )
  )

  ini({
    # ---- Structural parameters (Subrtova 2026 Table 2, Fixed effects) ----
    # Vd is a single population value; no covariate was retained on it.
    lvc <- log(11.47)
    label("Volume of distribution (L)")                                            # Table 2: Vd_pop = 11.47 L (R.S.E. 6.88%; bootstrap 11.56, 11.45-11.65)

    # CL_pop is the clearance of a patient at the reference eGFR of 90 mL/min,
    # not a covariate-free intercept. Results: "for an individual with normal
    # renal function status (eGFR = 90 mL/min), daptomycin Vd and CL were
    # 11.47 L and 0.69 L/h ... which corresponds to a t1/2 of 11.5 h".
    lcl <- log(0.69)
    label("Clearance at the reference eGFR of 90 mL/min (L/h)")                     # Table 2: CL_pop = 0.69 L/h (R.S.E. 5.13%; bootstrap 0.70, 0.69-0.70)

    # Power exponent on (eGFR / 90) for CL. Log-transformed here (a) so the
    # exponent stays positive and (b) so its log-normal IIV is mu-referenced;
    # the untransformed canonical covariate-effect name e_crcl_cl is recovered
    # in model(). Precedent for the `le_` form: Kong_2025_piperacillin_tazobactam.R.
    le_crcl_cl <- log(0.40)
    label("Log power exponent on (CRCL/90) for CL (unitless)")                      # Table 2: beta_CL_eGFR = 0.40 (R.S.E. 32.40%; bootstrap 0.41, 0.39-0.43)

    # ---- Inter-individual variability (Subrtova 2026 Table 2, "Between
    # subject variability (%)") ----
    # SCALE. The column is headed "(%)" and reports 18.0 / 20.0 / 39.0. It is
    # read here as CV%, so omega^2 = log(CV^2 + 1). The column is definitely
    # NOT a variance: with N = 31 subjects the Cramer-Rao floor on the RSE of
    # a variance is sqrt(2/31) = 25.4%, and the CL row's reported R.S.E. of
    # 22.3% sits below it, which is arithmetically impossible for a variance
    # (the SD floor is sqrt(1/62) = 12.7%). Separating CV% from a bare SD is
    # not possible from the paper, but the two readings differ by at most 4%
    # relative on omega at these magnitudes (e.g. 39% CV -> omega = 0.376 vs
    # SD 0.39); see the vignette Errata.
    etalvc ~ 0.0319
    # Table 2: BSV Vd = 18.0% (R.S.E. 29.8%); omega^2 = log(1 + 0.18^2) = 0.03189

    etalcl ~ 0.0392
    # Table 2: BSV CL = 20.0% (R.S.E. 22.3%); omega^2 = log(1 + 0.20^2) = 0.03922

    # IIV on the eGFR exponent. Unusual for a popPK covariate coefficient, but
    # it follows directly from the estimation setup: Methods state the
    # time-varying covariates "were tested by incorporating them as
    # regressors", and a Monolix regressor forces its coefficient to be
    # declared as an individual parameter, which by default carries a random
    # effect. Methods also state "the model parameters were assumed to be
    # log-normally distributed", which fixes the distribution.
    #
    # CONFIRMED BY PARAMETER COUNT, not just by inference. Supplement Table S1
    # reports OFV, AIC and BIC for every candidate model, and Monolix computes
    # AIC = OFV + 2p and BIC = OFV + p*log(N_subjects) with N = 31. Solving for
    # p on the selected final model row (OFV 1173.24, AIC 1187.24, BIC 1197.28)
    # gives p = (1187.24 - 1173.24)/2 = 7 and (1197.28 - 1173.24)/log(31) =
    # 7.00. Seven estimated parameters is Vd_pop + CL_pop + beta_CL_eGFR (3
    # fixed effects) + THREE random effects + 1 proportional error. Six would
    # be the count if Table 2's third "Between subject variability" row were a
    # mislabelled duplicate rather than a real omega. The same arithmetic
    # reproduces every other Table S1 row exactly (base 1-cmt p = 5: AIC
    # 1211.98, BIC 1219.15; 2-cmt and split renal/non-renal CL both p = 9),
    # which confirms the AIC/BIC convention being solved against.
    etale_crcl_cl ~ 0.1416
    # Table 2: BSV beta_CL_eGFR = 39.0% (R.S.E. 29.8%); omega^2 = log(1 + 0.39^2) = 0.14157

    # ---- Residual variability ----
    # Monolix proportional error model y = f * (1 + b * eps), b = 0.28, which
    # is exactly nlmixr2's prop(propSd). A combined and a constant error model
    # were both tested and rejected (Methods; Results, first paragraph of
    # "Population pharmacokinetic analysis").
    propSd <- 0.28
    label("Proportional residual error (fraction)")                                 # Table 2: Error model parameter, Proportional = 0.28 (R.S.E. 7.72%)
  })

  model({
    # 1. Covariate effect. Recover the canonical untransformed covariate-effect
    #    parameter from its log-scale, mu-referenced counterpart.
    e_crcl_cl <- exp(le_crcl_cl + etale_crcl_cl)

    # 2. Individual parameters. Subrtova 2026 Table 2, Fixed effects block:
    #      Vd = Vd_pop
    #      CL = CL_pop * (eGFR / 90)^beta_CL_eGFR
    #    CRCL is time-varying, so cl is re-evaluated at every record.
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (CRCL / 90)^e_crcl_cl

    # 3. Micro-constant. Written out rather than using linCmt() because a
    #    time-varying regressor on CL makes kel time-varying, which the
    #    linCmt() solution cannot represent.
    kel <- cl / vc

    # 4. ODE system. Intravenous administration only (30- or 60-min infusion);
    #    there is no absorption compartment and no bioavailability term.
    d/dt(central) <- -kel * central

    # 5. Observation and error. Total (not unbound) daptomycin serum
    #    concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
