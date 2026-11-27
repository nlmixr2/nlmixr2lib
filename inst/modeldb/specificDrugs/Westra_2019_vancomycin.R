Westra_2019_vancomycin <- function() {
  description <- "Two-compartment population pharmacokinetic model for intravenous vancomycin in adults with end-stage kidney disease on intermittent online hemodiafiltration (HDF), developed by iterative two-stage Bayesian analysis in MW\\Pharm 3.83 and evaluated in low-flux hemodialysis patients. Total clearance is a non-renal arm plus a renal arm proportional to residual creatinine clearance plus, during a dialysis session only, an extracorporeal arm equal to blood flow times a dialyser extraction ratio. Both volumes scale with a fat-corrected lean body mass."
  reference <- paste(
    "Westra N, Proost JH, Franssen CFM, Wilms EB, van Buren M, Touw DJ.",
    "Vancomycin pharmacokinetic model development in patients on intermittent",
    "online hemodiafiltration. PLoS One. 2019;14(5):e0216801.",
    "doi:10.1371/journal.pone.0216801. PMCID: PMC6516654.",
    sep = " "
  )
  vignette <- "Westra_2019_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The quadratic coefficient of the assay-error polynomial has no canonical
  # residual-error name; see the residual-error block in ini().
  paper_specific_residual_sds <- c("quadSd")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the fat-corrected lean body mass",
        "LBMc = LBM + ffat * (BW - LBM) with ffat = 0.4 (Methods, Eq 3),",
        "which multiplies both per-kg volumes. Clearance does not scale with",
        "weight. Development cohort median 72.9 kg, range 53-95.7 kg (Table 1)."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Input to the sex-specific lean body mass printed in Methods below",
        "Eq 3: LBM = 50.0 + 0.9 * (Height - 152) for men and",
        "45.5 + 0.9 * (Height - 152) for women. The paper prints no cap of LBM",
        "at total body weight, so none is applied; for a short-and-heavy or",
        "tall-and-light subject LBM can exceed WT and then LBMc > WT.",
        "Development cohort median 165 cm, range 145-189 cm (Table 1)."
      ),
      source_name = "Height"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the intercept of the lean-body-mass formula only (50.0 men,",
        "45.5 women). Development cohort 9 men / 8 women (Table 1)."
      ),
      source_name = "Gender"
    ),
    CRCL = list(
      description = paste(
        "Residual renal function, measured as the mean of the 24-h urinary",
        "creatinine and urea clearances (U*V/P), in ABSOLUTE mL/min (not",
        "BSA-normalised)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods Eq 1: renal function (mL/min) = 0.5 * (creatinine clearance",
        "+ urea clearance), each computed from a 24-h urine collection and",
        "the plasma level. Table 1 'CLcr (ml/min)' (median 2.9, range",
        "0.4-16.5 mL/min, mean 4.4) matches the supplement S1 'Residual renal",
        "clearance ml/min' column, not its '/1.73 m2' column, so the",
        "un-normalised value is the one the model used. Enters Eq 2 as",
        "fr * CLcr with fr = 0.1 fixed; converted to L/h inside model() by",
        "0.06. The dosing tables assume CLcr = 4.4 mL/min."
      ),
      source_name = "CLcr"
    ),
    BFR = list(
      description = "Extracorporeal blood flow through the dialyser during a session",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Qec of Methods Eq 2, CL = CLm + fr * CLcr + Qec * Eec. Multiplied by",
        "RRT_HEMODIAL_ACTIVE, so its value between sessions is irrelevant",
        "('With the exception of dialysis sessions the Qec is zero'). Converted",
        "to L/h by 0.06 inside model(). Per-patient blood flows are not",
        "reported; the dosing tables (Tables 3-4) use 200, 250 and 300 mL/min,",
        "and the Discussion's 74 L/day dialysis clearance at Eec = 0.212",
        "implies about 242 mL/min."
      ),
      source_name = "Qec"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemo(dia)filtration session running indicator, 1 during a session, 0 otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no session running)",
      notes = paste(
        "Time-varying within subject; gates the extracorporeal arm BFR * Eec.",
        "The development cohort was dialysed with online HDF on FX80",
        "high-flux filters and the evaluation cohort with low-flux HD; the",
        "paper reports that one model describes both. All sessions lasted",
        "4 h (Discussion) and vancomycin was infused during the last 30 or",
        "60 min of a session (Methods, Patient data)."
      ),
      source_name = "Qec > 0"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 17,
    n_studies = 1,
    n_observations = 159,
    age_range = "46-85 years",
    age_median = "63 years",
    weight_range = "53-95.7 kg",
    weight_median = "72.9 kg",
    height_range = "145-189 cm",
    height_median = "165 cm",
    sex_female_pct = 47,
    race_ethnicity = "Not reported.",
    disease_state = "End-stage kidney disease on intermittent online hemodiafiltration, receiving intravenous vancomycin",
    renal_function = "Residual renal function (24-h urine creatinine/urea clearance) median 2.9 mL/min, range 0.4-16.5 mL/min",
    dose_range = "Intravenous vancomycin infused during the last 30 or 60 min of a 4-h dialysis session; doses not tabulated",
    regions = "Netherlands (HagaZiekenhuis, The Hague)",
    notes = paste(
      "Development cohort: 17 online-HDF patients (FX80 high-flux filters)",
      "treated January 2002 - January 2007 with more than one extracorporeal session during the",
      "vancomycin course, 159 samples (median 8 per patient, range 2-24),",
      "drawn at the start of dialysis, before the vancomycin infusion and",
      "after the end of dialysis (Methods; Table 1). 135 of 152 screened",
      "patients were excluded for missing data. Evaluation cohort (not used",
      "for estimation): 21 low-flux HD patients at UMCG Groningen, 2009-2016,",
      "132 samples, median weight 70 kg (38-133), residual CLcr median",
      "2.8 mL/min (Table 1)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Table 2 'Population parameters for vancomycin during online HDF
    # dialysis', Final model. MW\Pharm KinPop reports each parameter as a
    # 'Mean' and an inter-individual 'SD' on the natural scale. IIV was
    # 'assumed to be log-normally distributed' (Methods), so the Mean is
    # carried as the typical value and the SD is converted to a log-scale
    # variance with omega^2 = log(1 + (SD / Mean)^2). fr and CL12 were
    # 'fixed on respectively 0.1 and 9.96 L/h' (Results) with SD 0.
    # ---------------------------------------------------------------------
    lcl_nonren <- log(0.473)
    label("Non-renal clearance CLm (L/h)") # Table 2 CLm Mean = 0.473 L/h (95% CI 0.291-0.633); Abstract CLm = 0.473 (0.271) L/h

    e_crcl_cl_renal <- fixed(0.1)
    label("Renal vancomycin clearance as a fraction fr of creatinine clearance (unitless)") # Table 2 Fr Mean = 0.100 (fixed), SD 0; Methods Eq 2

    lvc <- log(0.278)
    label("Central volume V1 per kg fat-corrected lean body mass (L/kg LBMc)") # Table 2 V1 Mean = 0.278 L/kgLBMc (95% CI 0.232-0.349)

    lq <- fixed(log(9.96))
    label("Intercompartmental clearance CL12 (L/h)") # Table 2 CL12 Mean = 9.960 (fixed) L/h, SD 0; Results 'fixed ... at a literature value'

    lvp <- log(0.686)
    label("Peripheral volume V2 per kg fat-corrected lean body mass (L/kg LBMc)") # Table 2 V2 Mean = 0.686 L/kgLBMc (95% CI 0.482-1.087)

    ledia <- log(0.212)
    label("Dialyser extraction ratio Eec applied to extracorporeal blood flow (fraction)") # Table 2 Eec Mean = 0.212 (95% CI 0.126-0.281)

    ffat_v <- fixed(0.4)
    label("Fat-distribution factor ffat in LBMc = LBM + ffat * (WT - LBM) (unitless)") # Methods below Eq 3: 'ffat distribution over fat factor and was set to 0.4'

    # IIV: omega^2 = log(1 + CV^2), CV = SD / Mean from Table 2.
    etalcl_nonren ~ 0.283869 # Table 2 CLm SD 0.271 / Mean 0.473 = CV 57.29 pct; eta-shrinkage 9.1 pct
    etalvc ~ 0.103926 # Table 2 V1 SD 0.092 / Mean 0.278 = CV 33.09 pct; eta-shrinkage 25.3 pct
    etalvp ~ 0.213880 # Table 2 V2 SD 0.335 / Mean 0.686 = CV 48.83 pct; eta-shrinkage 30.4 pct
    etaledia ~ 0.100688 # Table 2 Eec SD 0.069 / Mean 0.212 = CV 32.55 pct; eta-shrinkage 63.7 pct

    # Residual error. MW\Pharm weighted the fit by the assay error:
    # 'SD = 1.3842 + 0.0626 x C + 0.0018 x C^2, where C is the vancomycin
    # concentration in mg/L' (Methods, Vancomycin assay). Supplement S2's
    # 'SD' column reproduces this polynomial at the measured concentration
    # (e.g. C = 11 mg/L -> 2.2906). The three terms are added linearly, not
    # in quadrature, and are assay characteristics rather than fitted
    # values, so all are fixed. No additional residual term is reported.
    addSd <- fixed(1.3842)
    label("Assay-error SD intercept (mg/L)") # Methods 'SD = 1.3842 + 0.0626 x C + 0.0018 x C2'
    propSd <- fixed(0.0626)
    label("Assay-error SD linear coefficient (fraction of C)") # Methods 'SD = 1.3842 + 0.0626 x C + 0.0018 x C2'
    quadSd <- fixed(0.0018)
    label("Assay-error SD quadratic coefficient (L/mg)") # Methods 'SD = 1.3842 + 0.0626 x C + 0.0018 x C2'
  })

  model({
    # Lean body mass (kg), Methods below Eq 3: 50.0 + 0.9 * (Height - 152)
    # for men and 45.5 + 0.9 * (Height - 152) for women.
    lbm <- (50.0 - 4.5 * SEXF) + 0.9 * (HT - 152)
    # Fat-corrected lean body mass (kg), Methods Eq 3.
    lbmc <- lbm + ffat_v * (WT - lbm)

    # Methods Eq 2: CL = CLm + fr * CLcr + Qec * Eec, with Qec = 0 outside a
    # dialysis session. CLcr and Qec are in mL/min; 0.06 converts to L/h.
    # The gated total is written into `cl` itself.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- e_crcl_cl_renal * CRCL * 0.06
    edia <- exp(ledia + etaledia)
    cl_hemodialysis <- BFR * 0.06 * edia
    cl <- cl_nonren + cl_renal + RRT_HEMODIAL_ACTIVE * cl_hemodialysis

    vc <- exp(lvc + etalvc) * lbmc
    vp <- exp(lvp + etalvp) * lbmc
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    sdCc <- addSd + propSd * Cc + quadSd * Cc^2
    Cc ~ add(sdCc)
  })
}
