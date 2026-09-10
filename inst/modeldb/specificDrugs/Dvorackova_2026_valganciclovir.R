Dvorackova_2026_valganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for ganciclovir in adult lung",
    "transplant recipients receiving intravenous ganciclovir and/or oral",
    "valganciclovir for cytomegalovirus prophylaxis or treatment",
    "(Dvorackova 2026). Both routes are described by one joint fit of 379",
    "serum ganciclovir concentrations from 110 patients: intravenous",
    "ganciclovir doses enter the central compartment directly, and oral",
    "valganciclovir doses enter a first-order depot (ka 0.334 1/h) with a",
    "0.563 h lag time and bioavailability 0.575, held on the logit scale",
    "(theta_F 0.304) so it cannot leave (0, 1). Clearance is a LINEAR (not",
    "power) function of CKD-EPI 2021 estimated glomerular filtration rate,",
    "CL = 2.05 + 4.96 * (eGFR / 85.2) L/h, giving 7.01 L/h at the cohort",
    "median eGFR of 85.2 mL/min/1.73 m^2 and rising by 0.058 L/h per",
    "mL/min/1.73 m^2; eGFR was the only covariate retained, and body weight",
    "was NOT a covariate on either volume. Central volume is 43.1 L,",
    "peripheral volume 140 L and intercompartmental clearance 2.1 L/h.",
    "Interindividual variability is log-normal on CL (variance 0.165) and on",
    "central volume (variance 0.431) and normal on logit-F (variance 2.63);",
    "the data did not support variability on ka, peripheral volume or",
    "intercompartmental clearance, and inter-occasion variability did not",
    "improve the fit. Residual variability is proportional (variance 0.236,",
    "i.e. SD 0.486). The oral dose record is in administered valganciclovir",
    "mg: the model applies the paper's molar prodrug conversion",
    "(255.23 / 354.362 = 0.720 g ganciclovir per g valganciclovir) inside",
    "f(depot), so `depot` and `central` both hold ganciclovir mg."
  )
  reference <- paste(
    "Dvorackova E, Michalickova D, Petrus J, Klapkova E, Dutkova A,",
    "Kotowski T, Krekels EHJ, Havlin J, Lischke R, Slanar O (2026).",
    "Population pharmacokinetics and dose optimization of valganciclovir",
    "and ganciclovir in lung transplant recipients.",
    "Med Princ Pract 35:169-180. doi:10.1159/000548942",
    sep = " "
  )
  vignette <- "Dvorackova_2026_valganciclovir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Both routes the source study used are dosable, and they are NOT
  # interchangeable in dose units:
  #   * `central` takes an intravenous GANCICLOVIR dose in mg (the study
  #     protocol gave 5 mg/kg every 12 h as a 60-min infusion, so use a
  #     rate / duration on the dose record).
  #   * `depot` takes an oral VALGANCICLOVIR dose in mg as prescribed
  #     (e.g. 900 mg). The molar prodrug conversion to ganciclovir
  #     equivalents is applied inside f(depot), so the dose record must NOT
  #     be pre-converted.
  dosing <- c("depot", "central")

  # Issue #482. `depot` holds ganciclovir mg, not valganciclovir mg: f(depot)
  # applies both the bioavailability and the molar prodrug conversion at dose
  # time, so the amount actually entering the state is
  # amt_valganciclovir * F * (255.23 / 354.362). The assayed matrix is serum
  # -- Dvorackova 2026 Methods 'Study Design': "Serum blood concentrations of
  # GCV for PK analysis were obtained just before administration of the next
  # dose (Ctrough) and at least two other times on the same day".
  compartmentData <- list(
    depot = list(
      analyte = "ganciclovir", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "ganciclovir", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ganciclovir", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate, body-surface-area-normalised",
        "to 1.73 m^2, calculated with the race-free CKD-EPI 2021",
        "creatinine-based equation (Inker et al. 2021, the paper's",
        "reference [15]). NOT a measured creatinine clearance and NOT a",
        "cystatin-C-based estimate."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Enters clearance",
        "LINEARLY, not as a power term:",
        "CL = CLp + theta_eGFR * (CRCL / 85.2) with CLp = 2.05 L/h and",
        "theta_eGFR = 4.96 L/h (Dvorackova 2026 Table 3). The 85.2",
        "normalising constant is the cohort median eGFR (Table 2). Because",
        "the relationship is linear with a POSITIVE intercept, clearance",
        "does not go to zero as eGFR does, and the model must not be used",
        "outside the observed eGFR range: Discussion 'the linear decline in",
        "CL with eGFR was estimated based on data with a lowest eGFR of",
        "16.8 mL/min/1.73 m^2 and extrapolations with this model below this",
        "range should be performed with extreme caution'. Observed range",
        "16.8-153 mL/min/1.73 m^2. Time-varying if the source dataset",
        "supplied repeated creatinine measurements; the paper does not",
        "state whether the eGFR column was baseline-only or updated per",
        "sample."
      ),
      source_name        = "eGFR"
    )
  )

  # Covariates Dvorackova 2026 screened by stepwise covariate modelling but
  # did NOT retain in the final model -- Results 'Population PK Model': "No
  # other covariates were found to be statistically significant for any other
  # PK parameter." The paper reports no point estimate for any of these, so
  # they are documentation of the covariate screen only and are deliberately
  # absent from model().
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested as a covariate on CL and on central volume V1 (Methods",
        "'PK Model Development'); not significant. Cohort median 75 kg",
        "(IQR 67-86, Table 2). Discussion: 'Like our study, other studies",
        "describing VGCV/GCV PK in solid organ transplant patients also did",
        "not identify body weight as a covariate for V1 and V2'. There is",
        "therefore NO allometric term in this model -- do not add one."
      ),
      source_name        = "body weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tested as a covariate on CL, V1 and F; not significant. Cohort",
        "70 male / 40 female (64% / 36%), Table 2. The source recorded",
        "'sex' without stating its coding direction; the canonical",
        "female-indicator orientation is used here for documentation only."
      ),
      source_name        = "sex"
    ),
    AGE = list(
      description        = "Age at the time of the PK sampling occasion.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded as a potential covariate (Methods 'Study Design') but not",
        "retained. Cohort median 55 years (IQR 46-62), Table 2. Note that",
        "age also enters the model indirectly through the CKD-EPI 2021",
        "equation used to compute CRCL."
      ),
      source_name        = "age"
    ),
    HT = list(
      description        = "Body height.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded as a potential covariate (Methods 'Study Design') but not",
        "retained and not summarised in Table 2. No point estimate",
        "reported."
      ),
      source_name        = "height"
    ),
    CREAT = list(
      description        = "Serum creatinine.",
      units              = "not stated by the source",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded as a potential covariate (Methods 'Study Design') but not",
        "retained as a covariate in its own right; it enters the model only",
        "through the CKD-EPI 2021 eGFR carried as CRCL. The source does not",
        "state whether creatinine was recorded in umol/L (the usual Czech",
        "laboratory unit) or mg/dL, and does not summarise it in Table 2."
      ),
      source_name        = "serum creatinine"
    ),
    DIS_CF = list(
      description        = "Cystic fibrosis as the indication for transplantation, 1 = yes, 0 = no.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cystic fibrosis)",
      notes              = paste(
        "Tested as a covariate on CL, V1 and F (Methods 'PK Model",
        "Development'); not significant. 9 of 110 patients (8%), Table 2."
      ),
      source_name        = "cystic fibrosis"
    ),
    CONMED_AZOLE = list(
      description        = paste(
        "Co-treatment with an azole antifungal (voriconazole, posaconazole",
        "or fluconazole), 1 = yes, 0 = no."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no azole antifungal co-treatment)",
      notes              = paste(
        "Tested as a covariate on CL (Methods 'PK Model Development'); not",
        "significant. 31 of 110 patients (28%), Table 2. The three azoles",
        "were pooled into a single indicator by the source; it does not",
        "report per-drug effects."
      ),
      source_name        = "co-medication with antimycotics"
    ),
    T_TRANSPLANT = list(
      description        = "Time elapsed since lung transplantation.",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded as a potential covariate (Methods 'Study Design') but not",
        "retained. Median 14.5 days from transplantation to the first drawn",
        "concentration (range 1-1,936), Table 2. This variable also defined",
        "the three inter-occasion-variability strata that were tested and",
        "rejected (up to month 1, 1-2 months, more than 2 months after",
        "transplantation)."
      ),
      source_name        = "time after transplantation"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 110L,
    n_studies      = 1L,
    n_observations = 379L,
    age_median     = "55 years",
    age_range      = "IQR 46-62 years (inclusion required age > 18 years)",
    weight_median  = "75 kg",
    weight_range   = "IQR 67-86 kg",
    sex_female_pct = 36,
    disease_state  = paste(
      "Adult lung transplant recipients receiving ganciclovir or",
      "valganciclovir as prophylaxis against, or treatment of,",
      "cytomegalovirus infection. All received basiliximab or",
      "antithymocyte globulin induction plus triple maintenance",
      "immunosuppression with tacrolimus, mycophenolate mofetil and",
      "prednisone/methylprednisolone. Patients with a second solid-organ",
      "transplant or on renal replacement therapy were excluded. Cystic",
      "fibrosis was the transplant indication in 9 patients (8%)."
    ),
    renal_function = paste(
      "eGFR (CKD-EPI 2021) median 85.2 mL/min/1.73 m^2, range",
      "16.8-153 mL/min/1.73 m^2. No patient was on renal replacement",
      "therapy (exclusion criterion)."
    ),
    co_medication  = "Azole antifungal (voriconazole, posaconazole or fluconazole) in 31 patients (28%).",
    dose_range     = paste(
      "Standard protocol: intravenous ganciclovir 5 mg/kg every 12 h as a",
      "60-min infusion at concentrations not exceeding 10 mg/mL for the",
      "first 14 days, then oral valganciclovir for at least 90 days to",
      "12 months, renally adjusted per Table 1 (900 mg q12h at eGFR",
      "> 90 and 60-90; 900 mg morning plus 450 mg evening at eGFR 30-60;",
      "450 mg q12h at eGFR 15-30 mL/min/1.73 m^2). Intravenous ganciclovir",
      "in Table 1 is 500 / 400 / 200 / 100 mg q12h over the same eGFR",
      "strata."
    ),
    regions        = "Czech Republic (single-centre: Prague Lung Transplant Program, Motol University Hospital and General University Hospital in Prague)",
    notes          = paste(
      "Demographics from Dvorackova 2026 Table 2; values are median",
      "(interquartile range) except eGFR and time since transplantation,",
      "which are median (min-max). Prospective, open-label",
      "(laboratory-blinded) PK study conducted January 2020 - July 2024.",
      "379 ganciclovir concentrations (185 after intravenous ganciclovir,",
      "194 after oral valganciclovir), median 3 samples per patient",
      "(range 1-12), concentrations 0.1-19.2 mg/L. Samples were typically",
      "collected at steady state (regimen unchanged for at least 5 days).",
      "11 patients (10%) contributed concentrations on more than one",
      "occasion. Model built in NONMEM 7.4.0 with FOCE-I."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Dvorackova 2026 Table 3, "Final model
    # (RSE %)" column. The paper estimated these on the natural scale;
    # they are carried log-transformed here per library convention so
    # they cannot go negative.
    # ------------------------------------------------------------------
    lka <- log(0.334)
    label("Absorption rate constant of oral valganciclovir ka (1/h)")
    # Table 3 kA = 0.334 (RSE 27%), bootstrap 0.337 (0.202-0.502).
    # Table 3 prints the units of kA as "L/h", which is a typo: a
    # first-order absorption rate constant has units 1/h. See the vignette
    # Errata.

    ltlag <- log(0.563)
    label("Absorption lag time of oral valganciclovir Tlag (h)")
    # Table 3 Tlag = 0.563 h (RSE 34%), bootstrap 0.615 (0.105-0.881).
    # Results: "Incorporation of lag time decreased OFV by 6.476 points and
    # was therefore retained in the model."

    # Clearance is LINEAR in eGFR, not a power function:
    #   CL = CLp + theta_eGFR * (eGFR / 85.2)   [L/h]
    # `lcl` therefore carries CLp -- the eGFR-independent INTERCEPT of that
    # line -- and NOT the typical clearance. Typical clearance at the cohort
    # median eGFR of 85.2 mL/min/1.73 m^2 is 2.05 + 4.96 = 7.01 L/h, which
    # is what the Discussion compares against the 2.2-15.8 L/h reported by
    # the six other transplant popPK studies for "a typical male individual
    # from our study with an eGFR of 85.2 mL/min/1.73 m^2 and body weight of
    # 75 kg".
    lcl <- log(2.05)
    label("Intercept CLp of the linear eGFR clearance model (L/h)")
    # Table 3 CLp = 2.05 (RSE 31%), bootstrap 2.017 (0.760-3.882).

    lvc <- log(43.1)
    label("Central compartment volume V1 (L)")
    # Table 3 V1 = 43.1 L (RSE 19%), bootstrap 42.5 (10.9-57.8).

    lvp <- log(140)
    label("Peripheral compartment volume V2 (L)")
    # Table 3 V2 = 140 L (RSE 37%), bootstrap 118 (12-960). Results: the
    # bootstrap median for V2 is the one parameter outside the 10% agreement
    # band, "within 15% of the final estimate".

    lq <- log(2.1)
    label("Intercompartmental clearance Q (L/h)")
    # Table 3 Q = 2.1 (RSE 20%), bootstrap 2.11 (0.77-4.17).

    # Bioavailability of oral valganciclovir, relative to the
    # ganciclovir-equivalent dose. Held on the logit scale by the source:
    # Results, "For F, a logit transformation was applied and a normal
    # distribution for IIV was incorporated in the logit domain."
    logitfdepot <- 0.304
    label("Logit of oral valganciclovir bioavailability F (unitless; F = 0.575)")
    # Table 3 theta_F = 0.304 (RSE 48%), bootstrap 0.334 (-0.438 to 1.111).
    # Table 3 footnote a: 'Derived from F = e^theta_F / (1 + e^theta_F)',
    # and the F row prints 0.575 (derived); expit(0.304) = 0.5754, which
    # reproduces it. The estimated theta_F is used here rather than
    # back-transforming the rounded 0.575.

    # ------------------------------------------------------------------
    # Covariate effect. Additive on the natural scale, so it is a
    # clearance in L/h rather than a unitless exponent.
    # ------------------------------------------------------------------
    e_crcl_cl <- 4.96
    label("Additive eGFR slope theta_eGFR on CL, per unit of eGFR/85.2 (L/h)")
    # Table 3 theta_eGFR = 4.96 (RSE 16%), bootstrap 4.93 (2.70-6.78),
    # inside the equation 'CL = CLp + theta_eGFR x (eGFR/85.2), L/h'. The
    # Table 3 footnote glosses theta_eGFR as the 'increase in CL per
    # mL/min/1.73 m^2 eGFR', which does not match the printed equation --
    # per unit eGFR the slope is 4.96 / 85.2 = 0.0582 L/h, and the Results
    # text confirms the equation: "For every 1 mL/min/1.73 m^2 decrease in
    # eGFR, there was a 0.06 L/h decrease in GCV/VGCV CL." See the vignette
    # Errata.

    # ------------------------------------------------------------------
    # Interindividual variability. Table 3 block header reads
    # 'Interindividual variability (variance)', so every value below is a
    # VARIANCE and is used as-is. Results: "lognormally distributed IIV on
    # CL and V1. For F, a logit transformation was applied and a normal
    # distribution for IIV was incorporated in the logit domain."
    # The source reports no off-diagonal covariances, so the matrix is
    # diagonal here. Results also state the data were insufficient to
    # estimate IIV on kA, V2 and Q, and that inter-occasion variability on
    # CL, V1 and F did not improve the fit; neither is implemented.
    # ------------------------------------------------------------------
    etalcl ~ 0.165   # Table 3, IIV 'CL' variance 0.165 (RSE 23%); bootstrap 0.152 (0.076-0.233)
    etalvc ~ 0.431   # Table 3, IIV 'V1' variance 0.431 (RSE 48%); bootstrap 0.391 (0.103-2.301)
    etalogitfdepot ~ 2.63  # Table 3, IIV 'F' variance 2.63 (RSE 39%) in the logit domain; bootstrap 2.53 (0.70-4.85)

    # ------------------------------------------------------------------
    # Residual error. Table 3 block header reads 'Residual unexplained
    # variability (variance)', so the printed 0.236 is a variance and the
    # standard deviation nlmixr2 wants is sqrt(0.236) = 0.485798.
    # Results: "Residual variability was best explained by a proportional
    # residual error model."
    # ------------------------------------------------------------------
    propSd <- 0.485798
    label("Proportional residual error SD (fraction)")
    # Table 3 'Proportional error' variance 0.236 (RSE 11%), bootstrap
    # 0.229 (0.172-0.297); sqrt(0.236) = 0.485798.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # Clearance: the covariate enters ADDITIVELY on the natural scale, so
    # the eGFR term is built first and the log-normal interindividual
    # variability is then applied multiplicatively to the TOTAL clearance
    # (Dvorackova 2026 Table 3 equation
    # 'CL = CLp + theta_eGFR x (eGFR/85.2), L/h' with lognormal IIV on CL).
    # This deliberately does not mu-reference etalcl -- an additive
    # covariate model cannot be written in mu-referenced form.
    # ------------------------------------------------------------------
    cl <- (exp(lcl) + e_crcl_cl * (CRCL / 85.2)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    # Bioavailability is estimated in the logit domain with normally
    # distributed variability there, so the inverse logit is taken after
    # adding the eta. This keeps every individual F strictly inside (0, 1).
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Oral valganciclovir enters `depot`; intravenous
    # ganciclovir is dosed straight into `central`.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Dose scaling and lag time.
    #
    # The oral dose record is in administered VALGANCICLOVIR mg. The source
    # converted it to ganciclovir equivalents as a data step before fitting
    # -- Methods 'PK Model Development': "As GCV was given as its prodrug,
    # VGCV, we converted the VGCV dosage to a GCV equivalent using molecular
    # weights (354.362 g/mol for VGCV and 255.23 g/mol for GCV)." That
    # conversion is applied here instead, inside f(depot), so a user doses
    # the prescribed valganciclovir amount (e.g. 900 mg) and both states
    # still hold ganciclovir mg. F = 0.575 is bioavailability relative to
    # the ganciclovir-equivalent dose, so the two factors multiply.
    mw_gcv_per_vgcv <- 255.23 / 354.362
    f(depot) <- fdepot * mw_gcv_per_vgcv

    # The lag time applies to oral absorption only -- Methods: "Inclusion of
    # a lag time on the oral absorption was tested." Attaching it to `depot`
    # leaves intravenous doses into `central` unlagged.
    alag(depot) <- tlag

    # 5. Observation and residual error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
