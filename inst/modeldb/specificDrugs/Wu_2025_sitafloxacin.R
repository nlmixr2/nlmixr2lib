Wu_2025_sitafloxacin <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral sitafloxacin built by Wu 2025 from 12 pooled ",
    "clinical trials (342 Japanese and Chinese healthy volunteers, elderly volunteers, subjects ",
    "with renal impairment and patients with respiratory-tract infection; 3,294 plasma ",
    "concentrations) to drive Monte Carlo PK/PD cut-off (clinical-breakpoint) determination ",
    "against Streptococcus pneumoniae, Staphylococcus aureus, Escherichia coli, Klebsiella ",
    "pneumoniae and Pseudomonas aeruginosa. Absorption is sequential: zero-order input into the ",
    "depot over duration D1, then first-order transfer to central at ka, behind an absorption lag ",
    "TLAG. Creatinine clearance is a power covariate on apparent clearance; body weight and age ",
    "are power covariates on the apparent central volume; the fed state lengthens the zero-order ",
    "absorption duration. All disposition parameters are apparent (/F): only oral data were ",
    "analysed and absolute bioavailability was not estimated."
  )
  reference <- paste(
    "Wu H, Li Y, Li X, Fan Y, Guo B, Liu X, Li W, Chen M, Chen Y, Zhang J.",
    "Model-guided development of pharmacokinetic/pharmacodynamic cut-offs and evaluation of",
    "sitafloxacin dosing regimens against target pathogens.",
    "Front Pharmacol. 2025;16:1476158. doi:10.3389/fphar.2025.1476158.",
    "See also modellib('Rodjun_2023_sitafloxacin') for the independent one-compartment",
    "sitafloxacin model that Rodjun 2023 carried over from Tanigawara 2013.",
    sep = " "
  )
  vignette <- "Wu_2025_sitafloxacin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "sitafloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "sitafloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sitafloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters apparent clearance as the power term (CRCL / 106.88)^0.460, exactly as printed in ",
        "the Wu 2025 Section 2.2 display equation 'CL = 14.7 x (CRCL/106.88)^0.46 x exp(eta1)'. ",
        "The normalizing value 106.88 mL/min is the reference the authors printed inside the ",
        "equation; it is close to but distinct from the Table 1 cohort MEAN of 115.1 +/- 53.5 ",
        "mL/min, so it is most likely the cohort median (Wu 2025 does not label it). Units are raw ",
        "mL/min, matching the Table 1 unit column; no BSA normalization is applied anywhere in the ",
        "paper. The PK/PD simulations in Wu 2025 Section 4.2 were restricted to normal renal ",
        "function and mild renal impairment (CRCL >= 50 mL/min)."
      ),
      source_name        = "CRCL"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters the apparent central volume as the power term (WT / 58.55)^0.966, per the Wu 2025 ",
        "Section 2.2 display equation 'V2 = 89.8 x (WT/58.55)^0.966 x (AGE/31)^0.286 x exp(eta2)'. ",
        "The exponent is estimated (0.966, RSE 13.5%), not fixed at an allometric 1. The ",
        "normalizing value 58.55 kg is printed inside the equation and is close to the Table 1 ",
        "cohort mean of 58.2 +/- 10.2 kg. Weight does NOT scale clearance in this model."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters the apparent central volume as the power term (AGE / 31)^0.286, per the Wu 2025 ",
        "Section 2.2 display equation for V2. The normalizing value 31 years is printed inside the ",
        "equation and is well below the Table 1 cohort mean of 43.0 +/- 21.7 years, consistent ",
        "with a median dominated by the 183 healthy volunteers. The positive exponent means older ",
        "subjects have a larger apparent central volume."
      ),
      source_name        = "AGE"
    ),
    FED = list(
      description        = "Fed state at the time of dosing",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (fasting)",
      notes              = paste0(
        "1 = dose taken in the fed state, 0 = fasting. Wu 2025 calls the column FOOD and applies ",
        "it to the zero-order absorption duration in the unusual multiplicative POWER-OF-(1+FOOD) ",
        "form printed in the Section 2.2 display equation: 'D1 = 0.281 x (1 + FOOD)^1.59 x ",
        "exp(eta3)'. This is NOT the usual exp(theta * FOOD) indicator form, and it is encoded ",
        "literally in model() rather than being rewritten. The paper's own worked values confirm ",
        "both the form and the orientation of the indicator: the same display equation states ",
        "'in fasting state D1(h) = 0.281' and 'in fed state D1(h) = 0.846', and ",
        "0.281 * 2^1.59 = 0.8457, which rounds to the printed 0.846. Food therefore lengthens the ",
        "zero-order absorption window roughly 3-fold. Wu 2025 does not report the meal ",
        "composition or caloric content, so the generic FED indicator is used rather than the ",
        "FED_HIGHFAT / FED_LOWFAT refinements."
      ),
      source_name        = "FOOD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 342L,
    n_observations = 3294L,
    n_studies      = 12L,
    age_range      = "Mean (SD) 43.0 (21.7) years (Wu 2025 Table 1). Individual minimum and maximum are not reported; the studies included a dedicated elderly-volunteer arm.",
    weight_range   = "Mean (SD) 58.2 (10.2) kg (Wu 2025 Table 1). Mean (SD) height 164.7 (9.6) cm and BMI 21.4 (3.0) kg/m^2.",
    sex_female_pct = 28.1,
    race_ethnicity = "Not reported as such. 294 subjects were enrolled in Japan (147 patients with respiratory-system infection, 12 subjects with renal insufficiency, 135 healthy subjects) and 48 healthy subjects were enrolled in China (Wu 2025 Section 2.1).",
    disease_state  = "Pooled: 147 patients with respiratory-system infection, 12 subjects with renal insufficiency, and 183 healthy subjects (Wu 2025 Section 2.1). The paper separately notes 24 subjects with moderate renal impairment and 4 with severe renal impairment in the dataset (Wu 2025 Discussion).",
    dose_range     = "Single doses of 3-200 mg or multiple doses of 50 or 100 mg q12h or 100 mg q8h for 7-14 days, all oral (Wu 2025 Section 4.1). The three regimens carried forward into the PK/PD simulations were 50 mg q12h, 100 mg q24h and 100 mg q12h.",
    regions        = "Japan (11 studies) and China (1 study).",
    renal_function = "Mean (SD) creatinine clearance 115.1 (53.5) mL/min (Wu 2025 Table 1). The dataset spans normal renal function through severe impairment, but the published PK/PD cut-offs apply only to CRCL >= 50 mL/min.",
    notes          = "Fitted in NONMEM 7.4 with FOCE-I (Wu 2025 Section 4.1). 1.02% of concentration records were missing or excluded. Model qualification: 1,000-replicate bootstrap (945 successful, 94.5% convergence) and a 1,000-replicate visual predictive check."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural PK -- Wu 2025 Table 2, 'Model estimates Mean (RSE%)' column,
    # cross-checked against the Section 2.2 display equations. All disposition
    # parameters are apparent (/F): Wu 2025 analysed oral data only and did not
    # estimate absolute bioavailability, so F is implicitly 1 and no lfdepot
    # term is used. NONMEM compartment numbering in the source is ADVAN4-style
    # (1 = depot, 2 = central, 3 = peripheral), which is what fixes V2 as the
    # CENTRAL volume, V3 as the PERIPHERAL volume, and D1 / TLAG as the
    # zero-order duration and lag time of the DEPOT compartment.
    # -----------------------------------------------------------------------
    lcl   <- log(14.7)  ; label("Apparent total clearance CL/F (L/h)")                              # Wu 2025 Table 2 CL = 14.7 L/h (RSE 1.8%; bootstrap mean 14.6, 95% CI 14.0-15.2)
    lvc   <- log(89.8)  ; label("Apparent central volume of distribution V2/F (L)")                 # Wu 2025 Table 2 V2 = 89.8 L (RSE 2.8%; bootstrap mean 89.6, 95% CI 84.2-95.6)
    lq    <- log(5.81)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Wu 2025 Table 2 Q = 5.81 L/h (RSE 6.9%; bootstrap mean 5.73, 95% CI 3.83-7.64)
    lvp   <- log(33.8)  ; label("Apparent peripheral volume of distribution V3/F (L)")              # Wu 2025 Table 2 V3 = 33.8 L (RSE 3.6%; bootstrap mean 33.6, 95% CI 30.5-36.5)
    lka   <- log(5.19)  ; label("First-order absorption rate constant ka (1/h)")                    # Wu 2025 Table 2 Ka = 5.19 1/h (RSE 16.3%; bootstrap mean 5.16, 95% CI 4.02-6.44)
    ld1   <- log(0.281) ; label("Zero-order absorption input duration D1, fasting (h)")             # Wu 2025 Table 2 D1 = 0.281 h (RSE 4.2%; bootstrap mean 0.278, 95% CI 0.211-0.309); Section 2.2 'in fasting state D1(h) = 0.281'
    ltlag <- log(0.205) ; label("Absorption lag time TLAG (h)")                                     # Wu 2025 Table 2 TLAG = 0.205 h (RSE 1.1%; bootstrap mean 0.207, 95% CI 0.196-0.219)

    # -----------------------------------------------------------------------
    # Covariate effects. Wu 2025 prints the final covariate model as display
    # equations in Section 2.2, immediately above the Table 2 abbreviation
    # note:
    #
    #   CL = 14.7 * (CRCL/106.88)^0.46  * exp(eta1)                (L/h)
    #   V2 = 89.8 * (WT/58.55)^0.966 * (AGE/31)^0.286 * exp(eta2)  (L)
    #   V3 = 33.8                                                  (L)
    #   Q  = 5.81                                                  (L/h)
    #   D1 = 0.281 * (1 + FOOD)^1.59 * exp(eta3)                   (h)
    #   TLAG = 0.205                                               (h)
    #
    # Every coefficient below is the printed equation coefficient, and each
    # matches its Table 2 row to the digits Table 2 prints (Table 2 gives
    # 0.460 where the equation rounds to 0.46). The continuous covariates are
    # already in power form in the source, so no exp()/power rewriting is
    # needed. The food effect is left in its literal (1 + FOOD)^theta form --
    # see the FED entry in covariateData for why, and for the arithmetic that
    # confirms it against the paper's own printed fed-state value of 0.846 h.
    # -----------------------------------------------------------------------
    e_crcl_cl <- 0.460 ; label("Power effect of creatinine clearance on CL/F (unitless, reference 106.88 mL/min)") # Wu 2025 Table 2 'CRCL on CL' = 0.460 (RSE 6.5%; bootstrap mean 0.486, 95% CI 0.417-0.570); Section 2.2 equation '(CRCL/106.88)^0.46'
    e_wt_vc   <- 0.966 ; label("Power effect of body weight on V2/F (unitless, reference 58.55 kg)")               # Wu 2025 Table 2 'WT on V2' = 0.966 (RSE 13.5%; bootstrap mean 0.973, 95% CI 0.711-1.27); Section 2.2 equation '(WT/58.55)^0.966'
    e_age_vc  <- 0.286 ; label("Power effect of age on V2/F (unitless, reference 31 years)")                       # Wu 2025 Table 2 'Age on V2' = 0.286 (RSE 17.9%; bootstrap mean 0.260, 95% CI 0.114-0.372); Section 2.2 equation '(AGE/31)^0.286'
    e_fed_d1  <- 1.59  ; label("Fed-state exponent on D1, applied as (1 + FED)^theta (unitless)")                  # Wu 2025 Table 2 'Food on D1' = 1.59 (RSE 3.8%; bootstrap mean 1.52, 95% CI 0.898-1.83); Section 2.2 equation '(1 + FOOD)^1.59'

    # -----------------------------------------------------------------------
    # Inter-individual variability. Wu 2025 Table 2 reports these under the
    # heading "Interindividual variability (omega^2)", so the tabulated
    # numbers are VARIANCES and are used directly (no CV%-to-variance
    # conversion). The Table 2 abbreviation footnote assigns them:
    # "eta1, interindividual variability of CL; eta2, interindividual
    # variability of V2; eta3, interindividual variability of D1;
    # eta4, interindividual variability of IOV."
    #
    # eta3 on D1 is genuinely enormous as published (omega^2 = 3.46, i.e.
    # omega = 1.86 on the log scale) with 24.9% shrinkage. It is encoded here
    # exactly as reported rather than being tempered: a 3-fold food effect on
    # D1 plus pooled fasted/fed dosing across 12 studies is consistent with an
    # absorption-duration term that is poorly identified in individual
    # subjects. The consequence for simulation is discussed in the vignette
    # Errata -- because D1 only shapes the absorption window and the dose mass
    # is conserved, even extreme eta3 draws leave AUC untouched and perturb
    # only Cmax and Tmax.
    #
    # ka, Q/F, V3/F and TLAG carry no IIV term in Wu 2025 Table 2 and none is
    # invented here (the model is encoded exactly as published).
    #
    # !! eta4 IS NOT ENCODED -- OPEN ITEM. Wu 2025 Table 2 reports a fourth
    # random effect with omega^2 = 0.0279 (RSE 8.7%, shrinkage 29.1%,
    # bootstrap mean 0.0284, 95% CI 0.0184-0.0394), described in the
    # abbreviation footnote only as "interindividual variability of IOV" --
    # i.e. an inter-occasion variability term. The paper never states which
    # structural parameter carries it, and never defines what an occasion is.
    # The Section 2.2 display equations account for eta1, eta2 and eta3 only;
    # eta4 appears nowhere in them. No supplement or control stream is
    # available (EuropePMC returns HTTP 500 for this PMCID's
    # supplementaryFiles endpoint and the Frontiers file endpoints 404), and
    # the supplement's documented contents are a study-design table and mean
    # concentration-time figures, neither of which would place an eta.
    # Assigning it to CL, D1 or F would be inventing model structure, so the
    # term is omitted and its reported value preserved here instead. See the
    # vignette Errata.
    # -----------------------------------------------------------------------
    etalcl ~ 0.0515 ; label("IIV on apparent clearance (variance, log scale)")                    # Wu 2025 Table 2 eta1 omega^2 = 0.0515 (RSE 5.5%, shrinkage 14.7%; bootstrap mean 0.0475, 95% CI 0.0289-0.0621)
    etalvc ~ 0.0635 ; label("IIV on apparent central volume (variance, log scale)")               # Wu 2025 Table 2 eta2 omega^2 = 0.0635 (RSE 9.3%, shrinkage 29.6%; bootstrap mean 0.0585, 95% CI 0.0359-0.0741)
    etald1 ~ 3.46   ; label("IIV on zero-order absorption duration (variance, log scale)")        # Wu 2025 Table 2 eta3 omega^2 = 3.46 (RSE 7.8%, shrinkage 24.9%; bootstrap mean 3.46, 95% CI 2.62-4.43)

    # -----------------------------------------------------------------------
    # Residual error: combined proportional plus additive. Wu 2025 Table 2
    # reports both under the heading "Residual variability (sigma^2)", so the
    # tabulated numbers are VARIANCES and the square roots are entered here as
    # standard deviations:
    #   propSd = sqrt(0.0361)  = 0.19        (an exact square root, which
    #                                         independently confirms the
    #                                         variance reading -- 19% CV)
    #   addSd  = sqrt(2.35e-5) = 0.0048477   mg/L (= 4.85 ng/mL)
    # The Table 2 abbreviation footnote assigns them: "epsilon1, proportional
    # residual error; epsilon2, additive residuals". Concentration units are
    # mg/L because dose is in mg and V2/F is in L.
    # -----------------------------------------------------------------------
    propSd <- 0.19      ; label("Proportional residual error (fraction)")   # Wu 2025 Table 2 epsilon1 sigma^2 = 0.0361 (RSE 1.5%, shrinkage 10.3%; bootstrap mean 0.0362, 95% CI 0.0290-0.0429) -> sqrt = 0.19
    addSd  <- 0.0048477 ; label("Additive residual error (mg/L)")           # Wu 2025 Table 2 epsilon2 sigma^2 = 2.35e-5 (RSE 6.5%, shrinkage 10.3%; bootstrap mean 2.24e-5, 95% CI 8.20e-6-3.57e-5) -> sqrt = 0.0048477
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # CL/F carries the creatinine-clearance power term; V2/F carries the
    # body-weight and age power terms; D1 carries the fed-state term. Q/F,
    # V3/F, ka and TLAG are covariate-free. At the reference covariate vector
    # (CRCL 106.88 mL/min, WT 58.55 kg, AGE 31 years, FED 0) every multiplier
    # collapses to 1 and the parameters reduce to the Table 2 typical values.
    # ---------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 106.88)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 58.55)^e_wt_vc * (AGE / 31)^e_age_vc
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1) * (1 + FED)^e_fed_d1
    tlag <- exp(ltlag)

    # 2. Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # 3. Two-compartment ODE system with sequential zero-order then
    #    first-order absorption (Wu 2025 Section 2.2: "The absorption of
    #    sitafloxacin involves zero-order and first-order kinetics", with the
    #    display equation listing both a zero-order D1 and a first-order Ka
    #    under a single "Absorption" brace). The dose enters `depot` at a
    #    constant rate over the zero-order window of duration D1, starting
    #    TLAG after the dose record, and `depot` then drains into `central`
    #    first-order at ka.
    #
    #    Because D1 is a MODELLED duration, dose records must carry
    #    rate = -2 so that rxode2 uses dur(depot); a plain bolus would
    #    collapse the zero-order phase and bias Cmax upward.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot)  <- d1
    alag(depot) <- tlag

    # ---------------------------------------------------------------------
    # 4. Observation. Dose is in mg and vc in L, so central / vc is mg/L
    #    (= ug/mL), which is the unit Wu 2025 uses for the MIC and PK/PD
    #    cut-off scale throughout (Tables 3 and 4).
    #
    #    Cc is the TOTAL plasma concentration. Wu 2025's PK/PD index is the
    #    UNBOUND fAUC24h/MIC, but the paper never reports a sitafloxacin
    #    unbound fraction, so no fu term is introduced here. The validation
    #    vignette applies fu = 0.388 from the separate, already-packaged
    #    Rodjun_2023_sitafloxacin model (Rodjun 2023, sourced there to
    #    Tanigawara 2013) purely to reproduce the published PK/PD cut-offs,
    #    and flags it as non-Wu provenance.
    # ---------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
