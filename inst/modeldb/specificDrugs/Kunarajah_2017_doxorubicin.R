Kunarajah_2017_doxorubicin <- function() {
  description <- paste(
    "Population PK/PD model for IV doxorubicin (3-compartment) with",
    "first-order metabolism to doxorubicinol (1-compartment) and a",
    "cardiac troponin I (cTnI) turnover sub-model in paediatric oncology",
    "patients (Kunarajah 2017). Body surface area enters as a linear",
    "factor on every clearance and volume parameter; age enters as an",
    "additional power factor on doxorubicin clearance. The cTnI turnover",
    "sub-model is driven by a saturable Emax stimulation of cTnI",
    "synthesis by the combined doxorubicin + doxorubicinol plasma",
    "concentration, with the cTnI baseline shifted linearly by the prior",
    "cumulative anthracyclines dose received by the patient before the",
    "first dose analysed."
  )
  reference <- "Kunarajah K, Hennig S, Norris RLG, Lobb M, Charles BG, Pinkerton R, Moore AS. Population pharmacokinetic modelling of doxorubicin and doxorubicinol in children with cancer: is there a relationship with cardiac troponin profiles? Cancer Chemother Pharmacol. 2017;79(6):1209-1217. doi:10.1007/s00280-017-3309-6"
  vignette <- "Kunarajah_2017_doxorubicin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at the analysed dose. Enters as a linear factor `1 + (BSA - 1.8) * 0.465` applied identically to every clearance and volume parameter (CL, V1, Q2, V2, Q3, V3, CLm, V4, Qm) in the doxorubicin-doxorubicinol popPK model.",
      units              = "m^2",
      type               = "continuous",
      reference_category = "n/a -- linearised around BSA = 1.8 m^2 (the typical-adult reference used by the source NM-TRAN .ctl). The paper does not specify a BSA computation formula; in the validation vignette the Mosteller formula is used.",
      notes              = "Time-fixed per subject in the source dataset. Children in the cohort had BSA well below 1.8 m^2, so the linear factor evaluates to less than 1 across the population (paediatric typical CL / V are smaller than the adult-reference parameter values).",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Subject age at the analysed dose. Enters as an additional power factor `1 + (AGE/8.4)^0.736` on doxorubicin total clearance (CL only); does not modify any volume or doxorubicinol parameter.",
      units              = "years",
      type               = "continuous",
      reference_category = "n/a -- normalised by AGE / 8.4 (years).",
      notes              = "The 8.4-year reference and the 0.736 exponent were carried over from the Kontny 2013 / Voller 2015 doxorubicin maturation literature (FIXED in the Kunarajah 2017 fit, .ctl THETA(9) FIX) rather than estimated from the Kunarajah 2017 cohort.",
      source_name        = "AGE"
    ),
    PRIOR_ANTHRACYCLINE_DOSE = list(
      description        = "Cumulative prior anthracycline dose received before the first dose analysed in the current cycle, expressed in doxorubicin-equivalent body-surface-area-normalised mg/m^2. Enters as a linear shift `1 + 0.00308 * (PRIOR_ANTHRACYCLINE_DOSE - 90)` on the typical baseline cardiac troponin I (cTnI) before the next dose; does not modify any PK parameter.",
      units              = "mg/m^2",
      type               = "continuous",
      reference_category = "n/a -- linearised around the cohort median 90 mg/m^2.",
      notes              = "Time-fixed per subject for the modelled cycle (the running cumulative anthracycline dose at the first observed dose; further doses received during the modelled window are accounted for through the PK and the Emax stimulation, not through this baseline shift). Cohort range 0-225 mg/m^2 (six anthracycline-naive patients with PRIOR_ANTHRACYCLINE_DOSE = 0; eleven previously-exposed patients with median prior dose 100 mg/m^2). See `inst/references/covariate-columns.md` for the canonical entry.",
      source_name        = "PCAMT"
    )
  )

  population <- list(
    n_subjects     = 17,
    n_studies      = 1,
    age_range      = "3.42-14.67 years (median 7.50)",
    weight_range   = "11.0-88.6 kg (Table 1)",
    height_range   = "0.9-1.8 m (Table 1)",
    sex_female_pct = 29.4,
    species        = "Human (paediatric oncology)",
    disease_state  = "Childhood malignancies treated with doxorubicin (Hodgkin lymphoma, non-Hodgkin lymphoma, hepatoblastoma, Wilms' tumour, T-ALL, Pre-B ALL, synovial sarcoma, neuroblastoma, osteosarcoma; Table 1).",
    dose_range     = "Single IV doxorubicin doses 25-75 mg/m^2 with infusion durations 0.25-72.5 h (Table 1; median 30 mg/m^2).",
    regions        = "Royal Children's Hospital, Brisbane, Australia (now Lady Cilento Children's Hospital).",
    prior_anthracycline_pct = 64.7,
    prior_anthracycline_dose_range = "0-225 mg/m^2 (median 100 mg/m^2 in the 11 previously-exposed patients).",
    notes          = "19 children enrolled, 17 had blood sampled for analysis. Sampling: pre-infusion; (post-infusion) 5-10 min, 1-2 h, 2-12 h, 24-120 h (1-5 days), 168 h (7 days). 99 doxorubicin and 119 doxorubicinol concentrations and 104 cTnI measurements available for modelling. Lower limit of quantification 4.7 ng/mL for both doxorubicin and doxorubicinol; cTnI lower limit 0.04 ug/L (= 99th percentile of the assay reference range)."
  )

  ini({
    # Doxorubicin disposition. Table 2 reports parameter estimates as
    # typical values for a 1.8 m^2 reference patient (the body-surface-
    # area-normalised parameterisation used by Kontny 2013); five of the
    # nine PK fixed-effects (Q2, V2, Q3, V3, Qm) are FIX-carried from
    # Kontny 2013 / Voller 2015, the four others (CL, V1, CLm, V4) are
    # estimated in this analysis.
    lcl       <- log(58.7)  ; label("Doxorubicin total clearance CL at BSA = 1.8 m^2, AGE = 0 (L/h)")           # Kunarajah 2017 Table 2 (estimated; .ctl THETA(1))
    lvc       <- log(32.2)  ; label("Doxorubicin central volume V1 at BSA = 1.8 m^2 (L)")                       # Table 2 (estimated; THETA(2))
    lq        <- fixed(log(35.8))  ; label("Doxorubicin inter-compartmental clearance Q to peripheral1 at BSA = 1.8 m^2 (L/h); paper Q2") # Table 2 FIX (THETA(5))
    lvp       <- fixed(log(3810))  ; label("Doxorubicin peripheral1 volume Vp at BSA = 1.8 m^2 (L); paper V2")  # Table 2 FIX (THETA(6))
    lq2       <- fixed(log(65.1))  ; label("Doxorubicin inter-compartmental clearance Q2 to peripheral2 at BSA = 1.8 m^2 (L/h); paper Q3") # Table 2 FIX (THETA(7))
    lvp2      <- fixed(log(705))   ; label("Doxorubicin peripheral2 volume Vp2 at BSA = 1.8 m^2 (L); paper V3") # Table 2 FIX (THETA(8))
    lqm       <- fixed(log(32.1))  ; label("Formation clearance Qm doxorubicin -> doxorubicinol at BSA = 1.8 m^2 (L/h)") # Table 2 FIX (THETA(13))

    # Doxorubicinol disposition (single central compartment).
    lcl_doxol <- log(19.9)  ; label("Doxorubicinol elimination clearance CLm at BSA = 1.8 m^2 (L/h)")           # Table 2 (estimated; THETA(11))
    lvc_doxol <- log(508)   ; label("Doxorubicinol central volume V4 at BSA = 1.8 m^2 (L)")                     # Table 2 (estimated; THETA(12))

    # Covariate effect coefficients. Both FIXED in this analysis.
    e_age_cl    <- fixed(0.736)   ; label("Age-on-CL exponent in `1 + (AGE/8.4)^e_age_cl` (unitless)")          # Table 2 FIX (THETA(9))
    e_bsa_clvc  <- fixed(0.465)   ; label("BSA-on-all-CL/V linear coefficient in `1 + (BSA - 1.8) * e_bsa_clvc` (1/m^2); shared identically across CL, V1, Q2, V2, Q3, V3, CLm, V4, Qm") # Table 2 FIX (THETA(10))

    # Cardiac troponin I (cTnI) turnover sub-model.
    lkdeg       <- log(0.6)   ; label("cTnI first-order degradation rate constant kdeg (1/h)")                 # Table 2 (estimated; THETA(16))
    lemax       <- log(0.15)  ; label("Maximum drug-driven fractional increase in cTnI synthesis rate (Emax, unitless)") # Table 2 (estimated; THETA(19))
    lec50       <- log(11.8)  ; label("Combined doxorubicin + doxorubicinol plasma concentration giving half-maximal cTnI synthesis stimulation (EC50, ug/L)") # Table 2 (estimated; THETA(20))
    lbl_ctni    <- log(0.021) ; label("Baseline cTnI for a patient with PRIOR_ANTHRACYCLINE_DOSE = 90 mg/m^2 (Cbase, ug/L)") # Table 2 (estimated; THETA(17) reported as 20.5 pg/mL = 0.021 ug/L)
    e_pcamt_bl_ctni <- fixed(0.00308) ; label("Linear coefficient on prior cumulative anthracycline dose in `1 + e_pcamt_bl_ctni * (PRIOR_ANTHRACYCLINE_DOSE - 90)` shift on baseline cTnI (per mg/m^2)") # Table 2 (estimated; THETA(21)). Paper text: 0.31% baseline shift per 1 mg/m^2 prior dose.

    # IIV. Doxorubicin CL and doxorubicinol CL share a 2x2 NONMEM
    # $OMEGA BLOCK; QM, Q3 (q2 here), V3 (vp2 here), and the cTnI
    # baseline carry independent diagonal $OMEGA blocks.
    etalcl + etalcl_doxol ~ c(0.0344, 0.0523, 0.116)   # NONMEM $OMEGA BLOCK(2): var(CL), cov(CL,CLm), var(CLm)
    etalqm     ~ 0.11      # NONMEM $OMEGA 0.11   ; BSV QM       (var on log scale, ~33% CV)
    etalq2     ~ 0.0497    # NONMEM $OMEGA 0.0497 ; BSV Q3       (~22% CV; Q3 -> q2 in nlmixr2 convention)
    etalvp2    ~ 0.0241    # NONMEM $OMEGA 0.0241 ; BSV V3       (~16% CV; V3 -> vp2)
    etalbl_ctni ~ 0.0607   # NONMEM $OMEGA 0.0607 ; BSV TBASE    (~25% CV)

    # Residual error. Doxorubicin and doxorubicinol use combined
    # additive + proportional; cTnI uses proportional only. NONMEM
    # SIGMA = 1 with SDs carried on the THETA line; in nlmixr2 the
    # values below are the residual SDs (proportional fractions are
    # CV-on-the-individual-prediction; additive SDs are in ug/L).
    propSd       <- 0.203 ; label("Doxorubicin proportional residual SD (fraction)")                           # Table 2 (THETA(3); .ctl ;Doxol label is a transcription typo)
    addSd        <- 0.24  ; label("Doxorubicin additive residual SD (ug/L)")                                   # Table 2 (THETA(4); paper rounds 0.238 to 0.24)
    propSd_doxol <- 0.229 ; label("Doxorubicinol proportional residual SD (fraction)")                         # Table 2 (THETA(14))
    addSd_doxol  <- 0.95  ; label("Doxorubicinol additive residual SD (ug/L)")                                 # Table 2 (THETA(15); paper rounds 0.954 to 0.95)
    propSd_cTnI  <- 0.582 ; label("cTnI proportional residual SD (fraction)")                                  # Table 2 (THETA(18))
  })

  model({
    # Covariate factors. The BSA factor centres the parameters on a
    # 1.8 m^2 reference and applies identically to every clearance and
    # volume parameter (Table 2 footnote: "Same parameterisation for
    # all other clearance and distribution parameters in the model").
    # The age factor applies only to total doxorubicin CL.
    fbsacl <- 1 + (BSA - 1.8) * e_bsa_clvc
    fage   <- 1 + (AGE / 8.4)^e_age_cl
    fpcamt <- 1 + e_pcamt_bl_ctni * (PRIOR_ANTHRACYCLINE_DOSE - 90)

    # Doxorubicin individual parameters.
    cl  <- exp(lcl  + etalcl)  * fbsacl * fage
    vc  <- exp(lvc)            * fbsacl
    q   <- exp(lq)             * fbsacl
    vp  <- exp(lvp)            * fbsacl
    q2  <- exp(lq2  + etalq2)  * fbsacl
    vp2 <- exp(lvp2 + etalvp2) * fbsacl
    qm  <- exp(lqm  + etalqm)  * fbsacl

    # Doxorubicinol individual parameters.
    cl_doxol <- exp(lcl_doxol + etalcl_doxol) * fbsacl
    vc_doxol <- exp(lvc_doxol)                * fbsacl

    # cTnI turnover individual parameters.
    kdeg    <- exp(lkdeg)
    emax    <- exp(lemax)
    ec50    <- exp(lec50)
    bl_ctni <- exp(lbl_ctni + etalbl_ctni) * fpcamt
    ksyn    <- bl_ctni * kdeg

    # Doxorubicin: 3-compartment IV disposition. The total clearance
    # CL already includes the doxorubicin -> doxorubicinol formation
    # pathway, so the parent-loss term in d/dt(central) uses cl/vc and
    # the doxorubicinol input term uses qm/vc separately (Kunarajah
    # 2017 Appendix $DES). Mass balance: cl >= qm.
    d/dt(central)     <- -cl / vc * central -
                          q  / vc * central + q  / vp  * peripheral1 -
                          q2 / vc * central + q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2

    # Doxorubicinol: 1-compartment fed by qm * central / vc.
    d/dt(central_doxol) <- qm / vc * central -
                           cl_doxol / vc_doxol * central_doxol

    # Plasma concentrations in ug/L. The compartmental amounts central
    # and central_doxol are in mg (matching the dose units), so dividing
    # by the corresponding volume (L) gives mg/L; the *1000 converts to
    # the paper's reporting units (ug/L = ng/mL).
    Cc       <- (central       / vc)       * 1000
    Cc_doxol <- (central_doxol / vc_doxol) * 1000

    # cTnI: indirect-response turnover with combined doxorubicin +
    # doxorubicinol Emax stimulation of zero-order synthesis.
    drug_effect <- emax * (Cc + Cc_doxol) / (ec50 + Cc + Cc_doxol)
    effect(0)   <- bl_ctni
    d/dt(effect) <- ksyn * (1 + drug_effect) - kdeg * effect
    cTnI <- effect

    Cc       ~ add(addSd)       + prop(propSd)
    Cc_doxol ~ add(addSd_doxol) + prop(propSd_doxol)
    cTnI     ~ prop(propSd_cTnI)
  })
}
