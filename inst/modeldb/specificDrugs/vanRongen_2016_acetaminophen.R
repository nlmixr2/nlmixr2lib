vanRongen_2016_acetaminophen <- function() {
  description <- "Parent-and-metabolites population PK model for intravenous acetaminophen (paracetamol) and its glucuronide, sulphate, and CYP2E1-oxidation (cysteine + mercapturate) metabolites in morbidly obese and non-obese adults (van Rongen 2016). One-compartment plasma disposition for parent acetaminophen with four parallel elimination pathways from the central compartment (glucuronidation, sulphation, CYP2E1 oxidation, and unchanged renal); one-compartment plasma disposition for glucuronide and cysteine + mercapturate metabolites each fed via a single-transit-compartment delay; two-compartment plasma disposition for sulphate (central + peripheral, fixed equal volumes 5.66 L each). Lean body weight (LBW; Janmahasatian et al. 2005 equation) enters as a power-law covariate on parent V, all three formation clearances, the CYP2E1 transit rate constant, and glucuronide elimination CL. Total body weight enters on the glucuronide volume of distribution."
  reference <- paste(
    "van Rongen A, Valitalo PAJ, Peeters MYM, Boerma D, Huisman FW,",
    "van Ramshorst B, van Dongen EPA, van den Anker JN, Knibbe CAJ (2016).",
    "Morbidly obese patients exhibit increased CYP2E1-mediated",
    "oxidation of acetaminophen.",
    "Clin Pharmacokinet 55(7):833-847.",
    "doi:10.1007/s40262-015-0357-0.",
    sep = " "
  )
  vignette <- "vanRongen_2016_acetaminophen"
  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling on glucuronide volume of distribution with reference TBW = 130.9 kg (pooled population median). Time-fixed at baseline. Source column 'TBW' renamed to canonical 'WT' on input.",
      source_name        = "TBW"
    ),
    LBM = list(
      description        = "Lean body mass at baseline (Janmahasatian et al. 2005 formula).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-law scaling with reference LBW = 65.2 kg (pooled population median) on parent V, the three formation clearances (CL_gluc, CL_sulph, CL_CYP2E1), the CYP2E1 transit rate constant, and the glucuronide elimination CL. Source column 'LBW' renamed to canonical 'LBM' on input (same biological quantity, no value transformation).",
      source_name        = "LBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 1L,
    age_range      = "18-58 years",
    age_median     = "41 years (pooled)",
    weight_range   = "53.4-193.1 kg (TBW); 36.0-96.2 kg (LBW)",
    weight_median  = "130.9 kg TBW / 65.2 kg LBW (pooled medians)",
    sex_female_pct = 67.9,
    race_ethnicity = NA_character_,
    disease_state  = "20 morbidly obese patients (BMI 40-55.2 kg/m^2) undergoing bariatric surgery + 8 non-obese patients undergoing other elective surgery; all received a single 2 g intravenous acetaminophen study dose followed by standard postoperative pain protocol (1 g intravenous acetaminophen every 6 h up to 24 h).",
    dose_range     = "Single 2 g (~13231 umol) intravenous acetaminophen infusion over 15 minutes. Optional follow-on standard-of-care dosing 1 g every 6 h up to 24 h; data after the 2 g infusion period inform the IIV but not the structural model.",
    regions        = "Netherlands (St Antonius Hospital, Nieuwegein)",
    notes          = "Demographics from Table 1 of van Rongen 2016. Pooled medians of LBW (65.2 kg) and TBW (130.9 kg) are the references used in the power-law covariate equations of the final model (Equations on page 839 and Table 2)."
  )

  ini({
    # Final population PK parameters from Table 2 of van Rongen 2016 (page
    # 839-840), final-model column. Bootstrap 95% CIs are reported in the
    # right-most column and confirm the point estimates. Concentrations were
    # logarithmically transformed before estimation, and a proportional
    # residual error model was used (Methods page 836); the reported CV% in
    # the residual-error section of Table 2 maps directly to nlmixr2's
    # propSd in linear space.

    # Structural disposition - parent acetaminophen central volume.
    lvc        <- log(67.2)   ; label("Acetaminophen central volume at LBW = 65.2 kg (L)")                                 # Table 2 final V_acetaminophen = V_{65.2 kg} = 67.2 L (RSE 2.8%)

    # Formation clearances from the central acetaminophen compartment to
    # each metabolite pool. Reference values are at the pooled population
    # median LBW (65.2 kg). The unchanged-renal arm CL_unchanged is held
    # at 5% of the typical metabolic CL sum at the reference LBW (paper
    # Methods page 836); the value 0.0163 L/min is reproduced from
    # 0.05/0.95 * (0.224 + 0.065 + 0.021).
    lcl_gluc   <- log(0.224)  ; label("Glucuronidation formation clearance at LBW = 65.2 kg (L/min)")                      # Table 2 final CL_{gluc,65.2 kg} = 0.224 L/min (RSE 5%)
    lcl_sulf   <- log(0.065)  ; label("Sulphation formation clearance at LBW = 65.2 kg (L/min)")                           # Table 2 final CL_{sulph,65.2 kg} = 0.065 L/min (RSE 6%)
    lcl_cysmer <- log(0.021)  ; label("CYP2E1 (cysteine + mercapturate) formation clearance at LBW = 65.2 kg (L/min)")     # Table 2 final CL_{CYP2E1,65.2 kg} = 0.021 L/min (RSE 14.6%)
    lcl_renal  <- fixed(log(0.05 / 0.95 * (0.224 + 0.065 + 0.021))) ; label("Renal CL of unchanged acetaminophen (L/min)") # Methods page 836: "CL_unchanged was assumed to be 5% of CL_tot (= CL_unchanged + CL_gluc + CL_sulph + CL_CYP2E1) of a 70 kg individual"

    # Glucuronide pathway - transit compartment + plasma metabolite.
    lvc_gluc   <- log(32.3)   ; label("Glucuronide volume of distribution at TBW = 130.9 kg (L)")                          # Table 2 final V_{glucuronide,130.9 kg} = 32.3 L (RSE 4.1%)
    lktr_gluc  <- log(0.095)  ; label("Glucuronide transit-compartment rate constant Ktr_gluc (1/min)")                    # Table 2 final Ktr_{gluc} = 0.095 1/min (RSE 11.5%), corresponding to MTT = 10.5 min for n = 1 transit compartment
    lcle_gluc  <- log(0.222)  ; label("Glucuronide elimination clearance at LBW = 65.2 kg (L/min)")                        # Table 2 final CLE_{gluc,65.2 kg} = 0.222 L/min (RSE 6.3%)

    # Sulphate pathway - two-compartment plasma disposition with equalized
    # central + peripheral volumes (5.66 L FIX per Methods page 836; cites
    # Liukas 2011 / ref [31]).
    lvc_sulf   <- fixed(log(5.66))  ; label("Sulphate central volume of distribution (L, fixed)")                          # Methods page 836: V_sulphate,central = V_sulphate,peripheral = 5.66 L FIX, citing Liukas 2011 (ref 31)
    lvp_sulf   <- fixed(log(5.66))  ; label("Sulphate peripheral volume of distribution (L, fixed equal to central)")      # Methods page 836: V_sulphate,central = V_sulphate,peripheral
    lq_sulf    <- log(0.339) ; label("Sulphate inter-compartmental clearance Q (L/min)")                                   # Table 2 final Q = 0.339 L/min (RSE 19.6%)
    lcle_sulf  <- log(0.096) ; label("Sulphate elimination clearance (L/min)")                                             # Table 2 final CLE_{sulph} = 0.096 L/min (RSE 3.4%)

    # CYP2E1 (cysteine + mercapturate) pathway - transit compartment +
    # plasma metabolite pool (fixed apparent volume per Liukas 2011).
    lvc_cysmer <- fixed(log(15.6))  ; label("Cysteine + mercapturate combined volume of distribution (L, fixed)")          # Methods page 836: V_{cysteine and mercapturate} = 15.6 L FIX, citing Liukas 2011 (ref 31)
    lktr_cysmer <- log(0.0057)      ; label("CYP2E1 transit-compartment rate constant Ktr_CYP2E1 at LBW = 65.2 kg (1/min)") # Table 2 final Ktr_{CYP2E1,65.2 kg} = 0.0057 1/min (RSE 12.2%), MTT = 175.4 min for n = 1 transit compartment
    lcle_cysmer <- log(0.329)       ; label("Cysteine + mercapturate elimination clearance (L/min)")                       # Table 2 final CLE_{cys and mercap} = 0.329 L/min (RSE 14.5%)

    # Covariate effects - power-law exponents from Table 2 of van Rongen
    # 2016. Six are on LBW (parent V, three formation CLs, CYP2E1 transit
    # rate constant, and glucuronide elimination CL); one is on TBW (the
    # glucuronide volume of distribution). All exponents are fitted with
    # RSE reported in Table 2, so they are estimated (no fixed()).
    e_lbm_vc          <- 0.90  ; label("Power exponent of LBW on parent V")                                                # Table 2 final S = 0.90 (RSE 17.4%)
    e_lbm_cl_gluc     <- 1.33  ; label("Power exponent of LBW on glucuronidation formation CL")                            # Table 2 final T = 1.33 (RSE 17%)
    e_lbm_cl_sulf     <- 0.92  ; label("Power exponent of LBW on sulphation formation CL")                                 # Table 2 final U = 0.92 (RSE 19.9%)
    e_lbm_cl_cysmer   <- 0.67  ; label("Power exponent of LBW on CYP2E1 formation CL")                                     # Table 2 final W = 0.67 (RSE 27.4%)
    e_wt_vc_gluc      <- 0.55  ; label("Power exponent of TBW on glucuronide volume of distribution")                      # Table 2 final X = 0.55 (RSE 23.3%)
    e_lbm_ktr_cysmer  <- 1.1   ; label("Power exponent of LBW on CYP2E1 transit rate constant")                            # Table 2 final Y = 1.1 (RSE 33%)
    e_lbm_cle_gluc    <- 0.89  ; label("Power exponent of LBW on glucuronide elimination CL")                              # Table 2 final Z = 0.89 (RSE 31%)

    # Inter-individual variability. Table 2 reports IIV as CV% on the
    # log-normal scale; the internal variance is omega^2 = log(1 + CV^2).
    # No IIV is reported on Ktr_gluc, CL_E,sulph, Q_sulph, or fixed
    # volumes (Vc_sulf, Vp_sulf, Vc_cysmer); CL_unchanged is held fixed
    # and carries no IIV.
    etalvc         ~ 0.02053  # log(1 + 0.144^2) -- Table 2 final IIV V_acetaminophen = 14.4% (RSE 32.1%)
    etalcl_gluc    ~ 0.04642  # log(1 + 0.218^2) -- Table 2 final IIV CL_gluc = 21.8% (RSE 32.5%)
    etalcl_sulf    ~ 0.05738  # log(1 + 0.243^2) -- Table 2 final IIV CL_sulph = 24.3% (RSE 30.1%)
    etalcl_cysmer  ~ 0.05286  # log(1 + 0.233^2) -- Table 2 final IIV CL_CYP2E1 = 23.3% (RSE 37.4%)
    etalvc_gluc    ~ 0.04938  # log(1 + 0.225^2) -- Table 2 final IIV V_glucuronide = 22.5% (RSE 29.5%)
    etalcle_gluc   ~ 0.08785  # log(1 + 0.303^2) -- Table 2 final IIV CLE_gluc = 30.3% (RSE 23.9%)
    etalcle_cysmer ~ 0.11491  # log(1 + 0.349^2) -- Table 2 final IIV CLE_cys&mercap = 34.9% (RSE 33.4%)

    # Residual error - proportional in linear space for every output. The
    # paper estimated on log-transformed concentrations with proportional
    # residual error (Methods page 836); the CV% reported in Table 2's
    # residual-variability section is the proportional residual SD on the
    # linear scale.
    propSd        <- 0.171  ; label("Acetaminophen proportional residual SD (fraction)")             # Table 2 final proportional error for acetaminophen = 17.1% (RSE 27%)
    propSd_gluc   <- 0.197  ; label("Glucuronide proportional residual SD (fraction)")                # Table 2 final proportional error for glucuronide = 19.7% (RSE 27.9%)
    propSd_sulf   <- 0.185  ; label("Sulphate proportional residual SD (fraction)")                   # Table 2 final proportional error for sulphate = 18.5% (RSE 20.6%)
    propSd_cysmer <- 0.250  ; label("Cysteine + mercapturate proportional residual SD (fraction)")    # Table 2 final proportional error for cys and mercap = 25.0% (RSE 8.7%)
  })

  model({
    # Individual structural parameters. LBW enters as a power-law
    # multiplier on the parent V, three formation CLs, the CYP2E1
    # transit rate, and the glucuronide elimination CL with reference
    # 65.2 kg; TBW enters on the glucuronide V with reference 130.9 kg.
    # CL_unchanged is held at a fixed typical value (no covariate, no
    # IIV) per Methods page 836.
    vc         <- exp(lvc          + etalvc)         * (LBM / 65.2)^e_lbm_vc
    cl_gluc    <- exp(lcl_gluc     + etalcl_gluc)    * (LBM / 65.2)^e_lbm_cl_gluc
    cl_sulf    <- exp(lcl_sulf     + etalcl_sulf)    * (LBM / 65.2)^e_lbm_cl_sulf
    cl_cysmer  <- exp(lcl_cysmer   + etalcl_cysmer)  * (LBM / 65.2)^e_lbm_cl_cysmer
    cl_renal   <- exp(lcl_renal)
    vc_gluc    <- exp(lvc_gluc     + etalvc_gluc)    * (WT  / 130.9)^e_wt_vc_gluc
    ktr_gluc   <- exp(lktr_gluc)
    cle_gluc   <- exp(lcle_gluc    + etalcle_gluc)   * (LBM / 65.2)^e_lbm_cle_gluc
    vc_sulf    <- exp(lvc_sulf)
    vp_sulf    <- exp(lvp_sulf)
    q_sulf     <- exp(lq_sulf)
    cle_sulf   <- exp(lcle_sulf)
    vc_cysmer  <- exp(lvc_cysmer)
    ktr_cysmer <- exp(lktr_cysmer) * (LBM / 65.2)^e_lbm_ktr_cysmer
    cle_cysmer <- exp(lcle_cysmer  + etalcle_cysmer)

    # Micro-constants. Each formation rate is the partial CL of the
    # parent divided by the parent V; the elimination rates of each
    # metabolite are the metabolite CL_E divided by the metabolite V.
    # The sulphate distribution rates use the inter-compartmental CL
    # Q and the (equal) central and peripheral volumes.
    k_form_gluc   <- cl_gluc   / vc
    k_form_sulf   <- cl_sulf   / vc
    k_form_cysmer <- cl_cysmer / vc
    k_renal       <- cl_renal  / vc
    kel_gluc      <- cle_gluc  / vc_gluc
    kel_sulf      <- cle_sulf  / vc_sulf
    kel_cysmer    <- cle_cysmer/ vc_cysmer
    k_sulf_cp     <- q_sulf    / vc_sulf
    k_sulf_pc     <- q_sulf    / vp_sulf

    # ODE system matching Figure 1 of van Rongen 2016 (Methods page 837).
    # Parent acetaminophen has four parallel elimination arms; the
    # glucuronide and CYP2E1 (cysteine + mercapturate) arms route
    # through a single transit compartment each to capture the metabolite
    # formation delay reported in Methods page 836 ("To capture eventual
    # delay in formation ... a varying number of transit compartments
    # was tested ... n = 1 transit compartment").
    d/dt(central)          <- -(k_form_gluc + k_form_sulf + k_form_cysmer + k_renal) * central
    d/dt(transit1_gluc)    <-  k_form_gluc   * central          - ktr_gluc   * transit1_gluc
    d/dt(central_gluc)     <-  ktr_gluc      * transit1_gluc    - kel_gluc   * central_gluc
    d/dt(central_sulf)     <-  k_form_sulf   * central          - kel_sulf   * central_sulf -
                               k_sulf_cp     * central_sulf     + k_sulf_pc  * peripheral1_sulf
    d/dt(peripheral1_sulf) <-  k_sulf_cp     * central_sulf     - k_sulf_pc  * peripheral1_sulf
    d/dt(transit1_cysmer)  <-  k_form_cysmer * central          - ktr_cysmer * transit1_cysmer
    d/dt(central_cysmer)   <-  ktr_cysmer    * transit1_cysmer  - kel_cysmer * central_cysmer

    # Plasma concentrations. Dose is in micromoles (umol), volumes in
    # litres, so amount / volume returns umol/L, matching van Rongen
    # 2016 which expressed all concentrations in umol/L (Methods page
    # 836: "Concentrations were expressed in micromoles per litre,
    # using the molecular weights of acetaminophen, acetaminophen
    # glucuronide, acetaminophen sulphate, acetaminophen cysteine and
    # acetaminophen mercapturate 151.16, 327.29, 231.23, 270.30 and
    # 312.24 g/mol, respectively"). Users dosing in mg should convert:
    # dose_umol = dose_mg / 0.15116.
    Cc        <- central        / vc
    Cc_gluc   <- central_gluc   / vc_gluc
    Cc_sulf   <- central_sulf   / vc_sulf
    Cc_cysmer <- central_cysmer / vc_cysmer

    Cc        ~ prop(propSd)
    Cc_gluc   ~ prop(propSd_gluc)
    Cc_sulf   ~ prop(propSd_sulf)
    Cc_cysmer ~ prop(propSd_cysmer)
  })
}
