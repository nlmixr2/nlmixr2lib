Smythe_2012_rifampicin <- function() {
  description <- "Semimechanistic population PK / enzyme-turnover autoinduction model for oral rifampicin in adult tuberculosis patients (Smythe 2012). One-compartment disposition with single-transit absorption (N = 1 FIX) feeds the central compartment; rifampicin plasma concentration drives a nonlinear Emax production-rate increase on a unitary-baseline enzyme pool, which in turn multiplies apparent oral clearance. CL/F and V/F are Anderson-Holford normal-fat-mass (NFM) allometrically scaled to a 70-kg patient with separate estimated Ffat contributions on CL/F and V/F. HIV-positive status increases V/F by 29.6%. IIV is on CL/F (correlated with V/F at 91.1%), V/F, and EC50; interoccasion variability is on MTT and bioavailability F. Residual error is combined additive + proportional. The same one-compartment + single-transit + autoinduction structural backbone (and the autoinduction parameters MTT, N, Emax, EC50, kENZ) is inherited verbatim by Clewe 2015 and Svensson 2016 -- see modellib('Clewe_2015_rifampicin') and modellib('Svensson_2016_rifampicin')."
  reference <- paste(
    "Smythe W, Khandelwal A, Merle C, Rustomjee R, Gninafon M, Lo MB,",
    "Sow OB, Olliaro PL, Lienhardt C, Horton J, Smith P, McIlleron H,",
    "Simonsson USH. (2012).",
    "A semimechanistic pharmacokinetic-enzyme turnover model for rifampin",
    "autoinduction in adult tuberculosis patients.",
    "Antimicrob Agents Chemother 56(4):2091-2098.",
    "doi:10.1128/AAC.05792-11.",
    sep = " "
  )
  vignette <- "Smythe_2012_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline (time-fixed in the Smythe 2012 analysis).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used together with FFM in the Anderson-Holford normal-fat-mass (NFM) scaling: NFM_param = FFM + Ffat_param * (WT - FFM); then CL/F is proportional to (NFM_CL / 70)^0.75 and V/F to (NFM_V / 70)^1.0. Cohort median 55 kg (interquartile range 50-61 kg across the four study sites; Smythe 2012 Table 1). Reference 70 kg for the allometric standardization. Range 38-80 kg (Methods 'Patients' paragraph 1).",
      source_name        = "WT"
    ),
    FFM = list(
      description        = "Fat-free mass derived from body weight, height, and sex via the Janmahasatian / Anderson-Holford WHSmax / WHS50 formula (Smythe 2012 Eqs. 8-10).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sex-specific FFM formula (Smythe 2012 Methods 'Covariate analysis' Eq. 8): men WHSmax = 42.92 kg/m^2 and WHS50 = 30.93 kg/m^2; women WHSmax = 37.99 kg/m^2 and WHS50 = 35.98 kg/m^2; FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT) with HT in meters and WT in kg. Used together with WT in the Anderson-Holford NFM scaling (paper Eqs. 9-10 give the NFM expressions for CL/F and V/F respectively). Cohort median FFM 45 kg (range 37-50 kg across sites; Smythe 2012 Table 1). The user data set should supply FFM directly as a covariate column; the vignette includes a chunk demonstrating the FFM derivation.",
      source_name        = "FFM"
    ),
    HIV_POS = list(
      description        = "HIV-1 antibody positive comorbidity indicator at study entry (1 = HIV-positive, 0 = HIV-negative). Time-fixed per subject; all HIV-positive subjects in the Smythe 2012 cohort were antiretroviral naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = "Multiplicative effect on V/F: V/F_HIV = V/F * (1 + 0.296 * HIV_POS); HIV-positive patients have a 29.6% larger apparent volume of distribution than HIV-negative reference subjects (Smythe 2012 Table 3 row 'V/F-HIV (%)' = 29.6, RSE 17.2%). Sex, age, study site, and creatinine clearance were tested but not retained in the final NFM-scaled covariate model (Methods 'Covariate analysis' paragraph 1; Results paragraphs 2-3). 54 of 174 patients were HIV-positive (49 South Africa, 0 Senegal, 2 Benin, 3 Guinea per Smythe 2012 Table 1).",
      source_name        = "HIV"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator: 1 = preinduced state (first dose, day 1), 2 = induced state (~day 29 after daily dosing, when autoinduction is at or near steady state). Smythe 2012 fit IOV across these two occasions on MTT and F.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Time-varying within subject; constant within an occasion. Decomposed inside model() into binary indicators oc1, oc2 that multiplex the per-occasion IOV etas on log(MTT) and log(F). Smythe 2012 reports a single shared IOV variance per parameter across the two occasions (Table 3 rows 'IOV MTT (%)' = 68 and 'IOV F (%)' = 16.2); encoded with the first occasion's eta estimated and the second occasion's eta fixed to the same variance (the Wilkins 2008 / Barnett 2018 OMEGA BLOCK(1) SAME idiom; nlmixr2 has no SAME shortcut). The two-occasion design is what the source data set supports; users wanting to simulate IOV at additional doses can extend the OCC range and add corresponding etaiov_*_<k> ~ fix() entries.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 174L,
    n_studies      = 1L,
    age_range      = "18-65 years (Methods 'Patients' paragraph 1)",
    age_median     = "28 years (South Africa) / 26 (Senegal) / 25 (Benin) / 34 (Guinea); pooled interquartile range 24-37 (South Africa, the largest cohort; Smythe 2012 Table 1)",
    weight_range   = "38-80 kg (Methods 'Patients' paragraph 1)",
    weight_median  = "55 kg (interquartile range 50-61 kg in South Africa, similar across sites; Smythe 2012 Table 1)",
    sex_female_pct = 34.5,
    race_ethnicity = "African; pooled across South Africa, Senegal, Benin, Guinea (Smythe 2012 Methods 'Patients' paragraph 1 and Table 1).",
    disease_state  = "Adults with newly diagnosed pulmonary tuberculosis enrolled in the control arm of the OFLOTUB phase III multicenter trial (NCT00216385). 54 of 174 patients (31%) were HIV-positive and antiretroviral naive (Smythe 2012 Table 1).",
    dose_range     = "Rifampicin 450 mg orally daily (patients < 50 kg body weight, n = 71) or 600 mg orally daily (patients >= 50 kg body weight, n = 103) co-administered with isoniazid, pyrazinamide, and ethambutol as fixed-dose combination tablets, six days per week, by directly observed therapy during the 2-month intensive phase of standard short-course antituberculosis treatment.",
    regions        = "South Africa (n = 101), Senegal (n = 27), Benin (n = 19), Guinea (n = 27).",
    notes          = "Baseline demographics from Smythe 2012 Table 1. 946 rifampicin plasma concentrations across the two occasions were modeled. Three sparse samples per occasion per patient: 1-2 h and 2.5-3.5 h post-dose for all patients, with the third sample block-randomized to 4-6 h or 8-10 h post-dose on occasion 1; on occasion 2 (median sampling day 29, range 26-50 days) the third sample was predose, 4-6 h, or 8-10 h post-dose. The clinical trial identifier is NCT00216385."
  )

  ini({
    # =========================================================================
    # Smythe 2012 Table 3 final model = Model 3 (NFM allometric scaling) with
    # HIV included as a 29.6% increase on V/F. Time unit is HOURS throughout
    # (matching the paper's parameter table). CL/F and V/F are reported in
    # L/h and L respectively, standardized to a 70-kg patient via NFM.
    # =========================================================================

    # --- Structural disposition (preinduced, allometrically standardized) ---
    lcl       <- log(10.0)              ; label("Apparent oral clearance CL/F (L/h) at the preinduced state, NFM-standardized to a 70-kg patient") # Smythe 2012 Table 3: TV(CL/F)STD = 10.0 L/h (RSE 3.7%)
    lvc       <- log(86.7)              ; label("Apparent central volume of distribution V/F (L), NFM-standardized to a 70-kg patient")             # Smythe 2012 Table 3: TV(V/F)STD = 86.7 L (RSE 2.3%)

    # --- Anderson-Holford NFM (normal fat mass) fractional fat contributions ---
    # NFM_CL = FFM + e_fat_cl * (WT - FFM); allometric exponent 0.75 on CL/F.
    # NFM_VC = FFM + e_fat_vc * (WT - FFM); allometric exponent 1.0  on V/F.
    e_fat_cl  <- 0.311                  ; label("Fractional fat contribution to NFM for CL/F (Ffat_CL/F, unitless)") # Smythe 2012 Table 3: (Ffat)CL/F = 0.311 (RSE 40.2%)
    e_fat_vc  <- 0.188                  ; label("Fractional fat contribution to NFM for V/F  (Ffat_V/F,  unitless)") # Smythe 2012 Table 3: (Ffat)V/F  = 0.188 (RSE 49.1%)

    # --- HIV covariate effect on V/F: multiplicative (1 + e_hiv_vc * HIV_POS) ---
    e_hiv_pos_vc <- 0.296                ; label("Fractional increase in V/F for HIV-positive patients (unitless)")  # Smythe 2012 Table 3: V/F-HIV = 29.6% (RSE 17.2%)

    # --- Transit absorption (Savic / Wilkins parameterization) ---
    lmtt      <- log(0.713)             ; label("Mean transit time MTT (h)")                                          # Smythe 2012 Table 3: MTT = 0.713 h (RSE 1.6%)
    nn_fix    <- fixed(1)               ; label("Number of transit compartments N (unitless, integer)")               # Smythe 2012 Table 3: 'No. of transit compartments = 1 FIX'

    # --- Autoinduction enzyme-pool turnover (Hassan-type, unity baseline) ---
    # d(enz_pool)/dt = kENZ * (1 + EFF) - kENZ * enz_pool, with the baseline
    # production rate of kENZ chosen so that enz_pool = 1 at steady state when
    # rifampicin concentration is zero (Smythe 2012 Methods 'Structural model'
    # paragraph 5).
    lemax     <- log(1.04)              ; label("Maximal fractional increase in enzyme production rate Emax (unitless)") # Smythe 2012 Table 3: Emax = 1.04 (RSE 2.6%)
    lec50     <- log(0.0705)            ; label("Rifampicin plasma concentration producing half Emax (mg/L)")            # Smythe 2012 Table 3: EC50 = 0.0705 mg/L (RSE 6.3%)
    lkenz     <- log(0.00369)           ; label("First-order degradation rate constant of the enzyme pool kENZ (1/h)")   # Smythe 2012 Table 3: kENZ = 0.00369/h (RSE 5.6%) -> turnover half-life log(2)/kENZ = 187.8 h ~ 7.8 days

    # --- Bioavailability anchor (CL/F and V/F are apparent F-relative) ---
    lfdepot   <- fixed(log(1))          ; label("Oral bioavailability F (fixed at 1 because CL and V are apparent F-relative)") # Implicit anchor: Smythe 2012 reports CL/F and V/F (no separate F estimate); see Methods 'Stochastic model' paragraph 1

    # =========================================================================
    # IIV. Smythe 2012 Table 3 reports IIV (%CV) on CL/F, V/F, and EC50, with
    # a 91.1% correlation between CL/F and V/F. omega^2 = log(1 + CV^2) for a
    # log-normal random effect.
    # =========================================================================
    # var_CL = log(1 + 0.30^2)  = 0.08618
    # var_VC = log(1 + 0.192^2) = 0.03625
    # cov(CL,VC) = 0.911 * sqrt(var_CL * var_VC)
    #            = 0.911 * sqrt(0.08618 * 0.03625)
    #            = 0.911 * sqrt(0.003124) = 0.911 * 0.05589 = 0.05091
    etalcl + etalvc ~ c(0.08618,
                        0.05091, 0.03625)   # Smythe 2012 Table 3: IIV CL/F = 30.0% CV (RSE 12.3%), IIV V/F = 19.2% CV (RSE 14.8%), CL-V correlation = 91.1% (RSE 20.7%)

    # var_EC50 = log(1 + 4.93^2) = log(25.305) = 3.231
    etalec50 ~ 3.231                       # Smythe 2012 Table 3: IIV EC50 = 493% CV (RSE 19%)

    # =========================================================================
    # IOV across the two occasions (preinduced and induced states). Smythe 2012
    # reports a single shared variance per parameter, encoded here with the
    # first occasion's eta estimated and the second occasion fixed to the
    # same variance (mirrors the Barnett_2018_rifampicin / Wilkins_2008
    # OMEGA BLOCK(1) SAME idiom -- nlmixr2 has no SAME shortcut).
    # =========================================================================
    # var_iov_MTT = log(1 + 0.68^2)  = log(1.4624) = 0.38014
    etaiov_mtt_1 ~ 0.38014                 # Smythe 2012 Table 3: IOV MTT = 68% CV (RSE 7%)
    etaiov_mtt_2 ~ fix(0.38014)            # same IOV variance across occasions (Smythe 2012 BLOCK SAME idiom)
    # var_iov_F   = log(1 + 0.162^2) = log(1.02624) = 0.02591
    etaiov_f_1   ~ 0.02591                 # Smythe 2012 Table 3: IOV F = 16.2% CV (RSE 11.2%)
    etaiov_f_2   ~ fix(0.02591)            # same IOV variance across occasions

    # =========================================================================
    # Residual error -- combined additive (mg/L) + proportional (fraction).
    # =========================================================================
    addSd     <- 0.965                    ; label("Additive residual error (mg/L)")        # Smythe 2012 Table 3: Additive error = 0.965 mg/L (RSE 2.8%)
    propSd    <- 0.099                    ; label("Proportional residual error (fraction)") # Smythe 2012 Table 3: Proportional error = 9.9% (RSE 4.7%)
  })

  model({
    # --- 1. Decompose OCC into binary indicators for IOV multiplexing on -----
    #        log(MTT) and log(F). Occasion 1 = preinduced (first dose);
    #        occasion 2 = induced state (~day 29). Users supplying additional
    #        occasions can extend etaiov_*_<k> ~ fix(.) entries in ini().
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_f   <- oc1 * etaiov_f_1   + oc2 * etaiov_f_2

    # --- 2. Anderson-Holford NFM allometric scaling (Smythe 2012 Eqs. 9-10) --
    #        NFM_CL = FFM + Ffat_CL * (WT - FFM)
    #        NFM_VC = FFM + Ffat_VC * (WT - FFM)
    #        CL/F = TV(CL/F)_STD * (NFM_CL / 70)^0.75
    #        V/F  = TV(V/F)_STD  * (NFM_VC / 70)^1.0 * (1 + e_hiv_pos_vc * HIV_POS)
    nfm_cl <- FFM + e_fat_cl * (WT - FFM)
    nfm_vc <- FFM + e_fat_vc * (WT - FFM)
    cl_base <- exp(lcl + etalcl) * (nfm_cl / 70) ^ 0.75
    vc      <- exp(lvc + etalvc) * (nfm_vc / 70) ^ 1.0 * (1 + e_hiv_pos_vc * HIV_POS)

    # --- 3. Autoinduction enzyme-pool individual parameters and transit chain.
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    kenz <- exp(lkenz)
    mtt  <- exp(lmtt + iov_mtt)
    nn   <- nn_fix
    ktr  <- (nn + 1) / mtt

    # --- 4. Plasma concentration (mg/L) and effective time-varying clearance.
    #        The enz_pool state starts at 1 (Smythe 2012 normalization: the
    #        zero-order enzyme production rate is set to kENZ so that the
    #        pool is in steady-state equal to unity when rifampicin is
    #        absent). When rifampicin enters the central compartment, EFF
    #        > 0 raises the production rate above kENZ; the pool grows above
    #        1, and the effective clearance cl_base * enz_pool exceeds the
    #        preinduced cl_base by the same factor.
    Cc <- central / vc
    cl <- cl_base * enz_pool

    # --- 5. PK ODE system: depot -> transit1 -> central with concentration- ---
    #        driven autoinduction of the enzyme pool. With N = 1 transit
    #        compartment the Savic chain reduces to a single transit-rate
    #        constant ktr = (N + 1) / MTT = 2 / MTT into both depot loss and
    #        transit1 inflow. The elimination from central is written as
    #        cl * Cc (= cl * central / vc = kel * central) so the
    #        autoinduction multiplier is read directly off cl.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1
    d/dt(central)  <-  ktr * transit1 - cl * Cc
    d/dt(enz_pool) <-  kenz * (1 + emax * Cc / (ec50 + Cc)) - kenz * enz_pool

    # Enzyme pool starts at unity by construction.
    enz_pool(0) <- 1.0

    # --- 6. Bioavailability. F is anchored at 1 (CL and V are apparent      ---
    #        F-relative) and carries log-scale IOV across the two occasions
    #        with variance log(1 + 0.162^2).
    f(depot) <- exp(lfdepot + iov_f)

    # --- 7. Observation: rifampicin plasma concentration with combined     ---
    #        additive + proportional residual error.
    Cc ~ add(addSd) + prop(propSd)
  })
}
