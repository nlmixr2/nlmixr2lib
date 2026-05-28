Svensson_2016_rifampicin <- function() {
  description <- "Combined population PK/PD model for rifampicin in adults with drug-susceptible pulmonary tuberculosis: a one-compartment, single-transit, oral PK model with first-order plasma-concentration-driven autoinduction of clearance via an enzyme-pool turnover (structure from Smythe 2012) linked to the Multistate Tuberculosis Pharmacometric (MTP) three-state bacterial disease model (fast-, slow-, and nonmultiplying Mycobacterium tuberculosis states; structure from Clewe 2016) with rifampicin drug effects as fixed-at-100% on/off inhibition of fast-multiplying bacterial growth plus second-order plasma-concentration-driven death of slow- and nonmultiplying bacteria; all PK parameters and all MTP transfer/growth rates are fixed to the upstream-paper estimates, while the system carrying capacity Bmax (with 152% CV IIV) and the two second-order death rates SDk and NDk are re-estimated against 19 patients from a 1966-1977 Kenyan rifampicin monotherapy trial."
  reference <- paste(
    "Svensson R. J., Simonsson U. S. H. (2016).",
    "Application of the Multistate Tuberculosis Pharmacometric Model in",
    "Patients With Rifampicin-Treated Pulmonary Tuberculosis.",
    "CPT: Pharmacometrics & Systems Pharmacology 5(5):264-273.",
    "doi:10.1002/psp4.12079.",
    "PK structure (one-compartment + single transit + enzyme-pool",
    "autoinduction + Anderson-Holford NFM allometric scaling) adapted",
    "from Smythe et al. (2012) Antimicrob Agents Chemother 56(4):2091-2098",
    "doi:10.1128/AAC.05792-11.",
    "MTP disease model structure (three bacterial substates with",
    "time-dependent fast-to-slow transfer) from Clewe et al. (2016)",
    "J Antimicrob Chemother 71(4):964-974 doi:10.1093/jac/dkv478.",
    sep = " "
  )
  vignette <- "Svensson_2016_rifampicin"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight. Time-fixed per subject in the Svensson 2016 analysis (set to the Smythe 2012 cohort mean of 56 kg for every patient because no individual covariate values were available for the 1966-1977 sputum dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used together with FFM in the Anderson-Holford normal-fat-mass (NFM) allometric scaling of CL/F and V/F (Smythe 2012 Table 3 model 3): NFM_param = FFM + Ffat_param * (WT - FFM), then CL/F propto NFM_CL^0.75 and V/F propto NFM_V^1.0, both standardized to a 70-kg patient. Reference value 56 kg. Svensson 2016 Methods 'Population pharmacokinetic model' paragraph 2 reports the cohort body weight assumption.",
      source_name        = "WT"
    ),
    FFM = list(
      description        = "Fat-free mass. Time-fixed per subject in the Svensson 2016 analysis (set to the Smythe 2012 cohort mean of 45 kg for every patient because no individual covariate values were available for the 1966-1977 sputum dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed from body weight, height, and sex via the Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065) formula in the Smythe 2012 cohort; for this PD model, FFM is supplied directly as a covariate column (default 45 kg). Drives the Anderson-Holford NFM allometric scaling alongside WT. Reference value 45 kg. Svensson 2016 Methods 'Population pharmacokinetic model' paragraph 2 reports the cohort FFM assumption.",
      source_name        = "FFM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19L,
    n_studies      = 1L,
    age_range      = "adult (Jindani et al. 1980 source trial)",
    weight_range   = "set to the cohort-mean 56 kg for all patients (no individual demographics available)",
    weight_median  = "56 kg (cohort-mean, fixed for all patients)",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in the source paper (Kenyan patients)",
    disease_state  = "Treatment-naive drug-susceptible pulmonary tuberculosis (Mycobacterium tuberculosis); all patients assumed HIV-negative and in stationary phase of infection at trial entry.",
    dose_range     = "Rifampicin oral monotherapy at 5 mg/kg (n=3), 10 mg/kg (n=8), or 20 mg/kg (n=8) once daily for 14 days at 08:00, plus a no-treatment negative control arm (n=4) for the structural disease model fit (the negative-control 4 patients are not counted in n_subjects=19).",
    regions        = "Kenya (1966-1977 trial reanalysed by Svensson 2016).",
    notes          = "Baseline demographics from Svensson 2016 Methods 'Patients and study design'. Sputum CFU sampled every 2 days during the 12-hour overnight window 8 PM-8 AM, including two baseline measurements before the first rifampicin dose. The original trial pre-dated formal IRB review; data reuse was approved by the UK National Research Ethics Service via the PreDiCT-TB consortium. The PD model also performs external validation against three retrospective rifampicin-monotherapy trials (Sirgel 2005 n=14-16 per arm at 300 / 600 mg / 20 mg/kg in South Africa; Diacon 2007 n=13 at 20 mg/kg in South Africa; Rustomjee 2008 n=15 at 600 mg in South Africa; Svensson 2016 Table 1)."
  )

  ini({
    # =========================================================================
    # Rifampicin PK-enzyme turnover model (structure: Smythe 2012 model 3,
    # parameters: Svensson 2016 Table 2 'Drug PK parameters' all FIX).
    #
    # Time unit in this model is DAYS (chosen to match the MTP model rates).
    # Smythe 2012 reports CL/F in L/h and kENZ in 1/h; multiplying by 24
    # converts to L/day and 1/day, and dividing MTT by 24 converts h to days.
    # All ten PK parameters are fixed exactly as Svensson 2016 reports them
    # (no covariates, no IIV, no residual error on the PK; this is the
    # 'population PK parameter approach' of Zhang et al. 2003 used because
    # the Svensson 2016 sputum dataset contains no rifampicin concentration
    # observations).
    # =========================================================================
    lcl     <- fixed(log(10.0 * 24))     ; label("Apparent oral clearance CL/F (L/day) at the preinduced state, standardized to 70 kg")  # Svensson 2016 Table 2 CL/F = 10.0 L/h * 24 = 240 L/day; carried from Smythe 2012 Table 3
    lvc     <- fixed(log(86.7))          ; label("Apparent central volume of distribution V/F (L), standardized to 70 kg")               # Svensson 2016 Table 2 V/F = 86.7 L; carried from Smythe 2012 Table 3
    lmtt    <- fixed(log(0.713 / 24))    ; label("Mean transit time (day, on log scale)")                                                  # Svensson 2016 Table 2 MTT = 0.713 h / 24 = 0.0297 day; carried from Smythe 2012 Table 3
    nn_fix  <- fixed(1)                  ; label("Number of transit compartments (integer, unitless)")                                    # Svensson 2016 Table 2 NN = 1.00 FIX; carried from Smythe 2012 Table 3 'No. of transit compartments = 1 FIX'
    lemax   <- fixed(log(1.04))          ; label("Maximal fractional increase in the enzyme production rate (unitless)")                  # Svensson 2016 Table 2 Emax = 1.04; carried from Smythe 2012 Table 3
    lec50   <- fixed(log(0.0705))        ; label("Rifampicin plasma concentration producing half Emax (mg/L)")                            # Svensson 2016 Table 2 EC50 = 0.0705 mg/L; carried from Smythe 2012 Table 3
    lkenz   <- fixed(log(0.00369 * 24))  ; label("Rate constant for first-order degradation of the enzyme pool (1/day)")                  # Svensson 2016 Table 2 kENZ = 0.00369/h * 24 = 0.0886/day; carried from Smythe 2012 Table 3
    lfdepot <- fixed(log(1))             ; label("Oral bioavailability (fixed at 1 because CL and V are reported apparent F-relative)")   # Svensson 2016 Table 2 F = 1.00 FIX
    e_fat_cl <- fixed(0.311)             ; label("Ffat on CL/F: fractional contribution of fat mass to NFM in the CL allometric scaling") # Svensson 2016 Table 2 (Ffat)CL/F = 0.311 FIX; carried from Smythe 2012 Table 3 model 3
    e_fat_vc <- fixed(0.188)             ; label("Ffat on V/F:  fractional contribution of fat mass to NFM in the V  allometric scaling") # Svensson 2016 Table 2 (Ffat)V/F  = 0.188 FIX; carried from Smythe 2012 Table 3 model 3

    # =========================================================================
    # MTP disease model (Clewe 2016 / Svensson 2016 Table 2 'Multistate
    # Tuberculosis Pharmacometric model'). All transfer / growth rates and
    # initial bacterial numbers F0 and S0 are fixed to the in vitro
    # estimates of Clewe 2016 because they were unidentifiable from the
    # short-duration clinical dataset (Svensson 2016 Results paragraph 4:
    # 'Parameters kG, kFSlin, kSF, kSN, kNS, and kFN were unidentifiable.
    # This was probably because of the short study duration and that data
    # were obtained at a stationary phase of infection.'). Only Bmax was
    # re-estimated. Time unit is DAYS throughout.
    # =========================================================================
    lkg     <- fixed(log(0.206))        ; label("Fast-multiplying bacterial growth rate kG (1/day)")                                    # Svensson 2016 Table 2 kG = 0.206/day FIX; carried from Clewe 2016
    lkfn    <- fixed(log(8.97e-7))      ; label("Transfer rate from fast- to nonmultiplying state kFN (1/day)")                          # Svensson 2016 Table 2 kFN = 8.97e-7/day FIX
    lksn    <- fixed(log(0.186))        ; label("Transfer rate from slow- to nonmultiplying state kSN (1/day)")                          # Svensson 2016 Table 2 kSN = 0.186/day FIX
    lksf    <- fixed(log(0.0145))       ; label("Transfer rate from slow- to fast-multiplying state kSF (1/day)")                        # Svensson 2016 Table 2 kSF = 0.0145/day FIX
    lkns    <- fixed(log(0.00123))      ; label("Transfer rate from non- to slow-multiplying state kNS (1/day)")                         # Svensson 2016 Table 2 kNS = 0.00123/day FIX
    lkfslin <- fixed(log(0.00166))      ; label("Time-dependent slope of the fast-to-slow transfer rate kFSlin (1/day^2)")               # Svensson 2016 Table 2 kFSlin = 0.00166/day^2 FIX
    lf0     <- fixed(log(4.10))         ; label("Initial bacterial number of fast-multiplying state at infection F0 (1/mL)")             # Svensson 2016 Table 2 F0 = 4.10/mL FIX
    lrbase     <- fixed(log(9770))         ; label("Initial bacterial number of slow-multiplying state at infection S0 (1/mL)")             # Svensson 2016 Table 2 S0 = 9770/mL FIX
    lbmax   <- log(2.61e9)              ; label("System carrying capacity Bmax (1/mL); estimated typical value")                         # Svensson 2016 Table 2 Bmax = 2.61e9/mL (RSE 30.5%, 95% CI 1.51e9-4.52e9); only MTP parameter re-estimated in the clinical fit

    # =========================================================================
    # Rifampicin drug effects on the bacterial states (Svensson 2016
    # Table 2 'Exposure-response parameters'). All three effects act on
    # plasma concentration Crif = central / vc (in mg/L).
    # =========================================================================
    fg_on_off <- fixed(1.00)            ; label("Fractional inhibition of fast-multiplying growth when rifampicin Crif > 0 (unitless, on/off)") # Svensson 2016 Table 2 FGon/off = 1.00 FIX (estimated close to 1 and then fixed)
    lsdk      <- log(0.200)             ; label("Second-order slow-multiplying death rate SDk (L/mg/day); estimated typical value")              # Svensson 2016 Table 2 SDk = 0.200 L/mg/day (RSE 41.6%, 95% CI 0.0854-0.390)
    lndk      <- log(0.106)             ; label("Second-order nonmultiplying death rate NDk (L/mg/day); estimated typical value")                # Svensson 2016 Table 2 NDk = 0.106 L/mg/day (RSE 19.0%, 95% CI 0.0643-0.188)

    # =========================================================================
    # IIV (Svensson 2016 Table 2 IIV row). Only Bmax carries IIV; the other
    # MTP parameters were fixed to in vitro estimates and the PK has no IIV
    # in the population-PK-parameter approach.
    # =========================================================================
    etalbmax  ~ 1.197                   # Svensson 2016 Table 2 IIV Bmax = 152% CV; omega^2 = log(1 + 1.52^2) = 1.197

    # =========================================================================
    # Residual error. Svensson 2016 Table 2 reports TWO additive-on-log-scale
    # components: an across-replicate component (110% CV) and a
    # between-replicate-within-same-sputum-sample component (23.1% CV),
    # implemented via the NONMEM L2 data item to absorb the correlation
    # between technical replicates of the same sputum.  nlmixr2lib has no
    # idiomatic encoding for L2-style nested residuals, so the two
    # components are folded into a single additive (on log-CFU/mL) SD via
    # sigma_total = sqrt(sigma_repl^2 + sigma_betw^2); see vignette
    # Assumptions and deviations.
    # =========================================================================
    addSd <- 1.124                       ; label("Combined additive residual SD on log(CFU/mL); folds Svensson 2016 e=1.10 + e_repl=0.231 into one component") # Svensson 2016 Table 2 e = 110% + e_repl = 23.1% combined as sqrt(1.10^2 + 0.231^2) = 1.124
  })

  model({
    # --- 1. NFM allometric scaling (Anderson and Holford 2008, Smythe 2012 -----
    #        model 3). NFM contains an estimated fractional contribution of
    #        fat mass (Ffat) on each disposition parameter, then the standard
    #        allometric exponents 0.75 on clearance and 1.0 on volume.
    nfm_cl <- FFM + e_fat_cl * (WT - FFM)
    nfm_vc <- FFM + e_fat_vc * (WT - FFM)

    # --- 2. Individual PK parameters (no IIV; PK fixed in the
    #        population-PK-parameter approach). cl_base is the preinduced
    #        clearance; the effective clearance is cl_base * enz_pool, with
    #        the enzyme-pool state starting at 1 by construction.
    cl_base <- exp(lcl) * (nfm_cl / 70) ^ 0.75
    vc      <- exp(lvc) * (nfm_vc / 70) ^ 1.0
    mtt     <- exp(lmtt)
    kenz    <- exp(lkenz)
    emax    <- exp(lemax)
    ec50    <- exp(lec50)
    nn      <- nn_fix
    ktr     <- (nn + 1) / mtt
    cl      <- cl_base * enz_pool

    # --- 3. Disease-model individual parameters.
    kg     <- exp(lkg)
    kfn    <- exp(lkfn)
    ksn    <- exp(lksn)
    ksf    <- exp(lksf)
    kns    <- exp(lkns)
    kfslin <- exp(lkfslin)
    bmax   <- exp(lbmax + etalbmax)
    sdk    <- exp(lsdk)
    ndk    <- exp(lndk)

    # --- 4. Time-dependent fast-to-slow transfer kFS(t) = kFSlin * t. The
    #        variable `t` is rxode2's elapsed simulation time in the model's
    #        time unit (here, days since infection). Svensson 2016 assumes
    #        t = 0 is infection and treatment begins at t = 150 days, the
    #        time point at which the in vitro MTP model predicts stationary
    #        phase (Svensson 2016 Results paragraph 2).
    kfs <- kfslin * t

    # --- 5. Rifampicin plasma concentration (mg/L) and the on/off drug
    #        indicator. Once dosing has begun, Crif stays above zero
    #        between doses, so FGon/off = 1 effectively shuts off
    #        fast-multiplying growth for the entire treatment window
    #        (Svensson 2016 Discussion paragraph 7: 'As the model included
    #        no concentration-dependent effect on fast-multiplying bacteria,
    #        but only complete inhibition of the growth of fast-multiplying
    #        bacteria, no circadian change is seen for this state.').
    Crif <- central / vc
    ind_drug <- (Crif > 0)
    fg_inhib <- 1 - fg_on_off * ind_drug

    # --- 6. Total live bacteria (used in the Gompertz growth term). The
    #        small additive epsilon avoids log(0) when all states have been
    #        eradicated; it is negligible at biological scales (CFU/mL).
    total_bact <- fast + slow + nonm + 1e-6

    # --- 7. PK ODE system (single transit absorption + central + ENZ pool).
    d/dt(depot)     <- -ktr * depot
    d/dt(transit1)  <-  ktr * depot - ktr * transit1
    d/dt(central)   <-  ktr * transit1 - cl_base * enz_pool * Crif
    d/dt(enz_pool)  <-  kenz * (1 + emax * Crif / (ec50 + Crif)) - kenz * enz_pool

    # Enzyme pool starts at unity by construction (Smythe 2012 Eq 1
    # 'To normalize the enzyme concentrations to unity at baseline, the
    # zero-order production rate of the enzyme was set to kENZ.').
    enz_pool(0) <- 1.0

    # --- 8. MTP disease ODE system.
    d/dt(fast) <- kg * fast * log(bmax / total_bact) * fg_inhib -
                  kfs * fast - kfn * fast + ksf * slow
    d/dt(slow) <- kfs * fast - ksf * slow - ksn * slow + kns * nonm -
                  sdk * Crif * slow
    d/dt(nonm) <- kfn * fast + ksn * slow - kns * nonm -
                  ndk * Crif * nonm

    # Initial bacterial loads at infection (t = 0). N0 = 0 is the standard
    # MTP initial condition (no nonmultiplying bacteria at the moment of
    # infection; the nonmultiplying pool fills up over time via the kFN
    # and kSN transfers).
    fast(0) <- exp(lf0)
    slow(0) <- exp(lrbase)
    nonm(0) <- 0.0

    # --- 9. Bioavailability.
    f(depot) <- exp(lfdepot)

    # --- 10. Observation: natural log of (fast + slow)/mL of sputum, modelled
    #         as additive on the log scale (LTBS in Svensson 2016 NONMEM
    #         implementation). Only the fast- and slow-multiplying states
    #         contribute to CFU; the nonmultiplying state is nonculturable
    #         by definition (Svensson 2016 Methods paragraph 5: 'Only the
    #         bacterial numbers in the fast- and slow-multiplying states
    #         were visible as CFU, the nonmultiplying state was considered
    #         nonculturable.'). The published model uses an additional
    #         sputum-sample compartment that averages (fast + slow) over
    #         the 12-h overnight collection window; nlmixr2lib uses the
    #         instantaneous mid-/last-time-point readout instead because
    #         Svensson 2016 Results paragraph 6 reports that 'The last
    #         time point and midtime point methods provided similar OFV'
    #         to the sample-compartment method. See vignette Assumptions
    #         and deviations.
    Sputum_lnCFU <- log(fast + slow + 1e-6)
    Sputum_lnCFU ~ add(addSd)
  })
}
