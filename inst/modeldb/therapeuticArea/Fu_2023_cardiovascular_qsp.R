Fu_2023_cardiovascular_qsp <- function() {
  description <- "QSP (rat). Cardiovascular systems (CVS) hemodynamic model developed by Snelder et al. (2013, 2014) and re-used as a fixed system by Fu et al. (2023) for a stochastic-simulation and re-estimation (SSE) identifiability analysis. The model links five interdependent hemodynamic variables in rats -- heart rate (HR, beats/min), stroke volume (SV, mL/beat), total peripheral resistance (TPR, mmHg*min/mL), cardiac output (CO = HR*SV, mL/min) and mean arterial pressure (MAP = CO*TPR, mmHg) -- through three coupled indirect-response ODEs on HR, SV-turnover (SVT), and TPR, a shared MAP-negative-feedback term on the production rate of each state, cosine circadian rhythms on HR and TPR, and a paper-specific direct HR-on-SV coupling (SV = SVT*(1 - HRSV*log(HR/BSLHR))). Drug PK is a hypothetical one-compartment IV-bolus model with k = 0.173/h (half-life 4 h). In the on-disk NONMEM control stream, the drug effect (Emax model on the drug amount) acts on HR-production, matching the mode-of-action-on-HR scenario of Fu 2023 Section 2.3.3; the vignette shows how to relocate the effect to SV or TPR for the SSE identifiability sweep. All 12 system-specific parameters (baselines, k_out rate constants, MAP-feedback constant, HR-on-SV coupling constant, four circadian rhythm parameters) and all 3 IIV variances (on the three baselines) are FIXED to the Snelder-derived rat estimates transcribed in Fu 2023 Supplemental Materials S1."

  reference <- paste(
    "Fu Y, Snelder N, Guo T, van der Graaf PH, van Hasselt JGC.",
    "Evaluation of a Cardiovascular Systems Model for Design and Analysis",
    "of Hemodynamic Safety Studies.",
    "Pharmaceutics. 2023 Apr 7;15(4):1175.",
    "doi:10.3390/pharmaceutics15041175.",
    "System-specific structure and parameter values inherited from",
    "Snelder et al. (2013, 2014); Fu 2023 fixes the system parameters",
    "as reported in Supplemental Materials S1 (NONMEM control stream).",
    sep = " "
  )
  vignette <- "Fu_2023_cardiovascular"

  # Paper-mechanistic ODE states of the CVS system (Fu 2023 Suppl. S1 $MODEL:
  # COMP(HR), COMP(SVT), COMP(TPR)). Declared here so checkModelConventions()
  # accepts them as legitimate paper-specific compartments without requiring
  # generic canonicalisation -- see compartment-names.md section "Paper-specific
  # compartments" and the Bizzotto 2016 / Yuan 2019 / FehlingKaschek 2019
  # precedents. If additional cardiovascular-systems models are added later
  # (Snelder 2013, Snelder 2014, van Rijn 2011 dog CVS, etc.) these three
  # hemodynamic states may be promoted to first-class canonicals.
  paper_specific_compartments <- c("hr", "svt", "tpr")

  units <- list(
    time          = "hour",
    dosing        = "amount unit (implicit V = 1 in the hypothetical-drug PK; the drug-effect Emax uses the central amount directly so any amount unit is valid provided EC50 is on the same scale)",
    concentration = "beats/min (HR), mL/beat (SV), mmHg*min/mL (TPR), mL/min (CO), mmHg (MAP), and hypothetical amount units for the drug (Cc)"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Wistar-Kyoto, per Snelder 2013 J Pharmacokinet Pharmacodyn 40(3):357-373 and 2014 Br J Pharmacol 171(22):5076-5092, from which Fu 2023 inherits the CVS system parameters)",
    n_subjects     = "3, 4, or 5 per SSE scenario (Fu 2023 Section 2.3.5; 100 simulated datasets per scenario)",
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy rat (baseline hemodynamics with pre-clinical iv-bolus drug perturbation; hypothetical drug with 4-h half-life applied at 0.1, 0.3, 1, 3, or 10 mg/kg iv bolus across the SSE scenarios of Fu 2023 Section 2.2)",
    dose_range     = "0.1, 1, or 10 mg/kg iv bolus on days 1-3 (three ascending doses) or 0.1, 0.3, 1, 3, 10 mg/kg on days 1-5 (five ascending doses); the hypothetical drug is a placeholder for pre-clinical safety pharmacology compounds in general",
    regions        = "Rat (Leiden Academic Centre for Drug Research and LAP&P Consultants BV, Netherlands; Snelder-cohort data)",
    n_ode_states   = 4L,
    n_fixed_system_params = 12L,
    n_fixed_iiv    = 3L,
    n_drug_params  = 2L,
    fixed_params_source = "All 12 system-specific parameters and 3 IIV variances are fixed to the Snelder-derived rat estimates transcribed in Fu 2023 Supplemental Materials S1 (NONMEM control stream, $THETA and $OMEGA blocks with the corresponding FIX flags). Fu 2023 states (Section 2.2): 'For all the scenarios focusing on identification of the MoA, drug-specific parameters (Emax and EC50), interindividual variability (IIV), and residual errors of each type of measurement were estimated; the system-specific parameters were fixed to the parameter values from the published CVS model [7].'",
    notes          = "Fu 2023 is a simulation-and-re-estimation identifiability study, not a data-fit paper. The rat physiologic backbone (baselines, k_out rate constants, MAP feedback, circadian rhythms, HR-on-SV coupling) is inherited from Snelder et al. and held constant across all SSE scenarios. The drug-specific parameters (Emax, EC50) are the ones the paper evaluates for identifiability."
  )

  ini({
    # ============================================================
    # Structural parameters
    #
    # All values are transcribed from Fu 2023 Supplemental Materials
    # S1 (NONMEM control stream, $THETA block). Every system-specific
    # parameter (TH1, TH4-TH15) carries a FIX flag in the source
    # control stream; the drug-specific parameters (TH2 EC50, TH3
    # EMAX) are the identifiability targets of the paper and are
    # left as initial-estimate anchors here.
    # ============================================================

    # ---- Hypothetical drug PK -------------------------------------
    # First-order elimination rate constant of the one-compartment
    # iv-bolus PK model (Fu 2023 Section 2.1 and Supplemental S1
    # $THETA TH1 K = 0.17325 FIX; half-life = ln(2)/k = 4.00 h,
    # matching the 4-h realistic-half-life claim in the same section).
    lkel <- fixed(log(0.17325))  ; label("Hypothetical-drug elimination rate constant, k (1/h)")  # Fu 2023 Suppl. S1 $THETA TH1 (0.17325 FIX)

    # ---- Drug-specific PD (identifiability targets) --------------
    # EC50 on the same scale as the central-compartment amount
    # (Supplemental S1 uses A(1) directly in the Emax term, so EC50
    # is expressed in the same "amount = concentration" scale used by
    # the hypothetical-drug PK). Fu 2023 SSE scenarios test EC50 at
    # 100, 1000, or 100000 ng/mL and Emax at -1, 1, or 10 (Table S1).
    # The 100-ng/mL initial estimate here matches the "quantification
    # of system- and drug-specific parameters" scenario at the top of
    # Table S1.
    lec50 <- log(100)            ; label("Half-maximal drug concentration for HR-production effect, EC50 (concentration units, same scale as A(1))")  # Fu 2023 Suppl. S1 $THETA TH2 (initial 100); Table S1 test values 100, 1000, 100000
    emax  <- -1                  ; label("Maximum fractional drug effect on HR production, Emax (signed; -1 = inhibition of HR-production rate, +1 or +10 = stimulation)")  # Fu 2023 Suppl. S1 $THETA TH3 (initial 1 with bounds -1 to 1); Table S1 test values -1, 1, 10; Emax kept on linear scale because it can be signed

    # ---- MAP negative-feedback constant --------------------------
    # Log-transformed for consistency with the positive-constrained
    # parameter convention. Represents the strength of MAP feedback
    # on all three primary hemodynamic states (HR, SV, TPR): the
    # production-rate multiplier is (1 - fb * MAP), so higher fb
    # means stronger MAP-driven suppression of production.
    lfb <- fixed(log(0.0029))    ; label("MAP negative-feedback constant, fb (1/mmHg)")  # Fu 2023 Suppl. S1 $THETA TH4 (0.0029 FIX)

    # ---- Rat hemodynamic baselines (Snelder-derived, fixed) ------
    # Log-transformed so the log-normal IIV expression matches the
    # supplement's THETA(i)*EXP(ETA(j)) form (Fu 2023 Suppl. S1 $PK:
    # BSLHR = THETA(5)*EXP(ETA(1)), BSLMAP = THETA(6)*EXP(ETA(2)),
    # BSLCO = THETA(7)*EXP(ETA(3)); BSLSV and BSLTPR are derived
    # from these three via BSLSV = BSLCO/BSLHR and BSLTPR =
    # BSLMAP/BSLCO).
    lrbase_HR  <- fixed(log(310)) ; label("Baseline heart rate, BSLHR (beats/min)")                     # Fu 2023 Suppl. S1 $THETA TH5 (310 FIX)
    lrbase_MAP <- fixed(log(155)) ; label("Baseline mean arterial pressure, BSLMAP (mmHg)")             # Fu 2023 Suppl. S1 $THETA TH6 (155 FIX)
    lrbase_CO  <- fixed(log(69))  ; label("Baseline cardiac output, BSLCO (mL/min)")                    # Fu 2023 Suppl. S1 $THETA TH7 (69 FIX)

    # ---- k_out dissipation rate constants (Snelder-derived, fixed)
    # Each of HR, SV, and TPR is an indirect-response state with its
    # own first-order dissipation rate constant. The k_in production
    # rate is derived at steady state as k_in_X = k_out_X * BSL_X /
    # (1 - fb * BSLMAP), so no k_in parameters are declared here.
    lkout_HR  <- fixed(log(11.6))  ; label("First-order dissipation rate for HR, KOUTHR (1/h)")   # Fu 2023 Suppl. S1 $THETA TH8 (11.6 FIX)
    lkout_SV  <- fixed(log(0.126)) ; label("First-order dissipation rate for SV, KOUTSV (1/h)")   # Fu 2023 Suppl. S1 $THETA TH9 (0.126 FIX)
    lkout_TPR <- fixed(log(3.58))  ; label("First-order dissipation rate for TPR, KOUTTPR (1/h)") # Fu 2023 Suppl. S1 $THETA TH10 (3.58 FIX)

    # ---- HR-on-SV coupling constant ------------------------------
    # Describes the magnitude of the direct HR-to-SV coupling in the
    # $DES definition SV = A(3)*(1 - HRSV*LOG(A(2)/BSLHR)). At
    # baseline (A(2) = BSLHR) the term is 1 and SV = A(3) = BSLSV.
    # When HR rises above baseline, actual SV decreases proportional
    # to log(HR/BSLHR) with slope HRSV.
    lhrsv <- fixed(log(0.312))    ; label("HR-on-SV direct coupling constant, HRSV (unitless)")  # Fu 2023 Suppl. S1 $THETA TH11 (0.312 FIX)

    # ---- Circadian rhythm parameters (fixed) ---------------------
    # HR circadian rhythm: cosine with amplitude AMPHR and horizontal
    # displacement HORHR (period 24 h). TPR circadian rhythm: cosine
    # with amplitude AMPHR * AMPTPRratio and horizontal displacement
    # HORTPR (Fu 2023 Suppl. S1 $PK: AMP2 = AMP1 * THETA(15) where
    # THETA(15) = 1 FIX; the ratio is 1 in the published model).
    # Displacements and amplitudes are kept on the linear scale
    # because phase offsets can (in principle) be signed and cosine
    # amplitudes below the mean production rate are permissible.
    hor_HR        <- fixed(8.73)   ; label("Horizontal displacement of HR circadian rhythm, HORHR (h)")           # Fu 2023 Suppl. S1 $THETA TH12 (8.73 FIX)
    amp_HR        <- fixed(0.0918) ; label("Amplitude (fraction of production) of HR circadian rhythm, AMPHR")    # Fu 2023 Suppl. S1 $THETA TH13 (0.0918 FIX)
    hor_TPR       <- fixed(19.3)   ; label("Horizontal displacement of TPR circadian rhythm, HORTPR (h)")         # Fu 2023 Suppl. S1 $THETA TH14 (19.3 FIX)
    amp_TPR_ratio <- fixed(1)      ; label("Ratio of TPR-circadian amplitude to HR-circadian amplitude, AMPTPRratio (unitless)")  # Fu 2023 Suppl. S1 $THETA TH15 (1 FIX)

    # ---- Inter-individual variability (log-normal, fixed) --------
    # Snelder-derived IIV on the three baseline hemodynamic
    # variables, encoded as log-normal variance ($OMEGA in
    # Supplemental S1, BSLX = THETAX * EXP(ETAX) in $PK). Fu 2023
    # Section 2.2 states that "the system-specific parameters were
    # fixed to the parameter values from the published CVS model"
    # in the identifiability-focused SSE scenarios; the IIVs are
    # part of that system parameter set (they belong to the rat
    # physiologic baseline distribution, not to the drug-specific
    # parameters that were estimated). Wrapped in fixed() to
    # preserve the "held constant during SSE re-estimation"
    # provenance.
    etalrbase_HR  ~ fixed(0.00372) # Fu 2023 Suppl. S1 $OMEGA (0.00372) IIV_BSLHR  ; log-scale variance
    etalrbase_MAP ~ fixed(0.00137) # Fu 2023 Suppl. S1 $OMEGA (0.00137) IIV_BSLMAP ; log-scale variance
    etalrbase_CO  ~ fixed(0.0515)  # Fu 2023 Suppl. S1 $OMEGA (0.0515)  IIV_BSLCO  ; log-scale variance

    # ---- Residual error (proportional, per output) ---------------
    # $SIGMA in the supplement gives NONMEM VARIANCE for the
    # proportional residual EPS on each of HR, CO, MAP (Y = IPRED *
    # (1 + EPS(k))). nlmixr2 propSd_<output> is expressed as the
    # standard deviation, so propSd = sqrt(sigma^2). All three are
    # fixed in the SSE reference set (system-specific residual
    # errors from Snelder that Fu 2023 uses as the simulation
    # reference).
    propSd_HR  <- fixed(sqrt(0.006084)) ; label("Proportional residual SD on HR (fraction)")  # Fu 2023 Suppl. S1 $SIGMA 0.006084 (variance); sqrt = 0.0780
    propSd_CO  <- fixed(sqrt(0.004761)) ; label("Proportional residual SD on CO (fraction)")  # Fu 2023 Suppl. S1 $SIGMA 0.004761 (variance); sqrt = 0.0690
    propSd_MAP <- fixed(sqrt(0.0036))   ; label("Proportional residual SD on MAP (fraction)") # Fu 2023 Suppl. S1 $SIGMA 0.0036   (variance); sqrt = 0.0600
  })

  model({
    # ---- Derived scalar constants --------------------------------
    # Circadian rhythm period is 24 h (Fu 2023 Suppl. S1 $PK: PER =
    # 24). The cosine argument is (2*pi/PER) * (t + hor) so that a
    # positive hor shifts the peak to earlier time.
    circ_per <- 24

    # ---- Bare (linear-scale) hemodynamic parameters --------------
    # Log-normal IIV attaches only to the three baselines that
    # $PK explicitly randomizes (BSLHR, BSLMAP, BSLCO); all other
    # system parameters have no IIV in the supplement.
    fb            <- exp(lfb)
    kel           <- exp(lkel)
    ec50          <- exp(lec50)
    bsl_HR        <- exp(lrbase_HR  + etalrbase_HR)
    bsl_MAP       <- exp(lrbase_MAP + etalrbase_MAP)
    bsl_CO        <- exp(lrbase_CO  + etalrbase_CO)
    kout_HR       <- exp(lkout_HR)
    kout_SV       <- exp(lkout_SV)
    kout_TPR      <- exp(lkout_TPR)
    hrsv          <- exp(lhrsv)

    # ---- Derived hemodynamic baselines ---------------------------
    # SV baseline = CO / HR (Fu 2023 Suppl. S1 $PK: BSLSV =
    # BSLCO/BSLHR). TPR baseline = MAP / CO (BSLTPR = BSLMAP/BSLCO).
    bsl_SV  <- bsl_CO  / bsl_HR
    bsl_TPR <- bsl_MAP / bsl_CO

    # ---- Steady-state kin (production) rate constants ------------
    # Solving the ODE at steady state under baseline conditions
    # gives kin_X = kout_X * BSL_X / (1 - fb * BSLMAP) for each of
    # the three primary states (Fu 2023 Suppl. S1 $PK: KINHR =
    # KOUTHR*BSLHR/(1 - FB*BSLMAP), KINSV = KOUTSV*BSLSV/(1 -
    # FB*BSLMAP), KINTPR = KOUTTPR*BSLTPR/(1 - FB*BSLMAP)).
    fb_bsl  <- 1 - fb * bsl_MAP
    kin_HR  <- kout_HR  * bsl_HR  / fb_bsl
    kin_SV  <- kout_SV  * bsl_SV  / fb_bsl
    kin_TPR <- kout_TPR * bsl_TPR / fb_bsl

    # ---- Circadian rhythm signals --------------------------------
    # Cosine functions with amplitude relative to the production
    # rate (multiplicative on kin_HR and kin_TPR, respectively).
    cs_HR  <- amp_HR                 * cos(2 * pi * (t + hor_HR)  / circ_per)
    cs_TPR <- amp_HR * amp_TPR_ratio * cos(2 * pi * (t + hor_TPR) / circ_per)

    # ---- Derived hemodynamic outputs from ODE states -------------
    # Actual SV incorporates the direct HR-on-SV coupling: SV =
    # SVT * (1 - HRSV * log(HR/BSLHR)) (Fu 2023 Suppl. S1 $DES).
    # CO = HR * SV and MAP = CO * TPR.
    sv_actual <- svt * (1 - hrsv * log(hr / bsl_HR))
    co        <- hr * sv_actual
    map       <- co * tpr

    # ---- Drug effect (Emax on HR-production) ---------------------
    # The on-disk NONMEM control stream places the drug effect on
    # the HR-production term as (1 + emax * A(1)/(EC50 + A(1))).
    # emax is signed (-1 = inhibition, +1 or +10 = stimulation)
    # so the factor spans (0, 2] for stimulation or [0, 1) for
    # inhibition depending on the sign. To relocate the effect to
    # SV or TPR (Fu 2023 SSE scenarios of Section 2.3), move the
    # e_drug factor into the corresponding d/dt(svt) or d/dt(tpr)
    # equation.
    e_drug <- emax * central / (ec50 + central)

    # ---- ODE system ---------------------------------------------
    # 1. Hypothetical-drug PK (one-compartment iv bolus). The
    #    central compartment carries the drug amount; because the
    #    supplement uses A(1) directly in the Emax expression, the
    #    "concentration" in the drug-effect term is numerically
    #    equal to the central-compartment amount (implicit V = 1).
    d/dt(central) <- -kel * central

    # 2. Heart-rate indirect-response ODE (state units: beats/min).
    #    Production combines the circadian rhythm, MAP negative
    #    feedback, and the drug effect on HR-production; dissipation
    #    is first-order at kout_HR.
    d/dt(hr) <- kin_HR * (1 + cs_HR) * (1 - fb * map) * (1 + e_drug) - kout_HR * hr

    # 3. SV-turnover indirect-response ODE (state units: mL/beat).
    #    No circadian rhythm on SV-production (the supplement does
    #    not include one); MAP feedback is applied.
    d/dt(svt) <- kin_SV * (1 - fb * map) - kout_SV * svt

    # 4. TPR indirect-response ODE (state units: mmHg*min/mL).
    #    Circadian rhythm on TPR-production and MAP feedback.
    d/dt(tpr) <- kin_TPR * (1 + cs_TPR) * (1 - fb * map) - kout_TPR * tpr

    # ---- Initial conditions (baseline hemodynamics) -------------
    central(0) <- 0
    hr(0)      <- bsl_HR
    svt(0)     <- bsl_SV
    tpr(0)     <- bsl_TPR

    # ---- Observations and residual error -------------------------
    # HR, CO, MAP are the three measured hemodynamic readouts of the
    # SSE scenarios (Fu 2023 Section 2.2); each has its own
    # proportional residual SD. Cc (drug amount in central) is
    # exposed for downstream plotting but is not a paper observation.
    HR  <- hr
    CO  <- co
    MAP <- map
    Cc  <- central

    HR  ~ prop(propSd_HR)
    CO  ~ prop(propSd_CO)
    MAP ~ prop(propSd_MAP)
  })
}
