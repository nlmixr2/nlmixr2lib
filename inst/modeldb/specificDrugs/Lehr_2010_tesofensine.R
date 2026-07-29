Lehr_2010_tesofensine <- function() {
  description <- "Joint parent (tesofensine) + metabolite (M1, CYP3A4-formed) population PK and effect-compartment PK/PD model in mild Alzheimer's disease (Lehr 2010 Phase IIa fit; 62 patients across two 4-week placebo-controlled studies). Parent is one-compartment with first-order absorption (ka FIXED from upstream Phase I popPK) and parallel elimination through a metabolite-formation arm (CL_met = parent -> M1 flux) and a non-formation arm (CL_non-met = elimination via routes other than M1 formation). M1 is one-compartment with apparent volume FIXED at 0.768-fold of the parent apparent volume (mouse-derived ratio, Lehr 2010 ref 17). Tesofensine and M1 each drive their own effect compartment (shared keo FIXED at a small value, equivalent to a long effect-equilibration half-life); the combined drug effect on ADAS-Cog uses an extended Emax with competitive interaction in which the M1 effect-compartment concentration is divided by 5 to reflect the in-vivo M1 potency one-fifth that of the parent (Lehr 2010 Methods, ref 17). The ADAS-Cog observation equals the sum of drug, placebo, and disease-progression contributions (change from each subject's baseline). The placebo bi-exponential model (onset rate keq, offset rate kel_pla, scaling beta_pla) is fully FIXED to literature values from a published large-AD-cohort placebo model (Lehr 2010 ref 34); the linear disease-progression slope is FIXED at 6 ADAS-Cog points/year (Lehr 2010 ref 26). Emax is negative (a clinically meaningful ADAS-Cog improvement is a score reduction); the sign is applied inside model() while |Emax| is the ini-scale magnitude carried with multiplicative IIV."
  reference <- paste(
    "Lehr T, Staab A, Trommeshauser D, Schaefer HG, Kloft C.",
    "Quantitative Pharmacology Approach in Alzheimer's Disease:",
    "Efficacy Modeling of Early Clinical Data to Predict Clinical",
    "Outcome of Tesofensine. AAPS J. 2010;12(2):117-129.",
    "doi:10.1208/s12248-009-9164-6.",
    sep = " "
  )
  vignette <- "Lehr_2010_tesofensine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species         = "human",
    n_subjects      = 62L,
    n_studies       = 2L,
    age_range       = "70-80 years (median 70)",
    age_median      = "70 years",
    weight_range    = "54-129 kg (median 80)",
    weight_median   = "80 kg",
    sex_female_pct  = 44,
    race_ethnicity  = c(White = 96.8, Black = 3.2),
    disease_state   = "mild Alzheimer's disease (median ADAS-Cog at baseline 8.5; range 3.3-23.7)",
    dose_range      = "0.125-1.0 mg/day oral; loading of 0.5-2.0 mg/day for the first 3 days then maintenance for 25 days",
    renal_function  = "normal hepatic function; median creatinine clearance 77 mL/min (range 41-136)",
    regions         = "(not reported)",
    n_observations  = "357 tesofensine and 341 M1 plasma concentrations from 44 active-treated patients; 176 active and 72 placebo ADAS-Cog measurements",
    notes           = "Demographics from Lehr 2010 Results, Data Base section: 44 active + 18 placebo treated patients with mild AD, median ADAS-Cog at baseline 8.5 points (range 3.3-23.7). The race composition is 60 Caucasians + 2 African-Americans (paper text). The model was developed on the Phase IIa dataset (this file uses the Phase IIa column of Table I) and subsequently re-estimated on the combined Phase IIa+IIb dataset (430 additional patients); see the validation vignette for the combined-dataset re-estimation comparison."
  )

  ini({
    # Structural PK parameters - reference values are the final Phase IIa
    # NONMEM estimates from Lehr 2010 Table I (column 'Phase IIa studies
    # dataset / NONMEM / Value'). Fixed-vs-estimated encoding mirrors the
    # Table I 'RSE,%' column: FIX entries become fixed() wrappers.

    lka        <- fixed(log(0.69))  ; label("First-order tesofensine absorption rate ka (1/h)")             # Lehr 2010 Table I: ka = 0.69 1/h FIX (Phase IIa); FIXED to a value from an unpublished Phase I popPK analysis (Methods, Population PK Model section)
    lcl_nonmet <- log(1.31)         ; label("Apparent non-formation clearance CL_non-met/F of tesofensine (L/h); elimination via routes other than M1 formation, carries the parent CL IIV") # Lehr 2010 Table I: CL_non-met/F = 1.31 L/h, RSE 5.5%
    lcl_met    <- log(0.416)        ; label("Apparent formation clearance CL_met/F of tesofensine to M1 (L/h); fixed effect, no IIV per the source paper") # Lehr 2010 Table I: CL_met/F = 0.416 L/h, RSE 6.9%
    lvc        <- log(720)          ; label("Apparent tesofensine central volume V2/F (L)")                 # Lehr 2010 Table I: V2/F = 720 L, RSE 4.9%
    lcl_m1     <- log(1.17)         ; label("Apparent M1 elimination clearance CL_M1/F (L/h)")              # Lehr 2010 Table I: CL_M1/F = 1.17 L/h, RSE 8.2%
    lvc_m1     <- fixed(log(553))   ; label("Apparent M1 central volume V3/F (L); FIXED at 0.768-fold of V2/F per Lehr 2010 Table I footnote d (mouse-derived ratio, ref 17)") # Lehr 2010 Table I: V3/F = 553 L FIX (Phase IIa)

    # PD parameters
    lkeo       <- fixed(log(0.0001)); label("Effect-compartment equilibration rate constant keo (1/h); FIXED at the upper boundary of the operator-insensitive plateau identified in the sensitivity analysis (Methods, Population PK/PD Model section)") # Lehr 2010 Table I: keo = 0.0001 1/h FIX
    lemax      <- log(1.46)         ; label("Log absolute Emax for ADAS-Cog improvement (paper Emax = -1.46 ADAS-Cog points; sign applied in model())") # Lehr 2010 Table I: E_MAX = -1.46, RSE 29.3%
    lec50      <- log(0.0139)       ; label("Log EC50 (ng/mL) producing 50% of |Emax|; tesofensine effect-compartment scale (M1 effect-compartment concentration is scaled by 1/5 to share this EC50, per the competitive-interaction extended Emax structure)") # Lehr 2010 Table I: EC50 = 0.0139 ng/mL, RSE 49.6%

    # Placebo model - all parameters FIXED to literature values from a
    # published large-AD-cohort placebo model (Lehr 2010 ref 34); the
    # paper notes the parameters could not be estimated from the sparse
    # 4-week Phase IIa placebo data (only 72 ADAS-Cog measurements from
    # 18 placebo patients).
    lkeq       <- fixed(log(0.00183))    ; label("Log placebo onset rate constant keq (1/h); FIXED from literature placebo model (Lehr 2010 ref 34)") # Lehr 2010 Table I: keq = 0.00183 1/h FIX
    lkel_pla   <- fixed(log(0.000473))   ; label("Log placebo offset rate constant kel_pla (1/h); paper symbol is 'kel' but renamed here to lkel_pla to avoid collision with the K-PD canonical lkel; FIXED from literature placebo model") # Lehr 2010 Table I: kel = 0.000473 1/h FIX
    lbeta_pla  <- fixed(log(1.42))       ; label("Log absolute placebo scaling magnitude (paper beta = -1.42 ADAS-Cog points scaling; sign applied in model()); FIXED from literature placebo model") # Lehr 2010 Table I: beta = -1.42 FIX

    # Disease progression - linear, FIXED at the literature value used
    # by Lehr 2010 (6 ADAS-Cog points/year per Lehr 2010 ref 26).
    # Stored in (ADAS-Cog points per hour) for unit consistency with
    # the rest of the model.
    ldp_rate   <- fixed(log(6 / (365.25 * 24))) ; label("Log linear disease-progression rate (ADAS-Cog points/hour); FIXED at literature 6 points/year (Lehr 2010 ref 26)") # Lehr 2010 Disease Progression and Placebo Effect Model section

    # IIV - Lehr 2010 Table I (Phase IIa column) reports each random
    # effect as a percent CV on the parameter. The log-normal variance
    # is omega^2 = log(1 + (CV%/100)^2). The Table I footnote a states
    # 'RSE is given on the variance scale' (i.e., the RSE column applies
    # to omega^2, not CV%); the CV% values themselves are the random-
    # effect magnitudes used here.
    # Block correlation: NONMEM $OMEGA BLOCK on (etalcl_nonmet, etalvc)
    # with rho = 0.72 (Lehr 2010 Table I: 'Corr CL non-met /F_V2/F').
    # Covariance = rho * sqrt(omega_cl_nonmet^2 * omega_vc^2)
    #            = 0.72 * sqrt(0.1640 * 0.0878)
    #            = 0.72 * 0.1200 = 0.0863
    etalcl_nonmet + etalvc ~ c(0.1640, 0.0863, 0.0878)   # CL_non-met CV%=42.2 -> omega^2=log(1+0.422^2)=0.1640; V2 CV%=30.3 -> 0.0878; rho=0.72 -> cov=0.0863
    etalcl_m1 ~ 0.0444                                    # CL_M1 CV%=21.3 -> omega^2=log(1+0.213^2)=0.0444
    etalvc_m1 ~ 0.1664                                    # V3 CV%=42.5 -> omega^2=log(1+0.425^2)=0.1664
    etalemax  ~ 0.6700                                    # Emax CV%=97.7 -> omega^2=log(1+0.977^2)=0.6700; multiplicative IIV on |Emax|, sign preserved by the lognormal form
    etalbeta_pla ~ fixed(0.9698)                          # beta CV%=128 FIX -> omega^2=log(1+1.28^2)=0.9698; carried from the literature placebo model (Lehr 2010 ref 34)

    # Residual error
    propSd          <- 0.149  ; label("Proportional residual error on tesofensine plasma concentration (fraction)") # Lehr 2010 Table I: prop err tesofensine = 14.9%, RSE 23.5%
    propSd_m1       <- 0.173  ; label("Proportional residual error on M1 plasma concentration (fraction)")          # Lehr 2010 Table I: prop err M1 = 17.3%, RSE 36.9%
    addSd_ADAS_cog  <- 2.3    ; label("Additive residual SD on the ADAS-Cog change-from-baseline observation (ADAS-Cog points)") # Lehr 2010 Table I: add err ADAS-Cog = +/-2.3, RSE 19.5%
  })

  model({
    # Individual structural parameters
    ka        <- exp(lka)
    cl_nonmet <- exp(lcl_nonmet + etalcl_nonmet)
    cl_met    <- exp(lcl_met)
    cl_m1     <- exp(lcl_m1 + etalcl_m1)
    vc        <- exp(lvc + etalvc)
    vc_m1     <- exp(lvc_m1 + etalvc_m1)
    keo       <- exp(lkeo)
    # Sign applied here (paper Emax = -1.46); the lognormal IIV on
    # |Emax| preserves the sign across individuals.
    emax      <- -exp(lemax + etalemax)
    ec50      <- exp(lec50)
    keq       <- exp(lkeq)
    kel_pla   <- exp(lkel_pla)
    beta_pla  <- -exp(lbeta_pla + etalbeta_pla)
    dp_rate   <- exp(ldp_rate)

    # Plasma concentrations (ng/mL). Dose is in mg, vc in L, so
    # central/vc has units mg/L = ug/mL; multiply by 1000 to express
    # in ng/mL (the source paper's reporting scale).
    Cc    <- (central    / vc)    * 1000
    Cc_m1 <- (central_m1 / vc_m1) * 1000

    # PK ODE system - parent is one-compartment with first-order absorption
    # and parallel elimination (formation of M1 + non-formation routes);
    # M1 is one-compartment with linear elimination, fed by the
    # cl_met/vc * central formation flux out of the parent's central
    # compartment.
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  ka * depot - (cl_met + cl_nonmet) / vc * central
    d/dt(central_m1) <-  cl_met / vc * central - cl_m1 / vc_m1 * central_m1

    # Effect compartments (concentrations directly modeled in ng/mL,
    # paired with the parent and metabolite plasma compartments).
    d/dt(effect)    <- keo * (Cc    - effect)
    d/dt(effect_m1) <- keo * (Cc_m1 - effect_m1)

    # Drug effect on ADAS-Cog (extended Emax with competitive
    # interaction; M1 potency one-fifth of tesofensine per Lehr 2010
    # Methods, ref 17). The numerator Ce/EC50 contributions add (parent
    # plus M1 weighted by 1/5); the denominator is the standard
    # Emax(1+sum) form.
    ce_norm <- effect / ec50 + effect_m1 / (ec50 * 5)
    drug    <- emax * ce_norm / (1 + ce_norm)

    # Placebo bi-exponential model (Lehr 2010 Eq. 1, ref 34):
    # PLACEBO = (beta * keq) / (keq - kel) * (exp(-kel*time) - exp(-keq*time)).
    # Numerically stable because keq (0.00183) > kel_pla (0.000473)
    # so the denominator never vanishes. At time = 0 the placebo
    # contribution is exactly 0.
    placebo <- (beta_pla * keq) / (keq - kel_pla) *
               (exp(-kel_pla * time) - exp(-keq * time))

    # Linear disease progression added to the overall ADAS-Cog
    # change-from-baseline. dp_rate is in ADAS-Cog points/hour;
    # multiplying by time (hours) gives the cumulative progression.
    dp <- dp_rate * time

    # ADAS-Cog change from baseline = drug + placebo + DP. The source
    # paper writes ADAS-Cog = Baseline + Drug + Placebo + Disease
    # Progression; this model emits the change-from-baseline component
    # so that downstream users can add an individual baseline value
    # explicitly (the baseline is per-subject and not estimated in the
    # source paper).
    ADAS_cog <- drug + placebo + dp

    Cc       ~ prop(propSd)
    Cc_m1    ~ prop(propSd_m1)
    ADAS_cog ~ add(addSd_ADAS_cog)
  })
}
