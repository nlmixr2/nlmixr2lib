Wright_2013_allopurinol <- function() {
  description <- "Sequential parent-metabolite population pharmacokinetic model for oral and intravenous allopurinol and its principal long-lived active metabolite oxypurinol in patients with chronic gout and healthy volunteers, pooled from three Christchurch (New Zealand) clinical studies (dose-escalation, furosemide-interaction, vitamin-C interaction) and two digitally-extracted historic intravenous-allopurinol studies. Structural model: 2-compartment first-order absorption disposition for allopurinol (parent) with 80 percent molar conversion to a 1-compartment oxypurinol (metabolite). Fat-free-mass (FFM) allometric scaling on all clearance parameters (exponent 0.75) and volume parameters (exponent 1.0) referenced at FFM = 70 kg. Oxypurinol clearance is decomposed into non-renal (0.178 L/h/70 kg FFM) plus renal (0.777 L/h/70 kg FFM per creatinine clearance of 6 L/h) components, with a 39 percent reduction of the renal component in patients on concomitant diuretic therapy. All allopurinol structural parameters, allometric exponents, oral bioavailability, and the molar conversion fraction were held FIXED during the sequential fit; oxypurinol structural parameters, IIV, and residual error were estimated. Wright 2013 Table 2 reports residual epsilon components under a '(CV%)' column label whose numeric values are interpreted here as SD on fraction scale (see the model's validation vignette Assumptions and deviations section for the rationale)."
  reference <- paste(
    "Wright DFB, Stamp LK, Merriman TR, Barclay ML, Duffull SB, Holford NHG.",
    "The population pharmacokinetics of allopurinol and oxypurinol in patients with gout.",
    "Eur J Clin Pharmacol. 2013;69(7):1411-1421.",
    "doi:10.1007/s00228-013-1478-8.",
    sep = " "
  )
  vignette <- "Wright_2013_allopurinol"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass computed via the Janmahasatian 2005 formula from total body weight, height and sex",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Reference FFM = 70 kg per the Anderson & Holford 2008 allometric convention (Wright 2013 Methods, 'Covariate model'). Allometric exponent 0.75 on CL_allo, Q_allo, CL_oxy and 1.0 on V1_allo, V2_allo, V_oxy. Range 35 to 99 kg with median 70 kg across the pooled 104 subjects (Wright 2013 Table 1).",
      source_name        = "FFM"
    ),
    CRCL = list(
      description        = "Creatinine clearance standardised to 70 kg body weight (Cockcroft-Gault formula with Anderson-Holford allometric normalisation), NOT BSA-normalised to 1.73 m^2",
      units              = "L/h per 70 kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Reference CLcr_STD = 6 L/h per 70 kg (equivalent to 100 mL/min per 70 kg body weight), Wright 2013 Methods 'Covariate model' equation for RF. Range 1.0 to 8.5 L/h (median 4.3 L/h) across the pooled dataset (Wright 2013 Table 1). The renal component of oxypurinol clearance is proportional to CRCL / 6. Wright's CLcr is scaled per 70 kg body weight rather than to 1.73 m^2 body surface area, in keeping with the Anderson-Holford PK-scaling convention.",
      source_name        = "CLcr"
    ),
    CONMED_DIURETIC = list(
      description        = "Concomitant diuretic therapy indicator (any diuretic class)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant diuretic use)",
      notes              = "Time-fixed at baseline. Coded 1 for subjects on any diuretic class (loop, thiazide, thiazide-like, or potassium-sparing). Wright 2013 Table 1 reports 29 of 104 pooled subjects (28 percent) on diuretics. Applied as a multiplicative 0.61 factor on the RENAL component of oxypurinol clearance only (Wright 2013 page 6 final-model equation for CL_oxy), corresponding to a 39 percent reduction in renal CL_oxy in diuretic users.",
      source_name        = "Diuretic"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104,
    n_studies      = 5,
    age_range      = "adult",
    weight_range   = "49 - 171 kg (total body weight; Wright 2013 Table 1)",
    ffm_range      = "35 - 99 kg (fat-free mass; Wright 2013 Table 1)",
    sex_ratio      = "93 male : 11 female (Wright 2013 Table 1)",
    disease_state  = "Chronic gout in 92 patients from three Christchurch (New Zealand) rheumatology cohorts led by author LKS, plus 12 healthy volunteers from two digitally-extracted historic intravenous-allopurinol studies (Appelbaum et al. and Tittel & Breithaupt).",
    dose_range     = "Oral allopurinol 100 to 600 mg once daily in the Christchurch cohorts; digitally-extracted intravenous-allopurinol single-dose schedules in the two historic cohorts. All concentrations reported as umol/L (Wright 2013 ESM 1).",
    regions        = "New Zealand (Christchurch); historic IV cohorts contribute pharmacokinetic profiles digitally extracted from published figures (Appelbaum, Tittel-Breithaupt).",
    notes          = "Pooled patient-level data from five studies. Christchurch (n = 92 gout patients): dose-escalation study (n = 74, 442 obs), furosemide-interaction study (n = 10, 111 obs), vitamin-C interaction study (n = 16, 141 obs). Digitally-extracted intravenous cohorts (n = 12 healthy volunteers). Median CLcr 4.3 L/h (range 1.0 - 8.5); median FFM 70 kg (range 35 - 99). 29 of 104 subjects on concomitant diuretic therapy. 40 of 104 on ACE inhibitors, 43 on beta blockers, 18 on NSAIDs, 24 on urate-lowering therapy (Wright 2013 Table 1); of these concomitant categories only diuretic use was retained in the final model."
  )

  ini({
    # --- Allopurinol structural parameters (all FIXED in the sequential
    #     parent-metabolite step from Wright 2013's step-2 allopurinol-only
    #     fit; Wright 2013 Table 2 lists identical values in the 'Base
    #     sequential model' and 'Final sequential model' columns because the
    #     allopurinol layer was held fixed while oxypurinol parameters were
    #     estimated).
    lka       <- fixed(log(1.6));  label("Allopurinol first-order absorption rate constant ka (1/h)")                          # Wright 2013 Table 2, ka = 1.6 (fixed)
    lcl_allo  <- fixed(log(49.6)); label("Allopurinol clearance CL_allo (L/h) at FFM = 70 kg")                                 # Wright 2013 Table 2, CL_allo = 49.6 (fixed)
    lvc_allo  <- fixed(log(11.4)); label("Allopurinol central volume V1_allo (L) at FFM = 70 kg")                              # Wright 2013 Table 2, V1_allo = 11.4 (fixed)
    lq_allo   <- fixed(log(142));  label("Allopurinol inter-compartmental clearance Q_allo (L/h) at FFM = 70 kg")              # Wright 2013 Table 2, Q = 142 (fixed)
    lvp_allo  <- fixed(log(90.7)); label("Allopurinol peripheral volume V2_allo (L) at FFM = 70 kg")                           # Wright 2013 Table 2, V2_allo = 90.7 (fixed)
    lfdepot   <- fixed(log(0.85)); label("Oral bioavailability F of allopurinol (unitless; 0.85 per prior literature)")       # Wright 2013 Table 2, F = 0.85 (fixed)

    # --- Oxypurinol structural parameters (final sequential estimates from
    #     the Wright 2013 page-6 CL_oxy and V_oxy equations and Table 2
    #     'Final sequential model' column). The renal + non-renal
    #     decomposition and the diuretic fractional-effect multiplier are
    #     encoded on linear scale to preserve the additive form of the CL
    #     equation; log-normal IIV is applied multiplicatively to the total.
    cl_oxy_nonren <- 0.178; label("Oxypurinol non-renal clearance CL_oxy,nonren (L/h) at FFM = 70 kg")                                   # Wright 2013 page 6 equation, "0.178"
    cl_oxy_renref <- 0.777; label("Oxypurinol renal clearance CL_oxy,renal (L/h) at FFM = 70 kg and CLcr = 6 L/h reference")             # Wright 2013 page 6 equation, "0.777"
    lv_oxy        <- log(41.4); label("Oxypurinol central volume V_oxy (L) at FFM = 70 kg")                                              # Wright 2013 Table 2, V_oxy = 41.4 (final)

    # --- Allometric exponents and metabolic-fraction parameters (all FIXED).
    e_ffm_cl <- fixed(0.75);       label("Allometric exponent on CL_allo, Q_allo, CL_oxy versus FFM / 70 kg")                    # Wright 2013 Methods, "Clearance was allometrically scaled to an exponent of 0.75"
    e_ffm_vc <- fixed(1.0);        label("Allometric exponent on V1_allo, V2_allo, V_oxy versus FFM / 70 kg")                    # Wright 2013 Methods, "volume to an exponent of 1"
    fmet_oxy <- fixed(0.80);       label("Molar fraction of allopurinol converted to oxypurinol (unitless; 0.80 per literature)")  # Wright 2013 Results, "fixed at 80% based on values previously reported [3, 14]"

    # --- Diuretic covariate effect on the RENAL component of oxypurinol
    #     clearance. Encoded as a fractional multiplier (0.61 = 39 percent
    #     reduction of CL_oxy,renal in diuretic users; reference = 1.0).
    e_conmed_diuretic_cl_oxy_renal <- 0.61; label("Diuretic fractional-effect multiplier on renal CL_oxy (unitless; reference 1.0; diuretic user 0.61)")  # Wright 2013 Table 2, "Diuretic fractional effect 0.61"; page 6 equation

    # --- Between-subject variability (log-normal on parameter scale). Wright
    #     2013 reports BSV as approximate CV percent; converted here to the
    #     log-normal omega variance via omega^2 = log(CV^2 + 1) with CV as a
    #     fraction (standard log-normal moment-match).
    etalcl_allo  ~ fix(0.1029)  # BSV_CL_allo = 32.9% CV (fixed from Wright 2013 step 2). omega^2 = log(0.329^2 + 1) = 0.1029
    etalvp_allo  ~ fix(0.0800)  # BSV_V2_allo = 28.9% CV (fixed). omega^2 = log(0.289^2 + 1) = 0.0800
    # Oxypurinol IIV block: correlated CL_oxy - V_oxy with Corr(CL,V) = 0.148
    # (Wright 2013 Table 2 final). Diagonal entries are omega^2:
    #   CL_oxy: 28.2% CV -> omega^2 = log(0.282^2 + 1) = 0.07648
    #   V_oxy:  14.9% CV -> omega^2 = log(0.149^2 + 1) = 0.02195
    # Off-diagonal: corr x sqrt(omega1^2 x omega2^2) = 0.148 x sqrt(0.07648 x 0.02195) = 0.006064
    etacl_oxy + etalv_oxy ~ c(
      0.07648,
      0.006064, 0.02195
    )

    # --- Residual error (combined additive + proportional per analyte).
    #     Wright 2013 Table 2's '(CV%)' header column values are interpreted
    #     here as the SD on FRACTION scale, not literal percent, in keeping
    #     with the Fig 3b prediction-corrected VPC spread and with the
    #     convention of the same group's Wright 2016 follow-up paper (BJCP
    #     doi:10.1111/bcp.12799), which reports oxypurinol proportional SD
    #     as 0.20 = 20% CV under the identical Table-column heading style.
    #     See vignette Assumptions and deviations for the full rationale.
    addSd      <- fixed(0.0001);  label("Allopurinol additive residual SD (umol/L)")                                             # Wright 2013 Table 2, epsilon_add allo = 0.0001 (fixed)
    propSd     <- fixed(0.74);    label("Allopurinol proportional residual SD (fraction; interpreted as 74% CV)")                # Wright 2013 Table 2, epsilon_prop allo = 0.74 (fixed)
    addSd_oxy  <- 0.21;           label("Oxypurinol additive residual SD (umol/L)")                                                     # Wright 2013 Table 2, epsilon_add oxy = 0.21
    propSd_oxy <- 0.00004;        label("Oxypurinol proportional residual SD (fraction; near-zero after fit)")                          # Wright 2013 Table 2, epsilon_prop oxy = 0.00004
  })

  model({
    # Molecular weights (g/mol) used to convert the mass flux of allopurinol
    # leaving its central compartment (mg/h) to the mass flux of oxypurinol
    # entering the metabolite central compartment (mg/h). fmet_oxy = 0.80 is
    # a MOLAR fraction, so the mass conversion is 0.80 * (mwOxy / mwAllo).
    mwAllo      <- 136.11
    mwOxy       <- 152.11
    molarFactor <- mwOxy / mwAllo

    # Individual PK parameters. Allometric FFM scaling per Anderson & Holford
    # (FFM / 70 kg)^0.75 on clearance and (FFM / 70 kg)^1.0 on volume.
    ka       <- exp(lka)
    cl_allo  <- exp(lcl_allo + etalcl_allo) * (FFM / 70)^e_ffm_cl
    vc_allo  <- exp(lvc_allo)                * (FFM / 70)^e_ffm_vc
    q_allo   <- exp(lq_allo)                 * (FFM / 70)^e_ffm_cl
    vp_allo  <- exp(lvp_allo + etalvp_allo) * (FFM / 70)^e_ffm_vc

    # Oxypurinol clearance: linear sum of non-renal + renal components. The
    # renal component is proportional to (CRCL / 6), scaled by the diuretic
    # multiplier raised to the CONMED_DIURETIC indicator. The whole clearance
    # is allometrically scaled by (FFM / 70)^0.75 (Wright 2013 page 6 final
    # equation). Log-normal BSV is applied multiplicatively to the total.
    cl_oxy_renal <- cl_oxy_renref * (CRCL / 6) * (e_conmed_diuretic_cl_oxy_renal^CONMED_DIURETIC)
    cl_oxy_typ   <- (cl_oxy_nonren + cl_oxy_renal) * (FFM / 70)^e_ffm_cl
    cl_oxy       <- cl_oxy_typ * exp(etacl_oxy)
    v_oxy        <- exp(lv_oxy + etalv_oxy) * (FFM / 70)^e_ffm_vc

    # Micro-rate constants (1/h).
    kel_allo <- cl_allo / vc_allo
    k12_allo <- q_allo  / vc_allo
    k21_allo <- q_allo  / vp_allo
    kel_oxy  <- cl_oxy  / v_oxy

    # ODE system. All compartment amounts are in mg (allopurinol) and mg
    # (oxypurinol). Oxypurinol formation flux = fmet_oxy (mole fraction) x
    # kel_allo x central (allopurinol mass leaving) x molarFactor (mg-of-oxy
    # per mg-of-allo).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel_allo * central -
                          k12_allo * central + k21_allo * peripheral1
    d/dt(peripheral1) <-  k12_allo * central - k21_allo * peripheral1
    d/dt(central_oxy) <-  fmet_oxy * kel_allo * central * molarFactor -
                          kel_oxy * central_oxy

    # Oral bioavailability of allopurinol (fixed at 0.85 per prior literature;
    # Wright 2013 Table 2). Bioavailability targets the depot compartment;
    # intravenous dosing bypasses the depot and enters `central` directly.
    f(depot) <- exp(lfdepot)

    # Plasma observations in umol/L. Compartment amount (mg) divided by
    # volume (L) gives mg/L; multiply by 1000 / MW to convert to umol/L.
    Cc     <- 1000 * central     / (vc_allo * mwAllo)
    Cc_oxy <- 1000 * central_oxy / (v_oxy   * mwOxy)

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_oxy ~ add(addSd_oxy) + prop(propSd_oxy)
  })
}
