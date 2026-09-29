Bulitta_2019_pefloxacin <- function() {
  description <- "Three-compartment population PK model for pefloxacin and its two urinary metabolites (norfloxacin, pefloxacin N-oxide) in 8 adult cystic-fibrosis patients and 10 healthy volunteers, each given 400 mg pefloxacin as a 30 min IV infusion and 400 mg orally in a randomised two-way crossover (Bulitta 2019). Disposition is central plus one peripheral compartment plus an intestinal recirculation pool: pefloxacin is exsorbed from the central compartment into the gut lumen by a saturable clearance CLEX = CLGUT * KmEX / (KmEX + C1), with CLGUT fixed to gut blood flow (66 L/h) and all of the exsorbed drug subsequently reabsorbed at a first-order rate, so the pool acts as a distribution site rather than an elimination route. Parent elimination is split into renal and non-renal arms; the non-renal arm is partitioned by logit-transformed formation fractions into norfloxacin and pefloxacin N-oxide, each tracked as a cumulative urinary amount. Fat-free mass is the size descriptor with allometric scaling (exponents 0.75 on all clearances, 1.0 on all volumes; FFM_STD = 53 kg). Healthy volunteers are the reference group: cystic fibrosis multiplicatively scales renal clearance (1.53), non-renal clearance (0.861) and both volumes (0.916) via fcyf_*^DIS_CF, and separately estimated typical values apply to bioavailability, absorption half-life and reabsorption half-life."
  reference <- paste(
    "Bulitta JB, Jiao Y, Landersdorfer CB, Sutaria DS, Tao X, Shin E,",
    "Hohl R, Holzgrabe U, Stephan U, Sorgel F. Comparable bioavailability",
    "and disposition of pefloxacin in patients with cystic fibrosis and",
    "healthy volunteers assessed via population pharmacokinetics.",
    "Pharmaceutics. 2019;11(7):323. doi:10.3390/pharmaceutics11070323."
  )
  vignette <- "Bulitta_2019_pefloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "pefloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pefloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pefloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    gut_lumen = list(analyte = "pefloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "pefloxacin", units = "mg", specimen = "urine", verified = TRUE),
    urine_norflox = list(analyte = "norfloxacin", units = "mg", specimen = "urine", verified = TRUE),
    urine_noxpeflox = list(analyte = "pefloxacin N-oxide", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. Computed by the Janmahasatian formula",
        "(Bulitta 2019 Table 1 footnote a, citing ref 26). Allometric scaling",
        "on all clearance terms (exponent 0.75 fixed) and on all volume terms",
        "(exponent 1.0 fixed) with reference FFM_STD = 53 kg (Bulitta 2019",
        "Methods 'Body size and composition', eqs. 9 and 10, and the",
        "supplement Figure S3 lines 75-76: FWTCL = (FFM/53)**0.75,",
        "FWTV = (FFM/53)). Table 1 medians [ranges]: 33.3 [27.3-46.4] kg in",
        "patients with CF and 52.1 [37.7-64.0] kg in healthy volunteers."
      ),
      source_name = "FFM"
    ),
    DIS_CF = list(
      description = "Cystic-fibrosis disease-state indicator (1 = CF patient, 0 = healthy volunteer)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-fixed per subject. Corresponds to the GRP flag in the supplement",
        "Figure S3 estimation code (GRP == 1 selects the CF branch). Acts on the",
        "model in two distinct ways, exactly as the source code does. (1) Power-form",
        "disease scale factors fcyf_clr^DIS_CF, fcyf_clnr^DIS_CF and fcyf_vss^DIS_CF",
        "multiply renal clearance, non-renal clearance and both volumes respectively",
        "(Bulitta 2019 Table 4, FFM-allometric row; Figure S3 lines 41-44, 54-57,",
        "79-84). (2) Bioavailability, absorption rate, reabsorption rate and both",
        "metabolite formation fractions carry separately estimated typical values",
        "and separately estimated between-subject variances for the two groups",
        "(Figure S3 lines 37-39, 46-47, 50-52, 59-60); healthy volunteers keep the",
        "bare parameter name and CF patients take the _cf suffix."
      ),
      source_name = "GRP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "CF patients 19 [17-24] years (median [range]); healthy volunteers 24 [18-27] years",
    weight_range = "CF patients 46.3 [35.5-63.5] kg total body weight; healthy volunteers 77.5 [55.0-82.0] kg",
    ffm_range = "CF patients 33.3 [27.3-46.4] kg; healthy volunteers 52.1 [37.7-64.0] kg",
    lbm_range = "CF patients 38.0 [31.0-47.1] kg; healthy volunteers 55.6 [43.0-64.5] kg",
    height_range = "CF patients 166 [158-175] cm; healthy volunteers 174 [168-191] cm",
    bmi_range = "CF patients 17.6 [13.4-22.2] kg/m^2; healthy volunteers 21.6 [19.5-27.3] kg/m^2",
    sex_female_pct = 61.1,
    race_ethnicity = "Caucasian (all 18 subjects)",
    disease_state = "8 adult cystic-fibrosis patients (one aged 17 years, enrolled with consent from his legal representative) and 10 healthy volunteers",
    dose_range = "400 mg pefloxacin as a 30 min IV infusion and 400 mg orally, randomised two-way crossover with a 10-day washout",
    regions = "Germany (single centre; University Hospital Essen ethics approval, 1984)",
    notes = paste(
      "Demographics from Bulitta 2019 Table 1. Plasma sampled to 48 h after each",
      "dose; urine collected over eleven intervals to 48 h, but only the cumulative",
      "amount excreted to the last interval was available for modelling, which is",
      "why the urinary states are observed as cumulative amounts. Pefloxacin,",
      "norfloxacin and pefloxacin N-oxide were quantified by reversed-phase HPLC",
      "with fluorescence detection (linear range 0.078-20 mg/L for pefloxacin in",
      "plasma). Plasma concentrations of the two metabolites were NOT measured, so",
      "their clearances and volumes are not identifiable and are not part of this",
      "model; only their formation fractions are estimated. All population modelling",
      "used importance sampling (pmethod = 4) in S-ADAPT 1.57 via SADAPT-TRAN.",
      "Structural parameters from Table 3 and supplementary Table S1 (original",
      "dataset column); disease scale factors from Table 4, FFM-allometric row;",
      "model code from supplementary Figure S3."
    )
  )

  ini({
    # ========================================================================
    # Reference subject: HEALTHY VOLUNTEER with FFM = FFM_STD = 53 kg.
    # The 'Patients with CF' columns of Table 3 are reproduced from these
    # healthy-volunteer values by the disease factors and _cf parameters
    # below; see the vignette source-trace table for the arithmetic.
    # ========================================================================

    # ---- Absorption (Bulitta 2019 Table 3) ---------------------------------
    # The source estimates absorption and reabsorption as HALF-LIVES in
    # minutes; Figure S3 lines 38-39 and 51-52 convert them to first-order
    # rate constants per hour as KA = LOG(2)/(Tabs/60). The log-scale rate
    # constants below are that conversion applied to the published half-lives,
    # so the eta variances (which are on a log scale) carry over unchanged.
    lfdepot <- log(1.03)
    label("Oral bioavailability F_BIO, healthy volunteers (fraction)") # Table 3, row 'Oral bioavailability', healthy column: 1.03 (SE 4.70%)
    lfdepot_cf <- log(1.00)
    label("Oral bioavailability F_BIO, patients with CF (fraction)") # Table 3, row 'Oral bioavailability', CF column: 1.00 (SE 5.60%)

    ltlag <- log(13.3 / 60)
    label("Oral absorption lag time T_lag, both groups (h)") # Table 3, row 'Absorption lag-time': 13.3 min (SE 4.00%), same estimate in both groups

    lka <- log(3.713288)
    label("Oral absorption rate constant, healthy volunteers (1/h)") # Table 3, row 'Absorption half-life', healthy: T_abs = 11.2 min -> ka = log(2)/(11.2/60) = 3.7133 /h
    lka_cf <- log(2.100446)
    label("Oral absorption rate constant, patients with CF (1/h)") # Table 3, row 'Absorption half-life', CF: T_abs = 19.8 min -> ka = log(2)/(19.8/60) = 2.1004 /h

    lkr <- log(1.999463)
    label("Reabsorption rate constant from the gut lumen, healthy volunteers (1/h)") # Table 3, row 'Reabsorption half-life from intestine', healthy: T_reabs = 20.8 min -> kr = log(2)/(20.8/60) = 1.9995 /h
    lkr_cf <- log(0.6359148)
    label("Reabsorption rate constant from the gut lumen, patients with CF (1/h)") # Table 3, row 'Reabsorption half-life from intestine', CF: T_reabs = 65.4 min -> kr = log(2)/(65.4/60) = 0.6359 /h

    # ---- Disposition (Bulitta 2019 Table 3, healthy volunteers at FFM 53 kg)
    lvc <- log(40.8)
    label("Central volume of distribution V1 (L)") # Table 3, row 'Volume of distribution for central compartment', healthy: 40.8 L (SE 12.9%)
    lvp <- log(65.4)
    label("Peripheral volume of distribution V2 (L)") # Table 3, row 'Volume of distribution for peripheral compartment', healthy: 65.4 L (SE 5%)
    lcl_nonren <- log(8.56)
    label("Non-renal clearance CL_NR (L/h)") # Table 3, row 'Non-renal clearance', healthy: 8.56 L/h (SE 7.70%)
    lcl_renal <- log(0.705)
    label("Renal clearance CL_R (L/h)") # Table 3, row 'Renal clearance', healthy: 0.705 L/h (SE 6.50%)
    lq <- log(406)
    label("Distribution clearance CL_D (L/h)") # Table 3, row 'Distribution clearance': 406 L/h (SE 33.7%), same estimate in both groups

    # ---- Saturable intestinal exsorption (Bulitta 2019 eq. 4) --------------
    # CLEX = VmaxEX / (KmEX + C1) with VmaxEX = CLGUT * KmEX, i.e.
    # CLEX = CLGUT * KmEX / (KmEX + C1) as coded in Figure S3 line 21.
    lclgut <- fixed(log(66))
    label("Maximum exsorption clearance CL_GUT, set to blood flow to gut (L/h)") # Table 3, row 'Gut clearance for enterohepatic circulation': 66 L/h (fixed); Results 'Population Pharmacokinetic Modeling' explains the fixing to gut blood flow
    lkm <- log(1.44)
    label("Plasma concentration at half-maximal exsorption clearance Km_EX (mg/L)") # Table 3, row 'Plasma concentration associated with half-maximal CL GUT': 1.44 mg/L (SE 12.0%)

    # ---- Metabolite formation fractions (Bulitta 2019 Table 3) -------------
    # The source holds these on a logit scale (Figure S3 lines 65-71) so the
    # fractions stay bounded; Table 3 footnote d reports the median of the
    # individual subject estimates for each group. The logit values below are
    # log(fm / (1 - fm)) of those medians.
    logitfm_norflox <- -1.557539
    label("Logit of the non-renal clearance fraction forming norfloxacin, healthy volunteers (logit scale, unitless)") # Table 3, row 'Norfloxacin / Formation fraction', healthy: fm_NOR = 0.174; logit(0.174) = -1.5575
    logitfm_norflox_cf <- -1.380056
    label("Logit of the non-renal clearance fraction forming norfloxacin, patients with CF (logit scale, unitless)") # Table 3, row 'Norfloxacin / Formation fraction', CF: fm_NOR = 0.201; logit(0.201) = -1.3801
    logitfm_noxpeflox <- -1.529957
    label("Logit of the non-renal clearance fraction forming pefloxacin N-oxide, healthy volunteers (logit scale, unitless)") # Table 3, row 'Pefloxacin N-oxide / Formation fraction', healthy: fm_NOX = 0.178; logit(0.178) = -1.5300
    logitfm_noxpeflox_cf <- -1.130873
    label("Logit of the non-renal clearance fraction forming pefloxacin N-oxide, patients with CF (logit scale, unitless)") # Table 3, row 'Pefloxacin N-oxide / Formation fraction', CF: fm_NOX = 0.244; logit(0.244) = -1.1309

    # ---- Disease-state scale factors (Bulitta 2019 Table 4, FFM-allometric row)
    fcyf_clr <- 1.53
    label("CF / healthy ratio on renal clearance (unitless)") # Table 4, 'FFM allometric' row, F_CYF,CLR = 1.53 (SE 12.5%); Table S1 bootstrap median 1.53 (95% CI 1.21-1.78)
    fcyf_clnr <- 0.861
    label("CF / healthy ratio on non-renal clearance (unitless)") # Table 4, 'FFM allometric' row, F_CYF,CLNR = 0.861 (SE 12.3%); Table S1 bootstrap median 0.850 (95% CI 0.694-1.07)
    fcyf_vss <- 0.916
    label("CF / healthy ratio on both volumes of distribution (unitless)") # Table 4, 'FFM allometric' row, F_CYF,VSS = 0.916 (SE 14.7%); Table S1 bootstrap median 0.931 (95% CI 0.739-1.15)

    # ---- Allometric scaling by fat-free mass (Bulitta 2019 eqs. 9 and 10) --
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on every clearance term (unitless)") # Methods 'Body size and composition': 'We fixed the allometric exponent to 1.0 for all volumes and to 0.75 for all clearances'; Figure S3 line 75
    e_ffm_vc <- fixed(1.00)
    label("Allometric exponent on the central volume (unitless)") # Methods 'Body size and composition'; Figure S3 line 76
    e_ffm_vp <- fixed(1.00)
    label("Allometric exponent on the peripheral volume (unitless)") # Methods 'Body size and composition'; Figure S3 line 76
    ffm_std <- fixed(53)
    label("Reference fat-free mass FFM_STD (kg)") # Methods 'Body size and composition': 'a standard fat-free mass FFMSTD of 53 kg'; Figure S3 lines 75-76 divide by 53

    # ========================================================================
    # Between-subject variability (Bulitta 2019 Table 3 / Table S1).
    # Table 3 footnote a: the tabulated BSV values are 'apparent coefficients
    # of variation of a normal distribution on natural logarithmic scale',
    # defined in Methods 'Between-subject variability model' as 'the square
    # roots of the estimated variances'. The omega entries below are therefore
    # the SQUARE of the tabulated value, NOT log(1 + CV^2). The distinction is
    # immaterial for the small BSVs but decisive for CL_D (1.19^2 = 1.4161 vs
    # log(1 + 1.19^2) = 0.8819).
    #
    # V1, V2, CL_NR, CL_R and CL_D carry one shared variance across both
    # groups (a single line in Table S1); F_BIO, T_abs and T_reabs carry
    # separate per-group variances (two lines, footnotes a and b). The
    # group-specific etas below are selected by DIS_CF inside model(), so
    # exactly one of each pair contributes to any given subject.
    # ========================================================================
    etalfdepot ~ 0.016129 # Table 3, 'Oral bioavailability' healthy BSV 0.127 -> 0.127^2
    etalfdepot_cf ~ 0.021316 # Table 3, 'Oral bioavailability' CF BSV 0.146 -> 0.146^2
    etaltlag ~ 0.00506944 # Table 3, 'Absorption lag-time' BSV 0.0712 -> 0.0712^2 (same in both groups)
    etalka ~ 0.966289 # Table 3, 'Absorption half-life' healthy BSV 0.983 -> 0.983^2
    etalka_cf ~ 1.21 # Table 3, 'Absorption half-life' CF BSV 1.10 -> 1.10^2
    etalkr ~ 0.646416 # Table 3, 'Reabsorption half-life from intestine' healthy BSV 0.804 -> 0.804^2
    etalkr_cf ~ 0.111556 # Table 3, 'Reabsorption half-life from intestine' CF BSV 0.334 -> 0.334^2
    etalvc ~ 0.189225 # Table 3, 'Volume of distribution for central compartment' BSV 0.435 -> 0.435^2 (shared)
    etalvp ~ 0.017161 # Table 3, 'Volume of distribution for peripheral compartment' BSV 0.131 -> 0.131^2 (shared)
    etalcl_nonren ~ 0.056644 # Table 3, 'Non-renal clearance' BSV 0.238 -> 0.238^2 (shared)
    etalcl_renal ~ 0.028224 # Table 3, 'Renal clearance' BSV 0.168 -> 0.168^2 (shared)
    etalq ~ 1.4161 # Table 3, 'Distribution clearance' BSV 1.19 -> 1.19^2 (shared)
    etalkm ~ fixed(0.01) # Table 3, 'Plasma concentration associated with half-maximal CL GUT', BSV of 0.1 held constant by the authors -> 0.1^2

    # ---- Residual unexplained variability (Bulitta 2019 Table 3 footnote) --
    addSd <- 0.00984
    label("Additive residual error on plasma pefloxacin concentrations (mg/L)") # Table 3 footnote: 'The additive and proportional residual errors of plasma concentrations were 0.00984 mg/L and 15.1% for pefloxacin'
    propSd <- 0.151
    label("Proportional residual error on plasma pefloxacin concentrations (fraction)") # Table 3 footnote: proportional residual error 15.1%
    # The source fixed the additive residual error of each urinary output to
    # 1% OF THE DOSE (Methods 'Residual error model and uncertainty'; Figure S3
    # lines 102-104 output 100 * X / 400, i.e. percent of the 400 mg dose).
    # The outputs here are cumulative AMOUNTS in mg, so the equivalent SD for
    # this study's fixed 400 mg dose is 0.01 * 400 = 4 mg.
    addSd_Aurine <- fixed(4)
    label("Additive residual error on cumulative urinary pefloxacin (mg)") # Methods 'Residual error model and uncertainty': additive error fixed to 1% of dose; 0.01 * 400 mg = 4 mg
    addSd_Aurine_norflox <- fixed(4)
    label("Additive residual error on cumulative urinary norfloxacin (mg)") # Methods 'Residual error model and uncertainty': additive error fixed to 1% of dose; 0.01 * 400 mg = 4 mg
    addSd_Aurine_noxpeflox <- fixed(4)
    label("Additive residual error on cumulative urinary pefloxacin N-oxide (mg)") # Methods 'Residual error model and uncertainty': additive error fixed to 1% of dose; 0.01 * 400 mg = 4 mg
  })

  model({
    # ---- Allometric size factors (Bulitta 2019 eqs. 9-10; Figure S3 75-76) --
    size_cl <- (FFM / ffm_std)^e_ffm_cl
    size_vc <- (FFM / ffm_std)^e_ffm_vc
    size_vp <- (FFM / ffm_std)^e_ffm_vp

    # ---- Group-selected parameters -----------------------------------------
    # Figure S3 lines 35-61 branch on GRP: the CF branch uses its own THETA
    # and its own eta, the healthy branch uses the reference pair. Writing the
    # selection as a DIS_CF-weighted sum of the two (fixed effect + eta) terms
    # reproduces that branch exactly, including the group-specific variances.
    fdepot <- exp((lfdepot + etalfdepot) * (1 - DIS_CF) +
      (lfdepot_cf + etalfdepot_cf) * DIS_CF)
    ka <- exp((lka + etalka) * (1 - DIS_CF) +
      (lka_cf + etalka_cf) * DIS_CF)
    kr <- exp((lkr + etalkr) * (1 - DIS_CF) +
      (lkr_cf + etalkr_cf) * DIS_CF)
    tlag <- exp(ltlag + etaltlag)

    # ---- Individual disposition parameters (Figure S3 lines 79-84) ---------
    cl_renal <- exp(lcl_renal + etalcl_renal) * size_cl * fcyf_clr^DIS_CF
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * size_cl * fcyf_clnr^DIS_CF
    q <- exp(lq + etalq) * size_cl
    clgut <- exp(lclgut) * size_cl
    vc <- exp(lvc + etalvc) * size_vc * fcyf_vss^DIS_CF
    vp <- exp(lvp + etalvp) * size_vp * fcyf_vss^DIS_CF
    km <- exp(lkm + etalkm)

    # ---- Metabolite formation fractions (Figure S3 lines 65-71) ------------
    # Shares of the NON-RENAL clearance; the complement (1 - fm_norflox -
    # fm_noxpeflox) leaves the system as unmeasured non-renal elimination.
    logitfm_norflox_ind <- logitfm_norflox * (1 - DIS_CF) + logitfm_norflox_cf * DIS_CF
    logitfm_noxpeflox_ind <- logitfm_noxpeflox * (1 - DIS_CF) + logitfm_noxpeflox_cf * DIS_CF
    fm_norflox <- expit(logitfm_norflox_ind)
    fm_noxpeflox <- expit(logitfm_noxpeflox_ind)

    # ---- Concentrations ----------------------------------------------------
    Cc <- central / vc
    Cp <- peripheral1 / vp

    # ---- Saturable exsorption clearance into the gut lumen (eq. 4) ---------
    cl_exs <- clgut * km / (km + Cc)

    # ---- ODE system (Bulitta 2019 eqs. 1-3, 5-8; Figure S3 lines 23-31) ----
    # The intravenous infusion enters `central` directly (the R(1) term of
    # Figure S3 line 24) and is supplied through the event table rather than
    # as a structural parameter. All exsorbed pefloxacin is reabsorbed, so
    # `gut_lumen` is a distribution site with no elimination arm of its own.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (cl_renal + cl_nonren + cl_exs) * Cc -
      q * (Cc - Cp) +
      kr * gut_lumen
    d/dt(peripheral1) <- q * (Cc - Cp)
    d/dt(gut_lumen) <- cl_exs * Cc - kr * gut_lumen
    d/dt(urine) <- cl_renal * Cc
    d/dt(urine_norflox) <- fm_norflox * cl_nonren * Cc
    d/dt(urine_noxpeflox) <- fm_noxpeflox * cl_nonren * Cc

    # Oral bioavailability and lag time act on the depot (Figure S3 line 63
    # BOLUSF(1) = FBIO, and lines 15-19 which withhold absorption until the
    # lag time has elapsed).
    f(depot) <- fdepot
    alag(depot) <- tlag

    # ---- Observations ------------------------------------------------------
    # Cumulative urinary amounts in mg; the source reports them as a percent
    # of the 400 mg dose (Figure S3 lines 102-104).
    Aurine <- urine
    Aurine_norflox <- urine_norflox
    Aurine_noxpeflox <- urine_noxpeflox

    Cc ~ add(addSd) + prop(propSd)
    Aurine ~ add(addSd_Aurine)
    Aurine_norflox ~ add(addSd_Aurine_norflox)
    Aurine_noxpeflox ~ add(addSd_Aurine_noxpeflox)
  })
}
