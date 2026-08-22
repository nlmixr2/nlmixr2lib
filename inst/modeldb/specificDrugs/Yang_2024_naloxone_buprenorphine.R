Yang_2024_naloxone_buprenorphine <- function() {
  description <- paste(
    "QSP. Mechanistic PK-PD model of buprenorphine-induced respiratory",
    "depression and its reversal by intramuscular naloxone auto-injector",
    "(NAI). Couples a three-compartment buprenorphine PK model (the",
    "parent, unsuffixed) to the Yang 2024 naloxone PK model (depot plus",
    "three transit absorption compartments and a two-compartment",
    "disposition, all _naloxone-suffixed). Both drugs equilibrate with",
    "their own biophase (effect) compartment and then compete for the",
    "mu-opioid (MOP) receptor: a single receptor-occupancy state RL_op",
    "carries the apparent fractional buprenorphine-receptor occupancy,",
    "with naloxone entering through the quasi-steady-state competitive",
    "term (1 - Ce_naloxone / (KD_naloxone + Ce_naloxone)). Minute",
    "ventilation follows a linear transduction function",
    "VE = V0 * (1 - alpha * RL_op). Buprenorphine's very slow receptor",
    "dissociation (koff 0.0172 /min, half-life ~40 min) is what makes it",
    "the hardest of the four Yang 2024 opioids for naloxone to reverse.",
    "Buprenorphine PK-PD parameters are inherited unchanged from Yassen",
    "2007; the naloxone PK layer is Yang 2024's own fit."
  )
  reference <- paste(
    "Yang TE, Del Bene F, Lavezzi SM, Iavarone L, Zhang J, Kim J,",
    "Gjurich B, Kessler C. Mechanistic pharmacokinetic-pharmacodynamic",
    "modeling and simulations of naloxone auto-injector 10 mg reversal",
    "of opioid-induced respiratory depression.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1722-1733.",
    "doi:10.1002/psp4.13215. PMCID PMC11494827.",
    "Parameter source: Appendix S1 Table S5 (Buprenorphine Population",
    "PK-PD Parameter Estimates) and Appendix S1 example NONMEM control",
    "stream (A) Buprenorphine, whose $THETA / $OMEGA blocks and $DES",
    "equations this file reproduces one-for-one.",
    "Buprenorphine PK-PD parameters originate from Yassen A, Olofsen E,",
    "van Dorp E, Sarton E, Teppema L, Danhof M, Dahan A.",
    "Mechanism-based pharmacokinetic-pharmacodynamic modelling of the",
    "reversal of buprenorphine-induced respiratory depression by",
    "naloxone: a study in healthy volunteers.",
    "Clin Pharmacokinet. 2007;46(11):965-980.",
    "doi:10.2165/00003088-200746110-00004.",
    "Naloxone PK parameters are Yang 2024 Appendix S1 Table S4."
  )
  vignette <- "Yang_2024_naloxone_opioid_reversal"
  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ng/mL"
  )
  # The absorption-chain heuristic in buildModelDb() recognises only
  # `depot` and `central`; naloxone is dosed into the _naloxone-suffixed
  # depot, so both dosing compartments are declared explicitly.
  dosing <- c("central", "depot_naloxone")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Appendix S1 control stream
  # (A) $MODEL block, which names every compartment, and its unit header
  # ("; AMT in ug ... ; PK volumes = L ... ventilation = L/min").
  compartmentData <- list(
    central              = list(analyte = "buprenorphine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "buprenorphine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2          = list(analyte = "buprenorphine", units = "ug", specimen = "plasma", verified = TRUE),
    effect               = list(analyte = "buprenorphine biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    depot_naloxone       = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit1_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit2_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit3_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    central_naloxone     = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_naloxone = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    effect_naloxone      = list(analyte = "naloxone biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    RL_op                = list(analyte = "apparent fractional buprenorphine occupancy of the mu-opioid receptor pool", units = "fraction", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Acts only on naloxone apparent clearance:",
        "CL/F = 3.26 * (WT/70)^0.538 (Yang 2024 Table S4). The reference",
        "weight of 70 kg is explicit in control stream (A)",
        "('WT= 70 ; kg'). Yang 2024's own simulations fixed the virtual",
        "population at WT = 70 kg because the weight effect on naloxone",
        "AUC was judged not clinically significant; WT is kept exposed",
        "here so subject-level weight can vary. The buprenorphine PK",
        "parameters carry no weight scaling in Yassen 2007."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Healthy adults (naloxone layer 23-54 years)",
    weight_range   = "Simulations fix WT = 70 kg; naloxone PK data spanned 57.2-100.2 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy adult volunteers. The buprenorphine PK-PD layer comes",
      "from the Yassen 2007 healthy-volunteer study of",
      "buprenorphine-induced respiratory depression reversed by a",
      "naloxone infusion; the naloxone PK layer is Yang 2024's own",
      "auto-injector population PK in 48 healthy adults."
    ),
    dose_range     = paste(
      "Yang 2024 simulated IV buprenorphine at 0.9, 9.9 and",
      "18.9 ug/kg (the lowest corresponding to 50 percent ventilation",
      "suppression, the highest 21-fold higher), each with naloxone",
      "auto-injector 2 mg or 10 mg IM/SC. Naloxone was given when",
      "ventilation had fallen 30 percent from baseline (a lower",
      "trigger than the 60 percent used for the other three opioids,",
      "because buprenorphine's respiratory depressant effect has a",
      "slow onset)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Yang 2024 validated this constructed model by digitising the",
      "ventilation-time profiles for buprenorphine 0.2 and 0.4 mg/70 kg",
      "from Figure 3 of Yassen 2007 and confirming that most observed",
      "points fell inside the 90 percent prediction interval of 1000",
      "simulations (Appendix S1 Figure S7).",
      "Baseline ventilation V0 was carried as a per-subject data column",
      "(E0) in the control stream rather than as a THETA; it is encoded",
      "here as the canonical le0 parameter with the Table S5 typical",
      "value and IIV so that the model is self-contained.",
      "Buprenorphine is a partial agonist: alpha = 0.67 caps the",
      "attainable ventilation suppression at 67 percent even at full",
      "receptor occupancy."
    )
  )

  ini({
    # Every value below is held FIX in Appendix S1 control stream (A):
    # the buprenorphine PK-PD block is inherited unchanged from Yassen
    # 2007 and the naloxone PK block from Yang 2024's own Table S4, and
    # neither was re-fitted for this analysis. They are therefore
    # wrapped in fixed().

    # --- Naloxone PK (Yang 2024 Table S4; control stream (A) THETA 1-6)
    lktr_naloxone <- fixed(log(0.696))
    label("Naloxone transit absorption rate constant (KTR, 1/min)")            # (A) THETA 1 = 0.696; Table S4
    lcl_naloxone <- fixed(log(3.26))
    label("Naloxone apparent clearance at 70 kg (CL/F, L/min)")                # (A) THETA 2 = 3.26; Table S4
    e_wt_cl_naloxone <- fixed(0.538)
    label("Allometric exponent of body weight on naloxone CL/F (unitless)")    # (A) THETA 3 = 0.538; Table S4
    lvc_naloxone <- fixed(log(404))
    label("Naloxone apparent central volume (V2/F, L)")                        # (A) THETA 4 = 404; Table S4
    lq_naloxone <- fixed(log(0.0847))
    label("Naloxone apparent intercompartmental clearance (Q/F, L/min)")       # (A) THETA 5 = 0.0847; Table S4
    lvp_naloxone <- fixed(log(81.8))
    label("Naloxone apparent peripheral volume (V3/F, L)")                     # (A) THETA 6 = 81.8; Table S4

    # --- Buprenorphine PK (Table S5; control stream (A) THETA 7-12).
    # NONMEM V1/V2/V3 map to the canonical vc/vp/vp2 and Q2/Q3 to q/q2.
    lcl <- fixed(log(1.43))
    label("Buprenorphine clearance (CL, L/min)")                               # (A) THETA 7 = 1.43; Table S5
    lvc <- fixed(log(5.71))
    label("Buprenorphine central volume (V1, L)")                              # (A) THETA 8 = 5.71; Table S5
    lq <- fixed(log(1.45))
    label("Buprenorphine intercompartmental clearance to peripheral1 (Q2, L/min)")  # (A) THETA 9 = 1.45; Table S5
    lvp <- fixed(log(22.1))
    label("Buprenorphine first peripheral volume (V2, L)")                     # (A) THETA 10 = 22.1; Table S5
    lq2 <- fixed(log(1.00))
    label("Buprenorphine intercompartmental clearance to peripheral2 (Q3, L/min)")  # (A) THETA 11 = 1.00; Table S5
    lvp2 <- fixed(log(135))
    label("Buprenorphine second peripheral volume (V3, L)")                    # (A) THETA 12 = 135; Table S5

    # --- Biophase equilibration (equation 1 of Appendix S1)
    lke0_naloxone <- fixed(log(0.106))
    label("Naloxone biophase equilibration rate constant (ke0, 1/min)")        # (A) THETA 13 = 0.106; Table S5
    lke0 <- fixed(log(0.00404))
    label("Buprenorphine biophase equilibration rate constant (ke0, 1/min)")   # (A) THETA 14 = 0.00404; Table S5

    # --- MOP receptor association / dissociation (equations 4-6)
    lkon <- fixed(log(0.203))
    label("Buprenorphine receptor association rate constant (kon, mL/ng/min)") # (A) THETA 15 = 0.203; Table S5
    lkoff <- fixed(log(0.0172))
    label("Buprenorphine receptor dissociation rate constant (koff, 1/min)")   # (A) THETA 16 = 0.0172; Table S5
    lec50_naloxone <- fixed(log(1.87))
    label("Naloxone equilibrium dissociation constant (KD = koff/kon, ng/mL)") # (A) THETA 17 = 1.87; Table S5

    # --- Linear transduction function (equation 7)
    lalpha <- fixed(log(0.67))
    label("Buprenorphine intrinsic activity (alpha, unitless 0-1)")            # (A) THETA 18 = 0.67; Table S5
    le0 <- fixed(log(24.0))
    label("Baseline minute ventilation (V0, L/min)")                           # Table S5: V0 24.0 L/min (carried as the E0 data column in control stream (A))

    # --- IIV. Control stream (A) $OMEGA entries 1-13, in order. Table S5
    # prints these as "IIV (%)"; each tabulated percentage is the square
    # root of the corresponding $OMEGA variance (e.g. buprenorphine CL
    # 19.1% -> 0.036481 = 0.191^2), so the two sources agree exactly.
    etalktr_naloxone ~ fixed(0.111)
    etalcl_naloxone ~ fixed(0.0129)
    etalvc_naloxone ~ fixed(0.0658)
    etalcl ~ fixed(0.036481)
    etalvc ~ fixed(0.033856)
    etalq ~ fixed(0.133225)
    etalq2 ~ fixed(0.078961)
    etalvp2 ~ fixed(0.123904)
    etalke0 ~ fixed(0.697225)
    etalkon ~ fixed(0.744769)
    etalkoff ~ fixed(0.931225)
    etalec50_naloxone ~ fixed(0.799236)
    etalalpha ~ fixed(0.061504)
    # Baseline ventilation IIV. Control stream (A) has no $OMEGA entry
    # for it because E0 arrived as a per-subject data column already
    # carrying the variability; Table S5 reports V0 IIV = 21.8%, giving
    # the variance 0.218^2 = 0.047524 on the same convention as every
    # other row of that table.
    etale0 ~ fixed(0.047524)

    # Buprenorphine plasma-concentration residual error, Table S5
    # "Proportional error (%) 13.5". Control stream (A) does not fit the
    # PK endpoint (it simulates ventilation only, with the PK held
    # fixed), so this term is carried here for completeness of the
    # PK endpoint rather than reproduced from that stream.
    propSd <- fixed(0.135)
    label("Buprenorphine proportional residual error SD (fraction)")           # Table S5: proportional error 13.5%
    # Ventilation residual error, control stream (A) THETA 19 with
    # $SIGMA 1 FIX, so THETA(19) is the additive SD directly.
    addSd <- fixed(2.24)
    label("Ventilation additive residual error SD (L/min)")                    # (A) THETA 19 = 2.24; Table S5
  })

  model({
    # ---- Naloxone individual parameters
    ktr_naloxone <- exp(lktr_naloxone + etalktr_naloxone)
    cl_naloxone <- exp(lcl_naloxone + etalcl_naloxone) * (WT / 70)^e_wt_cl_naloxone
    vc_naloxone <- exp(lvc_naloxone + etalvc_naloxone)
    q_naloxone <- exp(lq_naloxone)
    vp_naloxone <- exp(lvp_naloxone)
    ke0_naloxone <- exp(lke0_naloxone)
    ec50_naloxone <- exp(lec50_naloxone + etalec50_naloxone)

    # ---- Buprenorphine individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)
    ke0 <- exp(lke0 + etalke0)
    kon <- exp(lkon + etalkon)
    koff <- exp(lkoff + etalkoff)
    alpha <- exp(lalpha + etalalpha)
    e0 <- exp(le0 + etale0)

    # ---- Naloxone PK: depot + three transit compartments in series
    # (all governed by KTR) feeding a two-compartment disposition.
    # Control stream (A) $DES DADT(1)-DADT(6).
    d/dt(depot_naloxone) <- -ktr_naloxone * depot_naloxone
    d/dt(transit1_naloxone) <- ktr_naloxone * depot_naloxone - ktr_naloxone * transit1_naloxone
    d/dt(transit2_naloxone) <- ktr_naloxone * transit1_naloxone - ktr_naloxone * transit2_naloxone
    d/dt(transit3_naloxone) <- ktr_naloxone * transit2_naloxone - ktr_naloxone * transit3_naloxone
    d/dt(central_naloxone) <- ktr_naloxone * transit3_naloxone -
      (cl_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vc_naloxone) * central_naloxone +
      (q_naloxone / vp_naloxone) * peripheral1_naloxone
    d/dt(peripheral1_naloxone) <- (q_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vp_naloxone) * peripheral1_naloxone

    # Naloxone biophase equilibration (equation 1). Amounts are ug and
    # volumes L, so central_naloxone/vc_naloxone is ug/L = ng/mL, the
    # scale on which KD_naloxone (1.87 ng/mL) is expressed.
    d/dt(effect_naloxone) <- ke0_naloxone * (central_naloxone / vc_naloxone - effect_naloxone)

    # ---- Buprenorphine PK: three compartments, control stream (A)
    # $DES DADT(8)-DADT(10) written with micro-constants K80 = CL/V1,
    # K89 = Q2/V1, K98 = Q2/V2, K810 = Q3/V1, K108 = Q3/V3.
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(peripheral2) <- (q2 / vc) * central - (q2 / vp2) * peripheral2

    # Buprenorphine biophase equilibration (equation 1).
    d/dt(effect) <- ke0 * (central / vc - effect)

    # ---- MOP receptor kinetics with competitive naloxone inhibition.
    # Appendix S1 equation 6 / control stream (A) DADT(12):
    #   dRL/dt = kon*Ce_op*(1 - RL)*(1 - Ce_nal/(KD + Ce_nal)) - koff*RL
    # The naloxone term is the quasi-steady-state occupancy of
    # equation 5, valid because naloxone unbinds far faster than
    # buprenorphine. RL_op is a fraction, so it is dimensionless and
    # takes no dose.
    d/dt(RL_op) <- kon * effect * (1 - RL_op) *
      (1 - effect_naloxone / (ec50_naloxone + effect_naloxone)) -
      koff * RL_op

    # ---- Outputs
    Cc <- central / vc
    Cnal <- central_naloxone / vc_naloxone
    # Linear transduction function, equation 7 / control stream (A)
    # $ERROR E = E0*(1 - ALPHA*A(12)).
    VE <- e0 * (1 - alpha * RL_op)
    # Ventilation as a fraction of each subject's own baseline, the
    # quantity Yang 2024 reports its 40 / 70 / 85 percent recovery
    # thresholds against (control stream ERATIO = IPRED/E0).
    ERATIO <- VE / e0

    Cc ~ prop(propSd)
    VE ~ add(addSd)
  })
}
