Yang_2024_naloxone_morphine <- function() {
  description <- paste(
    "QSP. Mechanistic PK-PD model of morphine-induced respiratory",
    "depression and its reversal by intramuscular naloxone auto-injector",
    "(NAI). Couples a three-compartment morphine PK model parameterised",
    "directly in micro-constants (the parent, unsuffixed) to the",
    "Yang 2024 naloxone PK model (depot plus three transit absorption",
    "compartments and a two-compartment disposition, all",
    "_naloxone-suffixed). Both drugs equilibrate with their own biophase",
    "(effect) compartment; a single receptor-occupancy state RL_op",
    "carries the apparent fractional morphine-receptor occupancy, and",
    "naloxone displaces morphine through a steep sigmoidal competitive",
    "term with its own shape parameter gamma. Minute ventilation follows",
    "a linear transduction function VE = V0 * (1 - alpha * RL_op) with",
    "alpha fixed at 1, so full receptor occupancy abolishes ventilation.",
    "Morphine PK-PD parameters are inherited unchanged from Olofsen",
    "2010; the naloxone PK layer is Yang 2024's own fit."
  )
  reference <- paste(
    "Yang TE, Del Bene F, Lavezzi SM, Iavarone L, Zhang J, Kim J,",
    "Gjurich B, Kessler C. Mechanistic pharmacokinetic-pharmacodynamic",
    "modeling and simulations of naloxone auto-injector 10 mg reversal",
    "of opioid-induced respiratory depression.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1722-1733.",
    "doi:10.1002/psp4.13215. PMCID PMC11494827.",
    "Parameter source: Appendix S1 Table S6 (Morphine Population PK-PD",
    "Parameter Estimates) and Appendix S1 example NONMEM control stream",
    "(B) Morphine, whose $THETA / $OMEGA blocks and $DES equations this",
    "file reproduces one-for-one. Where the two disagree in precision",
    "the control stream is used, since it carries more digits",
    "(kon 0.002989416 vs the tabulated 0.00299; C50 naloxone 0.6021768",
    "vs the tabulated 0.602).",
    "Morphine PK-PD parameters originate from Olofsen E, van Dorp E,",
    "Teppema L, Aarts L, Smith TW, Dahan A, Sarton E. Naloxone reversal",
    "of morphine- and morphine-6-glucuronide-induced respiratory",
    "depression in healthy volunteers: a mechanism-based",
    "pharmacokinetic-pharmacodynamic modeling study.",
    "Anesthesiology. 2010;112(6):1417-1427.",
    "doi:10.1097/ALN.0b013e3181d5e29d.",
    "Naloxone PK parameters are Yang 2024 Appendix S1 Table S4."
  )
  vignette <- "Yang_2024_naloxone_opioid_reversal"
  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ng/mL"
  )
  dosing <- c("central", "depot_naloxone")

  compartmentData <- list(
    central              = list(analyte = "morphine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "morphine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2          = list(analyte = "morphine", units = "ug", specimen = "plasma", verified = TRUE),
    effect               = list(analyte = "morphine biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    depot_naloxone       = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit1_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit2_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit3_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    central_naloxone     = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_naloxone = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    effect_naloxone      = list(analyte = "naloxone biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    RL_op                = list(analyte = "apparent fractional morphine occupancy of the mu-opioid receptor pool", units = "fraction", specimen = "not applicable", verified = TRUE)
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
        "weight of 70 kg is explicit in control stream (B)",
        "('WT= 70 ; kg'). The morphine micro-constants carry no weight",
        "scaling in Olofsen 2010."
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
      "Healthy adult volunteers. The morphine PK-PD layer comes from",
      "the Olofsen 2010 healthy-volunteer study of naloxone reversal of",
      "morphine- and morphine-6-glucuronide-induced respiratory",
      "depression; the naloxone PK layer is Yang 2024's own",
      "auto-injector population PK in 48 healthy adults."
    ),
    dose_range     = paste(
      "Yang 2024 simulated IV morphine at 0.2, 2.2 and 4.2 mg/kg (the",
      "lowest corresponding to 50 percent ventilation suppression, the",
      "highest 21-fold higher and representing a lethal dose), each",
      "with naloxone auto-injector 2 mg or 10 mg IM/SC given once",
      "ventilation had fallen 60 percent from baseline."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Yang 2024 validated this constructed model by digitising the",
      "ventilation-time profile for morphine 10.5 mg/70 kg from",
      "Figure 1 of Olofsen 2010 and confirming that most observed",
      "points fell inside the 90 percent prediction interval of 1000",
      "simulations (Appendix S1 Figure S8).",
      "Baseline ventilation V0 was carried as a per-subject data column",
      "(E0) in the control stream rather than as a THETA; it is encoded",
      "here as the canonical le0 parameter with the Table S6 typical",
      "value and IIV so that the model is self-contained.",
      "Morphine disposition is parameterised in micro-constants rather",
      "than clearances because that is how Olofsen 2010 reported it and",
      "because the IIV terms sit on the micro-constants; the equivalent",
      "macro-constants are CL = k10*vc = 1.574 L/min, Q2 = k12*vc =",
      "0.905 L/min, V2 = Q2/k21 = 9.97 L, Q3 = k13*vc = 1.528 L/min and",
      "V3 = Q3/k31 = 114.0 L (the same conversion control stream (B)",
      "performs in $PK for reporting)."
    )
  )

  ini({
    # Every value below is held FIX in Appendix S1 control stream (B):
    # the morphine PK-PD block is inherited unchanged from Olofsen 2010
    # and the naloxone PK block from Yang 2024's own Table S4, and
    # neither was re-fitted for this analysis.

    # --- Naloxone PK (Yang 2024 Table S4; control stream (B) THETA 1-6)
    lktr_naloxone <- fixed(log(0.696))
    label("Naloxone transit absorption rate constant (KTR, 1/min)")            # (B) THETA 1 = 0.696; Table S4
    lcl_naloxone <- fixed(log(3.26))
    label("Naloxone apparent clearance at 70 kg (CL/F, L/min)")                # (B) THETA 2 = 3.26; Table S4
    e_wt_cl_naloxone <- fixed(0.538)
    label("Allometric exponent of body weight on naloxone CL/F (unitless)")    # (B) THETA 3 = 0.538; Table S4
    lvc_naloxone <- fixed(log(404))
    label("Naloxone apparent central volume (V2/F, L)")                        # (B) THETA 4 = 404; Table S4
    lq_naloxone <- fixed(log(0.0847))
    label("Naloxone apparent intercompartmental clearance (Q/F, L/min)")       # (B) THETA 5 = 0.0847; Table S4
    lvp_naloxone <- fixed(log(81.8))
    label("Naloxone apparent peripheral volume (V3/F, L)")                     # (B) THETA 6 = 81.8; Table S4

    # --- Morphine PK, micro-constant parameterisation (Table S6;
    # control stream (B) THETA 7-12). Table S6 labels k10 / k12 / k21 /
    # k13 / k31 as "L/min", which is a units typo: control stream (B)
    # annotates each of them "[1/min]", and k10 * V1 = 0.308 * 5.11 =
    # 1.57 L/min recovers the textbook morphine clearance, which the
    # tabulated label could not.
    lvc <- fixed(log(5.11))
    label("Morphine central volume (V1, L)")                                   # (B) THETA 7 = 5.11; Table S6
    lkel <- fixed(log(0.308))
    label("Morphine elimination micro-constant (k10, 1/min)")                  # (B) THETA 8 = 0.308 [1/min]; Table S6
    lk12 <- fixed(log(0.177))
    label("Morphine central-to-peripheral1 micro-constant (k12, 1/min)")       # (B) THETA 9 = 0.177 [1/min]; Table S6
    lk21 <- fixed(log(0.0907))
    label("Morphine peripheral1-to-central micro-constant (k21, 1/min)")       # (B) THETA 10 = 0.0907 [1/min]; Table S6
    lk13 <- fixed(log(0.299))
    label("Morphine central-to-peripheral2 micro-constant (k13, 1/min)")       # (B) THETA 11 = 0.299 [1/min]; Table S6
    lk31 <- fixed(log(0.0134))
    label("Morphine peripheral2-to-central micro-constant (k31, 1/min)")       # (B) THETA 12 = 0.0134 [1/min]; Table S6

    # --- Biophase equilibration (equation 1). Olofsen 2010 reports the
    # plasma-effect-site equilibration HALF-LIFE; control stream (B)
    # converts it in $DES as ke0 = 0.693 / t(1/2)ke0. The conversion is
    # done here so that ke0 is the stored parameter.
    lke0_naloxone <- fixed(log(0.693 / 11.2))
    label("Naloxone biophase equilibration rate constant (ke0, 1/min)")        # (B) THETA 13 = 11.2 min half-life -> 0.06188 1/min; Table S6
    lke0 <- fixed(log(0.693 / 74.4))
    label("Morphine biophase equilibration rate constant (ke0, 1/min)")        # (B) THETA 14 = 74.4 min half-life -> 0.009315 1/min; Table S6

    # --- MOP receptor association / dissociation (equations 4, 8, 9)
    lkon <- fixed(log(0.002989416))
    label("Morphine receptor association rate constant (kon, mL/ng/min)")      # (B) THETA 15 = 0.002989416 (Table S6 rounds to 0.00299)
    lkoff <- fixed(log(0.138))
    label("Morphine receptor dissociation rate constant (koff, 1/min)")        # (B) THETA 16 = 0.138; Table S6
    lec50_naloxone <- fixed(log(0.6021768))
    label("Naloxone half-maximal effect-site concentration (C50, ng/mL)")      # (B) THETA 17 = 0.6021768 (Table S6 rounds to 0.602)
    lhill_naloxone <- fixed(log(4.18))
    label("Naloxone competitive-displacement shape parameter (gamma, unitless)")  # (B) THETA 18 = 4.18; Table S6

    # --- Linear transduction function (equation 10)
    lalpha <- fixed(log(1))
    label("Morphine intrinsic activity (alpha, unitless 0-1)")                 # Table S6: alpha = 1 (control stream (B) $ERROR writes E = E0*(1 - A(12)) with no alpha term, i.e. alpha = 1)
    le0 <- fixed(log(26.5))
    label("Baseline minute ventilation (V0, L/min)")                           # Table S6: V0 26.5 L/min (carried as the E0 data column in control stream (B))

    # --- IIV. Control stream (B) $OMEGA entries 1-13, in order. Each
    # equals the square of the Table S6 "IIV (%)" column (e.g. morphine
    # k10 11.3% -> 0.0128 = 0.113^2). k12 has no IIV, matching the "-"
    # in Table S6.
    etalktr_naloxone ~ fixed(0.111)
    etalcl_naloxone ~ fixed(0.0129)
    etalvc_naloxone ~ fixed(0.0658)
    etalvc ~ fixed(0.0179)
    etalkel ~ fixed(0.0128)
    etalk21 ~ fixed(0.0446)
    etalk13 ~ fixed(0.0207)
    etalk31 ~ fixed(0.0312)
    # Control stream (B) places this eta on the morphine equilibration
    # HALF-LIFE, whereas the parameter stored above is the rate
    # constant. Because ke0 = 0.693/t(1/2), an exp(eta) multiplier on
    # the half-life is exactly an exp(-eta) multiplier on ke0; the
    # log-normal distribution is symmetric in the exponent, so the same
    # variance on lke0 reproduces the identical population spread.
    etalke0 ~ fixed(0.16)
    etalkon ~ fixed(0.192)
    etalkoff ~ fixed(0.112)
    etalec50_naloxone ~ fixed(0.141)
    etalhill_naloxone ~ fixed(0.928)
    # Baseline ventilation IIV. Control stream (B) has no $OMEGA entry
    # for it because E0 arrived as a per-subject data column already
    # carrying the variability; Table S6 reports V0 IIV = 12%, giving
    # the variance 0.12^2 = 0.0144.
    etale0 ~ fixed(0.0144)

    # Table S6 reports no residual error on the morphine PK endpoint
    # (the row is absent, and control stream (B) fits only ventilation
    # with the PK held fixed). Encoded as an explicit zero rather than
    # invented; see the vignette Assumptions and deviations section.
    propSd <- fixed(0)
    label("Morphine proportional residual error SD (fraction; not reported)")  # Table S6 reports no PK residual-error row
    addSd <- fixed(2.18)
    label("Ventilation additive residual error SD (L/min)")                    # (B) THETA 19 = 2.18; Table S6
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
    hill_naloxone <- exp(lhill_naloxone + etalhill_naloxone)

    # ---- Morphine individual parameters (micro-constants)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21 + etalk21)
    k13 <- exp(lk13 + etalk13)
    k31 <- exp(lk31 + etalk31)
    ke0 <- exp(lke0 + etalke0)
    kon <- exp(lkon + etalkon)
    koff <- exp(lkoff + etalkoff)
    alpha <- exp(lalpha)
    e0 <- exp(le0 + etale0)

    # ---- Naloxone PK, control stream (B) $DES DADT(1)-DADT(6).
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
    d/dt(effect_naloxone) <- ke0_naloxone * (central_naloxone / vc_naloxone - effect_naloxone)

    # ---- Morphine PK, control stream (B) $DES DADT(8)-DADT(10),
    # written directly in micro-constants as the source does.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(effect) <- ke0 * (central / vc - effect)

    # ---- MOP receptor kinetics with competitive naloxone displacement.
    # Control stream (B) DADT(12):
    #   NRATIO = Ce_nal / C50_nal
    #   dRL/dt = kon*Ce_op*(1 - RL)*(1 - NRATIO^gamma/(1 + NRATIO)^gamma)
    #            - koff*RL
    # Note the parenthesisation: the exponent gamma is applied to the
    # WHOLE (1 + NRATIO) denominator, so the naloxone occupancy term is
    # (NRATIO/(1 + NRATIO))^gamma, not the textbook sigmoid
    # NRATIO^gamma/(1 + NRATIO^gamma). Appendix S1 equation 8 is printed
    # the same way and the control stream implements it literally, so
    # this form -- not the textbook one -- is what produced the
    # published results. See the vignette Assumptions and deviations
    # section, which notes that the Table S6 footnote calling C50 "the
    # concentration causing 50% effect" does not hold under this form.
    d/dt(RL_op) <- kon * effect * (1 - RL_op) *
      (1 - (effect_naloxone / ec50_naloxone)^hill_naloxone /
        (1 + effect_naloxone / ec50_naloxone)^hill_naloxone) -
      koff * RL_op

    # ---- Outputs
    Cc <- central / vc
    Cnal <- central_naloxone / vc_naloxone
    VE <- e0 * (1 - alpha * RL_op)
    ERATIO <- VE / e0

    Cc ~ prop(propSd)
    VE ~ add(addSd)
  })
}
