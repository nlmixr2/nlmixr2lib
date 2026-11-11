Ynineb_2025_propofol_remifentanil <- function() {
  description <- paste(
    "Augmented propofol pharmacokinetic-pharmacodynamic model with a",
    "body-mass-index-driven fat 'trap' compartment, plus classical",
    "remifentanil pharmacokinetics and a Greco response-surface model for the",
    "Bispectral Index (BIS) under total intravenous anesthesia (Ynineb 2025).",
    "Propofol follows the Schnider three-compartment model (blood, muscle, fat;",
    "covariates age, weight, height and James lean body mass) in which the fat",
    "compartment exchanges with a trap sub-compartment of volume",
    "V3 * BMI / 100 at clearance Cl3 / Risk, where Risk is a three-term",
    "sum-of-sines function of BMI. Each drug drives its own effect",
    "compartment; BIS = E0 - Emax * I^gamma / (1 + I^gamma) with",
    "I = Ce_propofol / C50p + Ce_remifentanil / C50r + interaction term. No",
    "parameter was estimated in the source: the PK parameters come from the",
    "Schnider (propofol) and Minto (remifentanil) models and the PD",
    "parameters are patient-specific values identified in earlier work by the",
    "same group. Defaults are those of the paper's clinical case (VitalDB",
    "patient 5428), including that patient's propofol ke0 of 1.2 1/min",
    "(0.456 1/min for the paper's other patient); the remifentanil PK is",
    "fixed at that patient's printed values because the source does not",
    "print the Minto covariate equations.",
    "The paper's alternative fractional-order (partial Caputo) fat-release",
    "model is not representable as an ordinary differential equation and is",
    "not encoded.",
    sep = " "
  )
  reference <- paste(
    "Ynineb AR, Yumuk E, Copot D, Ayvaz B, Khoumeri B, Ben Othman G, Mihai MD,",
    "Birs IR, Muresan CI, Ionescu CM. (2025). Modeling drug retention as",
    "memory effects in obese patients using fractional and augmented models.",
    "Scientific Reports 15:26234. doi:10.1038/s41598-025-10172-1.",
    "PMCID: PMC12276283. Schnider propofol covariate equations and the James",
    "lean-body-mass formulas as printed in the companion paper: Ynineb AR,",
    "Yumuk E, Copot D, et al. (2026). Pharmacokinetic modelling during",
    "long-term anesthesia: minimizing the gap. Journal of Advanced Research",
    "82:535-561. doi:10.1016/j.jare.2025.06.047. PMCID: PMC13001178.",
    sep = " "
  )
  vignette <- "Ynineb_2025_propofol_remifentanil"
  units <- list(
    time = "min",
    dosing = "mg (propofol, cmt central); ug (remifentanil, cmt central_remifentanil)",
    concentration = "mg/L = ug/mL (propofol); ug/L = ng/mL (remifentanil)"
  )
  # Time is minutes because every clearance in the source is in L/min and
  # every rate constant in 1/min (Tables 5 and 8). The propofol input u(t) is
  # in mg/min and V1 in L, so propofol concentrations are mg/L (= ug/mL, the
  # unit of C50p in Tables 4 and 7). The remifentanil dose unit is not stated
  # in the source; ug is assumed so that its concentration is in ng/mL, the
  # conventional unit for remifentanil C50.

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effects on the propofol muscle volume V2 and muscle",
        "clearance Cl2, centred at 53 years (Schnider model, companion paper",
        "Eq. 6). Patient 21 was 62 years (Table 4), patient 5428 61 years",
        "(Table 7)."
      ),
      source_name = "Age"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on propofol Cl1 centred at 77 kg (Schnider model,",
        "companion paper Eq. 6); also enters the James lean body mass and",
        "the body mass index BMI = WT / (HT / 100)^2 that drives the trap",
        "compartment."
      ),
      source_name = "Weight"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on propofol Cl1 centred at 177 cm (Schnider model,",
        "companion paper Eq. 6); also enters the James lean body mass and",
        "BMI (converted to metres inside model())."
      ),
      source_name = "Height"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the James lean-body-mass formula: female",
        "1.07 * WT - 148 * (WT / HT)^2, male 1.1 * WT - 128 * (WT / HT)^2",
        "(companion paper Eqs. 7-8). Tables 4 and 7 code sex as F / M."
      ),
      source_name = "Gender"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "not applicable", verified = FALSE),
    peripheral2 = list(analyte = "propofol", units = "mg", specimen = "not applicable", verified = FALSE),
    peripheral2_deep = list(analyte = "propofol", units = "mg", specimen = "not applicable", verified = FALSE),
    effect = list(analyte = "propofol", units = "mg/L", specimen = "not applicable", verified = FALSE),
    central_remifentanil = list(analyte = "remifentanil", units = "ug", specimen = "plasma", verified = FALSE),
    peripheral1_remifentanil = list(
      analyte = "remifentanil",
      units = "ug",
      specimen = "not applicable",
      verified = FALSE
    ),
    peripheral2_remifentanil = list(
      analyte = "remifentanil",
      units = "ug",
      specimen = "not applicable",
      verified = FALSE
    ),
    effect_remifentanil = list(analyte = "remifentanil", units = "ug/L", specimen = "not applicable", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 2L,
    n_studies = 2L,
    age_range = "61-62 years",
    weight_range = "81-88 kg",
    height_range = "166-168 cm",
    bmi_range = "29.4-31.2 kg/m^2",
    sex_female_pct = 50,
    disease_state = paste(
      "Adults with overweight or obesity under total intravenous anesthesia:",
      "patient 21 of a 70-patient Ghent University Hospital database",
      "(female, BMI 31.2, hysterectomy; propofol only) and patient 5428 of the",
      "public VitalDB database (male, BMI 29.4, abdominal aortic aneurysm",
      "repair; propofol plus remifentanil)."
    ),
    dose_range = paste(
      "Simulated propofol protocols for patient 21: 120 mg over 0.6 min then",
      "a maintenance infusion from 2 to 6.4 min; a 100 mg bolus over 1 min;",
      "four 48 mg boluses over 20 s every 2 min; and a continuous infusion",
      "from 2 to 20 min (371.52 mg in total). Patient 5428 received the",
      "propofol and remifentanil infusions recorded in VitalDB over a",
      "350 min surgery."
    ),
    regions = "Belgium (Ghent University Hospital); South Korea (VitalDB, Seoul National University Hospital)",
    notes = paste(
      "Simulation and single-patient evaluation study; no population fit.",
      "Biometrics from Tables 4 and 7. Patient 21's PD values match the",
      "70-patient database of the companion paper (Table B.9, index 21).",
      "BIS was recorded at 1 s intervals for patient 5428 and compared with",
      "the simulated BIS during induction (0-30 min) and awakening",
      "(330 min to end)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Propofol: Schnider three-compartment model. The source paper
    # computes these 'based on the Schnider model' (Results) without
    # printing the equations; they are taken from the same group's
    # companion paper (J Adv Res 2026, Eq. 6) and reproduce every Cl, V
    # and k in Tables 5 and 8 of the source exactly. Nothing was
    # estimated in the source, so every value is fixed().
    # Cl1 = 1.89 + 0.0456 (WT - 77) - 0.0681 (LBM - 59) + 0.0264 (HT - 177)
    # -----------------------------------------------------------------
    lcl <- fixed(log(1.89))
    label("Propofol elimination clearance Cl1 at WT 77 kg, LBM 59 kg, HT 177 cm (L/min)") # Companion Eq. 6: Cl1 intercept 1.89
    e_wt_cl <- fixed(0.0456)
    label("Linear slope of Cl1 on weight (L/min per kg)") # Companion Eq. 6: +0.0456 (weight - 77)
    e_lbm_cl <- fixed(-0.0681)
    label("Linear slope of Cl1 on lean body mass (L/min per kg)") # Companion Eq. 6: -0.0681 (lbm - 59)
    e_ht_cl <- fixed(0.0264)
    label("Linear slope of Cl1 on height (L/min per cm)") # Companion Eq. 6: +0.0264 (height - 177)
    lvc <- fixed(log(4.27))
    label("Propofol blood (central) volume V1 (L)") # Companion Eq. 6: V1 = 4.27; Table 5 V1 = 4.27
    lq <- fixed(log(1.29))
    label("Propofol blood-muscle clearance Cl2 at age 53 (L/min)") # Companion Eq. 6: Cl2 intercept 1.29
    e_age_q <- fixed(-0.024)
    label("Linear slope of Cl2 on age (L/min per year)") # Companion Eq. 6: -0.024 (age - 53)
    lvp <- fixed(log(18.9))
    label("Propofol muscle volume V2 at age 53 (L)") # Companion Eq. 6: V2 intercept 18.9
    e_age_vp <- fixed(-0.391)
    label("Linear slope of V2 on age (L per year)") # Companion Eq. 6: -0.391 (age - 53)
    lq2 <- fixed(log(0.836))
    label("Propofol blood-fat clearance Cl3 (L/min)") # Companion Eq. 6: Cl3 = 0.836; Table 5 Cl3 = 0.836
    lvp2 <- fixed(log(238))
    label("Propofol fat volume V3 (L)") # Companion Eq. 6: V3 = 238; Table 5 V3 = 238
    # Propofol k1e = ke0 is patient-specific in the source: 0.456 1/min
    # (the Schnider value) for patient 21 in Table 5, but 1.2 1/min for
    # patient 5428 in Table 8. The defaults are patient 5428's.
    lke0 <- fixed(log(1.2))
    label("Propofol effect-site rate constants k1e = ke0, patient 5428 (1/min)") # Table 8 row 1: k1e = ke0 = 1.2 (Table 5, patient 21: 0.456)

    # -----------------------------------------------------------------
    # Fat trap. Eq. 4: k3t = Clt / Vt, kt3 = Clt / V3; Eq. 5:
    # Vt = V3 * BMI / 100, Clt = Cl3 / R; Eq. 6:
    # R = sum_i a_i sin(b_i BMI + c_i), coefficients from Table 3 (the
    # 3-sine fit, also companion Table A.6 'Approximation 1').
    # -----------------------------------------------------------------
    risk_a1 <- fixed(27.4472)
    label("Trap-risk sum-of-sines amplitude a1 (unitless)") # Table 3: a1 = 27.4472
    risk_b1 <- fixed(0.0033)
    label("Trap-risk sum-of-sines frequency b1 (per kg/m^2)") # Table 3: b1 = 0.0033
    risk_c1 <- fixed(-0.0570)
    label("Trap-risk sum-of-sines phase c1 (rad)") # Table 3: c1 = -0.0570
    risk_a2 <- fixed(0.8002)
    label("Trap-risk sum-of-sines amplitude a2 (unitless)") # Table 3: a2 = 0.8002
    risk_b2 <- fixed(0.5390)
    label("Trap-risk sum-of-sines frequency b2 (per kg/m^2)") # Table 3: b2 = 0.5390
    risk_c2 <- fixed(-4.2180)
    label("Trap-risk sum-of-sines phase c2 (rad)") # Table 3: c2 = -4.2180
    risk_a3 <- fixed(0.4815)
    label("Trap-risk sum-of-sines amplitude a3 (unitless)") # Table 3: a3 = 0.4815
    risk_b3 <- fixed(0.7595)
    label("Trap-risk sum-of-sines frequency b3 (per kg/m^2)") # Table 3: b3 = 0.7595
    risk_c3 <- fixed(-6.2212)
    label("Trap-risk sum-of-sines phase c3 (rad)") # Table 3: c3 = -6.2212

    # -----------------------------------------------------------------
    # Remifentanil: classical three-compartment model (Eq. 1), Minto
    # parameters for patient 5428 as printed in Table 8, row 2. The
    # Minto covariate equations are not printed in the source or the
    # companion paper, so these are that patient's values. Table 8
    # prints V1 = 4.27 for remifentanil, which is the propofol value
    # repeated: the same row's k10 = 0.4716 = Cl1 / V1 needs
    # V1 = 2.329 / 0.4716 = 4.9385 L, as do k12 = Cl2 / V1 = 0.2871 and
    # k13 = Cl3 / V1 = 0.0106, so V1 is taken from the rate constants.
    # -----------------------------------------------------------------
    lcl_remifentanil <- fixed(log(2.329))
    label("Remifentanil elimination clearance Cl1, patient 5428 (L/min)") # Table 8 row 2: Cl1 = 2.329
    lvc_remifentanil <- fixed(log(4.9385))
    label("Remifentanil central volume V1, patient 5428 (L)") # Table 8 row 2: k10 = 0.4716 with Cl1 = 2.329 gives V1 = 4.9385 (printed 4.27 is the propofol value)
    lq_remifentanil <- fixed(log(1.4179))
    label("Remifentanil fast-peripheral clearance Cl2, patient 5428 (L/min)") # Table 8 row 2: Cl2 = 1.4179
    lvp_remifentanil <- fixed(log(8.5082))
    label("Remifentanil fast-peripheral volume V2, patient 5428 (L)") # Table 8 row 2: V2 = 8.5082
    lq2_remifentanil <- fixed(log(0.0523))
    label("Remifentanil slow-peripheral clearance Cl3, patient 5428 (L/min)") # Table 8 row 2: Cl3 = 0.0523
    lvp2_remifentanil <- fixed(log(5.42))
    label("Remifentanil slow-peripheral volume V3, patient 5428 (L)") # Table 8 row 2: V3 = 5.42
    lke0_remifentanil <- fixed(log(0.448))
    label("Remifentanil effect-site rate constants k1e = ke0, patient 5428 (1/min)") # Table 8 row 2: k1e = ke0 = 0.448

    # -----------------------------------------------------------------
    # BIS response surface (Eqs. 25 and 27), patient 5428 (Table 7).
    # Emax is not printed for either patient. Pairing the effect-site
    # concentration and BIS curves for patient 21 (Figs. 4-5) and
    # re-running the paper's three protocols (Table 6 BIS nadirs) both
    # require Emax = E0, so Emax is set equal to E0 (see the vignette).
    # -----------------------------------------------------------------
    le0 <- fixed(log(90.72))
    label("Baseline BIS E0, patient 5428 (BIS units)") # Table 7: E0 = 90.72
    lemax <- fixed(log(90.72))
    label("Maximum BIS reduction Emax, patient 5428 (BIS units)") # Not printed; set equal to E0 (Table 7: E0 = 90.72), back-solved from Figs. 4-5 and Table 6
    lec50 <- fixed(log(4.8053))
    label("Propofol effect-site concentration at half-maximal effect C50p, patient 5428 (ug/mL)") # Table 7: C50p = 4.8053 ug/mL
    lec50_remifentanil <- fixed(log(6))
    label("Remifentanil effect-site concentration at half-maximal effect C50r, patient 5428 (ng/mL)") # Table 7: C50r = 6
    lhill <- fixed(log(1.6305))
    label("Hill coefficient gamma of the response surface, patient 5428 (unitless)") # Table 7: gamma = 1.6305
    greco_int <- fixed(0)
    label("Greco interaction parameter Gamma, patient 5428 (unitless; 0 = additive)") # Table 7: Gamma = 0
  })

  model({
    # Body-size quantities. BMI = WT / HT^2 with HT in metres (Eq. 5 and
    # text); James lean body mass with HT in cm (companion Eqs. 7-8).
    bmi <- WT / (HT / 100)^2
    lbm <- SEXF * (1.07 * WT - 148 * (WT / HT)^2) +
      (1 - SEXF) * (1.1 * WT - 128 * (WT / HT)^2)

    # Propofol Schnider clearances and volumes (companion Eq. 6), additive
    # linear covariate effects.
    cl <- exp(lcl) +
      e_wt_cl * (WT - 77) +
      e_lbm_cl * (lbm - 59) +
      e_ht_cl * (HT - 177)
    vc <- exp(lvc)
    q <- exp(lq) + e_age_q * (AGE - 53)
    vp <- exp(lvp) + e_age_vp * (AGE - 53)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    ke0 <- exp(lke0)

    # Trap volume and clearance (Eqs. 4-6).
    risk <- risk_a1 * sin(risk_b1 * bmi + risk_c1) +
      risk_a2 * sin(risk_b2 * bmi + risk_c2) +
      risk_a3 * sin(risk_b3 * bmi + risk_c3)
    vp2_deep <- vp2 * bmi / 100
    cl_deep <- q2 / risk

    # Micro-constants (Eq. 2 and Eq. 4).
    kel <- cl / vc
    k12 <- q / vc
    k13 <- q2 / vc
    k21 <- q / vp
    k31 <- q2 / vp2
    k3t <- cl_deep / vp2_deep
    kt3 <- cl_deep / vp2

    # Propofol states (Eqs. 7a-7d). The paper writes the states as
    # concentrations x_i, but with the micro-constants above Eq. 7 is
    # mass-balanced only for x_i = amount_i / V1, so the states here are
    # amounts and the paper's x2, x3 and xt are peripheral1 / vc,
    # peripheral2 / vc and peripheral2_deep / vc.
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2 - k3t * peripheral2 + kt3 * peripheral2_deep
    d/dt(peripheral2_deep) <- k3t * peripheral2 - kt3 * peripheral2_deep

    # Remifentanil (Eq. 1 with Minto parameters).
    cl_remifentanil <- exp(lcl_remifentanil)
    vc_remifentanil <- exp(lvc_remifentanil)
    q_remifentanil <- exp(lq_remifentanil)
    vp_remifentanil <- exp(lvp_remifentanil)
    q2_remifentanil <- exp(lq2_remifentanil)
    vp2_remifentanil <- exp(lvp2_remifentanil)
    ke0_remifentanil <- exp(lke0_remifentanil)
    kel_remifentanil <- cl_remifentanil / vc_remifentanil
    k12_remifentanil <- q_remifentanil / vc_remifentanil
    k13_remifentanil <- q2_remifentanil / vc_remifentanil
    k21_remifentanil <- q_remifentanil / vp_remifentanil
    k31_remifentanil <- q2_remifentanil / vp2_remifentanil
    d/dt(central_remifentanil) <- -(kel_remifentanil + k12_remifentanil + k13_remifentanil) * central_remifentanil +
      k21_remifentanil * peripheral1_remifentanil + k31_remifentanil * peripheral2_remifentanil
    d/dt(peripheral1_remifentanil) <- k12_remifentanil * central_remifentanil - k21_remifentanil * peripheral1_remifentanil
    d/dt(peripheral2_remifentanil) <- k13_remifentanil * central_remifentanil - k31_remifentanil * peripheral2_remifentanil

    Cc <- central / vc
    Cc_remifentanil <- central_remifentanil / vc_remifentanil

    # Effect sites (Eq. 24) with k1e = ke0 for each drug.
    d/dt(effect) <- ke0 * Cc - ke0 * effect
    d/dt(effect_remifentanil) <- ke0_remifentanil * Cc_remifentanil - ke0_remifentanil * effect_remifentanil

    # Greco response surface (Eq. 27) and BIS (Eq. 25). Eq. 27 as printed
    # divides BOTH effect-site concentrations by C50 of drug 2; the
    # interaction term and the Greco model divide each drug by its own C50,
    # which is used here. With no remifentanil, I = effect / ec50 and BIS
    # reduces to the single-drug Hill model of Eq. 26.
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    ec50_remifentanil <- exp(lec50_remifentanil)
    hill <- exp(lhill)
    u_prop <- effect / ec50
    u_remi <- effect_remifentanil / ec50_remifentanil
    inter <- u_prop + u_remi + greco_int * u_prop * u_remi
    BIS <- e0 - emax * inter^hill / (1 + inter^hill)
  })
}
