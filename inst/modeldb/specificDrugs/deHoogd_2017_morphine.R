deHoogd_2017_morphine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for morphine and its two",
    "glucuronide metabolites (M3G, M6G) in 20 morbidly obese adults",
    "(post-gastric-bypass) and 20 healthy adult volunteers (de Hoogd 2017).",
    "Morphine: three-compartment IV model with total body weight (TBW)",
    "covariate on the second peripheral volume V5M. Non-glucuronide",
    "morphine clearance is structurally fixed at 35% of total morphine CL",
    "in a 70-kg healthy adult. M3G and M6G are each one-compartment models",
    "fed by formation-delay transit chains (n = 5 for M3G, n = 2 for M6G);",
    "VM3G = VM6G is a structural equality. TBW covariates apply to CLF M6G,",
    "the M3G transit rate Ktr, M3G elimination CL, and M6G elimination CL,",
    "all power-form normalised to a reference of 98.5 kg (population median).",
    "Proportional residual error is reported separately for the healthy-",
    "volunteer cohort and the morbidly obese cohort, selected via the",
    "binary indicator DIS_MORBOBESE."
  )
  reference <- paste(
    "de Hoogd S, Valitalo PAJ, Dahan A, van Kralingen S,",
    "Coughtrie MMW, van Dongen EPA, van Ramshorst B, Knibbe CAJ.",
    "Influence of morbid obesity on the pharmacokinetics of morphine,",
    "morphine-3-glucuronide, and morphine-6-glucuronide.",
    "Clin Pharmacokinet. 2017;56(12):1577-1587.",
    "doi:10.1007/s40262-017-0544-2.",
    "ClinicalTrials.gov NCT01097148.",
    sep = " "
  )
  vignette <- "deHoogd_2017_morphine"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "nmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed total body weight (TBW). Power-form covariate on M6G",
        "formation clearance (CLF M6G, exponent -0.329), morphine peripheral",
        "volume V5M (exponent 0.483), M3G transit rate Ktr (exponent -0.701),",
        "M3G elimination clearance (exponent -1.08), and M6G elimination",
        "clearance (exponent -1.03). All covariates normalised to a reference",
        "TBW of 98.5 kg (population median across pooled cohort). The paper",
        "reports no significant TBW effect on morphine clearance itself."
      ),
      source_name        = "TBW"
    ),
    DIS_MORBOBESE = list(
      description        = "Morbidly obese cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = paste(
        "1 = morbidly obese surgical patient (post gastric bypass / banding /",
        "sleeve; BMI > 40 kg/m^2), 0 = non-obese healthy volunteer. Selects",
        "the cohort-specific proportional residual error for each of the",
        "three observed species (morphine, M3G, M6G); the paper estimated",
        "separate residual variabilities for the healthy-volunteer and",
        "morbidly obese cohorts (Table 2, residual variability rows)."
      ),
      source_name        = "morbid obesity indicator (paper Table 2 cohort split)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40,
    n_studies      = 3,
    age_range      = "20-59 years",
    age_median     = "morbidly obese 44.1 +/- 10.6 years; healthy volunteers 25.5 +/- 4.1 years",
    weight_range   = "56-251.9 kg",
    weight_median  = "98.5 kg (population median across pooled cohort)",
    sex_female_pct = 52.5,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Pooled cohort of 20 morbidly obese adults (BMI 37.9-78.6 kg/m^2,",
      "weight 112-251.9 kg) undergoing laparoscopic gastric bypass / banding /",
      "sleeve surgery, plus 20 historical-control healthy adult volunteers",
      "(weight 56-85 kg) from two prior morphine PK studies. All subjects",
      "had normal renal and liver function. ASA II/III."
    ),
    dose_range     = paste(
      "Morbidly obese: 10 mg IV bolus of morphine HCl at end of surgery",
      "(additional postoperative boluses as needed; mean total 15.7 mg).",
      "Healthy volunteers: 0.10 mg/kg IV bolus followed by 0.030 mg/kg/h",
      "infusion for 1 h (mean total 9.2 mg)."
    ),
    regions        = "The Netherlands (St Antonius Hospital + Leiden University Medical Center)",
    notes          = paste(
      "Demographics from de Hoogd 2017 Table 1. NONMEM 7.2 with FOCE-INTER;",
      "modeled concentrations expressed in nmol/L using morphine MW 285.33",
      "g/mol and M3G/M6G MW 461.46 g/mol; administered dose corrected from",
      "morphine HCl (MW 321.8 g/mol). Sampling 0-420 minutes post first",
      "morphine dose; 11 samples per obese subject, 15 samples per healthy",
      "volunteer. Identifiability of the full structural model was verified",
      "via COMBOS (Electronic Supplementary Material 1)."
    )
  )

  ini({
    # ============================================================
    # Morphine structural parameters -- de Hoogd 2017 Table 2 (Final model)
    # ============================================================
    lcl_form_m3g <- log(0.748)
    label("Morphine -> M3G formation clearance (L/min)")              # Table 2 Final: CLF M3G = 0.748 L/min (no TBW covariate)
    lcl_form_m6g <- log(0.129)
    label("Morphine -> M6G formation clearance at TBW = 98.5 kg (L/min)") # Table 2 Final: CLF M6G, 98.5 kg = 0.129 L/min
    e_wt_cl_form_m6g <- -0.329
    label("TBW power exponent on CLF M6G (unitless)")                  # Table 2 Final: K = -0.329 (RSE 36%)

    lvc <- log(4.62)
    label("Morphine central volume V1M (L)")                           # Table 2 Final: V1M = 4.62 L
    lvp <- log(9.52)
    label("Morphine first peripheral volume V4M (L)")                  # Table 2 Final: V4M = 9.52 L
    lvp2 <- log(118)
    label("Morphine second peripheral volume V5M at TBW = 98.5 kg (L)") # Table 2 Final: V98.5kg = 118 L
    e_wt_vp2 <- 0.483
    label("TBW power exponent on V5M (unitless)")                      # Table 2 Final: L = 0.483 (RSE 48%)

    lq <- log(0.814)
    label("Morphine inter-compartmental clearance Q2 (V1M <-> V4M) (L/min)")  # Table 2 Final: Q2 = 0.814 L/min
    lq2 <- log(1.29)
    label("Morphine inter-compartmental clearance Q3 (V1M <-> V5M) (L/min)")  # Table 2 Final: Q3 = 1.29 L/min

    # Non-glucuronide CL is structurally fixed (Methods 2.4): 35% of
    # total morphine clearance in a 70-kg healthy subject. Total CL at
    # 70 kg = (CLF_M3G + CLF_M6G_70kg) / 0.65, where CLF_M6G_70kg =
    # 0.129 * (70/98.5)^-0.329 = 0.1444. Total CL_70kg = 0.8924/0.65 =
    # 1.3729 L/min; CL_nongluc = 0.35 * 1.3729 = 0.4805 L/min. The
    # paper reports no TBW effect on morphine CL, so this constant is
    # used unchanged for every subject.
    lcl_nongluc <- fixed(log(0.4805))
    label("Morphine non-glucuronide elimination clearance (L/min, FIXED structural)")  # Methods 2.4: 35% of CL_total(70 kg)

    # ============================================================
    # M3G formation-delay transit chain
    # ============================================================
    lktr_m3g <- log(1.68)
    label("Morphine -> M3G transit-chain rate Ktr at TBW = 98.5 kg (1/min)") # Table 2 Final: Ktr_98.5kg = 1.68 /min (n = 5 transits, MTT = 5/Ktr ~ 2.98 min at 98.5 kg)
    e_wt_ktr_m3g <- -0.701
    label("TBW power exponent on Ktr (M3G transit rate) (unitless)")   # Table 2 Final: M = -0.701 (RSE 30%)

    # ============================================================
    # M6G formation-delay transit chain
    # ============================================================
    lktr_m6g <- log(0.159)
    label("Morphine -> M6G transit-chain rate Ktr2 (1/min)")           # Table 2 Final: Ktr2 = 0.159 /min (n = 2 transits, MTT = 2/Ktr2 ~ 12.6 min; no TBW covariate)

    # ============================================================
    # Metabolite elimination
    # ============================================================
    lvc_m3g <- log(5.29)
    label("M3G central volume = M6G central volume (L; VM3G = VM6G structural equality)") # Table 2 Final: VM3G = VM6G = 5.29 L
    lcl_m3g <- log(0.134)
    label("M3G elimination clearance at TBW = 98.5 kg (L/min)")        # Table 2 Final: CLE M3G_98.5kg = 0.134 L/min
    e_wt_cl_m3g <- -1.08
    label("TBW power exponent on M3G elimination clearance (unitless)") # Table 2 Final: N = -1.08 (RSE 22%)
    lcl_m6g <- log(0.149)
    label("M6G elimination clearance at TBW = 98.5 kg (L/min)")        # Table 2 Final: CLE M6G_98.5kg = 0.149 L/min
    e_wt_cl_m6g <- -1.03
    label("TBW power exponent on M6G elimination clearance (unitless)") # Table 2 Final: O = -1.03 (RSE 31%)

    # ============================================================
    # IIV - log-normal; omega^2 = log(1 + CV^2). Reported CV% from
    # Table 2 Final model column.
    # ============================================================
    etalcl_form_m3g ~ log(1 + 0.208^2)
    # Table 2 Final IIV: CLF M3G = 20.8% CV -> omega^2 = log(1 + 0.208^2)
    etalcl_m3g ~ log(1 + 0.659^2)
    # Table 2 Final IIV: CLE M3G = 65.9% CV -> omega^2 = log(1 + 0.659^2)
    etalvc_m3g ~ log(1 + 0.297^2)
    # Table 2 Final IIV: VM3G = VM6G = 29.7% CV -> omega^2 = log(1 + 0.297^2); single ETA applied to the shared metabolite volume
    etalktr_m6g ~ log(1 + 0.368^2)
    # Table 2 Final IIV: Ktr2 = 36.8% CV -> omega^2 = log(1 + 0.368^2)

    # ============================================================
    # Residual error - proportional, separate magnitudes for the
    # healthy-volunteer (DIS_MORBOBESE = 0) and morbidly obese
    # (DIS_MORBOBESE = 1) cohorts (Table 2 Final model column).
    # The canonical residual-error parameters propSd, propSd_m3g, and
    # propSd_m6g are derived in model() as cohort-weighted
    # combinations of the four per-cohort SDs below.
    # ============================================================
    propSd_morphine_hv <- 0.140
    label("Morphine proportional residual SD in healthy volunteers (fraction)") # Table 2 Final: HV morphine 14.0%
    propSd_morphine_mo <- 0.379
    label("Morphine proportional residual SD in morbidly obese (fraction)")     # Table 2 Final: MO morphine 37.9%
    propSd_m3g_hv <- 0.179
    label("M3G proportional residual SD in healthy volunteers (fraction)")      # Table 2 Final: HV M3G 17.9%
    propSd_m3g_mo <- 0.171
    label("M3G proportional residual SD in morbidly obese (fraction)")          # Table 2 Final: MO M3G 17.1%
    propSd_m6g_hv <- 0.295
    label("M6G proportional residual SD in healthy volunteers (fraction)")      # Table 2 Final: HV M6G 29.5%
    propSd_m6g_mo <- 0.281
    label("M6G proportional residual SD in morbidly obese (fraction)")          # Table 2 Final: MO M6G 28.1%
  })

  model({
    # ------------------------------------------------------------
    # Body-weight covariate factors (reference TBW = 98.5 kg).
    # ------------------------------------------------------------
    wt_ratio <- WT / 98.5

    # ------------------------------------------------------------
    # Individual structural parameters.
    # ------------------------------------------------------------
    cl_form_m3g <- exp(lcl_form_m3g + etalcl_form_m3g)
    cl_form_m6g <- exp(lcl_form_m6g) * wt_ratio^e_wt_cl_form_m6g
    cl_nongluc  <- exp(lcl_nongluc)

    vc  <- exp(lvc)
    vp  <- exp(lvp)
    vp2 <- exp(lvp2) * wt_ratio^e_wt_vp2
    q   <- exp(lq)
    q2  <- exp(lq2)

    ktr_m3g <- exp(lktr_m3g) * wt_ratio^e_wt_ktr_m3g
    ktr_m6g <- exp(lktr_m6g + etalktr_m6g)

    vc_m3g <- exp(lvc_m3g + etalvc_m3g)
    vc_m6g <- vc_m3g                                   # VM3G = VM6G structural equality
    cl_m3g <- exp(lcl_m3g + etalcl_m3g) * wt_ratio^e_wt_cl_m3g
    cl_m6g <- exp(lcl_m6g) * wt_ratio^e_wt_cl_m6g

    # ------------------------------------------------------------
    # Internal state quantities are mg of morphine hydrochloride
    # equivalents in every compartment (morphine HCl MW 321.8 g/mol).
    # Concentrations Cc, Cc_m3g, Cc_m6g are output in nmol/L (the
    # paper's unit) by multiplying mg/L by the molar conversion
    # factor 1e6/321.8 in the observation equations. This is done at
    # observation time (not via f(central)) because rxode2's f()
    # rescales the infusion DURATION by the factor f -- using
    # f(central) > 1 for unit conversion breaks infusion dosing.
    # 1:1 molar stoichiometry across morphine HCl -> morphine free
    # base -> M3G -> M6G means the same molar conversion factor
    # applies to every output.
    # ------------------------------------------------------------
    mw_conv <- 1e6 / 321.8

    # ------------------------------------------------------------
    # Morphine three-compartment IV disposition. Total morphine
    # elimination = cl_nongluc + cl_form_m3g + cl_form_m6g (the two
    # formation clearances appear as both elimination from morphine
    # and source flux for the metabolite transit chains).
    # ------------------------------------------------------------
    d/dt(central) <- q  * peripheral1 / vp  + q2 * peripheral2 / vp2 -
                     (q + q2 + cl_nongluc + cl_form_m3g + cl_form_m6g) * central / vc
    d/dt(peripheral1) <- q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <- q2 * central / vc - q2 * peripheral2 / vp2

    # ------------------------------------------------------------
    # M3G formation-delay chain (5 transit compartments) feeding the
    # M3G central compartment. The formation flux entering transit1
    # is cl_form_m3g * Cc_morphine, in nmol/min of morphine
    # equivalents (= nmol/min of M3G by 1:1 molar stoichiometry).
    # ------------------------------------------------------------
    d/dt(transit1_m3g) <- cl_form_m3g * central / vc - ktr_m3g * transit1_m3g
    d/dt(transit2_m3g) <- ktr_m3g * transit1_m3g - ktr_m3g * transit2_m3g
    d/dt(transit3_m3g) <- ktr_m3g * transit2_m3g - ktr_m3g * transit3_m3g
    d/dt(transit4_m3g) <- ktr_m3g * transit3_m3g - ktr_m3g * transit4_m3g
    d/dt(transit5_m3g) <- ktr_m3g * transit4_m3g - ktr_m3g * transit5_m3g
    d/dt(central_m3g)  <- ktr_m3g * transit5_m3g - cl_m3g * central_m3g / vc_m3g

    # ------------------------------------------------------------
    # M6G formation-delay chain (2 transit compartments).
    # ------------------------------------------------------------
    d/dt(transit1_m6g) <- cl_form_m6g * central / vc - ktr_m6g * transit1_m6g
    d/dt(transit2_m6g) <- ktr_m6g * transit1_m6g - ktr_m6g * transit2_m6g
    d/dt(central_m6g)  <- ktr_m6g * transit2_m6g - cl_m6g * central_m6g / vc_m6g

    # ------------------------------------------------------------
    # Observations. Concentrations expressed in nmol/L (paper units)
    # via the molar conversion factor mw_conv = 1e6/321.8.
    # ------------------------------------------------------------
    Cc     <- (central     / vc)     * mw_conv
    Cc_m3g <- (central_m3g / vc_m3g) * mw_conv
    Cc_m6g <- (central_m6g / vc_m6g) * mw_conv

    # ------------------------------------------------------------
    # Cohort-conditional residual error. DIS_MORBOBESE = 1 selects
    # the morbidly obese magnitude, DIS_MORBOBESE = 0 selects the
    # healthy-volunteer magnitude (Table 2 Final, "Residual
    # variability (%)" block).
    # ------------------------------------------------------------
    propSd     <- propSd_morphine_hv * (1 - DIS_MORBOBESE) + propSd_morphine_mo * DIS_MORBOBESE
    propSd_m3g <- propSd_m3g_hv      * (1 - DIS_MORBOBESE) + propSd_m3g_mo      * DIS_MORBOBESE
    propSd_m6g <- propSd_m6g_hv      * (1 - DIS_MORBOBESE) + propSd_m6g_mo      * DIS_MORBOBESE

    Cc     ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
  })
}
