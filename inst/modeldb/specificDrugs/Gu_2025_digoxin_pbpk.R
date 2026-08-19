Gu_2025_digoxin_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 24-state, custom Phoenix WinNonlin 8.4",
    "implementation). digoxin disposition in healthy adults and in",
    "chronic heart failure across NYHA classes II, III and IV. digoxin is",
    "a cardiac glycoside eliminated primarily by the kidney (Supplementary Section 2.1).",
    "The circuit is lung, liver, kidney, spleen, heart, brain, muscle,",
    "skin, adipose and rest-of-body perfusion-limited tissues plus",
    "arterial and venous blood, a six-segment gastrointestinal lumen",
    "(stomach, duodenum, jejunum, ileum, caecum, colon) and the six",
    "matching gut-wall tissues, with first-order absorption from the",
    "duodenum, jejunum and ileum only. Heart failure is applied by",
    "rescaling organ blood flows (Table 1): splanchnic x 0.76 / 0.54 /",
    "0.46, renal x 0.78 / 0.55 / 0.63 and skin, adipose and muscle x",
    "0.57 / 0.44 / 0.28 for NYHA II / III / IV; intrinsic hepatic and",
    "renal clearances are unchanged. Deterministic: the paper virtual",
    "populations are uniform 80-120% draws on CLli,int, CLk,int, fu,b,",
    "Peff and ka, not lognormal random effects, so no IIV is encoded."
  )
  reference <- paste(
    "Gu W, Shao Q, Jiang L (2025). Predicting Pharmacokinetics of Drugs",
    "in Patients with Heart Failure and Optimizing Their Dosing",
    "Strategies Using a Physiologically Based Pharmacokinetic Model.",
    "Pharmaceutics 17(11):1394. doi:10.3390/pharmaceutics17111394.",
    "ODE system from Supplementary Material Eq S1-S15; system physiology",
    "from Supplementary Table S1; digoxin parameters from Supplementary",
    "Table S2; heart-failure blood flows from Table 1; validation targets",
    "from Table 2 (heart failure) and Supplementary Table S4 (healthy).",
    sep = " "
  )
  vignette <- "Gu_2025_heart_failure_pbpk"
  units    <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  # Gut-WALL tissue states and the caecum / colon lumen segments have no
  # canonical entry yet. wall_duodenum / wall_jejunum / wall_ileum follow
  # Luo_2024_benazepril_pbpk.R; wall_stomach, wall_cecum and wall_colon
  # extend that same wall_<segment> family, and cecum / colon extend the
  # stomach / duodenum / jejunum / ileum lumen family.
  paper_specific_compartments <- c(
    "wall_stomach", "wall_duodenum", "wall_jejunum", "wall_ileum",
    "wall_cecum", "wall_colon", "cecum", "colon"
  )

  covariateData <- list(
    DIS_CHF_NYHA2 = list(
      description        = "NYHA functional class II heart-failure indicator (1 = mild chronic heart failure).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes              = paste(
        "The three DIS_CHF_NYHA* indicators are mutually exclusive; all",
        "three 0 selects the healthy 70 kg physiology of Supplementary",
        "Table S1. DIS_CHF_NYHA2 = 1 selects the II Grade column of Gu",
        "2025 Table 1: splanchnic flows (hepatic artery, spleen, stomach",
        "and the five intestinal-wall segments) x 0.76, renal flow x 0.78,",
        "and skin, adipose and muscle flows x 0.57. Heart, brain and",
        "rest-of-body flows are unchanged and cardiac output is the sum",
        "of the rescaled flows. Intrinsic clearances are NOT rescaled",
        "(Section 2.3)."
      ),
      source_name        = "HF-II"
    ),
    DIS_CHF_NYHA3 = list(
      description        = "NYHA functional class III heart-failure indicator (1 = moderate chronic heart failure).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes              = paste(
        "The three DIS_CHF_NYHA* indicators are mutually exclusive; all",
        "three 0 selects the healthy 70 kg physiology of Supplementary",
        "Table S1. DIS_CHF_NYHA3 = 1 selects the III Grade column of Gu",
        "2025 Table 1: splanchnic flows (hepatic artery, spleen, stomach",
        "and the five intestinal-wall segments) x 0.54, renal flow x 0.55,",
        "and skin, adipose and muscle flows x 0.44. Heart, brain and",
        "rest-of-body flows are unchanged and cardiac output is the sum",
        "of the rescaled flows. Intrinsic clearances are NOT rescaled",
        "(Section 2.3)."
      ),
      source_name        = "HF-III"
    ),
    DIS_CHF_NYHA4 = list(
      description        = "NYHA functional class IV heart-failure indicator (1 = severe chronic heart failure).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes              = paste(
        "The three DIS_CHF_NYHA* indicators are mutually exclusive; all",
        "three 0 selects the healthy 70 kg physiology of Supplementary",
        "Table S1. DIS_CHF_NYHA4 = 1 selects the IV Grade column of Gu",
        "2025 Table 1: splanchnic flows (hepatic artery, spleen, stomach",
        "and the five intestinal-wall segments) x 0.46, renal flow x 0.63,",
        "and skin, adipose and muscle flows x 0.28. Heart, brain and",
        "rest-of-body flows are unchanged and cardiac output is the sum",
        "of the rescaled flows. Intrinsic clearances are NOT rescaled",
        "(Section 2.3)."
      ),
      source_name        = "HF-IV"
    )
  )

  compartmentData <- list(
    stomach        = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum       = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum        = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    ileum          = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    cecum          = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    colon          = list(analyte = "digoxin", units = "mg", specimen = "administration site", verified = TRUE),
    wall_stomach   = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    wall_duodenum  = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum   = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum     = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    wall_cecum     = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    wall_colon     = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    liver          = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    kidney         = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    spleen         = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    lung           = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    heart          = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    brain          = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    muscle         = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    skin           = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    adipose        = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    other          = list(analyte = "digoxin", units = "mg", specimen = "tissue", verified = TRUE),
    arterial       = list(analyte = "digoxin", units = "mg", specimen = "whole blood", verified = TRUE),
    venous         = list(analyte = "digoxin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 115L,
    n_studies      = 9L,
    age_range      = "adults (18-90 years across the pooled reports)",
    weight_median  = "70 kg (the model is parameterised for a 70 kg adult; no weight covariate)",
    disease_state  = paste(
      "Pooled healthy volunteers and chronic heart-failure patients",
      "(NYHA classes II-IV) from the literature reports tabulated in Gu",
      "2025 Supplementary Table S3."
    ),
    dose_range     = "0.1-1 mg single oral (healthy 0.25-1 mg; HF 0.1 and 0.25 mg)",
    notes          = paste(
      "Literature-digitised clinical data; subject counts are the sum of",
      "the per-report n in Supplementary Table S3. The authors simulated",
      "1000 virtual individuals per population by drawing CLli,int,",
      "CLk,int, fu,b, Peff and ka uniformly over 80-120% of their",
      "baseline values (Section 2.3), so the reported 5th-95th percentile",
      "bands are a parameter-uncertainty sweep rather than estimated",
      "between-subject variability."
    )
  )

  ini({
    # ---- Physicochemical / binding properties (Supplementary Table S2) ----
    fu_p       <- fixed(0.71);    label("Fraction unbound in plasma (unitless)")  # Table S2 fu,plasma [6] Neuhoff 2013
    bpr        <- fixed(1);    label("Blood:plasma concentration ratio Rb (unitless)")  # Table S2 Rbp assumed (footnote c)
    # Table S2 gives fu,g = 1 assumed (footnote c), but digoxin has no
    # tabulated intestinal-wall intrinsic clearance (CLGWi,int is '/' in
    # Table S2), so Eq S14 carries no gut-metabolism term and fu,g never
    # enters the model. It is therefore not declared as a parameter.

    # ---- Organ clearances (Supplementary Table S2), mL/min -> L/min. ----
    # These are the literature PLASMA clearances; Eq S9 converts them to
    # blood clearances and Eq S7 / S11 invert the well-stirred relation to
    # the intrinsic clearances the ODEs use.
    lcl_renal  <- fixed(log(0.0984));  label("Renal plasma clearance CLrenal (log L/min)")  # Table S2 CLrenal 98.4 mL/min [29] Tsutsumi 2002
    lcl_liver  <- fixed(log(0.0416));  label("Hepatic plasma clearance CLliver (log L/min)")  # Table S2 CLliver 41.6 mL/min [29] Tsutsumi 2002

    # ---- Absorption ----
    # Peff is tabulated, so Eq S15 (ka,i = 2 x Peff / ri) gives a distinct
    # absorption rate constant per intestinal segment.
    lpeff      <- fixed(log(0.0282));  label("Effective permeability coefficient Peff,A-B (log cm/min)")  # Table S2 Peff,A-B 0.0282 cm/min [6] Neuhoff 2013

    # ---- Tissue:plasma partition coefficients (Supplementary Table S2) ----
    # footnote e (method of [37]); KLiver:p and KMuscle:p from [24] / [36]
    kp_adipose   <- fixed(142.15  ); label("Adipose:plasma partition coefficient (unitless)")  # Table S2 KAdipose:plasma
    kp_liver     <- fixed(10.83   ); label("Liver:plasma partition coefficient (unitless)")  # Table S2 KLiver:plasma
    kp_muscle    <- fixed(7.35    ); label("Muscle:plasma partition coefficient (unitless)")  # Table S2 KMuscle:plasma
    kp_lung      <- fixed(4.83    ); label("Lung:plasma partition coefficient (unitless)")  # Table S2 KLung:plasma
    kp_kidney    <- fixed(3.73    ); label("Kidney:plasma partition coefficient (unitless)")  # Table S2 KKidney:plasma
    kp_stomach   <- fixed(2.55    ); label("Stomach:plasma partition coefficient (unitless)")  # Table S2 Kstomach:plasma
    kp_other     <- fixed(0.01    ); label("Rest-of-body:plasma partition coefficient (unitless)")  # Table S2 KROB:plasma
    kp_brain     <- fixed(7.11    ); label("Brain:plasma partition coefficient (unitless)")  # Table S2 KBrain:plasma
    kp_heart     <- fixed(3.42    ); label("Heart:plasma partition coefficient (unitless)")  # Table S2 KHeart:plasma
    kp_spleen    <- fixed(2.36    ); label("Spleen:plasma partition coefficient (unitless)")  # Table S2 KSpleen:plasma
    kp_skin      <- fixed(10.64   ); label("Skin:plasma partition coefficient (unitless)")  # Table S2 KSkin:plasma
    kp_gut       <- fixed(7.88    ); label("Intestinal-wall:plasma partition coefficient (unitless)")  # Table S2 Kgut:plasma

    # ---- Residual error ----
    # Gu 2025 report no residual-error model: the PBPK model is a forward
    # predictor assessed by 0.5-2.0-fold AUC / Cmax agreement and 5th-95th
    # percentile coverage (Section 2.4), not by a fitted sigma. Non-fitted
    # placeholder per the in-repo PBPK convention (Luo_2024_benazepril_pbpk.R).
    propSd     <- fixed(0.10); label("Proportional residual error placeholder (fraction)")  # not reported by Gu 2025
  })

  model({
    # =========================================================
    # 1. NYHA decode. The three indicators are mutually exclusive;
    #    hf_none = 1 selects the healthy Table S1 physiology.
    # =========================================================
    hf_2    <- DIS_CHF_NYHA2
    hf_3    <- DIS_CHF_NYHA3
    hf_4    <- DIS_CHF_NYHA4
    hf_none <- 1 - hf_2 - hf_3 - hf_4

    # Table 1 blood-flow multipliers relative to Table S1 (Section 2.3).
    f_splanchnic <- hf_none * 1 + hf_2 * 0.76 + hf_3 * 0.54 + hf_4 * 0.46
    f_renal      <- hf_none * 1 + hf_2 * 0.78 + hf_3 * 0.55 + hf_4 * 0.63
    f_periph     <- hf_none * 1 + hf_2 * 0.57 + hf_3 * 0.44 + hf_4 * 0.28

    # =========================================================
    # 2. Tissue volumes, Supplementary Table S1 (mL -> L). Unchanged by
    #    heart failure. The 18 tissue volumes sum to exactly 70 L.
    # =========================================================
    v_lung      <- 1.17     # Table S1 1170 mL
    v_kidney    <- 0.28     # Table S1 280 mL
    v_heart     <- 0.31     # Table S1 310 mL
    v_liver     <- 1.69     # Table S1 1690 mL
    v_muscle    <- 35       # Table S1 35000 mL
    v_skin      <- 7.8      # Table S1 7800 mL
    v_brain     <- 1.45     # Table S1 1450 mL
    v_adipose   <- 10       # Table S1 10000 mL
    v_other     <- 5.1      # Table S1 5100 mL
    v_spleen    <- 0.19     # Table S1 190 mL
    v_venous    <- 3.47     # Table S1 3470 mL
    v_arterial  <- 1.73     # Table S1 1730 mL
    v_stomach   <- 0.16     # Table S1 160 mL
    v_duodenum  <- 0.07     # Table S1 70 mL
    v_jejunum   <- 0.209    # Table S1 209 mL
    v_ileum     <- 0.139    # Table S1 139 mL
    v_cecum     <- 0.116    # Table S1 116 mL
    v_colon     <- 1.116    # Table S1 1116 mL

    # =========================================================
    # 3. Blood flows, Supplementary Table S1 (mL/min -> L/min) rescaled by
    #    the Table 1 heart-failure multipliers. Reproduces every entry of
    #    Table 1 exactly, including the cardiac-output (Lung) row.
    # =========================================================
    q_hepatic_artery <- 0.3 * f_splanchnic     # Table S1 hepatic arterial 300 mL/min (QLiver 1518 minus the splanchnic inflows)
    q_stomach        <- 0.038  * f_splanchnic  # Table S1 38 mL/min
    q_duodenum       <- 0.118  * f_splanchnic  # Table S1 118 mL/min
    q_jejunum        <- 0.413  * f_splanchnic  # Table S1 413 mL/min
    q_ileum          <- 0.244  * f_splanchnic  # Table S1 244 mL/min
    q_cecum          <- 0.044  * f_splanchnic  # Table S1 44 mL/min
    q_colon          <- 0.281  * f_splanchnic  # Table S1 281 mL/min
    q_spleen         <- 0.08   * f_splanchnic  # Table S1 80 mL/min
    q_kidney         <- 1.24 * f_renal        # Table S1 1240 mL/min
    q_muscle         <- 0.75 * f_periph       # Table S1 750 mL/min
    q_skin           <- 0.3 * f_periph        # Table S1 300 mL/min
    q_adipose        <- 0.26 * f_periph       # Table S1 260 mL/min
    q_heart          <- 0.24                  # Table S1 240 mL/min (unchanged, Table 1 footnote a)
    q_brain          <- 0.7                   # Table S1 700 mL/min (unchanged, Table 1 footnote a)
    q_other          <- 0.592                 # Table S1 592 mL/min (unchanged, Table 1 footnote a)

    # Eq S5: QLi = QLa + QSt + QSp + sum(QGwi).
    q_liver <- q_hepatic_artery + q_stomach + q_spleen + q_duodenum +
               q_jejunum + q_ileum + q_cecum + q_colon
    # Cardiac output = sum of the arterial outflows (5600 mL/min healthy).
    q_total <- q_liver + q_kidney + q_muscle + q_skin + q_adipose +
               q_heart + q_brain + q_other

    # =========================================================
    # 4. Blood binding and intrinsic clearances.
    #    Eq S6:  fu,b   = fu,p / Rb
    #    Eq S9:  CL,b   = CL,p / (1 - Hct + Rb x Hct),  Hct = 0.43
    #    Eq S7:  CLli,b = QLi x fu,b x CLli,int / (QLi + fu,b x CLli,int)
    #    Eq S11: CLK    = QK  x fu,b x CLK,int  / (QK  + fu,b x CLK,int)
    #    S7 / S11 are inverted here to recover the intrinsic clearances.
    #    The inversion uses the HEALTHY reference flows because heart
    #    failure does not alter intrinsic clearance (Section 2.3); only the
    #    q_liver / q_kidney in the ODEs below carry the HF rescaling.
    # =========================================================
    hct    <- 0.43                            # Supplementary text after Eq S9 [1]
    fu_b   <- fu_p / bpr
    q_liver_ref  <- 1.518                     # Table S1 healthy hepatic flow 1518 mL/min
    q_kidney_ref <- 1.24                      # Table S1 healthy renal flow 1240 mL/min
    cl_liver_b <- exp(lcl_liver) / (1 - hct + bpr * hct)
    cl_renal_b <- exp(lcl_renal) / (1 - hct + bpr * hct)
    cl_liver_int <- q_liver_ref  * cl_liver_b / ((q_liver_ref  - cl_liver_b) * fu_b)
    cl_renal_int <- q_kidney_ref * cl_renal_b / ((q_kidney_ref - cl_renal_b) * fu_b)

    # =========================================================
    # 5. Segment absorption rate constants.
    # =========================================================
    # Eq S15: ka,i = 2 x Peff,A-B / ri, with the Table S1 intestinal radii.
    ka_duodenum <- 2 * exp(lpeff) / 2      # Table S1 duodenal radius 2 cm
    ka_jejunum  <- 2 * exp(lpeff) / 1.63   # Table S1 jejunal radius 1.63 cm
    ka_ileum    <- 2 * exp(lpeff) / 1.45   # Table S1 ileal radius 1.45 cm

    # Gastric-emptying and intestinal-transit rate constants, Table S1.
    kt_stomach   <- 0.0462    # Table S1 0.0462 1/min
    kt_duodenum  <- 0.0462    # Table S1 0.0462 1/min
    kt_jejunum   <- 0.012     # Table S1 0.012 1/min
    kt_ileum     <- 0.0058    # Table S1 0.0058 1/min
    kt_cecum     <- 0.0025    # Table S1 0.0025 1/min
    kt_colon     <- 0.00085   # Table S1 0.00085 1/min

    # =========================================================
    # 6. Tissue concentrations (mg/L = ug/mL) and the blood-side effluent
    #    concentration Ct / (Kt:p / Rb) that leaves each tissue.
    # =========================================================
    c_liver          <- liver / v_liver
    c_kidney         <- kidney / v_kidney
    c_spleen         <- spleen / v_spleen
    c_lung           <- lung / v_lung
    c_heart          <- heart / v_heart
    c_brain          <- brain / v_brain
    c_muscle         <- muscle / v_muscle
    c_skin           <- skin / v_skin
    c_adipose        <- adipose / v_adipose
    c_other          <- other / v_other
    c_arterial       <- arterial / v_arterial
    c_venous         <- venous / v_venous
    c_wall_stomach   <- wall_stomach / v_stomach
    c_wall_duodenum  <- wall_duodenum / v_duodenum
    c_wall_jejunum   <- wall_jejunum / v_jejunum
    c_wall_ileum     <- wall_ileum / v_ileum
    c_wall_cecum     <- wall_cecum / v_cecum
    c_wall_colon     <- wall_colon / v_colon

    eff_liver          <- c_liver * bpr / kp_liver
    eff_kidney         <- c_kidney * bpr / kp_kidney
    eff_spleen         <- c_spleen * bpr / kp_spleen
    eff_lung           <- c_lung * bpr / kp_lung
    eff_heart          <- c_heart * bpr / kp_heart
    eff_brain          <- c_brain * bpr / kp_brain
    eff_muscle         <- c_muscle * bpr / kp_muscle
    eff_skin           <- c_skin * bpr / kp_skin
    eff_adipose        <- c_adipose * bpr / kp_adipose
    eff_other          <- c_other * bpr / kp_other
    eff_wall_stomach   <- c_wall_stomach * bpr / kp_stomach
    eff_wall_duodenum  <- c_wall_duodenum * bpr / kp_gut
    eff_wall_jejunum   <- c_wall_jejunum * bpr / kp_gut
    eff_wall_ileum     <- c_wall_ileum * bpr / kp_gut
    eff_wall_cecum     <- c_wall_cecum * bpr / kp_gut
    eff_wall_colon     <- c_wall_colon * bpr / kp_gut

    # =========================================================
    # 7. Gastrointestinal lumen. Eq S12 (stomach) and Eq S13 (segments).
    #    Eq S13 is printed as dAi/dt = Kt,i-1 Ai-1 - (Kt,i-1 + ka,i) Ai, i.e.
    #    with the PRECEDING segment's transit constant on the efflux term.
    #    That is a typographical slip: it leaves the colon transit constant
    #    of Table S1 unused and shifts every other one by a segment. The
    #    standard compartmental-absorption-and-transit form of Yu & Amidon
    #    (Supplementary reference [2]), with each segment effluxing at its
    #    OWN Kt,i, is used instead -- it reproduces the Table S4 / Table 2
    #    predictions to a 1% median error versus 26% for the printed form.
    #    See the vignette's Assumptions and deviations section.
    # =========================================================
    d/dt(stomach)  <- -kt_stomach * stomach
    d/dt(duodenum) <-  kt_stomach  * stomach  - kt_duodenum * duodenum - ka_duodenum * duodenum
    d/dt(jejunum)  <-  kt_duodenum * duodenum - kt_jejunum  * jejunum  - ka_jejunum  * jejunum
    d/dt(ileum)    <-  kt_jejunum  * jejunum  - kt_ileum    * ileum    - ka_ileum    * ileum
    d/dt(cecum)    <-  kt_ileum    * ileum    - kt_cecum    * cecum
    d/dt(colon)    <-  kt_cecum    * cecum    - kt_colon    * colon

    # =========================================================
    # 8. Gut-wall tissues, Eq S14. Absorption enters the duodenal, jejunal
    #    and ileal walls only; every wall drains into the liver (Eq S5).
    # =========================================================
    d/dt(wall_stomach) <- q_stomach * c_arterial - q_stomach * eff_wall_stomach
    d/dt(wall_duodenum) <- q_duodenum * c_arterial + ka_duodenum * duodenum - q_duodenum * eff_wall_duodenum
    d/dt(wall_jejunum) <- q_jejunum * c_arterial + ka_jejunum * jejunum - q_jejunum * eff_wall_jejunum
    d/dt(wall_ileum) <- q_ileum * c_arterial + ka_ileum * ileum - q_ileum * eff_wall_ileum
    d/dt(wall_cecum) <- q_cecum * c_arterial - q_cecum * eff_wall_cecum
    d/dt(wall_colon) <- q_colon * c_arterial - q_colon * eff_wall_colon

    # =========================================================
    # 9. Perfusion-limited tissues, Eq S1; spleen also drains to liver.
    # =========================================================
    d/dt(heart   ) <- q_heart * (c_arterial - eff_heart)
    d/dt(brain   ) <- q_brain * (c_arterial - eff_brain)
    d/dt(muscle  ) <- q_muscle * (c_arterial - eff_muscle)
    d/dt(skin    ) <- q_skin * (c_arterial - eff_skin)
    d/dt(adipose ) <- q_adipose * (c_arterial - eff_adipose)
    d/dt(other   ) <- q_other * (c_arterial - eff_other)
    d/dt(spleen  ) <- q_spleen * (c_arterial - eff_spleen)

    # Kidney, Eq S10.
    d/dt(kidney) <- q_kidney * c_arterial -
      (q_kidney + fu_b * cl_renal_int) * eff_kidney

    # Liver, Eq S5: hepatic artery plus stomach, spleen and all five
    # intestinal-wall effluents in, well-stirred hepatic extraction out.
    d/dt(liver) <- q_hepatic_artery * c_arterial + q_stomach * eff_wall_stomach +
      q_spleen * eff_spleen + q_duodenum * eff_wall_duodenum +
      q_jejunum * eff_wall_jejunum + q_ileum * eff_wall_ileum +
      q_cecum * eff_wall_cecum + q_colon * eff_wall_colon -
      (q_liver + fu_b * cl_liver_int) * eff_liver

    # Venous blood, Eq S2; lung, Eq S3; arterial blood, Eq S4.
    d/dt(venous) <- q_liver * eff_liver + q_kidney * eff_kidney +
      q_heart * eff_heart + q_brain * eff_brain + q_muscle * eff_muscle +
      q_skin * eff_skin + q_adipose * eff_adipose + q_other * eff_other -
      q_total * c_venous
    d/dt(lung)     <- q_total * (c_venous - eff_lung)
    d/dt(arterial) <- q_total * (eff_lung - c_arterial)

    # =========================================================
    # 10. Observation. States are mg and volumes L, so concentrations are
    #     mg/L = ug/mL. c_venous is a BLOOD concentration, so dividing by
    #     Rb gives the venous PLASMA concentration the clinical reports and
    #     Table 2 / Table S4 are expressed in.
    # =========================================================
    Cc <- c_venous / bpr

    Cc ~ prop(propSd)
  })
}

