Kim_2025_infliximab_dubinsky <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in children with inflammatory bowel disease (Dubinsky 2017 model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Dubinsky MC, Phan BL, Singh N, Rabizadeh S, Mould DR. Pharmacokinetic",
    "dashboard-recommended dosing is different than standard of care dosing in",
    "infliximab-treated pediatric IBD patients. AAPS J. 2017;19(1):215-222.",
    "doi:10.1208/s12248-016-9994-y (reference 21 of Kim 2025).",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #5) with MAXEVAL=0,",
    "i.e. every THETA / OMEGA / SIGMA is held at the originally published value.",
    "Values here are transcribed from that control stream and cross-checked",
    "against Kim 2025 Table S2.",
    sep = " "
  )
  vignette <- "Kim_2025_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on all four disposition parameters, normalized to a 70 kg reference: CL *= (WT/70)^0.612, Vc *= (WT/70)^0.696, Vp *= (WT/70)^0.604 and Q *= (WT/70)^1.15. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement) and identified this as a probable contributor to the poor predictive performance of the paediatric models.",
      source_name        = "WGT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, normalized to a 4 g/dL reference: CL *= (ALB/4)^(-2.3), i.e. lower albumin gives markedly higher clearance. The source control stream states the covariate is in g/dL; the canonical ALB column is in SI g/L, so model() converts inline (g/L * 0.1 = g/dL). Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies toward infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Linear fractional effect on CL: CL *= (1 + 0.231 * ADA_POS), i.e. +23.1% clearance when ADA-positive. Source paper labels this covariate 'ATI' (antibodies toward infliximab); renamed to canonical ADA_POS per covariate-columns.md. Time-varying in the Kim 2025 validation dataset. Kim 2025 showed (Figure 3) that assuming every patient to be ATI-positive improved predictive performance for the five ATI-carrying models, which they attributed to a positive bias inherent in the population PK models.",
      source_name        = "ATI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "Paediatric; mean age 13 years (Kim 2025 Table S1).",
    weight_range   = "Not reported; mean body weight 41 kg (Kim 2025 Table S1). The model's reference weight is 70 kg, well above the development cohort's mean.",
    sex_female_pct = 44,
    race_ethnicity = "Not specified.",
    disease_state  = "Inflammatory bowel disease: Crohn's disease (n = 41) and ulcerative colitis (n = 9).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "United States.",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Dubinsky model: Crohn's disease (n = 41) and ulcerative",
      "colitis (n = 9), paediatric patients, induction and maintenance phases,",
      "trough sampling, 28% ATI positive, mean age 13 years, 44% female, mean",
      "body weight 41 kg, mean albumin 3.9 g/dL, 22% on an immunomodulator.",
      "Mean values are reported for this model rather than medians.",
      "The external-validation cohort in Kim 2025 itself (57 Korean paediatric",
      "IBD patients with 158 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model gave MPE 31.7%, MAE 1.8 mg/L and rRMSE 198.1%",
      "(Kim 2025 Table 2). Kim 2025 noted that this model and the Wojciechowski",
      "model carry the largest inter-individual variability of the paediatric",
      "models evaluated (approximately 110% on Q) and cautioned that such large",
      "IIV can mask structural model deficiencies while still failing to",
      "improve future-concentration prediction.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- typical values for the reference subject
    # (WT = 70 kg, ALB = 4 g/dL, ADA-negative).
    lcl <- log(0.296);  label("Typical clearance for the reference subject (CL, L/day)")                    # Data S1 model #5 $THETA1 = 0.296 L/day; Table S2 "0.296"
    lvc <- log(3.30);   label("Typical central volume for the reference subject (Vc, L)")                   # Data S1 model #5 $THETA2 = 3.30 L; Table S2 "3.3"
    lvp <- log(1.16);   label("Typical peripheral volume for the reference subject (Vp, L)")                # Data S1 model #5 $THETA3 = 1.16 L; Table S2 "1.16"
    lq  <- log(0.0781); label("Typical inter-compartmental clearance for the reference subject (Q, L/day)") # Data S1 model #5 $THETA4 = 0.0781 L/day; Table S2 "0.0781"

    # Covariate effects
    e_wt_cl  <-  0.612; label("Power exponent of body weight on CL ((WT/70)^e_wt_cl)")             # Data S1 model #5 $THETA5 = 0.612
    e_alb_cl <- -2.3;   label("Power exponent of serum albumin on CL ((ALB/4)^e_alb_cl)")          # Data S1 model #5 $THETA6 = -2.3
    e_ada_cl <-  0.231; label("Fractional increase in CL when ADA-positive (1 + e_ada_cl)")        # Data S1 model #5 $THETA7 = 0.231
    e_wt_vc  <-  0.696; label("Power exponent of body weight on Vc ((WT/70)^e_wt_vc)")             # Data S1 model #5 $THETA8 = 0.696
    e_wt_vp  <-  0.604; label("Power exponent of body weight on Vp ((WT/70)^e_wt_vp)")             # Data S1 model #5 $THETA9 = 0.604
    e_wt_q   <-  1.15;  label("Power exponent of body weight on Q ((WT/70)^e_wt_q)")               # Data S1 model #5 $THETA10 = 1.15

    # Inter-individual variability. Kim 2025's control streams use
    # OMEGA = (CV)^2 (see the annotation on the Aubourg stream); each value's
    # %CV is given in the stream comment and reproduced in Table S2.
    etalcl ~ 0.098   # Data S1 model #5 $OMEGA ETA1 = 0.098 (31.3% CV); Table S2 '(31.3%)'
    etalvc ~ 0.010   # Data S1 model #5 $OMEGA ETA2 = 0.010 (9.85% CV); Table S2 '(9.85%)'
    etalvp ~ 0.579   # Data S1 model #5 $OMEGA ETA3 = 0.579 (76.1% CV); Table S2 '(76.1%)'
    etalq  ~ 1.232   # Data S1 model #5 $OMEGA ETA4 = 1.232 (111% CV); Table S2 '(111%)'

    # Residual error -- proportional only.
    propSd <- 0.419; label("Proportional residual error (fraction)")  # Data S1 model #5 $SIGMA 0.176, comment "proportional error (41.9 %CV)"
  })
  model({
    # Canonical ALB is in SI g/L; this model was calibrated on g/dL.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL

    # Individual PK parameters for the reference subject WT = 70 kg,
    # alb_gdL = 4 g/dL, ADA-negative:
    #   CL (L/day) = 0.296  * (WT/70)^0.612 * (alb_gdL/4)^(-2.3) * (1 + 0.231*ADA)
    #   Vc (L)     = 3.30   * (WT/70)^0.696
    #   Vp (L)     = 1.16   * (WT/70)^0.604
    #   Q  (L/day) = 0.0781 * (WT/70)^1.15
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (alb_gdL / 4)^e_alb_cl *
      (1 + e_ada_cl * ADA_POS)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model; infliximab is given as an IV infusion into central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L, numerically identical to ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
