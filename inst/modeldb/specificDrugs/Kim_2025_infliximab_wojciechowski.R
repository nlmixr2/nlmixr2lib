Kim_2025_infliximab_wojciechowski <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in patients with inflammatory bowel disease (Wojciechowski 2017 model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Wojciechowski J, Upton RN, Mould DR, Wiese MD, Foster DJR. Infliximab",
    "maintenance dosing in inflammatory bowel disease: an example for in silico",
    "assessment of adaptive dosing strategies. AAPS J. 2017;19(4):1136-1147.",
    "doi:10.1208/s12248-017-0082-8 (reference 24 of Kim 2025).",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #10) with MAXEVAL=0,",
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
      notes              = "Power effect on all four disposition parameters, normalized to a 70 kg reference: CL *= (WT/70)^0.614, Vc *= (WT/70)^0.691, Vp *= (WT/70)^0.59 and Q *= (WT/70)^1.1. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement).",
      source_name        = "WGT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, normalized to a 4 g/dL reference: CL *= (ALB/4)^(-1.17), i.e. lower albumin gives higher clearance. The source control stream states the covariate is in g/dL; the canonical ALB column is in SI g/L, so model() converts inline (g/L * 0.1 = g/dL). Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies toward infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Exponential effect on CL: CL *= exp(0.257 * ADA_POS), i.e. a 1.293-fold (+29.3%) increase in clearance when ADA-positive. Note this is an exponential form, unlike the linear (1 + theta) form used by the Dubinsky and Fasanmade models. Source paper labels this covariate 'ATI' (antibodies toward infliximab); renamed to canonical ADA_POS per covariate-columns.md. Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ATI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 655L,
    n_studies      = 1L,
    age_range      = "Pooled paediatric and adult cohort; age not specified in Kim 2025 Table S1 for this model.",
    weight_range   = "Not reported; median body weight 70 kg (Kim 2025 Table S1), which is the reference weight used by the model.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not specified.",
    disease_state  = "Inflammatory bowel disease: Crohn's disease (n = 112) and ulcerative colitis (n = 543).",
    dose_range     = "Intravenous infliximab; treatment phase and sampling times not specified in Kim 2025 Table S1.",
    regions        = "Not specified.",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Wojciechowski model: Crohn's disease (n = 112) and",
      "ulcerative colitis (n = 543), pooled paediatric and adult patients,",
      "median body weight 70 kg, median albumin 4.0 g/dL. Treatment phase,",
      "sampling time, ATI positivity, age, sex and immunomodulator use were not",
      "specified.",
      "Although the development cohort was predominantly adult, Kim 2025",
      "evaluated this model as one of the seven paediatric models because that",
      "is how it is deployed in TDMx.",
      "The external-validation cohort in Kim 2025 itself (57 Korean paediatric",
      "IBD patients with 158 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model had the lowest bias of the paediatric models (MPE 30.4%,",
      "MAE 1.8 mg/L, rRMSE 150.0%; Kim 2025 Table 2) but still did not reach",
      "clinical acceptability. Kim 2025 noted that this model and the Dubinsky",
      "model carry the largest inter-individual variability of the paediatric",
      "models evaluated (approximately 110% on Q).",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- typical values for the reference subject
    # (WT = 70 kg, ALB = 4 g/dL, ADA-negative).
    lcl <- log(0.294);  label("Typical clearance for the reference subject (CL, L/day)")                    # Data S1 model #10 $THETA1 = 0.294 L/day; Table S2 "0.294"
    lvc <- log(3.33);   label("Typical central volume for the reference subject (Vc, L)")                   # Data S1 model #10 $THETA2 = 3.33 L; Table S2 "3.33"
    lvp <- log(1.14);   label("Typical peripheral volume for the reference subject (Vp, L)")                # Data S1 model #10 $THETA3 = 1.14 L; Table S2 "1.14"
    lq  <- log(0.0719); label("Typical inter-compartmental clearance for the reference subject (Q, L/day)") # Data S1 model #10 $THETA4 = 0.0719 L/day; Table S2 "0.0719"

    # Covariate effects
    e_wt_cl  <-  0.614; label("Power exponent of body weight on CL ((WT/70)^e_wt_cl)")                     # Data S1 model #10 $THETA5 = 0.614
    e_alb_cl <- -1.17;  label("Power exponent of serum albumin on CL ((ALB/4)^e_alb_cl)")                  # Data S1 model #10 $THETA6 = -1.17
    e_ada_cl <-  0.257; label("Log-scale ADA-positive effect on CL (multiplier exp(e_ada_cl * ADA_POS))")  # Data S1 model #10 $THETA7 = 0.257
    e_wt_vc  <-  0.691; label("Power exponent of body weight on Vc ((WT/70)^e_wt_vc)")                     # Data S1 model #10 $THETA8 = 0.691
    e_wt_vp  <-  0.59;  label("Power exponent of body weight on Vp ((WT/70)^e_wt_vp)")                     # Data S1 model #10 $THETA9 = 0.59
    e_wt_q   <-  1.1;   label("Power exponent of body weight on Q ((WT/70)^e_wt_q)")                       # Data S1 model #10 $THETA10 = 1.1

    # Inter-individual variability; OMEGA = (CV)^2 per Kim 2025's control-stream
    # annotation, with each %CV given in the stream comment and Table S2.
    etalcl ~ 0.107   # Data S1 model #10 $OMEGA ETA1 = 0.107 (32.7% CV); Table S2 '(32.7%)'
    etalvc ~ 0.023   # Data S1 model #10 $OMEGA ETA2 = 0.023 (15% CV); Table S2 '(15.2%)'
    etalvp ~ 0.638   # Data S1 model #10 $OMEGA ETA3 = 0.638 (79.9% CV); Table S2 '(79.9%)'
    etalq  ~ 1.210   # Data S1 model #10 $OMEGA ETA4 = 1.210 (110% CV); Table S2 '(110%)'

    # Residual error -- proportional only.
    propSd <- 0.419; label("Proportional residual error (fraction)")  # Data S1 model #10 $SIGMA 0.176, comment "proportional error (41.9 %CV)"
  })
  model({
    # Canonical ALB is in SI g/L; this model was calibrated on g/dL.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL

    # Individual PK parameters for the reference subject WT = 70 kg,
    # alb_gdL = 4 g/dL, ADA-negative:
    #   CL (L/day) = 0.294  * (WT/70)^0.614 * (alb_gdL/4)^(-1.17) * exp(0.257*ADA)
    #   Vc (L)     = 3.33   * (WT/70)^0.691
    #   Vp (L)     = 1.14   * (WT/70)^0.59
    #   Q  (L/day) = 0.0719 * (WT/70)^1.1
    cl <- exp(lcl + e_ada_cl * ADA_POS + etalcl) * (WT / 70)^e_wt_cl *
      (alb_gdL / 4)^e_alb_cl
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
