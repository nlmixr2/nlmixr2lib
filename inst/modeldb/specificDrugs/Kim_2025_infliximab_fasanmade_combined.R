Kim_2025_infliximab_fasanmade_combined <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in children and adults with Crohn's disease, with inter-occasion variability on clearance (Fasanmade 2011 pooled REACH + ACCENT I model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Fasanmade AA, Adedokun OJ, Blank M, Zhou H, Davis HM. Pharmacokinetic",
    "properties of infliximab in children and adults with Crohn's disease: a",
    "retrospective analysis of data from 2 phase III clinical trials.",
    "Clin Ther. 2011;33(7):946-964. doi:10.1016/j.clinthera.2011.06.002",
    "(reference 22 of Kim 2025), pooling the paediatric REACH and adult",
    "ACCENT I trials.",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #6) with MAXEVAL=0,",
    "i.e. every THETA / OMEGA / SIGMA is held at the originally published value.",
    "Values here are transcribed from that control stream and cross-checked",
    "against Kim 2025 Table S2.",
    "NOTE: modellib('Frymoyer_2017_infliximab') is an independent transcription",
    "of this same underlying Fasanmade model, taken from Frymoyer 2017 rather",
    "than from Kim 2025. The two transcriptions agree on every typical value and",
    "on the albumin, ADA and immunomodulator effects, but disagree on the SIGN",
    "of the three body-weight exponents; see the e_wt_* comments in ini() below.",
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
      notes              = "The source equations are written on per-kilogram CL / V values with reference 65 kg, so the total-parameter exponent is 1 plus the published per-kg exponent. Kim 2025's control stream gives per-kg exponents of -0.313 (CL), -0.233 (Vc) and -0.588 (Vp), i.e. total-parameter exponents of 0.687, 0.767 and 0.412; Q is a constant per-kg value, giving a total-Q exponent of 1.0. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement).",
      source_name        = "WGT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, normalized to a 4.1 g/dL reference: CL *= (ALB/4.1)^(-0.855). The source control stream states the covariate is in g/dL; the canonical ALB column is in SI g/L, so model() converts inline (g/L * 0.1 = g/dL). Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies toward infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Linear fractional effect on CL: CL *= (1 + 0.291 * ADA_POS), i.e. +29.1% clearance when ADA-positive. Source paper labels this covariate 'ATI' (antibodies toward infliximab); renamed to canonical ADA_POS per covariate-columns.md. Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ATI"
    ),
    CONMED_IMMUNOMOD = list(
      description        = "Concomitant immunomodulator therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant immunomodulator)",
      notes              = "Linear fractional effect on CL: CL *= (1 - 0.137 * CONMED_IMMUNOMOD), i.e. -13.7% clearance while on an immunomodulator. In Kim 2025's validation dataset the immunomodulator class comprised azathioprine, mercaptopurine, methotrexate, cyclosporine and tacrolimus (Table 1 footnote). Source paper labels this covariate 'IMM'; renamed to canonical CONMED_IMMUNOMOD per covariate-columns.md. Time-varying in the Kim 2025 validation dataset.",
      source_name        = "IMM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 692L,
    n_studies      = 2L,
    age_range      = "Pooled paediatric and adult cohort; median age 33 years (Kim 2025 Table S1).",
    weight_range   = "Not reported; median body weight 65 kg (Kim 2025 Table S1), which is the reference weight used by the model.",
    sex_female_pct = 56,
    race_ethnicity = "Not specified.",
    disease_state  = "Crohn's disease (n = 692).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "Multi-regional pivotal trials of infliximab in Crohn's disease (REACH and ACCENT I).",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Fasanmade_combined model: Crohn's disease (n = 692),",
      "pooled paediatric and adult patients, induction and maintenance phases,",
      "peak / intermediate / trough sampling, 10% ATI positive, median age",
      "33 years, 56% female, median body weight 65 kg, median albumin 4.1 g/dL,",
      "100% on an immunomodulator.",
      "The external-validation cohort in Kim 2025 itself (57 Korean paediatric",
      "IBD patients with 158 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model was the most precise of the seven paediatric models",
      "(MPE 32.4%, MAE 1.7 mg/L, rRMSE 96.3%; Kim 2025 Table 2) but, like all of",
      "them, did not reach clinical acceptability.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. The source equations are stated per kilogram with a
    # 65 kg reference; the typical values below are the per-kg values multiplied
    # by 65 kg and converted from mL to L:
    #   lcl = log(5.42 * 65 / 1000) = log(0.3523 L/day)
    #   lvc = log(52.4 * 65 / 1000) = log(3.406 L)
    #   lvp = log(19.6 * 65 / 1000) = log(1.274 L)
    #   lq  = log(2.26 * 65 / 1000) = log(0.1469 L/day)
    lcl <- log(0.3523); label("Typical clearance for the reference subject (CL, L/day)")                    # Data S1 model #6 $THETA1 = 5.42 mL/kg/day * 65 kg / 1000; Table S2 "5.42"
    lvc <- log(3.406);  label("Typical central volume for the reference subject (Vc, L)")                   # Data S1 model #6 $THETA2 = 52.4 mL/kg * 65 kg / 1000; Table S2 "52.4"
    lvp <- log(1.274);  label("Typical peripheral volume for the reference subject (Vp, L)")                # Data S1 model #6 $THETA3 = 19.6 mL/kg * 65 kg / 1000; Table S2 "19.6"
    lq  <- log(0.1469); label("Typical inter-compartmental clearance for the reference subject (Q, L/day)") # Data S1 model #6 $THETA4 = 2.26 mL/kg/day * 65 kg / 1000; Table S2 "2.26"

    # Covariate effects. The three e_wt_* values are the exponents on the PER-KG
    # parameters; model() adds the +1 that converts them to total-parameter
    # exponents (0.687 on CL, 0.767 on Vc, 0.412 on Vp).
    #
    # DISCREPANCY: modellib('Frymoyer_2017_infliximab') transcribes the same
    # underlying Fasanmade model with these three exponents POSITIVE (+0.313,
    # +0.233, +0.588), giving total-parameter exponents of 1.313 / 1.233 / 1.588.
    # Kim 2025's control stream is unambiguous that they are negative, and the
    # negative signs are the physiologically expected ones: they give a total-CL
    # weight exponent of 0.687, close to the canonical allometric 0.75, whereas
    # the positive signs imply clearance rising faster than proportionally with
    # body weight. The two encodings agree exactly at the 65 kg reference and
    # diverge away from it (about 62% apart in CL at 30 kg). This file follows
    # Kim 2025; the divergence is flagged in the validation vignette.
    e_wt_cl  <- -0.313; label("Power exponent of body weight on per-kg CL ((WT/65)^(1 + e_wt_cl) on total CL)")  # Data S1 model #6 $THETA8 = -0.313
    e_wt_vc  <- -0.233; label("Power exponent of body weight on per-kg Vc ((WT/65)^(1 + e_wt_vc) on total Vc)")  # Data S1 model #6 $THETA9 = -0.233
    e_wt_vp  <- -0.588; label("Power exponent of body weight on per-kg Vp ((WT/65)^(1 + e_wt_vp) on total Vp)")  # Data S1 model #6 $THETA10 = -0.588
    e_alb_cl <- -0.855; label("Power exponent of serum albumin on CL ((ALB/4.1)^e_alb_cl)")                      # Data S1 model #6 $THETA5 = -0.855
    e_ada_cl <-  0.291; label("Fractional increase in CL when ADA-positive (1 + e_ada_cl)")                      # Data S1 model #6 $THETA6 = 0.291
    e_imm_cl <- -0.137; label("Fractional change in CL on a concomitant immunomodulator (1 + e_imm_cl)")         # Data S1 model #6 $THETA7 = -0.137

    # Inter-individual variability; OMEGA = (CV)^2 per Kim 2025's control-stream
    # annotation, with each %CV given in the stream comment and Table S2.
    etalcl ~ 0.0942   # Data S1 model #6 $OMEGA IIV_CL = 0.0942 (30.7% CV); Table S2 '(30.7%)'
    etalvc ~ 0.0159   # Data S1 model #6 $OMEGA IIV_V1 = 0.0159 (12.6% CV); Table S2 '(12.6%)'
    etalvp ~ 0.306    # Data S1 model #6 $OMEGA IIV_V2 = 0.306 (55.3% CV); Table S2 '(55.3%)'

    # Inter-occasion variability on CL. The source declares two occasion blocks
    # (OMEGA BLOCK(1) 0.0335 each) selected by FLAG1 / FLAG2, but hard-codes
    # FLAG1 = 1 and FLAG2 = 0 with the comment "In this cohort, all data were
    # from the maintenance phase", so exactly one occasion is ever active. Only
    # that occasion's random effect is carried here; the second, structurally
    # inert occasion is not declared.
    etaiov_cl_1 ~ 0.0335   # Data S1 model #6 $OMEGA BLOCK(1) IOV_CL = 0.0335 (18.3% CV), occasion 1

    # Residual error -- combined proportional + additive; SDs from the stream's
    # own comments (the $SIGMA entries are the corresponding variances).
    propSd <- 0.292; label("Proportional residual error (fraction)")  # Data S1 model #6 $SIGMA 0.0853, comment "proportional error (29.2 %CV)"
    addSd  <- 0.371; label("Additive residual error (mg/L)")         # Data S1 model #6 $SIGMA 0.138, comment "additive error (0.371 mg/L)"
  })
  model({
    # Canonical ALB is in SI g/L; this model was calibrated on g/dL.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL

    # Individual PK parameters. The published per-kg forms are
    #   CL (mL/kg/day) = 5.42 * (WT/65)^-0.313 * (alb_gdL/4.1)^-0.855
    #                         * (1 + 0.291*ADA) * (1 - 0.137*IMM)
    #   Vc (mL/kg)     = 52.4 * (WT/65)^-0.233
    #   Vp (mL/kg)     = 19.6 * (WT/65)^-0.588
    #   Q  (mL/kg/day) = 2.26  (constant per-kg)
    # Multiplying by WT to obtain totals carries an implicit exponent of +1 on
    # each parameter, which is what the (1 + e_wt_*) terms below express.
    # Inter-occasion variability on CL. The source applies
    # EXP(ETA(1) + FLAG1*ETA(4) + FLAG2*ETA(5)) with FLAG1 = 1 and FLAG2 = 0, so
    # only occasion 1 contributes. Routing it through a named intermediate keeps
    # the mu-referenced exponent to a single eta, which rxode2 requires.
    iov_cl <- etaiov_cl_1

    cl <- exp(lcl + etalcl + iov_cl) *
      (WT / 65)^(1 + e_wt_cl) *
      (alb_gdL / 4.1)^e_alb_cl *
      (1 + e_ada_cl * ADA_POS) *
      (1 + e_imm_cl * CONMED_IMMUNOMOD)
    vc <- exp(lvc + etalvc) * (WT / 65)^(1 + e_wt_vc)
    vp <- exp(lvp + etalvp) * (WT / 65)^(1 + e_wt_vp)
    q  <- exp(lq)           * (WT / 65)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model; infliximab is given as an IV infusion into central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L, numerically identical to ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
