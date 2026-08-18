Kim_2025_infliximab_fasanmade_pediatric <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in children with Crohn's disease, with inter-occasion variability on clearance (Fasanmade 2011 paediatric REACH model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
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
    "(reference 22 of Kim 2025).",
    "This is the paediatric-only parameterization fitted to the REACH trial",
    "alone; the pooled paediatric + adult parameterization from the same",
    "publication is modellib('Kim_2025_infliximab_fasanmade_combined').",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #7) with MAXEVAL=0,",
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
      notes              = "The source equations are written on per-kilogram CL / V values with reference 42 kg, so the total-parameter exponent is 1 plus the published per-kg exponent. Kim 2025's control stream gives per-kg exponents of -0.34 (CL), -0.171 (Vc) and -0.414 (Vp), i.e. total-parameter exponents of 0.66, 0.829 and 0.586; Q is a constant per-kg value, giving a total-Q exponent of 1.0. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement) and identified this as a probable contributor to the poor predictive performance of the paediatric models.",
      source_name        = "WGT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, normalized to a 3.8 g/dL reference: CL *= (ALB/3.8)^(-1.22). The source control stream states the covariate is in g/dL; the canonical ALB column is in SI g/L, so model() converts inline (g/L * 0.1 = g/dL). Time-varying in the Kim 2025 validation dataset.",
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112L,
    n_studies      = 1L,
    age_range      = "Paediatric; median age 13 years (Kim 2025 Table S1).",
    weight_range   = "Not reported; median body weight 42 kg (Kim 2025 Table S1), which is the reference weight used by the model.",
    sex_female_pct = 41,
    race_ethnicity = "Not specified.",
    disease_state  = "Crohn's disease (n = 112).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "Multi-regional paediatric Crohn's disease trial (REACH).",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Fasanmade_pediatric model: Crohn's disease (n = 112),",
      "paediatric patients, induction and maintenance phases, peak /",
      "intermediate / trough sampling, 3% ATI positive, median age 13 years,",
      "41% female, median body weight 42 kg, median albumin 3.8 g/dL, 100% on an",
      "immunomodulator.",
      "Unlike the pooled parameterization, this paediatric-only model retains no",
      "ADA or immunomodulator effect on clearance -- only albumin and body",
      "weight.",
      "The external-validation cohort in Kim 2025 itself (57 Korean paediatric",
      "IBD patients with 158 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model gave MPE 50.6%, MAE 1.7 mg/L and rRMSE 196.0%",
      "(Kim 2025 Table 2) and did not reach clinical acceptability.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. The source equations are stated per kilogram with a
    # 42 kg reference; the typical values below are the per-kg values multiplied
    # by 42 kg and converted from mL to L:
    #   lcl = log(5.43 * 42 / 1000) = log(0.22806 L/day)
    #   lvc = log(54.2 * 42 / 1000) = log(2.2764 L)
    #   lvp = log(29.2 * 42 / 1000) = log(1.2264 L)
    #   lq  = log(3.52 * 42 / 1000) = log(0.14784 L/day)
    lcl <- log(0.22806); label("Typical clearance for the reference subject (CL, L/day)")                    # Data S1 model #7 $THETA1 = 5.43 mL/kg/day * 42 kg / 1000; Table S2 "5.43"
    lvc <- log(2.2764);  label("Typical central volume for the reference subject (Vc, L)")                   # Data S1 model #7 $THETA2 = 54.2 mL/kg * 42 kg / 1000; Table S2 "54.2"
    lvp <- log(1.2264);  label("Typical peripheral volume for the reference subject (Vp, L)")                # Data S1 model #7 $THETA3 = 29.2 mL/kg * 42 kg / 1000; Table S2 "29.2"
    lq  <- log(0.14784); label("Typical inter-compartmental clearance for the reference subject (Q, L/day)") # Data S1 model #7 $THETA4 = 3.52 mL/kg/day * 42 kg / 1000; Table S2 "3.52"

    # Covariate effects. The three e_wt_* values are the exponents on the PER-KG
    # parameters; model() adds the +1 that converts them to total-parameter
    # exponents (0.66 on CL, 0.829 on Vc, 0.586 on Vp).
    e_wt_cl  <- -0.34;  label("Power exponent of body weight on per-kg CL ((WT/42)^(1 + e_wt_cl) on total CL)")  # Data S1 model #7 $THETA6 = -0.34
    e_wt_vc  <- -0.171; label("Power exponent of body weight on per-kg Vc ((WT/42)^(1 + e_wt_vc) on total Vc)")  # Data S1 model #7 $THETA7 = -0.171
    e_wt_vp  <- -0.414; label("Power exponent of body weight on per-kg Vp ((WT/42)^(1 + e_wt_vp) on total Vp)")  # Data S1 model #7 $THETA8 = -0.414
    e_alb_cl <- -1.22;  label("Power exponent of serum albumin on CL ((ALB/3.8)^e_alb_cl)")                      # Data S1 model #7 $THETA5 = -1.22

    # Inter-individual variability; OMEGA = (CV)^2 per Kim 2025's control-stream
    # annotation, with each %CV given in the stream comment and Table S2.
    etalcl ~ 0.064   # Data S1 model #7 $OMEGA IIV_CL = 0.064 (25.2% CV); Table S2 "(25.2%)"
    etalvc ~ 0.027   # Data S1 model #7 $OMEGA IIV_V1 = 0.027 (16.3% CV); Table S2 "(16.3%)"
    etalvp ~ 0.122   # Data S1 model #7 $OMEGA IIV_V2 = 0.122 (34.9% CV); Table S2 "(34.9%)"

    # Inter-occasion variability on CL. As in the pooled parameterization, the
    # source declares two occasion blocks (OMEGA BLOCK(1) 0.048 each) selected by
    # FLAG1 / FLAG2 but hard-codes FLAG1 = 1 and FLAG2 = 0 with the comment "In
    # this cohort, all data were from the maintenance phase", so exactly one
    # occasion is ever active. Only that occasion's random effect is carried
    # here; the second, structurally inert occasion is not declared.
    etaiov_cl_1 ~ 0.048   # Data S1 model #7 $OMEGA BLOCK(1) IOV_CL = 0.048 (21.9% CV), occasion 1

    # Residual error -- combined proportional + additive; SDs from the stream's
    # own comments (the $SIGMA entries are the corresponding variances).
    propSd <- 0.275; label("Proportional residual error (fraction)")  # Data S1 model #7 $SIGMA 0.076, comment "proportional error (27.5 %CV)"
    addSd  <- 0.244; label("Additive residual error (mg/L)")         # Data S1 model #7 $SIGMA 0.060, comment "additive error (0.244 mg/L)"
  })
  model({
    # Canonical ALB is in SI g/L; this model was calibrated on g/dL.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL

    # Individual PK parameters. The published per-kg forms are
    #   CL (mL/kg/day) = 5.43 * (WT/42)^-0.34 * (alb_gdL/3.8)^-1.22
    #   Vc (mL/kg)     = 54.2 * (WT/42)^-0.171
    #   Vp (mL/kg)     = 29.2 * (WT/42)^-0.414
    #   Q  (mL/kg/day) = 3.52  (constant per-kg)
    # Multiplying by WT to obtain totals carries an implicit exponent of +1 on
    # each parameter, which is what the (1 + e_wt_*) terms below express.
    # Inter-occasion variability on CL. The source applies
    # EXP(ETA(1) + FLAG1*ETA(4) + FLAG2*ETA(5)) with FLAG1 = 1 and FLAG2 = 0, so
    # only occasion 1 contributes. Routing it through a named intermediate keeps
    # the mu-referenced exponent to a single eta, which rxode2 requires.
    iov_cl <- etaiov_cl_1

    cl <- exp(lcl + etalcl + iov_cl) *
      (WT / 42)^(1 + e_wt_cl) *
      (alb_gdL / 3.8)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT / 42)^(1 + e_wt_vc)
    vp <- exp(lvp + etalvp) * (WT / 42)^(1 + e_wt_vp)
    q  <- exp(lq)           * (WT / 42)

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
