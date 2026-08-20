Kim_2025_infliximab_passot_pediatric <- function() {
  description <- "One-compartment population PK model of intravenous infliximab in children with inflammatory bowel disease (Passot 2016 paediatric model with an age effect on central volume, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Passot C, Mulleman D, Bejan-Angoulvant T, et al. The underlying",
    "inflammatory chronic disease influences infliximab pharmacokinetics.",
    "MAbs. 2016;8(7):1407-1416. doi:10.1080/19420862.2016.1216741",
    "(reference 18 of Kim 2025).",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #8) with MAXEVAL=0,",
    "i.e. every THETA / OMEGA / SIGMA is held at the originally published value.",
    "Values here are transcribed from that control stream and cross-checked",
    "against Kim 2025 Table S2.",
    "Kim 2025 Table S1 footnote states that this paediatric variant is identical",
    "to the adult variant (modellib('Kim_2025_infliximab_passot_adult')) except",
    "for the addition of age as a covariate on the central volume.",
    sep = " "
  )
  vignette <- "Kim_2025_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on both CL and Vc, normalized to a 67 kg reference: CL *= (WT/67)^0.603 and Vc *= (WT/67)^0.277. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement) and identified this as a probable contributor to the poor predictive performance of the paediatric models, because weight changes materially over a long paediatric follow-up.",
      source_name        = "WGT"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Exponential effect on CL and Vc with female as the reference level: the multiplier is 1 for females and exp(0.181) = 1.198 on CL and exp(0.209) = 1.232 on Vc for males. The source control stream codes SEX = 1 for female and SEX = 0 for male, matching the canonical SEXF polarity directly (no value inversion needed). Time-invariant.",
      source_name        = "SEX"
    ),
    IBD_CD = list(
      description        = "Inflammatory bowel disease subtype (1 = Crohn's disease, 0 = ulcerative colitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "no reference level within IBD: Crohn's disease and ulcerative colitis each carry their own multiplier",
      notes              = "The source paper (Passot 2016) fitted infliximab PK across several chronic inflammatory diseases, so neither IBD subtype is the model's global reference category; both Crohn's disease and ulcerative colitis receive an explicit exponential multiplier on CL and on Vc, applied on top of the typical values THETA1 / THETA2. On CL: exp(0.384) = 1.468 for Crohn's disease and exp(0.472) = 1.603 for ulcerative colitis. On Vc: exp(0.399) = 1.490 for Crohn's disease and exp(0.417) = 1.517 for ulcerative colitis. The source control stream codes DX = 0 for Crohn's disease and DX = 1 for ulcerative colitis, which is the inverse of the canonical IBD_CD polarity; derive IBD_CD = 1 - DX. Time-invariant.",
      source_name        = "DX"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = "age >= 15 years (multiplier 1)",
      notes              = "Entered as a dichotomy at 15 years rather than as a continuous effect, and only on Vc: the multiplier is 1 for age >= 15 years and exp(-0.396) = 0.673 for age < 15 years. This is the only covariate that distinguishes the paediatric variant of the Passot model from the adult variant. Time-invariant in the Kim 2025 validation dataset.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 79L,
    n_studies      = 1L,
    age_range      = "Pooled paediatric and adult cohort; age not specified in Kim 2025 Table S1 for this model.",
    weight_range   = "Not specified in Kim 2025 Table S1; the model's reference weight is 67 kg.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not specified; developed in a French cohort.",
    disease_state  = "Inflammatory bowel disease: Crohn's disease (n = 63) and ulcerative colitis (n = 16).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "France.",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Passot model: Crohn's disease (n = 63) and ulcerative",
      "colitis (n = 16), pooled paediatric and adult patients, induction and",
      "maintenance phases, trough sampling only, 0% ATI positive. Age, sex,",
      "body weight, albumin, ESR and immunomodulator use were not specified.",
      "The external-validation cohort in Kim 2025 itself (57 Korean paediatric",
      "IBD patients with 158 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters. In that cohort",
      "this model gave MPE 45.2%, MAE 1.8 mg/L and rRMSE 112.4%",
      "(Kim 2025 Table 2). Kim 2025 concluded that none of the seven paediatric",
      "models reached clinically acceptable predictive performance in Korean",
      "children, including after supplying one or two previous trough",
      "concentrations, and recommended that paediatric-specific models be",
      "developed for this population.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- typical values before any covariate multiplier.
    # Neither IBD subtype is a unit multiplier (see e_ibd_* below), so these are
    # not the values for a "reference patient" on their own.
    lcl <- log(0.23); label("Typical clearance before covariate effects (CL, L/day)")   # Data S1 model #8 $THETA1 = 0.23 L/day; Table S2 "0.23"
    lvc <- log(5.2);  label("Typical central volume before covariate effects (Vc, L)")  # Data S1 model #8 $THETA2 = 5.2 L; Table S2 "5.2"

    # Covariate effects on CL -- identical to the adult variant.
    e_wt_cl     <- 0.603; label("Power exponent of body weight on CL ((WT/67)^e_wt_cl)")                     # Data S1 model #8 $THETA3 = 0.603
    e_sex_cl    <- 0.181; label("Log-scale male effect on CL (multiplier exp(e_sex_cl) when SEXF = 0)")      # Data S1 model #8 $THETA4 = 0.181
    e_ibd_cd_cl <- 0.384; label("Log-scale Crohn's disease effect on CL (multiplier exp(e_ibd_cd_cl))")      # Data S1 model #8 $THETA5 = 0.384
    e_ibd_uc_cl <- 0.472; label("Log-scale ulcerative colitis effect on CL (multiplier exp(e_ibd_uc_cl))")   # Data S1 model #8 $THETA6 = 0.472

    # Covariate effects on Vc -- identical to the adult variant, plus the age
    # effect that is unique to this paediatric variant.
    e_wt_vc     <-  0.277; label("Power exponent of body weight on Vc ((WT/67)^e_wt_vc)")                    # Data S1 model #8 $THETA7 = 0.277
    e_sex_vc    <-  0.209; label("Log-scale male effect on Vc (multiplier exp(e_sex_vc) when SEXF = 0)")     # Data S1 model #8 $THETA8 = 0.209
    e_ibd_cd_vc <-  0.399; label("Log-scale Crohn's disease effect on Vc (multiplier exp(e_ibd_cd_vc))")     # Data S1 model #8 $THETA9 = 0.399
    e_ibd_uc_vc <-  0.417; label("Log-scale ulcerative colitis effect on Vc (multiplier exp(e_ibd_uc_vc))")  # Data S1 model #8 $THETA10 = 0.417
    e_age_vc    <- -0.396; label("Log-scale effect on Vc for age < 15 years (multiplier exp(e_age_vc))")     # Data S1 model #8 $THETA11 = -0.396

    # Inter-individual variability. The control stream gives the SD in each
    # comment: 0.092 = 0.304^2 (CL); 0.050 = 0.224^2 (Vc).
    etalcl ~ 0.092   # Data S1 model #8 $OMEGA ETA1 = 0.092 (SD 0.304); Table S2 '(30.3%)'
    etalvc ~ 0.050   # Data S1 model #8 $OMEGA ETA2 = 0.050 (SD 0.224); Table S2 '(22.4%)'

    # Residual error -- combined proportional + additive; SDs from the stream's
    # own comments (the $SIGMA entries are the corresponding variances).
    propSd <- 0.223; label("Proportional residual error (fraction)")  # Data S1 model #8 $SIGMA 0.05, comment "proportional error (0.223)"
    addSd  <- 0.72;  label("Additive residual error (mg/L)")         # Data S1 model #8 $SIGMA 0.518, comment "additive error (0.72)"
  })
  model({
    # Individual PK parameters. Exactly one of the two e_ibd_* terms is active
    # per subject; the age term is a step at 15 years applied to Vc only:
    #   CL (L/day) = 0.23 * (WT/67)^0.603 * exp(0.181)^(male)
    #                     * exp(0.384)^(CD) * exp(0.472)^(UC)
    #   Vc (L)     = 5.2  * (WT/67)^0.277 * exp(0.209)^(male)
    #                     * exp(0.399)^(CD) * exp(0.417)^(UC)
    #                     * exp(-0.396)^(AGE < 15)
    cl <- exp(lcl +
                e_sex_cl * (1 - SEXF) +
                e_ibd_cd_cl * IBD_CD +
                e_ibd_uc_cl * (1 - IBD_CD) +
                etalcl) * (WT / 67)^e_wt_cl
    vc <- exp(lvc +
                e_sex_vc * (1 - SEXF) +
                e_ibd_cd_vc * IBD_CD +
                e_ibd_uc_vc * (1 - IBD_CD) +
                e_age_vc * (AGE < 15) +
                etalvc) * (WT / 67)^e_wt_vc

    kel <- cl / vc

    # One-compartment model; infliximab is given as an IV infusion into central.
    # The source control stream declares a single compartment but also carries a
    # stray "DADT(2) = A(1)/V1" line, which is an unused cumulative-AUC
    # accumulator (compartment 2 is never declared in $MODEL and never enters
    # $ERROR); it is intentionally not reproduced here.
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L, numerically identical to ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
