Chudasama_2026_enfortumab <- function() {
  description <- paste(
    "Population pharmacokinetic / tumor growth inhibition (PK-TGI) model for",
    "enfortumab vedotin monotherapy in locally advanced or metastatic",
    "urothelial cancer (Chudasama 2026; 603 patients from EV-101, EV-102,",
    "EV-201 and EV-301, 0.5-1.25 mg/kg on days 1, 8 and 15 of a 28-day",
    "cycle). Claret-type TGI: first-order (exponential) tumor growth minus a",
    "tumor shrinkage rate driven by the conjugated-ADC plasma concentration",
    "through an Emax term, damped over time by a resistance rate and held",
    "above zero by a residual FLOOR effect: d/dt(tumor_size) = p *",
    "tumor_size - kkill * (kill_floor + exp(-res * t)) * Cc / (ic50 + Cc) *",
    "tumor_size, with tumor_size(0) = observed baseline sum of diameters.",
    "Covariates on the growth rate are baseline tumor size, albumin,",
    "hemoglobin and >= 3 prior lines of therapy. The ADC concentration comes",
    "from the linear three-compartment enfortumab vedotin population PK",
    "model of Zuo 2024 (final model, all covariates), which Chudasama 2026",
    "used through its individual post hoc estimates; its parameters are held",
    "fixed here."
  )
  reference <- paste(
    "Chudasama VL, Patilea-Vrana G, Mu H, Hamed S, Poondru S, Garg A, Chia YL.",
    "Evaluating the Impact of Enfortumab Vedotin Dose Modifications on",
    "Clinical Efficacy Using Tumor Growth Inhibition Modeling in Patients",
    "With Advanced Urothelial Cancer. J Clin Pharmacol. 2026.",
    "doi:10.1002/jcph.70241. PMID 42466532; PMCID PMC13377668.",
    "PK layer: Zuo P, Bonate P, Garg A, Matsangou M, Tang M. Population",
    "Pharmacokinetic Modeling and Exposure-Response Analysis for the",
    "Antibody-Drug Conjugate Enfortumab Vedotin in Locally Advanced or",
    "Metastatic Urothelial Carcinoma. Clin Pharmacol Ther.",
    "2024;116(5):1278-1288. doi:10.1002/cpt.3383."
  )
  vignette <- "Chudasama_2026_enfortumab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight. Allometric power covariate on all three clearances and all three volumes of the ADC PK layer (Zuo 2024), reference 75 kg.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Zuo 2024 Table 1 shares one exponent across CLA, Q2 and Q3 and another across V1, V2 and V3. The enfortumab vedotin dose is 1.25 mg/kg capped at 125 mg (body weight > 100 kg); the cap is applied in the dose, not in the model. TGI analysis median 74.2 kg (range 36.9-158; Chudasama 2026 Table S10).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at baseline. Power covariate on ADC clearance (Zuo 2024), reference 68 years.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "TGI analysis median 68 years (range 24-90; Chudasama 2026 Table S10).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male. Proportional-change covariate on ADC clearance, central volume and second peripheral volume (Zuo 2024).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Zuo 2024 Table 1 'Sex (female)' effects; the reference population is male. 25.4% of the TGI analysis set was female (Chudasama 2026 Table S10).",
      source_name = "SEX"
    ),
    ALB = list(
      description = "Baseline serum albumin. Power covariate on ADC clearance (Zuo 2024, reference 39 g/L) and on the tumor growth rate p (Chudasama 2026, reference 39 g/L).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Chudasama 2026 imputed missing values (n = 4) with the population median (control-stream column ALBMissAsMed). Zuo 2024's reference population text prints '39 g/mL', a unit typo for g/L. TGI analysis median 39.0 g/L (range 25.0-51.0; Table S10).",
      source_name = "ALBMissAsMed"
    ),
    HGB = list(
      description = "Baseline hemoglobin. Power covariate on the ADC central volume (Zuo 2024, reference 112 g/L) and on the tumor growth rate p (Chudasama 2026, reference 112 g/L).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Chudasama 2026 imputed missing values (n = 11) with the population median (control-stream column HGBMissAsMed). TGI analysis median 112 g/L (range 79-172; Table S10).",
      source_name = "HGBMissAsMed"
    ),
    TUMSZ = list(
      description = "Baseline sum of diameters (SOD) of RECIST 1.1 target lesions. Initial condition of the tumor_size state, power covariate on the tumor growth rate p (reference 64.1 mm) and on ADC clearance (Zuo 2024, reference 65 mm).",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Chudasama 2026 control stream: A_0(4) = SODB and (SODB2/64.1)^THETA on KG; the TGI dataset median was 64.1 mm (range 11-331; Table S10). Zuo 2024 standardises the same baseline SOD to its PK-population median of 65 mm.",
      source_name = "SODB"
    ),
    LINE_4L_PLUS = list(
      description = "1 = three or more prior lines of anticancer therapy in all settings (so enfortumab vedotin is fourth line or later), 0 = zero to two prior lines. Proportional-change covariate on the tumor growth rate p.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Chudasama 2026 control-stream column PTHRALL2N (1 = 0-2 prior lines, 3 = 3+ prior lines; the mrgsolve code sets PTHRALL2N3 = 1 only for PTHRALL2N == 3), Table S9 'Prior therapies (3+) on KG'. The Results text says '>3 prior lines', but the Methods ('3+ lines vs all others'), Table S9 and the code all define the category as three or more. Prior lines are counted, so LINE_4L_PLUS = as.integer(n_prior_lines >= 3). 29.9% of the TGI analysis set (180/603).",
      source_name = "PTHRALL2N"
    ),
    FORM_ENFORTUMAB_PROCESSA = list(
      description = "1 = enfortumab vedotin drug product made by manufacturing process A, 0 = the other manufacturing process. Proportional-change covariate on ADC clearance (Zuo 2024).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Zuo 2024 Table 1 'Manufacturer process (A) ~ CLA' = -0.110; the reference population in Zuo 2024 Figure 2 is the non-A process. Zuo 2024 does not name the other process; set 0 to simulate the reference material.",
      source_name = "Manufacturer process"
    ),
    ASSAY_LAB_INTERTEK = list(
      description = "1 = enfortumab vedotin concentrations quantified by the Intertek analytical laboratory, 0 = the other laboratory. Proportional-change covariate on the ADC second peripheral volume (Zuo 2024).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Zuo 2024 Table 1 'Analytic laboratory (Intertek) ~ V3' = -0.355. A bioanalytical-site effect, not a patient property: simulations of new patients set it to 0.",
      source_name = "Analytic laboratory"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "enfortumab vedotin (conjugated antibody-drug conjugate)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "enfortumab vedotin (conjugated antibody-drug conjugate)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "enfortumab vedotin (conjugated antibody-drug conjugate)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_size = list(
      analyte = "tumor (sum of RECIST 1.1 target-lesion diameters)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 603L,
    n_studies = 4L,
    age_range = "24-90 years (median 68)",
    weight_range = "36.9-158 kg (median 74.2)",
    sex_female_pct = 25.4,
    race_ethnicity = "Non-Asian 70.0%, Asian 23.4%, missing 6.6% (Chudasama 2026 Table S10)",
    disease_state = "Locally advanced or metastatic urothelial cancer with a measurable baseline tumor assessment and at least one non-zero post-baseline assessment; most patients previously treated with platinum chemotherapy and/or a PD-1/PD-L1 inhibitor.",
    dose_range = "Enfortumab vedotin 0.5, 0.75, 1.0 or 1.25 mg/kg (91.4% at 1.25 mg/kg) as a 30-minute intravenous infusion on days 1, 8 and 15 of a 28-day cycle, with protocol dose reductions to 1.0 and 0.75 mg/kg.",
    regions = "Global (EV-101 NCT02091999, EV-102 NCT03070990 Japan, EV-201 NCT03219333, EV-301 NCT03474107).",
    notes = "2167 tumor SOD observations. Baseline tumor size median 64.1 mm (11-331); albumin < 35 g/L in 15.4% and hemoglobin < 100 g/L in 22.4%; 0 / 1 / 2 / >= 3 prior lines of therapy in 0.5% / 15.8% / 53.9% / 29.9% (Results; Table S10). Seven patients with a post-baseline SOD of 0 mm were excluded. The ADC PK layer was developed by Zuo 2024 on 748 patients from five studies (median 75 kg, 68 years)."
  )

  ini({
    # ADC pharmacokinetics: Zuo 2024 Table 1, final model. Chudasama 2026 drove
    # the TGI with individual post hoc PK parameters from this model (control
    # stream $INPUT CLD V1 Q2D V2 Q3D V3), so every PK value is fixed here.
    # Zuo 2024 reports clearances in L/h; model() converts them to L/day.
    lcl <- fixed(log(0.110)); label("ADC clearance CLA (L/h)") # Zuo 2024 Table 1: CLA = 0.110 L/h (RSE 1.3%)
    lvc <- fixed(log(3.63)); label("ADC central volume V1 (L)") # Zuo 2024 Table 1: V1 = 3.63 L (RSE 1.1%)
    lq <- fixed(log(0.00446)); label("ADC intercompartmental clearance to the first peripheral compartment Q2 (L/h)") # Zuo 2024 Table 1: Q2 = 0.00446 L/h (RSE 7.0%)
    lvp <- fixed(log(5.79)); label("ADC first peripheral volume V2 (L)") # Zuo 2024 Table 1: V2 = 5.79 L (RSE 12%)
    lq2 <- fixed(log(0.0397)); label("ADC intercompartmental clearance to the second peripheral compartment Q3 (L/h)") # Zuo 2024 Table 1: Q3 = 0.0397 L/h (RSE 2.6%)
    lvp2 <- fixed(log(3.42)); label("ADC second peripheral volume V3 (L)") # Zuo 2024 Table 1: V3 = 3.42 L (RSE 2.9%)

    e_wt_cl <- fixed(0.572); label("Power exponent of body weight on CLA, Q2 and Q3 (unitless)") # Zuo 2024 Table 1: Weight ~ CLA, Q2, Q3 = 0.572 (RSE 9.0%)
    e_wt_vc <- fixed(0.607); label("Power exponent of body weight on V1, V2 and V3 (unitless)") # Zuo 2024 Table 1: Weight ~ V1, V2, V3 = 0.607 (RSE 8.0%)
    e_age_cl <- fixed(-0.325); label("Power exponent of age on CLA (unitless)") # Zuo 2024 Table 1: Age ~ CLA = -0.325 (RSE 17%)
    e_alb_cl <- fixed(-0.507); label("Power exponent of albumin on CLA (unitless)") # Zuo 2024 Table 1: Albumin ~ CLA = -0.507 (RSE 16%)
    e_form_processa_cl <- fixed(-0.110); label("Proportional change in CLA for manufacturing process A (unitless)") # Zuo 2024 Table 1: Manufacturer process (A) ~ CLA = -0.110 (RSE 14%)
    e_sexf_cl <- fixed(-0.138); label("Proportional change in CLA for female sex (unitless)") # Zuo 2024 Table 1: Sex (female) ~ CLA = -0.138 (RSE 14%)
    e_tumsz_cl <- fixed(0.0555); label("Power exponent of baseline tumor SOD on CLA (unitless)") # Zuo 2024 Table 1: SOD ~ CLA = 0.0555 (RSE 22%)
    e_hgb_vc <- fixed(-0.380); label("Power exponent of hemoglobin on V1 (unitless)") # Zuo 2024 Table 1: Hemoglobin ~ V1 = -0.380 (RSE 17%)
    e_sexf_vc <- fixed(-0.148); label("Proportional change in V1 for female sex (unitless)") # Zuo 2024 Table 1: Sex (female) ~ V1 = -0.148 (RSE 14%)
    e_assay_intertek_vp2 <- fixed(-0.355); label("Proportional change in V3 for the Intertek analytical laboratory (unitless)") # Zuo 2024 Table 1: Analytic laboratory (Intertek) ~ V3 = -0.355 (RSE 7.8%)
    e_sexf_vp2 <- fixed(-0.139); label("Proportional change in V3 for female sex (unitless)") # Zuo 2024 Table 1: Sex (female) ~ V3 = -0.139 (RSE 28%)

    # Zuo 2024 Table 1 reports IIV as CV%; omega^2 = log(1 + CV^2).
    etalcl ~ fixed(0.040383) # Zuo 2024 Table 1: omega CLA 20.3 %CV
    etalvc ~ fixed(0.044370) # Zuo 2024 Table 1: omega V1 21.3 %CV
    etalvp2 ~ fixed(0.130920) # Zuo 2024 Table 1: omega V3 37.4 %CV

    propSd <- fixed(0.204); label("Proportional residual error on ADC concentration (fraction)") # Zuo 2024 Table 1: proportional error 20.4 %CV (RSE 2.3%)
    addSd <- fixed(0.0459); label("Additive residual error on ADC concentration (ug/mL)") # Zuo 2024 Table 1: additive error 0.0459 ug/mL (RSE 9.1%)

    # Tumor growth inhibition: Chudasama 2026 Table S9, final model.
    lp <- log(0.00137); label("Apparent first-order tumor growth rate KG (1/day)") # Chudasama 2026 Table S9 final model: KG = 0.00137 1/day
    lkkill <- log(0.0327); label("Maximum ADC-driven tumor shrinkage rate KS (1/day)") # Chudasama 2026 Table S9 final model: KS = 0.0327 1/day
    lres <- log(0.0278); label("Resistance rate damping the shrinkage effect over time, lambda (1/day)") # Chudasama 2026 Table S9 final model: LAMBDA = 0.0278 1/day
    lic50 <- log(3.40); label("ADC concentration at half the maximum shrinkage rate IC50 (ug/mL)") # Chudasama 2026 Table S9 final model: IC50 = 3.40 ug/mL
    lkill_floor <- log(0.149); label("Residual shrinkage effect at infinite time, FLOOR, as a multiple of kkill (unitless)") # Chudasama 2026 Table S9 final model: FLOOR = 0.149

    e_tumsz_p <- 0.664; label("Power exponent of baseline tumor SOD on KG (unitless)") # Chudasama 2026 Table S9 final model: Baseline tumor SOD on KG = 0.664
    e_alb_p <- -0.593; label("Power exponent of baseline albumin on KG (unitless)") # Chudasama 2026 Table S9 final model: Baseline albumin on KG = -0.593
    e_hgb_p <- -0.421; label("Power exponent of baseline hemoglobin on KG (unitless)") # Chudasama 2026 Table S9 final model: Baseline hemoglobin on KG = -0.421
    e_line_4l_plus_p <- 0.354; label("Proportional change in KG for >= 3 prior lines of therapy (unitless)") # Chudasama 2026 Table S9 final model: Prior therapies (3+) on KG = 0.354

    # Table S9 lists the IIV as 'Var[...]', i.e. variances of exponential etas.
    etalp ~ 1.01 # Chudasama 2026 Table S9 final model: IIV on KG variance 1.01
    etalic50 ~ 3.51 # Chudasama 2026 Table S9 final model: IIV on IC50 variance 3.51

    addSd_tumor_size <- 7.85; label("Additive residual error on tumor SOD (mm)") # Chudasama 2026 Table S9 final model: ADD RE = 7.85 mm
  })

  model({
    # ADC PK (Zuo 2024): continuous covariates are power functions of the
    # covariate over its reference value, categorical covariates multiply by
    # (1 + theta * indicator). The factor 24 converts L/h to L/day.
    cl <- exp(lcl + etalcl) * 24 *
      (WT / 75)^e_wt_cl * (AGE / 68)^e_age_cl * (ALB / 39)^e_alb_cl *
      (TUMSZ / 65)^e_tumsz_cl *
      (1 + e_form_processa_cl * FORM_ENFORTUMAB_PROCESSA) *
      (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc * (HGB / 112)^e_hgb_vc *
      (1 + e_sexf_vc * SEXF)
    q <- exp(lq) * 24 * (WT / 75)^e_wt_cl
    vp <- exp(lvp) * (WT / 75)^e_wt_vc
    q2 <- exp(lq2) * 24 * (WT / 75)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 75)^e_wt_vc *
      (1 + e_assay_intertek_vp2 * ASSAY_LAB_INTERTEK) *
      (1 + e_sexf_vp2 * SEXF)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # TGI (Chudasama 2026 final model, Table S3 control stream and Table S4
    # mrgsolve code): KG carries the covariates, IC50 the IIV only.
    p <- exp(lp + etalp) *
      (TUMSZ / 64.1)^e_tumsz_p * (ALB / 39)^e_alb_p * (HGB / 112)^e_hgb_p *
      (1 + e_line_4l_plus_p * LINE_4L_PLUS)
    kkill <- exp(lkkill)
    res <- exp(lres)
    ic50 <- exp(lic50 + etalic50)
    kill_floor <- exp(lkill_floor)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ADC concentration in mg/L = ug/mL, the units of IC50.
    Cc <- central / vc

    # t is time since the first enfortumab vedotin dose (days); the resistance
    # term exp(-res * t) is not reset by later doses.
    drug_effect <- Cc / (ic50 + Cc)
    tumor_size(0) <- TUMSZ
    d/dt(tumor_size) <- p * tumor_size - kkill * (kill_floor + exp(-res * t)) * drug_effect * tumor_size

    Cc ~ add(addSd) + prop(propSd)
    tumor_size ~ add(addSd_tumor_size)
  })
}
