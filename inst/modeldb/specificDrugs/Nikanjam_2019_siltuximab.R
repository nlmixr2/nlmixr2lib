Nikanjam_2019_siltuximab <- function() {
  description <- "Two-compartment population PK model for siltuximab (anti-IL-6) in adults pooled across healthy volunteers and oncology cohorts including Castleman's disease, smoldering multiple myeloma, and other tumor types (Nikanjam 2019)"
  reference <- "Nikanjam M, Cho HJ, Capparelli EV. Population pharmacokinetics of siltuximab: impact of disease state. Cancer Chemother Pharmacol. 2019;84(5):993-1001. doi:10.1007/s00280-019-03939-7"
  vignette <- "Nikanjam_2019_siltuximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on V1 and V2 (combined Vss); normalized as WT/73 per Nikanjam 2019 (reference: 73 kg, overall median from Table 1).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on CL and on V1/V2; normalized as ALB/4.1 per Nikanjam 2019 (reference: 4.1 g/dL, overall median from Table 1).",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALT/19 per Nikanjam 2019 (reference: 19 U/L, overall median from Table 1).",
      source_name        = "ALT"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on V1 and V2; normalized as CREAT/0.9 per Nikanjam 2019 (reference: 0.9 mg/dL, overall median from Table 1). Source paper uses the column label 'SCR' (serum creatinine); renamed to the canonical CREAT per covariate-columns.md.",
      source_name        = "SCR"
    ),
    DIS_HV = list(
      description        = "Healthy-volunteer cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient subject from any of the pooled disease cohorts: Castleman's disease, smoldering multiple myeloma, multiple myeloma, MGUS, renal cell carcinoma, ovarian cancer, KRAS-mutated tumors, or other solid tumors)",
      notes              = "Multiplicative effects in Nikanjam 2019: 0.77 on CL (a 23% reduction in healthy volunteers vs. the pooled patient reference) and 0.83 on V1/V2 (a 17% reduction).",
      source_name        = "HV"
    ),
    DIS_CASTLEMAN = list(
      description        = "Castleman's disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Castleman subject; pooled across healthy volunteers and other indications)",
      notes              = "Multiplicative effect in Nikanjam 2019: 1.24 on CL (a 24% increase in Castleman's disease vs. the pooled non-Castleman reference). No effect on V1/V2.",
      source_name        = "CD"
    ),
    DIS_SMM = list(
      description        = "Smoldering multiple myeloma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-SMM subject; pooled across healthy volunteers and other indications)",
      notes              = "Multiplicative effect in Nikanjam 2019: 0.77 on V1/V2 (a 23% reduction in SMM vs. the pooled non-SMM reference). No effect on CL.",
      source_name        = "SMM"
    )
  )

  population <- list(
    n_subjects     = 460L,
    n_studies      = 7L,
    age_range      = "18-84 years (overall, Table 1)",
    age_median     = "51 years",
    weight_range   = "40-170 kg",
    weight_median  = "73 kg",
    sex_female_pct = 47.8,
    race_ethnicity = "White, Asian, African-American, and Other (per-study composition reported in source; race was tested as a covariate but not retained in the final model).",
    disease_state  = "Pooled cohort: healthy volunteers (T08, n=141) plus oncology and plasma-cell-disorder patients across renal cell carcinoma (T01), non-Hodgkin's lymphoma / multiple myeloma / Castleman's disease (T03), multicentric Castleman's disease (MCD2001), MGUS / smoldering and active multiple myeloma (SMM1001), high-risk smoldering multiple myeloma (SMM2001), and advanced cancer including ovarian and KRAS-mutated tumors (STM2001).",
    dose_range     = "0.15-15.6 mg/kg IV infusion (1-2 h) every 1-4 weeks; pooled across single-dose (T08) and multiple-dose (other) regimens.",
    regions        = "Multi-center; geographic detail not specified in the source.",
    notes          = "Baseline demographics and laboratory values pooled from Nikanjam 2019 Table 1. Reference covariate values used in the final-model parameterization: WT = 73 kg, ALB = 4.1 g/dL, ALT = 19 U/L, CREAT = 0.9 mg/dL, with HV = CD = SMM = 0 (i.e., the reference subject is a patient with a non-Castleman, non-SMM tumor type). 7,761 siltuximab concentrations contributed to the analysis."
  )

  ini({
    # Structural parameters (typical values for the reference subject:
    # 73-kg patient with ALB 4.1 g/dL, ALT 19 U/L, CREAT 0.9 mg/dL, and a
    # non-Castleman / non-SMM / non-HV oncology indication). Source: Nikanjam
    # 2019 Table 2 (final population PK model parameters).
    lcl <- log(0.214); label("Clearance (CL, L/day)")                       # Nikanjam 2019 Table 2: theta2
    lvc <- log(3.66);  label("Central volume of distribution (V1, L)")      # Nikanjam 2019 Table 2: theta1
    lvp <- log(3.13);  label("Peripheral volume of distribution (V2, L)")   # Nikanjam 2019 Table 2: theta3
    lq  <- log(0.624); label("Inter-compartmental clearance (Q, L/day)")    # Nikanjam 2019 Table 2: theta4

    # Covariate effects (Nikanjam 2019 Table 2). The continuous-covariate
    # exponents act as power-form effects normalized to the population
    # medians; the disease-state coefficients act as multiplicative scalars
    # raised to the {0,1} indicator (e.g., HV=1 multiplies CL by e_hv_cl).
    e_alb_cl   <- -0.84;   label("Power exponent of albumin on CL (unitless)")              # Nikanjam 2019 Table 2: theta5
    e_wt_v     <-  0.65;   label("Power exponent of weight on V1 and V2 (unitless)")        # Nikanjam 2019 Table 2: theta6
    e_alb_v    <- -0.35;   label("Power exponent of albumin on V1 and V2 (unitless)")       # Nikanjam 2019 Table 2: theta7
    e_hv_v     <-  0.83;   label("Multiplicative effect of healthy volunteer on V1/V2")      # Nikanjam 2019 Table 2: theta8
    e_smm_v    <-  0.77;   label("Multiplicative effect of smoldering MM on V1/V2")          # Nikanjam 2019 Table 2: theta9
    e_hv_cl    <-  0.77;   label("Multiplicative effect of healthy volunteer on CL")         # Nikanjam 2019 Table 2: theta10
    e_alt_cl   <- -0.096;  label("Power exponent of ALT on CL (unitless)")                   # Nikanjam 2019 Table 2: theta11
    e_creat_v  <-  0.16;   label("Power exponent of serum creatinine on V1 and V2")          # Nikanjam 2019 Table 2: theta12
    e_cd_cl    <-  1.24;   label("Multiplicative effect of Castleman's disease on CL")        # Nikanjam 2019 Table 2: theta13

    # Inter-individual variability: a single eta on CL (BSV_CL = 20.0% CV)
    # and a single eta shared between V1 and V2 (BSV_Vss = 42.0% CV).
    # omega^2 = log(CV^2 + 1) for log-normal IIV.
    etalcl ~ 0.0392   # log(1 + 0.20^2);   Nikanjam 2019 Table 2: BSV CL = 20.0% CV
    etalvc ~ 0.1625   # log(1 + 0.42^2);   Nikanjam 2019 Table 2: BSV Vss = 42.0% CV (shared on V1 and V2)

    # Residual error (combined additive + proportional). Source: Nikanjam 2019 Table 2.
    propSd <- 0.224;  label("Proportional residual error (fraction)")        # Nikanjam 2019 Table 2: 22.4%
    addSd  <- 0.0217; label("Additive residual error (ug/mL)")                # Nikanjam 2019 Table 2: 0.0217
  })
  model({
    # Individual PK parameters. Reference subject: 73-kg patient with
    # ALB 4.1 g/dL, ALT 19 U/L, CREAT 0.9 mg/dL, and a non-HV, non-Castleman,
    # non-SMM oncology indication. Continuous covariates enter as power-form
    # effects normalized to the population medians; disease-state indicators
    # enter as multiplicative scalars raised to the binary indicator.
    cl <- exp(lcl + etalcl) *
      (ALB / 4.1)^e_alb_cl *
      (ALT / 19)^e_alt_cl *
      e_hv_cl^DIS_HV *
      e_cd_cl^DIS_CASTLEMAN

    vc <- exp(lvc + etalvc) *
      (WT  / 73)^e_wt_v *
      (ALB / 4.1)^e_alb_v *
      (CREAT / 0.9)^e_creat_v *
      e_hv_v^DIS_HV *
      e_smm_v^DIS_SMM

    vp <- exp(lvp + etalvc) *
      (WT  / 73)^e_wt_v *
      (ALB / 4.1)^e_alb_v *
      (CREAT / 0.9)^e_creat_v *
      e_hv_v^DIS_HV *
      e_smm_v^DIS_SMM

    q <- exp(lq)

    # Two-compartment IV model. All siltuximab dosing in the source is IV
    # infusion; no extravascular depot is needed.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
