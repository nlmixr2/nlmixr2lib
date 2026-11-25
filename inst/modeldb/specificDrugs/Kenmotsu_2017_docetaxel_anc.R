Kenmotsu_2017_docetaxel_anc <- function() {
  description <- "Exposure-toxicity model for the percent decrease in absolute neutrophil count (ANC) at nadir during the first 3-week course of docetaxel 60 mg/m^2 in Japanese adults with non-small cell lung cancer (Kenmotsu 2017). Sigmoid Emax function of the unbound docetaxel AUC (fu * AUC, with AUC = Dose / CL): Emax 84.5%, Hill coefficient 3.72, and an EC50 of 0.101 ug*h/mL when baseline ANC <= 4341 /uL versus 0.142 ug*h/mL above it. The individual clearance comes from the embedded companion population PK model (modellib('Kenmotsu_2017_docetaxel'): two compartments, CL proportional to BSA and decreasing with alpha-1 acid glycoprotein, Vc decreasing with albumin) and the per-subject measured unbound fraction is supplied as the FU covariate. One observation per subject and course, no IIV on the PD parameters, exponential (proportional) residual error."
  reference <- "Kenmotsu H, Imamura CK, Ono A, Omori S, Nakashima K, Wakuda K, Taira T, Naito T, Murakami H, Takahashi T, Tanigawara Y. The effects of advanced age and serum alpha1-acid glycoprotein on docetaxel unbound exposure and dose-limiting toxicity in cancer patients. Br J Clin Pharmacol. 2017;83(11):2416-2425. doi:10.1111/bcp.13354. PMCID: PMC5651322. Correction: Br J Clin Pharmacol. 2018;84(10):2449. doi:10.1111/bcp.13694 (relabels the Figure 2 y-axes as clearance; no model value changed)."
  vignette <- "Kenmotsu_2017_docetaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL", response = "% decrease from pretreatment ANC")

  compartmentData <- list(
    central = list(analyte = "docetaxel", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "docetaxel", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Clearance is directly proportional to BSA (Table 2: CL = theta1 * BSA * (AAG/121)^theta2). BSA formula not stated. Cohort range 1.19-2.09 m^2 (Table 1).",
      source_name = "BSA"
    ),
    AAG = list(
      description = "Serum alpha-1 acid glycoprotein concentration at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Source unit mg/dL with clearance reference 121 mg/dL (Table 2); 1.21 g/L in canonical units gives the identical ratio. Cohort range 0.50-2.64 g/L (Table 1). AAG enters only through clearance (and hence AUC); the measured fu is supplied separately through FU.",
      source_name = "AAG"
    ),
    ALB = list(
      description = "Serum albumin concentration at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Source unit g/dL with central-volume reference 3.7 g/dL (Table 2); 37 g/L in canonical units gives the identical ratio. Affects Vc only, so it does not change AUC = Dose / CL or the nadir prediction.",
      source_name = "ALB"
    ),
    FU = list(
      description = "Measured unbound fraction of docetaxel in serum",
      units = "unitless fraction",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject value measured by equilibrium dialysis of a serum sample drawn 10-60 min after the end of the infusion (Methods). Supplied as a unitless FRACTION (the paper reports percent: median 8.1%, range 3.1-15.6%, Figure 3A; mean 7.8 +/- 2.7% for age < 75 and 8.9 +/- 2.1% for age >= 75 years, Table 3). The paper does not model fu from AAG or albumin (fu correlated only weakly with AAG, R^2 = 0.149, Figure 3C), so it must be supplied as data.",
      source_name = "fu"
    ),
    NEUT = list(
      description = "Absolute neutrophil count at baseline (pretreatment)",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Enters only as the dichotomy baseline ANC <= 4341 /uL (the cohort median, Table S2 footnote) versus > 4341 /uL, which selects the EC50 (Table 2). Cohort medians 4461 /uL (range 2132-9650, age < 75) and 4001 /uL (range 2088-9911, age >= 75), Table 1.",
      source_name = "Baseline ANC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "34-84 years",
    sex_female_pct = 15.7,
    race_ethnicity = c(Asian = 100),
    disease_state = "Non-small cell lung cancer, docetaxel monotherapy",
    dose_range = "60 mg/m^2 as a 1-h intravenous infusion every 3 weeks (first course)",
    regions = "Japan (Shizuoka Cancer Center; UMIN000014701)",
    notes = "49 of the 51 PK patients were evaluable for neutropenia (33 aged < 75 years and 16 aged >= 75 years); two were excluded for prophylactic G-CSF (Results). Percent decrease in ANC at nadir = (pretreatment ANC - nadir ANC) / pretreatment ANC * 100 (Equation 3) over the first 3-week course. Age and sex were tested on Emax and EC50 and not retained (Table S2)."
  )

  ini({
    # Pharmacokinetics -- identical to modellib('Kenmotsu_2017_docetaxel').
    lcl <- log(14.5)
    label("Total docetaxel clearance per m^2 body surface area at AAG = 121 mg/dL, theta1 (L/h/m^2)") # Table 2, theta1 = 14.5 (SE 0.516)
    e_aag_cl <- -0.495
    label("Power exponent of AAG/121 mg/dL on CL, theta2 (unitless)") # Table 2, theta2 = -0.495 (SE 0.0816)
    lvc <- log(8.89)
    label("Central volume at ALB = 3.7 g/dL, theta3 (L)") # Table 2, theta3 = 8.89 (SE 0.403)
    e_alb_vc <- -0.550
    label("Power exponent of ALB/3.7 g/dL on Vc, theta4 (unitless)") # Table 2, theta4 = -0.550 (SE 0.185)
    lvss <- log(129)
    label("Volume of distribution at steady state, Vss (L)") # Table 2, Vss = 129 (SE 8.28)
    lq <- log(11.0)
    label("Intercompartmental clearance, Q (L/h)") # Table 2, Q = 11.0 (SE 0.542)

    # Exposure-toxicity relationship for neutropenia (Equation 4).
    lemax <- log(84.5)
    label("Maximum percent decrease in ANC at nadir, Emax (%)") # Table 2, Emax = 84.5% (SE 2.13)
    lec50 <- log(0.142)
    label("fu*AUC at half-maximal effect when baseline ANC > 4341 /uL, EC50 (ug*h/mL)") # Table 2, EC50 baseline ANC > 4341 = 0.142 (SE 0.00267)
    lec50_neutlow <- log(0.101)
    label("fu*AUC at half-maximal effect when baseline ANC <= 4341 /uL, EC50 (ug*h/mL)") # Table 2, EC50 baseline ANC <= 4341 = 0.101 (SE 0.0159)
    lhill <- log(3.72)
    label("Sigmoidicity (shape) coefficient, gamma (unitless)") # Table 2, shape parameter gamma = 3.72 (SE 0.518)

    etalcl ~ 0.0426 # Table 2, omega2 CL = 0.0426 (SE 0.00903)
    etalvc ~ 0.00885 # Table 2, omega2 Vc = 0.00885 (SE 0.00923)
    etalvss ~ 0.0176 # Table 2, omega2 Vss = 0.0176 (SE 0.00987)

    propSd <- 0.12806
    label("Proportional (exponential) residual error on percent ANC decrease at nadir (fraction)") # Table 2, exposure-toxicity residual variability 0.0164 read as a variance -> sqrt(0.0164) = 0.12806
  })

  model({
    cl <- exp(lcl + etalcl) * BSA * (AAG / 1.21)^e_aag_cl
    vc <- exp(lvc + etalvc) * (ALB / 37)^e_alb_vc
    vss <- exp(lvss + etalvss)
    q <- exp(lq)
    vp <- vss - vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Methods: AUC = Dose / CL for the course's dose (mg / (L/h) = ug*h/mL),
    # and the exposure driving neutropenia is the unbound AUC fu * AUC.
    auc_course <- podo(central) / cl
    fu_auc <- FU * auc_course

    emax <- exp(lemax)
    hill <- exp(lhill)
    # Table 2 / Figure 4A: separate EC50 for baseline ANC <= 4341 /uL.
    ec50 <- exp(lec50)
    if (NEUT <= 4341) ec50 <- exp(lec50_neutlow)

    # Equation 4
    nadir_anc_pct_decrease <- emax * fu_auc^hill / (ec50^hill + fu_auc^hill)
    nadir_anc_pct_decrease ~ prop(propSd)
  })
}
