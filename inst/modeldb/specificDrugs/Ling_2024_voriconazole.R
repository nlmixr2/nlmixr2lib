Ling_2024_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for intravenous and oral voriconazole in Chinese adults with invasive fungal infections (Ling 2024); C-reactive protein enters clearance as an exponential inflammation effect that is switched off in CYP2C19 poor metabolizers, alongside CYP2C19 phenotype, age, serum albumin and sex effects on clearance and a power-form body-weight effect on volume of distribution."
  reference <- "Ling J, Yang X, Dong L, Jiang Y, Zou S, Hu N. Influence of C-reactive protein on the pharmacokinetics of voriconazole in relation to the CYP2C19 genotype: a population pharmacokinetics analysis. Front Pharmacol. 2024;15:1455721. doi:10.3389/fphar.2024.1455721"
  vignette <- "Ling_2024_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot   = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRP = list(
      description        = "C-reactive protein concentration (standard clinical assay), time-varying: Ling 2024 recorded CRP on the same day as each therapeutic-drug-monitoring sample",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL as an exponential effect scaled by the cohort-median CRP of 59 mg/L",
        "(Table 1 reports median 58.85 mg/L, rounded to 59 in the published final-model equation):",
        "exp(e_crp_cl * CRP / 59). Note this is a median-SCALED but not median-CENTERED form -- the",
        "printed typical CL of 3.83 L/h is therefore the value at CRP = 0, not at the median CRP.",
        "At the median CRP the typical CL is 3.83 * exp(-0.155) = 3.28 L/h. The uncentered reading is",
        "confirmed by the paper's own Monte Carlo target-attainment percentages; see the model vignette.",
        "The effect is multiplied by (1 - CYP2C19_PM) because Ling 2024 estimated the CRP effect only in",
        "CYP2C19 normal and intermediate metabolizers -- the separately-fitted PM exponent was -0.0172",
        "(a 3% CL change over CRP 1 to 100 mg/L) and the term was dropped for PM patients.",
        "Cohort CRP: mean 77.10, SD 68.74, median 58.85, range 0.9-306.6 mg/L (Ling 2024 Table 1)."
      ),
      source_name        = "CRP"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = paste(
        "Ling 2024 assigned CPIC phenotypes from fluorescence in-situ hybridization genotyping.",
        "IM genotypes pooled by Ling 2024 were *1/*2, *1/*3 and *2/*17 (Methods, 'Genotyping and phenotype",
        "assignment'). 72 of 167 patients (43.1%) were IM (Table 1). The paper's reference category is the",
        "normal metabolizer (*1/*1) phenotype, which coincides with the canonical EM/UM implicit reference",
        "because no ultrarapid (*17/*17) or rapid (*1/*17) metabolizers were enrolled -- so no",
        "reparameterization was needed. The published effect is a multiplicative ratio applied as",
        "e_cyp2c19_im_cl^CYP2C19_IM (Ling 2024 final-model equation)."
      ),
      source_name        = "IM"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal metabolizer; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = paste(
        "Companion to CYP2C19_IM; see those notes for the reference-category rationale. PM genotypes pooled",
        "by Ling 2024 were *2/*2, *2/*3 and *3/*3. 29 of 167 patients (17.4%) were PM (Table 1).",
        "CYP2C19_PM appears twice in model(): once as the multiplicative phenotype ratio on CL",
        "(e_cyp2c19_pm_cl^CYP2C19_PM) and once as the (1 - CYP2C19_PM) switch that removes the CRP",
        "inflammation effect for poor metabolizers."
      ),
      source_name        = "PM"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL scaled to the cohort-median age of 71 years (Ling 2024 Table 1):",
        "(AGE / 71)^e_age_cl. Cohort age mean 68.87, SD 14.87, median 71, range 16-97 years. The abstract",
        "states patients aged >= 16 years while the Methods state >= 18 years; Table 1's range starts at 16."
      ),
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported by Ling 2024 in SI units (g/L), matching the canonical unit -- no conversion is applied",
        "in model(). Power-form effect on CL scaled to the cohort-median albumin of 34.8 g/L",
        "(Ling 2024 Table 1): (ALB / 34.8)^e_alb_cl. Cohort albumin mean 36.36, SD 8.54, median 34.8,",
        "range 18.3-75.1 g/L. The positive exponent means low albumin lowers clearance, consistent with",
        "the paper's Discussion recommendation to monitor for toxicity in hypoalbuminaemic patients."
      ),
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Sex indicator, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Ling 2024 codes gender as M = 0, F = 1 in the final-model equation, matching the canonical SEXF",
        "orientation exactly -- no value transformation is needed. The effect is a multiplicative ratio",
        "applied as e_sexf_cl^SEXF, so women have 1.41-fold higher voriconazole clearance than men.",
        "Cohort composition 119 male / 48 female (Ling 2024 Table 1)."
      ),
      source_name        = "gender"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on V scaled to the cohort-median weight of 65 kg (Ling 2024 Table 1):",
        "(WT / 65)^e_wt_vc. Cohort weight mean 64.48, SD 12.24, median 65, range 37-100 kg. The estimated",
        "exponent 2.21 is far above the allometric-theory value of 1 and is imprecisely estimated",
        "(RSE 28.3%, bootstrap 95% CI 0.705-3.408); Ling 2024 attributes the weak identifiability of V to",
        "the trough-only sampling design, which gave 65.7% eta-shrinkage on V. Do not extrapolate this",
        "exponent outside the observed 37-100 kg range."
      ),
      source_name        = "WT"
    )
  )

  # Covariates Ling 2024 screened during stepwise covariate modelling but did
  # not retain in the final model. Documentation only -- these names are
  # deliberately absent from model(). Uric acid (mean 196.97, median 156.8,
  # range 27.4-1917 umol/L; Ling 2024 Table 1) was also screened and dropped
  # but has no canonical register entry, so it is recorded in population$notes
  # rather than here.
  covariatesDataExcluded <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor use",
      units       = "(binary)",
      type        = "binary",
      notes       = "117 of 167 patients (70.1%) received a PPI (Ling 2024 Table 1). No significant effect on voriconazole PK was found; the Discussion attributes this to most patients receiving rabeprazole or pantoprazole rather than omeprazole."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic corticosteroid use",
      units       = "(binary)",
      type        = "binary",
      notes       = "43 of 167 patients (25.7%) received a glucocorticoid (Ling 2024 Table 1). Screened as a co-medication covariate and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Mean 63.98, median 34.7, range 7.8-1211.6 U/L (Ling 2024 Table 1). Screened as a liver-function covariate and not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Mean 44.42, median 25.35, range 2.6-629.1 U/L (Ling 2024 Table 1). Screened as a liver-function covariate and not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Mean 26.65, median 13.35, range 1.7-514.8 umol/L (Ling 2024 Table 1). Screened as a liver-function covariate and not retained."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Mean 97.96, median 96, range 64-152 (Ling 2024 Table 1). Table 1 prints the unit as mmol/L, but the magnitudes are unambiguously g/L (the SI mmol/L scale for hemoglobin runs about 4-10). Screened as a complete-blood-count covariate and not retained."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Mean 174.25, median 165, range 4-624 x 10^9/L (Ling 2024 Table 1). Screened as a complete-blood-count covariate and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Mean 106.50, median 78, range 2.77-641 umol/L (Ling 2024 Table 1). Screened as a renal-function covariate and not retained, consistent with voriconazole being cleared by hepatic metabolism."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 167L,
    n_studies      = 1L,
    n_observations = 232L,
    age_range      = "16-97 years",
    age_median     = "71 years",
    age_mean       = "68.87 +/- 14.87 years",
    weight_range   = "37-100 kg",
    weight_median  = "65 kg",
    weight_mean    = "64.48 +/- 12.24 kg",
    sex_female_pct = 28.7,
    race_ethnicity = c(Chinese = 100),
    cyp2c19_phenotype = c(NM_pct = 39.5, IM_pct = 43.1, PM_pct = 17.4, RM_pct = 0, UM_pct = 0),
    co_medication  = c(ProtonPumpInhibitor_pct = 70.1, Corticosteroid_pct = 25.7),
    disease_state  = "Adult inpatients with proven, probable or possible invasive fungal infection treated with voriconazole and undergoing therapeutic drug monitoring. Patients receiving other antifungals or co-medications known to alter voriconazole PK were excluded.",
    dose_range     = "Intravenous or oral voriconazole twice daily. 57 patients received an intravenous loading dose of 600 or 400 mg twice daily followed by an intravenous maintenance dose of 300 or 200 mg twice daily; 4 patients received an oral loading dose of 400 or 300 mg twice daily followed by 200 mg twice daily orally; 106 patients received 200 mg twice daily intravenously or orally with no loading dose. Intravenous infusion rate was kept below 3 mg/kg/h. Oral doses were taken 1 h before or after a meal.",
    regions        = "Single center: The First People's Hospital of Changzhou / The Third Affiliated Hospital of Soochow University, Changzhou, Jiangsu, China.",
    notes          = paste(
      "Retrospective single-center study, October 2020 - June 2023 (the Methods 'Patients' section states",
      "March 2020 - August 2023 for the study period; the abstract gives the October 2020 - June 2023",
      "data-collection window). Ethics approval No. 2023-038. 232 steady-state trough concentrations from",
      "167 patients, measured by HPLC-MS/MS over a validated 0.1-20 ug/mL range. Samples were drawn within",
      "30 min before the next dose, on treatment day 3 with a loading dose or day 5 without. Fewer than 5%",
      "of concentrations were below the quantification limit and were discarded (M1 method). Observed",
      "trough concentrations: mean 4.67, SD 2.66, median 4.45, range 0.25-16.75 ug/mL; 59.1% (137/232) fell",
      "in the 0.5-5.0 ug/mL therapeutic window of the Chinese voriconazole practice guidelines. Uric acid",
      "(mean 196.97, median 156.8, range 27.4-1917 umol/L) was screened as a covariate and not retained.",
      "Baseline demographics per Ling 2024 Table 1; final-model estimates and 1000-replicate bootstrap",
      "(898 successful) per Ling 2024 Table 2."
    )
  )

  ini({
    # Structural parameters. The reference subject for the typical-value
    # equation is a CYP2C19 normal metabolizer (both phenotype indicators
    # 0) male (SEXF = 0) at AGE = 71 years, ALB = 34.8 g/L, WT = 65 kg --
    # the cohort medians of Ling 2024 Table 1 -- and, because the CRP term
    # is median-scaled rather than median-centered, at CRP = 0 mg/L.

    # Absorption: ka fixed at 1.1/h because the retrospective dataset held
    # trough concentrations only and contained no absorption-phase data.
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h)")  # Ling 2024 Methods 'Population pharmacokinetic modeling': "ka was fixed at 1.1 h-1, as previously reported by Pascual et al. (2012)"; Table 2 lists ka as 1.1 (fixed) in both the base and final model

    lcl <- log(3.83); label("Clearance at the reference covariate values (L/h)")  # Ling 2024 Table 2 final model: CL = 3.83 L/h (RSE 9.5%, bootstrap median 3.84, 95% CI 3.17-4.77)

    lvc <- log(134); label("Volume of distribution at WT = 65 kg (L)")  # Ling 2024 Table 2 final model: V = 134 L (RSE 12.0%, bootstrap median 134, 95% CI 99.7-169.9)

    # Oral bioavailability. Estimated, not fixed; interindividual
    # variability on F was tested and rejected (shrinkage 90.63%).
    lfdepot <- log(0.965); label("Oral bioavailability (fraction)")  # Ling 2024 Table 2 final model: F = 0.965 (RSE 7.4%, bootstrap median 0.965, 95% CI 0.829-1.110)

    # Covariate effects on CL. The CYP2C19 phenotype and sex effects are
    # multiplicative ratios applied as ratio^indicator; age and albumin
    # are power terms on the median-normalized covariate; CRP is an
    # exponential term on the median-scaled covariate.
    e_cyp2c19_im_cl <- 0.794; label("CL ratio for CYP2C19 intermediate vs normal metabolizer (unitless)")  # Ling 2024 Table 2 final model, "IM on CL" = 0.794 (RSE 8.0%, bootstrap 95% CI 0.657-0.905)
    e_cyp2c19_pm_cl <- 0.635; label("CL ratio for CYP2C19 poor vs normal metabolizer (unitless)")  # Ling 2024 Table 2 final model, "PM on CL" = 0.635 (RSE 12.0%, bootstrap 95% CI 0.472-0.781)
    e_crp_cl        <- -0.155; label("Exponential CRP coefficient on CL, per unit of CRP/59 (unitless)")  # Ling 2024 Table 2 final model, "CRP on CL" = -0.155 (RSE 25.2%, bootstrap 95% CI -0.273 to -0.076)
    e_age_cl        <- -0.582; label("Power exponent for AGE on CL (unitless)")  # Ling 2024 Table 2 final model, "Age on CL" = -0.582 (RSE 25.3%, bootstrap 95% CI -0.928 to -0.307)
    e_alb_cl        <- 0.644; label("Power exponent for ALB on CL (unitless)")  # Ling 2024 Table 2 final model, "ALB on CL" = 0.644 (RSE 26.4%, bootstrap 95% CI 0.300-1.100)
    e_sexf_cl       <- 1.410; label("CL ratio for female vs male (unitless)")  # Ling 2024 Table 2 final model, "Gender on CL" = 1.410 (RSE 7.2%, bootstrap 95% CI 1.203-1.619)

    # Covariate effect on V.
    e_wt_vc <- 2.210; label("Power exponent for WT on V (unitless)")  # Ling 2024 Table 2 final model, "WT on V" = 2.210 (RSE 28.3%, bootstrap 95% CI 0.705-3.408)

    # IIV. Ling 2024 Methods states an exponential interindividual model,
    # Pj = P_typical * exp(eta_j) with eta ~ N(0, omega^2). Table 2 reports
    # the final-model IIV as 38.9% for CL and 45.2% for V. Using the usual
    # NONMEM reporting convention CV% ~= sqrt(omega^2) * 100, the internal
    # variances are 0.389^2 = 0.151321 and 0.452^2 = 0.204304. This
    # SD-not-variance reading is corroborated by reproducing the paper's
    # published Monte Carlo target-attainment percentages; see the vignette.
    etalcl ~ 0.151321  # Ling 2024 Table 2 final model: IIV on CL = 38.9% (RSE 20.4%, eta-shrinkage 16.2%, bootstrap median 37.7%, 95% CI 28.4-46.7%); var = 0.389^2
    etalvc ~ 0.204304  # Ling 2024 Table 2 final model: IIV on V  = 45.2% (RSE 31.1%, eta-shrinkage 65.7%, bootstrap median 46.4%, 95% CI 16.5-72.2%); var = 0.452^2

    # Residual error. Ling 2024 Methods specifies a combined exponential-
    # plus-additive residual model, Cij = Chat_ij * exp(eps1) + eps2, but
    # Table 2 labels and reports the first component as a proportional
    # error ("Prop (%)"). The two forms agree to first order at this
    # magnitude (exp(0.147) - 1 = 0.158), and the combined additive-plus-
    # proportional form is the one nlmixr2 supports directly.
    propSd <- 0.147; label("Proportional residual error (fraction)")  # Ling 2024 Table 2 final model: Prop = 14.7% (RSE 27.4%, eps-shrinkage 35.3%, bootstrap median 14.3%, 95% CI 7.82-18.6%)
    addSd  <- 0.58; label("Additive residual error (ug/mL)")  # Ling 2024 Table 2 final model: Add = 0.58 ug/mL (RSE 30.2%, eps-shrinkage 35.3%, bootstrap median 0.55, 95% CI 0.23-0.74)
  })

  model({
    # Individual PK parameters. The clearance equation reproduces the
    # Ling 2024 published final model:
    #
    #   CL = 3.83 * 0.794^(IM) * 0.635^(PM)
    #        * [exp(-0.155 * CRP/59)]^(NM=1, IM=1, PM=0)
    #        * (AGE/71)^-0.582 * (ALB/34.8)^0.644 * 1.41^(gender: M=0, F=1)
    #
    # The paper writes the CRP switch as an indicator exponent that is 1
    # for normal and intermediate metabolizers and 0 for poor
    # metabolizers. Because the cohort contained only NM, IM and PM
    # phenotypes, that indicator is exactly (1 - CYP2C19_PM).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      e_cyp2c19_im_cl^CYP2C19_IM *
      e_cyp2c19_pm_cl^CYP2C19_PM *
      exp(e_crp_cl * (CRP / 59) * (1 - CYP2C19_PM)) *
      (AGE / 71)^e_age_cl *
      (ALB / 34.8)^e_alb_cl *
      e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc

    # One-compartment disposition with first-order oral absorption.
    # Intravenous doses bypass the depot and enter central directly.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Oral bioavailability applies only to doses entering via the depot.
    f(depot) <- exp(lfdepot)

    # Observation. Amounts in mg and volume in L give mg/L = ug/mL, the
    # scale on which the paper reports trough concentrations.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
