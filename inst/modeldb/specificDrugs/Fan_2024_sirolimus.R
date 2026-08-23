Fan_2024_sirolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral sirolimus (rapamycin) in Chinese children with vascular anomalies, developed from routine therapeutic-drug-monitoring trough concentrations (Fan 2024). Body weight is the only covariate retained in the final model, entering as separately estimated allometric exponents on apparent clearance (1.23) and on apparent volume of distribution (1.62), each normalized to the 16 kg cohort median body weight. The absorption rate constant was held at 0.485 per hour from prior sirolimus studies in Chinese children because no absorption-phase samples were available; inter-individual variability was estimated on CL/F only (V/F variability was not estimable from trough-only data) and residual variability is additive."
  reference <- paste(
    "Fan L, Guo HL, Zhao YT, Li Y, Wang WJ, Huang J, Hu YH, Zou JJ, Chen F.",
    "Population pharmacokinetic study in children with vascular anomalies:",
    "body weight as a key variable in predicting the initial dose and dosing",
    "frequency of sirolimus. Front Pharmacol. 2024;15:1457614.",
    "doi:10.3389/fphar.2024.1457614",
    sep = " "
  )
  vignette <- "Fan_2024_sirolimus"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Fan 2024 Methods 2.2: whole blood
  # (1-2 mL, EDTA K2) was assayed by EMIT for sirolimus, so the observed
  # matrix is whole blood rather than plasma (the paper's prose uses
  # "plasma" loosely in places; the sampling description is unambiguous).
  compartmentData <- list(
    depot   = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the Fan 2024 final model. Enters as a power function normalized to the 16 kg cohort median body weight (BWmedian), with separately estimated exponents on each parameter: CL/F = 4.06 * (WT/16)^1.23 (Equation 22) and V/F = 155 * (WT/16)^1.62 (Equation 23). Both exponents were estimated rather than fixed at the theoretical allometric 0.75 / 1 values; the paper's Model I (simple exponential, Equations 6-7) achieved a lower objective function value than the fixed-exponent Model II (Equations 8-9) and three maturation-model alternatives (Supplementary Table S4). Age was screened but excluded because of collinearity with body weight (r = 0.86). Studied range 3.3-65 kg.",
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the size / maturation model series (Fan 2024 Methods 2.3.2, Models III and V) but not retained: body weight and age were strongly correlated (r = 0.86), so age was dropped to avoid collinearity and parameter-estimation instability (Discussion). Studied range 0.08-12 years."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in the demographic table (24 male / 25 female) but not screened as a covariate on CL/F in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened during stepwise forward inclusion / backward exclusion but not significant on CL/F (Fan 2024 Results 3.2; Supplementary Table S4)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a marker of hepatic function but not significant on CL/F. The Discussion attributes this to an insufficient number of patients with liver impairment in the retrospective cohort."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a marker of hepatic function but not significant on CL/F (Fan 2024 Results 3.2)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened because sirolimus distributes predominantly into red blood cells, but not significant on CL/F. The Discussion attributes this to blood-cell indices fluctuating only within normal ranges in this cohort."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened alongside the other red-cell indices but not significant on CL/F."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49L,
    n_studies      = 1L,
    n_observations = 134L,
    age_range      = "0.08-12 years",
    age_median     = "3.5 years",
    weight_range   = "3.3-65 kg",
    weight_median  = "16 kg",
    sex_female_pct = 51.0,
    race_ethnicity = "Chinese (single-centre cohort in Nanjing, China); not reported by subgroup.",
    disease_state  = "Children with vascular anomalies (including tufted angioma, kaposiform hemangioendothelioma, and lymphatic and venous malformations) receiving oral sirolimus.",
    dose_range     = "0.018-0.152 mg/kg/day oral sirolimus; the common initial regimen was 0.08 mg/kg/day given at a dosing interval of either 12 or 24 h.",
    regions        = "China (single centre: Children's Hospital of Nanjing Medical University).",
    assay          = "Enzyme multiplied immunoassay technique (EMIT 2000, Siemens) on whole blood, calibration range 3.5-30 ng/mL; quality-control deviations over the collection period ranged from -13.2 to 14.8 percent.",
    sampling       = "All 134 observations are steady-state trough concentrations drawn 30 min before the next maintenance dose, at least 7 days after starting sirolimus. No absorption-phase or peak samples were collected, which is why Ka was not estimable and V/F carries no inter-individual variability.",
    target_range   = "Trough concentration target of 5-15 ng/mL used for the Monte Carlo dose-optimization simulations.",
    notes          = "Retrospective TDM cohort collected between July 2017 and April 2022 (Fan 2024 Table 1). Exclusion criteria were concentrations beyond the assay detection limit and ongoing serious infections or multiple organ injuries. Sex: 24 male / 25 female. Genotyping for CYP3A4, CYP3A5, mTOR, ABCB1, ABCC2, CYP3A7, POR, IL10, IL18, SUMO4, NR1I2 and TCF7L2 variants was performed and screened, but no single nucleotide polymorphism was retained in the final model."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, reported at the 16 kg cohort-median body
    # weight (the BWmedian normalization constant in Equations 22-23).
    # Final-model point estimates from Fan 2024 Table 2.
    # ------------------------------------------------------------------
    lcl <- log(4.06); label("Apparent clearance CL/F at the 16 kg reference body weight (L/h)")                  # Table 2: CL/F = 4.06 L/h (4% RSE); Equation 22
    lvc <- log(155);  label("Apparent volume of distribution V/F at the 16 kg reference body weight (L)")        # Table 2: V/F = 155 L (26% RSE); Equation 23

    # Ka was not estimated. Methods 2.3.1: "Due to the absence of
    # observations during the absorption phase, absorption rate constant
    # (Ka) was established at 0.485 h-1 according to previously reports in
    # the literature". Table 2 reports Ka = 0.485 with "fixed" in place of
    # an RSE and no bootstrap confidence interval.
    lka <- fixed(log(0.485)); label("First-order absorption rate constant Ka (1/h; literature value from prior sirolimus studies in Chinese children)")  # Table 2: Ka = 0.485 1/h

    # ------------------------------------------------------------------
    # Allometric exponents on body weight. Both were ESTIMATED (Model I,
    # the "simple exponential model" of Equations 6-7), not fixed at the
    # theoretical 0.75 / 1 of Model II, and both carry an RSE and a
    # bootstrap confidence interval in Table 2 -- hence no fixed() wrapper.
    # ------------------------------------------------------------------
    e_wt_cl <- 1.23; label("Allometric exponent on body weight for CL/F (unitless)")  # Table 2: m = 1.23 (8% RSE, bootstrap 0.98-1.43); Equation 22
    e_wt_vc <- 1.62; label("Allometric exponent on body weight for V/F (unitless)")   # Table 2: n = 1.62 (13% RSE, bootstrap 1.09-2.07); Equation 23

    # ------------------------------------------------------------------
    # Inter-individual variability. Methods 2.3.1 Equation 1 specifies the
    # exponential model Pi = TV(P) * exp(eta_i), so the reported omega is
    # a standard deviation on the natural-log scale. Table 2 gives
    # omega_CL = 25.94% (32% RSE, 17% shrinkage), encoded here as the
    # variance 0.2594^2 = 0.06729.
    #
    # Results 3.2 is explicit that there is no eta on V/F: "Since only
    # Ctrough were applied in the study, the inter-individual variability
    # of V/F was not estimated." No etalvc is declared, matching the
    # published OMEGA structure.
    # ------------------------------------------------------------------
    etalcl ~ 0.06729  # omega_CL = 25.94% -> 0.2594^2; Table 2

    # ------------------------------------------------------------------
    # Residual unexplained variability. Results 3.2: "The residual
    # variability was best described by an additive error model
    # (Equation 2)", i.e. Y = IPRED + eps. Table 2 reports sigma in
    # ng/mL, so the value is a standard deviation on the concentration
    # scale rather than a variance.
    # ------------------------------------------------------------------
    addSd <- 3.41; label("Additive residual error (ng/mL)")  # Table 2: sigma = 3.41 ng/mL (14% RSE, 13% shrinkage)
  })

  model({
    # Reference body weight: the cohort median (BWmedian) that Fan 2024
    # normalizes both power functions to (Table 1 median WT = 16 kg;
    # Equations 22-23 are written with a literal 16 in the denominator).
    ref_wt <- 16

    # Individual parameters. Both allometric exponents were estimated
    # separately, so CL/F and V/F do not share a common exponent.
    #   CL/F = 4.06 * (WT/16)^1.23   (Equation 22)
    #   V/F  = 155  * (WT/16)^1.62   (Equation 23)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc)          * (WT / ref_wt)^e_wt_vc

    kel <- cl / vc

    # ODE system: one-compartment disposition with first-order absorption.
    # Bioavailability is not separately identifiable from oral-only data,
    # so CL and V are the apparent quantities CL/F and V/F and no f(depot)
    # term is applied (equivalent to F fixed at 1).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation: whole-blood sirolimus concentration in ng/mL. Doses are
    # in mg and volumes in L, so central / vc is in mg/L; multiply by 1000
    # to convert to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
