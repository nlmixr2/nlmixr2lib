Jin_2020_amikacin <- function() {
  description <- "Two-compartment population PK model for intravenous amikacin in Korean adults with nontuberculous mycobacterial pulmonary disease (Jin 2020), with a power effect of MDRD eGFR on clearance (reference 91.1 mL/min/1.73 m^2) and a power effect of body weight on central volume (reference 51.1 kg)."
  reference <- paste(
    "Jin X, Oh J, Cho JY, Lee S, Rhee SJ.",
    "Population Pharmacokinetic Analysis of Amikacin for Optimal Pharmacotherapy",
    "in Korean Patients with Nontuberculous Mycobacterial Pulmonary Disease.",
    "Antibiotics (Basel). 2020;9(11):784.",
    "doi:10.3390/antibiotics9110784.",
    sep = " "
  )
  vignette <- "Jin_2020_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the Modification of Diet in Renal Disease (MDRD) equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column eGFR. Jin 2020 Table 1 footnote and Methods 4.1: eGFR =",
        "175 x SCr^-1.154 x Age^-0.203 (x 0.742 if female), i.e. the IDMS-traceable",
        "4-variable MDRD equation, in mL/min/1.73 m^2. Enters CL as the power",
        "(CRCL / 91.1)^0.229, where 91.1 mL/min/1.73 m^2 is the population median",
        "(Results 2.2). Cohort mean 95.7 (SD 36.4), range 26.8-297.8 (Table 1).",
        "Recorded immediately before TDM sampling; hemodialysis patients were",
        "excluded."
      ),
      source_name = "eGFR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V1 as the power (WT / 51.1)^0.702, where 51.1 kg is the population",
        "median (Results 2.2). Cohort mean 51.9 kg (SD 11.2), range 29.9-79.8 kg",
        "(Table 1). Recorded immediately before TDM sampling."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a continuous covariate (linear and power forms; Methods 4.2) but not retained in the final model. Cohort range 25-85 years (Table 1)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened as a continuous covariate (Methods 4.2) but not retained in the final model. Cohort range 146.1-177.0 cm (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened as a continuous covariate (Methods 4.2) but not retained; renal function enters only through MDRD eGFR. Cohort range 0.3-1.9 mg/dL (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened as a continuous covariate (Methods 4.2) but not retained in the final model. Cohort range 1.6-4.8 g/dL (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a categorical covariate (Methods 4.2; source coded IND = 1 for male) but not retained in the final model. 51 of 70 patients (73%) were female."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 1L,
    n_observations = 848L,
    age_range = "25-85 years (mean 66.0, SD 11.5)",
    weight_range = "29.9-79.8 kg (mean 51.9, SD 11.2; median 51.1)",
    height_range = "146.1-177.0 cm (mean 160.5, SD 8.2)",
    bmi_range = "12.8-33.3 kg/m^2 (mean 20.2, SD 3.7)",
    sex_female_pct = 72.9,
    race_ethnicity = "Korean (single-centre, Seoul National University Hospital)",
    disease_state = "Nontuberculous mycobacterial pulmonary disease (NTM-PD) treated with intravenous amikacin with therapeutic drug monitoring; hemodialysis patients excluded.",
    renal_function = "MDRD eGFR 26.8-297.8 mL/min/1.73 m^2 (mean 95.7, SD 36.4; median 91.1); serum creatinine 0.3-1.9 mg/dL (mean 0.7)",
    albumin_range = "1.6-4.8 g/dL (mean 3.6)",
    dose_range = "Intravenous infusion over 30 min or 1 h at dosing intervals of 8, 12, 24 or 48 h (routine clinical dosing; guideline 10-15 mg/kg once daily).",
    regions = "Republic of Korea (Seoul National University Hospital), retrospective, 1 December 2009 to 1 December 2019.",
    notes = paste(
      "Demographics from Jin 2020 Table 1. Sparse TDM sampling (mostly a peak",
      "within 60 min after the end of infusion and a trough within 30 min before",
      "the next dose); 55 of the 848 samples fell outside those windows.",
      "Concentrations below the limit of quantitation were excluded. Fit in",
      "NONMEM 7.4.3 with FOCE-I on log-transformed concentrations; evaluated by",
      "GOF, pcVPC and a 1000-replicate bootstrap."
    )
  )

  ini({
    # Structural parameters: Jin 2020 Table 2, 'Final Model' column, and the
    # final-model equations printed in the Table 2 footnote:
    #   CL = 3.52 x (eGFR/91.1)^0.229 x exp(eta); V1 = 14.4 x (WT/51.1)^0.702 x exp(eta);
    #   V2 = 14.2; Q = 0.464.
    lcl <- log(3.52); label("Clearance at eGFR 91.1 mL/min/1.73 m^2 (L/h)") # Table 2 final model: CL = 3.52 L/h (RSE 5%)
    lvc <- log(14.4); label("Central volume V1 at WT 51.1 kg (L)") # Table 2 final model: V1 = 14.4 L (RSE 4%)
    lvp <- log(14.2); label("Peripheral volume V2 (L)") # Table 2 final model: V2 = 14.2 L (RSE 34%)
    lq <- log(0.464); label("Intercompartmental clearance Q (L/h)") # Table 2 final model: Q = 0.464 L/h (RSE 25%)

    # Covariate effects, power form (Methods 4.2 Eq. 2, P = theta1 x (COV/MED)^theta2).
    e_crcl_cl <- 0.229; label("Power exponent of (eGFR / 91.1) on CL (unitless)") # Table 2 final model: eGFR effect on CL = 0.229 (RSE 32%)
    e_wt_vc <- 0.702; label("Power exponent of (WT / 51.1) on V1 (unitless)") # Table 2 final model: WT effect on V1 = 0.702 (RSE 17%)

    # Between-subject variability, exponential model (Methods 4.2), estimated as
    # an OMEGA BLOCK on CL and V1 (Results 2.2). Table 2 reports the IIV as CV%
    # (CL 27.9%, V1 18%), converted to log-scale variances as omega^2 =
    # log(CV^2 + 1); the covariance is taken as printed, -0.0135 (bootstrap
    # median -0.012). Implied correlation -0.0135 / sqrt(0.07494 x 0.03189) = -0.276.
    etalcl + etalvc ~ c(0.07494, -0.0135, 0.03189) # Table 2 final model: IIV CL 27.9% -> log(0.279^2 + 1) = 0.07494; covariance -0.0135; IIV V1 18% -> log(0.18^2 + 1) = 0.03189

    # Residual error: 'additive' on log-transformed concentrations (Methods 4.2),
    # i.e. exponential/log-normal on the linear scale. Table 2 prints 0.299 without
    # a scale; it is read as a standard deviation because it reproduces the
    # interquartile ranges of Figure 4 (a variance reading, SD 0.547, overshoots
    # the published upper quartiles by 10-20 mg/L). See the vignette.
    expSd <- 0.299; label("Additive residual SD on log-transformed concentration") # Table 2 final model: additive residual error = 0.299 (RSE 8%)
  })

  model({
    # Individual PK parameters (Jin 2020 Table 2 footnote).
    cl <- exp(lcl + etalcl) * (CRCL / 91.1)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 51.1)^e_wt_vc
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the central
    # compartment. Amikacin was given as a 30-min or 1-h IV infusion; users
    # specify the infusion duration per dose in the event table.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Serum concentration: dose in mg, vc in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
