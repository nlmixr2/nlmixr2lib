Purwonugroho_2012_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in hospitalised Thai adults (Purwonugroho 2012). Clearance is proportional to Cockcroft-Gault creatinine clearance (CL = 0.044 * CRCL) and central volume is proportional to age (V1 = 0.542 * AGE); intercompartmental clearance and peripheral volume have no covariates. Exponential IIV on all four parameters and an additive residual error."
  reference <- "Purwonugroho TA, Chulavatnatol S, Preechagoon Y, Chindavijak B, Malathum K, Bunuparadah P. Population pharmacokinetics of vancomycin in Thai patients. ScientificWorldJournal. 2012;2012:762649. doi:10.1100/2012/762649"
  vignette <- "Purwonugroho_2012_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation; raw mL/min, not BSA-normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Source name CLCr (Cockcroft-Gault). Purwonugroho 2012 Table 1: modeling group mean 35.07 +/- 29.83 mL/min. Stored under the canonical CRCL column per the inst/references/covariate-columns.md precedent that accepts raw mL/min when the source does not BSA-normalise (Revilla 2010 vancomycin, Jeon 2014 piperacillin). Enters CL as a through-origin linear term CL (L/h) = 0.044 * CRCL (Table 2, theta1); there is no non-renal intercept.",
      source_name = "CLCr"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Purwonugroho 2012 Table 1: modeling group mean 66.62 +/- 18.38 years (adults > 18 years by the inclusion criterion). Enters the central volume as a through-origin linear term V1 (L) = 0.542 * AGE (Table 2, theta2). The authors note in the Discussion that age is probably a surrogate of illness severity rather than a direct physiological determinant of V1.",
      source_name = "Age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Results, covariate screening) but not retained; the Discussion attributes this to the narrow weight range (86 percent of patients 40-70 kg). Table 1 modeling-group mean 57.64 +/- 11.62 kg."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Significant on CL in univariate screening but dropped before forward selection because its OFV contribution was much smaller than Cockcroft-Gault CLCr."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      notes = "Significant on CL in univariate screening but excluded from forward selection because it shares the physiological basis of CLCr and contributed less."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 212L,
    n_studies = 1L,
    age_range = "> 18 years (adults)",
    age_median = "66.62 years (mean, SD 18.38)",
    weight_range = "86 percent of patients 40-70 kg",
    weight_median = "57.64 kg (mean, SD 11.62)",
    sex_female_pct = 47.17,
    race_ethnicity = "Thai",
    disease_state = "Hospitalised adults in any ward receiving vancomycin under therapeutic drug monitoring (about one third ICU patients; 34.43 percent of the modeling group on dialysis)",
    dose_range = "Vancomycin IV infusion; dosing regimens varied between and within patients (doses and infusion durations not tabulated)",
    regions = "Thailand (Ramathibodi Hospital, Mahidol University, Bangkok; January 2007 to June 2010)",
    renal_function = "Cockcroft-Gault CLCr mean 35.07 +/- 29.83 mL/min; MDRD4 eGFR mean 44.08 +/- 36.71 mL/min/1.73 m^2; serum creatinine mean 3.12 +/- 2.92 mg/dL",
    n_concentrations = 391L,
    notes = "Retrospective TDM cohort randomised into a modeling group (212 patients, 391 serum concentrations, 1-8 per patient) and a validation group (34 patients, 68 concentrations). Serum vancomycin by FPIA (Axsym; LLOQ 2.00 mg/L). Fit in NONMEM VII (ADVAN3 TRANS4). External validation: mean prediction error -1.43 mg/L (95% CI -5.82 to 2.99), RMSE 12.28 mg/L. Demographics in Table 1."
  )

  ini({
    # Purwonugroho 2012 Table 2 (final model). Both covariate relationships are
    # through-origin linear terms, so lcl and lvc are the log of the slope:
    #   CL (L/h) = theta1 * CLCr (mL/min);  V1 (L) = theta2 * Age (years).
    lcl <- log(0.044); label("Clearance per unit Cockcroft-Gault CRCL (L/h per mL/min)") # Table 2 theta1 = 0.044 (estimation error 22.34%)
    lvc <- log(0.542); label("Central volume per year of age (L/year)") # Table 2 theta2 = 0.542 (estimation error 16.00%)
    lq <- log(6.95); label("Intercompartmental clearance Q (L/h)") # Table 2 theta3 = 6.950 L/h (estimation error 45.17%)
    lvp <- log(44.2); label("Peripheral volume V2 (L)") # Table 2 theta4 = 44.200 L (estimation error 19.58%)

    # Table 2 reports IIV as percent CV under the 'proportional' interindividual
    # error model (Results). NONMEM's CV for that model is sqrt(omega^2), so the
    # variance is (CV/100)^2 and is applied here as an exponential eta.
    etalcl ~ 0.128 # Table 2 omega CL = 35.78 percent; 0.3578^2 = 0.1280
    etalvc ~ 0.0438 # Table 2 omega V1 = 20.93 percent; 0.2093^2 = 0.0438
    etalq ~ 0.156 # Table 2 omega Q = 39.50 percent; 0.3950^2 = 0.1560
    etalvp ~ 0.328 # Table 2 omega V2 = 57.27 percent; 0.5727^2 = 0.3280

    addSd <- 4.51; label("Additive residual error (mg/L)") # Table 2 sigma = 4.51 mg/L (estimation error 38.07%)
  })

  model({
    cl <- exp(lcl + etalcl) * CRCL
    vc <- exp(lvc + etalvc) * AGE
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L: Cc in mg/L. Vancomycin is given as an IV
    # infusion; set the rate or duration on the dosing record.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
