Li_2021_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Chinese infants younger than one year with septicemia (Li 2021). CL scales with body weight (reference 70 kg, estimated exponent 1.06), as a power of serum creatinine (reference 20 umol/L, exponent -0.315) and is multiplied by 1.46 with concomitant ceftriaxone. V scales linearly with body weight (reference 70 kg). IIV is on CL only; residual error is exponential (log-normal)."
  reference <- "Li Z, Li H, Wang C, Jiao Z, Xu F, Sun H. Establishment of a population pharmacokinetics model of vancomycin in 94 infants with septicemia and its application in individualized therapy. BMC Pharmacol Toxicol. 2021;22(1):26. doi:10.1186/s40360-021-00489-8"
  vignette <- "Li_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 1: mean 4.686 kg (SD 2.57), median 4 kg (range 1.4-18). Normalised to 70 kg in both the CL and V equations of the final model (Li 2021 Results 'Final model').",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 1: mean 19.91 umol/L (SD 7.45), median 18.25 (range 5.5-50). Enters CL as (SCR/20)^-0.315, so higher SCR lowers CL (Li 2021 Results 'Final model').",
      source_name = "SCR"
    ),
    CONMED_CEFTRIAXONE = list(
      description = "Concomitant ceftriaxone indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Li 2021 Results 'Final model': 'DC = 1 when incorporated with ceftriaxone, then DC = 0'. Enters CL as 1.46^DC (46% higher CL with ceftriaxone). The paper does not report how many infants received ceftriaxone.",
      source_name = "DC"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Age screened on CL and not retained (Li 2021 Methods 'PopPK modeling'); highly correlated with WT and height (Fig. 1). Table 1 reports it in days (median 88.5, range 1-345)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Li 2021 Methods 'PopPK modeling')."
    ),
    CRCL = list(
      description = "Creatinine clearance (Schwartz)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Tested as an alternative to SCR on CL; gave a similar OFV and SCR was preferred as easier to obtain (Li 2021 Results 'Covariate Models')."
    ),
    CONMED_MEROPENEM = list(
      description = "Concomitant meropenem indicator",
      units = "(binary)",
      type = "binary",
      notes = "Included in forward selection and removed at backward elimination (dOFV < 6.635; effect on CL < 20%) (Li 2021 Results 'Covariate Models')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 1L,
    age_range = "Postnatal age 1-345 days (infants younger than one year); gestational age at birth 25.7-41.4 weeks; corrected gestational age 31.14-83.07 weeks",
    age_median = "88.5 days (Table 1; the mean is 67.14 days)",
    weight_range = "1.4-18 kg",
    weight_median = "4 kg",
    sex_female_pct = 38.3,
    race_ethnicity = "Chinese (single centre, Shanghai Children's Hospital, Shanghai Jiao Tong University)",
    disease_state = "Infants with septicemia treated with intravenous vancomycin",
    dose_range = "Daily vancomycin dose 20-200 mg/day (median 60 mg/day)",
    regions = "China (Shanghai)",
    renal_function = "Serum creatinine median 18.25 umol/L (range 5.5-50); creatinine clearance median 120 mL/min/1.73 m^2 (range 21.46-280); BUN median 2.75 mmol/L.",
    n_concentrations = 205L,
    notes = "Li 2021 Table 1; enrolment January 2009 to December 2015; 58 male and 36 female infants. Routine therapeutic-drug-monitoring trough and peak samples only, measured by fluorescence polarization immunoassay (AxSYM, LLOQ 2 mg/L). NONMEM 7.2, FOCE with interaction. V IIV was not estimated because of its very large RSE with the sparse design (Table 3 model 2 and Results)."
  )

  ini({
    # Structural parameters (Li 2021 Table 2 'Parameter Value of the Final
    # Model'; reference subject WT = 70 kg, SCR = 20 umol/L, no ceftriaxone).
    lcl <- log(10.3); label("Typical CL at WT = 70 kg, SCR = 20 umol/L, no ceftriaxone (L/h)") # Li 2021 Table 2: CL = 10.3 (RSE 29.6%)
    lvc <- log(50.6); label("Typical V at WT = 70 kg (L)") # Li 2021 Table 2: V = 50.6 (RSE 7.5%)

    # Covariate effects on CL (Li 2021 Results 'Final model':
    #   CL = 10.3 * (WT/70)^1.06 * (SCR/20)^-0.315 * 1.46^DC )
    e_wt_cl <- 1.06; label("Power exponent of WT/70 on CL (unitless)") # Li 2021 Table 2: theta1 = 1.06 (RSE 9.4%)
    e_creat_cl <- -0.315; label("Power exponent of CREAT/20 on CL (unitless)") # Li 2021 Table 2: theta2 = -0.315 (RSE 20.7%)
    e_conmed_ceftriaxone_cl <- 1.46; label("Multiplicative ceftriaxone effect on CL (applied as ratio^CONMED_CEFTRIAXONE)") # Li 2021 Table 2: theta3 = 1.46 (RSE 16.7%)

    # IIV on CL only. Table 2 eta1 = 0.145 (RSE 27.7%) is read as the NONMEM
    # omega^2 (variance, CV about 39%); the paper's Figure 3 Monte Carlo
    # percentile bands are reproduced only on the variance reading (see the
    # vignette).
    etalcl ~ 0.145

    # Residual error. Table 2 labels eps1 = 0.194 (RSE 16.1%) 'Proportional
    # within-subject variability'; Methods 'PopPK modeling' states the
    # residual 'was better described by the exponential model'. It is read as
    # the NONMEM sigma^2 of Y = F * exp(eps) (FOCE-I estimates this as
    # proportional), so expSd = sqrt(0.194). Only this log-normal form
    # reproduces the Figure 3 Monte Carlo percentiles (see the vignette).
    expSd <- 0.4405; label("Exponential residual error (log-scale SD)") # Li 2021 Table 2: eps1 = 0.194 variance -> SD sqrt(0.194) = 0.4405
  })
  model({
    # Individual PK parameters (Li 2021 Results 'Final model').
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (CREAT / 20)^e_creat_cl *
      e_conmed_ceftriaxone_cl^CONMED_CEFTRIAXONE
    # V scales linearly with WT (no exponent printed in the final-model equation).
    vc <- exp(lvc) * (WT / 70)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
