Cheng_2020_sirolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for oral sirolimus in Chinese children with refractory immune cytopenia, developed from routine steady-state therapeutic-drug-monitoring trough concentrations (Cheng 2020). Apparent clearance scales with body weight (power exponent 0.50, normalized to the 28.5 kg cohort median) and total bilirubin (power exponent -0.32, normalized to the 11.29 umol/L cohort median); apparent volume carries no covariate. The absorption rate constant was fixed at 0.7521 per hour from the literature because only trough samples were available. Inter-individual variability is exponential on CL/F and V/F and the residual error is proportional."
  reference <- paste(
    "Cheng X, Zhao Y, Gu H, Zhao L, Zang Y, Wang X, Wu R.",
    "The first study in pediatric: Population pharmacokinetics of sirolimus",
    "and its application in Chinese children with immune cytopenia.",
    "Int J Immunopathol Pharmacol. 2020;34:2058738420934936.",
    "doi:10.1177/2058738420934936",
    sep = " "
  )
  vignette <- "Cheng_2020_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Sirolimus was assayed in whole blood by fluorescence polarization
  # immunoassay (Methods, 'Assay of sirolimus': 'Sirolimus blood
  # concentrations were determined by ... FPIA').
  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F normalized to the 28.50 kg cohort median body weight, exponent 0.50 (Results, final-model equation; Table 3 'f CL-WT'). Forward inclusion dropped the OFV by 16.88 (Table 2). Studied range 7-43 kg (mean 27.03 +/- 10.87 kg, Table 1).",
      source_name = "WT"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F normalized to the 11.29 umol/L cohort median total bilirubin, exponent -0.32 (Results, final-model equation; Table 3 'f CL-TBIL'). Forward inclusion dropped the OFV by 11.81 (Table 2). Cohort mean 12.13 +/- 7.35 umol/L (Table 1); Table 4 splits dosing at 20.50 umol/L.",
      source_name = "TBIL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods, 19 covariates) but not retained. Cohort 8.16 +/- 3.60 years, range 1-15."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as the one categorical covariate but not retained. 18 male / 9 female."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    ALB = list(description = "Serum albumin", units = "g/L", type = "continuous", notes = "Screened but not retained."),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained (indirect bilirubin, IBIL, was also screened and not retained)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (source column SCr) but not retained."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    TRIG = list(
      description = "Triglycerides",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened (source column TG) but not retained."
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened (source column TCHO) but not retained."
    ),
    LDLC = list(
      description = "Low-density lipoprotein cholesterol",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened (source column LDL) but not retained."
    ),
    RBC = list(
      description = "Red blood cell count",
      units = "10^12 cells/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (source column HB) but not retained."
    ),
    HCT = list(description = "Hematocrit", units = "%", type = "continuous", notes = "Screened but not retained.")
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    n_observations = 107L,
    age_range = "1-15 years",
    age_mean = "8.16 +/- 3.60 years",
    weight_range = "7-43 kg",
    weight_mean = "27.03 +/- 10.87 kg",
    weight_median = "28.50 kg",
    sex_female_pct = 33.3,
    race_ethnicity = "Chinese (single centre, Beijing Children's Hospital).",
    disease_state = "Children with refractory acquired or congenital single- or multiple-lineage autoimmune cytopenia in whom first- and second-line treatment was ineffective.",
    dose_range = "Oral sirolimus (gelatin capsule) once daily, initial dose 1.5 mg/m^2, then adjusted by TDM to a trough target of 5-15 ng/mL.",
    regions = "China (Beijing Children's Hospital, Capital Medical University).",
    assay = "Fluorescence polarization immunoassay (TDx/FLx, Abbott) on whole blood.",
    sampling = "107 steady-state trough concentrations drawn 0.5 h before the next dose, at least 7 days after starting sirolimus (168-7368 h after the first dose).",
    target_range = "Trough concentration 5-15 ng/mL.",
    notes = "Pilot study January 2016 to December 2017 (Table 1 demographics). Nineteen covariates were screened; only body weight and total bilirubin were retained, both on CL/F. Analysed with Phoenix NLME 1.3."
  )

  ini({
    # Final-model estimates, Cheng 2020 Table 3 (Model estimate column) and
    # the final-model equation in Results, 'PK model building'.
    lcl <- log(5.63); label("Apparent clearance CL/F at 28.5 kg and total bilirubin 11.29 umol/L (L/h)") # Table 3: CL = 5.63 L/h (RSE 9.00%)
    lvc <- log(144.16); label("Apparent volume of distribution V/F (L)") # Table 3: V = 144.16 L (RSE 23.89%)

    # Methods, 'PPK analysis': 'absorption rate constant (Ka) was fixed at
    # 0.7521/h according to relevant literature'; Table 3 prints 0.75 with RSE 0.
    lka <- fixed(log(0.7521)); label("First-order absorption rate constant Ka (1/h)") # Methods: Ka = 0.7521/h (Table 3: 0.75, RSE 0)

    e_tbili_cl <- -0.32; label("Power exponent of total bilirubin on CL/F (unitless)") # Table 3: f CL-TBIL = -0.32 (RSE 31.12%)
    e_wt_cl <- 0.50; label("Power exponent of body weight on CL/F (unitless)") # Table 3: f CL-WT = 0.50 (RSE 18.28%)

    # Abstract: 'Inter-individual variabilities for CL/F and V/F were 3.53%
    # and 7.27%'. Table 3 prints no omega rows. Read as the Phoenix omega
    # diagonal (a log-scale variance) times 100, the same x100 convention
    # the abstract applies to the residual stdev (Table 3 0.22 = '22.45%').
    etalcl ~ 0.0353 # Abstract: IIV CL/F = 3.53% (omega^2 = 0.0353)
    etalvc ~ 0.0727 # Abstract: IIV V/F = 7.27% (omega^2 = 0.0727)

    propSd <- 0.2245; label("Proportional residual error (fraction)") # Abstract: 22.45%; Table 3: Sigma = 0.22 (RSE 11.42%)
  })

  model({
    # Cohort medians used as the reference values in the final-model equation
    # (Results: '11.29 is median TBIL (umol/L)'; 'The median weight of
    # children was 28.50 kg').
    ref_wt <- 28.5
    ref_tbili <- 11.29

    # CL/F = 5.63 * (TBIL/11.29)^-0.32 * (WT/28.50)^0.5 * exp(eta_CL)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (TBILI / ref_tbili)^e_tbili_cl * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L; x1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
