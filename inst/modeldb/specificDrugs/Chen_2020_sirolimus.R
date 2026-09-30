Chen_2020_sirolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral sirolimus in Chinese children with kaposiform haemangioendothelioma, developed from routine therapeutic-drug-monitoring trough concentrations (Chen 2020). Apparent clearance CL/F scales with body weight by a fixed allometric exponent of 0.75 and is 1.999-fold higher in CYP3A5*1 carriers (expressers) than in CYP3A5*3/*3 nonexpressers; apparent volume V/F scales linearly with body weight. Both are normalized to a 70 kg standard weight. The absorption rate constant is fixed at 0.485 per hour from the group's earlier paediatric sirolimus model. Inter-individual variability is on CL/F only and the residual error is proportional."
  reference <- paste(
    "Chen X, Wang DD, Xu H, Li ZP. Initial dose recommendation for sirolimus",
    "in paediatric kaposiform haemangioendothelioma patients based on",
    "population pharmacokinetics and pharmacogenomics. J Int Med Res.",
    "2020;48(8):300060520947627. doi:10.1177/0300060520947627",
    sep = " "
  )
  vignette <- "Chen_2020_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Chen 2020 Methods 'Sirolimus administration': concentrations were
  # measured with the Emit 2000 Sirolimus Assay (Siemens), a whole-blood
  # immunoassay; the paper itself does not name the matrix, so verified =
  # FALSE on the central specimen.
  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F and V/F as a power function normalized to the 70 kg standard weight (Chen 2020 Equation iii, WTstd = 70 kg) with FIXED allometric exponents 0.75 on CL/F and 1 on V/F (Equations vi-vii). Cohort mean +/- SD 8.87 +/- 4.12 kg (Table 1); the range and median were not reported. The paper's Monte Carlo dose recommendations extrapolate to 5-60 kg.",
      source_name = "WT"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = carries at least one CYP3A5*1 allele (*1/*1 or *1/*3), 0 = CYP3A5*3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = "Chen 2020 Results below Equation vii: 'For individuals who carry the CYP3A5*3/*3 allele, the value of CYP3A5 is 0, and for individuals who carry the CYP3A5*1 allele, the value of CYP3A5 is 1' -- identical orientation to the canonical, so no value transformation is needed. Enters CL/F through the categorical form of Equation v as (1 - theta * CYP3A5) with theta = -0.999 (Table 3), i.e. a 1.999-fold CL/F in expressers (Results 'Simulation': 'the sirolimus clearance rates in individuals with CYP3A5*3/*3 and CYP3A5*1 were 1:1.999'). Cohort genotype counts (Table 2): *1/*1 n = 1, *1/*3 n = 6, *3/*3 n = 7, so 7 of 14 (50%) are expressers.",
      source_name = "CYP3A5"
    )
  )

  # Covariates screened in the stepwise search (Chen 2020 Methods 'Covariate
  # model') but not retained in the final model. Documentation only.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened ('gender'), not retained. Cohort 9 male / 5 female (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 1.53 +/- 1.40 years (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 40.41 +/- 2.86 g/L (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 29.96 +/- 18.92 IU/L (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 44.65 +/- 14.42 IU/L (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 23.52 +/- 5.38 umol/L (Table 1)."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened ('urea'), not retained. Mean +/- SD 17.57 +/- 46.25 mmol/L as printed in Table 1 (the SD larger than the mean suggests an outlier or a transcription slip in the source)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 60.71 +/- 3.92 g/L (Table 1)."
    ),
    TBA = list(
      description = "Total bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 5.36 +/- 3.11 umol/L (Table 1)."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 1.89 +/- 3.02 umol/L (Table 1)."
    ),
    TBIL = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 5.35 +/- 3.76 umol/L (Table 1)."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 33.72 +/- 2.77% (Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. Mean +/- SD 110.65 +/- 9.63 g/L (Table 1)."
    ),
    CONMED_PHENOBARBITAL = list(
      description = "Concomitant phenobarbital",
      units = "(binary)",
      type = "binary",
      notes = "Screened ('phenobarbitone'), not retained. 1 of 14 patients (Table 1)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 2 of 14 patients (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "not reported; mean +/- SD 1.53 +/- 1.40 years",
    weight_range = "not reported; mean +/- SD 8.87 +/- 4.12 kg",
    sex_female_pct = 35.7,
    race_ethnicity = "Chinese (single centre, Shanghai).",
    disease_state = "Paediatric kaposiform haemangioendothelioma treated with oral sirolimus.",
    dose_range = "0.16-1.5 mg/day oral sirolimus, adjusted by efficacy, adverse effects and TDM.",
    regions = "China (Children's Hospital of Fudan University, Shanghai), March 2016 - July 2019.",
    assay = "Emit 2000 Sirolimus Assay (Siemens), linear range 3.5-30 ng/mL.",
    sampling = "Retrospective real-world TDM; every concentration is a pre-dose trough. The number of observations is not reported.",
    target_range = "Trough 10-15 ng/mL used for the Monte Carlo initial-dose simulations.",
    notes = "Retrospective single-centre cohort; some patients' clinical data overlap the group's earlier paediatric sirolimus popPK study (Wang 2019, Oncol Lett 18:2412-2419). Pharmacogenomic panel (PGxOne, Table 2): ABCB1 rs1045642, ABCC4 rs1751034, ABCC8 rs757110, CYP2C19, CYP3A4, CYP3A5, UGT1A1 and UGT1A8 rs1042597 were screened; only CYP3A5 was retained."
  )

  ini({
    # Final-model point estimates, Chen 2020 Table 3, at the 70 kg
    # standard weight of Equation iii and in a CYP3A5*3/*3 subject.
    lcl <- log(7.55); label("Apparent clearance CL/F at 70 kg, CYP3A5*3/*3 (L/h)")        # Table 3: CL/F = 7.55 L/h (SE 15.2%); Equation vi
    lvc <- log(1840); label("Apparent volume of distribution V/F at 70 kg (L)")            # Table 3: V/F = 1840 L (SE 12.7%); Equation vii

    # Ka was not estimated: Methods 'Population pharmacokinetics model'
    # states Ka 'was fixed at 0.485/hour' (citing the earlier Wang 2019
    # paediatric sirolimus model); Table 3 lists '0.485 (fixed)'.
    lka <- fixed(log(0.485)); label("Absorption rate constant Ka (1/h)")                  # Table 3: Ka = 0.485 1/h (fixed)

    # Allometric exponents are theory-fixed (Methods 'Covariate model':
    # 'the index is the allometric coefficient 0.75 for the CL/F and 1 for
    # the V/F'); no SE is reported for either.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Methods Equation iii; Equation vi
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on V/F (unitless)")   # Methods Equation iii; Equation vii

    # CYP3A5 effect in the categorical form of Equation v,
    # CL = TV(CL) * (1 - theta * CYP3A5), with theta = -0.999 as printed, so
    # expressers have 1 + 0.999 = 1.999 times the clearance of
    # nonexpressers (Results 'Simulation': ratio 1:1.999).
    e_cyp3a5_cl <- -0.999; label("Coefficient of CYP3A5 expresser status on CL/F in (1 - theta*CYP3A5_EXPR) (unitless)")  # Table 3: theta CYP3A5 = -0.999 (SE 40.3%); Equation vi

    # IIV, exponential (Equation i), on CL/F only. Table 3 row 'omega CL/F'
    # = 0.348. Methods define the variance as omega^2, so the printed omega
    # is a standard deviation; this is independently confirmed by the
    # paper's own Monte Carlo (Figures 4 and 5), whose 95% bands and
    # target-attainment probabilities are reproduced with SD 0.348 and not
    # with variance 0.348 (see vignette). Variance = 0.348^2 = 0.121104.
    etalcl ~ 0.121104 # Table 3: omega CL/F = 0.348 (SE 17.0%), an SD; squared

    # Residual error, Equation ii: C = (1 + eps1) * Y, i.e. proportional.
    # Table 3 row 'sigma 1' = 0.390, printed with the same notation as the
    # omega row and so read as an SD (39% proportional).
    propSd <- 0.390; label("Proportional residual error (fraction)")  # Table 3: sigma 1 = 0.390 (SE 7.3%), proportional error
  })

  model({
    # Standard weight of Equation iii.
    ref_wt <- 70

    # Chen 2020 Equations vi-vii:
    #   CL/F = 7.55 * (WT/70)^0.75 * (1 - (-0.999) * CYP3A5)
    #   V/F  = 1840 * (WT/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl * (1 - e_cyp3a5_cl * CYP3A5_EXPR)
    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition with first-order absorption; F is not
    # identifiable from oral-only data and is absorbed into CL/F and V/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Doses in mg and volume in L give mg/L; x1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
