Kirwan_2021_vancomycin <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for vancomycin given by intermittent infusion",
    "to critically ill adults receiving continuous venovenous haemodiafiltration (CVVHDF) in a",
    "single Irish tertiary ICU. Fitted non-parametrically with the NPAG algorithm in Pmetrics to",
    "peak and trough therapeutic-drug-monitoring levels only. Clearance (total, i.e. residual",
    "native plus CVVHDF) and volume of distribution carry inter-individual variability; no",
    "covariate was retained -- the base structural model is the final model, and it is the model",
    "the authors used for their probability-of-target-attainment dosing simulations. Residual",
    "unexplained variability is carried as fixed(0) because the Pmetrics assay-error polynomial",
    "and the fitted gamma/lambda term were never published.",
    sep = " "
  )
  reference <- paste(
    "Kirwan M, Munshi R, O'Keeffe H, Judge C, Coyle M, Deasy E, Kelly YP, Lavin PJ, Donnelly M,",
    "D'Arcy DM. Exploring population pharmacokinetic models in patients treated with vancomycin",
    "during continuous venovenous haemodiafiltration (CVVHDF). Crit Care. 2021;25(1):443.",
    "doi:10.1186/s13054-021-03863-4. PMCID: PMC8691013.",
    sep = " "
  )
  vignette <- "Kirwan_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened on CL and on V (Methods, 'Continuous covariates'); weight on V did not give a",
        "statistically significant improvement over the base model (Results) and no parameter in",
        "this model is scaled by body size. Cohort mean 81.8 kg (SD 24), median 79.9 kg (Table 1).",
        sep = " "
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL (Methods, 'Continuous covariates'); not retained in the final (base) model."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened on CL (Methods, 'Continuous covariates'); not retained in the final (base) model."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened on CL (Methods, 'Continuous covariates', 'white cell count'); not retained in the final (base) model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_dosing_intervals = 106L,
    age_range = "mean 65.5 (SD 12.3) years; median 67 (IQR 56.8-75.3)",
    weight_range = "mean 81.8 (SD 24) kg; median 79.9 (IQR 61.3-91.5)",
    bmi = "mean 28.6 (SD 8.2) kg/m^2; median 26.8 (IQR 22.8-31)",
    sex_female_pct = 25,
    race_ethnicity = "Not reported (single-centre Irish cohort)",
    disease_state = paste(
      "Critically ill adults in a tertiary ICU treated concurrently with intravenous vancomycin",
      "and CVVHDF (Prismaflex, polyarylethersulfone haemofilter; blood flow 120-220 mL/min),",
      "with CVVHDF running for the majority of each included dosing interval. APACHE II mean",
      "23 (SD 5); SOFA mean 10.5 (SD 3.2); in-hospital mortality 25%.",
      sep = " "
    ),
    renal_function = paste(
      "Acute kidney injury on continuous renal replacement therapy; urine output on study day 1",
      "mean 238.6 (SD 344.2) mL/24 h, median 77.9 mL/24 h (Table 1).",
      sep = " "
    ),
    dose_range = paste(
      "ICU policy 25 mg/kg loading dose then 15-20 mg/kg once daily guided by levels, infused at",
      "10 mg/min. Observed maintenance doses mean 1098 (SD 249.4) mg, median 1000 mg; dosing",
      "interval mean 18.5 (SD 7) h (Table 2).",
      sep = " "
    ),
    regions = "Ireland (Tallaght University Hospital, Dublin; single centre)",
    notes = paste(
      "Retrospective electronic-health-record study, 22 January 2015 to July 2016. Peak levels",
      "drawn 60 min after the end of infusion and trough levels immediately before the next dose;",
      "total serum vancomycin by the Architect iVancomycin chemiluminescent microparticle",
      "immunoassay. The Results text reports 74 peak and 96 trough levels (170 in total) while",
      "the Abstract reports 155 levels; both are transcribed as printed.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural values are the MEANS of the NPAG non-parametric parameter
    # distribution for the base model, Kirwan 2021 Table 3. The Results text
    # quotes the means ('the mean CL was 2.59 L/h +/- 0.49 with a mean V of
    # 80.98 L'), and the means (not the medians 2.70 L/h / 73.72 L) reproduce
    # the paper's own probability-of-target-attainment results -- see vignette.
    # ------------------------------------------------------------------------
    lcl <- log(2.59); label("Total vancomycin clearance on CVVHDF (L/h)")
    # Table 3, base model CL mean = 2.59 L/h (SD 0.49, CV 18.99%, median 2.70)
    lvc <- log(80.98); label("Volume of distribution (L)")
    # Table 3, base model V mean = 80.98 L (SD 16.89, CV 20.86%, median 73.72)

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # joint density, not a parametric omega. Table 3 prints a CV% (= SD/mean,
    # table footnote) per parameter, carried here as a LOG-NORMAL approximation
    # via omega^2 = log(CV^2 + 1). No CL-V covariance is reported, so the two
    # etas are independent.
    # ------------------------------------------------------------------------
    # Table 3 base model CL CV% = 18.99 -> log(0.1899^2 + 1)
    etalcl ~ 0.035427
    # Table 3 base model V CV% = 20.86 -> log(0.2086^2 + 1)
    etalvc ~ 0.042594

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported: neither the Pmetrics
    # assay-error polynomial (C0-C3) nor the fitted gamma/lambda appear in the
    # paper or in Additional file 1. Carried as fixed(0) rather than invented.
    # ------------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    # 1. Individual parameters (no covariates in the final base model).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment ODE; vancomycin is given as an intermittent
    #    intravenous infusion into central (ICU policy 10 mg/min).
    d/dt(central) <- -kel * central

    # 4. Observation: total serum vancomycin. Pmetrics weights observations by
    #    a linear assay-error polynomial, hence the combined1() form.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
