Chen_2019_risperidone_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog). Three-compartment intravenous disposition",
    "model for risperidone, fitted by the authors to the mean plasma",
    "profile after a 1 mg intravenous bolus in four beagle dogs (DAS 2.0)",
    "and used unaltered as the disposition layer of their GastroPlus",
    "ACAT/OCCAT absorption model for a risperidone orodispersible film.",
    "Disposition-only: the oral / supralingual / sublingual absorption was",
    "simulated in GastroPlus and is not reproduced here.",
    "Typical-value model (no IIV or residual error reported)."
  )
  reference <- paste(
    "Chen F, Liu H, Wang B, Yang L, Cai W, Jiao Z, Yang Z, Chen Y, Quan Y,",
    "Xiang X, Wang H. (2020). Physiologically Based Pharmacokinetic Modeling",
    "to Understand the Absorption of Risperidone Orodispersible Film.",
    "Front Pharmacol 10:1692. doi:10.3389/fphar.2019.01692."
  )
  vignette <- "Chen_2019_risperidone_dog"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Volumes and clearance are reported per kg (Table 2) and scale",
        "linearly with body weight, as GastroPlus applies them. The fitted",
        "cohort mean was 9.0425 kg (Table 2, 'Body weight')."
      ),
      source_name = "Body weight"
    )
  )

  compartmentData <- list(
    central = list(analyte = "risperidone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "risperidone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "risperidone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 4L,
    n_studies = 1L,
    weight_range = "9.04 +/- 1.88 kg (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy beagle dogs, fasted 12 h before dosing",
    dose_range = paste(
      "1 mg/body single dose in a four-period crossover (1-week washout):",
      "i.v. solution (0.2 mg/mL in 0.05% tartaric acid), i.g. dissolved",
      "ODF, supralingual ODF and sublingual ODF. Only the i.v. period",
      "informs this disposition model."
    ),
    regions = "China (Shanghai)",
    notes = paste(
      "Plasma sampled at 0.167, 0.333, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8,",
      "12, 24 and 32 h; HPLC-MS/MS, LLOQ 0.2 ng/mL. The three-compartment",
      "parameters were fitted to the i.v. profile in DAS 2.0 and entered in",
      "GastroPlus 9.7 without further alteration (Methods, Model",
      "Development)."
    )
  )

  ini({
    lcl <- log(0.5903); label("Clearance per kg body weight (L/h/kg)") # Table 2, 'CL (L/h/kg)' = 0.5903
    lvc <- log(0.3139); label("Central volume per kg body weight (L/kg)") # Table 2, 'Vc (L/kg)' = 0.3139
    lk12 <- log(16.352); label("Central to peripheral1 rate constant (1/h)") # Table 2, 'K12 (h-1)' = 16.352
    lk21 <- log(9.007); label("Peripheral1 to central rate constant (1/h)") # Table 2, 'K21 (h-1)' = 9.007
    lk13 <- log(0.4625); label("Central to peripheral2 rate constant (1/h)") # Table 2, 'K13 (h-1)' = 0.4625
    lk31 <- log(0.4103); label("Peripheral2 to central rate constant (1/h)") # Table 2, 'K31 (h-1)' = 0.4103

    propSd <- fixed(0); label("Proportional residual error (fraction); not reported") # no residual error reported; typical-value model
  })

  model({
    # Table 2 per-kg values scale linearly with body weight (GastroPlus
    # convention; the fit cohort weighed 9.0425 kg).
    cl <- exp(lcl) * WT
    vc <- exp(lvc) * WT
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    kel <- cl / vc

    # Table 2 also prints V2 = 0.56988 and V3 = 0.35383 L/kg; these are
    # implied by vc * k12 / k21 and vc * k13 / k31 and are not separate
    # parameters.
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL (= ug/L, the
    # unit of Table 1).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
