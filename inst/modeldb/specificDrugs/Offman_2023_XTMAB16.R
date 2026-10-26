Offman_2023_XTMAB16 <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "XTMAB-16, a chimeric anti-TNF-alpha monoclonal antibody, after a single",
    "2-h intravenous infusion of 2 or 4 mg/kg in 19 healthy adults of the",
    "first-in-human study XTMAB-16-101 (NCT04971395; Offman 2023",
    "Supplementary Table S4). Central and peripheral volumes share one",
    "estimated body-weight power exponent (reference 72.4 kg); clearance is",
    "multiplied by 1.274 in subjects who were ever anti-drug-antibody",
    "positive. IIV on CL, V and V2; combined additive plus proportional",
    "residual error.",
    sep = " "
  )
  reference <- paste(
    "Offman E, Singh N, Julian MW, Locke LW, Bicer S, Mitchell J,",
    "Matthews T, Anderson K, Crouser ED. Leveraging in vitro and",
    "pharmacokinetic models to support bench to bedside investigation of",
    "XTMAB-16 as a novel pulmonary sarcoidosis treatment.",
    "Front Pharmacol. 2023;14:1066454. doi:10.3389/fphar.2023.1066454.",
    "Final parameter estimates from Supplementary Table S4; structural and",
    "covariate model from the Results section 'Model development and",
    "qualification'.",
    sep = " "
  )
  vignette <- "Offman_2023_XTMAB16"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "XTMAB-16", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "XTMAB-16", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight, normalised to 72.4 kg (the study median, Offman 2023 Table 1). One estimated exponent of 0.426 scales both the central (V) and peripheral (V2) volumes (Supplementary Table S4). Table S4 carries no weight term on CL, and the CL-per-kg and AUCinf values of Supplementary Table S7 for a 70 kg subject are reproduced only with CL independent of weight, although the Discussion lists weight among the covariates on CL.",
      source_name = "Weight"
    ),
    ADA_POS = list(
      description = "Subject-level anti-drug-antibody status (1 = at least one ADA-positive sample)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "Time-fixed per subject: a subject with at least one positive ADA sample is ADA-positive (Offman 2023 Table 1 footnote). In the study all ADA-positive subjects were also neutralising-ADA-positive (Table 1). Clearance is multiplied by 1.274 for ADA-positive subjects (Supplementary Table S4). The study was single-dose, and the paper's simulations used ADA-negative subjects only.",
      source_name = "subject-level ADA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19,
    n_studies = 1,
    n_observations = 277,
    age_range = "19-45 years",
    age_median = "32 years (mean 31.8)",
    weight_range = "53.2-97.2 kg",
    weight_median = "72.4 kg (mean 72.3)",
    sex_female_pct = 68.4,
    race_ethnicity = c(`Black or African American` = 84.2, White = 10.5, `American Indian or Alaska Native` = 5.3),
    disease_state = "healthy",
    dose_range = "single 2 or 4 mg/kg IV infusion over 2 h",
    regions = "not reported",
    ada_positive_pct = 42.1,
    notes = paste(
      "Randomised, double-blind, placebo-controlled first-in-human study",
      "XTMAB-16-101 (NCT04971395) in healthy non-smoking adults aged 18-45",
      "years with BMI 18-30 kg/m^2. 25 participants were enrolled (10 at",
      "2 mg/kg, 9 at 4 mg/kg and 6 on placebo, Supplementary Table S1); the",
      "PK analysis used the 19 XTMAB-16-treated subjects with 277 serum",
      "samples, of which 226 (81.6%) were quantifiable; BLQ samples were set",
      "to missing (Offman 2023 Table 1 and Results 'PK measurements').",
      "Height 158-185 cm, BMI 20.9-29.8 kg/m^2. Fit in Phoenix NLME 8.4.",
      sep = " "
    )
  )

  ini({
    # Structural parameters at the reference subject: 72.4 kg, ADA-negative
    # (Supplementary Table S4).
    lcl <- log(0.00902) ; label("Clearance CL in ADA-negative subjects (L/h)")                  # Table S4: CL = 0.00902 L/h (RSE 6.2%, 95% CI 0.008-0.010)
    lvc <- log(2.30)    ; label("Central volume of distribution V at 72.4 kg (L)")              # Table S4: V = 2.30 x (Weight/72.4)^0.426 L (RSE 2.2%, 95% CI 2.20-2.40)
    lvp <- log(0.715)   ; label("Peripheral volume of distribution V2 at 72.4 kg (L)")          # Table S4: V2 = 0.715 x (Weight/72.4)^0.426 L (RSE 13.1%, 95% CI 0.53-0.90)
    lq  <- log(0.0119)  ; label("Inter-compartmental clearance Q (L/h)")                        # Table S4: Q = 0.0119 L/h (RSE 22.5%, 95% CI 0.007-0.017)

    # One weight exponent is printed for both volumes in Table S4.
    e_wt_vc <- 0.426 ; label("Power exponent on (WT/72.4) for V and V2 (unitless)")             # Table S4: (Weight/72.4)^0.426 on both V and V2

    e_ada_pos_cl <- 1.274 ; label("Multiplicative factor on CL for ADA-positive subjects (unitless)")  # Table S4: CL x 1.274 if subject is ADA positive

    # IIV, reported as variances (Table S4 footnote a).
    etalcl ~ 0.0246   # Table S4: IIV on CL = 0.0246 (RSE 6.7%, shrinkage 7.1%)
    etalvc ~ 0.00936  # Table S4: IIV on V = 0.00936 (RSE 2.7%, shrinkage 4.2%)
    etalvp ~ 0.203    # Table S4: IIV on V2 = 0.203 (RSE 16.8%, shrinkage 17.7%)

    # Combined additive plus proportional residual error.
    propSd <- 0.0689 ; label("Proportional residual error (fraction)")   # Table S4: proportional error = 6.89%
    addSd  <- 0.819  ; label("Additive residual error (mg/L)")           # Table S4: additive error = 819 ng/mL = 0.819 mg/L
  })

  model({
    # Individual parameters. Weight scales only the volumes; ADA_POS is the
    # subject-level ever-positive indicator.
    cl <- exp(lcl + etalcl) * e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 72.4)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 72.4)^e_wt_vc
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion into central (2 h in the study).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # central in mg and vc in L, so Cc is in mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
