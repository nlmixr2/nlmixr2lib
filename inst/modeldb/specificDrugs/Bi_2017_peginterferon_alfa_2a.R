Bi_2017_peginterferon_alfa_2a <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for",
    "peginterferon alfa-2a in adult patients with chronic hepatitis B",
    "(Bi 2017). Creatinine clearance (Cockcroft-Gault, mL/min, not",
    "BSA-normalized) modifies clearance via a power form, and body mass",
    "index modifies central volume via a power form. Exponential IIV on",
    "CL, V, and Ka; combined proportional + additive residual error on",
    "plasma concentration."
  )
  reference <- paste(
    "Bi J, Li X, Liu J, Chen D, Li S, Hou J, Zhou Y, Zhu S, Zhao Z,",
    "Qin E, Wei Z. Population pharmacokinetics of peginterferon alfa-2a",
    "in patients with chronic hepatitis B.",
    "Sci Rep. 2017;7. doi:10.1038/s41598-017-08205-5"
  )
  vignette <- "Bi_2017_peginterferon_alfa_2a"
  units <- list(time = "hour", dosing = "ng", concentration = "ng/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance in raw mL/min (NOT BSA-normalized). Source column CCR in Bi 2017.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used with a power form on clearance ((CRCL / 91.39)^e_crcl_cl); reference 91.39 mL/min is the cohort mean (Bi 2017 Table 1). Source paper writes the column as CCR; canonical column is CRCL with the raw-Cockcroft-Gault flavour (cf. Delattre 2010 amikacin entry).",
      source_name        = "CCR"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used with a power form on central volume ((BMI / 23.41)^e_bmi_vc); reference 23.41 kg/m^2 is the cohort mean (Bi 2017 Table 1).",
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 178L,
    n_studies      = 1L,
    n_observations = 208L,
    age_range      = "15-75 years (median 50.5; mean 48.40 +/- 12.91; Bi 2017 Table 1)",
    weight_range   = "42.5-100 kg (median 64; mean 65.52 +/- 11.74; Bi 2017 Table 1)",
    bmi_range      = "15.43-33.80 kg/m^2 (median 23.33; mean 23.41 +/- 3.48; Bi 2017 Table 1)",
    height_range   = "145-191 cm (median 168; mean 166.89 +/- 7.95; Bi 2017 Table 1)",
    crcl_range     = "44.80-166.87 mL/min (median 91.66; mean 91.39 +/- 24.33; Bi 2017 Table 1)",
    sex_female_pct = 44.38,
    race_ethnicity = "Not explicitly tabulated; trial conducted at 302 Military Hospital, Beijing, China, so the cohort is Chinese.",
    disease_state  = "Chronic hepatitis B (HBsAg+, HBeAg+, Anti-HBeAg-, HBV DNA >= 1e5 copies/mL, 2x ULN <= ALT <= 10x ULN, serum total bilirubin <= 2x ULN). Disease grade either hepatitis (APRI <= 2) or compensated cirrhosis (APRI > 2).",
    dose_range     = "Subcutaneous peginterferon alfa-2a once weekly; doses in source data 50000-180000 ng (50-180 ug). Standard regimen per Bi 2017 Introduction is 180 ug SC weekly for 48 weeks.",
    regions        = "China (single-center, 302 Military Hospital, Beijing).",
    trial          = "ChiCTR-RO-13004320 (Chinese Clinical Trial Registry); enrolment October 2013 - June 2016",
    notes          = "Sparse-sampling design: 208 observations from 178 patients (1-4 per patient). Blood samples randomised within 0-48 h, 48-96 h, and >96 h post-dose (target Tmax ~72 h). Serum peginterferon alfa-2a was assayed with the Pestka Biomedical Laboratories Human IFN-alpha Multi-Subtype ELISA Kit (product #41105), LLOQ 15 pg/mL. The original final-model parameter set is reported in Bi 2017 Table 2."
  )

  ini({
    # Structural parameters (Bi 2017 Table 2 final model)
    lcl <- log(0.094); label("Typical clearance at CCR 91.39 mL/min (L/h)")             # Bi 2017 Table 2: CL = 0.094 L/h
    lvc <- log(15.6);  label("Typical central volume at BMI 23.41 kg/m^2 (L)")          # Bi 2017 Table 2: V  = 15.6 L
    lka <- log(0.028); label("First-order SC absorption rate constant Ka (1/h)")        # Bi 2017 Table 2: Ka = 0.028 1/h

    # Covariate effects (Bi 2017 Equations 4-5, Table 2)
    e_crcl_cl <- 0.31; label("Exponent of (CRCL / 91.39 mL/min) on clearance (unitless)")    # Bi 2017 Table 2: CCR-CL = 0.31
    e_bmi_vc  <- 1.81; label("Exponent of (BMI / 23.41 kg/m^2) on central volume (unitless)") # Bi 2017 Table 2: BMI-V  = 1.81

    # IIV - exponential, log-normal (Bi 2017 Methods: "IIV is assumed to follow
    # a log-normal distribution"; CV% reported in Table 2). Internal variance
    # is omega^2 = log(CV^2 + 1):
    #   omega^2_CL = log(0.295^2 + 1) = 0.08344
    #   omega^2_V  = log(1.010^2 + 1) = 0.70321
    #   omega^2_Ka = log(0.640^2 + 1) = 0.34340
    etalcl ~ 0.08344  # Bi 2017 Table 2: IIV CL CV = 29.5% -> omega^2 = log(0.295^2 + 1)
    etalvc ~ 0.70321  # Bi 2017 Table 2: IIV V  CV = 101%  -> omega^2 = log(1.010^2 + 1)
    etalka ~ 0.34340  # Bi 2017 Table 2: IIV Ka CV = 64.0% -> omega^2 = log(0.640^2 + 1)

    # Residual error - combined proportional + additive on linear-scale Cc
    # (Bi 2017 Methods, "Combined error model"; Table 2 reports CV% and SD).
    propSd <- 0.194; label("Proportional residual error (fraction)")                    # Bi 2017 Table 2: residual proportional CV = 19.4%
    addSd  <- 0.32;  label("Additive residual error (ng/L)")                            # Bi 2017 Table 2: residual additive SD = 0.32 ng/L
  })

  model({
    # Individual PK parameters (Bi 2017 Equations 4-6)
    cl <- exp(lcl + etalcl) * (CRCL / 91.39)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (BMI  / 23.41)^e_bmi_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma peginterferon alfa-2a concentration: doses entered in ng, vc in
    # L -> Cc in ng/L (numerically equal to pg/mL). Bi 2017 Table 2 reports
    # the additive residual SD in ng/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
