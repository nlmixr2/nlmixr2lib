Yin_2021_mitotane_nogenotype <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral mitotane in adults with adrenocortical carcinoma (Yin 2021), reduced model without genotype covariates. This is the authors' alternative model for patients whose CYP2C19 / SLCO1B3 / SLCO1B1 genotype is unknown (Online Resource 1 Table S3, built into the authors' Shiny app). Apparent central volume carries a power effect of fat amount (total body weight minus lean body weight by the Boer formula); apparent clearance has no covariate. Log-normal interindividual variability on CL/F, Vc/F, Vp/F and Q/F, interoccasion variability on CL/F with one occasion per 200 days of treatment, and combined additive plus proportional residual error. The absorption rate constant is fixed. Time in days. Yin_2021_mitotane is the companion pharmacogenetic model with genotype effects on CL/F."
  reference <- paste(
    "Yin A, Ettaieb MHT, Swen JJ, van Deun L, Kerkhofs TMA,",
    "van der Straaten RJHM, Corssmit EPM, Gelderblom H, Kerstens MN,",
    "Feelders RA, Eekhoff M, Timmers HJLM, D'Avolio A, Cusato J,",
    "Guchelaar HJ, Haak HR, Moes DJAR. Population pharmacokinetic and",
    "pharmacogenetic analysis of mitotane in patients with adrenocortical",
    "carcinoma: towards individualized dosing. Clin Pharmacokinet.",
    "2021;60(1):89-102. doi:10.1007/s40262-020-00913-y.",
    "Parameter estimates from Online Resource 1 Table S3 (reduced model",
    "without genotype covariates); the ODE system, the standard-deviation",
    "scale of the printed CV% values and the Boer lean-body-weight",
    "coefficients are taken from the authors' published Shiny app script",
    "(github.com/AnyueYin/Shiny-app-script-for-model-simulation---Population-PK-and-PG-analysis-of-mitotane).",
    sep = " "
  )
  vignette <- "Yin_2021_mitotane"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mitotane", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mitotane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mitotane", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at the start of mitotane treatment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; Table 1 mean 80.0 kg (SD 15.9, range 52.5-120). Enters only through the derived fat amount (FAT = WT - LBW).",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height at the start of mitotane treatment",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; Table 1 mean 172 cm (SD 10.0, range 154-193). Used only in the Boer lean-body-weight formula that yields fat amount.",
      source_name = "HT"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the sex-specific Boer lean-body-weight formula: male LBW = 0.407*WT + 0.267*HT - 19.2, female LBW = 0.252*WT + 0.473*HT - 48.3. Used only to derive fat amount for the Vc/F effect.",
      source_name = "SEX"
    ),
    OCC = list(
      description = "Occasion index for interoccasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Every 200 days of treatment defines an occasion (Table S3 footnote a). Decomposed inside model() into occ1..occ8 selecting per-occasion IOV etas on CL. Eight occasions cover 0-1600 days. See the vignette Errata.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    disease_state = "adrenocortical carcinoma (ENSAT I-IV)",
    age_range = "22.6-76.8 years (mean 52.0)",
    weight_range = "52.5-120 kg (mean 80.0)",
    sex_female_pct = 56.3,
    dose_range = "0.5-16 g/day total daily dose",
    regions = "Netherlands (Dutch Adrenal Network Registry)",
    notes = "Same cohort as the companion pharmacogenetic model Yin_2021_mitotane; this reduced model omits the genotype covariates for use when genotype is unavailable."
  )

  ini({
    lka <- fixed(log(15.0))
    label("Absorption rate constant KA (1/day); Table S3, held at the value estimated in the pharmacogenetic model")
    lcl <- log(217)
    label("Apparent clearance CL/F (L/day); Table S3 = 217")
    lvc <- log(8450)
    label("Apparent central volume Vc/F (L); Table S3 = 8450")
    lq <- log(609)
    label("Apparent intercompartmental clearance Q/F (L/day); Table S3 = 609")
    lvp <- log(15500)
    label("Apparent peripheral volume Vp/F (L); Table S3 = 15500")

    e_fat_vc <- 1.12
    label("Power exponent of fat amount on Vc/F, (FAT/23.6)^e (unitless); Table S3 Vc_FAT = 1.12, reference FAT 23.6 kg")

    # IIV variances = (printed CV%/100)^2 on the SD scale (authors' Shiny app).
    etalcl ~ 0.439569 # Table S3 IIV CL/F 66.3% -> SD 0.663
    etalvc ~ 0.403225 # Table S3 IIV Vc/F 63.5% -> SD 0.635
    etalq ~ 1.010025 # Table S3 IIV Q/F 100.5% -> SD 1.005
    etalvp ~ 0.646416 # Table S3 IIV Vp/F 80.4% -> SD 0.804

    # IOV on CL/F, one occasion per 200 days; Table S3 IOV 31.2% -> SD 0.312.
    etaiov_cl_1 ~ 0.097344
    etaiov_cl_2 ~ fixed(0.097344)
    etaiov_cl_3 ~ fixed(0.097344)
    etaiov_cl_4 ~ fixed(0.097344)
    etaiov_cl_5 ~ fixed(0.097344)
    etaiov_cl_6 ~ fixed(0.097344)
    etaiov_cl_7 ~ fixed(0.097344)
    etaiov_cl_8 ~ fixed(0.097344)

    propSd <- 0.167
    label("Proportional residual error (Table S3 PRO CV% = 16.7, SD scale)")
    addSd <- 0.907
    label("Additive residual error (mg/L) (Table S3 ADD = 0.907)")
  })

  model({
    lbw <- SEXF * (0.252 * WT + 0.473 * HT - 48.3) +
      (1 - SEXF) * (0.407 * WT + 0.267 * HT - 19.2)
    fat <- WT - lbw

    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    occ7 <- (OCC == 7)
    occ8 <- (OCC == 8)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 +
      occ3 * etaiov_cl_3 + occ4 * etaiov_cl_4 +
      occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7 + occ8 * etaiov_cl_8

    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc) * (fat / 23.6)^e_fat_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
