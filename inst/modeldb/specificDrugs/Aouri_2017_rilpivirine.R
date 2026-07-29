Aouri_2017_rilpivirine <- function() {
  description <- "One-compartment population PK model for oral rilpivirine in HIV-1-infected adults (Aouri 2017), with zero-order absorption from the gastrointestinal tract directly into the central compartment (duration D1 = 4 h; derived mean absorption time D1/2 = 2 h), apparent clearance CL/F = 11.7 L/h, apparent volume of distribution V/F = 401 L, combined proportional plus additive residual error (21.6% and 9.8 ng/mL), and inter-individual variability on CL/F only (33% CV). No demographic, clinical, or genetic covariates (sex, body weight, height, age, race, AST, ALT, HCV, HBV, comedications, CYP3A4*22, CYP3A5*3, CYP2C19*2, CYP2C19*17, UGT1A1*28, UGT1A4*2) were retained in the final covariate model."
  reference <- paste(
    "Aouri M, Barcelo C, Guidi M, Rotger M, Cavassini M, Hizrel C,",
    "Buclin T, Decosterd LA, Csajka C, the Swiss HIV Cohort Study. (2017).",
    "Population pharmacokinetics and pharmacogenetics analysis of",
    "rilpivirine in HIV-1-infected individuals.",
    "Antimicrob Agents Chemother 61(1):e00899-16.",
    "doi:10.1128/AAC.00899-16"
  )
  vignette <- "Aouri_2017_rilpivirine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Male",
      notes              = "Screened in the covariate analysis but not retained in the final model. Aouri 2017 Results page (after Table 1): a 13% decrease in CL/F in females compared with males was observed (95% CI 3-25%) but did not reach statistical significance (dOFV = -5.3, below the multiple-testing significance threshold of -7.88).",
      source_name        = "SEX"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Median 74 kg (range 42-112) in the study cohort (Aouri 2017 Table 1).",
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Median 175 cm (range 150-198) in the study cohort (Aouri 2017 Table 1).",
      source_name        = "HT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Median 46 years (range 22-80) in the study cohort (Aouri 2017 Table 1).",
      source_name        = "AGE"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Recoded to a dichotomous indicator using a cutoff of 1.5 x ULN in the source covariate analysis (Aouri 2017 Methods, Covariate model section).",
      source_name        = "AST"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Recoded to a dichotomous indicator using a cutoff of 1.5 x ULN in the source covariate analysis (Aouri 2017 Methods, Covariate model section).",
      source_name        = "ALT"
    ),
    HCV_POS = list(
      description        = "Chronic hepatitis C virus coinfection indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not coinfected",
      notes              = "Screened but not retained (Aouri 2017 Methods, Covariate model section).",
      source_name        = "HCV"
    ),
    HBV_POS = list(
      description        = "Hepatitis B virus coinfection indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Not coinfected",
      notes              = "Screened but not retained (Aouri 2017 Methods, Covariate model section).",
      source_name        = "HBV"
    ),
    RACE_BLACK = list(
      description        = "Black / African ethnicity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Non-Black",
      notes              = "Race screened but not retained. Aouri 2017 Table 1 cohort: 51.4% Caucasian, 11.6% African, 0.4% Asian, 1.2% Other, 35.3% Unknown.",
      source_name        = "RACE"
    ),
    CYP3A4_STAR22 = list(
      description        = "CYP3A4*22 (rs35599367) loss-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained (dOFV less than -3.05). Genotype available for 119 Caucasian subjects in the Swiss HIV Cohort Study; MAF 0.04.",
      source_name        = "CYP3A4_22"
    ),
    CYP3A5_STAR3 = list(
      description        = "CYP3A5*3 (rs776746) loss-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained. MAF 0.90 in the Caucasian subset (Aouri 2017 Genotyping section).",
      source_name        = "CYP3A5_3"
    ),
    CYP2C19_STAR2 = list(
      description        = "CYP2C19*2 (rs4244285) loss-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained. MAF 0.16 in the Caucasian subset (Aouri 2017 Genotyping section).",
      source_name        = "CYP2C19_2"
    ),
    CYP2C19_STAR17 = list(
      description        = "CYP2C19*17 (rs12248560) gain-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained. MAF 0.24 in the Caucasian subset (Aouri 2017 Genotyping section).",
      source_name        = "CYP2C19_17"
    ),
    UGT1A1_STAR28 = list(
      description        = "UGT1A1*28 (rs8175347) loss-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained. MAF 0.39 in the Caucasian subset (Aouri 2017 Genotyping section).",
      source_name        = "UGT1A1_28"
    ),
    UGT1A4_STAR2 = list(
      description        = "UGT1A4*2 (rs6755571) loss-of-function allele indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Homozygous wild-type",
      notes              = "Tested but not retained. MAF 0.04 in the Caucasian subset (Aouri 2017 Genotyping section).",
      source_name        = "UGT1A4_2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 249L,
    n_observations = 325L,
    n_studies      = 1L,
    age_range      = "22-80 years",
    age_median     = "46 years",
    weight_range   = "42-112 kg",
    weight_median  = "74 kg",
    height_range   = "150-198 cm",
    height_median  = "175 cm",
    sex_female_pct = 27.3,
    race_ethnicity = c(Caucasian = 51.4, African = 11.6, Asian = 0.4, Other = 1.2, Unknown = 35.3),
    disease_state  = "HIV-1-infected adults receiving therapeutic drug monitoring (87% antiretroviral-experienced, 13% antiretroviral-naive)",
    dose_range     = "25 mg orally once daily, fixed-dose single-tablet regimen with emtricitabine (FTC) and tenofovir disoproxil fumarate (TDF)",
    regions        = "Switzerland (Lausanne, Zurich, Basel, Geneva, Saint Gall, Bern)",
    study_period   = "April 2013 - January 2015",
    notes          = "Demographic and genetic characteristics from Aouri 2017 Table 1. Samples were drawn 1.25-36 h after the last drug intake under steady-state conditions (at least 4 weeks after rilpivirine regimen initiation). Plasma rilpivirine concentrations ranged from 12 to 255 ng/mL; the LC-MS/MS assay LLOQ was 5 ng/mL. Genotyping was performed in 119 Swiss HIV Cohort Study participants; 130 of 249 patients had unknown genotypes."
  )

  ini({
    # Structural parameters -- Aouri 2017 Table 2 final-model estimates.
    # One-compartment open model with zero-order absorption directly into
    # the central compartment (NONMEM ADVAN1-like). Apparent oral clearance
    # CL/F and apparent volume V/F are reported with bioavailability F
    # folded into the parameters (rilpivirine was administered only orally,
    # so F is not separately identifiable). Time in hours; concentrations
    # converted from internal mg/L to ng/mL inside model() to match the
    # source paper.
    lcl <- log(11.7) ; label("Apparent clearance CL/F (L/h)")                              # Aouri 2017 Table 2 final CL/F = 11.7 L/h, RSE 2.8%
    lvc <- log(401)  ; label("Apparent volume of distribution V/F (L)")                    # Aouri 2017 Table 2 final V/F = 401 L, RSE 14.1%
    ld1 <- log(4)    ; label("Zero-order absorption duration D1 (h)")                      # Aouri 2017 Table 2 final D1 = 4 h, RSE 24.7%; mean absorption time = D1/2 = 2 h (Aouri 2017 Results, "Population pharmacokinetic analysis")

    # Inter-individual variability -- Aouri 2017 Table 2 reports IIV on CL/F
    # only (33% CV, RSE ~6%). IIV on V/F or D1 did not improve the model fit
    # (dOFV = -0.0; Aouri 2017 Results, "Population pharmacokinetic analysis").
    # Convert CV to log-scale variance: omega^2 = log(CV^2 + 1)
    #   CL/F: 33% CV -> log(0.33^2 + 1) = 0.10337
    etalcl ~ 0.10337                                                                       # Aouri 2017 Table 2 IIV CL/F = 33% CV

    # Residual unexplained variability -- Aouri 2017 Table 2 reports a
    # combined proportional + additive model on the linear concentration
    # scale (Aouri 2017 Methods, "Statistical model"). Both residual-error
    # magnitudes carry 30% shrinkage.
    addSd  <- 9.8   ; label("Additive residual error (ng/mL)")                             # Aouri 2017 Table 2: additive residual error = 9.8 ng/mL, RSE 13.6%
    propSd <- 0.216 ; label("Proportional residual error (fraction)")                      # Aouri 2017 Table 2: proportional residual error = 21.6%, RSE 7%
  })

  model({
    # Individual PK parameters -- no covariates retained in the final model.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    d1 <- exp(ld1)

    # Micro-constants
    kel <- cl / vc

    # One-compartment ODE with zero-order input directly into the central
    # compartment. The dose is delivered at constant rate amt / d1 between
    # t = 0 and t = d1 (NONMEM ADVAN1 with rate = -2 or with a duration on
    # the dose record). Mean absorption time is d1 / 2.
    d/dt(central) <- -kel * central

    # Zero-order absorption: duration d1 applies to dose records that
    # target the central compartment.
    dur(central) <- d1

    # Plasma rilpivirine concentration: dose in mg, vc in L gives mg/L;
    # multiply by 1000 to report ng/mL, matching Aouri 2017 Figure 1 axes
    # and the published Cmin / Cmax / target cutoff values (50 ng/mL).
    Cc <- (central / vc) * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
