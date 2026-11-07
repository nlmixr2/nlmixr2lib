FernandezAlarcon_2024_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model for once-daily extended-release",
    "tacrolimus (Advagraf) in de novo adult kidney transplant recipients during",
    "the first 15 days post-transplant (Fernandez-Alarcon 2024), with",
    "first-order absorption (ka fixed at 2 1/h) after a lag time and linear",
    "elimination. Apparent clearance CL/F depends on recipient CYP3A5",
    "expresser status (non-expressers 0.666-fold) and decreases linearly with",
    "age above 60 years; CL/F, CLD/F and Vc/F are allometrically scaled to 70",
    "kg body weight. Disposition parameters refer to whole-blood",
    "concentrations standardized to a hematocrit of 45%; the predicted",
    "whole-blood concentration is the standardized concentration times",
    "HCT/45. Between-subject variability on CL/F and Vc/F and a proportional",
    "residual error.",
    sep = " "
  )
  reference <- paste(
    "Fernandez-Alarcon B, Nolberger O, Vidal-Alabro A, Rigo-Bonnin R,",
    "Grinyo JM, Melilli E, Montero N, Manonelles A, Coloma A, Fava A,",
    "Codina S, Cruzado JM, Colom H, Lloberas N.",
    "Guiding the starting dose of the once-daily formulation of tacrolimus in",
    "'de novo' adult renal transplant patients: a population approach.",
    "Front Pharmacol. 2024;15:1456565. doi:10.3389/fphar.2024.1456565.",
    sep = " "
  )
  vignette <- "FernandezAlarcon_2024_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling to 70 kg with fixed exponents: 0.75 on CL/F and",
        "CLD/F, 1 on Vc/F; Vp/F is not scaled (Table 4 footnote). Cohort",
        "median 70.5 kg (range 44.0-109.0; Table 1).",
        sep = " "
      ),
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F only for patients older than 60 years, through the",
        "paper's agecat indicator (agecat = 1 for age > 60 years, 0",
        "otherwise; Equation 3): CL/F loses 0.0562 L/h per year of AGE (the",
        "whole age, not the years above 60) times (WT/70)^0.75. Cohort median",
        "54 years (range 19-82; Table 1).",
        sep = " "
      ),
      source_name = "AGE"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (rs776746)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (CYP3A5*1 carrier, *1/*1 or *1/*3)",
      notes = paste(
        "1 = CYP3A5*1 carrier (*1/*1 or *1/*3); 0 = CYP3A5*3/*3 non-carrier.",
        "The paper's reference group is the carrier group: theta1 is the",
        "carrier CL/F and theta2 = 0.666 multiplies it for non-carriers",
        "(Table 4 footnote), so the model applies 0.666^(1 - CYP3A5_EXPR).",
        "Cohort: *1/*1 3 (2.2%), *1/*3 17 (12.3%), *3/*3 118 (85.5%);",
        "Table 1.",
        sep = " "
      ),
      source_name = "CYP3A5 genotype"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. The model's disposition parameters describe",
        "whole-blood concentrations standardized to a hematocrit of 45%",
        "(Methods section 2.6.2, following Storset 2014); the observed",
        "whole-blood concentration is the standardized prediction times",
        "HCT/45, with HCT in percent. Median 33.2% overall (range",
        "20.9-46.6); 29.7% on day 5, 31.3% on day 10, 32.9% on day 15",
        "(Table 2).",
        sep = " "
      ),
      source_name = "hematocrit"
    )
  )

  # Covariates screened by Fernandez-Alarcon 2024 but not retained in the
  # final model. Documented for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on CL/F alone and inside the CYP3A4/CYP3A5 cluster phenotype;",
        "neither improved the model over CYP3A5 alone (Results section 3.2).",
        "Cohort: *1/*1 126 (91.4%), *1/*22 12 (8.6%); Table 1.",
        sep = " "
      )
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 c.3435C>T (rs1045642) variant carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on CL/F (Methods section 2.6.2); not retained. Cohort:",
        "*C/*C 27.5%, *C/*T 52.2%, *T/*T 20.3%; Table 1.",
        sep = " "
      )
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested for allometric scaling of all disposition parameters; not",
        "significant and under-predicted peak concentrations (Results",
        "section 3.2). Cohort median 51.9 kg (range 30.5-73.9); Table 1.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a categorical covariate (Methods section 2.6.2); not retained. Cohort 48 of 138 female; Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested as a power covariate (Methods section 2.6.2); not retained. Median 128 umol/L (range 57-768); Table 2."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 138,
    n_studies = 1,
    n_observations = 1067,
    age_range = "54 years median (range 19-82)",
    weight_range = "70.5 kg median (range 44.0-109.0)",
    sex_female_pct = 34.8,
    disease_state = paste(
      "de novo adult kidney transplant recipients (living donor 61%,",
      "deceased donor 39%) on once-daily extended-release tacrolimus",
      "(Advagraf) started 24 h after transplantation at 0.1 mg/kg/day and",
      "titrated to a 5-10 ng/mL trough, with mycophenolic acid and",
      "corticosteroids",
      sep = " "
    ),
    dose_range = "6 mg/day median (range 0.5-23); 0.091 mg/kg/day median (0.0068-0.211)",
    regions = "Spain (single centre, Hospital Universitari de Bellvitge, Barcelona)",
    hematocrit = "33.2% median (range 20.9-46.6), rising from 29.7% on day 5 to 32.9% on day 15",
    renal_function = "eGFR (CKD-EPI) 53 mL/min median (range 5-111)",
    genotypes = paste(
      "CYP3A5 *1/*1 3 (2.2%), *1/*3 17 (12.3%), *3/*3 118 (85.5%);",
      "CYP3A4 *1/*1 126 (91.4%), *1/*22 12 (8.6%);",
      "ABCB1 *C/*C 38 (27.5%), *C/*T 72 (52.2%), *T/*T 28 (20.3%)",
      sep = " "
    ),
    notes = paste(
      "Fernandez-Alarcon 2024 Table 1 and Methods section 2.2. 29 patients",
      "had a full 24-h profile on day 5 post-transplant (0, 0.5, 1, 2, 4, 8,",
      "12, 20, 24 h) plus troughs on days 10 and 15; the other 109 had",
      "troughs on days 5, 10 and 15. Whole-blood tacrolimus by LC-MS/MS",
      "(LLOQ 1.0 ng/mL).",
      sep = " "
    )
  )

  ini({
    # Apparent clearance: CL/F = (theta1 * theta2 - theta3 * AGE * agecat) *
    # (BW/70)^0.75 (Table 4), theta2 applying to CYP3A5*1 non-carriers.
    lcl <- log(26.5); label("Apparent clearance CL/F for a 70-kg CYP3A5*1 carrier aged 60 years or younger (L/h)") # Table 4 theta1 = 26.5 L/h/70 kg (RSE 7.8%)
    e_cyp3a5_nonexpr_cl <- 0.666; label("Multiplier on CL/F for a CYP3A5*1 non-carrier (*3/*3) relative to a carrier (unitless)") # Table 4 theta2 = 0.666 (RSE 8.1%)
    e_age_cl <- 0.0562; label("Linear decrease of CL/F per year of age, applied only above 60 years (L/h/y per 70 kg)") # Table 4 theta3 = 0.0562 L/h/70 kg/y (RSE 26.0%)

    lvc <- log(327); label("Apparent central volume Vc/F for 70 kg (L)") # Table 4 Vc/F = 327 L/70 kg (RSE 14.0%)
    lvp <- log(298); label("Apparent peripheral volume Vp/F (L)") # Table 4 Vp/F = 298 L (RSE 17.4%)
    lq <- log(51.9); label("Apparent intercompartmental clearance CLD/F for 70 kg (L/h)") # Table 4 CLD/F = 51.9 L/h/70 kg (RSE 27.0%)

    lka <- fixed(log(2)); label("First-order absorption rate constant ka (1/h)") # Table 4 Ka = 2.0 h-1 FIX
    ltlag <- log(0.341); label("Absorption lag time (h)") # Table 4 Lag-time = 0.341 h (RSE 9.4%)

    # Table 4 reports the between-patient variances as CV%; converted with
    # omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.321^2) # Table 4 omega2 CL = 32.1% (RSE 15.1%)
    etalvc ~ log(1 + 0.482^2) # Table 4 omega2 Vc = 48.2% (RSE 75.4%)

    propSd <- 0.282; label("Proportional residual error (fraction)") # Table 4 Residual variability = 28.2% (RSE 3.8%)
  })

  model({
    # agecat = 1 for patients older than 60 years (Equation 3)
    agecat <- AGE > 60

    cl_typ <- (exp(lcl) * e_cyp3a5_nonexpr_cl^(1 - CYP3A5_EXPR) -
      e_age_cl * AGE * agecat) * (WT / 70)^0.75
    cl <- cl_typ * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 70)
    vp <- exp(lvp)
    q <- exp(lq) * (WT / 70)^0.75
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Disposition parameters were estimated on whole-blood concentrations
    # standardized to a hematocrit of 45%; the observable whole-blood
    # concentration scales linearly with hematocrit. mg / L * 1000 = ng/mL.
    Cstd <- 1000 * central / vc
    Cc <- Cstd * HCT / 45
    Cc ~ prop(propSd)
  })
}
