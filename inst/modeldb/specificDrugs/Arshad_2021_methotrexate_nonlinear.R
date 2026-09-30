Arshad_2021_methotrexate_nonlinear <- function() {
  description <- "Three-compartment population PK model for high-dose intravenous methotrexate in adults with haematological malignancies or solid tumours (Arshad 2021, supplementary combined model), with parallel linear and Michaelis-Menten elimination from the central compartment. The linear clearance carries power effects of baseline serum creatinine (centred at 0.74 mg/dL) and age (58 years) and a fractional decrease in women; correlated inter-individual variability on linear CL and central volume; combined additive and exponential residual error. The authors preferred the linear-elimination model (Arshad_2021_methotrexate) for covariate and dosing work because this model's estimation was unstable."
  reference <- "Arshad U, Taubert M, Seeger-Nukpezah T, Ullah S, Spindeldreier KC, Jaehde U, Hallek M, Fuhr U, Vehreschild JJ, Jakob C. Evaluation of body-surface-area adjusted dosing of high-dose methotrexate by population pharmacokinetics in a large cohort of cancer patients. BMC Cancer. 2021;21:719. doi:10.1186/s12885-021-08443-x (Additional file 1, Supplementary Table)"
  vignette <- "Arshad_2021_methotrexate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CREAT = list(
      description = "Baseline serum creatinine at the start of treatment",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the linear clearance. The Supplementary Table prints the coefficient (-0.91) but not the functional form; the power form centred at the cohort median 0.74 mg/dL of the paper's final linear model (Eq. 4) is assumed.",
      source_name = "SCr"
    ),
    AGE = list(
      description = "Age at the start of treatment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the linear clearance centred at the cohort median 58 years, the form of Eq. 4 (assumed; the Supplementary Table prints the coefficient -0.23 only).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Supplementary Table 'Sex (fractional decrease in females)' = -0.28, applied as (1 + Sex x -0.28) following the categorical form of the Methods and Eq. 4 (sex coded 0 = male, 1 = female).",
      source_name = "Sex"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area (Du Bois formula)",
      units = "m^2",
      type = "continuous",
      notes = "Retained on CL in the linear-elimination final model (Arshad_2021_methotrexate) but printed as '-' (not included) in the Supplementary Table of this combined linear + nonlinear model."
    )
  )

  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "umol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 229,
    n_studies = 1,
    n_observations = 2182,
    age_range = "19-82 years",
    age_median = "58 years",
    weight_range = "41.5-227 kg",
    weight_median = "78.4 kg",
    bsa_range = "1.34-3.42 m^2",
    bsa_median = "1.96 m^2",
    sex_female_pct = 36.2,
    disease_state = "Adults with haematological malignancies (predominantly low-aggressive non-Hodgkin lymphoma and acute lymphoblastic leukaemia) or solid tumours treated with high-dose methotrexate; neutropenic patients from the Cologne Cohort of Neutropenic Patients (CoCoNut).",
    dose_range = "High-dose methotrexate by 4 h or 24 h intravenous infusion (18 patients occasionally 12 h or 48 h, one patient 72 h); median 3 cycles per patient (range 1-9).",
    regions = "Germany (University Hospital of Cologne, 2005-2018)",
    notes = "Same cohort as Arshad_2021_methotrexate (Arshad 2021 Tables 1 and 2). The combined model improved the OFV by 70 points over the linear model but took about 60 h per run and was unstable, so the authors did not carry it forward (Results, 'PK model')."
  )

  ini({
    # Structural parameters: bootstrap MEANS of the Supplementary Table (the
    # only values printed for this model).
    lcl <- log(4.77); label("Linear clearance for a 58-year-old man with SCr 0.74 mg/dL (L/h)") # Arshad 2021 Supplementary Table LCL = 4.77 L/h (RSE 14.1%); Results: 'linear component in the combined model was 4.77 L/h'
    lvmax <- log(2.46); label("Maximum Michaelis-Menten elimination rate (umol/h)") # Arshad 2021 Supplementary Table Vmax = 2.46 umol/h (RSE 31.6%)
    lkm <- log(1.02); label("Michaelis-Menten constant (umol/L)") # Arshad 2021 Supplementary Table Km = 1.02 umol/L (RSE 31.9%)
    lvc <- log(1.12); label("Central volume of distribution V1 (L)") # Arshad 2021 Supplementary Table V1 = 1.12 L (RSE 32.5%)
    lvp <- log(3.87); label("First peripheral volume V2 (L)") # Arshad 2021 Supplementary Table V2 = 3.87 L (RSE 24.5%)
    lvp2 <- log(5.08); label("Second peripheral volume V3 (L)") # Arshad 2021 Supplementary Table V3 = 5.08 L (RSE 30.2%)
    lq <- log(0.52); label("Intercompartmental clearance Q1, central-V2 (L/h)") # Arshad 2021 Supplementary Table Q1 = 0.52 L/h (RSE 26.7%)
    lq2 <- log(0.04); label("Intercompartmental clearance Q2, central-V3 (L/h)") # Arshad 2021 Supplementary Table Q2 = 0.04 L/h (RSE 24.1%)

    # Covariate effects on the linear clearance
    e_creat_cl <- -0.91; label("Power exponent of serum creatinine (/0.74 mg/dL) on linear CL (unitless)") # Arshad 2021 Supplementary Table SCr = -0.91 (RSE 12.0%)
    e_age_cl <- -0.23; label("Power exponent of age (/58 years) on linear CL (unitless)") # Arshad 2021 Supplementary Table Age = -0.23 (RSE 37.0%)
    e_sexf_cl <- -0.28; label("Fractional change in linear CL for women (unitless)") # Arshad 2021 Supplementary Table Sex = -0.28 (RSE 21.9%)

    # Inter-individual variability (Supplementary Table omega^2). The table's
    # IOV block repeats these three rows digit for digit, so the IOV values are
    # not recoverable and IOV is not encoded; see the vignette.
    etalcl + etalvc ~ c(0.07, 0.10, 2.127) # Arshad 2021 Supplementary Table IIV LCL 0.07, COV(LCL,V1) 0.10, V1 2.127 (omega^2)

    # Residual error: combined additive and exponential; exponential arm as its
    # first-order proportional equivalent (rxode2 cannot simulate lnorm() + add()).
    addSd <- 0.1414; label("Additive residual error SD (umol/L)") # Arshad 2021 Supplementary Table additive sigma^2 = 0.02; sqrt(0.02) = 0.1414
    propSd <- 0.4690; label("Proportional residual error, first-order equivalent of the exponential arm (fraction)") # Arshad 2021 Supplementary Table exponential sigma^2 = 0.22; sqrt(0.22) = 0.4690
  })
  model({
    # 1. Individual parameters. The linear clearance is named cllin rather
    #    than cl so rxode2 never treats the cl/vc pair as a closed-form linear
    #    system that would drop the Michaelis-Menten arm.
    cllin <- exp(lcl + etalcl) *
      (CREAT / 0.74)^e_creat_cl *
      (AGE / 58)^e_age_cl *
      (1 + e_sexf_cl * SEXF)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q <- exp(lq)
    q2 <- exp(lq2)

    # 2. Micro-constants
    kel <- cllin / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODEs (amounts in umol); parallel first-order and Michaelis-Menten
    #    elimination from the central compartment
    Cc <- central / vc
    d / dt(central) <- -kel * central - vmax * Cc / (km + Cc) -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    d / dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Observation (umol/L) and residual error
    Cc ~ add(addSd) + prop(propSd)
  })
}
