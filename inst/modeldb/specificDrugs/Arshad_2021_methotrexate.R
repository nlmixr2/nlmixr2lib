Arshad_2021_methotrexate <- function() {
  description <- "Three-compartment population PK model for high-dose intravenous methotrexate (4 h or 24 h infusion) in adults with haematological malignancies or solid tumours (Arshad 2021), with linear elimination. Clearance carries power effects of baseline serum creatinine (centred at 0.74 mg/dL), age (58 years) and body surface area (1.73 m^2) and a fractional decrease in women; correlated inter-individual variability on CL and central volume and inter-occasion (per treatment cycle) variability on CL; combined additive and exponential residual error."
  reference <- "Arshad U, Taubert M, Seeger-Nukpezah T, Ullah S, Spindeldreier KC, Jaehde U, Hallek M, Fuhr U, Vehreschild JJ, Jakob C. Evaluation of body-surface-area adjusted dosing of high-dose methotrexate by population pharmacokinetics in a large cohort of cancer patients. BMC Cancer. 2021;21:719. doi:10.1186/s12885-021-08443-x"
  vignette <- "Arshad_2021_methotrexate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CREAT = list(
      description = "Baseline serum creatinine at the start of treatment",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred at the cohort median 0.74 mg/dL (Arshad 2021 Eq. 4; Table 1 median 0.74, range 0.36-1.66 mg/dL). Covariates were evaluated at their baseline (start-of-treatment) values because of missing on-treatment data (Methods).",
      source_name = "SCr"
    ),
    AGE = list(
      description = "Age at the start of treatment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred at the cohort median 58 years (Arshad 2021 Eq. 4; Table 1 median 58, range 19-82 years).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Arshad 2021 Eq. 4 codes sex as 0 for males and 1 for females and applies (1 + Sex x -0.16) to CL, i.e. a 16% lower CL in women; the source encoding already matches SEXF.",
      source_name = "Sex"
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula) at the start of treatment",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred at 1.73 m^2 ('BSA effect was centered around the typical value of 1.73 m2', Methods). Du Bois formula per the Table 1 footnote. The exponent is +0.23 (Table 3), not the -0.23 printed in Eq. 4; see the vignette 'Assumptions and deviations'.",
      source_name = "BSA"
    ),
    OCC = list(
      description = "High-dose methotrexate treatment cycle (occasion) index, 1-9",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Inter-occasion variability was 'defined as the variability between individual cycles of MTX therapy' (Arshad 2021 Methods). Patients contributed a median of 3 cycles (range 1-9, Results), so nine occasions are encoded; occasions 2-9 share the occasion-1 variance (the analogue of NONMEM $OMEGA BLOCK(1) SAME). Set OCC to the cycle number of each dosing interval; a single-cycle simulation may use OCC = 1 throughout. Records with OCC outside 1-9 carry no IOV.",
      source_name = "OCC"
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
    notes = "Retrospective therapeutic-drug-monitoring data (Arshad 2021 Tables 1 and 2). Serum creatinine median 0.74 mg/dL (0.36-1.66); 83 female / 146 male. Methotrexate measured by competitive immunoassay, LLOQ 0.009 umol/L."
  )

  ini({
    # Structural parameters. Only CL has a printed final point estimate (the
    # 4.52 L/h of Eq. 4); every other value is the bootstrap median of Table 3
    # (1000 samples), which is the only parameter table the paper prints.
    lcl <- log(4.52); label("Clearance for a 58-year-old man with SCr 0.74 mg/dL and BSA 1.73 m^2 (L/h)") # Arshad 2021 Eq. 4 typical value 4.52 L/h (Table 3 bootstrap median 4.33, 95% CI 2.95-5.92)
    lvc <- log(4.29); label("Central volume of distribution V1 (L)") # Arshad 2021 Table 3 V1 = 4.29 L (RSE 52.5%)
    lvp <- log(2.51); label("First peripheral volume V2 (L)") # Arshad 2021 Table 3 V2 = 2.51 L (RSE 61.1%)
    lvp2 <- log(2.36); label("Second peripheral volume V3 (L)") # Arshad 2021 Table 3 V3 = 2.36 L (RSE 35.5%)
    lq <- log(0.37); label("Intercompartmental clearance Q1, central-V2 (L/h)") # Arshad 2021 Table 3 Q1 = 0.37 L/h (RSE 38.2%)
    lq2 <- log(0.02); label("Intercompartmental clearance Q2, central-V3 (L/h)") # Arshad 2021 Table 3 Q2 = 0.02 L/h (RSE 51.38%)

    # Covariate effects on CL (Eq. 4 power form; values agree with Table 3)
    e_creat_cl <- -0.49; label("Power exponent of serum creatinine (/0.74 mg/dL) on CL (unitless)") # Arshad 2021 Eq. 4 and Table 3 SCr = -0.49 (95% CI -0.31 to -0.08 as printed)
    e_age_cl <- -0.18; label("Power exponent of age (/58 years) on CL (unitless)") # Arshad 2021 Eq. 4 and Table 3 Age = -0.18 (95% CI -0.30 to 0.05)
    e_sexf_cl <- -0.16; label("Fractional change in CL for women (unitless)") # Arshad 2021 Eq. 4 and Table 3 Sex = -0.16 (95% CI -0.25 to -0.07)
    e_bsa_cl <- 0.23; label("Power exponent of BSA (/1.73 m^2) on CL (unitless)") # Arshad 2021 Table 3 BSA = 0.23 (95% CI -0.33 to 0.67); Eq. 4 prints -0.23 (sign typo, see vignette)

    # Inter-individual variability (Table 3 omega^2, bootstrap medians)
    etalcl + etalvc ~ c(0.11, 0.29, 1.34) # Arshad 2021 Table 3 IIV CL 0.11, COV(CL,V1) 0.29, V1 1.34 (omega^2)

    # Inter-occasion variability on CL, one eta per treatment cycle; occasions
    # 2-9 repeat the occasion-1 variance and are therefore fixed().
    etaiov_cl_1 ~ 0.09 # Arshad 2021 Table 3 IOV CL omega^2 = 0.09 (RSE 15.7%)
    etaiov_cl_2 ~ fixed(0.09)
    etaiov_cl_3 ~ fixed(0.09)
    etaiov_cl_4 ~ fixed(0.09)
    etaiov_cl_5 ~ fixed(0.09)
    etaiov_cl_6 ~ fixed(0.09)
    etaiov_cl_7 ~ fixed(0.09)
    etaiov_cl_8 ~ fixed(0.09)
    etaiov_cl_9 ~ fixed(0.09)

    # Residual error: 'combined (additive and exponential) error model'. The
    # exponential arm is carried as its first-order proportional equivalent
    # because rxode2 cannot simulate lnorm() + add(); see the vignette.
    addSd <- 0.1414; label("Additive residual error SD (umol/L)") # Arshad 2021 Table 3 additive sigma^2 = 0.02; sqrt(0.02) = 0.1414
    propSd <- 0.5099; label("Proportional residual error, first-order equivalent of the exponential arm (fraction)") # Arshad 2021 Table 3 exponential sigma^2 = 0.26; sqrt(0.26) = 0.5099
  })
  model({
    # 1. Occasion indicators (treatment cycles 1-9)
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    occ7 <- (OCC == 7)
    occ8 <- (OCC == 8)
    occ9 <- (OCC == 9)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7 + occ8 * etaiov_cl_8 + occ9 * etaiov_cl_9

    # 2. Individual parameters. Arshad 2021 Eq. 4:
    #    CL = 4.52 (SCr/0.74)^-0.49 (Age/58)^-0.18 (BSA/1.73)^0.23 (1 + Sex x -0.16)
    cl <- exp(lcl + etalcl + iov_cl) *
      (CREAT / 0.74)^e_creat_cl *
      (AGE / 58)^e_age_cl *
      (BSA / 1.73)^e_bsa_cl *
      (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q <- exp(lq)
    q2 <- exp(lq2)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. ODEs (amounts in umol)
    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    d / dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 5. Observation (umol/L) and residual error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
