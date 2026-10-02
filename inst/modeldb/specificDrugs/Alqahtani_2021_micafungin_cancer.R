Alqahtani_2021_micafungin_cancer <- function() {
  description <- "Two-compartment population PK model for IV micafungin in hospitalised adult patients with cancer (Alqahtani 2021, cancer-group parameter set). Linear elimination from the central compartment, log-normal IIV on CL, V1, Q and V2, and a combined additive + proportional residual error. The paper reports that AST, ALT and body weight influenced CL and that BMI, body weight, total bilirubin and albumin influenced V1, but prints no coefficient or functional form for any of them, so this file carries the cancer-group typical values of Table 2 without covariate effects. The companion file Alqahtani_2021_micafungin_noncancer holds the non-cancer group."
  reference <- paste(
    "Alqahtani S, Alfarhan A, Alsultan A, Alsarhani E, Alsubaie A, Asiri Y.",
    "Assessment of Micafungin Dosage Regimens in Patients with Cancer Using",
    "Pharmacokinetic/Pharmacodynamic Modeling and Monte Carlo Simulation.",
    "Antibiotics (Basel). 2021;10(11):1363.",
    "doi:10.3390/antibiotics10111363.",
    sep = " "
  )
  vignette <- "Alqahtani_2021_micafungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Reported as a significant covariate on CL and on V1 in both groups (Alqahtani 2021 Section 3.2 and Discussion), and varied (50, 70, 100 kg) in the Monte Carlo simulations of Table 3, but no coefficient, functional form or reference weight is printed in the article or in Supplementary File S1. Not encoded; see the vignette Assumptions and deviations.",
      source_name = "body weight"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Reported as a significant covariate on CL (Alqahtani 2021 Section 3.2); no coefficient, functional form or centring value is printed anywhere. The paper states that varying AST in the simulations did not change the PTA values (Section 3.3). Not encoded.",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Reported as a significant covariate on CL (Alqahtani 2021 Section 3.2); no coefficient, functional form or centring value is printed anywhere. The paper states that varying ALT in the simulations did not change the PTA values (Section 3.3). Not encoded.",
      source_name = "ALT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Reported as influencing V1 (Alqahtani 2021 Section 3.2); no coefficient or functional form printed. Not encoded.",
      source_name = "BMI"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Reported as influencing V1 (Alqahtani 2021 Section 3.2); no coefficient or functional form printed. Table 1 gives no unit for bilirubin; the cancer-group mean 26.5 is consistent with umol/L. Not encoded.",
      source_name = "total bilirubin"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Reported as influencing V1 (Alqahtani 2021 Section 3.2); no coefficient or functional form printed. Table 1 gives no unit for albumin; the cancer-group mean 25.6 is consistent with g/L. Not encoded.",
      source_name = "albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "mean 47.3 years (SD 12.3); adults >= 18 years",
    weight_range = "mean 63.4 kg (SD 18.2)",
    height_range = "mean 162.2 cm (SD 9.9)",
    sex_female_pct = 40,
    race_ethnicity = "Not reported (single-centre cohort in Riyadh, Saudi Arabia).",
    disease_state = "Hospitalised adults with cancer (tumour types not reported) receiving micafungin empirically or for a confirmed fungal infection. Mean serum creatinine 74.7 (printed as mmol/L; umol/L intended), CLcr 103 mL/min, albumin 25.6, AST 34.2, ALT 26.3, total bilirubin 26.5, SOFA score 7 (Table 1).",
    dose_range = "100 mg micafungin once daily as a 60-min IV infusion (all 10 cancer patients); at least two doses before sampling.",
    regions = "Saudi Arabia (King Saud University Medical City, Riyadh).",
    n_observations = 70L,
    notes = "Cancer arm of a prospective single-centre PK study of 19 patients (10 with cancer, 9 without); 7 samples per patient at 1, 2, 4, 6, 8, 12 and 24 h after the start of the infusion (133 samples in total). Demographics from Alqahtani 2021 Table 1. The companion non-cancer parameter set is Alqahtani_2021_micafungin_noncancer."
  )

  ini({
    # Structural PK parameters: Alqahtani 2021 Table 2, 'Patients with Cancer'
    # column (Monolix 4.4, SAEM).
    lcl <- log(1.2)
    label("Clearance CL (L/h)") # Table 2: CL = 1.2 L/h (RSE 11.6%)
    lvc <- log(10.7)
    label("Central volume of distribution V1 (L)") # Table 2: V1 = 10.7 L (RSE 23.6%)
    lq <- log(0.144)
    label("Intercompartmental clearance Q (L/h)") # Table 2: Q = 0.144 L/h (RSE 14%)
    lvp <- log(3.5)
    label("Peripheral volume of distribution V2 (L)") # Table 2: V2 = 3.5 L (RSE 16%)

    # IIV: log-normal (Supplementary File S1, CL_j = CL_pop x exp(eta_j)).
    # Table 2 footnote: IIV 'Expressed as coefficient of variation';
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.11 # Table 2: IIV for CL 34.1% CV; log(1 + 0.341^2)
    etalvc ~ 0.0057594 # Table 2: IIV for V1 7.6% CV; log(1 + 0.076^2)
    etalq ~ 0.098654 # Table 2: IIV for Q 32.2% CV; log(1 + 0.322^2)
    etalvp ~ 0.12701 # Table 2: IIV for V2 36.8% CV; log(1 + 0.368^2)

    # Residual error: Supplementary File S1,
    # Cobs = Cpred x (1 + eps_prop) + eps_const with independent epsilons,
    # i.e. combined2. Table 2 rows 'a' (constant) and 'b' (proportional).
    addSd <- 0.21
    label("Additive residual error (mg/L)") # Table 2: residual error a = 0.21 (RSE 10.7%)
    propSd <- 0.22
    label("Proportional residual error (fraction)") # Table 2: residual error b = 0.22 (RSE 4.5%)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, V1 in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
