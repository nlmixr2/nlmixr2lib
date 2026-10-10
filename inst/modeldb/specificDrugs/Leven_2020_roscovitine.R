Leven_2020_roscovitine <- function() {
  description <- "Two-compartment parent-metabolite population PK model for oral roscovitine (seliciclib) and its carboxylate metabolite M3 in adults with cystic fibrosis chronically infected with Pseudomonas aeruginosa (Leven 2020), with a dose-dependent saturable first-pass effect (Emax-type fraction of the dose escaping first-pass conversion) and a sum-of-inverse-Gaussian-densities absorption input (two parallel roscovitine inputs, one pre-systemic M3 input sharing the first roscovitine input's shape). The roscovitine and M3 central volumes share one apparent V/F, and roscovitine has no elimination route other than systemic conversion to M3."
  reference <- "Leven C, Schutz S, Audrezet MP, Nowak E, Meijer L, Montier T. Non-Linear Pharmacokinetics of Oral Roscovitine (Seliciclib) in Cystic Fibrosis Patients Chronically Infected with Pseudomonas aeruginosa: A Study on Population Pharmacokinetics with Monte Carlo Simulations. Pharmaceutics. 2020;12(11):1087. doi:10.3390/pharmaceutics12111087"
  vignette <- "Leven_2020_roscovitine"
  units <- list(time = "h", dosing = "umol", concentration = "nmol/L")

  covariateData <- list(
    HT = list(
      description = "Body height at baseline",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the shared apparent central volume V/F, (HT / 170)^e_ht_vc, with the reference height of 170 cm stated in the Results (Section 3.2). Observed range 158-182 cm (Table 2).",
      source_name = "height"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no proton-pump inhibitor)",
      notes = "Subject-level indicator of concomitant esomeprazole, omeprazole or pantoprazole during the study (15 of 23 patients, Table 2). Enters as ln(T1max_i) = ln(T1max) + beta * PPI (Equation 9), lengthening the time to peak of the first inverse-Gaussian input (and, through TM3max = T1max, of the pre-systemic M3 input) by exp(0.680) = 1.97-fold. Table 3 and its footnote place the effect on T1max; the Results prose says 'dT2max', which Figure 9 (the M3 peak moving later under PPI) rules out.",
      source_name = "PPI"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "roscovitine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "roscovitine", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "roscovitine", units = "umol", specimen = "plasma", verified = TRUE),
    central_m3 = list(
      analyte = "roscovitine carboxylate metabolite M3",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "22-51 years",
    age_median = "33 years",
    weight_range = "48-90 kg",
    weight_median = "58 kg",
    sex_female_pct = 43.5,
    disease_state = "Adults with cystic fibrosis carrying two CF-causing mutations including at least one F508del-CFTR mutation, FEV1 >= 40%, chronically infected with Pseudomonas aeruginosa",
    dose_range = "200, 400 or 800 mg roscovitine orally (hard gelatin capsules, 200 mg each); PK sampled over 12 h after the first dose",
    regions = "France (12 CF clinical trial centres)",
    notes = "ROSCO-CF phase IIa dose-ranging trial (NCT02649751): 23 of 34 randomised subjects received roscovitine (9 at 200 mg, 7 at 400 mg, 7 at 800 mg); 138 roscovitine (19 BLQ) and 138 M3 (9 BLQ) concentrations (Table 1). Height median 166 cm (158-182); MDRD GFR median 113 mL/min (61-202); 15 of 23 on proton-pump inhibitors (Table 2). Estimation in Monolix 2019R1 (SAEM), BLQ data handled by the M3 method."
  )

  ini({
    # Dose-dependent first-pass effect: fdose = Dose / (Dose + D50), Equation 7
    # and Appendix A 'Fr = amtDose/(amtDose + D50)'; dose in umol of roscovitine
    # (MW 354.45 g/mol, Section 2.3).
    led50 <- log(1190); label("Oral roscovitine dose at which half the dose escapes first-pass conversion to M3, D50 (umol)") # Table 3 'D50 (umol) 1190'

    # Sum-of-inverse-Gaussian absorption input, Equations 4-6 and Appendix A.
    logitfrel <- log(0.285 / (1 - 0.285)); label("Logit of the weight pi1 of the first roscovitine inverse-Gaussian input (unitless)") # Table 3 'pi1 0.285' (logit-normal, Section 2.5; no IIV estimated)
    ltpeak_invgauss1 <- log(0.676); label("Time of the peak of the first roscovitine inverse-Gaussian input density, T1max (h)") # Table 3 'T1max (h) 0.676'
    e_conmed_ppi_tpeak_invgauss1 <- 0.680; label("Log-scale shift in T1max with concomitant proton-pump inhibitor, beta (unitless)") # Table 3 'beta (PPI = 1) on T1max 0.680'
    ldtpeak_invgauss2 <- log(1.04); label("Delay of the second roscovitine inverse-Gaussian input peak after the first, dT2max (h)") # Table 3 'dT2max (h) 1.04'; T2max = T1max + dT2max (Section 3.2)
    lcv_invgauss1 <- log(0.542); label("Coefficient of variation of the first inverse-Gaussian input density, CV1 (unitless)") # Table 3 'CV1 0.542'
    lcv_invgauss2 <- log(0.354); label("Coefficient of variation of the second inverse-Gaussian input density, CV2 (unitless)") # Table 3 'CV2 0.354'

    # Disposition
    lvc <- log(62.2); label("Apparent central volume V/F shared by roscovitine and M3 at 170 cm height (L)") # Table 3 'V/F (liters) 62.2'
    e_ht_vc <- 6.47; label("Power exponent of height (HT / 170 cm) on V/F (unitless)") # Table 3 'beta height (cm) on V 6.47'; Equation 10, reference 170 cm (Section 3.2)
    lk12 <- log(1.82); label("Roscovitine central-to-peripheral rate constant k12 (1/h)") # Table 3 'k12 (h-1) 1.82'
    lk21 <- log(0.768); label("Roscovitine peripheral-to-central rate constant k21 (1/h)") # Table 3 'k21 (h-1) 0.768'
    lkmet <- log(2.07); label("Systemic conversion rate constant of roscovitine to M3, kmet (1/h)") # Table 3 'kmet (h-1) 2.07'
    lkel_m3 <- log(2.58); label("Elimination rate constant of M3, ke (1/h)") # Table 3 'ke (h-1) 2.58'

    # IIV: Table 3 reports standard deviations of log-normal random effects;
    # variances are omega^2, covariances corr * omega_a * omega_b.
    etaled50 ~ 0.474721 # Table 3 'omega D50 0.689' -> 0.689^2
    etaltpeak_invgauss1 + etalcv_invgauss2 ~ c(0.204304, -0.207438, 0.299209) # Table 3 'omega T1max 0.452', 'omega CV2 0.547', 'Corr. T1max CV2 -0.839' -> -0.839 * 0.452 * 0.547
    etaldtpeak_invgauss2 + etalcv_invgauss1 ~ c(0.418609, 0.171988, 0.181476) # Table 3 'omega dT2max 0.647', 'omega CV1 0.426', 'Corr. dT2max CV1 0.624' -> 0.624 * 0.647 * 0.426
    etalvc ~ 0.459684 # Table 3 'omega V/F 0.678' -> 0.678^2
    etalkel_m3 ~ 0.154449 # Table 3 'omega ke 0.393' -> 0.393^2

    # Residual error (Monolix proportional for roscovitine; combined
    # g = sqrt(a^2 + b^2 f^2) for M3, Section 2.5)
    propSd <- 0.297; label("Proportional residual error for roscovitine (fraction)") # Table 3 'Roscovitine b1 0.297'
    addSd_m3 <- 9.20; label("Additive residual error for M3 (nmol/L)") # Table 3 'M3 a2 (nmol L-1) 9.20'
    propSd_m3 <- 0.271; label("Proportional residual error for M3 (fraction)") # Table 3 'M3 b2 0.271'
  })

  model({
    # Individual parameters
    ed50 <- exp(led50 + etaled50)
    frel <- expit(logitfrel)
    tpeak_invgauss1 <- exp(ltpeak_invgauss1 + e_conmed_ppi_tpeak_invgauss1 * CONMED_PPI + etaltpeak_invgauss1)
    dtpeak_invgauss2 <- exp(ldtpeak_invgauss2 + etaldtpeak_invgauss2)
    cv_invgauss1 <- exp(lcv_invgauss1 + etalcv_invgauss1)
    cv_invgauss2 <- exp(lcv_invgauss2 + etalcv_invgauss2)
    vc <- exp(lvc + etalvc) * (HT / 170)^e_ht_vc
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kmet <- exp(lkmet)
    kel_m3 <- exp(lkel_m3 + etalkel_m3)

    # Mean absorption times from the peak time and CV (Equation 6; Appendix A
    # MAT_1 and MAT_2). The M3 input shares T1max and CV1 (TM3max = T1max,
    # CVM3 = CV1, Section 3.2), so its density is the first roscovitine one.
    tpeak_invgauss2 <- tpeak_invgauss1 + dtpeak_invgauss2
    mat_invgauss1 <- tpeak_invgauss1 / (sqrt(1 + 9 / 4 * cv_invgauss1^4) - 3 / 2 * cv_invgauss1^2)
    mat_invgauss2 <- tpeak_invgauss2 / (sqrt(1 + 9 / 4 * cv_invgauss2^4) - 3 / 2 * cv_invgauss2^2)

    # Dose amount (umol) and time since the dose: Monolix amtDose and t in
    # Appendix A. The whole dose lands in depot and depot empties at exactly the
    # prescribed input rate, so mass is conserved.
    dose_umol <- podo(depot)
    tdose <- tad(depot)
    fdose <- dose_umol / (dose_umol + ed50)

    # Inverse-Gaussian input densities (Equation 5; Appendix A inv_gauss_1..3),
    # zero at and before the dose time where the density's limit is 0.
    invgauss1 <- 0
    invgauss2 <- 0
    if (tdose > 0) {
      invgauss1 <- sqrt(mat_invgauss1 / (2 * pi * cv_invgauss1^2 * tdose^3)) * exp(-(tdose - mat_invgauss1)^2 / (2 * cv_invgauss1^2 * mat_invgauss1 * tdose))
      invgauss2 <- sqrt(mat_invgauss2 / (2 * pi * cv_invgauss2^2 * tdose^3)) * exp(-(tdose - mat_invgauss2)^2 / (2 * cv_invgauss2^2 * mat_invgauss2 * tdose))
    }
    input_parent <- fdose * dose_umol * (frel * invgauss1 + (1 - frel) * invgauss2)
    input_m3 <- (1 - fdose) * dose_umol * invgauss1

    d/dt(depot) <- -input_parent - input_m3
    d/dt(central) <- input_parent - kmet * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_m3) <- kmet * central + input_m3 - kel_m3 * central_m3

    # Concentrations in nmol/L (umol / L * 1000; Appendix A C_rosco, C_M3)
    Cc <- central / vc * 1000
    Cc_m3 <- central_m3 / vc * 1000

    Cc ~ prop(propSd)
    Cc_m3 ~ add(addSd_m3) + prop(propSd_m3)
  })
}
