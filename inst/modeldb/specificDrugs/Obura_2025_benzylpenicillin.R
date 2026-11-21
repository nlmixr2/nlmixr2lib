Obura_2025_benzylpenicillin <- function() {
  description <- paste(
    "Two-compartment IV population PK model of intrapartum benzylpenicillin (penicillin G) in",
    "12 labouring women at risk of early-onset neonatal Group B Streptococcus disease (Obura",
    "2025), with a third, umbilical cord compartment linked to the maternal central compartment",
    "by first-order rate constants. Maternal plasma and umbilical cord blood concentrations are",
    "both observed. Pmetrics NPAG non-parametric fit; typical values are the weighted means of",
    "the joint density (Table 2) and inter-individual variability is a log-normal approximation",
    "built from the reported CV% and the CL-Vc correlation. No covariates were retained.",
    "Concentrations are total drug; the paper converts to unbound with a fixed 60% protein",
    "binding for target attainment.",
    sep = " "
  )
  reference <- paste(
    "Obura B, Unsworth J, Jimenez-Valverde A, Waitt C, Das S, Hope W, Navaratnam K.",
    "Pharmacokinetics of Intrapartum Benzylpenicillin: Insights Into Candidate Regimens to",
    "Prevent Early Onset Neonatal Group B Streptococcus Disease.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(9):1504-1514.",
    "doi:10.1002/psp4.70072. PMID: 40629558. PMCID: PMC12439288.",
    sep = " "
  )
  vignette <- "Obura_2025_benzylpenicillin"

  # `cord` is the umbilical cord compartment of Obura 2025 Figure 2 and
  # Equations 6 and 8 (X3, volume Vu). It is an observed physiologic matrix
  # with its own output equation, not a generic `peripheral2`. Declared
  # paper-specific following Nie_2025_sufentanil.R, which carries the same
  # `cord` state, `lvcord` volume and `Ccord` output.
  paper_specific_compartments <- c("cord")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "benzylpenicillin", units = "mg", specimen = "tissue", verified = TRUE),
    cord = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Maternal age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods, Covariate Modeling); not retained. Table 1: median 30.5, range 19-37 years."
    ),
    WT = list(
      description = "Total body weight at the booking visit (mother plus unborn baby)",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened with allometric exponents fixed at 1 (volume) and 0.75 (clearance); not retained,",
        "although the weight-Vc correlation 'barely approached significance'. Table 1: median 75,",
        "range 52-132 kg."
      )
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw mL/min, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened; not retained. The CrCl distribution is not printed; Table 1 gives serum",
        "creatinine 0.56 +/- 0.12 mg/dL (mean +/- SD)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "19-37 years",
    age_median = "30.5 years",
    weight_range = "52-132 kg",
    weight_median = "75 kg",
    bmi_mean = "30.7 +/- 10.45 kg/m^2",
    sex_female_pct = 100,
    race_ethnicity = c(White = 75, Asian = 25),
    disease_state = paste(
      "Pregnant women in labour requiring intrapartum antibiotic prophylaxis against early-onset",
      "neonatal GBS disease: preterm labour (n = 5), GBS-positive antenatal swab (n = 4) or GBS",
      "bacteriuria (n = 3). Planned elective caesarean section, penicillin allergy and intrapartum",
      "pyrexia (>= 38 C) were exclusions. Gravidity median 2 (1-4); parity median 0 (0-3)."
    ),
    renal_function = "Serum creatinine 0.56 +/- 0.12 mg/dL (mean +/- SD).",
    dose_range = paste(
      "Benzylpenicillin 3 g IV loading dose at onset of labour, then 1.5 g every 4 h until",
      "delivery, each as an approximately 30-min intermittent infusion (RCOG regimen). Two women",
      "received no maintenance dose, four received 1, three received 2, two received 3 and one",
      "received 4."
    ),
    regions = "United Kingdom (Liverpool Women's Hospital), January 2022 to June 2023.",
    n_observations = paste(
      "45 maternal plasma and 10 umbilical cord blood concentrations. Maternal samples before and",
      "0.25, 0.5, 1, 2 and 4 h after the first dose; one cord sample at delivery. One cord and two",
      "plasma samples were excluded. LC-MS/MS, LLOQ 0.1 mg/L; all samples quantifiable."
    ),
    notes = paste(
      "Pmetrics 2.1.1 NPAG (R 4.3.1). Total body weight, age and Cockcroft-Gault creatinine",
      "clearance were screened and none retained, so the base model is final."
    )
  )

  ini({
    # Pmetrics NPAG final-model estimates, Obura 2025 Table 2 (columns Mean /
    # Median / 95% CI / CV%). The weighted MEANS are used as typical values:
    # the Abstract and Results quote them ('mean volume of the central
    # compartment and first-order clearance were 16.55 L and 41.24 L/h'), and
    # with them a typical-value solve reproduces the paper's own simulated
    # cord concentrations (Figure 4 continuous-infusion plateau, Table S1
    # cord AUC20-24) far better than the medians do. Table 2's 95% CI brackets
    # the mean, so it is precision of the mean, not between-subject spread.
    lvc <- log(16.55); label("Maternal central volume of distribution Vc (L)") # Table 2 Vc mean 16.55 (median 13.69, 95% CI 13.08-20.82, CV 40.15%)
    lcl <- log(41.24); label("Maternal clearance CL (L/h)") # Table 2 CL mean 41.24 (median 40.55, 95% CI 32.33-52.3, CV 24.06%)
    lk12 <- log(1.49); label("Central to peripheral transfer rate constant kcp (1/h)") # Table 2 kcp mean 1.49 (median 1.22, 95% CI 0.46-1.94, CV 89.35%)
    lk21 <- log(0.98); label("Peripheral to central transfer rate constant kpc (1/h)") # Table 2 kpc mean 0.98 (median 1.02, 95% CI 0.33-1.5, CV 52.72%)
    lk_central_cord <- log(0.09); label("Central to umbilical cord transfer rate constant kcu (1/h)") # Table 2 kcu mean 0.09 (median 0.07, 95% CI 0.04-0.13, CV 89.89%)
    lk_cord_central <- log(0.5); label("Umbilical cord to central transfer rate constant kuc (1/h)") # Table 2 kuc mean 0.5 (median 0.48, 95% CI 0.2-0.84, CV 70.16%)
    lvcord <- log(7.96); label("Umbilical cord compartment volume Vu (L)") # Table 2 Vu mean 7.96 (median 7.05, 95% CI 6.2-9.37, CV 51.92%)

    # Inter-individual variability. Methods: CV% = parameter SD / weighted
    # mean of the discrete NPAG distribution. Approximated here as log-normal,
    # omega^2 = log(1 + CV^2), which reproduces each printed CV exactly.
    # Table 2 'CL ~ Vc correlation' 0.67 is a correlation of the support
    # points on the linear scale; it is carried as the log-scale correlation
    # of etalcl and etalvc. No other correlation is reported.
    etalcl + etalvc ~ c(
      log(1 + 0.2406^2),
      0.67 * sqrt(log(1 + 0.2406^2) * log(1 + 0.4015^2)),
      log(1 + 0.4015^2)
    ) # Table 2 CV 24.06% (CL), 40.15% (Vc), CL ~ Vc correlation 0.67
    etalk12 ~ log(1 + 0.8935^2) # Table 2 kcp CV 89.35%
    etalk21 ~ log(1 + 0.5272^2) # Table 2 kpc CV 52.72%
    etalk_central_cord ~ log(1 + 0.8989^2) # Table 2 kcu CV 89.89%
    etalk_cord_central ~ log(1 + 0.7016^2) # Table 2 kuc CV 70.16%
    etalvcord ~ log(1 + 0.5192^2) # Table 2 Vu CV 51.92%

    # Residual error. Methods: the assay SD is the polynomial
    # SD = C0 + C1*C + C2*C^2 + C3*C^3 with C0 = 0.1 and C1 = 0.02 from the
    # assay validation and C2 = C3 = 0; a multiplicative gamma model
    # (residual error = SD x gamma) was selected (Results 'A proportional
    # error model best described the residual variability'). Data S1 Text S2
    # gives one gamma and the same polynomial row for both outputs, and
    # Table 2 reports gamma = 4.72. The resulting SD = 4.72 * (0.1 + 0.02*C)
    # is a LINEAR sum of an additive and a proportional term, i.e.
    # combined1(), with addSd = 4.72 * 0.1 mg/L and propSd = 4.72 * 0.02.
    addSd <- 4.72 * 0.1; label("Maternal plasma additive residual SD, gamma x C0 (mg/L)") # Table 2 gamma 4.72 x Methods C0 0.1
    propSd <- 4.72 * 0.02; label("Maternal plasma proportional residual SD, gamma x C1 (fraction)") # Table 2 gamma 4.72 x Methods C1 0.02
    addSd_Ccord <- 4.72 * 0.1; label("Umbilical cord additive residual SD, gamma x C0 (mg/L)") # Table 2 gamma 4.72 x Text S2 cord polynomial C0 0.10
    propSd_Ccord <- 4.72 * 0.02; label("Umbilical cord proportional residual SD, gamma x C1 (fraction)") # Table 2 gamma 4.72 x Text S2 cord polynomial C1 0.02
  })
  model({
    # No covariate effects: none were retained (Results).
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    k_central_cord <- exp(lk_central_cord + etalk_central_cord)
    k_cord_central <- exp(lk_cord_central + etalk_cord_central)
    vcord <- exp(lvcord + etalvcord)

    # Text S2 secondary variable Ke = Cl/Vc
    kel <- cl / vc

    # Equations 4-6 (Text S2): X1 central, X2 peripheral, X3 umbilical cord.
    # The cord compartment returns drug to the central compartment, so it is
    # an observed distribution space, not an elimination route.
    d/dt(central) <- -(kel + k12 + k_central_cord) * central + k21 * peripheral1 +
      k_cord_central * cord
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(cord) <- k_central_cord * central - k_cord_central * cord

    # Equations 7-8: maternal plasma (Y1) and umbilical cord blood (Y2) total
    # concentrations.
    Cc <- central / vc
    Ccord <- cord / vcord

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Ccord ~ add(addSd_Ccord) + prop(propSd_Ccord) + combined1()
  })
}
