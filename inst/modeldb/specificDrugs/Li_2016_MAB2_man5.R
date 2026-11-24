Li_2016_MAB2_man5 <- function() {
  description <- paste(
    "Mechanistic two-variant, two-compartment PK model for the Man5",
    "(mannose-5) Fc glycoform quality attribute of MAB2, a non-disclosed",
    "humanized IgG1 antibody (Biogen), in 3 patients after a single high IV",
    "dose (Li 2016). Man5-free and Man5-containing MAB2 each follow a linear",
    "two-compartment model with shared distribution rate constants and",
    "volume; Man5 MAB2 is cleared 2.3 times faster than the rest and is not",
    "formed in circulation. The initial Man5 fraction of the dose is a",
    "fixed model input. Typical-value model: no IIV or residual error was",
    "reported."
  )
  reference <- paste(
    "Li Y, Monine M, Huang Y, Swann P, Nestorov I, Lyubarskaya Y.",
    "Quantitation and pharmacokinetic modeling of therapeutic antibody",
    "quality attributes in human studies. MAbs. 2016;8(6):1079-1087.",
    "doi:10.1080/19420862.2016.1186322. Model equations 1-4 and the dose",
    "split from the Methods section 'Pharmacokinetic modeling'; parameter",
    "values from Table 1 (Man5 model column).",
    sep = " "
  )
  vignette <- "Li_2016_MAB2_quality_attributes"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "MAB2, forms without the Man5 Fc glycan",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "MAB2, forms without the Man5 Fc glycan",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_man5 = list(analyte = "MAB2, Man5 Fc glycoform", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_man5 = list(analyte = "MAB2, Man5 Fc glycoform", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 3L,
    n_studies = 1L,
    disease_state = paste(
      "Patients in a single-ascending-dose PK study of MAB2 (indication not",
      "disclosed)."
    ),
    dose_range = paste(
      "Single high IV dose, the same in all 3 patients; the dose amount is",
      "not reported."
    ),
    regions = "Not reported",
    notes = paste(
      "Li 2016 Methods 'Clinical sample information': 3 patients, serum",
      "sampled pre-dose, 10 min, 1 h, 4 h, 12 h, and days 2, 4, 7, 14, 21",
      "and 42 (samples to 24 weeks were collected, only the first 6 weeks",
      "were profiled). One drug-product batch, so every patient received",
      "the same initial Man5 level. Antibody recovered by anti-idiotypic",
      "affinity purification and quantified by Lys-C peptide mapping LC-MS",
      "against a heavy-isotope internal calibrant. Parameters estimated in",
      "ADAPT; the CVs printed in Table 1 are estimate precision, not",
      "between-subject variability. No demographics are reported."
    )
  )

  ini({
    lk12 <- log(0.18); label("Central-to-peripheral rate constant kCP, both forms (1/day)") # Table 1 Man5 model: kCP = 0.18 /day (CV 57.4%)
    lk21 <- log(0.096); label("Peripheral-to-central rate constant kPC, both forms (1/day)") # Table 1 Man5 model: kPC = 0.096 /day (CV 7.3%)
    lkel <- log(0.127); label("Elimination rate constant kcl0 of the Man5-free form (1/day)") # Table 1 Man5 model: kcl0 = 0.127 /day (CV 4.6%)
    lkel_man5 <- log(0.295); label("Elimination rate constant kclmod of the Man5 glycoform (1/day)") # Table 1 Man5 model: kclmod = 0.295 /day (CV 8.1%)
    lvc <- log(3.32); label("Central volume V, both forms (L)") # Table 1 Man5 model: V = 3.32 L (CV 28.6%)
    fdose_man5 <- fixed(0.030); label("Man5 fraction f of the administered dose (fraction)") # Table 2 and Figure 6B: initial level 3.0 percent, the base of the 1x/2x/3x/5x series
  })

  model({
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    kel_man5 <- exp(lkel_man5)
    vc <- exp(lvc)

    # Equations 1-4 with ktr = 0 (Table 1; Methods assumption 2: Man5 does
    # not form in circulation), so the two forms do not exchange
    d/dt(central) <- -(k12 + kel) * central + k21 * peripheral1
    d/dt(peripheral1) <- -k21 * peripheral1 + k12 * central
    d/dt(central_man5) <- -(k12 + kel_man5) * central_man5 + k21 * peripheral1_man5
    d/dt(peripheral1_man5) <- -k21 * peripheral1_man5 + k12 * central_man5

    # A dose is split between the two forms: AC0(0) = (1 - f) * Dose and
    # ACmod(0) = f * Dose. Give the full dose amount twice at the same
    # time, once into central and once into central_man5.
    f(central) <- 1 - fdose_man5
    f(central_man5) <- fdose_man5

    # Total MAB2, Man5 MAB2, and the relative Man5 level (%)
    Cc <- (central + central_man5) / vc
    Cc_man5 <- central_man5 / vc
    pct_man5 <- 100 * central_man5 / (central + central_man5)
  })
}
