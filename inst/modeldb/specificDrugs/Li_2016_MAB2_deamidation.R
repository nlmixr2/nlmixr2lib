Li_2016_MAB2_deamidation <- function() {
  description <- paste(
    "Mechanistic two-variant, two-compartment PK model for the Fc-site",
    "deamidation quality attribute of MAB2, a non-disclosed humanized IgG1",
    "antibody (Biogen), in 3 patients after a single high IV dose (Li 2016).",
    "Original (non-deamidated) and deamidated MAB2 each follow a linear",
    "two-compartment model with shared distribution rate constants and",
    "volume; the original form converts first-order to the deamidated form",
    "in both the central and the peripheral compartment, and both forms are",
    "cleared at the same rate. The initial deamidated fraction of the dose",
    "is a fixed model input. Typical-value model: no IIV or residual error",
    "was reported."
  )
  reference <- paste(
    "Li Y, Monine M, Huang Y, Swann P, Nestorov I, Lyubarskaya Y.",
    "Quantitation and pharmacokinetic modeling of therapeutic antibody",
    "quality attributes in human studies. MAbs. 2016;8(6):1079-1087.",
    "doi:10.1080/19420862.2016.1186322. Model equations 1-4 and the dose",
    "split from the Methods section 'Pharmacokinetic modeling'; parameter",
    "values from Table 1 (Deamidation model column).",
    sep = " "
  )
  vignette <- "Li_2016_MAB2_quality_attributes"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "MAB2, original (non-deamidated) form",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "MAB2, original (non-deamidated) form",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_deamid = list(analyte = "MAB2, Fc-deamidated form", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_deamid = list(analyte = "MAB2, Fc-deamidated form", units = "mg", specimen = "tissue", verified = TRUE)
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
      "the same initial deamidation level. Antibody recovered by",
      "anti-idiotypic affinity purification and quantified by Lys-C",
      "peptide mapping LC-MS against a heavy-isotope internal calibrant.",
      "Parameters estimated in ADAPT; the CVs printed in Table 1 are",
      "estimate precision, not between-subject variability. No",
      "demographics are reported."
    )
  )

  ini({
    lk12 <- log(0.17); label("Central-to-peripheral rate constant kCP, both forms (1/day)") # Table 1 Deamidation model: kCP = 0.17 /day (CV 35.9%)
    lk21 <- log(0.093); label("Peripheral-to-central rate constant kPC, both forms (1/day)") # Table 1 Deamidation model: kPC = 0.093 /day (CV 1.4%)
    lkmet_deamid <- log(0.01); label("Conversion rate constant ktr from original to deamidated form (1/day)") # Table 1 Deamidation model: ktr = 0.01 /day (CV 3.1%)
    lkel <- log(0.13); label("Elimination rate constant kcl0 of both forms (1/day)") # Table 1 Deamidation model: kcl0 = 0.13 /day (CV 10.8%); kclmod = kcl0 by Methods assumption 1
    lvc <- log(3.34); label("Central volume V, both forms (L)") # Table 1 Deamidation model: V = 3.34 L (CV 26.2%)
    fdose_deamid <- fixed(0.056); label("Deamidated fraction f of the administered dose (fraction)") # Table 2 and Figure 6A: initial level 5.6 percent, the base of the 1x/2x/3x/5x series
  })

  model({
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kmet_deamid <- exp(lkmet_deamid)
    kel <- exp(lkel)
    # Methods assumption 1: Fc-site deamidation does not change clearance
    kel_deamid <- kel
    vc <- exp(lvc)

    # Equations 1-4; conversion occurs in both the central and the
    # peripheral compartment at the same rate constant
    d/dt(central) <- -(k12 + kmet_deamid + kel) * central + k21 * peripheral1
    d/dt(peripheral1) <- -(k21 + kmet_deamid) * peripheral1 + k12 * central
    d/dt(central_deamid) <- -(k12 + kel_deamid) * central_deamid + kmet_deamid * central + k21 * peripheral1_deamid
    d/dt(peripheral1_deamid) <- -k21 * peripheral1_deamid + kmet_deamid * peripheral1 + k12 * central_deamid

    # A dose is split between the two forms: AC0(0) = (1 - f) * Dose and
    # ACmod(0) = f * Dose. Give the full dose amount twice at the same
    # time, once into central and once into central_deamid.
    f(central) <- 1 - fdose_deamid
    f(central_deamid) <- fdose_deamid

    # Total MAB2, deamidated MAB2, and the relative deamidation level (%)
    Cc <- (central + central_deamid) / vc
    Cc_deamid <- central_deamid / vc
    pct_deamid <- 100 * central_deamid / (central + central_deamid)
  })
}
