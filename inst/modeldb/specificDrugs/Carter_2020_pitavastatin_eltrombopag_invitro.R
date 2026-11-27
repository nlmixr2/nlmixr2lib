Carter_2020_pitavastatin_eltrombopag_invitro <- function() {
  description <- paste(
    "In vitro (cryopreserved suspended human hepatocytes, 10-donor pool).",
    "Micro-rate-constant mechanistic model of hepatocyte uptake of",
    "pitavastatin (substrate) and eltrombopag (inhibitor) with mutual",
    "competitive inhibition at a shared uptake transporter (Carter 2020",
    "model 1, the final model). Each drug has three states holding an amount",
    "per 1e6 cells: incubation medium, bound to the uptake transporter, and",
    "intracellular. Drug in the medium binds the free transporter by",
    "second-order association (free transporter = total transporter minus",
    "both drug-bound amounts), dissociates back to the medium, or is",
    "translocated into the cell; both drugs also cross the membrane passively",
    "in both directions, and intracellular pitavastatin is eliminated",
    "(metabolised) first order. The observed cellular concentration is the",
    "transporter-bound plus intracellular amount divided by the cell volume.",
    "Parameters are the Table 2 median individual estimates; no",
    "between-experiment variance or residual-error magnitude is reported."
  )
  reference <- paste(
    "Carter SJ, Chouhan B, Sharma P, Chappell MJ. Prediction of Clinical",
    "Transporter-Mediated Drug-Drug Interactions via Comeasurement of",
    "Pitavastatin and Eltrombopag in Human Hepatocyte Models. CPT",
    "Pharmacometrics Syst Pharmacol. 2020;9(4):211-221.",
    "doi:10.1002/psp4.12505 (PMID 32142598, PMC7179958).",
    "ODEs: Supplemental Table S1 (Eqs. S10-S17) and the authors' Monolix",
    "model code 'Pita_Etm_Micro_Model_Comp' (Model Codes supplement).",
    "Parameter estimates: Table 2 ('Model 1: Micro-rate constant parameter",
    "estimates')."
  )
  vignette <- "Carter_2020_pitavastatin_eltrombopag"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/mL")

  dosing <- c("medium", "medium_eltrombopag")

  # Per-drug states of the hepatocyte incubation. Pitavastatin (the
  # substrate) keeps the bare names and eltrombopag (the inhibitor) carries
  # the `_eltrombopag` suffix, following the victim / perpetrator convention
  # of the transporter-DDI models. All six hold amounts per 1e6 cells.
  paper_specific_compartments <- c(
    "medium",
    "transporter_bound",
    "intracellular",
    "medium_eltrombopag",
    "transporter_bound_eltrombopag",
    "intracellular_eltrombopag"
  )

  compartmentData <- list(
    medium = list(analyte = "pitavastatin", units = "nmol", specimen = "administration site", verified = TRUE),
    transporter_bound = list(analyte = "pitavastatin", units = "nmol", specimen = "tissue", verified = TRUE),
    intracellular = list(analyte = "pitavastatin", units = "nmol", specimen = "tissue", verified = TRUE),
    medium_eltrombopag = list(
      analyte = "eltrombopag",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transporter_bound_eltrombopag = list(analyte = "eltrombopag", units = "nmol", specimen = "tissue", verified = TRUE),
    intracellular_eltrombopag = list(analyte = "eltrombopag", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (cryopreserved suspended human hepatocytes)",
    n_subjects = 35L,
    n_studies = 1L,
    donors = "BioIVT LiverPool lot LYB: 10 donors (8 Caucasian, 1 African American, 1 Hispanic)",
    system = "Suspended hepatocytes (viability 84-87 percent) at 1e6 cells/mL in Leibovitz L15 medium, 37 C; cells separated from the medium by an oil-spin method",
    dose_range = paste(
      "Pitavastatin 0.3, 1, 3, 10, 30 and 100 nmol/mL with or without",
      "eltrombopag; eltrombopag pre-incubated for 15 min at 45 nmol/mL",
      "(1.5e6 cells/mL) and diluted to 30 nmol/mL (1e6 cells/mL) by the",
      "pitavastatin addition. Per 1e6 cells both phases hold 30 nmol of",
      "eltrombopag."
    ),
    duration = "Sampling 0.25-30 min after pitavastatin addition",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Each pitavastatin concentration in each of three experiments was",
      "fitted as an individual data set in Monolix 2018R2 (35 data sets in",
      "the deposited data file: 18 pitavastatin-only and 17 with",
      "eltrombopag, the 100 nmol/mL eltrombopag arm of one experiment being",
      "excluded). Table 2 gives the conditional mode with the",
      "minimum-maximum over the individual data sets; log-normal parameters",
      "and a proportional residual error were used, but the variances were",
      "not reported. The initial eltrombopag amounts at the time of",
      "pitavastatin addition were fitted as free parameters per data set",
      "and are not reported."
    )
  )

  ini({
    # Pitavastatin (substrate)
    lkon <- log(0.17); label("Pitavastatin association rate constant to the uptake transporter, k_aP (1/nmol/min)") # Table 2, k_aX for pitavastatin 0.17 (0.14-0.25)
    lkoff <- log(2.2); label("Pitavastatin dissociation rate constant from the transporter back to the medium, k_dP (1/min)") # Table 2, k_dX for pitavastatin 2.2 (1.97-2.37)
    lktrans <- log(1.65); label("Pitavastatin translocation rate constant from the transporter into the cell, k_tP (1/min)") # Table 2, k_tX for pitavastatin 1.65 (1.57-1.74)
    lkmet <- log(0.22); label("Intracellular pitavastatin elimination (metabolism) rate constant, k_eP (1/min)") # Table 2, k_eP 0.22 (0.2-0.24)
    lkpass_in <- log(5.5e-4); label("Passive rate constant of pitavastatin from the medium into the cell, k_fP (1/min)") # Table 2, k_fX for pitavastatin 5.5e-4 (4.6e-4 to 6.1e-4)
    lkpass_out <- log(0.21); label("Passive rate constant of pitavastatin from the cell into the medium, k_bP (1/min)") # Table 2, k_bX for pitavastatin 0.21 (0.18-0.22)

    # Shared uptake transporter
    ltransporter_tot <- log(0.18); label("Total amount of uptake transporter per 1e6 cells, T_o (nmol)") # Table 2, T_o 0.18 (0.11-0.37)

    # Eltrombopag (inhibitor)
    lkon_eltrombopag <- log(0.26); label("Eltrombopag association rate constant to the uptake transporter, k_aE (1/nmol/min)") # Table 2, k_aX for eltrombopag 0.26 (0.23-0.31)
    lkoff_eltrombopag <- log(1.57); label("Eltrombopag dissociation rate constant from the transporter back to the medium, k_dE (1/min)") # Table 2, k_dX for eltrombopag 1.57 (1.42-2)
    lktrans_eltrombopag <- log(0.27); label("Eltrombopag translocation rate constant from the transporter into the cell, k_tE (1/min)") # Table 2, k_tX for eltrombopag 0.27 (0.24-0.32)
    lkpass_in_eltrombopag <- log(0.05); label("Passive rate constant of eltrombopag from the medium into the cell, k_fE (1/min)") # Table 2, k_fX for eltrombopag 0.05 (0.04-0.06)
    lkpass_out_eltrombopag <- log(0.43); label("Passive rate constant of eltrombopag from the cell into the medium, k_bE (1/min)") # Table 2, k_bX for eltrombopag 0.43 (0.35-0.65)

    # Experimental constant
    vcell <- fixed(0.003); label("Cell volume per 1e6 hepatocytes, V_cell (mL)") # Methods 'Data analysis': V_cell 3 uL per 1e6 human hepatocytes; Eqs. S16-S17

    # Residual error: proportional in Monolix, magnitude not reported
    propSd_Ccell <- fixed(0); label("Proportional residual error on the cellular pitavastatin concentration (fraction)") # Methods 'Parameter estimation': proportional residual error; value not reported
    propSd_Ccell_eltrombopag <- fixed(0); label("Proportional residual error on the cellular eltrombopag concentration (fraction)") # Methods 'Parameter estimation': proportional residual error; value not reported
  })

  model({
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ktrans <- exp(lktrans)
    kmet <- exp(lkmet)
    kpass_in <- exp(lkpass_in)
    kpass_out <- exp(lkpass_out)
    transporter_tot <- exp(ltransporter_tot)
    kon_eltrombopag <- exp(lkon_eltrombopag)
    koff_eltrombopag <- exp(lkoff_eltrombopag)
    ktrans_eltrombopag <- exp(lktrans_eltrombopag)
    kpass_in_eltrombopag <- exp(lkpass_in_eltrombopag)
    kpass_out_eltrombopag <- exp(lkpass_out_eltrombopag)

    # Free transporter by conservation (Eq. S3): T_f = T_o - S2 - I2
    transporter_free <- transporter_tot - transporter_bound - transporter_bound_eltrombopag

    # Pitavastatin, Eqs. S10, S12, S14
    d/dt(medium) <- -kon * medium * transporter_free - kpass_in * medium + koff * transporter_bound + kpass_out * intracellular
    d/dt(transporter_bound) <- kon * medium * transporter_free - (koff + ktrans) * transporter_bound
    d/dt(intracellular) <- kpass_in * medium + ktrans * transporter_bound - (kpass_out + kmet) * intracellular

    # Eltrombopag, Eqs. S11, S13, S15
    d/dt(medium_eltrombopag) <- -kon_eltrombopag * medium_eltrombopag * transporter_free - kpass_in_eltrombopag * medium_eltrombopag + koff_eltrombopag * transporter_bound_eltrombopag + kpass_out_eltrombopag * intracellular_eltrombopag
    d/dt(transporter_bound_eltrombopag) <- kon_eltrombopag * medium_eltrombopag * transporter_free - (koff_eltrombopag + ktrans_eltrombopag) * transporter_bound_eltrombopag
    d/dt(intracellular_eltrombopag) <- kpass_in_eltrombopag * medium_eltrombopag + ktrans_eltrombopag * transporter_bound_eltrombopag - kpass_out_eltrombopag * intracellular_eltrombopag

    # Observed cellular concentrations, Eqs. S16-S17 (nmol/mL)
    Ccell <- (transporter_bound + intracellular) / vcell
    Ccell_eltrombopag <- (transporter_bound_eltrombopag + intracellular_eltrombopag) / vcell

    Ccell ~ prop(propSd_Ccell)
    Ccell_eltrombopag ~ prop(propSd_Ccell_eltrombopag)
  })
}
