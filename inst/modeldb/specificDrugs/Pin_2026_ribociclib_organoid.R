Pin_2026_ribociclib_organoid <- function() {
  description <- paste(
    "In vitro (human small-intestinal organoids, hSIO). QSP. Pin 2026",
    "organoid cell-dynamics model of ribociclib-induced intestinal toxicity:",
    "proliferative (MCM2+), non-proliferative and dead hSIO cells plus toxic",
    "cell-produced waste accumulating in an unrefreshed medium. Ribociclib",
    "(CDK4/6 inhibitor) kills proliferative cells through a Hill function",
    "of the free drug concentration in the hSIO cells, supplied as the",
    "time-varying covariate CEFFECT (uM; taken equal to the measured free",
    "medium concentration). Deterministic typical-value fit to one 96 h",
    "experiment (no IIV, no residual error). Time is in hours.",
    "The same kill function drives the human epithelium model",
    "Pin_2026_ribociclib_human."
  )
  reference <- paste(
    "Pin C, Maheshvare M D, Gall L, Hall A, Coen M, Bianco A, Liu S,",
    "Pryde K. (2026). Intestinal Organoid-Based Mathematical Modeling",
    "Predicts Clinical Gastrointestinal Toxicity of Oral Oncology Drugs.",
    "CPT Pharmacometrics Syst Pharmacol 15(5):e70256.",
    "doi:10.1002/psp4.70256.",
    "Model code deposited as Supporting Information Data S2 and at BioModels",
    "(MODEL2212120003)."
  )
  vignette <- "Pin_2026_intestinal_organoid_gi_toxicity"

  paper_specific_compartments <- c("prolif", "nonprolif", "dead", "waste")

  units <- list(
    time = "h",
    dosing = "n/a (the drug exposure is supplied externally via the time-varying covariate CEFFECT; the model has no dosing compartment)",
    concentration = "uM for CEFFECT; cells for the hSIO cell states; arbitrary units for waste"
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Time-varying free ribociclib concentration in the hSIO cells, the site of action of the kill (and arrest) functions.",
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Pin 2026 Eqs 5-6 take the free drug in the hSIO cells equal to the",
        "free drug in the medium (rapid equilibration; negligible protein",
        "binding in the serum-free medium). The medium concentration was",
        "measured over the 96 h experiment; supply it at each record and",
        "interpolate linearly. Set to 0 for a vehicle control. The paper's",
        "symbol is Drug."
      ),
      source_name = "Drug"
    )
  )

  compartmentData <- list(
    prolif = list(analyte = "proliferative (MCM2+) hSIO cells", units = "cells", specimen = "tissue", verified = TRUE),
    nonprolif = list(
      analyte = "non-proliferative (MCM2-) hSIO cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    dead = list(analyte = "dead hSIO cells not yet detached", units = "cells", specimen = "tissue", verified = TRUE),
    waste = list(
      analyte = "toxic cell-produced waste in the culture medium",
      units = "arbitrary units",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (human small-intestinal organoids)",
    n_subjects = "One adult stem-cell-derived hSIO line from one donor (HUB Organoid Technology biobank); 3-4 replicate wells per concentration.",
    n_studies = 1L,
    disease_state = "Healthy human small-intestinal epithelium (organoid culture)",
    dose_range = "four concentrations covering clinically relevant ranges (Figure S1B); vehicle control on the same plates",
    regions = "In vitro (AstraZeneca)",
    notes = paste(
      "hSIOs in serum-free organoid growth medium, sampled at 0, 24, 48, 72",
      "and 96 h without medium refreshment; MCM2+ proliferative and MCM2-",
      "non-proliferative cells counted by immunofluorescence image analysis",
      "(Pin 2026 Methods and Supporting Information). Parameters fitted",
      "sequentially in MATLAB (nlinfit, log-transformed counts): control",
      "cell dynamics first, then the toxicity function."
    )
  )

  ini({
    # hSIO cell dynamics -- Pin 2026 Table 1, Ribociclib column.
    lambda_p <- 0.0334
    label("Proliferation rate of proliferative cells lambda_P (1/h)") # Table 1, Ribociclib column: 0.0334
    delta_pn <- 0.0076
    label("Transfer rate proliferative to non-proliferative delta_PN (1/h)") # Table 1, Ribociclib column: 0.0076
    delta_nd <- 0.0152
    label("Transfer rate non-proliferative to dead delta_ND (1/h)") # Table 1, Ribociclib column: 0.0152
    tau_d <- fixed(log(2) / 12)
    label("Detachment rate of dead cells tau_D (1/h)") # Table 1: ln(2)/12, fixed [ref 16]
    kwaste <- fixed(0.1)
    label("Cell waste production rate k_waste (a.u. per cell per h)") # Table 1: 0.1, fixed
    kkillwaste <- 0.02
    label("Waste-induced cell killing rate k_killwaste (1/h)") # Table 1: 0.02
    c50waste <- 59874.14
    label("Waste level inducing half-maximal killing C_50waste (a.u.)") # Table 1: 59874.14 (exp(11) in the deposited code)
    phi <- fixed(0)
    label("Clearance rate of waste from the medium phi (1/h)") # Table 1: 0, fixed (no medium refreshment)
    pref <- 4316
    label("Reference number of proliferative MCM2+ cells p_ref (cells)") # Table 1, Ribociclib column: 4316
    nref <- 1.3e7
    label("Reference number of non-proliferative cells n_ref (cells)") # Table 1, Ribociclib column: 1.3 x 10^7

    # Drug effect -- Pin 2026 Table 1 'Drug effects', Eqs 5-6.
    karrest <- 0
    label("Maximal drug-induced cell-cycle arrest k_arrest (unitless)") # Table 1, Ribociclib column: 0 (no arrest detected)
    c50arrest <- fixed(1)
    label("Free drug concentration for half-maximal arrest C_50arrest (uM)") # inactive while k_arrest = 0; value 1 is the deposited-code placeholder (HIOs_DynamicModel_Abema.m), not reported in Table 1
    hillarrest <- fixed(1)
    label("Hill coefficient of the arrest function (unitless)") # inactive while k_arrest = 0; deposited-code placeholder 1
    kkill <- 0.8969
    label("Maximal drug-induced killing of proliferative cells k_kill (unitless)") # Table 1, Ribociclib column: 0.8969
    c50kill <- 10.0385
    label("Free drug concentration for half-maximal killing C_50kill (uM)") # Table 1, Ribociclib column: 10.0385
    hillkill <- fixed(1)
    label("Hill coefficient of the kill function (unitless)") # Table 1: 1, assumed

    # Initial conditions -- Pin 2026 Table 1 'Initial values', Ribociclib column.
    p0 <- 997
    label("Proliferative MCM2+ cells on day 0 P(0) (cells)") # Table 1, Ribociclib column: 997
    n0 <- 222
    label("Non-proliferative cells on day 0 N(0) (cells)") # Table 1, Ribociclib column: 222
    d0 <- 1
    label("Dead cells on day 0 D(0) (cells)") # Table 1: 1
    w0 <- 0
    label("Cell waste in the medium on day 0 W(0) (a.u.)") # Table 1: 0
  })
  model({
    # Drug-induced arrest and killing of proliferative cells (Eqs 5-6).
    farrest <- karrest * CEFFECT^hillarrest / (CEFFECT^hillarrest + c50arrest^hillarrest)
    fkill <- kkill * CEFFECT^hillkill / (CEFFECT^hillkill + c50kill^hillkill)

    # Waste-induced killing, a Hill function of the accumulated waste.
    wastekill <- kkillwaste * waste / (waste + c50waste)

    # hSIO cell dynamics (Eqs 1-4).
    d/dt(prolif) <- lambda_p * prolif * (1 - farrest - fkill) -
      delta_pn * prolif * prolif / pref - wastekill * prolif
    d/dt(nonprolif) <- delta_pn * prolif * prolif / pref + lambda_p * prolif * farrest -
      delta_nd * nonprolif * nonprolif / nref - wastekill * nonprolif
    d/dt(dead) <- delta_nd * nonprolif * nonprolif / nref + lambda_p * prolif * fkill -
      tau_d * dead + wastekill * (prolif + nonprolif)
    d/dt(waste) <- kwaste * (prolif + nonprolif) - phi * waste

    prolif(0) <- p0
    nonprolif(0) <- n0
    dead(0) <- d0
    waste(0) <- w0

    # Figure 1 plots log10(1 + cells).
    lprolif <- log10(1 + prolif)
    lnonprolif <- log10(1 + nonprolif)
  })
}
