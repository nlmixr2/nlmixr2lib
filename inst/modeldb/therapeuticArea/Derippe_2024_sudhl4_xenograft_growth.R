Derippe_2024_sudhl4_xenograft_growth <- function() {
  description <- "Preclinical (mouse). Untreated exponential tumor growth of a subcutaneous SU-DHL-4 diffuse large B-cell lymphoma xenograft, fit in Monolix to the digitized vehicle-control arm of Phillips 2020. Single first-order growth rate constant with an estimated initial tumor volume and purely additive residual error. Derippe 2024 reinterprets the growth rate constant as the per-time-step probability that a cell divides in the paper's minimal agent-based model, so this is the drug-free reference arm against which the BH3-mimetic apoptosis QSP model is bridged in vivo."
  reference <- paste(
    "Derippe T, Fouliard S, Decleves X, Mager DE.",
    "Quantitative systems pharmacology modeling of tumor heterogeneity in response to",
    "BH3-mimetics using virtual tumors calibrated with cell viability assays.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1252-1263. doi:10.1002/psp4.13158.",
    "Parameter values are from the Supporting Information section 'Control PD modeling'",
    "(supplementary file PSP4-13-1252-s001.docx).",
    "Source tumor volume data digitized from Phillips DC, Jin S, Gregory GP, et al.",
    "A novel CDK9 inhibitor increases the efficacy of venetoclax (ABT-199) in multiple",
    "models of hematologic malignancies. Leukemia. 2020;34(6):1646-1657.",
    sep = " "
  )
  vignette <- "Derippe_2024_BH3_mimetics_virtual_tumors"
  units <- list(time = "day", dosing = "mg/kg", concentration = "mm^3")

  covariateData <- list()

  covariatesDataExcluded <- list()


  compartmentData <- list(
    tumor_size = list(analyte = "SU-DHL-4 xenograft tumor volume", units = "mm^3", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "mouse (SU-DHL-4 subcutaneous xenograft)",
    n_subjects     = 1L,
    n_studies      = 1L,
    disease_state  = "Established subcutaneous SU-DHL-4 diffuse large B-cell lymphoma xenograft",
    dose_range     = "Vehicle control (no drug)",
    regions        = "Preclinical (literature-digitized)",
    notes          = "The digitized control arm is a single mean tumor volume time course (6 observations over days 0-15), so n_subjects = 1 refers to that one group-mean profile rather than to individual animals; the per-animal group size is not reported by Derippe 2024. The companion treated arms (venetoclax 50 mg/kg PO QD, A-1592668 1.5 mg/kg PO three times weekly, and the combination) are simulated by the paper's agent-based model rather than by this ODE model. No inter-individual variability was reported."
  )

  ini({
    # Supplement "Control PD modeling": "The exponential growth rate constant was
    # estimated to 0.1313 h-1 (RSE = 2.2%), the initial tumor volume was estimated
    # to 258 mm3 (RSE = 3.7%). The error model was purely additive (24.28 mm3, RSE = 29%)."
    #
    # The printed "h-1" unit is a typo: 0.1313 1/h would be a 5.3-hour tumor
    # doubling time. The rate is per DAY. This is settled by the deposited source
    # data (data/mice_SU_DHL4_full.csv), whose control arm goes from 255.9 to
    # 1819.9 mm3 over 15 days -- an ordinary least squares fit of log(volume) on
    # day gives 0.136 1/day and an intercept of 243 mm3, matching the reported
    # Monolix values. The agent-based model in the same Supplement also consumes
    # the constant as "probability = kgrowth x step" with step = 0.1 DAYS.
    # See vignette Errata.
    lkge   <- log(0.1313); label("First-order exponential tumor growth rate constant (1/day)")  # Supplement "Control PD modeling" (2.2% RSE); unit corrected from the printed 1/h
    lrbase <- log(258);    label("Initial tumor volume (mm^3)")                                  # Supplement "Control PD modeling" (3.7% RSE)

    addSd <- 24.28; label("Additive residual SD on tumor volume (mm^3)")                         # Supplement "Control PD modeling" (29% RSE)
  })

  model({
    kge   <- exp(lkge)
    rbase <- exp(lrbase)

    tumor_size(0) <- rbase
    d/dt(tumor_size) <- kge * tumor_size

    tumor_size ~ add(addSd)
  })
}
