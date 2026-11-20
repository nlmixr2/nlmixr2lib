Singley_2026_gepotidacin_Ecoli_ALL <- function() {
  description <- paste(
    "Preclinical (mouse, neutropenic male CD-1; thigh infection with",
    "Escherichia coli ALL). Gepotidacin PK/PD model of Singley 2026 for",
    "E. coli isolate ALL (Figure S8). The typical-value one-compartment",
    "subcutaneous PK of Singley_2026_gepotidacin gives the free-drug",
    "AUC over the first 24 h of treatment (fauc_0_24; the PK model",
    "predicts unbound concentrations directly); a three-parameter",
    "inhibitory Imax model (Hill coefficient fixed at 1, Table S4",
    "footnote) of fAUC0-24/MIC then gives the change in thigh bacterial",
    "burden between the start of treatment and the end of the study",
    "(dlog10cfu, log10 CFU per thigh, positive = net growth), and the",
    "observed output log_cfu adds that change to the nominal 6.5 log10",
    "CFU starting burden. Both are the paper's 24 h endpoint (doses",
    "every 6 h from 1 h to 19 h after infection, thighs harvested at 25",
    "h) and are meant to be read at t = 24 h after the first dose;",
    "before 24 h they reflect only the exposure accumulated so far.",
    "Gepotidacin MIC of this isolate 4 ug/mL (Table 3); resistance:",
    "cephems, FQ, carbapenems; NDM-1. Mean growth of vehicle-treated",
    "controls over the study +1.52 log10 CFU (Tables 2 and S5/S6).",
    "Published fAUC/MIC targets for stasis, 1-log10 and 2-log10 kill:",
    "2.9, 6.7 and 14.8. No residual error was reported, so addSd is",
    "held at 0."
  )
  reference <- paste(
    "Singley CM, Kurumaddali A, Hoover JL. Pharmacokinetic and",
    "pharmacodynamic characterization of gepotidacin against Escherichia",
    "coli and Klebsiella pneumoniae in a neutropenic mouse thigh infection",
    "model. Antimicrob Agents Chemother. 2026;70(1):e01176-25.",
    "doi:10.1128/aac.01176-25. PK from Supplemental Material Table S1; PD",
    "parameters from Figure S8 panel 'ALL' fAUC/MIC table."
  )
  vignette <- "Singley_2026_gepotidacin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL unbound (Cc); log10 CFU per thigh (log_cfu)"
  )

  paper_specific_compartments <- c("fauc_0_24")

  covariateData <- list(
    DOSE = list(
      description = "Gepotidacin dose given at each subcutaneous administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Dose per administration in mg per mouse (the unit of amt),",
        "entering ka as (DOSE / 2.7)^-0.437 (Figure S1). Convert a mg/kg",
        "dose by multiplying by the body weight in kg (27-28 g mice; 2.7 mg",
        "is 100 mg/kg in a 27 g mouse). For the every-6-h regimen it is the",
        "size of each dose, not the daily total. For an untreated control",
        "give no dose records and any positive DOSE (e.g. 2.7): the dose",
        "term is undefined at DOSE = 0, and with no drug fauc_0_24 stays 0",
        "so dlog10cfu equals e0. Place the event columns before DOSE in the",
        "event table."
      ),
      source_name = "Dose (Figure S1 'Individual Ka = Kapop * (Dose/2.7)^Deff * e^eta')"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "gepotidacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gepotidacin (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    fauc_0_24 = list(analyte = "gepotidacin (unbound)", units = "ug*h/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, outbred, male; neutropenic)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_range = "27-28 g",
    sex_female_pct = 0,
    organism = "Escherichia coli ALL; gepotidacin MIC 4 ug/mL by CLSI broth microdilution",
    disease_state = paste(
      "Neutropenic (cyclophosphamide 150 mg/kg on day -4 and 100 mg/kg on",
      "day -1) mice infected intramuscularly in one thigh; inoculum tailored",
      "per isolate for a 1 h baseline of about 6.5 log10 CFU per thigh"
    ),
    dose_range = "1-600 mg/kg subcutaneously every 6 h for 24 h (four doses), from 1 h after infection",
    regions = "United States (GSK)",
    notes = paste(
      "Five mice per efficacy group; one group sacrificed at 1 h after",
      "infection for the baseline count and the rest at 25 h. The",
      "fAUC/MIC of each group was derived by the authors from",
      "typical-value mrgsolve simulations of the population PK model,",
      "which this file reproduces. Fitted model diagnostics: AIC 74.1 and",
      "R^2 0.9. Lower and upper limits of quantification of the thigh",
      "counts 1.2 and 9.4 log10 CFU."
    )
  )

  ini({
    # Population PK typical values (Table S1); the authors computed the
    # PK/PD index from typical-value simulations, so no between-subject
    # variability is carried.
    lka <- log(0.886)
    label("Absorption rate constant ka at the 2.7 mg reference dose (1/h)") # Table S1 'Kapop (1/h)' = 0.886
    lcl <- log(0.101)
    label("Clearance CL (L/h)") # Table S1 'CLpop (L/h)' = 0.101
    lvc <- log(0.018)
    label("Central volume of distribution V (L)") # Table S1 'Vpop (L)' = 0.018
    e_dose_ka <- -0.437
    label("Power exponent of dose per administration (DOSE/2.7 mg) on ka (unitless)") # Table S1 'Dose effect on Kapop' = -0.437

    # Isolate MIC, a measured input to the index.
    mic <- fixed(4)
    label("Broth-microdilution MIC of E. coli ALL (ug/mL)") # Table 3 'Gepotidacin MIC (ug/mL)' = 4

    # Inhibitory Imax model of fAUC0-24/MIC (Figure S5):
    # E = E0 - Imax * X / (IC50 + X), Hill coefficient fixed at 1.
    # Published targets: stasis 2.9, 1-log 6.7, 2-log 14.8; AIC 74.1, R^2 0.9.
    e0 <- 1.5
    label("Change in thigh burden of vehicle-treated controls, e0 (log10 CFU)") # Figure S8 panel 'ALL' fAUC/MIC table 'E0' = 1.5 (SE 0.18)
    limax <- log(5.2)
    label("Maximum change from e0, imax (log10 CFU)") # Figure S8 panel 'ALL' fAUC/MIC table 'Imax' = 5.2 (SE 0.34)
    lic50 <- log(7.3)
    label("fAUC0-24/MIC at half the maximum effect, ic50 (h)") # Figure S8 panel 'ALL' fAUC/MIC table 'IC50' = 7.3 (SE 1.73)
    lhill <- fixed(log(1))
    label("Hill coefficient of the fAUC0-24/MIC effect curve (unitless)") # Table S4 footnote a: slope (gamma) set as a constant value of 1

    # Nominal burden at the start of treatment, a design input.
    log10_cfu0 <- fixed(6.5)
    label("Bacterial burden at the start of treatment (log10 CFU per thigh)") # Materials and Methods 'Infection': inoculum tailored to a 1 h baseline of approximately 6.5 log10 CFU/mouse

    addSd <- fixed(0)
    label("Additive residual SD on log_cfu (log10 CFU per thigh)") # not reported (AIC and R^2 only)
  })

  model({
    ka <- exp(lka) * (DOSE / 2.7)^e_dose_ka
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc

    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Unbound plasma concentration (ug/mL).
    Cc <- central / vc

    # Free-drug AUC over the first 24 h after the first dose; the state is
    # constant after 24 h, so the index holds its 24 h value from then on.
    d/dt(fauc_0_24) <- Cc * (t < 24)
    fauc_mic <- fauc_0_24 / mic

    # Change in log10 CFU per thigh from the start of treatment to the end
    # of the study; e0 at zero exposure, approaching e0 - imax at high
    # exposure.
    dlog10cfu <- e0 - imax * fauc_mic^hill / (fauc_mic^hill + ic50^hill)
    log_cfu <- log10_cfu0 + dlog10cfu
    log_cfu ~ add(addSd)
  })
}
