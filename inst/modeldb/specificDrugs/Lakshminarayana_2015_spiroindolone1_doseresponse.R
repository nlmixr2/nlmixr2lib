Lakshminarayana_2015_spiroindolone1_doseresponse <- function() {
  description <- paste(
    "Preclinical (mouse, NMRI, Plasmodium berghei-infected). Sigmoid dose-response",
    "model for the day-4 percent parasitemia after a single oral dose of",
    "spiroindolone analog (+)-1 of the Lakshminarayana 2015 series (an active (+)-enantiomer from the Novartis cipargamin lead-optimisation campaign; structure in supplementary Figure S1),",
    "given 24 h after infection (Lakshminarayana 2015 Table 3). Fitted in NONMEM VI",
    "(FOCE) on log10 parasitemia against log10 dose with equation 1,",
    "log10 P = log10 Pmin + (log10 Pmax - log10 Pmin) / (1 + 10^((log10 EZ50 - log10 x) * s)),",
    "where x is the dose in mg/kg. Table 3 reports ED90 rather than EZ50; the model",
    "derives EZ50 = ED90 / 9^(1 / |s|), the reading under which ED90 is the dose",
    "giving 90% of the maximal log10-parasitemia reduction (checked against the",
    "Figure 1 curves). There is no PK component: the dose enters as the covariate",
    "DOSE (mg/kg). Sibling dose-response fits for the other analogs and for cipargamin share this structure (Lakshminarayana_2015_<compound>_doseresponse)."
  )
  reference <- paste(
    "Lakshminarayana SB, Freymond C, Fischli C, Yu J, Weber S, Goh A, Yeung BKS,",
    "Ho PC, Dartois V, Diagana TT, Rottmann M, Blasco F. (2015).",
    "Pharmacokinetic-pharmacodynamic analysis of spiroindolone analogs and KAE609",
    "in a murine malaria model. Antimicrob Agents Chemother 59(2):1200-1210.",
    "doi:10.1128/AAC.03274-14. PMCID: PMC4335872."
  )
  vignette <- "Lakshminarayana_2015_spiroindolones"
  units <- list(
    time = "h",
    dosing = "mg/kg (single oral dose, supplied as the DOSE covariate)",
    concentration = "% parasitemia (observation; parasitized erythrocytes per 100 erythrocytes)"
  )

  covariateData <- list(
    DOSE = list(
      description = "Single oral dose of (+)-1 given 24 h post-infection",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE register entry: the animal's assigned dose level,",
        "here per kilogram of body weight. Set to 0 for vehicle controls so the",
        "predicted parasitemia equals Pmax. The parasitemia readout is 48 h after",
        "the dose (72 h post-infection), so the model is a single-time-point",
        "dose-response and does not describe a time course."
      ),
      source_name = "x (equation 1; dose in mg/kg)"
    )
  )

  compartmentData <- list(
    parasites = list(
      analyte = "Plasmodium berghei GFP ANKA-infected erythrocytes",
      units = "% parasitemia",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (NMRI, female, specific-pathogen-free)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Plasmodium berghei GFP ANKA (PbGFPCON) blood-stage infection, intravenous",
      "inoculation with 2 x 10^7 parasitized erythrocytes on day 0"
    ),
    weight_range = "20-22 g",
    sex_female_pct = 100,
    dose_range = "Single oral doses of 2.5, 5, 10, 30 and 100 mg/kg, in 10% ethanol, 30% PEG 400, 60% (10% vitamin E TPGS), 10 mL/kg",
    sampling = "Parasitemia by flow cytometry at 72 h post-infection (48 h after the dose)",
    regions = "Switzerland (Swiss Tropical and Public Health Institute, Basel)",
    notes = paste(
      "Groups of five mice per dose level; Figure 1 pools several studies (one",
      "symbol per study), so the total number of animals is not stated. Untreated",
      "controls typically died on day 6-7. In vitro IC50 (Table 2): 13 +/- 2.2 nM",
      "against P. falciparum NF54 and 242 +/- 124 nM against P. berghei."
    )
  )

  ini({
    lpmin <- log(0.03)
    label("Minimum percent parasitemia Pmin (%)") # Table 3, (+)-1 row, Pmin = 0.03
    lpmax <- log(33.11)
    label("Maximum percent parasitemia Pmax (%)") # Table 3, (+)-1 row, Pmax = 33.11
    led90 <- log(15.3)
    label("Dose giving 90% of the maximal log10-parasitemia reduction ED90 (mg/kg)") # Table 3, (+)-1 row, ED90 = 15.3 mg/kg
    hill <- -3.1
    label("Hill slope s on the log10 dose scale (unitless; negative = parasitemia falls with dose)") # Table 3, (+)-1 row, Hill slope = -3.1

    # Table 3 reports no residual error magnitude for the NONMEM fit.
    expSd <- fixed(0)
    label("Log-normal residual SD on parasitemia (natural-log scale; not reported, held at 0)")
  })

  model({
    pmin <- exp(lpmin)
    pmax <- exp(lpmax)
    ed90 <- exp(led90)
    h <- -hill

    # EZ50 from ED90: 90% of the maximal effect is reached where
    # (x / EZ50)^h = 9, so EZ50 = ED90 / 9^(1 / h).
    ez50 <- ed90 / 9^(1 / h)

    # Equation 1, written in the algebraically equivalent Hill form
    # EZ50^h / (EZ50^h + x^h) so that it stays defined at DOSE = 0.
    log10_par <- log10(pmin) + (log10(pmax) - log10(pmin)) * ez50^h / (ez50^h + DOSE^h)
    parasites <- 10^log10_par

    parasites ~ lnorm(expSd)
  })
}
