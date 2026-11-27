Lakshminarayana_2015_cipargamin_auctre <- function() {
  description <- paste(
    "Preclinical (mouse, NMRI, Plasmodium berghei-infected). Sigmoid (variable-slope)",
    "PK-PD index model relating day-4 percent parasitemia to the AUC0-48/TRE index of",
    "oral cipargamin (KAE609) in the Lakshminarayana 2015 dose-fractionation study,",
    "where TRE = 2 x the ex vivo P. berghei IC99 = 48 ng/mL. Fitted in GraphPad Prism",
    "on log10 parasitemia against log10 index:",
    "log10 P = log10 Pmin + (log10 Pmax - log10 Pmin) / (1 + 10^((log10 EC50 - log10 x) * s)).",
    "There is no PK component: the 0-48 h plasma AUC enters as the covariate",
    "AUC_CIPARGAMIN and is divided by the threshold parameter tre. The paper's own AUC values",
    "were simulated with the PK model Lakshminarayana_2015_cipargamin_mouse. Sibling",
    "fits of the same data to Cmax/TRE and %T>TRE are",
    "Lakshminarayana_2015_cipargamin_cmaxtre and Lakshminarayana_2015_cipargamin_ttre."
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
    dosing = "ug*h/mL (cipargamin 0-48 h plasma AUC, supplied as the AUC_CIPARGAMIN covariate)",
    concentration = "% parasitemia (observation; parasitized erythrocytes per 100 erythrocytes)"
  )

  covariateData <- list(
    AUC_CIPARGAMIN = list(
      description = "Cipargamin plasma AUC from the first dose to 48 h (AUC0-48), the numerator of the AUC0-48/TRE index",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma concentration. The window is fixed at 0-48 h after the first",
        "dose, the interval between the start of treatment (24 h post-infection) and",
        "the parasitemia readout (72 h post-infection), whatever the number of doses.",
        "In the source the value was a typical-value PK-model simulation for each",
        "dose-fractionation regimen (Lakshminarayana_2015_cipargamin_mouse). Set to 0",
        "for vehicle controls so the predicted parasitemia equals Pmax."
      ),
      source_name = "AUC0-48 (Table 1 tabulates AUC0-48/TRE)"
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
    n_subjects = 60L,
    n_studies = 1L,
    disease_state = paste(
      "Plasmodium berghei GFP ANKA (PbGFPCON) blood-stage infection, intravenous",
      "inoculation with 2 x 10^7 parasitized erythrocytes on day 0"
    ),
    weight_range = "20-22 g",
    sex_female_pct = 100,
    dose_range = paste(
      "Total daily oral doses of 0.5, 1, 2, 4 and 8 mg/kg given once, as two doses",
      "12 h apart, or as three doses 8 h apart, starting 24 h post-infection"
    ),
    sampling = "Parasitemia by flow cytometry at 72 h post-infection (48 h after the first dose)",
    regions = "Switzerland (Swiss Tropical and Public Health Institute, Basel)",
    notes = paste(
      "15 regimens x 4 mice (Table 1); vehicle-treated controls had a mean",
      "parasitemia of 35.51% +/- 3.86%. The index values were computed with",
      "TRE = 48 ng/mL (2 x IC99 of 61.6 nM, i.e. 123 nM)."
    )
  )

  ini({
    lpmin <- log(0.13)
    label("Minimum percent parasitemia Pmin (%)") # Table 5, AUC0-48/TRE column, Pmin = 0.13 (0.09 to 0.20)
    lpmax <- log(32.66)
    label("Maximum percent parasitemia Pmax (%)") # Table 5, AUC0-48/TRE column, Pmax = 32.66 (25.06 to 42.56)
    lec50 <- log(52.72)
    label("AUC0-48/TRE index producing half the maximal log10-parasitemia reduction EC50 (h)") # Table 5, AUC0-48/TRE column, EC50 = 52.72 (44.16 to 62.95)
    hill <- -1.63
    label("Hill slope s on the log10 index scale (unitless; negative = parasitemia falls with exposure)") # Table 5, AUC0-48/TRE column, Hill slope = -1.63 (-2.04 to -1.22)

    tre <- fixed(0.048)
    label("PK-PD threshold TRE = 2 x P. berghei IC99 (ug/mL)") # Results: IC99 = 61.6 nM; 2 x IC99 = 123 nM or 48 ng/ml used as TRE

    # Results: standard deviation of the regression Sy.x = 0.22 (log10
    # parasitemia units); converted to the natural-log scale for lnorm().
    expSd <- 0.22 * log(10)
    label("Log-normal residual SD on parasitemia (natural-log scale)")
  })

  model({
    pmin <- exp(lpmin)
    pmax <- exp(lpmax)
    ec50 <- exp(lec50)

    # Equation 1 with x = AUC0-48/TRE. With s < 0,
    # 1 / (1 + 10^((log10 EC50 - log10 x) * s)) = EC50^h / (EC50^h + x^h), h = -s,
    # which stays defined at x = 0 (vehicle control).
    auc_tre <- AUC_CIPARGAMIN / tre
    h <- -hill
    log10_par <- log10(pmin) + (log10(pmax) - log10(pmin)) * ec50^h / (ec50^h + auc_tre^h)
    parasites <- 10^log10_par

    parasites ~ lnorm(expSd)
  })
}
