Lakshminarayana_2015_cipargamin_ttre <- function() {
  description <- paste(
    "Preclinical (mouse, NMRI, Plasmodium berghei-infected). Sigmoid (variable-slope)",
    "PK-PD index model relating day-4 percent parasitemia to the %T>TRE index of",
    "oral cipargamin (KAE609) in the Lakshminarayana 2015 dose-fractionation study,",
    "where TRE = 2 x the ex vivo P. berghei IC99 = 48 ng/mL. Fitted in GraphPad Prism",
    "on log10 parasitemia against the untransformed index x (percent of the 48 h",
    "window above TRE): log10 P = log10 Pmin + (log10 Pmax - log10 Pmin) / (1 + 10^((EC50 - x) * s)).",
    "There is no PK component: the percentage of time above TRE enters as the covariate",
    "FTTRE_CIPARGAMIN, already formed against the threshold. The paper's own values",
    "were simulated with the PK model Lakshminarayana_2015_cipargamin_mouse. Sibling",
    "fits of the same data to Cmax/TRE and AUC0-48/TRE are",
    "Lakshminarayana_2015_cipargamin_cmaxtre and Lakshminarayana_2015_cipargamin_auctre."
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
    dosing = "percent (time above TRE over 0-48 h, supplied as the FTTRE_CIPARGAMIN covariate)",
    concentration = "% parasitemia (observation; parasitized erythrocytes per 100 erythrocytes)"
  )

  covariateData <- list(
    FTTRE_CIPARGAMIN = list(
      description = "Percentage of the 48 h after the first dose during which the total cipargamin plasma concentration exceeded TRE (2 x the P. berghei IC99 = 48 ng/mL)",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Percent scale, not a fraction: the fitted EC50 (32.06) is on the same",
        "scale. The window is fixed at 0-48 h after the first dose. In the source",
        "the value was computed from typical-value PK-model simulations for each",
        "dose-fractionation regimen (Lakshminarayana_2015_cipargamin_mouse) against",
        "TRE = 48 ng/mL. Set to 0 for vehicle controls."
      ),
      source_name = "%T>TRE (Table 1)"
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
    lpmin <- log(0.15)
    label("Minimum percent parasitemia Pmin (%)") # Table 5, %T>TRE column, Pmin = 0.15 (0.12 to 0.20)
    lpmax <- log(53.09)
    label("Maximum percent parasitemia Pmax (%)") # Table 5, %T>TRE column, Pmax = 53.09 (32.28 to 87.50)
    lec50 <- log(32.06)
    label("Percent time above TRE producing half the maximal log10-parasitemia reduction EC50 (%)") # Table 5, %T>TRE column, EC50 = 32.06 (28.04 to 36.07)
    hill <- -0.031
    label("Hill slope s on the untransformed %T>TRE scale (1/%; negative = parasitemia falls with exposure)") # Table 5, %T>TRE column, Hill slope = -0.031 (-0.037 to -0.024)

    # Results: standard deviation of the regression Sy.x = 0.17 (log10
    # parasitemia units); converted to the natural-log scale for lnorm().
    expSd <- 0.17 * log(10)
    label("Log-normal residual SD on parasitemia (natural-log scale)")
  })

  model({
    pmin <- exp(lpmin)
    pmax <- exp(lpmax)
    ec50 <- exp(lec50)

    # Equation 1 with x = %T>TRE entered untransformed (Methods: the
    # percentage index is not log-transformed), so the curve is a logistic
    # in x rather than a Hill function of x.
    log10_par <- log10(pmin) + (log10(pmax) - log10(pmin)) / (1 + 10^((ec50 - FTTRE_CIPARGAMIN) * hill))
    parasites <- 10^log10_par

    parasites ~ lnorm(expSd)
  })
}
