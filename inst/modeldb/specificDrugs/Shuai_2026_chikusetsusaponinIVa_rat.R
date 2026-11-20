Shuai_2026_chikusetsusaponinIVa_rat <- function() {
  description <- paste(
    "Preclinical (rat). One-compartment pharmacokinetic model with first-order",
    "absorption (no lag time) for chikusetsusaponin IVa (CS-IVa) given as a",
    "single 15 mg/kg intraperitoneal dose at the onset of reperfusion to male",
    "Sprague-Dawley rats with myocardial ischemia/reperfusion injury (30 min",
    "left-anterior-descending ligation) (Shuai 2026). The model was fitted by",
    "unweighted least squares in Phoenix WinNonlin to the 13-point mean",
    "concentration-time profile of 8 rats, so it carries no between-animal",
    "variability. An exploratory effect compartment with the representative",
    "ke0 = 0.001 1/min selected by the authors' ke0 sensitivity grid is",
    "attached; the concentration-effect function linking the effect-site",
    "concentration to the CK-MB response was not reported, so the model outputs",
    "the effect-site concentration only."
  )
  reference <- paste(
    "Shuai X, Wang H, Luo J, Zeng Y, Wang Y, Zhu L, Liu Z, Cheng Y (2026).",
    "Pharmacokinetics and Exploratory Exposure-Response Analysis of",
    "Chikusetsusaponin IVa in Myocardial Ischemia/Reperfusion-Injured Rats.",
    "Pharmaceuticals (Basel) 19(5):749. doi:10.3390/ph19050749"
  )
  vignette <- "Shuai_2026_chikusetsusaponinIVa_rat"
  units <- list(time = "min", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "chikusetsusaponin IVa",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "chikusetsusaponin IVa",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "chikusetsusaponin IVa",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 8,
    n_studies = 1,
    sex_female_pct = 0,
    weight_range = "200 +/- 20 g",
    disease_state = paste(
      "myocardial ischemia/reperfusion injury (30 min left anterior descending",
      "coronary artery ligation followed by reperfusion)"
    ),
    dose_range = "15 mg/kg single intraperitoneal dose at the onset of reperfusion",
    regions = "China (Guangzhou University of Chinese Medicine)",
    n_observations = "13 mean plasma concentrations, 5-1440 min post-dose",
    notes = paste(
      "Shuai 2026 Section 4.2 and 4.5: specific-pathogen-free male",
      "Sprague-Dawley rats weighing 200 +/- 20 g; MI/RI-15 mg/kg group n = 8.",
      "Normal rats (7.5, 15 and 30 mg/kg, n = 8 per group) were analysed by",
      "non-compartmental analysis only (Table 2); no compartmental model was",
      "fitted to them. Plasma CS-IVa was measured by UHPLC-MS/MS (LLOQ 1.22 nM,",
      "about 0.00097 mg/L at MW 795.0 g/mol). The structural fit used the group",
      "mean profile: the AIC and SBC in Table S6 are reproduced only with",
      "n = 13 observations and 3 parameters."
    )
  )

  ini({
    # One-compartment, first-order absorption, no lag, unweighted (Shuai 2026
    # Section 2.4; Supplementary Table S6 model 1 and Table S7 'Full data' row).
    lka <- log(0.01518); label("First-order absorption rate constant K01 (1/min)") # Table S7, full data, K01 = 0.01518 1/min (CV 40.57)
    lcl <- log(0.03589); label("Apparent clearance CL/F (L/min/kg)") # Table S7, full data, CL/F = 0.03589 L/min/kg (CV 11.33)
    lvc <- log(5.4017); label("Apparent volume of distribution V/F (L/kg)") # Table S7, full data, V/F = 5.4017 L/kg (CV 32.18)

    # Exploratory effect compartment (Section 2.5; Supplementary S7 and Table S8):
    # ke0 was not estimated but chosen from the grid 0.1-0.00001 1/min as the
    # best-fitting converged candidate.
    lke0 <- fixed(log(0.001)); label("Effect-compartment equilibration rate constant ke0 (1/min)") # Table S8 and Figure 6C, ke0 = 0.001 1/min

    # Unweighted least squares = additive residual error. Not printed; derived
    # as sqrt(SSR / (n - p)) = sqrt(0.154555 / (13 - 3)), with SSR from Table S6
    # model 1 and n = 13 implied by AIC = n * ln(SSR) + 2p = -18.2736.
    addSd <- 0.1243; label("Additive residual SD on CS-IVa plasma concentration (mg/L)") # derived from Table S6 model 1 SSR and AIC
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ke0 <- exp(lke0)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(effect) <- ke0 * (central / vc - effect)

    Cc <- central / vc
    Ce <- effect

    Cc ~ add(addSd)
  })
}
