Bastida_2018_tocilizumab <- function() {
  description <- "One-compartment population PK model for intravenous tocilizumab in adults with rheumatoid arthritis (Bastida 2018), with parallel first-order linear and Michaelis-Menten elimination from the central compartment; total body weight and time-varying C-reactive protein on linear CL."
  reference <- "Bastida C, Ruiz-Esquide V, Pascal M, de Vries Schultink AHM, Yague J, Sanmarti R, Huitema ADR, Soy D. Fixed dosing of intravenous tocilizumab in rheumatoid arthritis. Results from a population pharmacokinetic analysis. Br J Clin Pharmacol. 2018;84(4):716-725. doi:10.1111/bcp.13500"
  vignette <- "Bastida_2018_tocilizumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 62 kg per the Bastida 2018 final-model equation CL = 0.0104 * (WT/62)^0.360 * (1 + 0.131 * (CRP - 0.484)) (Results, p719). The reference 62 kg is close to the cohort mean weight of 63.5 kg (Table 1). Treated as baseline-fixed; the source paper does not describe time-varying weight handling.",
      source_name        = "WT"
    ),
    CRP = list(
      description        = "C-reactive protein concentration (time-varying)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear additive-offset effect on linear CL with reference 0.484 mg/dL per the Bastida 2018 final-model equation CL = 0.0104 * (WT/62)^0.360 * (1 + 0.131 * (CRP - 0.484)) (Results, p719; equation in paper uses the Spanish abbreviation PCR for C-reactive protein, identical to CRP in the English text). The paper assesses CRP as a time-varying covariate (Methods, p717; Discussion, p722) so this column is intended to be supplied at every observation time. Assay type is standard CRP (not high-sensitivity). The reference value 0.484 mg/dL is the cohort median; the cohort baseline mean from Table 1 is 0.29 mg/dL (this lower mean reflects that most patients in the cohort were already in remission at inclusion). The canonical CRP register entry documents the unit as mg/L; this model carries the source unit mg/dL (1 mg/dL = 10 mg/L) to match the paper's equation exactly without unit conversion, consistent with the Wang_2020_ontamalimab.R precedent.",
      source_name        = "PCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_observations = 109L,
    n_studies      = 1L,
    age_range      = "18 years and older",
    age_mean_sd    = "54.1 +/- 12.3 years",
    weight_range   = "estimated 40-120 kg (simulated dose-range; cohort mean 63.5 +/- 13.8 kg)",
    weight_mean_sd = "63.5 +/- 13.8 kg",
    sex_female_pct = 88.6,
    race_ethnicity = c("White/Caucasian" = 80.0, Hispanic = 17.1, "Afro-American" = 2.9),
    disease_state  = "Adult rheumatoid arthritis (74% in DAS28 remission, 11% low and 14% moderate disease activity at inclusion). Erosive RA in 82.9% of subjects; 71.4% RF-positive; 77.1% anti-CCP-positive.",
    dose_range     = "4, 6, or 8 mg/kg by 1-hour intravenous infusion every 28 days at the discretion of the treating rheumatologist.",
    regions        = "Single-center prospective observational study, Hospital Clinic Barcelona, Spain (protocol HCB/2015/0533).",
    notes          = paste(
      "Baseline demographics from Bastida 2018 Table 1 (p718-719).",
      "Sample distribution: 54/109 trough concentrations and 55/109 intermediate samples (approximately 7, 14, 21 days post-infusion).",
      "16/35 patients contributed a single sample; the remaining 19 contributed an average of five samples each.",
      "ADA-positive: 0/17 tested (assay limited to samples with tocilizumab < 1 ug/mL).",
      "Concomitant therapy: corticosteroids 57.1%, conventional DMARDs 68.6% (methotrexate 75% of DMARD users, leflunomide 20.8%, sulfasalazine 4.2%); 65.7% previously treated with biologic agents.",
      "Smoking status: 62.9% non-smokers, 25.7% ex-smokers, 11.4% active smokers.",
      "Baseline laboratory mean (SD): serum creatinine 0.69 (0.29) mg/mL; haemoglobin 134 (14) g/L; total protein 68.1 (4.59) g/L; albumin 43.5 (1.98) g/L; HDL-cholesterol 61.8 (13.8) mg/dL; CRP 0.29 (0.76) mg/dL; ESR 4.69 (6.95) mm.",
      sep = " "
    )
  )

  ini({
    # Structural PK parameters - Bastida 2018 Table 2 final-model estimates (p720).
    # Reference covariate values for the typical subject: total body weight 62 kg,
    # CRP 0.484 mg/dL. The paper reports time in hours; CL is L/h and Vmax is mg/h.
    lcl   <- log(0.0104); label("Linear clearance CL (L/h)")                                # Bastida 2018 Table 2, CL row
    lvc   <- log(4.83);   label("Central volume of distribution V (L)")                     # Bastida 2018 Table 2, V row
    lvmax <- log(0.239);  label("Maximum Michaelis-Menten elimination rate Vmax (mg/h)")    # Bastida 2018 Table 2, VM row
    lkm   <- log(4.22);   label("Michaelis-Menten constant Km (ug/mL)")                     # Bastida 2018 Table 2, KM row

    # Covariate effects on linear CL - Bastida 2018 Table 2 ("Covariate effects" section)
    # combined with the closed-form equation reproduced in Results (p719):
    #   CL = 0.0104 * (WT/62)^0.360 * (1 + 0.131 * (CRP - 0.484))
    # Weight enters as a power-of-ratio (WT/62)^e_wt_cl with reference 62 kg;
    # CRP enters as a linear additive offset (1 + e_crp_cl * (CRP - 0.484)) with
    # reference 0.484 mg/dL. The reference 62 kg and 0.484 mg/dL are read directly
    # off the equation in the paper.
    e_wt_cl  <- 0.360; label("Power exponent of WT/62 on linear CL (unitless)")             # Bastida 2018 Table 2, WT(kg) on CL
    e_crp_cl <- 0.131; label("Linear coefficient of (CRP - 0.484 mg/dL) on linear CL (1/(mg/dL))")  # Bastida 2018 Table 2, CRP (mg dL-1) on CL

    # Inter-individual variability - Bastida 2018 Table 2 ("Interindividual variability"
    # section) reports CV% on the linear-parameter scale for CL and V; no off-diagonal
    # correlation is reported. Convert each CV% to log-normal variance via
    # omega^2 = log(1 + CV^2):
    #   CL CV 17.0% -> omega^2 = log(1 + 0.170^2) = 0.02850
    #   V  CV 30.8% -> omega^2 = log(1 + 0.308^2) = 0.09065
    # Bastida 2018 Table 2 VAR(eta_CL) = 17.0% CV
    etalcl ~ 0.02850
    # Bastida 2018 Table 2 VAR(eta_V) = 30.8% CV
    etalvc ~ 0.09065

    # Residual error - Bastida 2018 Table 2 (Residual error section). Combined
    # additive + proportional model: Cobs = Cpred * (1 + eps_prop) + eps_add.
    propSd <- 0.255; label("Proportional residual error (fraction)")                        # Bastida 2018 Table 2 Proportional row, 25.5%
    addSd  <- 0.161; label("Additive residual error (ug/mL)")                               # Bastida 2018 Table 2 Additive row, 0.161 ug/mL
  })
  model({
    # Individual PK parameters with the Bastida 2018 final-model covariate equations.
    # CL = 0.0104 * (WT/62)^0.360 * (1 + 0.131 * (CRP - 0.484))   [Results, p719]
    cl   <- exp(lcl + etalcl) *
            (WT / 62)^e_wt_cl *
            (1 + e_crp_cl * (CRP - 0.484))
    vc   <- exp(lvc + etalvc)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # One-compartment IV PK with parallel linear and Michaelis-Menten elimination
    # from the central compartment. Concentration in ug/mL; central in mg
    # (1 mg/L = 1 ug/mL). Dose enters central directly via the 1-hour IV infusion;
    # there is no depot or peripheral compartment in the Bastida 2018 model.
    Cc <- central / vc

    d/dt(central) <- -(cl / vc) * central -
                      vmax * Cc / (km + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
