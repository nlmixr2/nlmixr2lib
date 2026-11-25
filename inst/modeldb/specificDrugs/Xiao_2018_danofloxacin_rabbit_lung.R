Xiao_2018_danofloxacin_rabbit_lung <- function() {
  description <- paste0(
    "Preclinical (rabbit, New Zealand). Inhibitory sigmoid Emax PK/PD-index model for the in vivo ",
    "antibacterial effect of oral danofloxacin against Salmonella Typhimurium strain HJ0310J in an ",
    "intravenous-challenge systemic infection model, driven by the AUC24h/MIC index. This file carries ",
    "the fit for the LUNG bacterial load; the other three organ fits are packaged as sibling ",
    "models. Xiao 2018 fitted E = E0 + Emax * Ce^N / (EC50^N + Ce^N) in WinNonlin 6.1, where E is the ",
    "change in log10 CFU/g of tissue over the treatment course and Ce is AUC24h/MIC. Table 2 reports Emax as ",
    "the positive magnitude of the maximal kill (4.74) while Figure 4 plots E falling from 0 ",
    "and Table 2 calls the 228.66 target a '3 Log10 cfu/mL reduction', so the effect is encoded in ",
    "the inhibitory form E = e0 - imax * Ce^N / (EC50^N + Ce^N) with imax = Table 2 Emax; this ",
    "reproduces the Table 2 target (a 3.00 log10 reduction from e0 at AUC24h/MIC = 228.66) and the Figure 4 plateau. ",
    "There is NO PK component: Xiao 2018 analysed serum danofloxacin non-compartmentally (Table 1) ",
    "and published no structural PK model, so exposure enters as the covariate AUC_DANOFLOXACIN_SERUM (serum ",
    "AUC0-24h after one oral dose) divided by the parameter mic. Xiao 2018 reports no residual SD ",
    "or between-animal variability for the fit, so no eta parameters are present and addSd is FIXED ",
    "at 0 (typical-value simulation only). Sibling models: Xiao_2018_danofloxacin_rabbit_liver, ",
    "Xiao_2018_danofloxacin_rabbit_spleen, Xiao_2018_danofloxacin_rabbit_lung."
  )
  reference <- paste(
    "Xiao X, Pei L, Jiang L-J, Lan W-X, Xiao J-Y, Jiang Y-J, Wang Z-Q. (2018).",
    "In Vivo Pharmacokinetic/Pharmacodynamic Profiles of Danofloxacin in Rabbits",
    "Infected With Salmonella typhimurium After Oral Administration.",
    "Frontiers in Pharmacology 9:391.",
    "doi:10.3389/fphar.2018.00391. PMCID: PMC5913287.",
    "Model equation: Materials and Methods, 'Pharmacokinetics and Pharmacodynamics",
    "Analysis'. Parameter estimates and 3-log10 targets: Table 2, Lung column.",
    "Non-compartmental PK: Table 1. AUC24h/MIC per dose: Results, 'In Vivo",
    "Pharmacokinetic and Pharmacodynamic Profiles of Danofloxacin'. MIC: Results,",
    "'In Vitro Pharmacodynamics'. Fitted curve and data: Figure 4.",
    sep = " "
  )
  vignette <- "Xiao_2018_danofloxacin_rabbit"

  units <- list(
    time = "h",
    dosing = "ug*h/mL (danofloxacin serum AUC0-24h, supplied as a covariate)",
    concentration = "log10 CFU/g (observation; the signed change in lung S. Typhimurium load over the 3-day treatment course)"
  )

  depends <- c("AUC_DANOFLOXACIN_SERUM")

  covariateData <- list(
    AUC_DANOFLOXACIN_SERUM = list(
      description = "Danofloxacin serum area under the concentration-time curve from 0 to 24 h after a single oral dose, the numerator of the AUC24h/MIC PK/PD index driving the in vivo effect",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Xiao 2018 measured AUC24h non-compartmentally at 1, 10 and 30 mg/kg (Table 1: 1.55, 11.70 ",
        "and 28.10 ug*h/mL) and obtained the other efficacy-study doses from a linear regression of ",
        "AUC24h on dose (r^2 = 0.9935; least-squares refit AUC24h = 1.4885 + 0.8996 * dose). The ",
        "AUC24h/MIC values listed in the Results text (76.42 to 911.24 for 1 to 30 mg/kg) equal that ",
        "regression divided by mic = 0.03125 ug/mL at all eight doses, including 1, 10 and 30 mg/kg; ",
        "Figure 4, however, plots the 1, 10 and 30 mg/kg groups at the MEASURED AUC24h/MIC (about 50, ",
        "374 and 899). See the vignette. Set to 0 for the untreated control so the sigmoid term ",
        "vanishes and the predicted change equals e0."
      ),
      source_name = "AUC24 h (Xiao 2018 Table 1; Materials and Methods, 'Pharmacokinetics and Pharmacodynamics Analysis')"
    )
  )

  population <- list(
    species = "rabbit (New Zealand), infected with Salmonella Typhimurium",
    n_subjects = 45L,
    n_studies = 1L,
    age_range = "5-6 weeks",
    weight_range = "0.8-1.2 kg",
    organism = "Salmonella enterica serovar Typhimurium strain HJ0310J (virulent field strain, Yangzhou University). Danofloxacin MIC 0.03 ug/mL (two-fold dilution step 0.03125) and MBC 0.06 ug/mL in Mueller-Hinton broth; 0.06 and 0.12 ug/mL in rabbit serum; MPC 0.25 ug/mL. The broth MIC is the denominator of the AUC24h/MIC index",
    system = "Systemic infection by two intravenous (ear marginal vein) injections of 0.5 mL of a 10^8 CFU/mL suspension 12 h apart. Pre-treatment loads 12 h after the last inoculation (infection-model study): 6.51 +/- 0.32 (blood), 6.56 +/- 0.26 (liver), 7.73 +/- 0.18 (spleen) and 7.32 +/- 0.32 (lung) log10 CFU/mL or CFU/g. Limit of detection 2 x 10^3 CFU/mL or CFU/g",
    disease_state = "experimental systemic Salmonella Typhimurium infection (sepsis model)",
    dose_range = "0 (0.85% NaCl control), 1, 2.5, 5, 7.5, 10, 15, 20 and 30 mg/kg danofloxacin by gavage once daily for 3 days, starting 12 h after infection; blood, liver, spleen and lung sampled 24 h after the last dose",
    regions = "China (Yangzhou University)",
    notes = "115 rabbits in total: 20 for the inoculum-finding infection model, 20 for the mortality study, 30 for pharmacokinetics (single oral dose of 1, 10 or 30 mg/kg, n = 10 per dose, serum sampled 0.25-24 h) and 45 for the efficacy study (9 groups of 5) that this model was fitted to. Sex not reported. The four organ fits (blood, liver, spleen, lung) are reported side by side in Table 2 and are packaged as four model files."
  )

  ini({
    # =================================================================
    # Xiao 2018 Table 2, Lung column
    # =================================================================
    # Printed equation (Materials and Methods, 'Pharmacokinetics and
    # Pharmacodynamics Analysis'):
    #   E = E0 + Emax * Ce^N / (EC50^N + Ce^N)
    # Table 2 reports Emax as a POSITIVE magnitude, but Figure 4 plots
    # E falling from about 0, and Table 2 labels the AUC24h/MIC
    # of 228.66 the target for a '3 Log10 cfu/mL reduction'. The fit
    # is therefore the WinNonlin inhibitory form
    #   E = E0 - Imax * Ce^N / (EC50^N + Ce^N),  Imax = Table 2 Emax,
    # whose inhibitory term equals 3.00 at Ce = 228.66 (vignette gate).
    e0 <- 0.01
    label("Change in lung S. Typhimurium load at zero exposure E0 (log10 CFU/g)") # Xiao 2018 Table 2, Lung column, E0 = 0.01
    limax <- log(4.74)
    label("Log maximal reduction in lung S. Typhimurium load Imax (log10 CFU/g; Emax in Xiao 2018)") # Xiao 2018 Table 2, Lung column, Emax = 4.74
    lec50 <- log(148.18)
    label("Log AUC24h/MIC producing half the maximal reduction EC50 (h)") # Xiao 2018 Table 2, Lung column, EC50 = 148.18
    lhill <- log(1.26)
    label("Log Hill coefficient N of the index-effect curve (unitless)") # Xiao 2018 Table 2, Lung column, Slope (N) = 1.26

    # =================================================================
    # Minimum inhibitory concentration (index denominator)
    # =================================================================
    # Printed as 0.03 ug/mL, the rounded label of the 1/32 two-fold
    # dilution step. 0.03125 is the value the paper divided by: the
    # AUC24h regression on dose divided by 0.03125 reproduces all eight
    # AUC24h/MIC values in the Results to within 0.03, whereas 0.03
    # overshoots them by 3.2-38.0. FIXED: a measured strain property.
    mic <- fixed(0.03125)
    label("Danofloxacin MIC against S. Typhimurium HJ0310J in MH broth (ug/mL)") # Xiao 2018 Results, 'In Vitro Pharmacodynamics': MIC 0.03 ug/mL in MH broth (dilution step 0.03125; see vignette)

    # =================================================================
    # Residual error
    # =================================================================
    # Xiao 2018 reports r^2 = 0.8708 for this fit and no residual SD or
    # standard errors.
    addSd <- fixed(0)
    label("Additive residual SD on the change in log10 CFU/g (0; not reported in Xiao 2018)") # Xiao 2018 Results reports r^2 only
  })

  model({
    imax <- exp(limax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC24h/MIC (units: h)
    aucmic <- AUC_DANOFLOXACIN_SERUM / mic

    # Signed change in lung log10 CFU/g over the treatment course
    # (negative = reduction): e0 at zero exposure, e0 - imax at
    # saturating exposure.
    Cc <- e0 - imax * aucmic^hill / (ec50^hill + aucmic^hill)
    Cc ~ add(addSd)
  })
}
