Jin_2025_dorzagliatin_wx390 <- function() {
  description <- paste(
    "Preclinical (mouse). Sigmoid Emax model of the percent reduction in fasting",
    "plasma glucose produced by oral dorzagliatin (a glucokinase activator) given",
    "with the PI3K/mTOR inhibitor WX390, relative to the time-matched",
    "WX390-alone control group, driven by an",
    "externally supplied dorzagliatin plasma concentration",
    "(CP_DORZAGLIATIN_NGML). Fit in female BALB/c nude mice after 28 days of",
    "dorzagliatin 10 mg/kg once daily by oral gavage, given 30 min before",
    "WX390 0.2 mg/kg (Jin 2025 Table 5, row 'Dorz + WX390'). PD-only model:",
    "Jin 2025 fitted a two-compartment first-order",
    "absorption PK model in Phoenix but published none of its parameter values,",
    "so no PK layer is packaged and users must supply the concentration",
    "trajectory. Companion files Jin_2025_dorzagliatin (dorzagliatin alone) and",
    "Jin_2025_dorzagliatin_alpelisib (with alpelisib, BYL719) carry the other",
    "two separately fitted curves.",
    sep = " "
  )

  reference <- paste(
    "Jin G, Zheng K, Liu S, Yi H, Wei W, Xu C, Xiang X, Kang Y. In Vivo PK-PD and",
    "Drug-Drug Interaction Study of Dorzagliatin for the Management of PI3Kalpha Inhibitor-Induced",
    "Hyperglycemia. Pharmaceuticals (Basel). 2025;18(6):927.",
    "doi:10.3390/ph18060927.",
    sep = " "
  )

  vignette <- "Jin_2025_dorzagliatin_pi3ki_hyperglycemia"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external dorzagliatin plasma-concentration covariate)",
    concentration = "(observation fpg_reduction is a percent, 0-100; driving covariate CP_DORZAGLIATIN_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_DORZAGLIATIN_NGML = list(
      description = paste(
        "Instantaneous dorzagliatin plasma concentration at the time of each PD",
        "observation, supplied as a time-varying covariate from observed plasma",
        "samples or an upstream PK source."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the sigmoid Emax expression",
        "fpg_reduction = e0 + emax * C^hill / (ec50^hill + C^hill) (Jin 2025",
        "Equation 6, where Ce is 'the plasma concentration of the drug'; no",
        "effect-compartment delay is reported). At C = 0 the output equals e0,",
        "not 0."
      ),
      source_name = "C / Ce (Jin 2025 Equation 6 and Table 5)"
    )
  )

  population <- list(
    species = "mouse (female BALB/c nude)",
    n_subjects = 3L,
    n_studies = 1,
    age_range = "6-8 weeks at purchase",
    weight_range = "18-20 g at purchase",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state = paste(
      "Nude mice of the OVCAR3 ovarian-cancer xenograft programme; the glucose",
      "pharmacodynamic study was fasted from 6 h before the last dose."
    ),
    dose_range = "Dorzagliatin 10 mg/kg orally once daily for 28 days, 30 min before WX390 0.2 mg/kg",
    regions = "China",
    notes = paste(
      "PD study (Jin 2025 Section 4.5): 18 female BALB/c nude mice in six groups",
      "of n = 3 (vehicle, WX390, BYL719, dorzagliatin, dorzagliatin + WX390,",
      "dorzagliatin + BYL719); this file is the Dorz + WX390 group. On",
      "day 28 tail-vein glucose was measured by glucometer at 0, 1, 2, 4, 6 and",
      "8 h (Section 4.6.2 lists 0.5 h as well). The response is the percent",
      "reduction relative to the WX390-alone control group (Jin 2025",
      "Equation 5),",
      "and the plasma concentrations come from the separate PK study (n = 6",
      "per group, sparse 3 samples per mouse; Section 4.2.3). The model was",
      "fitted to group-mean data in Phoenix 8.3.5; no between-animal variability",
      "or residual error was reported, so this file is typical-value only."
    )
  )

  ini({
    le0 <- log(12); label("Percent reduction in fasting plasma glucose at zero dorzagliatin concentration (%)") # Jin 2025 Table 5, row 'Dorz + WX390': E = 12 + ...
    lemax <- log(56); label("Maximum additional percent reduction in fasting plasma glucose (%)") # Jin 2025 Table 5, row 'Dorz + WX390': ... + 56 x C^gamma / ...
    lec50 <- log(490); label("Dorzagliatin plasma concentration giving half the maximum effect (ng/mL)") # Jin 2025 Results 2.5: EC50 = 490 ng/mL for the Dorz + WX390 group; Table 5 '490^1.35'
    lhill <- log(1.35); label("Hill coefficient (unitless)") # Jin 2025 Results 2.5: gamma = 1.35 for the Dorz + WX390 group; Table 5 exponents 1.35
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Jin 2025 Equation 6: sigmoid Emax on the plasma concentration, with the
    # response the percent reduction in fasting plasma glucose versus the
    # time-matched WX390-alone control group (Equation 5).
    fpg_reduction <- e0 + emax * CP_DORZAGLIATIN_NGML^hill / (ec50^hill + CP_DORZAGLIATIN_NGML^hill)
  })
}
