Dahan_2026_naloxone <- function() {
  description <- paste(
    "Mechanistic PK/PD simulation (utility) model of intravenous naloxone",
    "reversing fentanyl-induced respiratory depression versus precipitating",
    "opioid withdrawal in opioid-tolerant adults with opioid use disorder.",
    "One-compartment naloxone PK with a hysteresis effect compartment;",
    "naloxone competitively lowers the effective fentanyl effect-site",
    "concentration (CEF = CE / (1 + CNALE / Ki)). CEF drives an inhibitory",
    "sigmoid-Emax model of isohypercapnic minute ventilation and the",
    "normalized mu-opioid receptor activation R = CEF / (CEF + C50).",
    "Respiratory reversal is R < 0.5; precipitated withdrawal is R falling",
    "below (1 - w) times its pre-naloxone baseline. The fentanyl effect-site",
    "concentration is held constant at the value giving the overdose",
    "severity R0 for a typical subject (CE = R0 / (1 - R0) * C50). The only",
    "random effect is log-normal inter-individual variability in C50;",
    "population probabilities and the utility P(reversal) - k * P(withdrawal)",
    "are the means of the per-subject indicators across a simulated cohort."
  )
  reference <- paste(
    "Dahan A, Dahan JDC, Olofsen E, Van Lemmen M, Van Velzen M, Sarton E,",
    "Dungan G, Miller TL, Raffa RB, Martini C, Niesters M.",
    "Naloxone reversal of fentanyl-induced respiratory depression versus",
    "naloxone-induced precipitated withdrawal: an in silico utility model.",
    "Front Pharmacol. 2026;17:1944044. doi:10.3389/fphar.2026.1944044."
  )
  vignette <- "Dahan_2026_naloxone"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")
  dosing <- "central"

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "naloxone",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "naloxone",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Virtual opioid-tolerant adults with opioid use disorder (OUD) in a",
      "fentanyl overdose treated with intravenous naloxone"
    ),
    dose_range = paste(
      "Steady-state naloxone effect-site concentrations 0.1-300 ng/mL;",
      "intravenous naloxone bolus doses about 0.06-7 mg in the calibration",
      "against clinical withdrawal data"
    ),
    notes = paste(
      "In silico simulation study; no new human data were collected (Dahan",
      "2026 Methods). Monte Carlo cohorts of 40,000 virtual subjects per",
      "scenario, with only C50 varying between subjects. Parameter values are",
      "taken from the literature (Table 1): ventilatory model from Boom 2013,",
      "OUD-specific C50 from Algera 2020, naloxone Ki from Moss 2020 /",
      "Stolbach 2023, naloxone ke0 from van Lemmen 2025b, naloxone volume and",
      "elimination rate constant from the FDA 2022 review. The withdrawal",
      "threshold w was calibrated against emergency-department (Yugar 2023)",
      "and prospective naloxone (van Lemmen 2025b, 2026) dose-withdrawal data",
      "(Figure 5)."
    )
  )

  ini({
    # Naloxone pharmacokinetics (Table 1, 'Naloxone pharmacokinetic and antagonism model parameters')
    lvc <- fixed(log(128))
    label("Naloxone apparent volume of distribution Vd (L)") # Table 1: Vd = 128 L (FDA 2022)
    lkel <- fixed(log(0.0077))
    label("Naloxone elimination rate constant ke (1/min)") # Table 1: ke = 0.0077 min-1 (FDA 2022)
    lke0 <- fixed(log(0.087))
    label("Naloxone blood-effect-site equilibration rate constant ke0N (1/min)") # Table 1: ke0 = 0.087 min-1 (van Lemmen 2025b)
    lkd <- fixed(log(1.5))
    label("Naloxone competitive inhibition constant Ki (ng/mL)") # Table 1: Ki = 1.5 ng/mL (Moss 2020; Stolbach 2023); Equation 1

    # Fentanyl ventilatory pharmacodynamics (Table 1, 'Fentanyl pharmacokinetic and pharmacodynamic model parameters')
    le0 <- fixed(log(19.9))
    label("Baseline isohypercapnic minute ventilation EMAX (L/min)") # Table 1: EMAX = 19.9 L/min (Boom 2013); Equation 2
    lec50 <- fixed(log(4.08))
    label("Fentanyl effect-site concentration causing 50% respiratory depression C50 in OUD (ng/mL)") # Methods text: C50 = 4.08 ng/mL (Table 1 rounds to 4.1; Algera 2020)
    hill <- fixed(1)
    label("Steepness of the ventilation concentration-response gamma (unitless)") # Table 1: gamma = 1 (Boom 2013); Equation 3

    # Overdose severity, withdrawal threshold and utility weight (Table 1, 'Range of withdrawal and utility model parameters')
    rbase <- fixed(0.9)
    label("Baseline (pre-naloxone) normalized mu-opioid receptor activation R0 for a typical subject (fraction)") # Table 1 range 0-0.97; Figure 1 and Figures 3-4 reference value R0 = 0.9
    wd_thres <- 0.3
    label("Withdrawal threshold w: fractional drop in receptor activation that precipitates withdrawal (fraction)") # Abstract and Table 1: calibrated w = 0.30 (95% CI 0.16-0.40); Equation 6
    wd_weight <- fixed(1)
    label("Weight k on withdrawal probability in the utility function (unitless)") # Table 1 range 0.25-4; Equation 7, k = 1 is the reference value

    # Inter-individual variability
    etalec50 ~ fixed(0.31) # Table 1: omega2 (C50) = 0.31, log-domain variance (Boom 2013)
  })

  model({
    # Naloxone disposition: one compartment, IV bolus into central
    vc <- exp(lvc)
    kel <- exp(lkel)
    ke0 <- exp(lke0)
    kd <- exp(lkd)

    # Ventilatory and receptor parameters; C50 is the only parameter with IIV
    e0 <- exp(le0)
    ec50 <- exp(lec50 + etalec50)

    # Fentanyl effect-site concentration held at the level giving the
    # overdose severity rbase for the typical C50 (Equation 4 rearranged;
    # Figure 1: CE = R0 / (1 - R0) * C50 with C50 = 4.08 ng/mL)
    ce_fent <- rbase / (1 - rbase) * exp(lec50)

    d/dt(central) <- -kel * central
    # Naloxone plasma concentration: mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc
    # Naloxone effect-site concentration (ng/mL)
    d/dt(effect) <- ke0 * (Cc - effect)

    # Competitive antagonism: effective fentanyl concentration at the receptor (Equation 1)
    ceff_fent <- ce_fent / (1 + effect / kd)

    # Isohypercapnic minute ventilation, inhibitory sigmoid Emax (Equations 2-3)
    vent <- e0 / (1 + (ceff_fent / ec50)^hill)

    # Normalized mu-opioid receptor activation before (Equation 4) and after
    # (Equation 5) naloxone, for this subject's C50
    ract0 <- ce_fent / (ce_fent + ec50)
    ract <- ceff_fent / (ceff_fent + ec50)

    # Model-defined respiratory reversal: residual respiratory depression
    # below 50% of baseline, equivalent to CEF < C50 and to ract < 0.5
    reversal <- 0
    if (ract < 0.5) {
      reversal <- 1
    }
    # Precipitated withdrawal (Equation 6)
    withdrawal <- 0
    if (ract < (1 - wd_thres) * ract0) {
      withdrawal <- 1
    }
    # Per-subject utility; its cohort mean is Equation 7
    utility <- reversal - wd_weight * withdrawal
  })
}
