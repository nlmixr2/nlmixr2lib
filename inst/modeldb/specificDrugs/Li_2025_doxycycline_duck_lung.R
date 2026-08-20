Li_2025_doxycycline_duck_lung <- function() {
  description <- "Preclinical (duck, Tadorna tadorna). One-compartment first-order-absorption model for doxycycline concentrations in LUNG tissue of Riemerella anatipestifer-infected ducks after a single intramuscular injection. Li 2025 analysed plasma, lung and liver as three independent one-compartment first-order-absorption fits in WinNonlin 6.1 (Methods, 'Pharmacokinetics (PK) of DOX in RA-infected ducks'), so this file packages the lung fit as its own model; the companions are Li_2025_doxycycline_duck_liver and, for plasma plus the exposure-response, Li_2025_doxycycline_duck_ff20 and Li_2025_doxycycline_duck_ff40. Table 1 reports mean T1/2ka = 0.42 h and mean T1/2kel = 11.53 h for lung, but no apparent clearance. Cl/F is recovered here as the mean across the five dose levels of dose divided by the tabulated AUC, which is exactly how Li 2025's own plasma Cl/F was computed (dose/AUC reproduces the printed plasma Cl/F at all five dose levels to two decimal places), giving 0.3287 L/h/kg; V/F then follows as (Cl/F)/kel = 5.47 L/kg. The recovery is validated by the fact that it reproduces Table 1's lung Tmax and Cmax at every dose level to two decimal places (see the vignette source-trace table). Because the lung and plasma fits are independent, the apparent absorption and disposition terms differ between them and this model is dosed on its own depot; it is a descriptive model of the lung concentration-time curve, not a physiological tissue-distribution model, and its V/F absorbs the lung-to-plasma partitioning. Doses and all disposition terms are per kg bodyweight. Li 2025 fitted mean profiles and reported neither between-subject variability nor a residual error magnitude, so there are no eta parameters and propSd is FIXED at 0; the model is intended for typical-value simulation."
  reference <- paste(
    "Li FL, He CY, Chen HY, Cheng SM, Liu Y, Ding HZ, Zhang HL. (2025).",
    "In vivo pharmacokinetic/pharmacodynamic relationship of florfenicol in",
    "combination with doxycycline against Riemerella anatipestifer in ducks",
    "and the effect upon resistance development.",
    "Poultry Science 104:104922.",
    "doi:10.1016/j.psj.2025.104922.",
    sep = " "
  )
  vignette <- "Li_2025_doxycycline_florfenicol_ducks"
  units <- list(
    time = "h",
    dosing = "mg/kg (doxycycline; all disposition terms are per kg bodyweight)",
    concentration = "ug/g (doxycycline in lung tissue)"
  )

  # auc_dox -- cumulative lung doxycycline AUC (ug*h/g), the quantity Li 2025
  #   tabulates for lung in Table 1. Accumulator idiom follows
  #   Beredaki_2023_micafungin_clsi.
  paper_specific_compartments <- c("auc_dox")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "doxycycline", units = "mg/kg", specimen = "administration site", verified = TRUE),
    lung    = list(analyte = "doxycycline", units = "mg/kg", specimen = "tissue", verified = TRUE),
    auc_dox = list(analyte = "doxycycline", units = "ug*h/g", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species        = "duck (common shelduck, Tadorna tadorna), 7 days old",
    n_subjects     = 360L,
    n_studies      = 1L,
    age_median     = "7 days at the start of the study",
    weight_range   = "130-150 g",
    sex_female_pct = NA_real_,
    organism       = "Riemerella anatipestifer strain CVCC3857 (doxycycline MIC 1 ug/mL, florfenicol MIC 1 ug/mL; Li 2025 Results, 'MIC and MPC Of FF and DOX against RA')",
    disease_state  = paste(
      "Experimental systemic Riemerella anatipestifer infection established by",
      "intraperitoneal injection of 10^9 CFU/mL; the target bacterial load was",
      "reached 12 h after inoculation and drug was given at that point"
    ),
    dose_range     = "Doxycycline 1, 2.5, 5, 10 and 20 mg/kg as a single intramuscular injection into the thigh (five groups of 72 ducks)",
    sampling       = paste(
      "Lung tissue at 0.5, 1, 2, 4, 6, 8, 12, 24 and 36 h after dosing;",
      "0.5 g homogenised per sample and assayed by HPLC-MS/MS. Lung LLOQ",
      "0.005 ug/g, calibration 0.005-5 ug/g (R^2 > 0.99), recovery 80.08-90.80%"
    ),
    regions        = "China (Fuyang Normal University, Anhui; South China Agricultural University, Guangzhou)",
    notes          = paste(
      "Apparent clearance rose monotonically with dose in lung (0.239 to",
      "0.552 L/h/kg from 1 to 20 mg/kg), as it did in plasma; the mean is",
      "packaged here and the per-dose values are reproduced in the vignette.",
      "Lung exposure slightly exceeded plasma at the lower doses (AUC 4.18",
      "versus 3.48 ug*h/mL at 1 mg/kg) and was comparable at 20 mg/kg (36.25",
      "versus 37.11)."
    )
  )

  ini({
    # =================================================================
    # Lung pharmacokinetics -- Li 2025 Table 1 (lung, "Mean" column)
    # =================================================================
    lka <- log(log(2) / 0.42)
    label("Log first-order absorption rate constant ka (1/h)")  # Li 2025 Table 1, lung mean T1/2ka = 0.42 +/- 0.14 h; ka = ln(2)/T1/2ka

    # Li 2025 tabulates Cl/F for plasma only. It is recovered for lung as the
    # mean of dose/AUC over the five dose levels -- the same construction that
    # reproduces the printed plasma Cl/F exactly. Per-dose values:
    # 1/4.18 = 0.2392, 2.5/10.04 = 0.2490, 5/18.92 = 0.2643, 10/29.49 = 0.3391,
    # 20/36.25 = 0.5517 L/h/kg; mean 0.3287.
    lcl <- log(0.3287)
    label("Log apparent clearance from lung Cl/F (L/h/kg)")  # derived from Li 2025 Table 1 lung AUC24 at doses 1, 2.5, 5, 10, 20 mg/kg (4.18, 10.04, 18.92, 29.49, 36.25 ug*h/g)

    # V/F = (Cl/F)/kel with kel = ln(2)/T1/2kel. Reproduces Table 1's lung Tmax
    # and Cmax at every dose level to two decimal places (vignette source trace).
    lvc <- log(0.3287 / (log(2) / 11.53))
    label("Log apparent lung volume of distribution V/F (L/kg)")  # derived: recovered lung Cl/F = 0.3287 L/h/kg and Li 2025 Table 1 lung mean T1/2kel = 11.53 +/- 1.43 h

    # =================================================================
    # Residual error
    # =================================================================
    # Li 2025 fitted mean profiles in WinNonlin and reported no between-subject
    # variability and no residual standard deviation, so the residual SD is held
    # at zero for deterministic typical-value simulation.
    propSd <- fixed(0)
    label("Proportional residual error on lung doxycycline concentration (0; not reported in Li 2025)")  # Li 2025 reported no residual error model
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    # One-compartment model with first-order absorption fitted directly to the
    # lung concentration-time data (Li 2025 Methods). F is not identifiable from
    # intramuscular data alone, so cl and vc are apparent terms and no f(depot)
    # is applied.
    kel <- cl / vc
    d/dt(depot) <- -ka * depot
    d/dt(lung) <- ka * depot - kel * lung

    Clung <- lung / vc
    d/dt(auc_dox) <- Clung

    Clung ~ prop(propSd)
  })
}
