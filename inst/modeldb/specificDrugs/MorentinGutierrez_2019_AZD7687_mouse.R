MorentinGutierrez_2019_AZD7687_mouse <- function() {
  description <- "Preclinical (mouse, ICR). PK/PD model for AZD7687, a diacylglycerol acyltransferase 1 (DGAT1) inhibitor, and plasma triglyceride (TAG) during an oral lipid tolerance test with a 20% soybean-oil emulsion (Morentin Gutierrez 2019). One-compartment first-order-absorption PK of the FREE (unbound) plasma AZD7687 concentration per kg body weight (ka fixed); the soybean-oil fat enters a lagged first-order TAG depot whose input to plasma is inhibited by AZD7687 (1 - C/(C + IC50)), on top of an endogenous TAG turnover (baseline R0, loss rate kout)."
  reference <- "Morentin Gutierrez P, Yates J, Nilsson C, Birtles S. Evolving data analysis of an Oral Lipid Tolerance Test toward the standard for the Oral Glucose Tolerance Test: Cross species modeling effects of AZD7687 on plasma triacylglycerol. Pharmacol Res Perspect. 2019;7(2):e00465. doi:10.1002/prp2.465"
  vignette <- "MorentinGutierrez_2019_AZD7687"

  # `tg` (plasma triglyceride) and `depot_tg` (gavaged fat in the gut) are not
  # yet canonical compartments; `tg` is declared paper-specific in the same way
  # by Ahn_2023_vupanorsen.R and Ousey_2026_plozasiran.R.
  paper_specific_compartments <- c("depot_tg", "tg")

  # The lipid gavage is given into `depot_tg`, not into the drug depot.
  dosing <- c("depot", "depot_tg")

  units <- list(
    time = "h",
    dosing = "umol/kg (AZD7687 into depot; convert mg/kg with the molecular weight 367.4 g/mol) and mg per animal (soybean-oil triglyceride into depot_tg; 10 mL/kg x 20% oil x 0.917 g/mL x body weight)",
    concentration = "umol/L (free, unbound AZD7687, Cc) and mg/dL (plasma triglyceride, tg)"
  )

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "AZD7687",
      units = "umol/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "AZD7687",
      units = "umol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_tg = list(
      analyte = "triglycerides",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    tg = list(
      analyte = "triglycerides",
      units = "mg/dL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (ICR)",
    n_subjects = 41,
    n_studies = 1,
    weight_range = "approximately 25 g",
    sex_female_pct = 0,
    disease_state = "Healthy male mice given an oral lipid tolerance test",
    dose_range = "AZD7687 0.1, 1 or 3 mg/kg by oral gavage (n = 11 per group) or vehicle (n = 8)",
    meal_challenge = "Intralipos 20% (soybean-oil emulsion) 10 mL/kg by oral gavage 0.5 h after the AZD7687 dose",
    notes = "ICR mouse OLTT (Methods, Data: Mouse). PK modelled as free concentrations, fu = 0.1; 3 mice per active group were serially bled for PK at 1, 2, 3 and 4 h post dose, TAG sampled at 0.5, 1.5, 2.5, 3.5 and 4.5 h."
  )

  ini({
    # AZD7687 PK (free concentration, per kg) -- Supplementary Table 3
    lka <- fixed(log(3)); label("Absorption rate constant ka (1/h)") # Suppl Table 3 tvKa = 3 1/hr FIXED
    lcl <- log(2.17); label("Apparent clearance of free AZD7687 CL (L/h/kg)") # Suppl Table 3 tvCl = 2.17 L/(kg*hr)
    lvc <- log(13.14); label("Apparent volume of free AZD7687 V (L/kg)") # Suppl Table 3 tvV = 13.14 L/kg

    # Plasma TAG PK/PD -- main Table 1 and Supplementary Table 6
    ltlag_tg <- log(0.98); label("Lag time of soybean-oil TAG absorption TlagTAG (h)") # Suppl Table 6 tvTlagTAG = 0.98 hr
    lka_tg <- log(2.10); label("First-order absorption rate of soybean-oil TAG kaTAG (1/h)") # Suppl Table 6 tvKaTAG = 2.10 1/hr
    lvc_tg <- log(0.012); label("Volume of distribution of exogenous TAG VTAG (dL)") # Suppl Table 6 tvVTAG = 0.012 dL
    lrbase_tg <- log(178.07); label("Baseline plasma TAG R0 (mg/dL)") # Suppl Table 6 tvR0 = 178.07 mg/dL (main Table 1: 178.0)
    lkout_tg <- log(4.41); label("Fractional turnover rate of plasma TAG kout (1/h)") # Suppl Table 6 tvkout = 4.41 1/hr
    lic50_tg <- log(0.081); label("Free AZD7687 concentration giving 50% inhibition of meal TAG input IC50 (umol/L)") # main Table 1 IC50 = 0.081 umol/L (Suppl Table 6 prints 0.08)

    # IIV: log-normal; reported as %CV, converted with omega^2 = log(CV^2 + 1)
    etalcl + etalvc ~ c(
      0.12123,
      0.039666, 0.073928
    ) # Suppl Table 3 IIV CL 35.9%, IIV V 27.7%, Cl-V correlation 0.419
    etaltlag_tg ~ 0.00000361 # Suppl Table 6 IIV lagTAG 0.19%
    etalka_tg ~ 0.11744 # Suppl Table 6 IIV KaTAG 35.3%
    etalrbase_tg ~ 0.030848 # Suppl Table 6 IIV R0 17.7%
    etalic50_tg ~ 0.32527 # Suppl Table 6 IIV IC50 62.0%

    # Residual error
    propSd <- 0.112; label("Proportional residual error, free AZD7687 (fraction)") # Suppl Table 3 Conc Residual Multiplicative = 0.112
    propSd_tg <- 0.198; label("Proportional residual error, plasma TAG (fraction)") # Suppl Table 6 TAG Residual Multiplicative = 0.198
  })

  model({
    # AZD7687 PK per kg body weight (Suppl section S1, mouse)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Plasma TAG (Suppl section S2, mouse; same structure as human)
    tlag_tg <- exp(ltlag_tg + etaltlag_tg)
    ka_tg <- exp(lka_tg + etalka_tg)
    vc_tg <- exp(lvc_tg)
    rbase_tg <- exp(lrbase_tg + etalrbase_tg)
    kout_tg <- exp(lkout_tg)
    ic50_tg <- exp(lic50_tg + etalic50_tg)

    # Eq 1: fraction of meal TAG input not blocked by AZD7687
    inh_tg <- 1 - Cc / (Cc + ic50_tg)

    # Eqs 2-4; the lag acts on the meal depot
    d/dt(depot_tg) <- -ka_tg * depot_tg
    alag(depot_tg) <- tlag_tg
    d/dt(tg) <- rbase_tg * kout_tg + ka_tg * depot_tg / vc_tg * inh_tg - kout_tg * tg
    tg(0) <- rbase_tg

    Cc ~ prop(propSd)
    tg ~ prop(propSd_tg)
  })
}
