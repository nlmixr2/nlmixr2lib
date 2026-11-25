MorentinGutierrez_2019_AZD7687_rat <- function() {
  description <- "Preclinical (rat, Han-Wistar). PK/PD model for AZD7687, a diacylglycerol acyltransferase 1 (DGAT1) inhibitor, and plasma triglyceride (TAG) during an oral lipid tolerance test with corn oil (Morentin Gutierrez 2019). One-compartment first-order-absorption PK of the FREE (unbound) plasma AZD7687 concentration per kg body weight (ka fixed); the corn-oil fat enters a lagged first-order TAG depot whose input to plasma is inhibited by AZD7687 (1 - C/(C + IC50)), on top of an endogenous TAG turnover (baseline R0, loss rate kout)."
  reference <- "Morentin Gutierrez P, Yates J, Nilsson C, Birtles S. Evolving data analysis of an Oral Lipid Tolerance Test toward the standard for the Oral Glucose Tolerance Test: Cross species modeling effects of AZD7687 on plasma triacylglycerol. Pharmacol Res Perspect. 2019;7(2):e00465. doi:10.1002/prp2.465"
  vignette <- "MorentinGutierrez_2019_AZD7687"

  # `tg` (plasma triglyceride) and `depot_tg` (gavaged fat in the gut) are not
  # yet canonical compartments; `tg` is declared paper-specific in the same way
  # by Ahn_2023_vupanorsen.R and Ousey_2026_plozasiran.R.
  paper_specific_compartments <- c("depot_tg", "tg")

  # The oil gavage is given into `depot_tg`, not into the drug depot.
  dosing <- c("depot", "depot_tg")

  units <- list(
    time = "h",
    dosing = "umol/kg (AZD7687 into depot; convert mg/kg with the molecular weight 367.4 g/mol) and mg per animal (corn-oil triglyceride into depot_tg; 5 mL/kg x 0.9 g/mL x body weight)",
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
    species = "rat (Han-Wistar)",
    n_subjects = 53,
    n_studies = 1,
    weight_range = "approximately 230 g",
    sex_female_pct = 0,
    disease_state = "Healthy chow-fed male rats fasted 14 h before an oral lipid tolerance test",
    dose_range = "AZD7687 0.1, 0.3 or 3 mg/kg by oral gavage (n = 12 per group), vehicle (n = 11), or naive untreated (n = 6)",
    meal_challenge = "Corn oil 5 mL/kg by oral gavage 2 h after the AZD7687 dose (no oil in the naive group)",
    notes = "Han-Wistar rat OLTT (Methods, Data: Rat). PK modelled as free concentrations, fu = 0.03; 3 rats per active group were serially bled for PK at 2, 3.5, 5, 8 and 11 h post dose, TAG sampled at 2, 2.5, 3.5, 5, 8 and 11 h."
  )

  ini({
    # AZD7687 PK (free concentration, per kg) -- Supplementary Table 2
    lka <- fixed(log(3)); label("Absorption rate constant ka (1/h)") # Suppl Table 2 tvKa = 3 1/hr FIXED
    lcl <- log(0.73); label("Apparent clearance of free AZD7687 CL (L/h/kg)") # Suppl Table 2 tvCl = 0.73 L/(kg*hr)
    lvc <- log(8.66); label("Apparent volume of free AZD7687 V (L/kg)") # Suppl Table 2 tvV = 8.66 L/kg

    # Plasma TAG PK/PD -- main Table 1 and Supplementary Table 5
    ltlag_tg <- log(1.24); label("Lag time of corn-oil TAG absorption TlagTAG (h)") # Suppl Table 5 tvTlagTAG = 1.24 hr
    lka_tg <- log(0.37); label("First-order absorption rate of corn-oil TAG kaTAG (1/h)") # Suppl Table 5 tvKaTAG = 0.37 1/hr
    lvc_tg <- log(2.19); label("Volume of distribution of exogenous TAG VTAG (dL)") # Suppl Table 5 tvVTAG = 2.19 dL
    lrbase_tg <- log(109.2); label("Baseline plasma TAG R0 (mg/dL)") # main Table 1 R0 = 109.2 mg/dL (Suppl Table 5 prints 109.3)
    lkout_tg <- log(0.28); label("Fractional turnover rate of plasma TAG kout (1/h)") # Suppl Table 5 tvkout = 0.28 1/hr
    lic50_tg <- log(0.162); label("Free AZD7687 concentration giving 50% inhibition of meal TAG input IC50 (umol/L)") # Suppl Table 5 tvIC50 = 0.162 umol/L

    # IIV: log-normal; reported as %CV, converted with omega^2 = log(CV^2 + 1)
    etalcl ~ 0.18741 # Suppl Table 2 IIV CL 45.4%
    etalvc ~ 0.013596 # Suppl Table 2 IIV V 11.7%
    etaltlag_tg ~ 0.00205 # Suppl Table 5 IIV lagTAG 4.53%
    etalka_tg ~ 0.68515 # Suppl Table 5 IIV KaTAG 99.2%
    etalvc_tg ~ 0.27192 # Suppl Table 5 IIV VTAG 55.9%
    etalrbase_tg ~ 0.033652 # Suppl Table 5 IIV R0 18.5%

    # Residual error
    addSd <- 0.0023; label("Additive residual error, free AZD7687 (umol/L)") # Suppl Table 2 Conc Residual Additive = 0.0023 uM
    propSd <- 0.15; label("Proportional residual error, free AZD7687 (fraction)") # Suppl Table 2 Conc Residual Multiplicative = 0.15
    propSd_tg <- 0.3; label("Proportional residual error, plasma TAG (fraction)") # Suppl Table 5 TAG Residual Multiplicative = 0.3
  })

  model({
    # AZD7687 PK per kg body weight (Suppl section S1, rat)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Plasma TAG (Suppl section S2, rat; same structure as human)
    tlag_tg <- exp(ltlag_tg + etaltlag_tg)
    ka_tg <- exp(lka_tg + etalka_tg)
    vc_tg <- exp(lvc_tg + etalvc_tg)
    rbase_tg <- exp(lrbase_tg + etalrbase_tg)
    kout_tg <- exp(lkout_tg)
    ic50_tg <- exp(lic50_tg)

    # Eq 1: fraction of meal TAG input not blocked by AZD7687
    inh_tg <- 1 - Cc / (Cc + ic50_tg)

    # Eqs 2-4; the lag acts on the meal depot
    d/dt(depot_tg) <- -ka_tg * depot_tg
    alag(depot_tg) <- tlag_tg
    d/dt(tg) <- rbase_tg * kout_tg + ka_tg * depot_tg / vc_tg * inh_tg - kout_tg * tg
    tg(0) <- rbase_tg

    Cc ~ add(addSd) + prop(propSd)
    tg ~ prop(propSd_tg)
  })
}
