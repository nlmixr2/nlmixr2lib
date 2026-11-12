Santos_2024_amphotericinB_rat <- function() {
  description <- "Preclinical (rat). Two-compartment population PK model for intravenous amphotericin B deoxycholate in healthy and Candida albicans-infected Wistar rats, fitted jointly to total plasma and microdialysis free kidney concentrations, with the peripheral compartment representing the kidney. The ODEs are the paper's printed Equations 5-6, which carry a 0.56 factor on the kidney inflow only; with the published estimates this structure does not conserve mass and does not reproduce the paper's figures or NCA AUCs (see the vignette) (Santos 2024)"
  reference <- "Santos VV, Pereira LC, de Araujo JMS, Borges MAH, Brandao CM, Santos LO, Villarreal CF, Azeredo FJ. Amphotericin B tissue penetration and pharmacokinetics in healthy and Candida albicans-infected rats: insights from microdialysis and population modeling. Front Pharmacol. 2024;15:1515462. doi:10.3389/fphar.2024.1515462"
  vignette <- "Santos_2024_amphotericinB_rat"
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  # Santos 2024 Results, Eqs. 5-8: A1 is the drug amount in the central
  # (plasma) compartment and A2 the amount in the kidney compartment; Cplasma =
  # A1/V1 is the total plasma concentration and Ckidney = A2/V2 the free
  # concentration measured by microdialysis in the kidney cortex. Volumes are
  # body-weight-normalised (mL/kg), so amounts are per kg of body weight.
  compartmentData <- list(
    central = list(analyte = "amphotericin B (total)", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amphotericin B (unbound)", units = "ug/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  # Infection was tested as a categorical covariate and not retained.
  covariatesDataExcluded <- list(
    DIS_INFECT_ACTIVE = list(
      description = "Experimental systemic Candida albicans infection, 1 = infected rat, 0 = healthy rat",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy)",
      notes = "Experimental, time-fixed per rat: infected animals received 0.1 mL of 2.5 x 10^6 CFU/mL C. albicans ATCC 10231 into the lateral tail vein 2 days before the PK experiment (Santos 2024 Methods 2.2). Tested on the PK parameters by likelihood-ratio and Wald tests (P < 0.05) and not retained: 'The categorical covariate infection by C. albicans did not influence the model performance (p > 0.05)' (Results). The non-compartmental AUCs did not differ significantly between groups either.",
      source_name = "infection by C. albicans"
    )
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 12L,
    n_studies = 1L,
    weight_range = "250-350 g",
    sex_female_pct = 0,
    disease_state = "Healthy (n = 6) and systemically Candida albicans-infected (n = 6) male Wistar rats; infection by intravenous inoculation of C. albicans ATCC 10231 two days before dosing. Rats were anaesthetised with urethane 1.25 g/kg IP with the kidney surgically exposed for the microdialysis probe.",
    dose_range = "Single 2.5 mg/kg intravenous bolus of amphotericin B deoxycholate (Anforicin B, Cristalia) prepared in water for injection.",
    regions = "Brazil",
    notes = "266 observations (116 total plasma, 150 free kidney) from 12 rats (6 healthy, 6 infected) fitted simultaneously in Monolix 2024R1 (Santos 2024 Methods 2.5). Plasma sampled at 0.083, 0.25, 0.5, 1, 1.5, 2, 4, 8, 12, 24 and 36 h; kidney-cortex microdialysate (CMA/20 probe, 1.5 uL/min Ringer's) collected over 20 h in 1-h intervals and corrected by the in vivo retrodialysis recovery of 10.9% +/- 3.7%. HPLC-DAD assay range 0.25-10 ug/mL. In vitro plasma protein binding was concentration-dependent: 56.4% (2 ug/mL) and 25.5% (8 ug/mL) in healthy rat plasma, 60.9% and 41.8% in infected rat plasma."
  )

  ini({
    # Santos 2024 Table 1, final model. Volumes are mL/kg and the dose is
    # entered in ug/kg, so A/V is in ug/mL.
    lvc <- log(3.3)
    label("Central (plasma) volume V1 (mL/kg)") # Table 1, 'V1 (mL/kg)' = 3.3 (RSE 46.4%)
    lvp <- log(328.58)
    label("Kidney compartment volume V2 (mL/kg)") # Table 1, 'V2 (mL/kg)' = 328.58 (RSE 21.6%)
    lkel <- log(8.53)
    label("Elimination rate constant from the central compartment k10 (1/h)") # Table 1, 'k10 (h-1)' = 8.53 (RSE 33.4%)
    lk12 <- log(91.58)
    label("Plasma-to-kidney distribution rate constant k12 (1/h)") # Table 1, 'k12 (h-1)' = 91.58 (RSE 31.3%)
    lk21 <- log(1.07)
    label("Kidney-to-plasma distribution rate constant k21 (1/h)") # Table 1, 'k21 (h-1)' = 1.07 (RSE 15.1%)

    # IIV (Table 1). Printed under an 'omega2' header, but the printed CVs only
    # match if the values are Monolix omega SDs: sqrt(exp(0.87^2) - 1) = 106%
    # against the printed 106.18%, whereas reading 0.87 as a variance gives
    # 118%. The variances below are the squared SDs. k21 carries no IIV.
    etalvc ~ 0.87^2 # Table 1, 'omega2 (V1)' = 0.87 read as SD (CV 106.18%)
    etalvp ~ 0.67^2 # Table 1, 'omega2 (V2)' = 0.67 read as SD (CV 74.71%)
    etalkel ~ 0.54^2 # Table 1, 'omega2 (k10)' = 0.54 read as SD (CV 58.74%)
    etalk12 ~ 0.57^2 # Table 1, 'omega2 (k12)' = 0.57 read as SD (CV 62.34%)

    # Proportional residual error per output (Table 1, 'Error model
    # parameters'). Printed as 'sigma2' but encoded as the Monolix
    # proportional coefficient b (an SD on the proportional scale), the same
    # reading as the IIV block.
    propSd <- 0.76
    label("Proportional residual error, total plasma concentration (fraction)") # Table 1, 'proportional error, plasma' = 0.76 (RSE 11.3%)
    propSd_Ckidney <- 0.65
    label("Proportional residual error, free kidney concentration (fraction)") # Table 1, 'proportional error, kidney' = 0.65 (RSE 8.26%)
  })
  model({
    # Exponential IIV, Santos 2024 Eq. 4.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)

    # Santos 2024 Eqs. 5-6, transcribed as printed. The 0.56 factor scales
    # only the kidney's gain from plasma, so 0.44 * k12 * central leaves the
    # system as a second elimination route; the paper does not explain the
    # factor (the measured plasma protein binding at 2 ug/mL was 56.4%).
    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- 0.56 * k12 * central - k21 * peripheral1

    # Eqs. 7-8: total plasma and free kidney concentrations (ug/mL).
    Cc <- central / vc
    Ckidney <- peripheral1 / vp

    Cc ~ prop(propSd)
    Ckidney ~ prop(propSd_Ckidney)
  })
}
