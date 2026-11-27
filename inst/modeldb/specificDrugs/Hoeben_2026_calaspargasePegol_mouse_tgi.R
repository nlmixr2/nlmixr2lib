Hoeben_2026_calaspargasePegol_mouse_tgi <- function() {
  description <- paste(
    "Preclinical (mouse, female BALB/c nude bearing subcutaneous JHH-5",
    "hepatocellular-carcinoma xenografts). Sequential PK-PD tumor-growth-inhibition",
    "model of calaspargase pegol: a one-compartment plasma asparaginase activity",
    "(PAA) model with Michaelis-Menten elimination scaled by a relative",
    "proteolytic-enzyme pool whose production PAA inhibits (fixed from an earlier",
    "mouse PK analysis) drives a Simeoni 2004 tumor-growth-inhibition model with",
    "exponential-then-linear growth, a linear PAA kill term and three transit",
    "compartments of dying cells."
  )
  reference <- paste(
    "Hoeben E, Madelain V, Acharya C, Barbier Saint Hilaire P, Tyagi E, Okour M,",
    "Fouliard S. PKPD-Based Translational Modeling of Calaspargase Pegol",
    "Preclinical Activity in Hepatocellular Carcinoma. Eur J Drug Metab",
    "Pharmacokinet. 2026. doi:10.1007/s13318-026-01010-4.",
    "Structural TGI formulation from Simeoni M, Magni P, Cammia C, De Nicolao G,",
    "Croci V, Pesenti E, et al. Predictive pharmacokinetic-pharmacodynamic",
    "modeling of tumor growth kinetics in xenograft models after administration",
    "of anticancer agents. Cancer Res. 2004;64(3):1094-1101.",
    "doi:10.1158/0008-5472.CAN-03-2524.",
    sep = " "
  )
  vignette <- "Hoeben_2026_calaspargasePegol"

  units <- list(time = "h", dosing = "IU/kg", concentration = "IU/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "calaspargase pegol (asparaginase activity)",
      units = "IU/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    enzyme = list(
      analyte = "proteolytic enzyme eliminating calaspargase pegol (relative amount, 1 at baseline)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "non-replicating dying tumor cells (transit state 1)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "non-replicating dying tumor cells (transit state 2)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "non-replicating dying tumor cells (transit state 3)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female BALB/c nude, JHH-5 hepatocellular-carcinoma xenograft)",
    n_subjects = 40L,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = paste(
      "subcutaneous JHH-5 hepatocellular-carcinoma xenograft (ASNS",
      "promoter-hypermethylated), mean tumor volume 175 mm^3 at the start of",
      "treatment"
    ),
    dose_range = "vehicle, 20, 40 and 100 IU/kg calaspargase pegol IV every 4 days for 28 days (8 doses)",
    regions = "preclinical in-vivo xenograft study",
    notes = paste(
      "10 mice per group (Methods 'Data Collection'; Supplementary Table S1);",
      "527 tumor volume and 87 PAA observations. Tumor volumes below 1 mm^3",
      "(36% of observations) were imputed at the 1 mm^3 limit for the final fit.",
      "The PK sub-model (and its IIV) was estimated earlier on 344 PAA samples",
      "from 90 non-grafted and grafted female BALB/c nude mice dosed 6-300 IU/kg",
      "(Supplementary Table S2) and was held fixed in this PK-PD fit."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Mouse PK of plasma asparaginase activity (PAA), fixed from the earlier
    # mouse PK analysis (Methods 'PK Model and Prediction of PAA'; Table 1
    # 'Mice PK Model'). Supplementary Figure S1B gives the equations.
    # ---------------------------------------------------------------------
    lvc <- fixed(log(0.0636)); label("Central volume of distribution V (L/kg)") # Table 1 Mice PK: V = 0.0636 L/kg
    lvmax <- fixed(log(5.07)); label("Maximum elimination capacity Vm (IU/kg/h)") # Table 1 Mice PK: Vm = 5.07 IU/kg/h
    lkm <- fixed(log(4000)); label("PAA at half-maximum elimination capacity Km (IU/L)") # Table 1 Mice PK: Km = 4000 IU/L
    lkout <- fixed(log(0.000813)); label("First-order loss rate of the eliminating enzyme kout (1/h)") # Table 1 Mice PK: kout = 0.000813 1/h
    lec50 <- fixed(log(107)); label("PAA reducing the enzyme production to 50% EC50 (IU/L)") # Table 1 Mice PK: EC50 = 107 IU/L

    # ---------------------------------------------------------------------
    # Simeoni tumor-growth-inhibition parameters (Table 1 'Mice PKPD TGI
    # Model', final model with imputation of values below 1 mm^3).
    # ---------------------------------------------------------------------
    ltumorExpGrowth <- log(0.00708); label("Exponential tumor growth rate kg (1/h)") # Table 1 Mice PKPD TGI: kg = 0.00708 1/h
    ltumorLinGrowth <- log(70.0); label("Linear tumor growth rate kgl (mm^3/h)") # Table 1 Mice PKPD TGI: kgl = 70.0 mm3/h
    ldrugSlope <- log(3.22e-5); label("Slope of the linear PAA kill effect slp (L/IU/h)") # Table 1 Mice PKPD TGI: slp = 3.22 x 10^-5 L/IU/h
    ldamageTransit <- log(1 / 41.1); label("Transit rate through the dying-cell stages 1/tau (1/h)") # Table 1 Mice PKPD TGI: tau = 41.1 h
    lrbase_tumor <- log(152); label("Tumor volume at time zero W0 (mm^3)") # Table 1 Mice PKPD TGI: W0 = 152 mm3
    psi <- fixed(20); label("Switch sharpness between exponential and linear growth psi (unitless)") # Table 1 Mice PKPD TGI: psi = 20 FIX

    # IIV. Table 1 reports IIV as SD of the log-normal random effect, so the
    # variance is SD^2. The PK IIV was fixed with the rest of the PK model.
    etalvc ~ fixed(0.00585225) # Table 1 Mice PK: IIV SD V = 0.0765
    etalvmax ~ fixed(0.026896) # Table 1 Mice PK: IIV SD Vm = 0.164
    etalkm ~ fixed(0.00390625) # Table 1 Mice PK: IIV SD Km = 0.0625
    etaltumorExpGrowth ~ 0.050176 # Table 1 Mice PKPD TGI: IIV SD kg = 0.224
    etaldrugSlope ~ 0.128881 # Table 1 Mice PKPD TGI: IIV SD slp = 0.359
    etalrbase_tumor ~ 0.012996 # Table 1 Mice PKPD TGI: IIV SD W0 = 0.114

    # Residual error. The PAA model's proportional residual error is named in
    # Methods but its magnitude is not reported, so it is carried as fixed(0).
    propSd <- fixed(0); label("Proportional residual error on PAA (fraction; ZERO -- not reported in the source)")
    addSd_tumor_vol <- 0.827; label("Additive residual error on tumor volume (mm^3)") # Table 1: RUV1 SD additive component = 0.827
    propSd_tumor_vol <- 0.437; label("Proportional residual error on tumor volume (fraction)") # Table 1: RUV2 SD proportional component = 0.437
  })

  model({
    vc <- exp(lvc + etalvc)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)
    kout <- exp(lkout)
    ec50 <- exp(lec50)
    # Supplementary Figure S1B: Rin = Kout, so the enzyme pool sits at 1 in the
    # absence of drug and can only fall below 1 while PAA is present.
    kin <- kout

    tumorExpGrowth <- exp(ltumorExpGrowth + etaltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    drugSlope <- exp(ldrugSlope + etaldrugSlope)
    damageTransit <- exp(ldamageTransit)
    rbase_tumor <- exp(lrbase_tumor + etalrbase_tumor)

    # PAA (IU/L) = amount (IU/kg) / V (L/kg)
    Cc <- central / vc

    # Supplementary Figure S1B: dAc/dt = -Vmax*Ac*Ae/(Km*V + Ac), which is
    # Vmax*Ae*Cc/(Km + Cc) after dividing numerator and denominator by V.
    d/dt(central) <- -vmax * enzyme * Cc / (km + Cc)
    d/dt(enzyme) <- kin * (1 - Cc / (Cc + ec50)) - kout * enzyme
    enzyme(0) <- 1

    # Hoeben 2026 Equations 1-2. The growth term uses the replicating cells
    # (TS) in both the numerator and the exponential-to-linear switch, as
    # printed; the measured tumor volume is the sum of the replicating and the
    # three dying-cell stages (Simeoni 2004).
    drugEffectCyclingCells <- drugSlope * Cc
    d/dt(cycling_cells) <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth / tumorLinGrowth * cycling_cells)^psi)^(1 / psi) -
      drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)
    cycling_cells(0) <- rbase_tumor

    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    Cc ~ prop(propSd)
    tumor_vol ~ add(addSd_tumor_vol) + prop(propSd_tumor_vol)
  })
}
