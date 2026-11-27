Tapani_2014_nicotinicAcid_rat_obese <- function() {
  description <- paste(
    "Preclinical (rat). Separate-analysis PK/PD feedback model for",
    "nicotinic acid (NiAc) and plasma non-esterified fatty acids (NEFA) in",
    "male obese Zucker (fa/fa) rats after IV constant-rate NiAc infusions;",
    "the obese-rat comparator of the Tapani 2014 joint normal-obese model.",
    "NiAc PK (fixed from earlier analyses) is a one-compartment model with",
    "one Michaelis-Menten elimination and endogenous synthesis. NEFA",
    "turnover: NiAc inhibits NEFA formation (Hill-Imax, Imax = 1), with an",
    "eight-compartment moderator feedback chain (moderator1 inhibits",
    "formation, moderator8 stimulates loss) and a NiAc-independent formation",
    "term kcap that is estimated at practically zero."
  )
  reference <- paste(
    "Tapani S, Almquist J, Leander J, Ahlstrom C, Peletier LA, Jirstrand M,",
    "Gabrielsson J. Joint feedback analysis modeling of nonesterified fatty",
    "acids in obese Zucker rats and normal Sprague-Dawley rats after",
    "different routes of administration of nicotinic acid.",
    "J Pharm Sci. 2014;103(8):2571-2584. doi:10.1002/jps.24077"
  )
  vignette <- "Tapani_2014_nicotinicAcid_rat"
  units <- list(
    time = "min",
    dosing = "umol/kg",
    concentration = "umol/L (NiAc); mmol/L (NEFA)"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "nicotinic acid", units = "umol/kg", specimen = "plasma", verified = TRUE),
    nefa = list(analyte = "non-esterified fatty acids", units = "mmol/L", specimen = "plasma", verified = TRUE),
    moderator1 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator2 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator3 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator4 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator5 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator6 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator7 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator8 = list(
      analyte = "NEFA moderator (feedback signal)",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (male obese Zucker fa/fa)",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "16 weeks at experiment",
    weight_range = "473-547 g",
    sex_female_pct = 0,
    disease_state = paste(
      "Obese, insulin-resistant Zucker (fa/fa) rats; fasted 14 h before",
      "and throughout the experiment."
    ),
    dose_range = paste(
      "NiAc 20 umol/kg IV over 30 min (n = 8) or 51 umol/kg IV over 300 min",
      "(n = 7); vehicle controls (n = 2 + 2)."
    ),
    regions = "Sweden (AstraZeneca R&D Molndal)",
    notes = paste(
      "Obese Zucker subset of the pooled Tapani 2014 data set; n_subjects",
      "is the sum of the group sizes printed in Experimental Design.",
      "11-15 arterial samples per rat."
    )
  )

  ini({
    # ---- NiAc PK, obese Zucker rats (Table 1, fixed from Refs 4-5) ----
    lvmax <- fixed(log(1.59)); label("Vmax of NiAc elimination (umol/min/kg)") # Table 1 obese: Vmax1 = 1.59 (RSE 13.9%)
    lkm <- fixed(log(18.9)); label("Km of NiAc elimination (umol/L)") # Table 1 obese: Km1 = 18.9 (RSE 21.5%)
    lvc <- fixed(log(0.323)); label("Central volume Vc (L/kg)") # Table 1 obese: Vc = 0.323 (RSE 12.4%)
    lsynt <- fixed(log(0.00280)); label("Endogenous NiAc synthesis rate Synt (umol/min/kg)") # Table 1 obese: Synt = 0.00280 (RSE 10.1%)

    # ---- NEFA feedback model, separate analysis (Table 4 'Separate Analysis', obese rats) ----
    lrbase <- log(1.04); label("Baseline NEFA R0 (mmol/L)") # Table 4 separate obese: R0 = 1.04 (RSE 6.31%)
    lkout <- log(0.0683); label("Fractional NEFA turnover kout (L/mmol/min)") # Table 4 separate obese: kout = 0.0683 (RSE 24.0%)
    lktol <- log(0.0697); label("Moderator turnover rate ktol (1/min)") # Table 4 separate obese: ktol = 0.0697 (RSE 9.87%)
    lkcap <- log(4.94e-8); label("NiAc-independent NEFA formation kcap (mmol/L/min)") # Table 4 separate obese: kcap = 4.94 x 10^-8 (RSE not reported)
    lpmod <- log(3.19); label("Amplification factor p on moderator1 (unitless)") # Table 4 separate obese: p = 3.19 (RSE 15.3%)
    lic50 <- log(0.0433); label("NiAc IC50 of NEFA-formation inhibition (umol/L)") # Table 4 separate obese: IC50 = 0.0433 (RSE 49.6%)
    lgam <- log(0.563); label("Sigmoidicity gamma of NiAc inhibition (unitless)") # Table 4 separate obese: gamma = 0.563 (RSE 8.99%)
    imax <- fixed(1); label("Maximal NiAc inhibition of NEFA formation Imax (fraction)") # Table 4: Imax = 1 (fixed, Initial Parameter Estimates)

    # ---- IIV: the paper's IIV % is taken as 100 * sqrt(omega) (CV^2 below); see vignette ----
    etalvmax ~ fixed(0.214^2) # Table 1 obese: IIV Vmax1 21.4%
    etalsynt ~ fixed(0.953^2) # Table 1 obese: IIV Synt 95.3%
    etalrbase ~ 0.158^2 # Table 4 separate obese: IIV R0 15.8%
    etalkout ~ 0.297^2 # Table 4 separate obese: IIV kout 29.7%
    etalic50 ~ 0.747^2 # Table 4 separate obese: IIV IC50 74.7%

    # ---- Residual error ----
    propSd <- fixed(0.400); label("Proportional residual SD on NiAc (fraction)") # Table 1 obese: sigma1 = 40.0% (RSE 26.3%)
    addSd_NEFA <- sqrt(0.0111); label("Additive residual SD on NEFA (mmol/L)") # Table 4 separate obese: sigma2 = 0.0111 (RSE 2.56%), read as a variance
  })

  model({
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)
    vc <- exp(lvc)
    synt <- exp(lsynt + etalsynt)

    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    ktol <- exp(lktol)
    kcap <- exp(lkcap)
    pmod <- exp(lpmod)
    ic50 <- exp(lic50 + etalic50)
    gam <- exp(lgam)

    # Eq. 10: kin from the baseline steady state R = M1 = M8 = R0.
    kin <- (kout * rbase^2 - kcap) * rbase^pmod

    Cp <- central / vc
    # Eq. 4 written as an amount per kg (Vc * dCp/dt); IV infusions dose central.
    d/dt(central) <- synt - vmax * Cp / (km + Cp)
    # Endogenous NiAc steady state of Eq. 4 with no dose.
    central(0) <- synt * km / (vmax - synt) * vc

    # Eq. 6: fraction of NEFA formation remaining, driven by total NiAc.
    icp <- 1 - imax * Cp^gam / (ic50^gam + Cp^gam)

    # Eqs. 5 and 7-9.
    d/dt(nefa) <- kin * icp / moderator1^pmod + kcap - kout * nefa * moderator8
    d/dt(moderator1) <- ktol * (nefa - moderator1)
    d/dt(moderator2) <- ktol * (moderator1 - moderator2)
    d/dt(moderator3) <- ktol * (moderator2 - moderator3)
    d/dt(moderator4) <- ktol * (moderator3 - moderator4)
    d/dt(moderator5) <- ktol * (moderator4 - moderator5)
    d/dt(moderator6) <- ktol * (moderator5 - moderator6)
    d/dt(moderator7) <- ktol * (moderator6 - moderator7)
    d/dt(moderator8) <- ktol * (moderator7 - moderator8)

    # Mi,0 = R0 (text after Eq. 9).
    nefa(0) <- rbase
    moderator1(0) <- rbase
    moderator2(0) <- rbase
    moderator3(0) <- rbase
    moderator4(0) <- rbase
    moderator5(0) <- rbase
    moderator6(0) <- rbase
    moderator7(0) <- rbase
    moderator8(0) <- rbase

    Cc <- Cp
    NEFA <- nefa
    Cc ~ prop(propSd)
    NEFA ~ add(addSd_NEFA)
  })
}
