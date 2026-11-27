Tapani_2014_nicotinicAcid_rat_normal <- function() {
  description <- paste(
    "Preclinical (rat). Separate-analysis PK/PD feedback model for",
    "nicotinic acid (NiAc) and plasma non-esterified fatty acids (NEFA) in",
    "normal male Sprague-Dawley rats after IV constant-rate, stepwise IV and",
    "oral NiAc; the normal-rat comparator of the Tapani 2014 joint",
    "normal-obese model. NiAc PK (fixed from earlier analyses) is a",
    "two-compartment model with two parallel Michaelis-Menten eliminations,",
    "endogenous synthesis and parallel linear + saturable gut absorption.",
    "NEFA turnover: NiAc inhibits NEFA formation (Hill-Imax, Imax = 1), with",
    "an eight-compartment moderator feedback chain (moderator1 inhibits",
    "formation, moderator8 stimulates loss) and a NiAc-independent formation",
    "term kcap."
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
    depot = list(analyte = "nicotinic acid", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nicotinic acid", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotinic acid", units = "umol/kg", specimen = "tissue", verified = TRUE),
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
    species = "rat (male Sprague-Dawley)",
    n_subjects = 99L,
    n_studies = 3L,
    age_range = "16 weeks at experiment",
    weight_range = "220-367 g",
    sex_female_pct = 0,
    disease_state = "Healthy; fasted 14 h before and throughout the experiment.",
    dose_range = paste(
      "NiAc 1, 5 or 20 umol/kg IV over 30 min; 5, 10 or 51 umol/kg IV over",
      "300 min; 24.4, 81.2 or 812 umol/kg oral gavage; 5 umol/kg over",
      "30 min followed by stepwise-decreasing infusion (with or without a",
      "second 5 umol/kg 30-min infusion at 210 min); vehicle controls."
    ),
    regions = "Sweden (AstraZeneca R&D Molndal)",
    notes = paste(
      "Sprague-Dawley subset of the pooled Tapani 2014 data set; n_subjects",
      "is the sum of the group sizes printed in Experimental Design.",
      "11-15 arterial samples per rat."
    )
  )

  ini({
    # ---- NiAc PK, normal Sprague-Dawley rats (Table 1, fixed from Refs 4-5) ----
    lvmax_hi <- fixed(log(0.0871)); label("Vmax of high-affinity NiAc elimination (umol/min/kg)") # Table 1 SD: Vmax1 = 0.0871 (RSE 22.8%)
    lkm_hi <- fixed(log(0.235)); label("Km of high-affinity NiAc elimination (umol/L)") # Table 1 SD: Km1 = 0.235 (RSE 29.2%)
    lvmax_lo <- fixed(log(7.09)); label("Vmax of low-affinity NiAc elimination (umol/min/kg)") # Table 1 SD: Vmax2 = 7.09 (RSE 39.6%)
    lkm_lo <- fixed(log(74.5)); label("Km of low-affinity NiAc elimination (umol/L)") # Table 1 SD: Km2 = 74.5 (RSE 43.4%)
    lvc <- fixed(log(0.393)); label("Central volume Vc (L/kg)") # Table 1 SD: Vc = 0.393 (RSE 5.29%)
    lvp <- fixed(log(0.172)); label("Peripheral volume Vt (L/kg)") # Table 1 SD: Vt = 0.172 (RSE 35.2%)
    lq <- fixed(log(0.000852)); label("Intercompartmental clearance Cld (L/min/kg)") # Table 1 SD: Cld = 0.000852 (RSE 27.8%)
    lsynt <- fixed(log(0.00355)); label("Endogenous NiAc synthesis rate Synt (umol/min/kg)") # Table 1 SD: Synt = 0.00355 (RSE 23.3%)
    lka <- fixed(log(0.00477)); label("First-order gut absorption rate ka (1/min)") # Table 1 SD: ka = 0.00477 (RSE 33.5%)
    ltmax_abs <- fixed(log(2.96)); label("Maximum saturable gut absorption rate Vmax,g (umol/min/kg)") # Table 1 SD: Vmax,g = 2.96 (RSE 21.3%)
    lkt_abs <- fixed(log(20.5)); label("Gut amount at half-maximal saturable absorption Km,g (umol/kg)") # Table 1 SD: Km,g = 20.5 (RSE 29.6%)

    # ---- NEFA feedback model, separate analysis (Table 4 'Separate Analysis', SD rats) ----
    lrbase <- log(0.739); label("Baseline NEFA R0 (mmol/L)") # Table 4 separate SD: R0 = 0.739 (RSE 3.50%)
    lkout <- log(0.316); label("Fractional NEFA turnover kout (L/mmol/min)") # Table 4 separate SD: kout = 0.316 (RSE 9.83%)
    lktol <- log(0.0239); label("Moderator turnover rate ktol (1/min)") # Table 4 separate SD: ktol = 0.0239 (RSE 1.92%)
    lkcap <- log(0.0281); label("NiAc-independent NEFA formation kcap (mmol/L/min)") # Table 4 separate SD: kcap = 0.0281 (RSE 9.85%)
    lpmod <- log(1.15); label("Amplification factor p on moderator1 (unitless)") # Table 4 separate SD: p = 1.15 (RSE 2.61%)
    lic50 <- log(0.0757); label("NiAc IC50 of NEFA-formation inhibition (umol/L)") # Table 4 separate SD: IC50 = 0.0757 (RSE 14.9%)
    lgam <- log(2.19); label("Sigmoidicity gamma of NiAc inhibition (unitless)") # Table 4 separate SD: gamma = 2.19 (RSE 4.33%)
    imax <- fixed(1); label("Maximal NiAc inhibition of NEFA formation Imax (fraction)") # Table 4: Imax = 1 (fixed, Initial Parameter Estimates)

    # ---- IIV: the paper's IIV % is taken as 100 * sqrt(omega) (CV^2 below); see vignette ----
    etalvmax_hi ~ fixed(0.927^2) # Table 1 SD: IIV Vmax1 92.7%
    etalvmax_lo ~ fixed(0.291^2) # Table 1 SD: IIV Vmax2 29.1%
    etalsynt ~ fixed(1.09^2) # Table 1 SD: IIV Synt 109%
    etalka ~ fixed(0.108^2) # Table 1 SD: IIV ka 10.8%
    etaltmax_abs ~ fixed(0.108^2) # Table 1 SD: IIV Vmax,g 10.8%
    etalrbase ~ 0.199^2 # Table 4 separate SD: IIV R0 19.9%
    etalkout ~ 0.474^2 # Table 4 separate SD: IIV kout 47.4%
    etalic50 ~ 1.34^2 # Table 4 separate SD: IIV IC50 134%

    # ---- Residual error ----
    propSd <- fixed(0.428); label("Proportional residual SD on NiAc (fraction)") # Table 1 SD: sigma1 = 42.8% (RSE 5.16%)
    addSd_NEFA <- sqrt(0.0102); label("Additive residual SD on NEFA (mmol/L)") # Table 4 separate SD: sigma2 = 0.0102 (RSE 2.45%), read as a variance
  })

  model({
    vmax_hi <- exp(lvmax_hi + etalvmax_hi)
    km_hi <- exp(lkm_hi)
    vmax_lo <- exp(lvmax_lo + etalvmax_lo)
    km_lo <- exp(lkm_lo)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    synt <- exp(lsynt + etalsynt)
    ka <- exp(lka + etalka)
    tmax_abs <- exp(ltmax_abs + etaltmax_abs)
    kt_abs <- exp(lkt_abs)

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
    Ct <- peripheral1 / vp

    # Eq. 1: gut amount with parallel linear (A1) and saturable (A2) loss.
    absLin <- ka * depot
    absSat <- tmax_abs * depot / (kt_abs + depot)
    d/dt(depot) <- -absLin - absSat
    # Eqs. 2-3 written as amounts per kg (Vc * dCp/dt); IV infusions dose central.
    d/dt(central) <- absLin + absSat + synt - vmax_hi * Cp / (km_hi + Cp) - vmax_lo * Cp / (km_lo + Cp) - q * Cp + q * Ct
    d/dt(peripheral1) <- q * Cp - q * Ct

    # Endogenous NiAc steady state of Eqs. 2-3 with no dose: the single
    # positive root of a*Cp^2 + b*Cp + c = 0 (a < 0 < c).
    aq <- synt - vmax_hi - vmax_lo
    bq <- synt * (km_hi + km_lo) - vmax_hi * km_lo - vmax_lo * km_hi
    cq <- synt * km_hi * km_lo
    Cp0 <- (-bq - sqrt(bq * bq - 4 * aq * cq)) / (2 * aq)
    central(0) <- Cp0 * vc
    peripheral1(0) <- Cp0 * vp

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
