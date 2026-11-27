Tapani_2014_nicotinicAcid_rat <- function() {
  description <- paste(
    "Preclinical (rat). Joint normal-obese PK/PD feedback model for",
    "nicotinic acid (NiAc) and plasma non-esterified fatty acids (NEFA) in",
    "male Sprague-Dawley and obese Zucker (fa/fa) rats after IV",
    "constant-rate, stepwise IV and oral NiAc. NiAc PK (fixed from earlier",
    "analyses) is a two-compartment model with two parallel Michaelis-Menten",
    "eliminations, endogenous synthesis and parallel linear + saturable gut",
    "absorption in normal rats, and a one-compartment model with one",
    "Michaelis-Menten elimination and endogenous synthesis in obese rats.",
    "NEFA turnover: NiAc inhibits NEFA formation (Hill-Imax, Imax = 1), with",
    "an eight-compartment moderator feedback chain (moderator1 inhibits",
    "formation, moderator8 stimulates loss) and a NiAc-independent",
    "formation term kcap. Baseline NEFA R0, fractional turnover kout and the",
    "sigmoidicity gamma differ between the strains (DIS_OBESE); the other",
    "NEFA parameters are shared."
  )
  reference <- paste(
    "Tapani S, Almquist J, Leander J, Ahlstrom C, Peletier LA, Jirstrand M,",
    "Gabrielsson J. Joint feedback analysis modeling of nonesterified fatty",
    "acids in obese Zucker rats and normal Sprague-Dawley rats after",
    "different routes of administration of nicotinic acid.",
    "J Pharm Sci. 2014;103(8):2571-2584. doi:10.1002/jps.24077"
  )
  vignette <- "Tapani_2014_nicotinicAcid_rat"
  dosing <- c("depot", "central")
  units <- list(
    time = "min",
    dosing = "umol/kg",
    concentration = "umol/L (NiAc); mmol/L (NEFA)"
  )

  covariateData <- list(
    DIS_OBESE = list(
      description = "Obese Zucker (fa/fa) rat indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal Sprague-Dawley rat)",
      notes = paste(
        "1 = obese Zucker (fa/fa) rat, 0 = normal Sprague-Dawley rat.",
        "The paper's indicator Z_i (Eq. 11) is 1 for normal and 0 for",
        "obese animals, so Z_i = 1 - DIS_OBESE. Group-specific values",
        "enter linearly as in Eq. 11, theta = (1 - DIS_OBESE) * theta_N +",
        "DIS_OBESE * theta_D, for R0, kout and gamma and for the separate",
        "normal and obese NiAc PK parameter sets. The fitted model only",
        "uses 0 and 1; the paper's 'degree of disease' simulations",
        "(Figure 7) set Z between 0 and 1 and weight the solutions of the",
        "two PK models, which this single-PK-structure implementation does",
        "not reproduce (see the vignette)."
      ),
      source_name = "Z (Eq. 11; Z = 1 - DIS_OBESE)"
    )
  )

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
    species = "rat (male Sprague-Dawley and male obese Zucker fa/fa)",
    n_subjects = 118L,
    n_studies = 3L,
    age_range = "16 weeks at experiment",
    weight_range = "220-367 g (Sprague-Dawley); 473-547 g (obese Zucker)",
    sex_female_pct = 0,
    disease_state = paste(
      "Normal Sprague-Dawley rats and obese, insulin-resistant Zucker",
      "(fa/fa) rats; fasted 14 h before and throughout the experiment."
    ),
    dose_range = paste(
      "Sprague-Dawley: NiAc 1, 5 or 20 umol/kg IV over 30 min; 5, 10 or",
      "51 umol/kg IV over 300 min; 24.4, 81.2 or 812 umol/kg oral gavage;",
      "5 umol/kg over 30 min followed by stepwise-decreasing infusion",
      "(with or without a second 5 umol/kg 30-min infusion at 210 min);",
      "vehicle controls. Obese Zucker: 20 umol/kg IV over 30 min or",
      "51 umol/kg IV over 300 min; vehicle controls."
    ),
    regions = "Sweden (AstraZeneca R&D Molndal)",
    notes = paste(
      "Pooled data from earlier NiAc/NEFA studies (Experimental Design).",
      "The group sizes printed for the Sprague-Dawley arms sum to about 99",
      "and for the obese Zucker arms to 19; the abstract says 'more than",
      "100' rats. 11-15 arterial samples per rat; 1440 NEFA observations",
      "(Results). NEFA by enzymatic colorimetric assay (Wako); NiAc by",
      "LC-MS/MS (LLOQ 1 nmol/L)."
    )
  )

  ini({
    # ---- NiAc PK, normal Sprague-Dawley rats (Table 1, fixed from Refs 4-5) ----
    lvmax_hi <- fixed(log(0.0871)); label("Vmax of high-affinity NiAc elimination, normal rats (umol/min/kg)") # Table 1 SD: Vmax1 = 0.0871 (RSE 22.8%)
    lkm_hi <- fixed(log(0.235)); label("Km of high-affinity NiAc elimination, normal rats (umol/L)") # Table 1 SD: Km1 = 0.235 (RSE 29.2%)
    lvmax_lo <- fixed(log(7.09)); label("Vmax of low-affinity NiAc elimination, normal rats (umol/min/kg)") # Table 1 SD: Vmax2 = 7.09 (RSE 39.6%)
    lkm_lo <- fixed(log(74.5)); label("Km of low-affinity NiAc elimination, normal rats (umol/L)") # Table 1 SD: Km2 = 74.5 (RSE 43.4%)
    lvc <- fixed(log(0.393)); label("Central volume Vc, normal rats (L/kg)") # Table 1 SD: Vc = 0.393 (RSE 5.29%)
    lvp <- fixed(log(0.172)); label("Peripheral volume Vt, normal rats (L/kg)") # Table 1 SD: Vt = 0.172 (RSE 35.2%)
    lq <- fixed(log(0.000852)); label("Intercompartmental clearance Cld, normal rats (L/min/kg)") # Table 1 SD: Cld = 0.000852 (RSE 27.8%)
    lsynt <- fixed(log(0.00355)); label("Endogenous NiAc synthesis rate Synt, normal rats (umol/min/kg)") # Table 1 SD: Synt = 0.00355 (RSE 23.3%)
    lka <- fixed(log(0.00477)); label("First-order gut absorption rate ka (1/min)") # Table 1 SD: ka = 0.00477 (RSE 33.5%)
    ltmax_abs <- fixed(log(2.96)); label("Maximum saturable gut absorption rate Vmax,g (umol/min/kg)") # Table 1 SD: Vmax,g = 2.96 (RSE 21.3%)
    lkt_abs <- fixed(log(20.5)); label("Gut amount at half-maximal saturable absorption Km,g (umol/kg)") # Table 1 SD: Km,g = 20.5 (RSE 29.6%)

    # ---- NiAc PK, obese Zucker rats (Table 1, fixed from Refs 4-5) ----
    lvmax_hi_obese <- fixed(log(1.59)); label("Vmax of NiAc elimination, obese rats (umol/min/kg)") # Table 1 obese: Vmax1 = 1.59 (RSE 13.9%)
    lkm_hi_obese <- fixed(log(18.9)); label("Km of NiAc elimination, obese rats (umol/L)") # Table 1 obese: Km1 = 18.9 (RSE 21.5%)
    lvc_obese <- fixed(log(0.323)); label("Central volume Vc, obese rats (L/kg)") # Table 1 obese: Vc = 0.323 (RSE 12.4%)
    lsynt_obese <- fixed(log(0.00280)); label("Endogenous NiAc synthesis rate Synt, obese rats (umol/min/kg)") # Table 1 obese: Synt = 0.00280 (RSE 10.1%)

    # ---- NEFA feedback model, joint analysis (Table 4 'Joint Analysis') ----
    lrbase <- log(0.741); label("Baseline NEFA R0, normal rats (mmol/L)") # Table 4 joint SD: R0 = 0.741 (RSE 4.42%)
    lrbase_obese <- log(1.07); label("Baseline NEFA R0, obese rats (mmol/L)") # Table 4 joint obese: R0 = 1.07 (RSE 6.34%)
    lkout <- log(0.290); label("Fractional NEFA turnover kout, normal rats (L/mmol/min)") # Table 4 joint SD: kout = 0.290 (RSE 10.1%)
    lkout_obese <- log(0.343); label("Fractional NEFA turnover kout, obese rats (L/mmol/min)") # Table 4 joint obese: kout = 0.343 (RSE 31.1%)
    lgam <- log(2.16); label("Sigmoidicity gamma of NiAc inhibition, normal rats (unitless)") # Table 4 joint SD: gamma = 2.16 (RSE 4.87%)
    lgam_obese <- log(0.298); label("Sigmoidicity gamma of NiAc inhibition, obese rats (unitless)") # Table 4 joint obese: gamma = 0.298 (RSE 9.63%)
    lktol <- log(0.0245); label("Moderator turnover rate ktol, both strains (1/min)") # Table 4 joint: ktol = 0.0245 (RSE 2.00%)
    lkcap <- log(0.0245); label("NiAc-independent NEFA formation kcap, both strains (mmol/L/min)") # Table 4 joint: kcap = 0.0245 (RSE 10.2%)
    lpmod <- log(1.20); label("Amplification factor p on moderator1, both strains (unitless)") # Table 4 joint: p = 1.20 (RSE 4.02%)
    lic50 <- log(0.0820); label("NiAc IC50 of NEFA-formation inhibition, both strains (umol/L)") # Table 4 joint: IC50 = 0.0820 (RSE 15.5%)
    imax <- fixed(1); label("Maximal NiAc inhibition of NEFA formation Imax (fraction)") # Table 4: Imax = 1 (fixed, Initial Parameter Estimates)

    # ---- IIV: the paper's IIV % is taken as 100 * sqrt(omega) (CV^2 below); see vignette ----
    etalvmax_hi ~ fixed(0.927^2) # Table 1 SD: IIV Vmax1 92.7%
    etalvmax_lo ~ fixed(0.291^2) # Table 1 SD: IIV Vmax2 29.1%
    etalsynt ~ fixed(1.09^2) # Table 1 SD: IIV Synt 109%
    etalka ~ fixed(0.108^2) # Table 1 SD: IIV ka 10.8%
    etaltmax_abs ~ fixed(0.108^2) # Table 1 SD: IIV Vmax,g 10.8%
    etalvmax_hi_obese ~ fixed(0.214^2) # Table 1 obese: IIV Vmax1 21.4%
    etalsynt_obese ~ fixed(0.953^2) # Table 1 obese: IIV Synt 95.3%
    etalrbase ~ 0.197^2 # Table 4 joint SD: IIV R0 19.7%
    etalrbase_obese ~ 0.240^2 # Table 4 joint obese: IIV R0 24.0%
    etalkout ~ 0.714^2 # Table 4 joint SD: IIV kout 71.4%
    etalkout_obese ~ 1.51^2 # Table 4 joint obese: IIV kout 151%
    etalic50 ~ 1.29^2 # Table 4 joint: IIV IC50 129%, shared by both strains

    # ---- Residual error ----
    propSd_normal <- fixed(0.428); label("Proportional residual SD on NiAc, normal rats (fraction)") # Table 1 SD: sigma1 = 42.8% (RSE 5.16%)
    propSd_obese <- fixed(0.400); label("Proportional residual SD on NiAc, obese rats (fraction)") # Table 1 obese: sigma1 = 40.0% (RSE 26.3%)
    addSd_NEFA <- sqrt(0.0110); label("Additive residual SD on NEFA (mmol/L)") # Table 4 joint: sigma2 = 0.0110 (RSE 2.84%), read as a variance
  })

  model({
    nrm <- 1 - DIS_OBESE

    # NiAc PK. Obese rats use the one-compartment, single-Michaelis-Menten
    # model of Eq. 4: the low-affinity pathway and the peripheral exchange
    # are switched off by DIS_OBESE.
    vmax_hi_n <- exp(lvmax_hi + etalvmax_hi)
    vmax_hi_d <- exp(lvmax_hi_obese + etalvmax_hi_obese)
    vmax_hi <- nrm * vmax_hi_n + DIS_OBESE * vmax_hi_d
    km_hi <- nrm * exp(lkm_hi) + DIS_OBESE * exp(lkm_hi_obese)
    vmax_lo <- nrm * exp(lvmax_lo + etalvmax_lo)
    km_lo <- exp(lkm_lo)
    vc <- nrm * exp(lvc) + DIS_OBESE * exp(lvc_obese)
    vp <- exp(lvp)
    q <- nrm * exp(lq)
    synt_n <- exp(lsynt + etalsynt)
    synt_d <- exp(lsynt_obese + etalsynt_obese)
    synt <- nrm * synt_n + DIS_OBESE * synt_d
    ka <- exp(lka + etalka)
    tmax_abs <- exp(ltmax_abs + etaltmax_abs)
    kt_abs <- exp(lkt_abs)

    # NEFA parameters; group-specific values combine linearly as in Eq. 11.
    rbase_n <- exp(lrbase + etalrbase)
    rbase_d <- exp(lrbase_obese + etalrbase_obese)
    rbase <- nrm * rbase_n + DIS_OBESE * rbase_d
    kout_n <- exp(lkout + etalkout)
    kout_d <- exp(lkout_obese + etalkout_obese)
    kout <- nrm * kout_n + DIS_OBESE * kout_d
    gam <- nrm * exp(lgam) + DIS_OBESE * exp(lgam_obese)
    ktol <- exp(lktol)
    kcap <- exp(lkcap)
    pmod <- exp(lpmod)
    ic50 <- exp(lic50 + etalic50)

    # Eq. 10: kin from the baseline steady state R = M1 = M8 = R0.
    kin <- (kout * rbase^2 - kcap) * rbase^pmod

    Cp <- central / vc
    Ct <- peripheral1 / vp

    # Eq. 1: gut amount with parallel linear (A1) and saturable (A2) loss.
    absLin <- ka * depot
    absSat <- tmax_abs * depot / (kt_abs + depot)
    d/dt(depot) <- -absLin - absSat
    # Eqs. 2-4 written as amounts per kg (Vc * dCp/dt); IV infusions dose central.
    d/dt(central) <- absLin + absSat + synt - vmax_hi * Cp / (km_hi + Cp) - vmax_lo * Cp / (km_lo + Cp) - q * Cp + q * Ct
    d/dt(peripheral1) <- q * Cp - q * Ct

    # Endogenous NiAc steady state of Eqs. 2-4 with no dose. Clearing the
    # denominators gives a*Cp^2 + b*Cp + c = 0 with a < 0 < c, whose single
    # positive root is taken; with vmax_lo = 0 (obese) it reduces to
    # synt * km_hi / (vmax_hi - synt).
    aq <- synt - vmax_hi - vmax_lo
    bq <- synt * (km_hi + km_lo) - vmax_hi * km_lo - vmax_lo * km_hi
    cq <- synt * km_hi * km_lo
    Cp0 <- (-bq - sqrt(bq * bq - 4 * aq * cq)) / (2 * aq)
    central(0) <- Cp0 * vc
    peripheral1(0) <- Cp0 * vp

    # Eq. 6: I(Cp) is the fraction of NEFA formation remaining. It is driven
    # by total (endogenous + dosed) NiAc.
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

    propSd <- nrm * propSd_normal + DIS_OBESE * propSd_obese

    Cc <- Cp
    NEFA <- nefa
    Cc ~ prop(propSd)
    NEFA ~ add(addSd_NEFA)
  })
}
