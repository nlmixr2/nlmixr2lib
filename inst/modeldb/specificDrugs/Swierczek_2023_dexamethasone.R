Swierczek_2023_dexamethasone <- function() {
  description <- paste(
    "Translational population PK/PD model of the anti-inflammatory effects",
    "of oral dexamethasone in hospitalized adults with COVID-19. PK: two-",
    "compartment model with first-order absorption after a lag time and",
    "linear elimination (apparent clearances and volumes, oral F unknown),",
    "fitted to 167 blood concentrations from 27 patients after a single 6 mg",
    "oral dose. PD: three indirect-response turnover pools -- TNF-alpha and",
    "IL-6 with zero-order production inhibited by dexamethasone (Imax/IC50",
    "on the central blood concentration) and CRP whose production is",
    "linearly stimulated by IL-6. The PD parameters were not fitted to",
    "COVID-19 data: Imax/IC50 were translated from the authors' rat LPS",
    "model, kout values from human LPS endotoxemia, baselines from COVID-19",
    "non-survivors, and the PD IIV was set by assumption or trial-and-error",
    "simulation. No covariates."
  )
  reference <- paste(
    "Swierczek A, Jusko WJ. Anti-inflammatory effects of dexamethasone in",
    "COVID-19 patients: Translational population PK/PD modeling and",
    "simulation. Clin Transl Sci. 2023;16(9):1667-1679.",
    "doi:10.1111/cts.13577.",
    sep = " "
  )
  vignette <- "Swierczek_2023_dexamethasone"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (dexamethasone); pg/mL (tnf, il6); mg/L (crp)"
  )

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "dexamethasone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dexamethasone", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "dexamethasone", units = "mg", specimen = "tissue", verified = TRUE),
    tnf = list(analyte = "TNF-alpha", units = "pg/mL", specimen = "plasma", verified = TRUE),
    il6 = list(analyte = "IL-6", units = "pg/mL", specimen = "plasma", verified = TRUE),
    crp = list(analyte = "C-reactive protein", units = "mg/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27,
    n_studies = 1,
    age_mean = "63.5 years (30 enrolled patients)",
    weight_mean = "81 kg (30 enrolled patients)",
    disease_state = "hospitalized adults with COVID-19",
    dose_range = "PK: single 6 mg oral dexamethasone tablet, fasting. Simulations: 1.5, 3, 6 or 12 mg orally once daily for 10 days",
    regions = "Switzerland (Geneva University Hospitals)",
    notes = paste(
      "PK data: the prospective observational study of Abouir et al.",
      "(NCT04996784; 30 hospitalized patients with COVID-19, male and",
      "female) with capillary blood samples at 0, 0.5, 1, 2, 4, 6 and 8 h",
      "after a single 6 mg dose. Three profiles with a pre-dose",
      "concentration > 50 ng/mL were excluded, leaving 167 concentrations",
      "from 27 patients. No individual covariates were available. The PD",
      "layer has no fitted population: Imax and IC50 come from the authors'",
      "rat LPS model, kout values and the IL-6 -> CRP stimulation from",
      "human LPS endotoxemia (Liu et al.), and the baseline TNF-alpha,",
      "IL-6 and CRP from the median of COVID-19 non-survivors at Tongji",
      "Hospital, Wuhan."
    )
  )

  ini({
    # PK -- Table 1 (Monolix SAEM). The Supplement 1 Simulx file lists the
    # same values in its 'mlx_Pop' block.
    ltlag <- log(0.29)
    label("Absorption lag time tlag (h)") # Table 1, tlag = 0.29 h (RSE 20.5%)
    lka <- log(2.24)
    label("First-order absorption rate constant ka (1/h)") # Table 1, ka = 2.24 1/h (RSE 42.2%)
    lcl <- log(6.27)
    label("Apparent clearance CL/F (L/h)") # Table 1, CL/F = 6.27 L/h (RSE 11.7%)
    lvc <- log(20.83)
    label("Apparent central volume Vc/F (L)") # Table 1, Vc/F = 20.83 L (RSE 25.9%)
    lq <- log(19.72)
    label("Apparent distributional clearance CLd/F (L/h)") # Table 1, CLd/F = 19.72 L/h (RSE 48.6%)
    lvp <- log(13.02)
    label("Apparent peripheral volume Vp/F (L)") # Table 1, Vp/F = 13.02 L (RSE 26.7%)

    # PD -- Table 2. kin values are secondary parameters (kin = R0 * kout;
    # Table 2 footnote a) and carry no IIV.
    lkin_tnf <- fixed(log(3.890))
    label("Zero-order TNF-alpha production rate kinTNF (pg/mL/h)") # Table 2, kinTNF = 3.890 (secondary; = R0TNF 10.9 x koutTNF 0.357)
    lkout_tnf <- fixed(log(0.357))
    label("First-order TNF-alpha degradation rate constant koutTNF (1/h)") # Table 2, koutTNF = 0.357 (from Liu et al., ref 30)
    lic50_tnf <- fixed(log(2.25))
    label("Dexamethasone blood concentration for 50% of Imax on TNF-alpha production IC50DEX(TNF) (ng/mL)") # Table 2, IC50DEX(TNF) = 2.25 (from rat LPS model, ref 29)
    logitimax_tnf <- fixed(qlogis(0.805))
    label("Logit of maximum inhibition of TNF-alpha production ImaxDEX(TNF) (unitless)") # Table 2, ImaxDEX(TNF) = 0.805 (ref 29)
    lkin_il6 <- fixed(log(18.21))
    label("Zero-order IL-6 production rate kinIL6 (pg/mL/h)") # Table 2, kinIL6 = 18.21 (secondary; = R0IL6 59.7 x koutIL6 0.305)
    lkout_il6 <- fixed(log(0.305))
    label("First-order IL-6 degradation rate constant koutIL6 (1/h)") # Table 2, koutIL6 = 0.305 (ref 30)
    lic50_il6 <- fixed(log(3.35))
    label("Dexamethasone blood concentration for 50% of Imax on IL-6 production IC50DEX(IL6) (ng/mL)") # Table 2, IC50DEX(IL6) = 3.35 (ref 29)
    logitimax_il6 <- fixed(qlogis(0.83))
    label("Logit of maximum inhibition of IL-6 production ImaxDEX(IL6) (unitless)") # Table 2, ImaxDEX(IL6) = 0.83 (ref 29)
    lkin_crp <- fixed(log(0.0016))
    label("Zero-order CRP production rate kinCRP (mg/L/h)") # Table 2, kinCRP = 0.0016 (secondary)
    lkout_crp <- fixed(log(0.0183))
    label("First-order CRP degradation rate constant koutCRP (1/h)") # Table 2, koutCRP = 0.0183 (ref 30)
    lsstim_il6_crp <- fixed(log(19.2))
    label("Linear stimulation of CRP production by IL-6 SIL6onCRP (mL/pg)") # Table 2, SIL6onCRP = 19.2 (ref 30)

    # IIV -- Monolix 'SD of random effects' (Tables 1-2), squared to
    # variances. PK IIV was estimated; PD IIV was assumed (IC50/Imax) or set
    # by trial-and-error simulation (kout) and is therefore fixed.
    etaltlag ~ 0.0256 # Table 1, SD 0.16 (normal IIV: added to tlag in h, not to ltlag)
    etalka ~ 2.8561 # Table 1, SD 1.69
    etalcl ~ 0.3025 # Table 1, SD 0.55
    etalvc ~ 0.3844 # Table 1, SD 0.62
    etalq ~ 0.4489 # Table 1, SD 0.67
    etalvp ~ 0.16 # Table 1, SD 0.4
    etalkout_tnf ~ fixed(0.01147041) # Table 2, SD 0.1071 (trial-and-error)
    etalic50_tnf ~ fixed(0.2025) # Table 2, SD 0.45 (assumed)
    etalogitimax_tnf ~ fixed(0.027556) # Table 2, SD 0.166 on the logit scale (assumed)
    etalkout_il6 ~ fixed(0.033489) # Table 2, SD 0.183 (trial-and-error)
    etalic50_il6 ~ fixed(0.49) # Table 2, SD 0.7 (assumed)
    etalogitimax_il6 ~ fixed(0.025921) # Table 2, SD 0.161 on the logit scale (assumed)
    etalkout_crp ~ fixed(0.00001369) # Table 2, SD 0.0037 (trial-and-error)

    addSd <- 14.76
    label("Additive residual error on blood dexamethasone (ng/mL)") # Table 1, a = 14.76 ng/mL (RSE 8.42%)
  })
  model({
    # Supplement 1 Monolix code: Tlag is normally distributed
    # (Tlag = Tlag_pop + eta), all other PK and kout / IC50 parameters are
    # log-normal, and Imax is logit-normal. etaltlag keeps the canonical
    # eta name but is added to tlag in hours, not to ltlag.
    tlag <- exp(ltlag) + etaltlag
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kin_tnf <- exp(lkin_tnf)
    kout_tnf <- exp(lkout_tnf + etalkout_tnf)
    ic50_tnf <- exp(lic50_tnf + etalic50_tnf)
    imax_tnf <- expit(logitimax_tnf + etalogitimax_tnf)
    kin_il6 <- exp(lkin_il6)
    kout_il6 <- exp(lkout_il6 + etalkout_il6)
    ic50_il6 <- exp(lic50_il6 + etalic50_il6)
    imax_il6 <- expit(logitimax_il6 + etalogitimax_il6)
    kin_crp <- exp(lkin_crp)
    kout_crp <- exp(lkout_crp + etalkout_crp)
    sstim_il6_crp <- exp(lsstim_il6_crp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg, volumes in L: mg/L x 1000 = ng/mL (blood; the authors
    # assume a blood-to-plasma ratio of 1).
    Cc <- 1000 * central / vc

    # Eqs 4-6. The drug effect is driven by the central blood concentration
    # Cc, as in the Supplement 1 Monolix code and the Table 2 definition of
    # IC50 as a 'DEX plasma concentration'. Each pool starts at its own
    # pre-dose steady state, kin / kout, so the IIV on kout also spreads the
    # baselines (Methods, 'Verification of the model').
    d/dt(tnf) <- kin_tnf * (1 - imax_tnf * Cc / (ic50_tnf + Cc)) - kout_tnf * tnf
    d/dt(il6) <- kin_il6 * (1 - imax_il6 * Cc / (ic50_il6 + Cc)) - kout_il6 * il6
    d/dt(crp) <- kin_crp * (1 + sstim_il6_crp * il6) - kout_crp * crp
    tnf(0) <- kin_tnf / kout_tnf
    il6(0) <- kin_il6 / kout_il6
    crp(0) <- kin_crp * (1 + sstim_il6_crp * kin_il6 / kout_il6) / kout_crp

    Cc ~ add(addSd)
  })
}
