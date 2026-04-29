PerezRuixo_2025_posdinemab <- function() {
  description <- "Mechanism-based population PK-PD model with full TMDD for the anti-tau monoclonal antibody posdinemab in serum, CSF, and ISF (Perez-Ruixo 2025): two-compartment serum disposition with linear elimination, distribution into a CSF compartment and a downstream ISF compartment, explicit second-order binding of free posdinemab to free p217+tau in CSF and to tau seeds in ISF, internalization of free target and drug-target complex, and Alzheimer's-disease-vs-healthy effect on baseline p217+tau."
  reference <- "Perez-Ruixo C, Liu L, Galpern WR, Perez-Ruixo JJ. Mechanistic Population Pharmacokinetic-Pharmacodynamic Model of the Tau-Targeted Antibody Posdinemab in Healthy Participants and Participants with Alzheimer's Disease. Clin Pharmacol Ther. 2026;119(4):979-990. doi:10.1002/cpt.70173"
  vignette <- "PerezRuixo_2025_posdinemab"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "pmol/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Allometric scaling exponents fixed at 0.75 for clearances (CL, Q, QCSF, QISF) and 1.00 for volumes (V1, V2, VCSF, VISF), normalised to 70 kg.",
      source_name        = "WT"
    ),
    DIS_AD = list(
      description        = "Alzheimer's disease patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant)",
      notes              = "Multiplicative log-shift on baseline free p217+tau in CSF (R0): exp(e_ad_r0) = 5.995 / 0.793 = 7.56-fold higher in AD vs healthy. No PK-parameter effects.",
      source_name        = "STATUS"
    )
  )

  population <- list(
    n_subjects     = 69L,
    n_studies      = 1L,
    study          = "Phase 1 first-in-human dose-escalation, NCT03375697 (Galpern 2024)",
    age_range      = "55-78 years (median 67)",
    height_range   = "150-192 cm (median 169)",
    weight_range   = "51-106 kg (median 75)",
    sex_female_pct = 46.4,
    race_ethnicity = "98.6% White (98.5% Caucasian, 1.5% Hispanic of White), 1.4% Asian (Table 1)",
    disease_state  = "Healthy adults (n = 56, 81.2%) and adults with Alzheimer's disease (n = 13, 18.8%)",
    dose_range     = "Single ascending dose 1, 3, 10, 30, 60 mg/kg IV (healthy); multiple ascending dose 5, 15, 30, or 50 mg/kg IV every 28 days for 13 weeks (5 and 50 mg/kg in healthy; 15 and 30 mg/kg in AD); placebo n = 16",
    regions        = "Multinational Phase 1 (Galpern 2024)",
    notes          = "53 active-treatment participants contributed 907 serum and 169 CSF posdinemab concentration observations; all 69 contributed 294 CSF p217+tau observations (free + total combined). Analysis used FOCE in NONMEM 7.3.0."
  )

  ini({
    # ----- Posdinemab serum PK (Table 2; allometrically scaled to 70 kg) -----
    lcl    <- log(9.21e-3);  label("Posdinemab clearance from central (CL, L/h, 70 kg ref)")            # Table 2: 9.21 x 10^-3 L/h
    lv1    <- log(3.14);     label("Posdinemab central volume of distribution (V1, L, 70 kg ref)")      # Table 2: 3.14 L
    lq     <- log(24.9e-3);  label("Posdinemab intercompartmental clearance central<->peripheral1 (Q, L/h, 70 kg ref)") # Table 2: 24.9 x 10^-3 L/h
    lv2    <- log(2.87);     label("Posdinemab peripheral1 volume of distribution (V2, L, 70 kg ref)")   # Table 2: 2.87 L

    # ----- Posdinemab CSF / ISF distribution (Table 2; allometric) -----------
    lqcsf  <- log(4.02e-6);  label("Posdinemab intercompartmental clearance central<->CSF (QCSF, L/h, 70 kg ref)") # Table 2: 4.02 x 10^-6 L/h
    lvcsf  <- log(229e-3);   label("Posdinemab CSF volume of distribution (VCSF, L, 70 kg ref)")        # Table 2: 229 x 10^-3 L = 0.229 L
    lqisf  <- log(1.83e-3);  label("Posdinemab intercompartmental clearance CSF<->ISF (QISF, L/h, 70 kg ref)") # Table 2: 1.83 x 10^-3 L/h
    lvisf  <- log(43.4e-3);  label("Posdinemab ISF volume of distribution (VISF, L, 70 kg ref)")        # Table 2: 43.4 x 10^-3 L = 0.0434 L

    # ----- Allometric exponents (Results: fixed; Germovsek 2021 priors) ------
    allo_cl <- 0.75;  label("Allometric exponent on clearances (fixed; Results: estimating exponents did not improve MOFV)")
    allo_v  <- 1.00;  label("Allometric exponent on volumes (fixed; Results)")

    # ----- Mechanistic p217+tau / tau-seed parameters (Table 2) --------------
    lr0    <- log(0.793);                  label("Baseline free p217+tau in CSF, healthy (R0_HV, pmol/L)") # Table 2: R0 healthy 0.793 pmol/L
    e_ad_r0 <- log(5.995/0.793);           label("Log-shift on R0 for AD (additive on log-scale; exp = 7.56-fold higher)") # Table 2: R0 AD 5.995 pmol/L => log(5.995/0.793) = 2.023
    lkc    <- log(0.040);                  label("First-order elimination of free p217+tau / tau seeds (kc, 1/h)") # Table 2: 0.040 1/h
    lkint  <- log(0.299);                  label("First-order elimination of free posdinemab in CSF and posdinemab-target complex (kint, 1/h)") # Table 2: 0.299 1/h

    # kon reported in Table 2 as 264 nmol/mL^-1 / hour. Converting to (pmol/L)^-1 / h:
    # 1 nmol/mL = 10^6 pmol/L => kon[(pmol/L)^-1/h] = 264 * 10^-6 = 2.64e-4
    # Sanity check: kd = koff/kon = 0.224 / 2.64e-4 = 848.5 pmol/L (matches Discussion: 848.5 pM).
    lkon   <- log(2.64e-4);                label("Posdinemab-p217+tau second-order association rate in CSF (kon, (pmol/L)^-1 h^-1; converted from 264 (nmol/mL)^-1 h^-1)") # Table 2: 264 (nmol/mL)^-1 h^-1
    lkoff  <- log(0.224);                  label("Posdinemab-p217+tau dissociation rate in CSF (koff, 1/h)") # Table 2: 0.224 1/h

    # Fixed mechanistic assumptions (Methods, "Mechanism-based popPK-PD model"):
    aff_isf_ratio <- 20; label("Affinity ratio kd(CSF)/kd(ISF), fixed at 20 (Methods: ISF binds tau seeds with 20-fold higher affinity)")
    seed_ratio    <- 10; label("Baseline ISF tau seed / CSF p217+tau ratio, fixed at 10 (Methods: ISF tau seed levels 10-fold higher than CSF)")

    # MW for mg-to-pmol dose conversion. Posdinemab is a humanized IgG1 ~148 kDa
    # (Discussion: "13,805 pmol of posdinemab (148 kDa)").
    mw_da <- 148000; label("Posdinemab molecular weight (Da); used for dose mg-to-pmol conversion only")

    # ----- IIV (Table 2; CV% -> omega^2 = log(1 + CV^2)) ---------------------
    etalcl   ~ 0.05181  # CL CV 23.4%
    etalv1   ~ 0.03012  # V1 CV 17.5%
    etalq    ~ 0.18525  # Q  CV 44.9%
    etalv2   ~ 0.04764  # V2 CV 22.1%
    etalqcsf ~ 0.08300  # QCSF CV 29.4%
    etalvcsf ~ 0.06244  # VCSF CV 25.5%
    etalvisf ~ 0.60466  # VISF CV 91.0%
    etalr0   ~ 0.37113  # R0   CV 67.7%
    etalkc   ~ 0.26156  # kc   CV 54.7%
    etalkint ~ 0.13510  # kint CV 38.2%

    # ----- Residual unexplained variability (Table 2) ------------------------
    # Paper Methods: "an additive error model after natural logarithmic
    # transformation was used", which is proportional in linear space. Table 2
    # sigma values match the IIV CV% column convention (bare percentages).
    CcpropSd       <- 0.0873; label("Proportional residual error on serum posdinemab (sigma_1 = 8.73)")        # Table 2: sigma_1 8.73
    CcsfpropSd     <- 0.164;  label("Proportional residual error on CSF posdinemab (sigma_2 = 16.4)")          # Table 2: sigma_2 16.4
    TotalTaupropSd <- 0.112;  label("Proportional residual error on CSF total p217+tau (sigma_3 = 11.2)")      # Table 2: sigma_3 11.2
    FreeTaupropSd  <- 0.133;  label("Proportional residual error on CSF free p217+tau (sigma_4 = 13.3)")       # Table 2: sigma_4 13.3
  })

  model({
    # Allometric scaling to 70 kg (Results: fixed exponents).
    wt_cl <- (WT / 70)^allo_cl
    wt_v  <- (WT / 70)^allo_v

    # Individual structural parameters.
    cl    <- exp(lcl   + etalcl)   * wt_cl
    v1    <- exp(lv1   + etalv1)   * wt_v
    q     <- exp(lq    + etalq)    * wt_cl
    v2    <- exp(lv2   + etalv2)   * wt_v
    qcsf  <- exp(lqcsf + etalqcsf) * wt_cl
    vcsf  <- exp(lvcsf + etalvcsf) * wt_v
    qisf  <- exp(lqisf)            * wt_cl
    visf  <- exp(lvisf + etalvisf) * wt_v

    # Mechanistic parameters.
    kc    <- exp(lkc   + etalkc)
    kint  <- exp(lkint + etalkint)
    kon   <- exp(lkon)
    koff  <- exp(lkoff)
    R0    <- exp(lr0 + e_ad_r0 * DIS_AD + etalr0)

    # ISF binding rate constants implementing the 20-fold higher ISF affinity
    # (kd_ISF = kd_CSF / aff_isf_ratio). Decomposition: hold kon constant
    # (mAb-target on-rates are typically near the diffusion limit) and reduce
    # koff by aff_isf_ratio. See vignette "Assumptions and deviations".
    kon_isf  <- kon
    koff_isf <- koff / aff_isf_ratio

    # Endogenous p217+tau / tau-seed production rates (pmol/L/h). The paper's
    # printed equation uses one symbol "ksyn" for both compartments but the
    # initial condition RISF(0) = 10 * R(0) requires that ksyn for ISF be 10x
    # larger than for CSF so that the system is in baseline steady state.
    # See vignette "Errata".
    ksyn_csf <- kc * R0
    R0_isf   <- seed_ratio * R0
    ksyn_isf <- kc * R0_isf

    # Concentrations (pmol/L). Drug compartments are tracked in pmol so that
    # central / v1 yields pmol/L directly.
    Cc   <- central    / v1
    Cp   <- peripheral1 / v2
    Ccsf <- ACSF       / vcsf
    Cisf <- AISF       / visf

    # ODE system. Drug compartments (central, peripheral1, ACSF, AISF) are in
    # pmol; target compartments (R, RC, RISF, RCISF) are concentrations in
    # pmol/L. Equations 4 and 7 of the supplement are implemented in their
    # mass-balanced form (see vignette "Errata"); equations 1-3, 5, 6, 8 are
    # printed correctly.
    d/dt(central)    <- -cl * Cc - q * Cc + q * Cp - qcsf * Cc + qcsf * Ccsf
    d/dt(peripheral1) <-  q * Cc - q * Cp
    d/dt(ACSF)       <-  qcsf * Cc - qcsf * Ccsf - qisf * Ccsf + qisf * Cisf -
                          kon * R * ACSF + koff * RC * vcsf
    d/dt(AISF)       <-  qisf * Ccsf - qisf * Cisf -
                          kon_isf * RISF * AISF + koff_isf * RCISF * visf

    # Free p217+tau in CSF (pmol/L).
    d/dt(R)          <-  ksyn_csf - kc * R - kon * R * Ccsf + koff * RC

    # Posdinemab-p217+tau complex in CSF (pmol/L).
    d/dt(RC)         <-  kon * R * Ccsf - (koff + kint) * RC

    # Free tau seeds in ISF (pmol/L).
    d/dt(RISF)       <-  ksyn_isf - kc * RISF - kon_isf * RISF * Cisf + koff_isf * RCISF

    # Posdinemab-tau seed complex in ISF (pmol/L).
    d/dt(RCISF)      <-  kon_isf * RISF * Cisf - (koff_isf + kint) * RCISF

    # Initial conditions for the target species. Drug compartments default to 0;
    # IV bolus places the dose into central via f(central) below.
    R(0)    <- R0
    RISF(0) <- R0_isf

    # Bioavailability used as a unit-conversion factor: input dose in mg,
    # central state in pmol. 1 mg of posdinemab = 1e9 / mw_da pmol.
    f(central) <- 1e9 / mw_da

    # Outputs (pmol/L).
    TotalTau <- R + RC
    FreeTau  <- R

    Cc       ~ prop(CcpropSd)
    Ccsf     ~ prop(CcsfpropSd)
    TotalTau ~ prop(TotalTaupropSd)
    FreeTau  ~ prop(FreeTaupropSd)
  })
}
