Ngo_2020_HL2351 <- function() {
  description <- "Population PK model for HL2351 (hIL-1Ra-hyFc, ~97 kDa) in healthy adult Korean men: a quasi-steady-state target-mediated drug disposition (QSS-TMDD) model coupled with FcRn-mediated recycling. The injection-site depot feeds a separate distribution space where free drug equilibrates with FcRn (QSS dissociation constant AKSS1, total FcRn AFcRn_t); free drug moves to the central compartment either directly (Ka2) or by FcRn-mediated recycling of the FcRn-drug complex (Krec). In the central compartment free drug equilibrates with IL1R (QSS dissociation constant KSS2, total IL1R CIL1R_t), is taken up back to the distribution space (Kup), exchanged with one peripheral compartment (Q/F), and eliminated linearly (CL/F). The IL1R-drug complex degrades at Kdeg2. All drug amounts and concentrations are in nmol / nmol/L; convert mg dosing using molecular weight 97 kDa (1 mg HL2351 = approximately 10306 nmol)."
  reference <- "Ngo L, Lee J, Lim L, Lim H, Bae KS, Hong T, Bae S, Hong Y. Development of a Pharmacokinetic Model Describing Neonatal Fc Receptor-Mediated Recycling of HL2351, a Novel Hybrid Fc-Fused Interleukin-1 Receptor Antagonist, to Optimize Dosage Regimen. CPT Pharmacometrics Syst Pharmacol. 2020 Oct;9(10):584-595. doi:10.1002/psp4.12552. PMID: 32945613."
  vignette <- "Ngo_2020_HL2351"
  units <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  population <- list(
    n_subjects     = 40,
    n_studies      = 1,
    study_id       = "NCT02175056 (Phase I single-ascending-dose, Seoul National University Hospital)",
    age_range      = "20-45 years",
    weight_range   = "not reported",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult Korean men",
    dose_range     = "Single SC dose: HL2351 1, 2, 4, 8, or 12 mg/kg (n = 8 per dose group)",
    regions        = "Republic of Korea",
    notes          = "Phase I single-ascending-dose first-in-human study NCT02175056 (Ngo 2020 Methods, Materials and PK data collection). The full study had six cohorts (HL2351 1, 2, 4, 8, 12 mg/kg and anakinra 100 mg, n = 8 each); the HL2351 model in this file is fit to the five HL2351 cohorts only (n = 40). Anakinra parameters (Ngo 2020 Table S2 in the supplement) are not implemented here. PK observations were collected to 672 hours postdose for HL2351 (Ngo 2020 Methods)."
  )

  ini({
    # Structural parameters - typical adult, healthy Korean men. All values
    # are the observed-data medians from Ngo 2020 Table 1; bootstrap medians
    # are reported for context but are not implemented as a separate model.
    lka1      <- log(1.21);     label("Absorption rate from injection-site depot to distribution space (Ka1, 1/h)")  # Ngo 2020 Table 1: Ka1 = 1.21 1/h
    lka2      <- log(0.0171);   label("Direct rate of free drug from distribution space to central (Ka2, 1/h)")     # Ngo 2020 Table 1: Ka2 = 0.0171 1/h
    lkrec     <- log(0.0338);   label("FcRn-mediated recycling rate of FcRn-drug complex into central (Krec, 1/h)") # Ngo 2020 Table 1: Krec = 0.0338 1/h
    lkdeg1    <- log(0.0264);   label("Degradation rate of free drug at the distribution space (Kdeg1, 1/h)")       # Ngo 2020 Table 1: Kdeg1 = 0.0264 1/h
    lkdeg2    <- fixed(log(0.206));  label("Degradation rate of IL1R-drug complex at the central compartment (Kdeg2, 1/h; FIXED)") # Ngo 2020 Table 1: Kdeg2 = 0.206 1/h (FIX)
    lkup      <- fixed(log(0.00952)); label("Uptake rate of free drug from central back to distribution space (Kup, 1/h; FIXED)") # Ngo 2020 Table 1: Kup = 0.00952 1/h (FIX)
    lcl       <- log(0.208);    label("Apparent clearance of free drug from central (CL/F, L/h)")                   # Ngo 2020 Table 1: CL/F = 0.208 L/h
    lvc       <- log(11.3);     label("Apparent central volume of distribution (Vc/F, L)")                          # Ngo 2020 Table 1: Vc/F = 11.3 L
    lq        <- log(0.0288);   label("Apparent inter-compartmental clearance (Q/F, L/h)")                          # Ngo 2020 Table 1: Q/F = 0.0288 L/h
    lvp       <- fixed(log(5.06));  label("Apparent peripheral volume of distribution (Vp/F, L; FIXED)")            # Ngo 2020 Table 1: Vp/F = 5.06 L (FIX)
    la_kss1   <- log(237);      label("Drug-FcRn QSS dissociation constant in distribution space (AKSS1, nmol)")    # Ngo 2020 Table 1: AKSS1 = 237 nmol
    la_fcrn_t <- log(749);      label("Total active amount of FcRn in distribution space (AFcRn_t, nmol)")          # Ngo 2020 Table 1: AFcRn_t = 749 nmol
    lkss2     <- log(14.5);     label("Drug-IL1R QSS dissociation constant in central (KSS2, nmol/L)")              # Ngo 2020 Table 1: KSS2 = 14.5 nmol/L
    lc_il1r_t <- log(2.23);     label("Total active concentration of IL1R in central (CIL1R_t, nmol/L)")            # Ngo 2020 Table 1: CIL1R_t = 2.23 nmol/L
    lalag     <- log(0.312);    label("Subcutaneous absorption lag time (Alag, h)")                                 # Ngo 2020 Table 1: Alag = 0.312 h

    # IIV - Ngo 2020 Table 1 IIV section reports CV%; omega^2 = log(CV^2 + 1)
    etalka1   ~ 0.6672   # 97.4% CV  -> log(1 + 0.974^2)  = 0.6672
    etalka2   ~ 0.4215   # 72.4% CV  -> log(1 + 0.724^2)  = 0.4215
    etalkrec  ~ 0.1046   # 33.2% CV  -> log(1 + 0.332^2)  = 0.1046
    etalkdeg1 ~ 0.0729   # 27.5% CV  -> log(1 + 0.275^2)  = 0.0729
    etalkup   ~ 0.5973   # 90.4% CV  -> log(1 + 0.904^2)  = 0.5973 (typical Kup is FIXED but IIV is estimated; Ngo 2020 Table 1)
    etalcl    ~ 0.0511   # 22.9% CV  -> log(1 + 0.229^2)  = 0.0511
    etalalag  ~ 0.2685   # 55.5% CV  -> log(1 + 0.555^2)  = 0.2685

    # Residual error - Ngo 2020 Table 1 residual variability section. Observed
    # serum concentrations are in nmol/L (additive units must match concentration units).
    addSd  <- 0.184; label("Additive residual error (nmol/L)")           # Ngo 2020 Table 1: Additive = 0.184 nmol/L
    propSd <- 0.115; label("Proportional residual error (fraction)")     # Ngo 2020 Table 1: Proportional = 11.5%
  })

  model({
    # Individual parameters
    ka1      <- exp(lka1   + etalka1)
    ka2      <- exp(lka2   + etalka2)
    krec     <- exp(lkrec  + etalkrec)
    kdeg1    <- exp(lkdeg1 + etalkdeg1)
    kdeg2    <- exp(lkdeg2)
    kup      <- exp(lkup   + etalkup)
    cl       <- exp(lcl    + etalcl)
    vc       <- exp(lvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    a_kss1   <- exp(la_kss1)
    a_fcrn_t <- exp(la_fcrn_t)
    kss2     <- exp(lkss2)
    c_il1r_t <- exp(lc_il1r_t)
    alag_t   <- exp(lalag  + etalalag)

    # Linear-disposition micro-rate constants from CL/F, Vc/F, Q/F, Vp/F
    kel <- cl / vc
    kcp <- q  / vc
    kpc <- q  / vp

    # QSS approximation - Eqs. (4) and (5) of Ngo 2020.
    # The state variables are TOTAL drug amounts: dist holds free drug + FcRn-drug
    # complex (Atot1, nmol); central holds free drug + IL1R-drug complex
    # (total amount = Ctot2 * Vc, nmol). Free amounts/concentrations are derived
    # from the totals via the QSS quadratic.
    diff1     <- dist - a_fcrn_t - a_kss1
    a_dfree1  <- 0.5 * (diff1 + sqrt(diff1 * diff1 + 4 * a_kss1 * dist))
    a_fcrn_d  <- dist - a_dfree1                                     # FcRn-drug complex amount (nmol)

    c_tot2    <- central / vc
    diff2     <- c_tot2 - c_il1r_t - kss2
    c_dfree2  <- 0.5 * (diff2 + sqrt(diff2 * diff2 + 4 * kss2 * c_tot2))
    a_dfree2  <- c_dfree2 * vc                                       # free drug amount in central (nmol)
    a_il1r_d  <- (c_tot2 - c_dfree2) * vc                            # IL1R-drug complex amount (nmol)

    # ODEs - mass balance on TOTAL drug amounts in each space (Ngo 2020 Figure 1a)
    d/dt(depot)       <- -ka1 * depot
    d/dt(dist)        <-  ka1 * depot - krec * a_fcrn_d - ka2 * a_dfree1 - kdeg1 * a_dfree1 + kup * a_dfree2
    d/dt(central)     <-  krec * a_fcrn_d + ka2 * a_dfree1 - kup * a_dfree2 - kel * a_dfree2 - kdeg2 * a_il1r_d - kcp * a_dfree2 + kpc * peripheral1
    d/dt(peripheral1) <-  kcp * a_dfree2 - kpc * peripheral1

    # Subcutaneous absorption lag time
    alag(depot) <- alag_t

    # Observation: total drug concentration in the central compartment (free + IL1R-bound)
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
