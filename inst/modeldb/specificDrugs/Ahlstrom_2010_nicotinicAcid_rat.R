Ahlstrom_2010_nicotinicAcid_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "PK/PD feedback model for nicotinic acid (NiAc) and non-esterified",
    "fatty acids (NEFA) in male Sprague-Dawley rats following IV",
    "infusions. NiAc disposition is a two-compartment model with two",
    "parallel capacity-limited (Michaelis-Menten) elimination processes",
    "(likely glycine conjugation and amidation) plus endogenous synthesis.",
    "NEFA turnover is described by an inhibitory drug-mechanism function",
    "(Hill-Imax, with Imax fixed at 1) acting on the formation of NEFA,",
    "coupled to a moderator feedback chain of 8 transit compartments",
    "(precursor1..precursor8): the first compartment inhibits NEFA",
    "formation amplified by exponent p and the last compartment",
    "stimulates NEFA loss. A NiAc-independent capillary release rate",
    "kcap sets the lower physiological limit of NEFA in plasma. All",
    "structural parameters are body-weight-normalized (per kg).",
    "Parameter values from Ahlstrom 2010 Tables 1 (NiAc PK) and 2",
    "(NEFA PD)."
  )
  reference <- paste(
    "Ahlstrom C, Peletier LA, Jansson-Lofmark R, Gabrielsson J. (2011).",
    "Feedback modeling of non-esterified fatty acids in rats after",
    "nicotinic acid infusions.",
    "J Pharmacokinet Pharmacodyn 38(1):1-24.",
    "doi:10.1007/s10928-010-9172-2.",
    sep = " "
  )
  vignette <- "Ahlstrom_2010_nicotinicAcid_rat"
  units <- list(
    time          = "min",
    dosing        = "umol/kg",
    concentration = "umol/L"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (male Sprague-Dawley)",
    n_subjects     = 63L,
    n_studies      = 1L,
    age_range      = "(not reported in the source publication)",
    weight_range   = "220-367 g",
    sex_female_pct = 0,
    disease_state  = "healthy; fasted 14 h prior to dosing and throughout the experiment",
    dose_range     = paste(
      "NiAc IV constant-rate infusion (jugular vein). 30 min infusion arm:",
      "vehicle (n=10), 1 umol/kg (n=4), 5 umol/kg (n=8), 20 umol/kg (n=9).",
      "300 min infusion arm: vehicle (n=8), 5 umol/kg (n=9), 10 umol/kg",
      "(n=8), 51 umol/kg (n=7). 8 groups total."
    ),
    regions        = "Sweden (AstraZeneca R&D Molndal)",
    notes          = paste(
      "63 male Sprague-Dawley rats from Harlan Nederlands B.V., housed and",
      "acclimatised for at least 1 week prior to surgery. Catheters were",
      "implanted in the left carotid artery for blood sampling and in the",
      "right external jugular vein for drug administration under",
      "isoflurane anaesthesia. Animals had 5-8 days post-surgery recovery,",
      "were fasted 14 h prior to dosing and throughout the experiment, and",
      "received NiAc dissolved in 0.9% NaCl as an IV constant-rate",
      "infusion. 11-13 arterial blood samples per rat were collected",
      "during and after infusion. NiAc was measured by LC-MS (LLOQ 1",
      "nmol/L, linear over 0.001-28 umol/L) and NEFA by enzymatic",
      "colorimetric assay (Wako Chemicals). Endogenous NiAc concentration",
      "was estimated to 6.8 nmol/L (population range 4.8-15 nmol/L) per",
      "Eq. 12 / Appendix A; baseline NEFA ranged 0.3-1.1 mmol/L. See",
      "Ahlstrom 2010 Methods (Animals and surgical procedures,",
      "Experimental design, Analytical assays) and Tables 1-2 for the",
      "fitted model output. Body weight was not used as a model covariate",
      "in the paper -- all PK parameters are reported on a per-kg basis,",
      "and the model file follows that parameterisation (per-kg amounts,",
      "per-kg volumes); the user can scale by absolute body weight",
      "downstream if needed."
    )
  )

  ini({
    # --------------------------------------------------------------
    # NiAc PK parameters -- Ahlstrom 2010 Table 1.
    # All structural parameters per kg of body weight.
    # RSE% from Table 1 are noted in trailing comments.
    # --------------------------------------------------------------
    lvmax_hi <- log(0.0573)  ; label("Vmax of high-affinity NiAc elimination, pathway 1 (umol/min/kg)")    # Table 1: Vmax1 = 0.0573 (RSE 1.57%)
    lkm_hi   <- log(0.00468) ; label("Km of high-affinity NiAc elimination, pathway 1 (umol/L)")           # Table 1: Km1 = 0.00468 (RSE 1.97%)
    lvmax_lo <- log(1.46)    ; label("Vmax of low-affinity NiAc elimination, pathway 2 (umol/min/kg)")     # Table 1: Vmax2 = 1.46 (RSE 4.20%)
    lkm_lo   <- log(16.6)    ; label("Km of low-affinity NiAc elimination, pathway 2 (umol/L)")            # Table 1: Km2 = 16.6 (RSE 6.20%)
    lvc      <- log(0.345)   ; label("Central volume of distribution Vc (L/kg)")                            # Table 1: Vc = 0.345 (RSE 0.559%)
    lq       <- log(0.0203)  ; label("Inter-compartmental clearance Q (paper Cld) (L/min/kg)")              # Table 1: Cld = 0.0203 (RSE 3.08%)
    lvp      <- log(3.54)    ; label("Peripheral volume of distribution Vp (paper Vt) (L/kg)")              # Table 1: Vt = 3.54 (RSE 7.66%)
    lsynt    <- log(0.0346)  ; label("Endogenous synthesis rate of NiAc Synt (umol/min/kg)")                # Table 1: Synt = 0.0346 (RSE 4.10%)

    # --------------------------------------------------------------
    # NEFA PD parameters -- Ahlstrom 2010 Table 2.
    # --------------------------------------------------------------
    lr0      <- log(0.606)   ; label("Baseline NEFA concentration R0 (mmol/L)")                              # Table 2: R0 = 0.606 (RSE 3.51%)
    lkout    <- log(0.411)   ; label("Fractional turnover rate of NEFA kout (L/mmol/min)")                   # Table 2: kout = 0.411 (RSE 9.22%)
    lktol    <- log(0.0267)  ; label("Turnover rate of moderator chain ktol (1/min)")                        # Table 2: ktol = 0.0267 (RSE 3.40%)
    lic50    <- log(0.0446)  ; label("Potency: NiAc Cp at 50% maximal NEFA-formation inhibition IC50 (umol/L)") # Table 2: IC50 = 0.0446 (RSE 7.20%)
    lgam     <- log(1.48)    ; label("Sigmoidicity (Hill exponent) of inhibitory drug function gamma (unitless)") # Table 2: gamma = 1.48 (RSE 3.79%)
    lpmod    <- log(1.21)    ; label("Amplification factor p on moderator M1 in NEFA formation (unitless)")  # Table 2: p = 1.21 (RSE 6.64%)
    lkcap    <- log(0.0318)  ; label("NiAc-independent capillary NEFA release rate kcap (mmol/L/min)")        # Table 2: kcap = 0.0318 (RSE 6.51%)
    imax     <- fixed(1.0)   ; label("Maximum drug-induced inhibition of NEFA formation Imax (fixed)")        # Table 2: Imax = 1 (fixed per Discussion)

    # --------------------------------------------------------------
    # IIV -- omega^2 = log(CV^2 + 1) conversion from Tables 1 and 2.
    # Only the parameters Tables 1/2 reported IIV for are encoded;
    # all others are population-typical (no eta).
    #   NiAc PK (Table 1):
    #     Vmax2:  CV 98.2% -> log(1 + 0.982^2) = 0.6754
    #     Cld:    CV 54.5% -> log(1 + 0.545^2) = 0.2606
    #     Synt:   CV 15.8% -> log(1 + 0.158^2) = 0.02472
    #   NEFA PD (Table 2):
    #     R0:     CV 29.0% -> log(1 + 0.290^2) = 0.08047
    #     kout:   CV 48.2% -> log(1 + 0.482^2) = 0.2089
    #     IC50:   CV 101%  -> log(1 + 1.01^2)  = 0.7039
    # --------------------------------------------------------------
    etalvmax_lo ~ 0.6754   # Table 1 IIV Vmax2 CV 98.2% (RSE 24.8)
    etalq       ~ 0.2606   # Table 1 IIV Cld   CV 54.5% (RSE 88.6)
    etalsynt    ~ 0.02472  # Table 1 IIV Synt  CV 15.8% (RSE 62.5)
    etalr0      ~ 0.08047  # Table 2 IIV R0    CV 29.0% (RSE 43.1)
    etalkout    ~ 0.2089   # Table 2 IIV kout  CV 48.2% (RSE 35.6)
    etalic50    ~ 0.7039   # Table 2 IIV IC50  CV 101%  (RSE 46.1)

    # --------------------------------------------------------------
    # Residual error.
    # NiAc PK (Table 1): proportional 34.4% + additive 0.800.
    #   Table 1 lists r2 with a "(%)" units header that appears to be a
    #   copy artefact from r1 (the value 0.800 cannot be a CV%; additive
    #   residual error has to carry concentration units). The most
    #   plausible reading consistent with the NONMEM data scale (NiAc
    #   plasma concentrations are reported in umol/L throughout the
    #   paper, range 0.001-20 umol/L) is that 0.800 is the additive SD
    #   in umol/L. Documented as an interpretation in the vignette
    #   Assumptions and deviations section. Encoded literally as 0.800.
    # NEFA PD (Table 2): proportional only (additive was originally
    #   modelled but estimated as negligible -- only proportional retained).
    # --------------------------------------------------------------
    propSd      <- 0.344 ; label("Proportional residual error on NiAc plasma concentration (fraction)")     # Table 1: r1 = 34.4% (RSE 18.6)
    addSd       <- 0.800 ; label("Additive residual SD on NiAc plasma concentration (umol/L; unit inferred)") # Table 1: r2 = 0.800 (RSE 5.80); see vignette Assumptions
    propSd_NEFA <- 0.221 ; label("Proportional residual error on NEFA plasma concentration (fraction)")     # Table 2: r = 22.1% (RSE 8.51)
  })

  model({
    # --------------------------------------------------------------
    # 1. Individual structural parameters.
    # Exponential IIV on log-transformed typical values per Methods:
    # "Interindividual variability was modeled as exponential models
    # for all disposition parameters of NiAc and NEFA".
    # --------------------------------------------------------------
    vmax_hi <- exp(lvmax_hi)
    km_hi   <- exp(lkm_hi)
    vmax_lo <- exp(lvmax_lo + etalvmax_lo)
    km_lo   <- exp(lkm_lo)
    vc      <- exp(lvc)
    q       <- exp(lq + etalq)
    vp      <- exp(lvp)
    synt    <- exp(lsynt + etalsynt)

    r0      <- exp(lr0 + etalr0)
    kout    <- exp(lkout + etalkout)
    ktol    <- exp(lktol)
    ic50    <- exp(lic50 + etalic50)
    gam     <- exp(lgam)
    pmod    <- exp(lpmod)
    kcap    <- exp(lkcap)

    # --------------------------------------------------------------
    # 2. NiAc disposition: 2-cmt with 2 parallel MM elimination plus
    #    endogenous zero-order synthesis (Ahlstrom 2010 Eqs. 1a-1b).
    # Per-kg formulation: state-amounts in umol/kg, volumes in L/kg,
    # concentrations in umol/L; dose amount is umol/kg.
    # --------------------------------------------------------------
    Cp <- central     / vc                                          # umol/L
    Ct <- peripheral1 / vp                                          # umol/L

    d/dt(central)     <- synt -
                         vmax_hi * Cp / (km_hi + Cp) -
                         vmax_lo * Cp / (km_lo + Cp) -
                         q * (Cp - Ct)
    d/dt(peripheral1) <- q * (Cp - Ct)

    # --------------------------------------------------------------
    # 3. Endogenous-baseline initialisation of NiAc compartments.
    # At baseline (Inf = 0, dCp/dt = dCt/dt = 0) the central-compartment
    # balance reduces to a quadratic in Cp_baseline (paper Appendix A
    # Eq. 12). With A*Cp^2 + B*Cp + C = 0 and
    #   A = synt - vmax_hi - vmax_lo               (negative under physiological values)
    #   B = synt*(km_hi + km_lo) - vmax_hi*km_lo - vmax_lo*km_hi
    #   C = synt * km_hi * km_lo
    # the unique positive root is (-B - sqrt(B^2 - 4AC)) / (2A).
    # For the typical-value parameter set this gives Cp_baseline
    # ~ 0.0068 umol/L (= 6.8 nmol/L), matching Ahlstrom 2010 Results.
    # Per-individual draws of (synt, vmax_lo) produce per-individual
    # baselines spanning ~4.8-15 nmol/L (Results, NiAc disposition
    # section), reproducing the paper's population baseline range.
    # --------------------------------------------------------------
    Acoef <- synt - vmax_hi - vmax_lo
    Bcoef <- synt * (km_hi + km_lo) - vmax_hi * km_lo - vmax_lo * km_hi
    Ccoef <- synt * km_hi * km_lo
    Cp0   <- (-Bcoef - sqrt(Bcoef * Bcoef - 4 * Acoef * Ccoef)) / (2 * Acoef)

    central(0)     <- Cp0 * vc
    peripheral1(0) <- Cp0 * vp

    # --------------------------------------------------------------
    # 4. NEFA feedback model (Ahlstrom 2010 Eqs. 2-4 + N=8 chain).
    #
    # Inhibitory drug-mechanism function (Eq. 2):
    #   I(Cp) = Imax * Cp^gamma / (IC50^gamma + Cp^gamma)
    #
    # NEFA balance (Eq. 3, with kin solved from steady-state R = M_i = R0):
    #   d/dt(NEFA) = kin*(1 - I(Cp)) / M1^p + kcap - kout * M8 * NEFA
    #   kin        = R0^p * (kout*R0^2 - kcap)
    #
    # Moderator feedback chain (Eq. 4 generalised to N = 8 per Results
    # "the optimal number of moderator transit compartments N was
    # estimated to be 8"):
    #   d/dt(M1) = ktol*(NEFA - M1)
    #   d/dt(Mi) = ktol*(M(i-1) - Mi)   for i = 2,...,8
    #
    # nlmixr2lib mapping: the moderator chain is encoded with the
    # blessed precursor<n> chain prefix (closest existing convention
    # for a first-order delay cascade representing biosignal
    # transduction -- see Friberg 2002 paclitaxel for the precedent).
    # The nefa compartment is registered in R/conventions.R alongside
    # this model.
    # --------------------------------------------------------------
    inh <- imax * Cp^gam / (ic50^gam + Cp^gam)
    kin <- r0^pmod * (kout * r0 * r0 - kcap)

    d/dt(nefa)       <- kin * (1 - inh) / precursor1^pmod +
                        kcap - kout * precursor8 * nefa
    d/dt(precursor1) <- ktol * (nefa       - precursor1)
    d/dt(precursor2) <- ktol * (precursor1 - precursor2)
    d/dt(precursor3) <- ktol * (precursor2 - precursor3)
    d/dt(precursor4) <- ktol * (precursor3 - precursor4)
    d/dt(precursor5) <- ktol * (precursor4 - precursor5)
    d/dt(precursor6) <- ktol * (precursor5 - precursor6)
    d/dt(precursor7) <- ktol * (precursor6 - precursor7)
    d/dt(precursor8) <- ktol * (precursor7 - precursor8)

    nefa(0)       <- r0
    precursor1(0) <- r0
    precursor2(0) <- r0
    precursor3(0) <- r0
    precursor4(0) <- r0
    precursor5(0) <- r0
    precursor6(0) <- r0
    precursor7(0) <- r0
    precursor8(0) <- r0

    # --------------------------------------------------------------
    # 5. Observations.
    # Cc   -- NiAc plasma concentration (umol/L), canonical PK output.
    # NEFA -- response (mmol/L), uppercase alias for the nefa state to
    #         match the paper's published symbol.
    # --------------------------------------------------------------
    Cc   <- Cp
    NEFA <- nefa

    Cc   ~ prop(propSd) + add(addSd)
    NEFA ~ prop(propSd_NEFA)
  })
}
