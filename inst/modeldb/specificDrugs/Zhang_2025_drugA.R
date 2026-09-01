Zhang_2025_drugA <- function() {
  description <- paste(
    "Semi-mechanistic dose-exposure-response (D-E-R) model for the",
    "anti-psoriatic 'Drug A' of Zhang 2025, linking subcutaneous PK to the",
    "Psoriasis Area and Severity Index (PASI) score, with a covariate effect",
    "for a region of interest in a multi-regional Phase 2 trial.",
    "PK: two-compartment disposition with first-order SC absorption and",
    "parallel linear (kel 0.0724 /day) and Michaelis-Menten (Vmax 7.4 mg/day,",
    "Km 0.01 ug/mL) elimination from a 3.9 L central compartment; F 57.6%.",
    "PD: central drug drives an effect ('signaling') compartment whose",
    "concentration inhibits plaque synthesis in an indirect-response PASI",
    "model, ksyn = ksyn0 * (1 - Imax * Ce / (IC50 + Ce)) with Imax fixed at 1",
    "and IC50 2.86 ug/mL; PASI starts at its baseline 13.5 and degrades at",
    "kdeg = ksyn0 / baseline. A placebo effect MULTIPLIES the degradation",
    "rate, kdegx = kdeg * (1 + Pmax * exp(-kpl * t)), and is maximal at t = 0,",
    "so the placebo arm improves to about an 18% PASI reduction near week 2-3",
    "and then recovers (Zhang 2025 Figure 1B dashed lines).",
    "The inter-regional difference is a multiplicative effect on IC50,",
    "IC50 = IC50 * R_regionX^REGION_CHINA; R_regionX = 2.6 is the paper's",
    "borderline scenario, at which Region X patients reach a 50% median PASI",
    "reduction at week 12 on 210 mg versus 68.4% in typical patients.",
    "'Drug A' is anonymized in the article; the D-E-R model is cited to",
    "Salinger 2014 and Papp 2012, both of which are brodalumab (AMG 827),",
    "and the parameter values are the authors' lightly adjusted set for the",
    "planned Phase 2 trial rather than a brodalumab fit -- so the model is",
    "registered under the paper's own compound label."
  )
  reference <- paste(
    "Zhang X, Xiao Y, Wu J, Marshall S, Zhou X.",
    "Pharmacometric Model-Based Sample Size Allocation for a Region of",
    "Interest in a Multi-Regional Phase 2 Trial: A Case Study of an",
    "Anti-Psoriatic Drug.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1673-1682.",
    "doi:10.1002/psp4.70090.",
    "Parameter values from Table S1; structural equations and variance terms",
    "from the Data S4 NONMEM control streams.",
    sep = " "
  )
  vignette <- "Zhang_2025_drugA"

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot = list(
      analyte = "drug A (anti-IL-17 receptor monoclonal antibody; brodalumab per the cited source models)",
      units = "mg", specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "drug A", units = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "drug A", units = "mg", specimen = "plasma", verified = TRUE
    ),
    pasi = list(
      analyte = "none", units = "PASI units (0-72 clinical score)",
      specimen = "not applicable", verified = TRUE
    ),
    effect = list(
      analyte = "drug A", units = "mg", specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    REGION_CHINA = list(
      description        = "1 = subject enrolled at a study site in the region of interest, 0 = enrolled elsewhere.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Region X / rest of world)",
      notes              = paste(
        "The published article refers to this region only as 'Region X' and",
        "stresses the generality of the method. The region is identified as",
        "China by the authors' own deposited Data S4 R code comments",
        "('sample size for cn (Region X)' and 'add label for Chinese and ROW",
        "populations (POP=2 is Region X, POP=1 is non-Region X)'), not by the",
        "article body. Encoded from the control stream's POP column as",
        "REGION_CHINA = POP - 1, matching COV_CN**(POP-1) in $PK."
      ),
      source_name        = "POP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 175L,
    n_studies      = 1L,
    age_range      = "not reported (adult psoriasis population)",
    weight_range   = "not reported",
    disease_state  = "moderate-to-severe plaque psoriasis",
    dose_range     = paste(
      "Drug A 70, 140 or 210 mg SC at day 1 and weeks 1, 2, 4, 6, 8 and 10;",
      "or 280 mg SC at day 1 and weeks 4 and 8; or placebo.",
      "Randomization 1:1:1:1:1."
    ),
    regions        = "multi-regional; a designated 'Region X' (China, per the Data S4 code comments) versus the rest of the world",
    notes          = paste(
      "This is a SIMULATION model, not a fit to the Phase 2 data. The N = 175",
      "is the planned size of the hypothetical multi-regional Phase 2",
      "dose-ranging trial that Zhang 2025 simulates (Methods 2.1); the",
      "underlying D-E-R model was developed on Phase 1 data for Drug A",
      "(Zhang 2025 refs 10-11) and carried into this paper 'with minor",
      "adjustments to parameter values' (Methods 2.2). No fitted-cohort",
      "demographics are reported. The primary endpoint is percentage",
      "improvement from baseline in PASI at week 12; PASI evaluation visits",
      "were day 0 and weeks 1, 2, 4, 6, 8, 10, 12 and 16, and PK sampling",
      "was at days 7, 8, 11, 14, 28, 56, 60, 64, 70, 84, 112 and 154",
      "(Data S4). NONMEM 7.5.0 (ADVAN13) with PsN 5.3.1 and R 4.3.2."
    )
  )

  ini({
    # ================================================================
    # Pharmacokinetics. Data S4 "NONMEM model code for PK simulation":
    # every $THETA, $OMEGA and $SIGMA in that stream carries FIX,
    # because the PK layer was fixed from the Phase 1 Drug A model and
    # simulated sequentially ahead of the PD layer -- hence fixed()
    # throughout this block. Values cross-checked against Table S1.
    # ================================================================
    lkel <- fixed(log(0.0724))
    label("First-order elimination rate constant from central (1/day)")
    # Table S1 Kel = 0.0724 /day; Data S4 PK $THETA(1) 0.0724 FIX.

    lvc <- fixed(log(3.9))
    label("Central compartment volume (L)")
    # Table S1 V = 3.9 L; Data S4 PK $THETA(2) 3.9 FIX. (The control
    # stream's inline comment says "V2 mL"; Table S1's L is the
    # dimensionally consistent unit -- dose mg / 3.9 L = ug/mL.)

    lvmax <- fixed(log(7.4))
    label("Maximum rate of Michaelis-Menten elimination (mg/day)")
    # Table S1 Vmax = 7.4 mg/day; Data S4 PK $THETA(3) 7.4 FIX.

    lkm <- fixed(log(0.01))
    label("Concentration at half-maximal nonlinear elimination (ug/mL)")
    # Table S1 Km = 0.01 ug/mL; Data S4 PK $THETA(4) 0.01 FIX.

    lka <- fixed(log(0.255))
    label("First-order absorption rate constant (1/day)")
    # Table S1 KA = 0.255 /day; Data S4 PK $THETA(5) 0.255 FIX.
    # (The stream's inline comment says "ka 1/hr"; Table S1 and the
    # day-based dosing/sampling grid both make it 1/day.)

    lk12 <- fixed(log(0.26))
    label("Central-to-peripheral rate constant (1/day)")
    # Table S1 K12 = 0.26 /day; Data S4 PK $THETA(6) 0.26 FIX.

    lk21 <- fixed(log(0.351))
    label("Peripheral-to-central rate constant (1/day)")
    # Table S1 K21 = 0.351 /day; Data S4 PK $THETA(7) 0.351 FIX.

    lfdepot <- fixed(log(0.576))
    label("Subcutaneous bioavailability (fraction)")
    # Table S1 F1 = 57.6%; Data S4 PK $THETA(8) 0.576 FIX.

    # ================================================================
    # Pharmacodynamics: indirect-response PASI model driven by an
    # effect ("signaling") compartment. Data S4 "NONMEM model code for
    # simulation of treatment effects based on post-hoc PK parameters".
    # These $THETA entries carry no FIX except Imax and the regional
    # effect, so only those two are wrapped in fixed().
    # ================================================================
    lrbase <- log(13.5)
    label("Baseline PASI score (PASI units)")
    # Table S1 BL = 13.5 PASI units; Data S4 PD $THETA(1) 13.5.

    lksyn <- log(0.862)
    label("Zero-order rate of plaque formation, ksyn0 (PASI units/day)")
    # Table S1 Ksyn = 0.862 PASI units/day; Data S4 PD $THETA(2) 0.862.

    limax <- fixed(log(1))
    label("Maximum fractional inhibition of plaque formation (unitless)")
    # Table S1 Imax = 1; Data S4 PD $THETA(3) "1 FIX".

    lic50 <- log(2.86)
    label("Effect-compartment concentration at half-maximal inhibition (ug/mL)")
    # Table S1 IC50 = 2.86 ug/mL; Data S4 PD $THETA(4) 2.86.

    lpmax <- log(0.439)
    label("Maximum fractional increase of the PASI degradation rate by placebo (unitless)")
    # Table S1 PBOmax = 0.439; Data S4 PD $THETA(5) 0.439. Table S1
    # labels this "1/day", but it enters as the dimensionless
    # multiplier PBOMX in KDEGX = KDEG * (1 + PBOMX * EXP(-KPBO*T)).

    lkpl <- log(0.0372)
    label("Decay rate constant of the placebo effect (1/day)")
    # Table S1 KPBO = 0.0372 /day; Data S4 PD $THETA(6) 0.0372.

    lke0 <- log(0.0469)
    label("Effect-compartment (signaling) equilibration rate constant (1/day)")
    # Table S1 Keo = 0.0469 /day; Data S4 PD $THETA(7) 0.0469 (K50).

    e_region_china_ic50 <- fixed(2.6)
    label("Ratio of IC50 in Region X patients to IC50 in typical patients (R_regionX, unitless)")
    # Data S4 PD $THETA(8) COV_CN, entered as COV_CN**(POP-1) and set to
    # the scenario value via "COV_CN_VALUE FIX"; the MCMP generating
    # script sets COV_CN_VALUE_PD <- 2.6. Zhang 2025 Table 1 and the
    # Figure 1 caption give 2.6 as the BORDERLINE clinically relevant
    # difference (Region X reaches exactly a 50% median PASI reduction
    # at week 12 on 210 mg, versus 68.4% for typical patients). Table 1
    # maps the whole range: 1.0 = no difference, 1.0-2.6 = inferior but
    # still clinically meaningful, > 2.6 = clinically relevant
    # difference. The paper's SSE scenarios sweep 1.0 to 5.4.

    # ================================================================
    # Between-subject variability. Exponential on all six parameters;
    # the $OMEGA entries are variances and each square-roots to the
    # %CV printed in Table S1. The three PK omegas carry FIX in the PK
    # stream (fixed from the Phase 1 model); the three PD omegas do not.
    # ================================================================
    etalvc ~ fixed(0.0841)
    # Data S4 PK $OMEGA(1) 0.0841 FIX; sqrt = 0.290 = Table S1 omega(V) 29.0% CV.

    etalka ~ fixed(0.564)
    # Data S4 PK $OMEGA(2) 0.564 FIX; sqrt = 0.751 = Table S1 omega(KA) 75.1% CV.

    etalvmax ~ fixed(0.0992)
    # Data S4 PK $OMEGA(3) 0.0992 FIX; sqrt = 0.315 = Table S1 omega(Vmax) 31.5% CV.

    etalrbase ~ 0.046656
    # Data S4 PD $OMEGA(1) 0.046656; sqrt = 0.216 = Table S1 omega(BL) 21.6% CV.

    etalic50 ~ 1.8496
    # Data S4 PD $OMEGA(2) 1.8496; sqrt = 1.36 = Table S1 omega(IC50) 136% CV.

    etalke0 ~ 0.327184
    # Data S4 PD $OMEGA(3) 0.327184; sqrt = 0.572 = Table S1 omega(Keo) 57.2% CV.

    # ================================================================
    # Residual error. Both endpoints use Y = IPRED*(1+EPS1) + EPS2, so
    # combined proportional + additive. The $SIGMA entries are
    # VARIANCES; the SDs below are their square roots.
    #
    # ERRATUM: Table S1 reports the PK additive residual as
    # "SD (ug/mL) 0.09", but 0.09 is the $SIGMA VARIANCE, so the SD is
    # sqrt(0.09) = 0.3 ug/mL. The other three rows prove the rule:
    # sqrt(0.2683) = 0.518 = Table S1's stated 51.8% CV,
    # sqrt(0.00404496) = 0.0636 = Table S1's stated 6.36%, and
    # sqrt(1.31) = 1.145 = Table S1's stated 1.145 PASI units.
    # Three consistent rows falsify the fourth; the control stream wins.
    # ================================================================
    propSd <- fixed(0.518)
    label("Proportional residual error for drug concentration (fraction)")
    # Data S4 PK $SIGMA(1) 0.2683 FIX; sqrt = 0.518 = Table S1 51.8% CV.

    addSd <- fixed(0.3)
    label("Additive residual error for drug concentration (ug/mL)")
    # Data S4 PK $SIGMA(2) 0.09 FIX; sqrt = 0.3 ug/mL. See the erratum note above.

    propSd_pasi <- 0.0636
    label("Proportional residual error for PASI score (fraction)")
    # Data S4 PD $SIGMA(1) 0.00404496; sqrt = 0.0636 = Table S1 6.36% CV.

    addSd_pasi <- 1.145
    label("Additive residual error for PASI score (PASI units)")
    # Data S4 PD $SIGMA(2) 1.31; sqrt = 1.145 = Table S1 1.145 PASI units.
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters (Data S4 PK $PK block).
    # ------------------------------------------------------------
    kel  <- exp(lkel)
    vc   <- exp(lvc + etalvc)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)
    ka   <- exp(lka + etalka)
    k12  <- exp(lk12)
    k21  <- exp(lk21)

    # ------------------------------------------------------------
    # Individual PD parameters (Data S4 PD $PK block).
    # The regional effect multiplies IC50 in the power form
    # COV_CN**(POP-1), i.e. IC50 * R_regionX^REGION_CHINA with
    # REGION_CHINA = POP - 1 in {0, 1}.
    # ------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    ksyn0 <- exp(lksyn)
    imax  <- exp(limax)
    ic50  <- exp(lic50 + etalic50) * e_region_china_ic50^REGION_CHINA
    ke0   <- exp(lke0 + etalke0)
    pmax  <- exp(lpmax)
    kpl   <- exp(lkpl)

    # PASI degradation is pinned so that the drug- and placebo-free
    # system sits exactly at baseline: KDEG = KSYN0 / BL.
    kdeg <- ksyn0 / rbase

    # Effect-compartment volume, V5 = K25 * V2 / K50 with K25 = K
    # (Data S4 PD $PK). This scaling makes the effect-compartment
    # concentration equal the central concentration at equilibrium.
    ve <- kel * vc / ke0

    Cc <- central / vc
    Ce <- effect / ve

    # Imax inhibition of plaque formation, and the placebo effect as a
    # multiplicative increase of the degradation rate that is maximal
    # at t = 0 and decays with kpl. t is absolute time, matching the
    # NONMEM $DES variable T. The stream's "IF (TRT.EQ.0) KSYN=KSYN0"
    # line is not needed here: it only exists because the placebo arm
    # was dosed with AMT = 1e-10 mg to make NONMEM generate individual
    # PK parameters. With a genuinely undosed placebo arm, Ce is
    # identically zero and ksyn reduces to ksyn0 on its own.
    ksyn  <- ksyn0 * (1 - imax * Ce / (ic50 + Ce))
    kdegx <- kdeg * (1 + pmax * exp(-kpl * t))

    # ------------------------------------------------------------
    # ODE system (Data S4 $DES). State order matches the control
    # stream's $MODEL: DEPOT, CENTRAL, PERIPH, PASI, SIGCMP.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - vmax * Cc / (km + Cc) -
                           kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pasi)        <- ksyn - kdegx * pasi
    d/dt(effect)      <- kel * central - ke0 * effect

    pasi(0) <- rbase

    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------
    # Observations. The PK stream observes A(2)/V2; the PD stream
    # observes A(4) itself (COMP=(PASI,DEFOBS), IPRED=A(4)).
    # ------------------------------------------------------------
    Cc   ~ add(addSd) + prop(propSd)
    pasi ~ add(addSd_pasi) + prop(propSd_pasi)
  })
}
