Ande_2018_paclitaxel_everolimus_dasatinib <- function() {
  description <- "In vitro (JIMT-1 HER2-positive trastuzumab-resistant breast cancer cell line) three-dimensional and dynamic (3DD) BelloCell-bioreactor PK/PD model for the sequential triple combination of paclitaxel (PAC) + everolimus (EVE) + dasatinib (DAS). PAC follows a two-compartment mammillary PK model (CL/Q/Vc=Vp_paper/Vp=Vt_paper) driven by a 3-h IV infusion. EVE and DAS share a single combined concentration covariate (their molar concentrations are equal throughout the experiment) constructed inside model() as a step-function-plus-exponential-decline: 0 nM for t<24 h, 50 nM during 24<=t<96 h, then exponential decline with rate kde from t>=96 h derived from the BelloCell perfusion mechanics (0.53 mL/min / 500 mL). PAC and DE concentrations jointly stimulate caspase-3 production via two paper-named slope parameters (spac, sde) entering a 5-state Lobo-Balthasar signal-transduction chain (4 transit relays + 1 active caspase-3 pool) with a Mager-Jusko-style feedback loop in which the production rate of the first transit is divided by the active caspase-3 level. The change in active caspase-3 from baseline drives the kill term in a modified Simeoni (2004) tumor-growth model (exponential-then-linear growth with psi=20 FIXED) on JIMT-1 cell count. No IIV or residual error is reported in the paper (Table 3 reports % RSE on point estimates only); a small fixed additive residual error is included per output as a placeholder so the model is fittable, with values flagged in the vignette Assumptions and deviations section."
  reference <- paste(
    "Ande A, Vaidya TR, Tran BN, Vicchiarelli M, Brown AN, Ait-Oudhia S (2018).",
    "Utility of a Novel Three-Dimensional and Dynamic (3DD) Cell Culture System",
    "for PK/PD Studies: Evaluation of a Triple Combination Therapy at Overcoming",
    "Anti-HER2 Treatment Resistance in Breast Cancer.",
    "Frontiers in Pharmacology 9:403.",
    "doi:10.3389/fphar.2018.00403.",
    sep = " "
  )
  vignette <- "Ande_2018_paclitaxel_everolimus_dasatinib"
  units <- list(
    time          = "hour",
    dosing        = "ug (PAC IV infusion); 50 nM target bath concentration (DE)",
    concentration = "nmol/L (PAC PD driver Cc); fold/L (caspase3Act, unitless relative-to-baseline); cells/L (tumorCells, total bioreactor count)"
  )

  covariateData <- list()

  population <- list(
    species             = "in vitro (JIMT-1 HER2-positive trastuzumab-resistant breast cancer cell line)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Homo sapiens JIMT-1 cell line (HER2-positive, MUC4-expressing, refractory to trastuzumab, pertuzumab, T-DM1 and lapatinib; AddexBio)",
    system              = "BelloCell HD continuous cell culture system (Cesco Bioengineering) with 865 BioNOC II PET fabric carriers (bed volume 100 cm^3); two BelloCell-500AP perfusion bottles per arm (one control, one treatment); 500 mL working volume; oscillating up-down speed 1.5 mm/s upper-hold 2 min bottom-hold 1.5 min; 0.53 mL/min perfusion in/out to maintain isovolumetric medium replacement.",
    medium              = "DMEM + 5% FBS + 1% penicillin-streptomycin; 1% DMSO maintained throughout the 5-day treatment period; assay medium with and without drugs made fresh and added daily.",
    temperature         = "37 C, 5% CO2",
    duration            = "21 days (5-day treatment of PAC + DAS + EVE followed by 2-day washout and 14-day regrowth observation); cell counts every 2 days; caspase-3 measured daily for the first 5 days; PAC PK sampling at 2, 3, 6, 24, 27, 30, 48, 56, 72, 78, 96 h.",
    inoculum            = "30 mL of JIMT-1 cell suspension at 2.6e6 cells/mL (8e7 cells per bottle) seeded onto BioNOC II carriers; 72 h pre-attachment incubation before PAC infusion at experiment t = 0; observed initial post-attachment count approximately 1.5e8 cells per Figure 6B starting points.",
    regimens            = "Treatment arm: PAC 3.75 ug/min (225.14 ug/h) constant-rate IV infusion into the bioreactor for 3 h starting at t = 0; DAS and EVE spiked at t = 24 h to 50 nM each and maintained at 50 nM via continuous perfusion infusion for the following 72 h (DE 'on' during 24 <= t < 96 h), then washout via daily medium replacement.",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. The paper fits in Monolix and reports point estimates with % RSE for each parameter (Tables 1 / 2 / 3) but does NOT tabulate a residual-error SD, between-bottle eta, or assay precision. The packaged 3DD model accordingly contains no etas; per-output additive residual SDs are encoded as small placeholders (flagged with fixed() and operator-derived inline comments) so the model parses but should not be interpreted as Monolix-reported uncertainty. See Ande 2018 Methods (3DD Cell Culture Experiments, Drugs Dosing Schedule, PK/PD Modeling of Data From 3DD Cell Culture Settings) and Table 3."
  )

  ini({
    # ===================================================================
    # PACLITAXEL PHARMACOKINETICS
    # ===================================================================
    # 2-compartment mammillary model. Paper writes the central volume as
    # Vp (Cp = central concentration) and the peripheral volume as Vt
    # (Ct = peripheral concentration); the canonical mapping is
    # vc <- paper Vp and vp <- paper Vt (per parameter-names.md, V1 / V2
    # always map to vc / vp regardless of the paper's local naming).
    lcl <- log(0.541)
    label("Paclitaxel clearance CL (L/h)")  # Ande 2018 Table 3
    lvc <- log(9.8)
    label("Paclitaxel central volume Vc (L; paper's Vp)")  # Ande 2018 Table 3
    lq  <- log(0.998)
    label("Paclitaxel inter-compartmental clearance Q (L/h)")  # Ande 2018 Table 3
    lvp <- log(25)
    label("Paclitaxel peripheral volume Vp (L; paper's Vt)")  # Ande 2018 Table 3

    # ===================================================================
    # DE (DAS + EVE) BIOREACTOR-WASHOUT KINETICS
    # ===================================================================
    # The paper does not tabulate the EVE/DAS washout decay rate; the
    # profile is described as 'a simple step-function followed by an
    # exponential decline' simulated from BelloCell mass balance
    # (Methods, PK/PD Modeling of Data From 3DD Cell Culture Settings).
    # The perfusion rate is 0.53 mL/min over a 500 mL working volume,
    # giving a first-order washout of 0.53 / 500 * 60 = 0.0636 /h. Held
    # FIXED at this Methods-derived value.
    lkde <- fixed(log(0.0636))
    label("DE bioreactor washout rate kde (1/h; FIXED from BelloCell perfusion mechanics)")  # Ande 2018 Methods (3DD Cell Culture System Calibration)

    # ===================================================================
    # CASPASE-3 SIGNAL-TRANSDUCTION CHAIN
    # ===================================================================
    # Lobo-Balthasar 5-state transit-compartment model with Mager-Jusko
    # feedback. Kin / Kout / Ktr collapse to a single rate constant Ktr
    # because the model is parameterised so the no-treatment baseline
    # caspase-3 level is unity (Methods: 'In control... Kin and Kout
    # become numerically equivalent to Ktr. Hence, a single rate constant
    # Ktr was estimated.').
    lktr <- log(0.146)
    label("Caspase-3 signal-transduction transit rate constant Ktr (1/h)")  # Ande 2018 Table 3

    # Slope parameters quantifying drug effect on caspase-3 production.
    # Paper uses S_PAC and S_DE; canonical rendering is lspac / lsde
    # (lowercase paper-symbol + log prefix). The 'PAC' driver is the
    # central-compartment PAC concentration in nM; the 'DE' driver is
    # the bath EVE / DAS concentration (which is equal throughout the
    # experiment, per paper Methods 'CONC_DE represents the concentration
    # of DAS or EVE, which is equivalent throughout the duration of the
    # study').
    lspac <- log(0.312)
    label("Paclitaxel slope on caspase-3 production S_PAC (L/nmol)")  # Ande 2018 Table 3
    lsde  <- log(0.306)
    label("Dasatinib+everolimus combined slope on caspase-3 production S_DE (L/nmol)")  # Ande 2018 Table 3

    # ===================================================================
    # MODIFIED SIMEONI TUMOR GROWTH (JIMT-1 cell count)
    # ===================================================================
    # Simeoni 2004 exponential-to-linear growth: dR/dt =
    # lambda0*R / (1 + (lambda0*R/lambda1)^psi)^(1/psi) - kd*deltaC3*R.
    # psi controls the smoothness of the exponential -> linear switch
    # (psi = 20 makes it nearly a hard threshold).
    llambda0 <- log(0.0077)
    label("Exponential tumor-growth rate constant lambda0 (1/h)")  # Ande 2018 Table 3
    llambda1 <- log(7.3)
    label("Linear tumor-growth rate constant lambda1 (cells / h)")  # Ande 2018 Table 3
    lpsi     <- fixed(log(20))
    label("Simeoni exponential-to-linear transition exponent psi (unitless; FIXED)")  # Ande 2018 Table 3
    lkd      <- log(0.0096)
    label("JIMT-1 cell-death rate constant kd (cells / h per unit caspase-3 change)")  # Ande 2018 Table 3

    # Initial JIMT-1 cell count at the start of PK/PD experiment.
    # Paper Methods seeds 8e7 cells per BelloCell bottle; after the 72-h
    # pre-attachment / growth period the cell number visible at the
    # start of observation (Figure 6B initial point) is approximately
    # 1.51e8 (treatment arm) and 1.557e8 (control arm). Not reported in
    # any table -- operator-derived from Figure 6B initial points.
    # operator-derived from Figure 6B initial cell count (~151 million cells)
    lr0 <- log(1.51e8)
    label("Initial JIMT-1 cell count at start of PK/PD experiment R0 (cells)")  # Ande 2018 Figure 6B (initial point of treatment arm)

    # ===================================================================
    # PACLITAXEL MOLAR-MASS CONVERSION
    # ===================================================================
    # PAC central-compartment amount is tracked in ug; the PD slope
    # parameter S_PAC is in L/nmol so the PD driver Cc must be in nM.
    # The conversion factor is pac_mw_ug_per_nmol = PAC MW (g/mol) /
    # 1000 = 853.91 / 1000 = 0.85391 ug/nmol. Held FIXED at the
    # tabulated PAC molar mass (PubChem CID 36314).
    lpac_mw <- fixed(log(0.85391))
    label("Paclitaxel molar-mass conversion (ug per nmol; FIXED literature)")  # PubChem CID 36314 (paclitaxel MW 853.91 g/mol)

    # ===================================================================
    # RESIDUAL ERROR (PLACEHOLDERS)
    # ===================================================================
    # Ande 2018 reports % RSE on point estimates only; no residual SD
    # is tabulated for PAC concentration, caspase-3 activity, or cell
    # count. The values below are small operator-chosen placeholders so
    # the multi-output model parses cleanly; see vignette Assumptions
    # and deviations.
    # operator-chosen placeholder (paper reports no residual SD)
    addSd <- fixed(1)
    label("Additive residual SD on paclitaxel concentration Cc (nM; FIXED placeholder)")
    # operator-chosen placeholder
    addSd_caspase3Act <- fixed(0.1)
    label("Additive residual SD on caspase-3 activity (relative-to-baseline units; FIXED placeholder)")
    # operator-chosen placeholder
    addSd_tumorCells <- fixed(5e6)
    label("Additive residual SD on JIMT-1 cell count (cells; FIXED placeholder)")
  })

  model({
    # ===================================================================
    # 1. Typical-value parameters (no IIV)
    # ===================================================================
    cl       <- exp(lcl)
    vc       <- exp(lvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    kde      <- exp(lkde)
    ktr      <- exp(lktr)
    spac     <- exp(lspac)
    sde      <- exp(lsde)
    lambda0  <- exp(llambda0)
    lambda1  <- exp(llambda1)
    psi      <- exp(lpsi)
    kd       <- exp(lkd)
    r0       <- exp(lr0)
    pac_mw   <- exp(lpac_mw)

    # ===================================================================
    # 2. Paclitaxel PK micro-constants (2-cmt mammillary)
    # ===================================================================
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ===================================================================
    # 3. PAC central / peripheral ODEs and PD driver
    # ===================================================================
    # central holds the amount in ug; Cc (the PD driver) is converted
    # to nM by dividing the per-volume concentration by pac_mw (ug per
    # nmol). PAC dosing is a 3-h IV infusion delivered via the user
    # event table (amt = 675.42, rate = 225.14, cmt = central).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    Cc <- (central / vc) / pac_mw

    # ===================================================================
    # 4. DE bath concentration (paper-described step-plus-exponential)
    # ===================================================================
    # Per Ande 2018 Methods (Drugs Dosing Schedule + PK/PD Modeling of
    # Data From 3DD Cell Culture Settings):
    #   conc_de = 0       for 0  <= t < 24
    #   conc_de = 50 nM   for 24 <= t < 96
    #   conc_de = 50 nM * exp(-kde * (t - 96))   for t >= 96
    # Encoded via indicator multiplication so rxode2 evaluates a single
    # continuous expression at every solver step (no piecewise branch
    # discontinuities).
    ind_dose  <- (time >= 24.0) * (time < 96.0)
    ind_post  <- (time >= 96.0)
    conc_de   <- 50 * ind_dose + 50 * ind_post * exp(-kde * (time - 96.0))

    # ===================================================================
    # 5. Caspase-3 signal-transduction chain (Ande 2018 Eqs 14-18)
    # ===================================================================
    # Stimulation: 1 + spac * Cc + sde * conc_de.
    # First transit has Mager-Jusko feedback: production rate is
    # (ktr / transit5) * stim. Subsequent relays propagate at rate ktr.
    # The 'fifth' state (transit5) accumulates with Kout = ktr (= Kin =
    # Ktr at baseline) so the no-drug steady-state level of every state
    # is unity. delta_casp = transit5 - 1 is the drug-perturbation
    # signal that drives tumor-cell death (zero at baseline).
    stim <- 1 + spac * Cc + sde * conc_de
    d/dt(transit1) <-  (ktr / transit5) * stim - ktr * transit1
    d/dt(transit2) <-  ktr * (transit1 - transit2)
    d/dt(transit3) <-  ktr * (transit2 - transit3)
    d/dt(transit4) <-  ktr * (transit3 - transit4)
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5

    # ===================================================================
    # 6. Modified Simeoni JIMT-1 tumor growth with caspase-driven death
    # ===================================================================
    # dR/dt = lambda0*R / (1 + (lambda0*R/lambda1)^psi)^(1/psi)
    #         - kd * (transit5 - 1) * R
    delta_casp <- transit5 - 1
    grow_num   <- lambda0 * tumor
    grow_den   <- (1 + (lambda0 * tumor / lambda1) ^ psi) ^ (1 / psi)
    d/dt(tumor) <- grow_num / grow_den - kd * delta_casp * tumor

    # ===================================================================
    # 7. Initial conditions (no-treatment baseline = 1 for caspase chain)
    # ===================================================================
    central(0)     <- 0
    peripheral1(0) <- 0
    transit1(0)    <- 1
    transit2(0)    <- 1
    transit3(0)    <- 1
    transit4(0)    <- 1
    transit5(0)    <- 1
    tumor(0)       <- r0

    # ===================================================================
    # 8. Observation variables (multi-output)
    # ===================================================================
    # Cc          -- PAC central-compartment concentration in nM (the
    #                PD driver; PK observations from Figure 5A in the
    #                paper's reported nM units after MW conversion).
    # caspase3Act -- relative-to-baseline active caspase-3 (paper's P5
    #                state; baseline level = 1, fold-increase from
    #                treatment).
    # tumorCells  -- JIMT-1 cell count (cells; Figure 6B observation).
    caspase3Act <- transit5
    tumorCells  <- tumor

    Cc          ~ add(addSd)
    caspase3Act ~ add(addSd_caspase3Act)
    tumorCells  ~ add(addSd_tumorCells)
  })
}
