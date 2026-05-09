Guo_2022_PF_06939999 <- function() {
  description <- "Population PK/PD model for PF-06939999 (a small-molecule PRMT5 inhibitor) in 28 adults with advanced solid tumors enrolled in the dose-escalation part of NCT03854227. PK is a two-compartment model with first-order absorption (CL/F, V1/F, Q/F, V2/F, Ka). Plasma SDMA (the PD biomarker for PRMT5 inhibition) is modelled by an indirect-response model with saturable Imax inhibition on zero-order SDMA production (Kin/Kout), the log-transformed SDMA observation taking an additive (log-normal) residual error. Platelet count is described by the Friberg semi-mechanistic myelosuppression model (proliferating cells plus three transit compartments feeding a circulating compartment) with a linear drug effect Slope*Cc on the proliferation rate and feedback (Circ0/circ)^gamma."
  reference <- paste(
    "Guo C, Liao KH, Li M, Wang I-M, Shaik N, Yin D.",
    "PK/PD model-informed dose selection for oncology phase I expansion:",
    "Case study based on PF-06939999, a PRMT5 inhibitor.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12:1619-1625.",
    "doi:10.1002/psp4.12882.",
    sep = " "
  )
  vignette <- "Guo_2022_PF_06939999"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL",
    sdma          = "ng/mL",
    platelet      = "10^9/L"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 28L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adults with advanced or metastatic solid tumors enrolled in the dose-escalation part (Part 1) of the first-in-patient study NCT03854227 of PF-06939999, a PRMT5 inhibitor. Two confirmed partial responses were observed (one each in the 2 mg b.i.d. and 4 mg b.i.d. cohorts). Four subjects experienced dose-limiting toxicities: thrombocytopenia (n = 2) in the 6 mg b.i.d. cohort, anemia (n = 1) in the 8 mg q.d. cohort, and neutropenia (n = 1) in the 6 mg q.d. cohort.",
    dose_range     = "Oral PF-06939999 once daily (q.d.) at 0.5 mg (n = 1), 4 mg (n = 5), 6 mg (n = 6), or 8 mg (n = 3); or twice daily (b.i.d.) at 0.5 mg (n = 1), 1 mg (n = 2), 2 mg (n = 3), 4 mg (n = 3), or 6 mg (n = 4). Recommended dose for expansion (RDE) selected from the simulations is 6 mg q.d.",
    regions        = NA_character_,
    notes          = "Baseline body weight, age, and hepatic function were tested as PK covariates and were not retained in the final model (drop in objective function value < 3.84). Demographic detail (age, weight, sex, race) is summarised in Appendix Table S1, which is not on disk in this worktree; the population fields above therefore record disease state, dosing, and dose-limiting toxicities verbatim from the main text (Analysis Plan and Results). Trial registration: NCT03854227."
  )

  ini({
    # PK structural parameters - Table 1 PK model rows. Apparent (CL/F, V/F, Q/F)
    # because PF-06939999 is dosed orally and absolute bioavailability was not
    # determined.
    lka <- log(2.31);  label("First-order absorption rate Ka (1/h)")          # Table 1, Ka = 2.31 1/h
    lcl <- log(9.53);  label("Apparent clearance CL/F (L/h)")                 # Table 1, CL/F = 9.53 L/h
    lvc <- log(160);   label("Apparent central volume V1/F (L)")              # Table 1, V1/F = 160 L
    lq  <- log(26.2);  label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table 1, Q/F = 26.2 L/h
    lvp <- log(285);   label("Apparent peripheral volume V2/F (L)")           # Table 1, V2/F = 285 L

    # SDMA indirect-response PD parameters - Table 1 SDMA PD model rows.
    # Imax is bounded [0,1] and reported without IIV; carried as a bare
    # parameter to keep the typical-value mapping transparent.
    imax    <- 0.823;        label("Maximum fractional inhibition of SDMA production Imax (unitless, 0-1)")  # Table 1, Imax = 0.823
    lic50   <- log(0.425);   label("PF-06939999 concentration giving 50% of maximum SDMA inhibition IC50 (ng/mL)")  # Table 1, IC50 = 0.425 ng/mL
    lkout   <- log(0.00708); label("First-order SDMA elimination rate Kout (1/h)")     # Table 1, Kout = 0.00708 1/h
    lblsdma <- log(113);     label("Baseline plasma SDMA Kin/Kout (ng/mL)")            # Table 1, Baseline SDMA = 113 ng/mL

    # Friberg semi-mechanistic platelet PD parameters - Table 1 Platelet count PD model rows.
    lmtt   <- log(134);     label("Mean transit time MTT through proliferation -> circulation chain (h)")  # Table 1, MTT = 134 h
    lslope <- log(0.00496); label("Linear drug-effect slope Slope on PF-06939999 plasma concentration (mL/ng)")  # Table 1, Slope = 0.00496 per ng/mL
    lgamma <- log(0.217);   label("Feedback exponent gamma on (Circ0/circ) (unitless)")  # Table 1, gamma = 0.217
    lblplt <- log(232);     label("Baseline circulating platelet count Circ0 (10^9/L)")  # Table 1, Baseline PLT = 232

    # Inter-individual variability - Table 1 IIV (%) column. The paper reports
    # IIV as %CV; the variance entered here uses the exact log-normal
    # relationship omega^2 = log(1 + CV^2). Etas modify the log-scale typical
    # values (eta on lcl, etc.).
    etalcl     ~ 0.141   # CV 38.9% -> log(1 + 0.389^2) = 0.141; Table 1 PK CL/F IIV
    etalvc     ~ 0.317   # CV 61.1% -> log(1 + 0.611^2) = 0.317; Table 1 PK V1/F IIV
    etalblsdma ~ 0.0812  # CV 29.1% -> log(1 + 0.291^2) = 0.0812; Table 1 SDMA Baseline IIV
    etalslope  ~ 0.241   # CV 52.2% -> log(1 + 0.522^2) = 0.241; Table 1 Platelet Slope IIV
    etalgamma  ~ 0.199   # CV 46.9% -> log(1 + 0.469^2) = 0.199; Table 1 Platelet gamma IIV
    etalblplt  ~ 0.0771  # CV 28.3% -> log(1 + 0.283^2) = 0.0771; Table 1 Platelet Baseline IIV

    # Residual error - Table 1 residual error rows. PK and platelet models use
    # an exponential residual error (Methods, Analysis Plan), which maps
    # onto a log-normal residual in nlmixr2 (Cc ~ lnorm(SD)). SDMA was fit on
    # log-transformed observations with additive residual on the log scale, which
    # is mathematically equivalent (log(F) + EPS == log(F * exp(EPS))).
    # The Table 1 estimates are NONMEM SIGMA values (variances on the
    # log / proportional scale); the SDs entered here are sqrt of those values.
    propSd      <- 0.335;  label("PK residual SD (proportional / log-scale)")            # Table 1 PK residual error 0.112 -> sqrt(0.112) = 0.335
    propSd_SDMA <- 0.121;  label("SDMA residual SD (additive on log-transformed SDMA)")  # Table 1 SDMA residual error 0.0146 -> sqrt(0.0146) = 0.121
    propSd_PLT  <- 0.153;  label("Platelet residual SD (proportional / log-scale)")      # Table 1 PLT residual error 0.0235 -> sqrt(0.0235) = 0.153
  })

  model({
    # Individual parameters (typical-value PK and PD parameters with IIV
    # introduced on the log scale per Methods).
    ka      <- exp(lka)
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc + etalvc)
    q       <- exp(lq)
    vp      <- exp(lvp)
    ic50    <- exp(lic50)
    kout    <- exp(lkout)
    blsdma  <- exp(lblsdma + etalblsdma)
    mtt     <- exp(lmtt)
    slope   <- exp(lslope + etalslope)
    gamma   <- exp(lgamma + etalgamma)
    blplt   <- exp(lblplt + etalblplt)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    kin <- blsdma * kout                    # zero-order SDMA production rate at baseline
    ktr <- 4 / mtt                          # transit-rate constant; 4 = (3 transit + 1 proliferation) compartments

    # Plasma PF-06939999 concentration drives both PD models. Doses are
    # administered in mg and volumes are in L, so central / vc has units mg/L
    # (== ug/mL). Multiply by 1000 to express Cc in ng/mL (the unit the paper
    # reports for both PK observations and the SDMA IC50 / platelet Slope
    # parameters; see units$concentration).
    Cc <- (central / vc) * 1000

    # PK ODEs - 2-compartment with first-order oral absorption (Figure 1, PK for PF-06939999)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # SDMA indirect-response (Figure 1, PK/PD for plasma SDMA): saturable Imax
    # inhibition on Kin (zero-order SDMA production), first-order Kout elimination.
    d/dt(sdma) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * sdma

    # Friberg semi-mechanistic platelet model (Figure 1, PK/PD for platelet):
    # proliferating cells (precursor1) feed three transit compartments (precursor2..4)
    # feeding the circulating compartment. Linear drug effect 1 - slope*Cc on
    # proliferation; feedback (Circ0 / circ)^gamma.
    d/dt(precursor1) <-  ktr * precursor1 * (1 - slope * Cc) * (blplt / circ)^gamma - ktr * precursor1
    d/dt(precursor2) <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <-  ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <-  ktr * precursor4 - ktr * circ

    # Initial conditions: SDMA at endogenous baseline; all platelet
    # compartments at the individual baseline platelet count.
    sdma(0)       <- blsdma
    precursor1(0) <- blplt
    precursor2(0) <- blplt
    precursor3(0) <- blplt
    precursor4(0) <- blplt
    circ(0)       <- blplt

    # Multi-output observations and residual error (Methods, Analysis Plan)
    SDMA <- sdma
    PLT  <- circ
    Cc   ~ lnorm(propSd)        # PK exponential residual error
    SDMA ~ lnorm(propSd_SDMA)   # additive residual on log-transformed SDMA == log-normal
    PLT  ~ lnorm(propSd_PLT)    # platelet exponential residual error
  })
}
