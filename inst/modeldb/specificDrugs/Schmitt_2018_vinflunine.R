Schmitt_2018_vinflunine <- function() {
  description <- "Combined population PK / PD model for IV vinflunine in adult cancer patients (Schmitt 2018, 18 phase I/II trials, n=372). Four-compartment IV-infusion popPK with creatinine clearance, body surface area, body weight, and PEGylated liposomal doxorubicin co-administration covariates, plus a five-compartment Friberg-style semi-mechanistic myelosuppression PD model for absolute neutrophil count (proliferation + 3 transit + circulation; linear drug effect 1 - slope*Cc on proliferation; (circ0/circ)^gamma feedback)."
  reference <- paste(
    "Schmitt A, Nguyen L, Zorza G, Ferre P, Petain A (2018).",
    "Better characterization of vinflunine pharmacokinetics variability",
    "and exposure/toxicity relationship to improve its use:",
    "Analyses from 18 trials.",
    "Br J Clin Pharmacol 84(7):1506-1517.",
    "doi:10.1111/bcp.13518.",
    sep = " "
  )
  vignette <- "Schmitt_2018_vinflunine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL", ANC = "10^9/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault, raw mL/min, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cockcroft-Gault creatinine clearance in raw mL/min (no BSA normalization). Reference value 82 mL/min (population median per Schmitt 2018 Table 1). Used as power scaling (CRCL / 82)^0.134 on vinflunine CL.",
      source_name        = "CLCR"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 1.8 m^2 (population median per Schmitt 2018 Table 1). Used as power scaling (BSA / 1.8)^0.542 on vinflunine CL.",
      source_name        = "BSA"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 69 kg (population median per Schmitt 2018 Table 1). Used as power scaling (WT / 69)^0.498 on V3 (peripheral2 volume) and (WT / 69)^0.650 on V4 (peripheral3 volume). The 69 kg reference is paper-consistent (the CL formula in the publication uses Table 1 population medians for CRCL and BSA); see vignette Assumptions and deviations.",
      source_name        = "WT"
    ),
    CONMED_PLDH = list(
      description        = "PEGylated liposomal doxorubicin (PLDH) combination indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant PLDH)",
      notes              = "1 if subject is receiving vinflunine in combination with PEGylated liposomal doxorubicin (PLDH; Doxil / Caelyx); 0 otherwise. Used as power-form factor 0.865^CONMED_PLDH on vinflunine CL (i.e. CL is 86.5% of single-agent value under PLDH co-administration).",
      source_name        = "PLDH"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 372L,
    n_studies      = 18L,
    age_range      = "19-82 years",
    age_median     = "59 years",
    weight_range   = "42-114 kg",
    weight_median  = "69 kg",
    sex_female_pct = 44.4,
    disease_state  = "Adult cancer patients across multiple solid tumour indications (bladder, breast, NSCLC and others) enrolled in 8 phase I, 2 phase I/II, 6 phase II and 2 phase I special-population studies (renal impairment, liver dysfunction). 141 of 372 patients had liver metastases. 43% received vinflunine in combination with one of: cisplatin, gemcitabine, carboplatin, PEGylated liposomal doxorubicin (PLDH), capecitabine, epirubicin, or doxorubicin.",
    dose_range     = "Intravenous vinflunine (Javlor) administered as a short zero-order infusion. SPC doses: 320 mg/m^2 once every 3 weeks (q3w) for normal renal function; 280 mg/m^2 q3w for moderate renal impairment (CrCl 40-60 mL/min); 250 mg/m^2 q3w for severe renal impairment (CrCl 20-<40 mL/min). 4980 vinflunine concentrations across 656 PK profiles after outlier handling (4154 concentrations retained).",
    notes          = "PK analysis: 372 subjects, 4154 vinflunine concentrations (after outlier handling), 656 PK profiles across the 18 studies. PK/PD analysis: 210 vinflunine-monotherapy subjects, 423 administrations, 1871 absolute neutrophil count (ANC) observations from 5 phase I and 6 phase II studies. Cockcroft-Gault CrCl median 82 mL/min (range 29-199). See Schmitt 2018 Table 1 for the full baseline covariate distribution."
  )

  ini({
    # Structural PK parameters - Final covariate model column of Schmitt 2018 Table 2
    # (page 1609 of the published article; matches the explicit CL formula on p.1610).
    lcl  <- log(40.5);  label("Vinflunine clearance CL (L/h)")             # Table 2, final covariate model
    lvc  <- log(20.7);  label("Central volume Vc / V (L)")                 # Table 2, final covariate model
    lvp  <- log(115);   label("Peripheral1 volume Vp / V2 (L)")            # Table 2, final covariate model
    lvp2 <- log(384);   label("Peripheral2 volume Vp2 / V3 (L)")           # Table 2, final covariate model
    lvp3 <- log(593);   label("Peripheral3 volume Vp3 / V4 (L)")           # Table 2, final covariate model

    # Intercompartmental clearances Q1, Q2, Q3 computed from the source paper's
    # back-rate constants K21, K31, K41 via Q = k_xy * V_y (mass-balance):
    #   Q  = K21 * V2 = 0.946  * 115 = 108.79 L/h
    #   Q2 = K31 * V3 = 0.141  * 384 = 54.144 L/h
    #   Q3 = K41 * V4 = 0.0229 * 593 = 13.580 L/h
    # The paper parameterises the 4-compartment model via the back-rate constants
    # rather than Q (NONMEM rate-constant idiom). The Q-based encoding here is
    # mathematically equivalent and follows nlmixr2lib's canonical convention
    # (vc + vp/q + vp2/q2 + vp3/q3); the source-trace anchor remains the
    # K21 / K31 / K41 values in Table 2.
    lq   <- log(108.79);  label("Intercompartmental clearance Q1 = K21*V2 (L/h)")   # Table 2 (K21 = 0.946, V2 = 115)
    lq2  <- log(54.144);  label("Intercompartmental clearance Q2 = K31*V3 (L/h)")   # Table 2 (K31 = 0.141, V3 = 384)
    lq3  <- log(13.580);  label("Intercompartmental clearance Q3 = K41*V4 (L/h)")   # Table 2 (K41 = 0.0229, V4 = 593)

    # Covariate effects on CL - exponents from Schmitt 2018 Table 2 final covariate
    # model and reproduced verbatim in the explicit CL formula on p.1610:
    #   CLi (L/h) = 40.5 * (CLcr / 82)^0.134 * (BSA / 1.8)^0.542 * 0.865^PLDH
    e_crcl_cl <- 0.134;   label("Power exponent of (CRCL / 82) on CL")              # Table 2
    e_bsa_cl  <- 0.542;   label("Power exponent of (BSA / 1.8) on CL")              # Table 2
    e_pldh_cl <- 0.865;   label("Power-form factor on CL under PLDH co-admin (CL *= 0.865^CONMED_PLDH)")  # Table 2

    # Covariate effects on V3 and V4 (peripheral2 / peripheral3) - power form on
    # WT with reference 69 kg (population median per Table 1). The source paper
    # does not print an explicit (WT / ref)^exp formula for V3 / V4 in the way
    # it does for CL; the reference 69 kg is the paper-consistent population
    # median (matching the CL formula's use of CRCL_ref = 82 mL/min and
    # BSA_ref = 1.8 m^2, both Table 1 medians). See vignette Assumptions and
    # deviations.
    e_wt_vp2 <- 0.498;    label("Power exponent of (WT / 69) on V3 / peripheral2")  # Table 2 (WT on V3)
    e_wt_vp3 <- 0.650;    label("Power exponent of (WT / 69) on V4 / peripheral3")  # Table 2 (WT on V4)

    # IIV variances (omega^2 on log-eta scale) computed from CV% via the standard
    # log-normal mapping omega^2 = log(CV^2 + 1):
    #   CL:  25.3% CV -> omega^2 = log(0.253^2 + 1) = 0.06205
    #   V:   73.9% CV -> omega^2 = log(0.739^2 + 1) = 0.4358
    #   V2:  81.3% CV -> omega^2 = log(0.813^2 + 1) = 0.5076
    #   V3:  40.6% CV -> omega^2 = log(0.406^2 + 1) = 0.1526
    #   V4:  39.1% CV -> omega^2 = log(0.391^2 + 1) = 0.1423
    # Correlations from Table 2 are converted to covariances via
    # cov(eta_a, eta_b) = corr * sqrt(omega^2_a * omega^2_b):
    #   corr(CL, V4) = 0.726 -> cov(etalcl, etalvp3) = 0.726 * sqrt(0.06205 * 0.1423) = 0.06823
    #   corr(V, V2)  = 0.644 -> cov(etalvc, etalvp)  = 0.644 * sqrt(0.4358 * 0.5076) = 0.30287
    # IIV on K21, K31, K41 was fixed to 0 in Schmitt 2018 (Table 2 footnote);
    # we do not declare etas on the Q parameters because they inherit the
    # underlying rate-constant IIV-fixed-at-zero assumption.
    etalcl + etalvp3 ~ c(0.06205,
                         0.06823, 0.1423)                                 # 2x2 block (Table 2 corr CL-V4 = 0.726)
    etalvc + etalvp  ~ c(0.4358,
                         0.30287, 0.5076)                                 # 2x2 block (Table 2 corr V-V2 = 0.644)
    etalvp2          ~ 0.1526                                             # Table 2 (IIV V3 = 40.6%)

    # Residual error on plasma vinflunine concentration: proportional, 20.3% CV
    # per Schmitt 2018 Table 2 final covariate model ("Residual variability (%)").
    propSd <- 0.203;    label("Proportional residual error on vinflunine concentration (fraction)")  # Table 2

    # ----- Friberg-style myelosuppression PD layer (Table 3) -----
    # Final PK/PD model parameters from Schmitt 2018 Table 3 (page 1611).
    # Five-compartment structure (1 proliferation + 3 transit + 1 circulating).
    lcirc0 <- log(4.53);  label("Baseline circulating ANC Circ0 (10^9 cells/L)")       # Table 3 (Base / Circ0)
    lmtt   <- log(124);   label("Mean transit time MTT through proliferation -> circ chain (h)")  # Table 3 (MTT)
    lslope <- log(0.425); label("Linear drug-effect slope on vinflunine concentration ((ng/mL)^-1)")  # Table 3 (Slope)
    gamma  <- fixed(0.170); label("Feedback exponent gamma on (Circ0 / circ) (unitless)")  # Table 3 (no IIV in source)

    # PD IIV - omega^2 computed from CV% as above:
    #   Circ0: 46.5% CV -> omega^2 = log(0.465^2 + 1) = 0.1957
    #   Slope: 35.9% CV -> omega^2 = log(0.359^2 + 1) = 0.1213
    # MTT IIV is fixed at 0 in Schmitt 2018 (Table 3); IOV on MTT (24%) and
    # Slope (23.8%) is not represented here (nlmixr2 does not natively model
    # interoccasion variability without an OCC data column; see vignette
    # Assumptions and deviations).
    etalcirc0 ~ 0.1957                                                    # Table 3 (IIV Circ0 = 46.5%)
    etalslope ~ 0.1213                                                    # Table 3 (IIV Slope = 35.9%)

    # PD residual error: Schmitt 2018 Table 3 reports "additive model on
    # log-transformed ANC" with CV% = 42.4%. This is the NONMEM
    # Y = IPRED * EXP(EPS(1)) form, i.e. log-normal residual on ANC.
    # We encode it as expSd = 0.424 (log-scale SD), consistent with the
    # standard CV%-as-log-SD reporting convention for log-additive errors.
    expSd_ANC <- 0.424;  label("Log-normal residual SD on ANC (log-scale)")  # Table 3 (Residual log-additive CV%)
  })

  model({
    # ----- Individual PK parameters with covariate effects -----
    # CLi (L/h) = 40.5 * (CRCL/82)^0.134 * (BSA/1.8)^0.542 * 0.865^CONMED_PLDH
    cl  <- exp(lcl  + etalcl)  * (CRCL / 82)^e_crcl_cl *
                                  (BSA  / 1.8)^e_bsa_cl *
                                  e_pldh_cl^CONMED_PLDH
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)
    # V3 and V4 (peripheral2 / peripheral3) carry the WT effect with reference 69 kg
    vp2 <- exp(lvp2 + etalvp2) * (WT / 69)^e_wt_vp2
    vp3 <- exp(lvp3 + etalvp3) * (WT / 69)^e_wt_vp3
    q   <- exp(lq)
    q2  <- exp(lq2)
    q3  <- exp(lq3)

    # Micro-constants for the linear 4-compartment ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    k14 <- q3 / vc
    k41 <- q3 / vp3

    # ----- Linear 4-compartment IV PK -----
    # IV-infusion doses enter `central` directly (no depot). The data set
    # supplies dose duration / rate via RATE or DUR columns; the model code
    # does not impose a structural infusion duration.
    d/dt(central)     <- -kel * central - k12 * central - k13 * central -
                          k14 * central +
                          k21 * peripheral1 +
                          k31 * peripheral2 +
                          k41 * peripheral3
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(peripheral3) <-  k14 * central - k41 * peripheral3

    # Plasma vinflunine concentration in ng/mL.
    # central is carried in mg (matching the mg dosing convention);
    # Cc [ng/mL] = (central [mg] / vc [L]) * 1000.
    Cc <- 1000 * central / vc

    # ----- Friberg-style myelosuppression PD -----
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt)
    slope <- exp(lslope + etalslope)

    # Transit-rate constant ktr = (n_transit + 1) / MTT for a chain of
    # 3 transit compartments plus a proliferation compartment, per the
    # Friberg 2002 / Schmitt 2018 5-compartment myelosuppression structure
    # ("the model consists of a series of compartments (i.e. five)").
    ktr <- 4 / mtt

    # Drug effect on proliferation rate and feedback from circulating cells
    edrug <- 1 - slope * Cc
    feed  <- (circ0 / circ)^gamma

    # Five-compartment Friberg myelosuppression chain:
    # precursor1 = proliferation; precursor2..precursor4 = transit; circ = circulating ANC.
    d/dt(circ)       <- ktr * precursor4 - ktr * circ
    d/dt(precursor1) <- ktr * precursor1 * edrug * feed - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4

    # Initial conditions: all five myelosuppression compartments start at
    # the individual baseline circulating ANC.
    circ(0)       <- circ0
    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0

    # Observation: absolute neutrophil count (10^9 cells/L)
    ANC <- circ

    # Two observation variables; per-output residual error.
    Cc  ~ prop(propSd)
    ANC ~ lnorm(expSd_ANC)
  })
}
