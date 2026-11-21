Sou_2026_siremadlin <- function() {
  description <- "Population PK/PD model for oral siremadlin (HDM201, an HDM2-p53 inhibitor) and platelet counts in adults with WT-p53 solid tumours or haematological malignancies (AML, ALL, high-risk MDS). PK: one-compartment disposition with linear clearance and Savic transit-compartment absorption (Monolix Ktr / Mtt / ka parameterisation), body weight on V/F. PD: Friberg-type platelet maturation chain (proliferating pool, three transit compartments, circulating platelets) with a linear siremadlin effect increasing apoptosis of drug-sensitive proliferating precursors, and an irreversible first-order transition of precursors from a drug-sensitive to a non-responsive subpopulation. Baseline platelet count and mean maturation time take separate typical values and separate IIV for haematological-malignancy and solid-tumour patients."
  reference <- paste(
    "Sou T, Malek K, Weber HJ.",
    "Population PK/PD Modeling for Evaluation of Intertwining Effects of",
    "Drug and Disease on Thrombocytopenia in Acute Leukemias.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(3):e70210.",
    "doi:10.1002/psp4.70210.",
    "Model equations from Data S1 (Monolix model code); parameter values",
    "from Table 1.",
    sep = " "
  )
  vignette <- "Sou_2026_siremadlin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    plt = "10^9/L"
  )
  dosing <- "depot"

  # The proliferating precursor pool is split into a drug-sensitive
  # subpopulation (the canonical Friberg `prol` state, which holds every
  # precursor at baseline) and a drug-non-responsive subpopulation that
  # fills only through the irreversible S -> R transition.
  paper_specific_compartments <- c("prol_res")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on V/F, (WT / 70)^0.883. The paper's covariate is Monolix 'tBWKG' (a transformed body weight) and the reference weight is not stated; 70 kg is assumed.",
      source_name = "BWKG"
    ),
    DIS_HEME_MALIG = list(
      description = "1 = haematological malignancy (AML, ALL or high-risk MDS); 0 = solid tumour",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the group-specific typical value AND inter-individual variance of the baseline platelet count (PLTZ) and the platelet mean maturation time (MMTP). The source code carries the complementary indicator SOL (1 = solid tumour), so DIS_HEME_MALIG = 1 - SOL. Both groups share all other parameters.",
      source_name = "SOL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "siremadlin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "siremadlin", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(
      analyte = "drug-sensitive proliferating platelet precursor cells (bone marrow)",
      units = "10^9 cells/L (scaled to the circulating platelet count)",
      specimen = "not applicable",
      verified = TRUE
    ),
    prol_res = list(
      analyte = "drug-non-responsive proliferating platelet precursor cells (bone marrow)",
      units = "10^9 cells/L (scaled to the circulating platelet count)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "maturing platelet precursors (stage 1)",
      units = "10^9 cells/L (scaled to the circulating platelet count)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing platelet precursors (stage 2)",
      units = "10^9 cells/L (scaled to the circulating platelet count)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing platelet precursors (stage 3)",
      units = "10^9 cells/L (scaled to the circulating platelet count)",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "platelets", units = "10^9 cells/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 153L,
    n_studies = 1L,
    age_range = "adults; not reported in the source",
    weight_range = "not reported in the source",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in the source",
    disease_state = "WT-p53 advanced solid tumours (colon, melanoma, sarcoma; 88 subjects, 1638 platelet observations) and haematological malignancies (AML, ALL, high-risk MDS; 65 subjects, 598 platelet observations)",
    dose_range = "Oral siremadlin in four regimens: 12.5-350 mg on day 1 of a 21-day cycle; 120-200 mg on days 1 and 8 of a 28-day cycle; 15-25 mg once daily on days 1-7 of a 28-day cycle; 1-20 mg once daily on days 1-14 of a 28-day cycle",
    regions = "multi-centre; not reported in the source",
    notes = "First-in-human phase I study CHDM201X2101 (NCT02143635). Patients who received any blood-product transfusion were excluded from the analysis. The 153 subjects are those contributing platelet observations; the PK subset is not counted separately in the source. External validation used 8 post-allogeneic-HSCT AML patients in remission from study CHDM201K12201 (NCT05447663), whose observations were not used in fitting. No baseline-demographics table is published."
  )

  ini({
    # PK (Table 1, 'PK' block). Monolix reports V/F in 'L/h', a typo for L.
    lktr <- log(6.69); label("Transit absorption rate constant Ktr (1/h)") # Table 1 'Ktr' = 6.69 1/h
    lmtt <- log(0.793); label("Mean transit time for absorption Mtt (h)") # Table 1 'MTT' = 0.793 h
    lka <- log(3.69); label("First-order absorption rate constant ka (1/h)") # Table 1 'Ka' = 3.69 1/h
    lvc <- log(115); label("Apparent volume of distribution V/F at 70 kg (L)") # Table 1 'V/F' = 115
    lcl <- log(5.9); label("Apparent clearance CL/F (L/h)") # Table 1 'Cl/F' = 5.9 L/h
    e_wt_vc <- 0.883; label("Power exponent of body weight on V/F (unitless)") # Table 1 'beta_V_tBWKG' = 0.883

    # PD (Table 1, 'PD' block). The reference (unsuffixed) typical values
    # are the solid-tumour ones; the _heme values apply when
    # DIS_HEME_MALIG = 1.
    lcirc0 <- log(229); label("Baseline platelet count Circ0 in solid-tumour patients (10^9/L)") # Table 1 'PLTZ_SOL' = 229 G/L
    lcirc0_heme <- log(30.8); label("Baseline platelet count Circ0 in haematological-malignancy patients (10^9/L)") # Table 1 'PLTZ_HEM' = 30.8 G/L
    lmtt_plt <- log(844); label("Platelet mean maturation time MMTP in solid-tumour patients (h)") # Table 1 'MMTP_SOL' = 844 h
    lmtt_plt_heme <- log(142); label("Platelet mean maturation time MMTP in haematological-malignancy patients (h)") # Table 1 'MMTP_HEM' = 142 h
    lgamma <- log(0.113); label("Exponent gamma of the (Circ0/circ) feedback on proliferation (unitless)") # Table 1 'gam' = 0.113
    lslope <- log(0.00184); label("Linear slope SLP of the siremadlin effect on sensitive-precursor apoptosis (mL/ng)") # Table 1 'SLP' = 0.00184
    lksr <- log(1.31e-5); label("First-order transition rate constant kSR of precursors from drug-sensitive to non-responsive (1/h)") # Table 1 'kSR' = 1.31e-05 1/h

    # IIV: Table 1 reports the IIV column 'in SD scale' (Monolix omega);
    # each is squared to a variance. All parameters are log-normal.
    etalktr ~ 0.383161 # Table 1 'Ktr' IIV SD 0.619
    etalmtt ~ 0.103684 # Table 1 'MTT' IIV SD 0.322
    etalka ~ 2.0736 # Table 1 'Ka' IIV SD 1.44
    # Table 1 'corr_V_Cl' = 0.631; cov = 0.631 * 0.558 * 0.315
    etalcl + etalvc ~ c(0.311364, 0.1109109, 0.099225) # Table 1 'Cl/F' IIV SD 0.558, 'V/F' IIV SD 0.315
    etalcirc0 ~ 0.121104 # Table 1 'PLTZ_SOL' IIV SD 0.348
    etalcirc0_heme ~ 0.606841 # Table 1 'PLTZ_HEM' IIV SD 0.779
    etalmtt_plt ~ 0.4356 # Table 1 'MMTP_SOL' IIV SD 0.66
    etalmtt_plt_heme ~ 1 # Table 1 'MMTP_HEM' IIV SD 1
    etalgamma ~ 1.1449 # Table 1 'gam' IIV SD 1.07
    etalslope ~ 0.835396 # Table 1 'SLP' IIV SD 0.914
    etalksr ~ 23.04 # Table 1 'kSR' IIV SD 4.8

    # Residual error (Table 1). Monolix 'combined' PK error, encoded as
    # combined1 (SD = a + b * f).
    addSd <- 0.559; label("Additive residual error, siremadlin (ng/mL)") # Table 1 'a1' = 0.559 ng/mL
    propSd <- 0.337; label("Proportional residual error, siremadlin (fraction)") # Table 1 'b1' = 0.337
    propSd_PLT <- 0.275; label("Proportional residual error, platelet count (fraction)") # Table 1 'b2' = 0.275
  })

  model({
    # PK individual parameters
    ktr <- exp(lktr + etalktr)
    mtt <- exp(lmtt + etalmtt)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl <- exp(lcl + etalcl)
    kel <- cl / vc

    # PD individual parameters. Data S1: 'if SOL==0 Circ0 = PLTZ,
    # MMTPx = MMTP elseif SOL==1 Circ0 = PLTZsol, MMTPx = MMTPsol'. Each
    # group has its own typical value and its own IIV, so both etas are
    # declared and the covariate selects one.
    circ0 <- exp(
      (lcirc0 + etalcirc0) * (1 - DIS_HEME_MALIG) +
        (lcirc0_heme + etalcirc0_heme) * DIS_HEME_MALIG
    )
    mtt_plt <- exp(
      (lmtt_plt + etalmtt_plt) * (1 - DIS_HEME_MALIG) +
        (lmtt_plt_heme + etalmtt_plt_heme) * DIS_HEME_MALIG
    )
    gamma <- exp(lgamma + etalgamma)
    slope <- exp(lslope + etalslope)
    ksr <- exp(lksr + etalksr)

    # Data S1: n = 3; ktrP = (n + 1) / MMTPx; kprol = kcirc = ktrP
    ktr_plt <- 4 / mtt_plt
    kprol <- ktr_plt
    kcirc <- ktr_plt

    # Monolix pkmodel(Ktr, Mtt, ka, V, Cl): Savic transit absorption with a
    # non-integer number of transit compartments nn = Mtt * Ktr - 1 feeding
    # an absorption compartment emptied at ka. Written in the integrated
    # form: depot holds all drug not yet absorbed (transit chain plus
    # absorption compartment); the amount still in the chain is
    # dose * (1 - P(nn + 1, Ktr * tad)), with P the regularised lower
    # incomplete gamma function, and only the remainder is absorbed at ka.
    # This stays bounded when IIV drives nn below 0. The chain term follows
    # the most recent dose only, which is exact for once-daily dosing
    # because Mtt (0.79 h) is far shorter than the dosing interval.
    nn <- mtt * ktr - 1
    tad_depot <- tad(depot)
    chain <- depot
    if (tad_depot > 0) {
      chain <- podo(depot) * (1 - gammap(nn + 1, ktr * tad_depot))
    }
    abs_cmt <- depot - chain

    d/dt(depot) <- -ka * abs_cmt
    d/dt(central) <- ka * abs_cmt - kel * central

    # Dose in mg over V/F in L gives mg/L; x 1000 gives ng/mL, the unit of
    # the additive residual error 'a1'.
    Cc <- 1000 * central / vc

    # Data S1: 'Edrug = SLP * Cc'
    edrug <- slope * Cc

    # Data S1 ODEs (ProlS -> prol, ProlR -> prol_res, T1..T3 ->
    # transit1..3, Circ -> circ)
    d/dt(prol) <- kprol * prol * (circ0 / circ)^gamma -
      ktr_plt * (1 + edrug) * prol - ksr * prol
    d/dt(prol_res) <- kprol * prol_res * (circ0 / circ)^gamma -
      ktr_plt * prol_res + ksr * prol
    d/dt(transit1) <- ktr_plt * (prol + prol_res) - ktr_plt * transit1
    d/dt(transit2) <- ktr_plt * transit1 - ktr_plt * transit2
    d/dt(transit3) <- ktr_plt * transit2 - ktr_plt * transit3
    d/dt(circ) <- ktr_plt * transit3 - kcirc * circ

    # Data S1 initial values: every precursor starts in the sensitive pool
    prol(0) <- circ0
    prol_res(0) <- 0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0

    PLT <- circ

    Cc ~ add(addSd) + prop(propSd) + combined1()
    PLT ~ prop(propSd_PLT)
  })
}
