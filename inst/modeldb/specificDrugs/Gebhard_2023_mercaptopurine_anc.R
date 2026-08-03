Gebhard_2023_mercaptopurine_anc <- function() {
  description <- "Joint PK/PD model of oral 6-mercaptopurine (6MP) myelosuppression during maintenance therapy for childhood acute lymphoblastic leukemia (Gebhard 2023, final model PKPD^6MP_{lin,mm}). This is the paper's headline model. A FIXED one-compartment plasma layer with first-order absorption feeds a red-cell compartment holding 6-thioguanine nucleotides (E-TGN) via SATURABLE Michaelis-Menten influx; the red-cell concentration then drives a linear effect function Edrug = slope * E-TGN that suppresses the proliferation rate of a Friberg-style myelosuppression chain (proliferating pool, three transit compartments, circulating neutrophils) with power-law rebound feedback (base / ANC)^gamma. Because every patient's record begins mid-maintenance-therapy rather than at drug-free baseline, the red-cell compartment is initialised at the patient's first observed E-TGN (covariate BL_TGN_RBC) and the Friberg chain is initialised at a treatment steady state scaled by an estimated fraction-of-baseline parameter inieff. Estimated on 452 children with 4624 E-TGN and 9808 ANC observations (Table 5, whole-data-set column). Two outputs: E-TGN (umol/L) and ANC (G/L). Volumes are per m^2 body-surface area and doses are supplied in umol/m^2, so BSA cancels and no BSA covariate is needed."
  reference <- paste(
    "Gebhard A, Lilienthal P, Metzler M, Rauh M, Sager S, Schmiegelow K, Toksvang LN, Zierk J. (2023).",
    "Pharmacokinetic-pharmacodynamic modeling of maintenance therapy for childhood acute lymphoblastic leukemia.",
    "Sci Rep 13:11749. doi:10.1038/s41598-023-38414-0.",
    "Fixed plasma-PK parameters are reproduced in Gebhard 2023 Supplementary Table S11 from",
    "Lennard L, Keen D, Lilleyman JS (1986) and Brunton LL, Lazo JS, Parker K (2005);",
    "k_circ is fixed from Jost F, Zierk J, Le TT, et al. (2020).",
    "The myelosuppression structure is Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO. (2002).",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 20(24):4713-4721. doi:10.1200/JCO.2002.02.140.",
    "The 6MP-only PK arm is modellib('Gebhard_2023_mercaptopurine') and the methotrexate arm is",
    "modellib('Gebhard_2023_methotrexate').",
    sep = " "
  )
  vignette <- "Gebhard_2023_leukemia_maintenance_therapy"
  units <- list(
    time          = "day",
    dosing        = "umol/m^2",
    concentration = "umol/L",
    anc           = "G/L"
  )
  # Unit note. Gebhard 2023 states that the E-TGN OBSERVATIONS were converted to
  # umol/L (assuming a hemoglobin molecular weight of 64458 g/mol and 330 g Hb/L
  # erythrocytes) but never states the units of the PLASMA concentrations. They
  # are fixed by dimensional analysis of equation 6MP3b: V_mm is reported in
  # umol/L/day and K_mm in umol/L (Table 5), and K_mm is added directly to the
  # plasma concentration X_C / V_C in the denominator, so the plasma
  # concentration must also be in umol/L. Consistently, slope^6MP is reported in
  # L/umol so that Edrug = slope * E-TGN is dimensionless. Volumes are L/m^2
  # (Supplementary Table S11), so doses must be supplied in umol/m^2 and BSA
  # cancels out of the model entirely. Convert a clinical mg/m^2 dose with the
  # 6-mercaptopurine molecular weight 152.18 g/mol:
  #   amt [umol/m^2] = dose [mg/m^2] * 1000 / 152.18.
  # That conversion factor is a derived quantity, not printed in the paper.

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot      = list(analyte = "6-mercaptopurine (6MP)", units = NA_character_, specimen = "administration site", verified = FALSE),
    central    = list(analyte = "6-mercaptopurine (6MP)", units = NA_character_, specimen = "plasma", verified = FALSE),
    rbc_tgn    = list(analyte = "6-thioguanine nucleotides (E-TGN)", units = NA_character_, specimen = "blood cell", verified = FALSE),
    precursor1 = list(analyte = "neutrophils", units = NA_character_, specimen = "not applicable", verified = FALSE),
    precursor2 = list(analyte = "neutrophils", units = NA_character_, specimen = "not applicable", verified = FALSE),
    precursor3 = list(analyte = "neutrophils", units = NA_character_, specimen = "not applicable", verified = FALSE),
    precursor4 = list(analyte = "neutrophils", units = NA_character_, specimen = "not applicable", verified = FALSE),
    circ       = list(analyte = "neutrophils", units = NA_character_, specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    BL_TGN_RBC = list(
      description        = "Baseline (first observed) erythrocyte 6-thioguanine-nucleotide concentration, used as the initial condition of the rbc_tgn compartment.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Gebhard 2023 Supplementary Information section I: X_E^6MP(0) = INITGN. Every patient's record begins after weeks of maintenance therapy, 28 days after the last high-dose intravenous MTX course, so the red-cell pool has already accumulated and cannot be initialised at zero. Cohort median 0.83 umol/L, range 0-7.6 umol/L (Table 1). Originally measured in nmol/mmol hemoglobin and converted to umol/L assuming Hb MW 64458 g/mol and 330 g Hb/L erythrocytes (Gebhard 2023 Methods). Note the model's mid-therapy Friberg initialisation uses the ESTIMATED inieff parameter rather than an observed baseline ANC, so no baseline-neutrophil covariate is required.",
      source_name        = "INITGN"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 452L,
    n_studies      = 1L,
    age_range      = "2.4-16.9 years (median 5.9)",
    weight_range   = "10.3-105.5 kg (median 21.5)",
    height_range   = "81.5-180.0 cm (median 114.0)",
    disease_state  = "Precursor-B-cell acute lymphoblastic leukemia in children, during oral maintenance therapy with 6-mercaptopurine plus low-dose oral methotrexate. Treated on the NOPHO ALL-92 protocol (Nordic Society for Paediatric Haematology and Oncology); Copenhagen ethics approval V.200.2080/91.",
    dose_range     = "6-Mercaptopurine 5.4-175.0 mg/m^2 daily by mouth (median 57.1 mg/m^2); the concurrent methotrexate dose was 1.3-45.0 mg/m^2 weekly (median 15.0 mg/m^2). Doses were titrated to a target white-blood-cell count of 1.5-3.5 G/L; the model's own ANC target range is 0.5-2.0 G/L, following Jost et al.",
    regions        = "Nordic countries (NOPHO ALL-92)",
    observations   = "9808 ANC observations and 4624 E-TGN observations across 452 patients; the data set also holds 4192 E-MTX observations, fitted by the decoupled MTX arm (Table 1).",
    baseline_anc   = "Cohort median ANC 1.6 G/L, range 0-22.5 G/L (Table 1). The model's estimated drug-free baseline is base = 2.17 G/L (Table 5).",
    notes          = "Each patient's record begins with the first paired E-TGN / E-MTX measurement taken 28 days after the last high-dose intravenous methotrexate course, so the high-dose-MTX period is excluded and the record starts mid-therapy. Patients were dropped when parameter estimation would have been impossible: fewer than two E-TGN or E-MTX observations, no 6MP or MTX dose, or no height, weight or ANC observation. Estimation used NONMEM 7.5 with IMP followed by SAEM for this full PKPD model. Population-prediction ANC errors for this model (Table 4): median RMSE 1.12 G/L, mean RMSE 1.27 G/L, median MAE 0.90 G/L -- roughly half those of the PKPD_Jost reference model (median RMSE 2.12 G/L)."
  )

  ini({
    # =========================================================================
    # Plasma layer -- FIXED one-compartment model, Gebhard 2023 Supplementary
    # Table S11 (identical to Suppl Table S2 of the 6MP-only PK model). The
    # paper reports micro-constants (ka, ke, V_C); they are reparameterised
    # here to the canonical lka / lcl / lvc. The product is written out inline
    # so the arithmetic is auditable.
    # =========================================================================
    lka     <- fixed(log(21.07));          label("ka: first-order absorption rate constant (1/day)")                 # Gebhard 2023 Suppl Table S11: ka^6MP = 21.07 1/day
    lvc     <- fixed(log(20.1911));        label("V_C: central volume per body-surface area (L/m^2)")                # Gebhard 2023 Suppl Table S11: V_C^6MP = 20.1911 L/m^2
    lcl     <- fixed(log(15.4 * 20.1911)); label("CL: total clearance per body-surface area (L/day/m^2) = ke * V_C") # Gebhard 2023 Suppl Table S11: ke^6MP = 15.4 1/day, V_C^6MP = 20.1911 L/m^2
    lfdepot <- fixed(log(0.12));           label("F: oral bioavailability of 6-mercaptopurine (fraction)")           # Gebhard 2023 Suppl Table S11: F^6MP = 0.12

    # Circulating-neutrophil elimination rate, FIXED following Jost et al.
    # (Gebhard 2023 main text: "Similar to Jost et al. we fixed both the
    # proliferating rate of the first compartment and k_circ, with
    # k_prol = k_tr and k_circ = 2.3765 1/day"). Cross-check: the implied
    # mature-neutrophil half-life ln(2) / 2.3765 = 0.2917 day is exactly the
    # t_0.5 = 0.29 quoted in the Supplementary sensitivity analysis (Table S12).
    lkout <- fixed(log(2.3765)); label("k_circ: circulating-neutrophil elimination rate constant (1/day)") # Gebhard 2023 Suppl Table S11: k_circ = 2.3765 1/day (Jost et al. 2020)

    # =========================================================================
    # ESTIMATED on the whole data set, Gebhard 2023 Table 5 (final model
    # PKPD^6MP_{lin,mm}). Values come from the whole-data-set column, not the
    # cross-validation column and not initial estimates. Note the values differ
    # from the standalone 6MP PK model (Table 3) -- the paper reports that
    # V_mm^6MP and K_mm^6MP "increase distinctly" when the ANC data are fitted
    # jointly -- so these are NOT interchangeable with Gebhard_2023_mercaptopurine.
    # =========================================================================
    # -- 6MP red-cell layer --
    lvmax_rbc <- log(0.21);  label("V_mm: maximum saturable influx rate, plasma -> red cells (umol/L/day)")          # Gebhard 2023 Table 5: V_mm^6MP  = 0.21 umol/L/day (RSE 4%)
    lkm_rbc   <- log(0.14);  label("K_mm: plasma 6MP concentration at half-maximal red-cell influx (umol/L)")        # Gebhard 2023 Table 5: K_mm^6MP  = 0.14 umol/L (RSE 0.1%)
    lkeff_rbc <- log(0.050); label("K_eff: first-order efflux / elimination rate constant out of red cells (1/day)") # Gebhard 2023 Table 5: K_eff^6MP = 0.050 1/day (RSE 2%)

    # -- Friberg myelosuppression layer --
    lrbase  <- log(2.17); label("base: drug-free baseline absolute neutrophil count (G/L)")                       # Gebhard 2023 Table 5: base = 2.17 G/L (RSE 3%)
    lktr    <- log(0.15); label("k_tr: transit rate constant through the maturation chain (1/day)")               # Gebhard 2023 Table 5: k_tr = 0.15 1/day (RSE 2%); k_prol is fixed equal to k_tr per the main text
    lslope  <- log(0.16); label("slope: linear drug-effect coefficient on red-cell TGN (L/umol)")                 # Gebhard 2023 Table 5: slope^6MP = 0.16 L/umol (RSE 4%)
    lgamma  <- log(0.79); label("gamma: exponent of the (base / ANC) rebound-feedback term (dimensionless)")      # Gebhard 2023 Table 5: gamma = 0.79 (RSE 5%)
    linieff <- log(0.87); label("inieff: fraction of baseline at which the chain is initialised (dimensionless)") # Gebhard 2023 Table 5: inieff = 0.87 (RSE 22%)

    # IIV. Table 5 reports the coefficient of variation, and its footnote gives
    # the definition CV = sqrt(exp(omega^2) - 1); inverting, omega^2 =
    # log(CV^2 + 1). The inversion is written out inline so it is auditable.
    # There is deliberately NO eta on K_mm: Table 5 leaves the K_mm^6MP row of
    # the CV block empty (printed as "-"), matching the 'pop' parameterisation
    # carried over from the standalone 6MP PK model.
    etalvmax_rbc ~ log(0.56^2 + 1)  # Gebhard 2023 Table 5: CV(V_mm^6MP)  = 56% [46, 65], eta-shrinkage 27% -> omega^2 = 0.27273
    etalkeff_rbc ~ log(0.54^2 + 1)  # Gebhard 2023 Table 5: CV(K_eff^6MP) = 54% [45, 63], eta-shrinkage 29% -> omega^2 = 0.25581
    etalrbase    ~ log(0.27^2 + 1)  # Gebhard 2023 Table 5: CV(base)      = 27% [24, 31], eta-shrinkage 22% -> omega^2 = 0.07037
    etalktr      ~ log(0.72^2 + 1)  # Gebhard 2023 Table 5: CV(k_tr)      = 72% [60, 82], eta-shrinkage 35% -> omega^2 = 0.41769
    etalslope    ~ log(0.81^2 + 1)  # Gebhard 2023 Table 5: CV(slope^6MP) = 81% [65, 97], eta-shrinkage 40% -> omega^2 = 0.50447
    etalgamma    ~ log(0.11^2 + 1)  # Gebhard 2023 Table 5: CV(gamma)     = 11% [8, 12],  eta-shrinkage 49% -> omega^2 = 0.01203
    etalinieff   ~ log(0.54^2 + 1)  # Gebhard 2023 Table 5: CV(inieff)    = 54% [48, 59], eta-shrinkage 20% -> omega^2 = 0.25581

    # -------------------------------------------------------------------------
    # Residual unexplained variability. Gebhard 2023 does not state whether the
    # tabulated numbers are variances or standard deviations; they are
    # VARIANCES. The 6MP arm decides it: read as variances the combined
    # residual SD at the cohort median E-TGN of 0.83 umol/L is
    # sqrt(0.024 + (0.24698 * 0.83)^2) = 0.257 umol/L, matching the observed
    # median / mean individual-prediction RMSE of 0.21 / 0.24 reported for the
    # same arm in Table 3; read as standard deviations it would be 0.058,
    # about 3.6 times too small. The ANC reading is corroborated
    # independently: as a variance, 0.25 gives a 50% proportional residual CV,
    # so at the estimated baseline of 2.17 G/L the residual SD is 1.09 G/L
    # against Table 4's population-prediction median / mean ANC RMSE of
    # 1.12 / 1.27 G/L -- a direct match; as a standard deviation it would be
    # 25%, giving 0.54 G/L, half the observed error. The square roots are
    # written out inline so the reading is auditable. Table 5 also reports the
    # decoupled MTX arm's residual errors (additive 0.000018, proportional
    # 0.024); those belong to Gebhard_2023_methotrexate.R -- see the note on
    # the MTX arm in the model() body below.
    # -------------------------------------------------------------------------
    addSd_Crbc_tgn  <- sqrt(0.024); label("Additive residual SD on E-TGN (umol/L)")        # Gebhard 2023 Table 5: additive residual error, 6MP = 0.024 (RSE 9%), read as a variance
    propSd_Crbc_tgn <- sqrt(0.061); label("Proportional residual SD on E-TGN (fraction)")  # Gebhard 2023 Table 5: proportional residual error, 6MP = 0.061 (RSE 5%), read as a variance
    propSd_ANC  <- sqrt(0.25);  label("Proportional residual SD on ANC (fraction)")    # Gebhard 2023 Table 5: proportional residual error, ANC = 0.25 (RSE 2%), read as a variance
  })

  model({
    # ---- 1. Individual parameters ----
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    vmax_rbc <- exp(lvmax_rbc + etalvmax_rbc)
    km_rbc   <- exp(lkm_rbc)                    # population-only: no eta (Table 5 CV column is "-")
    keff_rbc <- exp(lkeff_rbc + etalkeff_rbc)

    rbase  <- exp(lrbase  + etalrbase)
    ktr    <- exp(lktr    + etalktr)
    slope  <- exp(lslope  + etalslope)
    gamma  <- exp(lgamma  + etalgamma)
    inieff <- exp(linieff + etalinieff)
    kout   <- exp(lkout)                        # k_circ, fixed at 2.3765 1/day

    # ---- 2. Micro-constants ----
    kel <- cl / vc   # = ke^6MP = 15.4 1/day

    # ---- 3. 6MP PK arm: equations 6MP1, 6MP2 and 6MP3b ----
    # Reproduced verbatim from Gebhard 2023 Supplementary Information section I
    # ("The final PKPD model consists of the ODE system below"). Cc must be
    # defined before the red-cell ODE that consumes it.
    d/dt(depot)   <- -ka * depot                 # 6MP1
    d/dt(central) <-  ka * depot - kel * central # 6MP2
    Cc <- central / vc
    d/dt(rbc_tgn) <- vmax_rbc * Cc / (km_rbc + Cc) - keff_rbc * rbc_tgn  # 6MP3b

    # ---- 4. Effect function: equation PKPD7a ----
    # The final model's superscript (PKPD^6MP) denotes that the effect function
    # is driven by E-TGN ALONE. The rejected sibling PKPD^{6MP,MTX} (equation
    # PKPD7b, Suppl Table S10) adds a slope^MTX * E-MTX term. Table 5 does also
    # report K_in^MTX, K_eff^MTX and MTX residual errors because the MTX arm was
    # carried during estimation to fit the E-MTX observations, but it never
    # enters this effect function and so cannot influence ANC -- consistent with
    # the Discussion's "our final PKPD model does not include a MTX PK
    # submodel" and with Suppl section I, which lists no MTX states. That arm
    # is the separate model Gebhard_2023_methotrexate.R, so no endpoint is lost.
    edrug <- slope * rbc_tgn                     # PKPD7a

    # ---- 5. Friberg myelosuppression chain: equations PKPD1-PKPD6 ----
    # Compartment naming follows the Friberg-family precedent in nlmixr2lib
    # (Friberg_2002_paclitaxel, Ozawa_2007_docetaxel, Kawamura_2018_eribulin):
    #   X_prol      -> precursor1 (proliferating pool)
    #   X_tr1..X_tr3 -> precursor2..precursor4 (maturation transit chain)
    #   X_ma        -> circ (mature circulating neutrophils, observed as ANC)
    # k_prol is fixed equal to k_tr by the authors (main text, following Jost
    # et al.), so a single ktr drives both the proliferation and transit terms.
    fb <- (rbase / circ)^gamma                                          # PKPD6
    d/dt(precursor1) <- ktr * precursor1 * (1 - edrug) * fb - ktr * precursor1 # PKPD1
    d/dt(precursor2) <- ktr * (precursor1 - precursor2)                 # PKPD2
    d/dt(precursor3) <- ktr * (precursor2 - precursor3)                 # PKPD3
    d/dt(precursor4) <- ktr * (precursor3 - precursor4)                 # PKPD4
    d/dt(circ)       <- ktr * precursor4 - kout * circ                  # PKPD5

    # ---- 6. Bioavailability ----
    f(depot) <- exp(lfdepot)

    # ---- 7. Initial conditions (Gebhard 2023 Suppl Information section I) ----
    # X_GI(0) = X_C(0) = 0 (rxode2 default); X_E^6MP(0) = INITGN.
    #
    # The PD chain is initialised at a treatment steady state rather than at
    # drug-free baseline, because the data set begins weeks into maintenance
    # therapy: X_ma(0) = inieff * base, and the upstream states at
    # inieff * base * k_ma / k_tr.
    #
    # DERIVATION OF k_ma. The symbol k_ma appears in the published initial
    # conditions but is never defined anywhere in the paper or supplement. It
    # is uniquely determined, not guessed: the authors state "we assume to have
    # reached a treatment steady state", and setting
    #   dX_ma/dt = k_tr * X_tr3 - k_circ * X_ma = 0
    # together with X_ma(0) = inieff * base forces
    #   X_tr3(0) = inieff * base * k_circ / k_tr,
    # i.e. k_ma == k_circ = 2.3765 1/day. That is the value used below (kout).
    rbc_tgn(0)    <- BL_TGN_RBC
    precursor1(0) <- inieff * rbase * kout / ktr
    precursor2(0) <- inieff * rbase * kout / ktr
    precursor3(0) <- inieff * rbase * kout / ktr
    precursor4(0) <- inieff * rbase * kout / ktr
    circ(0)       <- inieff * rbase

    # ---- 8. Observations ----
    # Two endpoints. Cc (defined above) is the plasma 6-mercaptopurine
    # concentration (umol/L); it is not observed in this study but is exposed
    # so downstream users and the validation vignette can inspect the fixed
    # plasma layer that drives the red-cell influx.
    Crbc_tgn <- rbc_tgn   # E-TGN (umol/L)
    ANC  <- circ      # absolute neutrophil count (G/L)
    Crbc_tgn ~ add(addSd_Crbc_tgn) + prop(propSd_Crbc_tgn)
    ANC  ~ prop(propSd_ANC)
  })
}
