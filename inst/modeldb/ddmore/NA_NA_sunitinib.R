NA_NA_sunitinib <- function() {
  description <- "Semi-mechanistic PK/PD/tumor-growth model for sunitinib in non-small cell lung cancer (NSCLC) patients (DDMORE Foundation Model Repository entry DDMODEL00000231; MDL/PharmML deposit; no linked publication identified). Parent and metabolite each follow a 2-compartment oral PK model with first-order absorption from a separate depot; both depots receive the dose, with effective bioavailability split (1 - fp) to the parent and fp to the metabolite (fp = 0.21 hard-coded). Four indirect-response PD biomarker compartments (biom1, biom2, biom3, biom4) are driven by parent (and optionally metabolite) plasma concentration through 1/(1 + pd*conc) inhibition factors; biom1 is on Kout-modulation, biom2-4 are on Kin-modulation. Tumor volume is described by a sphere-volume state (radius is the observation) with a doubling-time-capped exponential growth term modulated by a delayed parent-concentration memory and a lambda-feedback growth-rate state. A resistance-accumulator state, a parent-concentration integrator, and a delayed-signal compartment combine into the tumor's effective growth rate. The model preserves four hard-coded structural placeholders from the MDL: fp = 0.21 (metabolite-formation fraction), th1..4 = 1 (parent-only drive on biomarkers; metabolite drive disabled), thettum = 1 (parent-only drive on tumor; metabolite drive disabled), dres = 0 (no decay of resistance accumulator)."

  reference <- paste(
    "DDMORE Foundation Model Repository: DDMODEL00000231",
    "(MPD6: Sutent / sunitinib semi-mechanistic PK/PD model in non-small",
    "cell lung cancer; MDL/PharmML deposit, version 3 in the dpastoor",
    "scrape). No linked publication has been located on disk; the bundle",
    "ships only Sunitinib_MPD6_model.mdl and the auto-rendered",
    "Sunitinib_MPD6_model.xml. The deposit's metadata description reads:",
    "\"This is a semi-mechanistic PK/PD model for sunitinib therapy in",
    "non-small cell lung cancer patients. It was developed and validated",
    "using clinical trial data provided by the drug developer.\"",
    sep = " "
  )
  vignette <- "NA_NA_sunitinib"
  paper_specific_compartments <- c("biom1", "biom2", "biom3", "biom4", "resistance", "lat_signal", "parent_integ", "lam_feedback")
  paper_specific_residual_sds <- c("propSd_Cc_metab")


  units <- list(time = "day", dosing = "mg", concentration = "mg/mL")
  ddmore_id <- "DDMODEL00000231"
  replicate_of <- NULL

  covariateData <- list(
    # No covariates declared in the MDL parObj or mdlObj; the dataObj only
    # carries ID, TIME, AMT, DVID, DV, MDV. The model is covariate-free.
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "not reported in the DDMORE bundle (no .lst, no Model_Accomodations.text, no linked publication)",
    weight_range   = "not reported in the DDMORE bundle",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Non-small cell lung cancer (NSCLC), per the DDMORE bundle's RDF model-has-description: 'semi-mechanistic PK/PD model for sunitinib therapy in non-small cell lung cancer patients ... developed and validated using clinical trial data provided by the drug developer.'",
    dose_range     = "not reported in the DDMORE bundle. Sunitinib's standard oncology regimen at the time of MPD6 deposit was 50 mg PO QD on a 4-weeks-on / 2-weeks-off schedule (GIST/RCC label) or 37.5 mg PO QD continuously (later label expansion); the MDL is dose-input-agnostic and accepts any AMT.",
    regions        = "not reported in the DDMORE bundle",
    biomarkers     = "Four PD biomarker compartments labelled only as A5, A6, A7, A8 in the MDL. Their biological identity is not specified in the bundle; the standard sunitinib biomarker panel reported elsewhere (e.g. Hansson 2013 GIST / DDMODEL00000197) is VEGF, sVEGFR-2, sVEGFR-3, and sKIT, and the index ordering in this MDL is consistent with Kout-modulated biom1 (VEGF-like) and Kin-modulated biom2/3/4 (soluble-receptor-like), but no source-of-truth in the bundle confirms the mapping. Compartments are therefore named generically.",
    notes          = "DDMORE-source bundle without a linked publication, an Output_real_*.lst, an Output_simulated_*.lst, a Simulated_*.csv, or a Model_Accomodations.text. Subject counts, study counts, age/weight/sex/race distributions, and dose-level breakdowns are all unobtainable from the bundle alone. Population fields can be populated in a follow-up edit if the linked publication is located."
  )

  ini({
    # ----------------------------------------------------------------------
    # All values come from the MDL parObj STRUCTURAL{} and VARIABILITY{}
    # blocks (Sunitinib_MPD6_model.mdl). For DDMORE-source extractions the
    # final estimates normally come from Output_real_*.lst, but this bundle
    # ships no listing -- see vignette Errata. Per the operator decision
    # (sidecar response 030/response-001.json: extract_mdl), the parObj
    # values are treated as the deposited final estimates.
    #
    # Time unit: day (the MDL multiplies all rate constants by 24 inside
    # MODEL_PREDICTION, e.g. Kah = Ka*24, q = 24*Cl/V1).
    # Volume unit: bundle-internal (likely mL given the magnitudes); the
    # absolute concentration scale is not externally checkable.
    # ----------------------------------------------------------------------

    # ----- Parent PK (2-compartment with first-order oral absorption) -----
    lka      <- log(0.0715);  label("Parent absorption rate constant (1/h, multiplied by 24 in model())")  # MDL parObj POP_Ka
    lcl      <- log(26.9);    label("Parent total clearance (mL/h)")                                       # MDL parObj POP_Cl
    lvc      <- log(3220);    label("Parent central volume (mL)")                                          # MDL parObj POP_V1 = 3.22e3
    lq       <- log(17.5);    label("Parent inter-compartmental clearance (mL/h)")                         # MDL parObj POP_QQ
    lvp      <- log(127);     label("Parent peripheral volume (mL)")                                       # MDL parObj POP_V2

    # ----- Metabolite PK (2-compartment with first-order oral absorption) -----
    lka_metab <- log(0.177);  label("Metabolite absorption rate constant (1/h, multiplied by 24 in model())") # MDL parObj POP_Kam
    lcl_metab <- log(15.6);   label("Metabolite total clearance (mL/h)")                                      # MDL parObj POP_Clm
    lvc_metab <- log(3710);   label("Metabolite central volume (mL)")                                         # MDL parObj POP_Vm1 = 3.71e3
    lq_metab  <- log(159);    label("Metabolite inter-compartmental clearance (mL/h)")                        # MDL parObj POP_QQm
    lvp_metab <- log(156);    label("Metabolite peripheral volume (mL)")                                      # MDL parObj POP_Vm2

    # ----- PD biomarker degradation rates (Kout in indirect-response chains) -----
    lkout_biom1 <- log(0.111);   label("Biomarker 1 degradation / turnover rate (1/day)")  # MDL parObj POP_d1
    lkout_biom2 <- log(0.101);   label("Biomarker 2 degradation / turnover rate (1/day)")  # MDL parObj POP_d2
    lkout_biom3 <- log(0.169);   label("Biomarker 3 degradation / turnover rate (1/day)")  # MDL parObj POP_d3
    lkout_biom4 <- log(0.00659); label("Biomarker 4 degradation / turnover rate (1/day)")  # MDL parObj POP_d4

    # ----- Drug-effect (1/(1+pd*conc) inhibition factors) -----
    lpd1 <- log(144);   label("Inverse-EC50 drug-effect coefficient on biomarker 1 (1/(mg/mL))")  # MDL parObj POP_pd1
    lpd2 <- log(22);    label("Inverse-EC50 drug-effect coefficient on biomarker 2 (1/(mg/mL))")  # MDL parObj POP_pd2
    lpd3 <- log(36.4);  label("Inverse-EC50 drug-effect coefficient on biomarker 3 (1/(mg/mL))")  # MDL parObj POP_pd3
    lpd4 <- log(98.1);  label("Inverse-EC50 drug-effect coefficient on biomarker 4 (1/(mg/mL))")  # MDL parObj POP_pd4

    # ----- Biomarker-baseline log-setpoint (POP=0 in MDL; baseline = 1.0) -----
    # The MDL parObj POP_st1..st4 are exactly 0; biomarker steady-state baseline
    # is therefore exp(0) = 1 for every typical-value subject. Per-subject
    # baselines are exp(eta_st_*).
    lst1 <- 0;  label("Log-baseline setpoint for biomarker 1 (typical = 0, baseline = 1)")  # MDL parObj POP_st1
    lst2 <- 0;  label("Log-baseline setpoint for biomarker 2 (typical = 0, baseline = 1)")  # MDL parObj POP_st2
    lst3 <- 0;  label("Log-baseline setpoint for biomarker 3 (typical = 0, baseline = 1)")  # MDL parObj POP_st3
    lst4 <- 0;  label("Log-baseline setpoint for biomarker 4 (typical = 0, baseline = 1)")  # MDL parObj POP_st4

    # ----- Tumor-volume initial-radius log (POP=0; baseline tumor radius = 1) -----
    lx0 <- 0;  label("Log-baseline tumor radius (typical = 0, radius = 1; volume = 4/3*pi*r^3)")  # MDL parObj POP_x0

    # ----- Parent-concentration integrator decay rate -----
    lpdm <- log(0.129);  label("Decay rate of parent-concentration memory integrator (1/day)")  # MDL parObj POP_pdm

    # ----- Lambda-feedback growth-rate state (initial value, growth, offset) -----
    llam     <- log(0.000193);  label("Initial value of lambda-feedback growth-rate state (1/day)")          # MDL parObj POP_lam
    llam0    <- log(0.00189);   label("Lambda offset (scales the resistance feedback into tumor growth, 1/day)") # MDL parObj POP_lam0
    lalphres <- log(0.0102);    label("Lambda-state self-growth rate (1/day)")                                   # MDL parObj POP_alphres

    # ----- Resistance / memory dynamics rate -----
    lpdr <- log(0.0286);  label("Memory accumulation rate for delayed-signal and resistance compartments (1/day)")  # MDL parObj POP_pdr

    # ----------------------------------------------------------------------
    # Inter-individual variability -- MDL parObj VARIABILITY{} block.
    # MDL declares `type is sd` for every per-parameter omega: variance =
    # omega^2. The MDL also declares an OMEGA `type is corr` 14x14 lower-
    # triangular correlation matrix with mostly small (~|0.01|-|0.07|)
    # off-diagonal entries; those off-diagonals are NOT encoded here (see
    # vignette Errata for the simplification). All non-zero diagonal etas
    # are kept as independent IIV.
    #
    # Etas with omega = 0 in the MDL (QQm, pdr, st2, st3, st4, x0, lam0,
    # alphres, pd1, QQ, V2, Vm2, Ka, Kam) are dropped -- those parameters
    # carry no inter-individual variability in the MDL.
    # ----------------------------------------------------------------------
    etalkout_biom1 ~ 24.4;     label("IIV variance on log biomarker-1 Kout (variance = 4.94^2; very large -- see vignette Errata)") # MDL omega_d1 = 4.94 sd
    etalkout_biom2 ~ 0.235;    label("IIV variance on log biomarker-2 Kout (variance = 0.485^2)")  # MDL omega_d2 = 0.485 sd
    etalkout_biom3 ~ 0.267;    label("IIV variance on log biomarker-3 Kout (variance = 0.517^2)")  # MDL omega_d3 = 0.517 sd
    etalkout_biom4 ~ 0.219;    label("IIV variance on log biomarker-4 Kout (variance = 0.468^2)")  # MDL omega_d4 = 0.468 sd
    etalst1        ~ 0.121;    label("IIV variance on log biomarker-1 baseline setpoint (variance = 0.348^2)")  # MDL omega_st1 = 0.348 sd
    etalpdm        ~ 3.03;     label("IIV variance on log parent-integrator decay rate (variance = 1.74^2; very large -- see vignette Errata)")  # MDL omega_pdm = 1.74 sd
    etallam        ~ 4.45;     label("IIV variance on log lambda-feedback initial value (variance = 2.11^2; very large -- see vignette Errata)") # MDL omega_lam = 2.11 sd
    etalpd2        ~ 0.219;    label("IIV variance on log inverse-EC50 for biomarker 2 (variance = 0.468^2)")  # MDL omega_pd2 = 0.468 sd
    etalpd3        ~ 0.976;    label("IIV variance on log inverse-EC50 for biomarker 3 (variance = 0.988^2)")  # MDL omega_pd3 = 0.988 sd
    etalpd4        ~ 0.147;    label("IIV variance on log inverse-EC50 for biomarker 4 (variance = 0.384^2)")  # MDL omega_pd4 = 0.384 sd
    etalcl         ~ 0.0882;   label("IIV variance on log parent CL (variance = 0.297^2)")  # MDL omega_Cl = 0.297 sd
    etalvc         ~ 1.69;     label("IIV variance on log parent Vc (variance = 1.3^2; large -- see vignette Errata)")  # MDL omega_V1 = 1.3 sd
    etalcl_metab   ~ 0.197;    label("IIV variance on log metabolite CL (variance = 0.444^2)")  # MDL omega_Clm = 0.444 sd
    etalvc_metab   ~ 0.824;    label("IIV variance on log metabolite Vc (variance = 0.908^2; large -- see vignette Errata)") # MDL omega_Vm1 = 0.908 sd

    # ----------------------------------------------------------------------
    # Residual error -- MDL parObj b_1..b_6 (proportional) and a_7/b_7
    # (combined additive + proportional on Y7 = tumor radius).
    # MDL OBSERVATION block:
    #   Y1..Y6 = proportionalError(proportional = b_*, eps = EPS_Y, prediction = output*)
    #   Y7     = combinedError1(additive = a_7, proportional = b_7, eps = EPS_Y, prediction = output7)
    # ----------------------------------------------------------------------
    propSd          <- 0.512;  label("Proportional residual error on parent concentration (Y1, fraction)")    # MDL parObj b_1
    propSd_Cc_metab <- 0.429;  label("Proportional residual error on metabolite concentration (Y2, fraction)") # MDL parObj b_2
    propSd_biom1 <- 0.503;  label("Proportional residual error on biomarker 1 (Y3, fraction)")              # MDL parObj b_3
    propSd_biom2 <- 0.137;  label("Proportional residual error on biomarker 2 (Y4, fraction)")              # MDL parObj b_4
    propSd_biom3 <- 0.28;   label("Proportional residual error on biomarker 3 (Y5, fraction)")              # MDL parObj b_5
    propSd_biom4 <- 0.135;  label("Proportional residual error on biomarker 4 (Y6, fraction)")              # MDL parObj b_6
    addSd_tumorRadius  <- 0.241;  label("Additive residual error on tumor radius (Y7)")                     # MDL parObj a_7
    propSd_tumorRadius <- 0.0856; label("Proportional residual error on tumor radius (Y7, fraction)")       # MDL parObj b_7
  })

  model({
    # Hard-coded structural placeholders preserved verbatim from the MDL
    # MODEL_PREDICTION block. The MDL declares these as if they were
    # tunable but assigns each a numeric constant before the DEQ block;
    # they are documented in the vignette Errata.
    fp      <- 0.21  # MDL line `fp=0.21` -- fraction of dose entering the metabolite arm
    th1     <- 1     # MDL `th1=1` -- parent-only drive on biomarker 1 (metabolite drive disabled)
    th2     <- 1     # MDL `th2=1`
    th3     <- 1     # MDL `th3=1`
    th4     <- 1     # MDL `th4=1`
    thettum <- 1     # MDL `thettum=1` -- parent-only drive on tumor
    dres    <- 0     # MDL `dres=0` -- no decay of resistance accumulator
    pi_     <- 3.1416  # MDL hard-codes PI as 3.1416 in the tumor-radius output

    # Individual structural parameters (typical * exp(eta))
    ka        <- exp(lka)         # parent absorption rate (1/h before *24 conversion)
    cl        <- exp(lcl + etalcl)
    vc        <- exp(lvc + etalvc)
    q         <- exp(lq)
    vp        <- exp(lvp)

    ka_metab  <- exp(lka_metab)
    cl_metab  <- exp(lcl_metab + etalcl_metab)
    vc_metab  <- exp(lvc_metab + etalvc_metab)
    q_metab   <- exp(lq_metab)
    vp_metab  <- exp(lvp_metab)

    kout_biom1 <- exp(lkout_biom1 + etalkout_biom1)
    kout_biom2 <- exp(lkout_biom2 + etalkout_biom2)
    kout_biom3 <- exp(lkout_biom3 + etalkout_biom3)
    kout_biom4 <- exp(lkout_biom4 + etalkout_biom4)

    pd1 <- exp(lpd1)
    pd2 <- exp(lpd2 + etalpd2)
    pd3 <- exp(lpd3 + etalpd3)
    pd4 <- exp(lpd4 + etalpd4)

    st1 <- lst1 + etalst1
    st2 <- lst2
    st3 <- lst3
    st4 <- lst4
    x0  <- lx0

    pdm     <- exp(lpdm + etalpdm)
    pdr     <- exp(lpdr)
    lam     <- exp(llam + etallam)
    lam0    <- exp(llam0)
    alphres <- exp(lalphres)

    # MDL conversions (1/h -> 1/day) and micro-constants
    Kah  <- ka       * 24                # parent depot -> central absorption (1/day)
    Kamh <- ka_metab * 24                # metabolite depot -> central absorption (1/day)
    qel       <- 24 * cl       / vc      # parent kel (1/day)
    qel_metab <- 24 * cl_metab / vc_metab  # metabolite kel (1/day)
    k12       <- 24 * q        / vc      # parent central -> peripheral (1/day)
    k21       <- 24 * q        / vp      # parent peripheral -> central (1/day)
    km12      <- 24 * q_metab  / vc_metab  # metabolite central -> peripheral (1/day)
    km21      <- 24 * q_metab  / vp_metab  # metabolite peripheral -> central (1/day)

    # Biomarker baseline setpoints and Kin
    stst1 <- exp(st1)
    stst2 <- exp(st2)
    stst3 <- exp(st3)
    stst4 <- exp(st4)
    kin1 <- stst1 * kout_biom1
    kin2 <- stst2 * kout_biom2
    kin3 <- stst3 * kout_biom3
    kin4 <- stst4 * kout_biom4

    # Tumor doubling-time cap (30-day doubling at maximum growth)
    maxgr <- log(2) / 30

    # Plasma concentrations (MDL output1, output2)
    Cc       <- central       / vc
    Cc_metab <- central_metab / vc_metab

    # 1/(1 + pd*(th*Cp + (1-th)*Cm)) drug-effect inhibition factors per biomarker
    eff1 <- 1 + pd1 * (th1 * Cc + (1 - th1) * Cc_metab)
    eff2 <- 1 + pd2 * (th2 * Cc + (1 - th2) * Cc_metab)
    eff3 <- 1 + pd3 * (th3 * Cc + (1 - th3) * Cc_metab)
    eff4 <- 1 + pd4 * (th4 * Cc + (1 - th4) * Cc_metab)

    # Tumor sphere geometry: A9 (state) = volume, observation Y7 = radius
    tumor_safe <- (tumor + sqrt(tumor * tumor)) / 2  # max(tumor, 0) without branching
    tumorRadius <- ((3 / (4 * pi_)) * tumor_safe)^(1.0 / 3.0)

    # Tumor-growth-rate input. MDL: RATEIN = (min(maxgr, A13) - (lam0/pdr)*A12 + lam0*A11) * A9 if A9>0 else 0.
    # min(a,b) without branching = (a + b - |a - b|) / 2.
    growth_term <- (maxgr + lam_feedback - sqrt((maxgr - lam_feedback)^2)) / 2
    rate_in <- (growth_term - (lam0 / pdr) * parent_integ + lam0 * lat_signal) * tumor_safe

    # Concentration drives positive parts (max(., 0)) for parent and metabolite
    Cc_pos       <- (Cc       + sqrt(Cc       * Cc      )) / 2
    Cc_metab_pos <- (Cc_metab + sqrt(Cc_metab * Cc_metab)) / 2
    resistance_pos <- (resistance + sqrt(resistance * resistance)) / 2

    # ODEs (MDL DEQ block, A1..A15)
    d/dt(depot)             <- -Kah  * depot                                            # MDL A14: parent depot
    d/dt(depot_metab)       <- -Kamh * depot_metab                                      # MDL A15: metabolite depot
    d/dt(central)           <-  Kah  * depot * (1 - fp) - k12  * central       + k21  * peripheral1       - qel       * central        # MDL A1
    d/dt(peripheral1)       <-  k12  * central          - k21  * peripheral1                                                            # MDL A2
    d/dt(central_metab)     <-  Kamh * depot_metab * fp - km12 * central_metab + km21 * peripheral1_metab - qel_metab * central_metab   # MDL A3
    d/dt(peripheral1_metab) <-  km12 * central_metab    - km21 * peripheral1_metab                                                       # MDL A4
    d/dt(biom1) <- kin1                  - kout_biom1 * biom1 / eff1   # MDL A5: drug INCREASES Kout via 1/eff1 inhibition
    d/dt(biom2) <- kin2 / eff2 - kout_biom2 * biom2                    # MDL A6: drug DECREASES Kin
    d/dt(biom3) <- kin3 / eff3 - kout_biom3 * biom3                    # MDL A7
    d/dt(biom4) <- kin4 / eff4 - kout_biom4 * biom4                    # MDL A8
    d/dt(tumor)        <- rate_in                                                                                          # MDL A9
    d/dt(resistance)   <- pdr * (thettum * Cc_pos + (1 - thettum) * Cc_metab_pos) - dres * resistance_pos                  # MDL A10
    d/dt(lat_signal)   <- parent_integ - pdr * lat_signal                                                                  # MDL A11
    d/dt(parent_integ) <- Cc - pdm * parent_integ                                                                          # MDL A12
    d/dt(lam_feedback) <- alphres * lam_feedback                                                                           # MDL A13

    # Initial conditions (MDL DEQ init = ...). Depots A14/A15 receive the
    # dose via standard rxode2 dose events, not via init = D -- that MDL
    # convention is handled by rxode2's dose-event mechanism.
    biom1(0)        <- stst1
    biom2(0)        <- stst2
    biom3(0)        <- stst3
    biom4(0)        <- stst4
    tumor(0)        <- (4.0 / 3.0) * pi_ * exp(x0)^3   # MDL A9 init = (4/3)*3.1416*exp(x0)^3
    lam_feedback(0) <- lam                              # MDL A13 init = lam (parameter)
    # resistance(0), lat_signal(0), parent_integ(0) default to 0 (MDL init = 0).

    # Observations (MDL output1..output7)
    Cc          ~ prop(propSd)                                             # Y1
    Cc_metab    ~ prop(propSd_Cc_metab)                                    # Y2
    biom1       ~ prop(propSd_biom1)                                       # Y3
    biom2       ~ prop(propSd_biom2)                                       # Y4
    biom3       ~ prop(propSd_biom3)                                       # Y5
    biom4       ~ prop(propSd_biom4)                                       # Y6
    tumorRadius ~ add(addSd_tumorRadius) + prop(propSd_tumorRadius)        # Y7
  })
}
