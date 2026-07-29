Li_2012_dexamethasone_rat <- function() {
  description <- "Preclinical (rat). Mechanism-based PK/PD model for CYP3A1/2 induction by dexamethasone (DEX, single 100 mg/kg ip dose) in male Sprague-Dawley rats. PK is a two-compartment mammillary model with zero-order ip absorption of duration T0 directly into the central compartment (no first-order rate constant; CL/F, Q/F, Vc/F, Vp/F all reported as kg-normalised apparent values). PK BSV is exponential and is retained only on Q/F (all other PK BSVs were not significant). The PD cascade describes CYP3A1 and CYP3A2 induction at three molecular levels: (1) mRNA dynamics use an indirect-response (Dayneka-style) model in which a Hill-type fractional occupancy of CYP3A DNA-responsive elements by the DEX-PXR complex (FO = Cp^gamma / (SC50^gamma + Cp^gamma)) drives a stimulation signal Si,0 = Smax * FO that flows through a per-isoform chain of transit compartments with mean transit time tau (one compartment for CYP3A1, eight compartments for CYP3A2) before stimulating mRNA synthesis as d/dt(mRNAi) = kin,i * (1 + Si,ni) - kout,i * mRNAi. (2) Protein dynamics translate mRNA to CYP3A protein via d/dt(CYP3Ai) = ksyn,i * mRNAi^mi - kdeg,i * CYP3Ai, where the per-isoform amplification factor mi is a paper-mechanistic power exponent on mRNA. (3) Enzyme activity (rate of 6beta-hydroxytestosterone formation) is the algebraic linear combination EA = alpha * CYP3A1 + beta * CYP3A2 with per-isoform turnover-number rates alpha and beta (pmol 6beta-OHT / min / pmol CYP3A). The PK and PD layers were fit sequentially in NONMEM 7.1.2 with FOCE+I, the PK model first then the PD layers simultaneously with PK fixed. Three PD layers were fit by the naive pool approach (each animal contributed one PD observation per time point), so no PD IIVs are present. Numbers of transit compartments (n1 = 1, n2 = 8) are paper-mechanistic fixed structural integers."
  reference <- paste(
    "Li L, Li Z, Deng C, Ning M, Li H, Bi S, Zhou T, Lu W. (2012).",
    "A mechanism-based pharmacokinetic/pharmacodynamic model for CYP3A1/2 induction",
    "by dexamethasone in rats.",
    "Acta Pharmacologica Sinica 33(1):127-136.",
    "doi:10.1038/aps.2011.161"
  )
  vignette <- "Li_2012_dexamethasone_rat"
  units <- list(time = "hour", dosing = "mg/kg", concentration = "ug/mL")

  paper_specific_compartments <- c(
    "stim_cyp3a1_1",
    "stim_cyp3a2_1", "stim_cyp3a2_2", "stim_cyp3a2_3", "stim_cyp3a2_4",
    "stim_cyp3a2_5", "stim_cyp3a2_6", "stim_cyp3a2_7", "stim_cyp3a2_8",
    "mrna_cyp3a1", "mrna_cyp3a2",
    "prot_cyp3a1", "prot_cyp3a2"
  )
  paper_specific_residual_sds <- c(
    "propSd_mrna_cyp3a1", "propSd_mrna_cyp3a2",
    "propSd_prot_cyp3a1", "propSd_prot_cyp3a2",
    "propSd_EA"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley)",
    sex            = "male",
    n_subjects_pk  = 3L,
    n_subjects_pd  = 84L,
    n_studies      = 1L,
    weight_range   = "200-250 g",
    disease_state  = "healthy (normal rats; CYP3A induction by exogenous corticosteroid challenge)",
    dose_range     = "100 mg/kg DEX in 5 mL/kg corn oil, single intraperitoneal dose after 12 h fast",
    regions        = "China (Peking University Health Science Center, Beijing)",
    notes          = paste(
      "PK cohort: serial blood samples from 3 treated rats at 16 time points (0.083, 0.25, 0.5,",
      "0.75, 1, 2, 3, 4, 6, 8, 12, 16, 24, 30, 36, 48 h). PD cohort: 84 rats randomized to DEX or",
      "vehicle (corn oil) and sacrificed at 14 time points (0, 1, 2, 4, 8, 12, 16, 24, 30, 36, 42,",
      "48, 54, 60 h; n = 3 per time point). Liver microsomes prepared by differential",
      "centrifugation; CYP3A1/2 mRNA quantified by RT-PCR with absolute-quantification reference",
      "standards (units: attomol per ug total RNA); CYP3A1/2 protein by non-competitive ELISA",
      "(units: pmol per mg microsomal protein); CYP3A1/2 enzyme activity by testosterone substrate",
      "assay reporting 6beta-OHT formation rate (units: pmol per min per mg microsomal protein).",
      "Acclimatized 10 d at 22 C under 12 h/12 h light/dark before dosing; food and water ad",
      "libitum. Approved by the Peking University Committee on Animal Care and Use; conducted per",
      "European Community guidelines."
    )
  )

  ini({
    # ============================================================
    # PK -- two-compartment with zero-order ip absorption.
    # Apparent (CL/F, Q/F, Vc/F, Vp/F) values from Li 2012 Table 2.
    # ============================================================

    # Reported CL/F = 172.7 mL/kg/h; converted to L/kg/h via /1000 inside model().
    lcl   <- log(172.7)
    label("Apparent systemic clearance CL/F (mL/kg/h)")                            # Li 2012 Table 2: CL/F = 172.7 mL/kg/h (RSE 6.70%)

    lvc   <- log(657.4)
    label("Apparent central volume Vc/F (mL/kg)")                                  # Li 2012 Table 2: Vc/F = 657.4 mL/kg (RSE 12.55%)

    lq    <- log(14.32)
    label("Apparent inter-compartmental clearance Q/F (mL/kg/h)")                  # Li 2012 Table 2: Q/F = 14.32 mL/kg/h (RSE 36.52%)

    lvp   <- log(263.2)
    label("Apparent peripheral volume Vp/F (mL/kg)")                               # Li 2012 Table 2: Vp/F = 263.2 mL/kg (RSE 31.95%)

    ldur  <- log(10.47)
    label("Zero-order absorption duration T0 (h)")                                 # Li 2012 Table 2: T0 = 10.47 h (RSE 6.36%)

    # PK BSV. Li 2012 Methods page (Data analysis): "An exponential error model
    # was selected for modeling the between-subject variability (BSV)." Results
    # (Pharmacokinetics paragraph): "Only the BSV for inter-compartment clearance
    # (Q) was included in the final PK model." Table 2 reports BSV Q/F = 0.655
    # (RSE 36.66%) as the variance on the log scale; encoded as-is.
    etalq ~ 0.655                                                                  # Li 2012 Table 2: BSV Q/F variance = 0.655 (RSE 36.66%)

    # PK residual error. Li 2012 Methods: "The residual variability for both the
    # PK and PD models was modeled initially with a combined error model; if one
    # of the components (additive or proportional) of the residual was
    # negligible, it was deleted from the model." The paper does NOT report a
    # numeric magnitude for the retained residual component(s) of the PK or any
    # PD layer. The propSd value below is therefore a STRUCTURAL PLACEHOLDER
    # (fixed at 0.10 = 10% CV) chosen only to make the nlmixr2 model file
    # syntactically complete; it is NOT an estimate from the source. Forward
    # simulations should use rxode2::zeroRe() or override propSd to recover the
    # paper's typical-value trajectories. See vignette Errata.
    propSd <- fixed(0.10)
    label("Proportional residual SD on DEX plasma concentration -- PLACEHOLDER, not from Li 2012")  # not in source

    # ============================================================
    # PD layer 1: CYP3A1 / CYP3A2 mRNA dynamics (Eq 3-8).
    # ============================================================

    # Production rates kin,i (zero-order, in attomol/h/ug total RNA). The mRNA
    # baseline mRNAi(0) = kin,i / kout,i (Eq 8) -- not estimated separately.
    lkin_mrna_cyp3a1  <- log(7.47)
    label("CYP3A1 mRNA transcription rate kin,1 (attomol/h/ug total RNA)")          # Li 2012 Table 3: kin,1 = 7.47 (RSE 4.65%)
    lkin_mrna_cyp3a2  <- log(40.00)
    label("CYP3A2 mRNA transcription rate kin,2 (attomol/h/ug total RNA)")          # Li 2012 Table 3: kin,2 = 40.00 (RSE 19.56%)

    # First-order mRNA degradation rates kout,i (1/h).
    lkout_mrna_cyp3a1 <- log(0.23)
    label("CYP3A1 mRNA degradation rate kout,1 (1/h)")                              # Li 2012 Table 3: kout,1 = 0.23 (RSE 6.45%)
    lkout_mrna_cyp3a2 <- log(0.197)
    label("CYP3A2 mRNA degradation rate kout,2 (1/h)")                              # Li 2012 Table 3: kout,2 = 0.197 (RSE 21.00%)

    # Maximum stimulation of the transcription rate (unitless fold-change).
    lsmax_mrna_cyp3a1 <- log(22.42)
    label("CYP3A1 mRNA Smax,1 -- max fold-stimulation of kin,1 (unitless)")         # Li 2012 Table 3: Smax,1 = 22.42 (RSE 10.90%)
    lsmax_mrna_cyp3a2 <- log(8.55)
    label("CYP3A2 mRNA Smax,2 -- max fold-stimulation of kin,2 (unitless)")         # Li 2012 Table 3: Smax,2 = 8.55 (RSE 11.72%)

    # DEX plasma concentration at 50% Smax (ug/mL, same units as Cp).
    lsc50_mrna_cyp3a1 <- log(2.39)
    label("CYP3A1 mRNA SC50,1 -- DEX Cp for 50% Smax,1 (ug/mL)")                    # Li 2012 Table 3: SC50,1 = 2.39 ug/mL (RSE 15.89%)
    lsc50_mrna_cyp3a2 <- log(2.82)
    label("CYP3A2 mRNA SC50,2 -- DEX Cp for 50% Smax,2 (ug/mL)")                    # Li 2012 Table 3: SC50,2 = 2.82 ug/mL (RSE 38.50%)

    # Hill exponent on the fractional-occupancy FO = Cp^gamma / (SC50^gamma + Cp^gamma).
    lhill_mrna_cyp3a1 <- log(8.01)
    label("CYP3A1 Hill exponent gamma_1 on DEX-PXR fractional occupancy (unitless)") # Li 2012 Table 3: gamma_1 = 8.01 (RSE 6.92%)
    lhill_mrna_cyp3a2 <- log(5.00)
    label("CYP3A2 Hill exponent gamma_2 on DEX-PXR fractional occupancy (unitless)") # Li 2012 Table 3: gamma_2 = 5.00 (RSE 13.16%)

    # Mean transit time tau_i between each stimulation-chain transit
    # compartment (h). Number of compartments n_i is paper-mechanistic-fixed at
    # n_1 = 1 (CYP3A1) and n_2 = 8 (CYP3A2) and is encoded in the ODE structure
    # rather than as a numeric parameter.
    lmtt_mrna_cyp3a1  <- log(4.59)
    label("CYP3A1 mean transit time tau_1 -- 1 compartment in chain (h)")           # Li 2012 Table 3: tau_1 = 4.59 h (RSE 4.71%); n_1 = 1 FIXED
    lmtt_mrna_cyp3a2  <- log(2.58)
    label("CYP3A2 mean transit time tau_2 -- 8 compartments in chain (h)")          # Li 2012 Table 3: tau_2 = 2.58 h (RSE 12.79%); n_2 = 8 FIXED

    # ============================================================
    # PD layer 2: CYP3A1 / CYP3A2 protein dynamics (Eq 9-10).
    # d/dt(prot_i) = ksyn,i * mRNAi^m_i - kdeg,i * prot_i
    # Initial condition: prot_i(0) = (ksyn,i / kdeg,i) * mRNAi(0)^m_i (Eq 10).
    # ============================================================

    # Protein synthesis rate constants ksyn,i (pmol/h/mg protein per (attomol mRNA/ug RNA)^m_i).
    lksyn_cyp3a1   <- log(0.0359)
    label("CYP3A1 protein synthesis rate ksyn,1 (pmol/h/mg prot / (attomol/ug)^m1)") # Li 2012 Table 4: ksyn,1 = 0.0359 (RSE 18.68%)
    lksyn_cyp3a2   <- log(0.486)
    label("CYP3A2 protein synthesis rate ksyn,2 (pmol/h/mg prot / (attomol/ug)^m2)") # Li 2012 Table 4: ksyn,2 = 0.486 (RSE 11.77%)

    # Protein degradation rate constants kdeg,i (1/h).
    lkdeg_cyp3a1   <- log(0.0268)
    label("CYP3A1 protein degradation rate kdeg,1 (1/h)")                            # Li 2012 Table 4: kdeg,1 = 0.0268 (RSE 15.99%)
    lkdeg_cyp3a2   <- log(0.0567)
    label("CYP3A2 protein degradation rate kdeg,2 (1/h)")                            # Li 2012 Table 4: kdeg,2 = 0.0567 (RSE 11.79%)

    # Paper-mechanistic amplification factor m_i -- a power exponent on mRNA in
    # the protein synthesis term (Eq 9: ksyn,i * mRNAi^m_i). The paper describes
    # m_i as "the amplification factor, indicating that one copy of the mRNA can
    # be translated into multiple copies of the protein", but the reported values
    # m_1 = 0.911 and m_2 = 0.287 are < 1, so the power-law form is sublinear in
    # mRNA. lgamma_<output> follows the parameter-names register precedent for
    # paper-mechanistic power exponents (e.g., Ait-Oudhia 2012 CRP-transit
    # amplification); the lhill canonical is reserved for the Hill exponent in
    # the FO equation above, so the two roles do not collide.
    lgamma_cyp3a1  <- log(0.911)
    label("CYP3A1 protein synthesis power exponent m_1 on mRNA (unitless)")          # Li 2012 Table 4: m_1 = 0.911 (RSE 3.99%)
    lgamma_cyp3a2  <- log(0.287)
    label("CYP3A2 protein synthesis power exponent m_2 on mRNA (unitless)")          # Li 2012 Table 4: m_2 = 0.287 (RSE 9.68%)

    # ============================================================
    # PD layer 3: CYP3A1/2 enzyme activity (Eq 11).
    # EA = alpha * prot_cyp3a1 + beta * prot_cyp3a2
    # alpha, beta are per-isoform catalytic-rate constants (turnover numbers)
    # for 6beta-hydroxytestosterone formation -- canonical kcat role.
    # ============================================================
    lkcat_cyp3a1   <- log(0.463)
    label("CYP3A1 catalytic rate alpha -- 6beta-OHT formation (pmol/min/pmol CYP3A1)") # Li 2012 Table 5: alpha = 0.463 (RSE 22.73%)
    lkcat_cyp3a2   <- log(7.49)
    label("CYP3A2 catalytic rate beta -- 6beta-OHT formation (pmol/min/pmol CYP3A2)")  # Li 2012 Table 5: beta = 7.49 (RSE 5.51%)

    # ============================================================
    # PD residual error -- placeholders. Li 2012 reports neither a residual
    # variance nor a residual SD for any of the four PD outputs (mRNA CYP3A1,
    # mRNA CYP3A2, protein CYP3A1, protein CYP3A2, EA). PD data were fit naive
    # pool (single observation per animal per time point), so individual-level
    # variability is folded into the residual; magnitudes are not in Tables 3,
    # 4, or 5 nor the prose. The four propSd_* fixed placeholders below are
    # ONLY for nlmixr2 syntactic completeness and do not represent estimates
    # from the source. Forward simulations using rxode2::rxSolve()
    # +rxode2::zeroRe() give the paper's typical-value trajectories.
    # See vignette Errata.
    propSd_mrna_cyp3a1 <- fixed(0.10)
    label("CYP3A1 mRNA proportional residual SD -- PLACEHOLDER, not from Li 2012")     # not in source
    propSd_mrna_cyp3a2 <- fixed(0.10)
    label("CYP3A2 mRNA proportional residual SD -- PLACEHOLDER, not from Li 2012")     # not in source
    propSd_prot_cyp3a1 <- fixed(0.10)
    label("CYP3A1 protein proportional residual SD -- PLACEHOLDER, not from Li 2012")  # not in source
    propSd_prot_cyp3a2 <- fixed(0.10)
    label("CYP3A2 protein proportional residual SD -- PLACEHOLDER, not from Li 2012")  # not in source
    propSd_EA          <- fixed(0.10)
    label("Enzyme activity proportional residual SD -- PLACEHOLDER, not from Li 2012") # not in source
  })

  model({
    # ============================================================
    # PK individual parameters and micro-constants.
    # ============================================================
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    q      <- exp(lq + etalq)
    vp     <- exp(lvp)
    t0_dur <- exp(ldur)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment PK with zero-order ip absorption directly into central
    # (Li 2012 Eq 1-2). The dose is delivered to `central` as a zero-order
    # infusion of duration T0 = t0_dur via dur(central); rxode2 expects rate = -2
    # on the dose event to invoke the modeled duration. F is implicitly 1.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(central) <- t0_dur

    # DEX plasma concentration (ug/mL). Dose in mg/kg, volumes in mL/kg, so
    # mass/volume = mg/kg / (mL/kg) = mg/mL; multiplying by 1000 converts to
    # ug/mL to match the assay (LLQ 0.25 ug/mL, Methods).
    Cc <- central / vc * 1000

    # ============================================================
    # PD layer 1: CYP3A1 and CYP3A2 mRNA dynamics (Eq 3-8).
    # ============================================================

    # mRNA-related individual parameters.
    kin_cyp3a1  <- exp(lkin_mrna_cyp3a1)
    kout_cyp3a1 <- exp(lkout_mrna_cyp3a1)
    smax_cyp3a1 <- exp(lsmax_mrna_cyp3a1)
    sc50_cyp3a1 <- exp(lsc50_mrna_cyp3a1)
    hill_cyp3a1 <- exp(lhill_mrna_cyp3a1)
    mtt_cyp3a1  <- exp(lmtt_mrna_cyp3a1)

    kin_cyp3a2  <- exp(lkin_mrna_cyp3a2)
    kout_cyp3a2 <- exp(lkout_mrna_cyp3a2)
    smax_cyp3a2 <- exp(lsmax_mrna_cyp3a2)
    sc50_cyp3a2 <- exp(lsc50_mrna_cyp3a2)
    hill_cyp3a2 <- exp(lhill_mrna_cyp3a2)
    mtt_cyp3a2  <- exp(lmtt_mrna_cyp3a2)

    # Fractional occupancy of the DNA-responsive elements by the DEX-PXR
    # complex (Eq 3). Computed from the DEX plasma concentration Cc with a
    # per-isoform Hill exponent.
    fo_cyp3a1 <- Cc^hill_cyp3a1 / (sc50_cyp3a1^hill_cyp3a1 + Cc^hill_cyp3a1)
    fo_cyp3a2 <- Cc^hill_cyp3a2 / (sc50_cyp3a2^hill_cyp3a2 + Cc^hill_cyp3a2)

    # Stimulation signal entering the first transit (Eq 4).
    s0_cyp3a1 <- smax_cyp3a1 * fo_cyp3a1
    s0_cyp3a2 <- smax_cyp3a2 * fo_cyp3a2

    # Stimulation-signal transit chains (Eq 5-6).
    # CYP3A1: n_1 = 1 (single compartment; output S_1,1 feeds mRNA equation).
    d/dt(stim_cyp3a1_1) <- (s0_cyp3a1 - stim_cyp3a1_1) / mtt_cyp3a1

    # CYP3A2: n_2 = 8 (eight serial compartments; output S_2,8 feeds mRNA equation).
    d/dt(stim_cyp3a2_1) <- (s0_cyp3a2    - stim_cyp3a2_1) / mtt_cyp3a2
    d/dt(stim_cyp3a2_2) <- (stim_cyp3a2_1 - stim_cyp3a2_2) / mtt_cyp3a2
    d/dt(stim_cyp3a2_3) <- (stim_cyp3a2_2 - stim_cyp3a2_3) / mtt_cyp3a2
    d/dt(stim_cyp3a2_4) <- (stim_cyp3a2_3 - stim_cyp3a2_4) / mtt_cyp3a2
    d/dt(stim_cyp3a2_5) <- (stim_cyp3a2_4 - stim_cyp3a2_5) / mtt_cyp3a2
    d/dt(stim_cyp3a2_6) <- (stim_cyp3a2_5 - stim_cyp3a2_6) / mtt_cyp3a2
    d/dt(stim_cyp3a2_7) <- (stim_cyp3a2_6 - stim_cyp3a2_7) / mtt_cyp3a2
    d/dt(stim_cyp3a2_8) <- (stim_cyp3a2_7 - stim_cyp3a2_8) / mtt_cyp3a2

    # mRNA dynamics (Eq 7): drug-stimulated zero-order production minus
    # first-order degradation. Baseline mRNAi(0) = kin,i / kout,i (Eq 8) is
    # set as the compartment initial condition.
    d/dt(mrna_cyp3a1) <- kin_cyp3a1 * (1 + stim_cyp3a1_1) - kout_cyp3a1 * mrna_cyp3a1
    d/dt(mrna_cyp3a2) <- kin_cyp3a2 * (1 + stim_cyp3a2_8) - kout_cyp3a2 * mrna_cyp3a2

    mrna_cyp3a1(0) <- kin_cyp3a1 / kout_cyp3a1
    mrna_cyp3a2(0) <- kin_cyp3a2 / kout_cyp3a2

    # ============================================================
    # PD layer 2: CYP3A1 and CYP3A2 protein dynamics (Eq 9-10).
    # d/dt(prot_i) = ksyn,i * mRNAi^m_i - kdeg,i * prot_i
    # ============================================================
    ksyn_cyp3a1 <- exp(lksyn_cyp3a1)
    kdeg_cyp3a1 <- exp(lkdeg_cyp3a1)
    m_cyp3a1    <- exp(lgamma_cyp3a1)

    ksyn_cyp3a2 <- exp(lksyn_cyp3a2)
    kdeg_cyp3a2 <- exp(lkdeg_cyp3a2)
    m_cyp3a2    <- exp(lgamma_cyp3a2)

    d/dt(prot_cyp3a1) <- ksyn_cyp3a1 * mrna_cyp3a1^m_cyp3a1 - kdeg_cyp3a1 * prot_cyp3a1
    d/dt(prot_cyp3a2) <- ksyn_cyp3a2 * mrna_cyp3a2^m_cyp3a2 - kdeg_cyp3a2 * prot_cyp3a2

    # Protein initial conditions (Eq 10): prot_i(0) = ksyn,i/kdeg,i * mRNAi(0)^m_i.
    prot_cyp3a1(0) <- (ksyn_cyp3a1 / kdeg_cyp3a1) * (kin_cyp3a1 / kout_cyp3a1)^m_cyp3a1
    prot_cyp3a2(0) <- (ksyn_cyp3a2 / kdeg_cyp3a2) * (kin_cyp3a2 / kout_cyp3a2)^m_cyp3a2

    # ============================================================
    # PD layer 3: CYP3A1/2 total enzyme activity (Eq 11).
    # Algebraic linear combination -- alpha and beta are turnover-number
    # rates of 6beta-OHT formation per pmol of each CYP3A isoform.
    # ============================================================
    alpha <- exp(lkcat_cyp3a1)
    beta  <- exp(lkcat_cyp3a2)
    EA <- alpha * prot_cyp3a1 + beta * prot_cyp3a2

    # ============================================================
    # Observation and error models.
    # ============================================================
    Cc          ~ prop(propSd)
    mrna_cyp3a1 ~ prop(propSd_mrna_cyp3a1)
    mrna_cyp3a2 ~ prop(propSd_mrna_cyp3a2)
    prot_cyp3a1 ~ prop(propSd_prot_cyp3a1)
    prot_cyp3a2 ~ prop(propSd_prot_cyp3a2)
    EA          ~ prop(propSd_EA)
  })
}
