Lindauer_2017_pembrolizumab <- function() {
  description <- paste(
    "QSP / mini-PBPK. Translational semi-mechanistic PK/PD/TGI model for the",
    "anti-PD-1 monoclonal antibody pembrolizumab in advanced melanoma. Couples a",
    "two-compartment plasma PK (parallel linear + Michaelis-Menten clearance,",
    "human PK substituted from Elassaiss-Schaap 2017 KEYNOTE-001) to a",
    "Shah-Betts (2012) physiologic tumor tissue compartment (vascular, endosomal,",
    "interstitial sub-spaces with FcRn recycling), mechanistic pembrolizumab-PD-1",
    "binding in both blood and tumor, an indirect-response positive feedback that",
    "upregulates tumor PD-1 expression when the complex forms, and a Simeoni-type",
    "tumor-growth model in which the antitumor effect is a power function of the",
    "tumor receptor occupancy. Mouse-derived parameter estimates plus three human",
    "melanoma growth-rate scenarios (slow/medium/fast) and two kill-rate scaling",
    "options (allometric / growth-proportional) are tabulated in Lindauer 2017",
    "Table 1 and Table S3; the default human parameterisation here is medium",
    "growth with allometric kill-rate scaling (the central reference scenario)."
  )
  reference <- paste(
    "Lindauer A, Valiathan CR, Mehta K, Sriram V, de Greef R, Elassaiss-Schaap J,",
    "de Alwis DP. Translational Pharmacokinetic/Pharmacodynamic Modeling of Tumor",
    "Growth Inhibition Supports Dose-Range Selection of the Anti-PD-1 Antibody",
    "Pembrolizumab. CPT Pharmacometrics Syst Pharmacol. 2017;6(1):11-20.",
    "doi:10.1002/psp4.12130. Human plasma PK adapted from Elassaiss-Schaap J et",
    "al. (2017) CPT Pharmacometrics Syst Pharmacol 6(1):21-28;",
    "see modellib('Elassaiss-Schaap_2017_pembrolizumab'). Tumor-tissue physiologic",
    "structure follows Shah DK, Betts AM. J Pharmacokinet Pharmacodyn.",
    "2012;39:67-86. doi:10.1007/s10928-011-9232-2. Tumor-growth backbone from",
    "Simeoni M et al. Cancer Res. 2004;64(3):1094-1101.",
    sep = " "
  )
  vignette <- "Lindauer_2017_pembrolizumab"
  paper_specific_compartments <- c(
    "tumor_vs", "tumor_is", "tumor_es_ub", "tumor_es_b",
    "complex_blood", "complex_tumor", "target_tumor", "tumor_vol"
  )

  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species         = "human (translational projection from preclinical C57BL/6 mouse MC38 colon-adenocarcinoma allograft)",
    n_subjects      = NA_integer_,
    n_studies       = 0L,
    age_range       = "not applicable (simulated typical adult, 70 kg reference)",
    age_median      = "not applicable",
    weight_range    = "70 kg reference body weight; allometric scaling assumes a single typical adult",
    weight_median   = "70 kg",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "not applicable (translational simulation)",
    disease_state   = "Advanced / metastatic cutaneous melanoma; KEYNOTE-001 expansion-cohort dose-selection rationale",
    dose_range      = "0.1-10 mg/kg IV Q2W or Q3W in the published dose-response simulations; 2 mg/kg Q3W is the recommended lowest maximally efficacious dose",
    regions         = "Translational simulation (no clinical trial population)",
    notes           = paste(
      "The model is fit on mouse PK + receptor-occupancy + tumor-volume data",
      "from MC38-bearing C57BL/6 syngeneic allograft mice receiving the mouse or",
      "rat DX400 surrogate anti-mouse-PD-1 antibody (216 mouse + 100 rat",
      "antibody PK samples from 316 mice; 466 tumor-volume + 139 receptor-occupancy",
      "PD samples; see Lindauer 2017 Tables S1-S2). For human dose-response",
      "simulations the plasma-PK sub-model is replaced by the Elassaiss-Schaap",
      "2017 KEYNOTE-001 PK; the tumor-tissue mini-PBPK and the PD-1 feedback",
      "parameters are kept constant across species, and KdegPD-1 is allometrically",
      "scaled. Three melanoma growth-rate scenarios (fast / medium / slow) and two",
      "kill-rate scaling assumptions (allometric / growth-proportional) are listed",
      "in Table 1 and Table S3; the default coefficients in this file encode the",
      "central reference (medium growth, allometric scaling). Run the validation",
      "vignette for the other five scenarios. FcRn is treated as a conserved",
      "species in the endosomal space (FcRn_free + complex = constant total),",
      "which is the standard Shah & Betts 2012 implementation and resolves",
      "apparent typos in the Lindauer 2017 supplement transcription of the",
      "FcRn dFcRn/dt equation. The bimolecular Kon * (free PD-1) consumption",
      "term in the supplement's central-compartment equation lacks the C1",
      "antibody multiplier; the implementation here applies the canonical mass",
      "balance Kon * C_antibody * (free target). See the validation vignette",
      "Assumptions and deviations section for details."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Pembrolizumab molecular weight (used only inside model() to convert
    # between mg dosing and nmol state). Fixed at the published apparent
    # molecular weight of an IgG4 mAb (~149 kg/mol).
    # ----------------------------------------------------------------------
    mw_mab <- fixed(149000); label("Pembrolizumab molecular weight (g/mol)")                                                  # IgG4 typical MW ~149 kDa (textbook value)

    # ----------------------------------------------------------------------
    # Plasma PK -- HUMAN values from Lindauer 2017 Table 1 "Value in man"
    # column (= Elassaiss-Schaap 2017 KEYNOTE-001 popPK Table 2). All
    # structural plasma-PK parameters are TYPICAL VALUES from a previously
    # fit popPK, so they are wrapped in fixed() here.
    # ----------------------------------------------------------------------
    lcl   <- fixed(log(0.167));   label("Linear clearance CL_lin (L/day)")                                                    # Lindauer 2017 Table 1 row CL human = 167 mL/day = 0.167 L/day
    lvc   <- fixed(log(2.877));   label("Central volume V1 (L)")                                                              # Lindauer 2017 Table 1 row V1 human = 2877 mL = 2.877 L
    lq    <- fixed(log(0.384));   label("Inter-compartmental clearance Q (L/day)")                                            # Lindauer 2017 Table 1 row Q human = 384 mL/day = 0.384 L/day
    lvp   <- fixed(log(2.854));   label("Peripheral volume V2 (L)")                                                           # Lindauer 2017 Table 1 row V2 human = 2854 mL = 2.854 L
    lvmax <- fixed(log(0.114));   label("Maximum non-linear elimination rate Vmax (mg/day)")                                  # Lindauer 2017 Table 1 row Vmax human = 114 ug/day = 0.114 mg/day
    lkm   <- fixed(log(0.078));   label("Michaelis-Menten constant Km (ug/mL)")                                               # Lindauer 2017 Table 1 row Km human = 0.078 ug/mL

    # ----------------------------------------------------------------------
    # Tumor mini-PBPK (Shah & Betts 2012) -- "constant across species" per
    # Lindauer 2017 Table 1 "Value in mouse" column (the human column is
    # marked '-' for these rows, meaning the mouse value carries over).
    # ----------------------------------------------------------------------
    f_v_es     <- fixed(0.005); label("Endosomal-space volume as fraction of total tumor volume (unitless)")                  # Lindauer 2017 Table 1 row V_es = 0.5% of V_tot (Shah & Betts)
    f_v_is     <- fixed(0.55);  label("Interstitial-space volume as fraction of total tumor volume (unitless)")               # Lindauer 2017 Table 1 row V_is = 55% of V_tot (Shah & Betts)
    f_v_vs     <- fixed(0.07);  label("Vascular-space volume as fraction of total tumor volume (unitless)")                   # Lindauer 2017 Table 1 row V_vs = 7% of V_tot (Shah & Betts)
    plq_norm   <- fixed(304.8); label("Tumor plasma flow per unit tissue volume (1/day; = 12.7 1/h * 24)")                    # Lindauer 2017 Table 1 row PLQ = 12.7 L/h/L tissue (Shah & Betts)
    f_lymph    <- fixed(0.002); label("Lymph flow as fraction of plasma flow (unitless)")                                     # Lindauer 2017 Table 1 row L = 0.2% of plasma flow (Shah & Betts)
    clup_norm  <- fixed(0.8784); label("Endosomal pinocytosis per unit endosomal-space volume (1/day; = 0.0366 1/h * 24)")    # Lindauer 2017 Table 1 row CLup = 0.0366 L/h/L endosomal (Shah & Betts)
    kdeg_endo  <- fixed(1029.6); label("Endosomal degradation rate constant of free antibody (1/day; = 42.9 1/h * 24)")       # Lindauer 2017 Table 1 row Kdeg = 42.9 1/h (Shah & Betts)
    v_ref      <- fixed(0.842); label("Vascular reflection coefficient (unitless)")                                           # Lindauer 2017 Table 1 row v_ref = 0.842 (Shah & Betts)
    v_ref_is   <- fixed(0.2);   label("Lymph / interstitial reflection coefficient (unitless)")                               # Lindauer 2017 Table 1 row v_ref_is = 0.2 (Shah & Betts)
    fcrn_init  <- fixed(49800); label("Initial endosomal FcRn concentration (nM; = 49.8 uM)")                                 # Lindauer 2017 Table 1 row FcRni = 49.8 uM (Shah & Betts)
    fr_recycle <- fixed(0.715); label("Fraction of endosomal FcRn-bound antibody recycled to vascular space (unitless)")      # Lindauer 2017 Table 1 row FR = 0.715 (Shah & Betts)

    # ----------------------------------------------------------------------
    # FcRn-antibody binding kinetics -- HUMAN values from Shah & Betts 2012.
    # Both rates fixed. Converted from 1E6/M/h to 1/(nM*day): * 24 / 1e9.
    # ----------------------------------------------------------------------
    kon_fcrn   <- fixed(19.008); label("FcRn-antibody association rate constant (1/(nM*day); = 792 (1E6/M/h) * 24 / 1e9)")    # Lindauer 2017 Table 1 row Kon_FcRn human = 792e6 /M/h (Shah & Betts)
    koff_fcrn  <- fixed(573.6);  label("FcRn-antibody dissociation rate constant (1/day; = 23.9 1/h * 24)")                   # Lindauer 2017 Table 1 row Koff_FcRn human = 23.9 1/h (Shah & Betts)

    # ----------------------------------------------------------------------
    # Pembrolizumab-PD-1 binding (HUMAN in vitro, Merck data on file).
    # Mouse fit used rat DX400 surrogate (340e6 /M/h, 0.106 1/h); replaced
    # here by pembrolizumab.
    # ----------------------------------------------------------------------
    kon_pd1    <- fixed(69.12); label("Pembrolizumab-PD-1 association rate constant (1/(nM*day); = 2880 (1E6/M/h) * 24 / 1e9)")  # Lindauer 2017 Table 1 row Kon_PD-1 human = 2880e6 /M/h (Merck)
    koff_pd1   <- fixed(3.456); label("Pembrolizumab-PD-1 dissociation rate constant (1/day; = 0.144 1/h * 24)")                # Lindauer 2017 Table 1 row Koff_PD-1 human = 0.144 1/h (Merck)

    # ----------------------------------------------------------------------
    # PD-1 receptor expression and feedback -- estimated in mouse, kept
    # constant across species (Lindauer 2017 Table 1 "Constant across
    # species" source rows). N_Tcell uses the human Merck Manuals reference
    # value (792/uL); the mouse model used 1000/uL.
    # ----------------------------------------------------------------------
    n_tcell    <- fixed(792);   label("T-cell concentration in blood (cells per uL of blood)")                                # Lindauer 2017 Table 1 row N_Tcell human = 792 / uL
    n_pd1_tc   <- fixed(10000); label("PD-1 receptors per T cell (receptors/cell)")                                           # Lindauer 2017 Table 1 row N_PD-1_TC = 10000 (assumed)
    tmulti     <- fixed(4.32);  label("Initial ratio of total PD-1 concentration tumor:blood (unitless)")                     # Lindauer 2017 Table S2 mouse PK/PD T_multi = 4.32 (constant across species)
    emax_tp    <- fixed(94.7);  label("Maximal fold-increase of PD-1 production by complex feedback (unitless)")              # Lindauer 2017 Table S2 mouse PK/PD EMAXTP = 94.7 (constant across species)
    ec50_tp    <- fixed(1.46);  label("Pembrolizumab-PD-1 complex concentration at half-maximal feedback (nM)")               # Lindauer 2017 Table S2 mouse PK/PD EC50TP = 1.46 nM (constant across species)

    # ----------------------------------------------------------------------
    # Complex degradation (internalisation, T-cell turnover); mouse estimate
    # 0.0194 1/h allometrically scaled to human: 0.0194 * (70000/20)^-0.25
    # = 0.00246 1/h.
    # ----------------------------------------------------------------------
    kdeg_pd1   <- fixed(0.0590); label("Pembrolizumab-PD-1 complex degradation rate (1/day; = 0.00246 1/h * 24)")             # Lindauer 2017 Table 1 row KdegPD-1 human = 0.00246 1/h (allometric from mouse)

    # ----------------------------------------------------------------------
    # Tumor growth -- HUMAN melanoma reference scenario (medium growth,
    # allometric scaling of the kill rate; Lindauer 2017 Table S3 row
    # "Medium growth / allometric scaling"). Mouse-fit PSI and gamma kept.
    #
    # Other scenarios (Table S3):
    #   slow,   allometric:           l0 = 0.0017, sltg = 2.575E-06
    #   slow,   growth-proportional:  l0 = 0.0017, sltg = 2.886E-07
    #   medium, allometric (default): l0 = 0.0036, sltg = 2.575E-06
    #   medium, growth-proportional:  l0 = 0.0036, sltg = 6.247E-07
    #   fast,   allometric:           l0 = 0.0088, sltg = 2.575E-06
    #   fast,   growth-proportional:  l0 = 0.0088, sltg = 1.542E-06
    # ----------------------------------------------------------------------
    ll0        <- log(0.0036);    label("Tumor exponential growth-rate constant L0 (1/day) -- medium-growth scenario")        # Lindauer 2017 Table S3 medium growth = 0.0036/day (doubling time 193 d)
    ll1        <- fixed(log(1e6)); label("Tumor linear growth-rate constant L1 (mL/day) -- effectively disabled for human melanoma (exponential growth only)")  # Lindauer 2017 Table S3 "L1 NA" (linear growth not used for human; "Methods")
    lrbase     <- log(41.5);      label("Initial tumor volume W0 (mL) at start of treatment")                                 # Lindauer 2017 Table S3 baseline = 41.5 mL (= 64 mm SLD per Chiu/Ouellet, sphere conversion)
    lsltg      <- fixed(log(2.575e-6)); label("Drug-effect slope SLtg on tumor (1/day; allometric scaling)")                  # Lindauer 2017 Table S3 medium/allometric SLtg = 2.575E-06
    lgamma     <- fixed(log(2.28)); label("Exponent of the drug-effect power function (unitless)")                            # Lindauer 2017 Table S2 mouse PK/PD gamma = 2.28 (kept for human)
    psi        <- fixed(20);      label("Simeoni shape parameter (unitless) -- standard PSI = 20")                            # Simeoni 2004 PSI = 20 (cited by Lindauer 2017)

    # ----------------------------------------------------------------------
    # Inter-individual variability -- carried unchanged from the mouse
    # PK/PD fit (Lindauer 2017 Table S2). Only L0 and W0 have IIV;
    # omega^2 = log(CV^2 + 1).
    #   L0:  CV 59.4% -> omega^2 = log(0.594^2 + 1) = 0.30237
    #   W0:  CV 37.5% -> omega^2 = log(0.375^2 + 1) = 0.13164
    # ----------------------------------------------------------------------
    etall0    ~ 0.30237                                                                                                       # Lindauer 2017 Table S2 IIV L0 = 59.4% CV
    etalrbase ~ 0.13164                                                                                                       # Lindauer 2017 Table S2 IIV W0 = 37.5% CV

    # ----------------------------------------------------------------------
    # Residual error -- mouse-fit values reused as placeholders for human
    # simulations (the paper does not refit residuals on human data).
    # ----------------------------------------------------------------------
    propSd    <- fixed(0.197);  label("Proportional residual error on plasma pembrolizumab concentration (fraction; mouse-fit placeholder)")  # Lindauer 2017 Table S2 mouse PK proportional = 19.7% CV
    addSd     <- fixed(0.0658); label("Additive residual error on plasma pembrolizumab concentration (ug/mL; mouse-fit placeholder)")          # Lindauer 2017 Table S2 mouse PK additive = 0.0658 ug/mL
    propSd_tumor_vol <- fixed(0.20); label("Proportional residual error on tumor volume (fraction; placeholder)")                              # Not directly reported for human; placeholder ~Simeoni 2004 order
    expSd_R0_blood   <- fixed(0.627); label("Exponential residual error on blood receptor occupancy (CV; mouse-fit placeholder)")              # Lindauer 2017 Table S2 RO blood = 62.7% CV (exponential)
    expSd_R0_tumor   <- fixed(0.508); label("Exponential residual error on tumor receptor occupancy (CV; mouse-fit placeholder)")              # Lindauer 2017 Table S2 RO tumor = 50.8% CV (exponential)
  })

  model({
    # ============================================================
    # Internal unit conventions
    # ---------------------------------------------------------------
    # Time:                       day
    # Volumes:                    L
    # Antibody amounts:           nmol
    # Antibody concentrations:    nM (= nmol/L)
    # FcRn concentration:         nM (initial 49800 nM = 49.8 uM)
    # All Kon constants:          1/(nM*day)
    # All Koff and Kdeg:          1/day
    # PLQ_norm, CLup_norm:        1/day per L of tissue / endosomal vol
    # ============================================================

    # ---- Plasma PK structural parameters ----
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # ---- Tumor growth structural parameters (with IIV on L0 and W0) ----
    l0     <- exp(ll0 + etall0)
    l1     <- exp(ll1)
    rbase  <- exp(lrbase + etalrbase)
    sltg   <- exp(lsltg)
    gamma  <- exp(lgamma)

    # ---- Unit conversion factors ----
    # Convert mg dose to nmol via molecular weight: 1 mg = (1e6 / MW) nmol.
    mg_to_nmol <- 1e6 / mw_mab
    # Convert ug/mL to nM via molecular weight: 1 ug/mL = mg/L; mg/L / MW
    # (g/mol) * 1e6 = nM; equivalently 1e6 / MW = nM per ug/mL.
    ugmL_to_nM <- 1e6 / mw_mab
    # f() bioavailability so a "1 mg" dose lands as mg_to_nmol nmol in central.
    f(central) <- mg_to_nmol

    # ---- Tumor volume and the three Shah-Betts sub-volumes (in L) ----
    # tumor_vol is tracked in mL; multiply by 1e-3 to get L.
    tv_L  <- tumor_vol * 1e-3
    v_vs  <- f_v_vs * tv_L
    v_is  <- f_v_is * tv_L
    v_es  <- f_v_es * tv_L

    # ---- Tumor flow rates (L/day) ----
    plq_total  <- plq_norm  * tv_L           # plasma flow to tumor (L/day)
    l_total    <- f_lymph   * plq_total      # lymph flow (L/day)
    clup_total <- clup_norm * v_es           # endosomal pinocytosis (L/day)

    # ---- Antibody concentrations in nM (derived from amounts) ----
    C1   <- central     / vc
    C2   <- peripheral1 / vp
    Cvs  <- tumor_vs    / v_vs
    Cis  <- tumor_is    / v_is
    Ceub <- tumor_es_ub / v_es
    Ceb  <- tumor_es_b  / v_es

    # ---- FcRn conservation (FcRn_free + bound = fcrn_init total) ----
    # Total FcRn amount in endosome = fcrn_init * v_es (constant
    # conservation pool; equivalent to the Shah & Betts steady-state
    # assumption). Free FcRn = total minus the FcRn-bound antibody
    # amount per endosomal volume.
    fcrn_free <- fcrn_init - Ceb

    # ---- Pembrolizumab-PD-1 complex concentrations (nM) ----
    PD1b <- complex_blood / vc
    PD1t <- complex_tumor / v_is

    # ---- Total PD-1 receptor concentrations (nM) ----
    # Blood: constant. N_Tcell [cells/uL] * 1e6 [uL/L] * N_PD-1_TC
    # [receptors/cell] / 6.022e23 [/mol] * 1e9 [nM per M] collapses to
    # N_Tcell * N_PD-1_TC * 1.6605e-9 nM.
    C_PD1_b <- n_tcell * n_pd1_tc * 1.6605e-9
    # Tumor: dynamic, follows the target_tumor state amount in nmol
    # divided by the interstitial volume in L.
    C_PD1_t <- target_tumor / v_is

    # ---- Receptor occupancies (%) ----
    R0_blood <- 100 * PD1b / C_PD1_b
    R0_tumor <- 100 * PD1t / C_PD1_t

    # ---- Drug effect on tumor (Lindauer 2017 supplement) ----
    DE <- sltg * R0_tumor^gamma

    # ============================================================
    # ODE system. Equations follow Lindauer 2017 supplement
    # "Mathematical representation of the model as differential
    # equation system", with two implementation notes:
    #   (1) The supplement's central-compartment binding consumption
    #       term Kon * (C_PD1_b - PD1_b) is missing the C1 antibody
    #       multiplier and the V1 scaling; the canonical bimolecular
    #       mass balance vc * kon_pd1 * C1 * (C_PD1_b - PD1b) and
    #       dissociation vc * koff_pd1 * PD1b are used here.
    #   (2) FcRn is treated as a conserved species (FcRn_free + bound
    #       = constant) rather than carrying a separate FcRn ODE; this
    #       is the standard Shah & Betts implementation and avoids
    #       the supplement's apparent transcription typo in the
    #       dFcRn/dt equation.
    # All ODEs in nmol/day (= V * dC/dt from the supplement).
    # ============================================================

    # Plasma central (amount, nmol/day)
    d/dt(central) <-
      - cl   * C1                                                # linear elimination (L/day * nM = nmol/day)
      - vmax * mg_to_nmol * C1 / (km * ugmL_to_nM + C1)          # MM elim: Vmax in mg/day -> nmol/day, Km in ug/mL -> nM
      - q    * C1 + q * C2                                       # peripheral exchange
      - plq_total * C1 + plq_total * Cvs                         # plasma <-> tumor vascular exchange
      - vc * kon_pd1  * C1 * (C_PD1_b - PD1b)                    # PD-1 binding (consumes antibody)
      + vc * koff_pd1 * PD1b                                     # PD-1 dissociation

    # Plasma peripheral (amount, nmol/day)
    d/dt(peripheral1) <- q * C1 - q * C2

    # Tumor vascular space (amount, nmol/day; Lindauer 2017 supplement
    # "Vascular space tumor")
    d/dt(tumor_vs) <-
        plq_total                * C1
      - (plq_total - l_total)    * Cvs
      - (1 - v_ref)   * l_total  * Cvs
      - clup_total               * Cvs
      + clup_total * fr_recycle  * Ceb

    # Tumor endosomal unbound (amount, nmol/day; Lindauer 2017 supplement
    # "Endosomal space mAb unbound to FcRn")
    d/dt(tumor_es_ub) <-
        clup_total * (Cvs + Cis)
      - v_es * kon_fcrn  * Ceub * fcrn_free
      + v_es * koff_fcrn * Ceb
      - v_es * kdeg_endo * Ceub

    # Tumor endosomal FcRn-bound (amount, nmol/day; Lindauer 2017
    # supplement "Endosomal space mAb bound to FcRn", with the inflow
    # term corrected to a proper FcRn-binding mass balance and the
    # outflow CLup recycling explicit)
    d/dt(tumor_es_b) <-
        v_es * kon_fcrn  * Ceub * fcrn_free
      - v_es * koff_fcrn * Ceb
      - clup_total * Ceb

    # Tumor interstitial space (amount, nmol/day; Lindauer 2017
    # supplement "Interstitial compartment")
    d/dt(tumor_is) <-
        (1 - v_ref)    * l_total * Cvs
      - (1 - v_ref_is) * l_total * Cis
      - clup_total * Cis
      + clup_total * (1 - fr_recycle) * Ceb
      - v_is * kon_pd1  * Cis * (C_PD1_t - PD1t)
      + v_is * koff_pd1 * PD1t

    # Pembrolizumab-PD-1 complex in blood (amount, nmol/day; Lindauer
    # 2017 supplement "Drug receptor binding in blood")
    d/dt(complex_blood) <-
        vc * kon_pd1  * C1 * (C_PD1_b - PD1b)
      - vc * koff_pd1 * PD1b
      -     kdeg_pd1  * complex_blood

    # Pembrolizumab-PD-1 complex in tumor (amount, nmol/day; Lindauer
    # 2017 supplement "Drug receptor binding in the tumor")
    d/dt(complex_tumor) <-
        v_is * kon_pd1  * Cis * (C_PD1_t - PD1t)
      - v_is * koff_pd1 * PD1t
      -        kdeg_pd1 * complex_tumor

    # Total tumor PD-1 receptor amount (nmol; Lindauer 2017 supplement
    # "Tumor PD-1 receptor upregulation and elimination"). At PD1t = 0
    # (no drug) the steady-state amount equals MPD1_t(0). Setting
    # k_out = kdeg_pd1 and choosing k_in = kdeg_pd1 * MPD1_t(0) gives
    # the right baseline.
    kin_tumor <- kdeg_pd1 * tmulti * C_PD1_b * (rbase * 1e-3 * f_v_is)
    d/dt(target_tumor) <-
        kin_tumor * (1 + emax_tp * PD1t / (ec50_tp + PD1t))
      - kdeg_pd1  * target_tumor

    # Tumor volume (mL/day; Lindauer 2017 supplement "Tumor volume",
    # Simeoni-style exp-to-linear growth with the drug effect DE).
    d/dt(tumor_vol) <-
        l0 * tumor_vol / (1 + (l0 / l1 * tumor_vol)^psi)^(1 / psi)
      - DE * tumor_vol

    # ---- Initial conditions ----
    central(0)       <- 0
    peripheral1(0)   <- 0
    tumor_vs(0)      <- 0
    tumor_is(0)      <- 0
    tumor_es_ub(0)   <- 0
    tumor_es_b(0)    <- 0
    complex_blood(0) <- 0
    complex_tumor(0) <- 0
    target_tumor(0)  <- tmulti * C_PD1_b * (rbase * 1e-3 * f_v_is)
    tumor_vol(0)     <- rbase

    # ---- Observation variables and residual error ----
    # Plasma pembrolizumab concentration in ug/mL: central(nmol) / vc(L) / ugmL_to_nM
    Cc <- (central / vc) / ugmL_to_nM

    Cc        ~ add(addSd) + prop(propSd)
    tumor_vol ~ prop(propSd_tumor_vol)
    R0_blood  ~ lnorm(expSd_R0_blood)
    R0_tumor  ~ lnorm(expSd_R0_tumor)
  })
}
