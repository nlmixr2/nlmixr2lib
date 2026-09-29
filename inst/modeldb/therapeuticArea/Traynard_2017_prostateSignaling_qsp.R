Traynard_2017_prostateSignaling_qsp <- function() {
  description <- paste(
    "QSP. In vitro (LNCaP human prostate cancer cell line). Logic-based ODE",
    "model of the MAPK / PI3K / JAK-STAT / IKK signaling network that drives",
    "castration-resistant prostate cancer survival (Traynard 2017, the worked",
    "example of a CPT:PSP logic-modeling tutorial). Twenty network nodes are",
    "carried as ODE states holding a normalized activity on [0, 1]; five",
    "ligand stimuli (EGF, IGF-1, IL-6, TNF-alpha, DHT) and five kinase",
    "inhibitors (PI3K, MEK, mTOR, p38, IKK) enter as per-condition",
    "covariates. Each state relaxes toward the continuous (HillCube)",
    "homologue of its Boolean rule at a node-specific rate tau, with each",
    "regulatory edge passed through a normalized Hill transfer function of",
    "strength k. Multi-input nodes (AR, PI3K, RAS, JNK) use continuous OR.",
    "The 45 estimated parameters (20 tau, 25 k) were fitted with CellNOptR /",
    "CNORode2017 to 44 perturbation conditions of the Lescarbeau 2014 LNCaP",
    "phosphoproteome data at 30 and 240 min. All activities start at the",
    "normalized basal value 0.5; a fully inhibited node is held at basal.",
    "Deterministic: no IIV and no residual error are reported. The paper's",
    "companion MaBoSS continuous-time Boolean model (which adds the Survival,",
    "Cell_cycle, MYC, Caspase8, Caspase9, p53, NFkB and beta-catenin nodes)",
    "is a stochastic Markov process over Boolean states, not an ODE system,",
    "and is therefore not part of this file.",
    sep = " "
  )
  reference <- paste(
    "Traynard P, Tobalina L, Eduati F, Calzone L, Saez-Rodriguez J (2017).",
    "Logic Modeling in Quantitative Systems Pharmacology.",
    "CPT Pharmacometrics Syst Pharmacol 6(8):499-511.",
    "doi:10.1002/psp4.12225. PMCID: PMC5572374.",
    "The trained model is published only in the supplementary GitHub",
    "repository https://github.com/saezlab/CPT_QSPtutorial (GPLv3), cited in",
    "the paper's Supplementary Materials section; the main text prints the",
    "structural ODE form and two of the 25 edge strengths.",
    "Node lifetimes tau are supplement file logicODEparameters_nodes.txt;",
    "edge strengths k are supplement file logicODEparameters_edges.txt;",
    "the prior knowledge network is supplement file",
    "PriorKnowledgeNetwork.sif; the fitting script is supplement file",
    "CellNOptR_optimisation.R.",
    "Training data: Lescarbeau RM, Kaplan DL (2014). Quantitative analysis of",
    "castration resistant prostate cancer progression through phosphoproteome",
    "signaling. BMC Cancer 14:325. doi:10.1186/1471-2407-14-325.",
    "Transfer-function and right-hand-side definitions were read from the",
    "C sources of CNORode2017 (github.com/saezlab/CNORode2017), the package",
    "the supplement's fitting script installs and calls.",
    sep = " "
  )
  vignette <- "Traynard_2017_prostateSignaling_qsp"

  # Every state is a signaling-network node carrying a normalized protein
  # activity, not a drug amount, so none of them map onto a canonical PK
  # compartment role. Node symbols are the source network's (supplement
  # PriorKnowledgeNetwork.sif), which are standard protein symbols whose
  # capitalization is meaningful (p38, mTOR, GSK3a, Stat3).
  paper_specific_compartments <- c(
    "EGFR",
    "IGF1_R",
    "IL6R",
    "TNFR",
    "AR",
    "PI3K",
    "AKT",
    "mTOR",
    "RPS6",
    "GSK3a",
    "Jak",
    "Stat3",
    "RAS",
    "MEK",
    "ERK1_2",
    "Rac",
    "p38",
    "HSP27",
    "JNK",
    "IKKa"
  )

  units <- list(
    time = "min",
    dosing = paste(
      "(no dosing events; ligands and inhibitors are applied to the culture",
      "at time 0 and enter as the STIM_* and *_INHIBITED covariates)",
      sep = " "
    ),
    concentration = paste(
      "(normalized node activity on [0, 1]; 0.5 is the unperturbed basal",
      "state, i.e. a log2 fold change of zero)",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every state is a normalized activity of a protein (or
  # protein family) measured in, or inferred for, cultured LNCaP cells, so
  # the specimen is the cultured cell monolayer ("tissue"). The eight nodes
  # flagged `measured` in the analyte string are the phosphosites actually
  # assayed by Lescarbeau 2014 and used for training; the other twelve are
  # unmeasured network states.
  compartmentData <- list(
    EGFR = list(
      analyte = "EGF receptor activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    IGF1_R = list(
      analyte = "IGF-1 receptor activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    IL6R = list(
      analyte = "IL-6 receptor activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    TNFR = list(
      analyte = "TNF receptor activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    AR = list(
      analyte = "androgen receptor activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    PI3K = list(
      analyte = "phosphoinositide 3-kinase activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    AKT = list(
      analyte = "phospho-AKT (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    mTOR = list(
      analyte = "mTOR activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    RPS6 = list(
      analyte = "phospho-ribosomal protein S6 (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    GSK3a = list(
      analyte = "phospho-GSK3 alpha (measured phosphosite, normalized; the assayed site is inhibitory, so the sign of the GSK3 edges was inverted by the authors)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    Jak = list(
      analyte = "Janus kinase activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    Stat3 = list(
      analyte = "phospho-Stat3 (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    RAS = list(
      analyte = "RAS activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    MEK = list(
      analyte = "MEK activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    ERK1_2 = list(
      analyte = "phospho-ERK1/2 (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    Rac = list(
      analyte = "Rac GTPase activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    p38 = list(
      analyte = "phospho-p38 MAPK (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    HSP27 = list(
      analyte = "phospho-HSP27 (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    JNK = list(
      analyte = "phospho-JNK (measured phosphosite, normalized)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    ),
    IKKa = list(
      analyte = "IKK alpha activity (normalized, unmeasured network state)",
      units = "fraction of maximal activity",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    STIM_EGF_NORM = list(
      description = paste(
        "Normalized activity of the EGF input node. 0.5 = no EGF added (the",
        "unstimulated basal level that CellNOptR assigns to an unstimulated",
        "stimulus node), 1 = EGF applied to the culture.",
        sep = " "
      ),
      units = "fraction of maximal activity (unitless, [0, 1])",
      type = "continuous",
      reference_category = "0.5 (no EGF)",
      notes = paste(
        "Per-condition covariate, constant in time. The supplement fitting",
        "script sets unstimulated stimulus nodes to 0.5 rather than 0",
        "(CellNOptR_optimisation.R line",
        "'cnolist$valueStimuli[cnolist$valueStimuli==0]=0.5'), matching the",
        "log2-fold-change normalization in which 0.5 is the basal state. The",
        "applied ligand concentration is reported in the upstream data paper",
        "(Lescarbeau 2014), not in Traynard 2017, so only the binary",
        "applied / not-applied contrast is recoverable here and the covariate",
        "is expressed on the model's own normalized activity scale.",
        sep = " "
      ),
      source_name = "TR:EGF"
    ),
    STIM_IGF1_NORM = list(
      description = paste(
        "Normalized activity of the IGF-1 input node. 0.5 = no IGF-1 added,",
        "1 = IGF-1 applied to the culture.",
        sep = " "
      ),
      units = "fraction of maximal activity (unitless, [0, 1])",
      type = "continuous",
      reference_category = "0.5 (no IGF-1)",
      notes = paste(
        "Per-condition covariate, constant in time. Traynard 2017 notes that",
        "the fitted tau_IGF1_R is exactly 0, so IGF1_R is frozen at basal and",
        "this covariate has no effect on any state in the trained model",
        "('IGF1-R, TNFR, and IKK are associated with null tau parameters ...",
        "The system therefore does not depend on the inputs IGF-1 and TNF",
        "alpha'). It is retained because it is part of the published model",
        "structure.",
        sep = " "
      ),
      source_name = "TR:IGF_1"
    ),
    STIM_IL6_NORM = list(
      description = paste(
        "Normalized activity of the IL-6 input node. 0.5 = no IL-6 added,",
        "1 = IL-6 applied to the culture.",
        sep = " "
      ),
      units = "fraction of maximal activity (unitless, [0, 1])",
      type = "continuous",
      reference_category = "0.5 (no IL-6)",
      notes = paste(
        "Per-condition covariate, constant in time. tau_IL6R is very small",
        "(3.14e-04 per min), so the IL-6 arm moves only slightly over the",
        "240 min of the experiment.",
        sep = " "
      ),
      source_name = "TR:IL6"
    ),
    STIM_TNFA_NORM = list(
      description = paste(
        "Normalized activity of the TNF-alpha input node. 0.5 = no",
        "TNF-alpha added, 1 = TNF-alpha applied to the culture.",
        sep = " "
      ),
      units = "fraction of maximal activity (unitless, [0, 1])",
      type = "continuous",
      reference_category = "0.5 (no TNF-alpha)",
      notes = paste(
        "Per-condition covariate, constant in time. tau_TNFR is exactly 0 in",
        "the trained model, so TNFR is frozen at basal and this covariate has",
        "no effect on any state (Traynard 2017, paragraph following Table 1).",
        sep = " "
      ),
      source_name = "TR:TNFa"
    ),
    STIM_DHT_NORM = list(
      description = paste(
        "Normalized activity of the dihydrotestosterone (DHT) input node.",
        "0.5 = no DHT added, 1 = DHT applied to the culture.",
        sep = " "
      ),
      units = "fraction of maximal activity (unitless, [0, 1])",
      type = "continuous",
      reference_category = "0.5 (no DHT)",
      notes = paste(
        "Per-condition covariate, constant in time. DHT is the androgen",
        "stimulus whose survival advantage in LNCaP cells, and the abrogation",
        "of that advantage by PI3K inhibition, the paper reproduces.",
        sep = " "
      ),
      source_name = "TR:DHT"
    ),
    PI3K_INHIBITED = list(
      description = "1 = the culture was treated with a PI3K inhibitor, 0 = untreated.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PI3K inhibitor)",
      notes = paste(
        "Per-condition covariate, constant in time. In CNORode2017 an",
        "inhibitor multiplies the whole right-hand side of its target node by",
        "(1 - INH), so a fully inhibited node cannot move away from its",
        "initial value and is held at the normalized basal activity 0.5. The",
        "covariate is therefore usable on [0, 1] for a partial inhibition,",
        "although the Lescarbeau 2014 design used 0 / 1 only. The inhibitor",
        "compound identity and concentration are reported in the upstream data",
        "paper, not in Traynard 2017.",
        sep = " "
      ),
      source_name = "TR:PI3Ki"
    ),
    MEK_INHIBITED = list(
      description = "1 = the culture was treated with a MEK inhibitor, 0 = untreated.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no MEK inhibitor)",
      notes = "Per-condition covariate, constant in time. Acts on the MEK node exactly as PI3K_INHIBITED acts on PI3K.",
      source_name = "TR:MEKi"
    ),
    MTOR_INHIBITED = list(
      description = "1 = the culture was treated with an mTOR inhibitor, 0 = untreated.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mTOR inhibitor)",
      notes = "Per-condition covariate, constant in time. Acts on the mTOR node exactly as PI3K_INHIBITED acts on PI3K.",
      source_name = "TR:mTORi"
    ),
    P38_INHIBITED = list(
      description = "1 = the culture was treated with a p38 MAPK inhibitor, 0 = untreated.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no p38 inhibitor)",
      notes = "Per-condition covariate, constant in time. Acts on the p38 node exactly as PI3K_INHIBITED acts on PI3K.",
      source_name = "TR:p38i"
    ),
    IKK_INHIBITED = list(
      description = "1 = the culture was treated with an IKK inhibitor, 0 = untreated.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no IKK inhibitor)",
      notes = paste(
        "Per-condition covariate, constant in time. Acts on the IKKa node.",
        "tau_IKKa is exactly 0 in the trained model, so IKKa is frozen at",
        "basal regardless of this covariate.",
        sep = " "
      ),
      source_name = "TR:IKKai"
    )
  )

  population <- list(
    species = "in vitro (LNCaP human prostate cancer cell line)",
    n_subjects = 44,
    n_studies = 1,
    age_range = "(not applicable; cultured cell line)",
    weight_range = "(not applicable; cultured cell line)",
    sex_female_pct = NA_real_,
    disease_state = "castration-resistant prostate cancer (LNCaP, an androgen-sensitive line used to study the transition to androgen independence)",
    dose_range = paste(
      "44 perturbation conditions combining the five ligands",
      "(EGF, IGF-1, IL-6, TNF-alpha, DHT) with the five kinase inhibitors",
      "(PI3K, MEK, IKK, mTOR, p38). Ligand and inhibitor concentrations are",
      "reported in Lescarbeau 2014, not in Traynard 2017.",
      sep = " "
    ),
    regions = "(not applicable; cell culture)",
    n_conditions = 44,
    n_timepoints_used = "2 of 3 (30 min and 240 min; the 24 h timepoint was dropped because phosphorylation signaling is expected to reach a semi-steady state within a few hours and a 24 h window would require modeling transcriptional rewiring)",
    n_measured_nodes = 8,
    notes = paste(
      "The 'subject' in this model is a perturbation condition, not an",
      "organism. Training data are the LNCaP arm of Lescarbeau 2014, an",
      "antibody-based phosphoproteomic panel of eight phosphosites",
      "(AKT, RPS6, GSK3, ERK1/2, p38, JNK, HSP27, Stat3), normalized per",
      "species as the log2 fold change versus the unperturbed basal state and",
      "then linearly rescaled onto [0, 1] with 0.5 at basal. Traynard 2017",
      "excluded the docetaxel arm of Lescarbeau 2014 because it produced",
      "little phosphoproteome variation and its target (beta-tubulin) is",
      "outside the network. Optimization was repeated 10 times from random",
      "starts (coefficient of variation across runs 0.009); the best run is",
      "the parameter set carried here. Bootstrap (300x), data randomization",
      "(300x) and network randomization (100x) all showed the trained model",
      "performs significantly better than chance (Figure 3d).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Logic-based ODE model in the formalism of Wittmann 2009, as
    # implemented in CellNOptR / CNORode2017 (Terfve 2012). For each state
    # node j,
    #
    #   dx_j/dt = tau_j * ( B_j(f(x_i1, k_i1->j), ...) - x_j ) * (1 - inh_j)
    #
    # where f() is the edge transfer function, B_j is the continuous
    # (multilinear / HillCube) homologue of node j's Boolean rule, tau_j is
    # the node's inverse lifetime, and inh_j is the inhibitor covariate for
    # that node. Traynard 2017 prints the AKT instance of this equation in
    # the section "Training the logic-based ODE model with CellNOpt"; the
    # general right-hand side is CNORode2017 src/rhsODE.c.
    #
    # The supplement's fitting script (CellNOptR_optimisation.R) selects
    # transfer_function = 4, which CNORode2017 src/sim_logic_ode.c dispatches
    # to FG_transfer_function (src/FG_transfer_function.c):
    #
    #   f(x, n, k) = 1 - (1 - x)^n / ( (1 - x)^n + k^n ) * ( 1 + k^n )
    #
    # This is a monotone sigmoid on [0, 1] with f(0) = 0 and f(1) = 1; the
    # paper describes it as "a monotonically increasing sigmoid in the
    # variable ... (in the range {0,1})" whose "increase rate depends on the
    # parameter k ... defining the strength of the regulatory interaction",
    # and notes that k = 0 makes the target independent of that regulator
    # (f is identically 0 there), which is exactly this expression.
    #
    # SCALE. tau and k are bounded on [0, 1] by the fitting script
    # (LB_tau = 0, UB_tau = 1, LB_k = 0, UB_k = 1) and several estimates sit
    # exactly on a bound (tau_TNFR = tau_IGF1_R = tau_IKKa = 0; six k = 1).
    # They are therefore carried on the natural scale, NOT log-transformed:
    # log(0) is undefined and a log transform is the wrong parameterization
    # for a parameter that is bounded above.
    #
    # ESTIMATED vs FIXED. tau and k were estimated (opt_tau = TRUE,
    # opt_k = TRUE) by enhanced scatter search (MEIGO essm) against the
    # 44-condition LNCaP dataset, so they are NOT wrapped in fixed(). The
    # Hill exponent n was held at its default (default_n = 3,
    # opt_n = FALSE) and the basal normalization level is a data-processing
    # constant; both are fixed().
    #
    # No standard errors, confidence intervals, IIV or residual error are
    # reported anywhere in the paper or its supplement: this is a single
    # deterministic fit whose uncertainty is reported graphically as a
    # bootstrap performance distribution (Figure 3d), not as per-parameter
    # intervals.
    # =====================================================================

    nhill <- fixed(3) ; label("Hill exponent n of every edge transfer function (unitless)") # supplement CellNOptR_optimisation.R: createLBodeContPars(..., default_n = 3, opt_n = FALSE)
    basal <- fixed(0.5) ; label("Normalized basal node activity, used as the initial condition of every state (unitless)") # Traynard 2017 'Training the logic-based ODE model with CellNOpt' step 3: 'all initial conditions are set to 0.5 (which is the basal state, or the state at time 0)'

    # ---------------------------------------------------------------------
    # Node lifetimes tau (1/min). Supplement logicODEparameters_nodes.txt.
    # The same 20 values appear again as the MaBoSS activation rates
    # $u_<node> in supplement trainedmodel.cfg, which is an independent
    # transcription of the same fit (the paper states that the tau values
    # 'obtained in the optimized logic-based ODE model' were assigned to the
    # MaBoSS node-activation transition rates).
    # ---------------------------------------------------------------------
    tau_EGFR <- 1 ; label("Lifetime parameter tau of EGFR (1/min)") # logicODEparameters_nodes.txt tau_EGFR; cfg $u_EGFR = 1
    tau_IGF1_R <- 0 ; label("Lifetime parameter tau of IGF1_R (1/min); zero, so IGF1_R is frozen at basal") # logicODEparameters_nodes.txt tau_IGF1_R; cfg $u_IGF1_R = 0
    tau_IL6R <- 0.000313976514259364 ; label("Lifetime parameter tau of IL6R (1/min)") # logicODEparameters_nodes.txt tau_IL6R; cfg $u_IL6R
    tau_TNFR <- 0 ; label("Lifetime parameter tau of TNFR (1/min); zero, so TNFR is frozen at basal") # logicODEparameters_nodes.txt tau_TNFR; cfg $u_TNFR = 0
    tau_AR <- 0.0372178932277208 ; label("Lifetime parameter tau of AR (1/min)") # logicODEparameters_nodes.txt tau_AR; cfg $u_AR
    tau_PI3K <- 0.0021076346277392 ; label("Lifetime parameter tau of PI3K (1/min)") # logicODEparameters_nodes.txt tau_PI3K; cfg $u_PI3K
    tau_AKT <- 0.0142451029797276 ; label("Lifetime parameter tau of AKT (1/min)") # logicODEparameters_nodes.txt tau_AKT; cfg $u_AKT
    tau_mTOR <- 0.0262595544549834 ; label("Lifetime parameter tau of mTOR (1/min)") # logicODEparameters_nodes.txt tau_mTOR; cfg $u_mTOR
    tau_RPS6 <- 0.023113503254164 ; label("Lifetime parameter tau of RPS6 (1/min)") # logicODEparameters_nodes.txt tau_RPS6; cfg $u_RPS6
    tau_GSK3a <- 0.0377782342041935 ; label("Lifetime parameter tau of GSK3a (1/min)") # logicODEparameters_nodes.txt tau_GSK3a; cfg $u_GSK3a
    tau_Jak <- 0.00801935194391314 ; label("Lifetime parameter tau of Jak (1/min)") # logicODEparameters_nodes.txt tau_Jak; cfg $u_Jak
    tau_Stat3 <- 0.00481896239850702 ; label("Lifetime parameter tau of Stat3 (1/min)") # logicODEparameters_nodes.txt tau_Stat3; cfg $u_Stat3
    tau_RAS <- 0.000929779281639963 ; label("Lifetime parameter tau of RAS (1/min)") # logicODEparameters_nodes.txt tau_RAS; cfg $u_RAS
    tau_MEK <- 0.0866293601256323 ; label("Lifetime parameter tau of MEK (1/min)") # logicODEparameters_nodes.txt tau_MEK; cfg $u_MEK
    tau_ERK1_2 <- 0.0277052045529414 ; label("Lifetime parameter tau of ERK1_2 (1/min)") # logicODEparameters_nodes.txt tau_ERK1_2; cfg $u_ERK1_2
    tau_Rac <- 0.0603798602331782 ; label("Lifetime parameter tau of Rac (1/min)") # logicODEparameters_nodes.txt tau_Rac; cfg $u_Rac
    tau_p38 <- 0.11875303258739 ; label("Lifetime parameter tau of p38 (1/min)") # logicODEparameters_nodes.txt tau_p38; cfg $u_p38
    tau_HSP27 <- 0.0715007312382617 ; label("Lifetime parameter tau of HSP27 (1/min)") # logicODEparameters_nodes.txt tau_HSP27; cfg $u_HSP27
    tau_JNK <- 0.015887038493103 ; label("Lifetime parameter tau of JNK (1/min)") # logicODEparameters_nodes.txt tau_JNK; cfg $u_JNK
    tau_IKKa <- 0 ; label("Lifetime parameter tau of IKKa (1/min); zero, so IKKa is frozen at basal") # logicODEparameters_nodes.txt tau_IKKa; cfg $u_IKKa = 0

    # ---------------------------------------------------------------------
    # Edge strengths k (unitless). Supplement logicODEparameters_edges.txt,
    # whose naming scheme is '<source>_k_<target>'. Twenty-five edges: the
    # 41-edge prior knowledge network minus the non-observable and
    # non-controllable branch removed by CellNOptR's compression step
    # (Cell_cycle, MYC, Caspase8, Caspase9, p53, NFkB, beta_catenin, the two
    # AND gates and the Stress input), with the compressed AND gate
    # 'and2 = TNFR AND NOT p53' collapsing to the single edge TNFR -> IKKa.
    # ---------------------------------------------------------------------
    k_EGF_EGFR <- 0.550337329553758 ; label("Edge strength k of EGF -> EGFR (unitless)") # logicODEparameters_edges.txt EGF_k_EGFR
    k_IGF1_IGF1R <- 1 ; label("Edge strength k of IGF-1 -> IGF1_R (unitless)") # logicODEparameters_edges.txt IGF_1_k_IGF1_R
    k_IL6_IL6R <- 0.78134194522735 ; label("Edge strength k of IL-6 -> IL6R (unitless)") # logicODEparameters_edges.txt IL6_k_IL6R
    k_TNFa_TNFR <- 1 ; label("Edge strength k of TNF-alpha -> TNFR (unitless)") # logicODEparameters_edges.txt TNFa_k_TNFR
    k_AKT_AR <- 0.000835340126624366 ; label("Edge strength k of AKT -> AR (unitless)") # logicODEparameters_edges.txt AKT_k_AR; also printed in the main text as 'kAKT -> AR = 8.4e-04'
    k_DHT_AR <- 0.590287658825809 ; label("Edge strength k of DHT -> AR (unitless)") # logicODEparameters_edges.txt DHT_k_AR
    k_AR_PI3K <- 1 ; label("Edge strength k of AR -> PI3K (unitless)") # logicODEparameters_edges.txt AR_k_PI3K
    k_EGFR_PI3K <- 1 ; label("Edge strength k of EGFR -> PI3K (unitless)") # logicODEparameters_edges.txt EGFR_k_PI3K
    k_IGF1R_PI3K <- 1 ; label("Edge strength k of IGF1_R -> PI3K (unitless)") # logicODEparameters_edges.txt IGF1_R_k_PI3K
    k_PI3K_AKT <- 0.352193224711085 ; label("Edge strength k of PI3K -> AKT (unitless)") # logicODEparameters_edges.txt PI3K_k_AKT
    k_AKT_mTOR <- 1 ; label("Edge strength k of AKT -> mTOR (unitless)") # logicODEparameters_edges.txt AKT_k_mTOR
    k_mTOR_RPS6 <- 0.303674761072917 ; label("Edge strength k of mTOR -> RPS6 (unitless)") # logicODEparameters_edges.txt mTOR_k_RPS6
    k_AKT_GSK3a <- 0.558225196843653 ; label("Edge strength k of AKT -> GSK3a (unitless)") # logicODEparameters_edges.txt AKT_k_GSK3a
    k_IL6R_Jak <- 1 ; label("Edge strength k of IL6R -> Jak (unitless)") # logicODEparameters_edges.txt IL6R_k_Jak
    k_Jak_Stat3 <- 0.28320076784032 ; label("Edge strength k of Jak -> Stat3 (unitless)") # logicODEparameters_edges.txt Jak_k_Stat3
    k_Jak_RAS <- 0.674347577514922 ; label("Edge strength k of Jak -> RAS (unitless)") # logicODEparameters_edges.txt Jak_k_RAS
    k_EGFR_RAS <- 0.25869659979145 ; label("Edge strength k of EGFR -> RAS (unitless)") # logicODEparameters_edges.txt EGFR_k_RAS
    k_RAS_MEK <- 0.490290337695327 ; label("Edge strength k of RAS -> MEK (unitless)") # logicODEparameters_edges.txt RAS_k_MEK
    k_MEK_ERK1_2 <- 0.462172528684214 ; label("Edge strength k of MEK -> ERK1_2 (unitless)") # logicODEparameters_edges.txt MEK_k_ERK1_2
    k_RAS_Rac <- 0.284328823821596 ; label("Edge strength k of RAS -> Rac (unitless)") # logicODEparameters_edges.txt RAS_k_Rac
    k_Rac_p38 <- 1 ; label("Edge strength k of Rac -> p38 (unitless)") # logicODEparameters_edges.txt Rac_k_p38
    k_p38_HSP27 <- 0.993983089528085 ; label("Edge strength k of p38 -> HSP27 (unitless)") # logicODEparameters_edges.txt p38_k_HSP27
    k_TNFR_JNK <- 0.361795499785572 ; label("Edge strength k of TNFR -> JNK (unitless)") # logicODEparameters_edges.txt TNFR_k_JNK
    k_Rac_JNK <- 0.0000673323536405412 ; label("Edge strength k of Rac -> JNK (unitless)") # logicODEparameters_edges.txt Rac_k_JNK; also printed in the main text as 'kRac -> JNK = 6.7e-05'
    k_TNFR_IKKa <- 0.554961626612893 ; label("Edge strength k of TNFR -> IKKa (unitless)") # logicODEparameters_edges.txt TNFR_k_IKKa
  })

  model({
    # =====================================================================
    # Edge transfer functions. CNORode2017 transfer_function = 4, i.e.
    # src/FG_transfer_function.c:
    #   f(x, n, k) = 1 - (1 - x)^n / ( (1 - x)^n + k^n ) * ( 1 + k^n )
    # evaluated at the CURRENT activity of the source node (or, for the five
    # ligand inputs, at the constant stimulus covariate).
    # =====================================================================

    # --- Receptor layer: one stimulus covariate into one receptor node ----
    fEgfEgfr <- 1 - (1 - STIM_EGF_NORM)^nhill / ((1 - STIM_EGF_NORM)^nhill + k_EGF_EGFR^nhill) * (1 + k_EGF_EGFR^nhill)
    fIgf1Igf1r <- 1 - (1 - STIM_IGF1_NORM)^nhill / ((1 - STIM_IGF1_NORM)^nhill + k_IGF1_IGF1R^nhill) * (1 + k_IGF1_IGF1R^nhill)
    fIl6Il6r <- 1 - (1 - STIM_IL6_NORM)^nhill / ((1 - STIM_IL6_NORM)^nhill + k_IL6_IL6R^nhill) * (1 + k_IL6_IL6R^nhill)
    fTnfaTnfr <- 1 - (1 - STIM_TNFA_NORM)^nhill / ((1 - STIM_TNFA_NORM)^nhill + k_TNFa_TNFR^nhill) * (1 + k_TNFa_TNFR^nhill)

    # --- Androgen receptor: OR(AKT, DHT) ---------------------------------
    fAktAr <- 1 - (1 - AKT)^nhill / ((1 - AKT)^nhill + k_AKT_AR^nhill) * (1 + k_AKT_AR^nhill)
    fDhtAr <- 1 - (1 - STIM_DHT_NORM)^nhill / ((1 - STIM_DHT_NORM)^nhill + k_DHT_AR^nhill) * (1 + k_DHT_AR^nhill)

    # --- PI3K: OR(AR, EGFR, IGF1_R) --------------------------------------
    fArPi3k <- 1 - (1 - AR)^nhill / ((1 - AR)^nhill + k_AR_PI3K^nhill) * (1 + k_AR_PI3K^nhill)
    fEgfrPi3k <- 1 - (1 - EGFR)^nhill / ((1 - EGFR)^nhill + k_EGFR_PI3K^nhill) * (1 + k_EGFR_PI3K^nhill)
    fIgf1rPi3k <- 1 - (1 - IGF1_R)^nhill / ((1 - IGF1_R)^nhill + k_IGF1R_PI3K^nhill) * (1 + k_IGF1R_PI3K^nhill)

    # --- Single-input edges of the PI3K / mTOR and JAK-STAT arms ---------
    fPi3kAkt <- 1 - (1 - PI3K)^nhill / ((1 - PI3K)^nhill + k_PI3K_AKT^nhill) * (1 + k_PI3K_AKT^nhill)
    fAktMtor <- 1 - (1 - AKT)^nhill / ((1 - AKT)^nhill + k_AKT_mTOR^nhill) * (1 + k_AKT_mTOR^nhill)
    fMtorRps6 <- 1 - (1 - mTOR)^nhill / ((1 - mTOR)^nhill + k_mTOR_RPS6^nhill) * (1 + k_mTOR_RPS6^nhill)
    fAktGsk3a <- 1 - (1 - AKT)^nhill / ((1 - AKT)^nhill + k_AKT_GSK3a^nhill) * (1 + k_AKT_GSK3a^nhill)
    fIl6rJak <- 1 - (1 - IL6R)^nhill / ((1 - IL6R)^nhill + k_IL6R_Jak^nhill) * (1 + k_IL6R_Jak^nhill)
    fJakStat3 <- 1 - (1 - Jak)^nhill / ((1 - Jak)^nhill + k_Jak_Stat3^nhill) * (1 + k_Jak_Stat3^nhill)

    # --- RAS: OR(Jak, EGFR) ----------------------------------------------
    fJakRas <- 1 - (1 - Jak)^nhill / ((1 - Jak)^nhill + k_Jak_RAS^nhill) * (1 + k_Jak_RAS^nhill)
    fEgfrRas <- 1 - (1 - EGFR)^nhill / ((1 - EGFR)^nhill + k_EGFR_RAS^nhill) * (1 + k_EGFR_RAS^nhill)

    # --- MAPK and stress arms --------------------------------------------
    fRasMek <- 1 - (1 - RAS)^nhill / ((1 - RAS)^nhill + k_RAS_MEK^nhill) * (1 + k_RAS_MEK^nhill)
    fMekErk <- 1 - (1 - MEK)^nhill / ((1 - MEK)^nhill + k_MEK_ERK1_2^nhill) * (1 + k_MEK_ERK1_2^nhill)
    fRasRac <- 1 - (1 - RAS)^nhill / ((1 - RAS)^nhill + k_RAS_Rac^nhill) * (1 + k_RAS_Rac^nhill)
    fRacP38 <- 1 - (1 - Rac)^nhill / ((1 - Rac)^nhill + k_Rac_p38^nhill) * (1 + k_Rac_p38^nhill)
    fP38Hsp27 <- 1 - (1 - p38)^nhill / ((1 - p38)^nhill + k_p38_HSP27^nhill) * (1 + k_p38_HSP27^nhill)

    # --- JNK: OR(TNFR, Rac); IKKa: single input --------------------------
    fTnfrJnk <- 1 - (1 - TNFR)^nhill / ((1 - TNFR)^nhill + k_TNFR_JNK^nhill) * (1 + k_TNFR_JNK^nhill)
    fRacJnk <- 1 - (1 - Rac)^nhill / ((1 - Rac)^nhill + k_Rac_JNK^nhill) * (1 + k_Rac_JNK^nhill)
    fTnfrIkka <- 1 - (1 - TNFR)^nhill / ((1 - TNFR)^nhill + k_TNFR_IKKa^nhill) * (1 + k_TNFR_IKKa^nhill)

    # =====================================================================
    # Continuous Boolean targets B_j. For a single activating input the
    # HillCube target is just the transfer value; for the four multi-input
    # nodes the Boolean rule is OR (Traynard 2017 Figure 2b caption: "All
    # other nodes in the model with more than one input edge are modeled
    # with a simple OR gate"), whose multilinear homologue is
    # 1 - prod(1 - f_i). Neither of the prior network's two AND gates
    # survives compression, so no product term is needed here.
    # =====================================================================
    bAr <- 1 - (1 - fAktAr) * (1 - fDhtAr)
    bPi3k <- 1 - (1 - fArPi3k) * (1 - fEgfrPi3k) * (1 - fIgf1rPi3k)
    bRas <- 1 - (1 - fJakRas) * (1 - fEgfrRas)
    bJnk <- 1 - (1 - fTnfrJnk) * (1 - fRacJnk)

    # =====================================================================
    # ODE system. Each node relaxes toward its continuous Boolean target at
    # rate tau. The five inhibitor covariates multiply the whole right-hand
    # side of their target node (CNORode2017 src/rhsODE.c, factor
    # (1 - inhibitor_array[j])), so a fully inhibited node cannot leave its
    # initial condition.
    # =====================================================================
    d/dt(EGFR) <- tau_EGFR * (fEgfEgfr - EGFR)
    d/dt(IGF1_R) <- tau_IGF1_R * (fIgf1Igf1r - IGF1_R)
    d/dt(IL6R) <- tau_IL6R * (fIl6Il6r - IL6R)
    d/dt(TNFR) <- tau_TNFR * (fTnfaTnfr - TNFR)
    d/dt(AR) <- tau_AR * (bAr - AR)
    d/dt(PI3K) <- tau_PI3K * (bPi3k - PI3K) * (1 - PI3K_INHIBITED)
    d/dt(AKT) <- tau_AKT * (fPi3kAkt - AKT)
    d/dt(mTOR) <- tau_mTOR * (fAktMtor - mTOR) * (1 - MTOR_INHIBITED)
    d/dt(RPS6) <- tau_RPS6 * (fMtorRps6 - RPS6)
    d/dt(GSK3a) <- tau_GSK3a * (fAktGsk3a - GSK3a)
    d/dt(Jak) <- tau_Jak * (fIl6rJak - Jak)
    d/dt(Stat3) <- tau_Stat3 * (fJakStat3 - Stat3)
    d/dt(RAS) <- tau_RAS * (bRas - RAS)
    d/dt(MEK) <- tau_MEK * (fRasMek - MEK) * (1 - MEK_INHIBITED)
    d/dt(ERK1_2) <- tau_ERK1_2 * (fMekErk - ERK1_2)
    d/dt(Rac) <- tau_Rac * (fRasRac - Rac)
    d/dt(p38) <- tau_p38 * (fRacP38 - p38) * (1 - P38_INHIBITED)
    d/dt(HSP27) <- tau_HSP27 * (fP38Hsp27 - HSP27)
    d/dt(JNK) <- tau_JNK * (bJnk - JNK)
    d/dt(IKKa) <- tau_IKKa * (fTnfrIkka - IKKa) * (1 - IKK_INHIBITED)

    # =====================================================================
    # Initial conditions. Every state starts at the normalized basal
    # activity (CNORode2017 src/simulateODE.c initializes the whole state
    # array to 0.5, and the MIDAS training file carries 0.5 for every
    # readout at time 0).
    # =====================================================================
    EGFR(0) <- basal
    IGF1_R(0) <- basal
    IL6R(0) <- basal
    TNFR(0) <- basal
    AR(0) <- basal
    PI3K(0) <- basal
    AKT(0) <- basal
    mTOR(0) <- basal
    RPS6(0) <- basal
    GSK3a(0) <- basal
    Jak(0) <- basal
    Stat3(0) <- basal
    RAS(0) <- basal
    MEK(0) <- basal
    ERK1_2(0) <- basal
    Rac(0) <- basal
    p38(0) <- basal
    HSP27(0) <- basal
    JNK(0) <- basal
    IKKa(0) <- basal

    # =====================================================================
    # Observations. The eight trained readouts are the states themselves --
    # the model is written on the same normalized [0, 1] scale as the
    # MIDAS training data, so no transformation is applied. No residual
    # error is reported: Traynard 2017 minimizes a residual sum of squares
    # plus a steady-state penalty and reports no error model, no IIV and no
    # parameter standard errors.
    # =====================================================================
  })
}
