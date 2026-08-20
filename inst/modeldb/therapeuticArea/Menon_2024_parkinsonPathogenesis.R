Menon_2024_parkinsonPathogenesis <- function() {
  description <- paste0(
    "QSP. Combined ROS + proteasome-sequestration mechanistic model of ",
    "Parkinson's disease pathogenesis in dopaminergic neurons, from the ",
    "Menon-Bakshi-Krishnan 2024 CPT:PSP Perspective. Five ODEs describe ",
    "the intracellular time-evolution (hours) of misfolded alpha-synuclein ",
    "monomer (alpha_syn), the larger alpha-syn aggregate species F (fibrils ",
    "and oligomers lumped, fagg), free proteasome components (pfree), the ",
    "proteasome-aggregate complex CP (pcomplex), and cytoplasmic reactive ",
    "oxygen species (ros). The ROS module (positive feedback: misfolded ",
    "alpha-syn drives cytoplasmic ROS via mitochondrial + DA-vesicle leak; ",
    "ROS drives further alpha-syn misfolding) is combined with the ",
    "proteasomal module (mutual inhibition: proteasome degrades aggregates, ",
    "aggregates sequester proteasome). Combining these two nonlinear ",
    "feedback modules yields the emergent landscape mapped in the paper: ",
    "monostable or bistable alpha-syn (from the ROS positive feedback) x ",
    "steady-state or oscillatory proteasome (supercritical or subcritical ",
    "Hopf, from the proteasome negative feedback). All concentrations are ",
    "non-dimensionalised (scale factor 50 pM); time is in hours. The basal ",
    "cytoplasmic ROS production rate k0 is the bifurcation-input parameter ",
    "the paper varies to demonstrate pulse- and step-triggered switching to ",
    "the diseased state (elevated misfolded alpha-syn, sustained or ",
    "transient proteasome-sequestration oscillations, timescale ~1 year). ",
    "Default ini() values are the Figure 2(A) scenario (supercritical Hopf ",
    "in the proteasomal module, no overlap with the bistable regime); the ",
    "vignette walks through Figures 1 (individual modules), 2 (combined ",
    "landscape), and 2(I) / S2 (proteasome-availability intervention). The ",
    "proteasome-degradation rate kP is not tabulated in the paper's ",
    "parameter tables and is inferred here from the steady-state balance ",
    "at healthy [P]_ss = 1 (see vignette Assumptions and deviations). This ",
    "is a mechanistic ODE model without drug input, without IIV, and ",
    "without residual error (typical-value simulation). A companion DA ",
    "compartmentalization variant is packaged as ",
    "Menon_2024_parkinsonPathogenesis_DA."
  )
  reference <- paste(
    "Menon G, Bakshi S, Krishnan J. (2024).",
    "The interaction of core modules as a basis for elucidating network",
    "behavior determining Parkinson's disease pathogenesis.",
    "CPT Pharmacometrics Syst Pharmacol 13(3):335-340.",
    "doi:10.1002/psp4.13108.",
    "Structural components inherited from Cloutier & Wellstead (2012)",
    "IET Syst Biol 6(3):86-93 (ROS positive-feedback module) and Sneppen",
    "et al. (2009) Phys Biol 6:036005 (proteasome sequestration module).",
    sep = " "
  )
  vignette <- "Menon_2024_parkinsonPathogenesis"
  paper_specific_compartments <- c("alpha_syn", "fagg", "pfree", "pcomplex", "ros")

  units <- list(
    time          = "h",
    dosing        = "(none; mechanistic disease-pathogenesis model, no drug input)",
    concentration = "(dimensionless; all state variables scaled by 50 pM per Supinfo2 Section 1.4)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    fagg      = list(analyte = "alpha-syn aggregate species", units = NA_character_, specimen = "administration site", verified = FALSE),
    pcomplex  = list(analyte = "proteasome-aggregate complex", units = NA_character_, specimen = "administration site", verified = FALSE),
    pfree     = list(analyte = "free proteasome components", units = NA_character_, specimen = "administration site", verified = FALSE),
    alpha_syn = list(analyte = "misfolded alpha-synuclein monomer", units = NA_character_, specimen = "administration site", verified = FALSE),
    ros       = list(analyte = "cytoplasmic reactive oxygen species", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "in silico (dopaminergic neuron cell-level mechanistic model)",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "(not applicable; typical-value cellular mechanistic model, not fit to individual subjects)",
    weight_range   = "(not applicable)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Parkinson's disease pathogenesis (dopaminergic-neuron intracellular network). The healthy vs. diseased distinction is a bistable-switch state of the model itself: at low k0 the system rests in a low-alpha_syn healthy attractor; a pulse or step increase in k0 (or in other network parameters) triggers a switch to a high-alpha_syn diseased attractor characterised by elevated misfolded alpha-synuclein and (depending on the proteasomal-module parameters) sustained proteasome-sequestration oscillations.",
    dose_range     = "(not applicable)",
    regions        = "(not applicable)",
    notes          = "Menon 2024 is explicitly labelled a 'Perspective' paper. It presents an original combined ODE system assembled from two published upstream modules -- Cloutier & Wellstead 2012 (ROS positive-feedback module) and Sneppen et al. 2009 (proteasome sequestration module) -- and analyses the emergent systems landscape (bistability, Hopf bifurcations, pulse/step-triggered switching, intervention reversibility). Parameter values in Supinfo3 span 10 figure scenarios (Fig 1B, 1CD, 2A-I, S1-S8) that vary k3, kagg, krem, sigma, gamma, nu, and m to expose different qualitative behaviours. The Fig 2(A) values used here as the default correspond to the supercritical-Hopf-no-bistable-overlap scenario used in the paper's main-text analysis of pulse-triggered switching. The paper reports simulations in MATLAB via ode15s; MATLAB code is at https://github.com/gm1613/parkinsons_disease_models. This nlmixr2lib encoding uses rxode2 as the ODE engine; qualitative behaviours (bistable switch, sustained oscillations, intervention reversal) are reproduced in the validation vignette. All concentrations non-dimensionalised (scale factor 50 pM per Supinfo2 Section 1.4)."
  )

  ini({
    # ------------------------------------------------------------------
    # Default parameter set: Figure 2(A) scenario (Supinfo3 Table
    # 'Figure 2 (A), Figure S1(C,D)'). This scenario is supercritical
    # Hopf in the proteasomal module with NO overlap between the
    # bistable and oscillatory regimes -- a pulse in k0 triggers a
    # switch to the diseased state with transient oscillations; a step
    # in k0 triggers a switch to the diseased state with sustained
    # low-amplitude oscillations (paper Fig S1(C,D)).
    #
    # k0 (basal cytoplasmic ROS production) is the paper's bifurcation
    # parameter -- Supinfo3 does not tabulate a specific k0 value, but
    # rather varies k0 across each figure's bifurcation diagram. The
    # value 0.001 (dimensionless) placed here is chosen small enough
    # to hold the system in the low-alpha_syn healthy attractor at
    # steady state (see vignette 'Healthy baseline' chunk); the
    # vignette also shows the diseased-state simulation obtained by
    # stepping k0 up.
    #
    # All rate constants are h^-1 unless dimensionless (d_asyn,
    # K_asyn, k0 -- see comments). Concentrations are
    # non-dimensionalised (scale factor 50 pM per Supinfo2 Sec 1.4).
    # ------------------------------------------------------------------

    k0     <- 0.001;   label("Basal cytoplasmic ROS production rate (dimensionless bifurcation input; see vignette 'Bifurcation input k0')")
    # Supinfo2 Eq (5): k1*(k0 + d_asyn*(alpha_syn^4/(1+alpha_syn^4))) is the total ROS production term; k0 is the alpha_syn-independent component. Supinfo3 does not tabulate a specific k0 -- it is the bifurcation-input parameter varied per figure. Default 0.001 chosen small enough to hold the system in the healthy attractor.

    m      <- 0.00042; label("Basal misfolded alpha-synuclein formation rate (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' m = 0.00042 h^-1. Supinfo2 Eq (4): d[alpha_syn]/dt has +m as constant basal formation.

    k1     <- 35;      label("Cytoplasmic ROS input scaling rate constant (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' k1 = 35 h^-1. Constant across all scenarios (Cloutier 2012 module inheritance).

    k2     <- 35;      label("Cytoplasmic ROS removal rate constant (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' k2 = 35 h^-1. Constant across all scenarios.

    k3     <- 0.0087;  label("Rate constant of alpha-synuclein misfolding driven by cytoplasmic ROS (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' k3 = 0.0087 h^-1. Varies across scenarios (0.0047-0.0087) as the ROS-driven misfolding strength that modulates bistable/oscillatory regimes.

    k4     <- 0.0042;  label("Rate constant of misfolded alpha-synuclein removal by proteolysis (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' k4 = 0.0042 h^-1. Supinfo2 Eq (4): -k4*[P]*[alpha_syn] is the proteolytic removal term. NOTE the isolated ROS module in Fig 1(B) uses k4 = 0.0087 (h^-1) instead; see vignette 'ROS module in isolation' for the isolated-module reproduction.

    kagg   <- 0.0042;  label("Rate constant of formation of larger alpha-synuclein aggregates (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' kagg = 0.0042 h^-1. Supinfo2 Eq (1)/(4): kagg*[alpha_syn] transfers monomer into the aggregate compartment F.

    krem   <- 0.0045;  label("Background removal rate of misfolded alpha-synuclein by other mechanisms (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' krem = 0.0045 h^-1. Supinfo2 Eq (4): -krem*[alpha_syn] is a lumped removal term for lysosomal degradation, cell leakage, and any non-proteasomal pathway.

    sigma  <- 0.0042;  label("Background production rate of free proteasome components (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' sigma = 0.0042 h^-1. Supinfo2 Eq (3): +sigma is the constant zeroth-order proteasome-input term. Increased sigma models proteasome-availability therapy (Fig 2(I): sigma = 0.0416).

    kP     <- 0.0042;  label("Background degradation rate of free proteasome components (h^-1; INFERRED - see vignette Assumptions and deviations)")
    # Supinfo2 Eq (3): -kP*[P] is the linear background degradation term. NOT tabulated in Supinfo3 for any figure. Inferred here from the steady-state balance at F=CP~0: [P]_ss = sigma/kP; the paper's bifurcation diagrams (Fig 2) show healthy-state [P] ~ 1 (dimensionless), which requires kP = sigma = 0.0042. See vignette 'Assumptions and deviations' for the derivation. Under this inference the Fig 2(I) 'increased proteasome availability' scenario (sigma = 0.0416, kP unchanged) yields [P]_ss = 0.0416/0.0042 ~ 9.9 -- an ~10-fold rise consistent with the paper's Fig 2(I) description.

    gamma  <- 0.0042;  label("Rate constant of proteasome binding to larger alpha-synuclein aggregates (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' gamma = 0.0042 h^-1. Supinfo2 Eq (1)/(2)/(3): gamma*[F]*[P] is the bimolecular binding rate producing CP.

    nu     <- 0.0042;  label("Rate constant of proteolysis of the proteasome-aggregate complex (h^-1)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' nu = 0.0042 h^-1. Supinfo2 Eq (2)/(3): -nu*[CP] releases proteasome back to free pool while removing the aggregate.

    d_asyn <- 10;      label("Susceptibility of ROS production to misfolded alpha-synuclein (dimensionless)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' d_asyn = 10. Constant across all combined-model scenarios (matches Cloutier 2012 module). Multiplies the Hill activation in Eq (5) for the alpha_syn-driven ROS-production term.

    K_asyn <- 1;       label("Hill activation half-max concentration for alpha_syn-driven ROS production (dimensionless)")
    # Supinfo3 Table 'Figure 2 (A), Figure S1(C,D)' K_asyn = 1. Constant across all combined-model scenarios (matches Cloutier 2012 module). Enters Eq (5) as ([alpha_syn]/K_asyn)^4 / (1 + ([alpha_syn]/K_asyn)^4). With K_asyn = 1 the Hill function simplifies to alpha_syn^4/(1+alpha_syn^4) as printed in the paper; K_asyn is retained explicitly here so its value is source-traceable.

    # ------------------------------------------------------------------
    # Initial-condition parameters (dimensionless; scale factor 50 pM).
    # The paper simulates from a healthy attractor at low alpha_syn and
    # [P] ~ 1; F, CP, ROS all initialise near zero and equilibrate to
    # the healthy steady state within ~one simulated hour. Baseline
    # values placed here as free ini() parameters (rather than hardcoded
    # in model()) so a user can start from the diseased attractor by
    # over-riding bl_alpha_syn, bl_pfree, etc. See vignette 'Healthy
    # baseline' chunk.
    # ------------------------------------------------------------------
    bl_alpha_syn <- 0.1;   label("Initial misfolded alpha-synuclein monomer (dimensionless)")
    # Chosen to place the system near the low-alpha_syn healthy attractor at default parameters; equilibrates within simulated hours.

    bl_fagg      <- 0.1;   label("Initial alpha-synuclein larger-aggregate concentration F (dimensionless)")
    # As above; equilibrates via Eq (1) balance kagg*[alpha_syn] = gamma*[F]*[P].

    bl_pfree     <- 1.0;   label("Initial free proteasome component concentration P (dimensionless)")
    # Matches the paper's healthy [P]_ss ~ 1 (visible on Fig 2 bifurcation diagrams at low k0). Under the inferred kP = sigma this is exactly sigma/kP = 1.

    bl_pcomplex  <- 0.0;   label("Initial proteasome-aggregate complex CP (dimensionless)")
    # Equilibrates rapidly via Eq (2) balance gamma*[F]*[P] = nu*[CP].

    bl_ros       <- 0.0;   label("Initial cytoplasmic ROS (dimensionless)")
    # Equilibrates rapidly via Eq (5); healthy [ROS]_ss = (k1/k2)*k0 = k0 for k1 = k2.
  })

  model({
    # ------------------------------------------------------------------
    # Combined ROS + proteasome-sequestration ODE system.
    # Menon 2024 Supinfo2 Section 1.4, Equations (1)-(5).
    # State variables:
    #   alpha_syn : misfolded alpha-synuclein monomer [alpha*syn]
    #   fagg      : larger aggregates F (fibrils + oligomers, lumped)
    #   pfree     : free proteasome components [P]
    #   pcomplex  : proteasome-aggregate complex [CP]
    #   ros       : cytoplasmic ROS [ROS]
    # All concentrations dimensionless (scale factor 50 pM); time in hours.
    # ------------------------------------------------------------------

    # Hill activation: alpha_syn^4 / (K_asyn^4 + alpha_syn^4).
    # Reproduces the sharp switching signature that makes the ROS
    # positive feedback bistable. With K_asyn = 1 as printed in the
    # paper, this simplifies to alpha_syn^4/(1+alpha_syn^4).
    hill_asyn <- alpha_syn^4 / (K_asyn^4 + alpha_syn^4)

    # Eq (1): d[F]/dt = kagg*[alpha_syn] - gamma*[F]*[P]
    d/dt(fagg)     <- kagg * alpha_syn - gamma * fagg * pfree

    # Eq (2): d[CP]/dt = gamma*[F]*[P] - nu*[CP]
    d/dt(pcomplex) <- gamma * fagg * pfree - nu * pcomplex

    # Eq (3): d[P]/dt = sigma - kP*[P] - gamma*[F]*[P] + nu*[CP]
    d/dt(pfree)    <- sigma - kP * pfree - gamma * fagg * pfree + nu * pcomplex

    # Eq (4): d[alpha_syn]/dt = m + k3*[ROS] - kagg*[alpha_syn]
    #                           - k4*[P]*[alpha_syn] - krem*[alpha_syn]
    d/dt(alpha_syn) <- m + k3 * ros - kagg * alpha_syn -
                       k4 * pfree * alpha_syn - krem * alpha_syn

    # Eq (5): d[ROS]/dt = k1*(k0 + d_asyn * hill_asyn) - k2*[ROS]
    d/dt(ros)       <- k1 * (k0 + d_asyn * hill_asyn) - k2 * ros

    # Initial conditions from the bl_* ini() parameters.
    alpha_syn(0)  <- bl_alpha_syn
    fagg(0)       <- bl_fagg
    pfree(0)      <- bl_pfree
    pcomplex(0)   <- bl_pcomplex
    ros(0)        <- bl_ros
  })
}
