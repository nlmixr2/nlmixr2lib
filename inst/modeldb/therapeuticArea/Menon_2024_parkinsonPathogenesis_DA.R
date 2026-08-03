Menon_2024_parkinsonPathogenesis_DA <- function() {
  description <- paste0(
    "QSP. Dopamine compartmentalization + ROS-feedback variant of the ",
    "Menon-Bakshi-Krishnan 2024 Parkinson's-disease pathogenesis model. ",
    "Four ODEs describe the intracellular time-evolution (hours) of ",
    "misfolded alpha-synuclein monomer (alpha_syn), cytoplasmic reactive ",
    "oxygen species (ros), vesicular dopamine (da_ves), and cytoplasmic ",
    "dopamine (da_cyto). This variant REPLACES the proteasomal-module of ",
    "the combined model (Menon_2024_parkinsonPathogenesis) with an explicit ",
    "three-compartment DA metabolism representation (vesicular / cytoplasmic ",
    "/ extracellular; extracellular DA fixed as a boundary condition D_ex). ",
    "Cytoplasmic DA contributes to the ROS-production term (parameter kDA), ",
    "which now enters the same Hill-activated ROS-positive-feedback that ",
    "drives alpha-synuclein misfolding. Misfolded alpha-synuclein aggregates ",
    "in turn drive DA leakage from vesicles into the cytoplasm (parameter ",
    "tr2 applied through the same Hill activation). The paper uses this ",
    "variant to demonstrate that inhibiting DA transport from extracellular ",
    "space into the cytoplasm (parameter tr3) raises the bistable-switching ",
    "threshold in alpha_syn -- a potential DA-transporter-targeting ",
    "therapeutic mechanism. Default ini() values are the Figure S8(B) ",
    "scenario (constant k_synth DA production, tr3 = 40000 h^-1 baseline ",
    "high-transport regime). The alternative Figure S8(C) scenario with ",
    "extracellular-DA-mediated feedback inhibition of DA synthesis ",
    "(k_synth = alpha/(beta + D_ex)) is documented in the vignette. ",
    "DA concentrations are in uM (paper: 'DA concentrations are in uM', ",
    "Supinfo2 Section 1.6); alpha_syn and ROS remain dimensionless (scale ",
    "factor 50 pM). This is a mechanistic ODE model without drug input, ",
    "without IIV, and without residual error. Companion combined-model ",
    "variant is packaged as Menon_2024_parkinsonPathogenesis."
  )
  reference <- paste(
    "Menon G, Bakshi S, Krishnan J. (2024).",
    "The interaction of core modules as a basis for elucidating network",
    "behavior determining Parkinson's disease pathogenesis.",
    "CPT Pharmacometrics Syst Pharmacol 13(3):335-340.",
    "doi:10.1002/psp4.13108.",
    "DA-compartmentalization structure inherited from Best, Nijhout & Reed",
    "(2009) Theor Biol Med Model 6:21 (see modellib('Menon_2024_parkinsonPathogenesis')",
    "for the companion combined ROS + proteasome-sequestration model).",
    sep = " "
  )
  vignette <- "Menon_2024_parkinsonPathogenesis"
  paper_specific_compartments <- c("alpha_syn", "ros", "da_ves", "da_cyto")

  units <- list(
    time          = "h",
    dosing        = "(none; mechanistic DA-compartmentalization variant, no drug input)",
    concentration = "(alpha_syn and ROS dimensionless with 50 pM scale factor per Supinfo2 Section 1.4; DA states da_ves, da_cyto in uM per Supinfo2 Section 1.6)"
  )

  covariateData <- list()

  population <- list(
    species        = "in silico (dopaminergic neuron cell-level mechanistic model)",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "(not applicable)",
    weight_range   = "(not applicable)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Parkinson's disease pathogenesis, DA-compartmentalization variant. The paper demonstrates that reducing DA transport from extracellular space into the cytoplasm (parameter tr3) lowers steady-state cytoplasmic DA and raises the bistable-switching threshold in alpha_syn -- a mechanism for targeting DA transporters therapeutically.",
    dose_range     = "(not applicable)",
    regions        = "(not applicable)",
    notes          = "Extracted from Menon 2024 Supinfo2 Section 1.6 (Equations 6-9) and Supinfo3 Figure S8(B) parameter table. The DA-compartmentalization variant does NOT include the proteasomal-sequestration module (the combined 5-ODE model is the companion Menon_2024_parkinsonPathogenesis). Extracellular DA (D_ex) is a boundary condition held constant per the paper's simplification ('For simplicity, in our analysis we assume that D_ex is maintained at a fixed level of the order of 1 nM'). The tr3 parameter's Supinfo3 entry '40, 20, 0.4 (x10^3)' encodes three scenarios in the Fig S8(B) sweep -- default here is the high-transport 40000 h^-1 baseline; the vignette also runs the low-transport 400 h^-1 scenario to show the switching-threshold effect. Not extracted here: the extracellular DA feedback inhibition variant (Fig S8(C), k_synth = alpha/(beta + D_ex)) is described in the vignette Errata for downstream users to construct if needed."
  )

  ini({
    # ------------------------------------------------------------------
    # Default parameter set: Figure S8(B) scenario (Supinfo3 Table
    # 'Figure S8 (B)'). Constant DA-synthesis rate k_synth = 20
    # uM/h; tr3 default = 40000 h^-1 (the highest of the three
    # scenario values '40, 20, 0.4 (x10^3)').
    #
    # k0 (basal cytoplasmic ROS production) is the bifurcation-input
    # parameter varied across the S8(B) sweep. Default 0.001 chosen
    # to hold the system in the healthy attractor at steady state.
    # ------------------------------------------------------------------

    k0     <- 0.001;   label("Basal cytoplasmic ROS production rate (dimensionless bifurcation input; see vignette)")
    # Supinfo2 Eq (7): k1*(k0 + d_asyn*hill_asyn + kDA*da_cyto) is the total ROS production term. k0 is the alpha_syn- and DA-independent basal component.

    m      <- 0;       label("Basal misfolded alpha-synuclein formation rate (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' m = 0 h^-1. The Fig S8 DA-focus scenarios use m = 0 so that alpha_syn is entirely driven by the ROS positive-feedback (no basal misfolding); this isolates the DA-transporter mechanism.

    k1     <- 35;      label("Cytoplasmic ROS input scaling rate constant (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k1 = 35 h^-1.

    k2     <- 35;      label("Cytoplasmic ROS removal rate constant (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k2 = 35 h^-1.

    k3     <- 0.0087;  label("Rate constant of alpha-synuclein misfolding driven by cytoplasmic ROS (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k3 = 0.0087 h^-1.

    k4     <- 0.0087;  label("Rate constant of misfolded alpha-synuclein removal (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k4 = 0.0087 h^-1. NOTE the DA variant Eq (6) removes alpha_syn by k4*[alpha_syn] (linear, no proteasome factor) -- the proteasomal machinery is not modelled in this variant.

    d_asyn <- 10;      label("Susceptibility of ROS production to misfolded alpha-synuclein (dimensionless)")
    # Supinfo3 Table 'Figure S8 (B)' d_asyn = 10.

    K_asyn <- 1;       label("Hill activation half-max concentration for alpha_syn (dimensionless)")
    # Supinfo3 Table 'Figure S8 (B)' K_asyn = 1.

    tr1    <- 40;      label("Rate constant of DA transport from cytoplasm into vesicles (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' tr1 = 40 h^-1. Supinfo2 Eq (8)/(9): tr1*[D_C] moves DA from cytoplasm to vesicle.

    tr2    <- 10;      label("Rate constant of alpha_syn-mediated DA leak from vesicles into cytoplasm (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' tr2 = 10 h^-1. Supinfo2 Eq (8)/(9): tr2 * hill_asyn * [D_V] moves DA from vesicle to cytoplasm; the Hill activation captures the alpha_syn-aggregate-mediated vesicle disruption.

    tr3    <- 40000;   label("Rate constant of DA transport from extracellular space into cytoplasm (h^-1; scenario high-transport)")
    # Supinfo3 Table 'Figure S8 (B)' tr3 = '40, 20, 0.4 (x10^3)' h^-1 -- three scenario values. Default 40000 (= 40 x 10^3) is the highest / baseline high-transport case; low-transport 400 h^-1 is the DA-transporter-inhibition target (see vignette). Supinfo2 Eq (9): tr3*[D_ex] is the extracellular-to-cytoplasm DA flux.

    D_ex   <- 0.002;   label("Extracellular DA concentration (uM; boundary condition, held constant)")
    # Supinfo3 Table 'Figure S8 (B)' D_ex = 0.002 uM. Supinfo2 Sec 1.6: 'For simplicity, in our analysis we assume that D_ex is maintained at a fixed level (of the order of 1 nM)'; 0.002 uM = 2 nM.

    k_synth <- 20;     label("Cytoplasmic DA synthesis rate (uM h^-1; constant for Fig S8(B))")
    # Supinfo3 Table 'Figure S8 (B)' k_synth = 20 uM h^-1. Fig S8(C) instead uses k_synth = alpha/(beta + D_ex) with alpha = 0.04 uM h^-1 and beta = 0.002 uM (extracellular-DA feedback inhibition of DA synthesis); see vignette Errata.

    k_deg   <- 40;     label("Cytoplasmic DA degradation rate (h^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k_deg = 40 h^-1.

    kDA     <- 0.05;   label("Rate constant of ROS production driven by cytoplasmic DA (h^-1 uM^-1)")
    # Supinfo3 Table 'Figure S8 (B)' k_DA = 0.05 h^-1 uM^-1. Supinfo2 Eq (7): k1*(... + kDA*[D_C]) adds cytoplasmic-DA-driven ROS to the total production term.

    # Initial conditions.
    bl_alpha_syn <- 0.05;  label("Initial misfolded alpha-synuclein monomer (dimensionless)")
    # Low initial value places the system near the low-alpha_syn healthy attractor.

    bl_ros       <- 0.0;   label("Initial cytoplasmic ROS (dimensionless)")

    bl_da_ves    <- 1.0;   label("Initial vesicular DA (uM)")
    # Order-of-magnitude healthy vesicular DA per Best 2009 estimates (Supinfo2 Sec 1.6: 'we choose parameter values consistent with the estimated reaction fluxes and steady state concentrations reported in [10]'); equilibrates rapidly.

    bl_da_cyto   <- 0.5;   label("Initial cytoplasmic DA (uM)")
    # Order-of-magnitude healthy cytoplasmic DA; equilibrates rapidly.
  })

  model({
    # ------------------------------------------------------------------
    # DA-compartmentalization + ROS-feedback ODE system.
    # Menon 2024 Supinfo2 Section 1.6, Equations (6)-(9).
    # State variables:
    #   alpha_syn : misfolded alpha-synuclein monomer [alpha*syn] (dimensionless)
    #   ros       : cytoplasmic ROS [ROS] (dimensionless)
    #   da_ves    : vesicular DA [D_V] (uM)
    #   da_cyto   : cytoplasmic DA [D_C] (uM)
    # Extracellular DA [D_ex] is a boundary condition (constant D_ex above).
    # ------------------------------------------------------------------

    # Hill activation: alpha_syn^4 / (K_asyn^4 + alpha_syn^4).
    hill_asyn <- alpha_syn^4 / (K_asyn^4 + alpha_syn^4)

    # Eq (6): d[alpha_syn]/dt = m + k3*[ROS] - k4*[alpha_syn]
    d/dt(alpha_syn) <- m + k3 * ros - k4 * alpha_syn

    # Eq (7): d[ROS]/dt = k1*(k0 + d_asyn*hill_asyn + kDA*[D_C]) - k2*[ROS]
    d/dt(ros)       <- k1 * (k0 + d_asyn * hill_asyn + kDA * da_cyto) - k2 * ros

    # Eq (8): d[D_V]/dt = tr1*[D_C] - tr2 * hill_asyn * [D_V]
    d/dt(da_ves)    <- tr1 * da_cyto - tr2 * hill_asyn * da_ves

    # Eq (9): d[D_C]/dt = k_synth + tr3*[D_ex] + tr2 * hill_asyn * [D_V]
    #                     - tr1*[D_C] - k_deg*[D_C]
    d/dt(da_cyto)   <- k_synth + tr3 * D_ex + tr2 * hill_asyn * da_ves -
                       tr1 * da_cyto - k_deg * da_cyto

    # Initial conditions.
    alpha_syn(0) <- bl_alpha_syn
    ros(0)       <- bl_ros
    da_ves(0)    <- bl_da_ves
    da_cyto(0)   <- bl_da_cyto
  })
}
