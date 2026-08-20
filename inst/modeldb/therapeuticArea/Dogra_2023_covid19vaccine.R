Dogra_2023_covid19vaccine <- function() {
  description <- "QSP. Mechanistic ODE model (17 states) of the adaptive immune response to COVID-19 mRNA vaccination in healthy and immunocompromised adults. States include naive and activated antigen-presenting cells (apc, apc_active); naive and effector CD4+ T-cells (cd4, cd4_effector); naive and effector CD8+ T-cells (cd8, cd8_effector); naive and activated B-cells (bcell, bcell_active); antibody-secreting plasma cells (plasma_cell); neutralizing antibody titer (antibody, U/mL); type-I and type-II interferons and interleukin-6 (ifn1, ifn2, il6, pg/mL); and a target-cell-limited SARS-CoV-2 infection subsystem (healthy_cell, infected_cell, virus). Antigen from mRNA vaccine doses enters via a nanoparticle depot (depot_antigen); the immune-status modulation factor f_immune in [0,1] scales the naive-cell homeostasis of CD4, CD8, and B cells to represent healthy (f_immune=1), mildly (~0.75), or severely (~0.55) immunocompromised subjects. Deterministic typical-value mechanism (no IIV, no residual error) calibrated jointly against Lucas 2020 infection kinetics (n=80) and Collier 2021 healthy-adult (n=31) and Peeters 2022 cancer-patient vaccine kinetics (n=63 chemotherapy, n=16 immunotherapy); validated against Bayart 2021 (n=158) and Papazisis 2022 (n=110). Parameter values are the healthy-adult vaccine calibration (Table 1 VCH column); IC (infection-calibration) alternates for the six calibration-dependent parameters are documented in comments."
  reference <- paste(
    "Dogra P, Schiavone C, Wang Z, Ruiz-Ramirez J, Caserta S, Staquicini DI,",
    "Markosian C, Wang J, Sostman HD, Pasqualini R, Arap W, Cristini V. (2023).",
    "A modeling-based approach to optimize COVID-19 vaccine dosing schedules",
    "for improved protection. JCI Insight 8(13):e169860.",
    "doi:10.1172/jci.insight.169860.",
    sep = " "
  )
  vignette <- "Dogra_2023_covid19vaccine"
  paper_specific_compartments <- c(
    "depot_antigen",
    "apc", "apc_active",
    "cd4", "cd4_effector",
    "cd8", "cd8_effector",
    "bcell", "bcell_active",
    "plasma_cell",
    "antibody",
    "ifn1", "ifn2", "il6",
    "healthy_cell", "infected_cell", "virus"
  )

  units <- list(
    time = "day",
    dosing = "dimensionless (unit dose)",
    concentration = "U/mL (neutralizing antibody titer)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot_antigen = list(analyte = "antigen", units = NA_character_, specimen = "administration site", verified = FALSE),
    apc           = list(analyte = "naive antigen-presenting cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    apc_active    = list(analyte = "activated antigen-presenting cells", units = NA_character_, specimen = "lymph", verified = FALSE),
    cd4           = list(analyte = "naive CD4+ T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    cd4_effector  = list(analyte = "effector CD4+ T-cells", units = NA_character_, specimen = "tissue", verified = FALSE),
    cd8           = list(analyte = "naive CD8+ T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    cd8_effector  = list(analyte = "effector CD8+ T-cells", units = NA_character_, specimen = "tissue", verified = FALSE),
    bcell         = list(analyte = "naive B-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    bcell_active  = list(analyte = "activated B-cells", units = NA_character_, specimen = "lymph", verified = FALSE),
    plasma_cell   = list(analyte = "antibody-secreting plasma cells", units = NA_character_, specimen = "plasma", verified = FALSE),
    antibody      = list(analyte = "neutralizing antibody", units = NA_character_, specimen = "serum", verified = FALSE),
    ifn1          = list(analyte = "type-I interferon", units = NA_character_, specimen = "blood cell", verified = FALSE),
    ifn2          = list(analyte = "type-II interferon", units = NA_character_, specimen = "tissue", verified = FALSE),
    il6           = list(analyte = "interleukin-6", units = NA_character_, specimen = "serum", verified = FALSE),
    healthy_cell  = list(analyte = "healthy cells", units = NA_character_, specimen = "tissue", verified = FALSE),
    infected_cell = list(analyte = "infected cells", units = NA_character_, specimen = "tissue", verified = FALSE),
    virus         = list(analyte = "SARS-CoV-2 virus", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 458L,
    n_studies = 5L,
    age_range = "adult (individual studies span 18-89 years)",
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Healthy adults (Collier 2021, Bayart 2021, Papazisis 2022); moderately infected COVID-19 patients (Lucas 2020 calibration set); cancer patients undergoing chemotherapy (Peeters 2022 n=63) or immunotherapy (Peeters 2022 n=16)",
    dose_range = "Pfizer-BioNTech BNT162b2 or Moderna mRNA-1273 COVID-19 mRNA vaccines, intramuscular; canonical two-dose primary series 21-28 days apart, boosters 5-9 months after second dose. Dose encoded as dimensionless amt=1 in event tables.",
    regions = "United States (Lucas 2020 Yale, Collier 2021 BIDMC, Papazisis 2022 Greece pooled), Belgium (Peeters 2022 UZA), multi-site (Bayart 2021)",
    notes = "Model was calibrated jointly against five clinical data sets: (1) Lucas et al. 2020 Nature 584:463-469 -- viral load and immune kinetics in 80 moderately infected COVID-19 patients (infection calibration IC); (2) Collier et al. 2021 Nature 596:417-422 -- neutralizing antibody and T-cell kinetics in 31 healthy adults receiving 2 doses of Pfizer-BioNTech (vaccine calibration healthy VCH); (3) Peeters et al. 2022 ESMO Open 6:100274 -- antibody kinetics in 63 chemotherapy and 16 immunotherapy cancer patients (VCC / VCI f-parameter calibration); (4) Bayart et al. 2021 Vaccines 9:1092 -- 158 healthy adults, Pfizer-BioNTech 2 doses (external validation); (5) Papazisis et al. 2022 Emerg Microbes Infect 11:2016-2023 -- 110 healthy adults, Pfizer-BioNTech 3 doses (external validation). Total n=458 across calibration and validation cohorts."
  )

  ini({
    # ============================================================================
    # All values are the healthy-adult vaccine calibration (Table 1 VCH column).
    # Six parameters differ between the vaccine (VCH) and infection (IC) calibrations;
    # the IC value is recorded in a comment alongside the VCH value for reference.
    # This model is deterministic (no IIV, no residual error) -- values are typical
    # mechanistic constants (endogenous / QSP model conventions).
    # ============================================================================

    # ---- Vaccine (antigen) parameters -------------------------------------------
    t_np       <- 1;        label("Characteristic time of vaccine-nanoparticle clearance (day)")                                                                                     # Dogra 2023 Table 1 T_NP row (ref 42; NP diameter 100 nm)
    k_eff      <- 18.95;    label("Michaelis constant for vaccine efficacy vs antibody titer (U/mL)")                                                                                # Dogra 2023 Table 1 K_eff row (estimated; Methods Eq 1)
    v_eff_max  <- 92.47;    label("Maximum theoretical vaccine efficacy (%)")                                                                                                        # Dogra 2023 Table 1 V_eff_max row (estimated; Methods Eq 1)
    ab_escape  <- 1;        label("Antibody-strain binding score (1 = wild-type, 0.2 = Omicron)")                                                                                    # Dogra 2023 Table 1 Ab_escape row (ref 61; user overrides for VOCs)

    # ---- Infection-related parameters ------------------------------------------
    # Drive the healthy_cell / infected_cell / virus subsystem. For a pure vaccine
    # simulation, healthy_cell(0) = infected_cell(0) = virus(0) = 0 and this
    # subsystem stays quiescent.
    beta_inf   <- 0.004;    label("Viral infectivity rate of healthy respiratory cells (mL/GE/day)")                                                                                 # Dogra 2023 Table 1 beta row (estimated)
    delta_ic   <- 0.15;     label("Cytopathic death rate of infected cells (1/day)")                                                                                                 # Dogra 2023 Table 1 delta row (ref 57)
    delta_c    <- 4.51e-5;  label("Infected-cell death rate mediated by effector CD8+ T-cells (mL/cell/day)")                                                                        # Dogra 2023 Table 1 delta_C row (estimated)
    p_v        <- 3.39;     label("Virion production rate by infected cells (GE/cell/day)")                                                                                          # Dogra 2023 Table 1 P_v row (estimated)
    k_ifn1     <- 4.86;     label("Michaelis constant for type-I IFN suppression of virion production (pg/mL)")                                                                      # Dogra 2023 Table 1 K_IFN1 row (estimated)
    k_apc      <- 1.16;     label("Rate of virus clearance by naive APCs (mL/cell/day)")                                                                                             # Dogra 2023 Table 1 k_APC row (estimated)
    k_ab       <- 0.11;     label("Antibody-mediated virus neutralization rate (mL/U/day)")                                                                                          # Dogra 2023 Table 1 k_Ab row (estimated)

    # ---- Innate immunity parameters (naive and activated APCs) -----------------
    gamma_apc  <- 0.4873;   label("Growth (regeneration) rate of naive APCs (1/day)")                                                                                                # Dogra 2023 Table 1 gamma_APC row (ref 66)
    t_apc      <- 0.159;    label("APC activation rate per antigen (1/day); VCH -- IC calibration reports 36.3")                                                                     # Dogra 2023 Table 1 T_APC row VCH column (estimated); IC value 36.3 documented but VCH used here
    k_v        <- 0.50;     label("Michaelis constant for antigen-induced APC activation (1/mL); VCH -- IC calibration reports 0.0625")                                              # Dogra 2023 Table 1 K_v row VCH column (estimated); IC value 0.0625 documented but VCH used here
    k_ifn2     <- 0.0835;   label("Michaelis constant for type-II IFN amplification of T-cell activation (pg/mL)")                                                                   # Dogra 2023 Table 1 K_IFN2 row (estimated)
    apc_bar    <- 1e6;      label("Carrying capacity of naive APCs (cell/mL)")                                                                                                       # Dogra 2023 Table 1 APC_bar row (ref 57; 10^6)
    delta_apc  <- 0.2;      label("Death rate of activated APCs (1/day)")                                                                                                            # Dogra 2023 Table 1 delta_APC row (ref 57)

    # ---- Cellular immunity parameters (CD4 and CD8 T-cells) --------------------
    gamma_cd4  <- 1.5122;   label("Growth rate of naive CD4+ T-cells (1/day)")                                                                                                       # Dogra 2023 Table 1 gamma_CD4 row (ref 67)
    cd4_bar    <- 10^5.8;   label("Carrying capacity of naive CD4+ T-cells (cell/mL); 10^5.8 ~ 6.31e5")                                                                               # Dogra 2023 Table 1 CD4_bar row (ref 57; 10^5.8)
    k_il6      <- 18.93;    label("Michaelis constant for IL-6-induced naive T-cell exhaustion (pg/mL)")                                                                             # Dogra 2023 Table 1 K_IL6 row (estimated)
    t_cd4      <- 0.0223;   label("Activation rate of naive CD4+ T-cells (mL/cell/day)")                                                                                             # Dogra 2023 Table 1 T_CD4 row (estimated)
    delta_t    <- 0.00039;  label("Death rate of effector CD4+ and CD8+ T-cells (1/day); VCH -- IC calibration reports 0.0075")                                                      # Dogra 2023 Table 1 delta_T row VCH column (estimated); IC value 0.0075 documented but VCH used here
    gamma_cd8  <- 2.0794;   label("Growth rate of naive CD8+ T-cells (1/day)")                                                                                                       # Dogra 2023 Table 1 gamma_CD8 row (ref 67)
    cd8_bar    <- 1e5;      label("Carrying capacity of naive CD8+ T-cells (cell/mL)")                                                                                               # Dogra 2023 Table 1 CD8_bar row (ref 57; 10^5)
    t_cd8      <- 0.023;    label("Activation rate of naive CD8+ T-cells (mL/cell/day)")                                                                                             # Dogra 2023 Table 1 T_CD8 row (estimated)

    # ---- Humoral immunity parameters (B-cells, plasma cells, antibody) --------
    gamma_b    <- 0.462;    label("Growth rate of naive B-cells (1/day)")                                                                                                            # Dogra 2023 Table 1 gamma_B row (ref 68)
    t_b        <- 0.4965;   label("Activation rate of naive B-cells (mL/cell/day)")                                                                                                  # Dogra 2023 Table 1 T_B row (estimated)
    bcell_bar  <- 1e5;      label("Carrying capacity of naive B-cells (cell/mL)")                                                                                                    # Dogra 2023 Table 1 B_bar row (ref 57; 10^5)
    t_bc       <- 0.36;     label("Differentiation rate of activated B-cells into plasma cells (mL/cell/day)")                                                                       # Dogra 2023 Table 1 T_BC row (ref 57)
    delta_p    <- 1.23;     label("Death rate of plasma cells (1/day); VCH -- IC calibration reports 0.0083")                                                                        # Dogra 2023 Table 1 delta_P row VCH column (ref 57 for IC; estimated for VCH); VCH used here
    p_ab       <- 0.35;     label("Antibody production rate per plasma cell (U/cell/day); VCH -- IC calibration reports 0.167")                                                      # Dogra 2023 Table 1 P_Ab row VCH column (estimated); IC value 0.167 documented but VCH used here
    cl_ab      <- 0.0027;   label("Antibody first-order clearance rate (1/day); VCH -- IC calibration reports 0.254")                                                                # Dogra 2023 Table 1 Cl_Ab row VCH column (estimated); IC value 0.254 documented but VCH used here

    # ---- Immunity mediator (cytokine) parameters -------------------------------
    p_ifn1     <- 4.2;      label("Production rate of type-I interferon per infected or vaccine-affected cell (pg/cell/day)")                                                        # Dogra 2023 Table 1 P_IFN1 row (estimated)
    delta_cyt  <- 1.71;     label("First-order degradation rate of cytokines IFN1, IFN2, IL-6 (1/day)")                                                                              # Dogra 2023 Table 1 delta_cyt row (estimated)
    p_ifn2     <- 0.174;    label("Production rate of type-II interferon per effector CD4+ or CD8+ T-cell (pg/cell/day)")                                                            # Dogra 2023 Table 1 P_IFN2 row (estimated)
    p_il6      <- 0.273;    label("Production rate of IL-6 per effector CD4+ / CD8+ T-cell or activated APC (pg/cell/day)")                                                          # Dogra 2023 Table 1 P_IL6 row (estimated)

    # ---- Immune-status modulation factor (patient-specific parameter) ---------
    # Paper Table 1 patient-specific parameter f: dimensionless factor in [0,1]
    # scaling naive CD4, CD8, and B-cell carrying capacities and initial conditions.
    # f=1 healthy adult (VCH); 0.5172 chemotherapy (VCC); 0.5885 immunotherapy (VCI);
    # virtual cohorts use f=0.9-1.0 healthy, 0.7-0.9 mildly immunocompromised,
    # 0.5-0.7 highly immunocompromised. Encoded here as a paper-specific
    # mechanistic parameter with the healthy default; users override at simulation
    # time via rxode2::rxSolve(mod, params = c(f_immune = 0.55), ...).
    f_immune   <- 1;        label("Immune-status modulation factor (0-1); 1 = healthy adult (VCH default)")                                                                          # Dogra 2023 Table 1 f row VCH (1); VCC 0.5172; VCI 0.5885
  })

  model({
    # ============================================================================
    # 17-state mechanistic QSP model of adaptive immune response to COVID-19 mRNA
    # vaccination (Dogra 2023 Supplemental Methods equations S1-S17).
    # ============================================================================

    # ---- Antigen kinetics (Ca(t), Eq S1) --------------------------------------
    # Paper Eq S1 defines Ca(t) as a sum of Gaussians centred at each vaccination
    # time tau_i with characteristic width t_np. For rxode2 / nlmixr2 event-table
    # simulation this is implemented as a first-order-decay depot compartment
    # (impulse response exp(-t/t_np), identical characteristic time t_np). Both
    # peak at Ca(tau) = Dose and share the ~1-day pulse width; the Gaussian and
    # exponential decays differ only in the sub-day tail shape and are
    # pharmacologically indistinguishable on the multi-week / multi-month
    # vaccination time scale relevant to the paper's optimisation analysis. See
    # vignette Errata for detail.
    d/dt(depot_antigen) <- -depot_antigen / t_np
    ag <- depot_antigen                          # Ag(t) in the vaccine case (paper: Ag=Ca(t))

    # ---- Naive and activated APCs (Eq S2, S3) ---------------------------------
    activation_apc <- t_apc * apc * ag / (k_v + ag)
    d/dt(apc)        <- gamma_apc * apc * (1 - apc / apc_bar) - activation_apc
    d/dt(apc_active) <- activation_apc - delta_apc * apc_active

    # ---- Naive and effector CD4+ T-cells (Eq S4, S5) --------------------------
    il6_exhaustion <- 1 - il6 / (k_il6 + il6)
    ifn2_amp       <- 1 + ifn2 / (k_ifn2 + ifn2)
    cd4_capacity   <- f_immune * cd4_bar * il6_exhaustion
    activation_cd4 <- t_cd4 * apc_active * cd4 * ifn2_amp
    d/dt(cd4)          <- gamma_cd4 * cd4 * (1 - cd4 / cd4_capacity) - activation_cd4
    d/dt(cd4_effector) <- activation_cd4 - delta_t * cd4_effector

    # ---- Naive and effector CD8+ T-cells (Eq S6, S7) --------------------------
    cd8_capacity   <- f_immune * cd8_bar * il6_exhaustion
    activation_cd8 <- t_cd8 * apc_active * cd8 * ifn2_amp
    d/dt(cd8)          <- gamma_cd8 * cd8 * (1 - cd8 / cd8_capacity) - activation_cd8
    d/dt(cd8_effector) <- activation_cd8 - delta_t * cd8_effector

    # ---- Naive and activated B-cells, plasma cells, antibody (Eq S8-S11) ------
    activation_b      <- t_b * apc_active * bcell
    differentiation_b <- t_bc * cd4_effector * bcell_active
    d/dt(bcell)         <- gamma_b * bcell * (1 - bcell / (f_immune * bcell_bar)) - activation_b
    d/dt(bcell_active)  <- activation_b - differentiation_b
    d/dt(plasma_cell)   <- differentiation_b - delta_p * plasma_cell
    d/dt(antibody)      <- p_ab * plasma_cell - cl_ab * antibody

    # ---- Cytokines (Eq S12-S14) -----------------------------------------------
    d/dt(ifn1) <- p_ifn1 * (infected_cell + depot_antigen) - delta_cyt * ifn1
    d/dt(ifn2) <- p_ifn2 * (cd4_effector + cd8_effector) - delta_cyt * ifn2
    d/dt(il6)  <- p_il6  * (cd4_effector + cd8_effector + apc_active) - delta_cyt * il6

    # ---- SARS-CoV-2 infection subsystem (Eq S15-S17) --------------------------
    # For a pure vaccine simulation set healthy_cell(0) = infected_cell(0) =
    # virus(0) = 0 in the event table and this subsystem stays quiescent;
    # setting virus(0) > 0 activates infection dynamics.
    infection_flux <- beta_inf * healthy_cell * virus
    d/dt(healthy_cell)  <- -infection_flux
    d/dt(infected_cell) <- infection_flux - delta_ic * infected_cell - delta_c * infected_cell * cd8_effector
    d/dt(virus)         <- p_v * infected_cell * (1 - ifn1 / (k_ifn1 + ifn1)) - k_ab * virus * antibody - k_apc * virus * apc

    # ---- Initial conditions (Eq S2-S17) ---------------------------------------
    # Naive lymphocyte pools scaled by f_immune (immunosuppression). Infection
    # subsystem starts quiescent; user overrides healthy_cell(0) / virus(0) for
    # an infection simulation.
    depot_antigen(0)  <- 0
    apc(0)            <- apc_bar
    apc_active(0)     <- 0
    cd4(0)            <- f_immune * cd4_bar
    cd4_effector(0)   <- 0
    cd8(0)            <- f_immune * cd8_bar
    cd8_effector(0)   <- 0
    bcell(0)          <- f_immune * bcell_bar
    bcell_active(0)   <- 0
    plasma_cell(0)    <- 0
    antibody(0)       <- 0
    ifn1(0)           <- 0
    ifn2(0)           <- 0
    il6(0)            <- 0
    healthy_cell(0)   <- 0
    infected_cell(0)  <- 0
    virus(0)          <- 0

    # ---- Vaccine efficacy outputs (Eqs 1 and 2 in main paper) -----------------
    # V_eff(t) is an algebraic transform of antibody titer against a chosen
    # variant. Wild-type uses ab_escape=1; Omicron uses ab_escape=0.2 -- the
    # paper's Eq 2 penalises the effective antibody concentration by the
    # strain-specific binding score. The (antibody * ab_escape) form is confirmed
    # by the paper's own protection-threshold arithmetic: 154 U/mL * 1 = 770 U/mL
    # * 0.2 both yield 82.3% vaccine efficacy on the Michaelis-Menten curve
    # (Methods paragraph "for the VOCs, the protective threshold was corrected
    # for by using the binding score"). See vignette Errata for the OCR
    # ambiguity in the trimmed supplement.
    ab_effective <- antibody * ab_escape
    vaccine_efficacy <- v_eff_max * ab_effective / (k_eff + ab_effective)
  })
}
