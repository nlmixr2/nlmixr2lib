DoldanMartelli_2013_EGF_IFN_chimera <- function() {
  description <- "In vitro (Daudi human Burkitt lymphoma cell line). Mechanistic kinetic model of an EGF-IFNalpha-2a chimeric ligand binding to EGFR and IFN receptor on the cell membrane: sequential two-subunit engagement, receptor lateral diffusion, and internalization (Doldan-Martelli 2013). Default parameters are wild-type IFN chimera in Daudi-EGFR cells (overexpressing EGFR ~300x parental); k2on / k2off can be overridden for K133A and R144A IFN mutants, and R1_0 / R2_0 for parental Daudi cells (see vignette)."
  reference <- "Doldan-Martelli V, Guantes R, Miguez DG. A mathematical model for the rational design of chimeric ligands in selective drug therapies. CPT Pharmacometrics Syst Pharmacol. 2013 Feb 13;2(2):e26. doi:10.1038/psp.2013.2."
  vignette <- "DoldanMartelli_2013_EGF_IFN_chimera"
  paper_specific_compartments <- c("egfr", "ifnr", "c_egf", "c_ifn", "c_full")

  units <- list(
    time = "min",
    dosing = "(none; constant extracellular ligand)",
    concentration = "molecules per cell"
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Daudi human Burkitt lymphoma cell line; parental Daudi and Daudi-EGFR engineered to overexpress EGFR ~300x)",
    n_subjects = NA_integer_,
    n_studies = 1,
    disease_state = "In-vitro cell-based assay; chimeric design targets cells overexpressing EGFR (proposed model for cancers and other diseases with EGFR upregulation).",
    notes = "Mechanistic deterministic ODE model -- no clinical subjects, no IIV, no residual error. Initial-condition receptor counts: parental Daudi cells have R1(0)=22 EGFR and R2(0)=2800 IFNR molecules per cell (Table 1, ref 14); Daudi-EGFR cells have R1(0)=5640 EGFR and R2(0)=3600 IFNR per cell. The Methods (paper p.6) treat extracellular ligand concentration L as constant during the 8-h simulated assay window because internalized ligand is a small fraction of bulk supply. Three IFN variants are reported in the paper (wild-type, K133A, R144A) with progressively lower IFNR affinity; this file defaults to wild-type and the vignette shows how to swap parameters for the mutants."
  )

  # All parameter values come from paper Table 1 (page 4); references in trailing
  # comments are the original-source citations as numbered in the paper.
  # Every parameter is held fixed -- this is a deterministic mechanistic model with
  # no estimation step; fixed() is used so downstream tooling (and a user who
  # later wraps the model into an estimation context) sees the provenance clearly.
  ini({
    # ----- EGF / EGFR system (paper Table 1, top half) -----
    k1on  <- fixed(0.09);    label("EGF-to-EGFR 3D association rate constant (1/(nmol/L)/min)")        # paper Table 1, ref 30
    k1off <- fixed(0.24);    label("EGF-EGFR dissociation rate constant (1/min)")                       # paper Table 1, ref 30
    ke1   <- fixed(0.15);    label("Internalization rate of EGF-bound complex C1 (1/min)")              # paper Table 1, ref 31
    h1    <- fixed(90);      label("EGFR extracellular-domain height above membrane (Angstrom)")        # paper Table 1, ref 32-33
    D1    <- fixed(2.2e-10); label("EGFR membrane lateral diffusion coefficient (cm^2/s)")              # paper Table 1, ref 32-33; midpoint of 2 - 2.4e-10

    # ----- IFN / IFNR system: WILD-TYPE defaults (paper Table 1, middle block) -----
    # Override k2on / k2off / ke2 to switch to K133A or R144A IFN mutants (see vignette).
    k2on  <- fixed(0.22);    label("IFN-to-IFNR 3D association rate constant, WT IFNalpha-2a (1/(nmol/L)/min); K133A=0.041, R144A=0.021")  # paper Table 1, ref 29
    k2off <- fixed(0.66);    label("IFN-IFNR dissociation rate constant, WT IFNalpha-2a (1/min); K133A=1.08, R144A=2.58")                  # paper Table 1, ref 29
    ke2   <- fixed(0.046);   label("Internalization rate of IFN-bound complex C2 (1/min)")              # paper Table 1, ref 28; same value for WT and mutants
    h2    <- fixed(50);      label("IFNR extracellular-domain height above membrane (Angstrom)")        # paper Table 1, ref 34
    D2    <- fixed(1e-10);   label("IFNR membrane lateral diffusion coefficient (cm^2/s)")              # paper Table 1, ref 21

    # ----- Cell-membrane geometry and chimera linker (paper Table 1, bottom block) -----
    A     <- fixed(900);     label("Typical mammalian-cell surface area (square micrometers)")          # paper Table 1, ref 23
    a     <- fixed(48.5);    label("Chimera linker effective end-to-end distance (Angstrom)")           # paper Table 1, ref 14; computed via Eq. 1 (worm-like chain, 7 Gly4-Ser subunits, residue length 3.8 A)

    # ----- Initial receptor counts: Daudi-EGFR defaults (paper Table 1, ref 14) -----
    # Override to R1_0=22, R2_0=2800 for parental Daudi cells (see vignette).
    R1_0  <- fixed(5640);    label("Initial EGFR molecules per cell, Daudi-EGFR default (parental Daudi=22)")     # paper Table 1, ref 14
    R2_0  <- fixed(3600);    label("Initial IFNR molecules per cell, Daudi-EGFR default (parental Daudi=2800)")  # paper Table 1, ref 14

    # ----- Extracellular ligand stimulation (held constant per paper Methods) -----
    L     <- fixed(1);       label("Extracellular chimeric-ligand concentration (nmol/L); held constant per paper Methods")  # paper Methods Eq. 11-15 explanation
  })

  model({
    # Alias ini() THETAs to model-local variables. This step is required because
    # nlmixr2est's mu-reference check forbids more than one THETA in a single
    # expression, and the ODE right-hand sides legitimately reference many
    # THETAs each. The aliased variables (kk1on, kk1off, ...) carry the same
    # values; structural-model arithmetic uses them directly below.
    kk1on  <- k1on
    kk1off <- k1off
    kke1   <- ke1
    hh1    <- h1
    DD1    <- D1
    kk2on  <- k2on
    kk2off <- k2off
    kke2   <- ke2
    hh2    <- h2
    DD2    <- D2
    AA     <- A
    aa     <- a
    RR1_0  <- R1_0
    RR2_0  <- R2_0
    LL     <- L

    # Physical constants
    Nav <- 6.022e23  # Avogadro's number (molecules / mol)

    # Effective reaction volumes V_i = A * (h_i + a) (paper Methods, text after Eq. 9)
    # A is in um^2, h_i and a are in A; 1 um^2 * 1 A = 1e-22 m^3 = 1e-19 L.
    V1_L <- AA * (hh1 + aa) * 1e-19
    V2_L <- AA * (hh2 + aa) * 1e-19

    # Sum of receptor diffusion coefficients, converted from cm^2/s to um^2/min.
    # 1 cm^2/s = 1e8 um^2/s = 6e9 um^2/min
    D_total <- (DD1 + DD2) * 6e9

    # Average half-distance between receptors on the cell surface (paper Eq. 8); um.
    # Computed from INITIAL receptor counts so b_i is a structural property (the
    # 2D rate constants k_diff_i, k'_on_i are constants of the model, not state-dependent).
    b1 <- sqrt(AA / (pi * RR1_0))
    b2 <- sqrt(AA / (pi * RR2_0))

    # Convert linker length a from Angstrom to micrometers for the b_i / a ratio.
    a_um <- aa * 1e-4

    # Diffusive rate constants (paper Eq. 7); paper's "# molecules^-1 min^-1" convention.
    kdiff1 <- 2 * pi * D_total / (AA * log(b1 / a_um))
    kdiff2 <- 2 * pi * D_total / (AA * log(b2 / a_um))

    # Effective 2D association-rate constants (paper Eq. 9). k_on is given in (nmol/L)^-1 min^-1
    # in Table 1; multiply by 1e9 to convert to (mol/L)^-1 min^-1 before dividing by N_av * V_i (L).
    k1on_prime <- kk1on * 1e9 / (Nav * V1_L)
    k2on_prime <- kk2on * 1e9 / (Nav * V2_L)

    # Coupling rate constants k_c (paper Eq. 6) and capture probabilities gamma (paper Eq. 10 text).
    kc1 <- 1 / (1 / kdiff1 + 1 / k1on_prime)
    kc2 <- 1 / (1 / kdiff2 + 1 / k2on_prime)
    gamma1 <- k1on_prime / (kdiff1 + k1on_prime)
    gamma2 <- k2on_prime / (kdiff2 + k2on_prime)

    # Uncoupling rates k_u (paper Eq. 10) and full-complex internalization rate (paper text).
    ku1 <- (1 - gamma1) * kk1off
    ku2 <- (1 - gamma2) * kk2off
    ke3 <- kke1 + kke2

    # ODE system (paper Eqs. 11-15). State variables (molecules per cell):
    #   egfr   = R1 (free EGF receptor)
    #   ifnr   = R2 (free IFN receptor)
    #   c_egf  = C1 (EGF subunit bound; IFN subunit free)
    #   c_ifn  = C2 (IFN subunit bound; EGF subunit free)
    #   c_full = C3 (chimera bound via both subunits to one EGFR + one IFNR)
    d/dt(egfr)   <-  kk1off * c_egf  + ku1 * c_full - kk1on * egfr * LL - kc1 * egfr * c_ifn
    d/dt(ifnr)   <-  kk2off * c_ifn  + ku2 * c_full - kk2on * ifnr * LL - kc2 * ifnr * c_egf
    d/dt(c_egf)  <-  kk1on * egfr * LL + ku2 * c_full - kk1off * c_egf  - kc2 * c_egf * ifnr - kke1 * c_egf
    d/dt(c_ifn)  <-  kk2on * ifnr * LL + ku1 * c_full - kk2off * c_ifn  - kc1 * c_ifn * egfr - kke2 * c_ifn
    d/dt(c_full) <-  kc1 * egfr * c_ifn + kc2 * ifnr * c_egf - (ku1 + ku2 + ke3) * c_full

    # Initial conditions: receptors at baseline, no complexes formed.
    egfr(0)   <- RR1_0
    ifnr(0)   <- RR2_0
    c_egf(0)  <- 0
    c_ifn(0)  <- 0
    c_full(0) <- 0

    # Reported readouts (paper Results, Figure 2 caption):
    # total IFN-active complexes = C2 + C3 = number of receptor-engaged IFN molecules,
    # which correlates with pSTAT1 signaling per ref. 15 cited in the paper.
    ifn_active <- c_ifn + c_full
  })
}
