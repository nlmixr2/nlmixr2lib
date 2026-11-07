Kim_2024_carmrna_mouse_qsp <- function() {
  description <- paste(
    "In vivo (mouse, immunodeficient NSG).",
    "QSP (mechanism-based cell-level PK-PD, Monolix 2023R1) model of transient",
    "surface chimeric antigen receptor (CAR) expression on circulating CD8+ T",
    "cells after intravenous dosing of a CD8-targeted polymeric nanoparticle",
    "(NP) carrying in-vitro-transcribed 1928z CAR mRNA. A two-compartment linear",
    "NP disposition (knp, k12, k21; parameters fixed from a separate NP PK fit)",
    "delivers mRNA to the central blood compartment, where an anti-CD8 antibody",
    "on the NP binds a constant pool of free CD8 receptors (kon, koff); the",
    "complex is internalised (kint) and dilutes with untransduced T-cell",
    "proliferation (kg_cell); the internalised NP releases free intracellular",
    "mRNA that degrades (kdeg_trans); and the free mRNA drives surface CAR",
    "protein formation through a sigmoidal Emax (Hill) process opposed by",
    "first-order protein degradation (kdeg). States are concentrations (pM); the",
    "output is percent surface CAR expression on CD8+ T cells.",
    sep = " "
  )
  reference <- paste(
    "Kim SJ, Mugundu GM, Singh AP (2024).",
    "Translational PK-PD model for in vivo CAR-T-cell therapy delivered using",
    "CAR mRNA-loaded polymeric nanoparticle vector.",
    "Clin Transl Sci 17(12):e70101.",
    "doi:10.1111/cts.70101.",
    "In vivo cell-level model: Table 2 ('In vivo' column), Appendix S1",
    "Eqs. 13-17, and the simulation code in Appendix S2. The NP PK parameters",
    "(knp, k12, k21, V1) were fixed from a two-compartment fit of a separate",
    "nanoparticle's published disposition (Appendix S1 Table S1, Eqs. 7-12).",
    sep = " "
  )
  vignette <- "Kim_2024_car_mrna_nanoparticle"
  units <- list(
    time = "h",
    dosing = "pM",
    concentration = "percent surface CAR expression"
  )

  # No patient covariates: typical-value mechanism with IIV on the CAR-formation
  # parameters only. Dose amount and timing are supplied through the event table.
  covariateData <- list()

  # Every state is a concentration (pM), except carexpr which is the percent
  # surface CAR expression (a PD readout). The NP central/peripheral mRNA
  # concentrations and the binding cascade are paper-mechanistic states with no
  # analogue in the canonical compartment register.
  paper_specific_compartments <- c(
    "mrna_central",
    "mrna_peripheral",
    "rc_complex",
    "mrna_free",
    "carexpr"
  )

  compartmentData <- list(
    mrna_central = list(
      analyte = "CAR mRNA (in anti-CD8 nanoparticle)",
      units = "pM",
      specimen = "whole blood",
      verified = TRUE
    ),
    mrna_peripheral = list(
      analyte = "CAR mRNA (in anti-CD8 nanoparticle)",
      units = "pM",
      specimen = "tissue",
      verified = TRUE
    ),
    rc_complex = list(
      analyte = "nanoparticle-CD8 receptor complex",
      units = "pM",
      specimen = "blood cell",
      verified = TRUE
    ),
    mrna_free = list(analyte = "free intracellular CAR mRNA", units = "pM", specimen = "blood cell", verified = TRUE),
    carexpr = list(
      analyte = "surface 1928z CAR protein (percent of CD8+ T cells)",
      units = "percent",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vivo mouse (immunodeficient NSG, human T-cell reconstituted)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = paste(
      "NSG mice inoculated with 1e6 CD19+ Raji cells (firefly luciferase),",
      "reconstituted 5 days later with 1e7 CD3+ human T cells"
    ),
    dose_range = paste(
      "Six weekly IV infusions of a CD8-targeted nanoparticle loaded with",
      "50 ug of 1928z CAR mRNA, starting 2 days after T-cell reconstitution.",
      "50 ug / MW 1e6 Da in the fitted NP central volume V1 = 7.87 mL =",
      "6353 pM, the per-dose mRNA_central increment."
    ),
    notes = paste(
      "Surface-CAR time course (percent 1928z CAR+ among CD8+ T cells in",
      "peripheral blood, sampled 2 days after each dose) digitised by the",
      "original authors from Parayath et al. (the study Kim 2024 cites as",
      "reference 23). The observed mean at the original regimen was",
      "10.7 +/- 3.6 percent CAR+. NP PK parameters were fixed from a separate",
      "two-compartment fit of a physicochemically similar nanoparticle; EC50,",
      "the Hill coefficient and the protein degradation rate were estimated on",
      "the in vivo surface-CAR data."
    )
  )

  ini({
    # ================================================================
    # Fixed NP disposition constants (Kim 2024 Table 2, 'In vivo' column;
    # derived in Appendix S1 Eqs. 7-12 from Table S1 and then fixed).
    # V1 = 7.87 mL (RSE 2.07%) is the fitted NP central volume. It does not
    # appear in the ODE: the states are concentrations, and V1 enters only in
    # converting the 50 ug dose to its 6353 pM central increment and in
    # deriving rtot. It is therefore documented here rather than carried as an
    # unused ini() parameter.
    # ================================================================
    lknp <- fixed(log(0.086)) ; label("NP elimination rate from the central compartment (1/h)")  # Table 2, 'k'
    lk12 <- fixed(log(0.64)) ; label("NP central-to-peripheral distribution rate (1/h)")         # Table 2
    lk21 <- fixed(log(0.61)) ; label("NP peripheral-to-central distribution rate (1/h)")         # Table 2

    # ================================================================
    # Fixed mechanistic constants (Kim 2024 Table 2, 'In vivo' column).
    # Binding / internalisation / receptor / proliferation constants fixed to
    # literature values (references 43-46); kdeg_trans fixed from the in vitro
    # relative-mRNA fit; Emax fixed to 1.
    # ================================================================
    kon <- fixed(0.00202) ; label("Association rate of anti-CD8 Ab with CD8 receptor (1/pM/h)") # Table 2, ref 44
    koff <- fixed(3.60) ; label("Dissociation rate of anti-CD8 Ab with CD8 receptor (1/h)")     # Table 2, ref 44
    kint <- fixed(0.0133) ; label("Internalisation rate of the NP-receptor complex (1/h)")      # Table 2, ref 43
    kg_cell <- fixed(0.0134) ; label("Proliferation (growth) rate of untransduced CD8+ T cells (1/h)") # Table 2, ref 45
    rtot <- fixed(2483) ; label("Concentration of free total CD8 receptors (pM)")               # Table 2, ref 46
    emax <- fixed(1) ; label("Maximum effect on CAR formation (unitless)")                      # Table 2
    kdeg_trans <- fixed(0.037) ; label("Free intracellular mRNA degradation rate (1/h)")        # Table 2, fixed from in vitro

    # ================================================================
    # Estimated surface-CAR formation / degradation parameters
    # (Kim 2024 Table 2, 'In vivo' column). IIV is log-normal.
    # omega^2 = log(1 + CV^2) with CV the apparent IIV coefficient of variation.
    # ================================================================
    lec50 <- log(124) ; label("Free intracellular mRNA concentration for 50% of max CAR formation (pM)") # Table 2, RSE 29.6%
    lhill <- log(1.71) ; label("Hill coefficient of CAR formation (unitless)")                 # Table 2, 'n1', RSE 11.2%
    lkdeg <- log(0.078) ; label("Surface CAR protein degradation rate (1/h)")                  # Table 2, 'k_deg', RSE 9.22%

    etalec50 ~ 0.0160003  # Table 2, IIV on EC50 = 12.7 percent CV (RSE 91.2%); log(1 + 0.127^2)
    etalhill ~ 0.00678319 # Table 2, IIV on n1 = 8.25 percent CV (RSE 79.9%); log(1 + 0.0825^2)
    etalkdeg ~ 0.00115533 # Table 2, IIV on k_deg = 3.40 percent CV (RSE 93.5%); log(1 + 0.034^2)

    # Proportional residual error on percent surface CAR expression (Table 2).
    propSd <- 0.395 ; label("Proportional residual error on percent CAR (fraction)")           # Table 2, 'b', RSE 8.55%
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual CAR-formation parameters.
    # ----------------------------------------------------------------------
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)
    kdeg <- exp(lkdeg + etalkdeg)

    # NP disposition micro-constants (fixed).
    knp <- exp(lknp)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # ----------------------------------------------------------------------
    # 2. Derived terms: sigmoidal (Hill) CAR-formation rate driven by free
    #    intracellular mRNA, and the proliferation dilution flux.
    # ----------------------------------------------------------------------
    carprod <- emax * mrna_free^hill / (ec50^hill + mrna_free^hill)
    dilflux <- kg_cell * rc_complex

    # ----------------------------------------------------------------------
    # 3. ODE system (Appendix S1 Eqs. 13-17; Appendix S2 code). States are
    #    concentrations (pM); the free receptor concentration (rtot) is held
    #    constant. Dose the mrna_central compartment (6353 pM for 50 ug) at
    #    each administration.
    # ----------------------------------------------------------------------
    d/dt(mrna_central) <- -knp * mrna_central - k12 * mrna_central +
      k21 * mrna_peripheral - kon * mrna_central * rtot + koff * rc_complex
    d/dt(mrna_peripheral) <- k12 * mrna_central - k21 * mrna_peripheral
    d/dt(rc_complex) <- kon * mrna_central * rtot - koff * rc_complex - kint * rc_complex - dilflux
    d/dt(mrna_free) <- kint * rc_complex - kdeg_trans * mrna_free
    d/dt(carexpr) <- carprod - kdeg * carexpr

    # ----------------------------------------------------------------------
    # 4. Observation and residual error. The percent surface CAR expression
    #    state (carexpr) is observed directly.
    # ----------------------------------------------------------------------
    carexpr ~ prop(propSd)
  })
}
