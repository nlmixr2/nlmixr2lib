Kim_2024_carmrna_invitro_qsp <- function() {
  description <- paste(
    "In vitro (human CD8+ T cells).",
    "QSP (mechanism-based cell-level PK-PD, Monolix 2023R1) model of transient",
    "surface chimeric antigen receptor (CAR) expression on naive T cells after",
    "transfection with a CD8-targeted polymeric nanoparticle (NP) carrying",
    "in-vitro-transcribed 1928z CAR mRNA. An anti-CD8 antibody on the NP binds a",
    "constant pool of free CD8 receptors (kon, koff), the NP-receptor complex is",
    "internalised (kint) and dilutes with untransduced T-cell proliferation",
    "(kg_cell), the internalised NP releases free intracellular mRNA that",
    "degrades (kdeg_trans), and the free mRNA drives surface CAR protein",
    "formation through a sigmoidal Emax (Hill) process opposed by first-order",
    "protein degradation (kdeg). States are concentrations (pM); the output is",
    "percent surface CAR expression.",
    sep = " "
  )
  reference <- paste(
    "Kim SJ, Mugundu GM, Singh AP (2024).",
    "Translational PK-PD model for in vivo CAR-T-cell therapy delivered using",
    "CAR mRNA-loaded polymeric nanoparticle vector.",
    "Clin Transl Sci 17(12):e70101.",
    "doi:10.1111/cts.70101.",
    "In vitro cell-level model: Table 2 ('In vitro' column) and Appendix S1",
    "Eqs. 2-5; kdeg_trans estimated separately from the relative intracellular",
    "mRNA data with a one-compartment model (Appendix S1 Eq. 6).",
    sep = " "
  )
  vignette <- "Kim_2024_car_mrna_nanoparticle"
  units <- list(
    time = "h",
    dosing = "pM",
    concentration = "percent surface CAR expression"
  )

  # No patient/experimental covariates: the model is a typical-value mechanism
  # with IIV on the CAR-formation parameters only. Dose amount and timing are
  # supplied through the event table, not as covariates.
  covariateData <- list()

  # Every state is a concentration (pM) in the cell-culture system, except
  # carexpr which is the percent surface CAR expression (a PD readout). These
  # are paper-mechanistic states with no analogue in the canonical compartment
  # register; see paper_specific_compartments.
  paper_specific_compartments <- c("mrna_media", "rc_complex", "mrna_free", "carexpr")

  compartmentData <- list(
    mrna_media = list(
      analyte = "CAR mRNA (in anti-CD8 nanoparticle)",
      units = "pM",
      specimen = "administration site",
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
    species = "in vitro (human CD8+ T cells)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "not applicable (cell culture)",
    dose_range = paste(
      "Single transfection with 3 ug of 1928z CAR mRNA in a CD8-targeted",
      "polymeric nanoparticle, in a suspension of 1e6 CD8+ T cells/mL at an",
      "assumed media volume of 100 uL (nine replicates). 3 ug / MW 1e6 Da in",
      "100 uL = 30000 pM, the mRNA_media initial concentration."
    ),
    notes = paste(
      "Normalised qPCR CAR mRNA and flow-cytometry surface-CAR time courses",
      "(days 1, 2, 5, 8, 10 post-transfection) digitised by the original",
      "authors from Parayath et al. (the study Kim 2024 cites as reference 23).",
      "Binding and internalisation parameters (kon, koff, kint), the free",
      "receptor concentration (rtot) and the T-cell proliferation rate",
      "(kg_cell) were fixed to literature values; the mRNA degradation rate",
      "(kdeg_trans) was estimated separately from the relative intracellular",
      "mRNA data and then fixed; EC50, the Hill coefficient and the protein",
      "degradation rate (kdeg) were estimated on the surface-CAR data."
    )
  )

  ini({
    # ================================================================
    # Fixed mechanistic constants (Kim 2024 Table 2, 'In vitro' column).
    # Binding / internalisation / receptor / proliferation constants were
    # fixed to literature values (references 43-46); Emax was fixed to 1 to
    # stay within physiologically plausible values.
    # ================================================================
    kon <- fixed(0.00202) ; label("Association rate of anti-CD8 Ab with CD8 receptor (1/pM/h)") # Table 2, ref 44
    koff <- fixed(3.60) ; label("Dissociation rate of anti-CD8 Ab with CD8 receptor (1/h)")     # Table 2, ref 44
    kint <- fixed(0.0133) ; label("Internalisation rate of the NP-receptor complex (1/h)")      # Table 2, ref 43
    kg_cell <- fixed(0.0134) ; label("Proliferation (growth) rate of untransduced CD8+ T cells (1/h)") # Table 2, ref 45
    rtot <- fixed(2483) ; label("Concentration of free total CD8 receptors (pM)")               # Table 2, ref 46
    emax <- fixed(1) ; label("Maximum effect on CAR formation (unitless)")                      # Table 2
    # kdeg_trans was estimated (0.037, RSE 10.5%; IIV 30%, RSE 28.9%) from the
    # relative intracellular mRNA data with a one-compartment model (Appendix S1
    # Eq. 6), then fixed in the surface-CAR model; the 30% IIV belongs to that
    # separate mRNA-only fit, so it is not carried here.
    kdeg_trans <- fixed(0.037) ; label("Free intracellular mRNA degradation rate (1/h)")        # Table 2 / Appendix S1 Eq. 6

    # ================================================================
    # Estimated surface-CAR formation / degradation parameters
    # (Kim 2024 Table 2, 'In vitro' column). IIV is log-normal.
    # omega^2 = log(1 + CV^2) with CV the apparent IIV coefficient of variation.
    # ================================================================
    lec50 <- log(941) ; label("Free intracellular mRNA concentration for 50% of max CAR formation (pM)") # Table 2, RSE 1.73%
    lhill <- log(13.09) ; label("Hill coefficient of CAR formation (unitless)")                 # Table 2, 'n1', RSE 8.06%
    lkdeg <- log(1.45) ; label("Surface CAR protein degradation rate (1/h)")                    # Table 2, 'k_deg', RSE 1.99%

    etalec50 ~ 0.00148115 # Table 2, IIV on EC50 = 3.85 percent CV (RSE 30.1%); log(1 + 0.0385^2)
    etalhill ~ 0.0403046  # Table 2, IIV on n1 = 20.28 percent CV (RSE 30.8%); log(1 + 0.2028^2)
    etalkdeg ~ 0.000184943 # Table 2, IIV on k_deg = 1.36 percent CV (RSE 87.6%); log(1 + 0.0136^2)

    # Combined residual error on percent surface CAR expression (Table 2).
    addSd <- 0.041 ; label("Additive (constant) residual error on percent CAR (percent)")       # Table 2, 'a', RSE 11.7%
    propSd <- 0.18 ; label("Proportional residual error on percent CAR (fraction)")             # Table 2, 'b', RSE 14.2%
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual CAR-formation parameters.
    # ----------------------------------------------------------------------
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)
    kdeg <- exp(lkdeg + etalkdeg)

    # ----------------------------------------------------------------------
    # 2. Derived terms.
    #    Sigmoidal (Hill) CAR-formation rate driven by free intracellular
    #    mRNA, and the dilution flux from untransduced T-cell proliferation.
    # ----------------------------------------------------------------------
    carprod <- emax * mrna_free^hill / (ec50^hill + mrna_free^hill)
    dilflux <- kg_cell * rc_complex

    # ----------------------------------------------------------------------
    # 3. ODE system (Appendix S1 Eqs. 2-5). States are concentrations (pM);
    #    the free receptor concentration (rtot) is held constant. Dose the
    #    mRNA_media compartment (30000 pM for 3 ug in 100 uL) at time 0.
    # ----------------------------------------------------------------------
    d/dt(mrna_media) <- -kon * mrna_media * rtot + koff * rc_complex
    d/dt(rc_complex) <- kon * mrna_media * rtot - koff * rc_complex - kint * rc_complex - dilflux
    d/dt(mrna_free) <- kint * rc_complex - kdeg_trans * mrna_free
    d/dt(carexpr) <- carprod - kdeg * carexpr

    # ----------------------------------------------------------------------
    # 4. Observation and residual error. The percent surface CAR expression
    #    state (carexpr) is observed directly.
    # ----------------------------------------------------------------------
    carexpr ~ add(addSd) + prop(propSd)
  })
}
