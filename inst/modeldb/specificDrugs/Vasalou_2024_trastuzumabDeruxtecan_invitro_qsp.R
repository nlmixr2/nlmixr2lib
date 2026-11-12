Vasalou_2024_trastuzumabDeruxtecan_invitro_qsp <- function() {
  description <- "In vitro (NCI-N87 and NCI-H322 cell monolayers). QSP. Cellular disposition model of DXd delivered by trastuzumab deruxtecan (T-DXd), adapted from Singh and Shah: ADC in the culture medium binds cell-surface antigen, the complex is internalized and degraded in endosomes/lysosomes releasing DXd into the cytosol, DXd crosses the cell membrane in both directions at one rate, and the ADC deconjugates in the medium and on the cell surface. Media and per-cell intracellular DXd concentrations are the outputs. Defaults are T-DXd in NCI-N87 cells; the same structure with Table S1 overrides describes NCI-H322 cells, datopotamab deruxtecan, isotype-control ADCs and unconjugated DXd. Deterministic: no IIV or residual error reported."
  reference <- "Vasalou C, Proia TA, Kazlauskas L, Przybyla A, Sung M, Mamidi S, Maratea K, Griffin M, Sargeant R, Urosevic J, Rosenbaum AI, Yuan J, Aluri KC, Ramsden D, Hariparsad N, Jones RDO, Mettetal JT. Quantitative evaluation of trastuzumab deruxtecan pharmacokinetics and pharmacodynamics in mouse models of varying degrees of HER2 expression. CPT Pharmacometrics Syst Pharmacol. 2024;13(6):994. doi:10.1002/psp4.13133. In-vitro model Equations 2-9 and Table S1 of Appendix S1."
  vignette <- "Vasalou_2024_trastuzumabDeruxtecan"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM"
  )
  # ADC and unconjugated DXd are added to the culture medium, which is not
  # named central or depot, so buildModelDb() cannot infer it.
  dosing <- c("media_adc", "media_dxd")

  # In-vitro media / single-cell states (Appendix S1 Eqs 2-7) have no
  # canonical compartment.
  paper_specific_compartments <- c(
    "media_adc",
    "cell_bound_adc",
    "cell_endo_adc",
    "cell_dxd",
    "media_dxd",
    "dar"
  )

  compartmentData <- list(
    media_adc = list(analyte = "trastuzumab deruxtecan", units = "nmol", specimen = "not applicable", verified = TRUE),
    cell_bound_adc = list(
      analyte = "trastuzumab deruxtecan-HER2 complex",
      units = "molecules per cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    cell_endo_adc = list(
      analyte = "trastuzumab deruxtecan (endosomal/lysosomal)",
      units = "molecules per cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    cell_dxd = list(
      analyte = "DXd (released payload, intracellular)",
      units = "molecules per cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    media_dxd = list(analyte = "DXd (released payload)", units = "nmol", specimen = "not applicable", verified = TRUE),
    dar = list(analyte = "drug-to-antibody ratio", units = "(ratio)", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (NCI-N87 and NCI-H322 human cancer cell lines)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "NCI-N87 gastric carcinoma (HER2 high) and NCI-H322 non-small-cell lung carcinoma cell monolayers in 6-well plates",
    dose_range = "Continuous exposure to 10 ug/mL T-DXd, datopotamab deruxtecan, DAR4 or DAR8 isotype-control ADC for up to 48 h; 50 ng/mL unconjugated DXd continuous (N87) or 2 h followed by washout (N87)",
    notes = "1.5 mL medium per well; 1.2e6 (N87) or 2e6 (H322) cells; three replicate wells per condition and time point. Fitted in Phoenix 8.3.3. Cell doubling is omitted because of slow growth."
  )

  ini({
    # Table S1 values for T-DXd in NCI-N87 cells
    ncell <- fixed(1.2e6); label("Number of cells in culture, N87") # Table S1: N cell N87 1200000, Fixed
    mv <- fixed(0.0015); label("Medium volume MV (L)") # Table S1: 0.0015, Fixed
    vcell <- fixed(3e-12); label("Volume of a single cell (L)") # Table S2: Vcell 3*10^-12, 'Fixed [3]' (Singh and Shah, the in-vitro model source)
    ag_cell <- fixed(3.5e6); label("Antigen (HER2) receptors per cell, N87") # Table S1: Ag cell N87 3.5*10^6
    lkoff <- fixed(log(0.014)); label("ADC/antigen dissociation rate koff (1/h)") # Table S1: Koff = 0.014 (Maas), Fixed
    kd <- fixed(0.047); label("ADC/antigen equilibrium dissociation constant KD (nM)") # Table S1: KD 0.047 for T-DXd, Fixed; Kon = koff/KD = 0.3
    lkint <- log(0.015); label("Internalization rate of ADC/antigen complex kint (1/h)") # Table S1: N87 0.015 (15%CV), model fit
    lkdeg <- fixed(log(0.3)); label("Endosomal/lysosomal ADC degradation rate kdegrad (1/h)") # Table S1: Kdegrad 0.3, Fixed [3]
    lkdec <- log(5.9e-5); label("ADC deconjugation rate kdec, N87 (1/h)") # Table S1: Kdec N87 5.9e-5 (15%CV), model fit
    lkdiff <- log(0.7); label("DXd influx = efflux rate across the cell membrane kinD = keffD (1/h)") # Table S1: kinD/keffD 0.7 (16%CV), model fit
    dar0 <- fixed(8); label("Drug-to-antibody ratio at the start of incubation (unitless)") # Table S2: DAR0 8 for T-DXd
  })

  model({
    koff <- exp(lkoff)
    kon <- koff / kd
    kint <- exp(lkint)
    kdeg <- exp(lkdeg)
    kdec <- exp(lkdec)
    kdiff <- exp(lkdiff)

    # Scaling factor from molecules to nmol (Table S1: 10^9 / 6.023e23)
    sf <- 1e9 / 6.023e23

    # Medium ADC concentration (nM); the state holds the medium amount (nmol)
    cm_adc <- media_adc / mv

    dar(0) <- dar0

    # Eq 2 (multiplied through by MV so the state is an amount)
    d/dt(media_adc) <- (-kon * cm_adc * (ag_cell - cell_bound_adc) +
      koff * cell_bound_adc) * ncell * sf - kdec * media_adc
    # Eqs 3-5: molecules per cell
    d/dt(cell_bound_adc) <- kon * cm_adc * (ag_cell - cell_bound_adc) -
      (koff + kint + kdec) * cell_bound_adc
    d/dt(cell_endo_adc) <- kint * cell_bound_adc - kdeg * cell_endo_adc
    d/dt(cell_dxd) <- kdeg * cell_endo_adc * dar +
      kdiff * vcell / mv * media_dxd / sf - kdiff * cell_dxd
    # Eq 6: medium DXd amount (nmol)
    d/dt(media_dxd) <- kdec * cm_adc * dar * mv +
      (kdiff * cell_dxd + kdec * cell_bound_adc * dar) * ncell * sf -
      kdiff * vcell / mv * ncell * media_dxd
    # Eq 7
    d/dt(dar) <- -kdec * dar

    # Eqs 8-9: medium and intracellular DXd concentrations (nM)
    Cmedia_dxd <- media_dxd / mv
    Ccell_dxd <- cell_dxd * sf / vcell
  })
}
