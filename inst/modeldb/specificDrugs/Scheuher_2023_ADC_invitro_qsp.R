Scheuher_2023_ADC_invitro_qsp <- function() {
  description <- "QSP. In vitro cellular ADC processing model for HER2-targeting antibody-drug conjugates (T-DM1 default; T-DXd variant supported via parameter overrides). 13 ODE states describing extracellular ADC / Ab / free payload, cell-surface HER2 with reversible ADC and antibody binding, endosomal HER2 species, endosomal / cytosolic payload, and cytosolic intracellular target (tubulin for DM1, TOPO-1 for DXd). Amounts in nmol per paper Modeling Convention; concentrations are amount/volume. Parameter set defaults to SK-BR-3 cell line with T-DM1 (Erickson 2012 in vitro incubation)."
  reference   <- "Scheuher B, Ghusinga KR, McGirr K, Nowak M, Panday S, Apgar J, Subramanian K, Betts A. Towards a platform quantitative systems pharmacology (QSP) model for preclinical to clinical translation of antibody drug conjugates (ADCs). J Pharmacokinet Pharmacodyn. 2023;51(1):5-30. doi:10.1007/s10928-023-09884-6. In vitro cellular model = Tables S1a, S2a, S3a-b."
  vignette    <- "Scheuher_2023_ADC_platform_qsp"
  units       <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    adc_ext       = list(analyte = "ADC", units = "nmol", specimen = "plasma", verified = FALSE),
    ab_ext        = list(analyte = "Ab", units = "nmol", specimen = "plasma", verified = FALSE),
    pl_ext        = list(analyte = "payload", units = "nmol", specimen = "plasma", verified = FALSE),
    her2          = list(analyte = "HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    her2_adc      = list(analyte = "ADC-HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    her2_ab       = list(analyte = "Ab-HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    her2_endo     = list(analyte = "endosomal HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    her2_adc_endo = list(analyte = "ADC-endosomal HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    her2_ab_endo  = list(analyte = "Ab-endosomal HER2", units = "nmol", specimen = "administration site", verified = FALSE),
    pl_endo       = list(analyte = "payload", units = "nmol", specimen = "endosome", verified = FALSE),
    pl_cyto       = list(analyte = "payload", units = "nmol", specimen = "plasma", verified = FALSE),
    t_cyto        = list(analyte = "target", units = "nmol", specimen = "plasma", verified = FALSE),
    tpl_cyto      = list(analyte = "payload-target complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "in vitro (SK-BR-3 cell line; BT-474 and MCF-7-neo/HER2 alternates supported)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = "HER2-expressing breast cancer cell lines under 2-hour ice incubation with radiolabelled T-DM1, washed, then measured for intracellular and extracellular DM1 catabolites over ~144 hours",
    dose_range     = "T-DM1 spike at experimental start (SK-BR-3 default from Erickson 2012 Figure 3); no repeated dosing",
    regions        = "Preclinical in vitro (Erickson et al. 2012)",
    notes          = "Cell-line-specific alternatives: V_cell = 3.82e-12 L (SK-BR-3, default), 3.96e-12 L (BT-474), 3.65e-12 L (MCF-7-neo/HER2); RPC_HER2 = 1.12e6 (SK-BR-3), 4.7e5 (BT-474), 1.3e6 (MCF-7-neo/HER2). Payload identity: DM1 (default; T_per_cell=65 nM tubulin, K_D_PL=930 nM). T-DXd alternative: k_in_PL=0.0128 /s, k_out_PL=0.00898 /s, T_per_cell=297 nM TOPO-1, K_D_PL=307 nM, DAR=8, kdec (mouse-level)=1.7e-7 /s. See vignette Errata for full parameter substitutions."
  )

  ini({
    # -------------------------------------------------------------------
    # Structural parameters - T-DM1 in SK-BR-3 (default in vitro parameterisation)
    # Source: Scheuher 2023 Supplement Table S2a
    # -------------------------------------------------------------------

    # Deconjugation rate (1/s), converted to 1/hour for consistency with time units.
    # Table S2a: kdec = 6.41e-7 /s -> 6.41e-7 * 3600 = 2.308e-3 /hour
    lkdec <- log(2.308e-3); label("Deconjugation rate constant (1/h)")             # Table S2a: 6.41e-7 /s (T-DM1); fit to Erickson 2012 Figure 3

    # HER2 binding (association rate constant, in nM^-1 hour^-1; converted from /s)
    # Table S2a: k_on_Ab = 1e-4 /nM/s -> 1e-4 * 3600 = 0.36 /nM/hour
    lkonab <- fixed(log(0.36)); label("HER2 binding association rate constant (1/nM/h)")  # Table S2a: 1e-4 /nM/s (Schlosshauer & Baker 2004 typical)

    # HER2 equilibrium binding constant (nM)
    lkdab <- log(0.314); label("HER2:ADC equilibrium binding constant K_D_Ab (nM)")   # Table S2a: 0.314 nM; fit to Austin 2004

    # HER2 receptor kinetics (1/s, converted to 1/hour)
    # Table S2a: k_endo_HER2 = 4.27e-5 /s -> 4.27e-5 * 3600 = 0.15372 /hour
    lkendoher2 <- log(0.15372); label("HER2 endocytosis rate constant (1/h)")      # Table S2a: 4.27e-5 /s; fit to Austin 2004
    lkrecher2  <- log(2.4e-5 * 3600); label("HER2 recycling rate constant (1/h)")  # Table S2a: 2.4e-5 /s = 0.0864 /hour
    lkdegher2  <- log(1.27e-4 * 3600); label("HER2 endosomal degradation rate constant (1/h)")  # Table S2a: 1.27e-4 /s = 0.4572 /hour

    # HER2:Ab complex kinetics assumed equal to free HER2 kinetics per Table S2a
    # (endocytosis, recycling, and degradation identical between HER2 and HER2:Ab)

    # Linker cleavage rate constant (1/hour); T-DM1 has non-cleavable linker so 0
    lkcleave <- fixed(log(1e-12)); label("Endosomal linker cleavage rate (1/h) - ~0 for T-DM1 non-cleavable linker")  # Table S2a: 0 /s

    # Payload target binding
    # k_on_PL = 1e-3 /nM/s -> 3.6 /nM/hour (typical binding, per Frost 2016 / Guan 2015)
    lkonpl <- fixed(log(3.6)); label("Payload:target association rate constant (1/nM/h)")  # Table S2a: 1e-3 /nM/s
    lkdpl  <- log(930); label("Payload:target equilibrium binding K_D_PL (nM)")               # Table S2a: 930 nM DM1:tubulin (Lopus 2010)

    # Cytosolic target concentration per cell (nM)
    lTpercell <- log(65); label("Intracellular target concentration T_per_cell (nM)")  # Table S2a: 65 nM tubulin (Shah 2012 estimate)

    # Payload cellular flux (1/s -> 1/hour)
    lkinpl <- log(5.95e-5 * 3600); label("Payload influx rate constant into cell (1/h)")  # Table S2a: 5.95e-5 /s = 0.2142 /hour (Khera 2018)
    lkoutpl <- log(3.95e-5 * 3600); label("Payload efflux rate constant from cell (1/h)") # Table S2a: 3.95e-5 /s = 0.1422 /hour (Khera 2018)

    # Cell/volume parameters (SK-BR-3 in vitro, T-DM1 with Erickson 2012 experiment)
    lVcell <- fixed(log(3.82e-12)); label("Cellular volume (L)")                        # Table S2a: 3.82e-12 L SK-BR-3
    lVmedia <- fixed(log(1e-4)); label("In vitro media volume (L)")                     # Table S2a: 1e-4 L (Erickson 2012 experiment)
    lNcell <- fixed(log(1e7)); label("Number of cells in in vitro incubation")          # Table S2a: 1e7 cells (Erickson 2012)
    lRPCher2 <- fixed(log(1.12e6)); label("HER2 receptors per cell RPC_HER2_tumor")     # Table S2a: 1.12e6 receptors/cell SK-BR-3

    # DAR (drug-to-antibody ratio)
    lDAR <- fixed(log(3.5)); label("Drug-to-antibody ratio (DAR)")                       # Table S2a: 3.5 for T-DM1 (Kadcyla BLA 2013)

    # Residual error placeholders - the in vitro model has no fitted subject-level
    # residual error in the paper (deterministic simulation vs. observed profiles),
    # but nlmixr2 conventions require at least one error term. Set a small
    # proportional error to allow simulation with observed noise levels.
    propSd <- fixed(0.1); label("Proportional residual error (fraction) - not from paper; 10% for simulation")  # Assumption: paper reports deterministic fits, no residual error term
  })

  model({
    # -------------------------------------------------------------------
    # Individual parameter transformations
    # -------------------------------------------------------------------
    kdec      <- exp(lkdec)
    konab     <- exp(lkonab)
    kdab      <- exp(lkdab)
    koffab    <- konab * kdab                       # Table S2a inferred: k_off,Ab = k_on,Ab * K_D,Ab
    kendoher2 <- exp(lkendoher2)
    krecher2  <- exp(lkrecher2)
    kdegher2  <- exp(lkdegher2)
    kcleave   <- exp(lkcleave)
    konpl     <- exp(lkonpl)
    kdpl      <- exp(lkdpl)
    koffpl    <- konpl * kdpl                       # Table S2a inferred: k_off,PL = k_on,PL * K_D,PL
    Tpercell  <- exp(lTpercell)
    kinpl     <- exp(lkinpl)
    koutpl    <- exp(lkoutpl)
    Vcell     <- exp(lVcell)
    Vmedia    <- exp(lVmedia)
    Ncell     <- exp(lNcell)
    RPCher2   <- exp(lRPCher2)
    DAR       <- exp(lDAR)

    # -------------------------------------------------------------------
    # Derived initial conditions and synthesis rate for HER2 steady state
    # Source: Table S2a inferred parameters
    # -------------------------------------------------------------------
    # Cell-surface HER2 amount per cell at steady state (before drug):
    #   HER2_0_per_cell = RPCher2 / 6.022e23 * 1e9  (receptors -> nmol)
    HER2_0_per_cell_nmol <- RPCher2 / 6.022e23 * 1e9

    # Zeroth-order synthesis rate (Table S2a inferred):
    #   k_synth = RPC * k_endo * k_deg / (k_deg + k_rec) - here in nmol/hour for all cells
    ksynth_nmol_per_hour <- HER2_0_per_cell_nmol * kendoher2 * kdegher2 / (kdegher2 + krecher2) * Ncell

    # Endosomal HER2 amount per cell at steady state (Table S2a inferred):
    #   RPC_HER2_endo = k_endo / (k_rec + k_deg) * RPCher2
    HER2_endo_0_per_cell_nmol <- kendoher2 / (krecher2 + kdegher2) * HER2_0_per_cell_nmol

    # -------------------------------------------------------------------
    # ODE system (Table S3a-b, in vitro cellular reactions and ODEs)
    # States are AMOUNTS in nmol (paper Modeling Convention).
    # -------------------------------------------------------------------

    # Deconjugation (v1): extracellular ADC -> antibody + DAR*payload
    v1 <- kdec * adc_ext

    # HER2 synthesis (v2): zeroth order into cell-surface HER2 pool
    v2 <- ksynth_nmol_per_hour

    # HER2 + ADC binding (v3): forward and reverse
    v3f <- konab / Vmedia * her2 * adc_ext
    v3r <- koffab * her2_adc

    # HER2 + Ab binding (v4): forward and reverse
    v4f <- konab / Vmedia * her2 * ab_ext
    v4r <- koffab * her2_ab

    # HER2 internalization / recycling (v5)
    v5f <- kendoher2 * her2
    v5r <- krecher2  * her2_endo

    # HER2:ADC internalization / recycling (v6)
    v6f <- kendoher2 * her2_adc
    v6r <- krecher2  * her2_adc_endo

    # HER2:Ab internalization / recycling (v7)
    v7f <- kendoher2 * her2_ab
    v7r <- krecher2  * her2_ab_endo

    # Endosomal HER2 degradation (v8)
    v8 <- kdegher2 * her2_endo

    # Endosomal HER2:ADC degradation -> DAR * payload_endo (v9)
    v9 <- kdegher2 * her2_adc_endo

    # Endosomal HER2:Ab degradation (v10)
    v10 <- kdegher2 * her2_ab_endo

    # Endosomal linker cleavage: HER2:ADC_endo -> HER2:Ab_endo + DAR*PL_endo (v11)
    v11 <- kcleave * her2_adc_endo

    # Payload endosome -> cytosol (v12)
    v12 <- kinpl * pl_endo

    # Payload:target binding in cytosol (v13); in vitro Pc_PL = 1 per Table S3a footnote
    v13f <- konpl / Vcell * t_cyto * pl_cyto        # Pc_PL = 1 in vitro (Table S3a footnote)
    v13r <- koffpl * tpl_cyto

    # Payload cytosol <-> extracellular flux (v14)
    v14f <- koutpl * pl_cyto
    v14r <- kinpl  * pl_ext

    # -------------------------------------------------------------------
    # ODEs (Table S3b)
    # -------------------------------------------------------------------
    d/dt(adc_ext)       <- -v1 - v3f + v3r                                # (1) Table S3b
    d/dt(ab_ext)        <-  v1 - v4f + v4r                                # (2)
    d/dt(pl_ext)        <-  v1 * DAR + v14f - v14r                        # (3)
    d/dt(her2)          <-  v2 - v5f + v5r - v3f + v3r - v4f + v4r        # (4)
    d/dt(her2_adc)      <-  v3f - v3r - v6f + v6r                         # (5)
    d/dt(her2_ab)       <-  v4f - v4r - v7f + v7r                         # (6)
    d/dt(her2_endo)     <-  v5f - v5r - v8                                # (7)
    d/dt(her2_adc_endo) <-  v6f - v6r - v9 - v11                          # (8)
    d/dt(her2_ab_endo)  <-  v7f - v7r - v10 + v11                         # (9)
    d/dt(pl_endo)       <-  (v9 + v11) * DAR - v12                        # (10)
    d/dt(pl_cyto)       <-  v12 - v13f + v13r - v14f + v14r               # (11)
    d/dt(t_cyto)        <- -v13f + v13r                                   # (12)
    d/dt(tpl_cyto)      <-  v13f - v13r                                   # (13; Table S3b Eqn 13 typo: printed dT_cyto/dt=v13f-v13r on this row but state is T:PL_cyto)

    # -------------------------------------------------------------------
    # Initial conditions
    # -------------------------------------------------------------------
    her2(0)       <- HER2_0_per_cell_nmol * Ncell
    her2_endo(0)  <- HER2_endo_0_per_cell_nmol * Ncell
    t_cyto(0)     <- Tpercell * Vcell * Ncell           # amount over all cells (Vcell*Ncell = total intracellular volume)

    # -------------------------------------------------------------------
    # Observations
    # -------------------------------------------------------------------
    # Intracellular unconjugated payload concentration (nM) in cells
    Cpl_intra <- (pl_cyto + tpl_cyto) / (Vcell * Ncell)
    # Extracellular ADC concentration (nM) in media
    Cc <- adc_ext / Vmedia
    # Extracellular Ab (naked trastuzumab) concentration (nM) in media
    Cab <- ab_ext / Vmedia
    # Extracellular unconjugated payload (nM) in media
    Cpl_ext <- pl_ext / Vmedia
    # Total intracellular DM1 (bound + free) in nmol / cell (compare to Erickson 2012 Figure 3)
    intra_dm1_per_cell <- (pl_cyto + tpl_cyto + pl_endo) / Ncell

    Cc ~ prop(propSd)
  })
}
