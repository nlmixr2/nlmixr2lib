Derippe_2024_venetoclax_apoptosis_qsp <- function() {
  description <- "QSP. In vitro (SU-DHL-4 and KARPAS-422 B-cell lymphoma cell lines) with mouse xenograft extension. Modified Lindner BCL2-family apoptosis network (57 ODE states) driving single-cell mitochondrial outer membrane permeabilisation (MOMP) in response to the BH3-mimetics venetoclax (anti-BCL2), A-1155463 (anti-BCL-XL) and A-1210477 (anti-MCL-1). Antiapoptotic proteins BCL2 / BCL-XL / MCL-1 sequester the BH3-only activators BIM / tBID / PUMA / NOXA and cap activated BAX / BAK; activated BAX and BAK oligomerise to dodecamers and apoptosis is triggered once more than 10 percent of total BAX + BAK sits in hexamers or larger (state Pore). Each drug removes its target by a second-order kill process. The eight per-cell initial protein levels define a virtual cell; a virtual tumor is 100 such cells calibrated to cell viability assays. In vivo extension adds mouse oral venetoclax PK and a kinetic-pharmacodynamic CDK9-inhibitor (A-1592668) arm that shuts down MCL-1 production."
  reference <- paste(
    "Derippe T, Fouliard S, Decleves X, Mager DE.",
    "Quantitative systems pharmacology modeling of tumor heterogeneity in response to",
    "BH3-mimetics using virtual tumors calibrated with cell viability assays.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1252-1263. doi:10.1002/psp4.13158.",
    "Every ODE and rate constant below is transcribed from the authors' own deposited RxODE",
    "implementation, file 0_Lindner_model_PaSM_config.R in the repository named by the paper's",
    "Data Availability Statement (https://github.com/Thibaudpmx/Virtual_tumor_publication;",
    "Zenodo snapshot doi:10.5281/zenodo.10826315).",
    "The underlying network topology and rate constants originate from",
    "Lindner AU, Concannon CG, Boukes GJ, et al. Systems analysis of BCL2 protein family",
    "interactions establishes a model to predict responses to chemotherapy.",
    "Cancer Res. 2013;73(2):519-528. doi:10.1158/0008-5472.CAN-12-2269.",
    sep = " "
  )
  vignette <- "Derippe_2024_BH3_mimetics_virtual_tumors"

  # Units are mixed by construction in the source implementation and are kept as
  # published: protein species are in nM, BH3-mimetic concentrations are in uM
  # (the kill constants k2_*_I are 1/uM/h), and time is in hours. The in vitro
  # drug exposures are NOT dosing events -- they are set through the parameters
  # Bcl2_I0 / Bclxl_I0 / Mcl1_I0. Only the in vivo arm uses dosing records
  # (oral venetoclax into Veneto_gut in mg/kg, and A-1592668 into A15).
  units <- list(time = "h", dosing = "mg/kg", concentration = "nM")

  # buildModelDb() auto-detects dosing compartments only when they are literally
  # named `depot` or `central`. This model's dosing states are `Veneto_gut` (oral
  # venetoclax) and `A15` (the A-1592668 kinetic-pharmacodynamic arm), so they are
  # declared explicitly -- without this the registry would advertise the model as
  # having no dosing compartment at all.
  dosing <- c("Veneto_gut", "A15")

  covariateData <- list()

  covariatesDataExcluded <- list()

  paper_specific_compartments <- c(
    "Veneto_gut", "Veneto_central", "A15", "Bcl2_I", "Bclxl_I", "Mcl1_I",
    "Bcl2", "Bclxl", "Mcl1", "BIM", "tBID", "PUMA", "NOXA",
    "Bcl2_BIM", "Bclxl_BIM", "Mcl1_BIM", "Bcl2_tBID", "Bclxl_tBID", "Mcl1_tBID",
    "Bcl2_PUMA", "Bclxl_PUMA", "Mcl1_PUMA", "Bcl2_NOXA", "Bclxl_NOXA", "Mcl1_NOXA",
    "Bcl2_BAKa", "BAKa", "Bcl2_BAXma", "BAXma", "Bclxl_BAKa", "Bclxl_BAXma",
    "Mcl1_BAKa", "Mcl1_BAXma", "BIM_BAXc", "BAXc", "tBID_BAXc", "PUMA_BAXc",
    "BIM_BAK", "BAK", "tBID_BAK", "PUMA_BAK", "BAXca", "BAK_VDAC2", "VDAC2",
    "BAK2", "BAK4", "BAK6", "BAK8", "BAK10", "BAK12",
    "BAX2", "BAX4", "BAX6", "BAX8", "BAX10", "BAX12", "TimeAbove"
  )

  compartmentData <- list(
    Veneto_gut     = list(analyte = "venetoclax", units = "mg/kg", specimen = "administration site", verified = TRUE),
    Veneto_central = list(analyte = "venetoclax", units = "mg/kg", specimen = "plasma", verified = TRUE),
    A15            = list(analyte = "A-1592668", units = "mg/kg", specimen = "not applicable", verified = TRUE),
    Bcl2_I         = list(analyte = "venetoclax", units = "uM", specimen = "tumor", verified = TRUE),
    Bclxl_I        = list(analyte = "A-1155463", units = "uM", specimen = "tumor", verified = TRUE),
    Mcl1_I         = list(analyte = "A-1210477", units = "uM", specimen = "tumor", verified = TRUE),
    Bcl2           = list(analyte = "BCL2 protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl          = list(analyte = "BCL-XL protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1           = list(analyte = "MCL-1 protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    BIM            = list(analyte = "BIM protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    tBID           = list(analyte = "truncated BID protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    PUMA           = list(analyte = "PUMA protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    NOXA           = list(analyte = "NOXA protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_BIM       = list(analyte = "BCL2:BIM complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_BIM      = list(analyte = "BCL-XL:BIM complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_BIM       = list(analyte = "MCL-1:BIM complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_tBID      = list(analyte = "BCL2:tBID complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_tBID     = list(analyte = "BCL-XL:tBID complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_tBID      = list(analyte = "MCL-1:tBID complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_PUMA      = list(analyte = "BCL2:PUMA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_PUMA     = list(analyte = "BCL-XL:PUMA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_PUMA      = list(analyte = "MCL-1:PUMA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_NOXA      = list(analyte = "BCL2:NOXA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_NOXA     = list(analyte = "BCL-XL:NOXA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_NOXA      = list(analyte = "MCL-1:NOXA complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_BAKa      = list(analyte = "BCL2:activated-BAK complex", units = "nM", specimen = "tumor", verified = TRUE),
    BAKa           = list(analyte = "activated BAK protein", units = "nM", specimen = "tumor", verified = TRUE),
    Bcl2_BAXma     = list(analyte = "BCL2:activated-membrane-BAX complex", units = "nM", specimen = "tumor", verified = TRUE),
    BAXma          = list(analyte = "activated membrane-bound BAX protein", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_BAKa     = list(analyte = "BCL-XL:activated-BAK complex", units = "nM", specimen = "tumor", verified = TRUE),
    Bclxl_BAXma    = list(analyte = "BCL-XL:activated-membrane-BAX complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_BAKa      = list(analyte = "MCL-1:activated-BAK complex", units = "nM", specimen = "tumor", verified = TRUE),
    Mcl1_BAXma     = list(analyte = "MCL-1:activated-membrane-BAX complex", units = "nM", specimen = "tumor", verified = TRUE),
    BIM_BAXc       = list(analyte = "BIM:cytosolic-BAX activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    BAXc           = list(analyte = "cytosolic (inactive) BAX protein", units = "nM", specimen = "tumor", verified = TRUE),
    tBID_BAXc      = list(analyte = "tBID:cytosolic-BAX activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    PUMA_BAXc      = list(analyte = "PUMA:cytosolic-BAX activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    BIM_BAK        = list(analyte = "BIM:BAK activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    BAK            = list(analyte = "inactive BAK protein", units = "nM", specimen = "tumor", verified = TRUE),
    tBID_BAK       = list(analyte = "tBID:BAK activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    PUMA_BAK       = list(analyte = "PUMA:BAK activation complex", units = "nM", specimen = "tumor", verified = TRUE),
    BAXca          = list(analyte = "activated cytosolic BAX protein", units = "nM", specimen = "tumor", verified = TRUE),
    BAK_VDAC2      = list(analyte = "BAK:VDAC2 complex", units = "nM", specimen = "tumor", verified = TRUE),
    VDAC2          = list(analyte = "VDAC2 protein (free)", units = "nM", specimen = "tumor", verified = TRUE),
    BAK2           = list(analyte = "BAK dimer", units = "nM", specimen = "tumor", verified = TRUE),
    BAK4           = list(analyte = "BAK tetramer", units = "nM", specimen = "tumor", verified = TRUE),
    BAK6           = list(analyte = "BAK hexamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAK8           = list(analyte = "BAK octamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAK10          = list(analyte = "BAK decamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAK12          = list(analyte = "BAK dodecamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAX2           = list(analyte = "BAX dimer", units = "nM", specimen = "tumor", verified = TRUE),
    BAX4           = list(analyte = "BAX tetramer", units = "nM", specimen = "tumor", verified = TRUE),
    BAX6           = list(analyte = "BAX hexamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAX8           = list(analyte = "BAX octamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAX10          = list(analyte = "BAX decamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    BAX12          = list(analyte = "BAX dodecamer (pore-forming)", units = "nM", specimen = "tumor", verified = TRUE),
    TimeAbove      = list(analyte = "cumulative time with Pore above threshold", units = "h", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "in vitro (SU-DHL-4 and KARPAS-422 diffuse large B-cell lymphoma cell lines); in vivo extension in mouse (SU-DHL-4 subcutaneous xenograft)",
    n_subjects     = 100L,
    n_studies      = 3L,
    disease_state  = "B-cell non-Hodgkin lymphoma (SU-DHL-4 and KARPAS-422 germinal-centre DLBCL lines)",
    dose_range     = "In vitro 48 h exposure: venetoclax or A-1155463 at 0, 0.08, 0.16, 0.32, 0.64, 1.3, 2.6, 5, 10, 20 uM, each crossed with A-1210477 at 0, 5, 10, 15 uM. In vivo: venetoclax 50 mg/kg PO QD x 21 days; A-1592668 1.5 mg/kg PO three times weekly x 3 weeks.",
    regions        = "Preclinical (literature-digitized in vitro and in vivo data)",
    notes          = "n_subjects = 100 refers to the number of virtual cells per calibrated virtual tumor (arbitrarily fixed at 100 by the authors), not to animals or patients. Cell viability data were digitized from Phillips DC et al., Blood Cancer J. 2015;5:e368 (80 calibration points per cell line: 40 venetoclax +/- A-1210477 and 40 A-1155463 +/- A-1210477). Mouse tumor growth inhibition data were digitized from Phillips DC et al., Leukemia. 2020;34(6):1646-1657. Mouse venetoclax PK was fit to data digitized from Eisenmann 2020 (J Chromatogr B 1152:122176) and Salem 2021 (Mol Cancer Ther 20(6):999-1008); see modellib('Derippe_2024_venetoclax_mouse'). 7,703,029 virtual cells were generated in total, of which ~2.4 million showed spontaneous apoptosis and were discarded. The default per-cell protein levels shipped below are virtual cell 4272, the modal cell of BOTH published three-drug calibrated virtual tumors (11/100 cells in SU-DHL-4, 9/100 in KARPAS-422)."
  )

  ini({
    # =====================================================================
    # 1. Per-cell protein levels -- the eight "parameters of interest".
    #    These define a virtual cell and are the intended user-varying
    #    inputs; a virtual tumor is 100 cells drawn from a calibrated set.
    #    Shipped defaults = virtual cell 4272 from the deposited
    #    celltheque (calibrated_VT/VT_both_cell_line_1.RDS, the SU-DHL-4
    #    three-drug virtual tumor), the modal cell of both published VTs.
    #    NOT wrapped in fixed(): the VT calibration selects these per cell,
    #    and downstream users are expected to replace them.
    # =====================================================================
    Bcl20  <- 1356.45; label("Initial free BCL2 concentration (nM)")    # deposited celltheque, cellid 4272 (paper Figure 2a value ranges: 20-2000 nM)
    Bclxl0 <- 688.19;  label("Initial free BCL-XL concentration (nM)")  # deposited celltheque, cellid 4272
    Mcl10  <- 424.63;  label("Initial free MCL-1 concentration (nM)")   # deposited celltheque, cellid 4272
    BIM0   <- 303.71;  label("Initial free BIM concentration (nM)")     # deposited celltheque, cellid 4272
    PUMA0  <- 133.76;  label("Initial free PUMA concentration (nM)")    # deposited celltheque, cellid 4272
    NOXA0  <- 141.20;  label("Initial free NOXA concentration (nM)")    # deposited celltheque, cellid 4272
    BAK0   <- 29.62;   label("Initial total BAK concentration (nM); also sets initial VDAC2")  # deposited celltheque, cellid 4272
    BAXc0  <- 652.49;  label("Initial cytosolic BAX concentration (nM)")                       # deposited celltheque, cellid 4272

    # =====================================================================
    # 2. In vitro BH3-mimetic exposures (uM). Set these to the assay
    #    concentration; they are held constant because ke_*_I = 0.
    # =====================================================================
    Bcl2_I0  <- 0; label("Venetoclax concentration in the in vitro assay (uM)")   # Methods: 0-20 uM, doubling
    Bclxl_I0 <- 0; label("A-1155463 concentration in the in vitro assay (uM)")    # Methods: 0-20 uM, doubling
    Mcl1_I0  <- 0; label("A-1210477 concentration in the in vitro assay (uM)")    # deposited data: 0, 5, 10, 15 uM

    # =====================================================================
    # 3. Drug action. Second-order kill constants set to 10 /uM/h for all
    #    three BH3-mimetics; no intrinsic drug elimination in vitro.
    # =====================================================================
    k2_Bcl2_I  <- fixed(10); label("Second-order kill rate of venetoclax on BCL2 (1/uM/h)")     # Supplement "QSP model: Modification from original version"; config parameters_default_values
    k2_Bclxl_I <- fixed(10); label("Second-order kill rate of A-1155463 on BCL-XL (1/uM/h)")    # Supplement, as above
    k2_Mcl1_I  <- fixed(10); label("Second-order kill rate of A-1210477 on MCL-1 (1/uM/h)")     # Supplement, as above
    ke_BCl2_I  <- fixed(0);  label("Elimination rate of venetoclax in the in vitro system (1/h)")  # config parameters_default_values (no drug elimination in vitro)
    ke_BClxl_I <- fixed(0);  label("Elimination rate of A-1155463 in the in vitro system (1/h)")   # config parameters_default_values
    ke_Mcl1_I  <- fixed(0);  label("Elimination rate of A-1210477 in vitro; also the A-1592668 K-PD elimination rate (1/h)")  # config parameters_default_values; reused for d/dt(A15) -- see vignette Errata

    # =====================================================================
    # 4. System / simulation controls
    # =====================================================================
    switchE       <- fixed(1);    label("Homeostasis switch: 1 enables BH3-only production and BAX/BAK turnover")  # Supplement: the Derippe modification of the Lindner model
    tburn         <- fixed(700);  label("Burn-in duration before drug exposure (h)")            # Supplement: "This burn-in phase is performed for 700 h"
    poreThreshold <- fixed(10);   label("Percent of total BAX + BAK in hexamers or larger that triggers MOMP (%)") # Supplement + Methods: apoptosis once 10% of BAX/BAK are in pores
    kelimBAXBAK   <- fixed(0.014); label("First-order turnover rate of every BAX/BAK-containing species (1/h)")    # config parameters_default_values; see vignette Errata (supplement prose says half-life 22 h, i.e. 0.0315 1/h)

    # =====================================================================
    # 5. Protein degradation rate constants (Lindner 2013, transcribed
    #    from the deposited RxODE model). All 1/h.
    # =====================================================================
    kdeg_Bcl2  <- fixed(0.139);        label("Degradation rate of free BCL2 (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_Bclxl <- fixed(0.139);        label("Degradation rate of free BCL-XL (1/h)")  # 0_Lindner_model_PaSM_config.R
    kdeg_Mcl1  <- fixed(0.00068 * 3600); label("Degradation rate of free MCL-1 (1/h)") # 0_Lindner_model_PaSM_config.R: 0.00068 1/s x 3600
    kdeg_BIM   <- fixed(0.173);        label("Degradation rate of free BIM (1/h)")     # 0_Lindner_model_PaSM_config.R
    kdeg_tBID  <- fixed(0.554);        label("Degradation rate of free tBID (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_PUMA  <- fixed(0.204);        label("Degradation rate of free PUMA (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_NOXA  <- fixed(0.695);        label("Degradation rate of free NOXA (1/h)")    # 0_Lindner_model_PaSM_config.R

    # Complex degradation rate constants (1/h)
    kdeg_Bcl2_BIM   <- fixed(0.554); label("Degradation rate of the BCL2:BIM complex (1/h)")     # 0_Lindner_model_PaSM_config.R
    kdeg_Bclxl_BIM  <- fixed(0.554); label("Degradation rate of the BCL-XL:BIM complex (1/h)")   # 0_Lindner_model_PaSM_config.R
    kdeg_Mcl1_BIM   <- fixed(0.277); label("Degradation rate of the MCL-1:BIM complex (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_Bcl2_tBID  <- fixed(0.554); label("Degradation rate of the BCL2:tBID complex (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_Bclxl_tBID <- fixed(0.554); label("Degradation rate of the BCL-XL:tBID complex (1/h)")  # 0_Lindner_model_PaSM_config.R
    kdeg_Mcl1_tBID  <- fixed(0.554); label("Degradation rate of the MCL-1:tBID complex (1/h)")   # 0_Lindner_model_PaSM_config.R
    kdeg_Bcl2_PUMA  <- fixed(0.554); label("Degradation rate of the BCL2:PUMA complex (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_Bclxl_PUMA <- fixed(0.554); label("Degradation rate of the BCL-XL:PUMA complex (1/h)")  # 0_Lindner_model_PaSM_config.R
    kdeg_Mcl1_PUMA  <- fixed(0.277); label("Degradation rate of the MCL-1:PUMA complex (1/h)")   # 0_Lindner_model_PaSM_config.R
    kdeg_Bcl2_NOXA  <- fixed(0.554); label("Degradation rate of the BCL2:NOXA complex (1/h)")    # 0_Lindner_model_PaSM_config.R
    kdeg_Bclxl_NOXA <- fixed(0.554); label("Degradation rate of the BCL-XL:NOXA complex (1/h)")  # 0_Lindner_model_PaSM_config.R
    kdeg_Mcl1_NOXA  <- fixed(0.925); label("Degradation rate of the MCL-1:NOXA complex (1/h)")   # 0_Lindner_model_PaSM_config.R

    # =====================================================================
    # 6. Sequestration association rate constants (1/nM/h). Each
    #    antiapoptotic protein has a single dissociation rate constant
    #    shared across all of its partners (kbackward_Bcl2 / _Bclxl /
    #    _Mcl1), exactly as in the deposited model.
    # =====================================================================
    kforward_Bcl2_BIM   <- fixed(0.108);    label("Association rate of BCL2 with BIM (1/nM/h)")       # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_BIM  <- fixed(1.98);     label("Association rate of BCL-XL with BIM (1/nM/h)")     # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_BIM   <- fixed(4.68);     label("Association rate of MCL-1 with BIM (1/nM/h)")      # 0_Lindner_model_PaSM_config.R
    kforward_Bcl2_tBID  <- fixed(0.126);    label("Association rate of BCL2 with tBID (1/nM/h)")      # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_tBID <- fixed(0.0583);   label("Association rate of BCL-XL with tBID (1/nM/h)")    # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_tBID  <- fixed(0.0947);   label("Association rate of MCL-1 with tBID (1/nM/h)")     # 0_Lindner_model_PaSM_config.R
    kforward_Bcl2_PUMA  <- fixed(0.028);    label("Association rate of BCL2 with PUMA (1/nM/h)")      # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_PUMA <- fixed(0.311);    label("Association rate of BCL-XL with PUMA (1/nM/h)")    # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_PUMA  <- fixed(0.493);    label("Association rate of MCL-1 with PUMA (1/nM/h)")     # 0_Lindner_model_PaSM_config.R
    kforward_Bcl2_NOXA  <- fixed(0.00262);  label("Association rate of BCL2 with NOXA (1/nM/h)")      # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_NOXA <- fixed(0.000158); label("Association rate of BCL-XL with NOXA (1/nM/h)")    # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_NOXA  <- fixed(0.0237);   label("Association rate of MCL-1 with NOXA (1/nM/h)")     # 0_Lindner_model_PaSM_config.R

    kforward_Bcl2_BAKa   <- fixed(5.04e-05); label("Association rate of BCL2 with activated BAK (1/nM/h)")        # 0_Lindner_model_PaSM_config.R
    kforward_Bcl2_BAXma  <- fixed(0.0336);   label("Association rate of BCL2 with membrane BAX (1/nM/h)")         # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_BAKa  <- fixed(0.0198);   label("Association rate of BCL-XL with activated BAK (1/nM/h)")      # 0_Lindner_model_PaSM_config.R
    kforward_Bclxl_BAXma <- fixed(0.00186);  label("Association rate of BCL-XL with membrane BAX (1/nM/h)")       # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_BAKa   <- fixed(0.117);    label("Association rate of MCL-1 with activated BAK (1/nM/h)")       # 0_Lindner_model_PaSM_config.R
    kforward_Mcl1_BAXma  <- fixed(9.36e-05); label("Association rate of MCL-1 with membrane BAX (1/nM/h)")        # 0_Lindner_model_PaSM_config.R

    kbackward_Bcl2  <- fixed(0.504); label("Dissociation rate of every BCL2 complex (1/h)")    # 0_Lindner_model_PaSM_config.R: kbackward_Bcl2_* all 0.504
    kbackward_Bclxl <- fixed(1.58);  label("Dissociation rate of every BCL-XL complex (1/h)")  # 0_Lindner_model_PaSM_config.R: kbackward_Bclxl_* all 1.58
    kbackward_Mcl1  <- fixed(0.936); label("Dissociation rate of every MCL-1 complex (1/h)")   # 0_Lindner_model_PaSM_config.R: kbackward_Mcl1_* all 0.936

    # =====================================================================
    # 7. BAX / BAK activation. The deposited model uses one association
    #    rate (0.00925), one dissociation rate (0.925) and one catalytic
    #    activation rate (418) for all six activator:effector pairs.
    # =====================================================================
    kforward_act  <- fixed(0.00925); label("Association rate of a BH3-only activator with BAX or BAK (1/nM/h)")  # 0_Lindner_model_PaSM_config.R: kforward_{BIM,tBID,PUMA}_{BAXc,BAK} all 0.00925
    kbackward_act <- fixed(0.925);   label("Dissociation rate of an activator:effector complex (1/h)")           # 0_Lindner_model_PaSM_config.R: kbackward_{BIM,tBID,PUMA}_{BAXc,BAK} all 0.925
    k_activation  <- fixed(418);     label("Catalytic activation rate of BAX or BAK from the activator complex (1/h)")  # 0_Lindner_model_PaSM_config.R: k_{BIM,tBID,PUMA}_{BAK,BAXc} and k_BAXca all 418

    kforward_BAK_VDAC2  <- fixed(0.00832); label("Association rate of BAK with VDAC2 (1/nM/h)")  # 0_Lindner_model_PaSM_config.R
    kbackward_BAK_VDAC2 <- fixed(8.32);    label("Dissociation rate of the BAK:VDAC2 complex (1/h)")  # 0_Lindner_model_PaSM_config.R

    # =====================================================================
    # 8. Oligomerisation. The deposited model assigns the same pair of
    #    constants to every BAK and BAX oligomerisation step (it literally
    #    re-assigns kforward_BAK8 / kforward_BAX12 etc. more than once with
    #    the same value), so they are carried here as one shared pair.
    # =====================================================================
    kforward_olig  <- fixed(0.0461); label("Association rate for every BAX/BAK oligomerisation step (1/nM/h)")  # 0_Lindner_model_PaSM_config.R: kforward_BAK2..BAK12 and kforward_BAX2..BAX12 all 0.0461
    kbackward_olig <- fixed(0.695);  label("Dissociation rate for every BAX/BAK oligomer (1/h)")                # 0_Lindner_model_PaSM_config.R: kbackward_BAK2..BAK12 and kbackward_BAX2..BAX12 all 0.695

    # =====================================================================
    # 9. In vivo extension: mouse oral venetoclax PK and the A-1592668
    #    K-PD arm. PK values are transferred from the paper's own mouse PK
    #    sub-model without re-fitting -- see
    #    modellib("Derippe_2024_venetoclax_mouse") for the full covariate
    #    model. The values hardcoded in the deposited QSP file are the
    #    male / Eisenmann column of Supplement Table "Mice PK modeling".
    # =====================================================================
    ka_Veneto  <- fixed(0.856); label("Oral absorption rate of venetoclax in mouse (1/h); taken from the paper's mouse PK sub-model")  # Supplement Table "Mice PK modeling" (male, Eisenmann)
    Cl_Veneto  <- fixed(0.449); label("Venetoclax clearance in mouse (L/h/kg); taken from the paper's mouse PK sub-model")             # Supplement Table "Mice PK modeling"; header prints L/h, see vignette Errata
    Vd_Veneto  <- fixed(3.54);  label("Venetoclax volume of distribution in mouse (L/kg); taken from the paper's mouse PK sub-model")  # Supplement Table "Mice PK modeling" (male, Eisenmann)
    ratioTumor <- fixed(1);     label("Plasma-to-tumor venetoclax conversion factor (unitless); assumed, set to 0.3 for the in vivo ABM") # Supplement "Minimal ABM limitations": 0.3 assumed with very high uncertainty; config default 1 (in vitro)
    hillA15    <- fixed(5);     label("Hill coefficient of A-1592668 inhibition of MCL-1 production (unitless)")  # Results: "sigmoidal Emax model with a Hill coefficient equal to 5"
    EC50A15    <- fixed(1e-6);  label("EC50 of A-1592668 inhibition of MCL-1 production (mg/kg); assumed small enough for complete inhibition") # Results: "an EC50 small enough such that drug concentrations were always above this value"; config value 1e-6
  })

  model({
    # -------------------------------------------------------------------
    # 1. Initial conditions. VDAC2 is initialised to the same total as BAK
    #    (one VDAC2 binding site per BAK molecule) exactly as deposited.
    # -------------------------------------------------------------------
    Bcl2_I(0)  <- Bcl2_I0
    Bclxl_I(0) <- Bclxl_I0
    Mcl1_I(0)  <- Mcl1_I0
    VDAC2(0)   <- BAK0
    BAK(0)     <- BAK0
    BAXc(0)    <- BAXc0
    NOXA(0)    <- NOXA0
    PUMA(0)    <- PUMA0
    BIM(0)     <- BIM0
    Mcl1(0)    <- Mcl10
    Bclxl(0)   <- Bclxl0
    Bcl2(0)    <- Bcl20

    # -------------------------------------------------------------------
    # 2. Zero-order production of the antiapoptotic proteins, set so that
    #    the drug-free system sits at its initial condition (the Derippe
    #    homeostasis modification of the Lindner model).
    # -------------------------------------------------------------------
    kprod_Bcl2  <- Bcl20 * kdeg_Bcl2
    kprod_Bclxl <- Bclxl0 * kdeg_Bclxl
    kprod_Mcl1  <- Mcl10 * kdeg_Mcl1

    # -------------------------------------------------------------------
    # 3. Reaction fluxes (numbering follows the deposited model R1-R74).
    # -------------------------------------------------------------------
    R1_deg_Bcl2  <- -kdeg_Bcl2 * Bcl2
    R2_prod_Bcl2 <- kprod_Bcl2
    R3_deg_Bclxl  <- -kdeg_Bclxl * Bclxl
    R4_prod_Bclxl <- kprod_Bclxl
    R5_deg_Mcl1  <- -kdeg_Mcl1 * Mcl1
    # A-1592668 (CDK9 inhibitor) shuts down MCL-1 production via a sigmoidal Imax term
    R6_prod_Mcl1 <- kprod_Mcl1 * (1 - A15^hillA15 / (A15^hillA15 + EC50A15^hillA15))
    R7_deg_BIM  <- -kdeg_BIM * BIM
    R8_deg_tBID <- -kdeg_tBID * tBID
    R9_deg_PUMA <- -kdeg_PUMA * PUMA
    R10_deg_NOXA <- -kdeg_NOXA * NOXA

    R11_deg_Bcl2_BIM     <- -kdeg_Bcl2_BIM * Bcl2_BIM
    R12_complex_Bcl2_BIM <- kforward_Bcl2_BIM * Bcl2 * BIM - kbackward_Bcl2 * Bcl2_BIM
    R13_deg_Bclxl_BIM     <- -kdeg_Bclxl_BIM * Bclxl_BIM
    R14_complex_Bclxl_BIM <- kforward_Bclxl_BIM * Bclxl * BIM - kbackward_Bclxl * Bclxl_BIM
    R15_deg_Mcl1_BIM     <- -kdeg_Mcl1_BIM * Mcl1_BIM
    R16_complex_Mcl1_BIM <- kforward_Mcl1_BIM * Mcl1 * BIM - kbackward_Mcl1 * Mcl1_BIM
    R17_deg_Bcl2_tBID     <- -kdeg_Bcl2_tBID * Bcl2_tBID
    R18_complex_Bcl2_tBID <- kforward_Bcl2_tBID * Bcl2 * tBID - kbackward_Bcl2 * Bcl2_tBID
    R19_deg_Bclxl_tBID     <- -kdeg_Bclxl_tBID * Bclxl_tBID
    R20_complex_Bclxl_tBID <- kforward_Bclxl_tBID * Bclxl * tBID - kbackward_Bclxl * Bclxl_tBID
    R21_deg_Mcl1_tBID     <- -kdeg_Mcl1_tBID * Mcl1_tBID
    R22_complex_Mcl1_tBID <- kforward_Mcl1_tBID * Mcl1 * tBID - kbackward_Mcl1 * Mcl1_tBID
    R23_deg_Bcl2_PUMA     <- -kdeg_Bcl2_PUMA * Bcl2_PUMA
    R24_complex_Bcl2_PUMA <- kforward_Bcl2_PUMA * Bcl2 * PUMA - kbackward_Bcl2 * Bcl2_PUMA
    R25_deg_Bclxl_PUMA     <- -kdeg_Bclxl_PUMA * Bclxl_PUMA
    R26_complex_Bclxl_PUMA <- kforward_Bclxl_PUMA * Bclxl * PUMA - kbackward_Bclxl * Bclxl_PUMA
    R27_deg_Mcl1_PUMA     <- -kdeg_Mcl1_PUMA * Mcl1_PUMA
    R28_complex_Mcl1_PUMA <- kforward_Mcl1_PUMA * Mcl1 * PUMA - kbackward_Mcl1 * Mcl1_PUMA
    R29_deg_Bcl2_NOXA     <- -kdeg_Bcl2_NOXA * Bcl2_NOXA
    R30_complex_Bcl2_NOXA <- kforward_Bcl2_NOXA * Bcl2 * NOXA - kbackward_Bcl2 * Bcl2_NOXA
    R31_deg_Bclxl_NOXA     <- -kdeg_Bclxl_NOXA * Bclxl_NOXA
    R32_complex_Bclxl_NOXA <- kforward_Bclxl_NOXA * Bclxl * NOXA - kbackward_Bclxl * Bclxl_NOXA
    R33_deg_Mcl1_NOXA     <- -kdeg_Mcl1_NOXA * Mcl1_NOXA
    R34_complex_Mcl1_NOXA <- kforward_Mcl1_NOXA * Mcl1 * NOXA - kbackward_Mcl1 * Mcl1_NOXA

    R35_complex_Bcl2_BAKa  <- kforward_Bcl2_BAKa * Bcl2 * BAKa - kbackward_Bcl2 * Bcl2_BAKa
    R36_complex_Bcl2_BAXma <- kforward_Bcl2_BAXma * Bcl2 * BAXma - kbackward_Bcl2 * Bcl2_BAXma
    R37_complex_Bclxl_BAKa  <- kforward_Bclxl_BAKa * Bclxl * BAKa - kbackward_Bclxl * Bclxl_BAKa
    R38_complex_Bclxl_BAXma <- kforward_Bclxl_BAXma * Bclxl * BAXma - kbackward_Bclxl * Bclxl_BAXma
    R39_complex_Mcl1_BAKa  <- kforward_Mcl1_BAKa * Mcl1 * BAKa - kbackward_Mcl1 * Mcl1_BAKa
    R40_complex_Mcl1_BAXma <- kforward_Mcl1_BAXma * Mcl1 * BAXma - kbackward_Mcl1 * Mcl1_BAXma

    R41_complex_BIM_BAXc  <- kforward_act * BIM * BAXc - kbackward_act * BIM_BAXc
    R42_complex_tBID_BAXc <- kforward_act * tBID * BAXc - kbackward_act * tBID_BAXc
    R43_complex_PUMA_BAXc <- kforward_act * PUMA * BAXc - kbackward_act * PUMA_BAXc
    R44_complex_BIM_BAK  <- kforward_act * BIM * BAK - kbackward_act * BIM_BAK
    R45_complex_tBID_BAK <- kforward_act * tBID * BAK - kbackward_act * tBID_BAK
    R46_complex_PUMA_BAK <- kforward_act * PUMA * BAK - kbackward_act * PUMA_BAK

    R47_disso_BIM_BAK  <- -k_activation * BIM_BAK
    R48_disso_tBID_BAK <- -k_activation * tBID_BAK
    R49_disso_PUMA_BAK <- -k_activation * PUMA_BAK
    R50_disso_BIM_BAXc  <- -k_activation * BIM_BAXc
    R51_disso_tBID_BAXc <- -k_activation * tBID_BAXc
    R52_disso_PUMA_BAXc <- -k_activation * PUMA_BAXc
    R53_disso_BAXca <- -k_activation * BAXca

    R54_complex_BAK_VDAC2 <- kforward_BAK_VDAC2 * BAK * VDAC2 - kbackward_BAK_VDAC2 * BAK_VDAC2

    R55_complex_BAK2  <- kforward_olig * BAKa * BAKa - kbackward_olig * BAK2
    R56_complex_BAK4  <- kforward_olig * BAK2 * BAK2 - kbackward_olig * BAK4
    R57_complex_BAK6  <- kforward_olig * BAK2 * BAK4 - kbackward_olig * BAK6
    R58_complex_BAK8  <- kforward_olig * BAK2 * BAK6 - kbackward_olig * BAK8
    R59_complex_BAK10 <- kforward_olig * BAK2 * BAK8 - kbackward_olig * BAK10
    R60_complex_BAK12 <- kforward_olig * BAK2 * BAK10 - kbackward_olig * BAK12
    R61_complex_BAK8  <- kforward_olig * BAK4 * BAK4 - kbackward_olig * BAK8
    R62_complex_BAK10 <- kforward_olig * BAK4 * BAK6 - kbackward_olig * BAK10
    R63_complex_BAK12 <- kforward_olig * BAK4 * BAK8 - kbackward_olig * BAK12
    R64_complex_BAK12 <- kforward_olig * BAK6 * BAK6 - kbackward_olig * BAK12
    R65_complex_BAX2  <- kforward_olig * BAXma * BAXma - kbackward_olig * BAX2
    R66_complex_BAX4  <- kforward_olig * BAX2 * BAX2 - kbackward_olig * BAX4
    R67_complex_BAX6  <- kforward_olig * BAX2 * BAX4 - kbackward_olig * BAX6
    R68_complex_BAX8  <- kforward_olig * BAX2 * BAX6 - kbackward_olig * BAX8
    R69_complex_BAX10 <- kforward_olig * BAX2 * BAX8 - kbackward_olig * BAX10
    R70_complex_BAX12 <- kforward_olig * BAX2 * BAX10 - kbackward_olig * BAX12
    R71_complex_BAX8  <- kforward_olig * BAX4 * BAX4 - kbackward_olig * BAX8
    R72_complex_BAX10 <- kforward_olig * BAX4 * BAX6 - kbackward_olig * BAX10
    R73_complex_BAX12 <- kforward_olig * BAX4 * BAX8 - kbackward_olig * BAX12
    R74_complex_BAX12 <- kforward_olig * BAX6 * BAX6 - kbackward_olig * BAX12

    # -------------------------------------------------------------------
    # 4. Drug switch and the in vivo mouse venetoclax PK. Drugs enter the
    #    system only after the burn-in; before that the system relaxes to
    #    its homeostatic equilibrium.
    # -------------------------------------------------------------------
    drug_act <- (t >= tburn)
    Ke_Veneto <- Cl_Veneto / Vd_Veneto
    d/dt(Veneto_gut)     <- -ka_Veneto * Veneto_gut * drug_act
    d/dt(Veneto_central) <- (ka_Veneto * Veneto_gut - Ke_Veneto * Veneto_central) * drug_act
    Veneto_plasma <- Veneto_central / Vd_Veneto
    # Deposited code multiplies the mg/L plasma value (not the uM value it also
    # computes) by ratioTumor before feeding the 1/uM/h kill term -- see vignette Errata.
    Veneto_tumor <- Veneto_plasma * ratioTumor
    # A-1592668 K-PD: amount decays at ke_Mcl1_I (0 by default, so a dose persists)
    d/dt(A15) <- -A15 * ke_Mcl1_I

    # In vitro drug states: constant while ke_*_I = 0
    d/dt(Bcl2_I)  <- -ke_BCl2_I * Bcl2_I * drug_act
    d/dt(Bclxl_I) <- -ke_BClxl_I * Bclxl_I * drug_act
    d/dt(Mcl1_I)  <- -ke_Mcl1_I * Mcl1_I * drug_act

    # -------------------------------------------------------------------
    # 5. Antiapoptotic proteins. BCL2 is attacked by both the in vitro
    #    venetoclax concentration (Bcl2_I) and the in vivo tumor exposure
    #    (Veneto_tumor); only one of the two is non-zero in practice.
    # -------------------------------------------------------------------
    d/dt(Bcl2) <- -Veneto_tumor * k2_Bcl2_I * Bcl2 * drug_act + R1_deg_Bcl2 + R2_prod_Bcl2 -
      R12_complex_Bcl2_BIM - R18_complex_Bcl2_tBID - R24_complex_Bcl2_PUMA - R30_complex_Bcl2_NOXA -
      R35_complex_Bcl2_BAKa - R36_complex_Bcl2_BAXma - Bcl2_I * k2_Bcl2_I * Bcl2 * drug_act
    d/dt(Bclxl) <- R3_deg_Bclxl + R4_prod_Bclxl - R14_complex_Bclxl_BIM - R20_complex_Bclxl_tBID -
      R26_complex_Bclxl_PUMA - R32_complex_Bclxl_NOXA - R37_complex_Bclxl_BAKa -
      R38_complex_Bclxl_BAXma - Bclxl_I * k2_Bclxl_I * Bclxl * drug_act
    d/dt(Mcl1) <- R5_deg_Mcl1 + R6_prod_Mcl1 - R16_complex_Mcl1_BIM - R22_complex_Mcl1_tBID -
      R28_complex_Mcl1_PUMA - R34_complex_Mcl1_NOXA - R39_complex_Mcl1_BAKa -
      R40_complex_Mcl1_BAXma - Mcl1_I * k2_Mcl1_I * Mcl1 * drug_act

    # -------------------------------------------------------------------
    # 6. BH3-only proteins. BIM, PUMA and NOXA carry the zero-order
    #    endogenous production that replaces Lindner's cytotoxic "stress
    #    dose"; tBID has none (it is the cytotoxic-stress species).
    # -------------------------------------------------------------------
    d/dt(BIM) <- R7_deg_BIM - R12_complex_Bcl2_BIM - R14_complex_Bclxl_BIM - R16_complex_Mcl1_BIM -
      R41_complex_BIM_BAXc - R44_complex_BIM_BAK - R47_disso_BIM_BAK - R50_disso_BIM_BAXc +
      BIM0 * kdeg_BIM * switchE
    d/dt(tBID) <- R8_deg_tBID - R18_complex_Bcl2_tBID - R20_complex_Bclxl_tBID - R22_complex_Mcl1_tBID -
      R42_complex_tBID_BAXc - R45_complex_tBID_BAK - R48_disso_tBID_BAK - R51_disso_tBID_BAXc
    d/dt(PUMA) <- R9_deg_PUMA - R24_complex_Bcl2_PUMA - R26_complex_Bclxl_PUMA - R28_complex_Mcl1_PUMA -
      R43_complex_PUMA_BAXc - R46_complex_PUMA_BAK - R49_disso_PUMA_BAK - R52_disso_PUMA_BAXc +
      PUMA0 * kdeg_PUMA * switchE
    d/dt(NOXA) <- R10_deg_NOXA - R30_complex_Bcl2_NOXA - R32_complex_Bclxl_NOXA - R34_complex_Mcl1_NOXA +
      NOXA0 * kdeg_NOXA * switchE

    # -------------------------------------------------------------------
    # 7. Sequestration complexes
    # -------------------------------------------------------------------
    d/dt(Bcl2_BIM)   <- R11_deg_Bcl2_BIM + R12_complex_Bcl2_BIM
    d/dt(Bclxl_BIM)  <- R13_deg_Bclxl_BIM + R14_complex_Bclxl_BIM
    d/dt(Mcl1_BIM)   <- R15_deg_Mcl1_BIM + R16_complex_Mcl1_BIM
    d/dt(Bcl2_tBID)  <- R17_deg_Bcl2_tBID + R18_complex_Bcl2_tBID
    d/dt(Bclxl_tBID) <- R19_deg_Bclxl_tBID + R20_complex_Bclxl_tBID
    d/dt(Mcl1_tBID)  <- R21_deg_Mcl1_tBID + R22_complex_Mcl1_tBID
    d/dt(Bcl2_PUMA)  <- R23_deg_Bcl2_PUMA + R24_complex_Bcl2_PUMA
    d/dt(Bclxl_PUMA) <- R25_deg_Bclxl_PUMA + R26_complex_Bclxl_PUMA
    d/dt(Mcl1_PUMA)  <- R27_deg_Mcl1_PUMA + R28_complex_Mcl1_PUMA
    d/dt(Bcl2_NOXA)  <- R29_deg_Bcl2_NOXA + R30_complex_Bcl2_NOXA
    d/dt(Bclxl_NOXA) <- R31_deg_Bclxl_NOXA + R32_complex_Bclxl_NOXA
    d/dt(Mcl1_NOXA)  <- R33_deg_Mcl1_NOXA + R34_complex_Mcl1_NOXA

    # -------------------------------------------------------------------
    # 8. BAX / BAK activation, capping and turnover. Every BAX/BAK-bearing
    #    species is degraded at kelimBAXBAK; BAXc, BAK and VDAC2 carry the
    #    matching zero-order production so total BAX and BAK are stationary.
    #
    #    NOTE: d/dt(Bclxl_BAKa) reproduces the deposited implementation
    #    verbatim, including the six BAXma flux terms it carries. Those
    #    terms have no mechanistic basis in the BCL-XL:BAK complex balance
    #    and are almost certainly a copy-paste artefact in the authors'
    #    code, but they are what generated every published result, so they
    #    are retained here. See the vignette Errata for the corrected form
    #    and a quantification of the difference.
    # -------------------------------------------------------------------
    d/dt(Bcl2_BAKa)  <- R35_complex_Bcl2_BAKa - Bcl2_BAKa * kelimBAXBAK * switchE
    d/dt(BAKa) <- -R35_complex_Bcl2_BAKa - R37_complex_Bclxl_BAKa - R39_complex_Mcl1_BAKa -
      R47_disso_BIM_BAK - R48_disso_tBID_BAK - R49_disso_PUMA_BAK - R55_complex_BAK2 -
      R55_complex_BAK2 - BAKa * kelimBAXBAK * switchE
    d/dt(Bcl2_BAXma) <- R36_complex_Bcl2_BAXma - Bcl2_BAXma * kelimBAXBAK * switchE
    d/dt(BAXma) <- -R36_complex_Bcl2_BAXma - R38_complex_Bclxl_BAXma - R40_complex_Mcl1_BAXma -
      R53_disso_BAXca - R65_complex_BAX2 - R65_complex_BAX2 - BAXma * kelimBAXBAK * switchE
    d/dt(Bclxl_BAKa) <- R37_complex_Bclxl_BAKa - R36_complex_Bcl2_BAXma - R38_complex_Bclxl_BAXma -
      R40_complex_Mcl1_BAXma - R53_disso_BAXca - R65_complex_BAX2 - R65_complex_BAX2 -
      Bclxl_BAKa * kelimBAXBAK * switchE
    d/dt(Bclxl_BAXma) <- R38_complex_Bclxl_BAXma - Bclxl_BAXma * kelimBAXBAK * switchE
    d/dt(Mcl1_BAKa)  <- R39_complex_Mcl1_BAKa - Mcl1_BAKa * kelimBAXBAK * switchE
    d/dt(Mcl1_BAXma) <- R40_complex_Mcl1_BAXma - Mcl1_BAXma * kelimBAXBAK * switchE
    d/dt(BIM_BAXc)  <- R41_complex_BIM_BAXc + R50_disso_BIM_BAXc - BIM_BAXc * kelimBAXBAK * switchE
    d/dt(BAXc) <- -R41_complex_BIM_BAXc - R42_complex_tBID_BAXc - R43_complex_PUMA_BAXc -
      BAXc * kelimBAXBAK * switchE + BAXc0 * kelimBAXBAK * switchE
    d/dt(tBID_BAXc) <- R42_complex_tBID_BAXc + R51_disso_tBID_BAXc - tBID_BAXc * kelimBAXBAK * switchE
    d/dt(PUMA_BAXc) <- R43_complex_PUMA_BAXc + R52_disso_PUMA_BAXc - PUMA_BAXc * kelimBAXBAK * switchE
    d/dt(BIM_BAK)  <- R44_complex_BIM_BAK + R47_disso_BIM_BAK - BIM_BAK * kelimBAXBAK * switchE
    d/dt(BAK) <- -R44_complex_BIM_BAK - R45_complex_tBID_BAK - R46_complex_PUMA_BAK -
      R54_complex_BAK_VDAC2 - BAK * kelimBAXBAK * switchE + BAK0 * kelimBAXBAK * switchE
    d/dt(tBID_BAK) <- R45_complex_tBID_BAK + R48_disso_tBID_BAK - tBID_BAK * kelimBAXBAK * switchE
    d/dt(PUMA_BAK) <- R46_complex_PUMA_BAK + R49_disso_PUMA_BAK - PUMA_BAK * kelimBAXBAK * switchE
    d/dt(BAXca) <- -R50_disso_BIM_BAXc - R51_disso_tBID_BAXc - R52_disso_PUMA_BAXc + R53_disso_BAXca -
      BAXca * kelimBAXBAK * switchE
    d/dt(BAK_VDAC2) <- R54_complex_BAK_VDAC2 - BAK_VDAC2 * kelimBAXBAK * switchE
    d/dt(VDAC2) <- -R54_complex_BAK_VDAC2 - VDAC2 * kelimBAXBAK * switchE + BAK0 * kelimBAXBAK * switchE

    # -------------------------------------------------------------------
    # 9. Oligomerisation cascades (dimer through dodecamer)
    # -------------------------------------------------------------------
    d/dt(BAK2) <- R55_complex_BAK2 - R56_complex_BAK4 - R56_complex_BAK4 - R57_complex_BAK6 -
      R58_complex_BAK8 - R59_complex_BAK10 - R60_complex_BAK12 - BAK2 * kelimBAXBAK * switchE
    d/dt(BAK4) <- R56_complex_BAK4 - R57_complex_BAK6 - R61_complex_BAK8 - R61_complex_BAK8 -
      R62_complex_BAK10 - R63_complex_BAK12 - BAK4 * kelimBAXBAK * switchE
    d/dt(BAK6) <- R57_complex_BAK6 - R58_complex_BAK8 - R62_complex_BAK10 - R64_complex_BAK12 -
      R64_complex_BAK12 - BAK6 * kelimBAXBAK * switchE
    d/dt(BAK8) <- R58_complex_BAK8 - R59_complex_BAK10 + R61_complex_BAK8 - R63_complex_BAK12 -
      BAK8 * kelimBAXBAK * switchE
    d/dt(BAK10) <- R59_complex_BAK10 - R60_complex_BAK12 + R62_complex_BAK10 -
      BAK10 * kelimBAXBAK * switchE
    d/dt(BAK12) <- R60_complex_BAK12 + R63_complex_BAK12 + R64_complex_BAK12 -
      BAK12 * kelimBAXBAK * switchE
    d/dt(BAX2) <- R65_complex_BAX2 - R66_complex_BAX4 - R66_complex_BAX4 - R67_complex_BAX6 -
      R68_complex_BAX8 - R69_complex_BAX10 - R70_complex_BAX12 - BAX2 * kelimBAXBAK * switchE
    d/dt(BAX4) <- R66_complex_BAX4 - R67_complex_BAX6 - R71_complex_BAX8 - R71_complex_BAX8 -
      R72_complex_BAX10 - R73_complex_BAX12 - BAX4 * kelimBAXBAK * switchE
    d/dt(BAX6) <- R67_complex_BAX6 - R68_complex_BAX8 - R72_complex_BAX10 - R74_complex_BAX12 -
      R74_complex_BAX12 - BAX6 * kelimBAXBAK * switchE
    d/dt(BAX8) <- R68_complex_BAX8 - R69_complex_BAX10 + R71_complex_BAX8 - R73_complex_BAX12 -
      BAX8 * kelimBAXBAK * switchE
    d/dt(BAX10) <- R69_complex_BAX10 - R70_complex_BAX12 + R72_complex_BAX10 -
      BAX10 * kelimBAXBAK * switchE
    d/dt(BAX12) <- R70_complex_BAX12 + R73_complex_BAX12 + R74_complex_BAX12 -
      BAX12 * kelimBAXBAK * switchE

    # -------------------------------------------------------------------
    # 10. MOMP readout. Pore is the percent of total BAX + BAK sitting in
    #     hexamers or larger. TimeAbove accumulates the hours spent above
    #     the threshold; the paper's cell-death criterion is
    #     max(TimeAbove) > 0.1 h between drug administration and 748 h.
    # -------------------------------------------------------------------
    pctBAX <- 100 * (6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12) / BAXc0
    pctBAK <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12) / BAK0
    # The denominator is written as ((BAK0) + (BAXc0)) rather than (BAK0 + BAXc0):
    # rxode2's mu-referencing check rejects two bare population parameters in one
    # expression, and the extra parentheses are the standard workaround. (Its error
    # message names the first two ini() entries, not the offending ones.)
    Pore <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12 +
                     6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12) / ((BAK0) + (BAXc0))
    d/dt(TimeAbove) <- (Pore > poreThreshold)
  })
}
