Dasti_2025_bnt162b2_qsp <- function() {
  description <- paste(
    "QSP (multiscale mRNA-vaccine immune response, tissue layer). Dasti 2025",
    "model of the humoral response to the Pfizer-BioNTech BNT162b2 COVID-19",
    "mRNA vaccine in the general adult population. Three physiological",
    "compartments (intramuscular injection site, draining lymph node, blood)",
    "track LNP-encapsulated mRNA uptake and degradation, four antigen-",
    "presenting-cell populations (neutrophils, monocytes, myeloid and",
    "plasmacytoid dendritic cells) with a reversible low/medium/high",
    "MHC-II-antigen exposure chain on each dendritic-cell population, naive /",
    "activated / memory / functional helper T cells, and 17 B-cell affinity",
    "subclones spanning a 2-fold Ka ladder that each activate, enter the",
    "germinal centre, and differentiate into memory B cells and short- and",
    "long-lived plasma cells. Plasma cells migrate to blood and secrete",
    "anti-RBD IgG, the single observable. 220 ODEs. The dose-dependent",
    "saturable mRNA-degradation term reproduces the non-linear dose /",
    "immunogenicity relationship. Deterministic mechanism model: the authors",
    "quantified uncertainty by Fisher-information sampling of the parameter",
    "covariance rather than by fitting IIV, so no etas are encoded and the",
    "residual SD is fixed at 0. The molecular layer that supplies the",
    "dendritic-cell transition rates is a separate model file",
    "(Dasti_2025_bnt162b2_molecular_qsp).",
    sep = " "
  )
  reference <- paste(
    "Dasti L, Giampiccolo S, Pettina E, Fiandaca G, Zangani N, Leonardelli L,",
    "De Lima Hedayioglu F, Campanile E, Marchetti L. A Multiscale Quantitative",
    "Systems Pharmacology Model for the Development and Optimization of mRNA",
    "Vaccines. CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1213-1224.",
    "doi:10.1002/psp4.70041.",
    "Tissue-layer state variables from Supporting Information Table S1,",
    "initial values from Table S2, parameters from Table S3 for the general adult population,",
    "equations from Supporting Information Section S1.4 (equations 1-66).",
    "Full-precision parameter values and the resolution of three Supporting",
    "Information typographical errors were taken from the authors' published",
    "MATLAB implementation at https://github.com/cosbi-research/QSPmRNAVaccines",
    "(Figure_reproduction/Figure 4/{model_equations,Parameters,simulation}.m);",
    "see the vignette Errata section.",
    sep = " "
  )
  vignette <- "Dasti_2025_mrna_vaccines"
  units <- list(
    time          = "day",
    dosing        = "pmol (LNP-encapsulated mRNA at the injection site)",
    concentration = "ng/mL (total anti-RBD IgG in serum)"
  )
  # The vaccine dose is given into the mRNA state, not depot/central.
  dosing <- "mRNA"
  # Every state is a paper-mechanistic immunological species; the j=1..17
  # suffixes index the B-cell affinity subclones of SI Table S1, not
  # arbitrary numbering.
  paper_specific_compartments <- c(
    "mRNA", "NP_IS", "NPL_IS", "NPAg_IS", "MN_IS", "MNL_IS",
    "MNAg_IS", "mDC_IS", "mDCL_IS", "mDCAgLon_IS", "mDCAgMon_IS", "mDCAgH_IS",
    "mDCAgMoff_IS", "mDCAgLoff_IS", "mDCoff_IS", "pDC_IS", "pDCL_IS", "pDCAgLon_IS",
    "pDCAgMon_IS", "pDCAgH_IS", "pDCAgMoff_IS", "pDCAgLoff_IS", "pDCoff_IS", "NPL_LN",
    "NPAg_LN", "MNL_LN", "MNAg_LN", "mDCL_LN", "mDCAgLon_LN", "mDCAgMon_LN",
    "mDCAgH_LN", "mDCAgMoff_LN", "mDCAgLoff_LN", "mDC_LN", "pDCL_LN", "pDCAgLon_LN",
    "pDCAgMon_LN", "pDCAgH_LN", "pDCAgMoff_LN", "pDCAgLoff_LN", "pDC_LN", "NT",
    "AT_N", "MT", "AT_M", "FT", "NB1", "NB2",
    "NB3", "NB4", "NB5", "NB6", "NB7", "NB8",
    "NB9", "NB10", "NB11", "NB12", "NB13", "NB14",
    "NB15", "NB16", "NB17", "ABN1", "ABN2", "ABN3",
    "ABN4", "ABN5", "ABN6", "ABN7", "ABN8", "ABN9",
    "ABN10", "ABN11", "ABN12", "ABN13", "ABN14", "ABN15",
    "ABN16", "ABN17", "GCB1", "GCB2", "GCB3", "GCB4",
    "GCB5", "GCB6", "GCB7", "GCB8", "GCB9", "GCB10",
    "GCB11", "GCB12", "GCB13", "GCB14", "GCB15", "GCB16",
    "GCB17", "MB1", "MB2", "MB3", "MB4", "MB5",
    "MB6", "MB7", "MB8", "MB9", "MB10", "MB11",
    "MB12", "MB13", "MB14", "MB15", "MB16", "MB17",
    "ABM1", "ABM2", "ABM3", "ABM4", "ABM5", "ABM6",
    "ABM7", "ABM8", "ABM9", "ABM10", "ABM11", "ABM12",
    "ABM13", "ABM14", "ABM15", "ABM16", "ABM17", "SP1",
    "SP2", "SP3", "SP4", "SP5", "SP6", "SP7",
    "SP8", "SP9", "SP10", "SP11", "SP12", "SP13",
    "SP14", "SP15", "SP16", "SP17", "LP1", "LP2",
    "LP3", "LP4", "LP5", "LP6", "LP7", "LP8",
    "LP9", "LP10", "LP11", "LP12", "LP13", "LP14",
    "LP15", "LP16", "LP17", "NP_BL", "MN_BL", "mDC_BL",
    "pDC_BL", "SPBL1", "SPBL2", "SPBL3", "SPBL4", "SPBL5",
    "SPBL6", "SPBL7", "SPBL8", "SPBL9", "SPBL10", "SPBL11",
    "SPBL12", "SPBL13", "SPBL14", "SPBL15", "SPBL16", "SPBL17",
    "LPBL1", "LPBL2", "LPBL3", "LPBL4", "LPBL5", "LPBL6",
    "LPBL7", "LPBL8", "LPBL9", "LPBL10", "LPBL11", "LPBL12",
    "LPBL13", "LPBL14", "LPBL15", "LPBL16", "LPBL17", "Ab1",
    "Ab2", "Ab3", "Ab4", "Ab5", "Ab6", "Ab7",
    "Ab8", "Ab9", "Ab10", "Ab11", "Ab12", "Ab13",
    "Ab14", "Ab15", "Ab16", "Ab17"
  )
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Descriptions transcribed from SI Table S1;
  # verified = TRUE means checked against that table. Sub-anatomical sites
  # (injection-site muscle, lymph node, blood) are mapped onto the
  # controlled specimen vocabulary; the analyte text carries the detail.
  compartmentData <- list(
    mRNA         = list(analyte = "LNP-encapsulated mRNA product at the intramuscular injection site", units = "pmol", specimen = "administration site", verified = TRUE),
    NP_IS        = list(analyte = "resident naive neutrophils at the injection site", units = "cells", specimen = "tissue", verified = TRUE),
    NPL_IS       = list(analyte = "neutrophils that internalized the LNPs", units = "cells", specimen = "tissue", verified = TRUE),
    NPAg_IS      = list(analyte = "neutrophils expressing the antigen protein", units = "cells", specimen = "tissue", verified = TRUE),
    MN_IS        = list(analyte = "resident naive monocytes at the injection site", units = "cells", specimen = "tissue", verified = TRUE),
    MNL_IS       = list(analyte = "monocytes that internalized the LNPs", units = "cells", specimen = "tissue", verified = TRUE),
    MNAg_IS      = list(analyte = "monocytes expressing the antigen protein", units = "cells", specimen = "tissue", verified = TRUE),
    mDC_IS       = list(analyte = "resident naive myeloid dendritic cells at the injection site", units = "cells", specimen = "tissue", verified = TRUE),
    mDCL_IS      = list(analyte = "myeloid dendritic cells that internalized the LNPs", units = "cells", specimen = "tissue", verified = TRUE),
    mDCAgLon_IS  = list(analyte = "myeloid dendritic cells with low antigen expression (binding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    mDCAgMon_IS  = list(analyte = "myeloid dendritic cells with medium antigen expression (binding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    mDCAgH_IS    = list(analyte = "myeloid dendritic cells with high antigen expression", units = "cells", specimen = "tissue", verified = TRUE),
    mDCAgMoff_IS = list(analyte = "myeloid dendritic cells with medium antigen expression (unbinding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    mDCAgLoff_IS = list(analyte = "myeloid dendritic cells with low antigen expression (unbinding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    mDCoff_IS    = list(analyte = "mature myeloid dendritic cells with no antigen expression (unbinding-phase endpoint)", units = "cells", specimen = "tissue", verified = TRUE),
    pDC_IS       = list(analyte = "resident naive plasmacytoid dendritic cells at the injection site", units = "cells", specimen = "tissue", verified = TRUE),
    pDCL_IS      = list(analyte = "plasmacytoid dendritic cells that internalized the LNPs", units = "cells", specimen = "tissue", verified = TRUE),
    pDCAgLon_IS  = list(analyte = "plasmacytoid dendritic cells with low antigen expression (binding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    pDCAgMon_IS  = list(analyte = "plasmacytoid dendritic cells with medium antigen expression (binding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    pDCAgH_IS    = list(analyte = "plasmacytoid dendritic cells with high antigen expression", units = "cells", specimen = "tissue", verified = TRUE),
    pDCAgMoff_IS = list(analyte = "plasmacytoid dendritic cells with medium antigen expression (unbinding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    pDCAgLoff_IS = list(analyte = "plasmacytoid dendritic cells with low antigen expression (unbinding phase)", units = "cells", specimen = "tissue", verified = TRUE),
    pDCoff_IS    = list(analyte = "mature plasmacytoid dendritic cells with no antigen expression (unbinding-phase endpoint)", units = "cells", specimen = "tissue", verified = TRUE),
    NPL_LN       = list(analyte = "neutrophils that internalized the LNPs", units = "cells", specimen = "lymph", verified = TRUE),
    NPAg_LN      = list(analyte = "neutrophils producing antigen protein from translated mRNA", units = "cells", specimen = "lymph", verified = TRUE),
    MNL_LN       = list(analyte = "monocytes that internalized the LNPs", units = "cells", specimen = "lymph", verified = TRUE),
    MNAg_LN      = list(analyte = "monocytes producing antigen protein from translated mRNA", units = "cells", specimen = "lymph", verified = TRUE),
    mDCL_LN      = list(analyte = "myeloid dendritic cells that internalized the LNPs", units = "cells", specimen = "lymph", verified = TRUE),
    mDCAgLon_LN  = list(analyte = "myeloid dendritic cells with low antigen expression (binding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    mDCAgMon_LN  = list(analyte = "myeloid dendritic cells with medium antigen expression (binding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    mDCAgH_LN    = list(analyte = "myeloid dendritic cells with high antigen expression", units = "cells", specimen = "lymph", verified = TRUE),
    mDCAgMoff_LN = list(analyte = "myeloid dendritic cells with medium antigen expression (unbinding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    mDCAgLoff_LN = list(analyte = "myeloid dendritic cells with low antigen expression (unbinding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    mDC_LN       = list(analyte = "mature myeloid dendritic cells with no antigen expression (unbinding-phase endpoint)", units = "cells", specimen = "lymph", verified = TRUE),
    pDCL_LN      = list(analyte = "plasmacytoid dendritic cells that internalized the LNPs", units = "cells", specimen = "lymph", verified = TRUE),
    pDCAgLon_LN  = list(analyte = "plasmacytoid dendritic cells with low antigen expression (binding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    pDCAgMon_LN  = list(analyte = "plasmacytoid dendritic cells with medium antigen expression (binding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    pDCAgH_LN    = list(analyte = "plasmacytoid dendritic cells with high antigen expression", units = "cells", specimen = "lymph", verified = TRUE),
    pDCAgMoff_LN = list(analyte = "plasmacytoid dendritic cells with medium antigen expression (unbinding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    pDCAgLoff_LN = list(analyte = "plasmacytoid dendritic cells with low antigen expression (unbinding phase)", units = "cells", specimen = "lymph", verified = TRUE),
    pDC_LN       = list(analyte = "mature plasmacytoid dendritic cells with no antigen expression (unbinding-phase endpoint)", units = "cells", specimen = "lymph", verified = TRUE),
    NT           = list(analyte = "naive T-helper cells", units = "cells", specimen = "lymph", verified = TRUE),
    AT_N         = list(analyte = "activated T-helper cells derived from naive T cells", units = "cells", specimen = "lymph", verified = TRUE),
    MT           = list(analyte = "memory T cells", units = "cells", specimen = "lymph", verified = TRUE),
    AT_M         = list(analyte = "activated T-helper cells derived from memory T cells", units = "cells", specimen = "lymph", verified = TRUE),
    FT           = list(analyte = "functional T cells", units = "cells", specimen = "lymph", verified = TRUE),
    NB1          = list(analyte = "naive B cells, affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB2          = list(analyte = "naive B cells, affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB3          = list(analyte = "naive B cells, affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB4          = list(analyte = "naive B cells, affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB5          = list(analyte = "naive B cells, affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB6          = list(analyte = "naive B cells, affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB7          = list(analyte = "naive B cells, affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB8          = list(analyte = "naive B cells, affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB9          = list(analyte = "naive B cells, affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB10         = list(analyte = "naive B cells, affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB11         = list(analyte = "naive B cells, affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB12         = list(analyte = "naive B cells, affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB13         = list(analyte = "naive B cells, affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB14         = list(analyte = "naive B cells, affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB15         = list(analyte = "naive B cells, affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB16         = list(analyte = "naive B cells, affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NB17         = list(analyte = "naive B cells, affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN1         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN2         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN3         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN4         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN5         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN6         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN7         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN8         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN9         = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN10        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN11        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN12        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN13        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN14        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN15        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN16        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABN17        = list(analyte = "activated B cells derived from naive B cells (SI Table S1 ANB), affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB1         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB2         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB3         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB4         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB5         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB6         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB7         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB8         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB9         = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB10        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB11        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB12        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB13        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB14        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB15        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB16        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    GCB17        = list(analyte = "activated B cells that entered the germinal centre, affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB1          = list(analyte = "memory B cells, affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB2          = list(analyte = "memory B cells, affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB3          = list(analyte = "memory B cells, affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB4          = list(analyte = "memory B cells, affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB5          = list(analyte = "memory B cells, affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB6          = list(analyte = "memory B cells, affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB7          = list(analyte = "memory B cells, affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB8          = list(analyte = "memory B cells, affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB9          = list(analyte = "memory B cells, affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB10         = list(analyte = "memory B cells, affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB11         = list(analyte = "memory B cells, affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB12         = list(analyte = "memory B cells, affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB13         = list(analyte = "memory B cells, affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB14         = list(analyte = "memory B cells, affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB15         = list(analyte = "memory B cells, affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB16         = list(analyte = "memory B cells, affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    MB17         = list(analyte = "memory B cells, affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM1         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM2         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM3         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM4         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM5         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM6         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM7         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM8         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM9         = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM10        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM11        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM12        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM13        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM14        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM15        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM16        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    ABM17        = list(analyte = "activated B cells derived from memory B cells (SI Table S1 AMB), affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP1          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP2          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP3          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP4          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP5          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP6          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP7          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP8          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP9          = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP10         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP11         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP12         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP13         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP14         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP15         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP16         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    SP17         = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP1          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 1 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP2          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 2 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP3          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 3 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP4          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 4 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP5          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 5 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP6          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 6 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP7          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 7 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP8          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 8 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP9          = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 9 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP10         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 10 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP11         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 11 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP12         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 12 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP13         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 13 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP14         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 14 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP15         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 15 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP16         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 16 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    LP17         = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 17 of 17", units = "cells", specimen = "lymph", verified = TRUE),
    NP_BL        = list(analyte = "naive neutrophils", units = "cells", specimen = "blood cell", verified = TRUE),
    MN_BL        = list(analyte = "naive monocytes", units = "cells", specimen = "blood cell", verified = TRUE),
    mDC_BL       = list(analyte = "naive myeloid dendritic cells", units = "cells", specimen = "blood cell", verified = TRUE),
    pDC_BL       = list(analyte = "naive plasmacytoid dendritic cells", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL1        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 1 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL2        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 2 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL3        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 3 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL4        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 4 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL5        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 5 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL6        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 6 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL7        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 7 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL8        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 8 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL9        = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 9 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL10       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 10 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL11       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 11 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL12       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 12 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL13       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 13 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL14       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 14 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL15       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 15 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL16       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 16 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    SPBL17       = list(analyte = "short-lived antibody-secreting plasma cells, affinity subclone 17 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL1        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 1 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL2        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 2 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL3        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 3 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL4        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 4 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL5        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 5 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL6        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 6 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL7        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 7 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL8        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 8 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL9        = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 9 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL10       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 10 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL11       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 11 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL12       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 12 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL13       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 13 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL14       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 14 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL15       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 15 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL16       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 16 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    LPBL17       = list(analyte = "long-lived antibody-secreting plasma cells, affinity subclone 17 of 17", units = "cells", specimen = "blood cell", verified = TRUE),
    Ab1          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 1 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab2          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 2 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab3          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 3 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab4          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 4 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab5          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 5 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab6          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 6 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab7          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 7 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab8          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 8 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab9          = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 9 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab10         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 10 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab11         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 11 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab12         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 12 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab13         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 13 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab14         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 14 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab15         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 15 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab16         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 16 of 17", units = "pmol", specimen = "serum", verified = TRUE),
    Ab17         = list(analyte = "anti-spike (anti-RBD) IgG antibody in blood, affinity subclone 17 of 17", units = "pmol", specimen = "serum", verified = TRUE)
  )
  covariateData <- list()
  population <- list(
    species       = "human",
    disease_state = "healthy adults receiving prophylactic COVID-19 mRNA vaccination",
    age_range     = "general adult population (not age-stratified; a separate over-60 calibration is provided in Dasti_2025_bnt162b2_over60_qsp)",
    vaccine       = "Pfizer-BioNTech BNT162b2",
    dose_range    = "1, 10, 20 and 30 ug mRNA per injection; 21-day and 70-day prime-boost intervals and a third dose at 9 months were simulated",
    n_subjects    = NA_integer_,
    notes         = paste(
      "Calibrated against anti-RBD IgG serum concentrations from Sahin 2021",
      "(BNT162b2 phase 1/2, 1 and 30 ug), Keshavarz 2022 and Goel 2021, and",
      "validated against Sahin 2021 (10 and 20 ug), Kontopoulou 2022,",
      "Naaber 2021, Takeuchi 2022 and Payne 2021 (PITCH, extended interval).",
      "Antigen-presenting-cell parameters were calibrated on rhesus-macaque",
      "data from Liang 2017 (an anti-H10N8 influenza mRNA vaccine), so the",
      "early injection-site and lymph-node dynamics are non-human-primate",
      "derived while the T-cell, B-cell and antibody parameters are human.",
      "Sample sizes of the source studies are listed in SI Table S8.",
      sep = " "
    )
  )

  ini({
    kdmRNA       <- fixed(log(2) / (10 / 24)); label("mRNA degradation rate (1/day; from the 10-h mRNA half-life)") # Table S3 general parameters (1.6636/day); ref 3
    VLN          <- fixed(5e-4); label("Draining lymph node volume (L)") # Table S3 general parameters; ref 5
    VBL          <- fixed(5); label("Blood volume (L)") # Table S2 blood initial values (all four are per-5-L); Table S3 prints 5.5 L -- see vignette Errata
    NAvo         <- fixed(6.022e23); label("Avogadro constant (molecules/mole)") # Table S3 general parameters; ref 7
    kdeg         <- fixed(0.118692187); label("LNP/mRNA degradation rate by muscle uptake (1/day)") # Table S3 (0.1186); fitted on Liang 2017 data
    KmRNA        <- fixed(0.2072438573); label("mRNA amount giving 50% naive-APC recruitment rate (pmol)") # Table S3 (0.2072); fitted on Liang 2017 data
    SFmRNA       <- fixed(2958.763891); label("Scaling factor for the mRNA uptake rate (1/pmol)") # Table S3 (2.9587e+03); fitted on Liang 2017 data
    mRNAmax      <- fixed(3.141600294); label("mRNA amount triggering the regulatory response at the injection site (pmol)") # Table S3 (3.1416); fitted on antibody data
    ksat         <- fixed(55.83299964); label("mRNA degradation rate from the regulatory mechanism (1/day)") # Table S3 (55.8330); fitted on antibody data
    kslope       <- fixed(1); label("Slope of the regulatory (tanh) switch (pmol)") # Table S3; model assumption to adimensionalise the term
    expAg        <- fixed(337132.5865); label("Scaling factor for antigen quantity in the lymph node (unitless)") # Table S3 (3.3713e+05); fitted on antibody data, held constant across all three calibrations (SI Section S5)
    kdtNP        <- fixed(24 * log(2) / 7); label("Neutrophil death rate (1/day; 7-h half-life)") # Table S3 (2.3765/day); ref 4
    krcNP        <- fixed(0.00086040978); label("Neutrophil recruitment rate, mRNA-driven (1/day)") # Table S3 (8.6040e-04); fitted on Liang 2017 data
    kupNP        <- fixed(2.321695675e-06); label("Neutrophil LNP internalisation rate (1/day)") # Table S3 (2.3216e-06); fitted on Liang 2017 data
    kexpNP       <- fixed(0.1822345729); label("Neutrophil antigen expression rate (1/day)") # Table S3 (0.1822); fitted on Liang 2017 data
    kis2lnNPlnp  <- fixed(0.02263828077); label("LNP-loaded neutrophil migration rate, injection site to lymph node (1/day)") # published code Parameters.m Liang_fit(8); Table S3 misprints this as 0.1822 (a duplicate of kexpNP) -- see vignette Errata
    kis2lnNPag   <- fixed(0.1112715793); label("Antigen-expressing neutrophil migration rate, injection site to lymph node (1/day)") # Table S3 (0.1112); fitted on Liang 2017 data
    kis2blNP     <- fixed(0.0004102812887); label("Naive neutrophil migration rate, injection site to blood (1/day)") # Table S3 (4.1028e-04); fitted on Liang 2017 data
    NP0IS        <- fixed(0); label("Initial naive neutrophils at the injection site (cells)") # Table S2 (0 cells); Liang 2017 first timepoint
    NP0BL        <- fixed(1.595e10); label("Initial naive neutrophils in blood (cells)") # Table S2 (1.5950e+10 cells); ref 4
    kdtMN        <- fixed(24 * log(2) / 24); label("Monocyte death rate (1/day; 24-h half-life)") # Table S3 (0.6931/day); ref 4
    krcMN        <- fixed(0.003809771376); label("Monocyte recruitment rate, mRNA-driven (1/day)") # Table S3 (0.0038); fitted on Liang 2017 data
    kupMN        <- fixed(4.059326032e-06); label("Monocyte LNP internalisation rate (1/day)") # Table S3 (4.0593e-06); fitted on Liang 2017 data
    kexpMN       <- fixed(4.30403771); label("Monocyte antigen expression rate (1/day)") # Table S3 (4.3040); fitted on Liang 2017 data
    kis2lnMNlnp  <- fixed(0.5219812732); label("LNP-loaded monocyte migration rate, injection site to lymph node (1/day)") # Table S3 (0.5219); fitted on Liang 2017 data
    kis2lnMNag   <- fixed(0.0001194687646); label("Antigen-expressing monocyte migration rate, injection site to lymph node (1/day)") # Table S3 (1.1946e-04); fitted on Liang 2017 data
    kis2blMN     <- fixed(0.006083687683); label("Naive monocyte migration rate, injection site to blood (1/day)") # Table S3 (0.0060); fitted on Liang 2017 data
    MN0IS        <- fixed(300); label("Initial naive monocytes at the injection site (cells)") # Table S2 (300 cells); Liang 2017 first timepoint
    MN0BL        <- fixed(1.88e9); label("Initial naive monocytes in blood (cells)") # Table S2 (1.8800e+09 cells); ref 5
    kdtImDC      <- fixed(0.0924); label("Naive (immature) myeloid dendritic cell death rate (1/day)") # Table S3 (0.0924); ref 4
    kdtmDC       <- fixed(1.140252427); label("Mature myeloid dendritic cell death rate (1/day)") # Table S3 (1.1402); fitted on Liang 2017 data
    krcmDC       <- fixed(0.005220380831); label("Myeloid dendritic cell recruitment rate, mRNA-driven (1/day)") # Table S3 (0.0052); fitted on Liang 2017 data
    kupmDC       <- fixed(5.085514757e-06); label("Myeloid dendritic cell LNP internalisation rate (1/day)") # Table S3 (5.0855e-06); fitted on Liang 2017 data
    kexpmDC      <- fixed(4.9219018); label("Myeloid dendritic cell maturation rate, mDC_LNP to mDC_Ag_Lon (1/day)") # Table S3 (4.9219); fitted on Liang 2017 data
    kis2lnmDClnp <- fixed(0.1468493298); label("LNP-loaded myeloid dendritic cell migration rate, injection site to lymph node (1/day)") # Table S3 (0.1468); fitted on Liang 2017 data
    kis2lnmDCag  <- fixed(0.2351721546); label("Antigen-presenting myeloid dendritic cell migration rate, injection site to lymph node (1/day)") # Table S3 (0.2351); fitted on Liang 2017 data
    kis2blmDC    <- fixed(0.2449181227); label("Naive myeloid dendritic cell migration rate, injection site to blood (1/day)") # Table S3 (0.2449); fitted on Liang 2017 data
    mDC0IS       <- fixed(1864.37); label("Initial naive myeloid dendritic cells at the injection site (cells)") # Table S2 (1.8644e+03 cells); Liang 2017 first timepoint
    mDC0BL       <- fixed(6.1925e7); label("Initial naive myeloid dendritic cells in blood (cells)") # Table S2 (6.1925e+07 cells); ref 5
    kdtIpDC      <- fixed(0.0924); label("Naive (immature) plasmacytoid dendritic cell death rate (1/day)") # Table S3 (0.0924); ref 4
    kdtpDC       <- fixed(1.350173307); label("Mature plasmacytoid dendritic cell death rate (1/day)") # Table S3 (1.3501); fitted on Liang 2017 data
    krcpDC       <- fixed(0.0004999447019); label("Plasmacytoid dendritic cell recruitment rate, mRNA-driven (1/day)") # Table S3 (4.9994e-04); fitted on Liang 2017 data
    kuppDC       <- fixed(4.430919993e-06); label("Plasmacytoid dendritic cell LNP internalisation rate (1/day)") # Table S3 (4.4309e-06); fitted on Liang 2017 data
    kexppDC      <- fixed(29.15987199); label("Plasmacytoid dendritic cell maturation rate, pDC_LNP to pDC_Ag_Lon (1/day)") # Table S3 (29.1598); fitted on Liang 2017 data
    kis2lnpDClnp <- fixed(27.07264527); label("LNP-loaded plasmacytoid dendritic cell migration rate, injection site to lymph node (1/day)") # Table S3 (27.0726); fitted on Liang 2017 data
    kis2lnpDCag  <- fixed(0.003917403408); label("Antigen-presenting plasmacytoid dendritic cell migration rate, injection site to lymph node (1/day)") # Table S3 (0.0039); fitted on Liang 2017 data
    kis2blpDC    <- fixed(0.1548673895); label("Naive plasmacytoid dendritic cell migration rate, injection site to blood (1/day)") # Table S3 (0.1548); fitted on Liang 2017 data
    pDC0IS       <- fixed(0); label("Initial naive plasmacytoid dendritic cells at the injection site (cells)") # Table S2 (0 cells); Liang 2017 first timepoint
    pDC0BL       <- fixed(3.5e7); label("Initial naive plasmacytoid dendritic cells in blood (cells)") # Table S2 (3.5e+07 cells); ref 5
    pL           <- fixed(0.2); label("Low antigen-expression maturation weight (unitless)") # Table S3 (0.2); model assumption
    pM           <- fixed(0.65); label("Medium antigen-expression maturation weight (unitless)") # Table S3 (0.65); model assumption
    pH           <- fixed(1); label("High antigen-expression maturation weight (unitless)") # Table S3 (1); model assumption
    ktrMmDC      <- fixed(5.984386009847674); label("mDC transition rate, Ag_Lon to Ag_Mon (1/day)") # Table S3 (5.9867); derived from the molecular layer, SI Section S3 -- full precision from published code Parameters.m
    ktrHmDC      <- fixed(4.023456291503178); label("mDC transition rate, Ag_Mon to Ag_H (1/day)") # Table S3 (4.0242); derived from the molecular layer, SI Section S3
    katrHmDC     <- fixed(2.00114454171865); label("mDC transition rate, Ag_H to Ag_Moff (1/day)") # Table S3 (2.0006); derived from the molecular layer, SI Section S3
    katrMmDC     <- fixed(1.18252548339979); label("mDC transition rate, Ag_Moff to Ag_Loff (1/day)") # Table S3 (1.1825); derived from the molecular layer, SI Section S3
    katrLmDC     <- fixed(0.393936426260313); label("mDC transition rate, Ag_Loff to mDC_off (1/day)") # Table S3 (0.3939); derived from the molecular layer, SI Section S3
    ktrMpDC      <- fixed(5.985073797735384); label("pDC transition rate, Ag_Lon to Ag_Mon (1/day)") # Table S3 (5.9887); derived from the molecular layer, SI Section S3
    ktrHpDC      <- fixed(4.025034796217001); label("pDC transition rate, Ag_Mon to Ag_H (1/day)") # Table S3 (4.0250); derived from the molecular layer, SI Section S3
    katrHpDC     <- fixed(2.000405122440113); label("pDC transition rate, Ag_H to Ag_Moff (1/day)") # Table S3 (2.0004); derived from the molecular layer, SI Section S3
    katrMpDC     <- fixed(1.182561965521271); label("pDC transition rate, Ag_Moff to Ag_Loff (1/day)") # Table S3 (1.1824); derived from the molecular layer, SI Section S3
    katrLpDC     <- fixed(0.393949613591618); label("pDC transition rate, Ag_Loff to pDC_off (1/day)") # Table S3 (0.3939); derived from the molecular layer, SI Section S3
    NmaxMHCmDC   <- fixed(12264.33347931808); label("Maximum MHC-II/antigen complexes on a mature mDC membrane (molecules)") # derived from the molecular layer, SI Section S3; value from published code Parameters.m
    NmaxMHCpDC   <- fixed(12264.33347931808); label("Maximum MHC-II/antigen complexes on a mature pDC membrane (molecules)") # derived from the molecular layer, SI Section S3; value from published code Parameters.m
    NT0          <- fixed(1445); label("Initial naive helper T cells in the lymph node (cells)") # Table S2 (1445 cells); ref 3
    KNT          <- fixed(400); label("MHC-II/peptide complexes for 50% naive T-cell activation (epitopes)") # Table S3 (400); ref 4
    KMT          <- fixed(40); label("MHC-II/peptide complexes for 50% memory T-cell activation (epitopes)") # Table S3 (40); ref 4
    kdtNT        <- fixed(0.0029); label("Naive helper T-cell death rate (1/day)") # Table S3 (0.0029); ref 4
    kdtAT        <- fixed(0.18); label("Activated helper T-cell death rate (1/day)") # Table S3 (0.18); ref 4
    kdtMT        <- fixed(2.7397e-4); label("Memory helper T-cell death rate (1/day)") # Table S3 (2.7397e-04); ref 4
    kdtFT        <- fixed(0.18); label("Functional helper T-cell death rate (1/day)") # Table S3 (0.1800); ref 4
    fAT          <- fixed(0.5); label("Fraction of activated T cells differentiating into memory T cells (unitless)") # Table S3 (0.5); model assumption
    kactNT       <- fixed(294.6588); label("Maximum naive helper T-cell activation rate (1/day)") # Table S3 (294.6588); fitted on BNT162b2 antibody data
    kactMT       <- fixed(935.8437988); label("Maximum memory helper T-cell activation rate (1/day)") # Table S3 (935.8438); fitted on BNT162b2 antibody data
    kprolAT      <- fixed(4.980899828); label("Maximum activated helper T-cell proliferation rate (1/day)") # Table S3 (4.9809); fitted on BNT162b2 antibody data
    BRN          <- fixed(75000); label("B-cell receptors per B cell (receptors/cell)") # Table S3 (75000); ref 4
    KaMid        <- fixed(1e-6); label("Ag/BCR association constant of the middle affinity subclone j=9 (1/pM); subclone j uses KaMid*2^(j-9)") # Table S3 (Ka,i from 3.9063e-09 to 2.56e-04 1/pM over the 17 subclones); ref 4
    KR           <- fixed(1); label("Occupied BCRs for 50% naive B-cell activation (unitless)") # Table S3 (1); ref 4
    kactNB       <- fixed(2.48); label("Maximum naive B-cell activation rate (1/day)") # Table S3 (2.48); ref 13
    kdtAB        <- fixed(0.2518); label("Activated B-cell death rate (1/day)") # Table S3 (0.2518); ref 4
    kdtNB        <- fixed(0.029); label("Naive B-cell death rate (1/day)") # Table S3 (0.029; the row is mislabelled k_dt^MB -- the description reads 'naive B-cells'); refs 13-15
    kdtMB        <- fixed(7.8278e-5); label("Memory B-cell death rate (1/day)") # Table S3 (7.8278e-05); ref 4
    g1           <- fixed(0.5); label("Fraction of activated B cells differentiating into memory B cells (unitless)") # Table S3 (0.5); ref 4
    g2           <- fixed(0.4); label("Fraction of activated B cells differentiating into short-lived plasma cells (unitless)") # Table S3 (0.4); ref 4
    CCN          <- fixed(66.17859977); label("Carrying capacity of 1 functional T cell for naive B-cell activation (unitless)") # Table S3 (66.1786); fitted on BNT162b2 antibody data
    kactMB       <- fixed(8.621199792); label("Maximum memory B-cell activation rate (1/day)") # Table S3 (8.6212); fitted on BNT162b2 antibody data
    kprolABN     <- fixed(6.061200085); label("Maximum proliferation rate of activated B cells derived from naive B cells (1/day)") # Table S3 (6.0612); fitted on BNT162b2 antibody data
    kprolABM     <- fixed(6.344299179); label("Maximum proliferation rate of activated B cells derived from memory B cells (1/day)") # Table S3 (6.3443); fitted on BNT162b2 antibody data
    kln2blPC     <- fixed(40.48209984); label("Plasma-cell migration rate, lymph node to blood (1/day)") # Table S3 (40.4821); fitted on BNT162b2 antibody data
    kdtSP        <- fixed(0.1383000154); label("Short-lived plasma-cell death rate (1/day)") # Table S3 (0.1383); fitted on BNT162b2 antibody data
    kdtLP        <- fixed(0.01430000263); label("Long-lived plasma-cell death rate (1/day)") # Table S3 (0.0143); fitted on BNT162b2 antibody data
    delayB       <- fixed(7.950590808e-07); label("Activated B-cell migration rate into the germinal centre (1/day)") # Table S3 (7.9506e-07); fitted on BNT162b2 antibody data
    NB0_1        <- fixed(8.653744563715161); label("Initial naive B cells, affinity subclone 1 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_2        <- fixed(24.845573920059596); label("Initial naive B cells, affinity subclone 2 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_3        <- fixed(61.97568463471193); label("Initial naive B cells, affinity subclone 3 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_4        <- fixed(134.31390362766157); label("Initial naive B cells, affinity subclone 4 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_5        <- fixed(252.899486900331); label("Initial naive B cells, affinity subclone 5 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_6        <- fixed(413.7159489978759); label("Initial naive B cells, affinity subclone 6 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_7        <- fixed(588.0089240916109); label("Initial naive B cells, affinity subclone 7 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_8        <- fixed(726.0941029474715); label("Initial naive B cells, affinity subclone 8 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_9        <- fixed(778.9852606331245); label("Initial naive B cells, affinity subclone 9 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_10       <- fixed(726.0941029474715); label("Initial naive B cells, affinity subclone 10 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_11       <- fixed(588.0089240916109); label("Initial naive B cells, affinity subclone 11 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_12       <- fixed(413.7159489978759); label("Initial naive B cells, affinity subclone 12 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_13       <- fixed(252.899486900331); label("Initial naive B cells, affinity subclone 13 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_14       <- fixed(134.31390362766157); label("Initial naive B cells, affinity subclone 14 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_15       <- fixed(61.97568463471193); label("Initial naive B cells, affinity subclone 15 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_16       <- fixed(24.845573920059596); label("Initial naive B cells, affinity subclone 16 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    NB0_17       <- fixed(8.653744563715161); label("Initial naive B cells, affinity subclone 17 (cells)") # Table S2 (NB0 spans 8.6537 to 778.9853 cells); normal-density allocation of 5200 naive B cells across the 17 subclones per published code Parameters.m
    kprodAb      <- fixed(8.64e8); label("Antibody secretion rate per plasma cell (molecules/(cell*day))") # Table S3 (8.64e+08); ref 4
    kdegAb       <- fixed(0.09179990135); label("Antibody elimination rate in blood (1/day)") # Table S3 (0.0918); fitted on BNT162b2 antibody data
    MWAb         <- fixed(150000); label("IgG molecular weight (g/mol)") # published code Parameters.m (MW_Ab = 150000 Da); used only for the pmol-to-ng/mL output conversion
    expSd        <- fixed(0); label("Exponential (log-scale) residual SD -- not reported by the paper") # SI Section S5 states Step 3 was fitted 'assuming an exponential error model' but reports no residual SD; fixed at 0 so the model simulates the deterministic typical-value trajectory the paper plots. See vignette Errata.
  })

  model({
    # ==== Homeostatic generation rates ================================
    # Table S3 reports these as 'computed to preserve steady state
    # without vaccination'. They are derived here from that constraint
    # (as the published code does) so the no-vaccine steady state is
    # exact; see vignette Errata on the two printed IS values.
    kbrNPIS <- NP0IS * (kdtNP + kis2blNP)
    kbrNPBL <- kdtNP * NP0BL - kis2blNP * NP0IS
    kbrMNIS <- MN0IS * (kdtMN + kis2blMN)
    kbrMNBL <- kdtMN * MN0BL - kis2blMN * MN0IS
    kbrmDCIS <- mDC0IS * (kdtImDC + kis2blmDC)
    kbrmDCBL <- kdtImDC * mDC0BL - kis2blmDC * mDC0IS
    kbrpDCIS <- pDC0IS * (kdtIpDC + kis2blpDC)
    kbrpDCBL <- kdtIpDC * pDC0BL - kis2blpDC * pDC0IS

    # ==== Antigen/BCR affinity ladder =================================
    # Ka_j = KaMid * 2^(j - 9), j = 1..17 (SI Table S3)
    Ka1 <- KaMid /  256
    Ka2 <- KaMid /  128
    Ka3 <- KaMid /  64
    Ka4 <- KaMid /  32
    Ka5 <- KaMid /  16
    Ka6 <- KaMid /  8
    Ka7 <- KaMid /  4
    Ka8 <- KaMid /  2
    Ka9 <- KaMid *  1
    Ka10 <- KaMid *  2
    Ka11 <- KaMid *  4
    Ka12 <- KaMid *  8
    Ka13 <- KaMid *  16
    Ka14 <- KaMid *  32
    Ka15 <- KaMid *  64
    Ka16 <- KaMid *  128
    Ka17 <- KaMid *  256

    # ==== Injection site: LNP-encapsulated mRNA (SI eq. 1) ============
    mRNAreg <- ksat * mRNA * 0.5 * (1 + tanh((mRNA - mRNAmax) / kslope))
    d/dt(mRNA) <- -kuppDC * mRNA * pDC_IS - kupmDC * mRNA * mDC_IS -
      kupMN * mRNA * MN_IS - kupNP * mRNA * NP_IS -
      kdmRNA * mRNA - kdeg * mRNA - mRNAreg
    # Note: the uptake terms in the mRNA balance carry no SFmRNA factor,
    # whereas the cell-side uptake terms below do. This asymmetry is
    # exactly as published (model_equations.m dydt(1) vs dydt(2..)).
    recFrac <- mRNA / (KmRNA + mRNA)

    # ==== Injection site: neutrophils (SI eqs. 2-4) ===================
    d/dt(NP_IS) <- kbrNPIS + krcNP * recFrac * NP_BL -
      kupNP * SFmRNA * mRNA * NP_IS - kdtNP * NP_IS - kis2blNP * NP_IS
    d/dt(NPL_IS) <- kupNP * SFmRNA * mRNA * NP_IS - kexpNP * NPL_IS -
      kdtNP * NPL_IS - kis2lnNPlnp * NPL_IS
    d/dt(NPAg_IS) <- kexpNP * NPL_IS - kdtNP * NPAg_IS - kis2lnNPag * NPAg_IS

    # ==== Injection site: monocytes (SI eqs. 5-7) =====================
    d/dt(MN_IS) <- kbrMNIS + krcMN * recFrac * MN_BL -
      kupMN * SFmRNA * mRNA * MN_IS - kdtMN * MN_IS - kis2blMN * MN_IS
    d/dt(MNL_IS) <- kupMN * SFmRNA * mRNA * MN_IS - kexpMN * MNL_IS -
      kdtMN * MNL_IS - kis2lnMNlnp * MNL_IS
    d/dt(MNAg_IS) <- kexpMN * MNL_IS - kdtMN * MNAg_IS - kis2lnMNag * MNAg_IS

    # ==== Injection site: mDC antigen-exposure chain (SI eqs. 8-15) ====
    d/dt(mDC_IS) <- kbrmDCIS + krcmDC * recFrac * mDC_BL -
      kupmDC * SFmRNA * mRNA * mDC_IS - kdtImDC * mDC_IS - kis2blmDC * mDC_IS
    d/dt(mDCL_IS) <- kupmDC * SFmRNA * mRNA * mDC_IS - kexpmDC * mDCL_IS -
      kdtmDC * mDCL_IS - kis2lnmDClnp * mDCL_IS
    d/dt(mDCAgLon_IS) <- kexpmDC * mDCL_IS - ktrMmDC * mDCAgLon_IS -
      kdtmDC * mDCAgLon_IS - kis2lnmDCag * mDCAgLon_IS
    d/dt(mDCAgMon_IS) <- ktrMmDC * mDCAgLon_IS - ktrHmDC * mDCAgMon_IS -
      kdtmDC * mDCAgMon_IS - kis2lnmDCag * mDCAgMon_IS
    d/dt(mDCAgH_IS) <- ktrHmDC * mDCAgMon_IS - katrHmDC * mDCAgH_IS -
      kdtmDC * mDCAgH_IS - kis2lnmDCag * mDCAgH_IS
    d/dt(mDCAgMoff_IS) <- katrHmDC * mDCAgH_IS - katrMmDC * mDCAgMoff_IS -
      kdtmDC * mDCAgMoff_IS - kis2lnmDCag * mDCAgMoff_IS
    d/dt(mDCAgLoff_IS) <- katrMmDC * mDCAgMoff_IS - katrLmDC * mDCAgLoff_IS -
      kdtmDC * mDCAgLoff_IS - kis2lnmDCag * mDCAgLoff_IS
    d/dt(mDCoff_IS) <- katrLmDC * mDCAgLoff_IS - kdtmDC * mDCoff_IS - kis2lnmDCag * mDCoff_IS

    # ==== Injection site: pDC antigen-exposure chain (SI eqs. 16-23) ====
    d/dt(pDC_IS) <- kbrpDCIS + krcpDC * recFrac * pDC_BL -
      kuppDC * SFmRNA * mRNA * pDC_IS - kdtIpDC * pDC_IS - kis2blpDC * pDC_IS
    d/dt(pDCL_IS) <- kuppDC * SFmRNA * mRNA * pDC_IS - kexppDC * pDCL_IS -
      kdtpDC * pDCL_IS - kis2lnpDClnp * pDCL_IS
    d/dt(pDCAgLon_IS) <- kexppDC * pDCL_IS - ktrMpDC * pDCAgLon_IS -
      kdtpDC * pDCAgLon_IS - kis2lnpDCag * pDCAgLon_IS
    d/dt(pDCAgMon_IS) <- ktrMpDC * pDCAgLon_IS - ktrHpDC * pDCAgMon_IS -
      kdtpDC * pDCAgMon_IS - kis2lnpDCag * pDCAgMon_IS
    d/dt(pDCAgH_IS) <- ktrHpDC * pDCAgMon_IS - katrHpDC * pDCAgH_IS -
      kdtpDC * pDCAgH_IS - kis2lnpDCag * pDCAgH_IS
    d/dt(pDCAgMoff_IS) <- katrHpDC * pDCAgH_IS - katrMpDC * pDCAgMoff_IS -
      kdtpDC * pDCAgMoff_IS - kis2lnpDCag * pDCAgMoff_IS
    d/dt(pDCAgLoff_IS) <- katrMpDC * pDCAgMoff_IS - katrLpDC * pDCAgLoff_IS -
      kdtpDC * pDCAgLoff_IS - kis2lnpDCag * pDCAgLoff_IS
    d/dt(pDCoff_IS) <- katrLpDC * pDCAgLoff_IS - kdtpDC * pDCoff_IS - kis2lnpDCag * pDCoff_IS

    # ==== Draining lymph node: neutrophils and monocytes (SI eqs. 24-27)
    d/dt(NPL_LN) <- kis2lnNPlnp * NPL_IS - kexpNP * NPL_LN - kdtNP * NPL_LN
    d/dt(NPAg_LN) <- kis2lnNPag * NPAg_IS + kexpNP * NPL_LN - kdtNP * NPAg_LN
    d/dt(MNL_LN) <- kis2lnMNlnp * MNL_IS - kexpMN * MNL_LN - kdtMN * MNL_LN
    d/dt(MNAg_LN) <- kis2lnMNag * MNAg_IS + kexpMN * MNL_LN - kdtMN * MNAg_LN

    # ==== Draining lymph node: mDC antigen-exposure chain (SI eqs. 28-34)
    d/dt(mDCL_LN) <- kis2lnmDClnp * mDCL_IS - kexpmDC * mDCL_LN - kdtmDC * mDCL_LN
    d/dt(mDCAgLon_LN) <- kis2lnmDCag * mDCAgLon_IS + kexpmDC * mDCL_LN -
      ktrMmDC * mDCAgLon_LN - kdtmDC * mDCAgLon_LN
    d/dt(mDCAgMon_LN) <- kis2lnmDCag * mDCAgMon_IS + ktrMmDC * mDCAgLon_LN -
      ktrHmDC * mDCAgMon_LN - kdtmDC * mDCAgMon_LN
    d/dt(mDCAgH_LN) <- kis2lnmDCag * mDCAgH_IS + ktrHmDC * mDCAgMon_LN -
      katrHmDC * mDCAgH_LN - kdtmDC * mDCAgH_LN
    d/dt(mDCAgMoff_LN) <- kis2lnmDCag * mDCAgMoff_IS + katrHmDC * mDCAgH_LN -
      katrMmDC * mDCAgMoff_LN - kdtmDC * mDCAgMoff_LN
    d/dt(mDCAgLoff_LN) <- kis2lnmDCag * mDCAgLoff_IS + katrMmDC * mDCAgMoff_LN -
      katrLmDC * mDCAgLoff_LN - kdtmDC * mDCAgLoff_LN
    d/dt(mDC_LN) <- kis2lnmDCag * mDCoff_IS + katrLmDC * mDCAgLoff_LN - kdtmDC * mDC_LN

    # ==== Draining lymph node: pDC antigen-exposure chain (SI eqs. 35-41)
    d/dt(pDCL_LN) <- kis2lnpDClnp * pDCL_IS - kexppDC * pDCL_LN - kdtpDC * pDCL_LN
    d/dt(pDCAgLon_LN) <- kis2lnpDCag * pDCAgLon_IS + kexppDC * pDCL_LN -
      ktrMpDC * pDCAgLon_LN - kdtpDC * pDCAgLon_LN
    d/dt(pDCAgMon_LN) <- kis2lnpDCag * pDCAgMon_IS + ktrMpDC * pDCAgLon_LN -
      ktrHpDC * pDCAgMon_LN - kdtpDC * pDCAgMon_LN
    d/dt(pDCAgH_LN) <- kis2lnpDCag * pDCAgH_IS + ktrHpDC * pDCAgMon_LN -
      katrHpDC * pDCAgH_LN - kdtpDC * pDCAgH_LN
    d/dt(pDCAgMoff_LN) <- kis2lnpDCag * pDCAgMoff_IS + katrHpDC * pDCAgH_LN -
      katrMpDC * pDCAgMoff_LN - kdtpDC * pDCAgMoff_LN
    d/dt(pDCAgLoff_LN) <- kis2lnpDCag * pDCAgLoff_IS + katrMpDC * pDCAgMoff_LN -
      katrLpDC * pDCAgLoff_LN - kdtpDC * pDCAgLoff_LN
    d/dt(pDC_LN) <- kis2lnpDCag * pDCoff_IS + katrLpDC * pDCAgLoff_LN - kdtpDC * pDC_LN

    # ==== Antigen presentation summaries (SI Section S1.4.2) ==========
    mDCAgL_LN <- mDCAgLon_LN + mDCAgLoff_LN
    mDCAgM_LN <- mDCAgMon_LN + mDCAgMoff_LN
    pDCAgL_LN <- pDCAgLon_LN + pDCAgLoff_LN
    pDCAgM_LN <- pDCAgMon_LN + pDCAgMoff_LN
    mDCmat <- mDCAgL_LN * pL + mDCAgM_LN * pM + mDCAgH_LN * pH
    pDCmat <- pDCAgL_LN * pL + pDCAgM_LN * pM + pDCAgH_LN * pH
    mDCtot <- mDCAgL_LN + mDCAgM_LN + mDCAgH_LN
    pDCtot <- pDCAgL_LN + pDCAgM_LN + pDCAgH_LN
    # Below one low-exposure-equivalent mature DC the published model
    # switches antigen presentation off entirely (a hard switch in
    # model_equations.m; it also guards the 0/0 at mDCtot = 0).
    if (mDCmat < pL) {
      NmedioMDC <- 0
    } else {
      NmedioMDC <- mDCmat * NmaxMHCmDC / mDCtot
    }
    if (pDCmat < pL) {
      NmedioPDC <- 0
    } else {
      NmedioPDC <- pDCmat * NmaxMHCpDC / pDCtot
    }
    tPool <- mDCtot + pDCtot + NT + AT_N + AT_M + MT
    sumMDCT <- mDCtot / (mDCtot + NT + AT_N + AT_M + MT)
    sumPDCT <- pDCtot / (pDCtot + NT + AT_N + AT_M + MT)
    if (mDCtot + pDCtot <= 0) {
      wMDC <- 0
      wPDC <- 0
    } else {
      wMDC <- mDCtot / (mDCtot + pDCtot)
      wPDC <- pDCtot / (mDCtot + pDCtot)
    }

    # ==== T-cell activation / proliferation functions (SI eqs. 42-46) =
    DmdcN <- sumMDCT * NmedioMDC / (NmedioMDC + KNT)
    DpdcN <- sumPDCT * NmedioPDC / (NmedioPDC + KNT)
    DmdcM <- sumMDCT * NmedioMDC / (NmedioMDC + KMT)
    DpdcM <- sumPDCT * NmedioPDC / (NmedioPDC + KMT)
    DN <- DmdcN * wMDC + DpdcN * wPDC
    DM <- DmdcM * wMDC + DpdcM * wPDC
    EmdcN <- sumMDCT * (NmedioMDC - KNT) / (NmedioMDC + KNT)
    EmdcM <- sumMDCT * (NmedioMDC - KMT) / (NmedioMDC + KMT)
    EpdcN <- sumPDCT * (NmedioPDC - KNT) / (NmedioPDC + KNT)
    EpdcM <- sumPDCT * (NmedioPDC - KMT) / (NmedioPDC + KMT)
    ENraw <- EmdcN * wMDC + EpdcN * wPDC
    EMraw <- EmdcM * wMDC + EpdcM * wPDC
    # The published code floors E at zero for the first 0.1 day, which
    # prevents a spurious negative proliferation signal before any
    # antigen-presenting DC has reached the lymph node.
    if (t < 0.1) {
      EN <- max(ENraw, 0)
      EM <- max(EMraw, 0)
    } else {
      EN <- ENraw
      EM <- EMraw
    }
    d/dt(NT) <- kdtNT * (NT0 - NT) - kactNT * DN * NT
    d/dt(AT_N) <- kactNT * DN * NT + kprolAT * EN * AT_N - kdtAT * AT_N
    d/dt(MT) <- kprolAT * (1 - EN) * fAT * AT_N +
      kprolAT * (1 - EM) * fAT * AT_M - kactMT * DM * MT - kdtMT * MT
    d/dt(AT_M) <- kactMT * DM * MT + kprolAT * EM * AT_M - kdtAT * AT_M
    d/dt(FT) <- kprolAT * (1 - EN) * (1 - fAT) * AT_N +
      kprolAT * (1 - EM) * (1 - fAT) * AT_M - kdtFT * FT

    # ==== Free antigen in the lymph node (SI eq. 47) ==================
    # Total antigen available to B cells, from the MHC-II/antigen
    # exposure of the mature DC pools.
    wLM <- abs(pL - pM) / 2 + pL
    wMH <- abs(pM - pH) / 2 + pM
    wHH <- abs(pH - 1) / 2 + pH
    maxMHCmDC <- NmaxMHCmDC / NAvo * 1e12
    maxMHCpDC <- NmaxMHCpDC / NAvo * 1e12
    AgLN <- expAg * (maxMHCmDC * (wLM * mDCAgL_LN + wMH * mDCAgM_LN +
      wHH * mDCAgH_LN) + maxMHCpDC * (wLM * pDCAgL_LN + wMH * pDCAgM_LN +
      wHH * pDCAgH_LN))
    AgTot <- AgLN / VLN
    bcrSc <- BRN / NAvo * 1e12 / VLN
    B1 <- (NB1 + ABN1 + GCB1 + ABM1 + MB1) * bcrSc
    B2 <- (NB2 + ABN2 + GCB2 + ABM2 + MB2) * bcrSc
    B3 <- (NB3 + ABN3 + GCB3 + ABM3 + MB3) * bcrSc
    B4 <- (NB4 + ABN4 + GCB4 + ABM4 + MB4) * bcrSc
    B5 <- (NB5 + ABN5 + GCB5 + ABM5 + MB5) * bcrSc
    B6 <- (NB6 + ABN6 + GCB6 + ABM6 + MB6) * bcrSc
    B7 <- (NB7 + ABN7 + GCB7 + ABM7 + MB7) * bcrSc
    B8 <- (NB8 + ABN8 + GCB8 + ABM8 + MB8) * bcrSc
    B9 <- (NB9 + ABN9 + GCB9 + ABM9 + MB9) * bcrSc
    B10 <- (NB10 + ABN10 + GCB10 + ABM10 + MB10) * bcrSc
    B11 <- (NB11 + ABN11 + GCB11 + ABM11 + MB11) * bcrSc
    B12 <- (NB12 + ABN12 + GCB12 + ABM12 + MB12) * bcrSc
    B13 <- (NB13 + ABN13 + GCB13 + ABM13 + MB13) * bcrSc
    B14 <- (NB14 + ABN14 + GCB14 + ABM14 + MB14) * bcrSc
    B15 <- (NB15 + ABN15 + GCB15 + ABM15 + MB15) * bcrSc
    B16 <- (NB16 + ABN16 + GCB16 + ABM16 + MB16) * bcrSc
    B17 <- (NB17 + ABN17 + GCB17 + ABM17 + MB17) * bcrSc
    # The published implementation solves the multi-affinity binding
    # equilibrium  AgTot = x + sum_j B_j*Ka_j*x/(1 + Ka_j*x)  for the free
    # concentration x with a numerical root finder (MATLAB fzero).
    # rxode2 has no in-RHS root finder, so the identical equation is
    # solved by Newton iteration written out in closed form. g(x) is
    # increasing and concave with g(0) < 0, so Newton started from x = 0
    # converges monotonically from below. The first step reduces to
    # AgTot/(1 + sum_j B_j*Ka_j); four further steps follow. Three steps
    # already reproduce the activation function F to machine precision
    # over the whole simulated trajectory (see vignette Errata).
    sBK <- B1 * Ka1 + B2 * Ka2 + B3 * Ka3 + B4 * Ka4 + B5 * Ka5 + B6 * Ka6 + 
      B7 * Ka7 + B8 * Ka8 + B9 * Ka9 + B10 * Ka10 + B11 * Ka11 + B12 * Ka12 + 
      B13 * Ka13 + B14 * Ka14 + B15 * Ka15 + B16 * Ka16 + B17 * Ka17
    agf <- AgTot / (1 + sBK)
    gS1 <- B1 * Ka1 * agf / (1 + Ka1 * agf) + B2 * Ka2 * agf / (1 + Ka2 * agf) + B3 * Ka3 * agf / (1 + Ka3 * agf) + B4 * Ka4 * agf / (1 + Ka4 * agf) + B5 * Ka5 * agf / (1 + Ka5 * agf) + B6 * Ka6 * agf / (1 + Ka6 * agf) + 
      B7 * Ka7 * agf / (1 + Ka7 * agf) + B8 * Ka8 * agf / (1 + Ka8 * agf) + B9 * Ka9 * agf / (1 + Ka9 * agf) + B10 * Ka10 * agf / (1 + Ka10 * agf) + B11 * Ka11 * agf / (1 + Ka11 * agf) + B12 * Ka12 * agf / (1 + Ka12 * agf) + 
      B13 * Ka13 * agf / (1 + Ka13 * agf) + B14 * Ka14 * agf / (1 + Ka14 * agf) + B15 * Ka15 * agf / (1 + Ka15 * agf) + B16 * Ka16 * agf / (1 + Ka16 * agf) + B17 * Ka17 * agf / (1 + Ka17 * agf)
    gD1 <- B1 * Ka1 / (1 + Ka1 * agf)^2 + B2 * Ka2 / (1 + Ka2 * agf)^2 + B3 * Ka3 / (1 + Ka3 * agf)^2 + B4 * Ka4 / (1 + Ka4 * agf)^2 + B5 * Ka5 / (1 + Ka5 * agf)^2 + B6 * Ka6 / (1 + Ka6 * agf)^2 + 
      B7 * Ka7 / (1 + Ka7 * agf)^2 + B8 * Ka8 / (1 + Ka8 * agf)^2 + B9 * Ka9 / (1 + Ka9 * agf)^2 + B10 * Ka10 / (1 + Ka10 * agf)^2 + B11 * Ka11 / (1 + Ka11 * agf)^2 + B12 * Ka12 / (1 + Ka12 * agf)^2 + 
      B13 * Ka13 / (1 + Ka13 * agf)^2 + B14 * Ka14 / (1 + Ka14 * agf)^2 + B15 * Ka15 / (1 + Ka15 * agf)^2 + B16 * Ka16 / (1 + Ka16 * agf)^2 + B17 * Ka17 / (1 + Ka17 * agf)^2
    agf <- agf - (agf + gS1 - AgTot) / (1 + gD1)
    gS2 <- B1 * Ka1 * agf / (1 + Ka1 * agf) + B2 * Ka2 * agf / (1 + Ka2 * agf) + B3 * Ka3 * agf / (1 + Ka3 * agf) + B4 * Ka4 * agf / (1 + Ka4 * agf) + B5 * Ka5 * agf / (1 + Ka5 * agf) + B6 * Ka6 * agf / (1 + Ka6 * agf) + 
      B7 * Ka7 * agf / (1 + Ka7 * agf) + B8 * Ka8 * agf / (1 + Ka8 * agf) + B9 * Ka9 * agf / (1 + Ka9 * agf) + B10 * Ka10 * agf / (1 + Ka10 * agf) + B11 * Ka11 * agf / (1 + Ka11 * agf) + B12 * Ka12 * agf / (1 + Ka12 * agf) + 
      B13 * Ka13 * agf / (1 + Ka13 * agf) + B14 * Ka14 * agf / (1 + Ka14 * agf) + B15 * Ka15 * agf / (1 + Ka15 * agf) + B16 * Ka16 * agf / (1 + Ka16 * agf) + B17 * Ka17 * agf / (1 + Ka17 * agf)
    gD2 <- B1 * Ka1 / (1 + Ka1 * agf)^2 + B2 * Ka2 / (1 + Ka2 * agf)^2 + B3 * Ka3 / (1 + Ka3 * agf)^2 + B4 * Ka4 / (1 + Ka4 * agf)^2 + B5 * Ka5 / (1 + Ka5 * agf)^2 + B6 * Ka6 / (1 + Ka6 * agf)^2 + 
      B7 * Ka7 / (1 + Ka7 * agf)^2 + B8 * Ka8 / (1 + Ka8 * agf)^2 + B9 * Ka9 / (1 + Ka9 * agf)^2 + B10 * Ka10 / (1 + Ka10 * agf)^2 + B11 * Ka11 / (1 + Ka11 * agf)^2 + B12 * Ka12 / (1 + Ka12 * agf)^2 + 
      B13 * Ka13 / (1 + Ka13 * agf)^2 + B14 * Ka14 / (1 + Ka14 * agf)^2 + B15 * Ka15 / (1 + Ka15 * agf)^2 + B16 * Ka16 / (1 + Ka16 * agf)^2 + B17 * Ka17 / (1 + Ka17 * agf)^2
    agf <- agf - (agf + gS2 - AgTot) / (1 + gD2)
    gS3 <- B1 * Ka1 * agf / (1 + Ka1 * agf) + B2 * Ka2 * agf / (1 + Ka2 * agf) + B3 * Ka3 * agf / (1 + Ka3 * agf) + B4 * Ka4 * agf / (1 + Ka4 * agf) + B5 * Ka5 * agf / (1 + Ka5 * agf) + B6 * Ka6 * agf / (1 + Ka6 * agf) + 
      B7 * Ka7 * agf / (1 + Ka7 * agf) + B8 * Ka8 * agf / (1 + Ka8 * agf) + B9 * Ka9 * agf / (1 + Ka9 * agf) + B10 * Ka10 * agf / (1 + Ka10 * agf) + B11 * Ka11 * agf / (1 + Ka11 * agf) + B12 * Ka12 * agf / (1 + Ka12 * agf) + 
      B13 * Ka13 * agf / (1 + Ka13 * agf) + B14 * Ka14 * agf / (1 + Ka14 * agf) + B15 * Ka15 * agf / (1 + Ka15 * agf) + B16 * Ka16 * agf / (1 + Ka16 * agf) + B17 * Ka17 * agf / (1 + Ka17 * agf)
    gD3 <- B1 * Ka1 / (1 + Ka1 * agf)^2 + B2 * Ka2 / (1 + Ka2 * agf)^2 + B3 * Ka3 / (1 + Ka3 * agf)^2 + B4 * Ka4 / (1 + Ka4 * agf)^2 + B5 * Ka5 / (1 + Ka5 * agf)^2 + B6 * Ka6 / (1 + Ka6 * agf)^2 + 
      B7 * Ka7 / (1 + Ka7 * agf)^2 + B8 * Ka8 / (1 + Ka8 * agf)^2 + B9 * Ka9 / (1 + Ka9 * agf)^2 + B10 * Ka10 / (1 + Ka10 * agf)^2 + B11 * Ka11 / (1 + Ka11 * agf)^2 + B12 * Ka12 / (1 + Ka12 * agf)^2 + 
      B13 * Ka13 / (1 + Ka13 * agf)^2 + B14 * Ka14 / (1 + Ka14 * agf)^2 + B15 * Ka15 / (1 + Ka15 * agf)^2 + B16 * Ka16 / (1 + Ka16 * agf)^2 + B17 * Ka17 / (1 + Ka17 * agf)^2
    agf <- agf - (agf + gS3 - AgTot) / (1 + gD3)
    gS4 <- B1 * Ka1 * agf / (1 + Ka1 * agf) + B2 * Ka2 * agf / (1 + Ka2 * agf) + B3 * Ka3 * agf / (1 + Ka3 * agf) + B4 * Ka4 * agf / (1 + Ka4 * agf) + B5 * Ka5 * agf / (1 + Ka5 * agf) + B6 * Ka6 * agf / (1 + Ka6 * agf) + 
      B7 * Ka7 * agf / (1 + Ka7 * agf) + B8 * Ka8 * agf / (1 + Ka8 * agf) + B9 * Ka9 * agf / (1 + Ka9 * agf) + B10 * Ka10 * agf / (1 + Ka10 * agf) + B11 * Ka11 * agf / (1 + Ka11 * agf) + B12 * Ka12 * agf / (1 + Ka12 * agf) + 
      B13 * Ka13 * agf / (1 + Ka13 * agf) + B14 * Ka14 * agf / (1 + Ka14 * agf) + B15 * Ka15 * agf / (1 + Ka15 * agf) + B16 * Ka16 * agf / (1 + Ka16 * agf) + B17 * Ka17 * agf / (1 + Ka17 * agf)
    gD4 <- B1 * Ka1 / (1 + Ka1 * agf)^2 + B2 * Ka2 / (1 + Ka2 * agf)^2 + B3 * Ka3 / (1 + Ka3 * agf)^2 + B4 * Ka4 / (1 + Ka4 * agf)^2 + B5 * Ka5 / (1 + Ka5 * agf)^2 + B6 * Ka6 / (1 + Ka6 * agf)^2 + 
      B7 * Ka7 / (1 + Ka7 * agf)^2 + B8 * Ka8 / (1 + Ka8 * agf)^2 + B9 * Ka9 / (1 + Ka9 * agf)^2 + B10 * Ka10 / (1 + Ka10 * agf)^2 + B11 * Ka11 / (1 + Ka11 * agf)^2 + B12 * Ka12 / (1 + Ka12 * agf)^2 + 
      B13 * Ka13 / (1 + Ka13 * agf)^2 + B14 * Ka14 / (1 + Ka14 * agf)^2 + B15 * Ka15 / (1 + Ka15 * agf)^2 + B16 * Ka16 / (1 + Ka16 * agf)^2 + B17 * Ka17 / (1 + Ka17 * agf)^2
    agf <- agf - (agf + gS4 - AgTot) / (1 + gD4)

    # ==== B-cell receptor occupancy functions (SI eqs. 48-50) =========
    ro1 <- Ka1 * agf / (1 + Ka1 * agf)
    ro2 <- Ka2 * agf / (1 + Ka2 * agf)
    ro3 <- Ka3 * agf / (1 + Ka3 * agf)
    ro4 <- Ka4 * agf / (1 + Ka4 * agf)
    ro5 <- Ka5 * agf / (1 + Ka5 * agf)
    ro6 <- Ka6 * agf / (1 + Ka6 * agf)
    ro7 <- Ka7 * agf / (1 + Ka7 * agf)
    ro8 <- Ka8 * agf / (1 + Ka8 * agf)
    ro9 <- Ka9 * agf / (1 + Ka9 * agf)
    ro10 <- Ka10 * agf / (1 + Ka10 * agf)
    ro11 <- Ka11 * agf / (1 + Ka11 * agf)
    ro12 <- Ka12 * agf / (1 + Ka12 * agf)
    ro13 <- Ka13 * agf / (1 + Ka13 * agf)
    ro14 <- Ka14 * agf / (1 + Ka14 * agf)
    ro15 <- Ka15 * agf / (1 + Ka15 * agf)
    ro16 <- Ka16 * agf / (1 + Ka16 * agf)
    ro17 <- Ka17 * agf / (1 + Ka17 * agf)
    R1 <- ro1 * BRN
    R2 <- ro2 * BRN
    R3 <- ro3 * BRN
    R4 <- ro4 * BRN
    R5 <- ro5 * BRN
    R6 <- ro6 * BRN
    R7 <- ro7 * BRN
    R8 <- ro8 * BRN
    R9 <- ro9 * BRN
    R10 <- ro10 * BRN
    R11 <- ro11 * BRN
    R12 <- ro12 * BRN
    R13 <- ro13 * BRN
    R14 <- ro14 * BRN
    R15 <- ro15 * BRN
    R16 <- ro16 * BRN
    R17 <- ro17 * BRN
    Ff1 <- R1 / (KR + R1)
    Gf1 <- (1 - ro1) * Ff1
    Hf1 <- (R1 - KR) / (R1 + KR)
    Ff2 <- R2 / (KR + R2)
    Gf2 <- (1 - ro2) * Ff2
    Hf2 <- (R2 - KR) / (R2 + KR)
    Ff3 <- R3 / (KR + R3)
    Gf3 <- (1 - ro3) * Ff3
    Hf3 <- (R3 - KR) / (R3 + KR)
    Ff4 <- R4 / (KR + R4)
    Gf4 <- (1 - ro4) * Ff4
    Hf4 <- (R4 - KR) / (R4 + KR)
    Ff5 <- R5 / (KR + R5)
    Gf5 <- (1 - ro5) * Ff5
    Hf5 <- (R5 - KR) / (R5 + KR)
    Ff6 <- R6 / (KR + R6)
    Gf6 <- (1 - ro6) * Ff6
    Hf6 <- (R6 - KR) / (R6 + KR)
    Ff7 <- R7 / (KR + R7)
    Gf7 <- (1 - ro7) * Ff7
    Hf7 <- (R7 - KR) / (R7 + KR)
    Ff8 <- R8 / (KR + R8)
    Gf8 <- (1 - ro8) * Ff8
    Hf8 <- (R8 - KR) / (R8 + KR)
    Ff9 <- R9 / (KR + R9)
    Gf9 <- (1 - ro9) * Ff9
    Hf9 <- (R9 - KR) / (R9 + KR)
    Ff10 <- R10 / (KR + R10)
    Gf10 <- (1 - ro10) * Ff10
    Hf10 <- (R10 - KR) / (R10 + KR)
    Ff11 <- R11 / (KR + R11)
    Gf11 <- (1 - ro11) * Ff11
    Hf11 <- (R11 - KR) / (R11 + KR)
    Ff12 <- R12 / (KR + R12)
    Gf12 <- (1 - ro12) * Ff12
    Hf12 <- (R12 - KR) / (R12 + KR)
    Ff13 <- R13 / (KR + R13)
    Gf13 <- (1 - ro13) * Ff13
    Hf13 <- (R13 - KR) / (R13 + KR)
    Ff14 <- R14 / (KR + R14)
    Gf14 <- (1 - ro14) * Ff14
    Hf14 <- (R14 - KR) / (R14 + KR)
    Ff15 <- R15 / (KR + R15)
    Gf15 <- (1 - ro15) * Ff15
    Hf15 <- (R15 - KR) / (R15 + KR)
    Ff16 <- R16 / (KR + R16)
    Gf16 <- (1 - ro16) * Ff16
    Hf16 <- (R16 - KR) / (R16 + KR)
    Ff17 <- R17 / (KR + R17)
    Gf17 <- (1 - ro17) * Ff17
    Hf17 <- (R17 - KR) / (R17 + KR)

    # ==== Helper-T support functions (SI eqs. 51-52) ==================
    bPool <- NB1 + ABN1 + GCB1 + ABM1 + MB1 + NB2 + ABN2 + GCB2 + ABM2 + MB2 + NB3 + ABN3 + GCB3 + ABM3 + MB3 + 
      NB4 + ABN4 + GCB4 + ABM4 + MB4 + NB5 + ABN5 + GCB5 + ABM5 + MB5 + NB6 + ABN6 + GCB6 + ABM6 + MB6 + 
      NB7 + ABN7 + GCB7 + ABM7 + MB7 + NB8 + ABN8 + GCB8 + ABM8 + MB8 + NB9 + ABN9 + GCB9 + ABM9 + MB9 + 
      NB10 + ABN10 + GCB10 + ABM10 + MB10 + NB11 + ABN11 + GCB11 + ABM11 + MB11 + NB12 + ABN12 + GCB12 + ABM12 + MB12 + 
      NB13 + ABN13 + GCB13 + ABM13 + MB13 + NB14 + ABN14 + GCB14 + ABM14 + MB14 + NB15 + ABN15 + GCB15 + ABM15 + MB15 + 
      NB16 + ABN16 + GCB16 + ABM16 + MB16 + NB17 + ABN17 + GCB17 + ABM17 + MB17
    PN <- CCN * FT / (CCN * FT + bPool)
    # CC_M is fixed at 10 x CC_N (Table S3 model assumption)
    PM <- (CCN * 10) * FT / ((CCN * 10) * FT + bPool)

    # ==== B-cell subclones in the lymph node (SI eqs. 53-59) ==========
    # 17 affinity subclones, indexed j = 1..17 by increasing Ka.
    d/dt(NB1) <- kdtNB * (NB0_1 - NB1) - kactNB * Ff1 * PN * NB1
    d/dt(NB2) <- kdtNB * (NB0_2 - NB2) - kactNB * Ff2 * PN * NB2
    d/dt(NB3) <- kdtNB * (NB0_3 - NB3) - kactNB * Ff3 * PN * NB3
    d/dt(NB4) <- kdtNB * (NB0_4 - NB4) - kactNB * Ff4 * PN * NB4
    d/dt(NB5) <- kdtNB * (NB0_5 - NB5) - kactNB * Ff5 * PN * NB5
    d/dt(NB6) <- kdtNB * (NB0_6 - NB6) - kactNB * Ff6 * PN * NB6
    d/dt(NB7) <- kdtNB * (NB0_7 - NB7) - kactNB * Ff7 * PN * NB7
    d/dt(NB8) <- kdtNB * (NB0_8 - NB8) - kactNB * Ff8 * PN * NB8
    d/dt(NB9) <- kdtNB * (NB0_9 - NB9) - kactNB * Ff9 * PN * NB9
    d/dt(NB10) <- kdtNB * (NB0_10 - NB10) - kactNB * Ff10 * PN * NB10
    d/dt(NB11) <- kdtNB * (NB0_11 - NB11) - kactNB * Ff11 * PN * NB11
    d/dt(NB12) <- kdtNB * (NB0_12 - NB12) - kactNB * Ff12 * PN * NB12
    d/dt(NB13) <- kdtNB * (NB0_13 - NB13) - kactNB * Ff13 * PN * NB13
    d/dt(NB14) <- kdtNB * (NB0_14 - NB14) - kactNB * Ff14 * PN * NB14
    d/dt(NB15) <- kdtNB * (NB0_15 - NB15) - kactNB * Ff15 * PN * NB15
    d/dt(NB16) <- kdtNB * (NB0_16 - NB16) - kactNB * Ff16 * PN * NB16
    d/dt(NB17) <- kdtNB * (NB0_17 - NB17) - kactNB * Ff17 * PN * NB17
    d/dt(ABN1) <- kactNB * Gf1 * PN * NB1 - delayB * ABN1 - kdtAB * ABN1
    d/dt(ABN2) <- kactNB * Gf2 * PN * NB2 - delayB * ABN2 - kdtAB * ABN2
    d/dt(ABN3) <- kactNB * Gf3 * PN * NB3 - delayB * ABN3 - kdtAB * ABN3
    d/dt(ABN4) <- kactNB * Gf4 * PN * NB4 - delayB * ABN4 - kdtAB * ABN4
    d/dt(ABN5) <- kactNB * Gf5 * PN * NB5 - delayB * ABN5 - kdtAB * ABN5
    d/dt(ABN6) <- kactNB * Gf6 * PN * NB6 - delayB * ABN6 - kdtAB * ABN6
    d/dt(ABN7) <- kactNB * Gf7 * PN * NB7 - delayB * ABN7 - kdtAB * ABN7
    d/dt(ABN8) <- kactNB * Gf8 * PN * NB8 - delayB * ABN8 - kdtAB * ABN8
    d/dt(ABN9) <- kactNB * Gf9 * PN * NB9 - delayB * ABN9 - kdtAB * ABN9
    d/dt(ABN10) <- kactNB * Gf10 * PN * NB10 - delayB * ABN10 - kdtAB * ABN10
    d/dt(ABN11) <- kactNB * Gf11 * PN * NB11 - delayB * ABN11 - kdtAB * ABN11
    d/dt(ABN12) <- kactNB * Gf12 * PN * NB12 - delayB * ABN12 - kdtAB * ABN12
    d/dt(ABN13) <- kactNB * Gf13 * PN * NB13 - delayB * ABN13 - kdtAB * ABN13
    d/dt(ABN14) <- kactNB * Gf14 * PN * NB14 - delayB * ABN14 - kdtAB * ABN14
    d/dt(ABN15) <- kactNB * Gf15 * PN * NB15 - delayB * ABN15 - kdtAB * ABN15
    d/dt(ABN16) <- kactNB * Gf16 * PN * NB16 - delayB * ABN16 - kdtAB * ABN16
    d/dt(ABN17) <- kactNB * Gf17 * PN * NB17 - delayB * ABN17 - kdtAB * ABN17
    d/dt(GCB1) <- delayB * ABN1 + kprolABN * Hf1 * PN * GCB1 - kdtAB * GCB1
    d/dt(GCB2) <- delayB * ABN2 + kprolABN * Hf2 * PN * GCB2 - kdtAB * GCB2
    d/dt(GCB3) <- delayB * ABN3 + kprolABN * Hf3 * PN * GCB3 - kdtAB * GCB3
    d/dt(GCB4) <- delayB * ABN4 + kprolABN * Hf4 * PN * GCB4 - kdtAB * GCB4
    d/dt(GCB5) <- delayB * ABN5 + kprolABN * Hf5 * PN * GCB5 - kdtAB * GCB5
    d/dt(GCB6) <- delayB * ABN6 + kprolABN * Hf6 * PN * GCB6 - kdtAB * GCB6
    d/dt(GCB7) <- delayB * ABN7 + kprolABN * Hf7 * PN * GCB7 - kdtAB * GCB7
    d/dt(GCB8) <- delayB * ABN8 + kprolABN * Hf8 * PN * GCB8 - kdtAB * GCB8
    d/dt(GCB9) <- delayB * ABN9 + kprolABN * Hf9 * PN * GCB9 - kdtAB * GCB9
    d/dt(GCB10) <- delayB * ABN10 + kprolABN * Hf10 * PN * GCB10 - kdtAB * GCB10
    d/dt(GCB11) <- delayB * ABN11 + kprolABN * Hf11 * PN * GCB11 - kdtAB * GCB11
    d/dt(GCB12) <- delayB * ABN12 + kprolABN * Hf12 * PN * GCB12 - kdtAB * GCB12
    d/dt(GCB13) <- delayB * ABN13 + kprolABN * Hf13 * PN * GCB13 - kdtAB * GCB13
    d/dt(GCB14) <- delayB * ABN14 + kprolABN * Hf14 * PN * GCB14 - kdtAB * GCB14
    d/dt(GCB15) <- delayB * ABN15 + kprolABN * Hf15 * PN * GCB15 - kdtAB * GCB15
    d/dt(GCB16) <- delayB * ABN16 + kprolABN * Hf16 * PN * GCB16 - kdtAB * GCB16
    d/dt(GCB17) <- delayB * ABN17 + kprolABN * Hf17 * PN * GCB17 - kdtAB * GCB17
    d/dt(MB1) <- kprolABN * (1 - Hf1) * PN * g1 * GCB1 +
      kprolABM * (1 - Hf1) * PM * g1 * ABM1 -
      kactMB * Ff1 * PM * MB1 - kdtMB * MB1
    d/dt(MB2) <- kprolABN * (1 - Hf2) * PN * g1 * GCB2 +
      kprolABM * (1 - Hf2) * PM * g1 * ABM2 -
      kactMB * Ff2 * PM * MB2 - kdtMB * MB2
    d/dt(MB3) <- kprolABN * (1 - Hf3) * PN * g1 * GCB3 +
      kprolABM * (1 - Hf3) * PM * g1 * ABM3 -
      kactMB * Ff3 * PM * MB3 - kdtMB * MB3
    d/dt(MB4) <- kprolABN * (1 - Hf4) * PN * g1 * GCB4 +
      kprolABM * (1 - Hf4) * PM * g1 * ABM4 -
      kactMB * Ff4 * PM * MB4 - kdtMB * MB4
    d/dt(MB5) <- kprolABN * (1 - Hf5) * PN * g1 * GCB5 +
      kprolABM * (1 - Hf5) * PM * g1 * ABM5 -
      kactMB * Ff5 * PM * MB5 - kdtMB * MB5
    d/dt(MB6) <- kprolABN * (1 - Hf6) * PN * g1 * GCB6 +
      kprolABM * (1 - Hf6) * PM * g1 * ABM6 -
      kactMB * Ff6 * PM * MB6 - kdtMB * MB6
    d/dt(MB7) <- kprolABN * (1 - Hf7) * PN * g1 * GCB7 +
      kprolABM * (1 - Hf7) * PM * g1 * ABM7 -
      kactMB * Ff7 * PM * MB7 - kdtMB * MB7
    d/dt(MB8) <- kprolABN * (1 - Hf8) * PN * g1 * GCB8 +
      kprolABM * (1 - Hf8) * PM * g1 * ABM8 -
      kactMB * Ff8 * PM * MB8 - kdtMB * MB8
    d/dt(MB9) <- kprolABN * (1 - Hf9) * PN * g1 * GCB9 +
      kprolABM * (1 - Hf9) * PM * g1 * ABM9 -
      kactMB * Ff9 * PM * MB9 - kdtMB * MB9
    d/dt(MB10) <- kprolABN * (1 - Hf10) * PN * g1 * GCB10 +
      kprolABM * (1 - Hf10) * PM * g1 * ABM10 -
      kactMB * Ff10 * PM * MB10 - kdtMB * MB10
    d/dt(MB11) <- kprolABN * (1 - Hf11) * PN * g1 * GCB11 +
      kprolABM * (1 - Hf11) * PM * g1 * ABM11 -
      kactMB * Ff11 * PM * MB11 - kdtMB * MB11
    d/dt(MB12) <- kprolABN * (1 - Hf12) * PN * g1 * GCB12 +
      kprolABM * (1 - Hf12) * PM * g1 * ABM12 -
      kactMB * Ff12 * PM * MB12 - kdtMB * MB12
    d/dt(MB13) <- kprolABN * (1 - Hf13) * PN * g1 * GCB13 +
      kprolABM * (1 - Hf13) * PM * g1 * ABM13 -
      kactMB * Ff13 * PM * MB13 - kdtMB * MB13
    d/dt(MB14) <- kprolABN * (1 - Hf14) * PN * g1 * GCB14 +
      kprolABM * (1 - Hf14) * PM * g1 * ABM14 -
      kactMB * Ff14 * PM * MB14 - kdtMB * MB14
    d/dt(MB15) <- kprolABN * (1 - Hf15) * PN * g1 * GCB15 +
      kprolABM * (1 - Hf15) * PM * g1 * ABM15 -
      kactMB * Ff15 * PM * MB15 - kdtMB * MB15
    d/dt(MB16) <- kprolABN * (1 - Hf16) * PN * g1 * GCB16 +
      kprolABM * (1 - Hf16) * PM * g1 * ABM16 -
      kactMB * Ff16 * PM * MB16 - kdtMB * MB16
    d/dt(MB17) <- kprolABN * (1 - Hf17) * PN * g1 * GCB17 +
      kprolABM * (1 - Hf17) * PM * g1 * ABM17 -
      kactMB * Ff17 * PM * MB17 - kdtMB * MB17
    d/dt(ABM1) <- kactMB * Gf1 * PM * MB1 +
      kprolABM * Hf1 * PM * ABM1 - kdtAB * ABM1
    d/dt(ABM2) <- kactMB * Gf2 * PM * MB2 +
      kprolABM * Hf2 * PM * ABM2 - kdtAB * ABM2
    d/dt(ABM3) <- kactMB * Gf3 * PM * MB3 +
      kprolABM * Hf3 * PM * ABM3 - kdtAB * ABM3
    d/dt(ABM4) <- kactMB * Gf4 * PM * MB4 +
      kprolABM * Hf4 * PM * ABM4 - kdtAB * ABM4
    d/dt(ABM5) <- kactMB * Gf5 * PM * MB5 +
      kprolABM * Hf5 * PM * ABM5 - kdtAB * ABM5
    d/dt(ABM6) <- kactMB * Gf6 * PM * MB6 +
      kprolABM * Hf6 * PM * ABM6 - kdtAB * ABM6
    d/dt(ABM7) <- kactMB * Gf7 * PM * MB7 +
      kprolABM * Hf7 * PM * ABM7 - kdtAB * ABM7
    d/dt(ABM8) <- kactMB * Gf8 * PM * MB8 +
      kprolABM * Hf8 * PM * ABM8 - kdtAB * ABM8
    d/dt(ABM9) <- kactMB * Gf9 * PM * MB9 +
      kprolABM * Hf9 * PM * ABM9 - kdtAB * ABM9
    d/dt(ABM10) <- kactMB * Gf10 * PM * MB10 +
      kprolABM * Hf10 * PM * ABM10 - kdtAB * ABM10
    d/dt(ABM11) <- kactMB * Gf11 * PM * MB11 +
      kprolABM * Hf11 * PM * ABM11 - kdtAB * ABM11
    d/dt(ABM12) <- kactMB * Gf12 * PM * MB12 +
      kprolABM * Hf12 * PM * ABM12 - kdtAB * ABM12
    d/dt(ABM13) <- kactMB * Gf13 * PM * MB13 +
      kprolABM * Hf13 * PM * ABM13 - kdtAB * ABM13
    d/dt(ABM14) <- kactMB * Gf14 * PM * MB14 +
      kprolABM * Hf14 * PM * ABM14 - kdtAB * ABM14
    d/dt(ABM15) <- kactMB * Gf15 * PM * MB15 +
      kprolABM * Hf15 * PM * ABM15 - kdtAB * ABM15
    d/dt(ABM16) <- kactMB * Gf16 * PM * MB16 +
      kprolABM * Hf16 * PM * ABM16 - kdtAB * ABM16
    d/dt(ABM17) <- kactMB * Gf17 * PM * MB17 +
      kprolABM * Hf17 * PM * ABM17 - kdtAB * ABM17
    d/dt(SP1) <- kprolABN * (1 - Hf1) * PN * g2 * GCB1 +
      kprolABM * (1 - Hf1) * PM * g2 * ABM1 -
      (kdtSP + kln2blPC) * SP1
    d/dt(SP2) <- kprolABN * (1 - Hf2) * PN * g2 * GCB2 +
      kprolABM * (1 - Hf2) * PM * g2 * ABM2 -
      (kdtSP + kln2blPC) * SP2
    d/dt(SP3) <- kprolABN * (1 - Hf3) * PN * g2 * GCB3 +
      kprolABM * (1 - Hf3) * PM * g2 * ABM3 -
      (kdtSP + kln2blPC) * SP3
    d/dt(SP4) <- kprolABN * (1 - Hf4) * PN * g2 * GCB4 +
      kprolABM * (1 - Hf4) * PM * g2 * ABM4 -
      (kdtSP + kln2blPC) * SP4
    d/dt(SP5) <- kprolABN * (1 - Hf5) * PN * g2 * GCB5 +
      kprolABM * (1 - Hf5) * PM * g2 * ABM5 -
      (kdtSP + kln2blPC) * SP5
    d/dt(SP6) <- kprolABN * (1 - Hf6) * PN * g2 * GCB6 +
      kprolABM * (1 - Hf6) * PM * g2 * ABM6 -
      (kdtSP + kln2blPC) * SP6
    d/dt(SP7) <- kprolABN * (1 - Hf7) * PN * g2 * GCB7 +
      kprolABM * (1 - Hf7) * PM * g2 * ABM7 -
      (kdtSP + kln2blPC) * SP7
    d/dt(SP8) <- kprolABN * (1 - Hf8) * PN * g2 * GCB8 +
      kprolABM * (1 - Hf8) * PM * g2 * ABM8 -
      (kdtSP + kln2blPC) * SP8
    d/dt(SP9) <- kprolABN * (1 - Hf9) * PN * g2 * GCB9 +
      kprolABM * (1 - Hf9) * PM * g2 * ABM9 -
      (kdtSP + kln2blPC) * SP9
    d/dt(SP10) <- kprolABN * (1 - Hf10) * PN * g2 * GCB10 +
      kprolABM * (1 - Hf10) * PM * g2 * ABM10 -
      (kdtSP + kln2blPC) * SP10
    d/dt(SP11) <- kprolABN * (1 - Hf11) * PN * g2 * GCB11 +
      kprolABM * (1 - Hf11) * PM * g2 * ABM11 -
      (kdtSP + kln2blPC) * SP11
    d/dt(SP12) <- kprolABN * (1 - Hf12) * PN * g2 * GCB12 +
      kprolABM * (1 - Hf12) * PM * g2 * ABM12 -
      (kdtSP + kln2blPC) * SP12
    d/dt(SP13) <- kprolABN * (1 - Hf13) * PN * g2 * GCB13 +
      kprolABM * (1 - Hf13) * PM * g2 * ABM13 -
      (kdtSP + kln2blPC) * SP13
    d/dt(SP14) <- kprolABN * (1 - Hf14) * PN * g2 * GCB14 +
      kprolABM * (1 - Hf14) * PM * g2 * ABM14 -
      (kdtSP + kln2blPC) * SP14
    d/dt(SP15) <- kprolABN * (1 - Hf15) * PN * g2 * GCB15 +
      kprolABM * (1 - Hf15) * PM * g2 * ABM15 -
      (kdtSP + kln2blPC) * SP15
    d/dt(SP16) <- kprolABN * (1 - Hf16) * PN * g2 * GCB16 +
      kprolABM * (1 - Hf16) * PM * g2 * ABM16 -
      (kdtSP + kln2blPC) * SP16
    d/dt(SP17) <- kprolABN * (1 - Hf17) * PN * g2 * GCB17 +
      kprolABM * (1 - Hf17) * PM * g2 * ABM17 -
      (kdtSP + kln2blPC) * SP17
    d/dt(LP1) <- kprolABN * (1 - Hf1) * PN * (1 - g1 - g2) * GCB1 +
      kprolABM * (1 - Hf1) * PM * (1 - g1 - g2) * ABM1 -
      (kdtLP + kln2blPC) * LP1
    d/dt(LP2) <- kprolABN * (1 - Hf2) * PN * (1 - g1 - g2) * GCB2 +
      kprolABM * (1 - Hf2) * PM * (1 - g1 - g2) * ABM2 -
      (kdtLP + kln2blPC) * LP2
    d/dt(LP3) <- kprolABN * (1 - Hf3) * PN * (1 - g1 - g2) * GCB3 +
      kprolABM * (1 - Hf3) * PM * (1 - g1 - g2) * ABM3 -
      (kdtLP + kln2blPC) * LP3
    d/dt(LP4) <- kprolABN * (1 - Hf4) * PN * (1 - g1 - g2) * GCB4 +
      kprolABM * (1 - Hf4) * PM * (1 - g1 - g2) * ABM4 -
      (kdtLP + kln2blPC) * LP4
    d/dt(LP5) <- kprolABN * (1 - Hf5) * PN * (1 - g1 - g2) * GCB5 +
      kprolABM * (1 - Hf5) * PM * (1 - g1 - g2) * ABM5 -
      (kdtLP + kln2blPC) * LP5
    d/dt(LP6) <- kprolABN * (1 - Hf6) * PN * (1 - g1 - g2) * GCB6 +
      kprolABM * (1 - Hf6) * PM * (1 - g1 - g2) * ABM6 -
      (kdtLP + kln2blPC) * LP6
    d/dt(LP7) <- kprolABN * (1 - Hf7) * PN * (1 - g1 - g2) * GCB7 +
      kprolABM * (1 - Hf7) * PM * (1 - g1 - g2) * ABM7 -
      (kdtLP + kln2blPC) * LP7
    d/dt(LP8) <- kprolABN * (1 - Hf8) * PN * (1 - g1 - g2) * GCB8 +
      kprolABM * (1 - Hf8) * PM * (1 - g1 - g2) * ABM8 -
      (kdtLP + kln2blPC) * LP8
    d/dt(LP9) <- kprolABN * (1 - Hf9) * PN * (1 - g1 - g2) * GCB9 +
      kprolABM * (1 - Hf9) * PM * (1 - g1 - g2) * ABM9 -
      (kdtLP + kln2blPC) * LP9
    d/dt(LP10) <- kprolABN * (1 - Hf10) * PN * (1 - g1 - g2) * GCB10 +
      kprolABM * (1 - Hf10) * PM * (1 - g1 - g2) * ABM10 -
      (kdtLP + kln2blPC) * LP10
    d/dt(LP11) <- kprolABN * (1 - Hf11) * PN * (1 - g1 - g2) * GCB11 +
      kprolABM * (1 - Hf11) * PM * (1 - g1 - g2) * ABM11 -
      (kdtLP + kln2blPC) * LP11
    d/dt(LP12) <- kprolABN * (1 - Hf12) * PN * (1 - g1 - g2) * GCB12 +
      kprolABM * (1 - Hf12) * PM * (1 - g1 - g2) * ABM12 -
      (kdtLP + kln2blPC) * LP12
    d/dt(LP13) <- kprolABN * (1 - Hf13) * PN * (1 - g1 - g2) * GCB13 +
      kprolABM * (1 - Hf13) * PM * (1 - g1 - g2) * ABM13 -
      (kdtLP + kln2blPC) * LP13
    d/dt(LP14) <- kprolABN * (1 - Hf14) * PN * (1 - g1 - g2) * GCB14 +
      kprolABM * (1 - Hf14) * PM * (1 - g1 - g2) * ABM14 -
      (kdtLP + kln2blPC) * LP14
    d/dt(LP15) <- kprolABN * (1 - Hf15) * PN * (1 - g1 - g2) * GCB15 +
      kprolABM * (1 - Hf15) * PM * (1 - g1 - g2) * ABM15 -
      (kdtLP + kln2blPC) * LP15
    d/dt(LP16) <- kprolABN * (1 - Hf16) * PN * (1 - g1 - g2) * GCB16 +
      kprolABM * (1 - Hf16) * PM * (1 - g1 - g2) * ABM16 -
      (kdtLP + kln2blPC) * LP16
    d/dt(LP17) <- kprolABN * (1 - Hf17) * PN * (1 - g1 - g2) * GCB17 +
      kprolABM * (1 - Hf17) * PM * (1 - g1 - g2) * ABM17 -
      (kdtLP + kln2blPC) * LP17

    # ==== Blood: innate populations (SI eqs. 60-63) ===================
    d/dt(NP_BL) <- kbrNPBL - krcNP * recFrac * NP_BL - kdtNP * NP_BL +
      kis2blNP * NP_IS
    d/dt(MN_BL) <- kbrMNBL - krcMN * recFrac * MN_BL - kdtMN * MN_BL +
      kis2blMN * MN_IS
    d/dt(mDC_BL) <- kbrmDCBL - krcmDC * recFrac * mDC_BL - kdtImDC * mDC_BL +
      kis2blmDC * mDC_IS
    d/dt(pDC_BL) <- kbrpDCBL - krcpDC * recFrac * pDC_BL - kdtIpDC * pDC_BL +
      kis2blpDC * pDC_IS

    # ==== Blood: plasma cells and antibody (SI eqs. 64-66) ============
    # kprodAb is molecules/(cell*day); /NAvo*1e12 converts to pmol.
    abSc <- kprodAb / NAvo * 1e12
    d/dt(SPBL1) <- kln2blPC * SP1 - kdtSP * SPBL1
    d/dt(SPBL2) <- kln2blPC * SP2 - kdtSP * SPBL2
    d/dt(SPBL3) <- kln2blPC * SP3 - kdtSP * SPBL3
    d/dt(SPBL4) <- kln2blPC * SP4 - kdtSP * SPBL4
    d/dt(SPBL5) <- kln2blPC * SP5 - kdtSP * SPBL5
    d/dt(SPBL6) <- kln2blPC * SP6 - kdtSP * SPBL6
    d/dt(SPBL7) <- kln2blPC * SP7 - kdtSP * SPBL7
    d/dt(SPBL8) <- kln2blPC * SP8 - kdtSP * SPBL8
    d/dt(SPBL9) <- kln2blPC * SP9 - kdtSP * SPBL9
    d/dt(SPBL10) <- kln2blPC * SP10 - kdtSP * SPBL10
    d/dt(SPBL11) <- kln2blPC * SP11 - kdtSP * SPBL11
    d/dt(SPBL12) <- kln2blPC * SP12 - kdtSP * SPBL12
    d/dt(SPBL13) <- kln2blPC * SP13 - kdtSP * SPBL13
    d/dt(SPBL14) <- kln2blPC * SP14 - kdtSP * SPBL14
    d/dt(SPBL15) <- kln2blPC * SP15 - kdtSP * SPBL15
    d/dt(SPBL16) <- kln2blPC * SP16 - kdtSP * SPBL16
    d/dt(SPBL17) <- kln2blPC * SP17 - kdtSP * SPBL17
    d/dt(LPBL1) <- kln2blPC * LP1 - kdtLP * LPBL1
    d/dt(LPBL2) <- kln2blPC * LP2 - kdtLP * LPBL2
    d/dt(LPBL3) <- kln2blPC * LP3 - kdtLP * LPBL3
    d/dt(LPBL4) <- kln2blPC * LP4 - kdtLP * LPBL4
    d/dt(LPBL5) <- kln2blPC * LP5 - kdtLP * LPBL5
    d/dt(LPBL6) <- kln2blPC * LP6 - kdtLP * LPBL6
    d/dt(LPBL7) <- kln2blPC * LP7 - kdtLP * LPBL7
    d/dt(LPBL8) <- kln2blPC * LP8 - kdtLP * LPBL8
    d/dt(LPBL9) <- kln2blPC * LP9 - kdtLP * LPBL9
    d/dt(LPBL10) <- kln2blPC * LP10 - kdtLP * LPBL10
    d/dt(LPBL11) <- kln2blPC * LP11 - kdtLP * LPBL11
    d/dt(LPBL12) <- kln2blPC * LP12 - kdtLP * LPBL12
    d/dt(LPBL13) <- kln2blPC * LP13 - kdtLP * LPBL13
    d/dt(LPBL14) <- kln2blPC * LP14 - kdtLP * LPBL14
    d/dt(LPBL15) <- kln2blPC * LP15 - kdtLP * LPBL15
    d/dt(LPBL16) <- kln2blPC * LP16 - kdtLP * LPBL16
    d/dt(LPBL17) <- kln2blPC * LP17 - kdtLP * LPBL17
    d/dt(Ab1) <- abSc * (SPBL1 + LPBL1) - kdegAb * Ab1
    d/dt(Ab2) <- abSc * (SPBL2 + LPBL2) - kdegAb * Ab2
    d/dt(Ab3) <- abSc * (SPBL3 + LPBL3) - kdegAb * Ab3
    d/dt(Ab4) <- abSc * (SPBL4 + LPBL4) - kdegAb * Ab4
    d/dt(Ab5) <- abSc * (SPBL5 + LPBL5) - kdegAb * Ab5
    d/dt(Ab6) <- abSc * (SPBL6 + LPBL6) - kdegAb * Ab6
    d/dt(Ab7) <- abSc * (SPBL7 + LPBL7) - kdegAb * Ab7
    d/dt(Ab8) <- abSc * (SPBL8 + LPBL8) - kdegAb * Ab8
    d/dt(Ab9) <- abSc * (SPBL9 + LPBL9) - kdegAb * Ab9
    d/dt(Ab10) <- abSc * (SPBL10 + LPBL10) - kdegAb * Ab10
    d/dt(Ab11) <- abSc * (SPBL11 + LPBL11) - kdegAb * Ab11
    d/dt(Ab12) <- abSc * (SPBL12 + LPBL12) - kdegAb * Ab12
    d/dt(Ab13) <- abSc * (SPBL13 + LPBL13) - kdegAb * Ab13
    d/dt(Ab14) <- abSc * (SPBL14 + LPBL14) - kdegAb * Ab14
    d/dt(Ab15) <- abSc * (SPBL15 + LPBL15) - kdegAb * Ab15
    d/dt(Ab16) <- abSc * (SPBL16 + LPBL16) - kdegAb * Ab16
    d/dt(Ab17) <- abSc * (SPBL17 + LPBL17) - kdegAb * Ab17

    # ==== Initial conditions (SI Table S2) ============================
    # The mRNA state starts empty; the vaccine dose is given as a bolus
    # into mRNA (pmol) from the event table, matching the published
    # code, which adds mRNA0_pmol to the mRNA state at each injection.
    NP_IS(0) <- NP0IS
    MN_IS(0) <- MN0IS
    mDC_IS(0) <- mDC0IS
    pDC_IS(0) <- pDC0IS
    NT(0) <- NT0
    NP_BL(0) <- NP0BL
    MN_BL(0) <- MN0BL
    mDC_BL(0) <- mDC0BL
    pDC_BL(0) <- pDC0BL
    NB1(0) <- NB0_1
    NB2(0) <- NB0_2
    NB3(0) <- NB0_3
    NB4(0) <- NB0_4
    NB5(0) <- NB0_5
    NB6(0) <- NB0_6
    NB7(0) <- NB0_7
    NB8(0) <- NB0_8
    NB9(0) <- NB0_9
    NB10(0) <- NB0_10
    NB11(0) <- NB0_11
    NB12(0) <- NB0_12
    NB13(0) <- NB0_13
    NB14(0) <- NB0_14
    NB15(0) <- NB0_15
    NB16(0) <- NB0_16
    NB17(0) <- NB0_17
    # every other state starts at 0 (Table S2)

    # ==== Observation ================================================
    # Total anti-RBD IgG in serum. Cc is the canonical single-output
    # observation name; here it is the antibody concentration in blood
    # (ng/mL), NOT a concentration of the dosed mRNA. The conversion is
    # pmol -> ng/mL: sum(Ab)*1e-12 mol * MWAb g/mol * 1e9 ng/g / (VBL*1e3 mL).
    AbTot <- Ab1 + Ab2 + Ab3 + Ab4 + Ab5 + Ab6 + 
      Ab7 + Ab8 + Ab9 + Ab10 + Ab11 + Ab12 + 
      Ab13 + Ab14 + Ab15 + Ab16 + Ab17
    Cc <- AbTot * 1e-12 * MWAb * 1e9 / (VBL * 1e3)
    Cc ~ lnorm(expSd)
  })
}
