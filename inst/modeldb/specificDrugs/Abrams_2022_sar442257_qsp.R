Abrams_2022_sar442257_qsp <- function() {
  description <- paste(
    "In vitro (human peripheral-blood T cells, CD38+ PBMCs and multiple myeloma",
    "cells). QSP. Rule-based quantitative systems pharmacology model of the",
    "CD38xCD28xCD3 trispecific T-cell engager SAR442257 in multiple myeloma: the",
    "dose-prediction form of the in vitro model (760 species plus 13 cumulative",
    "flux trackers, 773 ODE states). Naive, effector memory and active CD4+",
    "and CD8+ T cells, multiple myeloma (MM) cells, a lumped CD38+ PBMC population",
    "and shed soluble CD38 bind free drug through the CD3, CD28 and CD38 arms; any",
    "pair of cells whose arms can be bridged by drug forms a synapse with a",
    "collision-driven rate (Figure 1C). Effector memory cells activate in a",
    "CD3-arm synapse at a constant rate; naive cells need drug on both CD3 and CD28",
    "(a Michaelis-Menten AND gate). Active T cells proliferate and kill MM cells",
    "and PBMCs in synapse, T cell-MM synapses become killing-resistant, and active",
    "T cells and synapsed PBMCs/MM cells release IFN-gamma, TNF-alpha, IL-6 and",
    "IL-10 (outputs only). Deterministic mechanism model with no IIV and no error",
    "model. Drug is dosed as molecules into tsAb (1 nM is Vol * 6.022e14",
    "molecules); setting kon_CD28 = 0 gives the comparator CD38xCD3 bispecific.",
    sep = " "
  )
  reference <- paste(
    "Abrams RE, Pierre K, El-Murr N, Seung E, Wu L, Luna E, Mehta R, Li J,",
    "Larabi K, Ahmed M, Pelekanou V, Yang ZY, van de Velde H, Stamatelos SK.",
    "Quantitative systems pharmacology modeling sheds light into the dose",
    "response relationship of a trispecific T cell engager in multiple myeloma.",
    "Sci Rep. 2022;12:10976. doi:10.1038/s41598-022-14726-5",
    sep = " "
  )
  vignette <- "Abrams_2022_sar442257_qsp"
  units <- list(
    time = "h",
    dosing = "molecules (SAR442257 molecules per well; 1 nM = Vol * 6.022e14 molecules)",
    concentration = "nM (free SAR442257 in Cc); cells and molecules per well for all states"
  )

  # The ODE system is the deposited Matlab model-generation code
  # (Supplementary files MOESM1-3: the main script with
  # Nm = 'tsAb_FullODEs_DosePred', DefineSpecified_2cellSyn and
  # Receptors_Def_perCell) executed unchanged, with its output
  # tsAb_FullODEs_DosePred_Original.m translated term by term: X(idx.<name>)
  # becomes the state <name> and ps.<name> the parameter <name>. Four
  # generator-level substitutions follow the printed sources and are listed in
  # the vignette: koffc_<R> -> koff_<R> (Figure 1C writes the bridge
  # dissociation as koffR1 + koffR2); one antigen density per lineage
  # (Table S2 'CD3per_CD4 (A, EM, N)'); kprod_MM_<cytokine> -> the
  # kprod_TRGT_<cytokine> value (Table S2 'Assume same as TRGT rates'); and
  # the logistic CD38+ PBMC proliferation term, whose rate kpr_TRGT is 0 in
  # Table S2, is dropped together with its unreported carrying capacity
  # C_M_TRGT. The generator's constant input ps.Dose is replaced by dose
  # events into tsAb. The repeated synapse-formation, bridge-formation and
  # naive-activation sub-expressions are named once (kf_*, kb1_*, kb2_*, ka_*)
  # instead of being inlined at every use; the vignette shows the factored and
  # inlined systems give identical solutions.
  paper_specific_compartment_pattern <- paste0(
    "^(tsAb|sCD38|MM|TRGT|IFNg|TNFa|IL6|IL10|SynFormation|CD[48]_(N|EM|A)|",
    "(Tcell|MMcell|TRGTcell|CD3|CD28|CD38)(Syn|Deg)|[RS]_[A-Za-z0-9_]+)$"
  )

  covariateData <- list()

  compartmentData <- list(
    tsAb = list(analyte = "SAR442257 (free)", units = "molecules", specimen = "whole blood", verified = TRUE),
    CD8_N = list(analyte = "naive CD8+ T cells (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    CD8_EM = list(
      analyte = "effector memory CD8+ T cells (free)",
      units = "cells",
      specimen = "blood cell",
      verified = TRUE
    ),
    CD8_A = list(analyte = "active CD8+ T cells (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    CD4_N = list(analyte = "naive CD4+ T cells (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    CD4_EM = list(
      analyte = "effector memory CD4+ T cells (free)",
      units = "cells",
      specimen = "blood cell",
      verified = TRUE
    ),
    CD4_A = list(analyte = "active CD4+ T cells (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    MM = list(analyte = "multiple myeloma cells (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    TRGT = list(analyte = "CD38+ PBMCs (free)", units = "cells", specimen = "blood cell", verified = TRUE),
    sCD38 = list(analyte = "soluble CD38", units = "molecules", specimen = "whole blood", verified = TRUE),
    S_CD4_A_CD28_MM_CD38 = list(
      analyte = "active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38 = list(
      analyte = "active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28 = list(
      analyte = "active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38 = list(
      analyte = "active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38 = list(
      analyte = "active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38 = list(
      analyte = "effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38 = list(
      analyte = "effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38 = list(
      analyte = "effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38 = list(
      analyte = "naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38 = list(
      analyte = "naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38 = list(
      analyte = "naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38 = list(
      analyte = "active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38 = list(
      analyte = "active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28 = list(
      analyte = "active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38 = list(
      analyte = "active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38 = list(
      analyte = "active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38 = list(
      analyte = "effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38 = list(
      analyte = "effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38 = list(
      analyte = "effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38 = list(
      analyte = "naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38 = list(
      analyte = "naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38 = list(
      analyte = "naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38 = list(
      analyte = "multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38 = list(
      analyte = "multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT = list(
      analyte = "killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT = list(
      analyte = "killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT = list(
      analyte = "killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT = list(
      analyte = "killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT = list(
      analyte = "killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT = list(
      analyte = "killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "synapses",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_N_CD3 = list(
      analyte = "CD3 on free naive CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_EM_CD3 = list(
      analyte = "CD3 on free effector memory CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_A_CD3 = list(
      analyte = "CD3 on free active CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_N_CD3 = list(
      analyte = "CD3 on free naive CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_EM_CD3 = list(
      analyte = "CD3 on free effector memory CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_A_CD3 = list(
      analyte = "CD3 on free active CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_MM_CD38 = list(
      analyte = "CD38 on free multiple myeloma cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_TRGT_CD38 = list(
      analyte = "CD38 on free CD38+ PBMCs, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_sCD38_CD38 = list(
      analyte = "soluble CD38 (free)",
      units = "molecules",
      specimen = "whole blood",
      verified = TRUE
    ),
    R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on free naive CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on free effector memory CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on free active CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on free naive CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on free effector memory CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on free active CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_MM_CD38_tsAb = list(
      analyte = "CD38 on free multiple myeloma cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on free CD38+ PBMCs, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_sCD38_CD38_tsAb = list(
      analyte = "soluble CD38 bound to SAR442257",
      units = "molecules",
      specimen = "whole blood",
      verified = TRUE
    ),
    R_CD8_N_CD28 = list(
      analyte = "CD28 on free naive CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_EM_CD28 = list(
      analyte = "CD28 on free effector memory CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_A_CD28 = list(
      analyte = "CD28 on free active CD8+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_N_CD28 = list(
      analyte = "CD28 on free naive CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_EM_CD28 = list(
      analyte = "CD28 on free effector memory CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_A_CD28 = list(
      analyte = "CD28 on free active CD4+ T cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_MM_CD28 = list(
      analyte = "CD28 on free multiple myeloma cells, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_TRGT_CD28 = list(
      analyte = "CD28 on free CD38+ PBMCs, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on free naive CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on free effector memory CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on free active CD8+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on free naive CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on free effector memory CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on free active CD4+ T cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_MM_CD28_tsAb = list(
      analyte = "CD28 on free multiple myeloma cells, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on free CD38+ PBMCs, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 = list(
      analyte = "CD3 on the naive CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3 = list(
      analyte = "CD3 on the active CD4+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD4+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3 = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3 = list(
      analyte = "CD3 on the active CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 = list(
      analyte = "CD3 on the effector memory CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_TRGT_CD38 = list(
      analyte = "CD38 on the CD38+ PBMC within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38 = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD4+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD4+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb = list(
      analyte = "CD3 on the naive CD4+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb = list(
      analyte = "CD3 on the active CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb = list(
      analyte = "CD3 on the effector memory CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb = list(
      analyte = "CD38 on the CD38+ PBMC within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb = list(
      analyte = "CD38 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 = list(
      analyte = "CD28 on the naive CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28 = list(
      analyte = "CD28 on the active CD4+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD4+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28 = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28 = list(
      analyte = "CD28 on the active CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 = list(
      analyte = "CD28 on the effector memory CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_TRGT_CD28 = list(
      analyte = "CD28 on the CD38+ PBMC within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28 = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, unbound",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD8+ T cell within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD4+ T cell within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD4+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb = list(
      analyte = "CD28 on the naive CD4+ T cell within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb = list(
      analyte = "CD28 on the active CD8+ T cell within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb = list(
      analyte = "CD28 on the effector memory CD8+ T cell within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb = list(
      analyte = "CD28 on the CD38+ PBMC within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb = list(
      analyte = "CD28 on the multiple myeloma cell within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses, bound to SAR442257",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_EM_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_N_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD4+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within active CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_EM_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within effector memory CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_A_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-active CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-effector memory CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD4_N_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-naive CD4+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_A_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-active CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_EM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-effector memory CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_CD8_N_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-naive CD8+ T cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD28_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_N_CD3_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within naive CD8+ T cell (CD3 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_MM_CD38_Br = list(
      analyte = "SAR442257 bridges within multiple myeloma cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_MM_CD28_TRGT_CD38_Br = list(
      analyte = "SAR442257 bridges within multiple myeloma cell (CD28 arm)-CD38+ PBMC (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD28_MM_CD38_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD4+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD28_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD4_A_CD3_MM_CD38_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD4+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD28_MM_CD38_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD8+ T cell (CD28 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD28_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD28 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    S_CD8_A_CD3_MM_CD38_MUT_Br = list(
      analyte = "SAR442257 bridges within killing-resistant active CD8+ T cell (CD3 arm)-multiple myeloma cell (CD38 arm) synapses",
      units = "molecules",
      specimen = "blood cell",
      verified = TRUE
    ),
    IFNg = list(analyte = "IFN-gamma", units = "molecules", specimen = "whole blood", verified = TRUE),
    TNFa = list(analyte = "TNF-alpha", units = "molecules", specimen = "whole blood", verified = TRUE),
    IL6 = list(analyte = "IL-6", units = "molecules", specimen = "whole blood", verified = TRUE),
    IL10 = list(analyte = "IL-10", units = "molecules", specimen = "whole blood", verified = TRUE),
    TcellSyn = list(
      analyte = "cumulative cell flux tracker TcellSyn",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    TcellDeg = list(
      analyte = "cumulative cell flux tracker TcellDeg",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    MMcellSyn = list(
      analyte = "cumulative cell flux tracker MMcellSyn",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    MMcellDeg = list(
      analyte = "cumulative cell flux tracker MMcellDeg",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    TRGTcellSyn = list(
      analyte = "cumulative cell flux tracker TRGTcellSyn",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    TRGTcellDeg = list(
      analyte = "cumulative cell flux tracker TRGTcellDeg",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD3Syn = list(
      analyte = "cumulative receptor flux tracker CD3Syn",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD3Deg = list(
      analyte = "cumulative receptor flux tracker CD3Deg",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD28Syn = list(
      analyte = "cumulative receptor flux tracker CD28Syn",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD28Deg = list(
      analyte = "cumulative receptor flux tracker CD28Deg",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD38Syn = list(
      analyte = "cumulative receptor flux tracker CD38Syn",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    CD38Deg = list(
      analyte = "cumulative receptor flux tracker CD38Deg",
      units = "molecules",
      specimen = "not applicable",
      verified = TRUE
    ),
    SynFormation = list(
      analyte = "cumulative cell flux tracker SynFormation",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (human peripheral-blood T cells, CD38+ PBMCs and multiple myeloma cells)",
    n_subjects = paste(
      "Not applicable: the model was calibrated to in vitro T-cell activation, cytotoxicity",
      "(RPMI-8226 and KMS-11 cell lines) and MIMIC cytokine-assay data, and this parameter",
      "set is the dose-prediction module representing a peripheral-blood sample of a human",
      "multiple myeloma patient (Table S2 initial values)."
    ),
    disease_state = "Relapsed and refractory multiple myeloma (in vitro representation)",
    dose_range = paste(
      "Fixed drug concentrations from 8.4e-4 to 0.672 nM in the dose-prediction simulations",
      "(Figure 3), starting at the MABEL dose; 72-h simulations."
    ),
    notes = paste(
      "The paper's in vitro virtual population (sampled from the MIMIC calibration) is not",
      "deposited; this file carries the single parameter set of Table S2. The model informed",
      "the design of the first-in-human trial NCT04401020 of SAR442257."
    )
  )

  ini({
    # T-cell synthesis (zero in vitro)
    ksyn_CD8_N <- fixed(0); label("Synthesis rate of naive CD8+ T cells (cells/h)") # Table S2 'ksyn_CD8_N' = 0, 'Assume no synthesis in vitro'
    ksyn_CD8_EM <- fixed(0); label("Synthesis rate of effector memory CD8+ T cells (cells/h)") # Table S2 'ksyn_CD8_EM' = 0, 'Assume no synthesis in vitro'
    ksyn_CD4_N <- fixed(0); label("Synthesis rate of naive CD4+ T cells (cells/h)") # Table S2 'ksyn_CD4_N' = 0, 'Assume no synthesis in vitro'
    ksyn_CD4_EM <- fixed(0); label("Synthesis rate of effector memory CD4+ T cells (cells/h)") # Table S2 'ksyn_CD4_EM' = 0, 'Assume no synthesis in vitro'

    # Apoptotic (degradation) rates
    kdeg_CD8_N <- fixed(0.0064); label("In vitro apoptotic rate of naive CD8+ T cells (1/h)") # Table S2 'kdeg_CD8_N' = 0.0064, literature [74]
    kdeg_CD8_EM <- fixed(0.0076); label("In vitro apoptotic rate of effector memory CD8+ T cells (1/h)") # Table S2 'kdeg_CD8_EM' = 0.0076, literature [74]
    kdeg_CD8_A <- fixed(0.051); label("In vitro apoptotic rate of active CD8+ T cells (1/h)") # Table S2 'kdeg_CD8_A' = 0.051, literature [75]; checked against the MIMIC assay
    kdeg_CD4_N <- fixed(0.0085); label("In vitro apoptotic rate of naive CD4+ T cells (1/h)") # Table S2 'kdeg_CD4_N' = 0.0085, literature [76]
    kdeg_CD4_EM <- fixed(0.0076); label("In vitro apoptotic rate of effector memory CD4+ T cells (1/h)") # Table S2 'kdeg_CD4_EM' = 0.0076, literature [76]
    kdeg_CD4_A <- fixed(0.01); label("In vitro apoptotic rate of active CD4+ T cells (1/h)") # Table S2 'kdeg_CD4_A' = 0.01, literature [75]; checked against the MIMIC assay
    kdeg_MM <- fixed(0.0025); label("In vitro apoptotic rate of multiple myeloma cells (1/h)") # Table S2 'kdeg_MM' = 0.0025, literature [77]
    kdeg_TRGT <- fixed(0.0036); label("In vitro apoptotic rate of CD38+ PBMCs (1/h)") # Table S2 'kdeg_TRGT' = 0.0036, literature [78]; checked against the MIMIC assay

    # Proliferation
    kpr_CD8 <- fixed(0.0452); label("In vitro proliferation rate of active CD8+ T cells (1/h)") # Table S2 'kpr_CD8' = 0.0452, calibrated to the activation assay
    kpr_CD4 <- fixed(0.0099); label("In vitro proliferation rate of active CD4+ T cells (1/h)") # Table S2 'kpr_CD4' = 0.0099, calibrated to the activation assay
    kpr_MM <- fixed(1.00e-4); label("In vitro proliferation rate of multiple myeloma cells (1/h)") # Table S2 'kpr_MM' = 1.00e-4, literature [77]
    C_M_MM <- fixed(1.00e6); label("In vitro multiple myeloma carrying capacity (cells)") # Table S2 'C_M_MM' = 1.00e6, literature [77]

    # Killing and resistance
    kkillMM_CD8 <- fixed(64.2663); label("Tumour-cell killing rate by an active CD8+ T cell in synapse (1/h)") # Table S2 'kkillMM_CD8' = 64.2663, calibrated to the cytotoxicity assay
    kkillMM_CD4 <- fixed(6.4266); label("Tumour-cell killing rate by an active CD4+ T cell in synapse (1/h)") # Table S2 'kkillMM_CD4' = 6.4266, one order of magnitude below CD8
    kmut_SYN <- fixed(14.2651); label("Rate at which T cell-MM synapses become resistant to killing (1/h)") # Table S2 'kmut_SYN' = 14.2651, calibrated to the cytotoxicity assay
    kkillTRGT_CD8 <- fixed(1.17e-4); label("CD38+ PBMC killing rate by an active CD8+ T cell in synapse (1/h)") # Table S2 'kkillTRGT_CD8' = 1.17e-4, calibrated to the MIMIC assay
    kkillTRGT_CD4 <- fixed(3.71e-6); label("CD38+ PBMC killing rate by an active CD4+ T cell in synapse (1/h)") # Table S2 'kkillTRGT_CD4' = 3.71e-6, calibrated to the MIMIC assay

    # Cytokine production (molecules/cell/h); MM-cell rates equal the CD38+ PBMC rates
    kprod_CD8_A_IFNg <- fixed(1.56e3); label("Production rate of IFN-gamma by active CD8+ T cells (molecules/cell/h)") # Table S2 'kprod_CD8_A_IFNg' = 1.56e3, calibrated to the MIMIC assay
    kprod_CD4_A_IFNg <- fixed(1.32e4); label("Production rate of IFN-gamma by active CD4+ T cells (molecules/cell/h)") # Table S2 'kprod_CD4_A_IFNg' = 1.32e4, calibrated to the MIMIC assay
    kprod_TRGT_IFNg <- fixed(0); label("Production rate of IFN-gamma by CD38+ PBMCs and MM cells in synapse (molecules/cell/h)") # Table S2 'kprod_TRGT_IFNg' = 0, 'Assume PBMCs only produce IL-6 and IL-10'
    kprod_CD8_A_TNFa <- fixed(147.39); label("Production rate of TNF-alpha by active CD8+ T cells (molecules/cell/h)") # Table S2 'kprod_CD8_A_TNFa' = 147.39, calibrated to the MIMIC assay
    kprod_CD4_A_TNFa <- fixed(1.94e3); label("Production rate of TNF-alpha by active CD4+ T cells (molecules/cell/h)") # Table S2 'kprod_CD4_A_TNFa' = 1.94e3, calibrated to the MIMIC assay
    kprod_TRGT_TNFa <- fixed(0); label("Production rate of TNF-alpha by CD38+ PBMCs and MM cells in synapse (molecules/cell/h)") # Table S2 'kprod_TRGT_TNFa' = 0, 'Assume PBMCs only produce IL-6 and IL-10'
    kprod_CD8_A_IL6 <- fixed(896.169); label("Production rate of IL-6 by active CD8+ T cells (molecules/cell/h)") # Table S2 'kprod_CD8_A_IL6' = 896.169, calibrated to the MIMIC assay
    kprod_CD4_A_IL6 <- fixed(566.361); label("Production rate of IL-6 by active CD4+ T cells (molecules/cell/h)") # Table S2 'kprod_CD4_A_IL6' = 566.361, calibrated to the MIMIC assay
    kprod_TRGT_IL6 <- fixed(3.07e4); label("Production rate of IL-6 by CD38+ PBMCs and MM cells in synapse (molecules/cell/h)") # Table S2 'kprod_TRGT_IL6' = 3.07e4, calibrated to the MIMIC assay
    kprod_CD8_A_IL10 <- fixed(0); label("Production rate of IL-10 by active CD8+ T cells (molecules/cell/h)") # Table S2 'kprod_CD8_A_IL10' = 0, 'Assume T-cells do not produce IL-10'
    kprod_CD4_A_IL10 <- fixed(0); label("Production rate of IL-10 by active CD4+ T cells (molecules/cell/h)") # Table S2 'kprod_CD4_A_IL10' = 0, 'Assume T-cells do not produce IL-10'
    kprod_TRGT_IL10 <- fixed(4.42e4); label("Production rate of IL-10 by CD38+ PBMCs and MM cells in synapse (molecules/cell/h)") # Table S2 'kprod_TRGT_IL10' = 4.42e4, calibrated to the MIMIC assay
    kprod_CD8_N_IFNg <- fixed(0); label("Production rate of IFN-gamma by naive CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_EM_IFNg <- fixed(0); label("Production rate of IFN-gamma by effector memory CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_N_IFNg <- fixed(0); label("Production rate of IFN-gamma by naive CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_EM_IFNg <- fixed(0); label("Production rate of IFN-gamma by effector memory CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_N_TNFa <- fixed(0); label("Production rate of TNF-alpha by naive CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_EM_TNFa <- fixed(0); label("Production rate of TNF-alpha by effector memory CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_N_TNFa <- fixed(0); label("Production rate of TNF-alpha by naive CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_EM_TNFa <- fixed(0); label("Production rate of TNF-alpha by effector memory CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_N_IL6 <- fixed(0); label("Production rate of IL-6 by naive CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_EM_IL6 <- fixed(0); label("Production rate of IL-6 by effector memory CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_N_IL6 <- fixed(0); label("Production rate of IL-6 by naive CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_EM_IL6 <- fixed(0); label("Production rate of IL-6 by effector memory CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_N_IL10 <- fixed(0); label("Production rate of IL-10 by naive CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD8_EM_IL10 <- fixed(0); label("Production rate of IL-10 by effector memory CD8+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_N_IL10 <- fixed(0); label("Production rate of IL-10 by naive CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kprod_CD4_EM_IL10 <- fixed(0); label("Production rate of IL-10 by effector memory CD4+ T cells in synapse (molecules/cell/h)") # not listed in Table S2; zero per Methods 'active T-cells produce TNF-alpha, IFN-gamma, and IL-6'
    kdeg_IFNg <- fixed(0); label("Degradation rate of IFN-gamma (1/h)") # Table S2 'kdeg_IFNg/TNFa/IL6/IL10' = 0, no significant in vitro degradation
    kdeg_TNFa <- fixed(0); label("Degradation rate of TNF-alpha (1/h)") # Table S2 'kdeg_IFNg/TNFa/IL6/IL10' = 0, no significant in vitro degradation
    kdeg_IL6 <- fixed(0); label("Degradation rate of IL-6 (1/h)") # Table S2 'kdeg_IFNg/TNFa/IL6/IL10' = 0, no significant in vitro degradation
    kdeg_IL10 <- fixed(0); label("Degradation rate of IL-10 (1/h)") # Table S2 'kdeg_IFNg/TNFa/IL6/IL10' = 0, no significant in vitro degradation

    # Soluble CD38
    kshedMM_s38 <- fixed(0.005); label("Shedding rate of soluble CD38 per MM-cell CD38 receptor (1/h)") # Table S2 'kshedMM_s38' = 0.005, calibrated to the cytotoxicity assay
    kshedTRGT_s38 <- fixed(0); label("Shedding rate of soluble CD38 per CD38+ PBMC CD38 receptor (1/h)") # Table S2 'kshedTRGT_s38' = 0, only MM shedding assumed significant
    kdeg_sCD38 <- fixed(0); label("Degradation rate of soluble CD38 (1/h)") # Table S2 'kdeg_sCD38' = 0, no significant in vitro degradation

    # T-cell activation
    kact_N <- fixed(101); label("Maximum activation rate of naive T cells in synapse (1/h)") # Table S2 'kact_N' = 101, calibrated to the activation assay
    kact_EM <- fixed(2.14); label("Activation rate of effector memory T cells in synapse (1/h)") # Table S2 'kact_EM' = 2.14, calibrated to the activation assay
    EC50_CD3 <- fixed(0.0169); label("Fraction of synapse CD3 bound by drug giving half-maximal naive activation (fraction)") # Table S2 'EC50_CD3' = 0.0169, calibrated to the activation assay
    EC50_CD28 <- fixed(0.0137); label("Fraction of synapse CD28 bound by drug giving half-maximal naive activation (fraction)") # Table S2 'EC50_CD28' = 0.0137, calibrated to the activation assay
    E <- fixed(1); label("Small offset preventing division by zero in the naive activation term (molecules)") # Table S2 'E' = 1, 'Set to small error value'

    # Antigen densities
    CD3per_CD4 <- fixed(124000); label("CD3 antigen density on CD4+ T cells, all states (molecules/cell)") # Table S2 'CD3per_CD4 (A, EM, N)' = 124000, literature [62]
    CD3per_CD8 <- fixed(124000); label("CD3 antigen density on CD8+ T cells, all states (molecules/cell)") # Table S2 'CD3per_CD8 (A, EM, N)' = 124000, literature [62]
    CD28per_CD4 <- fixed(19000); label("CD28 antigen density on CD4+ T cells, all states (molecules/cell)") # Table S2 'CD28per_CD4 (A, EM, N)' = 19000, literature [51]
    CD28per_CD8 <- fixed(12500); label("CD28 antigen density on CD8+ T cells, all states (molecules/cell)") # Table S2 'CD28per_CD8 (A, EM, N)' = 12500, literature [51]
    CD28per_MM <- fixed(50000); label("CD28 antigen density on MM cells (molecules/cell)") # Table S2 'CD28per_MM' = 50000, mean of KMS-11 and RPMI-8226
    CD38per_MM <- fixed(23000); label("CD38 antigen density on MM cells (molecules/cell)") # Table S2 'CD38per_MM' = 23000, literature [50, 66]
    CD28per_TRGT <- fixed(0); label("CD28 antigen density on CD38+ PBMCs (molecules/cell)") # Table S2 'CD28per_TRGT' = 0, no significant CD28 on CD38+ PBMCs
    CD38per_TRGT <- fixed(3.32e3); label("CD38 antigen density on CD38+ PBMCs (molecules/cell)") # Table S2 'CD38per_TRGT' = 3.32e3, B/NK/monocyte mean with daratumumab downregulation

    # Drug binding (per molecule; second-order rates already scaled to the well volume)
    kon_CD3 <- fixed(6.7e-13); label("Association rate of drug with CD3 (1/(molecule*h))") # Table S2 'kon_CD3' = 6.7e-13, internal data [9]
    koff_CD3 <- fixed(1.6524); label("Dissociation rate of drug from CD3 (1/h)") # Table S2 'koff_CD3' = 1.6524, internal data [9]
    kon_CD28 <- fixed(2.8e-12); label("Association rate of drug with CD28 (1/(molecule*h))") # Table S2 'kon_CD28' = 2.8e-12, internal data [9]; 0 gives the CD38xCD3 bispecific
    koff_CD28 <- fixed(0.7344); label("Dissociation rate of drug from CD28 (1/h)") # Table S2 'koff_CD28' = 0.7344, internal data [9]
    kon_CD38 <- fixed(1.36e-11); label("Association rate of drug with CD38 (1/(molecule*h))") # Table S2 'kon_CD38' = 1.36e-11, internal data [9]
    koff_CD38 <- fixed(7.2); label("Dissociation rate of drug from CD38 (1/h)") # Table S2 'koff_CD38' = 7.2, internal data [9]
    Vol <- fixed(2.00e-4); label("In vitro well volume used to convert drug molecules to nM (L)") # Table S2 'Vol' = 2.00e-4

    # Synapse formation
    kcoll <- fixed(82.7); label("Collision factor, propensity of two cells to collide (unitless)") # Table S2 'kcoll' = 82.7, calibrated to the activation assay
    Br_perS <- fixed(14.2); label("Number of receptor bridges per synapse (unitless)") # Table S2 'Br_perS' = 14.2, calibrated to the activation assay
    kDis <- fixed(2.68e-4); label("Dissociation rate of synapses that can dissociate (1/h)") # Table S2 'kDis' = 2.68e-4, calibrated to the activation assay
    CD8_N_0 <- fixed(20274); label("Scaling cell number for CD8_N in the synapse-formation term (cells)") # Table S2 'CD8_N_0' = 20274, 'Set to initial cell number'
    CD8_EM_0 <- fixed(20274); label("Scaling cell number for CD8_EM in the synapse-formation term (cells)") # Table S2 'CD8_EM_0' = 20274, 'Set to initial cell number'
    CD8_A_0 <- fixed(20274); label("Scaling cell number for CD8_A in the synapse-formation term (cells)") # Table S2 'CD8_A_0' = 20274, 'Set to initial cell number'
    CD4_N_0 <- fixed(10488); label("Scaling cell number for CD4_N in the synapse-formation term (cells)") # Table S2 'CD4_N_0' = 10488, 'Set to initial cell number'
    CD4_EM_0 <- fixed(10488); label("Scaling cell number for CD4_EM in the synapse-formation term (cells)") # Table S2 'CD4_EM_0' = 10488, 'Set to initial cell number'
    CD4_A_0 <- fixed(10488); label("Scaling cell number for CD4_A in the synapse-formation term (cells)") # Table S2 'CD4_A_0' = 10488, 'Set to initial cell number'
    TRGT_0 <- fixed(4960000); label("Scaling cell number for TRGT in the synapse-formation term (cells)") # Table S2 'TRGT_0' = 4960000, 'Set to initial cell number'
    MM_0 <- fixed(2540800); label("Scaling cell number for MM in the synapse-formation term (cells)") # Table S2 'MM_0' = 2540800, 'Set to initial cell number'

    # Initial cell numbers (dose-prediction module: peripheral blood of an MM patient)
    bl_CD8_N <- fixed(1.55e4); label("Initial number of naive CD8+ T cells (cells)") # Table S2 initial value 'CD8_N' = 1.55e4, [39, 67]
    bl_CD8_EM <- fixed(8.12e3); label("Initial number of effector memory CD8+ T cells (cells)") # Table S2 initial value 'CD8_EM' = 8.12e3, [39, 67], MIMIC assay
    bl_CD8_A <- fixed(2.36e2); label("Initial number of active CD8+ T cells (cells)") # Table S2 initial value 'CD8_A' = 2.36e2, [39, 67], MIMIC assay
    bl_CD4_N <- fixed(1.20e4); label("Initial number of naive CD4+ T cells (cells)") # Table S2 initial value 'CD4_N' = 1.20e4, [39, 67]
    bl_CD4_EM <- fixed(3.11e3); label("Initial number of effector memory CD4+ T cells (cells)") # Table S2 initial value 'CD4_EM' = 3.11e3, [39, 67], MIMIC assay
    bl_CD4_A <- fixed(2.11e2); label("Initial number of active CD4+ T cells (cells)") # Table S2 initial value 'CD4_A' = 2.11e2, [39, 67], MIMIC assay
    bl_MM <- fixed(3.92e3); label("Initial number of multiple myeloma cells (cells)") # Table S2 initial value 'MM' = 3.92e3, [79]
    bl_TRGT <- fixed(2.90e4); label("Initial number of CD38+ PBMCs (cells)") # Table S2 initial value 'TRGT' = 2.90e4, [39, 67]
    bl_sCD38 <- fixed(4.96e6); label("Initial number of soluble CD38 molecules (molecules)") # Table S2 initial value 'sCD38' = 4.96e6, literature [80]
  })

  model({
    # Synapse-formation flux per cell pair (Figure 1C konOFF template):
    # kcoll/(CELL1_0*CELL2_0)*(konR1*fR1_CELL1*bR2_CELL2 + konR2*fR2_CELL2*bR1_CELL1 -
    #   (koffR1 + koffR2)*BR_SYN). Multiplied by CELL1*CELL2 in the ODEs below.
    kf_CD8_N_CD3_CD8_N_CD28 <- kcoll/(CD8_N_0*CD8_N_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD8_N_CD28_Br)
    kf_CD8_N_CD3_CD8_EM_CD28 <- kcoll/(CD8_N_0*CD8_EM_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_EM_CD28_tsAb + kon_CD28*
      R_CD8_EM_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD8_EM_CD28_Br)
    kf_CD8_N_CD3_CD8_A_CD28 <- kcoll/(CD8_N_0*CD8_A_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD8_A_CD28_Br)
    kf_CD8_N_CD3_CD4_N_CD28 <- kcoll/(CD8_N_0*CD4_N_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD4_N_CD28_Br)
    kf_CD8_N_CD3_CD4_EM_CD28 <- kcoll/(CD8_N_0*CD4_EM_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_EM_CD28_tsAb + kon_CD28*
      R_CD4_EM_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD4_EM_CD28_Br)
    kf_CD8_N_CD3_CD4_A_CD28 <- kcoll/(CD8_N_0*CD4_A_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_CD4_A_CD28_Br)
    kf_CD8_N_CD3_MM_CD28 <- kcoll/(CD8_N_0*MM_0)*(kon_CD3*R_CD8_N_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_N_CD3_MM_CD28_Br)
    kf_CD8_EM_CD3_CD8_N_CD28 <- kcoll/(CD8_EM_0*CD8_N_0)*(kon_CD3*R_CD8_EM_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD8_N_CD28_Br)
    kf_CD8_EM_CD3_CD8_EM_CD28 <- kcoll/(CD8_EM_0*CD8_EM_0)*(kon_CD3*R_CD8_EM_CD3*R_CD8_EM_CD28_tsAb +
       kon_CD28*R_CD8_EM_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD8_EM_CD28_Br)
    kf_CD8_EM_CD3_CD8_A_CD28 <- kcoll/(CD8_EM_0*CD8_A_0)*(kon_CD3*R_CD8_EM_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD8_A_CD28_Br)
    kf_CD8_EM_CD3_CD4_N_CD28 <- kcoll/(CD8_EM_0*CD4_N_0)*(kon_CD3*R_CD8_EM_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD4_N_CD28_Br)
    kf_CD8_EM_CD3_CD4_EM_CD28 <- kcoll/(CD8_EM_0*CD4_EM_0)*(kon_CD3*R_CD8_EM_CD3*R_CD4_EM_CD28_tsAb +
       kon_CD28*R_CD4_EM_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD4_EM_CD28_Br)
    kf_CD8_EM_CD3_CD4_A_CD28 <- kcoll/(CD8_EM_0*CD4_A_0)*(kon_CD3*R_CD8_EM_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_CD4_A_CD28_Br)
    kf_CD8_EM_CD3_MM_CD28 <- kcoll/(CD8_EM_0*MM_0)*(kon_CD3*R_CD8_EM_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_EM_CD3_MM_CD28_Br)
    kf_CD8_A_CD3_CD8_N_CD28 <- kcoll/(CD8_A_0*CD8_N_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD8_N_CD28_Br)
    kf_CD8_A_CD3_CD8_EM_CD28 <- kcoll/(CD8_A_0*CD8_EM_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_EM_CD28_tsAb + kon_CD28*
      R_CD8_EM_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD8_EM_CD28_Br)
    kf_CD8_A_CD3_CD8_A_CD28 <- kcoll/(CD8_A_0*CD8_A_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD8_A_CD28_Br)
    kf_CD8_A_CD3_CD4_N_CD28 <- kcoll/(CD8_A_0*CD4_N_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD4_N_CD28_Br)
    kf_CD8_A_CD3_CD4_EM_CD28 <- kcoll/(CD8_A_0*CD4_EM_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_EM_CD28_tsAb + kon_CD28*
      R_CD4_EM_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD4_EM_CD28_Br)
    kf_CD8_A_CD3_CD4_A_CD28 <- kcoll/(CD8_A_0*CD4_A_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_CD4_A_CD28_Br)
    kf_CD8_A_CD3_MM_CD28 <- kcoll/(CD8_A_0*MM_0)*(kon_CD3*R_CD8_A_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD8_A_CD3_MM_CD28_Br)
    kf_CD4_N_CD3_CD8_N_CD28 <- kcoll/(CD4_N_0*CD8_N_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD8_N_CD28_Br)
    kf_CD4_N_CD3_CD8_EM_CD28 <- kcoll/(CD4_N_0*CD8_EM_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_EM_CD28_tsAb + kon_CD28*
      R_CD8_EM_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD8_EM_CD28_Br)
    kf_CD4_N_CD3_CD8_A_CD28 <- kcoll/(CD4_N_0*CD8_A_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD8_A_CD28_Br)
    kf_CD4_N_CD3_CD4_N_CD28 <- kcoll/(CD4_N_0*CD4_N_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD4_N_CD28_Br)
    kf_CD4_N_CD3_CD4_EM_CD28 <- kcoll/(CD4_N_0*CD4_EM_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_EM_CD28_tsAb + kon_CD28*
      R_CD4_EM_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD4_EM_CD28_Br)
    kf_CD4_N_CD3_CD4_A_CD28 <- kcoll/(CD4_N_0*CD4_A_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_CD4_A_CD28_Br)
    kf_CD4_N_CD3_MM_CD28 <- kcoll/(CD4_N_0*MM_0)*(kon_CD3*R_CD4_N_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_N_CD3_MM_CD28_Br)
    kf_CD4_EM_CD3_CD8_N_CD28 <- kcoll/(CD4_EM_0*CD8_N_0)*(kon_CD3*R_CD4_EM_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD8_N_CD28_Br)
    kf_CD4_EM_CD3_CD8_EM_CD28 <- kcoll/(CD4_EM_0*CD8_EM_0)*(kon_CD3*R_CD4_EM_CD3*R_CD8_EM_CD28_tsAb +
       kon_CD28*R_CD8_EM_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD8_EM_CD28_Br)
    kf_CD4_EM_CD3_CD8_A_CD28 <- kcoll/(CD4_EM_0*CD8_A_0)*(kon_CD3*R_CD4_EM_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD8_A_CD28_Br)
    kf_CD4_EM_CD3_CD4_N_CD28 <- kcoll/(CD4_EM_0*CD4_N_0)*(kon_CD3*R_CD4_EM_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD4_N_CD28_Br)
    kf_CD4_EM_CD3_CD4_EM_CD28 <- kcoll/(CD4_EM_0*CD4_EM_0)*(kon_CD3*R_CD4_EM_CD3*R_CD4_EM_CD28_tsAb +
       kon_CD28*R_CD4_EM_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD4_EM_CD28_Br)
    kf_CD4_EM_CD3_CD4_A_CD28 <- kcoll/(CD4_EM_0*CD4_A_0)*(kon_CD3*R_CD4_EM_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_CD4_A_CD28_Br)
    kf_CD4_EM_CD3_MM_CD28 <- kcoll/(CD4_EM_0*MM_0)*(kon_CD3*R_CD4_EM_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_EM_CD3_MM_CD28_Br)
    kf_CD4_A_CD3_CD8_N_CD28 <- kcoll/(CD4_A_0*CD8_N_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_N_CD28_tsAb + kon_CD28*
      R_CD8_N_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD8_N_CD28_Br)
    kf_CD4_A_CD3_CD8_EM_CD28 <- kcoll/(CD4_A_0*CD8_EM_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_EM_CD28_tsAb + kon_CD28*
      R_CD8_EM_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD8_EM_CD28_Br)
    kf_CD4_A_CD3_CD8_A_CD28 <- kcoll/(CD4_A_0*CD8_A_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_A_CD28_tsAb + kon_CD28*
      R_CD8_A_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD8_A_CD28_Br)
    kf_CD4_A_CD3_CD4_N_CD28 <- kcoll/(CD4_A_0*CD4_N_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_N_CD28_tsAb + kon_CD28*
      R_CD4_N_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD4_N_CD28_Br)
    kf_CD4_A_CD3_CD4_EM_CD28 <- kcoll/(CD4_A_0*CD4_EM_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_EM_CD28_tsAb + kon_CD28*
      R_CD4_EM_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD4_EM_CD28_Br)
    kf_CD4_A_CD3_CD4_A_CD28 <- kcoll/(CD4_A_0*CD4_A_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_A_CD28_tsAb + kon_CD28*
      R_CD4_A_CD28*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_CD4_A_CD28_Br)
    kf_CD4_A_CD3_MM_CD28 <- kcoll/(CD4_A_0*MM_0)*(kon_CD3*R_CD4_A_CD3*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*
      R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD28)*S_CD4_A_CD3_MM_CD28_Br)
    kf_CD8_N_CD3_TRGT_CD38 <- kcoll/(CD8_N_0*TRGT_0)*(kon_CD3*R_CD8_N_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_N_CD3_TRGT_CD38_Br)
    kf_CD8_N_CD3_MM_CD38 <- kcoll/(CD8_N_0*MM_0)*(kon_CD3*R_CD8_N_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD8_N_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_N_CD3_MM_CD38_Br)
    kf_CD8_EM_CD3_TRGT_CD38 <- kcoll/(CD8_EM_0*TRGT_0)*(kon_CD3*R_CD8_EM_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_EM_CD3_TRGT_CD38_Br)
    kf_CD8_EM_CD3_MM_CD38 <- kcoll/(CD8_EM_0*MM_0)*(kon_CD3*R_CD8_EM_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD8_EM_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_EM_CD3_MM_CD38_Br)
    kf_CD8_A_CD3_TRGT_CD38 <- kcoll/(CD8_A_0*TRGT_0)*(kon_CD3*R_CD8_A_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_A_CD3_TRGT_CD38_Br)
    kf_CD8_A_CD3_MM_CD38 <- kcoll/(CD8_A_0*MM_0)*(kon_CD3*R_CD8_A_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD8_A_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD8_A_CD3_MM_CD38_Br)
    kf_CD4_N_CD3_TRGT_CD38 <- kcoll/(CD4_N_0*TRGT_0)*(kon_CD3*R_CD4_N_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_N_CD3_TRGT_CD38_Br)
    kf_CD4_N_CD3_MM_CD38 <- kcoll/(CD4_N_0*MM_0)*(kon_CD3*R_CD4_N_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD4_N_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_N_CD3_MM_CD38_Br)
    kf_CD4_EM_CD3_TRGT_CD38 <- kcoll/(CD4_EM_0*TRGT_0)*(kon_CD3*R_CD4_EM_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_EM_CD3_TRGT_CD38_Br)
    kf_CD4_EM_CD3_MM_CD38 <- kcoll/(CD4_EM_0*MM_0)*(kon_CD3*R_CD4_EM_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD4_EM_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_EM_CD3_MM_CD38_Br)
    kf_CD4_A_CD3_TRGT_CD38 <- kcoll/(CD4_A_0*TRGT_0)*(kon_CD3*R_CD4_A_CD3*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_A_CD3_TRGT_CD38_Br)
    kf_CD4_A_CD3_MM_CD38 <- kcoll/(CD4_A_0*MM_0)*(kon_CD3*R_CD4_A_CD3*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD4_A_CD3_tsAb - (koff_CD3 + koff_CD38)*S_CD4_A_CD3_MM_CD38_Br)
    kf_CD8_N_CD28_TRGT_CD38 <- kcoll/(CD8_N_0*TRGT_0)*(kon_CD28*R_CD8_N_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_N_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_N_CD28_TRGT_CD38_Br)
    kf_CD8_N_CD28_MM_CD38 <- kcoll/(CD8_N_0*MM_0)*(kon_CD28*R_CD8_N_CD28*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD8_N_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_N_CD28_MM_CD38_Br)
    kf_CD8_EM_CD28_TRGT_CD38 <- kcoll/(CD8_EM_0*TRGT_0)*(kon_CD28*R_CD8_EM_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_EM_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_EM_CD28_TRGT_CD38_Br)
    kf_CD8_EM_CD28_MM_CD38 <- kcoll/(CD8_EM_0*MM_0)*(kon_CD28*R_CD8_EM_CD28*R_MM_CD38_tsAb + kon_CD38*
      R_MM_CD38*R_CD8_EM_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_EM_CD28_MM_CD38_Br)
    kf_CD8_A_CD28_TRGT_CD38 <- kcoll/(CD8_A_0*TRGT_0)*(kon_CD28*R_CD8_A_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD8_A_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_A_CD28_TRGT_CD38_Br)
    kf_CD8_A_CD28_MM_CD38 <- kcoll/(CD8_A_0*MM_0)*(kon_CD28*R_CD8_A_CD28*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD8_A_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD8_A_CD28_MM_CD38_Br)
    kf_CD4_N_CD28_TRGT_CD38 <- kcoll/(CD4_N_0*TRGT_0)*(kon_CD28*R_CD4_N_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_N_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_N_CD28_TRGT_CD38_Br)
    kf_CD4_N_CD28_MM_CD38 <- kcoll/(CD4_N_0*MM_0)*(kon_CD28*R_CD4_N_CD28*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD4_N_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_N_CD28_MM_CD38_Br)
    kf_CD4_EM_CD28_TRGT_CD38 <- kcoll/(CD4_EM_0*TRGT_0)*(kon_CD28*R_CD4_EM_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_EM_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_EM_CD28_TRGT_CD38_Br)
    kf_CD4_EM_CD28_MM_CD38 <- kcoll/(CD4_EM_0*MM_0)*(kon_CD28*R_CD4_EM_CD28*R_MM_CD38_tsAb + kon_CD38*
      R_MM_CD38*R_CD4_EM_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_EM_CD28_MM_CD38_Br)
    kf_CD4_A_CD28_TRGT_CD38 <- kcoll/(CD4_A_0*TRGT_0)*(kon_CD28*R_CD4_A_CD28*R_TRGT_CD38_tsAb + kon_CD38*
      R_TRGT_CD38*R_CD4_A_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_A_CD28_TRGT_CD38_Br)
    kf_CD4_A_CD28_MM_CD38 <- kcoll/(CD4_A_0*MM_0)*(kon_CD28*R_CD4_A_CD28*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*
      R_CD4_A_CD28_tsAb - (koff_CD28 + koff_CD38)*S_CD4_A_CD28_MM_CD38_Br)
    kf_MM_CD28_TRGT_CD38 <- kcoll/(MM_0*TRGT_0)*(kon_CD28*R_MM_CD28*R_TRGT_CD38_tsAb + kon_CD38*R_TRGT_CD38*
      R_MM_CD28_tsAb - (koff_CD28 + koff_CD38)*S_MM_CD28_TRGT_CD38_Br)
    kf_MM_CD28_MM_CD38 <- kcoll/(MM_0*MM_0)*(kon_CD28*R_MM_CD28*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*R_MM_CD28_tsAb -
       (koff_CD28 + koff_CD38)*S_MM_CD28_MM_CD38_Br)

    # Bridge formation through each of the two drug arms of a synapse (generator
    # terms kbr_R1 / kbr_R2), scaled by the bridges-per-synapse factor Br_perS.
    kb1_CD8_N_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_N*CD8_N/(CD8_N_0*CD8_N_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_N_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD8_N_CD28_Br)
    kb1_CD8_N_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_N*CD8_EM/(CD8_N_0*CD8_EM_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_EM_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD8_EM_CD28_Br)
    kb1_CD8_N_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_N*CD8_A/(CD8_N_0*CD8_A_0)*(kon_CD3*R_CD8_N_CD3*R_CD8_A_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD8_A_CD28_Br)
    kb1_CD8_N_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_N*CD4_N/(CD8_N_0*CD4_N_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_N_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD4_N_CD28_Br)
    kb1_CD8_N_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_N*CD4_EM/(CD8_N_0*CD4_EM_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_EM_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD4_EM_CD28_Br)
    kb1_CD8_N_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_N*CD4_A/(CD8_N_0*CD4_A_0)*(kon_CD3*R_CD8_N_CD3*R_CD4_A_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_CD4_A_CD28_Br)
    kb1_CD8_N_CD3_MM_CD28 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD3*R_CD8_N_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD8_N_CD3_MM_CD28_Br)
    kb1_CD8_EM_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_EM*CD8_N/(CD8_EM_0*CD8_N_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD8_N_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD8_N_CD28_Br)
    kb1_CD8_EM_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_EM*CD8_EM/(CD8_EM_0*CD8_EM_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD8_EM_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD8_EM_CD28_Br)
    kb1_CD8_EM_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_EM*CD8_A/(CD8_EM_0*CD8_A_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD8_A_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD8_A_CD28_Br)
    kb1_CD8_EM_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_EM*CD4_N/(CD8_EM_0*CD4_N_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD4_N_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD4_N_CD28_Br)
    kb1_CD8_EM_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_EM*CD4_EM/(CD8_EM_0*CD4_EM_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD4_EM_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD4_EM_CD28_Br)
    kb1_CD8_EM_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_EM*CD4_A/(CD8_EM_0*CD4_A_0)*(kon_CD3*R_CD8_EM_CD3*
      R_CD4_A_CD28_tsAb - koff_CD3*S_CD8_EM_CD3_CD4_A_CD28_Br)
    kb1_CD8_EM_CD3_MM_CD28 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD3*R_CD8_EM_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD8_EM_CD3_MM_CD28_Br)
    kb1_CD8_A_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_A*CD8_N/(CD8_A_0*CD8_N_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_N_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD8_N_CD28_Br)
    kb1_CD8_A_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_A*CD8_EM/(CD8_A_0*CD8_EM_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_EM_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_Br)
    kb1_CD8_A_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_A*CD8_A/(CD8_A_0*CD8_A_0)*(kon_CD3*R_CD8_A_CD3*R_CD8_A_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD8_A_CD28_Br)
    kb1_CD8_A_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_A*CD4_N/(CD8_A_0*CD4_N_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_N_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD4_N_CD28_Br)
    kb1_CD8_A_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_A*CD4_EM/(CD8_A_0*CD4_EM_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_EM_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_Br)
    kb1_CD8_A_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_A*CD4_A/(CD8_A_0*CD4_A_0)*(kon_CD3*R_CD8_A_CD3*R_CD4_A_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_CD4_A_CD28_Br)
    kb1_CD8_A_CD3_MM_CD28 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD3*R_CD8_A_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD8_A_CD3_MM_CD28_Br)
    kb1_CD4_N_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_N*CD8_N/(CD4_N_0*CD8_N_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_N_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD8_N_CD28_Br)
    kb1_CD4_N_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_N*CD8_EM/(CD4_N_0*CD8_EM_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_EM_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD8_EM_CD28_Br)
    kb1_CD4_N_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_N*CD8_A/(CD4_N_0*CD8_A_0)*(kon_CD3*R_CD4_N_CD3*R_CD8_A_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD8_A_CD28_Br)
    kb1_CD4_N_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_N*CD4_N/(CD4_N_0*CD4_N_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_N_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD4_N_CD28_Br)
    kb1_CD4_N_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_N*CD4_EM/(CD4_N_0*CD4_EM_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_EM_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD4_EM_CD28_Br)
    kb1_CD4_N_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_N*CD4_A/(CD4_N_0*CD4_A_0)*(kon_CD3*R_CD4_N_CD3*R_CD4_A_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_CD4_A_CD28_Br)
    kb1_CD4_N_CD3_MM_CD28 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD3*R_CD4_N_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD4_N_CD3_MM_CD28_Br)
    kb1_CD4_EM_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_EM*CD8_N/(CD4_EM_0*CD8_N_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD8_N_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD8_N_CD28_Br)
    kb1_CD4_EM_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_EM*CD8_EM/(CD4_EM_0*CD8_EM_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD8_EM_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD8_EM_CD28_Br)
    kb1_CD4_EM_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_EM*CD8_A/(CD4_EM_0*CD8_A_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD8_A_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD8_A_CD28_Br)
    kb1_CD4_EM_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_EM*CD4_N/(CD4_EM_0*CD4_N_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD4_N_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD4_N_CD28_Br)
    kb1_CD4_EM_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_EM*CD4_EM/(CD4_EM_0*CD4_EM_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD4_EM_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD4_EM_CD28_Br)
    kb1_CD4_EM_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_EM*CD4_A/(CD4_EM_0*CD4_A_0)*(kon_CD3*R_CD4_EM_CD3*
      R_CD4_A_CD28_tsAb - koff_CD3*S_CD4_EM_CD3_CD4_A_CD28_Br)
    kb1_CD4_EM_CD3_MM_CD28 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD3*R_CD4_EM_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD4_EM_CD3_MM_CD28_Br)
    kb1_CD4_A_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_A*CD8_N/(CD4_A_0*CD8_N_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_N_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD8_N_CD28_Br)
    kb1_CD4_A_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_A*CD8_EM/(CD4_A_0*CD8_EM_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_EM_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_Br)
    kb1_CD4_A_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_A*CD8_A/(CD4_A_0*CD8_A_0)*(kon_CD3*R_CD4_A_CD3*R_CD8_A_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD8_A_CD28_Br)
    kb1_CD4_A_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_A*CD4_N/(CD4_A_0*CD4_N_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_N_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD4_N_CD28_Br)
    kb1_CD4_A_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_A*CD4_EM/(CD4_A_0*CD4_EM_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_EM_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_Br)
    kb1_CD4_A_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_A*CD4_A/(CD4_A_0*CD4_A_0)*(kon_CD3*R_CD4_A_CD3*R_CD4_A_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_CD4_A_CD28_Br)
    kb1_CD4_A_CD3_MM_CD28 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD3*R_CD4_A_CD3*R_MM_CD28_tsAb -
       koff_CD3*S_CD4_A_CD3_MM_CD28_Br)
    kb1_CD8_N_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_N*TRGT/(CD8_N_0*TRGT_0)*(kon_CD3*R_CD8_N_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD8_N_CD3_TRGT_CD38_Br)
    kb1_CD8_N_CD3_MM_CD38 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD3*R_CD8_N_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD8_N_CD3_MM_CD38_Br)
    kb1_CD8_EM_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_EM*TRGT/(CD8_EM_0*TRGT_0)*(kon_CD3*R_CD8_EM_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD8_EM_CD3_TRGT_CD38_Br)
    kb1_CD8_EM_CD3_MM_CD38 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD3*R_CD8_EM_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD8_EM_CD3_MM_CD38_Br)
    kb1_CD8_A_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_A*TRGT/(CD8_A_0*TRGT_0)*(kon_CD3*R_CD8_A_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD8_A_CD3_TRGT_CD38_Br)
    kb1_CD8_A_CD3_MM_CD38 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD3*R_CD8_A_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD8_A_CD3_MM_CD38_Br)
    kb1_CD4_N_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_N*TRGT/(CD4_N_0*TRGT_0)*(kon_CD3*R_CD4_N_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD4_N_CD3_TRGT_CD38_Br)
    kb1_CD4_N_CD3_MM_CD38 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD3*R_CD4_N_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD4_N_CD3_MM_CD38_Br)
    kb1_CD4_EM_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_EM*TRGT/(CD4_EM_0*TRGT_0)*(kon_CD3*R_CD4_EM_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD4_EM_CD3_TRGT_CD38_Br)
    kb1_CD4_EM_CD3_MM_CD38 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD3*R_CD4_EM_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD4_EM_CD3_MM_CD38_Br)
    kb1_CD4_A_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_A*TRGT/(CD4_A_0*TRGT_0)*(kon_CD3*R_CD4_A_CD3*R_TRGT_CD38_tsAb -
       koff_CD3*S_CD4_A_CD3_TRGT_CD38_Br)
    kb1_CD4_A_CD3_MM_CD38 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD3*R_CD4_A_CD3*R_MM_CD38_tsAb -
       koff_CD3*S_CD4_A_CD3_MM_CD38_Br)
    kb1_CD8_N_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_N*TRGT/(CD8_N_0*TRGT_0)*(kon_CD28*R_CD8_N_CD28*R_TRGT_CD38_tsAb -
       koff_CD28*S_CD8_N_CD28_TRGT_CD38_Br)
    kb1_CD8_N_CD28_MM_CD38 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD28*R_CD8_N_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD8_N_CD28_MM_CD38_Br)
    kb1_CD8_EM_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_EM*TRGT/(CD8_EM_0*TRGT_0)*(kon_CD28*R_CD8_EM_CD28*
      R_TRGT_CD38_tsAb - koff_CD28*S_CD8_EM_CD28_TRGT_CD38_Br)
    kb1_CD8_EM_CD28_MM_CD38 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD28*R_CD8_EM_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD8_EM_CD28_MM_CD38_Br)
    kb1_CD8_A_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_A*TRGT/(CD8_A_0*TRGT_0)*(kon_CD28*R_CD8_A_CD28*R_TRGT_CD38_tsAb -
       koff_CD28*S_CD8_A_CD28_TRGT_CD38_Br)
    kb1_CD8_A_CD28_MM_CD38 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD28*R_CD8_A_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD8_A_CD28_MM_CD38_Br)
    kb1_CD4_N_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_N*TRGT/(CD4_N_0*TRGT_0)*(kon_CD28*R_CD4_N_CD28*R_TRGT_CD38_tsAb -
       koff_CD28*S_CD4_N_CD28_TRGT_CD38_Br)
    kb1_CD4_N_CD28_MM_CD38 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD28*R_CD4_N_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD4_N_CD28_MM_CD38_Br)
    kb1_CD4_EM_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_EM*TRGT/(CD4_EM_0*TRGT_0)*(kon_CD28*R_CD4_EM_CD28*
      R_TRGT_CD38_tsAb - koff_CD28*S_CD4_EM_CD28_TRGT_CD38_Br)
    kb1_CD4_EM_CD28_MM_CD38 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD28*R_CD4_EM_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD4_EM_CD28_MM_CD38_Br)
    kb1_CD4_A_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_A*TRGT/(CD4_A_0*TRGT_0)*(kon_CD28*R_CD4_A_CD28*R_TRGT_CD38_tsAb -
       koff_CD28*S_CD4_A_CD28_TRGT_CD38_Br)
    kb1_CD4_A_CD28_MM_CD38 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD28*R_CD4_A_CD28*R_MM_CD38_tsAb -
       koff_CD28*S_CD4_A_CD28_MM_CD38_Br)
    kb1_MM_CD28_TRGT_CD38 <- Br_perS*kcoll*MM*TRGT/(MM_0*TRGT_0)*(kon_CD28*R_MM_CD28*R_TRGT_CD38_tsAb -
       koff_CD28*S_MM_CD28_TRGT_CD38_Br)
    kb1_MM_CD28_MM_CD38 <- Br_perS*kcoll*MM*MM/(MM_0*MM_0)*(kon_CD28*R_MM_CD28*R_MM_CD38_tsAb - koff_CD28*
      S_MM_CD28_MM_CD38_Br)
    kb2_CD8_N_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_N*CD8_N/(CD8_N_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*R_CD8_N_CD3_tsAb -
       koff_CD28*S_CD8_N_CD3_CD8_N_CD28_Br)
    kb2_CD8_N_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_N*CD8_EM/(CD8_N_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD8_N_CD3_tsAb - koff_CD28*S_CD8_N_CD3_CD8_EM_CD28_Br)
    kb2_CD8_N_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_N*CD8_A/(CD8_N_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*R_CD8_N_CD3_tsAb -
       koff_CD28*S_CD8_N_CD3_CD8_A_CD28_Br)
    kb2_CD8_N_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_N*CD4_N/(CD8_N_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*R_CD8_N_CD3_tsAb -
       koff_CD28*S_CD8_N_CD3_CD4_N_CD28_Br)
    kb2_CD8_N_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_N*CD4_EM/(CD8_N_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD8_N_CD3_tsAb - koff_CD28*S_CD8_N_CD3_CD4_EM_CD28_Br)
    kb2_CD8_N_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_N*CD4_A/(CD8_N_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*R_CD8_N_CD3_tsAb -
       koff_CD28*S_CD8_N_CD3_CD4_A_CD28_Br)
    kb2_CD8_N_CD3_MM_CD28 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD8_N_CD3_tsAb -
       koff_CD28*S_CD8_N_CD3_MM_CD28_Br)
    kb2_CD8_EM_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_EM*CD8_N/(CD8_EM_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD8_N_CD28_Br)
    kb2_CD8_EM_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_EM*CD8_EM/(CD8_EM_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD8_EM_CD28_Br)
    kb2_CD8_EM_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_EM*CD8_A/(CD8_EM_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD8_A_CD28_Br)
    kb2_CD8_EM_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_EM*CD4_N/(CD8_EM_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD4_N_CD28_Br)
    kb2_CD8_EM_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_EM*CD4_EM/(CD8_EM_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_Br)
    kb2_CD8_EM_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_EM*CD4_A/(CD8_EM_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*
      R_CD8_EM_CD3_tsAb - koff_CD28*S_CD8_EM_CD3_CD4_A_CD28_Br)
    kb2_CD8_EM_CD3_MM_CD28 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD8_EM_CD3_tsAb -
       koff_CD28*S_CD8_EM_CD3_MM_CD28_Br)
    kb2_CD8_A_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD8_A*CD8_N/(CD8_A_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*R_CD8_A_CD3_tsAb -
       koff_CD28*S_CD8_A_CD3_CD8_N_CD28_Br)
    kb2_CD8_A_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD8_A*CD8_EM/(CD8_A_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD8_A_CD3_tsAb - koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_Br)
    kb2_CD8_A_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD8_A*CD8_A/(CD8_A_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*R_CD8_A_CD3_tsAb -
       koff_CD28*S_CD8_A_CD3_CD8_A_CD28_Br)
    kb2_CD8_A_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD8_A*CD4_N/(CD8_A_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*R_CD8_A_CD3_tsAb -
       koff_CD28*S_CD8_A_CD3_CD4_N_CD28_Br)
    kb2_CD8_A_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD8_A*CD4_EM/(CD8_A_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD8_A_CD3_tsAb - koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_Br)
    kb2_CD8_A_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD8_A*CD4_A/(CD8_A_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*R_CD8_A_CD3_tsAb -
       koff_CD28*S_CD8_A_CD3_CD4_A_CD28_Br)
    kb2_CD8_A_CD3_MM_CD28 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD8_A_CD3_tsAb -
       koff_CD28*S_CD8_A_CD3_MM_CD28_Br)
    kb2_CD4_N_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_N*CD8_N/(CD4_N_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*R_CD4_N_CD3_tsAb -
       koff_CD28*S_CD4_N_CD3_CD8_N_CD28_Br)
    kb2_CD4_N_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_N*CD8_EM/(CD4_N_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD4_N_CD3_tsAb - koff_CD28*S_CD4_N_CD3_CD8_EM_CD28_Br)
    kb2_CD4_N_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_N*CD8_A/(CD4_N_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*R_CD4_N_CD3_tsAb -
       koff_CD28*S_CD4_N_CD3_CD8_A_CD28_Br)
    kb2_CD4_N_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_N*CD4_N/(CD4_N_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*R_CD4_N_CD3_tsAb -
       koff_CD28*S_CD4_N_CD3_CD4_N_CD28_Br)
    kb2_CD4_N_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_N*CD4_EM/(CD4_N_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD4_N_CD3_tsAb - koff_CD28*S_CD4_N_CD3_CD4_EM_CD28_Br)
    kb2_CD4_N_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_N*CD4_A/(CD4_N_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*R_CD4_N_CD3_tsAb -
       koff_CD28*S_CD4_N_CD3_CD4_A_CD28_Br)
    kb2_CD4_N_CD3_MM_CD28 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD4_N_CD3_tsAb -
       koff_CD28*S_CD4_N_CD3_MM_CD28_Br)
    kb2_CD4_EM_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_EM*CD8_N/(CD4_EM_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD8_N_CD28_Br)
    kb2_CD4_EM_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_EM*CD8_EM/(CD4_EM_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_Br)
    kb2_CD4_EM_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_EM*CD8_A/(CD4_EM_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD8_A_CD28_Br)
    kb2_CD4_EM_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_EM*CD4_N/(CD4_EM_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD4_N_CD28_Br)
    kb2_CD4_EM_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_EM*CD4_EM/(CD4_EM_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD4_EM_CD28_Br)
    kb2_CD4_EM_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_EM*CD4_A/(CD4_EM_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*
      R_CD4_EM_CD3_tsAb - koff_CD28*S_CD4_EM_CD3_CD4_A_CD28_Br)
    kb2_CD4_EM_CD3_MM_CD28 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD4_EM_CD3_tsAb -
       koff_CD28*S_CD4_EM_CD3_MM_CD28_Br)
    kb2_CD4_A_CD3_CD8_N_CD28 <- Br_perS*kcoll*CD4_A*CD8_N/(CD4_A_0*CD8_N_0)*(kon_CD28*R_CD8_N_CD28*R_CD4_A_CD3_tsAb -
       koff_CD28*S_CD4_A_CD3_CD8_N_CD28_Br)
    kb2_CD4_A_CD3_CD8_EM_CD28 <- Br_perS*kcoll*CD4_A*CD8_EM/(CD4_A_0*CD8_EM_0)*(kon_CD28*R_CD8_EM_CD28*
      R_CD4_A_CD3_tsAb - koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_Br)
    kb2_CD4_A_CD3_CD8_A_CD28 <- Br_perS*kcoll*CD4_A*CD8_A/(CD4_A_0*CD8_A_0)*(kon_CD28*R_CD8_A_CD28*R_CD4_A_CD3_tsAb -
       koff_CD28*S_CD4_A_CD3_CD8_A_CD28_Br)
    kb2_CD4_A_CD3_CD4_N_CD28 <- Br_perS*kcoll*CD4_A*CD4_N/(CD4_A_0*CD4_N_0)*(kon_CD28*R_CD4_N_CD28*R_CD4_A_CD3_tsAb -
       koff_CD28*S_CD4_A_CD3_CD4_N_CD28_Br)
    kb2_CD4_A_CD3_CD4_EM_CD28 <- Br_perS*kcoll*CD4_A*CD4_EM/(CD4_A_0*CD4_EM_0)*(kon_CD28*R_CD4_EM_CD28*
      R_CD4_A_CD3_tsAb - koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_Br)
    kb2_CD4_A_CD3_CD4_A_CD28 <- Br_perS*kcoll*CD4_A*CD4_A/(CD4_A_0*CD4_A_0)*(kon_CD28*R_CD4_A_CD28*R_CD4_A_CD3_tsAb -
       koff_CD28*S_CD4_A_CD3_CD4_A_CD28_Br)
    kb2_CD4_A_CD3_MM_CD28 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD28*R_MM_CD28*R_CD4_A_CD3_tsAb -
       koff_CD28*S_CD4_A_CD3_MM_CD28_Br)
    kb2_CD8_N_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_N*TRGT/(CD8_N_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_N_CD3_tsAb -
       koff_CD38*S_CD8_N_CD3_TRGT_CD38_Br)
    kb2_CD8_N_CD3_MM_CD38 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_N_CD3_tsAb -
       koff_CD38*S_CD8_N_CD3_MM_CD38_Br)
    kb2_CD8_EM_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_EM*TRGT/(CD8_EM_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_EM_CD3_tsAb -
       koff_CD38*S_CD8_EM_CD3_TRGT_CD38_Br)
    kb2_CD8_EM_CD3_MM_CD38 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_EM_CD3_tsAb -
       koff_CD38*S_CD8_EM_CD3_MM_CD38_Br)
    kb2_CD8_A_CD3_TRGT_CD38 <- Br_perS*kcoll*CD8_A*TRGT/(CD8_A_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_A_CD3_tsAb -
       koff_CD38*S_CD8_A_CD3_TRGT_CD38_Br)
    kb2_CD8_A_CD3_MM_CD38 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_A_CD3_tsAb -
       koff_CD38*S_CD8_A_CD3_MM_CD38_Br)
    kb2_CD4_N_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_N*TRGT/(CD4_N_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_N_CD3_tsAb -
       koff_CD38*S_CD4_N_CD3_TRGT_CD38_Br)
    kb2_CD4_N_CD3_MM_CD38 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_N_CD3_tsAb -
       koff_CD38*S_CD4_N_CD3_MM_CD38_Br)
    kb2_CD4_EM_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_EM*TRGT/(CD4_EM_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_EM_CD3_tsAb -
       koff_CD38*S_CD4_EM_CD3_TRGT_CD38_Br)
    kb2_CD4_EM_CD3_MM_CD38 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_EM_CD3_tsAb -
       koff_CD38*S_CD4_EM_CD3_MM_CD38_Br)
    kb2_CD4_A_CD3_TRGT_CD38 <- Br_perS*kcoll*CD4_A*TRGT/(CD4_A_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_A_CD3_tsAb -
       koff_CD38*S_CD4_A_CD3_TRGT_CD38_Br)
    kb2_CD4_A_CD3_MM_CD38 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_A_CD3_tsAb -
       koff_CD38*S_CD4_A_CD3_MM_CD38_Br)
    kb2_CD8_N_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_N*TRGT/(CD8_N_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_N_CD28_tsAb -
       koff_CD38*S_CD8_N_CD28_TRGT_CD38_Br)
    kb2_CD8_N_CD28_MM_CD38 <- Br_perS*kcoll*CD8_N*MM/(CD8_N_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_N_CD28_tsAb -
       koff_CD38*S_CD8_N_CD28_MM_CD38_Br)
    kb2_CD8_EM_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_EM*TRGT/(CD8_EM_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_EM_CD28_tsAb -
       koff_CD38*S_CD8_EM_CD28_TRGT_CD38_Br)
    kb2_CD8_EM_CD28_MM_CD38 <- Br_perS*kcoll*CD8_EM*MM/(CD8_EM_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_EM_CD28_tsAb -
       koff_CD38*S_CD8_EM_CD28_MM_CD38_Br)
    kb2_CD8_A_CD28_TRGT_CD38 <- Br_perS*kcoll*CD8_A*TRGT/(CD8_A_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD8_A_CD28_tsAb -
       koff_CD38*S_CD8_A_CD28_TRGT_CD38_Br)
    kb2_CD8_A_CD28_MM_CD38 <- Br_perS*kcoll*CD8_A*MM/(CD8_A_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD8_A_CD28_tsAb -
       koff_CD38*S_CD8_A_CD28_MM_CD38_Br)
    kb2_CD4_N_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_N*TRGT/(CD4_N_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_N_CD28_tsAb -
       koff_CD38*S_CD4_N_CD28_TRGT_CD38_Br)
    kb2_CD4_N_CD28_MM_CD38 <- Br_perS*kcoll*CD4_N*MM/(CD4_N_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_N_CD28_tsAb -
       koff_CD38*S_CD4_N_CD28_MM_CD38_Br)
    kb2_CD4_EM_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_EM*TRGT/(CD4_EM_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_EM_CD28_tsAb -
       koff_CD38*S_CD4_EM_CD28_TRGT_CD38_Br)
    kb2_CD4_EM_CD28_MM_CD38 <- Br_perS*kcoll*CD4_EM*MM/(CD4_EM_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_EM_CD28_tsAb -
       koff_CD38*S_CD4_EM_CD28_MM_CD38_Br)
    kb2_CD4_A_CD28_TRGT_CD38 <- Br_perS*kcoll*CD4_A*TRGT/(CD4_A_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_CD4_A_CD28_tsAb -
       koff_CD38*S_CD4_A_CD28_TRGT_CD38_Br)
    kb2_CD4_A_CD28_MM_CD38 <- Br_perS*kcoll*CD4_A*MM/(CD4_A_0*MM_0)*(kon_CD38*R_MM_CD38*R_CD4_A_CD28_tsAb -
       koff_CD38*S_CD4_A_CD28_MM_CD38_Br)
    kb2_MM_CD28_TRGT_CD38 <- Br_perS*kcoll*MM*TRGT/(MM_0*TRGT_0)*(kon_CD38*R_TRGT_CD38*R_MM_CD28_tsAb -
       koff_CD38*S_MM_CD28_TRGT_CD38_Br)
    kb2_MM_CD28_MM_CD38 <- Br_perS*kcoll*MM*MM/(MM_0*MM_0)*(kon_CD38*R_MM_CD38*R_MM_CD28_tsAb - koff_CD38*
      S_MM_CD28_MM_CD38_Br)

    # Naive T-cell activation rate in synapse (Methods: AND gate of Michaelis-Menten
    # terms in drug-bound CD3 and drug-bound CD28 on the naive cell, times kact_N).
    ka_CD4_N_CD3_CD4_A_CD28 <- kact_N*((S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_A_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3 + S_CD4_N_CD3_CD4_A_CD28_Br) +
       (S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_A_CD28_Br))*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28) +
       S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_CD4_EM_CD28 <- kact_N*((S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_EM_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3 +
       S_CD4_N_CD3_CD4_EM_CD28_Br) + (S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_EM_CD28_Br))*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb/(EC50_CD28*(E + S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb +
       S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28) + S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_CD4_N_CD28 <- kact_N*((S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_N_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 + S_CD4_N_CD3_CD4_N_CD28_Br) +
       (S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD4_N_CD28_Br))*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28) +
       S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_CD8_A_CD28 <- kact_N*((S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_A_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3 + S_CD4_N_CD3_CD8_A_CD28_Br) +
       (S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_A_CD28_Br))*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28) +
       S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_CD8_EM_CD28 <- kact_N*((S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_EM_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3 +
       S_CD4_N_CD3_CD8_EM_CD28_Br) + (S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_EM_CD28_Br))*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb/(EC50_CD28*(E + S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb +
       S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28) + S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_CD8_N_CD28 <- kact_N*((S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_N_CD28_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3 + S_CD4_N_CD3_CD8_N_CD28_Br) +
       (S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_CD8_N_CD28_Br))*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28) +
       S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_MM_CD28 <- kact_N*((S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD28_Br)/(EC50_CD3*
      (E + S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3 + S_CD4_N_CD3_MM_CD28_Br) +
       (S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD28_Br))*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28) + S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_MM_CD38 <- kact_N*((S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD38_Br)/(EC50_CD3*
      (E + S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3 + S_CD4_N_CD3_MM_CD38_Br) +
       (S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_MM_CD38_Br))*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28) + S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb))
    ka_CD4_N_CD3_TRGT_CD38 <- kact_N*((S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_TRGT_CD38_Br)/
      (EC50_CD3*(E + S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3 + S_CD4_N_CD3_TRGT_CD38_Br) +
       (S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + S_CD4_N_CD3_TRGT_CD38_Br))*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb + S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28) +
       S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb))
    ka_CD8_N_CD3_CD4_A_CD28 <- kact_N*((S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_A_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3 + S_CD8_N_CD3_CD4_A_CD28_Br) +
       (S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_A_CD28_Br))*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28) +
       S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_CD4_EM_CD28 <- kact_N*((S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_EM_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3 +
       S_CD8_N_CD3_CD4_EM_CD28_Br) + (S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_EM_CD28_Br))*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb/(EC50_CD28*(E + S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb +
       S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28) + S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_CD4_N_CD28 <- kact_N*((S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_N_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3 + S_CD8_N_CD3_CD4_N_CD28_Br) +
       (S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD4_N_CD28_Br))*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28) +
       S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_CD8_A_CD28 <- kact_N*((S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_A_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3 + S_CD8_N_CD3_CD8_A_CD28_Br) +
       (S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_A_CD28_Br))*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28) +
       S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_CD8_EM_CD28 <- kact_N*((S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_EM_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3 +
       S_CD8_N_CD3_CD8_EM_CD28_Br) + (S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_EM_CD28_Br))*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb/(EC50_CD28*(E + S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb +
       S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28) + S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_CD8_N_CD28 <- kact_N*((S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_N_CD28_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 + S_CD8_N_CD3_CD8_N_CD28_Br) +
       (S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_CD8_N_CD28_Br))*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28) +
       S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_MM_CD28 <- kact_N*((S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD28_Br)/(EC50_CD3*
      (E + S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3 + S_CD8_N_CD3_MM_CD28_Br) +
       (S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD28_Br))*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28) + S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_MM_CD38 <- kact_N*((S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD38_Br)/(EC50_CD3*
      (E + S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3 + S_CD8_N_CD3_MM_CD38_Br) +
       (S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_MM_CD38_Br))*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28) + S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb))
    ka_CD8_N_CD3_TRGT_CD38 <- kact_N*((S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_TRGT_CD38_Br)/
      (EC50_CD3*(E + S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3 + S_CD8_N_CD3_TRGT_CD38_Br) +
       (S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb + S_CD8_N_CD3_TRGT_CD38_Br))*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb/
      (EC50_CD28*(E + S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb + S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28) +
       S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb))

    d/dt(tsAb) <- -kon_CD3*R_CD8_N_CD3*tsAb + koff_CD3*R_CD8_N_CD3_tsAb - kon_CD28*R_CD8_N_CD28*tsAb +
       koff_CD28*R_CD8_N_CD28_tsAb - kon_CD3*R_CD8_EM_CD3*tsAb + koff_CD3*R_CD8_EM_CD3_tsAb - kon_CD28*
      R_CD8_EM_CD28*tsAb + koff_CD28*R_CD8_EM_CD28_tsAb - kon_CD3*R_CD8_A_CD3*tsAb + koff_CD3*R_CD8_A_CD3_tsAb -
       kon_CD28*R_CD8_A_CD28*tsAb + koff_CD28*R_CD8_A_CD28_tsAb - kon_CD3*R_CD4_N_CD3*tsAb + koff_CD3*
      R_CD4_N_CD3_tsAb - kon_CD28*R_CD4_N_CD28*tsAb + koff_CD28*R_CD4_N_CD28_tsAb - kon_CD3*R_CD4_EM_CD3*
      tsAb + koff_CD3*R_CD4_EM_CD3_tsAb - kon_CD28*R_CD4_EM_CD28*tsAb + koff_CD28*R_CD4_EM_CD28_tsAb -
       kon_CD3*R_CD4_A_CD3*tsAb + koff_CD3*R_CD4_A_CD3_tsAb - kon_CD28*R_CD4_A_CD28*tsAb + koff_CD28*
      R_CD4_A_CD28_tsAb - kon_CD38*R_MM_CD38*tsAb + koff_CD38*R_MM_CD38_tsAb - kon_CD28*R_MM_CD28*tsAb +
       koff_CD28*R_MM_CD28_tsAb - kon_CD38*R_TRGT_CD38*tsAb + koff_CD38*R_TRGT_CD38_tsAb - kon_CD28*R_TRGT_CD28*
      tsAb + koff_CD28*R_TRGT_CD28_tsAb - kon_CD38*R_sCD38_CD38*tsAb + koff_CD38*R_sCD38_CD38_tsAb - kon_CD3*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb - kon_CD28*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb - kon_CD38*
      S_CD4_A_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb - kon_CD28*
      S_CD4_A_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb - kon_CD3*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb - kon_CD38*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38*
      tsAb + koff_CD38*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28*
      tsAb + koff_CD28*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kon_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb - kon_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb - kon_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28*tsAb +
       koff_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb - kon_CD3*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb - kon_CD28*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb - kon_CD38*S_CD4_A_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*
      S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb -
       kon_CD38*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb -
       kon_CD38*S_CD4_EM_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_EM_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb -
       kon_CD38*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD38*S_CD4_EM_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_EM_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb -
       kon_CD38*S_CD4_EM_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_EM_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb -
       kon_CD38*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb -
       kon_CD38*S_CD4_N_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_N_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb -
       kon_CD38*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb -
       kon_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb - kon_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb - kon_CD3*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3*
      tsAb + koff_CD3*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb - kon_CD28*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28*
      tsAb + koff_CD28*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb - kon_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28*tsAb +
       koff_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb - kon_CD3*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3*tsAb +
       koff_CD3*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb - kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28*
      tsAb + koff_CD28*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb - kon_CD38*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38*
      tsAb + koff_CD38*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28*
      tsAb + koff_CD28*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kon_CD3*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb - kon_CD38*S_CD8_A_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_CD8_A_CD28_MM_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb - kon_CD3*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb - kon_CD38*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38*
      tsAb + koff_CD38*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28*
      tsAb + koff_CD28*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb - kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kon_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3*
      tsAb + koff_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb - kon_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb - kon_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28*tsAb +
       koff_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb - kon_CD3*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb - kon_CD28*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb - kon_CD38*S_CD8_A_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*
      S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb -
       kon_CD38*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb -
       kon_CD38*S_CD8_EM_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_EM_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb -
       kon_CD38*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD38*S_CD8_EM_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_EM_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb -
       kon_CD38*S_CD8_EM_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_EM_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb -
       kon_CD38*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb -
       kon_CD38*S_CD8_N_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_N_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb -
       kon_CD38*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kon_CD3*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb -
       kon_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb -
       kon_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb - kon_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb - kon_CD3*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3*
      tsAb + koff_CD3*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb - kon_CD28*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28*
      tsAb + koff_CD28*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb - kon_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28*tsAb +
       koff_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb - kon_CD3*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3*tsAb +
       koff_CD3*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb - kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28*
      tsAb + koff_CD28*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb - kon_CD38*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38*
      tsAb + koff_CD38*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28*
      tsAb + koff_CD28*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kon_CD38*S_MM_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb - kon_CD28*S_MM_CD28_MM_CD38_R_MM_CD28*tsAb +
       koff_CD28*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb - kon_CD38*S_MM_CD28_TRGT_CD38_R_MM_CD38*tsAb + koff_CD38*
      S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb - kon_CD28*S_MM_CD28_TRGT_CD38_R_MM_CD28*tsAb + koff_CD28*S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb -
       kon_CD38*S_MM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kon_CD28*S_MM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kon_CD3*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kon_CD38*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb -
       kon_CD38*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb -
       kon_CD3*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb -
       kon_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kon_CD38*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kon_CD38*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb -
       kon_CD38*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb -
       kon_CD3*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb -
       kon_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kon_CD38*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38*tsAb + koff_CD38*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb -
       kon_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(CD8_N) <- +ksyn_CD8_N - kdeg_CD8_N*CD8_N - kf_CD8_N_CD3_CD8_N_CD28*CD8_N*CD8_N - kf_CD8_N_CD3_CD8_N_CD28*
      CD8_N*CD8_N - kf_CD8_N_CD3_CD8_EM_CD28*CD8_N*CD8_EM - kf_CD8_N_CD3_CD8_A_CD28*CD8_N*CD8_A - kf_CD8_N_CD3_CD4_N_CD28*
      CD8_N*CD4_N - kf_CD8_N_CD3_CD4_EM_CD28*CD8_N*CD4_EM - kf_CD8_N_CD3_CD4_A_CD28*CD8_N*CD4_A - kf_CD8_N_CD3_MM_CD28*
      CD8_N*MM - kf_CD8_EM_CD3_CD8_N_CD28*CD8_EM*CD8_N - kf_CD8_A_CD3_CD8_N_CD28*CD8_A*CD8_N - kf_CD4_N_CD3_CD8_N_CD28*
      CD4_N*CD8_N - kf_CD4_EM_CD3_CD8_N_CD28*CD4_EM*CD8_N - kf_CD4_A_CD3_CD8_N_CD28*CD4_A*CD8_N - kf_CD8_N_CD3_TRGT_CD38*
      CD8_N*TRGT - kf_CD8_N_CD3_MM_CD38*CD8_N*MM - kf_CD8_N_CD28_TRGT_CD38*CD8_N*TRGT - kf_CD8_N_CD28_MM_CD38*
      CD8_N*MM + kDis*S_CD4_A_CD3_CD8_N_CD28 + kDis*S_CD4_EM_CD3_CD8_N_CD28 + kDis*S_CD4_N_CD3_CD8_N_CD28 +
       kDis*S_CD8_A_CD3_CD8_N_CD28 + kDis*S_CD8_EM_CD3_CD8_N_CD28 + kDis*S_CD8_N_CD28_MM_CD38 + kDis*
      S_CD8_N_CD28_TRGT_CD38 + kDis*S_CD8_N_CD3_CD4_A_CD28 + kDis*S_CD8_N_CD3_CD4_EM_CD28 + kDis*S_CD8_N_CD3_CD4_N_CD28 +
       kDis*S_CD8_N_CD3_CD8_A_CD28 + kDis*S_CD8_N_CD3_CD8_EM_CD28 + kDis*S_CD8_N_CD3_CD8_N_CD28 + kDis*
      S_CD8_N_CD3_CD8_N_CD28 + kDis*S_CD8_N_CD3_MM_CD28 + kDis*S_CD8_N_CD3_MM_CD38 + kDis*S_CD8_N_CD3_TRGT_CD38
    d/dt(CD8_EM) <- +ksyn_CD8_EM - kdeg_CD8_EM*CD8_EM - kf_CD8_N_CD3_CD8_EM_CD28*CD8_N*CD8_EM - kf_CD8_EM_CD3_CD8_N_CD28*
      CD8_EM*CD8_N - kf_CD8_EM_CD3_CD8_EM_CD28*CD8_EM*CD8_EM - kf_CD8_EM_CD3_CD8_EM_CD28*CD8_EM*CD8_EM -
       kf_CD8_EM_CD3_CD8_A_CD28*CD8_EM*CD8_A - kf_CD8_EM_CD3_CD4_N_CD28*CD8_EM*CD4_N - kf_CD8_EM_CD3_CD4_EM_CD28*
      CD8_EM*CD4_EM - kf_CD8_EM_CD3_CD4_A_CD28*CD8_EM*CD4_A - kf_CD8_EM_CD3_MM_CD28*CD8_EM*MM - kf_CD8_A_CD3_CD8_EM_CD28*
      CD8_A*CD8_EM - kf_CD4_N_CD3_CD8_EM_CD28*CD4_N*CD8_EM - kf_CD4_EM_CD3_CD8_EM_CD28*CD4_EM*CD8_EM -
       kf_CD4_A_CD3_CD8_EM_CD28*CD4_A*CD8_EM - kf_CD8_EM_CD3_TRGT_CD38*CD8_EM*TRGT - kf_CD8_EM_CD3_MM_CD38*
      CD8_EM*MM - kf_CD8_EM_CD28_TRGT_CD38*CD8_EM*TRGT - kf_CD8_EM_CD28_MM_CD38*CD8_EM*MM + kDis*S_CD4_A_CD3_CD8_EM_CD28 +
       kDis*S_CD4_EM_CD3_CD8_EM_CD28 + kDis*S_CD4_N_CD3_CD8_EM_CD28 + kDis*S_CD8_A_CD3_CD8_EM_CD28 + kDis*
      S_CD8_EM_CD28_MM_CD38 + kDis*S_CD8_EM_CD28_TRGT_CD38 + kDis*S_CD8_EM_CD3_CD4_A_CD28 + kDis*S_CD8_EM_CD3_CD4_EM_CD28 +
       kDis*S_CD8_EM_CD3_CD4_N_CD28 + kDis*S_CD8_EM_CD3_CD8_A_CD28 + kDis*S_CD8_EM_CD3_CD8_EM_CD28 + kDis*
      S_CD8_EM_CD3_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_CD8_N_CD28 + kDis*S_CD8_EM_CD3_MM_CD28 + kDis*S_CD8_EM_CD3_MM_CD38 +
       kDis*S_CD8_EM_CD3_TRGT_CD38 + kDis*S_CD8_N_CD3_CD8_EM_CD28
    d/dt(CD8_A) <- -kdeg_CD8_A*CD8_A + kpr_CD8*CD8_A - kf_CD8_N_CD3_CD8_A_CD28*CD8_N*CD8_A - kf_CD8_EM_CD3_CD8_A_CD28*
      CD8_EM*CD8_A - kf_CD8_A_CD3_CD8_N_CD28*CD8_A*CD8_N - kf_CD8_A_CD3_CD8_EM_CD28*CD8_A*CD8_EM - kf_CD8_A_CD3_CD8_A_CD28*
      CD8_A*CD8_A - kf_CD8_A_CD3_CD8_A_CD28*CD8_A*CD8_A - kf_CD8_A_CD3_CD4_N_CD28*CD8_A*CD4_N - kf_CD8_A_CD3_CD4_EM_CD28*
      CD8_A*CD4_EM - kf_CD8_A_CD3_CD4_A_CD28*CD8_A*CD4_A - kf_CD8_A_CD3_MM_CD28*CD8_A*MM - kf_CD4_N_CD3_CD8_A_CD28*
      CD4_N*CD8_A - kf_CD4_EM_CD3_CD8_A_CD28*CD4_EM*CD8_A - kf_CD4_A_CD3_CD8_A_CD28*CD4_A*CD8_A - kf_CD8_A_CD3_TRGT_CD38*
      CD8_A*TRGT - kf_CD8_A_CD3_MM_CD38*CD8_A*MM - kf_CD8_A_CD28_TRGT_CD38*CD8_A*TRGT - kf_CD8_A_CD28_MM_CD38*
      CD8_A*MM + kkillMM_CD8*S_CD8_A_CD28_MM_CD38 + kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38 + kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28 + kkillMM_CD8*S_CD8_A_CD3_MM_CD38 + kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38 + kDis*
      S_CD4_A_CD3_CD8_A_CD28 + kDis*S_CD4_EM_CD3_CD8_A_CD28 + kDis*S_CD4_N_CD3_CD8_A_CD28 + kDis*S_CD8_A_CD28_TRGT_CD38 +
       kDis*S_CD8_A_CD3_CD4_A_CD28 + kDis*S_CD8_A_CD3_CD4_EM_CD28 + kDis*S_CD8_A_CD3_CD4_N_CD28 + kDis*
      S_CD8_A_CD3_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD8_EM_CD28 + kDis*S_CD8_A_CD3_CD8_N_CD28 +
       kDis*S_CD8_A_CD3_TRGT_CD38 + kDis*S_CD8_EM_CD3_CD8_A_CD28 + kDis*S_CD8_N_CD3_CD8_A_CD28 + kDis*
      S_CD8_A_CD28_MM_CD38_MUT + kDis*S_CD8_A_CD3_MM_CD28_MUT + kDis*S_CD8_A_CD3_MM_CD38_MUT
    d/dt(CD4_N) <- +ksyn_CD4_N - kdeg_CD4_N*CD4_N - kf_CD8_N_CD3_CD4_N_CD28*CD8_N*CD4_N - kf_CD8_EM_CD3_CD4_N_CD28*
      CD8_EM*CD4_N - kf_CD8_A_CD3_CD4_N_CD28*CD8_A*CD4_N - kf_CD4_N_CD3_CD8_N_CD28*CD4_N*CD8_N - kf_CD4_N_CD3_CD8_EM_CD28*
      CD4_N*CD8_EM - kf_CD4_N_CD3_CD8_A_CD28*CD4_N*CD8_A - kf_CD4_N_CD3_CD4_N_CD28*CD4_N*CD4_N - kf_CD4_N_CD3_CD4_N_CD28*
      CD4_N*CD4_N - kf_CD4_N_CD3_CD4_EM_CD28*CD4_N*CD4_EM - kf_CD4_N_CD3_CD4_A_CD28*CD4_N*CD4_A - kf_CD4_N_CD3_MM_CD28*
      CD4_N*MM - kf_CD4_EM_CD3_CD4_N_CD28*CD4_EM*CD4_N - kf_CD4_A_CD3_CD4_N_CD28*CD4_A*CD4_N - kf_CD4_N_CD3_TRGT_CD38*
      CD4_N*TRGT - kf_CD4_N_CD3_MM_CD38*CD4_N*MM - kf_CD4_N_CD28_TRGT_CD38*CD4_N*TRGT - kf_CD4_N_CD28_MM_CD38*
      CD4_N*MM + kDis*S_CD4_A_CD3_CD4_N_CD28 + kDis*S_CD4_EM_CD3_CD4_N_CD28 + kDis*S_CD4_N_CD28_MM_CD38 +
       kDis*S_CD4_N_CD28_TRGT_CD38 + kDis*S_CD4_N_CD3_CD4_A_CD28 + kDis*S_CD4_N_CD3_CD4_EM_CD28 + kDis*
      S_CD4_N_CD3_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD8_A_CD28 + kDis*S_CD4_N_CD3_CD8_EM_CD28 +
       kDis*S_CD4_N_CD3_CD8_N_CD28 + kDis*S_CD4_N_CD3_MM_CD28 + kDis*S_CD4_N_CD3_MM_CD38 + kDis*S_CD4_N_CD3_TRGT_CD38 +
       kDis*S_CD8_A_CD3_CD4_N_CD28 + kDis*S_CD8_EM_CD3_CD4_N_CD28 + kDis*S_CD8_N_CD3_CD4_N_CD28
    d/dt(CD4_EM) <- +ksyn_CD4_EM - kdeg_CD4_EM*CD4_EM - kf_CD8_N_CD3_CD4_EM_CD28*CD8_N*CD4_EM - kf_CD8_EM_CD3_CD4_EM_CD28*
      CD8_EM*CD4_EM - kf_CD8_A_CD3_CD4_EM_CD28*CD8_A*CD4_EM - kf_CD4_N_CD3_CD4_EM_CD28*CD4_N*CD4_EM -
       kf_CD4_EM_CD3_CD8_N_CD28*CD4_EM*CD8_N - kf_CD4_EM_CD3_CD8_EM_CD28*CD4_EM*CD8_EM - kf_CD4_EM_CD3_CD8_A_CD28*
      CD4_EM*CD8_A - kf_CD4_EM_CD3_CD4_N_CD28*CD4_EM*CD4_N - kf_CD4_EM_CD3_CD4_EM_CD28*CD4_EM*CD4_EM -
       kf_CD4_EM_CD3_CD4_EM_CD28*CD4_EM*CD4_EM - kf_CD4_EM_CD3_CD4_A_CD28*CD4_EM*CD4_A - kf_CD4_EM_CD3_MM_CD28*
      CD4_EM*MM - kf_CD4_A_CD3_CD4_EM_CD28*CD4_A*CD4_EM - kf_CD4_EM_CD3_TRGT_CD38*CD4_EM*TRGT - kf_CD4_EM_CD3_MM_CD38*
      CD4_EM*MM - kf_CD4_EM_CD28_TRGT_CD38*CD4_EM*TRGT - kf_CD4_EM_CD28_MM_CD38*CD4_EM*MM + kDis*S_CD4_A_CD3_CD4_EM_CD28 +
       kDis*S_CD4_EM_CD28_MM_CD38 + kDis*S_CD4_EM_CD28_TRGT_CD38 + kDis*S_CD4_EM_CD3_CD4_A_CD28 + kDis*
      S_CD4_EM_CD3_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD4_N_CD28 + kDis*
      S_CD4_EM_CD3_CD8_A_CD28 + kDis*S_CD4_EM_CD3_CD8_EM_CD28 + kDis*S_CD4_EM_CD3_CD8_N_CD28 + kDis*S_CD4_EM_CD3_MM_CD28 +
       kDis*S_CD4_EM_CD3_MM_CD38 + kDis*S_CD4_EM_CD3_TRGT_CD38 + kDis*S_CD4_N_CD3_CD4_EM_CD28 + kDis*
      S_CD8_A_CD3_CD4_EM_CD28 + kDis*S_CD8_EM_CD3_CD4_EM_CD28 + kDis*S_CD8_N_CD3_CD4_EM_CD28
    d/dt(CD4_A) <- -kdeg_CD4_A*CD4_A + kpr_CD4*CD4_A - kf_CD8_N_CD3_CD4_A_CD28*CD8_N*CD4_A - kf_CD8_EM_CD3_CD4_A_CD28*
      CD8_EM*CD4_A - kf_CD8_A_CD3_CD4_A_CD28*CD8_A*CD4_A - kf_CD4_N_CD3_CD4_A_CD28*CD4_N*CD4_A - kf_CD4_EM_CD3_CD4_A_CD28*
      CD4_EM*CD4_A - kf_CD4_A_CD3_CD8_N_CD28*CD4_A*CD8_N - kf_CD4_A_CD3_CD8_EM_CD28*CD4_A*CD8_EM - kf_CD4_A_CD3_CD8_A_CD28*
      CD4_A*CD8_A - kf_CD4_A_CD3_CD4_N_CD28*CD4_A*CD4_N - kf_CD4_A_CD3_CD4_EM_CD28*CD4_A*CD4_EM - kf_CD4_A_CD3_CD4_A_CD28*
      CD4_A*CD4_A - kf_CD4_A_CD3_CD4_A_CD28*CD4_A*CD4_A - kf_CD4_A_CD3_MM_CD28*CD4_A*MM - kf_CD4_A_CD3_TRGT_CD38*
      CD4_A*TRGT - kf_CD4_A_CD3_MM_CD38*CD4_A*MM - kf_CD4_A_CD28_TRGT_CD38*CD4_A*TRGT - kf_CD4_A_CD28_MM_CD38*
      CD4_A*MM + kkillMM_CD4*S_CD4_A_CD28_MM_CD38 + kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38 + kkillMM_CD4*
      S_CD4_A_CD3_MM_CD28 + kkillMM_CD4*S_CD4_A_CD3_MM_CD38 + kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38 + kDis*
      S_CD4_A_CD28_TRGT_CD38 + kDis*S_CD4_A_CD3_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD4_EM_CD28 +
       kDis*S_CD4_A_CD3_CD4_N_CD28 + kDis*S_CD4_A_CD3_CD8_A_CD28 + kDis*S_CD4_A_CD3_CD8_EM_CD28 + kDis*
      S_CD4_A_CD3_CD8_N_CD28 + kDis*S_CD4_A_CD3_TRGT_CD38 + kDis*S_CD4_EM_CD3_CD4_A_CD28 + kDis*S_CD4_N_CD3_CD4_A_CD28 +
       kDis*S_CD8_A_CD3_CD4_A_CD28 + kDis*S_CD8_EM_CD3_CD4_A_CD28 + kDis*S_CD8_N_CD3_CD4_A_CD28 + kDis*
      S_CD4_A_CD28_MM_CD38_MUT + kDis*S_CD4_A_CD3_MM_CD28_MUT + kDis*S_CD4_A_CD3_MM_CD38_MUT
    d/dt(MM) <- -kdeg_MM*MM + kpr_MM*(1 - MM/C_M_MM)*MM - kf_CD8_N_CD3_MM_CD28*CD8_N*MM - kf_CD8_EM_CD3_MM_CD28*
      CD8_EM*MM - kf_CD8_A_CD3_MM_CD28*CD8_A*MM - kf_CD4_N_CD3_MM_CD28*CD4_N*MM - kf_CD4_EM_CD3_MM_CD28*
      CD4_EM*MM - kf_CD4_A_CD3_MM_CD28*CD4_A*MM - kf_CD8_N_CD3_MM_CD38*CD8_N*MM - kf_CD8_EM_CD3_MM_CD38*
      CD8_EM*MM - kf_CD8_A_CD3_MM_CD38*CD8_A*MM - kf_CD4_N_CD3_MM_CD38*CD4_N*MM - kf_CD4_EM_CD3_MM_CD38*
      CD4_EM*MM - kf_CD4_A_CD3_MM_CD38*CD4_A*MM - kf_CD8_N_CD28_MM_CD38*CD8_N*MM - kf_CD8_EM_CD28_MM_CD38*
      CD8_EM*MM - kf_CD8_A_CD28_MM_CD38*CD8_A*MM - kf_CD4_N_CD28_MM_CD38*CD4_N*MM - kf_CD4_EM_CD28_MM_CD38*
      CD4_EM*MM - kf_CD4_A_CD28_MM_CD38*CD4_A*MM - kf_MM_CD28_TRGT_CD38*MM*TRGT - kf_MM_CD28_MM_CD38*
      MM*MM - kf_MM_CD28_MM_CD38*MM*MM + kDis*S_CD4_EM_CD28_MM_CD38 + kDis*S_CD4_EM_CD3_MM_CD28 + kDis*
      S_CD4_EM_CD3_MM_CD38 + kDis*S_CD4_N_CD28_MM_CD38 + kDis*S_CD4_N_CD3_MM_CD28 + kDis*S_CD4_N_CD3_MM_CD38 +
       kDis*S_CD8_EM_CD28_MM_CD38 + kDis*S_CD8_EM_CD3_MM_CD28 + kDis*S_CD8_EM_CD3_MM_CD38 + kDis*S_CD8_N_CD28_MM_CD38 +
       kDis*S_CD8_N_CD3_MM_CD28 + kDis*S_CD8_N_CD3_MM_CD38 + kDis*S_MM_CD28_MM_CD38 + kDis*S_MM_CD28_MM_CD38 +
       kDis*S_MM_CD28_TRGT_CD38 + kDis*S_CD4_A_CD28_MM_CD38_MUT + kDis*S_CD4_A_CD3_MM_CD28_MUT + kDis*
      S_CD4_A_CD3_MM_CD38_MUT + kDis*S_CD8_A_CD28_MM_CD38_MUT + kDis*S_CD8_A_CD3_MM_CD28_MUT + kDis*S_CD8_A_CD3_MM_CD38_MUT
    d/dt(TRGT) <- -kdeg_TRGT*TRGT - kf_CD8_N_CD3_TRGT_CD38*CD8_N*TRGT - kf_CD8_EM_CD3_TRGT_CD38*CD8_EM*
      TRGT - kf_CD8_A_CD3_TRGT_CD38*CD8_A*TRGT - kf_CD4_N_CD3_TRGT_CD38*CD4_N*TRGT - kf_CD4_EM_CD3_TRGT_CD38*
      CD4_EM*TRGT - kf_CD4_A_CD3_TRGT_CD38*CD4_A*TRGT - kf_CD8_N_CD28_TRGT_CD38*CD8_N*TRGT - kf_CD8_EM_CD28_TRGT_CD38*
      CD8_EM*TRGT - kf_CD8_A_CD28_TRGT_CD38*CD8_A*TRGT - kf_CD4_N_CD28_TRGT_CD38*CD4_N*TRGT - kf_CD4_EM_CD28_TRGT_CD38*
      CD4_EM*TRGT - kf_CD4_A_CD28_TRGT_CD38*CD4_A*TRGT - kf_MM_CD28_TRGT_CD38*MM*TRGT + kDis*S_CD4_A_CD28_TRGT_CD38 +
       kDis*S_CD4_A_CD3_TRGT_CD38 + kDis*S_CD4_EM_CD28_TRGT_CD38 + kDis*S_CD4_EM_CD3_TRGT_CD38 + kDis*
      S_CD4_N_CD28_TRGT_CD38 + kDis*S_CD4_N_CD3_TRGT_CD38 + kDis*S_CD8_A_CD28_TRGT_CD38 + kDis*S_CD8_A_CD3_TRGT_CD38 +
       kDis*S_CD8_EM_CD28_TRGT_CD38 + kDis*S_CD8_EM_CD3_TRGT_CD38 + kDis*S_CD8_N_CD28_TRGT_CD38 + kDis*
      S_CD8_N_CD3_TRGT_CD38 + kDis*S_MM_CD28_TRGT_CD38
    d/dt(sCD38) <- +kshedMM_s38*R_MM_CD38 + kshedTRGT_s38*R_TRGT_CD38 - kdeg_sCD38*sCD38
    d/dt(S_CD4_A_CD28_MM_CD38) <- +kf_CD4_A_CD28_MM_CD38*CD4_A*MM - kkillMM_CD4*S_CD4_A_CD28_MM_CD38 -
       kmut_SYN*S_CD4_A_CD28_MM_CD38
    d/dt(S_CD4_A_CD28_TRGT_CD38) <- +kf_CD4_A_CD28_TRGT_CD38*CD4_A*TRGT - kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38 -
       kDis*S_CD4_A_CD28_TRGT_CD38
    d/dt(S_CD4_A_CD3_CD4_A_CD28) <- +kf_CD4_A_CD3_CD4_A_CD28*CD4_A*CD4_A + kact_EM*S_CD4_EM_CD3_CD4_A_CD28 +
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28 - kDis*S_CD4_A_CD3_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD4_EM_CD28) <- +kf_CD4_A_CD3_CD4_EM_CD28*CD4_A*CD4_EM + kact_EM*S_CD4_EM_CD3_CD4_EM_CD28 +
       ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28 - kDis*S_CD4_A_CD3_CD4_EM_CD28
    d/dt(S_CD4_A_CD3_CD4_N_CD28) <- +kf_CD4_A_CD3_CD4_N_CD28*CD4_A*CD4_N + kact_EM*S_CD4_EM_CD3_CD4_N_CD28 +
       ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28 - kDis*S_CD4_A_CD3_CD4_N_CD28
    d/dt(S_CD4_A_CD3_CD8_A_CD28) <- +kf_CD4_A_CD3_CD8_A_CD28*CD4_A*CD8_A + kact_EM*S_CD4_EM_CD3_CD8_A_CD28 +
       ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28 - kDis*S_CD4_A_CD3_CD8_A_CD28
    d/dt(S_CD4_A_CD3_CD8_EM_CD28) <- +kf_CD4_A_CD3_CD8_EM_CD28*CD4_A*CD8_EM + kact_EM*S_CD4_EM_CD3_CD8_EM_CD28 +
       ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28 - kDis*S_CD4_A_CD3_CD8_EM_CD28
    d/dt(S_CD4_A_CD3_CD8_N_CD28) <- +kf_CD4_A_CD3_CD8_N_CD28*CD4_A*CD8_N + kact_EM*S_CD4_EM_CD3_CD8_N_CD28 +
       ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28 - kDis*S_CD4_A_CD3_CD8_N_CD28
    d/dt(S_CD4_A_CD3_MM_CD28) <- +kf_CD4_A_CD3_MM_CD28*CD4_A*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD28 - kmut_SYN*
      S_CD4_A_CD3_MM_CD28 + kact_EM*S_CD4_EM_CD3_MM_CD28 + ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28
    d/dt(S_CD4_A_CD3_MM_CD38) <- +kf_CD4_A_CD3_MM_CD38*CD4_A*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD38 - kmut_SYN*
      S_CD4_A_CD3_MM_CD38 + kact_EM*S_CD4_EM_CD3_MM_CD38 + ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38
    d/dt(S_CD4_A_CD3_TRGT_CD38) <- +kf_CD4_A_CD3_TRGT_CD38*CD4_A*TRGT - kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38 +
       kact_EM*S_CD4_EM_CD3_TRGT_CD38 + ka_CD4_N_CD3_TRGT_CD38*S_CD4_N_CD3_TRGT_CD38 - kDis*S_CD4_A_CD3_TRGT_CD38
    d/dt(S_CD4_EM_CD28_MM_CD38) <- +kf_CD4_EM_CD28_MM_CD38*CD4_EM*MM - kDis*S_CD4_EM_CD28_MM_CD38
    d/dt(S_CD4_EM_CD28_TRGT_CD38) <- +kf_CD4_EM_CD28_TRGT_CD38*CD4_EM*TRGT - kDis*S_CD4_EM_CD28_TRGT_CD38
    d/dt(S_CD4_EM_CD3_CD4_A_CD28) <- +kf_CD4_EM_CD3_CD4_A_CD28*CD4_EM*CD4_A - kact_EM*S_CD4_EM_CD3_CD4_A_CD28 -
       kDis*S_CD4_EM_CD3_CD4_A_CD28
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD4_EM_CD28*CD4_EM*CD4_EM - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28 -
       kDis*S_CD4_EM_CD3_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD4_N_CD28) <- +kf_CD4_EM_CD3_CD4_N_CD28*CD4_EM*CD4_N - kact_EM*S_CD4_EM_CD3_CD4_N_CD28 -
       kDis*S_CD4_EM_CD3_CD4_N_CD28
    d/dt(S_CD4_EM_CD3_CD8_A_CD28) <- +kf_CD4_EM_CD3_CD8_A_CD28*CD4_EM*CD8_A - kact_EM*S_CD4_EM_CD3_CD8_A_CD28 -
       kDis*S_CD4_EM_CD3_CD8_A_CD28
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28) <- +kf_CD4_EM_CD3_CD8_EM_CD28*CD4_EM*CD8_EM - kact_EM*S_CD4_EM_CD3_CD8_EM_CD28 -
       kDis*S_CD4_EM_CD3_CD8_EM_CD28
    d/dt(S_CD4_EM_CD3_CD8_N_CD28) <- +kf_CD4_EM_CD3_CD8_N_CD28*CD4_EM*CD8_N - kact_EM*S_CD4_EM_CD3_CD8_N_CD28 -
       kDis*S_CD4_EM_CD3_CD8_N_CD28
    d/dt(S_CD4_EM_CD3_MM_CD28) <- +kf_CD4_EM_CD3_MM_CD28*CD4_EM*MM - kact_EM*S_CD4_EM_CD3_MM_CD28 - kDis*
      S_CD4_EM_CD3_MM_CD28
    d/dt(S_CD4_EM_CD3_MM_CD38) <- +kf_CD4_EM_CD3_MM_CD38*CD4_EM*MM - kact_EM*S_CD4_EM_CD3_MM_CD38 - kDis*
      S_CD4_EM_CD3_MM_CD38
    d/dt(S_CD4_EM_CD3_TRGT_CD38) <- +kf_CD4_EM_CD3_TRGT_CD38*CD4_EM*TRGT - kact_EM*S_CD4_EM_CD3_TRGT_CD38 -
       kDis*S_CD4_EM_CD3_TRGT_CD38
    d/dt(S_CD4_N_CD28_MM_CD38) <- +kf_CD4_N_CD28_MM_CD38*CD4_N*MM - kDis*S_CD4_N_CD28_MM_CD38
    d/dt(S_CD4_N_CD28_TRGT_CD38) <- +kf_CD4_N_CD28_TRGT_CD38*CD4_N*TRGT - kDis*S_CD4_N_CD28_TRGT_CD38
    d/dt(S_CD4_N_CD3_CD4_A_CD28) <- +kf_CD4_N_CD3_CD4_A_CD28*CD4_N*CD4_A - ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28 -
       kDis*S_CD4_N_CD3_CD4_A_CD28
    d/dt(S_CD4_N_CD3_CD4_EM_CD28) <- +kf_CD4_N_CD3_CD4_EM_CD28*CD4_N*CD4_EM - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28 - kDis*S_CD4_N_CD3_CD4_EM_CD28
    d/dt(S_CD4_N_CD3_CD4_N_CD28) <- +kf_CD4_N_CD3_CD4_N_CD28*CD4_N*CD4_N - ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28 -
       kDis*S_CD4_N_CD3_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD8_A_CD28) <- +kf_CD4_N_CD3_CD8_A_CD28*CD4_N*CD8_A - ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28 -
       kDis*S_CD4_N_CD3_CD8_A_CD28
    d/dt(S_CD4_N_CD3_CD8_EM_CD28) <- +kf_CD4_N_CD3_CD8_EM_CD28*CD4_N*CD8_EM - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28 - kDis*S_CD4_N_CD3_CD8_EM_CD28
    d/dt(S_CD4_N_CD3_CD8_N_CD28) <- +kf_CD4_N_CD3_CD8_N_CD28*CD4_N*CD8_N - ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28 -
       kDis*S_CD4_N_CD3_CD8_N_CD28
    d/dt(S_CD4_N_CD3_MM_CD28) <- +kf_CD4_N_CD3_MM_CD28*CD4_N*MM - ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28 -
       kDis*S_CD4_N_CD3_MM_CD28
    d/dt(S_CD4_N_CD3_MM_CD38) <- +kf_CD4_N_CD3_MM_CD38*CD4_N*MM - ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38 -
       kDis*S_CD4_N_CD3_MM_CD38
    d/dt(S_CD4_N_CD3_TRGT_CD38) <- +kf_CD4_N_CD3_TRGT_CD38*CD4_N*TRGT - ka_CD4_N_CD3_TRGT_CD38*S_CD4_N_CD3_TRGT_CD38 -
       kDis*S_CD4_N_CD3_TRGT_CD38
    d/dt(S_CD8_A_CD28_MM_CD38) <- +kf_CD8_A_CD28_MM_CD38*CD8_A*MM - kkillMM_CD8*S_CD8_A_CD28_MM_CD38 -
       kmut_SYN*S_CD8_A_CD28_MM_CD38
    d/dt(S_CD8_A_CD28_TRGT_CD38) <- +kf_CD8_A_CD28_TRGT_CD38*CD8_A*TRGT - kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38 -
       kDis*S_CD8_A_CD28_TRGT_CD38
    d/dt(S_CD8_A_CD3_CD4_A_CD28) <- +kf_CD8_A_CD3_CD4_A_CD28*CD8_A*CD4_A + kact_EM*S_CD8_EM_CD3_CD4_A_CD28 +
       ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28 - kDis*S_CD8_A_CD3_CD4_A_CD28
    d/dt(S_CD8_A_CD3_CD4_EM_CD28) <- +kf_CD8_A_CD3_CD4_EM_CD28*CD8_A*CD4_EM + kact_EM*S_CD8_EM_CD3_CD4_EM_CD28 +
       ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28 - kDis*S_CD8_A_CD3_CD4_EM_CD28
    d/dt(S_CD8_A_CD3_CD4_N_CD28) <- +kf_CD8_A_CD3_CD4_N_CD28*CD8_A*CD4_N + kact_EM*S_CD8_EM_CD3_CD4_N_CD28 +
       ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28 - kDis*S_CD8_A_CD3_CD4_N_CD28
    d/dt(S_CD8_A_CD3_CD8_A_CD28) <- +kf_CD8_A_CD3_CD8_A_CD28*CD8_A*CD8_A + kact_EM*S_CD8_EM_CD3_CD8_A_CD28 +
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28 - kDis*S_CD8_A_CD3_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD8_EM_CD28) <- +kf_CD8_A_CD3_CD8_EM_CD28*CD8_A*CD8_EM + kact_EM*S_CD8_EM_CD3_CD8_EM_CD28 +
       ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28 - kDis*S_CD8_A_CD3_CD8_EM_CD28
    d/dt(S_CD8_A_CD3_CD8_N_CD28) <- +kf_CD8_A_CD3_CD8_N_CD28*CD8_A*CD8_N + kact_EM*S_CD8_EM_CD3_CD8_N_CD28 +
       ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28 - kDis*S_CD8_A_CD3_CD8_N_CD28
    d/dt(S_CD8_A_CD3_MM_CD28) <- +kf_CD8_A_CD3_MM_CD28*CD8_A*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD28 - kmut_SYN*
      S_CD8_A_CD3_MM_CD28 + kact_EM*S_CD8_EM_CD3_MM_CD28 + ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28
    d/dt(S_CD8_A_CD3_MM_CD38) <- +kf_CD8_A_CD3_MM_CD38*CD8_A*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD38 - kmut_SYN*
      S_CD8_A_CD3_MM_CD38 + kact_EM*S_CD8_EM_CD3_MM_CD38 + ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38
    d/dt(S_CD8_A_CD3_TRGT_CD38) <- +kf_CD8_A_CD3_TRGT_CD38*CD8_A*TRGT - kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38 +
       kact_EM*S_CD8_EM_CD3_TRGT_CD38 + ka_CD8_N_CD3_TRGT_CD38*S_CD8_N_CD3_TRGT_CD38 - kDis*S_CD8_A_CD3_TRGT_CD38
    d/dt(S_CD8_EM_CD28_MM_CD38) <- +kf_CD8_EM_CD28_MM_CD38*CD8_EM*MM - kDis*S_CD8_EM_CD28_MM_CD38
    d/dt(S_CD8_EM_CD28_TRGT_CD38) <- +kf_CD8_EM_CD28_TRGT_CD38*CD8_EM*TRGT - kDis*S_CD8_EM_CD28_TRGT_CD38
    d/dt(S_CD8_EM_CD3_CD4_A_CD28) <- +kf_CD8_EM_CD3_CD4_A_CD28*CD8_EM*CD4_A - kact_EM*S_CD8_EM_CD3_CD4_A_CD28 -
       kDis*S_CD8_EM_CD3_CD4_A_CD28
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28) <- +kf_CD8_EM_CD3_CD4_EM_CD28*CD8_EM*CD4_EM - kact_EM*S_CD8_EM_CD3_CD4_EM_CD28 -
       kDis*S_CD8_EM_CD3_CD4_EM_CD28
    d/dt(S_CD8_EM_CD3_CD4_N_CD28) <- +kf_CD8_EM_CD3_CD4_N_CD28*CD8_EM*CD4_N - kact_EM*S_CD8_EM_CD3_CD4_N_CD28 -
       kDis*S_CD8_EM_CD3_CD4_N_CD28
    d/dt(S_CD8_EM_CD3_CD8_A_CD28) <- +kf_CD8_EM_CD3_CD8_A_CD28*CD8_EM*CD8_A - kact_EM*S_CD8_EM_CD3_CD8_A_CD28 -
       kDis*S_CD8_EM_CD3_CD8_A_CD28
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD8_EM_CD28*CD8_EM*CD8_EM - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28 -
       kDis*S_CD8_EM_CD3_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD8_N_CD28) <- +kf_CD8_EM_CD3_CD8_N_CD28*CD8_EM*CD8_N - kact_EM*S_CD8_EM_CD3_CD8_N_CD28 -
       kDis*S_CD8_EM_CD3_CD8_N_CD28
    d/dt(S_CD8_EM_CD3_MM_CD28) <- +kf_CD8_EM_CD3_MM_CD28*CD8_EM*MM - kact_EM*S_CD8_EM_CD3_MM_CD28 - kDis*
      S_CD8_EM_CD3_MM_CD28
    d/dt(S_CD8_EM_CD3_MM_CD38) <- +kf_CD8_EM_CD3_MM_CD38*CD8_EM*MM - kact_EM*S_CD8_EM_CD3_MM_CD38 - kDis*
      S_CD8_EM_CD3_MM_CD38
    d/dt(S_CD8_EM_CD3_TRGT_CD38) <- +kf_CD8_EM_CD3_TRGT_CD38*CD8_EM*TRGT - kact_EM*S_CD8_EM_CD3_TRGT_CD38 -
       kDis*S_CD8_EM_CD3_TRGT_CD38
    d/dt(S_CD8_N_CD28_MM_CD38) <- +kf_CD8_N_CD28_MM_CD38*CD8_N*MM - kDis*S_CD8_N_CD28_MM_CD38
    d/dt(S_CD8_N_CD28_TRGT_CD38) <- +kf_CD8_N_CD28_TRGT_CD38*CD8_N*TRGT - kDis*S_CD8_N_CD28_TRGT_CD38
    d/dt(S_CD8_N_CD3_CD4_A_CD28) <- +kf_CD8_N_CD3_CD4_A_CD28*CD8_N*CD4_A - ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28 -
       kDis*S_CD8_N_CD3_CD4_A_CD28
    d/dt(S_CD8_N_CD3_CD4_EM_CD28) <- +kf_CD8_N_CD3_CD4_EM_CD28*CD8_N*CD4_EM - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28 - kDis*S_CD8_N_CD3_CD4_EM_CD28
    d/dt(S_CD8_N_CD3_CD4_N_CD28) <- +kf_CD8_N_CD3_CD4_N_CD28*CD8_N*CD4_N - ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28 -
       kDis*S_CD8_N_CD3_CD4_N_CD28
    d/dt(S_CD8_N_CD3_CD8_A_CD28) <- +kf_CD8_N_CD3_CD8_A_CD28*CD8_N*CD8_A - ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28 -
       kDis*S_CD8_N_CD3_CD8_A_CD28
    d/dt(S_CD8_N_CD3_CD8_EM_CD28) <- +kf_CD8_N_CD3_CD8_EM_CD28*CD8_N*CD8_EM - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28 - kDis*S_CD8_N_CD3_CD8_EM_CD28
    d/dt(S_CD8_N_CD3_CD8_N_CD28) <- +kf_CD8_N_CD3_CD8_N_CD28*CD8_N*CD8_N - ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28 -
       kDis*S_CD8_N_CD3_CD8_N_CD28
    d/dt(S_CD8_N_CD3_MM_CD28) <- +kf_CD8_N_CD3_MM_CD28*CD8_N*MM - ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28 -
       kDis*S_CD8_N_CD3_MM_CD28
    d/dt(S_CD8_N_CD3_MM_CD38) <- +kf_CD8_N_CD3_MM_CD38*CD8_N*MM - ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38 -
       kDis*S_CD8_N_CD3_MM_CD38
    d/dt(S_CD8_N_CD3_TRGT_CD38) <- +kf_CD8_N_CD3_TRGT_CD38*CD8_N*TRGT - ka_CD8_N_CD3_TRGT_CD38*S_CD8_N_CD3_TRGT_CD38 -
       kDis*S_CD8_N_CD3_TRGT_CD38
    d/dt(S_MM_CD28_MM_CD38) <- +kf_MM_CD28_MM_CD38*MM*MM - kDis*S_MM_CD28_MM_CD38
    d/dt(S_MM_CD28_TRGT_CD38) <- +kf_MM_CD28_TRGT_CD38*MM*TRGT - kDis*S_MM_CD28_TRGT_CD38
    d/dt(S_CD4_A_CD28_MM_CD38_MUT) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38 - kDis*S_CD4_A_CD28_MM_CD38_MUT
    d/dt(S_CD4_A_CD3_MM_CD28_MUT) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28 - kDis*S_CD4_A_CD3_MM_CD28_MUT
    d/dt(S_CD4_A_CD3_MM_CD38_MUT) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38 - kDis*S_CD4_A_CD3_MM_CD38_MUT
    d/dt(S_CD8_A_CD28_MM_CD38_MUT) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38 - kDis*S_CD8_A_CD28_MM_CD38_MUT
    d/dt(S_CD8_A_CD3_MM_CD28_MUT) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28 - kDis*S_CD8_A_CD3_MM_CD28_MUT
    d/dt(S_CD8_A_CD3_MM_CD38_MUT) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38 - kDis*S_CD8_A_CD3_MM_CD38_MUT
    d/dt(R_CD8_N_CD3) <- +ksyn_CD8_N*CD3per_CD8 - kdeg_CD8_N*R_CD8_N_CD3 - kon_CD3*R_CD8_N_CD3*tsAb +
       koff_CD3*R_CD8_N_CD3_tsAb - kb1_CD8_N_CD3_CD8_N_CD28 - kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_N -
       kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_N - kb1_CD8_N_CD3_CD8_EM_CD28 - kf_CD8_N_CD3_CD8_EM_CD28*
      R_CD8_N_CD3*CD8_EM - kb1_CD8_N_CD3_CD8_A_CD28 - kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD3*CD8_A - kb1_CD8_N_CD3_CD4_N_CD28 -
       kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD3*CD4_N - kb1_CD8_N_CD3_CD4_EM_CD28 - kf_CD8_N_CD3_CD4_EM_CD28*
      R_CD8_N_CD3*CD4_EM - kb1_CD8_N_CD3_CD4_A_CD28 - kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD3*CD4_A - kb1_CD8_N_CD3_MM_CD28 -
       kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD3*MM - kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_EM - kf_CD8_A_CD3_CD8_N_CD28*
      R_CD8_N_CD3*CD8_A - kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD3*CD4_N - kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD3*
      CD4_EM - kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD3*CD4_A - kb1_CD8_N_CD3_TRGT_CD38 - kf_CD8_N_CD3_TRGT_CD38*
      R_CD8_N_CD3*TRGT - kb1_CD8_N_CD3_MM_CD38 - kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD3*MM - kf_CD8_N_CD28_TRGT_CD38*
      R_CD8_N_CD3*TRGT - kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD3*MM + kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3 +
       kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3 + kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3 +
       kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3 + kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3 +
       kDis*S_CD8_N_CD3_CD4_A_CD28_Br + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_CD4_EM_CD28_Br +
       kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_CD4_N_CD28_Br + kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3 +
       kDis*S_CD8_N_CD3_CD8_A_CD28_Br + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_CD8_EM_CD28_Br +
       kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_CD8_N_CD28_Br + kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 +
       kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_MM_CD28_Br + kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3 +
       kDis*S_CD8_N_CD3_MM_CD38_Br + kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3 + kDis*S_CD8_N_CD3_TRGT_CD38_Br +
       kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3
    d/dt(R_CD8_EM_CD3) <- +ksyn_CD8_EM*CD3per_CD8 - kdeg_CD8_EM*R_CD8_EM_CD3 - kon_CD3*R_CD8_EM_CD3*tsAb +
       koff_CD3*R_CD8_EM_CD3_tsAb - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_N - kb1_CD8_EM_CD3_CD8_N_CD28 -
       kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD3*CD8_N - kb1_CD8_EM_CD3_CD8_EM_CD28 - kf_CD8_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD3*CD8_EM - kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_EM - kb1_CD8_EM_CD3_CD8_A_CD28 -
       kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD3*CD8_A - kb1_CD8_EM_CD3_CD4_N_CD28 - kf_CD8_EM_CD3_CD4_N_CD28*
      R_CD8_EM_CD3*CD4_N - kb1_CD8_EM_CD3_CD4_EM_CD28 - kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD3*CD4_EM -
       kb1_CD8_EM_CD3_CD4_A_CD28 - kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD3*CD4_A - kb1_CD8_EM_CD3_MM_CD28 -
       kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD3*MM - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_A - kf_CD4_N_CD3_CD8_EM_CD28*
      R_CD8_EM_CD3*CD4_N - kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD4_EM - kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3*
      CD4_A - kb1_CD8_EM_CD3_TRGT_CD38 - kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD3*TRGT - kb1_CD8_EM_CD3_MM_CD38 -
       kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD3*MM - kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD3*TRGT - kf_CD8_EM_CD28_MM_CD38*
      R_CD8_EM_CD3*MM + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 +
       kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3 +
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_CD4_A_CD28_Br + kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3 +
       kDis*S_CD8_EM_CD3_CD4_EM_CD28_Br + kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_CD4_N_CD28_Br +
       kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_CD8_A_CD28_Br + kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3 +
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_Br + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 +
       kDis*S_CD8_EM_CD3_CD8_N_CD28_Br + kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_MM_CD28_Br +
       kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3 + kDis*S_CD8_EM_CD3_MM_CD38_Br + kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3 +
       kDis*S_CD8_EM_CD3_TRGT_CD38_Br + kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3 + kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(R_CD8_A_CD3) <- -kdeg_CD8_A*R_CD8_A_CD3 + kpr_CD8*CD8_A*CD3per_CD8 - kon_CD3*R_CD8_A_CD3*tsAb +
       koff_CD3*R_CD8_A_CD3_tsAb - kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD3*CD8_N - kf_CD8_EM_CD3_CD8_A_CD28*
      R_CD8_A_CD3*CD8_EM - kb1_CD8_A_CD3_CD8_N_CD28 - kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD3*CD8_N - kb1_CD8_A_CD3_CD8_EM_CD28 -
       kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD3*CD8_EM - kb1_CD8_A_CD3_CD8_A_CD28 - kf_CD8_A_CD3_CD8_A_CD28*
      R_CD8_A_CD3*CD8_A - kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD3*CD8_A - kb1_CD8_A_CD3_CD4_N_CD28 - kf_CD8_A_CD3_CD4_N_CD28*
      R_CD8_A_CD3*CD4_N - kb1_CD8_A_CD3_CD4_EM_CD28 - kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD3*CD4_EM - kb1_CD8_A_CD3_CD4_A_CD28 -
       kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD3*CD4_A - kb1_CD8_A_CD3_MM_CD28 - kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD3*
      MM - kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD3*CD4_N - kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD3*CD4_EM - kf_CD4_A_CD3_CD8_A_CD28*
      R_CD8_A_CD3*CD4_A - kb1_CD8_A_CD3_TRGT_CD38 - kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD3*TRGT - kb1_CD8_A_CD3_MM_CD38 -
       kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD3*MM - kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD3*TRGT - kf_CD8_A_CD28_MM_CD38*
      R_CD8_A_CD3*MM + kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3 + kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3 +
       kkillMM_CD8*S_CD8_A_CD3_MM_CD28_Br + kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3 + kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_Br + kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3 + kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_Br +
       kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3 + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3 +
       kDis*S_CD8_A_CD3_CD4_A_CD28_Br + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_CD4_EM_CD28_Br +
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_CD4_N_CD28_Br + kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3 +
       kDis*S_CD8_A_CD3_CD8_A_CD28_Br + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3 +
       kDis*S_CD8_A_CD3_CD8_EM_CD28_Br + kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_CD8_N_CD28_Br +
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_TRGT_CD38_Br + kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3 +
       kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3 + kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3 +
       kDis*S_CD8_A_CD3_MM_CD28_MUT_Br + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3 + kDis*S_CD8_A_CD3_MM_CD38_MUT_Br +
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3
    d/dt(R_CD4_N_CD3) <- +ksyn_CD4_N*CD3per_CD4 - kdeg_CD4_N*R_CD4_N_CD3 - kon_CD3*R_CD4_N_CD3*tsAb +
       koff_CD3*R_CD4_N_CD3_tsAb - kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD3*CD8_N - kf_CD8_EM_CD3_CD4_N_CD28*
      R_CD4_N_CD3*CD8_EM - kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD3*CD8_A - kb1_CD4_N_CD3_CD8_N_CD28 - kf_CD4_N_CD3_CD8_N_CD28*
      R_CD4_N_CD3*CD8_N - kb1_CD4_N_CD3_CD8_EM_CD28 - kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD3*CD8_EM - kb1_CD4_N_CD3_CD8_A_CD28 -
       kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD3*CD8_A - kb1_CD4_N_CD3_CD4_N_CD28 - kf_CD4_N_CD3_CD4_N_CD28*
      R_CD4_N_CD3*CD4_N - kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3*CD4_N - kb1_CD4_N_CD3_CD4_EM_CD28 - kf_CD4_N_CD3_CD4_EM_CD28*
      R_CD4_N_CD3*CD4_EM - kb1_CD4_N_CD3_CD4_A_CD28 - kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD3*CD4_A - kb1_CD4_N_CD3_MM_CD28 -
       kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD3*MM - kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD3*CD4_EM - kf_CD4_A_CD3_CD4_N_CD28*
      R_CD4_N_CD3*CD4_A - kb1_CD4_N_CD3_TRGT_CD38 - kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD3*TRGT - kb1_CD4_N_CD3_MM_CD38 -
       kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD3*MM - kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD3*TRGT - kf_CD4_N_CD28_MM_CD38*
      R_CD4_N_CD3*MM + kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3 + kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 +
       kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3 + kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_CD4_A_CD28_Br +
       kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_CD4_EM_CD28_Br + kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3 +
       kDis*S_CD4_N_CD3_CD4_N_CD28_Br + kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 +
       kDis*S_CD4_N_CD3_CD8_A_CD28_Br + kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_CD8_EM_CD28_Br +
       kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_CD8_N_CD28_Br + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3 +
       kDis*S_CD4_N_CD3_MM_CD28_Br + kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_MM_CD38_Br +
       kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3 + kDis*S_CD4_N_CD3_TRGT_CD38_Br + kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3 +
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3 + kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 + kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(R_CD4_EM_CD3) <- +ksyn_CD4_EM*CD3per_CD4 - kdeg_CD4_EM*R_CD4_EM_CD3 - kon_CD3*R_CD4_EM_CD3*tsAb +
       koff_CD3*R_CD4_EM_CD3_tsAb - kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD8_N - kf_CD8_EM_CD3_CD4_EM_CD28*
      R_CD4_EM_CD3*CD8_EM - kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD8_A - kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3*
      CD4_N - kb1_CD4_EM_CD3_CD8_N_CD28 - kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD3*CD8_N - kb1_CD4_EM_CD3_CD8_EM_CD28 -
       kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD3*CD8_EM - kb1_CD4_EM_CD3_CD8_A_CD28 - kf_CD4_EM_CD3_CD8_A_CD28*
      R_CD4_EM_CD3*CD8_A - kb1_CD4_EM_CD3_CD4_N_CD28 - kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_EM_CD3*CD4_N - kb1_CD4_EM_CD3_CD4_EM_CD28 -
       kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_EM - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_EM -
       kb1_CD4_EM_CD3_CD4_A_CD28 - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD3*CD4_A - kb1_CD4_EM_CD3_MM_CD28 -
       kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD3*MM - kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_A - kb1_CD4_EM_CD3_TRGT_CD38 -
       kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD3*TRGT - kb1_CD4_EM_CD3_MM_CD38 - kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD3*
      MM - kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD3*TRGT - kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD3*MM + kDis*
      S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3 +
       kDis*S_CD4_EM_CD3_CD4_A_CD28_Br + kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_CD4_EM_CD28_Br +
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + kDis*
      S_CD4_EM_CD3_CD4_N_CD28_Br + kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_CD8_A_CD28_Br +
       kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_CD8_EM_CD28_Br + kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3 +
       kDis*S_CD4_EM_CD3_CD8_N_CD28_Br + kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_MM_CD28_Br +
       kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3 + kDis*S_CD4_EM_CD3_MM_CD38_Br + kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3 +
       kDis*S_CD4_EM_CD3_TRGT_CD38_Br + kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3 + kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 +
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + kDis*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(R_CD4_A_CD3) <- -kdeg_CD4_A*R_CD4_A_CD3 + kpr_CD4*CD4_A*CD3per_CD4 - kon_CD3*R_CD4_A_CD3*tsAb +
       koff_CD3*R_CD4_A_CD3_tsAb - kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD3*CD8_N - kf_CD8_EM_CD3_CD4_A_CD28*
      R_CD4_A_CD3*CD8_EM - kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD3*CD8_A - kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD3*
      CD4_N - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_EM - kb1_CD4_A_CD3_CD8_N_CD28 - kf_CD4_A_CD3_CD8_N_CD28*
      R_CD4_A_CD3*CD8_N - kb1_CD4_A_CD3_CD8_EM_CD28 - kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD3*CD8_EM - kb1_CD4_A_CD3_CD8_A_CD28 -
       kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD3*CD8_A - kb1_CD4_A_CD3_CD4_N_CD28 - kf_CD4_A_CD3_CD4_N_CD28*
      R_CD4_A_CD3*CD4_N - kb1_CD4_A_CD3_CD4_EM_CD28 - kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD3*CD4_EM - kb1_CD4_A_CD3_CD4_A_CD28 -
       kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_A - kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_A - kb1_CD4_A_CD3_MM_CD28 -
       kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD3*MM - kb1_CD4_A_CD3_TRGT_CD38 - kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD3*
      TRGT - kb1_CD4_A_CD3_MM_CD38 - kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD3*MM - kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD3*
      TRGT - kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD3*MM + kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3 + kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3 + kkillMM_CD4*S_CD4_A_CD3_MM_CD28_Br + kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3 +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD38_Br + kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3 + kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_Br + kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3 + kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3 +
       kDis*S_CD4_A_CD3_CD4_A_CD28_Br + kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3 +
       kDis*S_CD4_A_CD3_CD4_EM_CD28_Br + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_CD4_N_CD28_Br +
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_CD8_A_CD28_Br + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3 +
       kDis*S_CD4_A_CD3_CD8_EM_CD28_Br + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_CD8_N_CD28_Br +
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_TRGT_CD38_Br + kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3 +
       kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 + kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3 + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3 +
       kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3 + kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3 +
       kDis*S_CD4_A_CD3_MM_CD28_MUT_Br + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3 + kDis*S_CD4_A_CD3_MM_CD38_MUT_Br +
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3
    d/dt(R_MM_CD38) <- -kdeg_MM*R_MM_CD38 + kpr_MM*(1 - MM/C_M_MM)*MM*CD38per_MM - kshedMM_s38*R_MM_CD38 -
       kon_CD38*R_MM_CD38*tsAb + koff_CD38*R_MM_CD38_tsAb - kf_CD8_N_CD3_MM_CD28*R_MM_CD38*CD8_N - kf_CD8_EM_CD3_MM_CD28*
      R_MM_CD38*CD8_EM - kf_CD8_A_CD3_MM_CD28*R_MM_CD38*CD8_A - kf_CD4_N_CD3_MM_CD28*R_MM_CD38*CD4_N -
       kf_CD4_EM_CD3_MM_CD28*R_MM_CD38*CD4_EM - kf_CD4_A_CD3_MM_CD28*R_MM_CD38*CD4_A - kb2_CD8_N_CD3_MM_CD38 -
       kf_CD8_N_CD3_MM_CD38*R_MM_CD38*CD8_N - kb2_CD8_EM_CD3_MM_CD38 - kf_CD8_EM_CD3_MM_CD38*R_MM_CD38*
      CD8_EM - kb2_CD8_A_CD3_MM_CD38 - kf_CD8_A_CD3_MM_CD38*R_MM_CD38*CD8_A - kb2_CD4_N_CD3_MM_CD38 -
       kf_CD4_N_CD3_MM_CD38*R_MM_CD38*CD4_N - kb2_CD4_EM_CD3_MM_CD38 - kf_CD4_EM_CD3_MM_CD38*R_MM_CD38*
      CD4_EM - kb2_CD4_A_CD3_MM_CD38 - kf_CD4_A_CD3_MM_CD38*R_MM_CD38*CD4_A - kb2_CD8_N_CD28_MM_CD38 -
       kf_CD8_N_CD28_MM_CD38*R_MM_CD38*CD8_N - kb2_CD8_EM_CD28_MM_CD38 - kf_CD8_EM_CD28_MM_CD38*R_MM_CD38*
      CD8_EM - kb2_CD8_A_CD28_MM_CD38 - kf_CD8_A_CD28_MM_CD38*R_MM_CD38*CD8_A - kb2_CD4_N_CD28_MM_CD38 -
       kf_CD4_N_CD28_MM_CD38*R_MM_CD38*CD4_N - kb2_CD4_EM_CD28_MM_CD38 - kf_CD4_EM_CD28_MM_CD38*R_MM_CD38*
      CD4_EM - kb2_CD4_A_CD28_MM_CD38 - kf_CD4_A_CD28_MM_CD38*R_MM_CD38*CD4_A - kf_MM_CD28_TRGT_CD38*
      R_MM_CD38*TRGT - kb2_MM_CD28_MM_CD38 - kf_MM_CD28_MM_CD38*R_MM_CD38*MM - kf_MM_CD28_MM_CD38*R_MM_CD38*
      MM + kDis*S_CD4_EM_CD28_MM_CD38_Br + kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD38 + kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD38 +
       kDis*S_CD4_EM_CD3_MM_CD38_Br + kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD38 + kDis*S_CD4_N_CD28_MM_CD38_Br +
       kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD38 + kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD38 + kDis*S_CD4_N_CD3_MM_CD38_Br +
       kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD38 + kDis*S_CD8_EM_CD28_MM_CD38_Br + kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD38 +
       kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD38 + kDis*S_CD8_EM_CD3_MM_CD38_Br + kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD38 +
       kDis*S_CD8_N_CD28_MM_CD38_Br + kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD38 + kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD38 +
       kDis*S_CD8_N_CD3_MM_CD38_Br + kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD38 + kDis*S_MM_CD28_MM_CD38_Br +
       kDis*S_MM_CD28_MM_CD38_R_MM_CD38 + kDis*S_MM_CD28_MM_CD38_R_MM_CD38 + kDis*S_MM_CD28_TRGT_CD38_R_MM_CD38 +
       kDis*S_CD4_A_CD28_MM_CD38_MUT_Br + kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38 + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38 +
       kDis*S_CD4_A_CD3_MM_CD38_MUT_Br + kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38 + kDis*S_CD8_A_CD28_MM_CD38_MUT_Br +
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38 + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38 + kDis*S_CD8_A_CD3_MM_CD38_MUT_Br +
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38
    d/dt(R_TRGT_CD38) <- -kdeg_TRGT*R_TRGT_CD38 - kshedTRGT_s38*R_TRGT_CD38 - kon_CD38*R_TRGT_CD38*tsAb +
       koff_CD38*R_TRGT_CD38_tsAb - kb2_CD8_N_CD3_TRGT_CD38 - kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD38*CD8_N -
       kb2_CD8_EM_CD3_TRGT_CD38 - kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD38*CD8_EM - kb2_CD8_A_CD3_TRGT_CD38 -
       kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD38*CD8_A - kb2_CD4_N_CD3_TRGT_CD38 - kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD38*
      CD4_N - kb2_CD4_EM_CD3_TRGT_CD38 - kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD38*CD4_EM - kb2_CD4_A_CD3_TRGT_CD38 -
       kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD38*CD4_A - kb2_CD8_N_CD28_TRGT_CD38 - kf_CD8_N_CD28_TRGT_CD38*
      R_TRGT_CD38*CD8_N - kb2_CD8_EM_CD28_TRGT_CD38 - kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD38*CD8_EM - kb2_CD8_A_CD28_TRGT_CD38 -
       kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD38*CD8_A - kb2_CD4_N_CD28_TRGT_CD38 - kf_CD4_N_CD28_TRGT_CD38*
      R_TRGT_CD38*CD4_N - kb2_CD4_EM_CD28_TRGT_CD38 - kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD38*CD4_EM - kb2_CD4_A_CD28_TRGT_CD38 -
       kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD38*CD4_A - kb2_MM_CD28_TRGT_CD38 - kf_MM_CD28_TRGT_CD38*R_TRGT_CD38*
      MM + kDis*S_CD4_A_CD28_TRGT_CD38_Br + kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD4_A_CD3_TRGT_CD38_Br +
       kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD4_EM_CD28_TRGT_CD38_Br + kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38 +
       kDis*S_CD4_EM_CD3_TRGT_CD38_Br + kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD4_N_CD28_TRGT_CD38_Br +
       kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD4_N_CD3_TRGT_CD38_Br + kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38 +
       kDis*S_CD8_A_CD28_TRGT_CD38_Br + kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD8_A_CD3_TRGT_CD38_Br +
       kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD8_EM_CD28_TRGT_CD38_Br + kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38 +
       kDis*S_CD8_EM_CD3_TRGT_CD38_Br + kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD8_N_CD28_TRGT_CD38_Br +
       kDis*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38 + kDis*S_CD8_N_CD3_TRGT_CD38_Br + kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38 +
       kDis*S_MM_CD28_TRGT_CD38_Br + kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(R_sCD38_CD38) <- +kshedMM_s38*R_MM_CD38 + kshedTRGT_s38*R_TRGT_CD38 - kdeg_sCD38*R_sCD38_CD38 -
       kon_CD38*R_sCD38_CD38*tsAb + koff_CD38*R_sCD38_CD38_tsAb
    d/dt(R_CD8_N_CD3_tsAb) <- -kdeg_CD8_N*R_CD8_N_CD3_tsAb + kon_CD3*R_CD8_N_CD3*tsAb - koff_CD3*R_CD8_N_CD3_tsAb -
       kb2_CD8_N_CD3_CD8_N_CD28 - kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_N - kf_CD8_N_CD3_CD8_N_CD28*
      R_CD8_N_CD3_tsAb*CD8_N - kb2_CD8_N_CD3_CD8_EM_CD28 - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD3_tsAb*
      CD8_EM - kb2_CD8_N_CD3_CD8_A_CD28 - kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD3_tsAb*CD8_A - kb2_CD8_N_CD3_CD4_N_CD28 -
       kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD3_tsAb*CD4_N - kb2_CD8_N_CD3_CD4_EM_CD28 - kf_CD8_N_CD3_CD4_EM_CD28*
      R_CD8_N_CD3_tsAb*CD4_EM - kb2_CD8_N_CD3_CD4_A_CD28 - kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD3_tsAb*CD4_A -
       kb2_CD8_N_CD3_MM_CD28 - kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD3_tsAb*MM - kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*
      CD8_EM - kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_A - kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*
      CD4_N - kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD4_EM - kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*
      CD4_A - kb2_CD8_N_CD3_TRGT_CD38 - kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD3_tsAb*TRGT - kb2_CD8_N_CD3_MM_CD38 -
       kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD3_tsAb*MM - kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD3_tsAb*TRGT - kf_CD8_N_CD28_MM_CD38*
      R_CD8_N_CD3_tsAb*MM + kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb +
       kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb +
       kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb +
       kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb + kDis*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb + kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb +
       kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb
    d/dt(R_CD8_EM_CD3_tsAb) <- -kdeg_CD8_EM*R_CD8_EM_CD3_tsAb + kon_CD3*R_CD8_EM_CD3*tsAb - koff_CD3*
      R_CD8_EM_CD3_tsAb - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_N - kb2_CD8_EM_CD3_CD8_N_CD28 -
       kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD3_tsAb*CD8_N - kb2_CD8_EM_CD3_CD8_EM_CD28 - kf_CD8_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD3_tsAb*CD8_EM - kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_EM - kb2_CD8_EM_CD3_CD8_A_CD28 -
       kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD3_tsAb*CD8_A - kb2_CD8_EM_CD3_CD4_N_CD28 - kf_CD8_EM_CD3_CD4_N_CD28*
      R_CD8_EM_CD3_tsAb*CD4_N - kb2_CD8_EM_CD3_CD4_EM_CD28 - kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD3_tsAb*
      CD4_EM - kb2_CD8_EM_CD3_CD4_A_CD28 - kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD3_tsAb*CD4_A - kb2_CD8_EM_CD3_MM_CD28 -
       kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD3_tsAb*MM - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_A -
       kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD4_N - kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*
      CD4_EM - kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD4_A - kb2_CD8_EM_CD3_TRGT_CD38 - kf_CD8_EM_CD3_TRGT_CD38*
      R_CD8_EM_CD3_tsAb*TRGT - kb2_CD8_EM_CD3_MM_CD38 - kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD3_tsAb*MM - kf_CD8_EM_CD28_TRGT_CD38*
      R_CD8_EM_CD3_tsAb*TRGT - kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD3_tsAb*MM + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb +
       kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb + kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb + kDis*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb + kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(R_CD8_A_CD3_tsAb) <- -kdeg_CD8_A*R_CD8_A_CD3_tsAb + kon_CD3*R_CD8_A_CD3*tsAb - koff_CD3*R_CD8_A_CD3_tsAb -
       kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_N - kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_EM -
       kb2_CD8_A_CD3_CD8_N_CD28 - kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD3_tsAb*CD8_N - kb2_CD8_A_CD3_CD8_EM_CD28 -
       kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD3_tsAb*CD8_EM - kb2_CD8_A_CD3_CD8_A_CD28 - kf_CD8_A_CD3_CD8_A_CD28*
      R_CD8_A_CD3_tsAb*CD8_A - kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_A - kb2_CD8_A_CD3_CD4_N_CD28 -
       kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD3_tsAb*CD4_N - kb2_CD8_A_CD3_CD4_EM_CD28 - kf_CD8_A_CD3_CD4_EM_CD28*
      R_CD8_A_CD3_tsAb*CD4_EM - kb2_CD8_A_CD3_CD4_A_CD28 - kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD3_tsAb*CD4_A -
       kb2_CD8_A_CD3_MM_CD28 - kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD3_tsAb*MM - kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*
      CD4_N - kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD4_EM - kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*
      CD4_A - kb2_CD8_A_CD3_TRGT_CD38 - kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD3_tsAb*TRGT - kb2_CD8_A_CD3_MM_CD38 -
       kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD3_tsAb*MM - kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD3_tsAb*TRGT - kf_CD8_A_CD28_MM_CD38*
      R_CD8_A_CD3_tsAb*MM + kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb + kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb +
       kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb + kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb +
       kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb +
       kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb + kDis*
      S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb + kDis*
      S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb +
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb + kDis*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kDis*
      S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb + kDis*
      S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb
    d/dt(R_CD4_N_CD3_tsAb) <- -kdeg_CD4_N*R_CD4_N_CD3_tsAb + kon_CD3*R_CD4_N_CD3*tsAb - koff_CD3*R_CD4_N_CD3_tsAb -
       kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_N - kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_EM -
       kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_A - kb2_CD4_N_CD3_CD8_N_CD28 - kf_CD4_N_CD3_CD8_N_CD28*
      R_CD4_N_CD3_tsAb*CD8_N - kb2_CD4_N_CD3_CD8_EM_CD28 - kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD3_tsAb*
      CD8_EM - kb2_CD4_N_CD3_CD8_A_CD28 - kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD3_tsAb*CD8_A - kb2_CD4_N_CD3_CD4_N_CD28 -
       kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_N - kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_N -
       kb2_CD4_N_CD3_CD4_EM_CD28 - kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD3_tsAb*CD4_EM - kb2_CD4_N_CD3_CD4_A_CD28 -
       kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD3_tsAb*CD4_A - kb2_CD4_N_CD3_MM_CD28 - kf_CD4_N_CD3_MM_CD28*
      R_CD4_N_CD3_tsAb*MM - kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_EM - kf_CD4_A_CD3_CD4_N_CD28*
      R_CD4_N_CD3_tsAb*CD4_A - kb2_CD4_N_CD3_TRGT_CD38 - kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD3_tsAb*TRGT -
       kb2_CD4_N_CD3_MM_CD38 - kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD3_tsAb*MM - kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD3_tsAb*
      TRGT - kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD3_tsAb*MM + kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb + kDis*
      S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb +
       kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kDis*
      S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb + kDis*
      S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb +
       kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kDis*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(R_CD4_EM_CD3_tsAb) <- -kdeg_CD4_EM*R_CD4_EM_CD3_tsAb + kon_CD3*R_CD4_EM_CD3*tsAb - koff_CD3*
      R_CD4_EM_CD3_tsAb - kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD8_N - kf_CD8_EM_CD3_CD4_EM_CD28*
      R_CD4_EM_CD3_tsAb*CD8_EM - kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD8_A - kf_CD4_N_CD3_CD4_EM_CD28*
      R_CD4_EM_CD3_tsAb*CD4_N - kb2_CD4_EM_CD3_CD8_N_CD28 - kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD3_tsAb*
      CD8_N - kb2_CD4_EM_CD3_CD8_EM_CD28 - kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD3_tsAb*CD8_EM - kb2_CD4_EM_CD3_CD8_A_CD28 -
       kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD3_tsAb*CD8_A - kb2_CD4_EM_CD3_CD4_N_CD28 - kf_CD4_EM_CD3_CD4_N_CD28*
      R_CD4_EM_CD3_tsAb*CD4_N - kb2_CD4_EM_CD3_CD4_EM_CD28 - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*
      CD4_EM - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD4_EM - kb2_CD4_EM_CD3_CD4_A_CD28 - kf_CD4_EM_CD3_CD4_A_CD28*
      R_CD4_EM_CD3_tsAb*CD4_A - kb2_CD4_EM_CD3_MM_CD28 - kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD3_tsAb*MM -
       kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD4_A - kb2_CD4_EM_CD3_TRGT_CD38 - kf_CD4_EM_CD3_TRGT_CD38*
      R_CD4_EM_CD3_tsAb*TRGT - kb2_CD4_EM_CD3_MM_CD38 - kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD3_tsAb*MM - kf_CD4_EM_CD28_TRGT_CD38*
      R_CD4_EM_CD3_tsAb*TRGT - kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD3_tsAb*MM + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb +
       kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb + kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb + kDis*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kDis*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(R_CD4_A_CD3_tsAb) <- -kdeg_CD4_A*R_CD4_A_CD3_tsAb + kon_CD3*R_CD4_A_CD3*tsAb - koff_CD3*R_CD4_A_CD3_tsAb -
       kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_N - kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_EM -
       kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_A - kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_N -
       kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_EM - kb2_CD4_A_CD3_CD8_N_CD28 - kf_CD4_A_CD3_CD8_N_CD28*
      R_CD4_A_CD3_tsAb*CD8_N - kb2_CD4_A_CD3_CD8_EM_CD28 - kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD3_tsAb*
      CD8_EM - kb2_CD4_A_CD3_CD8_A_CD28 - kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD3_tsAb*CD8_A - kb2_CD4_A_CD3_CD4_N_CD28 -
       kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD3_tsAb*CD4_N - kb2_CD4_A_CD3_CD4_EM_CD28 - kf_CD4_A_CD3_CD4_EM_CD28*
      R_CD4_A_CD3_tsAb*CD4_EM - kb2_CD4_A_CD3_CD4_A_CD28 - kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_A -
       kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_A - kb2_CD4_A_CD3_MM_CD28 - kf_CD4_A_CD3_MM_CD28*
      R_CD4_A_CD3_tsAb*MM - kb2_CD4_A_CD3_TRGT_CD38 - kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD3_tsAb*TRGT - kb2_CD4_A_CD3_MM_CD38 -
       kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD3_tsAb*MM - kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD3_tsAb*TRGT - kf_CD4_A_CD28_MM_CD38*
      R_CD4_A_CD3_tsAb*MM + kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb + kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb + kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb +
       kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb +
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*
      S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb + kDis*
      S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb + kDis*
      S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb + kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb +
       kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kDis*
      S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb + kDis*
      S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb
    d/dt(R_MM_CD38_tsAb) <- -kdeg_MM*R_MM_CD38_tsAb + kon_CD38*R_MM_CD38*tsAb - koff_CD38*R_MM_CD38_tsAb -
       kf_CD8_N_CD3_MM_CD28*R_MM_CD38_tsAb*CD8_N - kf_CD8_EM_CD3_MM_CD28*R_MM_CD38_tsAb*CD8_EM - kf_CD8_A_CD3_MM_CD28*
      R_MM_CD38_tsAb*CD8_A - kf_CD4_N_CD3_MM_CD28*R_MM_CD38_tsAb*CD4_N - kf_CD4_EM_CD3_MM_CD28*R_MM_CD38_tsAb*
      CD4_EM - kf_CD4_A_CD3_MM_CD28*R_MM_CD38_tsAb*CD4_A - kb1_CD8_N_CD3_MM_CD38 - kf_CD8_N_CD3_MM_CD38*
      R_MM_CD38_tsAb*CD8_N - kb1_CD8_EM_CD3_MM_CD38 - kf_CD8_EM_CD3_MM_CD38*R_MM_CD38_tsAb*CD8_EM - kb1_CD8_A_CD3_MM_CD38 -
       kf_CD8_A_CD3_MM_CD38*R_MM_CD38_tsAb*CD8_A - kb1_CD4_N_CD3_MM_CD38 - kf_CD4_N_CD3_MM_CD38*R_MM_CD38_tsAb*
      CD4_N - kb1_CD4_EM_CD3_MM_CD38 - kf_CD4_EM_CD3_MM_CD38*R_MM_CD38_tsAb*CD4_EM - kb1_CD4_A_CD3_MM_CD38 -
       kf_CD4_A_CD3_MM_CD38*R_MM_CD38_tsAb*CD4_A - kb1_CD8_N_CD28_MM_CD38 - kf_CD8_N_CD28_MM_CD38*R_MM_CD38_tsAb*
      CD8_N - kb1_CD8_EM_CD28_MM_CD38 - kf_CD8_EM_CD28_MM_CD38*R_MM_CD38_tsAb*CD8_EM - kb1_CD8_A_CD28_MM_CD38 -
       kf_CD8_A_CD28_MM_CD38*R_MM_CD38_tsAb*CD8_A - kb1_CD4_N_CD28_MM_CD38 - kf_CD4_N_CD28_MM_CD38*R_MM_CD38_tsAb*
      CD4_N - kb1_CD4_EM_CD28_MM_CD38 - kf_CD4_EM_CD28_MM_CD38*R_MM_CD38_tsAb*CD4_EM - kb1_CD4_A_CD28_MM_CD38 -
       kf_CD4_A_CD28_MM_CD38*R_MM_CD38_tsAb*CD4_A - kf_MM_CD28_TRGT_CD38*R_MM_CD38_tsAb*TRGT - kb1_MM_CD28_MM_CD38 -
       kf_MM_CD28_MM_CD38*R_MM_CD38_tsAb*MM - kf_MM_CD28_MM_CD38*R_MM_CD38_tsAb*MM + kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb +
       kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb + kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb + kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb +
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb + kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb + kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb +
       kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb + kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb + kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb +
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb + kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb + kDis*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb +
       kDis*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb + kDis*S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb + kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb +
       kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb + kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb + kDis*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb + kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb
    d/dt(R_TRGT_CD38_tsAb) <- -kdeg_TRGT*R_TRGT_CD38_tsAb + kon_CD38*R_TRGT_CD38*tsAb - koff_CD38*R_TRGT_CD38_tsAb -
       kb1_CD8_N_CD3_TRGT_CD38 - kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_N - kb1_CD8_EM_CD3_TRGT_CD38 -
       kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_EM - kb1_CD8_A_CD3_TRGT_CD38 - kf_CD8_A_CD3_TRGT_CD38*
      R_TRGT_CD38_tsAb*CD8_A - kb1_CD4_N_CD3_TRGT_CD38 - kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_N -
       kb1_CD4_EM_CD3_TRGT_CD38 - kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_EM - kb1_CD4_A_CD3_TRGT_CD38 -
       kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_A - kb1_CD8_N_CD28_TRGT_CD38 - kf_CD8_N_CD28_TRGT_CD38*
      R_TRGT_CD38_tsAb*CD8_N - kb1_CD8_EM_CD28_TRGT_CD38 - kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*
      CD8_EM - kb1_CD8_A_CD28_TRGT_CD38 - kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_A - kb1_CD4_N_CD28_TRGT_CD38 -
       kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_N - kb1_CD4_EM_CD28_TRGT_CD38 - kf_CD4_EM_CD28_TRGT_CD38*
      R_TRGT_CD38_tsAb*CD4_EM - kb1_CD4_A_CD28_TRGT_CD38 - kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_A -
       kb1_MM_CD28_TRGT_CD38 - kf_MM_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*MM + kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb +
       kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*
      S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb +
       kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*
      S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*
      S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(R_sCD38_CD38_tsAb) <- -kdeg_sCD38*R_sCD38_CD38_tsAb + kon_CD38*R_sCD38_CD38*tsAb - koff_CD38*
      R_sCD38_CD38_tsAb
    d/dt(R_CD8_N_CD28) <- +ksyn_CD8_N*CD28per_CD8 - kdeg_CD8_N*R_CD8_N_CD28 - kon_CD28*R_CD8_N_CD28*tsAb +
       koff_CD28*R_CD8_N_CD28_tsAb - kb2_CD8_N_CD3_CD8_N_CD28 - kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28*
      CD8_N - kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_N - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD28*CD8_EM -
       kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD28*CD8_A - kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD28*CD4_N - kf_CD8_N_CD3_CD4_EM_CD28*
      R_CD8_N_CD28*CD4_EM - kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD28*CD4_A - kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD28*
      MM - kb2_CD8_EM_CD3_CD8_N_CD28 - kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_EM - kb2_CD8_A_CD3_CD8_N_CD28 -
       kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_A - kb2_CD4_N_CD3_CD8_N_CD28 - kf_CD4_N_CD3_CD8_N_CD28*
      R_CD8_N_CD28*CD4_N - kb2_CD4_EM_CD3_CD8_N_CD28 - kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD28*CD4_EM -
       kb2_CD4_A_CD3_CD8_N_CD28 - kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD28*CD4_A - kf_CD8_N_CD3_TRGT_CD38*
      R_CD8_N_CD28*TRGT - kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD28*MM - kb1_CD8_N_CD28_TRGT_CD38 - kf_CD8_N_CD28_TRGT_CD38*
      R_CD8_N_CD28*TRGT - kb1_CD8_N_CD28_MM_CD38 - kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD28*MM + kDis*S_CD4_A_CD3_CD8_N_CD28_Br +
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28 + kDis*S_CD4_EM_CD3_CD8_N_CD28_Br + kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 +
       kDis*S_CD4_N_CD3_CD8_N_CD28_Br + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28 + kDis*S_CD8_A_CD3_CD8_N_CD28_Br +
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28 + kDis*S_CD8_EM_CD3_CD8_N_CD28_Br + kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 +
       kDis*S_CD8_N_CD28_MM_CD38_Br + kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28 + kDis*S_CD8_N_CD28_TRGT_CD38_Br +
       kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28 +
       kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28 +
       kDis*S_CD8_N_CD3_CD8_N_CD28_Br + kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 +
       kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28 + kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28
    d/dt(R_CD8_EM_CD28) <- +ksyn_CD8_EM*CD28per_CD8 - kdeg_CD8_EM*R_CD8_EM_CD28 - kon_CD28*R_CD8_EM_CD28*
      tsAb + koff_CD28*R_CD8_EM_CD28_tsAb - kb2_CD8_N_CD3_CD8_EM_CD28 - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28*
      CD8_N - kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD28*CD8_N - kb2_CD8_EM_CD3_CD8_EM_CD28 - kf_CD8_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD28*CD8_EM - kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_EM - kf_CD8_EM_CD3_CD8_A_CD28*
      R_CD8_EM_CD28*CD8_A - kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD28*CD4_N - kf_CD8_EM_CD3_CD4_EM_CD28*
      R_CD8_EM_CD28*CD4_EM - kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD28*CD4_A - kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD28*
      MM - kb2_CD8_A_CD3_CD8_EM_CD28 - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_A - kb2_CD4_N_CD3_CD8_EM_CD28 -
       kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD4_N - kb2_CD4_EM_CD3_CD8_EM_CD28 - kf_CD4_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD28*CD4_EM - kb2_CD4_A_CD3_CD8_EM_CD28 - kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD4_A -
       kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD28*TRGT - kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD28*MM - kb1_CD8_EM_CD28_TRGT_CD38 -
       kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD28*TRGT - kb1_CD8_EM_CD28_MM_CD38 - kf_CD8_EM_CD28_MM_CD38*
      R_CD8_EM_CD28*MM + kDis*S_CD4_A_CD3_CD8_EM_CD28_Br + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28 +
       kDis*S_CD4_EM_CD3_CD8_EM_CD28_Br + kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + kDis*S_CD4_N_CD3_CD8_EM_CD28_Br +
       kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + kDis*S_CD8_A_CD3_CD8_EM_CD28_Br + kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28 +
       kDis*S_CD8_EM_CD28_MM_CD38_Br + kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD28_TRGT_CD38_Br +
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28 + kDis*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28 +
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_Br + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 +
       kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28 + kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28 +
       kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28 + kDis*S_CD8_N_CD3_CD8_EM_CD28_Br + kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(R_CD8_A_CD28) <- -kdeg_CD8_A*R_CD8_A_CD28 + kpr_CD8*CD8_A*CD28per_CD8 - kon_CD28*R_CD8_A_CD28*
      tsAb + koff_CD28*R_CD8_A_CD28_tsAb - kb2_CD8_N_CD3_CD8_A_CD28 - kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD28*
      CD8_N - kb2_CD8_EM_CD3_CD8_A_CD28 - kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD28*CD8_EM - kf_CD8_A_CD3_CD8_N_CD28*
      R_CD8_A_CD28*CD8_N - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD28*CD8_EM - kb2_CD8_A_CD3_CD8_A_CD28 - kf_CD8_A_CD3_CD8_A_CD28*
      R_CD8_A_CD28*CD8_A - kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28*CD8_A - kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD28*
      CD4_N - kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD28*CD4_EM - kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD28*CD4_A -
       kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD28*MM - kb2_CD4_N_CD3_CD8_A_CD28 - kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD28*
      CD4_N - kb2_CD4_EM_CD3_CD8_A_CD28 - kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD28*CD4_EM - kb2_CD4_A_CD3_CD8_A_CD28 -
       kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD28*CD4_A - kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD28*TRGT - kf_CD8_A_CD3_MM_CD38*
      R_CD8_A_CD28*MM - kb1_CD8_A_CD28_TRGT_CD38 - kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD28*TRGT - kb1_CD8_A_CD28_MM_CD38 -
       kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD28*MM + kkillMM_CD8*S_CD8_A_CD28_MM_CD38_Br + kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28 +
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_Br + kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28 + kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28 + kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28 + kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28 + kDis*S_CD4_A_CD3_CD8_A_CD28_Br + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28 +
       kDis*S_CD4_EM_CD3_CD8_A_CD28_Br + kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 + kDis*S_CD4_N_CD3_CD8_A_CD28_Br +
       kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD28_TRGT_CD38_Br + kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28 +
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28 +
       kDis*S_CD8_A_CD3_CD8_A_CD28_Br + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28 +
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28 +
       kDis*S_CD8_EM_CD3_CD8_A_CD28_Br + kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 + kDis*S_CD8_N_CD3_CD8_A_CD28_Br +
       kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28 + kDis*S_CD8_A_CD28_MM_CD38_MUT_Br + kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28 +
       kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28 + kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28
    d/dt(R_CD4_N_CD28) <- +ksyn_CD4_N*CD28per_CD4 - kdeg_CD4_N*R_CD4_N_CD28 - kon_CD28*R_CD4_N_CD28*tsAb +
       koff_CD28*R_CD4_N_CD28_tsAb - kb2_CD8_N_CD3_CD4_N_CD28 - kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD28*
      CD8_N - kb2_CD8_EM_CD3_CD4_N_CD28 - kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD28*CD8_EM - kb2_CD8_A_CD3_CD4_N_CD28 -
       kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD28*CD8_A - kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD28*CD8_N - kf_CD4_N_CD3_CD8_EM_CD28*
      R_CD4_N_CD28*CD8_EM - kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD28*CD8_A - kb2_CD4_N_CD3_CD4_N_CD28 - kf_CD4_N_CD3_CD4_N_CD28*
      R_CD4_N_CD28*CD4_N - kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD28*CD4_N - kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD28*
      CD4_EM - kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD28*CD4_A - kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD28*MM - kb2_CD4_EM_CD3_CD4_N_CD28 -
       kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD28*CD4_EM - kb2_CD4_A_CD3_CD4_N_CD28 - kf_CD4_A_CD3_CD4_N_CD28*
      R_CD4_N_CD28*CD4_A - kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD28*TRGT - kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD28*
      MM - kb1_CD4_N_CD28_TRGT_CD38 - kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD28*TRGT - kb1_CD4_N_CD28_MM_CD38 -
       kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD28*MM + kDis*S_CD4_A_CD3_CD4_N_CD28_Br + kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28 +
       kDis*S_CD4_EM_CD3_CD4_N_CD28_Br + kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD28_MM_CD38_Br +
       kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28 + kDis*S_CD4_N_CD28_TRGT_CD38_Br + kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28 +
       kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD4_N_CD28_Br +
       kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28 +
       kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28 +
       kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28 + kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28 + kDis*S_CD8_A_CD3_CD4_N_CD28_Br +
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28 + kDis*S_CD8_EM_CD3_CD4_N_CD28_Br + kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 +
       kDis*S_CD8_N_CD3_CD4_N_CD28_Br + kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(R_CD4_EM_CD28) <- +ksyn_CD4_EM*CD28per_CD4 - kdeg_CD4_EM*R_CD4_EM_CD28 - kon_CD28*R_CD4_EM_CD28*
      tsAb + koff_CD28*R_CD4_EM_CD28_tsAb - kb2_CD8_N_CD3_CD4_EM_CD28 - kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28*
      CD8_N - kb2_CD8_EM_CD3_CD4_EM_CD28 - kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD8_EM - kb2_CD8_A_CD3_CD4_EM_CD28 -
       kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD8_A - kb2_CD4_N_CD3_CD4_EM_CD28 - kf_CD4_N_CD3_CD4_EM_CD28*
      R_CD4_EM_CD28*CD4_N - kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD28*CD8_N - kf_CD4_EM_CD3_CD8_EM_CD28*
      R_CD4_EM_CD28*CD8_EM - kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD28*CD8_A - kf_CD4_EM_CD3_CD4_N_CD28*
      R_CD4_EM_CD28*CD4_N - kb2_CD4_EM_CD3_CD4_EM_CD28 - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_EM -
       kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_EM - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD28*CD4_A -
       kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD28*MM - kb2_CD4_A_CD3_CD4_EM_CD28 - kf_CD4_A_CD3_CD4_EM_CD28*
      R_CD4_EM_CD28*CD4_A - kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD28*TRGT - kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD28*
      MM - kb1_CD4_EM_CD28_TRGT_CD38 - kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD28*TRGT - kb1_CD4_EM_CD28_MM_CD38 -
       kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD28*MM + kDis*S_CD4_A_CD3_CD4_EM_CD28_Br + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28 +
       kDis*S_CD4_EM_CD28_MM_CD38_Br + kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD28_TRGT_CD38_Br +
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28 + kDis*
      S_CD4_EM_CD3_CD4_EM_CD28_Br + kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 +
       kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28 + kDis*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28 +
       kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28 + kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28 + kDis*S_CD4_N_CD3_CD4_EM_CD28_Br +
       kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + kDis*S_CD8_A_CD3_CD4_EM_CD28_Br + kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28 +
       kDis*S_CD8_EM_CD3_CD4_EM_CD28_Br + kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + kDis*S_CD8_N_CD3_CD4_EM_CD28_Br +
       kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(R_CD4_A_CD28) <- -kdeg_CD4_A*R_CD4_A_CD28 + kpr_CD4*CD4_A*CD28per_CD4 - kon_CD28*R_CD4_A_CD28*
      tsAb + koff_CD28*R_CD4_A_CD28_tsAb - kb2_CD8_N_CD3_CD4_A_CD28 - kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD28*
      CD8_N - kb2_CD8_EM_CD3_CD4_A_CD28 - kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD28*CD8_EM - kb2_CD8_A_CD3_CD4_A_CD28 -
       kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD28*CD8_A - kb2_CD4_N_CD3_CD4_A_CD28 - kf_CD4_N_CD3_CD4_A_CD28*
      R_CD4_A_CD28*CD4_N - kb2_CD4_EM_CD3_CD4_A_CD28 - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD28*CD4_EM -
       kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD28*CD8_N - kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD28*CD8_EM - kf_CD4_A_CD3_CD8_A_CD28*
      R_CD4_A_CD28*CD8_A - kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD28*CD4_N - kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD28*
      CD4_EM - kb2_CD4_A_CD3_CD4_A_CD28 - kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD28*CD4_A - kf_CD4_A_CD3_CD4_A_CD28*
      R_CD4_A_CD28*CD4_A - kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD28*MM - kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD28*
      TRGT - kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD28*MM - kb1_CD4_A_CD28_TRGT_CD38 - kf_CD4_A_CD28_TRGT_CD38*
      R_CD4_A_CD28*TRGT - kb1_CD4_A_CD28_MM_CD38 - kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD28*MM + kkillMM_CD4*
      S_CD4_A_CD28_MM_CD38_Br + kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28 + kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_Br +
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28 + kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28 +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28 + kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28 +
       kDis*S_CD4_A_CD28_TRGT_CD38_Br + kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD4_A_CD28_Br +
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28 +
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28 +
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28 + kDis*S_CD4_EM_CD3_CD4_A_CD28_Br +
       kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 + kDis*S_CD4_N_CD3_CD4_A_CD28_Br + kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28 +
       kDis*S_CD8_A_CD3_CD4_A_CD28_Br + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28 + kDis*S_CD8_EM_CD3_CD4_A_CD28_Br +
       kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 + kDis*S_CD8_N_CD3_CD4_A_CD28_Br + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28 +
       kDis*S_CD4_A_CD28_MM_CD38_MUT_Br + kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28 + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28 +
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28
    d/dt(R_MM_CD28) <- -kdeg_MM*R_MM_CD28 + kpr_MM*(1 - MM/C_M_MM)*MM*CD28per_MM - kon_CD28*R_MM_CD28*
      tsAb + koff_CD28*R_MM_CD28_tsAb - kb2_CD8_N_CD3_MM_CD28 - kf_CD8_N_CD3_MM_CD28*R_MM_CD28*CD8_N -
       kb2_CD8_EM_CD3_MM_CD28 - kf_CD8_EM_CD3_MM_CD28*R_MM_CD28*CD8_EM - kb2_CD8_A_CD3_MM_CD28 - kf_CD8_A_CD3_MM_CD28*
      R_MM_CD28*CD8_A - kb2_CD4_N_CD3_MM_CD28 - kf_CD4_N_CD3_MM_CD28*R_MM_CD28*CD4_N - kb2_CD4_EM_CD3_MM_CD28 -
       kf_CD4_EM_CD3_MM_CD28*R_MM_CD28*CD4_EM - kb2_CD4_A_CD3_MM_CD28 - kf_CD4_A_CD3_MM_CD28*R_MM_CD28*
      CD4_A - kf_CD8_N_CD3_MM_CD38*R_MM_CD28*CD8_N - kf_CD8_EM_CD3_MM_CD38*R_MM_CD28*CD8_EM - kf_CD8_A_CD3_MM_CD38*
      R_MM_CD28*CD8_A - kf_CD4_N_CD3_MM_CD38*R_MM_CD28*CD4_N - kf_CD4_EM_CD3_MM_CD38*R_MM_CD28*CD4_EM -
       kf_CD4_A_CD3_MM_CD38*R_MM_CD28*CD4_A - kf_CD8_N_CD28_MM_CD38*R_MM_CD28*CD8_N - kf_CD8_EM_CD28_MM_CD38*
      R_MM_CD28*CD8_EM - kf_CD8_A_CD28_MM_CD38*R_MM_CD28*CD8_A - kf_CD4_N_CD28_MM_CD38*R_MM_CD28*CD4_N -
       kf_CD4_EM_CD28_MM_CD38*R_MM_CD28*CD4_EM - kf_CD4_A_CD28_MM_CD38*R_MM_CD28*CD4_A - kb1_MM_CD28_TRGT_CD38 -
       kf_MM_CD28_TRGT_CD38*R_MM_CD28*TRGT - kb1_MM_CD28_MM_CD38 - kf_MM_CD28_MM_CD38*R_MM_CD28*MM - kf_MM_CD28_MM_CD38*
      R_MM_CD28*MM + kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD28 + kDis*S_CD4_EM_CD3_MM_CD28_Br + kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD28 +
       kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD28 + kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD28 + kDis*S_CD4_N_CD3_MM_CD28_Br +
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD28 + kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD28 + kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD28 +
       kDis*S_CD8_EM_CD3_MM_CD28_Br + kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD28 + kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD28 +
       kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD28 + kDis*S_CD8_N_CD3_MM_CD28_Br + kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD28 +
       kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD28 + kDis*S_MM_CD28_MM_CD38_Br + kDis*S_MM_CD28_MM_CD38_R_MM_CD28 +
       kDis*S_MM_CD28_MM_CD38_R_MM_CD28 + kDis*S_MM_CD28_TRGT_CD38_Br + kDis*S_MM_CD28_TRGT_CD38_R_MM_CD28 +
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28 + kDis*S_CD4_A_CD3_MM_CD28_MUT_Br + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28 +
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28 + kDis*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28 + kDis*S_CD8_A_CD3_MM_CD28_MUT_Br +
       kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28 + kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28
    d/dt(R_TRGT_CD28) <- -kdeg_TRGT*R_TRGT_CD28 - kon_CD28*R_TRGT_CD28*tsAb + koff_CD28*R_TRGT_CD28_tsAb -
       kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD28*CD8_N - kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD28*CD8_EM - kf_CD8_A_CD3_TRGT_CD38*
      R_TRGT_CD28*CD8_A - kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD28*CD4_N - kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD28*
      CD4_EM - kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD28*CD4_A - kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD28*CD8_N -
       kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD28*CD8_EM - kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD28*CD8_A - kf_CD4_N_CD28_TRGT_CD38*
      R_TRGT_CD28*CD4_N - kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD28*CD4_EM - kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD28*
      CD4_A - kf_MM_CD28_TRGT_CD38*R_TRGT_CD28*MM + kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28 +
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28 +
       kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28 +
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28 + kDis*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28 +
       kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28 + kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(R_CD8_N_CD28_tsAb) <- -kdeg_CD8_N*R_CD8_N_CD28_tsAb + kon_CD28*R_CD8_N_CD28*tsAb - koff_CD28*
      R_CD8_N_CD28_tsAb - kb1_CD8_N_CD3_CD8_N_CD28 - kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_N -
       kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_N - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD28_tsAb*CD8_EM -
       kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD28_tsAb*CD8_A - kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD28_tsAb*CD4_N -
       kf_CD8_N_CD3_CD4_EM_CD28*R_CD8_N_CD28_tsAb*CD4_EM - kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD28_tsAb*
      CD4_A - kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD28_tsAb*MM - kb1_CD8_EM_CD3_CD8_N_CD28 - kf_CD8_EM_CD3_CD8_N_CD28*
      R_CD8_N_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_CD8_N_CD28 - kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*
      CD8_A - kb1_CD4_N_CD3_CD8_N_CD28 - kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_CD8_N_CD28 -
       kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD4_EM - kb1_CD4_A_CD3_CD8_N_CD28 - kf_CD4_A_CD3_CD8_N_CD28*
      R_CD8_N_CD28_tsAb*CD4_A - kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD28_tsAb*TRGT - kf_CD8_N_CD3_MM_CD38*
      R_CD8_N_CD28_tsAb*MM - kb2_CD8_N_CD28_TRGT_CD38 - kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD28_tsAb*TRGT -
       kb2_CD8_N_CD28_MM_CD38 - kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD28_tsAb*MM + kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb + kDis*
      S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb + kDis*
      S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb + kDis*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kDis*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb + kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb +
       kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb
    d/dt(R_CD8_EM_CD28_tsAb) <- -kdeg_CD8_EM*R_CD8_EM_CD28_tsAb + kon_CD28*R_CD8_EM_CD28*tsAb - koff_CD28*
      R_CD8_EM_CD28_tsAb - kb1_CD8_N_CD3_CD8_EM_CD28 - kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*CD8_N -
       kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD28_tsAb*CD8_N - kb1_CD8_EM_CD3_CD8_EM_CD28 - kf_CD8_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD28_tsAb*CD8_EM - kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*CD8_EM - kf_CD8_EM_CD3_CD8_A_CD28*
      R_CD8_EM_CD28_tsAb*CD8_A - kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD28_tsAb*CD4_N - kf_CD8_EM_CD3_CD4_EM_CD28*
      R_CD8_EM_CD28_tsAb*CD4_EM - kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD28_tsAb*CD4_A - kf_CD8_EM_CD3_MM_CD28*
      R_CD8_EM_CD28_tsAb*MM - kb1_CD8_A_CD3_CD8_EM_CD28 - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD8_A - kb1_CD4_N_CD3_CD8_EM_CD28 - kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_CD8_EM_CD28 -
       kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*CD4_EM - kb1_CD4_A_CD3_CD8_EM_CD28 - kf_CD4_A_CD3_CD8_EM_CD28*
      R_CD8_EM_CD28_tsAb*CD4_A - kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD28_tsAb*TRGT - kf_CD8_EM_CD3_MM_CD38*
      R_CD8_EM_CD28_tsAb*MM - kb2_CD8_EM_CD28_TRGT_CD38 - kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD28_tsAb*
      TRGT - kb2_CD8_EM_CD28_MM_CD38 - kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD28_tsAb*MM + kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb + kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb + kDis*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb + kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(R_CD8_A_CD28_tsAb) <- -kdeg_CD8_A*R_CD8_A_CD28_tsAb + kon_CD28*R_CD8_A_CD28*tsAb - koff_CD28*
      R_CD8_A_CD28_tsAb - kb1_CD8_N_CD3_CD8_A_CD28 - kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_N -
       kb1_CD8_EM_CD3_CD8_A_CD28 - kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_EM - kf_CD8_A_CD3_CD8_N_CD28*
      R_CD8_A_CD28_tsAb*CD8_N - kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_CD8_A_CD28 -
       kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_A - kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_A -
       kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD28_tsAb*CD4_N - kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD28_tsAb*CD4_EM -
       kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD28_tsAb*CD4_A - kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD28_tsAb*MM - kb1_CD4_N_CD3_CD8_A_CD28 -
       kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_CD8_A_CD28 - kf_CD4_EM_CD3_CD8_A_CD28*
      R_CD8_A_CD28_tsAb*CD4_EM - kb1_CD4_A_CD3_CD8_A_CD28 - kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*
      CD4_A - kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD28_tsAb*TRGT - kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD28_tsAb*
      MM - kb2_CD8_A_CD28_TRGT_CD38 - kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD28_tsAb*TRGT - kb2_CD8_A_CD28_MM_CD38 -
       kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD28_tsAb*MM + kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb +
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb + kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb +
       kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb + kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb +
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb
    d/dt(R_CD4_N_CD28_tsAb) <- -kdeg_CD4_N*R_CD4_N_CD28_tsAb + kon_CD28*R_CD4_N_CD28*tsAb - koff_CD28*
      R_CD4_N_CD28_tsAb - kb1_CD8_N_CD3_CD4_N_CD28 - kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_N -
       kb1_CD8_EM_CD3_CD4_N_CD28 - kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_CD4_N_CD28 -
       kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_A - kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD28_tsAb*CD8_N -
       kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD28_tsAb*CD8_EM - kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD28_tsAb*
      CD8_A - kb1_CD4_N_CD3_CD4_N_CD28 - kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_N - kf_CD4_N_CD3_CD4_N_CD28*
      R_CD4_N_CD28_tsAb*CD4_N - kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD28_tsAb*CD4_EM - kf_CD4_N_CD3_CD4_A_CD28*
      R_CD4_N_CD28_tsAb*CD4_A - kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD28_tsAb*MM - kb1_CD4_EM_CD3_CD4_N_CD28 -
       kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_EM - kb1_CD4_A_CD3_CD4_N_CD28 - kf_CD4_A_CD3_CD4_N_CD28*
      R_CD4_N_CD28_tsAb*CD4_A - kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD28_tsAb*TRGT - kf_CD4_N_CD3_MM_CD38*
      R_CD4_N_CD28_tsAb*MM - kb2_CD4_N_CD28_TRGT_CD38 - kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD28_tsAb*TRGT -
       kb2_CD4_N_CD28_MM_CD38 - kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD28_tsAb*MM + kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb +
       kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb + kDis*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb + kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kDis*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(R_CD4_EM_CD28_tsAb) <- -kdeg_CD4_EM*R_CD4_EM_CD28_tsAb + kon_CD28*R_CD4_EM_CD28*tsAb - koff_CD28*
      R_CD4_EM_CD28_tsAb - kb1_CD8_N_CD3_CD4_EM_CD28 - kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*CD8_N -
       kb1_CD8_EM_CD3_CD4_EM_CD28 - kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_CD4_EM_CD28 -
       kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*CD8_A - kb1_CD4_N_CD3_CD4_EM_CD28 - kf_CD4_N_CD3_CD4_EM_CD28*
      R_CD4_EM_CD28_tsAb*CD4_N - kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD28_tsAb*CD8_N - kf_CD4_EM_CD3_CD8_EM_CD28*
      R_CD4_EM_CD28_tsAb*CD8_EM - kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD28_tsAb*CD8_A - kf_CD4_EM_CD3_CD4_N_CD28*
      R_CD4_EM_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_CD4_EM_CD28 - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD4_EM - kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*CD4_EM - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD28_tsAb*
      CD4_A - kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD28_tsAb*MM - kb1_CD4_A_CD3_CD4_EM_CD28 - kf_CD4_A_CD3_CD4_EM_CD28*
      R_CD4_EM_CD28_tsAb*CD4_A - kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD28_tsAb*TRGT - kf_CD4_EM_CD3_MM_CD38*
      R_CD4_EM_CD28_tsAb*MM - kb2_CD4_EM_CD28_TRGT_CD38 - kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD28_tsAb*
      TRGT - kb2_CD4_EM_CD28_MM_CD38 - kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD28_tsAb*MM + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb + kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(R_CD4_A_CD28_tsAb) <- -kdeg_CD4_A*R_CD4_A_CD28_tsAb + kon_CD28*R_CD4_A_CD28*tsAb - koff_CD28*
      R_CD4_A_CD28_tsAb - kb1_CD8_N_CD3_CD4_A_CD28 - kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_N -
       kb1_CD8_EM_CD3_CD4_A_CD28 - kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_CD4_A_CD28 -
       kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_A - kb1_CD4_N_CD3_CD4_A_CD28 - kf_CD4_N_CD3_CD4_A_CD28*
      R_CD4_A_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_CD4_A_CD28 - kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*
      CD4_EM - kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD28_tsAb*CD8_N - kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD28_tsAb*
      CD8_EM - kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD28_tsAb*CD8_A - kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD28_tsAb*
      CD4_N - kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD28_tsAb*CD4_EM - kb1_CD4_A_CD3_CD4_A_CD28 - kf_CD4_A_CD3_CD4_A_CD28*
      R_CD4_A_CD28_tsAb*CD4_A - kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD4_A - kf_CD4_A_CD3_MM_CD28*
      R_CD4_A_CD28_tsAb*MM - kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD28_tsAb*TRGT - kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD28_tsAb*
      MM - kb2_CD4_A_CD28_TRGT_CD38 - kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD28_tsAb*TRGT - kb2_CD4_A_CD28_MM_CD38 -
       kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD28_tsAb*MM + kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb +
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb + kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb + kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb + kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb + kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb + kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb + kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb + kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb +
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb
    d/dt(R_MM_CD28_tsAb) <- -kdeg_MM*R_MM_CD28_tsAb + kon_CD28*R_MM_CD28*tsAb - koff_CD28*R_MM_CD28_tsAb -
       kb1_CD8_N_CD3_MM_CD28 - kf_CD8_N_CD3_MM_CD28*R_MM_CD28_tsAb*CD8_N - kb1_CD8_EM_CD3_MM_CD28 - kf_CD8_EM_CD3_MM_CD28*
      R_MM_CD28_tsAb*CD8_EM - kb1_CD8_A_CD3_MM_CD28 - kf_CD8_A_CD3_MM_CD28*R_MM_CD28_tsAb*CD8_A - kb1_CD4_N_CD3_MM_CD28 -
       kf_CD4_N_CD3_MM_CD28*R_MM_CD28_tsAb*CD4_N - kb1_CD4_EM_CD3_MM_CD28 - kf_CD4_EM_CD3_MM_CD28*R_MM_CD28_tsAb*
      CD4_EM - kb1_CD4_A_CD3_MM_CD28 - kf_CD4_A_CD3_MM_CD28*R_MM_CD28_tsAb*CD4_A - kf_CD8_N_CD3_MM_CD38*
      R_MM_CD28_tsAb*CD8_N - kf_CD8_EM_CD3_MM_CD38*R_MM_CD28_tsAb*CD8_EM - kf_CD8_A_CD3_MM_CD38*R_MM_CD28_tsAb*
      CD8_A - kf_CD4_N_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_N - kf_CD4_EM_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_EM -
       kf_CD4_A_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_A - kf_CD8_N_CD28_MM_CD38*R_MM_CD28_tsAb*CD8_N - kf_CD8_EM_CD28_MM_CD38*
      R_MM_CD28_tsAb*CD8_EM - kf_CD8_A_CD28_MM_CD38*R_MM_CD28_tsAb*CD8_A - kf_CD4_N_CD28_MM_CD38*R_MM_CD28_tsAb*
      CD4_N - kf_CD4_EM_CD28_MM_CD38*R_MM_CD28_tsAb*CD4_EM - kf_CD4_A_CD28_MM_CD38*R_MM_CD28_tsAb*CD4_A -
       kb2_MM_CD28_TRGT_CD38 - kf_MM_CD28_TRGT_CD38*R_MM_CD28_tsAb*TRGT - kb2_MM_CD28_MM_CD38 - kf_MM_CD28_MM_CD38*
      R_MM_CD28_tsAb*MM - kf_MM_CD28_MM_CD38*R_MM_CD28_tsAb*MM + kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb +
       kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb + kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb + kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb +
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb + kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb + kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb +
       kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb + kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb + kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb +
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb + kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb + kDis*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb +
       kDis*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb + kDis*S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb + kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb +
       kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb + kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb + kDis*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb + kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb + kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(R_TRGT_CD28_tsAb) <- -kdeg_TRGT*R_TRGT_CD28_tsAb + kon_CD28*R_TRGT_CD28*tsAb - koff_CD28*R_TRGT_CD28_tsAb -
       kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_N - kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_EM -
       kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_A - kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_N -
       kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_EM - kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_A -
       kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_N - kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_EM -
       kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_A - kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_N -
       kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_EM - kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_A -
       kf_MM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*MM + kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb +
       kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*
      S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb +
       kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3) <- +kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD3*MM - kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3 -
       kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3 - kon_CD3*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3) <- +kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD3*TRGT - kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3 - kon_CD3*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_A + kf_CD4_A_CD3_CD4_A_CD28*
      R_CD4_A_CD3*CD4_A + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3 + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 +
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3 + ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3 -
       kon_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3 - kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD3*CD4_EM + kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD3*CD4_N + kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD3*CD8_A + kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD3*CD8_EM + kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD3*CD8_N + kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3) <- +kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD3*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3 -
       kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3 + kact_EM*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3) <- +kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD3*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3 -
       kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3 + kact_EM*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3 + ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3 - kon_CD3*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3) <- +kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD3*TRGT - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3 + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3 + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3 - kon_CD3*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3
    d/dt(S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3) <- +kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD3*MM - kon_CD3*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3*
      tsAb + koff_CD3*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3) <- +kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD3*TRGT - kon_CD3*
      S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD3*CD4_A - kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_EM + kf_CD4_EM_CD3_CD4_EM_CD28*
      R_CD4_EM_CD3*CD4_EM - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 -
       kon_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 - kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_EM_CD3*CD4_N - kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD3*CD8_A - kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD3*CD8_EM - kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD3*CD8_N - kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD3*MM - kact_EM*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3 -
       kon_CD3*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD3*MM - kact_EM*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3 -
       kon_CD3*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3) <- +kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD3*TRGT - kact_EM*
      S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3 - kon_CD3*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3
    d/dt(S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3) <- +kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD3*MM - kon_CD3*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3*
      tsAb + koff_CD3*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3) <- +kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD3*TRGT - kon_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3*
      tsAb + koff_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD3*CD4_A - ka_CD4_N_CD3_CD4_A_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD3*CD4_EM - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3*CD4_N + kf_CD4_N_CD3_CD4_N_CD28*
      R_CD4_N_CD3*CD4_N - ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 - ka_CD4_N_CD3_CD4_N_CD28*
      S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD3*CD8_A - ka_CD4_N_CD3_CD8_A_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD3*CD8_EM - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD3*CD8_N - ka_CD4_N_CD3_CD8_N_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3) <- +kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD3*MM - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3) <- +kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD3*MM - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3) <- +kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD3*TRGT - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3 - kon_CD3*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3
    d/dt(S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3) <- +kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD3*MM - kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3 -
       kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3 - kon_CD3*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3) <- +kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD3*TRGT - kkillTRGT_CD8*
      S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3 - kon_CD3*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD3*CD4_A + kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD3*CD4_EM + kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD3*CD4_N + kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD3*CD8_A + kf_CD8_A_CD3_CD8_A_CD28*
      R_CD8_A_CD3*CD8_A + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3 + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 +
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3 + ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3 -
       kon_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3 - kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD3*CD8_EM + kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD3*CD8_N + kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3) <- +kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD3*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3 -
       kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3 + kact_EM*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3) <- +kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD3*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3 -
       kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3 + kact_EM*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3 + ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3 - kon_CD3*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3) <- +kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD3*TRGT - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3 + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3 + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3 - kon_CD3*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3
    d/dt(S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3) <- +kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD3*MM - kon_CD3*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3*
      tsAb + koff_CD3*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3) <- +kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD3*TRGT - kon_CD3*
      S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD3*CD4_A - kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD3*CD4_EM - kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD3*CD4_N - kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD3*CD8_A - kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_EM + kf_CD8_EM_CD3_CD8_EM_CD28*
      R_CD8_EM_CD3*CD8_EM - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 -
       kon_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 - kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD3*CD8_N - kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD3*MM - kact_EM*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3 -
       kon_CD3*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD3*MM - kact_EM*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3 -
       kon_CD3*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3) <- +kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD3*TRGT - kact_EM*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3 - kon_CD3*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3
    d/dt(S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3) <- +kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD3*MM - kon_CD3*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3*
      tsAb + koff_CD3*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3) <- +kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD3*TRGT - kon_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3*
      tsAb + koff_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD3*CD4_A - ka_CD8_N_CD3_CD4_A_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD8_N_CD3*CD4_EM - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD3*CD4_N - ka_CD8_N_CD3_CD4_N_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD3*CD8_A - ka_CD8_N_CD3_CD8_A_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD3*CD8_EM - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_N + kf_CD8_N_CD3_CD8_N_CD28*
      R_CD8_N_CD3*CD8_N - ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 - ka_CD8_N_CD3_CD8_N_CD28*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3) <- +kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD3*MM - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3) <- +kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD3*MM - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3) <- +kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD3*TRGT - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3 - kon_CD3*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3
    d/dt(S_MM_CD28_MM_CD38_R_MM_CD38) <- +kf_MM_CD28_MM_CD38*R_MM_CD38*MM + kf_MM_CD28_MM_CD38*R_MM_CD38*
      MM - kon_CD38*S_MM_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*
      S_MM_CD28_MM_CD38_R_MM_CD38 - kDis*S_MM_CD28_MM_CD38_R_MM_CD38
    d/dt(S_MM_CD28_TRGT_CD38_R_MM_CD38) <- +kf_MM_CD28_TRGT_CD38*R_MM_CD38*TRGT - kon_CD38*S_MM_CD28_TRGT_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_MM_CD38
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3 - kon_CD3*
      S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3 - kon_CD3*
      S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb - kDis*
      S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3 - kon_CD3*
      S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb - kDis*
      S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3 - kon_CD3*
      S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3 - kon_CD3*
      S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb - kDis*
      S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3 - kon_CD3*
      S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb - kDis*
      S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3
    d/dt(S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD3_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb - kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb + kon_CD3*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD3_tsAb*TRGT -
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb + kon_CD3*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_A +
       kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_A + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb +
       kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb +
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kon_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD3_tsAb*CD4_EM +
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD3_tsAb*CD4_N +
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD3_tsAb*CD8_A +
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD3_tsAb*CD8_EM +
       kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD3_tsAb*CD8_N +
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD3_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb + kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD3_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb + kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb) <- +kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD3_tsAb*TRGT - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3*tsAb - koff_CD3*
      S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD3_tsAb*MM + kon_CD3*
      S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3*tsAb - koff_CD3*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb - kDis*
      S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD3_tsAb*TRGT +
       kon_CD3*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3*tsAb - koff_CD3*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD3_tsAb*CD4_A -
       kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*
      CD4_EM + kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD4_EM - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_EM_CD3_tsAb*CD4_N -
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD3_tsAb*CD8_A -
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD3_tsAb*
      CD8_EM - kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD3_tsAb*CD8_N -
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD3_tsAb*MM - kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3*tsAb - koff_CD3*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD3_tsAb*MM - kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3*tsAb - koff_CD3*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD3_tsAb*TRGT -
       kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD3_tsAb*MM + kon_CD3*
      S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3*tsAb - koff_CD3*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD3_tsAb*TRGT +
       kon_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3*tsAb - koff_CD3*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD3_tsAb*CD4_A -
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD3_tsAb*CD4_EM -
       ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_N +
       kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_N - ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD3_tsAb*CD8_A -
       ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD3_tsAb*CD8_EM -
       ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD3_tsAb*CD8_N -
       ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD3_tsAb*MM - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3*tsAb - koff_CD3*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD3_tsAb*MM - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3*tsAb - koff_CD3*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb) <- +kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD3_tsAb*TRGT - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3*tsAb - koff_CD3*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD3_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD3_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb - kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb + kon_CD3*
      S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD3_tsAb*TRGT -
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb + kon_CD3*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD3_tsAb*CD4_A +
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD3_tsAb*CD4_EM +
       kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD3_tsAb*CD4_N +
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_A +
       kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_A + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb +
       kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb +
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kon_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD3_tsAb*CD8_EM +
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD3_tsAb*CD8_N +
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD3_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb + kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD3_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb + kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb) <- +kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD3_tsAb*TRGT - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3*tsAb - koff_CD3*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD3_tsAb*MM + kon_CD3*
      S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3*tsAb - koff_CD3*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb - kDis*
      S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD3_tsAb*TRGT +
       kon_CD3*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3*tsAb - koff_CD3*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD3_tsAb*CD4_A -
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD3_tsAb*
      CD4_EM - kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD3_tsAb*CD4_N -
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD3_tsAb*CD8_A -
       kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*
      CD8_EM + kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_EM - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD3_tsAb*CD8_N -
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD3_tsAb*MM - kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3*tsAb - koff_CD3*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD3_tsAb*MM - kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3*tsAb - koff_CD3*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD3_tsAb*TRGT -
       kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD3_tsAb*MM + kon_CD3*
      S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3*tsAb - koff_CD3*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD3_tsAb*TRGT +
       kon_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3*tsAb - koff_CD3*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD3_tsAb*CD4_A -
       ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD8_N_CD3_tsAb*CD4_EM -
       ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD3_tsAb*CD4_N -
       ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD3_tsAb*CD8_A -
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD3_tsAb*CD8_EM -
       ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_N +
       kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_N - ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD3_tsAb*MM - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3*tsAb - koff_CD3*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD3_tsAb*MM - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3*tsAb - koff_CD3*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb) <- +kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD3_tsAb*TRGT - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3*tsAb - koff_CD3*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD3_tsAb
    d/dt(S_MM_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_MM_CD28_MM_CD38*R_MM_CD38_tsAb*MM + kf_MM_CD28_MM_CD38*
      R_MM_CD38_tsAb*MM + kon_CD38*S_MM_CD28_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_MM_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb) <- +kf_MM_CD28_TRGT_CD38*R_MM_CD38_tsAb*TRGT + kon_CD38*
      S_MM_CD28_TRGT_CD38_R_MM_CD38*tsAb - koff_CD38*S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD3_tsAb +
       kon_CD3*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD3_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD3_tsAb +
       kon_CD3*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD3_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28) <- +kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD28*MM - kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28 -
       kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28 - kon_CD28*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28*tsAb +
       koff_CD28*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28) <- +kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD28*TRGT - kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28 - kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD28*CD4_A + kf_CD4_A_CD3_CD4_A_CD28*
      R_CD4_A_CD28*CD4_A + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28 + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 +
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28 + ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28 -
       kon_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28 - kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD28*CD4_EM + kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD28*CD4_N + kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD28*CD8_A + kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD28*CD8_EM + kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD28*CD8_N + kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28) <- +kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD28*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28 -
       kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28 + kact_EM*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28) <- +kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD28*MM - kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28 -
       kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28 + kact_EM*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28 + ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28) <- +kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD28*TRGT - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28 + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28 + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28 - kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28
    d/dt(S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28) <- +kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD28*MM - kon_CD28*
      S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28) <- +kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD28*TRGT - kon_CD28*
      S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD28*CD4_A - kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_EM +
       kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_EM - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 -
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb + koff_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 -
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_EM_CD28*CD4_N - kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD28*CD8_A - kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD28*CD8_EM -
       kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28*
      tsAb + koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD28*CD8_N - kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD28*MM - kact_EM*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28 -
       kon_CD28*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD28*MM - kact_EM*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28 -
       kon_CD28*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28) <- +kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD28*TRGT - kact_EM*
      S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28 - kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28
    d/dt(S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28) <- +kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD28*MM - kon_CD28*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28*
      tsAb + koff_CD28*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28) <- +kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD28*TRGT - kon_CD28*
      S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD28*CD4_A - ka_CD4_N_CD3_CD4_A_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD28*CD4_EM - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD28*CD4_N + kf_CD4_N_CD3_CD4_N_CD28*
      R_CD4_N_CD28*CD4_N - ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 - ka_CD4_N_CD3_CD4_N_CD28*
      S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD28*CD8_A - ka_CD4_N_CD3_CD8_A_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD28*CD8_EM - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD28*CD8_N - ka_CD4_N_CD3_CD8_N_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28) <- +kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD28*MM - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28) <- +kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD28*MM - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28) <- +kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD28*TRGT - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28 - kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28
    d/dt(S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28) <- +kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD28*MM - kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28 -
       kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28 - kon_CD28*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28*tsAb +
       koff_CD28*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28) <- +kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD28*TRGT - kkillTRGT_CD8*
      S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28 - kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD28*CD4_A + kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD28*CD4_EM + kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD28*CD4_N + kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28*CD8_A + kf_CD8_A_CD3_CD8_A_CD28*
      R_CD8_A_CD28*CD8_A + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28 + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 +
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28 + ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28 -
       kon_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28 - kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD28*CD8_EM + kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD28*CD8_N + kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28) <- +kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD28*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28 -
       kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28 + kact_EM*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28) <- +kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD28*MM - kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28 -
       kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28 + kact_EM*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28 + ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28) <- +kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD28*TRGT - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28 + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28 + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28 - kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28
    d/dt(S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28) <- +kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD28*MM - kon_CD28*
      S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28) <- +kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD28*TRGT - kon_CD28*
      S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD28*CD4_A - kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD28*CD4_EM -
       kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28*
      tsAb + koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD28*CD4_N - kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD28*CD8_A - kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_EM +
       kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_EM - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 -
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb + koff_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 -
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD28*CD8_N - kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD28*MM - kact_EM*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28 -
       kon_CD28*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD28*MM - kact_EM*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28 -
       kon_CD28*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28) <- +kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD28*TRGT - kact_EM*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28 - kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28
    d/dt(S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28) <- +kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD28*MM - kon_CD28*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28*
      tsAb + koff_CD28*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28) <- +kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD28*TRGT - kon_CD28*
      S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD28*CD4_A - ka_CD8_N_CD3_CD4_A_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD8_N_CD28*CD4_EM - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD28*CD4_N - ka_CD8_N_CD3_CD4_N_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD28*CD8_A - ka_CD8_N_CD3_CD8_A_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD28*CD8_EM - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_N + kf_CD8_N_CD3_CD8_N_CD28*
      R_CD8_N_CD28*CD8_N - ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 - ka_CD8_N_CD3_CD8_N_CD28*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28) <- +kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD28*MM - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28) <- +kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD28*MM - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28) <- +kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD28*TRGT - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28 - kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28
    d/dt(S_MM_CD28_MM_CD38_R_MM_CD28) <- +kf_MM_CD28_MM_CD38*R_MM_CD28*MM + kf_MM_CD28_MM_CD38*R_MM_CD28*
      MM - kon_CD28*S_MM_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*
      S_MM_CD28_MM_CD38_R_MM_CD28 - kDis*S_MM_CD28_MM_CD38_R_MM_CD28
    d/dt(S_MM_CD28_TRGT_CD38_R_MM_CD28) <- +kf_MM_CD28_TRGT_CD38*R_MM_CD28*TRGT - kon_CD28*S_MM_CD28_TRGT_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_MM_CD28
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28 - kon_CD28*
      S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28 - kon_CD28*
      S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28 - kon_CD28*
      S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28 - kon_CD28*
      S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28 - kon_CD28*
      S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28 - kon_CD28*
      S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28
    d/dt(S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD28_MM_CD38*R_CD4_A_CD28_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb - kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb + kon_CD28*
      S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD28_TRGT_CD38*R_CD4_A_CD28_tsAb*TRGT -
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb + kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD4_A +
       kf_CD4_A_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD4_A + kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb +
       kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb +
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kon_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_A_CD28_tsAb*CD4_EM +
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_A_CD28_tsAb*CD4_N +
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD4_A_CD28_tsAb*CD8_A +
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD4_A_CD28_tsAb*CD8_EM +
       kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD4_A_CD28_tsAb*CD8_N +
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_MM_CD28*R_CD4_A_CD28_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb + kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_MM_CD38*R_CD4_A_CD28_tsAb*MM - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb + kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb) <- +kf_CD4_A_CD3_TRGT_CD38*R_CD4_A_CD28_tsAb*TRGT -
       kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb +
       ka_CD4_N_CD3_TRGT_CD38*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD28_MM_CD38*R_CD4_EM_CD28_tsAb*MM +
       kon_CD28*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28*tsAb - koff_CD28*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD28_MM_CD38_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD28_TRGT_CD38*R_CD4_EM_CD28_tsAb*
      TRGT + kon_CD28*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28*tsAb - koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_EM_CD28_tsAb*
      CD4_A - kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD4_EM + kf_CD4_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*CD4_EM - kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_EM_CD28_tsAb*
      CD4_N - kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD4_EM_CD28_tsAb*
      CD8_A - kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD8_EM - kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD4_EM_CD28_tsAb*
      CD8_N - kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_MM_CD28*R_CD4_EM_CD28_tsAb*MM - kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28*tsAb - koff_CD28*
      S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_MM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_MM_CD38*R_CD4_EM_CD28_tsAb*MM - kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28*tsAb - koff_CD28*
      S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_MM_CD38_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_TRGT_CD38*R_CD4_EM_CD28_tsAb*TRGT -
       kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD28_MM_CD38*R_CD4_N_CD28_tsAb*MM + kon_CD28*
      S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28*tsAb - koff_CD28*S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb - kDis*
      S_CD4_N_CD28_MM_CD38_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD28_TRGT_CD38*R_CD4_N_CD28_tsAb*TRGT +
       kon_CD28*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28*tsAb - koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_N_CD28_TRGT_CD38_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_N_CD28_tsAb*CD4_A -
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_N_CD28_tsAb*CD4_EM -
       ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_N +
       kf_CD4_N_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_N - ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD4_N_CD28_tsAb*CD8_A -
       ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD4_N_CD28_tsAb*CD8_EM -
       ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD4_N_CD28_tsAb*CD8_N -
       ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_MM_CD28*R_CD4_N_CD28_tsAb*MM - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28*tsAb - koff_CD28*
      S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_MM_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_MM_CD38*R_CD4_N_CD28_tsAb*MM - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28*tsAb - koff_CD28*
      S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_MM_CD38_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb) <- +kf_CD4_N_CD3_TRGT_CD38*R_CD4_N_CD28_tsAb*TRGT -
       ka_CD4_N_CD3_TRGT_CD38*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_CD4_N_CD28_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD28_MM_CD38*R_CD8_A_CD28_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb - kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb + kon_CD28*
      S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD28_TRGT_CD38*R_CD8_A_CD28_tsAb*TRGT -
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb + kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD8_A_CD28_tsAb*CD4_A +
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD8_A_CD28_tsAb*CD4_EM +
       kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD8_A_CD28_tsAb*CD4_N +
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_A +
       kf_CD8_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_A + kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb +
       kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb +
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kon_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_A_CD28_tsAb*CD8_EM +
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_A_CD28_tsAb*CD8_N +
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_MM_CD28*R_CD8_A_CD28_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb + kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_MM_CD38*R_CD8_A_CD28_tsAb*MM - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb + kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb) <- +kf_CD8_A_CD3_TRGT_CD38*R_CD8_A_CD28_tsAb*TRGT -
       kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb +
       ka_CD8_N_CD3_TRGT_CD38*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD28_MM_CD38*R_CD8_EM_CD28_tsAb*MM +
       kon_CD28*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28*tsAb - koff_CD28*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD28_MM_CD38_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD28_TRGT_CD38*R_CD8_EM_CD28_tsAb*
      TRGT + kon_CD28*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28*tsAb - koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD8_EM_CD28_tsAb*
      CD4_A - kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD4_EM - kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD8_EM_CD28_tsAb*
      CD4_N - kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_EM_CD28_tsAb*
      CD8_A - kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD8_EM + kf_CD8_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*CD8_EM - kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_EM_CD28_tsAb*
      CD8_N - kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_MM_CD28*R_CD8_EM_CD28_tsAb*MM - kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28*tsAb - koff_CD28*
      S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_MM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_MM_CD38*R_CD8_EM_CD28_tsAb*MM - kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28*tsAb - koff_CD28*
      S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_MM_CD38_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_TRGT_CD38*R_CD8_EM_CD28_tsAb*TRGT -
       kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD28_MM_CD38*R_CD8_N_CD28_tsAb*MM + kon_CD28*
      S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28*tsAb - koff_CD28*S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb - kDis*
      S_CD8_N_CD28_MM_CD38_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD28_TRGT_CD38*R_CD8_N_CD28_tsAb*TRGT +
       kon_CD28*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28*tsAb - koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_N_CD28_TRGT_CD38_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD8_N_CD28_tsAb*CD4_A -
       ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD8_N_CD28_tsAb*CD4_EM -
       ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD8_N_CD28_tsAb*CD4_N -
       ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_N_CD28_tsAb*CD8_A -
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_N_CD28_tsAb*CD8_EM -
       ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_N +
       kf_CD8_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_N - ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_MM_CD28*R_CD8_N_CD28_tsAb*MM - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28*tsAb - koff_CD28*
      S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_MM_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_MM_CD38*R_CD8_N_CD28_tsAb*MM - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28*tsAb - koff_CD28*
      S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_MM_CD38_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb) <- +kf_CD8_N_CD3_TRGT_CD38*R_CD8_N_CD28_tsAb*TRGT -
       ka_CD8_N_CD3_TRGT_CD38*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_CD8_N_CD28_tsAb
    d/dt(S_MM_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_MM_CD28_MM_CD38*R_MM_CD28_tsAb*MM + kf_MM_CD28_MM_CD38*
      R_MM_CD28_tsAb*MM + kon_CD28*S_MM_CD28_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_MM_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb) <- +kf_MM_CD28_TRGT_CD38*R_MM_CD28_tsAb*TRGT + kon_CD28*
      S_MM_CD28_TRGT_CD38_R_MM_CD28*tsAb - koff_CD28*S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_CD4_A_CD28_tsAb +
       kon_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD28_MM_CD38_MUT_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_CD4_A_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_MM_CD28_MUT_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_CD4_A_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_MM_CD38_MUT_R_CD4_A_CD28_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_CD8_A_CD28_tsAb +
       kon_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD28_MM_CD38_MUT_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_CD8_A_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_MM_CD28_MUT_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_CD8_A_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_MM_CD38_MUT_R_CD8_A_CD28_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_R_MM_CD38) <- +kf_CD4_A_CD28_MM_CD38*R_MM_CD38*CD4_A - kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_MM_CD38 -
       kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD38 - kon_CD38*S_CD4_A_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*
      S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD38*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 -
       kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_N_CD3*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD3*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD3*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD4_A_CD3_MM_CD28_R_MM_CD38) <- +kf_CD4_A_CD3_MM_CD28*R_MM_CD38*CD4_A - kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_MM_CD38 -
       kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD38 + kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD38 + ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_MM_CD38) <- +kf_CD4_A_CD3_MM_CD38*R_MM_CD38*CD4_A - kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_MM_CD38 -
       kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD38 + kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD38 + ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD4_A_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD38*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38 + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38 + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_EM_CD28_MM_CD38_R_MM_CD38) <- +kf_CD4_EM_CD28_MM_CD38*R_MM_CD38*CD4_EM - kon_CD38*S_CD4_EM_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD38
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD38*CD4_EM - kon_CD38*
      S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 - kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD3*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 - kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD3*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD3*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 - kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD4_EM_CD3_MM_CD28_R_MM_CD38) <- +kf_CD4_EM_CD3_MM_CD28*R_MM_CD38*CD4_EM - kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD38 -
       kon_CD38*S_CD4_EM_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD38
    d/dt(S_CD4_EM_CD3_MM_CD38_R_MM_CD38) <- +kf_CD4_EM_CD3_MM_CD38*R_MM_CD38*CD4_EM - kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD38 -
       kon_CD38*S_CD4_EM_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD38
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD38*CD4_EM - kact_EM*
      S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_N_CD28_MM_CD38_R_MM_CD38) <- +kf_CD4_N_CD28_MM_CD38*R_MM_CD38*CD4_N - kon_CD38*S_CD4_N_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD38
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD38*CD4_N - kon_CD38*
      S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*
      S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD3*CD4_N - ka_CD4_N_CD3_CD4_A_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3 - kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD4_N - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD3*CD4_N - ka_CD4_N_CD3_CD8_A_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3 - kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD4_N - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD3*CD4_N - ka_CD4_N_CD3_CD8_N_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3 - kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD4_N_CD3_MM_CD28_R_MM_CD38) <- +kf_CD4_N_CD3_MM_CD28*R_MM_CD38*CD4_N - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD38
    d/dt(S_CD4_N_CD3_MM_CD38_R_MM_CD38) <- +kf_CD4_N_CD3_MM_CD38*R_MM_CD38*CD4_N - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD38
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD38*CD4_N - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_A_CD28_MM_CD38_R_MM_CD38) <- +kf_CD8_A_CD28_MM_CD38*R_MM_CD38*CD8_A - kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_MM_CD38 -
       kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD38 - kon_CD38*S_CD8_A_CD28_MM_CD38_R_MM_CD38*tsAb + koff_CD38*
      S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD38*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD3*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD3*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3 + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 -
       kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3 -
       kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD8_A_CD3_MM_CD28_R_MM_CD38) <- +kf_CD8_A_CD3_MM_CD28*R_MM_CD38*CD8_A - kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_MM_CD38 -
       kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD38 + kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD38 + ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_MM_CD38) <- +kf_CD8_A_CD3_MM_CD38*R_MM_CD38*CD8_A - kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_MM_CD38 -
       kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD38 + kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD38 + ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD8_A_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD38*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38 + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38 + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_EM_CD28_MM_CD38_R_MM_CD38) <- +kf_CD8_EM_CD28_MM_CD38*R_MM_CD38*CD8_EM - kon_CD38*S_CD8_EM_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD38
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD38*CD8_EM - kon_CD38*
      S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD3*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD3*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3 - kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD3*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3 - kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD3*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3 - kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb + koff_CD3*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3
    d/dt(S_CD8_EM_CD3_MM_CD28_R_MM_CD38) <- +kf_CD8_EM_CD3_MM_CD28*R_MM_CD38*CD8_EM - kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD38 -
       kon_CD38*S_CD8_EM_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD38
    d/dt(S_CD8_EM_CD3_MM_CD38_R_MM_CD38) <- +kf_CD8_EM_CD3_MM_CD38*R_MM_CD38*CD8_EM - kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD38 -
       kon_CD38*S_CD8_EM_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD38
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD38*CD8_EM - kact_EM*
      S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_N_CD28_MM_CD38_R_MM_CD38) <- +kf_CD8_N_CD28_MM_CD38*R_MM_CD38*CD8_N - kon_CD38*S_CD8_N_CD28_MM_CD38_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD38
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD38*CD8_N - kon_CD38*
      S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*
      S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD3*CD8_N - ka_CD8_N_CD3_CD4_A_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3 - kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3*CD8_N - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3 - kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD3*CD8_N - ka_CD8_N_CD3_CD4_N_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3 - kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD3*CD8_N - ka_CD8_N_CD3_CD8_A_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3 - kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3*CD8_N - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3 - kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb + koff_CD3*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3
    d/dt(S_CD8_N_CD3_MM_CD28_R_MM_CD38) <- +kf_CD8_N_CD3_MM_CD28*R_MM_CD38*CD8_N - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38*tsAb + koff_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD38
    d/dt(S_CD8_N_CD3_MM_CD38_R_MM_CD38) <- +kf_CD8_N_CD3_MM_CD38*R_MM_CD38*CD8_N - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38*tsAb + koff_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD38
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38) <- +kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD38*CD8_N - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38 - kon_CD38*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38*tsAb + koff_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38
    d/dt(S_MM_CD28_TRGT_CD38_R_TRGT_CD38) <- +kf_MM_CD28_TRGT_CD38*R_TRGT_CD38*MM - kon_CD38*S_MM_CD28_TRGT_CD38_R_TRGT_CD38*
      tsAb + koff_CD38*S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD38
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD38 - kon_CD38*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb - kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD38 - kon_CD38*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD38 - kon_CD38*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb - kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD38 - kon_CD38*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38*
      tsAb + koff_CD38*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38
    d/dt(S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_A_CD28_MM_CD38*R_MM_CD38_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb - kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD4_A_CD28_MM_CD38_R_MM_CD38*
      tsAb - koff_CD38*S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_A -
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38*
      tsAb - koff_CD38*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb - koff_CD3*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD4_A_CD3_MM_CD28*R_MM_CD38_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb + kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb +
       ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38*
      tsAb - koff_CD38*S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_A_CD3_MM_CD38*R_MM_CD38_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb + kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb +
       ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD4_A_CD3_MM_CD38_R_MM_CD38*
      tsAb - koff_CD38*S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_EM_CD28_MM_CD38*R_MM_CD38_tsAb*CD4_EM + kon_CD38*
      S_CD4_EM_CD28_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_EM +
       kon_CD38*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*
      CD4_EM - kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD4_EM_CD3_MM_CD28*R_MM_CD38_tsAb*CD4_EM - kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD4_EM_CD3_MM_CD28_R_MM_CD38*tsAb - koff_CD38*
      S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb - kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_EM_CD3_MM_CD38*R_MM_CD38_tsAb*CD4_EM - kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD4_EM_CD3_MM_CD38_R_MM_CD38*tsAb - koff_CD38*
      S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38*
      tsAb - koff_CD38*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_N_CD28_MM_CD38*R_MM_CD38_tsAb*CD4_N + kon_CD38*
      S_CD4_N_CD28_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_N +
       kon_CD38*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD4_N -
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD4_N -
       ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD4_N -
       ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD4_N -
       ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD4_N -
       ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD4_N_CD3_MM_CD28*R_MM_CD38_tsAb*CD4_N - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38*tsAb - koff_CD38*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD4_N_CD3_MM_CD38*R_MM_CD38_tsAb*CD4_N - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD4_N - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_A_CD28_MM_CD38*R_MM_CD38_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb - kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD8_A_CD28_MM_CD38_R_MM_CD38*
      tsAb - koff_CD38*S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_A -
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38*
      tsAb - koff_CD38*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb +
       kon_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3*tsAb - koff_CD3*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD8_A_CD3_MM_CD28*R_MM_CD38_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb + kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb +
       ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38*
      tsAb - koff_CD38*S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_A_CD3_MM_CD38*R_MM_CD38_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb + kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb +
       ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD8_A_CD3_MM_CD38_R_MM_CD38*
      tsAb - koff_CD38*S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_EM_CD28_MM_CD38*R_MM_CD38_tsAb*CD8_EM + kon_CD38*
      S_CD8_EM_CD28_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_EM +
       kon_CD38*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*
      CD8_EM - kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD3_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb + kon_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3*
      tsAb - koff_CD3*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD3_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD8_EM_CD3_MM_CD28*R_MM_CD38_tsAb*CD8_EM - kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD8_EM_CD3_MM_CD28_R_MM_CD38*tsAb - koff_CD38*
      S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb - kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_EM_CD3_MM_CD38*R_MM_CD38_tsAb*CD8_EM - kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD8_EM_CD3_MM_CD38_R_MM_CD38*tsAb - koff_CD38*
      S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38*
      tsAb - koff_CD38*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_N_CD28_MM_CD38*R_MM_CD38_tsAb*CD8_N + kon_CD38*
      S_CD8_N_CD28_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_N +
       kon_CD38*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb -
       kDis*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD3_tsAb*CD8_N -
       ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD3_tsAb*CD8_N -
       ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD3_tsAb*CD8_N -
       ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD3_tsAb*CD8_N -
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD3_tsAb
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD3_tsAb*CD8_N -
       ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb + kon_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3*
      tsAb - koff_CD3*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD3_tsAb
    d/dt(S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb) <- +kf_CD8_N_CD3_MM_CD28*R_MM_CD38_tsAb*CD8_N - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38*tsAb - koff_CD38*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD38_tsAb
    d/dt(S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb) <- +kf_CD8_N_CD3_MM_CD38*R_MM_CD38_tsAb*CD8_N - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38*tsAb - koff_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD38_tsAb
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD38_tsAb*CD8_N - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb + kon_CD38*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) <- +kf_MM_CD28_TRGT_CD38*R_TRGT_CD38_tsAb*MM + kon_CD38*
      S_MM_CD28_TRGT_CD38_R_TRGT_CD38*tsAb - koff_CD38*S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD38_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb + kon_CD38*
      S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*
      S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*
      S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb - kDis*
      S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*
      S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*
      S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb + kon_CD38*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb + kon_CD38*
      S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb - kDis*
      S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD38_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb + kon_CD38*
      S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38*tsAb - koff_CD38*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb - kDis*
      S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD38_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_R_MM_CD28) <- +kf_CD4_A_CD28_MM_CD38*R_MM_CD28*CD4_A - kkillMM_CD4*S_CD4_A_CD28_MM_CD38_R_MM_CD28 -
       kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD28 - kon_CD28*S_CD4_A_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*
      S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD28*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 -
       kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_N_CD28*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD28*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD28*CD4_A + kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD4_A_CD3_MM_CD28_R_MM_CD28) <- +kf_CD4_A_CD3_MM_CD28*R_MM_CD28*CD4_A - kkillMM_CD4*S_CD4_A_CD3_MM_CD28_R_MM_CD28 -
       kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD28 + kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD28 + ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_MM_CD28) <- +kf_CD4_A_CD3_MM_CD38*R_MM_CD28*CD4_A - kkillMM_CD4*S_CD4_A_CD3_MM_CD38_R_MM_CD28 -
       kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD28 + kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD28 + ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD28*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28 + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28 + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_EM_CD28_MM_CD38_R_MM_CD28) <- +kf_CD4_EM_CD28_MM_CD38*R_MM_CD28*CD4_EM - kon_CD28*S_CD4_EM_CD28_MM_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD28
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD28*CD4_EM - kon_CD28*
      S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD28*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 - kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD28*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 - kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD28*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb + koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD28*CD4_EM - kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 - kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD4_EM_CD3_MM_CD28_R_MM_CD28) <- +kf_CD4_EM_CD3_MM_CD28*R_MM_CD28*CD4_EM - kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD28 -
       kon_CD28*S_CD4_EM_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD28
    d/dt(S_CD4_EM_CD3_MM_CD38_R_MM_CD28) <- +kf_CD4_EM_CD3_MM_CD38*R_MM_CD28*CD4_EM - kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD28 -
       kon_CD28*S_CD4_EM_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD28
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD28*CD4_EM - kact_EM*
      S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_N_CD28_MM_CD38_R_MM_CD28) <- +kf_CD4_N_CD28_MM_CD38*R_MM_CD28*CD4_N - kon_CD28*S_CD4_N_CD28_MM_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD28
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD28*CD4_N - kon_CD28*
      S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*
      S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD28*CD4_N - ka_CD4_N_CD3_CD4_A_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28 - kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD4_N - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD28*CD4_N - ka_CD4_N_CD3_CD8_A_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28 - kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD4_N - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD28*CD4_N - ka_CD4_N_CD3_CD8_N_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28 - kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD4_N_CD3_MM_CD28_R_MM_CD28) <- +kf_CD4_N_CD3_MM_CD28*R_MM_CD28*CD4_N - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD28
    d/dt(S_CD4_N_CD3_MM_CD38_R_MM_CD28) <- +kf_CD4_N_CD3_MM_CD38*R_MM_CD28*CD4_N - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD28
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD28*CD4_N - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_A_CD28_MM_CD38_R_MM_CD28) <- +kf_CD8_A_CD28_MM_CD38*R_MM_CD28*CD8_A - kkillMM_CD8*S_CD8_A_CD28_MM_CD38_R_MM_CD28 -
       kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD28 - kon_CD28*S_CD8_A_CD28_MM_CD38_R_MM_CD28*tsAb + koff_CD28*
      S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD28*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD28*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD28*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28 + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 -
       kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_A + kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28 -
       kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD8_A_CD3_MM_CD28_R_MM_CD28) <- +kf_CD8_A_CD3_MM_CD28*R_MM_CD28*CD8_A - kkillMM_CD8*S_CD8_A_CD3_MM_CD28_R_MM_CD28 -
       kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD28 + kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD28 + ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_MM_CD28) <- +kf_CD8_A_CD3_MM_CD38*R_MM_CD28*CD8_A - kkillMM_CD8*S_CD8_A_CD3_MM_CD38_R_MM_CD28 -
       kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD28 + kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD28 + ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD28*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28 + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28 + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_EM_CD28_MM_CD38_R_MM_CD28) <- +kf_CD8_EM_CD28_MM_CD38*R_MM_CD28*CD8_EM - kon_CD28*S_CD8_EM_CD28_MM_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD28
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD28*CD8_EM - kon_CD28*
      S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD28*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb + koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD28*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28 - kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD28*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28 - kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD28*CD8_EM - kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28 - kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28
    d/dt(S_CD8_EM_CD3_MM_CD28_R_MM_CD28) <- +kf_CD8_EM_CD3_MM_CD28*R_MM_CD28*CD8_EM - kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD28 -
       kon_CD28*S_CD8_EM_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD28
    d/dt(S_CD8_EM_CD3_MM_CD38_R_MM_CD28) <- +kf_CD8_EM_CD3_MM_CD38*R_MM_CD28*CD8_EM - kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD28 -
       kon_CD28*S_CD8_EM_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD28
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD28*CD8_EM - kact_EM*
      S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_N_CD28_MM_CD38_R_MM_CD28) <- +kf_CD8_N_CD28_MM_CD38*R_MM_CD28*CD8_N - kon_CD28*S_CD8_N_CD28_MM_CD38_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD28
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD28*CD8_N - kon_CD28*
      S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*
      S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD28*CD8_N - ka_CD8_N_CD3_CD4_A_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28 - kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28*CD8_N - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28 - kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD28*CD8_N - ka_CD8_N_CD3_CD4_N_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28 - kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD28*CD8_N - ka_CD8_N_CD3_CD8_A_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28 - kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28*CD8_N - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28 - kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28
    d/dt(S_CD8_N_CD3_MM_CD28_R_MM_CD28) <- +kf_CD8_N_CD3_MM_CD28*R_MM_CD28*CD8_N - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD28
    d/dt(S_CD8_N_CD3_MM_CD38_R_MM_CD28) <- +kf_CD8_N_CD3_MM_CD38*R_MM_CD28*CD8_N - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28*tsAb + koff_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD28
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28) <- +kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD28*CD8_N - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28 - kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28*tsAb + koff_CD28*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28
    d/dt(S_MM_CD28_TRGT_CD38_R_TRGT_CD28) <- +kf_MM_CD28_TRGT_CD38*R_TRGT_CD28*MM - kon_CD28*S_MM_CD28_TRGT_CD38_R_TRGT_CD28*
      tsAb + koff_CD28*S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD28
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD28 - kon_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb - kDis*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD28 - kon_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb - kDis*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD28 - kon_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28*
      tsAb + koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28
    d/dt(S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_A_CD28_MM_CD38*R_MM_CD28_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb - kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD4_A_CD28_MM_CD38_R_MM_CD28*
      tsAb - koff_CD28*S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_A_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_A -
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28*
      tsAb - koff_CD28*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD4_A + kact_EM*S_CD4_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_A_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD4_A + kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD4_A_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD4_A +
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb - koff_CD28*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD4_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD4_A_CD3_MM_CD28*R_MM_CD28_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb + kact_EM*S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb +
       ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28*
      tsAb - koff_CD28*S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_A_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_A - kkillMM_CD4*
      S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb - kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb + kact_EM*S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb +
       ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28*
      tsAb - koff_CD28*S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_A_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_A - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*
      S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_EM_CD28_MM_CD38*R_MM_CD28_tsAb*CD4_EM + kon_CD28*
      S_CD4_EM_CD28_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD4_EM_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_EM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_EM +
       kon_CD28*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD4_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_EM_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD4_EM_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD4_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD4_EM - kact_EM*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD4_EM_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD4_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD4_EM_CD3_MM_CD28*R_MM_CD28_tsAb*CD4_EM - kact_EM*
      S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_MM_CD28_R_MM_CD28*tsAb - koff_CD28*
      S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb - kDis*S_CD4_EM_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_EM_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_EM - kact_EM*
      S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_MM_CD38_R_MM_CD28*tsAb - koff_CD28*
      S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD4_EM_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_EM_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_EM -
       kact_EM*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28*
      tsAb - koff_CD28*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_N_CD28_MM_CD38*R_MM_CD28_tsAb*CD4_N + kon_CD28*
      S_CD4_N_CD28_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD4_N_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_N_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_N +
       kon_CD28*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD4_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD4_N_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD4_N -
       ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD4_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD4_N - ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD4_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD4_N -
       ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD4_N - ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD4_N_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD4_N -
       ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD4_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD4_N_CD3_MM_CD28*R_MM_CD28_tsAb*CD4_N - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28*tsAb - koff_CD28*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD4_N_CD3_MM_CD38*R_MM_CD28_tsAb*CD4_N - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD4_N_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD4_N_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD4_N - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*
      S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD4_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_A_CD28_MM_CD38*R_MM_CD28_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb - kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD8_A_CD28_MM_CD38_R_MM_CD28*
      tsAb - koff_CD28*S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_A_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_A -
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28*
      tsAb - koff_CD28*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD8_A + kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD8_A_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_A_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD8_A + kact_EM*S_CD8_EM_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_A_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_A +
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb +
       kon_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28*tsAb - koff_CD28*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb -
       kDis*S_CD8_A_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD8_A_CD3_MM_CD28*R_MM_CD28_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb + kact_EM*S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb +
       ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28*
      tsAb - koff_CD28*S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_A_CD3_MM_CD38*R_MM_CD28_tsAb*CD8_A - kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb - kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb + kact_EM*S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb +
       ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28*
      tsAb - koff_CD28*S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_A_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_A - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*
      S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_EM_CD28_MM_CD38*R_MM_CD28_tsAb*CD8_EM + kon_CD28*
      S_CD8_EM_CD28_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD8_EM_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_EM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_EM +
       kon_CD28*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD8_EM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD8_EM - kact_EM*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD8_EM_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD8_EM_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_EM_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb) <- +kf_CD8_EM_CD3_CD8_N_CD28*R_CD8_N_CD28_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb - kDis*S_CD8_EM_CD3_CD8_N_CD28_R_CD8_N_CD28_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD8_EM_CD3_MM_CD28*R_MM_CD28_tsAb*CD8_EM - kact_EM*
      S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_MM_CD28_R_MM_CD28*tsAb - koff_CD28*
      S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb - kDis*S_CD8_EM_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_EM_CD3_MM_CD38*R_MM_CD28_tsAb*CD8_EM - kact_EM*
      S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_MM_CD38_R_MM_CD28*tsAb - koff_CD28*
      S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD8_EM_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_EM_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_EM -
       kact_EM*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28*
      tsAb - koff_CD28*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_EM_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_N_CD28_MM_CD38*R_MM_CD28_tsAb*CD8_N + kon_CD28*
      S_CD8_N_CD28_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb - kDis*S_CD8_N_CD28_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_N_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_N +
       kon_CD28*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb -
       kDis*S_CD8_N_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_A_CD28*R_CD4_A_CD28_tsAb*CD8_N -
       ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_A_CD28_R_CD4_A_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_EM_CD28*R_CD4_EM_CD28_tsAb*
      CD8_N - ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_EM_CD28_R_CD4_EM_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb) <- +kf_CD8_N_CD3_CD4_N_CD28*R_CD4_N_CD28_tsAb*CD8_N -
       ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb - kDis*S_CD8_N_CD3_CD4_N_CD28_R_CD4_N_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb) <- +kf_CD8_N_CD3_CD8_A_CD28*R_CD8_A_CD28_tsAb*CD8_N -
       ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_A_CD28_R_CD8_A_CD28_tsAb
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb) <- +kf_CD8_N_CD3_CD8_EM_CD28*R_CD8_EM_CD28_tsAb*
      CD8_N - ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb + kon_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28*
      tsAb - koff_CD28*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb - kDis*S_CD8_N_CD3_CD8_EM_CD28_R_CD8_EM_CD28_tsAb
    d/dt(S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb) <- +kf_CD8_N_CD3_MM_CD28*R_MM_CD28_tsAb*CD8_N - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28*tsAb - koff_CD28*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD28_R_MM_CD28_tsAb
    d/dt(S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb) <- +kf_CD8_N_CD3_MM_CD38*R_MM_CD28_tsAb*CD8_N - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28*tsAb - koff_CD28*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb -
       kDis*S_CD8_N_CD3_MM_CD38_R_MM_CD28_tsAb
    d/dt(S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_CD8_N_CD3_TRGT_CD38*R_TRGT_CD28_tsAb*CD8_N - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb + kon_CD28*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*
      S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_CD8_N_CD3_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) <- +kf_MM_CD28_TRGT_CD38*R_TRGT_CD28_tsAb*MM + kon_CD28*
      S_MM_CD28_TRGT_CD38_R_TRGT_CD28*tsAb - koff_CD28*S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb - kDis*S_MM_CD28_TRGT_CD38_R_TRGT_CD28_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb + kon_CD28*
      S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*
      S_CD4_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*
      S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb - kDis*
      S_CD4_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*
      S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*
      S_CD4_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb + kon_CD28*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*
      S_CD8_A_CD28_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb + kon_CD28*
      S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb - kDis*
      S_CD8_A_CD3_MM_CD28_MUT_R_MM_CD28_tsAb
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb + kon_CD28*
      S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28*tsAb - koff_CD28*S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb - kDis*
      S_CD8_A_CD3_MM_CD38_MUT_R_MM_CD28_tsAb
    d/dt(S_CD4_A_CD28_MM_CD38_Br) <- +kb1_CD4_A_CD28_MM_CD38 + kb2_CD4_A_CD28_MM_CD38 - kkillMM_CD4*S_CD4_A_CD28_MM_CD38_Br -
       kmut_SYN*S_CD4_A_CD28_MM_CD38_Br
    d/dt(S_CD4_A_CD28_TRGT_CD38_Br) <- +kb1_CD4_A_CD28_TRGT_CD38 + kb2_CD4_A_CD28_TRGT_CD38 - kkillTRGT_CD4*
      S_CD4_A_CD28_TRGT_CD38_Br - kDis*S_CD4_A_CD28_TRGT_CD38_Br
    d/dt(S_CD4_A_CD3_CD4_A_CD28_Br) <- +kb1_CD4_A_CD3_CD4_A_CD28 + kb2_CD4_A_CD3_CD4_A_CD28 + kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_Br + ka_CD4_N_CD3_CD4_A_CD28*S_CD4_N_CD3_CD4_A_CD28_Br - kDis*S_CD4_A_CD3_CD4_A_CD28_Br
    d/dt(S_CD4_A_CD3_CD4_EM_CD28_Br) <- +kb1_CD4_A_CD3_CD4_EM_CD28 + kb2_CD4_A_CD3_CD4_EM_CD28 + kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_Br + ka_CD4_N_CD3_CD4_EM_CD28*S_CD4_N_CD3_CD4_EM_CD28_Br - kDis*S_CD4_A_CD3_CD4_EM_CD28_Br
    d/dt(S_CD4_A_CD3_CD4_N_CD28_Br) <- +kb1_CD4_A_CD3_CD4_N_CD28 + kb2_CD4_A_CD3_CD4_N_CD28 + kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_Br + ka_CD4_N_CD3_CD4_N_CD28*S_CD4_N_CD3_CD4_N_CD28_Br - kDis*S_CD4_A_CD3_CD4_N_CD28_Br
    d/dt(S_CD4_A_CD3_CD8_A_CD28_Br) <- +kb1_CD4_A_CD3_CD8_A_CD28 + kb2_CD4_A_CD3_CD8_A_CD28 + kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_Br + ka_CD4_N_CD3_CD8_A_CD28*S_CD4_N_CD3_CD8_A_CD28_Br - kDis*S_CD4_A_CD3_CD8_A_CD28_Br
    d/dt(S_CD4_A_CD3_CD8_EM_CD28_Br) <- +kb1_CD4_A_CD3_CD8_EM_CD28 + kb2_CD4_A_CD3_CD8_EM_CD28 + kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_Br + ka_CD4_N_CD3_CD8_EM_CD28*S_CD4_N_CD3_CD8_EM_CD28_Br - kDis*S_CD4_A_CD3_CD8_EM_CD28_Br
    d/dt(S_CD4_A_CD3_CD8_N_CD28_Br) <- +kb1_CD4_A_CD3_CD8_N_CD28 + kb2_CD4_A_CD3_CD8_N_CD28 + kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_Br + ka_CD4_N_CD3_CD8_N_CD28*S_CD4_N_CD3_CD8_N_CD28_Br - kDis*S_CD4_A_CD3_CD8_N_CD28_Br
    d/dt(S_CD4_A_CD3_MM_CD28_Br) <- +kb1_CD4_A_CD3_MM_CD28 + kb2_CD4_A_CD3_MM_CD28 - kkillMM_CD4*S_CD4_A_CD3_MM_CD28_Br -
       kmut_SYN*S_CD4_A_CD3_MM_CD28_Br + kact_EM*S_CD4_EM_CD3_MM_CD28_Br + ka_CD4_N_CD3_MM_CD28*S_CD4_N_CD3_MM_CD28_Br
    d/dt(S_CD4_A_CD3_MM_CD38_Br) <- +kb1_CD4_A_CD3_MM_CD38 + kb2_CD4_A_CD3_MM_CD38 - kkillMM_CD4*S_CD4_A_CD3_MM_CD38_Br -
       kmut_SYN*S_CD4_A_CD3_MM_CD38_Br + kact_EM*S_CD4_EM_CD3_MM_CD38_Br + ka_CD4_N_CD3_MM_CD38*S_CD4_N_CD3_MM_CD38_Br
    d/dt(S_CD4_A_CD3_TRGT_CD38_Br) <- +kb1_CD4_A_CD3_TRGT_CD38 + kb2_CD4_A_CD3_TRGT_CD38 - kkillTRGT_CD4*
      S_CD4_A_CD3_TRGT_CD38_Br + kact_EM*S_CD4_EM_CD3_TRGT_CD38_Br + ka_CD4_N_CD3_TRGT_CD38*S_CD4_N_CD3_TRGT_CD38_Br -
       kDis*S_CD4_A_CD3_TRGT_CD38_Br
    d/dt(S_CD4_EM_CD28_MM_CD38_Br) <- +kb1_CD4_EM_CD28_MM_CD38 + kb2_CD4_EM_CD28_MM_CD38 - kDis*S_CD4_EM_CD28_MM_CD38_Br
    d/dt(S_CD4_EM_CD28_TRGT_CD38_Br) <- +kb1_CD4_EM_CD28_TRGT_CD38 + kb2_CD4_EM_CD28_TRGT_CD38 - kDis*
      S_CD4_EM_CD28_TRGT_CD38_Br
    d/dt(S_CD4_EM_CD3_CD4_A_CD28_Br) <- +kb1_CD4_EM_CD3_CD4_A_CD28 + kb2_CD4_EM_CD3_CD4_A_CD28 - kact_EM*
      S_CD4_EM_CD3_CD4_A_CD28_Br - kDis*S_CD4_EM_CD3_CD4_A_CD28_Br
    d/dt(S_CD4_EM_CD3_CD4_EM_CD28_Br) <- +kb1_CD4_EM_CD3_CD4_EM_CD28 + kb2_CD4_EM_CD3_CD4_EM_CD28 - kact_EM*
      S_CD4_EM_CD3_CD4_EM_CD28_Br - kDis*S_CD4_EM_CD3_CD4_EM_CD28_Br
    d/dt(S_CD4_EM_CD3_CD4_N_CD28_Br) <- +kb1_CD4_EM_CD3_CD4_N_CD28 + kb2_CD4_EM_CD3_CD4_N_CD28 - kact_EM*
      S_CD4_EM_CD3_CD4_N_CD28_Br - kDis*S_CD4_EM_CD3_CD4_N_CD28_Br
    d/dt(S_CD4_EM_CD3_CD8_A_CD28_Br) <- +kb1_CD4_EM_CD3_CD8_A_CD28 + kb2_CD4_EM_CD3_CD8_A_CD28 - kact_EM*
      S_CD4_EM_CD3_CD8_A_CD28_Br - kDis*S_CD4_EM_CD3_CD8_A_CD28_Br
    d/dt(S_CD4_EM_CD3_CD8_EM_CD28_Br) <- +kb1_CD4_EM_CD3_CD8_EM_CD28 + kb2_CD4_EM_CD3_CD8_EM_CD28 - kact_EM*
      S_CD4_EM_CD3_CD8_EM_CD28_Br - kDis*S_CD4_EM_CD3_CD8_EM_CD28_Br
    d/dt(S_CD4_EM_CD3_CD8_N_CD28_Br) <- +kb1_CD4_EM_CD3_CD8_N_CD28 + kb2_CD4_EM_CD3_CD8_N_CD28 - kact_EM*
      S_CD4_EM_CD3_CD8_N_CD28_Br - kDis*S_CD4_EM_CD3_CD8_N_CD28_Br
    d/dt(S_CD4_EM_CD3_MM_CD28_Br) <- +kb1_CD4_EM_CD3_MM_CD28 + kb2_CD4_EM_CD3_MM_CD28 - kact_EM*S_CD4_EM_CD3_MM_CD28_Br -
       kDis*S_CD4_EM_CD3_MM_CD28_Br
    d/dt(S_CD4_EM_CD3_MM_CD38_Br) <- +kb1_CD4_EM_CD3_MM_CD38 + kb2_CD4_EM_CD3_MM_CD38 - kact_EM*S_CD4_EM_CD3_MM_CD38_Br -
       kDis*S_CD4_EM_CD3_MM_CD38_Br
    d/dt(S_CD4_EM_CD3_TRGT_CD38_Br) <- +kb1_CD4_EM_CD3_TRGT_CD38 + kb2_CD4_EM_CD3_TRGT_CD38 - kact_EM*
      S_CD4_EM_CD3_TRGT_CD38_Br - kDis*S_CD4_EM_CD3_TRGT_CD38_Br
    d/dt(S_CD4_N_CD28_MM_CD38_Br) <- +kb1_CD4_N_CD28_MM_CD38 + kb2_CD4_N_CD28_MM_CD38 - kDis*S_CD4_N_CD28_MM_CD38_Br
    d/dt(S_CD4_N_CD28_TRGT_CD38_Br) <- +kb1_CD4_N_CD28_TRGT_CD38 + kb2_CD4_N_CD28_TRGT_CD38 - kDis*S_CD4_N_CD28_TRGT_CD38_Br
    d/dt(S_CD4_N_CD3_CD4_A_CD28_Br) <- +kb1_CD4_N_CD3_CD4_A_CD28 + kb2_CD4_N_CD3_CD4_A_CD28 - ka_CD4_N_CD3_CD4_A_CD28*
      S_CD4_N_CD3_CD4_A_CD28_Br - kDis*S_CD4_N_CD3_CD4_A_CD28_Br
    d/dt(S_CD4_N_CD3_CD4_EM_CD28_Br) <- +kb1_CD4_N_CD3_CD4_EM_CD28 + kb2_CD4_N_CD3_CD4_EM_CD28 - ka_CD4_N_CD3_CD4_EM_CD28*
      S_CD4_N_CD3_CD4_EM_CD28_Br - kDis*S_CD4_N_CD3_CD4_EM_CD28_Br
    d/dt(S_CD4_N_CD3_CD4_N_CD28_Br) <- +kb1_CD4_N_CD3_CD4_N_CD28 + kb2_CD4_N_CD3_CD4_N_CD28 - ka_CD4_N_CD3_CD4_N_CD28*
      S_CD4_N_CD3_CD4_N_CD28_Br - kDis*S_CD4_N_CD3_CD4_N_CD28_Br
    d/dt(S_CD4_N_CD3_CD8_A_CD28_Br) <- +kb1_CD4_N_CD3_CD8_A_CD28 + kb2_CD4_N_CD3_CD8_A_CD28 - ka_CD4_N_CD3_CD8_A_CD28*
      S_CD4_N_CD3_CD8_A_CD28_Br - kDis*S_CD4_N_CD3_CD8_A_CD28_Br
    d/dt(S_CD4_N_CD3_CD8_EM_CD28_Br) <- +kb1_CD4_N_CD3_CD8_EM_CD28 + kb2_CD4_N_CD3_CD8_EM_CD28 - ka_CD4_N_CD3_CD8_EM_CD28*
      S_CD4_N_CD3_CD8_EM_CD28_Br - kDis*S_CD4_N_CD3_CD8_EM_CD28_Br
    d/dt(S_CD4_N_CD3_CD8_N_CD28_Br) <- +kb1_CD4_N_CD3_CD8_N_CD28 + kb2_CD4_N_CD3_CD8_N_CD28 - ka_CD4_N_CD3_CD8_N_CD28*
      S_CD4_N_CD3_CD8_N_CD28_Br - kDis*S_CD4_N_CD3_CD8_N_CD28_Br
    d/dt(S_CD4_N_CD3_MM_CD28_Br) <- +kb1_CD4_N_CD3_MM_CD28 + kb2_CD4_N_CD3_MM_CD28 - ka_CD4_N_CD3_MM_CD28*
      S_CD4_N_CD3_MM_CD28_Br - kDis*S_CD4_N_CD3_MM_CD28_Br
    d/dt(S_CD4_N_CD3_MM_CD38_Br) <- +kb1_CD4_N_CD3_MM_CD38 + kb2_CD4_N_CD3_MM_CD38 - ka_CD4_N_CD3_MM_CD38*
      S_CD4_N_CD3_MM_CD38_Br - kDis*S_CD4_N_CD3_MM_CD38_Br
    d/dt(S_CD4_N_CD3_TRGT_CD38_Br) <- +kb1_CD4_N_CD3_TRGT_CD38 + kb2_CD4_N_CD3_TRGT_CD38 - ka_CD4_N_CD3_TRGT_CD38*
      S_CD4_N_CD3_TRGT_CD38_Br - kDis*S_CD4_N_CD3_TRGT_CD38_Br
    d/dt(S_CD8_A_CD28_MM_CD38_Br) <- +kb1_CD8_A_CD28_MM_CD38 + kb2_CD8_A_CD28_MM_CD38 - kkillMM_CD8*S_CD8_A_CD28_MM_CD38_Br -
       kmut_SYN*S_CD8_A_CD28_MM_CD38_Br
    d/dt(S_CD8_A_CD28_TRGT_CD38_Br) <- +kb1_CD8_A_CD28_TRGT_CD38 + kb2_CD8_A_CD28_TRGT_CD38 - kkillTRGT_CD8*
      S_CD8_A_CD28_TRGT_CD38_Br - kDis*S_CD8_A_CD28_TRGT_CD38_Br
    d/dt(S_CD8_A_CD3_CD4_A_CD28_Br) <- +kb1_CD8_A_CD3_CD4_A_CD28 + kb2_CD8_A_CD3_CD4_A_CD28 + kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_Br + ka_CD8_N_CD3_CD4_A_CD28*S_CD8_N_CD3_CD4_A_CD28_Br - kDis*S_CD8_A_CD3_CD4_A_CD28_Br
    d/dt(S_CD8_A_CD3_CD4_EM_CD28_Br) <- +kb1_CD8_A_CD3_CD4_EM_CD28 + kb2_CD8_A_CD3_CD4_EM_CD28 + kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_Br + ka_CD8_N_CD3_CD4_EM_CD28*S_CD8_N_CD3_CD4_EM_CD28_Br - kDis*S_CD8_A_CD3_CD4_EM_CD28_Br
    d/dt(S_CD8_A_CD3_CD4_N_CD28_Br) <- +kb1_CD8_A_CD3_CD4_N_CD28 + kb2_CD8_A_CD3_CD4_N_CD28 + kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_Br + ka_CD8_N_CD3_CD4_N_CD28*S_CD8_N_CD3_CD4_N_CD28_Br - kDis*S_CD8_A_CD3_CD4_N_CD28_Br
    d/dt(S_CD8_A_CD3_CD8_A_CD28_Br) <- +kb1_CD8_A_CD3_CD8_A_CD28 + kb2_CD8_A_CD3_CD8_A_CD28 + kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_Br + ka_CD8_N_CD3_CD8_A_CD28*S_CD8_N_CD3_CD8_A_CD28_Br - kDis*S_CD8_A_CD3_CD8_A_CD28_Br
    d/dt(S_CD8_A_CD3_CD8_EM_CD28_Br) <- +kb1_CD8_A_CD3_CD8_EM_CD28 + kb2_CD8_A_CD3_CD8_EM_CD28 + kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_Br + ka_CD8_N_CD3_CD8_EM_CD28*S_CD8_N_CD3_CD8_EM_CD28_Br - kDis*S_CD8_A_CD3_CD8_EM_CD28_Br
    d/dt(S_CD8_A_CD3_CD8_N_CD28_Br) <- +kb1_CD8_A_CD3_CD8_N_CD28 + kb2_CD8_A_CD3_CD8_N_CD28 + kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_Br + ka_CD8_N_CD3_CD8_N_CD28*S_CD8_N_CD3_CD8_N_CD28_Br - kDis*S_CD8_A_CD3_CD8_N_CD28_Br
    d/dt(S_CD8_A_CD3_MM_CD28_Br) <- +kb1_CD8_A_CD3_MM_CD28 + kb2_CD8_A_CD3_MM_CD28 - kkillMM_CD8*S_CD8_A_CD3_MM_CD28_Br -
       kmut_SYN*S_CD8_A_CD3_MM_CD28_Br + kact_EM*S_CD8_EM_CD3_MM_CD28_Br + ka_CD8_N_CD3_MM_CD28*S_CD8_N_CD3_MM_CD28_Br
    d/dt(S_CD8_A_CD3_MM_CD38_Br) <- +kb1_CD8_A_CD3_MM_CD38 + kb2_CD8_A_CD3_MM_CD38 - kkillMM_CD8*S_CD8_A_CD3_MM_CD38_Br -
       kmut_SYN*S_CD8_A_CD3_MM_CD38_Br + kact_EM*S_CD8_EM_CD3_MM_CD38_Br + ka_CD8_N_CD3_MM_CD38*S_CD8_N_CD3_MM_CD38_Br
    d/dt(S_CD8_A_CD3_TRGT_CD38_Br) <- +kb1_CD8_A_CD3_TRGT_CD38 + kb2_CD8_A_CD3_TRGT_CD38 - kkillTRGT_CD8*
      S_CD8_A_CD3_TRGT_CD38_Br + kact_EM*S_CD8_EM_CD3_TRGT_CD38_Br + ka_CD8_N_CD3_TRGT_CD38*S_CD8_N_CD3_TRGT_CD38_Br -
       kDis*S_CD8_A_CD3_TRGT_CD38_Br
    d/dt(S_CD8_EM_CD28_MM_CD38_Br) <- +kb1_CD8_EM_CD28_MM_CD38 + kb2_CD8_EM_CD28_MM_CD38 - kDis*S_CD8_EM_CD28_MM_CD38_Br
    d/dt(S_CD8_EM_CD28_TRGT_CD38_Br) <- +kb1_CD8_EM_CD28_TRGT_CD38 + kb2_CD8_EM_CD28_TRGT_CD38 - kDis*
      S_CD8_EM_CD28_TRGT_CD38_Br
    d/dt(S_CD8_EM_CD3_CD4_A_CD28_Br) <- +kb1_CD8_EM_CD3_CD4_A_CD28 + kb2_CD8_EM_CD3_CD4_A_CD28 - kact_EM*
      S_CD8_EM_CD3_CD4_A_CD28_Br - kDis*S_CD8_EM_CD3_CD4_A_CD28_Br
    d/dt(S_CD8_EM_CD3_CD4_EM_CD28_Br) <- +kb1_CD8_EM_CD3_CD4_EM_CD28 + kb2_CD8_EM_CD3_CD4_EM_CD28 - kact_EM*
      S_CD8_EM_CD3_CD4_EM_CD28_Br - kDis*S_CD8_EM_CD3_CD4_EM_CD28_Br
    d/dt(S_CD8_EM_CD3_CD4_N_CD28_Br) <- +kb1_CD8_EM_CD3_CD4_N_CD28 + kb2_CD8_EM_CD3_CD4_N_CD28 - kact_EM*
      S_CD8_EM_CD3_CD4_N_CD28_Br - kDis*S_CD8_EM_CD3_CD4_N_CD28_Br
    d/dt(S_CD8_EM_CD3_CD8_A_CD28_Br) <- +kb1_CD8_EM_CD3_CD8_A_CD28 + kb2_CD8_EM_CD3_CD8_A_CD28 - kact_EM*
      S_CD8_EM_CD3_CD8_A_CD28_Br - kDis*S_CD8_EM_CD3_CD8_A_CD28_Br
    d/dt(S_CD8_EM_CD3_CD8_EM_CD28_Br) <- +kb1_CD8_EM_CD3_CD8_EM_CD28 + kb2_CD8_EM_CD3_CD8_EM_CD28 - kact_EM*
      S_CD8_EM_CD3_CD8_EM_CD28_Br - kDis*S_CD8_EM_CD3_CD8_EM_CD28_Br
    d/dt(S_CD8_EM_CD3_CD8_N_CD28_Br) <- +kb1_CD8_EM_CD3_CD8_N_CD28 + kb2_CD8_EM_CD3_CD8_N_CD28 - kact_EM*
      S_CD8_EM_CD3_CD8_N_CD28_Br - kDis*S_CD8_EM_CD3_CD8_N_CD28_Br
    d/dt(S_CD8_EM_CD3_MM_CD28_Br) <- +kb1_CD8_EM_CD3_MM_CD28 + kb2_CD8_EM_CD3_MM_CD28 - kact_EM*S_CD8_EM_CD3_MM_CD28_Br -
       kDis*S_CD8_EM_CD3_MM_CD28_Br
    d/dt(S_CD8_EM_CD3_MM_CD38_Br) <- +kb1_CD8_EM_CD3_MM_CD38 + kb2_CD8_EM_CD3_MM_CD38 - kact_EM*S_CD8_EM_CD3_MM_CD38_Br -
       kDis*S_CD8_EM_CD3_MM_CD38_Br
    d/dt(S_CD8_EM_CD3_TRGT_CD38_Br) <- +kb1_CD8_EM_CD3_TRGT_CD38 + kb2_CD8_EM_CD3_TRGT_CD38 - kact_EM*
      S_CD8_EM_CD3_TRGT_CD38_Br - kDis*S_CD8_EM_CD3_TRGT_CD38_Br
    d/dt(S_CD8_N_CD28_MM_CD38_Br) <- +kb1_CD8_N_CD28_MM_CD38 + kb2_CD8_N_CD28_MM_CD38 - kDis*S_CD8_N_CD28_MM_CD38_Br
    d/dt(S_CD8_N_CD28_TRGT_CD38_Br) <- +kb1_CD8_N_CD28_TRGT_CD38 + kb2_CD8_N_CD28_TRGT_CD38 - kDis*S_CD8_N_CD28_TRGT_CD38_Br
    d/dt(S_CD8_N_CD3_CD4_A_CD28_Br) <- +kb1_CD8_N_CD3_CD4_A_CD28 + kb2_CD8_N_CD3_CD4_A_CD28 - ka_CD8_N_CD3_CD4_A_CD28*
      S_CD8_N_CD3_CD4_A_CD28_Br - kDis*S_CD8_N_CD3_CD4_A_CD28_Br
    d/dt(S_CD8_N_CD3_CD4_EM_CD28_Br) <- +kb1_CD8_N_CD3_CD4_EM_CD28 + kb2_CD8_N_CD3_CD4_EM_CD28 - ka_CD8_N_CD3_CD4_EM_CD28*
      S_CD8_N_CD3_CD4_EM_CD28_Br - kDis*S_CD8_N_CD3_CD4_EM_CD28_Br
    d/dt(S_CD8_N_CD3_CD4_N_CD28_Br) <- +kb1_CD8_N_CD3_CD4_N_CD28 + kb2_CD8_N_CD3_CD4_N_CD28 - ka_CD8_N_CD3_CD4_N_CD28*
      S_CD8_N_CD3_CD4_N_CD28_Br - kDis*S_CD8_N_CD3_CD4_N_CD28_Br
    d/dt(S_CD8_N_CD3_CD8_A_CD28_Br) <- +kb1_CD8_N_CD3_CD8_A_CD28 + kb2_CD8_N_CD3_CD8_A_CD28 - ka_CD8_N_CD3_CD8_A_CD28*
      S_CD8_N_CD3_CD8_A_CD28_Br - kDis*S_CD8_N_CD3_CD8_A_CD28_Br
    d/dt(S_CD8_N_CD3_CD8_EM_CD28_Br) <- +kb1_CD8_N_CD3_CD8_EM_CD28 + kb2_CD8_N_CD3_CD8_EM_CD28 - ka_CD8_N_CD3_CD8_EM_CD28*
      S_CD8_N_CD3_CD8_EM_CD28_Br - kDis*S_CD8_N_CD3_CD8_EM_CD28_Br
    d/dt(S_CD8_N_CD3_CD8_N_CD28_Br) <- +kb1_CD8_N_CD3_CD8_N_CD28 + kb2_CD8_N_CD3_CD8_N_CD28 - ka_CD8_N_CD3_CD8_N_CD28*
      S_CD8_N_CD3_CD8_N_CD28_Br - kDis*S_CD8_N_CD3_CD8_N_CD28_Br
    d/dt(S_CD8_N_CD3_MM_CD28_Br) <- +kb1_CD8_N_CD3_MM_CD28 + kb2_CD8_N_CD3_MM_CD28 - ka_CD8_N_CD3_MM_CD28*
      S_CD8_N_CD3_MM_CD28_Br - kDis*S_CD8_N_CD3_MM_CD28_Br
    d/dt(S_CD8_N_CD3_MM_CD38_Br) <- +kb1_CD8_N_CD3_MM_CD38 + kb2_CD8_N_CD3_MM_CD38 - ka_CD8_N_CD3_MM_CD38*
      S_CD8_N_CD3_MM_CD38_Br - kDis*S_CD8_N_CD3_MM_CD38_Br
    d/dt(S_CD8_N_CD3_TRGT_CD38_Br) <- +kb1_CD8_N_CD3_TRGT_CD38 + kb2_CD8_N_CD3_TRGT_CD38 - ka_CD8_N_CD3_TRGT_CD38*
      S_CD8_N_CD3_TRGT_CD38_Br - kDis*S_CD8_N_CD3_TRGT_CD38_Br
    d/dt(S_MM_CD28_MM_CD38_Br) <- +kb1_MM_CD28_MM_CD38 + kb2_MM_CD28_MM_CD38 - kDis*S_MM_CD28_MM_CD38_Br
    d/dt(S_MM_CD28_TRGT_CD38_Br) <- +kb1_MM_CD28_TRGT_CD38 + kb2_MM_CD28_TRGT_CD38 - kDis*S_MM_CD28_TRGT_CD38_Br
    d/dt(S_CD4_A_CD28_MM_CD38_MUT_Br) <- +kmut_SYN*S_CD4_A_CD28_MM_CD38_Br - kDis*S_CD4_A_CD28_MM_CD38_MUT_Br
    d/dt(S_CD4_A_CD3_MM_CD28_MUT_Br) <- +kmut_SYN*S_CD4_A_CD3_MM_CD28_Br - kDis*S_CD4_A_CD3_MM_CD28_MUT_Br
    d/dt(S_CD4_A_CD3_MM_CD38_MUT_Br) <- +kmut_SYN*S_CD4_A_CD3_MM_CD38_Br - kDis*S_CD4_A_CD3_MM_CD38_MUT_Br
    d/dt(S_CD8_A_CD28_MM_CD38_MUT_Br) <- +kmut_SYN*S_CD8_A_CD28_MM_CD38_Br - kDis*S_CD8_A_CD28_MM_CD38_MUT_Br
    d/dt(S_CD8_A_CD3_MM_CD28_MUT_Br) <- +kmut_SYN*S_CD8_A_CD3_MM_CD28_Br - kDis*S_CD8_A_CD3_MM_CD28_MUT_Br
    d/dt(S_CD8_A_CD3_MM_CD38_MUT_Br) <- +kmut_SYN*S_CD8_A_CD3_MM_CD38_Br - kDis*S_CD8_A_CD3_MM_CD38_MUT_Br
    d/dt(IFNg) <- +kprod_CD8_A_IFNg*CD8_A + kprod_CD4_A_IFNg*CD4_A - kdeg_IFNg*IFNg + (kprod_CD4_A_IFNg +
       kprod_TRGT_IFNg)*S_CD4_A_CD28_MM_CD38 + (kprod_CD4_A_IFNg + kprod_TRGT_IFNg)*S_CD4_A_CD28_TRGT_CD38 +
       (kprod_CD4_A_IFNg + kprod_CD4_A_IFNg)*S_CD4_A_CD3_CD4_A_CD28 + (kprod_CD4_A_IFNg + kprod_CD4_EM_IFNg)*
      S_CD4_A_CD3_CD4_EM_CD28 + (kprod_CD4_A_IFNg + kprod_CD4_N_IFNg)*S_CD4_A_CD3_CD4_N_CD28 + (kprod_CD4_A_IFNg +
       kprod_CD8_A_IFNg)*S_CD4_A_CD3_CD8_A_CD28 + (kprod_CD4_A_IFNg + kprod_CD8_EM_IFNg)*S_CD4_A_CD3_CD8_EM_CD28 +
       (kprod_CD4_A_IFNg + kprod_CD8_N_IFNg)*S_CD4_A_CD3_CD8_N_CD28 + (kprod_CD4_A_IFNg + kprod_TRGT_IFNg)*
      S_CD4_A_CD3_MM_CD28 + (kprod_CD4_A_IFNg + kprod_TRGT_IFNg)*S_CD4_A_CD3_MM_CD38 + (kprod_CD4_A_IFNg +
       kprod_TRGT_IFNg)*S_CD4_A_CD3_TRGT_CD38 + (kprod_CD4_EM_IFNg + kprod_TRGT_IFNg)*S_CD4_EM_CD28_MM_CD38 +
       (kprod_CD4_EM_IFNg + kprod_TRGT_IFNg)*S_CD4_EM_CD28_TRGT_CD38 + (kprod_CD4_EM_IFNg + kprod_CD4_A_IFNg)*
      S_CD4_EM_CD3_CD4_A_CD28 + (kprod_CD4_EM_IFNg + kprod_CD4_EM_IFNg)*S_CD4_EM_CD3_CD4_EM_CD28 + (kprod_CD4_EM_IFNg +
       kprod_CD4_N_IFNg)*S_CD4_EM_CD3_CD4_N_CD28 + (kprod_CD4_EM_IFNg + kprod_CD8_A_IFNg)*S_CD4_EM_CD3_CD8_A_CD28 +
       (kprod_CD4_EM_IFNg + kprod_CD8_EM_IFNg)*S_CD4_EM_CD3_CD8_EM_CD28 + (kprod_CD4_EM_IFNg + kprod_CD8_N_IFNg)*
      S_CD4_EM_CD3_CD8_N_CD28 + (kprod_CD4_EM_IFNg + kprod_TRGT_IFNg)*S_CD4_EM_CD3_MM_CD28 + (kprod_CD4_EM_IFNg +
       kprod_TRGT_IFNg)*S_CD4_EM_CD3_MM_CD38 + (kprod_CD4_EM_IFNg + kprod_TRGT_IFNg)*S_CD4_EM_CD3_TRGT_CD38 +
       (kprod_CD4_N_IFNg + kprod_TRGT_IFNg)*S_CD4_N_CD28_MM_CD38 + (kprod_CD4_N_IFNg + kprod_TRGT_IFNg)*
      S_CD4_N_CD28_TRGT_CD38 + (kprod_CD4_N_IFNg + kprod_CD4_A_IFNg)*S_CD4_N_CD3_CD4_A_CD28 + (kprod_CD4_N_IFNg +
       kprod_CD4_EM_IFNg)*S_CD4_N_CD3_CD4_EM_CD28 + (kprod_CD4_N_IFNg + kprod_CD4_N_IFNg)*S_CD4_N_CD3_CD4_N_CD28 +
       (kprod_CD4_N_IFNg + kprod_CD8_A_IFNg)*S_CD4_N_CD3_CD8_A_CD28 + (kprod_CD4_N_IFNg + kprod_CD8_EM_IFNg)*
      S_CD4_N_CD3_CD8_EM_CD28 + (kprod_CD4_N_IFNg + kprod_CD8_N_IFNg)*S_CD4_N_CD3_CD8_N_CD28 + (kprod_CD4_N_IFNg +
       kprod_TRGT_IFNg)*S_CD4_N_CD3_MM_CD28 + (kprod_CD4_N_IFNg + kprod_TRGT_IFNg)*S_CD4_N_CD3_MM_CD38 +
       (kprod_CD4_N_IFNg + kprod_TRGT_IFNg)*S_CD4_N_CD3_TRGT_CD38 + (kprod_CD8_A_IFNg + kprod_TRGT_IFNg)*
      S_CD8_A_CD28_MM_CD38 + (kprod_CD8_A_IFNg + kprod_TRGT_IFNg)*S_CD8_A_CD28_TRGT_CD38 + (kprod_CD8_A_IFNg +
       kprod_CD4_A_IFNg)*S_CD8_A_CD3_CD4_A_CD28 + (kprod_CD8_A_IFNg + kprod_CD4_EM_IFNg)*S_CD8_A_CD3_CD4_EM_CD28 +
       (kprod_CD8_A_IFNg + kprod_CD4_N_IFNg)*S_CD8_A_CD3_CD4_N_CD28 + (kprod_CD8_A_IFNg + kprod_CD8_A_IFNg)*
      S_CD8_A_CD3_CD8_A_CD28 + (kprod_CD8_A_IFNg + kprod_CD8_EM_IFNg)*S_CD8_A_CD3_CD8_EM_CD28 + (kprod_CD8_A_IFNg +
       kprod_CD8_N_IFNg)*S_CD8_A_CD3_CD8_N_CD28 + (kprod_CD8_A_IFNg + kprod_TRGT_IFNg)*S_CD8_A_CD3_MM_CD28 +
       (kprod_CD8_A_IFNg + kprod_TRGT_IFNg)*S_CD8_A_CD3_MM_CD38 + (kprod_CD8_A_IFNg + kprod_TRGT_IFNg)*
      S_CD8_A_CD3_TRGT_CD38 + (kprod_CD8_EM_IFNg + kprod_TRGT_IFNg)*S_CD8_EM_CD28_MM_CD38 + (kprod_CD8_EM_IFNg +
       kprod_TRGT_IFNg)*S_CD8_EM_CD28_TRGT_CD38 + (kprod_CD8_EM_IFNg + kprod_CD4_A_IFNg)*S_CD8_EM_CD3_CD4_A_CD28 +
       (kprod_CD8_EM_IFNg + kprod_CD4_EM_IFNg)*S_CD8_EM_CD3_CD4_EM_CD28 + (kprod_CD8_EM_IFNg + kprod_CD4_N_IFNg)*
      S_CD8_EM_CD3_CD4_N_CD28 + (kprod_CD8_EM_IFNg + kprod_CD8_A_IFNg)*S_CD8_EM_CD3_CD8_A_CD28 + (kprod_CD8_EM_IFNg +
       kprod_CD8_EM_IFNg)*S_CD8_EM_CD3_CD8_EM_CD28 + (kprod_CD8_EM_IFNg + kprod_CD8_N_IFNg)*S_CD8_EM_CD3_CD8_N_CD28 +
       (kprod_CD8_EM_IFNg + kprod_TRGT_IFNg)*S_CD8_EM_CD3_MM_CD28 + (kprod_CD8_EM_IFNg + kprod_TRGT_IFNg)*
      S_CD8_EM_CD3_MM_CD38 + (kprod_CD8_EM_IFNg + kprod_TRGT_IFNg)*S_CD8_EM_CD3_TRGT_CD38 + (kprod_CD8_N_IFNg +
       kprod_TRGT_IFNg)*S_CD8_N_CD28_MM_CD38 + (kprod_CD8_N_IFNg + kprod_TRGT_IFNg)*S_CD8_N_CD28_TRGT_CD38 +
       (kprod_CD8_N_IFNg + kprod_CD4_A_IFNg)*S_CD8_N_CD3_CD4_A_CD28 + (kprod_CD8_N_IFNg + kprod_CD4_EM_IFNg)*
      S_CD8_N_CD3_CD4_EM_CD28 + (kprod_CD8_N_IFNg + kprod_CD4_N_IFNg)*S_CD8_N_CD3_CD4_N_CD28 + (kprod_CD8_N_IFNg +
       kprod_CD8_A_IFNg)*S_CD8_N_CD3_CD8_A_CD28 + (kprod_CD8_N_IFNg + kprod_CD8_EM_IFNg)*S_CD8_N_CD3_CD8_EM_CD28 +
       (kprod_CD8_N_IFNg + kprod_CD8_N_IFNg)*S_CD8_N_CD3_CD8_N_CD28 + (kprod_CD8_N_IFNg + kprod_TRGT_IFNg)*
      S_CD8_N_CD3_MM_CD28 + (kprod_CD8_N_IFNg + kprod_TRGT_IFNg)*S_CD8_N_CD3_MM_CD38 + (kprod_CD8_N_IFNg +
       kprod_TRGT_IFNg)*S_CD8_N_CD3_TRGT_CD38
    d/dt(TNFa) <- +kprod_CD8_A_TNFa*CD8_A + kprod_CD4_A_TNFa*CD4_A - kdeg_TNFa*TNFa + (kprod_CD4_A_TNFa +
       kprod_TRGT_TNFa)*S_CD4_A_CD28_MM_CD38 + (kprod_CD4_A_TNFa + kprod_TRGT_TNFa)*S_CD4_A_CD28_TRGT_CD38 +
       (kprod_CD4_A_TNFa + kprod_CD4_A_TNFa)*S_CD4_A_CD3_CD4_A_CD28 + (kprod_CD4_A_TNFa + kprod_CD4_EM_TNFa)*
      S_CD4_A_CD3_CD4_EM_CD28 + (kprod_CD4_A_TNFa + kprod_CD4_N_TNFa)*S_CD4_A_CD3_CD4_N_CD28 + (kprod_CD4_A_TNFa +
       kprod_CD8_A_TNFa)*S_CD4_A_CD3_CD8_A_CD28 + (kprod_CD4_A_TNFa + kprod_CD8_EM_TNFa)*S_CD4_A_CD3_CD8_EM_CD28 +
       (kprod_CD4_A_TNFa + kprod_CD8_N_TNFa)*S_CD4_A_CD3_CD8_N_CD28 + (kprod_CD4_A_TNFa + kprod_TRGT_TNFa)*
      S_CD4_A_CD3_MM_CD28 + (kprod_CD4_A_TNFa + kprod_TRGT_TNFa)*S_CD4_A_CD3_MM_CD38 + (kprod_CD4_A_TNFa +
       kprod_TRGT_TNFa)*S_CD4_A_CD3_TRGT_CD38 + (kprod_CD4_EM_TNFa + kprod_TRGT_TNFa)*S_CD4_EM_CD28_MM_CD38 +
       (kprod_CD4_EM_TNFa + kprod_TRGT_TNFa)*S_CD4_EM_CD28_TRGT_CD38 + (kprod_CD4_EM_TNFa + kprod_CD4_A_TNFa)*
      S_CD4_EM_CD3_CD4_A_CD28 + (kprod_CD4_EM_TNFa + kprod_CD4_EM_TNFa)*S_CD4_EM_CD3_CD4_EM_CD28 + (kprod_CD4_EM_TNFa +
       kprod_CD4_N_TNFa)*S_CD4_EM_CD3_CD4_N_CD28 + (kprod_CD4_EM_TNFa + kprod_CD8_A_TNFa)*S_CD4_EM_CD3_CD8_A_CD28 +
       (kprod_CD4_EM_TNFa + kprod_CD8_EM_TNFa)*S_CD4_EM_CD3_CD8_EM_CD28 + (kprod_CD4_EM_TNFa + kprod_CD8_N_TNFa)*
      S_CD4_EM_CD3_CD8_N_CD28 + (kprod_CD4_EM_TNFa + kprod_TRGT_TNFa)*S_CD4_EM_CD3_MM_CD28 + (kprod_CD4_EM_TNFa +
       kprod_TRGT_TNFa)*S_CD4_EM_CD3_MM_CD38 + (kprod_CD4_EM_TNFa + kprod_TRGT_TNFa)*S_CD4_EM_CD3_TRGT_CD38 +
       (kprod_CD4_N_TNFa + kprod_TRGT_TNFa)*S_CD4_N_CD28_MM_CD38 + (kprod_CD4_N_TNFa + kprod_TRGT_TNFa)*
      S_CD4_N_CD28_TRGT_CD38 + (kprod_CD4_N_TNFa + kprod_CD4_A_TNFa)*S_CD4_N_CD3_CD4_A_CD28 + (kprod_CD4_N_TNFa +
       kprod_CD4_EM_TNFa)*S_CD4_N_CD3_CD4_EM_CD28 + (kprod_CD4_N_TNFa + kprod_CD4_N_TNFa)*S_CD4_N_CD3_CD4_N_CD28 +
       (kprod_CD4_N_TNFa + kprod_CD8_A_TNFa)*S_CD4_N_CD3_CD8_A_CD28 + (kprod_CD4_N_TNFa + kprod_CD8_EM_TNFa)*
      S_CD4_N_CD3_CD8_EM_CD28 + (kprod_CD4_N_TNFa + kprod_CD8_N_TNFa)*S_CD4_N_CD3_CD8_N_CD28 + (kprod_CD4_N_TNFa +
       kprod_TRGT_TNFa)*S_CD4_N_CD3_MM_CD28 + (kprod_CD4_N_TNFa + kprod_TRGT_TNFa)*S_CD4_N_CD3_MM_CD38 +
       (kprod_CD4_N_TNFa + kprod_TRGT_TNFa)*S_CD4_N_CD3_TRGT_CD38 + (kprod_CD8_A_TNFa + kprod_TRGT_TNFa)*
      S_CD8_A_CD28_MM_CD38 + (kprod_CD8_A_TNFa + kprod_TRGT_TNFa)*S_CD8_A_CD28_TRGT_CD38 + (kprod_CD8_A_TNFa +
       kprod_CD4_A_TNFa)*S_CD8_A_CD3_CD4_A_CD28 + (kprod_CD8_A_TNFa + kprod_CD4_EM_TNFa)*S_CD8_A_CD3_CD4_EM_CD28 +
       (kprod_CD8_A_TNFa + kprod_CD4_N_TNFa)*S_CD8_A_CD3_CD4_N_CD28 + (kprod_CD8_A_TNFa + kprod_CD8_A_TNFa)*
      S_CD8_A_CD3_CD8_A_CD28 + (kprod_CD8_A_TNFa + kprod_CD8_EM_TNFa)*S_CD8_A_CD3_CD8_EM_CD28 + (kprod_CD8_A_TNFa +
       kprod_CD8_N_TNFa)*S_CD8_A_CD3_CD8_N_CD28 + (kprod_CD8_A_TNFa + kprod_TRGT_TNFa)*S_CD8_A_CD3_MM_CD28 +
       (kprod_CD8_A_TNFa + kprod_TRGT_TNFa)*S_CD8_A_CD3_MM_CD38 + (kprod_CD8_A_TNFa + kprod_TRGT_TNFa)*
      S_CD8_A_CD3_TRGT_CD38 + (kprod_CD8_EM_TNFa + kprod_TRGT_TNFa)*S_CD8_EM_CD28_MM_CD38 + (kprod_CD8_EM_TNFa +
       kprod_TRGT_TNFa)*S_CD8_EM_CD28_TRGT_CD38 + (kprod_CD8_EM_TNFa + kprod_CD4_A_TNFa)*S_CD8_EM_CD3_CD4_A_CD28 +
       (kprod_CD8_EM_TNFa + kprod_CD4_EM_TNFa)*S_CD8_EM_CD3_CD4_EM_CD28 + (kprod_CD8_EM_TNFa + kprod_CD4_N_TNFa)*
      S_CD8_EM_CD3_CD4_N_CD28 + (kprod_CD8_EM_TNFa + kprod_CD8_A_TNFa)*S_CD8_EM_CD3_CD8_A_CD28 + (kprod_CD8_EM_TNFa +
       kprod_CD8_EM_TNFa)*S_CD8_EM_CD3_CD8_EM_CD28 + (kprod_CD8_EM_TNFa + kprod_CD8_N_TNFa)*S_CD8_EM_CD3_CD8_N_CD28 +
       (kprod_CD8_EM_TNFa + kprod_TRGT_TNFa)*S_CD8_EM_CD3_MM_CD28 + (kprod_CD8_EM_TNFa + kprod_TRGT_TNFa)*
      S_CD8_EM_CD3_MM_CD38 + (kprod_CD8_EM_TNFa + kprod_TRGT_TNFa)*S_CD8_EM_CD3_TRGT_CD38 + (kprod_CD8_N_TNFa +
       kprod_TRGT_TNFa)*S_CD8_N_CD28_MM_CD38 + (kprod_CD8_N_TNFa + kprod_TRGT_TNFa)*S_CD8_N_CD28_TRGT_CD38 +
       (kprod_CD8_N_TNFa + kprod_CD4_A_TNFa)*S_CD8_N_CD3_CD4_A_CD28 + (kprod_CD8_N_TNFa + kprod_CD4_EM_TNFa)*
      S_CD8_N_CD3_CD4_EM_CD28 + (kprod_CD8_N_TNFa + kprod_CD4_N_TNFa)*S_CD8_N_CD3_CD4_N_CD28 + (kprod_CD8_N_TNFa +
       kprod_CD8_A_TNFa)*S_CD8_N_CD3_CD8_A_CD28 + (kprod_CD8_N_TNFa + kprod_CD8_EM_TNFa)*S_CD8_N_CD3_CD8_EM_CD28 +
       (kprod_CD8_N_TNFa + kprod_CD8_N_TNFa)*S_CD8_N_CD3_CD8_N_CD28 + (kprod_CD8_N_TNFa + kprod_TRGT_TNFa)*
      S_CD8_N_CD3_MM_CD28 + (kprod_CD8_N_TNFa + kprod_TRGT_TNFa)*S_CD8_N_CD3_MM_CD38 + (kprod_CD8_N_TNFa +
       kprod_TRGT_TNFa)*S_CD8_N_CD3_TRGT_CD38
    d/dt(IL6) <- +kprod_CD8_A_IL6*CD8_A + kprod_CD4_A_IL6*CD4_A - kdeg_IL6*IL6 + (kprod_CD4_A_IL6 + kprod_TRGT_IL6)*
      S_CD4_A_CD28_MM_CD38 + (kprod_CD4_A_IL6 + kprod_TRGT_IL6)*S_CD4_A_CD28_TRGT_CD38 + (kprod_CD4_A_IL6 +
       kprod_CD4_A_IL6)*S_CD4_A_CD3_CD4_A_CD28 + (kprod_CD4_A_IL6 + kprod_CD4_EM_IL6)*S_CD4_A_CD3_CD4_EM_CD28 +
       (kprod_CD4_A_IL6 + kprod_CD4_N_IL6)*S_CD4_A_CD3_CD4_N_CD28 + (kprod_CD4_A_IL6 + kprod_CD8_A_IL6)*
      S_CD4_A_CD3_CD8_A_CD28 + (kprod_CD4_A_IL6 + kprod_CD8_EM_IL6)*S_CD4_A_CD3_CD8_EM_CD28 + (kprod_CD4_A_IL6 +
       kprod_CD8_N_IL6)*S_CD4_A_CD3_CD8_N_CD28 + (kprod_CD4_A_IL6 + kprod_TRGT_IL6)*S_CD4_A_CD3_MM_CD28 +
       (kprod_CD4_A_IL6 + kprod_TRGT_IL6)*S_CD4_A_CD3_MM_CD38 + (kprod_CD4_A_IL6 + kprod_TRGT_IL6)*S_CD4_A_CD3_TRGT_CD38 +
       (kprod_CD4_EM_IL6 + kprod_TRGT_IL6)*S_CD4_EM_CD28_MM_CD38 + (kprod_CD4_EM_IL6 + kprod_TRGT_IL6)*
      S_CD4_EM_CD28_TRGT_CD38 + (kprod_CD4_EM_IL6 + kprod_CD4_A_IL6)*S_CD4_EM_CD3_CD4_A_CD28 + (kprod_CD4_EM_IL6 +
       kprod_CD4_EM_IL6)*S_CD4_EM_CD3_CD4_EM_CD28 + (kprod_CD4_EM_IL6 + kprod_CD4_N_IL6)*S_CD4_EM_CD3_CD4_N_CD28 +
       (kprod_CD4_EM_IL6 + kprod_CD8_A_IL6)*S_CD4_EM_CD3_CD8_A_CD28 + (kprod_CD4_EM_IL6 + kprod_CD8_EM_IL6)*
      S_CD4_EM_CD3_CD8_EM_CD28 + (kprod_CD4_EM_IL6 + kprod_CD8_N_IL6)*S_CD4_EM_CD3_CD8_N_CD28 + (kprod_CD4_EM_IL6 +
       kprod_TRGT_IL6)*S_CD4_EM_CD3_MM_CD28 + (kprod_CD4_EM_IL6 + kprod_TRGT_IL6)*S_CD4_EM_CD3_MM_CD38 +
       (kprod_CD4_EM_IL6 + kprod_TRGT_IL6)*S_CD4_EM_CD3_TRGT_CD38 + (kprod_CD4_N_IL6 + kprod_TRGT_IL6)*
      S_CD4_N_CD28_MM_CD38 + (kprod_CD4_N_IL6 + kprod_TRGT_IL6)*S_CD4_N_CD28_TRGT_CD38 + (kprod_CD4_N_IL6 +
       kprod_CD4_A_IL6)*S_CD4_N_CD3_CD4_A_CD28 + (kprod_CD4_N_IL6 + kprod_CD4_EM_IL6)*S_CD4_N_CD3_CD4_EM_CD28 +
       (kprod_CD4_N_IL6 + kprod_CD4_N_IL6)*S_CD4_N_CD3_CD4_N_CD28 + (kprod_CD4_N_IL6 + kprod_CD8_A_IL6)*
      S_CD4_N_CD3_CD8_A_CD28 + (kprod_CD4_N_IL6 + kprod_CD8_EM_IL6)*S_CD4_N_CD3_CD8_EM_CD28 + (kprod_CD4_N_IL6 +
       kprod_CD8_N_IL6)*S_CD4_N_CD3_CD8_N_CD28 + (kprod_CD4_N_IL6 + kprod_TRGT_IL6)*S_CD4_N_CD3_MM_CD28 +
       (kprod_CD4_N_IL6 + kprod_TRGT_IL6)*S_CD4_N_CD3_MM_CD38 + (kprod_CD4_N_IL6 + kprod_TRGT_IL6)*S_CD4_N_CD3_TRGT_CD38 +
       (kprod_CD8_A_IL6 + kprod_TRGT_IL6)*S_CD8_A_CD28_MM_CD38 + (kprod_CD8_A_IL6 + kprod_TRGT_IL6)*S_CD8_A_CD28_TRGT_CD38 +
       (kprod_CD8_A_IL6 + kprod_CD4_A_IL6)*S_CD8_A_CD3_CD4_A_CD28 + (kprod_CD8_A_IL6 + kprod_CD4_EM_IL6)*
      S_CD8_A_CD3_CD4_EM_CD28 + (kprod_CD8_A_IL6 + kprod_CD4_N_IL6)*S_CD8_A_CD3_CD4_N_CD28 + (kprod_CD8_A_IL6 +
       kprod_CD8_A_IL6)*S_CD8_A_CD3_CD8_A_CD28 + (kprod_CD8_A_IL6 + kprod_CD8_EM_IL6)*S_CD8_A_CD3_CD8_EM_CD28 +
       (kprod_CD8_A_IL6 + kprod_CD8_N_IL6)*S_CD8_A_CD3_CD8_N_CD28 + (kprod_CD8_A_IL6 + kprod_TRGT_IL6)*
      S_CD8_A_CD3_MM_CD28 + (kprod_CD8_A_IL6 + kprod_TRGT_IL6)*S_CD8_A_CD3_MM_CD38 + (kprod_CD8_A_IL6 +
       kprod_TRGT_IL6)*S_CD8_A_CD3_TRGT_CD38 + (kprod_CD8_EM_IL6 + kprod_TRGT_IL6)*S_CD8_EM_CD28_MM_CD38 +
       (kprod_CD8_EM_IL6 + kprod_TRGT_IL6)*S_CD8_EM_CD28_TRGT_CD38 + (kprod_CD8_EM_IL6 + kprod_CD4_A_IL6)*
      S_CD8_EM_CD3_CD4_A_CD28 + (kprod_CD8_EM_IL6 + kprod_CD4_EM_IL6)*S_CD8_EM_CD3_CD4_EM_CD28 + (kprod_CD8_EM_IL6 +
       kprod_CD4_N_IL6)*S_CD8_EM_CD3_CD4_N_CD28 + (kprod_CD8_EM_IL6 + kprod_CD8_A_IL6)*S_CD8_EM_CD3_CD8_A_CD28 +
       (kprod_CD8_EM_IL6 + kprod_CD8_EM_IL6)*S_CD8_EM_CD3_CD8_EM_CD28 + (kprod_CD8_EM_IL6 + kprod_CD8_N_IL6)*
      S_CD8_EM_CD3_CD8_N_CD28 + (kprod_CD8_EM_IL6 + kprod_TRGT_IL6)*S_CD8_EM_CD3_MM_CD28 + (kprod_CD8_EM_IL6 +
       kprod_TRGT_IL6)*S_CD8_EM_CD3_MM_CD38 + (kprod_CD8_EM_IL6 + kprod_TRGT_IL6)*S_CD8_EM_CD3_TRGT_CD38 +
       (kprod_CD8_N_IL6 + kprod_TRGT_IL6)*S_CD8_N_CD28_MM_CD38 + (kprod_CD8_N_IL6 + kprod_TRGT_IL6)*S_CD8_N_CD28_TRGT_CD38 +
       (kprod_CD8_N_IL6 + kprod_CD4_A_IL6)*S_CD8_N_CD3_CD4_A_CD28 + (kprod_CD8_N_IL6 + kprod_CD4_EM_IL6)*
      S_CD8_N_CD3_CD4_EM_CD28 + (kprod_CD8_N_IL6 + kprod_CD4_N_IL6)*S_CD8_N_CD3_CD4_N_CD28 + (kprod_CD8_N_IL6 +
       kprod_CD8_A_IL6)*S_CD8_N_CD3_CD8_A_CD28 + (kprod_CD8_N_IL6 + kprod_CD8_EM_IL6)*S_CD8_N_CD3_CD8_EM_CD28 +
       (kprod_CD8_N_IL6 + kprod_CD8_N_IL6)*S_CD8_N_CD3_CD8_N_CD28 + (kprod_CD8_N_IL6 + kprod_TRGT_IL6)*
      S_CD8_N_CD3_MM_CD28 + (kprod_CD8_N_IL6 + kprod_TRGT_IL6)*S_CD8_N_CD3_MM_CD38 + (kprod_CD8_N_IL6 +
       kprod_TRGT_IL6)*S_CD8_N_CD3_TRGT_CD38
    d/dt(IL10) <- -kdeg_IL10*IL10 + (kprod_CD4_A_IL10 + kprod_TRGT_IL10)*S_CD4_A_CD28_MM_CD38 + (kprod_CD4_A_IL10 +
       kprod_TRGT_IL10)*S_CD4_A_CD28_TRGT_CD38 + (kprod_CD4_A_IL10 + kprod_CD4_A_IL10)*S_CD4_A_CD3_CD4_A_CD28 +
       (kprod_CD4_A_IL10 + kprod_CD4_EM_IL10)*S_CD4_A_CD3_CD4_EM_CD28 + (kprod_CD4_A_IL10 + kprod_CD4_N_IL10)*
      S_CD4_A_CD3_CD4_N_CD28 + (kprod_CD4_A_IL10 + kprod_CD8_A_IL10)*S_CD4_A_CD3_CD8_A_CD28 + (kprod_CD4_A_IL10 +
       kprod_CD8_EM_IL10)*S_CD4_A_CD3_CD8_EM_CD28 + (kprod_CD4_A_IL10 + kprod_CD8_N_IL10)*S_CD4_A_CD3_CD8_N_CD28 +
       (kprod_CD4_A_IL10 + kprod_TRGT_IL10)*S_CD4_A_CD3_MM_CD28 + (kprod_CD4_A_IL10 + kprod_TRGT_IL10)*
      S_CD4_A_CD3_MM_CD38 + (kprod_CD4_A_IL10 + kprod_TRGT_IL10)*S_CD4_A_CD3_TRGT_CD38 + (kprod_CD4_EM_IL10 +
       kprod_TRGT_IL10)*S_CD4_EM_CD28_MM_CD38 + (kprod_CD4_EM_IL10 + kprod_TRGT_IL10)*S_CD4_EM_CD28_TRGT_CD38 +
       (kprod_CD4_EM_IL10 + kprod_CD4_A_IL10)*S_CD4_EM_CD3_CD4_A_CD28 + (kprod_CD4_EM_IL10 + kprod_CD4_EM_IL10)*
      S_CD4_EM_CD3_CD4_EM_CD28 + (kprod_CD4_EM_IL10 + kprod_CD4_N_IL10)*S_CD4_EM_CD3_CD4_N_CD28 + (kprod_CD4_EM_IL10 +
       kprod_CD8_A_IL10)*S_CD4_EM_CD3_CD8_A_CD28 + (kprod_CD4_EM_IL10 + kprod_CD8_EM_IL10)*S_CD4_EM_CD3_CD8_EM_CD28 +
       (kprod_CD4_EM_IL10 + kprod_CD8_N_IL10)*S_CD4_EM_CD3_CD8_N_CD28 + (kprod_CD4_EM_IL10 + kprod_TRGT_IL10)*
      S_CD4_EM_CD3_MM_CD28 + (kprod_CD4_EM_IL10 + kprod_TRGT_IL10)*S_CD4_EM_CD3_MM_CD38 + (kprod_CD4_EM_IL10 +
       kprod_TRGT_IL10)*S_CD4_EM_CD3_TRGT_CD38 + (kprod_CD4_N_IL10 + kprod_TRGT_IL10)*S_CD4_N_CD28_MM_CD38 +
       (kprod_CD4_N_IL10 + kprod_TRGT_IL10)*S_CD4_N_CD28_TRGT_CD38 + (kprod_CD4_N_IL10 + kprod_CD4_A_IL10)*
      S_CD4_N_CD3_CD4_A_CD28 + (kprod_CD4_N_IL10 + kprod_CD4_EM_IL10)*S_CD4_N_CD3_CD4_EM_CD28 + (kprod_CD4_N_IL10 +
       kprod_CD4_N_IL10)*S_CD4_N_CD3_CD4_N_CD28 + (kprod_CD4_N_IL10 + kprod_CD8_A_IL10)*S_CD4_N_CD3_CD8_A_CD28 +
       (kprod_CD4_N_IL10 + kprod_CD8_EM_IL10)*S_CD4_N_CD3_CD8_EM_CD28 + (kprod_CD4_N_IL10 + kprod_CD8_N_IL10)*
      S_CD4_N_CD3_CD8_N_CD28 + (kprod_CD4_N_IL10 + kprod_TRGT_IL10)*S_CD4_N_CD3_MM_CD28 + (kprod_CD4_N_IL10 +
       kprod_TRGT_IL10)*S_CD4_N_CD3_MM_CD38 + (kprod_CD4_N_IL10 + kprod_TRGT_IL10)*S_CD4_N_CD3_TRGT_CD38 +
       (kprod_CD8_A_IL10 + kprod_TRGT_IL10)*S_CD8_A_CD28_MM_CD38 + (kprod_CD8_A_IL10 + kprod_TRGT_IL10)*
      S_CD8_A_CD28_TRGT_CD38 + (kprod_CD8_A_IL10 + kprod_CD4_A_IL10)*S_CD8_A_CD3_CD4_A_CD28 + (kprod_CD8_A_IL10 +
       kprod_CD4_EM_IL10)*S_CD8_A_CD3_CD4_EM_CD28 + (kprod_CD8_A_IL10 + kprod_CD4_N_IL10)*S_CD8_A_CD3_CD4_N_CD28 +
       (kprod_CD8_A_IL10 + kprod_CD8_A_IL10)*S_CD8_A_CD3_CD8_A_CD28 + (kprod_CD8_A_IL10 + kprod_CD8_EM_IL10)*
      S_CD8_A_CD3_CD8_EM_CD28 + (kprod_CD8_A_IL10 + kprod_CD8_N_IL10)*S_CD8_A_CD3_CD8_N_CD28 + (kprod_CD8_A_IL10 +
       kprod_TRGT_IL10)*S_CD8_A_CD3_MM_CD28 + (kprod_CD8_A_IL10 + kprod_TRGT_IL10)*S_CD8_A_CD3_MM_CD38 +
       (kprod_CD8_A_IL10 + kprod_TRGT_IL10)*S_CD8_A_CD3_TRGT_CD38 + (kprod_CD8_EM_IL10 + kprod_TRGT_IL10)*
      S_CD8_EM_CD28_MM_CD38 + (kprod_CD8_EM_IL10 + kprod_TRGT_IL10)*S_CD8_EM_CD28_TRGT_CD38 + (kprod_CD8_EM_IL10 +
       kprod_CD4_A_IL10)*S_CD8_EM_CD3_CD4_A_CD28 + (kprod_CD8_EM_IL10 + kprod_CD4_EM_IL10)*S_CD8_EM_CD3_CD4_EM_CD28 +
       (kprod_CD8_EM_IL10 + kprod_CD4_N_IL10)*S_CD8_EM_CD3_CD4_N_CD28 + (kprod_CD8_EM_IL10 + kprod_CD8_A_IL10)*
      S_CD8_EM_CD3_CD8_A_CD28 + (kprod_CD8_EM_IL10 + kprod_CD8_EM_IL10)*S_CD8_EM_CD3_CD8_EM_CD28 + (kprod_CD8_EM_IL10 +
       kprod_CD8_N_IL10)*S_CD8_EM_CD3_CD8_N_CD28 + (kprod_CD8_EM_IL10 + kprod_TRGT_IL10)*S_CD8_EM_CD3_MM_CD28 +
       (kprod_CD8_EM_IL10 + kprod_TRGT_IL10)*S_CD8_EM_CD3_MM_CD38 + (kprod_CD8_EM_IL10 + kprod_TRGT_IL10)*
      S_CD8_EM_CD3_TRGT_CD38 + (kprod_CD8_N_IL10 + kprod_TRGT_IL10)*S_CD8_N_CD28_MM_CD38 + (kprod_CD8_N_IL10 +
       kprod_TRGT_IL10)*S_CD8_N_CD28_TRGT_CD38 + (kprod_CD8_N_IL10 + kprod_CD4_A_IL10)*S_CD8_N_CD3_CD4_A_CD28 +
       (kprod_CD8_N_IL10 + kprod_CD4_EM_IL10)*S_CD8_N_CD3_CD4_EM_CD28 + (kprod_CD8_N_IL10 + kprod_CD4_N_IL10)*
      S_CD8_N_CD3_CD4_N_CD28 + (kprod_CD8_N_IL10 + kprod_CD8_A_IL10)*S_CD8_N_CD3_CD8_A_CD28 + (kprod_CD8_N_IL10 +
       kprod_CD8_EM_IL10)*S_CD8_N_CD3_CD8_EM_CD28 + (kprod_CD8_N_IL10 + kprod_CD8_N_IL10)*S_CD8_N_CD3_CD8_N_CD28 +
       (kprod_CD8_N_IL10 + kprod_TRGT_IL10)*S_CD8_N_CD3_MM_CD28 + (kprod_CD8_N_IL10 + kprod_TRGT_IL10)*
      S_CD8_N_CD3_MM_CD38 + (kprod_CD8_N_IL10 + kprod_TRGT_IL10)*S_CD8_N_CD3_TRGT_CD38
    d/dt(TcellSyn) <- +ksyn_CD8_N + ksyn_CD8_EM + kpr_CD8*CD8_A + ksyn_CD4_N + ksyn_CD4_EM + kpr_CD4*
      CD4_A
    d/dt(TcellDeg) <- +kdeg_CD8_N*CD8_N + kdeg_CD8_EM*CD8_EM + kdeg_CD8_A*CD8_A + kdeg_CD4_N*CD4_N + kdeg_CD4_EM*
      CD4_EM + kdeg_CD4_A*CD4_A
    d/dt(MMcellSyn) <- +kpr_MM*(1 - MM/C_M_MM)*MM
    d/dt(MMcellDeg) <- +kdeg_MM*MM + kkillMM_CD4*S_CD4_A_CD28_MM_CD38 + kkillMM_CD4*S_CD4_A_CD3_MM_CD28 +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD38 + kkillMM_CD8*S_CD8_A_CD28_MM_CD38 + kkillMM_CD8*S_CD8_A_CD3_MM_CD28 +
       kkillMM_CD8*S_CD8_A_CD3_MM_CD38
    d/dt(TRGTcellSyn) <- 0
    d/dt(TRGTcellDeg) <- +kdeg_TRGT*TRGT + kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38 + kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38 +
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38 + kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38
    d/dt(CD3Syn) <- +ksyn_CD8_N*CD3per_CD8 + ksyn_CD8_EM*CD3per_CD8 + kpr_CD8*CD8_A*CD3per_CD8 + ksyn_CD4_N*
      CD3per_CD4 + ksyn_CD4_EM*CD3per_CD4 + kpr_CD4*CD4_A*CD3per_CD4
    d/dt(CD3Deg) <- +kdeg_CD8_N*(R_CD8_N_CD3 + R_CD8_N_CD3_tsAb) + kdeg_CD8_EM*(R_CD8_EM_CD3 + R_CD8_EM_CD3_tsAb) +
       kdeg_CD8_A*(R_CD8_A_CD3 + R_CD8_A_CD3_tsAb) + kdeg_CD4_N*(R_CD4_N_CD3 + R_CD4_N_CD3_tsAb) + kdeg_CD4_EM*
      (R_CD4_EM_CD3 + R_CD4_EM_CD3_tsAb) + kdeg_CD4_A*(R_CD4_A_CD3 + R_CD4_A_CD3_tsAb)
    d/dt(CD28Syn) <- +ksyn_CD8_N*CD28per_CD8 + ksyn_CD8_EM*CD28per_CD8 + kpr_CD8*CD8_A*CD28per_CD8 + ksyn_CD4_N*
      CD28per_CD4 + ksyn_CD4_EM*CD28per_CD4 + kpr_CD4*CD4_A*CD28per_CD4 + kpr_MM*(1 - MM/C_M_MM)*MM*CD28per_MM
    d/dt(CD28Deg) <- +kdeg_CD8_N*(R_CD8_N_CD28 + R_CD8_N_CD28_tsAb) + kdeg_CD8_EM*(R_CD8_EM_CD28 + R_CD8_EM_CD28_tsAb) +
       kdeg_CD8_A*(R_CD8_A_CD28 + R_CD8_A_CD28_tsAb) + kdeg_CD4_N*(R_CD4_N_CD28 + R_CD4_N_CD28_tsAb) +
       kdeg_CD4_EM*(R_CD4_EM_CD28 + R_CD4_EM_CD28_tsAb) + kdeg_CD4_A*(R_CD4_A_CD28 + R_CD4_A_CD28_tsAb) +
       kdeg_MM*(R_MM_CD28 + R_MM_CD28_tsAb) + kdeg_TRGT*(R_TRGT_CD28 + R_TRGT_CD28_tsAb) + kkillMM_CD4*
      (S_CD4_A_CD28_MM_CD38_R_MM_CD28 + S_CD4_A_CD28_MM_CD38_R_MM_CD28_tsAb) + kkillTRGT_CD4*(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28 +
       S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) + kkillMM_CD4*(S_CD4_A_CD3_MM_CD28_R_MM_CD28 + S_CD4_A_CD3_MM_CD28_R_MM_CD28_tsAb) +
       kkillMM_CD4*S_CD4_A_CD3_MM_CD28_Br + kkillMM_CD4*(S_CD4_A_CD3_MM_CD38_R_MM_CD28 + S_CD4_A_CD3_MM_CD38_R_MM_CD28_tsAb) +
       kkillTRGT_CD4*(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28 + S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb) + kkillMM_CD8*
      (S_CD8_A_CD28_MM_CD38_R_MM_CD28 + S_CD8_A_CD28_MM_CD38_R_MM_CD28_tsAb) + kkillTRGT_CD8*(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28 +
       S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD28_tsAb) + kkillMM_CD8*(S_CD8_A_CD3_MM_CD28_R_MM_CD28 + S_CD8_A_CD3_MM_CD28_R_MM_CD28_tsAb) +
       kkillMM_CD8*S_CD8_A_CD3_MM_CD28_Br + kkillMM_CD8*(S_CD8_A_CD3_MM_CD38_R_MM_CD28 + S_CD8_A_CD3_MM_CD38_R_MM_CD28_tsAb) +
       kkillTRGT_CD8*(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28 + S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD28_tsAb)
    d/dt(CD38Syn) <- +kpr_MM*(1 - MM/C_M_MM)*MM*CD38per_MM
    d/dt(CD38Deg) <- +kdeg_MM*(R_MM_CD38 + R_MM_CD38_tsAb) + kdeg_TRGT*(R_TRGT_CD38 + R_TRGT_CD38_tsAb) +
       kdeg_sCD38*(R_sCD38_CD38 + R_sCD38_CD38_tsAb) + kkillMM_CD4*(S_CD4_A_CD28_MM_CD38_R_MM_CD38 + S_CD4_A_CD28_MM_CD38_R_MM_CD38_tsAb) +
       kkillMM_CD4*S_CD4_A_CD28_MM_CD38_Br + kkillTRGT_CD4*(S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38 + S_CD4_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) +
       kkillTRGT_CD4*S_CD4_A_CD28_TRGT_CD38_Br + kkillMM_CD4*(S_CD4_A_CD3_MM_CD28_R_MM_CD38 + S_CD4_A_CD3_MM_CD28_R_MM_CD38_tsAb) +
       kkillMM_CD4*(S_CD4_A_CD3_MM_CD38_R_MM_CD38 + S_CD4_A_CD3_MM_CD38_R_MM_CD38_tsAb) + kkillMM_CD4*
      S_CD4_A_CD3_MM_CD38_Br + kkillTRGT_CD4*(S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38 + S_CD4_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) +
       kkillTRGT_CD4*S_CD4_A_CD3_TRGT_CD38_Br + kkillMM_CD8*(S_CD8_A_CD28_MM_CD38_R_MM_CD38 + S_CD8_A_CD28_MM_CD38_R_MM_CD38_tsAb) +
       kkillMM_CD8*S_CD8_A_CD28_MM_CD38_Br + kkillTRGT_CD8*(S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38 + S_CD8_A_CD28_TRGT_CD38_R_TRGT_CD38_tsAb) +
       kkillTRGT_CD8*S_CD8_A_CD28_TRGT_CD38_Br + kkillMM_CD8*(S_CD8_A_CD3_MM_CD28_R_MM_CD38 + S_CD8_A_CD3_MM_CD28_R_MM_CD38_tsAb) +
       kkillMM_CD8*(S_CD8_A_CD3_MM_CD38_R_MM_CD38 + S_CD8_A_CD3_MM_CD38_R_MM_CD38_tsAb) + kkillMM_CD8*
      S_CD8_A_CD3_MM_CD38_Br + kkillTRGT_CD8*(S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38 + S_CD8_A_CD3_TRGT_CD38_R_TRGT_CD38_tsAb) +
       kkillTRGT_CD8*S_CD8_A_CD3_TRGT_CD38_Br
    d/dt(SynFormation) <- +kf_CD8_N_CD3_CD8_N_CD28*CD8_N*CD8_N + kf_CD8_N_CD3_CD8_EM_CD28*CD8_N*CD8_EM +
       kf_CD8_N_CD3_CD8_A_CD28*CD8_N*CD8_A + kf_CD8_N_CD3_CD4_N_CD28*CD8_N*CD4_N + kf_CD8_N_CD3_CD4_EM_CD28*
      CD8_N*CD4_EM + kf_CD8_N_CD3_CD4_A_CD28*CD8_N*CD4_A + kf_CD8_N_CD3_MM_CD28*CD8_N*MM + kf_CD8_EM_CD3_CD8_N_CD28*
      CD8_EM*CD8_N + kf_CD8_EM_CD3_CD8_EM_CD28*CD8_EM*CD8_EM + kf_CD8_EM_CD3_CD8_A_CD28*CD8_EM*CD8_A +
       kf_CD8_EM_CD3_CD4_N_CD28*CD8_EM*CD4_N + kf_CD8_EM_CD3_CD4_EM_CD28*CD8_EM*CD4_EM + kf_CD8_EM_CD3_CD4_A_CD28*
      CD8_EM*CD4_A + kf_CD8_EM_CD3_MM_CD28*CD8_EM*MM + kf_CD8_A_CD3_CD8_N_CD28*CD8_A*CD8_N + kf_CD8_A_CD3_CD8_EM_CD28*
      CD8_A*CD8_EM + kf_CD8_A_CD3_CD8_A_CD28*CD8_A*CD8_A + kf_CD8_A_CD3_CD4_N_CD28*CD8_A*CD4_N + kf_CD8_A_CD3_CD4_EM_CD28*
      CD8_A*CD4_EM + kf_CD8_A_CD3_CD4_A_CD28*CD8_A*CD4_A + kf_CD8_A_CD3_MM_CD28*CD8_A*MM + kf_CD4_N_CD3_CD8_N_CD28*
      CD4_N*CD8_N + kf_CD4_N_CD3_CD8_EM_CD28*CD4_N*CD8_EM + kf_CD4_N_CD3_CD8_A_CD28*CD4_N*CD8_A + kf_CD4_N_CD3_CD4_N_CD28*
      CD4_N*CD4_N + kf_CD4_N_CD3_CD4_EM_CD28*CD4_N*CD4_EM + kf_CD4_N_CD3_CD4_A_CD28*CD4_N*CD4_A + kf_CD4_N_CD3_MM_CD28*
      CD4_N*MM + kf_CD4_EM_CD3_CD8_N_CD28*CD4_EM*CD8_N + kf_CD4_EM_CD3_CD8_EM_CD28*CD4_EM*CD8_EM + kf_CD4_EM_CD3_CD8_A_CD28*
      CD4_EM*CD8_A + kf_CD4_EM_CD3_CD4_N_CD28*CD4_EM*CD4_N + kf_CD4_EM_CD3_CD4_EM_CD28*CD4_EM*CD4_EM +
       kf_CD4_EM_CD3_CD4_A_CD28*CD4_EM*CD4_A + kf_CD4_EM_CD3_MM_CD28*CD4_EM*MM + kf_CD4_A_CD3_CD8_N_CD28*
      CD4_A*CD8_N + kf_CD4_A_CD3_CD8_EM_CD28*CD4_A*CD8_EM + kf_CD4_A_CD3_CD8_A_CD28*CD4_A*CD8_A + kf_CD4_A_CD3_CD4_N_CD28*
      CD4_A*CD4_N + kf_CD4_A_CD3_CD4_EM_CD28*CD4_A*CD4_EM + kf_CD4_A_CD3_CD4_A_CD28*CD4_A*CD4_A + kf_CD4_A_CD3_MM_CD28*
      CD4_A*MM + kf_CD8_N_CD3_TRGT_CD38*CD8_N*TRGT + kf_CD8_N_CD3_MM_CD38*CD8_N*MM + kf_CD8_EM_CD3_TRGT_CD38*
      CD8_EM*TRGT + kf_CD8_EM_CD3_MM_CD38*CD8_EM*MM + kf_CD8_A_CD3_TRGT_CD38*CD8_A*TRGT + kf_CD8_A_CD3_MM_CD38*
      CD8_A*MM + kf_CD4_N_CD3_TRGT_CD38*CD4_N*TRGT + kf_CD4_N_CD3_MM_CD38*CD4_N*MM + kf_CD4_EM_CD3_TRGT_CD38*
      CD4_EM*TRGT + kf_CD4_EM_CD3_MM_CD38*CD4_EM*MM + kf_CD4_A_CD3_TRGT_CD38*CD4_A*TRGT + kf_CD4_A_CD3_MM_CD38*
      CD4_A*MM + kf_CD8_N_CD28_TRGT_CD38*CD8_N*TRGT + kf_CD8_N_CD28_MM_CD38*CD8_N*MM + kf_CD8_EM_CD28_TRGT_CD38*
      CD8_EM*TRGT + kf_CD8_EM_CD28_MM_CD38*CD8_EM*MM + kf_CD8_A_CD28_TRGT_CD38*CD8_A*TRGT + kf_CD8_A_CD28_MM_CD38*
      CD8_A*MM + kf_CD4_N_CD28_TRGT_CD38*CD4_N*TRGT + kf_CD4_N_CD28_MM_CD38*CD4_N*MM + kf_CD4_EM_CD28_TRGT_CD38*
      CD4_EM*TRGT + kf_CD4_EM_CD28_MM_CD38*CD4_EM*MM + kf_CD4_A_CD28_TRGT_CD38*CD4_A*TRGT + kf_CD4_A_CD28_MM_CD38*
      CD4_A*MM + kf_MM_CD28_TRGT_CD38*MM*TRGT + kf_MM_CD28_MM_CD38*MM*MM

    # Initial conditions (Table S2). Free receptor numbers are cell number x
    # antigen density; each soluble CD38 molecule carries one CD38 site. All
    # other states (drug, synapses, bound receptors, cytokines, flux
    # trackers) start at zero.
    CD8_N(0) <- bl_CD8_N
    R_CD8_N_CD3(0) <- bl_CD8_N * CD3per_CD8
    R_CD8_N_CD28(0) <- bl_CD8_N * CD28per_CD8
    CD8_EM(0) <- bl_CD8_EM
    R_CD8_EM_CD3(0) <- bl_CD8_EM * CD3per_CD8
    R_CD8_EM_CD28(0) <- bl_CD8_EM * CD28per_CD8
    CD8_A(0) <- bl_CD8_A
    R_CD8_A_CD3(0) <- bl_CD8_A * CD3per_CD8
    R_CD8_A_CD28(0) <- bl_CD8_A * CD28per_CD8
    CD4_N(0) <- bl_CD4_N
    R_CD4_N_CD3(0) <- bl_CD4_N * CD3per_CD4
    R_CD4_N_CD28(0) <- bl_CD4_N * CD28per_CD4
    CD4_EM(0) <- bl_CD4_EM
    R_CD4_EM_CD3(0) <- bl_CD4_EM * CD3per_CD4
    R_CD4_EM_CD28(0) <- bl_CD4_EM * CD28per_CD4
    CD4_A(0) <- bl_CD4_A
    R_CD4_A_CD3(0) <- bl_CD4_A * CD3per_CD4
    R_CD4_A_CD28(0) <- bl_CD4_A * CD28per_CD4
    MM(0) <- bl_MM
    R_MM_CD38(0) <- bl_MM * CD38per_MM
    R_MM_CD28(0) <- bl_MM * CD28per_MM
    TRGT(0) <- bl_TRGT
    R_TRGT_CD38(0) <- bl_TRGT * CD38per_TRGT
    R_TRGT_CD28(0) <- bl_TRGT * CD28per_TRGT
    sCD38(0) <- bl_sCD38
    R_sCD38_CD38(0) <- bl_sCD38

    # Outputs. Free drug concentration in nM.
    Cc <- tsAb / (Vol * 6.02214076e14)
    # Total MM cells, free plus synapsed (an MM-MM synapse holds two).
    mm_total <- MM + S_CD4_A_CD28_MM_CD38 + S_CD4_A_CD3_MM_CD28 + S_CD4_A_CD3_MM_CD38 +
      S_CD4_EM_CD28_MM_CD38 + S_CD4_EM_CD3_MM_CD28 + S_CD4_EM_CD3_MM_CD38 + S_CD4_N_CD28_MM_CD38 +
      S_CD4_N_CD3_MM_CD28 + S_CD4_N_CD3_MM_CD38 + S_CD8_A_CD28_MM_CD38 + S_CD8_A_CD3_MM_CD28 +
      S_CD8_A_CD3_MM_CD38 + S_CD8_EM_CD28_MM_CD38 + S_CD8_EM_CD3_MM_CD28 + S_CD8_EM_CD3_MM_CD38 +
      S_CD8_N_CD28_MM_CD38 + S_CD8_N_CD3_MM_CD28 + S_CD8_N_CD3_MM_CD38 + 2 * S_MM_CD28_MM_CD38 +
      S_MM_CD28_TRGT_CD38 + S_CD4_A_CD28_MM_CD38_MUT + S_CD4_A_CD3_MM_CD28_MUT +
      S_CD4_A_CD3_MM_CD38_MUT + S_CD8_A_CD28_MM_CD38_MUT + S_CD8_A_CD3_MM_CD28_MUT +
      S_CD8_A_CD3_MM_CD38_MUT
    # MM cell killing as reported in Figure 3A: loss from the initial count.
    pct_mm_killed <- 100 * (1 - mm_total / bl_MM)
    # Active T cells, free plus synapsed (Figure 3B) and free only (Figure 3C).
    tact_total <- CD8_A + CD4_A + S_CD4_A_CD28_MM_CD38 + S_CD4_A_CD28_TRGT_CD38 +
      2 * S_CD4_A_CD3_CD4_A_CD28 + S_CD4_A_CD3_CD4_EM_CD28 + S_CD4_A_CD3_CD4_N_CD28 +
      2 * S_CD4_A_CD3_CD8_A_CD28 + S_CD4_A_CD3_CD8_EM_CD28 + S_CD4_A_CD3_CD8_N_CD28 +
      S_CD4_A_CD3_MM_CD28 + S_CD4_A_CD3_MM_CD38 + S_CD4_A_CD3_TRGT_CD38 + S_CD4_EM_CD3_CD4_A_CD28 +
      S_CD4_EM_CD3_CD8_A_CD28 + S_CD4_N_CD3_CD4_A_CD28 + S_CD4_N_CD3_CD8_A_CD28 +
      S_CD8_A_CD28_MM_CD38 + S_CD8_A_CD28_TRGT_CD38 + 2 * S_CD8_A_CD3_CD4_A_CD28 +
      S_CD8_A_CD3_CD4_EM_CD28 + S_CD8_A_CD3_CD4_N_CD28 + 2 * S_CD8_A_CD3_CD8_A_CD28 +
      S_CD8_A_CD3_CD8_EM_CD28 + S_CD8_A_CD3_CD8_N_CD28 + S_CD8_A_CD3_MM_CD28 + S_CD8_A_CD3_MM_CD38 +
      S_CD8_A_CD3_TRGT_CD38 + S_CD8_EM_CD3_CD4_A_CD28 + S_CD8_EM_CD3_CD8_A_CD28 +
      S_CD8_N_CD3_CD4_A_CD28 + S_CD8_N_CD3_CD8_A_CD28 + S_CD4_A_CD28_MM_CD38_MUT +
      S_CD4_A_CD3_MM_CD28_MUT + S_CD4_A_CD3_MM_CD38_MUT + S_CD8_A_CD28_MM_CD38_MUT +
      S_CD8_A_CD3_MM_CD28_MUT + S_CD8_A_CD3_MM_CD38_MUT
    tact_free <- CD8_A + CD4_A
    # Ineffective MM synapses (the generator's SYN_MMiec list) as a percentage
    # of the initial MM cell number (Figure 3D).
    mm_ineff <- S_CD4_EM_CD28_MM_CD38 + S_CD4_N_CD28_MM_CD38 + S_CD8_EM_CD28_MM_CD38 +
      S_CD8_N_CD28_MM_CD38 + S_MM_CD28_MM_CD38 + S_MM_CD28_TRGT_CD38
    pct_mm_ineff <- 100 * mm_ineff / bl_MM
  })
}
