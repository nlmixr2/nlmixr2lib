Nakamura_2024_linagliptin_human_pbpk <- function() {
  description <- paste(
    "PBPK-TMDD (whole-body, 29-tissue, 5-compartment tandem liver,",
    "3-segment segregated-flow intestine). Linagliptin in healthy adults",
    "after single intravenous (0.5-10 mg, 1.5 h infusion) and single oral",
    "(1-10 mg) doses, fitted to plasma concentration, cumulative urinary",
    "and faecal excretion with the cluster Gauss-Newton method (CGNM).",
    "The model explains linagliptins target-mediated, nonlinear PK by",
    "reversible time-dependent binding (kon = koff/Kd) to soluble plasma",
    "DPP-4 (AsDPP4) and to transmembrane DPP-4 (AtDPP4) distributed over",
    "the organs by the tissue DPP-4 concentration ratios estimated from",
    "rat autoradiography, plus a DPP-4-complex-mediated saturable renal",
    "reabsorption clearance (CLreabDPP4) of the drug-tDPP4 complex in the",
    "renal tubule. Each tissue carries free drug, free sDPP-4, free",
    "tDPP-4 and the two drug-DPP-4 complexes as separate states (112 ODEs",
    "total). The 16 parameters without fixed() are the human CGNM medians",
    "of Table 2; every physiological volume, flow, partition coefficient",
    "and the rat-derived tissue DPP-4 concentrations are fixed() from",
    "Table S2 and the deposited model code (Data S2). CGNM is a",
    "fixed-effects nonlinear-least-squares method, so the model has no",
    "between-subject variability and the paper reports no residual-error",
    "model -- the propSd term is a placeholder. See the vignette Errata."
  )
  reference <- paste(
    "Nakamura R, Yoshikado T, Aoki Y, Sugiyama Y, Chiba K.",
    "Elucidation of DPP-4 involvement in systemic distribution and renal",
    "reabsorption of linagliptin by PBPK modeling with a cluster",
    "Gauss-Newton method. Clin Transl Sci. 2024;17(10):e70047.",
    "doi:10.1111/cts.70047.",
    "The ODE system is transcribed from the deposited human model code",
    "(Data S2, Program/pbpk_lng_human/lng_pbpk_model_human.r); the",
    "estimated parameters are the median column of Table 2 and the fixed",
    "physiological parameters are Table S2 and the deposited paramdata.",
    sep = " "
  )
  vignette <- "Nakamura_2024_linagliptin_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # GIMI (gastrointestinal-motility-inhibition) indicator. The paper
  # slows intestinal transit by the factor rAlpha until the colon mean
  # residence time MRTcolon after dosing (Supplementary Text, "Drug effect
  # on gastrointestinal motility"). The deposited event table keys the
  # window on the dose level, not the route: GIMI = 1 for t < MRTcolon in
  # every 1, 2.5, 5 and 10 mg arm, intravenous and oral alike, and GIMI = 0
  # throughout the 0.5 mg arm. It is a time-varying structural switch
  # rather than a fitted covariate effect, and is not in the canonical
  # covariate register; see the vignette Assumptions and deviations section.
  covariateData <- list(
    GIMI = list(
      description = paste(
        "1 = within the post-dose colon mean-residence-time window, where",
        "linagliptin inhibits gastrointestinal motility and every intestinal",
        "transit rate constant is scaled by rAlpha; 0 = normal motility."
      ),
      units = "binary",
      type = "categorical",
      reference_category = 0,
      notes = paste(
        "Time-varying per record. Set GIMI = 1 for records with time <",
        "MRTcolon, which depends on the dose level (2.00, 5.23, 10.6 and",
        "18.4 h for 1, 2.5, 5 and 10 mg) and applies to intravenous and oral",
        "arms alike; GIMI = 0 after MRTcolon and throughout a 0.5 mg dose,",
        "for which the authors define no window. Paper-specific motility",
        "switch, not a canonical register name."
      ),
      source_name = "GIMI"
    )
  )

  # PBPK organ-amount and DPP-4 (free / complex) states, kept with the
  # deposited Data S2 names (uppercase A-prefixed amounts) so the ODE
  # block can be verified line-by-line against the authors' code. They
  # are declared here so checkModelConventions() does not flag them as
  # non-canonical compartment names.
  paper_specific_compartments <- c(
    "Aehc_1",
    "Aehc_2",
    "Aehc_3",
    "Aserosa",
    "AsDPP4c_serosa",
    "AsDPP4f_serosa",
    "AtDPP4c_serosa",
    "AtDPP4f_serosa",
    "Amucosal_blood1",
    "AsDPP4c_mucosal_blood1",
    "AsDPP4f_mucosal_blood1",
    "AtDPP4c_mucosal_blood1",
    "AtDPP4f_mucosal_blood1",
    "Aduodenum_enterocyte",
    "Aduodenum_lumen",
    "AtDPP4c_duodenum_lumen",
    "AtDPP4f_duodenum_lumen",
    "Amucosal_blood2",
    "AsDPP4c_mucosal_blood2",
    "AsDPP4f_mucosal_blood2",
    "AtDPP4c_mucosal_blood2",
    "AtDPP4f_mucosal_blood2",
    "Ajejunum_enterocyte",
    "Ajejunum_lumen",
    "AtDPP4c_jejunum_lumen",
    "AtDPP4f_jejunum_lumen",
    "Amucosal_blood3",
    "AsDPP4c_mucosal_blood3",
    "AsDPP4f_mucosal_blood3",
    "AtDPP4c_mucosal_blood3",
    "AtDPP4f_mucosal_blood3",
    "Aileum_enterocyte",
    "Aileum_lumen",
    "AtDPP4c_ileum_lumen",
    "AtDPP4f_ileum_lumen",
    "ABlood",
    "AsDPP4c_blood",
    "AsDPP4f_blood",
    "AtDPP4c_blood",
    "AtDPP4f_blood",
    "Alung",
    "AsDPP4c_lung",
    "AsDPP4f_lung",
    "AtDPP4c_lung",
    "AtDPP4f_lung",
    "Aliver_1",
    "AsDPP4c_hepat_1",
    "AsDPP4f_hepat_1",
    "AtDPP4c_hepat_1",
    "AtDPP4f_hepat_1",
    "Aliver_2",
    "AsDPP4c_hepat_2",
    "AsDPP4f_hepat_2",
    "AtDPP4c_hepat_2",
    "AtDPP4f_hepat_2",
    "Aliver_3",
    "AsDPP4c_hepat_3",
    "AsDPP4f_hepat_3",
    "AtDPP4c_hepat_3",
    "AtDPP4f_hepat_3",
    "Aliver_4",
    "AsDPP4c_hepat_4",
    "AsDPP4f_hepat_4",
    "AtDPP4c_hepat_4",
    "AtDPP4f_hepat_4",
    "Aliver_5",
    "AsDPP4c_hepat_5",
    "AsDPP4f_hepat_5",
    "AtDPP4c_hepat_5",
    "AtDPP4f_hepat_5",
    "Areproductiveorgans",
    "AsDPP4c_reproductiveorgans",
    "AsDPP4f_reproductiveorgans",
    "AtDPP4c_reproductiveorgans",
    "AtDPP4f_reproductiveorgans",
    "Athymus",
    "AsDPP4c_thymus",
    "AsDPP4f_thymus",
    "AtDPP4c_thymus",
    "AtDPP4f_thymus",
    "Asalivarygland",
    "AsDPP4c_salivarygland",
    "AsDPP4f_salivarygland",
    "AtDPP4c_salivarygland",
    "AtDPP4f_salivarygland",
    "Apancreas",
    "AsDPP4c_pancreas",
    "AsDPP4f_pancreas",
    "Aspleen",
    "AsDPP4c_spleen",
    "AsDPP4f_spleen",
    "AtDPP4c_spleen",
    "AtDPP4f_spleen",
    "Askin",
    "AsDPP4c_skin",
    "AsDPP4f_skin",
    "AtDPP4c_skin",
    "AtDPP4f_skin",
    "Aremainder",
    "AsDPP4c_remainder",
    "AsDPP4f_remainder",
    "Akves",
    "AsDPP4c_kidneyVessel",
    "AsDPP4f_kidneyVessel",
    "AtDPP4c_kidneyVessel",
    "AtDPP4f_kidneyVessel",
    "Aktubule",
    "AtDPP4c_kidneytubule",
    "AtDPP4f_kidneytubule",
    "Xmetab",
    "Xurine",
    "Xfeces"
  )

  compartmentData <- list(
    Aehc_1 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aehc_2 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aehc_3 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aserosa = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_serosa = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_serosa = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_serosa = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_serosa = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Amucosal_blood1 = list(analyte = "linagliptin", units = "nmol", specimen = "whole blood", verified = FALSE),
    AsDPP4c_mucosal_blood1 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_mucosal_blood1 = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4c_mucosal_blood1 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_mucosal_blood1 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Aduodenum_enterocyte = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aduodenum_lumen = list(analyte = "linagliptin", units = "nmol", specimen = "administration site", verified = FALSE),
    AtDPP4c_duodenum_lumen = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    AtDPP4f_duodenum_lumen = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    Amucosal_blood2 = list(analyte = "linagliptin", units = "nmol", specimen = "whole blood", verified = FALSE),
    AsDPP4c_mucosal_blood2 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_mucosal_blood2 = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4c_mucosal_blood2 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_mucosal_blood2 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Ajejunum_enterocyte = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Ajejunum_lumen = list(analyte = "linagliptin", units = "nmol", specimen = "administration site", verified = FALSE),
    AtDPP4c_jejunum_lumen = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    AtDPP4f_jejunum_lumen = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    Amucosal_blood3 = list(analyte = "linagliptin", units = "nmol", specimen = "whole blood", verified = FALSE),
    AsDPP4c_mucosal_blood3 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_mucosal_blood3 = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4c_mucosal_blood3 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_mucosal_blood3 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Aileum_enterocyte = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aileum_lumen = list(analyte = "linagliptin", units = "nmol", specimen = "administration site", verified = FALSE),
    AtDPP4c_ileum_lumen = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    AtDPP4f_ileum_lumen = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    ABlood = list(analyte = "linagliptin", units = "nmol", specimen = "whole blood", verified = FALSE),
    AsDPP4c_blood = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_blood = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "whole blood", verified = FALSE),
    AtDPP4c_blood = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_blood = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Alung = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_lung = list(analyte = "linagliptin-sDPP-4 complex", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4f_lung = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_lung = list(analyte = "linagliptin-tDPP-4 complex", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4f_lung = list(analyte = "transmembrane DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    Aliver_1 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_hepat_1 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_hepat_1 = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_hepat_1 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_hepat_1 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Aliver_2 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_hepat_2 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_hepat_2 = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_hepat_2 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_hepat_2 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Aliver_3 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_hepat_3 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_hepat_3 = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_hepat_3 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_hepat_3 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Aliver_4 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_hepat_4 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_hepat_4 = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_hepat_4 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_hepat_4 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Aliver_5 = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_hepat_5 = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_hepat_5 = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_hepat_5 = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_hepat_5 = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Areproductiveorgans = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_reproductiveorgans = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_reproductiveorgans = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4c_reproductiveorgans = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_reproductiveorgans = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Athymus = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_thymus = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_thymus = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_thymus = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_thymus = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Asalivarygland = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_salivarygland = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_salivarygland = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4c_salivarygland = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_salivarygland = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Apancreas = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_pancreas = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_pancreas = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    Aspleen = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_spleen = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_spleen = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_spleen = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_spleen = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Askin = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_skin = list(analyte = "linagliptin-sDPP-4 complex", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4f_skin = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_skin = list(analyte = "linagliptin-tDPP-4 complex", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4f_skin = list(analyte = "transmembrane DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    Aremainder = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_remainder = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AsDPP4f_remainder = list(analyte = "soluble DPP-4 (free)", units = "nmol", specimen = "tissue", verified = FALSE),
    Akves = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AsDPP4c_kidneyVessel = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_kidneyVessel = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4c_kidneyVessel = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_kidneyVessel = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Aktubule = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    AtDPP4c_kidneytubule = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    AtDPP4f_kidneytubule = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    Xmetab = list(analyte = "linagliptin metabolites", units = "nmol", specimen = "not applicable", verified = FALSE),
    Xurine = list(analyte = "linagliptin", units = "nmol", specimen = "urine", verified = FALSE),
    Xfeces = list(analyte = "linagliptin", units = "nmol", specimen = "faeces", verified = FALSE)
  )

  population <- list(
    species = "human",
    disease_state = "healthy adults",
    dose_range = paste(
      "single intravenous 0.5, 2.5, 5 and 10 mg (1.5 h infusion) and",
      "single oral 1, 2.5, 5 and 10 mg"
    ),
    notes = paste(
      "The human clinical data (plasma concentration, urinary and faecal",
      "excretion after intravenous and oral linagliptin) are the pooled",
      "healthy-volunteer data of Sarashina 2010 / Retlich 2010 / Blech",
      "2010, the same data analysed in Sarashina 2020. Individual subject",
      "counts, ages, weights and sex are not reported by Nakamura 2024;",
      "all physiological parameters are scaled to a single typical 70 kg",
      "adult. Doses are supplied in nmol (linagliptin MW 472.54; the",
      "deposited code converts mg to nmol with the factor 2116.22 nmol/mg)."
    )
  )

  ini({
    lVBlood <- fixed(log(5.2)); label("Blood volume (L)")  # Table S2 blood volume 5.2
    lVduodenum_enterocyte <- fixed(log(0.0588)); label("Duodenum Enterocyte volume (L)")  # Table S2 duodenum enterocyte volume 0.0588
    lVduodenum_lumen <- fixed(log(0.0651)); label("Duodenum Lumen volume (L)")  # Table S2 duodenum lumen volume 0.0651
    lVeh <- fixed(log(0.469)); label("Extra-Hepatic (Disse/Interstitial) volume (L)")  # Table S2 extra hepatic volume 0.469
    lVh <- fixed(log(1.689)); label("Liver volume (L)")  # Table S2 liver volume 1.689
    lVileum_enterocyte <- fixed(log(0.2499)); label("Ileum Enterocyte volume (L)")  # Table S2 ileum enterocyte volume 0.2499
    lVileum_lumen <- fixed(log(0.2765)); label("Ileum Lumen volume (L)")  # Table S2 ileum lumen volume 0.2765
    lVjejunum_enterocyte <- fixed(log(0.2093)); label("Jejunum Enterocyte volume (L)")  # Table S2 jejunum enterocyte volume 0.2093
    lVjejunum_lumen <- fixed(log(0.2317)); label("Jejunum Lumen volume (L)")  # Table S2 jejunum lumen volume 0.2317
    lVktubule <- fixed(log(0.127)); label("Kidney Tubule (Duct) volume (L)")  # Table S2 kidney tubule(duct) volume 0.127
    lVkves <- fixed(log(0.153)); label("Kidney Vessel+cell volume (L)")  # Table S2 kidney vessel volume 0.153
    lVl <- fixed(log(1.17)); label("Lung volume (L)")  # Table S2 lung volume 1.17
    lVmucosal_blood1 <- fixed(log(0.00791)); label("Duodenum Mucosal Blood volume (L)")  # Table S2 duodenum mucosal blood volume 0.00791
    lVmucosal_blood2 <- fixed(log(0.028)); label("Jejunum Mucosal Blood volume (L)")  # Table S2 jejunum mucosal blood volume 0.028
    lVmucosal_blood3 <- fixed(log(0.03339)); label("Ileum Mucosal Blood volume (L)")  # Table S2 ileum mucosal blood volume 0.03339
    lVpan <- fixed(log(0.134)); label("Pancreas volume (L)")  # Table S2 pancreas volume 0.134
    lVrem <- fixed(log(65.744)); label("Remainder volume (L)")  # Table S2 remainder volume 65.744
    lVrep <- fixed(log(0.0374)); label("Reproductive Organs volume (L)")  # Table S2 reproductive organs volume 0.0374
    lVs <- fixed(log(7.8)); label("Skin volume (L)")  # Table S2 skin volume 7.8
    lVsal <- fixed(log(0.081)); label("Salivary Gland volume (L)")  # Table S2 salivary gland volume 0.081
    lVserosa <- fixed(log(0.4914)); label("Serosa volume (L)")  # Table S2 serosa volume 0.4914
    lVspl <- fixed(log(0.192)); label("Spleen volume (L)")  # Table S2 spleen volume 0.192
    lVthy <- fixed(log(0.024)); label("Thymus volume (L)")  # Table S2 thymus volume 0.024
    fVkvesv <- fixed(0.105); label("Kidney fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of kidney vessel volume 0.105
    fVlv <- fixed(0.262); label("Lung fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of lung vessel volume 0.262
    fVpanv <- fixed(0.1); label("Pancreas fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of pancreas vessel volume 0.1
    fVremv <- fixed(0.1); label("Remainder fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of remainder vessel volume 0.1
    fVrepv <- fixed(0.02); label("Reproductive Organs fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of reproductive organs vessel volume 0.02
    fVsalv <- fixed(0.1); label("Salivary Gland fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of salivary gland vessel volume 0.1
    fVserosav <- fixed(0.1); label("Serosa fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of serosa vessel volume 0.1
    fVsplv <- fixed(0.282); label("Spleen fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of spleen vessel volume 0.282
    fVsv <- fixed(0.019); label("Skin fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of skin vessel volume 0.019
    fVthyv <- fixed(0.03); label("Thymus fractional capillary (vessel) volume (unitless)")  # Table S2 fraction of thymus vessel volume 0.03
    lGFR <- fixed(log(7.5)); label("Glomerular filtration rate (L/h)")  # Table S2 glomerular filtration rate 7.5
    lQduode_muc <- fixed(log(1.148)); label("Duodenum Mucosal blood flow (L/h)")  # Table S2 duodenum mucosal blood flow 1.148
    lQh <- fixed(log(87)); label("Liver blood flow (L/h)")  # Table S2 liver blood flow 87
    lQha_proper <- fixed(log(13.72)); label("Proper Hepatic Artery blood flow (L/h)")  # Table S2 hepatic artery blood flow 13.72
    lQile_muc <- fixed(log(5.572)); label("Ileum Mucosal blood flow (L/h)")  # Table S2 ileum mucosal blood flow 5.572
    lQjeju_muc <- fixed(log(6.3)); label("Jejunum Mucosal blood flow (L/h)")  # Table S2 jejunum mucosal blood flow 6.3
    lQk <- fixed(log(74.4)); label("Kidney blood flow (L/h)")  # Table S2 kidney blood flow 74.4
    lQl <- fixed(log(262)); label("Lung blood flow (L/h)")  # Table S2 lung blood flow 262
    lQpan <- fixed(log(3.36)); label("Pancreas blood flow (L/h)")  # Table S2 pancreas blood flow 3.36
    lQrem <- fixed(log(155.1)); label("Remainder blood flow (L/h)")  # Table S2 remainder blood flow 155.1
    lQrep <- fixed(log(1.04)); label("Reproductive Organs blood flow (L/h)")  # Table S2 reproductive organs blood flow 1.04
    lQs <- fixed(log(18)); label("Skin blood flow (L/h)")  # Table S2 skin blood flow 18
    lQsal <- fixed(log(0.104)); label("Salivary Gland blood flow (L/h)")  # Table S2 salivary gland blood flow 0.104
    lQserosa <- fixed(log(52.08)); label("Serosa blood flow (L/h)")  # Table S2 serosa blood flow 52.08
    lQspl <- fixed(log(4.62)); label("Spleen blood flow (L/h)")  # Table S2 spleen blood flow 4.62
    lQthy <- fixed(log(0.201)); label("Thymus blood flow (L/h)")  # Table S2 thymus blood flow 0.201
    lQu <- fixed(log(0.05833333)); label("Urine blood flow (L/h)")  # Table S2 urin flow 0.05833333
    lKpa_calc <- fixed(log(0.201)); label("Adipose tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of adipose 0.201
    lKpgut_calc <- fixed(log(1.181)); label("Gut tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of gut 1.181
    lKph_calc <- fixed(log(1.915)); label("Liver tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of liver 1.915
    lKphrt_calc <- fixed(log(1.00752594533395)); label("Heart tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of heart 1.007526
    lKpk_calc <- fixed(log(2.058)); label("Kidney tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of kidney 2.058
    lKpl_calc <- fixed(log(1.59)); label("Lung tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of lung 1.59
    lKpm_calc <- fixed(log(0.913)); label("Muscle tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of muscle 0.913
    lKppan_calc <- fixed(log(0.92286251315842)); label("Pancreas tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of pancreas 0.9228625
    lKps_calc <- fixed(log(0.726)); label("Skin tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of skin 0.726
    lKpspl_calc <- fixed(log(1.36651167316662)); label("Spleen tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of spleen 1.366512
    lKpthy_calc <- fixed(log(1.1888252721636)); label("Thymus tissue:plasma partition coefficient (unitless)")  # Table S2 Kp value of thymus 1.188825
    fgut <- fixed(0.232); label("Free fraction in enterocyte (unitless)")  # Table S2 free fraction in enterocyte 0.232
    fupAlb <- fixed(0.227); label("Non-saturable plasma free fraction (albumin binding) (unitless)")  # Table S2 free fraction of 2nd binding protein (not saturable) 0.227
    Ht <- fixed(0.44); label("Haematocrit (unitless)")  # Table S2 bloodcell plasma volume ratio, hematocrit 0.44
    lkbile <- fixed(log(0.605)); label("Biliary transit rate constant (1/h)")  # Table S2 biliary transit rate constant 0.605
    lKmPgp <- fixed(log(187000)); label("Michaelis constant for intestinal P-gp (nmol/L)")  # Table S2 Km for P-gp in human 187000
    RB <- fixed(0.691); label("Blood-to-plasma concentration ratio (unitless)")  # Table S2 blood to plasma ratio 0.691
    lCDPP4_skin <- fixed(log(33.98647)); label("Transmembrane DPP-4 concentration in skin (nmol/L organ)")  # Data S2 human Fixed_param_temporary; rat-derived tissue DPP-4 concentration (nmol/L organ), Figure 1e
    lCDPP4_lung <- fixed(log(477.07926)); label("Transmembrane DPP-4 concentration in lung (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_liver <- fixed(log(715.28117)); label("Transmembrane DPP-4 concentration in liver (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_kidney <- fixed(log(12003.79348)); label("Transmembrane DPP-4 concentration in kidney (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_spleen <- fixed(log(281.01567)); label("Transmembrane DPP-4 concentration in spleen (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_salivarygland <- fixed(log(120.95336)); label("Transmembrane DPP-4 concentration in salivary gland (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_thymus <- fixed(log(174.63852)); label("Transmembrane DPP-4 concentration in thymus (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lCDPP4_reproductiveorgans <- fixed(log(176.00888)); label("Transmembrane DPP-4 concentration in reproductive organs (nmol/L organ)")  # Data S2; rat-derived tissue DPP-4 concentration (Figure 1e)
    lKprep_vivo <- fixed(log(2.128587)); label("Reproductive Organs tissue:plasma partition coefficient (unitless)")  # Table S2 Kp reproductive organ 2.13 (estimated in rat analysis)
    lKpsal_vivo <- fixed(log(1.634672)); label("Salivary Gland tissue:plasma partition coefficient (unitless)")  # Table S2 Kp salivary gland 1.63 (estimated in rat analysis)
    AR <- fixed(1); label("Apical reabsorption-ratio switch (unitless)")  # Data S2 AR=1 (apical reabsorption-ratio switch, off)
    Fdif_duodenum <- fixed(0.076); label("Duodenal fraction of small-intestinal diffusion area (unitless)")  # Data S2; duodenal fraction of small-intestinal passive-diffusion area (Pang segmented-flow)
    Fdif_jejunum <- fixed(0.395); label("Jejunal fraction of small-intestinal diffusion area (unitless)")  # Data S2; jejunal fraction of small-intestinal passive-diffusion area
    Fdif_ileum <- fixed(0.529); label("Ileal fraction of small-intestinal diffusion area (unitless)")  # Data S2; ileal fraction of small-intestinal passive-diffusion area
    Fpgp_duodenum <- fixed(0.057); label("Duodenal fraction of intestinal P-gp (unitless)")  # Data S2; duodenal fraction of intestinal P-gp
    Fpgp_jejunum <- fixed(0.273); label("Jejunal fraction of intestinal P-gp (unitless)")  # Data S2; jejunal fraction of intestinal P-gp
    Fpgp_ileum <- fixed(0.67); label("Ileal fraction of intestinal P-gp (unitless)")  # Data S2; ileal fraction of intestinal P-gp
    Gammaent <- fixed(1); label("Enterocyte basolateral:apical diffusion ratio (unitless)")  # Data S2 Gammaent=1 (enterocyte basolateral:apical diffusion ratio)
    lkfeces_duodenum <- fixed(log(5.12)); label("Duodenal luminal transit rate constant (1/h)")  # Data S2; duodenal luminal transit rate constant (1/h)
    lkfeces_ileum <- fixed(log(1.21)); label("Ileal luminal transit rate constant (1/h)")  # Data S2; ileal luminal transit rate constant (1/h)
    lkfeces_jejunum <- fixed(log(1.47)); label("Jejunal luminal transit rate constant (1/h)")  # Data S2; jejunal luminal transit rate constant (1/h)
    lKph_vivo <- fixed(log(13.07857)); label("Liver tissue:plasma partition coefficient (unitless)")  # Table S2 Kp liver 13.07857 (estimated in rat analysis)
    lCLinth <- log(407); label("Hepatic intrinsic clearance (L/h)")  # Table 2 CLint,h median 407 L/h
    lCLsecint <- log(61.5); label("Renal tubular secretion clearance (L/h)")  # Table 2 CLsec,int median 61.5 L/h
    lCLreabDPP4 <- log(0.0228); label("DPP-4-complex-mediated renal reabsorption clearance (L/h)")  # Table 2 CLreab,DPP4 median 0.0228 L/h
    lAsDPP4 <- log(31); label("Total amount of soluble DPP-4 (nmol)")  # Table 2 AsDPP4 median 31.0 nmol
    lAtDPP4 <- log(2610); label("Total amount of transmembrane DPP-4 (nmol)")  # Table 2 AtDPP4 median 2.61e3 nmol
    fbile <- 0.97; label("Fraction of hepatic intrinsic clearance that is biliary (unitless)")  # Table 2 fbile median 0.970
    lkoff <- log(0.229); label("LNG-DPP-4 dissociation rate constant (1/h)")  # Table 2 koff median 0.229 1/h
    lKd <- log(0.0259); label("LNG-DPP-4 dissociation constant (nmol/L)")  # Table 2 Kd median 0.0259 nM
    flumen_kidney <- 0.447; label("Fraction of kidney DPP-4 on the luminal side (unitless)")  # Table 2 flumen,kidney median 0.447
    lVmaxPgp <- log(10000000); label("Maximum intestinal P-gp transport rate (nmol/h)")  # Table 2 Vmax,P-gp median 1e7 nmol/h
    lPSdifentBE <- log(1.91); label("Mucosal-blood-to-enterocyte permeability clearance (L/h)")  # Table 2 PSdif,ent,BE median 1.91 L/h
    lLR <- log(6.22); label("Ratio of lumen-to-enterocyte vs blood-to-enterocyte permeability (unitless)")  # Table 2 LR median 6.22
    rAlpha <- 0.694; label("Intestinal-transit scaling factor under gastrointestinal-motility inhibition (unitless)")  # Table 2 rAlpha median 0.694
    lKprem_vivo <- log(1.85); label("Remainder tissue:plasma partition coefficient (unitless)")  # Table 2 Kp,rem.vivo median 1.85
    flumen_gut <- 0.56; label("Fraction of intestinal DPP-4 on the luminal side (unitless)")  # Table 2 flumen,gut median 0.56
    lrGL <- log(59.3); label("Intestine-to-liver DPP-4 amount ratio (unitless)")  # Table 2 rGL median 59.3
    propSd <- fixed(0.1); label("Proportional residual error (placeholder; not estimated)")  # placeholder: CGNM is fixed-effects least-squares; the paper reports no residual-error model
  })

  model({
    # Back-transform the log-scale parameters.
    VBlood <- exp(lVBlood)
    Vduodenum_enterocyte <- exp(lVduodenum_enterocyte)
    Vduodenum_lumen <- exp(lVduodenum_lumen)
    Veh <- exp(lVeh)
    Vh <- exp(lVh)
    Vileum_enterocyte <- exp(lVileum_enterocyte)
    Vileum_lumen <- exp(lVileum_lumen)
    Vjejunum_enterocyte <- exp(lVjejunum_enterocyte)
    Vjejunum_lumen <- exp(lVjejunum_lumen)
    Vktubule <- exp(lVktubule)
    Vkves <- exp(lVkves)
    Vl <- exp(lVl)
    Vmucosal_blood1 <- exp(lVmucosal_blood1)
    Vmucosal_blood2 <- exp(lVmucosal_blood2)
    Vmucosal_blood3 <- exp(lVmucosal_blood3)
    Vpan <- exp(lVpan)
    Vrem <- exp(lVrem)
    Vrep <- exp(lVrep)
    Vs <- exp(lVs)
    Vsal <- exp(lVsal)
    Vserosa <- exp(lVserosa)
    Vspl <- exp(lVspl)
    Vthy <- exp(lVthy)
    GFR <- exp(lGFR)
    Qduode_muc <- exp(lQduode_muc)
    Qh <- exp(lQh)
    Qha_proper <- exp(lQha_proper)
    Qile_muc <- exp(lQile_muc)
    Qjeju_muc <- exp(lQjeju_muc)
    Qk <- exp(lQk)
    Ql <- exp(lQl)
    Qpan <- exp(lQpan)
    Qrem <- exp(lQrem)
    Qrep <- exp(lQrep)
    Qs <- exp(lQs)
    Qsal <- exp(lQsal)
    Qserosa <- exp(lQserosa)
    Qspl <- exp(lQspl)
    Qthy <- exp(lQthy)
    Qu <- exp(lQu)
    Kpa_calc <- exp(lKpa_calc)
    Kpgut_calc <- exp(lKpgut_calc)
    Kph_calc <- exp(lKph_calc)
    Kphrt_calc <- exp(lKphrt_calc)
    Kpk_calc <- exp(lKpk_calc)
    Kpl_calc <- exp(lKpl_calc)
    Kpm_calc <- exp(lKpm_calc)
    Kppan_calc <- exp(lKppan_calc)
    Kps_calc <- exp(lKps_calc)
    Kpspl_calc <- exp(lKpspl_calc)
    Kpthy_calc <- exp(lKpthy_calc)
    kbile <- exp(lkbile)
    KmPgp <- exp(lKmPgp)
    CDPP4_skin <- exp(lCDPP4_skin)
    CDPP4_lung <- exp(lCDPP4_lung)
    CDPP4_liver <- exp(lCDPP4_liver)
    CDPP4_kidney <- exp(lCDPP4_kidney)
    CDPP4_spleen <- exp(lCDPP4_spleen)
    CDPP4_salivarygland <- exp(lCDPP4_salivarygland)
    CDPP4_thymus <- exp(lCDPP4_thymus)
    CDPP4_reproductiveorgans <- exp(lCDPP4_reproductiveorgans)
    Kprep_vivo <- exp(lKprep_vivo)
    Kpsal_vivo <- exp(lKpsal_vivo)
    kfeces_duodenum <- exp(lkfeces_duodenum)
    kfeces_ileum <- exp(lkfeces_ileum)
    kfeces_jejunum <- exp(lkfeces_jejunum)
    Kph_vivo <- exp(lKph_vivo)
    CLinth <- exp(lCLinth)
    CLsecint <- exp(lCLsecint)
    CLreabDPP4 <- exp(lCLreabDPP4)
    AsDPP4 <- exp(lAsDPP4)
    AtDPP4 <- exp(lAtDPP4)
    koff <- exp(lkoff)
    Kd <- exp(lKd)
    VmaxPgp <- exp(lVmaxPgp)
    PSdifentBE <- exp(lPSdifentBE)
    LR <- exp(lLR)
    Kprem_vivo <- exp(lKprem_vivo)
    rGL <- exp(lrGL)

  if (GIMI==0) {
    Alpha <- 1
  } else if (GIMI==1) {
    Alpha <- 1*rAlpha
  }
  
  # used in simulation
  Vplasma <- VBlood*(1-Ht)
  
  Vh_i <- Vh/5
  Veh_i <- Veh/5
  CLinth_i <- CLinth/5
  
  Vkvesv <- Vkves*fVkvesv
  Vsv <- Vs*fVsv
  Vlv <- Vl*fVlv
  
  Vpanv <- Vpan*fVpanv # pancreas
  Vsplv <- Vspl*fVsplv # spleen
  Vrepv <- Vrep*fVrepv # reproductiveorgans
  Vthyv <- Vthy*fVthyv # thymus
  Vsalv <- Vsal*fVsalv # salivarygland
  Vremv <- Vrem*fVremv # remainder
  Vserosav <- Vserosa*fVserosav # serosa
  
  # total gut volume
  Vgut <- Vduodenum_enterocyte+Vduodenum_lumen
  +Vjejunum_enterocyte+Vjejunum_lumen
  +Vileum_enterocyte+Vileum_lumen
  
  # DPP-4 binding
  kon <- koff/Kd
  
  # Kp values
  SFKp <- 1
  Kpk <- Kpk_calc*SFKp
  Kpspl <- Kpspl_calc*SFKp
  Kpm <- Kpm_calc*SFKp
  Kps <- Kps_calc*SFKp
  Kpa <- Kpa_calc*SFKp
  Kpl <- Kpl_calc*SFKp
  Kpthy <- Kpthy_calc*SFKp
  Kppan <- Kppan_calc*SFKp
  Kpgut <- Kpgut_calc*SFKp
  
  Kprep <- Kprep_vivo*SFKp
  Kprem <- Kprem_vivo*SFKp
  Kph <- Kph_vivo*SFKp
  
  Cplasma <- (ABlood/VBlood/RB)+(AsDPP4c_blood/Vplasma)
  
  #########################
  # Differential equation #
  #########################
  
  d/dt(Aehc_1) <- +CLinth_i*fbile*Aliver_1/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_2/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_3/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_4/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_5/Vh_i/Kph*fupAlb
  -kbile*Aehc_1
  
  d/dt(Aehc_2) <- +kbile*Aehc_1
  -kbile*Aehc_2
  
  d/dt(Aehc_3) <- +kbile*Aehc_2
  -kbile*Aehc_3
  
  
  d/dt(Aserosa) <- -kon*AsDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  +koff*AsDPP4c_serosa
  
  -kon*AtDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  +koff*AtDPP4c_serosa
  
  +Qserosa*Alung/Vl/Kpl*RB
  -Qserosa*Aserosa/Vserosa/Kpgut*RB
  
  
  d/dt(AsDPP4c_serosa) <- +kon*AsDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  -koff*AsDPP4c_serosa
  
  d/dt(AsDPP4f_serosa) <- -kon*AsDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  +koff*AsDPP4c_serosa
  
  d/dt(AtDPP4c_serosa) <- +kon*AtDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  -koff*AtDPP4c_serosa
  
  d/dt(AtDPP4f_serosa) <- -kon*AtDPP4f_serosa*Aserosa/Vserosa/Kpgut*fupAlb
  +koff*AtDPP4c_serosa
  
  
  d/dt(Amucosal_blood1) <- -kon*AsDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood1
  
  -kon*AtDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood1
  
  +Qduode_muc*Alung/Vl/Kpl*RB
  -Qduode_muc*Amucosal_blood1/Vmucosal_blood1
  +fgut*(PSdifentBE/Gammaent)*Fdif_duodenum*Aduodenum_enterocyte/Vduodenum_enterocyte
  -(fupAlb/RB)*PSdifentBE*Fdif_duodenum*Amucosal_blood1/Vmucosal_blood1
  
  d/dt(AsDPP4c_mucosal_blood1) <- +kon*AsDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  -koff*AsDPP4c_mucosal_blood1
  
  d/dt(AsDPP4f_mucosal_blood1) <- -kon*AsDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood1
  
  d/dt(AtDPP4c_mucosal_blood1) <- +kon*AtDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  -koff*AtDPP4c_mucosal_blood1
  
  d/dt(AtDPP4f_mucosal_blood1) <- -kon*AtDPP4f_mucosal_blood1*Amucosal_blood1/Vmucosal_blood1/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood1
  
  
  d/dt(Aduodenum_enterocyte) <- -fgut*(PSdifentBE/Gammaent)*Fdif_duodenum*Aduodenum_enterocyte/Vduodenum_enterocyte
  +(fupAlb/RB)*PSdifentBE*Fdif_duodenum*Amucosal_blood1/Vmucosal_blood1
  +LR*PSdifentBE*Fdif_duodenum*Aduodenum_lumen/Vduodenum_lumen
  -fgut*VmaxPgp*Fpgp_duodenum/(KmPgp+fgut*Aduodenum_enterocyte/Vduodenum_enterocyte)*Aduodenum_enterocyte/Vduodenum_enterocyte
  -fgut*AR*(PSdifentBE/Gammaent)*Fdif_duodenum*Aduodenum_enterocyte/Vduodenum_enterocyte
  
  
  d/dt(Aduodenum_lumen) <- -kon*AtDPP4f_duodenum_lumen*Aduodenum_lumen/Vduodenum_lumen
  +koff*AtDPP4c_duodenum_lumen
  
  +kbile*Aehc_3
  -Alpha*kfeces_duodenum*Aduodenum_lumen
  -LR*PSdifentBE*Fdif_duodenum*Aduodenum_lumen/Vduodenum_lumen
  +fgut*VmaxPgp*Fpgp_duodenum/(KmPgp+fgut*Aduodenum_enterocyte/Vduodenum_enterocyte)*Aduodenum_enterocyte/Vduodenum_enterocyte
  +fgut*AR*(PSdifentBE/Gammaent)*Fdif_duodenum*Aduodenum_enterocyte/Vduodenum_enterocyte
  
  d/dt(AtDPP4c_duodenum_lumen) <- +kon*AtDPP4f_duodenum_lumen*Aduodenum_lumen/Vduodenum_lumen
  -koff*AtDPP4c_duodenum_lumen
  
  d/dt(AtDPP4f_duodenum_lumen) <- -kon*AtDPP4f_duodenum_lumen*Aduodenum_lumen/Vduodenum_lumen
  +koff*AtDPP4c_duodenum_lumen
  
  
  d/dt(Amucosal_blood2) <- -kon*AsDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood2
  
  -kon*AtDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood2
  
  +Qjeju_muc*Alung/Vl/Kpl*RB
  -Qjeju_muc*Amucosal_blood2/Vmucosal_blood2
  +fgut*(PSdifentBE/Gammaent)*Fdif_jejunum*Ajejunum_enterocyte/Vjejunum_enterocyte
  -(fupAlb/RB)*PSdifentBE*Fdif_jejunum*Amucosal_blood2/Vmucosal_blood2
  
  d/dt(AsDPP4c_mucosal_blood2) <- +kon*AsDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  -koff*AsDPP4c_mucosal_blood2
  
  d/dt(AsDPP4f_mucosal_blood2) <- -kon*AsDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood2
  
  d/dt(AtDPP4c_mucosal_blood2) <- +kon*AtDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  -koff*AtDPP4c_mucosal_blood2
  
  d/dt(AtDPP4f_mucosal_blood2) <- -kon*AtDPP4f_mucosal_blood2*Amucosal_blood2/Vmucosal_blood2/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood2
  
  
  d/dt(Ajejunum_enterocyte) <- -fgut*(PSdifentBE/Gammaent)*Fdif_jejunum*Ajejunum_enterocyte/Vjejunum_enterocyte
  +(fupAlb/RB)*PSdifentBE*Fdif_jejunum*Amucosal_blood2/Vmucosal_blood2
  +LR*PSdifentBE*Fdif_jejunum*Ajejunum_lumen/Vjejunum_lumen
  -fgut*VmaxPgp*Fpgp_jejunum/(KmPgp+fgut*Ajejunum_enterocyte/Vjejunum_enterocyte)*Ajejunum_enterocyte/Vjejunum_enterocyte
  -fgut*AR*(PSdifentBE/Gammaent)*Fdif_jejunum*Ajejunum_enterocyte/Vjejunum_enterocyte
  
  
  d/dt(Ajejunum_lumen) <- -kon*AtDPP4f_jejunum_lumen*Ajejunum_lumen/Vjejunum_lumen
  +koff*AtDPP4c_jejunum_lumen
  
  +Alpha*kfeces_duodenum*Aduodenum_lumen
  -Alpha*kfeces_jejunum*Ajejunum_lumen
  -LR*PSdifentBE*Fdif_jejunum*Ajejunum_lumen/Vjejunum_lumen
  +fgut*VmaxPgp*Fpgp_jejunum/(KmPgp+fgut*Ajejunum_enterocyte/Vjejunum_enterocyte)*Ajejunum_enterocyte/Vjejunum_enterocyte
  +fgut*AR*(PSdifentBE/Gammaent)*Fdif_jejunum*Ajejunum_enterocyte/Vjejunum_enterocyte
  
  d/dt(AtDPP4c_jejunum_lumen) <- +kon*AtDPP4f_jejunum_lumen*Ajejunum_lumen/Vjejunum_lumen
  -koff*AtDPP4c_jejunum_lumen
  
  d/dt(AtDPP4f_jejunum_lumen) <- -kon*AtDPP4f_jejunum_lumen*Ajejunum_lumen/Vjejunum_lumen
  +koff*AtDPP4c_jejunum_lumen
  
  
  d/dt(Amucosal_blood3) <- -kon*AsDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood3
  
  -kon*AtDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood3
  
  +Qile_muc*Alung/Vl/Kpl*RB
  -Qile_muc*Amucosal_blood3/Vmucosal_blood3
  +fgut*(PSdifentBE/Gammaent)*Fdif_ileum*Aileum_enterocyte/Vileum_enterocyte
  -(fupAlb/RB)*PSdifentBE*Fdif_ileum*Amucosal_blood3/Vmucosal_blood3
  
  d/dt(AsDPP4c_mucosal_blood3) <- +kon*AsDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  -koff*AsDPP4c_mucosal_blood3
  
  d/dt(AsDPP4f_mucosal_blood3) <- -kon*AsDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood3
  
  d/dt(AtDPP4c_mucosal_blood3) <- +kon*AtDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  -koff*AtDPP4c_mucosal_blood3
  
  d/dt(AtDPP4f_mucosal_blood3) <- -kon*AtDPP4f_mucosal_blood3*Amucosal_blood3/Vmucosal_blood3/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood3
  
  
  d/dt(Aileum_enterocyte) <- -fgut*(PSdifentBE/Gammaent)*Fdif_ileum*Aileum_enterocyte/Vileum_enterocyte
  +(fupAlb/RB)*PSdifentBE*Fdif_ileum*Amucosal_blood3/Vmucosal_blood3
  +LR*PSdifentBE*Fdif_ileum*Aileum_lumen/Vileum_lumen
  -fgut*VmaxPgp*Fpgp_ileum/(KmPgp+fgut*Aileum_enterocyte/Vileum_enterocyte)*Aileum_enterocyte/Vileum_enterocyte
  -fgut*AR*(PSdifentBE/Gammaent)*Fdif_ileum*Aileum_enterocyte/Vileum_enterocyte
  
  
  d/dt(Aileum_lumen) <- -kon*AtDPP4f_ileum_lumen*Aileum_lumen/Vileum_lumen
  +koff*AtDPP4c_ileum_lumen
  
  +Alpha*kfeces_jejunum*Ajejunum_lumen
  -Alpha*kfeces_ileum*Aileum_lumen
  -LR*PSdifentBE*Fdif_ileum*Aileum_lumen/Vileum_lumen
  +fgut*VmaxPgp*Fpgp_ileum/(KmPgp+fgut*Aileum_enterocyte/Vileum_enterocyte)*Aileum_enterocyte/Vileum_enterocyte
  +fgut*AR*(PSdifentBE/Gammaent)*Fdif_ileum*Aileum_enterocyte/Vileum_enterocyte
  
  d/dt(AtDPP4c_ileum_lumen) <- +kon*AtDPP4f_ileum_lumen*Aileum_lumen/Vileum_lumen
  -koff*AtDPP4c_ileum_lumen
  
  d/dt(AtDPP4f_ileum_lumen) <- -kon*AtDPP4f_ileum_lumen*Aileum_lumen/Vileum_lumen
  +koff*AtDPP4c_ileum_lumen
  
  
  d/dt(ABlood) <- -kon*AsDPP4f_blood*ABlood/VBlood/RB*fupAlb
  +koff*AsDPP4c_blood
  
  -kon*AtDPP4f_blood*ABlood/VBlood/RB*fupAlb
  +koff*AtDPP4c_blood
  
  +Qh*Aliver_5/Vh_i/Kph*RB
  +Qrep*Areproductiveorgans/Vrep/Kprep*RB
  +Qthy*Athymus/Vthy/Kpthy*RB
  +Qsal*Asalivarygland/Vsal/Kpsal_vivo*RB
  +Qs*Askin/Vs/Kps*RB
  +Qrem*Aremainder/Vrem/Kprem*RB
  +(Qk-Qu)*Akves/Vkves/Kpk*RB
  -Ql*ABlood/VBlood
  
  d/dt(AsDPP4c_blood) <- +kon*AsDPP4f_blood*ABlood/VBlood/RB*fupAlb
  -koff*AsDPP4c_blood
  
  d/dt(AsDPP4f_blood) <- -kon*AsDPP4f_blood*ABlood/VBlood/RB*fupAlb
  +koff*AsDPP4c_blood
  
  d/dt(AtDPP4c_blood) <- +kon*AtDPP4f_blood*ABlood/VBlood/RB*fupAlb
  -koff*AtDPP4c_blood
  
  d/dt(AtDPP4f_blood) <- -kon*AtDPP4f_blood*ABlood/VBlood/RB*fupAlb
  +koff*AtDPP4c_blood
  
  
  d/dt(Alung) <- -kon*AsDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AsDPP4c_lung
  
  -kon*AtDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AtDPP4c_lung
  
  +Ql*ABlood/VBlood
  
  -Qserosa*Alung/Vl/Kpl*RB
  -Qduode_muc*Alung/Vl/Kpl*RB
  -Qjeju_muc*Alung/Vl/Kpl*RB
  -Qile_muc*Alung/Vl/Kpl*RB
  -Qha_proper*Alung/Vl/Kpl*RB
  -Qpan*Alung/Vl/Kpl*RB
  -Qspl*Alung/Vl/Kpl*RB
  -Qrep*Alung/Vl/Kpl*RB
  -Qthy*Alung/Vl/Kpl*RB
  -Qsal*Alung/Vl/Kpl*RB
  -Qs*Alung/Vl/Kpl*RB
  -Qrem*Alung/Vl/Kpl*RB
  -(Qk-GFR)*Alung/Vl/Kpl*RB
  -GFR*Alung/Vl/Kpl*fupAlb
  
  d/dt(AsDPP4c_lung) <- +kon*AsDPP4f_lung*Alung/Vl/Kpl*fupAlb
  -koff*AsDPP4c_lung
  
  d/dt(AsDPP4f_lung) <- -kon*AsDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AsDPP4c_lung
  
  d/dt(AtDPP4c_lung) <- +kon*AtDPP4f_lung*Alung/Vl/Kpl*fupAlb
  -koff*AtDPP4c_lung
  
  d/dt(AtDPP4f_lung) <- -kon*AtDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AtDPP4c_lung
  
  
  d/dt(Aliver_1) <- -kon*AsDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_1
  
  -kon*AtDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_1
  
  +Qserosa*Aserosa/Vserosa/Kpgut*RB
  +Qduode_muc*Amucosal_blood1/Vmucosal_blood1
  +Qjeju_muc*Amucosal_blood2/Vmucosal_blood2
  +Qile_muc*Amucosal_blood3/Vmucosal_blood3
  
  +Qha_proper*Alung/Vl/Kpl*RB
  +Qpan*Apancreas/Vpan/Kppan*RB
  +Qspl*Aspleen/Vspl/Kpspl*RB
  
  -Qh*Aliver_1/Vh_i/Kph*RB
  -CLinth_i*Aliver_1/Vh_i/Kph*fupAlb
  
  
  d/dt(AsDPP4c_hepat_1) <- +kon*AsDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  -koff*AsDPP4c_hepat_1
  
  d/dt(AsDPP4f_hepat_1) <- -kon*AsDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_1
  
  d/dt(AtDPP4c_hepat_1) <- +kon*AtDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  -koff*AtDPP4c_hepat_1
  
  d/dt(AtDPP4f_hepat_1) <- -kon*AtDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_1
  
  
  d/dt(Aliver_2) <- -kon*AsDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_2
  
  -kon*AtDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_2
  
  +Qh*Aliver_1/Vh_i/Kph*RB
  -Qh*Aliver_2/Vh_i/Kph*RB
  -CLinth_i*Aliver_2/Vh_i/Kph*fupAlb
  
  d/dt(AsDPP4c_hepat_2) <- +kon*AsDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  -koff*AsDPP4c_hepat_2
  
  d/dt(AsDPP4f_hepat_2) <- -kon*AsDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_2
  
  d/dt(AtDPP4c_hepat_2) <- +kon*AtDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  -koff*AtDPP4c_hepat_2
  
  d/dt(AtDPP4f_hepat_2) <- -kon*AtDPP4f_hepat_2*Aliver_2/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_2
  
  
  d/dt(Aliver_3) <- -kon*AsDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_3
  
  -kon*AtDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_3
  
  +Qh*Aliver_2/Vh_i/Kph*RB
  -Qh*Aliver_3/Vh_i/Kph*RB
  -CLinth_i*Aliver_3/Vh_i/Kph*fupAlb
  
  d/dt(AsDPP4c_hepat_3) <- +kon*AsDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  -koff*AsDPP4c_hepat_3
  
  d/dt(AsDPP4f_hepat_3) <- -kon*AsDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_3
  
  d/dt(AtDPP4c_hepat_3) <- +kon*AtDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  -koff*AtDPP4c_hepat_3
  
  d/dt(AtDPP4f_hepat_3) <- -kon*AtDPP4f_hepat_3*Aliver_3/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_3
  
  
  d/dt(Aliver_4) <- -kon*AsDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_4
  
  -kon*AtDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_4
  
  +Qh*Aliver_3/Vh_i/Kph*RB
  -Qh*Aliver_4/Vh_i/Kph*RB
  -CLinth_i*Aliver_4/Vh_i/Kph*fupAlb
  
  d/dt(AsDPP4c_hepat_4) <- +kon*AsDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  -koff*AsDPP4c_hepat_4
  
  d/dt(AsDPP4f_hepat_4) <- -kon*AsDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_4
  
  d/dt(AtDPP4c_hepat_4) <- +kon*AtDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  -koff*AtDPP4c_hepat_4
  
  d/dt(AtDPP4f_hepat_4) <- -kon*AtDPP4f_hepat_4*Aliver_4/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_4
  
  
  d/dt(Aliver_5) <- -kon*AsDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_5
  
  -kon*AtDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_5
  
  +Qh*Aliver_4/Vh_i/Kph*RB
  -Qh*Aliver_5/Vh_i/Kph*RB
  -CLinth_i*Aliver_5/Vh_i/Kph*fupAlb
  
  d/dt(AsDPP4c_hepat_5) <- +kon*AsDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  -koff*AsDPP4c_hepat_5
  
  d/dt(AsDPP4f_hepat_5) <- -kon*AsDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_5
  
  d/dt(AtDPP4c_hepat_5) <- +kon*AtDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  -koff*AtDPP4c_hepat_5
  
  d/dt(AtDPP4f_hepat_5) <- -kon*AtDPP4f_hepat_5*Aliver_5/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_5
  
  # Reproductive organs
  d/dt(Areproductiveorgans) <- -kon*AsDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  +koff*AsDPP4c_reproductiveorgans
  
  -kon*AtDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  +koff*AtDPP4c_reproductiveorgans
  
  -Qrep*Areproductiveorgans/Vrep/Kprep*RB
  +Qrep*Alung/Vl/Kpl*RB
  
  
  d/dt(AsDPP4c_reproductiveorgans) <- +kon*AsDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  -koff*AsDPP4c_reproductiveorgans
  
  d/dt(AsDPP4f_reproductiveorgans) <- -kon*AsDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  +koff*AsDPP4c_reproductiveorgans
  
  d/dt(AtDPP4c_reproductiveorgans) <- +kon*AtDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  -koff*AtDPP4c_reproductiveorgans
  
  d/dt(AtDPP4f_reproductiveorgans) <- -kon*AtDPP4f_reproductiveorgans*Areproductiveorgans/Vrep/Kprep*fupAlb
  +koff*AtDPP4c_reproductiveorgans
  
  # Thymus
  d/dt(Athymus) <- -kon*AsDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  +koff*AsDPP4c_thymus
  
  -kon*AtDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  +koff*AtDPP4c_thymus
  
  -Qthy*Athymus/Vthy/Kpthy*RB
  +Qthy*Alung/Vl/Kpl*RB
  
  
  d/dt(AsDPP4c_thymus) <- +kon*AsDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  -koff*AsDPP4c_thymus
  
  d/dt(AsDPP4f_thymus) <- -kon*AsDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  +koff*AsDPP4c_thymus
  
  d/dt(AtDPP4c_thymus) <- +kon*AtDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  -koff*AtDPP4c_thymus
  
  d/dt(AtDPP4f_thymus) <- -kon*AtDPP4f_thymus*Athymus/Vthy/Kpthy*fupAlb
  +koff*AtDPP4c_thymus
  
  # Salivary gland
  d/dt(Asalivarygland) <- -kon*AsDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  +koff*AsDPP4c_salivarygland
  
  -kon*AtDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  +koff*AtDPP4c_salivarygland
  
  -Qsal*Asalivarygland/Vsal/Kpsal_vivo*RB
  +Qsal*Alung/Vl/Kpl*RB
  
  
  d/dt(AsDPP4c_salivarygland) <- +kon*AsDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  -koff*AsDPP4c_salivarygland
  
  d/dt(AsDPP4f_salivarygland) <- -kon*AsDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  +koff*AsDPP4c_salivarygland
  
  d/dt(AtDPP4c_salivarygland) <- +kon*AtDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  -koff*AtDPP4c_salivarygland
  
  d/dt(AtDPP4f_salivarygland) <- -kon*AtDPP4f_salivarygland*Asalivarygland/Vsal/Kpsal_vivo*fupAlb
  +koff*AtDPP4c_salivarygland
  
  # Pancreas
  d/dt(Apancreas) <- -kon*AsDPP4f_pancreas*Apancreas/Vpan/Kppan*fupAlb
  +koff*AsDPP4c_pancreas
  
  -Qpan*Apancreas/Vpan/Kppan*RB
  +Qpan*Alung/Vl/Kpl*RB
  
  d/dt(AsDPP4c_pancreas) <- +kon*AsDPP4f_pancreas*Apancreas/Vpan/Kppan*fupAlb
  -koff*AsDPP4c_pancreas
  
  d/dt(AsDPP4f_pancreas) <- -kon*AsDPP4f_pancreas*Apancreas/Vpan/Kppan*fupAlb
  +koff*AsDPP4c_pancreas
  
  d/dt(Aspleen) <- -kon*AsDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AsDPP4c_spleen
  
  -kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AtDPP4c_spleen
  
  -Qspl*Aspleen/Vspl/Kpspl*RB
  +Qspl*Alung/Vl/Kpl*RB
  
  
  d/dt(AsDPP4c_spleen) <- +kon*AsDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  -koff*AsDPP4c_spleen
  
  d/dt(AsDPP4f_spleen) <- -kon*AsDPP4f_spleen *Aspleen/Vspl/Kpspl*fupAlb
  +koff*AsDPP4c_spleen
  
  d/dt(AtDPP4c_spleen) <- +kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  -koff*AtDPP4c_spleen
  
  d/dt(AtDPP4f_spleen) <- -kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AtDPP4c_spleen
  
  
  d/dt(Askin) <- -kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AsDPP4c_skin
  
  -kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AtDPP4c_skin
  
  -Qs*Askin/Vs/Kps*RB
  +Qs*Alung/Vl/Kpl*RB
  
  d/dt(AsDPP4c_skin) <- +kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  -koff*AsDPP4c_skin
  
  d/dt(AsDPP4f_skin) <- -kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AsDPP4c_skin
  
  d/dt(AtDPP4c_skin) <- +kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  -koff*AtDPP4c_skin
  
  d/dt(AtDPP4f_skin) <- -kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AtDPP4c_skin
  
  
  d/dt(Aremainder) <- -kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  +koff*AsDPP4c_remainder
  
  -Qrem*Aremainder/Vrem/Kprem*RB
  +Qrem*Alung/Vl/Kpl*RB
  
  d/dt(AsDPP4c_remainder) <- +kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  -koff*AsDPP4c_remainder
  
  d/dt(AsDPP4f_remainder) <- -kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  +koff*AsDPP4c_remainder
  
  
  d/dt(Akves) <- -kon*AsDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  +koff*AsDPP4c_kidneyVessel
  
  -kon*AtDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  +koff*AtDPP4c_kidneyVessel
  
  +(Qk-GFR)*Alung/Vl/Kpl*RB
  -(Qk-Qu)*Akves/Vkves/Kpk*RB
  
  -CLsecint*Akves/Vkves/Kpk*fupAlb
  +CLreabDPP4*AtDPP4c_kidneytubule/Vktubule
  
  
  d/dt(AsDPP4c_kidneyVessel) <- +kon*AsDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  -koff*AsDPP4c_kidneyVessel
  
  d/dt(AsDPP4f_kidneyVessel) <- -kon*AsDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  +koff*AsDPP4c_kidneyVessel
  
  d/dt(AtDPP4c_kidneyVessel) <- +kon*AtDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  -koff*AtDPP4c_kidneyVessel
  
  d/dt(AtDPP4f_kidneyVessel) <- -kon*AtDPP4f_kidneyVessel*Akves/Vkves/Kpk*fupAlb
  +koff*AtDPP4c_kidneyVessel
  
  
  d/dt(Aktubule) <- -kon*AtDPP4f_kidneytubule*Aktubule/Vktubule
  +koff*AtDPP4c_kidneytubule
  
  +GFR*Alung/Vl/Kpl*fupAlb # "GFR" is not a mistake
  +CLsecint*Akves/Vkves/Kpk*fupAlb
  -Qu*Aktubule/Vktubule
  
  d/dt(AtDPP4c_kidneytubule) <- +kon*AtDPP4f_kidneytubule*Aktubule/Vktubule
  -koff*AtDPP4c_kidneytubule
  -CLreabDPP4*AtDPP4c_kidneytubule/Vktubule
  
  d/dt(AtDPP4f_kidneytubule) <- -kon*AtDPP4f_kidneytubule*Aktubule/Vktubule
  +koff*AtDPP4c_kidneytubule
  +CLreabDPP4*AtDPP4c_kidneytubule/Vktubule
  
  
  d/dt(Xmetab) <- +CLinth_i*(1-fbile)*Aliver_1/Vh_i/Kph*fupAlb
  +CLinth_i*(1-fbile)*Aliver_2/Vh_i/Kph*fupAlb
  +CLinth_i*(1-fbile)*Aliver_3/Vh_i/Kph*fupAlb
  +CLinth_i*(1-fbile)*Aliver_4/Vh_i/Kph*fupAlb
  +CLinth_i*(1-fbile)*Aliver_5/Vh_i/Kph*fupAlb
  
  d/dt(Xurine) <- +Qu*Aktubule/Vktubule
  
  d/dt(Xfeces) <- +Alpha*kfeces_ileum*Aileum_lumen
  
  #################
  # soluble DPP-4 #
  #################
  
  Vtotal <- VBlood+Veh+Vsv+Vlv+Vkvesv+Vsplv+Vsalv+Vthyv+Vrepv+Vpanv+Vremv
  +Vmucosal_blood1+Vmucosal_blood2+Vmucosal_blood3+Vserosav
  
  AsDPP4f_blood(0) <- AsDPP4*VBlood/Vtotal
  AsDPP4f_hepat_1(0) <- AsDPP4*Veh/Vtotal/5
  AsDPP4f_hepat_2(0) <- AsDPP4*Veh/Vtotal/5
  AsDPP4f_hepat_3(0) <- AsDPP4*Veh/Vtotal/5
  AsDPP4f_hepat_4(0) <- AsDPP4*Veh/Vtotal/5
  AsDPP4f_hepat_5(0) <- AsDPP4*Veh/Vtotal/5
  AsDPP4f_skin(0) <- AsDPP4*Vsv/Vtotal
  AsDPP4f_lung(0) <- AsDPP4*Vlv/Vtotal
  AsDPP4f_kidneyVessel(0) <- AsDPP4*Vkvesv/Vtotal
  AsDPP4f_spleen(0) <- AsDPP4*Vsplv/Vtotal
  AsDPP4f_salivarygland(0) <- AsDPP4*Vsalv/Vtotal
  AsDPP4f_thymus(0) <- AsDPP4*Vthyv/Vtotal
  AsDPP4f_reproductiveorgans(0) <- AsDPP4*Vrepv/Vtotal
  AsDPP4f_pancreas(0) <- AsDPP4*Vpanv/Vtotal
  AsDPP4f_remainder(0) <- AsDPP4*Vremv/Vtotal
  AsDPP4f_mucosal_blood1(0) <- AsDPP4*Vmucosal_blood1/Vtotal
  AsDPP4f_mucosal_blood2(0) <- AsDPP4*Vmucosal_blood2/Vtotal
  AsDPP4f_mucosal_blood3(0) <- AsDPP4*Vmucosal_blood3/Vtotal
  AsDPP4f_serosa(0) <- AsDPP4*Vserosav/Vtotal
  
  ################
  # tissue DPP-4 #
  ################
  
  # Corrected for volume and converted back to DPP-4 expression per organ(nMol/organ).
  OrganDPP4_skin <- CDPP4_skin*Vs*(1-fVsv)
  OrganDPP4_lung <- CDPP4_lung*Vl*(1-fVlv)
  OrganDPP4_liver <- CDPP4_liver*(Vh-Veh)
  OrganDPP4_kidney <- CDPP4_kidney*(Vkves*(1-fVkvesv)+Vktubule)
  OrganDPP4_spleen <- CDPP4_spleen*Vspl*(1-fVsplv)
  OrganDPP4_salivarygland <- CDPP4_salivarygland*Vsal*(1-fVsalv)
  OrganDPP4_thymus <- CDPP4_thymus*Vthy*(1-fVthyv)
  OrganDPP4_reproductiveorgans <- CDPP4_reproductiveorgans*Vrep*(1-fVrepv)
  OrganDPP4_intestine <- OrganDPP4_liver*rGL
  
  # Convert to fraction
  ADPP4_total <- +OrganDPP4_skin+OrganDPP4_lung+OrganDPP4_liver+OrganDPP4_kidney
  +OrganDPP4_spleen+OrganDPP4_salivarygland+OrganDPP4_thymus
  +OrganDPP4_reproductiveorgans+OrganDPP4_intestine
  FAtDPP4_skin <- OrganDPP4_skin/ADPP4_total
  FAtDPP4_lung <- OrganDPP4_lung/ADPP4_total
  FAtDPP4_liver <- OrganDPP4_liver/ADPP4_total
  FAtDPP4_kidney <- OrganDPP4_kidney/ADPP4_total
  FAtDPP4_spleen <- OrganDPP4_spleen/ADPP4_total
  FAtDPP4_salivarygland <- OrganDPP4_salivarygland/ADPP4_total
  FAtDPP4_thymus <- OrganDPP4_thymus/ADPP4_total
  FAtDPP4_reproductiveorgans <- OrganDPP4_reproductiveorgans/ADPP4_total
  FAtDPP4_intestine <- OrganDPP4_intestine/ADPP4_total
  
  
  AtDPP4f_skin(0) <- AtDPP4*FAtDPP4_skin
  AtDPP4f_lung(0) <- AtDPP4*FAtDPP4_lung
  AtDPP4f_blood(0) <- 0
  AtDPP4f_hepat_1(0) <- AtDPP4*FAtDPP4_liver/5
  AtDPP4f_hepat_2(0) <- AtDPP4*FAtDPP4_liver/5
  AtDPP4f_hepat_3(0) <- AtDPP4*FAtDPP4_liver/5
  AtDPP4f_hepat_4(0) <- AtDPP4*FAtDPP4_liver/5
  AtDPP4f_hepat_5(0) <- AtDPP4*FAtDPP4_liver/5
  AtDPP4f_kidneyVessel(0) <- AtDPP4*FAtDPP4_kidney*(1-flumen_kidney)
  AtDPP4f_kidneytubule(0) <- AtDPP4*FAtDPP4_kidney*(flumen_kidney)
  AtDPP4f_spleen(0) <- AtDPP4*FAtDPP4_spleen
  AtDPP4f_salivarygland(0) <- AtDPP4*FAtDPP4_salivarygland
  AtDPP4f_thymus(0) <- AtDPP4*FAtDPP4_thymus
  AtDPP4f_reproductiveorgans(0) <- AtDPP4*FAtDPP4_reproductiveorgans
  
  AtDPP4_intestine <- AtDPP4*FAtDPP4_intestine
  Vmucosero <- Vmucosal_blood1*75.1+Vmucosal_blood2*105+Vmucosal_blood3*127+Vserosa*102.3667
  Vgutlumen <- Vduodenum_lumen*75.1+Vjejunum_lumen*105+Vileum_lumen*127
  
  AtDPP4f_duodenum_lumen(0) <- AtDPP4_intestine*flumen_gut*Vduodenum_lumen/Vgutlumen*75.1
  AtDPP4f_jejunum_lumen(0) <- AtDPP4_intestine*flumen_gut*Vjejunum_lumen/Vgutlumen*105
  AtDPP4f_ileum_lumen(0) <- AtDPP4_intestine*flumen_gut*Vileum_lumen/Vgutlumen*127
  
  AtDPP4f_mucosal_blood1(0) <- AtDPP4_intestine*(1-flumen_gut)*Vmucosal_blood1/Vmucosero*75.1
  AtDPP4f_mucosal_blood2(0) <- AtDPP4_intestine*(1-flumen_gut)*Vmucosal_blood2/Vmucosero*105
  AtDPP4f_mucosal_blood3(0) <- AtDPP4_intestine*(1-flumen_gut)*Vmucosal_blood3/Vmucosero*127
  AtDPP4f_serosa(0) <- AtDPP4_intestine*(1-flumen_gut)*Vserosa/Vmucosero*102.3667

    # Observation: plasma linagliptin concentration (nmol/L), total of
    # unbound drug and the sDPP-4-bound complex in plasma (deposited
    # Cplasma = ABlood/VBlood/RB + AsDPP4c_blood/Vplasma).
    Cc <- Cplasma
    Cc ~ prop(propSd)
  })
}
