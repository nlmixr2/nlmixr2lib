Nakamura_2024_linagliptin_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat). PBPK-TMDD (whole-body, single-segment",
    "segregated-flow intestine, 5-compartment tandem liver). Linagliptin",
    "in wild-type and Dpp-4-deficient Fischer/Sprague-Dawley rats after",
    "intravenous (0.01-50 mg/kg) and oral (2 mg/kg/day) dosing, fitted to",
    "plasma concentration and urinary/faecal excretion with the cluster",
    "Gauss-Newton method (CGNM). The nonlinear PK arises from reversible",
    "time-dependent binding (kon = koff/Kd) of linagliptin to soluble",
    "(AsDPP4) and transmembrane (AtDPP4) DPP-4, distributed over the",
    "organs by the tissue DPP-4 fractions estimated from rat",
    "autoradiography, plus a saturable DPP-4-complex-mediated renal",
    "reabsorption clearance (CLreabDPP4). The 16 parameters without",
    "fixed() are the CGNM medians of the plasma-concentration column of",
    "Table 1; every physiological constant and DPP-4 organ fraction is",
    "fixed() from Table S1 and the deposited code (Data S2). The",
    "Dpp-4-deficient rat is reproduced by setting the total DPP-4 amount",
    "ADPP4 to 0 (see the vignette). CGNM is a fixed-effects method, so the",
    "model has no between-subject variability and no residual-error model",
    "-- the propSd term is a placeholder. See the vignette Errata."
  )
  reference <- paste(
    "Nakamura R, Yoshikado T, Aoki Y, Sugiyama Y, Chiba K.",
    "Elucidation of DPP-4 involvement in systemic distribution and renal",
    "reabsorption of linagliptin by PBPK modeling with a cluster",
    "Gauss-Newton method. Clin Transl Sci. 2024;17(10):e70047.",
    "doi:10.1111/cts.70047.",
    "The ODE system is transcribed from the deposited rat plasma model",
    "code (Data S2, Program/pbpk_lng_rat_plasma/); estimated parameters",
    "are the plasma-concentration column of Table 1 and fixed parameters",
    "are Table S1 and the deposited paramdata.",
    sep = " "
  )
  vignette <- "Nakamura_2024_linagliptin_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  # PBPK organ-amount and DPP-4 (free / complex) states, kept with the
  # deposited Data S2 names (uppercase A-prefixed amounts) so the ODE
  # block can be verified line-by-line against the authors' code. They
  # are declared here so checkModelConventions() does not flag them as
  # non-canonical compartment names.
  paper_specific_compartments <- c(
    "Aehc",
    "Aserosa",
    "AsDPP4c_serosa",
    "AsDPP4f_serosa",
    "AtDPP4c_serosa",
    "AtDPP4f_serosa",
    "Amucosal_blood",
    "AsDPP4c_mucosal_blood",
    "AsDPP4f_mucosal_blood",
    "AtDPP4c_mucosal_blood",
    "AtDPP4f_mucosal_blood",
    "Aenterocyte",
    "Aintestinal_lumen",
    "AtDPP4c_intestinal_lumen",
    "AtDPP4f_intestinal_lumen",
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
    Aehc = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
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
    Amucosal_blood = list(analyte = "linagliptin", units = "nmol", specimen = "whole blood", verified = FALSE),
    AsDPP4c_mucosal_blood = list(
      analyte = "linagliptin-sDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AsDPP4f_mucosal_blood = list(
      analyte = "soluble DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4c_mucosal_blood = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    AtDPP4f_mucosal_blood = list(
      analyte = "transmembrane DPP-4 (free)",
      units = "nmol",
      specimen = "whole blood",
      verified = FALSE
    ),
    Aenterocyte = list(analyte = "linagliptin", units = "nmol", specimen = "tissue", verified = FALSE),
    Aintestinal_lumen = list(
      analyte = "linagliptin",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    AtDPP4c_intestinal_lumen = list(
      analyte = "linagliptin-tDPP-4 complex",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    AtDPP4f_intestinal_lumen = list(
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
    species = "rat (wild-type Fischer/Sprague-Dawley and Dpp-4-deficient Fischer)",
    disease_state = "healthy; wild-type vs Dpp-4-deficient",
    dose_range = "intravenous 0.01-50 mg/kg; oral 2 mg/kg/day",
    notes = paste(
      "Plasma concentration-time data (wild-type and Dpp-4-deficient rats,",
      "0.01-50 mg/kg i.v.) are from Retlich 2009; urinary and faecal",
      "excretion (wild-type, 2 mg/kg/day p.o.) and the autoradiography",
      "tissue data that fix the DPP-4 organ fractions are from Fuchs 2009.",
      "Physiological parameters are scaled to a 0.25 kg rat. Doses are",
      "supplied in nmol (linagliptin MW 472.54; the deposited code uses",
      "2116.22 nmol/mg, so a 2 mg/kg dose to a 0.25 kg rat is",
      "2116.22 * 2 * 0.25 = 1058 nmol)."
    )
  )

  ini({
    lVBlood <- fixed(log(0.0178571428571429)); label("Blood volume (L)")  # Table S1 rat physiological parameter VBlood = 0.01785714
    lVeh <- fixed(log(0.00328)); label("Extra-Hepatic (Disse/Interstitial) volume (L)")  # Table S1 rat physiological parameter Veh = 0.00328
    lVenterocyte <- fixed(log(0.0002083325)); label("Enterocyte volume (L)")  # Table S1 rat physiological parameter Venterocyte = 0.0002083325
    lVh <- fixed(log(0.0118)); label("Liver volume (L)")  # Table S1 rat physiological parameter Vh = 0.0118
    lVintestinal_lumen <- fixed(log(0.00138889)); label("Intestinal Lumen volume (L)")  # Table S1 rat physiological parameter Vintestinal_lumen = 0.00138889
    lVktubule <- fixed(log(0.000907142857142857)); label("Kidney Tubule (Duct) volume (L)")  # Table S1 rat physiological parameter Vktubule = 0.0009071429
    lVkves <- fixed(log(0.00109285714285714)); label("Kidney Vessel+cell volume (L)")  # Table S1 rat physiological parameter Vkves = 0.001092857
    lVl <- fixed(log(0.00155)); label("Lung volume (L)")  # Table S1 rat physiological parameter Vl = 0.00155
    lVmucosal_blood <- fixed(log(0.0001125)); label("Mucosal Blood volume (L)")  # Table S1 rat physiological parameter Vmucosal_blood = 0.0001125
    lVpan <- fixed(log(0.00125)); label("Pancreas volume (L)")  # Table S1 rat physiological parameter Vpan = 0.00125
    lVplasma <- fixed(log(0.01)); label("Plasma volume (L)")  # Table S1 rat physiological parameter Vplasma = 0.01
    lVrem <- fixed(log(0.17908)); label("Remainder volume (L)")  # Table S1 rat physiological parameter Vrem = 0.17908
    lVrep <- fixed(log(0.00288)); label("Reproductive Organs volume (L)")  # Table S1 rat physiological parameter Vrep = 0.00288
    lVs <- fixed(log(0.04)); label("Skin volume (L)")  # Table S1 rat physiological parameter Vs = 0.04
    lVsal <- fixed(log(0.00043)); label("Salivary Gland volume (L)")  # Table S1 rat physiological parameter Vsal = 0.00043
    lVserosa <- fixed(log(0.0028875)); label("Serosa volume (L)")  # Table S1 rat physiological parameter Vserosa = 0.0028875
    lVspl <- fixed(log(0.00075)); label("Spleen volume (L)")  # Table S1 rat physiological parameter Vspl = 0.00075
    lVthy <- fixed(log(0.00025)); label("Thymus volume (L)")  # Table S1 rat physiological parameter Vthy = 0.00025
    fVkvesv <- fixed(0.105); label("Kidney fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVkvesv = 0.105
    fVlv <- fixed(0.262); label("Lung fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVlv = 0.262
    fVpanv <- fixed(0.1); label("Pancreas fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVpanv = 0.1
    fVremv <- fixed(0.1); label("Remainder fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVremv = 0.1
    fVrepv <- fixed(0.02); label("Reproductive Organs fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVrepv = 0.02
    fVsalv <- fixed(0.1); label("Salivary Gland fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVsalv = 0.1
    fVserosav <- fixed(0.1); label("Serosa fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVserosav = 0.1
    fVsplv <- fixed(0.282); label("Spleen fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVsplv = 0.282
    fVsv <- fixed(0.019); label("Skin fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVsv = 0.019
    fVthyv <- fixed(0.03); label("Thymus fractional capillary (vessel) volume (unitless)")  # Table S1 rat physiological parameter fVthyv = 0.03
    lGFR <- fixed(log(0.0785)); label("Glomerular filtration rate (L/h)")  # Table S1 rat physiological parameter GFR = 0.0785
    lQh <- fixed(log(0.612)); label("Liver blood flow (L/h)")  # Table S1 rat physiological parameter Qh = 0.612
    lQha_proper <- fixed(log(0.2125)); label("Proper Hepatic Artery blood flow (L/h)")  # Table S1 rat physiological parameter Qha_proper = 0.2125
    lQk <- fixed(log(0.398)); label("Kidney blood flow (L/h)")  # Table S1 rat physiological parameter Qk = 0.398
    lQl <- fixed(log(3.06)); label("Lung blood flow (L/h)")  # Table S1 rat physiological parameter Ql = 3.06
    lQmucosal <- fixed(log(0.03325)); label("Mucosal blood flow (L/h)")  # Table S1 rat physiological parameter Qmucosal = 0.03325
    lQpan <- fixed(log(0.0428)); label("Pancreas blood flow (L/h)")  # Table S1 rat physiological parameter Qpan = 0.0428
    lQrem <- fixed(log(1.63573)); label("Remainder blood flow (L/h)")  # Table S1 rat physiological parameter Qrem = 1.63573
    lQrep <- fixed(log(0.026)); label("Reproductive Organs blood flow (L/h)")  # Table S1 rat physiological parameter Qrep = 0.026
    lQs <- fixed(log(0.315)); label("Skin blood flow (L/h)")  # Table S1 rat physiological parameter Qs = 0.315
    lQsal <- fixed(log(0.052)); label("Salivary Gland blood flow (L/h)")  # Table S1 rat physiological parameter Qsal = 0.052
    lQserosa <- fixed(log(0.3)); label("Serosa blood flow (L/h)")  # Table S1 rat physiological parameter Qserosa = 0.3
    lQspl <- fixed(log(0.0245)); label("Spleen blood flow (L/h)")  # Table S1 rat physiological parameter Qspl = 0.0245
    lQthy <- fixed(log(0.0214)); label("Thymus blood flow (L/h)")  # Table S1 rat physiological parameter Qthy = 0.0214
    lQu <- fixed(log(0.00208333)); label("Urine blood flow (L/h)")  # Table S1 rat physiological parameter Qu = 0.00208333
    lKpa_calc <- fixed(log(0.201)); label("Adipose tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpa_calc = 0.201
    lKpgut_calc <- fixed(log(1.181)); label("Gut tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpgut_calc = 1.181
    lKph_calc <- fixed(log(1.915)); label("Liver tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kph_calc = 1.915
    lKphrt_calc <- fixed(log(1.00752594533395)); label("Heart tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kphrt_calc = 1.007526
    lKpk_calc <- fixed(log(2.058)); label("Kidney tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpk_calc = 2.058
    lKpl_calc <- fixed(log(1.59)); label("Lung tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpl_calc = 1.59
    lKpm_calc <- fixed(log(0.913)); label("Muscle tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpm_calc = 0.913
    lKppan_calc <- fixed(log(0.92286251315842)); label("Pancreas tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kppan_calc = 0.9228625
    lKps_calc <- fixed(log(0.726)); label("Skin tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kps_calc = 0.726
    lKpspl_calc <- fixed(log(1.36651167316662)); label("Spleen tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpspl_calc = 1.366512
    lKpthy_calc <- fixed(log(1.1888252721636)); label("Thymus tissue:plasma partition coefficient (unitless)")  # Table S1 rat physiological parameter Kpthy_calc = 1.188825
    fgut <- fixed(0.232); label("Free fraction in enterocyte (unitless)")  # Table S1 rat physiological parameter fgut = 0.232
    fupAlb <- fixed(0.225); label("Non-saturable plasma free fraction (albumin binding) (unitless)")  # Table S1 rat physiological parameter fupAlb = 0.225
    lKmPgp <- fixed(log(187000)); label("Michaelis constant for intestinal P-gp (nmol/L)")  # Table S1 rat physiological parameter KmPgp = 187000
    RB <- fixed(0.691); label("Blood-to-plasma concentration ratio (unitless)")  # Table S1 rat physiological parameter RB = 0.691
    Alpha <- fixed(1); label("Gastrointestinal-motility transit factor (unitless)")  # Data S2 rat Alpha=1 (gastrointestinal-motility factor; not used in rat model)
    AR <- fixed(1); label("Apical reabsorption-ratio switch (unitless)")  # Data S2 AR=1 (apical reabsorption-ratio switch, off)
    Gammaent <- fixed(1); label("Enterocyte basolateral:apical diffusion ratio (unitless)")  # Data S2 Gammaent=1
    lrGL <- fixed(log(0.291)); label("Intestine-to-liver DPP-4 amount ratio (unitless)")  # Supplementary Text: rGL 0.291, DPP-4 gut:liver amount ratio from the Darmoul 1994 activity ratio
    lkfeces <- fixed(log(0.637)); label("Intestinal luminal transit rate constant (1/h)")  # Data S2 rat; intestinal luminal transit rate constant (1/h)
    fsDPP4 <- fixed(0.001975627); label("Fraction of total DPP-4 that is soluble (unitless)")  # Data S2 rat plasma; fraction of total DPP-4 that is soluble (from autoradiography analysis)
    ftDPP4 <- fixed(0.998024373); label("Fraction of total DPP-4 that is transmembrane (unitless)")  # Data S2 rat plasma; fraction of total DPP-4 that is transmembrane (1 - fsDPP4)
    fDPP4_skin <- fixed(0.04019386); label("Skin fraction of transmembrane DPP-4 (unitless)")  # Data S2; skin fraction of transmembrane DPP-4 (autoradiography analysis, Figure 1e)
    fDPP4_lung <- fixed(0.016447626); label("Lung fraction of transmembrane DPP-4 (unitless)")  # Data S2; lung fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_liver <- fixed(0.183671181); label("Liver fraction of transmembrane DPP-4 (unitless)")  # Data S2; liver fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_kidney <- fixed(0.682043533); label("Kidney fraction of transmembrane DPP-4 (unitless)")  # Data S2; kidney fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_spleen <- fixed(0.004560798); label("Spleen fraction of transmembrane DPP-4 (unitless)")  # Data S2; spleen fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_saliva <- fixed(0.001410761); label("Salivary Gland fraction of transmembrane DPP-4 (unitless)")  # Data S2; salivary-gland fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_thymus <- fixed(0.001276369); label("Thymus fraction of transmembrane DPP-4 (unitless)")  # Data S2; thymus fraction of transmembrane DPP-4 (Figure 1e)
    fDPP4_reprod <- fixed(0.01497193); label("Reproductive Organs fraction of transmembrane DPP-4 (unitless)")  # Data S2; reproductive-organ fraction of transmembrane DPP-4 (Figure 1e)
    lKph_vivo <- fixed(log(12.9439)); label("Liver tissue:plasma partition coefficient (unitless)")  # Data S2 rat plasma; Kp liver 12.94390 carried from rat autoradiography analysis
    lCLinth <- log(5.16); label("Hepatic intrinsic clearance (L/h)")  # Table 1 (plasma) CLint,h median 5.16 L/h
    lCLsecint <- log(0.695); label("Renal tubular secretion clearance (L/h)")  # Table 1 (plasma) CLsec,int median 0.695 L/h
    lCLreabDPP4 <- log(0.00000716); label("DPP-4-complex-mediated renal reabsorption clearance (L/h)")  # Table 1 (plasma) CLreab,DPP4 median 7.16e-6 L/h
    lADPP4 <- log(15.79); label("Total amount of DPP-4, soluble + transmembrane (nmol)")  # Table 1 (plasma) AsDPP4 0.0312 nmol / fsDPP4 0.001975627 = 15.79 nmol total DPP-4
    fbile <- 0.395; label("Fraction of hepatic intrinsic clearance that is biliary (unitless)")  # Table 1 (plasma) f bile median 0.395
    lkoff <- log(0.165); label("LNG-DPP-4 dissociation rate constant (1/h)")  # Table 1 (plasma) k off median 0.165 1/h
    lKd <- log(0.0238); label("LNG-DPP-4 dissociation constant (nmol/L)")  # Table 1 (plasma) K d median 0.0238 nM
    flumen_kidney <- 0.0419; label("Fraction of kidney DPP-4 on the luminal side (unitless)")  # Table 1 (plasma) f lumen,kidney median 0.0419
    flumen_gut <- 0.605; label("Fraction of intestinal DPP-4 on the luminal side (unitless)")  # Table 1 (plasma) f lumen,gut median 0.605
    lKprep_vivo <- log(2.31); label("Reproductive Organs tissue:plasma partition coefficient (unitless)")  # Table 1 (plasma) K p,rep.vivo median 2.31
    lKpsal_vivo <- log(1.76); label("Salivary Gland tissue:plasma partition coefficient (unitless)")  # Table 1 (plasma) K p,sal.vivo median 1.76
    lKprem_vivo <- log(3.13); label("Remainder tissue:plasma partition coefficient (unitless)")  # Table 1 (plasma) K p,rem.vivo median 3.13
    lkbile <- log(0.0623); label("Biliary transit rate constant (1/h)")  # Table 1 (plasma) k bile median 0.0623 1/h
    lVmaxPgp <- log(4170); label("Maximum intestinal P-gp transport rate (nmol/h)")  # Table 1 (plasma) V max,Pgp median 4.17e3 nmol/h
    lLR <- log(4.29); label("Ratio of lumen-to-enterocyte vs blood-to-enterocyte permeability (unitless)")  # Table 1 (plasma) LR median 4.29
    lPSdifentBE <- log(0.00141); label("Mucosal-blood-to-enterocyte permeability clearance (L/h)")  # Table 1 (plasma) PS dif,ent,BE median 0.00141 L/h
    propSd <- fixed(0.1); label("Proportional residual error (placeholder; not estimated)")  # placeholder: CGNM is fixed-effects least-squares; the paper reports no residual-error model
  })

  model({
    # Back-transform the log-scale parameters.
    VBlood <- exp(lVBlood)
    Veh <- exp(lVeh)
    Venterocyte <- exp(lVenterocyte)
    Vh <- exp(lVh)
    Vintestinal_lumen <- exp(lVintestinal_lumen)
    Vktubule <- exp(lVktubule)
    Vkves <- exp(lVkves)
    Vl <- exp(lVl)
    Vmucosal_blood <- exp(lVmucosal_blood)
    Vpan <- exp(lVpan)
    Vplasma <- exp(lVplasma)
    Vrem <- exp(lVrem)
    Vrep <- exp(lVrep)
    Vs <- exp(lVs)
    Vsal <- exp(lVsal)
    Vserosa <- exp(lVserosa)
    Vspl <- exp(lVspl)
    Vthy <- exp(lVthy)
    GFR <- exp(lGFR)
    Qh <- exp(lQh)
    Qha_proper <- exp(lQha_proper)
    Qk <- exp(lQk)
    Ql <- exp(lQl)
    Qmucosal <- exp(lQmucosal)
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
    KmPgp <- exp(lKmPgp)
    rGL <- exp(lrGL)
    kfeces <- exp(lkfeces)
    Kph_vivo <- exp(lKph_vivo)
    CLinth <- exp(lCLinth)
    CLsecint <- exp(lCLsecint)
    CLreabDPP4 <- exp(lCLreabDPP4)
    ADPP4 <- exp(lADPP4)
    koff <- exp(lkoff)
    Kd <- exp(lKd)
    Kprep_vivo <- exp(lKprep_vivo)
    Kpsal_vivo <- exp(lKpsal_vivo)
    Kprem_vivo <- exp(lKprem_vivo)
    kbile <- exp(lkbile)
    VmaxPgp <- exp(lVmaxPgp)
    LR <- exp(lLR)
    PSdifentBE <- exp(lPSdifentBE)

  
  # For Simulation
  
  Vh_i <- Vh/5          # Vh : Total volume of liver
  Veh_i <- Veh/5        # Veh: Extra hepatic (Disse + Interstitial) volume
  CLinth_i <- CLinth/5  # CLinth : Hepatic intrinsic clearance
  
  # Capillary Vessel Volume
  Vkvesv <- Vkves*fVkvesv       # Kidney
  Vsv <- Vs*fVsv                # Skin
  Vlv <- Vl*fVlv                # Lung
  Vpanv <- Vpan*fVpanv          # pancreas
  Vsplv <- Vspl*fVsplv          # spleen
  Vrepv <- Vrep*fVrepv          # reproductiveorgans
  Vthyv <- Vthy*fVthyv          # thymus
  Vsalv <- Vsal*fVsalv          # salivarygland
  Vremv <- Vrem*fVremv          # remainder
  Vserosav <- Vserosa*fVserosav # serosa
  
  # DPP-4 binding
  kon <- koff/Kd
  
  # Tissue to Plasma concentration ratio (Kp)
  SFKp <- 1
  
  Kpk <- Kpk_calc*SFKp   # Kidney
  Kphrt <- Kphrt_calc*SFKp # Heart
  Kpspl <- Kpspl_calc*SFKp # Spleen
  Kpm <- Kpm_calc*SFKp   # Muscle
  Kps <- Kps_calc*SFKp   # Skin
  Kpa <- Kpa_calc*SFKp   # Adipose
  Kpl <- Kpl_calc*SFKp   # Lung
  Kpthy <- Kpthy_calc*SFKp # Thymus
  Kppan <- Kppan_calc*SFKp # Pancreas
  Kpgut <- Kpgut_calc*SFKp # Intestine
  
  Kprep <- Kprep_vivo*SFKp # Reproductive Organ
  Kprem <- Kprem_vivo*SFKp # Remainder (Rest of Body)
  Kph <- Kph_vivo*SFKp   # Liver
  
  # Concentration for each organ
  Blood <- (ABlood+AsDPP4c_blood)/VBlood
  Plasma <- (ABlood/VBlood/RB)+(AsDPP4c_blood/Vplasma)
  Liver <- (Aliver_1+Aliver_2+Aliver_3+Aliver_4+Aliver_5+
             AsDPP4c_hepat_1+AsDPP4c_hepat_2+AsDPP4c_hepat_3+AsDPP4c_hepat_4+AsDPP4c_hepat_5+
             AtDPP4c_hepat_1+AtDPP4c_hepat_2+AtDPP4c_hepat_3+AtDPP4c_hepat_4+AtDPP4c_hepat_5
  )/Vh
  Pancreas <- (Apancreas+AsDPP4c_pancreas)/Vpan
  Spleen <- (Aspleen+AsDPP4c_spleen+AtDPP4c_spleen)/Vspl
  Lung <- (Alung+AsDPP4c_lung+AtDPP4c_lung)/Vl
  Reproductiveorgans <- (Areproductiveorgans+AsDPP4c_reproductiveorgans+AtDPP4c_reproductiveorgans)/Vrep
  Thymus <- (Athymus+AsDPP4c_thymus+AtDPP4c_thymus)/Vthy
  Skin <- (Askin+AsDPP4c_skin+AtDPP4c_skin)/Vs
  Salivarygland <- (Asalivarygland+AsDPP4c_salivarygland+AtDPP4c_salivarygland)/Vsal
  Kidney <- (Akves+AsDPP4c_kidneyVessel+AtDPP4c_kidneyVessel+Aktubule+AtDPP4c_kidneytubule)/(Vkves+Vktubule)
  
  #########################
  # Differential equation #
  #########################
  
  # EHC
  d/dt(Aehc) <- +CLinth_i*fbile*Aliver_1/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_2/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_3/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_4/Vh_i/Kph*fupAlb
  +CLinth_i*fbile*Aliver_5/Vh_i/Kph*fupAlb
  -kbile*Aehc
  
  # Serosa
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
  
  # Mucosal Blood
  d/dt(Amucosal_blood) <- -kon*AsDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood
  
  -kon*AtDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood
  
  +Qmucosal*Alung/Vl/Kpl*RB
  -Qmucosal*Amucosal_blood/Vmucosal_blood
  +fgut*(PSdifentBE/Gammaent)*Aenterocyte/Venterocyte
  -(fupAlb/RB)*PSdifentBE*Amucosal_blood/Vmucosal_blood
  
  d/dt(AsDPP4c_mucosal_blood) <- +kon*AsDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  -koff*AsDPP4c_mucosal_blood
  
  d/dt(AsDPP4f_mucosal_blood) <- -kon*AsDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  +koff*AsDPP4c_mucosal_blood
  
  d/dt(AtDPP4c_mucosal_blood) <- +kon*AtDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  -koff*AtDPP4c_mucosal_blood
  
  d/dt(AtDPP4f_mucosal_blood) <- -kon*AtDPP4f_mucosal_blood*Amucosal_blood/Vmucosal_blood/RB*fupAlb
  +koff*AtDPP4c_mucosal_blood
  
  # Enterocyte
  d/dt(Aenterocyte) <- -fgut*(PSdifentBE/Gammaent)*Aenterocyte/Venterocyte
  +(fupAlb/RB)*PSdifentBE*Amucosal_blood/Vmucosal_blood
  +LR*PSdifentBE*Aintestinal_lumen/Vintestinal_lumen
  -fgut*VmaxPgp/(KmPgp+fgut*Aenterocyte/Venterocyte)*Aenterocyte/Venterocyte
  -fgut*AR*(PSdifentBE/Gammaent)*Aenterocyte/Venterocyte
  
  # Intestinal Lumen
  d/dt(Aintestinal_lumen) <- -kon*AtDPP4f_intestinal_lumen*Aintestinal_lumen/Vintestinal_lumen
  +koff*AtDPP4c_intestinal_lumen
  
  +kbile*Aehc
  -Alpha*kfeces*Aintestinal_lumen
  -LR*PSdifentBE*Aintestinal_lumen/Vintestinal_lumen
  +fgut*VmaxPgp/(KmPgp+fgut*Aenterocyte/Venterocyte)*Aenterocyte/Venterocyte
  +fgut*AR*(PSdifentBE/Gammaent)*Aenterocyte/Venterocyte
  
  d/dt(AtDPP4c_intestinal_lumen) <- +kon*AtDPP4f_intestinal_lumen*Aintestinal_lumen/Vintestinal_lumen
  -koff*AtDPP4c_intestinal_lumen
  
  d/dt(AtDPP4f_intestinal_lumen) <- -kon*AtDPP4f_intestinal_lumen*Aintestinal_lumen/Vintestinal_lumen
  +koff*AtDPP4c_intestinal_lumen
  
  # Systemic Blood
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
  
  # Lung
  d/dt(Alung) <- -kon*AsDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AsDPP4c_lung
  
  -kon*AtDPP4f_lung*Alung/Vl/Kpl*fupAlb
  +koff*AtDPP4c_lung
  
  +Ql*ABlood/VBlood
  -Qserosa*Alung/Vl/Kpl*RB
  -Qmucosal*Alung/Vl/Kpl*RB
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
  
  # Liver 1
  d/dt(Aliver_1) <- -kon*AsDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AsDPP4c_hepat_1
  
  -kon*AtDPP4f_hepat_1*Aliver_1/Vh_i/Kph*fupAlb
  +koff*AtDPP4c_hepat_1
  
  +Qha_proper*Alung/Vl/Kpl*RB
  +Qpan*Apancreas/Vpan/Kppan*RB
  +Qspl*Aspleen/Vspl/Kpspl*RB
  +Qserosa*Aserosa/Vserosa/Kpgut*RB
  +Qmucosal*Amucosal_blood/Vmucosal_blood
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
  
  # Liver 2
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
  
  # Liver 3
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
  
  # Liver 4
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
  
  # Liver 5
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
  
  +Qrep*Alung/Vl/Kpl*RB
  -Qrep*Areproductiveorgans/Vrep/Kprep*RB
  
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
  
  +Qthy*Alung/Vl/Kpl*RB
  -Qthy*Athymus/Vthy/Kpthy*RB
  
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
  
  +Qsal*Alung/Vl/Kpl*RB
  -Qsal*Asalivarygland/Vsal/Kpsal_vivo*RB
  
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
  
  +Qpan*Alung/Vl/Kpl*RB
  -Qpan*Apancreas/Vpan/Kppan*RB
  
  d/dt(AsDPP4c_pancreas) <- +kon*AsDPP4f_pancreas*Apancreas/Vpan/Kppan*fupAlb
  -koff*AsDPP4c_pancreas
  
  d/dt(AsDPP4f_pancreas) <- -kon*AsDPP4f_pancreas*Apancreas/Vpan/Kppan*fupAlb
  +koff*AsDPP4c_pancreas
  
  
  d/dt(Aspleen) <- -kon*AsDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AsDPP4c_spleen
  
  -kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AtDPP4c_spleen
  
  +Qspl*Alung/Vl/Kpl*RB
  -Qspl*Aspleen/Vspl/Kpspl*RB
  
  d/dt(AsDPP4c_spleen) <- +kon*AsDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  -koff*AsDPP4c_spleen
  
  d/dt(AsDPP4f_spleen) <- -kon*AsDPP4f_spleen *Aspleen/Vspl/Kpspl*fupAlb
  +koff*AsDPP4c_spleen
  
  d/dt(AtDPP4c_spleen) <- +kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  -koff*AtDPP4c_spleen
  
  d/dt(AtDPP4f_spleen) <- -kon*AtDPP4f_spleen*Aspleen/Vspl/Kpspl*fupAlb
  +koff*AtDPP4c_spleen
  
  # Skin
  d/dt(Askin) <- -kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AsDPP4c_skin
  
  -kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AtDPP4c_skin
  
  +Qs*Alung/Vl/Kpl*RB
  -Qs*Askin/Vs/Kps*RB
  
  d/dt(AsDPP4c_skin) <- +kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  -koff*AsDPP4c_skin
  
  d/dt(AsDPP4f_skin) <- -kon*AsDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AsDPP4c_skin
  
  d/dt(AtDPP4c_skin) <- +kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  -koff*AtDPP4c_skin
  
  d/dt(AtDPP4f_skin) <- -kon*AtDPP4f_skin*Askin/Vs/Kps*fupAlb
  +koff*AtDPP4c_skin
  
  # Remainder
  d/dt(Aremainder) <- -kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  +koff*AsDPP4c_remainder
  
  +Qrem*Alung/Vl/Kpl*RB
  -Qrem*Aremainder/Vrem/Kprem*RB
  
  d/dt(AsDPP4c_remainder) <- +kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  -koff*AsDPP4c_remainder
  
  d/dt(AsDPP4f_remainder) <- -kon*AsDPP4f_remainder*Aremainder/Vrem/Kprem*fupAlb
  +koff*AsDPP4c_remainder
  
  # Kidney (Cell + Vessel)
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
  
  # Kidney (Tubule)
  d/dt(Aktubule) <- -kon*AtDPP4f_kidneytubule*Aktubule/Vktubule
  +koff*AtDPP4c_kidneytubule
  
  +GFR*Alung/Vl/Kpl*fupAlb # GFR is intentional, not a typo
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
  
  d/dt(Xfeces) <- +Alpha*kfeces*Aintestinal_lumen
  
  ######################################
  # Allocating Amount of Soluble DPP-4 #
  ######################################
  AsDPP4 <- ADPP4*fsDPP4
  Vtotal <- +VBlood+Veh+Vsv+Vlv+Vkvesv+Vsplv
  +Vsalv+Vthyv+Vrepv+Vpanv+Vremv
  +Vmucosal_blood+Vserosav
  
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
  AsDPP4f_mucosal_blood(0) <- AsDPP4*Vmucosal_blood/Vtotal
  AsDPP4f_serosa(0) <- AsDPP4*Vserosav/Vtotal
  
  ############################################
  # Allocating Amount of Transmembrane DPP-4 #
  ############################################
  AtDPP4 <- ADPP4*ftDPP4
  
  AtDPP4f_skin(0) <- AtDPP4*fDPP4_skin
  AtDPP4f_lung(0) <- AtDPP4*fDPP4_lung
  AtDPP4f_blood(0) <- 0
  AtDPP4f_hepat_1(0) <- AtDPP4*fDPP4_liver/5
  AtDPP4f_hepat_2(0) <- AtDPP4*fDPP4_liver/5
  AtDPP4f_hepat_3(0) <- AtDPP4*fDPP4_liver/5
  AtDPP4f_hepat_4(0) <- AtDPP4*fDPP4_liver/5
  AtDPP4f_hepat_5(0) <- AtDPP4*fDPP4_liver/5
  AtDPP4f_kidneyVessel(0) <- AtDPP4*fDPP4_kidney*(1-flumen_kidney)
  AtDPP4f_kidneytubule(0) <- AtDPP4*fDPP4_kidney*(flumen_kidney)
  AtDPP4f_spleen(0) <- AtDPP4*fDPP4_spleen
  AtDPP4f_salivarygland(0) <- AtDPP4*fDPP4_saliva
  AtDPP4f_thymus(0) <- AtDPP4*fDPP4_thymus
  AtDPP4f_reproductiveorgans(0) <- AtDPP4*fDPP4_reprod
  
  AtDPP4_intestine <- AtDPP4*fDPP4_liver*rGL
  AtDPP4f_intestinal_lumen(0) <- AtDPP4_intestine*(flumen_gut)
  AtDPP4f_mucosal_blood(0) <- AtDPP4_intestine*(1-flumen_gut)*Vmucosal_blood/(Vmucosal_blood+Vserosa)
  AtDPP4f_serosa(0) <- AtDPP4_intestine*(1-flumen_gut)*Vserosa/(Vmucosal_blood+Vserosa)

    # Observation: plasma linagliptin concentration (nmol/L), total of
    # unbound drug and sDPP-4-bound complex (deposited Plasma =
    # ABlood/VBlood/RB + AsDPP4c_blood/Vplasma).
    Cc <- Plasma
    Cc ~ prop(propSd)
  })
}
