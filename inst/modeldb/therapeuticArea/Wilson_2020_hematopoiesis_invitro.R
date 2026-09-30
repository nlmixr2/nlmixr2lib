Wilson_2020_hematopoiesis_invitro <- function() {
  description <- paste(
    "In vitro (human bone-marrow CD34+ cells, multilineage hematopoietic",
    "toxicity assay). QSP. Deterministic ODE model of in vitro",
    "hematopoiesis used to deconvolve drug-induced multilineage cytopenia",
    "into anti-proliferative and cell-killing mechanisms (Wilson 2020).",
    "Fifteen cell-count states: hematopoietic stem cells (HSC) renew",
    "themselves and feed multipotent progenitors (MPP), which branch into",
    "granulocyte-macrophage progenitors (GMP), megakaryocytes (MK), early",
    "erythroid cells and lymphoid progenitors; GMP branch into monocyte",
    "progenitors and granulocyte progenitors; granulocyte progenitors",
    "mature to granulocytes and then neutrophils (a proliferating and a",
    "quiescent pool); early erythroid cells mature to late erythroid cells;",
    "lymphoid progenitors differentiate to B cells. Each division makes two",
    "daughter cells that either renew (fraction rho) or differentiate",
    "(1 - rho). Only the most mature cell type of each lineage dies, at a",
    "uniform rate kDeath, into a total-dead-cells pool. A constant in vitro",
    "drug concentration (nM, dosed into the drug state) acts on every cell",
    "type through one total Emax (EmaxT) and one log(EC50) per cell type:",
    "min(1, EmaxT) inhibits the basal proliferation rate and",
    "max(0, EmaxT - 1) per day is an added first-order kill rate. Default",
    "ini() values are the system parameters calibrated to drug-free cell",
    "kinetics, with every drug effect switched off (EmaxT = 0); the fitted",
    "EmaxT and log(EC50) sets for the named compounds of the source",
    "supplement are carried in the drugEffectParameters metadata list and",
    "are applied with ini() or rxSolve(params = ). No IIV and no residual",
    "error are reported (the model was fitted to donor-mean data by a",
    "genetic algorithm with local search).",
    sep = " "
  )
  reference <- paste(
    "Wilson JL, Lu D, Corr N, Fullerton A, Lu J (2020).",
    "An in vitro quantitative systems pharmacology approach for",
    "deconvolving mechanisms of drug-induced, multilineage cytopenias.",
    "PLoS Comput Biol 16(7):e1007620. doi:10.1371/journal.pcbi.1007620.",
    "PMCID: PMC7402526.",
    "Fluxes, ODEs, rules, events and all system parameter values: S1 Text",
    "(S1Text_in_vitro_equations.pdf, the SimBiology export of the model).",
    "Parameter-name mapping: S3 Text. Drug-effect parameters: S2 Table",
    "(sheet 'log_EC50_parameters_from_fittin'; Table 1 of the paper prints",
    "the same values rounded). The SimBiology project file is deposited at",
    "https://github.com/jenwilson521/Multilineage_InvitroHem_Model.",
    sep = " "
  )
  vignette <- "Wilson_2020_hematopoiesis_invitro"
  units <- list(
    time = "day",
    dosing = "nM (in vitro drug concentration placed in the drug state at time 0)",
    concentration = "cells/mL (cell states); nM (drug)"
  )

  # The in vitro culture is a static system: the drug is added once at day 0
  # and nothing removes it (the SimBiology species 'drug' takes part in no
  # reaction). Following the in vitro convention of
  # Nielsen_2007_semimechanistic_antibiotic_pd.R, the drug state IS the bath
  # concentration and a bolus of amt = concentration (nM) sets it.
  dosing <- c("drug")

  # Cell-type states of Fig 2A / S1 Text. The SimBiology compartment
  # 'InVitro' is 1 mL, so its species amounts (unit 'molecule' = cells) are
  # numerically cells/mL, the unit of Fig 3.
  paper_specific_compartments <- c(
    "drug",
    "hsc",
    "mpp",
    "gmp",
    "mono_prog",
    "mono",
    "gran_prog",
    "gran",
    "neut_prolif",
    "neut_quies",
    "lym_prog",
    "bcell",
    "mk",
    "eryth1",
    "eryth2",
    "dead_cells"
  )

  compartmentData <- list(
    drug = list(analyte = "test compound", units = "nM", specimen = "administration site", verified = TRUE),
    hsc = list(analyte = "hematopoietic stem cells", units = "cells/mL", specimen = "not applicable", verified = TRUE),
    mpp = list(analyte = "multipotent progenitors", units = "cells/mL", specimen = "not applicable", verified = TRUE),
    gmp = list(
      analyte = "granulocyte-macrophage progenitors",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    mono_prog = list(
      analyte = "monocyte progenitors",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    mono = list(analyte = "monocyte-lineage cells", units = "cells/mL", specimen = "not applicable", verified = TRUE),
    gran_prog = list(
      analyte = "granulocyte-lineage progenitors",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    gran = list(
      analyte = "granulocyte-lineage cells",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    neut_prolif = list(
      analyte = "neutrophil-lineage cells, proliferating pool",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    neut_quies = list(
      analyte = "neutrophil-lineage cells, quiescent pool",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    lym_prog = list(analyte = "lymphoid progenitors", units = "cells/mL", specimen = "not applicable", verified = TRUE),
    bcell = list(analyte = "B-lineage cells", units = "cells/mL", specimen = "not applicable", verified = TRUE),
    mk = list(
      analyte = "megakaryocyte-lineage cells",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    eryth1 = list(
      analyte = "early erythroid cells (erythroid-lin I)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    eryth2 = list(
      analyte = "late erythroid cells (erythroid-lin II)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    dead_cells = list(
      analyte = "total dead cells (all cell types)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (human bone-marrow CD34+ cells)",
    n_subjects = 6L,
    n_studies = 1L,
    system = paste(
      "Freeze-thawed primary human bone-marrow-derived CD34+ cells in SFEM II",
      "medium with a custom multilineage cytokine cocktail, in ultra-low",
      "attachment 96-well plates at 37 C, 85% relative humidity, 5% CO2;",
      "cell types counted by 17-parameter flow cytometry with TruCount beads."
    ),
    disease_state = "not applicable (in vitro)",
    dose_range = paste(
      "Drug-free kinetics: six donor samples, counts on days 0, 2, 3, 4, 5, 6",
      "and 9 (system parameters fitted to days 2-6). Drug treatment: 51",
      "compounds, each at 0.2, 1, 5, 25, 100, 500 and 2500 nM added on day 0,",
      "counted on day 6 and normalized to vehicle-control wells; 2-7 donor",
      "samples per compound (S1 Table)."
    ),
    notes = paste(
      "n_subjects counts the six donor samples of the drug-free kinetics",
      "experiment that calibrated the system parameters. The system",
      "parameters and the initial cell counts were fitted to the donor mean,",
      "so the model describes a typical donor and carries no",
      "between-donor variability."
    )
  )

  # Drug-effect parameter sets for the 14 named compounds of S2 Table (5-FU is
  # keyed as fluorouracil), sheet
  # 'log_EC50_parameters_from_fittin' (full precision; Table 1 of the paper
  # prints the eight reference compounds rounded to 3 decimals, with EC50 =
  # exp(log EC50) in nM). Names match the ini() entries, so a set is applied
  # with rxode2::rxSolve(mod, events, params = drugEffectParameters$docetaxel)
  # or `mod |> ini(drugEffectParameters$docetaxel)`. The 37 anonymized
  # development compounds of S2 Table (class2_drug1 ... Pi3K_drug6) are not
  # carried. Thalidomide is the negative control: the paper notes an Emax
  # model could not describe its essentially flat responses (S7 Fig).
  drugEffectParameters <- list(
    abemaciclib = c(
      lec50_hsc = 4.711324803,
      lec50_mpp = 3.831302656,
      lec50_gmp = 3.619698744,
      lec50_gran_prog = 3.501970636,
      lec50_gran = 3.660765516,
      lec50_mono_prog = 3.749888807,
      lec50_mono = 3.832687789,
      lec50_neut = 4.192878639,
      lec50_eryth1 = 3.373606071,
      lec50_eryth2 = 3.80847646,
      lec50_mk = 6.186255302,
      lec50_lym_prog = 3.59750578,
      lec50_bcell = 4.237421368,
      emaxt_hsc = 0.9101713694,
      emaxt_mpp = 0.001563390276,
      emaxt_gmp = 1.03420858,
      emaxt_gran_prog = 0.002294951158,
      emaxt_gran = 0.001633925821,
      emaxt_mono_prog = 0.001661055877,
      emaxt_mono = 0.2307115982,
      emaxt_neut = 0.2345890627,
      emaxt_eryth1 = 1.024122326,
      emaxt_eryth2 = 0.03326860856,
      emaxt_mk = 0.0679720349,
      emaxt_lym_prog = 0.001599354253,
      emaxt_bcell = 0.6494149932
    ),
    dinaciclib = c(
      lec50_hsc = 4.258027775,
      lec50_mpp = 2.81424342,
      lec50_gmp = 3.894607668,
      lec50_gran_prog = 6.307529987,
      lec50_gran = 5.123352624,
      lec50_mono_prog = 3.444975375,
      lec50_mono = 3.345723456,
      lec50_neut = 2.653191726,
      lec50_eryth1 = 3.249216204,
      lec50_eryth2 = 3.489836297,
      lec50_mk = 3.144705193,
      lec50_lym_prog = 3.244750439,
      lec50_bcell = 3.727009999,
      emaxt_hsc = 1.030573267,
      emaxt_mpp = 0.693657522,
      emaxt_gmp = 0.0007410083598,
      emaxt_gran_prog = 0.0001354584935,
      emaxt_gran = 2.009771969,
      emaxt_mono_prog = 0.0003875483206,
      emaxt_mono = 0.4821666235,
      emaxt_neut = 0.2814839354,
      emaxt_eryth1 = 1.115485403,
      emaxt_eryth2 = 0.07095615749,
      emaxt_mk = 0.3310756894,
      emaxt_lym_prog = 1.149330798,
      emaxt_bcell = 0.7735883934
    ),
    docetaxel = c(
      lec50_hsc = 3.51844587,
      lec50_mpp = 1.548902756,
      lec50_gmp = 4.634167496,
      lec50_gran_prog = 4.850748535,
      lec50_gran = 2.067777973,
      lec50_mono_prog = 2.064921858,
      lec50_mono = 2.200515904,
      lec50_neut = 2.617350931,
      lec50_eryth1 = 2.096297529,
      lec50_eryth2 = 1.840319855,
      lec50_mk = 2.05267893,
      lec50_lym_prog = 2.86979661,
      lec50_bcell = 2.340466075,
      emaxt_hsc = 1.300402455,
      emaxt_mpp = 2.013795883,
      emaxt_gmp = 0.0004846578997,
      emaxt_gran_prog = 0.0004131229214,
      emaxt_gran = 2.762421737,
      emaxt_mono_prog = 1.45925997,
      emaxt_mono = 0.197808827,
      emaxt_neut = 0.231903643,
      emaxt_eryth1 = 1.725857969,
      emaxt_eryth2 = 0.01593846324,
      emaxt_mk = 0.1260644309,
      emaxt_lym_prog = 0.001600646328,
      emaxt_bcell = 1.018507819
    ),
    paclitaxel = c(
      lec50_hsc = 4.1319688,
      lec50_mpp = 3.746622945,
      lec50_gmp = 3.339667199,
      lec50_gran_prog = 4.587095141,
      lec50_gran = 2.758229231,
      lec50_mono_prog = 3.209525935,
      lec50_mono = 3.40131544,
      lec50_neut = 3.222542508,
      lec50_eryth1 = 2.731510432,
      lec50_eryth2 = 3.004688233,
      lec50_mk = 3.946080437,
      lec50_lym_prog = 2.459235153,
      lec50_bcell = 3.602646635,
      emaxt_hsc = 1.188440711,
      emaxt_mpp = 0.001756840393,
      emaxt_gmp = 0.001389823796,
      emaxt_gran_prog = 0.0004103865417,
      emaxt_gran = 2.360056111,
      emaxt_mono_prog = 1.555944937,
      emaxt_mono = 0.2582104984,
      emaxt_neut = 0.2139422239,
      emaxt_eryth1 = 1.466405598,
      emaxt_eryth2 = 0.07097545357,
      emaxt_mk = 0.2593265162,
      emaxt_lym_prog = 0.001599782424,
      emaxt_bcell = 0.9143363023
    ),
    palbociclib = c(
      lec50_hsc = 5.109941454,
      lec50_mpp = 5.869697062,
      lec50_gmp = 4.000351341,
      lec50_gran_prog = -2.250044504,
      lec50_gran = 3.031695734,
      lec50_mono_prog = 3.510945528,
      lec50_mono = 3.672134079,
      lec50_neut = 4.316898882,
      lec50_eryth1 = 3.749637007,
      lec50_eryth2 = 4.6005249,
      lec50_mk = 6.631800484,
      lec50_lym_prog = 3.177080281,
      lec50_bcell = 4.505910301,
      emaxt_hsc = 0.7761049786,
      emaxt_mpp = 0.0002640117243,
      emaxt_gmp = 0.8842506996,
      emaxt_gran_prog = 0.04704907626,
      emaxt_gran = 0.0003352269579,
      emaxt_mono_prog = 0.0003507590011,
      emaxt_mono = 0.1420773716,
      emaxt_neut = 0.1837259603,
      emaxt_eryth1 = 0.0003824303239,
      emaxt_eryth2 = 0.1176796618,
      emaxt_mk = 0.07370733456,
      emaxt_lym_prog = 0.0003199935932,
      emaxt_bcell = 0.4993894187
    ),
    pictilisib = c(
      lec50_hsc = 5.667583054,
      lec50_mpp = 5.909983376,
      lec50_gmp = 0.3408345984,
      lec50_gran_prog = 4.533105026,
      lec50_gran = -1.788901954,
      lec50_mono_prog = 4.758983041,
      lec50_mono = 6.931443303,
      lec50_neut = 4.881083128,
      lec50_eryth1 = 5.250513079,
      lec50_eryth2 = 5.252265871,
      lec50_mk = 6.975668458,
      lec50_lym_prog = 1.336279649,
      lec50_bcell = 4.514191861,
      emaxt_hsc = 0.4860296008,
      emaxt_mpp = 0.0002957530756,
      emaxt_gmp = 0.07081620335,
      emaxt_gran_prog = 0.6053350371,
      emaxt_gran = 0.2318310913,
      emaxt_mono_prog = 1.399145178,
      emaxt_mono = 0.03930542451,
      emaxt_neut = 0.1646835097,
      emaxt_eryth1 = 0.0003105398408,
      emaxt_eryth2 = 0.1356470715,
      emaxt_mk = 0.1226939349,
      emaxt_lym_prog = 0.0003199844664,
      emaxt_bcell = 0.5142336741
    ),
    ribociclib = c(
      lec50_hsc = 5.920833686,
      lec50_mpp = 6.551951895,
      lec50_gmp = 4.946387916,
      lec50_gran_prog = 3.06424402,
      lec50_gran = 3.098988983,
      lec50_mono_prog = 4.330381843,
      lec50_mono = 5.194360396,
      lec50_neut = 5.180577642,
      lec50_eryth1 = 4.484212334,
      lec50_eryth2 = 5.404239333,
      lec50_mk = 7.623036458,
      lec50_lym_prog = 2.512795293,
      lec50_bcell = 3.578737904,
      emaxt_hsc = 0.4888850628,
      emaxt_mpp = 0.0002366220938,
      emaxt_gmp = 0.7603567627,
      emaxt_gran_prog = 0.06432245126,
      emaxt_gran = 0.0003594020558,
      emaxt_mono_prog = 0.0003519262533,
      emaxt_mono = 0.1953403616,
      emaxt_neut = 0.1390651472,
      emaxt_eryth1 = 0.0005102879083,
      emaxt_eryth2 = 0.1078867609,
      emaxt_mk = 0.04812250075,
      emaxt_lym_prog = 0.0003199875775,
      emaxt_bcell = 0.3049300614
    ),
    thalidomide = c(
      lec50_hsc = 5.244680806,
      lec50_mpp = 6.577029215,
      lec50_gmp = 5.176112465,
      lec50_gran_prog = 4.856744036,
      lec50_gran = 3.60300118,
      lec50_mono_prog = 5.720450352,
      lec50_mono = 6.542986115,
      lec50_neut = -1.198050466,
      lec50_eryth1 = 6.051312058,
      lec50_eryth2 = 7.229812296,
      lec50_mk = 5.773881904,
      lec50_lym_prog = 2.999561729,
      lec50_bcell = 4.867978351,
      emaxt_hsc = 3.736757144e-06,
      emaxt_mpp = 1.114973012e-05,
      emaxt_gmp = 0.0001084520654,
      emaxt_gran_prog = 0.000120653332,
      emaxt_gran = 0.000103387226,
      emaxt_mono_prog = 4.935536354e-05,
      emaxt_mono = 6.079660939e-05,
      emaxt_neut = 0.01134197641,
      emaxt_eryth1 = 3.434588858e-05,
      emaxt_eryth2 = 3.511905878e-05,
      emaxt_mk = 3.881856916e-05,
      emaxt_lym_prog = 6.400051294e-05,
      emaxt_bcell = 0.002415304979
    ),
    fluorouracil = c(
      lec50_hsc = 7.552914559,
      lec50_mpp = 6.844588145,
      lec50_gmp = 6.623942605,
      lec50_gran_prog = 5.711427647,
      lec50_gran = -2.145573184,
      lec50_mono_prog = 5.180950027,
      lec50_mono = 5.699093879,
      lec50_neut = 6.940275133,
      lec50_eryth1 = 6.825795047,
      lec50_eryth2 = 7.485672475,
      lec50_mk = 7.014329374,
      lec50_lym_prog = 3.609935322,
      lec50_bcell = 6.993547245,
      emaxt_hsc = 0.7981970609,
      emaxt_mpp = 0.0005283051904,
      emaxt_gmp = 1.56865383,
      emaxt_gran_prog = 0.01620764915,
      emaxt_gran = 0.1807130656,
      emaxt_mono_prog = 0.001582229857,
      emaxt_mono = 0.304772599,
      emaxt_neut = 0.2667800036,
      emaxt_eryth1 = 2.33242665,
      emaxt_eryth2 = 0.03523655741,
      emaxt_mk = 0.322006228,
      emaxt_lym_prog = 0.001599600691,
      emaxt_bcell = 0.8248941531
    ),
    `A_1331852` = c(
      lec50_hsc = 7.578785988,
      lec50_mpp = 6.911570731,
      lec50_gmp = 6.814736578,
      lec50_gran_prog = 6.855782136,
      lec50_gran = -0.376039689,
      lec50_mono_prog = 5.202159337,
      lec50_mono = 2.738341539,
      lec50_neut = 5.990266201,
      lec50_eryth1 = 3.476507397,
      lec50_eryth2 = 7.060296563,
      lec50_mk = 3.182102396,
      lec50_lym_prog = 3.053231516,
      lec50_bcell = 2.570328626,
      emaxt_hsc = 0.08276186517,
      emaxt_mpp = 2.147275184e-06,
      emaxt_gmp = 1.135137006e-06,
      emaxt_gran_prog = 1.710665614e-06,
      emaxt_gran = 2.250390026,
      emaxt_mono_prog = 1.47314964e-05,
      emaxt_mono = 0.06103161557,
      emaxt_neut = 1.509588883e-07,
      emaxt_eryth1 = 2.575039902,
      emaxt_eryth2 = 1.513588302e-06,
      emaxt_mk = 0.1880420977,
      emaxt_lym_prog = 1.599966892e-05,
      emaxt_bcell = 0.1187395267
    ),
    bortezomib = c(
      lec50_hsc = 4.655311562,
      lec50_mpp = 3.500696861,
      lec50_gmp = 3.704799016,
      lec50_gran_prog = 4.88212604,
      lec50_gran = 4.228663124,
      lec50_mono_prog = 4.261293129,
      lec50_mono = 4.426473468,
      lec50_neut = 4.684742938,
      lec50_eryth1 = 4.107743431,
      lec50_eryth2 = 4.034117686,
      lec50_mk = 3.645379369,
      lec50_lym_prog = 4.718449023,
      lec50_bcell = 3.5885737,
      emaxt_hsc = 1.344935764,
      emaxt_mpp = 0.003021913156,
      emaxt_gmp = 0.1674680378,
      emaxt_gran_prog = 0.0006095296934,
      emaxt_gran = 1.838076299,
      emaxt_mono_prog = 0.001657737662,
      emaxt_mono = 0.3901188557,
      emaxt_neut = 0.2391037989,
      emaxt_eryth1 = 0.001616387582,
      emaxt_eryth2 = 0.1906344222,
      emaxt_mk = 0.3531905253,
      emaxt_lym_prog = 0.001599822091,
      emaxt_bcell = 0.8983701715
    ),
    cytarabine = c(
      lec50_hsc = 2.85462881,
      lec50_mpp = 4.806373547,
      lec50_gmp = 1.215426253,
      lec50_gran_prog = 1.940300506,
      lec50_gran = 3.486078726,
      lec50_mono_prog = 3.356087735,
      lec50_mono = 4.146236296,
      lec50_neut = 1.944067978,
      lec50_eryth1 = 1.461811539,
      lec50_eryth2 = -2.288262751,
      lec50_mk = 1.925389289,
      lec50_lym_prog = 4.243483903,
      lec50_bcell = 1.99437477,
      emaxt_hsc = 0.6131959251,
      emaxt_mpp = 0.000181030175,
      emaxt_gmp = 1.159696803,
      emaxt_gran_prog = 1.041032835,
      emaxt_gran = 0.000324332537,
      emaxt_mono_prog = 0.0003029573662,
      emaxt_mono = 0.2012096203,
      emaxt_neut = 0.1827840815,
      emaxt_eryth1 = 2.039666345,
      emaxt_eryth2 = 0.02919630388,
      emaxt_mk = 0.4281777896,
      emaxt_lym_prog = 0.0003199807065,
      emaxt_bcell = 0.8386419923
    ),
    lenalidomide = c(
      lec50_hsc = 7.028976261,
      lec50_mpp = 6.667272345,
      lec50_gmp = 8.026583352,
      lec50_gran_prog = 6.843051366,
      lec50_gran = 8.064721942,
      lec50_mono_prog = 6.958656151,
      lec50_mono = 5.675876286,
      lec50_neut = 7.897127742,
      lec50_eryth1 = 4.342340353,
      lec50_eryth2 = 5.568502482,
      lec50_mk = 6.636779676,
      lec50_lym_prog = 3.283696329,
      lec50_bcell = 6.910416573,
      emaxt_hsc = 1.333942771e-08,
      emaxt_mpp = 4.523947226e-07,
      emaxt_gmp = 1.843553003e-07,
      emaxt_gran_prog = 5.977547073e-07,
      emaxt_gran = 9.936703079e-07,
      emaxt_mono_prog = 3.831498518e-07,
      emaxt_mono = 3.682950572e-08,
      emaxt_neut = 4.067798177e-08,
      emaxt_eryth1 = 1.170390116e-07,
      emaxt_eryth2 = 1.482644305e-08,
      emaxt_mk = 2.229662148e-08,
      emaxt_lym_prog = 3.999973068e-06,
      emaxt_bcell = 1.486749347e-07
    ),
    `OTX015` = c(
      lec50_hsc = 5.48380272,
      lec50_mpp = 6.056825397,
      lec50_gmp = 6.745932801,
      lec50_gran_prog = 5.172069569,
      lec50_gran = 5.898000344,
      lec50_mono_prog = 4.733632849,
      lec50_mono = 5.577114608,
      lec50_neut = 4.48080188,
      lec50_eryth1 = 7.076801992,
      lec50_eryth2 = 5.153311214,
      lec50_mk = 8.326127969,
      lec50_lym_prog = 5.486671299,
      lec50_bcell = 6.052038567,
      emaxt_hsc = 0.4012463454,
      emaxt_mpp = 0.0005080572106,
      emaxt_gmp = 2.619652959,
      emaxt_gran_prog = 0.001559797013,
      emaxt_gran = 1.259986464,
      emaxt_mono_prog = 0.001792161047,
      emaxt_mono = 0.3745732323,
      emaxt_neut = 0.1672317832,
      emaxt_eryth1 = 1.947525747,
      emaxt_eryth2 = 0.1203206989,
      emaxt_mk = 0.2692261732,
      emaxt_lym_prog = 1.154153431,
      emaxt_bcell = 0.7675086009
    )
  )

  ini({
    # ---- System parameters: S1 Text, 'Parameters (Model Scoped)' --------
    # Basal proliferation rates kpro_<cell>_0 (1/day). S3 Text maps them to
    # the kappa of Fig 2B and the Methods.
    kpro0_hsc <- 5.2597; label("Basal proliferation rate of HSC (1/day)")                                      # S1 Text: kpro_HSC_0 = 5.2597 1/day
    kpro0_mpp <- 5.4285; label("Basal proliferation rate of MPP (1/day)")                                      # S1 Text: kpro_MPP_0 = 5.4285 1/day
    kpro0_gmp <- 0.27694; label("Basal proliferation rate of GMP (1/day)")                                     # S1 Text: kpro_GMP_0 = 0.27694 1/day
    kpro0_mono_prog <- 0.15714; label("Basal proliferation rate of monocyte progenitors (1/day)")              # S1 Text: kpro_MonoP_0 = 0.15714 1/day
    kpro0_mono <- 1.3532; label("Basal proliferation rate of monocyte-lineage cells (1/day)")                  # S1 Text: kpro_Mono_0 = 1.3532 1/day
    kpro0_gran_prog <- 5.7012; label("Basal proliferation rate of granulocyte progenitors (1/day)")            # S1 Text: kpro_GranP_0 = 5.7012 1/day
    kpro0_gran <- 1.1726e-4; label("Basal proliferation rate of granulocyte-lineage cells (1/day)")            # S1 Text: kpro_Gran_0 = 1.1726E-4 1/day
    kpro0_neut <- 2.6679; label("Basal proliferation rate of proliferating neutrophils (1/day)")               # S1 Text: kpro_Neut_0 = 2.6679 1/day
    kpro0_eryth1 <- 0.0077462; label("Basal proliferation rate of early erythroid cells (1/day)")              # S1 Text: kpro_ErythroidI_0 = 0.0077462 1/day
    kpro0_eryth2 <- 1.9612; label("Basal proliferation rate of late erythroid cells (1/day)")                  # S1 Text: kpro_ErythroidII_0 = 1.9612 1/day
    kpro0_mk <- 1.072; label("Basal proliferation rate of megakaryocyte-lineage cells (1/day)")                # S1 Text: kpro_MK_0 = 1.072 1/day
    kpro0_bcell <- 0.84171; label("Basal proliferation rate of B-lineage cells (1/day)")                       # S1 Text: kpro_B_0 = 0.84171 1/day

    # Renewal fractions rho_<cell> (fraction of daughter cells that renew
    # rather than differentiate). S1 Text also lists renewal_Mono = 0.26505,
    # which no flux uses (monocyte-lineage cells are terminal), so it is not
    # carried.
    renew_hsc <- 0.54755; label("Renewal fraction of HSC (unitless)")                                          # S1 Text: renewal_HSC = 0.54755
    renew_mpp <- 0.21625; label("Renewal fraction of MPP (unitless)")                                          # S1 Text: renewal_MPP = 0.21625
    renew_gmp <- 0.47546; label("Renewal fraction of GMP (unitless)")                                          # S1 Text: renewal_GMP = 0.47546
    renew_mono_prog <- 0.41977; label("Renewal fraction of monocyte progenitors (unitless)")                   # S1 Text: renewal_MonoP = 0.41977
    renew_gran_prog <- 0.49822; label("Renewal fraction of granulocyte progenitors (unitless)")                # S1 Text: renewal_GranP = 0.49822
    renew_gran <- 0.36766; label("Renewal fraction of granulocyte-lineage cells (unitless)")                   # S1 Text: renewal_Gran = 0.36766
    renew_eryth1 <- 0.33835; label("Renewal fraction of early erythroid cells (unitless)")                     # S1 Text: renewal_ErythroidI = 0.33835

    # Branching fractions beta. The MPP lymphoid branch is
    # max(0, 1 - kbranch_Erythroid - kbranch_MK - kbranch_GMP) and the GMP
    # granulocyte branch is 1 - kbranch_Mono (S1 Text fluxes 5 and 21).
    kbranch_gmp <- 0.30481; label("Fraction of differentiating MPP that become GMP (unitless)")                # S1 Text: kbranch_GMP = 0.30481
    kbranch_eryth <- 0.67317; label("Fraction of differentiating MPP that become early erythroid cells (unitless)") # S1 Text: kbranch_Erythroid = 0.67317
    kbranch_mk <- 0.021653; label("Fraction of differentiating MPP that become megakaryocytes (unitless)")     # S1 Text: kbranch_MK = 0.021653
    kbranch_mono <- 0.46395; label("Fraction of differentiating GMP that become monocyte progenitors (unitless)") # S1 Text: kbranch_Mono = 0.46395

    kdeath <- 0.45887; label("Death rate of the terminal cell types (1/day)")                                  # S1 Text: kDeath = 0.45887 1/day
    kdiff_lym <- 0.076318; label("Differentiation rate of lymphoid progenitors to B cells (1/day)")            # S1 Text: kdiff_lym = 0.076318 1/day
    qf_neut <- 0.99974; label("Fraction of the initial neutrophils placed in the quiescent pool (unitless)")   # S1 Text: QF_Neutrophil = 0.99974

    # ---- Initial cell counts (fitted; S1 Text, 'Species - InVitro') -----
    # Methods: the initial counts were fitted within [0.1, 1] x the day-0
    # measurement. Amounts in a 1 mL compartment, i.e. cells/mL.
    bl_hsc <- 275.67; label("Initial HSC count (cells/mL)")                                                    # S1 Text: Hematopoietic Stem Cell = 275.67
    bl_mpp <- 487.71; label("Initial MPP count (cells/mL)")                                                    # S1 Text: MPP = 487.71
    bl_gmp <- 754.32; label("Initial GMP count (cells/mL)")                                                    # S1 Text: GMP = 754.32
    bl_mono_prog <- 86.607; label("Initial monocyte-progenitor count (cells/mL)")                              # S1 Text: Monocyte prog = 86.607
    bl_mono <- 32.419; label("Initial monocyte-lineage count (cells/mL)")                                      # S1 Text: Monocyte-lin = 32.419
    bl_gran_prog <- 34.342; label("Initial granulocyte-progenitor count (cells/mL)")                           # S1 Text: Gran-lin prog = 34.342
    bl_gran <- 251.64; label("Initial granulocyte-lineage count (cells/mL)")                                   # S1 Text: Gran-lin = 251.64
    bl_neut <- 154.07; label("Initial neutrophil-lineage count before the quiescent split (cells/mL)")         # S1 Text: Prolif:Neutrophil-lin = 154.07 (= initNeutrophil)
    bl_lym_prog <- 562.67; label("Initial lymphoid-progenitor count (cells/mL)")                               # S1 Text: Lymphoid prog = 562.67
    bl_bcell <- 55.595; label("Initial B-lineage count (cells/mL)")                                            # S1 Text: B-lin = 55.595
    bl_mk <- 11.616; label("Initial megakaryocyte-lineage count (cells/mL)")                                   # S1 Text: MK-lin = 11.616
    bl_eryth1 <- 497.83; label("Initial early erythroid count (cells/mL)")                                     # S1 Text: Erythroid-lin I = 497.83
    bl_eryth2 <- 6.0643; label("Initial late erythroid count (cells/mL)")                                      # S1 Text: Erythroid-lin II = 6.0643
    bl_dead_cells <- 875.55; label("Initial total dead-cell count (cells/mL)")                                 # S1 Text: totalDeadCells = 875.55

    # ---- Drug-effect parameters (per compound) ---------------------------
    # One total Emax (EmaxT, unitless, bounded [0, 2] in the fit) and one
    # natural-log EC50 (log of nM; bounded [-2.3, 8.5]) per cell type
    # (Methods Eqs 2-3; S1 Text rules 5-29). EmaxT = 0 is the drug-free
    # system, matching the S1 Text defaults (Emax_drug_<cell> = 0.0,
    # log_EC50_<cell> = 1.0). Per-compound fitted sets are in the
    # drugEffectParameters metadata above (S2 Table).
    emaxt_hsc <- fixed(0); label("Total drug Emax on HSC (unitless)")                                          # S1 Text: Emax_drug_HSC = 0.0 (drug-free default)
    emaxt_mpp <- fixed(0); label("Total drug Emax on MPP (unitless)")                                          # S1 Text: Emax_drug_MPP = 0.0
    emaxt_gmp <- fixed(0); label("Total drug Emax on GMP (unitless)")                                          # S1 Text: Emax_drug_GMP = 0.0
    emaxt_mono_prog <- fixed(0); label("Total drug Emax on monocyte progenitors (unitless)")                   # S1 Text: Emax_drug_MonoP = 0.0
    emaxt_mono <- fixed(0); label("Total drug Emax on monocyte-lineage cells (unitless)")                      # S1 Text: Emax_drug_Mono = 0.0
    emaxt_gran_prog <- fixed(0); label("Total drug Emax on granulocyte progenitors (unitless)")                # S1 Text: Emax_drug_GranP = 0.0
    emaxt_gran <- fixed(0); label("Total drug Emax on granulocyte-lineage cells (unitless)")                   # S1 Text: Emax_drug_Gran = 0.0
    emaxt_neut <- fixed(0); label("Total drug Emax on neutrophils, both pools (unitless)")                     # S1 Text: Emax_drug_Neut = 0.0
    emaxt_lym_prog <- fixed(0); label("Total drug Emax on lymphoid progenitors (unitless)")                    # S1 Text: Emax_drug_LymP = 0.0
    emaxt_bcell <- fixed(0); label("Total drug Emax on B-lineage cells (unitless)")                            # S1 Text: Emax_drug_B = 0.0
    emaxt_mk <- fixed(0); label("Total drug Emax on megakaryocyte-lineage cells (unitless)")                   # S1 Text: Emax_drug_MK = 0.0
    emaxt_eryth1 <- fixed(0); label("Total drug Emax on early erythroid cells (unitless)")                     # S1 Text: Emax_drug_ErythroidI = 0.0
    emaxt_eryth2 <- fixed(0); label("Total drug Emax on late erythroid cells (unitless)")                      # S1 Text: Emax_drug_ErythroidII = 0.0

    lec50_hsc <- fixed(1); label("Natural log of the drug EC50 on HSC (log nM)")                               # S1 Text: log_EC50_HSC = 1.0 (placeholder; inert while EmaxT = 0)
    lec50_mpp <- fixed(1); label("Natural log of the drug EC50 on MPP (log nM)")                               # S1 Text: log_EC50_MPP = 1.0
    lec50_gmp <- fixed(1); label("Natural log of the drug EC50 on GMP (log nM)")                               # S1 Text: log_EC50_GMP = 1.0
    lec50_mono_prog <- fixed(1); label("Natural log of the drug EC50 on monocyte progenitors (log nM)")        # S1 Text: log_EC50_MonoP = 1.0
    lec50_mono <- fixed(1); label("Natural log of the drug EC50 on monocyte-lineage cells (log nM)")           # S1 Text: log_EC50_Mono = 1.0
    lec50_gran_prog <- fixed(1); label("Natural log of the drug EC50 on granulocyte progenitors (log nM)")     # S1 Text: log_EC50_GranP = 1.0
    lec50_gran <- fixed(1); label("Natural log of the drug EC50 on granulocyte-lineage cells (log nM)")        # S1 Text: log_EC50_Gran = 1.0
    lec50_neut <- fixed(1); label("Natural log of the drug EC50 on neutrophils (log nM)")                      # S1 Text: log_EC50_Neut = 1.0
    lec50_lym_prog <- fixed(1); label("Natural log of the drug EC50 on lymphoid progenitors (log nM)")         # S1 Text: log_EC50_LymP = 1.0
    lec50_bcell <- fixed(1); label("Natural log of the drug EC50 on B-lineage cells (log nM)")                 # S1 Text: log_EC50_B = 1.0
    lec50_mk <- fixed(1); label("Natural log of the drug EC50 on megakaryocyte-lineage cells (log nM)")        # S1 Text: log_EC50_MK = 1.0
    lec50_eryth1 <- fixed(1); label("Natural log of the drug EC50 on early erythroid cells (log nM)")          # S1 Text: log_EC50_ErythroidI = 1.0
    lec50_eryth2 <- fixed(1); label("Natural log of the drug EC50 on late erythroid cells (log nM)")           # S1 Text: log_EC50_ErythroidII = 1.0
  })

  model({
    # ---- 1. Drug-effect fractions (EC50_unit = 1 nM in S1 Text) ---------
    fdrug_hsc <- drug / (exp(lec50_hsc) + drug)
    fdrug_mpp <- drug / (exp(lec50_mpp) + drug)
    fdrug_gmp <- drug / (exp(lec50_gmp) + drug)
    fdrug_mono_prog <- drug / (exp(lec50_mono_prog) + drug)
    fdrug_mono <- drug / (exp(lec50_mono) + drug)
    fdrug_gran_prog <- drug / (exp(lec50_gran_prog) + drug)
    fdrug_gran <- drug / (exp(lec50_gran) + drug)
    fdrug_neut <- drug / (exp(lec50_neut) + drug)
    fdrug_lym_prog <- drug / (exp(lec50_lym_prog) + drug)
    fdrug_bcell <- drug / (exp(lec50_bcell) + drug)
    fdrug_mk <- drug / (exp(lec50_mk) + drug)
    fdrug_eryth1 <- drug / (exp(lec50_eryth1) + drug)
    fdrug_eryth2 <- drug / (exp(lec50_eryth2) + drug)

    # ---- 2. Anti-proliferation: S1 Text rules 5-16 (Methods Eq 2) -------
    # kpro_<cell> = kpro_<cell>_0 * (1 - min(1, EmaxT) * C / (EC50 + C)).
    # Lymphoid progenitors have no proliferation term (Methods).
    kpro_hsc <- kpro0_hsc * (1 - min(1, emaxt_hsc) * fdrug_hsc)
    kpro_mpp <- kpro0_mpp * (1 - min(1, emaxt_mpp) * fdrug_mpp)
    kpro_gmp <- kpro0_gmp * (1 - min(1, emaxt_gmp) * fdrug_gmp)
    kpro_mono_prog <- kpro0_mono_prog * (1 - min(1, emaxt_mono_prog) * fdrug_mono_prog)
    kpro_mono <- kpro0_mono * (1 - min(1, emaxt_mono) * fdrug_mono)
    kpro_gran_prog <- kpro0_gran_prog * (1 - min(1, emaxt_gran_prog) * fdrug_gran_prog)
    kpro_gran <- kpro0_gran * (1 - min(1, emaxt_gran) * fdrug_gran)
    kpro_neut <- kpro0_neut * (1 - min(1, emaxt_neut) * fdrug_neut)
    kpro_bcell <- kpro0_bcell * (1 - min(1, emaxt_bcell) * fdrug_bcell)
    kpro_mk <- kpro0_mk * (1 - min(1, emaxt_mk) * fdrug_mk)
    kpro_eryth1 <- kpro0_eryth1 * (1 - min(1, emaxt_eryth1) * fdrug_eryth1)
    kpro_eryth2 <- kpro0_eryth2 * (1 - min(1, emaxt_eryth2) * fdrug_eryth2)

    # ---- 3. Cell killing: S1 Text rules 17-29 and fluxes 25-38 ----------
    # Emax_cellkill_<cell> = max(0, EmaxT - 1) * one_over_day, applied as
    # Emax_cellkill * C / (C + EC50) (Methods Eq 3), in 1/day.
    kkill_hsc <- max(0, emaxt_hsc - 1) * fdrug_hsc
    kkill_mpp <- max(0, emaxt_mpp - 1) * fdrug_mpp
    kkill_gmp <- max(0, emaxt_gmp - 1) * fdrug_gmp
    kkill_mono_prog <- max(0, emaxt_mono_prog - 1) * fdrug_mono_prog
    kkill_mono <- max(0, emaxt_mono - 1) * fdrug_mono
    kkill_gran_prog <- max(0, emaxt_gran_prog - 1) * fdrug_gran_prog
    kkill_gran <- max(0, emaxt_gran - 1) * fdrug_gran
    kkill_neut <- max(0, emaxt_neut - 1) * fdrug_neut
    kkill_lym_prog <- max(0, emaxt_lym_prog - 1) * fdrug_lym_prog
    kkill_bcell <- max(0, emaxt_bcell - 1) * fdrug_bcell
    kkill_mk <- max(0, emaxt_mk - 1) * fdrug_mk
    kkill_eryth1 <- max(0, emaxt_eryth1 - 1) * fdrug_eryth1
    kkill_eryth2 <- max(0, emaxt_eryth2 - 1) * fdrug_eryth2

    # ---- 4. Fluxes (S1 Text 'Fluxes'; reaction numbers in comments) ------
    r1 <- kpro_hsc * hsc                                          # reaction_1: HSC proliferation
    r2 <- kpro_hsc * (1 - renew_hsc) * hsc                        # reaction_2: HSC -> MPP
    r3 <- kpro_mpp * mpp                                          # reaction_3: MPP proliferation
    mpp_diff <- kpro_mpp * (1 - renew_mpp) * mpp
    r5 <- kbranch_gmp * mpp_diff                                  # reaction_5: MPP -> GMP
    r6 <- max(0, 1 - kbranch_eryth - kbranch_mk - kbranch_gmp) * mpp_diff # reaction_6: MPP -> lymphoid prog
    r13 <- kbranch_mk * mpp_diff                                  # reaction_13: MPP -> MK
    r14 <- kbranch_eryth * mpp_diff                               # reaction_14: MPP -> erythroid I
    gmp_diff <- kpro_gmp * (1 - renew_gmp) * gmp
    r8 <- kbranch_mono * gmp_diff                                 # reaction_8: GMP -> monocyte prog
    r15 <- (1 - kbranch_mono) * gmp_diff                          # reaction_15: GMP -> gran-lin prog
    r21 <- kpro_gmp * gmp                                         # reaction_21: GMP proliferation
    r9 <- kpro_mono_prog * (1 - renew_mono_prog) * mono_prog      # reaction_9: monocyte prog -> monocyte-lin
    r22 <- kpro_mono_prog * mono_prog                             # reaction_22: monocyte prog proliferation
    r10 <- kpro_gran_prog * (1 - renew_gran_prog) * gran_prog     # reaction_10: gran-lin prog -> gran-lin
    r4 <- kpro_gran_prog * gran_prog                              # reaction_4: gran-lin prog proliferation
    r12 <- kpro_eryth1 * (1 - renew_eryth1) * eryth1              # reaction_12: erythroid I -> erythroid II
    r25 <- kpro_eryth1 * eryth1                                   # reaction_25: erythroid I proliferation
    r19 <- kdiff_lym * lym_prog                                   # reaction_19: lymphoid prog -> B-lin
    r24 <- kpro_mk * mk                                           # reaction_24: MK proliferation
    r26 <- kpro_eryth2 * eryth2                                   # reaction_26: erythroid II proliferation
    r28 <- kpro_gran * gran                                       # reaction_28: gran-lin proliferation
    r31 <- kpro_gran * (1 - renew_gran) * gran                    # reaction_31: gran-lin -> neutrophil
    r30 <- kpro_mono * mono                                       # reaction_30: monocyte-lin proliferation
    r32 <- kpro_bcell * bcell                                     # reaction_32: B-lin proliferation
    r16 <- kpro_neut * neut_prolif                                # reaction_16: proliferating neutrophil proliferation

    # Deaths of the terminal cell types (basal kDeath plus drug kill)
    r7 <- (kdeath + kkill_eryth2) * eryth2                        # reaction_7
    r11 <- (kdeath + kkill_mk) * mk                               # reaction_11
    r17 <- (kdeath + kkill_mono) * mono                           # reaction_17
    r18 <- (kdeath + kkill_neut) * neut_prolif                    # reaction_18
    r20 <- (kdeath + kkill_bcell) * bcell                         # reaction_20
    r23 <- (kdeath + kkill_neut) * neut_quies                     # reaction_23

    # Drug-induced killing of the non-terminal cell types
    kill_hsc <- kkill_hsc * hsc                                   # Kill_HSC
    kill_mpp <- kkill_mpp * mpp                                   # Kill_MPP
    kill_gmp <- kkill_gmp * gmp                                   # Kill_GMP
    kill_mono_prog <- kkill_mono_prog * mono_prog                 # Kill_MonoP
    kill_gran_prog <- kkill_gran_prog * gran_prog                 # Kill_GranP
    kill_gran <- kkill_gran * gran                                # Kill_Gran
    kill_lym_prog <- kkill_lym_prog * lym_prog                    # Kill_LymP
    kill_eryth1 <- kkill_eryth1 * eryth1                          # Kill_EryI

    # ---- 5. ODEs (S1 Text 'ODEs', stoichiometry as printed) --------------
    d/dt(drug) <- 0
    d/dt(hsc) <- r1 - 2 * r2 - kill_hsc
    d/dt(mpp) <- 2 * r2 + r3 - 2 * r5 - 2 * r6 - 2 * r13 - 2 * r14 - kill_mpp
    d/dt(gmp) <- 2 * r5 - 2 * r8 + r21 - 2 * r15 - kill_gmp
    d/dt(mono_prog) <- 2 * r8 - 2 * r9 + r22 - kill_mono_prog
    d/dt(mono) <- 2 * r9 + r30 - r17
    d/dt(gran_prog) <- -2 * r10 + 2 * r15 + r4 - kill_gran_prog
    d/dt(gran) <- 2 * r10 + r28 - 2 * r31 - kill_gran
    d/dt(neut_prolif) <- 2 * r31 + r16 - r18
    d/dt(neut_quies) <- -r23
    d/dt(lym_prog) <- 2 * r6 - r19 - kill_lym_prog
    d/dt(bcell) <- r19 + r32 - r20
    d/dt(mk) <- 2 * r13 + r24 - r11
    d/dt(eryth1) <- -2 * r12 + 2 * r14 + r25 - kill_eryth1
    d/dt(eryth2) <- 2 * r12 + r26 - r7
    d/dt(dead_cells) <- r7 + r11 + r17 + r18 + r20 + r23 +
      kill_eryth1 + kill_mpp + kill_hsc + kill_mono_prog + kill_gmp +
      kill_gran_prog + kill_gran + kill_lym_prog

    # ---- 6. Initial conditions -------------------------------------------
    # S1 Text 'Events': at time >= 0.01 day the neutrophils are split into a
    # quiescent (QF_Neutrophil) and a proliferating (1 - QF_Neutrophil) pool
    # of initNeutrophil = 154.07 cells. The split is applied here at t = 0
    # (see the vignette for the size of that approximation).
    hsc(0) <- bl_hsc
    mpp(0) <- bl_mpp
    gmp(0) <- bl_gmp
    mono_prog(0) <- bl_mono_prog
    mono(0) <- bl_mono
    gran_prog(0) <- bl_gran_prog
    gran(0) <- bl_gran
    neut_prolif(0) <- (1 - qf_neut) * bl_neut
    neut_quies(0) <- qf_neut * bl_neut
    lym_prog(0) <- bl_lym_prog
    bcell(0) <- bl_bcell
    mk(0) <- bl_mk
    eryth1(0) <- bl_eryth1
    eryth2(0) <- bl_eryth2
    dead_cells(0) <- bl_dead_cells

    # ---- 7. Outputs (S1 Text 'Repeated Assignments') ---------------------
    neut <- neut_prolif + neut_quies
    totalViableCells <- neut + lym_prog + mono + eryth1 + eryth2 + mk +
      gran_prog + mpp + hsc + gmp + mono_prog + bcell + gran
    viability <- totalViableCells / (totalViableCells + dead_cells)
  })
}
