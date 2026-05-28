#' nlmixr2lib convention standards
#'
#' Internal register of canonical parameter, compartment, covariate, and
#' residual-error names used by [checkModelConventions()]. The canonical
#' covariate list is parsed at runtime from
#' `inst/references/covariate-columns.md` (installed as
#' `system.file("references", "covariate-columns.md", package = "nlmixr2lib")`),
#' so that register remains the single authoritative source. Remaining fields
#' mirror the `extract-literature-model` skill's `naming-conventions.md` and
#' `vignettes/create-model-library.Rmd`.
#'
#' @keywords internal
#' @noRd
.nlmixr2libConventionsStatic <- list(
  pkParams = c(
    "lka", "lcl", "lvc", "lvp", "lvp2", "lq", "lq2", "lfdepot",
    "lvmax", "lcl_ss", "lcl_time", "lcl_renal", "lcl_nonren",
    # K-PD / single-rate-constant elimination (primary `ini()` form
    # used when no explicit `vc` is estimated). Canonical lkel adopted
    # 2026-05-28 per the naming audit (replaces `lke`).
    "lkel",
    # Canonical lag-time name. Replaces the legacy `lalag`, `llag`
    # forms per the 2026-05-28 naming audit.
    "ltlag",
    # Canonical acrophase (peak-time) name for circadian-IDR templates.
    # Replaces the legacy `ltz` form; semantically distinct from a
    # lag-time even though both denote a time-shift parameter
    # (operator clarification, 2026-05-28).
    "ltacro",
    # Influx (plasma -> tissue ECF/CSF) and efflux (tissue -> plasma)
    # clearances used by physiological CNS-distribution popPK models
    # in which the brain/tumor extracellular fluid is parameterised as
    # a separate compartment connected to plasma central via two
    # asymmetric clearances driven by unbound drug. Used in Campagne
    # 2019 cyclophosphamide mouse CNS penetration popPK.
    "lclin", "lclef"
  ),
  pkBareParams = c(
    "ka", "cl", "vc", "vp", "vp2", "q", "q2", "kel",
    "k12", "k21", "k13", "k31", "fdepot",
    "vmax", "cl_ss", "cl_time", "cl_renal", "cl_nonren",
    # Bare form of the canonical lag-time and acrophase parameters.
    "tlag", "tacro",
    # Bare counterparts of `lclin` / `lclef` (see pkParams).
    "clin", "clef"
  ),
  compartments = c(
    "depot", "central", "peripheral1", "peripheral2", "peripheral3", "effect",
    "target", "complex", "total_target",
    # Semi-physiological liver and kidney compartments used by paper-specific
    # extraction-ratio first-pass models (Xie_2019_agomelatine) and
    # whole-organ PBPK extractions (Ayyar_2024_givosiran,
    # Gilkey_2015_DiRnanoparticle). Always use the full English name -
    # never `liv` / `kid`.
    "liver", "kidney",
    # Cumulative-hazard state used by time-to-event / dropout sub-models
    # (Girard_2012_pimasertib). The state integrates the instantaneous
    # hazard so that survival = exp(-cumhaz); the source NONMEM idiom is
    # `$MODEL COMP=(CUMHAZ)` with `DADT(<cumhaz>) = HAZARD`.
    "cumhaz",
    # Renal-cortex accumulation compartment used by aminoglycoside
    # nephrotoxicity models (Llanos-Paez_2017_gentamicin). Tracks drug
    # amount sequestered in the renal cortex via saturable uptake from
    # the central compartment plus first-order tubular reabsorption back
    # out (Rougier 2003 / Croes 2011 mechanism).
    "renal_cortex",
    # Cerebrospinal-fluid (CSF) and interstitial-fluid (ISF) physiologic
    # compartments used by mechanistic mAb / target-mediated disposition
    # models with multiple body-fluid distribution volumes
    # (Perez-Ruixo_2025_posdinemab).
    "csf", "isf",
    # Brain / tumor extracellular fluid compartment used by cerebral-
    # microdialysis-based CNS-distribution popPK models in which the
    # ECF compartment carries unbound drug delivered via influx /
    # efflux clearances (Campagne 2019 cyclophosphamide mouse).
    "ecf",
    # Anatomic brain-region compartments used by mAb brain-distribution
    # PK models (Grimm_2023_trontinemab, Grimm_2023_gantenerumab) and
    # other mechanistic brain-PBPK extractions (Xie_2000_m3g_rat,
    # Stevens_2012_remoxipride). Each state holds the extracellular
    # drug concentration in the named region; total brain concentration
    # including residual plasma is derived as `Cbrain_<region>` in
    # model(). The `brain_<region>` namespace was adopted 2026-05 to
    # disambiguate brain-anatomical compartments from same-named
    # non-brain compartments (e.g., renal cortex). The bare region
    # names (`cerebellum`, `hippocampus`, `striatum`, `choroid_plexus`,
    # `brain_ecf`) are deprecated in favour of the prefixed forms;
    # `brain_csf` replaces the older `brain_ecf` for the cerebrospinal-
    # fluid compartment per the 2026-05-28 naming audit.
    "brain_cerebellum", "brain_hippocampus", "brain_striatum",
    "brain_cortex", "brain_choroid_plexus", "brain_csf", "brain_deep",
    # Physiological brain sub-compartments used by hybrid physiology-
    # based PK-PD models that resolve the blood-brain barrier transport
    # as two coupled states: drug in cerebral capillary blood
    # (brain_vascular, volume Vbv, fed by cerebral blood flow CLbv from
    # systemic central) and drug in brain tissue beyond the BBB
    # (brain_extravascular, volume Vbev, fed via the BBB-clearance term
    # CLbev which scales the unbound concentration on each side via the
    # fixed fu_plasma and fu_brain fractions). Used by
    # Johnson_2011_olanzapine_rat.
    "brain_vascular", "brain_extravascular",
    # Friberg-style myelosuppression circulating-cell compartment
    # (Friberg 2002 paclitaxel and derivatives). The terminal
    # compartment of a `precursor1 ... precursorN -> circ` maturation
    # chain; replaces a paper-naming `central` for circulating
    # neutrophils / platelets / lymphocytes when the model is a
    # maturation chain rather than a classical-PK central compartment.
    # Suffix-form `circ_<celltype>` (e.g., `circ_anc`, `circ_plt`) is
    # accepted for paired-output multi-cell-type models.
    "circ",
    # Urinary-excretion compartment used by renally cleared small
    # molecules (single bare `urine` for the parent drug, and per-
    # metabolite `urine_<metab>` for parent-plus-metabolite renal-
    # elimination models like Allegaert 2015 / Cook 2016 paracetamol
    # phase-II conjugates).
    "urine",
    # Cao 2013 mAb mPBPK family (Cao_2013_*, Yuan_2019_concizumab) uses
    # paper-anatomical compartment names that are an explicit exception
    # to the standard `central` / `peripheral1` / `peripheral2`
    # convention. The physiological meaning of plasma / tight / leaky /
    # lymph is load-bearing and would be lost under the generic
    # `peripheralN` renaming. Codified 2026-05-28 per the naming audit.
    "plasma", "tight", "leaky", "lymph",
    # Gallbladder / biliary recirculation compartment used by
    # enterohepatic-circulation (EHC) popPK models (Ide_2009_pravastatin
    # and similar). Drug accumulates from the central compartment via
    # biliary excretion (k12) and re-enters central after a delay
    # (k21 gated by gallbladder-emptying time tg), producing the
    # characteristic second-peak phenomenon.
    "gallbladder",
    # Soluble vascular endothelial growth factor receptor 2 (sVEGFR2)
    # plasma compartment used by indirect-response biomarker PD models
    # for angiogenesis inhibitors (Ait-Oudhia 2016, Hansson 2013a).
    # `svegfr3` is the sibling sVEGFR-3 turnover compartment used
    # alongside sVEGFR-2 in the Hansson 2013 sunitinib biomarker
    # panel; registered 2026-05-28 per the naming audit.
    "svegfr2", "svegfr3",
    # Tumor volume / size compartment used by tumor growth inhibition
    # (TGI) models in oncology (Ait-Oudhia 2016, NA_NA_sunitinib,
    # Schindler 2016, Wilbaux 2015 and similar). `tumor_size` is the
    # snake-case canonical output-state name registered by the 2026-
    # 05-28 audit for the TGI template family (tgi_no_sat_*, tgi_sat_*,
    # Ouerdani 2015 pazopanib, Mazzocco 2015, Zecchin 2016, Wilson
    # 2015 sunitinib_irinotecan_mouse). `carrying_capacity` is the
    # Gompertz / generalised-logistic ceiling K used alongside
    # tumor_size in saturable-growth variants.
    "tumor", "tumor_size", "carrying_capacity",
    # Oncology TGI cell-cycle decomposition states (Simeoni 2004 /
    # Wilson 2015 sunitinib_irinotecan_mouse): drug-driven killing
    # routes proliferating cells through three damaged-cell
    # transit compartments before the cells are cleared. Codified
    # 2026-05-28 per the naming audit.
    "cycling_cells", "damaged_cells1", "damaged_cells2", "damaged_cells3",
    # Endogenous plasma metabolic species used by glucose / lactate
    # turnover sub-models with drug-stimulated production (Oualha 2014
    # epinephrine: glucose zero-order production is stimulated by Ep
    # via an Emax function, and lactate is produced at the rate of
    # glucose elimination and itself first-order eliminated). Each
    # state holds a concentration (mmol/L) rather than an amount,
    # mirroring the source paper's mass-balance parameterisation.
    "glucose", "lactate",
    # Plasma non-esterified fatty acids (NEFA / free fatty acids) used
    # by lipid-turnover PD models with feedback control (Ahlstrom 2010:
    # NiAc inhibits hydrolysis of TG to NEFA; NEFA formation is also
    # suppressed by a moderator transit chain (precursor1..precursor8)
    # representing insulin-like delayed feedback, with a NiAc-independent
    # capillary release term setting the lower physiological limit).
    # State holds a concentration (mmol/L) rather than an amount.
    "nefa",
    # Purine metabolism PD compartments used by semi-mechanistic
    # xanthine / uric-acid turnover models (Hill-McManus 2017
    # doi:10.1111/bcp.13427). `xanthine` and `urate` hold serum amounts
    # (mg); `xanthine_urine` and `urate_urine` hold cumulative urinary
    # excretion amounts (mg) integrated from CLX / CLUA renal-clearance
    # outflows for direct comparison with 24-h urinary collection data.
    "xanthine", "urate", "xanthine_urine", "urate_urine",
    # Multistate Tuberculosis Pharmacometric (MTP) bacterial-state
    # compartments used by TB time-kill / popPK-PD models
    # (Clewe 2016 rifampicin, Chen 2017 mouse, Clewe 2018 TB MTP GPDI
    # in vitro, Svensson 2016 / Wicha 2018 rifampicin). The bacterial
    # population is partitioned into three states: `fast` / `fbugs`
    # (fast-multiplying), `slow` / `sbugs` (slow-multiplying), and
    # `nonm` / `nbugs` (non-multiplying). The original Clewe series
    # uses the `*bugs` form; later Svensson / Wicha rifampicin papers
    # use the bare `fast` / `slow` / `nonm` form. Both are canonical
    # under the MTP exception, registered 2026-05-28 per the naming
    # audit.
    "fast", "slow", "nonm", "fbugs", "sbugs", "nbugs",
    # Enzyme-induction reservoirs used by autoinduction popPK models
    # (Clewe 2015 / Svensson 2016 rifampicin: enz_pool drives
    # time-varying CL via an indirect-response loop on the central
    # enzyme pool; Wicha 2018 / Svensson 2018 rifampicin: bare
    # `enzyme` compartment for the autoinduction mass-action term).
    # Registered 2026-05-28 per the naming audit.
    "enzyme", "enz_pool",
    # DAS28 disease-activity score output compartment used by
    # rheumatoid-arthritis PD models (Frey 2013 tocilizumab,
    # Ma 2020 sarilumab DAS28-CRP). Single PD output, paper-named.
    # Registered 2026-05-28 per the naming audit.
    "das28",
    # Body-weight PD output compartment used by drug-induced weight-
    # change models (Han 2015 sibutramine, Thorsted 2016 somatropin
    # rat: rhGH-driven bodyweight gain). The state is the kg / g body-
    # weight value with first-order turnover driven by drug-modulated
    # production. Registered 2026-05-28 per the naming audit.
    "bw",
    # IGF-1 (insulin-like growth factor 1) plasma biomarker
    # compartment used by somatropin / GH PK/PD models (Thorsted 2016
    # somatropin rat + human). Stimulated by central GH via an Emax
    # function; drives downstream body-weight dynamics. Registered
    # 2026-05-28 per the naming audit.
    "igf1",
    # Gastric / stomach compartment used by gastric-emptying transit
    # models (Guiastrennec 2016 paracetamol, Back 2018 fenofibrate)
    # where the gastric mass-balance is resolved as a distinct state
    # ahead of the duodenal absorption depot. Registered 2026-05-28
    # per the naming audit.
    "stomach",
    # PBPK organ-amount compartments (a_<organ>) used by mass-balance
    # whole-body PBPK extractions (Zurlinden 2016 paracetamol). Each
    # state holds the drug amount in the named organ. The organ
    # vocabulary is spelled out (never abbreviated) per the operator
    # clarification on the 2026-05-28 naming audit:
    #   a_liver       = liver organ (Q_L * (CA - CVL))
    #   a_hepatic     = hepatic intermediate metabolite pool (distinct
    #                   from a_liver; Zurlinden treats hepatic
    #                   metabolite formation as a separate
    #                   compartment from the liver organ for the
    #                   APAP-sulfate and APAP-glucuronide species)
    #   a_fat         = adipose tissue compartment
    #   a_muscle      = muscle tissue
    #   a_kidney      = kidney organ
    #   a_rapidly_perfused = lumped rapidly perfused tissues (highly
    #                   perfused organs not modelled individually)
    #   a_slowly_perfused  = lumped slowly perfused tissues
    #   a_venous      = venous blood compartment
    #   a_arterial    = arterial blood compartment
    #   a_urine       = urinary excretion compartment (note this is
    #                   distinct from the bare `urine` compartment
    #                   above; PBPK models track urine on the
    #                   a_<organ> namespace alongside the other
    #                   amount-tracking compartments, while
    #                   non-PBPK renal-clearance models use the bare
    #                   `urine` form. Both are acceptable for their
    #                   respective conventions.)
    #   a_gut         = gut absorption / intestinal compartment
    "a_liver", "a_hepatic", "a_fat", "a_muscle", "a_kidney",
    "a_rapidly_perfused", "a_slowly_perfused", "a_venous",
    "a_arterial", "a_urine", "a_gut",
    # PBPK organ-vascular concentration compartments (vp_<organ>)
    # used by membrane-limited PBPK extractions where each organ
    # vascular volume is a distinct state (Parhiz 2024 mRNA-LNP,
    # Shah 2012 mAb PBPK). Spelled-out organ vocabulary; abbreviated
    # forms (vp_li / vp_lu / vp_ki / vp_sp / vp_he / vp_ht / vp_mu /
    # vp_sk / vp_ad / vp_bo / vp_br / vp_si / vp_lr / vp_pa / vp_th /
    # vp_po / vp_re / vp_ot) are deprecated per the 2026-05-28 audit.
    "vp_liver", "vp_lung", "vp_kidney", "vp_spleen", "vp_heart",
    "vp_muscle", "vp_skin", "vp_adipose", "vp_bone", "vp_brain",
    "vp_small_intestine", "vp_large_intestine", "vp_pancreas",
    "vp_thymus", "vp_portal", "vp_remainder", "vp_other",
    # Bimodal disease-progression state compartments used by the
    # Delor 2013 Alzheimer mixture-of-fast-and-slow-growths PD model
    # (per-subject mixture weight selects between a fast-progression
    # arm `a_fast` and a slow-progression arm `a_slow`). Distinct
    # from the PBPK perfusion compartments a_rapidly_perfused /
    # a_slowly_perfused above (different mechanistic role).
    "a_fast", "a_slow",
    # Whole-body blood compartments and helper states used by
    # membrane-limited PBPK models:
    #   blood   = whole-body central blood (Parhiz 2024 mRNA-LNP)
    #   bldeg   = blood-pool LNP degradation reservoir (Parhiz 2024)
    #   bcc     = central blood cells (Shah 2012 mAb PBPK)
    #   lnode   = lymph-node return compartment (Shah 2012)
    "blood", "bldeg", "bcc", "lnode",
    # Standard clinical PD-output biomarkers registered as canonical
    # compartments / output-state names so single-output PD models
    # using them pass the relaxed `Cc` rule. Each name is the
    # internationally standardised clinical abbreviation:
    #   ANC, PLT, WBC, RBC = cell counts (absolute neutrophil count,
    #     platelet, white blood cell, red blood cell)
    #   INR, PT, aPTT      = coagulation tests (international
    #     normalised ratio, prothrombin time, activated partial
    #     thromboplastin time)
    #   hb                  = hemoglobin
    #   PSA                 = prostate-specific antigen
    #   P1NP                = procollagen type I N-terminal propeptide
    #     (bone formation biomarker)
    #   OC                  = osteocalcin (bone turnover biomarker)
    #   TT                  = total testosterone / thrombin time
    #     (paper-dependent; both share the TT abbreviation in the
    #     contexts where it appears)
    # Registered 2026-05-28 per the naming audit (operator decision:
    # spell out paper-mechanistic names but the standard clinical
    # abbreviations are themselves canonical and need not be expanded).
    "ANC", "PLT", "WBC", "RBC", "INR", "PT", "aPTT", "hb",
    "PSA", "P1NP", "OC", "TT",
    # Bacterial-count PD output canonicals for TB / antibiotic models.
    # `cfu` is the linear colony-forming-unit count (Clewe 2016
    # rifampicin: cfu = fbugs + sbugs, proportional residual error).
    # `log_cfu` is the log-transformed sputum / culture CFU output
    # (universal TB-PK/PD endpoint, used by the Clewe 2018 / Khan 2015
    # / Mohamed 2016 / Sadouki 2025 / Svensson 2016 / Wicha 2018
    # rifampicin and combination-antibiotic models). The transform
    # base (ln vs log10) is paper-dependent and documented in each
    # source file. Registered 2026-05-28 per the naming audit.
    "cfu", "log_cfu", "MBL",
    # Paper-specific PD-endpoint output state names registered 2026-
    # 05-28 per the naming audit as canonical-compartment / output-
    # state names so single-output PD models that use them pass the
    # relaxed Cc rule.
    #
    # Cognitive / dementia scores:
    #   ADAS_cog  = Alzheimer Disease Assessment Scale - cognitive subscale
    #   ADAS_NORM = ADAS normalised (per-paper rescaling)
    #   cdr_mix   = Clinical Dementia Rating mixture-of-progression-rates output
    # Oncology / tumour endpoints:
    #   tumor_vol = TGI tumour volume output (Lobo 2002, Simeoni 2004)
    #   aescore   = composite adverse-event score (Girard 2012 pimasertib)
    # Ophthalmology:
    #   bcva      = best-corrected visual acuity (Mulyukov 2018 ranibizumab)
    # Pain / symptom scales:
    #   score          = generic pain score (Plan 2012)
    #   vas_pred       = visual-analog-scale prediction (Valitalo 2017 morphine)
    #   fatigue_grade  = fatigue grade (Hansson 2013c sunitinib)
    #   walkDist       = 6-minute walk test distance (Hamuro 2017 DMD)
    #   fev1pp         = FEV1 percent predicted (Harun 2019 cystic fibrosis)
    #   msHeadacheDays = monthly headache-day count (FiedlerKelly 2020 fremanezumab)
    #   migraineDays   = monthly migraine-day count (FiedlerKelly 2020 fremanezumab)
    # Virology / infection:
    #   viralLoad = viral load (Koloskoff 2025 ganciclovir)
    # Imaging / observation models:
    #   prob_roc  = probability output for ROC-style logistic PD models (Shin 2014 sevoflurane)
    # Endocrinology:
    #   prolactin = serum prolactin output (Stevens 2012 remoxipride)
    # Other tracked clinical states:
    #   aaaSize   = abdominal aortic aneurysm size (Sherer 2012 AAA)
    #   cel_count = cell counts in MS lesions (VelezdeMendizabal 2013 multiple sclerosis)
    #   G         = endogenous glucose output (Bizzotto 2016 glucose)
    "ADAS_cog", "ADAS_NORM", "cdr_mix", "tumor_vol", "aescore",
    "bcva", "score", "vas_pred", "fatigue_grade", "walkDist",
    "fev1pp", "msHeadacheDays", "migraineDays", "viralLoad",
    "prob_roc", "prolactin", "aaaSize", "cel_count", "G"
  ),
  # Bare numbered chains (transit / effect / precursor / lat / dar /
  # depot) and metabolite-suffixed compartments are validated
  # separately via .matchesCompartment() so that the registered
  # metabolite list can be honored at runtime; this static regex
  # covers only the numbered-chain patterns. `depot[0-9]+` accommodates
  # parallel-absorption models with two or more depots.
  compartmentRegex = "^(transit|effect|precursor|lat|depot)[0-9]+$",
  # Membrane-limited PBPK sub-compartment pattern: paper-prefix +
  # spelled-out organ name. Recognises the recurring
  # `<sub>_<organ>` shape used in Shah 2012 mAb PBPK and Parhiz 2024
  # mRNA-LNP extractions:
  #   bc   = vascular blood cells
  #   eu   = endosomal unbound
  #   eb   = endosomal FcRn-bound (mAb-FcRn complex)
  #   fr   = endosomal free FcRn
  #   is   = interstitial space
  #   int  = intracellular (Parhiz)
  #   mrna = mRNA pool (Parhiz)
  #   luc  = luciferase reporter (Parhiz)
  # The organ vocabulary mirrors the spelled-out a_<organ> / vp_<organ>
  # canonicals registered above. Registered 2026-05-28 per the naming
  # audit (operator clarification on spelled-out PBPK names).
  pbpkSubCompartmentRegex = "^(bc|eu|eb|fr|is|int|mrna|luc)_(liver|lung|kidney|spleen|heart|muscle|skin|adipose|bone|brain|small_intestine|large_intestine|pancreas|thymus|portal|remainder|other|hepatic|fat|rapidly_perfused|slowly_perfused|venous|arterial|urine|gut)$",
  darCompartmentRegex = "^dar[0-9]+_(central|peripheral[0-9]?)$",
  # Target species in physiologic body-fluid or named peripheral
  # compartments (e.g., target_csf, target_isf, target_peripheral,
  # target_peripheral1, complex_csf, complex_isf, complex_peripheral).
  # Used by mechanistic mAb / TMDD models with multiple distribution
  # volumes; the `_peripheral` variant covers extracellular target /
  # complex states tracked in a numbered peripheral compartment
  # (NA_NA_miridesap, Aguiar 2021 ustekinumab, Sahota 2015 miridesap).
  # Extended 2026-05-28 per the naming audit.
  targetLocationRegex = "^(target|complex)_(csf|isf|peripheral[0-9]?)$",
  observationVar = "Cc",
  # propSd and addSd are the canonical proportional and additive
  # residual-error SDs used with `~ prop(...)`, `~ add(...)`, and the
  # combined `~ prop(...) + add(...)` forms. expSd is the canonical
  # log-scale residual SD used with `~ lnorm(...)` (a distinct error
  # structure where the SD applies on the exponentiated scale), used by
  # popPK papers reporting log-transformed proportional error directly
  # (Cirincione_2017_exenatide, Wu_2024_inotuzumab).
  residualError = c("propSd", "addSd", "expSd"),
  transformPrefixes = c("l", "logit", "probit"),
  # Covariate-effect names match e_<cov>(_<continuation>)+_<param>
  # (canonical) or have an additional trailing token for metabolite /
  # shared / CL-component (e_<cov>_<param>_<suffix>). The pattern
  # accepts up to 6 underscore-separated tokens after the leading "e_"
  # to accommodate compound covariates (RACE_BLACK, ADA_POSITIVE,
  # FORM_CHO_PHASE2). Semantic interpretation of the trailing tokens
  # is handled by .classifyCovEffect().
  covEffectPattern = "^e_[A-Za-z0-9]+(_[A-Za-z0-9]+){1,5}$",
  # Lowercase paper / payload names allowed as a third-token suffix on
  # parameters and compartments for a non-parent species. The set is
  # extended whenever a new ADC payload, target binding partner, or
  # secondary-analyte appears in a model. Keep mutually disjoint from
  # pkBareParams and clComponents to preserve covariate-effect
  # disambiguation.
  registeredMetabolites = c(
    "mmae", "dxd", "sn38", "dm4", "medm4", "mcmmaf",
    "complex", "ige", "il1b", "tab", "nab",
    "dar0", "dar1", "dar2", "dar3", "dar4", "dar5", "dar6", "dar7", "dar8",
    # Small-molecule metabolites of agomelatine (Xie 2019): 3-hydroxy
    # and 7-desmethyl. Suffixes start with a digit; this is fine
    # because the convention check matches on `endsWith(name, "_<metab>")`
    # rather than treating the metabolite name itself as an R identifier.
    "3oh", "7dm",
    # N-desmethyl-bedaquiline metabolite (M2) of bedaquiline
    # (Svensson 2016 DDMODEL00000219).
    "m2",
    # Endoxifen (4-hydroxy-N-desmethyltamoxifen), the major active
    # metabolite of tamoxifen -- Ter Heine 2014.
    "endx",
    # Lidocaine sequential metabolites (DDMODEL00000281, NA_NA_lidocaine):
    # MEGX = monoethylglycinexylidide (LID -> MEGX via CYP1A2/3A4),
    # GX   = glycinexylidide (MEGX -> GX), and 2,6-XYL = 2,6-xylidide
    # (LID -> 2,6-XYL minor pathway). Each metabolite is a separate
    # central compartment with its own apparent volume in the source's
    # ADVAN5 parent + 3-metabolite structure.
    "megx", "gx", "xyl",
    # Morphine-3-glucuronide and morphine-6-glucuronide, the two major
    # glucuronide metabolites of morphine -- Knibbe 2009 DDMODEL00000248.
    "m3g", "m6g",
    # Phase-II conjugates: glucuronide (gluc) and sulphate (sulf).
    # Used for paracetamol-glucuronide / paracetamol-sulphate plasma
    # metabolite compartments in Allegaert 2015 (DDMODEL00000267).
    "gluc", "sulf",
    # Paracetamol (APAP) phase-II conjugate metabolites -- APAP-glucuronide
    # ("apapg") and APAP-sulphate ("apaps") used in the Cook 2016 newborn
    # model (DDMODEL00000271).
    "apapg", "apaps",
    # Colistin, the active polymyxin generated in vivo by hydrolysis of
    # the prodrug colistimethate sodium (CMS). Used as a metabolite
    # suffix in parent-prodrug CMS / metabolite-active-drug colistin
    # popPK models (Leuppi-Taegtmeyer 2019 DDMODEL00000295).
    "col",
    # Dihydroartemisinin, the active metabolite of artesunate
    # (Birgersson 2019 DDMODEL00000297).
    "dha",
    # Hydroxy-itraconazole (OH-ITZ), the major active metabolite of
    # itraconazole produced by CYP3A4 hydroxylation. Used as a metabolite
    # suffix in parent + metabolite simultaneous popPK models
    # (Hennig 2006 Clin Pharmacokinet 45(11):1099-1114; Hennig 2007 BJCP
    # 63(4):438-450).
    "ohi",
    # Doxorubicinol, the C-13 alcohol metabolite of doxorubicin
    # (Kunarajah 2017 paediatric oncology popPK/PD model).
    "doxol",
    # Daunorubicinol (DOL), the C-13 alcohol metabolite of
    # daunorubicin formed primarily by carbonyl reductase 1 (CBR1)
    # in adult AML patients (Varatharajan 2016 Cancer Chemother
    # Pharmacol 78(5):1051-1058 doi:10.1007/s00280-016-3166-8).
    "dol",
    # 25-O-desacetyl rifabutin, the primary active metabolite of
    # rifabutin formed by arylacetamide deacetylase (Hennig 2015
    # AAC doi:10.1128/AAC.01195-15).
    "desrbn",
    # AZ5104 (N-desmethyl osimertinib), an active EGFR-inhibitor
    # metabolite of osimertinib formed predominantly via CYP3A4/5
    # (Brown 2017 BJCP 83(6):1216-1226 doi:10.1111/bcp.13223).
    "az5104",
    # N-desmethyl-selumetinib (also reported as the active selumetinib
    # metabolite, ~3-5-fold more potent for MEK1 inhibition than the
    # parent), formed by oxidative N-demethylation of selumetinib
    # (Patel 2017 CPT Pharmacometrics Syst Pharmacol 6(5):305-314
    # doi:10.1002/psp4.12175).
    "ndsel",
    # Capecitabine sequential metabolites (Urien 2005
    # doi:10.1007/s10928-005-0018-2): 5'-DFCR (5'-deoxy-5-
    # fluorocytidine, formed in the liver by carboxylesterase from
    # capecitabine), 5'-DFUR (5'-deoxy-5-fluorouridine, formed from
    # 5'-DFCR by cytidine deaminase in liver and tumour cells), and
    # 5-FU (5-fluorouracil, formed from 5'-DFUR by thymidine
    # phosphorylase preferentially in tumour tissue). Each metabolite
    # has its own central compartment with apparent volume fixed to
    # 1 L (only output rate constants K23, K34, K40 are identifiable
    # in the source NONMEM ADVAN6 fit).
    "dfcr", "dfur", "5fu",
    # AS(N-1)3' truncated antisense strand of GalNAc-conjugated
    # siRNAs (givosiran and other galnac-siRNA conjugates), formed
    # by removal of the 3'-terminal nucleotide from the antisense
    # strand. Treated as the active metabolite that is equipotent
    # with the parent in terms of RISC loading and target mRNA
    # silencing (Ayyar 2024 doi:10.1016/j.xphs.2023.10.026).
    "asn1",
    # 10-monohydroxy derivative (MHD, also "10-hydroxy-carbazepine"),
    # the primary active metabolite of oxcarbazepine produced by
    # cytosolic arylketone reductases (Rodrigues 2017 BJCP
    # doi:10.1111/bcp.13392).
    "mhd",
    # Stereo-isomer (R / S) suffixes for enantiomer-resolved popPK
    # models in which both enantiomers are followed in plasma but no
    # interconversion is modelled (e.g. Valitalo 2017 ketorolac BJCP
    # doi:10.1111/bcp.13311). Treated as "non-parent analyte" suffixes
    # under the same registry as metabolites; neither enantiomer is the
    # parent.
    "r", "s",
    # Roflumilast N-oxide, the active metabolite of roflumilast that
    # contributes about 90% of total PDE4 inhibitory activity (tPDE4i).
    # Used as a metabolite suffix in parent-plus-metabolite popPK models
    # where parent (roflumilast) and metabolite (N-oxide) are fitted
    # jointly but with independent apparent absorption parameters
    # (Lahu 2010 doi:10.2165/11536600-000000000-00000).
    "noxide",
    # Combined acetaminophen cysteine + mercapturate compartment used by
    # CYP2E1-oxidation popPK models that lump the two oxidation
    # metabolites (acetaminophen cysteine and acetaminophen mercapturate)
    # into a single observation compartment because the two species are
    # in rapid equilibrium with overlapping disposition (van Rongen 2016
    # Clin Pharmacokinet doi:10.1007/s40262-015-0357-0).
    "cysmer",
    # Glucarpidase (CPG2), the bacterial carboxypeptidase G2 enzyme
    # given as a rescue therapy after high-dose methotrexate (MTX). Not
    # a metabolite of MTX but a co-administered perpetrator that
    # enzymatically hydrolyzes the substrate; the suffix marks the
    # non-parent species in 2-drug PK/PD models (Kimura 2023
    # doi:10.21873/anticanres.16351).
    "cpg2",
    # 3-N-acetyl-3,4-diaminopyridine (3-Ac DAP), the inactive N-acetyl
    # metabolite of 3,4-diaminopyridine (amifampridine) free base
    # produced by N-acetyltransferases. Used in parent-plus-metabolite
    # popPK models for amifampridine in patients with Lambert-Eaton
    # myasthenia (Thakkar 2017 doi:10.1002/psp4.12218).
    "acdap",
    # H4 (active thiol) metabolite of clopidogrel: the pharmacologically
    # active species responsible for P2Y12 receptor inhibition, formed
    # via sequential CYP-mediated oxidation of clopidogrel (CYP2C19 is
    # the dominant isoform for the second oxidation step). The "H4"
    # designator refers to the H4 stereoisomer specifically, which is
    # the antiplatelet-active diastereomer. Used in parent-plus-
    # metabolite popPK models (Danielak 2017 doi:10.1007/s00228-017-2334-z).
    "h4",
    # Mycophenolic acid glucuronide (MPAG, the 7-O-glucuronide phase II
    # metabolite of mycophenolic acid produced by UGT1A9 and UGT2B7).
    # MPAG is the major plasma metabolite of mycophenolic acid after
    # mycophenolate mofetil (MMF) dosing in renal transplant recipients.
    # Used in parent-plus-metabolite popPK models with explicit competitive
    # protein binding and enterohepatic recirculation
    # (de Winter 2009 doi:10.1007/s10928-009-9136-6).
    "mpag",
    # Cell-type suffixes used with Friberg-style `circ_<celltype>`
    # myelosuppression compartments and `precursor1_<celltype>` ...
    # `precursorN_<celltype>` maturation chains for paired-output
    # multi-cell models (Han_2015_decitabine and similar):
    #   anc = absolute neutrophil count, plt = platelet, wbc = white
    #   blood cell. Registered 2026-05-28 per the naming audit.
    "anc", "plt", "wbc",
    # Parent-drug suffixes used with `urine_<X>` (or `central_<X>`)
    # excretion compartments in parent + metabolite renal-elimination
    # models where the parent itself is also tracked in urine
    # (Allegaert_2015_paracetamol: urine_apap; Pierre_2017_morphine:
    # urine_morphine). Registered 2026-05-28 per the naming audit.
    "apap", "morphine",
    # Sibling-drug suffixes for the Hill-McManus 2017 dual-urate-lowering-
    # therapy PKPD model (doi:10.1111/bcp.13427), where febuxostat (`febx`,
    # xanthine oxidase inhibitor) and lesinurad (`lesn`, URAT1 uricosuric)
    # are co-administered and neither is the "parent"; both PK subsystems
    # use canonical compartment / PK-param names with the drug suffix
    # (`central_febx`, `lcl_febx`, `central_lesn`, `lcl_lesn`, etc.).
    # Same precedent as the existing co-administered-perpetrator (`cpg2`)
    # and stereoisomer (`r`, `s`) entries: registered for the
    # `<canonical>_<sibling>` pattern, not as chemical metabolites.
    "febx", "lesn",
    # Sibling-drug suffixes for fixed-combination antimalarial and
    # antibiotic models. `pyra` = pyrimethamine (paired with
    # sulfadoxine in Odongo 2015 / deKock 2017 sulfadoxine-
    # pyrimethamine models, with `depot_pyra` / `central_pyra` /
    # `peripheral1_pyra` PK subsystem). `mer` = meropenem (paired
    # with gentamicin / ciprofloxacin in Sadouki 2025 and with
    # linezolid / vancomycin in Wicha 2017). Registered 2026-05-28
    # per the naming audit.
    "pyra", "mer",
    # Zurlinden 2016 paracetamol PBPK metabolite shorthand:
    # `as` = APAP-sulfate (the same chemical species as the existing
    # `apaps` Cook 2016 suffix but using the Zurlinden notation), and
    # `ag` = APAP-glucuronide (paired with `apapg`). Registered as a
    # separate suffix here so the Zurlinden a_<organ>_as /
    # a_<organ>_ag compartment names pass the metabolite-suffix check
    # without rewriting the source-paper notation. Registered 2026-
    # 05-28 per the naming audit.
    "as", "ag"
  ),
  # Suffixes allowed for multi-component CL parameters. `_ss` denotes
  # the steady-state arm; `_time` denotes the time-varying decay arm.
  # `_renal` denotes the glomerular-filtration / tubular-secretion arm,
  # `_nonren` the non-renal (hepatic / metabolic / extra-renal) arm,
  # used by additive renal-plus-non-renal popPK models for renally
  # cleared small molecules (e.g. Jonckheere 2019 cefepime, where
  # CL_total = CL_renal + CL_nonren is the structural form).
  clComponents = c("ss", "time", "renal", "nonren"),
  # Paper-named mechanistic parameters that don't fit any canonical PK
  # naming pattern but recur across published models. Treated as
  # acceptable bare names (with the usual `l<name>` convention if
  # log-transformed). Add to this list rather than introducing a new
  # ad-hoc pattern.
  paperNamedParams = c(
    "kd", "kd0", "kdes", "kdecay", "krel", "kss", "kint",
    "frac", "alfm", "ksyn", "p", "vd", "kcat", "kpro", "krmr",
    # Transit-absorption naming used by published popPK models that
    # parameterize via mean-absorption-time / fraction-of-MAT
    # (Svensson 2016 bedaquiline DDMODEL00000219, Kovalenko 2020
    # dupilumab, etc.). `mat` = mean absorption time (hours / days);
    # `mtt` = mean transit time; `fr` = fraction of MAT in the transit
    # delay; `ktr` = first-order transit rate constant (= n_transit /
    # MTT for a chain of length n).
    "mat", "mtt", "fr", "ktr",
    # First-order rate constant for in-vivo formation of an active
    # metabolite from the parent central compartment, when the source
    # paper parameterises metabolite formation independently of the
    # parent's total clearance (i.e. fraction-metabolised is implicit
    # in kmet * Vc_parent / CL_parent rather than estimated as Fm).
    # Used in Krause 2017 selexipag (-> ACT-333679) popPK.
    "kmet",
    # Fraction of parent clearance routed to an active metabolite
    # (FM = "fraction metabolised"), estimated as a structural
    # parameter in parent-metabolite joint popPK models. Distinct
    # from `kmet` (formation rate constant): FM is unitless and
    # bounded in (0, 1] while kmet has rate units. Used in
    # Danielak 2017 clopidogrel -> H4 (doi:10.1007/s00228-017-2334-z).
    "fm",
    # Indirect-response / turnover rate constants (Dayneka 1993;
    # Jusko traditions). `kin` = zero-order production; `kout` =
    # first-order elimination of the turnover pool; `kdeg` = first-
    # order degradation (paper-synonym for elimination); `kpin` /
    # `kpout` = precursor-pool production / loss in `indirect_prec_*`
    # templates. Codified 2026-05-28 per the naming audit.
    "kin", "kout", "kdeg", "kpin", "kpout",
    # Canonical Hill / sigmoid-shape coefficient in sigmoidal Emax /
    # Imax functions (Cc^hill / (ec50^hill + Cc^hill)). Codified 2026-
    # 05-28 per the naming audit. Distinct from `gamma` for Friberg
    # myelosuppression feedback / TGI power-law growth exponents,
    # which retain `gamma` as a mechanistic-role designator.
    "hill",
    # Canonical baseline-value parameter for IDR / turnover state
    # initial conditions and TGI initial tumour sizes. Codified 2026-
    # 05-28 per the naming audit (replaces `r0`, `bl`, `base`, `s0`,
    # `ts0`).
    "rbase",
    # Fraction unbound in plasma, used as a fixed unitless multiplier
    # in cerebral-microdialysis-style CNS-distribution models where
    # only free drug crosses the BBB and the BBB transfer term is
    # CLin * fu * Cp. Reported in [0, 1]; usually held fixed at the
    # in-vitro equilibrium-dialysis-derived value (Campagne 2019).
    "fu"
  ),
  requiredUnits = c("time", "dosing", "concentration"),
  requiredMetadata = c("description", "reference", "units"),
  deprecatedResidualError = c(
    "prop.err", "add.err", "propErr", "addErr",
    "err.prop", "err.add"
  ),
  deprecatedIivPrefixes = c("iiv_", "IIV_", "bsv_", "BSV_"),
  # Bare volume names that should be replaced with vc / vp / vp2.
  deprecatedVolumeNames = c("v", "v1", "v2", "v3", "lv", "lv1", "lv2", "lv3"),
  # Deprecated Michaelis-Menten Vmax names.
  deprecatedVmaxNames = c("vm", "lvm"),
  # Deprecated parent-suffix marker. A model that names a parent-side
  # parameter `<base>_adc` should drop the `_adc` suffix; the parent
  # uses the canonical name unsuffixed.
  deprecatedParentSuffix = "_adc"
)

.covariateRegisterCache <- new.env(parent = emptyenv())

.covariateColumnsPath <- function() {
  p <- system.file("references", "covariate-columns.md",
                   package = "nlmixr2lib")
  if (nzchar(p)) return(p)
  stop("Could not locate inst/references/covariate-columns.md in the ",
       "nlmixr2lib package. Check the installation.", call. = FALSE)
}

#' Parse the canonical covariate register from covariate-columns.md.
#'
#' Walks the Markdown register and extracts one entry per H3 heading
#' (`### NAME (**...**)`). For each entry, captures the `Units`, `Type`,
#' `Scope`, `Source aliases`, and `Example models` fields. Aliases whose
#' backticked content is not a bare R identifier (e.g. `DVID = "study1"`)
#' are skipped. Example-model tokens are accepted as backticked file names
#' ending in `.R`; the `.R` suffix is stripped so the value matches the
#' bare model function name used throughout the rest of the package.
#'
#' @param path Path to the markdown file.
#' @return A named list keyed by canonical name. Each entry is a list with
#'   `units`, `type`, `scope` (one of `"general"` / `"specific"` / `NA`),
#'   `aliases` (character vector of alias names), and `example_models`
#'   (character vector of model function names).
#' @keywords internal
#' @noRd
.parseCovariateColumns <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- list()
  current <- NULL
  state <- "idle"

  flush <- function() {
    if (is.null(current)) return(invisible())
    for (nm in current$names) {
      entries[[nm]] <<- list(
        units = current$units %||% "",
        type = current$type %||% "",
        scope = current$scope %||% NA_character_,
        aliases = current$aliases %||% character(),
        example_models = current$example_models %||% character()
      )
    }
  }

  aliasRegex <- "^\\s*-\\s*`([^`]+)`"
  identRegex <- "^[A-Za-z_][A-Za-z0-9_]*$"
  modelFileRegex <- "^[A-Za-z_][A-Za-z0-9_-]*\\.R$"

  extractBacktickedModels <- function(text) {
    toks <- regmatches(text, gregexpr("`([^`]+)`", text))[[1]]
    models <- character()
    for (tok in toks) {
      inner <- gsub("`", "", tok)
      if (grepl(modelFileRegex, inner)) {
        models <- c(models, sub("\\.R$", "", inner))
      }
    }
    models
  }

  for (line in lines) {
    if (startsWith(line, "## ") && !startsWith(line, "### ")) {
      flush()
      current <- NULL
      state <- "idle"
      next
    }
    if (startsWith(line, "### ")) {
      flush()
      heading <- sub("^###\\s+", "", line)
      heading <- sub("\\s*\\(\\*\\*.*\\*\\*\\)\\s*$", "", heading)
      nms <- trimws(strsplit(heading, ",")[[1]])
      nms <- nms[grepl(identRegex, nms)]
      current <- list(names = nms, aliases = character(),
                      example_models = character())
      state <- "header"
      next
    }
    if (is.null(current)) next

    m <- regmatches(line, regexec("^- \\*\\*Units:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$units <- trimws(m[[2]])
      state <- "header"
      next
    }
    m <- regmatches(line, regexec("^- \\*\\*Type:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$type <- trimws(m[[2]])
      state <- "header"
      next
    }
    m <- regmatches(line, regexec("^- \\*\\*Scope:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      scope_raw <- tolower(trimws(sub("\\.$", "", m[[2]])))
      if (scope_raw %in% c("general", "specific")) {
        current$scope <- scope_raw
      }
      state <- "header"
      next
    }
    if (grepl("^- \\*\\*Source aliases:\\*\\*", line)) {
      state <- "aliases"
      after <- sub("^- \\*\\*Source aliases:\\*\\*\\s*", "", line)
      # "none", "none known", "none;"-style declarations have no aliases.
      if (grepl("^none\\b", after, ignore.case = TRUE)) next
      # Capture inline aliases up to the first em-dash prose separator.
      after <- strsplit(after, "\\s+\u2014\\s+", perl = TRUE)[[1]][1]
      inline <- regmatches(after, gregexpr("`([^`]+)`", after))[[1]]
      for (tok in inline) {
        inner <- gsub("`", "", tok)
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
      }
      next
    }
    if (grepl("^- \\*\\*Example models:\\*\\*", line)) {
      state <- "example_models"
      after <- sub("^- \\*\\*Example models:\\*\\*\\s*", "", line)
      current$example_models <- c(current$example_models,
                                  extractBacktickedModels(after))
      next
    }

    if (state == "aliases") {
      m <- regmatches(line, regexec(aliasRegex, line))[[1]]
      if (length(m) == 2) {
        inner <- m[[2]]
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
        next
      }
      if (grepl("^- \\*\\*", line)) {
        state <- "header"
      }
    }

    if (state == "example_models") {
      # Continuation bullet lines in a multi-line Example-models list.
      if (grepl("^\\s+-\\s", line)) {
        current$example_models <- c(current$example_models,
                                    extractBacktickedModels(line))
        next
      }
      if (grepl("^- \\*\\*", line)) {
        state <- "header"
      }
    }
  }
  flush()
  entries
}

.loadCanonicalCovariates <- function(force = FALSE) {
  if (!force && !is.null(.covariateRegisterCache$canonical)) {
    return(.covariateRegisterCache$canonical)
  }
  entries <- .parseCovariateColumns(.covariateColumnsPath())
  .covariateRegisterCache$canonical <- entries
  entries
}

.nlmixr2libConventions <- function() {
  out <- .nlmixr2libConventionsStatic
  out$canonicalCovariates <- .loadCanonicalCovariates()
  out
}

.nlmixr2libCovariateAliasMap <- function() {
  conv <- .loadCanonicalCovariates()
  out <- character()
  for (canon in names(conv)) {
    for (a in conv[[canon]]$aliases) {
      out[a] <- canon
    }
  }
  out
}

# Return TRUE when `name` is a canonical log-transformed PK parameter or
# a metabolite-suffixed PK parameter (`l<base>_<metab>`).
.isPkParam <- function(name, conv) {
  if (name %in% conv$pkParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkParams) return(TRUE)
    }
  }
  FALSE
}

# Return TRUE when `name` is a canonical bare PK parameter or a
# metabolite-suffixed bare PK parameter (`<base>_<metab>`).
.isPkBareParam <- function(name, conv) {
  if (name %in% conv$pkBareParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkBareParams) return(TRUE)
    }
  }
  FALSE
}

# Compartment name validator. Recognizes:
#   - canonical names from conv$compartments
#   - numbered chains via conv$compartmentRegex (transit/effect/precursor/lat/depot)
#   - DAR-numbered ADC isoforms via conv$darCompartmentRegex
#   - target species in physiologic compartments via conv$targetLocationRegex
#   - metabolite-suffixed compartments: <canonical>_<metab>
.matchesCompartment <- function(name, conv) {
  if (name %in% conv$compartments) return(TRUE)
  if (grepl(conv$compartmentRegex, name)) return(TRUE)
  if (grepl(conv$darCompartmentRegex, name)) return(TRUE)
  if (grepl(conv$targetLocationRegex, name)) return(TRUE)
  if (!is.null(conv$pbpkSubCompartmentRegex) &&
      grepl(conv$pbpkSubCompartmentRegex, name)) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$compartments) return(TRUE)
      # Compositions of a numbered-chain prefix with a metabolite
      # suffix are canonical: e.g., `transit1_m3g`, `precursor2_dxd`,
      # `lat1_complex`. Used by formation-delay transit chains feeding
      # a metabolite central compartment (deHoogd 2017 morphine model:
      # transit1_m3g..transit5_m3g and transit1_m6g..transit2_m6g).
      if (grepl(conv$compartmentRegex, base)) return(TRUE)
    }
  }
  FALSE
}

# TRUE when `name` ends with `_<metab>` for any registered metabolite.
.endsWithMetabolite <- function(name, conv) {
  for (metab in conv$registeredMetabolites) {
    if (endsWith(name, paste0("_", metab))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<component>` for any registered CL component.
.endsWithClComponent <- function(name, conv) {
  for (comp in conv$clComponents) {
    if (endsWith(name, paste0("_", comp))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<param>` for any bare PK parameter.
# Used to detect shared-exponent covariate effects like e_wt_cl_q.
.endsWithBarePkParam <- function(name, conv) {
  for (p in conv$pkBareParams) {
    if (endsWith(name, paste0("_", p))) return(TRUE)
  }
  FALSE
}

# Classify a covariate-effect name by its trailing suffix. Returns one of:
#   "two_token"  - matches e_<cov>_<param> with no third-token suffix
#   "metabolite" - matches e_<cov>_<param>_<metab>
#   "shared"     - matches e_<cov>_<param>_<param2> (shared exponent)
#   "component"  - matches e_<cov>_<param>_<component> (multi-CL arm)
#   "unknown"    - has a third-token suffix that doesn't match any
#                  registered category
.classifyCovEffect <- function(name, conv) {
  if (!startsWith(name, "e_")) return(NA_character_)
  if (!grepl(conv$covEffectPattern, name)) return(NA_character_)
  if (.endsWithMetabolite(name, conv)) return("metabolite")
  if (.endsWithClComponent(name, conv)) return("component")
  if (.endsWithBarePkParam(name, conv)) return("shared")
  # Strip the leading `e_` and check whether the rest is a single
  # `<cov>_<param>` pair (no third-token suffix). If yes, two_token;
  # otherwise the name has an unrecognized trailing suffix.
  rest <- substr(name, 3, nchar(name))
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  if (length(parts) == 2) return("two_token")
  "unknown"
}
