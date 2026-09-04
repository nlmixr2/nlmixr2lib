Campbell_2023_manganese_monkey_pbpk <- function() {
  description <- paste(
    "Preclinical (monkey). PBPK (whole-body, 65-ODE, tissue-binding).",
    "Manganese (Mn) disposition in the non-human primate with rapid",
    "association/dissociation to saturable tissue binding sites, recast from the",
    "slow-binding Schroeter 2011 monkey model following the Yoon 2019 adult-rat",
    "structure. Each of 8 tissues (liver, lung, bone, lumped other, globus",
    "pallidus, cerebellum, olfactory bulb, pituitary) carries a vascular",
    "sub-pool, a free-Mn pool and a capacity-limited bound-Mn pool. Includes",
    "regional inhaled-particle deposition (pulmonary, nasal respiratory, nasal",
    "olfactory with direct nose-to-olfactory-bulb transport), a multi-segment",
    "gut with enterocyte sequestration and faecal excretion, and Hill-type",
    "dose-dependent induction of biliary clearance driven by arterial Mn. A",
    "parallel 54Mn radiotracer system (`_mn54`) competes for the same binding",
    "sites, so tracer and bulk Mn are coupled and not separable. Deterministic:",
    "the publication reports no IIV and no residual error. Mn enters as",
    "covariate-driven zero-order diet and inhalation, NOT as dose events; the",
    "54Mn tracer is dosed as events. States start at zero, so a burn-in to",
    "steady state on the background diet is REQUIRED before any study exposure",
    "(the supplement's DSTART).",
    sep = " "
  )
  reference <- paste(
    "Campbell JL, Clewell HJ 3rd, Van Landingham C, Gentry PR, Keene AM,",
    "Taylor MD, Andersen ME. Incorporation of rapid association/dissociation",
    "processes in tissues into the monkey and human physiologically based",
    "pharmacokinetic models for manganese. Toxicol Sci. 2023 Feb",
    "17;191(2):212-226. doi:10.1093/toxsci/kfac123. PMCID: PMC9936208.",
    "Physiological parameters Table 1; chemical-specific parameters Table 2;",
    "the complete ODE listing is in the Supplementary Material.",
    sep = " "
  )
  vignette <- "Campbell_2023_manganese"
  units <- list(
    time          = "h",
    dosing        = "ug (54Mn tracer events; bulk Mn enters as zero-order diet/inhalation covariates)",
    concentration = "ug/L"
  )

  # Bulk (cold) Mn is NOT dosed as events -- it enters as the zero-order
  # diet and inhalation rates built from the exposure-scenario parameters.
  # These are the 54Mn tracer routes Campbell 2023 simulates in the monkey.
  # An oral tracer dose SPLITS: amt * f_diet_uptake into `liver_mn54` and
  # amt * (1 - f_diet_uptake) into `gut_lumen_mn54`, mirroring the
  # supplement's XRADI / XRAGL partition of KDIETR.
  dosing <- c(
    "venous_mn54",      # iv bolus (Furchner 1966, Figure 4B)
    "liver_mn54",       # oral, absorbed fraction (Furchner 1966, Figure 4C)
    "gut_lumen_mn54",   # oral, unabsorbed fraction (Furchner 1966, Figure 4C)
    "depot_ip_mn54",    # intraperitoneal (Dastur 1971, Figure 4A)
    "depot_sc_mn54",    # subcutaneous infusion (Newland 1987, Figure 5A)
    "depot_lung_mn54"   # nebulised inhalation (Newland 1987, Figure 5B)
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every tissue carries three states: `<organ>_vas` (the
  # organ's vascular sub-pool, taken as 3% of organ volume), `<organ>` (FREE
  # Mn, the species that diffuses and drives toxicity) and `bound_<organ>`
  # (Mn on saturable intracellular binding sites). Each state has an
  # `_mn54` twin holding the 54Mn radiotracer; the two compete for the same
  # binding sites. The `depot_*` states hold deposited or injected material
  # sitting at an administration site, not yet absorbed, so they carry no
  # volume and no concentration.
  compartmentData <- list(
    depot_lung = list(
      analyte = "Manganese deposited on pulmonary/tracheobronchial epithelium, not yet absorbed",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    lung_vas = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    lung = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_lung = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    depot_nasal = list(
      analyte = "Manganese deposited on nasal respiratory epithelium, not yet absorbed",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    nasal_respiratory = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    depot_brain = list(
      analyte = "Manganese deposited on nasal olfactory epithelium (direct nose-to-brain route)",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    liver_vas = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    liver = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_liver = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    gut_lumen = list(
      analyte = "Manganese in the gastrointestinal lumen",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    enterocyte = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    gut_lumen_lower = list(
      analyte = "Manganese in the gastrointestinal lumen",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    a_feces = list(
      analyte = "Manganese, cumulative amount excreted in faeces",
      units = "ug", specimen = "faeces", verified = TRUE
    ),
    brain_vascular = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    brain_globus_pallidus = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_globus_pallidus = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    brain_olfactory_bulb = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_olfactory_bulb = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    brain_cerebellum = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_cerebellum = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    pituitary = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_pituitary = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bone_vas = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    bone = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_bone = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    other_vas = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    other = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_other = list(
      analyte = "Manganese, bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    arterial = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    venous = list(
      analyte = "Manganese, free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    depot_sc = list(
      analyte = "Manganese at the subcutaneous injection site",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    a_oral = list(
      analyte = "Manganese, cumulative amount ingested",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    a_inhaled = list(
      analyte = "Manganese, cumulative amount inhaled",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    depot_lung_mn54 = list(
      analyte = "Manganese-54 (radiotracer) deposited on pulmonary/tracheobronchial epithelium, not yet absorbed",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    lung_vas_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    lung_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_lung_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    nasal_respiratory_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    depot_ip_mn54 = list(
      analyte = "Manganese-54 (radiotracer) in the peritoneal space",
      units = "ug", specimen = "administration site", verified = TRUE
    ),
    liver_vas_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    liver_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_liver_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    gut_lumen_mn54 = list(
      analyte = "Manganese-54 (radiotracer) in the gastrointestinal lumen",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    enterocyte_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    gut_lumen_lower_mn54 = list(
      analyte = "Manganese-54 (radiotracer) in the gastrointestinal lumen",
      units = "ug", specimen = "not applicable", verified = TRUE
    ),
    a_feces_mn54 = list(
      analyte = "Manganese-54 (radiotracer), cumulative amount excreted in faeces",
      units = "ug", specimen = "faeces", verified = TRUE
    ),
    brain_vascular_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    brain_globus_pallidus_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_globus_pallidus_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    brain_olfactory_bulb_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_olfactory_bulb_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    brain_cerebellum_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_brain_cerebellum_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    pituitary_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_pituitary_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bone_vas_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    bone_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_bone_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    other_vas_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    other_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    bound_other_mn54 = list(
      analyte = "Manganese-54 (radiotracer), bound to saturable tissue binding sites",
      units = "ug", specimen = "tissue", verified = TRUE
    ),
    arterial_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    venous_mn54 = list(
      analyte = "Manganese-54 (radiotracer), free (unbound)",
      units = "ug", specimen = "whole blood", verified = TRUE
    ),
    depot_sc_mn54 = list(
      analyte = "Manganese-54 (radiotracer) at the subcutaneous injection site",
      units = "ug", specimen = "administration site", verified = TRUE
    )
  )

  population <- list(
    species     = "monkey (rhesus/cynomolgus, Macaca)",
    n_subjects  = NA_integer_,
    n_studies   = 5L,
    weight_range = paste(
      "2.5 kg (Dorman 2006a inhalation cohort, n = 4-6 per exposure group, and",
      "Dastur 1971 ip cohort, n = 12); Table 1's generic monkey reference",
      "weight is 5.0 kg",
      sep = " "
    ),
    disease_state = "healthy",
    dose_range = paste(
      "Inhalation 0, 0.06, 0.3 and 1.5 mg Mn/m3 as MnSO4 aerosol (MMAD 2.0 um,",
      "GSD 1.5, density 2.95 g/cm3), 6 h/day, 5 days/week, 90 days;",
      "54MnCl2 tracer 200 uCi ip, 0.6 uCi iv, 0.6 uCi oral, 200 uCi sc,",
      "24 and 60 uCi nebulised",
      sep = " "
    ),
    notes = paste(
      "NOT a population fit -- a deterministic PBPK calibrated to pooled",
      "study means. Only the tissue influx (kin_*), efflux (kout_*) and",
      "maximal binding capacities (bmax_*) were re-estimated here, by",
      "simultaneous least squares on log model minus log data using the",
      "Subplex variant (NLOPT_LN_SBPLX) of Nelder-Mead in nloptr. The",
      "calibration set was the Dorman 2006a baseline-diet tissue levels and",
      "1.5 mg/m3 90-day inhalation time course (Figure 3A-J) plus the",
      "whole-body 54Mn clearance studies of Dastur 1971 (ip) and Furchner 1966",
      "(iv and oral) (Figure 4A-C). Association/dissociation constants came",
      "from Yoon et al. 2019 (adult rat) and the GI/biliary/deposition-clearance",
      "parameters from Schroeter et al. 2011. Notably, the dose-dependent",
      "brain-uptake term that earlier monkey models needed in the globus",
      "pallidus and pituitary was NOT required by this structure.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Physiological parameters -- Table 1, monkey column
    # =====================================================================
    bw <- fixed(2.5)
    label("Body weight (kg)")  # Results, Monkey model parameterization: "The body weight of 2.5 kg corresponded to the average weight of the monkeys in Dorman et al. (2006a)"; Table 1 footnote a directs use of study-specific weights, its generic monkey value being 5.0 kg

    # ---- tissue volumes as fraction of body weight (Table 1) ------------
    f_blood <- fixed(0.0734)
    label("Blood volume as fraction of body weight (unitless)")  # Table 1, row Blood, monkey
    f_bone <- fixed(0.144)
    label("Bone volume as fraction of body weight (unitless)")  # Table 1, row Bone, monkey (footnote c: set to human value)
    f_brain <- fixed(0.018)
    label("Brain volume as fraction of body weight (unitless)")  # Table 1, row Brain, monkey
    f_liver <- fixed(0.030)
    label("Liver volume as fraction of body weight (unitless)")  # Table 1, row Liver, monkey
    f_lung <- fixed(0.0066)
    label("Lung volume as fraction of body weight (unitless)")  # Table 1, row Lung, monkey

    # ---- brain-region volumes as fraction of BRAIN (Table 1) ------------
    # Footnote d: Dorman 2006a, with cerebellum, olfactory bulb and globus
    # pallidus taken as twice the reported right-hemisphere volume
    # (Supplementary Table S1).
    f_globus_pallidus <- fixed(0.0022)
    label("Globus pallidus volume as fraction of brain (unitless)")  # Table 1, row Globus Pallidus, monkey
    f_cerebellum <- fixed(0.085)
    label("Cerebellum volume as fraction of brain (unitless)")  # Table 1, row Cerebellum, monkey
    f_olfactory_bulb <- fixed(0.00056)
    label("Olfactory bulb volume as fraction of brain (unitless)")  # Table 1, row Olfactory bulb, monkey
    f_pituitary <- fixed(0.00037)
    label("Pituitary volume as fraction of brain (unitless)")  # Table 1, row Pituitary, monkey

    # ---- nasal geometry (Table 1) ---------------------------------------
    sa_nasal_olfactory <- fixed(1.54)
    label("Surface area of nasal olfactory epithelium (cm2/BW^0.75)")  # Table 1, monkey (footnote h: Menache 1997)
    sa_nasal_respiratory <- fixed(13.42)
    label("Surface area of nasal respiratory epithelium (cm2/BW^0.75)")  # Table 1, monkey (footnote h: Menache 1997)
    th_nasal <- fixed(375)
    label("Average nasal tissue thickness (um)")  # Table 1 (footnote j: Conolly 2000)

    # ---- blood flows (Table 1) ------------------------------------------
    qcc <- fixed(19.5)
    label("Cardiac output (L/h/kg^0.75)")  # Table 1, row Cardiac output, monkey
    vent_alveolar <- fixed(30.0)
    label("Alveolar ventilation (L/h/kg^0.75)")  # Table 1, row Pulmonary ventilation, monkey
    f_q_bone <- fixed(0.042)
    label("Bone blood flow as fraction of cardiac output (unitless)")  # Table 1, monkey
    f_q_brain <- fixed(0.065)
    label("Brain blood flow as fraction of cardiac output (unitless)")  # Table 1, monkey
    f_q_liver <- fixed(0.201)
    label("Liver blood flow as fraction of cardiac output (unitless)")  # Table 1, monkey
    f_q_nasal <- fixed(0.01)
    label("Nasal blood flow as fraction of cardiac output (unitless)")  # Table 1, monkey (footnote c)

    # =====================================================================
    # Chemical-specific parameters -- Table 2, monkey column
    # =====================================================================
    # Cellular transport: influx/efflux diffusion rate constants, each
    # scaled as constant * BW^-0.25. All were optimised to the monkey data
    # in this work and are carried over unchanged into the human model.
    kin_liver <- fixed(127.94)
    label("Liver influx rate constant, KINLIVC (BW^0.25/h)")  # Table 2, KINLIVC
    kin_lung <- fixed(6.61)
    label("Lung influx rate constant, KINLUNGC (BW^0.25/h)")  # Table 2, KINLUNGC
    kin_bone <- fixed(37.58)
    label("Bone influx rate constant, KINBONEC (BW^0.25/h)")  # Table 2, KINBONEC
    kin_other <- fixed(0.76)
    label("Rest-of-body influx rate constant, KINOTHC (BW^0.25/h)")  # Table 2, KINOTHC
    kin_globus_pallidus <- fixed(0.092)
    label("Globus pallidus influx rate constant, KINSTC (BW^0.25/h)")  # Table 2, KINSTC (labelled "striatum"; see description)
    kin_cerebellum <- fixed(0.74)
    label("Cerebellum influx rate constant, KINCBC (BW^0.25/h)")  # Table 2, KINCBC
    kin_olfactory_bulb <- fixed(2.04)
    label("Olfactory bulb influx rate constant, KINOBC (BW^0.25/h)")  # Table 2, KINOBC
    kin_pituitary <- fixed(0.016)
    label("Pituitary influx rate constant, KINPTC (BW^0.25/h)")  # Table 2, KINPTC

    kout_liver <- fixed(2.86)
    label("Liver efflux rate constant, KOUTLIVC (BW^0.25/h)")  # Table 2, KOUTLIVC
    kout_lung <- fixed(0.071)
    label("Lung efflux rate constant, KOUTLUNGC (BW^0.25/h)")  # Table 2, KOUTLUNGC
    kout_bone <- fixed(0.51)
    label("Bone efflux rate constant, KOUTBONEC (BW^0.25/h)")  # Table 2, KOUTBONEC
    kout_other <- fixed(0.0049)
    label("Rest-of-body efflux rate constant, KOUTOTHC (BW^0.25/h)")  # Table 2, KOUTOTHC
    kout_globus_pallidus <- fixed(0.027)
    label("Globus pallidus efflux rate constant, KOUTSTC (BW^0.25/h)")  # Table 2, KOUTSTC
    kout_cerebellum <- fixed(0.022)
    label("Cerebellum efflux rate constant, KOUTCBC (BW^0.25/h)")  # Table 2, KOUTCBC
    kout_olfactory_bulb <- fixed(16.69)
    label("Olfactory bulb efflux rate constant, KOUTOBC (BW^0.25/h)")  # Table 2, KOUTOBC
    kout_pituitary <- fixed(0.012)
    label("Pituitary efflux rate constant, KOUTPTC (BW^0.25/h)")  # Table 2, KOUTPTC

    # ---- tissue binding: association rate constants (Table 2) -----------
    # All eight tissues share one value, from Yoon et al. (2019).
    ka_bind <- fixed(0.182)
    label("Mn association rate constant for all tissues, KAST..KAOTH (per ug/L per h)")  # Table 2, KAST/KAOB/KACB/KAPT/KALIV/KALUNG/KABONE/KAOTH all = 0.182 (Yoon et al. 2019)

    # ---- tissue binding: equilibrium dissociation constants -------------
    # Table 2 gives the dissociation RATE constant as ka * KD, with KD equal
    # to 20 or 25 ug/L depending on tissue. Discussion: "Our fitted
    # dissociation rate constants, set at either 0.37 or 0.46 uM (20 or
    # 25 ug/l) depending on the tissue". Encoding KD (not kd) because KD is
    # the quantity the paper argues transfers across tissues and species
    # ("approximately 0.5 uM"); kd is recovered as ka_bind * kd_<organ>.
    kd_liver <- fixed(25)
    label("Mn equilibrium dissociation constant, liver (ug/L)")  # Table 2, KDLIV = KALIV*25 ug/l
    kd_lung <- fixed(25)
    label("Mn equilibrium dissociation constant, lung (ug/L)")  # Table 2, KDLUNG = KALUNG*25 ug/l
    kd_bone <- fixed(25)
    label("Mn equilibrium dissociation constant, bone (ug/L)")  # Table 2, KDBONE = KABONE*25 ug/l
    kd_other <- fixed(20)
    label("Mn equilibrium dissociation constant, rest of body (ug/L)")  # Table 2, KDOTH = KAOTH*20 ug/l
    kd_globus_pallidus <- fixed(20)
    label("Mn equilibrium dissociation constant, globus pallidus (ug/L)")  # Table 2, KDST = KAST*20 ug/l
    kd_cerebellum <- fixed(20)
    label("Mn equilibrium dissociation constant, cerebellum (ug/L)")  # Table 2, KDCB = KACB*20 ug/l
    kd_olfactory_bulb <- fixed(25)
    label("Mn equilibrium dissociation constant, olfactory bulb (ug/L)")  # Table 2, KDOB = KAOB*25 ug/l
    kd_pituitary <- fixed(20)
    label("Mn equilibrium dissociation constant, pituitary (ug/L)")  # Table 2, KDPT = KAPT*20 ug/l

    # ---- maximal binding capacities (Table 2), per L of tissue ----------
    bmax_liver <- fixed(7906.44)
    label("Maximal Mn binding capacity, liver (ug/L tissue)")  # Table 2, BMAXLIVC
    bmax_lung <- fixed(165.46)
    label("Maximal Mn binding capacity, lung (ug/L tissue)")  # Table 2, BMAXLUNGC
    bmax_bone <- fixed(189.76)
    label("Maximal Mn binding capacity, bone (ug/L tissue)")  # Table 2, BMAXBONEC
    bmax_other <- fixed(227.86)
    label("Maximal Mn binding capacity, rest of body (ug/L tissue)")  # Table 2, BMAXBODC
    bmax_globus_pallidus <- fixed(102.38)
    label("Maximal Mn binding capacity, globus pallidus (ug/L tissue)")  # Table 2, BMAXSTC
    bmax_cerebellum <- fixed(437.04)
    label("Maximal Mn binding capacity, cerebellum (ug/L tissue)")  # Table 2, BMAXCBC
    bmax_olfactory_bulb <- fixed(329.89)
    label("Maximal Mn binding capacity, olfactory bulb (ug/L tissue)")  # Table 2, BMAXOBC
    bmax_pituitary <- fixed(0.55)
    label("Maximal Mn binding capacity, pituitary (ug/L tissue)")  # Table 2, BMAXPTC -- deliberately near zero, so pituitary Mn is essentially all FREE; this is what makes predicted pituitary Mn track blood Mn ~108-fold and reproduces Figure 3H

    # ---- regional deposition fractions (Table 2) ------------------------
    # NOTE the olfactory/respiratory airflow split below is TRANSPOSED
    # relative to Table 2 as printed. Table 2 gives FDEPNO = 0.4313*0.91 and
    # FDEPNR = 0.4313*0.09, but Methods (Monkey model parameterization)
    # states the split is "9% olfactory, 91% respiratory (Kepler et al.,
    # 1998)". Scoring both readings against the paper's own published model
    # curves in Figure 3J: the text reading gives olfactory bulb 2.87 ug/g
    # against a published 3.00 (-4.4%), Table 2 as printed gives 21.07
    # (+602%). All nine Figure 3 tissues improve (median 9.3% -> 4.3%). The
    # text reading is used; see the vignette Errata.
    f_dep_lung <- fixed(0.1396)
    label("Fraction of inhaled Mn deposited in lung, FDEPLU (unitless)")  # Table 2, monkey (MPPD ver. 3.04)
    f_dep_nasal_olfactory <- fixed(0.038817)
    label("Fraction of inhaled Mn deposited in nasal olfactory region, FDEPNO (unitless)")  # 0.4313 * 0.09; MPPD head deposition 0.4313 x the 9% olfactory airflow share of Methods (Table 2 prints 0.4313*0.91 -- transposed, see Errata)
    f_dep_nasal_respiratory <- fixed(0.392483)
    label("Fraction of inhaled Mn deposited in nasal respiratory region, FDEPNR (unitless)")  # 0.4313 * 0.91; MPPD head deposition 0.4313 x the 91% respiratory airflow share of Methods (Table 2 prints 0.4313*0.09 -- transposed, see Errata)
    kdep_lung_deep <- fixed(10.0)
    label("Lung epithelium clearance to deep lung tissue, KDEPLUC (BW^0.25/h)")  # Table 2, KDEPLUC (Schroeter et al. 2011)
    kdep_lung_shallow <- fixed(100.0)
    label("Lung epithelium clearance directly to blood, KSHALLUC (BW^0.25/h)")  # Table 2, KSHALLUC (Schroeter et al. 2011)
    kdep_nasal <- fixed(0.0012)
    label("Nasal respiratory epithelium uptake to blood, KDEPNRC (BW^0.25/h)")  # Table 2, KDEPNRC (Schroeter et al. 2011)
    knose_to_bulb <- fixed(0.0035)
    label("Transport from nasal olfactory epithelium to olfactory bulb, KNPOBC (BW^0.25/h)")  # Table 2, KNPOBC (Schroeter et al. 2011)

    # ---- gastrointestinal parameters (Table 2) --------------------------
    f_diet_uptake <- fixed(0.002)
    label("Fraction of dietary Mn bioavailable, FDIETUP (unitless)")  # Table 2, FDIETUP monkey = 2.00E-03. The Results text for the Dastur 1971 simulation instead says "the baseline dietary uptake (0.0002) was used"; the printed parameter table takes authority over the prose (see Errata).
    kgi <- fixed(0.06)
    label("Rate constant of Mn uptake from GI lumen to epithelium, KGI (1/h)")  # Table 2, KGI monkey
    f_enterocyte <- fixed(0.011)
    label("Fraction of gut-lumen Mn taken into the GI epithelium, FENT (unitless)")  # Table 2, FENT monkey
    kent <- fixed(0.0022)
    label("Rate constant of Mn release from GI epithelium by sloughing, KENT (1/h)")  # Table 2, KENT
    kfeces <- fixed(0.6)
    label("Faecal excretion rate constant, KFECES (1/h)")  # Table 2, KFECES

    # ---- biliary elimination and its induction (Table 2, continued) -----
    cl_bile <- fixed(0.051)
    label("Basal biliary clearance of Mn, KBILEC (L/h/BW^0.75)")  # Table 2 (continued), KBILEC (Schroeter et al. 2011)
    q_bile <- fixed(0.00078)
    label("Bile flow, QBILE (L/h/BW)")  # Table 2 (continued), QBILE monkey; used only to report bile Mn concentration (NA for human)
    emax_bile <- fixed(2.50)
    label("Maximal fractional increase in biliary clearance, KBINDUC (unitless)")  # Table 2 (continued), KBINDUC
    km_bile <- fixed(0.027)
    label("Arterial Mn giving half-maximal biliary induction, KM (ug/g)")  # Table 2 (continued), KM. Table 2 labels the units ug/l, but the supplement compares KM against CART = arterial ug/kg / 1000, i.e. ug/g; the ug/g reading is the one consistent with the code (see Errata).
    hill_bile <- fixed(3.00)
    label("Hill coefficient for biliary induction, SLOPE (unitless)")  # Table 2 (continued), SLOPE

    # =====================================================================
    # Scenario inputs -- supplement "Dosing Controls" block. These are set
    # per simulation, not fitted. Defaults give an unexposed animal.
    # =====================================================================
    f_intake_diet <- fixed(0)
    label("Dietary intake factor, INFAC (kg diet/day/BW^0.75) -- NOT REPORTED by Campbell 2023")  # supplement, "INFAC = 99.9; # Dietary Intake Factors (kg/day/BW) (from EPA 1986)". No value appears in the paper, its supplement, or any on-disk source; it is needed ONLY to convert a ppm diet (diet_mn_ppm) to an intake rate. Default 0 leaves that route inert -- check the Rdiet output. Use diet_mn_mgd instead.
    ksc <- fixed(0)
    label("Absorption rate constant from the subcutaneous injection site, KSBQ (1/h) -- NOT REPORTED")  # supplement "Dosing Controls"; needed only for the Newland 1987 sc-infusion simulation (Figure 5A)
    kip <- fixed(0)
    label("Transfer rate constant from the peritoneal space to liver, KIPSLOW (1/h) -- NOT REPORTED")  # supplement "Dosing Controls"; needed only for the Dastur 1971 ip simulation (Figure 4A)
    dose_mn54 <- fixed(1)
    label("Administered 54Mn tracer dose, for expressing whole-body retention as a percentage (ug)")  # scenario input; retention is a ratio, so any positive value works provided it matches the dosed amount
    # =====================================================================
    # Exposure-scenario inputs. Following the inhalation-PBPK precedent of
    # `Boone_2025_vinylchloride_pbpk.R`, these are ini() parameters rather
    # than covariate columns: they are simulation levers, not subject
    # attributes. Defaults give an unexposed subject on the background diet.
    # Any of them may still be made TIME-VARYING by supplying a data column
    # of the same name to rxSolve() -- which is how the Freeland-Graves
    # stepped-diet protocol of Figure 10 is reproduced.
    # =====================================================================
    diet_mn_mgd <- fixed(0)
    label("Dietary Mn intake (mg Mn/day)")  # supplement "Dosing Controls", BDIET. Campbell 2023 does not report the monkey dietary Mn intake numerically, so the default is 0 and a value must be supplied per scenario; see the vignette Errata.
    diet_mn_ppm <- fixed(0)
    label("Dietary Mn concentration (mg Mn/kg diet, i.e. ppm)")  # supplement "Dosing Controls", DDIET. Converted via diet_mn_ppm * f_intake_diet * bw^0.75; REQUIRES f_intake_diet, which the paper does not report, so this route is inert by default. Check the Rdiet output for the realised intake rate.
    diet_on <- fixed(1)
    label("Dietary intake switch: 1 = diet on, 0 = diet withheld (unitless)")  # supplement "Dosing Controls", XDIN
    conc_air_mgm3 <- fixed(0)
    label("Inhaled Mn air concentration during an exposure window (mg Mn/m3)")  # supplement "Dosing Controls", DINH. Campbell 2023's monkey inhalation study used 0, 0.06, 0.3 and 1.5 mg/m3.
    expo_hours_per_day <- fixed(6)
    label("Duration of the daily inhalation exposure window (h)")  # Figures 3 and 6 protocol: "6 h/day, 5 days/week"
    expo_days_per_week <- fixed(5)
    label("Number of consecutive exposure days per 7-day week (days)")  # Figures 3 and 6 protocol: "6 h/day, 5 days/week"
    expo_duration_days <- fixed(0)
    label("Total length of the inhalation exposure phase (days); 0 switches inhalation off")  # Figures 3 and 6 protocol: 90 days
  })

  model({
    # =====================================================================
    # Scaling (supplement, Dynamics "Scaling calculations")
    # =====================================================================
    # ---- tissue volumes (L) ---------------------------------------------
    v_blood <- f_blood * bw
    v_bone  <- f_bone * bw
    v_liver <- f_liver * bw
    v_lung  <- f_lung * bw
    v_brain <- f_brain * bw
    v_gp    <- f_globus_pallidus * v_brain
    v_cb    <- f_cerebellum * v_brain
    v_ob    <- f_olfactory_bulb * v_brain
    v_pit   <- f_pituitary * v_brain
    # Nasal volumes come from area x thickness: cm2 * um -> L needs /1e4
    # (cm2 -> m2 is not the path; the supplement divides by 1e4 then 1e3).
    v_nasal_olf  <- sa_nasal_olfactory * bw^0.75 * th_nasal / 1e4 / 1e3
    v_nasal_resp <- sa_nasal_respiratory * bw^0.75 * th_nasal / 1e4 / 1e3
    v_nose <- v_nasal_olf + v_nasal_resp
    # "Remaining body" (Table 1). Note this subtracts only the four
    # RESOLVED brain regions, so unresolved brain is part of `other`.
    v_other <- bw - v_gp - v_ob - v_cb - v_liver - v_bone - v_lung -
      v_blood - v_nose - v_pit

    # ---- blood flows (L/h) ----------------------------------------------
    qc <- qcc * bw^0.75
    vent_alv <- vent_alveolar * bw^0.75
    q_liver <- f_q_liver * qc
    # Brain flow is apportioned to the resolved regions. The published
    # expression divides the summed region fractions by the brain fraction
    # of body weight; that is dimensionally odd (the region fractions are
    # already fractions OF BRAIN) but it is retained because it is what
    # reproduces Figure 3: the printed form gives globus pallidus, pituitary
    # and cerebellum within 8-13% of the published curves, while dropping
    # the divisor moves them to 42-61% off. See the vignette Errata.
    q_brain <- f_q_brain * qc *
      ((f_globus_pallidus + f_olfactory_bulb + f_cerebellum + f_pituitary) / f_brain)
    q_bone <- f_q_bone * qc
    q_nose <- f_q_nasal * qc
    q_other <- qc - q_bone - q_brain - q_liver - q_nose
    q_bile_abs <- q_bile * bw
    cl_bile_abs <- cl_bile * bw^0.75

    # ---- allometric scaling of first-order rate constants ---------------
    kdep_lu   <- kdep_lung_deep * bw^-0.25
    kshal_lu  <- kdep_lung_shallow * bw^-0.25
    kdep_nr   <- kdep_nasal * bw^-0.25
    knp_ob    <- knose_to_bulb * bw^-0.25
    kin_li <- kin_liver * bw^-0.25
    kin_lu <- kin_lung * bw^-0.25
    kin_bo <- kin_bone * bw^-0.25
    kin_ot <- kin_other * bw^-0.25
    kin_gp <- kin_globus_pallidus * bw^-0.25
    kin_cb <- kin_cerebellum * bw^-0.25
    kin_ob <- kin_olfactory_bulb * bw^-0.25
    kin_pt <- kin_pituitary * bw^-0.25
    kout_li <- kout_liver * bw^-0.25
    kout_lu <- kout_lung * bw^-0.25
    kout_bo <- kout_bone * bw^-0.25
    kout_ot <- kout_other * bw^-0.25
    kout_gp <- kout_globus_pallidus * bw^-0.25
    kout_cb <- kout_cerebellum * bw^-0.25
    kout_ob <- kout_olfactory_bulb * bw^-0.25
    kout_pt <- kout_pituitary * bw^-0.25

    # ---- dissociation RATE constants, kd = ka * KD (Figure 2) -----------
    kdr_li <- ka_bind * kd_liver
    kdr_lu <- ka_bind * kd_lung
    kdr_bo <- ka_bind * kd_bone
    kdr_ot <- ka_bind * kd_other
    kdr_gp <- ka_bind * kd_globus_pallidus
    kdr_cb <- ka_bind * kd_cerebellum
    kdr_ob <- ka_bind * kd_olfactory_bulb
    kdr_pt <- ka_bind * kd_pituitary

    # ---- absolute binding capacities (ug) -------------------------------
    bcap_li <- bmax_liver * v_liver
    bcap_lu <- bmax_lung * v_lung
    bcap_bo <- bmax_bone * v_bone
    bcap_ot <- bmax_other * v_other
    bcap_gp <- bmax_globus_pallidus * v_gp
    bcap_cb <- bmax_cerebellum * v_cb
    bcap_ob <- bmax_olfactory_bulb * v_ob
    bcap_pt <- bmax_pituitary * v_pit

    # ---- FREE binding sites: bulk Mn and 54Mn compete for the same sites
    bf_li <- bcap_li - bound_liver - bound_liver_mn54
    bf_lu <- bcap_lu - bound_lung - bound_lung_mn54
    bf_bo <- bcap_bo - bound_bone - bound_bone_mn54
    bf_ot <- bcap_ot - bound_other - bound_other_mn54
    bf_gp <- bcap_gp - bound_brain_globus_pallidus - bound_brain_globus_pallidus_mn54
    bf_cb <- bcap_cb - bound_brain_cerebellum - bound_brain_cerebellum_mn54
    bf_ob <- bcap_ob - bound_brain_olfactory_bulb - bound_brain_olfactory_bulb_mn54
    bf_pt <- bcap_pt - bound_pituitary - bound_pituitary_mn54

    # =====================================================================
    # Concentrations. Tissue vascular sub-pools are 3% of organ volume.
    # =====================================================================
    c_lung_vas <- lung_vas / (v_lung * 0.03)
    c_lung <- lung / v_lung
    c_nasal_vas <- nasal_respiratory / (v_nasal_resp * 0.03)
    c_liver_vas <- liver_vas / (v_liver * 0.03)
    c_liver <- liver / v_liver
    c_brain_vas <- brain_vascular / (v_brain * 0.03)
    c_gp <- brain_globus_pallidus / v_gp
    c_ob <- brain_olfactory_bulb / v_ob
    c_cb <- brain_cerebellum / v_cb
    c_pit <- pituitary / v_pit
    c_bone_vas <- bone_vas / (v_bone * 0.03)
    c_bone <- bone / v_bone
    c_other_vas <- other_vas / (v_other * 0.03)
    c_other <- other / v_other
    c_art <- arterial / v_blood
    c_ven <- venous / v_blood

    t_c_lung_vas <- lung_vas_mn54 / (v_lung * 0.03)
    t_c_lung <- lung_mn54 / v_lung
    t_c_nasal_vas <- nasal_respiratory_mn54 / (v_nasal_resp * 0.03)
    t_c_liver_vas <- liver_vas_mn54 / (v_liver * 0.03)
    t_c_liver <- liver_mn54 / v_liver
    t_c_brain_vas <- brain_vascular_mn54 / (v_brain * 0.03)
    t_c_gp <- brain_globus_pallidus_mn54 / v_gp
    t_c_ob <- brain_olfactory_bulb_mn54 / v_ob
    t_c_cb <- brain_cerebellum_mn54 / v_cb
    t_c_pit <- pituitary_mn54 / v_pit
    t_c_bone_vas <- bone_vas_mn54 / (v_bone * 0.03)
    t_c_bone <- bone_mn54 / v_bone
    t_c_other_vas <- other_vas_mn54 / (v_other * 0.03)
    t_c_other <- other_mn54 / v_other
    t_c_art <- arterial_mn54 / v_blood
    t_c_ven <- venous_mn54 / v_blood

    # =====================================================================
    # Inputs
    # =====================================================================
    # ---- inhalation square wave (idiom from Boone_2025_vinylchloride_pbpk)
    # MOD(t, P) is written as t - P*floor(t/P), and IF/THEN/ELSE as the
    # complement of a comparison, which rxode2 evaluates to 1 or 0. The wave
    # is ON for the first `expo_hours_per_day` hours of each of the first
    # `expo_days_per_week` days of each 7-day week, until
    # `expo_duration_days` have elapsed.
    hour_of_day <- time - 24 * floor(time / 24)
    day_index <- floor(time / 24)
    day_of_week <- day_index - 7 * floor(day_index / 7)
    expo_on <- (1 - (hour_of_day >= expo_hours_per_day)) *
      (1 - (day_of_week >= expo_days_per_week)) *
      (1 - (time >= expo_duration_days * 24))

    # ---- dietary input (ug/h) -------------------------------------------
    # Either mg/day directly, or ppm x intake factor x BW^0.75. The ppm arm
    # needs f_intake_diet, which the paper does not report.
    Rdiet <- (diet_mn_mgd + diet_mn_ppm * f_intake_diet * bw^0.75) * 1000 / 24
    # Inhalation (ug/h). DINH in mg/m3 is numerically ug/L of air.
    Rinh_lung <- vent_alv * conc_air_mgm3 * f_dep_lung * expo_on
    Rinh_nasal <- vent_alv * conc_air_mgm3 * f_dep_nasal_respiratory * expo_on
    Rinh_olf <- vent_alv * conc_air_mgm3 * f_dep_nasal_olfactory * expo_on

    # Biliary clearance with Hill-type induction by arterial Mn (ug/g).
    # Bulk and tracer arterial Mn both drive induction.
    cart_ugg <- (c_art + t_c_art) / 1000
    cl_bile_ind <- cl_bile_abs *
      (1 + emax_bile * cart_ugg^hill_bile / (km_bile^hill_bile + cart_ugg^hill_bile))
    Rbile <- cl_bile_ind * c_liver
    t_Rbile <- cl_bile_ind * t_c_liver
    # Dietary Mn and reabsorbed biliary Mn both enter the liver.
    Rabs_liver <- (Rdiet * diet_on + Rbile) * f_diet_uptake
    t_Rabs_liver <- t_Rbile * f_diet_uptake

    # =====================================================================
    # Bulk Mn ODEs (supplement, "Model differential Equations")
    # =====================================================================
    # ---- lung ------------------------------------------------------------
    d/dt(depot_lung) <- Rinh_lung - kdep_lu * depot_lung - kshal_lu * depot_lung
    d/dt(lung_vas) <- qc * (c_ven - c_lung_vas) + kout_lu * lung - kin_lu * lung_vas
    d/dt(lung) <- kin_lu * lung_vas - kout_lu * lung +
      kdr_lu * bound_lung - ka_bind * bf_lu * c_lung + kdep_lu * depot_lung
    d/dt(bound_lung) <- -kdr_lu * bound_lung + ka_bind * bf_lu * c_lung
    # ---- nose ------------------------------------------------------------
    d/dt(depot_nasal) <- Rinh_nasal - kdep_nr * depot_nasal
    d/dt(nasal_respiratory) <- q_nose * (c_art - c_nasal_vas) + kdep_nr * depot_nasal
    d/dt(depot_brain) <- Rinh_olf - knp_ob * depot_brain
    # ---- liver -----------------------------------------------------------
    d/dt(liver_vas) <- q_liver * (c_art - c_liver_vas) + kout_li * liver - kin_li * liver_vas
    d/dt(liver) <- kin_li * liver_vas - kout_li * liver +
      kdr_li * bound_liver - ka_bind * bf_li * c_liver - Rbile + Rabs_liver
    d/dt(bound_liver) <- -kdr_li * bound_liver + ka_bind * bf_li * c_liver
    # ---- gut -------------------------------------------------------------
    d/dt(gut_lumen) <- Rdiet * diet_on * (1 - f_diet_uptake) - kgi * gut_lumen
    d/dt(enterocyte) <- f_enterocyte * kgi * gut_lumen - kent * enterocyte
    d/dt(gut_lumen_lower) <- (1 - f_enterocyte) * kgi * gut_lumen -
      kfeces * gut_lumen_lower + kent * enterocyte + Rbile * (1 - f_diet_uptake)
    d/dt(a_feces) <- kfeces * gut_lumen_lower
    # ---- brain -----------------------------------------------------------
    d/dt(brain_vascular) <- q_brain * (c_art - c_brain_vas) +
      kout_gp * brain_globus_pallidus + kout_ob * brain_olfactory_bulb +
      kout_cb * brain_cerebellum + kout_pt * pituitary -
      (kin_gp + kin_ob + kin_cb + kin_pt) * brain_vascular
    d/dt(brain_globus_pallidus) <- kin_gp * brain_vascular - kout_gp * brain_globus_pallidus +
      kdr_gp * bound_brain_globus_pallidus - ka_bind * bf_gp * c_gp
    d/dt(bound_brain_globus_pallidus) <- -kdr_gp * bound_brain_globus_pallidus +
      ka_bind * bf_gp * c_gp
    # Olfactory bulb uniquely also receives the direct nose-to-brain input.
    d/dt(brain_olfactory_bulb) <- kin_ob * brain_vascular - kout_ob * brain_olfactory_bulb +
      kdr_ob * bound_brain_olfactory_bulb - ka_bind * bf_ob * c_ob + knp_ob * depot_brain
    d/dt(bound_brain_olfactory_bulb) <- -kdr_ob * bound_brain_olfactory_bulb +
      ka_bind * bf_ob * c_ob
    d/dt(brain_cerebellum) <- kin_cb * brain_vascular - kout_cb * brain_cerebellum +
      kdr_cb * bound_brain_cerebellum - ka_bind * bf_cb * c_cb
    d/dt(bound_brain_cerebellum) <- -kdr_cb * bound_brain_cerebellum + ka_bind * bf_cb * c_cb
    d/dt(pituitary) <- kin_pt * brain_vascular - kout_pt * pituitary +
      kdr_pt * bound_pituitary - ka_bind * bf_pt * c_pit
    d/dt(bound_pituitary) <- -kdr_pt * bound_pituitary + ka_bind * bf_pt * c_pit
    # ---- bone ------------------------------------------------------------
    d/dt(bone_vas) <- q_bone * (c_art - c_bone_vas) + kout_bo * bone - kin_bo * bone_vas
    d/dt(bone) <- kin_bo * bone_vas - kout_bo * bone +
      kdr_bo * bound_bone - ka_bind * bf_bo * c_bone
    d/dt(bound_bone) <- -kdr_bo * bound_bone + ka_bind * bf_bo * c_bone
    # ---- rest of body ----------------------------------------------------
    d/dt(other_vas) <- q_other * (c_art - c_other_vas) + kout_ot * other - kin_ot * other_vas
    d/dt(other) <- kin_ot * other_vas - kout_ot * other +
      kdr_ot * bound_other - ka_bind * bf_ot * c_other
    d/dt(bound_other) <- -kdr_ot * bound_other + ka_bind * bf_ot * c_other
    # ---- blood -----------------------------------------------------------
    d/dt(arterial) <- qc * c_lung_vas - qc * c_art + kshal_lu * depot_lung
    d/dt(venous) <- q_liver * c_liver_vas + q_brain * c_brain_vas +
      q_other * c_other_vas + q_bone * c_bone_vas + q_nose * c_nasal_vas -
      qc * c_ven + ksc * depot_sc
    d/dt(depot_sc) <- -ksc * depot_sc
    # ---- cumulative intake, for the mass-balance check -------------------
    d/dt(a_oral) <- Rdiet * diet_on
    d/dt(a_inhaled) <- Rinh_lung + Rinh_nasal + Rinh_olf

    # =====================================================================
    # 54Mn tracer ODEs (supplement, "Tracer Kinetics"). Structurally
    # identical, sharing every parameter and the free-site terms above.
    # =====================================================================
    d/dt(depot_lung_mn54) <- -kdep_lu * depot_lung_mn54 - kshal_lu * depot_lung_mn54
    d/dt(lung_vas_mn54) <- qc * (t_c_ven - t_c_lung_vas) + kout_lu * lung_mn54 -
      kin_lu * lung_vas_mn54
    d/dt(lung_mn54) <- kin_lu * lung_vas_mn54 - kout_lu * lung_mn54 +
      kdr_lu * bound_lung_mn54 - ka_bind * bf_lu * t_c_lung + kdep_lu * depot_lung_mn54
    d/dt(bound_lung_mn54) <- -kdr_lu * bound_lung_mn54 + ka_bind * bf_lu * t_c_lung
    d/dt(nasal_respiratory_mn54) <- q_nose * (t_c_art - t_c_nasal_vas)
    d/dt(depot_ip_mn54) <- -kip * depot_ip_mn54
    d/dt(liver_vas_mn54) <- q_liver * (t_c_art - t_c_liver_vas) + kout_li * liver_mn54 -
      kin_li * liver_vas_mn54
    d/dt(liver_mn54) <- kin_li * liver_vas_mn54 - kout_li * liver_mn54 +
      kdr_li * bound_liver_mn54 - ka_bind * bf_li * t_c_liver - t_Rbile +
      t_Rabs_liver + kip * depot_ip_mn54
    d/dt(bound_liver_mn54) <- -kdr_li * bound_liver_mn54 + ka_bind * bf_li * t_c_liver
    d/dt(gut_lumen_mn54) <- -kgi * gut_lumen_mn54
    d/dt(enterocyte_mn54) <- f_enterocyte * kgi * gut_lumen_mn54 - kent * enterocyte_mn54
    d/dt(gut_lumen_lower_mn54) <- (1 - f_enterocyte) * kgi * gut_lumen_mn54 -
      kfeces * gut_lumen_lower_mn54 + kent * enterocyte_mn54 + t_Rbile * (1 - f_diet_uptake)
    d/dt(a_feces_mn54) <- kfeces * gut_lumen_lower_mn54
    d/dt(brain_vascular_mn54) <- q_brain * (t_c_art - t_c_brain_vas) +
      kout_gp * brain_globus_pallidus_mn54 + kout_ob * brain_olfactory_bulb_mn54 +
      kout_cb * brain_cerebellum_mn54 + kout_pt * pituitary_mn54 -
      (kin_gp + kin_ob + kin_cb + kin_pt) * brain_vascular_mn54
    d/dt(brain_globus_pallidus_mn54) <- kin_gp * brain_vascular_mn54 -
      kout_gp * brain_globus_pallidus_mn54 +
      kdr_gp * bound_brain_globus_pallidus_mn54 - ka_bind * bf_gp * t_c_gp
    d/dt(bound_brain_globus_pallidus_mn54) <- -kdr_gp * bound_brain_globus_pallidus_mn54 +
      ka_bind * bf_gp * t_c_gp
    d/dt(brain_olfactory_bulb_mn54) <- kin_ob * brain_vascular_mn54 -
      kout_ob * brain_olfactory_bulb_mn54 +
      kdr_ob * bound_brain_olfactory_bulb_mn54 - ka_bind * bf_ob * t_c_ob
    d/dt(bound_brain_olfactory_bulb_mn54) <- -kdr_ob * bound_brain_olfactory_bulb_mn54 +
      ka_bind * bf_ob * t_c_ob
    d/dt(brain_cerebellum_mn54) <- kin_cb * brain_vascular_mn54 -
      kout_cb * brain_cerebellum_mn54 +
      kdr_cb * bound_brain_cerebellum_mn54 - ka_bind * bf_cb * t_c_cb
    d/dt(bound_brain_cerebellum_mn54) <- -kdr_cb * bound_brain_cerebellum_mn54 +
      ka_bind * bf_cb * t_c_cb
    d/dt(pituitary_mn54) <- kin_pt * brain_vascular_mn54 - kout_pt * pituitary_mn54 +
      kdr_pt * bound_pituitary_mn54 - ka_bind * bf_pt * t_c_pit
    d/dt(bound_pituitary_mn54) <- -kdr_pt * bound_pituitary_mn54 + ka_bind * bf_pt * t_c_pit
    d/dt(bone_vas_mn54) <- q_bone * (t_c_art - t_c_bone_vas) + kout_bo * bone_mn54 -
      kin_bo * bone_vas_mn54
    d/dt(bone_mn54) <- kin_bo * bone_vas_mn54 - kout_bo * bone_mn54 +
      kdr_bo * bound_bone_mn54 - ka_bind * bf_bo * t_c_bone
    d/dt(bound_bone_mn54) <- -kdr_bo * bound_bone_mn54 + ka_bind * bf_bo * t_c_bone
    d/dt(other_vas_mn54) <- q_other * (t_c_art - t_c_other_vas) + kout_ot * other_mn54 -
      kin_ot * other_vas_mn54
    d/dt(other_mn54) <- kin_ot * other_vas_mn54 - kout_ot * other_mn54 +
      kdr_ot * bound_other_mn54 - ka_bind * bf_ot * t_c_other
    d/dt(bound_other_mn54) <- -kdr_ot * bound_other_mn54 + ka_bind * bf_ot * t_c_other
    d/dt(arterial_mn54) <- qc * t_c_lung_vas - qc * t_c_art + kshal_lu * depot_lung_mn54
    d/dt(venous_mn54) <- q_liver * t_c_liver_vas + q_brain * t_c_brain_vas +
      q_other * t_c_other_vas + q_bone * t_c_bone_vas + q_nose * t_c_nasal_vas -
      qc * t_c_ven + ksc * depot_sc_mn54
    d/dt(depot_sc_mn54) <- -ksc * depot_sc_mn54

    # =====================================================================
    # Observations (supplement, CalcOutputs). Tissue concentrations are
    # reported as ug Mn per g tissue, the units of Figures 3 and 6.
    # =====================================================================
    # Cc is whole-blood Mn in ug/L -- the measured circulating quantity.
    Cc <- c_art
    # Freeland-Graves comparison (Figure 10): "Simulated plasma Mn
    # concentrations were calculated as 20% of simulated blood Mn."
    Cplasma <- c_art * 0.20
    Ctot_liver <- (liver + liver_vas + bound_liver) / (v_liver * 1000)
    Ctot_lung <- (lung + lung_vas + bound_lung + depot_lung) / (v_lung * 1000)
    Ctot_bone <- (bone + bone_vas + bound_bone) / (v_bone * 1000)
    Ctot_other <- (other + other_vas + bound_other) / (v_other * 1000)
    Ctot_nose <- (nasal_respiratory + depot_nasal + depot_brain) / (v_nose * 1000)
    # Brain regions carry a share of the brain vascular pool, apportioned
    # by region volume fraction (supplement: ast = ABRST + BMNST + ABRBD*VSTMC).
    Ctot_globus_pallidus <- (brain_globus_pallidus + bound_brain_globus_pallidus +
      brain_vascular * f_globus_pallidus) / v_gp / 1000
    Ctot_olfactory_bulb <- (brain_olfactory_bulb + bound_brain_olfactory_bulb +
      brain_vascular * f_olfactory_bulb) / v_ob / 1000
    Ctot_cerebellum <- (brain_cerebellum + bound_brain_cerebellum +
      brain_vascular * f_cerebellum) / v_cb / 1000
    Ctot_pituitary <- (pituitary + bound_pituitary +
      brain_vascular * f_pituitary) / v_pit / 1000
    Cbile <- Rbile / q_bile_abs / 1000
    # Percent of maximal binding capacity still free, per tissue -- the
    # saturation diagnostic the paper's mechanism turns on.
    Pfree_globus_pallidus <- bf_gp / bcap_gp * 100
    Pfree_liver <- bf_li / bcap_li * 100

    # ---- mass balance (supplement `balance`), an author-provided gate ----
    body_mn <- (liver + liver_vas + bound_liver) +
      (brain_globus_pallidus + brain_olfactory_bulb + brain_cerebellum + pituitary +
         bound_brain_globus_pallidus + bound_brain_olfactory_bulb +
         bound_brain_cerebellum + bound_pituitary + brain_vascular) +
      (other + other_vas + bound_other) + arterial + venous +
      (bone + bone_vas + bound_bone) + (lung + lung_vas + bound_lung + depot_lung) +
      (nasal_respiratory + depot_nasal + depot_brain) +
      (gut_lumen + enterocyte + gut_lumen_lower) + depot_sc
    Mn_balance <- (body_mn + a_feces) / (a_oral + a_inhaled + 1e-30)

    # ---- 54Mn whole-body retention (supplement `ret54`) ------------------
    body_mn54 <- (liver_mn54 + liver_vas_mn54 + bound_liver_mn54) +
      (brain_globus_pallidus_mn54 + brain_olfactory_bulb_mn54 + brain_cerebellum_mn54 +
         pituitary_mn54 + bound_brain_globus_pallidus_mn54 +
         bound_brain_olfactory_bulb_mn54 + bound_brain_cerebellum_mn54 +
         bound_pituitary_mn54 + brain_vascular_mn54) +
      (other_mn54 + other_vas_mn54 + bound_other_mn54) + arterial_mn54 + venous_mn54 +
      (bone_mn54 + bone_vas_mn54 + bound_bone_mn54) +
      (lung_mn54 + lung_vas_mn54 + bound_lung_mn54 + depot_lung_mn54) +
      nasal_respiratory_mn54 +
      (gut_lumen_mn54 + enterocyte_mn54 + gut_lumen_lower_mn54) +
      depot_ip_mn54 + depot_sc_mn54
    Retention_mn54 <- body_mn54 / dose_mn54 * 100
    # Tracer mass balance. The supplement writes dt(XLOSS) = XFECES, which
    # integrates the already-cumulative faecal state a second time and so
    # cannot satisfy the authors' own stated diagnostic ("value close to 1
    # means the model is mass balanced for 54Mn"). The cold-side mirror is
    # dt(LOSS) = RFECES, a rate; cumulative faecal tracer is therefore the
    # intended term and is what is used here. See the vignette Errata.
    Mn54_balance <- (body_mn54 + a_feces_mn54) / dose_mn54
    # No residual-error or IIV model: Campbell 2023 is a deterministic
    # mechanistic PBPK and reports neither.
  })
}
