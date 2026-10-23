Bihorel_2026_aflibercept <- function() {
  description <- "Integrated semi-mechanistic population PK model for free and VEGF-bound (adjusted bound) aflibercept in plasma after intravitreal (AFL-2 2 mg and AFL-HD 8 mg drug products), intravenous and subcutaneous administration (Bihorel 2026). Study and fellow eyes (fixed 4 mL vitreous volumes) exchange free aflibercept bidirectionally with the central compartment at the ocular clearance QE; free aflibercept distributes to two peripheral compartments, is eliminated first-order (K20), is taken up by a capacity-limited platelet pool, and binds VEGF in the central compartment by a Michaelis-Menten process that feeds the bound-aflibercept compartment (equal volume to central, first-order elimination K40). SC doses are absorbed first-order. Allometric body weight on all volumes and on K20/K40, albumin on K40, age and drug product (AFL-HD 39.4% slower) on QE, and healthy-volunteer / oncology indicators on Vmax. Pooled 2744 participants from 16 trials (healthy volunteers, nAMD, DME, solid tumors / lymphoma)."
  reference <- "Bihorel S, Chiu J, Chittenden J, Eissing T, Turner KC, Hochel J, Zimmermann T, He ZW, Davis JD, DiCioccio AT, Harnisch LO. Integrated Population Pharmacokinetic Analysis of Intravitreal Aflibercept Drug Products. CPT Pharmacometrics Syst Pharmacol. 2026. doi:10.1002/psp4.70287. Adapts the IV-only free/bound aflibercept model of Thai et al. (2011); see modellib('Thai_2011_aflibercept')."
  vignette <- "Bihorel_2026_aflibercept"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power scaling (WT/76) on V2 = V4, V3, V8, K20 and K40 with estimated exponents; 76 kg is the reference patient of the Figure 4 caption.",
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE/71)^-1.36 on the ocular clearance QE; 71 years is the reference patient of the Figure 4 caption.",
      source_name = "age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. The paper reports albumin in g/dL with reference 4.3 g/dL; the model converts the canonical g/L column inline (alb_gdL <- ALB * 0.1) before applying (ALB/4.3)^-0.767 on K40.",
      source_name = "ALB"
    ),
    FORM_AFL_HD = list(
      description = "Aflibercept 8 mg drug product (AFL-HD; 70 uL of a 114.3 mg/mL formulation) vs the 2 mg drug product (AFL-2; 50 uL of a 40 mg/mL formulation)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (AFL-2)",
      notes = "Paper symbol HD. Multiplies QE by 0.606^HD (39.4% slower ocular clearance for AFL-HD). The model carries one QE per subject, so the indicator describes the product injected in the study eye; set 0 for IV, SC and AFL-2 IVT dosing.",
      source_name = "HD"
    ),
    DIS_HEALTHY = list(
      description = "Healthy participant",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient: nAMD, DME, or solid tumor / lymphoma)",
      notes = "Paper symbol HV. Multiplies Vmax by 0.613^HV. The reference Vmax is the combined nAMD + DME population.",
      source_name = "HV"
    ),
    DIS_CANCER = list(
      description = "Participant with a solid tumor or lymphoma",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant, nAMD or DME)",
      notes = "Paper symbol ONCO ('1 in participants with tumors'). The oncology cohort (n = 38) pooled solid tumors and lymphoma; mapped to DIS_CANCER because the cohort is predominantly the advanced-solid-tumor IV studies of the upstream model. Multiplies Vmax by 1.24^ONCO.",
      source_name = "ONCO"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Listed among the evaluated covariates (Abstract) but not retained in the final model. Racial classification, Japanese origin and renal and hepatic function were likewise evaluated and assessed only as exposure subgroups (Figure S5)."
    )
  )

  compartmentData <- list(
    vitreous = list(analyte = "free aflibercept", units = "mg", specimen = "vitreous", verified = TRUE),
    central = list(analyte = "free aflibercept", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "free aflibercept", units = "mg", specimen = "tissue", verified = TRUE),
    complex = list(
      analyte = "VEGF:aflibercept complex (adjusted bound aflibercept, aflibercept-mass equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    vitreous_fellow = list(analyte = "free aflibercept", units = "mg", specimen = "vitreous", verified = TRUE),
    depot = list(analyte = "free aflibercept", units = "mg", specimen = "administration site", verified = TRUE),
    plt_afl = list(analyte = "platelet-associated aflibercept", units = "mg", specimen = "blood cell", verified = TRUE),
    peripheral2 = list(analyte = "free aflibercept", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2744L,
    n_studies = 16L,
    age_range = NA_character_,
    age_median = "71 years (reference patient, Figure 4 caption)",
    weight_range = NA_character_,
    weight_median = "76 kg (reference patient, Figure 4 caption)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Healthy participants (n = 76), neovascular age-related macular degeneration (nAMD, n = 1662), diabetic macular edema (DME, n = 968), and solid tumors or lymphoma (n = 38).",
    dose_range = "IVT 0.05-4 mg (AFL-2) or 8 mg (AFL-HD); IV 0.3-4 mg/kg; SC 0.025-2 mg/kg.",
    regions = "Multinational (16 Regeneron / Bayer trials including CANDELA, PULSAR and PHOTON).",
    notes = "31,326 plasma samples (15,657 free aflibercept, 62.9% BLQ; 15,669 adjusted bound aflibercept, 19.6% BLQ) analysed with Beal's M3 method in NONMEM 7.5. Baseline demographics are in Table S5 (supplement not on file); the age and weight medians above are the Figure 4 reference patient. Adjusted bound aflibercept = measured VEGF:aflibercept complex x 0.717 (aflibercept mass fraction of the hVEGF165:aflibercept complex)."
  )

  ini({
    # Ocular disposition (Table 1, stage 4 estimates)
    lf_vitreous <- log(0.701); label("Bioavailability of an intravitreal dose in either eye, F_IVT (fraction)") # Table 1 'F IVT' = 0.701
    lq_vitreous <- log(0.000608); label("Ocular clearance QE between each eye and central, AFL-2 (L/day)") # Table 1 'QE for AFL-2' = 0.608 mL/day
    e_form_afl_hd_q_vitreous <- 0.606; label("Proportional effect of AFL-HD drug product on QE (unitless, 0.606^HD)") # Table 1 = 0.606
    e_age_q_vitreous <- -1.36; label("Power exponent of age on QE (unitless, reference 71 years)") # Table 1 = -1.36; Section 3.1 equation
    lv_vitreous <- fixed(log(0.004)); label("Volume of study eye V1 and fellow eye V5 (L)") # Table 1 'V1 and V5' = 4 mL, fixed

    # Free aflibercept systemic disposition (Table 1, stage 1 or stage 3 estimates)
    lkel <- log(0.0807); label("First-order elimination rate constant of free aflibercept K20 (1/day)") # Table 1 'K20' = 0.0807
    e_wt_kel <- -0.192; label("Power exponent of body weight on K20 (unitless, reference 76 kg)") # Table 1 = -0.192
    lvc <- log(4.99); label("Central volume V2 (= bound-aflibercept volume V4) (L)") # Table 1 'V2' = 4.99 (V4 set equal)
    e_wt_vc <- 0.872; label("Power exponent of body weight on V2 and V4 (unitless)") # Table 1 = 0.872
    lvp <- log(1.08); label("First peripheral volume V3 (L)") # Table 1 'V3' = 1.08
    e_wt_vp <- 1.08; label("Power exponent of body weight on V3 (unitless)") # Table 1 = 1.08
    lvp2 <- log(1.18); label("Second peripheral volume V8 (L)") # Table 1 'V8' = 1.18
    e_wt_vp2 <- 1.16; label("Power exponent of body weight on V8 (unitless)") # Table 1 = 1.16
    lq <- log(0.849); label("First distribution clearance QF1 (L/day)") # Table 1 'QF1' = 0.849
    lq2 <- log(0.0763); label("Second distribution clearance QF2 (L/day)") # Table 1 'QF2' = 0.0763

    # Platelet uptake (Table 1, stage 1 estimates)
    lvmax_plt <- log(0.0310); label("Maximum platelet uptake rate Vmax,p (1/day/L)") # Table 1 'Vmax,p' = 0.0310
    lkm_plt <- log(42.7); label("Platelet-pool amount at half-maximal uptake Km,p (mg)") # Table 1 'Km,p' = 42.7
    lkel_plt <- log(0.0265); label("Elimination rate constant from the platelet pool K70 (1/day)") # Table 1 'K70' = 0.0265

    # Subcutaneous absorption (Table 1, stage 1 estimates)
    logitfdepot <- log(0.536 / (1 - 0.536)); label("Logit of SC bioavailability F_SC (logit scale, unitless)") # Table 1 'F SC' = 0.536; IIV on the logit scale (Table 1 note)
    lka <- log(0.368); label("SC absorption rate constant K62 (1/day)") # Table 1 'K62' = 0.368

    # VEGF binding and bound aflibercept
    lvmax <- log(0.224); label("Maximum binding rate to VEGF in nAMD or DME, Vmax (mg/day/L)") # Table 1 'Vmax' = 0.224 (95% CI 0.213, 0.227); the Section 3.1 equation prints 0.0224, see vignette
    e_dis_healthy_vmax <- 0.613; label("Proportional effect of healthy participant on Vmax (unitless, 0.613^HV)") # Table 1 = 0.613
    e_dis_cancer_vmax <- 1.24; label("Proportional effect of tumor participant on Vmax (unitless, 1.24^ONCO)") # Table 1 = 1.24
    lkm <- log(0.411); label("Free aflibercept concentration at half Vmax, Km (mg/L)") # Table 1 'Km' = 0.411
    lkint <- log(0.0350); label("First-order elimination rate constant of bound aflibercept K40 (1/day)") # Table 1 'K40' = 0.0350
    e_wt_kint <- -0.153; label("Power exponent of body weight on K40 (unitless)") # Table 1 = -0.153
    e_alb_kint <- -0.767; label("Power exponent of serum albumin on K40 (unitless, reference 4.3 g/dL)") # Table 1 = -0.767

    # IIV: variances on the log scale (logit scale for F_SC), Table 1 note
    etalq_vitreous ~ 0.597 # Table 1 'omega2 in QE' = 0.597 (CV 90.4%)
    etalkel + etalvc ~ c(0.207, -0.0727, 0.0618) # Table 1 'omega2 in K20' 0.207, covariance -0.0727, 'omega2 in V2 and V4' 0.0618
    etalvmax ~ 0.158 # Table 1 'omega2 in Vmax' = 0.158 (CV 41.3%)
    etalogitfdepot ~ 0.629 # Table 1 'omega2 in FSC' = 0.629 (logit scale)
    etalka ~ 0.852 # Table 1 'omega2 in K62' = 0.852 (CV 116%)
    etalkint ~ 0.0452 # Table 1 'omega2 in K40' = 0.0452 (CV 21.5%)

    # Residual error (SDs) after IVT dosing at the LLOQs of the pivotal
    # CANDELA / PULSAR / PHOTON assays (free 0.0156 mg/L, bound 0.0224 mg/L;
    # Figure 3 caption). Table 1 also reports IV/SC and other-LLOQ sets:
    # free IV/SC LLOQ 0.0313: add 0.0250, prop 0.403; free IV/SC LLOQ 0.0156:
    # add 0.00786, prop 0.357; bound IV+SC LLOQ 0.0315: add 0.0206, prop 0.167;
    # bound IVT LLOQ 0.0315: add 0.0214, prop 0.161.
    addSd <- 0.00644; label("Additive residual SD, free aflibercept after IVT (mg/L)") # Table 1 'sigma additive after IVT dosing when LLOQ is 0.0156' = 0.00644
    propSd <- 0.415; label("Proportional residual SD, free aflibercept after IVT (fraction)") # Table 1 'sigma proportional after IVT dosing when LLOQ is 0.0156' = 0.415
    addSd_complex <- 0.0180; label("Additive residual SD, adjusted bound aflibercept after IVT (mg/L)") # Table 1 'sigma additive after IVT dosing when LLOQ is 0.0224' = 0.0180
    propSd_complex <- 0.227; label("Proportional residual SD, adjusted bound aflibercept after IVT (fraction)") # Table 1 'sigma proportional after IVT dosing when LLOQ is 0.0224' = 0.227
  })

  model({
    # Covariate terms (Section 3.1 typical-value equations)
    alb_gdL <- ALB * 0.1 # SI g/L -> US-convention g/dL used by the paper
    wt_ratio <- WT / 76

    # Individual parameters
    f_vitreous <- exp(lf_vitreous)
    q_vitreous <- exp(lq_vitreous + etalq_vitreous) *
      e_form_afl_hd_q_vitreous^FORM_AFL_HD *
      (AGE / 71)^e_age_q_vitreous
    v_vitreous <- exp(lv_vitreous)
    kel <- exp(lkel + etalkel) * wt_ratio^e_wt_kel
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    vp <- exp(lvp) * wt_ratio^e_wt_vp
    vp2 <- exp(lvp2) * wt_ratio^e_wt_vp2
    q <- exp(lq)
    q2 <- exp(lq2)
    vmax_plt <- exp(lvmax_plt)
    km_plt <- exp(lkm_plt)
    kel_plt <- exp(lkel_plt)
    logit_fdepot <- logitfdepot + etalogitfdepot
    fdepot <- expit(logit_fdepot)
    ka <- exp(lka + etalka)
    vmax <- exp(lvmax + etalvmax) *
      e_dis_healthy_vmax^DIS_HEALTHY *
      e_dis_cancer_vmax^DIS_CANCER
    km <- exp(lkm)
    kint <- exp(lkint + etalkint) * wt_ratio^e_wt_kint * (alb_gdL / 4.3)^e_alb_kint
    v_complex <- vc # V4 = V2 (Figure 1 caption)

    # Saturable processes (Figure 1 ODEs). VEGF binding is Michaelis-Menten
    # in the free central concentration; platelet uptake is first order in
    # the free central amount and shuts off as the platelet pool fills.
    cfree <- central / vc
    binding <- vmax * vc * cfree / (km + cfree)
    uptake_plt <- vmax_plt * vc * (1 - plt_afl / (km_plt + plt_afl)) * central

    # ODE system, states in the paper's A1..A8 order (Figure 1 caption)
    d/dt(vitreous) <- q_vitreous / vc * central - q_vitreous / v_vitreous * vitreous
    d/dt(central) <- q_vitreous / v_vitreous * vitreous +
      q_vitreous / v_vitreous * vitreous_fellow +
      q / vp * peripheral1 + ka * depot + q2 / vp2 * peripheral2 -
      ((2 * q_vitreous + q + q2) / vc + kel) * central -
      binding - uptake_plt
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(complex) <- binding - kint * complex
    d/dt(vitreous_fellow) <- q_vitreous / vc * central - q_vitreous / v_vitreous * vitreous_fellow
    d/dt(depot) <- -ka * depot
    d/dt(plt_afl) <- uptake_plt - kel_plt * plt_afl
    d/dt(peripheral2) <- q2 / vc * central - q2 / vp2 * peripheral2

    # Bioavailability: F_IVT on doses into either eye, F_SC on SC doses
    # (Figure 1 diagram). IV doses go straight to central.
    f(vitreous) <- f_vitreous
    f(vitreous_fellow) <- f_vitreous
    f(depot) <- fdepot

    # Free aflibercept in the study eye (not an observation in the source)
    C_vitreous <- vitreous / v_vitreous

    Cc <- cfree
    Cc_complex <- complex / v_complex

    Cc ~ add(addSd) + prop(propSd)
    Cc_complex ~ add(addSd_complex) + prop(propSd_complex)
  })
}
