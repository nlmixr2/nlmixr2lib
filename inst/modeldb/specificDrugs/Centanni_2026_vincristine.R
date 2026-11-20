Centanni_2026_vincristine <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous vincristine in",
    "Dutch and Kenyan children with cancer, linked to a proportional-odds",
    "Markov model for vincristine-induced peripheral neuropathy (VIPN) with",
    "a separate logistic model for informative censoring (Centanni 2026,",
    "final models of Table 2, estimated on the pooled Dutch and Kenyan data).",
    "Linear elimination and linear distribution to one peripheral",
    "compartment, plus saturable (Michaelis-Menten) uptake into a binding",
    "compartment with first-order return; allometric weight scaling",
    "(exponents fixed at 0.75 and 1, 70 kg reference); a shared exponential",
    "Kenyan effect on CL and Vp; an infusion-only lag time; log-additive",
    "residual error with separate magnitudes for bolus and infusion samples",
    "and IIV on the infusion residual. The neuropathy model gives the",
    "probability of CTCAE grade >= 1 and >= 2 for each of four subscores",
    "(sensory, motor, constipation, neuralgia) driven by an Emax function of",
    "the occasion AUC = dose / CL, with a Markov term for a previous subscore",
    "above grade 0 and a Kenyan Emax fixed to effectively zero; the",
    "censoring model gives the per-assessment probability of dropping out,",
    "rising with the dosing occasion and after a previous grade >= 2."
  )
  reference <- paste(
    "Centanni M, van de Velde ME, Huitema ADR, Kaspers GJL, Karlsson MO,",
    "Friberg LE. Extrapolating Vincristine-Induced Peripheral Neuropathy From",
    "Caucasian to Kenyan Populations: Impact of Type I and Type II Selection",
    "Bias. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70303.",
    "doi:10.1002/psp4.70303"
  )
  vignette <- "Centanni_2026_vincristine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling with fixed exponents 0.75 on CL and Q and 1 on",
        "Vc and Vp, normalised to 70 kg (Methods 'Pharmacokinetic Model';",
        "Table 2 footnote a; supplement 2 $PK). The binding parameters are",
        "not weight-scaled."
      ),
      source_name = "WT"
    ),
    TINF = list(
      description = "Duration of the intravenous vincristine infusion; 0 for a bolus",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the administration method and sets the modelled infusion",
        "duration. TINF > 0 marks an infusion: the dose enters after the",
        "infusion lag time and over TINF hours (supplement 2 $PK",
        "'D1 = INFT' and 'ALAG1 = THETA(6)*EXP(ETA(3))' when DOSET >= 1),",
        "and the sample takes the infusion residual error with its IIV.",
        "TINF = 0 marks a bolus (no lag, bolus residual error); boluses",
        "lasting 2-5 min were modelled as instantaneous. Infusion dose",
        "records must carry rate = -2 so that rxode2 takes the duration from",
        "the model; supply TINF as a data column because rxode2 does not",
        "expose a dose record's duration to model(). The Dutch protocol",
        "infusion was 60 min (38 min of vincristine, 22 min of flushing),",
        "96 min in four and 15 min in two patients (Supporting Methods S1).",
        "The paper does not say whether INFT was the total or the",
        "vincristine-only time; all Kenyan doses were boluses."
      ),
      source_name = "INFT (with DOSET >= 1 marking an infusion)"
    ),
    REGION_KENYA = list(
      description = "Kenyan cohort indicator (1 = Kenyan child, 0 = Dutch child)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper's nationality covariate NAT (0 = Dutch, 1 = Kenyan).",
        "Multiplies CL and Vp by exp(-0.786) (one shared coefficient,",
        "Equation 4 and Table 2) and replaces Emax of the neuropathy model",
        "by exp(-100 x 20), which is effectively zero (Table 2 'Kenyan Emax",
        "parameter'). The Kenyan children were treated in Kenya and the",
        "Dutch children in the Netherlands, so nationality, enrollment",
        "country and study coincide."
      ),
      source_name = "NAT"
    ),
    CYCLE = list(
      description = "Vincristine dose occasion of the neuropathy assessment (0 = baseline before the first dose, 1 = after the first dose, ...)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Enters only the censoring model, as the step indicators",
        "CYCLE > 0, CYCLE > 1 and CYCLE > 2 (supplement 2 $PRED",
        "'OCCT1/OCCT2/OCCT3'; Table 2 'Censoring | OCC > 0/1/2'). The",
        "assessment after the first dose takes CYCLE = 1; the baseline",
        "assessment before any dose takes CYCLE = 0, where the censoring",
        "logit is the baseline -144.7. The VIPN logits do not depend on the",
        "occasion."
      ),
      source_name = "OCC"
    ),
    PREV_AE_SCORE = list(
      description = "CTCAE grade (0, 1, >= 2 coded 2) of the same VIPN subscore at the previous assessment",
      units = "(ordinal score)",
      type = "count",
      reference_category = 0,
      notes = paste(
        "Markov state on a three-level scale: CTCAE 0, 1 and >= 2 (Methods",
        "'Pharmacodynamic Model'). Only PREV_AE_SCORE > 0 versus 0 enters",
        "the VIPN logits (Table 2 'B1|0, if Grade > 0'; supplement 2",
        "'PREVk = 1 if PDVk > 0'). Each of the four subscores (peripheral",
        "sensory neuropathy, peripheral motor neuropathy, constipation,",
        "neuralgia) is its own chain, so the value on a row is the previous",
        "grade of that row's subscore; 0 at the first assessment. The",
        "censoring model instead conditions on any subscore having been",
        ">= 2 at the previous assessment, which spans the four chains and",
        "so is returned as two outputs (p_cens_vipn and p_cens_prev_vipn)",
        "rather than read from this column."
      ),
      source_name = "PDV1-PDV4 (previous score per subscore)"
    )
  )

  covariatesDataExcluded <- list(
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = paste(
        "Tested on clearance (Equation 1, linear between the observed",
        "extremes 22 and 862) and on vincristine binding, with values",
        "measured in the Kenyan children only; not retained (Results",
        "'Pharmacokinetic Model')."
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "vincristine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vincristine", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "vincristine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 2L,
    n_observations = "513 vincristine plasma concentrations; 800 neuropathy assessments",
    age_range = "2-17 years (Dutch 2-17, Kenyan 5-14)",
    disease_state = paste(
      "Pediatric cancer treated with vincristine: acute lymphoblastic",
      "leukemia (64% of the Dutch cohort), Hodgkin lymphoma,",
      "nephroblastoma, medulloblastoma, rhabdomyosarcoma, low-grade glioma",
      "and neuroblastoma"
    ),
    dose_range = paste(
      "Dutch: 1.5 or 2 mg/m^2 (cap 2 mg) as an IV bolus (20 patients) or",
      "infusion (15 patients); Kenyan: 2 mg/m^2 (cap 2.5 mg) as an IV bolus,",
      "with a 20% dose increase after low exposure without neuropathy;",
      "one to five dosing occasions per patient"
    ),
    regions = "Netherlands (35 children) and Kenya (15 children)",
    sampling = paste(
      "Dutch: 10, 20, 30, 40, 60, 75 and 140 min after the bolus or start",
      "of infusion and at 24 h when still hospitalised; Kenyan: 60, 90 and",
      "240 min after the bolus on two occasions; LC-MS/MS, LLOQ 0.1 ng/mL"
    ),
    notes = paste(
      "Pooled re-analysis of the van de Velde 2020 Dutch trial (35",
      "children) and the CHAPATI feasibility study in Kenya (15 children;",
      "Uittenboogaard 2024). Each dosing occasion was treated as a separate",
      "individual when the PK model was estimated, so the IIV on CL also",
      "carries occasion-to-occasion variability. NONMEM 7.5.1, importance",
      "sampling. The paper's demographic table is not part of the article;",
      "only the age ranges are reported."
    )
  )

  ini({
    # PK, Table 2 (allometric scaling to 70 kg, Table 2 footnote a)
    lcl <- log(42.2); label("Clearance at 70 kg, Dutch (L/h)") # Table 2 CL = 42.2 L/h (RSE 14%)
    lvc <- log(14.5); label("Central volume at 70 kg (L)") # Table 2 Vc = 14.5 L (RSE 19%)
    lq <- log(71.6); label("Intercompartmental clearance at 70 kg (L/h)") # Table 2 Q = 71.6 L/h (RSE 9%)
    lvp <- log(541.5); label("Peripheral volume at 70 kg, Dutch (L)") # Table 2 Vp = 541.5 L (RSE 13%)
    e_region_kenya_cl_vp <- -0.786; label("Exponential effect of Kenyan cohort on CL and Vp (unitless)") # Table 2 'Kenyan nationality on CL' and 'on Vp' = -0.786 (RSE 26%), one shared THETA(7) in supplement 2
    e_wt_cl <- fixed(0.75); label("Allometric exponent of weight on CL and Q (unitless)") # Methods 'Pharmacokinetic Model': 0.75 for clearance parameters
    e_wt_vc <- fixed(1); label("Allometric exponent of weight on Vc and Vp (unitless)") # Methods 'Pharmacokinetic Model': 1 for distribution parameters
    ltlag <- fixed(log(0.02)); label("Lag time of an infusion (h)") # supplement 2 $THETA(6) ALAG1 = 0.02 FIX (h); Table 2 prints 'LAG (min) 2 (FIX)'
    lvmax_bind <- log(7.17); label("Maximum rate of saturable uptake into the binding compartment, Bmax (mg/h)") # Table 2 Bmax = 7.17 'mg' (RSE 38%); amount per hour in supplement 2 $DES
    lkm_bind <- log(5.03); label("Central amount giving half-maximal uptake into the binding compartment, Kb (ug)") # Table 2 Kb = 5.03 ug (RSE 44%)
    lk31 <- log(28.3); label("Rate constant out of the binding compartment, kp (1/h)") # Table 2 kp = 28.3 /h (RSE 33%)

    # IIV, Table 2; CV% = 100 x sqrt(exp(omega^2) - 1) (footnote b), so omega^2 = log(1 + CV^2)
    etalcl ~ 0.77705 # Table 2 'IIV CL' = 108.4% -> log(1 + 1.084^2)
    etalvc ~ 1.16386 # Table 2 'IIV Vc' = 148.4% -> log(1 + 1.484^2)
    etaltlag ~ fixed(0.419665) # Table 2 'IIV LAG' = 72.2%, held constant; supplement 2 $OMEGA ETA(3) = 0.419665
    etaRUV ~ 0.46846 # Table 2 'IIV RUVprop, infusion' = 77.3% -> log(1 + 0.773^2)

    # Residual error: additive on log-transformed concentrations (Table 2 footnote c)
    expSd_bolus <- 0.36; label("Log-scale residual SD, bolus samples") # Table 2 'RUVprop,bolus' = 0.36 (RSE 7%)
    expSd_infusion <- 0.38; label("Log-scale residual SD, infusion samples") # Table 2 'RUVprop,infusion' = 0.38 (RSE 12%)

    # Neuropathy model, Table 2 'PKPD model'. Logits of P(CTCAE >= k) per
    # subscore, structure from the supplement 2 $PRED code.
    b1_vipn <- fixed(-10); label("Logit intercept for grade >= 1 given previous grade 0") # Table 2 'B1|0' = -10 (FIX)
    b_prev_vipn <- 2.989; label("Logit shift when the previous grade of the subscore is above 0") # Table 2 'B1|0, if Grade > 0' = 2.989
    b2b1_nm_vipn <- -1.102; label("Logit increment from grade >= 1 to grade >= 2, neuralgia and motor subscores") # Table 2 'B2|0 (neuralgia/motor)' = -1.102
    b2b1_cs_vipn <- -3.583; label("Logit increment from grade >= 1 to grade >= 2, constipation and sensory subscores") # Table 2 'B2|0 (constipation/sensory)' = -3.583
    emax_nm_vipn <- 8.697; label("Maximum effect of the occasion AUC on the logit scale, neuralgia and motor") # Table 2 'Emax (neuralgia/motor)' = 8.697
    emax_cs_vipn <- 9.149; label("Maximum effect of the occasion AUC on the logit scale, constipation and sensory") # Table 2 'Emax (constipation/sensory)' = 9.149
    e_region_kenya_emax_vipn <- fixed(20); label("Kenyan Emax parameter; Kenyan Emax = exp(-100 x value) (unitless)") # Table 2 'Kenyan Emax parameter' = 20 (FIX)
    auc50_vipn <- 2.587; label("Occasion AUC giving half the maximum effect (ng*h/mL)") # Table 2 'AUC50 (ng*h/mL)' = 2.587

    # Censoring model, Table 2 (Equation 2; supplement 2 'LG_P3')
    b_cens_vipn <- -144.7; label("Baseline logit of censoring, before the first occasion") # Table 2 'PX0 B3 | Baseline' = -144.7
    b_cens_occ1_vipn <- 136.9; label("Logit shift of censoring from the first dose occasion") # Table 2 'Censoring | OCC > 0' = 136.9
    b_cens_occ2_vipn <- 5.764; label("Additional logit shift of censoring from the second dose occasion") # Table 2 'Censoring | OCC > 1' = 5.764
    b_cens_occ3_vipn <- 1.502; label("Additional logit shift of censoring from the third dose occasion") # Table 2 'Censoring | OCC > 2' = 1.502
    b_cens_prev_vipn <- 3.05; label("Logit shift of censoring after a previous grade >= 2 in any subscore") # Table 2 'Censoring | Grade > 0' = 3.05; supplement 2 applies it when any previous subscore is above 1
  })

  model({
    # Administration method: TINF > 0 is an infusion, TINF = 0 a bolus
    inf <- 0
    if (TINF > 0) inf <- 1

    cl <- exp(lcl + etalcl + e_region_kenya_cl_vp * REGION_KENYA) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp + e_region_kenya_cl_vp * REGION_KENYA) * (WT / 70)^e_wt_vc
    tlag <- exp(ltlag + etaltlag)
    vmax_bind <- exp(lvmax_bind)
    # Kb is reported in ug; the model amounts are in mg
    km_bind <- exp(lkm_bind) / 1000
    k31 <- exp(lk31)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    bind <- vmax_bind * central / (km_bind + central)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - bind + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- bind - k31 * peripheral2

    # Infusion records carry rate = -2; the bolus records ignore dur()
    alag(central) <- tlag * inf
    dur(central) <- TINF

    # mg / L = ug/mL; x 1000 for ng/mL
    Cc <- central / vc * 1000

    # Neuropathy: the drug effect is the AUC of the current occasion,
    # AUC = dose / CL, in ng*h/mL. podo() is NA before the first dose
    # (baseline assessment): no drug effect.
    dose_occ <- podo(central)
    if (is.na(dose_occ)) dose_occ <- 0
    auc_occ <- dose_occ / cl * 1000
    emax_nm <- emax_nm_vipn
    emax_cs <- emax_cs_vipn
    if (REGION_KENYA == 1) {
      emax_nm <- exp(-100 * e_region_kenya_emax_vipn)
      emax_cs <- exp(-100 * e_region_kenya_emax_vipn)
    }
    prev_vipn <- 0
    if (PREV_AE_SCORE > 0) prev_vipn <- 1
    lgt_ge1_nm <- b1_vipn + emax_nm * auc_occ / (auc50_vipn + auc_occ) + b_prev_vipn * prev_vipn
    lgt_ge1_cs <- b1_vipn + emax_cs * auc_occ / (auc50_vipn + auc_occ) + b_prev_vipn * prev_vipn
    # Cumulative VIPN probabilities, not yet scaled for censoring
    p_vipn_ge1_nm <- expit(lgt_ge1_nm)
    p_vipn_ge2_nm <- expit(lgt_ge1_nm + b2b1_nm_vipn)
    p_vipn_ge1_cs <- expit(lgt_ge1_cs)
    p_vipn_ge2_cs <- expit(lgt_ge1_cs + b2b1_cs_vipn)

    # Censoring (Equation 2): no previous grade >= 2 in any subscore, and
    # with one. Equation 3 then gives the four outcome probabilities of a
    # subscore as P(0) = (1 - Pc)(1 - P(>=1)), P(1) = (1 - Pc)(P(>=1) - P(>=2)),
    # P(>=2) = (1 - Pc) P(>=2) and P(censored) = Pc.
    occ1 <- 0
    if (CYCLE > 0) occ1 <- 1
    occ2 <- 0
    if (CYCLE > 1) occ2 <- 1
    occ3 <- 0
    if (CYCLE > 2) occ3 <- 1
    lgt_cens <- b_cens_vipn + b_cens_occ1_vipn * occ1 + b_cens_occ2_vipn * occ2 + b_cens_occ3_vipn * occ3
    p_cens_vipn <- expit(lgt_cens)
    p_cens_prev_vipn <- expit(lgt_cens + b_cens_prev_vipn)

    expSd_i <- (1 - inf) * expSd_bolus + inf * expSd_infusion * exp(etaRUV)
    Cc ~ lnorm(expSd_i)
  })
}
