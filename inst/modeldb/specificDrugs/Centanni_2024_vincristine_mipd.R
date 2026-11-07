Centanni_2024_vincristine_mipd <- function() {
  description <- paste(
    "Simulation variant of the Centanni 2024 vincristine PK and neuropathy",
    "model for children with cancer, as used for the paper's",
    "precision-dosing trial: clearance of the Table 2 two-compartment model",
    "replaced by the two-group CYP3A5 mixture (42.3 L/h at 70 kg, 1.87-fold",
    "higher in the 10% high expressers) with 32.3% IIV and 36.3%",
    "inter-occasion variability (one occasion per dose); all other",
    "parameters as in Centanni_2024_vincristine. Linear disposition with",
    "allometric weight scaling (exponents fixed at 0.75 and 1, 70 kg",
    "reference), an infusion-only lag time and a multiplier on the",
    "documented infusion duration; log-additive residual error with separate",
    "magnitudes for bolus and infusion samples and IIV on the infusion",
    "residual. The neuropathy model gives the probability of CTCAE grade",
    ">= 1 and >= 2 for each of four subscores (sensory, motor, constipation,",
    "neuralgia) at each dosing occasion, driven by an Emax function of the",
    "occasion AUC = dose / CL, with a higher-risk intercept from the second",
    "occasion and a Markov term for a previous subscore above grade 0."
  )
  reference <- paste(
    "Centanni M, van de Velde ME, Uittenboogaard A, Kaspers GJL,",
    "Karlsson MO, Friberg LE. Model-Informed Precision Dosing to Reduce",
    "Vincristine-Induced Peripheral Neuropathy in Pediatric Patients: A",
    "Pharmacokinetic and Pharmacodynamic Modeling and Simulation Analysis.",
    "Clin Pharmacokinet. 2024;63(2):197-209. doi:10.1007/s40262-023-01336-1"
  )
  vignette <- "Centanni_2024_vincristine"
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
        "Table 2 footnote a). Cohort median 34 kg (range 12-95; Table 1)."
      ),
      source_name = "weight"
    ),
    TINF = list(
      description = "Documented duration of the intravenous vincristine infusion; 0 for a bolus",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the administration method and sets the modelled infusion",
        "duration. TINF > 0 marks an infusion: the dose enters after the",
        "infusion lag time, over TINF times the estimated duration multiplier,",
        "and the sample takes the infusion residual error. TINF = 0 marks a",
        "bolus (no lag, bolus residual error). The documented time is read",
        "as the vincristine bag time, which the multiplier stretches for the",
        "flush: the protocol infusion was 60 min, 38 min of vincristine and",
        "22 min of flushing (96 min = 60 + 36 min in four patients; 15 min in",
        "two; Methods 'Patients and Study Design'), so TINF = 38/60 h for the",
        "standard infusion. That reading reproduces the timing of the",
        "observed infusion profile of Figure 1, where concentrations fall",
        "from about 1 h. Infusion dose records must carry rate = -2 so",
        "that rxode2 takes the duration from the model; supply TINF as a data",
        "column because rxode2 does not expose a dose record's duration to",
        "model()."
      ),
      source_name = "administration method (bolus / infusion) and documented infusion time"
    ),
    CYCLE = list(
      description = "Vincristine dose occasion (1 = first dose, 2 = second, ...)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Only CYCLE > 1 versus CYCLE = 1 matters: from the second dosing",
        "occasion the logit intercept for grade >= 1 VIPN is b1_occ_vipn",
        "instead of b1_vipn (Table 3 'B1|0, if OCC > 1'). Set it on each",
        "neuropathy row to the occasion of the dose that preceded it.",
        "Patients contributed one to five occasions."
      ),
      source_name = "OCC"
    ),
    OCC = list(
      description = "Dosing occasion for the inter-occasion variability on CL (1-5, one occasion per dose)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Selects the IOV eta on CL: etalcl_oc1 ... etalcl_oc5 for occasions",
        "1-5 (the paper's simulated course has five doses; Methods",
        "'Generation of the Target Populations'). Equal to CYCLE in this",
        "design; carried separately because CYCLE drives the occasion",
        "effect of the neuropathy model and OCC multiplexes the IOV.",
        "Records with an OCC outside 1-5 take no IOV."
      ),
      source_name = "OCC"
    ),
    CYP3A5_EXPR = list(
      description = "Membership of the high-clearance (CYP3A5 high-expresser) mixture group, 1 = high, 0 = low",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper did not genotype the patients; it fitted a two-group",
        "mixture on CL with the CL ratio fixed to 1.87 and the probability",
        "of the high-clearance group fixed to 10%, from earlier CYP3A5",
        "genotype work (Methods 'Approximation of Non-estimated Factors'),",
        "and described the 10% group as likely fast metabolizers (Results).",
        "The mixture group is carried here as this indicator; for a",
        "population simulation draw it per subject as Bernoulli(0.10), as",
        "the paper did ('extensive/poor metabolizer: 10/90%')."
      ),
      source_name = "mixture subpopulation"
    ),
    PREV_AE_SCORE = list(
      description = "CTCAE grade (0, 1, >= 2 coded 2) of the same VIPN subscore at the previous assessment",
      units = "(ordinal score)",
      type = "count",
      reference_category = 0,
      notes = paste(
        "Markov state on a three-level scale: CTCAE 0, 1 and >= 2 (the three",
        "grade 3-4 scores in the data were pooled with grade 2; Methods",
        "'Pharmacodynamic Model'). Only PREV_AE_SCORE > 0 versus 0 enters",
        "the model (Table 3 'B > 0|Grade > 0'). Each of the four subscores",
        "(peripheral sensory neuropathy, peripheral motor neuropathy,",
        "constipation, neuralgia) is its own chain, so the value on a row is",
        "the previous grade of that row's subscore; 0 at the baseline",
        "assessment."
      ),
      source_name = "previous VIPN subscore (PDV)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on the PK parameters and on the neuropathy model; not retained (Results 'Pharmacokinetic Model' and 'Pharmacodynamic Model'). Cohort median 9 years (range 2-17; Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested on the PK parameters and on the neuropathy model; not retained. Used only to compute the 1.5 or 2 mg/m^2 dose (2 mg cap). Cohort median 1.2 m^2 (range 0.5-2.2; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on the neuropathy model; not retained. 16 of 35 patients (46%) male (Table 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "vincristine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vincristine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    n_observations = "425 vincristine concentrations; 276 neuropathy assessments",
    age_range = "2-17 years",
    age_median = "9 years",
    weight_range = "12-95 kg",
    weight_median = "34 kg",
    bsa_range = "0.5-2.2 m^2 (median 1.2)",
    sex_female_pct = 54,
    disease_state = "Pediatric cancer treated with vincristine; 64% acute lymphoblastic leukemia",
    dose_range = "Vincristine 1.5 or 2 mg/m^2 (2 mg cap) per protocol, one to five dosing occasions per patient, as an IV bolus (20 patients) or an IV infusion (15 patients; 60 min standard, 96 min in four and 15 min in two patients)",
    regions = "Not stated; the patients are the PK substudy of the randomized bolus-versus-infusion trial reported by van de Velde et al. 2020 (reference 9)",
    sampling = "10, 20, 30, 40, 60, 75 and 140 min after the bolus or the start of the infusion, and at 24 h when still hospitalized; LC-MS/MS, LLOQ 0.1 ng/mL",
    notes = paste(
      "Demographics from Table 1 (n = 35). The text gives 20 bolus and 15",
      "infusion patients, while Table 1 labels the 20 (57%) row 'Patients",
      "with infusion'. Each dosing occasion was handled as a separate",
      "individual when the Table 2 model was estimated; this variant takes",
      "its CL terms from the second-stage estimation that split the",
      "per-occasion empirical Bayes CL estimates into IIV and IOV and fitted",
      "the CYP3A5 mixture (Methods 'Approximation of Non-estimated Factors',",
      "Eq. 1; Results 'Pharmacokinetic Model'), which",
      "the paper used to generate true exposures in its precision-dosing",
      "simulations; the estimation model is Centanni_2024_vincristine."
    )
  )

  ini({
    # PK, Table 2 (allometric scaling to 70 kg, Table 2 footnote a)
    lcl <- log(42.3); label("Clearance at 70 kg, low-clearance mixture group (L/h)") # Results 'Pharmacokinetic Model': 42.3 L/h/70 kg for the 90% group
    e_cyp3a5_expr_cl <- fixed(1.87); label("CL ratio of the high-clearance (CYP3A5 high-expresser) group (unitless)") # Methods 'Approximation of Non-estimated Factors': ratio fixed to 1.87 from reference 13; 42.3 x 1.87 = 79.1 L/h/70 kg as reported
    lvc <- log(23.3); label("Central volume at 70 kg (L)") # Table 2 Vc = 23.3 L (RSE 16%)
    lq <- log(114); label("Intercompartmental clearance at 70 kg (L/h)") # Table 2 Q = 114 L/h (RSE 14%)
    lvp <- log(576); label("Peripheral volume at 70 kg (L)") # Table 2 Vp = 576 L (RSE 20%)
    ltlag <- log(0.058); label("Lag time of an infusion, bag to bloodstream (h)") # Table 2 LAGT = 0.058 h (RSE 18%)
    tinf_mult <- 1.30; label("Multiplier on the documented infusion duration (unitless)") # Table 2 INFT = 1.30 (RSE 4%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of weight on CL and Q (unitless)") # Methods 'Pharmacokinetic Model': exponent 0.75 for CL
    e_wt_vc <- fixed(1); label("Allometric exponent of weight on Vc and Vp (unitless)") # Methods 'Pharmacokinetic Model': exponent 1 for distribution parameters

    # IIV and IOV on CL from the second-stage estimation (Eq. 1; Results
    # 'Pharmacokinetic Model'). CV% = sqrt(exp(omega^2) - 1) as in Table 2
    # footnote b, so omega^2 = log(1 + CV^2). One IOV variance shared by the
    # five occasions.
    etalcl ~ 0.0992 # IIV CL CV 32.3% -> log(1 + 0.323^2)
    etalcl_oc1 ~ 0.1238 # IOV CL CV 36.3% -> log(1 + 0.363^2)
    etalcl_oc2 ~ fixed(0.1238) # equal to the occasion-1 IOV variance
    etalcl_oc3 ~ fixed(0.1238) # equal to the occasion-1 IOV variance
    etalcl_oc4 ~ fixed(0.1238) # equal to the occasion-1 IOV variance
    etalcl_oc5 ~ fixed(0.1238) # equal to the occasion-1 IOV variance

    # Other IIV, Table 2; omega^2 = log(1 + CV^2)
    etalvc ~ 0.4782 # Table 2 'IIV Vc' = 78.3% -> log(1 + 0.783^2)
    etaltlag ~ 0.4016 # Table 2 'IIV LAGT' = 70.3% -> log(1 + 0.703^2)
    etaRUV ~ 0.0765 # Table 2 'IIV RUV' = 28.2% (infusion samples only) -> log(1 + 0.282^2)

    # Residual error: additive on log-transformed concentrations (Table 2 footnote c)
    expSd_bolus <- 0.137; label("Log-scale residual SD, bolus samples") # Table 2 'RUVadd,bolus' = 0.137 (RSE 18%)
    expSd_infusion <- 0.204; label("Log-scale residual SD, infusion samples") # Table 2 'RUVadd,infusion' = 0.204 (RSE 27%)

    # Neuropathy model, Table 3. Logit of P(CTCAE >= k) per subscore (Eqs. 2, 4, 5).
    # The occasion-1 intercept and Emax nearly cancel (-121 + 119), so the net
    # intercept and the occasion effect need more digits than Table 3 prints;
    # b1_occ_vipn and auc50_vipn come from the expit curves of Online Resource 2
    # (Figure S6), which the authors computed from the parameter estimates.
    b1_vipn <- -121; label("Logit intercept for grade >= 1 given previous grade 0, first occasion") # Table 3 'B1|0' = -121 (RSE 45%)
    b1_occ_vipn <- -119.54; label("Logit intercept for grade >= 1 given previous grade 0, occasion > 1") # Table 3 'B1|0, if OCC > 1' = -120 (RSE 36%); -119.54 fitted to the Figure S6 curves (rounds to -120)
    b2b1_nm_vipn <- -1.22; label("Logit increment from grade >= 1 to grade >= 2, neuralgia and motor subscores") # Table 3 'B2|0 (neuropathy/motor)' = -1.22 (RSE 32%)
    b2b1_cs_vipn <- -3.59; label("Logit increment from grade >= 1 to grade >= 2, constipation and sensory subscores") # Table 3 'B2|0 (constipation/sensory)' = -3.59 (RSE 16%)
    b_prev_vipn <- 2.81; label("Logit shift when the previous grade of the subscore is above 0") # Table 3 'B > 0|Grade > 0' = 2.81 (RSE 17%)
    emax_vipn <- 119; label("Maximum effect of the occasion AUC on the logit scale") # Table 3 Emax = 119 (RSE 46%)
    auc50_vipn <- 0.164; label("Occasion AUC giving half the maximum effect (ng*h/mL)") # Table 3 prints 1.62 ng*h/mL (RSE 63%); 0.164 fitted to the Figure S6 curves (decimal slip)
  })

  model({
    # Administration method: TINF > 0 is an infusion, TINF = 0 a bolus
    inf <- 0
    if (TINF > 0) inf <- 1

    iov_cl <- (OCC == 1) * etalcl_oc1 + (OCC == 2) * etalcl_oc2 +
      (OCC == 3) * etalcl_oc3 + (OCC == 4) * etalcl_oc4 + (OCC == 5) * etalcl_oc5
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * e_cyp3a5_expr_cl^CYP3A5_EXPR
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Infusion records carry rate = -2; the bolus records ignore dur()
    alag(central) <- tlag * inf
    dur(central) <- TINF * tinf_mult

    # mg / L = ug/mL; x 1000 for ng/mL
    Cc <- central / vc * 1000

    # Neuropathy: the drug effect is the AUC of the current occasion,
    # AUC = dose / CL (Results 'Pharmacodynamic Model'), in ng*h/mL; CL
    # includes the occasion's IOV, as the true AUC of the simulations did
    # podo() is NA before the first dose (baseline assessment): no drug effect
    dose_occ <- podo(central)
    if (is.na(dose_occ)) dose_occ <- 0
    auc_occ <- dose_occ / cl * 1000
    drug_vipn <- emax_vipn * auc_occ / (auc50_vipn + auc_occ)
    b1 <- b1_vipn
    if (CYCLE > 1) b1 <- b1_occ_vipn
    prev_vipn <- 0
    if (PREV_AE_SCORE > 0) prev_vipn <- 1
    lgt_ge1 <- b1 + drug_vipn + b_prev_vipn * prev_vipn
    p_vipn_ge1 <- expit(lgt_ge1)
    p_vipn_ge2_nm <- expit(lgt_ge1 + b2b1_nm_vipn)
    p_vipn_ge2_cs <- expit(lgt_ge1 + b2b1_cs_vipn)

    expSd_i <- (1 - inf) * expSd_bolus + inf * expSd_infusion * exp(etaRUV)
    Cc ~ lnorm(expSd_i)
  })
}
