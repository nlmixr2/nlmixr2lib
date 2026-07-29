Bonate_2004_apomine <- function() {
  description <- paste(
    "Two-compartment population PK model for oral apomine (a synthetic",
    "bisphosphonate-ester anti-cancer agent) in 38 subjects -- 19 healthy",
    "adult males and 19 male and female patients with advanced solid",
    "tumours -- pooled from four model-development studies (Bonate 2004",
    "Studies 1, 2, 5, 6) with three validation studies (Studies 3, 4, 7).",
    "Apomine is administered orally with or without food in single and",
    "multiple-dose regimens over 30 to 2100 mg total daily dose.",
    "Disposition is a two-compartment model with linear elimination and a",
    "first-order absorption mixture: a dominant Group 1 subpopulation",
    "(97 %, paper P1 logit) with an estimable lag time and faster",
    "absorption rate, and a minority Group 2 subpopulation (3 %) with no",
    "lag time and slower absorption (Table 3). Apparent oral clearance is",
    "time-dependent via an empirical sigmoid Emax auto-induction model",
    "in elapsed time (CL = CL0 + CLmax * time^n / (t50^n + time^n);",
    "Table 3), reaching 50 % of the maximal induction-driven increment",
    "in about two days. Relative bioavailability F1 is dose-saturable",
    "(F1 = D50 / (Dose + D50)) with an additional fractional food effect",
    "(1 + theta_food * FED), where the F1max anchor is structurally fixed",
    "at 1 (Table 3). Cancer patients have lower baseline CL/F and lower",
    "central volume than healthy males (encoded via the DIS_CANCER",
    "indicator with log-additive effects e_cancer_cl and e_cancer_vc",
    "back-derived from Table 3); intercompartmental clearance and",
    "peripheral volume are common to both populations. Central volume",
    "scales proportionally with body weight at a fixed allometric",
    "exponent of 1.0 (Table 3). A bimodal high-Vp subpopulation",
    "observed in the four healthy-male multiple-dose subjects",
    "(Bonate 2004 Study 2) is encoded as a binary indicator MIX_HIGH_VP",
    "multiplying the typical peripheral volume by 23.5. Inter-occasion",
    "variability (18 % on CL and Vc, Table 3) is omitted from this file",
    "because occasion definitions are study-design-specific (per the",
    "standard nlmixr2lib practice; see vignette Assumptions and",
    "deviations).",
    sep = " "
  )
  reference <- paste(
    "Bonate PL, Floret S, Bentzen C (2004).",
    "Population pharmacokinetics of APOMINE: A meta-analysis in cancer",
    "patients and healthy males.",
    "Br J Clin Pharmacol 58(2):142-155.",
    "doi:10.1111/j.1365-2125.2004.02111.x",
    sep = " "
  )
  vignette <- "Bonate_2004_apomine"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Bonate 2004. Cohort weight range",
        "55-115 kg across the development studies (Table 2; healthy males",
        "55-94 kg, cancer patients 55-115 kg). Reference 75 kg with fixed",
        "allometric exponent 1.0 on Vc per Table 3 ('Power term for weight",
        "on central compartment volume = 1.00 Fixed')."
      ),
      source_name        = "WT"
    ),
    DIS_CANCER = list(
      description        = paste(
        "Advanced solid-tumour patient indicator: 1 = cancer patient,",
        "0 = healthy adult male (reference)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy adult male)",
      notes              = paste(
        "Encodes the patient-type stratification used in Bonate 2004 for",
        "the two distinct typical CL/F and Vc values reported in Table 3.",
        "Healthy males (reference): CL0 = 40.7 mL/h, Vc = 12.3 L; cancer",
        "patients: CL0 = 10.2 mL/h, Vc = 7.11 L. Encoded as log-additive",
        "shifts e_cancer_cl = log(10.2 / 40.7) = -1.384 and",
        "e_cancer_vc = log(7.11 / 12.3) = -0.548. The auto-induction",
        "parameters (CLmax, t50, n) and distribution parameters (Q, Vp)",
        "are common to both populations per Table 3 (asterisked footnote",
        "'common variability for both cancer patients and healthy males')."
      ),
      source_name        = "POP (paper used a binary patient-type indicator)"
    ),
    FED = list(
      description        = paste(
        "Fed-vs-fasted dose-record indicator: 1 = oral dose administered",
        "with food, 0 = fasted."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Time-varying at the dose-record level. In the Bonate 2004 cohort",
        "Studies 1, 2, and 5 were fasted (FED = 0) and Studies 3, 4, 6, 7",
        "were administered with food (FED = 1). Enters relative",
        "bioavailability as a fractional shift: F1 = D50 /",
        "(Dose + D50) * (1 + theta_food * FED) with theta_food = -0.360",
        "(Table 3), i.e. food reduces F1 by 36 %."
      ),
      source_name        = "Food (Table 1 Food column)"
    ),
    MIX_LAGGED_ABS = list(
      description        = paste(
        "Latent absorption-mixture class indicator: 1 = subject in the",
        "lagged-absorption Group 1 subpopulation (Bonate 2004 dominant",
        "class, ka = 1.77 /h with a 0.821 h lag time); 0 = subject in the",
        "no-lag Group 2 subpopulation (slower ka = 0.361 /h with lag",
        "time fixed at 0)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Group 2, no-lag rare-class minority)",
      notes              = paste(
        "Per-subject latent class assignment from the published mixture",
        "absorption model (Bonate 2004 equation (1.4) and Results). The",
        "population probability of Group 1 is P(1) = 1 / (1 + exp(P1))",
        "with P1 = -3.47 (Table 3), giving P(Group 1) = 0.970 and",
        "P(Group 2) = 0.030. The reference category (0 = Group 2) is",
        "chosen so the binary numerically maps onto the paper's class",
        "indicator; the dominant Group 1 is the non-reference category",
        "by intent (same encoding as Tsuji 2017 MIX_PDI). For typical-",
        "value simulation set MIX_LAGGED_ABS = 1 to reproduce the",
        "dominant 97 % phenotype; set MIX_LAGGED_ABS = 0 for the rarer",
        "no-lag phenotype. For population simulation, draw",
        "MIX_LAGGED_ABS ~ Bernoulli(0.970) per subject."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment)"
    ),
    MIX_HIGH_VP = list(
      description        = paste(
        "Latent peripheral-volume-class indicator: 1 = subject in the",
        "high-peripheral-volume subpopulation observed only in Bonate",
        "2004 Study 2 (n = 4 healthy-male multiple-dose subjects with",
        "Vp 23.5-fold higher than the rest of the cohort); 0 = subject",
        "in the typical Vp subpopulation (the other 34 of 38 model-",
        "development subjects)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (typical low-Vp majority)",
      notes              = paste(
        "Per-subject latent class assignment from Bonate 2004 Results:",
        "examination of the empirical Bayesian peripheral-volume",
        "estimates showed four subjects in Study 2 (healthy-male",
        "multiple-dose) with Vp many-fold higher than all other",
        "subjects, modelled as a multiplicative shift V3i = TVV3 *",
        "theta_Study2 for that subgroup with theta_Study2 = 23.5",
        "(Table 3). Cohort prevalence is 4 / 38 = 10.5 %. The paper",
        "notes that simulations with and without this multiplier showed",
        "minimal differences in concentration-time profiles, so for",
        "typical-value forward simulations the default MIX_HIGH_VP = 0",
        "is recommended; set MIX_HIGH_VP = 1 only to reproduce the",
        "Study 2 healthy-male multiple-dose subgroup."
      ),
      source_name        = "Study indicator combined with the Bayesian Vp histogram"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 38L,
    n_studies       = 4L,
    age_range       = "18-78 years (development studies; Table 2)",
    weight_range    = "55-115 kg (development studies; Table 2)",
    sex_female_pct  = 13.2,
    race_ethnicity  = "Not reported by race in Bonate 2004 Table 2.",
    disease_state   = paste(
      "Healthy adult males (Studies 1 and 2; n = 19) and adult male and",
      "female patients with advanced solid tumours (Studies 5 and 6;",
      "n = 19). Cancer patients received concomitant supportive-care",
      "medications but were excluded if on chemotherapy, immunotherapy,",
      "hormonal cancer therapy, radiation, CYP3A inhibitors, or other",
      "experimental medications."
    ),
    dose_range      = paste(
      "Oral apomine 30-2100 mg administered once or twice daily, with or",
      "without food. Study 1: single dose (fasting); Study 2: multiple",
      "dose qd (fasting); Study 5: once weekly (fasting); Study 6:",
      "multiple dose bd (food)."
    ),
    regions         = paste(
      "United Kingdom (Studies 1, 2: Tayside, Dundee, Scotland), United",
      "States (Studies 3-6: Arizona Cancer Center Tucson AZ; Cancer",
      "Therapy and Research Center San Antonio TX), and United Kingdom",
      "(Study 7: Charing Cross Hospital, London)."
    ),
    n_observations  = 801L,
    notes           = paste(
      "Model-development dataset 38 subjects with 801 plasma apomine",
      "concentration measurements (Bonate 2004 Results, paragraph 1).",
      "Demographics by study in Bonate 2004 Table 2; the development",
      "studies are 1, 2, 5, 6 -- Studies 3, 4, 7 are validation sets",
      "and not included in the parameter estimates. Two influential",
      "subjects (1.18 and 5.11) identified by jack-knife principal-",
      "components analysis were removed before the final parameter",
      "estimates in Table 3 were obtained (Results, paragraph 1)."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Structural PK -- healthy adult male reference (Bonate 2004 Table 3)
    # CL and Vc values converted from the paper's mL/h and L units to
    # L/h and L so dosing in mg and Vc in L give Cc = central/vc in
    # mg/L = ug/mL directly (matching the assay range 0.01-50 ug/mL).
    # -------------------------------------------------------------------
    lcl <- log(0.0407)   ; label("Baseline apparent oral clearance CL0/F in healthy males before auto-induction (L/h)")          # Bonate 2004 Table 3: CL0 = 40.7 mL/h (SE 8.07; BSV 68%)
    lvc <- log(12.3)     ; label("Apparent central volume of distribution Vc/F in healthy males at 75 kg reference (L)")        # Bonate 2004 Table 3: TVV2 = 12.3 L (SE 1.86; BSV 30%)

    # Cancer-vs-healthy disease effects on CL0 and Vc (log-additive
    # shifts back-derived from Table 3 paired values).
    e_cancer_cl <- log(10.2 / 40.7) ; label("Log-shift in CL0/F for cancer patients vs healthy males (unitless)")               # Bonate 2004 Table 3: CL0 cancer = 10.2 mL/h vs healthy 40.7 mL/h
    e_cancer_vc <- log(7.11 / 12.3) ; label("Log-shift in Vc/F for cancer patients vs healthy males (unitless)")                # Bonate 2004 Table 3: TVV2 cancer = 7.11 L vs healthy 12.3 L

    # -------------------------------------------------------------------
    # Auto-induction model (sigmoid Emax in elapsed time, common between
    # populations per Table 3 footnote *; no covariate effects retained
    # on t50 per paper Results paragraph 2).
    # CL(t) = CL0 + CLmax * t^n / (t50^n + t^n)
    # -------------------------------------------------------------------
    lclmax <- log(0.320) ; label("Maximum induction-driven increment in CL/F (L/h)")                                            # Bonate 2004 Table 3: CLmax = 320 mL/h (SE 85.3; variance too small to estimate)
    lt50   <- log(46.4)  ; label("Time to 50% maximal auto-induction (h)")                                                      # Bonate 2004 Table 3: t50 = 46.4 h (SE 17.6; BSV 88%)
    lhill  <- log(6.40)  ; label("Hill (shape) exponent for the auto-induction sigmoid (unitless)")                             # Bonate 2004 Table 3: Shape parameter (n) = 6.40 (SE 1.76; variance too small to estimate)

    # -------------------------------------------------------------------
    # Distribution parameters (common to both populations per Table 3)
    # -------------------------------------------------------------------
    lq           <- log(0.198) ; label("Inter-compartmental clearance Q/F (L/h)")                                               # Bonate 2004 Table 3: Q = 198 mL/h (SE 40.9; variance too small to estimate)
    lvp          <- log(1.83)  ; label("Apparent peripheral volume of distribution Vp/F (L; typical-Vp subgroup)")              # Bonate 2004 Table 3: V3 = 1.83 L (SE 0.679; BSV 141%)
    e_wt_vc      <- fixed(1.0) ; label("Allometric exponent of WT on Vc/F (unitless, fixed at 1)")                              # Bonate 2004 Table 3: Power term = 1.00 Fixed
    e_high_vp_vp <- log(23.5)  ; label("Log multiplier on Vp/F for the high-Vp subgroup (Study 2 healthy-male MD)")             # Bonate 2004 Table 3: Peripheral compartment multiplier for study 2 = 23.5 (SE 7.42)

    # -------------------------------------------------------------------
    # Absorption mixture (Group 1 = lagged dominant 97%; Group 2 = no
    # lag rare 3%; BSV common between groups per Table 3 footnote *)
    # -------------------------------------------------------------------
    lka1  <- log(1.77)  ; label("Absorption rate constant ka for Group 1 (lagged-absorption, 1/h)")                             # Bonate 2004 Table 3: ka Group 1 = 1.77 /h (SE 0.431; BSV 145%)
    ltlag <- log(0.821) ; label("Absorption lag time for Group 1 (h; Group 2 lag is fixed at 0)")                               # Bonate 2004 Table 3: Lag time Group 1 = 0.821 h (SE 0.0386); Lag time Group 2 = 0 fixed
    lka2  <- log(0.361) ; label("Absorption rate constant ka for Group 2 (no-lag absorption, 1/h)")                             # Bonate 2004 Table 3: ka Group 2 = 0.361 /h (SE 0.0362; BSV 145%)
    # Mixing logit P1 = -3.47 -> P(Group 1) = 1/(1 + exp(P1)) = 0.970;
    # encoded structurally via the MIX_LAGGED_ABS binary covariate, not
    # as an ini() parameter (the typical-value library does not estimate
    # mixture proportions; see covariateData[[MIX_LAGGED_ABS]]$notes).

    # -------------------------------------------------------------------
    # Bioavailability (dose-saturable with food effect; F1max fixed at 1)
    # F1 = (1 - F1max * Dose / (Dose + D50)) * (1 + theta_food * FED)
    #    = D50 / (Dose + D50) * (1 + theta_food * FED)            (F1max = 1)
    # -------------------------------------------------------------------
    ld50     <- log(128) ; label("Dose for 50% reduction in F1 (mg)")                                                           # Bonate 2004 Table 3: D50 = 128 mg (SE 39.4)
    e_food_f <- -0.360   ; label("Food effect on F1 (linear shift on the dose-saturable F1 expression)")                        # Bonate 2004 Table 3: Food effect on F1 (theta_food) = -0.360 (SE 0.119)

    # -------------------------------------------------------------------
    # Inter-individual variability (omega^2 = log(CV^2 + 1) on the log
    # scale; CV% from Table 3 BSV column)
    # -------------------------------------------------------------------
    etalcl  ~ 0.4267                                                                                                            # Bonate 2004 Table 3: BSV CL = 68% -> omega^2 = log(0.68^2 + 1) = 0.4267 (common to healthy + cancer)
    etalvc  ~ 0.08619                                                                                                           # Bonate 2004 Table 3: BSV Vc = 30% -> omega^2 = log(0.30^2 + 1) = 0.08619 (common to healthy + cancer)
    etalvp  ~ 0.7945                                                                                                            # Bonate 2004 Table 3: BSV Vp = 141% -> omega^2 = log(1.41^2 + 1) = 0.7945
    etalt50 ~ 0.5594                                                                                                            # Bonate 2004 Table 3: BSV t50 = 88% -> omega^2 = log(0.88^2 + 1) = 0.5594
    etalka  ~ 0.8126                                                                                                            # Bonate 2004 Table 3: BSV ka = 145% -> omega^2 = log(1.45^2 + 1) = 0.8126 (common to Group 1 and Group 2)

    # -------------------------------------------------------------------
    # Residual error (paper form: Cobs = Cpred * exp(eps1) + eps2, an
    # additive plus exponential-proportional model; for small CV this is
    # well-approximated by the nlmixr2 add(addSd) + prop(propSd) idiom).
    # -------------------------------------------------------------------
    propSd <- 0.11  ; label("Proportional residual SD (fraction)")                                                              # Bonate 2004 Table 3: Proportional error = 11%
    addSd  <- 0.168 ; label("Additive residual SD (ug/mL)")                                                                     # Bonate 2004 Table 3: Additive error = 0.168
  })

  model({
    # 1. Individual baseline CL0 (healthy reference + cancer log-shift)
    #    and the time-dependent total CL via the auto-induction sigmoid.
    #    Elapsed simulation time `time` corresponds to the paper's "time
    #    relative to first exposure to apomine" (Methods, "Base model
    #    development"). For the paper's single-dose Study 1, the authors
    #    held CL constant; with the formula below, induction at 24 h is
    #    ~1.4% of CLmax for the typical t50 = 46.4 h and n = 6.40, so
    #    the time-dependent CL has only a minor effect on single-dose
    #    simulations.
    cl0   <- exp(lcl + e_cancer_cl * DIS_CANCER + etalcl)
    clmax <- exp(lclmax)
    t50   <- exp(lt50 + etalt50)
    hill  <- exp(lhill)
    cl    <- cl0 + clmax * time^hill / (t50^hill + time^hill)

    # 2. Distribution parameters. Vc carries the cancer effect plus an
    #    allometric WT/(75 kg) scaling at fixed exponent 1.0; Vp carries
    #    the high-Vp subgroup multiplier; Q is common.
    vc <- exp(lvc + e_cancer_vc * DIS_CANCER + etalvc) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + e_high_vp_vp * MIX_HIGH_VP + etalvp)
    q  <- exp(lq)

    # 3. Absorption mixture: indicator-gated ka and lag-time. Shared
    #    etalka scales whichever group's ka is active for the subject;
    #    Group 2 lag time is structurally 0 (Table 3 fixed).
    ka1  <- exp(lka1 + etalka)
    ka2  <- exp(lka2 + etalka)
    ka   <- ka1 * MIX_LAGGED_ABS + ka2 * (1 - MIX_LAGGED_ABS)
    tlag <- exp(ltlag) * MIX_LAGGED_ABS

    # 4. Micro-constants for the two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. ODE system (two-compartment with first-order absorption).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 6. Bioavailability: F = D50/(Dose + D50) * (1 + e_food_f * FED).
    #    podo(depot) returns the most recent dose amount entering the
    #    depot compartment (in the units$dosing = mg declared above).
    d50      <- exp(ld50)
    f(depot) <- d50 / (podo(depot) + d50) * (1 + e_food_f * FED)

    # 7. Absorption lag time on the depot compartment.
    alag(depot) <- tlag

    # 8. Observation. Dose in mg, Vc in L -> central / vc is in mg/L =
    #    ug/mL, matching the GC/NPD assay range 0.01-50 ug/mL (Methods,
    #    "Analysis of apomine").
    Cc <- central / vc

    # 9. Residual error -- additive + proportional. See ini() comment
    #    for the paper's exponential-proportional + additive form.
    Cc ~ prop(propSd) + add(addSd)
  })
}
