Du_2019_magnesiumSulfate <- function() {
  description <- "Two-compartment population PK model of magnesium sulfate (MgSO4-7H2O) in women with preeclampsia, fitted to the CHANGE FROM BASELINE in serum magnesium, with allometric body-weight scaling on all four disposition parameters, a serum-creatinine power effect on clearance, and interoccasion variability on clearance gated by antepartum status; first-order intramuscular absorption parameters are fixed from Salinger 2013 (Du 2019)."
  reference <- "Du L, Wenning L, Migoya E, Xu Y, Carvalho B, Brookfield K, Witjes H, de Greef R, Lumbiganon P, Sangkomkamhang U, Titapant V, Duley L, Long Q, Oladapo OT. Population pharmacokinetic modeling to evaluate standard magnesium sulfate treatments and alternative dosing regimens for women with preeclampsia. J Clin Pharmacol 2019;59(3):374-385. doi:10.1002/jcph.1328"
  vignette <- "Du_2019_magnesiumSulfate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "magnesium", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "magnesium", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "magnesium", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Maternal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power effect on all four disposition parameters, centered on 85 kg. Du 2019 Table 3 fixes the exponent at 0.75 for CL and Q and at 1 for Vc and Vp; the final-model equations printed in Results normalise by '85 kg' verbatim. The 85 kg reference is the cohort median (Table 2 reports mean 90.3 +/- 20.2 kg, range 57-157, but the Model Simulations section identifies 85 kg as the middle/median weight and the Discussion cross-checks 'Vss = 32.4 L/85 kg', which equals the tabulated Vc 15.4 + Vp 17.0).",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Maternal serum creatinine concentration (baseline value, carried forward/backward for missing records)",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL only, encoded as (CREAT/0.8)^e_creat_cl with e_creat_cl = -0.731 from Du 2019 Table 3 ('Serum creatinine exponent for CL, theta'). The 0.8 mg/dL reference is printed verbatim in the final-model equation as the denominator '0.8 mg/dL' and is the cohort median (Table 2 mean 0.82 +/- 0.29 mg/dL, range 0.4-2.1). The negative exponent gives the conventional renal direction: higher serum creatinine lowers magnesium clearance. Du 2019 Methods state that the last observation was carried forward to impute missing serum creatinine during treatment, and that a missing baseline value was filled by carrying the first post-dose value backward.",
      source_name = "Cr"
    ),
    PREG = list(
      description = "Antepartum (still pregnant) versus postpartum status of the sampling occasion; 1 = antepartum, 0 = postpartum",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (postpartum)",
      notes = "Time-varying WITHIN a subject: Du 2019 continued magnesium for 24 h after delivery, so a woman contributes antepartum records (PREG = 1) before delivery and postpartum records (PREG = 0) afterwards. 270/623 samples (43.3%) were drawn before birth and 353 (56.4%) after. IMPORTANT -- this indicator carries NO fixed-effect coefficient in Du 2019. It gates the interoccasion-variability eta on clearance: the printed final-model equation is exp(eta_CL,i + AP_ij * eta_IOV,i), so the IOV deviation applies only to antepartum records and the postpartum occasion takes the reference (no IOV deviation). See the ini() note on etaiov_cl_1 and the vignette's Assumptions and deviations section.",
      source_name = "AP"
    )
  )

  # Covariates screened by Du 2019 but NOT retained in the final model. Table S1
  # of the supplement reports the univariate runs: age on CL (run 8, OFV
  # unchanged at 3211.7), gestational age on CL (run 9, 3206.5 vs 3211.7 --
  # dOFV 5.2, short of the chi-square(1) = 6.63 threshold for p < .01), age on
  # Vc (run 10, unchanged) and gestational age on Vc (run 11, unchanged).
  # Results states: "Effects of age and gestational age on CL were not
  # statistically significant." No point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Maternal age at baseline",
      units = "years",
      type = "continuous",
      notes = "Mean 30.0 +/- 7.3 years, range 19-44 (Du 2019 Table 2). Screened on CL and Vc; not retained and no point estimate reported."
    ),
    EGA = list(
      description = "Gestational age of the fetus at baseline",
      units = "weeks",
      type = "continuous",
      notes = "Mean 34.73 +/- 4.31 weeks, range 21.0-40.3 (Du 2019 Table 2). Screened on CL and Vc; not retained and no point estimate reported."
    ),
    BMI = list(
      description = "Maternal body mass index at baseline",
      units = "kg/m^2",
      type = "continuous",
      notes = "Mean 34.8 +/- 6.5 kg/m^2, range 20.9-52.3 (Du 2019 Table 2). Listed in Methods among the available characteristics assessed; body size was ultimately carried by WT allometry alone."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_studies = 1L,
    age_range = "19-44 years (mean 30.0 +/- 7.3)",
    age_median = "~30 years",
    weight_range = "57-157 kg (mean 90.3 +/- 20.2)",
    weight_median = "85 kg (the reference weight used in the covariate model)",
    sex_female_pct = 100,
    race_ethnicity = "not reported",
    disease_state = "Preeclampsia. Gestational age at baseline 21.0-40.3 weeks (mean 34.73 +/- 4.31). Serum creatinine 0.4-2.1 mg/dL (mean 0.82 +/- 0.29). Baseline serum magnesium 14-25 mg/L (mean 18.3 +/- 2.2).",
    dose_range = "All women in the analysis dataset received the same intravenous regimen: a 4 g MgSO4-7H2O loading infusion over 20 minutes followed by a 2 g/h continuous intravenous maintenance infusion, continued for 24 hours after delivery. The intramuscular regimens in Du 2019 Table 1 are SIMULATION scenarios only -- no intramuscular data were fitted.",
    regions = "United States (Stanford University IRB); source PK data from Brookfield et al.",
    notes = "Secondary population-PK analysis of the 92 preeclamptic women from the 111-woman cohort of Brookfield KF et al., Am J Obstet Gynecol 2016 (the 19 non-preeclamptic women were excluded). 623 serum magnesium concentrations; 370 (59.4%) drawn during intravenous treatment and 253 (40.6%) after discontinuation; 270 (43.3%) antepartum and 353 (56.4%) postpartum, the latter including 8 samples taken at delivery. Serial sampling at baseline, 0.5, 1, 2, 4 h and every 6 h during administration, then 1, 3, 6, 9 and 12 h after discontinuation. Assay LLOQ 0.2 mmol/L on a Siemens Dimension RxL Max. Estimation in NONMEM 7.3 by FOCE with eta-epsilon interaction; model development summarised in supplemental Table S1 (final model = run 12, OFV 3188.4).\n\nUNITS AND THE BASELINE. Doses are administered as MgSO4-7H2O (heptahydrate, MW 246.47) while serum concentrations are reported as elemental Mg (MW 24.305). The PK parameters in this file use dose units of mg of elemental Mg and concentration units of mg/L of elemental Mg, matching the sibling models Salinger_2013_magnesiumSulfate.R, Easterling_2018_magnesium_sulfate.R and Deng_2024_magnesiumSulfate.R. Convert administered MgSO4-7H2O grams to mg Mg by multiplying by 24.305/246.47 = 0.0986 (4 g = 394.4 mg Mg; a 1 g/h infusion = 98.6 mg Mg/h; a 10 g intramuscular dose = 986.2 mg Mg).\n\nCRITICALLY, Du 2019 modelled the CHANGE FROM BASELINE in serum magnesium, not the absolute concentration, because endogenous magnesium is present before dosing (Discussion, first limitation). `Cc` in this file is therefore the change from baseline in mg/L, and it is the quantity the combined residual error was fitted to. The paper's own simulations, Figures 2 and 3, Table 4 and the therapeutic (1.5-2.5 mmol/L) and toxicity (3.5 mmol/L) thresholds are all expressed as TOTAL serum magnesium, obtained by adding a constant baseline of 0.74 mmol/L = 18 mg/L (the observed median baseline, Model Simulations section). That constant is supplied here as the fixed parameter `lrbase` and the total is exposed as the derived variable `CcTotal`. Divide mg/L by 24.305 to obtain mmol/L.\n\nINTRAMUSCULAR ABSORPTION IS NOT FROM THIS COHORT. Every woman in the dataset was dosed intravenously, so Ka and F could not be estimated. Du 2019 took Ka = 0.317 /h and F = 0.862 from its reference 14, which is Salinger DH et al., BJOG 2013;120:894-900 -- the same paper extracted in this package as Salinger_2013_magnesiumSulfate.R, whose lka and lfdepot carry these identical values. Both are wrapped in fixed() here. The paper flags this as its second limitation."
  )

  ini({
    # Structural parameters from Du 2019 Table 3 (Parameter Estimates of the
    # Final Population Pharmacokinetic Model). Reference subject: WT = 85 kg,
    # CREAT = 0.8 mg/dL.
    # Dose units used here: mg of elemental Mg (multiply g of MgSO4-7H2O by 98.6).
    # Concentration units: mg/L of elemental Mg, as CHANGE FROM BASELINE.
    lcl <- log(3.72); label("Clearance for the reference subject (L/h)") # Du 2019 Table 3, CL (L/h), 3.5% RSE
    lvc <- log(15.4); label("Central volume of distribution for the reference subject (L)") # Du 2019 Table 3, Vc (L), 11.6% RSE
    lq <- log(3.66); label("Intercompartmental clearance for the reference subject (L/h)") # Du 2019 Table 3, Q (L/h), 24.5% RSE
    lvp <- log(17.0); label("Peripheral volume of distribution for the reference subject (L)") # Du 2019 Table 3, Vp (L), 9.8% RSE

    # Intramuscular absorption parameters. Du 2019 could not estimate these --
    # every woman in the analysis dataset was dosed intravenously -- and took
    # them from its reference 14 (Salinger 2013 BJOG, extracted in this package
    # as Salinger_2013_magnesiumSulfate.R, which carries the identical values).
    # Fixed, not estimated; see Du 2019 Model Simulations and the second
    # limitation in the Discussion.
    lka <- fixed(log(0.317)); label("Intramuscular first-order absorption rate constant (1/h)") # Du 2019 Model Simulations, from reference 14 (Salinger 2013)
    lfdepot <- fixed(log(0.862)); label("Intramuscular bioavailability (fraction)") # Du 2019 Model Simulations, from reference 14 (Salinger 2013)

    # Endogenous magnesium baseline. NOT a parameter of the fitted model -- Du
    # 2019 fitted the change from baseline -- but the paper fixes this constant
    # for every simulation it reports, so it is supplied here as a fixed value
    # to build the total-magnesium output CcTotal. 0.74 mmol/L * 24.305 = 18 mg/L.
    lrbase <- fixed(log(18)); label("Assumed endogenous baseline serum magnesium (mg/L)") # Du 2019 Model Simulations, 'Baseline magnesium concentration was assumed to be 0.74 mmol/L (18 mg/L)'

    # Covariate effects from Du 2019 Table 3 and the final-model equations
    # printed in Results:
    #   CL_i = CL * (Cr_i/0.8 mg/dL)^theta * (WT_i/85 kg)^0.75
    #              * exp(eta_CL,i + AP_ij * eta_IOV,i)
    #   Vc_i = Vc * (WT_i/85 kg)          * exp(eta_Vc,i)
    #   Q_i  = Q  * (WT_i/85 kg)^0.75
    #   Vp_i = Vp * (WT_i/85 kg)
    # The two weight exponents are single values shared across two disposition
    # parameters each, exactly as Table 3 reports them ('WT exponent for CL and
    # Q', 'WT exponent for Vc and Vp'), so they take the shared-exponent
    # e_<cov>_<param1>_<param2> form. Both were fixed to allometric theory
    # rather than estimated (Covariate Analysis: 'Allometric scaling factors
    # for body weight were fixed to values of 0.75 and 1 for CL and volume of
    # distribution (Vd), respectively').
    e_creat_cl <- -0.731; label("Power exponent of (CREAT/0.8) on CL (unitless)") # Du 2019 Table 3, serum creatinine exponent for CL, 14.2% RSE
    e_wt_cl_q <- fixed(0.75); label("Allometric power exponent of (WT/85) on CL and Q (unitless)") # Du 2019 Table 3, WT exponent for CL and Q
    e_wt_vc_vp <- fixed(1.0); label("Allometric power exponent of (WT/85) on Vc and Vp (unitless)") # Du 2019 Table 3, WT exponent for Vc and Vp

    # Random effects from Du 2019 Table 3. The table reports each as a VARIANCE
    # on the log scale with the %CV in parentheses; the Methods define
    # CV(%) = sqrt(exp(omega^2) - 1) * 100. Each tabulated CV reproduces from
    # the tabulated variance under that formula, confirming the scale:
    #   0.0749 -> 27.9%, 0.241 -> 52.2%, 0.056 -> 24.0% (table prints 23.9%).
    # Supplement Table S1 run 2 tested a CL-Vc eta correlation (OFV 3310.9 vs
    # 3315.0, dOFV 4.1, short of chi-square(1) = 6.63 for p < .01) and the final
    # model lineage (run 12 <- 7 <- 6 <- 4 <- 1) does not carry it, so the two
    # IIV terms are coded independently. No IIV was estimated on Q or Vp.
    etalcl ~ 0.0749; label("IIV on clearance (variance, log scale)") # Du 2019 Table 3, IIV on CL, 21.2% RSE, 14.3% shrinkage
    etalvc ~ 0.241; label("IIV on central volume (variance, log scale)") # Du 2019 Table 3, IIV on Vc, 30.7% RSE, 33.3% shrinkage

    # Interoccasion variability on clearance, antepartum versus postpartum.
    # Du 2019 prints the final-model clearance equation as
    #   exp(eta_CL,i + AP_ij * eta_IOV,i)
    # -- a SINGLE IOV deviation per subject, switched on by the antepartum
    # indicator, with the postpartum occasion taking the reference (no
    # deviation). The generic Methods formula P_i = TVP * exp(eta_Pi + kappa_ij)
    # would instead imply a separate kappa for each of the two occasions. Where
    # the generic Methods text and the printed final-model equation disagree,
    # the printed equation governs, so occasion 1 (antepartum) carries this eta
    # and occasion 2 (postpartum) carries none. See the vignette's Assumptions
    # and deviations section.
    etaiov_cl_1 ~ 0.056; label("IOV on clearance for the antepartum occasion (variance, log scale)") # Du 2019 Table 3, IOV on CL, 42.2% RSE, 48.9% shrinkage

    # Residual error from Du 2019 Table 3. The Methods print the combined form
    # verbatim as y_ij = yhat_ij * (1 + eps_prop,ij) + eps_add,ij, and Results
    # confirm 'The residual error was best described by a combined residual
    # error structure' (supplement Table S1 run 12, OFV 3188.4 vs 3211.7).
    # Both act on the CHANGE FROM BASELINE, which is the dependent variable.
    propSd <- 0.12; label("Proportional residual error (fraction)") # Du 2019 Table 3, proportional, 13.5% RSE
    addSd <- 4.97; label("Additive residual error (mg/L)") # Du 2019 Table 3, additive (mg/L), 7.2% RSE, 11.3% shrinkage
  })
  model({
    # Individual PK parameters, replicating the final-model equations printed in
    # Du 2019 Results. Reference subject: WT = 85 kg, CREAT = 0.8 mg/dL.
    # PREG is the antepartum indicator (1 = antepartum, 0 = postpartum) and
    # enters ONLY as the switch on the IOV eta -- Du 2019 estimates no
    # antepartum fixed effect.
    cl <- exp(lcl + etalcl + PREG * etaiov_cl_1) *
      (CREAT / 0.8)^e_creat_cl *
      (WT / 85)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 85)^e_wt_vc_vp
    q <- exp(lq) * (WT / 85)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 85)^e_wt_vc_vp
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with linear elimination from central and
    # first-order absorption from an intramuscular depot. Du 2019 Model
    # Structure: "MgSO4 was set to be dosed into the central compartment for
    # intravenous administration", so IV loading infusions and IV maintenance
    # infusions go to central (cmt = "central", with rate or dur set per record)
    # and intramuscular doses go to depot (cmt = "depot"). No intramuscular data
    # were fitted; the depot exists to run the paper's Table 1 simulation
    # scenarios.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Intramuscular bioavailability applies to depot doses only; intravenous
    # doses into central are unaffected.
    f(depot) <- exp(lfdepot)

    # Observation. Cc is the CHANGE FROM BASELINE in serum magnesium (mg/L of
    # elemental Mg) -- the dependent variable Du 2019 actually fitted, and the
    # quantity the combined residual error below was estimated on.
    Cc <- central / vc

    # Total serum magnesium, reconstructed the way Du 2019 reconstructs it for
    # every figure and table it reports: "the predicted magnesium change from
    # baseline was added to the baseline concentration (taking a population
    # average for magnesium baseline levels)". Compare CcTotal, not Cc, against
    # the paper's therapeutic (1.5-2.5 mmol/L) and toxicity (3.5 mmol/L)
    # thresholds; divide by 24.305 to convert mg/L to mmol/L.
    CcTotal <- Cc + exp(lrbase)

    Cc ~ prop(propSd) + add(addSd)
  })
}
