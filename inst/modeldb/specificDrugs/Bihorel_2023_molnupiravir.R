Bihorel_2023_molnupiravir <- function() {
  description <- "Two-compartment population PK model with Savic transit-compartment absorption and linear elimination for plasma beta-D-N4-hydroxycytidine (NHC), the circulating active nucleoside of the orally administered prodrug molnupiravir (MK-4482, EIDD-2801), in 1207 healthy adults and adults with COVID-19 pooled across one phase I, two phase II and one phase II/III trial (Bihorel 2023). Doses and concentrations are in molar units (molnupiravir 329.31 Da, NHC 259.2 Da). Absorption is an Erlang transit chain (NN = 7.84 compartments, mean transit time MTT = 0.435 h in the fasted-capsule reference) feeding a depot that empties first-order at ka = 0.797 1/h; MTT is raised 422% by a high-fat meal and lowered 61.6% for the oral solution relative to the capsule, neither of which alters the extent of absorption (relative bioavailability F1 fixed at 1). Apparent elimination clearance CL/F = 70.6 L/h in an 80 kg participant rises less-than-proportionally with body weight (power 0.412); apparent central volume Vc/F = 63.9 L in a man of BMI 28 kg/m2 rises with BMI (power 0.997) and is 33% lower in women. Interindividual variability is carried on CL/F (43.4% CV) and Vc/F (62.9% CV); the mean transit time carries inter-occasion variability (39.8% CV) over two occasions rather than IIV. Proportional residual error is stratified by trial phase (25.5% CV for the densely sampled phase I trial, 49.7% CV for the sparsely sampled phase II/III trials). No covariate effect moved the AUC(0-12) geometric mean ratio outside the 0.7-2.0 clinical comparability bounds, so no dose adjustment is recommended for any subpopulation studied."
  reference <- "Bihorel S, Cao Y, Chawla A, Birger R, Maas BM, Gao W, Roepcke S, Sardella S, Humphrey R, Kondragunta S, Jayaraman B, Martinho M, Painter W, Painter G, Holman W, De Anda C, Brown ML, Johnson MG, Paschke A, Rizk ML, Stone JA. Population pharmacokinetics of molnupiravir in adults with COVID-19: Lack of clinically important exposure variation across individuals. CPT Pharmacometrics Syst Pharmacol. 2023 Dec;12(12):1859-1871. doi:10.1002/psp4.13031. PMCID: PMC10725262."
  vignette <- "Bihorel_2023_molnupiravir"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  #
  # Molnupiravir is a prodrug that is essentially undetectable in plasma
  # (Discussion paragraph 1: "generally not detectable in plasma due to rapid
  # hydrolysis to NHC during absorption and first-pass metabolism"). Bihorel
  # 2023 therefore fits apparent parameters directly to plasma NHC while
  # dosing molnupiravir. Methods "Bioanalytical methods" converts dose amounts
  # with the molnupiravir molecular mass (329.31 Da) and concentrations with
  # the NHC molecular mass (259.2 Da), i.e. the model assumes 1 mol dosed
  # molnupiravir presents 1 mol NHC to the systemic circulation, with the
  # actual conversion fraction absorbed into the apparent (/F) scaling. The
  # depot accordingly holds dosed prodrug on a molar basis and the systemic
  # states hold NHC.
  compartmentData <- list(
    depot       = list(analyte = "molnupiravir", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "beta-D-N4-hydroxycytidine (NHC)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "beta-D-N4-hydroxycytidine (NHC)", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Less-than-proportional power effect on apparent elimination clearance, referenced to 80 kg: cl <- cl_typ * (WT / 80)^0.412 (Bihorel 2023 Table 2, 'Apparent central clearance in 80kg participants' and 'Power of body weight effect'). The 80 kg reference is stated in the Table 2 row label; the cohort median body weight was 85 kg (Table 1). The paper deliberately applies no allometry to the apparent distribution clearance Q/F or the apparent peripheral volume Vp/F, and flags that omission in Limitations as the reason the model should not be extrapolated to children.",
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on apparent central volume, referenced to 28 kg/m^2: vc <- vc_typ * (BMI / 28)^0.997 (Bihorel 2023 Table 2, 'Apparent central volume in 28kg/m 2 BMI male participants' and 'Power of BMI effect'). The 28 kg/m^2 reference is stated in the Table 2 row label and matches the reference individual used for the Discussion typical-profile simulation; the cohort median BMI was 30.4 kg/m^2 (Table 1). The estimated power of 0.997 is indistinguishable from 1, i.e. Vc/F is effectively proportional to BMI, but it is carried as estimated because Table 2 reports it with an RSE of 13.1% and no FIX flag.",
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Linear proportional shift on apparent central volume: vc <- vc_typ * (1 + (-0.330) * SEXF), i.e. a 33% lower Vc/F in women than in men at the same BMI (Bihorel 2023 Table 2, 'Proportional shift in female participants' = -0.330, RSE 11.8%; Results 'Base model refinement using phase III data' states '33% decrease in V C /F in women compared with men'). Male is the reference because the Table 2 typical value row is labelled 'male participants'. The effect raises Cmax in women by about 15% (Discussion) but does not change AUC(0-12), which depends only on CL/F.",
      source_name        = "SEX"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted or standard meal)",
      notes              = "Linear proportional shift on the absorption mean transit time: mtt <- mtt_typ * (1 + 4.22 * FED_HIGHFAT), i.e. a 422% increase in MTT after a high-fat meal (Bihorel 2023 Table 2, 'Proportional shift due to high-fat meal' = 4.22, RSE 6.29%; Results states '422% increase in MTT following a high-fat meal compared with fasting or a standard meal'). Time-varying per dose record: the food-effect part of MK-4482-004 dosed the same participants fasted and after a high-fat breakfast on two separate occasions in a balanced crossover. Food status was NOT collected in the phase II and III trials (Table 1, 'Food status collection' = No for all 1107 participants with COVID-19); Bihorel 2023 used mixture modelling at analysis stage 2 to assign a food status and, after a small-scale sensitivity analysis, FIXED the fraction of unknown-food-status participants treated as having eaten a high-fat meal at 25%, then hard-coded that assignment into the analysis dataset. Downstream users simulating a COVID-19 cohort should therefore draw FED_HIGHFAT ~ Bernoulli(0.25); see the vignette Assumptions and deviations. The high-fat meal changes Tmax and Cmax but not AUC, so molnupiravir may be given without regard to food (Discussion).",
      source_name        = "FOOD"
    ),
    FORM_SOLUTION = list(
      description        = "Oral-solution formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule or suspension)",
      notes              = "Linear proportional shift on the absorption mean transit time: mtt <- mtt_typ * (1 + (-0.616) * FORM_SOLUTION), i.e. a 61.6% shorter MTT for the oral solution than for the capsule reference (Bihorel 2023 Table 2, 'Proportional shift in oral solution' = -0.616, RSE 5.49%; Results states '61.6% decrease in MTT for oral solution compared with capsule or suspension'). The oral solution was used only in 36 of the 100 phase I MK-4482-004 participants (Table 1, 'Formulation'); all 1171 remaining participants, including every participant in MOVe-IN, MOVe-OUT and MK-4482-006, received capsules, so this indicator is 0 throughout any phase II/III simulation. Formulation affected the rate but not the extent of absorption (Discussion: 'Food status and formulation were not found to influence the extent of molnupiravir bioavailability').",
      source_name        = "FORM"
    ),
    OCC = list(
      description        = "Occasion index for inter-occasion variability on the absorption mean transit time",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1 and 2. Bihorel 2023 Table 2 note: 'The different occasions for IOV in MTT were labeled as occasion 1 and occasion 2. Participants were assumed to have only one occasion for IOV in MTT, except those enrolled in the food effect and MAD part of MK-4482-004.' The occasion count of exactly two is therefore stated by the source and is not an extraction-side construct. Decomposed inside model() into binary indicators oc1 and oc2 that multiplex the two IOV etas on log-MTT; the second variance is fixed equal to the first, encoding the shared NONMEM $OMEGA BLOCK(1) SAME. Records for participants with a single dosing occasion (every participant with COVID-19, and the single-ascending-dose part of MK-4482-004) take OCC = 1. Registered idiom; closest precedent is the two-occasion Chen_2023_nemonoxacin.R.",
      source_name        = "OCC"
    ),
    STUDY_MOV_PHASE23 = list(
      description        = "Phase II / phase II-III trial residual-error stratum indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase I trial MK-4482-004 observation)",
      notes              = "Record-level study-design property, not a subject-level covariate: it selects only the proportional residual-error magnitude. 1 = the NHC observation comes from MK-4482-006, MOVe-IN or MOVe-OUT (sparse sampling: 1-5 samples per participant, drawn around the ninth or tenth dose); 0 = the observation comes from the densely sampled phase I trial MK-4482-004 (7-26 samples per participant). Bihorel 2023 Table 2 reports two residual-variability rows, 'Phase I trials' = 0.0652 (25.5% CV) and 'Phase II/III trials' = 0.247 (49.7% CV); Results 'Base structural model development' notes that 'distinct residual variability models being estimated for the phase I and II studies' was one of the stage-2 refinements. Follows the paper-specific STUDY_<drug>_<phase> convention; the closest analogue is STUDY_NMV_PHASE23 for the Chan 2023 nirmatrelvir COVID-19 antiviral pooled analysis.",
      source_name        = "STUDY"
    )
  )

  # Covariates that Bihorel 2023 screened but did NOT retain in the final
  # model. Documentation only -- these names are deliberately absent from
  # model(). Each was either rejected during the formal forward-selection /
  # backward-elimination search (alpha = 0.01 forward, 0.001 backward), or
  # retained through stage 2 and then eliminated at stage 3.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not statistically significant on any PK parameter. Discussion: 'Age was not a statistically significant covariate and predicted NHC exposures were consistent even in the oldest participants (Figure 3).' Cohort median (range) 46 (18-91) years."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI / MDRD-style, BSA-normalized)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Screened; not retained. Discussion: 'Renal impairment was not identified as a statistically significant covariate in the analysis, which included 48% and 7% of participants with mild or moderate renal impairment, respectively. The predicted modest 18% AUC(0-12) increase with moderate renal impairment was not considered to be clinically relevant.' Reported only as impairment strata (normal / mild 60-89 / moderate 30-59 mL/min/1.73 m^2, Table 1), never as a continuous coefficient, so no usable point estimate exists."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (modified Child-Pugh)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested after the formal covariate screen using modified Child-Pugh criteria (Table S2) approximated from total bilirubin and albumin, with encephalopathy, ascites and INR assumed normal because they were not collected. Not retained: 'the available data in participants with mild (n = 60) and moderate (n = 3) hepatic impairment did not identify a distinguishable trend in NHC exposures with hepatic dysfunction' (Discussion). No coefficient reported."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (modified Child-Pugh)",
      units       = "(binary)",
      type        = "binary",
      notes       = "See HEPIMP_MILD. Only 3 of 1207 participants (0.2%) were classified moderate; Limitations flags the small subgroup explicitly. No coefficient reported."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Reported only as a post hoc forest-plot AUC(0-12) geometric mean ratio of 0.825 (Figure 3), which Discussion attributes to higher body weight / BMI in this subpopulation rather than to race itself. A forest-plot GMR is a geometric mean over shrunken empirical Bayes estimates, not a typical-value contrast, so it cannot be inverted into a coefficient. n = 63 (5.2%)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Post hoc forest-plot AUC(0-12) GMR 1.01 (Figure 3), i.e. no meaningful difference. A separate sensitivity analysis (Figure S10) compared Japanese with weight-matched non-Japanese participants and found no clear difference after accounting for body weight. n = 43 (3.6%). Not invertible into a coefficient; see RACE_BLACK."
    ),
    DIS_HOSPITALIZED = list(
      description = "Hospitalized-with-COVID-19 indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "This is the one covariate that entered the model and was then removed. A shift in the absorption parameter for hospitalized patients was carried through analysis stage 2 (on the zero-order duration D1) but, after the stage-3 reversion to the transit-compartment absorption model, 'the shift in MTT in hospitalized patients was close to 0 and poorly estimated (RSE: 87.2%)' and it was the single relationship dropped by the stage-3 backward elimination (Results, 'Base model refinement using phase III data' and 'Final model'). No final estimate is reported, so it cannot be encoded. n = 196 (16.2%)."
    ),
    CONMED_REMDESIVIR = list(
      description = "Baseline remdesivir use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Recorded in Table 1 (48 of 196 hospitalized participants, 4.0% overall; no non-hospitalized participant received remdesivir) but not reported as a tested or retained covariate anywhere in the paper. No coefficient exists."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1207L,
    n_studies        = 4L,
    n_observations   = 4202L,
    age_range        = "18-91 years",
    age_median       = "46 years",
    weight_range     = "36.1-172 kg",
    weight_median    = "85 kg",
    bmi_range        = "14.3-68.6 kg/m^2",
    bmi_median       = "30.4 kg/m^2",
    sex_female_pct   = 48.3,
    race_ethnicity   = c(
      White = 66.7, Other = 19.0, NativeAmericanOrAlaskaNative = 5.5,
      BlackOrAfricanAmerican = 5.2, Asian = 3.6, HawaiianOrPacificIslander = 0.1
    ),
    ethnicity        = c(NonHispanicOrLatino = 59.2, HispanicOrLatino = 40.8),
    disease_state    = "Healthy adults (100; 8.3%) and adults with COVID-19 (1107; 91.7%), the latter comprising 911 non-hospitalized (75.5%) and 196 hospitalized (16.2%) participants",
    dose_range       = "MK-4482-004 (phase I, healthy): 50-1600 mg single oral dose fasted, a fasted-vs-high-fat-breakfast balanced crossover food-effect part, and 50-800 mg every 12 h for 6 days. MK-4482-006, MOVe-IN and MOVe-OUT part 1: 200, 400 or 800 mg every 12 h for 5 days. MOVe-OUT part 2 (phase II/III): 800 mg every 12 h for 5 days.",
    regions          = c(Europe = 41.0, AsiaPacific = 29.3, NorthAmerica = 20.5, Africa = 7.0, LatinAmerica = 2.2),
    renal_function   = "Normal 45.0%; mild impairment (eGFR 60-89 mL/min/1.73 m^2) 48.1%; moderate impairment (eGFR 30-59 mL/min/1.73 m^2) 7.0%",
    hepatic_function = "Normal 94.8%; mild impairment 5.0%; moderate impairment 0.2%. Hepatic function was not captured directly; the degree of impairment was approximated with a modified Child-Pugh score (Table S2) from total bilirubin and albumin, with encephalopathy, ascites and INR assumed normal because they were not collected.",
    formulation      = "Capsule 97.0% (1171); oral solution 3.0% (36, all in MK-4482-004)",
    co_medication    = "48 of 196 hospitalized participants (24.5%; 4.0% overall) had taken or were receiving remdesivir when molnupiravir was started; no non-hospitalized participant received remdesivir",
    notes            = "Demographics and baseline characteristics from Bihorel 2023 Table 1. Four randomized, double-blind, placebo-controlled trials: MK-4482-004 (phase I, NCT04392219, n = 100), MK-4482-006 (phase IIa non-hospitalized, NCT04405570, n = 66), MK-4482-001 / MOVe-IN (phase II hospitalized, NCT04575584, n = 196) and MK-4482-002 / MOVe-OUT (phase II/III non-hospitalized, NCT04575597, n = 845 across parts 1 and 2). 4847 plasma NHC records were collected and 4202 (88%) were analysed; exclusions were mostly below-LLOQ samples from the low-dose phase I arms (LLOQ 5 ng/mL in MK-4482-004 and MK-4482-006, 1 ng/mL in MOVe-IN and MOVe-OUT). Sampling was dense in healthy participants (7-26 samples each) and sparse in participants with COVID-19 (1-5 samples each), so 72.8% of the overall population contributed at most two samples. Estimation was FOCE with interaction in NONMEM 7 level 3."
  )

  ini({
    # ---- Absorption (Bihorel 2023 Table 2) -------------------------------
    # Stage-1/stage-3 structure of Figure 1a: an Erlang transit chain with NN
    # compartments and mean transit time MTT delivers drug into a depot, which
    # then empties into central at first-order rate ka.
    lka <- log(0.797)
    label("First-order absorption rate constant from depot to central (1/h)")
    # Table 2 row 'k a / First-order absorption rate constant, 1/h' = 0.797 (RSE 2.57%)

    lmtt <- log(0.435)
    label("Mean absorption transit time, fasted capsule reference (h)")
    # Table 2 row 'MTT / Mean absorption transit time, h' = 0.435 (RSE 5.39%)

    nn <- 7.84
    label("Number of Savic-style transit compartments (unitless)")
    # Table 2 row 'NN / Number of transit compartments' = 7.84 (RSE 16.5%).
    # Estimated, not fixed, and non-integer; rxode2's transit() uses the
    # continuous Savic 2007 analytical form (lgamma(nn + 1)), so the published
    # non-integer value is carried as-is rather than rounded.

    e_highfat_mtt <- 4.22
    label("Proportional shift in MTT after a high-fat meal (unitless)")
    # Table 2 row 'Proportional shift due to high-fat meal' = 4.22 (RSE 6.29%);
    # Results: "422% increase in MTT following a high-fat meal compared with
    # fasting or a standard meal"

    e_solution_mtt <- -0.616
    label("Proportional shift in MTT for oral solution vs capsule (unitless)")
    # Table 2 row 'Proportional shift in oral solution' = -0.616 (RSE 5.49%);
    # Results: "61.6% decrease in MTT for oral solution compared with capsule
    # or suspension"

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F1 (unitless)")
    # Table 2 row 'F1 / Relative bioavailability' = 1.00, %RSE column reads
    # "Fixed". Food status and formulation changed absorption rate only, not
    # extent (Discussion), so there is no covariate effect on F1.

    # ---- Disposition (Bihorel 2023 Table 2) ------------------------------
    lcl <- log(70.6)
    label("Apparent elimination clearance CL/F in an 80 kg participant (L/h)")
    # Table 2 row 'CL/F / Apparent central clearance in 80kg participants, L/h'
    # = 70.6 (RSE 1.97%)

    e_wt_cl <- 0.412
    label("Power of body weight on CL/F (unitless)")
    # Table 2 row 'Power of body weight effect' = 0.412 (RSE 14.0%)

    lvc <- log(63.9)
    label("Apparent central volume Vc/F in a man of BMI 28 kg/m2 (L)")
    # Table 2 row 'V C /F / Apparent central volume in 28kg/m 2 BMI male
    # participants, L' = 63.9 (RSE 5.07%)

    e_bmi_vc <- 0.997
    label("Power of BMI on Vc/F (unitless)")
    # Table 2 row 'Power of BMI effect' = 0.997 (RSE 13.1%)

    e_sexf_vc <- -0.330
    label("Proportional shift in Vc/F in women vs men (unitless)")
    # Table 2 row 'Proportional shift in female participants' = -0.330
    # (RSE 11.8%); Results: "33% decrease in V C /F in women compared with men"

    lq <- log(2.99)
    label("Apparent distribution clearance Q/F (L/h)")
    # Table 2 row 'Q/F / Apparent distribution clearance, L/h' = 2.99
    # (RSE 5.70%). No body-size covariate: Limitations states "the selected
    # model does not include any allometric scaling of apparent distribution
    # clearance and apparent peripheral volume".

    lvp <- log(68.3)
    label("Apparent peripheral volume Vp/F (L)")
    # Table 2 row 'V P /F / Apparent peripheral volume, L' = 68.3 (RSE 14.6%)

    # ---- Interindividual variability (Bihorel 2023 Table 2) --------------
    # Table 2 reports the "Magnitude of variability" column as %CV, so each
    # variance is recovered as omega^2 = log(1 + CV^2) for a log-normal eta.
    # Q/F, Vp/F, ka, MTT and NN are all marked "NE" (not estimated) in that
    # column and therefore carry no IIV.
    etalcl ~ 0.172573
    # Table 2 CL/F variability = 43.4 %CV (RSE 8.99%); log(1 + 0.434^2) = 0.172573.
    # Shrinkage 9.3% (Table 2 note).

    etalvc ~ 0.333371
    # Table 2 V C /F variability = 62.9 %CV (RSE 24.1%); log(1 + 0.629^2) = 0.333371.
    # Shrinkage 27.9% (Table 2 note).

    # ---- Inter-occasion variability on MTT (Bihorel 2023 Table 2) --------
    # Results, stage 3: "variability in MTT was modeled using IOV after testing
    # that IIV in MTT reduced to zero when associated with IOV". Exactly two
    # occasions per the Table 2 note. The second variance is fixed equal to the
    # first, encoding the NONMEM $OMEGA BLOCK(1) SAME that a single reported
    # IOV magnitude implies.
    etaiov_mtt_1 ~ 0.147049
    # Table 2 row 'IOV in MTT' = 39.8 %CV (RSE 17.6%); log(1 + 0.398^2) = 0.147049.
    # Shrinkage 27.6% (occasion 1) and 7.8% (occasion 2) per the Table 2 note.
    etaiov_mtt_2 ~ fixed(0.147049)

    # ---- Residual variability (Bihorel 2023 Table 2) ---------------------
    # Table 2 prints the variance in the "Final estimate" column and its square
    # root as %CV: sqrt(0.0652) = 0.2553 -> 25.5 %CV and sqrt(0.247) = 0.4970
    # -> 49.7 %CV. The model is a constant-CV (proportional) error model.
    propSdPhase1 <- sqrt(0.0652)
    label("Proportional residual error, phase I trial MK-4482-004 (fraction)")
    # Table 2 'Residual variability / Phase I trials' = 0.0652 (RSE 11.9%), 25.5 %CV

    propSdPhase23 <- sqrt(0.247)
    label("Proportional residual error, phase II and II/III trials (fraction)")
    # Table 2 'Residual variability / Phase II/III trials' = 0.247 (RSE 4.00%), 49.7 %CV
  })

  model({
    # --- 1. Individual disposition parameters ----------------------------
    # CL/F rises less-than-proportionally with body weight about the 80 kg
    # reference stated in the Table 2 row label.
    cl <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl

    # Vc/F rises with BMI about the 28 kg/m2 male reference and is shifted down
    # proportionally in women. Both effects enter as reported: a power function
    # for the continuous covariate and a linear proportional shift for sex.
    vc <- exp(lvc + etalvc) * (BMI / 28)^e_bmi_vc * (1 + e_sexf_vc * SEXF)

    # No covariate acts on Q/F or Vp/F (Limitations).
    q  <- exp(lq)
    vp <- exp(lvp)

    # --- 2. Individual absorption parameters ------------------------------
    # OCC is decomposed into mutually exclusive occasion indicators that
    # multiplex the two IOV etas on log-MTT. Single-occasion records pass
    # OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # MTT carries the food and formulation effects as linear proportional
    # shifts. FED_HIGHFAT and FORM_SOLUTION are mutually exclusive in the
    # source data: the oral solution was given only in the phase I trial, whose
    # high-fat arm used capsules.
    mtt <- exp(lmtt + iov_mtt) *
      (1 + e_highfat_mtt * FED_HIGHFAT) *
      (1 + e_solution_mtt * FORM_SOLUTION)

    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

    # --- 3. Micro-constants ----------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- 4. ODE system (Bihorel 2023 Figure 1a) --------------------------
    # transit() is rxode2's implementation of the Savic 2007 analytical Erlang
    # transit-chain input rate,
    #   ktr = (nn + 1) / mtt
    #   rate = fdepot * podo(depot) * ktr * (ktr * tad(depot))^nn *
    #            exp(-ktr * tad(depot)) / gamma(nn + 1),
    # which reproduces a chain of nn transit compartments without integrating
    # them. It feeds the depot, which empties first-order at ka into the
    # two-compartment linear disposition system.
    d/dt(depot)       <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # transit() already feeds the whole dose in through the analytical chain
    # (it reads the raw amount from podo(depot) and applies fdepot itself), so
    # the dose event must not additionally deposit a bolus into depot.
    f(depot) <- 0

    # --- 5. Observation and error ----------------------------------------
    # Doses were converted to molar units with the molnupiravir molecular mass
    # and concentrations with the NHC molecular mass (Methods, "Bioanalytical
    # methods"), so amounts in nmol divided by vc in L give nmol/L directly.
    Cc <- central / vc

    # Two proportional residual-error magnitudes, selected by trial phase.
    ruvProp <- propSdPhase1 * (1 - STUDY_MOV_PHASE23) + propSdPhase23 * STUDY_MOV_PHASE23
    Cc ~ prop(ruvProp)
  })
}
