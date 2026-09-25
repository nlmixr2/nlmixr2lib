# Population PK model for lidocaine and three metabolites (MEGX, GX and
# 2,6-xylidine) after repeated application of the lidocaine 5% medicated plaster
# (Versatis) in post-herpetic neuralgia patients, from Bursi et al. (2017).
#
# IDENTIFICATION NOTE (this extraction's main incidental finding):
# `inst/modeldb/ddmore/NA_NA_lidocaine.R` is the DDMORE Foundation Model
# Repository entry DDMODEL00000281, extracted from a bundle whose header states
# that no linked publication could be identified. That bundle IS this analysis.
# The two agree on every structural feature and on all 16 thetas / 3 omegas /
# 4 sigmas to ~3 significant figures; additional concordances are the exact
# observation count (1989), the dosing record (AMT 21600 / RATE 1800 = 1800
# per hour for 12 h = one plaster), the DDMORE licence being registered to BAST
# Inc. Ltd (Bursi 2017 co-author Joachim Grevel's affiliation), and the bundle's
# run date of 29/11/2016 against this paper's 12/01/2017 online publication.
# The two files are cross-linked with `replicate_of`. See the vignette's
# 'Relationship to the DDMORE sibling model' section for the full concordance
# table and for the four discrepancies (subject count, and the third significant
# digit of five estimates) that show run249 is a neighbouring run rather than
# the published final model.

Bursi_2017_lidocaine <- function() {
  description <- paste(
    "Population PK model for lidocaine and three metabolites",
    "(monoethylglycinexylidide [MEGX], glycinexylidide [GX] and 2,6-xylidine)",
    "after repeated 12-h application of up to three lidocaine 5% medicated",
    "plasters in post-herpetic neuralgia patients, followed for up to 14.5",
    "months. NONMEM ADVAN5 general-linear four-compartment parent-metabolite",
    "model, one compartment per chemical entity. Each plaster is assumed to",
    "deliver lidocaine into the central compartment at a constant 1800 ug/h",
    "for 12 h; lidocaine has no direct elimination pathway and leaves only by",
    "metabolism to MEGX (k_megx_form, fixed 0.03 /h) and to 2,6-xylidine",
    "(k_xyl_form, fixed 0.007 /h), so apparent lidocaine elimination is their",
    "sum, 0.037 /h. MEGX is metabolised onward to GX (k_gx_form = 1.93 /h);",
    "GX and 2,6-xylidine each carry an elimination rate constant. The number",
    "of plasters applied simultaneously (DLVL > 2) switches both the GX",
    "elimination rate constant and the lidocaine apparent central volume,",
    "giving apparent clearances of 48.8 L/h for <= 2 plasters and 67.0 L/h for",
    "3 plasters, i.e. exposure that rises less than proportionally with dose.",
    "Five further binary covariates (bilirubin, creatinine clearance, CYP1A2",
    "substrate co-medication, body mass index, ALT) act additively on the GX",
    "elimination rate constant and one (ALT) on the 2,6-xylidine rate",
    "constant, with lactate dehydrogenase switching the latter's baseline.",
    "Because the metabolites were never dosed alone, the fraction of lidocaine",
    "converted and the metabolite volumes are unidentifiable, so all three",
    "metabolite volumes are fixed at an arbitrary 100 L and the lidocaine",
    "volume is apparent (V/F); the paper back-calculates a topical",
    "bioavailability of roughly 5%."
  )
  reference <- "Bursi R, Piana C, Grevel J, Huntjens D, Boesl I. Evaluation of the Population Pharmacokinetic Properties of Lidocaine and its Metabolites After Long-Term Multiple Applications of a Lidocaine Plaster in Post-Herpetic Neuralgia Patients. Eur J Drug Metab Pharmacokinet. 2017;42(5):801-814. doi:10.1007/s13318-017-0400-7. PMCID: PMC5597703."
  vignette <- "Bursi_2017_lidocaine"
  replicate_of <- "inst/modeldb/ddmore/NA_NA_lidocaine.R"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Bursi 2017 Sect. 2.5 ('Blood Sampling') and Sect. 2.6
  # ('Bioanalysis') state that serum was separated and that lidocaine, MEGX, GX
  # and 2,6-xylidine were assayed in serum by LC-MS/MS, so the specimen is
  # serum rather than plasma. Amount units are ug because the plaster input
  # rate is given in ug/h (Sect. 2.7.7 assumption 1) and concentrations are
  # reported in ug/L (Table 3 footnote, Table 5 footnote).
  compartmentData <- list(
    central = list(analyte = "lidocaine", units = "ug", specimen = "serum", verified = TRUE),
    central_megx = list(analyte = "monoethylglycinexylidide (MEGX)", units = "ug", specimen = "serum", verified = TRUE),
    central_gx = list(analyte = "glycinexylidide (GX)", units = "ug", specimen = "serum", verified = TRUE),
    central_xyl = list(analyte = "2,6-xylidine", units = "ug", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    DLVL = list(
      description = "Number of lidocaine 5% medicated plasters applied simultaneously (Bursi 2017 Fig. 1 legend: 'DLVL dose level, i.e., number of plasters applied simultaneously'). Binarised in the model as DLVL > 2, which switches BOTH the typical-value GX elimination rate constant (Table 3: 1.44 -> 2.07 /h) and the lidocaine apparent central volume (Table 3: 1320 -> 1810 L).",
      units = "(count of plasters; 1, 2 or 3)",
      type = "count",
      reference_category = "DLVL <= 2 (one or two plasters; Table 3 rows 'k 30 for DLVL B2' and 'V 1 for DLVL B2').",
      notes = "Protocol allowed up to three plasters applied simultaneously for at most 12 h in each 24-h period (Sect. 1 and Sect. 2.4). Table 2 shows observations contributed at 1, 2 and 3 plasters. Results Sect. 3: 'The influence of DLVL can be interpreted as a decrease of systemic bioavailability when the dose, i.e., the number of plasters, increases.' DLVL also sets the input rate in the dataset (1800 ug/h per plaster), so it enters the simulation twice - once through the dose record and once as this covariate. The `NA_NA_lidocaine.R` sibling binarises the same column at the same threshold.",
      source_name = "DLVL"
    ),
    TBILI = list(
      description = "Total serum bilirubin. Binarised as TBILI > 0.53, which adds an additive modifier of -0.526 /h to the typical-value GX elimination rate constant (Table 3 row 'Effect of BIL [0.53 on k 30').",
      units = "mg/dL (see notes - the paper's own footnote says umol/L)",
      type = "continuous",
      reference_category = "TBILI <= 0.53 (additive modifier off).",
      notes = "UNIT DISCREPANCY, carried deliberately: the Table 3 footnote reads 'BIL bilirubin (lmol/L)' i.e. umol/L, but 0.53 umol/L is roughly an order of magnitude below any measurable total bilirubin (clinical reference range 5-21 umol/L), whereas 0.53 mg/dL sits squarely inside the normal range (0.3-1.2 mg/dL) and is a plausible cohort split. The threshold VALUE is transcribed verbatim from Table 3; the unit label is recorded here as mg/dL. Sect. 2.7.7 assumption 4 and the Discussion make clear these thresholds are distributional cohort splits ('laboratory safety parameters from the 30th percentile'), not clinical cutoffs - the sibling ALT threshold of 11 U/L is likewise below the clinical reference range. Flagged in the vignette Errata. Paper column name `BIL`; canonical `TBILI`.",
      source_name = "BIL"
    ),
    CRCL = list(
      description = "Creatinine clearance. Binarised as CRCL <= 52.7, which adds an additive modifier of -0.32 /h to the typical-value GX elimination rate constant (Table 3 row 'Effect of CL CR B52.7 on k 30').",
      units = "mL/min (not stated to be BSA-normalised)",
      type = "continuous",
      reference_category = "CRCL > 52.7 (additive modifier off).",
      notes = "Table 3 footnote gives 'CLCR creatinine clearance (mL/min)' with no BSA normalisation and no statement of the estimating equation (Cockcroft-Gault vs MDRD vs measured). Results Sect. 3 notes this covariate was at the limit of significance (dOFV = 6.72 against the 6.63 backward-deletion criterion) and was retained on mechanistic grounds: 'CLCR ... describes the relationship between excretion of GX and kidney function, which naturally decreases with increasing age.' The canonical CRCL register entry is BSA-normalised; this model's values are not, which matters if a user supplies mL/min/1.73m^2. Paper column name `CLCR`.",
      source_name = "CLCR"
    ),
    S1A2 = list(
      description = "Concomitant administration of a cytochrome P450 1A2 (CYP1A2) substrate. Binary. Adds an additive modifier of +0.852 /h to the typical-value GX elimination rate constant (Table 3 row 'Effect of CYP1A2 substrate on k 30').",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant CYP1A2 substrate (additive modifier off).",
      notes = "Sect. 2.7.4: 'The effect of concomitant medications, such as atenolol (ATEN), a cytochrome P450 1A2 (CYP1A2) substrate, and metoprolol and beta-blockers, as separate group, was explored ... Concomitant medication was expressed as binary data (concomitant medication yes or no).' Atenolol as a separate term was dropped during backward deletion; the pooled CYP1A2-substrate indicator was retained. Sect. 2.7.7 assumption 5: 'All co-medications, independent of start or stop dates, were assumed to be in use throughout the observation period', so the covariate is time-fixed per subject. The Discussion cautions the effect may be spurious: 'the confidence interval for the effect of CYP1A2, as well as for the effects of LDH and ALT, was found to be quite large.' This paper supplies the biological meaning that the `NA_NA_lidocaine.R` sibling could only guess at from the bundle's integer-coded `S1A2` column.",
      source_name = "CYP1A2 substrate"
    ),
    BMI = list(
      description = "Body mass index. Binarised as BMI > 27.9, which adds an additive modifier of +0.938 /h to the typical-value GX elimination rate constant (Table 3 row 'Effect of BMI [27.9 on k 30').",
      units = "kg/m^2",
      type = "continuous",
      reference_category = "BMI <= 27.9 (additive modifier off).",
      notes = "Cohort median BMI 26.3 kg/m^2 (Table 1), so the 27.9 threshold is near the middle of the distribution rather than at the WHO overweight (25) or obesity (30) cutoff. The Discussion states the direction explicitly - 'the effect of BMI on the elimination half-life of GX, k 30 , indicated a modest decrease for patients with a BMI [27.9 kg/m2' - a shorter half-life, i.e. a LARGER k30, which is the sign carried here - 'which could not be explained by any physiological mechanism.'",
      source_name = "BMI"
    ),
    ALT = list(
      description = "Serum alanine transaminase. Binarised as ALT > 11, which adds an additive modifier of -0.492 /h to the GX elimination rate constant AND +0.229 /h to the 2,6-xylidine elimination rate constant (Table 3 rows 'Effect of ALT [11 on k 30' and 'Effect of ALT [11 on k 40').",
      units = "U/L",
      type = "continuous",
      reference_category = "ALT <= 11 (both additive modifiers off).",
      notes = "The threshold of 11 U/L is below the lower end of the usual adult reference range (about 7-56 U/L), confirming these binarisations are distributional cohort splits rather than clinical hepatic-impairment cutoffs; Sect. 3 describes simulating 'subjects not affected by the pharmacokinetic covariates ... and by laboratory safety parameters from the 30th percentile'. The two ALT effects have OPPOSITE signs (GX down, 2,6-xylidine up), which is transcribed as printed. Paper column name `ALT`.",
      source_name = "ALT"
    ),
    LDH = list(
      description = "Serum lactate dehydrogenase. Binarised as LDH > 195, which switches the typical-value 2,6-xylidine elimination rate constant baseline (Table 3: 0.667 /h for LDH <= 195, 0.410 /h for LDH > 195).",
      units = "U/L",
      type = "continuous",
      reference_category = "LDH <= 195 (baseline rate constant 0.667 /h).",
      notes = "Encoded as a stratum switch rather than an additive modifier because Table 3 prints two parallel baseline estimates ('k 40 for LDH B195' and 'k 40 for LDH [195') rather than a reference value plus an increment - the same shape as the DLVL effects. 195 U/L is within the usual clinical reference range (about 140-280 U/L). The Discussion groups LDH with CYP1A2 and ALT as effects whose confidence intervals 'was found to be quite large' and which 'might be due to spurious effects emerging from a large number of tests.'",
      source_name = "LDH"
    )
  )

  # Covariates Bursi 2017 screened (Sect. 2.7.4) but did not retain in the final
  # model. Documented for provenance only; checkModelConventions() does not
  # require these to appear in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Sect. 2.7.4 tested weight 'on all model parameters'; not retained in the final model (Sect. 3 lists the retained set as DLVL, ATEN, BIL, ALT, CYP1A2 substrate, CLCR, BMI, LDH, of which only ATEN was subsequently deleted). Table 1: median 72.5 kg, range 38-114."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Sect. 2.7.4 tested age on all model parameters; not retained. Table 1: median 72 years, range 45-92."
    ),
    HEIGHT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Sect. 2.7.4 tested height on all model parameters; not retained. Table 1: median 165 cm, range 142-189."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Sect. 2.7.4 tested sex on all model parameters; not retained. Table 1: 119 of 212 patients female (56.1%)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units = "U/L",
      type = "continuous",
      notes = "Sect. 2.7.4 tested AST on elimination; not retained. No point estimate published."
    ),
    CK = list(
      description = "Creatine kinase",
      units = "U/L",
      type = "continuous",
      notes = "Sect. 2.7.4 tested CK on elimination; not retained. No point estimate published."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Sect. 2.7.4 tested gammaGT on elimination; not retained. No point estimate published."
    ),
    SCR = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Sect. 2.7.4 tested serum creatinine on elimination; not retained (the derived CLCR was). No point estimate published."
    ),
    SMOKER = list(
      description = "Current smoking status",
      units = "(binary)",
      type = "binary",
      notes = "Sect. 2.7.4 tested smoking status on all model parameters; not retained. Of interest because smoking induces CYP1A2, the isozyme whose substrate co-medication indicator WAS retained. No point estimate published."
    ),
    HR = list(
      description = "Heart rate",
      units = "beats/min",
      type = "continuous",
      notes = "Sect. 2.7.4 tested heart rate on absorption and elimination; not retained. No point estimate published."
    ),
    SBP = list(
      description = "Systolic blood pressure",
      units = "mmHg",
      type = "continuous",
      notes = "Sect. 2.7.4 tested systolic blood pressure on all model parameters; not retained. No point estimate published."
    ),
    DBP = list(
      description = "Diastolic blood pressure",
      units = "mmHg",
      type = "continuous",
      notes = "Sect. 2.7.4 tested diastolic blood pressure on all model parameters; not retained. No point estimate published."
    ),
    CONMED_ATENOLOL = list(
      description = "Concomitant atenolol (a CYP1A2 substrate) coadministration indicator",
      units = "(binary)",
      type = "binary",
      notes = "Sect. 3: atenolol (ATEN) was significant during forward inclusion on k30 but is the ONLY covariate removed during backward deletion - 'During the backward deletion, only ATEN was found of no significance in the final model.' The broader pooled CYP1A2-substrate indicator (S1A2) was retained in its place. No final-model point estimate published."
    ),
    URINEPH = list(
      description = "Urine pH",
      units = "(pH units)",
      type = "continuous",
      notes = "Sect. 2.7.4 tested urine pH on elimination; not retained. No point estimate published."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 212L,
    n_studies = 2L,
    age_range = "45-92 years",
    age_median = "72 years",
    weight_range = "38-114 kg",
    weight_median = "72.5 kg",
    sex_female_pct = 56.1,
    bmi_range = "13.5-44.6 kg/m^2",
    bmi_median = "26.3 kg/m^2",
    height_range = "142-189 cm",
    disease_state = "Post-herpetic neuralgia, defined as neuropathic pain persisting for at least 3 months after healing of a herpes zoster skin rash, with average pain intensity of at least 4 on the 11-point numeric rating scale at screening.",
    dose_range = "Up to three lidocaine 5% medicated plasters (700 mg lidocaine each) applied simultaneously but not overlapping, for a maximum of 12 h in each 24-h period. Each plaster is assumed to deliver lidocaine at a constant 1800 ug/h (Sect. 2.7.7 assumption 1), i.e. 21600 ug per plaster per 12-h application.",
    treatment_duration = "Up to 14.5 months overall: up to 10 weeks in the first (randomised-withdrawal) trial and up to 12 months in the second (open-label long-term) trial; the Discussion notes data extend to 18 months for some subjects.",
    notes = "Demographics from Table 1 ('Demographic characteristics of the patients in the pharmacokinetic population'), which reports male (N = 93), female (N = 119) and pooled (N = 212) columns; the values recorded here are the pooled 'All' column. Sect. 2.3: the PK population is every patient with at least one measurable lidocaine, MEGX, GX or 2,6-xylidine concentration plus adequate dosing and sampling history. Sparse sampling, up to five occasions per patient per trial (Sect. 2.5), contributing 1989 concentrations in total (Table 2: 513 lidocaine, 474 MEGX, 480 GX, 522 2,6-xylidine). Both trials enrolled patients aged 50 years and older, so the tabulated minimum ages of 53 (male) and 45 (female) reflect age at the PK sampling visit rather than an enrolment violation. Race / ethnicity and region are not reported. Note that the DDMORE sibling run (DDMODEL00000281) reports 325 individuals against the same 1989 observations; see the vignette Errata."
  )

  ini({
    # ---- Lidocaine metabolic-formation rate constants (Table 3) ----
    # Both fixed, so lidocaine's apparent elimination rate constant is their
    # sum. The Discussion confirms the sum arithmetically via the apparent
    # clearance: 'lidocaine K el was estimated at 0.0037 h -1 (sum of k 12 and
    # k 14 )' and then 'the apparent clearance was computed as 48.8 L/h after
    # the application of two or fewer plasters and 67.0 L/h after the
    # application of three plasters'. 48.8 / 1320 = 0.03697 and 67.0 / 1810 =
    # 0.03702, so the intended Kel is 0.037 /h = 0.03 + 0.007 and the printed
    # '0.0037' is a typographical slip for 0.037. This also settles the
    # conflicting k12 value in Sect. 4.1 Limitations, which says k12 and k14
    # 'were fixed to values similar to those estimated during covariate forward
    # inclusion (0.003 and 0.007 h -1 , respectively)': 0.003 + 0.007 = 0.010
    # would give an apparent clearance of 13.2 L/h, not 48.8. Table 3's 0.03 is
    # therefore the correct value and is what is encoded. See vignette Errata.
    lk_megx_form <- fixed(log(0.03))
    label("Lidocaine to MEGX formation rate constant (k12, 1/h)") # Table 3 row 'k 12 (h - 1 )': Fixed to 0.03, SEE 'Not applicable', 'Not estimated'
    lk_xyl_form <- fixed(log(0.007))
    label("Lidocaine to 2,6-xylidine formation rate constant (k14, 1/h)") # Table 3 row 'k 14 (h - 1 )': Fixed to 0.007, SEE 'Not applicable', 'Not estimated'

    # ---- MEGX to GX formation rate constant (Table 3) ----
    # Estimated; no covariate effects and no IIV in the final model.
    lk_gx_form <- log(1.93)
    label("MEGX to GX formation rate constant (k23, 1/h)") # Table 3 row 'k 23 (h - 1 )': 1.93, SEE 0.175, 95% CI (1.59; 2.27)

    # ---- GX elimination rate constant, DLVL strata (Table 3) ----
    # Table 3 prints two parallel baselines rather than a reference plus an
    # offset, which is the stratum-suffix case in parameter-names.md. Sect.
    # 2.7.4 eq. (3) gives the switch form directly: hTV = h1 * COV + h2 *
    # (1 - COV).
    lkel_gx_dlvlle2 <- log(1.44)
    label("GX elimination rate constant for <= 2 plasters (k30, 1/h)") # Table 3 row 'k 30 (h - 1 ) for DLVL B 2': 1.44, SEE 0.169, 95% CI (1.11; 1.77)
    lkel_gx_dlvlgt2 <- log(2.07)
    label("GX elimination rate constant for 3 plasters (k30, 1/h)") # Table 3 row 'k 30 (h - 1 ) for DLVL [ 2': 2.07, SEE 0.278, 95% CI (1.53; 2.61)

    # ---- Additive covariate modifiers on the GX elimination rate constant ----
    # These are ADDITIVE on the linear 1/h scale, not multiplicative: Table 3
    # gives each row the unit '(h - 1 )' (a multiplicative factor would be
    # dimensionless) and three of the five are NEGATIVE, which cannot be a
    # parameter value for a rate constant. They are the (h1 - h2) increments of
    # Sect. 2.7.4 eq. (3) with the DLVL baseline above supplying h2.
    e_tbili_kel_gx <- -0.526
    label("Additive change in k30 when bilirubin > 0.53 (1/h)") # Table 3 row 'Effect of BIL [0.53 on k 30 (h - 1 )': -0.526, SEE 0.148, 95% CI (-0.816; -0.236)
    e_crcl_kel_gx <- -0.32
    label("Additive change in k30 when creatinine clearance <= 52.7 mL/min (1/h)") # Table 3 row 'Effect of CL CR B 52.7 on k 30 (h - 1 )': -0.32, SEE 0.166, 95% CI (-0.645; 0.005)
    e_s1a2_kel_gx <- 0.852
    label("Additive change in k30 with concomitant CYP1A2 substrate (1/h)") # Table 3 row 'Effect of CYP1A2 substrate on k 30 (h - 1 )': 0.852, SEE 0.27, 95% CI (0.323; 1.381)
    e_bmi_kel_gx <- 0.938
    label("Additive change in k30 when body mass index > 27.9 kg/m^2 (1/h)") # Table 3 row 'Effect of BMI [27.9 on k 30 (h - 1 )': 0.938, SEE 0.309, 95% CI (0.332; 1.544)
    e_alt_kel_gx <- -0.492
    label("Additive change in k30 when ALT > 11 U/L (1/h)") # Table 3 row 'Effect of ALT [11 on k 30 (h - 1 )': -0.492, SEE 0.193, 95% CI (-0.87; -0.114)

    # ---- 2,6-xylidine elimination rate constant, LDH strata (Table 3) ----
    lkel_xyl_ldhle195 <- log(0.667)
    label("2,6-xylidine elimination rate constant for LDH <= 195 U/L (k40, 1/h)") # Table 3 row 'k 40 (h - 1 ) for LDH B 195': 0.667, SEE 0.0383, 95% CI (0.592; 0.742)
    lkel_xyl_ldhgt195 <- log(0.41)
    label("2,6-xylidine elimination rate constant for LDH > 195 U/L (k40, 1/h)") # Table 3 row 'k 40 (h - 1 ) for LDH [ 195': 0.41, SEE 0.0614, 95% CI (0.29; 0.53)

    e_alt_kel_xyl <- 0.229
    label("Additive change in k40 when ALT > 11 U/L (1/h)") # Table 3 row 'Effect of ALT [11 on k 40 (h - 1 )': 0.229, SEE 0.0975, 95% CI (0.038; 0.420)

    # ---- Lidocaine apparent central volume, DLVL strata (Table 3) ----
    # Apparent (V/F): the plaster input is the nominal delivered amount, not
    # the absorbed amount. Discussion: 1320 / 70 kg = 18.9 and 1810 / 70 kg =
    # 25.9, reported as 'about 19 and 26 L/kg (assuming a typical subject of
    # 70 kg)'.
    lvc_dlvlle2 <- log(1320)
    label("Lidocaine apparent central volume for <= 2 plasters (V1/F, L)") # Table 3 row 'V 1 (L) for DLVL B 2': 1320, SEE 99.5, 95% CI (1124; 1515)
    lvc_dlvlgt2 <- log(1810)
    label("Lidocaine apparent central volume for 3 plasters (V1/F, L)") # Table 3 row 'V 1 (L) for DLVL [ 2': 1810, SEE 184, 95% CI (1449; 2170)

    # ---- Metabolite apparent volumes (Table 3) ----
    # Sect. 3: 'Since lidocaine metabolites were not administered alone and the
    # true fraction of lidocaine converted to its metabolites is unknown, the
    # fraction of lidocaine converted to its metabolites and the apparent
    # volume of distribution of the metabolites are unidentifiable in the
    # model. Hence, the apparent volumes of distribution of the metabolites
    # were fixed to an arbitrary value (100 L).' Declared as three separate
    # fixed parameters (one per metabolite) rather than one shared name so a
    # later model can relax them independently; all three carry the single
    # Table 3 value.
    lvc_megx <- fixed(log(100))
    label("MEGX apparent volume of distribution (L)") # Table 3 row 'V 2 , V 3 , V 4 (L)': Fixed to 100, arbitrary (Sect. 3, unidentifiable)
    lvc_gx <- fixed(log(100))
    label("GX apparent volume of distribution (L)") # Table 3 row 'V 2 , V 3 , V 4 (L)': Fixed to 100, arbitrary (Sect. 3, unidentifiable)
    lvc_xyl <- fixed(log(100))
    label("2,6-xylidine apparent volume of distribution (L)") # Table 3 row 'V 2 , V 3 , V 4 (L)': Fixed to 100, arbitrary (Sect. 3, unidentifiable)

    # ---- Between-subject variability (Table 3) ----
    # Sect. 2.7.1 eq. (1): hi = hTV * exp(gi), i.e. exponential (log-normal)
    # IIV, which is what attaching the eta to the log-scale parameter gives.
    # Table 3 labels these 'Proportional on' and reports the VARIANCE in the
    # Estimate column with the CV% in the right-hand column; Sect. 2.7.1 says
    # 'The magnitude of IIV was expressed as coefficient of variation (%CV),
    # which was approximated by the square root of the variance estimate' -
    # confirmed arithmetically, sqrt(0.39) = 0.624, sqrt(0.2) = 0.447,
    # sqrt(0.312) = 0.559. The Estimate column is therefore the variance and is
    # used directly. Sect. 3: IIV was supported only on V1, k30 and k40.
    # The etas attach to the bare canonical stems because the paper estimates
    # ONE variance per parameter shared across both strata.
    etalkel_gx ~ 0.39 # Table 3 IIV 'Proportional on k 30': 0.39, SEE 0.127, 95% CI (0.141; 0.639), CV% 62.4
    etalkel_xyl ~ 0.2 # Table 3 IIV 'Proportional on k 40': 0.2, SEE 0.0424, 95% CI (0.117; 0.283), CV% 44.7
    etalvc ~ 0.312 # Table 3 IIV 'Proportional on V 1': 0.312, SEE 0.0757, 95% CI (0.164; 0.46), CV% 55.9

    # ---- Residual error (Table 3) ----
    # Sect. 2.7.1 eq. (2): Cij = Chat_ij + e_aij, a pure additive error on the
    # concentration scale, one epsilon per analyte. Table 3's Estimate column
    # holds the VARIANCE and the right-hand column the SD in ug/L; verified
    # arithmetically (sqrt(364) = 19.08, sqrt(53.3) = 7.30, sqrt(47.9) = 6.92,
    # sqrt(6.39) = 2.53 against the printed 19.1, 7.3, 6.9, 2.5). The exact
    # sqrt() of the variance is used rather than the rounded SD column.
    addSd <- sqrt(364)
    label("Additive residual SD for lidocaine (ug/L)") # Table 3 'Additive for Lidocaine': 364 (variance), SEE 81.30, 95% CI (204; 523), SD 19.1 ug/L
    addSd_megx <- sqrt(53.3)
    label("Additive residual SD for MEGX (ug/L)") # Table 3 'Additive for MEGX': 53.3 (variance), SEE 12.50, 95% CI (28.8; 77.8), SD 7.3 ug/L
    addSd_gx <- sqrt(47.9)
    label("Additive residual SD for GX (ug/L)") # Table 3 'Additive for GX': 47.9 (variance), SEE 23.80, 95% CI (1.2; 94.5), SD 6.9 ug/L
    addSd_xyl <- sqrt(6.39)
    label("Additive residual SD for 2,6-xylidine (ug/L)") # Table 3 'Additive for 2,6-xylidine': 6.39 (variance), SEE 1.41, 95% CI (3.63; 9.15), SD 2.5 ug/L
  })

  model({
    # 1. Binary covariate derivations. Each reproduces one Table 3 threshold
    # verbatim. Sect. 2.7.4: 'The laboratory safety parameters were used as
    # dichotomous covariates.'
    DLVL_HIGH <- (DLVL > 2)
    TBILI_HIGH <- (TBILI > 0.53)
    CRCL_LOW <- (CRCL <= 52.7)
    S1A2_IND <- (S1A2 == 1)
    BMI_HIGH <- (BMI > 27.9)
    ALT_HIGH <- (ALT > 11)
    LDH_HIGH <- (LDH > 195)

    # 2. Lidocaine metabolic-formation rate constants. Fixed, no IIV, no
    # covariates.
    k_megx_form <- exp(lk_megx_form)
    k_xyl_form <- exp(lk_xyl_form)
    # MEGX to GX. Estimated, no IIV, no covariates.
    k_gx_form <- exp(lk_gx_form)

    # 3. GX elimination rate constant. Sect. 2.7.4 eq. (3) selects the DLVL
    # stratum baseline, then the five binary covariate modifiers are ADDED on
    # the linear 1/h scale, then exponential IIV per Sect. 2.7.1 eq. (1).
    # No lower bound is imposed: the minimum attainable typical value over all
    # 32 covariate combinations is 1.44 - 0.526 - 0.32 - 0.492 = 0.102 /h,
    # which is positive. The vignette asserts this by enumeration.
    typ_kel_gx <-
      exp(lkel_gx_dlvlle2) * (1 - DLVL_HIGH) + exp(lkel_gx_dlvlgt2) * DLVL_HIGH +
      e_tbili_kel_gx * TBILI_HIGH +
      e_crcl_kel_gx * CRCL_LOW +
      e_s1a2_kel_gx * S1A2_IND +
      e_bmi_kel_gx * BMI_HIGH +
      e_alt_kel_gx * ALT_HIGH
    kel_gx <- typ_kel_gx * exp(etalkel_gx)

    # 4. 2,6-xylidine elimination rate constant. LDH stratum baseline plus the
    # single ALT modifier. Minimum attainable typical value is 0.41 /h.
    typ_kel_xyl <-
      exp(lkel_xyl_ldhle195) * (1 - LDH_HIGH) + exp(lkel_xyl_ldhgt195) * LDH_HIGH +
      e_alt_kel_xyl * ALT_HIGH
    kel_xyl <- typ_kel_xyl * exp(etalkel_xyl)

    # 5. Volumes. Lidocaine apparent central volume is DLVL-stratified and
    # carries exponential IIV; the three metabolite volumes are fixed.
    vc <- exp(lvc_dlvlle2 * (1 - DLVL_HIGH) + lvc_dlvlgt2 * DLVL_HIGH + etalvc)
    vc_megx <- exp(lvc_megx)
    vc_gx <- exp(lvc_gx)
    vc_xyl <- exp(lvc_xyl)

    # 6. ODE system. Sect. 3: 'This model is a linear model of four
    # compartments, one for each chemical entity. The ADVAN 5 subroutine in
    # NONMEM, which implements a user-defined general linear model, was used.'
    # The Fig. 1 legend enumerates exactly five transfers: k12 (1 -> 2), k14
    # (1 -> 4), k23 (2 -> 3), k30 (3 -> out) and k40 (4 -> out). There is no
    # k10 and no k20, so lidocaine leaves only by metabolism and MEGX leaves
    # only by conversion to GX. Dosing is a zero-order input into `central`
    # supplied through the event table's rate column (Sect. 2.7.7 assumption 1:
    # plasters 'deliver lidocaine at a constant rate (i.e., 1800 lg/h)'), so
    # there is no absorption compartment.
    d/dt(central) <- -(k_megx_form + k_xyl_form) * central
    d/dt(central_megx) <- k_megx_form * central - k_gx_form * central_megx
    d/dt(central_gx) <- k_gx_form * central_megx - kel_gx * central_gx
    d/dt(central_xyl) <- k_xyl_form * central - kel_xyl * central_xyl

    # 7. Observations. One additive residual per analyte (Sect. 2.7.1 eq. 2).
    Cc <- central / vc
    Cc_megx <- central_megx / vc_megx
    Cc_gx <- central_gx / vc_gx
    Cc_xyl <- central_xyl / vc_xyl

    Cc ~ add(addSd)
    Cc_megx ~ add(addSd_megx)
    Cc_gx ~ add(addSd_gx)
    Cc_xyl ~ add(addSd_xyl)
  })
}
