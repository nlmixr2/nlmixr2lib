Schlachter_2026_atogepant <- function() {
  description <- paste(
    "Three-compartment population pharmacokinetic model for oral atogepant, a",
    "calcitonin gene-related peptide (CGRP) receptor antagonist for the",
    "preventive treatment of migraine, pooled across 12 phase 1 studies, 1",
    "phase 2b/3 study and 1 phase 3 study in healthy participants and patients",
    "with episodic migraine (Schlachter 2026 'Phase 3 Model'). Absorption is",
    "sequential zero-order into the depot over a duration Tk0 followed by a",
    "first-order transfer ka into the central compartment, the whole input",
    "delayed by an absorption lag time; ka is not estimated independently but",
    "is derived from Tk0 and the estimated zero-order time fraction Fk0.",
    "Apparent clearance is lower in migraine patients than in healthy",
    "participants and is modified by severe hepatic impairment, itraconazole,",
    "quinidine and single- versus multiple-dose rifampicin; apparent central",
    "volume increases with body weight; relative bioavailability increases",
    "with dose and is modified by itraconazole and rifampicin; a high-fat meal",
    "lengthens the lag time; and the zero-order duration depends on dose and",
    "on the early phase 1 tablet formulation. Residual error is proportional",
    "with four study-specific magnitudes, and a blood-to-plasma ratio converts",
    "the plasma prediction to the dried-blood-sample matrix."
  )
  reference <- paste(
    "Schlachter L, Stodtmann S, Voelkner A, Jonsson F, Lagraauw HM,",
    "Boinpally RR. Population Pharmacokinetics of Atogepant for the",
    "Prevention of Migraine. Clin Pharmacokinet. 2026;65(2):151-166.",
    "doi:10.1007/s40262-025-01566-5."
  )
  vignette <- "Schlachter_2026_atogepant"

  # Non-canonical residual-SD names carried by this paper: the proportional
  # residual error magnitude is study-specific, with four distinct strata in
  # the final Phase 3 Model (Table 2 'sigma prop' rows). Same shape as
  # Pohl_2022_linzagolix.R, which carries two study-specific strata.
  paper_specific_residual_sds <- c(
    "propSdPhase1", "propSdCgpPk02", "propSdCgpMd01", "propSdPhase3"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-model effect on apparent central volume only:",
        "V1/F = 86.1 * (WT / 76.8)^0.411, printed as a display equation in",
        "Section 3.1. The 76.8 kg reference is the value printed in that",
        "equation; note it is NOT the Phase 3 Model median body weight",
        "(77.6 kg, Table 1) but the Phase 1 Model median/mean, carried",
        "forward as the centering constant. The exponent 0.411 is far below",
        "the allometric 1.0, which the Discussion calls out explicitly as a",
        "mild weight-volume relationship. Analysis-population range 40.7 to",
        "196 kg (Table 1). Body weight was screened on, but not retained",
        "for, CL/F."
      ),
      source_name        = "Body weight"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with episodic or chronic migraine)",
      notes              = paste(
        "Selects between the two apparent-clearance typical values printed",
        "side by side in Table 2 and in the Section 3.1 display equation:",
        "CL/F = 22.9 L/h in healthy participants and 17.4 L/h in patients.",
        "Patients are the reference because Section 2.3.1 states the typical",
        "value is parameterised for the level constituting the largest",
        "proportion of the population, and patients supply 1005 of the 1356",
        "Phase 3 Model subjects (Table S2: 463 phase 2b/3 plus 542 phase 3).",
        "Section 3.1 expresses the same contrast the other way round, as",
        "CL/F being 23.7% lower in patients than in healthy participants.",
        "No statistically significant PK difference was found between",
        "patients with episodic and with chronic migraine, so both carry 0."
      ),
      source_name        = "Healthy participant vs patient"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, or mild or moderate impairment)",
      notes              = paste(
        "Multiplicative effect on apparent clearance: CL/F * (1 - 0.366),",
        "a 36.6% reduction (Table 2 and the Section 3.1 display equation).",
        "Mild and moderate hepatic impairment were tested and were NOT",
        "retained, so those subjects sit in the reference category together",
        "with subjects having normal hepatic function; this is why the",
        "indicator is HEPIMP_SEV and not the graded HEPIMP. Eight subjects",
        "(0.6%) were in each of the mild, moderate and severe strata",
        "(Table 1), all from the dedicated hepatic-impairment study",
        "CGP-PK-01 (Table S1)."
      ),
      source_name        = "Severe hepatic impairment"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 / P-gp inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant itraconazole)",
      notes              = paste(
        "Two multiplicative effects, both printed in Table 2 and in the",
        "Section 3.1 display equations: CL/F * (1 - 0.662), a 66.2%",
        "clearance reduction, and Frel * (1 + 0.949), a 1.95-fold increase",
        "in relative bioavailability that the Discussion attributes to",
        "itraconazole-mediated CYP3A4 inhibition in enterocytes. Studied in",
        "the dedicated drug-drug-interaction study CGP-PK-02 (40 subjects,",
        "2.9% of the Phase 3 Model population)."
      ),
      source_name        = "Itraconazole"
    ),
    CONMED_RIFAMPICIN_SD = list(
      description        = "Concomitant rifampicin, single dose (OATP1B1 inhibition phase)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant rifampicin)",
      notes              = paste(
        "Rifampicin ('rifampin' in the source paper's US usage) enters this",
        "model as TWO mutually exclusive indicators because its single-dose",
        "and multiple-dose effects act through different mechanisms and in",
        "opposite directions. This indicator is the single-dose state, which",
        "the Discussion attributes to OATP1B1 inhibition: CL/F * (1 - 0.128)",
        "and Frel * (1 + 1.42), i.e. 12.8% lower clearance and a 2.4-fold",
        "higher relative bioavailability. See CONMED_RIFAMPICIN_MD for the",
        "induction state. A subject in the dedicated study CGP-PK-12 carries",
        "this indicator on the single-dose occasion and the MD indicator on",
        "the multiple-dose occasion, so both are per-dose-record rather than",
        "per-subject. 31 subjects (2.3%) received rifampicin (Table 1).",
        "NOTE: the Section 3.1 display equation for Frel prints the",
        "single-dose and multiple-dose labels swapped relative to Table 2,",
        "the Section 3.1 prose and the Abstract; the three-to-one majority",
        "is followed here. See the vignette Errata."
      ),
      source_name        = "Rifampin after first dose"
    ),
    CONMED_RIFAMPICIN_MD = list(
      description        = "Concomitant rifampicin, multiple doses (CYP3A4 induction phase)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant rifampicin)",
      notes              = paste(
        "The multiple-dose rifampicin state, which the Discussion attributes",
        "to CYP3A4 induction: CL/F * (1 + 0.818), a 1.82-fold increase, and",
        "Frel * (1 - 0.248), a 24.8% decrease. Opposite in direction to the",
        "single-dose state carried by CONMED_RIFAMPICIN_SD on both",
        "parameters, which is why a single CONMED_RIFAMPICIN indicator",
        "cannot express this model. The two indicators are mutually",
        "exclusive: a record with both set to 1 is not a state the source",
        "analysis contains."
      ),
      source_name        = "Rifampin following multiple doses"
    ),
    CONMED_QUINIDINE = list(
      description        = "Concomitant quinidine (P-glycoprotein inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant quinidine)",
      notes              = paste(
        "Multiplicative effect on apparent clearance only: CL/F * (1 -",
        "0.285), a 28.5% reduction (Table 2 and the Section 3.1 display",
        "equation). Unlike itraconazole and rifampicin, quinidine carries no",
        "effect on relative bioavailability in this model. Studied in the",
        "dedicated drug-drug-interaction study 3101-103-002 (25 subjects,",
        "1.8% of the Phase 3 Model population; Table 1)."
      ),
      source_name        = "Quinidine"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat meal at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Per-dose-record indicator. Single multiplicative effect on the",
        "absorption lag time: ALAG = 0.276 * (1 + 0.672), printed as a",
        "display equation in Section 3.1, which lengthens the lag from 0.276",
        "to 0.46 h. Section 3.1 states the high-fat meal 'otherwise had no",
        "impact', so food does not touch Tk0, Frel or clearance. The meal is",
        "specifically a high-fat meal in the dedicated food-effect study",
        "3101-105-002 (Table S1), so FED_HIGHFAT applies rather than the",
        "general FED."
      ),
      source_name        = "Food"
    ),
    DOSE_ATOGEPANT_MG = list(
      description        = "Administered atogepant dose level",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record dose amount in mg, driving two power-model effects",
        "printed as display equations in Section 3.1:",
        "Frel = (dose / 60 mg)^0.119 and Tk0 = 0.908 * (dose / 60 mg)^0.199.",
        "The 60 mg reference is the marketed once-daily dose. Section 3.1",
        "reports that Frel is approximately 1.24-fold higher at 60 mg than",
        "at 10 mg, which the exponent reproduces: (60/10)^0.119 = 1.24.",
        "Because the Tk0 effect propagates into ka through the derived",
        "relation ka = Fk0 / (Tk0 * (1 - Fk0)), Table 2 labels the same",
        "0.199 estimate as an 'Exponential dose effect on ka'; it is applied",
        "to Tk0. Studied dose range 10 to 300 mg (Table S2). Must be",
        "supplied as a data column separate from the event-table amt: a",
        "covariate column literally named DOSE (any casing) is consumed by",
        "rxode2's etTrans() and never reaches model()."
      ),
      source_name        = "Dose"
    ),
    FORM_ATOGEPANT_EARLYTAB = list(
      description        = "Early phase 1 atogepant tablet formulation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the Formulation 5 tablet used in the phase 2b/3 and phase 3 studies)",
      notes              = paste(
        "Per-dose-record indicator for the early phase 1 tablet",
        "presentations, with a single multiplicative effect on the",
        "zero-order absorption duration: Tk0 * (1 - 0.353), a 35% shorter",
        "duration (Section 3.1 display equation). The reference is",
        "unambiguous across all three source statements and is the",
        "Formulation 5 tablet, which supplies 1274 of 1356 subjects (94.0%)",
        "and 100% of the phase 2b/3, phase 3 and both external-validation",
        "cohorts (Table S2). The source is internally inconsistent about",
        "WHICH early formulations carry the effect: the display equation",
        "says 'phase 1 / formulations 3 and 4', the Section 3.1 prose says",
        "'the formulation 1 tablet used in early phase 1 studies', and the",
        "Table 2 row alias says 'Formulation 2 tablet / formulation 4",
        "tablet' (a label carried over from the Phase 1 Model column, where",
        "two separate estimates of -0.44 and -0.42 sat on ka). All three",
        "agree the effect belongs to the early phase 1 tablets and not to",
        "Formulation 5, and Formulations 1 to 4 together are only 82 of 1356",
        "subjects, so a single aggregate early-tablet indicator is used.",
        "See the vignette Errata. Because the Tk0 effect propagates into ka,",
        "Table 2 labels this estimate as a formulation effect 'on ka'; it is",
        "applied to Tk0."
      ),
      source_name        = "Formulation"
    ),
    STUDY_CGP_PK_02 = list(
      description        = "Phase 1 itraconazole drug-drug-interaction study CGP-PK-02 indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other contributing study)",
      notes              = paste(
        "Selects the study-specific proportional residual error magnitude",
        "0.228 (22.8% CV) reported in Table 2 as 'sigma prop (study",
        "CGP-PK-02)'. CGP-PK-02 is the 40-subject single-dose itraconazole",
        "drug-drug-interaction study and is the only contributing study that",
        "provided both plasma and dried-blood-sample concentrations",
        "(Table S1). Per-record study-fixed indicator; mutually exclusive",
        "with STUDY_CGP_MD_01 and STUDY_ATOGEPANT_PHASE3."
      ),
      source_name        = "Study"
    ),
    STUDY_CGP_MD_01 = list(
      description        = "Phase 2b/3 episodic-migraine study CGP-MD-01 indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other contributing study)",
      notes              = paste(
        "Selects the study-specific proportional residual error magnitude",
        "0.585 (58.5% CV) reported in Table 2 as 'sigma prop (study",
        "CGP-MD-01)'. CGP-MD-01 is the 463-subject phase 2b/3 dose-ranging",
        "efficacy study in patients with episodic migraine; it contributed",
        "1431 sparsely sampled dried-blood-sample observations (Tables S1",
        "and S2), which is why its residual error is the largest of the four",
        "strata. Per-record study-fixed indicator; mutually exclusive with",
        "STUDY_CGP_PK_02 and STUDY_ATOGEPANT_PHASE3."
      ),
      source_name        = "Study"
    ),
    STUDY_ATOGEPANT_PHASE3 = list(
      description        = "Phase 3 ADVANCE or PROGRESS study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other contributing study)",
      notes              = paste(
        "Selects the study-specific proportional residual error magnitude",
        "0.491 (49.1% CV) reported in Table 2 as 'sigma prop (ADVANCE and",
        "PROGRESS studies)'. One indicator covers both pivotal phase 3",
        "studies because Table 2 reports a single shared estimate for them:",
        "ADVANCE (3101-301-002, episodic migraine, in the model-development",
        "dataset) and PROGRESS (3101-303-002, chronic migraine, in the",
        "external-validation dataset). Per-record study-fixed indicator;",
        "mutually exclusive with STUDY_CGP_PK_02 and STUDY_CGP_MD_01. When",
        "all three study indicators are 0 the record falls in the reference",
        "stratum, 'all phase 1 studies except CGP-PK-02', whose residual",
        "error is 0.307 (30.7% CV)."
      ),
      source_name        = "Study"
    ),
    OCC = list(
      description        = "Dosing-occasion index for the inter-occasion variability on relative bioavailability",
      units              = "(integer)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer occasion column, 1 to 4, multiplexing the four",
        "inter-occasion-variability etas on relative bioavailability. The",
        "source reports the IOV variance on Frel (Table 2 'omega2 IOV Frel'",
        "= 0.0373, i.e. 19.3% CV) but does NOT state how many occasions the",
        "analysis dataset defined, so the occasion count of four here is an",
        "implementation choice, not a source value; one variance is reported",
        "and every occasion shares it. Pass OCC = 1 for single-occasion",
        "data so the first IOV eta applies. See the vignette Errata."
      ),
      source_name        = "Occasion"
    )
  )

  # Covariates the source screened but did not retain in the final Phase 3
  # Model. Figure 4 presents these as a forest plot of model-predicted
  # steady-state AUC24 and Cmax ratios versus their reference groups; every
  # one changed exposure by less than 20% and was judged not clinically
  # significant, so none appears in model(). Recorded here to preserve the
  # provenance of the covariate screen without carrying convention warnings.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and V1/F in the Phase 1 Model and evaluated in the",
        "Figure 4 forest plot with a reference group of < 65 years; the",
        "confidence interval included 1 (no effect) for both AUC24 and Cmax.",
        "Analysis-population range 18 to 78 years (Table 1)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Retained in the Phase 1 Model as a -0.19 fractional effect on CL/F",
        "(Table 2) but NOT retained in the Phase 3 Model. Evaluated in the",
        "Figure 4 forest plot against a male reference group with a < 20%",
        "exposure effect. 74.7% of the Phase 3 Model population was female",
        "(Table 1)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "The stepwise covariate search found Asian descent statistically",
        "significant on CL/F in the Phase 2 Model but it was removed for",
        "over-parameterisation (Supplement Section 1.2.3), and it was not",
        "retained in the Phase 3 Model. An exploratory analysis of predicted",
        "exposures for patients in Asian versus non-Asian countries in the",
        "PROGRESS validation showed no substantial regional difference",
        "(Section 3.3)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F in the Phase 1 Model. Mild and moderate renal",
        "impairment were predicted to have no relevant effect on atogepant",
        "PK (Figure 4 and the Discussion). Phase 3 Model median 127 mL/min",
        "(range 46.5 to 392); 87.2% normal, 11.9% mild and 1.0% moderate",
        "renal impairment (Table 1)."
      )
    ),
    CONMED_STATIN = list(
      description = "Concomitant statin use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on both CL/F and Frel in the Phase 3 Model covariate search",
        "(Table S3) and found significant on neither (Section 3.1).",
        "Figure 4 puts the effect at a 14% AUC24 increase and a 13% Cmax",
        "increase, below the 20% clinical-relevance threshold. 4.9% of the",
        "Phase 3 Model population (Table 1)."
      )
    ),
    CONMED_BCRP_INHIBITOR = list(
      description = "Concomitant breast cancer resistance protein (BCRP) inhibitor use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Evaluated only in the Figure 4 forest plot, against a no-BCRP-",
        "inhibitor reference. Expected to change AUC24 by a 1% decrease with",
        "no change in Cmax, and the confidence interval included 1."
      )
    ),
    CONMED_BCRP_SUBSTRATE = list(
      description = "Concomitant breast cancer resistance protein (BCRP) substrate use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Evaluated only in the Figure 4 forest plot, against a",
        "no-BCRP-substrate reference. Expected to change AUC24 by a 3%",
        "decrease and Cmax by a 3% decrease, and the confidence interval",
        "included 1."
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "atogepant", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "atogepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "atogepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "atogepant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1356,
    n_studies      = 14,
    age_range      = "18-78 years",
    age_median     = "39.0 years",
    weight_range   = "40.7-196 kg",
    weight_median  = "77.6 kg",
    height_median  = "166 cm (range 146-204)",
    sex_female_pct = 74.7,
    race_ethnicity = c(
      Caucasian = 76.0, `Black/African American` = 17.9, Asian = 3.5,
      Multiple = 1.9, `Native American/Alaska Native` = 0.3,
      `Pacific Islander` = 0.2, `Not used` = 0.1
    ),
    disease_state  = paste(
      "351 healthy participants and 1005 patients with episodic migraine.",
      "Hepatic function: 98.2% none, 0.6% mild, 0.6% moderate, 0.6% severe.",
      "Renal function per Cockcroft-Gault: 87.2% normal (> 90 mL/min),",
      "11.9% mild (60-89), 1.0% moderate (30-59); no severe impairment."
    ),
    renal_function = "median creatinine clearance 127 mL/min (range 46.5-392)",
    co_medication  = paste(
      "Itraconazole 40 subjects (2.9%), rifampicin 31 (2.3%), quinidine 25",
      "(1.8%), concomitant statin 67 (4.9%); 1260 (92.9%) had no",
      "cotreatment (Table 1)."
    ),
    dose_range     = "10-300 mg oral tablet, once or twice daily, single and multiple dose",
    notes          = paste(
      "Baseline characteristics from Schlachter 2026 Table 1, Phase 3 Model",
      "column. 11,763 observations from 1356 participants across 12 phase 1",
      "studies, the phase 2b/3 study CGP-MD-01 and the phase 3 ADVANCE study",
      "(Tables S1 and S2). Concentrations below the assay lower limit of",
      "quantitation (0.1, 1.0 or 10 ng/mL for plasma and 2.5 ng/mL for",
      "dried blood samples, by study) were treated as missing and excluded,",
      "as were observations with |CWRES| > 6 and concentrations flagged by",
      "the dose-normalised outlier rules in Section 2.2. Study MK-8031-P001,",
      "which used an oral solution, was in the Phase 2 Model but was dropped",
      "from the Phase 3 Model because of its formulation; this is why the",
      "Phase 2 Model's solution effects on Frel and on lag time have no",
      "Phase 3 counterpart. The model was externally validated without",
      "re-estimation against 350 patients with chronic migraine (PROGRESS,",
      "1542 observations) and 127 patients with episodic migraine for whom",
      "two to four classes of conventional oral preventives had failed",
      "(ELEVATE, 562 observations); those 477 subjects are not counted in",
      "n_subjects above."
    )
  )

  ini({
    # ---- Structural disposition parameters (Schlachter 2026 Table 2,
    #      'Phase 3 Model' Estimate column). Volumes and clearances are
    #      apparent (X/F) because every study was oral.
    lcl  <- log(17.4); label("Apparent clearance CL/F in patients with migraine (L/h)")       # Table 2 'Apparent clearance patients' 17.4 (RSE 2%, 95% CI 16.8-18.1)
    lvc  <- log(86.1); label("Apparent central volume V1/F at 76.8 kg (L)")                   # Table 2 'Apparent central volume of distribution' 86.1 (RSE 2.3%, 95% CI 82.2-89.9)
    lq   <- log(1.43); label("Apparent first intercompartmental clearance Q/F (L/h)")         # Table 2 'Apparent first intercompartmental clearance' 1.43 (RSE 10%, 95% CI 1.15-1.71)
    lvp  <- log(40.5); label("Apparent first peripheral volume V2/F (L)")                     # Table 2 'Apparent first peripheral volume of distribution' 40.5 (RSE 7.3%, 95% CI 34.7-46.3)
    lq2  <- log(1.68); label("Apparent second intercompartmental clearance Q2/F (L/h)")       # Table 2 'Apparent second intercompartmental clearance' 1.68 (RSE 7.4%, 95% CI 1.44-1.93)
    lvp2 <- log(13.0); label("Apparent second peripheral volume V3/F (L)")                    # Table 2 'Apparent second peripheral volume of distribution' 13.0 (RSE 11.9%, 95% CI 9.96-16.0)

    # ---- Absorption. The input is sequential: a zero-order release into the
    #      depot lasting Tk0, then a first-order transfer ka out of the depot,
    #      the whole thing delayed by the lag time. ka is NOT an independent
    #      parameter -- Section 3.1 states it 'was linked to the zero-order
    #      input parameters through ka = Fk0 / [Tk0 * (1 - Fk0)]', so it is
    #      derived in model() and moves with each subject's Tk0.
    ld1   <- log(0.908); label("Duration of the zero-order input into the depot, Tk0, at 60 mg (h)")  # Table 2 'Duration zero-order absorption' 0.908 (RSE 4.7%, 95% CI 0.824-0.992)
    ltlag <- log(0.276); label("Absorption lag time ALAG when fasted (h)")                            # Table 2 'Lag time' 0.276 (RSE 0.5%, 95% CI 0.273-0.279)

    # Fk0 is the fraction of the total absorption time constant (Tk0 + 1/ka)
    # contributed by the zero-order step, NOT a fraction of the dose: the
    # printed relation ka = Fk0 / [Tk0 * (1 - Fk0)] rearranges exactly to
    # Fk0 = Tk0 / (Tk0 + 1/ka). It is bounded in (0, 1) and so is not
    # log-transformed. Check against the source: 0.908 / (0.908 + 1/2.486)
    # = 0.693, and Fk0 / [Tk0 * (1 - Fk0)] = 0.693 / (0.908 * 0.307) = 2.486,
    # reproducing the derived ka of 2.48/h printed in Section 3.1.
    fk0 <- 0.693; label("Fraction of the total absorption time constant that is zero-order, Fk0 (fraction)")  # Table 2 'Fraction zero-order absorption (Fk0)' 0.693 (RSE 2.9%, 95% CI 0.654-0.732)

    # Relative bioavailability is purely relative and is anchored at 1 for the
    # 60 mg reference dose with no interacting comedication; the whole
    # dose-dependence and all comedication effects are carried by the
    # e_*_fdepot terms below. Section 3.1 prints Frel as (dose/60 mg)^0.119
    # times the comedication factors, with no separate estimated scale term.
    lfdepot <- fixed(log(1)); label("Relative bioavailability Frel at the 60 mg reference dose (fraction)")  # Schlachter 2026 Section 3.1 Frel display equation

    # Blood-to-plasma ratio converting the plasma prediction to the
    # dried-blood-sample matrix used by studies CGP-PK-02 and CGP-MD-01.
    bpratio <- 0.573; label("Blood-to-plasma concentration ratio (fraction)")  # Table 2 'Blood-plasma ratio' 0.573 (RSE 2%, 95% CI 0.550-0.596); Section 3.1 'blood concentrations were predicted to be 57.3% of those in plasma'

    # ---- Covariate effects on apparent clearance. Section 2.3.1 gives the
    #      categorical form as P_TV * (1 + theta_Xm), so each estimate is a
    #      fractional change: -0.662 is a 66.2% reduction.
    e_dis_healthy_cl          <-  0.3161; label("Fractional change in CL/F for a healthy participant vs a migraine patient (unitless)")  # Derived from the two Table 2 typical values: 22.9 / 17.4 - 1 = 0.3161, so 17.4 * 1.3161 = 22.9 L/h. Section 3.1 states the same contrast as patients being 23.7% lower than healthy
    e_hepimp_sev_cl           <- -0.366;  label("Fractional change in CL/F with severe hepatic impairment (unitless)")                   # Table 2 'Effect of severe hepatic impairment on CL/F' -0.366 (RSE 19%, 95% CI -0.503 to -0.230)
    e_conmed_itraconazole_cl  <- -0.662;  label("Fractional change in CL/F with concomitant itraconazole (unitless)")                    # Table 2 'Itraconazole effect on CL/F' -0.662 (RSE 0.6%, 95% CI -0.669 to -0.654)
    e_conmed_quinidine_cl     <- -0.285;  label("Fractional change in CL/F with concomitant quinidine (unitless)")                       # Table 2 'Quinidine effect on CL/F' -0.285 (RSE 3.4%, 95% CI -0.305 to -0.266)
    e_conmed_rifampicin_sd_cl <- -0.128;  label("Fractional change in CL/F with a single rifampicin dose (unitless)")                    # Table 2 'Rifampin effect on CL/F after first dose' -0.128 (RSE 14.2%, 95% CI -0.164 to -0.0924)
    e_conmed_rifampicin_md_cl <-  0.818;  label("Fractional change in CL/F with multiple rifampicin doses (unitless)")                   # Table 2 'Rifampin effect on CL/F following multiple doses' 0.818 (RSE 3.9%, 95% CI 0.756-0.881); Section 3.1 calls this a 1.82-fold increase

    # ---- Covariate effect on apparent central volume. Continuous covariates
    #      take the power form of Section 2.3.1, centred at the reference
    #      value printed in the Section 3.1 display equation.
    e_wt_vc <- 0.411; label("Power exponent on (WT / 76.8 kg) for V1/F (unitless)")  # Table 2 'Exponential weight effect on V1/F' 0.411 (RSE 11.7%, 95% CI 0.317-0.505)

    # ---- Covariate effects on relative bioavailability.
    e_dose_atogepant_mg_fdepot    <-  0.119; label("Power exponent on (dose / 60 mg) for Frel (unitless)")                  # Table 2 'Exponential dose effect on Frel' 0.119 (RSE 11.2%, 95% CI 0.0928-0.145)
    e_conmed_itraconazole_fdepot  <-  0.949; label("Fractional change in Frel with concomitant itraconazole (unitless)")    # Table 2 'Itraconazole effect on Frel' 0.949 (RSE 11.4%, 95% CI 0.737-1.16); Section 3.1 calls this a 1.95-fold increase
    e_conmed_rifampicin_sd_fdepot <-  1.42;  label("Fractional change in Frel with a single rifampicin dose (unitless)")    # Table 2 'Rifampin effect on Frel after first dose' 1.42 (RSE 11%, 95% CI 1.12-1.73); Section 3.1 and the Abstract both call this a 2.4-fold increase
    e_conmed_rifampicin_md_fdepot <- -0.248; label("Fractional change in Frel with multiple rifampicin doses (unitless)")   # Table 2 'Rifampin effect on Frel following multiple doses' -0.248 (RSE 15.7%, 95% CI -0.325 to -0.172); Section 3.1 calls this a 24.8% decrease

    # ---- Covariate effect on the absorption lag time.
    e_fed_highfat_tlag <- 0.672; label("Fractional change in the lag time when dosed after a high-fat meal (unitless)")  # Table 2 'Food effect on ALAG' 0.672 (RSE 2.1%, 95% CI 0.644-0.699); 0.276 * 1.672 = 0.46 h, matching Section 3.1

    # ---- Covariate effects on the zero-order duration. Table 2 labels both of
    #      these as effects 'on ka' because they propagate into the derived ka;
    #      the Section 3.1 display equation applies them to Tk0.
    e_dose_atogepant_mg_d1       <-  0.199; label("Power exponent on (dose / 60 mg) for Tk0 (unitless)")             # Table 2 'Exponential dose effect on ka' 0.199 (RSE 9%, 95% CI 0.164-0.234); Section 3.1 display equation applies it to Tk0
    e_form_atogepant_earlytab_d1 <- -0.353; label("Fractional change in Tk0 for the early phase 1 tablet (unitless)")  # Table 2 'Formulation 2 tablet effect on ka/formulation 4 tablet effect on ka' -0.353 (RSE 22.3%, 95% CI -0.507 to -0.198); Section 3.1 display equation applies it to Tk0 and the prose calls it a 35% shorter duration

    # ---- Interindividual variability. Table 2 reports VARIANCES; the Table 1
    #      footnote states 'IIV is derived from variance according to
    #      sqrt(omega^2) * 100', and each value below reproduces the percent
    #      CV quoted in Section 3.1, which pins the scale.
    etalcl     ~ 0.0243  # Table 2 'omega2 CL/F' 0.0243 (RSE 9.9%, 95% CI 0.0196-0.0290); sqrt = 0.156, the 15.6% CV of Section 3.1
    etald1     ~ 0.163   # Table 2 'omega2 tk0' 0.163 (RSE 5.8%, 95% CI 0.144-0.182); sqrt = 0.404, the 40.4% CV of Section 3.1
    etalfdepot ~ 0.234   # Table 2 'omega2 Frel' 0.234 (RSE 5%, 95% CI 0.211-0.257); sqrt = 0.484, the 48.4% CV of Section 3.1

    # Q/F and V2/F share an estimated covariance (Table 2 'Cov Q/F, V2/F').
    # The implied correlation is 0.422 / sqrt(0.471 * 0.601) = 0.79.
    etalq + etalvp ~ c(0.471,
                       0.422, 0.601)  # Table 2 'omega2 Q/F' 0.471 (RSE 17.1%), 'Cov Q/F, V2/F' 0.422 (RSE 14.6%), 'omega2 V2/F' 0.601 (RSE 14.5%); sqrt of the diagonals gives the 68.6% and 77.5% CVs of Section 3.1

    # Inter-occasion variability on relative bioavailability. One variance is
    # reported and it is shared by every occasion; nlmixr2 has no NONMEM
    # $OMEGA BLOCK(1) SAME shortcut, so occasion 1 carries the estimated
    # variance and the later occasions fix it to the same value. The source
    # does not state how many occasions its dataset defined -- see the
    # OCC covariate note and the vignette Errata.
    etaiov_fdepot_1 ~ 0.0373       # Table 2 'omega IOV Frel' 0.0373 (RSE 3.5%, 95% CI 0.0347-0.0399); sqrt = 0.193, the 19.3% CV of Section 3.1
    etaiov_fdepot_2 ~ fixed(0.0373)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_3 ~ fixed(0.0373)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_4 ~ fixed(0.0373)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual unexplained variability: proportional, with four
    #      study-specific magnitudes. Table 2 reports these on the SD scale
    #      directly -- each value equals the percent CV quoted in Section 3.1
    #      divided by 100, so no square root is taken.
    propSdPhase1  <- 0.307; label("Proportional residual SD for all phase 1 studies except CGP-PK-02 (fraction)")  # Table 2 'sigma prop' 0.307 (RSE 0.5%, 95% CI 0.304-0.310); Section 3.1 quotes 30.7% CV
    propSdCgpPk02 <- 0.228; label("Proportional residual SD for phase 1 study CGP-PK-02 (fraction)")               # Table 2 'sigma prop (study CGP-PK-02)' 0.228 (RSE 1.4%, 95% CI 0.222-0.235); Section 3.1 quotes 22.8% CV
    propSdCgpMd01 <- 0.585; label("Proportional residual SD for phase 2b/3 study CGP-MD-01 (fraction)")            # Table 2 'sigma prop (study CGP-MD-01)' 0.585 (RSE 3.1%, 95% CI 0.550-0.620); Section 3.1 quotes 58.5% CV
    propSdPhase3  <- 0.491; label("Proportional residual SD for the ADVANCE and PROGRESS phase 3 studies (fraction)")  # Table 2 'sigma prop (ADVANCE and PROGRESS studies)' 0.491 (RSE 2%, 95% CI 0.472-0.510); Section 3.1 quotes 49.1% CV
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators to
    #    multiplex the four inter-occasion-variability etas on relative
    #    bioavailability. For single-occasion data pass OCC = 1.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <-
      oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # 2. Apparent clearance. The typical value is the migraine-patient value;
    #    healthy participants scale up to the 22.9 L/h of Table 2. The five
    #    remaining terms are the intrinsic- and extrinsic-factor effects of
    #    the Section 3.1 CL/F display equation. The two rifampicin indicators
    #    are mutually exclusive and pull in opposite directions.
    cl <- exp(lcl + etalcl) *
      (1 + e_dis_healthy_cl * DIS_HEALTHY) *
      (1 + e_hepimp_sev_cl * HEPIMP_SEV) *
      (1 + e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE) *
      (1 + e_conmed_quinidine_cl * CONMED_QUINIDINE) *
      (1 + e_conmed_rifampicin_sd_cl * CONMED_RIFAMPICIN_SD) *
      (1 + e_conmed_rifampicin_md_cl * CONMED_RIFAMPICIN_MD)

    # 3. Remaining disposition parameters. Body weight touches V1/F only, and
    #    the IIV on Q/F and V2/F is correlated.
    vc  <- exp(lvc) * (WT / 76.8)^e_wt_vc
    q   <- exp(lq  + etalq)
    vp  <- exp(lvp + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # 4. Absorption. Dose and formulation act on the zero-order duration; the
    #    first-order rate constant is then derived from it, so both effects
    #    propagate into ka exactly as Table 2's row labels describe. A high-fat
    #    meal acts on the lag time alone.
    d1 <- exp(ld1 + etald1) *
      (DOSE_ATOGEPANT_MG / 60)^e_dose_atogepant_mg_d1 *
      (1 + e_form_atogepant_earlytab_d1 * FORM_ATOGEPANT_EARLYTAB)
    ka   <- fk0 / (d1 * (1 - fk0))
    tlag <- exp(ltlag) * (1 + e_fed_highfat_tlag * FED_HIGHFAT)

    # 5. Relative bioavailability: a power function of dose plus the
    #    itraconazole and rifampicin factors, carrying both IIV and IOV.
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) *
      (DOSE_ATOGEPANT_MG / 60)^e_dose_atogepant_mg_fdepot *
      (1 + e_conmed_itraconazole_fdepot * CONMED_ITRACONAZOLE) *
      (1 + e_conmed_rifampicin_sd_fdepot * CONMED_RIFAMPICIN_SD) *
      (1 + e_conmed_rifampicin_md_fdepot * CONMED_RIFAMPICIN_MD)

    # 6. Micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 7. ODE system: three disposition compartments with linear elimination
    #    from the central compartment. Dose the depot with rate = -2 so
    #    rxode2 applies the modelled dur().
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    dur(depot)  <- d1
    alag(depot) <- tlag
    f(depot)    <- fdepot

    # 8. Observation. Volumes are apparent (V/F) and doses are in mg, so
    #    central/vc is mg/L = ug/mL; the source reports atogepant in ng/mL
    #    (plasma assay LLOQ 0.1 to 10 ng/mL by study), hence the factor 1000.
    #    Cb is the same drug measured in the dried-blood-sample matrix used
    #    by studies CGP-PK-02 and CGP-MD-01. It is carried as a derived
    #    output column rather than as a second estimated endpoint: the source
    #    applies the SAME study-specific proportional error to plasma and to
    #    DBS records (Table 2 strata are by study, not by matrix), so a
    #    separate endpoint would add no information while forcing every
    #    observation row in every downstream simulation to carry a dvid.
    Cc <- central / vc * 1000
    Cb <- Cc * bpratio

    # The proportional residual error magnitude is study-specific. When all
    # three study indicators are 0 the record falls in the reference stratum,
    # 'all phase 1 studies except CGP-PK-02'.
    propSdInd <-
      propSdPhase1  * (1 - STUDY_CGP_PK_02 - STUDY_CGP_MD_01 - STUDY_ATOGEPANT_PHASE3) +
      propSdCgpPk02 * STUDY_CGP_PK_02 +
      propSdCgpMd01 * STUDY_CGP_MD_01 +
      propSdPhase3  * STUDY_ATOGEPANT_PHASE3

    Cc ~ prop(propSdInd)
  })
}
