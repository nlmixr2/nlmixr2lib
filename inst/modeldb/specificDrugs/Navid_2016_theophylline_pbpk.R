Navid_2016_theophylline_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 83-ODE, four coupled compounds).",
    "Theophylline (THP) disposition and its drug-drug interactions with",
    "caffeine (CAF) and ciprofloxacin (CIP) in adults, with paraxanthine",
    "(PX) carried as a fourth species because it is both the major",
    "caffeine metabolite and a competitive inhibitor of theophylline",
    "metabolism. Each compound moves through twelve perfusion-limited",
    "organs plus arterial and venous blood; the three orally dosed",
    "compounds each carry a nine-state Yu and Amidon compartmental",
    "absorption and transit chain (stomach, seven small-intestine transit",
    "compartments, colon). Hepatic metabolism is Michaelis-Menten over",
    "fifteen enzyme-substrate reactions on CYP1A2, CYP2E1 and CYP3A4,",
    "with every reaction competitively inhibited by every other reaction",
    "sharing its enzyme -- that coupling is what produces the paper's",
    "drug-drug interactions, including the caffeine-to-theophylline",
    "metabolic route that makes coffee a theophylline source. Renal",
    "elimination is a fixed extraction ratio on kidney blood flow.",
    "Population variability is deterministic and enters through four",
    "covariates rather than random effects: SEXF switches the whole organ",
    "volume and blood-flow table between the ICRP reference man and woman,",
    "HEPFUNC_REL is the paper's Delta(met) hepatic-activity multiplier",
    "(0.66 in the average elderly), RENALFUNC_REL is its Delta(ren) renal",
    "multiplier (0.31 in the average elderly), and BODYFAT_PCT carries the",
    "adipose scaling the paper uses for its race comparison. The model has",
    "no between-subject variability and no residual-error model because",
    "the source reports neither; it was solved in Mathematica and",
    "calibrated against published concentration-time data rather than",
    "fitted. All values are the authors' Supplementary Tables 1-3.",
    sep = " "
  )
  reference <- paste(
    "Navid A, Ng DM, Wong SE, Lightstone FC.",
    "Application of a Physiologically Based Pharmacokinetic Model to Study",
    "Theophylline Metabolism and Its Interactions With Ciprofloxacin and",
    "Caffeine. CPT Pharmacometrics Syst Pharmacol. 2016;5(2):74-81.",
    "doi:10.1002/psp4.12061. PMC4761233.",
    "Equations 1-4 (the organ mass balance, the dosing term, bulk renal",
    "extraction and the inhibited Michaelis-Menten hepatic term) are from",
    "the main-text Methods. Every numerical value is from the",
    "Supplementary Materials: section 1 and the absorption equations from",
    "PSP4-5-74-s001.docx, organ volumes and blood flows from Table 1",
    "(PSP4-5-74-s002.docx), tissue:plasma partition coefficients from",
    "Table 2 (PSP4-5-74-s003.docx), and renal extraction ratios, fractions",
    "absorbed, transit times and all Vmax / Km values from Table 3",
    "(PSP4-5-74-s004.docx).",
    sep = " "
  )
  vignette <- "Navid_2016_theophylline_pbpk"

  # The supplement's parameter names are embedded MathType images rather
  # than text, so the .docx converts to a table of values with no row
  # labels. The names were recovered by carving the PDF payload out of
  # each EMF and reading it; three independent internal checks confirm the
  # row-to-name alignment, and each is reproduced in the vignette:
  #   * Table 3's last six rows are Tge and Tsi for the three oral drugs,
  #     and come out 30, 30, 30, 199, 199, 199 -- exactly the values
  #     supplement section 1 states in prose.
  #   * The two ciprofloxacin Km values land on 4.6e-6 and 5.1e-5, which
  #     are Fuhr's 0.18 mM and McLellan's 2 mM each divided by 39, the
  #     operation supplement note 9 describes in prose.
  #   * Supplement note 1 says the caffeine renal extraction ratio was
  #     reused for paraxanthine, and rows 2 and 4 are both 0.008.
  # The units COLUMN of Table 3 is the one thing the conversion does not
  # preserve: its four blocks have sizes 14 / 15 / 6 / 7 where the
  # parameter types require 7 / 14 / 15 / 6, i.e. the column is rotated by
  # seven rows because the leading dimensionless (blank) cells were
  # dropped. Units below are therefore taken from the parameter type, not
  # from that column. See the vignette Errata.

  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    SEXF = list(
      description = paste(
        "Biological sex. Switches the entire organ volume and blood-flow",
        "table between the men and women columns of Supplementary Table 1,",
        "and switches the reference percent body fat used to scale adipose",
        "volume. This is the paper's dominant source of between-subject",
        "difference: women are about 30 percent lighter, so an unadjusted",
        "milligram dose gives a higher exposure."
      ),
      source_name = "sex",
      units = "(binary)",
      type = "binary",
      reference = "0 (male)"
    ),
    HEPFUNC_REL = list(
      description = paste(
        "Delta(met) of Methods Eq 4: hepatic metabolic activity as a",
        "fraction of the reference group (healthy adults in their 20s).",
        "Multiplies every hepatic Michaelis-Menten velocity. The paper uses",
        "0.66 for the average elderly patient and reports safe thresholds",
        "as low as 0.35 in Table 1."
      ),
      source_name = "Delta_met",
      units = "(dimensionless)",
      type = "continuous",
      reference = "1 (reference group)"
    ),
    RENALFUNC_REL = list(
      description = paste(
        "Delta(ren) of the paper's Results: renal clearance as a fraction",
        "of the reference group. Multiplies the renal extraction ratio of",
        "every compound. The paper uses 0.31 for the average elderly",
        "patient, which matters far more for ciprofloxacin (extraction",
        "ratio 0.39) than for theophylline (0.004)."
      ),
      source_name = "Delta_ren",
      units = "(dimensionless)",
      type = "continuous",
      reference = "1 (reference group)"
    ),
    BODYFAT_PCT = list(
      description = paste(
        "Percent total body fat, used to scale the adipose compartment",
        "volume by the paper's Phi factor (supplement section 5). At the",
        "reference value for the subject's sex (24.932 percent for men,",
        "37.5 percent for women, both derived from Supplementary Table 1)",
        "the adipose volume is exactly the tabulated one. The paper's own",
        "per-race values come from Carpenter 2013, which was not available",
        "when this model was built; see the vignette Errata."
      ),
      source_name = "Phi",
      units = "% (percent)",
      type = "continuous",
      reference = "24.932 (men) / 37.5 (women)"
    )
  )

  # Every perfused-organ and blood state holds a CONCENTRATION in mmol/mL
  # (= M), which is the unit Eq 1 is written in. The luminal transit states
  # and the urine states hold AMOUNTS in mmol: the Yu and Amidon equations
  # of supplement section 1 are written in concentration but the lumen
  # volumes are never reported, and because every term of that chain is
  # first order the amount form is identical and is what conserves mass.
  compartmentData <- list(
    stomach_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit1_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit2_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit3_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit4_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit5_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit6_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    transit7_thp = list(analyte = "theophylline", units = "mmol", specimen = "administration site", verified = TRUE),
    colon_thp = list(analyte = "theophylline", units = "mmol", specimen = "faeces", verified = TRUE),
    adipose_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    bone_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    brain_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    gut_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    heart_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    kidney_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    liver_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    lung_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    muscle_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    other_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    skin_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    spleen_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    arterial_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    venous_thp = list(analyte = "theophylline", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    urine_thp = list(analyte = "theophylline", units = "mmol", specimen = "urine", verified = TRUE),
    stomach_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit1_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit2_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit3_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit4_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit5_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit6_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    transit7_caf = list(analyte = "caffeine", units = "mmol", specimen = "administration site", verified = TRUE),
    colon_caf = list(analyte = "caffeine", units = "mmol", specimen = "faeces", verified = TRUE),
    adipose_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    bone_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    brain_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    gut_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    heart_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    kidney_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    liver_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    lung_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    muscle_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    other_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    skin_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    spleen_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    arterial_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    venous_caf = list(analyte = "caffeine", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    urine_caf = list(analyte = "caffeine", units = "mmol", specimen = "urine", verified = TRUE),
    adipose_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    bone_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    brain_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    gut_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    heart_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    kidney_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    liver_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    lung_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    muscle_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    other_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    skin_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    spleen_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    arterial_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    venous_px = list(analyte = "paraxanthine", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    urine_px = list(analyte = "paraxanthine", units = "mmol", specimen = "urine", verified = TRUE),
    stomach_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit1_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit2_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit3_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit4_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit5_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit6_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    transit7_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "administration site", verified = TRUE),
    colon_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "faeces", verified = TRUE),
    adipose_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    bone_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    brain_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    gut_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    heart_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    kidney_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    liver_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    lung_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    muscle_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    other_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    skin_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    spleen_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "tissue", verified = TRUE),
    arterial_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    venous_cip = list(analyte = "ciprofloxacin", units = "mmol/mL", specimen = "whole blood", verified = TRUE),
    urine_cip = list(analyte = "ciprofloxacin", units = "mmol", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "adults; a reference group in their 20s and an average elderly group",
    disease_state = "healthy adults; theophylline is indicated for obstructive airway disease",
    dose_range = paste(
      "theophylline 125-1600 mg/day oral (the dose initiation and titration",
      "protocol runs 300 -> 400 -> 600 mg/day divided every 6-8 h);",
      "caffeine 183-960 mg/day oral; ciprofloxacin 100-1500 mg/day oral"
    ),
    notes = paste(
      "No subjects were fitted. This is a deterministic PBPK calibrated",
      "against published concentration-time curves from the literature",
      "(supplement Figure S2 reproduces theophylline data from Rovei 1982,",
      "caffeine from six sources and ciprofloxacin from three), so there is",
      "no study population, no n and no estimated variance. The reference",
      "group is the ICRP Publication 89 reference adult male and female."
    )
  )

  ini({
    # ==================================================================
    # Physiological parameters -- Supplementary Table 1 (file
    # PSP4-5-74-s002.docx), sourced by that table to ICRP Publication 89
    # (Valentin 2002). The twelve organ volumes plus arterial and venous
    # blood sum to exactly 73,000 mL (men) and 60,000 mL (women), i.e. the
    # ICRP reference adult male and female; and the organ blood flows sum
    # to exactly the lung flow (6,500 and 5,900 mL/min) once the portal
    # inflow is removed. Both identities are used in the vignette as
    # transcription checks.
    # ==================================================================
    v_adipose_male <- fixed(18200)
    label("Adipose volume, men (mL)")  # Table 1, row Adipose, column 'Volume men'
    v_adipose_female <- fixed(22500)
    label("Adipose volume, women (mL)")  # Table 1, row Adipose, column 'Volume women'
    v_bone_male <- fixed(5500)
    label("Bone volume, men (mL)")  # Table 1, row Bone, column 'Volume men'
    v_bone_female <- fixed(4000)
    label("Bone volume, women (mL)")  # Table 1, row Bone, column 'Volume women'
    v_brain_male <- fixed(1450)
    label("Brain volume, men (mL)")  # Table 1, row Brain, column 'Volume men'
    v_brain_female <- fixed(1300)
    label("Brain volume, women (mL)")  # Table 1, row Brain, column 'Volume women'
    v_gut_male <- fixed(2070)
    label("Gut volume, men (mL)")  # Table 1, row Gut, column 'Volume men'
    v_gut_female <- fixed(1930)
    label("Gut volume, women (mL)")  # Table 1, row Gut, column 'Volume women'
    v_heart_male <- fixed(330)
    label("Heart volume, men (mL)")  # Table 1, row Heart, column 'Volume men'
    v_heart_female <- fixed(250)
    label("Heart volume, women (mL)")  # Table 1, row Heart, column 'Volume women'
    v_kidney_male <- fixed(310)
    label("Kidney volume, men (mL)")  # Table 1, row Kidney, column 'Volume men'
    v_kidney_female <- fixed(275)
    label("Kidney volume, women (mL)")  # Table 1, row Kidney, column 'Volume women'
    v_liver_male <- fixed(1800)
    label("Liver volume, men (mL)")  # Table 1, row Liver, column 'Volume men'
    v_liver_female <- fixed(1400)
    label("Liver volume, women (mL)")  # Table 1, row Liver, column 'Volume women'
    v_lung_male <- fixed(500)
    label("Lung volume, men (mL)")  # Table 1, row Lung, column 'Volume men'
    v_lung_female <- fixed(420)
    label("Lung volume, women (mL)")  # Table 1, row Lung, column 'Volume women'
    v_muscle_male <- fixed(29000)
    label("Muscle volume, men (mL)")  # Table 1, row Muscle, column 'Volume men'
    v_muscle_female <- fixed(17500)
    label("Muscle volume, women (mL)")  # Table 1, row Muscle, column 'Volume women'
    v_other_male <- fixed(4790)
    label("Other volume, men (mL)")  # Table 1, row Other, column 'Volume men'
    v_other_female <- fixed(3895)
    label("Other volume, women (mL)")  # Table 1, row Other, column 'Volume women'
    v_skin_male <- fixed(3300)
    label("Skin volume, men (mL)")  # Table 1, row Skin, column 'Volume men'
    v_skin_female <- fixed(2300)
    label("Skin volume, women (mL)")  # Table 1, row Skin, column 'Volume women'
    v_spleen_male <- fixed(150)
    label("Spleen volume, men (mL)")  # Table 1, row Spleen, column 'Volume men'
    v_spleen_female <- fixed(130)
    label("Spleen volume, women (mL)")  # Table 1, row Spleen, column 'Volume women'
    v_arterial_male <- fixed(1867)
    label("Arterial blood volume, men (mL)")  # Table 1, row Blood -Arterial, column 'Volume men'
    v_arterial_female <- fixed(1367)
    label("Arterial blood volume, women (mL)")  # Table 1, row Blood -Arterial, column 'Volume women'
    v_venous_male <- fixed(3733)
    label("Venous blood volume, men (mL)")  # Table 1, row Blood -Venous, column 'Volume men'
    v_venous_female <- fixed(2733)
    label("Venous blood volume, women (mL)")  # Table 1, row Blood -Venous, column 'Volume women'

    q_adipose_male <- fixed(325)
    label("Adipose blood flow, men (mL/min)")  # Table 1, row Adipose, column 'Blood flow men'
    q_adipose_female <- fixed(501.5)
    label("Adipose blood flow, women (mL/min)")  # Table 1, row Adipose, column 'Blood flow women'
    q_bone_male <- fixed(325)
    label("Bone blood flow, men (mL/min)")  # Table 1, row Bone, column 'Blood flow men'
    q_bone_female <- fixed(295)
    label("Bone blood flow, women (mL/min)")  # Table 1, row Bone, column 'Blood flow women'
    q_brain_male <- fixed(780)
    label("Brain blood flow, men (mL/min)")  # Table 1, row Brain, column 'Blood flow men'
    q_brain_female <- fixed(708)
    label("Brain blood flow, women (mL/min)")  # Table 1, row Brain, column 'Blood flow women'
    q_gut_male <- fixed(975)
    label("Gut blood flow, men (mL/min)")  # Table 1, row Gut, column 'Blood flow men'
    q_gut_female <- fixed(1003)
    label("Gut blood flow, women (mL/min)")  # Table 1, row Gut, column 'Blood flow women'
    q_heart_male <- fixed(260)
    label("Heart blood flow, men (mL/min)")  # Table 1, row Heart, column 'Blood flow men'
    q_heart_female <- fixed(295)
    label("Heart blood flow, women (mL/min)")  # Table 1, row Heart, column 'Blood flow women'
    q_kidney_male <- fixed(1235)
    label("Kidney blood flow, men (mL/min)")  # Table 1, row Kidney, column 'Blood flow men'
    q_kidney_female <- fixed(1003)
    label("Kidney blood flow, women (mL/min)")  # Table 1, row Kidney, column 'Blood flow women'
    q_liver_male <- fixed(1657.5)
    label("Liver blood flow, men (mL/min)")  # Table 1, row Liver, column 'Blood flow men'
    q_liver_female <- fixed(1593)
    label("Liver blood flow, women (mL/min)")  # Table 1, row Liver, column 'Blood flow women'
    q_lung_male <- fixed(6500)
    label("Lung blood flow, men (mL/min)")  # Table 1, row Lung, column 'Blood flow men'
    q_lung_female <- fixed(5900)
    label("Lung blood flow, women (mL/min)")  # Table 1, row Lung, column 'Blood flow women'
    q_muscle_male <- fixed(1105)
    label("Muscle blood flow, men (mL/min)")  # Table 1, row Muscle, column 'Blood flow men'
    q_muscle_female <- fixed(708)
    label("Muscle blood flow, women (mL/min)")  # Table 1, row Muscle, column 'Blood flow women'
    q_other_male <- fixed(487.5)
    label("Other blood flow, men (mL/min)")  # Table 1, row Other, column 'Blood flow men'
    q_other_female <- fixed(501.5)
    label("Other blood flow, women (mL/min)")  # Table 1, row Other, column 'Blood flow women'
    q_skin_male <- fixed(325)
    label("Skin blood flow, men (mL/min)")  # Table 1, row Skin, column 'Blood flow men'
    q_skin_female <- fixed(295)
    label("Skin blood flow, women (mL/min)")  # Table 1, row Skin, column 'Blood flow women'
    q_spleen_male <- fixed(195)
    label("Spleen blood flow, men (mL/min)")  # Table 1, row Spleen, column 'Blood flow men'
    q_spleen_female <- fixed(177)
    label("Spleen blood flow, women (mL/min)")  # Table 1, row Spleen, column 'Blood flow women'

    bodyfat_ref_male <- fixed(24.932)
    label("Reference percent body fat, men (percent)")  # derived: Table 1 adipose 18,200 mL / 73,000 mL total
    bodyfat_ref_female <- fixed(37.5)
    label("Reference percent body fat, women (percent)")  # derived: Table 1 adipose 22,500 mL / 60,000 mL total

    # ==================================================================
    # Tissue:plasma partition coefficients -- Supplementary Table 2
    # (PSP4-5-74-s003.docx). Supplement section 1 note 2: calculated with
    # the Poulin/Theil and Berezhkovskiy formulas. Note 4: the caffeine
    # and theophylline coefficients were then uniformly multiplied by 0.6
    # and 0.8 respectively so predicted Cmax matched observed data -- the
    # tabulated values are POST-adjustment (the 'Rest of body' entries are
    # exactly 0.8 for THP and 0.6 for CAF, i.e. an unadjusted value of 1).
    # ==================================================================
    lkp_adipose_thp <- fixed(log(0.088))
    label("Log adipose:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Adipose, column P(alpha:plasma) THP = 0.088
    lkp_bone_thp <- fixed(log(0.33))
    label("Log bone:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Bone, column P(alpha:plasma) THP = 0.33
    lkp_brain_thp <- fixed(log(0.59))
    label("Log brain:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Brain, column P(alpha:plasma) THP = 0.59
    lkp_gut_thp <- fixed(log(0.52))
    label("Log gut:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Gut, column P(alpha:plasma) THP = 0.52
    lkp_heart_thp <- fixed(log(0.52))
    label("Log heart:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Heart, column P(alpha:plasma) THP = 0.52
    lkp_kidney_thp <- fixed(log(0.55))
    label("Log kidney:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Kidney, column P(alpha:plasma) THP = 0.55
    lkp_liver_thp <- fixed(log(0.53))
    label("Log liver:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Liver, column P(alpha:plasma) THP = 0.53
    lkp_lung_thp <- fixed(log(0.55))
    label("Log lung:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Lung, column P(alpha:plasma) THP = 0.55
    lkp_muscle_thp <- fixed(log(0.53))
    label("Log muscle:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Muscle, column P(alpha:plasma) THP = 0.53
    lkp_other_thp <- fixed(log(0.8))
    label("Log other:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Other, column P(alpha:plasma) THP = 0.8
    lkp_skin_thp <- fixed(log(0.51))
    label("Log skin:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Skin, column P(alpha:plasma) THP = 0.51
    lkp_spleen_thp <- fixed(log(0.55))
    label("Log spleen:plasma partition coefficient, theophylline (unitless)")  # Table 2, row Spleen, column P(alpha:plasma) THP = 0.55

    lkp_adipose_caf <- fixed(log(0.13))
    label("Log adipose:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Adipose, column P(alpha:plasma) CAF = 0.13
    lkp_bone_caf <- fixed(log(0.29))
    label("Log bone:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Bone, column P(alpha:plasma) CAF = 0.29
    lkp_brain_caf <- fixed(log(0.49))
    label("Log brain:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Brain, column P(alpha:plasma) CAF = 0.49
    lkp_gut_caf <- fixed(log(0.44))
    label("Log gut:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Gut, column P(alpha:plasma) CAF = 0.44
    lkp_heart_caf <- fixed(log(0.44))
    label("Log heart:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Heart, column P(alpha:plasma) CAF = 0.44
    lkp_kidney_caf <- fixed(log(0.46))
    label("Log kidney:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Kidney, column P(alpha:plasma) CAF = 0.46
    lkp_liver_caf <- fixed(log(0.46))
    label("Log liver:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Liver, column P(alpha:plasma) CAF = 0.46
    lkp_lung_caf <- fixed(log(0.47))
    label("Log lung:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Lung, column P(alpha:plasma) CAF = 0.47
    lkp_muscle_caf <- fixed(log(0.45))
    label("Log muscle:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Muscle, column P(alpha:plasma) CAF = 0.45
    lkp_other_caf <- fixed(log(0.6))
    label("Log other:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Other, column P(alpha:plasma) CAF = 0.6
    lkp_skin_caf <- fixed(log(0.43))
    label("Log skin:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Skin, column P(alpha:plasma) CAF = 0.43
    lkp_spleen_caf <- fixed(log(0.47))
    label("Log spleen:plasma partition coefficient, caffeine (unitless)")  # Table 2, row Spleen, column P(alpha:plasma) CAF = 0.47

    lkp_adipose_px <- fixed(log(0.22))
    label("Log adipose:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Adipose, column P(alpha:plasma) PX = 0.22
    lkp_bone_px <- fixed(log(0.47))
    label("Log bone:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Bone, column P(alpha:plasma) PX = 0.47
    lkp_brain_px <- fixed(log(0.83))
    label("Log brain:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Brain, column P(alpha:plasma) PX = 0.83
    lkp_gut_px <- fixed(log(0.74))
    label("Log gut:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Gut, column P(alpha:plasma) PX = 0.74
    lkp_heart_px <- fixed(log(0.76))
    label("Log heart:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Heart, column P(alpha:plasma) PX = 0.76
    lkp_kidney_px <- fixed(log(0.79))
    label("Log kidney:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Kidney, column P(alpha:plasma) PX = 0.79
    lkp_liver_px <- fixed(log(0.77))
    label("Log liver:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Liver, column P(alpha:plasma) PX = 0.77
    lkp_lung_px <- fixed(log(0.8))
    label("Log lung:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Lung, column P(alpha:plasma) PX = 0.8
    lkp_muscle_px <- fixed(log(0.76))
    label("Log muscle:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Muscle, column P(alpha:plasma) PX = 0.76
    lkp_other_px <- fixed(log(1))
    label("Log other:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Other, column P(alpha:plasma) PX = 1
    lkp_skin_px <- fixed(log(0.73))
    label("Log skin:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Skin, column P(alpha:plasma) PX = 0.73
    lkp_spleen_px <- fixed(log(0.8))
    label("Log spleen:plasma partition coefficient, paraxanthine (unitless)")  # Table 2, row Spleen, column P(alpha:plasma) PX = 0.8

    lkp_adipose_cip <- fixed(log(0.25))
    label("Log adipose:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Adipose, column P(alpha:plasma) CIP = 0.25
    lkp_bone_cip <- fixed(log(0.36))
    label("Log bone:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Bone, column P(alpha:plasma) CIP = 0.36
    lkp_brain_cip <- fixed(log(0.64))
    label("Log brain:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Brain, column P(alpha:plasma) CIP = 0.64
    lkp_gut_cip <- fixed(log(0.59))
    label("Log gut:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Gut, column P(alpha:plasma) CIP = 0.59
    lkp_heart_cip <- fixed(log(0.62))
    label("Log heart:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Heart, column P(alpha:plasma) CIP = 0.62
    lkp_kidney_cip <- fixed(log(0.64))
    label("Log kidney:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Kidney, column P(alpha:plasma) CIP = 0.64
    lkp_liver_cip <- fixed(log(0.62))
    label("Log liver:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Liver, column P(alpha:plasma) CIP = 0.62
    lkp_lung_cip <- fixed(log(0.66))
    label("Log lung:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Lung, column P(alpha:plasma) CIP = 0.66
    lkp_muscle_cip <- fixed(log(0.61))
    label("Log muscle:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Muscle, column P(alpha:plasma) CIP = 0.61
    lkp_other_cip <- fixed(log(1))
    label("Log other:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Other, column P(alpha:plasma) CIP = 1
    lkp_skin_cip <- fixed(log(0.58))
    label("Log skin:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Skin, column P(alpha:plasma) CIP = 0.58
    lkp_spleen_cip <- fixed(log(0.64))
    label("Log spleen:plasma partition coefficient, ciprofloxacin (unitless)")  # Table 2, row Spleen, column P(alpha:plasma) CIP = 0.64

    # Blood-to-plasma ratio BP. Eq 1 and Eq 4 use the partition coefficient
    # ONLY through the composite P(alpha:plasma)/BP, and BP is not reported
    # anywhere in the paper or the supplement. It is fixed at 1 so the
    # composite equals the tabulated coefficient; the 'Rest of body'
    # entries of exactly 1.00 for CIP and PX are what an unadjusted
    # coefficient with BP = 1 gives. See the vignette Errata.
    bp_thp <- fixed(1)
    label("Blood-to-plasma ratio, theophylline (unitless)")  # not reported; see Errata
    bp_caf <- fixed(1)
    label("Blood-to-plasma ratio, caffeine (unitless)")  # not reported; see Errata
    bp_px <- fixed(1)
    label("Blood-to-plasma ratio, paraxanthine (unitless)")  # not reported; see Errata
    bp_cip <- fixed(1)
    label("Blood-to-plasma ratio, ciprofloxacin (unitless)")  # not reported; see Errata

    # ==================================================================
    # Renal extraction ratios and absorption -- Supplementary Table 3
    # (PSP4-5-74-s004.docx). Eq 3: dZ/dt = E * Q * Cab / V, so renal
    # clearance is E * Q(kidney); for THP that is 0.004 * 1235 = 4.94
    # mL/min = 0.30 L/h, about a tenth of total theophylline clearance.
    # ==================================================================
    eren_thp <- fixed(0.004)
    label("Renal extraction ratio, theophylline (unitless)")  # Table 3, row E(kidney) THP
    eren_caf <- fixed(0.008)
    label("Renal extraction ratio, caffeine (unitless)")  # Table 3, row E(kidney) CAF
    eren_px <- fixed(0.008)
    label("Renal extraction ratio, paraxanthine (unitless)")  # Table 3, row E(kidney) PX
    eren_cip <- fixed(0.39)
    label("Renal extraction ratio, ciprofloxacin (unitless)")  # Table 3, row E(kidney) CIP
    fabs_cip <- fixed(0.7)
    label("Fraction of dose absorbed, ciprofloxacin (unitless)")  # Table 3, row F(CIP)
    tge <- fixed(30)
    label("Gastric emptying time (min)")  # Table 3, rows Tge(CAF), Tge(THP), Tge(CIP), all 30; supplement section 1 'Tge for all four drugs were set to 30 minutes'
    tsi <- fixed(199)
    label("Small-intestinal transit time (min)")  # Table 3, rows Tsi(CAF), Tsi(THP), Tsi(CIP), all 199; supplement section 1
    nsi <- fixed(7)
    label("Number of small-intestine transit compartments (count)")  # supplement section 1, N = 7 per Yu and Amidon
    ka_limit_mult <- fixed(100)
    label("Absorption rate constant as a multiple of Kt for a fully absorbed drug (unitless)")  # numerical stand-in for the F = 1 singularity of Ka = Kt*((1-F)^(-1/N)-1); see Errata

    # ==================================================================
    # Hepatic Michaelis-Menten parameters -- Supplementary Table 3.
    # Vmax in M/min, Km in mmol/mL (= M). Every reaction of Table 3 is
    # carried; products other than PX and THP are terminal sinks. The
    # CIP-CYP3A4 Km has no partner Vmax because CIP is not a CYP3A4
    # substrate here -- it enters only as a competitive inhibitor.
    # ==================================================================
    vmax_thp_1mx_1a2 <- fixed(3.4e-07)
    label("Vmax, THP -> 1MX by CYP1A2 (M/min)")  # Table 3, row Vmax THP -> 1MX by CYP1A2
    km_thp_1mx_1a2 <- fixed(2.8e-05)
    label("Km, THP -> 1MX by CYP1A2 (mmol/mL)")  # Table 3, row Km THP -> 1MX by CYP1A2
    vmax_thp_3mx_1a2 <- fixed(1.6e-07)
    label("Vmax, THP -> 3MX by CYP1A2 (M/min)")  # Table 3, row Vmax THP -> 3MX by CYP1A2
    km_thp_3mx_1a2 <- fixed(1.3e-05)
    label("Km, THP -> 3MX by CYP1A2 (mmol/mL)")  # Table 3, row Km THP -> 3MX by CYP1A2
    vmax_thp_13u_1a2 <- fixed(5e-06)
    label("Vmax, THP -> 13U by CYP1A2 (M/min)")  # Table 3, row Vmax THP -> 13U by CYP1A2
    km_thp_13u_1a2 <- fixed(0.00062)
    label("Km, THP -> 13U by CYP1A2 (mmol/mL)")  # Table 3, row Km THP -> 13U by CYP1A2
    vmax_thp_13u_2e1 <- fixed(6e-05)
    label("Vmax, THP -> 13U by CYP2E1 (M/min)")  # Table 3, row Vmax THP -> 13U by CYP2E1
    km_thp_13u_2e1 <- fixed(0.015)
    label("Km, THP -> 13U by CYP2E1 (mmol/mL)")  # Table 3, row Km THP -> 13U by CYP2E1
    vmax_caf_px_1a2 <- fixed(8.4e-06)
    label("Vmax, CAF -> PX by CYP1A2 (M/min)")  # Table 3, row Vmax CAF -> PX by CYP1A2
    km_caf_px_1a2 <- fixed(0.00019)
    label("Km, CAF -> PX by CYP1A2 (mmol/mL)")  # Table 3, row Km CAF -> PX by CYP1A2
    vmax_caf_tb_1a2 <- fixed(8.3e-07)
    label("Vmax, CAF -> TB by CYP1A2 (M/min)")  # Table 3, row Vmax CAF -> TB by CYP1A2
    km_caf_tb_1a2 <- fixed(0.00016)
    label("Km, CAF -> TB by CYP1A2 (mmol/mL)")  # Table 3, row Km CAF -> TB by CYP1A2
    vmax_caf_tb_2e1 <- fixed(7e-08)
    label("Vmax, CAF -> TB by CYP2E1 (M/min)")  # Table 3, row Vmax CAF -> TB by CYP2E1
    km_caf_tb_2e1 <- fixed(0.0014)
    label("Km, CAF -> TB by CYP2E1 (mmol/mL)")  # Table 3, row Km CAF -> TB by CYP2E1
    vmax_caf_thp_1a2 <- fixed(3.1e-07)
    label("Vmax, CAF -> THP by CYP1A2 (M/min)")  # Table 3, row Vmax CAF -> THP by CYP1A2
    km_caf_thp_1a2 <- fixed(0.00025)
    label("Km, CAF -> THP by CYP1A2 (mmol/mL)")  # Table 3, row Km CAF -> THP by CYP1A2
    vmax_caf_thp_2e1 <- fixed(5.2e-08)
    label("Vmax, CAF -> THP by CYP2E1 (M/min)")  # Table 3, row Vmax CAF -> THP by CYP2E1
    km_caf_thp_2e1 <- fixed(0.00084)
    label("Km, CAF -> THP by CYP2E1 (mmol/mL)")  # Table 3, row Km CAF -> THP by CYP2E1
    vmax_caf_ta_1a2 <- fixed(5.4e-07)
    label("Vmax, CAF -> TA by CYP1A2 (M/min)")  # Table 3, row Vmax CAF -> TA by CYP1A2
    km_caf_ta_1a2 <- fixed(0.00027)
    label("Km, CAF -> TA by CYP1A2 (mmol/mL)")  # Table 3, row Km CAF -> TA by CYP1A2
    vmax_caf_ta_2e1 <- fixed(4.3e-07)
    label("Vmax, CAF -> TA by CYP2E1 (M/min)")  # Table 3, row Vmax CAF -> TA by CYP2E1
    km_caf_ta_2e1 <- fixed(0.001)
    label("Km, CAF -> TA by CYP2E1 (mmol/mL)")  # Table 3, row Km CAF -> TA by CYP2E1
    vmax_cip_1a2 <- fixed(8e-07)
    label("Vmax, CIP -> CIP by CYP1A2 (M/min)")  # Table 3, row Vmax CIP -> CIP by CYP1A2
    km_cip_1a2 <- fixed(4.6e-06)
    label("Km, CIP -> CIP by CYP1A2 (mmol/mL)")  # Table 3, row Km CIP -> CIP by CYP1A2
    vmax_px_1a2 <- fixed(1.3e-05)
    label("Vmax, PX -> PX by CYP1A2 (M/min)")  # Table 3, row Vmax PX -> PX by CYP1A2
    km_px_1a2 <- fixed(0.00017)
    label("Km, PX -> PX by CYP1A2 (mmol/mL)")  # Table 3, row Km PX -> PX by CYP1A2
    vmax_caf_ta_3a4 <- fixed(1.7e-05)
    label("Vmax, CAF -> TA by CYP3A4 (M/min)")  # Table 3, row Vmax CAF -> TA by CYP3A4
    km_caf_ta_3a4 <- fixed(0.046)
    label("Km, CAF -> TA by CYP3A4 (mmol/mL)")  # Table 3, row Km CAF -> TA by CYP3A4
    km_cip_3a4_inhib <- fixed(5.1e-05)
    label("Inhibition constant of ciprofloxacin at CYP3A4 (mmol/mL)")  # Table 3, row Km CIP CYP3A4; no partner Vmax. Supplement note 9: McLellan 2 mM divided by 39

    mw_thp <- fixed(180.164)
    label("Molar mass, theophylline (g/mol)")  # standard chemical constant, used only for mg <-> mmol conversion
    mw_caf <- fixed(194.191)
    label("Molar mass, caffeine (g/mol)")  # standard chemical constant, used only for mg <-> mmol conversion
    mw_px <- fixed(180.164)
    label("Molar mass, paraxanthine (g/mol)")  # standard chemical constant, used only for mg <-> mmol conversion
    mw_cip <- fixed(331.346)
    label("Molar mass, ciprofloxacin (g/mol)")  # standard chemical constant, used only for mg <-> mmol conversion
  })

  model({
    # ---- sex-specific physiology (Supplementary Table 1) -------------
    v_adipose <- v_adipose_male * (1 - SEXF) + v_adipose_female * SEXF
    q_adipose <- q_adipose_male * (1 - SEXF) + q_adipose_female * SEXF
    v_bone <- v_bone_male * (1 - SEXF) + v_bone_female * SEXF
    q_bone <- q_bone_male * (1 - SEXF) + q_bone_female * SEXF
    v_brain <- v_brain_male * (1 - SEXF) + v_brain_female * SEXF
    q_brain <- q_brain_male * (1 - SEXF) + q_brain_female * SEXF
    v_gut <- v_gut_male * (1 - SEXF) + v_gut_female * SEXF
    q_gut <- q_gut_male * (1 - SEXF) + q_gut_female * SEXF
    v_heart <- v_heart_male * (1 - SEXF) + v_heart_female * SEXF
    q_heart <- q_heart_male * (1 - SEXF) + q_heart_female * SEXF
    v_kidney <- v_kidney_male * (1 - SEXF) + v_kidney_female * SEXF
    q_kidney <- q_kidney_male * (1 - SEXF) + q_kidney_female * SEXF
    v_liver <- v_liver_male * (1 - SEXF) + v_liver_female * SEXF
    q_liver <- q_liver_male * (1 - SEXF) + q_liver_female * SEXF
    v_lung <- v_lung_male * (1 - SEXF) + v_lung_female * SEXF
    q_lung <- q_lung_male * (1 - SEXF) + q_lung_female * SEXF
    v_muscle <- v_muscle_male * (1 - SEXF) + v_muscle_female * SEXF
    q_muscle <- q_muscle_male * (1 - SEXF) + q_muscle_female * SEXF
    v_other <- v_other_male * (1 - SEXF) + v_other_female * SEXF
    q_other <- q_other_male * (1 - SEXF) + q_other_female * SEXF
    v_skin <- v_skin_male * (1 - SEXF) + v_skin_female * SEXF
    q_skin <- q_skin_male * (1 - SEXF) + q_skin_female * SEXF
    v_spleen <- v_spleen_male * (1 - SEXF) + v_spleen_female * SEXF
    q_spleen <- q_spleen_male * (1 - SEXF) + q_spleen_female * SEXF
    v_arterial <- v_arterial_male * (1 - SEXF) + v_arterial_female * SEXF
    v_venous <- v_venous_male * (1 - SEXF) + v_venous_female * SEXF

    # Supplement section 5: race differences are simulated by scaling the
    # adipose volume by Phi = (mean adipose body fraction of the group) /
    # (mean adipose body fraction of the reference group). Phi is supplied
    # here as BODYFAT_PCT relative to the reference-group percent body fat,
    # because the per-race fractions live in Carpenter 2013, which was not
    # available when this model was built. BODYFAT_PCT at its reference
    # value leaves the model at the Table 1 adipose volume.
    bodyfat_ref <- bodyfat_ref_male * (1 - SEXF) + bodyfat_ref_female * SEXF
    v_adipose <- v_adipose * BODYFAT_PCT / bodyfat_ref

    # Cardiac output is the lung flow; the hepatic artery is the liver
    # flow less the portal (gut + spleen) inflow.
    q_co <- q_lung
    q_ha <- q_liver - q_gut - q_spleen

    # ---- absorption rate constants (supplement section 1) ------------
    kge <- 1 / tge
    kt <- nsi / tsi
    # Ka = Kt * (1/((1-F)^(1/N)) - 1). F(CIP) = 0.7 evaluates directly.
    ka_cip <- kt * ((1 - fabs_cip)^(-1 / nsi) - 1)
    # Table 3 gives F(THP) = F(CAF) = 1 exactly, for which that expression
    # is singular: complete absorption requires Ka -> Inf, and the limit is
    # gastric-emptying-limited absorption. ka_limit_mult is the numerical
    # stand-in for that limit, not a fitted rate; the vignette shows Cmax
    # is unchanged to four significant figures across a tenfold range.
    ka_thp <- kt * ka_limit_mult
    ka_caf <- kt * ka_limit_mult

    # ---- partition coefficients back-transformed from log scale ------
    kp_adipose_thp <- exp(lkp_adipose_thp)
    kp_bone_thp <- exp(lkp_bone_thp)
    kp_brain_thp <- exp(lkp_brain_thp)
    kp_gut_thp <- exp(lkp_gut_thp)
    kp_heart_thp <- exp(lkp_heart_thp)
    kp_kidney_thp <- exp(lkp_kidney_thp)
    kp_liver_thp <- exp(lkp_liver_thp)
    kp_lung_thp <- exp(lkp_lung_thp)
    kp_muscle_thp <- exp(lkp_muscle_thp)
    kp_other_thp <- exp(lkp_other_thp)
    kp_skin_thp <- exp(lkp_skin_thp)
    kp_spleen_thp <- exp(lkp_spleen_thp)

    kp_adipose_caf <- exp(lkp_adipose_caf)
    kp_bone_caf <- exp(lkp_bone_caf)
    kp_brain_caf <- exp(lkp_brain_caf)
    kp_gut_caf <- exp(lkp_gut_caf)
    kp_heart_caf <- exp(lkp_heart_caf)
    kp_kidney_caf <- exp(lkp_kidney_caf)
    kp_liver_caf <- exp(lkp_liver_caf)
    kp_lung_caf <- exp(lkp_lung_caf)
    kp_muscle_caf <- exp(lkp_muscle_caf)
    kp_other_caf <- exp(lkp_other_caf)
    kp_skin_caf <- exp(lkp_skin_caf)
    kp_spleen_caf <- exp(lkp_spleen_caf)

    kp_adipose_px <- exp(lkp_adipose_px)
    kp_bone_px <- exp(lkp_bone_px)
    kp_brain_px <- exp(lkp_brain_px)
    kp_gut_px <- exp(lkp_gut_px)
    kp_heart_px <- exp(lkp_heart_px)
    kp_kidney_px <- exp(lkp_kidney_px)
    kp_liver_px <- exp(lkp_liver_px)
    kp_lung_px <- exp(lkp_lung_px)
    kp_muscle_px <- exp(lkp_muscle_px)
    kp_other_px <- exp(lkp_other_px)
    kp_skin_px <- exp(lkp_skin_px)
    kp_spleen_px <- exp(lkp_spleen_px)

    kp_adipose_cip <- exp(lkp_adipose_cip)
    kp_bone_cip <- exp(lkp_bone_cip)
    kp_brain_cip <- exp(lkp_brain_cip)
    kp_gut_cip <- exp(lkp_gut_cip)
    kp_heart_cip <- exp(lkp_heart_cip)
    kp_kidney_cip <- exp(lkp_kidney_cip)
    kp_liver_cip <- exp(lkp_liver_cip)
    kp_lung_cip <- exp(lkp_lung_cip)
    kp_muscle_cip <- exp(lkp_muscle_cip)
    kp_other_cip <- exp(lkp_other_cip)
    kp_skin_cip <- exp(lkp_skin_cip)
    kp_spleen_cip <- exp(lkp_spleen_cip)

    # ---- emergent (venous-equilibrium) concentrations, Eq 1 ----------
    # Eq 1 divides the tissue concentration by P(alpha:plasma)/BP.
    cv_adipose_thp <- adipose_thp * bp_thp / kp_adipose_thp
    cv_bone_thp <- bone_thp * bp_thp / kp_bone_thp
    cv_brain_thp <- brain_thp * bp_thp / kp_brain_thp
    cv_gut_thp <- gut_thp * bp_thp / kp_gut_thp
    cv_heart_thp <- heart_thp * bp_thp / kp_heart_thp
    cv_kidney_thp <- kidney_thp * bp_thp / kp_kidney_thp
    cv_liver_thp <- liver_thp * bp_thp / kp_liver_thp
    cv_lung_thp <- lung_thp * bp_thp / kp_lung_thp
    cv_muscle_thp <- muscle_thp * bp_thp / kp_muscle_thp
    cv_other_thp <- other_thp * bp_thp / kp_other_thp
    cv_skin_thp <- skin_thp * bp_thp / kp_skin_thp
    cv_spleen_thp <- spleen_thp * bp_thp / kp_spleen_thp

    cv_adipose_caf <- adipose_caf * bp_caf / kp_adipose_caf
    cv_bone_caf <- bone_caf * bp_caf / kp_bone_caf
    cv_brain_caf <- brain_caf * bp_caf / kp_brain_caf
    cv_gut_caf <- gut_caf * bp_caf / kp_gut_caf
    cv_heart_caf <- heart_caf * bp_caf / kp_heart_caf
    cv_kidney_caf <- kidney_caf * bp_caf / kp_kidney_caf
    cv_liver_caf <- liver_caf * bp_caf / kp_liver_caf
    cv_lung_caf <- lung_caf * bp_caf / kp_lung_caf
    cv_muscle_caf <- muscle_caf * bp_caf / kp_muscle_caf
    cv_other_caf <- other_caf * bp_caf / kp_other_caf
    cv_skin_caf <- skin_caf * bp_caf / kp_skin_caf
    cv_spleen_caf <- spleen_caf * bp_caf / kp_spleen_caf

    cv_adipose_px <- adipose_px * bp_px / kp_adipose_px
    cv_bone_px <- bone_px * bp_px / kp_bone_px
    cv_brain_px <- brain_px * bp_px / kp_brain_px
    cv_gut_px <- gut_px * bp_px / kp_gut_px
    cv_heart_px <- heart_px * bp_px / kp_heart_px
    cv_kidney_px <- kidney_px * bp_px / kp_kidney_px
    cv_liver_px <- liver_px * bp_px / kp_liver_px
    cv_lung_px <- lung_px * bp_px / kp_lung_px
    cv_muscle_px <- muscle_px * bp_px / kp_muscle_px
    cv_other_px <- other_px * bp_px / kp_other_px
    cv_skin_px <- skin_px * bp_px / kp_skin_px
    cv_spleen_px <- spleen_px * bp_px / kp_spleen_px

    cv_adipose_cip <- adipose_cip * bp_cip / kp_adipose_cip
    cv_bone_cip <- bone_cip * bp_cip / kp_bone_cip
    cv_brain_cip <- brain_cip * bp_cip / kp_brain_cip
    cv_gut_cip <- gut_cip * bp_cip / kp_gut_cip
    cv_heart_cip <- heart_cip * bp_cip / kp_heart_cip
    cv_kidney_cip <- kidney_cip * bp_cip / kp_kidney_cip
    cv_liver_cip <- liver_cip * bp_cip / kp_liver_cip
    cv_lung_cip <- lung_cip * bp_cip / kp_lung_cip
    cv_muscle_cip <- muscle_cip * bp_cip / kp_muscle_cip
    cv_other_cip <- other_cip * bp_cip / kp_other_cip
    cv_skin_cip <- skin_cip * bp_cip / kp_skin_cip
    cv_spleen_cip <- spleen_cip * bp_cip / kp_spleen_cip

    # ---- hepatic concentrations driving metabolism, Eq 4 -------------
    u_thp <- cv_liver_thp
    u_caf <- cv_liver_caf
    u_px <- cv_liver_px
    u_cip <- cv_liver_cip

    # ---- hepatic reaction velocities, Eq 4 ---------------------------
    # Competitive inhibition: Km(i,j) * (1 + sum_l U(l)/Ki(l,j)) + U(i).
    # Methods: 'we assume that the enzyme-inhibitor complex is in
    # equilibrium with the free enzyme and the inhibitor and we use Km
    # values instead of Ki', so Ki(l,j) = Km(l,j) throughout. The
    # inhibitor index l of Eq 4 runs over competing COMPOUNDS, not over
    # reactions: each other compound handled by the same enzyme
    # contributes ONE term, using the Km of its dominant (highest Vmax)
    # pathway on that enzyme as the Ki. A compound's own parallel
    # pathways are not inhibitors of one another, so with a single
    # compound present each pathway is plain Michaelis-Menten. Two rival
    # readings -- summing over every other REACTION, with or without the
    # substrate's own -- were implemented and are falsified by the
    # paper's own Table 1; see the vignette Errata.
    # HEPFUNC_REL is Delta(met) of Eq 4.
    v_thp_1mx_1a2 <- HEPFUNC_REL * vmax_thp_1mx_1a2 * u_thp /
      (km_thp_1mx_1a2 * (1 + u_caf / km_caf_px_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_thp)
    v_thp_3mx_1a2 <- HEPFUNC_REL * vmax_thp_3mx_1a2 * u_thp /
      (km_thp_3mx_1a2 * (1 + u_caf / km_caf_px_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_thp)
    v_thp_13u_1a2 <- HEPFUNC_REL * vmax_thp_13u_1a2 * u_thp /
      (km_thp_13u_1a2 * (1 + u_caf / km_caf_px_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_thp)
    v_thp_13u_2e1 <- HEPFUNC_REL * vmax_thp_13u_2e1 * u_thp /
      (km_thp_13u_2e1 * (1 + u_caf / km_caf_ta_2e1) + u_thp)
    v_caf_px_1a2 <- HEPFUNC_REL * vmax_caf_px_1a2 * u_caf /
      (km_caf_px_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_caf)
    v_caf_tb_1a2 <- HEPFUNC_REL * vmax_caf_tb_1a2 * u_caf /
      (km_caf_tb_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_caf)
    v_caf_tb_2e1 <- HEPFUNC_REL * vmax_caf_tb_2e1 * u_caf /
      (km_caf_tb_2e1 * (1 + u_thp / km_thp_13u_2e1) + u_caf)
    v_caf_thp_1a2 <- HEPFUNC_REL * vmax_caf_thp_1a2 * u_caf /
      (km_caf_thp_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_caf)
    v_caf_thp_2e1 <- HEPFUNC_REL * vmax_caf_thp_2e1 * u_caf /
      (km_caf_thp_2e1 * (1 + u_thp / km_thp_13u_2e1) + u_caf)
    v_caf_ta_1a2 <- HEPFUNC_REL * vmax_caf_ta_1a2 * u_caf /
      (km_caf_ta_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_cip / km_cip_1a2 + u_px / km_px_1a2) + u_caf)
    v_caf_ta_2e1 <- HEPFUNC_REL * vmax_caf_ta_2e1 * u_caf /
      (km_caf_ta_2e1 * (1 + u_thp / km_thp_13u_2e1) + u_caf)
    v_cip_1a2 <- HEPFUNC_REL * vmax_cip_1a2 * u_cip /
      (km_cip_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_caf / km_caf_px_1a2 + u_px / km_px_1a2) + u_cip)
    v_px_1a2 <- HEPFUNC_REL * vmax_px_1a2 * u_px /
      (km_px_1a2 * (1 + u_thp / km_thp_13u_1a2 + u_caf / km_caf_px_1a2 + u_cip / km_cip_1a2) + u_px)
    v_caf_ta_3a4 <- HEPFUNC_REL * vmax_caf_ta_3a4 * u_caf /
      (km_caf_ta_3a4 * (1 + u_cip / km_cip_3a4_inhib) + u_caf)

    # ================= THEOPHYLLINE =================
    # Compartmental absorption and transit model of Yu and Amidon
    # (supplement section 1), states carried as amounts in mmol.
    d/dt(stomach_thp) <- -kge * stomach_thp
    d/dt(transit1_thp) <- kge * stomach_thp - kt * transit1_thp - ka_thp * transit1_thp
    d/dt(transit2_thp) <- kt * transit1_thp - kt * transit2_thp - ka_thp * transit2_thp
    d/dt(transit3_thp) <- kt * transit2_thp - kt * transit3_thp - ka_thp * transit3_thp
    d/dt(transit4_thp) <- kt * transit3_thp - kt * transit4_thp - ka_thp * transit4_thp
    d/dt(transit5_thp) <- kt * transit4_thp - kt * transit5_thp - ka_thp * transit5_thp
    d/dt(transit6_thp) <- kt * transit5_thp - kt * transit6_thp - ka_thp * transit6_thp
    d/dt(transit7_thp) <- kt * transit6_thp - kt * transit7_thp - ka_thp * transit7_thp
    d/dt(colon_thp) <- kt * transit7_thp
    d/dt(adipose_thp) <- q_adipose * (arterial_thp - cv_adipose_thp) / v_adipose
    d/dt(bone_thp) <- q_bone * (arterial_thp - cv_bone_thp) / v_bone
    d/dt(brain_thp) <- q_brain * (arterial_thp - cv_brain_thp) / v_brain
    d/dt(gut_thp) <- q_gut * (arterial_thp - cv_gut_thp) / v_gut+
      (ka_thp *
      (transit1_thp + transit2_thp + transit3_thp + transit4_thp + transit5_thp + transit6_thp + transit7_thp)) / v_gut
    d/dt(heart_thp) <- q_heart * (arterial_thp - cv_heart_thp) / v_heart
    d/dt(kidney_thp) <- q_kidney * (arterial_thp - cv_kidney_thp) / v_kidney -
      RENALFUNC_REL * eren_thp * q_kidney * arterial_thp / v_kidney
    d/dt(liver_thp) <- (q_ha * arterial_thp + q_gut * cv_gut_thp +
      q_spleen * cv_spleen_thp - q_liver * cv_liver_thp) / v_liver -
      (v_thp_1mx_1a2 + v_thp_3mx_1a2 + v_thp_13u_1a2 + v_thp_13u_2e1) +
      (v_caf_thp_1a2 + v_caf_thp_2e1)
    d/dt(lung_thp) <- q_co * (venous_thp - cv_lung_thp) / v_lung
    d/dt(muscle_thp) <- q_muscle * (arterial_thp - cv_muscle_thp) / v_muscle
    d/dt(other_thp) <- q_other * (arterial_thp - cv_other_thp) / v_other
    d/dt(skin_thp) <- q_skin * (arterial_thp - cv_skin_thp) / v_skin
    d/dt(spleen_thp) <- q_spleen * (arterial_thp - cv_spleen_thp) / v_spleen
    d/dt(arterial_thp) <- q_co * (cv_lung_thp - arterial_thp) / v_arterial
    d/dt(venous_thp) <-
      (q_adipose * cv_adipose_thp +
       q_bone * cv_bone_thp +
       q_brain * cv_brain_thp +
       q_heart * cv_heart_thp +
       q_kidney * cv_kidney_thp +
       q_liver * cv_liver_thp +
       q_muscle * cv_muscle_thp +
       q_other * cv_other_thp +
       q_skin * cv_skin_thp -
       q_co * venous_thp) / v_venous
    d/dt(urine_thp) <- RENALFUNC_REL * eren_thp * q_kidney * arterial_thp
    Cthp <- venous_thp * 1000 * mw_thp
    f(stomach_thp) <- 1 / mw_thp   # dose in mg -> mmol
    f(venous_thp) <- 1 / (mw_thp * v_venous)   # IV dose in mg -> mmol/mL

    # ================= CAFFEINE =================
    # Compartmental absorption and transit model of Yu and Amidon
    # (supplement section 1), states carried as amounts in mmol.
    d/dt(stomach_caf) <- -kge * stomach_caf
    d/dt(transit1_caf) <- kge * stomach_caf - kt * transit1_caf - ka_caf * transit1_caf
    d/dt(transit2_caf) <- kt * transit1_caf - kt * transit2_caf - ka_caf * transit2_caf
    d/dt(transit3_caf) <- kt * transit2_caf - kt * transit3_caf - ka_caf * transit3_caf
    d/dt(transit4_caf) <- kt * transit3_caf - kt * transit4_caf - ka_caf * transit4_caf
    d/dt(transit5_caf) <- kt * transit4_caf - kt * transit5_caf - ka_caf * transit5_caf
    d/dt(transit6_caf) <- kt * transit5_caf - kt * transit6_caf - ka_caf * transit6_caf
    d/dt(transit7_caf) <- kt * transit6_caf - kt * transit7_caf - ka_caf * transit7_caf
    d/dt(colon_caf) <- kt * transit7_caf
    d/dt(adipose_caf) <- q_adipose * (arterial_caf - cv_adipose_caf) / v_adipose
    d/dt(bone_caf) <- q_bone * (arterial_caf - cv_bone_caf) / v_bone
    d/dt(brain_caf) <- q_brain * (arterial_caf - cv_brain_caf) / v_brain
    d/dt(gut_caf) <- q_gut * (arterial_caf - cv_gut_caf) / v_gut+
      (ka_caf *
      (transit1_caf + transit2_caf + transit3_caf + transit4_caf + transit5_caf + transit6_caf + transit7_caf)) / v_gut
    d/dt(heart_caf) <- q_heart * (arterial_caf - cv_heart_caf) / v_heart
    d/dt(kidney_caf) <- q_kidney * (arterial_caf - cv_kidney_caf) / v_kidney -
      RENALFUNC_REL * eren_caf * q_kidney * arterial_caf / v_kidney
    d/dt(liver_caf) <- (q_ha * arterial_caf + q_gut * cv_gut_caf +
      q_spleen * cv_spleen_caf - q_liver * cv_liver_caf) / v_liver -
      (v_caf_px_1a2 + v_caf_tb_1a2 + v_caf_tb_2e1 + v_caf_thp_1a2 + v_caf_thp_2e1 + v_caf_ta_1a2 + v_caf_ta_2e1 + v_caf_ta_3a4)
    d/dt(lung_caf) <- q_co * (venous_caf - cv_lung_caf) / v_lung
    d/dt(muscle_caf) <- q_muscle * (arterial_caf - cv_muscle_caf) / v_muscle
    d/dt(other_caf) <- q_other * (arterial_caf - cv_other_caf) / v_other
    d/dt(skin_caf) <- q_skin * (arterial_caf - cv_skin_caf) / v_skin
    d/dt(spleen_caf) <- q_spleen * (arterial_caf - cv_spleen_caf) / v_spleen
    d/dt(arterial_caf) <- q_co * (cv_lung_caf - arterial_caf) / v_arterial
    d/dt(venous_caf) <-
      (q_adipose * cv_adipose_caf +
       q_bone * cv_bone_caf +
       q_brain * cv_brain_caf +
       q_heart * cv_heart_caf +
       q_kidney * cv_kidney_caf +
       q_liver * cv_liver_caf +
       q_muscle * cv_muscle_caf +
       q_other * cv_other_caf +
       q_skin * cv_skin_caf -
       q_co * venous_caf) / v_venous
    d/dt(urine_caf) <- RENALFUNC_REL * eren_caf * q_kidney * arterial_caf
    Ccaf <- venous_caf * 1000 * mw_caf
    f(stomach_caf) <- 1 / mw_caf   # dose in mg -> mmol
    f(venous_caf) <- 1 / (mw_caf * v_venous)   # IV dose in mg -> mmol/mL

    # ================= PARAXANTHINE =================
    # Paraxanthine is never dosed; it is formed only from caffeine.
    d/dt(adipose_px) <- q_adipose * (arterial_px - cv_adipose_px) / v_adipose
    d/dt(bone_px) <- q_bone * (arterial_px - cv_bone_px) / v_bone
    d/dt(brain_px) <- q_brain * (arterial_px - cv_brain_px) / v_brain
    d/dt(gut_px) <- q_gut * (arterial_px - cv_gut_px) / v_gut
    d/dt(heart_px) <- q_heart * (arterial_px - cv_heart_px) / v_heart
    d/dt(kidney_px) <- q_kidney * (arterial_px - cv_kidney_px) / v_kidney -
      RENALFUNC_REL * eren_px * q_kidney * arterial_px / v_kidney
    d/dt(liver_px) <- (q_ha * arterial_px + q_gut * cv_gut_px +
      q_spleen * cv_spleen_px - q_liver * cv_liver_px) / v_liver -
      (v_px_1a2) +
      (v_caf_px_1a2)
    d/dt(lung_px) <- q_co * (venous_px - cv_lung_px) / v_lung
    d/dt(muscle_px) <- q_muscle * (arterial_px - cv_muscle_px) / v_muscle
    d/dt(other_px) <- q_other * (arterial_px - cv_other_px) / v_other
    d/dt(skin_px) <- q_skin * (arterial_px - cv_skin_px) / v_skin
    d/dt(spleen_px) <- q_spleen * (arterial_px - cv_spleen_px) / v_spleen
    d/dt(arterial_px) <- q_co * (cv_lung_px - arterial_px) / v_arterial
    d/dt(venous_px) <-
      (q_adipose * cv_adipose_px +
       q_bone * cv_bone_px +
       q_brain * cv_brain_px +
       q_heart * cv_heart_px +
       q_kidney * cv_kidney_px +
       q_liver * cv_liver_px +
       q_muscle * cv_muscle_px +
       q_other * cv_other_px +
       q_skin * cv_skin_px -
       q_co * venous_px) / v_venous
    d/dt(urine_px) <- RENALFUNC_REL * eren_px * q_kidney * arterial_px
    Cpx <- venous_px * 1000 * mw_px
    f(venous_px) <- 1 / (mw_px * v_venous)   # IV dose in mg -> mmol/mL

    # ================= CIPROFLOXACIN =================
    # Compartmental absorption and transit model of Yu and Amidon
    # (supplement section 1), states carried as amounts in mmol.
    d/dt(stomach_cip) <- -kge * stomach_cip
    d/dt(transit1_cip) <- kge * stomach_cip - kt * transit1_cip - ka_cip * transit1_cip
    d/dt(transit2_cip) <- kt * transit1_cip - kt * transit2_cip - ka_cip * transit2_cip
    d/dt(transit3_cip) <- kt * transit2_cip - kt * transit3_cip - ka_cip * transit3_cip
    d/dt(transit4_cip) <- kt * transit3_cip - kt * transit4_cip - ka_cip * transit4_cip
    d/dt(transit5_cip) <- kt * transit4_cip - kt * transit5_cip - ka_cip * transit5_cip
    d/dt(transit6_cip) <- kt * transit5_cip - kt * transit6_cip - ka_cip * transit6_cip
    d/dt(transit7_cip) <- kt * transit6_cip - kt * transit7_cip - ka_cip * transit7_cip
    d/dt(colon_cip) <- kt * transit7_cip
    d/dt(adipose_cip) <- q_adipose * (arterial_cip - cv_adipose_cip) / v_adipose
    d/dt(bone_cip) <- q_bone * (arterial_cip - cv_bone_cip) / v_bone
    d/dt(brain_cip) <- q_brain * (arterial_cip - cv_brain_cip) / v_brain
    d/dt(gut_cip) <- q_gut * (arterial_cip - cv_gut_cip) / v_gut+
      (ka_cip *
      (transit1_cip + transit2_cip + transit3_cip + transit4_cip + transit5_cip + transit6_cip + transit7_cip)) / v_gut
    d/dt(heart_cip) <- q_heart * (arterial_cip - cv_heart_cip) / v_heart
    d/dt(kidney_cip) <- q_kidney * (arterial_cip - cv_kidney_cip) / v_kidney -
      RENALFUNC_REL * eren_cip * q_kidney * arterial_cip / v_kidney
    d/dt(liver_cip) <- (q_ha * arterial_cip + q_gut * cv_gut_cip +
      q_spleen * cv_spleen_cip - q_liver * cv_liver_cip) / v_liver -
      (v_cip_1a2)
    d/dt(lung_cip) <- q_co * (venous_cip - cv_lung_cip) / v_lung
    d/dt(muscle_cip) <- q_muscle * (arterial_cip - cv_muscle_cip) / v_muscle
    d/dt(other_cip) <- q_other * (arterial_cip - cv_other_cip) / v_other
    d/dt(skin_cip) <- q_skin * (arterial_cip - cv_skin_cip) / v_skin
    d/dt(spleen_cip) <- q_spleen * (arterial_cip - cv_spleen_cip) / v_spleen
    d/dt(arterial_cip) <- q_co * (cv_lung_cip - arterial_cip) / v_arterial
    d/dt(venous_cip) <-
      (q_adipose * cv_adipose_cip +
       q_bone * cv_bone_cip +
       q_brain * cv_brain_cip +
       q_heart * cv_heart_cip +
       q_kidney * cv_kidney_cip +
       q_liver * cv_liver_cip +
       q_muscle * cv_muscle_cip +
       q_other * cv_other_cip +
       q_skin * cv_skin_cip -
       q_co * venous_cip) / v_venous
    d/dt(urine_cip) <- RENALFUNC_REL * eren_cip * q_kidney * arterial_cip
    Ccip <- venous_cip * 1000 * mw_cip
    f(stomach_cip) <- 1 / mw_cip   # dose in mg -> mmol
    f(venous_cip) <- 1 / (mw_cip * v_venous)   # IV dose in mg -> mmol/mL

    # Navid 2016 is a deterministic PBPK: it reports no between-subject
    # variability and no residual-error model, so neither is encoded.
  })
}
