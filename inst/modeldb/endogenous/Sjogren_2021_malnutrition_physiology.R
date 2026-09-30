Sjogren_2021_malnutrition_physiology <- function() {
  description <- paste(
    "PBPK system-parameter (physiology) model of malnutrition",
    "(Sjogren 2021 Pharmaceutics). This is the SYSTEM layer of a",
    "whole-body PBPK framework, not a drug model: it has no drug, no",
    "dosing, no compartments and no ODEs. It returns the paper's own",
    "physiological scaling parameters (PSPs; Table 2), the",
    "multiplicative factors that transform the physiology of a",
    "non-malnourished individual into that of a malnourished one at",
    "three levels of malnutrition (mild, intermediate, severe). The",
    "PSPs scale organ and tissue volumes (bone, brain, fat, heart,",
    "kidney, liver, muscle, pancreas, skin, spleen, and one shared",
    "factor for gonads, intestines, lung and stomach), the arterial,",
    "venous and portal blood volumes, plasma-protein (albumin and",
    "alpha-1-acid glycoprotein) concentration and hematocrit. All",
    "PSPs are 1 for a non-malnourished individual. The PSPs were",
    "derived from adult data and applied by the authors to virtual",
    "pediatric populations, i.e. they are relative changes to be",
    "multiplied onto an age-appropriate non-malnourished physiology",
    "(for example the adult Schlender_2016_aging_physiology model).",
    "The drug-level PK-Sim PBPK models the paper uses to evaluate the",
    "PSPs (caffeine, cefoxitin, ciprofloxacin, lumefantrine,",
    "pyrimethamine and sulfadoxine) are NOT included: the paper gives",
    "only compound-card inputs, names the PK-Sim Standard or Rodgers",
    "and Rowland partition methods without values, and deposited no",
    "project file, so the drug layer is not reproducible outside that",
    "platform.",
    sep = " "
  )
  reference <- paste(
    "Sjogren E, Tarning J, Barnes KI, Jonsson EN.",
    "A Physiologically-Based Pharmacokinetic Framework for Prediction",
    "of Drug Exposure in Malnourished Children.",
    "Pharmaceutics. 2021;13(2):204.",
    "doi:10.3390/pharmaceutics13020204.",
    sep = " "
  )
  vignette <- "Sjogren_2021_malnutrition_physiology"
  units <- list(
    time = "n/a (the scaling parameters are time-invariant)",
    dosing = "n/a (no exogenous dosing; PBPK system-physiology model)",
    concentration = paste(
      "n/a (no drug concentration). Every output is a dimensionless",
      "physiological scaling parameter: the ratio of the malnourished",
      "to the non-malnourished value of the organ volume, blood",
      "volume, plasma-protein concentration or hematocrit it scales."
    )
  )

  covariateData <- list(
    MAL_NOURISH_MILD = list(
      description = "Mild malnutrition indicator (1 = mild, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not malnourished; with MAL_NOURISH_MOD and MAL_NOURISH_SEV also 0)",
      notes = paste(
        "Mutually exclusive with MAL_NOURISH_MOD and MAL_NOURISH_SEV; all",
        "three 0 selects the non-malnourished physiology (every PSP = 1).",
        "Sjogren 2021 takes its malnutrition levels from Barac-Nieto et",
        "al. (ref 18): 'mild nutritional impairment' is a body weight /",
        "height ratio of 89.5 % of standard (Table 1). The paper scales",
        "whole populations to one level rather than classifying",
        "individuals, so the indicator is a simulation-scenario switch."
      ),
      source_name = "M (mild nutritional impairment)"
    ),
    MAL_NOURISH_MOD = list(
      description = "Intermediate (moderate) malnutrition indicator (1 = intermediate, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not malnourished; with MAL_NOURISH_MILD and MAL_NOURISH_SEV also 0)",
      notes = paste(
        "Mutually exclusive with MAL_NOURISH_MILD and MAL_NOURISH_SEV.",
        "The paper's 'intermediate nutritional impairment' (Barac-Nieto,",
        "ref 18), a body weight / height ratio of 82.7 % of standard",
        "(Table 1); encoded on the _MOD member of the MILD / MOD / SEV",
        "severity-suffix family."
      ),
      source_name = "I (intermediate nutritional impairment)"
    ),
    MAL_NOURISH_SEV = list(
      description = "Severe malnutrition indicator (1 = severe, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not malnourished; with MAL_NOURISH_MILD and MAL_NOURISH_MOD also 0)",
      notes = paste(
        "Mutually exclusive with MAL_NOURISH_MILD and MAL_NOURISH_MOD.",
        "The paper's 'severe nutritional impairment' (Barac-Nieto, ref",
        "18), a body weight / height ratio of 73.9 % of standard (Table",
        "1). Severe was the level the authors selected for every PK",
        "evaluation (Results 3.1): applied to a virtual 46-48-month-old",
        "population it gives a mean weight-for-height z-score of about",
        "-3, the WHO threshold for severe acute malnutrition."
      ),
      source_name = "S (severe nutritional impairment)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 2L,
    age_range = "Adults (source data for the PSPs); applied by the paper to virtual children aged 6 months to 10 years",
    weight_range = "Malnourished adults 42.5 to 52.0 kg (Barac-Nieto) and 46.3 kg underweight vs 70.9 kg reference (Bosy-Westphal) (Table 1)",
    height_range = "156 to 157 cm (Barac-Nieto) and 165 to 178 cm (Bosy-Westphal) (Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Colombian (Barac-Nieto, males, n = 49) and Dutch (Bosy-Westphal, males and females, n = 39)",
    disease_state = paste(
      "Chronic adult undernutrition (Barac-Nieto: mild, intermediate",
      "and severe nutritional impairment; Bosy-Westphal: underweight vs",
      "intermediate weight). No suitable data set for malnourished",
      "children was found (Methods 2.3.2)."
    ),
    dose_range = "n/a (no exogenous drug PK is modelled by this system layer)",
    regions = "Colombia and the Netherlands (PSP source data); PK evaluation against studies in Nigeria, South Africa, Kenya, Mali and other African sites",
    notes = paste(
      "n_subjects and n_studies describe the two adult studies whose",
      "measurements (Table 1) the PSPs were computed from: Barac-Nieto",
      "et al. (ref 18; hematology, fat mass, and the three malnutrition",
      "levels) and Bosy-Westphal et al. (ref 19; organ masses, bone",
      "mineral content, lean soft-tissue and skeletal-muscle mass). The",
      "sex split of the Dutch cohort is not reported, so",
      "sex_female_pct is NA. Organ PSPs are ratios of malnourished to",
      "non-malnourished measurements, corrected for the height",
      "difference between the study groups with ICRP height-typical",
      "organ values (Equation 1); the ICRP 30-year-old man (176 cm,",
      "73 kg) is the non-malnourished reference unless stated",
      "otherwise, and mild / intermediate values missing from the",
      "source data were extrapolated linearly in the body weight /",
      "height ratio. Fat and muscle PSPs were adjusted so the scaled",
      "organ weights add up to the target body weight. The paper",
      "evaluated the PSPs by comparing virtual-population weight-for-",
      "height and weight-for-age z-scores with 131 severely",
      "malnourished and 266 non-malnourished Malian and Nigerien",
      "children aged 6-59 months (ref 23)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Sjogren 2021 Table 2 ('Derived physiological scaling parameters
    # for translation of physiological changes at different levels of
    # malnutrition'), page 9. Each PSP is the paper's own derived
    # point value with no uncertainty, so every value is fixed(). The
    # 'Not Malnourished' column is 1 for every component and is
    # structural (see model()). Suffixes name the malnutrition
    # stratum: _mild, _mod (the paper's 'Intermediate'), _sev.
    # -----------------------------------------------------------------
    psp_bone_mild <- fixed(0.947)
    label("PSP for bone volume, mild malnutrition (ratio)") # Table 2, row 'Bone', column 'Mild'
    psp_bone_mod <- fixed(0.913)
    label("PSP for bone volume, intermediate malnutrition (ratio)") # Table 2, row 'Bone', column 'Intermediate'
    psp_bone_sev <- fixed(0.869)
    label("PSP for bone volume, severe malnutrition (ratio)") # Table 2, row 'Bone', column 'Severe'

    psp_brain_mild <- fixed(0.918)
    label("PSP for brain volume, mild malnutrition (ratio)") # Table 2, row 'Brain', column 'Mild'
    psp_brain_mod <- fixed(0.865)
    label("PSP for brain volume, intermediate malnutrition (ratio)") # Table 2, row 'Brain', column 'Intermediate'
    psp_brain_sev <- fixed(0.797)
    label("PSP for brain volume, severe malnutrition (ratio)") # Table 2, row 'Brain', column 'Severe'

    psp_fat_mild <- fixed(0.817)
    label("PSP for fat (adipose) volume, mild malnutrition (ratio)") # Table 2, row 'Fat', column 'Mild'
    psp_fat_mod <- fixed(0.822)
    label("PSP for fat (adipose) volume, intermediate malnutrition (ratio)") # Table 2, row 'Fat', column 'Intermediate' (above the Mild value as printed; fat and muscle were adjusted to hit the target body weight, Methods 2.3.2)
    psp_fat_sev <- fixed(0.624)
    label("PSP for fat (adipose) volume, severe malnutrition (ratio)") # Table 2, row 'Fat', column 'Severe'

    psp_leanorg_mild <- fixed(0.936)
    label("PSP for gonads, intestines, lung and stomach volumes, mild malnutrition (ratio)") # Table 2, row 'Gonads, intestines, lung, stomach', column 'Mild'
    psp_leanorg_mod <- fixed(0.894)
    label("PSP for gonads, intestines, lung and stomach volumes, intermediate malnutrition (ratio)") # Table 2, row 'Gonads, intestines, lung, stomach', column 'Intermediate'
    psp_leanorg_sev <- fixed(0.84)
    label("PSP for gonads, intestines, lung and stomach volumes, severe malnutrition (ratio)") # Table 2, row 'Gonads, intestines, lung, stomach', column 'Severe'

    psp_heart_mild <- fixed(0.902)
    label("PSP for heart volume, mild malnutrition (ratio)") # Table 2, row 'Heart', column 'Mild'
    psp_heart_mod <- fixed(0.839)
    label("PSP for heart volume, intermediate malnutrition (ratio)") # Table 2, row 'Heart', column 'Intermediate'
    psp_heart_sev <- fixed(0.758)
    label("PSP for heart volume, severe malnutrition (ratio)") # Table 2, row 'Heart', column 'Severe'

    psp_kidney_mild <- fixed(0.874)
    label("PSP for kidney volume, mild malnutrition (ratio)") # Table 2, row 'Kidney', column 'Mild'
    psp_kidney_mod <- fixed(0.792)
    label("PSP for kidney volume, intermediate malnutrition (ratio)") # Table 2, row 'Kidney', column 'Intermediate'
    psp_kidney_sev <- fixed(0.686)
    label("PSP for kidney volume, severe malnutrition (ratio)") # Table 2, row 'Kidney', column 'Severe'

    psp_liver_mild <- fixed(0.872)
    label("PSP for liver volume, mild malnutrition (ratio)") # Table 2, row 'Liver', column 'Mild'
    psp_liver_mod <- fixed(0.789)
    label("PSP for liver volume, intermediate malnutrition (ratio)") # Table 2, row 'Liver', column 'Intermediate'
    psp_liver_sev <- fixed(0.682)
    label("PSP for liver volume, severe malnutrition (ratio)") # Table 2, row 'Liver', column 'Severe'

    psp_muscle_mild <- fixed(0.893)
    label("PSP for muscle volume, mild malnutrition (ratio)") # Table 2, row 'Muscle', column 'Mild'
    psp_muscle_mod <- fixed(0.771)
    label("PSP for muscle volume, intermediate malnutrition (ratio)") # Table 2, row 'Muscle', column 'Intermediate'
    psp_muscle_sev <- fixed(0.715)
    label("PSP for muscle volume, severe malnutrition (ratio)") # Table 2, row 'Muscle', column 'Severe'

    psp_pancreas_mild <- fixed(0.936)
    label("PSP for pancreas volume, mild malnutrition (ratio)") # Table 2, row 'Pancreas', column 'Mild'
    psp_pancreas_mod <- fixed(0.894)
    label("PSP for pancreas volume, intermediate malnutrition (ratio)") # Table 2, row 'Pancreas', column 'Intermediate'
    psp_pancreas_sev <- fixed(0.84)
    label("PSP for pancreas volume, severe malnutrition (ratio)") # Table 2, row 'Pancreas', column 'Severe'

    psp_skin_mild <- fixed(0.954)
    label("PSP for skin volume, mild malnutrition (ratio)") # Table 2, row 'Skin', column 'Mild'
    psp_skin_mod <- fixed(0.922)
    label("PSP for skin volume, intermediate malnutrition (ratio)") # Table 2, row 'Skin', column 'Intermediate'
    psp_skin_sev <- fixed(0.879)
    label("PSP for skin volume, severe malnutrition (ratio)") # Table 2, row 'Skin', column 'Severe'

    psp_spleen_mild <- fixed(0.844)
    label("PSP for spleen volume, mild malnutrition (ratio)") # Table 2, row 'Spleen', column 'Mild'
    psp_spleen_mod <- fixed(0.743)
    label("PSP for spleen volume, intermediate malnutrition (ratio)") # Table 2, row 'Spleen', column 'Intermediate'
    psp_spleen_sev <- fixed(0.612)
    label("PSP for spleen volume, severe malnutrition (ratio)") # Table 2, row 'Spleen', column 'Severe'

    psp_blood_art_mild <- fixed(1.03)
    label("PSP for arterial blood volume, mild malnutrition (ratio)") # Table 2, row 'Blood / arterial', column 'Mild'
    psp_blood_art_mod <- fixed(0.992)
    label("PSP for arterial blood volume, intermediate malnutrition (ratio)") # Table 2, row 'Blood / arterial', column 'Intermediate'
    psp_blood_art_sev <- fixed(0.833)
    label("PSP for arterial blood volume, severe malnutrition (ratio)") # Table 2, row 'Blood / arterial', column 'Severe'

    psp_blood_ven_mild <- fixed(1.02)
    label("PSP for venous blood volume, mild malnutrition (ratio)") # Table 2, row 'Blood / venous', column 'Mild'
    psp_blood_ven_mod <- fixed(0.979)
    label("PSP for venous blood volume, intermediate malnutrition (ratio)") # Table 2, row 'Blood / venous', column 'Intermediate'
    psp_blood_ven_sev <- fixed(0.822)
    label("PSP for venous blood volume, severe malnutrition (ratio)") # Table 2, row 'Blood / venous', column 'Severe'

    psp_blood_portal_mild <- fixed(1.02)
    label("PSP for portal-vein blood volume, mild malnutrition (ratio)") # Table 2, row 'Blood / portal vein', column 'Mild'
    psp_blood_portal_mod <- fixed(0.982)
    label("PSP for portal-vein blood volume, intermediate malnutrition (ratio)") # Table 2, row 'Blood / portal vein', column 'Intermediate'
    psp_blood_portal_sev <- fixed(0.825)
    label("PSP for portal-vein blood volume, severe malnutrition (ratio)") # Table 2, row 'Blood / portal vein', column 'Severe'

    psp_plasmaprot_mild <- fixed(0.894)
    label("PSP for plasma-protein (albumin, AAG) concentration, mild malnutrition (ratio)") # Table 2, row 'Plasma proteins', column 'Mild' (= 3.8 / 4.25, Table 1 albumin over the 4.25 g/100 mL reference, Methods 2.3.2)
    psp_plasmaprot_mod <- fixed(0.706)
    label("PSP for plasma-protein (albumin, AAG) concentration, intermediate malnutrition (ratio)") # Table 2, row 'Plasma proteins', column 'Intermediate' (= 3 / 4.25)
    psp_plasmaprot_sev <- fixed(0.494)
    label("PSP for plasma-protein (albumin, AAG) concentration, severe malnutrition (ratio)") # Table 2, row 'Plasma proteins', column 'Severe' (= 2.1 / 4.25)

    psp_hct_mild <- fixed(0.945)
    label("PSP for hematocrit, mild malnutrition (ratio)") # Table 2, row 'Hematocrit', column 'Mild'
    psp_hct_mod <- fixed(0.791)
    label("PSP for hematocrit, intermediate malnutrition (ratio)") # Table 2, row 'Hematocrit', column 'Intermediate'
    psp_hct_sev <- fixed(0.681)
    label("PSP for hematocrit, severe malnutrition (ratio)") # Table 2, row 'Hematocrit', column 'Severe'
  })

  model({
    # The three indicators are mutually exclusive; all three 0 is the
    # non-malnourished physiology, for which Table 2 gives PSP = 1 for
    # every component. Each output is therefore 1 plus the selected
    # stratum's departure from 1.

    # ---- Organ and tissue volumes (Table 2) ----
    psp_bone <- 1 + MAL_NOURISH_MILD * (psp_bone_mild - 1) +
      MAL_NOURISH_MOD * (psp_bone_mod - 1) +
      MAL_NOURISH_SEV * (psp_bone_sev - 1)
    psp_brain <- 1 + MAL_NOURISH_MILD * (psp_brain_mild - 1) +
      MAL_NOURISH_MOD * (psp_brain_mod - 1) +
      MAL_NOURISH_SEV * (psp_brain_sev - 1)
    psp_fat <- 1 + MAL_NOURISH_MILD * (psp_fat_mild - 1) +
      MAL_NOURISH_MOD * (psp_fat_mod - 1) +
      MAL_NOURISH_SEV * (psp_fat_sev - 1)
    psp_heart <- 1 + MAL_NOURISH_MILD * (psp_heart_mild - 1) +
      MAL_NOURISH_MOD * (psp_heart_mod - 1) +
      MAL_NOURISH_SEV * (psp_heart_sev - 1)
    psp_kidney <- 1 + MAL_NOURISH_MILD * (psp_kidney_mild - 1) +
      MAL_NOURISH_MOD * (psp_kidney_mod - 1) +
      MAL_NOURISH_SEV * (psp_kidney_sev - 1)
    psp_liver <- 1 + MAL_NOURISH_MILD * (psp_liver_mild - 1) +
      MAL_NOURISH_MOD * (psp_liver_mod - 1) +
      MAL_NOURISH_SEV * (psp_liver_sev - 1)
    psp_muscle <- 1 + MAL_NOURISH_MILD * (psp_muscle_mild - 1) +
      MAL_NOURISH_MOD * (psp_muscle_mod - 1) +
      MAL_NOURISH_SEV * (psp_muscle_sev - 1)
    psp_pancreas <- 1 + MAL_NOURISH_MILD * (psp_pancreas_mild - 1) +
      MAL_NOURISH_MOD * (psp_pancreas_mod - 1) +
      MAL_NOURISH_SEV * (psp_pancreas_sev - 1)
    psp_skin <- 1 + MAL_NOURISH_MILD * (psp_skin_mild - 1) +
      MAL_NOURISH_MOD * (psp_skin_mod - 1) +
      MAL_NOURISH_SEV * (psp_skin_sev - 1)
    psp_spleen <- 1 + MAL_NOURISH_MILD * (psp_spleen_mild - 1) +
      MAL_NOURISH_MOD * (psp_spleen_mod - 1) +
      MAL_NOURISH_SEV * (psp_spleen_sev - 1)

    # Table 2 prints ONE row for 'Gonads, intestines, lung, stomach':
    # lean body mass was the surrogate for all four (Methods 2.3.2).
    # The per-organ outputs below share that single factor.
    psp_leanorg <- 1 + MAL_NOURISH_MILD * (psp_leanorg_mild - 1) +
      MAL_NOURISH_MOD * (psp_leanorg_mod - 1) +
      MAL_NOURISH_SEV * (psp_leanorg_sev - 1)
    psp_gonads <- psp_leanorg
    psp_intestine <- psp_leanorg
    psp_lung <- psp_leanorg
    psp_stomach <- psp_leanorg

    # ---- Blood volumes (Table 2, 'Blood' rows) ----
    # The arterial, venous and portal volumes were calculated from the
    # PSP for total blood volume (Methods 2.3.2); the paper prints the
    # three compartment PSPs, not the total.
    psp_blood_art <- 1 + MAL_NOURISH_MILD * (psp_blood_art_mild - 1) +
      MAL_NOURISH_MOD * (psp_blood_art_mod - 1) +
      MAL_NOURISH_SEV * (psp_blood_art_sev - 1)
    psp_blood_ven <- 1 + MAL_NOURISH_MILD * (psp_blood_ven_mild - 1) +
      MAL_NOURISH_MOD * (psp_blood_ven_mod - 1) +
      MAL_NOURISH_SEV * (psp_blood_ven_sev - 1)
    psp_blood_portal <- 1 + MAL_NOURISH_MILD * (psp_blood_portal_mild - 1) +
      MAL_NOURISH_MOD * (psp_blood_portal_mod - 1) +
      MAL_NOURISH_SEV * (psp_blood_portal_sev - 1)

    # ---- Hematology (Table 2) ----
    # The plasma-protein PSP is the serum-albumin ratio; the paper
    # assumes an equal effect of nutritional status on alpha-1-acid
    # glycoprotein (Methods 2.3.2), and uses albumin as the binding
    # surrogate for lumefantrine's lipoprotein binding (Results 3.2.4).
    psp_plasmaprot <- 1 + MAL_NOURISH_MILD * (psp_plasmaprot_mild - 1) +
      MAL_NOURISH_MOD * (psp_plasmaprot_mod - 1) +
      MAL_NOURISH_SEV * (psp_plasmaprot_sev - 1)
    psp_albumin <- psp_plasmaprot
    psp_aag <- psp_plasmaprot
    psp_hct <- 1 + MAL_NOURISH_MILD * (psp_hct_mild - 1) +
      MAL_NOURISH_MOD * (psp_hct_mod - 1) +
      MAL_NOURISH_SEV * (psp_hct_sev - 1)
  })
}
