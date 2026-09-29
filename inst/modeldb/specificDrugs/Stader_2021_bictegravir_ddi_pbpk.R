Stader_2021_bictegravir_ddi_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, Stader et al. Matlab 2017a framework, as deposited",
    "with the paper). Drug-drug interactions of oral bictegravir with",
    "inhibitors and inducers of CYP3A and UGT1A1 in adults aged 20 to 50",
    "years. Three drugs are carried simultaneously, each with the",
    "framework's full 63-state whole-body model (16 organs with vascular,",
    "interstitial and intracellular sub-compartments, venous and arterial",
    "blood, and a compartmental absorption and transit gut) plus 12",
    "enzyme-turnover states (hepatic CYP3A4, 2C19, 2D6, 2C8, 1A2, 2A6, 2B6,",
    "2J2 and UGT1A1, intestinal CYP3A4 per segment): 225 ODEs. The victim",
    "drug (bictegravir) takes the bare state names; the two co-administered",
    "drugs are parameterised slots (_perpetrator, _perpetrator2) that ship",
    "set to voriconazole and rifampicin, and the drugLibrary metadata holds",
    "the as-run inputs of all 13 drug models of the deposited library so that",
    "any of the paper's scenarios can be reproduced by substituting them.",
    "Interactions act on enzyme synthesis (induction), degradation",
    "(mechanism-based inactivation) and turnover rate (competitive",
    "inhibition), with the perpetrators' own concentrations taken from",
    "their unperturbed slots as in the framework. ddi_on = 0 turns the",
    "victim block into the framework's victim-alone control slot. Organ",
    "weights, blood and lymph flows, plasma proteins, GFR, microsomal protein",
    "and enzyme abundances are age-, sex-, height- and weight-dependent",
    "regressions of the deposited virtual-population generator, whose",
    "per-subject random draws are carried as fixed etas. Cc is the victim's",
    "reported plasma concentration, which the framework reads from the",
    "venous-blood state.",
    sep = " "
  )
  reference <- paste(
    "Stader F, Battegay M, Marzolini C. Physiologically-Based Pharmacokinetic",
    "Modeling to Support the Clinical Management of Drug-Drug Interactions",
    "With Bictegravir. Clin Pharmacol Ther. 2021;110(5):1231-1239.",
    "doi:10.1002/cpt.2221. Model code: Supplementary Material s002",
    "(CPT_Matlab_Code, the Matlab source of the PBPK framework with the drug",
    "library Drug/DrugLibrary/*.m); voriconazole and cobicistat inputs:",
    "Supplementary Table S1. The bictegravir drug model is described in",
    "Stader F et al. Clin Pharmacol Ther. 2021;109(4):1025-1029,",
    "doi:10.1002/cpt.2178.",
    sep = " "
  )
  vignette <- "Stader_2021_bictegravir_ddi_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the age regressions of the deposited virtual-population",
        "generator: height and weight means (supplied here as covariates),",
        "lung, adipose, brain, gonad, heart, kidney, muscle and skin weights,",
        "blood weight, albumin, brain vascular fraction, cardiac output,",
        "adipose, kidney and liver blood-flow fractions, GFR and microsomal",
        "protein per gram liver. The generator draws integer ages (rounded",
        "Weibull draws, 61.73 scale and 1.55 shape) restricted to 20-99",
        "years; predictions outside 20-99 years are not supported by the",
        "framework.",
        sep = " "
      ),
      source_name = "BODY.Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (0)",
      notes = paste(
        "The deposited code stores sex as BODY.Sex with 1 = female and",
        "0 = male (PBPK_Population_Demographics.m, Generate_Sex: 'assign a",
        "1 to females and 0 to males'), so it maps onto SEXF without",
        "transformation.",
        sep = " "
      ),
      source_name = "BODY.Sex"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual height. The generator draws it as",
        "normal(mean = -0.0039*AGE^2 + 0.238*AGE - 12.5*SEXF + 176,",
        "CV 3.8%) and resets it to the mean when the resulting BMI falls",
        "outside 18.5-30 kg/m^2 (PBPK_Population_Demographics.m); the",
        "validation vignette reproduces that draw. Enters BSA, lung, bone,",
        "brain and gut weights, adipose weight, heart and kidney blood",
        "flows and intestinal segment lengths.",
        sep = " "
      ),
      source_name = "BODY.Height"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual weight. The generator draws it as",
        "normal(mean = -0.0039*AGE^2 + 1.12*HT + 0.611*AGE - 0.424*SEXF - 137,",
        "CV 15.2%) with the same BMI reset as height. Enters BSA, adipose,",
        "gonad and lymph-node weights, total lymph flow, and the remaining",
        "tissue as the body-weight balance.",
        sep = " "
      ),
      source_name = "BODY.Weight"
    )
  )

  # Sub-compartment and enzyme states with no registered canonical (thymus and
  # gonad sub-compartments, the CAT uptake layer and enterocytes, per-tissue
  # isoform enzyme pools, and every state of the second perpetrator slot).
  paper_specific_compartments <- c(
    "gonads_vas",
    "gonads_ew",
    "gonads_iw",
    "thymus_vas",
    "thymus_ew",
    "thymus_iw",
    "duodenum_uptake",
    "jejunum_uptake",
    "ileum_uptake",
    "colon_uptake",
    "duodenum_enterocyte",
    "jejunum_enterocyte",
    "ileum_enterocyte",
    "colon_enterocyte",
    "enzyme_cyp3a4_liver",
    "enzyme_cyp2c19_liver",
    "enzyme_cyp2d6_liver",
    "enzyme_cyp2c8_liver",
    "enzyme_cyp1a2_liver",
    "enzyme_cyp2a6_liver",
    "enzyme_cyp2b6_liver",
    "enzyme_cyp2j2_liver",
    "enzyme_ugt1a1_liver",
    "enzyme_cyp3a4_duodenum",
    "enzyme_cyp3a4_jejunum",
    "enzyme_cyp3a4_ileum",
    "gonads_vas_perpetrator",
    "gonads_ew_perpetrator",
    "gonads_iw_perpetrator",
    "thymus_vas_perpetrator",
    "thymus_ew_perpetrator",
    "thymus_iw_perpetrator",
    "duodenum_uptake_perpetrator",
    "jejunum_uptake_perpetrator",
    "ileum_uptake_perpetrator",
    "colon_uptake_perpetrator",
    "duodenum_enterocyte_perpetrator",
    "jejunum_enterocyte_perpetrator",
    "ileum_enterocyte_perpetrator",
    "colon_enterocyte_perpetrator",
    "enzyme_cyp3a4_liver_perpetrator",
    "enzyme_cyp2c19_liver_perpetrator",
    "enzyme_cyp2d6_liver_perpetrator",
    "enzyme_cyp2c8_liver_perpetrator",
    "enzyme_cyp1a2_liver_perpetrator",
    "enzyme_cyp2a6_liver_perpetrator",
    "enzyme_cyp2b6_liver_perpetrator",
    "enzyme_cyp2j2_liver_perpetrator",
    "enzyme_ugt1a1_liver_perpetrator",
    "enzyme_cyp3a4_duodenum_perpetrator",
    "enzyme_cyp3a4_jejunum_perpetrator",
    "enzyme_cyp3a4_ileum_perpetrator",
    "lung_vas_perpetrator2",
    "lung_ew_perpetrator2",
    "lung_iw_perpetrator2",
    "adipose_vas_perpetrator2",
    "adipose_ew_perpetrator2",
    "adipose_iw_perpetrator2",
    "bone_vas_perpetrator2",
    "bone_ew_perpetrator2",
    "bone_iw_perpetrator2",
    "brain_vas_perpetrator2",
    "brain_ew_perpetrator2",
    "brain_iw_perpetrator2",
    "gonads_vas_perpetrator2",
    "gonads_ew_perpetrator2",
    "gonads_iw_perpetrator2",
    "heart_vas_perpetrator2",
    "heart_ew_perpetrator2",
    "heart_iw_perpetrator2",
    "kidney_vas_perpetrator2",
    "kidney_ew_perpetrator2",
    "kidney_iw_perpetrator2",
    "muscle_vas_perpetrator2",
    "muscle_ew_perpetrator2",
    "muscle_iw_perpetrator2",
    "skin_vas_perpetrator2",
    "skin_ew_perpetrator2",
    "skin_iw_perpetrator2",
    "thymus_vas_perpetrator2",
    "thymus_ew_perpetrator2",
    "thymus_iw_perpetrator2",
    "spleen_vas_perpetrator2",
    "spleen_ew_perpetrator2",
    "spleen_iw_perpetrator2",
    "pancreas_vas_perpetrator2",
    "pancreas_ew_perpetrator2",
    "pancreas_iw_perpetrator2",
    "other_vas_perpetrator2",
    "other_ew_perpetrator2",
    "other_iw_perpetrator2",
    "gut_vas_perpetrator2",
    "gut_ew_perpetrator2",
    "liver_vas_perpetrator2",
    "liver_ew_perpetrator2",
    "liver_iw_perpetrator2",
    "lnode_vas_perpetrator2",
    "lnode_ew_perpetrator2",
    "lnode_iw_perpetrator2",
    "venous_perpetrator2",
    "arterial_perpetrator2",
    "stomach_perpetrator2",
    "duodenum_perpetrator2",
    "jejunum_perpetrator2",
    "ileum_perpetrator2",
    "colon_perpetrator2",
    "duodenum_uptake_perpetrator2",
    "jejunum_uptake_perpetrator2",
    "ileum_uptake_perpetrator2",
    "colon_uptake_perpetrator2",
    "duodenum_enterocyte_perpetrator2",
    "jejunum_enterocyte_perpetrator2",
    "ileum_enterocyte_perpetrator2",
    "colon_enterocyte_perpetrator2",
    "a_feces_perpetrator2",
    "enzyme_cyp3a4_liver_perpetrator2",
    "enzyme_cyp2c19_liver_perpetrator2",
    "enzyme_cyp2d6_liver_perpetrator2",
    "enzyme_cyp2c8_liver_perpetrator2",
    "enzyme_cyp1a2_liver_perpetrator2",
    "enzyme_cyp2a6_liver_perpetrator2",
    "enzyme_cyp2b6_liver_perpetrator2",
    "enzyme_cyp2j2_liver_perpetrator2",
    "enzyme_ugt1a1_liver_perpetrator2",
    "enzyme_cyp3a4_duodenum_perpetrator2",
    "enzyme_cyp3a4_jejunum_perpetrator2",
    "enzyme_cyp3a4_ileum_perpetrator2"
  )

  # The etas are the per-subject random draws of the deposited
  # virtual-population generator, not estimated between-subject variances.
  paper_specific_etas <- c(
    "etahct",
    "etahsa",
    "etaw_adipose",
    "etaw_bone",
    "etaw_brain",
    "etaw_gonads",
    "etaw_heart",
    "etaw_kidney",
    "etaw_muscle",
    "etaw_skin",
    "etaw_thymus",
    "etaw_gut",
    "etaw_spleen",
    "etaw_pancreas",
    "etaw_liver",
    "etaw_lnode",
    "etaw_blood",
    "etaco",
    "etaltot",
    "etagfr",
    "etamppgl",
    "etamppgl_redraw",
    "etacyp3a4_liver",
    "etacyp3a4_liver_redraw",
    "etaugt1a1_liver",
    "etaugt1a1_liver_redraw",
    "etagastric",
    "etasitt",
    "etasitt_redraw",
    "etacolt",
    "etacolt_redraw",
    "etacyp3a4_gut",
    "etacyp3a4_gut_redraw",
    "etaaag",
    "etacyp2c19_liver",
    "etacyp2c19_liver_redraw",
    "etacyp2d6_liver",
    "etacyp2d6_liver_redraw",
    "etacyp2c8_liver",
    "etacyp2c8_liver_redraw",
    "etacyp1a2_liver",
    "etacyp1a2_liver_redraw",
    "etacyp2a6_liver",
    "etacyp2a6_liver_redraw",
    "etacyp2b6_liver",
    "etacyp2b6_liver_redraw",
    "etacyp2j2_liver",
    "etacyp2j2_liver_redraw",
    "etaugt1a4_liver",
    "etaugt1a4_liver_redraw",
    "etacyp2c19_gut",
    "etacyp2c19_gut_redraw",
    "etacyp2d6_gut",
    "etacyp2d6_gut_redraw"
  )

  compartmentData <- list(
    lung_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    bone_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    bone_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    bone_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    brain_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    gonads_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    heart_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    heart_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    heart_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    kidney_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    muscle_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    skin_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    skin_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    skin_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    thymus_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    spleen_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    pancreas_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    other_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    other_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    other_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    gut_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    gut_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    liver_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    liver_ew = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    liver_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    lnode_vas = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    lnode_ew = list(analyte = "bictegravir", units = "mg", specimen = "lymph", verified = TRUE),
    lnode_iw = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "bictegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    stomach = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    ileum = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    colon = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    colon_uptake = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    duodenum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    colon_enterocyte = list(analyte = "bictegravir", units = "mg", specimen = "tissue", verified = TRUE),
    a_feces = list(analyte = "bictegravir", units = "mg", specimen = "faeces", verified = TRUE),
    enzyme_cyp3a4_liver = list(
      analyte = "CYP3A4 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c19_liver = list(
      analyte = "CYP2C19 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2d6_liver = list(
      analyte = "CYP2D6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c8_liver = list(
      analyte = "CYP2C8 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp1a2_liver = list(
      analyte = "CYP1A2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2a6_liver = list(
      analyte = "CYP2A6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2b6_liver = list(
      analyte = "CYP2B6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2j2_liver = list(
      analyte = "CYP2J2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_ugt1a1_liver = list(
      analyte = "UGT1A1 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_duodenum = list(
      analyte = "CYP3A4 abundance (duodenum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_jejunum = list(
      analyte = "CYP3A4 abundance (jejunum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_ileum = list(
      analyte = "CYP3A4 abundance (ileum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    lung_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    lung_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    bone_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    bone_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    bone_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    brain_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    brain_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    gonads_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    heart_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    heart_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    heart_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    kidney_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    muscle_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    skin_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    skin_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    skin_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    thymus_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    spleen_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    pancreas_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    other_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    other_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    other_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    gut_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    gut_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    liver_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    liver_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    liver_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    lnode_vas_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    lnode_ew_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "lymph", verified = TRUE),
    lnode_iw_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    venous_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    stomach_perpetrator = list(
      analyte = "voriconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_perpetrator = list(
      analyte = "voriconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum_perpetrator = list(
      analyte = "voriconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    colon_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum_uptake_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_uptake_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_uptake_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    colon_uptake_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    duodenum_enterocyte_perpetrator = list(
      analyte = "voriconazole",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    jejunum_enterocyte_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_enterocyte_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    colon_enterocyte_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    a_feces_perpetrator = list(analyte = "voriconazole", units = "mg", specimen = "faeces", verified = TRUE),
    enzyme_cyp3a4_liver_perpetrator = list(
      analyte = "CYP3A4 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c19_liver_perpetrator = list(
      analyte = "CYP2C19 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2d6_liver_perpetrator = list(
      analyte = "CYP2D6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c8_liver_perpetrator = list(
      analyte = "CYP2C8 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp1a2_liver_perpetrator = list(
      analyte = "CYP1A2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2a6_liver_perpetrator = list(
      analyte = "CYP2A6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2b6_liver_perpetrator = list(
      analyte = "CYP2B6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2j2_liver_perpetrator = list(
      analyte = "CYP2J2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_ugt1a1_liver_perpetrator = list(
      analyte = "UGT1A1 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_duodenum_perpetrator = list(
      analyte = "CYP3A4 abundance (duodenum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_jejunum_perpetrator = list(
      analyte = "CYP3A4 abundance (jejunum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_ileum_perpetrator = list(
      analyte = "CYP3A4 abundance (ileum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    lung_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    bone_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    bone_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    bone_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    brain_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    brain_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    gonads_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    gonads_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    heart_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    heart_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    heart_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    kidney_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    muscle_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    skin_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    skin_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    skin_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    thymus_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    thymus_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    spleen_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    pancreas_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    other_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    other_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    other_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    gut_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    gut_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    liver_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    liver_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    liver_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    lnode_vas_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    lnode_ew_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "lymph", verified = TRUE),
    lnode_iw_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    venous_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    stomach_perpetrator2 = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_perpetrator2 = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum_perpetrator2 = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    colon_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum_uptake_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_uptake_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_uptake_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    colon_uptake_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    duodenum_enterocyte_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    jejunum_enterocyte_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    ileum_enterocyte_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    colon_enterocyte_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    a_feces_perpetrator2 = list(analyte = "rifampicin", units = "mg", specimen = "faeces", verified = TRUE),
    enzyme_cyp3a4_liver_perpetrator2 = list(
      analyte = "CYP3A4 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c19_liver_perpetrator2 = list(
      analyte = "CYP2C19 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2d6_liver_perpetrator2 = list(
      analyte = "CYP2D6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2c8_liver_perpetrator2 = list(
      analyte = "CYP2C8 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp1a2_liver_perpetrator2 = list(
      analyte = "CYP1A2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2a6_liver_perpetrator2 = list(
      analyte = "CYP2A6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2b6_liver_perpetrator2 = list(
      analyte = "CYP2B6 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp2j2_liver_perpetrator2 = list(
      analyte = "CYP2J2 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_ugt1a1_liver_perpetrator2 = list(
      analyte = "UGT1A1 abundance (liver)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_duodenum_perpetrator2 = list(
      analyte = "CYP3A4 abundance (duodenum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_jejunum_perpetrator2 = list(
      analyte = "CYP3A4 abundance (jejunum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_cyp3a4_ileum_perpetrator2 = list(
      analyte = "CYP3A4 abundance (ileum)",
      units = "fraction of baseline",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = NA_integer_,
    age_range = "20-50 years (virtual individuals)",
    sex_female_pct = 50,
    disease_state = "Healthy virtual adults; clinical DDI data from healthy volunteers",
    dose_range = "Bictegravir 50 mg once daily (prospective scenarios); perpetrators at the Table S3 doses",
    regions = NA_character_,
    notes = paste(
      "A PBPK model, not a fit to individual data. Each scenario was",
      "simulated in 100 virtual individuals (10 trials of 10, 50% women)",
      "aged 20-50 years. Perpetrators were given for 14 days at the",
      "Table S3 doses and bictegravir 50 mg once daily from day 7; DDI",
      "magnitudes are AUC ratios on the first and seventh day of",
      "bictegravir. The framework's DDI predictions were verified against",
      "clinical DDI studies of bictegravir with voriconazole,",
      "darunavir/cobicistat, atazanavir, atazanavir/cobicistat and",
      "rifampicin (Table 1; data from Gilead Sciences), and the new",
      "voriconazole and cobicistat models against published healthy-volunteer",
      "studies (Table S2).",
      sep = " "
    )
  )

  # As-run inputs of all 13 drug models of the deposited library
  # (Drug/DrugLibrary/<drug>.m), one column per drug in library order
  # (PBPK_DefineParameters.m; clarithromycin is not registered there and is
  # appended). Values are the ones the framework uses: a KM, Kapp, IC50 or
  # tissue scalar of 0 is replaced by 1, and a Ki of 0 means no inhibition.
  # Parameter names match the ini() stems; append _perpetrator or
  # _perpetrator2 to load a drug into a co-administered slot.
  drugLibrary <- data.frame(
    parameter = c(
      "mw",
      "logp",
      "dtype",
      "pka1",
      "pka2",
      "bp",
      "fu",
      "pb_aag",
      "papp",
      "kperup",
      "fgp",
      "fabscolon",
      "lagtime",
      "kpscalar",
      "jin_all",
      "jin_adipose",
      "jin_muscle",
      "jin_liver",
      "fin_all",
      "vmax1_cyp3a4",
      "km1_cyp3a4",
      "vmax2_cyp3a4",
      "km2_cyp3a4",
      "clint_cyp3a4",
      "vmax1_cyp2c19",
      "km1_cyp2c19",
      "clint_cyp2c19",
      "vmax1_cyp2d6",
      "km1_cyp2d6",
      "clint_cyp2d6",
      "clint_cyp2c8",
      "clint_cyp1a2",
      "clint_cyp2a6",
      "clint_cyp2b6",
      "clint_cyp2j2",
      "clint_ugt1a1",
      "vmax_ugt1a4",
      "km_ugt1a4",
      "clint_hep",
      "clrenal",
      "clbile",
      "ki_cyp3a4",
      "ki_cyp2c19",
      "ki_cyp2d6",
      "ki_cyp2c8",
      "ki_cyp1a2",
      "ki_ugt1a1",
      "kinact_cyp3a4",
      "kapp_cyp3a4",
      "kinact_cyp2j2",
      "kapp_cyp2j2",
      "indmax_cyp3a4",
      "ic50_cyp3a4",
      "indmax_cyp2b6",
      "ic50_cyp2b6",
      "indmax_ugt1a1",
      "ic50_ugt1a1"
    ),
    midazolam = c(
      325.8,
      3.89,
      1,
      6.15,
      0,
      0.6,
      0.032,
      0,
      210,
      0,
      0.005,
      1,
      0,
      1,
      0.67,
      1,
      1,
      3,
      1,
      5.23,
      2.16,
      5.2,
      31.8,
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      30,
      64,
      0,
      0.085,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    ketoconazole = c(
      531.43,
      4.04,
      2,
      2.94,
      6.51,
      0.62,
      0.029,
      0,
      495,
      2,
      1,
      1,
      0,
      1,
      0.5,
      2,
      1,
      1,
      1,
      0,
      1,
      0,
      1,
      0.5238,
      0,
      1,
      0,
      0,
      1,
      0.4296,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0.015,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    voriconazole = c(
      349.31,
      1.8,
      1,
      1.76,
      0,
      1.23,
      0.42,
      0,
      28.1,
      0,
      1,
      1,
      0,
      2.75,
      1,
      1,
      1,
      1,
      1,
      1212,
      15,
      0,
      1,
      0,
      4.19,
      3.5,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0.346,
      0.096,
      0,
      0.66,
      5.1,
      0,
      0,
      0,
      0,
      9.33,
      0.015,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    nilotinib = c(
      529.5,
      5,
      2,
      5.35,
      3.9,
      0.68,
      0.016,
      1,
      5.99,
      0.27,
      0.2,
      0,
      0,
      1,
      2,
      1,
      1,
      0.08333333,
      1,
      0,
      1,
      0,
      1,
      0.157,
      0,
      1,
      0,
      0,
      1,
      0,
      0.127,
      0.018,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0.448,
      0,
      0,
      0.236,
      0,
      0.19,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    rifampicin = c(
      823,
      3.28,
      6,
      1.7,
      7.9,
      0.9,
      0.15,
      0,
      1.472,
      3,
      1,
      1,
      0,
      0.22,
      0.05,
      1,
      1,
      10,
      1,
      0,
      1,
      0,
      1,
      0.0036,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      6.55,
      1.2,
      0,
      10.5,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      16.68,
      0.32,
      0,
      1,
      1.668,
      0.32
    ),
    bictegravir = c(
      449.39,
      1.28,
      3,
      9.81,
      0,
      0.64,
      0.0025,
      0,
      24.6,
      0,
      1,
      1,
      0,
      1,
      0.7,
      1,
      1,
      2,
      1,
      0,
      1,
      0,
      1,
      0.114,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0.292,
      0,
      1,
      3.993,
      0.0043,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    atazanavir = c(
      705,
      4.5,
      1,
      5.62,
      0,
      0.75,
      0.14,
      1,
      19.5,
      1,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      1.5,
      1,
      0,
      1,
      0,
      1,
      6.57,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      1.94,
      0,
      0,
      2.35,
      0,
      0,
      2.1,
      12.1,
      1.9,
      30,
      0.84,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    darunavir = c(
      547.7,
      1.8,
      1,
      2.39,
      0,
      0.64,
      0.06,
      1,
      5.5,
      0.5,
      1,
      0,
      0,
      1,
      0.7,
      0.1,
      1,
      1,
      1.5,
      0,
      1,
      0,
      1,
      3.35,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      4.88,
      0,
      0,
      0.44,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      2.2,
      0.18,
      0,
      1,
      0,
      1
    ),
    ritonavir = c(
      720.95,
      4.3,
      1,
      2,
      0,
      0.587,
      0.015,
      1,
      2.1,
      0,
      0.000009,
      1,
      1,
      1,
      1,
      0.4,
      0.4,
      1.3,
      1,
      1.37,
      0.068,
      0,
      1,
      0,
      0,
      1,
      0,
      0.93,
      1,
      0,
      0,
      0,
      0,
      0,
      2,
      0,
      0,
      1,
      0,
      0.32,
      0,
      0.02928,
      0.15,
      2.9,
      0,
      0,
      0,
      192,
      0.091,
      4.941176,
      0.4641,
      13.4,
      0.44,
      0,
      1,
      3.1,
      0.44
    ),
    cobicistat = c(
      776,
      4.36,
      2,
      6.58,
      3.62,
      0.589,
      0.025,
      0,
      7.61,
      0.3,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      3.5,
      1,
      0,
      1,
      0,
      1,
      33.235,
      0,
      1,
      0,
      0,
      1,
      1.558,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0.93,
      4.59,
      0,
      0,
      0,
      0,
      0,
      0,
      26.4,
      0.175,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    ),
    efavirenz = c(
      315.7,
      4.6,
      3,
      10.2,
      0,
      0.74,
      0.02,
      0,
      2.5,
      0.25,
      1,
      1,
      1,
      1,
      3,
      5,
      5,
      1,
      1,
      0,
      1,
      0,
      1,
      0.002,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0.07,
      0.08,
      0.55,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      20.6,
      0,
      0,
      4.8,
      0,
      0,
      0,
      1,
      0,
      1,
      14.5,
      3.9,
      5.7,
      0.8,
      6.5,
      3.9
    ),
    etravirine = c(
      435.3,
      5.2,
      1,
      3.75,
      0,
      0.7,
      0.02,
      0,
      0.1963,
      0.13,
      1,
      0,
      0,
      1,
      4,
      1,
      1,
      1,
      1,
      0,
      1,
      0,
      1,
      2.321,
      0,
      1,
      4.306,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      1,
      12.5,
      0.517,
      0,
      1,
      2.5,
      0.517
    ),
    clarithromycin = c(
      748,
      3.16,
      1,
      8.99,
      0,
      0.64,
      0.43,
      0,
      1.23,
      1,
      1,
      1,
      0,
      1,
      0.7,
      1,
      1,
      5,
      1,
      0,
      1,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      11.83,
      7.8,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      2.3,
      0.027,
      0,
      1,
      0,
      1,
      0,
      1,
      0,
      1
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Victim drug: bictegravir (Drug/DrugLibrary/bictegravir.m of the deposited Matlab code).
    # Fixed PBPK inputs, not estimates. Replace the whole set to simulate
    # another drug of the library (see the drugLibrary metadata).
    mw <- fixed(449.39)
    label("Victim drug: molecular weight (g/mol)") # bictegravir.m line 12: DRUG.MolW = 449.39
    logp <- fixed(1.28)
    label("Victim drug: octanol:water partition coefficient logP (unitless)") # bictegravir.m line 15: DRUG.logP = 1.28
    dtype <- fixed(3)
    label("Victim drug: compound class (1 monoprotic base; 2 diprotic base; 3 monoprotic acid; 4 diprotic acid; 5 neutral; 6 zwitterion) (code)") # bictegravir.m line 20: DRUG.type = mono_acid
    pka1 <- fixed(9.81)
    label("Victim drug: pKa 1 (unitless)") # bictegravir.m line 25: DRUG.pka1 = 9.81
    pka2 <- fixed(0)
    label("Victim drug: pKa 2 (0 when not used) (unitless)") # bictegravir.m line 26: DRUG.pka2 = 0.0
    bp <- fixed(0.64)
    label("Victim drug: blood:plasma concentration ratio (unitless)") # bictegravir.m line 29: DRUG.BP = 0.64
    fu <- fixed(0.0025)
    label("Victim drug: fraction unbound in plasma at the reference binding-protein concentration (unitless)") # bictegravir.m line 32: DRUG.fu = 0.0025
    pb_aag <- fixed(0)
    label("Victim drug: main plasma binding protein (0 albumin; 1 alpha1-acid glycoprotein) (flag)") # bictegravir.m line 34: DRUG.protein = albumin
    papp <- fixed(24.6)
    label("Victim drug: Caco-2 apparent permeability (1e-6 cm/s)") # bictegravir.m line 61: DRUG.Papp = 24.6
    kperup <- fixed(0)
    label("Victim drug: uptake-layer to enterocyte permeation rate (0 means 1000 1/h) (1/h)") # bictegravir.m line 64: DRUG.kPerUP = 0.0
    fgp <- fixed(1)
    label("Victim drug: enterocyte to gut-tissue permeability scalar (unitless)") # bictegravir.m line 70: DRUG.GUPermScalar = 1.0
    fabscolon <- fixed(1)
    label("Victim drug: colonic absorption switch (1 on; 0 off) (flag)") # bictegravir.m line 67: DRUG.FabsColon = 1.0
    lagtime <- fixed(0)
    label("Victim drug: dosing lag: an oral dose enters the stomach this long after its nominal time (h)") # bictegravir.m: not set (framework default 0)
    kpscalar <- fixed(1)
    label("Victim drug: global partition-coefficient scalar (unitless)") # bictegravir.m line 78: DRUG.KpScalarAll = 1.0
    jin_all <- fixed(0.7)
    label("Victim drug: cellular influx scalar for all tissues (unitless)") # bictegravir.m line 82: DRUG.JinScalarAll = 0.7
    jin_adipose <- fixed(1)
    label("Victim drug: adipose-specific cellular influx scalar (unitless)") # bictegravir.m: not set (framework default 1)
    jin_muscle <- fixed(1)
    label("Victim drug: muscle-specific cellular influx scalar (unitless)") # bictegravir.m: not set (framework default 1)
    jin_liver <- fixed(2)
    label("Victim drug: liver-specific cellular influx scalar (unitless)") # bictegravir.m line 83: DRUG.JinScalar(liver) = 2.0
    fin_all <- fixed(1)
    label("Victim drug: global vascular-to-interstitial scalar (unitless)") # bictegravir.m line 86: DRUG.FinScalarAll = 1.0
    vmax1_cyp3a4 <- fixed(0)
    label("Victim drug: CYP3A4 Vmax (pathway 1) (pmol/min/pmol)") # bictegravir.m: not set (framework default 0)
    km1_cyp3a4 <- fixed(1)
    label("Victim drug: CYP3A4 Km (pathway 1) (uM)") # bictegravir.m: not set (framework default 1)
    vmax2_cyp3a4 <- fixed(0)
    label("Victim drug: CYP3A4 Vmax (pathway 2) (pmol/min/pmol)") # bictegravir.m: not set (framework default 0)
    km2_cyp3a4 <- fixed(1)
    label("Victim drug: CYP3A4 Km (pathway 2) (uM)") # bictegravir.m: not set (framework default 1)
    clint_cyp3a4 <- fixed(0.114)
    label("Victim drug: CYP3A4 intrinsic clearance (uL/min/pmol)") # bictegravir.m line 94: DRUG.CLint_CYP_1(CYP3A4) = 0.114
    vmax1_cyp2c19 <- fixed(0)
    label("Victim drug: CYP2C19 Vmax (pmol/min/pmol)") # bictegravir.m: not set (framework default 0)
    km1_cyp2c19 <- fixed(1)
    label("Victim drug: CYP2C19 Km (uM)") # bictegravir.m: not set (framework default 1)
    clint_cyp2c19 <- fixed(0)
    label("Victim drug: CYP2C19 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    vmax1_cyp2d6 <- fixed(0)
    label("Victim drug: CYP2D6 Vmax (pmol/min/pmol)") # bictegravir.m: not set (framework default 0)
    km1_cyp2d6 <- fixed(1)
    label("Victim drug: CYP2D6 Km (uM)") # bictegravir.m: not set (framework default 1)
    clint_cyp2d6 <- fixed(0)
    label("Victim drug: CYP2D6 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_cyp2c8 <- fixed(0)
    label("Victim drug: CYP2C8 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_cyp1a2 <- fixed(0)
    label("Victim drug: CYP1A2 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_cyp2a6 <- fixed(0)
    label("Victim drug: CYP2A6 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_cyp2b6 <- fixed(0)
    label("Victim drug: CYP2B6 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_cyp2j2 <- fixed(0)
    label("Victim drug: CYP2J2 intrinsic clearance (uL/min/pmol)") # bictegravir.m: not set (framework default 0)
    clint_ugt1a1 <- fixed(0.292)
    label("Victim drug: UGT1A1 intrinsic clearance (uL/min/pmol)") # bictegravir.m line 95: DRUG.CLint_UGT_1(UGT1A1) = 0.292
    vmax_ugt1a4 <- fixed(0)
    label("Victim drug: UGT1A4 Vmax (pmol/min/pmol)") # bictegravir.m: not set (framework default 0)
    km_ugt1a4 <- fixed(1)
    label("Victim drug: UGT1A4 Km (uM)") # bictegravir.m: not set (framework default 1)
    clint_hep <- fixed(3.993)
    label("Victim drug: hepatic intrinsic clearance not assigned to an enzyme (uL/min/mg)") # bictegravir.m line 98: DRUG.CLint = 3.993
    clrenal <- fixed(0.0043)
    label("Victim drug: renal clearance at GFR 130 (men) or 120 (women) mL/min (L/h)") # bictegravir.m line 101: DRUG.CLrenal = 0.0043
    clbile <- fixed(0)
    label("Victim drug: biliary clearance of unbound intracellular liver drug (L/h)") # bictegravir.m line 105: DRUG.CLbile = 0.0
    ki_cyp3a4 <- fixed(0)
    label("Victim drug: CYP3A4 competitive inhibition constant (0 none) (uM)") # bictegravir.m line 134: DRUG.Ki_CYP(CYP3A4) = 0.0
    ki_cyp2c19 <- fixed(0)
    label("Victim drug: CYP2C19 competitive inhibition constant (0 none) (uM)") # bictegravir.m: not set (framework default 0)
    ki_cyp2d6 <- fixed(0)
    label("Victim drug: CYP2D6 competitive inhibition constant (0 none) (uM)") # bictegravir.m: not set (framework default 0)
    ki_cyp2c8 <- fixed(0)
    label("Victim drug: CYP2C8 competitive inhibition constant (0 none) (uM)") # bictegravir.m: not set (framework default 0)
    ki_cyp1a2 <- fixed(0)
    label("Victim drug: CYP1A2 competitive inhibition constant (0 none) (uM)") # bictegravir.m: not set (framework default 0)
    ki_ugt1a1 <- fixed(0)
    label("Victim drug: UGT1A1 competitive inhibition constant (0 none) (uM)") # bictegravir.m line 135: DRUG.Ki_UGT(UGT1A1) = 0.0
    kinact_cyp3a4 <- fixed(0)
    label("Victim drug: CYP3A4 maximal inactivation rate (mechanism-based) (1/h)") # bictegravir.m line 142: DRUG.kinact_CYP(CYP3A4) = 0.0
    kapp_cyp3a4 <- fixed(1)
    label("Victim drug: CYP3A4 half-maximal inactivation concentration (uM)") # bictegravir.m line 143: DRUG.Kapp_CYP(CYP3A4) = 0.0 (0 replaced by 1 as run)
    kinact_cyp2j2 <- fixed(0)
    label("Victim drug: CYP2J2 maximal inactivation rate (mechanism-based) (1/h)") # bictegravir.m: not set (framework default 0)
    kapp_cyp2j2 <- fixed(1)
    label("Victim drug: CYP2J2 half-maximal inactivation concentration (uM)") # bictegravir.m: not set (framework default 1)
    indmax_cyp3a4 <- fixed(0)
    label("Victim drug: CYP3A4 maximal induction (fold)") # bictegravir.m line 148: DRUG.IndMax_CYP(CYP3A4) = 0.0
    ic50_cyp3a4 <- fixed(1)
    label("Victim drug: CYP3A4 half-maximal induction concentration (uM)") # bictegravir.m line 149: DRUG.IC50_CYP(CYP3A4) = 0.0 (0 replaced by 1 as run)
    indmax_cyp2b6 <- fixed(0)
    label("Victim drug: CYP2B6 maximal induction (fold)") # bictegravir.m: not set (framework default 0)
    ic50_cyp2b6 <- fixed(1)
    label("Victim drug: CYP2B6 half-maximal induction concentration (uM)") # bictegravir.m: not set (framework default 1)
    indmax_ugt1a1 <- fixed(0)
    label("Victim drug: UGT1A1 maximal induction (fold)") # bictegravir.m: not set (framework default 0)
    ic50_ugt1a1 <- fixed(1)
    label("Victim drug: UGT1A1 half-maximal induction concentration (uM)") # bictegravir.m: not set (framework default 1)
    # -------------------------------------------------------------------
    # Perpetrator drug: voriconazole (Drug/DrugLibrary/voriconazole.m of the deposited Matlab code).
    # Fixed PBPK inputs, not estimates. Replace the whole set to simulate
    # another drug of the library (see the drugLibrary metadata).
    mw_perpetrator <- fixed(349.31)
    label("Perpetrator drug: molecular weight (g/mol)") # voriconazole.m line 12: DRUG.MolW = 349.31; Table S1 'MW' 349.3
    logp_perpetrator <- fixed(1.8)
    label("Perpetrator drug: octanol:water partition coefficient logP (unitless)") # voriconazole.m line 15: DRUG.logP = 1.8; Table S1 'logP' 1.8
    dtype_perpetrator <- fixed(1)
    label("Perpetrator drug: compound class (1 monoprotic base; 2 diprotic base; 3 monoprotic acid; 4 diprotic acid; 5 neutral; 6 zwitterion) (code)") # voriconazole.m line 20: DRUG.type = mono_base; Table S1 'drug type' mb
    pka1_perpetrator <- fixed(1.76)
    label("Perpetrator drug: pKa 1 (unitless)") # voriconazole.m line 25: DRUG.pka1 = 1.76; Table S1 'pKa 1' 1.76
    pka2_perpetrator <- fixed(0)
    label("Perpetrator drug: pKa 2 (0 when not used) (unitless)") # voriconazole.m line 26: DRUG.pka2 = 0.0
    bp_perpetrator <- fixed(1.23)
    label("Perpetrator drug: blood:plasma concentration ratio (unitless)") # voriconazole.m line 29: DRUG.BP = 1.23; Table S1 'BP' 1.23
    fu_perpetrator <- fixed(0.42)
    label("Perpetrator drug: fraction unbound in plasma at the reference binding-protein concentration (unitless)") # voriconazole.m line 32: DRUG.fu = 0.42; Table S1 'fup' 0.42
    pb_aag_perpetrator <- fixed(0)
    label("Perpetrator drug: main plasma binding protein (0 albumin; 1 alpha1-acid glycoprotein) (flag)") # voriconazole.m line 34: DRUG.protein = albumin; Table S1 'binding protein' HSA
    papp_perpetrator <- fixed(28.1)
    label("Perpetrator drug: Caco-2 apparent permeability (1e-6 cm/s)") # voriconazole.m line 61: DRUG.Papp = 28.1; Table S1 'Papp' 28.1
    kperup_perpetrator <- fixed(0)
    label("Perpetrator drug: uptake-layer to enterocyte permeation rate (0 means 1000 1/h) (1/h)") # voriconazole.m line 64: DRUG.kPerUP = 0.0
    fgp_perpetrator <- fixed(1)
    label("Perpetrator drug: enterocyte to gut-tissue permeability scalar (unitless)") # voriconazole.m line 73: DRUG.GUPermScalar = 1.0
    fabscolon_perpetrator <- fixed(1)
    label("Perpetrator drug: colonic absorption switch (1 on; 0 off) (flag)") # voriconazole.m line 70: DRUG.FabsColon = 1.0
    lagtime_perpetrator <- fixed(0)
    label("Perpetrator drug: dosing lag: an oral dose enters the stomach this long after its nominal time (h)") # voriconazole.m line 67: DRUG.LagTime = 0.0
    kpscalar_perpetrator <- fixed(2.75)
    label("Perpetrator drug: global partition-coefficient scalar (unitless)") # voriconazole.m line 81: DRUG.KpScalarAll = 2.75
    jin_all_perpetrator <- fixed(1)
    label("Perpetrator drug: cellular influx scalar for all tissues (unitless)") # voriconazole.m line 85: DRUG.JinScalarAll = 1.0
    jin_adipose_perpetrator <- fixed(1)
    label("Perpetrator drug: adipose-specific cellular influx scalar (unitless)") # voriconazole.m line 86: DRUG.JinScalar(adipose) = 1.0
    jin_muscle_perpetrator <- fixed(1)
    label("Perpetrator drug: muscle-specific cellular influx scalar (unitless)") # voriconazole.m line 87: DRUG.JinScalar(muscle) = 1.0
    jin_liver_perpetrator <- fixed(1)
    label("Perpetrator drug: liver-specific cellular influx scalar (unitless)") # voriconazole.m line 88: DRUG.JinScalar(liver) = 1.0
    fin_all_perpetrator <- fixed(1)
    label("Perpetrator drug: global vascular-to-interstitial scalar (unitless)") # voriconazole.m line 91: DRUG.FinScalarAll = 1.0
    vmax1_cyp3a4_perpetrator <- fixed(1212)
    label("Perpetrator drug: CYP3A4 Vmax (pathway 1) (pmol/min/pmol)") # voriconazole.m line 100: DRUG.Vmax_CYP_1(CYP3A4) = 1212; Table S1 'CYP3A4 Vmax' 1212
    km1_cyp3a4_perpetrator <- fixed(15)
    label("Perpetrator drug: CYP3A4 Km (pathway 1) (uM)") # voriconazole.m line 101: DRUG.KM_CYP_1(CYP3A4) = 15; Table S1 'CYP3A4 KM' 15
    vmax2_cyp3a4_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP3A4 Vmax (pathway 2) (pmol/min/pmol)") # voriconazole.m: not set (framework default 0)
    km2_cyp3a4_perpetrator <- fixed(1)
    label("Perpetrator drug: CYP3A4 Km (pathway 2) (uM)") # voriconazole.m: not set (framework default 1)
    clint_cyp3a4_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP3A4 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    vmax1_cyp2c19_perpetrator <- fixed(4.19)
    label("Perpetrator drug: CYP2C19 Vmax (pmol/min/pmol)") # voriconazole.m line 103: DRUG.Vmax_CYP_1(CYP2C19) = 4.19; Table S1 'CYP2C19 Vmax' 4.19
    km1_cyp2c19_perpetrator <- fixed(3.5)
    label("Perpetrator drug: CYP2C19 Km (uM)") # voriconazole.m line 104: DRUG.KM_CYP_1(CYP2C19) = 3.5; Table S1 'CYP2C19 KM' 3.5
    clint_cyp2c19_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2C19 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    vmax1_cyp2d6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2D6 Vmax (pmol/min/pmol)") # voriconazole.m: not set (framework default 0)
    km1_cyp2d6_perpetrator <- fixed(1)
    label("Perpetrator drug: CYP2D6 Km (uM)") # voriconazole.m: not set (framework default 1)
    clint_cyp2d6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2D6 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_cyp2c8_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2C8 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_cyp1a2_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP1A2 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_cyp2a6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2A6 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_cyp2b6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2B6 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_cyp2j2_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2J2 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    clint_ugt1a1_perpetrator <- fixed(0)
    label("Perpetrator drug: UGT1A1 intrinsic clearance (uL/min/pmol)") # voriconazole.m: not set (framework default 0)
    vmax_ugt1a4_perpetrator <- fixed(0)
    label("Perpetrator drug: UGT1A4 Vmax (pmol/min/pmol)") # voriconazole.m: not set (framework default 0)
    km_ugt1a4_perpetrator <- fixed(1)
    label("Perpetrator drug: UGT1A4 Km (uM)") # voriconazole.m: not set (framework default 1)
    clint_hep_perpetrator <- fixed(0.346)
    label("Perpetrator drug: hepatic intrinsic clearance not assigned to an enzyme (uL/min/mg)") # voriconazole.m line 110: DRUG.CLint = 0.346; Table S1 'Unspecified' 0.346
    clrenal_perpetrator <- fixed(0.096)
    label("Perpetrator drug: renal clearance at GFR 130 (men) or 120 (women) mL/min (L/h)") # voriconazole.m line 113: DRUG.CLrenal = 0.096; Table S1 'CLrenal' 0.096
    clbile_perpetrator <- fixed(0)
    label("Perpetrator drug: biliary clearance of unbound intracellular liver drug (L/h)") # voriconazole.m line 117: DRUG.CLbile = 0.0
    ki_cyp3a4_perpetrator <- fixed(0.66)
    label("Perpetrator drug: CYP3A4 competitive inhibition constant (0 none) (uM)") # voriconazole.m line 145: DRUG.Ki_CYP(CYP3A4) = 0.66; Table S1 'CYP3A4 Ki' 0.66
    ki_cyp2c19_perpetrator <- fixed(5.1)
    label("Perpetrator drug: CYP2C19 competitive inhibition constant (0 none) (uM)") # voriconazole.m line 146: DRUG.Ki_CYP(CYP2C19) = 5.1; Table S1 'CYP2C19 Ki' 5.1
    ki_cyp2d6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2D6 competitive inhibition constant (0 none) (uM)") # voriconazole.m: not set (framework default 0)
    ki_cyp2c8_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2C8 competitive inhibition constant (0 none) (uM)") # voriconazole.m: not set (framework default 0)
    ki_cyp1a2_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP1A2 competitive inhibition constant (0 none) (uM)") # voriconazole.m: not set (framework default 0)
    ki_ugt1a1_perpetrator <- fixed(0)
    label("Perpetrator drug: UGT1A1 competitive inhibition constant (0 none) (uM)") # voriconazole.m: not set (framework default 0)
    kinact_cyp3a4_perpetrator <- fixed(9.33)
    label("Perpetrator drug: CYP3A4 maximal inactivation rate (mechanism-based) (1/h)") # voriconazole.m line 151: DRUG.kinact_CYP(CYP3A4) = 9.33; Table S1 'CYP3A4 kinact' 9.33
    kapp_cyp3a4_perpetrator <- fixed(0.015)
    label("Perpetrator drug: CYP3A4 half-maximal inactivation concentration (uM)") # voriconazole.m line 152: DRUG.Kapp_CYP(CYP3A4) = 0.015; Table S1 'CYP3A4 Kapp' 0.015
    kinact_cyp2j2_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2J2 maximal inactivation rate (mechanism-based) (1/h)") # voriconazole.m: not set (framework default 0)
    kapp_cyp2j2_perpetrator <- fixed(1)
    label("Perpetrator drug: CYP2J2 half-maximal inactivation concentration (uM)") # voriconazole.m: not set (framework default 1)
    indmax_cyp3a4_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP3A4 maximal induction (fold)") # voriconazole.m line 157: DRUG.IndMax_CYP(CYP3A4) = 0.0
    ic50_cyp3a4_perpetrator <- fixed(1)
    label("Perpetrator drug: CYP3A4 half-maximal induction concentration (uM)") # voriconazole.m line 158: DRUG.IC50_CYP(CYP3A4) = 0.0 (0 replaced by 1 as run)
    indmax_cyp2b6_perpetrator <- fixed(0)
    label("Perpetrator drug: CYP2B6 maximal induction (fold)") # voriconazole.m: not set (framework default 0)
    ic50_cyp2b6_perpetrator <- fixed(1)
    label("Perpetrator drug: CYP2B6 half-maximal induction concentration (uM)") # voriconazole.m: not set (framework default 1)
    indmax_ugt1a1_perpetrator <- fixed(0)
    label("Perpetrator drug: UGT1A1 maximal induction (fold)") # voriconazole.m: not set (framework default 0)
    ic50_ugt1a1_perpetrator <- fixed(1)
    label("Perpetrator drug: UGT1A1 half-maximal induction concentration (uM)") # voriconazole.m: not set (framework default 1)
    # -------------------------------------------------------------------
    # Second perpetrator drug: rifampicin (Drug/DrugLibrary/rifampicin.m of the deposited Matlab code).
    # Fixed PBPK inputs, not estimates. Replace the whole set to simulate
    # another drug of the library (see the drugLibrary metadata).
    mw_perpetrator2 <- fixed(823)
    label("Second perpetrator drug: molecular weight (g/mol)") # rifampicin.m line 14: DRUG.MolW = 823.0
    logp_perpetrator2 <- fixed(3.28)
    label("Second perpetrator drug: octanol:water partition coefficient logP (unitless)") # rifampicin.m line 17: DRUG.logP = 3.28
    dtype_perpetrator2 <- fixed(6)
    label("Second perpetrator drug: compound class (1 monoprotic base; 2 diprotic base; 3 monoprotic acid; 4 diprotic acid; 5 neutral; 6 zwitterion) (code)") # rifampicin.m line 22: DRUG.type = zwitterion
    pka1_perpetrator2 <- fixed(1.7)
    label("Second perpetrator drug: pKa 1 (unitless)") # rifampicin.m line 27: DRUG.pka1 = 1.7
    pka2_perpetrator2 <- fixed(7.9)
    label("Second perpetrator drug: pKa 2 (0 when not used) (unitless)") # rifampicin.m line 28: DRUG.pka2 = 7.9
    bp_perpetrator2 <- fixed(0.9)
    label("Second perpetrator drug: blood:plasma concentration ratio (unitless)") # rifampicin.m line 31: DRUG.BP = 0.90
    fu_perpetrator2 <- fixed(0.15)
    label("Second perpetrator drug: fraction unbound in plasma at the reference binding-protein concentration (unitless)") # rifampicin.m line 34: DRUG.fu = 0.15
    pb_aag_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: main plasma binding protein (0 albumin; 1 alpha1-acid glycoprotein) (flag)") # rifampicin.m line 36: DRUG.protein = albumin
    papp_perpetrator2 <- fixed(1.472)
    label("Second perpetrator drug: Caco-2 apparent permeability (1e-6 cm/s)") # rifampicin.m line 63: DRUG.Papp = 1.472
    kperup_perpetrator2 <- fixed(3)
    label("Second perpetrator drug: uptake-layer to enterocyte permeation rate (0 means 1000 1/h) (1/h)") # rifampicin.m line 66: DRUG.kPerUP = 3.0
    fgp_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: enterocyte to gut-tissue permeability scalar (unitless)") # rifampicin.m line 75: DRUG.GUPermScalar = 1.0
    fabscolon_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: colonic absorption switch (1 on; 0 off) (flag)") # rifampicin.m line 72: DRUG.FabsColon = 1.0
    lagtime_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: dosing lag: an oral dose enters the stomach this long after its nominal time (h)") # rifampicin.m line 69: DRUG.LagTime = 0.0
    kpscalar_perpetrator2 <- fixed(0.22)
    label("Second perpetrator drug: global partition-coefficient scalar (unitless)") # rifampicin.m line 83: DRUG.KpScalarAll = 0.22
    jin_all_perpetrator2 <- fixed(0.05)
    label("Second perpetrator drug: cellular influx scalar for all tissues (unitless)") # rifampicin.m line 87: DRUG.JinScalarAll = 0.05
    jin_adipose_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: adipose-specific cellular influx scalar (unitless)") # rifampicin.m: not set (framework default 1)
    jin_muscle_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: muscle-specific cellular influx scalar (unitless)") # rifampicin.m: not set (framework default 1)
    jin_liver_perpetrator2 <- fixed(10)
    label("Second perpetrator drug: liver-specific cellular influx scalar (unitless)") # rifampicin.m line 88: DRUG.JinScalar(liver) = 1.0/0.1
    fin_all_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: global vascular-to-interstitial scalar (unitless)") # rifampicin.m line 91: DRUG.FinScalarAll = 1.0
    vmax1_cyp3a4_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP3A4 Vmax (pathway 1) (pmol/min/pmol)") # rifampicin.m: not set (framework default 0)
    km1_cyp3a4_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP3A4 Km (pathway 1) (uM)") # rifampicin.m: not set (framework default 1)
    vmax2_cyp3a4_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP3A4 Vmax (pathway 2) (pmol/min/pmol)") # rifampicin.m: not set (framework default 0)
    km2_cyp3a4_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP3A4 Km (pathway 2) (uM)") # rifampicin.m: not set (framework default 1)
    clint_cyp3a4_perpetrator2 <- fixed(0.0036)
    label("Second perpetrator drug: CYP3A4 intrinsic clearance (uL/min/pmol)") # rifampicin.m line 99: DRUG.CLint_CYP_1(CYP3A4) = 0.0036
    vmax1_cyp2c19_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2C19 Vmax (pmol/min/pmol)") # rifampicin.m: not set (framework default 0)
    km1_cyp2c19_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP2C19 Km (uM)") # rifampicin.m: not set (framework default 1)
    clint_cyp2c19_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2C19 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    vmax1_cyp2d6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2D6 Vmax (pmol/min/pmol)") # rifampicin.m: not set (framework default 0)
    km1_cyp2d6_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP2D6 Km (uM)") # rifampicin.m: not set (framework default 1)
    clint_cyp2d6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2D6 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_cyp2c8_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2C8 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_cyp1a2_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP1A2 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_cyp2a6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2A6 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_cyp2b6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2B6 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_cyp2j2_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2J2 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    clint_ugt1a1_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: UGT1A1 intrinsic clearance (uL/min/pmol)") # rifampicin.m: not set (framework default 0)
    vmax_ugt1a4_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: UGT1A4 Vmax (pmol/min/pmol)") # rifampicin.m: not set (framework default 0)
    km_ugt1a4_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: UGT1A4 Km (uM)") # rifampicin.m: not set (framework default 1)
    clint_hep_perpetrator2 <- fixed(6.55)
    label("Second perpetrator drug: hepatic intrinsic clearance not assigned to an enzyme (uL/min/mg)") # rifampicin.m line 105: DRUG.CLint = 6.55
    clrenal_perpetrator2 <- fixed(1.2)
    label("Second perpetrator drug: renal clearance at GFR 130 (men) or 120 (women) mL/min (L/h)") # rifampicin.m line 108: DRUG.CLrenal = 1.2
    clbile_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: biliary clearance of unbound intracellular liver drug (L/h)") # rifampicin.m line 112: DRUG.CLbile = 0.0
    ki_cyp3a4_perpetrator2 <- fixed(10.5)
    label("Second perpetrator drug: CYP3A4 competitive inhibition constant (0 none) (uM)") # rifampicin.m line 140: DRUG.Ki_CYP(CYP3A4) = 10.5
    ki_cyp2c19_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2C19 competitive inhibition constant (0 none) (uM)") # rifampicin.m: not set (framework default 0)
    ki_cyp2d6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2D6 competitive inhibition constant (0 none) (uM)") # rifampicin.m: not set (framework default 0)
    ki_cyp2c8_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2C8 competitive inhibition constant (0 none) (uM)") # rifampicin.m: not set (framework default 0)
    ki_cyp1a2_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP1A2 competitive inhibition constant (0 none) (uM)") # rifampicin.m: not set (framework default 0)
    ki_ugt1a1_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: UGT1A1 competitive inhibition constant (0 none) (uM)") # rifampicin.m line 141: DRUG.Ki_UGT(UGT1A1) = 0.0
    kinact_cyp3a4_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP3A4 maximal inactivation rate (mechanism-based) (1/h)") # rifampicin.m line 148: DRUG.kinact_CYP(CYP3A4) = 0.0
    kapp_cyp3a4_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP3A4 half-maximal inactivation concentration (uM)") # rifampicin.m line 149: DRUG.Kapp_CYP(CYP3A4) = 0.0 (0 replaced by 1 as run)
    kinact_cyp2j2_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2J2 maximal inactivation rate (mechanism-based) (1/h)") # rifampicin.m: not set (framework default 0)
    kapp_cyp2j2_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP2J2 half-maximal inactivation concentration (uM)") # rifampicin.m: not set (framework default 1)
    indmax_cyp3a4_perpetrator2 <- fixed(16.68)
    label("Second perpetrator drug: CYP3A4 maximal induction (fold)") # rifampicin.m line 154: DRUG.IndMax_CYP(CYP3A4) = 16.68
    ic50_cyp3a4_perpetrator2 <- fixed(0.32)
    label("Second perpetrator drug: CYP3A4 half-maximal induction concentration (uM)") # rifampicin.m line 155: DRUG.IC50_CYP(CYP3A4) = 0.32
    indmax_cyp2b6_perpetrator2 <- fixed(0)
    label("Second perpetrator drug: CYP2B6 maximal induction (fold)") # rifampicin.m: not set (framework default 0)
    ic50_cyp2b6_perpetrator2 <- fixed(1)
    label("Second perpetrator drug: CYP2B6 half-maximal induction concentration (uM)") # rifampicin.m: not set (framework default 1)
    indmax_ugt1a1_perpetrator2 <- fixed(1.668)
    label("Second perpetrator drug: UGT1A1 maximal induction (fold)") # rifampicin.m line 157: DRUG.IndMax_UGT(UGT1A1) = 16.68*10^-1
    ic50_ugt1a1_perpetrator2 <- fixed(0.32)
    label("Second perpetrator drug: UGT1A1 half-maximal induction concentration (uM)") # rifampicin.m line 158: DRUG.IC50_UGT(UGT1A1) = 0.32
    # -------------------------------------------------------------------
    # Scenario switches (structure of a framework run, not drug inputs).
    ddi_on <- fixed(1)
    label("Victim arm: 1 = with the co-administered drugs (framework DDI slot), 0 = victim-alone control slot (unitless)") # PBPK_PreProcessing.m Get_KinPar / Get_DDIPar
    on_perpetrator <- fixed(1)
    label("Perpetrator drug present in the simulated run (1) or not (0) (unitless)") # PBPK_UserChoice.m DRUG.No
    on_perpetrator2 <- fixed(1)
    label("Second perpetrator drug present in the simulated run (1) or not (0) (unitless)") # PBPK_UserChoice.m DRUG.No
    fupref_perpetrator <- fixed(1)
    label("Perpetrator is drug 1 of the run, whose fup enters every KaPR (1) or not (0) (unitless)") # PBPK_Drug_distribution.m DRUG.fup(d); voriconazole (library index 3) precedes rifampicin (5) and bictegravir (6)
    fupref_perpetrator2 <- fixed(0)
    label("Second perpetrator is drug 1 of the run (1) or not (0) (unitless)") # PBPK_Drug_distribution.m DRUG.fup(d)

    # -------------------------------------------------------------------
    # Virtual-population random draws. The deposited generator draws each
    # quantity as normrnd(Mean, (CV/100)*Mean), so the value is
    # Mean*(1 + eta) with eta ~ N(0, (CV/100)^2); the variances below are
    # (CV/100)^2 of the printed CVs. MPPGL_CV is 2.3 instead of 46 when
    # ritonavir is the last drug of the run (PBPK_Population_Liver.m).
    etahct ~ fixed(0.020736) # HCT_CV = 14.4 (PBPK_Population_Tissue.m)
    etahsa ~ fixed(0.006241) # HSA_CV = 7.9 (PBPK_Population_Tissue.m)
    etaw_adipose ~ fixed(0.087616) # WAD_CV = 29.6 (PBPK_Population_Tissue.m)
    etaw_bone ~ fixed(0.017424) # WBO_CV = 13.2 (PBPK_Population_Tissue.m)
    etaw_brain ~ fixed(0.0081) # WBR_CV = 9.0 (PBPK_Population_Tissue.m)
    etaw_gonads ~ fixed(0.121104) # WGO_CV = 34.8 (PBPK_Population_Tissue.m)
    etaw_heart ~ fixed(0.041209) # WHE_CV = 20.3 (PBPK_Population_Tissue.m)
    etaw_kidney ~ fixed(0.045369) # WKI_CV = 21.3 (PBPK_Population_Tissue.m)
    etaw_muscle ~ fixed(0.013924) # WMU_CV = 11.8 (PBPK_Population_Tissue.m)
    etaw_skin ~ fixed(0.006889) # WSK_CV = 8.3 (PBPK_Population_Tissue.m)
    etaw_thymus ~ fixed(0.200704) # WTH_CV = 44.8 (PBPK_Population_Tissue.m)
    etaw_gut ~ fixed(0.005329) # WGU_CV = 7.3 (PBPK_Population_Tissue.m)
    etaw_spleen ~ fixed(0.267289) # WSP_CV = 51.7 (PBPK_Population_Tissue.m)
    etaw_pancreas ~ fixed(0.077284) # WPA_CV = 27.8 (PBPK_Population_Tissue.m)
    etaw_liver ~ fixed(0.056169) # WLI_CV = 23.7 (PBPK_Population_Tissue.m)
    etaw_lnode ~ fixed(0.04) # WLN_CV = 20 (PBPK_Population_Tissue.m)
    etaw_blood ~ fixed(0.010816) # WBL_CV = 10.4 (PBPK_Population_Tissue.m)
    etaco ~ fixed(0.044521) # CO_CV = 21.1 (PBPK_Population_Tissue.m)
    etaltot ~ fixed(0.09) # LT_CV = 30 (PBPK_Population_Tissue.m)
    etagfr ~ fixed(0.021609) # GFR_CV = 14.7 (PBPK_Population_Tissue.m)
    etamppgl ~ fixed(0.2116) # MPPGL_CV = 46 (PBPK_Population_Liver.m; 46 for every drug but ritonavir)
    etamppgl_redraw ~ fixed(0.002116) # Gen_LIPar redraw normrnd(Mean, (CV/1000)*Mean) = 4.6% (PBPK_Population_Liver.m)
    etacyp3a4_liver ~ fixed(0.6561) # hepatic CYP3A4_CV = 81 (PBPK_Population_Liver.m)
    etacyp3a4_liver_redraw ~ fixed(0.006561) # Gen_LIPar redraw at CV/10 = 8.1% (PBPK_Population_Liver.m)
    etaugt1a1_liver ~ fixed(0.36) # UGT1A1_CV = 60 (PBPK_Population_Liver.m)
    etaugt1a1_liver_redraw ~ fixed(0.0036) # Gen_LIPar redraw at CV/10 = 6.0% (PBPK_Population_Liver.m)
    etagastric ~ fixed(1) # GET = 0.25 + (1.00 - 0.25)*rand; a standard normal mapped to U(0,1) by phi() (PBPK_Population_GIT.m)
    etasitt ~ fixed(0.161604) # SIT_CV = 40.2 (PBPK_Population_GIT.m)
    etasitt_redraw ~ fixed(0.00161604) # Gen_GUPar redraw at CV/10 = 4.02% (PBPK_Population_GIT.m)
    etacolt ~ fixed(0.04) # CNT_CV = 20.0 (PBPK_Population_GIT.m)
    etacolt_redraw ~ fixed(0.0004) # Gen_GUPar redraw at CV/10 = 2.0% (PBPK_Population_GIT.m)
    etacyp3a4_gut ~ fixed(0.36) # intestinal CYP3A4_CV = 60 (PBPK_Population_GIT.m)
    etacyp3a4_gut_redraw ~ fixed(0.0036) # Gen_GUPar redraw at CV/10 = 6.0% (PBPK_Population_GIT.m)
    etaaag ~ fixed(0.059049) # AAG_CV = 24.3 (PBPK_Population_Tissue.m)
    etacyp2c19_liver ~ fixed(0.6724) # CYP2C19_CV = 82 (PBPK_Population_Liver.m)
    etacyp2c19_liver_redraw ~ fixed(0.006724) # Gen_LIPar redraw at CV/10 = 8.2% (PBPK_Population_Liver.m)
    etacyp2d6_liver ~ fixed(0.5476) # CYP2D6_CV = 74 (PBPK_Population_Liver.m)
    etacyp2d6_liver_redraw ~ fixed(0.005476) # Gen_LIPar redraw at CV/10 = 7.4% (PBPK_Population_Liver.m)
    etacyp2c8_liver ~ fixed(0.4624) # CYP2C8_CV = 68 (PBPK_Population_Liver.m)
    etacyp2c8_liver_redraw ~ fixed(0.004624) # Gen_LIPar redraw at CV/10 = 6.8% (PBPK_Population_Liver.m)
    etacyp1a2_liver ~ fixed(0.6084) # CYP1A2_CV = 78 (PBPK_Population_Liver.m)
    etacyp1a2_liver_redraw ~ fixed(0.006084) # Gen_LIPar redraw at CV/10 = 7.8% (PBPK_Population_Liver.m)
    etacyp2a6_liver ~ fixed(0.7396) # CYP2A6_CV = 86 (PBPK_Population_Liver.m)
    etacyp2a6_liver_redraw ~ fixed(0.007396) # Gen_LIPar redraw at CV/10 = 8.6% (PBPK_Population_Liver.m)
    etacyp2b6_liver ~ fixed(1.5625) # CYP2B6_CV = 125 (PBPK_Population_Liver.m)
    etacyp2b6_liver_redraw ~ fixed(0.015625) # Gen_LIPar redraw at CV/10 = 12.5% (PBPK_Population_Liver.m)
    etacyp2j2_liver ~ fixed(0.3364) # CYP2J2_CV = 58 (PBPK_Population_Liver.m)
    etacyp2j2_liver_redraw ~ fixed(0.003364) # Gen_LIPar redraw at CV/10 = 5.8% (PBPK_Population_Liver.m)
    etaugt1a4_liver ~ fixed(0.1369) # UGT1A4_CV = 37 (PBPK_Population_Liver.m)
    etaugt1a4_liver_redraw ~ fixed(0.001369) # Gen_LIPar redraw at CV/10 = 3.7% (PBPK_Population_Liver.m)
    etacyp2c19_gut ~ fixed(0.36) # intestinal CYP2C19_CV = 60 (PBPK_Population_GIT.m)
    etacyp2c19_gut_redraw ~ fixed(0.0036) # Gen_GUPar redraw at CV/10 = 6.0% (PBPK_Population_GIT.m)
    etacyp2d6_gut ~ fixed(0.36) # intestinal CYP2D6_CV = 60 (PBPK_Population_GIT.m)
    etacyp2d6_gut_redraw ~ fixed(0.0036) # Gen_GUPar redraw at CV/10 = 6.0% (PBPK_Population_GIT.m)
  })

  model({
    # Paths below are files of the deposited s002 CPT_Matlab_Code folder.
    # ====================================================================
    # 1. Demographics (Population/PBPK_Population_Demographics.m)
    # ====================================================================
    bsa <- 0.007184 * HT^0.725 * WT^0.425 # BSA, DuBois & DuBois (m^2)

    # ====================================================================
    # 2. Blood (Population/PBPK_Population_Tissue.m, 'Blood parameters').
    # GenTisPar truncates at Min/Max only when length(Min) is 1 or 2,
    # i.e. for populations of one or two subjects; for the 100- and
    # 500-subject populations of the paper the draws are untruncated.
    # ====================================================================
    hct <- (0.443 - 0.033 * SEXF) * (1 + etahct) # HCT_Mean
    hsa <- (-0.0709 * AGE + 47.7) * (1 + etahsa) # HSA_Mean (g/L)
    c_aag <- 0.798 * (1 + etaaag) # AAG_Mean (g/L), AAG_CV = 24.3; untruncated like the other GenTisPar draws

    # ====================================================================
    # 3. Organ weights in kg (PBPK_Population_Tissue.m, 'Organ weights').
    # A negative draw is replaced by the block's alternative (WXX_Alt);
    # a negative mean makes normrnd return NaN, which GenTisPar replaces
    # by the (negative) mean and then by the alternative.
    # ====================================================================
    w_lung <- exp(0.00771 * AGE + 0.0279 * HT - 5.58) # WLU_Mean; WLU_CV = 0

    w_adipose_tv <- 0.68 * WT - 0.56 * HT + 6.1 * SEXF + 65 # WAD_Mean
    w_adipose_alt <- ((0.1875 - 0.0796 * SEXF) * AGE + (14.296 + 17.305 * SEXF)) * (WT / 100) # WAD_Alt
    w_adipose_0 <- w_adipose_tv * (1 + etaw_adipose)
    if (w_adipose_tv < 0 || w_adipose_0 < 0) {
      w_adipose_0 <- w_adipose_alt
    }

    w_bone_tv <- exp(0.024 * HT - 1.9) # WBO_Mean
    w_bone <- w_bone_tv * (1 + etaw_bone)
    if (w_bone < 0) {
      w_bone <- w_bone_tv
    }

    w_brain_tv <- exp(-0.00075 * AGE + 0.00778 * HT - 0.97) # WBR_Mean
    w_brain <- w_brain_tv * (1 + etaw_brain)
    if (w_brain < 0) {
      w_brain <- w_brain_tv
    }

    w_gonads_tv <- -0.00022 * AGE - 0.00034 * WT - 0.030 * SEXF + 0.072 # WGO_Mean
    w_gonads_alt <- (0.00049 - 0.00037 * SEXF) * WT # WGO_Alt
    w_gonads <- w_gonads_tv * (1 + etaw_gonads)
    if (w_gonads_tv < 0 || w_gonads < 0) {
      w_gonads <- w_gonads_alt
    }

    w_heart_tv <- 0.34 * bsa + 0.0018 * AGE - 0.36 # WHE_Mean
    w_heart <- w_heart_tv * (1 + etaw_heart)
    if (w_heart < 0) {
      w_heart <- w_heart_tv
    }

    w_kidney_tv <- -0.00038 * AGE - 0.056 * SEXF + 0.33 # WKI_Mean
    w_kidney <- w_kidney_tv * (1 + etaw_kidney)
    if (w_kidney < 0) {
      w_kidney <- w_kidney_tv
    }

    w_muscle_tv <- 17.9 * bsa - 0.0667 * AGE - 5.68 * SEXF - 1.22 # WMU_Mean
    w_muscle_0 <- w_muscle_tv * (1 + etaw_muscle)
    if (w_muscle_0 < 0) {
      w_muscle_0 <- w_muscle_tv
    }

    w_skin_tv <- exp(-0.0058 * AGE - 0.37 * SEXF + 1.13) # WSK_Mean
    w_skin <- w_skin_tv * (1 + etaw_skin)
    if (w_skin < 0) {
      w_skin <- w_skin_tv
    }

    w_thymus_tv <- 0.0221 # WTH_Mean
    w_thymus <- w_thymus_tv * (1 + etaw_thymus)
    if (w_thymus < 0) {
      w_thymus <- w_thymus_tv
    }

    w_gut_tv <- 3e-06 * HT^2.49 # WGU_Mean
    w_gut <- w_gut_tv * (1 + etaw_gut)
    if (w_gut < 0) {
      w_gut <- w_gut_tv
    }

    w_spleen_tv <- exp(1.13 * bsa - 3.93) # WSP_Mean
    w_spleen <- w_spleen_tv * (1 + etaw_spleen)
    if (w_spleen < 0) {
      w_spleen <- 0.00228 * WT # WSP_Alt
    }

    w_pancreas_tv <- 0.103 # WPA_Mean
    w_pancreas <- w_pancreas_tv * (1 + etaw_pancreas)
    if (w_pancreas < 0) {
      w_pancreas <- w_pancreas_tv
    }

    w_liver_tv <- 0.52 * bsa^1.73 # WLI_Mean
    w_liver <- w_liver_tv * (1 + etaw_liver)
    if (w_liver < 0) {
      w_liver <- w_liver_tv
    }

    w_lnode_tv <- 0.00386 * WT # WLN_Mean (Gill et al. 2016)
    w_lnode <- w_lnode_tv * (1 + etaw_lnode)
    if (w_lnode < 0) {
      w_lnode <- w_lnode_tv
    }

    w_blood <- (-0.0098 * AGE - 1.89 * SEXF + 6.06) * (1 + etaw_blood) # WBL_Mean
    w_plasma <- w_blood * (1 - 0.91 * hct)
    w_rbc <- w_blood - w_plasma

    # Remaining tissue: body weight balance, held at >= 2% of body weight
    # by taking the shortfall 75:25 from adipose and muscle, and adipose
    # held at >= 2% of body weight by taking the shortfall from muscle.
    w_other_0 <- WT - w_blood - (w_lung + w_adipose_0 + w_bone + w_brain + w_gonads +
      w_heart + w_kidney + w_muscle_0 + w_skin + w_thymus + w_gut + w_spleen +
      w_pancreas + w_liver + w_lnode)
    mirema <- 0
    if (w_other_0 < 0.02 * WT) {
      mirema <- 0.02 * WT - w_other_0
    }
    w_adipose_1 <- w_adipose_0 - 0.75 * mirema
    w_muscle_1 <- w_muscle_0 - 0.25 * mirema
    miadma <- 0
    if (w_adipose_1 < 0.02 * WT) {
      miadma <- 0.02 * WT - w_adipose_1
    }
    w_adipose <- w_adipose_1 + miadma
    w_muscle <- w_muscle_1 - miadma
    w_other <- WT - w_blood - (w_lung + w_adipose + w_bone + w_brain + w_gonads +
      w_heart + w_kidney + w_muscle + w_skin + w_thymus + w_gut + w_spleen +
      w_pancreas + w_liver + w_lnode)
    w_all <- w_lung + w_adipose + w_bone + w_brain + w_gonads + w_heart + w_kidney +
      w_muscle + w_skin + w_thymus + w_gut + w_spleen + w_pancreas + w_liver +
      w_lnode + w_other + w_plasma + w_rbc

    # ====================================================================
    # 4. Organ volumes in L (PBPK_Population_Tissue.m, 'Organ density').
    # The remaining-tissue entries of the density and composition tables
    # are initialised to 1 and then replaced by the Worg-weighted mean
    # over all 18 compartments, a mean that includes the remaining tissue
    # itself at that initial value of 1 (as-run).
    # ====================================================================
    den_other <- (w_lung * 1.0 + w_adipose * 0.916 + w_bone * 1.9 + w_brain * 1.04 +
      w_gonads * 1.045 + w_heart * 1.03 + w_kidney * 1.05 + w_muscle * 1.041 +
      w_skin * 1.1 + w_thymus * 1.025 + w_gut * 1.042 + w_spleen * 1.06 +
      w_pancreas * 1.045 + w_liver * 1.08 + w_lnode * 1.0 + w_other * 1 +
      w_plasma * 1.027 + w_rbc * 1.09) / w_all
    v_lung <- w_lung / 1.0
    v_adipose <- w_adipose / 0.916
    v_bone <- w_bone / 1.9
    v_brain <- w_brain / 1.04
    v_gonads <- w_gonads / 1.045
    v_heart <- w_heart / 1.03
    v_kidney <- w_kidney / 1.05
    v_muscle <- w_muscle / 1.041
    v_skin <- w_skin / 1.1
    v_thymus <- w_thymus / 1.025
    v_gut <- w_gut / 1.042
    v_spleen <- w_spleen / 1.06
    v_pancreas <- w_pancreas / 1.045
    v_liver <- w_liver / 1.08
    v_lnode <- w_lnode / 1.0
    v_other <- w_other / den_other
    v_venous <- (2 / 3) * w_blood / 1.06 # Vvein = Wvein/1.06
    v_arterial <- (1 / 3) * w_blood / 1.06 # Vartery = Wartery/1.06

    # ====================================================================
    # 5. Tissue composition of the remaining tissue (PBPK_Population_Tissue.m,
    # 'Tissue composition'): Worg-weighted means of the tabulated FraEW,
    # FraIW, FraNL, FraNP and KpHSA, and of FraVas.
    # ====================================================================
    fe_other <- (w_lung * 0.348 + w_adipose * 0.141 + w_bone * 0.098 + w_brain * 0.092 +
      w_gonads * 0.239 + w_heart * 0.313 + w_kidney * 0.283 + w_muscle * 0.091 +
      w_skin * 0.623 + w_thymus * 0.150 + w_gut * 0.267 + w_spleen * 0.208 +
      w_pancreas * 0.12 + w_liver * 0.165 + w_lnode * 0.208 + w_other * 1 +
      w_plasma * 0.945 + w_rbc * 0) / w_all
    fi_other <- (w_lung * 0.463 + w_adipose * 0.039 + w_bone * 0.341 + w_brain * 0.678 +
      w_gonads * 0.561 + w_heart * 0.445 + w_kidney * 0.5 + w_muscle * 0.669 +
      w_skin * 0.0947 + w_thymus * 0.626 + w_gut * 0.451 + w_spleen * 0.58 +
      w_pancreas * 0.664 + w_liver * 0.586 + w_lnode * 0.58 + w_other * 1 +
      w_plasma * 0 + w_rbc * 0.666) / w_all
    fnl_other <- (w_lung * 0.003 + w_adipose * 0.79 + w_bone * 0.074 + w_brain * 0.051 +
      w_gonads * 0.007 + w_heart * 0.015 + w_kidney * 0.0207 + w_muscle * 0.0238 +
      w_skin * 0.0248 + w_thymus * 0.017 + w_gut * 0.0487 + w_spleen * 0.0201 +
      w_pancreas * 0.041 + w_liver * 0.0348 + w_lnode * 0.0201 + w_other * 1 +
      w_plasma * 0.35 + w_rbc * 0.17) / w_all
    fnp_other <- (w_lung * 0.009 + w_adipose * 0.002 + w_bone * 0.0011 + w_brain * 0.0565 +
      w_gonads * 0.0077 + w_heart * 0.0166 + w_kidney * 0.0162 + w_muscle * 0.0072 +
      w_skin * 0.0111 + w_thymus * 0.0092 + w_gut * 0.0163 + w_spleen * 0.0198 +
      w_pancreas * 0.0093 + w_liver * 0.0252 + w_lnode * 0.0198 + w_other * 1 +
      w_plasma * 0.225 + w_rbc * 0.29) / w_all
    kphsa_other <- (w_lung * 0.212 + w_adipose * 0.021 + w_bone * 0.1 + w_brain * 0.048 +
      w_gonads * 0.048 + w_heart * 0.157 + w_kidney * 0.13 + w_muscle * 0.025 +
      w_skin * 0.277 + w_thymus * 0.075 + w_gut * 0.158 + w_spleen * 0.097 +
      w_pancreas * 0.06 + w_liver * 0.086 + w_lnode * 0.097 + w_other * 1 +
      w_plasma * 1.0 + w_rbc * 0) / w_all
    ap_other <- (w_lung * 0.5 + w_adipose * 0.4 + w_bone * 0.67 + w_brain * 0.4 +
      w_gonads * 1.23 + w_heart * 3.07 + w_kidney * 2.48 + w_muscle * 2.49 +
      w_skin * 1.32 + w_thymus * 2.3 + w_gut * 2.84 + w_spleen * 2.81 +
      w_pancreas * 1.67 + w_liver * 5.09 + w_lnode * 2.81 + w_other * 1 +
      w_plasma * 0.04 + w_rbc * 0.44) / w_all # AP (acidic phospholipids, mg/g)
    fvas_brain <- -0.000545 * AGE + 0.056 # FraVas(brain)
    fvas_other <- (w_lung * 0.185 + w_adipose * 0.031 + w_bone * 0.05 + w_brain * fvas_brain +
      w_gonads * 0.069 + w_heart * 0.042 + w_kidney * 0.07 + w_muscle * 0.027 +
      w_skin * 0.05 + w_thymus * 0.05 + w_gut * 0.05 + w_spleen * 0.05 +
      w_pancreas * 0.05 + w_liver * 0.05 + w_lnode * 0.05 + w_other * 1 +
      w_plasma * 1 + w_rbc * 1) / w_all

    # ====================================================================
    # 6. Sub-compartment volumes (PBPK_Population_Tissue.m, 'Subcompartment
    # volume'): Vvas = FraVas*Vorg, Vint = Vorg*FraEW - Vvas*(1 - HCT),
    # Vcel = Vorg - Vvas - Vint. FraVas and FraEW per organ as tabulated.
    # ====================================================================
    vv_lung <- 0.185 * v_lung
    ve_lung <- v_lung * 0.348 - vv_lung * (1 - hct)
    vi_lung <- v_lung - vv_lung - ve_lung
    vv_adipose <- 0.031 * v_adipose
    ve_adipose <- v_adipose * 0.141 - vv_adipose * (1 - hct)
    vi_adipose <- v_adipose - vv_adipose - ve_adipose
    vv_bone <- 0.05 * v_bone
    ve_bone <- v_bone * 0.098 - vv_bone * (1 - hct)
    vi_bone <- v_bone - vv_bone - ve_bone
    vv_brain <- fvas_brain * v_brain
    ve_brain <- v_brain * 0.092 - vv_brain * (1 - hct)
    vi_brain <- v_brain - vv_brain - ve_brain
    vv_gonads <- 0.069 * v_gonads
    ve_gonads <- v_gonads * 0.239 - vv_gonads * (1 - hct)
    vi_gonads <- v_gonads - vv_gonads - ve_gonads
    vv_heart <- 0.042 * v_heart
    ve_heart <- v_heart * 0.313 - vv_heart * (1 - hct)
    vi_heart <- v_heart - vv_heart - ve_heart
    vv_kidney <- 0.07 * v_kidney
    ve_kidney <- v_kidney * 0.283 - vv_kidney * (1 - hct)
    vi_kidney <- v_kidney - vv_kidney - ve_kidney
    vv_muscle <- 0.027 * v_muscle
    ve_muscle <- v_muscle * 0.091 - vv_muscle * (1 - hct)
    vi_muscle <- v_muscle - vv_muscle - ve_muscle
    vv_skin <- 0.05 * v_skin
    ve_skin <- v_skin * 0.623 - vv_skin * (1 - hct)
    vi_skin <- v_skin - vv_skin - ve_skin
    vv_thymus <- 0.05 * v_thymus
    ve_thymus <- v_thymus * 0.150 - vv_thymus * (1 - hct)
    vi_thymus <- v_thymus - vv_thymus - ve_thymus
    # Gut: Vint uses the tissue-level Vvas (0.05*Vorg); PBPK_Population_GIT.m
    # then overwrites the gut Vvas with the sum of the four segment vascular
    # spaces, 0.05 of the duodenum, jejunum, ileum and colon shares.
    ve_gut <- v_gut * 0.267 - (0.05 * v_gut) * (1 - hct)
    vv_gut <- 0.05 * (0.060 + 0.279 + 0.316 + 0.192) * v_gut
    vv_spleen <- 0.05 * v_spleen
    ve_spleen <- v_spleen * 0.208 - vv_spleen * (1 - hct)
    vi_spleen <- v_spleen - vv_spleen - ve_spleen
    vv_pancreas <- 0.05 * v_pancreas
    ve_pancreas <- v_pancreas * 0.12 - vv_pancreas * (1 - hct)
    vi_pancreas <- v_pancreas - vv_pancreas - ve_pancreas
    vv_liver <- 0.05 * v_liver
    ve_liver <- v_liver * 0.165 - vv_liver * (1 - hct)
    vi_liver <- v_liver - vv_liver - ve_liver
    vv_lnode <- 0.05 * v_lnode
    ve_lnode <- v_lnode * 0.208 - vv_lnode * (1 - hct)
    vi_lnode <- v_lnode - vv_lnode - ve_lnode
    vv_other <- fvas_other * v_other
    ve_other <- v_other * fe_other - vv_other * (1 - hct)
    vi_other <- v_other - vv_other - ve_other

    # ====================================================================
    # 7. Blood flows (PBPK_Population_Tissue.m, 'blood flows'): percent of
    # cardiac output, then L/h. The remaining-tissue flow is the balance of
    # the listed organs; the brain flow is not subtracted from it (as-run).
    # ====================================================================
    co <- (159 * bsa - 1.56 * AGE + 114) * (1 + etaco) # CO_Mean (L/h)
    fq_adipose_0 <- (0.044 + 0.027 * SEXF) * AGE + 2.4 * SEXF + 3.9
    fq_bone <- 5
    fq_brain <- exp(-0.48 * bsa + 0.04 * SEXF + 3.5)
    fq_gonads <- -0.03 * SEXF + 0.05
    fq_heart_0 <- -0.72 * HT - 10 * SEXF + 134
    if (fq_heart_0 < 0 || fq_heart_0 > 12.0) {
      fq_heart_0 <- 0.04 + 0.01 * SEXF # 'values from Valentin (2002)', entered in percent units (as-run)
    }
    fq_kidney_0 <- -8.7 * bsa + 0.29 * HT - 0.081 * AGE - 13
    fq_muscle_0 <- -6.4 * SEXF + 17.5
    fq_skin <- 5
    fq_thymus <- 1.5
    fq_liver <- -0.108 * AGE + 1.04 * SEXF + 27.9
    fq_ha <- 6.5 # hepatic artery
    fq_pv <- fq_liver - fq_ha # portal vein
    fq_gut <- ((2 * SEXF + 14) * fq_pv) / (1.5 * SEXF + 19)
    fq_spleen <- (3 * fq_pv) / (1.5 * SEXF + 19)
    fq_pancreas <- (1 * fq_pv) / (1.5 * SEXF + 19)
    fq_by <- fq_pv - fq_gut - fq_spleen - fq_pancreas # arterial flow bypassing the portal organs
    fq_lnode <- 1.65
    fq_other_0 <- 100 - fq_adipose_0 - fq_bone - fq_gonads - fq_heart_0 - fq_kidney_0 -
      fq_muscle_0 - fq_skin - fq_thymus - fq_liver - fq_lnode
    fqmis <- 0
    if (fq_other_0 < 1) {
      fqmis <- 1 - fq_other_0
    }
    fq_adipose <- fq_adipose_0 - 0.25 * fqmis
    fq_heart <- fq_heart_0 - 0.25 * fqmis
    fq_kidney <- fq_kidney_0 - 0.25 * fqmis
    fq_muscle <- fq_muscle_0 - 0.25 * fqmis
    fq_other <- 100 - fq_adipose - fq_bone - fq_gonads - fq_heart - fq_kidney -
      fq_muscle - fq_skin - fq_thymus - fq_liver - fq_lnode

    q_lung <- co
    q_adipose <- fq_adipose / 100 * co
    q_bone <- fq_bone / 100 * co
    q_brain <- fq_brain / 100 * co
    q_gonads <- fq_gonads / 100 * co
    q_heart <- fq_heart / 100 * co
    q_kidney <- fq_kidney / 100 * co
    q_muscle <- fq_muscle / 100 * co
    q_skin <- fq_skin / 100 * co
    q_thymus <- fq_thymus / 100 * co
    q_gut <- fq_gut / 100 * co
    q_spleen <- fq_spleen / 100 * co
    q_pancreas <- fq_pancreas / 100 * co
    q_liver <- fq_liver / 100 * co
    q_lnode <- fq_lnode / 100 * co
    q_other <- fq_other / 100 * co
    q_ha <- fq_ha / 100 * co
    q_by <- fq_by / 100 * co

    # ====================================================================
    # 8. Lymph flows in L/h (PBPK_Population_Tissue.m, 'lymph flows'; Gill
    # et al. 2016); the remaining tissue takes the unlisted fraction.
    # ====================================================================
    ltot <- 0.00386 * WT * (1 + etaltot) # LT_Mean
    if (ltot < 0) {
      ltot <- 0.00386 * WT
    }
    l_lung <- 0.03 * ltot
    l_adipose <- 0.128 * ltot
    l_bone <- 0 * ltot
    l_brain <- 0.0105 * ltot
    l_gonads <- 0.013 * ltot
    l_heart <- 0.01 * ltot
    l_kidney <- 0.085 * ltot
    l_muscle <- 0.16 * ltot
    l_skin <- 0.073 * ltot
    l_thymus <- 0.011 * ltot
    l_gut <- 0.12 * ltot
    l_spleen <- 0 * ltot
    l_pancreas <- 0.003 * ltot
    l_liver <- 0.33 * ltot
    l_other <- (1 - (0.03 + 0.128 + 0 + 0.0105 + 0.013 + 0.01 + 0.085 + 0.16 + 0.073 +
      0.011 + 0.12 + 0 + 0.003 + 0.33)) * ltot

    # Plasma flow Porg = Qorg*(1 - HCT), QL = Qorg - Lorg, PL = Porg - Lorg
    # (PBPK_ODE_solution.m, 'Blood flows').
    ql_lung <- q_lung - l_lung
    ql_adipose <- q_adipose - l_adipose
    ql_bone <- q_bone - l_bone
    ql_brain <- q_brain - l_brain
    ql_gonads <- q_gonads - l_gonads
    ql_heart <- q_heart - l_heart
    ql_kidney <- q_kidney - l_kidney
    ql_muscle <- q_muscle - l_muscle
    ql_skin <- q_skin - l_skin
    ql_thymus <- q_thymus - l_thymus
    ql_gut <- q_gut - l_gut
    ql_spleen <- q_spleen - l_spleen
    ql_pancreas <- q_pancreas - l_pancreas
    ql_liver <- q_liver - l_liver
    ql_other <- q_other - l_other
    p_lung <- q_lung * (1 - hct)
    p_adipose <- q_adipose * (1 - hct)
    p_bone <- q_bone * (1 - hct)
    p_brain <- q_brain * (1 - hct)
    p_gonads <- q_gonads * (1 - hct)
    p_heart <- q_heart * (1 - hct)
    p_kidney <- q_kidney * (1 - hct)
    p_muscle <- q_muscle * (1 - hct)
    p_skin <- q_skin * (1 - hct)
    p_thymus <- q_thymus * (1 - hct)
    p_gut <- q_gut * (1 - hct)
    p_spleen <- q_spleen * (1 - hct)
    p_pancreas <- q_pancreas * (1 - hct)
    p_liver <- q_liver * (1 - hct)
    p_other <- q_other * (1 - hct)

    # ====================================================================
    # 9. Kidney and liver scalars
    # ====================================================================
    gfr <- exp(-0.0079 * AGE + 0.5 * bsa + 4.2) * (1 + etagfr) # GFR_Mean (mL/min), PBPK_Population_Tissue.m

    # Gen_LIPar (PBPK_Population_Liver.m): a draw outside [Min, Max] is
    # redrawn once from normrnd(Mean, (CV/1000)*Mean).
    mppgl_tv <- 10^(0.0000024 * AGE^3 - 0.00038 * AGE^2 + 0.0158 * AGE + 1.407) # MPPGL_Mean (mg/g liver), Barter 2008
    mppgl <- mppgl_tv * (1 + etamppgl)
    if (mppgl < 10 || mppgl > 110) {
      mppgl <- mppgl_tv * (1 + etamppgl_redraw)
    }
    ab_cyp3a4_liver <- 93.0 * (1 + etacyp3a4_liver) # CYP3A4_Mean (pmol/mg)
    if (ab_cyp3a4_liver < 18.6 || ab_cyp3a4_liver > 601) {
      ab_cyp3a4_liver <- 93.0 * (1 + etacyp3a4_liver_redraw)
    }
    ab_ugt1a1_liver <- 41.0 * (1 + etaugt1a1_liver) # UGT1A1_Mean (pmol/mg)
    if (ab_ugt1a1_liver < 4.0 || ab_ugt1a1_liver > 138) {
      ab_ugt1a1_liver <- 41.0 * (1 + etaugt1a1_liver_redraw)
    }
    # Further hepatic abundances (pmol/mg), same Gen_LIPar draw and redraw.
    # CYP3A5 is omitted: with phenotypes switched off every virtual subject
    # is a CYP3A5 poor metaboliser (PBPK_Population_Demographics.m), so its
    # abundance is 0. CYP2C9 and UGT1A3 metabolise no drug in the library.
    ab_cyp2c19_liver <- 11 * (1 + etacyp2c19_liver) # CYP2C19_Mean, CV 82
    if (ab_cyp2c19_liver < 1 || ab_cyp2c19_liver > 38) {
      ab_cyp2c19_liver <- 11 * (1 + etacyp2c19_liver_redraw)
    }
    ab_cyp2d6_liver <- 12.6 * (1 + etacyp2d6_liver) # CYP2D6_Mean, CV 74
    if (ab_cyp2d6_liver < 4.2 || ab_cyp2d6_liver > 38) {
      ab_cyp2d6_liver <- 12.6 * (1 + etacyp2d6_liver_redraw)
    }
    ab_cyp2c8_liver <- 22.4 * (1 + etacyp2c8_liver) # CYP2C8_Mean, CV 68
    if (ab_cyp2c8_liver < 7.5 || ab_cyp2c8_liver > 67) {
      ab_cyp2c8_liver <- 22.4 * (1 + etacyp2c8_liver_redraw)
    }
    ab_cyp1a2_liver <- 39 * (1 + etacyp1a2_liver) # CYP1A2_Mean, CV 78
    if (ab_cyp1a2_liver < 13 || ab_cyp1a2_liver > 117) {
      ab_cyp1a2_liver <- 39 * (1 + etacyp1a2_liver_redraw)
    }
    ab_cyp2a6_liver <- 27 * (1 + etacyp2a6_liver) # CYP2A6_Mean, CV 86
    if (ab_cyp2a6_liver < 9 || ab_cyp2a6_liver > 81) {
      ab_cyp2a6_liver <- 27 * (1 + etacyp2a6_liver_redraw)
    }
    ab_cyp2b6_liver <- 16 * (1 + etacyp2b6_liver) # CYP2B6_Mean, CV 125
    if (ab_cyp2b6_liver < 2.5 || ab_cyp2b6_liver > 104) {
      ab_cyp2b6_liver <- 16 * (1 + etacyp2b6_liver_redraw)
    }
    ab_cyp2j2_liver <- 1.2 * (1 + etacyp2j2_liver) # CYP2J2_Mean, CV 58
    if (ab_cyp2j2_liver < 0.4 || ab_cyp2j2_liver > 3.6) {
      ab_cyp2j2_liver <- 1.2 * (1 + etacyp2j2_liver_redraw)
    }
    ab_ugt1a4_liver <- 55.4 * (1 + etaugt1a4_liver) # UGT1A4_Mean, CV 37
    if (ab_ugt1a4_liver < 4 || ab_ugt1a4_liver > 106) {
      ab_ugt1a4_liver <- 55.4 * (1 + etaugt1a4_liver_redraw)
    }

    # ====================================================================
    # 10. Gut transit and intestinal CYP3A4 (PBPK_Population_GIT.m).
    # Gen_GUPar redraws out-of-range values like Gen_LIPar.
    # ====================================================================
    t_stomach <- 0.25 + (1.00 - 0.25) * phi(etagastric) # GET ~ U(0.25, 1.00) h
    sitt <- 3.4 * (1 + etasitt) # SIT_Mean (h)
    if (sitt < 0.5 || sitt > 9.5) {
      sitt <- 3.4 * (1 + etasitt_redraw)
    }
    colt <- 17.2 * (1 + etacolt) # CNT_Mean (h)
    if (colt < 17.2 / 5 || colt > 17.2 * 5) {
      colt <- 17.2 * (1 + etacolt_redraw)
    }
    t_duodenum <- sitt * 0.091
    t_jejunum <- sitt * 0.426
    t_ileum <- sitt * 0.483
    t_colon <- colt
    ab_cyp3a4_gut <- 66.2 * (1 + etacyp3a4_gut) # intestinal CYP3A4_Mean (nmol)
    if (ab_cyp3a4_gut < 66.2 / 3 || ab_cyp3a4_gut > 66.2 * 3) {
      ab_cyp3a4_gut <- 66.2 * (1 + etacyp3a4_gut_redraw)
    }
    ab_cyp2c19_gut <- 1.5 * (1 + etacyp2c19_gut) # intestinal CYP2C19_Mean (nmol), CV 60
    if (ab_cyp2c19_gut < 1.5 / 3 || ab_cyp2c19_gut > 1.5 * 3) {
      ab_cyp2c19_gut <- 1.5 * (1 + etacyp2c19_gut_redraw)
    }
    ab_cyp2d6_gut <- 0.8 * (1 + etacyp2d6_gut) # intestinal CYP2D6_Mean (nmol), CV 60
    if (ab_cyp2d6_gut < 0.8 / 3 || ab_cyp2d6_gut > 0.8 * 3) {
      ab_cyp2d6_gut <- 0.8 * (1 + etacyp2d6_gut_redraw)
    }
    # Segment abundances in pmol (nmol * 1e3), Paine 1997 split 0.136 /
    # 0.544 / 0.320 (duodenum / jejunum / ileum).
    abg_cyp3a4_duodenum <- 0.136 * ab_cyp3a4_gut * 1000
    abg_cyp3a4_jejunum <- 0.544 * ab_cyp3a4_gut * 1000
    abg_cyp3a4_ileum <- 0.320 * ab_cyp3a4_gut * 1000
    abg_cyp2c19_duodenum <- 0.136 * ab_cyp2c19_gut * 1000
    abg_cyp2c19_jejunum <- 0.544 * ab_cyp2c19_gut * 1000
    abg_cyp2c19_ileum <- 0.320 * ab_cyp2c19_gut * 1000
    abg_cyp2d6_duodenum <- 0.136 * ab_cyp2d6_gut * 1000
    abg_cyp2d6_jejunum <- 0.544 * ab_cyp2d6_gut * 1000
    abg_cyp2d6_ileum <- 0.320 * ab_cyp2d6_gut * 1000

    # Segment volumes (L): VsegCAT, luminal fluid VfluCAT, segment vascular
    # VvasCAT = 0.05*Vseg, segment interstitial VintCAT (share of the gut
    # Vint), uptake layer VlumCAT = Vseg - Vflu - Vvas - Vint.
    vseg_duodenum <- 0.060 * v_gut
    vseg_jejunum <- 0.279 * v_gut
    vseg_ileum <- 0.316 * v_gut
    vseg_colon <- 0.192 * v_gut
    vflu_duodenum <- 0.136 * vseg_duodenum
    vflu_jejunum <- 0.136 * vseg_jejunum
    vflu_ileum <- 0.136 * vseg_ileum
    vflu_colon <- 0.057 * vseg_colon
    vup_duodenum <- vseg_duodenum - vflu_duodenum - 0.05 * vseg_duodenum - 0.060 * ve_gut
    vup_jejunum <- vseg_jejunum - vflu_jejunum - 0.05 * vseg_jejunum - 0.279 * ve_gut
    vup_ileum <- vseg_ileum - vflu_ileum - 0.05 * vseg_ileum - 0.316 * ve_gut
    vup_colon <- vseg_colon - vflu_colon - 0.05 * vseg_colon - 0.192 * ve_gut

    # Lengths (cm; ICRP), cylinder radius and surface, enterocyte volume
    # VentCAT = 3040*(Surface/0.0001)*770*1e-15 (L), and the absorptive
    # surface SAB = Surface*FPC*FVI.
    len_duodenum <- 0.091 * (1.6 * HT)
    len_jejunum <- 0.426 * (1.6 * HT)
    len_ileum <- 0.483 * (1.6 * HT)
    len_colon <- 0.52 * HT + 18.5
    len_total <- len_duodenum + len_jejunum + len_ileum + len_colon
    surf_duodenum <- 2 * pi * sqrt((vseg_duodenum * 1000) / (pi * len_duodenum)) * len_duodenum
    surf_jejunum <- 2 * pi * sqrt((vseg_jejunum * 1000) / (pi * len_jejunum)) * len_jejunum
    surf_ileum <- 2 * pi * sqrt((vseg_ileum * 1000) / (pi * len_ileum)) * len_ileum
    surf_colon <- 2 * pi * sqrt((vseg_colon * 1000) / (pi * len_colon)) * len_colon
    vent_duodenum <- (3040 * (surf_duodenum / 0.0001)) * 770 * 1e-15
    vent_jejunum <- (3040 * (surf_jejunum / 0.0001)) * 770 * 1e-15
    vent_ileum <- (3040 * (surf_ileum / 0.0001)) * 770 * 1e-15
    vent_colon <- (3040 * (surf_colon / 0.0001)) * 770 * 1e-15

    # ====================================================================
    # 11. Drug blocks. The victim drug takes the bare state names; the two
    # co-administered drugs carry the _perpetrator / _perpetrator2 suffix.
    # ====================================================================
    # --- victim drug: ionisation and plasma binding (PBPK_Drug_distribution.m) ---
    kperup_eff <- kperup
    if (kperup == 0) {
      kperup_eff <- 1000 # DRUG.kPerUP = 0 replaced by 1000 1/h (PBPK_Drug_PostProcessing.m)
    }
    ib1 <- 0
    ia1 <- 0
    ib2 <- 0
    ia2 <- 0
    if (dtype == 1 || dtype == 2 || dtype == 6) {
      ib1 <- 1 # base on pKa1: mono_base, di_base, zwitterion
    }
    if (dtype == 3 || dtype == 4) {
      ia1 <- 1 # acid on pKa1: mono_acid, di_acid
    }
    if (dtype == 2) {
      ib2 <- 1 # base on pKa2: di_base
    }
    if (dtype == 4 || dtype == 6) {
      ia2 <- 1 # acid on pKa2: di_acid, zwitterion
    }
    krio_plasma <- 1 + ib1 * 10^(pka1 - 7.4) + ia1 * 10^(7.4 - pka1) + ib2 * 10^(pka2 - 7.4) + ia2 * 10^(7.4 - pka2)
    krio_rbc <- 1 + ib1 * 10^(pka1 - 7.21) + ia1 * 10^(7.21 - pka1) + ib2 * 10^(pka2 - 7.21) + ia2 * 10^(7.21 - pka2)
    krio_660 <- 1 + ib1 * 10^(pka1 - 6.6) + ia1 * 10^(6.6 - pka1) + ib2 * 10^(pka2 - 6.6) + ia2 * 10^(6.6 - pka2)
    krio_710 <- 1 + ib1 * 10^(pka1 - 7.1) + ia1 * 10^(7.1 - pka1) + ib2 * 10^(pka2 - 7.1) + ia2 * 10^(7.1 - pka2)
    krio_700 <- 1 + ib1 * 10^(pka1 - 7.0) + ia1 * 10^(7.0 - pka1) + ib2 * 10^(pka2 - 7.0) + ia2 * 10^(7.0 - pka2)
    krio_722 <- 1 + ib1 * 10^(pka1 - 7.22) + ia1 * 10^(7.22 - pka1) + ib2 * 10^(pka2 - 7.22) + ia2 * 10^(7.22 - pka2)
    krio_723 <- 1 + ib1 * 10^(pka1 - 7.23) + ia1 * 10^(7.23 - pka1) + ib2 * 10^(pka2 - 7.23) + ia2 * 10^(7.23 - pka2)
    prot <- hsa
    protref <- 45.6
    if (pb_aag == 1) {
      prot <- c_aag
      protref <- 0.798
    }
    fup <- 1 / (1 + (((1 / fu) - 1) / protref) * prot)
    # --- perpetrator drug (perpetrator): ionisation and plasma binding (PBPK_Drug_distribution.m) ---
    kperup_eff_perpetrator <- kperup_perpetrator
    if (kperup_perpetrator == 0) {
      kperup_eff_perpetrator <- 1000 # DRUG.kPerUP = 0 replaced by 1000 1/h (PBPK_Drug_PostProcessing.m)
    }
    ib1_perpetrator <- 0
    ia1_perpetrator <- 0
    ib2_perpetrator <- 0
    ia2_perpetrator <- 0
    if (dtype_perpetrator == 1 || dtype_perpetrator == 2 || dtype_perpetrator == 6) {
      ib1_perpetrator <- 1 # base on pKa1: mono_base, di_base, zwitterion
    }
    if (dtype_perpetrator == 3 || dtype_perpetrator == 4) {
      ia1_perpetrator <- 1 # acid on pKa1: mono_acid, di_acid
    }
    if (dtype_perpetrator == 2) {
      ib2_perpetrator <- 1 # base on pKa2: di_base
    }
    if (dtype_perpetrator == 4 || dtype_perpetrator == 6) {
      ia2_perpetrator <- 1 # acid on pKa2: di_acid, zwitterion
    }
    krio_plasma_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.4) + ia1_perpetrator * 10^(7.4 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.4) + ia2_perpetrator * 10^(7.4 - pka2_perpetrator)
    krio_rbc_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.21) + ia1_perpetrator * 10^(7.21 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.21) + ia2_perpetrator * 10^(7.21 - pka2_perpetrator)
    krio_660_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 6.6) + ia1_perpetrator * 10^(6.6 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 6.6) + ia2_perpetrator * 10^(6.6 - pka2_perpetrator)
    krio_710_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.1) + ia1_perpetrator * 10^(7.1 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.1) + ia2_perpetrator * 10^(7.1 - pka2_perpetrator)
    krio_700_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.0) + ia1_perpetrator * 10^(7.0 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.0) + ia2_perpetrator * 10^(7.0 - pka2_perpetrator)
    krio_722_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.22) + ia1_perpetrator * 10^(7.22 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.22) + ia2_perpetrator * 10^(7.22 - pka2_perpetrator)
    krio_723_perpetrator <- 1 + ib1_perpetrator * 10^(pka1_perpetrator - 7.23) + ia1_perpetrator * 10^(7.23 - pka1_perpetrator) + ib2_perpetrator * 10^(pka2_perpetrator - 7.23) + ia2_perpetrator * 10^(7.23 - pka2_perpetrator)
    prot_perpetrator <- hsa
    protref_perpetrator <- 45.6
    if (pb_aag_perpetrator == 1) {
      prot_perpetrator <- c_aag
      protref_perpetrator <- 0.798
    }
    fup_perpetrator <- 1 / (1 + (((1 / fu_perpetrator) - 1) / protref_perpetrator) * prot_perpetrator)
    # --- perpetrator2 drug (perpetrator2): ionisation and plasma binding (PBPK_Drug_distribution.m) ---
    kperup_eff_perpetrator2 <- kperup_perpetrator2
    if (kperup_perpetrator2 == 0) {
      kperup_eff_perpetrator2 <- 1000 # DRUG.kPerUP = 0 replaced by 1000 1/h (PBPK_Drug_PostProcessing.m)
    }
    ib1_perpetrator2 <- 0
    ia1_perpetrator2 <- 0
    ib2_perpetrator2 <- 0
    ia2_perpetrator2 <- 0
    if (dtype_perpetrator2 == 1 || dtype_perpetrator2 == 2 || dtype_perpetrator2 == 6) {
      ib1_perpetrator2 <- 1 # base on pKa1: mono_base, di_base, zwitterion
    }
    if (dtype_perpetrator2 == 3 || dtype_perpetrator2 == 4) {
      ia1_perpetrator2 <- 1 # acid on pKa1: mono_acid, di_acid
    }
    if (dtype_perpetrator2 == 2) {
      ib2_perpetrator2 <- 1 # base on pKa2: di_base
    }
    if (dtype_perpetrator2 == 4 || dtype_perpetrator2 == 6) {
      ia2_perpetrator2 <- 1 # acid on pKa2: di_acid, zwitterion
    }
    krio_plasma_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.4) + ia1_perpetrator2 * 10^(7.4 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.4) + ia2_perpetrator2 * 10^(7.4 - pka2_perpetrator2)
    krio_rbc_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.21) + ia1_perpetrator2 * 10^(7.21 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.21) + ia2_perpetrator2 * 10^(7.21 - pka2_perpetrator2)
    krio_660_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 6.6) + ia1_perpetrator2 * 10^(6.6 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 6.6) + ia2_perpetrator2 * 10^(6.6 - pka2_perpetrator2)
    krio_710_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.1) + ia1_perpetrator2 * 10^(7.1 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.1) + ia2_perpetrator2 * 10^(7.1 - pka2_perpetrator2)
    krio_700_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.0) + ia1_perpetrator2 * 10^(7.0 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.0) + ia2_perpetrator2 * 10^(7.0 - pka2_perpetrator2)
    krio_722_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.22) + ia1_perpetrator2 * 10^(7.22 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.22) + ia2_perpetrator2 * 10^(7.22 - pka2_perpetrator2)
    krio_723_perpetrator2 <- 1 + ib1_perpetrator2 * 10^(pka1_perpetrator2 - 7.23) + ia1_perpetrator2 * 10^(7.23 - pka1_perpetrator2) + ib2_perpetrator2 * 10^(pka2_perpetrator2 - 7.23) + ia2_perpetrator2 * 10^(7.23 - pka2_perpetrator2)
    prot_perpetrator2 <- hsa
    protref_perpetrator2 <- 45.6
    if (pb_aag_perpetrator2 == 1) {
      prot_perpetrator2 <- c_aag
      protref_perpetrator2 <- 0.798
    }
    fup_perpetrator2 <- 1 / (1 + (((1 / fu_perpetrator2) - 1) / protref_perpetrator2) * prot_perpetrator2)
    # Drug 1's fup, read for every drug's KaPR (as-run DRUG.fup(d) linear
    # index into the subject-by-drug fup array, i.e. drug 1's fup of virtual
    # subject d; here each subject's own drug-1 fup, exact when the subjects
    # of a run are identical). Drug 1 is the co-administered drug with the
    # lowest index in the framework's drug library (PBPK_DefineParameters.m).
    fup_ref <- fup
    if (fupref_perpetrator == 1) {
      fup_ref <- fup_perpetrator
    }
    if (fupref_perpetrator2 == 1) {
      fup_ref <- fup_perpetrator2
    }
    # --- victim drug: Rodgers and Rowland partitioning (PBPK_Drug_distribution.m) ---
    strong <- 0
    if ((dtype == 1 || dtype == 2) && (pka1 > 7 || pka2 > 7)) {
      strong <- 1 # strong base: acidic-phospholipid binding branch
    }
    if (strong == 1) {
      logd <- 1.115 * abs(logp) - 1.35 - log10(krio_plasma) # vegetable oil:water, adipose only
      kpurbc <- ((bp * hct) + (1 - 0.91 * hct)) / fup
      kaap <- (kpurbc - ((krio_rbc / krio_plasma) * 0.666) - ((logp * 0.17 + (0.3 * logp + 0.7) * 0.29) / krio_plasma)) *
        (krio_rbc / (0.44 * (krio_rbc - 1)))
      kpu_lung <- abs((krio_660 * 0.463 / krio_plasma + 0.348 + kaap * 0.5 * (krio_660 - 1) / krio_plasma + (logp * 0.003 + (0.3 * logp + 0.7) * 0.009) / krio_plasma) * kpscalar)
      kpu_adipose <- abs((krio_710 * 0.039 / krio_plasma + 0.141 + kaap * 0.4 * (krio_710 - 1) / krio_plasma + (logd * 0.79 + (0.3 * logd + 0.7) * 0.002) / krio_plasma) * kpscalar)
      kpu_bone <- abs((krio_700 * 0.341 / krio_plasma + 0.098 + kaap * 0.67 * (krio_700 - 1) / krio_plasma + (logp * 0.074 + (0.3 * logp + 0.7) * 0.0011) / krio_plasma) * kpscalar)
      kpu_brain <- abs((krio_710 * 0.678 / krio_plasma + 0.092 + kaap * 0.4 * (krio_710 - 1) / krio_plasma + (logp * 0.051 + (0.3 * logp + 0.7) * 0.0565) / krio_plasma) * kpscalar)
      kpu_gonads <- abs((krio_700 * 0.561 / krio_plasma + 0.239 + kaap * 1.23 * (krio_700 - 1) / krio_plasma + (logp * 0.007 + (0.3 * logp + 0.7) * 0.0077) / krio_plasma) * kpscalar)
      kpu_heart <- abs((krio_710 * 0.445 / krio_plasma + 0.313 + kaap * 3.07 * (krio_710 - 1) / krio_plasma + (logp * 0.015 + (0.3 * logp + 0.7) * 0.0166) / krio_plasma) * kpscalar)
      kpu_kidney <- abs((krio_722 * 0.5 / krio_plasma + 0.283 + kaap * 2.48 * (krio_722 - 1) / krio_plasma + (logp * 0.0207 + (0.3 * logp + 0.7) * 0.0162) / krio_plasma) * kpscalar)
      kpu_muscle <- abs((krio_700 * 0.669 / krio_plasma + 0.091 + kaap * 2.49 * (krio_700 - 1) / krio_plasma + (logp * 0.0238 + (0.3 * logp + 0.7) * 0.0072) / krio_plasma) * kpscalar)
      kpu_skin <- abs((krio_700 * 0.0947 / krio_plasma + 0.623 + kaap * 1.32 * (krio_700 - 1) / krio_plasma + (logp * 0.0248 + (0.3 * logp + 0.7) * 0.0111) / krio_plasma) * kpscalar)
      kpu_thymus <- abs((krio_700 * 0.626 / krio_plasma + 0.15 + kaap * 2.3 * (krio_700 - 1) / krio_plasma + (logp * 0.017 + (0.3 * logp + 0.7) * 0.0092) / krio_plasma) * kpscalar)
      kpu_gut <- abs((krio_700 * 0.451 / krio_plasma + 0.267 + kaap * 2.84 * (krio_700 - 1) / krio_plasma + (logp * 0.0487 + (0.3 * logp + 0.7) * 0.0163) / krio_plasma) * kpscalar)
      kpu_spleen <- abs((krio_700 * 0.58 / krio_plasma + 0.208 + kaap * 2.81 * (krio_700 - 1) / krio_plasma + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma) * kpscalar)
      kpu_pancreas <- abs((krio_700 * 0.664 / krio_plasma + 0.12 + kaap * 1.67 * (krio_700 - 1) / krio_plasma + (logp * 0.041 + (0.3 * logp + 0.7) * 0.0093) / krio_plasma) * kpscalar)
      kpu_liver <- abs((krio_723 * 0.586 / krio_plasma + 0.165 + kaap * 5.09 * (krio_723 - 1) / krio_plasma + (logp * 0.0348 + (0.3 * logp + 0.7) * 0.0252) / krio_plasma) * kpscalar)
      kpu_lnode <- abs((krio_700 * 0.58 / krio_plasma + 0.208 + kaap * 2.81 * (krio_700 - 1) / krio_plasma + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma) * kpscalar)
      kpu_other <- abs((krio_700 * fi_other / krio_plasma + fe_other + kaap * ap_other * (krio_700 - 1) / krio_plasma + (logp * fnl_other + (0.3 * logp + 0.7) * fnp_other) / krio_plasma) * kpscalar)
    } else {
      # KaPR divides by DRUG.fup(d), a linear index that reads drug 1's fup
      # (fup_ref), not this drug's own (as run).
      kapr <- ((1 / fup_ref) - 1 - (logp * 0.35 + (0.3 * logp + 0.7) * 0.225) / krio_plasma) / prot
      kpu_lung <- abs((krio_660 * 0.463 / krio_plasma + 0.348 + (logp * 0.003 + (0.3 * logp + 0.7) * 0.009) / krio_plasma + kapr * 0.212 * prot) * kpscalar)
      kpu_adipose <- abs((krio_710 * 0.039 / krio_plasma + 0.141 + (logp * 0.79 + (0.3 * logp + 0.7) * 0.002) / krio_plasma + kapr * 0.021 * prot) * kpscalar)
      kpu_bone <- abs((krio_700 * 0.341 / krio_plasma + 0.098 + (logp * 0.074 + (0.3 * logp + 0.7) * 0.0011) / krio_plasma + kapr * 0.1 * prot) * kpscalar)
      kpu_brain <- abs((krio_710 * 0.678 / krio_plasma + 0.092 + (logp * 0.051 + (0.3 * logp + 0.7) * 0.0565) / krio_plasma + kapr * 0.048 * prot) * kpscalar)
      kpu_gonads <- abs((krio_700 * 0.561 / krio_plasma + 0.239 + (logp * 0.007 + (0.3 * logp + 0.7) * 0.0077) / krio_plasma + kapr * 0.048 * prot) * kpscalar)
      kpu_heart <- abs((krio_710 * 0.445 / krio_plasma + 0.313 + (logp * 0.015 + (0.3 * logp + 0.7) * 0.0166) / krio_plasma + kapr * 0.157 * prot) * kpscalar)
      kpu_kidney <- abs((krio_722 * 0.5 / krio_plasma + 0.283 + (logp * 0.0207 + (0.3 * logp + 0.7) * 0.0162) / krio_plasma + kapr * 0.13 * prot) * kpscalar)
      kpu_muscle <- abs((krio_700 * 0.669 / krio_plasma + 0.091 + (logp * 0.0238 + (0.3 * logp + 0.7) * 0.0072) / krio_plasma + kapr * 0.025 * prot) * kpscalar)
      kpu_skin <- abs((krio_700 * 0.0947 / krio_plasma + 0.623 + (logp * 0.0248 + (0.3 * logp + 0.7) * 0.0111) / krio_plasma + kapr * 0.277 * prot) * kpscalar)
      kpu_thymus <- abs((krio_700 * 0.626 / krio_plasma + 0.15 + (logp * 0.017 + (0.3 * logp + 0.7) * 0.0092) / krio_plasma + kapr * 0.075 * prot) * kpscalar)
      kpu_gut <- abs((krio_700 * 0.451 / krio_plasma + 0.267 + (logp * 0.0487 + (0.3 * logp + 0.7) * 0.0163) / krio_plasma + kapr * 0.158 * prot) * kpscalar)
      kpu_spleen <- abs((krio_700 * 0.58 / krio_plasma + 0.208 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + kapr * 0.097 * prot) * kpscalar)
      kpu_pancreas <- abs((krio_700 * 0.664 / krio_plasma + 0.12 + (logp * 0.041 + (0.3 * logp + 0.7) * 0.0093) / krio_plasma + kapr * 0.06 * prot) * kpscalar)
      kpu_liver <- abs((krio_723 * 0.586 / krio_plasma + 0.165 + (logp * 0.0348 + (0.3 * logp + 0.7) * 0.0252) / krio_plasma + kapr * 0.086 * prot) * kpscalar)
      kpu_lnode <- abs((krio_700 * 0.58 / krio_plasma + 0.208 + (logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + kapr * 0.097 * prot) * kpscalar)
      kpu_other <- abs((krio_700 * fi_other / krio_plasma + fe_other + (logp * fnl_other + (0.3 * logp + 0.7) * fnp_other) / krio_plasma + kapr * kphsa_other * prot) * kpscalar)
    }
    # Unbound fractions: fuint, fucel for albumin binders; 1 for AAG binders.
    fuint_lung <- 1 / ((1 - pb_aag) * (0.212 / 0.348) * ((1 / fup) - 1) + 1)
    fuint_adipose <- 1 / ((1 - pb_aag) * (0.021 / 0.141) * ((1 / fup) - 1) + 1)
    fuint_bone <- 1 / ((1 - pb_aag) * (0.1 / 0.098) * ((1 / fup) - 1) + 1)
    fuint_brain <- 1 / ((1 - pb_aag) * (0.048 / 0.092) * ((1 / fup) - 1) + 1)
    fuint_gonads <- 1 / ((1 - pb_aag) * (0.048 / 0.239) * ((1 / fup) - 1) + 1)
    fuint_heart <- 1 / ((1 - pb_aag) * (0.157 / 0.313) * ((1 / fup) - 1) + 1)
    fuint_kidney <- 1 / ((1 - pb_aag) * (0.13 / 0.283) * ((1 / fup) - 1) + 1)
    fuint_muscle <- 1 / ((1 - pb_aag) * (0.025 / 0.091) * ((1 / fup) - 1) + 1)
    fuint_skin <- 1 / ((1 - pb_aag) * (0.277 / 0.623) * ((1 / fup) - 1) + 1)
    fuint_thymus <- 1 / ((1 - pb_aag) * (0.075 / 0.15) * ((1 / fup) - 1) + 1)
    fuint_spleen <- 1 / ((1 - pb_aag) * (0.097 / 0.208) * ((1 / fup) - 1) + 1)
    fuint_pancreas <- 1 / ((1 - pb_aag) * (0.06 / 0.12) * ((1 / fup) - 1) + 1)
    fuint_liver <- 1 / ((1 - pb_aag) * (0.086 / 0.165) * ((1 / fup) - 1) + 1)
    fuint_lnode <- 1 / ((1 - pb_aag) * (0.097 / 0.208) * ((1 / fup) - 1) + 1)
    fuint_other <- 1 / ((1 - pb_aag) * (kphsa_other / fe_other) * ((1 / fup) - 1) + 1)
    fucel_lung <- 1 / (1 + (1 - pb_aag) * ((logp * 0.003 + (0.3 * logp + 0.7) * 0.009) / krio_plasma + 0.212 * hsa))
    fucel_adipose <- 1 / (1 + (1 - pb_aag) * ((logp * 0.79 + (0.3 * logp + 0.7) * 0.002) / krio_plasma + 0.021 * hsa))
    fucel_bone <- 1 / (1 + (1 - pb_aag) * ((logp * 0.074 + (0.3 * logp + 0.7) * 0.0011) / krio_plasma + 0.1 * hsa))
    fucel_brain <- 1 / (1 + (1 - pb_aag) * ((logp * 0.051 + (0.3 * logp + 0.7) * 0.0565) / krio_plasma + 0.048 * hsa))
    fucel_gonads <- 1 / (1 + (1 - pb_aag) * ((logp * 0.007 + (0.3 * logp + 0.7) * 0.0077) / krio_plasma + 0.048 * hsa))
    fucel_heart <- 1 / (1 + (1 - pb_aag) * ((logp * 0.015 + (0.3 * logp + 0.7) * 0.0166) / krio_plasma + 0.157 * hsa))
    fucel_kidney <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0207 + (0.3 * logp + 0.7) * 0.0162) / krio_plasma + 0.13 * hsa))
    fucel_muscle <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0238 + (0.3 * logp + 0.7) * 0.0072) / krio_plasma + 0.025 * hsa))
    fucel_skin <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0248 + (0.3 * logp + 0.7) * 0.0111) / krio_plasma + 0.277 * hsa))
    fucel_thymus <- 1 / (1 + (1 - pb_aag) * ((logp * 0.017 + (0.3 * logp + 0.7) * 0.0092) / krio_plasma + 0.075 * hsa))
    fucel_gut <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0487 + (0.3 * logp + 0.7) * 0.0163) / krio_plasma + 0.158 * hsa))
    fucel_spleen <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + 0.097 * hsa))
    fucel_pancreas <- 1 / (1 + (1 - pb_aag) * ((logp * 0.041 + (0.3 * logp + 0.7) * 0.0093) / krio_plasma + 0.06 * hsa))
    fucel_liver <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0348 + (0.3 * logp + 0.7) * 0.0252) / krio_plasma + 0.086 * hsa))
    fucel_lnode <- 1 / (1 + (1 - pb_aag) * ((logp * 0.0201 + (0.3 * logp + 0.7) * 0.0198) / krio_plasma + 0.097 * hsa))
    fucel_other <- 1 / (1 + (1 - pb_aag) * ((logp * fnl_other + (0.3 * logp + 0.7) * fnp_other) / krio_plasma + kphsa_other * hsa))
    # Membrane clearance CLin = CLout = ((Qorg - Lorg)/Kpu)/fup (L/h).
    clin_lung <- (ql_lung / kpu_lung) / fup
    clin_adipose <- (ql_adipose / kpu_adipose) / fup
    clin_bone <- (ql_bone / kpu_bone) / fup
    clin_brain <- (ql_brain / kpu_brain) / fup
    clin_gonads <- (ql_gonads / kpu_gonads) / fup
    clin_heart <- (ql_heart / kpu_heart) / fup
    clin_kidney <- (ql_kidney / kpu_kidney) / fup
    clin_muscle <- (ql_muscle / kpu_muscle) / fup
    clin_skin <- (ql_skin / kpu_skin) / fup
    clin_thymus <- (ql_thymus / kpu_thymus) / fup
    clin_gut <- (ql_gut / kpu_gut) / fup
    clin_spleen <- (ql_spleen / kpu_spleen) / fup
    clin_pancreas <- (ql_pancreas / kpu_pancreas) / fup
    clin_liver <- (ql_liver / kpu_liver) / fup
    clin_lnode <- (q_lnode / kpu_lnode) / fup
    clin_other <- (ql_other / kpu_other) / fup
    # Absorption (PBPK_Drug_absorption.m): Peff (Sun 2002) apportioned by length.
    peff <- 10^(0.6795 * log10(papp) - 0.3355)
    clab_duodenum <- surf_duodenum * 1.6 * 6.5 * peff * (len_duodenum / len_total) * 1e-4 * 3600 * 0.001
    clab_jejunum <- surf_jejunum * 1.6 * 8.6 * peff * (len_jejunum / len_total) * 1e-4 * 3600 * 0.001
    clab_ileum <- surf_ileum * 1.6 * 4.5 * peff * (len_ileum / len_total) * 1e-4 * 3600 * 0.001
    clab_colon <- surf_colon * 1.0 * 6.5 * peff * (len_colon / len_total) * 1e-4 * 3600 * 0.001 * fabscolon
    # Renal clearance on the kidney vascular (blood) concentration (PBPK_Drug_elimination.m).
    clr <- clrenal * (gfr / (130 - 10 * SEXF)) * (fup / fu)
    kiinv_cyp3a4 <- 0
    if (ki_cyp3a4 > 0) {
      kiinv_cyp3a4 <- 1 / ki_cyp3a4
    }
    kiinv_cyp2c19 <- 0
    if (ki_cyp2c19 > 0) {
      kiinv_cyp2c19 <- 1 / ki_cyp2c19
    }
    kiinv_cyp2d6 <- 0
    if (ki_cyp2d6 > 0) {
      kiinv_cyp2d6 <- 1 / ki_cyp2d6
    }
    kiinv_cyp2c8 <- 0
    if (ki_cyp2c8 > 0) {
      kiinv_cyp2c8 <- 1 / ki_cyp2c8
    }
    kiinv_cyp1a2 <- 0
    if (ki_cyp1a2 > 0) {
      kiinv_cyp1a2 <- 1 / ki_cyp1a2
    }
    kiinv_ugt1a1 <- 0
    if (ki_ugt1a1 > 0) {
      kiinv_ugt1a1 <- 1 / ki_ugt1a1
    }
    # --- perpetrator drug: Rodgers and Rowland partitioning (PBPK_Drug_distribution.m) ---
    strong_perpetrator <- 0
    if ((dtype_perpetrator == 1 || dtype_perpetrator == 2) && (pka1_perpetrator > 7 || pka2_perpetrator > 7)) {
      strong_perpetrator <- 1 # strong base: acidic-phospholipid binding branch
    }
    if (strong_perpetrator == 1) {
      logd_perpetrator <- 1.115 * abs(logp_perpetrator) - 1.35 - log10(krio_plasma_perpetrator) # vegetable oil:water, adipose only
      kpurbc_perpetrator <- ((bp_perpetrator * hct) + (1 - 0.91 * hct)) / fup_perpetrator
      kaap_perpetrator <- (kpurbc_perpetrator - ((krio_rbc_perpetrator / krio_plasma_perpetrator) * 0.666) - ((logp_perpetrator * 0.17 + (0.3 * logp_perpetrator + 0.7) * 0.29) / krio_plasma_perpetrator)) *
        (krio_rbc_perpetrator / (0.44 * (krio_rbc_perpetrator - 1)))
      kpu_lung_perpetrator <- abs((krio_660_perpetrator * 0.463 / krio_plasma_perpetrator + 0.348 + kaap_perpetrator * 0.5 * (krio_660_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.003 + (0.3 * logp_perpetrator + 0.7) * 0.009) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_adipose_perpetrator <- abs((krio_710_perpetrator * 0.039 / krio_plasma_perpetrator + 0.141 + kaap_perpetrator * 0.4 * (krio_710_perpetrator - 1) / krio_plasma_perpetrator + (logd_perpetrator * 0.79 + (0.3 * logd_perpetrator + 0.7) * 0.002) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_bone_perpetrator <- abs((krio_700_perpetrator * 0.341 / krio_plasma_perpetrator + 0.098 + kaap_perpetrator * 0.67 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.074 + (0.3 * logp_perpetrator + 0.7) * 0.0011) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_brain_perpetrator <- abs((krio_710_perpetrator * 0.678 / krio_plasma_perpetrator + 0.092 + kaap_perpetrator * 0.4 * (krio_710_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.051 + (0.3 * logp_perpetrator + 0.7) * 0.0565) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_gonads_perpetrator <- abs((krio_700_perpetrator * 0.561 / krio_plasma_perpetrator + 0.239 + kaap_perpetrator * 1.23 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.007 + (0.3 * logp_perpetrator + 0.7) * 0.0077) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_heart_perpetrator <- abs((krio_710_perpetrator * 0.445 / krio_plasma_perpetrator + 0.313 + kaap_perpetrator * 3.07 * (krio_710_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.015 + (0.3 * logp_perpetrator + 0.7) * 0.0166) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_kidney_perpetrator <- abs((krio_722_perpetrator * 0.5 / krio_plasma_perpetrator + 0.283 + kaap_perpetrator * 2.48 * (krio_722_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0207 + (0.3 * logp_perpetrator + 0.7) * 0.0162) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_muscle_perpetrator <- abs((krio_700_perpetrator * 0.669 / krio_plasma_perpetrator + 0.091 + kaap_perpetrator * 2.49 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0238 + (0.3 * logp_perpetrator + 0.7) * 0.0072) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_skin_perpetrator <- abs((krio_700_perpetrator * 0.0947 / krio_plasma_perpetrator + 0.623 + kaap_perpetrator * 1.32 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0248 + (0.3 * logp_perpetrator + 0.7) * 0.0111) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_thymus_perpetrator <- abs((krio_700_perpetrator * 0.626 / krio_plasma_perpetrator + 0.15 + kaap_perpetrator * 2.3 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.017 + (0.3 * logp_perpetrator + 0.7) * 0.0092) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_gut_perpetrator <- abs((krio_700_perpetrator * 0.451 / krio_plasma_perpetrator + 0.267 + kaap_perpetrator * 2.84 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0487 + (0.3 * logp_perpetrator + 0.7) * 0.0163) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_spleen_perpetrator <- abs((krio_700_perpetrator * 0.58 / krio_plasma_perpetrator + 0.208 + kaap_perpetrator * 2.81 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_pancreas_perpetrator <- abs((krio_700_perpetrator * 0.664 / krio_plasma_perpetrator + 0.12 + kaap_perpetrator * 1.67 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.041 + (0.3 * logp_perpetrator + 0.7) * 0.0093) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_liver_perpetrator <- abs((krio_723_perpetrator * 0.586 / krio_plasma_perpetrator + 0.165 + kaap_perpetrator * 5.09 * (krio_723_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0348 + (0.3 * logp_perpetrator + 0.7) * 0.0252) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_lnode_perpetrator <- abs((krio_700_perpetrator * 0.58 / krio_plasma_perpetrator + 0.208 + kaap_perpetrator * 2.81 * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator) * kpscalar_perpetrator)
      kpu_other_perpetrator <- abs((krio_700_perpetrator * fi_other / krio_plasma_perpetrator + fe_other + kaap_perpetrator * ap_other * (krio_700_perpetrator - 1) / krio_plasma_perpetrator + (logp_perpetrator * fnl_other + (0.3 * logp_perpetrator + 0.7) * fnp_other) / krio_plasma_perpetrator) * kpscalar_perpetrator)
    } else {
      # KaPR divides by DRUG.fup(d), a linear index that reads drug 1's fup
      # (fup_ref), not this drug's own (as run).
      kapr_perpetrator <- ((1 / fup_ref) - 1 - (logp_perpetrator * 0.35 + (0.3 * logp_perpetrator + 0.7) * 0.225) / krio_plasma_perpetrator) / prot_perpetrator
      kpu_lung_perpetrator <- abs((krio_660_perpetrator * 0.463 / krio_plasma_perpetrator + 0.348 + (logp_perpetrator * 0.003 + (0.3 * logp_perpetrator + 0.7) * 0.009) / krio_plasma_perpetrator + kapr_perpetrator * 0.212 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_adipose_perpetrator <- abs((krio_710_perpetrator * 0.039 / krio_plasma_perpetrator + 0.141 + (logp_perpetrator * 0.79 + (0.3 * logp_perpetrator + 0.7) * 0.002) / krio_plasma_perpetrator + kapr_perpetrator * 0.021 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_bone_perpetrator <- abs((krio_700_perpetrator * 0.341 / krio_plasma_perpetrator + 0.098 + (logp_perpetrator * 0.074 + (0.3 * logp_perpetrator + 0.7) * 0.0011) / krio_plasma_perpetrator + kapr_perpetrator * 0.1 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_brain_perpetrator <- abs((krio_710_perpetrator * 0.678 / krio_plasma_perpetrator + 0.092 + (logp_perpetrator * 0.051 + (0.3 * logp_perpetrator + 0.7) * 0.0565) / krio_plasma_perpetrator + kapr_perpetrator * 0.048 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_gonads_perpetrator <- abs((krio_700_perpetrator * 0.561 / krio_plasma_perpetrator + 0.239 + (logp_perpetrator * 0.007 + (0.3 * logp_perpetrator + 0.7) * 0.0077) / krio_plasma_perpetrator + kapr_perpetrator * 0.048 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_heart_perpetrator <- abs((krio_710_perpetrator * 0.445 / krio_plasma_perpetrator + 0.313 + (logp_perpetrator * 0.015 + (0.3 * logp_perpetrator + 0.7) * 0.0166) / krio_plasma_perpetrator + kapr_perpetrator * 0.157 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_kidney_perpetrator <- abs((krio_722_perpetrator * 0.5 / krio_plasma_perpetrator + 0.283 + (logp_perpetrator * 0.0207 + (0.3 * logp_perpetrator + 0.7) * 0.0162) / krio_plasma_perpetrator + kapr_perpetrator * 0.13 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_muscle_perpetrator <- abs((krio_700_perpetrator * 0.669 / krio_plasma_perpetrator + 0.091 + (logp_perpetrator * 0.0238 + (0.3 * logp_perpetrator + 0.7) * 0.0072) / krio_plasma_perpetrator + kapr_perpetrator * 0.025 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_skin_perpetrator <- abs((krio_700_perpetrator * 0.0947 / krio_plasma_perpetrator + 0.623 + (logp_perpetrator * 0.0248 + (0.3 * logp_perpetrator + 0.7) * 0.0111) / krio_plasma_perpetrator + kapr_perpetrator * 0.277 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_thymus_perpetrator <- abs((krio_700_perpetrator * 0.626 / krio_plasma_perpetrator + 0.15 + (logp_perpetrator * 0.017 + (0.3 * logp_perpetrator + 0.7) * 0.0092) / krio_plasma_perpetrator + kapr_perpetrator * 0.075 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_gut_perpetrator <- abs((krio_700_perpetrator * 0.451 / krio_plasma_perpetrator + 0.267 + (logp_perpetrator * 0.0487 + (0.3 * logp_perpetrator + 0.7) * 0.0163) / krio_plasma_perpetrator + kapr_perpetrator * 0.158 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_spleen_perpetrator <- abs((krio_700_perpetrator * 0.58 / krio_plasma_perpetrator + 0.208 + (logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator + kapr_perpetrator * 0.097 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_pancreas_perpetrator <- abs((krio_700_perpetrator * 0.664 / krio_plasma_perpetrator + 0.12 + (logp_perpetrator * 0.041 + (0.3 * logp_perpetrator + 0.7) * 0.0093) / krio_plasma_perpetrator + kapr_perpetrator * 0.06 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_liver_perpetrator <- abs((krio_723_perpetrator * 0.586 / krio_plasma_perpetrator + 0.165 + (logp_perpetrator * 0.0348 + (0.3 * logp_perpetrator + 0.7) * 0.0252) / krio_plasma_perpetrator + kapr_perpetrator * 0.086 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_lnode_perpetrator <- abs((krio_700_perpetrator * 0.58 / krio_plasma_perpetrator + 0.208 + (logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator + kapr_perpetrator * 0.097 * prot_perpetrator) * kpscalar_perpetrator)
      kpu_other_perpetrator <- abs((krio_700_perpetrator * fi_other / krio_plasma_perpetrator + fe_other + (logp_perpetrator * fnl_other + (0.3 * logp_perpetrator + 0.7) * fnp_other) / krio_plasma_perpetrator + kapr_perpetrator * kphsa_other * prot_perpetrator) * kpscalar_perpetrator)
    }
    # Unbound fractions: fuint, fucel for albumin binders; 1 for AAG binders.
    fuint_lung_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.212 / 0.348) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_adipose_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.021 / 0.141) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_bone_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.1 / 0.098) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_brain_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.048 / 0.092) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_gonads_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.048 / 0.239) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_heart_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.157 / 0.313) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_kidney_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.13 / 0.283) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_muscle_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.025 / 0.091) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_skin_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.277 / 0.623) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_thymus_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.075 / 0.15) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_spleen_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.097 / 0.208) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_pancreas_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.06 / 0.12) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_liver_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.086 / 0.165) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_lnode_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (0.097 / 0.208) * ((1 / fup_perpetrator) - 1) + 1)
    fuint_other_perpetrator <- 1 / ((1 - pb_aag_perpetrator) * (kphsa_other / fe_other) * ((1 / fup_perpetrator) - 1) + 1)
    fucel_lung_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.003 + (0.3 * logp_perpetrator + 0.7) * 0.009) / krio_plasma_perpetrator + 0.212 * hsa))
    fucel_adipose_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.79 + (0.3 * logp_perpetrator + 0.7) * 0.002) / krio_plasma_perpetrator + 0.021 * hsa))
    fucel_bone_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.074 + (0.3 * logp_perpetrator + 0.7) * 0.0011) / krio_plasma_perpetrator + 0.1 * hsa))
    fucel_brain_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.051 + (0.3 * logp_perpetrator + 0.7) * 0.0565) / krio_plasma_perpetrator + 0.048 * hsa))
    fucel_gonads_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.007 + (0.3 * logp_perpetrator + 0.7) * 0.0077) / krio_plasma_perpetrator + 0.048 * hsa))
    fucel_heart_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.015 + (0.3 * logp_perpetrator + 0.7) * 0.0166) / krio_plasma_perpetrator + 0.157 * hsa))
    fucel_kidney_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0207 + (0.3 * logp_perpetrator + 0.7) * 0.0162) / krio_plasma_perpetrator + 0.13 * hsa))
    fucel_muscle_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0238 + (0.3 * logp_perpetrator + 0.7) * 0.0072) / krio_plasma_perpetrator + 0.025 * hsa))
    fucel_skin_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0248 + (0.3 * logp_perpetrator + 0.7) * 0.0111) / krio_plasma_perpetrator + 0.277 * hsa))
    fucel_thymus_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.017 + (0.3 * logp_perpetrator + 0.7) * 0.0092) / krio_plasma_perpetrator + 0.075 * hsa))
    fucel_gut_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0487 + (0.3 * logp_perpetrator + 0.7) * 0.0163) / krio_plasma_perpetrator + 0.158 * hsa))
    fucel_spleen_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator + 0.097 * hsa))
    fucel_pancreas_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.041 + (0.3 * logp_perpetrator + 0.7) * 0.0093) / krio_plasma_perpetrator + 0.06 * hsa))
    fucel_liver_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0348 + (0.3 * logp_perpetrator + 0.7) * 0.0252) / krio_plasma_perpetrator + 0.086 * hsa))
    fucel_lnode_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * 0.0201 + (0.3 * logp_perpetrator + 0.7) * 0.0198) / krio_plasma_perpetrator + 0.097 * hsa))
    fucel_other_perpetrator <- 1 / (1 + (1 - pb_aag_perpetrator) * ((logp_perpetrator * fnl_other + (0.3 * logp_perpetrator + 0.7) * fnp_other) / krio_plasma_perpetrator + kphsa_other * hsa))
    # Membrane clearance CLin = CLout = ((Qorg - Lorg)/Kpu)/fup (L/h).
    clin_lung_perpetrator <- (ql_lung / kpu_lung_perpetrator) / fup_perpetrator
    clin_adipose_perpetrator <- (ql_adipose / kpu_adipose_perpetrator) / fup_perpetrator
    clin_bone_perpetrator <- (ql_bone / kpu_bone_perpetrator) / fup_perpetrator
    clin_brain_perpetrator <- (ql_brain / kpu_brain_perpetrator) / fup_perpetrator
    clin_gonads_perpetrator <- (ql_gonads / kpu_gonads_perpetrator) / fup_perpetrator
    clin_heart_perpetrator <- (ql_heart / kpu_heart_perpetrator) / fup_perpetrator
    clin_kidney_perpetrator <- (ql_kidney / kpu_kidney_perpetrator) / fup_perpetrator
    clin_muscle_perpetrator <- (ql_muscle / kpu_muscle_perpetrator) / fup_perpetrator
    clin_skin_perpetrator <- (ql_skin / kpu_skin_perpetrator) / fup_perpetrator
    clin_thymus_perpetrator <- (ql_thymus / kpu_thymus_perpetrator) / fup_perpetrator
    clin_gut_perpetrator <- (ql_gut / kpu_gut_perpetrator) / fup_perpetrator
    clin_spleen_perpetrator <- (ql_spleen / kpu_spleen_perpetrator) / fup_perpetrator
    clin_pancreas_perpetrator <- (ql_pancreas / kpu_pancreas_perpetrator) / fup_perpetrator
    clin_liver_perpetrator <- (ql_liver / kpu_liver_perpetrator) / fup_perpetrator
    clin_lnode_perpetrator <- (q_lnode / kpu_lnode_perpetrator) / fup_perpetrator
    clin_other_perpetrator <- (ql_other / kpu_other_perpetrator) / fup_perpetrator
    # Absorption (PBPK_Drug_absorption.m): Peff (Sun 2002) apportioned by length.
    peff_perpetrator <- 10^(0.6795 * log10(papp_perpetrator) - 0.3355)
    clab_duodenum_perpetrator <- surf_duodenum * 1.6 * 6.5 * peff_perpetrator * (len_duodenum / len_total) * 1e-4 * 3600 * 0.001
    clab_jejunum_perpetrator <- surf_jejunum * 1.6 * 8.6 * peff_perpetrator * (len_jejunum / len_total) * 1e-4 * 3600 * 0.001
    clab_ileum_perpetrator <- surf_ileum * 1.6 * 4.5 * peff_perpetrator * (len_ileum / len_total) * 1e-4 * 3600 * 0.001
    clab_colon_perpetrator <- surf_colon * 1.0 * 6.5 * peff_perpetrator * (len_colon / len_total) * 1e-4 * 3600 * 0.001 * fabscolon_perpetrator
    # Renal clearance on the kidney vascular (blood) concentration (PBPK_Drug_elimination.m).
    clr_perpetrator <- clrenal_perpetrator * (gfr / (130 - 10 * SEXF)) * (fup_perpetrator / fu_perpetrator)
    kiinv_cyp3a4_perpetrator <- 0
    if (ki_cyp3a4_perpetrator > 0) {
      kiinv_cyp3a4_perpetrator <- 1 / ki_cyp3a4_perpetrator
    }
    kiinv_cyp2c19_perpetrator <- 0
    if (ki_cyp2c19_perpetrator > 0) {
      kiinv_cyp2c19_perpetrator <- 1 / ki_cyp2c19_perpetrator
    }
    kiinv_cyp2d6_perpetrator <- 0
    if (ki_cyp2d6_perpetrator > 0) {
      kiinv_cyp2d6_perpetrator <- 1 / ki_cyp2d6_perpetrator
    }
    kiinv_cyp2c8_perpetrator <- 0
    if (ki_cyp2c8_perpetrator > 0) {
      kiinv_cyp2c8_perpetrator <- 1 / ki_cyp2c8_perpetrator
    }
    kiinv_cyp1a2_perpetrator <- 0
    if (ki_cyp1a2_perpetrator > 0) {
      kiinv_cyp1a2_perpetrator <- 1 / ki_cyp1a2_perpetrator
    }
    kiinv_ugt1a1_perpetrator <- 0
    if (ki_ugt1a1_perpetrator > 0) {
      kiinv_ugt1a1_perpetrator <- 1 / ki_ugt1a1_perpetrator
    }
    # --- perpetrator2 drug: Rodgers and Rowland partitioning (PBPK_Drug_distribution.m) ---
    strong_perpetrator2 <- 0
    if ((dtype_perpetrator2 == 1 || dtype_perpetrator2 == 2) && (pka1_perpetrator2 > 7 || pka2_perpetrator2 > 7)) {
      strong_perpetrator2 <- 1 # strong base: acidic-phospholipid binding branch
    }
    if (strong_perpetrator2 == 1) {
      logd_perpetrator2 <- 1.115 * abs(logp_perpetrator2) - 1.35 - log10(krio_plasma_perpetrator2) # vegetable oil:water, adipose only
      kpurbc_perpetrator2 <- ((bp_perpetrator2 * hct) + (1 - 0.91 * hct)) / fup_perpetrator2
      kaap_perpetrator2 <- (kpurbc_perpetrator2 - ((krio_rbc_perpetrator2 / krio_plasma_perpetrator2) * 0.666) - ((logp_perpetrator2 * 0.17 + (0.3 * logp_perpetrator2 + 0.7) * 0.29) / krio_plasma_perpetrator2)) *
        (krio_rbc_perpetrator2 / (0.44 * (krio_rbc_perpetrator2 - 1)))
      kpu_lung_perpetrator2 <- abs((krio_660_perpetrator2 * 0.463 / krio_plasma_perpetrator2 + 0.348 + kaap_perpetrator2 * 0.5 * (krio_660_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.003 + (0.3 * logp_perpetrator2 + 0.7) * 0.009) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_adipose_perpetrator2 <- abs((krio_710_perpetrator2 * 0.039 / krio_plasma_perpetrator2 + 0.141 + kaap_perpetrator2 * 0.4 * (krio_710_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logd_perpetrator2 * 0.79 + (0.3 * logd_perpetrator2 + 0.7) * 0.002) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_bone_perpetrator2 <- abs((krio_700_perpetrator2 * 0.341 / krio_plasma_perpetrator2 + 0.098 + kaap_perpetrator2 * 0.67 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.074 + (0.3 * logp_perpetrator2 + 0.7) * 0.0011) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_brain_perpetrator2 <- abs((krio_710_perpetrator2 * 0.678 / krio_plasma_perpetrator2 + 0.092 + kaap_perpetrator2 * 0.4 * (krio_710_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.051 + (0.3 * logp_perpetrator2 + 0.7) * 0.0565) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_gonads_perpetrator2 <- abs((krio_700_perpetrator2 * 0.561 / krio_plasma_perpetrator2 + 0.239 + kaap_perpetrator2 * 1.23 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.007 + (0.3 * logp_perpetrator2 + 0.7) * 0.0077) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_heart_perpetrator2 <- abs((krio_710_perpetrator2 * 0.445 / krio_plasma_perpetrator2 + 0.313 + kaap_perpetrator2 * 3.07 * (krio_710_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.015 + (0.3 * logp_perpetrator2 + 0.7) * 0.0166) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_kidney_perpetrator2 <- abs((krio_722_perpetrator2 * 0.5 / krio_plasma_perpetrator2 + 0.283 + kaap_perpetrator2 * 2.48 * (krio_722_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0207 + (0.3 * logp_perpetrator2 + 0.7) * 0.0162) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_muscle_perpetrator2 <- abs((krio_700_perpetrator2 * 0.669 / krio_plasma_perpetrator2 + 0.091 + kaap_perpetrator2 * 2.49 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0238 + (0.3 * logp_perpetrator2 + 0.7) * 0.0072) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_skin_perpetrator2 <- abs((krio_700_perpetrator2 * 0.0947 / krio_plasma_perpetrator2 + 0.623 + kaap_perpetrator2 * 1.32 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0248 + (0.3 * logp_perpetrator2 + 0.7) * 0.0111) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_thymus_perpetrator2 <- abs((krio_700_perpetrator2 * 0.626 / krio_plasma_perpetrator2 + 0.15 + kaap_perpetrator2 * 2.3 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.017 + (0.3 * logp_perpetrator2 + 0.7) * 0.0092) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_gut_perpetrator2 <- abs((krio_700_perpetrator2 * 0.451 / krio_plasma_perpetrator2 + 0.267 + kaap_perpetrator2 * 2.84 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0487 + (0.3 * logp_perpetrator2 + 0.7) * 0.0163) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_spleen_perpetrator2 <- abs((krio_700_perpetrator2 * 0.58 / krio_plasma_perpetrator2 + 0.208 + kaap_perpetrator2 * 2.81 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_pancreas_perpetrator2 <- abs((krio_700_perpetrator2 * 0.664 / krio_plasma_perpetrator2 + 0.12 + kaap_perpetrator2 * 1.67 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.041 + (0.3 * logp_perpetrator2 + 0.7) * 0.0093) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_liver_perpetrator2 <- abs((krio_723_perpetrator2 * 0.586 / krio_plasma_perpetrator2 + 0.165 + kaap_perpetrator2 * 5.09 * (krio_723_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0348 + (0.3 * logp_perpetrator2 + 0.7) * 0.0252) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_lnode_perpetrator2 <- abs((krio_700_perpetrator2 * 0.58 / krio_plasma_perpetrator2 + 0.208 + kaap_perpetrator2 * 2.81 * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
      kpu_other_perpetrator2 <- abs((krio_700_perpetrator2 * fi_other / krio_plasma_perpetrator2 + fe_other + kaap_perpetrator2 * ap_other * (krio_700_perpetrator2 - 1) / krio_plasma_perpetrator2 + (logp_perpetrator2 * fnl_other + (0.3 * logp_perpetrator2 + 0.7) * fnp_other) / krio_plasma_perpetrator2) * kpscalar_perpetrator2)
    } else {
      # KaPR divides by DRUG.fup(d), a linear index that reads drug 1's fup
      # (fup_ref), not this drug's own (as run).
      kapr_perpetrator2 <- ((1 / fup_ref) - 1 - (logp_perpetrator2 * 0.35 + (0.3 * logp_perpetrator2 + 0.7) * 0.225) / krio_plasma_perpetrator2) / prot_perpetrator2
      kpu_lung_perpetrator2 <- abs((krio_660_perpetrator2 * 0.463 / krio_plasma_perpetrator2 + 0.348 + (logp_perpetrator2 * 0.003 + (0.3 * logp_perpetrator2 + 0.7) * 0.009) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.212 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_adipose_perpetrator2 <- abs((krio_710_perpetrator2 * 0.039 / krio_plasma_perpetrator2 + 0.141 + (logp_perpetrator2 * 0.79 + (0.3 * logp_perpetrator2 + 0.7) * 0.002) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.021 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_bone_perpetrator2 <- abs((krio_700_perpetrator2 * 0.341 / krio_plasma_perpetrator2 + 0.098 + (logp_perpetrator2 * 0.074 + (0.3 * logp_perpetrator2 + 0.7) * 0.0011) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.1 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_brain_perpetrator2 <- abs((krio_710_perpetrator2 * 0.678 / krio_plasma_perpetrator2 + 0.092 + (logp_perpetrator2 * 0.051 + (0.3 * logp_perpetrator2 + 0.7) * 0.0565) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.048 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_gonads_perpetrator2 <- abs((krio_700_perpetrator2 * 0.561 / krio_plasma_perpetrator2 + 0.239 + (logp_perpetrator2 * 0.007 + (0.3 * logp_perpetrator2 + 0.7) * 0.0077) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.048 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_heart_perpetrator2 <- abs((krio_710_perpetrator2 * 0.445 / krio_plasma_perpetrator2 + 0.313 + (logp_perpetrator2 * 0.015 + (0.3 * logp_perpetrator2 + 0.7) * 0.0166) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.157 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_kidney_perpetrator2 <- abs((krio_722_perpetrator2 * 0.5 / krio_plasma_perpetrator2 + 0.283 + (logp_perpetrator2 * 0.0207 + (0.3 * logp_perpetrator2 + 0.7) * 0.0162) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.13 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_muscle_perpetrator2 <- abs((krio_700_perpetrator2 * 0.669 / krio_plasma_perpetrator2 + 0.091 + (logp_perpetrator2 * 0.0238 + (0.3 * logp_perpetrator2 + 0.7) * 0.0072) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.025 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_skin_perpetrator2 <- abs((krio_700_perpetrator2 * 0.0947 / krio_plasma_perpetrator2 + 0.623 + (logp_perpetrator2 * 0.0248 + (0.3 * logp_perpetrator2 + 0.7) * 0.0111) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.277 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_thymus_perpetrator2 <- abs((krio_700_perpetrator2 * 0.626 / krio_plasma_perpetrator2 + 0.15 + (logp_perpetrator2 * 0.017 + (0.3 * logp_perpetrator2 + 0.7) * 0.0092) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.075 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_gut_perpetrator2 <- abs((krio_700_perpetrator2 * 0.451 / krio_plasma_perpetrator2 + 0.267 + (logp_perpetrator2 * 0.0487 + (0.3 * logp_perpetrator2 + 0.7) * 0.0163) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.158 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_spleen_perpetrator2 <- abs((krio_700_perpetrator2 * 0.58 / krio_plasma_perpetrator2 + 0.208 + (logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.097 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_pancreas_perpetrator2 <- abs((krio_700_perpetrator2 * 0.664 / krio_plasma_perpetrator2 + 0.12 + (logp_perpetrator2 * 0.041 + (0.3 * logp_perpetrator2 + 0.7) * 0.0093) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.06 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_liver_perpetrator2 <- abs((krio_723_perpetrator2 * 0.586 / krio_plasma_perpetrator2 + 0.165 + (logp_perpetrator2 * 0.0348 + (0.3 * logp_perpetrator2 + 0.7) * 0.0252) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.086 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_lnode_perpetrator2 <- abs((krio_700_perpetrator2 * 0.58 / krio_plasma_perpetrator2 + 0.208 + (logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2 + kapr_perpetrator2 * 0.097 * prot_perpetrator2) * kpscalar_perpetrator2)
      kpu_other_perpetrator2 <- abs((krio_700_perpetrator2 * fi_other / krio_plasma_perpetrator2 + fe_other + (logp_perpetrator2 * fnl_other + (0.3 * logp_perpetrator2 + 0.7) * fnp_other) / krio_plasma_perpetrator2 + kapr_perpetrator2 * kphsa_other * prot_perpetrator2) * kpscalar_perpetrator2)
    }
    # Unbound fractions: fuint, fucel for albumin binders; 1 for AAG binders.
    fuint_lung_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.212 / 0.348) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_adipose_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.021 / 0.141) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_bone_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.1 / 0.098) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_brain_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.048 / 0.092) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_gonads_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.048 / 0.239) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_heart_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.157 / 0.313) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_kidney_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.13 / 0.283) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_muscle_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.025 / 0.091) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_skin_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.277 / 0.623) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_thymus_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.075 / 0.15) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_spleen_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.097 / 0.208) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_pancreas_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.06 / 0.12) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_liver_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.086 / 0.165) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_lnode_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (0.097 / 0.208) * ((1 / fup_perpetrator2) - 1) + 1)
    fuint_other_perpetrator2 <- 1 / ((1 - pb_aag_perpetrator2) * (kphsa_other / fe_other) * ((1 / fup_perpetrator2) - 1) + 1)
    fucel_lung_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.003 + (0.3 * logp_perpetrator2 + 0.7) * 0.009) / krio_plasma_perpetrator2 + 0.212 * hsa))
    fucel_adipose_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.79 + (0.3 * logp_perpetrator2 + 0.7) * 0.002) / krio_plasma_perpetrator2 + 0.021 * hsa))
    fucel_bone_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.074 + (0.3 * logp_perpetrator2 + 0.7) * 0.0011) / krio_plasma_perpetrator2 + 0.1 * hsa))
    fucel_brain_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.051 + (0.3 * logp_perpetrator2 + 0.7) * 0.0565) / krio_plasma_perpetrator2 + 0.048 * hsa))
    fucel_gonads_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.007 + (0.3 * logp_perpetrator2 + 0.7) * 0.0077) / krio_plasma_perpetrator2 + 0.048 * hsa))
    fucel_heart_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.015 + (0.3 * logp_perpetrator2 + 0.7) * 0.0166) / krio_plasma_perpetrator2 + 0.157 * hsa))
    fucel_kidney_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0207 + (0.3 * logp_perpetrator2 + 0.7) * 0.0162) / krio_plasma_perpetrator2 + 0.13 * hsa))
    fucel_muscle_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0238 + (0.3 * logp_perpetrator2 + 0.7) * 0.0072) / krio_plasma_perpetrator2 + 0.025 * hsa))
    fucel_skin_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0248 + (0.3 * logp_perpetrator2 + 0.7) * 0.0111) / krio_plasma_perpetrator2 + 0.277 * hsa))
    fucel_thymus_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.017 + (0.3 * logp_perpetrator2 + 0.7) * 0.0092) / krio_plasma_perpetrator2 + 0.075 * hsa))
    fucel_gut_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0487 + (0.3 * logp_perpetrator2 + 0.7) * 0.0163) / krio_plasma_perpetrator2 + 0.158 * hsa))
    fucel_spleen_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2 + 0.097 * hsa))
    fucel_pancreas_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.041 + (0.3 * logp_perpetrator2 + 0.7) * 0.0093) / krio_plasma_perpetrator2 + 0.06 * hsa))
    fucel_liver_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0348 + (0.3 * logp_perpetrator2 + 0.7) * 0.0252) / krio_plasma_perpetrator2 + 0.086 * hsa))
    fucel_lnode_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * 0.0201 + (0.3 * logp_perpetrator2 + 0.7) * 0.0198) / krio_plasma_perpetrator2 + 0.097 * hsa))
    fucel_other_perpetrator2 <- 1 / (1 + (1 - pb_aag_perpetrator2) * ((logp_perpetrator2 * fnl_other + (0.3 * logp_perpetrator2 + 0.7) * fnp_other) / krio_plasma_perpetrator2 + kphsa_other * hsa))
    # Membrane clearance CLin = CLout = ((Qorg - Lorg)/Kpu)/fup (L/h).
    clin_lung_perpetrator2 <- (ql_lung / kpu_lung_perpetrator2) / fup_perpetrator2
    clin_adipose_perpetrator2 <- (ql_adipose / kpu_adipose_perpetrator2) / fup_perpetrator2
    clin_bone_perpetrator2 <- (ql_bone / kpu_bone_perpetrator2) / fup_perpetrator2
    clin_brain_perpetrator2 <- (ql_brain / kpu_brain_perpetrator2) / fup_perpetrator2
    clin_gonads_perpetrator2 <- (ql_gonads / kpu_gonads_perpetrator2) / fup_perpetrator2
    clin_heart_perpetrator2 <- (ql_heart / kpu_heart_perpetrator2) / fup_perpetrator2
    clin_kidney_perpetrator2 <- (ql_kidney / kpu_kidney_perpetrator2) / fup_perpetrator2
    clin_muscle_perpetrator2 <- (ql_muscle / kpu_muscle_perpetrator2) / fup_perpetrator2
    clin_skin_perpetrator2 <- (ql_skin / kpu_skin_perpetrator2) / fup_perpetrator2
    clin_thymus_perpetrator2 <- (ql_thymus / kpu_thymus_perpetrator2) / fup_perpetrator2
    clin_gut_perpetrator2 <- (ql_gut / kpu_gut_perpetrator2) / fup_perpetrator2
    clin_spleen_perpetrator2 <- (ql_spleen / kpu_spleen_perpetrator2) / fup_perpetrator2
    clin_pancreas_perpetrator2 <- (ql_pancreas / kpu_pancreas_perpetrator2) / fup_perpetrator2
    clin_liver_perpetrator2 <- (ql_liver / kpu_liver_perpetrator2) / fup_perpetrator2
    clin_lnode_perpetrator2 <- (q_lnode / kpu_lnode_perpetrator2) / fup_perpetrator2
    clin_other_perpetrator2 <- (ql_other / kpu_other_perpetrator2) / fup_perpetrator2
    # Absorption (PBPK_Drug_absorption.m): Peff (Sun 2002) apportioned by length.
    peff_perpetrator2 <- 10^(0.6795 * log10(papp_perpetrator2) - 0.3355)
    clab_duodenum_perpetrator2 <- surf_duodenum * 1.6 * 6.5 * peff_perpetrator2 * (len_duodenum / len_total) * 1e-4 * 3600 * 0.001
    clab_jejunum_perpetrator2 <- surf_jejunum * 1.6 * 8.6 * peff_perpetrator2 * (len_jejunum / len_total) * 1e-4 * 3600 * 0.001
    clab_ileum_perpetrator2 <- surf_ileum * 1.6 * 4.5 * peff_perpetrator2 * (len_ileum / len_total) * 1e-4 * 3600 * 0.001
    clab_colon_perpetrator2 <- surf_colon * 1.0 * 6.5 * peff_perpetrator2 * (len_colon / len_total) * 1e-4 * 3600 * 0.001 * fabscolon_perpetrator2
    # Renal clearance on the kidney vascular (blood) concentration (PBPK_Drug_elimination.m).
    clr_perpetrator2 <- clrenal_perpetrator2 * (gfr / (130 - 10 * SEXF)) * (fup_perpetrator2 / fu_perpetrator2)
    kiinv_cyp3a4_perpetrator2 <- 0
    if (ki_cyp3a4_perpetrator2 > 0) {
      kiinv_cyp3a4_perpetrator2 <- 1 / ki_cyp3a4_perpetrator2
    }
    kiinv_cyp2c19_perpetrator2 <- 0
    if (ki_cyp2c19_perpetrator2 > 0) {
      kiinv_cyp2c19_perpetrator2 <- 1 / ki_cyp2c19_perpetrator2
    }
    kiinv_cyp2d6_perpetrator2 <- 0
    if (ki_cyp2d6_perpetrator2 > 0) {
      kiinv_cyp2d6_perpetrator2 <- 1 / ki_cyp2d6_perpetrator2
    }
    kiinv_cyp2c8_perpetrator2 <- 0
    if (ki_cyp2c8_perpetrator2 > 0) {
      kiinv_cyp2c8_perpetrator2 <- 1 / ki_cyp2c8_perpetrator2
    }
    kiinv_cyp1a2_perpetrator2 <- 0
    if (ki_cyp1a2_perpetrator2 > 0) {
      kiinv_cyp1a2_perpetrator2 <- 1 / ki_cyp1a2_perpetrator2
    }
    kiinv_ugt1a1_perpetrator2 <- 0
    if (ki_ugt1a1_perpetrator2 > 0) {
      kiinv_ugt1a1_perpetrator2 <- 1 / ki_ugt1a1_perpetrator2
    }

    # Concentrations (mg/L; vascular states hold blood concentrations).
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    cv_lung <- lung_vas / vv_lung
    ce_lung <- lung_ew / ve_lung
    ci_lung <- lung_iw / vi_lung
    cv_adipose <- adipose_vas / vv_adipose
    ce_adipose <- adipose_ew / ve_adipose
    ci_adipose <- adipose_iw / vi_adipose
    cv_bone <- bone_vas / vv_bone
    ce_bone <- bone_ew / ve_bone
    ci_bone <- bone_iw / vi_bone
    cv_brain <- brain_vas / vv_brain
    ce_brain <- brain_ew / ve_brain
    ci_brain <- brain_iw / vi_brain
    cv_gonads <- gonads_vas / vv_gonads
    ce_gonads <- gonads_ew / ve_gonads
    ci_gonads <- gonads_iw / vi_gonads
    cv_heart <- heart_vas / vv_heart
    ce_heart <- heart_ew / ve_heart
    ci_heart <- heart_iw / vi_heart
    cv_kidney <- kidney_vas / vv_kidney
    ce_kidney <- kidney_ew / ve_kidney
    ci_kidney <- kidney_iw / vi_kidney
    cv_muscle <- muscle_vas / vv_muscle
    ce_muscle <- muscle_ew / ve_muscle
    ci_muscle <- muscle_iw / vi_muscle
    cv_skin <- skin_vas / vv_skin
    ce_skin <- skin_ew / ve_skin
    ci_skin <- skin_iw / vi_skin
    cv_thymus <- thymus_vas / vv_thymus
    ce_thymus <- thymus_ew / ve_thymus
    ci_thymus <- thymus_iw / vi_thymus
    cv_gut <- gut_vas / vv_gut
    ce_gut <- gut_ew / ve_gut
    cv_spleen <- spleen_vas / vv_spleen
    ce_spleen <- spleen_ew / ve_spleen
    ci_spleen <- spleen_iw / vi_spleen
    cv_pancreas <- pancreas_vas / vv_pancreas
    ce_pancreas <- pancreas_ew / ve_pancreas
    ci_pancreas <- pancreas_iw / vi_pancreas
    cv_liver <- liver_vas / vv_liver
    ce_liver <- liver_ew / ve_liver
    ci_liver <- liver_iw / vi_liver
    cv_lnode <- lnode_vas / vv_lnode
    ce_lnode <- lnode_ew / ve_lnode
    ci_lnode <- lnode_iw / vi_lnode
    cv_other <- other_vas / vv_other
    ce_other <- other_ew / ve_other
    ci_other <- other_iw / vi_other
    cent_duodenum <- duodenum_enterocyte / vent_duodenum
    cent_jejunum <- jejunum_enterocyte / vent_jejunum
    cent_ileum <- ileum_enterocyte / vent_ileum
    cent_colon <- colon_enterocyte / vent_colon
    cli_um <- ci_liver * 1000 / mw # total intracellular liver, uM
    cul <- cli_um * fucel_liver # unbound, uM
    cdu_um <- cent_duodenum * 1000 / mw
    cudu <- cdu_um * fucel_gut
    cje_um <- cent_jejunum * 1000 / mw
    cuje <- cje_um * fucel_gut
    cil_um <- cent_ileum * 1000 / mw
    cuil <- cil_um * fucel_gut
    c_venous_perpetrator <- venous_perpetrator / v_venous
    c_arterial_perpetrator <- arterial_perpetrator / v_arterial
    cv_lung_perpetrator <- lung_vas_perpetrator / vv_lung
    ce_lung_perpetrator <- lung_ew_perpetrator / ve_lung
    ci_lung_perpetrator <- lung_iw_perpetrator / vi_lung
    cv_adipose_perpetrator <- adipose_vas_perpetrator / vv_adipose
    ce_adipose_perpetrator <- adipose_ew_perpetrator / ve_adipose
    ci_adipose_perpetrator <- adipose_iw_perpetrator / vi_adipose
    cv_bone_perpetrator <- bone_vas_perpetrator / vv_bone
    ce_bone_perpetrator <- bone_ew_perpetrator / ve_bone
    ci_bone_perpetrator <- bone_iw_perpetrator / vi_bone
    cv_brain_perpetrator <- brain_vas_perpetrator / vv_brain
    ce_brain_perpetrator <- brain_ew_perpetrator / ve_brain
    ci_brain_perpetrator <- brain_iw_perpetrator / vi_brain
    cv_gonads_perpetrator <- gonads_vas_perpetrator / vv_gonads
    ce_gonads_perpetrator <- gonads_ew_perpetrator / ve_gonads
    ci_gonads_perpetrator <- gonads_iw_perpetrator / vi_gonads
    cv_heart_perpetrator <- heart_vas_perpetrator / vv_heart
    ce_heart_perpetrator <- heart_ew_perpetrator / ve_heart
    ci_heart_perpetrator <- heart_iw_perpetrator / vi_heart
    cv_kidney_perpetrator <- kidney_vas_perpetrator / vv_kidney
    ce_kidney_perpetrator <- kidney_ew_perpetrator / ve_kidney
    ci_kidney_perpetrator <- kidney_iw_perpetrator / vi_kidney
    cv_muscle_perpetrator <- muscle_vas_perpetrator / vv_muscle
    ce_muscle_perpetrator <- muscle_ew_perpetrator / ve_muscle
    ci_muscle_perpetrator <- muscle_iw_perpetrator / vi_muscle
    cv_skin_perpetrator <- skin_vas_perpetrator / vv_skin
    ce_skin_perpetrator <- skin_ew_perpetrator / ve_skin
    ci_skin_perpetrator <- skin_iw_perpetrator / vi_skin
    cv_thymus_perpetrator <- thymus_vas_perpetrator / vv_thymus
    ce_thymus_perpetrator <- thymus_ew_perpetrator / ve_thymus
    ci_thymus_perpetrator <- thymus_iw_perpetrator / vi_thymus
    cv_gut_perpetrator <- gut_vas_perpetrator / vv_gut
    ce_gut_perpetrator <- gut_ew_perpetrator / ve_gut
    cv_spleen_perpetrator <- spleen_vas_perpetrator / vv_spleen
    ce_spleen_perpetrator <- spleen_ew_perpetrator / ve_spleen
    ci_spleen_perpetrator <- spleen_iw_perpetrator / vi_spleen
    cv_pancreas_perpetrator <- pancreas_vas_perpetrator / vv_pancreas
    ce_pancreas_perpetrator <- pancreas_ew_perpetrator / ve_pancreas
    ci_pancreas_perpetrator <- pancreas_iw_perpetrator / vi_pancreas
    cv_liver_perpetrator <- liver_vas_perpetrator / vv_liver
    ce_liver_perpetrator <- liver_ew_perpetrator / ve_liver
    ci_liver_perpetrator <- liver_iw_perpetrator / vi_liver
    cv_lnode_perpetrator <- lnode_vas_perpetrator / vv_lnode
    ce_lnode_perpetrator <- lnode_ew_perpetrator / ve_lnode
    ci_lnode_perpetrator <- lnode_iw_perpetrator / vi_lnode
    cv_other_perpetrator <- other_vas_perpetrator / vv_other
    ce_other_perpetrator <- other_ew_perpetrator / ve_other
    ci_other_perpetrator <- other_iw_perpetrator / vi_other
    cent_duodenum_perpetrator <- duodenum_enterocyte_perpetrator / vent_duodenum
    cent_jejunum_perpetrator <- jejunum_enterocyte_perpetrator / vent_jejunum
    cent_ileum_perpetrator <- ileum_enterocyte_perpetrator / vent_ileum
    cent_colon_perpetrator <- colon_enterocyte_perpetrator / vent_colon
    cli_um_perpetrator <- ci_liver_perpetrator * 1000 / mw_perpetrator # total intracellular liver, uM
    cul_perpetrator <- cli_um_perpetrator * fucel_liver_perpetrator # unbound, uM
    cdu_um_perpetrator <- cent_duodenum_perpetrator * 1000 / mw_perpetrator
    cudu_perpetrator <- cdu_um_perpetrator * fucel_gut_perpetrator
    cje_um_perpetrator <- cent_jejunum_perpetrator * 1000 / mw_perpetrator
    cuje_perpetrator <- cje_um_perpetrator * fucel_gut_perpetrator
    cil_um_perpetrator <- cent_ileum_perpetrator * 1000 / mw_perpetrator
    cuil_perpetrator <- cil_um_perpetrator * fucel_gut_perpetrator
    c_venous_perpetrator2 <- venous_perpetrator2 / v_venous
    c_arterial_perpetrator2 <- arterial_perpetrator2 / v_arterial
    cv_lung_perpetrator2 <- lung_vas_perpetrator2 / vv_lung
    ce_lung_perpetrator2 <- lung_ew_perpetrator2 / ve_lung
    ci_lung_perpetrator2 <- lung_iw_perpetrator2 / vi_lung
    cv_adipose_perpetrator2 <- adipose_vas_perpetrator2 / vv_adipose
    ce_adipose_perpetrator2 <- adipose_ew_perpetrator2 / ve_adipose
    ci_adipose_perpetrator2 <- adipose_iw_perpetrator2 / vi_adipose
    cv_bone_perpetrator2 <- bone_vas_perpetrator2 / vv_bone
    ce_bone_perpetrator2 <- bone_ew_perpetrator2 / ve_bone
    ci_bone_perpetrator2 <- bone_iw_perpetrator2 / vi_bone
    cv_brain_perpetrator2 <- brain_vas_perpetrator2 / vv_brain
    ce_brain_perpetrator2 <- brain_ew_perpetrator2 / ve_brain
    ci_brain_perpetrator2 <- brain_iw_perpetrator2 / vi_brain
    cv_gonads_perpetrator2 <- gonads_vas_perpetrator2 / vv_gonads
    ce_gonads_perpetrator2 <- gonads_ew_perpetrator2 / ve_gonads
    ci_gonads_perpetrator2 <- gonads_iw_perpetrator2 / vi_gonads
    cv_heart_perpetrator2 <- heart_vas_perpetrator2 / vv_heart
    ce_heart_perpetrator2 <- heart_ew_perpetrator2 / ve_heart
    ci_heart_perpetrator2 <- heart_iw_perpetrator2 / vi_heart
    cv_kidney_perpetrator2 <- kidney_vas_perpetrator2 / vv_kidney
    ce_kidney_perpetrator2 <- kidney_ew_perpetrator2 / ve_kidney
    ci_kidney_perpetrator2 <- kidney_iw_perpetrator2 / vi_kidney
    cv_muscle_perpetrator2 <- muscle_vas_perpetrator2 / vv_muscle
    ce_muscle_perpetrator2 <- muscle_ew_perpetrator2 / ve_muscle
    ci_muscle_perpetrator2 <- muscle_iw_perpetrator2 / vi_muscle
    cv_skin_perpetrator2 <- skin_vas_perpetrator2 / vv_skin
    ce_skin_perpetrator2 <- skin_ew_perpetrator2 / ve_skin
    ci_skin_perpetrator2 <- skin_iw_perpetrator2 / vi_skin
    cv_thymus_perpetrator2 <- thymus_vas_perpetrator2 / vv_thymus
    ce_thymus_perpetrator2 <- thymus_ew_perpetrator2 / ve_thymus
    ci_thymus_perpetrator2 <- thymus_iw_perpetrator2 / vi_thymus
    cv_gut_perpetrator2 <- gut_vas_perpetrator2 / vv_gut
    ce_gut_perpetrator2 <- gut_ew_perpetrator2 / ve_gut
    cv_spleen_perpetrator2 <- spleen_vas_perpetrator2 / vv_spleen
    ce_spleen_perpetrator2 <- spleen_ew_perpetrator2 / ve_spleen
    ci_spleen_perpetrator2 <- spleen_iw_perpetrator2 / vi_spleen
    cv_pancreas_perpetrator2 <- pancreas_vas_perpetrator2 / vv_pancreas
    ce_pancreas_perpetrator2 <- pancreas_ew_perpetrator2 / ve_pancreas
    ci_pancreas_perpetrator2 <- pancreas_iw_perpetrator2 / vi_pancreas
    cv_liver_perpetrator2 <- liver_vas_perpetrator2 / vv_liver
    ce_liver_perpetrator2 <- liver_ew_perpetrator2 / ve_liver
    ci_liver_perpetrator2 <- liver_iw_perpetrator2 / vi_liver
    cv_lnode_perpetrator2 <- lnode_vas_perpetrator2 / vv_lnode
    ce_lnode_perpetrator2 <- lnode_ew_perpetrator2 / ve_lnode
    ci_lnode_perpetrator2 <- lnode_iw_perpetrator2 / vi_lnode
    cv_other_perpetrator2 <- other_vas_perpetrator2 / vv_other
    ce_other_perpetrator2 <- other_ew_perpetrator2 / ve_other
    ci_other_perpetrator2 <- other_iw_perpetrator2 / vi_other
    cent_duodenum_perpetrator2 <- duodenum_enterocyte_perpetrator2 / vent_duodenum
    cent_jejunum_perpetrator2 <- jejunum_enterocyte_perpetrator2 / vent_jejunum
    cent_ileum_perpetrator2 <- ileum_enterocyte_perpetrator2 / vent_ileum
    cent_colon_perpetrator2 <- colon_enterocyte_perpetrator2 / vent_colon
    cli_um_perpetrator2 <- ci_liver_perpetrator2 * 1000 / mw_perpetrator2 # total intracellular liver, uM
    cul_perpetrator2 <- cli_um_perpetrator2 * fucel_liver_perpetrator2 # unbound, uM
    cdu_um_perpetrator2 <- cent_duodenum_perpetrator2 * 1000 / mw_perpetrator2
    cudu_perpetrator2 <- cdu_um_perpetrator2 * fucel_gut_perpetrator2
    cje_um_perpetrator2 <- cent_jejunum_perpetrator2 * 1000 / mw_perpetrator2
    cuje_perpetrator2 <- cje_um_perpetrator2 * fucel_gut_perpetrator2
    cil_um_perpetrator2 <- cent_ileum_perpetrator2 * 1000 / mw_perpetrator2
    cuil_perpetrator2 <- cil_um_perpetrator2 * fucel_gut_perpetrator2
    # --- interactions acting on the victim block's enzymes (PBPK_ODE_solution.m) ---
    xq_perpetrator <- (ddi_on * cli_um_perpetrator + (1 - ddi_on) * cli_um) * fucel_liver_perpetrator # concentration of drug perpetrator seen by the victim (uM x fucel)
    xq_perpetrator2 <- (ddi_on * cli_um_perpetrator2 + (1 - ddi_on) * cli_um) * fucel_liver_perpetrator2 # concentration of drug perpetrator2 seen by the victim (uM x fucel)
    ind_cyp3a4 <- (indmax_cyp3a4 - 1) * cli_um * fucel_liver / (ic50_cyp3a4 + cli_um * fucel_liver) +
      on_perpetrator * (ddi_on * indmax_cyp3a4_perpetrator - 1) * xq_perpetrator / (ddi_on * ic50_cyp3a4_perpetrator + (1 - ddi_on) + xq_perpetrator) +
      on_perpetrator2 * (ddi_on * indmax_cyp3a4_perpetrator2 - 1) * xq_perpetrator2 / (ddi_on * ic50_cyp3a4_perpetrator2 + (1 - ddi_on) + xq_perpetrator2)
    ind_cyp2c19 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    ind_cyp2d6 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    ind_cyp2c8 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    ind_cyp1a2 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    ind_cyp2a6 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    ind_cyp2b6 <- (indmax_cyp2b6 - 1) * cli_um * fucel_liver / (ic50_cyp2b6 + cli_um * fucel_liver) +
      on_perpetrator * (ddi_on * indmax_cyp2b6_perpetrator - 1) * xq_perpetrator / (ddi_on * ic50_cyp2b6_perpetrator + (1 - ddi_on) + xq_perpetrator) +
      on_perpetrator2 * (ddi_on * indmax_cyp2b6_perpetrator2 - 1) * xq_perpetrator2 / (ddi_on * ic50_cyp2b6_perpetrator2 + (1 - ddi_on) + xq_perpetrator2)
    ind_cyp2j2 <- -cli_um * fucel_liver / (1 + cli_um * fucel_liver) -
      on_perpetrator * xq_perpetrator / (1 + xq_perpetrator) -
      on_perpetrator2 * xq_perpetrator2 / (1 + xq_perpetrator2)
    mbi_cyp3a4 <- kinact_cyp3a4 * cli_um * fucel_liver / (kapp_cyp3a4 + cli_um * fucel_liver) +
      on_perpetrator * ddi_on * kinact_cyp3a4_perpetrator * xq_perpetrator / (kapp_cyp3a4_perpetrator + xq_perpetrator) +
      on_perpetrator2 * ddi_on * kinact_cyp3a4_perpetrator2 * xq_perpetrator2 / (kapp_cyp3a4_perpetrator2 + xq_perpetrator2)
    mbi_cyp2j2 <- kinact_cyp2j2 * cli_um * fucel_liver / (kapp_cyp2j2 + cli_um * fucel_liver) +
      on_perpetrator * ddi_on * kinact_cyp2j2_perpetrator * xq_perpetrator / (kapp_cyp2j2_perpetrator + xq_perpetrator) +
      on_perpetrator2 * ddi_on * kinact_cyp2j2_perpetrator2 * xq_perpetrator2 / (kapp_cyp2j2_perpetrator2 + xq_perpetrator2)
    com_cyp3a4 <- 1 + kiinv_cyp3a4 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_cyp3a4_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp3a4_perpetrator2 * xq_perpetrator2
    com_cyp2c19 <- 1 + kiinv_cyp2c19 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_cyp2c19_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2c19_perpetrator2 * xq_perpetrator2
    com_cyp2d6 <- 1 + kiinv_cyp2d6 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_cyp2d6_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2d6_perpetrator2 * xq_perpetrator2
    com_cyp2c8 <- 1 + kiinv_cyp2c8 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_cyp2c8_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2c8_perpetrator2 * xq_perpetrator2
    com_cyp1a2 <- 1 + kiinv_cyp1a2 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_cyp1a2_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp1a2_perpetrator2 * xq_perpetrator2
    indu_ugt1a1 <- indmax_ugt1a1 * cli_um * fucel_liver / (ic50_ugt1a1 + cli_um * fucel_liver) +
      on_perpetrator * ddi_on * indmax_ugt1a1_perpetrator * xq_perpetrator / (ic50_ugt1a1_perpetrator + xq_perpetrator) +
      on_perpetrator2 * ddi_on * indmax_ugt1a1_perpetrator2 * xq_perpetrator2 / (ic50_ugt1a1_perpetrator2 + xq_perpetrator2)
    com_ugt1a1 <- 1 + kiinv_ugt1a1 * cli_um * fucel_liver +
      on_perpetrator * ddi_on * kiinv_ugt1a1_perpetrator * xq_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_ugt1a1_perpetrator2 * xq_perpetrator2
    indg_du <- indmax_cyp3a4 * cudu / (ic50_cyp3a4 + cudu) +
      on_perpetrator * ddi_on * indmax_cyp3a4_perpetrator * cudu_perpetrator / (ic50_cyp3a4_perpetrator + cudu_perpetrator) +
      on_perpetrator2 * ddi_on * indmax_cyp3a4_perpetrator2 * cudu_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cudu_perpetrator2)
    mbig_du <- kinact_cyp3a4 * cli_um * fucel_gut / (kapp_cyp3a4 + cli_um * fucel_gut) +
      on_perpetrator * ddi_on * kinact_cyp3a4_perpetrator * cudu_perpetrator / (kapp_cyp3a4_perpetrator + cudu_perpetrator) +
      on_perpetrator2 * ddi_on * kinact_cyp3a4_perpetrator2 * cudu_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cudu_perpetrator2)
    comg_cyp3a4_du <- 1 + kiinv_cyp3a4 * cudu +
      on_perpetrator * ddi_on * kiinv_cyp3a4_perpetrator * cudu_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp3a4_perpetrator2 * cudu_perpetrator2
    comg_cyp2c19_du <- 1 + kiinv_cyp2c19 * cudu +
      on_perpetrator * ddi_on * kiinv_cyp2c19_perpetrator * cudu_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2c19_perpetrator2 * cudu_perpetrator2
    comg_cyp2d6_du <- 1 + kiinv_cyp2d6 * cudu +
      on_perpetrator * ddi_on * kiinv_cyp2d6_perpetrator * cudu_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2d6_perpetrator2 * cudu_perpetrator2
    indg_je <- indmax_cyp3a4 * cuje / (ic50_cyp3a4 + cuje) +
      on_perpetrator * ddi_on * indmax_cyp3a4_perpetrator * cuje_perpetrator / (ic50_cyp3a4_perpetrator + cuje_perpetrator) +
      on_perpetrator2 * ddi_on * indmax_cyp3a4_perpetrator2 * cuje_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cuje_perpetrator2)
    mbig_je <- kinact_cyp3a4 * cli_um * fucel_gut / (kapp_cyp3a4 + cli_um * fucel_gut) +
      on_perpetrator * ddi_on * kinact_cyp3a4_perpetrator * cuje_perpetrator / (kapp_cyp3a4_perpetrator + cuje_perpetrator) +
      on_perpetrator2 * ddi_on * kinact_cyp3a4_perpetrator2 * cuje_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cuje_perpetrator2)
    comg_cyp3a4_je <- 1 + kiinv_cyp3a4 * cuje +
      on_perpetrator * ddi_on * kiinv_cyp3a4_perpetrator * cuje_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp3a4_perpetrator2 * cuje_perpetrator2
    comg_cyp2c19_je <- 1 + kiinv_cyp2c19 * cuje +
      on_perpetrator * ddi_on * kiinv_cyp2c19_perpetrator * cuje_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2c19_perpetrator2 * cuje_perpetrator2
    comg_cyp2d6_je <- 1 + kiinv_cyp2d6 * cuje +
      on_perpetrator * ddi_on * kiinv_cyp2d6_perpetrator * cuje_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2d6_perpetrator2 * cuje_perpetrator2
    indg_il <- indmax_cyp3a4 * cuil / (ic50_cyp3a4 + cuil) +
      on_perpetrator * ddi_on * indmax_cyp3a4_perpetrator * cuil_perpetrator / (ic50_cyp3a4_perpetrator + cuil_perpetrator) +
      on_perpetrator2 * ddi_on * indmax_cyp3a4_perpetrator2 * cuil_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cuil_perpetrator2)
    mbig_il <- kinact_cyp3a4 * cli_um * fucel_gut / (kapp_cyp3a4 + cli_um * fucel_gut) +
      on_perpetrator * ddi_on * kinact_cyp3a4_perpetrator * cuil_perpetrator / (kapp_cyp3a4_perpetrator + cuil_perpetrator) +
      on_perpetrator2 * ddi_on * kinact_cyp3a4_perpetrator2 * cuil_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cuil_perpetrator2)
    comg_cyp3a4_il <- 1 + kiinv_cyp3a4 * cuil +
      on_perpetrator * ddi_on * kiinv_cyp3a4_perpetrator * cuil_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp3a4_perpetrator2 * cuil_perpetrator2
    comg_cyp2c19_il <- 1 + kiinv_cyp2c19 * cuil +
      on_perpetrator * ddi_on * kiinv_cyp2c19_perpetrator * cuil_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2c19_perpetrator2 * cuil_perpetrator2
    comg_cyp2d6_il <- 1 + kiinv_cyp2d6 * cuil +
      on_perpetrator * ddi_on * kiinv_cyp2d6_perpetrator * cuil_perpetrator +
      on_perpetrator2 * ddi_on * kiinv_cyp2d6_perpetrator2 * cuil_perpetrator2
    # --- interactions acting on the perpetrator block's enzymes (PBPK_ODE_solution.m) ---
    ind_cyp3a4_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) +
      (indmax_cyp3a4_perpetrator - 1) * cli_um_perpetrator * fucel_liver_perpetrator / (ic50_cyp3a4_perpetrator + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2c19_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2d6_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2c8_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp1a2_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2a6_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2b6_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) +
      (indmax_cyp2b6_perpetrator - 1) * cli_um_perpetrator * fucel_liver_perpetrator / (ic50_cyp2b6_perpetrator + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    ind_cyp2j2_perpetrator <- -1 * cli_um_perpetrator * fucel_liver / (1 + cli_um_perpetrator * fucel_liver) -
      cli_um_perpetrator * fucel_liver_perpetrator / (1 + cli_um_perpetrator * fucel_liver_perpetrator) -
      on_perpetrator2 * cli_um_perpetrator * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator * fucel_liver_perpetrator2)
    mbi_cyp3a4_perpetrator <- kinact_cyp3a4_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator / (kapp_cyp3a4_perpetrator + cli_um_perpetrator * fucel_liver_perpetrator)
    mbi_cyp2j2_perpetrator <- kinact_cyp2j2_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator / (kapp_cyp2j2_perpetrator + cli_um_perpetrator * fucel_liver_perpetrator)
    com_cyp3a4_perpetrator <- 1 + kiinv_cyp3a4_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    com_cyp2c19_perpetrator <- 1 + kiinv_cyp2c19_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    com_cyp2d6_perpetrator <- 1 + kiinv_cyp2d6_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    com_cyp2c8_perpetrator <- 1 + kiinv_cyp2c8_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    com_cyp1a2_perpetrator <- 1 + kiinv_cyp1a2_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    indu_ugt1a1_perpetrator <- indmax_ugt1a1_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator / (ic50_ugt1a1_perpetrator + cli_um_perpetrator * fucel_liver_perpetrator)
    com_ugt1a1_perpetrator <- 1 + kiinv_ugt1a1_perpetrator * cli_um_perpetrator * fucel_liver_perpetrator
    indg_du_perpetrator <- indmax_cyp3a4_perpetrator * cudu_perpetrator / (ic50_cyp3a4_perpetrator + cudu_perpetrator)
    mbig_du_perpetrator <- kinact_cyp3a4_perpetrator * cli_um_perpetrator * fucel_gut_perpetrator / (kapp_cyp3a4_perpetrator + cli_um_perpetrator * fucel_gut_perpetrator)
    comg_cyp3a4_du_perpetrator <- 1 + kiinv_cyp3a4_perpetrator * cudu_perpetrator
    comg_cyp2c19_du_perpetrator <- 1 + kiinv_cyp2c19_perpetrator * cudu_perpetrator
    comg_cyp2d6_du_perpetrator <- 1 + kiinv_cyp2d6_perpetrator * cudu_perpetrator
    indg_je_perpetrator <- indmax_cyp3a4_perpetrator * cuje_perpetrator / (ic50_cyp3a4_perpetrator + cuje_perpetrator)
    mbig_je_perpetrator <- kinact_cyp3a4_perpetrator * cli_um_perpetrator * fucel_gut_perpetrator / (kapp_cyp3a4_perpetrator + cli_um_perpetrator * fucel_gut_perpetrator)
    comg_cyp3a4_je_perpetrator <- 1 + kiinv_cyp3a4_perpetrator * cuje_perpetrator
    comg_cyp2c19_je_perpetrator <- 1 + kiinv_cyp2c19_perpetrator * cuje_perpetrator
    comg_cyp2d6_je_perpetrator <- 1 + kiinv_cyp2d6_perpetrator * cuje_perpetrator
    indg_il_perpetrator <- indmax_cyp3a4_perpetrator * cuil_perpetrator / (ic50_cyp3a4_perpetrator + cuil_perpetrator)
    mbig_il_perpetrator <- kinact_cyp3a4_perpetrator * cli_um_perpetrator * fucel_gut_perpetrator / (kapp_cyp3a4_perpetrator + cli_um_perpetrator * fucel_gut_perpetrator)
    comg_cyp3a4_il_perpetrator <- 1 + kiinv_cyp3a4_perpetrator * cuil_perpetrator
    comg_cyp2c19_il_perpetrator <- 1 + kiinv_cyp2c19_perpetrator * cuil_perpetrator
    comg_cyp2d6_il_perpetrator <- 1 + kiinv_cyp2d6_perpetrator * cuil_perpetrator
    # --- interactions acting on the perpetrator2 block's enzymes (PBPK_ODE_solution.m) ---
    ind_cyp3a4_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) +
      (indmax_cyp3a4_perpetrator2 - 1) * cli_um_perpetrator2 * fucel_liver_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2c19_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2d6_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2c8_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp1a2_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2a6_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2b6_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) +
      (indmax_cyp2b6_perpetrator2 - 1) * cli_um_perpetrator2 * fucel_liver_perpetrator2 / (ic50_cyp2b6_perpetrator2 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    ind_cyp2j2_perpetrator2 <- -1 * cli_um_perpetrator2 * fucel_liver / (1 + cli_um_perpetrator2 * fucel_liver) -
      on_perpetrator * cli_um_perpetrator2 * fucel_liver_perpetrator / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator) -
      cli_um_perpetrator2 * fucel_liver_perpetrator2 / (1 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    mbi_cyp3a4_perpetrator2 <- kinact_cyp3a4_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    mbi_cyp2j2_perpetrator2 <- kinact_cyp2j2_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2 / (kapp_cyp2j2_perpetrator2 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    com_cyp3a4_perpetrator2 <- 1 + kiinv_cyp3a4_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    com_cyp2c19_perpetrator2 <- 1 + kiinv_cyp2c19_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    com_cyp2d6_perpetrator2 <- 1 + kiinv_cyp2d6_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    com_cyp2c8_perpetrator2 <- 1 + kiinv_cyp2c8_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    com_cyp1a2_perpetrator2 <- 1 + kiinv_cyp1a2_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    indu_ugt1a1_perpetrator2 <- indmax_ugt1a1_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2 / (ic50_ugt1a1_perpetrator2 + cli_um_perpetrator2 * fucel_liver_perpetrator2)
    com_ugt1a1_perpetrator2 <- 1 + kiinv_ugt1a1_perpetrator2 * cli_um_perpetrator2 * fucel_liver_perpetrator2
    indg_du_perpetrator2 <- indmax_cyp3a4_perpetrator2 * cudu_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cudu_perpetrator2)
    mbig_du_perpetrator2 <- kinact_cyp3a4_perpetrator2 * cli_um_perpetrator2 * fucel_gut_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cli_um_perpetrator2 * fucel_gut_perpetrator2)
    comg_cyp3a4_du_perpetrator2 <- 1 + kiinv_cyp3a4_perpetrator2 * cudu_perpetrator2
    comg_cyp2c19_du_perpetrator2 <- 1 + kiinv_cyp2c19_perpetrator2 * cudu_perpetrator2
    comg_cyp2d6_du_perpetrator2 <- 1 + kiinv_cyp2d6_perpetrator2 * cudu_perpetrator2
    indg_je_perpetrator2 <- indmax_cyp3a4_perpetrator2 * cuje_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cuje_perpetrator2)
    mbig_je_perpetrator2 <- kinact_cyp3a4_perpetrator2 * cli_um_perpetrator2 * fucel_gut_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cli_um_perpetrator2 * fucel_gut_perpetrator2)
    comg_cyp3a4_je_perpetrator2 <- 1 + kiinv_cyp3a4_perpetrator2 * cuje_perpetrator2
    comg_cyp2c19_je_perpetrator2 <- 1 + kiinv_cyp2c19_perpetrator2 * cuje_perpetrator2
    comg_cyp2d6_je_perpetrator2 <- 1 + kiinv_cyp2d6_perpetrator2 * cuje_perpetrator2
    indg_il_perpetrator2 <- indmax_cyp3a4_perpetrator2 * cuil_perpetrator2 / (ic50_cyp3a4_perpetrator2 + cuil_perpetrator2)
    mbig_il_perpetrator2 <- kinact_cyp3a4_perpetrator2 * cli_um_perpetrator2 * fucel_gut_perpetrator2 / (kapp_cyp3a4_perpetrator2 + cli_um_perpetrator2 * fucel_gut_perpetrator2)
    comg_cyp3a4_il_perpetrator2 <- 1 + kiinv_cyp3a4_perpetrator2 * cuil_perpetrator2
    comg_cyp2c19_il_perpetrator2 <- 1 + kiinv_cyp2c19_perpetrator2 * cuil_perpetrator2
    comg_cyp2d6_il_perpetrator2 <- 1 + kiinv_cyp2d6_perpetrator2 * cuil_perpetrator2
    # --- victim block metabolism: Vmax/(Km*Com + Cu) + CLint/Com per enzyme (L/h/pmol) ---
    r_cyp3a4 <- (vmax1_cyp3a4 / (km1_cyp3a4 * com_cyp3a4 + cul) + vmax2_cyp3a4 / (km2_cyp3a4 * com_cyp3a4 + cul) +
      clint_cyp3a4 / com_cyp3a4) * 60e-6 * ab_cyp3a4_liver * enzyme_cyp3a4_liver
    r_cyp2c19 <- (vmax1_cyp2c19 / (km1_cyp2c19 * com_cyp2c19 + cul) + clint_cyp2c19 / com_cyp2c19) * 60e-6 * ab_cyp2c19_liver * enzyme_cyp2c19_liver
    r_cyp2d6 <- (vmax1_cyp2d6 / (km1_cyp2d6 * com_cyp2d6 + cul) + clint_cyp2d6 / com_cyp2d6) * 60e-6 * ab_cyp2d6_liver * enzyme_cyp2d6_liver
    r_cyp2c8 <- clint_cyp2c8 / com_cyp2c8 * 60e-6 * ab_cyp2c8_liver * enzyme_cyp2c8_liver
    r_cyp1a2 <- clint_cyp1a2 / com_cyp1a2 * 60e-6 * ab_cyp1a2_liver * enzyme_cyp1a2_liver
    r_cyp2a6 <- clint_cyp2a6 * 60e-6 * ab_cyp2a6_liver * enzyme_cyp2a6_liver
    r_cyp2b6 <- clint_cyp2b6 * 60e-6 * ab_cyp2b6_liver * enzyme_cyp2b6_liver
    r_cyp2j2 <- clint_cyp2j2 * 60e-6 * ab_cyp2j2_liver * enzyme_cyp2j2_liver
    r_ugt1a1 <- clint_ugt1a1 / com_ugt1a1 * 60e-6 * ab_ugt1a1_liver * enzyme_ugt1a1_liver
    r_ugt1a4 <- vmax_ugt1a4 / (km_ugt1a4 + cul) * 60e-6 * ab_ugt1a4_liver
    clmet_liver <- (r_cyp3a4 + r_cyp2c19 + r_cyp2d6 + r_cyp2c8 + r_cyp1a2 + r_cyp2a6 +
      r_cyp2b6 + r_cyp2j2 + r_ugt1a1 + r_ugt1a4 + clint_hep * 60e-6) * mppgl * (w_liver * 1000) + clbile
    clmet_duodenum <- ((vmax1_cyp3a4 / (km1_cyp3a4 * comg_cyp3a4_du + cudu) + vmax2_cyp3a4 / (km2_cyp3a4 * comg_cyp3a4_du + cudu) +
      clint_cyp3a4 / comg_cyp3a4_du) * abg_cyp3a4_duodenum * enzyme_cyp3a4_duodenum +
      (vmax1_cyp2c19 / (km1_cyp2c19 * comg_cyp2c19_du + cudu) + clint_cyp2c19 / comg_cyp2c19_du) * abg_cyp2c19_duodenum +
      (vmax1_cyp2d6 / (km1_cyp2d6 * comg_cyp2d6_du + cudu) + clint_cyp2d6 / comg_cyp2d6_du) * abg_cyp2d6_duodenum) * 60e-6
    clmet_jejunum <- ((vmax1_cyp3a4 / (km1_cyp3a4 * comg_cyp3a4_je + cuje) + vmax2_cyp3a4 / (km2_cyp3a4 * comg_cyp3a4_je + cuje) +
      clint_cyp3a4 / comg_cyp3a4_je) * abg_cyp3a4_jejunum * enzyme_cyp3a4_jejunum +
      (vmax1_cyp2c19 / (km1_cyp2c19 * comg_cyp2c19_je + cuje) + clint_cyp2c19 / comg_cyp2c19_je) * abg_cyp2c19_jejunum +
      (vmax1_cyp2d6 / (km1_cyp2d6 * comg_cyp2d6_je + cuje) + clint_cyp2d6 / comg_cyp2d6_je) * abg_cyp2d6_jejunum) * 60e-6
    clmet_ileum <- ((vmax1_cyp3a4 / (km1_cyp3a4 * comg_cyp3a4_il + cuil) + vmax2_cyp3a4 / (km2_cyp3a4 * comg_cyp3a4_il + cuil) +
      clint_cyp3a4 / comg_cyp3a4_il) * abg_cyp3a4_ileum * enzyme_cyp3a4_ileum +
      (vmax1_cyp2c19 / (km1_cyp2c19 * comg_cyp2c19_il + cuil) + clint_cyp2c19 / comg_cyp2c19_il) * abg_cyp2c19_ileum +
      (vmax1_cyp2d6 / (km1_cyp2d6 * comg_cyp2d6_il + cuil) + clint_cyp2d6 / comg_cyp2d6_il) * abg_cyp2d6_ileum) * 60e-6
    # --- perpetrator block metabolism: Vmax/(Km*Com + Cu) + CLint/Com per enzyme (L/h/pmol) ---
    r_cyp3a4_perpetrator <- (vmax1_cyp3a4_perpetrator / (km1_cyp3a4_perpetrator * com_cyp3a4_perpetrator + cul_perpetrator) + vmax2_cyp3a4_perpetrator / (km2_cyp3a4_perpetrator * com_cyp3a4_perpetrator + cul_perpetrator) +
      clint_cyp3a4_perpetrator / com_cyp3a4_perpetrator) * 60e-6 * ab_cyp3a4_liver * enzyme_cyp3a4_liver_perpetrator
    r_cyp2c19_perpetrator <- (vmax1_cyp2c19_perpetrator / (km1_cyp2c19_perpetrator * com_cyp2c19_perpetrator + cul_perpetrator) + clint_cyp2c19_perpetrator / com_cyp2c19_perpetrator) * 60e-6 * ab_cyp2c19_liver * enzyme_cyp2c19_liver_perpetrator
    r_cyp2d6_perpetrator <- (vmax1_cyp2d6_perpetrator / (km1_cyp2d6_perpetrator * com_cyp2d6_perpetrator + cul_perpetrator) + clint_cyp2d6_perpetrator / com_cyp2d6_perpetrator) * 60e-6 * ab_cyp2d6_liver * enzyme_cyp2d6_liver_perpetrator
    r_cyp2c8_perpetrator <- clint_cyp2c8_perpetrator / com_cyp2c8_perpetrator * 60e-6 * ab_cyp2c8_liver * enzyme_cyp2c8_liver_perpetrator
    r_cyp1a2_perpetrator <- clint_cyp1a2_perpetrator / com_cyp1a2_perpetrator * 60e-6 * ab_cyp1a2_liver * enzyme_cyp1a2_liver_perpetrator
    r_cyp2a6_perpetrator <- clint_cyp2a6_perpetrator * 60e-6 * ab_cyp2a6_liver * enzyme_cyp2a6_liver_perpetrator
    r_cyp2b6_perpetrator <- clint_cyp2b6_perpetrator * 60e-6 * ab_cyp2b6_liver * enzyme_cyp2b6_liver_perpetrator
    r_cyp2j2_perpetrator <- clint_cyp2j2_perpetrator * 60e-6 * ab_cyp2j2_liver * enzyme_cyp2j2_liver_perpetrator
    r_ugt1a1_perpetrator <- clint_ugt1a1_perpetrator / com_ugt1a1_perpetrator * 60e-6 * ab_ugt1a1_liver * enzyme_ugt1a1_liver_perpetrator
    r_ugt1a4_perpetrator <- vmax_ugt1a4_perpetrator / (km_ugt1a4_perpetrator + cul_perpetrator) * 60e-6 * ab_ugt1a4_liver
    clmet_liver_perpetrator <- (r_cyp3a4_perpetrator + r_cyp2c19_perpetrator + r_cyp2d6_perpetrator + r_cyp2c8_perpetrator + r_cyp1a2_perpetrator + r_cyp2a6_perpetrator +
      r_cyp2b6_perpetrator + r_cyp2j2_perpetrator + r_ugt1a1_perpetrator + r_ugt1a4_perpetrator + clint_hep_perpetrator * 60e-6) * mppgl * (w_liver * 1000) + clbile_perpetrator
    clmet_duodenum_perpetrator <- ((vmax1_cyp3a4_perpetrator / (km1_cyp3a4_perpetrator * comg_cyp3a4_du_perpetrator + cudu_perpetrator) + vmax2_cyp3a4_perpetrator / (km2_cyp3a4_perpetrator * comg_cyp3a4_du_perpetrator + cudu_perpetrator) +
      clint_cyp3a4_perpetrator / comg_cyp3a4_du_perpetrator) * abg_cyp3a4_duodenum * enzyme_cyp3a4_duodenum_perpetrator +
      (vmax1_cyp2c19_perpetrator / (km1_cyp2c19_perpetrator * comg_cyp2c19_du_perpetrator + cudu_perpetrator) + clint_cyp2c19_perpetrator / comg_cyp2c19_du_perpetrator) * abg_cyp2c19_duodenum +
      (vmax1_cyp2d6_perpetrator / (km1_cyp2d6_perpetrator * comg_cyp2d6_du_perpetrator + cudu_perpetrator) + clint_cyp2d6_perpetrator / comg_cyp2d6_du_perpetrator) * abg_cyp2d6_duodenum) * 60e-6
    clmet_jejunum_perpetrator <- ((vmax1_cyp3a4_perpetrator / (km1_cyp3a4_perpetrator * comg_cyp3a4_je_perpetrator + cuje_perpetrator) + vmax2_cyp3a4_perpetrator / (km2_cyp3a4_perpetrator * comg_cyp3a4_je_perpetrator + cuje_perpetrator) +
      clint_cyp3a4_perpetrator / comg_cyp3a4_je_perpetrator) * abg_cyp3a4_jejunum * enzyme_cyp3a4_jejunum_perpetrator +
      (vmax1_cyp2c19_perpetrator / (km1_cyp2c19_perpetrator * comg_cyp2c19_je_perpetrator + cuje_perpetrator) + clint_cyp2c19_perpetrator / comg_cyp2c19_je_perpetrator) * abg_cyp2c19_jejunum +
      (vmax1_cyp2d6_perpetrator / (km1_cyp2d6_perpetrator * comg_cyp2d6_je_perpetrator + cuje_perpetrator) + clint_cyp2d6_perpetrator / comg_cyp2d6_je_perpetrator) * abg_cyp2d6_jejunum) * 60e-6
    clmet_ileum_perpetrator <- ((vmax1_cyp3a4_perpetrator / (km1_cyp3a4_perpetrator * comg_cyp3a4_il_perpetrator + cuil_perpetrator) + vmax2_cyp3a4_perpetrator / (km2_cyp3a4_perpetrator * comg_cyp3a4_il_perpetrator + cuil_perpetrator) +
      clint_cyp3a4_perpetrator / comg_cyp3a4_il_perpetrator) * abg_cyp3a4_ileum * enzyme_cyp3a4_ileum_perpetrator +
      (vmax1_cyp2c19_perpetrator / (km1_cyp2c19_perpetrator * comg_cyp2c19_il_perpetrator + cuil_perpetrator) + clint_cyp2c19_perpetrator / comg_cyp2c19_il_perpetrator) * abg_cyp2c19_ileum +
      (vmax1_cyp2d6_perpetrator / (km1_cyp2d6_perpetrator * comg_cyp2d6_il_perpetrator + cuil_perpetrator) + clint_cyp2d6_perpetrator / comg_cyp2d6_il_perpetrator) * abg_cyp2d6_ileum) * 60e-6
    # --- perpetrator2 block metabolism: Vmax/(Km*Com + Cu) + CLint/Com per enzyme (L/h/pmol) ---
    r_cyp3a4_perpetrator2 <- (vmax1_cyp3a4_perpetrator2 / (km1_cyp3a4_perpetrator2 * com_cyp3a4_perpetrator2 + cul_perpetrator2) + vmax2_cyp3a4_perpetrator2 / (km2_cyp3a4_perpetrator2 * com_cyp3a4_perpetrator2 + cul_perpetrator2) +
      clint_cyp3a4_perpetrator2 / com_cyp3a4_perpetrator2) * 60e-6 * ab_cyp3a4_liver * enzyme_cyp3a4_liver_perpetrator2
    r_cyp2c19_perpetrator2 <- (vmax1_cyp2c19_perpetrator2 / (km1_cyp2c19_perpetrator2 * com_cyp2c19_perpetrator2 + cul_perpetrator2) + clint_cyp2c19_perpetrator2 / com_cyp2c19_perpetrator2) * 60e-6 * ab_cyp2c19_liver * enzyme_cyp2c19_liver_perpetrator2
    r_cyp2d6_perpetrator2 <- (vmax1_cyp2d6_perpetrator2 / (km1_cyp2d6_perpetrator2 * com_cyp2d6_perpetrator2 + cul_perpetrator2) + clint_cyp2d6_perpetrator2 / com_cyp2d6_perpetrator2) * 60e-6 * ab_cyp2d6_liver * enzyme_cyp2d6_liver_perpetrator2
    r_cyp2c8_perpetrator2 <- clint_cyp2c8_perpetrator2 / com_cyp2c8_perpetrator2 * 60e-6 * ab_cyp2c8_liver * enzyme_cyp2c8_liver_perpetrator2
    r_cyp1a2_perpetrator2 <- clint_cyp1a2_perpetrator2 / com_cyp1a2_perpetrator2 * 60e-6 * ab_cyp1a2_liver * enzyme_cyp1a2_liver_perpetrator2
    r_cyp2a6_perpetrator2 <- clint_cyp2a6_perpetrator2 * 60e-6 * ab_cyp2a6_liver * enzyme_cyp2a6_liver_perpetrator2
    r_cyp2b6_perpetrator2 <- clint_cyp2b6_perpetrator2 * 60e-6 * ab_cyp2b6_liver * enzyme_cyp2b6_liver_perpetrator2
    r_cyp2j2_perpetrator2 <- clint_cyp2j2_perpetrator2 * 60e-6 * ab_cyp2j2_liver * enzyme_cyp2j2_liver_perpetrator2
    r_ugt1a1_perpetrator2 <- clint_ugt1a1_perpetrator2 / com_ugt1a1_perpetrator2 * 60e-6 * ab_ugt1a1_liver * enzyme_ugt1a1_liver_perpetrator2
    r_ugt1a4_perpetrator2 <- vmax_ugt1a4_perpetrator2 / (km_ugt1a4_perpetrator2 + cul_perpetrator2) * 60e-6 * ab_ugt1a4_liver
    clmet_liver_perpetrator2 <- (r_cyp3a4_perpetrator2 + r_cyp2c19_perpetrator2 + r_cyp2d6_perpetrator2 + r_cyp2c8_perpetrator2 + r_cyp1a2_perpetrator2 + r_cyp2a6_perpetrator2 +
      r_cyp2b6_perpetrator2 + r_cyp2j2_perpetrator2 + r_ugt1a1_perpetrator2 + r_ugt1a4_perpetrator2 + clint_hep_perpetrator2 * 60e-6) * mppgl * (w_liver * 1000) + clbile_perpetrator2
    clmet_duodenum_perpetrator2 <- ((vmax1_cyp3a4_perpetrator2 / (km1_cyp3a4_perpetrator2 * comg_cyp3a4_du_perpetrator2 + cudu_perpetrator2) + vmax2_cyp3a4_perpetrator2 / (km2_cyp3a4_perpetrator2 * comg_cyp3a4_du_perpetrator2 + cudu_perpetrator2) +
      clint_cyp3a4_perpetrator2 / comg_cyp3a4_du_perpetrator2) * abg_cyp3a4_duodenum * enzyme_cyp3a4_duodenum_perpetrator2 +
      (vmax1_cyp2c19_perpetrator2 / (km1_cyp2c19_perpetrator2 * comg_cyp2c19_du_perpetrator2 + cudu_perpetrator2) + clint_cyp2c19_perpetrator2 / comg_cyp2c19_du_perpetrator2) * abg_cyp2c19_duodenum +
      (vmax1_cyp2d6_perpetrator2 / (km1_cyp2d6_perpetrator2 * comg_cyp2d6_du_perpetrator2 + cudu_perpetrator2) + clint_cyp2d6_perpetrator2 / comg_cyp2d6_du_perpetrator2) * abg_cyp2d6_duodenum) * 60e-6
    clmet_jejunum_perpetrator2 <- ((vmax1_cyp3a4_perpetrator2 / (km1_cyp3a4_perpetrator2 * comg_cyp3a4_je_perpetrator2 + cuje_perpetrator2) + vmax2_cyp3a4_perpetrator2 / (km2_cyp3a4_perpetrator2 * comg_cyp3a4_je_perpetrator2 + cuje_perpetrator2) +
      clint_cyp3a4_perpetrator2 / comg_cyp3a4_je_perpetrator2) * abg_cyp3a4_jejunum * enzyme_cyp3a4_jejunum_perpetrator2 +
      (vmax1_cyp2c19_perpetrator2 / (km1_cyp2c19_perpetrator2 * comg_cyp2c19_je_perpetrator2 + cuje_perpetrator2) + clint_cyp2c19_perpetrator2 / comg_cyp2c19_je_perpetrator2) * abg_cyp2c19_jejunum +
      (vmax1_cyp2d6_perpetrator2 / (km1_cyp2d6_perpetrator2 * comg_cyp2d6_je_perpetrator2 + cuje_perpetrator2) + clint_cyp2d6_perpetrator2 / comg_cyp2d6_je_perpetrator2) * abg_cyp2d6_jejunum) * 60e-6
    clmet_ileum_perpetrator2 <- ((vmax1_cyp3a4_perpetrator2 / (km1_cyp3a4_perpetrator2 * comg_cyp3a4_il_perpetrator2 + cuil_perpetrator2) + vmax2_cyp3a4_perpetrator2 / (km2_cyp3a4_perpetrator2 * comg_cyp3a4_il_perpetrator2 + cuil_perpetrator2) +
      clint_cyp3a4_perpetrator2 / comg_cyp3a4_il_perpetrator2) * abg_cyp3a4_ileum * enzyme_cyp3a4_ileum_perpetrator2 +
      (vmax1_cyp2c19_perpetrator2 / (km1_cyp2c19_perpetrator2 * comg_cyp2c19_il_perpetrator2 + cuil_perpetrator2) + clint_cyp2c19_perpetrator2 / comg_cyp2c19_il_perpetrator2) * abg_cyp2c19_ileum +
      (vmax1_cyp2d6_perpetrator2 / (km1_cyp2d6_perpetrator2 * comg_cyp2d6_il_perpetrator2 + cuil_perpetrator2) + clint_cyp2d6_perpetrator2 / comg_cyp2d6_il_perpetrator2) * abg_cyp2d6_ileum) * 60e-6
    # --- victim block ODEs (PBPK_ODE_solution.m rhs_function, amount form) ---
    d/dt(lung_vas) <- q_lung * c_venous - ql_lung * cv_lung - p_lung * fin_all * jin_adipose * cv_lung / bp + (p_lung - l_lung) * ce_lung
    d/dt(lung_ew) <- p_lung * fin_all * jin_adipose * cv_lung / bp - (p_lung - l_lung) * ce_lung - l_lung * ce_lung -
      clin_lung * jin_all * ce_lung * fuint_lung + clin_lung * ci_lung * fucel_lung
    d/dt(lung_iw) <- clin_lung * jin_all * ce_lung * fuint_lung - clin_lung * ci_lung * fucel_lung
    d/dt(adipose_vas) <- q_adipose * c_arterial - ql_adipose * cv_adipose - p_adipose * fin_all * cv_adipose / bp + (p_adipose - l_adipose) * ce_adipose
    d/dt(adipose_ew) <- p_adipose * fin_all * cv_adipose / bp - (p_adipose - l_adipose) * ce_adipose - l_adipose * ce_adipose -
      clin_adipose * jin_all * jin_adipose * ce_adipose * fuint_adipose + clin_adipose * ci_adipose * fucel_adipose
    d/dt(adipose_iw) <- clin_adipose * jin_all * jin_adipose * ce_adipose * fuint_adipose - clin_adipose * ci_adipose * fucel_adipose
    d/dt(bone_vas) <- q_bone * c_arterial - ql_bone * cv_bone - p_bone * fin_all * cv_bone / bp + (p_bone - l_bone) * ce_bone
    d/dt(bone_ew) <- p_bone * fin_all * cv_bone / bp - (p_bone - l_bone) * ce_bone - l_bone * ce_bone -
      clin_bone * jin_all * ce_bone * fuint_bone + clin_bone * ci_bone * fucel_bone
    d/dt(bone_iw) <- clin_bone * jin_all * ce_bone * fuint_bone - clin_bone * ci_bone * fucel_bone
    d/dt(brain_vas) <- q_brain * c_arterial - ql_brain * cv_brain - p_brain * fin_all * cv_brain / bp + (p_brain - l_brain) * ce_brain
    d/dt(brain_ew) <- p_brain * fin_all * cv_brain / bp - (p_brain - l_brain) * ce_brain - l_brain * ce_brain -
      clin_brain * jin_all * ce_brain * fuint_brain + clin_brain * ci_brain * fucel_brain
    d/dt(brain_iw) <- clin_brain * jin_all * ce_brain * fuint_brain - clin_brain * ci_brain * fucel_brain
    d/dt(gonads_vas) <- q_gonads * c_arterial - ql_gonads * cv_gonads - p_gonads * fin_all * cv_gonads / bp + (p_gonads - l_gonads) * ce_gonads
    d/dt(gonads_ew) <- p_gonads * fin_all * cv_gonads / bp - (p_gonads - l_gonads) * ce_gonads - l_gonads * ce_gonads -
      clin_gonads * jin_all * ce_gonads * fuint_gonads + clin_gonads * ci_gonads * fucel_gonads
    d/dt(gonads_iw) <- clin_gonads * jin_all * ce_gonads * fuint_gonads - clin_gonads * ci_gonads * fucel_gonads
    d/dt(heart_vas) <- q_heart * c_arterial - ql_heart * cv_heart - p_heart * fin_all * cv_heart / bp + (p_heart - l_heart) * ce_heart
    d/dt(heart_ew) <- p_heart * fin_all * cv_heart / bp - (p_heart - l_heart) * ce_heart - l_heart * ce_heart -
      clin_heart * jin_all * ce_heart * fuint_heart + clin_heart * ci_heart * fucel_heart
    d/dt(heart_iw) <- clin_heart * jin_all * ce_heart * fuint_heart - clin_heart * ci_heart * fucel_heart
    d/dt(kidney_vas) <- q_kidney * c_arterial - ql_kidney * cv_kidney - p_kidney * fin_all * jin_muscle * cv_kidney / bp + (p_kidney - l_kidney) * ce_kidney - clr * cv_kidney
    d/dt(kidney_ew) <- p_kidney * fin_all * jin_muscle * cv_kidney / bp - (p_kidney - l_kidney) * ce_kidney - l_kidney * ce_kidney -
      clin_kidney * jin_all * ce_kidney * fuint_kidney + clin_kidney * ci_kidney * fucel_kidney
    d/dt(kidney_iw) <- clin_kidney * jin_all * ce_kidney * fuint_kidney - clin_kidney * ci_kidney * fucel_kidney
    d/dt(muscle_vas) <- q_muscle * c_arterial - ql_muscle * cv_muscle - p_muscle * fin_all * cv_muscle / bp + (p_muscle - l_muscle) * ce_muscle
    d/dt(muscle_ew) <- p_muscle * fin_all * cv_muscle / bp - (p_muscle - l_muscle) * ce_muscle - l_muscle * ce_muscle -
      clin_muscle * jin_all * jin_muscle * ce_muscle * fuint_muscle + clin_muscle * ci_muscle * fucel_muscle
    d/dt(muscle_iw) <- clin_muscle * jin_all * jin_muscle * ce_muscle * fuint_muscle - clin_muscle * ci_muscle * fucel_muscle
    d/dt(skin_vas) <- q_skin * c_arterial - ql_skin * cv_skin - p_skin * fin_all * cv_skin / bp + (p_skin - l_skin) * ce_skin
    d/dt(skin_ew) <- p_skin * fin_all * cv_skin / bp - (p_skin - l_skin) * ce_skin - l_skin * ce_skin -
      clin_skin * jin_all * ce_skin * fuint_skin + clin_skin * ci_skin * fucel_skin
    d/dt(skin_iw) <- clin_skin * jin_all * ce_skin * fuint_skin - clin_skin * ci_skin * fucel_skin
    d/dt(thymus_vas) <- q_thymus * c_arterial - ql_thymus * cv_thymus - p_thymus * fin_all * cv_thymus / bp + (p_thymus - l_thymus) * ce_thymus
    d/dt(thymus_ew) <- p_thymus * fin_all * cv_thymus / bp - (p_thymus - l_thymus) * ce_thymus - l_thymus * ce_thymus -
      clin_thymus * jin_all * ce_thymus * fuint_thymus + clin_thymus * ci_thymus * fucel_thymus
    d/dt(thymus_iw) <- clin_thymus * jin_all * ce_thymus * fuint_thymus - clin_thymus * ci_thymus * fucel_thymus
    d/dt(spleen_vas) <- q_spleen * c_arterial - ql_spleen * cv_spleen - p_spleen * fin_all * cv_spleen / bp + (p_spleen - l_spleen) * ce_spleen
    d/dt(spleen_ew) <- p_spleen * fin_all * cv_spleen / bp - (p_spleen - l_spleen) * ce_spleen - l_spleen * ce_spleen -
      clin_spleen * jin_all * ce_spleen * fuint_spleen + clin_spleen * ci_spleen * fucel_spleen
    d/dt(spleen_iw) <- clin_spleen * jin_all * ce_spleen * fuint_spleen - clin_spleen * ci_spleen * fucel_spleen
    d/dt(pancreas_vas) <- q_pancreas * c_arterial - ql_pancreas * cv_pancreas - p_pancreas * fin_all * jin_liver * cv_pancreas / bp + (p_pancreas - l_pancreas) * ce_pancreas
    d/dt(pancreas_ew) <- p_pancreas * fin_all * jin_liver * cv_pancreas / bp - (p_pancreas - l_pancreas) * ce_pancreas - l_pancreas * ce_pancreas -
      clin_pancreas * jin_all * ce_pancreas * fuint_pancreas + clin_pancreas * ci_pancreas * fucel_pancreas
    d/dt(pancreas_iw) <- clin_pancreas * jin_all * ce_pancreas * fuint_pancreas - clin_pancreas * ci_pancreas * fucel_pancreas
    d/dt(other_vas) <- q_other * c_arterial - ql_other * cv_other - p_other * fin_all * cv_other / bp + (p_other - l_other) * ce_other
    d/dt(other_ew) <- p_other * fin_all * cv_other / bp - (p_other - l_other) * ce_other - l_other * ce_other -
      clin_other * jin_all * ce_other * fuint_other + clin_other * ci_other * fucel_other
    d/dt(other_iw) <- clin_other * jin_all * ce_other * fuint_other - clin_other * ci_other * fucel_other
    d/dt(gut_vas) <- q_gut * c_arterial - ql_gut * cv_gut - p_gut * fin_all * cv_gut / bp + (p_gut - l_gut) * ce_gut
    d/dt(gut_ew) <- p_gut * fin_all * cv_gut / bp - (p_gut - l_gut) * ce_gut - l_gut * ce_gut +
      clin_gut * fgp * fucel_gut * (cent_duodenum + cent_jejunum + cent_ileum + cent_colon)
    d/dt(liver_vas) <- ql_gut * cv_gut + ql_spleen * cv_spleen + ql_pancreas * cv_pancreas + (q_ha + q_by) * c_arterial -
      ql_liver * cv_liver - p_liver * fin_all * cv_liver / bp + (p_liver - l_liver) * ce_liver
    d/dt(liver_ew) <- p_liver * fin_all * cv_liver / bp - (p_liver - l_liver) * ce_liver - l_liver * ce_liver -
      clin_liver * jin_all * jin_liver * ce_liver * fuint_liver + clin_liver * ci_liver * fucel_liver
    d/dt(liver_iw) <- clin_liver * jin_all * jin_liver * ce_liver * fuint_liver - clin_liver * ci_liver * fucel_liver -
      clmet_liver * ci_liver * fucel_liver
    d/dt(lnode_vas) <- q_lnode * c_arterial - q_lnode * cv_lnode + ltot * ce_lnode - ltot * cv_lnode
    d/dt(lnode_ew) <- l_lung * ce_lung + l_adipose * ce_adipose + l_bone * ce_bone + l_brain * ce_brain +
      l_gonads * ce_gonads + l_heart * ce_heart + l_kidney * ce_kidney + l_muscle * ce_muscle +
      l_skin * ce_skin + l_thymus * ce_thymus + l_gut * ce_gut + l_spleen * ce_spleen +
      l_pancreas * ce_pancreas + l_liver * ce_liver + l_other * ce_other - ltot * ce_lnode -
      clin_lnode * jin_all * ce_lnode * fuint_lnode + clin_lnode * ci_lnode * fucel_lnode
    d/dt(lnode_iw) <- clin_lnode * jin_all * ce_lnode * fuint_lnode - clin_lnode * ci_lnode * fucel_lnode
    d/dt(venous) <- ql_adipose * cv_adipose + ql_bone * cv_bone + ql_brain * cv_brain +
      ql_gonads * cv_gonads + ql_heart * cv_heart + ql_kidney * cv_kidney +
      ql_muscle * cv_muscle + ql_skin * cv_skin + ql_thymus * cv_thymus +
      ql_liver * cv_liver + q_lnode * cv_lnode + ql_other * cv_other +
      ltot * cv_lnode - q_lung * c_venous
    d/dt(arterial) <- ql_lung * cv_lung - (q_adipose + q_bone + q_brain + q_gonads + q_heart +
      q_kidney + q_muscle + q_skin + q_thymus + q_gut + q_spleen + q_pancreas +
      q_ha + q_by + q_lnode + q_other) * c_arterial
    d/dt(stomach) <- -stomach / t_stomach
    alag(stomach) <- lagtime # dose rows are placed at start + LagTime (PBPK_StudyDesign.m RegimenFun)
    d/dt(duodenum) <- stomach / t_stomach - duodenum / t_duodenum - clab_duodenum * duodenum / vflu_duodenum
    d/dt(jejunum) <- duodenum / t_duodenum - jejunum / t_jejunum - clab_jejunum * jejunum / vflu_jejunum
    d/dt(ileum) <- jejunum / t_jejunum - ileum / t_ileum - clab_ileum * ileum / vflu_ileum
    d/dt(colon) <- ileum / t_ileum - colon / t_colon - clab_colon * colon / vflu_colon
    d/dt(duodenum_uptake) <- clab_duodenum * duodenum / vflu_duodenum - kperup_eff * duodenum_uptake
    d/dt(jejunum_uptake) <- clab_jejunum * jejunum / vflu_jejunum - kperup_eff * jejunum_uptake
    d/dt(ileum_uptake) <- clab_ileum * ileum / vflu_ileum - kperup_eff * ileum_uptake
    d/dt(colon_uptake) <- clab_colon * colon / vflu_colon - kperup_eff * colon_uptake
    d/dt(duodenum_enterocyte) <- kperup_eff * duodenum_uptake - clin_gut * fgp * cent_duodenum * fucel_gut -
      clmet_duodenum * cent_duodenum * fucel_gut
    d/dt(jejunum_enterocyte) <- kperup_eff * jejunum_uptake - clin_gut * fgp * cent_jejunum * fucel_gut -
      clmet_jejunum * cent_jejunum * fucel_gut
    d/dt(ileum_enterocyte) <- kperup_eff * ileum_uptake - clin_gut * fgp * cent_ileum * fucel_gut -
      clmet_ileum * cent_ileum * fucel_gut
    d/dt(colon_enterocyte) <- kperup_eff * colon_uptake - clin_gut * fgp * cent_colon * fucel_gut
    d/dt(a_feces) <- colon / t_colon
    # Enzyme turnover, relative to baseline abundance (1 = baseline).
    d/dt(enzyme_cyp3a4_liver) <- 0.0077 * (1 + ind_cyp3a4) - (0.0077 + mbi_cyp3a4) * enzyme_cyp3a4_liver
    d/dt(enzyme_cyp2c19_liver) <- 0.0267 * (1 + ind_cyp2c19) - (0.0267) * enzyme_cyp2c19_liver
    d/dt(enzyme_cyp2d6_liver) <- 0.0143 * (1 + ind_cyp2d6) - (0.0143) * enzyme_cyp2d6_liver
    d/dt(enzyme_cyp2c8_liver) <- 0.0310 * (1 + ind_cyp2c8) - (0.0310) * enzyme_cyp2c8_liver
    d/dt(enzyme_cyp1a2_liver) <- 0.0149 * (1 + ind_cyp1a2) - (0.0149) * enzyme_cyp1a2_liver
    d/dt(enzyme_cyp2a6_liver) <- 0.0267 * (1 + ind_cyp2a6) - (0.0267) * enzyme_cyp2a6_liver
    d/dt(enzyme_cyp2b6_liver) <- 0.0217 * (1 + ind_cyp2b6) - (0.0217) * enzyme_cyp2b6_liver
    d/dt(enzyme_cyp2j2_liver) <- 0.01 * (1 + ind_cyp2j2) - (0.01 + mbi_cyp2j2) * enzyme_cyp2j2_liver
    d/dt(enzyme_ugt1a1_liver) <- 0.0693 * (1 + indu_ugt1a1) - 0.0693 * enzyme_ugt1a1_liver
    d/dt(enzyme_cyp3a4_duodenum) <- 0.03 * (1 + indg_du) - (0.03 + mbig_du) * enzyme_cyp3a4_duodenum
    d/dt(enzyme_cyp3a4_jejunum) <- 0.03 * (1 + indg_je) - (0.03 + mbig_je) * enzyme_cyp3a4_jejunum
    d/dt(enzyme_cyp3a4_ileum) <- 0.03 * (1 + indg_il) - (0.03 + mbig_il) * enzyme_cyp3a4_ileum
    enzyme_cyp3a4_liver(0) <- 1
    enzyme_cyp2c19_liver(0) <- 1
    enzyme_cyp2d6_liver(0) <- 1
    enzyme_cyp2c8_liver(0) <- 1
    enzyme_cyp1a2_liver(0) <- 1
    enzyme_cyp2a6_liver(0) <- 1
    enzyme_cyp2b6_liver(0) <- 1
    enzyme_cyp2j2_liver(0) <- 1
    enzyme_ugt1a1_liver(0) <- 1
    enzyme_cyp3a4_duodenum(0) <- 1
    enzyme_cyp3a4_jejunum(0) <- 1
    enzyme_cyp3a4_ileum(0) <- 1
    # --- perpetrator block ODEs (PBPK_ODE_solution.m rhs_function, amount form) ---
    d/dt(lung_vas_perpetrator) <- q_lung * c_venous_perpetrator - ql_lung * cv_lung_perpetrator - p_lung * fin_all_perpetrator * jin_adipose_perpetrator * cv_lung_perpetrator / bp_perpetrator + (p_lung - l_lung) * ce_lung_perpetrator
    d/dt(lung_ew_perpetrator) <- p_lung * fin_all_perpetrator * jin_adipose_perpetrator * cv_lung_perpetrator / bp_perpetrator - (p_lung - l_lung) * ce_lung_perpetrator - l_lung * ce_lung_perpetrator -
      clin_lung_perpetrator * jin_all_perpetrator * ce_lung_perpetrator * fuint_lung_perpetrator + clin_lung_perpetrator * ci_lung_perpetrator * fucel_lung_perpetrator
    d/dt(lung_iw_perpetrator) <- clin_lung_perpetrator * jin_all_perpetrator * ce_lung_perpetrator * fuint_lung_perpetrator - clin_lung_perpetrator * ci_lung_perpetrator * fucel_lung_perpetrator
    d/dt(adipose_vas_perpetrator) <- q_adipose * c_arterial_perpetrator - ql_adipose * cv_adipose_perpetrator - p_adipose * fin_all_perpetrator * cv_adipose_perpetrator / bp_perpetrator + (p_adipose - l_adipose) * ce_adipose_perpetrator
    d/dt(adipose_ew_perpetrator) <- p_adipose * fin_all_perpetrator * cv_adipose_perpetrator / bp_perpetrator - (p_adipose - l_adipose) * ce_adipose_perpetrator - l_adipose * ce_adipose_perpetrator -
      clin_adipose_perpetrator * jin_all_perpetrator * jin_adipose_perpetrator * ce_adipose_perpetrator * fuint_adipose_perpetrator + clin_adipose_perpetrator * ci_adipose_perpetrator * fucel_adipose_perpetrator
    d/dt(adipose_iw_perpetrator) <- clin_adipose_perpetrator * jin_all_perpetrator * jin_adipose_perpetrator * ce_adipose_perpetrator * fuint_adipose_perpetrator - clin_adipose_perpetrator * ci_adipose_perpetrator * fucel_adipose_perpetrator
    d/dt(bone_vas_perpetrator) <- q_bone * c_arterial_perpetrator - ql_bone * cv_bone_perpetrator - p_bone * fin_all_perpetrator * cv_bone_perpetrator / bp_perpetrator + (p_bone - l_bone) * ce_bone_perpetrator
    d/dt(bone_ew_perpetrator) <- p_bone * fin_all_perpetrator * cv_bone_perpetrator / bp_perpetrator - (p_bone - l_bone) * ce_bone_perpetrator - l_bone * ce_bone_perpetrator -
      clin_bone_perpetrator * jin_all_perpetrator * ce_bone_perpetrator * fuint_bone_perpetrator + clin_bone_perpetrator * ci_bone_perpetrator * fucel_bone_perpetrator
    d/dt(bone_iw_perpetrator) <- clin_bone_perpetrator * jin_all_perpetrator * ce_bone_perpetrator * fuint_bone_perpetrator - clin_bone_perpetrator * ci_bone_perpetrator * fucel_bone_perpetrator
    d/dt(brain_vas_perpetrator) <- q_brain * c_arterial_perpetrator - ql_brain * cv_brain_perpetrator - p_brain * fin_all_perpetrator * cv_brain_perpetrator / bp_perpetrator + (p_brain - l_brain) * ce_brain_perpetrator
    d/dt(brain_ew_perpetrator) <- p_brain * fin_all_perpetrator * cv_brain_perpetrator / bp_perpetrator - (p_brain - l_brain) * ce_brain_perpetrator - l_brain * ce_brain_perpetrator -
      clin_brain_perpetrator * jin_all_perpetrator * ce_brain_perpetrator * fuint_brain_perpetrator + clin_brain_perpetrator * ci_brain_perpetrator * fucel_brain_perpetrator
    d/dt(brain_iw_perpetrator) <- clin_brain_perpetrator * jin_all_perpetrator * ce_brain_perpetrator * fuint_brain_perpetrator - clin_brain_perpetrator * ci_brain_perpetrator * fucel_brain_perpetrator
    d/dt(gonads_vas_perpetrator) <- q_gonads * c_arterial_perpetrator - ql_gonads * cv_gonads_perpetrator - p_gonads * fin_all_perpetrator * cv_gonads_perpetrator / bp_perpetrator + (p_gonads - l_gonads) * ce_gonads_perpetrator
    d/dt(gonads_ew_perpetrator) <- p_gonads * fin_all_perpetrator * cv_gonads_perpetrator / bp_perpetrator - (p_gonads - l_gonads) * ce_gonads_perpetrator - l_gonads * ce_gonads_perpetrator -
      clin_gonads_perpetrator * jin_all_perpetrator * ce_gonads_perpetrator * fuint_gonads_perpetrator + clin_gonads_perpetrator * ci_gonads_perpetrator * fucel_gonads_perpetrator
    d/dt(gonads_iw_perpetrator) <- clin_gonads_perpetrator * jin_all_perpetrator * ce_gonads_perpetrator * fuint_gonads_perpetrator - clin_gonads_perpetrator * ci_gonads_perpetrator * fucel_gonads_perpetrator
    d/dt(heart_vas_perpetrator) <- q_heart * c_arterial_perpetrator - ql_heart * cv_heart_perpetrator - p_heart * fin_all_perpetrator * cv_heart_perpetrator / bp_perpetrator + (p_heart - l_heart) * ce_heart_perpetrator
    d/dt(heart_ew_perpetrator) <- p_heart * fin_all_perpetrator * cv_heart_perpetrator / bp_perpetrator - (p_heart - l_heart) * ce_heart_perpetrator - l_heart * ce_heart_perpetrator -
      clin_heart_perpetrator * jin_all_perpetrator * ce_heart_perpetrator * fuint_heart_perpetrator + clin_heart_perpetrator * ci_heart_perpetrator * fucel_heart_perpetrator
    d/dt(heart_iw_perpetrator) <- clin_heart_perpetrator * jin_all_perpetrator * ce_heart_perpetrator * fuint_heart_perpetrator - clin_heart_perpetrator * ci_heart_perpetrator * fucel_heart_perpetrator
    d/dt(kidney_vas_perpetrator) <- q_kidney * c_arterial_perpetrator - ql_kidney * cv_kidney_perpetrator - p_kidney * fin_all_perpetrator * jin_muscle_perpetrator * cv_kidney_perpetrator / bp_perpetrator + (p_kidney - l_kidney) * ce_kidney_perpetrator - clr_perpetrator * cv_kidney_perpetrator
    d/dt(kidney_ew_perpetrator) <- p_kidney * fin_all_perpetrator * jin_muscle_perpetrator * cv_kidney_perpetrator / bp_perpetrator - (p_kidney - l_kidney) * ce_kidney_perpetrator - l_kidney * ce_kidney_perpetrator -
      clin_kidney_perpetrator * jin_all_perpetrator * ce_kidney_perpetrator * fuint_kidney_perpetrator + clin_kidney_perpetrator * ci_kidney_perpetrator * fucel_kidney_perpetrator
    d/dt(kidney_iw_perpetrator) <- clin_kidney_perpetrator * jin_all_perpetrator * ce_kidney_perpetrator * fuint_kidney_perpetrator - clin_kidney_perpetrator * ci_kidney_perpetrator * fucel_kidney_perpetrator
    d/dt(muscle_vas_perpetrator) <- q_muscle * c_arterial_perpetrator - ql_muscle * cv_muscle_perpetrator - p_muscle * fin_all_perpetrator * cv_muscle_perpetrator / bp_perpetrator + (p_muscle - l_muscle) * ce_muscle_perpetrator
    d/dt(muscle_ew_perpetrator) <- p_muscle * fin_all_perpetrator * cv_muscle_perpetrator / bp_perpetrator - (p_muscle - l_muscle) * ce_muscle_perpetrator - l_muscle * ce_muscle_perpetrator -
      clin_muscle_perpetrator * jin_all_perpetrator * jin_muscle_perpetrator * ce_muscle_perpetrator * fuint_muscle_perpetrator + clin_muscle_perpetrator * ci_muscle_perpetrator * fucel_muscle_perpetrator
    d/dt(muscle_iw_perpetrator) <- clin_muscle_perpetrator * jin_all_perpetrator * jin_muscle_perpetrator * ce_muscle_perpetrator * fuint_muscle_perpetrator - clin_muscle_perpetrator * ci_muscle_perpetrator * fucel_muscle_perpetrator
    d/dt(skin_vas_perpetrator) <- q_skin * c_arterial_perpetrator - ql_skin * cv_skin_perpetrator - p_skin * fin_all_perpetrator * cv_skin_perpetrator / bp_perpetrator + (p_skin - l_skin) * ce_skin_perpetrator
    d/dt(skin_ew_perpetrator) <- p_skin * fin_all_perpetrator * cv_skin_perpetrator / bp_perpetrator - (p_skin - l_skin) * ce_skin_perpetrator - l_skin * ce_skin_perpetrator -
      clin_skin_perpetrator * jin_all_perpetrator * ce_skin_perpetrator * fuint_skin_perpetrator + clin_skin_perpetrator * ci_skin_perpetrator * fucel_skin_perpetrator
    d/dt(skin_iw_perpetrator) <- clin_skin_perpetrator * jin_all_perpetrator * ce_skin_perpetrator * fuint_skin_perpetrator - clin_skin_perpetrator * ci_skin_perpetrator * fucel_skin_perpetrator
    d/dt(thymus_vas_perpetrator) <- q_thymus * c_arterial_perpetrator - ql_thymus * cv_thymus_perpetrator - p_thymus * fin_all_perpetrator * cv_thymus_perpetrator / bp_perpetrator + (p_thymus - l_thymus) * ce_thymus_perpetrator
    d/dt(thymus_ew_perpetrator) <- p_thymus * fin_all_perpetrator * cv_thymus_perpetrator / bp_perpetrator - (p_thymus - l_thymus) * ce_thymus_perpetrator - l_thymus * ce_thymus_perpetrator -
      clin_thymus_perpetrator * jin_all_perpetrator * ce_thymus_perpetrator * fuint_thymus_perpetrator + clin_thymus_perpetrator * ci_thymus_perpetrator * fucel_thymus_perpetrator
    d/dt(thymus_iw_perpetrator) <- clin_thymus_perpetrator * jin_all_perpetrator * ce_thymus_perpetrator * fuint_thymus_perpetrator - clin_thymus_perpetrator * ci_thymus_perpetrator * fucel_thymus_perpetrator
    d/dt(spleen_vas_perpetrator) <- q_spleen * c_arterial_perpetrator - ql_spleen * cv_spleen_perpetrator - p_spleen * fin_all_perpetrator * cv_spleen_perpetrator / bp_perpetrator + (p_spleen - l_spleen) * ce_spleen_perpetrator
    d/dt(spleen_ew_perpetrator) <- p_spleen * fin_all_perpetrator * cv_spleen_perpetrator / bp_perpetrator - (p_spleen - l_spleen) * ce_spleen_perpetrator - l_spleen * ce_spleen_perpetrator -
      clin_spleen_perpetrator * jin_all_perpetrator * ce_spleen_perpetrator * fuint_spleen_perpetrator + clin_spleen_perpetrator * ci_spleen_perpetrator * fucel_spleen_perpetrator
    d/dt(spleen_iw_perpetrator) <- clin_spleen_perpetrator * jin_all_perpetrator * ce_spleen_perpetrator * fuint_spleen_perpetrator - clin_spleen_perpetrator * ci_spleen_perpetrator * fucel_spleen_perpetrator
    d/dt(pancreas_vas_perpetrator) <- q_pancreas * c_arterial_perpetrator - ql_pancreas * cv_pancreas_perpetrator - p_pancreas * fin_all_perpetrator * jin_liver_perpetrator * cv_pancreas_perpetrator / bp_perpetrator + (p_pancreas - l_pancreas) * ce_pancreas_perpetrator
    d/dt(pancreas_ew_perpetrator) <- p_pancreas * fin_all_perpetrator * jin_liver_perpetrator * cv_pancreas_perpetrator / bp_perpetrator - (p_pancreas - l_pancreas) * ce_pancreas_perpetrator - l_pancreas * ce_pancreas_perpetrator -
      clin_pancreas_perpetrator * jin_all_perpetrator * ce_pancreas_perpetrator * fuint_pancreas_perpetrator + clin_pancreas_perpetrator * ci_pancreas_perpetrator * fucel_pancreas_perpetrator
    d/dt(pancreas_iw_perpetrator) <- clin_pancreas_perpetrator * jin_all_perpetrator * ce_pancreas_perpetrator * fuint_pancreas_perpetrator - clin_pancreas_perpetrator * ci_pancreas_perpetrator * fucel_pancreas_perpetrator
    d/dt(other_vas_perpetrator) <- q_other * c_arterial_perpetrator - ql_other * cv_other_perpetrator - p_other * fin_all_perpetrator * cv_other_perpetrator / bp_perpetrator + (p_other - l_other) * ce_other_perpetrator
    d/dt(other_ew_perpetrator) <- p_other * fin_all_perpetrator * cv_other_perpetrator / bp_perpetrator - (p_other - l_other) * ce_other_perpetrator - l_other * ce_other_perpetrator -
      clin_other_perpetrator * jin_all_perpetrator * ce_other_perpetrator * fuint_other_perpetrator + clin_other_perpetrator * ci_other_perpetrator * fucel_other_perpetrator
    d/dt(other_iw_perpetrator) <- clin_other_perpetrator * jin_all_perpetrator * ce_other_perpetrator * fuint_other_perpetrator - clin_other_perpetrator * ci_other_perpetrator * fucel_other_perpetrator
    d/dt(gut_vas_perpetrator) <- q_gut * c_arterial_perpetrator - ql_gut * cv_gut_perpetrator - p_gut * fin_all_perpetrator * cv_gut_perpetrator / bp_perpetrator + (p_gut - l_gut) * ce_gut_perpetrator
    d/dt(gut_ew_perpetrator) <- p_gut * fin_all_perpetrator * cv_gut_perpetrator / bp_perpetrator - (p_gut - l_gut) * ce_gut_perpetrator - l_gut * ce_gut_perpetrator +
      clin_gut_perpetrator * fgp_perpetrator * fucel_gut_perpetrator * (cent_duodenum_perpetrator + cent_jejunum_perpetrator + cent_ileum_perpetrator + cent_colon_perpetrator)
    d/dt(liver_vas_perpetrator) <- ql_gut * cv_gut_perpetrator + ql_spleen * cv_spleen_perpetrator + ql_pancreas * cv_pancreas_perpetrator + (q_ha + q_by) * c_arterial_perpetrator -
      ql_liver * cv_liver_perpetrator - p_liver * fin_all_perpetrator * cv_liver_perpetrator / bp_perpetrator + (p_liver - l_liver) * ce_liver_perpetrator
    d/dt(liver_ew_perpetrator) <- p_liver * fin_all_perpetrator * cv_liver_perpetrator / bp_perpetrator - (p_liver - l_liver) * ce_liver_perpetrator - l_liver * ce_liver_perpetrator -
      clin_liver_perpetrator * jin_all_perpetrator * jin_liver_perpetrator * ce_liver_perpetrator * fuint_liver_perpetrator + clin_liver_perpetrator * ci_liver_perpetrator * fucel_liver_perpetrator
    d/dt(liver_iw_perpetrator) <- clin_liver_perpetrator * jin_all_perpetrator * jin_liver_perpetrator * ce_liver_perpetrator * fuint_liver_perpetrator - clin_liver_perpetrator * ci_liver_perpetrator * fucel_liver_perpetrator -
      clmet_liver_perpetrator * ci_liver_perpetrator * fucel_liver_perpetrator
    d/dt(lnode_vas_perpetrator) <- q_lnode * c_arterial_perpetrator - q_lnode * cv_lnode_perpetrator + ltot * ce_lnode_perpetrator - ltot * cv_lnode_perpetrator
    d/dt(lnode_ew_perpetrator) <- l_lung * ce_lung_perpetrator + l_adipose * ce_adipose_perpetrator + l_bone * ce_bone_perpetrator + l_brain * ce_brain_perpetrator +
      l_gonads * ce_gonads_perpetrator + l_heart * ce_heart_perpetrator + l_kidney * ce_kidney_perpetrator + l_muscle * ce_muscle_perpetrator +
      l_skin * ce_skin_perpetrator + l_thymus * ce_thymus_perpetrator + l_gut * ce_gut_perpetrator + l_spleen * ce_spleen_perpetrator +
      l_pancreas * ce_pancreas_perpetrator + l_liver * ce_liver_perpetrator + l_other * ce_other_perpetrator - ltot * ce_lnode_perpetrator -
      clin_lnode_perpetrator * jin_all_perpetrator * ce_lnode_perpetrator * fuint_lnode_perpetrator + clin_lnode_perpetrator * ci_lnode_perpetrator * fucel_lnode_perpetrator
    d/dt(lnode_iw_perpetrator) <- clin_lnode_perpetrator * jin_all_perpetrator * ce_lnode_perpetrator * fuint_lnode_perpetrator - clin_lnode_perpetrator * ci_lnode_perpetrator * fucel_lnode_perpetrator
    d/dt(venous_perpetrator) <- ql_adipose * cv_adipose_perpetrator + ql_bone * cv_bone_perpetrator + ql_brain * cv_brain_perpetrator +
      ql_gonads * cv_gonads_perpetrator + ql_heart * cv_heart_perpetrator + ql_kidney * cv_kidney_perpetrator +
      ql_muscle * cv_muscle_perpetrator + ql_skin * cv_skin_perpetrator + ql_thymus * cv_thymus_perpetrator +
      ql_liver * cv_liver_perpetrator + q_lnode * cv_lnode_perpetrator + ql_other * cv_other_perpetrator +
      ltot * cv_lnode_perpetrator - q_lung * c_venous_perpetrator
    d/dt(arterial_perpetrator) <- ql_lung * cv_lung_perpetrator - (q_adipose + q_bone + q_brain + q_gonads + q_heart +
      q_kidney + q_muscle + q_skin + q_thymus + q_gut + q_spleen + q_pancreas +
      q_ha + q_by + q_lnode + q_other) * c_arterial_perpetrator
    d/dt(stomach_perpetrator) <- -stomach_perpetrator / t_stomach
    alag(stomach_perpetrator) <- lagtime_perpetrator # dose rows are placed at start + LagTime (PBPK_StudyDesign.m RegimenFun)
    d/dt(duodenum_perpetrator) <- stomach_perpetrator / t_stomach - duodenum_perpetrator / t_duodenum - clab_duodenum_perpetrator * duodenum_perpetrator / vflu_duodenum
    d/dt(jejunum_perpetrator) <- duodenum_perpetrator / t_duodenum - jejunum_perpetrator / t_jejunum - clab_jejunum_perpetrator * jejunum_perpetrator / vflu_jejunum
    d/dt(ileum_perpetrator) <- jejunum_perpetrator / t_jejunum - ileum_perpetrator / t_ileum - clab_ileum_perpetrator * ileum_perpetrator / vflu_ileum
    d/dt(colon_perpetrator) <- ileum_perpetrator / t_ileum - colon_perpetrator / t_colon - clab_colon_perpetrator * colon_perpetrator / vflu_colon
    d/dt(duodenum_uptake_perpetrator) <- clab_duodenum_perpetrator * duodenum_perpetrator / vflu_duodenum - kperup_eff_perpetrator * duodenum_uptake_perpetrator
    d/dt(jejunum_uptake_perpetrator) <- clab_jejunum_perpetrator * jejunum_perpetrator / vflu_jejunum - kperup_eff_perpetrator * jejunum_uptake_perpetrator
    d/dt(ileum_uptake_perpetrator) <- clab_ileum_perpetrator * ileum_perpetrator / vflu_ileum - kperup_eff_perpetrator * ileum_uptake_perpetrator
    d/dt(colon_uptake_perpetrator) <- clab_colon_perpetrator * colon_perpetrator / vflu_colon - kperup_eff_perpetrator * colon_uptake_perpetrator
    d/dt(duodenum_enterocyte_perpetrator) <- kperup_eff_perpetrator * duodenum_uptake_perpetrator - clin_gut_perpetrator * fgp_perpetrator * cent_duodenum_perpetrator * fucel_gut_perpetrator -
      clmet_duodenum_perpetrator * cent_duodenum_perpetrator * fucel_gut_perpetrator
    d/dt(jejunum_enterocyte_perpetrator) <- kperup_eff_perpetrator * jejunum_uptake_perpetrator - clin_gut_perpetrator * fgp_perpetrator * cent_jejunum_perpetrator * fucel_gut_perpetrator -
      clmet_jejunum_perpetrator * cent_jejunum_perpetrator * fucel_gut_perpetrator
    d/dt(ileum_enterocyte_perpetrator) <- kperup_eff_perpetrator * ileum_uptake_perpetrator - clin_gut_perpetrator * fgp_perpetrator * cent_ileum_perpetrator * fucel_gut_perpetrator -
      clmet_ileum_perpetrator * cent_ileum_perpetrator * fucel_gut_perpetrator
    d/dt(colon_enterocyte_perpetrator) <- kperup_eff_perpetrator * colon_uptake_perpetrator - clin_gut_perpetrator * fgp_perpetrator * cent_colon_perpetrator * fucel_gut_perpetrator
    d/dt(a_feces_perpetrator) <- colon_perpetrator / t_colon
    # Enzyme turnover, relative to baseline abundance (1 = baseline).
    d/dt(enzyme_cyp3a4_liver_perpetrator) <- 0.0077 * (1 + ind_cyp3a4_perpetrator) - (0.0077 + mbi_cyp3a4_perpetrator) * enzyme_cyp3a4_liver_perpetrator
    d/dt(enzyme_cyp2c19_liver_perpetrator) <- 0.0267 * (1 + ind_cyp2c19_perpetrator) - (0.0267) * enzyme_cyp2c19_liver_perpetrator
    d/dt(enzyme_cyp2d6_liver_perpetrator) <- 0.0143 * (1 + ind_cyp2d6_perpetrator) - (0.0143) * enzyme_cyp2d6_liver_perpetrator
    d/dt(enzyme_cyp2c8_liver_perpetrator) <- 0.0310 * (1 + ind_cyp2c8_perpetrator) - (0.0310) * enzyme_cyp2c8_liver_perpetrator
    d/dt(enzyme_cyp1a2_liver_perpetrator) <- 0.0149 * (1 + ind_cyp1a2_perpetrator) - (0.0149) * enzyme_cyp1a2_liver_perpetrator
    d/dt(enzyme_cyp2a6_liver_perpetrator) <- 0.0267 * (1 + ind_cyp2a6_perpetrator) - (0.0267) * enzyme_cyp2a6_liver_perpetrator
    d/dt(enzyme_cyp2b6_liver_perpetrator) <- 0.0217 * (1 + ind_cyp2b6_perpetrator) - (0.0217) * enzyme_cyp2b6_liver_perpetrator
    d/dt(enzyme_cyp2j2_liver_perpetrator) <- 0.01 * (1 + ind_cyp2j2_perpetrator) - (0.01 + mbi_cyp2j2_perpetrator) * enzyme_cyp2j2_liver_perpetrator
    d/dt(enzyme_ugt1a1_liver_perpetrator) <- 0.0693 * (1 + indu_ugt1a1_perpetrator) - 0.0693 * enzyme_ugt1a1_liver_perpetrator
    d/dt(enzyme_cyp3a4_duodenum_perpetrator) <- 0.03 * (1 + indg_du_perpetrator) - (0.03 + mbig_du_perpetrator) * enzyme_cyp3a4_duodenum_perpetrator
    d/dt(enzyme_cyp3a4_jejunum_perpetrator) <- 0.03 * (1 + indg_je_perpetrator) - (0.03 + mbig_je_perpetrator) * enzyme_cyp3a4_jejunum_perpetrator
    d/dt(enzyme_cyp3a4_ileum_perpetrator) <- 0.03 * (1 + indg_il_perpetrator) - (0.03 + mbig_il_perpetrator) * enzyme_cyp3a4_ileum_perpetrator
    enzyme_cyp3a4_liver_perpetrator(0) <- 1
    enzyme_cyp2c19_liver_perpetrator(0) <- 1
    enzyme_cyp2d6_liver_perpetrator(0) <- 1
    enzyme_cyp2c8_liver_perpetrator(0) <- 1
    enzyme_cyp1a2_liver_perpetrator(0) <- 1
    enzyme_cyp2a6_liver_perpetrator(0) <- 1
    enzyme_cyp2b6_liver_perpetrator(0) <- 1
    enzyme_cyp2j2_liver_perpetrator(0) <- 1
    enzyme_ugt1a1_liver_perpetrator(0) <- 1
    enzyme_cyp3a4_duodenum_perpetrator(0) <- 1
    enzyme_cyp3a4_jejunum_perpetrator(0) <- 1
    enzyme_cyp3a4_ileum_perpetrator(0) <- 1
    # --- perpetrator2 block ODEs (PBPK_ODE_solution.m rhs_function, amount form) ---
    d/dt(lung_vas_perpetrator2) <- q_lung * c_venous_perpetrator2 - ql_lung * cv_lung_perpetrator2 - p_lung * fin_all_perpetrator2 * jin_adipose_perpetrator2 * cv_lung_perpetrator2 / bp_perpetrator2 + (p_lung - l_lung) * ce_lung_perpetrator2
    d/dt(lung_ew_perpetrator2) <- p_lung * fin_all_perpetrator2 * jin_adipose_perpetrator2 * cv_lung_perpetrator2 / bp_perpetrator2 - (p_lung - l_lung) * ce_lung_perpetrator2 - l_lung * ce_lung_perpetrator2 -
      clin_lung_perpetrator2 * jin_all_perpetrator2 * ce_lung_perpetrator2 * fuint_lung_perpetrator2 + clin_lung_perpetrator2 * ci_lung_perpetrator2 * fucel_lung_perpetrator2
    d/dt(lung_iw_perpetrator2) <- clin_lung_perpetrator2 * jin_all_perpetrator2 * ce_lung_perpetrator2 * fuint_lung_perpetrator2 - clin_lung_perpetrator2 * ci_lung_perpetrator2 * fucel_lung_perpetrator2
    d/dt(adipose_vas_perpetrator2) <- q_adipose * c_arterial_perpetrator2 - ql_adipose * cv_adipose_perpetrator2 - p_adipose * fin_all_perpetrator2 * cv_adipose_perpetrator2 / bp_perpetrator2 + (p_adipose - l_adipose) * ce_adipose_perpetrator2
    d/dt(adipose_ew_perpetrator2) <- p_adipose * fin_all_perpetrator2 * cv_adipose_perpetrator2 / bp_perpetrator2 - (p_adipose - l_adipose) * ce_adipose_perpetrator2 - l_adipose * ce_adipose_perpetrator2 -
      clin_adipose_perpetrator2 * jin_all_perpetrator2 * jin_adipose_perpetrator2 * ce_adipose_perpetrator2 * fuint_adipose_perpetrator2 + clin_adipose_perpetrator2 * ci_adipose_perpetrator2 * fucel_adipose_perpetrator2
    d/dt(adipose_iw_perpetrator2) <- clin_adipose_perpetrator2 * jin_all_perpetrator2 * jin_adipose_perpetrator2 * ce_adipose_perpetrator2 * fuint_adipose_perpetrator2 - clin_adipose_perpetrator2 * ci_adipose_perpetrator2 * fucel_adipose_perpetrator2
    d/dt(bone_vas_perpetrator2) <- q_bone * c_arterial_perpetrator2 - ql_bone * cv_bone_perpetrator2 - p_bone * fin_all_perpetrator2 * cv_bone_perpetrator2 / bp_perpetrator2 + (p_bone - l_bone) * ce_bone_perpetrator2
    d/dt(bone_ew_perpetrator2) <- p_bone * fin_all_perpetrator2 * cv_bone_perpetrator2 / bp_perpetrator2 - (p_bone - l_bone) * ce_bone_perpetrator2 - l_bone * ce_bone_perpetrator2 -
      clin_bone_perpetrator2 * jin_all_perpetrator2 * ce_bone_perpetrator2 * fuint_bone_perpetrator2 + clin_bone_perpetrator2 * ci_bone_perpetrator2 * fucel_bone_perpetrator2
    d/dt(bone_iw_perpetrator2) <- clin_bone_perpetrator2 * jin_all_perpetrator2 * ce_bone_perpetrator2 * fuint_bone_perpetrator2 - clin_bone_perpetrator2 * ci_bone_perpetrator2 * fucel_bone_perpetrator2
    d/dt(brain_vas_perpetrator2) <- q_brain * c_arterial_perpetrator2 - ql_brain * cv_brain_perpetrator2 - p_brain * fin_all_perpetrator2 * cv_brain_perpetrator2 / bp_perpetrator2 + (p_brain - l_brain) * ce_brain_perpetrator2
    d/dt(brain_ew_perpetrator2) <- p_brain * fin_all_perpetrator2 * cv_brain_perpetrator2 / bp_perpetrator2 - (p_brain - l_brain) * ce_brain_perpetrator2 - l_brain * ce_brain_perpetrator2 -
      clin_brain_perpetrator2 * jin_all_perpetrator2 * ce_brain_perpetrator2 * fuint_brain_perpetrator2 + clin_brain_perpetrator2 * ci_brain_perpetrator2 * fucel_brain_perpetrator2
    d/dt(brain_iw_perpetrator2) <- clin_brain_perpetrator2 * jin_all_perpetrator2 * ce_brain_perpetrator2 * fuint_brain_perpetrator2 - clin_brain_perpetrator2 * ci_brain_perpetrator2 * fucel_brain_perpetrator2
    d/dt(gonads_vas_perpetrator2) <- q_gonads * c_arterial_perpetrator2 - ql_gonads * cv_gonads_perpetrator2 - p_gonads * fin_all_perpetrator2 * cv_gonads_perpetrator2 / bp_perpetrator2 + (p_gonads - l_gonads) * ce_gonads_perpetrator2
    d/dt(gonads_ew_perpetrator2) <- p_gonads * fin_all_perpetrator2 * cv_gonads_perpetrator2 / bp_perpetrator2 - (p_gonads - l_gonads) * ce_gonads_perpetrator2 - l_gonads * ce_gonads_perpetrator2 -
      clin_gonads_perpetrator2 * jin_all_perpetrator2 * ce_gonads_perpetrator2 * fuint_gonads_perpetrator2 + clin_gonads_perpetrator2 * ci_gonads_perpetrator2 * fucel_gonads_perpetrator2
    d/dt(gonads_iw_perpetrator2) <- clin_gonads_perpetrator2 * jin_all_perpetrator2 * ce_gonads_perpetrator2 * fuint_gonads_perpetrator2 - clin_gonads_perpetrator2 * ci_gonads_perpetrator2 * fucel_gonads_perpetrator2
    d/dt(heart_vas_perpetrator2) <- q_heart * c_arterial_perpetrator2 - ql_heart * cv_heart_perpetrator2 - p_heart * fin_all_perpetrator2 * cv_heart_perpetrator2 / bp_perpetrator2 + (p_heart - l_heart) * ce_heart_perpetrator2
    d/dt(heart_ew_perpetrator2) <- p_heart * fin_all_perpetrator2 * cv_heart_perpetrator2 / bp_perpetrator2 - (p_heart - l_heart) * ce_heart_perpetrator2 - l_heart * ce_heart_perpetrator2 -
      clin_heart_perpetrator2 * jin_all_perpetrator2 * ce_heart_perpetrator2 * fuint_heart_perpetrator2 + clin_heart_perpetrator2 * ci_heart_perpetrator2 * fucel_heart_perpetrator2
    d/dt(heart_iw_perpetrator2) <- clin_heart_perpetrator2 * jin_all_perpetrator2 * ce_heart_perpetrator2 * fuint_heart_perpetrator2 - clin_heart_perpetrator2 * ci_heart_perpetrator2 * fucel_heart_perpetrator2
    d/dt(kidney_vas_perpetrator2) <- q_kidney * c_arterial_perpetrator2 - ql_kidney * cv_kidney_perpetrator2 - p_kidney * fin_all_perpetrator2 * jin_muscle_perpetrator2 * cv_kidney_perpetrator2 / bp_perpetrator2 + (p_kidney - l_kidney) * ce_kidney_perpetrator2 - clr_perpetrator2 * cv_kidney_perpetrator2
    d/dt(kidney_ew_perpetrator2) <- p_kidney * fin_all_perpetrator2 * jin_muscle_perpetrator2 * cv_kidney_perpetrator2 / bp_perpetrator2 - (p_kidney - l_kidney) * ce_kidney_perpetrator2 - l_kidney * ce_kidney_perpetrator2 -
      clin_kidney_perpetrator2 * jin_all_perpetrator2 * ce_kidney_perpetrator2 * fuint_kidney_perpetrator2 + clin_kidney_perpetrator2 * ci_kidney_perpetrator2 * fucel_kidney_perpetrator2
    d/dt(kidney_iw_perpetrator2) <- clin_kidney_perpetrator2 * jin_all_perpetrator2 * ce_kidney_perpetrator2 * fuint_kidney_perpetrator2 - clin_kidney_perpetrator2 * ci_kidney_perpetrator2 * fucel_kidney_perpetrator2
    d/dt(muscle_vas_perpetrator2) <- q_muscle * c_arterial_perpetrator2 - ql_muscle * cv_muscle_perpetrator2 - p_muscle * fin_all_perpetrator2 * cv_muscle_perpetrator2 / bp_perpetrator2 + (p_muscle - l_muscle) * ce_muscle_perpetrator2
    d/dt(muscle_ew_perpetrator2) <- p_muscle * fin_all_perpetrator2 * cv_muscle_perpetrator2 / bp_perpetrator2 - (p_muscle - l_muscle) * ce_muscle_perpetrator2 - l_muscle * ce_muscle_perpetrator2 -
      clin_muscle_perpetrator2 * jin_all_perpetrator2 * jin_muscle_perpetrator2 * ce_muscle_perpetrator2 * fuint_muscle_perpetrator2 + clin_muscle_perpetrator2 * ci_muscle_perpetrator2 * fucel_muscle_perpetrator2
    d/dt(muscle_iw_perpetrator2) <- clin_muscle_perpetrator2 * jin_all_perpetrator2 * jin_muscle_perpetrator2 * ce_muscle_perpetrator2 * fuint_muscle_perpetrator2 - clin_muscle_perpetrator2 * ci_muscle_perpetrator2 * fucel_muscle_perpetrator2
    d/dt(skin_vas_perpetrator2) <- q_skin * c_arterial_perpetrator2 - ql_skin * cv_skin_perpetrator2 - p_skin * fin_all_perpetrator2 * cv_skin_perpetrator2 / bp_perpetrator2 + (p_skin - l_skin) * ce_skin_perpetrator2
    d/dt(skin_ew_perpetrator2) <- p_skin * fin_all_perpetrator2 * cv_skin_perpetrator2 / bp_perpetrator2 - (p_skin - l_skin) * ce_skin_perpetrator2 - l_skin * ce_skin_perpetrator2 -
      clin_skin_perpetrator2 * jin_all_perpetrator2 * ce_skin_perpetrator2 * fuint_skin_perpetrator2 + clin_skin_perpetrator2 * ci_skin_perpetrator2 * fucel_skin_perpetrator2
    d/dt(skin_iw_perpetrator2) <- clin_skin_perpetrator2 * jin_all_perpetrator2 * ce_skin_perpetrator2 * fuint_skin_perpetrator2 - clin_skin_perpetrator2 * ci_skin_perpetrator2 * fucel_skin_perpetrator2
    d/dt(thymus_vas_perpetrator2) <- q_thymus * c_arterial_perpetrator2 - ql_thymus * cv_thymus_perpetrator2 - p_thymus * fin_all_perpetrator2 * cv_thymus_perpetrator2 / bp_perpetrator2 + (p_thymus - l_thymus) * ce_thymus_perpetrator2
    d/dt(thymus_ew_perpetrator2) <- p_thymus * fin_all_perpetrator2 * cv_thymus_perpetrator2 / bp_perpetrator2 - (p_thymus - l_thymus) * ce_thymus_perpetrator2 - l_thymus * ce_thymus_perpetrator2 -
      clin_thymus_perpetrator2 * jin_all_perpetrator2 * ce_thymus_perpetrator2 * fuint_thymus_perpetrator2 + clin_thymus_perpetrator2 * ci_thymus_perpetrator2 * fucel_thymus_perpetrator2
    d/dt(thymus_iw_perpetrator2) <- clin_thymus_perpetrator2 * jin_all_perpetrator2 * ce_thymus_perpetrator2 * fuint_thymus_perpetrator2 - clin_thymus_perpetrator2 * ci_thymus_perpetrator2 * fucel_thymus_perpetrator2
    d/dt(spleen_vas_perpetrator2) <- q_spleen * c_arterial_perpetrator2 - ql_spleen * cv_spleen_perpetrator2 - p_spleen * fin_all_perpetrator2 * cv_spleen_perpetrator2 / bp_perpetrator2 + (p_spleen - l_spleen) * ce_spleen_perpetrator2
    d/dt(spleen_ew_perpetrator2) <- p_spleen * fin_all_perpetrator2 * cv_spleen_perpetrator2 / bp_perpetrator2 - (p_spleen - l_spleen) * ce_spleen_perpetrator2 - l_spleen * ce_spleen_perpetrator2 -
      clin_spleen_perpetrator2 * jin_all_perpetrator2 * ce_spleen_perpetrator2 * fuint_spleen_perpetrator2 + clin_spleen_perpetrator2 * ci_spleen_perpetrator2 * fucel_spleen_perpetrator2
    d/dt(spleen_iw_perpetrator2) <- clin_spleen_perpetrator2 * jin_all_perpetrator2 * ce_spleen_perpetrator2 * fuint_spleen_perpetrator2 - clin_spleen_perpetrator2 * ci_spleen_perpetrator2 * fucel_spleen_perpetrator2
    d/dt(pancreas_vas_perpetrator2) <- q_pancreas * c_arterial_perpetrator2 - ql_pancreas * cv_pancreas_perpetrator2 - p_pancreas * fin_all_perpetrator2 * jin_liver_perpetrator2 * cv_pancreas_perpetrator2 / bp_perpetrator2 + (p_pancreas - l_pancreas) * ce_pancreas_perpetrator2
    d/dt(pancreas_ew_perpetrator2) <- p_pancreas * fin_all_perpetrator2 * jin_liver_perpetrator2 * cv_pancreas_perpetrator2 / bp_perpetrator2 - (p_pancreas - l_pancreas) * ce_pancreas_perpetrator2 - l_pancreas * ce_pancreas_perpetrator2 -
      clin_pancreas_perpetrator2 * jin_all_perpetrator2 * ce_pancreas_perpetrator2 * fuint_pancreas_perpetrator2 + clin_pancreas_perpetrator2 * ci_pancreas_perpetrator2 * fucel_pancreas_perpetrator2
    d/dt(pancreas_iw_perpetrator2) <- clin_pancreas_perpetrator2 * jin_all_perpetrator2 * ce_pancreas_perpetrator2 * fuint_pancreas_perpetrator2 - clin_pancreas_perpetrator2 * ci_pancreas_perpetrator2 * fucel_pancreas_perpetrator2
    d/dt(other_vas_perpetrator2) <- q_other * c_arterial_perpetrator2 - ql_other * cv_other_perpetrator2 - p_other * fin_all_perpetrator2 * cv_other_perpetrator2 / bp_perpetrator2 + (p_other - l_other) * ce_other_perpetrator2
    d/dt(other_ew_perpetrator2) <- p_other * fin_all_perpetrator2 * cv_other_perpetrator2 / bp_perpetrator2 - (p_other - l_other) * ce_other_perpetrator2 - l_other * ce_other_perpetrator2 -
      clin_other_perpetrator2 * jin_all_perpetrator2 * ce_other_perpetrator2 * fuint_other_perpetrator2 + clin_other_perpetrator2 * ci_other_perpetrator2 * fucel_other_perpetrator2
    d/dt(other_iw_perpetrator2) <- clin_other_perpetrator2 * jin_all_perpetrator2 * ce_other_perpetrator2 * fuint_other_perpetrator2 - clin_other_perpetrator2 * ci_other_perpetrator2 * fucel_other_perpetrator2
    d/dt(gut_vas_perpetrator2) <- q_gut * c_arterial_perpetrator2 - ql_gut * cv_gut_perpetrator2 - p_gut * fin_all_perpetrator2 * cv_gut_perpetrator2 / bp_perpetrator2 + (p_gut - l_gut) * ce_gut_perpetrator2
    d/dt(gut_ew_perpetrator2) <- p_gut * fin_all_perpetrator2 * cv_gut_perpetrator2 / bp_perpetrator2 - (p_gut - l_gut) * ce_gut_perpetrator2 - l_gut * ce_gut_perpetrator2 +
      clin_gut_perpetrator2 * fgp_perpetrator2 * fucel_gut_perpetrator2 * (cent_duodenum_perpetrator2 + cent_jejunum_perpetrator2 + cent_ileum_perpetrator2 + cent_colon_perpetrator2)
    d/dt(liver_vas_perpetrator2) <- ql_gut * cv_gut_perpetrator2 + ql_spleen * cv_spleen_perpetrator2 + ql_pancreas * cv_pancreas_perpetrator2 + (q_ha + q_by) * c_arterial_perpetrator2 -
      ql_liver * cv_liver_perpetrator2 - p_liver * fin_all_perpetrator2 * cv_liver_perpetrator2 / bp_perpetrator2 + (p_liver - l_liver) * ce_liver_perpetrator2
    d/dt(liver_ew_perpetrator2) <- p_liver * fin_all_perpetrator2 * cv_liver_perpetrator2 / bp_perpetrator2 - (p_liver - l_liver) * ce_liver_perpetrator2 - l_liver * ce_liver_perpetrator2 -
      clin_liver_perpetrator2 * jin_all_perpetrator2 * jin_liver_perpetrator2 * ce_liver_perpetrator2 * fuint_liver_perpetrator2 + clin_liver_perpetrator2 * ci_liver_perpetrator2 * fucel_liver_perpetrator2
    d/dt(liver_iw_perpetrator2) <- clin_liver_perpetrator2 * jin_all_perpetrator2 * jin_liver_perpetrator2 * ce_liver_perpetrator2 * fuint_liver_perpetrator2 - clin_liver_perpetrator2 * ci_liver_perpetrator2 * fucel_liver_perpetrator2 -
      clmet_liver_perpetrator2 * ci_liver_perpetrator2 * fucel_liver_perpetrator2
    d/dt(lnode_vas_perpetrator2) <- q_lnode * c_arterial_perpetrator2 - q_lnode * cv_lnode_perpetrator2 + ltot * ce_lnode_perpetrator2 - ltot * cv_lnode_perpetrator2
    d/dt(lnode_ew_perpetrator2) <- l_lung * ce_lung_perpetrator2 + l_adipose * ce_adipose_perpetrator2 + l_bone * ce_bone_perpetrator2 + l_brain * ce_brain_perpetrator2 +
      l_gonads * ce_gonads_perpetrator2 + l_heart * ce_heart_perpetrator2 + l_kidney * ce_kidney_perpetrator2 + l_muscle * ce_muscle_perpetrator2 +
      l_skin * ce_skin_perpetrator2 + l_thymus * ce_thymus_perpetrator2 + l_gut * ce_gut_perpetrator2 + l_spleen * ce_spleen_perpetrator2 +
      l_pancreas * ce_pancreas_perpetrator2 + l_liver * ce_liver_perpetrator2 + l_other * ce_other_perpetrator2 - ltot * ce_lnode_perpetrator2 -
      clin_lnode_perpetrator2 * jin_all_perpetrator2 * ce_lnode_perpetrator2 * fuint_lnode_perpetrator2 + clin_lnode_perpetrator2 * ci_lnode_perpetrator2 * fucel_lnode_perpetrator2
    d/dt(lnode_iw_perpetrator2) <- clin_lnode_perpetrator2 * jin_all_perpetrator2 * ce_lnode_perpetrator2 * fuint_lnode_perpetrator2 - clin_lnode_perpetrator2 * ci_lnode_perpetrator2 * fucel_lnode_perpetrator2
    d/dt(venous_perpetrator2) <- ql_adipose * cv_adipose_perpetrator2 + ql_bone * cv_bone_perpetrator2 + ql_brain * cv_brain_perpetrator2 +
      ql_gonads * cv_gonads_perpetrator2 + ql_heart * cv_heart_perpetrator2 + ql_kidney * cv_kidney_perpetrator2 +
      ql_muscle * cv_muscle_perpetrator2 + ql_skin * cv_skin_perpetrator2 + ql_thymus * cv_thymus_perpetrator2 +
      ql_liver * cv_liver_perpetrator2 + q_lnode * cv_lnode_perpetrator2 + ql_other * cv_other_perpetrator2 +
      ltot * cv_lnode_perpetrator2 - q_lung * c_venous_perpetrator2
    d/dt(arterial_perpetrator2) <- ql_lung * cv_lung_perpetrator2 - (q_adipose + q_bone + q_brain + q_gonads + q_heart +
      q_kidney + q_muscle + q_skin + q_thymus + q_gut + q_spleen + q_pancreas +
      q_ha + q_by + q_lnode + q_other) * c_arterial_perpetrator2
    d/dt(stomach_perpetrator2) <- -stomach_perpetrator2 / t_stomach
    alag(stomach_perpetrator2) <- lagtime_perpetrator2 # dose rows are placed at start + LagTime (PBPK_StudyDesign.m RegimenFun)
    d/dt(duodenum_perpetrator2) <- stomach_perpetrator2 / t_stomach - duodenum_perpetrator2 / t_duodenum - clab_duodenum_perpetrator2 * duodenum_perpetrator2 / vflu_duodenum
    d/dt(jejunum_perpetrator2) <- duodenum_perpetrator2 / t_duodenum - jejunum_perpetrator2 / t_jejunum - clab_jejunum_perpetrator2 * jejunum_perpetrator2 / vflu_jejunum
    d/dt(ileum_perpetrator2) <- jejunum_perpetrator2 / t_jejunum - ileum_perpetrator2 / t_ileum - clab_ileum_perpetrator2 * ileum_perpetrator2 / vflu_ileum
    d/dt(colon_perpetrator2) <- ileum_perpetrator2 / t_ileum - colon_perpetrator2 / t_colon - clab_colon_perpetrator2 * colon_perpetrator2 / vflu_colon
    d/dt(duodenum_uptake_perpetrator2) <- clab_duodenum_perpetrator2 * duodenum_perpetrator2 / vflu_duodenum - kperup_eff_perpetrator2 * duodenum_uptake_perpetrator2
    d/dt(jejunum_uptake_perpetrator2) <- clab_jejunum_perpetrator2 * jejunum_perpetrator2 / vflu_jejunum - kperup_eff_perpetrator2 * jejunum_uptake_perpetrator2
    d/dt(ileum_uptake_perpetrator2) <- clab_ileum_perpetrator2 * ileum_perpetrator2 / vflu_ileum - kperup_eff_perpetrator2 * ileum_uptake_perpetrator2
    d/dt(colon_uptake_perpetrator2) <- clab_colon_perpetrator2 * colon_perpetrator2 / vflu_colon - kperup_eff_perpetrator2 * colon_uptake_perpetrator2
    d/dt(duodenum_enterocyte_perpetrator2) <- kperup_eff_perpetrator2 * duodenum_uptake_perpetrator2 - clin_gut_perpetrator2 * fgp_perpetrator2 * cent_duodenum_perpetrator2 * fucel_gut_perpetrator2 -
      clmet_duodenum_perpetrator2 * cent_duodenum_perpetrator2 * fucel_gut_perpetrator2
    d/dt(jejunum_enterocyte_perpetrator2) <- kperup_eff_perpetrator2 * jejunum_uptake_perpetrator2 - clin_gut_perpetrator2 * fgp_perpetrator2 * cent_jejunum_perpetrator2 * fucel_gut_perpetrator2 -
      clmet_jejunum_perpetrator2 * cent_jejunum_perpetrator2 * fucel_gut_perpetrator2
    d/dt(ileum_enterocyte_perpetrator2) <- kperup_eff_perpetrator2 * ileum_uptake_perpetrator2 - clin_gut_perpetrator2 * fgp_perpetrator2 * cent_ileum_perpetrator2 * fucel_gut_perpetrator2 -
      clmet_ileum_perpetrator2 * cent_ileum_perpetrator2 * fucel_gut_perpetrator2
    d/dt(colon_enterocyte_perpetrator2) <- kperup_eff_perpetrator2 * colon_uptake_perpetrator2 - clin_gut_perpetrator2 * fgp_perpetrator2 * cent_colon_perpetrator2 * fucel_gut_perpetrator2
    d/dt(a_feces_perpetrator2) <- colon_perpetrator2 / t_colon
    # Enzyme turnover, relative to baseline abundance (1 = baseline).
    d/dt(enzyme_cyp3a4_liver_perpetrator2) <- 0.0077 * (1 + ind_cyp3a4_perpetrator2) - (0.0077 + mbi_cyp3a4_perpetrator2) * enzyme_cyp3a4_liver_perpetrator2
    d/dt(enzyme_cyp2c19_liver_perpetrator2) <- 0.0267 * (1 + ind_cyp2c19_perpetrator2) - (0.0267) * enzyme_cyp2c19_liver_perpetrator2
    d/dt(enzyme_cyp2d6_liver_perpetrator2) <- 0.0143 * (1 + ind_cyp2d6_perpetrator2) - (0.0143) * enzyme_cyp2d6_liver_perpetrator2
    d/dt(enzyme_cyp2c8_liver_perpetrator2) <- 0.0310 * (1 + ind_cyp2c8_perpetrator2) - (0.0310) * enzyme_cyp2c8_liver_perpetrator2
    d/dt(enzyme_cyp1a2_liver_perpetrator2) <- 0.0149 * (1 + ind_cyp1a2_perpetrator2) - (0.0149) * enzyme_cyp1a2_liver_perpetrator2
    d/dt(enzyme_cyp2a6_liver_perpetrator2) <- 0.0267 * (1 + ind_cyp2a6_perpetrator2) - (0.0267) * enzyme_cyp2a6_liver_perpetrator2
    d/dt(enzyme_cyp2b6_liver_perpetrator2) <- 0.0217 * (1 + ind_cyp2b6_perpetrator2) - (0.0217) * enzyme_cyp2b6_liver_perpetrator2
    d/dt(enzyme_cyp2j2_liver_perpetrator2) <- 0.01 * (1 + ind_cyp2j2_perpetrator2) - (0.01 + mbi_cyp2j2_perpetrator2) * enzyme_cyp2j2_liver_perpetrator2
    d/dt(enzyme_ugt1a1_liver_perpetrator2) <- 0.0693 * (1 + indu_ugt1a1_perpetrator2) - 0.0693 * enzyme_ugt1a1_liver_perpetrator2
    d/dt(enzyme_cyp3a4_duodenum_perpetrator2) <- 0.03 * (1 + indg_du_perpetrator2) - (0.03 + mbig_du_perpetrator2) * enzyme_cyp3a4_duodenum_perpetrator2
    d/dt(enzyme_cyp3a4_jejunum_perpetrator2) <- 0.03 * (1 + indg_je_perpetrator2) - (0.03 + mbig_je_perpetrator2) * enzyme_cyp3a4_jejunum_perpetrator2
    d/dt(enzyme_cyp3a4_ileum_perpetrator2) <- 0.03 * (1 + indg_il_perpetrator2) - (0.03 + mbig_il_perpetrator2) * enzyme_cyp3a4_ileum_perpetrator2
    enzyme_cyp3a4_liver_perpetrator2(0) <- 1
    enzyme_cyp2c19_liver_perpetrator2(0) <- 1
    enzyme_cyp2d6_liver_perpetrator2(0) <- 1
    enzyme_cyp2c8_liver_perpetrator2(0) <- 1
    enzyme_cyp1a2_liver_perpetrator2(0) <- 1
    enzyme_cyp2a6_liver_perpetrator2(0) <- 1
    enzyme_cyp2b6_liver_perpetrator2(0) <- 1
    enzyme_cyp2j2_liver_perpetrator2(0) <- 1
    enzyme_ugt1a1_liver_perpetrator2(0) <- 1
    enzyme_cyp3a4_duodenum_perpetrator2(0) <- 1
    enzyme_cyp3a4_jejunum_perpetrator2(0) <- 1
    enzyme_cyp3a4_ileum_perpetrator2(0) <- 1

    # Reported concentration (ng/mL): the venous-blood state times 1000, as
    # PBPK_ExtractConcentration.m reports CONC.VB * MW with no blood:plasma
    # conversion.
    Cc <- 1000 * c_venous
    Cc_perpetrator <- 1000 * c_venous_perpetrator
    Cc_perpetrator2 <- 1000 * c_venous_perpetrator2
  })
}
