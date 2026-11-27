Nakamura_2018_tamoxifen_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, five-compartment dispersion liver, segregated-flow",
    "intestine, enterohepatic circulation). Tamoxifen and its four",
    "metabolites N-desmethyltamoxifen (NDMT), 4-hydroxytamoxifen (4-OHT),",
    "4-OHT-glucuronide (4-OHT-G) and endoxifen after oral tamoxifen in adult",
    "postmenopausal women with breast cancer. Five structurally identical",
    "PBPK modules -- one per chemical species -- are coupled by hepatic and",
    "intestinal metabolic clearances, so the whole system is one jointly",
    "fitted model with 88 ODE states. Each module carries a central (blood +",
    "rapidly equilibrating tissue) compartment, perfusion-limited muscle,",
    "skin, adipose and gut serosa, a liver resolved as five well-stirred",
    "sub-compartments in series (a dispersion-model approximation), a",
    "three-compartment biliary transit chain, a segregated-flow intestine",
    "(lumen / enterocyte / mucosal blood) and faecal and urinary sinks. The",
    "compound-specific kinetic and metabolic parameters were fitted with the",
    "Cluster Newton method (CNM), which returns a family of equally",
    "well-fitting parameter vectors rather than a single optimum; the values",
    "below are set ID 1, the authors' deposited virtual-clinical-study worked",
    "example, and the vignette shows how to substitute the other seven (set",
    "ID 3 was the authors' Grubbs outlier and was dropped from their",
    "analyses). The CYP2D6, CYP3A4 and CYP2C9 covariates are relative",
    "enzyme-activity multipliers (default 1 = wild-type) scaling the",
    "genotype-dependent formation clearances; the virtual TARGET-1 study",
    "varies CYP2D6 across diplotypes. Between-subject variability is the",
    "normal physiological variability and lognormal kinetic / enzyme",
    "variability of Supplementary Table S3. The paper reports no",
    "residual-error model (CNM is a fixed-effects least-squares method), so",
    "the propSd terms are placeholders. The clinical-response layer (the",
    "stepwise endoxifen x tumoral-ER probability model and the virtual",
    "TARGET-1 trial) is reproduced in the vignette as post-processing, not as",
    "ODE states. See the vignette Errata for the mass-balance corrections",
    "applied to the published equations and the single-chain-liver reading."
  )
  reference <- paste(
    "Nakamura T, Toshimoto K, Lee W, Imamura CK, Tanigawara Y, Sugiyama Y.",
    "Application of PBPK Modeling and Virtual Clinical Study Approaches to",
    "Predict the Outcomes of CYP2D6 Genotype-Guided Dosing of Tamoxifen.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(7):474-482.",
    "doi:10.1002/psp4.12307.",
    "The ordinary differential equations are transcribed from the",
    "Supplementary Methods (supplementary information 1). Fixed physiological",
    "and physicochemical constants are Supplementary Table S1; the eight CNM",
    "parameter sets are Supplementary Table S2 (set ID 1 used here); the",
    "between-subject variability is Supplementary Table S3; the CYP2D6",
    "activity-based diplotype distribution is Table 1. The deposited MATLAB",
    "CNM code and R virtual-clinical-study code (supplementary information)",
    "were read to settle the model structure and the ID 1 values.",
    sep = " "
  )
  vignette <- "Nakamura_2018_tamoxifen_pbpk"

  # Compartment vocabulary. `central`, `muscle`, `skin`, `adipose`,
  # `gut_lumen`, `a_feces` and `a_urine` are canonical
  # (inst/references/compartment-names.md); the metabolite suffixes
  # `ndmtam` (N-desmethyltamoxifen), `ohtam` (4-hydroxytamoxifen),
  # `ohtamg` (4-OHT-glucuronide) and `endox` (endoxifen) are registered in
  # that same file, and tamoxifen itself is the parent analyte so its states
  # carry no suffix. The remaining stems are declared paper-specific:
  # `serosa` is the gut serosal tissue, `liver1`..`liver5` the five
  # well-stirred hepatic dispersion sub-compartments (a SINGLE chain -- this
  # model has no separate hepatocyte/sinusoid split and no active hepatic
  # uptake, unlike the Toshimoto 2017 irinotecan sibling), `ehc1`..`ehc3`
  # the biliary transit chain, and `intestine_ent` / `intestine_muc` the
  # absorptive-enterocyte and mucosal-blood compartments of the
  # segregated-flow intestine.
  paper_specific_compartments <- c(
    "serosa",
    "liver1",
    "liver2",
    "liver3",
    "liver4",
    "liver5",
    "ehc1",
    "ehc2",
    "ehc3",
    "intestine_ent",
    "intestine_muc",
    "serosa_ndmtam",
    "liver1_ndmtam",
    "liver2_ndmtam",
    "liver3_ndmtam",
    "liver4_ndmtam",
    "liver5_ndmtam",
    "ehc1_ndmtam",
    "ehc2_ndmtam",
    "ehc3_ndmtam",
    "intestine_ent_ndmtam",
    "intestine_muc_ndmtam",
    "serosa_ohtam",
    "liver1_ohtam",
    "liver2_ohtam",
    "liver3_ohtam",
    "liver4_ohtam",
    "liver5_ohtam",
    "ehc1_ohtam",
    "ehc2_ohtam",
    "ehc3_ohtam",
    "intestine_ent_ohtam",
    "intestine_muc_ohtam",
    "serosa_ohtamg",
    "liver1_ohtamg",
    "liver2_ohtamg",
    "liver3_ohtamg",
    "liver4_ohtamg",
    "liver5_ohtamg",
    "ehc1_ohtamg",
    "ehc2_ohtamg",
    "ehc3_ohtamg",
    "intestine_ent_ohtamg",
    "intestine_muc_ohtamg",
    "serosa_endox",
    "liver1_endox",
    "liver2_endox",
    "liver3_endox",
    "liver4_endox",
    "liver5_endox",
    "ehc1_endox",
    "ehc2_endox",
    "ehc3_endox",
    "intestine_ent_endox",
    "intestine_muc_endox"
  )

  # Time in hours. Every volume, flow and clearance in Supplementary Tables
  # S1 and S2 is per kilogram of body weight and is multiplied by WT in
  # model(), so the tissue / blood states hold concentrations (nmol/L = nM)
  # and the luminal, biliary, faecal and urinary states hold amounts (nmol).
  # Tamoxifen doses are supplied in nmol (tamoxifen base MW 371.51 g/mol, so
  # a 20 mg daily dose is 20000 / 371.51 = 53834 nmol); the endoxifen blood
  # concentration Cc_endox is therefore directly comparable to the paper's
  # 14.3 and 35.1 nM response thresholds.
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Central and mucosal-blood states are whole-blood
  # concentrations (nM); the other tissue and liver states are tissue
  # concentrations (nM); the luminal, biliary, faecal and urinary states
  # are amounts (nmol).
  compartmentData <- list(
    central = list(analyte = "tamoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    muscle = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    serosa = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver1 = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver2 = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver3 = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver4 = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver5 = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    ehc1 = list(analyte = "tamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2 = list(analyte = "tamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3 = list(analyte = "tamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    gut_lumen = list(analyte = "tamoxifen", units = "nmol", specimen = "administration site", verified = TRUE),
    intestine_ent = list(analyte = "tamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    intestine_muc = list(analyte = "tamoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    a_feces = list(analyte = "tamoxifen", units = "nmol", specimen = "faeces", verified = TRUE),
    a_urine = list(analyte = "tamoxifen", units = "nmol", specimen = "urine", verified = TRUE),
    central_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    muscle_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    skin_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    adipose_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    serosa_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver1_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver2_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver3_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver4_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver5_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    ehc1_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    gut_lumen_ndmtam = list(
      analyte = "N-desmethyltamoxifen",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    intestine_ent_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    intestine_muc_ndmtam = list(
      analyte = "N-desmethyltamoxifen",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    a_feces_ndmtam = list(analyte = "N-desmethyltamoxifen", units = "nmol", specimen = "faeces", verified = TRUE),
    central_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    muscle_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    skin_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    adipose_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    serosa_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver1_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver2_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver3_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver4_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver5_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    ehc1_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    gut_lumen_ohtam = list(
      analyte = "4-hydroxytamoxifen",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    intestine_ent_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    intestine_muc_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    a_feces_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nmol", specimen = "faeces", verified = TRUE),
    a_urine_ohtam = list(analyte = "4-hydroxytamoxifen", units = "nmol", specimen = "urine", verified = TRUE),
    central_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    muscle_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    skin_ohtamg = list(analyte = "4-hydroxytamoxifen glucuronide", units = "nM", specimen = "tissue", verified = TRUE),
    adipose_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    serosa_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver1_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver2_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver3_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver4_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver5_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    ehc1_ohtamg = list(analyte = "4-hydroxytamoxifen glucuronide", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2_ohtamg = list(analyte = "4-hydroxytamoxifen glucuronide", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3_ohtamg = list(analyte = "4-hydroxytamoxifen glucuronide", units = "nmol", specimen = "bile", verified = TRUE),
    gut_lumen_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    intestine_ent_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    intestine_muc_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    a_urine_ohtamg = list(
      analyte = "4-hydroxytamoxifen glucuronide",
      units = "nmol",
      specimen = "urine",
      verified = TRUE
    ),
    central_endox = list(analyte = "endoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    muscle_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    skin_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    adipose_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    serosa_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver1_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver2_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver3_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver4_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    liver5_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    ehc1_endox = list(analyte = "endoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2_endox = list(analyte = "endoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3_endox = list(analyte = "endoxifen", units = "nmol", specimen = "bile", verified = TRUE),
    gut_lumen_endox = list(analyte = "endoxifen", units = "nmol", specimen = "administration site", verified = TRUE),
    intestine_ent_endox = list(analyte = "endoxifen", units = "nM", specimen = "tissue", verified = TRUE),
    intestine_muc_endox = list(analyte = "endoxifen", units = "nM", specimen = "whole blood", verified = TRUE),
    a_feces_endox = list(analyte = "endoxifen", units = "nmol", specimen = "faeces", verified = TRUE),
    a_urine_endox = list(analyte = "endoxifen", units = "nmol", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = paste(
        "Body weight. Scales every per-kilogram volume, flow and clearance",
        "to absolute units inside model(). Reference 67.58 kg (the mean of",
        "the US postmenopausal-women cohort the model was fitted to;",
        "Supplementary Table S1 gives Caucasian 66.39 / 67.58 kg and",
        "Japanese 53.85 / 54.44 kg)."
      ),
      units = "kg",
      type = "continuous",
      source_name = "BW"
    ),
    CYP2D6 = list(
      description = paste(
        "Relative CYP2D6 enzyme activity multiplier, 1 = wild-type *1/*1.",
        "Scales the CYP2D6-mediated formation clearances tamoxifen ->",
        "4-OHT and NDMT -> endoxifen. The virtual TARGET-1 study draws this",
        "per subject from the activity-based diplotype distribution of",
        "Table 1 (relative activity 1.13 increase/normal down to 0.004",
        "null/null, with diplotype-specific %CV)."
      ),
      units = "unitless (relative activity, reference 1 = CYP2D6 *1/*1)",
      type = "continuous",
      source_name = "CYP2D6",
      notes = paste(
        "Table 1 relative activity scores (reference CYP2D6*1/*1 = 1):",
        "increase/normal 1.13 (CV 64.3%), normal/normal 0.89 (47.9%),",
        "normal/decrease 0.71 (53.6%), normal/null 0.42 (99.0%),",
        "decrease/decrease 0.24 (66.0%), decrease/null 0.20 (109.0%),",
        "null/null 0.004 (103.0%). The %CV is the published between-subject",
        "variability of CYP2D6 activity WITHIN a diplotype and is applied in",
        "the vignette as the per-subject draw, so no separate IIV eta sits",
        "on the CYP2D6 pathways in the model."
      )
    ),
    CYP3A4 = list(
      description = paste(
        "Relative CYP3A4 enzyme activity multiplier, 1 = reference. Scales",
        "the CYP3A4-mediated formation and elimination clearances (tamoxifen",
        "-> NDMT, NDMT -> other, 4-OHT -> endoxifen, and the 4-OHT 3A",
        "elimination routes). Held at the reference value 1 throughout the",
        "virtual TARGET-1 study, which genotypes only CYP2D6; the 33% CV",
        "lognormal CYP3A variability of Supplementary Table S3 is carried by",
        "the model's IIV etas, not by this covariate."
      ),
      units = "unitless (relative activity, reference 1)",
      type = "continuous",
      source_name = "CYP3A4"
    ),
    CYP2C9 = list(
      description = paste(
        "Relative CYP2C9 enzyme activity multiplier, 1 = wild-type *1/*1.",
        "Scales the CYP2C9-mediated tamoxifen -> 4-OHT formation clearance.",
        "Held at the reference value 1 throughout the virtual TARGET-1",
        "study, which genotypes only CYP2D6."
      ),
      units = "unitless (relative activity, reference 1 = CYP2C9 *1/*1)",
      type = "continuous",
      source_name = "CYP2C9"
    )
  )

  population <- list(
    species = "human",
    n_subjects = paste(
      "The compound-specific kinetic and metabolic parameters were fitted",
      "by the Cluster Newton method to the mean blood concentration-time",
      "profiles of tamoxifen, NDMT, 4-OHT and endoxifen in 30 healthy",
      "postmenopausal women in the United States (Nichols et al. 2014),",
      "converted from plasma to blood with the tabulated blood-to-plasma",
      "ratios. The virtual TARGET-1 study simulated 180 virtual Japanese",
      "breast-cancer patients per run (200 runs), with sensitivity runs at",
      "n = 140, 220 and 260."
    ),
    disease_state = "hormone receptor-positive metastatic or recurrent breast cancer (virtual target population); healthy postmenopausal women (PK calibration cohort)",
    dose_range = paste(
      "Oral tamoxifen 20 mg/day (standard) or 40 mg/day (escalated) to",
      "steady state. 20 mg tamoxifen base (MW 371.51) is 53834 nmol."
    ),
    weight_range = "mean 67.58 kg (US / Caucasian) and 54.44 kg (Japanese); lognormal, CV 19 percent (Caucasian) / 15 percent (Japanese), Supplementary Tables S1 and S3",
    sex_female_pct = 100,
    race_ethnicity = "US / Caucasian calibration cohort; Japanese virtual target population (ethnic differences in body weight, CYP2D6 diplotype frequency and tumoral ER expression, Table 1)",
    regions = "United States (PK calibration); Japan (virtual TARGET-1 target population)",
    notes = paste(
      "Ages are 45-65 years (the SimCYP window used to derive the body",
      "weight distribution). The CYP2D6 activity-based diplotype",
      "frequencies differ by ethnicity (Table 1): e.g. null/null 6.2",
      "percent in Caucasians versus 0.5 percent in Japanese, and",
      "normal/decrease 11.2 versus 37.3 percent. Tumoral ER expression",
      "(moderate >=10 to <50 percent versus strong >=50 percent) enters the",
      "clinical-response model in the vignette, not the PBPK ODE system."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # 1. Cluster Newton method (CNM) parameters -- Supplementary Table S2,
    #    parameter set ID 1 (the authors' deposited virtual-clinical-study
    #    worked example). The CNM returns eight parameter vectors that fit
    #    the observed profiles equally well; the vignette tabulates all eight
    #    and shows how to substitute the other seven (ID 3 was the authors'
    #    Grubbs outlier). These are estimates, not fixed constants, so they
    #    carry no fixed() wrapper. SFKp is the per-compound scaling factor
    #    applied to every tissue partition coefficient.
    # ---------------------------------------------------------------

    # -- Tamoxifen (Supplementary Table S2a, row ID 1).
    lv_central <- log(0.410)
    label("Tamoxifen central (blood + rapid tissue) volume (L/kg)") # Suppl Table S2a ID 1: Vcentral 0.410
    lsfkp <- log(5.92)
    label("Tamoxifen Kp scaling factor SFKp (unitless)") # Suppl Table S2a ID 1: SFKp 5.92
    lka <- log(3.70)
    label("Tamoxifen lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table S2a ID 1: ka 3.70
    lkbile <- log(1.49)
    label("Tamoxifen biliary transit rate constant (1/h)") # Suppl Table S2a ID 1: kbile 1.49
    lps_dif <- log(0.527)
    label("Tamoxifen enterocyte passive-diffusion clearance PSdif (L/h/kg)") # Suppl Table S2a ID 1: PSdif 0.527
    lkfeces <- log(0.0224)
    label("Tamoxifen faecal transit rate constant (1/h)") # Suppl Table S2a ID 1: kfeces 0.0224
    lcl_renal <- log(0.00583)
    label("Tamoxifen renal clearance from blood (L/h/kg)") # Suppl Table S2a ID 1: CLR 0.00583
    lcl_bile <- log(7.02)
    label("Tamoxifen biliary intrinsic clearance (L/h/kg)") # Suppl Table S2a ID 1: CLbile 7.02
    lcl_ndmt_hep <- log(1.59)
    label("Hepatic intrinsic clearance, tamoxifen to NDMT via CYP3A (L/h/kg)") # Suppl Table S2a ID 1: CL NDMTby3A,hep 1.59
    lcl_4oht_2d6_hep <- log(0.233)
    label("Hepatic intrinsic clearance, tamoxifen to 4-OHT via CYP2D6 (L/h/kg)") # Suppl Table S2a ID 1: CL 4OHTby2D6,hep 0.233
    lcl_4oht_2c9_hep <- log(0.174)
    label("Hepatic intrinsic clearance, tamoxifen to 4-OHT via CYP2C9 (L/h/kg)") # Suppl Table S2a ID 1: CL 4OHTby2C9,hep 0.174
    lcl_4oht_other_hep <- log(0.0335)
    label("Hepatic intrinsic clearance, tamoxifen to 4-OHT via other enzymes (L/h/kg)") # Suppl Table S2a ID 1: CL 4OHTbyOther,hep 0.0335
    lcl_other_hep <- log(0.966)
    label("Hepatic intrinsic clearance, tamoxifen to other products (L/h/kg)") # Suppl Table S2a ID 1: CL Other,hep 0.966
    lcl_ndmt_ent <- log(0.0281)
    label("Intestinal intrinsic clearance, tamoxifen to NDMT via CYP3A (L/h/kg)") # Suppl Table S2a ID 1: CL NDMTby3A,ent 0.0281

    # -- NDMT (Supplementary Table S2b, row ID 1). Renal clearance fixed to
    #    zero (Supplementary Table S1: CLR '(0)'; renal excretion of NDMT
    #    deemed negligible from a human mass-balance study).
    lv_central_ndmtam <- log(0.129)
    label("NDMT central (blood + rapid tissue) volume (L/kg)") # Suppl Table S2b ID 1: Vcentral 0.129
    lsfkp_ndmtam <- log(1.01)
    label("NDMT Kp scaling factor SFKp (unitless)") # Suppl Table S2b ID 1: SFKp 1.01
    lka_ndmtam <- log(1.50)
    label("NDMT lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table S2b ID 1: ka 1.50
    lkbile_ndmtam <- log(0.655)
    label("NDMT biliary transit rate constant (1/h)") # Suppl Table S2b ID 1: kbile 0.655
    lps_dif_ndmtam <- log(1.06)
    label("NDMT enterocyte passive-diffusion clearance PSdif (L/h/kg)") # Suppl Table S2b ID 1: PSdif 1.06
    lkfeces_ndmtam <- log(0.0681)
    label("NDMT faecal transit rate constant (1/h)") # Suppl Table S2b ID 1: kfeces 0.0681
    lcl_bile_ndmtam <- log(0.0457)
    label("NDMT biliary intrinsic clearance (L/h/kg)") # Suppl Table S2b ID 1: CLbile 0.0457
    lcl_end_2d6_hep_ndmtam <- log(0.263)
    label("Hepatic intrinsic clearance, NDMT to endoxifen via CYP2D6 (L/h/kg)") # Suppl Table S2b ID 1: CL ENDby2D6,hep 0.263
    lcl_other_3a4_hep_ndmtam <- log(0.242)
    label("Hepatic intrinsic clearance, NDMT to other products via CYP3A (L/h/kg)") # Suppl Table S2b ID 1: CL Otherby3A,hep 0.242
    lcl_other_3a4_ent_ndmtam <- log(0.199)
    label("Intestinal intrinsic clearance, NDMT to other products via CYP3A (L/h/kg)") # Suppl Table S2b ID 1: CL Otherby3A,ent 0.199

    # -- 4-OHT (Supplementary Table S2c, row ID 1).
    lv_central_ohtam <- log(0.0876)
    label("4-OHT central (blood + rapid tissue) volume (L/kg)") # Suppl Table S2c ID 1: Vcentral 0.0876
    lsfkp_ohtam <- log(1.01)
    label("4-OHT Kp scaling factor SFKp (unitless)") # Suppl Table S2c ID 1: SFKp 1.01
    lka_ohtam <- log(0.651)
    label("4-OHT lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table S2c ID 1: ka 0.651
    lkbile_ohtam <- log(0.105)
    label("4-OHT biliary transit rate constant (1/h)") # Suppl Table S2c ID 1: kbile 0.105
    lps_dif_ohtam <- log(0.759)
    label("4-OHT enterocyte passive-diffusion clearance PSdif (L/h/kg)") # Suppl Table S2c ID 1: PSdif 0.759
    lkfeces_ohtam <- log(0.187)
    label("4-OHT faecal transit rate constant (1/h)") # Suppl Table S2c ID 1: kfeces 0.187
    lcl_renal_ohtam <- log(0.000132)
    label("4-OHT renal clearance from blood (L/h/kg)") # Suppl Table S2c ID 1: CLR 0.000132
    lcl_bile_ohtam <- log(5.68)
    label("4-OHT biliary intrinsic clearance (L/h/kg)") # Suppl Table S2c ID 1: CLbile 5.68
    lcl_end_3a4_hep_ohtam <- log(12.1)
    label("Hepatic intrinsic clearance, 4-OHT to endoxifen via CYP3A (L/h/kg)") # Suppl Table S2c ID 1: CL ENDby3A,hep 12.1
    lcl_4ohtg_ugt_hep_ohtam <- log(0.489)
    label("Hepatic intrinsic clearance, 4-OHT to 4-OHT-G via UGT (L/h/kg)") # Suppl Table S2c ID 1: CL 4OHTGbyUGT,hep 0.489
    lcl_other_3a4_hep_ohtam <- log(1.07)
    label("Hepatic intrinsic clearance, 4-OHT to other products via CYP3A (L/h/kg)") # Suppl Table S2c ID 1: CL Otherby3A,hep 1.07
    lcl_other_hep_ohtam <- log(0.219)
    label("Hepatic intrinsic clearance, 4-OHT to other products (L/h/kg)") # Suppl Table S2c ID 1: CL Other,hep 0.219
    lcl_end_3a4_ent_ohtam <- log(1.78)
    label("Intestinal intrinsic clearance, 4-OHT to endoxifen via CYP3A (L/h/kg)") # Suppl Table S2c ID 1: CL ENDby3A,ent 1.78
    lcl_4ohtg_ugt_ent_ohtam <- log(0.469)
    label("Intestinal intrinsic clearance, 4-OHT to 4-OHT-G via UGT (L/h/kg)") # Suppl Table S2c ID 1: CL 4OHTGbyUGT,ent 0.469
    lcl_other_3a4_ent_ohtam <- log(0.833)
    label("Intestinal intrinsic clearance, 4-OHT to other products via CYP3A (L/h/kg)") # Suppl Table S2c ID 1: CL Otherby3A,ent 0.833

    # -- 4-OHT-glucuronide (Supplementary Table S2d, row ID 1). No
    #    metabolism; kfeces fixed to zero (Supplementary Table S1: kfeces
    #    '(0)'). kdeg is the luminal deconjugation rate 4-OHT-G -> 4-OHT.
    lv_central_ohtamg <- log(0.0812)
    label("4-OHT-G central (blood + rapid tissue) volume (L/kg)") # Suppl Table S2d ID 1: Vcentral 0.0812
    lsfkp_ohtamg <- log(1.45)
    label("4-OHT-G Kp scaling factor SFKp (unitless)") # Suppl Table S2d ID 1: SFKp 1.45
    lka_ohtamg <- log(3.28)
    label("4-OHT-G lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table S2d ID 1: ka 3.28
    lkbile_ohtamg <- log(2.62)
    label("4-OHT-G biliary transit rate constant (1/h)") # Suppl Table S2d ID 1: kbile 2.62
    lps_dif_ohtamg <- log(1.12)
    label("4-OHT-G enterocyte passive-diffusion clearance PSdif (L/h/kg)") # Suppl Table S2d ID 1: PSdif 1.12
    lcl_renal_ohtamg <- log(0.00701)
    label("4-OHT-G renal clearance from blood (L/h/kg)") # Suppl Table S2d ID 1: CLR 0.00701
    lcl_bile_ohtamg <- log(36.4)
    label("4-OHT-G biliary intrinsic clearance (L/h/kg)") # Suppl Table S2d ID 1: CLbile 36.4
    lkdeg_ohtamg <- log(0.769)
    label("Luminal deconjugation rate constant, 4-OHT-G to 4-OHT (1/h)") # Suppl Table S2d ID 1: kdeg 0.769

    # -- Endoxifen (Supplementary Table S2e, row ID 1). Endoxifen is the one
    #    compound with active intestinal efflux (P-glycoprotein), PSact,ent.
    lv_central_endox <- log(0.0779)
    label("Endoxifen central (blood + rapid tissue) volume (L/kg)") # Suppl Table S2e ID 1: Vcentral 0.0779
    lsfkp_endox <- log(1.04)
    label("Endoxifen Kp scaling factor SFKp (unitless)") # Suppl Table S2e ID 1: SFKp 1.04
    lka_endox <- log(0.427)
    label("Endoxifen lumen-to-enterocyte absorption rate constant (1/h)") # Suppl Table S2e ID 1: ka 0.427
    lkbile_endox <- log(4.13)
    label("Endoxifen biliary transit rate constant (1/h)") # Suppl Table S2e ID 1: kbile 4.13
    lps_dif_endox <- log(0.0570)
    label("Endoxifen enterocyte passive-diffusion clearance PSdif (L/h/kg)") # Suppl Table S2e ID 1: PSdif 0.0570
    lps_act_ent_endox <- log(1.10)
    label("Endoxifen enterocyte active (Pgp) efflux clearance PSact,ent (L/h/kg)") # Suppl Table S2e ID 1: PSact,ent 1.10
    lcl_renal_endox <- log(0.00110)
    label("Endoxifen renal clearance from blood (L/h/kg)") # Suppl Table S2e ID 1: CLR 0.00110
    lkfeces_endox <- log(0.246)
    label("Endoxifen faecal transit rate constant (1/h)") # Suppl Table S2e ID 1: kfeces 0.246
    lcl_bile_endox <- log(4.20)
    label("Endoxifen biliary intrinsic clearance (L/h/kg)") # Suppl Table S2e ID 1: CLbile 4.20
    lcl_other_hep_endox <- log(0.134)
    label("Hepatic intrinsic clearance, endoxifen to other products (L/h/kg)") # Suppl Table S2e ID 1: CL Other,hep 0.134
    lcl_other_ent_endox <- log(31.4)
    label("Intestinal intrinsic clearance, endoxifen to other products (L/h/kg)") # Suppl Table S2e ID 1: CL Other,ent 31.4

    # ---------------------------------------------------------------
    # 2. Fixed physiological constants -- Supplementary Table S1, scaled per
    #    kilogram of body weight. Shared by all five compounds.
    # ---------------------------------------------------------------
    lq_liver <- fixed(log(1.243))
    label("Hepatic blood flow Qhep (L/h/kg)") # Suppl Table S1: Qhep 1.243
    lq_muscle <- fixed(log(0.643))
    label("Muscle blood flow (L/h/kg)") # Suppl Table S1: Qmus 0.643
    lq_skin <- fixed(log(0.257))
    label("Skin blood flow (L/h/kg)") # Suppl Table S1: Qskin 0.257
    lq_adipose <- fixed(log(0.223))
    label("Adipose blood flow (L/h/kg)") # Suppl Table S1: Qadi 0.223
    lq_mucosa <- fixed(log(0.257))
    label("Intestinal mucosal blood flow Qmuc,b (L/h/kg)") # Suppl Table S1: Qmuc,b 0.257
    lq_serosa <- fixed(log(0.274))
    label("Intestinal serosal blood flow Qser (L/h/kg)") # Suppl Table S1: Qser 0.274
    lv_liver <- fixed(log(0.0241))
    label("Total liver volume Vliver (L/kg)") # Suppl Table S1: Vliver 0.0241
    lv_muscle <- fixed(log(0.429))
    label("Muscle volume (L/kg)") # Suppl Table S1: Vmus 0.429
    lv_skin <- fixed(log(0.111))
    label("Skin volume (L/kg)") # Suppl Table S1: Vskin 0.111
    lv_adipose <- fixed(log(0.143))
    label("Adipose volume (L/kg)") # Suppl Table S1: Vadi 0.143
    lv_ent <- fixed(log(0.00739))
    label("Enterocyte volume Vent (L/kg)") # Suppl Table S1: Vent 0.00739
    lv_mucb <- fixed(log(0.00099))
    label("Mucosal-blood volume Vmuc,b (L/kg)") # Suppl Table S1: Vmuc,b 0.00099
    lv_serosa <- fixed(log(0.00893))
    label("Serosa volume Vser (L/kg)") # Suppl Table S1: Vser 0.00893
    ar <- fixed(20)
    label("Enterocyte apical-to-basolateral surface-area ratio AR (unitless)") # Suppl Table S1: AR 20

    # ---------------------------------------------------------------
    # 3. Fixed physicochemical constants -- Supplementary Table S1B. fp is
    #    the unbound fraction in plasma, RB the blood-to-plasma ratio, Kp the
    #    tissue-to-blood partition coefficients and fgut the enterocyte free
    #    fraction. The blood unbound fraction fb = fp / RB and the liver
    #    free fraction fliver = fb / (Kp,liver * SFKp) are derived in model().
    # ---------------------------------------------------------------
    fp <- fixed(0.0157)
    label("Tamoxifen unbound fraction in plasma (unitless)") # Suppl Table S1B: fp 0.0157
    rb <- fixed(0.803)
    label("Tamoxifen blood-to-plasma ratio RB (unitless)") # Suppl Table S1B: RB 0.803
    lkp_liver <- fixed(log(5.23))
    label("Tamoxifen liver-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,liver 5.23
    lkp_muscle <- fixed(log(3.03))
    label("Tamoxifen muscle-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,mus 3.03
    lkp_skin <- fixed(log(14.8))
    label("Tamoxifen skin-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,skin 14.8
    lkp_adipose <- fixed(log(36.2))
    label("Tamoxifen adipose-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,adi 36.2
    lkp_gut <- fixed(log(10.0))
    label("Tamoxifen gut-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,gut 10.0
    fgut <- fixed(0.00214)
    label("Tamoxifen enterocyte free fraction fgut (unitless)") # Suppl Table S1B: fgut 0.00214
    fp_ndmtam <- fixed(0.0198)
    label("NDMT unbound fraction in plasma (unitless)") # Suppl Table S1B: fp 0.0198
    rb_ndmtam <- fixed(0.815)
    label("NDMT blood-to-plasma ratio RB (unitless)") # Suppl Table S1B: RB 0.815
    lkp_liver_ndmtam <- fixed(log(9.52))
    label("NDMT liver-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,liver 9.52
    lkp_muscle_ndmtam <- fixed(log(3.36))
    label("NDMT muscle-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,mus 3.36
    lkp_skin_ndmtam <- fixed(log(4.19))
    label("NDMT skin-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,skin 4.19
    lkp_adipose_ndmtam <- fixed(log(4.44))
    label("NDMT adipose-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,adi 4.44
    lkp_gut_ndmtam <- fixed(log(5.82))
    label("NDMT gut-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,gut 5.82
    fgut_ndmtam <- fixed(0.00476)
    label("NDMT enterocyte free fraction fgut (unitless)") # Suppl Table S1B: fgut 0.00476
    fp_ohtam <- fixed(0.0135)
    label("4-OHT unbound fraction in plasma (unitless)") # Suppl Table S1B: fp 0.0135
    rb_ohtam <- fixed(0.880)
    label("4-OHT blood-to-plasma ratio RB (unitless)") # Suppl Table S1B: RB 0.880
    lkp_liver_ohtam <- fixed(log(11.2))
    label("4-OHT liver-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,liver 11.2
    lkp_muscle_ohtam <- fixed(log(4.03))
    label("4-OHT muscle-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,mus 4.03
    lkp_skin_ohtam <- fixed(log(6.06))
    label("4-OHT skin-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,skin 6.06
    lkp_adipose_ohtam <- fixed(log(7.45))
    label("4-OHT adipose-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,adi 7.45
    lkp_gut_ohtam <- fixed(log(7.45))
    label("4-OHT gut-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,gut 7.45
    fgut_ohtam <- fixed(0.00224)
    label("4-OHT enterocyte free fraction fgut (unitless)") # Suppl Table S1B: fgut 0.00224
    fp_ohtamg <- fixed(0.0365)
    label("4-OHT-G unbound fraction in plasma (unitless)") # Suppl Table S1B: fp 0.0365
    rb_ohtamg <- fixed(0.648)
    label("4-OHT-G blood-to-plasma ratio RB (unitless)") # Suppl Table S1B: RB 0.648
    lkp_liver_ohtamg <- fixed(log(0.185))
    label("4-OHT-G liver-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,liver 0.185
    lkp_muscle_ohtamg <- fixed(log(0.148))
    label("4-OHT-G muscle-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,mus 0.148
    lkp_skin_ohtamg <- fixed(log(0.507))
    label("4-OHT-G skin-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,skin 0.507
    lkp_adipose_ohtamg <- fixed(log(0.132))
    label("4-OHT-G adipose-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,adi 0.132
    lkp_gut_ohtamg <- fixed(log(0.316))
    label("4-OHT-G gut-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,gut 0.316
    fgut_ohtamg <- fixed(0.390)
    label("4-OHT-G enterocyte free fraction fgut (unitless)") # Suppl Table S1B: fgut 0.390
    fp_endox <- fixed(0.0200)
    label("Endoxifen unbound fraction in plasma (unitless)") # Suppl Table S1B: fp 0.0200
    rb_endox <- fixed(0.897)
    label("Endoxifen blood-to-plasma ratio RB (unitless)") # Suppl Table S1B: RB 0.897
    lkp_liver_endox <- fixed(log(12.3))
    label("Endoxifen liver-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,liver 12.3
    lkp_muscle_endox <- fixed(log(4.18))
    label("Endoxifen muscle-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,mus 4.18
    lkp_skin_endox <- fixed(log(4.02))
    label("Endoxifen skin-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,skin 4.02
    lkp_adipose_endox <- fixed(log(2.06))
    label("Endoxifen adipose-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,adi 2.06
    lkp_gut_endox <- fixed(log(6.74))
    label("Endoxifen gut-to-blood partition coefficient (unitless)") # Suppl Table S1B: Kp,gut 6.74
    fgut_endox <- fixed(0.0376)
    label("Endoxifen enterocyte free fraction fgut (unitless)") # Suppl Table S1B: fgut 0.0376

    # ---------------------------------------------------------------
    # 4. Between-subject variability -- Supplementary Table S3. Physiological
    #    volumes and flows (Table S3 "Physiological Parameters") are NORMAL,
    #    so their omega is CV^2 and the eta enters model() as (1 + eta).
    #    Kinetic / enzyme / transport parameters (Table S3 "Pharmacokinetic
    #    parameters") are LOGNORMAL, so omega is log(CV^2 + 1) and the eta
    #    enters as exp(eta). One random effect is drawn per PROCESS and
    #    shared across the five compounds, which matches the deposited
    #    parameter files' "same"-linkage of the kinetic constants: a subject
    #    with fast CYP3A metabolism is fast for every CYP3A step, so that
    #    shared draw is what makes the enzyme variability interpretable. The
    #    CYP2D6-pathway variability is NOT an eta here -- it is carried by
    #    the per-subject CYP2D6 activity covariate (Table 1 diplotype %CV).
    #    See the vignette Errata for the one eta (Vcentral) and the
    #    shared-draw simplification taken from the deposited parameter file.
    # ---------------------------------------------------------------
    etalv_central ~ 0.0861777 # deposited ID 1 parameter file 'Vb' CV 30 percent; omega = log(1 + 0.3^2); lognormal IIV on central volume (shared across compounds). Not in Suppl Table S3; see vignette Errata
    etalps_dif ~ 0.0099503 # Suppl Table S3 'PSdif,ent' CV 10 percent; omega = log(1 + 0.1^2); IIV on enterocyte passive diffusion
    etacl_3a4 ~ 0.1033685 # Suppl Table S3 'CLmet,hep by CYP3A' CV 33 percent; omega = log(1 + 0.33^2); IIV on hepatic CYP3A metabolism
    etacl_ugt ~ 0.0861777 # Suppl Table S3 'CLmet,hep by UGT' CV 30 percent; omega = log(1 + 0.3^2); IIV on UGT glucuronidation
    etalcl_other_hep ~ 0.0861777 # Suppl Table S3 'CLmet,hep by Other enzymes' CV 30 percent; omega = log(1 + 0.3^2); IIV on other-enzyme hepatic metabolism
    etalcl_bile ~ 0.0644423 # Suppl Table S3 'CLbile' CV 25.8 percent; omega = log(1 + 0.258^2); IIV on biliary intrinsic clearance
    etalkbile ~ 0.0861777 # Suppl Table S3 'kbile' CV 30 percent; omega = log(1 + 0.3^2); IIV on biliary transit rate constant
    etalka ~ 0.0861777 # Suppl Table S3 'ka' CV 30 percent; omega = log(1 + 0.3^2); IIV on absorption rate constant
    etalkfeces ~ 0.0861777 # Suppl Table S3 'kfeces' CV 30 percent; omega = log(1 + 0.3^2); IIV on faecal transit rate constant
    etaps_act_ent ~ 0.0644423 # Suppl Table S3 'PSact,ent' CV 25.8 percent; omega = log(1 + 0.258^2); IIV on endoxifen active intestinal efflux
    etacl_3a4_ent ~ 0.1033685 # Suppl Table S3 'CLmet,ent by CYP3A' CV 33 percent; omega = log(1 + 0.33^2); IIV on intestinal CYP3A metabolism
    etacl_other_ent ~ 0.0861777 # Suppl Table S3 'CLmet,ent by other enzymes' CV 30 percent; omega = log(1 + 0.3^2); IIV on other-enzyme intestinal metabolism
    etalcl_renal ~ 0.1105990 # Suppl Table S3 'CLint,sec' CV 34.2 percent; omega = log(1 + 0.342^2); IIV on renal clearance
    etalv_liver ~ 0.012996 # Suppl Table S3 'Vliver' CV 11.4 percent; omega = 0.114^2; normal IIV on liver volume
    etalv_muscle ~ 0.012996 # Suppl Table S3 'Vmus' CV 11.4 percent; omega = 0.114^2; normal IIV on muscle volume
    etalv_skin ~ 0.012996 # Suppl Table S3 'Vskin' CV 11.4 percent; omega = 0.114^2; normal IIV on skin volume
    etalv_adipose ~ 0.012996 # Suppl Table S3 'Vadi' CV 11.4 percent; omega = 0.114^2; normal IIV on adipose volume
    etalv_ent ~ 0.012996 # Suppl Table S3 footnote a (CV assumed same as liver) CV 11.4 percent; omega = 0.114^2; normal IIV on enterocyte volume
    etalv_mucb ~ 0.012996 # Suppl Table S3 'Vmuc,b' CV 11.4 percent; omega = 0.114^2; normal IIV on mucosal-blood volume
    etalv_serosa ~ 0.012996 # Suppl Table S3 'Vser' CV 11.4 percent; omega = 0.114^2; normal IIV on serosa volume
    etalq_liver ~ 0.014884 # Suppl Table S3 'Qhep' CV 12.2 percent; omega = 0.122^2; normal IIV on hepatic blood flow
    etalq_muscle ~ 0.014884 # Suppl Table S3 'Qmus' CV 12.2 percent; omega = 0.122^2; normal IIV on muscle blood flow
    etalq_skin ~ 0.014884 # Suppl Table S3 'Qskin' CV 12.2 percent; omega = 0.122^2; normal IIV on skin blood flow
    etalq_adipose ~ 0.014884 # Suppl Table S3 'Qadi' CV 12.2 percent; omega = 0.122^2; normal IIV on adipose blood flow
    etalq_mucosa ~ 0.014884 # Suppl Table S3 'Qmuc,b' CV 12.2 percent; omega = 0.122^2; normal IIV on mucosal blood flow
    etalq_serosa ~ 0.014884 # Suppl Table S3 'Qser' CV 12.2 percent; omega = 0.122^2; normal IIV on serosal blood flow

    # ---------------------------------------------------------------
    # 5. Residual error. The CNM minimises a weighted sum of squares (a
    #    fixed-effects criterion) and no residual-error model is reported, so
    #    the five propSd values are syntactic placeholders and must NOT be
    #    read as estimates. Same convention as Toshimoto_2017_irinotecan_pbpk.R.
    # ---------------------------------------------------------------
    propSd <- fixed(0.1)
    label("Tamoxifen proportional residual error placeholder (fraction)") # not reported in Nakamura 2018; placeholder only
    propSd_ndmtam <- fixed(0.1)
    label("NDMT proportional residual error placeholder (fraction)") # not reported in Nakamura 2018; placeholder only
    propSd_ohtam <- fixed(0.1)
    label("4-OHT proportional residual error placeholder (fraction)") # not reported in Nakamura 2018; placeholder only
    propSd_ohtamg <- fixed(0.1)
    label("4-OHT-G proportional residual error placeholder (fraction)") # not reported in Nakamura 2018; placeholder only
    propSd_endox <- fixed(0.1)
    label("Endoxifen proportional residual error placeholder (fraction)") # not reported in Nakamura 2018; placeholder only
  })

  model({
    # ===============================================================
    # 1. Physiological system, scaled to the individual body weight and
    #    perturbed by the NORMAL between-subject variability of
    #    Supplementary Table S3. Supplementary Tables S1 and S3 are
    #    per-kilogram, so every volume and flow is multiplied by WT.
    # ===============================================================
    v_liver <- exp(lv_liver) * WT * (1 + etalv_liver)
    v_muscle <- exp(lv_muscle) * WT * (1 + etalv_muscle)
    v_skin <- exp(lv_skin) * WT * (1 + etalv_skin)
    v_adipose <- exp(lv_adipose) * WT * (1 + etalv_adipose)
    v_ent <- exp(lv_ent) * WT * (1 + etalv_ent)
    v_mucb <- exp(lv_mucb) * WT * (1 + etalv_mucb)
    v_serosa <- exp(lv_serosa) * WT * (1 + etalv_serosa)
    q_liver <- exp(lq_liver) * WT * (1 + etalq_liver)
    q_muscle <- exp(lq_muscle) * WT * (1 + etalq_muscle)
    q_skin <- exp(lq_skin) * WT * (1 + etalq_skin)
    q_adipose <- exp(lq_adipose) * WT * (1 + etalq_adipose)
    q_mucosa <- exp(lq_mucosa) * WT * (1 + etalq_mucosa)
    q_serosa <- exp(lq_serosa) * WT * (1 + etalq_serosa)
    # Hepatic-artery flow: Qhep = Qha + Qser + Qmuc,b (Supplementary Methods
    # liver inflow term (Qhep - Qser - Qmuc) * Ccentral).
    q_art <- q_liver - q_serosa - q_mucosa
    # fdisp = 1/5: each of the five liver sub-compartments holds one fifth of
    # the liver volume and carries one fifth of every liver-wide clearance.
    fdisp <- 0.2
    v_liver_unit <- v_liver * fdisp

    # ===============================================================
    # 2. Compound-specific derived quantities. fb = fp / RB (blood unbound
    #    fraction); fliver = fb / (Kp,liver * SFKp) (the free-fraction factor
    #    multiplying every hepatic intrinsic clearance, Supplementary Methods
    #    "Optional equation"); the effective partition coefficients fold in
    #    the per-compound SFKp (every Kp appears as Kp * SFKp in the ODEs).
    # ===============================================================
    sfkp <- exp(lsfkp)
    fb <- fp / rb
    kpl <- exp(lkp_liver) * sfkp
    kpm <- exp(lkp_muscle) * sfkp
    kps <- exp(lkp_skin) * sfkp
    kpa <- exp(lkp_adipose) * sfkp
    kpg <- exp(lkp_gut) * sfkp
    fliver <- fb / kpl
    v_central <- exp(lv_central) * WT * exp(etalv_central)

    sfkp_ndmtam <- exp(lsfkp_ndmtam)
    fb_ndmtam <- fp_ndmtam / rb_ndmtam
    kpl_ndmtam <- exp(lkp_liver_ndmtam) * sfkp_ndmtam
    kpm_ndmtam <- exp(lkp_muscle_ndmtam) * sfkp_ndmtam
    kps_ndmtam <- exp(lkp_skin_ndmtam) * sfkp_ndmtam
    kpa_ndmtam <- exp(lkp_adipose_ndmtam) * sfkp_ndmtam
    kpg_ndmtam <- exp(lkp_gut_ndmtam) * sfkp_ndmtam
    fliver_ndmtam <- fb_ndmtam / kpl_ndmtam
    v_central_ndmtam <- exp(lv_central_ndmtam) * WT * exp(etalv_central)

    sfkp_ohtam <- exp(lsfkp_ohtam)
    fb_ohtam <- fp_ohtam / rb_ohtam
    kpl_ohtam <- exp(lkp_liver_ohtam) * sfkp_ohtam
    kpm_ohtam <- exp(lkp_muscle_ohtam) * sfkp_ohtam
    kps_ohtam <- exp(lkp_skin_ohtam) * sfkp_ohtam
    kpa_ohtam <- exp(lkp_adipose_ohtam) * sfkp_ohtam
    kpg_ohtam <- exp(lkp_gut_ohtam) * sfkp_ohtam
    fliver_ohtam <- fb_ohtam / kpl_ohtam
    v_central_ohtam <- exp(lv_central_ohtam) * WT * exp(etalv_central)

    sfkp_ohtamg <- exp(lsfkp_ohtamg)
    fb_ohtamg <- fp_ohtamg / rb_ohtamg
    kpl_ohtamg <- exp(lkp_liver_ohtamg) * sfkp_ohtamg
    kpm_ohtamg <- exp(lkp_muscle_ohtamg) * sfkp_ohtamg
    kps_ohtamg <- exp(lkp_skin_ohtamg) * sfkp_ohtamg
    kpa_ohtamg <- exp(lkp_adipose_ohtamg) * sfkp_ohtamg
    kpg_ohtamg <- exp(lkp_gut_ohtamg) * sfkp_ohtamg
    fliver_ohtamg <- fb_ohtamg / kpl_ohtamg
    v_central_ohtamg <- exp(lv_central_ohtamg) * WT * exp(etalv_central)

    sfkp_endox <- exp(lsfkp_endox)
    fb_endox <- fp_endox / rb_endox
    kpl_endox <- exp(lkp_liver_endox) * sfkp_endox
    kpm_endox <- exp(lkp_muscle_endox) * sfkp_endox
    kps_endox <- exp(lkp_skin_endox) * sfkp_endox
    kpa_endox <- exp(lkp_adipose_endox) * sfkp_endox
    kpg_endox <- exp(lkp_gut_endox) * sfkp_endox
    fliver_endox <- fb_endox / kpl_endox
    v_central_endox <- exp(lv_central_endox) * WT * exp(etalv_central)

    # ===============================================================
    # 3. Individual kinetic and metabolic parameters. Rate constants (ka,
    #    kbile, kfeces, kdeg) are 1/h and are NOT weight-scaled; every
    #    permeability-surface-area product and intrinsic clearance is per
    #    kilogram (L/h/kg) and is multiplied by WT. The CYP2D6 / CYP3A4 /
    #    CYP2C9 covariates (default 1) multiply the enzyme-specific pathways;
    #    the lognormal IIV etas of Supplementary Table S3 are applied once
    #    per process and shared across compounds.
    # ===============================================================
    ka <- exp(lka) * exp(etalka)
    kbile <- exp(lkbile) * exp(etalkbile)
    kfeces <- exp(lkfeces) * exp(etalkfeces)
    ps_dif <- exp(lps_dif) * WT * exp(etalps_dif)
    cl_renal <- exp(lcl_renal) * WT * exp(etalcl_renal)
    cl_bile <- exp(lcl_bile) * WT * exp(etalcl_bile)
    cl_tam_ndmt <- exp(lcl_ndmt_hep) * WT * CYP3A4 * exp(etacl_3a4)
    cl_tam_4oht_2d6 <- exp(lcl_4oht_2d6_hep) * WT * CYP2D6
    cl_tam_4oht_2c9 <- exp(lcl_4oht_2c9_hep) * WT * CYP2C9
    cl_tam_4oht_oth <- exp(lcl_4oht_other_hep) * WT * exp(etalcl_other_hep)
    cl_tam_other <- exp(lcl_other_hep) * WT * exp(etalcl_other_hep)
    cl_tam_ndmt_ent <- exp(lcl_ndmt_ent) * WT * CYP3A4 * exp(etacl_3a4_ent)
    cl_met_hep <- cl_tam_ndmt + cl_tam_4oht_2d6 + cl_tam_4oht_2c9 +
      cl_tam_4oht_oth + cl_tam_other
    cl_met_ent <- cl_tam_ndmt_ent

    ka_ndmtam <- exp(lka_ndmtam) * exp(etalka)
    kbile_ndmtam <- exp(lkbile_ndmtam) * exp(etalkbile)
    kfeces_ndmtam <- exp(lkfeces_ndmtam) * exp(etalkfeces)
    ps_dif_ndmtam <- exp(lps_dif_ndmtam) * WT * exp(etalps_dif)
    cl_bile_ndmtam <- exp(lcl_bile_ndmtam) * WT * exp(etalcl_bile)
    cl_ndmt_end <- exp(lcl_end_2d6_hep_ndmtam) * WT * CYP2D6
    cl_ndmt_other <- exp(lcl_other_3a4_hep_ndmtam) * WT * CYP3A4 * exp(etacl_3a4)
    cl_ndmt_other_ent <- exp(lcl_other_3a4_ent_ndmtam) * WT * CYP3A4 * exp(etacl_3a4_ent)
    cl_met_hep_ndmtam <- cl_ndmt_end + cl_ndmt_other
    cl_met_ent_ndmtam <- cl_ndmt_other_ent

    ka_ohtam <- exp(lka_ohtam) * exp(etalka)
    kbile_ohtam <- exp(lkbile_ohtam) * exp(etalkbile)
    kfeces_ohtam <- exp(lkfeces_ohtam) * exp(etalkfeces)
    ps_dif_ohtam <- exp(lps_dif_ohtam) * WT * exp(etalps_dif)
    cl_renal_ohtam <- exp(lcl_renal_ohtam) * WT * exp(etalcl_renal)
    cl_bile_ohtam <- exp(lcl_bile_ohtam) * WT * exp(etalcl_bile)
    cl_oht_end <- exp(lcl_end_3a4_hep_ohtam) * WT * CYP3A4 * exp(etacl_3a4)
    cl_oht_4ohtg <- exp(lcl_4ohtg_ugt_hep_ohtam) * WT * exp(etacl_ugt)
    cl_oht_other3a4 <- exp(lcl_other_3a4_hep_ohtam) * WT * CYP3A4 * exp(etacl_3a4)
    cl_oht_other <- exp(lcl_other_hep_ohtam) * WT * exp(etalcl_other_hep)
    cl_oht_end_ent <- exp(lcl_end_3a4_ent_ohtam) * WT * CYP3A4 * exp(etacl_3a4_ent)
    cl_oht_4ohtg_ent <- exp(lcl_4ohtg_ugt_ent_ohtam) * WT * exp(etacl_ugt)
    cl_oht_other_ent <- exp(lcl_other_3a4_ent_ohtam) * WT * CYP3A4 * exp(etacl_3a4_ent)
    cl_met_hep_ohtam <- cl_oht_end + cl_oht_4ohtg + cl_oht_other3a4 + cl_oht_other
    cl_met_ent_ohtam <- cl_oht_end_ent + cl_oht_4ohtg_ent + cl_oht_other_ent

    ka_ohtamg <- exp(lka_ohtamg) * exp(etalka)
    kbile_ohtamg <- exp(lkbile_ohtamg) * exp(etalkbile)
    ps_dif_ohtamg <- exp(lps_dif_ohtamg) * WT * exp(etalps_dif)
    cl_renal_ohtamg <- exp(lcl_renal_ohtamg) * WT * exp(etalcl_renal)
    cl_bile_ohtamg <- exp(lcl_bile_ohtamg) * WT * exp(etalcl_bile)
    kdeg <- exp(lkdeg_ohtamg)

    ka_endox <- exp(lka_endox) * exp(etalka)
    kbile_endox <- exp(lkbile_endox) * exp(etalkbile)
    kfeces_endox <- exp(lkfeces_endox) * exp(etalkfeces)
    ps_dif_endox <- exp(lps_dif_endox) * WT * exp(etalps_dif)
    ps_act_ent_endox <- exp(lps_act_ent_endox) * WT * exp(etaps_act_ent)
    cl_renal_endox <- exp(lcl_renal_endox) * WT * exp(etalcl_renal)
    cl_bile_endox <- exp(lcl_bile_endox) * WT * exp(etalcl_bile)
    cl_end_other <- exp(lcl_other_hep_endox) * WT * exp(etalcl_other_hep)
    cl_end_other_ent <- exp(lcl_other_ent_endox) * WT * exp(etacl_other_ent)
    cl_met_hep_endox <- cl_end_other
    cl_met_ent_endox <- cl_end_other_ent

    # ===============================================================
    # 4. ODE system -- Supplementary Methods, written out once per chemical
    #    species. The central, tissue and liver states hold CONCENTRATIONS
    #    (nM); the lumen, biliary, faecal and urinary states hold AMOUNTS
    #    (nmol). Two departures from the printed equations are made, both
    #    forced by mass balance and documented in the vignette Errata:
    #      (a) inter-liver transport uses the blood-equivalent concentration
    #          Qhep * Cliver / (Kp,liver * SFKp) at BOTH ends (the printed
    #          equation drops the partition divisor on the inflow); this is
    #          what the deposited MATLAB code does;
    #      (b) the metabolite formation source terms Xi(t) and Z(t) are
    #          printed without the parent's free-fraction factor and, in the
    #          liver, without the one-fifth dispersion weight that the
    #          matching loss term in the parent equation carries; both are
    #          restored so that what leaves the parent exactly arrives in the
    #          metabolite.
    # ===============================================================

    # ---------------- Tamoxifen (parent) ----------------
    d/dt(central) <- (q_liver * liver5 / kpl - q_liver * central +
      q_muscle * (muscle / kpm - central) +
      q_skin * (skin / kps - central) +
      q_adipose * (adipose / kpa - central) -
      cl_renal * central) / v_central
    d/dt(muscle) <- q_muscle * (central - muscle / kpm) / v_muscle
    d/dt(skin) <- q_skin * (central - skin / kps) / v_skin
    d/dt(adipose) <- q_adipose * (central - adipose / kpa) / v_adipose
    d/dt(serosa) <- q_serosa * (central - serosa / kpg) / v_serosa
    d/dt(liver1) <- (q_art * central + q_serosa * serosa / kpg +
      q_mucosa * intestine_muc - q_liver * liver1 / kpl -
      fliver * (cl_met_hep + cl_bile) * fdisp * liver1) / v_liver_unit
    d/dt(liver2) <- (q_liver * (liver1 - liver2) / kpl -
      fliver * (cl_met_hep + cl_bile) * fdisp * liver2) / v_liver_unit
    d/dt(liver3) <- (q_liver * (liver2 - liver3) / kpl -
      fliver * (cl_met_hep + cl_bile) * fdisp * liver3) / v_liver_unit
    d/dt(liver4) <- (q_liver * (liver3 - liver4) / kpl -
      fliver * (cl_met_hep + cl_bile) * fdisp * liver4) / v_liver_unit
    d/dt(liver5) <- (q_liver * (liver4 - liver5) / kpl -
      fliver * (cl_met_hep + cl_bile) * fdisp * liver5) / v_liver_unit
    d/dt(ehc1) <- fliver * cl_bile * fdisp *
      (liver1 + liver2 + liver3 + liver4 + liver5) - kbile * ehc1
    d/dt(ehc2) <- kbile * (ehc1 - ehc2)
    d/dt(ehc3) <- kbile * (ehc2 - ehc3)
    d/dt(gut_lumen) <- kbile * ehc3 +
      fgut * ar * ps_dif * intestine_ent -
      (ka + kfeces) * gut_lumen
    d/dt(intestine_ent) <- (ka * gut_lumen + fb * ps_dif * intestine_muc -
      fgut * ((ar + 1) * ps_dif + cl_met_ent) * intestine_ent) / v_ent
    d/dt(intestine_muc) <- (q_mucosa * central + fgut * ps_dif * intestine_ent -
      (q_mucosa + fb * ps_dif) * intestine_muc) / v_mucb
    d/dt(a_feces) <- kfeces * gut_lumen
    d/dt(a_urine) <- cl_renal * central

    # ---------------- NDMT ----------------
    d/dt(central_ndmtam) <- (q_liver * liver5_ndmtam / kpl_ndmtam - q_liver * central_ndmtam +
      q_muscle * (muscle_ndmtam / kpm_ndmtam - central_ndmtam) +
      q_skin * (skin_ndmtam / kps_ndmtam - central_ndmtam) +
      q_adipose * (adipose_ndmtam / kpa_ndmtam - central_ndmtam)) / v_central_ndmtam
    d/dt(muscle_ndmtam) <- q_muscle * (central_ndmtam - muscle_ndmtam / kpm_ndmtam) / v_muscle
    d/dt(skin_ndmtam) <- q_skin * (central_ndmtam - skin_ndmtam / kps_ndmtam) / v_skin
    d/dt(adipose_ndmtam) <- q_adipose * (central_ndmtam - adipose_ndmtam / kpa_ndmtam) / v_adipose
    d/dt(serosa_ndmtam) <- q_serosa * (central_ndmtam - serosa_ndmtam / kpg_ndmtam) / v_serosa
    d/dt(liver1_ndmtam) <- (q_art * central_ndmtam + q_serosa * serosa_ndmtam / kpg_ndmtam +
      q_mucosa * intestine_muc_ndmtam - q_liver * liver1_ndmtam / kpl_ndmtam -
      fliver_ndmtam * (cl_met_hep_ndmtam + cl_bile_ndmtam) * fdisp * liver1_ndmtam +
      fliver * cl_tam_ndmt * fdisp * liver1) / v_liver_unit
    d/dt(liver2_ndmtam) <- (q_liver * (liver1_ndmtam - liver2_ndmtam) / kpl_ndmtam -
      fliver_ndmtam * (cl_met_hep_ndmtam + cl_bile_ndmtam) * fdisp * liver2_ndmtam +
      fliver * cl_tam_ndmt * fdisp * liver2) / v_liver_unit
    d/dt(liver3_ndmtam) <- (q_liver * (liver2_ndmtam - liver3_ndmtam) / kpl_ndmtam -
      fliver_ndmtam * (cl_met_hep_ndmtam + cl_bile_ndmtam) * fdisp * liver3_ndmtam +
      fliver * cl_tam_ndmt * fdisp * liver3) / v_liver_unit
    d/dt(liver4_ndmtam) <- (q_liver * (liver3_ndmtam - liver4_ndmtam) / kpl_ndmtam -
      fliver_ndmtam * (cl_met_hep_ndmtam + cl_bile_ndmtam) * fdisp * liver4_ndmtam +
      fliver * cl_tam_ndmt * fdisp * liver4) / v_liver_unit
    d/dt(liver5_ndmtam) <- (q_liver * (liver4_ndmtam - liver5_ndmtam) / kpl_ndmtam -
      fliver_ndmtam * (cl_met_hep_ndmtam + cl_bile_ndmtam) * fdisp * liver5_ndmtam +
      fliver * cl_tam_ndmt * fdisp * liver5) / v_liver_unit
    d/dt(ehc1_ndmtam) <- fliver_ndmtam * cl_bile_ndmtam * fdisp *
      (liver1_ndmtam + liver2_ndmtam + liver3_ndmtam + liver4_ndmtam + liver5_ndmtam) -
      kbile_ndmtam * ehc1_ndmtam
    d/dt(ehc2_ndmtam) <- kbile_ndmtam * (ehc1_ndmtam - ehc2_ndmtam)
    d/dt(ehc3_ndmtam) <- kbile_ndmtam * (ehc2_ndmtam - ehc3_ndmtam)
    d/dt(gut_lumen_ndmtam) <- kbile_ndmtam * ehc3_ndmtam +
      fgut_ndmtam * ar * ps_dif_ndmtam * intestine_ent_ndmtam -
      (ka_ndmtam + kfeces_ndmtam) * gut_lumen_ndmtam
    d/dt(intestine_ent_ndmtam) <- (ka_ndmtam * gut_lumen_ndmtam +
      fb_ndmtam * ps_dif_ndmtam * intestine_muc_ndmtam -
      fgut_ndmtam * ((ar + 1) * ps_dif_ndmtam + cl_met_ent_ndmtam) * intestine_ent_ndmtam +
      fgut * cl_tam_ndmt_ent * intestine_ent) / v_ent
    d/dt(intestine_muc_ndmtam) <- (q_mucosa * central_ndmtam +
      fgut_ndmtam * ps_dif_ndmtam * intestine_ent_ndmtam -
      (q_mucosa + fb_ndmtam * ps_dif_ndmtam) * intestine_muc_ndmtam) / v_mucb
    d/dt(a_feces_ndmtam) <- kfeces_ndmtam * gut_lumen_ndmtam

    # ---------------- 4-OHT ----------------
    d/dt(central_ohtam) <- (q_liver * liver5_ohtam / kpl_ohtam - q_liver * central_ohtam +
      q_muscle * (muscle_ohtam / kpm_ohtam - central_ohtam) +
      q_skin * (skin_ohtam / kps_ohtam - central_ohtam) +
      q_adipose * (adipose_ohtam / kpa_ohtam - central_ohtam) -
      cl_renal_ohtam * central_ohtam) / v_central_ohtam
    d/dt(muscle_ohtam) <- q_muscle * (central_ohtam - muscle_ohtam / kpm_ohtam) / v_muscle
    d/dt(skin_ohtam) <- q_skin * (central_ohtam - skin_ohtam / kps_ohtam) / v_skin
    d/dt(adipose_ohtam) <- q_adipose * (central_ohtam - adipose_ohtam / kpa_ohtam) / v_adipose
    d/dt(serosa_ohtam) <- q_serosa * (central_ohtam - serosa_ohtam / kpg_ohtam) / v_serosa
    d/dt(liver1_ohtam) <- (q_art * central_ohtam + q_serosa * serosa_ohtam / kpg_ohtam +
      q_mucosa * intestine_muc_ohtam - q_liver * liver1_ohtam / kpl_ohtam -
      fliver_ohtam * (cl_met_hep_ohtam + cl_bile_ohtam) * fdisp * liver1_ohtam +
      fliver * (cl_tam_4oht_2d6 + cl_tam_4oht_2c9 + cl_tam_4oht_oth) * fdisp * liver1) / v_liver_unit
    d/dt(liver2_ohtam) <- (q_liver * (liver1_ohtam - liver2_ohtam) / kpl_ohtam -
      fliver_ohtam * (cl_met_hep_ohtam + cl_bile_ohtam) * fdisp * liver2_ohtam +
      fliver * (cl_tam_4oht_2d6 + cl_tam_4oht_2c9 + cl_tam_4oht_oth) * fdisp * liver2) / v_liver_unit
    d/dt(liver3_ohtam) <- (q_liver * (liver2_ohtam - liver3_ohtam) / kpl_ohtam -
      fliver_ohtam * (cl_met_hep_ohtam + cl_bile_ohtam) * fdisp * liver3_ohtam +
      fliver * (cl_tam_4oht_2d6 + cl_tam_4oht_2c9 + cl_tam_4oht_oth) * fdisp * liver3) / v_liver_unit
    d/dt(liver4_ohtam) <- (q_liver * (liver3_ohtam - liver4_ohtam) / kpl_ohtam -
      fliver_ohtam * (cl_met_hep_ohtam + cl_bile_ohtam) * fdisp * liver4_ohtam +
      fliver * (cl_tam_4oht_2d6 + cl_tam_4oht_2c9 + cl_tam_4oht_oth) * fdisp * liver4) / v_liver_unit
    d/dt(liver5_ohtam) <- (q_liver * (liver4_ohtam - liver5_ohtam) / kpl_ohtam -
      fliver_ohtam * (cl_met_hep_ohtam + cl_bile_ohtam) * fdisp * liver5_ohtam +
      fliver * (cl_tam_4oht_2d6 + cl_tam_4oht_2c9 + cl_tam_4oht_oth) * fdisp * liver5) / v_liver_unit
    d/dt(ehc1_ohtam) <- fliver_ohtam * cl_bile_ohtam * fdisp *
      (liver1_ohtam + liver2_ohtam + liver3_ohtam + liver4_ohtam + liver5_ohtam) -
      kbile_ohtam * ehc1_ohtam
    d/dt(ehc2_ohtam) <- kbile_ohtam * (ehc1_ohtam - ehc2_ohtam)
    d/dt(ehc3_ohtam) <- kbile_ohtam * (ehc2_ohtam - ehc3_ohtam)
    # 4-OHT additionally gains the luminal deconjugation of 4-OHT-G.
    d/dt(gut_lumen_ohtam) <- kbile_ohtam * ehc3_ohtam +
      fgut_ohtam * ar * ps_dif_ohtam * intestine_ent_ohtam -
      (ka_ohtam + kfeces_ohtam) * gut_lumen_ohtam +
      kdeg * gut_lumen_ohtamg
    d/dt(intestine_ent_ohtam) <- (ka_ohtam * gut_lumen_ohtam +
      fb_ohtam * ps_dif_ohtam * intestine_muc_ohtam -
      fgut_ohtam * ((ar + 1) * ps_dif_ohtam + cl_met_ent_ohtam) * intestine_ent_ohtam) / v_ent
    d/dt(intestine_muc_ohtam) <- (q_mucosa * central_ohtam +
      fgut_ohtam * ps_dif_ohtam * intestine_ent_ohtam -
      (q_mucosa + fb_ohtam * ps_dif_ohtam) * intestine_muc_ohtam) / v_mucb
    d/dt(a_feces_ohtam) <- kfeces_ohtam * gut_lumen_ohtam
    d/dt(a_urine_ohtam) <- cl_renal_ohtam * central_ohtam

    # ---------------- 4-OHT-glucuronide ----------------
    d/dt(central_ohtamg) <- (q_liver * liver5_ohtamg / kpl_ohtamg - q_liver * central_ohtamg +
      q_muscle * (muscle_ohtamg / kpm_ohtamg - central_ohtamg) +
      q_skin * (skin_ohtamg / kps_ohtamg - central_ohtamg) +
      q_adipose * (adipose_ohtamg / kpa_ohtamg - central_ohtamg) -
      cl_renal_ohtamg * central_ohtamg) / v_central_ohtamg
    d/dt(muscle_ohtamg) <- q_muscle * (central_ohtamg - muscle_ohtamg / kpm_ohtamg) / v_muscle
    d/dt(skin_ohtamg) <- q_skin * (central_ohtamg - skin_ohtamg / kps_ohtamg) / v_skin
    d/dt(adipose_ohtamg) <- q_adipose * (central_ohtamg - adipose_ohtamg / kpa_ohtamg) / v_adipose
    d/dt(serosa_ohtamg) <- q_serosa * (central_ohtamg - serosa_ohtamg / kpg_ohtamg) / v_serosa
    # 4-OHT-G has no metabolism (CLmet,hep = 0); it is formed from 4-OHT by
    # UGT glucuronidation in the liver.
    d/dt(liver1_ohtamg) <- (q_art * central_ohtamg + q_serosa * serosa_ohtamg / kpg_ohtamg +
      q_mucosa * intestine_muc_ohtamg - q_liver * liver1_ohtamg / kpl_ohtamg -
      fliver_ohtamg * cl_bile_ohtamg * fdisp * liver1_ohtamg +
      fliver_ohtam * cl_oht_4ohtg * fdisp * liver1_ohtam) / v_liver_unit
    d/dt(liver2_ohtamg) <- (q_liver * (liver1_ohtamg - liver2_ohtamg) / kpl_ohtamg -
      fliver_ohtamg * cl_bile_ohtamg * fdisp * liver2_ohtamg +
      fliver_ohtam * cl_oht_4ohtg * fdisp * liver2_ohtam) / v_liver_unit
    d/dt(liver3_ohtamg) <- (q_liver * (liver2_ohtamg - liver3_ohtamg) / kpl_ohtamg -
      fliver_ohtamg * cl_bile_ohtamg * fdisp * liver3_ohtamg +
      fliver_ohtam * cl_oht_4ohtg * fdisp * liver3_ohtam) / v_liver_unit
    d/dt(liver4_ohtamg) <- (q_liver * (liver3_ohtamg - liver4_ohtamg) / kpl_ohtamg -
      fliver_ohtamg * cl_bile_ohtamg * fdisp * liver4_ohtamg +
      fliver_ohtam * cl_oht_4ohtg * fdisp * liver4_ohtam) / v_liver_unit
    d/dt(liver5_ohtamg) <- (q_liver * (liver4_ohtamg - liver5_ohtamg) / kpl_ohtamg -
      fliver_ohtamg * cl_bile_ohtamg * fdisp * liver5_ohtamg +
      fliver_ohtam * cl_oht_4ohtg * fdisp * liver5_ohtam) / v_liver_unit
    d/dt(ehc1_ohtamg) <- fliver_ohtamg * cl_bile_ohtamg * fdisp *
      (liver1_ohtamg + liver2_ohtamg + liver3_ohtamg + liver4_ohtamg + liver5_ohtamg) -
      kbile_ohtamg * ehc1_ohtamg
    d/dt(ehc2_ohtamg) <- kbile_ohtamg * (ehc1_ohtamg - ehc2_ohtamg)
    d/dt(ehc3_ohtamg) <- kbile_ohtamg * (ehc2_ohtamg - ehc3_ohtamg)
    # 4-OHT-G has kfeces fixed to zero; it loses the luminal deconjugation
    # flux to 4-OHT (gained by gut_lumen_ohtam above).
    d/dt(gut_lumen_ohtamg) <- kbile_ohtamg * ehc3_ohtamg +
      fgut_ohtamg * ar * ps_dif_ohtamg * intestine_ent_ohtamg -
      ka_ohtamg * gut_lumen_ohtamg -
      kdeg * gut_lumen_ohtamg
    d/dt(intestine_ent_ohtamg) <- (ka_ohtamg * gut_lumen_ohtamg +
      fb_ohtamg * ps_dif_ohtamg * intestine_muc_ohtamg -
      fgut_ohtamg * ((ar + 1) * ps_dif_ohtamg) * intestine_ent_ohtamg +
      fgut_ohtam * cl_oht_4ohtg_ent * intestine_ent_ohtam) / v_ent
    d/dt(intestine_muc_ohtamg) <- (q_mucosa * central_ohtamg +
      fgut_ohtamg * ps_dif_ohtamg * intestine_ent_ohtamg -
      (q_mucosa + fb_ohtamg * ps_dif_ohtamg) * intestine_muc_ohtamg) / v_mucb
    d/dt(a_urine_ohtamg) <- cl_renal_ohtamg * central_ohtamg

    # ---------------- Endoxifen ----------------
    d/dt(central_endox) <- (q_liver * liver5_endox / kpl_endox - q_liver * central_endox +
      q_muscle * (muscle_endox / kpm_endox - central_endox) +
      q_skin * (skin_endox / kps_endox - central_endox) +
      q_adipose * (adipose_endox / kpa_endox - central_endox) -
      cl_renal_endox * central_endox) / v_central_endox
    d/dt(muscle_endox) <- q_muscle * (central_endox - muscle_endox / kpm_endox) / v_muscle
    d/dt(skin_endox) <- q_skin * (central_endox - skin_endox / kps_endox) / v_skin
    d/dt(adipose_endox) <- q_adipose * (central_endox - adipose_endox / kpa_endox) / v_adipose
    d/dt(serosa_endox) <- q_serosa * (central_endox - serosa_endox / kpg_endox) / v_serosa
    # Endoxifen is formed in the liver from NDMT (CYP2D6) and from 4-OHT (CYP3A).
    d/dt(liver1_endox) <- (q_art * central_endox + q_serosa * serosa_endox / kpg_endox +
      q_mucosa * intestine_muc_endox - q_liver * liver1_endox / kpl_endox -
      fliver_endox * (cl_met_hep_endox + cl_bile_endox) * fdisp * liver1_endox +
      fliver_ndmtam * cl_ndmt_end * fdisp * liver1_ndmtam +
      fliver_ohtam * cl_oht_end * fdisp * liver1_ohtam) / v_liver_unit
    d/dt(liver2_endox) <- (q_liver * (liver1_endox - liver2_endox) / kpl_endox -
      fliver_endox * (cl_met_hep_endox + cl_bile_endox) * fdisp * liver2_endox +
      fliver_ndmtam * cl_ndmt_end * fdisp * liver2_ndmtam +
      fliver_ohtam * cl_oht_end * fdisp * liver2_ohtam) / v_liver_unit
    d/dt(liver3_endox) <- (q_liver * (liver2_endox - liver3_endox) / kpl_endox -
      fliver_endox * (cl_met_hep_endox + cl_bile_endox) * fdisp * liver3_endox +
      fliver_ndmtam * cl_ndmt_end * fdisp * liver3_ndmtam +
      fliver_ohtam * cl_oht_end * fdisp * liver3_ohtam) / v_liver_unit
    d/dt(liver4_endox) <- (q_liver * (liver3_endox - liver4_endox) / kpl_endox -
      fliver_endox * (cl_met_hep_endox + cl_bile_endox) * fdisp * liver4_endox +
      fliver_ndmtam * cl_ndmt_end * fdisp * liver4_ndmtam +
      fliver_ohtam * cl_oht_end * fdisp * liver4_ohtam) / v_liver_unit
    d/dt(liver5_endox) <- (q_liver * (liver4_endox - liver5_endox) / kpl_endox -
      fliver_endox * (cl_met_hep_endox + cl_bile_endox) * fdisp * liver5_endox +
      fliver_ndmtam * cl_ndmt_end * fdisp * liver5_ndmtam +
      fliver_ohtam * cl_oht_end * fdisp * liver5_ohtam) / v_liver_unit
    d/dt(ehc1_endox) <- fliver_endox * cl_bile_endox * fdisp *
      (liver1_endox + liver2_endox + liver3_endox + liver4_endox + liver5_endox) -
      kbile_endox * ehc1_endox
    d/dt(ehc2_endox) <- kbile_endox * (ehc1_endox - ehc2_endox)
    d/dt(ehc3_endox) <- kbile_endox * (ehc2_endox - ehc3_endox)
    # Endoxifen has active (Pgp) apical efflux PSact,ent in addition to
    # passive diffusion.
    d/dt(gut_lumen_endox) <- kbile_endox * ehc3_endox +
      fgut_endox * (ar * ps_dif_endox + ps_act_ent_endox) * intestine_ent_endox -
      (ka_endox + kfeces_endox) * gut_lumen_endox
    d/dt(intestine_ent_endox) <- (ka_endox * gut_lumen_endox +
      fb_endox * ps_dif_endox * intestine_muc_endox -
      fgut_endox * ((ar + 1) * ps_dif_endox + cl_met_ent_endox + ps_act_ent_endox) *
        intestine_ent_endox +
      fgut_ohtam * cl_oht_end_ent * intestine_ent_ohtam) / v_ent
    d/dt(intestine_muc_endox) <- (q_mucosa * central_endox +
      fgut_endox * ps_dif_endox * intestine_ent_endox -
      (q_mucosa + fb_endox * ps_dif_endox) * intestine_muc_endox) / v_mucb
    d/dt(a_feces_endox) <- kfeces_endox * gut_lumen_endox
    d/dt(a_urine_endox) <- cl_renal_endox * central_endox

    # ===============================================================
    # 5. Observations. Each central state already is a blood concentration
    #    in nM. Tamoxifen is dosed orally into gut_lumen in the event table
    #    (no bioavailability transform: the absorbed fraction emerges from
    #    the ka / kfeces competition).
    # ===============================================================
    Cc <- central
    Cc_ndmtam <- central_ndmtam
    Cc_ohtam <- central_ohtam
    Cc_ohtamg <- central_ohtamg
    Cc_endox <- central_endox
    Cc ~ prop(propSd)
    Cc_ndmtam ~ prop(propSd_ndmtam)
    Cc_ohtam ~ prop(propSd_ohtam)
    Cc_ohtamg ~ prop(propSd_ohtamg)
    Cc_endox ~ prop(propSd_endox)
  })
}
