Yang_2016_methylphenidate_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited, acslX) with an ACAT-style physiological",
    "absorption model for oral immediate-release (IR) and extended-release (ER)",
    "methylphenidate (MPH) in adults. d- and l-MPH are carried separately; each",
    "has a nine-segment gut lumen (stomach, duodenum, two jejunum and three",
    "ileum segments, cecum, ascending colon) holding ER beads, undissolved IR",
    "drug and dissolved drug, an enterocyte state per intestinal segment",
    "(pH-dependent passive permeability, first-order non-specific gut loss),",
    "and seven perfused tissues plus plasma. Hepatic CES1A1 hydrolysis is",
    "Michaelis-Menten with competitive d/l inhibition; hepatic oxidation is",
    "linear. IR drug dissolves by a Noyes-Whitney diffusion-layer model with",
    "pH- and bile-dependent solubility; ER beads release through a lagged",
    "Weibull term (Tlag 4 h). Body weight scales flows and Michaelis-Menten",
    "Vmax as BW^0.75 and volumes as BW; sex selects the male or female",
    "fractional flows and volumes; fed state switches gastric emptying and",
    "proximal-gut pH and bile. Ritalinic acid formed in liver and gut is",
    "carried as a one-compartment metabolite. Deterministic typical-value",
    "model (the paper's Monte Carlo distributions are not encoded)."
  )
  reference <- paste(
    "Yang X, Duan J, Fisher J (2016). Application of Physiologically Based",
    "Absorption Modeling to Characterize the Pharmacokinetic Profiles of",
    "Oral Extended Release Methylphenidate Products in Adults. PLoS One",
    "11(10):e0164641. doi:10.1371/journal.pone.0164641.",
    "Structure and fixed constants from the acslX source code in S1 File",
    "(PBPK modeling code); calibrated and visually fitted values from",
    "Table 2; gut physiology from Table 1; Equations 1-6."
  )
  vignette <- "Yang_2016_methylphenidate_pbpk"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # IR drug is dosed into the stomach_ir_* states and ER beads into the
  # stomach_er_* states, half of the racemic dose into each enantiomer.
  dosing <- c("stomach_ir_d", "stomach_ir_l", "stomach_er_d", "stomach_er_l")

  paper_specific_compartments <- c(
    "stomach_er_d",
    "stomach_ir_d",
    "stomach_lumen_d",
    "duodenum_er_d",
    "duodenum_ir_d",
    "duodenum_lumen_d",
    "duodenum_ent_d",
    "jejunum1_er_d",
    "jejunum1_ir_d",
    "jejunum1_lumen_d",
    "jejunum1_ent_d",
    "jejunum2_er_d",
    "jejunum2_ir_d",
    "jejunum2_lumen_d",
    "jejunum2_ent_d",
    "ileum1_er_d",
    "ileum1_ir_d",
    "ileum1_lumen_d",
    "ileum1_ent_d",
    "ileum2_er_d",
    "ileum2_ir_d",
    "ileum2_lumen_d",
    "ileum2_ent_d",
    "ileum3_er_d",
    "ileum3_ir_d",
    "ileum3_lumen_d",
    "ileum3_ent_d",
    "cecum_er_d",
    "cecum_ir_d",
    "cecum_lumen_d",
    "cecum_ent_d",
    "colon_er_d",
    "colon_ir_d",
    "colon_lumen_d",
    "colon_ent_d",
    "central_ra_d",
    "stomach_er_l",
    "stomach_ir_l",
    "stomach_lumen_l",
    "duodenum_er_l",
    "duodenum_ir_l",
    "duodenum_lumen_l",
    "duodenum_ent_l",
    "jejunum1_er_l",
    "jejunum1_ir_l",
    "jejunum1_lumen_l",
    "jejunum1_ent_l",
    "jejunum2_er_l",
    "jejunum2_ir_l",
    "jejunum2_lumen_l",
    "jejunum2_ent_l",
    "ileum1_er_l",
    "ileum1_ir_l",
    "ileum1_lumen_l",
    "ileum1_ent_l",
    "ileum2_er_l",
    "ileum2_ir_l",
    "ileum2_lumen_l",
    "ileum2_ent_l",
    "ileum3_er_l",
    "ileum3_ir_l",
    "ileum3_lumen_l",
    "ileum3_ent_l",
    "cecum_er_l",
    "cecum_ir_l",
    "cecum_lumen_l",
    "cecum_ent_l",
    "colon_er_l",
    "colon_ir_l",
    "colon_lumen_l",
    "colon_ent_l",
    "central_ra_l"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cardiac output, enterocyte loss constants, hepatic Vmax and oxidation",
        "clearance scale with WT^0.75; tissue, plasma and enterocyte volumes",
        "scale with WT. The acslX source default is 80 kg; the paper set the",
        "body-weight distribution per study (S1 Table: 83.4 +/- 25.4 kg for",
        "the default population)."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; the acslX source defaults are the male values)",
      notes = paste(
        "Selects the female fractional blood flows to fat, gonads, heart,",
        "liver and each intestinal segment, and the female fractional volumes",
        "of fat, gonads, heart and the duodenal, ileal and ascending-colon",
        "walls (Table 1 footnote h and the female values in the S1 File",
        "comments / S1 Table)."
      ),
      source_name = "male/female"
    ),
    FED = list(
      description = "Fed state at dosing, 1 = fed, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Switches the gastric transit time (0.25 -> 1 h), the pH of the",
        "stomach, duodenum and both jejunum segments, and the bile-salt",
        "concentration of the duodenum to ileum 3 to the fed values of",
        "Table 1 footnote f. 'Fed' in the paper means a meal in general,",
        "not a high-fat meal. Only the Medikinet Retard studies were fed."
      ),
      source_name = "fasting/fed"
    ),
    FORM_METHYLPHENIDATE_MLR = list(
      description = "Multilayer-release methylphenidate (MLR-MPH) ER formulation indicator, 1 = MLR-MPH",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Ritalin LA, Metadate CD or Medikinet Retard ER beads)",
      notes = paste(
        "Selects the MLR-MPH Weibull release constants (A = 6 h^b, b = 4)",
        "instead of the shared Ritalin LA / Metadate CD / Medikinet Retard",
        "constants (A = 0.5 h^b, b = 2) (Table 2). Has no effect on the IR",
        "component or on an IR-only product. The IR:ER split is a dosing",
        "choice, not a covariate: MLR-MPH 40:60 (37:63 for the 80 mg",
        "capsule), Ritalin LA 50:50, Metadate CD 30:70, Medikinet Retard",
        "50:50."
      ),
      source_name = "MLR-MPH"
    )
  )

  compartmentData <- list(
    stomach_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    stomach_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    stomach_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_er_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_ir_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum1_er_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_ir_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum2_er_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_ir_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum1_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum1_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum1_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum1_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum2_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum2_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum2_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum2_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum3_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum3_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum3_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum3_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    cecum_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    cecum_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    cecum_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    cecum_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    colon_er_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    colon_ir_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    colon_lumen_d = list(
      analyte = "d-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    colon_ent_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    central_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "plasma", verified = TRUE),
    adipose_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    gonads_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    liver_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    brain_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    heart_d = list(analyte = "d-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    central_ra_d = list(analyte = "d-ritalinic acid", units = "ug", specimen = "plasma", verified = TRUE),
    stomach_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    stomach_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    stomach_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_er_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_ir_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    duodenum_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum1_er_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_ir_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum1_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum2_er_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_ir_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum2_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum1_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum1_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum1_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum1_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum2_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum2_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum2_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum2_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    ileum3_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum3_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    ileum3_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    ileum3_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    cecum_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    cecum_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    cecum_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    cecum_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    colon_er_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    colon_ir_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "administration site", verified = TRUE),
    colon_lumen_l = list(
      analyte = "l-methylphenidate",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    colon_ent_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    central_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "plasma", verified = TRUE),
    adipose_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    gonads_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    liver_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    brain_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    heart_l = list(analyte = "l-methylphenidate", units = "ug", specimen = "tissue", verified = TRUE),
    central_ra_l = list(analyte = "l-ritalinic acid", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 266,
    n_studies = 13,
    age_range = "18-45 years (study means 24-35 years)",
    weight_range = "study means 59.6-82.2 kg",
    sex_female_pct = NA,
    disease_state = "healthy adult volunteers",
    dose_range = paste(
      "IR MPH 0.3 mg/kg or 40 mg; MLR-MPH 20 or 80 mg; Ritalin LA 20 or",
      "40 mg; Metadate CD 20 mg; Medikinet Retard 20 or 40 mg; single oral doses"
    ),
    regions = "not reported",
    notes = paste(
      "Not a fit to individual data: the model was calibrated to digitized",
      "mean d- and l-MPH profiles from three IR studies (Patrick 2013, n = 24;",
      "Patrick 2007, n = 20; Wong 1998, n = 21) and evaluated against mean",
      "total-MPH profiles from ten ER datasets (Methods, 'Key pharmacokinetic",
      "studies'). n_subjects sums the 13 datasets, counting once the 28 men",
      "who received both Ritalin LA and Medikinet Retard (Haessler 2008).",
      "Mixed-sex cohorts except Wong 1998 and Haessler 2008 (men only)."
    )
  )

  ini({
    # All values are fixed: Yang 2016 is a PBPK model run by simulation, and
    # the three optimized constants (k5_dc, k5_lc, hpeff) and the visually
    # fitted Weibull constants were set by the authors, not estimated with
    # uncertainty. 'S1 File' is the acslX source code.

    # ---- Blood flows (fraction of cardiac output unless stated) ----
    qcc <- fixed(15.87); label("Cardiac output (L/h/kg^0.75)") # S1 File QCC = 15.87; S1 Table QCC 15.87
    fq_fat <- fixed(0.053); label("Fat blood flow, male (fraction of QC)") # S1 File QFatC = 0.053
    fq_fat_f <- fixed(0.091); label("Fat blood flow, female (fraction of QC)") # S1 File QFatC comment 'female: 0.091'; S1 Table
    fq_brain <- fixed(0.11); label("Brain blood flow (fraction of QC)") # S1 File QbrainC = 0.11
    fq_gonad <- fixed(0.00054); label("Gonad blood flow, male (fraction of QC)") # S1 File QgonadC = 0.00054
    fq_gonad_f <- fixed(0.00022); label("Gonad blood flow, female (fraction of QC)") # S1 File QgonadC comment 'female: 0.00022'
    fq_heart <- fixed(0.038); label("Heart blood flow, male (fraction of QC)") # S1 File QHeartC = 0.038
    fq_heart_f <- fixed(0.047); label("Heart blood flow, female (fraction of QC)") # S1 File QHeartC comment 'female: 0.047'
    fq_liver <- fixed(0.255); label("Total liver blood flow (portal + arterial), male (fraction of QC)") # S1 File QLiverC = 0.255; Methods 'Systemic uptake'
    fq_liver_f <- fixed(0.27); label("Total liver blood flow (portal + arterial), female (fraction of QC)") # S1 File QLiverC comment '0.27 for female'
    fq_seg <- fixed(0.02375); label("Blood flow to each of the 8 intestinal segments, male (fraction of QC)") # S1 File QDUOC..QASCENDINGC = 0.02375 ((0.255 - 0.065)/8); Table 1 0.024
    fq_seg_f <- fixed(0.025625); label("Blood flow to each of the 8 intestinal segments, female (fraction of QC)") # S1 File comment '0.025625 for females'; Table 1 0.026
    fq_slow_tot <- fixed(0.24); label("Slowly perfused + fat + heart blood flow (fraction of QC)") # S1 File 'QS = 0.24*QC - QFat - Qheart'
    fq_rich_tot <- fixed(0.76); label("Richly perfused + liver + gonad + brain blood flow (fraction of QC)") # S1 File 'QR = 0.76*QC - QLiver - Qgonad - Qbrain'

    # ---- Tissue volumes (fraction of body weight, L/kg) ----
    fv_plasma <- fixed(0.0435); label("Plasma volume (fraction of BW)") # S1 File VplasmaC = 0.0435
    fv_fat <- fixed(0.213); label("Fat volume, male (fraction of BW)") # S1 File VFatC = 0.213
    fv_fat_f <- fixed(0.327); label("Fat volume, female (fraction of BW)") # S1 File VFatC comment 'female: 0.327'
    fv_brain <- fixed(0.02); label("Brain volume (fraction of BW)") # S1 File VbrainC = 0.02
    fv_gonad <- fixed(0.0007); label("Gonad volume, male (fraction of BW)") # S1 File VgonadC = 0.0007
    fv_gonad_f <- fixed(0.0027); label("Gonad volume, female (fraction of BW)") # S1 File VgonadC comment 'female: 0.0027'
    fv_heart <- fixed(0.0045); label("Heart volume, male (fraction of BW)") # S1 File VHeartC = 0.0045
    fv_heart_f <- fixed(0.0042); label("Heart volume, female (fraction of BW)") # S1 File VHeartC comment 'female: 0.0042'
    fv_liver <- fixed(0.026); label("Liver volume (fraction of BW)") # S1 File VLiverC = 0.026
    fv_slow_tot <- fixed(0.60); label("Slowly perfused + fat + heart volume (fraction of BW)") # S1 File 'VS = 0.60*BW - VFat - Vheart'
    fv_rich_tot <- fixed(0.33); label("Richly perfused + liver + plasma + gonad + brain + gut-wall volume (fraction of BW)") # S1 File 'VR = 0.33*BW - VLiver - Vplasma - ...'
    fv_body <- fixed(0.6); label("Ritalinic acid distribution volume (fraction of BW)") # S1 File VbodyC = 0.6
    fv_ent_duodenum <- fixed(0.000767); label("Duodenal wall volume, male (fraction of BW)") # S1 File VMEMDUOC = 0.000767; Table 1 0.00077
    fv_ent_duodenum_f <- fixed(0.000877); label("Duodenal wall volume, female (fraction of BW)") # S1 File comment 'Female: 0.000877'; Table 1 0.00088
    fv_ent_jejunum <- fixed(0.00192); label("Wall volume of each jejunum segment (fraction of BW; both sexes)") # S1 File VMEMJEJ1C = VMEMJEJ2C = 0.00192 ('Female 0.00192')
    fv_ent_ileum <- fixed(0.00142); label("Wall volume of each ileum segment, male (fraction of BW)") # S1 File VMEMILL1C..3C = 0.00142; Table 1 0.0014
    fv_ent_ileum_f <- fixed(0.00149); label("Wall volume of each ileum segment, female (fraction of BW)") # S1 File comment 'Female 0.00149'; Table 1 0.0015
    fv_ent_cecum <- fixed(0.0004); label("Cecum wall volume (fraction of BW; both sexes)") # S1 File VMEMCECUMMEMC = 0.0004; Table 1 0.0004/0.0004
    fv_ent_colon <- fixed(0.000848); label("Ascending-colon wall volume, male (fraction of BW)") # S1 File VMEMASCENDINGMEMC = 0.000848; Table 1 0.00085
    fv_ent_colon_f <- fixed(0.00093); label("Ascending-colon wall volume, female (fraction of BW)") # S1 File comment '0.00093'; Table 1 0.00093

    # ---- Partition coefficients (tissue:blood) ----
    p_fat <- fixed(1.79); label("Fat:blood partition coefficient (unitless)") # S1 File PFat = 1.79
    p_brain <- fixed(6.07); label("Brain:blood partition coefficient (unitless)") # S1 File Pbrain = 6.07
    p_rich <- fixed(5.66); label("Richly perfused:blood partition coefficient (unitless)") # S1 File PRich = 5.66
    p_slow <- fixed(2.47); label("Slowly perfused:blood partition coefficient (unitless)") # S1 File PSlow = 2.47
    p_gonad <- fixed(3.12); label("Gonad:blood partition coefficient (unitless)") # S1 File PGonad = 3.12
    p_heart <- fixed(2.19); label("Heart:blood partition coefficient (unitless)") # S1 File Pheart = 2.19
    p_liver <- fixed(5.66); label("Liver:blood partition coefficient (unitless)") # S1 File PLiver = 5.66
    p_gut <- fixed(5.66); label("Small-intestine enterocyte:blood partition coefficient (unitless)") # S1 File KPGUT = 5.66; Table 2 'set to the value of the liver'
    p_colon <- fixed(5.66); label("Cecum / ascending-colon wall:blood partition coefficient (unitless)") # S1 File KPcolon = 5.66

    # ---- Gastrointestinal physiology (Table 1; S1 File) ----
    mtt_stomach <- fixed(0.25); label("Stomach transit time, fasted (h)") # Table 1 0.25 (fasting); S1 File TSTOMACH = 0.25
    mtt_stomach_fed <- fixed(1); label("Stomach transit time, fed (h)") # Table 1 0.25/1 (fasting/fed)
    mtt_duodenum <- fixed(0.26); label("Duodenum transit time (h)") # Table 1; S1 File TDUO = 0.26
    mtt_jejunum1 <- fixed(0.95); label("Jejunum 1 transit time (h)") # Table 1; S1 File TJEJ1 = 0.95
    mtt_jejunum2 <- fixed(0.76); label("Jejunum 2 transit time (h)") # Table 1; S1 File TJEJ2 = 0.76
    mtt_ileum1 <- fixed(0.59); label("Ileum 1 transit time (h)") # Table 1; S1 File TILL1 = 0.59
    mtt_ileum2 <- fixed(0.43); label("Ileum 2 transit time (h)") # Table 1; S1 File TILL2 = 0.43
    mtt_ileum3 <- fixed(0.31); label("Ileum 3 transit time (h)") # Table 1; S1 File TILL3 = 0.31
    mtt_cecum <- fixed(4.5); label("Cecum transit time (h)") # Table 1; S1 File Tcecum = 4.5
    mtt_colon <- fixed(13.5); label("Ascending-colon transit time (h)") # Table 1; S1 File Tascending = 13.5
    ph_stomach <- fixed(1.3); label("Stomach pH, fasted") # Table 1 1.3/4.9; S1 File PHSTOMACH = 1.3
    ph_stomach_fed <- fixed(4.9); label("Stomach pH, fed") # Table 1 1.3/4.9 (fasting/fed)
    ph_duodenum <- fixed(6.0); label("Duodenum pH, fasted") # Table 1 6.0/5.4; S1 File PHDUO = 6.0
    ph_duodenum_fed <- fixed(5.4); label("Duodenum pH, fed") # Table 1 6.0/5.4
    ph_jejunum1 <- fixed(6.2); label("Jejunum 1 pH, fasted") # Table 1 6.2/5.4; S1 File PHJEJ1 = 6.2
    ph_jejunum1_fed <- fixed(5.4); label("Jejunum 1 pH, fed") # Table 1 6.2/5.4
    ph_jejunum2 <- fixed(6.4); label("Jejunum 2 pH, fasted") # Table 1 6.4/6.0; S1 File PHJEJ2 = 6.4
    ph_jejunum2_fed <- fixed(6.0); label("Jejunum 2 pH, fed") # Table 1 6.4/6.0
    ph_ileum1 <- fixed(6.6); label("Ileum 1 pH") # Table 1; S1 File PHILL1 = 6.6
    ph_ileum2 <- fixed(6.9); label("Ileum 2 pH") # Table 1; S1 File PHILL2 = 6.9
    ph_ileum3 <- fixed(7.4); label("Ileum 3 pH") # Table 1; S1 File PHILL3 = 7.4
    ph_cecum <- fixed(6.4); label("Cecum pH") # Table 1; S1 File PHcecum = 6.4
    ph_colon <- fixed(6.8); label("Ascending-colon pH") # Table 1; S1 File PHascending = 6.8
    vlumen_stomach <- fixed(50); label("Stomach luminal volume (mL)") # Table 1; S1 File VSTOMACH = 50
    vlumen_duodenum <- fixed(48); label("Duodenum luminal volume (mL)") # Table 1; S1 File VDUO = 48
    vlumen_jejunum1 <- fixed(175); label("Jejunum 1 luminal volume (mL)") # Table 1; S1 File VJEJ1 = 175
    vlumen_jejunum2 <- fixed(140); label("Jejunum 2 luminal volume (mL)") # Table 1; S1 File VJEJ2 = 140
    vlumen_ileum1 <- fixed(109); label("Ileum 1 luminal volume (mL)") # Table 1; S1 File VILL1 = 109
    vlumen_ileum2 <- fixed(79); label("Ileum 2 luminal volume (mL)") # Table 1; S1 File VILL2 = 79
    vlumen_ileum3 <- fixed(56); label("Ileum 3 luminal volume (mL)") # Table 1; S1 File VILL3 = 56
    vlumen_cecum <- fixed(53); label("Cecum luminal volume (mL)") # Table 1; S1 File Vcecum = 53
    vlumen_colon <- fixed(57); label("Ascending-colon luminal volume (mL)") # Table 1; S1 File Vascending = 57
    esa_duodenum <- fixed(19995); label("Duodenum absorptive surface area (cm^2)") # Table 1; S1 File ESA_DUO = 19995
    esa_jejunum1 <- fixed(77482); label("Jejunum 1 absorptive surface area (cm^2)") # Table 1; S1 File ESA_JEJ1 = 77482
    esa_jejunum2 <- fixed(69217); label("Jejunum 2 absorptive surface area (cm^2)") # Table 1; S1 File ESA_JEJ2 = 69217
    esa_ileum1 <- fixed(60952); label("Ileum 1 absorptive surface area (cm^2)") # Table 1; S1 File ESA_ILL1 = 60952
    esa_ileum2 <- fixed(52171); label("Ileum 2 absorptive surface area (cm^2)") # Table 1; S1 File ESA_ILL2 = 52171
    esa_ileum3 <- fixed(43906); label("Ileum 3 absorptive surface area (cm^2)") # Table 1; S1 File ESA_ILL3 = 43906
    esa_cecum <- fixed(1964); label("Cecum absorptive surface area (cm^2)") # Table 1; S1 File ESA_CECUM = 1964
    esa_colon <- fixed(2961); label("Ascending-colon absorptive surface area (cm^2)") # Table 1; S1 File ESA_ASCENDING = 2961
    bile_duodenum <- fixed(2.8); label("Duodenum bile-salt concentration, fasted (mM)") # Table 1 2.8/14.44; S1 File BILE_DUO = 2.8
    bile_duodenum_fed <- fixed(14.44); label("Duodenum bile-salt concentration, fed (mM)") # Table 1 2.8/14.44
    bile_jejunum1 <- fixed(2.33); label("Jejunum 1 bile-salt concentration, fasted (mM)") # Table 1 2.33/12.02; S1 File BILE_JEJ1 = 2.33
    bile_jejunum1_fed <- fixed(12.02); label("Jejunum 1 bile-salt concentration, fed (mM)") # Table 1 2.33/12.02
    bile_jejunum2 <- fixed(2.03); label("Jejunum 2 bile-salt concentration, fasted (mM)") # Table 1 2.03/10.46; S1 File BILE_JEJ2 = 2.03
    bile_jejunum2_fed <- fixed(10.46); label("Jejunum 2 bile-salt concentration, fed (mM)") # Table 1 2.03/10.46
    bile_ileum1 <- fixed(1.41); label("Ileum 1 bile-salt concentration, fasted (mM)") # Table 1 1.41/7.28; S1 File BILE_ILL1 = 1.41
    bile_ileum1_fed <- fixed(7.28); label("Ileum 1 bile-salt concentration, fed (mM)") # Table 1 1.41/7.28
    bile_ileum2 <- fixed(1.16); label("Ileum 2 bile-salt concentration, fasted (mM)") # Table 1 1.16/5.99; S1 File BILE_ILL2 = 1.16
    bile_ileum2_fed <- fixed(5.99); label("Ileum 2 bile-salt concentration, fed (mM)") # Table 1 1.16/5.99
    bile_ileum3 <- fixed(0.14); label("Ileum 3 bile-salt concentration, fasted (mM)") # Table 1 0.14/0.73; S1 File BILE_ILL3 = 0.14
    bile_ileum3_fed <- fixed(0.73); label("Ileum 3 bile-salt concentration, fed (mM)") # Table 1 0.14/0.73

    # ---- Methylphenidate physicochemical properties (Table 2; Eqs 1-6) ----
    mw <- fixed(233.3062); label("MPH molecular weight (g/mol)") # Table 2; S1 File MW = 233.3062
    mw_water <- fixed(18.02); label("Water molecular weight (g/mol)") # S1 File MW_WATER = 18.02
    mw_ra <- fixed(220); label("Ritalinic acid molecular weight used for the MPH -> RA mass conversion (g/mol)") # S1 File '220/233.3062'
    pka <- fixed(8.77); label("MPH pKa (base)") # Table 2; S1 File PKABASE1 = 8.77
    logp <- fixed(0.20); label("MPH octanol/water logP (log10 units)") # Table 2; S1 File LOGP = 0.2
    sol_ref <- fixed(9250); label("MPH aqueous solubility at the reference pH (mg/L)") # Table 2 9250 mg/L at pH 10.58; S1 File SOL = 9250
    ph_sol_ref <- fixed(10.58); label("Reference pH of the solubility value") # Table 2; S1 File REFPHSOL = 10.58
    hpeff <- fixed(0.963e-4); label("Human effective jejunal permeability (cm/s)") # Table 2 0.963 x 10^-4 cm/s; Results 9.63e-5 cm/s (optimized)
    diffcoeff <- fixed(1e-4); label("Drug diffusion coefficient (cm^2/min)") # Methods Eq 1 text; S1 File DIFFCOEFF_IR = 1e-4
    rho <- fixed(1); label("Drug particle density (g/mL)") # Methods Eq 1 text; S1 File rho = 1
    rparticle <- fixed(5); label("Drug particle radius (um)") # Methods Eq 1 text; S1 File rparticle = 5.0
    dlt <- fixed(30); label("Diffusion layer thickness (um)") # Methods Eq 1 text; S1 File DLT = 30

    # ---- Gut non-specific loss and hepatic metabolism (Table 2) ----
    k5_dc <- fixed(0.79); label("d-MPH non-specific loss in enterocytes (1/h/kg^0.75)") # Table 2 K5dC = 0.79 (optimization); Results
    k5_lc <- fixed(37.8); label("l-MPH non-specific loss in enterocytes (1/h/kg^0.75)") # Table 2 K5lC = 37.8 (optimization); Results
    f_hydrol_gut <- fixed(0.8); label("Fraction of gut loss forming ritalinic acid (unitless)") # S1 File F = 0.8
    km_liver_d <- fixed(27600); label("d-MPH hepatic hydrolysis Km (ug/L)") # Table 2 Kmliverd = 27,600
    km_liver_l <- fixed(10172); label("l-MPH hepatic hydrolysis Km (ug/L)") # Table 2 Kmliverl = 10,172
    vmax_liver_dc <- fixed(25826); label("d-MPH hepatic hydrolysis Vmax (ug/h/kg^0.75)") # Table 2 VmaxliverdC = 25,826
    vmax_liver_lc <- fixed(52404); label("l-MPH hepatic hydrolysis Vmax (ug/h/kg^0.75)") # Table 2 VmaxliverlC = 52,404
    kmet_dc <- fixed(0.43); label("d-MPH hepatic oxidation clearance (L/h/kg^0.75)") # Table 2 KmetdC = 0.43
    kmet_lc <- fixed(0.43); label("l-MPH hepatic oxidation clearance (L/h/kg^0.75)") # Table 2 KmetlC = 0.43
    ku_ra_dc <- fixed(0.305); label("d-ritalinic acid urinary clearance (L/h/kg^0.75)") # S1 File Ku_RAdC = 0.305
    ku_ra_lc <- fixed(0.168); label("l-ritalinic acid urinary clearance (L/h/kg^0.75)") # S1 File Ku_RAlC = 0.168

    # ---- ER bead release (Weibull, Eq 2; Table 2 'Visual Fit') ----
    tlag_er <- fixed(4); label("ER release lag time (h)") # Table 2 Tlag = 4 h (all ER products)
    wb_a <- fixed(0.5); label("Weibull scale A, Ritalin LA / Metadate CD / Medikinet Retard (h^b)") # Table 2 A = 0.5
    wb_b <- fixed(2); label("Weibull shape b, Ritalin LA / Metadate CD / Medikinet Retard (unitless)") # Table 2 b = 2
    wb_a_mlr <- fixed(6); label("Weibull scale A, MLR-MPH (h^b)") # Table 2 A = 6
    wb_b_mlr <- fixed(4); label("Weibull shape b, MLR-MPH (unitless)") # Table 2 b = 4

    # ---- Residual error (not reported; deterministic model) ----
    propSd_d <- fixed(0); label("Proportional residual error, d-MPH plasma (fraction; not reported)") # not reported
    propSd_l <- fixed(0); label("Proportional residual error, l-MPH plasma (fraction; not reported)") # not reported
    propSd_rac <- fixed(0); label("Proportional residual error, total MPH plasma (fraction; not reported)") # not reported
  })

  model({
    # ---- Body-weight and sex scaling (S1 File DERIVATIVE block) ----
    bw075 <- WT^0.75
    qc <- qcc * bw075
    qfat <- (fq_fat * (1 - SEXF) + fq_fat_f * SEXF) * qc
    qbrain <- fq_brain * qc
    qgonad <- (fq_gonad * (1 - SEXF) + fq_gonad_f * SEXF) * qc
    qheart <- (fq_heart * (1 - SEXF) + fq_heart_f * SEXF) * qc
    qliver <- (fq_liver * (1 - SEXF) + fq_liver_f * SEXF) * qc
    q_seg <- (fq_seg * (1 - SEXF) + fq_seg_f * SEXF) * qc
    qslow <- fq_slow_tot * qc - qfat - qheart
    qrich <- fq_rich_tot * qc - qliver - qgonad - qbrain

    vplasma <- fv_plasma * WT
    vfat <- (fv_fat * (1 - SEXF) + fv_fat_f * SEXF) * WT
    vbrain <- fv_brain * WT
    vgonad <- (fv_gonad * (1 - SEXF) + fv_gonad_f * SEXF) * WT
    vheart <- (fv_heart * (1 - SEXF) + fv_heart_f * SEXF) * WT
    vliver <- fv_liver * WT
    vbody <- fv_body * WT
    vent_duodenum <- (fv_ent_duodenum * (1 - SEXF) + fv_ent_duodenum_f * SEXF) * WT
    vent_jejunum1 <- fv_ent_jejunum * WT
    vent_jejunum2 <- fv_ent_jejunum * WT
    vent_ileum1 <- (fv_ent_ileum * (1 - SEXF) + fv_ent_ileum_f * SEXF) * WT
    vent_ileum2 <- vent_ileum1
    vent_ileum3 <- vent_ileum1
    vent_cecum <- fv_ent_cecum * WT
    vent_colon <- (fv_ent_colon * (1 - SEXF) + fv_ent_colon_f * SEXF) * WT
    vslow <- fv_slow_tot * WT - vfat - vheart
    vrich <- fv_rich_tot * WT - vliver - vplasma - vgonad - vbrain - vent_duodenum -
      vent_jejunum1 - vent_jejunum2 - vent_ileum1 - vent_ileum2 - vent_ileum3 -
      vent_cecum - vent_colon

    k5_d <- k5_dc * bw075
    k5_l <- k5_lc * bw075
    vmax_liver_d <- vmax_liver_dc * bw075
    vmax_liver_l <- vmax_liver_lc * bw075
    kmet_d <- kmet_dc * bw075
    kmet_l <- kmet_lc * bw075
    ku_ra_d <- ku_ra_dc * bw075
    ku_ra_l <- ku_ra_lc * bw075

    # ---- Fed-state switches (Table 1 footnote f) ----
    mtt_stomach_i <- mtt_stomach * (1 - FED) + mtt_stomach_fed * FED
    ph_stomach_i <- ph_stomach * (1 - FED) + ph_stomach_fed * FED
    ph_duodenum_i <- ph_duodenum * (1 - FED) + ph_duodenum_fed * FED
    ph_jejunum1_i <- ph_jejunum1 * (1 - FED) + ph_jejunum1_fed * FED
    ph_jejunum2_i <- ph_jejunum2 * (1 - FED) + ph_jejunum2_fed * FED
    bile_duodenum_i <- bile_duodenum * (1 - FED) + bile_duodenum_fed * FED
    bile_jejunum1_i <- bile_jejunum1 * (1 - FED) + bile_jejunum1_fed * FED
    bile_jejunum2_i <- bile_jejunum2 * (1 - FED) + bile_jejunum2_fed * FED
    bile_ileum1_i <- bile_ileum1 * (1 - FED) + bile_ileum1_fed * FED
    bile_ileum2_i <- bile_ileum2 * (1 - FED) + bile_ileum2_fed * FED
    bile_ileum3_i <- bile_ileum3 * (1 - FED) + bile_ileum3_fed * FED

    # ---- Ionization and solubility (Eqs 3-5 as coded in S1 File) ----
    # Monoprotic base: HH = 1 + 10^(pKa - pH); unionized fraction NI = 1/HH.
    hh_ref <- 1 + 10^(pka - ph_sol_ref)
    sol_int <- sol_ref / hh_ref
    hh_stomach <- 1 + 10^(pka - ph_stomach_i)
    hh_duodenum <- 1 + 10^(pka - ph_duodenum_i)
    hh_jejunum1 <- 1 + 10^(pka - ph_jejunum1_i)
    hh_jejunum2 <- 1 + 10^(pka - ph_jejunum2_i)
    hh_ileum1 <- 1 + 10^(pka - ph_ileum1)
    hh_ileum2 <- 1 + 10^(pka - ph_ileum2)
    hh_ileum3 <- 1 + 10^(pka - ph_ileum3)
    hh_cecum <- 1 + 10^(pka - ph_cecum)
    hh_colon <- 1 + 10^(pka - ph_colon)
    # Bile solubilization ratio (Mithani 1996); the S1 File applies exp() to
    # log10 SR and halves the summed solubility, kept as run.
    sr_bile <- exp(0.606 * logp + 2.234)
    scaq_fac <- sr_bile * (mw_water / mw) * 1e-6 * mw
    solif_stomach <- sol_int * hh_stomach * 0.5
    solif_duodenum <- sol_int * hh_duodenum * (1 + scaq_fac * bile_duodenum_i) * 0.5
    solif_jejunum1 <- sol_int * hh_jejunum1 * (1 + scaq_fac * bile_jejunum1_i) * 0.5
    solif_jejunum2 <- sol_int * hh_jejunum2 * (1 + scaq_fac * bile_jejunum2_i) * 0.5
    solif_ileum1 <- sol_int * hh_ileum1 * (1 + scaq_fac * bile_ileum1_i) * 0.5
    solif_ileum2 <- sol_int * hh_ileum2 * (1 + scaq_fac * bile_ileum2_i) * 0.5
    solif_ileum3 <- sol_int * hh_ileum3 * (1 + scaq_fac * bile_ileum3_i) * 0.5
    solif_cecum <- sol_int * hh_cecum * 0.5
    solif_colon <- sol_int * hh_colon * 0.5
    ni_duodenum <- 1 / hh_duodenum
    ni_jejunum1 <- 1 / hh_jejunum1
    ni_jejunum2 <- 1 / hh_jejunum2
    ni_ileum1 <- 1 / hh_ileum1
    ni_ileum2 <- 1 / hh_ileum2
    ni_ileum3 <- 1 / hh_ileum3
    ni_cecum <- 1 / hh_cecum
    ni_colon <- 1 / hh_colon

    # ---- Permeability clearance (Eq 6) and dissolution constant (Eq 1) ----
    # diff_<seg> (mL/h) = 3600 s/h * Peff (cm/s) * surface area (cm^2).
    diff_duodenum <- 3600 * hpeff * esa_duodenum
    diff_jejunum1 <- 3600 * hpeff * esa_jejunum1
    diff_jejunum2 <- 3600 * hpeff * esa_jejunum2
    diff_ileum1 <- 3600 * hpeff * esa_ileum1
    diff_ileum2 <- 3600 * hpeff * esa_ileum2
    diff_ileum3 <- 3600 * hpeff * esa_ileum3
    diff_cecum <- 3600 * hpeff * esa_cecum
    diff_colon <- 3600 * hpeff * esa_colon
    # Eq 1 with its 60 x 100 factor; dissolution rate = kd_ir * X_solid * (Sol - Cdiss).
    kd_ir <- 3 * (diffcoeff / (rho * rparticle * dlt)) * 60 * 100

    # ---- ER bead release (S1 File) ----
    # The release clock is the time since the last ER dose (the S1 File uses
    # simulation time, identical for its single dose at time 0).
    # The S1 File releases X_ER(segment) * f(t) per hour, where f is the
    # Weibull DENSITY of Eq 2 (not its hazard), so at most 1 - exp(-1) = 63%
    # of the ER load is released; Table 3 reports 61-62% released.
    wb_a_i <- wb_a * (1 - FORM_METHYLPHENIDATE_MLR) + wb_a_mlr * FORM_METHYLPHENIDATE_MLR
    wb_b_i <- wb_b * (1 - FORM_METHYLPHENIDATE_MLR) + wb_b_mlr * FORM_METHYLPHENIDATE_MLR
    t_er <- tad(stomach_er_d) - tlag_er
    wb_er <- 0
    if (t_er > 0) {
      wb_er <- (wb_b_i / wb_a_i) * t_er^(wb_b_i - 1) * exp(-(t_er^wb_b_i) / wb_a_i)
    }

    # D-MPH arterial plasma and venous-blood tissue concentrations (ug/L)
    ca_d <- central_d / vplasma
    cv_liver_d <- liver_d / (vliver * p_liver)
    cv_fat_d <- adipose_d / (vfat * p_fat)
    cv_gonad_d <- gonads_d / (vgonad * p_gonad)
    cv_brain_d <- brain_d / (vbrain * p_brain)
    cv_rich_d <- a_rapidly_perfused_d / (vrich * p_rich)
    cv_slow_d <- a_slowly_perfused_d / (vslow * p_slow)
    cv_heart_d <- heart_d / (vheart * p_heart)
    cv_d <- (cv_liver_d * qliver + cv_fat_d * qfat + cv_rich_d * qrich +
      cv_slow_d * qslow + cv_gonad_d * qgonad + cv_brain_d * qbrain + cv_heart_d * qheart) / qc
    # L-MPH arterial plasma and venous-blood tissue concentrations (ug/L)
    ca_l <- central_l / vplasma
    cv_liver_l <- liver_l / (vliver * p_liver)
    cv_fat_l <- adipose_l / (vfat * p_fat)
    cv_gonad_l <- gonads_l / (vgonad * p_gonad)
    cv_brain_l <- brain_l / (vbrain * p_brain)
    cv_rich_l <- a_rapidly_perfused_l / (vrich * p_rich)
    cv_slow_l <- a_slowly_perfused_l / (vslow * p_slow)
    cv_heart_l <- heart_l / (vheart * p_heart)
    cv_l <- (cv_liver_l * qliver + cv_fat_l * qfat + cv_rich_l * qrich +
      cv_slow_l * qslow + cv_gonad_l * qgonad + cv_brain_l * qbrain + cv_heart_l * qheart) / qc

    # ---- D-MPH: stomach (no ER release in the stomach; code RELEASE_stomach = 0) ----
    diss_stomach_d <- kd_ir * stomach_ir_d * (solif_stomach - stomach_lumen_d / vlumen_stomach)
    d/dt(stomach_er_d) <- -stomach_er_d / mtt_stomach_i
    d/dt(stomach_ir_d) <- -diss_stomach_d - stomach_ir_d / mtt_stomach_i
    d/dt(stomach_lumen_d) <- diss_stomach_d - stomach_lumen_d / mtt_stomach_i

    # ---- D-MPH: duodenum ----
    diss_duodenum_d <- kd_ir * duodenum_ir_d * (solif_duodenum - duodenum_lumen_d / vlumen_duodenum)
    rel_duodenum_d <- duodenum_er_d * wb_er
    abs_duodenum_d <- ni_duodenum * diff_duodenum * (duodenum_lumen_d / vlumen_duodenum - duodenum_ent_d / vent_duodenum / 1000)
    d/dt(duodenum_er_d) <- stomach_er_d / mtt_stomach_i - duodenum_er_d / mtt_duodenum - rel_duodenum_d
    d/dt(duodenum_ir_d) <- stomach_ir_d / mtt_stomach_i - duodenum_ir_d / mtt_duodenum - diss_duodenum_d
    d/dt(duodenum_lumen_d) <- stomach_lumen_d / mtt_stomach_i - duodenum_lumen_d / mtt_duodenum + diss_duodenum_d + rel_duodenum_d - abs_duodenum_d
    out_duodenum_d <- q_seg * duodenum_ent_d / (vent_duodenum * p_gut)
    d/dt(duodenum_ent_d) <- q_seg * ca_d - out_duodenum_d - k5_d * duodenum_ent_d + abs_duodenum_d

    # ---- D-MPH: jejunum1 ----
    diss_jejunum1_d <- kd_ir * jejunum1_ir_d * (solif_jejunum1 - jejunum1_lumen_d / vlumen_jejunum1)
    rel_jejunum1_d <- jejunum1_er_d * wb_er
    abs_jejunum1_d <- ni_jejunum1 * diff_jejunum1 * (jejunum1_lumen_d / vlumen_jejunum1 - jejunum1_ent_d / vent_jejunum1 / 1000)
    d/dt(jejunum1_er_d) <- duodenum_er_d / mtt_duodenum - jejunum1_er_d / mtt_jejunum1 - rel_jejunum1_d
    d/dt(jejunum1_ir_d) <- duodenum_ir_d / mtt_duodenum - jejunum1_ir_d / mtt_jejunum1 - diss_jejunum1_d
    d/dt(jejunum1_lumen_d) <- duodenum_lumen_d / mtt_duodenum - jejunum1_lumen_d / mtt_jejunum1 + diss_jejunum1_d + rel_jejunum1_d - abs_jejunum1_d
    out_jejunum1_d <- q_seg * jejunum1_ent_d / (vent_jejunum1 * p_gut)
    d/dt(jejunum1_ent_d) <- q_seg * ca_d - out_jejunum1_d - k5_d * jejunum1_ent_d + abs_jejunum1_d

    # ---- D-MPH: jejunum2 ----
    diss_jejunum2_d <- kd_ir * jejunum2_ir_d * (solif_jejunum2 - jejunum2_lumen_d / vlumen_jejunum2)
    rel_jejunum2_d <- jejunum2_er_d * wb_er
    abs_jejunum2_d <- ni_jejunum2 * diff_jejunum2 * (jejunum2_lumen_d / vlumen_jejunum2 - jejunum2_ent_d / vent_jejunum2 / 1000)
    d/dt(jejunum2_er_d) <- jejunum1_er_d / mtt_jejunum1 - jejunum2_er_d / mtt_jejunum2 - rel_jejunum2_d
    d/dt(jejunum2_ir_d) <- jejunum1_ir_d / mtt_jejunum1 - jejunum2_ir_d / mtt_jejunum2 - diss_jejunum2_d
    d/dt(jejunum2_lumen_d) <- jejunum1_lumen_d / mtt_jejunum1 - jejunum2_lumen_d / mtt_jejunum2 + diss_jejunum2_d + rel_jejunum2_d - abs_jejunum2_d
    out_jejunum2_d <- q_seg * jejunum2_ent_d / (vent_jejunum2 * p_gut)
    d/dt(jejunum2_ent_d) <- q_seg * ca_d - out_jejunum2_d - k5_d * jejunum2_ent_d + abs_jejunum2_d

    # ---- D-MPH: ileum1 ----
    diss_ileum1_d <- kd_ir * ileum1_ir_d * (solif_ileum1 - ileum1_lumen_d / vlumen_ileum1)
    rel_ileum1_d <- ileum1_er_d * wb_er
    abs_ileum1_d <- ni_ileum1 * diff_ileum1 * (ileum1_lumen_d / vlumen_ileum1 - ileum1_ent_d / vent_ileum1 / 1000)
    d/dt(ileum1_er_d) <- jejunum2_er_d / mtt_jejunum2 - ileum1_er_d / mtt_ileum1 - rel_ileum1_d
    d/dt(ileum1_ir_d) <- jejunum2_ir_d / mtt_jejunum2 - ileum1_ir_d / mtt_ileum1 - diss_ileum1_d
    d/dt(ileum1_lumen_d) <- jejunum2_lumen_d / mtt_jejunum2 - ileum1_lumen_d / mtt_ileum1 + diss_ileum1_d + rel_ileum1_d - abs_ileum1_d
    out_ileum1_d <- q_seg * ileum1_ent_d / (vent_ileum1 * p_gut)
    d/dt(ileum1_ent_d) <- q_seg * ca_d - out_ileum1_d - k5_d * ileum1_ent_d + abs_ileum1_d

    # ---- D-MPH: ileum2 ----
    diss_ileum2_d <- kd_ir * ileum2_ir_d * (solif_ileum2 - ileum2_lumen_d / vlumen_ileum2)
    rel_ileum2_d <- ileum2_er_d * wb_er
    abs_ileum2_d <- ni_ileum2 * diff_ileum2 * (ileum2_lumen_d / vlumen_ileum2 - ileum2_ent_d / vent_ileum2 / 1000)
    d/dt(ileum2_er_d) <- ileum1_er_d / mtt_ileum1 - ileum2_er_d / mtt_ileum2 - rel_ileum2_d
    d/dt(ileum2_ir_d) <- ileum1_ir_d / mtt_ileum1 - ileum2_ir_d / mtt_ileum2 - diss_ileum2_d
    d/dt(ileum2_lumen_d) <- ileum1_lumen_d / mtt_ileum1 - ileum2_lumen_d / mtt_ileum2 + diss_ileum2_d + rel_ileum2_d - abs_ileum2_d
    out_ileum2_d <- q_seg * ileum2_ent_d / (vent_ileum2 * p_gut)
    d/dt(ileum2_ent_d) <- q_seg * ca_d - out_ileum2_d - k5_d * ileum2_ent_d + abs_ileum2_d

    # ---- D-MPH: ileum3 ----
    diss_ileum3_d <- kd_ir * ileum3_ir_d * (solif_ileum3 - ileum3_lumen_d / vlumen_ileum3)
    rel_ileum3_d <- ileum3_er_d * wb_er
    abs_ileum3_d <- ni_ileum3 * diff_ileum3 * (ileum3_lumen_d / vlumen_ileum3 - ileum3_ent_d / vent_ileum3 / 1000)
    d/dt(ileum3_er_d) <- ileum2_er_d / mtt_ileum2 - ileum3_er_d / mtt_ileum3 - rel_ileum3_d
    d/dt(ileum3_ir_d) <- ileum2_ir_d / mtt_ileum2 - ileum3_ir_d / mtt_ileum3 - diss_ileum3_d
    d/dt(ileum3_lumen_d) <- ileum2_lumen_d / mtt_ileum2 - ileum3_lumen_d / mtt_ileum3 + diss_ileum3_d + rel_ileum3_d - abs_ileum3_d
    out_ileum3_d <- q_seg * ileum3_ent_d / (vent_ileum3 * p_gut)
    d/dt(ileum3_ent_d) <- q_seg * ca_d - out_ileum3_d - k5_d * ileum3_ent_d + abs_ileum3_d

    # ---- D-MPH: cecum ----
    diss_cecum_d <- kd_ir * cecum_ir_d * (solif_cecum - cecum_lumen_d / vlumen_cecum)
    rel_cecum_d <- cecum_er_d * wb_er
    abs_cecum_d <- ni_cecum * diff_cecum * (cecum_lumen_d / vlumen_cecum - cecum_ent_d / vent_cecum / 1000)
    d/dt(cecum_er_d) <- ileum3_er_d / mtt_ileum3 - cecum_er_d / mtt_cecum - rel_cecum_d
    d/dt(cecum_ir_d) <- ileum3_ir_d / mtt_ileum3 - cecum_ir_d / mtt_cecum - diss_cecum_d
    d/dt(cecum_lumen_d) <- ileum3_lumen_d / mtt_ileum3 - cecum_lumen_d / mtt_cecum + diss_cecum_d + rel_cecum_d - abs_cecum_d
    out_cecum_d <- q_seg * cecum_ent_d / (vent_cecum * p_colon)
    d/dt(cecum_ent_d) <- q_seg * ca_d - out_cecum_d - k5_d * cecum_ent_d + abs_cecum_d

    # ---- D-MPH: colon ----
    diss_colon_d <- kd_ir * colon_ir_d * (solif_colon - colon_lumen_d / vlumen_colon)
    rel_colon_d <- colon_er_d * wb_er
    abs_colon_d <- ni_colon * diff_colon * (colon_lumen_d / vlumen_colon - colon_ent_d / vent_colon / 1000)
    d/dt(colon_er_d) <- cecum_er_d / mtt_cecum - colon_er_d / mtt_colon - rel_colon_d
    d/dt(colon_ir_d) <- cecum_ir_d / mtt_cecum - colon_ir_d / mtt_colon - diss_colon_d
    d/dt(colon_lumen_d) <- cecum_lumen_d / mtt_cecum - colon_lumen_d / mtt_colon + diss_colon_d + rel_colon_d - abs_colon_d
    out_colon_d <- q_seg * colon_ent_d / (vent_colon * p_colon)
    d/dt(colon_ent_d) <- q_seg * ca_d - out_colon_d - k5_d * colon_ent_d + abs_colon_d

    # ---- D-MPH: systemic PBPK (flow-limited; concentrations ug/L) ----
    rbs_d <- out_duodenum_d + out_jejunum1_d + out_jejunum2_d + out_ileum1_d + out_ileum2_d + out_ileum3_d + out_cecum_d + out_colon_d
    gutloss_d <- k5_d * (duodenum_ent_d + jejunum1_ent_d + jejunum2_ent_d + ileum1_ent_d + ileum2_ent_d + ileum3_ent_d + cecum_ent_d + colon_ent_d)
    rox_d <- kmet_d * cv_liver_d
    rhyd_d <- vmax_liver_d * cv_liver_d / (km_liver_d * (1 + cv_liver_l / km_liver_l) + cv_liver_d)
    d/dt(central_d) <- qc * (cv_d - ca_d)
    d/dt(adipose_d) <- qfat * (ca_d - cv_fat_d)
    d/dt(gonads_d) <- qgonad * (ca_d - cv_gonad_d)
    d/dt(liver_d) <- (qliver - 8 * q_seg) * ca_d - qliver * cv_liver_d - rox_d - rhyd_d + rbs_d
    d/dt(brain_d) <- qbrain * (ca_d - cv_brain_d)
    d/dt(a_rapidly_perfused_d) <- qrich * (ca_d - cv_rich_d)
    d/dt(a_slowly_perfused_d) <- qslow * (ca_d - cv_slow_d)
    d/dt(heart_d) <- qheart * (ca_d - cv_heart_d)
    # ritalinic acid formed in the liver and from the gut loss (fraction f_hydrol_gut)
    d/dt(central_ra_d) <- (rhyd_d + f_hydrol_gut * gutloss_d) * mw_ra / mw - ku_ra_d * central_ra_d / vbody

    # ---- L-MPH: stomach (no ER release in the stomach; code RELEASE_stomach = 0) ----
    diss_stomach_l <- kd_ir * stomach_ir_l * (solif_stomach - stomach_lumen_l / vlumen_stomach)
    d/dt(stomach_er_l) <- -stomach_er_l / mtt_stomach_i
    d/dt(stomach_ir_l) <- -diss_stomach_l - stomach_ir_l / mtt_stomach_i
    d/dt(stomach_lumen_l) <- diss_stomach_l - stomach_lumen_l / mtt_stomach_i

    # ---- L-MPH: duodenum ----
    diss_duodenum_l <- kd_ir * duodenum_ir_l * (solif_duodenum - duodenum_lumen_l / vlumen_duodenum)
    rel_duodenum_l <- duodenum_er_l * wb_er
    abs_duodenum_l <- ni_duodenum * diff_duodenum * (duodenum_lumen_l / vlumen_duodenum - duodenum_ent_l / vent_duodenum / 1000)
    d/dt(duodenum_er_l) <- stomach_er_l / mtt_stomach_i - duodenum_er_l / mtt_duodenum - rel_duodenum_l
    d/dt(duodenum_ir_l) <- stomach_ir_l / mtt_stomach_i - duodenum_ir_l / mtt_duodenum - diss_duodenum_l
    d/dt(duodenum_lumen_l) <- stomach_lumen_l / mtt_stomach_i - duodenum_lumen_l / mtt_duodenum + diss_duodenum_l + rel_duodenum_l - abs_duodenum_l
    out_duodenum_l <- q_seg * duodenum_ent_l / (vent_duodenum * p_gut)
    d/dt(duodenum_ent_l) <- q_seg * ca_l - out_duodenum_l - k5_l * duodenum_ent_l + abs_duodenum_l

    # ---- L-MPH: jejunum1 ----
    diss_jejunum1_l <- kd_ir * jejunum1_ir_l * (solif_jejunum1 - jejunum1_lumen_l / vlumen_jejunum1)
    rel_jejunum1_l <- jejunum1_er_l * wb_er
    abs_jejunum1_l <- ni_jejunum1 * diff_jejunum1 * (jejunum1_lumen_l / vlumen_jejunum1 - jejunum1_ent_l / vent_jejunum1 / 1000)
    d/dt(jejunum1_er_l) <- duodenum_er_l / mtt_duodenum - jejunum1_er_l / mtt_jejunum1 - rel_jejunum1_l
    d/dt(jejunum1_ir_l) <- duodenum_ir_l / mtt_duodenum - jejunum1_ir_l / mtt_jejunum1 - diss_jejunum1_l
    d/dt(jejunum1_lumen_l) <- duodenum_lumen_l / mtt_duodenum - jejunum1_lumen_l / mtt_jejunum1 + diss_jejunum1_l + rel_jejunum1_l - abs_jejunum1_l
    out_jejunum1_l <- q_seg * jejunum1_ent_l / (vent_jejunum1 * p_gut)
    d/dt(jejunum1_ent_l) <- q_seg * ca_l - out_jejunum1_l - k5_l * jejunum1_ent_l + abs_jejunum1_l

    # ---- L-MPH: jejunum2 ----
    diss_jejunum2_l <- kd_ir * jejunum2_ir_l * (solif_jejunum2 - jejunum2_lumen_l / vlumen_jejunum2)
    rel_jejunum2_l <- jejunum2_er_l * wb_er
    abs_jejunum2_l <- ni_jejunum2 * diff_jejunum2 * (jejunum2_lumen_l / vlumen_jejunum2 - jejunum2_ent_l / vent_jejunum2 / 1000)
    d/dt(jejunum2_er_l) <- jejunum1_er_l / mtt_jejunum1 - jejunum2_er_l / mtt_jejunum2 - rel_jejunum2_l
    d/dt(jejunum2_ir_l) <- jejunum1_ir_l / mtt_jejunum1 - jejunum2_ir_l / mtt_jejunum2 - diss_jejunum2_l
    d/dt(jejunum2_lumen_l) <- jejunum1_lumen_l / mtt_jejunum1 - jejunum2_lumen_l / mtt_jejunum2 + diss_jejunum2_l + rel_jejunum2_l - abs_jejunum2_l
    out_jejunum2_l <- q_seg * jejunum2_ent_l / (vent_jejunum2 * p_gut)
    d/dt(jejunum2_ent_l) <- q_seg * ca_l - out_jejunum2_l - k5_l * jejunum2_ent_l + abs_jejunum2_l

    # ---- L-MPH: ileum1 ----
    diss_ileum1_l <- kd_ir * ileum1_ir_l * (solif_ileum1 - ileum1_lumen_l / vlumen_ileum1)
    rel_ileum1_l <- ileum1_er_l * wb_er
    abs_ileum1_l <- ni_ileum1 * diff_ileum1 * (ileum1_lumen_l / vlumen_ileum1 - ileum1_ent_l / vent_ileum1 / 1000)
    d/dt(ileum1_er_l) <- jejunum2_er_l / mtt_jejunum2 - ileum1_er_l / mtt_ileum1 - rel_ileum1_l
    d/dt(ileum1_ir_l) <- jejunum2_ir_l / mtt_jejunum2 - ileum1_ir_l / mtt_ileum1 - diss_ileum1_l
    d/dt(ileum1_lumen_l) <- jejunum2_lumen_l / mtt_jejunum2 - ileum1_lumen_l / mtt_ileum1 + diss_ileum1_l + rel_ileum1_l - abs_ileum1_l
    out_ileum1_l <- q_seg * ileum1_ent_l / (vent_ileum1 * p_gut)
    d/dt(ileum1_ent_l) <- q_seg * ca_l - out_ileum1_l - k5_l * ileum1_ent_l + abs_ileum1_l

    # ---- L-MPH: ileum2 ----
    diss_ileum2_l <- kd_ir * ileum2_ir_l * (solif_ileum2 - ileum2_lumen_l / vlumen_ileum2)
    rel_ileum2_l <- ileum2_er_l * wb_er
    abs_ileum2_l <- ni_ileum2 * diff_ileum2 * (ileum2_lumen_l / vlumen_ileum2 - ileum2_ent_l / vent_ileum2 / 1000)
    d/dt(ileum2_er_l) <- ileum1_er_l / mtt_ileum1 - ileum2_er_l / mtt_ileum2 - rel_ileum2_l
    d/dt(ileum2_ir_l) <- ileum1_ir_l / mtt_ileum1 - ileum2_ir_l / mtt_ileum2 - diss_ileum2_l
    d/dt(ileum2_lumen_l) <- ileum1_lumen_l / mtt_ileum1 - ileum2_lumen_l / mtt_ileum2 + diss_ileum2_l + rel_ileum2_l - abs_ileum2_l
    out_ileum2_l <- q_seg * ileum2_ent_l / (vent_ileum2 * p_gut)
    d/dt(ileum2_ent_l) <- q_seg * ca_l - out_ileum2_l - k5_l * ileum2_ent_l + abs_ileum2_l

    # ---- L-MPH: ileum3 ----
    diss_ileum3_l <- kd_ir * ileum3_ir_l * (solif_ileum3 - ileum3_lumen_l / vlumen_ileum3)
    rel_ileum3_l <- ileum3_er_l * wb_er
    abs_ileum3_l <- ni_ileum3 * diff_ileum3 * (ileum3_lumen_l / vlumen_ileum3 - ileum3_ent_l / vent_ileum3 / 1000)
    d/dt(ileum3_er_l) <- ileum2_er_l / mtt_ileum2 - ileum3_er_l / mtt_ileum3 - rel_ileum3_l
    d/dt(ileum3_ir_l) <- ileum2_ir_l / mtt_ileum2 - ileum3_ir_l / mtt_ileum3 - diss_ileum3_l
    d/dt(ileum3_lumen_l) <- ileum2_lumen_l / mtt_ileum2 - ileum3_lumen_l / mtt_ileum3 + diss_ileum3_l + rel_ileum3_l - abs_ileum3_l
    out_ileum3_l <- q_seg * ileum3_ent_l / (vent_ileum3 * p_gut)
    d/dt(ileum3_ent_l) <- q_seg * ca_l - out_ileum3_l - k5_l * ileum3_ent_l + abs_ileum3_l

    # ---- L-MPH: cecum ----
    diss_cecum_l <- kd_ir * cecum_ir_l * (solif_cecum - cecum_lumen_l / vlumen_cecum)
    rel_cecum_l <- cecum_er_l * wb_er
    abs_cecum_l <- ni_cecum * diff_cecum * (cecum_lumen_l / vlumen_cecum - cecum_ent_l / vent_cecum / 1000)
    d/dt(cecum_er_l) <- ileum3_er_l / mtt_ileum3 - cecum_er_l / mtt_cecum - rel_cecum_l
    d/dt(cecum_ir_l) <- ileum3_ir_l / mtt_ileum3 - cecum_ir_l / mtt_cecum - diss_cecum_l
    d/dt(cecum_lumen_l) <- ileum3_lumen_l / mtt_ileum3 - cecum_lumen_l / mtt_cecum + diss_cecum_l + rel_cecum_l - abs_cecum_l
    out_cecum_l <- q_seg * cecum_ent_l / (vent_cecum * p_colon)
    d/dt(cecum_ent_l) <- q_seg * ca_l - out_cecum_l - k5_l * cecum_ent_l + abs_cecum_l

    # ---- L-MPH: colon ----
    diss_colon_l <- kd_ir * colon_ir_l * (solif_colon - colon_lumen_l / vlumen_colon)
    rel_colon_l <- colon_er_l * wb_er
    abs_colon_l <- ni_colon * diff_colon * (colon_lumen_l / vlumen_colon - colon_ent_l / vent_colon / 1000)
    d/dt(colon_er_l) <- cecum_er_l / mtt_cecum - colon_er_l / mtt_colon - rel_colon_l
    d/dt(colon_ir_l) <- cecum_ir_l / mtt_cecum - colon_ir_l / mtt_colon - diss_colon_l
    d/dt(colon_lumen_l) <- cecum_lumen_l / mtt_cecum - colon_lumen_l / mtt_colon + diss_colon_l + rel_colon_l - abs_colon_l
    out_colon_l <- q_seg * colon_ent_l / (vent_colon * p_colon)
    d/dt(colon_ent_l) <- q_seg * ca_l - out_colon_l - k5_l * colon_ent_l + abs_colon_l

    # ---- L-MPH: systemic PBPK (flow-limited; concentrations ug/L) ----
    rbs_l <- out_duodenum_l + out_jejunum1_l + out_jejunum2_l + out_ileum1_l + out_ileum2_l + out_ileum3_l + out_cecum_l + out_colon_l
    gutloss_l <- k5_l * (duodenum_ent_l + jejunum1_ent_l + jejunum2_ent_l + ileum1_ent_l + ileum2_ent_l + ileum3_ent_l + cecum_ent_l + colon_ent_l)
    rox_l <- kmet_l * cv_liver_l
    rhyd_l <- vmax_liver_l * cv_liver_l / (km_liver_l * (1 + cv_liver_d / km_liver_d) + cv_liver_l)
    d/dt(central_l) <- qc * (cv_l - ca_l)
    d/dt(adipose_l) <- qfat * (ca_l - cv_fat_l)
    d/dt(gonads_l) <- qgonad * (ca_l - cv_gonad_l)
    d/dt(liver_l) <- (qliver - 8 * q_seg) * ca_l - qliver * cv_liver_l - rox_l - rhyd_l + rbs_l
    d/dt(brain_l) <- qbrain * (ca_l - cv_brain_l)
    d/dt(a_rapidly_perfused_l) <- qrich * (ca_l - cv_rich_l)
    d/dt(a_slowly_perfused_l) <- qslow * (ca_l - cv_slow_l)
    d/dt(heart_l) <- qheart * (ca_l - cv_heart_l)
    # ritalinic acid formed in the liver and from the gut loss (fraction f_hydrol_gut)
    d/dt(central_ra_l) <- (rhyd_l + f_hydrol_gut * gutloss_l) * mw_ra / mw - ku_ra_l * central_ra_l / vbody

    # ---- Outputs: mixed venous plasma concentrations (ug/L = ng/mL) ----
    Cc_d <- cv_d
    Cc_l <- cv_l
    Cc_rac <- Cc_d + Cc_l
    C_ra_d <- central_ra_d / vbody
    C_ra_l <- central_ra_l / vbody

    Cc_d ~ prop(propSd_d)
    Cc_l ~ prop(propSd_l)
    Cc_rac ~ prop(propSd_rac)
  })
}
