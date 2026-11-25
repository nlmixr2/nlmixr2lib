Mohamed_2026_tacrolimus <- function() {
  description <- paste0(
    "Steady-state population PK model of tacrolimus clearance and oral ",
    "bioavailability in adult allogeneic hematopoietic cell transplant ",
    "recipients, with clinical, pharmacogenomic and gut-microbiome covariates ",
    "(Mohamed 2026; n = 44, 222 continuous-IV-infusion steady-state ",
    "concentrations and 436 oral steady-state troughs). The model is ",
    "algebraic -- Css = F * daily dose / (24 * CL), with the daily dose ",
    "supplied as the DOSE_TAC_MGD covariate and F = 1 for the continuous IV ",
    "infusion (ROUTE_ORAL = 0) -- because every observation was a steady-state ",
    "concentration and no volume or absorption parameter was estimated. CL is ",
    "57.5% higher in CYP3A5*1 carriers and 20.3% / 16.9% lower with ",
    "voriconazole / eGFR < 60 mL/min/1.73 m^2 (linear factors), with ",
    "exponential effects of three stool bacterial species (centred-log-ratio ",
    "abundances); logit(F) is shifted additively by corticosteroid and ",
    "laxative use and ten bacterial species. Log-normal IIV on CL, ",
    "logit-normal IIV on oral F, proportional residual error."
  )
  reference <- paste(
    "Mohamed ME, Cheng S, Staley C, Rashidi A, El Jurdi N, Holtan SG,",
    "Jacobson PA. Identification of Microbiome Associations with Tacrolimus",
    "Pharmacokinetics in Adult Hematopoietic Cell Transplantation Using",
    "Population Pharmacokinetic and Machine Learning. Pharm Res.",
    "2026;43(9):3125. doi:10.1007/s11095-026-04152-x"
  )
  vignette <- "Mohamed_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  clrNotes <- paste0(
    "Centred-log-ratio (CLR) transformed relative abundance of the species in ",
    "the stool sample collected within +/- 10 days of the first oral ",
    "steady-state tacrolimus trough (MetaPhlAn 4.1.1 profiles from Aviti ",
    "shotgun metagenomic sequencing). Mohamed 2026 Methods: species present ",
    "in < 30% of samples or with mean relative abundance < 0.1% were removed, ",
    "leaving 32 of 897 species; zeros were replaced by half the smallest ",
    "relative abundance and the 32 abundances were CLR-transformed. The ",
    "control-stream centring medians are reproduced (to within 0.006) from ",
    "the per-subject abundances in Supplementary file 3 only when the 32 ",
    "retained abundances are re-closed to proportions before zero ",
    "replacement, so the value depends on the whole 32-species panel and is ",
    "not portable to a different taxon panel or sequencing pipeline. One ",
    "sample per subject; time-fixed."
  )

  covariateData <- list(
    DOSE_TAC_MGD = list(
      description = "Total daily tacrolimus dose (oral) or 24 x continuous IV infusion rate",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "The dose-rate input of the steady-state equations (Mohamed 2026 ",
        "Eqs. 1-2). Oral: the total daily immediate-release dose (median ",
        "3 mg/day, range 0.5-15, given in two divided doses; Table I). IV: ",
        "the continuous infusion rate times 24 (median 86 ug/h = 2.064 ",
        "mg/day; Table I and Results). The control stream computes ",
        "CSS = F1 * RATE * 1e3 / CL with RATE in mg/h; DOSE_TAC_MGD / 24 is ",
        "that RATE."
      ),
      source_name = "RATE (x 24)"
    ),
    ROUTE_ORAL = list(
      description = "Oral immediate-release dosing (1) vs continuous IV infusion (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (continuous IV infusion; F fixed to 1)",
      notes = paste0(
        "Control-stream ROUTE column: ROUTE = 0 sets F1 = 1 (IV continuous ",
        "infusion steady-state concentration, Css), ROUTE = 1 applies the ",
        "logit-normal oral bioavailability (oral steady-state trough, ",
        "Ct,ss). ROUTE_ORAL = ROUTE. Per-observation record."
      ),
      source_name = "ROUTE"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser (at least one CYP3A5*1 allele, rs776746)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3)",
      notes = paste0(
        "Control stream column CYP3A5_3_2 is 1 for *3/*3 (the most common ",
        "level, factor 1) and 0 for *1 carriers (factor 1 + THETA(3)), so ",
        "CYP3A5_EXPR = 1 - CYP3A5_3_2. Table I: *3/*3 77.27%, *1/*3 18.18%, ",
        "*1/*1 4.55%. CYP3A5*6 and *7 were absent from the cohort."
      ),
      source_name = "CYP3A5_3_2"
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole at the tacrolimus concentration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole)",
      notes = "Time-varying per concentration; 296 of 658 concentrations (44.98%, Table II).",
      source_name = "VORICONAZOLE"
    ),
    RENALIMP = list(
      description = "eGFR < 60 mL/min/1.73 m^2 (2021 CKD-EPI race-free creatinine equation)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (eGFR >= 60 mL/min/1.73 m^2)",
      notes = paste0(
        "Threshold for this model: eGFR < 60 mL/min/1.73 m^2 by the 2021 ",
        "CKD-EPI race-free equation from serum creatinine (Mohamed 2026 ",
        "Methods; Table II footnote). Time-varying (laboratory value on the ",
        "day of the concentration, otherwise the closest one carried over). ",
        "95 of 658 concentrations (14.44%, Table II). Control stream column ",
        "eGFR_EPI_RF_grp (1 = < 60)."
      ),
      source_name = "eGFR_EPI_RF_grp"
    ),
    CONMED_STEROID = list(
      description = "Concomitant corticosteroid at the tacrolimus concentration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no corticosteroid)",
      notes = paste0(
        "Time-varying per concentration; 157 of 658 concentrations (23.86%, ",
        "Table II). The paper does not list the agents or routes. Acts on ",
        "logit(F) only, so it has no effect on IV records."
      ),
      source_name = "Corticosteroids"
    ),
    CONMED_LAXATIVE = list(
      description = "Concomitant laxative at the tacrolimus concentration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no laxative)",
      notes = paste0(
        "Time-varying per concentration; 66 of 658 concentrations (10.03%, ",
        "Table II). Table II footnote: sennosides/docusate sodium, ",
        "polyethylene glycol, sennosides, docusate sodium, psyllium. Tested ",
        "on logit(F) only (Methods)."
      ),
      source_name = "Laxatives"
    ),
    GUT_BACTEROIDES_UNIFORMIS_CLR = list(
      description = "CLR-transformed stool relative abundance of Bacteroides uniformis",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -0.7839; acts on CL and logit(F)."),
      source_name = "uniformis"
    ),
    GUT_BACTEROIDES_OVATUS_CLR = list(
      description = "CLR-transformed stool relative abundance of Bacteroides ovatus",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.4625; acts on CL and logit(F)."),
      source_name = "ovatus"
    ),
    GUT_LIMOSILACTOBACILLUS_FERMENTUM_CLR = list(
      description = "CLR-transformed stool relative abundance of Limosilactobacillus fermentum",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.6278; acts on CL and logit(F)."),
      source_name = "fermentum"
    ),
    GUT_RUTHENIBACTERIUM_LACTATIFORMANS_CLR = list(
      description = "CLR-transformed stool relative abundance of Ruthenibacterium lactatiformans",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -0.6759; acts on logit(F)."),
      source_name = "lactatiformans"
    ),
    GUT_VEILLONELLA_ATYPICA_CLR = list(
      description = "CLR-transformed stool relative abundance of Veillonella atypica",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.7002; acts on logit(F)."),
      source_name = "atypica"
    ),
    GUT_SELLIMONAS_INTESTINALIS_CLR = list(
      description = "CLR-transformed stool relative abundance of Sellimonas intestinalis",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.4625; acts on logit(F)."),
      source_name = "intestinalis"
    ),
    GUT_ENTEROCOCCUS_FAECIUM_CLR = list(
      description = "CLR-transformed stool relative abundance of Enterococcus faecium",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -0.3043; acts on logit(F)."),
      source_name = "faecium"
    ),
    GUT_EGGERTHELLA_LENTA_CLR = list(
      description = "CLR-transformed stool relative abundance of Eggerthella lenta",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.5057; acts on logit(F)."),
      source_name = "lenta"
    ),
    GUT_BLAUTIA_CAECIMURIS_CLR = list(
      description = "CLR-transformed stool relative abundance of Blautia caecimuris",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.7002; acts on logit(F)."),
      source_name = "caecimuris"
    ),
    GUT_STREPTOCOCCUS_SALIVARIUS_CLR = list(
      description = "CLR-transformed stool relative abundance of Streptococcus salivarius",
      units = "(log-ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(clrNotes, " Centred at the control-stream median -1.3441; acts on logit(F)."),
      source_name = "salivarius"
    )
  )

  covariatesDataExcluded <- list(
    AGE_GT65 = list(
      description = "Age > 65 years",
      units = "(binary)",
      type = "binary",
      notes = "Screened by SCM on CL and F (Methods); not retained. Table I: 11 of 44 (25%) > 65 years."
    ),
    WT = list(
      description = "Actual body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = "Screened as a power model (Methods); not retained. Table I median 86.7 kg (range 46.5-160.7)."
    ),
    IBW = list(
      description = "Ideal body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = "Screened as a power model in a separate SCM run (Methods); not retained. Table I median 69.6 kg (range 43.2-80)."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. Table I: 13 of 44 female (29.55%)."
    ),
    HCT_COND_RIC = list(
      description = "Reduced-intensity conditioning (vs myeloablative)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. Table I: RIC 29 (non-myeloablative n = 1 pooled with RIC), myeloablative 15."
    ),
    SNP_POR_RS1057868 = list(
      description = "POR*28 (rs1057868) genotype",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened (*1/*1 vs *1/*28 vs *28/*28); not retained. Table I: 18 / 24 / 2."
    ),
    CONMED_POSACONAZOLE = list(
      description = "Concomitant posaconazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 110 of 658 concentrations (Table II)."
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 133 of 658 concentrations (Table II)."
    ),
    CONMED_TMPSMX = list(
      description = "Concomitant sulfamethoxazole-trimethoprim",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 162 of 658 concentrations (Table II)."
    ),
    CONMED_CCB = list(
      description = "Concomitant dihydropyridine calcium channel blocker (amlodipine or nifedipine)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 276 of 658 concentrations (Table II)."
    ),
    ALB = list(
      description = "Serum albumin (screened as < 3.4 g/dL vs normal)",
      units = "g/dL",
      type = "continuous",
      notes = "Screened as a time-varying category; not retained. Low at 367 of 658 concentrations (Table II)."
    ),
    BUN = list(
      description = "Blood urea nitrogen (screened as high / normal / low)",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened as a time-varying category (> 30, 7-30, < 7 mg/dL); not retained (Table II)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    n_observations = 658L,
    age_range = "33 of 44 (75%) <= 65 years, 11 (25%) > 65 years (age not otherwise summarised)",
    weight_range = "46.5-160.7 kg (actual body weight)",
    weight_median = "86.7 kg",
    sex_female_pct = 29.55,
    race_ethnicity = c(Caucasian = 95.45, `Native American` = 4.55),
    disease_state = paste0(
      "Adult allogeneic hematopoietic cell transplant recipients (University ",
      "of Minnesota, MT2017-12 microbiome study) on continuous IV tacrolimus ",
      "from day -3 followed by oral immediate-release tacrolimus. Diagnoses: ",
      "AML 34.1%, MDS 29.5%, ALL 15.9%, CML 4.5%, CMML 4.5%, other 11.4%. ",
      "Excluded: liver transplant or total bilirubin > 3 x ULN at admission, ",
      "colostomy or ileostomy."
    ),
    dose_range = paste0(
      "IV continuous infusion started at 0.03 mg/kg/day (observed rate 30-138 ",
      "ug/h, median 86); oral immediate-release capsules 0.5-15 mg/day ",
      "(median 3 mg/day in two divided doses) from a median of day 16 after ",
      "HCT. Target troughs 12-15 ng/mL in week 1, then 5-10 ng/mL."
    ),
    co_medication = paste0(
      "Fractions of the 658 concentrations: voriconazole 45%, posaconazole ",
      "16.7%, fluconazole 20.2%, sulfamethoxazole-trimethoprim 24.6%, ",
      "corticosteroids 23.9%, dihydropyridine CCB 42%, laxatives 10% ",
      "(Table II)."
    ),
    genotypes = paste0(
      "CYP3A5 *3/*3 77.27%, *1/*3 18.18%, *1/*1 4.55%; CYP3A4*1B carriers ",
      "13.64%; CYP3A4*22 carriers 6.82%; POR*28 *1/*28 54.55%, *28/*28 4.55% ",
      "(Table I)."
    ),
    regions = "United States (University of Minnesota, single centre)",
    sampling_window = paste0(
      "Steady-state concentrations only (>= 48 h on unchanged daily doses): ",
      "continuous-IV-infusion Css (median 5 per subject) and oral troughs up ",
      "to 60 days after the first oral steady-state trough (median 11 per ",
      "subject). Whole-blood LC/MS assay."
    ),
    notes = paste0(
      "NONMEM 7.5 FOCE-I; SCM for clinical / PGx / co-medication covariates; ",
      "LASSO on the CL and logit(F) EBEs to pick microbiome species, which ",
      "were then added by full fixed-effects modelling. Uncertainty by SIR."
    )
  )

  ini({
    # Final estimates: Mohamed 2026 Supplementary file 1 NONMEM control
    # stream $THETA / $OMEGA / $SIGMA, which round to Table III.
    lcl <- log(6.9101); label("Clearance CL for a CYP3A5*3/*3 patient without voriconazole, eGFR >= 60 (L/h)") # $THETA T1 6.9101; Table III CL 6.91 L/hr
    logitfdepot <- 0.59240; label("Logit of typical oral bioavailability F (logit scale)") # $THETA T2 0.59240; Table III F 0.644 = expit(0.5924)

    # Categorical effects on CL: factor (1 + theta) when the covariate is 1.
    e_cyp3a5_expr_cl <- 0.57532; label("Fractional change in CL for CYP3A5*1 carriers (unitless)") # $THETA CLCYP3A5_3_21 0.57532; Table III 0.575
    e_conmed_voriconazole_cl <- -0.20317; label("Fractional change in CL with voriconazole (unitless)") # $THETA CLVORICONAZOLE1 -0.20317; Table III -0.203
    e_renalimp_cl <- -0.16902; label("Fractional change in CL with eGFR < 60 mL/min/1.73 m^2 (unitless)") # $THETA CLeGFR_EPI_RF_grp1 -0.16902; Table III -0.169

    # Categorical effects on logit(F): additive shift when the covariate is 1.
    e_conmed_steroid_fdepot <- -0.44743; label("Additive shift in logit(F) with corticosteroids (logit scale)") # $THETA F1Corticosteroids1 -0.44743; Table III -0.447
    e_conmed_laxative_fdepot <- 0.52704; label("Additive shift in logit(F) with laxatives (logit scale)") # $THETA F1Laxatives1 0.52704; Table III 0.527

    # Microbiome effects on CL: exp(theta * (CLR - median)) (Eq. 9).
    e_gut_bacteroides_uniformis_clr_cl <- 0.020259; label("Exponential CL coefficient per unit CLR of Bacteroides uniformis (1/log-ratio)") # $THETA T8 0.020259; Table III 0.0203
    e_gut_bacteroides_ovatus_clr_cl <- 0.011339; label("Exponential CL coefficient per unit CLR of Bacteroides ovatus (1/log-ratio)") # $THETA T9 0.011339; Table III 0.0113
    e_gut_limosilactobacillus_fermentum_clr_cl <- -0.034464; label("Exponential CL coefficient per unit CLR of Limosilactobacillus fermentum (1/log-ratio)") # $THETA T10 -0.034464; Table III -0.0345

    # Microbiome effects on logit(F): theta * (CLR - median) (Eq. 10).
    e_gut_ruthenibacterium_lactatiformans_clr_fdepot <- 0.11965; label("Logit(F) shift per unit CLR of Ruthenibacterium lactatiformans (logit per log-ratio)") # $THETA T11 0.11965; Table III 0.12
    e_gut_veillonella_atypica_clr_fdepot <- 0.050486; label("Logit(F) shift per unit CLR of Veillonella atypica (logit per log-ratio)") # $THETA T12 0.050486; Table III 0.0505
    e_gut_sellimonas_intestinalis_clr_fdepot <- 0.050774; label("Logit(F) shift per unit CLR of Sellimonas intestinalis (logit per log-ratio)") # $THETA T13 0.050774; Table III 0.0508
    e_gut_enterococcus_faecium_clr_fdepot <- -0.060493; label("Logit(F) shift per unit CLR of Enterococcus faecium (logit per log-ratio)") # $THETA T14 -0.060493; Table III -0.0605
    e_gut_bacteroides_uniformis_clr_fdepot <- -0.067766; label("Logit(F) shift per unit CLR of Bacteroides uniformis (logit per log-ratio)") # $THETA T15 -0.067766; Table III -0.0678
    e_gut_bacteroides_ovatus_clr_fdepot <- -0.059150; label("Logit(F) shift per unit CLR of Bacteroides ovatus (logit per log-ratio)") # $THETA T16 -0.059150; Table III -0.0591
    e_gut_limosilactobacillus_fermentum_clr_fdepot <- 0.055190; label("Logit(F) shift per unit CLR of Limosilactobacillus fermentum (logit per log-ratio)") # $THETA T17 0.055190; Table III 0.0552
    e_gut_eggerthella_lenta_clr_fdepot <- 0.13937; label("Logit(F) shift per unit CLR of Eggerthella lenta (logit per log-ratio)") # $THETA T18 0.13937; Table III 0.139
    e_gut_blautia_caecimuris_clr_fdepot <- -0.11461; label("Logit(F) shift per unit CLR of Blautia caecimuris (logit per log-ratio)") # $THETA T19 -0.11461; Table III -0.115
    e_gut_streptococcus_salivarius_clr_fdepot <- -0.072378; label("Logit(F) shift per unit CLR of Streptococcus salivarius (logit per log-ratio)") # $THETA T20 -0.072378; Table III -0.0724

    # Diagonal IIV; the CL-F covariance was dropped once covariates were added.
    etalcl ~ 0.026110 # $OMEGA 0.026110; Table III IIV-CL 0.0261 (CV 16.3%)
    etalogitfdepot ~ 0.20796 # $OMEGA 0.20796; Table III IIV-logit(F) 0.208

    propSd <- 0.34969; label("Proportional residual error (fraction)") # $SIGMA 0.12228 variance, sqrt = 0.34969; Table III 0.122 (CV 35%)
  })

  model({
    # CL (control-stream $PRED): linear categorical factors times exponential
    # microbiome factors centred at the cohort medians of the CLR values.
    cl <- exp(lcl + etalcl) *
      (1 + e_cyp3a5_expr_cl * CYP3A5_EXPR) *
      (1 + e_conmed_voriconazole_cl * CONMED_VORICONAZOLE) *
      (1 + e_renalimp_cl * RENALIMP) *
      exp(e_gut_bacteroides_uniformis_clr_cl * (GUT_BACTEROIDES_UNIFORMIS_CLR - (-0.7839))) *
      exp(e_gut_bacteroides_ovatus_clr_cl * (GUT_BACTEROIDES_OVATUS_CLR - (-1.4625))) *
      exp(e_gut_limosilactobacillus_fermentum_clr_cl * (GUT_LIMOSILACTOBACILLUS_FERMENTUM_CLR - (-1.6278)))

    # logit(F): additive shifts (Eqs. 4, 8, 10); medians from the control stream.
    logitf <- logitfdepot + etalogitfdepot +
      e_conmed_steroid_fdepot * CONMED_STEROID +
      e_conmed_laxative_fdepot * CONMED_LAXATIVE +
      e_gut_ruthenibacterium_lactatiformans_clr_fdepot * (GUT_RUTHENIBACTERIUM_LACTATIFORMANS_CLR - (-0.6759)) +
      e_gut_veillonella_atypica_clr_fdepot * (GUT_VEILLONELLA_ATYPICA_CLR - (-1.7002)) +
      e_gut_sellimonas_intestinalis_clr_fdepot * (GUT_SELLIMONAS_INTESTINALIS_CLR - (-1.4625)) +
      e_gut_enterococcus_faecium_clr_fdepot * (GUT_ENTEROCOCCUS_FAECIUM_CLR - (-0.3043)) +
      e_gut_bacteroides_uniformis_clr_fdepot * (GUT_BACTEROIDES_UNIFORMIS_CLR - (-0.7839)) +
      e_gut_bacteroides_ovatus_clr_fdepot * (GUT_BACTEROIDES_OVATUS_CLR - (-1.4625)) +
      e_gut_limosilactobacillus_fermentum_clr_fdepot * (GUT_LIMOSILACTOBACILLUS_FERMENTUM_CLR - (-1.6278)) +
      e_gut_eggerthella_lenta_clr_fdepot * (GUT_EGGERTHELLA_LENTA_CLR - (-1.5057)) +
      e_gut_blautia_caecimuris_clr_fdepot * (GUT_BLAUTIA_CAECIMURIS_CLR - (-1.7002)) +
      e_gut_streptococcus_salivarius_clr_fdepot * (GUT_STREPTOCOCCUS_SALIVARIUS_CLR - (-1.3441))

    # Eq. 5 inverse logit for oral records; F = 1 for the IV infusion.
    fdepot <- ROUTE_ORAL * expit(logitf) + (1 - ROUTE_ORAL)

    # Eqs. 1-2: Css = F * daily dose / (24 * CL). mg/day / 24 / (L/h) = mg/L;
    # x 1000 gives ng/mL (control stream: CSS = F1 * RATE * 1e3 / CL).
    Cc <- fdepot * DOSE_TAC_MGD / 24 * 1000 / cl

    Cc ~ prop(propSd)
  })
}
