Marathe_2023_belzutifan <- function() {
  description <- paste(
    "Two-compartment population PK model for oral belzutifan (Welireg, a",
    "hypoxia-inducible factor 2 alpha inhibitor) in 239 pooled subjects",
    "across four phase I studies and one pivotal phase II study: 83 healthy",
    "participants, 74 patients with advanced renal cell carcinoma, 21 with",
    "other advanced solid tumors, and 61 with von Hippel-Lindau",
    "disease-associated RCC. First-order absorption with an absorption lag",
    "and linear elimination from the central compartment. Body weight scales",
    "clearances (shared exponent on CL/F and Q/F) and volumes (shared",
    "exponent on V2/F and V3/F); age scales CL/F and V2/F; fed state and the",
    "final market formulation each slow the absorption rate constant; and the",
    "polymorphic UGT2B17 and CYP2C19 metabolizer phenotypes shift CL/F, with",
    "UGT2B17 poor metabolizers additionally showing 11% higher relative",
    "bioavailability. Residual variability is split into healthy-participant",
    "and patient strata."
  )
  reference <- paste(
    "Marathe DD, Jauslin PM, Kleijn HJ, de Miranda Silva C, Chain A,",
    "Bateman T, Shaw PM, Abraham AK, Kauh EA, Liu Y, Perini RF,",
    "de Alwis DP, Jain L. Population pharmacokinetic analyses for",
    "belzutifan to inform dosing considerations and labeling.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1499-1510.",
    "doi:10.1002/psp4.13028",
    sep = " "
  )
  vignette <- "Marathe_2023_belzutifan"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power (allometric) scaling with reference weight 73.64 kg (the",
        "pooled-cohort median, Marathe 2023 Table S2 reports median 73.6 kg;",
        "the control stream in supplement 2 uses WT/73.64). A single",
        "estimated exponent 0.65 is shared by CL/F and Q/F, and a single",
        "estimated exponent 1.06 is shared by V2/F and V3/F -- see Marathe",
        "2023 Results ('CL/F and apparent inter-compartmental clearance",
        "(Q/F) were scaled by the same estimated exponent, as were V2/F and",
        "V3/F') and the Table 2 caption equations. Baseline (time-fixed) in",
        "the source analysis. Cohort range 42.1-165.8 kg."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate normalized as (AGE / 55), where 55 years is the",
        "pooled-cohort median age (Marathe 2023 Table S2). Applied to CL/F",
        "with exponent -0.36 and to V2/F with exponent -0.20 (Table 2).",
        "V3/F carries no age effect. Time-fixed at baseline. Cohort range",
        "19-84 years. The paper concluded the age effect is moderate and not",
        "clinically relevant (Table S4, Discussion)."
      ),
      source_name        = "AGE"
    ),
    FED = list(
      description        = "Fed-versus-fasted state at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Per-dose-record indicator. Linear-deviation effect on the",
        "absorption rate constant: KA = KApop * (1 + (-0.88) * FED), i.e.",
        "an 87.6% reduction in KA when dosed fed (2.40 -> 0.30 1/h; Marathe",
        "2023 Table 3). Food had no effect on AUCss (Discussion). Food was",
        "assessed in a crossover food-effect study (Study 2), so a subject",
        "can contribute both FED = 0 and FED = 1 records; Table S2 reports",
        "239 (100%) fasted and 15 (6.3%) fed records. The food effect on lag",
        "time reached significance in the SCM (Table S3) but could not be",
        "retained in the final model because the coefficient's confidence",
        "interval included zero and the run terminated with rounding errors",
        "(Results, 'Final population PK model')."
      ),
      source_name        = "FED"
    ),
    FORM_BELZ_FMF = list(
      description        = paste(
        "Belzutifan final market formulation (film-coated tablet) indicator;",
        "0 = oral compressed fit-for-purpose (FFP) tablet"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (FFP oral compressed fit-for-purpose tablet)",
      notes              = paste(
        "Per-dose-record indicator. Linear-deviation effect on the",
        "absorption rate constant: KA = KApop * (1 + (-0.47) *",
        "FORM_BELZ_FMF), i.e. a 47.4% reduction in KA for the film-coated",
        "final market formulation relative to the FFP compressed tablet",
        "(2.40 -> 1.26 1/h; Marathe 2023 Table 3). Formulation had no effect",
        "on AUCss (Discussion). Assessed in a crossover relative-",
        "bioavailability study (Study 6), so a subject can contribute both",
        "levels; Table S2 reports 190 (79.5%) FFP and 67 (28.0%) FMF",
        "records. The source control stream (supplement 2) codes FORM = 1",
        "for FFP (the most common level, reference), FORM = 2 for FMF, and",
        "FORM = -99 for missing, with missing mapped onto the FFP reference;",
        "the canonical column is re-oriented as FORM_BELZ_FMF = 1 for FMF",
        "and 0 for FFP or missing. Formulation was tested on lag time in the",
        "SCM but not retained (Table S3)."
      ),
      source_name        = "FORM"
    ),
    UGT2B17_EM = list(
      description        = "UGT2B17 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (UGT2B17 intermediate or poor metabolizer); paired with",
        "UGT2B17_PM = 0 the reference is the UGT2B17 intermediate-metabolizer",
        "stratum"
      ),
      notes              = paste(
        "Time-fixed per subject (germline UGT2B17 copy-number / deletion",
        "genotype-derived phenotype). Paired with UGT2B17_PM to encode the",
        "three-level EM / IM (reference) / PM phenotype with two binary",
        "indicators, following the CYP2D6_PM / CYP2D6_EM and CYP2C19_IM /",
        "CYP2C19_PM register precedents. Linear-deviation effect on CL/F:",
        "(1 + 0.39 * UGT2B17_EM), i.e. +39.1% CL/F for extensive",
        "metabolizers versus the intermediate-metabolizer reference (5.63 ->",
        "7.83 L/h; Marathe 2023 Table 2 and Table 3). Table S2 phenotype",
        "distribution: 46 poor (19.2%), 98 intermediate (41.0%), 84",
        "extensive (35.1%), 11 missing (4.6%); the source control stream",
        "maps the missing sentinel UGT2B17P = -99 onto the intermediate",
        "reference. All three UGT2B17 categories were supported on CL/F, but",
        "only two could be maintained for the effect on bioavailability, so",
        "EM and IM were merged there (Results, 'Final population PK model')",
        "-- hence UGT2B17_EM has no effect on F."
      ),
      source_name        = "UGT2B17P"
    ),
    UGT2B17_PM = list(
      description        = "UGT2B17 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (UGT2B17 intermediate or extensive metabolizer); for the CL/F",
        "effect, paired with UGT2B17_EM = 0 the reference is the",
        "intermediate-metabolizer stratum, whereas for the bioavailability",
        "effect the reference is the pooled intermediate + extensive stratum"
      ),
      notes              = paste(
        "Time-fixed per subject (homozygous UGT2B17 gene deletion). Two",
        "distinct effects. (1) Linear-deviation effect on CL/F:",
        "(1 + (-0.24) * UGT2B17_PM), i.e. -24.2% CL/F versus the",
        "intermediate-metabolizer reference (5.63 -> 4.27 L/h; Marathe 2023",
        "Table 2 and Table 3). (2) Linear-deviation effect on relative",
        "bioavailability: F = 1 * (1 + 0.11 * UGT2B17_PM), i.e. +11.0% F",
        "versus the pooled intermediate + extensive reference (Table 2",
        "caption equation and Table 3). The reference categories for the two",
        "effects differ because only two UGT2B17 levels could be maintained",
        "on F (EM and IM merged) while all three were supported on CL/F",
        "(Results, 'Final population PK model'). The paper attributes the F",
        "effect to reduced first-pass glucuronidation in the gut when",
        "UGT2B17 activity is absent (Discussion)."
      ),
      source_name        = "UGT2B17P"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (CYP2C19 non-poor metabolizer: the pooled intermediate,",
        "extensive, rapid, and ultrarapid strata)"
      ),
      notes              = paste(
        "Time-fixed per subject (germline genotype-derived phenotype).",
        "Linear-deviation effect on CL/F: (1 + (-0.36) * CYP2C19_PM), i.e.",
        "-36.0% CL/F versus the non-poor-metabolizer reference (5.63 -> 3.60",
        "L/h; Marathe 2023 Table 2 and Table 3). The source control stream",
        "branches on CYP2C19P .EQ. 1 (poor) versus .NE. 1 (everything else),",
        "so all non-PM phenotypes are pooled into the reference -- the paper",
        "states 'CYP2C19 categories were merged, as only the PM category had",
        "a significantly different CL/F compared to all other categories'",
        "(Results, 'Final population PK model'), and Figure 2 confirms 'the",
        "category extensive CYP2C19 metabolizer includes all non-poor",
        "metabolizers'. Table S2 phenotype distribution: 19 poor (7.9%), 65",
        "intermediate (27.2%), 96 extensive (40.2%), 42 rapid (17.6%), 6",
        "ultrarapid (2.5%), 11 missing (4.6%). CYP2C19 phenotype was tested",
        "on bioavailability in the SCM but not retained (Table S3)."
      ),
      source_name        = "CYP2C19P"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (patient: the pooled advanced-RCC, other-advanced-solid-tumor,",
        "and VHL-RCC cohorts from Studies 1 and 4)"
      ),
      notes              = paste(
        "Stratifies the proportional residual error only; it has no effect",
        "on any structural or random-effect PK parameter. Healthy",
        "participants use propSd = 0.26 (RES HV) and patients use",
        "propSd = 0.29 (RES PAT) per Marathe 2023 Table 2 (Residual error",
        "rows). The source control stream (supplement 2) sets",
        "PROP = THETA(7) by default and PROP = THETA(8) IF (STUDY .EQ. 1",
        ".OR. STUDY .EQ. 4) -- Studies 1 and 4 are the two patient studies",
        "(phase I dose-escalation/expansion in advanced RCC and other solid",
        "tumors, and the pivotal phase II single-arm study in VHL-RCC), so",
        "the STUDY switch is exactly a patient-versus-healthy-participant",
        "switch and is re-expressed here as DIS_HEALTHY = 1 for the phase I",
        "healthy-participant studies 2, 6, and 7. Cohort split (Table S2):",
        "83 healthy participants (34.7%) and 156 patients (65.3%: 74 RCC,",
        "21 other solid tumors, 61 VHL-RCC). Disease status was also tested",
        "as a covariate on CL/F and V2/F in the SCM and was NOT retained",
        "(Table S3), so no structural disease effect appears in this model."
      ),
      source_name        = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "A sex effect on V2/F was identified by the stepwise covariate",
        "model (Table S3) but was deliberately NOT retained in the final",
        "model: the studies were unbalanced with respect to sex (the",
        "single-dose healthy-participant Studies 2 and 7 enrolled no male",
        "participants and Study 6 only 3, whereas most patients in the",
        "multiple-dose studies were men), so the authors concluded the",
        "apparent sex effect may represent a population or study effect.",
        "Exploratory boxplots of etas by sex showed no influence, and",
        "repeating the SCM without sex identified no replacement covariate",
        "(Results, 'Covariate selection'). Cohort: 105 male (43.9%),",
        "Table S2."
      )
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F for the renal-impairment assessment and not",
        "retained (Table S3): 'No effects of renal and hepatic impairment",
        "were found in the covariate analysis based on available data.'",
        "Cohort median 77.5 (range 19.6-171.2) mL/min/1.73m^2, with 80",
        "normal, 104 mild, 52 moderate, 1 severe, 2 missing (Table S2).",
        "Severe renal impairment could not be evaluated (n = 1)."
      )
    ),
    HEP_NCI = list(
      description = paste(
        "NCI Organ Dysfunction Working Group hepatic-dysfunction category",
        "(composite of aspartate aminotransferase and total bilirubin)"
      ),
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened on CL/F for the hepatic-impairment assessment and not",
        "retained (Table S3). Cohort: 226 normal (94.6%), 12 mild (5.0%), 1",
        "moderate (0.4%) (Table S2). Moderate and severe hepatic impairment",
        "could not be evaluated because too few or no such participants were",
        "enrolled (Results, 'Final population PK model')."
      )
    ),
    RACE = list(
      description = "Race category",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened on CL/F and V2/F and not retained (Table S3). Cohort: 171",
        "White (71.5%), 38 Asian (15.9%), 23 Black (9.6%), 1 Pacific",
        "Islander (0.4%), 4 Multiple/Other (1.7%), 2 missing (0.8%)",
        "(Table S2)."
      )
    ),
    ETHNIC = list(
      description = "Hispanic or Latino ethnicity category",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened on CL/F and V2/F and not retained (Table S3). Cohort: 202",
        "not Hispanic (84.5%), 34 Hispanic (14.2%), 3 missing (1.3%)",
        "(Table S2)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Listed among the investigated covariates (Methods, 'Covariates')",
        "but explicitly never tested: 'BMI was not tested, because it is",
        "correlated with body weight, which was already included in the base",
        "model' (Results, 'Covariate selection')."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 239L,
    n_studies      = 5L,
    n_observations = 5291L,
    age_range      = "19-84 years (median 55)",
    age_median     = "55 years",
    weight_range   = "42.1-165.8 kg (median 73.6)",
    weight_median  = "73.6 kg",
    sex_female_pct = 56.1,
    race_ethnicity = c(
      White = 71.5, Black = 9.6, Asian = 15.9,
      `Pacific Islander` = 0.4, `Multiple/Other` = 1.7, Missing = 0.8
    ),
    ethnicity      = c(
      `Not Hispanic` = 84.5, Hispanic = 14.2, Missing = 1.3
    ),
    disease_state  = paste(
      "Pooled cohort across five studies: 83 healthy participants (34.7%)",
      "from the phase I food-effect (Study 2, NCT03445169), relative-",
      "bioavailability (Study 6, PT2977-104 / MK-6482-006), and",
      "Japanese-versus-White genotype/race (Study 7, PT2977-104 /",
      "MK-6482-007) studies; 74 patients with advanced renal cell carcinoma",
      "(31.0%) and 21 with other advanced solid tumors (8.8%) from the phase",
      "I dose-escalation/expansion Study 1 (NCT02974738); and 61 patients",
      "with von Hippel-Lindau disease-associated RCC (25.5%) from the",
      "pivotal phase II single-arm Study 4 (NCT03401788). Renal function:",
      "80 normal, 104 mild, 52 moderate, 1 severe impairment. Hepatic",
      "function (NCI-ODWG): 226 normal, 12 mild, 1 moderate."
    ),
    dose_range     = paste(
      "Oral. Study 1 dose escalation 20, 120, 160, and 240 mg q.d. and",
      "120 mg b.i.d., expansion 120 mg q.d.; Study 2 120 mg single dose",
      "fasted and fed; Study 4 120 mg q.d.; Study 6 120 mg single dose",
      "(FFP versus FMF crossover); Study 7 40, 120, and 200 mg single dose.",
      "Recommended clinical dosage is 120 mg orally once daily."
    ),
    regions        = paste(
      "Not stated per study; Study 7 enrolled Japanese and White",
      "participants in a parallel-group design"
    ),
    cohort_split   = "83 healthy participants (34.7%) + 156 patients (65.3%) = 239 total",
    egfr_median    = "77.5 mL/min/1.73m^2 (range 19.6-171.2)",
    pharmacogenetics = paste(
      "UGT2B17 phenotype: 46 poor (19.2%), 98 intermediate (41.0%), 84",
      "extensive (35.1%), 11 missing (4.6%). CYP2C19 phenotype: 19 poor",
      "(7.9%), 65 intermediate (27.2%), 96 extensive (40.2%), 42 rapid",
      "(17.6%), 6 ultrarapid (2.5%), 11 missing (4.6%). Ten subjects were",
      "UGT2B17/CYP2C19 dual poor metabolizers, all from Study 7; no dual PM",
      "subjects were present in the patient studies (Studies 1 and 4)."
    ),
    notes          = paste(
      "Baseline demographics from Marathe 2023 Table S2 (All Studies",
      "column, N = 239). The full PK analysis dataset contained 5291",
      "measurable observations; one Study 1 patient was excluded for having",
      "no postdose observations and seven observations were excluded (one",
      "missing sampling time, six predefined outliers), leaving 5284",
      "observations from 239 participants. The data cutoff for Studies 1 and",
      "4 (the two studies still ongoing at the time of analysis) was",
      "June 1, 2020. Model fit in NONMEM 7.3 with FOCE-I via ADVAN4 TRANS4."
    )
  )

  ini({
    # Structural PK parameters, Marathe 2023 Table 2 (final model). Typical
    # values correspond to the reference subject: body weight 73.64 kg, age
    # 55 years, UGT2B17 intermediate metabolizer, CYP2C19 non-poor
    # metabolizer, dosed fasted with the FFP formulation.
    lka   <- log(2.40);  label("First-order absorption rate constant for a fasted FFP dose (ka, 1/h)")  # Marathe 2023 Table 2 (KA = 2.40 /h)
    lcl   <- log(5.63);  label("Apparent clearance for the reference subject (CL/F, L/h)")              # Marathe 2023 Table 2 (CL/F = 5.63 L/h)
    lvc   <- log(85.4);  label("Apparent central volume of distribution for the reference subject (V2/F, L)") # Marathe 2023 Table 2 (V2/F = 85.4 L)
    lq    <- log(5.37);  label("Apparent inter-compartmental clearance for the reference subject (Q/F, L/h)") # Marathe 2023 Table 2 (Q/F = 5.37 L/h)
    lvp   <- log(30.38); label("Apparent peripheral volume of distribution for the reference subject (V3/F, L)") # Marathe 2023 Table 2 (V3/F = 30.38 L)
    ltlag <- log(0.16);  label("Absorption lag time (ALAG, h)")                                         # Marathe 2023 Table 2 (ALAG = 0.16 h)

    # Relative bioavailability is a structural anchor fixed to 1 for the
    # pooled UGT2B17 intermediate + extensive reference stratum: the source
    # control stream (supplement 2) writes F1 = 1 * F1UGT2B17P, so only the
    # UGT2B17 poor-metabolizer deviation below is estimated.
    lfdepot <- fixed(log(1)); label("Relative bioavailability for the UGT2B17 IM/EM reference (F, fraction)") # Marathe 2023 Table 2 caption (F = 1 * (1 + F-UGT2B17P)); supplement 2 control stream F1 = 1 * F1UGT2B17P

    # Covariate effects, Marathe 2023 Table 2. Body weight and age enter as
    # power (allometric-style) exponents on normalized covariates; the
    # categorical effects enter as linear deviations (1 + coefficient), which
    # is how the source control stream codes them.
    #
    # A single weight exponent is shared by CL/F and Q/F, and a single weight
    # exponent is shared by V2/F and V3/F (Marathe 2023 Results: "CL/F and
    # apparent inter-compartmental clearance (Q/F) were scaled by the same
    # estimated exponent, as were V2/F and V3/F"); the shared-exponent
    # covariate-effect naming convention encodes this as e_wt_cl_q and
    # e_wt_vc_vp. The estimates approximately recover the theoretical
    # small-molecule allometric values 0.75 and 1.0 (Discussion).
    e_wt_cl_q  <-  0.65; label("Power exponent of (WT / 73.64 kg) on CL/F and Q/F (unitless)")          # Marathe 2023 Table 2 (CL-WT = 0.65)
    e_wt_vc_vp <-  1.06; label("Power exponent of (WT / 73.64 kg) on V2/F and V3/F (unitless)")         # Marathe 2023 Table 2 (V-WT = 1.06)

    # Age enters only on CL/F and V2/F; V3/F carries no age effect (Table 2
    # caption equation "V3/F = V3/Fpop * (WT/73.64)^V-WT * eta3").
    e_age_cl   <- -0.36; label("Power exponent of (AGE / 55 years) on CL/F (unitless)")                 # Marathe 2023 Table 2 (CL-AGE = -0.36)
    e_age_vc   <- -0.20; label("Power exponent of (AGE / 55 years) on V2/F (unitless)")                 # Marathe 2023 Table 2 (V-AGE = -0.20)

    # Absorption-rate covariates: food and formulation both slow absorption.
    e_fed_ka           <- -0.88; label("Linear-deviation effect of fed state on ka (fraction)")         # Marathe 2023 Table 2 (KA-FED = -0.88; Table 3 gives -87.6%, 2.40 -> 0.30 1/h)
    e_form_belz_fmf_ka <- -0.47; label("Linear-deviation effect of the final market formulation on ka (fraction)") # Marathe 2023 Table 2 (KA-FORM = -0.47; Table 3 gives -47.4%, 2.40 -> 1.26 1/h)

    # Metabolizer-phenotype effects on CL/F. UGT2B17 is encoded with two
    # paired binary indicators so that the intermediate-metabolizer stratum
    # (both indicators zero) is the reference, matching the source control
    # stream's three-way branch on UGT2B17P. CYP2C19 pools every non-poor
    # phenotype into the reference.
    e_ugt2b17_em_cl <-  0.39; label("Linear-deviation effect of UGT2B17 extensive-metabolizer phenotype on CL/F (fraction)") # Marathe 2023 Table 2 (CL-UGT2B17 extensive metabolizers = 0.39; Table 3 gives +39.1%, 5.63 -> 7.83 L/h)
    e_ugt2b17_pm_cl <- -0.24; label("Linear-deviation effect of UGT2B17 poor-metabolizer phenotype on CL/F (fraction)")      # Marathe 2023 Table 2 (CL-UGT2B17 poor metabolizers = -0.24; Table 3 gives -24.2%, 5.63 -> 4.27 L/h)
    e_cyp2c19_pm_cl <- -0.36; label("Linear-deviation effect of CYP2C19 poor-metabolizer phenotype on CL/F (fraction)")      # Marathe 2023 Table 2 (CL-CYP2C19 poor metabolizers = -0.36; Table 3 gives -36.0%, 5.63 -> 3.60 L/h)

    # UGT2B17 poor metabolizers additionally show higher relative
    # bioavailability; the reference for THIS effect is the pooled
    # intermediate + extensive stratum because only two UGT2B17 levels could
    # be maintained on F (Results, "Final population PK model").
    e_ugt2b17_pm_fdepot <- 0.11; label("Linear-deviation effect of UGT2B17 poor-metabolizer phenotype on F (fraction)")      # Marathe 2023 Table 2 (F-UGT2B17 poor metabolizers = 0.11; Table 3 gives +11.0%, F 1 -> 1.11)

    # IIV. Marathe 2023 Table 2 reports the omega^2 variances directly (the
    # source control stream uses exponential IIV, EXP(ETA(n)), on CL, V2, V3,
    # and KA, with an OMEGA BLOCK(3) over CL-V2-V3 and a separate diagonal
    # OMEGA for KA). Q/F carries no IIV (control stream: Q = THETA(3) * CLWT,
    # no ETA).
    #
    # Variances (Table 2 "Random effects"):
    #   omega^2 CL/F = 0.15   -> CV = sqrt(exp(0.15) - 1)   = 40.3%
    #   omega^2 V2/F = 0.013  -> CV = sqrt(exp(0.013) - 1)  = 11.4%
    #   omega^2 V3/F = 0.19   -> CV = sqrt(exp(0.19) - 1)   = 45.8%
    #   omega^2 KA   = 1.15   -> CV = sqrt(exp(1.15) - 1)   = 138.1%
    #
    # The Table 2 note reports the eta CORRELATIONS, not covariances:
    #   CL/F-V2/F: 0.40, CL/F-V3/F: 0.54, V2/F-V3/F: 0.38
    # Covariance = rho * sqrt(var_1 * var_2):
    #   cov(CL, V2) = 0.40 * sqrt(0.15  * 0.013) = 0.40 * 0.0441588 = 0.0176635
    #   cov(CL, V3) = 0.54 * sqrt(0.15  * 0.19 ) = 0.54 * 0.1688194 = 0.0911625
    #   cov(V2, V3) = 0.38 * sqrt(0.013 * 0.19 ) = 0.38 * 0.0496991 = 0.0188857
    # nlmixr2 block order is lower-triangular row-major:
    #   c(var_cl, cov_cl_vc, var_vc, cov_cl_vp, cov_vc_vp, var_vp)
    etalcl + etalvc + etalvp ~ c(
      0.15,
      0.0176635, 0.013,
      0.0911625, 0.0188857, 0.19
    )                                                                                                  # Marathe 2023 Table 2 (IIV on CL/F, V2/F, V3/F) + Table 2 note (eta correlations 0.40 / 0.54 / 0.38)
    etalka ~ 1.15                                                                                      # Marathe 2023 Table 2 (IIV on KA = 1.15)

    # Residual error. The source uses a proportional error model where
    # $SIGMA is fixed to 1 and the proportional magnitude is carried by a
    # THETA: the control stream (supplement 2) sets W = IPRED * PROP and
    # Y = IPRED + W * ERR(1) with $SIGMA 1 FIX, so the reported RES values
    # ARE the proportional residual standard deviations (no sqrt needed).
    # PROP switches by study: THETA(7) by default and THETA(8) for
    # STUDY .EQ. 1 .OR. STUDY .EQ. 4 (the two patient studies).
    propSd_hv <- 0.26; label("Proportional residual SD for healthy participants (fraction)")            # Marathe 2023 Table 2 (RESHV = 0.26; EPS = 1 FIX)
    propSd_pt <- 0.29; label("Proportional residual SD for patients (fraction)")                        # Marathe 2023 Table 2 (RES PAT = 0.29; EPS = 1 FIX)
  })

  model({
    # Derived covariate multipliers. Reference values are the pooled-cohort
    # medians used by the source control stream: WT / 73.64 kg and
    # AGE / 55 years.
    #
    # UGT2B17 is a three-level phenotype encoded with two mutually exclusive
    # binary indicators, so the additive linear-deviation form below
    # reproduces the source control stream's three-way branch exactly:
    #   UGT2B17_EM = 1, UGT2B17_PM = 0 -> 1 + 0.39  (extensive)
    #   UGT2B17_EM = 0, UGT2B17_PM = 0 -> 1         (intermediate, reference)
    #   UGT2B17_EM = 0, UGT2B17_PM = 1 -> 1 - 0.24  (poor)
    cl_ugt2b17 <- 1 + e_ugt2b17_em_cl * UGT2B17_EM + e_ugt2b17_pm_cl * UGT2B17_PM
    cl_cyp2c19 <- 1 + e_cyp2c19_pm_cl * CYP2C19_PM
    ka_cov     <- (1 + e_fed_ka * FED) * (1 + e_form_belz_fmf_ka * FORM_BELZ_FMF)

    # Individual PK parameters. Marathe 2023 Table 2 caption equations:
    #   CL/F = CL/Fpop * (WT/73.64)^CL-WT * (1 + CL-UGT2B17P) *
    #            (1 + CL-CYP2C19P) * (AGE/55)^CL-AGE * eta1
    #   V2/F = V2/Fpop * (WT/73.64)^V-WT * (AGE/55)^V-AGE * eta2
    #   Q/F  = Q/Fpop  * (WT/73.64)^CL-WT
    #   V3/F = V3/Fpop * (WT/73.64)^V-WT * eta3
    #   KA   = KApop   * (1 + KA-FED) * (1 + KA-FORM) * eta4
    ka <- exp(lka + etalka) * ka_cov
    cl <- exp(lcl + etalcl) * (WT / 73.64)^e_wt_cl_q * cl_ugt2b17 * cl_cyp2c19 *
      (AGE / 55)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 73.64)^e_wt_vc_vp * (AGE / 55)^e_age_vc
    q  <- exp(lq) * (WT / 73.64)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 73.64)^e_wt_vc_vp

    # Micro-constants for the 2-compartment central-disposition model
    # (NONMEM ADVAN4 TRANS4: depot -> central (V2) <-> peripheral (V3),
    # elimination from central).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Relative bioavailability and absorption lag on the depot. The source
    # control stream applies F1 and ALAG1 to compartment 1, the depot.
    f(depot)    <- exp(lfdepot) * (1 + e_ugt2b17_pm_fdepot * UGT2B17_PM)
    alag(depot) <- exp(ltlag)

    # Plasma belzutifan concentration. Dose in mg and vc in L give
    # central / vc in mg/L = ug/mL; the source control stream scales with
    # S2 = V2/1000, i.e. DV = 1000 * A2 / V2, reporting ng/mL (Table 5
    # geometric-mean Cmax 1362.54 ng/mL at 120 mg q.d. steady state).
    Cc <- 1000 * central / vc

    # Population-stratified proportional residual error: healthy
    # participants use propSd_hv (RESHV = 0.26) and patients use propSd_pt
    # (RES PAT = 0.29) per Marathe 2023 Table 2.
    propSdEff <- propSd_hv * DIS_HEALTHY + propSd_pt * (1 - DIS_HEALTHY)
    Cc ~ prop(propSdEff)
  })
}
