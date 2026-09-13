Inoue_2025_valemetostat_orr_bicr <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-efficacy model for objective ",
    "response rate (ORR, complete or partial response) assessed by ",
    "blinded independent central review (BICR) in adults with ",
    "relapsed/refractory peripheral T-cell lymphoma treated with the oral ",
    "EZH1/EZH2 dual inhibitor valemetostat (Inoue 2025, n = 119, 200 mg ",
    "once daily). The probability of the event is expit(a0 + <exposure ",
    "term> + x'b1 + x'b2 * zE), where x is the baseline covariate vector, ",
    "b1 are covariate effects on the logit intercept and b2 are covariate ",
    "effects on the exposure slope (supplementary Table S4). The exposure ",
    "term is linear in the standardized exposure zE = (CSSU_VALE - 13.9) ",
    "/ 12.2 ng/mL. Continuous covariates enter centred and scaled; binary ",
    "covariates enter raw. There is no PK layer and no ODE: the ",
    "individual exposure metric is supplied as the CSSU_VALE data column, ",
    "computed by the companion population PK model ",
    "Inoue_2025_valemetostat.R from each patient actual dosing history. ",
    "No between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). One of seven companion exposure-response ",
    "models in the Inoue_2025_valemetostat_* family."
  )
  reference <- paste(
    "Inoue H, Wang X, Garcia R, et al.",
    "Population Pharmacokinetics of Valemetostat and Exposure-Response Analyses of",
    "Efficacy and Safety in Patients with Relapsed/Refractory Peripheral T-Cell Lymphoma.",
    "J Clin Pharmacol. 2025;65(12):1699-1711. doi:10.1002/jcph.70100.",
    "Exposure-response parameters are from supplementary Table S4;",
    "the centring and scaling constants are from supplementary Figure S6.",
    sep = " "
  )
  vignette <- "Inoue_2025_valemetostat_ptcl"
  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CSSU_VALE covariate column)",
    concentration = "prob_orr_central (probability of objective response by blinded independent central review, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    CSSU_VALE = list(
      description        = "Unbound (free) valemetostat average plasma concentration over the on-treatment window up to the endpoint event (the paper Cavgtte). Supplied as data: this model has no PK layer, and the source analysis computed it from each patient post hoc parameters, actual dosing records and concomitant medications using the companion population PK model, so it reflects dose reductions and interruptions.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 13.9 ng/mL and scaled by 12.2 ng/mL inside model(). The centre is the Figure S6 reference subject exposure; the scale is printed in the Results narrative as the per-12.2 ng/mL increment.",
      source_name        = "Cavgtte"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase concentration. Enters the model as log(LDH), standardized on the log scale.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as log(LDH), centred at 5.5 (log units) and scaled by 0.579 inside model(); the centre corresponds to 244.7 untransformed. Scale source: solved from the Figure S6 panel; Table 1 reports a natural-scale SD only.",
      source_name        = "log baseline LDH, log U/L"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size, as the sum of the products of perpendicular diameters of target lesions. Enters the model as log(TUMSZ), standardized on the log scale.",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as log(TUMSZ), centred at 7.37 (log units) and scaled by 1.104 inside model(); the centre corresponds to 1588 untransformed. Scale source: solved from the Figure S6 panel; Table 1 reports a natural-scale SD only.",
      source_name        = "log baseline tumor size, mm2"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 69 and scaled by 12.6 inside model(). Centre is the Figure S6 reference subject value; scale source: Table 1 ER-efficacy SD 12.6 (the free-centre solve returns 12.37).",
      source_name        = "Age, years"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 74 and scaled by 13.4 inside model(). Centre is the Figure S6 reference subject value; scale source: Table 1 ER-efficacy SD 13.4 (solve 13.47).",
      source_name        = "Weight, kg"
    ),
    AAG = list(
      description        = "Baseline plasma alpha-1-acid glycoprotein concentration -- the binding protein that drives the saturable-binding component of the companion population PK model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 120 and scaled by 61.6 inside model(). Centre is the Figure S6 reference subject value; scale source: Table 1 ER-efficacy SD 61.6 (solve 62.16).",
      source_name        = "AAG, mg/dL"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Binary covariates are NOT centred or scaled.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the Figure S6 / S8 reference subject is male)",
      notes              = "Source row: Sex: Female.",
      source_name        = "Sex: Female"
    ),
    RACE_ASIAN_OTH = list(
      description        = "Asian non-Japanese indicator; 1 = Asian non-Japanese, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_WHITE and RACE_OTHER are also 0)",
      notes              = "Source row: Race/Country: Asian non-Japanese.",
      source_name        = "Race/Country: Asian non-Japanese"
    ),
    RACE_WHITE = list(
      description        = "White indicator; 1 = White, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_ASIAN_OTH and RACE_OTHER are also 0)",
      notes              = "Source row: Race/Country: White.",
      source_name        = "Race/Country: White"
    ),
    RACE_OTHER = list(
      description        = "Other race/country indicator; 1 = a race/country category outside Asian Japanese, Asian non-Japanese and White; 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_ASIAN_OTH and RACE_WHITE are also 0)",
      notes              = "Source row: Race/Country: Other.",
      source_name        = "Race/Country: Other"
    ),
    ECOG_GE1 = list(
      description        = "ECOG performance status indicator; 1 = ECOG PS 1 or greater, 0 = ECOG PS 0. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG PS 0)",
      notes              = "Source row: ECOG PS: 1+.",
      source_name        = "ECOG PS: 1+"
    ),
    LINE_3L = list(
      description        = "Third-line therapy indicator; 1 = exactly two prior anticancer regimens, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (paired with LINE_4L_PLUS; both 0 denotes one prior regimen, i.e. second-line therapy)",
      notes              = "The source column counts PRIOR REGIMENS, not lines: this indicator is 1 when the patient had exactly TWO prior regimens, which is third-line therapy. See the LINE_3L register entry for the off-by-one warning. Source row: Number of prior regimens: 2.",
      source_name        = "Number of prior regimens: 2"
    ),
    LINE_4L_PLUS = list(
      description        = "Fourth-line-or-later therapy indicator; 1 = more than two prior anticancer regimens, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (paired with LINE_3L; both 0 denotes one prior regimen, i.e. second-line therapy)",
      notes              = "The source column counts PRIOR REGIMENS, not lines: this indicator is 1 when the patient had MORE THAN TWO prior regimens, which is fourth-line therapy or later. Source row: Number of prior regimens: >2.",
      source_name        = "Number of prior regimens: >2"
    ),
    TUMTP_PTCL_NOS = list(
      description        = "PTCL not-otherwise-specified subtype indicator; 1 = PTCL-NOS, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (angioimmunoblastic T-cell lymphoma, when TUMTP_ALCL and TUMTP_PTCL_OTHER are also 0)",
      notes              = "Source row: PTCL subtypes: PTCL, NOS.",
      source_name        = "PTCL subtypes: PTCL, NOS"
    ),
    TUMTP_ALCL = list(
      description        = "Anaplastic large-cell lymphoma subtype indicator; 1 = ALCL, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (angioimmunoblastic T-cell lymphoma, when TUMTP_PTCL_NOS and TUMTP_PTCL_OTHER are also 0)",
      notes              = "Source row: PTCL subtypes: ALCL.",
      source_name        = "PTCL subtypes: ALCL"
    ),
    TUMTP_PTCL_OTHER = list(
      description        = "Residual PTCL subtype indicator; 1 = a PTCL subtype other than AITL, PTCL-NOS or ALCL; 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (angioimmunoblastic T-cell lymphoma, when TUMTP_PTCL_NOS and TUMTP_ALCL are also 0)",
      notes              = "Source row: PTCL subtypes: All other subtypes.",
      source_name        = "PTCL subtypes: All other subtypes"
    ),
    TX_HCT = list(
      description        = "Prior hematopoietic stem cell transplant indicator; 1 = yes, 0 = no. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior transplant)",
      notes              = "Source row: Prior transplants: Yes.",
      source_name        = "Prior transplants: Yes"
    ),
    HEPIMP = list(
      description        = "NCI-ODWG hepatic impairment indicator; 1 = mild or moderate impairment, 0 = normal hepatic function. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = "The source covariate is a two-level factor, normal versus mild-or-moderate impairment (Table S3). No patient with severe hepatic impairment was enrolled in either exposure-response analysis set (Table 1 reports only normal / mild / moderate), so the canonical HEPIMP mild-or-worse definition and the paper mild-or-moderate wording coincide exactly in this cohort. Do not substitute HEPIMP_MILD, which would silently exclude the moderate patients. Source row: NCI-ODWG hepatic function: Mild or Moderate impairment.",
      source_name        = "NCI-ODWG hepatic function: Mild or Moderate impairment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 119L,
    n_studies      = 1L,
    n_observations = "119 binary objective-response records, one per patient",
    age_range      = "mean 65.6 years, SD 12.6 (Inoue 2025 Table 1, ER efficacy population)",
    weight_range   = "mean 73.4 kg, SD 13.4 (Inoue 2025 Table 1, ER efficacy population)",
    sex_female_pct = 28.6,
    race_ethnicity = c(`Asian, Japanese` = 12.6, `Asian, non-Japanese` = 12.6,
                       White = 62.2, Other = 12.6),
    disease_state  = paste0(
        "relapsed or refractory peripheral T-cell lymphoma; subtypes AITL ",
        "35.3%, PTCL not otherwise specified 34.5%, ALCL 7.6%, all other ",
        "subtypes 22.7%"
      ),
    dose_range     = "valemetostat 200 mg orally once daily (a single dose level)",
    regions        = "multiregional (VALENTINE-PTCL01)",
    notes          = paste0(
        "Exposure-response efficacy analysis set: the 119 patients with PTCL ",
        "enrolled in the phase 2 VALENTINE-PTCL01 trial. Observed ORR by BICR ",
        "43.7%. Because only one dose level was studied, the exposure range ",
        "is narrow and the paper reports the exposure-response relationship ",
        "for ORR as flat (the exposure odds ratio credible interval spans 1). ",
        "ECOG PS 0 in 42.9%; one prior regimen 30.3%, two 25.2%, more than ",
        "two 44.5%; prior transplant 2.5%; NCI-ODWG hepatic function normal ",
        "80.7%, mild 16.8%, moderate 2.5%."
      )
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Inoue 2025 supplementary
    # Table S4. The table note states the reporting convention: 'main
    # effects and exposure effect exponentiated, with null value of 1;
    # interaction effects are on the original scale, with null value
    # of 0'. Each main effect and the exposure effect is therefore
    # written as log(printed odds ratio) so the published number stays
    # visible at the trace site, while each interaction effect is
    # written as the printed number itself.
    #
    # Reference subject (Figure S6 caption): 69 years, 74 kg, AAG 120
    # mg/dL, LDH 244.7 U/L, tumor size 1588 mm^2, Cavgtte 13.9 ng/mL,
    # male, Asian Japanese, ECOG PS 0, one prior regimen, AITL
    # subtype, no prior transplant, normal hepatic function. Centring
    # constants are that subject covariate values; scaling constants
    # are given per parameter below and in covariateData.
    # ==================================================================

    # ----- Logit intercept (reference-subject event probability) -----
    logit_ref <- log(0.597 / (1 - 0.597))
    label("Logit of the objective response by blinded independent central review probability for the reference subject (unitless logit)")
    # Inoue 2025 Table S4, 'Population mean' 0.597 (0.362, 0.802); Rhat 1.00, ESS-bulk 13,974

    # ----- Exposure effect -----
    e_cssu_logit <- log(1.10)
    label("Log-odds of objective response by blinded independent central review per 12.2 ng/mL increase in unbound valemetostat Cavgtte (unitless logit)")
    # Inoue 2025 Table S4, 'Unbound Cavg up to response (Cavgtte) by BICR' 1.10 (0.330, 4.05); Rhat 1.00, ESS-bulk 13,724

    # ----- Covariate effects on the logit intercept (b1) -----
    e_ldh_logit <- log(0.570)
    label("Log-odds shift on the objective response by blinded independent central review logit for a 0.579 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'log baseline LDH, log U/L' OR 0.570 (0.305, 0.960); Rhat 1.00, ESS-bulk 17,508
    e_tumsz_logit <- log(0.780)
    label("Log-odds shift on the objective response by blinded independent central review logit for a 1.104 log-unit increase in log baseline tumor size (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'log baseline tumor size, mm2' OR 0.780 (0.476, 1.25); Rhat 1.00, ESS-bulk 21,240
    e_age_logit <- log(1.32)
    label("Log-odds shift on the objective response by blinded independent central review logit for a 12.6 years increase in age (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Age, years' OR 1.32 (0.838, 2.12); Rhat 1.00, ESS-bulk 20,522
    e_wt_logit <- log(0.776)
    label("Log-odds shift on the objective response by blinded independent central review logit for a 13.4 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Weight, kg' OR 0.776 (0.484, 1.22); Rhat 1.00, ESS-bulk 19,972
    e_aag_logit <- log(0.692)
    label("Log-odds shift on the objective response by blinded independent central review logit for a 61.6 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'AAG, mg/dL' OR 0.692 (0.409, 1.14); Rhat 1.00, ESS-bulk 20,364
    e_sexf_logit <- log(1.61)
    label("Log-odds shift on the objective response by blinded independent central review logit for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Sex: Female' OR 1.61 (0.763, 4.00); Rhat 1.00, ESS-bulk 18,747
    e_race_asian_oth_logit <- log(0.799)
    label("Log-odds shift on the objective response by blinded independent central review logit for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Race/Country: Asian non-Japanese' OR 0.799 (0.315, 1.82); Rhat 1.00, ESS-bulk 18,597
    e_race_white_logit <- log(1.18)
    label("Log-odds shift on the objective response by blinded independent central review logit for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Race/Country: White' OR 1.18 (0.557, 2.51); Rhat 1.00, ESS-bulk 16,906
    e_race_other_logit <- log(0.724)
    label("Log-odds shift on the objective response by blinded independent central review logit for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Race/Country: Other' OR 0.724 (0.265, 1.70); Rhat 1.00, ESS-bulk 16,575
    e_ecog_ge1_logit <- log(0.575)
    label("Log-odds shift on the objective response by blinded independent central review logit for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'ECOG PS: 1+' OR 0.575 (0.258, 1.15); Rhat 1.00, ESS-bulk 19,687
    e_line_3l_logit <- log(1.51)
    label("Log-odds shift on the objective response by blinded independent central review logit for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Number of prior regimens: 2' OR 1.51 (0.711, 3.68); Rhat 1.00, ESS-bulk 16,258
    e_line_4l_plus_logit <- log(0.766)
    label("Log-odds shift on the objective response by blinded independent central review logit for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Number of prior regimens: >2' OR 0.766 (0.369, 1.55); Rhat 1.00, ESS-bulk 20,092
    e_tumtp_ptcl_nos_logit <- log(0.662)
    label("Log-odds shift on the objective response by blinded independent central review logit for PTCL-NOS subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'PTCL subtypes: PTCL, NOS' OR 0.662 (0.295, 1.37); Rhat 1.00, ESS-bulk 18,465
    e_tumtp_alcl_logit <- log(0.923)
    label("Log-odds shift on the objective response by blinded independent central review logit for ALCL subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'PTCL subtypes: ALCL' OR 0.923 (0.359, 2.25); Rhat 1.00, ESS-bulk 20,928
    e_tumtp_ptcl_other_logit <- log(0.879)
    label("Log-odds shift on the objective response by blinded independent central review logit for other PTCL subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'PTCL subtypes: All other subtypes' OR 0.879 (0.408, 1.81); Rhat 1.00, ESS-bulk 19,810
    e_tx_hct_logit <- log(0.945)
    label("Log-odds shift on the objective response by blinded independent central review logit for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'Prior transplants: Yes' OR 0.945 (0.297, 2.69); Rhat 1.00, ESS-bulk 19,110
    e_hepimp_logit <- log(0.812)
    label("Log-odds shift on the objective response by blinded independent central review logit for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S4 main effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' OR 0.812 (0.346, 1.81); Rhat 1.00, ESS-bulk 22,084

    # ----- Covariate effects on the exposure slope (b2) -----
    # These are printed on the logit scale already (null value 0), so
    # they are NOT wrapped in log().
    e_ldh_slope <- 0.0833
    label("Shift in the unbound-exposure slope for a 0.579 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'log baseline LDH, log U/L' 0.0833 (-0.551, 0.719); Rhat 1.00, ESS-bulk 18,146
    e_tumsz_slope <- 0.251
    label("Shift in the unbound-exposure slope for a 1.104 log-unit increase in log baseline tumor size (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'log baseline tumor size, mm2' 0.251 (-0.310, 0.817); Rhat 1.00, ESS-bulk 16,860
    e_age_slope <- -0.248
    label("Shift in the unbound-exposure slope for a 12.6 years increase in age (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Age, years' -0.248 (-0.803, 0.314); Rhat 1.00, ESS-bulk 21,466
    e_wt_slope <- -0.0289
    label("Shift in the unbound-exposure slope for a 13.4 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Weight, kg' -0.0289 (-0.593, 0.552); Rhat 1.00, ESS-bulk 18,455
    e_aag_slope <- -0.301
    label("Shift in the unbound-exposure slope for a 61.6 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'AAG, mg/dL' -0.301 (-0.920, 0.232); Rhat 1.00, ESS-bulk 18,079
    e_sexf_slope <- -0.00958
    label("Shift in the unbound-exposure slope for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Sex: Female' -0.00958 (-0.844, 0.817); Rhat 1.00, ESS-bulk 21,711
    e_race_asian_oth_slope <- -0.127
    label("Shift in the unbound-exposure slope for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Race/Country: Asian, non-Japanese' -0.127 (-1.12, 0.787); Rhat 1.00, ESS-bulk 19,869
    e_race_white_slope <- -0.568
    label("Shift in the unbound-exposure slope for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Race/Country: White' -0.568 (-1.82, 0.259); Rhat 1.00, ESS-bulk 13,478
    e_race_other_slope <- 0.0535
    label("Shift in the unbound-exposure slope for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Race/Country: Other' 0.0535 (-0.938, 1.07); Rhat 1.00, ESS-bulk 25,006
    e_ecog_ge1_slope <- -0.170
    label("Shift in the unbound-exposure slope for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'ECOG PS: 1+' -0.170 (-0.996, 0.622); Rhat 1.00, ESS-bulk 20,436
    e_line_3l_slope <- 0.382
    label("Shift in the unbound-exposure slope for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Number of prior regimens: 2' 0.382 (-0.520, 1.47); Rhat 1.00, ESS-bulk 19,312
    e_line_4l_plus_slope <- -0.126
    label("Shift in the unbound-exposure slope for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Number of prior regimens: >2' -0.126 (-0.931, 0.656); Rhat 1.00, ESS-bulk 21,071
    e_tumtp_ptcl_nos_slope <- 0.164
    label("Shift in the unbound-exposure slope for PTCL-NOS subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'PTCL subtypes: PTCL, NOS' 0.164 (-0.649, 1.00); Rhat 1.00, ESS-bulk 21,842
    e_tumtp_alcl_slope <- 0.0157
    label("Shift in the unbound-exposure slope for ALCL subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'PTCL subtypes: ALCL' 0.0157 (-1.00, 1.06); Rhat 1.00, ESS-bulk 20,572
    e_tumtp_ptcl_other_slope <- -0.155
    label("Shift in the unbound-exposure slope for other PTCL subtype vs AITL reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'PTCL subtypes: All other subtypes' -0.155 (-1.03, 0.684); Rhat 1.00, ESS-bulk 24,627
    e_tx_hct_slope <- -0.0217
    label("Shift in the unbound-exposure slope for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'Prior transplants: Yes' -0.0217 (-1.08, 1.05); Rhat 1.00, ESS-bulk 23,557
    e_hepimp_slope <- -0.304
    label("Shift in the unbound-exposure slope for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S4 interaction effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' -0.304 (-1.22, 0.528); Rhat 1.00, ESS-bulk 19,940

    # ----- No between-subject variability, no residual error -----
    # The source likelihood is Bernoulli on a single binary record per
    # patient, so there is no residual-error parameter and no random
    # effect to carry over. The placeholder additive SD below exists
    # only because rxode2 requires an endpoint definition; it is fixed
    # at a negligible value and is not from the source.
    addSd_prob_orr_central <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value probability; the source likelihood is Bernoulli (no source residual)")
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    zexpo <- (CSSU_VALE - 13.9) / 12.2
    zldh <- (log(LDH) - 5.5) / 0.579
    ztumsz <- (log(TUMSZ) - 7.37) / 1.104
    zage <- (AGE - 69) / 12.6
    zwt <- (WT - 74) / 13.4
    zaag <- (AAG - 120) / 61.6

    # ----- Covariate effect on the logit intercept (x' b1) -----
    cov_logit <- e_ldh_logit * zldh +
                 e_tumsz_logit * ztumsz +
                 e_age_logit * zage +
                 e_wt_logit * zwt +
                 e_aag_logit * zaag +
                 e_sexf_logit * SEXF +
                 e_race_asian_oth_logit * RACE_ASIAN_OTH +
                 e_race_white_logit * RACE_WHITE +
                 e_race_other_logit * RACE_OTHER +
                 e_ecog_ge1_logit * ECOG_GE1 +
                 e_line_3l_logit * LINE_3L +
                 e_line_4l_plus_logit * LINE_4L_PLUS +
                 e_tumtp_ptcl_nos_logit * TUMTP_PTCL_NOS +
                 e_tumtp_alcl_logit * TUMTP_ALCL +
                 e_tumtp_ptcl_other_logit * TUMTP_PTCL_OTHER +
                 e_tx_hct_logit * TX_HCT +
                 e_hepimp_logit * HEPIMP

    # ----- Covariate effect on the exposure slope (x' b2) -----
    cov_slope <- e_ldh_slope * zldh +
                 e_tumsz_slope * ztumsz +
                 e_age_slope * zage +
                 e_wt_slope * zwt +
                 e_aag_slope * zaag +
                 e_sexf_slope * SEXF +
                 e_race_asian_oth_slope * RACE_ASIAN_OTH +
                 e_race_white_slope * RACE_WHITE +
                 e_race_other_slope * RACE_OTHER +
                 e_ecog_ge1_slope * ECOG_GE1 +
                 e_line_3l_slope * LINE_3L +
                 e_line_4l_plus_slope * LINE_4L_PLUS +
                 e_tumtp_ptcl_nos_slope * TUMTP_PTCL_NOS +
                 e_tumtp_alcl_slope * TUMTP_ALCL +
                 e_tumtp_ptcl_other_slope * TUMTP_PTCL_OTHER +
                 e_tx_hct_slope * TX_HCT +
                 e_hepimp_slope * HEPIMP

    # ----- Linear exposure term on the standardized scale -----
    expo_effect <- e_cssu_logit * zexpo

    # ----- Linear predictor -----
    logit_orr_central <- logit_ref + expo_effect + cov_logit + cov_slope * zexpo

    prob_orr_central <- expit(logit_orr_central)

    # ----- Observation -----
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
