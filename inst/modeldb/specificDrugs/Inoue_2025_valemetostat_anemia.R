Inoue_2025_valemetostat_anemia <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-safety model for CTCAE grade ",
    ">= 3 anemia in adults with relapsed/refractory peripheral T-cell ",
    "lymphoma treated with the oral EZH1/EZH2 dual inhibitor valemetostat ",
    "(Inoue 2025, n = 251, 200 mg once daily). The probability of the ",
    "event is expit(a0 + <exposure term> + x'b1 + x'b2 * zE), where x is ",
    "the baseline covariate vector, b1 are covariate effects on the logit ",
    "intercept and b2 are covariate effects on the exposure slope ",
    "(supplementary Table S5). The exposure term is linear in the ",
    "standardized exposure zE = (CSSU_VALE - 17.7) / 13.7 ng/mL. ",
    "Continuous covariates enter centred and scaled; binary covariates ",
    "enter raw. There is no PK layer and no ODE: the individual exposure ",
    "metric is supplied as the CSSU_VALE data column, computed by the ",
    "companion population PK model Inoue_2025_valemetostat.R from each ",
    "patient actual dosing history. No between-subject random effect and ",
    "no residual error are estimated (Bernoulli likelihood). One of seven ",
    "companion exposure-response models in the Inoue_2025_valemetostat_* ",
    "family."
  )
  reference <- paste(
    "Inoue H, Wang X, Garcia R, et al.",
    "Population Pharmacokinetics of Valemetostat and Exposure-Response Analyses of",
    "Efficacy and Safety in Patients with Relapsed/Refractory Peripheral T-Cell Lymphoma.",
    "J Clin Pharmacol. 2025;65(12):1699-1711. doi:10.1002/jcph.70100.",
    "Exposure-response parameters are from supplementary Table S5;",
    "the centring and scaling constants are from supplementary Figure S8.",
    sep = " "
  )
  vignette <- "Inoue_2025_valemetostat_ptcl"
  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CSSU_VALE covariate column)",
    concentration = "prob_anemia (probability of CTCAE grade >= 3 anemia, 0-1; also logit_anemia)"
  )

  covariateData <- list(
    CSSU_VALE = list(
      description        = "Unbound (free) valemetostat average plasma concentration over the on-treatment window up to the endpoint event (the paper Cavgtte). Supplied as data: this model has no PK layer, and the source analysis computed it from each patient post hoc parameters, actual dosing records and concomitant medications using the companion population PK model, so it reflects dose reductions and interruptions.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 17.7 ng/mL and scaled by 13.7 ng/mL inside model(). The centre is the Figure S8 reference subject exposure; the scale is free-fit from the three Figure S8 panels whose net exposure slope is well determined (thrombocytopenia, any Grade >= 3 TEAE and dose reduction), which agree on a centre of 17.4-18.2 and a scale of 13.5-13.8.",
      source_name        = "Cavgtte"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 69 and scaled by 11.85 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit over the five linear-exposure Figure S8 panels, which return a centre of 68.95 and a scale of 11.85; Table 1 ER-safety SD 11.9 agrees.",
      source_name        = "Age, years"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 71 and scaled by 15.57 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit over the Figure S8 panels with a non-negligible weight effect (centre 71.06, scale 15.57); Table 1 ER-safety SD is 16.1.",
      source_name        = "Weight, kg"
    ),
    HGB_BL = list(
      description        = "Baseline blood hemoglobin concentration -- the laboratory value corresponding to this model anemia endpoint.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 112.07 and scaled by 18.45 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit from the Figure S8A panel; Table 1 ER-safety SD 18.5 agrees.",
      source_name        = "Baseline hemoglobin, g/L"
    ),
    AAG = list(
      description        = "Baseline plasma alpha-1-acid glycoprotein concentration -- the binding protein that drives the saturable-binding component of the companion population PK model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 110 and scaled by 62.16 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit over the five linear-exposure Figure S8 panels (centre 109.55, scale 62.16); Table 1 ER-safety SD is 64.6.",
      source_name        = "Alpha-1-acid glycoprotein, mg/dL"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase concentration. Enters the model as log(LDH), standardized on the log scale.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as log(LDH), centred at 5.56 (log units) and scaled by 0.5803 inside model(); the centre corresponds to 259.8 untransformed. Scale source: free-fit over the five linear-exposure Figure S8 panels (centre 5.556, scale 0.5803); Table 1 reports a natural-scale SD only.",
      source_name        = "log baseline LDH, log U/L"
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
    HEPIMP = list(
      description        = "NCI-ODWG hepatic impairment indicator; 1 = mild or moderate impairment, 0 = normal hepatic function. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = "The source covariate is a two-level factor, normal versus mild-or-moderate impairment (Table S3). No patient with severe hepatic impairment was enrolled in either exposure-response analysis set (Table 1 reports only normal / mild / moderate), so the canonical HEPIMP mild-or-worse definition and the paper mild-or-moderate wording coincide exactly in this cohort. Do not substitute HEPIMP_MILD, which would silently exclude the moderate patients. Source row: NCI-ODWG hepatic function: Mild or Moderate impairment.",
      source_name        = "NCI-ODWG hepatic function: Mild or Moderate impairment"
    ),
    ECOG_GE1 = list(
      description        = "ECOG performance status indicator; 1 = ECOG PS 1 or greater, 0 = ECOG PS 0. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG PS 0)",
      notes              = "Source row: ECOG PS: 1+.",
      source_name        = "ECOG PS: 1+"
    ),
    TUMTP_PTCL = list(
      description        = "Patient type indicator; 1 = peripheral T-cell lymphoma, 0 = adult T-cell leukemia/lymphoma. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult T-cell leukemia/lymphoma; ATLL is the model reference patient type)",
      notes              = "The contrast is inverted relative to the companion population PK model Inoue_2025_valemetostat.R, which carries TUMTP_ATLL against a PTCL reference. Both directions are as their source estimated them and are not interchangeable once interaction terms are present. Source row: Patient type: PTCL.",
      source_name        = "Patient type: PTCL"
    ),
    TX_HCT = list(
      description        = "Prior hematopoietic stem cell transplant indicator; 1 = yes, 0 = no. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior transplant)",
      notes              = "Source row: Prior transplants: Yes.",
      source_name        = "Prior transplants: Yes"
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
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 251L,
    n_studies      = 3L,
    n_observations = "251 binary event records, one per patient",
    age_range      = "mean 65.6 years, SD 11.9 (Inoue 2025 Table 1, ER safety population)",
    weight_range   = "mean 71.0 kg, SD 16.1 (Inoue 2025 Table 1, ER safety population)",
    sex_female_pct = 35.9,
    race_ethnicity = c(`Asian, Japanese` = 30.3, `Asian, non-Japanese` = 8.4,
                       White = 45.4, Other = 15.9),
    disease_state  = paste0(
        "relapsed or refractory peripheral T-cell lymphoma (75.7%) or adult ",
        "T-cell leukemia/lymphoma (24.3%)"
      ),
    dose_range     = "valemetostat 150-300 mg orally once daily across the pooled trials; 200 mg in VALENTINE-PTCL01",
    regions        = "Japan (DS3201-A-J101, DS3201-A-J201) and multiregional (VALENTINE-PTCL01)",
    notes          = paste0(
        "Exposure-response safety analysis set: 251 patients pooled from ",
        "three trials -- 71 from DS3201-A-J101, 25 from DS3201-A-J201, 133 ",
        "from the VALENTINE-PTCL01 PTCL cohort and 22 from the VALENTINE- ",
        "PTCL01 ATLL cohort. The same 251-patient dataset and the same ",
        "standardization constants serve all six safety endpoints, so the six ",
        "models differ only in their estimates. ECOG PS 0 in 42.2%; one prior ",
        "regimen 27.9%, two 24.7%, more than two 47.4%; prior transplant ",
        "22.7%; NCI-ODWG hepatic function normal 74.9%, mild 23.5%, moderate ",
        "1.6%."
      )
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Inoue 2025 supplementary
    # Table S5. The table note states the reporting convention: 'main
    # effects and exposure effect exponentiated, with null value of 1;
    # interaction effects are on the original scale, with null value
    # of 0'. Each main effect and the exposure effect is therefore
    # written as log(printed odds ratio) so the published number stays
    # visible at the trace site, while each interaction effect is
    # written as the printed number itself.
    #
    # Reference subject (Figure S8 caption): 69 years, 70.9 kg, AAG
    # 110 mg/dL, log LDH 5.56, Cavgtte 18.1 ng/mL, male, Asian
    # Japanese, ECOG PS 0, PTCL, one prior regimen, no prior
    # transplant, normal hepatic function. Centring constants are that
    # subject covariate values; scaling constants are given per
    # parameter below and in covariateData.
    # ==================================================================

    # ----- Logit intercept (reference-subject event probability) -----
    logit_ref <- log(0.163 / (1 - 0.163))
    label("Logit of the CTCAE grade >= 3 anemia probability for the reference subject (unitless logit)")
    # Inoue 2025 Table S5, 'Population mean' 0.163 (0.0699, 0.314); Rhat 1.00, ESS-bulk 12,409

    # ----- Exposure effect -----
    e_cssu_logit <- log(1.60)
    label("Log-odds of CTCAE grade >= 3 anemia per 13.7 ng/mL increase in unbound valemetostat Cavgtte (unitless logit)")
    # Inoue 2025 Table S5, 'Unbound Cavg up to anemia (Cavgtte)' 1.60 (0.600, 4.22); Rhat 1.00, ESS-bulk 12,460

    # ----- Covariate effects on the logit intercept (b1) -----
    e_age_logit <- log(1.15)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for a 11.85 years increase in age (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Age, years' OR 1.15 (0.786, 1.73); Rhat 1.00, ESS-bulk 18,341
    e_wt_logit <- log(1.36)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for a 15.57 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Weight, kg' OR 1.36 (0.910, 2.04); Rhat 1.00, ESS-bulk 16,789
    e_hgb_logit <- log(0.447)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for a 18.45 g/L increase in baseline hemoglobin (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Baseline hemoglobin, g/L' OR 0.447 (0.301, 0.647); Rhat 1.00, ESS-bulk 18,094
    e_aag_logit <- log(1.04)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for a 62.16 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Alpha-1-acid glycoprotein, mg/dL' OR 1.04 (0.710, 1.52); Rhat 1.00, ESS-bulk 15,954
    e_ldh_logit <- log(1.43)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for a 0.5803 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'log baseline LDH, log U/L' OR 1.43 (1.01, 2.06); Rhat 1.00, ESS-bulk 17,184
    e_sexf_logit <- log(0.955)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Sex: Female' OR 0.955 (0.497, 1.83); Rhat 1.00, ESS-bulk 19,388
    e_race_asian_oth_logit <- log(1.48)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Race/Country: Asian non-Japanese' OR 1.48 (0.657, 3.81); Rhat 1.00, ESS-bulk 19,073
    e_race_white_logit <- log(0.892)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Race/Country: White' OR 0.892 (0.444, 1.76); Rhat 1.00, ESS-bulk 15,888
    e_race_other_logit <- log(0.777)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Race/Country: Other' OR 0.777 (0.355, 1.67); Rhat 1.00, ESS-bulk 16,991
    e_hepimp_logit <- log(0.896)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' OR 0.896 (0.443, 1.74); Rhat 1.00, ESS-bulk 22,579
    e_ecog_ge1_logit <- log(0.935)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'ECOG PS: 1+' OR 0.935 (0.491, 1.77); Rhat 1.00, ESS-bulk 17,815
    e_tumtp_ptcl_logit <- log(0.829)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for PTCL vs ATLL reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Patient type: PTCL' OR 0.829 (0.407, 1.70); Rhat 1.00, ESS-bulk 14,822
    e_tx_hct_logit <- log(1.20)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Prior transplants: Yes' OR 1.20 (0.608, 2.39); Rhat 1.00, ESS-bulk 18,495
    e_line_3l_logit <- log(0.776)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Number of prior regimens: 2' OR 0.776 (0.364, 1.58); Rhat 1.00, ESS-bulk 19,045
    e_line_4l_plus_logit <- log(1.26)
    label("Log-odds shift on the CTCAE grade >= 3 anemia logit for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S5 main effect, 'Number of prior regimens: >2' OR 1.26 (0.680, 2.38); Rhat 1.00, ESS-bulk 15,241

    # ----- Covariate effects on the exposure slope (b2) -----
    # These are printed on the logit scale already (null value 0), so
    # they are NOT wrapped in log().
    e_age_slope <- -0.150
    label("Shift in the unbound-exposure slope for a 11.85 years increase in age (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Age, years' -0.150 (-0.577, 0.270); Rhat 1.00, ESS-bulk 16,575
    e_wt_slope <- -0.222
    label("Shift in the unbound-exposure slope for a 15.57 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Weight, kg' -0.222 (-0.685, 0.238); Rhat 1.00, ESS-bulk 14,879
    e_hgb_slope <- 0.0394
    label("Shift in the unbound-exposure slope for a 18.45 g/L increase in baseline hemoglobin (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Baseline hemoglobin, g/L' 0.0394 (-0.417, 0.507); Rhat 1.00, ESS-bulk 21,031
    e_aag_slope <- -0.157
    label("Shift in the unbound-exposure slope for a 62.16 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'AAG, mg/dL' -0.157 (-0.533, 0.191); Rhat 1.00, ESS-bulk 15,975
    e_ldh_slope <- -0.115
    label("Shift in the unbound-exposure slope for a 0.5803 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'log baseline LDH, log U/L' -0.115 (-0.454, 0.232); Rhat 1.00, ESS-bulk 16,137
    e_sexf_slope <- 0.372
    label("Shift in the unbound-exposure slope for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Sex: Female' 0.372 (-0.280, 1.04); Rhat 1.00, ESS-bulk 18,531
    e_race_asian_oth_slope <- -0.361
    label("Shift in the unbound-exposure slope for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Race/Country: Asian non-Japanese' -0.361 (-1.45, 0.532); Rhat 1.00, ESS-bulk 15,729
    e_race_white_slope <- -0.326
    label("Shift in the unbound-exposure slope for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Race/Country: White' -0.326 (-1.04, 0.369); Rhat 1.00, ESS-bulk 16,577
    e_race_other_slope <- 0.435
    label("Shift in the unbound-exposure slope for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Race/Country: Other' 0.435 (-0.377, 1.32); Rhat 1.00, ESS-bulk 16,455
    e_hepimp_slope <- 0.358
    label("Shift in the unbound-exposure slope for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' 0.358 (-0.316, 1.07); Rhat 1.00, ESS-bulk 18,145
    e_ecog_ge1_slope <- 0.444
    label("Shift in the unbound-exposure slope for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'ECOG PS: 1+' 0.444 (-0.213, 1.13); Rhat 1.00, ESS-bulk 17,451
    e_tumtp_ptcl_slope <- -0.561
    label("Shift in the unbound-exposure slope for PTCL vs ATLL reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Patient type: PTCL' -0.561 (-1.37, 0.162); Rhat 1.00, ESS-bulk 14,304
    e_tx_hct_slope <- 0.462
    label("Shift in the unbound-exposure slope for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Prior transplants: Yes' 0.462 (-0.205, 1.19); Rhat 1.00, ESS-bulk 17,136
    e_line_3l_slope <- 0.122
    label("Shift in the unbound-exposure slope for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Number of prior regimens: 2' 0.122 (-0.631, 0.895); Rhat 1.00, ESS-bulk 18,740
    e_line_4l_plus_slope <- -0.216
    label("Shift in the unbound-exposure slope for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S5 interaction effect, 'Number of prior regimens: >2' -0.216 (-0.882, 0.437); Rhat 1.00, ESS-bulk 18,078

    # ----- No between-subject variability, no residual error -----
    # The source likelihood is Bernoulli on a single binary record per
    # patient, so there is no residual-error parameter and no random
    # effect to carry over. The placeholder additive SD below exists
    # only because rxode2 requires an endpoint definition; it is fixed
    # at a negligible value and is not from the source.
    addSd_prob_anemia <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value probability; the source likelihood is Bernoulli (no source residual)")
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    zexpo <- (CSSU_VALE - 17.7) / 13.7
    zage <- (AGE - 69) / 11.85
    zwt <- (WT - 71) / 15.57
    zhgb <- (HGB_BL - 112.07) / 18.45
    zaag <- (AAG - 110) / 62.16
    zldh <- (log(LDH) - 5.56) / 0.5803

    # ----- Covariate effect on the logit intercept (x' b1) -----
    cov_logit <- e_age_logit * zage +
                 e_wt_logit * zwt +
                 e_hgb_logit * zhgb +
                 e_aag_logit * zaag +
                 e_ldh_logit * zldh +
                 e_sexf_logit * SEXF +
                 e_race_asian_oth_logit * RACE_ASIAN_OTH +
                 e_race_white_logit * RACE_WHITE +
                 e_race_other_logit * RACE_OTHER +
                 e_hepimp_logit * HEPIMP +
                 e_ecog_ge1_logit * ECOG_GE1 +
                 e_tumtp_ptcl_logit * TUMTP_PTCL +
                 e_tx_hct_logit * TX_HCT +
                 e_line_3l_logit * LINE_3L +
                 e_line_4l_plus_logit * LINE_4L_PLUS

    # ----- Covariate effect on the exposure slope (x' b2) -----
    cov_slope <- e_age_slope * zage +
                 e_wt_slope * zwt +
                 e_hgb_slope * zhgb +
                 e_aag_slope * zaag +
                 e_ldh_slope * zldh +
                 e_sexf_slope * SEXF +
                 e_race_asian_oth_slope * RACE_ASIAN_OTH +
                 e_race_white_slope * RACE_WHITE +
                 e_race_other_slope * RACE_OTHER +
                 e_hepimp_slope * HEPIMP +
                 e_ecog_ge1_slope * ECOG_GE1 +
                 e_tumtp_ptcl_slope * TUMTP_PTCL +
                 e_tx_hct_slope * TX_HCT +
                 e_line_3l_slope * LINE_3L +
                 e_line_4l_plus_slope * LINE_4L_PLUS

    # ----- Linear exposure term on the standardized scale -----
    expo_effect <- e_cssu_logit * zexpo

    # ----- Linear predictor -----
    logit_anemia <- logit_ref + expo_effect + cov_logit + cov_slope * zexpo

    prob_anemia <- expit(logit_anemia)

    # ----- Observation -----
    prob_anemia ~ add(addSd_prob_anemia)
  })
}
