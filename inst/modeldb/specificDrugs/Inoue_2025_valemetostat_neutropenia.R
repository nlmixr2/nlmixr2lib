Inoue_2025_valemetostat_neutropenia <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-safety model for CTCAE grade ",
    ">= 3 neutropenia in adults with relapsed/refractory peripheral ",
    "T-cell lymphoma treated with the oral EZH1/EZH2 dual inhibitor ",
    "valemetostat (Inoue 2025, n = 251, 200 mg once daily). The ",
    "probability of the event is expit(a0 + <exposure term> + x'b1 + x'b2 ",
    "* zE), where x is the baseline covariate vector, b1 are covariate ",
    "effects on the logit intercept and b2 are covariate effects on the ",
    "exposure slope (supplementary Table S6). This endpoint uses a ",
    "SATURATING exposure variable E = C / (ED50 + C) with ED50 = 39.7 ",
    "ng/mL in place of the linear standardized exposure, in both the main ",
    "exposure effect (coefficient Emax) and the interaction terms -- the ",
    "one endpoint of the seven whose exposure enters nonlinearly. ",
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
    "Exposure-response parameters are from supplementary Table S6;",
    "the centring and scaling constants are from supplementary Figure S8.",
    sep = " "
  )
  vignette <- "Inoue_2025_valemetostat_ptcl"
  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CSSU_VALE covariate column)",
    concentration = "prob_anc_decrease (probability of CTCAE grade >= 3 neutropenia, 0-1; also logit_anc_decrease)"
  )

  covariateData <- list(
    CSSU_VALE = list(
      description        = "Unbound (free) valemetostat average plasma concentration over the on-treatment window up to the endpoint event (the paper Cavgtte). Supplied as data: this model has no PK layer, and the source analysis computed it from each patient post hoc parameters, actual dosing records and concomitant medications using the companion population PK model, so it reflects dose reductions and interruptions.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "This model does NOT standardize the exposure. It enters through the saturating variable E = CSSU_VALE / (ED50 + CSSU_VALE) with ED50 = 39.7 ng/mL (Table S6), which is used both for the main exposure effect and as the multiplier of every interaction term. E = 0.313 at the reference exposure of 18.1 ng/mL, so unlike the other five safety models the interaction terms do NOT vanish at the reference subject -- that is what the Figure S8B binary rows show, and it is why this model needs no exposure centring or scaling constant.",
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
    NEUT = list(
      description        = "Baseline absolute neutrophil count -- the laboratory value corresponding to this model neutropenia endpoint.",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 3.643 and scaled by 3.104 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit from the Figure S8B panel under the saturating exposure variable; Table 1 ER-safety SD 3.29 agrees. This is the one covariate whose standardization constant is demonstrably NOT the Table 1 SD. The four printed Figure S8B percentile / probability pairs fall on a near-perfect straight line in the implied z, and that line returns a scale of 2.493 rather than the 3.29 printed in Table 1; using 3.29 does not reproduce the published probabilities. The solved value is therefore used and Table 1 is cited only where the two agree (hemoglobin, platelets, weight).",
      source_name        = "Baseline neutrophils, 109 cells/L"
    ),
    AAG = list(
      description        = "Baseline plasma alpha-1-acid glycoprotein concentration -- the binding protein that drives the saturable-binding component of the companion population PK model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 110 and scaled by 62.16 inside model(). Centre is the Figure S8 reference subject value; scale source: free-fit over the five linear-exposure Figure S8 panels (centre 109.55, scale 62.16); Table 1 ER-safety SD is 64.6.",
      source_name        = "AAG, mg/dL"
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
    # Table S6. The table note states the reporting convention: 'main
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
    # For this endpoint the intercept is the probability at ZERO
    # exposure, because the saturating exposure variable E = C / (ED50
    # + C) vanishes only at C = 0. That is why the Table S6 population
    # mean (0.241) differs from the Figure S8B ATLL reference row
    # (0.265): the gap is the exposure term evaluated at the reference
    # exposure, log(1.48) * 18.1 / (39.7 + 18.1) = 0.123 against an
    # observed logit gap of 0.127. For the other five safety endpoints
    # the exposure enters as a standardized term that vanishes AT the
    # reference exposure, so their population mean and their reference
    # row coincide exactly.
    #
    # Evidence that E also multiplies the INTERACTION terms, rather
    # than a standardized exposure that would be zero at the reference
    # subject: each binary covariate row of a forest panel satisfies
    # logit(p_k) - logit(p_ref) = a_k + b_k * zE_ref exactly, which
    # inverts to a scale-free estimate of zE_ref. The nine binary rows
    # of the anemia panel return zE_ref between -0.044 and +0.009 (the
    # interaction vanishes, as assumed for the linear models), while
    # the eight invertible binary rows of this panel return +0.35 to
    # +0.47 -- and C / (ED50 + C) = 0.313 at the reference exposure.
    # Substituting E for the linear term throughout drops the root-
    # mean-square residual over the 25 published rows of this panel
    # from 0.079 to 0.014 logit units with no fitted parameter. It
    # also reconciles the baseline neutrophil scale with Table 1:
    # solving under the old assumption returned 2.49 against a printed
    # SD of 3.29, whereas solving under E returns 3.10.
    logit_ref <- log(0.241 / (1 - 0.241))
    label("Logit of the CTCAE grade >= 3 neutropenia probability for the reference subject (unitless logit)")
    # Inoue 2025 Table S6, 'Population mean' 0.241 (0.0996, 0.459); Rhat 1.00, ESS-bulk 11,076

    # ----- Exposure effect -----
    emax_logit <- log(1.48)
    label("Maximum log-odds increase in grade >= 3 neutropenia at saturating unbound exposure (unitless logit)")
    # Inoue 2025 Table S6, 'Emax' 1.48 (0.173, 11.8); Rhat 1.00, ESS-bulk 13,556
    ed50 <- 39.7
    label("Unbound valemetostat concentration giving half the maximum neutropenia effect (ng/mL)")
    # Inoue 2025 Table S6, 'ED50' 39.7 (7.89, 84.8); Rhat 1.00, ESS-bulk 6172

    # ----- Covariate effects on the logit intercept (b1) -----
    e_age_logit <- log(1.26)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for a 11.85 years increase in age (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Age, years' OR 1.26 (0.802, 2.12); Rhat 1.00, ESS-bulk 10,914
    e_wt_logit <- log(0.989)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for a 15.57 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Weight, kg' OR 0.989 (0.626, 1.59); Rhat 1.00, ESS-bulk 13,106
    e_neut_logit <- log(0.538)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for a 3.104 10^9 cells/L increase in baseline neutrophils (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Baseline neutrophils, 109 cells/L' OR 0.538 (0.297, 1.07); Rhat 1.00, ESS-bulk 10,835
    e_aag_logit <- log(1.40)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for a 62.16 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'AAG, mg/dL' OR 1.40 (0.854, 2.33); Rhat 1.00, ESS-bulk 13,802
    e_ldh_logit <- log(1.07)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for a 0.5803 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'log baseline LDH, log U/L' OR 1.07 (0.683, 1.70); Rhat 1.00, ESS-bulk 13,879
    e_sexf_logit <- log(0.987)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Sex: Female' OR 0.987 (0.496, 1.87); Rhat 1.00, ESS-bulk 15,830
    e_race_asian_oth_logit <- log(1.25)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Race/Country: Asian non-Japanese' OR 1.25 (0.572, 2.83); Rhat 1.00, ESS-bulk 17,902
    e_race_white_logit <- log(0.572)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Race/Country: White' OR 0.572 (0.256, 1.18); Rhat 1.00, ESS-bulk 14,804
    e_race_other_logit <- log(0.788)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Race/Country: Other' OR 0.788 (0.354, 1.67); Rhat 1.00, ESS-bulk 16,901
    e_hepimp_logit <- log(1.36)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' OR 1.36 (0.696, 2.68); Rhat 1.00, ESS-bulk 19,188
    e_ecog_ge1_logit <- log(0.841)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'ECOG PS: 1+' OR 0.841 (0.452, 1.55); Rhat 1.00, ESS-bulk 17,412
    e_tumtp_ptcl_logit <- log(0.977)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for PTCL vs ATLL reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Patient type: PTCL' OR 0.977 (0.493, 1.94); Rhat 1.00, ESS-bulk 16,986
    e_tx_hct_logit <- log(1.02)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Prior transplants: Yes' OR 1.02 (0.516, 1.99); Rhat 1.00, ESS-bulk 18,901
    e_line_3l_logit <- log(1.09)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Number of prior regimens: 2' OR 1.09 (0.551, 2.21); Rhat 1.00, ESS-bulk 17,116
    e_line_4l_plus_logit <- log(1.22)
    label("Log-odds shift on the CTCAE grade >= 3 neutropenia logit for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S6 main effect, 'Number of prior regimens: >2' OR 1.22 (0.641, 2.29); Rhat 1.00, ESS-bulk 15,251

    # ----- Covariate effects on the exposure slope (b2) -----
    # These are printed on the logit scale already (null value 0), so
    # they are NOT wrapped in log().
    e_age_slope <- -0.298
    label("Shift in the unbound-exposure slope for a 11.85 years increase in age (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Age, years' -0.298 (-1.88, 0.638); Rhat 1.00, ESS-bulk 9548
    e_wt_slope <- -0.0697
    label("Shift in the unbound-exposure slope for a 15.57 kg increase in body weight (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Weight, kg' -0.0697 (-1.19, 0.962); Rhat 1.00, ESS-bulk 12,668
    e_neut_slope <- -0.485
    label("Shift in the unbound-exposure slope for a 3.104 10^9 cells/L increase in baseline neutrophils (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Baseline neutrophils, 109 cells/L' -0.485 (-2.98, 0.550); Rhat 1.00, ESS-bulk 6630
    e_aag_slope <- -0.154
    label("Shift in the unbound-exposure slope for a 62.16 mg/dL increase in AAG (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Alpha-1-acid glycoprotein, mg/dL' -0.154 (-1.35, 0.784); Rhat 1.00, ESS-bulk 12,377
    e_ldh_slope <- -0.166
    label("Shift in the unbound-exposure slope for a 0.5803 log-unit increase in log baseline LDH (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'log baseline LDH, log U/L' -0.166 (-1.35, 0.794); Rhat 1.00, ESS-bulk 12,839
    e_sexf_slope <- 0.255
    label("Shift in the unbound-exposure slope for female sex vs male reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Sex: Female' 0.255 (-0.800, 2.09); Rhat 1.00, ESS-bulk 11,304
    e_race_asian_oth_slope <- -0.0250
    label("Shift in the unbound-exposure slope for Asian non-Japanese vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Race/Country: Asian non-Japanese' -0.0250 (-1.44, 1.23); Rhat 1.00, ESS-bulk 17,167
    e_race_white_slope <- -0.215
    label("Shift in the unbound-exposure slope for White vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Race/Country: White' -0.215 (-2.02, 0.874); Rhat 1.00, ESS-bulk 12,102
    e_race_other_slope <- 0.00547
    label("Shift in the unbound-exposure slope for other race/country vs Asian Japanese reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Race/Country: Other' 0.00547 (-1.26, 1.29); Rhat 1.00, ESS-bulk 14,934
    e_hepimp_slope <- -0.103
    label("Shift in the unbound-exposure slope for mild or moderate hepatic impairment vs normal reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'NCI-ODWG hepatic function: Mild or Moderate impairment' -0.103 (-1.62, 1.00); Rhat 1.00, ESS-bulk 14,583
    e_ecog_ge1_slope <- -0.104
    label("Shift in the unbound-exposure slope for ECOG PS 1+ vs ECOG PS 0 reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'ECOG PS: 1+' -0.104 (-1.47, 1.00); Rhat 1.00, ESS-bulk 15,093
    e_tumtp_ptcl_slope <- -0.00689
    label("Shift in the unbound-exposure slope for PTCL vs ATLL reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Patient type: PTCL' -0.00689 (-1.23, 1.16); Rhat 1.00, ESS-bulk 18,022
    e_tx_hct_slope <- -0.0214
    label("Shift in the unbound-exposure slope for prior transplant vs none (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Prior transplants: Yes' -0.0214 (-1.31, 1.15); Rhat 1.00, ESS-bulk 15,897
    e_line_3l_slope <- -0.137
    label("Shift in the unbound-exposure slope for two prior regimens (third line) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Number of prior regimens: 2' -0.137 (-1.81, 0.953); Rhat 1.00, ESS-bulk 14,495
    e_line_4l_plus_slope <- 0.141
    label("Shift in the unbound-exposure slope for more than two prior regimens (fourth line or later) vs one prior regimen reference (unitless logit)")
    # Inoue 2025 Table S6 interaction effect, 'Number of prior regimens: >2' 0.141 (-0.928, 1.60); Rhat 1.00, ESS-bulk 15,039

    # ----- No between-subject variability, no residual error -----
    # The source likelihood is Bernoulli on a single binary record per
    # patient, so there is no residual-error parameter and no random
    # effect to carry over. The placeholder additive SD below exists
    # only because rxode2 requires an endpoint definition; it is fixed
    # at a negligible value and is not from the source.
    addSd_prob_anc_decrease <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value probability; the source likelihood is Bernoulli (no source residual)")
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    zage <- (AGE - 69) / 11.85
    zwt <- (WT - 71) / 15.57
    zneut <- (NEUT - 3.643) / 3.104
    zaag <- (AAG - 110) / 62.16
    zldh <- (log(LDH) - 5.56) / 0.5803

    # ----- Covariate effect on the logit intercept (x' b1) -----
    cov_logit <- e_age_logit * zage +
                 e_wt_logit * zwt +
                 e_neut_logit * zneut +
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
                 e_neut_slope * zneut +
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

    # ----- Saturating exposure variable -----
    # E rises from 0 at zero exposure to 1 at saturating exposure and
    # equals 0.313 at the reference exposure of 18.1 ng/mL. It
    # replaces the linear standardized exposure of the other six
    # models EVERYWHERE -- in the main exposure effect below and in
    # the interaction term of the linear predictor. See the ini()
    # intercept note for the evidence.
    expo_sat <- CSSU_VALE / (ed50 + CSSU_VALE)
    expo_effect <- emax_logit * expo_sat

    # ----- Linear predictor -----
    logit_anc_decrease <- logit_ref + expo_effect + cov_logit + cov_slope * expo_sat

    prob_anc_decrease <- expit(logit_anc_decrease)

    # ----- Observation -----
    prob_anc_decrease ~ add(addSd_prob_anc_decrease)
  })
}
