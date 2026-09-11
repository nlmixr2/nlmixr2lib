Inoue_2025_valemetostat <- function() {
  description <- paste0(
    "Three-compartment population pharmacokinetic model of TOTAL and ",
    "UNBOUND valemetostat, fitted simultaneously to both analytes, in ",
    "adults with relapsed/refractory non-Hodgkin lymphoma (adult T-cell ",
    "leukemia/lymphoma, peripheral T-cell lymphoma or other NHL) and in ",
    "healthy participants (Inoue 2025, n = 342 pooled across six trials: ",
    "J101, J201, VALENTINE-PTCL01, J107, J109 and U106). Valemetostat is ",
    "an oral EZH1/EZH2 dual inhibitor given 200 mg once daily. Absorption ",
    "is a sequential linked zero-order/first-order process: the dose is ",
    "released into the depot at a zero-order rate over duration D1 and ",
    "then absorbed first-order with KA. The central compartment holds ",
    "TOTAL valemetostat; the unbound concentration is recovered from it ",
    "by a saturable single-site binding sub-model with capacity BMAX ",
    "(the paper's RMAX) and dissociation constant KD, giving the closed ",
    "form Cu = (-(KD + BMAX - Ctot) + sqrt((KD + BMAX - Ctot)^2 + ",
    "4*KD*Ctot))/2. Elimination and BOTH inter-compartmental ",
    "distributions act on the UNBOUND concentration, so the model is ",
    "nonlinear in total drug even though unbound clearance is linear. ",
    "Alpha-1-acid glycoprotein drives the binding capacity ",
    "(BMAX ~ AAG^0.805) and additionally enters CL and F1 through a ",
    "single COMMON exponent (0.336), which is the device that makes ",
    "unbound exposure independent of AAG while total exposure rises with ",
    "it -- the paper's central pharmacological finding. Other covariates ",
    "(age, creatinine clearance, disease type, sex, race/country and ",
    "NCI-ODWG hepatic function on CL; P-gp and CYP3A inhibitor ",
    "comedication on F1) act on clearance or bioavailability, and body ",
    "weight enters with allometric exponents fixed at 0.75 and 1. A ",
    "study-specific factor (0.638) rescales unbound concentrations ",
    "measured by the J101 assay, which also carries its own residual ",
    "error. Fitted in NONMEM 7.5 by SAEM with interaction followed by ",
    "importance sampling. Companion exposure-response models for one ",
    "efficacy and six safety endpoints are the ",
    "Inoue_2025_valemetostat_* family, which consume the unbound ",
    "average concentration this model predicts."
  )
  reference <- paste(
    "Inoue H, Wang X, Garcia R, Reilly B, Tachibana M, Yoo Y, Lau Y, Chen Y.",
    "Population pharmacokinetics of valemetostat and exposure-response analyses",
    "of efficacy and safety in patients with relapsed/refractory peripheral",
    "T-cell lymphoma.",
    "J Clin Pharmacol. 2025;65(12):1699-1711. doi:10.1002/jcph.70100.",
    sep = " "
  )
  vignette <- "Inoue_2025_valemetostat_ptcl"

  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )
  # Unit convention. Inoue 2025 Table 2 prints KD and RMAX in nmol/L and the
  # paper nowhere states a valemetostat molar mass, so the binding sub-model
  # can only be written on a MOLAR scale without importing an off-source
  # constant. Doses are therefore in nmol and both outputs in nmol/L. The
  # volumes (L) and clearances (L/h) are unit-agnostic, so the unbound
  # steady-state average Cu = dose rate / CL is reproduced exactly in whatever
  # mass unit the dose is supplied in: 200 mg once daily over 24 h divided by
  # CL = 520 L/h gives 16.0 ng/mL, which brackets the paper's own reference
  # unbound Cavg values of 13.9 ng/mL (efficacy, Figure S6 caption) and
  # 18.1 ng/mL (safety, Figure S8 caption). See the vignette Errata.

  covariateData <- list(
    AAG = list(
      description        = "Baseline serum alpha-1-acid glycoprotein concentration. The sole binding partner in the central compartment and the only covariate the paper judged to materially affect TOTAL valemetostat exposure.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported in mg/dL by Inoue 2025 rather than in the register's canonical g/L; 100 mg/dL = 1 g/L. Reference value 100 mg/dL (Inoue 2025 Figure 2 caption reference individual). Enters in THREE places with only TWO estimated exponents: BMAX (exponent 0.805) and a single COMMON exponent (0.336) shared by CL and F1. Because CL and F1 carry the same exponent, unbound AUCss = dose * F1 / CL is algebraically independent of AAG while total exposure still rises with AAG through BMAX -- this is the paper's key mechanistic result (Results 'Final PPK Model Including Covariate Effects'; Discussion). Pooled PPK cohort mean (SD) 113 (62.3) mg/dL (Table 1); patients ran higher than non-patients (Figure S1).",
      source_name        = "AAG"
    ),
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric, with exponents FIXED at 0.75 for all clearance terms (CL, Q2, Q3) and 1 for all volume terms (V1, V2, V3) -- Inoue 2025 Results and Table 2 rows 'CL ~ WT' and 'V ~ WT', both marked FIXED. Reference weight 68.2 kg (Figure 2 caption reference individual). Pooled PPK cohort mean (SD) 69.3 (15.9) kg (Table 1).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on unbound CL with reference 65 years (Figure 2 caption reference individual). Exponent -0.205 with RSE 58.4% and a bootstrap 95% CI spanning zero (-0.553 to 0.0156, Table 2), i.e. retained in the full covariate model but not statistically resolved. Pooled PPK cohort mean (SD) 59.8 (17.2) years (Table 1).",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault calculated creatinine clearance at baseline.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on unbound CL with reference 83 mL/min (Figure 2 caption reference individual). Exponent -0.0107 with RSE 804% -- the least well-determined parameter in the model, consistent with valemetostat being eliminated hepatically. Pooled PPK cohort mean (SD) 89.2 (35.5) mL/min (Table 1).",
      source_name        = "CrCl"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the Figure 2 reference individual is male)",
      notes              = "Multiplicative effect on unbound CL, 1.06. Pooled PPK cohort 64.3% male (Table 1).",
      source_name        = "Female"
    ),
    TUMTP_ATLL = list(
      description        = "Adult T-cell leukemia/lymphoma indicator; 1 = ATLL, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (peripheral T-cell lymphoma, when DIS_BCELLNHL and DIS_HEALTHY are also 0)",
      notes              = "One of three indicators decomposing the paper's four-level 'population type' covariate (healthy participant / PTCL / ATLL / other NHL) against a PTCL reference; the Figure 2 caption reference individual is 'a male patient with PTCL'. Multiplicative effect on unbound CL, 0.828. Pooled PPK cohort 17.8% ATLL (Table 1).",
      source_name        = "ATLL"
    ),
    DIS_BCELLNHL = list(
      description        = "Other (non-ATLL, non-PTCL) non-Hodgkin lymphoma indicator; 1 = other NHL, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (peripheral T-cell lymphoma, when TUMTP_ATLL and DIS_HEALTHY are also 0)",
      notes              = "Inoue 2025 labels this stratum 'other NHL'. It is mapped onto the existing DIS_BCELLNHL canonical because the only study contributing non-ATLL non-PTCL lymphoma patients is J101, which the Methods describe as enrolling 'R/R non-Hodgkin lymphoma (NHL), including B-cell lymphomas, ATLL, and PTCL' -- so the residual NHL stratum is the B-cell lymphoma group. Multiplicative effect on unbound CL, 0.847. Pooled PPK cohort 5.6% (19 of 342, Table 1). See the vignette Errata: the paper's label is the histology-agnostic 'other NHL', so a future paper using an explicitly non-B-cell residual stratum should not reuse this mapping without checking.",
      source_name        = "OTHER NHL"
    ),
    DIS_HEALTHY = list(
      description        = "Non-patient indicator; 1 = healthy participant or non-cancer participant with hepatic impairment, 0 = lymphoma patient.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (peripheral T-cell lymphoma patient, when TUMTP_ATLL and DIS_BCELLNHL are also 0)",
      notes              = "Inoue 2025 calls this level 'non-patient'; it pools the healthy Japanese participants of the J107 DDI study and J109 food-effect study with the non-cancer hepatic-impairment participants of U106 (Figure S1 caption: 'The non-patient group includes both healthy participants and non-cancer patients with hepatic impairment'). Multiplicative effect on unbound CL, 0.942. Pooled PPK cohort 21.1% (72 of 342, Table 1).",
      source_name        = "NON-PATIENT"
    ),
    RACE_ASIAN_OTH = list(
      description        = "Asian non-Japanese indicator; 1 = Asian and enrolled outside Japan, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_WHITE and RACE_OTHER are also 0)",
      notes              = "The paper's covariate is a combined 'race/country enrolled' factor with four levels (Asian Japanese / Asian non-Japanese / White / Other) and an Asian-Japanese reference (Figure 2 caption reference individual is 'Asian Japanese'); the dominant Asian subgroup required by the RACE_ASIAN_OTH register entry is therefore Japanese. Multiplicative effect on unbound CL, 1.25 -- the largest non-AAG covariate effect in the model. Pooled PPK cohort 6.1% (Table 1).",
      source_name        = "ASIAN NON-J"
    ),
    RACE_WHITE = list(
      description        = "White race indicator; 1 = White, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_ASIAN_OTH and RACE_OTHER are also 0)",
      notes              = "Multiplicative effect on unbound CL, 0.967. Pooled PPK cohort 39.2% (Table 1).",
      source_name        = "WHITE"
    ),
    RACE_OTHER = list(
      description        = "Race category 'Other' indicator; 1 = other, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian Japanese, when RACE_ASIAN_OTH and RACE_WHITE are also 0)",
      notes              = "Multiplicative effect on unbound CL, 0.816. Pooled PPK cohort 12.9% (Table 1).",
      source_name        = "OTHER RACE"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator by NCI-ODWG criteria; 1 = mild, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, when HEPIMP_MOD is also 0)",
      notes              = "Multiplicative effect on unbound CL, 0.892. Pooled PPK cohort 21.1% mild (Table 1). Paired with HEPIMP_MOD; no severe stratum was enrolled.",
      source_name        = "HEPAT MILD"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator by NCI-ODWG criteria; 1 = moderate, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, when HEPIMP_MILD is also 0)",
      notes              = "Multiplicative effect on unbound CL, 0.752 -- a 33% higher unbound exposure. The paper quantifies the same effect from the post hoc simulations as a 29% (95% CI 1%-66%) increase in unbound AUCss. Pooled PPK cohort 3.5% moderate (12 of 342, Table 1).",
      source_name        = "HEPAT MOD"
    ),
    CONMED_PGP_INH = list(
      description        = "Concomitant P-glycoprotein inhibitor indicator; 1 = on a P-gp inhibitor and NOT on a CYP3A inhibitor, 0 = otherwise. Time-varying.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no interacting comedication, when CONMED_CYP3A4_INH is also 0)",
      notes              = "Inoue 2025 grouped its comedication categories 'due to small sample size' into none / P-gp inhibitor / CYP3A inhibitor with or without a P-gp inhibitor, so this indicator is the P-gp-ONLY arm and is superseded by CONMED_CYP3A4_INH when both are present -- the model() code enforces the precedence explicitly. Multiplicative effect on F1, 1.29. Treated as a TIME-VARYING covariate with immediate onset and immediate loss of effect (a limitation the Discussion flags). 14 of 342 participants (4.1%) took a P-gp inhibitor (Table 1).",
      source_name        = "DDI PGPi"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A inhibitor indicator; 1 = on a moderate or strong CYP3A inhibitor with or without a P-gp inhibitor, 0 = otherwise. Time-varying.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no interacting comedication, when CONMED_PGP_INH is also 0)",
      notes              = "Takes precedence over CONMED_PGP_INH per the paper's grouping ('CYP3Ai +/- P-gpi'). Multiplicative effect on F1, 1.23. 23 participants (6.7%) took a moderate CYP3A inhibitor, 1 (0.3%) a strong one and 1 (0.3%) a P-gp plus CYP3A inhibitor (Table 1). The Discussion notes the estimated magnitude is smaller than the dedicated DDI study found (4-fold with itraconazole, 1.6-fold with fluconazole) and attributes the gap to small numbers, missing comedication duration and the immediate-onset assumption.",
      source_name        = "DDI CYPi PGPi"
    ),
    STUDY_J101 = list(
      description        = "DS3201-A-J101 study indicator; 1 = the observation comes from the J101 phase 1 study, 0 = any of the other five studies.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (J201, VALENTINE-PTCL01, J107, J109 or U106)",
      notes              = "Applies ONLY to unbound valemetostat observations, and does two separate things: it rescales the predicted unbound concentration by the estimated assay factor 0.638 (Table 2 row 'ASSAY, DS3201-A-J101 unbound adjustment factor'), and it selects a different residual error (Sigma(3,3) = 0.327 rather than Sigma(2,2) = 0.404). Total valemetostat observations are unaffected. J101 contributed 71 of the 251 ER-safety patients (Table 1).",
      source_name        = "ASSAY"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "valemetostat", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "valemetostat", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "valemetostat", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "valemetostat", units = "nmol",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 342L,
    n_studies      = 6L,
    n_observations = "4635 total valemetostat concentrations from 342 participants and 3085 unbound valemetostat concentrations from 339 participants; 131 total (2.8%) and 43 unbound (1.4%) records were below the limit of quantification (Inoue 2025 Results)",
    age_range      = "mean (SD) 59.8 (17.2) years (Inoue 2025 Table 1)",
    weight_range   = "mean (SD) 69.3 (15.9) kg (Inoue 2025 Table 1)",
    sex_female_pct = 35.7,
    race_ethnicity = c(
      `Asian, Japanese`     = 41.8,
      `Asian, non-Japanese` = 6.1,
      White                 = 39.2,
      Other                 = 12.9
    ),
    disease_state  = "relapsed or refractory non-Hodgkin lymphoma -- peripheral T-cell lymphoma 55.6%, adult T-cell leukemia/lymphoma 17.8%, other NHL 5.6% -- pooled with 21.1% non-patients (healthy Japanese participants and non-cancer participants with hepatic impairment)",
    dose_range     = "valemetostat orally once daily; 200 mg is the approved and predominant dose (J201 and VALENTINE-PTCL01), with J101 contributing a multiple-ascending-dose escalation and J107/J109 single-dose healthy-participant data",
    regions        = "Japan and non-Japanese sites (VALENTINE-PTCL01 is multinational; J107 and J109 enrolled Japanese healthy participants; U106 enrolled non-Japanese participants with hepatic impairment)",
    hepatic_function = "normal 74.9%, mild NCI-ODWG impairment 21.1%, moderate 3.5% (Inoue 2025 Table 1)",
    renal_function = "creatinine clearance mean (SD) 89.2 (35.5) mL/min (Inoue 2025 Table 1)",
    co_medication  = "moderate CYP3A inhibitor 6.7%, strong CYP3A inhibitor 0.3%, P-gp inhibitor 4.1%, P-gp plus CYP3A inhibitor 0.3% (Inoue 2025 Table 1)",
    notes          = paste0(
      "Six trials: DS3201-A-J101 (NCT02732275, R/R NHL), J201 ",
      "(NCT04102150, R/R ATLL), VALENTINE-PTCL01 (NCT04703192, R/R PTCL ",
      "and ATLL), DS3201-A-J107 (jRCT2080225242, DDI in healthy Japanese ",
      "participants), DS3201-A-J109 (jRCT2071200043, food effect) and ",
      "DS3201-A-U106 (NCT04276662, hepatic impairment). Study-level ",
      "detail is in Inoue 2025 Table S1; demographics are in Table 1."
    )
  )

  ini({
    # ==================================================================
    # All structural and covariate values are the point estimates in
    # Inoue 2025 Table 2 ('Summary of PPK Parameter Estimates'), column
    # 'Estimate'. The Median and 95% CI columns of that table are
    # non-parametric bootstrap summaries (n = 462) and are NOT used here.
    #
    # Covariate reference values come from the Figure 2 caption, which
    # defines the reference individual as 'a male patient with PTCL,
    # Asian Japanese, no concomitant medication [CYP3Ai or P-gpi],
    # normal hepatic function, weighing 68.2 kg, 65 years of age, with a
    # AAG of 100 mg/dL, and a CrCl of 83 mL/min'.
    #
    # Scale conventions, both confirmed by the table's own printed
    # summary columns rather than assumed:
    #   * Continuous covariates enter as power terms (COV/ref)^theta --
    #     their estimates are centred on 0 (age -0.205, CrCl -0.0107)
    #     and sit alongside the two allometric WT exponents in the same
    #     block, which are unambiguously exponents.
    #   * Categorical covariates enter as multiplicative ratios --
    #     their estimates are centred on 1 (ATLL 0.828, female 1.06).
    # ==================================================================

    # ----- Disposition; CL, Q2 and Q3 act on the UNBOUND concentration -----
    lcl  <- log(520)   ; label("Apparent clearance of UNBOUND valemetostat (L/h)")            # Table 2, 'CL/F, L/h' = 520, RSE 10.8%
    lvc  <- log(42.5)  ; label("Apparent central volume of distribution (L)")                 # Table 2, 'V1/F, L' = 42.5, RSE 16.2%
    lq   <- log(40.0)  ; label("Apparent inter-compartmental clearance to peripheral1 (L/h)") # Table 2, 'Q2/F, L/h' = 40.0, RSE 26.9%
    lvp  <- log(3670)  ; label("Apparent first peripheral volume of distribution (L)")        # Table 2, 'V2/F, L' = 3670, RSE 14.7%
    lq2  <- log(230)   ; label("Apparent inter-compartmental clearance to peripheral2 (L/h)") # Table 2, 'Q3/F, L/h' = 230, RSE 9.46%
    lvp2 <- log(2950)  ; label("Apparent second peripheral volume of distribution (L)")       # Table 2, 'V3/F, L' = 2950, RSE 14.7%

    # ----- Sequential linked zero-order / first-order absorption -----
    lka <- log(0.374)  ; label("First-order absorption rate constant (1/h)")                  # Table 2, 'KA, 1/h' = 0.374, RSE 7.32%
    ld1 <- log(1.22)   ; label("Duration of the zero-order release into the depot (h)")       # Table 2, 'D1, h' = 1.22, RSE 11.8%

    # Relative bioavailability, FIXED to 1 as the estimation anchor; its
    # covariate effects (AAG and the comedication categories) are still
    # estimated and are applied on top of this value.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability of the depot input (fraction)") # Table 2, 'F1, Relative bioavailability' = 1.00, marked FIXED

    # ----- Saturable plasma-protein binding in the central compartment -----
    # The paper's RMAX is the canonical bmax (maximum binding capacity).
    lkd   <- log(221)  ; label("Equilibrium dissociation constant for valemetostat-AAG binding (nmol/L)") # Table 2, 'KD, nmol/L' = 221, RSE 11.4%
    lbmax <- log(8280) ; label("Total binding capacity at the reference AAG of 100 mg/dL (nmol/L)")       # Table 2, 'RMAX, nmol/L' = 8280, RSE 10.3%

    # ----- Study-specific unbound assay adjustment -----
    e_study_j101_cu <- 0.638 ; label("Multiplicative adjustment applied to UNBOUND concentrations measured by the DS3201-A-J101 assay (ratio)") # Table 2, 'ASSAY, DS3201-A-J101 unbound adjustment factor' = 0.638, RSE 6.76%

    # ----- Allometric body-weight exponents, both FIXED -----
    e_wt_cl <- fixed(0.750) ; label("Allometric exponent on CL, Q2 and Q3 (unitless)")  # Table 2, 'CL ~ WT' = 0.750, marked FIXED
    e_wt_vc  <- fixed(1.00)  ; label("Allometric exponent on V1, V2 and V3 (unitless)")  # Table 2, 'V ~ WT' = 1.00, marked FIXED

    # ----- Continuous covariate exponents -----
    e_aag_bmax <- 0.805   ; label("Exponent of AAG on the binding capacity bmax (unitless)")        # Table 2, 'RMAX ~ AAG' = 0.805, RSE 2.90%
    e_aag_cl   <- 0.336   ; label("COMMON exponent of AAG on both CL and F1 (unitless)")            # Table 2, 'CL/F ~ AAG' = 0.336 and 'F1 ~ AAG' = 0.336, both RSE 18.6% -- one estimated parameter reported on two rows
    e_age_cl   <- -0.205  ; label("Exponent of age on CL (unitless)")                               # Table 2, 'CL/F ~ AGE' = -0.205, RSE 58.4%
    e_crcl_cl  <- -0.0107 ; label("Exponent of creatinine clearance on CL (unitless)")              # Table 2, 'CL/F ~ CrCl' = -0.0107, RSE 804%

    # ----- Categorical covariate ratios on CL (reference: male, PTCL, Asian Japanese, normal hepatic function) -----
    e_tumtp_atll_cl     <- 0.828 ; label("CL ratio for adult T-cell leukemia/lymphoma versus PTCL (ratio)")      # Table 2, 'CL/F ~ ATLL' = 0.828, RSE 8.73%
    e_dis_bcellnhl_cl   <- 0.847 ; label("CL ratio for other non-Hodgkin lymphoma versus PTCL (ratio)")          # Table 2, 'CL/F ~ OTHER NHL' = 0.847, RSE 19.5%
    e_dis_healthy_cl    <- 0.942 ; label("CL ratio for non-patients versus PTCL patients (ratio)")               # Table 2, 'CL/F ~ NON-PATIENT' = 0.942, RSE 13.1%
    e_sexf_cl           <- 1.06  ; label("CL ratio for female versus male (ratio)")                              # Table 2, 'CL/F ~ FEMALE' = 1.06, RSE 5.87%
    e_race_asian_oth_cl <- 1.25  ; label("CL ratio for Asian non-Japanese versus Asian Japanese (ratio)")        # Table 2, 'CL/F ~ ASIAN NON-J' = 1.25, RSE 14.4%
    e_race_white_cl     <- 0.967 ; label("CL ratio for White versus Asian Japanese (ratio)")                     # Table 2, 'CL/F ~ WHITE' = 0.967, RSE 9.28%
    e_race_other_cl     <- 0.816 ; label("CL ratio for other race versus Asian Japanese (ratio)")                # Table 2, 'CL/F ~ OTHER RACE' = 0.816, RSE 8.81%
    e_hepimp_mild_cl    <- 0.892 ; label("CL ratio for mild NCI-ODWG hepatic impairment versus normal (ratio)")  # Table 2, 'CL/F ~ HEPAT MILD' = 0.892, RSE 6.20%
    e_hepimp_mod_cl     <- 0.752 ; label("CL ratio for moderate NCI-ODWG hepatic impairment versus normal (ratio)") # Table 2, 'CL/F ~ HEPAT MOD' = 0.752, RSE 18.8%

    # ----- Categorical covariate ratios on F1 (reference: no interacting comedication) -----
    e_conmed_pgp_inh_f1    <- 1.29 ; label("F1 ratio while taking a P-gp inhibitor without a CYP3A inhibitor (ratio)") # Table 2, 'F1 ~ DDI PGPi' = 1.29, RSE 30.3%
    e_conmed_cyp3a4_inh_f1 <- 1.23 ; label("F1 ratio while taking a CYP3A inhibitor with or without a P-gp inhibitor (ratio)") # Table 2, 'F1 ~ DDI CYPi PGPi' = 1.23, RSE 18.4%

    # ----- Inter-individual variability -----
    # Table 2 reports OMEGA VARIANCES on the log scale. The reporting
    # convention is settled by the table's own two derived columns, both
    # of which reproduce exactly:
    #   CV% = sqrt(exp(var) - 1)  -- 0.249 -> 53.2%, 2.49 -> 333%,
    #                               0.0563 -> 24.1%, 1.65 -> 205%,
    #                               0.0564 -> 24.1%, 0.485 -> 79.0%
    #   Corr = cov / sqrt(var1 * var2) -- 0.613/sqrt(0.249*2.49) = 0.778
    #                                     0.0629/sqrt(0.0563*1.65) = 0.207
    # Two block matrices, exactly as the paper describes: one for CL and
    # V1, another for KA and D1.
    etalcl + etalvc ~ c(0.249,
                        0.613, 2.49)     # Table 2, Omega(1,1) 0.249, Omega(2,1) 0.613 and Omega(2,2) 2.49; block for CL/F and V1/F
    etalka + etald1 ~ c(0.0563,
                        0.0629, 1.65)    # Table 2, Omega(3,3) 0.0563, Omega(4,3) 0.0629 and Omega(4,4) 1.65; block for KA and D1
    etalbmax ~ 0.0564                    # Table 2, Omega(5,5) 0.0564 [CV% = 24.1], RSE 10.8%; the paper's IIV-RMAX
    etalfdepot ~ 0.485                   # Table 2, Omega(6,6) 0.485 [CV% = 79.0], RSE 12.7%; the paper's IIV-F1

    # ----- Residual error -----
    # Table 2 reports three 'log-additive' SIGMA VARIANCES. NONMEM
    # log-additive error is proportional error in nlmixr2's linear
    # space, and the residual SD is sqrt(variance). Note that the
    # SIGMA rows use a DIFFERENT derived-column convention from the
    # OMEGA rows above: here the printed CV% is simply sqrt(variance)
    # (0.388 -> 62.3%, 0.404 -> 63.6%, 0.327 -> 57.2%), not
    # sqrt(exp(var) - 1), which would give 68.9%, 70.4% and 62.2%.
    propSd      <- sqrt(0.388) ; label("Proportional residual error for TOTAL valemetostat (fraction)")                       # Table 2, 'Log-additive - Total', Sigma(1,1) = 0.388 [CV% = 62.3], RSE 1.69%
    propSd_Cu      <- sqrt(0.404) ; label("Proportional residual error for UNBOUND valemetostat outside DS3201-A-J101 (fraction)") # Table 2, 'Log-additive - Unbound/non-DS3201-A-J101', Sigma(2,2) = 0.404 [CV% = 63.6], RSE 2.20%
    propSd_Cu_j101 <- sqrt(0.327) ; label("Proportional residual error for UNBOUND valemetostat within DS3201-A-J101 (fraction)")  # Table 2, 'Log-additive - Unbound/DS3201-A-J101', Sigma(3,3) = 0.327 [CV% = 57.2], RSE 14.7%
  })

  model({
    # ---- 1. Derived covariate multipliers -------------------------------
    # Continuous covariates are power terms on the Figure 2 reference
    # individual; categorical covariates are multiplicative ratios written
    # as (1 + (ratio - 1) * indicator) so that an indicator of 0 leaves the
    # reference value untouched.
    clCov <-
      (AGE / 65)^e_age_cl *
      (CRCL / 83)^e_crcl_cl *
      (1 + (e_tumtp_atll_cl - 1) * TUMTP_ATLL) *
      (1 + (e_dis_bcellnhl_cl - 1) * DIS_BCELLNHL) *
      (1 + (e_dis_healthy_cl - 1) * DIS_HEALTHY) *
      (1 + (e_sexf_cl - 1) * SEXF) *
      (1 + (e_race_asian_oth_cl - 1) * RACE_ASIAN_OTH) *
      (1 + (e_race_white_cl - 1) * RACE_WHITE) *
      (1 + (e_race_other_cl - 1) * RACE_OTHER) *
      (1 + (e_hepimp_mild_cl - 1) * HEPIMP_MILD) *
      (1 + (e_hepimp_mod_cl - 1) * HEPIMP_MOD)

    # The comedication categories are mutually exclusive: a participant on
    # both a CYP3A and a P-gp inhibitor belongs to the CYP3A arm, so the
    # P-gp-only term is switched off whenever CONMED_CYP3A4_INH is 1.
    f1Cov <-
      (1 + (e_conmed_pgp_inh_f1 - 1) * CONMED_PGP_INH * (1 - CONMED_CYP3A4_INH)) *
      (1 + (e_conmed_cyp3a4_inh_f1 - 1) * CONMED_CYP3A4_INH)

    # AAG enters CL and F1 through ONE shared exponent, which is what makes
    # unbound exposure (dose * F1 / CL) independent of AAG.
    aagCl <- (AAG / 100)^e_aag_cl

    # ---- 2. Individual parameters ---------------------------------------
    cl  <- exp(lcl + etalcl)  * (WT / 68.2)^e_wt_cl * clCov * aagCl
    vc  <- exp(lvc + etalvc)  * (WT / 68.2)^e_wt_vc
    q   <- exp(lq)            * (WT / 68.2)^e_wt_cl
    vp  <- exp(lvp)           * (WT / 68.2)^e_wt_vc
    q2  <- exp(lq2)           * (WT / 68.2)^e_wt_cl
    vp2 <- exp(lvp2)          * (WT / 68.2)^e_wt_vc

    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)
    f1 <- exp(lfdepot + etalfdepot) * f1Cov * aagCl

    bmax <- exp(lbmax + etalbmax) * (AAG / 100)^e_aag_bmax
    kd   <- exp(lkd)

    # ---- 3. Saturable binding: unbound concentration from total ---------
    # Ctot = Cu + bmax * Cu / (kd + Cu) rearranges to the quadratic
    # Cu^2 + (kd + bmax - Ctot) * Cu - kd * Ctot = 0, whose positive root is
    # taken below. At Cu << kd this collapses to Cu = Ctot / (1 + bmax/kd),
    # i.e. an unbound fraction of 1/(1 + 8280/221) = 2.60% at the reference
    # AAG -- the low-concentration free fraction implied by Table 2.
    ctot   <- central / vc
    bqterm <- kd + bmax - ctot
    cufree <- (-bqterm + sqrt(bqterm * bqterm + 4 * kd * ctot)) / 2

    # ---- 4. ODE system --------------------------------------------------
    # Elimination and BOTH distribution processes are driven by the unbound
    # concentration (Figure 1: CL leaves the "Unbound" species, and the
    # peripheral compartments equilibrate with it). This is what keeps the
    # terminal half-life near 11 h; driving them from the total
    # concentration instead would give roughly 14 days and could not reach
    # steady state under once-daily dosing.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * cufree -
                          q * (cufree - peripheral1 / vp) -
                          q2 * (cufree - peripheral2 / vp2)
    d/dt(peripheral1) <-  q * (cufree - peripheral1 / vp)
    d/dt(peripheral2) <-  q2 * (cufree - peripheral2 / vp2)

    # ---- 5. Dose input: zero-order release, then first-order absorption --
    # Requires rate = -2 in the event table so rxode2 applies dur().
    f(depot)   <- f1
    dur(depot) <- d1

    # ---- 6. Observations -------------------------------------------------
    # Unbound concentrations measured by the J101 assay are rescaled by the
    # estimated assay factor and carry their own residual error; total
    # concentrations are unaffected by study.
    Cc <- ctot
    Cu <- cufree * (1 - STUDY_J101 + e_study_j101_cu * STUDY_J101)

    sdCu <- propSd_Cu * (1 - STUDY_J101) + propSd_Cu_j101 * STUDY_J101

    Cc ~ prop(propSd)
    Cu ~ prop(sdCu)
  })
}
