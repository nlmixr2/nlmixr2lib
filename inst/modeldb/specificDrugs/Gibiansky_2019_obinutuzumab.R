Gibiansky_2019_obinutuzumab <- function() {
  description <- "Two-compartment population PK model of intravenous obinutuzumab (glycoengineered type II anti-CD20 IgG1) in adults with CD20+ B-cell malignancies (chronic lymphocytic leukaemia, indolent non-Hodgkin lymphoma, diffuse large B-cell lymphoma, mantle cell lymphoma) pooled from six trials including GADOLIN and GATHER; total clearance is a time-independent component cl_exp_inf plus a time-dependent component cl_exp_component * exp(-cl_exp_kdes * time) that decays as CD20 target is depleted. Covariates are body weight, albumin, age, sex, disease type and indolent lymphoma subtype, and baseline tumour size on the clearances; partner chemotherapy (CHOP, bendamustine, fludarabine + cyclophosphamide) and marginal zone lymphoma on the decay rate; weight and sex on the volumes. Log-scale residual error whose SD falls with concentration, with IIV on the residual SD (Gibiansky 2019, update of Gibiansky 2014)."
  reference <- "Gibiansky E, Gibiansky L, Buchheit V, Frey N, Brewster M, Fingerle-Rowson G, Jamois C. Pharmacokinetics, exposure, efficacy and safety of obinutuzumab in rituximab-refractory follicular lymphoma patients in the GADOLIN phase III study. Br J Clin Pharmacol. 2019;85(9):1935-1945. doi:10.1111/bcp.13974"
  vignette <- "Gibiansky_2019_obinutuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The IIV on the residual-error SD (Table S5 'omega^2 EPS') has no
  # structural typical-value partner; the residual SD is the
  # concentration-dependent expression built from sdL / sdH / sd50 in model().
  paper_specific_etas <- c("etaRUV")

  compartmentData <- list(
    central = list(analyte = "obinutuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "obinutuzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 75 kg, the median of the analysis data set (Table 1 footnote). Power exponents 0.639 on CL_inf, 0.378 on V1 and 1.08 on V2 (Table S5); no weight effect on CL_T or Q in the final model. Table 1 percentages at 52 and 115 kg (2.5th and 97.5th percentiles) are reproduced by (WT/75)^theta.",
      source_name = "BW"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 40 g/L (Table 1). Power exponent -0.682 on CL_inf (Table S5 'CL ALB'); CL_inf is 25.4% higher at 28.7 g/L (Table 1).",
      source_name = "BALB"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 65 years (Table 1). Power exponent -0.241 on CL_inf (Table S5 'CL AGE'); 13.8% higher CL_inf at 38 years (Table 1).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) is the paper's reference subject",
      notes = "The paper's reference subject is female (Results; Table 1 'Female' reference). The male factors in Table S5 (CL_T 1.45, CL_inf 1.18, V1 1.19) are applied as theta^(1 - SEXF). The source column coding is not printed; the male-vs-female direction is fixed by Table 1.",
      source_name = "SEX"
    ),
    TUMSZ = list(
      description = "Baseline tumour size (sum of the products of the diameters of the target lesions)",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column BSIZ (Table S4: overall mean 4900 (SD 5760) mm^2). Reference 3000 mm^2 (Table 1). Power exponents 0.351 on CL_T and 0.09 on CL_inf (Table S5 'CL T,BSIZ' and 'CL BSIZ'). The paper does not state how missing tumour size was handled; supply a value (for example the 3000 mm^2 reference) for every subject.",
      source_name = "BSIZ"
    ),
    TUMTP_CLL = list(
      description = "Chronic lymphocytic leukaemia disease-type indicator (1 = CLL, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FL or DLBCL is the paper's reference disease type)",
      notes = "Disease type was CLL, iNHL, DLBCL or MCL (Table S2). Factors 2.25 on CL_T and 1.47 on CL_inf (Table S5 'CL T,CLL' and 'CL CLL'). Small lymphocytic lymphoma is an iNHL subtype carried by TUMTP_SLL, not by this indicator.",
      source_name = "disease type"
    ),
    TUMTP_MCL = list(
      description = "Mantle cell lymphoma disease-type indicator (1 = MCL, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FL or DLBCL is the paper's reference disease type)",
      notes = "Factors 2.80 on CL_T and 2.07 on CL_inf (Table S5 'CL T,MCL' and 'CL MCL'); 20 MCL patients in the analysis.",
      source_name = "disease type"
    ),
    TUMTP_SLL = list(
      description = "Small lymphocytic lymphoma indicator (1 = SLL subtype of indolent NHL, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FL is the reference iNHL subtype)",
      notes = "iNHL subtype (FL, WM, MZL, SLL; Table S2). Factor 1.38 on CL_inf (Table S5 'CL SLL'); 22 SLL patients (Table S6).",
      source_name = "iNHL subtype"
    ),
    TUMTP_MZL = list(
      description = "Marginal zone lymphoma indicator (1 = MZL subtype of indolent NHL, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FL is the reference iNHL subtype)",
      notes = "Factor 0.299 on the CL_T decay rate k_des (Table S5 'K des,MZL'; CL_T half-life 21.1 days in the Results); 35 MZL patients (Table S6).",
      source_name = "iNHL subtype"
    ),
    CONMED_CHOP = list(
      description = "Co-administered CHOP chemotherapy (cyclophosphamide, doxorubicin, vincristine, prednisone) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CHOP)",
      notes = "Factor 0.309 on k_des (Table S5 'k des,CHOP'; Table 1 'Received CHOP' -69.1%). CHOP was given with obinutuzumab in GAUDI and GATHER (Table S1). The Results attribute the 20.4-day CL_T half-life (= the 0.309 factor) to 'CLL treated with chlorambucil', but no chlorambucil effect is tabulated; the model follows the tables (see the vignette).",
      source_name = "CHOP"
    ),
    CONMED_BENDAMUSTINE = list(
      description = "Co-administered bendamustine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no bendamustine)",
      notes = "Factor 0.591 on k_des (Table S5 'k des,bendamustine'); CL_T half-life 10.7 days in iNHL treated with bendamustine (Results). Bendamustine was the partner drug in GADOLIN and in a GAUDI arm (Table S1).",
      source_name = "bendamustine"
    ),
    CONMED_FC = list(
      description = "Co-administered fludarabine + cyclophosphamide (FC) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no FC)",
      notes = "Factor 2.64 on k_des (Table S5 'K des,fludarabine/cyclophosphamide'); CL_T half-life 2.4 days in FL treated with FC (Results). FC was given in GAUDI (Table S1).",
      source_name = "fludarabine + cyclophosphamide"
    )
  )

  covariatesDataExcluded <- list(
    ADA_POS = list(
      description = "Anti-drug antibody status",
      units = "(binary)",
      type = "binary",
      notes = "Tested on clearance (Table S2); antidrug antibodies were detected in one GADOLIN patient with no influence on PK (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 961L,
    n_studies = 6L,
    n_observations = 16301L,
    age_range = "mean 64.1 (SD 12.2) years; reference subject 65 years; 2.5th-97.5th percentiles 38-83 years (Table S4, Table 1)",
    weight_range = "mean 77.4 (SD 16.3) kg; median 75 kg; 2.5th-97.5th percentiles 52-115 kg (Table S4, Table 1)",
    albumin_range = "mean 40.3 (SD 5.18) g/L; median 40 g/L; 2.5th-97.5th percentiles 28.7-48.7 g/L",
    tumor_size_summary = "Baseline tumour size mean 4900 (SD 5760) mm^2; median 3000 mm^2; 2.5th-97.5th percentiles 304-22,400 mm^2",
    disease_state = "CD20+ B-cell malignancies: CLL 342, indolent NHL 469 (FL 406, MZL 35, SLL 22, WM 5), DLBCL 130, MCL 20 patients. GADOLIN (183 patients, rituximab-refractory iNHL, 145 FL) and GATHER (100 DLBCL) were added to the four studies of Gibiansky 2014.",
    dose_range = "Intravenous obinutuzumab 50-2000 mg; standard regimen 1000 mg on days 1, 8 and 15 of cycle 1 and day 1 of later cycles. GADOLIN: 1000 mg with bendamustine for six 28-day cycles, then 1000 mg every 2 months for up to 2 years.",
    regions = "International multicentre trials",
    studies = "GAUSS (BO21003, n = 105), GAUGUIN (BO20999, n = 131), GAUDI (BO21000, n = 134), CLL11 (BO21004, n = 308), GADOLIN (GAO4753g, n = 183), GATHER (GAO4915g, n = 100) (Table S1).",
    notes = "Sex distribution is not reported in the paper or supplement. Serum concentrations by ELISA, LLOQ 4.05 ng/mL. NONMEM 7.3, IMPMAP estimation."
  )

  ini({
    # Structural parameters (Table S5, final covariate model). Reference
    # subject: 65-year-old female with FL or DLBCL, 75 kg, albumin 40 g/L,
    # baseline tumour size 3000 mm^2, no partner chemotherapy (Results and
    # Table 1).
    lcl_exp_kdes <- log(0.11); label("Decay rate constant of the time-dependent clearance k_des (1/day)") # Table S5 k_des = 0.11 1/day (half-life 6.3 days, Results)
    lcl_exp_component <- log(0.154); label("Initial value of the time-dependent clearance CL_T (L/day)") # Table S5 CL_T = 0.154 L/day
    lcl_exp_inf <- log(0.0742); label("Time-independent clearance CL_inf (L/day)") # Table S5 CL_inf = 0.0742 L/day
    lvc <- log(2.72); label("Central volume of distribution V1 (L)") # Table S5 V1 = 2.72 L
    lvp <- log(1.23); label("Peripheral volume of distribution V2 (L)") # Table S5 V2 = 1.23 L
    lq <- log(1.32); label("Intercompartmental clearance Q (L/day)") # Table S5 Q = 1.32 L/day

    # Continuous covariates: power functions of covariate / reference
    # (Table 1 percentages reproduce (cov/ref)^theta).
    e_wt_cl_exp_inf <- 0.639; label("Power exponent of (WT/75) on CL_inf (unitless)") # Table S5 CL_inf,WT = 0.639
    e_wt_vc <- 0.378; label("Power exponent of (WT/75) on V1 (unitless)") # Table S5 V1,WT = 0.378
    e_wt_vp <- 1.08; label("Power exponent of (WT/75) on V2 (unitless)") # Table S5 V2,WT = 1.08
    e_alb_cl_exp_inf <- -0.682; label("Power exponent of (ALB/40) on CL_inf (unitless)") # Table S5 CL ALB = -0.682
    e_age_cl_exp_inf <- -0.241; label("Power exponent of (AGE/65) on CL_inf (unitless)") # Table S5 CL AGE = -0.241
    e_tumsz_cl_exp_inf <- 0.09; label("Power exponent of (TUMSZ/3000) on CL_inf (unitless)") # Table S5 CL BSIZ = 0.09
    e_tumsz_cl_exp_component <- 0.351; label("Power exponent of (TUMSZ/3000) on CL_T (unitless)") # Table S5 CL T,BSIZ = 0.351

    # Categorical covariates: multiplicative factors theta^indicator.
    e_sex_cl_exp_inf <- 1.18; label("Multiplicative factor on CL_inf for men, theta^(1 - SEXF) (unitless)") # Table S5 CL SEX = 1.18
    e_sex_cl_exp_component <- 1.45; label("Multiplicative factor on CL_T for men, theta^(1 - SEXF) (unitless)") # Table S5 CL T,SEX = 1.45
    e_sex_vc <- 1.19; label("Multiplicative factor on V1 for men, theta^(1 - SEXF) (unitless)") # Table S5 V1,SEX = 1.19
    e_tumtp_cll_cl_exp_inf <- 1.47; label("Multiplicative factor on CL_inf for CLL (unitless)") # Table S5 CL CLL = 1.47
    e_tumtp_sll_cl_exp_inf <- 1.38; label("Multiplicative factor on CL_inf for SLL (unitless)") # Table S5 CL SLL = 1.38
    e_tumtp_mcl_cl_exp_inf <- 2.07; label("Multiplicative factor on CL_inf for MCL (unitless)") # Table S5 CL MCL = 2.07
    e_tumtp_cll_cl_exp_component <- 2.25; label("Multiplicative factor on CL_T for CLL (unitless)") # Table S5 CL T,CLL = 2.25
    e_tumtp_mcl_cl_exp_component <- 2.80; label("Multiplicative factor on CL_T for MCL (unitless)") # Table S5 CL T,MCL = 2.80
    e_conmed_chop_cl_exp_kdes <- 0.309; label("Multiplicative factor on k_des with CHOP (unitless)") # Table S5 k des,CHOP = 0.309
    e_conmed_bendamustine_cl_exp_kdes <- 0.591; label("Multiplicative factor on k_des with bendamustine (unitless)") # Table S5 k des,bendamustine = 0.591
    e_conmed_fc_cl_exp_kdes <- 2.64; label("Multiplicative factor on k_des with fludarabine + cyclophosphamide (unitless)") # Table S5 K des,fludarabine/cyclophosphamide = 2.64
    e_tumtp_mzl_cl_exp_kdes <- 0.299; label("Multiplicative factor on k_des for MZL (unitless)") # Table S5 K des,MZL = 0.299

    # Residual error (Table S5). Log-scale ('exponential') residual SD
    # W = sdL - (sdL - sdH) * C / (sd50 + C), equal to (sdL + sdH)/2 at
    # C = sd50 (Table S5 footnote), times exp(etaRUV); sigma^2 fixed to 1.
    sdL <- 1.77; label("Log-scale residual SD at low concentrations (unitless)") # Table S5 sigma_L = 1.77
    sdH <- 0.111; label("Log-scale residual SD at high concentrations (unitless)") # Table S5 sigma_H = 0.111
    sd50 <- 8.26; label("Concentration at which the residual SD is (sdL + sdH)/2 (ug/mL)") # Table S5 sigma_50 = 8.26 ug/mL

    # Inter-individual variability (Table S5; omega^2, log scale; the
    # printed CV is 100 * sqrt(omega^2)). No covariances are reported.
    etalcl_exp_kdes ~ 0.82 # Table S5 omega^2 kdes = 0.82 (CV 90.6%)
    etalcl_exp_component ~ 1.21 # Table S5 omega^2 CLT = 1.21 (CV 110%)
    etalcl_exp_inf ~ 0.199 # Table S5 omega^2 CLinf = 0.199 (CV 44.6%)
    etalvc ~ 0.0341 # Table S5 omega^2 V1 = 0.0341 (CV 18.5%)
    etalvp ~ 0.23 # Table S5 omega^2 V2 = 0.23 (CV 47.9%)
    etalq ~ 0.956 # Table S5 omega^2 Q = 0.956 (CV 97.8%)
    etaRUV ~ 0.224 # Table S5 omega^2 EPS = 0.224 (CV 47.3%), IIV on the residual SD
  })

  model({
    # Covariate ratios and the male indicator (reference subject is female).
    sex_male <- 1 - SEXF
    iwt <- WT / 75
    ialb <- ALB / 40
    iage <- AGE / 65
    isiz <- TUMSZ / 3000

    # Individual parameters
    cl_exp_kdes <- exp(lcl_exp_kdes + etalcl_exp_kdes) *
      e_conmed_chop_cl_exp_kdes^CONMED_CHOP *
      e_conmed_bendamustine_cl_exp_kdes^CONMED_BENDAMUSTINE *
      e_conmed_fc_cl_exp_kdes^CONMED_FC *
      e_tumtp_mzl_cl_exp_kdes^TUMTP_MZL
    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component) *
      e_tumtp_cll_cl_exp_component^TUMTP_CLL *
      e_tumtp_mcl_cl_exp_component^TUMTP_MCL *
      isiz^e_tumsz_cl_exp_component *
      e_sex_cl_exp_component^sex_male
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) *
      iwt^e_wt_cl_exp_inf * ialb^e_alb_cl_exp_inf * iage^e_age_cl_exp_inf *
      isiz^e_tumsz_cl_exp_inf * e_sex_cl_exp_inf^sex_male *
      e_tumtp_cll_cl_exp_inf^TUMTP_CLL * e_tumtp_sll_cl_exp_inf^TUMTP_SLL *
      e_tumtp_mcl_cl_exp_inf^TUMTP_MCL
    vc <- exp(lvc + etalvc) * iwt^e_wt_vc * e_sex_vc^sex_male
    vp <- exp(lvp + etalvp) * iwt^e_wt_vp
    q <- exp(lq + etalq)

    # Total clearance decays from CL_T + CL_inf towards CL_inf (Methods:
    # 'clearance was a sum of time-dependent clearance (CLT) and
    # time-independent clearance (CLinf)'; 'CLT decreased with time with a
    # decay coefficient (Kdes)'). time is the time since the first
    # obinutuzumab dose.
    cl <- cl_exp_component * exp(-cl_exp_kdes * time) + cl_exp_inf

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L = ug/mL.
    Cc <- central / vc

    # Log-scale residual error: Y = log(IPRED) + W * EPS with EPS ~ N(0, 1).
    w <- (sdL - (sdL - sdH) * Cc / (sd50 + Cc)) * exp(etaRUV)
    Cc ~ lnorm(w)
  })
}
