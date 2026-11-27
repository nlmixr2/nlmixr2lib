Jamois_2019_obinutuzumab <- function() {
  description <- "Two-compartment population PK model of intravenous obinutuzumab (glycoengineered type II anti-CD20 mAb) in adults with CD20+ B-cell malignancies, updated with the GALLIUM first-line follicular lymphoma data; total clearance is a time-independent component cl_exp_inf plus a time-dependent component cl_exp_component * exp(-cl_exp_kdes * time) whose decay rate depends on the chemotherapy backbone (CHOP/CVP, bendamustine, FC). Covariates are body weight, sex, baseline tumor size, serum albumin, age, histology (CLL, SLL, MCL, MZL vs FL) and anti-drug antibodies; log-scale residual error whose SD falls with concentration, with IIV on the residual SD (Jamois 2019)."
  reference <- "Jamois C, Gibiansky E, Gibiansky L, Buchheit V, Sahin D, Cartron G, Marcus R, Hiddemann W, Seymour JF, Strefford JC, Hargreaves CE, Meneses-Lorente G, Frey N, Fingerle-Rowson G. Role of obinutuzumab exposure on clinical outcome of follicular lymphoma treated with first-line immunochemotherapy. Br J Clin Pharmacol. 2019;85(7):1495-1506. doi:10.1111/bcp.13920"
  vignette <- "Jamois_2019_obinutuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The IIV on the residual-error SD (Table S1 Omega(7,7), 'EPS') has no
  # structural typical-value partner; the residual SD is the
  # concentration-dependent expression built from sdL / sdH / sd50 in model().
  paper_specific_etas <- c("etaRUV")

  compartmentData <- list(
    central = list(analyte = "obinutuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "obinutuzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 75 kg (Table S2). Power exponents 0.721 on CL_inf, 0.437 on V1 and 0.776 on V2 (Table S1 theta7-theta9). Table S1 lists no weight effect on CL_T or Q, so none is applied.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female; Table S2 reference 'Female')",
      notes = "The paper reports male-vs-female ratios (Table S1 theta17 on CL_T = 1.48, theta22 on CL_inf = 1.18, theta25 on V1 = 1.19; Table S2 reference 'Female'), so the typical values are female values and the effects are applied through (1 - SEXF).",
      source_name = "SEX"
    ),
    TUMSZ = list(
      description = "Baseline tumor size, sum of the products of the perpendicular diameters of the target lesions (SPD)",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value (BSIZ), assessed by CT (Methods 'Tumour burden assessment'). Reference 3000 mm^2 (Table S2). Power exponents 0.345 on CL_T (theta16) and 0.0596 on CL_inf (theta21). Median 5109 mm^2 (range 235-80,720) in the GALLIUM PK/PFS population (Table S3). The paper does not say how missing values were handled.",
      source_name = "BSIZ"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 40 g/L (Table S2). Power exponent -0.493 on CL_inf (Table S1 theta24). The 2.5th-97.5th percentiles of the pooled dataset were 29.0-48.6 g/L (Table S2 footnote).",
      source_name = "ALB"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 65 years (Table S2). Power exponent -0.211 on CL_inf (Table S1 theta23). The 2.5th-97.5th percentiles of the pooled dataset were 38-82 years (Table S2).",
      source_name = "AGE"
    ),
    TUMTP_CLL = list(
      description = "Chronic lymphocytic leukemia indicator (1 = CLL, 0 = other histology)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (follicular lymphoma, also DLBCL; Table S2)",
      notes = "Multiplies CL_T by 2.28 (theta14) and CL_inf by 1.45 (theta19). Table S2 assigns the 'CL_CLL' effect to CL_inf only; its CL_T entry (+128.3%) is theta14 alone.",
      source_name = "disease type"
    ),
    TUMTP_SLL = list(
      description = "Small lymphocytic lymphoma indicator (1 = SLL, 0 = other histology)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (follicular lymphoma, MZL or DLBCL; Table S2)",
      notes = "Multiplies CL_inf by 1.38 (theta18). No SLL effect on CL_T in Table S1/S2. Mutually exclusive with TUMTP_CLL.",
      source_name = "disease type"
    ),
    TUMTP_MCL = list(
      description = "Mantle cell lymphoma indicator (1 = MCL, 0 = other histology)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (follicular lymphoma, also DLBCL; Table S2)",
      notes = "Multiplies CL_T by 2.98 (theta15) and CL_inf by 2.04 (theta20).",
      source_name = "disease type"
    ),
    TUMTP_MZL = list(
      description = "Marginal zone lymphoma indicator (1 = MZL, 0 = other histology)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (follicular lymphoma or DLBCL; Table S2)",
      notes = "Multiplies CL_T by 1.40 (theta13). MZL shares the CL_inf reference with FL and DLBCL (Table S2).",
      source_name = "disease type"
    ),
    CONMED_CHOP = list(
      description = "Concomitant CHOP or CVP chemotherapy backbone (1 = obinutuzumab given with CHOP or CVP, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chemotherapy: obinutuzumab monotherapy)",
      notes = "One shared effect for CHOP and CVP (Table S1 'k_des,CHOP-CVP'; Table S2 'CHOP or CVP'): multiplies k_des by 0.285 (theta10), giving a 21.6-day half-life of the decay of the time-dependent clearance in GALLIUM (Results). CVP shares three of the four CHOP drugs, which is the register's 'CHOP-like' definition. Mutually exclusive with CONMED_BENDAMUSTINE and CONMED_FC.",
      source_name = "CHOP-CVP"
    ),
    CONMED_BENDAMUSTINE = list(
      description = "Concomitant bendamustine (1 = obinutuzumab given with bendamustine, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chemotherapy: obinutuzumab monotherapy)",
      notes = "Multiplies k_des by 0.466 (Table S1 theta11), giving a 13.2-day half-life of the decay of the time-dependent clearance (Results).",
      source_name = "bendamustine"
    ),
    CONMED_FC = list(
      description = "Concomitant fludarabine plus cyclophosphamide (FC) (1 = obinutuzumab given with FC, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chemotherapy: obinutuzumab monotherapy)",
      notes = "Multiplies k_des by 2.88 (Table S1 theta12).",
      source_name = "FC"
    ),
    ADA_POS = list(
      description = "Anti-drug antibodies detected (1 = detected, 0 = not detected)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not detected)",
      notes = "Multiplies CL_T by 3.09 (Table S1 theta29) and k_des by 1.74 (theta30). The paper does not say whether the indicator was time-varying or ever-positive.",
      source_name = "ADA"
    )
  )

  population <- list(
    species = "human",
    disease_state = "CD20+ B-cell malignancies: the obinutuzumab population PK dataset of Gibiansky 2014 (678 patients with CLL, follicular and other indolent NHL, DLBCL or MCL from GAUGUIN, GAUDI, GAUSS and CLL11) plus 408 patients with previously untreated advanced follicular lymphoma from the phase III GALLIUM trial (Results 'PK analysis'). The pooled total is not stated in the paper.",
    age_range = "2.5th-97.5th percentile 38-82 years in the pooled dataset (Table S2); GALLIUM PK/PFS population median 59 (27-88) years (Table S3)",
    weight_range = "2.5th-97.5th percentile 50-115 kg in the pooled dataset (Table S2); GALLIUM PK/PFS population median 73.7 (35.3-155) kg (Table S3)",
    sex_female_pct = 52.6,
    tumor_size_summary = "GALLIUM PK/PFS population baseline tumour burden (SPD) median 5109 mm^2 (range 235.2-80,720); pooled-dataset 2.5th-97.5th percentiles 344-27,500 mm^2 (Tables S2 and S3)",
    albumin_summary = "Pooled-dataset 2.5th-97.5th percentiles 29.0-48.6 g/L (Table S2)",
    dose_range = "GALLIUM: intravenous obinutuzumab 1000 mg on days 1, 8 and 15 of cycle 1 and day 1 of later cycles, with bendamustine (6 x 28-day cycles, 8000 mg in total) or CHOP/CVP (8 x 21-day cycles, 10,000 mg in total), then 1000 mg every 2 months maintenance for up to 2 years in responders. Earlier studies: 200-2000 mg per infusion.",
    regions = "International multicentre trials",
    studies = "GAUGUIN, GAUDI, GAUSS and CLL11 (Gibiansky 2014) plus GALLIUM (NCT01332968).",
    notes = "sex_female_pct is from the 401-patient GALLIUM PK/PFS population (190 male, Table S3), not the pooled PK dataset. Serum obinutuzumab LLOQ 4.05 ng/mL. NONMEM 7.3.0."
  )

  ini({
    # Structural PK parameters (Table S1, final model). Reference subject:
    # female, 75 kg, 65 years, albumin 40 g/L, baseline tumour size 3000 mm^2,
    # follicular lymphoma, no chemotherapy, no anti-drug antibodies (Table S2).
    lcl_exp_kdes <- log(0.112); label("Decay rate constant of the time-dependent clearance k_des (1/day)") # Table S1 exp(theta1) k_des = 0.112 1/day
    lcl_exp_component <- log(0.145); label("Initial time-dependent clearance CL_T (L/day)") # Table S1 exp(theta2) CL_T = 0.145 L/day
    lcl_exp_inf <- log(0.0753); label("Time-independent clearance CL_inf (L/day)") # Table S1 exp(theta3) CL_inf = 0.0753 L/day
    lvc <- log(2.59); label("Central volume of distribution V1 (L)") # Table S1 exp(theta4) V1 = 2.59 L
    lvp <- log(1.47); label("Peripheral volume of distribution V2 (L)") # Table S1 exp(theta5) V2 = 1.47 L
    lq <- log(1.32); label("Intercompartmental clearance Q (L/day)") # Table S1 exp(theta6) Q = 1.32 L/day

    # Body weight (power on WT/75; Table S2 reference 75 kg)
    e_wt_cl_exp_inf <- 0.721; label("Power exponent of (WT/75) on CL_inf (unitless)") # Table S1 theta7 CL_inf,WT = 0.721
    e_wt_vc <- 0.437; label("Power exponent of (WT/75) on V1 (unitless)") # Table S1 theta8 V1,WT = 0.437
    e_wt_vp <- 0.776; label("Power exponent of (WT/75) on V2 (unitless)") # Table S1 theta9 V2,WT = 0.776

    # Chemotherapy backbone on k_des (log of the reported ratio exp(theta))
    e_conmed_chop_cl_exp_kdes <- log(0.285); label("Log ratio of k_des with CHOP or CVP vs no chemotherapy (unitless)") # Table S1 exp(theta10) k_des,CHOP-CVP = 0.285
    e_conmed_bendamustine_cl_exp_kdes <- log(0.466); label("Log ratio of k_des with bendamustine vs no chemotherapy (unitless)") # Table S1 exp(theta11) k_des,bendamustine = 0.466
    e_conmed_fc_cl_exp_kdes <- log(2.88); label("Log ratio of k_des with fludarabine + cyclophosphamide vs no chemotherapy (unitless)") # Table S1 exp(theta12) k_des,FC = 2.88

    # Effects on the time-dependent clearance CL_T
    e_mzl_cl_exp_component <- log(1.4); label("Log ratio of CL_T for MZL vs FL (unitless)") # Table S1 exp(theta13) CL_T,MZL = 1.4
    e_cll_cl_exp_component <- log(2.28); label("Log ratio of CL_T for CLL vs FL (unitless)") # Table S1 exp(theta14) CL_T,CLL = 2.28
    e_mcl_cl_exp_component <- log(2.98); label("Log ratio of CL_T for MCL vs FL (unitless)") # Table S1 exp(theta15) CL_T,MCL = 2.98
    e_tumsz_cl_exp_component <- 0.345; label("Power exponent of (TUMSZ/3000) on CL_T (unitless)") # Table S1 theta16 CL_T,BSIZ = 0.345
    e_sex_cl_exp_component <- log(1.48); label("Log ratio of CL_T for males vs females, applied via 1 - SEXF (unitless)") # Table S1 exp(theta17) CL_T,SEX = 1.48

    # Effects on the time-independent clearance CL_inf (Table S2 groups
    # theta18-theta24 under CL_inf)
    e_sll_cl_exp_inf <- log(1.38); label("Log ratio of CL_inf for SLL vs FL (unitless)") # Table S1 exp(theta18) CL_SLL = 1.38
    e_cll_cl_exp_inf <- log(1.45); label("Log ratio of CL_inf for CLL vs FL (unitless)") # Table S1 exp(theta19) CL_CLL = 1.45
    e_mcl_cl_exp_inf <- log(2.04); label("Log ratio of CL_inf for MCL vs FL (unitless)") # Table S1 exp(theta20) CL_MCL = 2.04
    e_tumsz_cl_exp_inf <- 0.0596; label("Power exponent of (TUMSZ/3000) on CL_inf (unitless)") # Table S1 theta21 CL_BSIZ = 0.0596
    e_sex_cl_exp_inf <- log(1.18); label("Log ratio of CL_inf for males vs females, applied via 1 - SEXF (unitless)") # Table S1 exp(theta22) CL_SEX = 1.18
    e_age_cl_exp_inf <- -0.211; label("Power exponent of (AGE/65) on CL_inf (unitless)") # Table S1 theta23 CL_AGE = -0.211
    e_alb_cl_exp_inf <- -0.493; label("Power exponent of (ALB/40) on CL_inf (unitless)") # Table S1 theta24 CL_ALB = -0.493

    # Sex on V1
    e_sex_vc <- log(1.19); label("Log ratio of V1 for males vs females, applied via 1 - SEXF (unitless)") # Table S1 exp(theta25) V1,SEX = 1.19

    # Anti-drug antibodies
    e_ada_cl_exp_component <- log(3.09); label("Log ratio of CL_T with anti-drug antibodies detected (unitless)") # Table S1 exp(theta29) CL_T,ADA = 3.09
    e_ada_cl_exp_kdes <- log(1.74); label("Log ratio of k_des with anti-drug antibodies detected (unitless)") # Table S1 exp(theta30) k_des,ADA = 1.74

    # Residual error: exponential (log-scale) error whose SD follows an Emax
    # function of the concentration, falling from sdL at low to sdH at high
    # concentrations, equal to (sdL + sdH)/2 at sd50; sigma^2 is fixed to 1,
    # so the thetas are SDs.
    sdL <- 1.59; label("Log-scale residual SD at low concentrations (unitless)") # Table S1 theta26 sigma_L = 1.59
    sdH <- 0.121; label("Log-scale residual SD at high concentrations (unitless)") # Table S1 theta27 sigma_H = 0.121
    sd50 <- 7.55; label("Concentration at which the residual SD is (sdL + sdH)/2 (ug/mL)") # Table S1 theta28 sigma_50 = 7.55 ug/mL

    # Inter-individual variability (Table S1; omega^2 on the log scale, diagonal)
    etalcl_exp_kdes ~ 0.821 # Table S1 Omega(1,1) k_des = 0.821 (CV 90.6%)
    etalcl_exp_component ~ 1.09 # Table S1 Omega(2,2) CL_T = 1.09 (CV 104%)
    etalcl_exp_inf ~ 0.163 # Table S1 Omega(3,3) CL_inf = 0.163 (CV 40.4%)
    etalvc ~ 0.0355 # Table S1 Omega(4,4) V1 = 0.0355 (CV 18.9%)
    etalvp ~ 0.143 # Table S1 Omega(5,5) V2 = 0.143 (CV 37.8%)
    etalq ~ 0.997 # Table S1 Omega(6,6) Q = 0.997 (CV 99.9%)
    etaRUV ~ 0.191 # Table S1 Omega(7,7) EPS = 0.191 (CV 43.6%), IIV on the residual SD
  })

  model({
    # Male indicator: the reported sex effects are male-vs-female ratios.
    sex_male <- 1 - SEXF

    # Individual parameters (power models on continuous covariates normalised
    # to the Table S2 reference values; exp(theta)^indicator on binary ones)
    cl_exp_kdes <- exp(
      lcl_exp_kdes + etalcl_exp_kdes +
        e_conmed_chop_cl_exp_kdes * CONMED_CHOP +
        e_conmed_bendamustine_cl_exp_kdes * CONMED_BENDAMUSTINE +
        e_conmed_fc_cl_exp_kdes * CONMED_FC +
        e_ada_cl_exp_kdes * ADA_POS
    )
    cl_exp_component <- exp(
      lcl_exp_component + etalcl_exp_component +
        e_mzl_cl_exp_component * TUMTP_MZL +
        e_cll_cl_exp_component * TUMTP_CLL +
        e_mcl_cl_exp_component * TUMTP_MCL +
        e_sex_cl_exp_component * sex_male +
        e_ada_cl_exp_component * ADA_POS
    ) * (TUMSZ / 3000)^e_tumsz_cl_exp_component
    cl_exp_inf <- exp(
      lcl_exp_inf + etalcl_exp_inf +
        e_sll_cl_exp_inf * TUMTP_SLL +
        e_cll_cl_exp_inf * TUMTP_CLL +
        e_mcl_cl_exp_inf * TUMTP_MCL +
        e_sex_cl_exp_inf * sex_male
    ) * (WT / 75)^e_wt_cl_exp_inf * (TUMSZ / 3000)^e_tumsz_cl_exp_inf *
      (AGE / 65)^e_age_cl_exp_inf * (ALB / 40)^e_alb_cl_exp_inf
    vc <- exp(lvc + etalvc + e_sex_vc * sex_male) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    q <- exp(lq + etalq)

    # Total clearance: CL(t) = CL_inf + CL_T * exp(-k_des * t), with t the
    # time since the first obinutuzumab dose (Gibiansky 2014 structure, which
    # this model updates; Discussion).
    cl <- cl_exp_component * exp(-cl_exp_kdes * time) + cl_exp_inf

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give ug/mL.
    Cc <- central / vc

    # Log-scale residual error: log(Y) = log(IPRED) + W * EPS, EPS ~ N(0, 1).
    w <- (sdL - (sdL - sdH) * Cc / (sd50 + Cc)) * exp(etaRUV)
    Cc ~ lnorm(w)
  })
}
