Ahamadi_2017_pembrolizumab <- function() {
  description <- "Two-compartment population PK model for pembrolizumab (humanized anti-PD-1 IgG4 monoclonal antibody) with allometric scaling and covariate effects of sex, albumin, tumor type, ECOG performance status, prior ipilimumab status, eGFR, and baseline tumor burden, in adults with advanced solid tumors (Ahamadi 2017, KEYNOTE-001/-002/-006)"
  reference <- "Ahamadi M, Freshwater T, Prohn M, Li CH, de Alwis DP, de Greef R, Elassaiss-Schaap J, Kondic A, Stone JA. Model-based characterization of the pharmacokinetics of pembrolizumab: a humanized anti-PD-1 monoclonal antibody in advanced solid tumors. CPT Pharmacometrics Syst Pharmacol. 2017;6(1):49-57. doi:10.1002/psp4.12139"
  vignette <- "Ahamadi_2017_pembrolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on shared CL/Q (exponent 0.595) and on shared Vc/Vp (exponent 0.489) with reference 76.8 kg (Ahamadi 2017 Table 3 footnote a/b equations: (WGT/76.8)^alpha). Reference weight is not listed in Table 2 demographics; the value 76.8 kg comes from the denominator in the footnote model equations.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL (exponent -0.907) and on Vc (exponent -0.208) with reference 39.6 g/L (Ahamadi 2017 Table 3 footnote a/b: (ALB/39.6)^theta). The cohort median in Table 2 is 40 g/L.",
      source_name        = "ALB"
    ),
    TUM_SLD = list(
      description        = "Baseline tumor burden (sum of longest diameters of target lesions per RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL (exponent 0.0872) with reference 89.6 mm (Ahamadi 2017 Table 3 footnote a: (BSLD/89.6)^theta). Source paper labels this 'baseline tumor burden (sum of longest dimensions of target lesions)' (BSLD), which maps to the canonical RECIST 1.1 sum-of-longest-diameters metric (TUM_SLD). The cohort median in Table 2 is 86 mm with 11.2% missing data imputed to median.",
      source_name        = "BSLD"
    ),
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL (exponent 0.135) with reference 88.47 mL/min/1.73 m^2 (Ahamadi 2017 Table 3 footnote a: (eGFR/88.47)^theta). Source column name is eGFR; stored under the canonical CRCL. The cohort median in Table 2 is 88.7 mL/min/1.73 m^2.",
      source_name        = "eGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional change on CL ((1 - 0.152) for female) and on Vc ((1 - 0.134) for female) per Ahamadi 2017 Table 3 footnote a/b. The paper's reference category is male, matching the canonical SEXF = 0 (male) convention.",
      source_name        = "SEX"
    ),
    TUMTP_NSCLC = list(
      description        = "Non-small cell lung cancer tumor-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma or other tumor type; melanoma is the implicit reference)",
      notes              = "Proportional change on CL ((1 + 0.145) for NSCLC) per Ahamadi 2017 Table 3 footnote a. Cancer type was tested as a three-level categorical (melanoma 73.7%, NSCLC 25.3%, other 1.01%); only NSCLC vs melanoma was retained in the final model, with the 'other' category pooled into the melanoma reference. Decompose the source TUMTP column into TUMTP_NSCLC = as.integer(TUMTP == 'NSCLC').",
      source_name        = "TUMTP"
    ),
    ECOG_GE1 = list(
      description        = "Eastern Cooperative Oncology Group performance-status indicator (>= 1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG-PS = 0, asymptomatic)",
      notes              = "Proportional change on CL ((1 - 0.0739) for ECOG-PS >= 1) per Ahamadi 2017 Table 3 row (Baseline ECOG-PS on CL = -0.0739). The Discussion confirms this direction: 'Relative to ECOG-PS 1, ECOG-PS 0 was associated with a 7.3% increase in clearance'. Table 3 footnote a prints (1 + 0.0739) for ECOG = 1, which is inconsistent with both the Table 3 sign and the Discussion narrative; treated as a typo in the footnote. ECOG-PS observed values were 0 (57.4%) and 1 (42.4%) with 0.2% missing imputed to mode.",
      source_name        = "ECOG"
    ),
    PRIOR_IPI = list(
      description        = "Prior ipilimumab treatment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ipilimumab-naive; 'missing' is pooled into the reference)",
      notes              = "Proportional change on CL ((1 + 0.140) for IPI-treated) and on Vc ((1 + 0.0736) for IPI-treated) per Ahamadi 2017 Table 3 footnote a/b. Source paper tested IPI status as a three-level categorical (naive 39.1%, treated 34.5%, missing 26.4%) and explicitly kept 'missing' as a separate category in the development dataset; the published Table 3 reports only the naive-vs-treated coefficient. For the packaged model the 'missing' subjects are treated like naive (PRIOR_IPI = 0), so the canonical column carries only the naive/treated contrast that the model coefficient describes.",
      source_name        = "IPI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 2188L,
    n_studies       = 3L,
    age_range       = "15-94 years",
    age_median      = "62 years",
    weight_range    = "not reported in Table 2 demographics; reference 76.8 kg used in the model equation",
    weight_median   = "not reported in Table 2 (only the equation denominator 76.8 kg is published)",
    sex_female_pct  = 40.9,
    race_ethnicity  = "not tabulated by category in the paper (race was tested as a covariate and not retained)",
    disease_state   = "Advanced / metastatic solid tumors (melanoma 73.7%, NSCLC 25.3%, other cancer type 1.01%)",
    dose_range      = "1-10 mg/kg IV infusion Q2W or Q3W (2 mg/kg Q3W is the approved regimen; 10 mg/kg Q3W and 10 mg/kg Q2W were the most populated cohorts)",
    regions         = "Multinational pooled KEYNOTE-001 (NCT01295827), KEYNOTE-002 (NCT01704287), KEYNOTE-006 (NCT01866319)",
    ecog_distribution = "ECOG 0 (asymptomatic) 57.4%, ECOG 1 (symptomatic) 42.4%, missing 0.2%",
    renal_function  = "Baseline eGFR 25.4-403.0 mL/min/1.73 m^2 (median 88.7); 1.2% missing imputed to median",
    notes           = "Baseline demographics per Ahamadi 2017 Table 2 (N = 2,188). KEYNOTE-001 contributed 1,223 patients, KEYNOTE-002 contributed 421 patients, and KEYNOTE-006 contributed 551 patients. Total observations: 12,171 (Table 1 sum). Continuous covariate medians: age 62 y, baseline tumor burden 86 mm, eGFR 88.7 mL/min/1.73 m^2, bilirubin 8.55 umol/L, AST 21 IU/L, albumin 40 g/L. Prior ipilimumab status: IPI-naive 39.1%, IPI-treated 34.5%, missing 26.4% (kept as a separate category during covariate selection per Methods). Coadministered systemic glucocorticoids: yes 14.9%, no 85.1%."
  )

  ini({
    # Structural parameters - typical values at the reference covariate vector
    # (WT 76.8 kg, ALB 39.6 g/L, TUM_SLD 89.6 mm, CRCL 88.47 mL/min/1.73 m^2,
    # male, melanoma, ECOG_GE1 = 0, IPI-naive). Values from Ahamadi 2017 Table 3
    # Estimate column. Q and Vp carry no covariates in the final model (they
    # only inherit body-weight allometry via shared exponents).
    lcl  <- log(0.22);  label("Baseline clearance at reference covariates CL_REF (L/day)")        # Ahamadi 2017 Table 3: CL = 0.22 L/day
    lvc  <- log(3.48);  label("Central volume of distribution Vc_REF (L)")                         # Ahamadi 2017 Table 3: Vc = 3.48 L
    lq   <- log(0.795); label("Intercompartmental clearance Q_REF (L/day)")                        # Ahamadi 2017 Table 3: Q  = 0.795 L/day
    lvp  <- log(4.06);  label("Peripheral volume of distribution Vp_REF (L)")                      # Ahamadi 2017 Table 3: Vp = 4.06 L

    # Allometric exponents (Ahamadi 2017 Table 3 rows 'a for CL and Q' and
    # 'a for Vc or Vp'). The paper uses the same exponent for CL and Q, and
    # the same exponent for Vc and Vp.
    e_wt_cl_q   <- 0.595; label("Power exponent of WT on CL and Q (unitless)")                     # Ahamadi 2017 Table 3: alpha for CL and Q = 0.595
    e_wt_vc_vp  <- 0.489; label("Power exponent of WT on Vc and Vp (unitless)")                    # Ahamadi 2017 Table 3: alpha for Vc and Vp = 0.489

    # Continuous covariate effects on CL - all power-form on (cov / ref)
    # (Ahamadi 2017 Table 3 footnote a equation; THETA values from Table 3
    # Estimate column).
    e_alb_cl    <- -0.907;  label("Power exponent of ALB on CL (unitless)")                        # Ahamadi 2017 Table 3: ALB on CL = -0.907
    e_tsld_cl   <-  0.0872; label("Power exponent of TUM_SLD on CL (unitless)")                    # Ahamadi 2017 Table 3: BSLD on CL = 0.0872
    e_crcl_cl   <-  0.135;  label("Power exponent of CRCL on CL (unitless)")                       # Ahamadi 2017 Table 3: eGFR on CL = 0.135

    # Categorical covariate effects on CL - proportional change
    # (1 + theta * I) per Ahamadi 2017 Methods ('categorical covariates were
    # described as a proportional change') and Table 3 footnote a equation.
    e_sexf_cl       <- -0.152;  label("Proportional change of female sex on CL (unitless)")        # Ahamadi 2017 Table 3: Sex on CL = -0.152
    e_nsclc_cl      <-  0.145;  label("Proportional change of NSCLC tumor type on CL (unitless)")  # Ahamadi 2017 Table 3: Cancer type on CL = 0.145
    e_ecog_ge1_cl   <- -0.0739; label("Proportional change of ECOG_GE1 = 1 on CL (unitless)")      # Ahamadi 2017 Table 3: Baseline ECOG-PS on CL = -0.0739
    e_ipi_cl        <-  0.140;  label("Proportional change of prior ipilimumab on CL (unitless)")  # Ahamadi 2017 Table 3: IPI status on CL = 0.140

    # Continuous covariate effects on Vc - power form (Ahamadi 2017 Table 3
    # footnote b equation).
    e_alb_vc        <- -0.208;  label("Power exponent of ALB on Vc (unitless)")                    # Ahamadi 2017 Table 3: ALB on Vc = -0.208

    # Categorical covariate effects on Vc - proportional change.
    e_sexf_vc       <- -0.134;  label("Proportional change of female sex on Vc (unitless)")        # Ahamadi 2017 Table 3: Sex on Vc = -0.134
    e_ipi_vc        <-  0.0736; label("Proportional change of prior ipilimumab on Vc (unitless)")  # Ahamadi 2017 Table 3: IPI status on Vc = 0.0736

    # Inter-individual variability (Ahamadi 2017 Table 3 Random effect section).
    # The paper describes a SHARED eta on CL and Q (called eta1 in the source)
    # and a SHARED eta on Vc and Vp (called eta2). Table 3 reports
    # omega^2_g1 = 0.134 (38% CV) and omega^2_g2 = 0.0417 (21% CV) but does
    # not tabulate a covariance between eta1 and eta2; that off-diagonal is
    # set to 0 here (see vignette Assumptions and deviations for the
    # development-narrative note about an estimated covariance whose value
    # is not published).
    etalcl + etalvc ~ c(0.134,
                        0.0,    0.0417)                                                            # Ahamadi 2017 Table 3: omega^2_CL/Q = 0.134, cov_CL:Vc set to 0 (not reported), omega^2_Vc/Vp = 0.0417

    # Residual error - the paper used 'additive residual error on the log
    # scale' (Y = LOG(F) + EPS) which is equivalent to a proportional error
    # in linear space with SD = sigma. Table 3 reports the residual-error
    # SD as 0.272 (27% CV).
    propSd <- 0.272; label("Proportional residual error (fraction)")                               # Ahamadi 2017 Table 3: residual error SD = 0.272
  })
  model({
    # Individual PK parameters (Ahamadi 2017 Table 3 footnotes a and b).
    # Reference covariate values: WT 76.8 kg, ALB 39.6 g/L, TUM_SLD 89.6 mm,
    # CRCL 88.47 mL/min/1.73 m^2, male (SEXF = 0), melanoma (TUMTP_NSCLC = 0),
    # ECOG 0 (ECOG_GE1 = 0), ipilimumab-naive (PRIOR_IPI = 0). Continuous
    # covariates are power-form on (cov / ref); categorical covariates are
    # proportional changes (1 + theta * I). The shared etalcl drives both
    # CL and Q; the shared etalvc drives both Vc and Vp.
    cl <- exp(lcl + etalcl) *
      (WT      / 76.8)^e_wt_cl_q *
      (ALB     / 39.6)^e_alb_cl *
      (TUM_SLD / 89.6)^e_tsld_cl *
      (CRCL    / 88.47)^e_crcl_cl *
      (1 + e_sexf_cl     * SEXF) *
      (1 + e_nsclc_cl    * TUMTP_NSCLC) *
      (1 + e_ecog_ge1_cl * ECOG_GE1) *
      (1 + e_ipi_cl      * PRIOR_IPI)

    vc <- exp(lvc + etalvc) *
      (WT  / 76.8)^e_wt_vc_vp *
      (ALB / 39.6)^e_alb_vc *
      (1 + e_sexf_vc * SEXF) *
      (1 + e_ipi_vc  * PRIOR_IPI)

    q  <- exp(lq + etalcl) * (WT / 76.8)^e_wt_cl_q
    vp <- exp(lvp + etalvc) * (WT / 76.8)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and Vc in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
