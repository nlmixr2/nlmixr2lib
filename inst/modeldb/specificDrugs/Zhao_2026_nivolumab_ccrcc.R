Zhao_2026_nivolumab_ccrcc <- function() {
  description <- "Two-compartment population PK model for subcutaneous and intravenous nivolumab (anti-PD-1 IgG4) in previously treated advanced or metastatic clear-cell renal cell carcinoma (CheckMate 67T), with first-order subcutaneous absorption, logit-scale bioavailability and time-varying clearance (sigmoid Emax); the Stage 2 model, estimated with NONMEM $PRIOR on the Stage 1 pooled model, from which the CheckMate 67T co-primary PK endpoints were derived (Zhao 2026)"
  reference <- "Zhao Y, Hu Z, Gibson E, Zhu L, Vezina H, Kondic AG. Model-based clinical pharmacology profiling of subcutaneous nivolumab in renal cell carcinoma: bridging to alternative dosing regimens and indications. Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70481"
  vignette <- "Zhao_2026_nivolumab"
  dosing <- c("depot", "central")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The paper describes a two-compartment model with zero-order intravenous
  # infusion into the sampled serum pool and first-order absorption from the
  # subcutaneous injection site (Discussion; Data S1 'Model development').
  compartmentData <- list(
    depot = list(analyte = "nivolumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL and Vc with reference weight 80 kg (Table S5 Note 1, Figure S1 Note 3). The Data S1 equation legend instead says 'WTB [75 kg]'; the 80 kg in the table notes is used because it is the reference the reported estimates are defined at (see the vignette Errata).",
      source_name = "WTB"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL with reference 90 mL/min/1.73 m^2 (Table S5 Note 1; Data S1 equation legend 'baseline eGFR (90 mL/min/1.73 m2)'). Source column GFRB / EGFR; stored under the canonical CRCL (BSA-normalized renal function). The estimating equation is not named.",
      source_name = "GFRB"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The reference subject is male (Table S5 Note 1). Female sex carries exponential effects on CL and Vc ('e^CL_FEMALE (if SEX=female)') and a multiplicative factor on the typical subcutaneous bioavailability ('F_TV = F_REF x F_FEMALE (if SEX=female)'), so the canonical SEXF orientation (1 = female) applies without sign inversion.",
      source_name = "SEX"
    ),
    ECOG_GE1 = list(
      description = "Baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator (1 if ECOG >= 1, else 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status = 0)",
      notes = "Exponential effect on CL 'e^CL_PS (if PS >= 1)' (Data S1). Table S5 also lists an F1_PS factor held at 1.00 (not estimated), carried here as a fixed() parameter so the effect is inert. Renamed from the source column PS to the canonical ECOG_GE1.",
      source_name = "PS"
    )
  )

  # The Data S1 CL0 equation carries a baseline-albumin power term, but the
  # final-model covariate list in Data S1 omits albumin and Table S5 reports
  # no albumin exponent, so the effect is not part of the reported model.
  covariatesDataExcluded <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Appears in the generic CL0 equation of Data S1 as (ALBB/ALBB_REF)^CL_ALBB with reference 4.0 g/dL, but is absent from the Data S1 list of final-model covariate effects and from Table S5, so no point estimate is available. Summarized in Table 1 (overall median 4 g/dL, range 2.1-5.1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 487L,
    n_studies = 1L,
    age_range = "median 65 years (range 20 - 93); mean 63.9 (SD 10.3)",
    age_median = "65 years",
    weight_range = "median 76.5 kg (range 35 - 157); mean 77.7 (SD 17.3)",
    weight_median = "76.5 kg",
    sex_female_pct = 32.4,
    race_ethnicity = c(
      White = 85.0,
      Other = 12.5,
      `American Indian/Alaska Native` = 1.0,
      Asian = 0.8,
      `Black/African American` = 0.4,
      Missing = 0.2
    ),
    disease_state = "Advanced or metastatic clear-cell renal cell carcinoma after prior systemic therapy (second-line)",
    dose_range = "Subcutaneous nivolumab 1200 mg co-formulated with rHuPH20 20,000 units Q4W (N = 242) or intravenous nivolumab 3 mg/kg Q2W as a 30-minute infusion (N = 245)",
    regions = "Global (multicenter phase III)",
    ecog_distribution = "ECOG 0 39.2%, ECOG 1 54.0%, ECOG 2 6.8%",
    renal_function = "Baseline eGFR median 61.2 (range 20.7 - 124) mL/min/1.73 m^2; mean 63.7 (SD 19.9)",
    albumin = "Baseline albumin median 4 (range 2.1 - 5.1) g/dL; 1.44% missing",
    n_subjects_sc = 242L,
    n_subjects_iv = 245L,
    notes = "CheckMate 67T (CA20967T, NCT04810078) Stage 2 analysis population, Table 1 (N = 487; 3734 PK samples per Table S2). Only CheckMate 67T data were fitted; information from the Stage 1 pooled model (17 historical intravenous studies plus CheckMate 8KX) entered through NONMEM $PRIOR on CL, VC, VP, Q, EMAX, HILL, T50, KA and F and on the OMEGAs of CL, VC, VP and EMAX (Data S1). Nivolumab ADA: 72.9% negative, 12.5% any positive, 14.6% missing."
  )

  ini({
    # Structural parameters, Table S5 ('Parameter estimates of the final
    # population pharmacokinetic model from stage 2'). Reference subject:
    # 80 kg, eGFR 90 mL/min/1.73 m^2, male, PS 0 (Table S5 Note 1). The
    # source time unit is hours; rates are converted to days (x 24).
    lcl <- log(0.0108 * 24); label("Baseline clearance CL0_REF at the reference covariates (L/day)")     # Table S5: CL0_REF = 0.0108 L/hr; x 24 = 0.2592 L/day
    lvc <- log(4.25);        label("Central volume of distribution VC_REF (L)")                          # Table S5: VC_REF = 4.25 L
    lq  <- log(0.0312 * 24); label("Intercompartmental clearance Q_REF (L/day)")                         # Table S5: Q_REF = 0.0312 L/hr; x 24 = 0.7488 L/day
    lvp <- log(2.67);        label("Peripheral volume of distribution VP_REF (L)")                       # Table S5: VP_REF = 2.67 L

    # Subcutaneous absorption.
    lka <- log(0.0123 * 24);                  label("First-order subcutaneous absorption rate constant KA_REF (1/day)")     # Table S5: KA = 0.0123 /hr; x 24 = 0.2952 /day (Results: '0.0123 hour-1 (or 0.295 day-1)')
    logitfdepot <- log(0.788 / (1 - 0.788));  label("Logit of typical subcutaneous bioavailability F_REF (fraction)")       # Table S5: F1 = 0.788; logit(0.788) = 1.3129

    # Covariate effects on CL0 (Data S1 CL0_TV equation): power on WT and
    # eGFR; exponential on female sex and PS >= 1.
    e_wt_cl       <-  0.644; label("Power exponent of WT on CL (unitless)")                                  # Table S5: CL_WTB = 0.644
    e_crcl_cl     <-  0.306; label("Power exponent of CRCL (eGFR) on CL (unitless)")                         # Table S5: CL_GFR = 0.306
    e_sex_cl      <- -0.149; label("Exponential coefficient of female sex on CL (unitless)")                 # Table S5: CL_SEX = -0.149
    e_ecog_ge1_cl <-  0.122; label("Exponential coefficient of ECOG_GE1 on CL (unitless)")                   # Table S5: CL_PS = 0.122

    # Covariate effects on Vc (Data S1 VC equation): power on WT; exponential
    # on female sex.
    e_wt_vc  <-  0.525; label("Power exponent of WT on Vc (unitless)")                                       # Table S5: VC_WTB = 0.525
    e_sex_vc <- -0.228; label("Exponential coefficient of female sex on Vc (unitless)")                      # Table S5: VC_SEX = -0.228

    # Covariate effects on the typical subcutaneous bioavailability, applied
    # multiplicatively on the fraction scale before the logit transform
    # (Data S1 'F_TV,i = F_REF x F_FEMALE (if SEX=female)').
    e_sex_fdepot      <- 0.942;    label("Multiplicative factor of female sex on subcutaneous bioavailability (unitless)")  # Table S5: F1_SEX = 0.942
    e_ecog_ge1_fdepot <- fixed(1); label("Multiplicative factor of ECOG_GE1 on subcutaneous bioavailability (unitless)")    # Table S5: F1_PS = 1.00, held constant by the authors (superscript f)

    # Time-varying clearance, sigmoid Emax in time since first dose (Data S1
    # CL_i(t) equation). EMAX is the maximal change in log CL: CL approaches
    # CL0 x exp(-0.260) = 0.771 x baseline. T50 is converted from hours to
    # days; the Hill ratio is unit-invariant when t and T50 share a unit.
    cl_time_max   <- -0.260;       label("Maximal change in log CL, EMAX_REF (unitless)")                    # Table S5: CL_EMAX = -0.260
    lcl_t50       <- log(1.41e3 / 24); label("Log of the time to 50%% of the maximal change in CL, T50 (log days)")         # Table S5: CL_T50 = 1.41E+03 hr; / 24 = 58.75 days
    lcl_time_hill <- log(3.39);       label("Log of the Hill / sigmoidicity exponent of time on CL (log unitless)")  # Table S5: CL_HILL = 3.39

    # Inter-individual variability (Table S5 'Random effects'), printed as
    # variance (SD) on the diagonal and covariance (correlation) off it; each
    # parenthetical reproduces (sqrt(0.135) = 0.367; 0.0395 / sqrt(0.135 x
    # 0.0827) = 0.374; 0.0744 / sqrt(0.211 x 0.615) = 0.207). CL, Vc, Vp and
    # Ka are log-normal; EMAX carries an additive eta; the bioavailability
    # eta acts on the logit scale (Data S1 F_i equation).
    etalcl + etalvc ~ c(0.135,
                        0.0395, 0.0827)                                                                      # Table S5: ZCL 0.135 (0.367), ZCL:ZVC 0.0395 (0.374), ZVC 0.0827 (0.288)
    etalvp ~ 0.211                                                                                           # Table S5: ZVP 0.211 (0.459)
    etacl_time_max ~ 0.0572                                                                                  # Table S5: ZEMAX 0.0572 (0.239); additive on EMAX
    etalka + etalogitfdepot ~ c(0.211,
                                0.0744, 0.615)                                                               # Table S5: ZKA 0.211 (0.459), ZKA:F1 0.0744 (0.207), F1 0.615 (0.784)

    # Residual error: proportional. Reported as a THETA (Table S5 symbol
    # theta10), i.e. the proportional standard deviation, as in the
    # sibling pre-specified model's control streams (Zhao 2025, where the
    # same THETA multiplies IPRED with $SIGMA 1 FIX).
    propSd <- 0.157; label("Proportional residual error (fraction)")                                         # Table S5: Proportional = 0.157 (theta10)
  })
  model({
    # Individual baseline CL and Vc (Data S1). Every covariate term equals 1
    # for the reference subject (80 kg, eGFR 90, male, PS 0).
    cl_base <- exp(lcl + etalcl) *
      (WT / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_sex_cl * SEXF) *
      exp(e_ecog_ge1_cl * ECOG_GE1)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sex_vc * SEXF)

    # Time-varying CL with additive IIV on EMAX; t is time since first dose.
    cl_t50 <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    ka <- exp(lka + etalka)

    # Subcutaneous bioavailability: covariate factors on the fraction scale,
    # IIV on the logit scale (Data S1 F_i and F_TV,i equations).
    fdepot_tv <- expit(logitfdepot) *
      e_sex_fdepot^SEXF *
      e_ecog_ge1_fdepot^ECOG_GE1
    fdepot <- expit(logit(fdepot_tv) + etalogitfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Subcutaneous doses go to depot (absorbed with ka, fraction fdepot);
    # intravenous doses are infused into central with F = 1.
    f(depot) <- fdepot

    # Dose in mg and volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
