Wang_2024_sugemalimab <- function() {
  description <- "Two-compartment population PK model with sigmoidal-emax time-varying clearance for intravenous sugemalimab (anti-PD-L1 IgG4) in adults with advanced solid tumours or lymphomas across nine Phase I-III trials (Wang 2024)"
  reference <- "Wang K, Pan C, Xu F, Tse AN, Sheng Y. Comprehensive population pharmacokinetic modelling of sugemalimab, an anti-programmed death-ligand 1 (PD-L1) human monoclonal antibody, in patients with solid tumours or lymphomas across multiple Phase I-III studies. Br J Clin Pharmacol. 2025;91(3):748-760. doi:10.1111/bcp.16276"
  vignette <- "Wang_2024_sugemalimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL and on Vc. Reference 61 kg from Wang 2024 Table 3 footnote (typical lung cancer male patient).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL and on Vc. Reference 41.5 g/L from Wang 2024 Table 3 footnote (typical lung cancer male patient). Source paper reports albumin in g/L (SI convention).",
      source_name        = "ALB"
    ),
    TUMSZ = list(
      description        = "Baseline tumour burden (sum of longest target-lesion diameters by RECIST, investigator-assessed)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL only (no Vc effect retained). Reference 47 mm from Wang 2024 Table 3 footnote (typical lung cancer male patient). Source column TUMORB; Wang 2024 reports tumour burden in mm directly. 14.7% of subjects had missing tumour burden in the source dataset (Wang 2024 Table 2 footnote b); imputation handling not stated in the paper.",
      source_name        = "TUMORB"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on baseline CL and on Vc; female reference value is lower for both parameters. Source paper Table 3 footnote: 'SEX: 0 for male, 1 for female', which matches canonical SEXF directly.",
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Anti-drug-antibody positivity status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Exponential effect on baseline CL only (no Vc effect retained). Source column ADA; Wang 2024 Table 3 footnote: 'ADA: 0 for negative, 1 for positive'. ADA-positive subjects had ~10.5% higher CL. ADA positivity rate in the source dataset was 8.8% (143 / 1628; Wang 2024 Table 2). Wang 2024 treats ADA as a time-fixed indicator (not time-varying) per the paper's Discussion: 'the ADA positivity rate of sugemalimab is low (8.7%), and therefore the time-dependent effect of ADA on PK was not considered'.",
      source_name        = "ADA"
    ),
    TUMTP_LYMPH = list(
      description        = "Lymphoma tumour-type indicator (heterogeneous lymphoma pool)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-lymphoma tumour type — when paired with TUMTP_OTH = TUMTP_GC = TUMTP_ESCC = 0, the reference is lung cancer / NSCLC)",
      notes              = "Exponential effect on baseline CL and on Vc. Wang 2024 pools two lymphoma histologies (extranodal NK/T-cell lymphoma from CS1001-201 / NCT03595657 and classical / relapsed-refractory Hodgkin lymphoma from CS1001-202 / NCT03505996) into a single indicator (n = 164 in the pooled dataset). Source column TTYPE level 1.",
      source_name        = "TTYPE1"
    ),
    TUMTP_OTH = list(
      description        = "'Other' tumour-type residual indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (one of the named groups: lung cancer reference, lymphoma, GCGEJ, ESCC)",
      notes              = "Exponential effect on baseline CL and on Vc. Heterogeneous solid-tumour residual group spanning miscellaneous histologies not captured by the named groups (n = 174 in Wang 2024). Source column TTYPE level 3.",
      source_name        = "TTYPE3"
    ),
    TUMTP_GC = list(
      description        = "Gastric / gastroesophageal-junction adenocarcinoma indicator (GCGEJ)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-GCGEJ tumour type)",
      notes              = "Exponential effect on baseline CL and on Vc. Wang 2024 GCGEJ group pools gastric adenocarcinoma (GC) and adenocarcinoma of the gastroesophageal junction (GEJ) into a single indicator (n = 275 in the pooled dataset; primarily from CS1001-303 / NCT03802591). Source column TTYPE level 4.",
      source_name        = "TTYPE4"
    ),
    TUMTP_ESCC = list(
      description        = "Oesophageal squamous cell carcinoma (ESCC) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-ESCC tumour type)",
      notes              = "Exponential effect on baseline CL and on Vc. n = 401 in the pooled dataset (primarily from CS1001-304 / NCT04187352). Source column TTYPE level 5.",
      source_name        = "TTYPE5"
    )
  )

  population <- list(
    n_subjects     = 1628L,
    n_studies      = 9L,
    n_observations = 11040L,
    age_range      = "18-78 years",
    age_median     = "60 years",
    weight_range   = "36.0-141.0 kg",
    weight_median  = "60 kg",
    sex_female_pct = 21.6,
    race_ethnicity = c(White = 2.15, Asian = 97.6, Other = 0.246),
    disease_state  = "Adults with advanced solid tumours or lymphomas. Tumor-type mix in the pooled dataset (n = 1628): lung cancer 614 (37.7%; primarily NSCLC), ESCC 401 (24.6%), GCGEJ 275 (16.9%), 'Other' 174 (10.7%), lymphoma 164 (10.1%; mix of extranodal NK/T-cell lymphoma and classical Hodgkin lymphoma).",
    dose_range     = "Sugemalimab 3-40 mg/kg or 1200 mg IV once every 3 weeks (Q3W); single Phase Ib/II combination cohort received 1800 mg Q4W. The vast majority of subjects received 1200 mg Q3W, reflecting the approved adult dose.",
    regions        = "Pooled global trials, predominantly Asia (Asian race = 97.6%); 35 White and 4 'Other' subjects.",
    ada_status     = "ADA-negative 1485 (91.2%); ADA-positive 143 (8.8%).",
    ecog_status    = "ECOG performance status 0 in 414 (25.4%); ECOG >= 1 in 1213 (74.5%); missing in 1.",
    tumour_burden  = "47 mm median (10-343 mm range); 239 (14.7%) subjects had missing tumour burden.",
    notes          = "Baseline demographics per Wang 2024 Table 2. Pooled dataset spans nine Phase I-III sugemalimab trials (NCT03312842 phase Ia/Ib, NCT03744403 phase I, NCT03595657 phase II ENKTL, NCT03505996 phase II r/r-cHL, NCT04200404 phase Ib/II + regorafenib, NCT03728556 phase III stage-III NSCLC, NCT03789604 phase III stage-IV NSCLC, NCT03802591 phase III GCGEJ, NCT04187352 phase III ESCC). The model was fit using NONMEM 7.5 with FOCE-I; final-model parameter estimates and the simulation-target reference-patient definition are in Wang 2024 Table 3."
  )

  ini({
    # Structural parameters - reference patient: typical lung cancer male
    # patient, WT 61 kg, ALB 41.5 g/L, TUMORB 47 mm, SEX = 0 (male),
    # ADA = 0, all TUMTP_* indicators = 0 (NSCLC reference), per the
    # Wang 2024 Table 3 footnote ("For a typical lung cancer (LC) male
    # patient with ALB of 41.5 g/L, body weight of 61 kg, tumour burden of
    # 47 mm and negative ADA, the estimated CL was 0.259 L/day ...").
    lcl     <- log(0.259); label("Baseline clearance CL0 at reference covariates (L/day)")     # Wang 2024 Table 3: exp(theta1) = 0.259 L/day
    lvc     <- log(3.42);  label("Central volume Vc at reference covariates (L)")              # Wang 2024 Table 3: exp(theta2) = 3.42 L
    lq      <- log(0.489); label("Intercompartmental clearance Q (L/day)")                     # Wang 2024 Table 3: exp(theta3) = 0.489 L/day
    lvp     <- log(1.57);  label("Peripheral volume Vp (L)")                                   # Wang 2024 Table 3: exp(theta4) = 1.57 L

    # Time-varying CL Hill-emax function (Wang 2024 Methods equation):
    #   CL(t) = exp(CL0 + emax * t^lambda / (T50^lambda + t^lambda))
    # so CL(0) = exp(CL0) and CL(t -> inf) = exp(CL0 + emax). With
    # emax = -0.528, the asymptotic CL is exp(-0.528) = 0.590 of baseline,
    # i.e., a 41.0% reduction at full saturation.
    emax     <-      -0.528;     label("Reference emax of time-varying CL (theta5; unitless, negative = CL decreases at full saturation)") # Wang 2024 Table 3: theta5 = -0.528
    lt50     <- log(53.6);       label("log T50 - time at which half of emax is reached (log days)") # Wang 2024 Table 3: theta6 = 53.6 day
    llambda  <- log(2.60);       label("log lambda - sigmoidicity (Hill coefficient) of the time-on-CL function (log unitless)") # Wang 2024 Table 3: theta7 = 2.60

    # Continuous covariate effects on CL and Vc - power form on log-transformed
    # covariates, equivalent to (cov / ref)^exponent in linear space (Wang 2024
    # Table 3 footnote final-model equations).
    e_wt_cl     <-  0.585;  label("Power exponent of body weight on baseline CL (unitless)")        # Wang 2024 Table 3: theta8 = 0.585
    e_wt_vc     <-  0.471;  label("Power exponent of body weight on Vc (unitless)")                 # Wang 2024 Table 3: theta9 = 0.471
    e_alb_cl    <- -0.836;  label("Power exponent of albumin on baseline CL (unitless)")            # Wang 2024 Table 3: theta10 = -0.836
    e_alb_vc    <- -0.364;  label("Power exponent of albumin on Vc (unitless)")                     # Wang 2024 Table 3: theta11 = -0.364
    e_tumsz_cl  <-  0.0421; label("Power exponent of baseline tumour burden on baseline CL (unitless)") # Wang 2024 Table 3: theta23 = 0.0421

    # Categorical covariate effects on CL (exponential form: exp(theta * indicator)).
    # Wang 2024 Table 3 reports the back-transformed multiplicative factor exp(theta_k);
    # the natural-log of that factor is the additive coefficient on log(CL).
    e_sexf_cl   <- log(0.865); label("Exponential coefficient of female sex on baseline CL (unitless; log(exp(theta12)))") # Wang 2024 Table 3: exp(theta12) = 0.865
    e_ada_cl    <- log(1.11);  label("Exponential coefficient of ADA-positive on baseline CL (unitless; log(exp(theta14)))") # Wang 2024 Table 3: exp(theta14) = 1.11
    e_lymph_cl  <- log(0.877); label("Exponential coefficient of lymphoma tumour type (TTYPE1) on baseline CL (unitless; log(exp(theta15)))") # Wang 2024 Table 3: exp(theta15) = 0.877
    e_oth_cl    <- log(0.885); label("Exponential coefficient of 'other' tumour type (TTYPE3) on baseline CL (unitless; log(exp(theta16)))") # Wang 2024 Table 3: exp(theta16) = 0.885
    e_gc_cl     <- log(1.13);  label("Exponential coefficient of GCGEJ tumour type (TTYPE4) on baseline CL (unitless; log(exp(theta17)))") # Wang 2024 Table 3: exp(theta17) = 1.13
    e_escc_cl   <- log(0.99);  label("Exponential coefficient of ESCC tumour type (TTYPE5) on baseline CL (unitless; log(exp(theta18)))") # Wang 2024 Table 3: exp(theta18) = 0.99

    # Categorical covariate effects on Vc (exponential form).
    e_sexf_vc   <- log(0.866); label("Exponential coefficient of female sex on Vc (unitless; log(exp(theta13)))") # Wang 2024 Table 3: exp(theta13) = 0.866
    e_lymph_vc  <- log(0.879); label("Exponential coefficient of lymphoma tumour type (TTYPE1) on Vc (unitless; log(exp(theta19)))") # Wang 2024 Table 3: exp(theta19) = 0.879
    e_oth_vc    <- log(0.926); label("Exponential coefficient of 'other' tumour type (TTYPE3) on Vc (unitless; log(exp(theta20)))") # Wang 2024 Table 3: exp(theta20) = 0.926
    e_gc_vc     <- log(1.14);  label("Exponential coefficient of GCGEJ tumour type (TTYPE4) on Vc (unitless; log(exp(theta21)))") # Wang 2024 Table 3: exp(theta21) = 1.14
    e_escc_vc   <- log(1.08);  label("Exponential coefficient of ESCC tumour type (TTYPE5) on Vc (unitless; log(exp(theta22)))") # Wang 2024 Table 3: exp(theta22) = 1.08

    # IIV. CL and Vc form a 2x2 log-normal block; Vp and T50 are independent
    # log-normal etas; emax has an independent additive eta on the linear
    # scale (Wang 2024 Table 3 footnote: "emax_i = theta5 + eta_Emax,i").
    # Q has no IIV reported. Source CV%-to-omega^2 conversions (log-normal):
    #   omega^2 = log(1 + CV%^2)
    #   CL  19.5%  -> 0.0373
    #   Vc  15.5%  -> 0.0237
    #   Vp  68.5%  -> 0.3847
    #   T50 64.2%  -> 0.3451
    # Source covariance Cov(CL, Vc) = 0.0161 reported directly on the
    # omega-block scale; correlation = 0.0161 / sqrt(0.0373 * 0.0237) = 0.541.
    # emax additive-eta variance: source CV% 18.5% interpreted as
    # SD(eta_Emax) / |emax| -> SD = 0.185 * 0.528 = 0.0977,
    # variance = 0.00955.
    etalcl + etalvc ~ c(0.0373,
                        0.0161, 0.0237)  # Wang 2024 Table 3: IIV CL 19.5%, IIV Vc 15.5%, Cov(CL,Vc) 0.0161
    etalvp   ~ 0.3847                    # Wang 2024 Table 3: IIV Vp 68.5%
    etaemax  ~ 0.00955                   # Wang 2024 Table 3: IIV emax 18.5% (additive eta on linear-scale emax)
    etalt50  ~ 0.3451                    # Wang 2024 Table 3: IIV T50 64.2%

    # Residual error. Source residual model: ln(y_ij) = ln(yhat_ij) + eps_ij
    # with var(eps) = sigma^2; sigma reported as 17.9% in Wang 2024 Table 3.
    # NONMEM "additive on log-scale" residual maps to nlmixr2 prop() with
    # propSd = sigma on the linear-fraction scale.
    propSd <- 0.179; label("Proportional residual error (fraction)")  # Wang 2024 Table 3: Residual error 17.9%
  })
  model({
    # Individual baseline CL (CL at t = 0). Power form for continuous covariates
    # (cov / ref)^exponent; exponential form for categorical covariates
    # exp(theta_k * indicator_k). Reference covariates: WT 61 kg, ALB 41.5 g/L,
    # TUMSZ 47 mm, SEXF 0, ADA_POS 0, all TUMTP_* indicators 0 (NSCLC).
    cl0 <- exp(lcl + etalcl) *
      (WT    / 61)^e_wt_cl *
      (ALB   / 41.5)^e_alb_cl *
      (TUMSZ / 47)^e_tumsz_cl *
      exp(e_sexf_cl  * SEXF) *
      exp(e_ada_cl   * ADA_POS) *
      exp(e_lymph_cl * TUMTP_LYMPH) *
      exp(e_oth_cl   * TUMTP_OTH) *
      exp(e_gc_cl    * TUMTP_GC) *
      exp(e_escc_cl  * TUMTP_ESCC)

    vc <- exp(lvc + etalvc) *
      (WT  / 61)^e_wt_vc *
      (ALB / 41.5)^e_alb_vc *
      exp(e_sexf_vc  * SEXF) *
      exp(e_lymph_vc * TUMTP_LYMPH) *
      exp(e_oth_vc   * TUMTP_OTH) *
      exp(e_gc_vc    * TUMTP_GC) *
      exp(e_escc_vc  * TUMTP_ESCC)

    q  <- exp(lq)              # Q has no IIV in Wang 2024
    vp <- exp(lvp + etalvp)

    # Time-varying CL: sigmoidal emax function of time since first dose.
    t50    <- exp(lt50 + etalt50)
    lambda <- exp(llambda)
    emax_i <- emax + etaemax
    cl     <- cl0 * exp(emax_i * t^lambda / (t50^lambda + t^lambda))

    # Two-compartment micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L => central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
