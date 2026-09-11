Wang_2025_serplulimab <- function() {
  description <- "Two-compartment population PK model with sigmoidal time-varying clearance for intravenous serplulimab (anti-PD-1 IgG4) in adults with advanced solid tumours, including extensive-stage small cell lung cancer, across eight Phase I-III trials (Wang 2025)"
  reference <- "Wang K, Shen Y, Hu C, Xu F, Wang Q, Gao Y, Zhou L. Population Pharmacokinetics and Exposure-Response Analysis of Serplulimab in Small Cell Lung Cancer Patients. Clin Transl Sci. 2025;18(9):e70322. doi:10.1111/cts.70322"
  vignette <- "Wang_2025_serplulimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2025 Section 3.2 (the final
  # model is a two-compartment linear-elimination model fit to serum
  # serplulimab concentrations; Section 3.1: "6677 serplulimab serum
  # concentration measurements").
  compartmentData <- list(
    central     = list(analyte = "serplulimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "serplulimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL (exponent 0.531) and on Vc (exponent 0.450). Reference 65 kg, taken from the printed final-model equations in Wang 2025 Section 3.2 ('ln(WT/65)'); the PK-dataset median is 64.5 kg (Table 1), so 65 kg is the rounded population median. Wang 2025 Supporting Information notes body weight was tested first because it is correlated with BMI and BSA, and was retained in the base model on the strength of the OFV drop.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on baseline CL only (exponent -0.783); no albumin effect on Vc was retained. Reference 41.3 g/L, taken from the printed final-model equation in Wang 2025 Section 3.2 ('ln(ALB/41.3)'), which equals the PK-dataset median in Table 1 exactly. Source paper reports albumin in g/L (SI convention), matching the canonical unit. Higher albumin lowers CL and therefore raises exposure (Wang 2025 Section 3.3: exposure ratios 0.818-1.17 across albumin quartiles).",
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on Vc only (coefficient -0.121); no sex effect on CL was retained. Wang 2025 Section 3.2 defines 'SEX = 0 for male, SEX = 1 for female', which matches canonical SEXF directly with no value transformation. Female subjects have exp(-0.121) = 0.886 of the male Vc, giving slightly higher exposure in females (Wang 2025 Section 3.3: exposure ratios 0.945-1.07). PK dataset was 19.84% female (227/1144; Table 1).",
      source_name        = "SEX"
    ),
    TUMTP_NONLUNG = list(
      description        = "Pooled non-lung-cancer tumour-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (lung cancer -- NSCLC or SCLC)",
      notes              = "Exponential effect on Vc only (coefficient -0.0887); no tumour-type effect on CL was retained. Wang 2025 Section 3.2 defines 'TUMTP = 0 for lung cancer, TUMTP = 1 for hepatocellular carcinoma, colorectal cancer, and other tumor types', i.e. a single pooled non-lung indicator rather than one indicator per histology. In the PK dataset (Table 1) the reference lung-cancer group is 817/1144 (71.42%) and the pooled non-lung group is 327/1144 (28.58%: hepatic cancer 125, colorectal cancer 86, other 116). Wang 2025 Section 3.3 forest plot (Figure 1) reports the lung / non-lung contrast as exposure ratios 1.01 vs 0.973.",
      source_name        = "TUMTP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1144L,
    n_studies      = 8L,
    n_observations = 6650L,
    age_range      = "23.0-83.0 years",
    age_median     = "61.0 years",
    weight_range   = "33.0-131 kg",
    weight_median  = "64.5 kg",
    sex_female_pct = 19.84,
    race_ethnicity = c(Asian = 78.41, `Non-Asian` = 21.59),
    disease_state  = "Adults with advanced solid tumours. Tumour-type mix in the PK dataset (n = 1144): lung cancer 817 (71.42%; NSCLC and SCLC pooled), hepatic cancer 125 (10.93%), colorectal cancer 86 (7.52%), other 116 (10.14%). The exposure-response efficacy dataset is the ES-SCLC subset from the Phase III ASTRUM-005 trial (HLX10-005-SCLC301 / NCT04063163; n = 389).",
    dose_range     = "Serplulimab 0.3-10 mg/kg IV across the Phase I dose-escalation trials; the recommended Phase II/III dose is 3 mg/kg Q2W or 4.5 mg/kg Q3W. ASTRUM-005 used 4.5 mg/kg Q3W with a 1-h infusion. Flat doses of 200 mg and 300 mg were also studied in Phase I.",
    regions        = "Predominantly China (Asian 78.41%); the Phase III ASTRUM-005 trial was multinational, contributing the 21.59% non-Asian subjects.",
    ada_status     = "ADA-negative 1078 (94.23%); ADA-positive 23 (2.01%); missing 43 (3.76%). Geometric-mean CL was 13.3% higher in ADA-positive subjects, which Wang 2025 judged not clinically meaningful; ADA was not retained in the final model.",
    ecog_status    = "ECOG performance status 0 in 306 (26.75%), 1 in 835 (72.99%), 2 in 3 (0.26%).",
    albumin        = "41.3 g/L median (23.9-67.9 g/L range).",
    tumour_burden  = "88.0 mm median (10.0-350 mm range) in the PK dataset; 117 mm median (13.8-323 mm) in the ER efficacy dataset.",
    renal_function = "Creatinine clearance 90.9 mL/min median (28.5-291 mL/min); serum creatinine 68.9 umol/L median (23.0-156 umol/L). Neither was a significant covariate.",
    notes          = "Baseline demographics per Wang 2025 Table 1 (PK dataset column). Pooled dataset spans eight serplulimab (HLX10) trials: two Phase I (HLX10-001 / NCT03952403, HLX10HLX04-001 / NCT04818359), four Phase II (HLX10-008-HCC201 / NCT05246164, HLX10-010-MSI201 / NCT04747236, HLX10-011-CC201 / NCT03973112, HLX10HLX07-001 / NCT04297995) and two Phase III (HLX10-004-NSCLC303 / NCT04778904, HLX10-005-SCLC301 / NCT04063163). 6677 serum concentrations were collected; 27 below the LLOQ were excluded, leaving 6650 in the analysis. Two subjects (0.18%) with missing weight were excluded from the covariate-effect simulations only. Fit with NONMEM 7.5.0 using FOCE-I; PsN 4.2.0 for diagnostics and a 1000-replicate bootstrap."
  )

  # Covariates screened by Wang 2025 (Supporting Information Section 1) but NOT
  # retained in the final model. Wang 2025 Section 4 states these "did not have
  # a statistically significant effect on the PK parameters of serplulimab", so
  # no point estimate exists to encode. Documented here to preserve the
  # provenance of the covariate screen without a "declared but not referenced"
  # convention warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Wang 2025 Supporting Information Section 1, demographic category); not significant at the p = 0.01 forward / p = 0.001 backward criteria."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not statistically significant. Wang 2025 nevertheless carried race into the Figure 1 forest plot 'due to its potential clinical interest': non-Asian subjects (21.6%, N = 247) showed exposure ratios 0.945-1.17 vs Asian, which the authors note 'may be confounded by differences in body weight between the groups'. No point estimate is published, so the effect cannot be encoded."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody positivity status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Wang 2025 Section 4 reports geometric-mean CL 13.3% higher in the 23 ADA-positive subjects (2.01% incidence), judged not clinically meaningful. No model coefficient is published."
    ),
    TUMBUR = list(
      description = "Baseline tumour burden (sum of longest target-lesion diameters, RECIST)",
      units       = "mm",
      type        = "continuous",
      notes       = "Screened as a PK covariate and not retained. Tumour burden IS a significant predictor of overall survival in the Cox exposure-response model (Wang 2025 Table S4: beta = 0.4554 on log(TUMBUR), Wald p = 0.0149), but that is a survival covariate, not a PK covariate; the Cox sub-model is not encoded here (see the vignette Errata)."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a PK covariate and not retained. LDH is the strongest predictor of overall survival in the Cox exposure-response model (Wang 2025 Table S4: beta = 1.0558 on log(LDH), Wald p < 0.0010), but that is a survival covariate, not a PK covariate."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the laboratory-parameter category; not retained in the final PK model."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the laboratory-parameter category; not retained in the final PK model."
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened in the laboratory-parameter category; not retained in the final PK model."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened in the laboratory-parameter category; not retained in the final PK model."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the laboratory-parameter category; not retained in the final PK model. Consistent with a 148 kDa IgG4 mAb not being renally cleared."
    ),
    ECOG = list(
      description = "Eastern Cooperative Oncology Group performance status",
      units       = "(score)",
      type        = "categorical",
      notes       = "Screened in the disease/treatment category; not retained in the final PK model."
    ),
    CONMED_CHEMO = list(
      description = "Concomitant chemotherapy indicator (source column COMB)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the disease/treatment category; not retained. 762/1144 (66.61%) received concomitant chemotherapy."
    )
  )

  ini({
    # ---- Structural parameters (Wang 2025 Table 2, "Final model estimate") ----
    # Reference patient for the covariate model: WT 65 kg, ALB 41.3 g/L,
    # SEX = 0 (male), TUMTP = 0 (lung cancer), per the printed final-model
    # equations in Wang 2025 Section 3.2.
    lcl <- log(0.225); label("Baseline clearance CL0 at reference covariates (L/day)")   # Wang 2025 Table 2: CL0 = 0.225 L/day (RSE 2.85%, 95% CI 0.213-0.238)
    lvc <- log(3.52);  label("Central volume Vc at reference covariates (L)")            # Wang 2025 Table 2: Vc = 3.52 L (RSE 0.849%, 95% CI 3.46-3.58)
    lq  <- log(0.463); label("Intercompartmental clearance Q (L/day)")                   # Wang 2025 Table 2: Q = 0.463 L/day (RSE 5.98%, 95% CI 0.412-0.52)
    lvp <- log(2.21);  label("Peripheral volume Vp (L)")                                 # Wang 2025 Table 2: Vp = 2.21 L (RSE 7.41%, 95% CI 1.91-2.55)

    # ---- Time-varying clearance (Wang 2025 Section 3.2 equation 1; the
    # algebraically identical form is given in the Supporting Information
    # Section 1 as CL = CL0 * exp(Emax * (t/T50)^lambda / (1 + (t/T50)^lambda)):
    #   CL_i(t) = CL0 * exp(Emax_i * t^lambda / (T50^lambda + t^lambda) + <covariates> + eta_CL,i)
    # so CL(0) = CL0 and CL(t -> inf) = CL0 * exp(Emax). With Emax = -0.364 the
    # asymptotic CL is exp(-0.364) = 0.695 of baseline, i.e. a 30.5% reduction
    # at full saturation -- which is exactly the "exp(Emax) = 0.695" row of
    # Wang 2025 Table 2, an independent confirmation of the sign and scale.
    cl_time_max   <- -0.364;     label("Maximum log-scale change in CL from baseline (Emax; unitless, negative = CL decreases over time)") # Wang 2025 Section 3.2: Emax_i = -0.364 + eta_Emax,i; Table 2 reports the back-transform exp(Emax) = 0.695 (RSE 4.52%, 95% CI 0.636-0.759)
    lcl_t50       <- log(106);   label("log T50 - time at which half of the maximum CL change is reached (log days)")                       # Wang 2025 Table 2: T50 = 106 day (RSE 4.52%, 95% CI 83.2-129); the same 106 appears in the Section 3.2 equation
    lcl_time_hill <- log(2.05);  label("log lambda - sigmoidicity (Hill coefficient) of the time-on-CL function (log unitless)")            # Wang 2025 Table 2: lambda = 2.05 (RSE 14.2%, 95% CI 1.48-2.62); the same 2.05 is the printed exponent in the Section 3.2 equation

    # ---- Covariate effects (Wang 2025 Table 2 + Section 3.2 equations) ----
    # Continuous covariates enter as coefficients on ln(cov / reference) inside
    # the exp(), which is identical to a power form (cov / reference)^exponent.
    e_wt_cl       <-  0.531;  label("Power exponent of body weight on baseline CL (unitless)")  # Wang 2025 Table 2 CLwt = 0.531 (RSE 12.2%, 95% CI 0.403-0.658); Section 3.2 equation: + 0.531 * ln(WT / 65)
    e_alb_cl      <- -0.783;  label("Power exponent of albumin on baseline CL (unitless)")      # Wang 2025 Table 2 CLalb = -0.783 (RSE 13.4%, 95% CI -0.893 to -0.674); Section 3.2 equation: - 0.783 * ln(ALB / 41.3)
    e_wt_vc       <-  0.450;  label("Power exponent of body weight on Vc (unitless)")           # Wang 2025 Table 2 Vwt = 0.450 (RSE 7.3%, 95% CI 0.385-0.514); Section 3.2 equation: + 0.450 * ln(WT / 65)

    # Categorical covariates enter as exp(coefficient * indicator).
    e_sexf_vc     <- -0.121;  label("Exponential coefficient of female sex on Vc (unitless)")            # Wang 2025 Table 2 Vcsex = -0.121 (95% CI -0.152 to -0.0891); Section 3.2 equation: - 0.121 * SEX with SEX = 1 for female
    e_nonlung_vc  <- -0.0887; label("Exponential coefficient of non-lung tumour type on Vc (unitless)")  # Wang 2025 Table 2 Vctumtp = -0.0887 (RSE 15.5%, 95% CI -0.116 to -0.0617); Section 3.2 equation: - 0.0887 * TUMTP with TUMTP = 1 for non-lung

    # ---- Inter-individual variability (Wang 2025 Table 2) ----
    # Table 2 note: "IIV for CL, Vc, Q, Vp, Emax, and residual are reported as
    # approximate CV%". Each reported percentage is omega itself (the SD on the
    # estimation scale) x 100, NOT the variance. This is settled arithmetically
    # by the published 95% CI columns, which are symmetric on the VARIANCE
    # scale for every single row:
    #   CL    25.8% -> 0.258^2 = 0.06656; CI (23.5, 27.9) -> (0.05523, 0.07784), midpoint 0.06654
    #   Vc    15.4% -> 0.154^2 = 0.02372; CI (13.7, 17.0) -> (0.01877, 0.02890), midpoint 0.02383
    #   Q     49.7% -> 0.497^2 = 0.24701; CI (17.0, 68.1) -> (0.02890, 0.46376), midpoint 0.24633
    #   Vp    51.5% -> 0.515^2 = 0.26523; CI (41.7, 59.7) -> (0.17389, 0.35641), midpoint 0.26515
    #   Emax  26.3% -> 0.263^2 = 0.06917; CI (21.1, 30.6) -> (0.04452, 0.09364), midpoint 0.06908
    #   sigma 16.6% -> 0.166^2 = 0.02756; CI (15.7, 17.4) -> (0.02465, 0.03028), midpoint 0.02746
    # A variance reading would put every midpoint far from the point estimate,
    # so the SD reading is the only one consistent with the published CIs. The
    # same conversion applies to Emax even though its eta is additive rather
    # than log-normal -- the reporting routine emitted sqrt(omega^2) x 100 for
    # every row, so omega_Emax = 0.263 (NOT 0.263 x |Emax|).
    #
    # Cov(CL, Vc) = 0.017 is reported directly on the omega-block scale (its own
    # 95% CI 0.0128-0.0212 is symmetric about 0.017), giving a CL-Vc correlation
    # of 0.017 / (0.258 * 0.154) = 0.428.
    #
    # CL, Vc, Q and Vp carry log-normal etas (Supporting Information Section 1:
    # "Inter-individual variability (IIV) of PK parameters was estimated on a
    # logarithmic scale"). Emax carries an additive eta on the linear scale,
    # per the printed Section 3.2 relation Emax_i = -0.364 + eta_Emax,i.
    etalcl + etalvc ~ c(0.06656,
                        0.017, 0.02372)  # Wang 2025 Table 2: IIV CL 25.8%, IIV Vc 15.4%, Cov(CL, Vc) 0.017
    etalq          ~ 0.24701             # Wang 2025 Table 2: IIV Q 49.7%
    etalvp         ~ 0.26523             # Wang 2025 Table 2: IIV Vp 51.5%
    etacl_time_max ~ 0.06917             # Wang 2025 Table 2: IIV Emax 26.3% (additive eta on the linear-scale Emax)

    # ---- Residual error (Wang 2025 Table 2) ----
    # Supporting Information Section 1 prints the residual model as
    # ln(y_ij) = ln(yhat_ij) + eps_ij with var(eps) = sigma^2 ("Since the
    # concentration data were estimated on a logarithmic scale"). Additive
    # error on the log scale is proportional error in nlmixr2's linear space,
    # so propSd = sigma on the linear-fraction scale.
    propSd <- 0.166; label("Proportional residual error (fraction)")  # Wang 2025 Table 2: Residual sigma = 16.6% (RSE 2.74%, 95% CI 15.7-17.4)
  })
  model({
    # Individual baseline CL (CL at t = 0) and Vc. Wang 2025 Section 3.2 writes
    # every covariate term inside a single exp(), i.e.
    #   CL_i = 0.225 * exp(... + 0.531 * ln(WT/65) - 0.783 * ln(ALB/41.3) + eta_CL,i)
    #   Vc_i = 3.52  * exp(0.450 * ln(WT/65) - 0.0887 * TUMTP - 0.121 * SEX + eta_Vc,i)
    # exp(beta * ln(cov / ref)) is identically (cov / ref)^beta, so the
    # continuous terms are written in the equivalent power form below.
    cl0 <- exp(lcl + etalcl) *
      (WT  / 65)^e_wt_cl *
      (ALB / 41.3)^e_alb_cl

    vc <- exp(lvc + etalvc) *
      (WT / 65)^e_wt_vc *
      exp(e_sexf_vc    * SEXF) *
      exp(e_nonlung_vc * TUMTP_NONLUNG)

    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Time-varying CL: sigmoidal function of time since the first dose.
    # t^lambda / (T50^lambda + t^lambda) rises from 0 at t = 0 to 1 as
    # t >> T50, so CL falls from cl0 to cl0 * exp(Emax_i).
    cl_t50        <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl0 * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # Two-compartment micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L => central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
