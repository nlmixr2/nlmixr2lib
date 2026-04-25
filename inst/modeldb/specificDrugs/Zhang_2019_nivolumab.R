Zhang_2019_nivolumab <- function() {
  description <- "Two-compartment population PK model with time-varying clearance for intravenous nivolumab (anti-PD-1 IgG4) in adults with advanced solid tumors, alone or in combination with ipilimumab or chemotherapy (Zhang 2019)"
  reference <- "Zhang J, Sanghavi K, Shen J, et al. Population Pharmacokinetics of Nivolumab in Combination With Ipilimumab in Patients With Advanced Malignancies. CPT Pharmacometrics Syst Pharmacol. 2019;8(12):962-970. doi:10.1002/psp4.12476"
  vignette <- "Zhang_2019_nivolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, Q (with the CL exponent) and on Vc, Vp (with the Vc exponent). Reference 80 kg from Zhang 2019 Figure 1 caption (reference patient).",
      source_name        = "BBWT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 90 mL/min/1.73 m^2 (Zhang 2019 Figure 1 caption). The renormalization formula (MDRD vs CKD-EPI) is not explicitly stated in Zhang 2019 but is inherited from the Bajaj 2017 base model that this analysis re-estimated.",
      source_name        = "eGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on CL and on Vc; female reference value lower for both parameters. Source paper denotes the column SEX (1 = female, 0 = male per the Figure 1 caption naming the male reference patient).",
      source_name        = "SEX"
    ),
    ECOG_PS_GT0 = list(
      description        = "ECOG performance status > 0 indicator (binary collapse of the four-level scale)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG PS 0; fully active)",
      notes              = "Multiplicative effect on baseline CL (exp(0.181) = 1.198 fold higher CL when PS>0) and additive effect on the time-varying-CL Emax parameter. Zhang 2019 uses the binary collapse PS=0 vs. PS>0 rather than the full ECOG scale.",
      source_name        = "PS"
    ),
    RACE_BLACK = list(
      description        = "Race indicator: African American",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (white or other; the Zhang 2019 reference combines white and 'other' into one reference group)",
      notes              = "Exponential effect on CL (CL_RAAA = 0.0374; not statistically significant in Zhang 2019, 95% CI -0.0308 to 0.111). Source column name RAAA, defined in the Zhang 2019 Table 2 footnote as 'African American race'. Renamed to canonical RACE_BLACK per covariate-columns.md.",
      source_name        = "RAAA"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator: Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (white or other; same composite reference as RACE_BLACK)",
      notes              = "Exponential effect on CL (CL_RAAS = -0.0354 = 3.5% lower CL than the reference group). Source column name RAAS, defined in the Zhang 2019 Table 2 footnote as 'Asian race'. Renamed to canonical RACE_ASIAN per covariate-columns.md.",
      source_name        = "RAAS"
    ),
    COADMIN_IPI_3Q3W = list(
      description        = "Coadministration regimen: nivolumab + ipilimumab 3 mg/kg every 3 weeks (4-dose induction)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-3Q3W regimen)",
      notes              = "Exponential effect on baseline CL (exp(0.227) = 1.255 fold higher CL than monotherapy). Per Zhang 2019 Methods, the unmodeled ipilimumab schedules (1 mg/kg q3w x 4 induction; 1 mg/kg q12w) showed no statistically significant effect on CL and were collapsed into the reference group along with monotherapy.",
      source_name        = "IPI3Q3W"
    ),
    COADMIN_IPI_1Q6W = list(
      description        = "Coadministration regimen: nivolumab + ipilimumab 1 mg/kg every 6 weeks (continuous maintenance)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-1Q6W regimen)",
      notes              = "Exponential effect on baseline CL (exp(0.159) = 1.172 fold higher CL than monotherapy).",
      source_name        = "IPI1Q6W"
    ),
    COADMIN_CHEMO = list(
      description        = "Coadministration regimen: nivolumab + platinum-based chemotherapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no chemotherapy coadministration)",
      notes              = "Exponential effect on baseline CL (exp(-0.104) = 0.901 fold, ~9.7% lower than monotherapy). Pooled across the four chemotherapy backbones (gemcitabine + cisplatin, pemetrexed + cisplatin, paclitaxel + carboplatin, platinum doublet).",
      source_name        = "CHEMO"
    ),
    COADMIN_IPI_ANY = list(
      description        = "Any-ipilimumab-coadministration indicator (regimen-agnostic)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ipilimumab coadministration)",
      notes              = "Additive effect on the time-varying-CL Emax parameter (Emax_IPICO = -0.0668; subjects on any nivolumab + ipilimumab combination show an additional 6.7% reduction of CL at full saturation relative to monotherapy or chemotherapy combination). Source column name IPICO. Coexists with COADMIN_IPI_3Q3W / COADMIN_IPI_1Q6W because those act on baseline CL rather than on Emax.",
      source_name        = "IPICO"
    )
  )

  population <- list(
    n_subjects     = 6468L,
    n_studies      = 25L,
    n_observations = 32835L,
    age_range      = "Adults with advanced solid tumors (specific age summaries not tabulated in the published main text)",
    weight_range   = "47.7-122.0 kg (5th-95th percentiles)",
    weight_mean    = "77.6 kg (SD 18.8 kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Composite reference group is white/other; Asian (RAAS) and African American (RAAA) coded as separate non-reference indicators. Per-category percentages not tabulated in the published main text.",
    disease_state  = "Advanced solid tumors: NSCLC 38.25%, melanoma 26.93%, RCC 19.25%, SCLC 6.03%, HCC 5.89%, CRC 3.65%.",
    dose_range     = "Nivolumab 0.1-10 mg/kg or 240/480 mg flat IV every 2-4 weeks. Validated dosing regimens for VPC: 3 mg/kg or 240 mg q2w monotherapy; 3 mg/kg q2w + ipilimumab 1 mg/kg q6w; 3 mg/kg + ipilimumab 1 mg/kg q3w x 4 followed by 3 mg/kg q2w; 1 mg/kg + ipilimumab 3 mg/kg q3w x 4 followed by 3 mg/kg q2w.",
    regions        = "Pooled global trials (7 phase I, 2 phase I/II, 6 phase II, 9 phase III, 1 phase IIIb/IV)",
    coadministration = "Monotherapy 55.12%; ipilimumab 1 mg/kg q12w 0.56%; ipilimumab 1 mg/kg q6w 11.75%; ipilimumab 1 mg/kg q3w x 4 induction 15.06%; ipilimumab 3 mg/kg q3w x 4 induction 13.84%; chemotherapy 3.68%.",
    performance_status = "ECOG PS 0 47.02%; PS 1 51.27%; PS 2 1.62%; PS 3 0.02%; missing 0.08%.",
    notes          = "Baseline demographics per Zhang 2019 Table 1. The PopPK model re-estimates the Bajaj 2017 nivolumab monotherapy base model (reference 7 in Zhang 2019) on a 6,468-subject pooled dataset that adds combination-therapy cohorts and several additional tumor types."
  )

  ini({
    # Structural parameters - reference patient: 80 kg, eGFR 90 mL/min/1.73 m^2,
    # male, PS 0, white/other race, monotherapy, NSCLC (Zhang 2019 Figure 1 caption).
    # Source values reported in mL/hour and hour; converted to L/day and days
    # so the time axis matches the rest of nlmixr2lib's mAb library.
    lcl     <- log(10.8 * 24 / 1000); label("Baseline clearance CL0 at reference covariates (L/day)")     # Zhang 2019 Table 2: CL0 REF = 10.8 mL/hour
    lvc     <- log(4.27);             label("Central volume Vc at reference covariates (L)")              # Zhang 2019 Table 2: VC REF = 4.27 L
    lq      <- log(34.9 * 24 / 1000); label("Intercompartmental clearance Q at reference covariates (L/day)") # Zhang 2019 Table 2: Q REF = 34.9 mL/hour
    lvp     <- log(2.70);             label("Peripheral volume Vp at reference covariates (L)")           # Zhang 2019 Table 2: VP REF = 2.70 L

    # Time-varying CL Hill-Emax function (Zhang 2019 final-model equation):
    #   CL(t) = CL0 * exp(Emax * t^HILL / (T50^HILL + t^HILL))
    # T50 reported in hours; converted to days. HILL is unitless. Emax is on
    # the linear scale (negative = CL decreases over time).
    Emax     <-      -0.240;            label("Reference Emax of time-varying CL (Zhang 2019 'Emax REF'; unitless, negative = CL decreases at full saturation)") # Zhang 2019 Table 2: Emax REF = -0.240
    lt50     <- log(2200 / 24);         label("log T50 - time at which half of Emax is reached (log days)") # Zhang 2019 Table 2: T50 = 2200 hour
    lhill    <- log(2.77);              label("log HILL - sigmoidicity of the time-on-CL function (log unitless)") # Zhang 2019 Table 2: HILL = 2.77

    # Continuous covariate effects on CL and Vc (power form on log-transformed covariates,
    # equivalent to (cov/ref)^exponent in linear space; Zhang 2019 final-model equation).
    e_wt_cl    <-  0.530; label("Power exponent of body weight on CL (unitless)")  # Zhang 2019 Table 2: CL BBWT = 0.530
    e_egfr_cl  <-  0.202; label("Power exponent of eGFR on CL (unitless)")         # Zhang 2019 Table 2: CL eGFR = 0.202
    e_wt_vc    <-  0.534; label("Power exponent of body weight on Vc (unitless)")  # Zhang 2019 Table 2: VC BBWT = 0.534

    # Categorical covariate effects on CL (exponential form: exp(theta * indicator);
    # indicators are 0 / 1).
    e_sexf_cl       <- -0.181; label("Exponential coefficient of female sex on CL (unitless)")                # Zhang 2019 Table 2: CL SEX = -0.181
    e_ecog_cl       <-  0.181; label("Exponential coefficient of ECOG PS > 0 on CL (unitless)")               # Zhang 2019 Table 2: CL PS = 0.181
    e_black_cl      <-  0.0374; label("Exponential coefficient of African American race on CL (unitless)")    # Zhang 2019 Table 2: CL RAAA = 0.0374 (not statistically significant)
    e_asian_cl      <- -0.0354; label("Exponential coefficient of Asian race on CL (unitless)")               # Zhang 2019 Table 2: CL RAAS = -0.0354
    e_ipi3q3w_cl    <-  0.227; label("Exponential coefficient of ipilimumab 3 mg/kg q3w coadministration on CL (unitless)") # Zhang 2019 Table 2: CL IPI3Q3W = 0.227
    e_ipi1q6w_cl    <-  0.159; label("Exponential coefficient of ipilimumab 1 mg/kg q6w coadministration on CL (unitless)") # Zhang 2019 Table 2: CL IPI1Q6W = 0.159
    e_chemo_cl      <- -0.104; label("Exponential coefficient of chemotherapy coadministration on CL (unitless)")           # Zhang 2019 Table 2: CL CHEMO = -0.104

    # Categorical covariate effects on Vc (exponential form).
    e_sexf_vc <- -0.161; label("Exponential coefficient of female sex on Vc (unitless)") # Zhang 2019 Table 2: VC SEX = -0.161

    # Additive covariate effects on Emax (Zhang 2019 Emax equation:
    #   Emax_i = Emax_REF + Emax_PS * I[PS>0] + Emax_IPICO * I[IPICO=1] + eta_Emax_i).
    e_ecog_emax  <- -0.138;  label("Additive effect of ECOG PS > 0 on Emax (unitless)")          # Zhang 2019 Table 2: Emax PS = -0.138
    e_ipico_emax <- -0.0668; label("Additive effect of any ipilimumab coadministration on Emax (unitless)") # Zhang 2019 Table 2: Emax IPICO = -0.0668

    # IIV. CL and Vc form a 2x2 log-normal block; Q and Vp are independent etas
    # constrained to share the variance of CL and Vc respectively per the
    # Methods ("the IIV random effect of Q follows the same distribution as
    # that of CL ... the IIV random effect of VP follows the same distribution
    # as that of VC"; Q and Vp do not appear in the Table 2 random-effects
    # block because their variances are pinned, not separately estimated).
    # Emax has an independent additive (not log-normal) eta on the linear scale.
    etalcl + etalvc ~ c(0.157,
                        0.0596, 0.152)  # Zhang 2019 Table 2: omega^2_CL, cov(CL,VC), omega^2_VC
    etalq    ~ 0.157   # constrained equal to omega^2_CL (Zhang 2019 Methods)
    etalvp   ~ 0.152   # constrained equal to omega^2_VC (Zhang 2019 Methods)
    etaEmax  ~ 0.0874  # Zhang 2019 Table 2: omega^2_Emax = 0.0874 (additive eta on linear-scale Emax)

    # Residual error (proportional only). Source value 0.245 is the standard
    # deviation on the linear scale per the table footnote naming the column
    # "Estimate" (not variance) and the proportional residual error model.
    propSd <- 0.245; label("Proportional residual error (fraction)")  # Zhang 2019 Table 2: Proportional = 0.245
  })
  model({
    # Individual baseline CL (CL at t=0). All categorical effects are exponential
    # of the form exp(theta * indicator); continuous effects are power form
    # (cov / ref)^exponent. Reference covariates: WT 80 kg, eGFR 90 mL/min/1.73 m^2,
    # male (SEXF=0), PS 0 (ECOG_PS_GT0=0), white/other race (RACE_BLACK=0 and
    # RACE_ASIAN=0), monotherapy (all COADMIN_*=0).
    cl0 <- exp(lcl + etalcl) *
      (WT   / 80)^e_wt_cl *
      (CRCL / 90)^e_egfr_cl *
      exp(e_sexf_cl    * SEXF) *
      exp(e_ecog_cl    * ECOG_PS_GT0) *
      exp(e_black_cl   * RACE_BLACK) *
      exp(e_asian_cl   * RACE_ASIAN) *
      exp(e_ipi3q3w_cl * COADMIN_IPI_3Q3W) *
      exp(e_ipi1q6w_cl * COADMIN_IPI_1Q6W) *
      exp(e_chemo_cl   * COADMIN_CHEMO)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sexf_vc * SEXF)

    # Q and Vp inherit the body-weight scaling of CL and Vc respectively
    # (Zhang 2019 Methods: "the effect of BBWT was added on Q and VP, and their
    # estimates were fixed to be similar to those of CL and VC").
    q  <- exp(lq  + etalq)  * (WT / 80)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 80)^e_wt_vc

    # Time-varying CL: Hill-Emax function of time since first dose.
    t50  <- exp(lt50)
    hill <- exp(lhill)
    Emax_i <- Emax + e_ecog_emax * ECOG_PS_GT0 + e_ipico_emax * COADMIN_IPI_ANY + etaEmax
    cl <- cl0 * exp(Emax_i * t^hill / (t50^hill + t^hill))

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
