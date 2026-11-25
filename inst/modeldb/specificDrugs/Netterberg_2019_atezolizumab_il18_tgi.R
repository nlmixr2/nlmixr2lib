Netterberg_2019_atezolizumab_il18_tgi <- function() {
  description <- paste(
    "Integrated PK / biomarker / tumor-size model for atezolizumab in",
    "relapsed or refractory non-small cell lung cancer (Netterberg 2019,",
    "phase I study PCD4989g). A two-compartment PK model with linear",
    "elimination drives an effect compartment that linearly stimulates the",
    "release of interleukin-18 (IL-18) from a depletable precursor pool into",
    "plasma (pool indirect-response model). Tumor size (sum of the longest",
    "diameters, SLD) follows a Claret-type model with zero-order growth and two",
    "first-order shrinkage terms: one proportional to the cycle-specific",
    "exposure AUC = dose / CL that wanes exponentially with time since the",
    "first dose, and one proportional to the individual relative change from",
    "baseline in IL-18 at day 21 (zero before day 21). A constant-hazard",
    "dropout model with the predicted relative change in tumor size as",
    "predictor is included for simulating dropout. Covariates: more than two",
    "metastatic sites (SLD0), liver metastases (growth rate) and never-smoking",
    "status (no AUC-driven shrinkage). Time is in days; the tumor-model",
    "parameters keep the paper's per-week units and are converted in model().",
    "Companion model: modellib('Netterberg_2019_atezolizumab_itac').",
    sep = " "
  )
  reference <- paste(
    "Netterberg I, Li CC, Molinero L, Budha N, Sukumaran S, Stroh M,",
    "Jonsson EN, Friberg LE. A PK/PD Analysis of Circulating Biomarkers and",
    "Their Relationship to Tumor Response in Atezolizumab-Treated non-small",
    "Cell Lung Cancer Patients. Clin Pharmacol Ther. 2019;105(2):486-495.",
    "doi:10.1002/cpt.1198.",
    "PK parameters are in Supplementary Table S3; biomarker, tumor and",
    "dropout parameters in Table 1.",
    sep = " "
  )
  vignette <- "Netterberg_2019_atezolizumab"

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # The IL-18 precursor pool, the IL-18 plasma state and the bookkeeping state
  # that holds the IL-18 value reached at day 21 are paper-specific states.
  paper_specific_compartments <- c("il18_pool", "il18", "il18_d21")

  covariateData <- list(
    DOSE_ATEZOLIZUMAB_MG = list(
      description = "Atezolizumab dose administered in the current 3-week cycle (mg); 0 before the first dose.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: carries the amount of the most recent dose, so it must agree with the dosing records' amt. Enters the cycle-specific exposure AUCcycle,n = Dose_n / CL_i (Eq. 2) that drives AUC-related tumor shrinkage. For mg/kg dosing it is the mg/kg dose times body weight. Median AUC in the analysis dataset was 5,600 mg*day/L (5th-95th percentile 2,900-9,500; Figure 4).",
      source_name = "Dose_n"
    ),
    MET_GE3 = list(
      description = "Baseline number of metastatic sites >= 3 indicator (1 = more than two metastatic sites at baseline, 0 = two or fewer).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (two or fewer metastatic sites at baseline)",
      notes = "Time-fixed. Multiplies the typical baseline SLD by (1 + 0.927): 4.75 cm vs 9.15 cm (Netterberg 2019 Results: 'For patients with two or fewer metastases at baseline, the typical SLD0 was 4.8 cm, whereas for patients with more than two metastases, the corresponding value was 9.1 cm'). Derived from the paper's integer count BMETS (median 2, range 0-4; Supplementary Table S5) as as.integer(BMETS > 2).",
      source_name = "BMETS"
    ),
    LMET = list(
      description = "Liver metastases at baseline (1 = present, 0 = absent).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no liver metastases at baseline)",
      notes = "Time-fixed. Multiplies the zero-order tumor growth rate by (1 + 8.48): 0.017 vs 0.16 cm/week (Netterberg 2019 Results). 26 of 88 patients (30%) had liver metastases (Supplementary Table S5).",
      source_name = "BLMETS"
    ),
    SMOKE_NEVER = list(
      description = "Never-smoker indicator (1 = never smoked, 0 = current or former smoker).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (current or former smoker)",
      notes = "Time-fixed, used unpaired (never vs ever). The AUC-related shrinkage rate constant KShr,Drug applies to current / former smokers; in never smokers it is fixed to 0 (Table 1 row 'theta KShr,Drug-Non-Smoking = 0 fixed'; Results: 'KShr,Drug approached zero in patients who had never smoked and was set to 0 in the final model'). 17 of 88 patients (19%) never smoked (Supplementary Table S5).",
      source_name = "SMOKE"
    )
  )

  covariatesDataExcluded <- list(
    PDL1_IC3_TC3 = list(
      description = "PD-L1-high indicator (immunohistochemistry IC3 or TC3) screened on the tumor-model parameters.",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search (Supplementary Table S6) but not retained: 'No statistically significant relationship was found between baseline expression of PD-L1 ... on any of the model parameters' (Netterberg 2019 Results). Baseline SLD, ECOG, histology, EGFR mutation, neutrophil/lymphocyte ratio and tumor gene expression (IL8, MICA, MICB, CCL2, IL1b) were screened likewise and not retained."
    )
  )

  compartmentData <- list(
    central = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "atezolizumab", units = "mg", specimen = "not applicable", verified = TRUE),
    effect = list(
      analyte = "atezolizumab (effect-compartment concentration)",
      units = "ug/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    il18_pool = list(
      analyte = "interleukin-18 (precursor pool)",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    il18 = list(analyte = "interleukin-18", units = "pg/mL", specimen = "plasma", verified = TRUE),
    il18_d21 = list(
      analyte = "interleukin-18 (value held from day 21)",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_size = list(
      analyte = "tumor size (sum of the longest diameters of RECIST target lesions)",
      units = "cm",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative dropout hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 1L,
    disease_state = "relapsed / refractory non-small cell lung cancer",
    dose_range = "atezolizumab 10, 15 or 20 mg/kg or 1,200 mg flat dose as a 30-60 min IV infusion every 3 weeks (one patient received 1 mg/kg)",
    regions = "not reported (phase Ia study PCD4989g, NCT01375842)",
    n_observations = "IL-18: 450 above LOQ; tumor SLD assessments: 485 (Supplementary Table S2)",
    covariate_summary = "More than two metastatic sites: count median 2 (range 0-4); liver metastases 30%; never smokers 19%, current 11%, former 69%; ECOG 0 28%, ECOG 1 72%; non-squamous histology 76% (Supplementary Table S5)",
    notes = "The PK model was developed on 167 patients from PCD4989g who received 10, 15 or 20 mg/kg every 3 weeks; its individual estimates were used for the 88 NSCLC patients analysed for biomarkers and tumor size. Median tumor-assessment duration was 18 weeks (range 0-153). Age, weight, sex and race are not tabulated in Netterberg 2019."
  )

  ini({
    # Two-compartment PK (Supplementary Table S3; the same values are quoted in
    # the Figure S2 legend).
    lcl <- log(0.233)
    label("Clearance CL (L/day)") # Supplementary Table S3: CL = 0.233 L/day (RSE 4.5%)
    lvc <- log(3.35)
    label("Central volume Vc (L)") # Supplementary Table S3: Vc = 3.35 L (RSE 2.3%)
    lq <- log(0.452)
    label("Intercompartmental clearance Q (L/day)") # Supplementary Table S3: Q = 0.452 L/day (RSE 8.6%)
    lvp <- log(3.36)
    label("Peripheral volume Vp (L)") # Supplementary Table S3: Vp = 3.36 L (RSE 7.8%)

    # IL-18 pool indirect-response model (Table 1, IL-18 model; Eq. 3-4)
    lkin_il18 <- log(23.0)
    label("IL-18 zero-order input rate Rin (pg/mL/day)") # Table 1: Rin = 23.0 pg/mL/day (RSE 12%)
    lkout_il18 <- log(0.106)
    label("IL-18 first-order turnover rate constant kout (1/day)") # Table 1: kout = 0.106 1/day (RSE 13%)
    lbl_il18_pool <- log(499)
    label("Baseline IL-18 concentration in the precursor pool Pool0 (pg/mL)") # Table 1: Pool0 = 499 pg/mL (RSE 18%)
    lslope_il18 <- log(0.168)
    label("Linear effect of effect-compartment atezolizumab on IL-18 release (mL/ug)") # Table 1: Slope = 0.168 mL/ug (RSE 26%)
    lke0 <- log(4.36e-3)
    label("Effect-compartment rate constant ke0 = ke1 (1/day)") # Table 1: ke0 = 4.36e-3 1/day (RSE 25%)

    # Tumor growth inhibition model (Table 1, TGI model; Eq. 1). Rates are in
    # the paper's per-week units.
    lrbase_tumor_size <- log(4.75)
    label("Baseline tumor size SLD0, two or fewer metastatic sites (cm)") # Table 1: SLD0 = 4.75 cm (RSE 11%)
    lrgrowth <- log(0.0169)
    label("Zero-order tumor growth rate RGrowth, no liver metastases (cm/week)") # Table 1: RGrowth = 0.0169 cm/week (RSE 17%)
    lkshr_drug <- log(2.75e-6)
    label("Tumor shrinkage rate constant per unit cycle AUC, current/former smokers (L/mg/day per week)") # Table 1: KShr,Drug = 2.75e-6 (RSE 23%)
    llambda_drug <- log(0.0609)
    label("Rate constant for the waning AUC effect on tumor shrinkage lambdaDrug (1/week)") # Table 1: lambdaDrug = 0.0609 1/week (RSE 8.7%)
    lkshr_il18 <- log(0.00282)
    label("Tumor shrinkage rate constant per unit IL-18 relative change at day 21 (1/week)") # Table 1: KShr,IL-18 = 0.00282 (RSE 14.1%)
    e_lmet_rgrowth <- 8.48
    label("Fractional change in RGrowth with liver metastases at baseline (1 + theta)") # Table 1: theta Rgrowth-LiverMets = 8.48 (RSE 8.5%)
    e_met_ge3_rbase_tumor_size <- 0.927
    label("Fractional change in SLD0 with more than two metastatic sites (1 + theta)") # Table 1: theta SLD0-Mets = 0.927 (RSE 27%)
    kshr_drug_smoke_never <- fixed(0)
    label("AUC-related tumor shrinkage rate constant in never smokers (L/mg/day per week)") # Table 1: theta KShr,Drug-Non-Smoking = 0 fixed

    # Dropout model (Table 1, dropout model; Eq. 5)
    lhaz_base <- log(0.0306)
    label("Constant baseline dropout hazard h0 (1/week)") # Table 1: h0 = 0.0306 1/week (RSE 1.8%)
    e_tumor_size_rcfb_hazard <- 1.56
    label("Coefficient of the predicted relative change from baseline tumor size on the log dropout hazard (unitless)") # Table 1: beta = 1.56 (RSE 8.8%)

    # IIV: Table 1 / Table S3 report IIV as CV% of log-normal etas; variances
    # are omega^2 = log(1 + CV^2).
    etalcl ~ 0.162430 # Supplementary Table S3: IIV CL 42% -> log(1 + 0.42^2)
    etalvc ~ 0.070363 # Supplementary Table S3: IIV Vc 27% -> log(1 + 0.27^2)
    etalkout_il18 ~ 0.199590 # Table 1: IIV kout (IL-18) 47% -> log(1 + 0.47^2)
    etalbl_il18_pool ~ 0.475240 # Table 1: IIV Pool0 (IL-18) 78% -> log(1 + 0.78^2)
    etalslope_il18 ~ 0.494696 # Table 1: IIV Slope (IL-18) 80% -> log(1 + 0.80^2)
    etalrbase_tumor_size ~ 0.325266 # Table 1: IIV SLD0 62% -> log(1 + 0.62^2)
    etalrgrowth ~ 0.623207 # Table 1: IIV RGrowth 93% -> log(1 + 0.93^2)
    etalkshr_drug ~ 0.931219 # Table 1: IIV KShr,Drug 124% -> log(1 + 1.24^2)
    etalkshr_il18 ~ 1.453117 # Table 1: IIV KShr,IL-18 181% -> log(1 + 1.81^2)

    # Residual error: log-transformed data with additive error on the log
    # scale, encoded as proportional error.
    propSd <- 0.187
    label("Proportional residual error, atezolizumab serum concentration (fraction)") # Supplementary Table S3: residual error 18.7% (RSE 6.6%)
    propSd_il18 <- 0.223
    label("Proportional residual error, plasma IL-18 (fraction)") # Table 1: IL-18 residual error 22.3% (RSE 2.5%)
    propSd_tumor_size <- 0.134
    label("Proportional residual error, tumor size SLD (fraction)") # Table 1: TGI residual error 13.4% (RSE 2.1%)
  })
  model({
    # Individual PK parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual IL-18 parameters. Baseline plasma IL-18 is Rin / kout and the
    # pool-to-plasma release rate constant is krel = kout * BioM0 / Pool0
    # (Eq. 3), i.e. krel = Rin / Pool0.
    kin_il18 <- exp(lkin_il18)
    kout_il18 <- exp(lkout_il18 + etalkout_il18)
    bl_il18_pool <- exp(lbl_il18_pool + etalbl_il18_pool)
    slope_il18 <- exp(lslope_il18 + etalslope_il18)
    ke0 <- exp(lke0)
    bl_il18 <- kin_il18 / kout_il18
    krel_il18 <- kout_il18 * bl_il18 / bl_il18_pool

    # Individual tumor-model parameters (per week)
    rbase_tumor_size <- exp(lrbase_tumor_size + etalrbase_tumor_size) *
      (1 + e_met_ge3_rbase_tumor_size * MET_GE3)
    rgrowth <- exp(lrgrowth + etalrgrowth) * (1 + e_lmet_rgrowth * LMET)
    kshr_drug <- (1 - SMOKE_NEVER) * exp(lkshr_drug + etalkshr_drug) +
      SMOKE_NEVER * kshr_drug_smoke_never
    lambda_drug <- exp(llambda_drug)
    kshr_il18 <- exp(lkshr_il18 + etalkshr_il18)
    haz_base <- exp(lhaz_base)

    # PK and effect compartment (Eq. 4, ke1 = ke0; the effect compartment
    # starts at zero)
    Cc <- central / vc
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(effect) <- ke0 * Cc - ke0 * effect

    # IL-18 pool indirect-response model (Eq. 3) with linear drug effect
    eff_il18 <- slope_il18 * effect
    rate_il18 <- krel_il18 * il18_pool * (1 + eff_il18) - kout_il18 * il18
    d/dt(il18_pool) <- kin_il18 - krel_il18 * il18_pool * (1 + eff_il18)
    d/dt(il18) <- rate_il18
    il18_pool(0) <- bl_il18_pool
    il18(0) <- bl_il18

    # il18_d21 follows IL-18 until day 21 and then holds the day-21 value, so
    # the relative change from baseline at day 21 is available to the tumor
    # model; it is set to zero before day 21 (Methods, Tumor growth inhibition
    # model).
    d/dt(il18_d21) <- (t < 21) * rate_il18
    il18_d21(0) <- bl_il18
    rcfb_il18_d21 <- (t >= 21) * (il18_d21 - bl_il18) / bl_il18

    # Cycle-specific exposure AUCcycle,n = Dose_n / CL_i (Eq. 2), in mg*day/L.
    # The cycle dose is a data column rather than podo(central): inside the
    # ODE right-hand side podo() returns the infusion rate, not the amount,
    # for an infused dose.
    auc_cycle <- DOSE_ATEZOLIZUMAB_MG / cl

    # Tumor size (Eq. 1). Time since the first dose (t, days) is converted to
    # weeks in the waning AUC effect; per-week rates are divided by 7 to give
    # per-day rates.
    d/dt(tumor_size) <- rgrowth / 7 -
      (kshr_drug * auc_cycle * exp(-lambda_drug * t / 7) +
        kshr_il18 * rcfb_il18_d21) / 7 * tumor_size
    tumor_size(0) <- rbase_tumor_size

    # Dropout (Eq. 5): constant baseline hazard and the predicted relative
    # change from baseline tumor size as predictor.
    tumor_size_rcfb <- (tumor_size - rbase_tumor_size) / rbase_tumor_size
    hazard <- haz_base / 7 * exp(e_tumor_size_rcfb_hazard * tumor_size_rcfb)
    d/dt(cumhaz) <- hazard
    surv_dropout <- exp(-cumhaz)

    Cc ~ prop(propSd)
    il18 ~ prop(propSd_il18)
    tumor_size ~ prop(propSd_tumor_size)
  })
}
