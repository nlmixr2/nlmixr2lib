Gastonguay_2005_efaproxiral <- function() {
  description <- "Linear 2-compartment IV popPK model for efaproxiral (RSR13) with an algebraic linear RBC:plasma proportionality (Crbc = SLPRBC * Cc) and a linear PD model relating partial pressure of oxygen at 50% hemoglobin saturation (p50, mmHg) to RBC efaproxiral concentration (p50 = INTp50 + SLPp50 * Crbc). Population of 451 cancer patients receiving radiation therapy pooled across six phase I-III trials. Full covariate model: BSA and AGE on CL; BSA, AGE and baseline albumin (BALB) on V1 and V2; BSA on Q; max administered dose (MDOS), AGE and BALB on SLPRBC; primary cancer type indicators (breast, glioma, other; reference = lung) on SLPp50. Three-way correlated IIV block on CL, V1, V2 plus diagonal IIV on Q, SLPRBC, INTp50 and SLPp50."
  reference <- paste(
    "Gastonguay MR, Venitz J, Steffen RP, Hackman J.",
    "Population pharmacokinetic-pharmacodynamic modeling of",
    "efaproxiral in cancer patients receiving radiation therapy.",
    "Clin Pharmacol Ther. 2005;77(2):P55 (PII-101).",
    "ASCPT 2005 annual meeting poster.",
    "doi:10.1016/j.clpt.2004.12.233.",
    sep = " "
  )
  vignette <- "Gastonguay_2005_efaproxiral"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area (baseline, per subject; m^2).",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 1.8 m^2 per Equation 1 of Gastonguay 2005; cohort median 1.82 m^2 (Table 1). Used as power-model effect on CL, V1, V2 and Q.",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Age at study entry (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 60 yrs per Equation 1; cohort median 57 yrs (Table 1). Used as power-model effect on CL, V1, V2 and SLPRBC.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin. Source paper reports baseline albumin (BALB) in g/dL with reference value 3.5 g/dL in Equation 1; canonical ALB is in g/L (SI), so model() applies an inline conversion alb_gdL <- ALB * 0.1 before the power-model power-of-(BALB/3.5) terms.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BALB in g/dL (Table 1 cohort median 3.6 g/dL). Reference value 3.5 g/dL per Equation 1 of Gastonguay 2005. Power-model effect on V1, V2, SLPRBC. Conversion factor 1 g/dL = 10 g/L applied in model().",
      source_name        = "BALB"
    ),
    DOSE_EFP_MAX_MG = list(
      description        = "Per-subject maximum administered efaproxiral dose (mg). Time-fixed; defined as the largest single-administration dose the subject received during the trial.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 6800 mg per Equation 1 of Gastonguay 2005 (close to median 100 mg/kg x ~72 kg). Used as power-model effect on SLPRBC. The paper uses MDOS as the variable symbol.",
      source_name        = "MDOS"
    ),
    TUMTP_BREAST = list(
      description        = "Breast-cancer indicator (1 = primary cancer type breast, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-breast; lung cancer is the implicit reference among the four-level CATP factor)",
      notes              = "Source factor CATP with levels 1=lung (reference), 2=breast, 3=cranial GBM, 4=other (Data section). Decompose as TUMTP_BREAST = as.integer(CATP == 2). Used as a categorical power-model multiplier on SLPp50.",
      source_name        = "(derived from CATP == 2)"
    ),
    TUMTP_GLIO = list(
      description        = "Cranial glioblastoma multiforme indicator (1 = primary cancer type cranial GBM, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-GBM; lung cancer is the implicit reference)",
      notes              = "Source factor CATP level 3 = cranial GBM. Decompose as TUMTP_GLIO = as.integer(CATP == 3). Used as a categorical power-model multiplier on SLPp50.",
      source_name        = "(derived from CATP == 3)"
    ),
    TUMTP_OTHER = list(
      description        = "Other primary cancer type indicator (1 = primary cancer type other, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-other; lung cancer is the implicit reference)",
      notes              = "Source factor CATP level 4 = other. Decompose as TUMTP_OTHER = as.integer(CATP == 4). Used as a categorical power-model multiplier on SLPp50.",
      source_name        = "(derived from CATP == 4)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 451L,
    n_studies      = 6L,
    age_range      = "28-87 years",
    age_median     = "57 years",
    weight_range   = "38.6-129 kg",
    weight_median  = "72.3 kg",
    sex_female_pct = 50.3,
    race_ethnicity = c(Caucasian = 90.2, Black = 5.8, NativeAmerican = 0.2,
                       Asian = 0.7, Hispanic = 2.0, Other = 1.1),
    disease_state  = "Adults with cancer receiving radiation therapy. Primary tumor types: lung (50.8%), breast (17.3%), cranial glioblastoma multiforme (16.0%), other (16.0%).",
    dose_range     = "Efaproxiral 75-100 mg/kg as a 30-minute IV infusion, 2-3 times per week (Data section).",
    regions        = "Pooled six phase I-III trials; specific region distribution not reported in the poster.",
    bsa_range      = "1.31-2.50 m^2 (median 1.82)",
    albumin_range  = "2.2-5.0 g/dL (median 3.6)",
    notes          = "Demographics from Table 1 (continuous: age, weight, height, baseline hemoglobin, baseline albumin, baseline creatinine, ideal body weight, body surface area, baseline and truncated creatinine clearance; categorical: sex, primary cancer type, race). Database contained 2582 plasma EFP concentrations, 2881 RBC EFP concentrations and 2483 p50 values (the latter from phase I-II studies only); total 7946 observations."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters - Gastonguay 2005 Table 3 typical values.
    # Reference covariate values from Equation 1: BSA = 1.8 m^2, AGE = 60 yrs,
    # BALB = 3.5 g/dL, MDOS = 6800 mg.
    # ------------------------------------------------------------------
    lcl <- log(1.88);   label("Population total clearance CL (L/h)")          # Table 3: CL = 1.88 L/h (RSE 8%)
    lvc <- log(10.5);   label("Population central volume V1 (L)")             # Table 3: V1 = 10.5 L (RSE 2%)
    lq  <- log(2.58);   label("Population inter-compartmental clearance Q (L/h)") # Table 3: Q = 2.58 L/h (RSE 12%)
    lvp <- log(18.1);   label("Population peripheral volume V2 (L)")          # Table 3: V2 = 18.1 L (RSE 10%)

    # ------------------------------------------------------------------
    # RBC:plasma proportionality (algebraic, not an ODE state).
    # Crbc = slprbc * Cc; SLPRBC near 1 indicates nearly equal RBC and
    # plasma concentrations (Crbc ~ Cc).
    # ------------------------------------------------------------------
    lslprbc <- log(0.982); label("RBC:plasma proportionality SLPRBC (dimensionless)") # Table 3: SLPRBC = 0.982 (RSE 1%)

    # ------------------------------------------------------------------
    # PD: linear p50 model.  p50 = INTp50 + SLPp50 * Crbc.
    # INTp50 is the baseline (drug-free) p50; SLPp50 is the linear PD
    # slope on RBC efaproxiral concentration.
    # ------------------------------------------------------------------
    lrbase_p50 <- log(26.9);   label("Baseline (drug-free) p50 INTp50 (mmHg)")                 # Table 3: INTp50 = 26.9 mmHg (RSE ~0%)
    lslope_p50 <- log(0.0193); label("PD linear slope SLPp50 of p50 on RBC EFP (mmHg per ug/mL)") # Table 3: SLPp50 = 0.0193 mmHg/(ug/mL) (RSE 5%)

    # ------------------------------------------------------------------
    # Covariate effects on PK parameters (power-model exponents).
    # Sign convention matches Equation 1: parameter_i = TVP *
    #   prod_k (cov_k / cov_ref_k)^theta_k.
    # ------------------------------------------------------------------
    e_bsa_cl  <-  0.528;  label("Power-model exponent of (BSA/1.8) on CL")    # Table 3: theta_CL~BSA = 0.528 (RSE 101%)
    e_age_cl  <- -1.0;    label("Power-model exponent of (AGE/60) on CL")    # Table 3: theta_CL~AGE = -1 (RSE 25%)
    e_bsa_vc  <-  1.15;   label("Power-model exponent of (BSA/1.8) on V1")    # Table 3: theta_V1~BSA = 1.15 (RSE 11%)
    e_age_vc  <- -0.282;  label("Power-model exponent of (AGE/60) on V1")    # Table 3: theta_V1~AGE = -0.282 (RSE 27%)
    e_alb_vc  <-  0.357;  label("Power-model exponent of (BALB/3.5) on V1")   # Table 3: theta_V1~BALB = 0.357 (RSE 47%)
    e_bsa_q   <-  0.199;  label("Power-model exponent of (BSA/1.8) on Q")    # Table 3: theta_Q~BSA = 0.199 (RSE 426%)
    e_bsa_vp  <-  2.62;   label("Power-model exponent of (BSA/1.8) on V2")    # Table 3: theta_V2~BSA = 2.62 (RSE 29%)
    e_age_vp  <-  0.524;  label("Power-model exponent of (AGE/60) on V2")    # Table 3: theta_V2~AGE = 0.524 (RSE 81%)
    e_alb_vp  <- -2.35;   label("Power-model exponent of (BALB/3.5) on V2")   # Table 3: theta_V2~BALB = -2.35 (RSE 57%)

    # ------------------------------------------------------------------
    # Covariate effects on SLPRBC (continuous power model).
    # ------------------------------------------------------------------
    e_dose_efp_max_slprbc <- -0.125;  label("Power-model exponent of (MDOS/6800) on SLPRBC") # Table 3: theta_SLPRBC~MDOS = -0.125 (RSE 42%)
    e_age_slprbc          <- -0.109;  label("Power-model exponent of (AGE/60) on SLPRBC")   # Table 3: theta_SLPRBC~AGE = -0.109 (RSE 49%)
    e_alb_slprbc          <-  0.367;  label("Power-model exponent of (BALB/3.5) on SLPRBC") # Table 3: theta_SLPRBC~BALB = 0.367 (RSE 38%)

    # ------------------------------------------------------------------
    # Covariate effects on SLPp50 (categorical power model with lung as
    # the implicit reference).  Equation 1 form: SLPp50 = theta *
    #   (theta_CATP2)^CATP2 * (theta_CATP3)^CATP3 * (theta_CATP4)^CATP4
    # is equivalent to the log-linear form
    #   log(SLPp50) = log(theta) + CATP2 * log(theta_CATP2) + ...
    # so the canonical e_<cov>_<param> effect parameters carry the
    # natural log of the paper's linear-scale multipliers.
    # ------------------------------------------------------------------
    e_breast_slope_p50 <- log(1.11);   label("Log-multiplier on SLPp50 for breast (CATP=2 vs lung)")        # Table 3: theta_SLPp50~CATP2 = 1.11 (RSE 8%)
    e_glio_slope_p50   <- log(1.28);   label("Log-multiplier on SLPp50 for cranial GBM (CATP=3 vs lung)")   # Table 3: theta_SLPp50~CATP3 = 1.28 (RSE 18%)
    e_other_slope_p50  <- log(0.854);  label("Log-multiplier on SLPp50 for other cancer type (CATP=4 vs lung)") # Table 3: theta_SLPp50~CATP4 = 0.854 (RSE 11%)

    # ------------------------------------------------------------------
    # IIV - Table 3 reports omega^2 directly (and a back-converted CV%);
    # the three-way correlated block is on CL, V1, V2 with the reported
    # covariances.
    # Diagonal variances:
    #   omega^2(CL)     = 0.283   (CV% 53%)
    #   omega^2(V1)     = 0.0324  (CV% 18%)
    #   omega^2(V2)     = 4.66    (CV% 216%)
    #   omega^2(Q)      = 0.154   (CV% 39%)
    #   omega^2(SLPRBC) = 0.0176  (CV% 13%)
    #   omega^2(INTp50) = 0.00217 (CV% 5%)
    #   omega^2(SLPp50) = 0.0568  (CV% 24%)
    # Off-diagonals (Table 3 footnotes a, b, c):
    #   cov(CL,V1) = 0.0264, cov(CL,V2) = -0.41, cov(V1,V2) = 0.0341
    # ------------------------------------------------------------------
    etalcl + etalvc + etalvp ~ c(0.283,
                                 0.0264, 0.0324,
                                 -0.41,  0.0341, 4.66) # Table 3 omega^2 + footnotes a,b,c
    etalq          ~ 0.154                              # Table 3 omega^2(Q) = 0.154
    etalslprbc     ~ 0.0176                             # Table 3 omega^2(SLPRBC) = 0.0176
    etalrbase_p50  ~ 0.00217                            # Table 3 omega^2(INTp50) = 0.00217
    etalslope_p50  ~ 0.0568                             # Table 3 omega^2(SLPp50) = 0.0568

    # ------------------------------------------------------------------
    # Residual error - combined additive + proportional on plasma and
    # RBC; proportional-only on p50.  Table 3 reports the additive
    # components as SD on the linear scale and the proportional / p50
    # components as CV%.
    # ------------------------------------------------------------------
    addSd       <- 24.74;     label("Additive residual SD on plasma EFP (ug/mL)")      # Table 3: SD = 24.74 ug/mL (RSE 26%)
    propSd      <- 0.1296;    label("Proportional residual SD on plasma EFP (fraction)") # Table 3: CV% = 12.96% (RSE 11%)
    addSd_Crbc  <- 15.62;     label("Additive residual SD on RBC EFP (ug/mL)")         # Table 3: SD = 15.62 ug/mL (RSE 27%)
    propSd_Crbc <- 0.1556;    label("Proportional residual SD on RBC EFP (fraction)")    # Table 3: CV% = 15.56% (RSE 15%)
    propSd_p50  <- 0.0686;    label("Proportional residual SD on p50 (fraction)")        # Table 3: CV% = 6.86% (RSE 18%)
  })

  model({
    # ---- 1. Derived covariate terms ----
    # Convert canonical SI albumin (g/L) to the g/dL units used by the
    # source equation's BALB / 3.5 g/dL power-of term.
    alb_gdL <- ALB * 0.1

    # ---- 2. Individual PK parameters (power-model covariate effects) ----
    cl <- exp(lcl + etalcl) *
            (BSA / 1.8)^e_bsa_cl *
            (AGE / 60)^e_age_cl
    vc <- exp(lvc + etalvc) *
            (BSA / 1.8)^e_bsa_vc *
            (AGE / 60)^e_age_vc *
            (alb_gdL / 3.5)^e_alb_vc
    q  <- exp(lq  + etalq) *
            (BSA / 1.8)^e_bsa_q
    vp <- exp(lvp + etalvp) *
            (BSA / 1.8)^e_bsa_vp *
            (AGE / 60)^e_age_vp *
            (alb_gdL / 3.5)^e_alb_vp

    # ---- 3. Individual RBC:plasma proportionality ----
    slprbc <- exp(lslprbc + etalslprbc) *
                (DOSE_EFP_MAX_MG / 6800)^e_dose_efp_max_slprbc *
                (AGE / 60)^e_age_slprbc *
                (alb_gdL / 3.5)^e_alb_slprbc

    # ---- 4. Individual PD parameters ----
    rbase_p50 <- exp(lrbase_p50 + etalrbase_p50)
    slope_p50 <- exp(lslope_p50 + etalslope_p50 +
                       e_breast_slope_p50 * TUMTP_BREAST +
                       e_glio_slope_p50   * TUMTP_GLIO   +
                       e_other_slope_p50  * TUMTP_OTHER)

    # ---- 5. Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 6. ODE system (linear 2-compartment, IV input via central) ----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 7. Observation variables ----
    Cc   <- central / vc          # plasma EFP concentration (ug/mL)
    Crbc <- slprbc  * Cc          # RBC EFP concentration (ug/mL)
    p50  <- rbase_p50 + slope_p50 * Crbc # p50 (mmHg)

    # ---- 8. Residual error per output ----
    Cc   ~ add(addSd)      + prop(propSd)
    Crbc ~ add(addSd_Crbc) + prop(propSd_Crbc)
    p50  ~                   prop(propSd_p50)
  })
}
