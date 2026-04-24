Bajaj_2017_nivolumab <- function() {
  description <- "Two-compartment population PK model for nivolumab (anti-PD-1 IgG4) with time-varying clearance (sigmoid Emax) in patients with advanced solid tumors (Bajaj 2017)"
  reference <- "Bajaj G, Wang X, Agrawal S, Gupta M, Roy A, Feng Y. Model-based population pharmacokinetic analysis of nivolumab in patients with solid tumors. CPT Pharmacometrics Syst Pharmacol. 2017;6(1):58-66. doi:10.1002/psp4.12143"
  vignette <- "Bajaj_2017_nivolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 80 kg (Bajaj 2017 Table 1 footnote and Eqs. 7 and 10).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Baseline CKD-EPI estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 90 mL/min/1.73 m^2 (Bajaj 2017 Table 1 footnote and Eq. 7). Source column name in Bajaj 2017 is eGFR; stored under the canonical CRCL. Bajaj 2017 Methods states eGFR was estimated using the CKD-EPI equation.",
      source_name        = "eGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The paper's own reference category is female (see notes).",
      notes              = "Bajaj 2017 encodes sex as a male-indicator (1 = male, 0 = female) with female as the reference category (Table 1 footnote: 'white female reference'). To store under the canonical SEXF (1 = female, 0 = male) while preserving Bajaj's female-reference CL_REF and VC_REF, the effect is applied in model() as exp(e_sex_cl * (1 - SEXF)) and exp(e_sex_vc * (1 - SEXF)), so SEXF = 1 yields factor 1 and SEXF = 0 yields the paper's male-vs-female exp-coefficient.",
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Indicator for Asian race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; pooled White, Black / African American, and Other)",
      notes              = "Exponential effect on CL for Asian race (Bajaj 2017 Table 1 row CL_RAAS and Eq. 7). Renamed from the source column RAAS to the canonical RACE_ASIAN per covariate-columns.md.",
      source_name        = "RAAS"
    ),
    ECOG_GE1 = list(
      description        = "Baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator (1 if ECOG >= 1, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status = 0, i.e., fully active)",
      notes              = "Exponential effect on CL for patients with ECOG >= 1 (Bajaj 2017 Table 1 row CL_BPS and Eq. 7). Bajaj 2017 ECOG values came directly from each study except CA209025, which collected Karnofsky Performance Status (KPS) and was mapped to ECOG via the Oken 1982 crosswalk before binarization. Renamed from the source column PS to the canonical ECOG_GE1 per covariate-columns.md.",
      source_name        = "PS"
    )
  )

  population <- list(
    n_subjects     = 1895L,
    n_studies      = 11L,
    age_range      = "mean 61.1 years (SD 11.1)",
    age_median     = "not reported (mean reported instead)",
    weight_range   = "mean 79.1 kg (SD 19.3); model-application range 34.1 - 168.2 kg",
    weight_median  = "not reported (mean reported instead)",
    sex_female_pct = 33.3,
    race_ethnicity = c(White = 88.92, Asian = 6.44, `Black/African American` = 2.8, Other = 1.74),
    disease_state  = "Advanced / metastatic solid tumors (melanoma 29.82%, NSCLC 34.78%, RCC 31.93%, other 3.48%)",
    dose_range     = "0.3 - 10.0 mg/kg IV infusion (1-hour) Q2W or Q3W across 11 trials",
    regions        = "Global (US, EU, Japan) across phase I / II / III studies",
    ecog_distribution = "ECOG 0 38.73%, ECOG 1 58.52%, ECOG 2 2.74%",
    renal_function = "Baseline CKD-EPI eGFR mean 78.5 (SD 21.6) mL/min/1.73 m^2",
    notes          = "Baseline demographics per Bajaj 2017 Table 3 (N = 1,895). Studies pooled: MDX1106-01, ONO-4538-01, MDX1106-03, CA209010, CA209063, ONO-4538-02, CA209017, CA209037, CA209025, CA209057, CA209066 (Bajaj 2017 Table 2). Baseline lactate dehydrogenase mean 350.9 (SD 397.6) U/L. Liver-dysfunction groups: normal 89.08%, mild 9.82%, moderate 0.11%, severe 0.05%, missing 0.95%."
  )

  ini({
    # Structural parameters - reference values for a white female patient
    # weighing 80 kg, with eGFR 90 mL/min/1.73 m^2, ECOG performance
    # status = 0, non-Asian race (Bajaj 2017 Table 1 footnote a).
    # CL and Q are reported in mL/h in Table 1; converted to L/day
    # (x 24 / 1000) because this model keeps time in days.
    lcl  <- log(9.4  * 24 / 1000); label("Baseline clearance CL_BASE,REF for the paper's reference covariates (L/day)") # Bajaj 2017 Table 1: CL_REF = 9.4 mL/h
    lvc  <- log(3.63);             label("Central volume of distribution VC_REF (L)")                                   # Bajaj 2017 Table 1: VC_REF = 3.63 L
    lq   <- log(32.1 * 24 / 1000); label("Intercompartmental clearance Q_REF (L/day)")                                  # Bajaj 2017 Table 1: Q_REF = 32.1 mL/h
    lvp  <- log(2.78);             label("Peripheral volume of distribution VP_REF (L)")                                # Bajaj 2017 Table 1: VP_REF = 2.78 L

    # Covariate effects on CL (Bajaj 2017 Table 1 and Eq. 7). Power on BW and
    # eGFR; exponential on ECOG_GE1, sex (applied to male-indicator), and RACE_ASIAN.
    e_wt_cl         <-  0.566;  label("Power exponent of WT on CL (unitless)")                                          # Bajaj 2017 Table 1: CL_BW   = 0.566
    e_crcl_cl       <-  0.186;  label("Power exponent of CRCL (eGFR) on CL (unitless)")                                 # Bajaj 2017 Table 1: CL_eGFR = 0.186
    e_ecog_ge1_cl   <-  0.172;  label("Exponential coefficient of ECOG_GE1 on CL (unitless)")                           # Bajaj 2017 Table 1: CL_BPS  = 0.172
    e_sex_cl        <-  0.165;  label("Exponential coefficient of male sex on CL (unitless; applied as (1 - SEXF))")    # Bajaj 2017 Table 1: CL_SEX  = 0.165
    e_race_asian_cl <- -0.125;  label("Exponential coefficient of RACE_ASIAN on CL (unitless)")                         # Bajaj 2017 Table 1: CL_RAAS = -0.125

    # Time-varying clearance (sigmoid Emax of time since first dose; Bajaj 2017
    # Table 1 and Eq. 8). Emax is the maximal fractional change in CL; at t >> T50,
    # CL approaches CL_base * exp(Emax). With Emax = -0.295 the steady-state CL is
    # 0.745 of baseline (a 25.5% reduction), matching the paper's reported mean
    # maximal reduction of ~24.5%. T50 is reported in hours in Table 1 (1.41e3 h);
    # converted to days (/24) because this model keeps time in days.
    cl_emax <- -0.295;           label("Maximal fractional change in CL_EMAX (unitless)")                               # Bajaj 2017 Table 1: CL_EMAX = -0.295
    t50     <-  1.41e3 / 24;     label("Time at which the change in CL is 50%% of CL_EMAX (days)")                      # Bajaj 2017 Table 1: CL_T50  = 1.41e3 h
    cl_hill <-  3.15;            label("Hill / sigmoidicity exponent of time on CL (unitless)")                         # Bajaj 2017 Table 1: CL_HILL = 3.15

    # Covariate effects on VC (Bajaj 2017 Table 1 and Eq. 10). Power on BW;
    # exponential on sex (applied to male-indicator).
    e_wt_vc  <- 0.597; label("Power exponent of WT on VC (unitless)")                                                   # Bajaj 2017 Table 1: VC_BW  = 0.597
    e_sex_vc <- 0.152; label("Exponential coefficient of male sex on VC (unitless; applied as (1 - SEXF))")             # Bajaj 2017 Table 1: VC_SEX = 0.152

    # Inter-individual variability (Bajaj 2017 Table 1). CL and VC are
    # log-normal with a single CL:VC covariance (correlation 0.352); VP is an
    # independent log-normal eta; Emax has an INDEPENDENT ADDITIVE eta
    # (Emax_i = cl_emax + etacl_emax per Bajaj 2017 Eq. 3).
    etalcl + etalvc ~ c(0.123,
                        0.0432, 0.123)                                                                                  # Bajaj 2017 Table 1: omega^2_CL = 0.123, cov_CL:VC = 0.0432, omega^2_VC = 0.123
    etalvp     ~ 0.258                                                                                                  # Bajaj 2017 Table 1: omega^2_VP = 0.258
    etacl_emax ~ 0.0719                                                                                                 # Bajaj 2017 Table 1: omega^2_EMAX = 0.0719

    # Residual error (proportional only; Bajaj 2017 Table 1).
    propSd <- 0.215; label("Proportional residual error (fraction)")                                                    # Bajaj 2017 Table 1: proportional error = 0.215
  })
  model({
    # Derived sex term: Bajaj 2017 encodes sex as a male-indicator with female as
    # the reference category, so (1 - SEXF) reproduces the paper's male = 1 column
    # while keeping SEXF (1 = female) as the canonical storage convention.
    sex_male <- 1 - SEXF

    # Individual baseline CL and VC with covariate adjustments
    # (Bajaj 2017 Eqs. 7 and 10; CL_REF and VC_REF are the paper's typical values
    # at WT = 80 kg, CRCL = 90 mL/min/1.73 m^2, ECOG_GE1 = 0, female, non-Asian).
    cl_base <- exp(lcl + etalcl) *
      (WT   / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_ecog_ge1_cl   * ECOG_GE1) *
      exp(e_sex_cl        * sex_male) *
      exp(e_race_asian_cl * RACE_ASIAN)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sex_vc * sex_male)

    # Time-varying clearance (Bajaj 2017 Eqs. 8 and 3). Additive IIV on Emax.
    emax_i <- cl_emax + etacl_emax
    cl     <- cl_base * exp(emax_i * t^cl_hill / (t50^cl_hill + t^cl_hill))

    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
