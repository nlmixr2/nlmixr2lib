Bajaj_2017_nivolumab_ddmore <- function() {
  description <- "Two-compartment population PK model for nivolumab (anti-PD-1 IgG4) with time-varying clearance (sigmoid Emax of time since first dose) in patients with advanced solid tumors (Bajaj 2017; DDMORE Foundation Model Repository entry DDMODEL00000284). DDMORE-source replicate of inst/modeldb/specificDrugs/Bajaj_2017_nivolumab.R; parameter values are taken from the DDMORE bundle's Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE block and time is kept in hours to mirror the bundle directly."
  reference <- paste(
    "Bajaj G, Wang X, Agrawal S, Gupta M, Roy A, Feng Y. (2017).",
    "Model-based population pharmacokinetic analysis of nivolumab in patients with solid tumors.",
    "CPT Pharmacometrics Syst Pharmacol 6(1):58-66.",
    "doi:10.1002/psp4.12143.",
    "DDMORE Foundation Model Repository: DDMODEL00000284.",
    sep = " "
  )
  vignette <- "Bajaj_2017_nivolumab_ddmore"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")
  ddmore_id    <- "DDMODEL00000284"
  replicate_of <- "inst/modeldb/specificDrugs/Bajaj_2017_nivolumab.R"

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 80 kg (Output_real_Nivo-PPK.lst $PK lines 100, 159: BBWT_R = 80 kg). DDMORE bundle source column is BBWT (Simulated_pkdata1_dataset.csv); stored under canonical WT.",
      source_name        = "BBWT"
    ),
    CRCL = list(
      description        = "Baseline CKD-EPI estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 90 mL/min/1.73 m^2 (Output_real_Nivo-PPK.lst $PK line 113: BGFR_R = 90). Bajaj 2017 Methods states eGFR was estimated using the CKD-EPI equation. DDMORE bundle source column is BGFR; stored under canonical CRCL.",
      source_name        = "BGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The bundle's own reference category is female (SEXN = 2), see notes.",
      notes              = "Bajaj 2017 (and the DDMORE bundle) encode sex as SEXN with 1 = male and 2 = female and use female as the reference category in CL and Vc (Output_real_Nivo-PPK.lst $PK lines 132-133, 168, 179: 'reference is SEX=2 (Female)'). To store under the canonical SEXF (1 = female, 0 = male) while preserving the bundle's female-reference CL_REF and VC_REF, the effect is applied in model() as exp(e_sex_cl * (1 - SEXF)) and exp(e_sex_vc * (1 - SEXF)), so SEXF = 1 yields factor 1 and SEXF = 0 yields the bundle's male-vs-female exp-coefficient. Source column SEXN; convert as SEXF = as.integer(SEXN == 2).",
      source_name        = "SEXN"
    ),
    RACE_ASIAN = list(
      description        = "Indicator for Asian race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; pooled White, Black / African American, and Other)",
      notes              = "Exponential effect on CL for Asian race (Output_real_Nivo-PPK.lst $THETA line 51 and $PK lines 151, 174: theta18 = CL_RAAS, applied as exp(theta18) when RACE_I == 3). DDMORE bundle source column is RACEN with codes 1 = White, 2 = African American, 3 = Asian, 4 = Other (.mod header line 13); canonical RACE_ASIAN = as.integer(RACEN == 3).",
      source_name        = "RACEN"
    ),
    ECOG_GE1 = list(
      description        = "Baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator (1 if ECOG >= 1, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status = 0, i.e., fully active)",
      notes              = "Exponential effect on CL for patients with ECOG >= 1 (Output_real_Nivo-PPK.lst $THETA line 46 and $PK lines 138, 169: theta13 = CL_PS_1, applied as exp(theta13) when PS >= 1). DDMORE bundle source column is PS with codes 0 / 1 already binarized (.mod header line 14); canonical ECOG_GE1 takes the same numeric values. Bajaj 2017 ECOG values came directly from each study except CA209025, which collected Karnofsky Performance Status (KPS) and was mapped to ECOG via the Oken 1982 crosswalk before binarization.",
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
    notes          = "Population summary mirrors the Bajaj 2017 publication (Table 3, N = 1,895; 11 pooled trials: MDX1106-01, ONO-4538-01, MDX1106-03, CA209010, CA209063, ONO-4538-02, CA209017, CA209037, CA209025, CA209057, CA209066). The DDMORE bundle's NONMEM run (Output_real_Nivo-PPK.lst lines 258-259) reports TOT. NO. OF OBS RECS = 12,292 across TOT. NO. OF INDIVIDUALS = 1,895, matching the publication's analysis dataset. The full FCM (full covariate model) with non-significant terms FIXED to 0 is what the bundle's .mod and .lst encode (24 thetas; 9 fixed at 0 for AGE, BLDH, BALB, melanoma, others-tumor, RCC-tumor, African-American race, hepatic dysfunction, plus TH5 = additive residual error fixed at 0); the model file below carries only the significant covariates because thetas fixed to zero have no effect on predictions and would clutter the source-trace."
  )

  ini({
    # Structural parameters - reference values for a white female patient
    # weighing 80 kg, with eGFR 90 mL/min/1.73 m^2, ECOG performance
    # status = 0, non-Asian race (Output_real_Nivo-PPK.lst $PK lines 100, 113,
    # 132-133, 137, 148-149: BBWT_R = 80, BGFR_R = 90, SEX_I default = 1
    # but reference category is SEX = 2, PS reference 0, RACE reference 1=White).
    # Time is kept in HOURS to mirror the .lst directly; values reproduce the
    # FINAL PARAMETER ESTIMATE block (lst lines 524-525):
    #   TH1 = 9.40E-03 (CL, L/h), TH2 = 3.63E+00 (V1, L), TH3 = 3.21E-02 (Q, L/h),
    #   TH4 = 2.78E+00 (V2, L).
    lcl  <- log(0.00940); label("Baseline clearance CL_REF for the bundle's reference covariates (L/h)") # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH1 = 9.40E-03 L/h
    lvc  <- log(3.63);    label("Central volume of distribution V1_REF (L)")                            # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH2 = 3.63 L
    lq   <- log(0.0321);  label("Intercompartmental clearance Q_REF (L/h)")                             # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH3 = 3.21E-02 L/h
    lvp  <- log(2.78);    label("Peripheral volume of distribution V2_REF (L)")                         # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH4 = 2.78 L

    # Covariate effects on CL (Output_real_Nivo-PPK.lst $PK lines 104, 116, 133, 138, 151).
    # Power on WT and CRCL; exponential on ECOG_GE1, sex (applied to male-indicator
    # to preserve the female-reference convention), and RACE_ASIAN. Non-significant
    # FCM thetas (TH8 = CL_AGE, TH10 = CL_BLDH, TH11 = CL_BALB, TH14 = CL_MEL,
    # TH15 = CL_OTH, TH16 = CL_RCC, TH17 = CL_RAAA, TH19 = CL_HEPA) are FIXED
    # at 0 in the bundle and omitted here.
    e_wt_cl         <-  0.566;  label("Power exponent of WT on CL (unitless)")                                       # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH7  = 5.66E-01 (CL_BBWT)
    e_crcl_cl       <-  0.186;  label("Power exponent of CRCL (eGFR) on CL (unitless)")                              # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH9  = 1.86E-01 (CL_GFR)
    e_sex_cl        <-  0.165;  label("Exponential coefficient of male sex on CL (unitless; applied as (1 - SEXF))") # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH12 = 1.65E-01 (CL_SEX)
    e_ecog_ge1_cl   <-  0.172;  label("Exponential coefficient of ECOG_GE1 on CL (unitless)")                        # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH13 = 1.72E-01 (CL_PS)
    e_race_asian_cl <- -0.125;  label("Exponential coefficient of RACE_ASIAN on CL (unitless)")                      # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH18 = -1.25E-01 (CL_RAAS)

    # Time-varying clearance (sigmoid Emax of time since first dose; Output_real_Nivo-PPK.lst
    # $PK lines 162-164 and 79-81):
    #   EMAX_i = AEMAX + ZEMAX(eta)
    #   CL_TIME = exp(EMAX_i * TIME^HILL / (T50^HILL + TIME^HILL))
    # AEMAX is the maximal log-fractional change in CL; at t >> T50, CL approaches
    # CL_base * exp(AEMAX). With AEMAX = -0.295 the steady-state CL is
    # exp(-0.295) = 0.745 of baseline (a ~25% reduction). T50 is in hours.
    cl_emax <- -0.295; label("Maximal fractional change in CL_EMAX (unitless)")           # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH22 = -2.95E-01 (CL_EMAX)
    t50     <-  1410;  label("Time at which the change in CL is 50%% of CL_EMAX (h)")     # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH23 =  1.41E+03 h (CL_T50)
    cl_hill <-  3.15;  label("Hill / sigmoidicity exponent of time on CL (unitless)")     # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH24 =  3.15E+00 (CL_HILL)

    # Covariate effects on VC (Output_real_Nivo-PPK.lst $PK lines 159-160, 178-179).
    # Power on WT; exponential on sex (applied to male-indicator).
    e_wt_vc  <- 0.597; label("Power exponent of WT on VC (unitless)")                                       # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH20 = 5.97E-01 (VC_BBWT)
    e_sex_vc <- 0.152; label("Exponential coefficient of male sex on VC (unitless; applied as (1 - SEXF))") # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH21 = 1.52E-01 (VC_SEX)

    # Inter-individual variability (Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE
    # OMEGA block lines 532-544). CL and V1 are correlated log-normal etas (BLOCK(2)
    # in $OMEGA: var_CL = 0.123, cov_CL_V1 = 0.0432, var_V1 = 0.123 -> correlation 0.352).
    # V2 (peripheral) and EMAX have independent diagonal etas. Note ETA on EMAX is
    # ADDITIVE: EMAX_i = cl_emax + etacl_emax (.lst $PK line 163 / source line 79).
    etalcl + etalvc ~ c(0.123,
                        0.0432, 0.123)                                                                # Output_real_Nivo-PPK.lst OMEGA BLOCK(2): var(ETA1)=1.23E-01, cov(ETA1,ETA2)=4.32E-02, var(ETA2)=1.23E-01
    etalvp     ~ 0.258                                                                                # Output_real_Nivo-PPK.lst OMEGA: var(ETA3)=2.58E-01
    etacl_emax ~ 0.0719                                                                               # Output_real_Nivo-PPK.lst OMEGA: var(ETA4)=7.19E-02

    # Residual error (proportional only; the bundle fixes the additive residual
    # error TH5 to 0 (.lst $THETA line 38), and SIGMA is fixed at 1 so the
    # proportional SD is theta(6) directly). FINAL PARAMETER ESTIMATE TH6 = 0.215.
    propSd <- 0.215; label("Proportional residual error (fraction)")                                  # Output_real_Nivo-PPK.lst FINAL PARAMETER ESTIMATE TH6 = 2.15E-01 (PERR)
  })
  model({
    # Derived sex term: the bundle encodes sex as a male-indicator with female as
    # the reference category (.lst $PK line 168: 'reference is SEX=2 (Female)'),
    # so (1 - SEXF) reproduces the bundle's male = 1 column while keeping the
    # canonical SEXF (1 = female) storage convention.
    sex_male <- 1 - SEXF

    # Individual baseline CL and VC with covariate adjustments (Output_real_Nivo-PPK.lst
    # $PK lines 167-179). CL_REF and VC_REF are the bundle's typical values at
    # WT = 80 kg, CRCL = 90 mL/min/1.73 m^2, ECOG_GE1 = 0, female, non-Asian.
    cl_base <- exp(lcl + etalcl) *
      (WT   / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_sex_cl        * sex_male) *
      exp(e_ecog_ge1_cl   * ECOG_GE1) *
      exp(e_race_asian_cl * RACE_ASIAN)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sex_vc * sex_male)

    # Time-varying clearance (.lst $PK lines 162-164). Additive eta on Emax.
    emax_i <- cl_emax + etacl_emax
    cl     <- cl_base * exp(emax_i * t^cl_hill / (t50^cl_hill + t^cl_hill))

    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L -> central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
