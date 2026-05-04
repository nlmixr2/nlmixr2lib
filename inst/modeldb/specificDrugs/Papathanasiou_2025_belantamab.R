Papathanasiou_2025_belantamab <- function() {
  description <- "Two-compartment population PK model for the antibody-drug conjugate (ADC) belantamab mafodotin in patients with relapsed/refractory multiple myeloma, with sigmoidal time-varying clearance and covariate effects of baseline body weight, BMI, albumin, soluble BCMA, serum IgG, race, and combination therapy (Papathanasiou 2025; ADC moiety only -- the cys-mcMMAF payload sub-model is not included; see vignette for rationale)"
  reference <- "Papathanasiou T, Strougo A, Roy A, Vakkalagadda B, Stein A, Jewell RC, Boer J, Dahmane E. Population pharmacokinetics for belantamab mafodotin monotherapy and combination therapies in patients with relapsed/refractory multiple myeloma. Clin Pharmacokinet. 2025;64(6):925-942. doi:10.1007/s40262-025-01508-1"
  vignette <- "Papathanasiou_2025_belantamab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponent 0.929 on Vc and Vp (shared via theta_V_WTBL) and 0.542 on CL and Q (shared via theta_CL_WTBL); reference 75 kg per Papathanasiou 2025 typical-patient definition (Methods, Sect. 2.3). Source column WTBL maps to canonical WT.",
      source_name        = "WTBL"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponent -0.459 on Vc; reference 27 kg/m^2 per Papathanasiou 2025 typical-patient definition. Source column IBMIBL maps to canonical BMI.",
      source_name        = "IBMIBL"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline (SI units, g/L). Power exponents -0.698 on CL, -0.302 on Vc, +0.567 on Vp; reference 40 g/L per Papathanasiou 2025 typical-patient definition. Source column ALBBL maps to canonical ALB.",
      source_name        = "ALBBL"
    ),
    IGG = list(
      description        = "Baseline serum immunoglobulin G concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponents +0.170 on CL and +0.192 on Imax; reference 15 g/L per Papathanasiou 2025 typical-patient definition. Source column IGGBL maps to canonical IGG. The published TI50 equation in Table 2 also references an IgG-on-TI50 exponent (theta_TI50_IGGBL), but the parameter table provides no value for it; that term is therefore omitted (see vignette Assumptions and deviations).",
      source_name        = "IGGBL"
    ),
    SBCMA = list(
      description        = "Baseline serum soluble B-cell maturation antigen (sBCMA) concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. 1 ng/mL is numerically equivalent to 1 ug/L (the unit reported in Papathanasiou 2025 Table 1). Power exponents +0.113 on CL, +0.0401 on Vc, +0.160 on Imax; reference 50 ng/mL per Papathanasiou 2025 typical-patient definition. Source column SBCMABL maps to canonical SBCMA.",
      source_name        = "SBCMABL"
    ),
    RACE_ASIAN = list(
      description        = "Indicator for Asian race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; pooled White, Black/African American, and Other)",
      notes              = "Multiplicative factor 0.913 on CL when RACE_ASIAN = 1 (Papathanasiou 2025 Table 2). Reference is non-Asian race; mutually exclusive with RACE_BLACK in the source dataset.",
      source_name        = "RACE"
    ),
    RACE_BLACK = list(
      description        = "Indicator for Black/African American race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; pooled White, Asian, and Other)",
      notes              = "Multiplicative factor 0.861 on CL when RACE_BLACK = 1 (Papathanasiou 2025 Table 2). Reference is non-Black/African American race; mutually exclusive with RACE_ASIAN in the source dataset.",
      source_name        = "RACE"
    ),
    COMBO_BELAMAF = list(
      description        = "Belantamab mafodotin combination therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (belantamab mafodotin monotherapy)",
      notes              = "Multiplicative factor 1.44 on the typical Imax when COMBO_BELAMAF = 1 (Papathanasiou 2025 Table 2: theta_IMAX_COMBO). Pools the Bor-Dex / Len-Dex combination backbones tested in DREAMM-6 and DREAMM-7; the Pom-Dex backbone (DREAMM-8) was used only for external validation and would also be encoded as COMBO_BELAMAF = 1 if simulating that regimen.",
      source_name        = "COMBO"
    )
  )

  population <- list(
    n_subjects     = 977L,
    n_studies      = 6L,
    age_range      = "32 - 89 years",
    age_median     = "66 years",
    weight_range   = "37 - 170 kg",
    weight_median  = "74.0 kg",
    sex_female_pct = 43.6,
    race_ethnicity = c(White = 78.1, `Black/African American` = 6.2, Asian = 13.6, Other = 0.8, Missing = 1.2),
    disease_state  = "Relapsed/refractory multiple myeloma (RRMM); patients had received >= 1 prior line of therapy depending on the study setting (DREAMM-2 4L+, DREAMM-3 3L+, DREAMM-6/-7 2L+).",
    dose_range     = "Belantamab mafodotin IV every 21 days; pivotal labelled regimen 2.5 mg/kg q3w (cycle-1 reference dose used for exposure simulations); studies also included 1.9, 2.5, and 3.4 mg/kg q3w cohorts.",
    regions        = "Global: Europe 42.3%, Northeast Asia 9.3%, North America 20.8%, Rest of world 27.6% (Table 1).",
    treatment_mix  = "Monotherapy 59.7%, bortezomib + dexamethasone 35.7%, lenalidomide + dexamethasone 4.6% (DREAMM-8 pomalidomide + dexamethasone cohort, n = 150, was held out for external validation only).",
    biomarker_summary = "Median (range): sBCMA 56.0 (2.08 - 2030) ng/mL; serum IgG 13.1 (0.350 - 119) g/L; albumin 39.0 (19.0 - 57.0) g/L; BMI 26.7 (14.0 - 48.4) kg/m^2; beta-2 microglobulin 297 (94.9 - 5190) nmol/L (Table 1).",
    renal_distribution = "Normal (eGFR >= 90 mL/min) 31.9%, mild 42.0%, moderate 23.0%, severe 2.8%, end-stage 0.3% (NCI-ODWG / MDRD; Table 1).",
    hepatic_distribution = "Normal 84.5%, mild 11.9%, moderate 0.5%, severe 0.1% (NCI-ODWG; Table 1).",
    external_validation = "DREAMM-8 (NCT04484623; n = 150 patients, 1221 ADC samples) belantamab mafodotin + pomalidomide + dexamethasone; not used for parameter estimation.",
    reference_subject = "65-year-old male, 75 kg WT, 27 kg/m^2 BMI, sBCMA 50 ng/mL, IgG 15 g/L, albumin 40 g/L, monotherapy, non-Asian non-Black race (Methods, Sect. 2.3).",
    notes          = "Pooled analysis of DREAMM-2 (NCT03525678; n = 218), DREAMM-3 (NCT04162210; n = 217), DREAMM-6 (NCT03544281; n = 152), DREAMM-7 (NCT04246047; n = 242), DREAMM-12 (NCT04398745; n = 23), and DREAMM-14 (NCT05064358; n = 125). 8880 ADC concentrations entered estimation. Sex was not retained as a covariate because of collinearity with body weight."
  )

  ini({
    # Structural parameters (typical values for the reference subject:
    # 65-year-old male, 75 kg, BMI 27 kg/m^2, sBCMA 50 ng/mL, IgG 15 g/L,
    # albumin 40 g/L, monotherapy, non-Asian / non-Black race;
    # Papathanasiou 2025 Methods Sect. 2.3 and Table 2).
    lcl   <- log(0.926); label("Initial typical clearance CL at reference covariates (L/day)")             # Papathanasiou 2025 Table 2: CL = 0.926 L/day
    lvc   <- log(4.21);  label("Central volume of distribution Vc at reference covariates (L)")           # Papathanasiou 2025 Table 2: ADC Vc = 4.21 L
    lq    <- log(0.711); label("Intercompartmental clearance Q at reference covariates (L/day)")          # Papathanasiou 2025 Table 2: Q = 0.711 L/day
    lvp   <- log(6.63);  label("Peripheral volume of distribution Vp at reference covariates (L)")        # Papathanasiou 2025 Table 2: ADC Vp = 6.63 L

    # Sigmoidal time-varying CL (Papathanasiou 2025 Table 2 PK parameter
    # estimation block: CL_Time = Imax_i * Time^Gamma / (TI50_i^Gamma +
    # Time^Gamma) and CL_i = theta_CL * exp(CL_Time) * ... covariates ...).
    # At t = 0 the multiplier is 1; as t -> infinity the multiplier approaches
    # exp(Imax_i). For monotherapy at reference covariates the typical Imax
    # is -0.403, so steady-state CL = 0.926 * exp(-0.403) = 0.619 L/day
    # (matches the paper's reported 33.2% reduction); for combination
    # therapy the Imax multiplier is multiplied by 1.44, giving an
    # effective Imax = -0.580 and steady-state CL = 0.519 L/day (matches
    # the reported 44.0% reduction).
    imax  <- -0.403;     label("Maximal log-fold change in CL (typical, monotherapy reference) (unitless)") # Papathanasiou 2025 Table 2: Imax
    lti50 <- log(66.4);  label("Log of time at which time-varying CL change is half-maximal (log days)")    # Papathanasiou 2025 Table 2: TI50 = 66.4 days
    gamma <- 2.87;       label("Hill / sigmoidicity exponent of time on CL (unitless)")                    # Papathanasiou 2025 Table 2: Gamma

    # Covariate effect parameters (Papathanasiou 2025 Table 2). The two
    # weight exponents (theta_V_WTBL on volumes, theta_CL_WTBL on
    # clearances) are each shared between two structural parameters as
    # written in the source paper.
    e_wt_vc_vp    <-  0.929;  label("Shared power exponent of WT on Vc and Vp (theta_V_WTBL)")             # Papathanasiou 2025 Table 2: theta_V_WTBL
    e_wt_cl_q     <-  0.542;  label("Shared power exponent of WT on CL and Q (theta_CL_WTBL)")             # Papathanasiou 2025 Table 2: theta_CL_WTBL
    e_alb_cl      <- -0.698;  label("Power exponent of ALB on CL")                                         # Papathanasiou 2025 Table 2: theta_CL_ALBBL
    e_alb_vc      <- -0.302;  label("Power exponent of ALB on Vc")                                         # Papathanasiou 2025 Table 2: theta_ADC_Vc_ALBBL
    e_alb_vp      <-  0.567;  label("Power exponent of ALB on Vp")                                         # Papathanasiou 2025 Table 2: theta_ADC_Vp_ALBBL
    e_sbcma_cl    <-  0.113;  label("Power exponent of SBCMA on CL")                                       # Papathanasiou 2025 Table 2: theta_CL_SBCMABL
    e_sbcma_vc    <-  0.0401; label("Power exponent of SBCMA on Vc")                                       # Papathanasiou 2025 Table 2: theta_ADC_Vc_SBCMABL
    e_igg_cl      <-  0.170;  label("Power exponent of IGG on CL")                                         # Papathanasiou 2025 Table 2: theta_CL_IGGBL
    e_bmi_vc      <- -0.459;  label("Power exponent of BMI on Vc")                                         # Papathanasiou 2025 Table 2: theta_ADC_Vc_IBMIBL
    e_race_asian_cl <- 0.913; label("Multiplicative factor of RACE_ASIAN on CL")                           # Papathanasiou 2025 Table 2: theta_CL_RACEA
    e_race_black_cl <- 0.861; label("Multiplicative factor of RACE_BLACK on CL")                           # Papathanasiou 2025 Table 2: theta_CL_RACEB
    e_combo_imax  <-  1.44;   label("Multiplicative factor of COMBO_BELAMAF on Imax")                      # Papathanasiou 2025 Table 2: theta_IMAX_COMBO
    e_igg_imax    <-  0.192;  label("Power exponent of IGG on Imax")                                       # Papathanasiou 2025 Table 2: theta_IMAX_IGGBL
    e_sbcma_imax  <-  0.160;  label("Power exponent of SBCMA on Imax")                                     # Papathanasiou 2025 Table 2: theta_IMAX_SBCMABL

    # Inter-individual variability (Papathanasiou 2025 Table 2). All etas
    # are log-normal (omega^2 = log(CV^2 + 1)) except etaimax, which is
    # additive on the linear Imax scale (the paper's footnote: normal
    # distribution CV% = SQRT(Omega) / theta * 100, so
    # Omega = (CV/100 * |theta|)^2). CL and Vc are correlated as a 2x2
    # block; Q, Vp, Imax, and TI50 are independent etas.
    etalcl + etalvc ~ c(0.06593,
                        0.0328, 0.03922)                  # Papathanasiou 2025 Table 2: CL CV 26.1%, Vc CV 20.0%, CL-Vc covariance 0.0328 (correlation 0.646)
    etalq    ~ 0.03546                                    # Papathanasiou 2025 Table 2: Q CV 19.0%
    etalvp   ~ 0.08945                                    # Papathanasiou 2025 Table 2: Vp CV 30.6%
    etaimax  ~ 0.01366                                    # Papathanasiou 2025 Table 2: Imax CV 29.0% (normal-distribution form, theta = -0.403)
    etalti50 ~ 0.39236                                    # Papathanasiou 2025 Table 2: TI50 CV 69.3%

    # Residual error (Papathanasiou 2025 Table 2): "Y = ln(IPRED) + eps"
    # with Var(eps) = 0.0633 (log(ng/mL))^2. Additive on the natural-log
    # scale is equivalent to proportional in linear space for small
    # errors, with proportional SD = sqrt(0.0633) ~ 0.2516.
    propSd <- 0.2516; label("Proportional residual error (SD; equivalent to additive-on-log-scale SD)")    # Papathanasiou 2025 Table 2: sqrt(0.0633) = 0.2516
  })

  model({
    # 1. Individual Imax (additive eta on the linear scale; covariate
    # multipliers per Papathanasiou 2025 Table 2 Imax equation). The
    # COMBO factor is applied as e_combo_imax^COMBO_BELAMAF so that
    # COMBO_BELAMAF = 0 yields factor 1 and COMBO_BELAMAF = 1 yields
    # factor 1.44.
    imax_typ <- imax *
      e_combo_imax^COMBO_BELAMAF *
      (IGG   / 15)^e_igg_imax *
      (SBCMA / 50)^e_sbcma_imax
    imax_i <- imax_typ + etaimax

    # 2. Individual TI50 (log-normal eta). The published Table 2
    # equation also references an IgG-on-TI50 exponent
    # (theta_TI50_IGGBL), but no numeric value is reported anywhere in
    # the paper for that parameter; the term is therefore omitted (see
    # vignette Assumptions and deviations).
    ti50_i <- exp(lti50 + etalti50)

    # 3. Individual baseline CL (covariate-adjusted typical CL at t = 0,
    # before applying the time-varying multiplier; Papathanasiou 2025
    # Table 2 CL_i equation). Race factors are encoded as
    # multiplier^indicator so that subjects with RACE_ASIAN = 0 and
    # RACE_BLACK = 0 contribute factor 1.
    cl_base <- exp(lcl + etalcl) *
      (WT    / 75)^e_wt_cl_q *
      (ALB   / 40)^e_alb_cl *
      (SBCMA / 50)^e_sbcma_cl *
      (IGG   / 15)^e_igg_cl *
      e_race_asian_cl^RACE_ASIAN *
      e_race_black_cl^RACE_BLACK

    # 4. Sigmoidal time-varying CL multiplier (Papathanasiou 2025 Table 2).
    # rxode2's t variable is time since simulation start; with the first
    # dose at t = 0 this matches the paper's "Time" since first dose.
    cl <- cl_base * exp(imax_i * t^gamma / (ti50_i^gamma + t^gamma))

    # 5. Other individual structural parameters (Papathanasiou 2025
    # Table 2 Vc, Vp, Q equations).
    vc <- exp(lvc + etalvc) *
      (WT    / 75)^e_wt_vc_vp *
      (ALB   / 40)^e_alb_vc *
      (BMI   / 27)^e_bmi_vc *
      (SBCMA / 50)^e_sbcma_vc

    vp <- exp(lvp + etalvp) *
      (WT  / 75)^e_wt_vc_vp *
      (ALB / 40)^e_alb_vp

    q  <- exp(lq + etalq) *
      (WT / 75)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 6. ODE system: linear two-compartment ADC kinetics with
    # time-varying CL. Dose in mg, volumes in L -> central / vc has
    # units mg/L = ug/mL.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
