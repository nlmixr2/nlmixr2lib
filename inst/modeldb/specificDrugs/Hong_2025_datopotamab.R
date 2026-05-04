Hong_2025_datopotamab <- function() {
  description <- "Coupled population PK model for datopotamab deruxtecan (Dato-DXd, anti-TROP2 antibody-drug conjugate) and its released payload DXd in adults with advanced solid tumors (Hong 2025). Dato-DXd disposition is a two-compartment model with parallel linear (CL_lin) and Michaelis-Menten (Vmax / Km) elimination from the central compartment. DXd is a one-compartment model whose formation rate equals the total Dato-DXd elimination rate (linear + nonlinear) scaled by the molecular-weight ratio (493.5 / 150000) and a time-and-cycle-dependent drug-to-antibody ratio DAR(tad, CYCLE) = 4 * (0.25 + 0.75 * exp(-beta * tad)) * (1 if CYCLE = 1 else Factor1). Body weight is included as a mechanistic covariate with a fixed allometric exponent of 0.75 on Dato-DXd linear clearance and estimated exponents on Dato-DXd volumes (paper Eq. 8-10) and on DXd CL/Vc (Eq. 14-15)."
  reference <- "Hong Y, Peigne S, Pan Y, Friberg Hietala S, McLaughlin A, Tajima N, Uema D, Zebger-Gong H, Tang Z, Zhou D, Abutarif M, Garimella T. Population Pharmacokinetic Analysis of Datopotamab Deruxtecan (Dato-DXd), a TROP2-Directed Antibody-Drug Conjugate, in Patients With Advanced Solid Tumors. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2149-2160. doi:10.1002/psp4.70118. PMID 41035281."
  vignette <- "Hong_2025_datopotamab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect with fixed exponent 0.75 on Dato-DXd CL_lin (Eq. 8); estimated exponents on Dato-DXd Vc (0.415, Eq. 9), Vp (0.311, Eq. 10), DXd CL (0.298, Eq. 14), and DXd Vc (0.530, Eq. 15). Reference 66 kg (Hong 2025 Table 1 / 2 reference subject).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on Dato-DXd CL_lin (exponent -0.306; Hong 2025 Table 1, Eq. 8). Reference 62 years (Hong 2025 reference subject).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effect on Dato-DXd CL_lin (1 + (-0.263) = 0.737 for females; Hong 2025 Table 1, Eq. 8) and Dato-DXd Vc (1 + (-0.160) = 0.840 for females; Eq. 9). Multiplicative effect on DXd Vc (1 + (-0.185) = 0.815 for females; Eq. 15). Males are the reference category.",
      source_name        = "SEX"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units; not g/dL). Power effect on Dato-DXd CL_lin (exponent -0.788, Hong 2025 Table 1 / Eq. 8) and DXd CL (exponent 0.343, Table 2 / Eq. 14). Reference 38 g/L (Hong 2025 reference subject).",
      source_name        = "ALB"
    ),
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on DXd CL (exponent -0.154; Hong 2025 Table 2 / Eq. 14). Reference 22 U/L (Hong 2025 reference subject). Note: paper text labels the unit 'g/L' which is implausible for AST; the equation form, the reference value 22, and the canonical clinical-PK convention all imply U/L (= IU/L) and that is what we use.",
      source_name        = "AST"
    ),
    TBILI = list(
      description        = "Baseline total serum bilirubin",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (US units; 1 mg/dL ~= 17.1 umol/L). Power effect on DXd CL (exponent -0.139; Hong 2025 Table 2 / Eq. 14). Reference 0.4 mg/dL (Hong 2025 reference subject).",
      source_name        = "TBIL"
    ),
    TUMSZ = list(
      description        = "Baseline sum of the longest dimension of target lesions (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on Dato-DXd Vmax (exponent 0.125; Hong 2025 Table 1 / Eq. 12). Reference 66 mm (Hong 2025 reference subject).",
      source_name        = "TUMSIZE"
    ),
    REGION_JAPAN = list(
      description        = "Study-site region indicator: 1 = Japan, 0 = otherwise (US, Europe, or Rest of World)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effect on Dato-DXd CL_lin only: 1 + (-0.219) = 0.781 for region Japan (Hong 2025 Table 1, Eq. 8). NOT used on DXd CL: Hong 2025 groups Japan with US into the DXd CL reference category, so REGION_JAPAN = 1 subjects fall into the DXd CL reference (REGION_EUROPE = 0 AND REGION_ROW = 0) regardless of REGION_JAPAN. Distinct from `RACE_JAPANESE` (subject ancestry).",
      source_name        = "REGJP"
    ),
    REGION_EUROPE = list(
      description        = "Study-site region indicator: 1 = Europe, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effect on DXd CL only: 1 + 0.240 = 1.240 for region Europe vs the US/Japan reference (Hong 2025 Table 2, Eq. 14). Not used on Dato-DXd parameters.",
      source_name        = "REGEU"
    ),
    REGION_ROW = list(
      description        = "Study-site region indicator: 1 = Rest of World (i.e., not US, Japan, or Europe), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative effect on DXd CL only: 1 + 0.196 = 1.196 for Rest of World vs the US/Japan reference (Hong 2025 Table 2, Eq. 14). Not used on Dato-DXd parameters.",
      source_name        = "REGRW"
    ),
    CYCLE = list(
      description        = "Treatment cycle number (1 = first dosing cycle, 2 = second, ..., integer count, time-varying across a multi-cycle Q3W treatment course)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = "Required for the DXd-formation DAR equation (Hong 2025 Eq. 13). DAR(tad, CYCLE) = 4 * (0.25 + 0.75 * exp(-beta * tad)) * (1 if CYCLE = 1, else Factor1 = 0.696). Increment CYCLE at the start of each new dosing cycle. Does not affect Dato-DXd disposition or DXd elimination.",
      source_name        = "CYCLE"
    )
  )

  population <- list(
    n_subjects     = 729L,
    n_studies      = 3L,
    age_range      = "Adults with advanced or metastatic solid tumors (NSCLC and breast cancer); reference subject 62 years per Hong 2025 covariate analysis.",
    weight_range   = "37.0 - 156 kg analysis range (Hong 2025 Methods / Body Weight-Based Dosing); reference 66 kg.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Multi-regional enrollment including Japan, US, Europe, and Rest of the World; race per se was not retained as a covariate in the final model.",
    disease_state  = "Advanced or metastatic non-small-cell lung cancer (NSCLC) and breast cancer (BC).",
    dose_range     = "Dato-DXd IV every 3 weeks. Phase 1 TROPION-PanTumor01 (TP01) dose range 0.27 - 10 mg/kg; Phase 2 TL05 and Phase 3 TL01 used the labeled 6 mg/kg Q3W regimen. Approved labeled dosing is 6 mg/kg Q3W with a flat dose cap of 540 mg for body weight >= 90 kg (Datroway US prescribing information; Hong 2025 Body Weight-Based Dosing section).",
    regions        = "Multi-regional: Japan, US, Europe, and Rest of the World.",
    n_observations = "9036 Dato-DXd PK observations + 9012 DXd PK observations from 729 patients across three studies (Hong 2025 Methods).",
    studies        = "Phase 1 TROPION-PanTumor01 (TP01), Phase 2 TL05, and Phase 3 TL01.",
    reference_subject_dato_dxd = "66 kg male, age 62 years, albumin 38 g/L, tumor size 66 mm, region != Japan (Hong 2025 Figure 4 caption / Table 1 reference).",
    reference_subject_dxd      = "66 kg male, region US, albumin 38 g/L, AST 22 U/L, total bilirubin 0.4 mg/dL (Hong 2025 Figure 4 caption / Table 2 reference).",
    notes          = "Pooled analysis dataset (Hong 2025 Methods). Final model fit jointly to all 729 subjects across the 3 studies, excluding observations before the first Dato-DXd administration. Anti-drug antibodies were tested but had only a small (6.9%) effect on CL_lin and were dropped from the final model. Renal-impairment (creatinine clearance) and broad hepatic-impairment (HI) classifications were not statistically significant covariates on Dato-DXd disposition; AST and total bilirubin were retained on DXd CL only."
  )

  ini({
    # ----- Dato-DXd structural parameters, typical values for the reference
    # subject (66 kg male, age 62, albumin 38 g/L, tumor 66 mm, non-Japan).
    # Concentrations are in mg/L (= ug/mL) so that dose in mg / V in L gives
    # ug/mL directly, which matches the bioanalytical assay scale.
    lcl       <- log(0.386);  label("Dato-DXd linear clearance at reference covariates (CL_linDatoDXd, L/day)")  # Hong 2025 Table 1, Eq. 8
    lvc       <- log(3.06);   label("Dato-DXd central volume at reference (V_cDatoDXd, L)")                       # Hong 2025 Table 1, Eq. 9
    lq        <- log(0.422);  label("Dato-DXd intercompartmental clearance (Q_DatoDXd, L/day)")                   # Hong 2025 Table 1
    lvp       <- log(2.88);   label("Dato-DXd peripheral volume at reference (V_pDatoDXd, L)")                    # Hong 2025 Table 1, Eq. 10
    lvmax     <- log(8.41);   label("Dato-DXd Michaelis-Menten Vmax at reference tumor size (mg/day)")           # Hong 2025 Table 1, Eq. 12: 8410 ug/day = 8.41 mg/day
    lkm       <- log(4.49);   label("Dato-DXd Michaelis-Menten Km (mg/L = ug/mL)")                                # Hong 2025 Table 1, Eq. 11: 4490 ng/mL = 4.49 ug/mL

    # ----- DXd structural parameters, typical values for the reference subject
    # (66 kg male, region US, albumin 38, AST 22 U/L, total bilirubin 0.4 mg/dL).
    # Paper reports CL_DXd in L/h (Table 2, Eq. 14); converted to L/day for
    # consistency with all other rate constants in this model.
    lcl_dxd   <- log(2.66 * 24);  label("DXd linear clearance at reference covariates (CL_DXd, L/day)")  # Hong 2025 Table 2, Eq. 14: 2.66 L/h * 24 h/day = 63.84 L/day
    lvc_dxd   <- log(25.1);       label("DXd central volume at reference (V_cDXd, L)")                    # Hong 2025 Table 2, Eq. 15
    lfactor1  <- log(0.696);      label("DAR scaling factor for cycle 2 and later (Factor1, unitless)")  # Hong 2025 Table 2, Eq. 13
    lbeta     <- log(0.259);      label("Within-cycle DAR exponential decay rate constant (beta, 1/day)") # Hong 2025 Table 2, Eq. 13

    # ----- Covariate effects on Dato-DXd parameters (Hong 2025 Table 1 / Eq. 8-12).
    # Continuous covariates enter as power models (CovEff = (Cov / Cov_ref)^theta;
    # paper Eq. 4); categorical covariates enter as fractional-difference models
    # (CovEff = 1 + theta for non-reference; paper Eq. 6).
    e_wt_cl          <- fix(0.75); label("Power exponent of WT on Dato-DXd CL_lin (unitless, fixed)")             # Hong 2025 Table 1: 0.750 fixed (estimated 0.80; Methods)
    e_wt_vc          <-  0.415;    label("Power exponent of WT on Dato-DXd Vc (unitless)")                        # Hong 2025 Table 1
    e_wt_vp          <-  0.311;    label("Power exponent of WT on Dato-DXd Vp (unitless)")                        # Hong 2025 Table 1
    e_age_cl         <- -0.306;    label("Power exponent of AGE on Dato-DXd CL_lin (unitless)")                   # Hong 2025 Table 1
    e_alb_cl         <- -0.788;    label("Power exponent of ALB on Dato-DXd CL_lin (unitless)")                   # Hong 2025 Table 1
    e_japan_cl       <- -0.219;    label("Fractional shift of region Japan on Dato-DXd CL_lin (unitless)")        # Hong 2025 Table 1: 1 + (-0.219) for Japan
    e_sexf_cl        <- -0.263;    label("Fractional shift of female sex on Dato-DXd CL_lin (unitless)")          # Hong 2025 Table 1: 1 + (-0.263) for female
    e_sexf_vc        <- -0.160;    label("Fractional shift of female sex on Dato-DXd Vc (unitless)")              # Hong 2025 Table 1: 1 + (-0.160) for female
    e_tumsz_vmax     <-  0.125;    label("Power exponent of TUMSZ on Dato-DXd Vmax (unitless)")                   # Hong 2025 Table 1, Eq. 12

    # ----- Covariate effects on DXd parameters (Hong 2025 Table 2 / Eq. 14-15).
    e_wt_cl_dxd      <- fix(0.298); label("Power exponent of WT on DXd CL (unitless, fixed)")                     # Hong 2025 Table 2: paper "estimated then fixed"
    e_wt_vc_dxd      <- fix(0.530); label("Power exponent of WT on DXd Vc (unitless, fixed)")                     # Hong 2025 Table 2: paper "estimated then fixed"
    e_alb_cl_dxd     <-  0.343;     label("Power exponent of ALB on DXd CL (unitless)")                           # Hong 2025 Table 2
    e_ast_cl_dxd     <- -0.154;     label("Power exponent of AST on DXd CL (unitless)")                           # Hong 2025 Table 2
    e_tbili_cl_dxd   <- -0.139;     label("Power exponent of TBILI on DXd CL (unitless)")                         # Hong 2025 Table 2
    e_eu_cl_dxd      <-  0.240;     label("Fractional shift of region Europe on DXd CL vs US/Japan reference (unitless)") # Hong 2025 Table 2: 1 + 0.240 for Europe
    e_row_cl_dxd     <-  0.196;     label("Fractional shift of region RoW on DXd CL vs US/Japan reference (unitless)")     # Hong 2025 Table 2: 1 + 0.196 for RoW
    e_sexf_vc_dxd    <- -0.185;     label("Fractional shift of female sex on DXd Vc (unitless)")                  # Hong 2025 Table 2: 1 + (-0.185) for female

    # ----- IIV. Paper reports IIV as %CV on log-normal parameters (Eq. 2);
    # convert to omega^2 = log(CV^2 + 1).
    etalcl      ~ 0.071375  # Dato-DXd CL_lin 27.2% CV; log(0.272^2 + 1) = 0.071375  -- Hong 2025 Table 1
    etalvc      ~ 0.020807  # Dato-DXd Vc     14.5% CV; log(0.145^2 + 1) = 0.020807  -- Hong 2025 Table 1
    etalq       ~ 0.092325  # Dato-DXd Q      31.1% CV; log(0.311^2 + 1) = 0.092325  -- Hong 2025 Table 1
    etalvp      ~ 0.092893  # Dato-DXd Vp     31.2% CV; log(0.312^2 + 1) = 0.092893  -- Hong 2025 Table 1
    etalvmax    ~ 0.036201  # Dato-DXd Vmax   19.2% CV; log(0.192^2 + 1) = 0.036201  -- Hong 2025 Table 1

    etalcl_dxd  ~ 0.094033  # DXd CL          31.4% CV; log(0.314^2 + 1) = 0.094033  -- Hong 2025 Table 2
    etalvc_dxd  ~ 0.123782  # DXd Vc          36.3% CV; log(0.363^2 + 1) = 0.123782  -- Hong 2025 Table 2

    # ----- Residual error. Paper uses an additive-on-log-scale RUV with an
    # exponential IIV term on RUV (Eq. 3). The additive-on-log-scale SD maps
    # directly to a proportional SD in nlmixr2's linear-space residual model.
    # The IIV-on-RUV term (Dato-DXd 45.8% CV; DXd 29.2% CV) is a non-standard
    # subject-level scaling of the residual SD that nlmixr2 does not natively
    # support; this implementation uses the typical-subject residual SD only.
    # See vignette Assumptions and deviations.
    propSd     <- 0.121;  label("Proportional residual error on Dato-DXd Cc (fraction)")  # Hong 2025 Table 1: additive RUV on log scale, CV 0.121
    propSd_dxd   <- 0.283;  label("Proportional residual error on DXd Cc_dxd (fraction)")   # Hong 2025 Table 2: additive RUV on log scale, CV 0.283
  })
  model({
    # ----- Individual Dato-DXd parameters (Hong 2025 Eq. 8-12). Continuous
    # covariates are normalized to the reference values 66 kg / 62 yr / 38 g/L
    # and enter as power terms; categorical covariates (Japan, female) enter
    # multiplicatively as 1 + theta for the non-reference category.
    cl <- exp(lcl + etalcl) *
      (WT  / 66)^e_wt_cl *
      (AGE / 62)^e_age_cl *
      (ALB / 38)^e_alb_cl *
      (1 + e_japan_cl * REGION_JAPAN) *
      (1 + e_sexf_cl  * SEXF)
    vc <- exp(lvc + etalvc) *
      (WT / 66)^e_wt_vc *
      (1 + e_sexf_vc * SEXF)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) *
      (WT / 66)^e_wt_vp
    vmax   <- exp(lvmax   + etalvmax) *
      (TUMSZ / 66)^e_tumsz_vmax
    km     <- exp(lkm)

    # ----- Individual DXd parameters (Hong 2025 Eq. 14-15). REGION_EUROPE and
    # REGION_ROW are mutually exclusive (a subject is in at most one); when
    # both are 0 the subject is in the US/Japan reference group.
    cl_dxd <- exp(lcl_dxd + etalcl_dxd) *
      (WT    / 66)^e_wt_cl_dxd *
      (ALB   / 38)^e_alb_cl_dxd *
      (AST   / 22)^e_ast_cl_dxd *
      (TBILI / 0.4)^e_tbili_cl_dxd *
      (1 + e_eu_cl_dxd  * REGION_EUROPE +
           e_row_cl_dxd * REGION_ROW)
    vc_dxd <- exp(lvc_dxd + etalvc_dxd) *
      (WT / 66)^e_wt_vc_dxd *
      (1 + e_sexf_vc_dxd * SEXF)

    factor1 <- exp(lfactor1)
    beta    <- exp(lbeta)

    # ----- Dato-DXd plasma concentration (mg/L = ug/mL) and instantaneous
    # elimination rates. The Michaelis-Menten elimination rate has units of
    # mass / time; with Cc in mg/L = ug/mL and km in the same units,
    # vmax * C / (km + C) is in mg/day.
    Cc        <- central / vc
    rate_lin  <- cl * Cc
    rate_mm   <- vmax * Cc / (km + Cc)
    rate_elim <- rate_lin + rate_mm

    # ----- Time-and-cycle dependent DAR (Hong 2025 Eq. 13). tad() is the time
    # since the most-recent dose (days). CYCLE = 1 for the first dose; for any
    # subsequent cycle, the DAR is scaled by Factor1 (= 0.696). The within-
    # cycle exponential decay e^(-beta * tad) is the same in every cycle.
    factor_dar <- factor1
    if (CYCLE < 2) factor_dar <- 1.0
    dar <- 4.0 * (0.25 + 0.75 * exp(-beta * tad())) * factor_dar

    # ----- DXd formation rate (mg/day): the total Dato-DXd elimination rate
    # is converted to a DXd mass rate using the per-antibody DAR (mol DXd /
    # mol Dato-DXd) and the molecular-weight ratio (493.5 g/mol DXd over
    # ~150,000 g/mol Dato-DXd). Hong 2025 DXd Base Model Development +
    # Methods (493.5 / 150000 = 3.29e-3).
    rate_dxd_in <- rate_elim * dar * (493.5 / 150000)

    # ----- ODE system. Dato-DXd dosing (IV bolus or short infusion) goes
    # directly into central; DXd has no exogenous dose.
    d/dt(central)     <- -rate_elim - q * Cc + (q / vp) * peripheral1
    d/dt(peripheral1) <-              q * Cc - (q / vp) * peripheral1
    d/dt(central_dxd) <-  rate_dxd_in - cl_dxd * central_dxd / vc_dxd

    # ----- Observations: Dato-DXd plasma concentration (mg/L = ug/mL) and
    # DXd plasma concentration (mg/L = ug/mL). DXd is conventionally reported
    # in ng/mL in the paper; multiply Cc_dxd by 1000 outside the model to
    # obtain ng/mL.
    Cc_dxd <- central_dxd / vc_dxd

    Cc     ~ prop(propSd)
    Cc_dxd ~ prop(propSd_dxd)
  })
}
