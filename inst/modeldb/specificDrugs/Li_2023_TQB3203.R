Li_2023_TQB3203 <- function() {
  description <- paste(
    "Three-compartment population PK model with first-order elimination for",
    "TQ-B3203, a novel semisynthetic camptothecin (topoisomerase I inhibitor)",
    "derivative, following a 90-min intravenous infusion of TQ-B3203 liposome",
    "injection (TLI) in Chinese patients with advanced solid tumors. Built on",
    "316 plasma concentrations from 15 patients over 25 treatment episodes in a",
    "2-45 mg/m^2 dose-escalation phase I trial (NCT03447145). Parameterized as",
    "central volume V1, shallow peripheral volume V2, deep peripheral volume",
    "V3, elimination clearance CL, and inter-compartmental clearances CL2",
    "(central <-> shallow) and CL3 (central <-> deep). Three covariates entered",
    "the final model as median-normalized power functions: body mass index BMI",
    "and direct bilirubin DBIL on CL, DBIL and lean body mass LBM on CL2, and",
    "LBM on V1 and V2. Higher BMI raises CL while higher DBIL lowers it,",
    "consistent with TQ-B3203 being excreted unchanged into faeces via bile so",
    "that impaired biliary excretion (raised DBIL) reduces clearance.",
    "Inter-individual variability is exponential on CL, V3, CL2 and CL3, and",
    "the residual error is log-additive (encoded as lnorm). Only parent",
    "TQ-B3203 was modelled; conversion to SN-38 was <5% (paper Discussion)."
  )
  reference <- paste(
    "Li X, Bo Y, Yin H, Liu X, Li X, Yang F. Population pharmacokinetic",
    "analysis of TQ-B3203 following intravenous administration of TQ-B3203",
    "liposome injection in Chinese patients with advanced solid tumors.",
    "Front Pharmacol. 2023 Jan 16;14:1102244. doi:10.3389/fphar.2023.1102244"
  )
  vignette <- "Li_2023_TQB3203"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Doses were prescribed in mg/m^2 (2-45 mg/m^2) and must be converted to an
  # absolute mg amount with the subject's BSA before use in an event record.
  # State amounts are therefore mg and Cc is scaled by 1000 to reach ng/mL.
  compartmentData <- list(
    central = list(
      analyte = "TQ-B3203", units = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "TQ-B3203", units = "mg", specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "TQ-B3203", units = "mg", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-normalized power effect on CL, (BMI / 23.44)^0.78, per Li 2023",
        "Eq. 10 and Eq. 7. The normalizing constant 23.44 kg/m^2 is the",
        "population median (Li 2023 Table 1 and Results 3.3 narrative).",
        "Positive exponent: higher BMI is associated with higher CL, which the",
        "authors attribute to obesity-enhanced phase I / phase II metabolism",
        "(paper Discussion). Trial eligibility restricted BMI to 18.5-26",
        "kg/m^2 at screening, and the observed range was 18.64-28.97 kg/m^2,",
        "so the effect is not supported outside roughly 19-29 kg/m^2.",
        "Time-fixed at baseline."
      ),
      source_name        = "BMI"
    ),
    DBIL = list(
      description        = paste(
        "Direct (conjugated) serum bilirubin at baseline, a marker of biliary",
        "excretion function"
      ),
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-normalized power effects on CL, (DBIL / 2.82)^-0.24 (Li 2023",
        "Eq. 10, Table 4), and on CL2, (DBIL / 2.82)^-1.77 (Li 2023 Table 4).",
        "The normalizing constant 2.82 umol/L is the population median (Li 2023",
        "Table 1 and Results 3.3 narrative). The paper reports DBIL in SI",
        "umol/L, which is already the canonical unit for this column, so no",
        "conversion is applied. Negative exponent on CL: TQ-B3203 is excreted",
        "unchanged into faeces through bile, so impaired biliary excretion",
        "raises DBIL and lowers TQ-B3203 clearance (paper Discussion).",
        "Observed range 1.5-7.9 umol/L. Time-fixed at baseline."
      ),
      source_name        = "DBIL"
    ),
    LBM = list(
      description        = "Lean body mass (the paper's lean body weight, LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-normalized power effects on V1, (LBM / 49.53)^1.18 (Li 2023",
        "Eq. 9, Table 4); on V2, (LBM / 49.53)^-1.41 (Table 4); and on CL2,",
        "(LBM / 49.53)^-2.55 (Table 4). The normalizing constant 49.53 kg is",
        "the population median (Li 2023 Table 1 and Results 3.3 narrative).",
        "The paper's column name is LBW (lean body weight), a registered alias",
        "of the canonical LBM; same quantity in kg with no value",
        "transformation. Li 2023 does NOT state which body-composition formula",
        "produced LBW -- it only cites Park et al. 2018 (Br J Anaesth",
        "121:559-566) as the reason LBW and body-fat percentage were screened,",
        "given the high lipophilicity of the drug. Because the James / Boer /",
        "Hume formulae differ materially at a given height and weight, a",
        "downstream user must supply LBM on a formula consistent with the",
        "cohort in Li 2023 Table 1 (median 49.53 kg, range 32.09-59.40 kg for",
        "a cohort of median weight 64.0 kg and height 164 cm). Time-fixed at",
        "baseline."
      ),
      source_name        = "LBW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = paste(
        "Screened (Li 2023 Methods 2.5.2, Table 1: median 57 years, range",
        "31-70) but not retained in the final model."
      )
    ),
    WT = list(
      description = "Total body weight", units = "kg", type = "continuous",
      notes = paste(
        "Screened (Table 1: median 64.0 kg, range 47.9-80.0) but dropped",
        "before stepwise selection because it correlated (r > 0.5) with the",
        "retained body-size covariates (Li 2023 Results 3.3 collinearity",
        "check)."
      )
    ),
    HT = list(
      description = "Height", units = "cm", type = "continuous",
      notes = "Screened (Table 1: median 164 cm, range 148-178); not retained."
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula)", units = "m^2",
      type = "continuous",
      notes = paste(
        "Screened as a covariate and not retained (Li 2023 Methods 2.5.2).",
        "BSA nevertheless remains operationally necessary to convert the",
        "prescribed mg/m^2 dose into an absolute mg amount; population median",
        "1.678 m^2, range 1.420-1.938 m^2 (Table 1)."
      )
    ),
    IBW = list(
      description = "Ideal body weight (Devine formula, Li 2023 Eq. 5 / 5a)",
      units = "kg", type = "continuous",
      notes = "Screened (Table 1: median 60.50 kg); not retained."
    ),
    WT_ADJUSTED = list(
      description = "Adjusted body weight (Li 2023 Eq. 6 / 6a)",
      units = "kg", type = "continuous",
      notes = "Screened (Table 1: median 61.90 kg); not retained."
    ),
    BODYFAT_PCT = list(
      description = "Body fat percentage", units = "%", type = "continuous",
      notes = paste(
        "Screened alongside LBW because of the drug's high lipophilicity",
        "(Li 2023 Methods 2.5.2, citing Park 2018); Table 1 median 23.39%.",
        "Not retained."
      )
    ),
    CRE = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      notes = "Screened (Table 1: median 60 umol/L); not retained."
    ),
    CRCL = list(
      description = paste(
        "Endogenous creatinine clearance (Cockcroft-Gault, Li 2023 Eq. 4 / 4a)",
        "and its adjusted-weight variant"
      ),
      units = "mL/min", type = "continuous",
      notes = paste(
        "Only the adjusted-CLcr variant survived the collinearity check and",
        "entered stepwise screening (Li 2023 Results 3.3); it was not",
        "retained in the final model. Table 1 reports CLcr median 109.2 and",
        "adjusted CLcr median 101.9 (the paper labels both 'mg/dL', which is",
        "a unit typo -- creatinine clearance is a flow, conventionally",
        "mL/min)."
      )
    ),
    TBILI = list(
      description = "Total serum bilirubin", units = "umol/L",
      type = "continuous",
      notes = paste(
        "Screened but dropped in the univariate collinearity check against",
        "DBIL (Li 2023 Discussion); Table 1 median 9 umol/L."
      )
    ),
    IBIL = list(
      description = "Indirect (unconjugated) serum bilirubin",
      units = "umol/L", type = "continuous",
      notes = paste(
        "Screened but dropped in the collinearity check against DBIL",
        "(Li 2023 Discussion); Table 1 median 6 umol/L."
      )
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase", units = "IU/L",
      type = "continuous",
      notes = paste(
        "Retained through the collinearity check and entered stepwise",
        "screening, but the authors 'did not find evidence of ALT as a",
        "significant covariate on CL' (Li 2023 Discussion); Table 1 median",
        "11.3 IU/L."
      )
    ),
    AST = list(
      description = "Baseline aspartate transaminase", units = "IU/L",
      type = "continuous",
      notes = paste(
        "Screened but dropped in the collinearity check (Li 2023 Discussion);",
        "Table 1 median 16.7 IU/L."
      )
    ),
    SEXF = list(
      description = "Sex", units = "(binary)", type = "binary",
      notes = paste(
        "Screened and found to influence the same PK parameters as LBW and",
        "DBIL, but dropped before stepwise selection because it was",
        "correlated with them (Li 2023 Figure 4, Results 3.3). Cohort was",
        "33.3% female (Table 1)."
      )
    ),
    SNP_UGT1A1_RS8175347 = list(
      description = "UGT1A1*28 promoter TA-repeat genotype", units = "(genotype)",
      type = "categorical",
      notes = paste(
        "Screened (Table 1: 14/15 TA(6)/TA(6), 1/15 TA(6)/TA(7)); not",
        "retained. Effectively non-informative at this sample size."
      )
    ),
    SNP_UGT1A1_RS4148323 = list(
      description = "UGT1A1*6 (211G>A) genotype", units = "(genotype)",
      type = "categorical",
      notes = paste(
        "Screened (Table 1: 10/15 211G/G, 5/15 211G/A). Entered the full",
        "model as an effect on V2 at forward step 6, but was eliminated in",
        "backward elimination (dOFV +4.27, p > 0.01) and is absent from the",
        "final model (Li 2023 Table 3 step 9, Table 4)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    n_observations = 316L,
    n_episodes     = 25L,
    age_range      = "31-70 years (median 57, IQR 44-65)",
    age_median     = "57 years",
    weight_range   = "47.9-80.0 kg (median 64.0, IQR 57.5-68.0)",
    weight_median  = "64.0 kg",
    height_range   = "148-178 cm (median 164, IQR 160-173)",
    bmi_range      = "18.64-28.97 kg/m^2 (median 23.44, IQR 21.05-24.91)",
    lbm_range      = "32.09-59.40 kg (median 49.53, IQR 38.36-55.45)",
    bsa_range      = "1.420-1.938 m^2 (median 1.678, IQR 1.600-1.841)",
    sex_female_pct = 33.3,
    race_ethnicity = "100% Chinese (single-country trial; Table 1 reports no race strata).",
    hepatic_function = paste(
      "Normal primary organ function was an eligibility criterion. Baseline",
      "total bilirubin 2.9-25.8 umol/L (median 9), direct bilirubin 1.5-7.9",
      "umol/L (median 2.82), ALT 3-24 IU/L (median 11.3), AST 11.2-29 IU/L",
      "(median 16.7) -- Table 1."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 27.35-157.3 (median 109.2);",
      "serum creatinine 38-204 umol/L (median 60) -- Table 1."
    ),
    disease_state  = paste(
      "Clearly diagnosed advanced solid tumors, ECOG performance status 0-1,",
      "life expectancy > 3 months, no prior camptothecin-analog therapy."
    ),
    dose_range     = paste(
      "2-45 mg/m^2 TQ-B3203 liposome injection as a single 90-min IV infusion",
      "on day 1 of cycle 1 and day 22 of cycle 2 (dose levels 2, 4, 6, 10, 14,",
      "30 and 45 mg/m^2; Table 2)."
    ),
    regions        = "China (multi-centre; Peking University Cancer Hospital lead site).",
    notes          = paste(
      "Baseline demographics from Li 2023 Table 1; NCA summaries by dose level",
      "and cycle from Table 2. Estimation used Phoenix NLME 8.3 with FOCE-ELS",
      "(not NONMEM). Sampling was pre-dose, 45 and 90 min after infusion start,",
      "then 0.5, 1, 2, 4, 8, 12, 24, 48, 72 and 96 h after the end of the",
      "infusion, so the last sample is 97.5 h after the start of dosing. The",
      "TQ-B3203 assay was linear over 0.5-500 ng/mL. Inter-occasion",
      "variability on CL and V1 was tested and rejected (Results 3.3), and no",
      "significant cycle-1 versus cycle-2 difference was found, so this model",
      "carries no IOV. Because the covariate ranges are narrow (a phase I",
      "dose-escalation cohort of 15 patients), the authors explicitly caution",
      "that extrapolation beyond them is not supported (Discussion)."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- Li 2023 Table 4, "Final model / Estimate
    # (%RSE)" column. Every value carries a bootstrap 95% CI that contains
    # the point estimate, and all RSEs are < 38%, so all six are estimated
    # (none fixed).
    # ------------------------------------------------------------------------
    lvc  <- log(4.81);   label("Central volume V1 (L)")
    # Li 2023 Table 4: V1 = 4.81 L (RSE 8.47%, bootstrap 95% CI 4.01-5.61);
    # restated in Results 3.3 and Eq. 9 and in the Discussion.
    lvp  <- log(24.44);  label("Shallow peripheral volume V2 (L)")
    # Li 2023 Table 4: V2 = 24.44 L (RSE 5.94%, bootstrap 95% CI 21.58-27.30)
    lvp2 <- log(27.98);  label("Deep peripheral volume V3 (L)")
    # Li 2023 Table 4: V3 = 27.98 L (RSE 6.69%, bootstrap 95% CI 24.30-31.66)
    lcl  <- log(3.97);   label("Elimination clearance CL (L/h)")
    # Li 2023 Table 4: CL = 3.97 L/h (RSE 4.85%, bootstrap 95% CI 3.59-4.35);
    # restated in Results 3.3 and Eq. 10 and in the Abstract and Discussion.
    lq   <- log(1.95);   label("Inter-compartmental CL central <-> shallow peripheral, CL2 (L/h)")
    # Li 2023 Table 4: CL2 = 1.95 L/h (RSE 20.48%, bootstrap 95% CI 1.17-2.74)
    lq2  <- log(10.58);  label("Inter-compartmental CL central <-> deep peripheral, CL3 (L/h)")
    # Li 2023 Table 4: CL3 = 10.58 L/h (RSE 11.28%, bootstrap 95% CI 8.23-12.92)

    # ------------------------------------------------------------------------
    # COVARIATE EFFECTS -- Li 2023 Table 4 covariate rows. All are exponents
    # of a median-normalized power function, (Cov / Cov_median)^theta, per
    # Li 2023 Eq. 7. Population medians used as the normalizing constants are
    # BMI 23.44 kg/m^2, DBIL 2.82 umol/L and LBM (paper's LBW) 49.53 kg
    # (Li 2023 Results 3.3 narrative, matching Table 1).
    #
    # Table 4 labels these rows with the units of the AFFECTED parameter
    # (e.g. "BMI on CL (L/h)"), which is a labelling artifact -- a power-model
    # exponent is unitless. The exponents are used here exactly as printed.
    # ------------------------------------------------------------------------
    e_bmi_cl  <-  0.78;  label("Exponent of (BMI / 23.44) on CL (unitless)")
    # Li 2023 Table 4 "BMI on CL" = 0.78 (RSE 22.08%, bootstrap 95% CI
    # 0.44-1.12); appears as the exponent 0.78 in Eq. 10 and is named in
    # Results 3.3 ("0.78 ... between BMI and CL").
    e_dbil_cl <- -0.24;  label("Exponent of (DBIL / 2.82) on CL (unitless)")
    # Li 2023 Table 4 "DBIL on CL" = -0.24 (RSE -36.01%, bootstrap 95% CI
    # -0.42 to -0.07); appears as the exponent -0.24 in Eq. 10 and is named
    # in Results 3.3 ("-.24 ... between DBIL and CL").
    e_dbil_q  <- -1.77;  label("Exponent of (DBIL / 2.82) on CL2 (unitless)")
    # Li 2023 Table 4 "DBIL on CL2" = -1.77 (RSE -19.55%, bootstrap 95% CI
    # -2.46 to -1.09); retained in the final model per Table 3 step 9
    # ("CL2-DBIL-LBW").
    e_lbm_q   <- -2.55;  label("Exponent of (LBM / 49.53) on CL2 (unitless)")
    # Li 2023 Table 4 "LBW on CL2" = -2.55 (RSE -27.03%, bootstrap 95% CI
    # -3.90 to -1.19); retained per Table 3 step 9 ("CL2-DBIL-LBW").
    e_lbm_vc  <-  1.18;  label("Exponent of (LBM / 49.53) on V1 (unitless)")
    # Li 2023 Table 4 "LBW on V1" = 1.18 (RSE 37.01%, bootstrap 95% CI
    # 0.32-2.05); appears as the exponent 1.18 in Eq. 9 and is named in
    # Results 3.3 ("1.18 ... between LBW and V1").
    e_lbm_vp  <- -1.41;  label("Exponent of (LBM / 49.53) on V2 (unitless)")
    # Li 2023 Table 4 "LBW on V2" = -1.41 (RSE -19.49%, bootstrap 95% CI
    # -1.95 to -0.87); retained per Table 3 step 9 ("V2-LBW").

    # ------------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Li 2023 Table 4, "Inter-individual
    # variability" section. The paper used an exponential IIV model
    # (Eq. 1: P_ij = theta_i * exp(eta_ij)) and the tabulated quantities are
    # labelled "omega^2", i.e. VARIANCES on the log scale, so they are used
    # directly with no CV back-transform. Table 4 footnote a calls them
    # "variance of inter-individual variability".
    #
    # Table 4 reports exactly four IIV terms -- on CL, V3, CL2 and CL3. No
    # omega^2 is reported for V1 or V2, and the paper reports no IIV
    # covariance block ("the correlation between IIV of pharmacokinetic
    # parameters should be clarified to construct a covariance model",
    # Methods 2.5.2, but no covariance is tabulated), so all four enter as
    # independent diagonal terms.
    #
    # DEVIATION: Eq. 9 prints V1 with an exp(eta_V1) term, but Table 4
    # reports no omega^2 for V1. Rather than invent a variance, etalvc is
    # declared and held at fixed(0) -- this keeps the printed Eq. 9 structure
    # visible and auditable while contributing no variability, so a
    # downstream re-fit can free it if the value is ever recovered. V2 is
    # not given an eta at all, because no printed equation or table row
    # anywhere in the paper shows one. See the vignette Assumptions and
    # deviations section.
    # ------------------------------------------------------------------------
    etalvc  ~ fixed(0);
    # Li 2023 Eq. 9 prints V1 = 4.81 * (LBW/49.53)^1.18 * exp(eta_V1), but
    # Table 4 reports no omega^2 for V1; held at 0 rather than invented.
    etalcl  ~ 0.043;  # Li 2023 Table 4: omega^2 CL = 0.043 (RSE 36.00%, bootstrap 95% CI 0.013-0.073)
    etalvp2 ~ 0.117;  # Li 2023 Table 4: omega^2 V3 = 0.117 (RSE 29.88%, bootstrap 95% CI 0.048-0.185)
    etalq   ~ 0.573;  # Li 2023 Table 4: omega^2 CL2 = 0.573 (RSE 32.60%, bootstrap 95% CI 0.203-0.944)
    etalq2  ~ 0.290;  # Li 2023 Table 4: omega^2 CL3 = 0.290 (RSE 32.50%, bootstrap 95% CI 0.103-0.477)

    # ------------------------------------------------------------------------
    # RESIDUAL ERROR -- Li 2023 Table 4 "Residual variability (sigma)
    # stdev0" = 0.200. Results 3.3 states that "a log-additive model" was
    # selected "to describe intra-individual variability, so the
    # pharmacokinetic parameters were estimated with the natural
    # logarithm-transformed (log) plasma TQ-B3203 concentration data".
    #
    # A Phoenix NLME log-additive residual is log(C_obs) = log(C_pred) + eps
    # with eps ~ N(0, sigma^2), i.e. exactly a log-normal residual. That maps
    # to nlmixr2's `lnorm(expSd)`, NOT to `prop(propSd)`: the SD applies on
    # the log scale. Same encoding decision as Sano_2023_fesoterodine.R.
    # Consistency check: Results 3.3 quotes the base-model log-additive error
    # as 22.4% and 21.4% with IOV, bracketing the final 0.200 = 20.0%.
    # ------------------------------------------------------------------------
    expSd <- 0.200;  label("Log-additive (log-normal) residual SD on the log-concentration scale")
    # Li 2023 Table 4: sigma stdev0 = 0.200 (RSE 10.30%, bootstrap 95% CI 0.159-0.240)
  })

  model({
    # 1. Individual parameters. Each covariate enters as a median-normalized
    #    power function (Li 2023 Eq. 7); the normalizing medians are the
    #    population medians quoted in Results 3.3.
    #
    #    V1 (L)  = 4.81 * (LBW / 49.53)^1.18                    -- Eq. 9
    #    CL (L/h) = 3.97 * (BMI / 23.44)^0.78 * (DBIL / 2.82)^-0.24 -- Eq. 10
    #    V2, CL2 covariate forms are from Table 4 / Table 3 step 9; the paper
    #    prints closed-form equations only for V1 and CL.
    vc  <- exp(lvc + etalvc) * (LBM / 49.53)^e_lbm_vc
    vp  <- exp(lvp) * (LBM / 49.53)^e_lbm_vp
    vp2 <- exp(lvp2 + etalvp2)
    cl  <- exp(lcl + etalcl) * (BMI / 23.44)^e_bmi_cl * (DBIL / 2.82)^e_dbil_cl
    q   <- exp(lq + etalq) * (DBIL / 2.82)^e_dbil_q * (LBM / 49.53)^e_lbm_q
    q2  <- exp(lq2 + etalq2)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system -- three compartments, first-order elimination from
    #    central (Li 2023 Results 3.3). Dosing is a 90-min IV infusion into
    #    `central`; supply it as dur = 1.5 (h) or rate = amt / 1.5.
    d/dt(central)     <- -(kel + k12 + k13) * central +
      k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Observation. States are mg and volumes are L, so central / vc is
    #    mg/L = ug/mL; multiply by 1000 to reach the paper's ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
