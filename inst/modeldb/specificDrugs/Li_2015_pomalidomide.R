Li_2015_pomalidomide <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag for oral pomalidomide in healthy participants and patients with relapsed/refractory multiple myeloma (Li 2015). Apparent clearance and central volume are comparable between the two populations, but multiple myeloma raises apparent peripheral volume 8.46-fold and apparent intercompartmental clearance 3.71-fold and shortens the lag time; covariates are body weight and total serum protein on V2/F and sex on CL/F, and the log-scale residual error is population-specific."
  reference <- "Li Y, Xu Y, Liu L, Wang X, Palmisano M, Zhou S. Population pharmacokinetics of pomalidomide. The Journal of Clinical Pharmacology. 2015;55(5):563-572. doi:10.1002/jcph.455"
  vignette <- "Li_2015_pomalidomide"
  # Two population-specific log-scale residual SDs (Table 2 'sigma^2 (HNP)' and
  # 'sigma^2 (MM patients)'); the canonical bare `expSd` has only one slot.
  paper_specific_residual_sds <- c("expSdHnp", "expSdMm")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "pomalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "pomalidomide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pomalidomide", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight; power effect on V2/F with reference 78.3 kg taken from the final covariate equation printed in Results 'Covariate Analysis'. Note the demographics table (Table 1) reports a cohort median of 78.1 kg; the equation's own 78.3 kg is the centering constant the reported V2/F of 58.3 L is conditioned on and is the value used here.",
      source_name        = "WT"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline total serum protein; linear centered effect on V2/F with centering constant 73.0 g/L taken from the final covariate equation printed in Results 'Covariate Analysis'. Table 1 reports a cohort median of 75.0 g/L and the Discussion quotes an observed range of 56 to 148 g/L; as with WT the equation's own constant is used. The paper argues the effect is a surrogate for disease stage rather than a binding-capacity mechanism, because multiple myeloma raises serum protein and higher stage correlates with higher protein.",
      source_name        = "TPT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effect on CL/F via (1 + e_sexf_cl * SEXF) with e_sexf_cl = -0.234, i.e. female CL/F is 23.4% lower than male. The paper's own final covariate equation is written as a two-branch expression, 8.45 for male participants and 8.45 * (1 - 0.234) for female participants, which is the same linear indicator form as its generic categorical-covariate Equation 5. Sex was the only demographic covariate retained on CL/F and the paper judges it not clinically relevant.",
      source_name        = "sex"
    ),
    DIS_MM = list(
      description        = "Relapsed / refractory multiple myeloma patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy normal participant)",
      notes              = "1 = patient with relapsed and refractory multiple myeloma (the four MM studies CC-4047-MM-001, -002, -003 and -005); 0 = healthy normal participant (the two healthy-volunteer studies CC-4047-CP-006 and -007). Reference category is the healthy participant, because Table 2 reports each MM effect as a 'MM patient / HNP' ratio on top of a healthy-participant typical value. Drives FIVE things at once: the CL/F ratio 0.913, the V2/F ratio 1.20, the V3/F ratio 8.46, the Q/F ratio 3.71, and a switch between two absolute absorption lag times (0.385 h healthy, 0.206 h MM). It also selects the log-scale residual SD, the paper having fitted a separate residual magnitude for the two populations on the grounds that the healthy-volunteer studies were the better controlled.",
      source_name        = "status of health (healthy participants vs. patients)"
    )
  )

  # Covariates the paper screened in its stepwise covariate search but did NOT
  # retain in the final model. Documentation only -- none is referenced in
  # model(). Listed because the paper's negative findings on renal function are
  # a headline result, not an omission.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; median 53.0 years (range 19.0-83.0, Table 1). Graphic analysis suggested a negative correlation with CL/F and a positive one with V2/F, neither of which reached significance in the forward-selection step."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Listed among the demographic covariates tested (Methods 'Covariate Analysis'); not retained. Table 1 reports height and body mass index but not BSA itself."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a two-level covariate: 78% of the cohort was white, so all non-white participants were pooled into a single group for the race analysis (Results 'Covariate Analysis'). Not retained; the paper reports no apparent relationship between V2/F and race."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; Hispanic or Latino participants (22.5%) were pooled with participants of unknown ethnicity (11.9%) into a single group. Graphic analysis showed a significant CL/F difference between that pooled group and non-Hispanic participants, but ethnicity was not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker; median 40.0 (range 17.0-52.0, Table 1, reported there as g/dL but numerically on the g/L scale). Positively correlated with CL/F and negatively with V2/F by graphic analysis; not significant in forward selection."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker; median 9.2 umol/L (range 1.9-52.3, Table 1). 12% of participants were above the 17 umol/L upper limit of normal and 3% had moderate-to-severe hepatic impairment. Negatively correlated with V2/F by graphic analysis; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker; median 22.0 U/L (range 9.0-73.0, Table 1). Not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault, NOT BSA-normalized",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as the renal-function marker; median 100.4 mL/min (range 20.8-188.2, Table 1), with 17% of participants having CLcr 30-60 mL/min and 1.3% below 30 mL/min. Positively correlated with CL/F by graphic analysis (correlation coefficient 0.2042) but did not reach the prespecified significance threshold in univariate analysis, and the geometric-mean CL/F was comparable across normal, mild and moderate renal-impairment strata (Figure 4A). A post hoc linear regression of CL/F on CLcr gave an intercept of 5.93 L/h (90%CI 4.45-7.41) and a slope of 0.019 (90%CI 0.0016-0.036), i.e. non-renal clearance accounts for about 77% of whole-body clearance. This absence of a renal effect is one of the paper's principal conclusions, so CRCL is documented here rather than silently dropped."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 240L,
    n_studies      = 6L,
    age_range      = "19.0-83.0 years",
    age_median     = "53.0 years",
    weight_range   = "44.4-127.0 kg",
    weight_median  = "78.1 kg",
    sex_female_pct = 25.8,
    race_ethnicity = c(White = 78.0, Black = 19.9, Asian = 0.8, `Native Hawaiian or other Pacific Islander` = 0.4, Other = 0.8),
    disease_state  = "Pooled healthy normal participants (n = 96, studies CC-4047-CP-006 and CP-007) and patients with relapsed and refractory multiple myeloma (n = 144, studies CC-4047-MM-001, -002, -003 and -005). The multiple myeloma cohort spanned disease stages I to III and included the comorbid condition of renal impairment: 40 patients had CLcr 30-60 mL/min and 3 had CLcr below 30 mL/min.",
    dose_range     = "0.5-10 mg oral solid dosage form, once daily or once every alternate day. Healthy participants received 0.5-2 mg QD for 5 days or a single 3 or 4 mg dose; MM patients received 1-10 mg QD for up to 4 weeks, or 4 mg QD on days 1-21 of a 28-day cycle, or 1-4 mg QD on days 1-14 of a 21-day cycle (Supplementary Table 1).",
    regions        = "Not reported.",
    notes          = "3,909 evaluable pomalidomide plasma concentration records. Baseline characteristics are Table 1; the per-study designs and PK sampling schedules are Supplementary Table 1. Plasma was assayed by LC-MS/MS with a lower limit of quantification of 0.25 ng/mL. The final model characterised concentrations spanning 1.0-179 ng/mL. NONMEM 7.2, FOCEI, on log-transformed concentrations; stability was confirmed by a 500-replicate nonparametric bootstrap, 483 (96.6%) of which minimised successfully."
  )

  ini({
    # ---- Structural parameters ----
    # All parameters are APPARENT (CL/F, V2/F, Q/F, V3/F): pomalidomide was
    # given only orally in the pooled dataset, so bioavailability is not
    # separately identifiable and no f(depot) term is carried.
    # Reference subject is a HEALTHY male participant at the covariate
    # equation's centering values, WT 78.3 kg and TPRO 73.0 g/L.

    lka       <- log(1.25)  ; label("First-order absorption rate constant (1/h)")                    # Table 2: ka = 1.25 (bootstrap 1.07-1.47)
    lcl       <- log(8.52)  ; label("Apparent clearance CL/F in healthy participants (L/h)")         # Table 2: CL/F = 8.52 (bootstrap 8.04-8.99). See the model file's Errata note: the final covariate equation prints 8.45, which the paper's own arithmetic falsifies.
    lvc       <- log(58.3)  ; label("Apparent central volume V2/F in healthy participants (L)")      # Table 2: V2/F = 58.3 (bootstrap 55.79-60.88)
    lvp       <- log(8.45)  ; label("Apparent peripheral volume V3/F in healthy participants (L)")   # Table 2: V3/F = 8.45 (bootstrap 7.47-9.40)
    lq        <- log(1.01)  ; label("Apparent intercompartmental clearance Q/F in healthy participants (L/h)")  # Table 2: Q/F = 1.01 (bootstrap 0.75-1.28)

    # Absorption lag time. The paper estimated an absolute value in each
    # population rather than a ratio, so both strata carry an explicit suffix
    # (parameter-names.md 'Stratum-suffixed parameters').
    ltlag_hnp <- log(0.385) ; label("Absorption lag time in healthy participants (h)")               # Table 2: Alag1 = 0.385 h (bootstrap 0.37-0.40)
    ltlag_mm  <- log(0.206) ; label("Absorption lag time in multiple myeloma patients (h)")          # Table 2: 'Alag1 MM patient, h' = 0.206 (bootstrap 0.177-0.231)

    # ---- Disease-state (multiple myeloma vs healthy) ratios ----
    # Each is a MULTIPLICATIVE RATIO of the MM typical value to the healthy
    # typical value, applied as ratio^DIS_MM (DIS_MM is 0/1).
    e_mm_cl   <- 0.913      ; label("Ratio of CL/F in multiple myeloma to healthy participants (unitless)")  # Table 2: 'CL/F MM patient/HNP' = 0.913 -> MM CL/F = 8.52 * 0.913 = 7.78 L/h (Discussion)
    e_mm_vc   <- 1.20       ; label("Ratio of V2/F in multiple myeloma to healthy participants (unitless)")  # Table 2: 'V2/F MM patient/HNP' = 1.20 -> MM V2/F = 58.3 * 1.20 = 69.9 L (Discussion)
    e_mm_vp   <- 8.46       ; label("Ratio of V3/F in multiple myeloma to healthy participants (unitless)")  # Table 2: 'V3/F MM patient/HNP' = 8.46 -> MM V3/F = 8.45 * 8.46 = 71.5 L (Results)
    e_mm_q    <- 3.71       ; label("Ratio of Q/F in multiple myeloma to healthy participants (unitless)")   # Table 2: 'Q/F MM patient/HNP' = 3.71 -> MM Q/F = 1.01 * 3.71 = 3.75 L/h (Results)

    # ---- Covariate effects ----
    e_wt_vc   <- 0.686      ; label("Power exponent on (WT / 78.3) for V2/F (unitless)")             # Table 2: 'Effect of weight on V2F' = 0.686; Results final covariate equation V2/F = 58.3 * (WT/78.3)^0.686 * (1 + 0.00609 * (TPT - 73.0))
    e_tpro_vc <- 0.00609    ; label("Linear coefficient on (TPRO - 73.0) for V2/F (per g/L)")        # Table 2: 'Effect of TPT on V2F' = 0.00609; same Results equation
    e_sexf_cl <- -0.234     ; label("Fractional change in CL/F for female vs male (unitless)")       # Table 2: 'Effect of sex on CL/F' = -0.234; Results: females have 23.4% lower CL/F

    # ---- Interindividual variability ----
    # Log-normal IIV, P_i = P * exp(eta_i) (Equation 1). Table 2 reports the
    # VARIANCES directly (labelled v^2), which the paper's own CV% figures
    # confirm: sqrt(exp(0.168) - 1) = 42.8% matches the reported 42.77% for
    # CL/F and sqrt(exp(0.0352) - 1) = 18.9% matches the reported 18.9% for
    # V2/F. No IIV was reported on V3/F, Q/F or the lag times.
    etalvc + etalcl ~ c(
      0.0352,   # Table 2: 'v2 V2/F' variance = 0.0352
      0.0599,   # Table 2: 'v V2/F : v CL/F' covariance = 0.0599 (correlation 0.78)
      0.168     # Table 2: 'v2 CL/F' variance = 0.168
    )
    etalka ~ 0.976                                                                                   # Table 2: 'v2 Ka' variance = 0.976 (shrinkage 14.8%)

    # ---- Residual variability ----
    # Concentrations were log-transformed and the residual was additive on the
    # log scale, ln(Cij) = ln(Cmij) + eij (Equation 2), which is exactly a
    # log-normal residual in nlmixr2. The paper fitted a separate magnitude for
    # each population; the SDs below are the square roots of the tabulated
    # variances.
    expSdHnp  <- 0.200      ; label("Log-scale residual SD in healthy participants")                  # Table 2: 'sigma2 (HNP)' = 0.04 -> SD = sqrt(0.04) = 0.200 (shrinkage 5.09%)
    expSdMm   <- 0.490      ; label("Log-scale residual SD in multiple myeloma patients")             # Table 2: 'sigma2 (MM patients)' = 0.240 -> SD = sqrt(0.240) = 0.4899 (shrinkage 7.65%)
  })

  model({
    # ---- 1. Individual parameters ----
    # Covariate effects and the multiple-myeloma ratios are multiplicative and
    # commute; DIS_MM is a 0/1 indicator so ratio^DIS_MM selects the healthy
    # value (exponent 0) or the MM value (exponent 1) exactly.
    ka <- exp(lka + etalka)

    cl <- exp(lcl + etalcl) *
      (1 + e_sexf_cl * SEXF) *
      e_mm_cl^DIS_MM

    vc <- exp(lvc + etalvc) *
      (WT / 78.3)^e_wt_vc *
      (1 + e_tpro_vc * (TPRO - 73.0)) *
      e_mm_vc^DIS_MM

    vp <- exp(lvp) * e_mm_vp^DIS_MM
    q  <- exp(lq)  * e_mm_q^DIS_MM

    # Absorption lag: two absolute estimates selected by disease state.
    tlag <- exp(ltlag_hnp) * (1 - DIS_MM) + exp(ltlag_mm) * DIS_MM

    # ---- 2. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ---- 4. Absorption lag ----
    alag(depot) <- tlag

    # ---- 5. Observation and population-specific residual error ----
    # central is in mg and vc in L, so central/vc is mg/L = ug/mL; the factor
    # 1000 converts to the ng/mL the paper reports.
    Cc <- 1000 * central / vc

    # Same idiom as Cirincione_2017_exenatide.R: build the population-specific
    # SD as a model() intermediate and pass it to lnorm().
    expSd <- expSdHnp * (1 - DIS_MM) + expSdMm * DIS_MM
    Cc ~ lnorm(expSd)
  })
}
