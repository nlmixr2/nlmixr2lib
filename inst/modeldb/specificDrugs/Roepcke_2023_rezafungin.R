Roepcke_2023_rezafungin <- function() {
  description <- "Three-compartment population PK model for rezafungin after weekly IV infusion in healthy subjects, hepatically impaired subjects, and patients with candidemia and/or invasive candidiasis (Roepcke 2023), with body-surface-area scaling on CL, V1, and the shared peripheral volume V23, a serum-albumin effect on V23, and a healthy-vs-diseased disease-state shift on CL and V1."
  reference   <- "Roepcke S, Passarell J, Walker H, Flanagan S. Population pharmacokinetic modeling and target attainment analyses of rezafungin for the treatment of candidemia and invasive candidiasis. Antimicrob Agents Chemother. 2023;67(12):e00916-23. doi:10.1128/aac.00916-23"
  vignette    <- "Roepcke_2023_rezafungin"
  units       <- list(time = "hr", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, V1, and the shared peripheral volume V23 with reference 1.9 m^2 (the overall median BSA reported in Roepcke 2023 Table 1). Exponents 0.882 on CL, 1.56 on V1, and 1.17 on V23 per Table 2 and its footnote a. The BSA computation formula (DuBois / Mosteller / Haycock) is unspecified in the source; the paper reports BSA only as a tabulated baseline characteristic. Body height was missing for 7 patients and body weight for 5 patients; these were imputed with study- and sex-specific median values before deriving BSA (Roepcke 2023 Results, 'Data used for PK modeling'). The paper notes that BMI and body weight were highly correlated with BSA and that model fit would be expected to be comparable with any of these body-size measures.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Serum albumin at baseline",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on the shared peripheral volume V23 only, exponent -0.708, with reference 3.2 g/dL (Roepcke 2023 Table 2 footnote a: V23 = 19.1 x (BSA/1.9)^1.17 x (ALB/3.2)^-0.708). Table 2 lists the exponent magnitude 0.708 without a sign; the printed equation in footnote a carries the explicit minus sign, and the negative orientation is confirmed by Table 3, where the low-albumin patient cohort (median 2.7 g/dL) has a larger total steady-state volume of distribution (62.4 L) than healthy subjects (median 4.5 g/dL, 42.3 L). Roepcke 2023 reports albumin in US-convention g/dL; the canonical register unit is SI g/L, so model() applies an inline conversion alb_gdL <- ALB * 0.1 before the power term to keep the structural exponent aligned with its original calibration. Reference 3.2 g/dL = 32 g/L. Albumin is lower in patients with candidemia / invasive candidiasis and in hepatically impaired subjects, so this covariate is partly collinear with DIS_HEALTHY.",
      source_name        = "ALB"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (disease state)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (diseased: patient with candidemia and/or invasive candidiasis, or subject with impaired hepatic function from Study CD101.IV.1.15)",
      notes              = "Roepcke 2023 Table 2 footnote a defines I_healthy as 'the indicator variable for healthy subjects who are not infected and have normal liver function (0 = no; 1 = yes)', which maps directly onto the canonical DIS_HEALTHY orientation with no re-expression. Linear fractional effects: CL x (1 - 0.276 x DIS_HEALTHY) and V1 x (1 - 0.222 x DIS_HEALTHY). The typical values 0.328 L/h and 17.7 L are therefore the diseased-state reference. Results, 'Clinical relevance' confirms the direction: 'healthy subjects had a 27.6% lower CL than infected subjects or those with hepatic impairment'. The composite category was constructed during covariate analysis by pooling (a) patients from the Phase 2 STRIVE and Phase 3 ReSTORE studies with (b) hepatically impaired subjects from Study CD101.IV.1.15, while the normal-hepatic-function matched controls from that same study were grouped with the healthy subjects; this replaced separate 'study on CL' and 'infection status on V1' effects with a single interpretable dichotomy (Roepcke 2023 Results, 'Covariate analysis')."
    )
  )

  # Covariates screened during the stepwise covariate analysis but NOT retained
  # in the final model. Documented for provenance only; these names are not
  # referenced in model(). Roepcke 2023 Materials and Methods, 'Population
  # pharmacokinetic analysis methodology' lists the full screened set, and
  # Results / Fig. 3 report geometric-mean-ratio forest plots quantifying the
  # (non-clinically-meaningful) exposure differences across several of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained. Roepcke 2023 Results: GMR (90% CI) of AUC0-168h was 0.95 (0.85, 1.05) for ages 65 to <75 years and 1.20 (1.06, 1.36) for ages 75 to 89 years, both relative to patients aged 24 to <65 years; not deemed clinically meaningful."
    ),
    WT = list(
      description = "Body weight at baseline",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained as an independent covariate; highly correlated with the retained BSA term."
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained; highly correlated with BSA. Roepcke 2023 Results reports GMR (90% CI) of AUC0-168h of 0.84 (0.75, 0.94) for overweight and 0.77 (0.69, 0.86) for obese patients relative to the optimum-BMI group; not deemed clinically meaningful."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Female subjects showed a 21% higher AUC0-168h, which the paper attributes to their lower average BSA rather than to an independent sex effect."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened but not retained. Roepcke 2023 reports CrCL in uncorrected mL/min (Cockcroft-Gault), not the canonical body-surface-area-normalised mL/min/1.73 m^2. GMR (90% CI) of AUC0-168h relative to normal renal function was 0.95 (0.85, 1.06), 1.23 (1.11, 1.36), 1.01 (0.88, 1.16), and 1.16 (0.96, 1.39) for mild, moderate, and severe renal impairment and kidney failure; not clinically meaningful."
    ),
    CREAT = list(
      description = "Serum creatinine at baseline",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained. Missing for one patient and imputed with the study- and sex-specific median."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation disease score",
      units       = "points",
      type        = "continuous",
      notes       = "Screened but not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 277L,
    n_studies      = 7L,
    n_healthy      = 94L,
    n_hepimp       = 16L,
    n_patients     = 167L,
    n_observations = 2512L,
    age_range      = "20-89 years",
    age_median     = "53.0 years",
    weight_range   = "34.0-154.5 kg",
    weight_median  = "74.70 kg",
    bmi_range      = "13.7-64.4 kg/m^2",
    bmi_median     = "26.39 kg/m^2",
    bsa_range      = "1.2-2.7 m^2",
    bsa_median     = "1.90 m^2",
    albumin_range  = "1.2-5.1 g/dL",
    albumin_median = "3.20 g/dL",
    sex_female_pct = 39.4,
    race_ethnicity = c(White = 76.5, Black = 9.7, Asian = 10.1, `American Indian or Alaska Native` = 0.7, Other = 0.7, Missing = 2.2),
    disease_state  = "Pooled cohort of 94 healthy subjects, 16 uninfected subjects with moderate or severe hepatic impairment, and 167 patients with candidemia and/or invasive candidiasis (in part critically ill). 16 patients (9.6%) required dialysis; 91 of 167 patients (54%) had renal impairment ranging from mild to kidney failure.",
    renal_function = "Among the 167 patients: 76 normal (CrCL >= 90 mL/min), 33 mild (60-89), 36 moderate (30-59), 14 severe (15-29), and 8 with kidney failure (< 15 mL/min)",
    dose_range     = "50-1,400 mg rezafungin as IV infusion; single and multiple (weekly) dosing. Phase 1 single doses 50/100/200/400 mg over 1 h, 600 mg over 1.5 h, and a divided 1,400 mg dose over 1.5 h + 2 h; Phase 1 multiple doses 100/200/400 mg weekly over 1 h; Phase 2/3 400 mg loading dose week 1 followed by 200 mg (or 400 mg in STRIVE Part A Group 1) weekly over 1 h",
    regions        = "North America, Europe, and Asia-Pacific (the Asia-Pacific sites were unique to Study CD101.IV.3.05, which enrolled 26 of the 28 Asian subjects)",
    studies        = "CD101.IV.1.01, CD101.IV.1.02, CD101.IV.1.06, CD101.IV.1.07, CD101.IV.1.15 (Phase 1); CD101.IV.2.03 / STRIVE / NCT02734862 (Phase 2); CD101.IV.3.05 / ReSTORE / NCT03667690 (Phase 3)",
    notes          = "Baseline demographics from Roepcke 2023 Table 1; study designs and dosing regimens from Table S1 of the supplemental material (aac.00916-23-s0001.docx). Dense PK sampling in the Phase 1 studies and a mixture of dense and sparse sampling in Phase 2 and Phase 3. Of 2,512 concentration records, 5 were below and 5 above the quantitation limits (< 0.5% of observations) and were excluded, along with 8 further records flagged as outliers by |CWRES| > 5 during model development. Model estimated with NONMEM 7.3 using FOCE with interaction."
  )

  ini({
    # Structural parameters at the diseased-state reference
    # (DIS_HEALTHY = 0), BSA = 1.9 m^2 and albumin = 3.2 g/dL.
    # Roepcke 2023 Table 2, 'Population mean' column.
    lcl  <- log(0.328); label("Clearance at the diseased reference, BSA = 1.9 m^2 (CL, L/h)")                                # Roepcke 2023 Table 2 (2.77 %RSE)
    lvc  <- log(17.7);  label("Central volume at the diseased reference, BSA = 1.9 m^2 (V1, L)")                             # Roepcke 2023 Table 2 (3.96 %RSE)
    lvp  <- log(19.1);  label("Shared peripheral volume at BSA = 1.9 m^2 and ALB = 3.2 g/dL (V23, L)")                       # Roepcke 2023 Table 2 (2.41 %RSE)
    lq   <- log(0.236); label("Intercompartmental clearance to peripheral1 (Q2, L/h)")                                       # Roepcke 2023 Table 2 (5.38 %RSE)
    lq2  <- log(12.4);  label("Intercompartmental clearance to peripheral2 (Q3, L/h)")                                       # Roepcke 2023 Table 2 (4.37 %RSE)

    # Body-surface-area power exponents, reference 1.9 m^2
    # (Roepcke 2023 Table 2 footnote a).
    e_bsa_cl <- 0.882; label("Exponent of (BSA / 1.9 m^2) on CL (unitless)")                                                 # Roepcke 2023 Table 2 (16.4 %RSE)
    e_bsa_vc <- 1.56;  label("Exponent of (BSA / 1.9 m^2) on V1 (unitless)")                                                 # Roepcke 2023 Table 2 (9.47 %RSE)
    e_bsa_vp <- 1.17;  label("Exponent of (BSA / 1.9 m^2) on V23 (unitless)")                                                # Roepcke 2023 Table 2 (14.9 %RSE)

    # Serum-albumin power exponent on the shared peripheral volume,
    # reference 3.2 g/dL. Table 2 prints the magnitude 0.708; the sign is
    # taken from the equation in Table 2 footnote a, which reads
    # V23 = 19.1 x (BSA/1.9)^1.17 x (ALB/3.2)^-0.708.
    e_alb_vp <- -0.708; label("Exponent of (ALB / 3.2 g/dL) on V23 (unitless)")                                              # Roepcke 2023 Table 2 footnote a (11.4 %RSE)

    # Disease-state (healthy-vs-diseased) linear fractional shifts.
    # Table 2 lists the magnitudes 0.276 and 0.222 as 'proportional shift
    # ... for healthy'; footnote a prints them with an explicit minus sign
    # inside the (1 + theta x I_healthy) form used in model() below.
    e_dis_healthy_cl <- -0.276; label("Proportional shift in CL for healthy subjects (fraction)")                            # Roepcke 2023 Table 2 footnote a (9.25 %RSE)
    e_dis_healthy_vc <- -0.222; label("Proportional shift in V1 for healthy subjects (fraction)")                            # Roepcke 2023 Table 2 footnote a (18.1 %RSE)

    # Correlated interindividual variability on CL, V1, and V23
    # (Roepcke 2023 Table 2). The paper reports the diagonal as apparent
    # %CV of an exponential (log-normal) IIV model, so the log-scale
    # variance is omega^2 = log(CV^2 + 1):
    #   CL  30.5 %CV -> log(1 + 0.305^2) = 0.088949
    #   V1  37.6 %CV -> log(1 + 0.376^2) = 0.132235
    #   V23 29.3 %CV -> log(1 + 0.293^2) = 0.082362
    # The off-diagonal covariances are taken verbatim from Table 2. This
    # transform is confirmed by Table 2 footnotes b-d: back-computing
    # r = cov / sqrt(var1 * var2) from these variances reproduces the
    # published correlations 0.516, 0.723, and 0.357 exactly, which the
    # naive omega^2 = CV^2 convention does not.
    etalcl + etalvc + etalvp ~ c(0.088949,
                                 0.0560, 0.132235,
                                 0.0619, 0.0373, 0.082362)                                                                  # Roepcke 2023 Table 2 (13.7 / 14.7 / 18.9 %RSE on the diagonals; 17.3 / 18.8 / 38.8 %RSE on the covariances)

    # Residual variability -- proportional / constant coefficient of
    # variation. Table 2 reports the variance 0.00949 and the equivalent
    # 9.74 %CV; sqrt(0.00949) = 0.09742.
    propSd <- 0.0974; label("Proportional residual error (fraction)")                                                        # Roepcke 2023 Table 2 (9.34 %RSE)
  })

  model({
    # 1. Derived covariate terms.
    # The canonical ALB column is SI g/L; Roepcke 2023 calibrated the V23
    # albumin exponent against US-convention g/dL with a 3.2 g/dL
    # reference, so convert before applying the power term.
    alb_gdL <- ALB * 0.1

    # 2. Individual parameters, following Roepcke 2023 Table 2 footnote a
    #    verbatim. Typical values are the diseased-state reference, so the
    #    (1 + e_dis_healthy_* x DIS_HEALTHY) factors shift CL and V1 down
    #    by 27.6% and 22.2% for healthy subjects. Q2 and Q3 carry neither
    #    covariate effects nor IIV (Table 2 lists them as 'NE').
    cl  <- exp(lcl + etalcl) * (BSA / 1.9)^e_bsa_cl * (1 + e_dis_healthy_cl * DIS_HEALTHY)
    vc  <- exp(lvc + etalvc) * (BSA / 1.9)^e_bsa_vc * (1 + e_dis_healthy_vc * DIS_HEALTHY)
    vp  <- exp(lvp + etalvp) * (BSA / 1.9)^e_bsa_vp * (alb_gdL / 3.2)^e_alb_vp
    q   <- exp(lq)
    q2  <- exp(lq2)

    # Roepcke 2023 collapsed the two peripheral volumes V2 and V3 into a
    # single estimated parameter V23 after finding that their typical
    # estimates were very similar and their individual IIVs were highly
    # correlated (Results, 'PK model results'). Both peripheral
    # compartments therefore share one volume and one eta; only the
    # intercompartmental clearances Q2 and Q3 distinguish them.
    vp2 <- vp

    # 3. Micro-constants. Q2 (0.236 L/h) feeds the slowly equilibrating
    #    peripheral1; Q3 (12.4 L/h) feeds the rapidly equilibrating
    #    peripheral2.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 4. Three-compartment IV disposition with first-order elimination.
    #    Doses enter the central compartment as a zero-order infusion via
    #    the event table (1 h for the 50-400 mg regimens; Table S1).
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 5. Observation. Dose in mg over volume in L gives mg/L = ug/mL,
    #    the unit rezafungin concentrations were reported in
    #    (Materials and Methods, 'Study design').
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
