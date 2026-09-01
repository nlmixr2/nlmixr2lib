Zhao_2026_venetoclax <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for oral venetoclax in ",
    "Chinese pediatric patients with hematological malignancy (Zhao 2026): ",
    "first-order absorption with ka fixed at 0.15 1/h; a power effect of body ",
    "surface area (exponent 1.4, reference 1.27 m^2) and an exponential effect ",
    "of concomitant triazole antifungal therapy (coefficient -1.1, a 67% ",
    "reduction) on CL/F; and a power effect of total serum protein (exponent ",
    "1.7, reference 62 g/L) on V/F. Developed from sparse real-world ",
    "therapeutic drug monitoring data at a single Beijing centre."
  )
  reference <- paste0(
    "Zhao Y, Song X, Zhang L, Zhu Y, Chen J, Gong Y, Luo X, He H, Zhang X, ",
    "Huang L. Population Pharmacokinetic Analyses and Exposure-Efficacy ",
    "Relationships of Venetoclax in Chinese Pediatric Patients with ",
    "Hematological Malignancy in a Real-World Setting. Drug Des Devel Ther. ",
    "2026;20:583847. doi:10.2147/DDDT.S583847."
  )
  vignette <- "Zhao_2026_venetoclax"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhao 2026 Methods (oral venetoclax,
  # plasma concentrations by UHPLC with carbamazepine internal standard).
  compartmentData <- list(
    depot   = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form effect on CL/F with exponent 1.4, mean-centred at 1.27 m^2 ",
        "(Zhao 2026 Eq. 6 and Table 3; 1.27 m^2 is the cohort mean BSA reported ",
        "in Table 1, so the centring constant equals the cohort mean and is a ",
        "free transcription check). BSA was computed by the Stevenson formula, ",
        "BSA = 0.0061 * height(cm) + 0.0128 * weight(kg) - 0.1529 ",
        "(Discussion); use that formula rather than DuBois or Mosteller when ",
        "reproducing the covariate. BSA was preferred over height and weight ",
        "entered separately because pediatric venetoclax is dosed on BSA ",
        "(120 mg/m^2 on day 1 and 240 mg/m^2 on days 2-28 when not ",
        "co-administered with azoles). Baseline value in the source dataset."
      ),
      source_name        = "BSA"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form effect on V/F with exponent 1.7, mean-centred at 62 g/L ",
        "(Zhao 2026 Eq. 5 and Table 3; 62 g/L is the cohort median and is ",
        "within rounding of the cohort mean 62.53 g/L reported in Table 1). ",
        "The source reports total protein in SI units (g/L) already, so no ",
        "conversion is needed; a US-convention g/dL value must be multiplied ",
        "by 10 before use. The positive exponent means lower total protein ",
        "gives a smaller apparent volume, which the authors attribute to ",
        "reduced plasma-protein binding raising the measured concentration ",
        "(Discussion). Venetoclax is >99% plasma-protein bound. Baseline value ",
        "in the source dataset."
      ),
      source_name        = "TP"
    ),
    CONMED_AZOLE = list(
      description        = "Concomitant triazole antifungal therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant triazole antifungal)",
      notes              = paste0(
        "Exponential effect on CL/F with coefficient -1.1 (Zhao 2026 Eq. 6 and ",
        "Table 3), i.e. CL/F is multiplied by exp(-1.1) = 0.333, a 67% ",
        "reduction. Pooled over posaconazole and voriconazole only -- both ",
        "strong CYP3A inhibitors and P-glycoprotein inhibitors; no ",
        "fluconazole, itraconazole or isavuconazole patients were in the ",
        "cohort (Discussion). 32 of 96 patients (33.3%) received a triazole ",
        "(Table 1). The source treats the indicator as a per-subject baseline ",
        "flag rather than an explicitly time-varying one; the paper does not ",
        "describe a post-cessation lag, so unlike Kirubakaran 2022 tacrolimus ",
        "no carry-forward window should be assumed."
      ),
      source_name        = "Triazole"
    )
  )

  # Screened in the stepwise covariate search (Methods, Covariate Model) but not
  # retained in the final model, so they are documentation only and are not
  # referenced in model(). Restricted here to concepts that already have a
  # canonical register entry; the source additionally screened platelet count,
  # lymphocyte percentage, eGFR (CKD-EPI), total carbon dioxide and the
  # malignancy-type categorical, none of which has a canonical column name.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a mean-centred power covariate; not retained. Zhao 2026 chose not to model CYP3A maturation explicitly (unlike Badawi 2022 pediatric venetoclax), citing Phoenix NLME limitations (Discussion)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a mean-centred power covariate; not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a mean-centred power covariate; not retained. Most enrolled patients had normal liver function, limiting the power to detect a hepatic effect (Limitations)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a mean-centred power covariate; not retained. Most enrolled patients had normal renal function (Limitations)."
    ),
    LDH = list(
      description = "Serum lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a mean-centred power covariate; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 1L,
    n_concentrations = 225L,
    age_range      = "0.3-18 years",
    age_median     = "12 years",
    age_mean       = "11 years (SD 4.6)",
    weight_range   = "8-112 kg",
    weight_median  = "45 kg",
    weight_mean    = "42 kg (SD 21.2)",
    height_range   = "63-190 cm",
    height_median  = "156 cm",
    bsa_range      = "0.33-2.35 m^2",
    bsa_median     = "1.38 m^2",
    bsa_mean       = "1.27 m^2 (SD 0.4; the model centring constant)",
    sex_female_pct = 51.0,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste0(
      "Pediatric hematological malignancy: acute myeloid leukemia 64/96 ",
      "(66.7%), acute lymphoblastic leukemia 20/96 (20.8%), non-Hodgkin ",
      "lymphoma 5/96 (5.2%), mixed lineage acute leukemia 4/96 (4.2%), ",
      "chronic myeloid leukemia 2/96 (2.1%), myelodysplastic syndrome 1/96 ",
      "(1.0%) (Table 1). Malignancy type was screened as a categorical ",
      "covariate but not retained."
    ),
    dose_range     = paste0(
      "Oral venetoclax once daily. Stepwise dose escalation over 2-3 days to a ",
      "maintenance dose of 100-400 mg once daily, adjusted for response, ",
      "adverse reactions and tolerance (Study Design). Clinical practice at ",
      "the centre doses on BSA: 120 mg/m^2 on day 1 and 240 mg/m^2 on days ",
      "2-28 when not co-administered with azole antifungals (Discussion)."
    ),
    regions        = "Single centre, Peking University People's Hospital, Beijing, China",
    conmed_triazole_pct = 33.3,
    total_protein_range = "39.2-99.4 g/L (median 62, mean 62.53, SD 7.50)",
    covariates_screened_not_retained = paste0(
      "The stepwise search screened age, BSA, lymphocyte percentage (LY%), ",
      "hemoglobin, platelet count, total protein, total bilirubin, ",
      "creatinine, eGFR (CKD-EPI), total carbon dioxide (TCO2) and lactate ",
      "dehydrogenase as mean-centred power covariates, plus malignancy type ",
      "and triazole use as categoricals (Methods, Covariate Model). Only ",
      "BSA and triazole on CL/F and total protein on V/F were retained. The ",
      "screened-but-unretained concepts that have a canonical covariate ",
      "column are listed in covariatesDataExcluded; lymphocyte percentage, ",
      "platelet count, eGFR, total carbon dioxide and the malignancy-type ",
      "categorical have no canonical column in ",
      "inst/references/covariate-columns.md and are recorded here only, ",
      "since none is referenced by model()."
    ),
    notes          = paste0(
      "Retrospective real-world therapeutic drug monitoring, December 2021 to ",
      "December 2024. Sparse sampling: a pre-dose trough drawn 30 min before ",
      "the next dose (22-24 h post-dose, C0) and a 6 h post-dose sample ",
      "(6 +/- 0.5 h, C6), with additional samples at 2, 4, 8 and 10 h; ",
      "outpatient samples were taken after at least 5 days of continuous ",
      "dosing at the prescribed dose. 16 observations were excluded for being ",
      "below the limit of quantification or lacking dosing information, ",
      "leaving 225. Bioanalysis by UHPLC with a Thermo Hypersil GOLD C18 ",
      "column and carbamazepine internal standard; LLOQ 0.1 ug/mL, calibration ",
      "range 0.1-10 ug/mL. Estimation in Phoenix NLME 8.3 by FOCE-ELS. ",
      "External evaluation on 84 trough concentrations from 39 further ",
      "patients (April-December 2024) gave a mean prediction error of 17.7% ",
      "and a mean absolute error of 28.0%. A separate exposure-efficacy ",
      "analysis in 52 pediatric AML patients receiving venetoclax with a ",
      "hypomethylating agent is reported in the same paper but is a ",
      "single-factor logistic regression fitted in SPSS, not part of this ",
      "pharmacokinetic model."
    )
  )

  ini({
    # ----- Structural PK (Zhao 2026 Table 3 final-model column, Eq. 5 and 6) -----
    # One-compartment with first-order oral absorption. The authors also fitted a
    # two-compartment model and found no meaningful difference in -2LL or AIC,
    # and retained the one-compartment structure (Discussion).
    # ka fixed at the literature value because the sparse TDM sampling scheme
    # carries few absorption-phase samples (Structural Model; Limitations).
    lka <- fixed(log(0.15)); label("Absorption rate constant (ka, 1/h)")  # Zhao 2026 Table 3, ka = 0.15 1/h, taken from published literature
    lcl <- log(4.8);         label("Apparent clearance at reference covariates (CL/F, L/h)")  # Zhao 2026 Table 3, CL/F = 4.8 L/h (RSE 26.5%, 95% CI 2.3-7.3); Eq. 6
    lvc <- log(124.7);       label("Apparent volume of distribution at reference total protein (V/F, L)")  # Zhao 2026 Table 3, V/F = 124.7 L (RSE 15.0%, 95% CI 87.8-161.7); Eq. 5

    # ----- Covariate effects (Zhao 2026 Eq. 5 and Eq. 6; Table 3) -----
    # Continuous covariates enter as a mean-centred power (Eq. 1,
    # Pi = P * (Cov / Mean_cov)^theta_cov); the binary triazole indicator enters
    # exponentially (Eq. 2, Pi = P * exp(theta_cov * Cov)).
    e_bsa_cl          <-  1.4;  label("Power exponent of body surface area on CL/F (unitless)")  # Zhao 2026 Table 3, BSA on CL/F = 1.4 (RSE 18.4%, 95% CI 0.9-2.0); Eq. 6
    e_conmed_azole_cl <- -1.1;  label("Exponential coefficient of concomitant triazole antifungal on CL/F (unitless)")  # Zhao 2026 Table 3, Triazole on CL/F = -1.1 (RSE -34.4%, 95% CI -1.9 to -0.4); Eq. 6
    e_tpro_vc         <-  1.7;  label("Power exponent of total serum protein on V/F (unitless)")  # Zhao 2026 Table 3, TP on V/F = 1.7 (RSE 24.4%, 95% CI 0.9-2.5); Eq. 5

    # ----- IIV (Zhao 2026 Table 3 final-model column) -----
    # Lognormal IIV on CL/F and V/F (Structural Model). Diagonal: the Methods
    # mention "interindividual variability and covariance" but Table 3 reports
    # only the two diagonal terms and no covariance, so no block is encoded.
    #
    # SCALE INFERENCE (not a printed statement). Table 3's "Interindividual
    # variability" rows print omega_CL/F = 0.996 and omega_V/F = 0.608 with no
    # scale label, and the table's own 95% CI column is internally inconsistent
    # (see below), so the SD-vs-variance reading had to be inferred. Two
    # independent tests agree that the printed values are log-scale SDs, hence
    # the variances below are the printed values SQUARED:
    #
    #   (a) Base-vs-final variance decomposition. The Results give IIV on CL/F
    #       as 1.168 in the base model and 0.996 in the final model, so the two
    #       retained CL/F covariates must account for the drop. From Table 1 and
    #       Eq. 6 those covariates inject
    #         1.4^2 * var(log BSA) + 1.1^2 * p(1-p)
    #         = 1.96 * log(1 + (0.40/1.27)^2) + 1.21 * (1/3)(2/3) = 0.454.
    #       Reading the printed numbers as SDs gives a drop of
    #       1.168^2 - 0.996^2 = 0.372 (ratio 0.82); reading them as variances
    #       gives 1.168 - 0.996 = 0.172 (ratio 0.38). The SD reading lands in
    #       the expected band, the variance reading is a factor of ~2.7 short.
    #   (b) The V/F confidence interval. Table 3's IIV 95% CI cells are
    #       transposed between the two rows (0.996 is printed with 0.2-0.5 and
    #       0.608 with 0.6-1.3, yet the bootstrap column pairs them the other
    #       way round). Reassigning them, the V/F interval 0.2-0.5 reproduces
    #       the printed 22.6% RSE only on the VARIANCE 0.608^2 = 0.370
    #       (implied RSE 20.7%), not on 0.608 itself (implied RSE 12.6%) --
    #       i.e. the estimate column and the CI column are on different scales,
    #       and the estimate column is the SD.
    #
    # The reading is immaterial for CL/F (0.996 as a variance is 131% CV and as
    # an SD is 130% CV) but material for V/F (67% CV as an SD vs 91% CV as a
    # variance). See the vignette Errata for the full reconciliation.
    etalcl ~ 0.992016  # 0.996^2; log-scale variance, 130% CV; Zhao 2026 Table 3 omega CL/F = 0.996 (RSE 17.6%); base model 1.168
    etalvc ~ 0.369664  # 0.608^2; log-scale variance, 67% CV; Zhao 2026 Table 3 omega V/F = 0.608 (RSE 22.6%)

    # ----- Residual error -----
    # "A proportional model described residual error" / "a multiplicative model
    # for residual error" (Structural Model; Population Pharmacokinetic
    # Modeling). sigma = 0.38 is on the SD scale, by two tests. (a) Its own
    # 95% CI: 0.38 +/- 1.96 * 0.085 * 0.38 = 0.32-0.44 reproduces the printed
    # 0.3-0.4, whereas the same RSE on the variance 0.38^2 = 0.144 would give
    # 0.12-0.17. (Note this makes the sigma row's CI SD-scale while the two IIV
    # rows' CIs are variance-scale -- the table mixes scales; each row is
    # nevertheless self-consistent, and both readings agree the ESTIMATE column
    # is SD-scale.) (b) An RSE of 8.5% is below the Cramer-Rao floor for a
    # variance over 225 observations, sqrt(2/225) = 9.4%, but above the SD
    # floor sqrt(1/450) = 4.7%.
    propSd <- 0.38;  label("Proportional residual error (fraction)")  # Zhao 2026 Table 3, sigma = 0.38 (RSE 8.5%, 95% CI 0.3-0.4)
  })

  model({
    # ----- 1. Derived covariate terms -----
    # Eq. 1 (mean-centred power) for the continuous covariates and Eq. 2
    # (exponential) for the binary triazole indicator. Centring constants are
    # the cohort means quoted in the Results: BSA 1.27 m^2, TP 62 g/L.
    bsa_cl   <- (BSA / 1.27)^e_bsa_cl
    azole_cl <- exp(e_conmed_azole_cl * CONMED_AZOLE)
    tpro_vc  <- (TPRO / 62)^e_tpro_vc

    # ----- 2. Individual PK parameters -----
    # Eq. 6: CL/F = 4.8 * (BSA/1.27)^1.4 * exp(-1.1 * Triazole)
    # Eq. 5: V/F  = 124.7 * (TP/62)^1.7
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * bsa_cl * azole_cl
    vc <- exp(lvc + etalvc) * tpro_vc

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system (1-cmt with first-order oral absorption) -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Bioavailability -----
    # CL and V are apparent (CL/F, V/F); F is absorbed into them and is not
    # separately identifiable from oral-only data. rxode2's default F = 1 is
    # correct here, so no f(depot) statement is needed.

    # ----- 6. Observation and error -----
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # report venetoclax plasma concentration in ng/mL, the unit used for the
    # C0 / C6 exposure metrics in Table 2 and the exposure-efficacy analysis.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
