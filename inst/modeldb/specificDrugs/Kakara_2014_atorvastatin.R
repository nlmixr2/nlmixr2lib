Kakara_2014_atorvastatin <- function() {
  description <- paste(
    "PD-only indirect-response Imax model for LDL-cholesterol lowering by",
    "atorvastatin (Kakara 2014). One LDL-C compartment with zero-order",
    "synthesis Kin inhibited by Imax * DOSE / (ID50 + DOSE), where DOSE",
    "is the current daily atorvastatin dose (mg/day) supplied as a",
    "time-varying covariate column. An additive 0.109 contribution to the",
    "inhibition fraction is applied when ezetimibe is coadministered",
    "(CONMED_EZE = 1). The LDL-C synthesis-elimination loop is set up at",
    "steady state by enforcing Kin = Baseline * Kout (Kout derived inside",
    "model() as Kin / Baseline). Baseline LDL-C is age-scaled as",
    "152 * (AGE/62)^(-0.240). Imax (0.567), Kin (32.8 mg/dL/day),",
    "Baseline (152 mg/dL), the age power exponent (-0.240), the ezetimibe",
    "INH contribution (0.109), and the IIV magnitudes are shared with",
    "Kakara_2014_pitavastatin and Kakara_2014_rosuvastatin (one joint",
    "NONMEM 7.2 FOCE-INTER fit across 378 patients). Atorvastatin ID50 =",
    "2.22 mg per Kakara 2014 Table 2."
  )
  reference <- paste(
    "Kakara M, Nomura H, Fukae M, Gotanda K, Hirota T, Matsubayashi S,",
    "Shimomura H, Hirakawa M, Ieiri I.",
    "Population pharmacodynamic analysis of LDL-cholesterol lowering",
    "effects by statins and co-medications based on electronic medical",
    "records.",
    "Br J Clin Pharmacol. 2014;78(4):824-835.",
    "doi:10.1111/bcp.12405.",
    sep = " "
  )
  vignette <- "Kakara_2014_statins_LDLC"
  units <- list(
    time          = "day",
    dosing        = "mg/day (DOSE covariate)",
    concentration = "mg/dL (LDL-C)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age in years (time-fixed; baseline value).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline LDL-C is scaled by (AGE/62)^(-0.240), with the reference",
        "age 62 years the median of the study cohort (Kakara 2014 Table 1).",
        "Older subjects therefore have lower predicted baseline LDL-C; the",
        "exponent value -0.240 (PW age for Baseline, Kakara 2014 Table 2)",
        "is shared across atorvastatin, pitavastatin, and rosuvastatin",
        "files."
      ),
      source_name        = "Age"
    ),
    DOSE = list(
      description        = paste(
        "Current administered atorvastatin daily dose (mg/day) carried as",
        "a time-varying data column. Set to 0 before treatment initiation",
        "or during drug holidays so the dose-driven inhibition INH becomes",
        "0 and LDL-C stays at the age-adjusted baseline."
      ),
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "K-PD-style usage (covariate-columns.md DOSE 'use case (b)'):",
        "the daily dose drives the inhibition function INH = Imax * DOSE",
        "/ (ID50 + DOSE) directly without a PK ODE. The Kakara 2014 study",
        "observed atorvastatin daily doses of 5, 10, 15, or 20 mg",
        "(Table 1); the model accepts any non-negative dose. The paper",
        "treats the daily dose as the immediate exposure driver because",
        "atorvastatin plasma concentrations were not measured.",
        "Multiple-day dose-change events are represented by stepping DOSE",
        "as a step function across records."
      ),
      source_name        = "Dose"
    ),
    CONMED_EZE = list(
      description        = paste(
        "Concomitant ezetimibe coadministration indicator (1 = on",
        "ezetimibe, 0 = atorvastatin monotherapy)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ezetimibe coadministration)",
      notes              = paste(
        "In the Kakara 2014 dataset no atorvastatin-treated patient",
        "received ezetimibe coadministration (Table 1: ezetimibe row =",
        "0 for atorvastatin); the indicator is retained here for",
        "structural completeness so simulations of hypothetical",
        "atorvastatin + ezetimibe combinations are supported. The",
        "additive 0.109 contribution to INH (Kakara 2014 Table 2,",
        "INH_EZT row) was estimated on the rosuvastatin + ezetimibe",
        "subgroup (n=12) and applied across all three statins in the",
        "joint model."
      ),
      source_name        = "Ezetimibe"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 149L,
    n_studies        = 1L,
    age_range        = "31-89 years",
    age_median       = "62 years",
    sex_female_pct   = 28.9,
    race_ethnicity   = c(Japanese = 100),
    disease_state    = paste(
      "Adults with hypercholesterolaemia requiring atorvastatin therapy,",
      "enrolled retrospectively from electronic medical records at the",
      "Fukuoka Tokushukai Medical Center, Japan."
    ),
    dose_range       = "5, 10, 15, or 20 mg atorvastatin once daily (Kakara 2014 Table 1)",
    regions          = "Japan",
    n_observations   = 2863L,
    biomarker        = "Low-density lipoprotein cholesterol (LDL-C, mg/dL)",
    notes            = paste(
      "Atorvastatin cohort of the Kakara 2014 retrospective EMR study.",
      "Across the full study 378 patients contributed 2863 LDL-C samples;",
      "the atorvastatin subcohort had n=149 patients (Table 1) with",
      "median dosing period 362 days (range 28-854) and median 7 LDL-C",
      "samples per patient (range 1-25). Baseline LDL-C in the cohort",
      "had median 151 mg/dL (range 71-359) and median total cholesterol",
      "232 mg/dL (Table 1). One atorvastatin patient received concomitant",
      "fenofibrate and one received tocopherol nicotinate; no atorvastatin",
      "patient received bezafibrate or ezetimibe. The total-cohort joint",
      "model in Table 2 was the only model published; per-cohort fits",
      "were not reported."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters from Kakara 2014 Table 2 (Fixed effects).
    # The paper reports the typical population values on the linear
    # scale; the log-transformed forms below match the exponential
    # IIV model the paper used ("The exponential error model was used
    # to describe interindividual and residual variabilities.",
    # Methods, Population pharmacodynamic analysis).
    # ---------------------------------------------------------------
    lbase  <- log(152);   label("Baseline LDL-C at reference age 62 y (mg/dL)")               # Kakara 2014 Table 2 (Baseline = 152, RSE 1.06%)
    lkin   <- log(32.8);  label("LDL-C synthesis rate constant Kin (mg/dL/day)")              # Kakara 2014 Table 2 (Kin = 32.8, RSE 8.81%)
    limax  <- log(0.567); label("Maximum statin-driven inhibition of Kin Imax (fraction)")    # Kakara 2014 Table 2 (Imax = 0.567, RSE 7.72%)
    lid50  <- log(2.22);  label("Atorvastatin daily dose at 50% of Imax ID50,atorva (mg/day)") # Kakara 2014 Table 2 (ID50,Atorvastatin = 2.22, RSE 32.3%)

    # ---------------------------------------------------------------
    # Covariate effects. Both shared across statins in the paper's
    # joint model.
    # ---------------------------------------------------------------
    e_age_base       <- -0.240; label("Power exponent of (AGE/62) on Baseline (unitless)")           # Kakara 2014 Table 2 (PW age for Baseline = -0.240, RSE 29.8%)
    e_conmed_eze_kin <-  0.109; label("Additive contribution of ezetimibe coadministration to INH (fraction)")  # Kakara 2014 Table 2 (INH_EZT = 0.109, RSE 18.9%)

    # ---------------------------------------------------------------
    # IIV -- exponential (log-normal multiplicative). The paper
    # reports CV% in Table 2; convert to log-scale variance via the
    # log-normal identity omega^2 = log(1 + CV^2). The ID50 IIV is
    # a single shared variance across all three statins in the joint
    # fit; the per-statin model files therefore each carry the same
    # 0.26611 value on etalid50.
    # ---------------------------------------------------------------
    etalbase ~ 0.02225  # log(1 + 0.150^2);  Kakara 2014 Table 2 Baseline IIV  CV 15.0%
    etalimax ~ 0.15716  # log(1 + 0.414^2);  Kakara 2014 Table 2 Imax     IIV  CV 41.4%
    etalid50 ~ 0.26611  # log(1 + 0.554^2);  Kakara 2014 Table 2 ID50     IIV  CV 55.4% (shared across atorvastatin, pitavastatin, rosuvastatin)

    # ---------------------------------------------------------------
    # Residual error -- proportional (Kakara 2014 Discussion: "...
    # residual variability was best described with a proportional
    # error model.", and Table 2 footnote: "sigma^2, proportional
    # residual variance"). Table 2 reports sigma^2 (%) = 15.2 as a
    # 15.2 % proportional SD, matching the AitOudhia_2012 convention
    # where the table-percentage value is used directly as propSd.
    # ---------------------------------------------------------------
    propSd <- 0.152; label("Proportional residual error on LDL-C (fraction)")  # Kakara 2014 Table 2 (sigma^2 = 15.2%, RSE 2.94%)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Derived covariate terms.
    #    Baseline LDL-C scales with age via a power function with
    #    reference age 62 years and exponent e_age_base = -0.240
    #    (Kakara 2014 Results: "Baseline (mg dl-1) = 152 * (age/62)^
    #    -0.240").
    # ---------------------------------------------------------------
    age_factor <- (AGE / 62)^e_age_base

    # ---------------------------------------------------------------
    # 2. Individual structural parameters.
    # ---------------------------------------------------------------
    base_i <- exp(lbase + etalbase) * age_factor
    imax_i <- exp(limax + etalimax)
    id50_i <- exp(lid50 + etalid50)
    kin    <- exp(lkin)

    # Kout follows from the steady-state relation Kin = Baseline *
    # Kout (Kakara 2014 Methods: "the relationships between Kin,
    # Kout and baseline are defined by Kin = Baseline * Kout").
    kout   <- kin / base_i

    # ---------------------------------------------------------------
    # 3. Drug effect.
    #    INH = Imax * DOSE / (ID50 + DOSE) + 0.109 * CONMED_EZE
    #    (Kakara 2014 Results, regression equations). The dose-driven
    #    Imax inhibition is the canonical Imax form on the daily
    #    dose; the ezetimibe contribution is additive and non-zero
    #    only when ezetimibe is coadministered.
    # ---------------------------------------------------------------
    inh <- imax_i * DOSE / (id50_i + DOSE) + e_conmed_eze_kin * CONMED_EZE

    # ---------------------------------------------------------------
    # 4. ODE -- single LDL-C compartment, zero-order synthesis Kin
    #    inhibited by INH and first-order elimination at Kout
    #    (Kakara 2014 Methods, indirect response model). Canonical
    #    compartment `ldl` per inst/references/compartment-names.md.
    # ---------------------------------------------------------------
    d/dt(ldl) <- kin * (1 - inh) - kout * ldl

    # Initial condition -- steady-state pre-treatment LDL-C equals
    # the (age-adjusted) Baseline.
    ldl(0) <- base_i

    # ---------------------------------------------------------------
    # 5. Observation and residual error.
    # ---------------------------------------------------------------
    Cc <- ldl
    Cc ~ prop(propSd)
  })
}
