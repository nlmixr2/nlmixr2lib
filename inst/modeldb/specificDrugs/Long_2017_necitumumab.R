Long_2017_necitumumab <- function() {
  description <- "Two-compartment population PK model for necitumumab in cancer patients (Long 2017), with IV infusion input and parallel linear plus Michaelis-Menten (target-mediated) elimination from the central compartment and allometric weight scaling on CL, Q, V1, and V2."
  reference <- "Long A, Chigutsa E, Wallin J. Population Pharmacokinetics of Necitumumab in Cancer Patients. Clin Pharmacokinet. 2017;56(5):505-514. doi:10.1007/s40262-016-0452-x"
  vignette <- "Long_2017_necitumumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying in the source analysis; last-observation-carried-forward imputation for missing)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL and Q (exponent 0.768) and on V1 and V2 (exponent 0.498), each normalized as (WT/70)^exponent per Long 2017 Table 4 footnotes c-f.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 807L,
    n_observations = 4920L,
    n_studies      = 5L,
    age_range      = "19-84 years",
    age_median     = "62 years",
    weight_range   = "35-181 kg",
    weight_median  = "71 kg",
    sex_female_pct = 25,
    race_ethnicity = c(Caucasian = 85, Asian = 7, Black = 2, Other = 5, AmericanIndian = 0, Multiple = 0),
    disease_state  = "Advanced solid tumors (squamous and non-squamous non-small-cell lung cancer, colorectal cancer, and other advanced solid malignancies).",
    dose_range     = "600-800 mg necitumumab IV as approximately 1-h infusion; administered weekly (5%), every 2 weeks (1%), or on days 1 and 8 of a 21-day cycle (94%). Phase III marketed regimen is 800 mg on days 1 and 8 of a 21-day cycle, given with or without concomitant gemcitabine / cisplatin chemotherapy.",
    regions        = "Multi-regional (five pooled clinical studies JFCA/B/C/I/J, three phase I/II with rich PK sampling and two phase III with trough-only sampling).",
    notes          = "Baseline demographics from Long 2017 Tables 2 and 3 (807 patients; 4920 quantifiable concentrations; 184 below-LOQ samples treated as missing). 75% male; ECOG 0/1/2: 33%/60%/6% (2% missing). Concomitant cisplatin 89%; concomitant gemcitabine 42%. Hispanic/Latino ethnicity 11%. Final covariate search found patient bodyweight as the only retained covariate; age, sex, race, ethnicity, liver function (AST, ALT, bilirubin), creatinine clearance, ECOG, and concomitant gemcitabine/cisplatin were tested and not retained."
  )

  ini({
    # Structural PK parameters - Long 2017 Table 4 final-model estimates.
    # Reference body weight = 70 kg (matches the covariate equations in Table 4
    # footnotes c-f: CLind = CL * (WT/70)^0.768, V1ind = V1 * (WT/70)^0.498,
    # etc.). The population median weight is 71 kg (Table 2).
    lcl   <- log(0.0114); label("Linear (non-target-mediated) clearance CL (L/h)")              # Long 2017 Table 4, CL row (final model)
    lvc   <- log(3.41);   label("Central volume of distribution V1 (L)")                         # Long 2017 Table 4, V1 row (final model)
    lq    <- log(0.0183); label("Intercompartmental clearance Q (L/h)")                          # Long 2017 Table 4, Q row (final model)
    lvp   <- log(3.29);   label("Peripheral volume of distribution V2 (L)")                      # Long 2017 Table 4, V2 row (final model)

    # Parallel Michaelis-Menten (target-mediated) elimination from central.
    # CLtot = CL + Vmax / (C + Km); saturation threshold Km is well below
    # the clinical trough range (~40 mg/L), so PK is approximately linear
    # in the studied dose range.
    lvmax <- log(0.565);  label("Maximum target-mediated elimination rate Vmax (mg/h)")          # Long 2017 Table 4, Vmax row (final model)
    lkm   <- log(7.97);   label("Michaelis-Menten constant Km (mg/L)")                           # Long 2017 Table 4, Km row (final model)

    # Body-weight allometric covariate exponents (reference 70 kg).
    e_wt_cl <- 0.768;  label("Power exponent of WT/70 on CL and Q (unitless)")                   # Long 2017 Table 4, footnotes c and d (shared exponent)
    e_wt_v  <- 0.498;  label("Power exponent of WT/70 on V1 and V2 (unitless)")                  # Long 2017 Table 4, footnotes e and f (shared exponent)

    # Inter-patient variability. Long 2017 places a single eta on total CL
    # (CLtot = CL + Vmax/(C+Km)); applying the same eta to both linear CL
    # and Vmax makes CLtot vary exponentially as a whole (28.8% CV). IIV
    # terms for V1 (21.1% CV) and V2 (55.4% CV) are separate. Correlation
    # between the CLtot and V1 etas is 0.609; V2 is uncorrelated.
    #   omega^2 = log(CV^2 + 1)
    #   CLtot CV 28.8% -> omega^2 = log(0.288^2 + 1) = 0.07969
    #   V1    CV 21.1% -> omega^2 = log(0.211^2 + 1) = 0.04356
    #   V2    CV 55.4% -> omega^2 = log(0.554^2 + 1) = 0.26767
    #   cov(CLtot, V1) = 0.609 * sqrt(0.07969 * 0.04356) = 0.03588
    etalcl + etalvc ~ c(0.07969, 0.03588, 0.04356)  # Long 2017 Table 4: CLtot IIV 28.8% CV, V1 IIV 21.1% CV, correlation 0.609 (final model)
    etalvp ~ 0.26767                                 # Long 2017 Table 4: V2 IIV 55.4% CV (final model)

    # Combined additive + proportional residual error (Long 2017 Table 4).
    # Paper reports additive SD = 10.8 ug/mL (= mg/L since concentration
    # units are equivalent) and proportional = 23.7% CV (SD = 0.237).
    addSd  <- 10.8;  label("Additive residual error (mg/L)")                                     # Long 2017 Table 4, additive residual (final model)
    propSd <- 0.237; label("Proportional residual error (fraction)")                             # Long 2017 Table 4, proportional residual (final model)
  })
  model({
    # Individual PK parameters. A single eta (etalcl) is applied to both
    # linear CL and Vmax so that the total clearance CLtot = CL + Vmax/(C+Km)
    # inherits the 28.8% CV reported on CLtot. Weight effects from Long 2017
    # Table 4: CL and Q share exponent 0.768; V1 and V2 share exponent 0.498.
    cl   <- exp(lcl   + etalcl) * (WT / 70)^e_wt_cl
    vmax <- exp(lvmax + etalcl)
    km   <- exp(lkm)
    vc   <- exp(lvc   + etalvc) * (WT / 70)^e_wt_v
    q    <- exp(lq)             * (WT / 70)^e_wt_cl
    vp   <- exp(lvp   + etalvp) * (WT / 70)^e_wt_v

    # Two-compartment IV input model with parallel linear and
    # Michaelis-Menten elimination from the central compartment.
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q  / vc) * central +
                          (q  / vp) * peripheral1
    d/dt(peripheral1) <-  (q  / vc) * central - (q / vp) * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
