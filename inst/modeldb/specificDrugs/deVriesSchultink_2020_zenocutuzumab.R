deVriesSchultink_2020_zenocutuzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten non-linear elimination from the central compartment for intravenous zenocutuzumab (MCLA-128), a bispecific IgG1 (anti-HER2 x anti-HER3) monoclonal antibody, in patients with various advanced solid tumors (de Vries Schultink 2020)"
  reference <- "de Vries Schultink AHM, Bol K, Doornbos RP, Murat A, Wasserman E, Dorlo TPC, Schellens JHM, Beijnen JH, Huitema ADR. Population Pharmacokinetics of MCLA-128, a HER2/HER3 Bispecific Monoclonal Antibody, in Patients with Solid Tumors. Clin Pharmacokinet. 2020;59(7):875-884. doi:10.1007/s40262-020-00858-2"
  vignette <- "deVriesSchultink_2020_zenocutuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Baseline fat-free mass (Janmahasatian et al. 2005 formula, sex-specific)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on linear CL (exponent 1.19) and central volume V1 (exponent 0.71); reference 43.1 kg per the population median reported in de Vries Schultink 2020 Table 1. Computed sex-specifically from body weight, height, and sex per the Janmahasatian formulae quoted in Section 2.3.3 of the paper. After the FFM effect was retained, sex showed no further significant impact on PK and was not retained as a covariate.",
      source_name        = "FFM"
    ),
    TUM_SLD = list(
      description        = "Baseline sum of longest diameters of target lesions (RECIST 1.1)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on the maximum non-linear (Michaelis-Menten) clearance capacity Vmax (exponent 0.447); reference 70.0 mm per the population median reported in de Vries Schultink 2020 Table 1. The paper labels this column 'Sum of lesions (SoL)'; the canonical column name is TUM_SLD (RECIST 1.1 sum of longest diameters). HER2 status was tested but not retained.",
      source_name        = "SoL"
    )
  )

  population <- list(
    n_subjects        = 116L,
    n_studies         = 1L,
    n_observations    = 1115L,
    phase_mix         = "Pooled phase I dose-escalation (40-900 mg flat doses q3w, n = 23) and phase II dose-expansion (750 mg flat dose q3w, n = 93) cohorts of NCT02912949.",
    age_median        = "59 years",
    age_range         = "25-83 years",
    weight_median     = "68.2 kg",
    weight_range      = "40.2-112 kg",
    height_median     = "164 cm",
    height_range      = "147-199 cm",
    bsa_median        = "1.78 m^2",
    bsa_range         = "1.30-2.33 m^2",
    ffm_median        = "43.1 kg",
    ffm_range         = "28.6-72.45 kg",
    tum_sld_median    = "70.0 mm",
    tum_sld_range     = "12.0-266 mm",
    sex_female_pct    = 65.5,
    her2_status       = c(positive = 29.3, negative = 42.2, unknown = 28.5),
    tumor_type        = c(breast = 14.7, colorectal = 7.8, endometrium = 11.2, gastric = 21.5, lung = 6.9, ovarian = 31.0, others = 6.9),
    disease_state     = "Advanced solid tumors (HER2+ gastric, HER2+ breast, endometrium, esophagus-gastric junction, colon, NSCLC, ovarian, and other) eligible for the dose-escalation or HER2-overexpressing dose-expansion cohorts.",
    dose_range        = "Flat IV doses of 40, 80, 160, 240, 360, 480, 600, 750, or 900 mg administered every 3 weeks (q3w). 1-h infusion at <= 360 mg, 2-h infusion at > 360 mg. Of 116 patients, 93 received the 750 mg q3w expansion-cohort dose.",
    regions           = "Seven participating centers (NCT02912949).",
    reference_subject = "Population-median FFM 43.1 kg, sum of longest diameters of target lesions 70.0 mm (per de Vries Schultink 2020 Table 1).",
    notes             = "Baseline demographics per de Vries Schultink 2020 Table 1 and Results Section 3.1. Height was missing for one female and one male patient; the population-median height (164 cm) was imputed for both. Continuous-covariate missing values were imputed by the population median. The drug is a humanized full-length bispecific IgG1 antibody (MCLA-128 / zenocutuzumab) targeting HER2 and HER3; bispecifics are previously out-of-scope in nlmixr2lib's modality table, and this model is the first bispecific in the library."
  )

  ini({
    # Structural parameters (de Vries Schultink 2020 Table 2, "Covariate
    # model" column / SIR final point estimates). Reference subject:
    # population-median FFM (43.1 kg) and population-median sum of longest
    # diameters of target lesions (70.0 mm). Time unit hour, dose unit mg,
    # concentration unit mg/L (= ug/mL); LLOQ 50 ng/mL = 0.05 mg/L per
    # Section 2.2. The Cc state is central / vc with central in mg and vc
    # in L, giving Cc in mg/L. Vmax is mg/h and Km is mg/L, so the
    # Michaelis-Menten flux Vmax * Cc / (Km + Cc) is mg/h and matches the
    # linear-elimination flux kel * central (mg/h).
    lcl   <- log(0.0304); label("Linear clearance for the reference subject (L/h)")          # de Vries Schultink 2020 Table 2 (CL, covariate model)
    lvc   <- log(3.52);   label("Central volume of distribution V1 for the reference subject (L)") # de Vries Schultink 2020 Table 2 (V1, covariate model)
    lq    <- log(0.0254); label("Intercompartmental clearance Q (L/h)")                      # de Vries Schultink 2020 Table 2 (Q, covariate model)
    lvp   <- log(1.63);   label("Peripheral volume of distribution V2 (L)")                  # de Vries Schultink 2020 Table 2 (V2, covariate model)
    lvmax <- log(0.114);  label("Maximum non-linear (Michaelis-Menten) elimination rate Vmax for the reference subject (mg/h)") # de Vries Schultink 2020 Table 2 (Vmax, covariate model)
    lkm   <- log(0.211);  label("Michaelis-Menten constant Km (mg/L)")                       # de Vries Schultink 2020 Table 2 (Km, covariate model)

    # Covariate effects on linear CL and central volume V1 (FFM) and on
    # the maximum non-linear clearance capacity Vmax (TUM_SLD). All three
    # enter as power effects of the form (COV / median(COV))^exponent
    # per the continuous-covariate equation in Section 2.3.3.
    e_ffm_cl     <- 1.19;  label("Power exponent of FFM on linear CL (unitless)")            # de Vries Schultink 2020 Table 2 (FFM on CL)
    e_ffm_vc     <- 0.71;  label("Power exponent of FFM on central volume V1 (unitless)")    # de Vries Schultink 2020 Table 2 (FFM on V1)
    e_tumsld_vmax <- 0.447; label("Power exponent of sum of longest diameters of target lesions on Vmax (unitless)") # de Vries Schultink 2020 Table 2 (SoL on Vmax)

    # Inter-individual variability. de Vries Schultink 2020 Table 2 reports
    # omega as %CV on log-normal parameters; convert via
    # omega^2 = log(CV^2 + 1). CL and V1 share a correlated block with
    # correlation 0.55 (paper "Correlation omega_CL ~ omega_V1"); the
    # block covariance is cov = 0.55 * sqrt(omega^2_CL * omega^2_V1).
    #   omega^2_CL    = log(0.379^2 + 1) = 0.134263
    #   omega^2_V1    = log(0.210^2 + 1) = 0.043147
    #   omega^2_Vmax  = log(0.631^2 + 1) = 0.335229
    #   cov(CL,V1)    = 0.55 * sqrt(0.134263 * 0.043147) = 0.041863
    etalcl + etalvc ~ c(0.134263,
                        0.041863, 0.043147)   # CL 37.9% CV, V1 21.0% CV, corr 0.55 -- de Vries Schultink 2020 Table 2
    etalvmax ~ 0.335229                       # Vmax 63.1% CV -- de Vries Schultink 2020 Table 2

    # Residual error. de Vries Schultink 2020 Section 2.3.2 specifies a
    # combined proportional + additive error model with the additive term
    # fixed to LLOQ/2 (0.025 mg/L) per Beal 2001 method (Table 2 reports
    # "Add error (SD mg/L) 0.025 fixed").
    propSd <- 0.187;        label("Proportional residual error (fraction)")  # de Vries Schultink 2020 Table 2 (Prop. error 18.7% CV)
    addSd  <- fixed(0.025); label("Additive residual error (mg/L)")          # de Vries Schultink 2020 Table 2 (Add error 0.025 mg/L FIXED, = LLOQ/2 per Section 2.3.2)
  })

  model({
    # Individual PK parameters with FFM power effects on CL and V1 and a
    # TUM_SLD (sum of longest diameters of target lesions) power effect on
    # Vmax (de Vries Schultink 2020 final covariate model, Table 2 and
    # Results Section 3.2). Reference values are the population medians:
    # FFM 43.1 kg, TUM_SLD 70.0 mm.
    cl   <- exp(lcl   + etalcl)   * (FFM     / 43.1)^e_ffm_cl
    vc   <- exp(lvc   + etalvc)   * (FFM     / 43.1)^e_ffm_vc
    q    <- exp(lq)
    vp   <- exp(lvp)
    vmax <- exp(lvmax + etalvmax) * (TUM_SLD / 70.0)^e_tumsld_vmax
    km   <- exp(lkm)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment IV model with parallel linear and Michaelis-Menten
    # non-linear elimination from the central compartment, per the ODE
    # system in de Vries Schultink 2020 Results Section 3.2:
    #   dA1/dt = -CL/V1 * A1 - Vmax * C1 / (Km + C1) - Q/V1 * A1 + Q/V2 * A2
    #   dA2/dt =  Q/V1 * A1 - Q/V2 * A2
    # with C1 = A1 / V1. Rewriting CL/V1 = kel, Q/V1 = k12, Q/V2 = k21:
    d/dt(central)     <- -kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                          k12 * central - k21 * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}
