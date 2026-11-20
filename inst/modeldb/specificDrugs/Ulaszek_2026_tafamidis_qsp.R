Ulaszek_2026_tafamidis_qsp <- function() {
  description <- paste(
    "QSP. Tafamidis-transthyretin (TTR) binding and TTR turnover model in",
    "healthy adult men (in silico; 30 virtual subjects aged 30-40 years).",
    "Six ODEs (Ulaszek 2026 Eqs 5-10): free TTR tetramer with zero-order",
    "synthesis and first-order degradation, reversible dissociation into",
    "monomers (fourth-order reassociation, first-order monomer loss),",
    "sequential binding of tafamidis to the two thyroxine sites of the",
    "tetramer (single- and double-bound complexes), and one-site binding of",
    "tafamidis to albumin. Free tafamidis is the total plasma concentration",
    "minus the drug held in the three complexes (Eq 11). The two",
    "tafamidis-bound TTR complexes are cleared at a fixed fraction (0.735) of",
    "the free-TTR degradation rate, the paper's calibrated hypothesis for the",
    "clinically observed ~33% rise in circulating TTR under tafamidis.",
    "The tafamidis concentration is NOT computed here: Ulaszek 2026 generated",
    "it with a Simcyp V23 whole-body PBPK model with ADAM absorption whose",
    "organ volumes, blood flows and absorption rate are Simcyp database",
    "outputs that the paper does not publish, so that layer is not",
    "reproducible from the published sources and is deliberately NOT",
    "extracted. Total plasma tafamidis is instead supplied per record as the",
    "time-varying covariate CP_TAFAMIDIS_UM (uM). Binding constants default",
    "to the Nelson-based set used for the paper's Figure 3 and its deposited",
    "Simcyp script; the alternative Corazza-based dissociation constants",
    "(koff1 28.8 1/h, koff2 216 1/h) are documented in the parameter labels.",
    "Deterministic mechanism model: variability enters through per-subject",
    "TTR synthesis/degradation rates and albumin, not through fitted IIV, so",
    "no etas and no residual-error model are encoded.",
    sep = " "
  )
  reference <- paste(
    "Ulaszek S, Lisowski B, Wisniowska B, Polak S. Exploring Tafamidis",
    "Effects Through PBPK-QSP Modelling. Pharmaceutics. 2026;18(3):367.",
    "doi:10.3390/pharmaceutics18030367. PMCID PMC13029652.",
    "Equations from Section 2.3.2 (Eqs 5-11) and Appendix A Eq A1; parameter",
    "values from Appendix A Table A2, confirmed against the deposited Simcyp",
    "Lua script and the 'PD Individual Parameters' sheet of the",
    "supplementary simulation workbook (Supplementary Materials).",
    sep = " "
  )
  vignette <- "Ulaszek_2026_tafamidis"

  units <- list(
    time = "h",
    dosing = "(not applicable; total plasma tafamidis is supplied as the time-varying covariate CP_TAFAMIDIS_UM)",
    concentration = "umol/L (uM) for every state and output"
  )

  # Paper-specific QSP species; none is a registered compartment. TTR is a
  # tetramer with two thyroxine sites, so the drug-bound forms are numbered
  # by the number of tafamidis molecules bound.
  paper_specific_compartments <- c(
    "complex_ttr1",
    "complex_ttr2",
    "ttr_free",
    "ttr_monomer",
    "complex_albumin",
    "albumin"
  )

  compartmentData <- list(
    complex_ttr1 = list(
      analyte = "TTR tetramer with one bound tafamidis molecule",
      units = "umol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    complex_ttr2 = list(
      analyte = "TTR tetramer with two bound tafamidis molecules",
      units = "umol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    ttr_free = list(analyte = "unbound TTR tetramer", units = "umol/L", specimen = "plasma", verified = TRUE),
    ttr_monomer = list(analyte = "TTR monomer", units = "umol/L", specimen = "plasma", verified = TRUE),
    complex_albumin = list(
      analyte = "tafamidis-albumin complex",
      units = "umol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    albumin = list(analyte = "unbound albumin", units = "umol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CP_TAFAMIDIS_UM = list(
      description = "Time-varying TOTAL (bound plus unbound) tafamidis plasma concentration driving TTR and albumin binding",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Total plasma concentration (the Simcyp 'PD Input Type: Total",
        "concentration, PKPD Input Compartment: Plasma' of the deposited",
        "workbook); free tafamidis is derived inside model() by Eq 11.",
        "Tafamidis free-acid molecular weight 308.12 g/mol (Table A1), so",
        "1 uM = 308.12 ng/mL; divide a concentration in ng/mL by 308.12.",
        "In Ulaszek 2026 the concentration came from a Simcyp V23 whole-body",
        "PBPK model that is not extracted (see description). Published",
        "predicted anchors (Table 2) for 80 mg tafamidis meglumine once daily",
        "on day 7: mean Cmax 9160 ng/mL (29.7 uM) and mean AUC144-168h",
        "180,557 ng*h/mL; after a single 80 mg dose (fasted) geometric-mean",
        "Cmax 3586 ng/mL (11.6 uM). Set to 0 for drug-free records."
      ),
      source_name = "xin (Simcyp Lua odeRateStep input); Ctotal in Eq 11"
    ),
    ALB = list(
      description = "Individual serum albumin concentration, the total albumin pool available for tafamidis binding",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converted inside model() to umol/L with the albumin molecular weight",
        "66,500 g/mol, exactly as the deposited Simcyp script does",
        "(getIndivHSA() / 66500 * 1e6). Sets the initial free-albumin state;",
        "the albumin pool is conserved (albumin + complex_albumin = ALB).",
        "Ulaszek 2026 Table A2: Simcyp healthy-volunteer albumin model, a",
        "baseline of 50.34 g/L adjusted for age (coefficient -0.0575) and BMI",
        "(coefficient -0.0738) with 10% inter-individual CV. The 30 simulated",
        "subjects of the supplementary workbook span 555-833 umol/L",
        "(36.9-55.4 g/L), mean 681 umol/L (45.3 g/L)."
      ),
      source_name = "HSA (Simcyp individual human serum albumin, g/L)"
    )
  )

  population <- list(
    species = "human (in silico; Simcyp V23 healthy-volunteer virtual population)",
    n_subjects = 30L,
    n_studies = NA_integer_,
    age_range = "30-40 years",
    sex_female_pct = 0,
    disease_state = "Healthy volunteers (virtual). The tafamidis PBPK component was evaluated against three phase I studies in healthy volunteers: the thorough-QT study NCT01775761 (n = 42, single 400 mg), the food-effect study B3461054 (n = 16, single 80 mg fasted and fed) and the bioequivalence study B3461056 (n = 30 men, 80 mg once daily for 7 days).",
    dose_range = "Tafamidis meglumine 80 mg (4 x 20 mg; 48.8 mg tafamidis) once daily is the reference regimen; dose-sensitivity simulations at 20, 40, 60, 80 and 160 mg once daily",
    regions = NA_character_,
    notes = paste(
      "Baseline TTR for each virtual subject was sampled from a normal",
      "distribution using the mean and SD that Ingenbleek et al. report for",
      "men aged 30-40 years, and split into a synthesis rate r and a",
      "degradation rate Tdeg (TTR = r / Tdeg) by an optimiser started from",
      "the Oppenheimer et al. turnover estimates (Section 2.4; supplementary",
      "R script). The resulting 30 baseline TTR values span 3.06-7.51 uM",
      "(mean 5.54 uM). Women were deliberately not simulated. The QSP",
      "component was not fitted to individual data; its only calibrated",
      "quantity is the bound-TTR clearance fraction, set to reproduce the",
      "population-average 33% rise in total TTR reported in the literature."
    )
  )

  ini({
    # System (patient) parameters. Ulaszek 2026 assigns every virtual subject
    # its own r and Tdeg (Table A2 'Patient-specific'); the values below are
    # the Oppenheimer et al. population averages that seeded that split
    # (supplementary R script '1_TTR sampling.R': 'avg mean conc 6.12 uM, avg
    # deg 0.015125/hour, avg synth 0.092565 uM/h'). Per-subject values from
    # the paper's deposited Simcyp script are tabulated in the vignette.
    ksyn <- fixed(0.092565)
    label("Zero-order TTR tetramer synthesis rate r (uM/h)") # supplement R script initial_params[1] = 0.092565 (Oppenheimer); Table A2 'TTR synthesis rate' = patient-specific
    kdeg <- fixed(0.015125)
    label("First-order degradation rate constant of free TTR tetramer Tdeg (1/h)") # Section 2.4 Oppenheimer 0.363 1/day = 0.015125 1/h; supplement R script initial_params[2] = 0.015125

    # TTR tetramer-monomer kinetics (Table A2)
    kdiss <- fixed(0.0024)
    label("Dissociation rate constant of unbound TTR tetramer into monomers kdiss (1/h)") # Table A2 'Dissociation rate of unbound TTR' = 0.0024 1/h
    kreassoc <- fixed(360000)
    label("Fourth-order reassociation rate constant of TTR monomers into tetramer, paper symbol ka (1/(uM^3*h))") # Table A2 'Association rate of unbound TTR' = 360,000 1/(uM^3*h)
    kdeg_monomer <- fixed(0.016)
    label("First-order irreversible loss rate constant of TTR monomers Mdeg (1/h)") # Table A2 'Monomer degradation rate' = 0.016 1/h (assumed equal to the average TTR degradation rate)

    # Tafamidis-TTR binding (Table A2). koff1 / koff2 are the Nelson-based
    # set (Nelson 2021 KD combined with Corazza 2019 kon) that drives the
    # deposited Simcyp script and Figure 3. The Corazza 2019 set replaces
    # them with koff1 = 28.8 and koff2 = 216 1/h (Figure A2).
    kon1 <- fixed(16200)
    label("Association rate constant of the first tafamidis molecule to TTR kon1 (1/(uM*h))") # Table A2 'Association rate constant for first tafamidis-TTR binding' = 16,200
    koff1 <- fixed(50.2)
    label("Dissociation rate constant of the single-bound tafamidis-TTR complex koff1 (1/h); Corazza-based alternative 28.8") # Table A2 'Dissociation rate constant for first tafamidis-TTR binding' = 50.2 OR 28.8
    kon2 <- fixed(12600)
    label("Association rate constant of the second tafamidis molecule to TTR kon2 (1/(uM*h))") # Table A2 'Association rate constant for second tafamidis-TTR binding' = 12,600
    koff2 <- fixed(2998.8)
    label("Dissociation rate constant of the double-bound tafamidis-TTR complex koff2 (1/h); Corazza-based alternative 216") # Table A2 'Dissociation rate constant for second tafamidis-TTR binding' = 2998.8 OR 216

    # Tafamidis-albumin binding (Table A2); KD = koff_alb / kon_alb = 1.8 uM
    kon_alb <- fixed(1000)
    label("Association rate constant of tafamidis to albumin kon,albumin (1/(uM*h))") # Table A2 'Association rate constant for tafamidis-albumin binding' = 1000
    koff_alb <- fixed(1800)
    label("Dissociation rate constant of the tafamidis-albumin complex koff,albumin (1/h)") # Table A2 'Dissociation rate constant for tafamidis-albumin binding' = 1800

    # Turnover-scaling factor s applied to the degradation of both bound
    # complexes. The text and Table A2 give the analytic approximation
    # s = 1/1.33 ~ 0.75; the deposited Simcyp script that produced every
    # published simulation uses 0.735 (T_bound_deg = P[2]*0.735), which is
    # the value that returns the 33% rise in this full model (vignette).
    f_kdeg_bound <- fixed(0.735)
    label("Fraction of the free-TTR degradation rate applied to tafamidis-bound TTR complexes s (unitless)") # Simcyp Lua odeRateStep 'P[2]*0.735'; Section 2.3 Eq 4 / Table A2 'Tdeg x 0.75'
  })

  model({
    # Albumin molecular weight used by the deposited Simcyp script to convert
    # individual serum albumin (g/L) to umol/L.
    mw_albumin <- 66500

    # Irreversible elimination rate of both bound complexes (Table A2:
    # 'Tdeg x s'; Eqs 5-6 T_bound,deg)
    kdeg_bound <- kdeg * f_kdeg_bound

    # Eq 11: free tafamidis = total - drug bound to TTR (one molecule in the
    # single complex, two in the double complex) - drug bound to albumin
    ctaf <- CP_TAFAMIDIS_UM - (complex_ttr1 + 2 * complex_ttr2 + complex_albumin)

    # Eq 5: single-bound tafamidis-TTR complex
    d/dt(complex_ttr1) <- -kdeg_bound * complex_ttr1 + kon1 * ctaf * ttr_free - koff1 * complex_ttr1 + koff2 * complex_ttr2 - kon2 * ctaf * complex_ttr1
    # Eq 6: double-bound tafamidis-TTR complex
    d/dt(complex_ttr2) <- -kdeg_bound * complex_ttr2 + kon2 * ctaf * complex_ttr1 - koff2 * complex_ttr2
    # Eq 7: free TTR tetramer (synthesis, monomer reassociation, dissociation,
    # degradation, drug binding)
    d/dt(ttr_free) <- ksyn + kreassoc * ttr_monomer^4 - (kdiss + kdeg) * ttr_free - kon1 * ctaf * ttr_free + koff1 * complex_ttr1
    # Eq 8: TTR monomer (four monomers per tetramer)
    d/dt(ttr_monomer) <- 4 * kdiss * ttr_free - 4 * kreassoc * ttr_monomer^4 - kdeg_monomer * ttr_monomer
    # Eq 9: tafamidis-albumin complex
    d/dt(complex_albumin) <- kon_alb * ctaf * albumin - koff_alb * complex_albumin
    # Eq 10: free albumin
    d/dt(albumin) <- -kon_alb * ctaf * albumin + koff_alb * complex_albumin

    # Initial conditions (Simcyp odeInitStep): drug-free baseline with free
    # TTR at ksyn / kdeg (Eq 2), monomer from Eq A1, all albumin unbound
    ttr_free(0) <- ksyn / kdeg
    ttr_monomer(0) <- (kdiss / kreassoc * ksyn / kdeg)^(1 / 4)
    albumin(0) <- ALB * 1e6 / mw_albumin

    # Outputs: bound and total TTR tetramer (the paper's surrogate of drug
    # action) and the unbound-to-total TTR ratio of Figure 4
    ttr_bound <- complex_ttr1 + complex_ttr2
    ttr_total <- ttr_free + complex_ttr1 + complex_ttr2
    ttr_unbound_ratio <- ttr_free / ttr_total
  })
}
