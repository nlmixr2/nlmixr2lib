FernandezRubio_2025_ceftriaxone <- function() {
  description <- "Two-compartment population PK model for unbound (free) intravenous ceftriaxone in elderly patients (>55 years) receiving high-dose ceftriaxone plus ampicillin for Enterococcus faecalis infective endocarditis. Fitted directly to ultrafiltrate-measured FREE ceftriaxone concentrations against the total administered dose, so clearance and both volumes are apparent unbound-drug parameters roughly an order of magnitude larger than the corresponding total-drug values. Creatinine clearance scales clearance and body mass index scales the central volume, both as median-normalised power terms. Estimated in Monolix 2024R1 by SAEM. Fernandez Rubio 2025, n = 16 patients / 24 treatment episodes, 3 samples per episode (pre-dose, +2 h, +4 h) at steady state."
  reference <- "Fernandez Rubio B, Docobo Perez F, Herrera Hidalgo L, Lopez-Cortes LE, Luque Marquez R, Lomas Cabezas JM, Lopez-Cortes LF, Mejias Trueba M, Guisado Gil AB, Gutierrez Valencia A, de Alarcon Gonzalez A, Gil Navarro MV. High-Dose Ceftriaxone in Elderly Patients with Enterococcal Infective Endocarditis: Population Pharmacokinetics of Free Ceftriaxone and Dose Optimization. Antibiotics. 2025;14(5):508. doi:10.3390/antibiotics14050508"
  vignette <- "FernandezRubio_2025_ceftriaxone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation and standardised to a body surface area of 1.73 m^2",
      units              = "mL/min/1.73m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per treatment episode in the source analysis. Fernandez Rubio 2025 Methods Sect. 4.3 states creatinine clearance was estimated with Cockcroft-Gault (reference [37], Cockcroft & Gault 1976) and 'standardized for 1.73 m2', so this column carries the BSA-normalised variant of the CRCL canonical rather than a raw mL/min value. Enters as the median-normalised power term (CRCL / 59.5)^0.62 on clearance (Results Sect. 2.2, first displayed equation). The 59.5 centring constant is the cohort median in Table 1 (IQR 48.5-88.4), so unlike several sibling entries in this register the printed equation constant and the tabulated median agree exactly. Eligibility excluded serum creatinine > 1.5 mg/dL and CrCl < 10 mL/min (Methods Sect. 4.2), so the model carries no information about severe renal impairment; the cohort is mildly-to-moderately impaired, consistent with its age (median 77 years). The bootstrap IQR of the exponent spans zero (0.5, -0.32-1.1), which the authors nonetheless retained on the -2LL / BIC criterion; the point estimate is carried here as printed.",
      source_name        = "CrCl"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per treatment episode. Enters as the median-normalised power term (BMI / 29.3)^2.52 on the central volume of distribution (Fernandez Rubio 2025 Results Sect. 2.2, second displayed equation). The 29.3 centring constant is the cohort median in Table 1 (IQR 26.7-33.3, an overweight-to-obese elderly cohort), matching the printed equation exactly. The exponent 2.52 is unusually steep for a body-size term on volume -- far above the allometric 1 -- and it is estimated on 24 episodes with 37.3% RSE and a bootstrap IQR that spans zero (2.83, -0.26-6.13); it is carried here exactly as printed but should not be extrapolated outside the observed 26.7-33.3 interquartile band without care. Body weight, age and sex were also screened as covariates (Methods Sect. 4.5) but were not retained; see covariatesDataExcluded.",
      source_name        = "BMI"
    )
  )

  # Screened in the source's forward-inclusion / backward-deletion covariate
  # search (Fernandez Rubio 2025 Methods Sect. 4.5: "Creatinine clearance
  # standardized for 1.73 m2 (CrCl), age, weight, body mass index (BMI), and
  # sex were explored as covariates for each structural model") but NOT
  # retained in the final model, which reports only the CrCl-on-CL and
  # BMI-on-V1 effects (Table 2). No point estimate is published for any of
  # these, so they are documentation only.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a covariate on each structural parameter (Methods Sect. 4.5) but not retained in the final model. Cohort median 90 kg (IQR 69.8-99.5), Table 1. Body size enters the final model only through BMI on V1."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a covariate (Methods Sect. 4.5) but not retained. Cohort median 77 years (IQR 71-78), Table 1; the enrolment floor was 55 years (Methods Sect. 4.2). Age is effectively constant across this cohort, which is itself the reason the study exists -- the paper's framing is that elderly-specific ceftriaxone PK is under-studied, not that age is a within-cohort covariate."
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "The only categorical covariate screened (Methods Sect. 4.5), tested through the exponential form theta_i = theta_pop * exp(beta_sex * SEX) with the source coding SEX 0 = female, 1 = male. Not retained in the final model, so no coefficient is published. Recorded here on the canonical SEXF orientation (1 = female), which is the inverse of the source's coding; had the effect been retained, the sign would need to be flipped. Cohort 18 of 24 episodes male (75%), Table 1."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "ceftriaxone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftriaxone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 16L,
    n_episodes       = 24L,
    n_studies        = 1L,
    n_samples        = 72L,
    age_median       = "77 years (IQR 71-78); enrolment floor 55 years",
    weight_median    = "90 kg (IQR 69.8-99.5)",
    height_median    = "168 cm (IQR 161-180)",
    bmi_median       = "29.3 kg/m^2 (IQR 26.7-33.3)",
    # Patient-level, to match n_subjects: 10 of 16 patients (63%) were male,
    # so 6 of 16 (37.5%) were female. Fernandez Rubio 2025 Table 1 reports the
    # EPISODE-level figure instead (18 of 24 episodes, 75%, male -> 25%
    # female); the two differ because a patient who received more than one
    # regimen contributed more than one episode.
    sex_female_pct   = 37.5,
    disease_state    = "Enterococcus faecalis infective endocarditis caused by ampicillin-susceptible strains (24/24 episodes), treated with ampicillin plus high-dose ceftriaxone. Baseline laboratory values (median, IQR): serum protein 5.8 g/dL (5.5-6.5), serum creatinine 1.1 mg/dL (0.8-1.2). Ten of the 16 patients (63%) were male, accounting for 18 of the 24 episodes (75%).",
    renal_function   = "Creatinine clearance (Cockcroft-Gault, standardised to 1.73 m^2) median 59.5 mL/min/1.73 m^2 (IQR 48.5-88.4). Exclusion criteria removed serum creatinine > 1.5 mg/dL and severe renal impairment (CrCl < 10 mL/min), so no dialysis or severely impaired patient is represented.",
    dose_range       = "Ceftriaxone at least 4 g/day for more than 48 h: 2 g every 12 h in 18 episodes (75%), 4 g every 24 h in 3 episodes (12.5%), and 6 g every 24 h in 3 episodes (12.5%). All episodes also received ampicillin, which is not part of this model.",
    protein_binding  = "Measured directly rather than assumed. Mean plasma protein binding was 85.7 +/- 6.4% pre-dose, 74.7 +/- 10.1% at +2 h, and 79.6 +/- 9.3% at +4 h -- markedly lower (i.e. a higher free fraction) than the 79-94% quoted for healthy young individuals, which the authors attribute to the hypoproteinaemia of ageing plus binding saturation at these high doses. Mean TOTAL ceftriaxone concentrations were 51.6 +/- 21.9 mg/L pre-dose, 129 +/- 53.2 mg/L at +2 h and 105.2 +/- 49.4 mg/L at +4 h; the corresponding mean FREE concentrations, which are what this model was fitted to, were 7.8 +/- 6.5, 34 +/- 26.5 and 22.7 +/- 19.7 mg/L.",
    regions          = "Spain (two tertiary teaching hospitals in Seville), 2021-2022",
    notes            = "Prospective observational PK study. Three blood samples were drawn per treatment episode at steady state (at least 48 h after treatment start): immediately pre-dose (Cmin), 2 +/- 0.5 h after the dose (C2) and 4 +/- 0.5 h after the dose (C4); a patient who received more than one regimen contributed one sample set per regimen, which is why 16 patients yield 24 episodes. Total and free ceftriaxone were quantified by LC-MS/MS, the free fraction isolated by 37 C ultrafiltration (Amicon Ultra 0.5 mL 30 K); the free standard curve was linear over 0.5-200 mg/L with a 0.5 mg/L lower limit of quantification. Albumin was NOT recorded (a stated limitation); total plasma protein was recorded instead and was not tested as a covariate."
  )

  ini({
    # Structural parameters: Fernandez Rubio 2025 Table 2, "Mean (%RSE)"
    # column of the final free-ceftriaxone model. Table 2 also reports a
    # 1000-iteration nonparametric bootstrap "Median (IQR)" column, which is
    # noted per line but is not the value carried here.
    #
    # These are UNBOUND-drug parameters fitted against the TOTAL administered
    # dose, so both volumes and the clearance carry the reciprocal of the free
    # fraction: they are approximately fu^-1 times the corresponding
    # total-drug values quoted in the Background (Vd 10.69-11.01 L,
    # CL 833-1023 mL/h). At the cohort's measured free fraction of roughly
    # 0.15-0.25 this places apparent unbound Vc near 44-73 L and apparent
    # unbound CL near 3-7 L/h, which brackets the estimates below.

    lcl <- log(11.57)
    label("Apparent unbound clearance CL at the reference CrCl of 59.5 mL/min/1.73m^2 (L/h)")
    # Fernandez Rubio 2025 Table 2: Cl 11.57 L/h (10.5% RSE);
    # bootstrap median 11.33 (IQR 9.04-14.07); eta-shrinkage 0.51%

    lvc <- log(43.6)
    label("Apparent unbound central volume V1 at the reference BMI of 29.3 kg/m^2 (L)")
    # Fernandez Rubio 2025 Table 2: V1 43.6 L (27.5% RSE);
    # bootstrap median 33.9 (IQR 15.01-58.03); eta-shrinkage 14.16%

    lq <- log(19.8)
    label("Apparent unbound intercompartmental clearance Q (L/h)")
    # Fernandez Rubio 2025 Table 2: Q 19.8 (32.4% RSE); bootstrap median
    # 20.97 (IQR 11.31-49.15); eta-shrinkage 26.38%. Table 2 prints the
    # unit for Q as "h-1"; that is a typographical error -- Q is defined in
    # the Results text and in the Table 2 footnote as the
    # INTERCOMPARTMENTAL CLEARANCE, which in a CL/V1/Q/V2 parameterisation
    # has units of L/h. Reading 19.8 as a rate constant instead would make
    # the model dimensionally inconsistent. See vignette Errata.

    lvp <- log(40.94)
    label("Apparent unbound peripheral volume V2 (L)")
    # Fernandez Rubio 2025 Table 2: V2 40.94 L (12.2% RSE);
    # bootstrap median 47.17 (IQR 39.94-58.67); eta-shrinkage 43.28%

    # Covariate effects. Fernandez Rubio 2025 Results Sect. 2.2 prints the
    # two final-model covariate relationships as median-normalised power
    # terms (the paper calls this a "log-linear model", i.e. linear in the
    # logarithm of the parameter):
    #     Cl_i = Cl_pop * (CrCl_i / 59.5)^beta_CrCl
    #     V1_i = V1_pop * (BMI_i  / 29.3)^beta_BMI
    # matching the general form given in Methods Sect. 4.5,
    #     theta_i = theta_pop * (Cov_i / Cov_median)^beta_cov,
    # where beta_cov "represents the estimated exponent describing the
    # covariate's effect". The two centring constants are the Table 1
    # cohort medians (CrCl 59.5 mL/min/1.73m^2, BMI 29.3 kg/m^2).
    e_crcl_cl <- 0.62
    label("Power exponent of creatinine clearance on CL (unitless)")
    # Fernandez Rubio 2025 Table 2, "Effect of CrCl on CL": 0.62 (44.8% RSE);
    # bootstrap median 0.5 (IQR -0.32-1.1)

    e_bmi_vc <- 2.52
    label("Power exponent of body mass index on V1 (unitless)")
    # Fernandez Rubio 2025 Table 2, "Effect of BMI on V1": 2.52 (37.3% RSE);
    # bootstrap median 2.83 (IQR -0.26-6.13)

    # Between-subject variability. Methods Sect. 4.5: "The between-subject
    # variability (BSV or omega) was ascribed to an exponential
    # distribution", i.e. theta_i = theta_pop * exp(eta_i) with
    # eta_i ~ N(0, omega^2). The Table 2 footnote glosses omega as
    # "coefficient of variation for between-subject variability", but the
    # tabulated quantity is Monolix's omega_<parameter> output, which is the
    # STANDARD DEVIATION of the random effect on the log scale, not a
    # variance and not a percentage CV. Two independent checks support the
    # SD reading, and both are reproduced in the vignette:
    #   (1) Software convention -- Monolix 2024R1 (the stated estimation
    #       tool) reports omega as the SD of eta; it has no variance-scale
    #       output. The loose "CV" gloss is the standard abuse of language
    #       for a log-normal random effect, where CV = sqrt(exp(omega^2)-1)
    #       approaches omega for small omega.
    #   (2) The paper's own Table 3 -- reconstructing the six Monte Carlo
    #       dosing scenarios reproduces the published probabilities of
    #       target attainment only under the SD reading. Under the
    #       variance reading (omega = sqrt(0.48) = 0.69 on CL) the
    #       continuous-infusion PTAs are badly under-predicted; the
    #       discriminating cell is 6 g/24 h by 24 h infusion at a 10 mg/L
    #       target, published at 92.0%.
    # Encoded here as variances, omega^2:
    #   Cl : 0.48^2 = 0.2304
    #   V1 : 0.77^2 = 0.5929
    #   Q  : 0.75^2 = 0.5625
    #   V2 : 0.09^2 = 0.0081
    # The variance-covariance structure of the BSV was itself a modelled
    # hypothesis (Methods Sect. 4.5) but Table 2 reports no off-diagonal
    # terms, so the etas are taken as uncorrelated.
    etalcl ~ 0.2304  # Fernandez Rubio 2025 Table 2: omega Cl 0.48 (15.9% RSE); bootstrap median 0.47 (IQR 0.31-0.6)
    etalvc ~ 0.5929  # Fernandez Rubio 2025 Table 2: omega V1 0.77 (20.6% RSE); bootstrap median 0.82 (IQR 0.12-1.26)
    etalq  ~ 0.5625  # Fernandez Rubio 2025 Table 2: omega Q  0.75 (40.8% RSE); bootstrap median 0.83 (IQR 0.33-1.38)
    etalvp ~ 0.0081  # Fernandez Rubio 2025 Table 2: omega V2 0.09 (68.5% RSE); bootstrap median 0.12 (IQR 0.044-0.26)

    # Residual error. Methods Sect. 4.5 states additive, proportional,
    # combined-normal and log-normal residual models were all tested; the
    # Table 2 footnote records the retained one as "sigma, constant error to
    # ceftriaxone observations". "Constant" is Monolix's name for the purely
    # ADDITIVE error model y = f + a * e, so sigma is the additive standard
    # deviation in the concentration units of the observations (mg/L of free
    # ceftriaxone). No proportional term was retained.
    addSd <- 1.41
    label("Additive residual error on free ceftriaxone concentration (mg/L)")
    # Fernandez Rubio 2025 Table 2: sigma 1.41 (20% RSE);
    # bootstrap median 1.24 (IQR 0.5-1.95)
  })

  model({
    # Covariate centring constants: the cohort medians used to normalise
    # each covariate in the Results Sect. 2.2 equations, equal to the
    # Table 1 medians.
    crcl_ref <- 59.5  # mL/min/1.73m^2, Fernandez Rubio 2025 Table 1 median
    bmi_ref  <- 29.3  # kg/m^2,         Fernandez Rubio 2025 Table 1 median

    # Individual parameters. Ceftriaxone is given intravenously, so there is
    # no absorption step: the Results Sect. 2.2 sentence "a two-compartment
    # model with first-order absorption and elimination" is a boilerplate
    # slip -- the model's own listed parameters are CL, V1, Q and V2 with no
    # absorption rate constant or bioavailability term anywhere in Table 2,
    # and every simulated regimen in Table 3 is an intravenous infusion.
    # See vignette Errata.
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (BMI / bmi_ref)^e_bmi_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # The observation is the UNBOUND (free) ceftriaxone plasma
    # concentration. The model was fitted to ultrafiltrate concentrations
    # against the total administered dose, so no protein-binding term
    # appears: the free fraction is absorbed into the apparent volumes.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
