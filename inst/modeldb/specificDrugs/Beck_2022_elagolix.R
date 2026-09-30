Beck_2022_elagolix <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order",
    "absorption, absorption lag time and first-order elimination for oral",
    "elagolix (a gonadotropin-releasing hormone receptor antagonist) in 2168",
    "premenopausal women pooled from six phase I studies in healthy women,",
    "four phase III studies in endometriosis and three phase III studies in",
    "uterine fibroids, where elagolix 300 mg BID was given alone or with",
    "estradiol 1 mg / norethindrone acetate 0.5 mg QD add-back therapy (Beck",
    "2022). OATP1B1 c.521T>C (rs4149056) transporter status acts on relative",
    "bioavailability (intermediate, poor and missing-genotype strata) and",
    "body weight acts on the apparent central volume by a power function.",
    "Combined proportional-plus-additive residual error with separate",
    "magnitudes for the phase I and phase III studies."
  )
  reference <- paste(
    "Beck D, Winzenborg I, Liu M, Degner J, Mostafa NM, Noertersheuser P,",
    "Shebley M. Population Pharmacokinetics of Elagolix in Combination with",
    "Low-Dose Estradiol/Norethindrone Acetate in Women with Uterine Fibroids.",
    "Clin Pharmacokinet. 2022;61(4):577-587. doi:10.1007/s40262-021-01096-w.",
    "Parameter estimates are from Table 3 and its footnotes c and d; the",
    "residual-error equation is Equation 2 of the Electronic Supplementary",
    "Material (40262_2021_1096_MOESM1_ESM.pdf)."
  )
  vignette <- "Beck_2022_elagolix"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the apparent central volume normalised to the",
        "population median of 76 kg (Table 3 footnote c:",
        "Vc/F = 279 * (WT/76)^0.160). Baseline value in the source."
      ),
      source_name = "Body weight"
    ),
    SNP_SLCO1B1_RS4149056_HET = list(
      description = "SLCO1B1 c.521T>C (rs4149056) heterozygous genotype, OATP1B1 'intermediate transporter'; 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive transporter, homozygous wild type, when SNP_SLCO1B1_RS4149056_HOM and SNP_SLCO1B1_RS4149056_MISSING are also 0)",
      notes = paste(
        "Beck 2022 Methods 2.4: intermediate transporter (IT) =",
        "heterozygous for 521T>C. Relative bioavailability 1.42 vs the",
        "extensive-transporter reference (Table 3 footnote d). 335 of 2168",
        "participants (15.5%, Table 2)."
      ),
      source_name = "OATP1B1 genotype status = intermediate transporter"
    ),
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = "SLCO1B1 c.521T>C (rs4149056) homozygous-variant genotype, OATP1B1 'poor transporter'; 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive transporter, homozygous wild type, when SNP_SLCO1B1_RS4149056_HET and SNP_SLCO1B1_RS4149056_MISSING are also 0)",
      notes = paste(
        "Beck 2022 Methods 2.4: poor transporter (PT) = homozygous variant",
        "521T>C. Relative bioavailability 1.96 vs the extensive-transporter",
        "reference (Table 3 footnote d). 32 of 2168 participants (1.48%,",
        "Table 2)."
      ),
      source_name = "OATP1B1 genotype status = poor transporter"
    ),
    SNP_SLCO1B1_RS4149056_MISSING = list(
      description = "SLCO1B1 c.521T>C (rs4149056) genotype not available; 1 = missing, 0 = genotyped",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotype known)",
      notes = paste(
        "Beck 2022 estimates a separate relative-bioavailability factor of",
        "1.10 for participants with no pharmacogenetic sample (Table 3",
        "footnote d); 545 of 2168 participants (25.1%, Table 2), mostly from",
        "non-consent to pharmacogenetic testing. When this indicator is 1,",
        "SNP_SLCO1B1_RS4149056_HET and SNP_SLCO1B1_RS4149056_HOM must be 0.",
        "The paper's covariate simulations did not use the missing stratum."
      ),
      source_name = "OATP1B1 genotype status = missing"
    ),
    STUDY_PHASE3 = list(
      description = "Record from a phase III study; 1 = phase III (sparse sampling), 0 = phase I (intensive sampling)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the six phase I studies in healthy women)",
      notes = paste(
        "Selects the residual-error magnitudes only (Table 3 'phase I",
        "studies' and 'phase III studies' rows); it changes no structural",
        "or covariate parameter, so the typical-value prediction is",
        "identical for both values. Set 1 to simulate observations like the",
        "endometriosis and uterine-fibroid trials, 0 for the intensively",
        "sampled phase I design."
      ),
      source_name = "Study phase (phase I vs phase III)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "elagolix", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "elagolix", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "elagolix", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2168L,
    n_studies = 13L,
    age_range = "18-53 years",
    age_median = "36 years",
    weight_range = "40-160 kg",
    weight_median = "76 kg",
    sex_female_pct = 100,
    race_ethnicity = c(Black = 30.4, `White and others` = 69.6),
    disease_state = paste(
      "Premenopausal women: 175 healthy (phase I), 1310 with",
      "endometriosis-associated pain and 683 with heavy menstrual bleeding",
      "associated with uterine fibroids (phase III)"
    ),
    dose_range = paste(
      "Oral elagolix single doses of 150-300 mg and multiple doses from",
      "150 mg QD to 400 mg BID for 9 days to 12 months; in the uterine-fibroid",
      "studies 300 mg BID alone or with estradiol 1 mg / norethindrone",
      "acetate 0.5 mg QD"
    ),
    regions = "Multinational (phase I studies at US sites)",
    oatp1b1_status = c(
      extensive = 57.9,
      intermediate = 15.5,
      poor = 1.48,
      missing = 25.1
    ),
    co_medication = "Estradiol/norethindrone acetate 1/0.5 mg QD add-back therapy in part of the uterine-fibroid cohort (not a significant covariate)",
    notes = paste(
      "Baseline demographics from Beck 2022 Table 2; studies and regimens",
      "from Table 1. 17,915 elagolix plasma concentrations were analysed",
      "(4511 phase I, 8685 endometriosis, 4719 uterine fibroids)."
    )
  )

  # Beck 2022 estimated a separate combined residual-error model for the
  # phase I and the phase III studies. nlmixr2 takes one residual-SD symbol
  # per term, so the stratum-specific SDs are separate ini() parameters
  # combined inside model() with the STUDY_PHASE3 indicator (the
  # Rich_2026_momelotinib.R pattern), declared here so checkModelConventions()
  # does not read the stratum suffixes as deviant residual-error names.
  paper_specific_residual_sds <- c(
    "propSdPh1",
    "addSdPh1",
    "propSdPh3",
    "addSdPh3"
  )

  ini({
    # Structural parameters (Table 3, 'Final pharmacokinetic model').
    # All are apparent (per unit oral bioavailability).
    lcl <- log(125); label("Apparent clearance CL/F (L/h)") # Table 3 'CL/F (L/h)' = 125
    lvc <- log(279); label("Apparent central volume Vc/F at 76 kg (L)") # Table 3 'Vc/F (L)' = 279
    lka <- log(2.46); label("First-order absorption rate constant ka (1/h)") # Table 3 'KA' = 2.46 (printed unit 'L/h' is a typo for 1/h)
    lq <- log(5.63); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (L/h)' = 5.63
    lvp <- log(51.7); label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp/F (L)' = 51.7
    ltlag <- log(0.207); label("Absorption lag time (h)") # Table 3 'Lag time (h)' = 0.207
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 in extensive transporters (unitless)") # Table 3 'F1' = 1.00 (fix)

    # Covariate effects
    e_wt_vc <- 0.160; label("Power exponent of body weight on Vc/F (unitless)") # Table 3 'Body weight on Vc/F' = 0.160; footnote c
    e_slco1b1_521_het_fdepot <- 0.421; label("Fractional increase in F1 for OATP1B1 intermediate transporters (unitless)") # Table 3 'Intermediate transporter on F1' = 0.421
    e_slco1b1_521_hom_fdepot <- 0.963; label("Fractional increase in F1 for OATP1B1 poor transporters (unitless)") # Table 3 'Poor transporter on F1' = 0.963
    e_slco1b1_521_missing_fdepot <- 0.101; label("Fractional increase in F1 for missing OATP1B1 genotype (unitless)") # Table 3 'Missing transporter on F1' = 0.101

    # IIV (variances, Table 3; %CV = 100*sqrt(exp(omega^2)-1) per footnote b).
    # The base model estimated a CL/F-Vc/F covariance, but Table 3 does not
    # report it; the covariance is set to zero here (vignette Errata).
    etalcl ~ 0.198 # Table 3 'IIV on CL/F' = 0.198 (46.8% CV)
    etalvc ~ 0.208 # Table 3 'IIV on Vc/F' = 0.208 (48.1% CV)

    # Residual error: Table 3 prints NONMEM SIGMA variances (ESM Equation 2,
    # eps ~ N(0, sigma^2)); nlmixr2 takes SDs, so each is square-rooted.
    propSdPh1 <- sqrt(0.145); label("Proportional residual SD, phase I studies (fraction)") # Table 3 'Proportional error (phase I studies)' = 0.145 (variance)
    addSdPh1 <- sqrt(5.26e-05); label("Additive residual SD, phase I studies (ng/mL)") # Table 3 'Additive error (phase I studies)' = 5.26e-05 (variance)
    propSdPh3 <- sqrt(0.284); label("Proportional residual SD, phase III studies (fraction)") # Table 3 'Proportional error (phase III studies)' = 0.284 (variance)
    addSdPh3 <- sqrt(0.266); label("Additive residual SD, phase III studies (ng/mL)") # Table 3 'Additive error (phase III studies)' = 0.266 (variance)
  })

  model({
    # Individual parameters (ESM Equation 3; Table 3 footnotes c and d)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (WT / 76)^e_wt_vc
    ka <- exp(lka)
    q <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag)

    # Relative bioavailability: F1 = 1 + theta_IT*IT + theta_PT*PT + theta_miss*missing
    # (Table 3 footnote d: 1.00 ET, 1.42 IT, 1.96 PT, 1.10 missing)
    fdepot <- exp(lfdepot) *
      (1 +
        e_slco1b1_521_het_fdepot * SNP_SLCO1B1_RS4149056_HET +
        e_slco1b1_521_hom_fdepot * SNP_SLCO1B1_RS4149056_HOM +
        e_slco1b1_521_missing_fdepot * SNP_SLCO1B1_RS4149056_MISSING)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x1000 gives ng/mL
    Cc <- 1000 * central / vc

    # Phase-stratified combined error (ESM Equation 2)
    propSd <- propSdPh3 * STUDY_PHASE3 + propSdPh1 * (1 - STUDY_PHASE3)
    addSd <- addSdPh3 * STUDY_PHASE3 + addSdPh1 * (1 - STUDY_PHASE3)
    Cc ~ add(addSd) + prop(propSd)
  })
}
