Karakitsios_2024_gevokizumab_mbma <- function() {
  description <- paste0(
    "mPBPK / MBMA. Five-dosage-group meta-analysis of gevokizumab ",
    "(anti-interleukin-1beta humanized IgG2) in adults with type 2 ",
    "diabetes, fitted to AGGREGATE (mean and SD versus time) plasma ",
    "concentration data from five dose arms simultaneously. The structural ",
    "model is the Cao 2013 second-generation minimal PBPK model; ",
    "Karakitsios and Dokoumetzidis estimated its three drug-specific ",
    "parameters (tight- and leaky-tissue vascular reflection coefficients ",
    "and plasma clearance) hierarchically, with an exponential ",
    "INTER-GROUP variability (IGV) term on each, following a log-Student's ",
    "t distribution with 5 degrees of freedom to down-weight outlying dose ",
    "arms. The random effects in this file are therefore DOSAGE-GROUP-level ",
    "(eta_study_*), not between-subject: the model simulates group-mean ",
    "concentration-time profiles and is NOT suitable for individual-subject ",
    "simulation. Between-subject variability was estimated separately for ",
    "each of the five arms (a semi-hierarchical design with no distribution ",
    "assumed across arms), so no single population IIV exists; the five ",
    "per-arm pairs are tabulated in population$notes. For an ",
    "individual-level version of the same structure, see the single ",
    "dosage-group fit modellib('Karakitsios_2024_gevokizumab')."
  )

  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "A Meta-Analysis Methodology in Stan to Estimate Population",
    "Pharmacokinetic Parameters from Multiple Aggregate Concentration-Time",
    "Datasets: Application to Gevokizumab mPBPK Model.",
    "Pharmaceutics. 2024 Aug 27;16(9):1129.",
    "doi:10.3390/pharmaceutics16091129.",
    "Structural mPBPK equations and physiologic constants are those of",
    "Cao Y, Balthasar JP, Jusko WJ. J Pharmacokinet Pharmacodyn.",
    "2013;40(5):597-607, reproduced in this paper's Supplementary Material",
    "equations (S1)-(S8); see modellib('Cao_2013_gevokizumab').",
    "Aggregate plasma data were digitised from",
    "Cavelti-Weder C et al. Diabetes Care. 2012;35(8):1654-1662",
    "(PMID 22699287).",
    sep = " "
  )
  vignette <- "Karakitsios_2024_gevokizumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. States hold antibody AMOUNTS (mg); the supplement
  # writes equations (S1)-(S4) in concentration form, which is the same
  # system divided through by each compartment volume. verified = TRUE:
  # analyte and specimen confirmed against the paper's Supplementary
  # Material (definitions following equation (S4)) and Figure S1.
  compartmentData <- list(
    plasma = list(analyte = "gevokizumab", units = "mg", specimen = "plasma", verified = TRUE),
    tight  = list(analyte = "gevokizumab", units = "mg", specimen = "tissue", verified = TRUE),
    leaky  = list(analyte = "gevokizumab", units = "mg", specimen = "tissue", verified = TRUE),
    lymph  = list(analyte = "gevokizumab", units = "mg", specimen = "lymph", verified = TRUE)
  )

  covariateData <- list()

  # Body weight underlies both the physiologic constants and the mg/kg-to-mg
  # dose conversion, but it never enters model() as a covariate: the digitised
  # aggregate data carry no individual weights. Dose is not a covariate either
  # -- it rides on the event table's amt (see population$dose_range).
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Not a data column. The physiologic parameters (L = 2.9 L/day, ",
        "VISF = 15.6 L, Vlymph = 5.2 L, Vplasma = 2.6 L) are quoted for a ",
        "70 kg person (Supplementary Material, text following equation ",
        "(S8)), and the five arms' mg/kg doses convert to mg on a 70 kg ",
        "basis (Sections 2.1 and 2.2). No individual weights were available ",
        "from the digitised aggregate data."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    n_groups       = 5L,
    age_range      = "adults (per Cavelti-Weder 2012 source study)",
    weight_range   = "70 kg reference body weight (Supplementary Material, text following equation (S8))",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Type 2 diabetes mellitus",
    dose_range     = "0.01, 0.03, 0.1, 0.3 and 1.0 mg/kg single IV dose (0.7, 2.1, 7, 21 and 70 mg for a 70 kg person; Section 2.1)",
    regions        = NA,
    data_origin    = paste0(
      "AGGREGATE data only. Mean plasma concentrations and their SDs versus ",
      "time were digitised with Digitizer software from the published graph ",
      "of Cavelti-Weder C et al. Diabetes Care 2012;35:1654-1662 (Section ",
      "2.1). No individual concentration records were available. The sixth ",
      "published arm (3.0 mg/kg, n = 5) was EXCLUDED by the authors because ",
      "5 patients were considered too few to give accurate means and SDs ",
      "(Section 2.4)."
    ),
    notes          = paste0(
      "Karakitsios 2024 Table 4 (posterior means), five dosage groups of 10 ",
      "patients each (50 patients total). Sampling times 4, 8, 24, 48, 72, ",
      "96, 168, 216, 264, 336, 504, 672, 1008 and 1344 h after a single IV ",
      "dose (Section 2.1). Fitted in RStan 2.19.2, 4 chains x 2000 samples ",
      "with 2000 warm-up iterations, run time 41,173.5 s; all Rhat < 1.01 ",
      "and bulk/tail ESS > 400 (Sections 2.4 and 3.2.2). Source code: ",
      "https://github.com/PMXathens/Gevokizumab (FIVE_DOSES.RData, ",
      "FIVE_DOSES_FIT.R, FIVE_DOSES_FIT.stan). ",
      "PER-ARM BETWEEN-SUBJECT VARIABILITY (Table 4; omega = SD on the log ",
      "scale, arms 1-5 = 0.01, 0.03, 0.1, 0.3, 1.0 mg/kg): omega_CLp = ",
      "0.1813, 0.1138, 0.0798, 0.0374, 0.2002; omega_V = 0.1676, 0.3213, ",
      "0.0733, 0.0960, 0.1000. These are NOT encoded in ini() because the ",
      "authors assumed no distribution across arms, so there is no ",
      "population IIV to carry; override etalcl / etalvc on ",
      "modellib('Karakitsios_2024_gevokizumab') to simulate a chosen arm. ",
      "PER-ARM PARAMETER ESTIMATES (Table 4, log scale, arms 1-5): ",
      "log_CLp = -5.0572, -4.9818, -5.0418, -5.1067, -5.1280; log_rc1 = ",
      "-0.0588, -0.0528, -0.0441, -0.0489, -0.0513; log_rc2 = -0.4421, ",
      "-0.2947, -0.2639, -0.1878, -0.2086. ",
      "The paper reports TWO residual error terms because the fitted ",
      "quantities are aggregate summaries: sigma_1 = 0.0734 on the mean ",
      "concentrations (encoded here as expSd) and sigma_2 = 0.2706 on the ",
      "SDs of the concentrations (no individual-level analogue; see ",
      "vignette Assumptions and deviations)."
    )
  )

  ini({
    # Mean-group (population-level) drug-specific parameters estimated by
    # Karakitsios 2024 (Table 4, posterior means). Physiologic restrictions
    # (S5)-(S6) require both reflection coefficients to lie in (0, 1); a
    # beta(1, 1) prior enforced these bounds on the mean-group values during
    # estimation (Section 2.4).
    sigma_tight <- 0.9504; label("Mean-group vascular reflection coefficient for tight tissues (rc1_mean, unitless)")  # Karakitsios 2024 Table 4: rc1_mean = 0.9504 (95% CrI 0.895-0.990)
    sigma_leaky <- 0.7674; label("Mean-group vascular reflection coefficient for leaky tissues (rc2_mean, unitless)")  # Karakitsios 2024 Table 4: rc2_mean = 0.7674 (95% CrI 0.647-0.896)
    lcl         <- log(0.0064); label("Mean-group plasma clearance (CLp_mean, L/h)")                                   # Karakitsios 2024 Table 4: CLp_mean = 0.0064 L/h (95% CrI 0.006-0.007)

    # Plasma volume is a fixed physiologic constant for a 70 kg person, not
    # an estimated parameter.
    lvc <- fixed(log(2.6)); label("Plasma volume (Vplasma, L; physiologic constant for a 70 kg person)")  # Supplementary Material, text following equation (S8): Vplasma = 2.6 L (Cao et al. 2013)

    # INTER-GROUP variability (IGV), not between-subject variability. Each
    # term is the SD of an exponential (log-scale) random effect across the
    # five dosage groups; the group-level parameters were assumed to follow
    # a log-Student's t distribution with df fixed at 5 to down-weight
    # outlying arms (Section 2.4). nlmixr2 draws these etas from a normal
    # rather than a t distribution -- see vignette Assumptions and
    # deviations. Variance entered here is the reported SD squared.
    eta_study_lcl         ~ 0.01572516  # Karakitsios 2024 Table 4: gamma_CLp = 0.1254 -> 0.1254^2 = 0.01572516
    eta_study_sigma_tight ~ 0.00114244  # Karakitsios 2024 Table 4: gamma_rc1 = 0.0338 -> 0.0338^2 = 0.00114244
    eta_study_sigma_leaky ~ 0.03276100  # Karakitsios 2024 Table 4: gamma_rc2 = 0.1810 -> 0.1810^2 = 0.032761

    # Residual error. The paper's sigma_1 is an EXPONENTIAL error between the
    # model-predicted and observed MEAN concentrations (an aggregate-data
    # residual), which is why it is mapped onto a lognormal observation
    # error here. See population$notes and the vignette Assumptions and
    # deviations for why this is not an individual-level RUV estimate.
    expSd <- 0.0734; label("Exponential residual error on the aggregate mean concentrations (sigma_1, fraction)")  # Karakitsios 2024 Table 4: sigma_1 = 0.0734 (95% CrI 0.061-0.089)
  })

  model({
    # Fixed physiologic constants for a 70 kg person. Karakitsios 2024
    # Supplementary Material, text following equation (S8), quoting Cao et
    # al. 2013: L = 2.9 L/day, VISF = 15.6 L, Vlymph = 5.2 L, Vplasma = 2.6
    # L. Lymph flow is converted to L/h because this model's time unit is
    # hours (the paper reports CLp in L/h and samples in h).
    sigmal    <- 0.2       # Karakitsios 2024 Supplementary Material, text following equation (S4): lymphatic capillary reflection coefficient rcL assumed 0.2
    kp        <- 0.8       # Karakitsios 2024 Supplementary Material, text following equation (S8): Kp (available ISF fraction) = 0.8 for gevokizumab
    lymphflow <- 2.9 / 24  # Karakitsios 2024 Supplementary Material: L = 2.9 L/day -> L/h

    vplasma <- exp(lvc)
    visf    <- 15.6  # Karakitsios 2024 Supplementary Material, text following equation (S8): VISF = 15.6 L
    vlymph  <- 5.2   # Karakitsios 2024 Supplementary Material, text following equation (S8): Vlymph = 5.2 L

    vtight <- 0.65 * visf * kp  # Karakitsios 2024 Supplementary Material equation (S7)
    vleaky <- 0.35 * visf * kp  # Karakitsios 2024 Supplementary Material equation (S8)
    l1     <- 0.33 * lymphflow  # Karakitsios 2024 Supplementary Material, text following equation (S4): L1 = 0.33 * L
    l2     <- 0.67 * lymphflow  # Karakitsios 2024 Supplementary Material, text following equation (S4): L2 = 0.67 * L

    # Dosage-group realisations of the three drug-specific parameters. The
    # exponential IGV is applied on the log scale, reproducing the paper's
    # log_CLp[i], log_rc1[i] and log_rc2[i] group parameters (Table 4).
    cl           <- exp(lcl + eta_study_lcl)
    sigma_tightI <- sigma_tight * exp(eta_study_sigma_tight)
    sigma_leakyI <- sigma_leaky * exp(eta_study_sigma_leaky)

    cp     <- plasma / vplasma
    ctight <- tight / vtight
    cleaky <- leaky / vleaky
    clymph <- lymph / vlymph

    # Karakitsios 2024 Supplementary Material equations (S1)-(S4), written
    # here in amounts (each concentration-form equation multiplied through
    # by its compartment volume).
    d/dt(plasma) <- clymph * lymphflow -
                    cp * l1 * (1 - sigma_tightI) -
                    cp * l2 * (1 - sigma_leakyI) -
                    cl * cp
    d/dt(tight)  <- l1 * (1 - sigma_tightI) * cp -
                    l1 * (1 - sigmal) * ctight
    d/dt(leaky)  <- l2 * (1 - sigma_leakyI) * cp -
                    l2 * (1 - sigmal) * cleaky
    d/dt(lymph)  <- l1 * (1 - sigmal) * ctight +
                    l2 * (1 - sigmal) * cleaky -
                    clymph * lymphflow

    Cc <- cp
    Cc ~ lnorm(expSd)
  })
}
