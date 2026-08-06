Karakitsios_2024_gevokizumab <- function() {
  description <- paste0(
    "mPBPK. Second-generation minimal physiologically-based population PK ",
    "model for gevokizumab (anti-interleukin-1beta humanized IgG2) in adults ",
    "with type 2 diabetes, with between-subject variability estimated from ",
    "AGGREGATE (mean and SD versus time) plasma concentration data rather ",
    "than individual records. Karakitsios and Dokoumetzidis re-estimated the ",
    "three drug-specific parameters of the Cao 2013 second-generation mPBPK ",
    "structure (the tight- and leaky-tissue vascular reflection coefficients ",
    "and plasma clearance) plus two IIV terms, using a Bayesian ",
    "reconstruct-the-likelihood-from-aggregate-data method implemented in ",
    "RStan: at each MCMC iteration a Latin-hypercube virtual population is ",
    "solved, its mean and SD versus time are computed, and those summaries ",
    "are fitted to the digitised published mean and SD profiles. This file ",
    "holds the SINGLE dosage-group fit (7 mg IV, n = 10; Table 3). A single ",
    "lognormal body-size random effect (etalvc) is shared by all four ",
    "physiologic volumes, per the authors' assumption that a larger patient ",
    "has proportionally larger plasma, tight-tissue ISF, leaky-tissue ISF ",
    "and lymph volumes. The five-dosage-group meta-analysis from the same ",
    "paper is a separate model; see ",
    "modellib('Karakitsios_2024_gevokizumab_mbma')."
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

  # Body weight is the implicit basis of the physiologic constants (all are
  # quoted for a 70 kg person) and of the authors' shared volume random
  # effect, but it is never a data column in this analysis: the aggregate
  # data carry no individual weights, so no weight term appears in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Not a data column. The physiologic parameters (L = 2.9 L/day, ",
        "VISF = 15.6 L, Vlymph = 5.2 L, Vplasma = 2.6 L) are quoted for a ",
        "70 kg person (Supplementary Material, text following equation ",
        "(S8)), and the 7 mg dose is the study's 0.1 mg/kg arm scaled to ",
        "70 kg (Section 2.2). The authors' shared volume IIV (etalvc) is ",
        "their proxy for unmeasured body-size differences, so weight is ",
        "represented as random variability rather than as a covariate."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "adults (per Cavelti-Weder 2012 source study)",
    weight_range   = "70 kg reference body weight (Supplementary Material, text following equation (S8))",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Type 2 diabetes mellitus",
    dose_range     = "7 mg single IV dose (the 0.1 mg/kg arm scaled to 70 kg; Sections 2.2 and 3.2.1)",
    regions        = NA,
    data_origin    = paste0(
      "AGGREGATE data only. Mean plasma concentrations and their SDs versus ",
      "time were digitised with Digitizer software from the published graph ",
      "of Cavelti-Weder C et al. Diabetes Care 2012;35:1654-1662 (Section ",
      "2.1). No individual concentration records were available."
    ),
    notes          = paste0(
      "Karakitsios 2024 Table 3 (posterior means), 7 mg dosage group, n = 10 ",
      "patients. Sampling times 4, 8, 24, 48, 72, 96, 168, 216, 264, 336, ",
      "504, 672, 1008 and 1344 h after a single IV dose (Section 2.1). ",
      "Fitted in RStan 2.19.2, 4 chains x 1000 samples with 1000 warm-up ",
      "iterations; all Rhat < 1.01 and bulk/tail ESS > 400 (Section 3.2.1). ",
      "Source code: https://github.com/PMXathens/Gevokizumab ",
      "(ONE_DOSE.RData, '7 mg dosage-group.R'). The paper reports TWO ",
      "residual error terms because the fitted quantities are aggregate ",
      "summaries: sigma_1 = 0.0758 on the mean concentrations (encoded here ",
      "as expSd) and sigma_2 = 0.2316 on the SDs of the concentrations ",
      "(no individual-level analogue; see vignette Assumptions and ",
      "deviations). The authors deliberately applied no residual error to ",
      "the individual virtual patients' concentrations (Section 3.3)."
    )
  )

  ini({
    # Drug-specific parameters estimated by Karakitsios 2024 (Table 3,
    # posterior means). Physiologic restrictions (S5)-(S6) require both
    # reflection coefficients to lie in (0, 1); a beta(1, 1) prior enforced
    # these bounds during estimation (Section 2.2).
    sigma_tight <- 0.9584; label("Vascular reflection coefficient for tight tissues (rc1, unitless)")  # Karakitsios 2024 Table 3: rc1_mean = 0.9584 (95% CrI 0.877-0.999)
    sigma_leaky <- 0.7645; label("Vascular reflection coefficient for leaky tissues (rc2, unitless)")  # Karakitsios 2024 Table 3: rc2_mean = 0.7645 (95% CrI 0.709-0.830)
    lcl         <- log(0.0065); label("Plasma clearance (CLp, L/h)")                                   # Karakitsios 2024 Table 3: CLp_mean = 0.0065 L/h (95% CrI 0.006-0.007)

    # Plasma volume is a fixed physiologic constant for a 70 kg person, not
    # an estimated parameter. It is carried in ini() rather than as a
    # model() local because it is the anchor for the shared volume random
    # effect etalvc below.
    lvc <- fixed(log(2.6)); label("Plasma volume (Vplasma, L; physiologic constant for a 70 kg person)")   # Supplementary Material, text following equation (S8): Vplasma = 2.6 L (Cao et al. 2013)

    # IIV. Both terms are the SD of a lognormal distribution (Section 2.2),
    # so omega = SD on the log scale and the variance entered here is that
    # SD squared. IIV was deliberately NOT applied to the two reflection
    # coefficients because doing so gave simulated SDs far from the observed
    # ones (Section 2.2).
    etalcl ~ 0.00600625  # Karakitsios 2024 Table 3: omega_CLp = 0.0775 -> 0.0775^2 = 0.00600625
    etalvc ~ 0.00488601  # Karakitsios 2024 Table 3: omega_V   = 0.0699 -> 0.0699^2 = 0.00488601

    # Residual error. The paper's sigma_1 is an EXPONENTIAL error between the
    # model-predicted and observed MEAN concentrations (an aggregate-data
    # residual), which is why it is mapped onto a lognormal observation
    # error here. See the population$notes and vignette Assumptions and
    # deviations for why this is not an individual-level RUV estimate.
    expSd <- 0.0758; label("Exponential residual error on the aggregate mean concentrations (sigma_1, fraction)")  # Karakitsios 2024 Table 3: sigma_1 = 0.0758 (95% CrI 0.049-0.123)
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

    # Shared body-size random effect. Section 2.2: "the same IIV term for
    # volume ... is applied to all 4 volumes of the mPBPK model ... a bigger
    # patient will obviously have more plasma and tight and leaky tissues as
    # well as larger lymph volumes than a smaller one."
    volscale <- exp(etalvc)

    vplasma <- exp(lvc) * volscale
    visf    <- 15.6 * volscale  # Karakitsios 2024 Supplementary Material, text following equation (S8): VISF = 15.6 L
    vlymph  <- 5.2 * volscale   # Karakitsios 2024 Supplementary Material, text following equation (S8): Vlymph = 5.2 L

    vtight <- 0.65 * visf * kp  # Karakitsios 2024 Supplementary Material equation (S7)
    vleaky <- 0.35 * visf * kp  # Karakitsios 2024 Supplementary Material equation (S8)
    l1     <- 0.33 * lymphflow  # Karakitsios 2024 Supplementary Material, text following equation (S4): L1 = 0.33 * L
    l2     <- 0.67 * lymphflow  # Karakitsios 2024 Supplementary Material, text following equation (S4): L2 = 0.67 * L

    cl <- exp(lcl + etalcl)

    cp     <- plasma / vplasma
    ctight <- tight / vtight
    cleaky <- leaky / vleaky
    clymph <- lymph / vlymph

    # Karakitsios 2024 Supplementary Material equations (S1)-(S4), written
    # here in amounts (each concentration-form equation multiplied through
    # by its compartment volume).
    d/dt(plasma) <- clymph * lymphflow -
                    cp * l1 * (1 - sigma_tight) -
                    cp * l2 * (1 - sigma_leaky) -
                    cl * cp
    d/dt(tight)  <- l1 * (1 - sigma_tight) * cp -
                    l1 * (1 - sigmal) * ctight
    d/dt(leaky)  <- l2 * (1 - sigma_leaky) * cp -
                    l2 * (1 - sigmal) * cleaky
    d/dt(lymph)  <- l1 * (1 - sigmal) * ctight +
                    l2 * (1 - sigmal) * cleaky -
                    clymph * lymphflow

    Cc <- cp
    Cc ~ lnorm(expSd)
  })
}
