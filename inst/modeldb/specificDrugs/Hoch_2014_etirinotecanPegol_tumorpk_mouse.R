Hoch_2014_etirinotecanPegol_tumorpk_mouse <- function() {
  description <- paste(
    "Preclinical (mouse with HT29 colorectal-cancer xenograft). Descriptive",
    "two-compartment model of etirinotecan pegol (NKTR-102) concentration in",
    "tumor after IV bolus doses, fitted by the authors directly to mean tumor",
    "concentrations (no plasma link). The paper reports no parameter values;",
    "the biexponential disposition was back-solved by the maintainers from",
    "the model-predicted tumor curve of Hoch 2014 Figure 2d."
  )
  reference <- paste(
    "Hoch U, Staschen C-M, Johnson RK, Eldon MA. Nonclinical",
    "pharmacokinetics and activity of etirinotecan pegol (NKTR-102), a",
    "long-acting topoisomerase 1 inhibitor, in multiple cancer models.",
    "Cancer Chemother Pharmacol. 2014;74(6):1125-1137.",
    "doi:10.1007/s00280-014-2577-7.",
    sep = " "
  )
  vignette <- "Hoch_2014_etirinotecanPegol"
  paper_specific_compartments <- c("tumor_ep")

  units <- list(
    time = "day",
    dosing = "mg/kg (irinotecan equivalents)",
    concentration = "ug/g (etirinotecan pegol in tumor homogenate, as irinotecan equivalents)"
  )

  compartmentData <- list(
    tumor_ep = list(
      analyte = "etirinotecan pegol",
      units = "mg/kg (irinotecan equivalents; apparent amount)",
      specimen = "tumor",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "etirinotecan pegol",
      units = "mg/kg (irinotecan equivalents; apparent amount)",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "mouse (female athymic nude, Ncr:Nu and Nu:Nu, with subcutaneous HT29 human colon-carcinoma xenograft)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "6-10 weeks (xenograft studies)",
    weight_range = "14-28 g (xenograft studies)",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state = "subcutaneous HT29 human colon-carcinoma xenograft; median tumor volume 100-172 mm^3 at the first dose (13 days after implant)",
    dose_range = "40 mg/kg irinotecan-equivalent IV bolus of etirinotecan pegol on days 0, 4 and 8 (q4d x 3)",
    regions = "preclinical (Nektar Therapeutics, USA)",
    notes = paste(
      "Groups of 4 mice per sampling time; tumors were sampled before dosing and at 0.5, 1, 4, 8",
      "and 12 h and 1, 2, 3, 4, 5, 10, 12, 15, 20, 30, 40, 50 and 60 days after the start of dosing.",
      "Mean concentrations were fitted in Berkeley Madonna to predict a continuous profile for",
      "noncompartmental analysis (Hoch 2014 Methods, 'Plasma and tumor pharmacokinetics...'); the",
      "fit has no between-animal variability and no residual-error estimate. The total number of",
      "animals is not reported."
    )
  )

  ini({
    # Not reported in the paper (Results: 'A three- and two-compartment PK model
    # well described the concentrations of etirinotecan pegol in plasma and tumor,
    # respectively'). Back-solved by the maintainers by least squares (log scale)
    # to the solid model-predicted tumor curve of Figure 2d, digitised at 86
    # points (log-RMSE 3.5%). The fitted biexponential per 40 mg/kg dose is
    # 19.3 * exp(-1.716 t) + 14.34 * exp(-0.0465 t) ug/g, converted to the
    # rate constants below for an IV bolus into the sampled tumor pool.
    lkel <- log(0.105); label("Elimination rate constant from the sampled tumor pool (1/day)") # back-solved from Figure 2d (maintainer-estimated)
    lk12 <- log(0.899); label("Transfer rate constant, sampled tumor pool to second compartment (1/day)") # back-solved from Figure 2d (maintainer-estimated)
    lk21 <- log(0.758); label("Transfer rate constant, second compartment to sampled tumor pool (1/day)") # back-solved from Figure 2d (maintainer-estimated)
    lv_tumor_ep <- log(1.19); label("Apparent dose-to-tumor-concentration ratio for etirinotecan pegol ((mg/kg)/(ug/g))") # back-solved from Figure 2d (maintainer-estimated)

    # Not reported: the fit was to group-mean concentrations. Fixed at zero.
    propSd <- fixed(0); label("Proportional residual error on tumor etirinotecan pegol (fraction; not reported)") # not reported in Hoch 2014
  })

  model({
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    v_tumor_ep <- exp(lv_tumor_ep)

    # Dose etirinotecan pegol (mg/kg irinotecan equivalents) into tumor_ep.
    d/dt(tumor_ep) <- -(kel + k12) * tumor_ep + k21 * peripheral1
    d/dt(peripheral1) <- k12 * tumor_ep - k21 * peripheral1

    Ctumor <- tumor_ep / v_tumor_ep
    Ctumor ~ prop(propSd)
  })
}
