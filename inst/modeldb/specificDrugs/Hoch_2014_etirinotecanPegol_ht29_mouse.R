Hoch_2014_etirinotecanPegol_ht29_mouse <- function() {
  description <- paste(
    "Preclinical (mouse with HT29 colorectal-cancer xenograft). Tumor PK/PD",
    "model for etirinotecan pegol (NKTR-102) and conventional irinotecan:",
    "intratumoral SN-38 drives an inhibitory Emax effect (Emax fixed at 1)",
    "on Gompertz growth of the relative tumor weight. The Gompertz and Emax",
    "parameters are the published estimates; the tumor SN-38 input functions",
    "(first-order release after etirinotecan pegol, bolus-type input after",
    "irinotecan) were not reported and were back-solved by the maintainers",
    "from the model-predicted SN-38 curves of Hoch 2014 Figure 3."
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
  paper_specific_compartments <- c("tumor_sn38", "tumor_sn38_irinotecan")

  units <- list(
    time = "day",
    dosing = "mg/kg (irinotecan equivalents for both etirinotecan pegol and irinotecan)",
    concentration = "ug/g (SN-38 in tumor homogenate)"
  )

  compartmentData <- list(
    depot = list(
      analyte = "etirinotecan pegol",
      units = "mg/kg (irinotecan equivalents)",
      specimen = "administration site",
      verified = TRUE
    ),
    tumor_sn38 = list(
      analyte = "SN-38 (formed from etirinotecan pegol)",
      units = "mg/kg (irinotecan equivalents; apparent amount)",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_sn38_irinotecan = list(
      analyte = "SN-38 (formed from conventional irinotecan)",
      units = "mg/kg (irinotecan equivalents; apparent amount)",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_vol = list(
      analyte = "none",
      units = "relative tumor weight (unitless; 1 at the first dose)",
      specimen = "not applicable",
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
    disease_state = "subcutaneous HT29 human colon-carcinoma xenograft (30-40 mg tumor fragments near the right axillary area); median tumor volume 100-172 mm^3 at the first dose (13 days after implant)",
    dose_range = "40 mg/kg irinotecan-equivalent IV bolus of etirinotecan pegol or of irinotecan on days 0, 4 and 8 (q4d x 3); vehicle control",
    regions = "preclinical (Nektar Therapeutics, USA)",
    notes = paste(
      "Groups of 4 mice per sampling time; blood and tumor were sampled before dosing and at",
      "0.5, 1, 4, 8 and 12 h and 1, 2, 3, 4, 5, 10, 12, 15, 20, 30, 40, 50 and 60 days after the start",
      "of dosing (Hoch 2014 Methods, 'Pharmacokinetic and pharmacodynamic study in tumor-bearing mice').",
      "The PD model was fitted by nonlinear regression to the MEAN relative tumor volume of the",
      "control, irinotecan and etirinotecan pegol groups simultaneously (Berkeley Madonna,",
      "Runge-Kutta 4), so the model has no between-animal variability and no residual-error",
      "estimate. The total number of animals is not reported."
    )
  )

  ini({
    # --- Tumor growth and drug effect: published estimates --------------------
    # Hoch 2014 Results, 'Etirinotecan pegol accumulation in tumor leads to high
    # tumor SN38 exposure...': 'Estimates for PK/PD parameters kgr, Limit, and
    # EC50 were as follows: 0.033 (day-1), 52, and 0.20 ug/g for HT29'.
    lkgrow <- log(0.033); label("Intrinsic tumor growth rate kgr (1/day)") # Results: kgr = 0.033 day-1 (HT29)
    ltsmax <- log(52); label("Gompertz limit of tumor growth Limit (relative tumor weight)") # Results: Limit = 52 (HT29)
    lec50 <- log(0.20); label("Tumor SN-38 concentration for half-maximal growth inhibition EC50 (ug/g)") # Results: EC50 = 0.20 ug/g (HT29)
    emax <- fixed(1); label("Maximal fractional inhibition of tumor growth Emax (unitless)") # Methods, 'Tumor PK/PD analysis': 'Emax was fixed at a value of 1'

    # --- Tumor SN-38 after etirinotecan pegol: back-solved -------------------
    # Not reported in the paper. Back-solved by the maintainers by least squares
    # (log scale) to the dashed model-predicted tumor SN-38 curve of Figure 3c,
    # digitised at 74 points (log-RMSE 3.5%). First-order release from a dose
    # depot into tumor SN-38; the slower of the two rate constants is assigned to
    # the release (flip-flop), consistent with the slow loss of etirinotecan pegol
    # from tumor (abstract: t1/2 = 17 days).
    lka <- log(0.0457); label("First-order release of SN-38 into tumor after etirinotecan pegol (1/day)") # back-solved from Figure 3c (maintainer-estimated)
    lkel_sn38 <- log(0.235); label("Elimination rate constant of tumor SN-38 after etirinotecan pegol (1/day)") # back-solved from Figure 3c (maintainer-estimated)
    lv_tumor_sn38 <- log(15.3); label("Apparent dose-to-tumor-concentration ratio for SN-38 after etirinotecan pegol ((mg/kg)/(ug/g))") # back-solved from Figure 3c (maintainer-estimated)

    # --- Tumor SN-38 after conventional irinotecan: back-solved --------------
    # Not reported in the paper. Bolus-type one-compartment input (Results: 'One-
    # and two-compartment PK models ... for conventional irinotecan and
    # etirinotecan pegol, respectively'). The predicted spike height in Figure 3a
    # is 140 ng/g at each dose; the decay is narrower than one pixel, so the
    # elimination rate constant is set from the printed tumor SN-38 AUC of the
    # irinotecan arm (Results: '0.100 ug day/g' for three doses):
    # k = 3 * 0.140 / 0.100 = 4.2 1/day.
    lkel_sn38_irinotecan <- log(4.2); label("Elimination rate constant of tumor SN-38 after irinotecan (1/day)") # derived: 3 doses x 0.140 ug/g (Figure 3a spike) / 0.100 ug*day/g (Results text AUC)
    lv_tumor_sn38_irinotecan <- log(286); label("Apparent dose-to-tumor-concentration ratio for SN-38 after irinotecan ((mg/kg)/(ug/g))") # 40 mg/kg / 0.140 ug/g; spike height digitised from Figure 3a (maintainer-estimated)

    # --- Residual error -------------------------------------------------------
    # Not reported: the PD model was fitted to group-mean relative tumor volume by
    # nonlinear regression and the paper gives no error model. Fixed at zero.
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on relative tumor weight (fraction; not reported)") # not reported in Hoch 2014
    propSd_Ctumor_sn38 <- fixed(0); label("Proportional residual error on tumor SN-38 concentration (fraction; not reported)") # not reported in Hoch 2014
  })

  model({
    kgrow <- exp(lkgrow)
    tsmax <- exp(ltsmax)
    ec50 <- exp(lec50)

    ka <- exp(lka)
    kel_sn38 <- exp(lkel_sn38)
    v_tumor_sn38 <- exp(lv_tumor_sn38)
    kel_sn38_irinotecan <- exp(lkel_sn38_irinotecan)
    v_tumor_sn38_irinotecan <- exp(lv_tumor_sn38_irinotecan)

    # Tumor SN-38 input functions (Online Resource 1 'Tumor Exposure Model').
    # Dose etirinotecan pegol into depot and irinotecan into
    # tumor_sn38_irinotecan, both in mg/kg irinotecan equivalents.
    d/dt(depot) <- -ka * depot
    d/dt(tumor_sn38) <- ka * depot - kel_sn38 * tumor_sn38
    d/dt(tumor_sn38_irinotecan) <- -kel_sn38_irinotecan * tumor_sn38_irinotecan

    # Total tumor SN-38 (ug/g). Only one drug is given per treatment arm.
    Ctumor_sn38 <- tumor_sn38 / v_tumor_sn38 + tumor_sn38_irinotecan / v_tumor_sn38_irinotecan

    # Online Resource 1: Effect = 1 - Emax * Ctu(t) / (EC50 + Ctu(t)), applied
    # to the intrinsic Gompertz growth model {kgr, Limit}.
    effect <- 1 - emax * Ctumor_sn38 / (ec50 + Ctumor_sn38)
    d/dt(tumor_vol) <- kgrow * tumor_vol * log(tsmax / tumor_vol) * effect
    tumor_vol(0) <- 1

    tumor_vol ~ prop(propSd_tumor_vol)
    Ctumor_sn38 ~ prop(propSd_Ctumor_sn38)
  })
}
