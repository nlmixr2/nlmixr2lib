HernandezLozano_2019_erlotinib_4c_aif <- function() {
  description <- paste(
    "Hepatobiliary PET kinetics of [11C]erlotinib in healthy adults, four-compartment (4C)",
    "model: arterial blood input (not a state), a hepatocyte compartment, an intrahepatic bile",
    "duct compartment (k3) that drains (k5) into an extrahepatic bile duct / gall bladder",
    "(eBD/GB) sink. The liver input function is the sampled radial-artery whole-blood",
    "radioactivity (arterial input function, AIF). Deterministic mean-parameter model: each",
    "rate constant is the arithmetic mean of the per-subject weighted least-squares estimates",
    "for one of four scan conditions (baseline and post-erlotinib scans of one six-subject",
    "study; baseline and post-rifampicin scans of a second six-subject study), selected by",
    "covariates. Outputs are the liver and eBD/GB region-of-interest radioactivity",
    "concentrations; the blood input is an exogenous forcing function.",
    sep = " "
  )
  reference <- paste(
    "Hernandez Lozano I, Karch R, Bauer M, Blaickner M, Matsuda A, Wulkersdorfer B, Hacker M,",
    "Zeitlinger M, Langer O. Towards Improved Pharmacokinetic Models for the Analysis of",
    "Transporter-Mediated Hepatic Disposition of Drug Molecules with Positron Emission Tomography.",
    "AAPS J. 2019;21(4):61. doi:10.1208/s12248-019-0323-0",
    sep = " "
  )
  vignette <- "HernandezLozano_2019_erlotinib"

  units <- list(
    time = "min",
    dosing = "n/a (no dosing event; the arterial blood radioactivity CP_ERLOTINIB_KBQML is a forcing function)",
    concentration = "kBq/mL"
  )

  compartmentData <- list(
    liver = list(analyte = "[11C]erlotinib", units = "kBq", specimen = "tissue", verified = TRUE),
    bile_transit = list(analyte = "[11C]erlotinib", units = "kBq", specimen = "bile", verified = TRUE),
    gallbladder = list(analyte = "[11C]erlotinib", units = "kBq", specimen = "bile", verified = TRUE)
  )

  covariateData <- list(
    CP_ERLOTINIB_KBQML = list(
      description = "Arterial whole-blood [11C]erlotinib radioactivity concentration (input-function forcing function)",
      units = "kBq/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying error-free regressor, declared linear() so rxode2 interpolates linearly between",
        "rows. Hernandez Lozano 2019 Methods: blood samples were drawn from the radial artery in",
        "parallel to the 90-min PET scan and the radial-artery concentration is assumed equal to the",
        "hepatic-artery concentration C_HA(t), which is used directly as the input function C_in(t).",
        "The per-subject arterial curves are not tabulated; supply an observed or assumed curve.",
        "Must be supplied on every record, starting with 0 at the tracer injection (time 0).",
        sep = " "
      ),
      source_name = "C_HA(t) / C_in(t)"
    ),
    ORGVOL_LIVER = list(
      description = "Physiological total liver volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Hernandez Lozano 2019 Eq. 9-11: V_liver, the physiological liver volume, calculated per",
        "subject by the body-surface-area / body-weight method of Vauthey et al. 2002 (reference 28).",
        "Sets the blood volume in the liver (fvasc_liver * V_liver) and converts the eBD/GB amount into",
        "the eBD/GB concentration; the liver ROI concentration is independent of it.",
        sep = " "
      ),
      source_name = "V_liver"
    ),
    ORGVOL_EBDGB = list(
      description = "Volume of the extrahepatic bile duct / gall bladder region of interest",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Hernandez Lozano 2019 Eq. 12: V_eBD/GB, derived per subject from the PET images; values not reported.",
      source_name = "V_eBD/GB"
    ),
    CONMED_ERLOTINIB = list(
      description = "Unlabelled therapeutic erlotinib given before the [11C]erlotinib PET scan",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = second scan of the [11C]erlotinib + erlotinib study, acquired 3 h after a 300 mg oral dose",
        "of unlabelled erlotinib; 0 = baseline scan. Selects the second-scan mean parameter set within",
        "STUDY_BAUER2018RIF = 0 and is ignored when STUDY_BAUER2018RIF = 1.",
        sep = " "
      ),
      source_name = "Second scan (after oral erlotinib)"
    ),
    CONMED_RIFAMPICIN_SD = list(
      description = "Single-dose rifampicin infused before the [11C]erlotinib PET scan",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = second scan of the [11C]erlotinib + rifampicin study, acquired after a 30-min i.v. infusion",
        "of 600 mg rifampicin (acute OATP/SLCO inhibition); 0 = baseline scan. Selects the second-scan",
        "mean parameter set within STUDY_BAUER2018RIF = 1 and is ignored when STUDY_BAUER2018RIF = 0.",
        sep = " "
      ),
      source_name = "Second scan (after i.v. rifampicin)"
    ),
    STUDY_BAUER2018RIF = list(
      description = "Subject belongs to the [11C]erlotinib + rifampicin study cohort (Bauer et al. 2018, Mol Pharm)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = [11C]erlotinib + rifampicin study (subjects p37-p42, Hernandez Lozano 2019 reference 19);",
        "0 = [11C]erlotinib + erlotinib study (subjects p12-p30, reference 18). The two studies were",
        "fitted and summarised separately, so each has its own baseline-scan mean parameter set.",
        sep = " "
      ),
      source_name = "[11C]erlotinib + rifampicin data set"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12,
    n_studies = 2,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "i.v. bolus microdose of [11C]erlotinib for each of two 90-min dynamic abdominal PET scans per subject;",
      "second scan 3 h after 300 mg oral unlabelled erlotinib (study 1) or after 600 mg i.v. rifampicin over 30 min (study 2)",
      sep = " "
    ),
    regions = "Austria (Medical University of Vienna; EudraCT 2015-001593-18)",
    notes = paste(
      "Six healthy volunteers per study, each with a baseline and a second [11C]erlotinib PET scan about 7 days apart",
      "(Hernandez Lozano 2019 Methods, Data Sets; data previously published by Bauer et al. 2018, references 18-19).",
      "Demographics are not reported in this paper. Each subject and scan was fitted individually",
      "(MATLAB lsqnonlin, weighted least squares); no population model was fitted.",
      sep = " "
    )
  )

  ini({
    # Mean of the six per-subject estimates printed in Hernandez Lozano 2019 Table II / Table IV,
    # 4C model, 'Arterial input function' columns, 'Mean +/- SD' rows. The per-subject fits are a
    # two-stage analysis, so these are summary means of individual estimates, held fixed
    # here; the SD across subjects is quoted in each comment but carries no random effect.

    # [11C]erlotinib + erlotinib study, baseline scan
    kinf_liver_erlbl <- fixed(2.391)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); erlotinib-study baseline scan") # Table II 'Baseline scan' Mean +/- SD = 2.391 +/- 0.372 1/min
    keff_liver_erlbl <- fixed(0.060)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); erlotinib-study baseline scan") # Table II 'Baseline scan' Mean +/- SD = 0.060 +/- 0.021 1/min
    kbile_erlbl <- fixed(0.077)
    label("Hepatocyte-to-intrahepatic-bile-duct transfer rate constant k3 (1/min); erlotinib-study baseline scan") # Table II 'Baseline scan' Mean +/- SD = 0.077 +/- 0.106 1/min
    ktr_bile_erlbl <- fixed(0.013)
    label("Intrahepatic-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k5 (1/min); erlotinib-study baseline scan") # Table II 'Baseline scan' Mean +/- SD = 0.013 +/- 0.013 1/min

    # [11C]erlotinib + erlotinib study, second scan 3 h after 300 mg oral erlotinib
    kinf_liver_erlpost <- fixed(2.438)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); scan after oral erlotinib") # Table II 'Second scan (after oral erlotinib)' Mean +/- SD = 2.438 +/- 0.323 1/min
    keff_liver_erlpost <- fixed(0.191)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); scan after oral erlotinib") # Table II 'Second scan (after oral erlotinib)' Mean +/- SD = 0.191 +/- 0.031 1/min
    kbile_erlpost <- fixed(0.009)
    label("Hepatocyte-to-intrahepatic-bile-duct transfer rate constant k3 (1/min); scan after oral erlotinib") # Table II 'Second scan (after oral erlotinib)' Mean +/- SD = 0.009 +/- 0.006 1/min
    ktr_bile_erlpost <- fixed(0.013)
    label("Intrahepatic-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k5 (1/min); scan after oral erlotinib") # Table II 'Second scan (after oral erlotinib)' Mean +/- SD = 0.013 +/- 0.007 1/min

    # [11C]erlotinib + rifampicin study, baseline scan
    kinf_liver_rifbl <- fixed(2.039)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); rifampicin-study baseline scan") # Table IV 'Baseline scan' Mean +/- SD = 2.039 +/- 0.467 1/min
    keff_liver_rifbl <- fixed(0.058)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); rifampicin-study baseline scan") # Table IV 'Baseline scan' Mean +/- SD = 0.058 +/- 0.056 1/min
    kbile_rifbl <- fixed(0.063)
    label("Hepatocyte-to-intrahepatic-bile-duct transfer rate constant k3 (1/min); rifampicin-study baseline scan") # Table IV 'Baseline scan' Mean +/- SD = 0.063 +/- 0.023 1/min
    ktr_bile_rifbl <- fixed(0.013)
    label("Intrahepatic-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k5 (1/min); rifampicin-study baseline scan") # Table IV 'Baseline scan' Mean +/- SD = 0.013 +/- 0.008 1/min

    # [11C]erlotinib + rifampicin study, second scan after 600 mg i.v. rifampicin
    kinf_liver_rifpost <- fixed(1.574)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); scan after i.v. rifampicin") # Table IV 'Second scan (after i.v. rifampicin)' Mean +/- SD = 1.574 +/- 0.440 1/min
    keff_liver_rifpost <- fixed(0.032)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); scan after i.v. rifampicin") # Table IV 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.032 +/- 0.027 1/min
    kbile_rifpost <- fixed(0.034)
    label("Hepatocyte-to-intrahepatic-bile-duct transfer rate constant k3 (1/min); scan after i.v. rifampicin") # Table IV 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.034 +/- 0.033 1/min
    ktr_bile_rifpost <- fixed(0.078)
    label("Intrahepatic-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k5 (1/min); scan after i.v. rifampicin") # Table IV 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.078 +/- 0.097 1/min

    # Fixed physiological constants (not estimated)
    fvasc_liver <- fixed(0.25)
    label("Blood volume fraction of the liver, Vb (mL blood / mL liver)") # Methods Eq. 7 and Eq. 11: Vb set to 0.25 (reference 24)

    # No residual-error model is reported (weighted least squares with frame-length and
    # decay-based weights, Eq. 23); zero residual SDs keep the model a typical-value predictor.
    propSd_Cliver <- fixed(0)
    label("Proportional residual SD of the liver ROI concentration (fraction); not reported") # not reported (Eqs. 21-23 WLS)
    propSd_Cgallbladder <- fixed(0)
    label("Proportional residual SD of the eBD/GB ROI concentration (fraction); not reported") # not reported (Eqs. 21-23 WLS)
  })

  model({
    # Arterial whole-blood radioactivity, an error-free forcing function interpolated linearly.
    linear(CP_ERLOTINIB_KBQML)

    # Scan-condition selection: exactly one weight is 1 for any valid covariate combination.
    sel_erlbl <- (1 - STUDY_BAUER2018RIF) * (1 - CONMED_ERLOTINIB)
    sel_erlpost <- (1 - STUDY_BAUER2018RIF) * CONMED_ERLOTINIB
    sel_rifbl <- STUDY_BAUER2018RIF * (1 - CONMED_RIFAMPICIN_SD)
    sel_rifpost <- STUDY_BAUER2018RIF * CONMED_RIFAMPICIN_SD
    kinf_liver <-
      sel_erlbl * kinf_liver_erlbl + sel_erlpost * kinf_liver_erlpost +
      sel_rifbl * kinf_liver_rifbl + sel_rifpost * kinf_liver_rifpost
    keff_liver <-
      sel_erlbl * keff_liver_erlbl + sel_erlpost * keff_liver_erlpost +
      sel_rifbl * keff_liver_rifbl + sel_rifpost * keff_liver_rifpost
    kbile <-
      sel_erlbl * kbile_erlbl + sel_erlpost * kbile_erlpost +
      sel_rifbl * kbile_rifbl + sel_rifpost * kbile_rifpost
    ktr_bile <-
      sel_erlbl * ktr_bile_erlbl + sel_erlpost * ktr_bile_erlpost +
      sel_rifbl * ktr_bile_rifbl + sel_rifpost * ktr_bile_rifpost

    # AIF approach: the input function is the sampled arterial curve (Eq. 7, C_in = C_HA).
    cin <- CP_ERLOTINIB_KBQML

    # Amount of tracer in the hepatic blood space, Eq. 7: X_blood = V_blood * C_in with
    # V_blood = Vb * V_liver.
    xblood <- fvasc_liver * ORGVOL_LIVER * cin

    # Eqs. 4-6 (amounts in kBq).
    d/dt(liver) <- kinf_liver * xblood - (keff_liver + kbile) * liver
    d/dt(bile_transit) <- kbile * liver - ktr_bile * bile_transit
    d/dt(gallbladder) <- ktr_bile * bile_transit

    # Liver ROI = hepatic blood + hepatocytes + intrahepatic bile duct, Eq. 11; eBD/GB ROI, Eq. 12.
    Cliver <- fvasc_liver * cin + (liver + bile_transit) / ORGVOL_LIVER
    Cgallbladder <- gallbladder / ORGVOL_EBDGB

    Cliver ~ prop(propSd_Cliver)
    Cgallbladder ~ prop(propSd_Cgallbladder)
  })
}
