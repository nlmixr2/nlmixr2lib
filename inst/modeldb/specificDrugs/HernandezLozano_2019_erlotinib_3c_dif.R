HernandezLozano_2019_erlotinib_3c_dif <- function() {
  description <- paste(
    "Hepatobiliary PET kinetics of [11C]erlotinib in healthy adults, three-compartment (3C)",
    "model: arterial blood input (not a state), one liver-tissue compartment, and an",
    "extrahepatic bile duct / gall bladder (eBD/GB) sink fed directly from the liver tissue",
    "(k3). The liver input function is a dual input function (DIF): 25% hepatic-artery blood",
    "(the sampled arterial curve) plus 75% portal-vein blood, the portal-vein curve being the",
    "convolution of the arterial curve with the impulse-response function beta/(beta+t)^2 (beta",
    "estimated per subject). Deterministic mean-parameter model: each rate constant is the",
    "arithmetic mean of the per-subject weighted least-squares estimates for one of four scan",
    "conditions (baseline and post-erlotinib scans of one six-subject study; baseline and post-",
    "rifampicin scans of a second six-subject study), selected by covariates. Outputs are the",
    "liver and eBD/GB region-of-interest radioactivity concentrations; the blood input is an",
    "exogenous forcing function.",
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

  # pvconv1..pvconv14 hold the sum-of-exponentials quadrature of the portal-vein
  # convolution integral (Eqs. 1-2); they have no anatomical meaning of their own.
  paper_specific_compartments <- paste0("pvconv", 1:14)

  units <- list(
    time = "min",
    dosing = "n/a (no dosing event; the arterial blood radioactivity CP_ERLOTINIB_KBQML is a forcing function)",
    concentration = "kBq/mL"
  )

  compartmentData <- list(
    liver = list(analyte = "[11C]erlotinib", units = "kBq", specimen = "tissue", verified = TRUE),
    gallbladder = list(analyte = "[11C]erlotinib", units = "kBq", specimen = "bile", verified = TRUE),
    pvconv1 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv2 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv3 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv4 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv5 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv6 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv7 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv8 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv9 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv10 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv11 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv12 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv13 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE),
    pvconv14 = list(analyte = "[11C]erlotinib", units = "kBq/mL", specimen = "not applicable", verified = TRUE)
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
        "hepatic-artery concentration C_HA(t); the model derives the portal-vein and dual input from it.",
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
    # Mean of the six per-subject estimates printed in Hernandez Lozano 2019 Table I / Table III,
    # 3C model, 'Dual input function' columns, 'Mean +/- SD' rows. The per-subject fits are a
    # two-stage analysis, so these are summary means of individual estimates, held fixed
    # here; the SD across subjects is quoted in each comment but carries no random effect.

    # [11C]erlotinib + erlotinib study, baseline scan
    kinf_liver_erlbl <- fixed(4.208)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); erlotinib-study baseline scan") # Table I 'Baseline scan' Mean +/- SD = 4.208 +/- 0.567 1/min
    keff_liver_erlbl <- fixed(0.012)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); erlotinib-study baseline scan") # Table I 'Baseline scan' Mean +/- SD = 0.012 +/- 0.012 1/min
    kbile_erlbl <- fixed(0.007)
    label("Liver-tissue-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k3 (1/min); erlotinib-study baseline scan") # Table I 'Baseline scan' Mean +/- SD = 0.007 +/- 0.003 1/min
    beta_pv_erlbl <- fixed(1.283)
    label("Portal-vein impulse-response parameter beta (min); erlotinib-study baseline scan") # Table I 'Baseline scan' Mean +/- SD = 1.283 +/- 0.490 min

    # [11C]erlotinib + erlotinib study, second scan 3 h after 300 mg oral erlotinib
    kinf_liver_erlpost <- fixed(3.805)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); scan after oral erlotinib") # Table I 'Second scan (after oral erlotinib)' Mean +/- SD = 3.805 +/- 0.668 1/min
    keff_liver_erlpost <- fixed(0.107)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); scan after oral erlotinib") # Table I 'Second scan (after oral erlotinib)' Mean +/- SD = 0.107 +/- 0.038 1/min
    kbile_erlpost <- fixed(0.003)
    label("Liver-tissue-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k3 (1/min); scan after oral erlotinib") # Table I 'Second scan (after oral erlotinib)' Mean +/- SD = 0.003 +/- 0.002 1/min
    beta_pv_erlpost <- fixed(1.311)
    label("Portal-vein impulse-response parameter beta (min); scan after oral erlotinib") # Table I 'Second scan (after oral erlotinib)' Mean +/- SD = 1.311 +/- 0.561 min

    # [11C]erlotinib + rifampicin study, baseline scan
    kinf_liver_rifbl <- fixed(3.943)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); rifampicin-study baseline scan") # Table III 'Baseline scan' Mean +/- SD = 3.943 +/- 0.697 1/min
    keff_liver_rifbl <- fixed(0.007)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); rifampicin-study baseline scan") # Table III 'Baseline scan' Mean +/- SD = 0.007 +/- 0.009 1/min
    kbile_rifbl <- fixed(0.010)
    label("Liver-tissue-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k3 (1/min); rifampicin-study baseline scan") # Table III 'Baseline scan' Mean +/- SD = 0.010 +/- 0.004 1/min
    beta_pv_rifbl <- fixed(1.080)
    label("Portal-vein impulse-response parameter beta (min); rifampicin-study baseline scan") # Table III 'Baseline scan' Mean +/- SD = 1.080 +/- 0.590 min

    # [11C]erlotinib + rifampicin study, second scan after 600 mg i.v. rifampicin
    kinf_liver_rifpost <- fixed(3.212)
    label("Blood-to-liver-tissue transfer rate constant k1 (1/min); scan after i.v. rifampicin") # Table III 'Second scan (after i.v. rifampicin)' Mean +/- SD = 3.212 +/- 1.091 1/min
    keff_liver_rifpost <- fixed(0.009)
    label("Liver-tissue-to-blood transfer rate constant k2 (1/min); scan after i.v. rifampicin") # Table III 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.009 +/- 0.005 1/min
    kbile_rifpost <- fixed(0.009)
    label("Liver-tissue-to-extrahepatic-bile-duct/gall-bladder transfer rate constant k3 (1/min); scan after i.v. rifampicin") # Table III 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.009 +/- 0.004 1/min
    beta_pv_rifpost <- fixed(0.971)
    label("Portal-vein impulse-response parameter beta (min); scan after i.v. rifampicin") # Table III 'Second scan (after i.v. rifampicin)' Mean +/- SD = 0.971 +/- 0.569 min

    # Fixed physiological constants (not estimated)
    fvasc_liver <- fixed(0.25)
    label("Blood volume fraction of the liver, Vb (mL blood / mL liver)") # Methods Eq. 7 and Eq. 11: Vb set to 0.25 (reference 24)
    fha <- fixed(0.25)
    label("Hepatic-artery fraction of liver blood supply, fHA (fraction)") # Methods Eq. 3: fHA fixed at 0.25
    fpv <- fixed(0.75)
    label("Portal-vein fraction of liver blood supply, fPV (fraction)") # Methods Eq. 3: fPV fixed at 0.75

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
    beta_pv <-
      sel_erlbl * beta_pv_erlbl + sel_erlpost * beta_pv_erlpost +
      sel_rifbl * beta_pv_rifbl + sel_rifpost * beta_pv_rifpost

    # Portal-vein input C_PV (Eqs. 1-2), the sum of the pvconv quadrature states solved below.
    cpv <- pvconv1 + pvconv2 + pvconv3 + pvconv4 + pvconv5 + pvconv6 + pvconv7 +
      pvconv8 + pvconv9 + pvconv10 + pvconv11 + pvconv12 + pvconv13 + pvconv14

    # Dual input function, Eq. 3.
    cin <- fha * CP_ERLOTINIB_KBQML + fpv * cpv

    # Amount of tracer in the hepatic blood space, Eq. 7: X_blood = V_blood * C_in with
    # V_blood = Vb * V_liver.
    xblood <- fvasc_liver * ORGVOL_LIVER * cin

    # Eqs. 15-16 (amounts in kBq).
    d/dt(liver) <- kinf_liver * xblood - (keff_liver + kbile) * liver
    d/dt(gallbladder) <- kbile * liver

    # Eqs. 1-2: C_PV(t) = integral_0^t h(t - tau) C_HA(tau) dtau with h(t) = beta / (beta + t)^2.
    # Writing h(t) = (1/beta) * integral_0^inf u exp(-u) exp(-u t / beta) du and applying the
    # trapezoid rule in x = log(u) on x = -9, -8, ..., 4 gives
    # h(t) ~= sum_k c_k / beta * exp(-u_k t / beta), u_k = exp(x_k), c_k = exp(2 x_k - exp(x_k)),
    # within 0.5% of the exact kernel for t = 0-120 min and any beta. Each exponential term of the
    # convolution is one linear ODE.
    d/dt(pvconv1) <- (exp(-18 - exp(-9)) * CP_ERLOTINIB_KBQML - exp(-9) * pvconv1) / beta_pv
    d/dt(pvconv2) <- (exp(-16 - exp(-8)) * CP_ERLOTINIB_KBQML - exp(-8) * pvconv2) / beta_pv
    d/dt(pvconv3) <- (exp(-14 - exp(-7)) * CP_ERLOTINIB_KBQML - exp(-7) * pvconv3) / beta_pv
    d/dt(pvconv4) <- (exp(-12 - exp(-6)) * CP_ERLOTINIB_KBQML - exp(-6) * pvconv4) / beta_pv
    d/dt(pvconv5) <- (exp(-10 - exp(-5)) * CP_ERLOTINIB_KBQML - exp(-5) * pvconv5) / beta_pv
    d/dt(pvconv6) <- (exp(-8 - exp(-4)) * CP_ERLOTINIB_KBQML - exp(-4) * pvconv6) / beta_pv
    d/dt(pvconv7) <- (exp(-6 - exp(-3)) * CP_ERLOTINIB_KBQML - exp(-3) * pvconv7) / beta_pv
    d/dt(pvconv8) <- (exp(-4 - exp(-2)) * CP_ERLOTINIB_KBQML - exp(-2) * pvconv8) / beta_pv
    d/dt(pvconv9) <- (exp(-2 - exp(-1)) * CP_ERLOTINIB_KBQML - exp(-1) * pvconv9) / beta_pv
    d/dt(pvconv10) <- (exp(0 - exp(0)) * CP_ERLOTINIB_KBQML - exp(0) * pvconv10) / beta_pv
    d/dt(pvconv11) <- (exp(2 - exp(1)) * CP_ERLOTINIB_KBQML - exp(1) * pvconv11) / beta_pv
    d/dt(pvconv12) <- (exp(4 - exp(2)) * CP_ERLOTINIB_KBQML - exp(2) * pvconv12) / beta_pv
    d/dt(pvconv13) <- (exp(6 - exp(3)) * CP_ERLOTINIB_KBQML - exp(3) * pvconv13) / beta_pv
    d/dt(pvconv14) <- (exp(8 - exp(4)) * CP_ERLOTINIB_KBQML - exp(4) * pvconv14) / beta_pv

    # Liver ROI = hepatic blood + liver tissue, Eq. 19; eBD/GB ROI, Eq. 12.
    Cliver <- fvasc_liver * cin + liver / ORGVOL_LIVER
    Cgallbladder <- gallbladder / ORGVOL_EBDGB

    Cliver ~ prop(propSd_Cliver)
    Cgallbladder ~ prop(propSd_Cgallbladder)
  })
}
