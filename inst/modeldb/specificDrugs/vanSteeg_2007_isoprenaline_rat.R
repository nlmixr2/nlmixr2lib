vanSteeg_2007_isoprenaline_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment micro-constant IV PK of the beta-adrenoceptor agonist isoprenaline in",
    "conscious male Wistar-Kyoto (WKY) rats, linked to a DIRECT (no effect compartment)",
    "sigmoid Emax model for the resulting increase in heart rate.",
    "Isoprenaline is so potent that a maximal chronotropic response occurs at plasma",
    "concentrations below the assay limit of quantification, so the PK and the PD were run as",
    "SEPARATE experiments and the PD-experiment plasma concentrations were predicted from this",
    "PK model rather than measured.",
    "This model supplies the tachycardia background against which the companion",
    "S(-)-atenolol models were built:",
    "modellib('vanSteeg_2007_atenolol_iso_rat') and modellib('vanSteeg_2007_atenolol_noniso_rat').",
    "Note the two non-standard variability scales, both taken verbatim from the paper:",
    "PK inter-individual variability is EXPONENTIAL (Equation 1) whereas PD inter-individual",
    "variability is ADDITIVE in beats/min (Equation 5), and the PD residual error is additive in",
    "beats/min despite the paper's Equation 6 printing a proportional form (see the ini() notes).",
    sep = " "
  )

  reference <- paste(
    "van Steeg TJ, Freijer J, Danhof M, de Lange ECM.",
    "Pharmacokinetic-pharmacodynamic modelling of S(-)-atenolol in rats:",
    "reduction of isoprenaline-induced tachycardia as a continuous pharmacodynamic endpoint.",
    "British Journal of Pharmacology. 2007;151(3):356-366.",
    "doi:10.1038/sj.bjp.0707234. PMCID: PMC2013984.",
    "Isoprenaline PK parameters are Table 1; isoprenaline PD parameters are Table 2.",
    "Companion model files from the same paper:",
    "modellib('vanSteeg_2007_atenolol_iso_rat'), modellib('vanSteeg_2007_atenolol_noniso_rat').",
    sep = " "
  )

  vignette <- "vanSteeg_2007_atenolol_isoprenaline_rat"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  # No covariates. The paper fits absolute (not weight-normalised) V1 and absolute
  # micro-constants, so body weight enters only when converting the reported ug/kg
  # dose into the absolute infused amount. See population$notes.
  covariateData <- list()

  compartmentData <- list(
    central     = list(analyte = "isoprenaline", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoprenaline", units = "ng", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "rat (Wistar-Kyoto, WKY; male)",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = "not reported; acclimatised at least 5 days, instrumented 7 days before the experiment",
    weight_range   = "294 g +/- 49 (mean +/- SD across the 39 WKY rats of this study)",
    sex_female_pct = 0,
    disease_state  = "normotensive (WKY) conscious rats; no disease model",
    dose_range     = paste(
      "PK experiment: isoprenaline 25 ug/kg (n = 7) or 50 ug/kg (n = 7) as a 10-minute IV infusion.",
      "PD experiment: continuous IV infusions of 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1 and 2.5 ug/kg/h",
      "(n = 6 rats, sequential rates), plus steady-state infusions of 5 ug/kg/h (n = 17) and",
      "10 ug/kg/h (n = 8) contributed from other experiments.",
      sep = " "
    ),
    regions        = "The Netherlands (Leiden University)",
    biomarkers     = paste(
      "Heart rate captured continuously from a femoral-artery pressure signal (P10EZ-1 transducer,",
      "Spike 2 acquisition) and used directly as the PD endpoint.",
      "Plasma isoprenaline by reversed-phase HPLC with electrochemical detection after ion-paired",
      "liquid-liquid extraction; calibration range 1-1000 ng/mL, LOQ 0.3 ng/mL (150 uL plasma);",
      "intra- and inter-assay variability 5% and 12%.",
      sep = " "
    ),
    notes = paste(
      "SUBJECT COUNT. n_subjects = 45 is the union of the two isoprenaline experiments:",
      "14 rats in the PK experiment (two dose groups of 7) and 31 rats in the PD experiment",
      "(6 in the sequential-infusion dose-ranging arm plus 17 and 8 at 5 and 10 ug/kg/h from",
      "other experiments). The 39 WKY rats stated in Methods/Animals cover only this study's own",
      "experiments (isoprenaline PK 14 + isoprenaline PD 6 + atenolol 9 + 8 = 37); the 17 and 8",
      "animals at 5 and 10 ug/kg/h are explicitly 'additional data from other experiments'.",
      "DOSE UNITS: the paper reports ug/kg doses but fits an absolute V1 (mL) and absolute",
      "micro-constants, so a body weight is needed to build an event table. The validation vignette",
      "doses the paper's own typical 294 g rat (Discussion), for which 25 and 50 ug/kg are 7350 and",
      "14700 ng. The body weight is NOT a model covariate.",
      "PK-PD SEPARATION: 'Isoprenaline is an extremely potent beta-adrenoceptor agonist and a maximal",
      "increase in heart rate in rats is obtained with plasma concentrations below the limit of",
      "quantification. Therefore, the PK and the PD of isoprenaline were determined in separate",
      "experiments and the concentrations in the PD experiments were predicted using the PK model.'",
      "The PD is therefore fitted against MODEL-PREDICTED, not measured, plasma concentrations.",
      "NO EFFECT COMPARTMENT: unlike the companion S(-)-atenolol ISO model, the isoprenaline PD was",
      "fitted as a direct-effect sigmoid Emax on the plasma concentration; the paper introduces the",
      "effect compartment only to resolve the hysteresis seen with S(-)-atenolol.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETICS -- van Steeg 2007 Table 1 ("Population estimates
    # of PK parameters (micro constants) for isoprenaline in the rat and
    # estimates of the BS replicates"), 'Value (CV)' column. A standard
    # two-compartment model parameterised directly on micro-constants
    # (Methods/Pharmacokinetics: "PK compartmental analysis for
    # isoprenaline was performed by fitting a standard two compartment
    # model to the concentration-time profiles").
    #
    # The '(CV)' column of every table in this paper is the RELATIVE
    # STANDARD ERROR of the estimate, not a between-animal CV. The paper
    # confirms this twice: Results/Isoprenaline says "The coefficient of
    # variation of the parameter estimates varied between 17 and 47%",
    # which brackets exactly the Table 1 range (17.6% to 47.1%), and the
    # S(-)-atenolol section says "coefficients of variation ranging
    # between 3 and 36%", bracketing exactly the Table 3 range
    # (3.4% to 36.7%).
    # ==================================================================
    lkel <- log(1.32)   ; label("Elimination rate constant k10 (1/min)")               # Table 1 K10 = 1.32 1/min (RSE 20.5%; bootstrap 1.33)
    lk12 <- log(0.391)  ; label("Central-to-peripheral rate constant k12 (1/min)")     # Table 1 K12 = 0.391 1/min (RSE 36.1%; bootstrap 0.412)
    lk21 <- log(0.177)  ; label("Peripheral-to-central rate constant k21 (1/min)")     # Table 1 K21 = 0.177 1/min (RSE 17.6%; bootstrap 0.177)
    lvc  <- log(79.6)   ; label("Central volume of distribution V1 (mL)")              # Table 1 V1 = 79.6 mL (RSE 24.4%; bootstrap 83.3)

    # ------------------------------------------------------------------
    # PK inter-individual variability. Exponential, per Equation 1:
    #   P_i = theta * exp(eta_i),  eta_i ~ N(0, omega^2)
    # Table 1 reports the VARIANCES omega^2 directly. "Estimation of
    # interindividual variability was possible for k12, k21 and V1" --
    # there is no IIV on k10, so lkel carries no eta.
    # ------------------------------------------------------------------
    etalk12 ~ 0.213     # Table 1 omega^2_K12 = 0.213 (RSE 46.8%; bootstrap 0.215) -> CV 48.9%
    etalk21 ~ 0.114     # Table 1 omega^2_K21 = 0.114 (RSE 43.8%; bootstrap 0.121). The bootstrap RSE is 103.0%, flagged in Results as an uncertainty in estimating the inter-animal variability on k21; the bootstrap point estimate is nevertheless nearly identical to the final one.
    etalvc  ~ 0.0982    # Table 1 omega^2_V1 = 0.0982 (RSE 47.1%; bootstrap 0.0900) -> CV 32.1%

    # ------------------------------------------------------------------
    # PK residual error. Proportional, per Equation 2:
    #   C_obs = C_pred * (1 + eps),  eps ~ N(0, sigma^2)
    # Table 1 reports the VARIANCE sigma^2 = 0.0796, so the SD that
    # nlmixr2's prop() wants is its square root (0.282, i.e. 28.2%).
    # The table row is labelled "sigma^2_PD" even though Table 1 is the
    # PK table; that subscript is a copy-paste slip carried across all
    # four tables of the paper (see the vignette Errata). Its placement
    # under the "Residual error" heading is unambiguous.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.0796) ; label("Proportional residual error on plasma isoprenaline (fraction)")   # Table 1 sigma^2 = 0.0796 (RSE 22.4%; bootstrap 0.0760)

    # ==================================================================
    # PHARMACODYNAMICS -- van Steeg 2007 Table 2 ("Population estimates
    # of PD parameters and variabilities for isoprenaline in the rat").
    # Sigmoid Emax on the plasma concentration, Equation 3:
    #   E = E0 + Emax * C^n / (EC50^n + C^n)
    # Emax is POSITIVE here: isoprenaline is an agonist and raises heart
    # rate. (The companion S(-)-atenolol models carry a negative Emax.)
    #
    # E0 and Emax are NOT log-transformed. The paper models the PD
    # parameters ADDITIVELY (Equation 5, P_i = theta + eta_i), unlike the
    # exponential PK Equation 1, so the eta must sit on the natural scale
    # of beats/min. The Table 2 magnitude confirms it: omega^2_E0 = 860
    # is an SD of 29.3 beats/min on a 374 beats/min baseline (7.8%), and
    # is uninterpretable as an exponential variance. The paper's own
    # Discussion corroborates the number -- "The variation in heart rate
    # owing to circadian rhythms, movement and stress is approximately
    # 30 b.p.m. in rats".
    #
    # EC50 and the Hill coefficient carry no IIV, are strictly positive,
    # and therefore take the ordinary log-transformed canonical form.
    # ==================================================================
    e0     <- 374       ; label("Baseline (drug-free) heart rate E0 (beats/min)")                                  # Table 2 E0 = 374 beats/min (RSE 1.9%); Results text gives 374.0 +/- 7.0
    emax   <- 130       ; label("Maximum isoprenaline effect on heart rate, Emax (beats/min; positive = increase)") # Table 2 Emax = 130 beats/min (RSE 5.9%); Results text gives 130 +/- 7.7
    lec50  <- log(0.0138) ; label("Plasma concentration at half-maximal effect, EC50 (ng/mL)")                      # Table 2 EC50 = 0.0138 ng/mL (RSE 31.9%); Results text rounds to 0.014 +/- 0.0044
    lhill  <- log(1.18) ; label("Hill coefficient n of the isoprenaline sigmoid (unitless)")                        # Table 2 n = 1.18 (RSE 19.3%); Results text gives 1.18 +/- 0.23

    # ------------------------------------------------------------------
    # PD inter-individual variability, ADDITIVE in beats/min per
    # Equation 5. "interanimal variability was observed for the baseline
    # only", so Emax, EC50 and n carry no eta.
    # ------------------------------------------------------------------
    etae0 ~ 860         # Table 2 omega^2_E0 = 860 (beats/min)^2 (RSE 29.8%) -> SD 29.3 beats/min

    # ------------------------------------------------------------------
    # PD residual error, ADDITIVE in beats/min. Table 2 reports the
    # VARIANCE sigma^2 = 409, so the SD is 20.2 beats/min.
    #
    # The paper's Equation 6 prints a PROPORTIONAL form for the PD
    # residual, but that equation is a verbatim copy-paste of the PK
    # Equation 2 -- its own explanatory sentence still reads "in which
    # C_obs,ij is the jth observed CONCENTRATION", inside the
    # Pharmacodynamics section. The reported magnitudes settle it:
    # sigma^2 = 409 read proportionally would be a residual CV of 2022%,
    # whereas read additively it is 20.2 beats/min on a 374-504
    # beats/min signal (4-5%), matching a continuously recorded
    # heart-rate trace. All three of this paper's PD residual variances
    # (409 here, 747 and 250 in Table 4) behave the same way. Encoded
    # additively; recorded in the vignette Errata.
    # ------------------------------------------------------------------
    addSd_hr <- sqrt(409) ; label("Additive residual error on heart rate (beats/min)")   # Table 2 sigma^2 = 409 (beats/min)^2 (RSE 18.3%)
  })

  model({
    # ---- 1. Individual parameters ----
    # PK: exponential etas (Equation 1). PD: additive etas (Equation 5).
    kel  <- exp(lkel)
    k12  <- exp(lk12 + etalk12)
    k21  <- exp(lk21 + etalk21)
    vc   <- exp(lvc  + etalvc)

    e0_i <- e0 + etae0
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---- 2. Two-compartment PK on micro-constants ----
    # Amount form, so rxode2's own zero-order infusion machinery supplies
    # the infusion rate; dose `central` with dur = 10 min for the PK
    # experiment, or with an open-ended `rate` for the PD infusions.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1

    Cc <- central / vc

    # ---- 3. Direct-effect sigmoid Emax on heart rate (Equation 3) ----
    # No effect compartment: the paper introduces one only for
    # S(-)-atenolol, whose concentration-effect plot showed hysteresis.
    hr <- e0_i + emax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- 4. Observations ----
    Cc ~ prop(propSd)
    hr ~ add(addSd_hr)
  })
}
