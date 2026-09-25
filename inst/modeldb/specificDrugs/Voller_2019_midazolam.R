Voller_2019_midazolam <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for midazolam in preterm",
    "neonates (55 neonates, gestational age 24.0-33.6 weeks, actual body weight",
    "0.6-4.3 kg) pooled from the prospective multicentre DINO study and the",
    "de Wildt 1998 data set. Actual body weight enters clearance as a power",
    "function with an estimated exponent of 1.69 and central volume as a linear",
    "function (exponent fixed at 1), both centred on the 1.12 kg cohort median.",
    "The steeply supra-allometric clearance exponent reflects CYP3A maturation",
    "over the preterm weight range and is the basis of the paper's finding that",
    "the registered flat 0.03 mg/(kg*h) infusion does not give equal exposure.",
    "No covariate other than actual body weight was retained."
  )
  reference <- paste(
    "Voller S, Flint RB, Beggah F, Reiss I, Andriessen P, Zimmermann LJI,",
    "van den Anker JN, Liem KD, Koch BCP, de Wildt S, Knibbe CAJ, Simons SHP",
    "(2019). Recently Registered Midazolam Doses for Preterm Neonates Do Not",
    "Lead to Equal Exposure: A Population Pharmacokinetic Model.",
    "J Clin Pharmacol 59(10):1300-1308. doi:10.1002/jcph.1429."
  )
  vignette <- "Voller_2019_midazolam"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Voller 2019 Methods ('Population Pharmacokinetic",
        "Analysis and Model Evaluation') records actual body weight at each day",
        "in the neonatal intensive care unit, with linear interpolation between",
        "available measurements and last-value-carried-forward for other",
        "covariates. Table 1 gives a combined-cohort median of 1.1 kg (range",
        "0.6-4.3 kg); the Table 2 footnote states the covariate-centring value",
        "explicitly as 'median WT = 1.12 kg', and 1.12 kg is the value used",
        "here. Enters with two distinct functional forms, both centred at",
        "1.12 kg: a power function on clearance,",
        "CL_i = 0.0737 * (WT_i / 1.12)^1.69 (Table 2), and a linear function on",
        "central volume, V_i = 1.03 * (WT_i / 1.12) (Table 2). Results state",
        "that a power exponent on volume was estimated close to 1, that the",
        "linear form fitted equally well (+0.2 points in objective function",
        "value) and was therefore carried forward -- so the volume exponent is",
        "encoded as fixed(1), not estimated. Actual body weight superseded",
        "birth weight, gestational age, postnatal age, postmenstrual age and",
        "sex in the stepwise covariate search; see covariatesDataExcluded."
      ),
      source_name = "WT"
    )
  )

  # Screened in the stepwise covariate modelling procedure (Voller 2019
  # Methods, forward inclusion p <= 0.01 / backward elimination p <= 0.001)
  # but NOT retained in the final model: Results states 'No other covariates
  # were found' after actual body weight entered on clearance and volume.
  # Documented here so the paper's covariate screen is preserved without
  # declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    WT_BIRTH = list(
      description = "Body weight at birth",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: combined-cohort median 0.95 kg (range 0.47-1.8 kg), reported",
        "for 29 of the 55 neonates. Screened but not retained; the Discussion",
        "attributes this to actual body weight being the better size descriptor",
        "once the observation period extends to 30 days."
      ),
      source_name = "BW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: combined-cohort median 27.3 weeks (range 24.0-33.6 weeks).",
        "Screened but not retained. Discussion attributes the null result to",
        "the correlation between age and weight and to the narrow gestational",
        "age range of the pooled cohort."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age at day of inclusion",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: combined-cohort median 9 days (range 1-88 days). Reported in",
        "days by the source; the PAGE/PNA register default is months. Screened",
        "but not retained."
      ),
      source_name = "PNA"
    ),
    PAGE = list(
      description = "Postmenstrual age, calculated from gestational age and postnatal age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Voller 2019 Methods lists postmenstrual age among the screened",
        "covariates, 'calculated using gestational age and postnatal age'. No",
        "summary statistics are tabulated for it. Screened but not retained.",
        "Reported on the weeks scale used throughout the neonatal literature,",
        "not the register-default months."
      ),
      source_name = "PMA"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table 1: 29 of 55 neonates female (53%). Voller 2019 Methods lists",
        "'gender' among the screened covariates. Screened but not retained."
      ),
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    central = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 2L,
    age_range = "Postnatal age 1-88 days at inclusion; gestational age at birth 24.0-33.6 weeks",
    age_median = "Postnatal age 9 days at inclusion; gestational age at birth 27.3 weeks",
    weight_range = "Actual body weight 0.6-4.3 kg; birth weight 0.47-1.8 kg",
    weight_median = "Actual body weight 1.1 kg (Table 1); 1.12 kg is the covariate-centring median given in the Table 2 footnote",
    sex_female_pct = 52.7,
    disease_state = paste(
      "Preterm neonates admitted to level III neonatal intensive care units and",
      "receiving intravenous midazolam for a stressful procedure or for",
      "continuous sedation. All were born before 32 weeks of gestation in study",
      "1; study 2 extended to 33.6 weeks. No target midazolam concentration has",
      "been defined for this population."
    ),
    dose_range = paste(
      "Study 1 (DINO): 123 intravenous doses in 29 neonates -- 59 given for",
      "stressful procedures as infusions of 30 minutes or less and 65 given as",
      "continuous sedation infusions longer than 30 minutes; dose chosen by the",
      "treating physician. Study 2 (de Wildt et al): a single intravenous dose",
      "of 0.1 mg/kg over 30 minutes in each of 26 neonates. 149 doses and 232",
      "plasma samples in total."
    ),
    regions = "Netherlands",
    ga_range = "24.0-33.6 weeks (study 1 median 26.7, study 2 median 28.1)",
    notes = paste(
      "Pooled from two studies (Voller 2019 Table 1). Study 1 is the",
      "prospective multicentre Drug dosage Improvement NeOnates (DINO) study",
      "(MEC-2014-067, NL47409.078.14, NCT02421068), September 2014 to September",
      "2017, in four Dutch level III neonatal intensive care units: 29 neonates",
      "(13 female, 16 male), 60 sparse opportunistic samples. Study 2 is the",
      "previously published de Wildt et al data set (MEC171.586/1998/125) from",
      "the Sophia Children's Hospital Rotterdam: 26 neonates (16 female, 10",
      "male), 172 densely sampled concentrations at 0.5, 1, 2, 4, 6, 12 and 24",
      "hours after a single 0.1 mg/kg intravenous dose. Two different assays",
      "were used -- gas chromatography with mass spectrometric detection for",
      "study 2 (lower limit of quantification 1 ug/L) and liquid",
      "chromatography-tandem mass spectrometry for study 1 (lower limit of",
      "quantification 4 ug/L) -- which the Discussion names as one likely",
      "contributor to the large unexplained interindividual variability.",
      "C-reactive protein and other inflammation markers were too sparsely",
      "sampled to be tested."
    )
  )

  ini({
    # Structural parameters from Voller 2019 Table 2, 'Final Model Estimate
    # (RSE%)' column. Both are the typical value for a neonate at the
    # covariate-centring median actual body weight of 1.12 kg (Table 2
    # footnote: 'p, population mean value of a parameter for an individual with
    # body weight of 1.12 kg').
    #
    # Clearance is printed in Table 2 as 0.0737 L/h; the Abstract restates the
    # same estimate as '73.7 mL/h for a neonate weighing 1.1 kg'. Bootstrap
    # median 0.0787 L/h (95%CI 0.0550-0.102), 1000 samples, 83.1% convergence.
    lcl <- log(0.0737)
    label("Clearance at the median actual body weight of 1.12 kg (L/h)")

    # Central volume 1.03 L at 1.12 kg; Abstract restates as '1.03 L for a
    # neonate weighing 1.1 kg'. Bootstrap median 1.11 L (95%CI 0.877-1.37).
    lvc <- log(1.03)
    label("Central volume of distribution at the median actual body weight of 1.12 kg (L)")

    # Covariate effects of actual body weight, both centred at 1.12 kg.
    #
    # Power exponent on clearance, estimated. Table 2 row
    # 'CL_i = CL_p x (WT_i/Median WT)^theta_WT'; theta_WT = 1.69 (RSE 10%),
    # bootstrap 1.67 (95%CI 0.406-2.47). Results: clearance was best predicted
    # by actual body weight using a power function (p < .001, -116 points in
    # objective function value).
    e_wt_cl <- 1.69
    label("Power exponent of (WT / 1.12 kg) on clearance (unitless)")

    # Exponent on central volume, FIXED at 1 (linear). Table 2 prints the
    # volume row as 'V_i = V_p x (WT_i/Median WT)' with no exponent symbol, and
    # Results states: 'Because the estimated exponent for volume was close to
    # 1, a linear influence of body weight on volume of distribution was tested
    # for comparison. This led to an equally good fit (+0.2 points in OFV) and
    # was therefore carried forward.' Encoded as fixed(1) to preserve the fact
    # that the linear form is a modelling choice, not an estimate.
    e_wt_vc <- fixed(1)
    label("Linear-form exponent of (WT / 1.12 kg) on central volume (unitless)")

    # Interindividual variability. Table 2 reports both etas as percentages:
    # 'On CL [%] 91.9 (15%) [8%]' and 'On V [%] 67.2 (16%) [17%]' (RSE in
    # round brackets, eta shrinkage in square brackets); bootstrap 96.6
    # (67.5-128) and 72.7 (47.4-96.8). The paper does not say whether the
    # percentage is omega (the square root of the OMEGA variance, the common
    # approximate-CV convention) or the exact log-normal coefficient of
    # variation sqrt(exp(omega^2) - 1). The two readings differ materially at
    # this magnitude: omega = 0.919 versus omega = sqrt(log(1 + 0.919^2)) =
    # 0.783.
    #
    # The Results section discriminates them. It reports, from 1000 simulations
    # of a 0.03 mg/(kg*h) infusion at 72 h, the proportion of individual
    # simulated concentrations above 1000 ug/L (27.8%, 10.6% and 5.4% at 0.5,
    # 1.25 and 2.5 kg) and below 200 ug/L (5.7%, 19.6% and 37.6%), plus pooled
    # figures across all nine simulated weights of 11.6% and 23.2%. Reproducing
    # that simulation with omega read directly as the printed percentage gives
    # 27.8 / 11.2 / 4.4 above and 6.8 / 21.0 / 38.7 below, pooling to 11.7% and
    # 24.1% -- agreeing with all eight published figures to within Monte Carlo
    # error. The exact-CV reading gives 26.1 / 8.3 / 2.5 and 3.9 / 17.0 / 36.6,
    # pooling to 9.3% and 20.9%, which misses the 2.5 kg tail by a factor of
    # two. The printed percentages are therefore omega x 100, so
    #   omega^2(CL) = 0.919^2 = 0.844561
    #   omega^2(V)  = 0.672^2 = 0.451584
    # The same simulation also confirms that the published proportions are
    # individual predictions WITHOUT residual error added. Table 2 reports no
    # correlation between the two etas, so they are entered as independent
    # diagonal elements. The full reproduction is the omega-scale gate in the
    # validation vignette.
    etalcl ~ 0.844561 # Table 2, row 'On CL [%]' = 91.9 (RSE 15%, shrinkage 8%); omega = 0.919
    etalvc ~ 0.451584 # Table 2, row 'On V [%]' = 67.2 (RSE 16%, shrinkage 17%); omega = 0.672

    # Combined proportional-plus-additive residual error, Table 2 'Residual
    # variability'. Proportional 33.8% (RSE 16%), bootstrap 34.0
    # (28.6-39.4) -- entered as the linear-scale fraction 0.338. Additive
    # 0.218 ug/L (RSE 56%), bootstrap 0.283 (0.0800-0.585) -- entered on the
    # ug/L concentration scale declared in `units`. The additive term is far
    # below both assay lower limits of quantification (1 ug/L for study 2,
    # 4 ug/L for study 1) and contributes negligibly over the observed
    # concentration range.
    propSd <- 0.338
    label("Proportional residual error (fraction)")
    addSd <- 0.218
    label("Additive residual error (ug/L)")
  })
  model({
    # Individual parameters. Both covariate relationships are centred on the
    # 1.12 kg median actual body weight given in the Voller 2019 Table 2
    # footnote:
    #   CL_i = CL_p * (WT_i / 1.12)^theta_WT   with theta_WT = 1.69 estimated
    #   V_i  = V_p  * (WT_i / 1.12)            i.e. exponent fixed at 1
    cl <- exp(lcl + etalcl) * (WT / 1.12)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 1.12)^e_wt_vc

    kel <- cl / vc

    # One-compartment model with intravenous input only. Results: 'A
    # 1-compartment model described the available data best.' Both studies
    # administered midazolam intravenously -- as bolus or short infusion for
    # stressful procedures and as continuous infusion for sedation -- so doses
    # target `central` directly and there is no depot state and no
    # bioavailability term.
    d/dt(central) <- -kel * central

    # Dose in ug and volume in L give ug/L, the concentration unit used
    # throughout Voller 2019 (target concentration 400 ug/L; assay lower limits
    # of quantification 1 and 4 ug/L). ug/L is numerically equal to ng/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
