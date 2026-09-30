Agema_2021_oxycodone <- function() {
  description <- paste(
    "Joint one-compartment population PK model for oral oxycodone and its",
    "sequential N-demethylated metabolites noroxycodone and noroxymorphone in",
    "hospitalised adults with cancer-related pain (Agema 2021). Immediate-",
    "release (IR) and extended-release (ER) tablets are dosed into two separate",
    "first-order absorption depots (depot = IR, depot2 = ER); oxycodone is",
    "eliminated by a fixed first-order route and converted by a first-order rate",
    "constant to noroxycodone, which is eliminated and further converted to",
    "noroxymorphone. The model is parameterised in micro-rate constants and",
    "apparent volumes, with IIV on the two formation/elimination constants of",
    "noroxycodone and on noroxymorphone elimination, and no covariates. DOSES",
    "ARE IN MMOL OF OXYCODONE and all three concentrations are in nmol/L, as in",
    "the authors' NONMEM dataset: 10 mg oxycodone (free base, MW 315.36 g/mol)",
    "is 0.03171 mmol."
  )
  reference <- paste(
    "Agema BC, Oosten AW, Sassen SDT, Rietdijk WJR, van der Rijt CCD, Koch BCP,",
    "Mathijssen RHJ, Koolen SLW. Population Pharmacokinetics of Oxycodone and",
    "Metabolites in Patients with Cancer-Related Pain. Cancers (Basel).",
    "2021;13(11):2768. doi:10.3390/cancers13112768"
  )
  vignette <- "Agema_2021_oxycodone"

  # Supplementary Document S2 (final-model control stream), $PK:
  # ';scaling from nmol/L (observations) to mmol (dosages)' with
  # S3 = V3/1000000, S4 = V4/1000000, S5 = V5/1000000. Dataset AMT is
  # therefore mmol of oxycodone and DV is nmol/L for every analyte. The
  # stream's LLOQ constants (0.6342, 3.318511 and 3.480561 nmol/L) equal the
  # assay LLOQs of Methods 2.3 (0.200, 1.00 and 1.00 ng/mL) divided by the
  # free-base molecular weights of oxycodone (315.36), noroxycodone (301.34)
  # and noroxymorphone (287.31 g/mol), which pins the concentration
  # conversion. The dose conversion (mg tablet -> mmol) is not printed.
  units <- list(time = "h", dosing = "mmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(
      analyte = "oxycodone (immediate-release tablet)",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte = "oxycodone (extended-release tablet)",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "oxycodone", units = "mmol", specimen = "plasma", verified = TRUE),
    central_noroxycod = list(analyte = "noroxycodone", units = "mmol", specimen = "plasma", verified = TRUE),
    central_noroxymor = list(analyte = "noroxymorphone", units = "mmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Results 3.2: 'None of the covariates significantly improved the model when
  # included and were therefore not incorporated into the model.' Methods 2.1
  # lists the baseline covariates collected; the control stream's $INPUT
  # carries HT WT AGE GNDR BMI ALBU CREAT CYP2D6 CYP3A4 UGT2B7. None is
  # referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (median-centred power / exponential models, Document S1); not retained. Cohort median 80.0 kg, range 46-135 (Table 1).",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Collected at baseline and carried in the dataset ($INPUT HT); not retained.",
      source_name = "HT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 27.7, range 18.9-42.6 (Table 1).",
      source_name = "BMI"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 62.5 years, range 39-81 (Table 1). The Discussion attributes the null result partly to the narrow age distribution.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a proportional (categorical) covariate; not retained. 12 of 28 (43%) female (Table 1). Dataset coding of GNDR is not printed.",
      source_name = "GNDR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 42.5 g/L, range 29-49 (Table 1).",
      source_name = "ALBU"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min",
      type = "continuous",
      notes = "Screened (dataset column CREAT); not retained. Cohort median eGFR 86.0 mL/min, range 37 to >90 (Table 1); only three patients below 60 mL/min.",
      source_name = "CREAT"
    ),
    CYP2D6_IM = list(
      description = "CYP2D6 intermediate-metabolizer phenotype indicator (reference: extensive metabolizer)",
      units = "(binary)",
      type = "binary",
      notes = "Screened with the genotype models of Document S1; not retained. 12 EM, 10 IM, 6 missing (Table 1); no poor or ultrarapid metabolizers. Discussion: no significant EM-vs-IM difference in the noroxycodone -> noroxymorphone conversion rate.",
      source_name = "CYP2D6"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 22 *1/*1, 6 *1/*22, 0 *22/*22 (Table 1). Discussion: no difference in the oxycodone -> noroxycodone conversion rate.",
      source_name = "CYP3A4"
    ),
    SNP_UGT2B7_RS7439366 = list(
      description = "UGT2B7*2 (rs7439366, H268Y) allele count or carrier indicator",
      units = "(count)",
      type = "count",
      notes = "Screened with exponential / dominant / recessive genotype models (Document S1); not retained. 6 wild type, 17 heterozygous, 5 variant (Table 1). Discussion reports a non-significant trend towards decreased oxycodone metabolism in T-allele carriers.",
      source_name = "UGT2B7"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28,
    n_observations = 1207,
    n_studies = 1,
    age_range = "39-81 years",
    age_median = "62.5 years",
    weight_range = "46-135 kg",
    weight_median = "80.0 kg",
    bmi_median = "27.7 kg/m^2 (range 18.9-42.6)",
    sex_female_pct = 43,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Hospitalised adults with moderate to severe nociceptive cancer-related",
      "pain (primary tumour urogenital 32%, breast 18%, GIST / soft-tissue",
      "sarcoma 14%, melanoma 11%, other 25%); WHO performance status 1-3; no",
      "Child-Pugh B or C liver dysfunction."
    ),
    renal_function = "eGFR median 86.0 mL/min (range 37 to >90); three patients below 60 mL/min",
    genotypes = "CYP2D6 12 EM / 10 IM / 6 missing; CYP3A4*22 22 *1/*1 / 6 *1/*22; UGT2B7*2 6 wild type / 17 heterozygous / 5 variant",
    dose_range = "Oral oxycodone ER tablets 5-100 mg twice daily (08:00 and 20:00) and IR tablets 5-30 mg as needed",
    regions = "The Netherlands (Erasmus MC Cancer Institute, Rotterdam)",
    notes = paste(
      "Observational study (MEC 09.332, NTR4369). 302 plasma samples from 28",
      "patients (one patient enrolled twice), each assayed for oxycodone,",
      "oxymorphone, noroxycodone and noroxymorphone (1207 measurements, 29.2%",
      "BLQ; M3 method). Oxymorphone was not modelled (179 of 309 samples BLQ).",
      "Samples twice daily at the ER dosing times plus, at most once daily, pre",
      "and 5, 15, 30 and 60 min after an IR tablet (Methods 2.2; Table 1)."
    )
  )

  ini({
    # ---- Structural parameters: Table 2 'Parameter Estimate' column. The
    # control stream (Document S2) names the NONMEM compartments 1 = ABSIR,
    # 2 = ABSER, 3 = CENTRAL (oxycodone), 4 = METAB1 (noroxycodone),
    # 5 = METAB3 (noroxymorphone); its $THETA block holds INITIAL estimates
    # (e.g. V3 566 vs the Table 2 final 619), so the Table 2 finals are
    # encoded, except for the fixed K30 (see below).
    lka_ir <- log(3.61); label("Absorption rate constant of the immediate-release tablet, Ka,IR (K13) (1/h)") # Table 2 Ka,IR = 3.61 1/h (RSE 39.3%)
    lka_er <- log(0.329); label("Absorption rate constant of the extended-release tablet, Ka,ER (K23) (1/h)") # Table 2 Ka,ER = 0.329 1/h (RSE 38.6%); the printed row label 'Ka, ER 23' carries the stream's K23 subscript
    lvc <- log(619); label("Apparent volume of oxycodone, V3/F (L)") # Table 2 V3/F = 619 L (RSE 9.2%)
    # K30 was fixed from an intermediate model (Results 3.2). Table 2 prints it
    # rounded to 0.012; the control stream's $THETA(3) '(0, 0.01224) FIX' is
    # the value the final run actually used, so it is encoded here.
    lkel <- fixed(log(0.01224)); label("Oxycodone elimination rate constant by routes other than noroxycodone formation, K30 (1/h)") # Table 2 K30 = 0.012 FIX; Document S2 $THETA(3) = 0.01224 FIX
    lkmet_noroxycod <- log(0.086); label("Oxycodone-to-noroxycodone conversion rate constant, K34 (1/h)") # Table 2 K34 = 0.086 1/h (RSE 10.9%)
    lvc_noroxycod <- fixed(log(16.3)); label("Apparent volume of noroxycodone, V4/F (L)") # Table 2 V4/F = 16.3 L FIX; Document S2 $THETA(7) = 16.3 FIX
    lkel_noroxycod <- log(3.28); label("Noroxycodone elimination rate constant, K40 (1/h)") # Table 2 K40 = 3.28 1/h (RSE 30.1%)
    lkmet_noroxymor <- log(1.36); label("Noroxycodone-to-noroxymorphone conversion rate constant, K45 (1/h)") # Table 2 K45 = 1.36 1/h (RSE 32.1%)
    lvc_noroxymor <- fixed(log(64.1)); label("Apparent volume of noroxymorphone, V5/F (L)") # Table 2 V5/F = 64.1 L FIX; Document S2 $THETA(12) = 64.1 FIX
    lkel_noroxymor <- log(1.97); label("Noroxymorphone elimination rate constant, K50 (1/h)") # Table 2 K50 = 1.97 1/h (RSE 35.9%)

    # ---- IIV: exponential etas (Methods 2.4; Document S2 'K34 = THETA(9) *
    # EXP(ETA(1))' etc.). Table 2 reports CV%, converted with
    # omega^2 = log(1 + CV^2). The stream estimates a full $OMEGA BLOCK(3), but
    # the paper prints no final covariances, so the block is encoded diagonal.
    etalkmet_noroxycod ~ 0.1199576 # Table 2 IIV K34 = 35.7 CV% (shrinkage 5.4%); log(1 + 0.357^2)
    etalkel_noroxycod ~ 0.6801476 # Table 2 IIV K40 = 98.7 CV% (shrinkage 3.4%); log(1 + 0.987^2)
    etalkel_noroxymor ~ 0.5113189 # Table 2 IIV K50 = 81.7 CV% (shrinkage 7.1%); log(1 + 0.817^2)

    # ---- Residual error. Document S2 $ERROR writes each analyte's SD as
    # W = IPRED * THETA_prop + THETA_add with $SIGMA 1 FIX, i.e. the
    # proportional and additive SDs are summed linearly (nlmixr2 combined1).
    # The oxycodone additive term was fixed to 0 (Results 3.2), leaving a
    # purely proportional error.
    propSd <- 0.397; label("Proportional residual error, oxycodone (fraction)") # Table 2 oxycodone proportional 39.7% (RSE 13.4%); additive 0 FIX
    propSd_noroxycod <- 0.167; label("Proportional residual error, noroxycodone (fraction)") # Table 2 noroxycodone proportional 16.7% (RSE 18.5%)
    addSd_noroxycod <- 3.34; label("Additive residual error, noroxycodone (nmol/L)") # Table 2 noroxycodone additive 3.34 nM (RSE 40.4%)
    propSd_noroxymor <- 0.156; label("Proportional residual error, noroxymorphone (fraction)") # Table 2 noroxymorphone proportional 15.6% (RSE 21.5%)
    addSd_noroxymor <- 1.09; label("Additive residual error, noroxymorphone (nmol/L)") # Table 2 noroxymorphone additive 1.09 nM (RSE 24.3%)
  })

  model({
    # 1. Individual parameters (Document S2 $PK).
    ka_ir <- exp(lka_ir)
    ka_er <- exp(lka_er)
    vc <- exp(lvc)
    kel <- exp(lkel)
    kmet_noroxycod <- exp(lkmet_noroxycod + etalkmet_noroxycod)
    vc_noroxycod <- exp(lvc_noroxycod)
    kel_noroxycod <- exp(lkel_noroxycod + etalkel_noroxycod)
    kmet_noroxymor <- exp(lkmet_noroxymor)
    vc_noroxymor <- exp(lvc_noroxymor)
    kel_noroxymor <- exp(lkel_noroxymor + etalkel_noroxymor)

    # 2. ODE system, Document S2 $DES (Figure 2). Amounts are mmol throughout;
    #    conversion moves amount one-for-one, so the metabolite volumes are
    #    apparent (V/F, per mmol of oxycodone dosed).
    d/dt(depot) <- -ka_ir * depot
    d/dt(depot2) <- -ka_er * depot2
    d/dt(central) <- ka_ir * depot + ka_er * depot2 - (kel + kmet_noroxycod) * central
    d/dt(central_noroxycod) <- kmet_noroxycod * central - (kel_noroxycod + kmet_noroxymor) * central_noroxycod
    d/dt(central_noroxymor) <- kmet_noroxymor * central_noroxycod - kel_noroxymor * central_noroxymor

    # 3. Observations: mmol / L * 1e6 = nmol/L (Document S2 S3 = V3/1000000).
    Cc <- 1e6 * central / vc
    Cc_noroxycod <- 1e6 * central_noroxycod / vc_noroxycod
    Cc_noroxymor <- 1e6 * central_noroxymor / vc_noroxymor

    Cc ~ prop(propSd)
    Cc_noroxycod ~ add(addSd_noroxycod) + prop(propSd_noroxycod) + combined1()
    Cc_noroxymor ~ add(addSd_noroxymor) + prop(propSd_noroxymor) + combined1()
  })
}
