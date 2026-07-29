Wilbaux_2014_ovarianCancer_ca125 <- function() {
  description <- "K-PD joint model of tumour size and CA-125 kinetics during platinum-based chemotherapy in recurrent ovarian cancer (CALYPSO phase III trial; carboplatin + paclitaxel or carboplatin + pegylated liposomal doxorubicin pooled)"
  reference <- "Wilbaux M, Henin E, Oza A, Colomban O, Pujade-Lauraine E, Freyer G, Tod M, You B. Prediction of tumour response induced by chemotherapy using modelling of CA-125 kinetics in recurrent ovarian cancer patients. Br J Cancer. 2014;110(6):1517-1524. doi:10.1038/bjc.2014.75"
  vignette <- "Wilbaux_2014_ovarianCancer_ca125"
  paper_specific_compartments <- c("ca125")
  paper_specific_etas <- c(
    "etalkprol", "etala50", "etalkreduc",
    "etalkprod1", "etalkprod2", "etalk2",
    "etalrbase_ca125"
  )
  paper_specific_residual_sds <- c("addSd_ca125", "lambda_ca125")
  units <- list(
    time = "day",
    dosing = "a.u. (arbitrary K-PD dose, 1 unit per chemotherapy cycle)",
    concentration = "U/mL for CA-125; mm for tumour size"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at baseline",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise forward-selection / backward-elimination; not retained (Methods + Results).",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "Weight"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "Height"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "Body surface area"
    ),
    CREAT = list(
      description = "Serum creatinine at baseline",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "Creatinine"
    ),
    PFS_PRIOR = list(
      description = "Progression-free survival following the first prior chemotherapy line",
      units = "months",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "PFS 1st chemo"
    ),
    THERAPY_FREE_INT = list(
      description = "Patient therapy-free interval before study entry",
      units = "months",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "Patient therapy free interval"
    ),
    TRT = list(
      description = "Chemotherapy regimen indicator (CP = carboplatin + paclitaxel; CD = carboplatin + pegylated liposomal doxorubicin)",
      units = "(binary)",
      type = "binary",
      reference_category = "CP",
      notes = "Separate K-PD kinetics were initially fit for CP and CD but the estimates were not significantly different, so a single set of drug-kinetic parameters was retained (Methods, Drug kinetics paragraph). The treatment indicator was therefore not entered into the covariate selection step.",
      source_name = "Treatment"
    ),
    SURGERY_28D = list(
      description = "Any debulking surgery within 28 days of study entry",
      units = "(binary)",
      type = "binary",
      reference_category = "No",
      notes = "Screened; not retained.",
      source_name = "Any surgery within 28 days"
    ),
    FIGO_STAGE = list(
      description = "FIGO ovarian-cancer stage at diagnosis (I, II, or III)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened; not retained.",
      source_name = "FIGO stage"
    ),
    TUMOUR_SITE = list(
      description = "Primary tumour anatomical site (Ovary, Peritoneal, or Fallopian)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened; not retained.",
      source_name = "Primary tumour site"
    ),
    WBC_HIGH = list(
      description = "Elevated white blood cell count indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "No",
      notes = "Screened; not retained.",
      source_name = "Elevated white blood cells"
    ),
    ASCITES = list(
      description = "Ascites involvement indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "No",
      notes = "Screened; not retained.",
      source_name = "Ascite involvement"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 535,
    n_studies = 1,
    age_range = "27-82 years",
    age_median = "61 years",
    weight_range = "41-150 kg",
    weight_median = "69 kg",
    sex_female_pct = 100,
    disease_state = "Platinum-sensitive recurrent ovarian, peritoneal, or fallopian cancer (ROC); platinum-sensitive defined by the CALYPSO inclusion criteria",
    dose_range = "Arbitrary K-PD dose = 1 a.u. per chemotherapy cycle (carboplatin + paclitaxel or carboplatin + pegylated liposomal doxorubicin; cycles every 3-4 weeks per the CALYPSO protocol)",
    regions = "Multicentre international (CALYPSO phase III trial)",
    notes = paste(
      "Modelling time horizon limited to the first 513 days after randomisation",
      "(50% of patients had discontinued by this day; Methods, Data management).",
      "Learning data set: 357 patients used to estimate parameters; validation",
      "data set: 178 patients used for advanced internal evaluation. 297",
      "additional patients with non-measurable disease (and therefore missing",
      "tumour-size data) were excluded from model building but used to predict",
      "latent tumour response from CA-125 alone. Patient characteristics are in",
      "Table 1 of Wilbaux 2014."
    )
  )

  ini({
    # Drug kinetics (K-PD): a single first-order rate constant K drives both
    # decay of the first K-PD compartment (depot) and decay of the transit
    # compartment that lags the drug effect (transit1). Reported as one
    # parameter K because separate kinetics for CP vs CD were not significantly
    # different (Methods, Drug kinetics paragraph).
    lkel        <- log(0.019) ; label("K-PD treatment kinetic rate constant K (1/day)")                       # Table 2 (K = 0.019)

    # Tumour dynamics: constant proliferation rate KPROL inhibited saturably
    # by the K-PD drug amount in transit1, with a linear size-proportional
    # decrease at rate KREDUC. Drug acts as inhibitor of growth in a saturable
    # manner per Methods and Figure 1.
    lkprol      <- log(0.869) ; label("Tumour growth rate KPROL (mm/day)")                                    # Table 2 (KPROL = 0.869)
    la50        <- log(0.162) ; label("K-PD drug amount producing 50% maximum growth inhibition A50 (a.u.)")  # Table 2 (A50 = 0.162)
    lkreduc     <- log(0.013) ; label("Tumour size decrease rate constant KREDUC (1/day)")                    # Table 2 (KREDUC = 0.013)

    # CA-125 indirect response: basal production KPROD1 plus stationary-tumour
    # production KPROD2 modulated by exp(K2 * VARTS) so growing tumours
    # exponentially upregulate CA-125 shedding and shrinking tumours
    # exponentially downregulate it (production stays strictly positive);
    # first-order elimination at KELIM. See Methods and Figure 1 caption
    # and the explicit equation block "dCA125/dt = KPROD1 + KPROD2 *
    # exp(K2*VARTS) - KELIM*CA125".
    lkprod1     <- log(0.452) ; label("CA-125 basal production rate KPROD1 (U/mL per day)")                   # Table 2 (KPROD1 = 0.452)
    lkprod2     <- log(0.615) ; label("CA-125 stationary-tumour production rate KPROD2 (U/mL per day)")        # Table 2 (KPROD2 = 0.615)
    lk2         <- log(21.4)  ; label("Exponential factor linking KPROD2 to dTS/dt K2 (day/mm)")               # Table 2 (K2 = 21.4)
    lkout       <- log(0.037) ; label("CA-125 elimination rate constant KELIM (1/day)")                       # Table 2 (KELIM = 0.037)

    # Estimated baseline initial conditions; logit-constrained TS0 to be less
    # than the no-drug steady state KPROL/KREDUC = 66.8 mm so the tumour can
    # grow (Methods, Supplementary Material 3).
    lrbase_tumor_size <- log(57.8) ; label("Baseline tumour size TS0 (mm)")                                   # Table 2 (TS0 = 57.8)
    lrbase_ca125      <- log(167)  ; label("Baseline CA-125 concentration CA0 (U/mL)")                        # Table 2 (CA0 = 167)

    # Inter-individual variability. Reported in Table 2 as CV%; converted to
    # omega^2 on the log scale via omega^2 = log(CV^2 + 1).
    etalkel              ~ log(0.721^2 + 1)   # Table 2 (K IIV 72.1 %CV)
    etalkprol            ~ log(1.536^2 + 1)   # Table 2 (KPROL IIV 153.6 %CV)
    etala50              ~ log(1.670^2 + 1)   # Table 2 (A50 IIV 167.0 %CV)
    etalkreduc           ~ log(1.404^2 + 1)   # Table 2 (KREDUC IIV 140.4 %CV)
    etalkprod1           ~ log(0.870^2 + 1)   # Table 2 (KPROD1 IIV 87.0 %CV)
    etalkprod2           ~ log(2.258^2 + 1)   # Table 2 (KPROD2 IIV 225.8 %CV)
    etalk2               ~ log(1.086^2 + 1)   # Table 2 (K2 IIV 108.6 %CV)
    etalkout             ~ log(0.636^2 + 1)   # Table 2 (KELIM IIV 63.6 %CV)
    etalrbase_tumor_size ~ log(0.756^2 + 1)   # Table 2 (TS0 IIV 75.6 %CV)
    etalrbase_ca125      ~ log(1.127^2 + 1)   # Table 2 (CA0 IIV 112.7 %CV)

    # Residual error.
    # Tumour size was log-transformed during fitting and an additive error
    # was applied on the log scale, i.e. log-normal residual on the linear
    # scale. The paper also reports a fixed below-LOQ residual SD of 0.576
    # (LOQ = 10 mm, M5 method with BLQ set to LOQ/2 and SD anchored at
    # LOQ/4); BLQ handling is done in the data, the model itself uses the
    # above-LOQ residual.
    expSd_tumor_size <- 0.273         ; label("Tumour size log-normal residual SD (above LOQ, log scale)")     # Table 2 (above-LOQ residual = 0.273)
    # CA-125 was Box-Cox transformed (lambda fixed at the R car::powerTransform
    # estimate of -0.16) and the residual was additive on the transformed
    # scale; the boxCox(lambda) factor below applies the same transformation
    # to predicted and observed CA-125 when forming the likelihood.
    addSd_ca125  <- 0.165             ; label("CA-125 additive residual SD on Box-Cox(lambda) scale")          # Table 2 (CA-125 residual = 0.165)
    lambda_ca125 <- fixed(-0.16)      ; label("CA-125 Box-Cox transformation parameter lambda (unitless)")     # Methods, Data management (R car::powerTransform estimate)
  })

  model({
    # Individual K-PD and PD parameters
    kel    <- exp(lkel    + etalkel)
    kprol  <- exp(lkprol  + etalkprol)
    a50    <- exp(la50    + etala50)
    kreduc <- exp(lkreduc + etalkreduc)
    kprod1 <- exp(lkprod1 + etalkprod1)
    kprod2 <- exp(lkprod2 + etalkprod2)
    k2     <- exp(lk2     + etalk2)
    kout   <- exp(lkout   + etalkout)

    # Initial conditions for the dynamic states
    tumor_size(0) <- exp(lrbase_tumor_size + etalrbase_tumor_size)
    ca125(0)      <- exp(lrbase_ca125      + etalrbase_ca125)

    # K-PD chain. depot receives the per-cycle bolus dose of 1 a.u.; transit1
    # holds the delayed drug amount that drives tumour-growth inhibition.
    # Both decay first-order at the same rate K (Methods).
    d/dt(depot)    <- -kel * depot
    d/dt(transit1) <-  kel * depot - kel * transit1

    # Tumour dynamics. Drug saturably inhibits the constant growth rate
    # KPROL through A2 / (A50 + A2) (paper's A2 = transit1 here). Decrease
    # is linear in tumour size at rate KREDUC, independent of drug. Without
    # drug, dTS/dt = KPROL - KREDUC*TS reaches a steady state at KPROL/KREDUC
    # = 66.8 mm; TS0 = 57.8 mm sits below this ceiling by construction.
    drugInhibition <- a50 / (a50 + transit1)
    varts          <- kprol * drugInhibition - kreduc * tumor_size
    d/dt(tumor_size) <- varts

    # CA-125 indirect response. Basal (normal-tissue) and stationary-tumour
    # production are additive; the stationary-tumour rate is multiplied by
    # exp(K2 * VARTS) so growing tumours (VARTS > 0) exponentially
    # upregulate CA-125 shedding and shrinking tumours (VARTS < 0)
    # exponentially downregulate it. Elimination is first-order at KELIM.
    d/dt(ca125) <- kprod1 + kprod2 * exp(k2 * varts) - kout * ca125

    # Observations and error model
    tumor_size ~ lnorm(expSd_tumor_size)
    ca125      ~ add(addSd_ca125) + boxCox(lambda_ca125)
  })
}
