Dickinson_2021_dolutegravir_motherinfant <- function() {
  description <- "Mother-to-infant population PK model for dolutegravir (DolPHIN-1), extending the maternal model of modellib('Dickinson_2021_dolutegravir') with a one-compartment breastfed-infant model fitted sequentially with the maternal parameters fixed. At delivery the infant is born with the transplacental amount on board, equal to the predicted umbilical-cord concentration times the mother's apparent central volume; thereafter the infant receives dolutegravir by first-order input from the mother's breast-milk compartment and eliminates it by a first-order infant elimination rate constant. The time-varying PREG indicator marks delivery: while PREG = 1 the infant state tracks the transplacental amount, and from the first PREG = 0 record the infant ODE runs."
  reference <- paste(
    "Dickinson L, Walimbwa S, Singh Y, Kaboggoza J, Kintu K, Sihlangu M,",
    "Coombs JA, Malaba TR, Byamugisha J, Pertinez H, Amara A, Gini J, Else L,",
    "Heiberg C, Hodel EM, Reynolds H, Myer L, Waitt C, Khoo S, Lamorde M,",
    "Orrell C; DolPHIN-1 Study Group.",
    "Infant exposure to dolutegravir through placental and breast milk",
    "transfer: a population pharmacokinetic analysis of DolPHIN-1.",
    "Clin Infect Dis. 2021;73(5):e1200-e1207. doi:10.1093/cid/ciaa1861.",
    "Structure from the Results text, Figure 3 and the Supplementary",
    "Material; parameter values from Table 2 and its footnote b.",
    "Maternal layer shared with modellib('Dickinson_2021_dolutegravir').",
    sep = " "
  )
  vignette <- "Dickinson_2021_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # See Dickinson_2021_dolutegravir.R: the non-depleting cord-blood state
  # follows the Fauchet_2015_lopinavir_placental.R precedent.
  paper_specific_compartments <- c("fetal")

  covariateData <- list(
    PREG = list(
      description = "Pregnancy status of the mother, 1 = pregnant (before delivery), 0 = postpartum (from delivery onward)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (postpartum)",
      notes = "Time-varying. Used here only as the DELIVERY SWITCH of the dyad, not as a covariate effect: pregnancy had no retained effect on any maternal parameter (Discussion). While PREG = 1 the infant state is driven so that it equals the transplacental amount fetal * vc exactly; from delivery (the first PREG = 0 record) that amount becomes the infant's initial condition -- the paper's 'initial dolutegravir dose ... estimated by conversion of the cord concentration at time of delivery to an amount in milligrams by multiplication of the concentration and maternal central volume of distribution' (Supplementary Material) -- and the infant model starts. Set PREG = 1 on every record before the delivery time and PREG = 0 on every record at and after it. To simulate an infant from a known transplacental amount instead, pass PREG = 0 throughout and dose that amount into infant_central at the delivery time.",
      source_name = "delivery time (time of membrane rupture, or 15 min before cord sampling, when delivery time was missing)"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for between-occasion variability on maternal CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "As in Dickinson_2021_dolutegravir.R: a single interoccasion variability of 21.0% on CL/F (Table 2 footnote b) with occasions undefined in the paper; three slots are encoded for the third-trimester, delivery and postpartum visits, and occasions 2 and 3 fix their variance to the occasion-1 value. For a single-occasion record pass OCC = 1.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT_BIRTH = list(
      description = "Infant birth weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on the infant elimination rate constant but not retained (Results: 'No infant covariates were associated with kINF'). Cohort median 3.3 kg (2.5-4.3), Table 1.",
      source_name = "birth weight"
    ),
    BSA = list(
      description = "Infant body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on kINF but not retained. Cohort median 0.22 m^2 (0.18-0.25), Table 1.",
      source_name = "BSA"
    ),
    PNA = list(
      description = "Infant postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on kINF but not retained. Cohort median 7 days (3-18), Table 1; the canonical PNA is in months.",
      source_name = "postnatal age"
    ),
    GA = list(
      description = "Infant gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on kINF but not retained. Cohort median 39 weeks (35-43), Table 1.",
      source_name = "gestational age"
    ),
    PAGE = list(
      description = "Infant postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on kINF but not retained. Cohort median 40 weeks (36-44), Table 1.",
      source_name = "postmenstrual age"
    ),
    SEXF = list(
      description = "Infant sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on kINF but not retained. 5 of 22 infants (23%) female, Table 1.",
      source_name = "sex"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    fetal = list(analyte = "dolutegravir", units = "mg/L", specimen = "whole blood", verified = FALSE),
    milk = list(analyte = "dolutegravir", units = "mg", specimen = "milk", verified = TRUE),
    infant_central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 22,
    n_studies = 1,
    age_range = "Infants: postnatal age 3-18 days; mothers: 19-42 years",
    age_median = "Infants: 7 days; mothers: 27 years",
    weight_range = "Infant birth weight 2.5-4.3 kg; maternal weight 44-160 kg",
    weight_median = "Infant birth weight 3.3 kg; maternal weight 67 kg",
    sex_female_pct = 23,
    disease_state = "Breastfed, HIV-exposed infants of women living with HIV who started dolutegravir late in pregnancy (DolPHIN-1); exposed in utero and through breast milk",
    dose_range = "No direct infant dose. Mothers: dolutegravir 50 mg orally once daily from the third trimester until a switch to efavirenz-based therapy 2-18 days postpartum. Transplacental infant 'dose' at delivery 39.9 mg (range 15.5-59.0 mg), i.e. 12.5 mg/kg birth weight (5.0-19.6 mg/kg)",
    regions = "Uganda (10 infants) and South Africa (12 infants)",
    n_observations = "65 infant plasma concentrations",
    mother_partner = "Maternal layer identical to Dickinson_2021_dolutegravir (n = 28 mothers). Maternal individual parameters were fixed to their Bayesian estimates while the infant parameters were estimated (sequential approach, Supplementary Material).",
    notes = "DolPHIN-1 (NCT02245022). 22 of 27 infants had a recorded delivery time and could be modelled. 17 male (77%), 5 female (23%). Infant demographics from Table 1. 1 of 65 infant samples was BLQ and was included as LLQ/2."
  )

  ini({
    # ---- Maternal layer: identical to Dickinson_2021_dolutegravir ----
    lcl <- log(1.50)
    label("Apparent maternal clearance CL/F (L/h)") # Table 2 'CL/F' = 1.50 L/h (RSE 1.9%)
    lvc <- log(24.6)
    label("Apparent maternal central volume Vc/F (L)") # Table 2 'Vc/F' = 24.6 L (RSE 4.7%)
    lq <- log(0.0138)
    label("Apparent maternal intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F' = 0.0138 L/h (RSE 3.7%)
    lvp <- log(2.01)
    label("Apparent maternal peripheral volume Vp/F (L)") # Table 2 'Vp/F' = 2.01 L (RSE 16.7%)
    lka <- log(0.75)
    label("Absorption rate constant ka (1/h)") # Table 2 'ka' = 0.75 1/h (RSE 4.9%)
    lk_central_fetal <- log(2.81)
    label("Mother-to-fetus transfer rate constant kM-F (1/h)") # Table 2 'kM-F' = 2.81 1/h (RSE 3.6%)
    lk_fetal_central <- log(2.20)
    label("Fetus-to-mother transfer rate constant kF-M (1/h)") # Table 2 'kF-M' = 2.20 1/h (RSE 40.9%)
    lk_central_milk <- log(0.0027)
    label("Mother-to-breast-milk transfer rate constant kM-BM (1/h)") # Table 2 'kM-BM' = 0.0027 1/h (RSE 32.3%)
    lk_milk_central <- log(16.3)
    label("Breast-milk-to-mother transfer rate constant kBM-M (1/h)") # Table 2 'kBM-M' = 16.3 1/h (RSE 4.6%)
    lvmilk <- fixed(log(0.125))
    label("Breast-milk compartment volume VBM (L)") # Table 2 'VBM' = 0.125 L (Fixed)

    # ---- Infant layer (Table 2, 'Infant (n = 22)') ----
    lkmilkinf <- log(3.22)
    label("Breast-milk-to-infant transfer rate constant kBM-INF (1/h)") # Table 2 'kBM-INF' = 3.22 1/h (RSE 16.2%)
    lkel_infant <- log(0.0162)
    label("Infant elimination rate constant kINF (1/h)") # Table 2 'kINF' = 0.0162 1/h (RSE 6.0%)
    lvc_infant <- log(30.1)
    label("Infant apparent volume of distribution VINF/F (L)") # Table 2 'VINF/F' = 30.1 L (RSE 21.3%)

    # ---- IIV / IOV (Table 2 footnote b), omega^2 = log(1 + CV^2) ----
    etalcl + etalvc ~ c(0.020240, 0.025647, 0.041965) # Table 2 footnote b: IIV CL/F 14.3%, Vc/F 20.7%; Table 2 correlation 0.88
    etaiov_cl_1 ~ 0.043155 # Table 2 footnote b: IOV CL/F 21.0%; log(1 + 0.210^2)
    etaiov_cl_2 ~ fix(0.043155) # SAME-equivalent: equal to the occasion-1 variance
    etaiov_cl_3 ~ fix(0.043155) # SAME-equivalent: equal to the occasion-1 variance
    etalkel_infant ~ 0.174034 # Table 2 footnote b: IIV kINF 43.6% (RSE 61.6%); log(1 + 0.436^2)

    # ---- Residual error (Table 2 'Residual error, %'), proportional ----
    propSd <- 0.365
    label("Proportional residual error, maternal plasma (fraction)") # Table 2 'Maternal plasma' = 36.5%
    propSd_Cfetal <- 0.330
    label("Proportional residual error, umbilical cord (fraction)") # Table 2 'Umbilical cord' = 33.0%
    propSd_Cmilk <- 0.584
    label("Proportional residual error, breast milk (fraction)") # Table 2 'Breast milk' = 58.4%
    propSd_Cinfant <- 0.344
    label("Proportional residual error, infant plasma (fraction)") # Table 2 'Infant' = 34.4% (RSE 22.1%)
  })

  model({
    # 1. Occasion multiplexing of the between-occasion eta on maternal CL/F.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # 2. Individual parameters (maternal, then infant).
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    k_central_fetal <- exp(lk_central_fetal)
    k_fetal_central <- exp(lk_fetal_central)
    k_central_milk <- exp(lk_central_milk)
    k_milk_central <- exp(lk_milk_central)
    vmilk <- exp(lvmilk)

    kmilkinf <- exp(lkmilkinf)
    kel_infant <- exp(lkel_infant + etalkel_infant)
    vc_infant <- exp(lvc_infant)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # 4. Maternal ODEs (as in Dickinson_2021_dolutegravir.R).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1 - k_central_milk * central + k_milk_central * milk
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(fetal) <- k_central_fetal * Cc - k_fetal_central * fetal
    d/dt(milk) <- k_central_milk * central - k_milk_central * milk

    # 5. Infant ODE. Before delivery (PREG = 1) the infant state is driven by
    #    vc * d(fetal)/dt, so that it equals the transplacental amount
    #    Cfetal * vc at every instant (both start at 0 and vc is constant
    #    within a subject); at delivery that amount is the infant's initial
    #    'dose' (Supplementary Material). From delivery (PREG = 0) the infant
    #    receives first-order input kBM-INF times the breast-milk AMOUNT and
    #    eliminates at kINF (Figure 3). The breast-milk compartment is not
    #    drained by this input: the maternal model, including breast milk,
    #    was fixed to its individual estimates when the infant model was fitted
    #    (Figure 3 caption), and no emptying of the milk compartment was
    #    retained (Supplementary Material).
    d/dt(infant_central) <- PREG * vc *
      (k_central_fetal * Cc - k_fetal_central * fetal) +
      (1 - PREG) * (kmilkinf * milk - kel_infant * infant_central)

    # 6. Observations.
    Cfetal <- fetal
    Cmilk <- milk / vmilk
    Cinfant <- infant_central / vc_infant

    Cc ~ prop(propSd)
    Cfetal ~ prop(propSd_Cfetal)
    Cmilk ~ prop(propSd_Cmilk)
    Cinfant ~ prop(propSd_Cinfant)
  })
}
