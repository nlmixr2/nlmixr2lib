Truong_2025_ivacaftor <- function() {
  description <- "One-compartment population PK model for oral ivacaftor in children with cystic fibrosis aged 2-18 years carrying at least one F508del allele (Truong 2025, MODUL-CF). Absorption lag time and first-order absorption, first-order elimination, with allometric scaling of CL/F (exponent 0.75 fixed) and V/F (exponent 1 fixed) on body weight normalised to 70 kg. The lag time and absorption rate constant were not estimable from the sparse therapeutic-drug-monitoring design and were fixed to previously published values. Between-subject variability was estimable on apparent clearance only; residual variability is combined additive plus proportional, the only one of the three drugs to need an additive term. One of three independent per-drug models the paper reports for the elexacaftor/tezacaftor/ivacaftor combination; ivacaftor is dosed every 12 h whereas the two correctors are dosed once daily."
  reference   <- "Truong NH, Benaboud S, Bouazza N, et al. Elexacaftor/Tezacaftor/Ivacaftor Population Pharmacokinetics in Pediatric Patients With Cystic Fibrosis. Clin Transl Sci. 2025;18(5):e70245. doi:10.1111/cts.70245"
  vignette    <- "Truong_2025_elexacaftor_tezacaftor_ivacaftor"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Truong 2025 Methods 2.3 (plasma
  # concentrations assayed by LC-MS/MS, mg/L) and 2.5 (oral mg doses).
  compartmentData <- list(
    depot   = list(analyte = "ivacaftor", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ivacaftor", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Body weight enters as classical allometric scaling on both CL/F and V/F with the power exponents fixed to the theoretical values 0.75 and 1 respectively (Truong 2025 Methods 2.4, final paragraph, and Results 3.2). Weight was compared against BMI and BSA as size descriptors and 'the latter showed the most significant effect' (Results 3.2); adding it reduced the BIC by 35.2 for ivacaftor - the smallest of the three drugs - and reduced the between-subject variability on CL/F from 0.48 to 0.37. The reference weight of 70 kg is the normalisation the paper itself uses when it tabulates CL/F in L/h/70kg and V/F in L/70kg (Table 2) and when it compares the standardised estimates against adult values (Discussion paragraph 1). The paper does not state whether the weight column was time-varying across the 2022-2024 therapeutic-drug-monitoring period or fixed at baseline; Table 1 reports baseline weight only.",
      source_name        = "WT"
    )
  )

  # Covariates the paper screened but did not retain in the final model.
  # Documented for provenance only; they are deliberately absent from
  # model() and so must not live in covariateData.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a main covariate and not retained (Truong 2025 Methods 2.4 and Results 3.2: 'None of the other covariates including AGE, SEX, and GAL were related to PK parameters'). Age-related maturation functions were additionally evaluated on CL/F and on bioavailability and 'did not enhance the fit', which the authors attribute to CYP3A4/5 reaching 80% of adult activity by age two (Discussion paragraph 2)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as covariate SEX and not retained (Truong 2025 Results 3.2). The cohort was 66.7% male (Results 3.1)."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as an alternative size descriptor against body weight and body surface area; body weight showed the most significant effect and BMI was not retained (Truong 2025 Results 3.2)."
    ),
    BSA = list(
      description = "Body surface area at baseline.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened as an alternative size descriptor against body weight and BMI; body weight showed the most significant effect and BSA was not retained (Truong 2025 Results 3.2). The BSA computation formula is not stated in the paper."
    ),
    FORM_ETI_GRANULES = list(
      description = "Elexacaftor/tezacaftor/ivacaftor formulation indicator (1 = oral granule packet, 0 = tablet).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as covariate GAL (type of ETI formulation, granules or tablet) and not retained (Truong 2025 Methods 2.4 and Results 3.2). Recorded here under the canonical FORM_<drug>_<formulation> family name for provenance; because it is not retained it is not registered in inst/references/covariate-columns.md."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96,
    n_studies      = 1,
    age_range      = "2-18 years",
    age_strata     = "47 children <6 years (median 4.0 years, IQR 3.2-5.0); 19 children 6 to <12 years (median 8.8, IQR 7.0-10.4); 30 children >=12 years (median 15.8, IQR 14.6-16.6)",
    weight_range   = "Median weight by age stratum (Truong 2025 Table 1): 16.0 kg (IQR 14.0-17.5) for <6 years; 25.0 kg (IQR 21.6-31.2) for 6 to <12 years; 51.1 kg (IQR 44.5-57.5) for >=12 years. Overall cohort median 21 kg (IQR 16-43) per Table S1.",
    weight_median  = "21 kg (IQR 16-43), overall cohort (Truong 2025 Table S1)",
    sex_female_pct = 33.3,
    race_ethnicity = "Not reported. Cystic fibrosis was 'initially identified in the Caucasian population' (Introduction), but the paper tabulates no race or ethnicity data and did not screen either as a covariate.",
    disease_state  = "Cystic fibrosis with at least one F508del CFTR allele, receiving elexacaftor/tezacaftor/ivacaftor per marketing authorisation and undergoing therapeutic drug monitoring. None of the children had renal failure, liver disease, or coadministration of strong CYP3A inducers or inhibitors (Discussion, limitations paragraph).",
    dose_range     = "Weight-banded label dosing of ivacaftor (Truong 2025 Methods 2.5): 60 mg each morning and 59.5 mg each evening for <14 kg; 75 mg every 12 h for 14-30 kg; 150 mg every 12 h for >30 kg or age >12 years.",
    regions        = "France; 20 pediatric university hospital centres.",
    n_observations = "150 ivacaftor plasma concentrations. Median 1 sample per child per compound (range 1-5). Exactly one sample across all three compounds fell below the limit of quantification, an ivacaftor measurement from a 3-year-old, handled as left-censored.",
    co_medication  = "Concomitant elexacaftor and tezacaftor as part of the fixed triple combination. Strong CYP3A inducers and inhibitors were absent from the cohort.",
    notes          = "Ancillary pharmacokinetic sub-study of the MODUL-CF French prospective multicentre cohort (EUDRACT 2018-002624-16, NCT03894657), enrolling March 2022 to March 2024. Sparse therapeutic-drug-monitoring design: a trough just before intake and/or a peak roughly 4 h after intake. Estimation by MCMC-SAEM in Monolix 2023R1; concentrations below the limit of quantification were handled as left-censored. Assay validated over 0.0525-14 mg/L for ivacaftor. Model qualified by prediction-corrected VPC (Figure 1) and NPDE (Figure S2). Dietary fat intake was not recorded, so the label's food effect on ivacaftor bioavailability (AUC 2.5- to 4-fold higher with a moderate-fat meal) is not represented in this model; all parameters are apparent (X/F) under real-world dosing conditions."
  )

  ini({
    # Structural parameters - final-model Monolix estimates (Truong 2025
    # Table 2, "Ivacaftor" column). All clearance and volume terms are
    # apparent (X/F): the study is oral-only with no intravenous reference
    # arm, so bioavailability is not separately identifiable.

    # Absorption parameters were NOT estimated. Truong 2025 Results 3.2:
    # "Given the limited number of plasma samples available during the
    # absorption phase in this study, the absorption parameters (i.e., lag
    # time [Tlag] and/or absorption constant [Ka]) were set to values
    # reported in the literature for the three molecules [27]." Both values
    # are marked "(fix)" in Table 2 and are traced to that table.
    #
    # PROVENANCE CAVEAT: the paper cites reference [27] - Tsai A, Wu SP,
    # Haseltine E, et al. Pulm Ther. 2020;6(2):275-286
    # (doi:10.1007/s41030-020-00124-7) - collectively for all three
    # molecules, and that source does confirm the elexacaftor
    # (ka 0.59 /h, Tlag 2.17 h) and tezacaftor (ka 0.894 /h) values used in
    # the sibling model files. It does NOT, however, publish a first-order
    # ka or Tlag for ivacaftor: its ivacaftor compound file uses an ADAM
    # absorption model with the input type "Predicted" from a Caco-2
    # permeability of 11.9e-6 cm/s, and its text states only that "the
    # absorption rate constant (ka) was predicted by the Simcyp" simulator.
    # The specific ivacaftor values 0.42 /h and 0.70 h are therefore
    # traceable to Truong 2025 Table 2 itself and to no upstream printed
    # table. They are used as printed.
    ltlag <- fixed(log(0.70)); label("Absorption lag time (h)")                    # Truong 2025 Table 2: Tlag = 0.70 (fix)
    lka   <- fixed(log(0.42)); label("First-order absorption rate constant (1/h)")   # Truong 2025 Table 2: Ka = 0.42 (fix)

    lcl <- log(13.4); label("Apparent oral clearance at 70 kg (L/h)")              # Truong 2025 Table 2: CL/F = 13.4 L/h/70kg (4.8% RSE)
    lvc <- log(183);  label("Apparent volume of distribution at 70 kg (L)")         # Truong 2025 Table 2: V/F = 183 L/70kg (12.6% RSE)

    # Allometric exponents on body weight, held fixed at the theoretical
    # values by the source authors. Truong 2025 Results 3.2: "According to
    # the allometric rule, the power exponents for bodyweight effect were
    # fixed to 0.75 and 1 for the apparent clearance and the apparent volume
    # of distribution parameters."
    e_wt_cl <- fixed(0.75); label("Allometric WT exponent on CL/F (unitless)")           # Truong 2025 Methods 2.4 and Results 3.2: exponent fixed to 0.75
    e_wt_vc <- fixed(1);    label("Allometric WT exponent on V/F (unitless)")            # Truong 2025 Methods 2.4 and Results 3.2: exponent fixed to 1

    # IIV. Truong 2025 Methods 2.4: "The between-subject variabilities (BSV)
    # were ascribed to an exponential model." The Table 2 footnote defines
    # omega as the "interindividual variability estimate expressed as
    # standard deviation", so the tabulated 0.37 is a log-scale SD and ini()
    # takes its square. Only CL/F carried an estimable BSV (Results 3.2);
    # no BSV on V/F is reported, so none is encoded.
    etalcl ~ 0.1369   # 0.37^2;  Truong 2025 Table 2: omega_CL = 0.37 (12.4% RSE)

    # Residual error. Truong 2025 Results 3.2: "a combined residual error
    # model for ivacaftor" - the only one of the three drugs to carry an
    # additive term. Monolix distinguishes combined1 (sd = a + b*f) from
    # combined2 (sd = sqrt(a^2 + (b*f)^2)) and the paper names neither, so
    # the nlmixr2 default add() + prop() form, which is combined2, is used
    # here. See the vignette's Assumptions and deviations section.
    addSd  <- 0.131; label("Additive residual error (mg/L)")                             # Truong 2025 Table 2: sigma_additive = 0.131 mg/L (43.0% RSE)
    propSd <- 0.278; label("Proportional residual error (fraction)")                     # Truong 2025 Table 2: sigma_proportional = 0.278 (12.5% RSE)
  })

  model({
    # Reference body weight for allometric scaling. Truong 2025 Table 2
    # tabulates CL/F in L/h/70kg and V/F in L/70kg, and the Discussion
    # compares "the main PK parameters of ETI, when standardized to a 70 kg
    # individual" against adult values.
    ref_wt <- 70

    # Individual parameters.
    tlag <- exp(ltlag)
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc   <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    # Micro-constant.
    kel <- cl / vc

    # ODE system: one compartment with lagged first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot)   <- tlag

    # Observation. Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
