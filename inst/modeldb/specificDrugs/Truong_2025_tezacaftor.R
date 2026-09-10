Truong_2025_tezacaftor <- function() {
  description <- "One-compartment population PK model for oral tezacaftor in children with cystic fibrosis aged 2-18 years carrying at least one F508del allele (Truong 2025, MODUL-CF). First-order absorption without a lag time and first-order elimination, with allometric scaling of CL/F (exponent 0.75 fixed) and V/F (exponent 1 fixed) on body weight normalised to 70 kg. The absorption rate constant was not estimable from the sparse therapeutic-drug-monitoring design and was fixed to a previously published value. Between-subject variability was estimable on apparent clearance only; residual variability is proportional. One of three independent per-drug models the paper reports for the elexacaftor/tezacaftor/ivacaftor combination; tezacaftor is the only one of the three that did not need an absorption lag time."
  reference   <- "Truong NH, Benaboud S, Bouazza N, et al. Elexacaftor/Tezacaftor/Ivacaftor Population Pharmacokinetics in Pediatric Patients With Cystic Fibrosis. Clin Transl Sci. 2025;18(5):e70245. doi:10.1111/cts.70245"
  vignette    <- "Truong_2025_elexacaftor_tezacaftor_ivacaftor"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Truong 2025 Methods 2.3 (plasma
  # concentrations assayed by LC-MS/MS, mg/L) and 2.5 (oral mg doses).
  compartmentData <- list(
    depot   = list(analyte = "tezacaftor", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tezacaftor", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Body weight enters as classical allometric scaling on both CL/F and V/F with the power exponents fixed to the theoretical values 0.75 and 1 respectively (Truong 2025 Methods 2.4, final paragraph, and Results 3.2). Weight was compared against BMI and BSA as size descriptors and 'the latter showed the most significant effect' (Results 3.2); adding it reduced the BIC by 84.2 for tezacaftor - the largest of the three drugs - and reduced the between-subject variability on CL/F from 0.39 to 0.26. The reference weight of 70 kg is the normalisation the paper itself uses when it tabulates CL/F in L/h/70kg and V/F in L/70kg (Table 2) and when it compares the standardised estimates against adult values (Discussion paragraph 1). The paper does not state whether the weight column was time-varying across the 2022-2024 therapeutic-drug-monitoring period or fixed at baseline; Table 1 reports baseline weight only.",
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
    dose_range     = "Weight-banded label dosing of tezacaftor once daily (Truong 2025 Methods 2.5): 40 mg for <14 kg; 50 mg for 14-30 kg; 100 mg for >30 kg or age >12 years.",
    regions        = "France; 20 pediatric university hospital centres.",
    n_observations = "150 tezacaftor plasma concentrations. Median 1 sample per child per compound (range 1-5). No tezacaftor sample fell below the limit of quantification.",
    co_medication  = "Concomitant elexacaftor and ivacaftor as part of the fixed triple combination. Strong CYP3A inducers and inhibitors were absent from the cohort.",
    notes          = "Ancillary pharmacokinetic sub-study of the MODUL-CF French prospective multicentre cohort (EUDRACT 2018-002624-16, NCT03894657), enrolling March 2022 to March 2024. Sparse therapeutic-drug-monitoring design: a trough just before intake and/or a peak roughly 4 h after intake. Estimation by MCMC-SAEM in Monolix 2023R1; concentrations below the limit of quantification were handled as left-censored. Assay validated over 0.075-20 mg/L for tezacaftor. Model qualified by prediction-corrected VPC (Figure 1) and NPDE (Figure S2). Dietary fat intake was not recorded; the label reports that food does not materially affect tezacaftor bioavailability (unlike elexacaftor and ivacaftor), so this omission is least consequential for tezacaftor. All parameters are apparent (X/F) under real-world dosing conditions."
  )

  ini({
    # Structural parameters - final-model Monolix estimates (Truong 2025
    # Table 2, "Tezacaftor" column). All clearance and volume terms are
    # apparent (X/F): the study is oral-only with no intravenous reference
    # arm, so bioavailability is not separately identifiable.

    # No absorption lag time for tezacaftor. Truong 2025 Results 3.2: "A
    # one-compartment model with first-order absorption and elimination
    # optimally described the pharmacokinetics of tezacaftor, while the
    # pharmacokinetics of elexacaftor and ivacaftor were best fitted by a
    # one-compartment model incorporating absorption lag time". Table 2
    # accordingly prints "-" in the tezacaftor Tlag cell.

    # Ka was NOT estimated. Truong 2025 Results 3.2: "the absorption
    # parameters (i.e., lag time [Tlag] and/or absorption constant [Ka])
    # were set to values reported in the literature for the three molecules
    # [27]." Reference [27] is Tsai A, Wu SP, Haseltine E, et al. Pulm Ther.
    # 2020;6(2):275-286 (doi:10.1007/s41030-020-00124-7), whose tezacaftor
    # compound file carries ka = 0.894 /h - the same value Truong 2025
    # prints. It is marked "(fix)" in Table 2.
    lka <- fixed(log(0.894)); label("First-order absorption rate constant (1/h)")     # Truong 2025 Table 2: Ka = 0.894 (fix)

    lcl <- log(1.32); label("Apparent oral clearance at 70 kg (L/h)")               # Truong 2025 Table 2: CL/F = 1.32 L/h/70kg (3.8% RSE)
    lvc <- log(29.3); label("Apparent volume of distribution at 70 kg (L)")          # Truong 2025 Table 2: V/F = 29.3 L/70kg (6.8% RSE)

    # Allometric exponents on body weight, held fixed at the theoretical
    # values by the source authors. Truong 2025 Results 3.2: "According to
    # the allometric rule, the power exponents for bodyweight effect were
    # fixed to 0.75 and 1 for the apparent clearance and the apparent volume
    # of distribution parameters."
    e_wt_cl <- fixed(0.75); label("Allometric WT exponent on CL/F (unitless)")            # Truong 2025 Methods 2.4 and Results 3.2: exponent fixed to 0.75
    e_wt_vc <- fixed(1);    label("Allometric WT exponent on V/F (unitless)")             # Truong 2025 Methods 2.4 and Results 3.2: exponent fixed to 1

    # IIV. Truong 2025 Methods 2.4: "The between-subject variabilities (BSV)
    # were ascribed to an exponential model." The Table 2 footnote defines
    # omega as the "interindividual variability estimate expressed as
    # standard deviation", so the tabulated 0.26 is a log-scale SD and ini()
    # takes its square. Only CL/F carried an estimable BSV (Results 3.2);
    # no BSV on V/F is reported, so none is encoded.
    etalcl ~ 0.0676   # 0.26^2;  Truong 2025 Table 2: omega_CL = 0.26 (14.8% RSE)

    # Residual error. Truong 2025 Results 3.2: "A proportional residual
    # error model was retained for elexacaftor/tezacaftor".
    propSd <- 0.302; label("Proportional residual error (fraction)")                      # Truong 2025 Table 2: sigma_proportional = 0.302 (8.5% RSE)
  })

  model({
    # Reference body weight for allometric scaling. Truong 2025 Table 2
    # tabulates CL/F in L/h/70kg and V/F in L/70kg, and the Discussion
    # compares "the main PK parameters of ETI, when standardized to a 70 kg
    # individual" against adult values.
    ref_wt <- 70

    # Individual parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    # Micro-constant.
    kel <- cl / vc

    # ODE system: one compartment with first-order absorption, no lag time.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
