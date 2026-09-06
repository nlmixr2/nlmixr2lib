Clegg_2024_tixagevimab_cilgavimab <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order intramuscular",
    "absorption and linear elimination for AZD7442 (Evusheld), the",
    "co-formulated combination of the YTE-modified SARS-CoV-2-neutralising",
    "monoclonal antibodies tixagevimab and cilgavimab, in 4,940 adults from",
    "eight COVID-19 prophylaxis and treatment trials (Clegg 2024 Table 1).",
    "AZD7442 is the summed analyte: serum AZD7442 concentration is the sum of",
    "the measured tixagevimab and cilgavimab concentrations. Body weight is",
    "carried on CL, Vc, Q and Vp with allometric exponents fixed at 0.75 and 1",
    "and a 70 kg reference. Absorption rate carries male sex, age > 65 years,",
    "BMI >= 30 kg/m2 and diabetes; clearance carries diabetes; central volume",
    "carries Black race; and intramuscular bioavailability carries the",
    "injection site (thigh versus gluteal). IIV is a 3x3 block on ka, CL and",
    "Vc plus an independent eta on Vp, and the residual error is additive on",
    "the log scale with separate magnitudes for intramuscular and intravenous",
    "records. Companion analyte-specific models:",
    "modellib('Clegg_2024_tixagevimab') and modellib('Clegg_2024_cilgavimab').",
    sep = " "
  )
  reference <- paste(
    "Clegg LE, Stepanov O, Schmidt H, Tang W, Zhang H, Webber C, Cohen TS,",
    "Esser MT, Nagard M. Accelerating therapeutics development during a",
    "pandemic: population pharmacokinetics of the long-acting antibody",
    "combination AZD7442 (tixagevimab/cilgavimab) in the prophylaxis and",
    "treatment of COVID-19. Antimicrob Agents Chemother. 2024;68(5):e01587-23.",
    "doi:10.1128/aac.01587-23. Parameter estimates from Table 1 and",
    "supplementary Table S3; model structure from supplementary Table S2 and",
    "Materials and Methods 'Population PK model development'.",
    sep = " "
  )
  vignette <- "Clegg_2024_tixagevimab_cilgavimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "AZD7442 (tixagevimab + cilgavimab)", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "AZD7442 (tixagevimab + cilgavimab)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "AZD7442 (tixagevimab + cilgavimab)", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight, normalised to a 70 kg reference. Supplementary Table S3 records every body-weight row as 'Baseline weight in kg on <parameter> (centered around: 70 kg)'; the exponents were fixed at 0.75 on the clearance-like parameters (CL, Q) and 1 on the volume-like parameters (Vc, Vp) rather than estimated (supplementary Table S2 'Covariates included in the base model'). Note that the 70 kg allometric reference is NOT the cohort median: the Figure 4 comparator subject weighs 80.6 kg and the cohort mean is 83.9 kg (supplementary Table S1).",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female). NOTE: the source model's reference sex is FEMALE, not the canonical column's 0 group -- see notes.",
      notes              = "Clegg 2024 codes sex as SEXM (1 = male) and reports beta_ka(SEXM_1) = 0.543 as a log-additive term on ka, so the typical ka = 0.11 /day belongs to a FEMALE subject. This file carries the canonical SEXF orientation (1 = female, 0 = male) and forms the paper's male indicator inside model() as (1 - SEXF), the same construction used by Bajaj_2017_nivolumab.R and Wada_2023_sparsentan.R. exp(0.543) - 1 = 72.1% faster absorption in males, matching the Results text.",
      source_name        = "SEXM"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model only through the paper's dichotomised AGECAT indicator (AGECAT_1 = age > 65 years, reference age <= 65 years), which model() forms as (AGE > 65). The continuous canonical column is carried rather than a pre-dichotomised one so that a user simulating the model supplies the demographic value itself and the cut point stays visible in model(). Cohort mean 50.5 years (SD 15.8, range 18.0-98.0); 871 of 4,940 participants (17.6%) were > 65 years (supplementary Table S1).",
      source_name        = "AGECAT"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model only through the paper's dichotomised BMICAT indicator (BMICAT_1 = BMI >= 30 kg/m2, reference BMI < 30 kg/m2), which model() forms as (BMI >= 30). Cohort mean 29.2 kg/m2 (SD 6.7, range 13.6-72.6); 1,959 of 4,940 participants (39.7%) had BMI >= 30 kg/m2 (supplementary Table S1). BMI is carried alongside WT because the paper fits both: WT drives the allometric scaling of the disposition parameters and BMI drives a separate absorption-rate effect.",
      source_name        = "BMICAT"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Time-fixed at study entry. Two log-additive effects: beta_ka(DIAB_1) = -0.318 on ka (exp(-0.318) - 1 = 27.2% slower absorption) and beta_CL(DIAB_1) = 0.198 on CL (exp(0.198) - 1 = 21.9% higher clearance), both matching the Results text. 638 of 4,940 participants (12.9%) had diabetes (supplementary Table S1).",
      source_name        = "DIAB"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black race: White, Asian, Other, not reported or unknown)",
      notes              = "Log-additive effect on the central volume of distribution: beta_Vc(RACEB_1) = -0.251, i.e. exp(-0.251) - 1 = 22.2% lower Vc in Black participants, matching the Results text. 668 of 4,940 participants (13.5%) were Black (supplementary Table S1). Kept as the plain RACE_BLACK dichotomy rather than the composite RACE_BLACK_OTH used in Clegg_2024_nirsevimab.R -- this paper's covariate is Black versus all other race groups, with no pooling.",
      source_name        = "RACEB"
    ),
    INJSITE_THIGH = list(
      description        = "Intramuscular injection into the anterolateral thigh",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (gluteal region)",
      notes              = "Per-dose-record indicator. Clegg 2024 carries this covariate under the study name ACTIV2 because the thigh was used only in ACTIV-2 while every other study injected into the gluteal region; supplementary Table S2 states explicitly that 'the covariate ACTIV-2 in this context is used to distinguish the site of IM administrations'. It is therefore registered here as the anatomical-site covariate rather than as a study indicator. Log-additive effect on bioavailability: beta_FIM(ACTIV2_1) = 0.416, i.e. exp(0.416) - 1 = 51.6% higher bioavailability for thigh injection, matching the Results text. Because the effect is log-additive and the gluteal FIM is already 0.671, the implied thigh bioavailability is 0.671 * exp(0.416) = 1.017, marginally above unity -- see the vignette Errata. Applies only to intramuscular doses (records routed through the depot); it has no effect on intravenous doses, which are placed directly into the central compartment.",
      source_name        = "ACTIV2"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous administration of AZD7442",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intramuscular). NOTE: the non-IV reference in this paper is IM, not the SC reference given in the register's general ROUTE_IV description.",
      notes              = "Per-dose-record dosing-route indicator (1 = IV infusion, 0 = IM injection) used to switch the log-scale additive residual SD between the IV value (expSdIv = 0.104) and the IM value (expSdIm = 0.24) per Table 1 / supplementary Table S3 rows error_ADD2 and error_ADD1. Distinct from the rxode2 cmt event column, which routes the dose (cmt = central for IV, cmt = depot for IM); ROUTE_IV is the covariate that picks the residual-error term. 342 of 4,940 participants (6.9%) received IV infusions and 4,598 (93.1%) received IM injections (supplementary Table S1). Same role as in Zierhut_2008_osteoprotegerin.R and Wang_2021_pertuzumab.R.",
      source_name        = "route of administration"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4940L,
    n_observations = 31895L,
    n_studies      = 8L,
    age_range      = "18.0-98.0 years",
    age_median     = "mean 50.5 years (SD 15.8); 17.6% were > 65 years",
    weight_range   = "36.0-216.0 kg",
    weight_median  = "mean 83.9 kg (SD 21.3); the Figure 4 typical-subject weight is the cohort median of 80.6 kg",
    bmi_range      = "13.6-72.6 kg/m^2 (mean 29.2, SD 6.7); 39.7% had BMI >= 30 kg/m^2",
    height_range   = "122.0-204.0 cm (mean 169, SD 10.1)",
    sex_female_pct = 46.9,
    race_ethnicity = c(White = 70.7, Black = 13.5, Asian = 8.7, `Other / not reported / unknown` = 7.0),
    ethnicity      = c(`Not Hispanic or Latino` = 71.4, `Hispanic or Latino` = 24.8, `Not reported / unknown` = 3.8),
    disease_state  = "Mixed prophylaxis and treatment populations: healthy adults (three phase I studies and one Chinese phase II study), adults at increased risk of inadequate response to COVID-19 vaccination or of SARS-CoV-2 exposure (PROVENT), adults exposed to SARS-CoV-2 within 8 days (STORM CHASER), and outpatients with mild-to-moderate COVID-19 (TACKLE, ACTIV-2). Comorbidities: cardiovascular disease 33.3%, diabetes 12.9%, chronic obstructive pulmonary disease 5.1%, chronic kidney disease 4.5%, chronic liver disease 4.1%, immunocompromised 0.8%.",
    dose_range     = "Single doses of 300 mg or 600 mg intramuscularly, and 300, 600, 1,000 or 3,000 mg intravenously. Some PROVENT participants received a second 300 mg IM dose 10-14 months after the first as part of a sub-study.",
    route          = "IM 4,598 participants (93.1%; gluteal region except ACTIV-2, which used the anterolateral thigh) and IV 342 participants (6.9%).",
    regions        = "North America, South America, Europe and Asia (including dedicated phase I studies in healthy Japanese and healthy Chinese adults and a phase II study in Chinese adults).",
    studies        = "Global phase I first-in-human NCT04507256 (n=50); PROVENT phase III pre-exposure prophylaxis NCT04625725 (n=3,288); STORM CHASER phase III post-exposure prophylaxis NCT04625972 (n=722); TACKLE phase III outpatient treatment NCT04723394 (n=442); ACTIV-2 phase II outpatient treatment NCT04518410 (n=157); Japanese phase I NCT04896541 (n=30); Chinese phase I NCT05437289 (n=49); Chinese phase II NCT05184062 (n=202).",
    reference_subject = "The typical values in Table 1 belong to a 70 kg FEMALE with age <= 65 years, BMI < 30 kg/m^2, non-Black race, no diabetes and a gluteal intramuscular injection. This is not the same as the Figure 4 covariate-comparison comparator, which is an 80.6 kg MALE with otherwise the same covariate levels.",
    notes          = "Baseline demographics from Clegg 2024 supplementary Table S1; study descriptions and sampling schemes from supplementary Table S6. Serum tixagevimab and cilgavimab were assayed by a validated hybrid ligand-binding LC-MS/MS method (LLOQ 0.3 ug/mL) and AZD7442 concentrations were computed as the sum of the two. Data below the limit of quantitation were handled with the NONMEM M3 method. Excluded from the analysis dataset: two ACTIV-2 participants injected at the wrong site, 170 participants with no observation records, 45 participants whose observations were all BLQ or missing, 13 implausibly high (> 100 ug/mL) STORM CHASER observations, and physiologically implausible intermittent BLQ observations from 79 participants."
  )

  ini({
    # Typical structural parameters at the reference covariate set defined in
    # population$reference_subject (70 kg, female, age <= 65 y, BMI < 30 kg/m2,
    # non-Black, no diabetes, gluteal IM injection). Time in days, dose in mg,
    # volumes in L, so Cc = central / vc is mg/L = ug/mL, matching the
    # published serum concentration units.
    lka     <- log(0.11);   label("First-order intramuscular absorption rate ka (1/day)")   # Clegg 2024 Table 1 / Table S3
    lcl     <- log(0.0504); label("Clearance CL at the reference covariates (L/day)")       # Clegg 2024 Table 1 / Table S3
    lvc     <- log(3.36);   label("Central volume of distribution Vc at the reference covariates (L)") # Clegg 2024 Table 1 / Table S3
    lq      <- log(0.395);  label("Inter-compartmental clearance Q at the reference covariates (L/day)") # Clegg 2024 Table 1 / Table S3
    lvp     <- log(1.83);   label("Peripheral volume of distribution Vp at the reference covariates (L)") # Clegg 2024 Table 1 / Table S3
    lfdepot <- log(0.671);  label("Absolute intramuscular bioavailability for a gluteal injection (fraction)") # Clegg 2024 Table 1 / Table S3

    # Body-weight allometry, exponents held at the standard 0.75 / 1 values
    # (supplementary Table S2: 'with fixed allometric exponents (0.75 for rates
    # and 1.0 for volumes)'); reported without RSE or CI in Table 1 and flagged
    # (FIX) in supplementary Table S3. Reference weight 70 kg.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")   # Clegg 2024 Table S3, beta_CL(BWT)
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on Vc (unitless)")   # Clegg 2024 Table S3, beta_Vc(BWT)
    e_wt_q  <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)")    # Clegg 2024 Table S3, beta_Q(BWT)
    e_wt_vp <- fixed(1);    label("Allometric exponent of body weight on Vp (unitless)")   # Clegg 2024 Table S3, beta_Vp(BWT)

    # Categorical covariate effects. Clegg 2024 Materials and Methods
    # ('Population PK model development'): "As log-transformed serum
    # concentrations were modeled, covariate-parameters relationships were
    # implemented as additive terms in the log-transformed domain of the
    # respective parameter." Every coefficient below is therefore a log-scale
    # additive shift, and exp(coefficient) - 1 reproduces the percentage
    # changes quoted in the Results.
    e_sexf_ka          <-  0.543; label("Log-scale effect of male sex on ka, applied to (1 - SEXF) (unitless)") # Clegg 2024 Table 1, beta_ka(SEXM_1)
    e_age_ka           <- -0.325; label("Log-scale effect of age > 65 years on ka (unitless)")                  # Clegg 2024 Table 1, beta_ka(AGECAT_1)
    e_bmi_ka           <- -0.254; label("Log-scale effect of BMI >= 30 kg/m2 on ka (unitless)")                 # Clegg 2024 Table 1, beta_ka(BMICAT_1)
    e_dis_diab_ka      <- -0.318; label("Log-scale effect of diabetes on ka (unitless)")                        # Clegg 2024 Table 1, beta_ka(DIAB_1)
    e_dis_diab_cl      <-  0.198; label("Log-scale effect of diabetes on CL (unitless)")                        # Clegg 2024 Table 1, beta_CL(DIAB_1)
    e_race_black_vc    <- -0.251; label("Log-scale effect of Black race on Vc (unitless)")                      # Clegg 2024 Table 1, beta_Vc(RACEB_1)
    e_injsite_thigh_fdepot <- 0.416; label("Log-scale effect of thigh (versus gluteal) intramuscular injection on bioavailability (unitless)") # Clegg 2024 Table 1, beta_FIM(ACTIV2_1)

    # Inter-individual variability. Table 1 and supplementary Table S3 report
    # log-normal IIV as %CV, converted here with omega^2 = log(1 + CV^2):
    #   ka  59.89% -> 0.3065145      cl  42.64% -> 0.1670531
    #   vc  39.79% -> 0.1469745      vp  35.88% -> 0.1210997
    # Supplementary Table S2 states that "Correlation between the random
    # effects of ka, CL, and Vc was considered", so those three form a 3x3
    # block and Vp carries an independent eta. Off-diagonals are
    # corr * sqrt(var_i * var_j) using Table 1 corr(ka,CL) = -0.387,
    # corr(ka,Vc) = -0.689 and corr(CL,Vc) = 0.588; the resulting Omega is
    # positive definite (smallest eigenvalue 0.0451).
    # FIM and Q IIV are reported as 0% (FIX) in supplementary Table S3, i.e.
    # no random effect was estimated on either, so no etalfdepot / etalq
    # appears here -- a fixed(0) variance would make Omega singular.
    etalka + etalcl + etalvc ~ c(
       0.3065145,
      -0.0875717,  0.1670531,
      -0.1462398,  0.0921352,  0.1469745
    )                                            # Clegg 2024 Table 1 IIV %CV and correlation-of-random-effects rows
    etalvp ~ 0.1210997                           # Clegg 2024 Table 1, Vp 35.9 %CV (Table S3: 35.88%)

    # Residual error. Clegg 2024 fit log-transformed serum concentrations with
    # "separate additive residual error models for log transformed IV and IM
    # data" (supplementary Table S2); supplementary Table S3 gives the units of
    # both terms as log(ug/mL). An additive error on the natural-log scale is
    # exactly a log-normal residual on the linear scale, so these map onto
    # nlmixr2's `~ lnorm(...)` with the SD selected per record by ROUTE_IV.
    expSdIm <- 0.24;  label("Log-scale additive residual SD for intramuscular records (log ug/mL)") # Clegg 2024 Table 1, error_ADD1
    expSdIv <- 0.104; label("Log-scale additive residual SD for intravenous records (log ug/mL)")   # Clegg 2024 Table 1, error_ADD2
  })

  model({
    # 1. Derived covariate terms. Clegg 2024 supplies sex as a male indicator
    #    and age / BMI as pre-dichotomised categories; the canonical columns
    #    are SEXF and the continuous AGE / BMI, so the paper's indicators are
    #    reconstructed here with the published cut points.
    sexm   <- 1 - SEXF
    agecat <- AGE > 65            # AGECAT_1 = age > 65 years
    bmicat <- BMI >= 30           # BMICAT_1 = BMI >= 30 kg/m2

    # 2. Individual parameters. All covariate effects are additive on the
    #    log-transformed parameter (Materials and Methods), and body weight
    #    enters as a power function normalised to 70 kg.
    ka <- exp(lka + etalka +
                e_sexf_ka * sexm +
                e_age_ka * agecat +
                e_bmi_ka * bmicat +
                e_dis_diab_ka * DIS_DIAB)
    cl <- exp(lcl + etalcl + e_dis_diab_cl * DIS_DIAB) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc + e_race_black_vc * RACE_BLACK) * (WT / 70)^e_wt_vc
    q  <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # Absolute intramuscular bioavailability; the thigh effect is log-additive
    # like every other covariate in this model.
    fdepot <- exp(lfdepot + e_injsite_thigh_fdepot * INJSITE_THIGH)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. IM doses enter the depot and are absorbed first-order;
    #    IV doses are placed directly into the central compartment (zero-order
    #    input handled by the event table's rate / dur, per supplementary
    #    Table S2 'Zero-order administration of IV doses into the central
    #    compartment').
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability applies to the depot (IM) route only.
    f(depot) <- fdepot

    # 6. Observation and route-specific log-scale residual error.
    Cc <- central / vc
    expSd <- expSdIm + (expSdIv - expSdIm) * ROUTE_IV
    Cc ~ lnorm(expSd)
  })
}
