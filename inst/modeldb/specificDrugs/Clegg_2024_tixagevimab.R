Clegg_2024_tixagevimab <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order intramuscular",
    "absorption and linear elimination for tixagevimab, the YTE-modified",
    "SARS-CoV-2-neutralising monoclonal antibody component of AZD7442",
    "(Evusheld), in 4,940 adults from eight COVID-19 prophylaxis and treatment",
    "trials (Clegg 2024 Table 1). This is the tixagevimab-analyte refit of the",
    "final AZD7442 model: same structure and same covariate model, parameters",
    "re-estimated against the tixagevimab serum concentrations. Body weight is",
    "carried on CL, Vc, Q and Vp with allometric exponents fixed at 0.75 and 1",
    "and a 70 kg reference. Absorption rate carries male sex, age > 65 years,",
    "BMI >= 30 kg/m2 and diabetes; clearance carries diabetes; central volume",
    "carries Black race; and intramuscular bioavailability carries the",
    "injection site (thigh versus gluteal). IIV is a 3x3 block on ka, CL and",
    "Vc plus an independent eta on Vp, and the residual error is additive on",
    "the log scale with separate magnitudes for intramuscular and intravenous",
    "records. Companion models:",
    "modellib('Clegg_2024_tixagevimab_cilgavimab') for the summed AZD7442",
    "analyte and modellib('Clegg_2024_cilgavimab') for the other component.",
    sep = " "
  )
  reference <- paste(
    "Clegg LE, Stepanov O, Schmidt H, Tang W, Zhang H, Webber C, Cohen TS,",
    "Esser MT, Nagard M. Accelerating therapeutics development during a",
    "pandemic: population pharmacokinetics of the long-acting antibody",
    "combination AZD7442 (tixagevimab/cilgavimab) in the prophylaxis and",
    "treatment of COVID-19. Antimicrob Agents Chemother. 2024;68(5):e01587-23.",
    "doi:10.1128/aac.01587-23. Parameter estimates from the Tixagevimab column",
    "of Table 1, cross-checked against supplementary Table S5 (whose caption is",
    "swapped with Table S4 -- see the vignette Errata); model structure from",
    "supplementary Table S2 and Materials and Methods 'Population PK model",
    "development'.",
    sep = " "
  )
  vignette <- "Clegg_2024_tixagevimab_cilgavimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "tixagevimab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tixagevimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tixagevimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight, normalised to a 70 kg reference. Supplementary Table S5 records every body-weight row as 'Baseline weight in kg on <parameter> (centered around: 70 kg)'; the exponents were fixed at 0.75 on the clearance-like parameters (CL, Q) and 1 on the volume-like parameters (Vc, Vp) rather than estimated (supplementary Table S2). The 70 kg allometric reference is not the cohort median: the Figure 4 comparator subject weighs 80.6 kg and the cohort mean is 83.9 kg (supplementary Table S1).",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female). NOTE: the source model's reference sex is FEMALE, not the canonical column's 0 group -- see notes.",
      notes              = "Clegg 2024 codes sex as SEXM (1 = male) and reports beta_ka(SEXM_1) = 0.632 as a log-additive term on ka, so the typical ka = 0.122 /day belongs to a FEMALE subject. This file carries the canonical SEXF orientation (1 = female, 0 = male) and forms the paper's male indicator inside model() as (1 - SEXF), the same construction used by Bajaj_2017_nivolumab.R and Wada_2023_sparsentan.R.",
      source_name        = "SEXM"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model only through the paper's dichotomised AGECAT indicator (AGECAT_1 = age > 65 years, reference age <= 65 years), which model() forms as (AGE > 65). Cohort mean 50.5 years (SD 15.8, range 18.0-98.0); 871 of 4,940 participants (17.6%) were > 65 years (supplementary Table S1).",
      source_name        = "AGECAT"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model only through the paper's dichotomised BMICAT indicator (BMICAT_1 = BMI >= 30 kg/m2, reference BMI < 30 kg/m2), which model() forms as (BMI >= 30). Cohort mean 29.2 kg/m2 (SD 6.7, range 13.6-72.6); 39.7% had BMI >= 30 kg/m2 (supplementary Table S1). Carried alongside WT because the paper fits both: WT drives the allometric scaling of disposition and BMI drives a separate absorption-rate effect.",
      source_name        = "BMICAT"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Time-fixed at study entry. Two log-additive effects: beta_ka(DIAB_1) = -0.262 on ka (23.1% slower absorption) and beta_CL(DIAB_1) = 0.167 on CL (18.2% higher clearance). 638 of 4,940 participants (12.9%) had diabetes (supplementary Table S1).",
      source_name        = "DIAB"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black race: White, Asian, Other, not reported or unknown)",
      notes              = "Log-additive effect on the central volume of distribution: beta_Vc(RACEB_1) = -0.242, i.e. 21.5% lower Vc in Black participants. 668 of 4,940 participants (13.5%) were Black (supplementary Table S1). Kept as the plain RACE_BLACK dichotomy rather than the composite RACE_BLACK_OTH used in Clegg_2024_nirsevimab.R -- this paper's covariate is Black versus all other race groups, with no pooling.",
      source_name        = "RACEB"
    ),
    INJSITE_THIGH = list(
      description        = "Intramuscular injection into the anterolateral thigh",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (gluteal region)",
      notes              = "Per-dose-record indicator. Clegg 2024 carries this covariate under the study name ACTIV2 because the thigh was used only in ACTIV-2 while every other study injected into the gluteal region; supplementary Table S2 states explicitly that 'the covariate ACTIV-2 in this context is used to distinguish the site of IM administrations'. Log-additive effect on bioavailability: beta_FIM(ACTIV2_1) = 0.378, i.e. 45.9% higher bioavailability for thigh injection (0.615 * exp(0.378) = 0.897). Applies only to intramuscular doses; intravenous doses bypass the depot.",
      source_name        = "ACTIV2"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intramuscular). NOTE: the non-IV reference in this paper is IM, not the SC reference given in the register's general ROUTE_IV description.",
      notes              = "Per-dose-record dosing-route indicator (1 = IV infusion, 0 = IM injection) used to switch the log-scale additive residual SD between the IV value (expSdIv = 0.108) and the IM value (expSdIm = 0.272) per Table 1 rows error_ADD2 and error_ADD1. Distinct from the rxode2 cmt event column, which routes the dose (cmt = central for IV, cmt = depot for IM). 342 of 4,940 participants (6.9%) received IV infusions (supplementary Table S1). Same role as in Zierhut_2008_osteoprotegerin.R and Wang_2021_pertuzumab.R.",
      source_name        = "route of administration"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4940L,
    n_studies      = 8L,
    age_range      = "18.0-98.0 years",
    age_median     = "mean 50.5 years (SD 15.8); 17.6% were > 65 years",
    weight_range   = "36.0-216.0 kg",
    weight_median  = "mean 83.9 kg (SD 21.3); the Figure 4 typical-subject weight is the cohort median of 80.6 kg",
    bmi_range      = "13.6-72.6 kg/m^2 (mean 29.2, SD 6.7); 39.7% had BMI >= 30 kg/m^2",
    sex_female_pct = 46.9,
    race_ethnicity = c(White = 70.7, Black = 13.5, Asian = 8.7, `Other / not reported / unknown` = 7.0),
    disease_state  = "Mixed prophylaxis and treatment populations: healthy adults (three phase I studies and one Chinese phase II study), adults at increased risk of inadequate response to COVID-19 vaccination or of SARS-CoV-2 exposure (PROVENT), adults exposed to SARS-CoV-2 within 8 days (STORM CHASER), and outpatients with mild-to-moderate COVID-19 (TACKLE, ACTIV-2).",
    dose_range     = "Single doses of 300 mg or 600 mg AZD7442 intramuscularly, and 300, 600, 1,000 or 3,000 mg AZD7442 intravenously. AZD7442 is a 1:1 combination, so a 300 mg AZD7442 dose delivers 150 mg tixagevimab.",
    route          = "IM 4,598 participants (93.1%; gluteal region except ACTIV-2, which used the anterolateral thigh) and IV 342 participants (6.9%).",
    regions        = "North America, South America, Europe and Asia.",
    reference_subject = "The typical values in Table 1 belong to a 70 kg FEMALE with age <= 65 years, BMI < 30 kg/m^2, non-Black race, no diabetes and a gluteal intramuscular injection.",
    notes          = "Same 4,940-participant pooled dataset as the AZD7442 model; the analyte differs. Clegg 2024 Materials and Methods: 'abbreviated model development for tixagevimab and cilgavimab was performed by refining the final AZD7442 model on an as-needed basis'. Serum tixagevimab was assayed by a validated hybrid ligand-binding LC-MS/MS method with LLOQ 0.3 ug/mL; BLQ data were handled with the NONMEM M3 method. Baseline demographics from supplementary Table S1. IMPORTANT: dose amounts supplied to this model must be the tixagevimab amount, i.e. half the labelled AZD7442 dose."
  )

  ini({
    # Typical structural parameters at the reference covariate set (70 kg,
    # female, age <= 65 y, BMI < 30 kg/m2, non-Black, no diabetes, gluteal IM
    # injection). Time in days, dose in mg, volumes in L, so Cc = central / vc
    # is mg/L = ug/mL.
    lka     <- log(0.122);  label("First-order intramuscular absorption rate ka (1/day)")   # Clegg 2024 Table 1, Tixagevimab column
    lcl     <- log(0.0457); label("Clearance CL at the reference covariates (L/day)")       # Clegg 2024 Table 1, Tixagevimab column
    lvc     <- log(3.17);   label("Central volume of distribution Vc at the reference covariates (L)") # Clegg 2024 Table 1, Tixagevimab column
    lq      <- log(0.432);  label("Inter-compartmental clearance Q at the reference covariates (L/day)") # Clegg 2024 Table 1, Tixagevimab column
    lvp     <- log(1.77);   label("Peripheral volume of distribution Vp at the reference covariates (L)") # Clegg 2024 Table 1, Tixagevimab column
    lfdepot <- log(0.615);  label("Absolute intramuscular bioavailability for a gluteal injection (fraction)") # Clegg 2024 Table 1, Tixagevimab column

    # Body-weight allometry, exponents held at the standard 0.75 / 1 values
    # (supplementary Table S2); reported without RSE or CI in Table 1 and
    # flagged (FIX) in supplementary Table S5. Reference weight 70 kg.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)")   # Clegg 2024 Table 1, beta_CL(BWT)
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on Vc (unitless)")   # Clegg 2024 Table 1, beta_Vc(BWT)
    e_wt_q  <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)")    # Clegg 2024 Table 1, beta_Q(BWT)
    e_wt_vp <- fixed(1);    label("Allometric exponent of body weight on Vp (unitless)")   # Clegg 2024 Table 1, beta_Vp(BWT)

    # Categorical covariate effects, all additive on the log-transformed
    # parameter (Clegg 2024 Materials and Methods, 'Population PK model
    # development').
    e_sexf_ka          <-  0.632; label("Log-scale effect of male sex on ka, applied to (1 - SEXF) (unitless)") # Clegg 2024 Table 1, beta_ka(SEXM_1)
    e_age_ka           <- -0.443; label("Log-scale effect of age > 65 years on ka (unitless)")                  # Clegg 2024 Table 1, beta_ka(AGECAT_1)
    e_bmi_ka           <- -0.217; label("Log-scale effect of BMI >= 30 kg/m2 on ka (unitless)")                 # Clegg 2024 Table 1, beta_ka(BMICAT_1)
    e_dis_diab_ka      <- -0.262; label("Log-scale effect of diabetes on ka (unitless)")                        # Clegg 2024 Table 1, beta_ka(DIAB_1)
    e_dis_diab_cl      <-  0.167; label("Log-scale effect of diabetes on CL (unitless)")                        # Clegg 2024 Table 1, beta_CL(DIAB_1)
    e_race_black_vc    <- -0.242; label("Log-scale effect of Black race on Vc (unitless)")                      # Clegg 2024 Table 1, beta_Vc(RACEB_1)
    e_injsite_thigh_fdepot <- 0.378; label("Log-scale effect of thigh (versus gluteal) intramuscular injection on bioavailability (unitless)") # Clegg 2024 Table 1, beta_FIM(ACTIV2_1)

    # Inter-individual variability, log-normal, reported as %CV and converted
    # with omega^2 = log(1 + CV^2):
    #   ka  78.56% -> 0.4806761      cl  40.81% -> 0.1540469
    #   vc  38.00% -> 0.1348805      vp  36.85% -> 0.1273304
    # ka, CL and Vc form a 3x3 block (supplementary Table S2: "Correlation
    # between the random effects of ka, CL, and Vc was considered") with
    # off-diagonals corr * sqrt(var_i * var_j) from Table 1
    # corr(ka,CL) = -0.496, corr(ka,Vc) = -0.835, corr(CL,Vc) = 0.697; the
    # resulting Omega is positive definite (smallest eigenvalue 0.0210).
    # FIM and Q IIV are reported as 0% (FIX), i.e. no random effect was
    # estimated on either, so no etalfdepot / etalq appears here.
    etalka + etalcl + etalvc ~ c(
       0.4806761,
      -0.1349691,  0.1540469,
      -0.2126118,  0.1004694,  0.1348805
    )                                            # Clegg 2024 Table 1 IIV %CV and correlation-of-random-effects rows
    etalvp ~ 0.1273304                           # Clegg 2024 Table 1, Vp 36.8 %CV (Table S5: 36.85%)

    # Residual error: additive on the natural-log scale (supplementary
    # Table S5 gives the units as log(ug/mL)), which is exactly a log-normal
    # residual on the linear scale, so `~ lnorm(...)` with the SD selected per
    # record by ROUTE_IV.
    expSdIm <- 0.272; label("Log-scale additive residual SD for intramuscular records (log ug/mL)") # Clegg 2024 Table 1, error_ADD1
    expSdIv <- 0.108; label("Log-scale additive residual SD for intravenous records (log ug/mL)")   # Clegg 2024 Table 1, error_ADD2
  })

  model({
    # 1. Derived covariate terms.
    sexm   <- 1 - SEXF
    agecat <- AGE > 65            # AGECAT_1 = age > 65 years
    bmicat <- BMI >= 30           # BMICAT_1 = BMI >= 30 kg/m2

    # 2. Individual parameters.
    ka <- exp(lka + etalka +
                e_sexf_ka * sexm +
                e_age_ka * agecat +
                e_bmi_ka * bmicat +
                e_dis_diab_ka * DIS_DIAB)
    cl <- exp(lcl + etalcl + e_dis_diab_cl * DIS_DIAB) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc + e_race_black_vc * RACE_BLACK) * (WT / 70)^e_wt_vc
    q  <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    fdepot <- exp(lfdepot + e_injsite_thigh_fdepot * INJSITE_THIGH)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. IM doses enter the depot; IV doses are placed directly
    #    into the central compartment (zero-order input handled by the event
    #    table's rate / dur).
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
