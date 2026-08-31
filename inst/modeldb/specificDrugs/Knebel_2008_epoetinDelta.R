Knebel_2008_epoetinDelta <- function() {
  description <- "One-compartment population PK model with first-order absorption and linear elimination for subcutaneous and intravenous epoetin delta in pediatric patients aged 1-17 years with chronic kidney disease (Knebel 2008). Serum erythropoietin is the sum of the drug contribution and an additive endogenous baseline concentration. Clearance and central volume are allometrically scaled by body weight (exponents 0.75 and 1) about a 35 kg reference, with additional power effects of age above 10 years, a sex factor on both, and dialysis-modality factors on volume. The same model was fit jointly to epoetin alfa data, which enters through a product indicator that shifts the absorption rate constant and the subcutaneous bioavailability."
  reference <- paste(
    "Knebel W, Palmen M, Dowell JA, Gastonguay M.",
    "Population pharmacokinetic modeling of epoetin delta in pediatric patients",
    "with chronic kidney disease.",
    "J Clin Pharmacol. 2008;48(7):837-848.",
    "doi:10.1177/0091270008318218",
    sep = " "
  )
  vignette <- "Knebel_2008_epoetinDelta"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "mIU/mL"
  )

  # Serum erythropoietin was measured by ELISA in mEU/mL and converted to
  # bioassay/potency international units at 1 IU = 1.3 EU (Knebel 2008
  # "Biological Methods"), so the reported concentrations and the doses share
  # the same activity unit. With doses in IU and volumes in L the state/volume
  # ratio is IU/L, which equals mIU/mL numerically.
  compartmentData <- list(
    depot   = list(analyte = "epoetin delta", units = "IU", specimen = "administration site", verified = TRUE),
    central = list(analyte = "epoetin delta", units = "IU", specimen = "serum",               verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Knebel 2008 Table II: median 34.5 kg, mean 35.6 kg, range 11.3-83 kg. Enters CL and Vc as the allometric power model of equation (3), TVP = theta * (WT/WTref)^theta_allo, with the exponents held at the physiologic values 0.75 (clearance) and 1 (volume). The Results section states the weight normalization for the full covariate model was WTref = 35 kg; the base model used a 70 kg reference (Table III footnote a), so base-model and final-model typical values are not directly comparable. Weight was chosen over BMI and BSA because the three body-size metrics correlated at r >= 0.71.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Knebel 2008 Table II: median 13 years, mean 11.8 years, range 1-17 years. Results: 'Age was included as a power function, normalized by the reference age of 10 years for all patients older than 10 years of age.' The effect is therefore a hinge - the (AGE/10)^theta term applies only above 10 years and equals 1 at or below it, which is why the abstract's reference individual is a '35-kg male <= 10 years'. Implemented as max(AGE, 10)/10. Acts on both CL (exponent 0.999) and Vc (exponent 2.89).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Knebel 2008 Table II: 38 of 60 subjects (63%) were male, so 22 (37%) were female. Results: 'sex entered the model as a power function, with a separate dichotomous (0,1) covariate serving as an on-off switch', i.e. the Table III rows *theta9^SEX[female] and *theta10^SEX[female]. Male is the reference (the abstract's reference individual is male). The paper's source column labels the female level; no value transformation is needed.",
      source_name        = "SEX"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis modality indicator (1 = hemodialysis, 0 = predialysis or peritoneal dialysis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (predialysis; see notes)",
      notes              = "Knebel 2008 Table II: 28 of 60 subjects (47%) were on hemodialysis, 15 (25%) on peritoneal dialysis, 17 (28%) predialysis. Results: 'dialysis type was entered multiplicatively, with the predialysis type as a reference and a switch to determine what dialysis type effect was being estimated', i.e. the Table III row *theta11^DIT1[hemodialysis]. Paired with PERIT_DIAL as the two switches of a three-level covariate: predialysis subjects carry RRT_HEMODIAL_STATUS = 0 and PERIT_DIAL = 0, so the reference category here is predialysis rather than the 'not on hemodialysis' pooling implied by the bare canonical definition. Acts on Vc only; a dialysis-type effect on CL was not retained.",
      source_name        = "DIT1"
    ),
    PERIT_DIAL = list(
      description        = "Peritoneal-dialysis modality indicator (1 = peritoneal dialysis, 0 = predialysis or hemodialysis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (predialysis; see notes)",
      notes              = "Knebel 2008 Table II: 15 of 60 subjects (25%) were on peritoneal dialysis. Table III row *theta12^DIT2[peritoneal dialysis]. The second of the two dialysis-modality switches described in Results; see the RRT_HEMODIAL_STATUS notes for the three-level encoding. Unlike Takama 2007 - which took hemodialysis as the reference and reports PERIT_DIAL relative to it - Knebel 2008 takes predialysis as the reference, so both indicators are 0 for the predialysis stratum and the two effects are read relative to predialysis, not relative to each other. Acts on Vc only.",
      source_name        = "DIT2"
    ),
    FORM_EPO_ALFA = list(
      description        = "Epoetin alfa product indicator (1 = epoetin alfa, CHO-cell-derived; 0 = epoetin delta, human-cell-line-derived)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (epoetin delta)",
      notes              = "Knebel 2008 Table I: 47 of 60 subjects received epoetin delta and 13 received epoetin alfa; patients were randomized approximately 3:1. Results: 'The base model was then modified to allow Ka and F1 to differ between epoetin delta and epoetin alfa SC administration', i.e. the Table III rows *theta6^TRT[epoetin alfa] on Ka and *theta7^TRT[epoetin alfa] on F1. Subject-level and time-fixed (each subject received one product for the whole 24-week study). Both products are recombinant human erythropoietin and are quantified by the same ELISA, so the covariate distinguishes the drug product rather than the analyte. The two effects are estimated only from the subcutaneous arms; the single intravenous epoetin alfa subject contributes no absorption information, and neither term reaches the central compartment for an intravenous dose.",
      source_name        = "TRT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 60L,
    n_studies        = 1L,
    n_concentrations = 261L,
    age_range        = "1-17 years",
    age_median       = "13 years",
    weight_range     = "11.3-83 kg",
    weight_median    = "34.5 kg",
    sex_female_pct   = 37,
    race_ethnicity   = c(White = 90, `African American` = 3, Multiracial = 7),
    disease_state    = "Pediatric chronic kidney disease with associated anemia. Knebel 2008 Table II: 28 subjects (47%) on hemodialysis, 15 (25%) on peritoneal dialysis, 17 (28%) predialysis. Primary diagnoses leading to CKD included glomerulonephritis (9, 15%), secondary glomerulonephritis/vasculitis (1, 2%), and interstitial nephritis/pyelonephritis (4, 7%). All patients were on a stable erythropoiesis-stimulating agent before entry with hemoglobin 10-13 g/dL, which is why a subject-specific baseline erythropoietin concentration was not estimable.",
    dose_range       = "Epoetin delta: 26-191 IU/kg subcutaneous (37 subjects), 54-769 IU/kg intravenous (10 subjects). Epoetin alfa: 24-190 IU/kg subcutaneous (12 subjects), 36-88 IU/kg intravenous (1 subject). Administered once weekly, twice weekly, or three times weekly for 24 weeks, titrated to keep hemoglobin in the 10-13 g/dL target range (Knebel 2008 Table I and Study Design).",
    regions          = "United States and Argentina (per the ethics-committee appendix of Knebel 2008)",
    body_size        = "Knebel 2008 Table II also reports body mass index (median 17.2, mean 18.1, range 12.3-28.6 kg/m^2) and body surface area (median 1.17, mean 1.15, range 0.506-1.95 m^2). Neither was retained; weight was chosen as the single body-size metric because weight, BSA, BMI, and age were correlated at r >= 0.71 except BSA-BMI (r = 0.59) and age-BMI (r = 0.33).",
    notes            = "Phase III, open-label, randomized, stratified, multicenter 24-week study of epoetin delta (Dynepo). Sparse sampling: most patients gave one PK sample per scheduled visit (week 2 predose, week 4 at 2-4 h, week 8 at 12-24 h, week 12 at 20-36 h, week 20 at 48 h postdose), averaging 4 samples per patient, with most samples within 50 h of a dose. Serum erythropoietin was measured by a validated ELISA linear over 2.50-250 mEU/mL with an LLQ of 2.50 mEU/mL. Fit in NONMEM VI level 1.1 with FOCE-INT. Race was recorded but no formal comparisons across racial categories were made because of the small African American and multiracial strata, so race is not a covariate in the model."
  )

  ini({
    # Structural parameters. All values are the Knebel 2008 Table III "Final
    # Model" column (NONMEM final estimates, not initial values); the
    # parenthesised numbers in the trailing comments are the %SE and the
    # nonparametric-bootstrap 95% CI reported in the same table. The reference
    # individual is a 35 kg male aged 10 years or less, predialysis, receiving
    # subcutaneous epoetin delta (Knebel 2008 abstract and Results). The paper
    # reports CL, V, Ka, F1, and Bepo on the linear scale; they are stored on
    # the log scale here per the nlmixr2lib convention.
    lcl     <- log(0.268);  label("Clearance at the reference covariate values (L/h)")                          # Knebel 2008 Table III final model: theta1 = 0.268 L/h (35% SE; bootstrap 95% CI 0.148, 0.827)
    lvc     <- log(1.03);   label("Central volume of distribution at the reference covariate values (L)")       # Knebel 2008 Table III final model: theta2 = 1.03 L (45% SE; bootstrap 95% CI 0.344, 6.44)
    lka     <- log(0.0554); label("First-order absorption rate constant for subcutaneous epoetin delta (1/h)")  # Knebel 2008 Table III final model: theta3 = 0.0554 1/h (16% SE; bootstrap 95% CI 0.0405, 0.199)
    lfdepot <- log(0.708);  label("Subcutaneous bioavailability of epoetin delta (unitless)")                   # Knebel 2008 Table III final model: theta4 = 0.708 (38% SE; bootstrap 95% CI 0.337, 2.12)
    lrbase  <- log(6.71);   label("Baseline endogenous erythropoietin serum concentration Bepo (mIU/mL)")       # Knebel 2008 Table III final model: theta5 = 6.71 mIU/mL (7% SE; bootstrap 95% CI 5.72, 7.70)

    # Allometric exponents. Knebel 2008 Results: "V and CL were allometrically
    # scaled by weight using a power model with the exponents fixed to 0.75 for
    # CL and 1 for V, as described in equation (3)". Table III shows these rows
    # as NA in every estimate column, confirming they were not estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/35) for CL (unitless)")  # Knebel 2008 Results and equation (3); Table III row *(WT/35)^0.75 reported NA
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/35) for Vc (unitless)")  # Knebel 2008 Results and equation (3); Table III row *(WT/35)^1 reported NA

    # Covariate effects, all estimated. Knebel 2008 equation (4) gives the
    # normalized-power form for continuous covariates and the
    # exponentiated-switch form theta^cov for categorical (0-1) covariates, so
    # the age coefficients are exponents on (AGE/10) while the sex, dialysis,
    # and product coefficients are multiplicative factors raised to the 0/1
    # indicator. None of these are held fixed by the paper.
    e_age_cl         <- 0.999; label("Power exponent on (AGE/10) for CL above the 10-year reference age (unitless)")  # Knebel 2008 Table III final model: theta8 = 0.999 (54% SE; bootstrap 95% CI -0.0690, 2.1)
    e_sexf_cl        <- 0.923; label("Multiplicative factor on CL when SEXF = 1 (unitless)")                          # Knebel 2008 Table III final model: theta9 = 0.923 (18% SE; bootstrap 95% CI 0.594, 1.32)
    e_age_vc         <- 2.89;  label("Power exponent on (AGE/10) for Vc above the 10-year reference age (unitless)")  # Knebel 2008 Table III final model: theta13 = 2.89 (29% SE; bootstrap 95% CI 0.0281, 4.95)
    e_sexf_vc        <- 0.994; label("Multiplicative factor on Vc when SEXF = 1 (unitless)")                          # Knebel 2008 Table III final model: theta10 = 0.994 (39% SE; bootstrap 95% CI 0.370, 2.31)
    e_hemodial_vc    <- 4.53;  label("Multiplicative factor on Vc when RRT_HEMODIAL_STATUS = 1 (unitless)")           # Knebel 2008 Table III final model: theta11 = 4.53 (45% SE; bootstrap 95% CI 1.65, 16.5)
    e_perit_dial_vc  <- 2.48;  label("Multiplicative factor on Vc when PERIT_DIAL = 1 (unitless)")                    # Knebel 2008 Table III final model: theta12 = 2.48 (39% SE; bootstrap 95% CI 1.08, 6.13)
    e_epoalfa_ka     <- 1.23;  label("Multiplicative factor on Ka when FORM_EPO_ALFA = 1 (unitless)")                 # Knebel 2008 Table III final model: theta6 = 1.23 (14% SE; bootstrap 95% CI 0.495, 1.62); Results quote it as relative Ka = 123% [50%, 162%]
    e_epoalfa_fdepot <- 0.544; label("Multiplicative factor on F1 when FORM_EPO_ALFA = 1 (unitless)")                 # Knebel 2008 Table III final model: theta7 = 0.544 (19% SE; bootstrap 95% CI 0.368, 0.824); Results quote it as relative bioavailability F1 = 54.4% [36.8%, 82.4%]

    # Inter-individual variability. Knebel 2008 equation (1) is the exponential
    # / lognormal form Pi = P * exp(eta_Pi) with eta ~ N(0, omega^2), so the
    # Table III "Interindividual Variance" entries are omega^2 on the log scale
    # and go into ini() unchanged. The paper's %CV column is the square root of
    # the variance rather than the lognormal sqrt(exp(omega^2) - 1): sqrt(0.387)
    # = 0.622 and sqrt(1.41) = 1.19 reproduce the printed 62.2% and 119%
    # exactly, which is what fixes these as variances. Table III reports only
    # the diagonal elements Omega(1,1) and Omega(2,2); although Methods says a
    # full block was attempted, no covariance is published, so the etas are
    # encoded as independent. Bepo carries no IIV - Results: "All patients were
    # stabilized on an ESA prior to study entry, resulting in ... the inability
    # to estimate a subject-specific baseline erythropoietin concentration."
    etalcl ~ 0.387  # Knebel 2008 Table III final model: Omega(1,1) on CL = 0.387 (29% SE; CV% = 62.2; bootstrap 95% CI 0.184, 0.671)
    etalvc ~ 1.41   # Knebel 2008 Table III final model: Omega(2,2) on V  = 1.41  (33% SE; CV% = 119;  bootstrap 95% CI 0.481, 2.51)

    # Residual error. Knebel 2008 equation (2) is Cij = Chat_ij * exp(eps_pij),
    # an exponential (log-normal) residual on the model prediction inclusive of
    # the endogenous baseline, which maps to nlmixr2's ~ lnorm(expSd). Table III
    # reports the residual VARIANCE sigma^2 = 0.245; expSd is the standard
    # deviation, sqrt(0.245) = 0.49497, and the paper's own CV% = 49.5 column is
    # that same square root, confirming the variance reading.
    expSd <- 0.49497; label("Exponential (log-scale) residual standard deviation (unitless)")  # Knebel 2008 Table III final model: sigma(1,1) = 0.245 (13% SE; CV% = 49.5; bootstrap 95% CI 0.184, 0.308); sqrt(0.245) = 0.49497
  })

  model({
    # Reference covariate values for the normalized power terms. Knebel 2008
    # Results: "The weight normalization for the allometric relationships for
    # the effect of weight on CL and V was 35 kg for the full covariate model"
    # and age was "normalized by the reference age of 10 years".
    wt_ref  <- 35  # kg
    age_ref <- 10  # years

    # Age hinge. Knebel 2008 Results applies the (AGE/10) power term only "for
    # all patients older than 10 years of age", so the term is identically 1 at
    # or below the reference age. That is what makes the abstract's reference
    # individual a "35-kg male <= 10 years" rather than one exactly 10 years old.
    age_term <- max(AGE, age_ref) / age_ref

    # Individual parameters. Continuous covariates use the normalized power
    # form of equation (4) and categorical (0-1) covariates the exponentiated
    # switch form, exactly as laid out row by row in Table III:
    #   CL = theta1 * (WT/35)^0.75 * (AGE/10)^theta8  * theta9^SEX
    #   V  = theta2 * (WT/35)^1    * theta10^SEX * theta11^DIT1 * theta12^DIT2 * (AGE/10)^theta13
    #   Ka = theta3 * theta6^TRT
    #   F1 = theta4 * theta7^TRT
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl * age_term^e_age_cl * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc * age_term^e_age_vc * e_sexf_vc^SEXF *
      e_hemodial_vc^RRT_HEMODIAL_STATUS * e_perit_dial_vc^PERIT_DIAL

    # Absorption rate and bioavailability carry no IIV in Table III and shift
    # only with the product indicator.
    ka     <- exp(lka)     * e_epoalfa_ka^FORM_EPO_ALFA
    fdepot <- exp(lfdepot) * e_epoalfa_fdepot^FORM_EPO_ALFA

    # Endogenous baseline erythropoietin, a population constant (no eta).
    rbase <- exp(lrbase)

    kel <- cl / vc

    # One-compartment model with first-order absorption and linear elimination
    # (Knebel 2008 Results, "Population Pharmacokinetic Modeling Results").
    # Subcutaneous doses go into depot; intravenous doses go directly into
    # central, where neither the bioavailability nor the absorption rate applies.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- fdepot

    # Observed serum erythropoietin is the drug contribution plus the
    # endogenous baseline; Knebel 2008 equation (2) states the model-predicted
    # value "includes the population-predicted baseline erythropoietin
    # concentration (Bepo)", so the residual error applies to the total.
    Cc <- central / vc + rbase
    Cc ~ lnorm(expSd)
  })
}
