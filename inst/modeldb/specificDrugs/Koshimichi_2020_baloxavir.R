Koshimichi_2020_baloxavir <- function() {
  description <- "Three-compartment population PK model with first-order absorption and an absorption lag time for baloxavir acid (the active form of the prodrug baloxavir marboxil) in healthy adults and in otherwise healthy and high-risk adult and adolescent influenza patients, with bodyweight, race (Asian vs non-Asian) and sex covariates (Koshimichi 2020)"
  reference <- "Koshimichi H, Retout S, Cosson V, Duval V, De Buck S, Tsuda Y, Ishibashi T, Wajima T. Population Pharmacokinetics and Exposure-Response Relationships of Baloxavir Marboxil in Influenza Patients at High Risk of Complications. Antimicrob Agents Chemother. 2020;64(7):e00119-20. doi:10.1128/AAC.00119-20"
  vignette <- "Koshimichi_2020_baloxavir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Doses are administered as baloxavir marboxil (the prodrug) and the model
  # describes plasma baloxavir acid. CL/F and V/F are apparent parameters
  # referenced to the administered baloxavir marboxil dose (no molecular-weight
  # correction appears in the paper), so compartment amounts are expressed in
  # mg of baloxavir marboxil dose equivalents.
  compartmentData <- list(
    depot = list(analyte = "baloxavir marboxil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "baloxavir acid",
      units = "mg (baloxavir marboxil dose equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "baloxavir acid",
      units = "mg (baloxavir marboxil dose equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "baloxavir acid",
      units = "mg (baloxavir marboxil dose equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling centered on 67.7 kg (Koshimichi 2020 Table 2 footnote; Figure S1 reference value). One exponent is shared by CL/F, Q1/F and Q2/F and a second exponent is shared by Vc/F, Vp1/F and Vp2/F.",
      source_name = "body weight"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Multiplicative factor 0.519^Asian on CL/F and 0.564^Asian on Vc/F (Koshimichi 2020 Table 2 footnote: 'Asian = 0 for nonAsian subject and 1 for Asian subject').",
      source_name = "Asian"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative factor 0.682^gender on ka (Koshimichi 2020 Table 2 footnote: 'Gender = 0 for male and 1 for female'); the source coding already matches SEXF.",
      source_name = "gender"
    )
  )

  covariatesDataExcluded <- list(
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Retained in the full model (no. 407) on CL/F but removed in the final model (no. 523) because its effect ratio was close to 1 (Koshimichi 2020 Results; Table S2f; Figure S1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Retained in the full model on CL/F and Vc/F but removed in the final model (Koshimichi 2020 Table S2f; Figure S1)."
    ),
    FED = list(
      description = "Fed dosing condition",
      units = "(binary)",
      type = "binary",
      notes = "Food (fed) effect on F retained in the full model but removed in the final model (Koshimichi 2020 Discussion; Table S2f; Figure S1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1827,
    n_studies = 13,
    n_observations = 11846,
    age_range = "12-85 years",
    weight_range = "36.0-217.3 kg",
    weight_median = "67.7 kg (reference weight of the covariate model)",
    sex_female_pct = 44.3,
    race_ethnicity = c(Asian = 58.8, `Non-Asian` = 41.2),
    disease_state = "277 healthy subjects from 10 phase 1 studies; 1,550 influenza patients from a phase 2 and a phase 3 study in otherwise healthy patients (886) and the CAPSTONE-2 phase 3 study in patients at high risk of influenza complications (664)",
    dose_range = "Single oral doses of baloxavir marboxil 6-80 mg; 40 mg (<80 kg) or 80 mg (>=80 kg) in the phase 3 studies",
    regions = "Japan, United States, United Kingdom and global phase 3 sites",
    notes = "Pooled from 13 clinical studies (Koshimichi 2020 Table S1). Baseline demographics by health status and race in Table 1: healthy Asian n = 231, healthy non-Asian n = 46, patient Asian n = 844, patient non-Asian n = 706. Phase 3 dosing was 40 mg for 40 to <80 kg and 80 mg for >=80 kg."
  )

  ini({
    # Structural parameters - typical values for a 67.7 kg non-Asian male
    lcl   <- log(10.8);  label("Apparent clearance (L/h)")                              # Koshimichi 2020 Table 2: CL/F = 10.8 L/h (RSE 1.8%)
    lvc   <- log(565);   label("Apparent central volume of distribution (L)")           # Koshimichi 2020 Table 2: Vc/F = 565 L (RSE 3.0%)
    lq    <- log(12.4);  label("Apparent first intercompartmental clearance (L/h)")     # Koshimichi 2020 Table 2: Q1/F = 12.4 L/h (RSE 6.9%)
    lvp   <- log(141);   label("Apparent first peripheral volume of distribution (L)")  # Koshimichi 2020 Table 2: Vp1/F = 141 L (RSE 3.1%)
    lq2   <- log(1.43);  label("Apparent second intercompartmental clearance (L/h)")    # Koshimichi 2020 Table 2: Q2/F = 1.43 L/h (RSE 4.3%)
    lvp2  <- log(139);   label("Apparent second peripheral volume of distribution (L)") # Koshimichi 2020 Table 2: Vp2/F = 139 L (RSE 2.4%)
    lka   <- log(1.03);  label("First-order absorption rate constant (1/h)")            # Koshimichi 2020 Table 2: Ka = 1.03 (RSE 6.3%; unit printed as liters/h, 1/h per the footnote)
    ltlag <- log(0.345); label("Absorption lag time (h)")                               # Koshimichi 2020 Table 2: lag time = 0.345 h (RSE 3.3%)

    # Covariate effects
    e_wt_cl_q       <- 0.362; label("Power exponent on (WT/67.7) shared by CL/F, Q1/F and Q2/F (unitless)")    # Koshimichi 2020 Table 2: effect of body wt on CL/F, Q1/F and Q2/F = 0.362 (RSE 10.4%)
    e_wt_vc_vp      <- 0.833; label("Power exponent on (WT/67.7) shared by Vc/F, Vp1/F and Vp2/F (unitless)")  # Koshimichi 2020 Table 2: effect of body wt on Vc/F, Vp1/F and Vp2/F = 0.833 (RSE 5.0%)
    e_race_asian_cl <- 0.519; label("Multiplicative factor on CL/F for Asian vs non-Asian subjects (unitless)")  # Koshimichi 2020 Table 2: effect of race (Asian) on CL/F = 0.519 (RSE 2.2%)
    e_race_asian_vc <- 0.564; label("Multiplicative factor on Vc/F for Asian vs non-Asian subjects (unitless)")  # Koshimichi 2020 Table 2: effect of race (Asian) on Vc/F = 0.564 (RSE 3.6%)
    e_sexf_ka       <- 0.682; label("Multiplicative factor on ka for female vs male subjects (unitless)")        # Koshimichi 2020 Table 2: effect of gender on Ka = 0.682 (RSE 9.1%)

    # Inter-individual variability. Table 2 prints each IIV as '% CV' and the
    # CL/F-Vc/F covariance as a raw value (0.209), so the diagonals are on the
    # same raw OMEGA scale: omega^2 = (CV/100)^2. The printed %RSEs agree with
    # this reading (e.g. Vc/F: bootstrap CI 59.7-65.4% implies an omega^2 RSE of
    # 4.6%, printed 4.6%; the log(1 + CV^2) reading implies 3.9%).
    # CL/F 41.1%, Vc/F 62.7% -> 0.168921, 0.393129; covariance 0.209 (corr 0.81).
    etalcl + etalvc ~ c(0.168921,
                        0.209, 0.393129)  # Koshimichi 2020 Table 2: CV 41.1 and 62.7 percent, covariance 0.209
    etalvp  ~ 0.085849  # Koshimichi 2020 Table 2: Vp1/F IIV CV 29.3 percent
    etalvp2 ~ 0.125316  # Koshimichi 2020 Table 2: Vp2/F IIV CV 35.4 percent
    etalka  ~ 1.530169  # Koshimichi 2020 Table 2: Ka IIV CV 123.7 percent

    # Residual error
    propSd <- 0.202; label("Proportional residual error (fraction)")  # Koshimichi 2020 Table 2: proportional residual error 20.2% CV (RSE 2.1%)
  })
  model({
    # Koshimichi 2020 Table 2 footnote (verbatim structure):
    #   CL/F  = 10.8 * (WT/67.7)^0.362 * 0.519^Asian
    #   Q1/F  = 12.4 * (WT/67.7)^0.362
    #   Q2/F  = 1.43 * (WT/67.7)^0.362
    #   Vc/F  = 565  * (WT/67.7)^0.833 * 0.564^Asian
    #   Vp1/F = 141  * (WT/67.7)^0.833
    #   Vp2/F = 139  * (WT/67.7)^0.833
    #   Ka    = 1.03 * 0.682^gender
    cl   <- exp(lcl + etalcl) * (WT / 67.7)^e_wt_cl_q * e_race_asian_cl^RACE_ASIAN
    vc   <- exp(lvc + etalvc) * (WT / 67.7)^e_wt_vc_vp * e_race_asian_vc^RACE_ASIAN
    q    <- exp(lq) * (WT / 67.7)^e_wt_cl_q
    vp   <- exp(lvp + etalvp) * (WT / 67.7)^e_wt_vc_vp
    q2   <- exp(lq2) * (WT / 67.7)^e_wt_cl_q
    vp2  <- exp(lvp2 + etalvp2) * (WT / 67.7)^e_wt_vc_vp
    ka   <- exp(lka + etalka) * e_sexf_ka^SEXF
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    alag(depot) <- tlag

    # Amounts are in mg and volumes in L, so central/vc is mg/L; x1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
