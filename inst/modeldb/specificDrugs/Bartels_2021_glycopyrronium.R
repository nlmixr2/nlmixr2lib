Bartels_2021_glycopyrronium <- function() {
  description <- "Two-compartment population PK model with bolus input for inhaled glycopyrronium in adults and adolescents with asthma receiving the indacaterol/glycopyrronium/mometasone furoate (IND/GLY/MF) fixed-dose combination or glycopyrronium monotherapy via the Breezhaler device (IRIDIUM Phase III study), with fixed allometric body-weight exponents on all clearance and volume terms and grouped-race (Japanese, other) effects on Vc/F (Bartels 2021)"
  reference <- paste(
    "Bartels C, Jain M, Yu J, Tillmann HC, Vaidya S. Population",
    "Pharmacokinetic Analysis of Indacaterol/Glycopyrronium/Mometasone Furoate",
    "After Administration of Combination Therapies Using the Breezhaler()",
    "Device in Patients with Asthma. Eur J Drug Metab Pharmacokinet.",
    "2021;46(4):487-504. doi:10.1007/s13318-021-00689-x.",
    sep = " "
  )
  vignette <- "Bartels_2021_indacaterol_glycopyrronium_mometasone"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")
  # Unit note: doses are in ug and volumes in L, so central / vc is in ug/L
  # (= ng/mL); the factor 1000 in the observation line converts to pg/mL, the
  # unit in which Bartels 2021 reports every concentration and the assay LLOQ
  # (0.250 pg/mL, Sect. 2.2).

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power model normalised to the Table 3 reference of 75 kg (Sect. 2.4 Eq. 1) with the default allometric exponents held fixed for GLY: 0.75 on CL/F and Q/F, 1 on Vc/F and Vp/F (Table 5; Sect. 3.4).",
      source_name = "Body weight"
    ),
    RACE_JAPANESE = list(
      description = "Japanese ethnicity indicator (1 = Japanese, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian/White, the reference level of the Bartels 2021 grouped-race covariate, Table 3 footnote a)",
      notes = "Level 'Japanese' of the three-level grouped-race covariate. Enters Vc/F as exp(-0.65 * RACE_JAPANESE) (Sect. 2.4 Eq. 2 form; Table 5 'Japanese ethnicity on Vc/F' -0.65). Mutually exclusive with RACE_OTHER; both 0 selects the Caucasian/White reference.",
      source_name = "Grouped race (Japanese)"
    ),
    RACE_OTHER = list(
      description = "Grouped-race level 'other': neither Caucasian/White nor Japanese (1 = other, 0 = Caucasian/White or Japanese)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian/White when RACE_JAPANESE is also 0)",
      notes = "The residual level of the Bartels 2021 three-level grouped-race covariate (Table 3 footnote a; Table 4 'Grouped race, Other'). It pools non-Japanese Asian, Black, Native American and other-race patients (68 of the 698 pooled patients). Enters Vc/F as exp(-0.063 * RACE_OTHER); the coefficient is reported only in the Sect. 3.4 text ('-0.063 [%RSE 320%]'), where the authors state it was kept only to allow estimation of the Japanese effect relative to Caucasian patients.",
      source_name = "Grouped race (other)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "glycopyrronium", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "glycopyrronium", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 698L,
    n_studies = 3L,
    age_range = "11-79 years (pooled Table 4 study ranges)",
    weight_range = "33.6-156 kg (pooled Table 4 study ranges; IRIDIUM mean 82 kg, range 44-136 kg)",
    sex_female_pct = 57.0,
    race_ethnicity = c(
      Caucasian = 74.9,
      Japanese = 15.3,
      `Other Asian` = 1.1,
      Black = 1.0,
      `Native American` = 3.3,
      Other = 4.3
    ),
    disease_state = "asthma (inadequately controlled on medium- or high-dose ICS/LABA)",
    dose_range = "glycopyrronium 50 ug once daily by oral inhalation as IND/GLY/MF 150/50/80 or 150/50/160 ug, or glycopyrronium 50 ug monotherapy, via the Breezhaler device",
    regions = "multinational (IRIDIUM multicentre study)",
    notes = "The 698 patients are the pooled pharmacokinetic analysis set of PALLADIUM (n = 273), IRIDIUM (n = 249) and E2201 (n = 176) (Table 4). Glycopyrronium was given only in IRIDIUM (Table 1), so the GLY model was fit to the IRIDIUM subset (249 patients: 140 female, mean age 52.7 years, mean weight 82 kg, 13 Japanese, 17 other grouped race; sparse sampling at -25 min, 2 min, 15 min and 1 h on Days 30 and 86). Race percentages quoted are for the pooled set."
  )

  ini({
    # Structural parameters: Bartels 2021 Table 5, GLY column (Monolix 2018R1,
    # SAEM). Reference patient: 75 kg, Caucasian/White (Table 3).
    lcl <- log(89); label("Apparent clearance CL/F (L/h)") # Table 5 GLY 'CL/F (L/h)' 89 (RSE 4.2%)
    lvc <- log(440); label("Apparent central volume Vc/F (L)") # Table 5 GLY 'Vc/F (L)' 440 (RSE 6.2%)
    lq <- log(350); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 5 GLY 'Q/F (L/h)' 350 (RSE 5.9%)
    # Sect. 2.3: when a two-compartment model could not be fitted, Vp/F was
    # fixed to the value found previously in COPD patients (ref. 17).
    lvp <- fixed(log(1300)); label("Apparent peripheral volume Vp/F (L)") # Table 5 GLY 'Vp/F (L)' 1300 (fixed)

    # Covariate effects: Table 5 GLY column and Sect. 3.4 text.
    e_wt_cl <- fixed(0.75); label("Power exponent of (WT/75) on CL/F (unitless)") # Table 5 GLY 'Body weight on CL/F' 0.75 (fixed)
    e_wt_vc <- fixed(1); label("Power exponent of (WT/75) on Vc/F (unitless)") # Table 5 GLY 'Body weight on Vc/F' 1 (fixed)
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Table 5 GLY 'Body weight on Q/F' 0.75 (fixed)
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Table 5 GLY 'Body weight on Vp/F' 1 (fixed)
    e_race_japanese_vc <- -0.65; label("Log-scale effect of Japanese ethnicity on Vc/F (unitless)") # Table 5 GLY 'Japanese ethnicity on Vc/F' -0.65 (RSE 32%)
    e_race_other_vc <- -0.063; label("Log-scale effect of grouped-race 'other' on Vc/F (unitless)") # Sect. 3.4 text: effect of non-Caucasian/Japanese patients -0.063 (RSE 320%)

    # Between-subject variability: Table 5 reports SDs of the random effects.
    # Variances are SD^2; the CL/F-Vc/F covariance is r * SD_CL * SD_Vc
    # = 0.5 * 0.39 * 0.53 = 0.10335.
    etalcl + etalvc ~ c(0.1521, 0.10335, 0.2809) # Table 5 GLY 'BSV on CL/F' SD 0.39, 'BSV on Vc/F' SD 0.53, 'Correlation between BSV CL/F and Vc/F' 0.5

    # Residual error: proportional only (Table 5).
    propSd <- 0.34; label("Proportional residual error (fraction)") # Table 5 GLY 'Proportional error, b (fraction)' 0.34 (RSE 3.4%)
  })

  model({
    # Individual parameters (Sect. 2.4 Eqs. 1-2; multiplicative exponential BSV,
    # Sect. 2.3).
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc + e_race_japanese_vc * RACE_JAPANESE + e_race_other_vc * RACE_OTHER) * (WT / 75)^e_wt_vc
    q <- exp(lq) * (WT / 75)^e_wt_q
    vp <- exp(lvp) * (WT / 75)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Bolus input (Sect. 3.2: 'a simpler model with bolus administration
    # described the data and did not require estimation of the absorption
    # rate'): each inhaled dose is a bolus on cmt = "central".
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
