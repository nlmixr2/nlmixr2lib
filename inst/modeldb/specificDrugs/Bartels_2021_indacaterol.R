Bartels_2021_indacaterol <- function() {
  description <- "Two-compartment population PK model with sequential zero-order/first-order absorption for inhaled indacaterol in adults and adolescents with asthma receiving the indacaterol/mometasone furoate (IND/MF) or indacaterol/glycopyrronium/mometasone furoate (IND/GLY/MF) fixed-dose combinations via the Breezhaler device (PALLADIUM and IRIDIUM Phase III studies), with estimated allometric body-weight exponents on CL/F and Vc/F, fixed allometric exponents on Q/F and Vp/F, a Japanese-ethnicity effect on Vc/F and an IRIDIUM study effect on Vc/F (Bartels 2021)"
  reference <- "Bartels C, Jain M, Yu J, Tillmann HC, Vaidya S. Population Pharmacokinetic Analysis of Indacaterol/Glycopyrronium/Mometasone Furoate After Administration of Combination Therapies Using the Breezhaler Device in Patients with Asthma. Eur J Drug Metab Pharmacokinet. 2021;46(4):489-506. doi:10.1007/s13318-021-00689-x"
  vignette <- "Bartels_2021_indacaterol_glycopyrronium_mometasone"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")
  # Unit note: doses are in ug and volumes in L, so central / vc is in ug/L
  # (= ng/mL); the factor 1000 in the observation line converts to pg/mL, the
  # unit in which Bartels 2021 reports every concentration (Online Resource 2,
  # Figs. 1, 2 and 4) and the assay LLOQ (5.00 pg/mL, Sect. 2.2).

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power model normalised to the Table 3 reference of 75 kg (Sect. 2.4 Eq. 1): estimated exponents on CL/F (0.28) and Vc/F (0.43); fixed allometric exponents 0.75 on Q/F and 1 on Vp/F (Table 5).",
      source_name = "Body weight"
    ),
    RACE_JAPANESE = list(
      description = "Japanese ethnicity indicator (1 = Japanese, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian/White, the reference level of the Bartels 2021 grouped-race covariate, Table 3 footnote a)",
      notes = "Level 'Japanese' of the three-level grouped-race covariate (Caucasian/White, Japanese, other). Enters Vc/F as exp(-0.29 * RACE_JAPANESE) (Sect. 2.4 Eq. 2 form; Table 5 'Japanese ethnicity on Vc/F'). The third level 'other' was also retained on Vc/F in the final IND model (Sect. 3.4) but its coefficient is not reported; see covariatesDataExcluded.",
      source_name = "Grouped race (Japanese)"
    ),
    STUDY_IRIDIUM = list(
      description = "IRIDIUM study indicator (1 = subject from IRIDIUM, NCT02571777; 0 = PALLADIUM, NCT02554786)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PALLADIUM)",
      notes = "Study effect on Vc/F, exp(0.25 * STUDY_IRIDIUM) (Table 5 'Study effect on Vc/F in IRIDIUM'); a larger IRIDIUM Vc/F lowers the IRIDIUM Cmax, matching the higher PALLADIUM Cmax noted in Sect. 3.1. Set to 0 to simulate the PALLADIUM reference.",
      source_name = "Study (IRIDIUM)"
    )
  )

  covariatesDataExcluded <- list(
    RACE_OTHER = list(
      description = "Grouped-race level 'other' (neither Caucasian/White nor Japanese)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian/White)",
      notes = "Retained on Vc/F in the final IND model (Sect. 3.4: 'grouped race (Caucasian/White, Japanese, other) on Vc/F (IND)'), but Table 5 and the text report only the Japanese coefficient. With no reported value the effect is not encoded, so subjects in this group are simulated as the Caucasian/White reference. Sect. 3.6 reports that such patients had a simulated mean Cmax only 5% above Caucasian patients, so the omission is small."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "indacaterol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "indacaterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "indacaterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 698L,
    n_studies = 3L,
    age_range = "11-79 years (pooled Table 4 study ranges)",
    weight_range = "33.6-156 kg (pooled Table 4 study ranges; study means 74.7-82 kg)",
    sex_female_pct = 57.0,
    race_ethnicity = c(
      Caucasian = 74.9,
      Japanese = 15.3,
      `Other Asian` = 1.1,
      Black = 1.0,
      `Native American` = 3.3,
      Other = 4.3
    ),
    disease_state = "asthma (inadequately controlled on medium- or high-dose ICS or ICS/LABA)",
    dose_range = "indacaterol 150 ug once daily by oral inhalation as IND/MF 150/160 or 150/320 ug, or IND/GLY/MF 150/50/80 or 150/50/160 ug, via the Breezhaler device",
    regions = "multinational (PALLADIUM, IRIDIUM, E2201 multicentre studies; 107 Japanese patients)",
    notes = "The 698 patients are the pooled pharmacokinetic analysis set of PALLADIUM (n = 273), IRIDIUM (n = 249) and E2201 (n = 176) (Table 4; 398 female). E2201 contributed mometasone furoate only, so the indacaterol model was fit to PALLADIUM and IRIDIUM data (sparse sampling to 1 h post dose on Days 30 and 84/86, Table 1). Race percentages are from the pooled Table 4 counts of the three analysis studies (Caucasian/White 523, Asian 115 of whom 107 Japanese, Black 7, Native American 23, Other 30); the grouped-race level 'other' is therefore 68 patients (9.7%). Mean baseline FEV1 1.9-2.1 L and mean eGFR 84-96 mL/min/1.73 m2 by study."
  )

  ini({
    # Structural parameters: Bartels 2021 Table 5, IND column (Monolix 2018R1,
    # SAEM). Reference patient: 75 kg, Caucasian/White, PALLADIUM (Table 3).
    lcl <- log(54); label("Apparent clearance CL/F (L/h)") # Table 5 IND 'CL/F (L/h)' 54 (RSE 3%)
    lvc <- log(600); label("Apparent central volume Vc/F (L)") # Table 5 IND 'Vc/F (L)' 600 (RSE 4.7%)
    lq <- log(380); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 5 IND 'Q/F (L/h)' 380 (RSE 3.5%)
    lvp <- log(5700); label("Apparent peripheral volume Vp/F (L)") # Table 5 IND 'Vp/F (L)' 5700 (RSE 15%)
    lka <- fixed(log(50)); label("First-order absorption rate constant Ka (1/h)") # Table 5 IND 'Ka (1/h)' 50 (fixed)
    ld1 <- log(0.055); label("Duration of the zero-order absorption D (h)") # Table 5 IND 'Duration of zero-order absorption (h)' 0.055 (RSE 8.2%)
    # Table 5 reports Fr, the fraction absorbed by the zero-order route. The
    # canonical logitffo is the complementary first-order fraction, so
    # logitffo = logit(1 - Fr) = -logit(Fr); see the eta line below.
    logitffo <- qlogis(1 - 0.33); label("Logit of the fraction of the dose absorbed by the first-order route, 1 - Fr (unitless)") # Table 5 IND 'Fraction absorbed via zero-order absorption, Fr' 0.33 (RSE 15%)

    # Covariate effects: Table 5 IND column. Continuous covariates are power
    # functions of WT / 75 kg (Sect. 2.4 Eq. 1); categorical ones are
    # exp(theta * indicator) (Sect. 2.4 Eq. 2).
    e_wt_cl <- 0.28; label("Power exponent of (WT/75) on CL/F (unitless)") # Table 5 IND 'Body weight on CL/F' 0.28 (RSE 46%)
    e_wt_vc <- 0.43; label("Power exponent of (WT/75) on Vc/F (unitless)") # Table 5 IND 'Body weight on Vc/F' 0.43 (RSE 29%)
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Table 5 IND 'Body weight on Q/F' 0.75 (fixed)
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Table 5 IND 'Body weight on Vp/F' 1 (fixed)
    e_race_japanese_vc <- -0.29; label("Log-scale effect of Japanese ethnicity on Vc/F (unitless)") # Table 5 IND 'Japanese ethnicity on Vc/F' -0.29 (RSE 32%)
    e_study_iridium_vc <- 0.25; label("Log-scale effect of the IRIDIUM study on Vc/F (unitless)") # Table 5 IND 'Study effect on Vc/F in IRIDIUM' 0.25 (RSE 21%)

    # Between-subject variability: Table 5 reports SDs of the random effects.
    # Variances are SD^2; the CL/F-Vc/F covariance is r * SD_CL * SD_Vc
    # = 0.54 * 0.48 * 0.38 = 0.098496.
    etalcl + etalvc ~ c(0.2304, 0.098496, 0.1444) # Table 5 IND 'BSV on CL/F' SD 0.48, 'BSV on Vc/F' SD 0.38, 'Correlation between BSV CL/F and Vc/F' 0.54
    etalq ~ 0.0225 # Table 5 IND 'BSV on Q/F' SD 0.15
    etalvp ~ 1.69 # Table 5 IND 'BSV on Vp/F' SD 1.3
    # BSV on Fr, taken as normal on the logit scale (see vignette). Because
    # logit(1 - Fr) = -logit(Fr) and the eta is symmetric, the variance on
    # logitffo equals the variance on logit(Fr).
    etalogitffo ~ 1.44 # Table 5 IND 'BSV on Fr' SD 1.2

    # Residual error: Table 5 lists only the proportional component for the
    # final IND model (the base-model additive term, Online Resource 4, was
    # dropped; Sect. 2.3).
    propSd <- 0.24; label("Proportional residual error (fraction)") # Table 5 IND 'Proportional error, b (fraction)' 0.24 (RSE 2.4%)
  })

  model({
    # Individual parameters (Sect. 2.4 Eqs. 1-2; multiplicative exponential BSV,
    # Sect. 2.3).
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc + e_race_japanese_vc * RACE_JAPANESE + e_study_iridium_vc * STUDY_IRIDIUM) * (WT / 75)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    ka <- exp(lka)
    d1 <- exp(ld1)
    ffo <- expit(logitffo + etalogitffo)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Sequential zero-order/first-order absorption (Sect. 3.2: 'a short
    # zero-order absorption of a fraction of the drug followed by a rapid
    # first-order absorption of the rest'). Each inhalation is TWO dose records
    # with the same amt: one on cmt = "central" with rate = -2 (the zero-order
    # fraction Fr = 1 - ffo, infused over d1) and one on cmt = "depot" (the
    # first-order fraction ffo, which starts once the zero-order input ends).
    f(central) <- 1 - ffo
    dur(central) <- d1
    f(depot) <- ffo
    alag(depot) <- d1

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
