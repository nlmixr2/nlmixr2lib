Bartels_2021_mometasoneFuroate <- function() {
  description <- "Two-compartment population PK model with mixed (simultaneous) zero-order/first-order absorption for inhaled mometasone furoate in adults and adolescents with asthma receiving the indacaterol/mometasone furoate (IND/MF) or indacaterol/glycopyrronium/mometasone furoate (IND/GLY/MF) fixed-dose combinations, or MF monotherapy, via the Breezhaler device (PALLADIUM, IRIDIUM and E2201 studies), with estimated allometric body-weight exponents on CL/F and Vc/F, baseline FEV1 on CL/F and Vc/F, and formulation (IND/GLY/MF; medium-dose MF strength) and IRIDIUM-study effects on Vc/F and relative bioavailability (Bartels 2021). The MF Twisthaler monotherapy comparator arms are not covered: their Vc/F, F and Vp/F formulation effects are not reported."
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
  # Unit note: doses are nominal MF ug and volumes in L, so central / vc is in
  # ug/L (= ng/mL); the factor 1000 in the observation line converts to pg/mL,
  # the unit in which Bartels 2021 reports every concentration and the assay
  # LLOQ (1.00 pg/mL, Sect. 2.2).

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power model normalised to the Table 3 reference of 75 kg (Sect. 2.4 Eq. 1): estimated exponents on CL/F (0.34) and Vc/F (0.33); fixed allometric exponents 0.75 on Q/F and 1 on Vp/F (Table 5).",
      source_name = "Body weight"
    ),
    FEV1 = list(
      description = "Baseline forced expiratory volume in 1 second (absolute, L)",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline FEV1 (Table 4 'FEV1 at baseline, L'; the spirometric manoeuvre is not stated). Power model normalised to the Table 3 reference of 2 L (Sect. 2.4 Eq. 1), with exponents -0.17 on CL/F and -0.23 on Vc/F (Table 5). The Vc/F term was added after the Wald-test step because of a residual random-effect correlation (Sect. 3.3). The model output is MF concentration, not FEV1, so the bare FEV1 canonical is used.",
      source_name = "Baseline FEV1"
    ),
    FORM_MF_INDGLYMF = list(
      description = "MF delivered as the indacaterol/glycopyrronium/mometasone furoate (IND/GLY/MF) Breezhaler fixed-dose combination (1) versus IND/MF or MF alone via Breezhaler (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IND/MF Breezhaler FDC, or MF monotherapy via Breezhaler)",
      notes = "Two effects. (1) A fixed dose-normalisation factor: Sect. 2.3 introduced 'a fixed multiplicative factor on the bioavailability ... 1.0 as part of IND/MF FDC or when delivered with the Breezhaler device, and 0.5 as part of IND/GLY/MF FDC' because the higher fine-particle mass of IND/GLY/MF makes 80/160 ug equivalent to 160/320 ug in IND/MF (Table 2). The factor is applied here as a divisor of the nominal dose, i.e. a relative bioavailability of 1/0.5 = 2; see the vignette for why the multiplier reading is falsified. (2) An estimated Vc/F effect, exp(-0.32 * FORM_MF_INDGLYMF) (Table 5 'IND/GLY/MF on Vc/F'). IND/GLY/MF was studied only in IRIDIUM, so in the source data FORM_MF_INDGLYMF = 1 implies STUDY_IRIDIUM = 1.",
      source_name = "Formulation (IND/GLY/MF)"
    ),
    FORM_MF_MEDIUM = list(
      description = "Medium-dose MF strength of the Breezhaler FDCs (1 = IND/MF 150/160 ug or IND/GLY/MF 150/50/80 ug; 0 = the high-dose strengths IND/MF 150/320 ug or IND/GLY/MF 150/50/160 ug, or MF monotherapy via Breezhaler)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (high-dose MF FDC strengths, IND/MF 150/320 ug and IND/GLY/MF 150/50/160 ug)",
      notes = "Relative bioavailability of the medium-dose MF formulations versus the corresponding high-dose ones, exp(0.18 * FORM_MF_MEDIUM) (Table 5 'IND/MF and IND/GLY/MF medium-dose MF on F'; Sect. 3.2 third bullet). The effect was estimated only for the medium-dose FDC strengths; the E2201 MF Breezhaler monotherapy arms (80 and 320 ug) are taken as the reference.",
      source_name = "Formulation (medium-dose MF)"
    ),
    STUDY_IRIDIUM = list(
      description = "IRIDIUM study indicator (1 = subject from IRIDIUM, NCT02571777; 0 = PALLADIUM, NCT02554786, or E2201, NCT01555151)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PALLADIUM or E2201)",
      notes = "Study effect on Vc/F, exp(0.19 * STUDY_IRIDIUM) (Table 5 'Study effect on Vc/F in IRIDIUM'). Sect. 3.6 reports that PALLADIUM patients had a 13% higher simulated mean Cmax than IRIDIUM patients.",
      source_name = "Study (IRIDIUM)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mometasone furoate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mometasone furoate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mometasone furoate", units = "ug", specimen = "plasma", verified = TRUE)
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
    disease_state = "asthma (inadequately controlled on ICS or ICS/LABA)",
    dose_range = "MF 80-320 ug once daily by oral inhalation via the Breezhaler device (IND/MF 150/160 or 150/320 ug; IND/GLY/MF 150/50/80 or 150/50/160 ug; MF 80 or 320 ug monotherapy in E2201); the analysis also included MF 200-800 ug per day via the Twisthaler device, which this model does not cover",
    regions = "multinational (PALLADIUM, IRIDIUM, E2201 multicentre studies; 107 Japanese patients)",
    notes = "Pooled pharmacokinetic analysis set of PALLADIUM (n = 273), IRIDIUM (n = 249) and E2201 (n = 176) (Table 4; 398 female). E2201 is a Phase II Breezhaler-versus-Twisthaler device-bridging study with 6 MF samples over 24 h on Days 1 and 28, which supports estimation of the MF distribution parameters; the Phase III studies sampled to 1 h post dose. Baseline FEV1 by study: mean 1.9-2.1 L, range 0.6-4.6 L (Table 4). 35 MF samples were below the LLOQ and were handled as censored in the likelihood (Sect. 3.1)."
  )

  ini({
    # Structural parameters: Bartels 2021 Table 5, MF column (Monolix 2018R1,
    # SAEM). Reference patient: 75 kg, baseline FEV1 2 L, PALLADIUM or E2201,
    # high-dose MF strength of IND/MF via Breezhaler (Table 3; Sect. 3.2).
    lcl <- log(210); label("Apparent clearance CL/F (L/h)") # Table 5 MF 'CL/F (L/h)' 210 (RSE 2.8%)
    lvc <- log(1800); label("Apparent central volume Vc/F (L)") # Table 5 MF 'Vc/F (L)' 1800 (RSE 5.2%)
    lq <- log(250); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 5 MF 'Q/F (L/h)' 250 (RSE 7.4%)
    lvp <- log(3700); label("Apparent peripheral volume Vp/F (L)") # Table 5 MF 'Vp/F (L)' 3700 (RSE 21%)
    lka <- log(2.4); label("First-order absorption rate constant Ka (1/h)") # Table 5 MF 'Ka (1/h)' 2.4 (RSE 9.8%)
    ld1 <- fixed(log(0.01)); label("Duration of the zero-order absorption D (h)") # Table 5 MF 'Duration of zero-order absorption (h)' 0.01 (fixed)
    # Table 5 reports Fr, the fraction absorbed by the zero-order route. The
    # canonical logitffo is the complementary first-order fraction, so
    # logitffo = logit(1 - Fr) = -logit(Fr); see the eta line below.
    logitffo <- qlogis(1 - 0.37); label("Logit of the fraction of the dose absorbed by the first-order route, 1 - Fr (unitless)") # Table 5 MF 'Fraction absorbed via zero-order absorption, Fr' 0.37 (RSE 3.3%)

    # Covariate effects: Table 5 MF column. Continuous covariates are power
    # functions of WT / 75 kg and FEV1 / 2 L (Sect. 2.4 Eq. 1); categorical
    # ones are exp(theta * indicator) (Sect. 2.4 Eq. 2).
    e_wt_cl <- 0.34; label("Power exponent of (WT/75) on CL/F (unitless)") # Table 5 MF 'Body weight on CL/F' 0.34 (RSE 28%)
    e_wt_vc <- 0.33; label("Power exponent of (WT/75) on Vc/F (unitless)") # Table 5 MF 'Body weight on Vc/F' 0.33 (RSE 27%)
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Table 5 MF 'Body weight on Q/F' 0.75 (fixed)
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Table 5 MF 'Body weight on Vp/F' 1 (fixed)
    e_fev1_cl <- -0.17; label("Power exponent of (FEV1/2) on CL/F (unitless)") # Table 5 MF 'Baseline FEV1 on CL/F' -0.17 (RSE 38%)
    e_fev1_vc <- -0.23; label("Power exponent of (FEV1/2) on Vc/F (unitless)") # Table 5 MF 'Baseline FEV1 on Vc/F' -0.23 (RSE 27%)
    e_study_iridium_vc <- 0.19; label("Log-scale effect of the IRIDIUM study on Vc/F (unitless)") # Table 5 MF 'Study effect on Vc/F in IRIDIUM' 0.19 (RSE 30%)
    e_form_mf_indglymf_vc <- -0.32; label("Log-scale effect of the IND/GLY/MF formulation on Vc/F (unitless)") # Table 5 MF 'IND/GLY/MF on Vc/F' -0.32 (RSE 19%)
    # Sect. 2.3 fixed factor 0.5 for IND/GLY/MF, applied as a dose divisor, so
    # the relative bioavailability multiplier is 1 / 0.5 = 2.
    e_form_mf_indglymf_f <- fixed(2); label("Relative bioavailability of MF in IND/GLY/MF versus IND/MF, from the fixed dose-normalisation factor 0.5 (unitless)") # Sect. 2.3 text: fixed multiplicative factor 0.5 as part of IND/GLY/MF FDC (1.0 for IND/MF or Breezhaler)
    e_form_mf_medium_f <- 0.18; label("Log-scale effect of the medium-dose MF strength on relative bioavailability (unitless)") # Table 5 MF 'IND/MF and IND/GLY/MF medium-dose MF on F' 0.18 (RSE 19%)

    # Between-subject variability: Table 5 reports SDs of the random effects.
    # Variances are SD^2; the CL/F-Vc/F covariance is r * SD_CL * SD_Vc
    # = 0.4 * 0.49 * 0.39 = 0.07644.
    etalcl + etalvc ~ c(0.2401, 0.07644, 0.1521) # Table 5 MF 'BSV on CL/F' SD 0.49, 'BSV on Vc/F' SD 0.39, 'Correlation between BSV CL/F and Vc/F' 0.4
    etalq ~ 0.1521 # Table 5 MF 'BSV on Q/F' SD 0.39
    etalvp ~ 4 # Table 5 MF 'BSV on Vp/F' SD 2
    # BSV on Fr, taken as normal on the logit scale (see vignette); the variance
    # on logitffo equals the variance on logit(Fr).
    etalogitffo ~ 0.0484 # Table 5 MF 'BSV on Fr' SD 0.22

    # Residual error: proportional only (Table 5).
    propSd <- 0.36; label("Proportional residual error (fraction)") # Table 5 MF 'Proportional error, b (fraction)' 0.36 (RSE 1.4%)
  })

  model({
    # Individual parameters (Sect. 2.4 Eqs. 1-2; multiplicative exponential BSV,
    # Sect. 2.3).
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl * (FEV1 / 2)^e_fev1_cl
    vc <- exp(lvc + etalvc + e_study_iridium_vc * STUDY_IRIDIUM + e_form_mf_indglymf_vc * FORM_MF_INDGLYMF) *
      (WT / 75)^e_wt_vc * (FEV1 / 2)^e_fev1_vc
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    ka <- exp(lka)
    d1 <- exp(ld1)
    ffo <- expit(logitffo + etalogitffo)

    # Relative bioavailability of the nominal MF dose (Sect. 2.3; Table 5).
    frel <- e_form_mf_indglymf_f^FORM_MF_INDGLYMF * exp(e_form_mf_medium_f * FORM_MF_MEDIUM)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Mixed zero-order/first-order absorption (Sect. 3.2: 'an initial very
    # rapid absorption of a fraction of the drug, overlaid by a slower
    # first-order absorption'). Each inhalation is TWO dose records with the
    # same amt, both at the dosing time: one on cmt = "central" with rate = -2
    # (the zero-order fraction Fr = 1 - ffo, infused over d1) and one on
    # cmt = "depot" (the first-order fraction ffo).
    f(central) <- frel * (1 - ffo)
    dur(central) <- d1
    f(depot) <- frel * ffo

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
