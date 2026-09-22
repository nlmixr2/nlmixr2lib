Yang_2026_pixavir <- function() {
  description <- paste(
    "Two-compartment population PK model with lagged first-order absorption for",
    "pixavir, the active metabolite of the oral cap-dependent endonuclease",
    "inhibitor pixavir marboxil (TG-1000), pooled over three Chinese clinical",
    "studies in 423 subjects: 56 healthy adults (phase I single ascending dose",
    "with a food-effect crossover) and 367 adults and adolescents aged 12 years",
    "and older with uncomplicated acute influenza (phase II and phase III).",
    "Body weight enters apparent clearance and apparent central volume as",
    "allometric power functions centred on the 61.30 kg population median, with",
    "both exponents estimated rather than fixed (0.961 on CL/F and 1.69 on",
    "Vc/F). Relative bioavailability is a hyperbolic function of dose,",
    "3.13 / (3.13 + DOSE/40), which reproduces the less-than-dose-proportional",
    "exposure observed over 10-160 mg. Prandial state acts on the absorption",
    "rate constant only, which falls from 0.56 /h fasted to 0.34 /h with a",
    "standard diet and 0.116 /h after a high-fat meal, leaving the extent of",
    "absorption unchanged. Inter-individual variability is carried on ka, CL/F",
    "and Vc/F, with CL/F and Vc/F correlated, and residual error is",
    "proportional."
  )

  reference <- paste(
    "Yang Y, Wang C, Liu X, Li P, Su C, Jiao Z.",
    "Population pharmacokinetics and exposure-response analysis of oral pixavir",
    "marboxil in adults and adolescents with influenza.",
    "Pharmaceutics. 2026;18(5):550.",
    "doi:10.3390/pharmaceutics18050550"
  )

  vignette <- "Yang_2026_pixavir"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "pixavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pixavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pixavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power effects on CL/F and Vc/F, each centred on the",
        "population median of 61.30 kg (the normalising constant printed in",
        "the Yang 2026 final-model equation block, Section 3.1.1). Both",
        "exponents were ESTIMATED, not fixed at the canonical 0.75 / 1.0:",
        "0.961 for CL/F (RSE 8.2%; bootstrap 95% CI 0.805-1.13) and 1.69 for",
        "Vc/F (RSE 8.5%; bootstrap 95% CI 1.43-1.96). The paper reports that",
        "the estimated exponents fitted significantly better than the fixed",
        "canonical pair (p < 0.01) but cautions that they are empirical",
        "approximations valid inside the studied weight range rather than",
        "generalizable physiological constants. Body weight was also formally",
        "tested on Vp/F and Q/F and NOT retained (dOFV < 3.84 for both,",
        "RSE > 50%), so those two parameters carry no size scaling. Cohort",
        "mean 63.8 kg (SD 11.9); the paper's own simulations span 40-120 kg.",
        "Time-fixed at baseline in the source analysis."
      ),
      source_name = "WT"
    ),
    DOSE = list(
      description = "Administered pixavir marboxil dose at the dose record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE canonical: the per-record administered dose",
        "level drives the dose-dependent relative bioavailability",
        "F = 3.13 / (3.13 + DOSE/40) (Yang 2026 Section 3.1.1 displayed",
        "equation, repeated in the final-model equation block). Note that this",
        "hyperbola is NOT anchored to 1 at the 40 mg reference dose: it gives",
        "F = 0.758 at 40 mg, 0.610 at 80 mg and 0.439 at 160 mg, so the",
        "tabulated CL/F, Vc/F, Vp/F and Q/F are the values that pair with this",
        "F, not apparent parameters at a unit-bioavailability reference. The",
        "40 in the denominator is the mg normalisation the paper states in",
        "prose ('DOSE is the administered dose normalized to the 40 mg",
        "reference dose'); DOSE itself is supplied in mg. Doses studied were",
        "single oral 10, 20, 40, 80, 120 and 160 mg, plus two 40 mg doses",
        "48 h apart (Table 1). Set DOSE on every dose record to the amount of",
        "that record; for the 40 mg + 40 mg regimen both records carry 40."
      ),
      source_name = "DOSE"
    ),
    FED = list(
      description = "Fed-versus-fasted state at the dose record, 1 = fed",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "First of the two indicators encoding the three prandial strata Yang",
        "2026 estimated a separate absorption rate constant for (Table 3):",
        "fasted (FED = 0, FED_HIGHFAT = 0, ka = 0.56 /h), standard diet",
        "(FED = 1, FED_HIGHFAT = 0, ka = 0.34 /h) and high-fat meal",
        "(FED = 1, FED_HIGHFAT = 1, ka = 0.116 /h). Food acted on the RATE and",
        "not the EXTENT of absorption - relative bioavailability carries no",
        "food term - so AUCinf is unaffected by prandial state in this model",
        "while Cmax and Tmax are. The 'standard diet' stratum is the ordinary",
        "state in this cohort rather than a special arm: 290 of 423 subjects",
        "(68.6%) are standard-diet, 133 (31.4%) fasted and only 12 (2.8%)",
        "high-fat (Table 2; the high-fat records come from the 12-subject",
        "two-way crossover food-effect part B of TG-1000-C-01, so those",
        "subjects also contribute fasted records). Dose-record level and",
        "therefore time-varying within a crossover subject."
      ),
      source_name = "Food intake"
    ),
    FED_HIGHFAT = list(
      description = "High-fat-meal-at-dosing indicator, 1 = dosed after a high-fat meal",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted or standard diet)",
      notes = paste(
        "Second indicator of the three-level prandial stratification described",
        "under FED; requires FED = 1 when set, matching the Goel 2016",
        "convention recorded in the covariate register. Only the part B",
        "crossover of the phase I study TG-1000-C-01 (n = 12) dosed after a",
        "high-fat meal. The paper reports that a high-fat meal reduced Cmax by",
        "about 41% and AUCinf by about 19% and delayed Tmax by about 1.5 h in",
        "the dedicated food-effect analysis, but that the effect retained in",
        "the popPK model is on ka alone and was not judged clinically",
        "meaningful in the target patient population."
      ),
      source_name = "Food intake"
    )
  )

  # Covariates Yang 2026 screened in the stepwise covariate search (Section
  # 2.2.1) but did not retain in the final model (Section 3.1.1). Documented
  # here so the provenance of the covariate screen survives without carrying
  # covariateData entries that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Demographic covariate screened; no statistically significant effect",
        "on pixavir PK (Section 3.1.1). Cohort mean 27.1 years (SD 8.08);",
        "adolescents mean 15.20 (SD 1.57), adults mean 27.03 (SD 8.41),",
        "Table 4. Adolescents are 12 years and older by protocol."
      )
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Demographic covariate screened; no statistically significant effect",
        "(Section 3.1.1). Cohort 190 of 423 female (44.9%), Table 2."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Evaluated as an alternative body-size descriptor alongside lean body",
        "mass; neither improved model performance nor model stability relative",
        "to total body weight (Section 3.1.1). Cohort mean 22.4 kg/m^2",
        "(SD 2.97), Table 2."
      )
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Evaluated as an alternative body-size descriptor; not retained in",
        "favour of total body weight (Section 3.1.1). No summary statistics",
        "are tabulated for this covariate."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Hepatic-function marker screened; not retained (Section 3.1.1).",
        "Cohort mean 21.4 U/L (SD 18.2), Table 2. Most participants had",
        "normal hepatic function."
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Hepatic-function marker screened; not retained (Section 3.1.1).",
        "Cohort mean 10.8 umol/L (SD 5.30), Table 2."
      )
    ),
    TPRO = list(
      description = "Total protein",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened as a hepatic-function marker; not retained (Section 3.1.1).",
        "Cohort mean 76.0 g/L (SD 5.31), Table 2. Of interest because pixavir",
        "is more than 96% plasma-protein bound."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Renal-function marker screened; not retained (Section 3.1.1). Cohort",
        "mean 131 umol/L (SD 35.3), Table 2. Most participants had normal",
        "renal function. Pixavir is eliminated predominantly via feces after",
        "UGT1A3-mediated glucuronidation, so a renal effect was not expected."
      )
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer versus influenza-patient indicator, 1 = healthy",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Disease status screened as a covariate; not retained (Section",
        "3.1.1 'No statistically significant effects of age, sex, or influenza",
        "virus type'; disease status is listed among the screened covariates",
        "in Section 2.2.1). Cohort 56 healthy (13.2%) and 367 influenza",
        "patients (86.8%), Table 2."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 423L,
    n_studies = 3L,
    n_observations = 3125L,
    age_range = "12 years and older (adolescents 12-17; overall mean 27.1 years, SD 8.08)",
    age_median = "not reported; mean 27.1 years (SD 8.08)",
    weight_range = "not reported; mean 63.8 kg (SD 11.9), model-centring median 61.30 kg",
    weight_median = "61.30 kg",
    sex_female_pct = 44.9,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Pooled: 56 healthy adults (13.2%) from the phase I single-ascending-dose",
      "and food-effect study TG-1000-C-01, and 367 adults and adolescents",
      "(86.8%) with uncomplicated acute influenza from the phase II",
      "dose-ranging study TG-1000-C-02 and the phase III study TG-1000-C-03"
    ),
    renal_function = "Most participants had normal renal function; mean serum creatinine 131 umol/L (SD 35.3)",
    hepatic_function = "Most participants had normal hepatic function; mean ALT 21.4 U/L (SD 18.2), mean total bilirubin 10.8 umol/L (SD 5.30)",
    dose_range = paste(
      "Single oral doses of 10, 20, 40, 80, 120 and 160 mg pixavir marboxil in",
      "phase I; single 40 mg or 80 mg, or two 40 mg doses 48 h apart, in phase",
      "II; single 40 mg (40-80 kg) or 80 mg (80 kg and above) in phase III"
    ),
    regions = "China (all studies conducted in China; all participants Chinese)",
    notes = paste(
      "Baseline demographics are Yang 2026 Table 2 and the study list is",
      "Table 1. Bioanalysis was validated LC-MS/MS with an LLOQ of 2.0 ng/mL",
      "for phase I and II samples and 1.0 ng/mL for phase III samples; BLQ",
      "records were 4.6% of all observations and, falling under the",
      "prespecified 5% threshold, were excluded from model fitting rather than",
      "handled by an M3-type likelihood. Only the prodrug's active metabolite",
      "pixavir was measured - plasma pixavir marboxil is generally below the",
      "LLOQ because hydrolysis by intestinal carboxylesterases is rapid - so",
      "the depot state here holds pixavir marboxil-equivalent amounts and the",
      "model carries no explicit prodrug compartment. Shrinkage was below 20%",
      "for all random effects. The companion exposure-response analysis in the",
      "same paper found no relationship between exposure and either efficacy",
      "or safety; it fits no model beyond descriptive linear regressions and",
      "is reproduced in the vignette rather than packaged as a model file."
    )
  )

  ini({
    # Structural parameters. Every value is the typical estimate for the
    # reference subject of Yang 2026 Table 3 and the final-model equation
    # block in Section 3.1.1: body weight 61.30 kg (the population median
    # used as the allometric centring constant). Because relative
    # bioavailability is NOT 1 at the reference dose (see e_dose_fdepot), the
    # clearance and volume values below pair with that bioavailability term
    # rather than being apparent parameters at F = 1.
    lcl <- log(9.14); label("Clearance CL/F at the reference weight (L/h)")                        # Yang 2026 Table 3: CL/F = 9.14 L/h (RSE 4.0%; bootstrap 9.12, 8.18-10.12)
    lvc <- log(268); label("Central volume of distribution Vc/F at the reference weight (L)")      # Yang 2026 Table 3: Vc/F = 268 L (RSE 4.9%; bootstrap 267, 231-304)
    lvp <- log(101); label("Peripheral volume of distribution Vp/F (L)")                           # Yang 2026 Table 3: Vp/F = 101 L (RSE 6.4%; bootstrap 101, 85.8-117)
    lq  <- log(5.28); label("Intercompartmental clearance Q/F (L/h)")                              # Yang 2026 Table 3: Q/F = 5.28 L/h (RSE 11.9%; bootstrap 5.26, 4.04-6.86)

    # Absorption rate constant. Yang 2026 estimated THREE values of ka in one
    # joint fit, one per prandial stratum, rather than a reference value plus
    # multiplicative food offsets - Table 3 lists three separate "k a (1/h)"
    # rows each with its own RSE and bootstrap interval, and the final-model
    # equation block writes ka as a three-branch case expression. They are
    # therefore encoded under the stratum-suffix convention (every stratum
    # carries a suffix; none keeps the bare canonical lka) so that each
    # published estimate is preserved verbatim and individually traceable.
    lka_fasted  <- log(0.56); label("Absorption rate constant ka, fasted (1/h)")                   # Yang 2026 Table 3: k a for fasting = 0.56 /h (RSE 10.7%; bootstrap 0.559, 0.497-0.63)
    lka_stddiet <- log(0.34); label("Absorption rate constant ka, standard diet (1/h)")            # Yang 2026 Table 3: k a for standard diet = 0.34 /h (RSE 6.1%; bootstrap 0.34, 0.297-0.390)
    lka_highfat <- log(0.116); label("Absorption rate constant ka, high-fat meal (1/h)")           # Yang 2026 Table 3: k a for high fat meal = 0.116 /h (RSE 10.1%; bootstrap 0.116, 0.08-0.154)

    ltlag <- log(0.333); label("Absorption lag time (h)")                                          # Yang 2026 Table 3: ALAG = 0.333 h (RSE 1.4%; bootstrap 0.333, 0.313-0.35)

    # Allometric exponents on body weight, both ESTIMATED rather than fixed at
    # the canonical 0.75 / 1.0 (Section 3.1.1 argues the estimated pair fits
    # significantly better, p < 0.01). No fixed() wrapper: Table 3 reports an
    # RSE and a bootstrap interval for each.
    e_wt_cl <- 0.961; label("Allometric exponent of body weight (/61.30 kg) on CL/F (unitless)")   # Yang 2026 Table 3: CL_WT = 0.961 (RSE 8.2%; bootstrap 0.961, 0.805-1.13)
    e_wt_vc <- 1.69; label("Allometric exponent of body weight (/61.30 kg) on Vc/F (unitless)")    # Yang 2026 Table 3: Vc_WT = 1.69 (RSE 8.5%; bootstrap 1.69, 1.43-1.96)

    # Dose effect on relative bioavailability. Yang 2026 Section 3.1.1
    # displays the hyperbolic form
    #   F = theta / (theta + DOSE/40)
    # and the final-model equation block repeats it with the estimate
    # substituted, F = 3.13 / (3.13 + DOSE/40). This descends from 0.758 at
    # 40 mg through 0.610 at 80 mg to 0.439 at 160 mg, reproducing the
    # less-than-dose-proportional exposure seen over 10-160 mg. The paper is
    # explicit that the function is empirical, not mechanistic: a
    # Michaelis-Menten alternative was tried and rejected for poor
    # identifiability.
    e_dose_fdepot <- 3.13; label("Hyperbolic dose coefficient on relative bioavailability (unitless)")  # Yang 2026 Table 3: F dose = 3.13 (RSE 13.5%; bootstrap 3.11, 2.19-4.73)

    # Inter-individual variability. Section 2.2.1: "Inter-individual
    # variability (IIV) was modeled assuming a log-normal distribution".
    # Table 3 reports the IIV as a percentage, so the percentages are read as
    # log-normal CVs and converted with omega^2 = log(CV^2 + 1):
    #   ka    75.8% -> log(0.758^2 + 1) = 0.4540
    #   CL/F  26.6% -> log(0.266^2 + 1) = 0.06836
    #   Vc/F  39.1% -> log(0.391^2 + 1) = 0.1423
    # The CL/F-Vc/F covariance is tabulated directly on the OMEGA scale as
    # 0.0928 and is used verbatim; with the converted variances it implies a
    # correlation of 0.94, which is what one expects for an oral drug whose
    # unmodelled bioavailability variability is shared by CL/F and Vc/F. The
    # resulting 2x2 block is positive definite (determinant 0.00111).
    # See the vignette's Assumptions section for the alternative reading in
    # which the tabulated percentages are sqrt(omega^2) instead.
    etalcl + etalvc ~ c(0.06836, 0.0928, 0.1423)  # Yang 2026 Table 3: eta CL/F 26.6% (RSE 3.9%), covariance CL/F-Vc/F 0.0928 (RSE 9.4%), eta Vc 39.1% (RSE 5.4%)
    etalka          ~ 0.4540                      # Yang 2026 Table 3: eta k a 75.8% (RSE 4.7%; bootstrap 75.59, 69.21-83.11)

    # Residual error. Section 3.1.1: "RUV was best described by a
    # proportional error model."
    propSd <- 0.2195; label("Proportional residual error SD (fraction)")                           # Yang 2026 Table 3: proportional error = 21.95% (RSE 0.7%; bootstrap 21.86, 19.97-23.56)
  })

  model({
    # ---- 1. Derived covariate terms ---------------------------------------
    # Reference body weight: the 61.30 kg population median printed in the
    # denominator of both allometric terms of the Yang 2026 final-model
    # equation block (Section 3.1.1).
    wtRef <- 61.30

    # Dose-dependent relative bioavailability (Yang 2026 Section 3.1.1). The
    # 40 is the mg reference dose the prose names; DOSE is supplied in mg.
    frel <- e_dose_fdepot / (e_dose_fdepot + DOSE / 40)

    # Prandial state selects one of the three estimated absorption rate
    # constants. Exactly one of the three products below is non-zero for any
    # admissible (FED, FED_HIGHFAT) pair, so this is a case selection written
    # arithmetically rather than a blend.
    lkaTyp <-
      lka_fasted  * (1 - FED) +
      lka_stddiet * FED * (1 - FED_HIGHFAT) +
      lka_highfat * FED * FED_HIGHFAT

    # ---- 2. Individual PK parameters --------------------------------------
    ka <- exp(lkaTyp + etalka)
    cl <- exp(lcl + etalcl) * (WT / wtRef)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wtRef)^e_wt_vc

    # Body weight was tested on Q/F and Vp/F and NOT retained (Section
    # 3.1.1), so neither carries an allometric term.
    q  <- exp(lq)
    vp <- exp(lvp)

    tlag <- exp(ltlag)

    # ---- 3. Micro-constants ----------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ---------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <-
      ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 5. Bioavailability and absorption lag ---------------------------
    f(depot) <- frel
    alag(depot) <- tlag

    # ---- 6. Observation and error model ----------------------------------
    # central is in mg and vc in L, so central/vc is mg/L; the factor 1000
    # converts to the ng/mL scale the paper reports throughout.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
