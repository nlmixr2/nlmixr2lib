Britz_2019_fluvoxamine_smoking <- function() {
  description <- "One-compartment population PK model for a single 50 mg oral dose of fluvoxamine in 24 healthy volunteers stratified by cigarette smoking (Britz 2019, Supplement S1 section 4, dataset of Spigset 1995). Zero-order absorption of duration D into the central compartment preceded by an absorption lag time ALAG, and first-order elimination; volume and clearance are apparent (Vc/F, CL/F). The single retained covariate is current cigarette smoking, which multiplies CL/F by 1.28 (a 28 percent increase) relative to non-smokers, reflecting induction of CYP1A2. All subjects were CYP2D6 extensive metabolizers. This is the NONMEM population-PK analysis reported alongside, and used to corroborate, the paper's PK-Sim whole-body PBPK model; the PBPK layer itself is not reproduced here."
  reference <- paste(
    "Britz H, Hanke N, Volz AK, Spigset O, Schwab M, Eissing T, Wendl T, Frechen S, Lehr T.",
    "Physiologically-Based Pharmacokinetic Models for CYP1A2 Drug-Drug Interaction Prediction:",
    "A Modeling Network of Fluvoxamine, Theophylline, Caffeine, Rifampicin, and Midazolam.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(5):296-307.",
    "doi:10.1002/psp4.12397.",
    "Population-PK parameters from Supplement S1, Table S11 (column 'Spigset 1995').",
    "Underlying clinical study: Spigset O, Carleborg L, Hedenmalm K, Dahlqvist R.",
    "Effect of cigarette smoking on fluvoxamine pharmacokinetics in humans.",
    "Clin Pharmacol Ther. 1995;58(4):399-403.",
    sep = " "
  )
  vignette <- "Britz_2019_fluvoxamine"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "fluvoxamine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SMOKE = list(
      description = "Current-smoker binary indicator; 1 = current cigarette smoker, 0 = non-smoker",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-smoker)",
      notes = paste(
        "Time-fixed per subject; the study enrolled two fixed strata (12 smokers and 12",
        "non-smokers) rather than following smoking-cessation changes within subjects, so the",
        "two-level SMOKE encoding is exact and the three-level SMOKE_CURRENT / SMOKE_NEVER",
        "pairing is not needed. Britz 2019 Supplement S1 section 4.3 reports the effect as",
        "'Fluvoxamine total clearance was approximately 28% higher in smokers compared to",
        "non-smokers (mean CL = 170 l/h vs. 133 l/h, p-value < 0.001)', which is the Table S11",
        "multiplier 1.28 applied to the typical CL/F of 133 L/h (133 * 1.28 = 170.2 L/h). The",
        "mechanism is induction of CYP1A2: the main text states that smoking is the strongest",
        "known inducer of CYP1A2 and that the companion PBPK model implemented it as a static",
        "1.38-fold increase in CYP1A2 activity. Note that the PBPK and population-PK effect",
        "sizes are not directly comparable -- 1.38 acts on the CYP1A2 enzyme activity alone,",
        "whereas 1.28 acts on total apparent clearance, which also includes the CYP2D6 and",
        "renal pathways.",
        sep = " "
      ),
      source_name = "smoking status (smoker / non-smoker); Britz 2019 Table S11 row 'Smoking on CL'"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    age_median = "36.5 years",
    weight_range = "51-95 kg",
    sex_female_pct = 41.7,
    disease_state = "healthy volunteers",
    dose_range = "50 mg oral single dose (tablet)",
    regions = "Sweden",
    smoking_status = "12 smokers, 12 non-smokers (5 female and 7 male in each group)",
    cyp2d6_phenotype = "all subjects CYP2D6 extensive metabolizers",
    n_observations = 311,
    notes = paste(
      "Britz 2019 Supplement S1 section 4.3 'Study population': 'The population of the study by",
      "Spigset et al. 1995 consisted of 24 young healthy volunteers, whereof 12 were",
      "non-smokers and 12 were smokers, with 5 females and 7 males in each group. All",
      "volunteers were characterized as CYP2D6 extensive metabolizers. Mean age was 36.5 years",
      "and bodyweight ranged from 51 to 95 kg.' Section 4.2 'Dataset': fluvoxamine was given",
      "orally as a single 50 mg dose; blood samples were taken pre-dose and 1, 2, 3, 4, 5, 6,",
      "7, 8, 10, 12, 24, 32 and 48 h after dosing; plasma levels were measured by HPLC with a",
      "lower limit of quantification of 0.5 nmol/L. One sample was missing, leaving 311",
      "concentrations for analysis.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # One-compartment model with zero-order absorption of duration D into the
    # central compartment, delayed by the absorption lag time ALAG, and linear
    # elimination from the central compartment (Britz 2019 Supplement S1
    # section 4.3, 'Population pharmacokinetic model'). Volume and clearance
    # are apparent (divided by the unknown oral bioavailability F), so no
    # separate bioavailability term is estimated.
    #
    # Every value below is from Britz 2019 Supplement S1 Table S11, column
    # 'Spigset 1995'; the bracketed percentages are the relative standard
    # errors printed in that table's adjacent RSE [%] column.
    # ========================================================================
    ld1 <- log(3.51); label("Zero-order absorption duration D (h)") # Table S11 'D (h)' = 3.51 h (RSE 9%); section 4.3 prose 'Zero-order input time differed slightly between both studies with 1.53 hours and 3.51 hours'
    ltlag <- log(1.79); label("Absorption lag time ALAG (h)") # Table S11 'ALAG (h)' = 1.79 h (RSE 6%); section 4.3 prose 'Absorption lag time was estimated at 2.75 hours and 1.79 hours'
    lvc <- log(3030); label("Apparent central volume of distribution Vc/F (L)") # Table S11 'VCentral (l/F)' = 3030 L (RSE 12%); section 4.3 prose 'The volumes of distribution were estimated at 2610 l/F and 3030 l/F'
    lcl <- log(133); label("Apparent clearance CL/F (L/h)") # Table S11 'CL (l/h/F)' = 133 L/h (RSE 18%); section 4.3 prose 'Fluvoxamine was cleared from the systemic circulation with 147 l/h/F and 133 l/h/F'

    # ---- Retained covariate effect on CL/F ---------------------------------
    # Multiplicative on the typical value, applied as a power of the binary
    # indicator so that CL/F is multiplied by 1.28 for current smokers and
    # left unchanged for non-smokers. Section 4.3 confirms the arithmetic:
    # 133 * 1.28 = 170.2 L/h, matching the quoted 'mean CL = 170 l/h vs. 133
    # l/h' and the quoted 28% increase.
    e_smoke_cl <- 1.28; label("Multiplicative effect of current smoking on CL/F (unitless)") # Table S11 'Smoking on CL' = 1.28 (RSE 21%)

    # ========================================================================
    # Inter-individual variability
    # Supplement S1 section 4.2 states that IIVs were 'modelled exponentially',
    # i.e. P_i = P_TV * exp(eta_i), and Table S11 reports each omega as a
    # coefficient of variation in percent, so the internal variance is
    #     omega^2 = log(1 + CV^2).
    # Table S11 reports no off-diagonal covariance, so the etas are diagonal.
    # ========================================================================
    etalvc ~ 0.2475630 # Table S11 'IIV VCentral (%CV)' = 53% (RSE 18%); log(1 + 0.53^2) = 0.2475630
    etalcl ~ 0.2151920 # Table S11 'IIV CL (%CV)' = 49% (RSE 22%); log(1 + 0.49^2) = 0.2151920

    # ========================================================================
    # Residual error: combined proportional plus additive. Section 4.3:
    # 'Residual variability was best described with a combined error model.
    # Although the additive error is very low, it was necessary to adequately
    # describe the data.' Unlike the Spigset 1997 column, this additive term
    # is printed without parentheses and carries an RSE, so it was estimated
    # rather than fixed.
    #
    # Table S11 reports the additive term in nmol/mL while this model reports
    # Cc in ng/mL, so it is converted with the fluvoxamine molecular weight of
    # 318.34 g/mol given in Supplement S1 Table S1b.
    # ========================================================================
    propSd <- 0.49; label("Proportional residual SD (fraction)") # Table S11 'Proportional (%)' = 49% (RSE 18%)
    addSd <- 9.5502e-04; label("Additive residual SD (ng/mL)") # Table S11 'Additive (nmol/ml)' = 3e-6 (RSE 50%); 3e-6 nmol/mL * 318.34 ng/nmol = 9.5502e-04 ng/mL
  })

  model({
    # ---- Individual PK parameters ------------------------------------------
    # SMOKE is a 0/1 indicator, so the power form reduces to a plain
    # multiplication by 1.28 for current smokers and by 1 for non-smokers.
    cl <- exp(lcl + etalcl) * e_smoke_cl^SMOKE
    vc <- exp(lvc + etalvc)
    d1 <- exp(ld1)
    tlag <- exp(ltlag)

    # ---- Micro-constant ----------------------------------------------------
    kel <- cl / vc

    # ---- ODE system --------------------------------------------------------
    # The oral dose enters the central compartment directly as a zero-order
    # input of duration d1 that starts after the lag time tlag; because Vc/F
    # and CL/F are apparent there is no depot compartment and no separate
    # bioavailability term. Dose records must carry rate = -2 so that rxode2
    # uses the modelled duration dur(central) = d1; without it the dose
    # collapses to an instantaneous bolus. Mean absorption time is
    # tlag + d1 / 2 = 3.545 h.
    d/dt(central) <- -kel * central

    dur(central) <- d1
    alag(central) <- tlag

    # ---- Observation and residual error ------------------------------------
    # central is in mg and vc in L, so central/vc is mg/L; multiply by 1000 to
    # report ng/mL, the unit used for the observed and predicted fluvoxamine
    # concentrations in Supplement S1 Table S1d.
    Cc <- (central / vc) * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
