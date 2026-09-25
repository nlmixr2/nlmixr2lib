Britz_2019_fluvoxamine_cyp2d6 <- function() {
  description <- "One-compartment population PK model for a single 50 mg oral dose of fluvoxamine in 10 healthy volunteers phenotyped for CYP2D6 (Britz 2019, Supplement S1 section 4, dataset of Spigset 1997). Zero-order absorption of duration D into the central compartment preceded by an absorption lag time ALAG, and first-order elimination; volume and clearance are apparent (Vc/F, CL/F). The single retained covariate is the CYP2D6 poor-metabolizer phenotype, which multiplies CL/F by 0.775 (a 22 percent reduction) relative to extensive metabolizers. All subjects were non-smokers. This is the NONMEM population-PK analysis reported alongside, and used to corroborate, the paper's PK-Sim whole-body PBPK model; the PBPK layer itself is not reproduced here."
  reference <- paste(
    "Britz H, Hanke N, Volz AK, Spigset O, Schwab M, Eissing T, Wendl T, Frechen S, Lehr T.",
    "Physiologically-Based Pharmacokinetic Models for CYP1A2 Drug-Drug Interaction Prediction:",
    "A Modeling Network of Fluvoxamine, Theophylline, Caffeine, Rifampicin, and Midazolam.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(5):296-307.",
    "doi:10.1002/psp4.12397.",
    "Population-PK parameters from Supplement S1, Table S11 (column 'Spigset 1997').",
    "Underlying clinical study: Spigset O, Granberg K, Hagg S, Norstrom A, Dahlqvist R.",
    "Relationship between fluvoxamine pharmacokinetics and CYP2D6/CYP2C19 phenotype polymorphisms.",
    "Eur J Clin Pharmacol. 1997;52(2):129-133.",
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
    CYP2D6_PM = list(
      description = "CYP2D6 poor-metabolizer phenotype indicator; 1 = CYP2D6 poor metabolizer, 0 = CYP2D6 extensive metabolizer",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2D6 extensive metabolizer)",
      notes = paste(
        "Time-fixed per subject. Britz 2019 Supplement S1 section 4.2 states that CYP2D6",
        "phenotype was tested both as a continuous covariate (the dextromethorphan metabolic",
        "ratio, on linear and logarithmic scales) and as a categorical covariate, and that the",
        "categorical encoding described the data best; the final model therefore carries a",
        "simple PM-versus-EM binary. Supplement S1 section 4.3 reports the resulting effect as",
        "'Volunteers phenotyped as CYP2D6 poor metabolizers had a 22% lower total clearance",
        "compared to extensive metabolizers (mean CL = 114 l/h vs. 147 l/h, p-value < 0.001)',",
        "which is the Table S11 multiplier 0.775 applied to the typical CL/F of 147 L/h",
        "(147 * 0.775 = 113.9 L/h). The study cohort contained 5 extensive and 5 poor",
        "metabolizers, so both levels are represented. No intermediate or ultrarapid",
        "metabolizers were present, so the paired CYP2D6_IM / CYP2D6_UM indicators are not",
        "needed and the reference category is unambiguously the extensive metabolizer.",
        sep = " "
      ),
      source_name = "CYP2D6 phenotype (extensive / poor metabolizer); Britz 2019 Table S11 row 'CYP2D6 poor metabolism on CL'"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10,
    n_studies = 1,
    age_median = "25 years",
    weight_range = "55-86 kg",
    sex_female_pct = 10,
    disease_state = "healthy volunteers",
    dose_range = "50 mg oral single dose (tablet)",
    regions = "Sweden",
    smoking_status = "all subjects non-smokers",
    cyp2d6_phenotype = "5 extensive metabolizers, 5 poor metabolizers",
    n_observations = 139,
    notes = paste(
      "Britz 2019 Supplement S1 section 4.3 'Study population': 'The population of the study",
      "by Spigset et al. 1997 consisted of 10 young healthy volunteers, that were categorized",
      "by CYP2D6 phenotype. The only female and 4 of the male volunteers were phenotyped as",
      "CYP2D6 extensive metabolizers, whereas 5 males were characterized as CYP2D6 poor",
      "metabolizers. All subjects were non-smokers. Mean age was 25 years and bodyweight ranged",
      "from 55 to 86 kg.' Section 4.2 'Dataset': fluvoxamine was given orally as a single 50 mg",
      "dose; blood samples were taken pre-dose and 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 24, 32 and",
      "48 h after dosing; plasma levels were measured by HPLC with a lower limit of",
      "quantification of 0.5 nmol/L. One sample was missing, leaving 139 concentrations for",
      "analysis.",
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
    # 'Spigset 1997'; the bracketed percentages are the relative standard
    # errors printed in that table's adjacent RSE [%] column.
    # ========================================================================
    ld1 <- log(1.53); label("Zero-order absorption duration D (h)") # Table S11 'D (h)' = 1.53 h (RSE 12%); section 4.3 prose 'Zero-order input time differed slightly between both studies with 1.53 hours and 3.51 hours'
    ltlag <- log(2.75); label("Absorption lag time ALAG (h)") # Table S11 'ALAG (h)' = 2.75 h (RSE 1%); section 4.3 prose 'Absorption lag time was estimated at 2.75 hours and 1.79 hours'
    lvc <- log(2610); label("Apparent central volume of distribution Vc/F (L)") # Table S11 'VCentral (l/F)' = 2610 L (RSE 11%); section 4.3 prose 'The volumes of distribution were estimated at 2610 l/F and 3030 l/F'
    lcl <- log(147); label("Apparent clearance CL/F (L/h)") # Table S11 'CL (l/h/F)' = 147 L/h (RSE 25%); section 4.3 prose 'Fluvoxamine was cleared from the systemic circulation with 147 l/h/F and 133 l/h/F'

    # ---- Retained covariate effect on CL/F ---------------------------------
    # Multiplicative on the typical value, applied as a power of the binary
    # indicator so that CL/F is multiplied by 0.775 for poor metabolizers and
    # left unchanged for extensive metabolizers. Section 4.3 confirms the
    # arithmetic: 147 * 0.775 = 113.9 L/h, matching the quoted 'mean CL = 114
    # l/h vs. 147 l/h' and the quoted 22% reduction.
    e_cyp2d6_pm_cl <- 0.775; label("Multiplicative effect of CYP2D6 poor metabolism on CL/F (unitless)") # Table S11 'CYP2D6 poor metabolism on CL' = 0.775 (RSE 33%)

    # ========================================================================
    # Inter-individual variability
    # Supplement S1 section 4.2 states that IIVs were 'modelled exponentially',
    # i.e. P_i = P_TV * exp(eta_i), and Table S11 reports each omega as a
    # coefficient of variation in percent, so the internal variance is
    #     omega^2 = log(1 + CV^2).
    # Table S11 reports no off-diagonal covariance, so the etas are diagonal.
    # ========================================================================
    etalvc ~ 0.0807501 # Table S11 'IIV VCentral (%CV)' = 29% (RSE 21%); log(1 + 0.29^2) = 0.0807501
    etalcl ~ 0.2151920 # Table S11 'IIV CL (%CV)' = 49% (RSE 17%); log(1 + 0.49^2) = 0.2151920

    # ========================================================================
    # Residual error: combined proportional plus additive. Section 4.3:
    # 'Residual variability was best described with a combined error model.
    # Although the additive error is very low, it was necessary to adequately
    # describe the data.' The additive term is printed in parentheses in
    # Table S11, and the table footnote states 'Parameter values in
    # parentheses were fixed' -- hence fixed() below.
    #
    # Table S11 reports the additive term in nmol/mL while this model reports
    # Cc in ng/mL, so it is converted with the fluvoxamine molecular weight of
    # 318.34 g/mol given in Supplement S1 Table S1b.
    # ========================================================================
    propSd <- 0.34; label("Proportional residual SD (fraction)") # Table S11 'Proportional (%)' = 34% (RSE 13%)
    addSd <- fixed(2.86506e-07); label("Additive residual SD (ng/mL)") # Table S11 'Additive (nmol/ml)' = (9e-10), fixed per the table footnote; 9e-10 nmol/mL * 318.34 ng/nmol = 2.86506e-07 ng/mL
  })

  model({
    # ---- Individual PK parameters ------------------------------------------
    # CYP2D6_PM is a 0/1 indicator, so the power form reduces to a plain
    # multiplication by 0.775 for poor metabolizers and by 1 for extensive
    # metabolizers.
    cl <- exp(lcl + etalcl) * e_cyp2d6_pm_cl^CYP2D6_PM
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
    # tlag + d1 / 2 = 3.515 h.
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
