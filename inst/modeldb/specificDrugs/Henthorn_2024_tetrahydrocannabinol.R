Henthorn_2024_tetrahydrocannabinol <- function() {
  description <- "Three-compartment population pharmacokinetic model for plasma delta-9-tetrahydrocannabinol (THC) after ad libitum inhalation of high-potency commercial cannabis flower and concentrate products in occasional and daily adult users, with the inhaled dose itself estimated as a fraction (Fi) of a nominal 15 mg fully bioavailable dose"
  reference <- paste(
    "Henthorn TK, Wang GS, Dooley G, Brooks-Russell A, Wrobel J, Limbacher S,",
    "Kosnett M. Dose Estimation Utility in a Population Pharmacokinetic",
    "Analysis of Inhaled Delta-9-Tetrahydrocannabinol Cannabis Market Products",
    "in Occasional and Daily Users. Ther Drug Monit. 2024;46(5):672-680.",
    "doi:10.1097/FTD.0000000000001224. PMCID: PMC11389879.",
    "The slowly equilibrating volume V3 and the elimination clearance CLe were",
    "fixed to the estimates of Sempio C, Huestis MA, Mikulich-Gilbertson SK,",
    "et al. Population pharmacokinetic modeling of plasma",
    "delta-9-tetrahydrocannabinol and an active and inactive metabolite",
    "following controlled smoked cannabis administration. Br J Clin Pharmacol.",
    "2020;86:611-619 (reference 22 of Henthorn 2024; see the Table 3 footnote).",
    sep = " "
  )
  vignette <- "Henthorn_2024_tetrahydrocannabinol"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CANNABIS_DAILY = list(
      description        = "Habitual cannabis use pattern: 1 = daily user, 0 = occasional user",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (occasional user)",
      notes              = paste(
        "Henthorn 2024 Methods (Participants): occasional use is 'an average of at",
        "least 2 days per month and no more than 3 days per week in the 90 days",
        "before enrollment' with typical use of flower products; daily use is daily",
        "consumption of either flower or concentrate products. The final model",
        "collapses the paper's three recruitment groups (occasional flower, daily",
        "flower, daily concentrate) to this two-level indicator: Methods states the",
        "analysis dichotomized pattern of use 'without regard to whether flower or",
        "concentrate was consumed in the dosing protocol'. Time-fixed per subject.",
        "The covariate acts only on the estimated inhaled-dose fraction Fi; it did",
        "not enter any disposition parameter.",
        sep = " "
      ),
      source_name        = "Covariate (daily user)"
    )
  )

  # Screened in the stepwise covariate search and NOT retained in the final
  # model (Henthorn 2024 Results, p. 676: "No demographic covariates (ie, age,
  # body weight, or sex), applied to any of the parameter estimates, reduced
  # the objective function by more than 3.84 and, therefore, did not reach
  # statistical significance in the forward segment of the stepwise covariate
  # selection process. ... Neither weighed dose nor inhalation device used
  # showed a significant reduction in -2*Loglikelihood in the forward step.").
  # Documentation only -- none of these is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate via the power form of Henthorn 2024 equation (2); dOFV < 3.84, not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a continuous covariate via the power form of Henthorn 2024 equation (2); dOFV < 3.84, not retained. Henthorn 2024 names 'body weight' in the Results sentence and 'body mass index' in the Methods list of continuous demographic covariates; BMI is the variable actually tabulated (Table 1)."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate via the exponential form of Henthorn 2024 equation (3); dOFV < 3.84, not retained. Cohort is 65.5% female (Table 1)."
    ),
    DOSE_THC_MG = list(
      description = "Weighed cannabis inhalation dose: product labeled THC concentration times the pre- to post-inhalation product weight difference",
      units       = "mg",
      type        = "continuous",
      notes       = "Screened via Henthorn 2024 equation (4), a power form normalised to the nominal 15 mg dose: theta_TV * (THCwt / 15 mg)^theta_cov. Not significant in the forward step and not retained. Figure 4 shows the weighed dose is a very poor predictor of the model-estimated inhaled dose."
    ),
    FORM_THC_CONCENTRATE = list(
      description = "Cannabis product form / inhalation device consumed on the study day: 1 = concentrate (vape pen, dab), 0 = flower (smoked)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Henthorn 2024 Results: 'Further delineation of usage pattern to include daily concentrate and daily flower use was significant at a lower level and did not further decrease the -2*Loglikelihood, so it was not included in the final model'; the inhalation device likewise showed no significant reduction. The retained CANNABIS_DAILY covariate is explicitly defined without regard to product form."
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "delta-9-tetrahydrocannabinol",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "delta-9-tetrahydrocannabinol",
      units    = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte  = "delta-9-tetrahydrocannabinol",
      units    = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29,
    n_studies      = 1,
    n_observations = 203,
    age_range      = "21-55 years (enrolled); 21-53 years observed",
    age_mean       = "32.8 years (SD 7.6)",
    bmi_mean       = "26.9 kg/m^2 (SD 6.2)",
    sex_female_pct = 65.5,
    race_ethnicity = c(White = 86.2, `Black/African American` = 10.3, `Other/no response` = 3.4),
    ethnicity      = c(`Hispanic/Latino` = 17.2, `Non-Hispanic/Latino` = 75.9, `Declined to respond` = 6.9),
    disease_state  = "healthy adult cannabis users; three recruitment groups: occasional flower (n = 9), daily flower (n = 10), daily concentrate (n = 10)",
    dose_range     = "ad libitum inhalation over a 15-minute session of participant-supplied Colorado commercial-market product; total THC 15-30% w/w for flower and 60-90% w/w for concentrate; mean weighed THC dose 53.0 mg (SD 58.5), mean model-estimated inhaled THC 11.1 mg (SD 11.1)",
    regions        = "United States (Colorado)",
    notes          = paste(
      "Henthorn 2024 Table 1 (demographics and cannabis use history) and Table 2",
      "(observed consumption). 30 adults completed the study; one occasional user",
      "with no detectable blood THC at any draw was dropped, leaving 29.",
      "Sampling ran 135 minutes from the inception of the first inhalation (blood",
      "draws ~1, 5 and 10 minutes after the final inhalation, then 30, 60, 90 and",
      "145 minutes after inception), so no terminal-phase data were collected.",
      "Each participant's pre-inhalation baseline blood THC was subtracted before",
      "analysis. Measured blood THC was divided by the blood:plasma ratio of 0.63",
      "so the model is fit to, and predicts, PLASMA THC. Participants were mostly",
      "young adult White females, which the authors flag as an extrapolation limit.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural disposition -- Henthorn 2024 Table 3, "Parameters of
    # the Final Pharmacokinetic Model". Volumes in L, clearances in
    # L/min, consistent with units$time = "min".
    # ---------------------------------------------------------------
    lvc  <- log(17.9)          ; label("Central volume of distribution VC (L)")                      # Table 3: VC = 17.9 +/- 1.20 L. The Abstract prints 19.9 +/- 1.2 L for the same parameter; Table 3 is the dedicated final-model table, its SEE (1.20) matches the Abstract's (1.2), and every other value in that Abstract sentence (V2 51.6, Q2 1.65, Q3 1.75) matches Table 3 exactly, so the Abstract value is a digit transposition. Operator ruling 2026-08-05: use Table 3. See vignette Errata.
    lvp  <- log(51.6)          ; label("Rapidly equilibrating peripheral volume V2 (L)")             # Table 3: V2 = 51.6 +/- 4.66 L (Abstract agrees: 51.6 +/- 4.7 L)
    lvp2 <- fixed(log(3327))   ; label("Slowly equilibrating peripheral volume V3 (L)")              # Table 3: V3 = 3327 L, marked "*Parameters fixed to values from Sempio et al." -- no terminal-phase data were collected in this 2.25 h study (Discussion)
    lq   <- log(1.65)          ; label("Intercompartmental clearance to peripheral1, Q2 (L/min)")    # Table 3: Q2 = 1.65 +/- 0.14 L/min (Abstract agrees)
    lq2  <- log(1.75)          ; label("Intercompartmental clearance to peripheral2, Q3 (L/min)")    # Table 3: Q3 = 1.75 +/- 0.10 L/min (Abstract agrees)
    lcl  <- fixed(log(0.72))   ; label("Elimination clearance CLe (L/min)")                          # Table 3: CLe = 0.72 L/min, marked "*Parameters fixed to values from Sempio et al."

    # ---------------------------------------------------------------
    # Inhaled-dose fraction. Fi is the fraction of the nominal 15 mg
    # fully bioavailable inhaled THC dose that the participant
    # actually absorbed (Methods, Pharmacokinetic Analyses: "The
    # inhalation dose was arbitrarily assumed, for modeling purposes,
    # to consist of a 60 mg THC joint smoked with 25% efficiency or
    # bioavailability (ie, 15 mg fully bioavailable by inhalation).
    # An individual's ad libitum inhaled dose was, therefore,
    # estimated as a multiple of this assumed dose, Finhaled").
    # Encoded as the canonical bioavailability parameter applied at
    # the central compartment, since Figure 2 shows drug inhaled
    # directly into VC with no absorption compartment.
    # ---------------------------------------------------------------
    lfdepot <- log(0.12)       ; label("Inhaled-dose fraction Fi of the nominal 15 mg bioavailable dose, occasional-user reference (fraction)")  # Table 3: Fi = 0.12 +/- 0.03 (typical value; the occasional-user reference level because Covariate_occasional user = 0 per the Table 3 footnote)

    # ---------------------------------------------------------------
    # Covariate effect. Henthorn 2024 equation (3), the categorical
    # form, is EXPONENTIAL: theta_i = theta_TV * e^(theta_cov) *
    # e^(eta_i). Arithmetic check against the paper's own Discussion:
    # 0.12 * exp(1.79) * 15 mg = 10.78 mg, and the Discussion states
    # "The typical inhaled THC dose for daily users estimated in this
    # study (10.78 mg)" -- an exact match to four significant figures,
    # which confirms the exponential reading. The paper's prose gloss
    # of this as "an approximate 5-fold difference" is loose;
    # exp(1.79) = 5.99.
    # ---------------------------------------------------------------
    e_cannabis_daily_fdepot <- 1.79 ; label("Exponential effect of daily (vs occasional) cannabis use on the inhaled-dose fraction Fi (log scale)")  # Table 3: Covariate (daily user) = 1.79 +/- 0.29. The Discussion prints the same point estimate with a different SEE ("daily = 1.79 +/- 0.67"); only the point estimate is encoded. See vignette Errata.

    # ---------------------------------------------------------------
    # Between-subject variability. Henthorn 2024 equation (1):
    # theta_j = theta_TV * e^(eta_j), so the tabulated omega^2 column
    # is directly the variance of a log-scale eta and is used as
    # printed with no CV%-to-variance conversion. Table 3 leaves the
    # omega^2 column blank for V3, Q3 and CLe, so no IIV is encoded on
    # those three; V3 and CLe were additionally fixed from Sempio et
    # al., while Q3 was estimated as a typical value only.
    # ---------------------------------------------------------------
    etalvc     ~ 0.044   # Table 3: VC omega^2 = 0.044 +/- 0.044 (shrinkage 0.29)
    etalvp     ~ 0.16    # Table 3: V2 omega^2 = 0.16 +/- 0.16 (shrinkage 0.10)
    etalq      ~ 0.13    # Table 3: Q2 omega^2 = 0.13 +/- 0.046 (shrinkage 0.25)
    etalfdepot ~ 0.82    # Table 3: Fi omega^2 = 0.82 +/- 0.02 (shrinkage 0.004). Results p. 676 gives the same point estimate with a different SEE, 0.82 +/- 0.23, and reports omega^2 = 1.92 for the model WITHOUT the usage-pattern covariate, so ~50% of the inhaled-dose variance is explained by CANNABIS_DAILY.

    # ---------------------------------------------------------------
    # Residual error. Table 3 row "d" = 0.15 +/- 0.02, footnoted
    # "d, proportional (relative) intrasubject variability"; Methods
    # confirms "The residual within-subject error was calculated using
    # the relative error method".
    # ---------------------------------------------------------------
    propSd <- 0.15 ; label("Proportional residual error (fraction)")  # Table 3: delta = 0.15 +/- 0.02
  })

  model({
    # ---------------------------------------------------------------
    # 1. Individual parameters. Log-normal per equation (1).
    # ---------------------------------------------------------------
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    vp2 <- exp(lvp2)
    q   <- exp(lq + etalq)
    q2  <- exp(lq2)
    cl  <- exp(lcl)

    # Inhaled-dose fraction with the daily-user covariate applied in
    # the exponential form of equation (3):
    #   Fi_i = Fi_TV * exp(theta_cov * CANNABIS_DAILY) * exp(eta_i)
    # CANNABIS_DAILY = 0 (occasional) reproduces the Table 3 typical
    # value Fi = 0.12; CANNABIS_DAILY = 1 gives 0.12 * exp(1.79) =
    # 0.719, i.e. 10.78 mg of the nominal 15 mg dose.
    fdepot <- exp(lfdepot + e_cannabis_daily_fdepot * CANNABIS_DAILY + etalfdepot)

    # ---------------------------------------------------------------
    # 2. Micro-constants (1/min).
    # ---------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---------------------------------------------------------------
    # 3. Three-compartment mammillary system (Henthorn 2024 Figure 2):
    #    drug is inhaled directly into VC, which exchanges with the
    #    rapidly equilibrating V2 (via Q2) and the slowly
    #    equilibrating V3 (via Q3), and is eliminated from VC by CLe.
    #    Amounts are in mg.
    # ---------------------------------------------------------------
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---------------------------------------------------------------
    # 4. Bioavailability. The dose record carries the nominal 15 mg
    #    fully bioavailable inhaled dose; fdepot scales it to the
    #    individual's estimated inhaled amount.
    # ---------------------------------------------------------------
    f(central) <- fdepot

    # ---------------------------------------------------------------
    # 5. Observation. central is in mg and vc in L, so central / vc is
    #    mg/L; the factor 1000 converts mg/L to the ng/mL declared in
    #    units$concentration and used throughout the paper. Cc is a
    #    PLASMA concentration: measured blood THC was divided by the
    #    blood:plasma ratio of 0.63 before fitting (Methods,
    #    Pharmacokinetic Analyses), so blood THC = Cc * 0.63.
    # ---------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
