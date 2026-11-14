Hanke_2026_iclepertin_hemoglobin <- function() {
  description <- paste(
    "Semi-mechanistic population PD model of the effect of steady-state iclepertin (GlyT1",
    "inhibitor) exposure on hemoglobin in healthy adults and in patients with cognitive",
    "impairment associated with schizophrenia (CIAS) or Alzheimer's disease. Two parallel",
    "four-compartment transit chains with a shared transfer rate 4 / lifespan (red blood",
    "cell lifespan fixed to 126 days) carry the mean corpuscular hemoglobin (MCH) and the",
    "red blood cell count (RBC). Iclepertin inhibits the MCH production rate through an",
    "Emax function of AUCtau,ss; total hemoglobin feeds back on the RBC production rate as a",
    "power function. Covariates: age, BMI and Black race on MCH production; sex, age, BMI,",
    "Asian race and ALT on RBC production. Outputs MCH, RBC and total hemoglobin."
  )
  reference <- paste(
    "Hanke N, Callisto SP, Nagy P, Desch M, Abdelwahab MT.",
    "Pharmacokinetic/pharmacodynamic analysis to characterize the effect of long-term",
    "GlyT1 inhibitor iclepertin exposure on hemoglobin levels.",
    "Front Pharmacol. 2026;17:1770980. doi:10.3389/fphar.2026.1770980.",
    "Parameter estimates from Table 2; the structural equations follow the NONMEM control",
    "stream in the article's Supplementary Material, section 4.",
    "The exposure input AUC_ICLEPERTIN can be generated from modellib('Hanke_2026_iclepertin')."
  )
  vignette <- "Hanke_2026_iclepertin"

  # Time is in days (lifespan and production rates are per day). There is no
  # drug compartment: iclepertin enters only through the AUC_ICLEPERTIN
  # covariate, so `dosing` is nominal.
  units <- list(time = "day", dosing = "mg", concentration = "g/L")

  covariateData <- list(
    AUC_ICLEPERTIN = list(
      description = paste(
        "Individual iclepertin steady-state AUC over the 24 h once-daily dosing interval",
        "(AUCtau,ss); 0 for placebo and outside the treatment window."
      ),
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results 'Final model': individual AUCtau,ss from the empirical Bayes estimates of",
        "the popPK model (Supplementary Table S2), applied from the first dose and set to zero",
        "48 h after the last dose. Step-wise time-varying; the user encodes the on/off time",
        "course. Typical CIAS value at 10 mg once daily: 4020 nmol*h/L (Table S3)."
      ),
      source_name = "AUCSS"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "female (SEXF = 1); males take the factor 1.07 on Kin_RBC",
      notes = "Control stream: KinRBC_0 = ... * SexEff**(1 - SEX), dataset SEX 1 = female (Table 2 'Effect of male sex' = 1.07).",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalized power effects (AGE / 59)^theta on Kin_MCH (0.0235) and Kin_RBC (-0.0592); Table 2 and control stream.",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalized power effects (BMI / 26.8)^theta on Kin_MCH (-0.0583) and Kin_RBC (0.0451); Table 2 and control stream.",
      source_name = "BMI"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalized power effect (ALT / 18)^0.0169 on Kin_RBC; Table 2 and control stream.",
      source_name = "ALT"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "non-Black (RACE_BLACK = 0)",
      notes = "Control stream RACEBL = 1 when RACE = 1; factor 0.948 on Kin_MCH (Table 2).",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "non-Asian (RACE_ASIAN = 0)",
      notes = "Control stream RACEASIA = 1 when RACE = 2; factor 0.970 on Kin_RBC (Table 2).",
      source_name = "RACE"
    )
  )

  compartmentData <- list(
    mch1 = list(analyte = "hemoglobin per erythrocyte", units = "pg", specimen = "blood cell", verified = TRUE),
    mch2 = list(analyte = "hemoglobin per erythrocyte", units = "pg", specimen = "blood cell", verified = TRUE),
    mch3 = list(analyte = "hemoglobin per erythrocyte", units = "pg", specimen = "blood cell", verified = TRUE),
    mch4 = list(analyte = "hemoglobin per erythrocyte", units = "pg", specimen = "blood cell", verified = TRUE),
    erythrocytes1 = list(analyte = "erythrocyte", units = "10^12/L", specimen = "blood cell", verified = TRUE),
    erythrocytes2 = list(analyte = "erythrocyte", units = "10^12/L", specimen = "blood cell", verified = TRUE),
    erythrocytes3 = list(analyte = "erythrocyte", units = "10^12/L", specimen = "blood cell", verified = TRUE),
    erythrocytes4 = list(analyte = "erythrocyte", units = "10^12/L", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1116,
    n_studies = 3,
    disease_state = paste(
      "Healthy volunteers (proof-of-mechanism study, n = 25), patients with cognitive",
      "impairment associated with schizophrenia (Phase II, n = 492) and patients with",
      "Alzheimer's disease (Phase II, n = 599)"
    ),
    dose_range = "Placebo or iclepertin 2, 5, 10, 25 or 50 mg once daily for 14 days or 12 weeks",
    covariate_medians = "Covariate reference values (Table 2): age 59 years, BMI 26.8 kg/m^2, ALT 18 U/L",
    notes = paste(
      "Table 1: proof of mechanism 5/10/25/50 mg for 14 days (6/6/8/5 subjects); Phase II",
      "CIAS placebo/2/5/10/25 mg for 12 weeks (170/78/81/83/80); Phase II AD",
      "placebo/2/5/10/25 mg for 12 weeks (118/120/118/122/121). Ten subjects with",
      "CWRES > |4| and out-of-range blood parameters were excluded, leaving 1116 subjects and",
      "12328 MCH and RBC observations (Supplementary Table S6). External validation on 821",
      "CIAS patients on 10 mg once daily for 26 weeks in the Phase III CONNEX-1/-2/-3 trials",
      "(34% female). Hemoglobin was not fitted; it is the product of MCH and RBC.",
      "NONMEM 7.4.3, FOCE-I."
    )
  )

  ini({
    lmtt_rbc <- fixed(log(126)); label("Red blood cell lifespan LS (day)") # Table 2 'Lifespan RBC' = 126 (fixed); control stream THETA(1) 126 FIX
    lkin_mch <- log(0.244); label("MCH production rate Kin_MCH (pg/cell/day)") # Table 2 = 0.244 (RSE 0.183%)
    lkin_rbc <- log(0.0351); label("Baseline RBC production rate Kin_RBC (10^12/L/day)") # Table 2 = 0.0351 (RSE 0.365%)
    lemax <- log(0.552); label("Maximum fractional inhibition of MCH production Emax (unitless)") # Table 2 = 0.552 (RSE 19.3%)
    lauc50 <- log(43500); label("Iclepertin AUCtau,ss giving half-maximal inhibition AUC50 (nmol*h/L)") # Table 2 = 43500 (RSE 22.7%)
    fb <- -1.41; label("Power exponent of the hemoglobin feedback on RBC production Feedbk (unitless)") # Table 2 'Feedback on RBC production rate' = -1.41 (RSE 9.46%)

    e_age_kin_mch <- 0.0235; label("Power exponent of AGE/59 on Kin_MCH (unitless)") # Table 2 = 0.0235 (RSE 20.3%); control stream AgeEff2 = THETA(15)
    e_bmi_kin_mch <- -0.0583; label("Power exponent of BMI/26.8 on Kin_MCH (unitless)") # Table 2 = -0.0583 (RSE 14.9%); control stream BMIEff = THETA(12)
    e_race_black_kin_mch <- 0.948; label("Kin_MCH multiplier for Black race (unitless)") # Table 2 = 0.948 (RSE 0.719%); control stream RaceEff = THETA(10)
    e_sexm_kin_rbc <- 1.07; label("Kin_RBC multiplier for male sex (unitless)") # Table 2 = 1.07 (RSE 0.484%); control stream SexEff = THETA(7)
    e_age_kin_rbc <- -0.0592; label("Power exponent of AGE/59 on Kin_RBC (unitless)") # Table 2 = -0.0592 (RSE 10.7%); control stream AgeEff = THETA(11)
    e_bmi_kin_rbc <- 0.0451; label("Power exponent of BMI/26.8 on Kin_RBC (unitless)") # Table 2 = 0.0451 (RSE 26.5%); control stream BMIEff2 = THETA(13)
    e_race_asian_kin_rbc <- 0.970; label("Kin_RBC multiplier for Asian race (unitless)") # Table 2 = 0.970 (RSE 0.610%); control stream RaceEff2 = THETA(14)
    e_alt_kin_rbc <- 0.0169; label("Power exponent of ALT/18 on Kin_RBC (unitless)") # Table 2 = 0.0169 (RSE 20.0%); control stream ALTEff = THETA(16)

    # IIV: Table 2 reports CV% = 100 * sqrt(exp(omega^2) - 1), so
    # omega^2 = log(1 + CV^2). The lifespan and Emax etas are $OMEGA 0 FIX in
    # the control stream and are omitted.
    etalkin_mch ~ 0.0033123 # Table 2 IIV MCH production rate 5.76 CV%: log(1 + 0.0576^2)
    etalkin_rbc ~ 0.0055943 # Table 2 IIV RBC production rate 7.49 CV%: log(1 + 0.0749^2)
    etalauc50 ~ 0.5240774 # Table 2 IIV AUC50 83.0 CV%: log(1 + 0.830^2)
    etafb ~ 0.4917706 # Table 2 IIV feedback 79.7 CV%: log(1 + 0.797^2)

    # Residual error: Y = IPRED + THETA * EPS(1) with $SIGMA 1 FIX, so the
    # THETAs are additive standard deviations.
    addSd_MCH <- 0.522; label("Additive residual error SD for MCH (pg/cell)") # Table 2 = 0.522 (RSE 1.44%); control stream THETA(8)
    addSd_RBC <- 0.159; label("Additive residual error SD for RBC (10^12/L)") # Table 2 = 0.159 (RSE 1.38%); control stream THETA(9)
  })

  model({
    # Individual parameters (control stream $PK)
    lifespan <- exp(lmtt_rbc)
    kinMch <- exp(lkin_mch + etalkin_mch) * e_race_black_kin_mch^RACE_BLACK *
      (BMI / 26.8)^e_bmi_kin_mch * (AGE / 59)^e_age_kin_mch
    kinRbc0 <- exp(lkin_rbc + etalkin_rbc) * e_sexm_kin_rbc^(1 - SEXF) *
      (AGE / 59)^e_age_kin_rbc * (BMI / 26.8)^e_bmi_kin_rbc *
      e_race_asian_kin_rbc^RACE_ASIAN * (ALT / 18)^e_alt_kin_rbc
    feedback <- fb * exp(etafb)
    emax <- exp(lemax)
    auc50 <- exp(lauc50 + etalauc50)

    # Shared transfer rate of both chains: Ktr = 4 / lifespan
    ktr <- 4 / lifespan

    # Each MCH compartment holds one quarter of the cell's hemoglobin content
    # (the observed MCH is their sum), so the pairwise sum below is one quarter
    # of total hemoglobin. Only its ratio to baseline enters the feedback.
    hbQuarter0 <- 4 * (kinMch / ktr) * (kinRbc0 / ktr)
    hbQuarter <- mch1 * erythrocytes1 + mch2 * erythrocytes2 +
      mch3 * erythrocytes3 + mch4 * erythrocytes4

    # Hemoglobin feedback on RBC production (power function) and iclepertin
    # Emax inhibition of MCH production (control stream $DES)
    kinRbc <- kinRbc0 * (hbQuarter / hbQuarter0)^feedback
    inh <- emax * AUC_ICLEPERTIN / (auc50 + AUC_ICLEPERTIN)

    d/dt(mch1) <- kinMch * (1 - inh) - ktr * mch1
    d/dt(mch2) <- ktr * mch1 - ktr * mch2
    d/dt(mch3) <- ktr * mch2 - ktr * mch3
    d/dt(mch4) <- ktr * mch3 - ktr * mch4

    d/dt(erythrocytes1) <- kinRbc - ktr * erythrocytes1
    d/dt(erythrocytes2) <- ktr * erythrocytes1 - ktr * erythrocytes2
    d/dt(erythrocytes3) <- ktr * erythrocytes2 - ktr * erythrocytes3
    d/dt(erythrocytes4) <- ktr * erythrocytes3 - ktr * erythrocytes4

    # Steady-state initial conditions (control stream A_0)
    mch1(0) <- kinMch / ktr
    mch2(0) <- kinMch / ktr
    mch3(0) <- kinMch / ktr
    mch4(0) <- kinMch / ktr
    erythrocytes1(0) <- kinRbc0 / ktr
    erythrocytes2(0) <- kinRbc0 / ktr
    erythrocytes3(0) <- kinRbc0 / ktr
    erythrocytes4(0) <- kinRbc0 / ktr

    # Observations: MCH (pg/cell) = sum of the MCH chain, RBC (10^12/L) = sum
    # of the RBC chain. Total hemoglobin (g/L) is Figure 1's
    # sum_i(MCH_i * RBC_i) with MCH_i the per-cell content of cohort i
    # (4 x mch_i); it equals MCH x RBC at steady state.
    MCH <- mch1 + mch2 + mch3 + mch4
    RBC <- erythrocytes1 + erythrocytes2 + erythrocytes3 + erythrocytes4
    thb <- 4 * hbQuarter

    MCH ~ add(addSd_MCH)
    RBC ~ add(addSd_RBC)
  })
}
