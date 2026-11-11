Upton_2025_cycloserine <- function() {
  description <- paste(
    "One-compartment population PK model for cycloserine in plasma and",
    "lumbar cerebrospinal fluid (CSF) after oral terizidone in adults with",
    "pulmonary drug-resistant tuberculosis without tuberculous meningitis",
    "(Upton 2025). Terizidone doses enter as cycloserine equivalents",
    "(two cycloserine molecules per terizidone: 750 mg terizidone = 506 mg",
    "cycloserine). Absorption is a Savic transit chain (mean transit time",
    "1.34 h, 2.38 transit compartments) feeding a first-order depot",
    "(ka 2.29 1/h); clearance 0.54 L/h and central volume 13.9 L at the",
    "cohort-typical fat-free mass of 37.21 kg, allometrically scaled on",
    "fat-free mass with fixed exponents 0.75 / 1, and clearance falls",
    "exponentially with age (-0.0145 per year from age 37). CSF is a",
    "Sheiner-style effect compartment holding a concentration that",
    "equilibrates with plasma at an equilibration half-life of 4.65 h toward",
    "a CSF-to-plasma pseudo-partition coefficient of 0.690 (total",
    "concentrations). Between-subject variability on clearance (25.6%) and",
    "two-occasion between-occasion variability on bioavailability (17.3%),",
    "ka (95.3%) and mean transit time (53.7%); residual error is combined",
    "proportional plus additive, separately for plasma and CSF."
  )
  reference <- paste(
    "Upton CM, Calderin JM, Diacon AH, Boeree MJ, Denti P, Wiesner L,",
    "Kellermann T, McCulloch M, Aarnoutse R, Wasserman S (2025).",
    "Cerebrospinal fluid penetration of cycloserine/terizidone and",
    "clofazimine in patients with pulmonary TB.",
    "Antimicrob Agents Chemother 69(12):e00931-25.",
    "doi:10.1128/aac.00931-25",
    sep = " "
  )
  vignette <- "Upton_2025_cycloserine_clofazimine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state holds a concentration, not an amount: the supplementary
  # CSF control stream reads CE = A(3) directly as the CSF concentration and
  # integrates DADT(3) = KE0_CSF*(PPC_CSF*C2 - A(3)), the Sheiner effect
  # compartment with negligible mass transfer and volume (Supplement S2,
  # 'Effect compartment implementation for CSF concentrations modelling').
  compartmentData <- list(
    depot = list(
      analyte = "cycloserine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "cycloserine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    csf = list(
      analyte = "cycloserine",
      units = "mg/L",
      specimen = "CSF",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric body-size descriptor on clearance (fixed exponent 0.75)",
        "and central volume (fixed exponent 1); FFM gave a larger drop in",
        "objective function than total body weight (dOFV -33.5 vs -20.4,",
        "Results, Cycloserine). Reference 37.21 kg is TVFFM in the",
        "supplementary plasma control stream ('TVFFM = 37.21'), the cohort",
        "typical value (Table S3 footnote c: 'fat-free mass of 37 kg').",
        "The control stream computes FFM with the Janmahasatian et al.",
        "(Clin Pharmacokinet 2005;44:1051-1065) formula from weight, height",
        "in metres and sex: FFM = WHSMAX * HT^2 * WT / (WHS50 * HT^2 + WT)",
        "with WHSMAX = 37.99, WHS50 = 35.98 for women and 42.92, 30.93 for",
        "men. Supply FFM directly."
      ),
      source_name = "FFM"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Exponential effect on clearance, exp(e_age_cl * (AGE - 37)),",
        "used by the authors as a proxy for age-related decline in renal",
        "function because serum creatinine was missing (Results and",
        "Discussion). The centring value 37 years is the control stream's",
        "EXP(THETA(9)*(AGE-37)); the Results text calls it 'the median age",
        "of the cohort', and Table 1 prints a median of 36 (range 20-60).",
        "The equation is used as printed."
      ),
      source_name = "AGE"
    ),
    OCC = list(
      description = "Dosing occasion for between-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions, matching the plasma control stream's IF(OCC==1) /",
        "IF(OCC==2) blocks with ETA(6), ETA(8), ETA(10), ETA(12) declared",
        "$OMEGA BLOCK(1) SAME. An occasion is a dosing event and its",
        "subsequent observations (Supplement S2). In the study each",
        "participant had a pre-dose sample (the previous, self-reported",
        "dose) and a 24-h profile after a directly observed dose, which is",
        "the natural two-occasion split. Any other OCC value switches the",
        "between-occasion variability off."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested as the allometric descriptor and lost to FFM (dOFV -20.4",
        "vs -33.5). Cohort median 52 kg (range 30-73), Table 1."
      )
    ),
    HIV_POS = list(
      description = "HIV-positive indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Screened with antiretroviral therapy, smoking and time on",
        "treatment (Supplement S2) and not retained: 'No other covariates",
        "were found to significantly impact the plasma pharmacokinetics of",
        "cycloserine'. 13 of 27 participants (48%) were living with HIV,",
        "9 (33%) on ART, Table 1."
      )
    ),
    CSF_ALB_RATIO = list(
      description = "CSF-to-serum albumin ratio",
      units = "(ratio)",
      type = "continuous",
      notes = paste(
        "Tested on the CSF parameters with CSF albumin and not retained",
        "('No covariates significantly affected the CSF pharmacokinetic",
        "parameters'). Median 0.005 (range 0.002-0.0097) in 22",
        "participants, Table 1."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    age_range = "20-60 years (median 36)",
    age_median = "36 years",
    weight_range = "30-73 kg (median 52)",
    weight_median = "52 kg",
    height_range = "1.44-1.81 m (median 1.63)",
    bmi_range = "14.9-26.5 kg/m^2 (median 19.2)",
    sex_female_pct = 59,
    race_ethnicity = "Not reported (South African cohort, Cape Town).",
    disease_state = paste(
      "Pulmonary drug-resistant tuberculosis without tuberculous",
      "meningitis; CSF:serum albumin ratio within age-related limits in",
      "every participant (intact blood-CSF barrier). 13 (48%) living with",
      "HIV (median CD4 127 cells/mm^3), 9 (33%) on antiretroviral therapy",
      "(6 protease-inhibitor-based, 3 dolutegravir-based)."
    ),
    dose_range = paste(
      "Terizidone 750 mg once daily (24 of 27) or 500 mg once daily",
      "(3 of 27) for at least 2 weeks before the PK day, equivalent to",
      "506 mg and 338 mg cycloserine. Doses were directly observed on the",
      "PK day; food was not standardised."
    ),
    regions = "South Africa",
    notes = paste(
      "TASK Cape Town observational CSF penetration study. 216 plasma",
      "samples (pre-dose and 2, 4, 6, 8, 10, 14, 24 h post-dose) and one",
      "time-matched lumbar CSF sample per participant at 2, 4, 6 or 8 h",
      "(27 CSF samples); all concentrations above the LLOQ (plasma",
      "0.16 mg/L, CSF 0.2 mg/L). Sequential analysis: the plasma model was",
      "fitted first, individual plasma parameters were then fixed and the",
      "CSF effect-compartment parameters estimated (NONMEM 7.5.1, FOCE-I,",
      "95% CI by SIR). Final estimates are Table S3; control streams are",
      "Supplement S1."
    )
  )

  ini({
    # --- Plasma structural parameters. Final estimates from Table S3
    # (95% CI by SIR). The Supplement S1 plasma control stream prints the
    # same values in $THETA to one more digit for some (CL 0.544, PROP
    # 0.0762, AGE_CL -0.0144); where the two differ beyond rounding, the
    # Table S3 value is used.
    lka <- log(2.29)
    label("First-order absorption rate constant from depot to central (1/h)")    # Table S3 Ka 2.29 (1.15-4.82); S1 $THETA 3 = 2.29
    lcl <- log(0.54)
    label("Clearance at FFM = 37.21 kg and age 37 years (L/h)")                  # Table S3 CL 0.54 (0.488-0.609); S1 $THETA 1 = 0.544
    lvc <- log(13.9)
    label("Central volume of distribution at FFM = 37.21 kg (L)")                # Table S3 Vc 13.9 (13.0-14.9); S1 $THETA 2 = 13.9
    lmtt <- log(1.34)
    label("Mean absorption transit time (h)")                                    # Table S3 MTT 1.34 (0.950-1.85); S1 $THETA 7 = 1.34
    lntr <- log(2.38)
    label("Number of absorption transit compartments (unitless)")                # Table S3 NN 2.38 (1.13-4.46); S1 $THETA 8 = 2.38
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (fraction)")                                     # Table S3 F 1 Fixed; S1 $THETA 4 (1) FIX

    # --- Covariate effects.
    e_age_cl <- -0.0145
    label("Exponential effect of age on clearance (1/year)")                     # Table S3 'Effect of age on CL (exponent)' -0.0145 (-0.0231 to -0.00622); Results -0.0145 1/year
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance (unitless)")        # Supplement S2 'exponents ... fixed to 0.75 and 1'; S1 ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on central volume (unitless)")   # Supplement S2; S1 ALLMV_FFM = (FFM/TVFFM)

    # --- CSF effect-compartment parameters (second, CSF estimation step).
    # Table S3 reports the equilibration half-life; the control stream
    # estimates KE0_CSF and derives HL_CSF = LOG(2)/KE0_CSF, so
    # ke0 = log(2) / 4.65 = 0.149 1/h, the S1 CSF $THETA 4 value.
    lke0 <- log(log(2) / 4.65)
    label("Plasma-to-CSF equilibration rate constant (1/h)")                     # Table S3 t1/2 Plasma-CSF 4.65 h (2.78-7.70); S1 CSF $THETA 4 KE0_CSF = 0.149
    lppc <- log(0.690)
    label("CSF-to-plasma pseudo-partition coefficient (fraction)")               # Table S3 PPC CSF-Plasma 0.690 (0.621-0.759); S1 CSF $THETA 3 = 0.692

    # --- Between-subject and between-occasion variability. The S1 $OMEGA
    # values are log-scale variances; Table S3 prints their square roots as
    # percentages: sqrt(0.0653) = 25.6%, sqrt(0.0299) = 17.3%,
    # sqrt(0.908) = 95.3%, sqrt(0.288) = 53.7%. BSV on V, ka and F is
    # $OMEGA 0 FIX in S1 and is omitted; BOV on CL is likewise 0 FIX.
    etalcl ~ 0.0653
    label("Between-subject variability in clearance (log-scale variance)")       # S1 $OMEGA 1 = 0.0653; Table S3 BSV CL 25.6%

    # Two occasions sharing one variance each ($OMEGA BLOCK(1) SAME).
    etaiov_fdepot_1 ~ 0.0299
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")  # S1 $OMEGA 7 = 0.0299; Table S3 BOV F 17.3%
    etaiov_fdepot_2 ~ fixed(0.0299)                                              # S1 $OMEGA 8 SAME as $OMEGA 7
    etaiov_ka_1 ~ 0.908
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)") # S1 $OMEGA 9 = 0.908; Table S3 BOV Ka 95.3%
    etaiov_ka_2 ~ fixed(0.908)                                                   # S1 $OMEGA 10 SAME as $OMEGA 9
    etaiov_mtt_1 ~ 0.288
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)") # S1 $OMEGA 11 = 0.288; Table S3 BOV MTT 53.7%
    etaiov_mtt_2 ~ fixed(0.288)                                                  # S1 $OMEGA 12 SAME as $OMEGA 11

    # --- Residual error, combined W = sqrt(ADD^2 + PROP^2) with $SIGMA 1 FIX.
    # The S1 $ERROR builds the additive SD as THETA + 0.2 * LLOQ with
    # LLOQ = 0.2 mg/L hard-coded in both control streams. For plasma,
    # Table S3 prints the estimated THETA (0.136, CI 0.0255-0.249), so the
    # total additive SD is 0.136 + 0.2 * 0.2 = 0.176 mg/L. For CSF the
    # THETA is 0 FIX and Table S3 prints the total, 0.2 * 0.2 = 0.04 mg/L
    # (footnote h).
    propSd <- 0.0760
    label("Proportional residual error for plasma cycloserine (fraction)")       # Table S3 plasma proportional 7.60% (6.66-8.81); S1 $THETA 5 = 0.0762
    addSd <- 0.176
    label("Additive residual error for plasma cycloserine (mg/L)")               # Table S3 plasma additive 0.136 (S1 $THETA 6) plus S1 $ERROR CENS_THR*0.2 with LLOQ = 0.2
    propSd_Ccsf <- 0.255
    label("Proportional residual error for CSF cycloserine (fraction)")          # Table S3 CSF proportional 25.5% (20.3-34.9); S1 CSF $THETA 5 = 0.255
    addSd_Ccsf <- fixed(0.04)
    label("Additive residual error for CSF cycloserine (mg/L)")                  # Table S3 CSF additive 4.00x10-2 Fixed, footnote h: 20% of the 0.2 mg/L CSF LLOQ
  })

  model({
    # 1. Occasion indicators for the two-occasion between-occasion
    # variability (S1 IF(OCC==1) / IF(OCC==2) blocks).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # 2. Individual parameters. S1: TVCL = THETA(1)*ALLMCL_FFM*
    # EXP(THETA(9)*(AGE-37)), TVV = THETA(2)*ALLMV_FFM, TVFFM = 37.21.
    cl <- exp(lcl + etalcl) * (FFM / 37.21)^e_ffm_cl * exp(e_age_cl * (AGE - 37))
    vc <- exp(lvc) * (FFM / 37.21)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    ke0 <- exp(lke0)
    ppc <- exp(lppc)
    kel <- cl / vc

    # 3. ODE system. The S1 $DES computes the Savic transit input
    # analytically, TRANSIT = EXP(LOG(BIO*PD*KTR) - GAMLN(NN+1) +
    # NN*LOG(KTR*TEMPO) - KTR*TEMPO) with KTR = (NN+1)/MTT and PD / TEMPO
    # the last dose and the time since it; rxode2's transit() is the same
    # closed form. F1 = 0 in S1 so the dose does not also enter the depot
    # as a bolus; f(depot) <- 0 reproduces that.
    d/dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc
    d/dt(csf) <- ke0 * (ppc * Cc - csf)
    f(depot) <- 0

    # 4. Observations. Doses are cycloserine-equivalent mg, so Cc and Ccsf
    # are cycloserine in mg/L.
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
