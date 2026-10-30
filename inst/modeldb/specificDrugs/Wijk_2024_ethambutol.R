Wijk_2024_ethambutol <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ethambutol at steady state",
    "in South African adults without HIV starting first-line treatment for",
    "drug-susceptible pulmonary tuberculosis (TRUST cohort, Wijk 2024).",
    "Absorption is a Savic analytical transit chain (NN = 2.01, i.e. 3.01",
    "transit compartments as Table 2 counts them, mean transit time 0.451 h)",
    "feeding a first-order absorption compartment (ka 1.39 1/h), with",
    "first-order elimination (CL 48.8 L/h), central volume 233 L, peripheral",
    "volume 510 L and intercompartmental clearance 34.4 L/h, allometrically",
    "scaled on fat-free mass with exponents 0.75 / 1. Clearance rises 0.264%",
    "per mL/min of Cockcroft-Gault creatinine clearance standardised to",
    "52 kg, centred on 127.1 mL/min. Doses are in mg of ethambutol base",
    "(0.74 times the mg of ethambutol dihydrochloride). Random effects are",
    "between-subject variability on CL (16.3%) and two-occasion",
    "between-occasion variability on bioavailability (19.8%), ka (112%) and",
    "mean transit time (94.2%), all omega standard deviations on the log",
    "scale. Residual error is combined proportional (10.5%) plus additive",
    "(0.0169 mg/L)."
  )
  reference <- paste(
    "Wijk M, Gausi K, Malatesta S, Weber SE, Court R, Myers B, Carney T,",
    "Parry CDH, Horsburgh CR, White LF, Wiesner L, Warren RM, Uren C,",
    "McIlleron H, Kloprogge F, Denti P, Jacobson KR. The impact of alcohol",
    "and illicit substance use on the pharmacokinetics of first-line TB",
    "drugs. J Antimicrob Chemother. 2024;79(8):2022.",
    "doi:10.1093/jac/dkae206.",
    "Parameter estimates from Table 2; model equations, the transit",
    "compartment count and the creatinine-clearance slope from the NONMEM",
    "control stream reproduced in the supplementary material.",
    sep = " "
  )
  vignette <- "Wijk_2024_first_line_tb_drugs"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass, computed from sex, total body weight and height by the Janmahasatian (2005) formula",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on CL, Vc, Vp and Q with exponents",
        "0.75 / 1, normalised to TVFFM = 43.4 kg as in the control stream",
        "(Table 2 footnote b names the 43.1 kg cohort median). Fat-free mass",
        "fitted better than total body weight (dOFV -30.0 versus -22.0).",
        "Cohort median 43.1 kg (IQR 35.8-47.0), Table 1."
      ),
      source_name = "FFM"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the Cockcroft-Gault creatinine clearance on CL,",
        "(140 - AGE) * WT / (72 * SCr in mg/dL), as in the control stream.",
        "Cohort median 37 years (IQR 27-48), Table 1."
      ),
      source_name = "AGE"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Control stream: CREATININE_MGDL = CREATININE/88.4; GFR =",
        "((140-AGE)*WT)/(CREATININE_MGDL*72) for males and females alike;",
        "GFR_STD = GFR*52/WT; TVGFR_STD = 1 + THETA(11)*(GFR_STD - 127.1).",
        "Weight cancels, so the standardised clearance is",
        "(140 - AGE) * 52 / (72 * CREAT / 88.4) mL/min and depends on age",
        "and serum creatinine only. Serum creatinine is not tabulated; the",
        "Table 1 Cockcroft-Gault creatinine clearance is 128 mL/min (IQR",
        "105-150). Source column CREATININE."
      ),
      source_name = "CREATININE"
    ),
    OCC = list(
      description = "Occasion index (1 or 2) for the between-occasion variability on bioavailability, ka and mean transit time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each dose and its following samples are one occasion (Methods); the",
        "control stream multiplexes the etas over OCC 1 and 2, each pair",
        "sharing one variance. Occasion 1 is the dose before the PK visit",
        "with its pre-dose sample and occasion 2 the observed dose on the PK",
        "visit. Any other OCC value switches the between-occasion variability",
        "off."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    ALCOHOL_USE_HIGH = list(
      description = "High alcohol use group from a latent class analysis of AUDIT, Timeline Follow-Back and phosphatidylethanol.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low or moderate alcohol use)",
      notes = paste(
        "The Results and Abstract report 13% higher ethambutol clearance in",
        "the high alcohol use group (dOFV -9.84, df = 1, P = 0.002). The",
        "effect is not part of the final model as published: the final",
        "control stream has no alcohol term although ALCOHOL_USE is in its",
        "$INPUT, Table 2 lists no alcohol row, and P = 0.002 does not meet",
        "the P < 0.001 backward-elimination criterion of the Methods. 29 of",
        "104 participants (28%) were in the high-use group."
      )
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the allometric descriptor and rejected in favour of fat-free mass; the control stream's Cockcroft-Gault clearance uses it but standardises it away (GFR_STD = GFR*52/WT). Cohort median 51.4 kg (IQR 45.5-56.9), Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 104L,
    n_studies = 1L,
    age_range = "median 37 years (IQR 27-48); 15 years or older",
    weight_range = "median 51.4 kg (IQR 45.5-56.9)",
    ffm_range = "median 43.1 kg (IQR 35.8-47.0)",
    sex_female_pct = 30,
    race_ethnicity = "Not reported; participants were enrolled in Worcester, South Africa",
    disease_state = "Drug-susceptible pulmonary tuberculosis, without HIV, starting first-line treatment; 42% moderate and 28% high alcohol use",
    renal_function = "Cockcroft-Gault creatinine clearance median 128 mL/min (IQR 105-150)",
    dose_range = paste(
      "Fixed-dose combination tablets once daily (ethambutol 275 mg per",
      "tablet, with rifampicin, isoniazid and pyrazinamide) with the number",
      "of tablets by WHO weight band; sampled at steady state 4-8 weeks into",
      "treatment under fasting conditions. The analysis expressed doses as",
      "ethambutol base, 0.74 times the dihydrochloride amount."
    ),
    regions = "South Africa (Worcester, Western Cape)",
    notes = paste(
      "Pharmacokinetic analysis nested in the TRUST prospective cohort.",
      "Samples pre-dose and 1.5, 3, 5 and 8 h post-dose; three ethambutol",
      "concentrations were below the 0.0844 mg/L LLOQ."
    )
  )

  ini({
    # Structural values are Wijk 2024 Table 2. The supplementary control
    # stream's $THETA block carries the same final values to more digits
    # (CL 48.7658, V2 232.883, KA 1.38704, MTT 0.451027, V3 510.098,
    # Q 34.3585, PROP 0.104762) and supplies the structure.
    lcl <- log(48.8)
    label("Clearance at the reference fat-free mass of 43.4 kg and standardised creatinine clearance of 127.1 mL/min (L/h)") # Table 2 CL = 48.8 L/h (44.9-52.9); control stream THETA(1) = 48.7658
    lvc <- log(233)
    label("Central volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vc = 233 L (204-261); control stream THETA(2) = 232.883
    lvp <- log(510)
    label("Peripheral volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vp = 510 L (368-731); control stream THETA(8) = 510.098
    lq <- log(34.4)
    label("Intercompartmental clearance at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 Q = 34.4 L/h (27.8-42.3); control stream THETA(9) = 34.3585
    lka <- log(1.39)
    label("First-order absorption rate constant (1/h)") # Table 2 Ka = 1.39 1/h (0.96-1.91); control stream THETA(3) = 1.38704
    lmtt <- log(0.451)
    label("Mean transit time through the absorption transit chain (h)") # Table 2 MTT = 0.451 h (0.422-0.483); control stream THETA(7) = 0.451027
    # Table 2 prints NN = 3.01 (2.07-4.00); the control stream estimate that
    # enters KTR = (NN+1)/MTT is THETA(10) = 2.01474, i.e. Table 2 reports
    # NN + 1 (see Wijk_2024_pyrazinamide.R). The control stream value is used.
    lntr <- log(2.01474)
    label("Shape parameter NN of the Savic absorption transit chain (unitless)") # control stream THETA(10) = 2.01474, KTR = (NN+1)/MTT; Table 2 NN = 3.01 = NN + 1
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (unitless)") # Table 2 F = 1, held constant; control stream THETA(4) = 1 FIX
    e_crcl_cl <- 0.00263624
    label("Fractional change in clearance per mL/min of standardised creatinine clearance (per mL/min)") # control stream THETA(11) = 0.00263624; Results 'each 10 units increase in creatinine clearance from the median leading to a 2.6% increase in clearance'

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance terms (unitless)") # Methods 'allometric exponents for clearances and volumes fixed to 0.75 and 1'; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume terms (unitless)") # Methods; control stream ALLMV_FFM = (FFM/TVFFM)

    # Table 2 percentages are omega standard deviations on the log scale;
    # sqrt() of each control stream $OMEGA reproduces the tabulated value.
    etalcl ~ 0.0266888
    label("Between-subject variability in CL (log-scale variance)") # Table 2 BSV CL = 16.3% (12.4-21.1); control stream $OMEGA 0.0266888
    etaiov_fdepot_1 ~ 0.039389
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)") # Table 2 BOV BIO = 19.8% (15.3-25.5); control stream $OMEGA 0.039389
    etaiov_fdepot_2 ~ fix(0.039389) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 1.2486
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)") # Table 2 BOV KA = 112% (88.1-138); control stream $OMEGA 1.2486
    etaiov_ka_2 ~ fix(1.2486) # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.887897
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)") # Table 2 BOV MTT = 94.2% (71.6-117); control stream $OMEGA 0.887897
    etaiov_mtt_2 ~ fix(0.887897) # control stream $OMEGA BLOCK(1) SAME

    # Control stream $ERROR: ADD = THETA(6) + 0.2 * LLOQ with THETA(6) = 0 FIX
    # and LLOQ = 0.0844 mg/L.
    propSd <- 0.105
    label("Proportional residual error (fraction)") # Table 2 proportional error = 10.5% (8.90-12.8); control stream THETA(5) = 0.104762
    addSd <- fixed(0.01688)
    label("Additive residual error (mg/L)") # Table 2 additive error 0.02 mg/L held constant; control stream ADD = 0 + 0.0844 * 0.2 = 0.01688
  })
  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # Cockcroft-Gault creatinine clearance without a sex factor, standardised
    # to 52 kg (control stream GFR_STD = GFR*52/WT, so weight cancels);
    # serum creatinine converted from umol/L to mg/dL by 88.4.
    crcl_std <- (140 - AGE) * 52 / (72 * CREAT / 88.4)
    f_crcl <- 1 + e_crcl_cl * (crcl_std - 127.1)

    cl <- exp(lcl + etalcl) * (FFM / 43.4)^e_ffm_cl * f_crcl
    vc <- exp(lvc) * (FFM / 43.4)^e_ffm_vc
    vp <- exp(lvp) * (FFM / 43.4)^e_ffm_vc
    q <- exp(lq) * (FFM / 43.4)^e_ffm_cl
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Savic transit input, as the control stream writes it:
    #   KTR = (NN+1)/MTT; PIZZA = LOG(BIO*PD*KTR) - GAMLN(NN+1)
    #   TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO) for TEMPO > 0
    # written out with podo() and tad() rather than through transit(), which
    # returns zero input once f(depot) is 0.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) +
        ntr * log(ktt) - ktt)
    }

    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # F1 = 0 in the control stream: the dose enters only through the transit
    # density.
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
