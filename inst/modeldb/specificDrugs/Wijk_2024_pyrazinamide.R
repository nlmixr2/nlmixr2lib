Wijk_2024_pyrazinamide <- function() {
  description <- paste(
    "One-compartment population PK model for oral pyrazinamide at steady",
    "state in South African adults without HIV starting first-line treatment",
    "for drug-susceptible pulmonary tuberculosis (TRUST cohort, Wijk 2024).",
    "Absorption is a Savic analytical transit chain (NN = 2.73, i.e. 3.73",
    "transit compartments as Table 2 counts them, mean transit time 0.289 h)",
    "feeding a first-order absorption compartment (ka 1.89 1/h), with",
    "first-order elimination (CL 4.48 L/h, V 41.4 L) allometrically scaled on",
    "fat-free mass with exponents 0.75 / 1. Random effects are",
    "between-subject variability on CL (20.3%) and two-occasion",
    "between-occasion variability on bioavailability (10.8%), ka (69.1%) and",
    "mean transit time (90.2%), all omega standard deviations on the log",
    "scale; the between-occasion etas are 2.49-fold larger for a dose whose",
    "intake was not observed by study staff. Residual error is combined",
    "proportional (2.66%) plus additive (0.63 mg/L)."
  )
  reference <- paste(
    "Wijk M, Gausi K, Malatesta S, Weber SE, Court R, Myers B, Carney T,",
    "Parry CDH, Horsburgh CR, White LF, Wiesner L, Warren RM, Uren C,",
    "McIlleron H, Kloprogge F, Denti P, Jacobson KR. The impact of alcohol",
    "and illicit substance use on the pharmacokinetics of first-line TB",
    "drugs. J Antimicrob Chemother. 2024;79(8):2022.",
    "doi:10.1093/jac/dkae206.",
    "Parameter estimates from Table 2; model equations and the transit",
    "compartment count from the NONMEM control stream reproduced in the",
    "supplementary material.",
    sep = " "
  )
  vignette <- "Wijk_2024_first_line_tb_drugs"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "pyrazinamide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "pyrazinamide",
      units = "mg",
      specimen = "plasma",
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
        "Allometric size descriptor on CL and V with exponents 0.75 / 1,",
        "normalised to TVFFM = 43.4 kg as in the control stream (Table 2",
        "footnote b names the 43.1 kg cohort median). Fat-free mass fitted",
        "better than total body weight (dOFV -107 versus -35.6). Cohort",
        "median 43.1 kg (IQR 35.8-47.0), Table 1."
      ),
      source_name = "FFM"
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
    ),
    SELFADMIN = list(
      description = "Self-administered (not directly observed) dosing-occasion indicator (1 = dose intake not observed by study staff, 0 = observed dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (observed dose)",
      notes = paste(
        "Carries the control stream's SPARSE column: IF (SPARSE==1) the",
        "occasion etas on bioavailability, ka and mean transit time are",
        "multiplied by THETA(9) = 2.49. The Results describe it as a",
        "2.49-fold larger BOV in the absorption parameters for the pre-dose",
        "concentrations, 'i.e. the concentrations resulting from doses whose",
        "intake was not observed by study staff' (dOFV -58.7, df = 1). In the",
        "study this is occasion 1. Set SELFADMIN = 0 to simulate observed",
        "dosing. Source column SPARSE."
      ),
      source_name = "SPARSE"
    )
  )

  covariatesDataExcluded <- list(
    ALCOHOL_USE_HIGH = list(
      description = "High alcohol use group from a latent class analysis of AUDIT, Timeline Follow-Back and phosphatidylethanol.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low or moderate alcohol use)",
      notes = "No statistically significant effect of alcohol use or illicit drug use on pyrazinamide PK (Results, Pyrazinamide pharmacokinetics)."
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
    dose_range = paste(
      "Fixed-dose combination tablets once daily (pyrazinamide 400 mg per",
      "tablet, with rifampicin, isoniazid and ethambutol) with the number of",
      "tablets by WHO weight band; sampled at steady state 4-8 weeks into",
      "treatment under fasting conditions."
    ),
    regions = "South Africa (Worcester, Western Cape)",
    notes = paste(
      "Pharmacokinetic analysis nested in the TRUST prospective cohort.",
      "Samples pre-dose and 1.5, 3, 5 and 8 h post-dose; one pyrazinamide",
      "concentration was below the 0.2 mg/L LLOQ."
    )
  )

  ini({
    # Structural values are Wijk 2024 Table 2. The supplementary control
    # stream's $THETA block carries the same final values to more digits
    # (CL 4.47911, V 41.3614, KA 1.88803, MTT 0.288724, PROP 0.0266402,
    # ADD 0.589092, BOV scaling 2.49166) and supplies the structure.
    lcl <- log(4.48)
    label("Clearance at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 CL = 4.48 L/h (4.27-4.70); control stream THETA(1) = 4.47911
    lvc <- log(41.4)
    label("Central volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vc = 41.4 L (40.2-42.6); control stream THETA(2) = 41.3614
    lka <- log(1.89)
    label("First-order absorption rate constant (1/h)") # Table 2 Ka = 1.89 1/h (1.57-2.36); control stream THETA(3) = 1.88803
    lmtt <- log(0.289)
    label("Mean transit time through the absorption transit chain (h)") # Table 2 MTT = 0.289 h (0.180-0.428); control stream THETA(7) = 0.288724
    # Table 2 prints NN = 3.73 (2.67-5.40). The control stream estimate that
    # enters KTR = (NN+1)/MTT and the gamma density is THETA(8) = 2.73017:
    # every other THETA in that block reproduces its Table 2 row, and Table 2
    # reports NN + 1 here (the same offset holds for ethambutol, 3.01 against
    # 2.01474, and isoniazid, 1 against 0). The control stream value is used.
    lntr <- log(2.73017)
    label("Shape parameter NN of the Savic absorption transit chain (unitless)") # control stream THETA(8) = 2.73017, KTR = (NN+1)/MTT; Table 2 NN = 3.73 = NN + 1
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (unitless)") # Table 2 F = 1, held constant; control stream THETA(4) = 1 FIX
    e_selfadmin_iov <- 2.49
    label("Multiplier on the between-occasion etas of F, ka and MTT for a dose not observed by study staff (unitless)") # Table 2 BOV scaling for pre-dose = 2.49-fold (1.84-3.04); control stream THETA(9) = 2.49166

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance (unitless)") # Methods 'allometric exponents for clearances and volumes fixed to 0.75 and 1'; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume (unitless)") # Methods; control stream ALLMV_FFM = (FFM/TVFFM)

    # Table 2 percentages are omega standard deviations on the log scale;
    # sqrt() of each control stream $OMEGA reproduces the tabulated value.
    etalcl ~ 0.0413833
    label("Between-subject variability in CL (log-scale variance)") # Table 2 BSV CL = 20.3% (17.4-24.6); control stream $OMEGA 0.0413833
    etaiov_fdepot_1 ~ 0.011694
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)") # Table 2 BOV BIO = 10.8% (9.01-13.2); control stream $OMEGA 0.011694
    etaiov_fdepot_2 ~ fix(0.011694) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 0.477405
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)") # Table 2 BOV KA = 69.1% (59.1-81.3); control stream $OMEGA 0.477405
    etaiov_ka_2 ~ fix(0.477405) # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.813801
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)") # Table 2 BOV MTT = 90.2% (66.0-115); control stream $OMEGA 0.813801
    etaiov_mtt_2 ~ fix(0.813801) # control stream $OMEGA BLOCK(1) SAME

    # Control stream $ERROR: ADD = THETA(6) + 0.2 * LLOQ with THETA(6)
    # estimated (0.589092, Table 2 0.59) and LLOQ = 0.2 mg/L, so the total
    # additive standard deviation is 0.59 + 0.04 = 0.63 mg/L.
    propSd <- 0.0266
    label("Proportional residual error (fraction)") # Table 2 proportional error = 2.66% (1.33-3.72); control stream THETA(5) = 0.0266402
    addSd <- 0.63
    label("Additive residual error (mg/L)") # Table 2 additive error 0.59 mg/L (0.38-0.82) plus the control stream's 0.2 * LLOQ = 0.04 mg/L term
  })
  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    # Etas of a dose not observed by study staff are inflated 2.49-fold
    # (control stream IF (SPARSE==1) BOVBIO = TVBOVBIO * THETA(9), and the
    # same for ka and MTT).
    iov_scale <- 1 + (e_selfadmin_iov - 1) * SELFADMIN
    iov_fdepot <- (oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2) * iov_scale
    iov_ka <- (oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2) * iov_scale
    iov_mtt <- (oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2) * iov_scale

    cl <- exp(lcl + etalcl) * (FFM / 43.4)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 43.4)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    kel <- cl / vc

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
    d/dt(central) <- ka * depot - kel * central

    # F1 = 0 in the control stream: the dose enters only through the transit
    # density.
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
