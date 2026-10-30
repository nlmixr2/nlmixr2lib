Wijk_2024_rifampicin <- function() {
  description <- paste(
    "One-compartment population PK model for oral rifampicin at steady state",
    "in South African adults without HIV starting first-line treatment for",
    "drug-susceptible pulmonary tuberculosis (TRUST cohort, Wijk 2024).",
    "Absorption is a Savic analytical transit chain (17 transit compartments,",
    "mean transit time 0.689 h) feeding a first-order absorption compartment",
    "(ka 1.82 1/h) that empties into a liver compartment, so first-pass",
    "extraction is structural. Elimination is a well-stirred liver with",
    "saturable intrinsic clearance, CLint = CLint,max * Km / (CH + Km) with",
    "Km 7 mg/L and CLint,max * fu 42.5 L/h, liver volume 1 L and hepatic",
    "blood flow 90 L/h for a 56.1 kg fat-free mass (both held constant).",
    "Central volume 55.1 L. Disposition parameters are allometrically scaled",
    "on fat-free mass with exponents 0.75 / 1. Bioavailability is 1 in",
    "participants with a positive urine test for methaqualone and 19.2%",
    "higher in those without. Random effects are between-subject variability",
    "on CLint,max (35.2%) and two-occasion between-occasion variability on",
    "bioavailability (5.40%), ka (136%) and mean transit time (79.9%), all",
    "omega standard deviations on the log scale. Residual error is combined",
    "proportional (13.7%) plus additive (0.0234 mg/L)."
  )
  reference <- paste(
    "Wijk M, Gausi K, Malatesta S, Weber SE, Court R, Myers B, Carney T,",
    "Parry CDH, Horsburgh CR, White LF, Wiesner L, Warren RM, Uren C,",
    "McIlleron H, Kloprogge F, Denti P, Jacobson KR. The impact of alcohol",
    "and illicit substance use on the pharmacokinetics of first-line TB",
    "drugs. J Antimicrob Chemother. 2024;79(8):2022.",
    "doi:10.1093/jac/dkae206.",
    "Parameter estimates from Table 2; model equations and the methaqualone",
    "effect from the NONMEM control stream reproduced in the supplementary",
    "material. Fat-free mass follows Janmahasatian et al. (2005) Clin",
    "Pharmacokinet 44(10):1051-1065 doi:10.2165/00003088-200544100-00004.",
    sep = " "
  )
  vignette <- "Wijk_2024_first_line_tb_drugs"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    liver = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central = list(
      analyte = "rifampicin",
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
        "Allometric size descriptor on every disposition parameter, with",
        "exponents 0.75 for clearances and 1 for volumes (Methods,",
        "Pharmacokinetic modelling analysis). Two references are used, as in",
        "the control stream $PK: CLint,max and V are normalised to",
        "TVFFM = 43.4 kg, and the hepatic blood flow and liver volume to",
        "56.1 kg, the fat-free mass of the standard 70 kg, 176 cm male for",
        "whom QH = 90 L/h and VH = 1 L (Table 2 footnote d). Table 2",
        "footnote b names 43.1 kg (the Table 1 cohort median) as the",
        "allometric reference; the control stream that produced the",
        "estimates carries 43.4 kg, which is used here. The control stream",
        "computes FFM = WHSMAX * HT^2 * WT / (WHS50 * HT^2 + WT) with HT in",
        "m, WHSMAX 42.92 and WHS50 30.93 for males and 37.99 and 35.98 for",
        "females. Fat-free mass fitted better than total body weight (dOFV",
        "-63.0 versus -12.0). Cohort median 43.1 kg (IQR 35.8-47.0), Table 1."
      ),
      source_name = "FFM"
    ),
    CONMED_METHAQUALONE = list(
      description = "1 = urine drug test positive for methaqualone (recreational use); 0 = negative",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (methaqualone positive), the reference level of the source control stream",
      notes = paste(
        "Control stream $PK: BIO_METHAQUALONE = 1 (comment 'reference,",
        "methaqualone positive'); IF (METHAQUALONE.EQ.0) BIO_METHAQUALONE =",
        "(1 + THETA(13)), THETA(13) = 0.191955. So the typical",
        "bioavailability of 1 belongs to methaqualone users and non-users",
        "have 19.2% higher bioavailability; the Results phrase this as 19%",
        "lower bioavailability with methaqualone use (dOFV -34.3, P < 0.001).",
        "Methaqualone (Mandrax) use was reported by 29 of 104 participants",
        "(28%, Table 1), ascertained by a urine drug test. Source column",
        "METHAQUALONE."
      ),
      source_name = "METHAQUALONE"
    ),
    OCC = list(
      description = "Occasion index (1 or 2) for the between-occasion variability on bioavailability, ka and mean transit time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each dose and its following samples are one occasion (Methods). The",
        "control stream multiplexes the etas with IF (OCC==1) / IF (OCC==2)",
        "over two occasions, each pair sharing one variance ($OMEGA BLOCK(1)",
        "followed by SAME). In the study, occasion 1 is the dose before the",
        "PK visit with its pre-dose sample and occasion 2 is the observed",
        "dose on the PK visit with the 1.5, 3, 5 and 8 h samples. Any other",
        "OCC value switches the between-occasion variability off."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the allometric descriptor and rejected in favour of fat-free mass; still needed upstream as an input to the Janmahasatian FFM formula. Cohort median 51.4 kg (IQR 45.5-56.9), Table 1."
    ),
    ALCOHOL_USE_HIGH = list(
      description = "High alcohol use group from a latent class analysis of AUDIT, Timeline Follow-Back and phosphatidylethanol.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low or moderate alcohol use)",
      notes = "High alcohol use gave 13% higher bioavailability (dOFV -6.11, P = 0.01), falling to 9.2% and losing significance (dOFV -1.44, P = 0.23) once methaqualone use was in the model, so it was not retained (Results, Rifampicin pharmacokinetics). PEth and AUDIT were also tested on their own. 29 of 104 participants (28%) were in the high-use group."
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
    disease_state = "Drug-susceptible pulmonary tuberculosis, without HIV, starting first-line treatment; 42% moderate and 28% high alcohol use, 28% methaqualone use",
    dose_range = paste(
      "Fixed-dose combination tablets once daily (rifampicin 150 mg,",
      "isoniazid 75 mg, pyrazinamide 400 mg, ethambutol 275 mg per tablet)",
      "with the number of tablets by WHO weight band; sampled at steady",
      "state 4-8 weeks into treatment under fasting conditions."
    ),
    regions = "South Africa (Worcester, Western Cape)",
    notes = paste(
      "Pharmacokinetic analysis nested in the TRUST prospective cohort",
      "(May 2017 to April 2022). Samples pre-dose and 1.5, 3, 5 and 8 h",
      "post-dose; 511 drug concentrations in total, of which 104 rifampicin",
      "concentrations were below the 0.117 mg/L LLOQ and handled by the M6",
      "method. Baseline characteristics are Table 1."
    )
  )

  ini({
    # Structural values are Wijk 2024 Table 2 (95% CI by sampling importance
    # resampling). The supplementary control stream's $THETA block agrees to
    # within 0.3% (exp(3.7471) = 42.4 for CLint,max * fu, exp(4.0075) = 55.0
    # for V), so Table 2 is used for the estimated values and the control
    # stream for structure and for the methaqualone effect, which Table 2
    # does not list.
    lclint_max <- log(42.5)
    label("Maximal unbound intrinsic hepatic clearance CLint,max * fu at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 CLint,max * fu = 42.5 L/h (38.4-47.1); control stream THETA(2) = 3.7471 on the log scale
    lkm <- fixed(log(7))
    label("Michaelis-Menten constant Km for saturable hepatic intrinsic clearance (mg/L)") # Table 2 Km = 7 mg/L held constant at the value the model estimated; control stream THETA(1) = 1.94591014906 FIX = log(7)
    lvc <- log(55.1)
    label("Central volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vc = 55.1 L (51.6-59.0); control stream THETA(3) = 4.0075 on the log scale
    lka <- log(1.82)
    label("First-order absorption rate constant from the absorption compartment into the liver (1/h)") # Table 2 Ka = 1.82 1/h (1.15-2.87); control stream THETA(4) = 0.597532 on the log scale
    lmtt <- log(0.689)
    label("Mean transit time through the absorption transit chain (h)") # Table 2 MTT = 0.689 h (0.566-0.838); control stream THETA(8) = -0.373117 on the log scale
    lntr <- fixed(log(17))
    label("Number of absorption transit compartments in the Savic chain (unitless)") # Table 2 NN = 17, held constant; control stream THETA(9) = 2.83321334406 FIX = log(17)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability in a methaqualone-positive participant (unitless)") # Table 2 F = 1, held constant; control stream THETA(5) = 0 FIX on the log scale
    e_no_methaqualone_fdepot <- 0.191955
    label("Fractional increase in bioavailability in participants without methaqualone use (fraction)") # control stream THETA(13) = 0.191955, BIO_METHAQUALONE = 1 + THETA(13) when METHAQUALONE = 0; Results 'Participants with reported methaqualone use had 19% lower bioavailability'

    # Hepatic physiology held constant (Table 2 footnote d: equivalent to
    # QH 90 L/h and VH 1 L for a standard 70 kg, 176 cm male of 56.1 kg FFM;
    # the table rows 73.9 L/h and 0.77 L are these scaled to the typical
    # fat-free mass). The control stream also carries FU = 0.2 FIX as
    # THETA(11) but never uses it: Table 2 footnote c reports CLint and fu
    # only as their product, which is the estimated parameter above.
    lqh <- fixed(log(90))
    label("Hepatic blood flow at the reference fat-free mass of 56.1 kg (L/h)") # Table 2 QH 73.9 L/h held constant with footnote d (90 L/h at 56.1 kg FFM); control stream THETA(10) = 90 FIX
    lvh <- fixed(log(1))
    label("Liver compartment volume at the reference fat-free mass of 56.1 kg (L)") # Table 2 VH 0.77 L held constant with footnote d (1 L at 56.1 kg FFM); control stream THETA(12) = 1 FIX

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance terms (unitless)") # Methods 'allometric exponents for clearances and volumes fixed to 0.75 and 1'; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume terms (unitless)") # Methods; control stream ALLMV_FFM = (FFM/TVFFM)

    # Table 2 percentages are omega standard deviations on the log scale:
    # sqrt(0.123619) = 0.352 reproduces the 35.2% BSV and sqrt(1.8383) = 1.356
    # the 136% BOV on ka, which the log-normal CV formula would not. The
    # variances are the control stream's $OMEGA values. The control stream
    # also carries zero-variance etas on V, ka, F, MTT, NN and Km, which are
    # omitted here.
    etalclint_max ~ 0.123619
    label("Between-subject variability in CLint,max (log-scale variance)") # Table 2 BSV CL = 35.2% (29.6-40.1); control stream $OMEGA 0.123619
    etaiov_fdepot_1 ~ 0.00292323
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)") # Table 2 BOV BIO = 5.40% (4.73-6.00); control stream $OMEGA 0.00292323
    etaiov_fdepot_2 ~ fix(0.00292323) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 1.8383
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)") # Table 2 BOV KA = 136% (91.2-169); control stream $OMEGA 1.8383
    etaiov_ka_2 ~ fix(1.8383) # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.638919
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)") # Table 2 BOV MTT = 79.9% (53.7-99.5); control stream $OMEGA 0.638919
    etaiov_mtt_2 ~ fix(0.638919) # control stream $OMEGA BLOCK(1) SAME

    # Control stream $ERROR: PROP = IPRED * THETA(6); ADD = THETA(7) + 0.2 * LLOQ
    # with THETA(7) = 0 FIX and LLOQ = 0.117 mg/L; W = SQRT(ADD**2 + PROP**2),
    # which is nlmixr2's combined add() + prop().
    propSd <- 0.137
    label("Proportional residual error (fraction)") # Table 2 proportional error = 13.7% (11.7-15.7); control stream THETA(6) = 0.136958
    addSd <- fixed(0.0234)
    label("Additive residual error (mg/L)") # Table 2 additive error 0.02 mg/L held constant; control stream ADD = 0 + 0.117 * 0.2 = 0.0234
  })
  model({
    # Occasion indicators for the between-occasion variability (control
    # stream IF (OCC==1) / IF (OCC==2) multiplexers).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # Allometric scaling on fat-free mass: 43.4 kg for the estimated
    # disposition parameters, 56.1 kg for the fixed hepatic physiology.
    clint_max <- exp(lclint_max + etalclint_max) * (FFM / 43.4)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 43.4)^e_ffm_vc
    qh <- exp(lqh) * (FFM / 56.1)^e_ffm_cl
    vh <- exp(lvh) * (FFM / 56.1)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    km <- exp(lkm)

    # Methaqualone-positive participants are the reference (F = 1).
    f_methaqualone <- 1 + e_no_methaqualone_fdepot * (1 - CONMED_METHAQUALONE)
    fdepot <- exp(lfdepot + iov_fdepot) * f_methaqualone

    # Savic transit input, as the control stream writes it:
    #   KTR = (NN+1)/MTT; PIZZA = LOG(BIO*PD*KTR) - GAMLN(NN+1)
    #   TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO) for TEMPO > 0
    # with PD the latest dose and TEMPO the time since it. Written out with
    # podo() and tad() rather than through transit(), which returns zero input
    # once f(depot) is 0.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) +
        ntr * log(ktt) - ktt)
    }

    # Well-stirred liver with saturable intrinsic clearance (control stream
    # $DES): VMAX = CLINT_MAX * KM; SAT_CL_FU = VMAX / (CH + KM);
    # EH = SAT_CL_FU / (SAT_CL_FU + QH); K20 = QH*EH/VH; K23 = QH*FH/VH;
    # K32 = QH/V. clint_max is already the unbound product CLint,max * fu.
    vmax <- clint_max * km
    c_liver <- liver / vh
    clint <- vmax / (c_liver + km)
    eh <- clint / (clint + qh)
    fh <- 1 - eh
    kliver_el <- qh * eh / vh
    kliver_c <- qh * fh / vh
    kc_liver <- qh / vc

    d/dt(depot) <- trin - ka * depot
    d/dt(liver) <- ka * depot - kliver_el * liver - kliver_c * liver + kc_liver * central
    d/dt(central) <- kliver_c * liver - kc_liver * central

    # F1 = 0 in the control stream: the whole dose enters through the transit
    # density, never as a bolus into the absorption compartment.
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
