Wijk_2024_isoniazid <- function() {
  description <- paste(
    "Two-compartment population PK model for oral isoniazid at steady state",
    "in South African adults without HIV starting first-line treatment for",
    "drug-susceptible pulmonary tuberculosis (TRUST cohort, Wijk 2024).",
    "Absorption is a single transit compartment (Savic chain with NN = 0,",
    "mean transit time 0.681 h) feeding a first-order absorption compartment",
    "(ka 3.28 1/h, held constant at a literature value). Elimination is a",
    "well-stirred liver with no liver volume: the extraction ratio",
    "EH = CLint * fu / (CLint * fu + QH) both removes drug from the central",
    "compartment at QH * EH / Vc and scales the absorbed amount by the",
    "first-pass fraction FH = 1 - EH. CLint * fu is trimodal in NAT2",
    "acetylator phenotype (14.8, 36.3 and 62.4 L/h for slow, intermediate",
    "and fast) and hepatic blood flow is 90 L/h for a 56.1 kg fat-free mass.",
    "Central volume 59.3 L, peripheral volume 38.2 L, intercompartmental",
    "clearance 3.45 L/h, all allometrically scaled on fat-free mass with",
    "exponents 0.75 / 1. Random effects are between-subject variability on",
    "CLint (29.3%) and two-occasion between-occasion variability on",
    "bioavailability (22.3%), ka (130%) and mean transit time (67.5%), all",
    "omega standard deviations on the log scale. Residual error is combined",
    "proportional (17.0%) plus additive (0.021 mg/L)."
  )
  reference <- paste(
    "Wijk M, Gausi K, Malatesta S, Weber SE, Court R, Myers B, Carney T,",
    "Parry CDH, Horsburgh CR, White LF, Wiesner L, Warren RM, Uren C,",
    "McIlleron H, Kloprogge F, Denti P, Jacobson KR. The impact of alcohol",
    "and illicit substance use on the pharmacokinetics of first-line TB",
    "drugs. J Antimicrob Chemother. 2024;79(8):2022.",
    "doi:10.1093/jac/dkae206.",
    "Parameter estimates from Table 2; model equations from the NONMEM",
    "control stream reproduced in the supplementary material.",
    sep = " "
  )
  vignette <- "Wijk_2024_first_line_tb_drugs"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "isoniazid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "isoniazid",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "isoniazid",
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
        "Allometric size descriptor on CLint, Vc, Vp and Q with exponents",
        "0.75 / 1, normalised to TVFFM = 43.4 kg as in the control stream",
        "(Table 2 footnote b names the 43.1 kg cohort median); hepatic blood",
        "flow is normalised to 56.1 kg. Fat-free mass fitted better than",
        "total body weight (dOFV -27.0 versus -2.51). Cohort median 43.1 kg",
        "(IQR 35.8-47.0), Table 1."
      ),
      source_name = "FFM"
    ),
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator indicator: 1 = slow, 0 = intermediate or fast. Paired with NAT2_RAPID; both 0 denotes an intermediate acetylator.",
      units = "(binary)",
      type = "binary",
      reference_category = "No single reference level: the phenotype selects one of three independent typical values of CLint * fu.",
      notes = paste(
        "Source column NAT2 coded 0 = fast, 1 = intermediate, 2 = slow,",
        "-99 = missing (control stream $PK). Phenotype was assigned from",
        "rs1799930, rs1799931, rs1801279 and rs1801280 following Hein and",
        "Doll. Genotype was available for 34 of 104 participants (33%): 32%",
        "slow, 41% intermediate and 27% fast. For the remaining participants",
        "the authors assigned the phenotype with a $MIX mixture model whose",
        "probabilities were held at 36% slow, 40% intermediate and 24% fast,",
        "the proportions reported for a South African population (Table 2).",
        "The mixture is a device for missing genotype data and is not part of",
        "this model; to simulate a participant of unknown genotype, draw the",
        "phenotype with those probabilities."
      ),
      source_name = "NAT2"
    ),
    NAT2_RAPID = list(
      description = "NAT2 rapid (fast) acetylator indicator: 1 = fast, 0 = intermediate or slow. Paired with NAT2_SLOW; both 0 denotes an intermediate acetylator.",
      units = "(binary)",
      type = "binary",
      reference_category = "No single reference level; see NAT2_SLOW.",
      notes = paste(
        "Fast acetylators have a typical CLint * fu 1.7-fold that of",
        "intermediate and 4.2-fold that of slow acetylators. Including NAT2",
        "with the mixture model improved the fit by dOFV -49.8 (df = 2)."
      ),
      source_name = "NAT2"
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
      notes = "No statistically significant effect of alcohol use (as groups, PEth or AUDIT) or of illicit drug use on isoniazid PK (Results, Isoniazid pharmacokinetics)."
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
      "Fixed-dose combination tablets once daily (isoniazid 75 mg per",
      "tablet, with rifampicin, pyrazinamide and ethambutol) with the number",
      "of tablets by WHO weight band; sampled at steady state 4-8 weeks into",
      "treatment under fasting conditions."
    ),
    regions = "South Africa (Worcester, Western Cape)",
    notes = paste(
      "Pharmacokinetic analysis nested in the TRUST prospective cohort.",
      "Samples pre-dose and 1.5, 3, 5 and 8 h post-dose; 113 isoniazid",
      "concentrations were below the 0.105 mg/L LLOQ and handled by the M6",
      "method. NAT2 genotype known for 34 participants (11 slow, 14",
      "intermediate, 9 fast by the reported 32/41/27% split)."
    )
  )

  ini({
    # Structural values are Wijk 2024 Table 2. The supplementary control
    # stream's $THETA block carries the same values (V2 59.3349, V3 38.1712,
    # Q 3.45432, MTT 0.681154), so it is the final vector; it supplies the
    # structure. Table 2 reports CLint and fu only as their product
    # (footnote c), and the control stream's EH uses that product directly
    # (its FU = 0.95 FIX THETA(15) is never used), so the estimated
    # parameters below are the products.
    lclint_slow <- log(14.8)
    label("Unbound intrinsic clearance CLint * fu in NAT2 slow acetylators at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 CLint slow = 14.8 L/h (12.7-17.3); control stream THETA(1)
    lclint_int <- log(36.3)
    label("Unbound intrinsic clearance CLint * fu in NAT2 intermediate acetylators at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 CLint intermediate = 36.3 L/h (31.1-43.1); control stream THETA(11)
    lclint_rapid <- log(62.4)
    label("Unbound intrinsic clearance CLint * fu in NAT2 fast acetylators at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 CLint fast = 62.4 L/h (50.9-77.5); control stream THETA(16)
    lvc <- log(59.3)
    label("Central volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vc = 59.3 L (53.1-65.1); control stream THETA(2) = 59.3349
    lvp <- log(38.2)
    label("Peripheral volume of distribution at the reference fat-free mass of 43.4 kg (L)") # Table 2 Vp = 38.2 L (23.5-61.8); control stream THETA(8) = 38.1712
    lq <- log(3.45)
    label("Intercompartmental clearance at the reference fat-free mass of 43.4 kg (L/h)") # Table 2 Q = 3.45 L/h (2.09-4.92); control stream THETA(9) = 3.45432
    lka <- fixed(log(3.28))
    label("First-order absorption rate constant (1/h)") # Table 2 Ka = 3.28 1/h, held constant at the value of Gausi et al. (Results); control stream THETA(3) = 3.28 FIX
    lmtt <- log(0.681)
    label("Mean transit time of the single absorption transit compartment (h)") # Table 2 MTT = 0.681 h (0.517-0.889); control stream THETA(7) = 0.681154
    lfdepot <- fixed(log(1))
    label("Oral bioavailability before first-pass hepatic extraction (unitless)") # Table 2 F = 1, held constant; control stream THETA(4) = 1 FIX
    lqh <- fixed(log(90))
    label("Hepatic blood flow at the reference fat-free mass of 56.1 kg (L/h)") # Table 2 QH row with footnote d (90 L/h at 56.1 kg FFM); control stream THETA(14) = 90 FIX

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance terms (unitless)") # Methods 'allometric exponents for clearances and volumes fixed to 0.75 and 1'; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume terms (unitless)") # Methods; control stream ALLMV_FFM = (FFM/TVFFM)

    # Table 2 percentages are omega standard deviations on the log scale;
    # sqrt() of each control stream $OMEGA reproduces the tabulated value.
    etalclint ~ 0.0860598
    label("Between-subject variability in CLint * fu (log-scale variance)") # Table 2 BSV CL = 29.3% (19.7-37.2); control stream $OMEGA 0.0860598
    etaiov_fdepot_1 ~ 0.049752
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)") # Table 2 BOV BIO = 22.3% (21.2-23.4); control stream $OMEGA 0.049752
    etaiov_fdepot_2 ~ fix(0.049752) # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 1.68227
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)") # Table 2 BOV KA = 130% (95.9-160); control stream $OMEGA 1.68227
    etaiov_ka_2 ~ fix(1.68227) # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.455921
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)") # Table 2 BOV MTT = 67.5% (47.3-84.4); control stream $OMEGA 0.455921
    etaiov_mtt_2 ~ fix(0.455921) # control stream $OMEGA BLOCK(1) SAME

    # Control stream $ERROR: ADD = THETA(6) + 0.2 * LLOQ with THETA(6) = 0 FIX
    # and LLOQ = 0.105 mg/L.
    propSd <- 0.170
    label("Proportional residual error (fraction)") # Table 2 proportional error = 17.0% (14.7-19.4); control stream THETA(5) = 0.169928
    addSd <- fixed(0.021)
    label("Additive residual error (mg/L)") # Table 2 additive error 0.02 mg/L held constant; control stream ADD = 0 + 0.105 * 0.2 = 0.021
  })
  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # NAT2 phenotype selects one of three typical CLint * fu values (control
    # stream TVCL_NAT = THETA(1) / THETA(11) / THETA(16)).
    nat2_int <- (1 - NAT2_SLOW) * (1 - NAT2_RAPID)
    lclint_typ <- NAT2_SLOW * lclint_slow + nat2_int * lclint_int +
      NAT2_RAPID * lclint_rapid

    clint <- exp(lclint_typ + etalclint) * (FFM / 43.4)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 43.4)^e_ffm_vc
    vp <- exp(lvp) * (FFM / 43.4)^e_ffm_vc
    q <- exp(lq) * (FFM / 43.4)^e_ffm_cl
    qh <- exp(lqh) * (FFM / 56.1)^e_ffm_cl
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    fdepot <- exp(lfdepot + iov_fdepot)

    # Well-stirred liver without a liver compartment (control stream $PK):
    # EH = CLINT/(CLINT+QH); FH = 1 - EH; K20 = QH*EH/V2.
    eh <- clint / (clint + qh)
    fh <- 1 - eh
    kel <- qh * eh / vc
    k12 <- q / vc
    k21 <- q / vp

    # Savic transit input with NN = 0 (control stream THETA(10) = 0 FIX), so
    # KTR = 1/MTT and the density reduces to BIO*PD*KTR*EXP(-KTR*TEMPO): one
    # transit compartment ahead of the absorption compartment. PD and TEMPO
    # are the latest dose and the time since it (podo() and tad()).
    ktr <- 1 / mtt
    tdos <- tad(depot)
    trin <- 0
    if (tdos > 0) {
      trin <- fdepot * podo(depot) * ktr * exp(-ktr * tdos)
    }

    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot * fh - k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # F1 = 0 in the control stream: the dose enters only through the transit
    # density.
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
