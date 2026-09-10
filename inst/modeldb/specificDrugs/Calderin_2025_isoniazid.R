Calderin_2025_isoniazid <- function() {
  description <- paste(
    "Two-compartment population PK model for oral isoniazid in plasma and",
    "lumbar cerebrospinal fluid (CSF) in South African adults with",
    "HIV-associated tuberculous meningitis (LASER-TBM PK substudy, Calderin",
    "2025). Absorption is a Savic transit chain (mean transit time 0.249 h,",
    "chain length fixed at 5) feeding a first-order depot (ka 2.21 1/h).",
    "Elimination is semi-mechanistic: a well-stirred liver model in which",
    "the estimated parameter is intrinsic clearance, so the hepatic",
    "extraction ratio EH = fu * CLint / (fu * CLint + Qh) both removes drug",
    "from the central compartment (rate Qh * EH / Vc) and reduces the",
    "absorbed fraction by first pass (FH = 1 - EH). Hepatic plasma flow is",
    "fixed at 90 L/h for a 56.1 kg fat-free mass (76.3 L/h at the cohort",
    "median) and the unbound fraction at 95%. Intrinsic clearance is",
    "trimodal in NAT2 acetylator phenotype: 14.6, 32.2 and 64.7 L/h for slow,",
    "intermediate and rapid acetylators at a fat-free mass of 45 kg. All",
    "disposition parameters are allometrically scaled on fat-free mass with",
    "fixed 0.75 / 1 exponents. CSF is a Sheiner-style effect compartment",
    "holding a concentration, equilibrating with plasma at 0.179 1/h",
    "(equilibration half-life 3.87 h) toward a CSF-to-plasma pseudo-partition",
    "coefficient of 1.04, i.e. CSF exposure matches plasma. Random effects",
    "are between-subject variability on intrinsic clearance (25.2%) and",
    "five-occasion between-occasion variability on bioavailability (32.1%),",
    "absorption rate constant (87.0%) and mean transit time (139%); all",
    "reported percentages are the omega standard deviation on the log scale.",
    "Residual error is combined proportional plus additive, separately for",
    "plasma (16.4%, 0.021 mg/L) and CSF (58.8%, 0.0117 mg/L). Neither",
    "high-dose rifampicin (35 mg/kg) nor study visit affected isoniazid PK."
  )
  reference <- paste(
    "Calderin JM, Wasserman S, Resendiz-Galvan JE, Abdelgawad N, Davis A,",
    "Stek C, Wiesner L, Wilkinson RJ, Denti P (2025).",
    "Population pharmacokinetics of pyrazinamide and isoniazid in plasma and",
    "cerebrospinal fluid from South African adults with tuberculous",
    "meningitis. Antimicrob Agents Chemother 69(8):e00099-25.",
    "doi:10.1128/aac.00099-25",
    sep = " "
  )
  vignette <- "Calderin_2025_pyrazinamide_isoniazid_tbm"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state is an exception to the usual "states hold an amount"
  # rule. Supplementary S2 writes the effect compartment directly in
  # concentration units, dC_CSF/dt = k_Plasma-CSF * (PPC * C_Plasma - C_CSF),
  # and the S10 control stream integrates that same equation as
  # DADT(4) = KE0*(PPC*C2 - A(4)) with the CSF observation read as
  # CE = A(4) (not A(4)/V). So `csf` carries mg/L, not mg.
  compartmentData <- list(
    depot = list(
      analyte = "isoniazid", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "isoniazid", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "isoniazid", units = "mg",
      specimen = "tissue", verified = TRUE
    ),
    csf = list(
      analyte = "isoniazid", units = "mg/L",
      specimen = "CSF", verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = paste(
        "Fat-free mass, computed with the Janmahasatian et al. formula from",
        "total body weight, height and sex. Allometric scaling on FFM was",
        "retained over total body weight (dOFV -16.95 vs -9.21)."
      ),
      units = "kg",
      type = "continuous",
      source_name = "FFM",
      reference_value = 45,
      notes = paste(
        "Cohort median FFM 45 kg (range 30-59), Table 1. The S10 control",
        "stream carries two allometric references: ALLMCL_FFM = (FFM/45)^0.75",
        "and ALLMV_FFM = (FFM/45) for the disposition parameters, and",
        "ALLMCL_FFM_HEP = (FFM/56.1)^0.75 for hepatic plasma flow, whose",
        "fixed 90 L/h value is defined for a 70 kg male (fat-free mass 56.1",
        "kg). Missing heights needed for the Janmahasatian formula were",
        "imputed with the Johansson and Karlsson regression given in",
        "supplementary S3."
      )
    ),
    NAT2_SLOW = list(
      description = paste(
        "NAT2 slow-acetylator indicator: 1 = slow acetylator, 0 =",
        "intermediate or rapid. Paired with NAT2_RAPID; the joint state",
        "(both 0) denotes an intermediate acetylator."
      ),
      units = "(binary)",
      type = "binary",
      source_name = "NAT2_phenotype",
      reference_category = paste(
        "The source model has no single reference level: it selects one of",
        "three independent typical-value intrinsic clearances from the",
        "phenotype (S10 IF (NAT2_phenotype.EQ.0/1/2) THEN TVCL = THETA(1/2/3))."
      ),
      notes = paste(
        "Phenotype assigned from rs1801279 (NAT2*14), rs1801280 (NAT2*5),",
        "rs1799930 (NAT2*6) and rs1799931 (NAT2*7). Genotype was available",
        "for 63% (n = 31) of participants: 19% slow, 55% intermediate, 26%",
        "rapid. The remaining participants were assigned a phenotype with a",
        "mixture model whose class probabilities were fixed to the observed",
        "frequencies; the mixture assignment is a fitting device and is not",
        "reproduced here, so the phenotype is supplied as a covariate."
      )
    ),
    NAT2_RAPID = list(
      description = paste(
        "NAT2 rapid (fast) acetylator indicator: 1 = rapid acetylator, 0 =",
        "intermediate or slow. Paired with NAT2_SLOW; the joint state",
        "(both 0) denotes an intermediate acetylator."
      ),
      units = "(binary)",
      type = "binary",
      source_name = "NAT2_phenotype",
      reference_category = paste(
        "The source model has no single reference level; see NAT2_SLOW."
      ),
      notes = paste(
        "Rapid acetylators had a typical intrinsic clearance 2.0-fold that",
        "of intermediate and 4.4-fold that of slow acetylators, the largest",
        "single source of isoniazid exposure variability in the cohort."
      )
    ),
    OCC = list(
      description = paste(
        "Occasion index for between-occasion variability. An occasion is a",
        "dosing event and its subsequent observations."
      ),
      units = "(count)",
      type = "categorical",
      source_name = "OCC",
      notes = paste(
        "The S10 control stream multiplexes the BOV etas with IF (OCC==k)",
        "blocks over five occasions (ETA 15-19 bioavailability, ETA 20-24",
        "ka, ETA 25-29 MTT), each declared $OMEGA BLOCK(1) SAME after the",
        "first. Unlike the companion pyrazinamide model, no study-visit",
        "effect was retained for isoniazid."
      )
    )
  )

  population <- list(
    n_subjects = 49,
    n_studies = 1,
    species = "human",
    age_range = "39 years (range 25-78)",
    weight_range = "60.0 kg (range 30.0-107)",
    sex_female_pct = 45,
    race_ethnicity = "South African; not further reported",
    disease_state = paste(
      "HIV-associated tuberculous meningitis; all participants received",
      "adjunctive dexamethasone"
    ),
    dose_range = paste(
      "Isoniazid 5 mg/kg once daily by WHO weight band, given as an oral",
      "fixed-dose combination"
    ),
    regions = "South Africa (Cape Town and Gqeberha)",
    notes = paste(
      "PK substudy nested in the open-label randomised phase 2A LASER-TBM",
      "trial. 414 plasma and 44 CSF isoniazid concentrations from 49",
      "participants at the day-3 visit and 34 at the day-28 visit. Arms:",
      "standard of care (rifampicin 10 mg/kg) versus high-dose rifampicin",
      "(35 mg/kg) plus linezolid, with or without aspirin; neither the",
      "rifampicin dose nor aspirin affected isoniazid PK. Intensive sampling",
      "on day 3 (pre-dose, 0.5, 1, 2, 3, 6, 8-10, 24 h) and sparse sampling",
      "on day 28 (pre-dose, 2, 4 h), with one lumbar CSF sample per visit",
      "randomised to a 1-3, 3-6, 6-10 or 24 h window. Plasma LLOQ 0.105 mg/L",
      "(27.5% BLQ), CSF LLOQ 0.0586 mg/L (2.3% BLQ). Median fat-free mass 45",
      "kg (range 30-59). NAT2 genotype available for 31 participants: 6 slow,",
      "17 intermediate, 8 rapid."
    )
  )

  ini({
    # --- Intrinsic clearance by NAT2 acetylator phenotype. The S10 control
    # stream comments that "CL now represents intrinsic clearance rather
    # than oral clearance", so Table 2's three clearance rows are CLint, not
    # CL/F. The well-stirred identity CL/F = fu * CLint reproduces the more
    # familiar oral clearances: 13.9, 30.6 and 61.5 L/h.
    lclint_slow <- log(14.6)
    label("Log of intrinsic clearance in NAT2 slow acetylators at FFM 45 kg (L/h)")               # Table 2 Slow acetylator 14.6 (12.1-17.1)
    lclint_int <- log(32.2)
    label("Log of intrinsic clearance in NAT2 intermediate acetylators at FFM 45 kg (L/h)")       # Table 2 Intermediate acetylator 32.2 (28.6-37.0)
    lclint_rapid <- log(64.7)
    label("Log of intrinsic clearance in NAT2 rapid acetylators at FFM 45 kg (L/h)")              # Table 2 Rapid acetylator 64.7 (54.1-77.7)

    # --- Distribution.
    lvc <- log(43.6)
    label("Log of central volume of distribution at FFM 45 kg (L)")                               # Table 2 Central volume of distribution 43.6 (39.7-48.7)
    lvp <- log(22.3)
    label("Log of peripheral volume of distribution at FFM 45 kg (L)")                            # Table 2 Peripheral volume of distribution 22.3 (15.5-30.3)
    lq <- log(5.02)
    label("Log of intercompartmental clearance at FFM 45 kg (L/h)")                               # Table 2 Intercompartmental clearance 5.02 (3.41-6.56)

    # Allometric exponents fixed a priori, not estimated: "the exponents for
    # clearance and volume were fixed to 0.75 and 1, respectively".
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on fat-free mass for clearance-type parameters")             # Methods, Pharmacokinetic modeling; S10 ALLMCL_FFM = (FFM/45)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent on fat-free mass for volume-type parameters")                # Methods, Pharmacokinetic modeling; S10 ALLMV_FFM = (FFM/45)

    # --- Well-stirred liver. Both values are fixed from the literature, not
    # estimated. Qh is defined at a 56.1 kg fat-free mass, which is the
    # fat-free mass of the 70 kg male the 90 L/h value refers to; Table 2
    # reports the value rescaled to the cohort median, 90*(45/56.1)^0.75 =
    # 76.3 L/h, so the control stream's 90 L/h at 56.1 kg is the parameter
    # and 76.3 L/h is its typical-individual consequence.
    lqh <- fixed(log(90))
    label("Log of hepatic plasma flow at FFM 56.1 kg (L/h)")                               # Methods "typical value of hepatic blood flow (Qh) was fixed at 90 L/h, which corresponds to a 70 kg male"; S10 $THETA 11 (90) FIX with ALLMCL_FFM_HEP = (FFM/56.1)**0.75; Table 2 footnote e
    fu <- fixed(0.95)
    label("Unbound fraction of isoniazid in plasma")                                      # Methods "the unbound fraction of isoniazid was fixed at 95%"; S10 $THETA 12 (0.95) FIX

    # --- Absorption. Savic transit chain feeding a first-order depot.
    lka <- log(2.21)
    label("Log of the first-order absorption rate constant (1/h)")                                # Table 2 First-order absorption rate constant 2.21 (1.59-3.08)
    lmtt <- log(0.249)
    label("Log of the mean absorption transit time (h)")                                          # Table 2 Mean absorption transit time 0.249 (0.162-0.328)
    lntr <- fixed(log(5))
    label("Log of the number of absorption transit compartments (unitless, 5)")                    # Table 2 Number of absorption transit compartments 5 Fixed; footnote d "fixed at 5, based on the previously estimated value, to improve model stability"
    lfdepot <- fixed(log(1))
    label("Log of pre-hepatic bioavailability (unitless, 1; all other parameters are relative to it)")  # Table 2 Bioavailability 1 Fixed; Methods "estimated relative to the pre-hepatic bioavailability, with the typical value of this parameter fixed to 1"; S10 $THETA 6 (1) FIX

    # --- CSF effect compartment. Table 2 reports the equilibration
    # half-life, so ke0 is recovered as log(2) / HL rather than read off the
    # control stream, whose $THETA block holds initial estimates only.
    lke0 <- log(log(2) / 3.87)
    label("Log of the plasma-to-CSF equilibration rate constant (1/h)")                           # Table 2 Plasma-to-CSF equilibrium half-life 3.87 h (2.47-7.59) -> ke0 = log(2)/3.87 = 0.179 1/h
    lppc <- log(1.04)
    label("Log of the CSF-to-plasma pseudo-partition coefficient (unitless)")                     # Table 2 CSF-to-plasma pseudo-partition coefficient 1.04 (0.76-1.38)

    # --- Random effects. As for the companion pyrazinamide model, Table 2
    # reports each variance component as the omega standard deviation on the
    # log scale times 100 (confirmed there by the S9 $OMEGA initials
    # reproducing the tabulated percentages exactly as sqrt(omega^2)).
    etalclint ~ 0.0635
    label("Between-subject variability in intrinsic clearance (log-scale variance)")              # Table 2 BSV in clearance 25.2% (17.0-30.5); 0.252^2 = 0.0635. Adding NAT2 dropped this from 47.7% to 25.2%
    # The control stream fixes BSV on V, ka, bioavailability, Vp, Q, MTT, ke0
    # and PPC to zero ($OMEGA 2-9 all "0 FIX"), so those etas are absent here.

    etaiov_fdepot_1 ~ 0.103
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")     # Table 2 BOV in bioavailability 32.1% (24.9-39.5); 0.321^2 = 0.103; S10 $OMEGA 15 0.103
    etaiov_fdepot_2 ~ fixed(0.103)                                                                # OMEGA 16 equal to OMEGA 15 per S10 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.103)                                                                # OMEGA 17 equal to OMEGA 15 per S10 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.103)                                                                # OMEGA 18 equal to OMEGA 15 per S10 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fixed(0.103)                                                                # OMEGA 19 equal to OMEGA 15 per S10 $OMEGA BLOCK(1) SAME

    etaiov_ka_1 ~ 0.757
    label("Between-occasion variability in absorption rate constant, occasion 1 (log-scale variance)")  # Table 2 BOV in absorption rate constant 87.0% (64.6-121); 0.870^2 = 0.757
    etaiov_ka_2 ~ fixed(0.757)                                                                    # OMEGA 21 equal to OMEGA 20 per S10 $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.757)                                                                    # OMEGA 22 equal to OMEGA 20 per S10 $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fixed(0.757)                                                                    # OMEGA 23 equal to OMEGA 20 per S10 $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fixed(0.757)                                                                    # OMEGA 24 equal to OMEGA 20 per S10 $OMEGA BLOCK(1) SAME

    etaiov_mtt_1 ~ 1.93
    label("Between-occasion variability in mean absorption transit time, occasion 1 (log-scale variance)")  # Table 2 BOV in mean absorption transit time 139% (114-178); 1.39^2 = 1.93
    etaiov_mtt_2 ~ fixed(1.93)                                                                    # OMEGA 26 equal to OMEGA 25 per S10 $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(1.93)                                                                    # OMEGA 27 equal to OMEGA 25 per S10 $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fixed(1.93)                                                                    # OMEGA 28 equal to OMEGA 25 per S10 $OMEGA BLOCK(1) SAME
    etaiov_mtt_5 ~ fixed(1.93)                                                                    # OMEGA 29 equal to OMEGA 25 per S10 $OMEGA BLOCK(1) SAME

    # --- Residual error. The S10 $ERROR builds each additive term as
    # THETA + 0.2 * LLOQ with both THETAs FIX 0, so the additive standard
    # deviations are exactly 20% of the matrix-specific LLOQ, which is what
    # Table 2 footnote c describes. Plasma: 0.2 * 0.105 = 0.021 mg/L
    # (Table 2 prints 0.02). CSF: 0.2 * 0.0586 = 0.0117 mg/L (Table 2 prints
    # 0.01). See the vignette Errata.
    propSd <- 0.164
    label("Proportional residual error for plasma isoniazid (fraction)")                          # Table 2 Proportional error for plasma 16.4% (14.5-19.2)
    addSd <- fixed(0.021)
    label("Additive residual error for plasma isoniazid (mg/L)")                                  # Footnote c rule = 0.2 * 0.105 mg/L plasma LLOQ (S10 LLOQ_P = 0.105); Table 2 prints 0.02
    propSd_Ccsf <- 0.588
    label("Proportional residual error for CSF isoniazid (fraction)")                             # Table 2 Proportional error for CSF 58.8% (46.9-82.7)
    addSd_Ccsf <- fixed(0.0117)
    label("Additive residual error for CSF isoniazid (mg/L)")                                     # Footnote c rule = 0.2 * 0.0586 mg/L CSF LLOQ (S10 LLOQ_E = 0.0586); Table 2 prints 0.01
  })

  model({
    # 1. Occasion indicators. The S10 control stream multiplexes the
    # between-occasion etas with IF (OCC==k) blocks over five occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5

    # 2. NAT2 acetylator phenotype selects one of three independent typical
    # intrinsic clearances (S10 IF (NAT2_phenotype.EQ.0/1/2) THEN
    # TVCL = THETA(1)/THETA(2)/THETA(3)). The canonical encoding carries the
    # slow and rapid indicators; intermediate is the joint (0, 0) state.
    nat2_int <- (1 - NAT2_SLOW) * (1 - NAT2_RAPID)
    lclint_typ <- NAT2_SLOW * lclint_slow + nat2_int * lclint_int +
      NAT2_RAPID * lclint_rapid

    # 3. Individual parameters, all allometrically scaled on fat-free mass
    # against the 45 kg cohort median. Hepatic plasma flow uses its own
    # 56.1 kg reference (S10 ALLMCL_FFM_HEP).
    clint <- exp(lclint_typ + etalclint) * (FFM / 45)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 45)^e_ffm_vc
    vp <- exp(lvp) * (FFM / 45)^e_ffm_vc
    q <- exp(lq) * (FFM / 45)^e_ffm_cl
    qh <- exp(lqh) * (FFM / 56.1)^e_ffm_cl
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    ke0 <- exp(lke0)
    ppc <- exp(lppc)

    # 4. Well-stirred liver. EH is the fraction extracted on a single pass;
    # FH is the fraction of the absorbed dose surviving first pass, and
    # qh * eh is the systemic hepatic clearance, so k20 = qh * eh / vc.
    #   S10: EH = (CLINT*FU)/((CLINT*FU)+QH)
    #        FH = 1 - EH
    #        K20 = QH*EH/V
    eh <- (clint * fu) / (clint * fu + qh)
    fh <- 1 - eh
    k20 <- qh * eh / vc
    k23 <- q / vc
    k32 <- q / vp

    # 5. Savic transit-compartment absorption, written out in closed form.
    # The S10 $DES computes the gamma-density input as
    #   KTR     = (NN + 1) / MTT
    #   PIZZA   = LOG(BIO*PD*KTR) - GAMLN(NN+1)
    #   TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO)
    #   DADT(1) = TRANSIT - KA*A(1)
    # with PD the most recent dose amount and TEMPO the time after that
    # dose, guarded by IF (PD > 0 AND TEMPO > 0). podo(depot) and tad(depot)
    # supply PD and TEMPO. The guard is reproduced because at TEMPO = 0 the
    # NN*LOG(0) term is -Inf.
    #
    # This is written out rather than delegated to rxode2's transit()
    # built-in: transit() combined with the f(depot) <- 0 that the control
    # stream's F1 = 0 requires evaluates to an identically zero input rate,
    # so the model would simulate flat zero concentrations. The closed form
    # below is unaffected because podo() and tad() are both still live under
    # f(depot) <- 0.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) +
                    ntr * log(ktt) - ktt)
    }

    # 6. ODE system. First-pass extraction scales the depot-to-central
    # transfer; systemic hepatic elimination leaves the central compartment.
    Cc <- central / vc

    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot * fh - k20 * central -
      k23 * central + k32 * peripheral1
    d/dt(peripheral1) <- k23 * central - k32 * peripheral1
    # The CSF state holds a concentration, not an amount (supplementary S2).
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    # The dose amount is delivered entirely through the transit density, so
    # the ordinary bolus into the depot is suppressed (S10 $PK sets F1 = 0).
    f(depot) <- 0

    # 7. Observations and residual error.
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
