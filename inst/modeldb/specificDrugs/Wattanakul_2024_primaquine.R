Wattanakul_2024_primaquine <- function() {
  description <- "Mechanistic joint parent-metabolite population PK model for oral primaquine and its carboxyprimaquine metabolite in lactating women with Plasmodium vivax infection, fitted simultaneously to maternal venous plasma, capillary plasma, and breast-milk concentrations. Four-transit absorption (NN = 4, ka set equal to ktr) into a one-compartment primaquine disposition model whose entire systemic clearance is a formation clearance into a one-compartment carboxyprimaquine disposition model, plus first-pass conversion of a fraction FM of the absorbed dose directly into carboxyprimaquine. Each central compartment exchanges with its own breast-milk compartment through a shared apparent inter-compartmental clearance and an analyte-specific milk:plasma partition coefficient; the exchange is gated by a square-wave breastfeeding function (10 feeds/day, 24-minute feeding window in a 2.4-hour cycle). Separate capillary:venous conversion factors are estimated for each analyte. Allometric body-weight scaling on clearance (0.75) and volume (1) referenced at WT = 51 kg; inter-occasion variability on bioavailability and mean transit time across four sampling occasions. All amounts are molar (nmol) because the primaquine to carboxyprimaquine conversion is 1:1 molar."
  reference <- paste(
    "Wattanakul T, Gilder ME, McGready R, Hanpithakpong W, Day NPJ,",
    "White NJ, Nosten F, Tarning J, Hoglund RM.",
    "Population pharmacokinetic modelling of primaquine exposures in",
    "lactating women and breastfed infants.",
    "Nat Commun. 2024;15:3851. doi:10.1038/s41467-024-47908-y.",
    "Parameter values taken from the NONMEM control stream in the",
    "Supplementary Information ('NONMEM code') and the identical Zenodo",
    "deposit doi:10.5281/zenodo.10925291 ('Mother-to-infant.mod'),",
    "cross-checked against Table 2 of the main paper.",
    sep = " "
  )
  vignette <- "Wattanakul_2024_primaquine_lactation"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All amounts are nmol because the source model is
  # parameterised on the molar scale (Wattanakul 2024 Methods,
  # "Pharmacokinetics in breastfeeding women": "The natural logarithm of
  # observed molar primaquine and carboxyprimaquine concentrations was
  # analysed"; control stream $INPUT "AMT ; DOSE AMOUNT OF PRIMAQUINE
  # (NMOL)"). verified = TRUE: analyte and specimen were read off the
  # $MODEL block comments of the source control stream.
  compartmentData <- list(
    depot       = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "primaquine", units = "nmol", specimen = "plasma", verified = TRUE),
    central_cpq = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "plasma", verified = TRUE),
    milk        = list(analyte = "primaquine", units = "nmol", specimen = "milk", verified = TRUE),
    milk_cpq    = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "milk", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline. Reference WT = 51 kg = the median maternal weight of the study cohort (Wattanakul 2024 Table 1, 51 kg (35-81)). Allometric exponents fixed a priori at 0.75 on CL/F for both analytes and 1 on V/F for both analytes (Results paragraph 3: 'Allometric scaling of body weight was implemented a priori on both clearance and volume of distribution ... Although the allometric scaling of body weight did not improve the model fit, it was retained'). Source column WT.",
      source_name        = "WT"
    ),
    WT_INFANT = list(
      description        = "Breastfed infant body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight of the mother's breastfed infant, i.e. of the DYAD PARTNER rather than of the modelled subject. It enters this maternal model only through the volume of the breast-milk compartments, which the paper derives from the amount of milk the infant ingests per feed: V_M = (0.15 L/kg/day * WT_INFANT) / feeds-per-day (Wattanakul 2024 Eq. 1 and Methods, using an average daily milk intake of 150 mL/kg infant body weight). Cohort median 6.8 kg (4.13-10.8) per Table 1; the resulting milk-compartment volumes ranged 0.062-0.162 L (Results, 'Pharmacokinetics in breastfeeding women'). Source column INFWT.",
      source_name        = "INFWT"
    ),
    OCC = list(
      description        = "Integer-valued sampling-occasion indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4 identify the four blood-sampling occasions -- day 0, day 3, day 7, and day 13 of the 14-day treatment course (Wattanakul 2024 Methods, 'Pharmacokinetics in breastfeeding women': 'The IOV ... was ... evaluated ... between sampling occasions i.e., day 0, day 3, day 7, and day 13'). Decomposed inside model() into binary indicators oc1 .. oc4 that multiplex the four IOV etas on relative bioavailability and on mean transit time, exactly as the source $PK block builds IOVF1 and IOVMTT from ETA(11)-ETA(18). For single-occasion simulations pass OCC = 1 so the first IOV eta applies. Source column OCC.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21,
    n_studies      = 1,
    age_range      = "18-40 years",
    age_median     = "23 years",
    weight_range   = "35-81 kg",
    weight_median  = "51 kg",
    sex_female_pct = 100,
    disease_state  = "Lactating women with a history of Plasmodium vivax infection and no previous primaquine radical-cure treatment; G6PD-normal by fluorescent spot test and Mahidol-variant genotyping in both mother and infant.",
    dose_range     = "Primaquine 0.5 mg base/kg once daily for 14 days, given orally under non-fasting conditions as directly observed therapy.",
    regions        = "Thai-Myanmar border (three Shoklo Malaria Research Unit clinics), enrolled 11 November 2012 to 24 June 2014; ClinicalTrials.gov NCT01780753.",
    infant_partner = "Each mother had one breastfeeding infant at least 28 days old: n = 21, age 0.42 years (0.13-1.81), weight 6.8 kg (4.13-10.8), 14 male / 7 female (Wattanakul 2024 Table 1). Infant weight enters this model only through the breast-milk compartment volume.",
    feeding_pattern = "Average 11 breastfeeds per day (6-18) and calculated infant daily milk intake 1020 mL (619-1620) per Table 1; the model's square-wave breastfeeding function uses the rounded average of 10 feeds/day (Results, 'Predicting infant concentrations').",
    notes          = "Baseline demographics from Wattanakul 2024 Table 1. Sampling: dense venous sampling on days 0 and 13 (0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 24 h post-dose), sparse venous sampling on days 3 and 7 (pre-dose and 2 h); capillary sampling at 0, 2, 6, 12 h on days 0 and 13; breast milk by manual expression in the 1-3, 3-7, 7-12, and 12-24 h windows on days 0 and 13 plus one 1-3 h sample on days 3 and 7. LLOQ 1.14 ng/mL (primaquine) and 4.88 ng/mL (carboxyprimaquine) in all matrices, with relative standard error below 10% for all drug measurements."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. All values are the FINAL estimates: the
    # supplement / Zenodo control stream is a `$SIM (123456) ONLYSIM
    # SUBPROBLEMS=1000` simulation deck, and its twelve $THETA entries
    # equal Table 2's population estimates exactly.
    #
    # Apparent (per-bioavailability) clearances and volumes refer to a
    # woman weighing WT = 51 kg, the cohort median.
    # ------------------------------------------------------------------
    lfdepot <- fixed(log(1));      label("Relative oral bioavailability F (unitless)")                                                    # Table 2 'F  1 fixed'; $THETA 1 '(1) FIX'
    lcl     <- log(17.1);          label("Apparent primaquine elimination clearance CL/F_PQ at WT = 51 kg (L/h)")                        # Table 2 'CL/F PQ (L/h) 17.1 (5.85% RSE)'; $THETA 2
    lvc     <- log(131);           label("Apparent primaquine central volume of distribution V/F_PQ at WT = 51 kg (L)")                  # Table 2 'V/F PQ (L) 131 (7.26% RSE)'; $THETA 3
    lcl_cpq <- log(0.967);         label("Apparent carboxyprimaquine elimination clearance CL/F_CPQ at WT = 51 kg (L/h)")                # Table 2 'CL/F CPQ (L/h) 0.967 (7.56% RSE)'; $THETA 4
    lvc_cpq <- log(22.7);          label("Apparent carboxyprimaquine central volume of distribution V/F_CPQ at WT = 51 kg (L)")          # Table 2 'V/F CPQ (L) 22.7 (6.82% RSE)'; $THETA 5
    lmtt    <- log(1.44);          label("Mean transit absorption time MTT (h)")                                                          # Table 2 'MTT (h) 1.44 (9.43% RSE)'; $THETA 6

    # First-pass metabolism is estimated on the logit scale, exactly as
    # the source $PK block does (LTF = LOG(TVFM/(1-TVFM)) followed by
    # FM = EXP(LTF+ETA(7))/(1+EXP(LTF+ETA(7)))). That logit transform is
    # also what proves Table 2's 'F_M (%)' column header is wrong: 0.282
    # is a FRACTION, not a percentage.
    logitfm <- log(0.282 / (1 - 0.282)); label("Logit of the fraction of primaquine converted to carboxyprimaquine by first-pass metabolism (FM = 0.282, unitless)")  # Table 2 'F_M 0.282 (8.14% RSE)'; $THETA 7

    # Breast-milk transfer. Q/F is shared: the source states and encodes
    # Q/F_CPQ = Q/F_PQ (QCPQ = QPQ in $PK; Table 2 lists the same
    # 0.400 (19.4% RSE) estimate on both rows with footnote c).
    lq_milk    <- log(0.400);      label("Apparent inter-compartmental clearance Q/F between the central and breast-milk compartments, shared by primaquine and carboxyprimaquine (L/h)")  # Table 2 'Q/F PQ (L/h) 0.400' and 'Q/F CPQ (L/h) 0.400' footnote c; $THETA 8
    pcmilk     <- 0.376;           label("Milk:plasma partition coefficient PC_PQ, the fraction of primaquine freely distributed to breast milk (unitless)")           # Table 2 'PC PQ 0.376 (3.70% RSE)'; $THETA 9
    pcmilk_cpq <- 0.00889;         label("Milk:plasma partition coefficient PC_CPQ, the fraction of carboxyprimaquine freely distributed to breast milk (unitless)")   # Table 2 'PC CPQ 0.00889 (3.13% RSE)'; $THETA 10

    # Capillary:venous conversion factors. Note the Table 2 value for
    # primaquine (0.898) is used, not the 0.902 quoted in the Results
    # narrative -- the control stream $THETA 11 carries 0.898 and the
    # narrative value comes from the earlier venous + capillary model
    # fitted before the breast-milk data were introduced.
    cfcap     <- 0.898;            label("Capillary:venous conversion factor CF_PQ for primaquine (unitless)")             # Table 2 'CF PQ 0.898 (2.69% RSE)'; $THETA 11
    cfcap_cpq <- 1.06;             label("Capillary:venous conversion factor CF_CPQ for carboxyprimaquine (unitless)")     # Table 2 'CF CPQ 1.06 (1.25% RSE)'; $THETA 12

    # Allometric exponents, fixed a priori (not estimated).
    e_wt_cl <- fixed(0.75);        label("Allometric exponent on CL/F_PQ and CL/F_CPQ (unitless)")   # Methods, 'Pharmacokinetics in breastfeeding women'; $PK '(WT/51)**0.75'
    e_wt_vc <- fixed(1);           label("Allometric exponent on V/F_PQ and V/F_CPQ (unitless)")     # Methods, 'Pharmacokinetics in breastfeeding women'; $PK '(WT/51)'

    # Breastfeeding pattern. These four constants define the square-wave
    # gate on the plasma <-> breast-milk exchange and the milk-compartment
    # volume; all are fixed design constants of the simulation, not
    # estimated parameters.
    feed_n      <- fixed(10);      label("Number of breastfeeds per day, giving a 24/10 = 2.4 h feeding cycle (feeds/day)")   # Results, 'Predicting infant concentrations'; $PK 'FEEDNO=10'
    feed_window <- fixed(0.4);     label("Duration of each breastfeeding window (h; 0.4 h = 24 minutes)")                     # Results, 'Predicting infant concentrations'; $DES 'FEEDWINDOW=0.4'
    feed_first  <- fixed(1);       label("Time of the first breastfeed after time zero (h)")                                  # $DES 'FIRSTFEEDT=1'
    milk_intake <- fixed(0.15);    label("Average daily breast-milk intake per kg of infant body weight (L/kg/day)")          # Methods, Eq. 1 (150 mL/kg infant body weight); $PK 'VMPQ=(0.15*INFWT)/FEEDNO'

    # ------------------------------------------------------------------
    # Inter-individual variability. The $OMEGA block of a simulation deck
    # carries FINAL variances. Cross-check: %CV = sqrt(exp(omega^2) - 1)
    # reproduces Table 2's IIV column to three significant figures for
    # nine of the ten terms.
    # ------------------------------------------------------------------
    etalfdepot ~ 0.0243   # $OMEGA 1 IIV_F1;  sqrt(exp(0.0243)-1) = 15.7% CV = Table 2 'F IIV: 15.7'
    etalcl     ~ 0.0214   # $OMEGA 2 IIV_CLP; sqrt(exp(0.0214)-1) = 14.7% CV vs Table 2 'CL/F PQ 15.1' -- see the model-file note below and vignette Errata
    etalvc     ~ 0.0362   # $OMEGA 3 IIV_V2;  sqrt(exp(0.0362)-1) = 19.2% CV = Table 2 'V/F PQ 19.2'
    etalcl_cpq ~ 0.0688   # $OMEGA 4 IIV_CLM; sqrt(exp(0.0688)-1) = 26.7% CV = Table 2 'CL/F CPQ 26.7'
    etalvc_cpq ~ 0.0307   # $OMEGA 5 IIV_V3;  sqrt(exp(0.0307)-1) = 17.7% CV = Table 2 'V/F CPQ 17.7'
    etalmtt    ~ 0.0411   # $OMEGA 6 IIV_MTT; sqrt(exp(0.0411)-1) = 20.5% CV = Table 2 'MTT IIV: 20.5'
    etalogitfm ~ 0.163    # $OMEGA 7 IIV_FM;  sqrt(exp(0.163)-1)  = 42.1% CV = Table 2 'F_M 42.1'
    etalq_milk ~ 0.590    # $OMEGA 8 IIV_QPRQ; sqrt(exp(0.590)-1) = 89.7% CV = Table 2 'Q/F PQ 89.7'

    # No IIV on PC_PQ, PC_CPQ, CF_PQ, or CF_CPQ. The authors estimated
    # IIV on the two conversion factors, found it <= 10% with poor
    # precision (%RSE 45 and 139), and fixed it to zero with no
    # systematic bias in the goodness-of-fit (Results paragraph 3);
    # $OMEGA 9 and 10 (IIV_PC1, IIV_PC2) are '0 FIX'. Table 2 reports '-'
    # in the IIV column for all four. They are omitted here rather than
    # written as `~ fixed(0)` because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve.

    # ------------------------------------------------------------------
    # Inter-occasion variability across the four sampling occasions
    # (day 0, day 3, day 7, day 13). The source uses $OMEGA BLOCK(1)
    # with SAME on occasions 2-4, i.e. one shared variance per parameter.
    # nlmixr2 has no SAME shortcut, so each occasion gets its own eta
    # with the variance fixed to the shared value after occasion 1
    # (the Jonsson_2011_ethambutol.R / Aregbe_2012_alvespimycin.R
    # convention).
    # ------------------------------------------------------------------
    etaiov_fdepot_1 ~ 0.0413        # $OMEGA 11 IOV_F1_OCC1; sqrt(exp(0.0413)-1) = 20.5% CV = Table 2 'F IOV: 20.5'
    etaiov_fdepot_2 ~ fix(0.0413)   # $OMEGA 12 'BLOCK (1) SAME'
    etaiov_fdepot_3 ~ fix(0.0413)   # $OMEGA 13 'BLOCK (1) SAME'
    etaiov_fdepot_4 ~ fix(0.0413)   # $OMEGA 14 'BLOCK (1) SAME'
    etaiov_mtt_1    ~ 0.280         # $OMEGA 15 IOV_MTT_OCC1; sqrt(exp(0.280)-1) = 56.8% CV = Table 2 'MTT IOV: 56.8'
    etaiov_mtt_2    ~ fix(0.280)    # $OMEGA 16 'BLOCK (1) SAME'
    etaiov_mtt_3    ~ fix(0.280)    # $OMEGA 17 'BLOCK (1) SAME'
    etaiov_mtt_4    ~ fix(0.280)    # $OMEGA 18 'BLOCK (1) SAME'

    # ------------------------------------------------------------------
    # Residual error. The source $ERROR block fits Y = IPRED + EPS with
    # IPRED subsequently log-transformed, i.e. additive error on the
    # natural-log scale -- which is nlmixr2's `~ lnorm(expSd)`, not
    # `~ prop(propSd)`. $SIGMA holds VARIANCES (NONMEM definition), so
    # expSd = sqrt(sigma^2); the values below are written as sqrt() of
    # the published variance so the Table 2 number stays visible.
    #
    # The variance reading is confirmed by an assay-plausibility check:
    # reading the $SIGMA entries as SDs would imply a 1.15% residual
    # error on capillary carboxyprimaquine, which no LC-MS/MS assay
    # achieves, and the paper reports RSE below 10% for all drug
    # measurements (Methods, 'Drug quantification').
    # ------------------------------------------------------------------
    expSd           <- sqrt(0.102);   label("Residual SD on the natural-log scale for primaquine in venous plasma (log units)")            # Table 2 'sigma PQ-venous 0.102'; $SIGMA 1
    expSd_Ccap      <- sqrt(0.0570);  label("Residual SD on the natural-log scale for primaquine in capillary plasma (log units)")          # Table 2 'sigma PQ-capillary 0.0570'; $SIGMA 2
    expSd_cpq    <- sqrt(0.0198);  label("Residual SD on the natural-log scale for carboxyprimaquine in venous plasma (log units)")      # Table 2 'sigma CPQ-venous 0.0198'; $SIGMA 3
    expSd_Ccap_cpq  <- sqrt(0.0115);  label("Residual SD on the natural-log scale for carboxyprimaquine in capillary plasma (log units)")   # Table 2 'sigma CPQ-capillary 0.0115'; $SIGMA 4
    expSd_Cmilk     <- sqrt(0.156);   label("Residual SD on the natural-log scale for primaquine in breast milk (log units)")               # Table 2 'sigma PQ-breast milk 0.156'; $SIGMA 5
    expSd_Cmilk_cpq <- sqrt(0.0911);  label("Residual SD on the natural-log scale for carboxyprimaquine in breast milk (log units)")        # Table 2 'sigma CPQ-breast milk 0.0911'; $SIGMA 6
  })

  model({
    # ---- 1. Occasion indicators and inter-occasion variability --------
    # Reproduces the source $PK block:
    #   IOVF1  = OCC1*ETA(11) + OCC2*ETA(12) + OCC3*ETA(13) + OCC4*ETA(14)
    #   IOVMTT = OCC1*ETA(15) + OCC2*ETA(16) + OCC3*ETA(17) + OCC4*ETA(18)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
                  oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2 +
                  oc3 * etaiov_mtt_3    + oc4 * etaiov_mtt_4

    # ---- 2. Individual parameters -------------------------------------
    fdepot  <- exp(lfdepot + etalfdepot + iov_fdepot)
    cl      <- exp(lcl     + etalcl)     * (WT / 51)^e_wt_cl
    vc      <- exp(lvc     + etalvc)     * (WT / 51)^e_wt_vc
    cl_cpq  <- exp(lcl_cpq + etalcl_cpq) * (WT / 51)^e_wt_cl
    vc_cpq  <- exp(lvc_cpq + etalvc_cpq) * (WT / 51)^e_wt_vc
    mtt     <- exp(lmtt    + etalmtt + iov_mtt)
    q_milk  <- exp(lq_milk + etalq_milk)

    # Inverse-logit back-transform of the first-pass metabolised fraction.
    fm <- exp(logitfm + etalogitfm) / (1 + exp(logitfm + etalogitfm))

    # Transit-chain rate constant. NN = 4 transit compartments and ka is
    # set equal to ktr, so KTR = (NN + 1) / MTT = 5 / MTT.
    ktr <- 5 / mtt

    # Volume of each breast-milk compartment (Eq. 1): the milk ingested
    # per feed = daily intake per kg * infant weight / feeds per day.
    # The same volume applies to primaquine and carboxyprimaquine.
    vmilk <- (milk_intake * WT_INFANT) / feed_n

    # ---- 3. Micro-rate constants (1/h) --------------------------------
    # kel is a pure FORMATION clearance: the source routes 100% of
    # systemic primaquine elimination into the carboxyprimaquine central
    # compartment (K23 appears with a + sign in DADT(3)), so no
    # primaquine leaves the system except as carboxyprimaquine.
    kel     <- cl     / vc         # K23, PQ central -> CPQ central
    kel_cpq <- cl_cpq / vc_cpq     # K30, CPQ central -> elimination

    k_central_milk         <- (q_milk / vc)     * pcmilk       # K28
    k_milk_central         <-  q_milk / vmilk                  # K82
    k_central_milk_cpq     <- (q_milk / vc_cpq) * pcmilk_cpq   # K39
    k_milk_central_cpq     <-  q_milk / vmilk                  # K93

    # ---- 4. Square-wave breastfeeding gate ----------------------------
    # Eqs. 5-8, transcribed from the $DES block of the control stream
    # (the display equations in the published PDF lose their minus signs
    # when the PDF text layer is extracted, so the control stream is the
    # authoritative form). The two waves are complementary:
    # sqwMilkToInfant = 1 during a feed, sqwVenousToMilk = 1 between
    # feeds. With feed_n = 10, feed_window = 0.4 and feed_first = 1 the
    # gate opens on t in [1.0, 1.4) h and repeats every 2.4 h.
    #
    # NOTE ON A SOURCE TYPO: the control stream writes
    # "SQW2=(SQW-1)*(-1)", referring to an undefined symbol SQW. It is
    # SQW1 that is meant -- Eq. 8 of the paper defines SQW_venous-to-milk
    # as the complement of SQW_milk-to-infant, and no other SQW symbol
    # exists. The typo is present identically in the Supplementary
    # Information listing and in the Zenodo deposit.
    feed_cycle <- 24 / feed_n
    feed_rest  <- feed_cycle - feed_window
    feed_psh1  <- 4 * pi * (feed_rest - feed_first) / feed_cycle
    feed_per1  <- 2 * pi *  feed_rest / feed_cycle
    sqwSin1    <- sin((pi - feed_per1) / 2)
    sqwSin2    <- sin(2 * pi * t / feed_cycle + (pi - feed_per1 + feed_psh1) / 2)
    # (|x| - x) / (2|x|) is 1 for x < 0 and 0 for x > 0. The published
    # form uses sqrt(x^2) for |x|, which is 0/0 exactly at the switching
    # instants; the 1e-12 offset keeps the expression finite there
    # without measurably softening the edge.
    sqwDiff          <- sqwSin2 - sqwSin1
    sqwAbs           <- sqrt(sqwDiff * sqwDiff + 1e-12)
    sqwMilkToInfant  <- (sqwAbs - sqwDiff) / (2 * sqwAbs)
    sqwVenousToMilk  <- 1 - sqwMilkToInfant

    # ---- 5. ODE system ------------------------------------------------
    # Compartment amounts are nmol; volumes are L; concentrations are
    # nmol/L. Absorption: depot -> transit1 -> ... -> transit4, then the
    # last transit splits by first-pass metabolism, K72 = KTR*(1-FM) into
    # primaquine central and K73 = KTR*FM into carboxyprimaquine central.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4

    d/dt(central) <- ktr * (1 - fm) * transit4 -
                     kel * central -
                     k_central_milk * central * sqwVenousToMilk +
                     k_milk_central * milk    * sqwVenousToMilk

    d/dt(central_cpq) <- ktr * fm * transit4 +
                         kel * central -
                         kel_cpq * central_cpq -
                         k_central_milk_cpq * central_cpq * sqwVenousToMilk +
                         k_milk_central_cpq * milk_cpq    * sqwVenousToMilk

    # Breast-milk compartments. In this maternal model milk exchanges
    # with plasma but is not drained by the infant: the milk-to-infant
    # transfer term of the mother-to-infant model
    # (Wattanakul_2024_primaquine_motherinfant) is absent here. That is
    # what makes the milk:plasma AUC ratio equal the partition
    # coefficient PC exactly, reproducing Table 3's 0.376 (0.375-0.377);
    # including the drain would drive the ratio to about 0.28.
    d/dt(milk) <- k_central_milk * central * sqwVenousToMilk -
                  k_milk_central * milk    * sqwVenousToMilk
    d/dt(milk_cpq) <- k_central_milk_cpq * central_cpq * sqwVenousToMilk -
                      k_milk_central_cpq * milk_cpq    * sqwVenousToMilk

    # ---- 6. Bioavailability -------------------------------------------
    f(depot) <- fdepot

    # ---- 7. Observations and residual error ---------------------------
    # Six simultaneously-fitted endpoints: each analyte in venous plasma,
    # capillary plasma, and breast milk. Capillary predictions are the
    # venous prediction scaled by the analyte's conversion factor, per
    # the $ERROR block (IPRED = (A(2)/S2)*CFPQ for CMT 2 / TYPE 2).
    Cc        <- central     / vc
    Ccap      <- Cc          * cfcap
    Cc_cpq    <- central_cpq / vc_cpq
    Ccap_cpq  <- Cc_cpq      * cfcap_cpq
    Cmilk     <- milk        / vmilk
    Cmilk_cpq <- milk_cpq    / vmilk

    Cc        ~ lnorm(expSd)
    Ccap      ~ lnorm(expSd_Ccap)
    Cc_cpq    ~ lnorm(expSd_cpq)
    Ccap_cpq  ~ lnorm(expSd_Ccap_cpq)
    Cmilk     ~ lnorm(expSd_Cmilk)
    Cmilk_cpq ~ lnorm(expSd_Cmilk_cpq)
  })
}
