Wattanakul_2024_primaquine_motherinfant <- function() {
  description <- "Mechanistic mother-to-infant transfer model for primaquine and carboxyprimaquine, extending the maternal population PK model of Wattanakul 2024 (see modellib('Wattanakul_2024_primaquine')) with a breastfed-infant compartment chain so that infant exposure can be predicted from maternal dosing. The maternal side is unchanged apart from the milk compartments now being drained by breastfeeding. During each feeding window a square-wave gate transfers the entire breast-milk content of each analyte into the infant's dose compartment; between feeds the gate closes and the milk compartments re-equilibrate with maternal plasma. The infant absorbs through a two-transit chain with mean transit time fixed to a paediatric literature value, and has no first-pass metabolism. Infant clearances and volumes are scaled from the mother's individual estimates by infant body weight (allometric exponents 0.75 and 1), with monoamine-oxidase-A maturation applied to primaquine clearance only. No random effects are placed on the infant parameters: all predicted infant variability is inherited from the mother. Simulation model -- no infant pharmacokinetic parameter was estimated, because all but one infant primaquine sample was below the limit of quantification."
  reference <- paste(
    "Wattanakul T, Gilder ME, McGready R, Hanpithakpong W, Day NPJ,",
    "White NJ, Nosten F, Tarning J, Hoglund RM.",
    "Population pharmacokinetic modelling of primaquine exposures in",
    "lactating women and breastfed infants.",
    "Nat Commun. 2024;15:3851. doi:10.1038/s41467-024-47908-y.",
    "Structure and parameter values taken from the NONMEM control stream in",
    "the Supplementary Information ('NONMEM code') and the identical Zenodo",
    "deposit doi:10.5281/zenodo.10925291 ('Mother-to-infant.mod'),",
    "cross-checked against Table 2 of the main paper.",
    "Maternal layer shared with modellib('Wattanakul_2024_primaquine').",
    sep = " "
  )
  vignette <- "Wattanakul_2024_primaquine_lactation"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All amounts are nmol because the source model is
  # parameterised on the molar scale (control stream $INPUT "AMT ; DOSE
  # AMOUNT OF PRIMAQUINE (NMOL)"). verified = TRUE: analyte and specimen
  # were read off the $MODEL block comments of the source control stream,
  # which names all seventeen compartments explicitly.
  compartmentData <- list(
    depot                = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1             = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2             = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3             = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4             = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    central              = list(analyte = "primaquine", units = "nmol", specimen = "plasma", verified = TRUE),
    central_cpq          = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "plasma", verified = TRUE),
    milk                 = list(analyte = "primaquine", units = "nmol", specimen = "milk", verified = TRUE),
    milk_cpq             = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "milk", verified = TRUE),
    infant_depot         = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_transit1      = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_transit2      = list(analyte = "primaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_central       = list(analyte = "primaquine", units = "nmol", specimen = "plasma", verified = TRUE),
    infant_depot_cpq     = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_transit1_cpq  = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_transit2_cpq  = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "administration site", verified = TRUE),
    infant_central_cpq   = list(analyte = "carboxyprimaquine", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline. Reference WT = 51 kg = the median maternal weight of the study cohort (Wattanakul 2024 Table 1). Allometric exponents fixed a priori at 0.75 on CL/F and 1 on V/F for both analytes. WT also appears in the denominator of the mother-to-infant allometric scaling (WT_INFANT / WT), so the infant's parameters are scaled from the mother's INDIVIDUAL estimates, not from the typical values. Source column WT.",
      source_name        = "WT"
    ),
    WT_INFANT = list(
      description        = "Breastfed infant body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight of the mother's breastfed infant, i.e. of the DYAD PARTNER rather than of the modelled subject. Enters twice: it sizes the breast-milk compartments via Eq. 1, V_M = (0.15 L/kg/day * WT_INFANT) / feeds-per-day, and it scales the infant's clearances by (WT_INFANT / WT)^0.75 and volumes by (WT_INFANT / WT). Cohort median 6.8 kg (4.13-10.8) per Table 1; the dosing-scenario simulations of Figs. 4-6 sweep 2-17 kg paired with infant ages of 0-24 months, from the WHO weight-for-age standard growth curves at z-scores -3 to +3. Source column INFWT.",
      source_name        = "INFWT"
    ),
    AGE_INFANT = list(
      description        = "Breastfed infant postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Postnatal age of the breastfed infant. Drives the monoamine-oxidase-A maturation factor on the infant's primaquine clearance ONLY (MAO-A is the enzyme converting primaquine to carboxyprimaquine); the infant's carboxyprimaquine clearance carries no maturation term. The model derives postmenstrual age internally as AGE_INFANT + 9.2 months, the paper's full-term-gestation assumption of 40 weeks (Methods, 'Predicting infant concentrations', Eq. 4), because individual gestational ages were not used. Cohort infants were at least 28 days old, median 0.42 years = 5.0 months (0.13-1.81 years = 1.6-21.7 months) per Table 1. Source column INFAGE.",
      source_name        = "INFAGE"
    ),
    OCC = list(
      description        = "Integer-valued sampling-occasion indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4 identify the four maternal blood-sampling occasions -- day 0, day 3, day 7, and day 13 of the 14-day treatment course. Decomposed inside model() into binary indicators oc1 .. oc4 that multiplex the four IOV etas on maternal relative bioavailability and on maternal mean transit time, exactly as the source $PK block builds IOVF1 and IOVMTT from ETA(11)-ETA(18). The infant carries no IOV of its own. For single-occasion simulations pass OCC = 1 so the first IOV eta applies. Source column OCC.",
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
    dose_range     = "Primaquine 0.5 mg base/kg once daily for 14 days, given orally under non-fasting conditions as directly observed therapy. The published simulations additionally evaluated 1.0 mg base/kg once daily for 7 days, 0.5 mg base/kg twice daily for 7 days, and a 0.25 mg base/kg single dose.",
    regions        = "Thai-Myanmar border (three Shoklo Malaria Research Unit clinics), enrolled 11 November 2012 to 24 June 2014; ClinicalTrials.gov NCT01780753.",
    infant_partner = "n = 21 breastfed infants at least 28 days old: age 0.42 years (0.13-1.81), weight 6.8 kg (4.13-10.8), 14 male / 7 female (Wattanakul 2024 Table 1). Infant capillary sampling at 0, 2, 6 h after the first breastfeed following maternal dosing on day 0, at 0 and 2 h on days 3 and 7, and at 0, 2, 6, 24 h on day 13.",
    feeding_pattern = "Average 11 breastfeeds per day (6-18) and calculated infant daily milk intake 1020 mL (619-1620) per Table 1; the model's square-wave breastfeeding function uses the rounded average of 10 feeds/day (Results, 'Predicting infant concentrations').",
    notes          = "No infant pharmacokinetic parameter could be estimated: every infant capillary primaquine concentration was below the LLOQ except one sample at 2.59 ng/mL, and 67.5% of infant carboxyprimaquine concentrations were below the LLOQ (Methods, 'Predicting infant concentrations'). The infant layer is therefore entirely predictive -- structure and scaling assumptions carried over from the mother plus a literature maturation function -- and was evaluated by overlaying the 95% prediction interval on the observed infant data (Fig. 3), not by fitting. Treat infant predictions from this model as extrapolation, not as a fitted description of infant pharmacokinetics."
  )

  ini({
    # ------------------------------------------------------------------
    # Maternal layer. Identical to Wattanakul_2024_primaquine; see that
    # file for the full source trace. All values are FINAL estimates:
    # the supplement / Zenodo control stream is a `$SIM (123456) ONLYSIM
    # SUBPROBLEMS=1000` simulation deck whose twelve $THETA entries equal
    # Table 2's population estimates exactly.
    # ------------------------------------------------------------------
    lfdepot <- fixed(log(1));      label("Maternal relative oral bioavailability F (unitless)")                                          # Table 2 'F  1 fixed'; $THETA 1 '(1) FIX'
    lcl     <- log(17.1);          label("Apparent maternal primaquine elimination clearance CL/F_PQ at WT = 51 kg (L/h)")               # Table 2 'CL/F PQ (L/h) 17.1'; $THETA 2
    lvc     <- log(131);           label("Apparent maternal primaquine central volume of distribution V/F_PQ at WT = 51 kg (L)")         # Table 2 'V/F PQ (L) 131'; $THETA 3
    lcl_cpq <- log(0.967);         label("Apparent maternal carboxyprimaquine elimination clearance CL/F_CPQ at WT = 51 kg (L/h)")       # Table 2 'CL/F CPQ (L/h) 0.967'; $THETA 4
    lvc_cpq <- log(22.7);          label("Apparent maternal carboxyprimaquine central volume of distribution V/F_CPQ at WT = 51 kg (L)") # Table 2 'V/F CPQ (L) 22.7'; $THETA 5
    lmtt    <- log(1.44);          label("Maternal mean transit absorption time MTT (h)")                                                 # Table 2 'MTT (h) 1.44'; $THETA 6
    logitfm <- log(0.282 / (1 - 0.282)); label("Logit of the fraction of primaquine converted to carboxyprimaquine by maternal first-pass metabolism (FM = 0.282, unitless)")  # Table 2 'F_M 0.282'; $THETA 7

    lq_milk    <- log(0.400);      label("Apparent inter-compartmental clearance Q/F between the maternal central and breast-milk compartments, shared by primaquine and carboxyprimaquine (L/h)")  # Table 2 'Q/F PQ' and 'Q/F CPQ' footnote c; $THETA 8
    pcmilk     <- 0.376;           label("Milk:plasma partition coefficient PC_PQ, the fraction of primaquine freely distributed to breast milk (unitless)")           # Table 2 'PC PQ 0.376'; $THETA 9
    pcmilk_cpq <- 0.00889;         label("Milk:plasma partition coefficient PC_CPQ, the fraction of carboxyprimaquine freely distributed to breast milk (unitless)")   # Table 2 'PC CPQ 0.00889'; $THETA 10

    cfcap     <- 0.898;            label("Capillary:venous conversion factor CF_PQ for primaquine (unitless)")            # Table 2 'CF PQ 0.898'; $THETA 11
    cfcap_cpq <- 1.06;             label("Capillary:venous conversion factor CF_CPQ for carboxyprimaquine (unitless)")    # Table 2 'CF CPQ 1.06'; $THETA 12

    e_wt_cl <- fixed(0.75);        label("Allometric exponent on CL/F for both analytes, mother and infant (unitless)")   # Methods; $PK '(WT/51)**0.75' and '(INFWT/WT)**0.75'
    e_wt_vc <- fixed(1);           label("Allometric exponent on V/F for both analytes, mother and infant (unitless)")    # Methods; $PK '(WT/51)' and '(INFWT/WT)'

    # ------------------------------------------------------------------
    # Breastfeeding pattern and milk-to-infant transfer.
    # ------------------------------------------------------------------
    feed_n      <- fixed(10);      label("Number of breastfeeds per day, giving a 24/10 = 2.4 h feeding cycle (feeds/day)")   # Results, 'Predicting infant concentrations'; $PK 'FEEDNO=10'
    feed_window <- fixed(0.4);     label("Duration of each breastfeeding window (h; 0.4 h = 24 minutes)")                     # Results, 'Predicting infant concentrations'; $DES 'FEEDWINDOW=0.4'
    feed_first  <- fixed(1);       label("Time of the first breastfeed after time zero (h)")                                  # $DES 'FIRSTFEEDT=1'
    milk_intake <- fixed(0.15);    label("Average daily breast-milk intake per kg of infant body weight (L/kg/day)")          # Methods, Eq. 1 (150 mL/kg infant body weight); $PK 'VMPQ=(0.15*INFWT)/FEEDNO'
    kmilkinf    <- fixed(100);     label("First-order transfer rate from each breast-milk compartment to the corresponding infant dose compartment during a feed (1/h)")  # Results: fixed high so >95% of the milk-compartment amount transfers within the feeding window; $PK 'MTINFPQ=100', 'MTINFCPQ=100'

    # ------------------------------------------------------------------
    # Infant layer. Every value is fixed: no infant parameter was
    # estimable (all but one infant primaquine sample was below the LLOQ).
    # ------------------------------------------------------------------
    lmtt_infant <- fixed(log(0.706)); label("Infant mean transit absorption time (h)")                                        # Methods, 'Predicting infant concentrations': fixed to a paediatric literature model (reference 13) to account for more rapid absorption; $PK 'INFMTT=0.706'
    pma_tm50    <- fixed(7.6);        label("Postmenstrual age at 50% mature MAO-A-mediated infant primaquine clearance (months)")  # Methods, Eq. 4 and 'Predicting infant concentrations': derived from the literature to reflect 55% enzyme activity at full-term birth; $PK 'TM50=7.6'

    # ------------------------------------------------------------------
    # Maternal inter-individual variability. The infant carries none:
    # 'The IIV was not incorporated into the infant's pharmacokinetic
    # parameters. Therefore, the variation in predicted infant
    # concentration ... is solely attributed to the variability in the
    # mother's pharmacokinetic parameters' (Methods).
    # ------------------------------------------------------------------
    etalfdepot ~ 0.0243   # $OMEGA 1 IIV_F1;   sqrt(exp(0.0243)-1) = 15.7% CV = Table 2 'F IIV: 15.7'
    etalcl     ~ 0.0214   # $OMEGA 2 IIV_CLP;  sqrt(exp(0.0214)-1) = 14.7% CV vs Table 2 'CL/F PQ 15.1' -- control stream wins; see vignette Errata
    etalvc     ~ 0.0362   # $OMEGA 3 IIV_V2;   sqrt(exp(0.0362)-1) = 19.2% CV = Table 2 'V/F PQ 19.2'
    etalcl_cpq ~ 0.0688   # $OMEGA 4 IIV_CLM;  sqrt(exp(0.0688)-1) = 26.7% CV = Table 2 'CL/F CPQ 26.7'
    etalvc_cpq ~ 0.0307   # $OMEGA 5 IIV_V3;   sqrt(exp(0.0307)-1) = 17.7% CV = Table 2 'V/F CPQ 17.7'
    etalmtt    ~ 0.0411   # $OMEGA 6 IIV_MTT;  sqrt(exp(0.0411)-1) = 20.5% CV = Table 2 'MTT IIV: 20.5'
    etalogitfm ~ 0.163    # $OMEGA 7 IIV_FM;   sqrt(exp(0.163)-1)  = 42.1% CV = Table 2 'F_M 42.1'
    etalq_milk ~ 0.590    # $OMEGA 8 IIV_QPRQ; sqrt(exp(0.590)-1)  = 89.7% CV = Table 2 'Q/F PQ 89.7'

    # No IIV on PC_PQ, PC_CPQ, CF_PQ, CF_CPQ ($OMEGA 9 and 10 are
    # '0 FIX'; Table 2 reports '-' for all four). Omitted rather than
    # written as `~ fixed(0)` because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve.

    # Maternal inter-occasion variability across the four sampling
    # occasions. $OMEGA BLOCK(1) with SAME on occasions 2-4; nlmixr2 has
    # no SAME shortcut, so the variance is fixed to the shared value
    # after occasion 1.
    etaiov_fdepot_1 ~ 0.0413        # $OMEGA 11 IOV_F1_OCC1; sqrt(exp(0.0413)-1) = 20.5% CV = Table 2 'F IOV: 20.5'
    etaiov_fdepot_2 ~ fix(0.0413)   # $OMEGA 12 'BLOCK (1) SAME'
    etaiov_fdepot_3 ~ fix(0.0413)   # $OMEGA 13 'BLOCK (1) SAME'
    etaiov_fdepot_4 ~ fix(0.0413)   # $OMEGA 14 'BLOCK (1) SAME'
    etaiov_mtt_1    ~ 0.280         # $OMEGA 15 IOV_MTT_OCC1; sqrt(exp(0.280)-1) = 56.8% CV = Table 2 'MTT IOV: 56.8'
    etaiov_mtt_2    ~ fix(0.280)    # $OMEGA 16 'BLOCK (1) SAME'
    etaiov_mtt_3    ~ fix(0.280)    # $OMEGA 17 'BLOCK (1) SAME'
    etaiov_mtt_4    ~ fix(0.280)    # $OMEGA 18 'BLOCK (1) SAME'

    # ------------------------------------------------------------------
    # Residual error, maternal endpoints only. Additive on the
    # natural-log scale (nlmixr2 `~ lnorm(expSd)`); $SIGMA holds
    # VARIANCES, so expSd = sqrt(sigma^2). The source $ERROR block
    # defines no infant residual error, so the infant observables below
    # are returned as predictions without a residual-error model.
    # ------------------------------------------------------------------
    expSd           <- sqrt(0.102);   label("Residual SD on the natural-log scale for maternal primaquine in venous plasma (log units)")          # Table 2 'sigma PQ-venous 0.102'; $SIGMA 1
    expSd_Ccap      <- sqrt(0.0570);  label("Residual SD on the natural-log scale for maternal primaquine in capillary plasma (log units)")        # Table 2 'sigma PQ-capillary 0.0570'; $SIGMA 2
    expSd_cpq       <- sqrt(0.0198);  label("Residual SD on the natural-log scale for maternal carboxyprimaquine in venous plasma (log units)")    # Table 2 'sigma CPQ-venous 0.0198'; $SIGMA 3
    expSd_Ccap_cpq  <- sqrt(0.0115);  label("Residual SD on the natural-log scale for maternal carboxyprimaquine in capillary plasma (log units)") # Table 2 'sigma CPQ-capillary 0.0115'; $SIGMA 4
    expSd_Cmilk     <- sqrt(0.156);   label("Residual SD on the natural-log scale for primaquine in breast milk (log units)")                      # Table 2 'sigma PQ-breast milk 0.156'; $SIGMA 5
    expSd_Cmilk_cpq <- sqrt(0.0911);  label("Residual SD on the natural-log scale for carboxyprimaquine in breast milk (log units)")               # Table 2 'sigma CPQ-breast milk 0.0911'; $SIGMA 6
  })

  model({
    # ---- 1. Occasion indicators and maternal IOV ----------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
                  oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2 +
                  oc3 * etaiov_mtt_3    + oc4 * etaiov_mtt_4

    # ---- 2. Maternal individual parameters ----------------------------
    fdepot  <- exp(lfdepot + etalfdepot + iov_fdepot)
    cl      <- exp(lcl     + etalcl)     * (WT / 51)^e_wt_cl
    vc      <- exp(lvc     + etalvc)     * (WT / 51)^e_wt_vc
    cl_cpq  <- exp(lcl_cpq + etalcl_cpq) * (WT / 51)^e_wt_cl
    vc_cpq  <- exp(lvc_cpq + etalvc_cpq) * (WT / 51)^e_wt_vc
    mtt     <- exp(lmtt    + etalmtt + iov_mtt)
    q_milk  <- exp(lq_milk + etalq_milk)
    fm      <- exp(logitfm + etalogitfm) / (1 + exp(logitfm + etalogitfm))
    ktr     <- 5 / mtt
    vmilk   <- (milk_intake * WT_INFANT) / feed_n

    # ---- 3. Infant parameters -----------------------------------------
    # Postmenstrual age assuming full-term gestation at 9.2 months
    # (40 weeks), per Eq. 4 -- individual gestational ages were not used.
    pma_infant <- AGE_INFANT + 9.2
    # Hyperbolic (Hill exponent 1) MAO-A maturation. At birth
    # (AGE_INFANT = 0) this gives 9.2 / (7.6 + 9.2) = 0.548, the "55%
    # of adult activity at full-term birth" the TM50 was chosen to match.
    mf_infant <- pma_infant / (pma_tm50 + pma_infant)

    # Infant clearances and volumes are scaled from the MOTHER'S
    # INDIVIDUAL values (cl, vc, cl_cpq, vc_cpq already carry her etas),
    # which is how the infant inherits all of its variability.
    # Maturation is applied to primaquine clearance only: MAO-A is the
    # enzyme that converts primaquine to carboxyprimaquine.
    cl_infant     <- cl     * (WT_INFANT / WT)^e_wt_cl * mf_infant
    vc_infant     <- vc     * (WT_INFANT / WT)^e_wt_vc
    cl_infant_cpq <- cl_cpq * (WT_INFANT / WT)^e_wt_cl
    vc_infant_cpq <- vc_cpq * (WT_INFANT / WT)^e_wt_vc

    # Infant absorption: 2 transit compartments, so ktr = (2 + 1) / MTT.
    mtt_infant <- exp(lmtt_infant)
    ktr_infant <- 3 / mtt_infant

    # ---- 4. Micro-rate constants (1/h) --------------------------------
    kel     <- cl     / vc         # K23, maternal PQ central -> maternal CPQ central
    kel_cpq <- cl_cpq / vc_cpq     # K30, maternal CPQ central -> elimination

    k_central_milk     <- (q_milk / vc)     * pcmilk       # K28
    k_milk_central     <-  q_milk / vmilk                  # K82
    k_central_milk_cpq <- (q_milk / vc_cpq) * pcmilk_cpq   # K39
    k_milk_central_cpq <-  q_milk / vmilk                  # K93

    # As in the mother, the infant's entire systemic primaquine clearance
    # is a formation clearance into infant carboxyprimaquine (K16T17
    # appears with a + sign in DADT(17)).
    kel_infant     <- cl_infant     / vc_infant       # K16T17
    kel_infant_cpq <- cl_infant_cpq / vc_infant_cpq   # K17T0

    # ---- 5. Square-wave breastfeeding gate ----------------------------
    # Eqs. 5-8, transcribed from the $DES block (the published PDF's
    # display equations lose their minus signs on text extraction, so the
    # control stream is the authoritative form). sqwMilkToInfant = 1
    # during a feed and gates milk -> infant; sqwVenousToMilk = 1 between
    # feeds and gates plasma <-> milk. With feed_n = 10, feed_window =
    # 0.4 and feed_first = 1 the feed occupies t in [1.0, 1.4) h,
    # repeating every 2.4 h.
    #
    # SOURCE TYPO: the control stream writes "SQW2=(SQW-1)*(-1)",
    # referring to an undefined symbol SQW; SQW1 is meant. Eq. 8 defines
    # SQW_venous-to-milk as the complement of SQW_milk-to-infant, and no
    # other SQW symbol exists. The typo is identical in the Supplementary
    # Information listing and in the Zenodo deposit.
    feed_cycle <- 24 / feed_n
    feed_rest  <- feed_cycle - feed_window
    feed_psh1  <- 4 * pi * (feed_rest - feed_first) / feed_cycle
    feed_per1  <- 2 * pi *  feed_rest / feed_cycle
    sqwSin1    <- sin((pi - feed_per1) / 2)
    sqwSin2    <- sin(2 * pi * t / feed_cycle + (pi - feed_per1 + feed_psh1) / 2)
    # (|x| - x) / (2|x|) is 1 for x < 0 and 0 for x > 0. The published
    # form uses sqrt(x^2) for |x|, which is 0/0 exactly at the switching
    # instants; the 1e-12 offset keeps the expression finite there.
    sqwDiff          <- sqwSin2 - sqwSin1
    sqwAbs           <- sqrt(sqwDiff * sqwDiff + 1e-12)
    sqwMilkToInfant  <- (sqwAbs - sqwDiff) / (2 * sqwAbs)
    sqwVenousToMilk  <- 1 - sqwMilkToInfant

    # ---- 6. ODE system ------------------------------------------------
    # Maternal absorption and disposition (identical to
    # Wattanakul_2024_primaquine).
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

    # Breast-milk compartments. Unlike the maternal-only model, these are
    # now drained into the infant during each feeding window. That drain
    # is what makes the milk:plasma AUC ratio of this model fall below
    # the partition coefficient PC (about 0.28 versus PC_PQ = 0.376):
    # use Wattanakul_2024_primaquine, not this model, to reproduce the
    # milk:plasma AUC ratio reported in Table 3.
    d/dt(milk) <- k_central_milk * central * sqwVenousToMilk -
                  k_milk_central * milk    * sqwVenousToMilk -
                  kmilkinf * milk * sqwMilkToInfant

    d/dt(milk_cpq) <- k_central_milk_cpq * central_cpq * sqwVenousToMilk -
                      k_milk_central_cpq * milk_cpq    * sqwVenousToMilk -
                      kmilkinf * milk_cpq * sqwMilkToInfant

    # Infant primaquine: milk -> dose compartment -> 2 transits ->
    # central. Infant relative bioavailability is 1 ($PK "INFF1=1"), and
    # the dose arrives as an ODE flux rather than as a dose event, so no
    # f() term applies.
    d/dt(infant_depot)    <- kmilkinf * milk * sqwMilkToInfant -
                             ktr_infant * infant_depot
    d/dt(infant_transit1) <- ktr_infant * infant_depot    - ktr_infant * infant_transit1
    d/dt(infant_transit2) <- ktr_infant * infant_transit1 - ktr_infant * infant_transit2

    # Infant carboxyprimaquine ingested pre-formed in the milk. The
    # infant has NO first-pass metabolism (Methods: "the infant model was
    # assumed to be identical to the primaquine-carboxyprimaquine model
    # developed in the mothers, but ignoring the first-pass metabolism"),
    # so the two chains stay separate until their central compartments.
    d/dt(infant_depot_cpq)    <- kmilkinf * milk_cpq * sqwMilkToInfant -
                                 ktr_infant * infant_depot_cpq
    d/dt(infant_transit1_cpq) <- ktr_infant * infant_depot_cpq    - ktr_infant * infant_transit1_cpq
    d/dt(infant_transit2_cpq) <- ktr_infant * infant_transit1_cpq - ktr_infant * infant_transit2_cpq

    d/dt(infant_central) <- ktr_infant * infant_transit2 -
                            kel_infant * infant_central

    d/dt(infant_central_cpq) <- kel_infant * infant_central +
                                ktr_infant * infant_transit2_cpq -
                                kel_infant_cpq * infant_central_cpq

    # ---- 7. Bioavailability -------------------------------------------
    f(depot) <- fdepot

    # ---- 8. Observations and residual error ---------------------------
    # Six maternal endpoints with residual error, plus two infant
    # predictions without one (no infant sigma was estimable).
    Cc          <- central            / vc
    Ccap        <- Cc                 * cfcap
    Cc_cpq      <- central_cpq        / vc_cpq
    Ccap_cpq    <- Cc_cpq             * cfcap_cpq
    Cmilk       <- milk               / vmilk
    Cmilk_cpq   <- milk_cpq           / vmilk
    Cinfant     <- infant_central     / vc_infant
    Cinfant_cpq <- infant_central_cpq / vc_infant_cpq

    Cc        ~ lnorm(expSd)
    Ccap      ~ lnorm(expSd_Ccap)
    Cc_cpq    ~ lnorm(expSd_cpq)
    Ccap_cpq  ~ lnorm(expSd_Ccap_cpq)
    Cmilk     ~ lnorm(expSd_Cmilk)
    Cmilk_cpq ~ lnorm(expSd_Cmilk_cpq)
  })
}
