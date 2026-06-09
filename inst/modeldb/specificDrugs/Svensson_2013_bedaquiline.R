Svensson_2013_bedaquiline <- function() {
  description <- "Three-compartment population PK model for bedaquiline (BDQ) with a two-compartment N-desmethyl metabolite M2 and a two-compartment N,N-bis-desmethyl metabolite M3 in healthy adult volunteers following single 400 mg oral doses, with Savic 2007 analytical transit-compartment absorption (non-integer NN feeding a first-order depot at rate ka) and an instantaneous-switch concomitant-efavirenz induction factor of 2.07 on apparent CL_BDQ and CL_M2 and 1.12 on apparent CL_M3, applied from 1 week after the start of 600 mg once-nightly efavirenz co-administration."
  reference <- paste(
    "Svensson E. M., Aweeka F., Park J.-G., Marzan F., Dooley K. E.,",
    "Karlsson M. O. (2013).",
    "Model-based estimates of the effects of efavirenz on bedaquiline",
    "pharmacokinetics and suggested dose adjustments for patients",
    "coinfected with HIV and tuberculosis.",
    "Antimicrobial Agents and Chemotherapy 57(6):2780-2787.",
    "doi:10.1128/AAC.00191-13.",
    sep = " "
  )
  vignette <- "Svensson_2013_bedaquiline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (used for allometric scaling around 70 kg)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline body weight. Allometric scaling applied with fixed exponents 0.75 on apparent clearances (CL/F, Q1/F, Q2/F, CL_M2, Q_M2, CL_M3, Q_M3) and 1 on apparent volumes (V/F, VP1/F, VP2/F, V_M2, VP_M2, V_M3, VP_M3) around a 70 kg reference adult. Svensson 2013 Methods 'Population modeling' states 'Allometric scaling of all disposition compartments with body weight as the size descriptor and fixed (0.75 for clearance and 1 for volume of distribution V) or estimated coefficients was evaluated' and Results 'Allometric scaling of disposition parameters with fixed coefficients improved the fit markedly (0.75 for clearances and 1 for volumes; estimation of the coefficients did not significantly improve the fit further)'. The reference weight is not explicitly stated in Svensson 2013, but 70 kg is the universally adopted allometric reference (Anderson and Holford 2008) and is the value used by the sibling Svensson 2014 bedaquiline DDI paper that extends this analysis; see vignette Assumptions and deviations.",
      source_name        = "WT"
    ),
    CONMED_EFV = list(
      description        = "Concomitant efavirenz co-administration at full CYP3A4 induction (1 = on daily 600 mg efavirenz for at least 1 week; 0 = not on efavirenz or pre-induction lag).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on efavirenz or first <1 week of co-administration)",
      notes              = "Subject- and time-varying indicator that the subject has reached the post-induction equilibrium with efavirenz co-administration. Svensson 2013 modelled the effect of EFV on BDQ, M2, and M3 apparent clearances as an instantaneous change in CL one week after initialization of EFV treatment (Svensson 2013 Methods 'Population modeling': 'the effect was modeled as an instantaneous change in clearance (CL) and/or bioavailability (F). Several time points between day 1 and day 14 of EFV treatment were tested for this change'; Results 'The impact of induction was described as an instantaneous change in clearance 1 week after initialization of EFV treatment'). Multiplicative factor on CL_BDQ and CL_M2: cl_eff = cl_base * 2.07^CONMED_EFV (Svensson 2013 Table 3 'EFV EFF BDQ and M2 = 2.07'). Multiplicative factor on CL_M3: cl_eff = cl_base * 1.12^CONMED_EFV (Svensson 2013 Table 3 'EFV EFF M3 = 1.12'). For simulation, set CONMED_EFV = 1 on observation rows that fall >= 1 week after the start of 600 mg once-nightly efavirenz co-administration and 0 otherwise (no EFV, or within the pre-induction lag). The same indicator drives the multiplicative effect on M2 CL because the increase in clearance with induction was not significantly different for BDQ and M2 (Svensson 2013 Results 'The increase in clearance with induction was not significantly different for BDQ and M2 and estimated to be about 2-fold').",
      source_name        = "EFV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_studies      = 1L,
    age_range      = "19-62 years",
    age_median     = "44 years",
    weight_range   = "57.3-118.8 kg",
    weight_median  = "82.3 kg",
    sex_female_pct = 8.6,
    race_ethnicity = c(
      `White, non-Hispanic`             = 68.6,
      `Black, non-Hispanic`             = 22.9,
      `Hispanic (regardless of race)`   = 5.7,
      `Asian or Pacific Islander`       = 2.9
    ),
    bmi_range     = "19.0-36.1 kg/m^2",
    bmi_median    = "25.9 kg/m^2",
    disease_state  = "Healthy adult volunteers (18-65 years) enrolled at four AIDS Clinical Trials Group (ACTG) sites with no clinical evidence of TB, negative HIV antibody test, normal standard blood tests, and normal corrected QT intervals. Women of reproductive potential were excluded.",
    dose_range     = "Single 400 mg oral dose of bedaquiline on Day 1 (period 1, alone) after a standard 670 kcal / 33% fat breakfast. From day 15 to 42, subjects received 600 mg efavirenz fasting each evening. On Day 29 (period 2, after 2 weeks of daily EFV), a second 400 mg oral dose of bedaquiline was administered. PK sampling: predose and at 1, 2, 3, 4, 5, 6, 8, 12, 24, 48, 72, 120, 168, 216, 264, and 336 h after each dose (i.e., up to day 14 / day 43).",
    regions        = "Four ACTG sites, United States",
    cyp2b6_status  = c(Slow = 8.6, Intermediate = 37.1, Extensive = 54.3),
    notes          = "Baseline demographics from Svensson 2013 Table 2. Thirty-seven healthy subjects were enrolled in ACTG study A5267; 35 completed the first PK sampling and were included in the modelling analysis (33 completed both PK sampling periods). The dataset comprises 1,152 observations each of BDQ and M2 and 560 observations of M3 (the M3 sampling subset). CYP2B6 composite 516/983 metabolizer genotype was tested as a covariate but was not significantly correlated with the magnitude of the EFV induction effect (Svensson 2013 Results). Few non-white subjects were enrolled; race was not investigated as a covariate."
  )

  ini({
    # Final parameter estimates from Svensson 2013 Table 3 'Parameter estimates
    # with relative standard errors of the final model'. Reported in apparent
    # F-relative form because bioavailability F was fixed to 1 in the original
    # analysis: structural CL and V are CL/F and V/F for the parent and
    # CL/(F*fm) / V/(F*fm) for the metabolites, where fm is the fraction of
    # parent metabolised to the next species in the cascade and is not
    # separately identifiable from the F-relative parameters.

    # Structural absorption: Savic 2007 dynamic transit-compartment model
    # (non-integer number of transit compartments NN feeding a depot at rate
    # KTR = (NN+1)/MTT, which then absorbs into central at rate KA).
    lmtt    <- log(1.31)    ; label("Mean transit time MTT (h, on log scale)")                                       # Svensson 2013 Table 3 'MTT = 1.31 h' (RSE 12.6%)
    lka     <- log(0.128)   ; label("First-order absorption rate constant ka (1/h, on log scale)")                    # Svensson 2013 Table 3 'KA = 0.128 1/h' (RSE 8.7%)
    lnn     <- log(5.21)    ; label("Number of transit compartments NN (unitless, on log scale)")                     # Svensson 2013 Table 3 'NN = 5.21' (RSE 20.5%)

    # Structural disposition: 3-compartment bedaquiline (parent),
    # 2-compartment M2 metabolite (formed from BDQ via CYP3A4 N-demethylation),
    # 2-compartment M3 metabolite (formed from M2 via further demethylation;
    # responsible enzyme(s) not directly identified in vitro but suspected
    # CYP3A4 by analogy with the BDQ -> M2 step). Apparent volumes and
    # clearances are CL/F, V/F for the parent and CL/(F*fm), V/(F*fm) for the
    # metabolites.
    lcl     <- log(2.96)    ; label("Apparent bedaquiline clearance CL/F at 70 kg (L/h)")                              # Svensson 2013 Table 3 'CL = 2.96 L/h' (RSE 9.5%)
    lvc     <- log(17.3)    ; label("Apparent bedaquiline central volume V/F at 70 kg (L)")                            # Svensson 2013 Table 3 'V = 17.3 L' (RSE 18.7%)
    lq      <- log(5.01)    ; label("Apparent bedaquiline first-peripheral inter-compartmental clearance Q1/F at 70 kg (L/h)") # Svensson 2013 Table 3 'Q1 = 5.01 L/h' (RSE 8.3%)
    lvp     <- log(2870)    ; label("Apparent bedaquiline first-peripheral volume VP1/F at 70 kg (L)")                  # Svensson 2013 Table 3 'VP1 = 2,870 L' (RSE 15.3%)
    lq2     <- log(4.16)    ; label("Apparent bedaquiline second-peripheral inter-compartmental clearance Q2/F at 70 kg (L/h)") # Svensson 2013 Table 3 'Q2 = 4.16 L/h' (RSE 10.2%)
    lvp2    <- log(136)     ; label("Apparent bedaquiline second-peripheral volume VP2/F at 70 kg (L)")                 # Svensson 2013 Table 3 'VP2 = 136 L' (RSE 9.0%)
    lcl_m2  <- log(12.3)    ; label("Apparent M2 metabolite clearance CL_M2/(F*fm_M2) at 70 kg (L/h)")                  # Svensson 2013 Table 3 'CL M2 = 12.3 L/h' (RSE 10.1%)
    lvc_m2  <- log(659)     ; label("Apparent M2 metabolite central volume V_M2/(F*fm_M2) at 70 kg (L)")                # Svensson 2013 Table 3 'V M2 = 659 L' (RSE 7.2%)
    lq_m2   <- log(103)     ; label("Apparent M2 metabolite inter-compartmental clearance Q1_M2/(F*fm_M2) at 70 kg (L/h)") # Svensson 2013 Table 3 'Q1 M2 = 103 L/h' (RSE 10.5%)
    lvp_m2  <- log(2840)    ; label("Apparent M2 metabolite peripheral volume VP1_M2/(F*fm_M2) at 70 kg (L)")           # Svensson 2013 Table 3 'VP1 M2 = 2,840 L' (RSE 6.0%)
    lcl_m3  <- log(39.2)    ; label("Apparent M3 metabolite clearance CL_M3/(F*fm_M2*fm_M3) at 70 kg (L/h)")            # Svensson 2013 Table 3 'CL M3 = 39.2 L/h' (RSE 9.0%)
    lvc_m3  <- log(11.2)    ; label("Apparent M3 metabolite central volume V_M3/(F*fm_M2*fm_M3) at 70 kg (L)")          # Svensson 2013 Table 3 'V M3 = 11.2 L' (RSE 44.7%)
    lq_m3   <- log(106)     ; label("Apparent M3 metabolite inter-compartmental clearance Q_M3/(F*fm_M2*fm_M3) at 70 kg (L/h)") # Svensson 2013 Table 3 'Q M3 = 106 L/h' (RSE 9.9%)
    lvp_m3  <- log(2680)    ; label("Apparent M3 metabolite peripheral volume VP_M3/(F*fm_M2*fm_M3) at 70 kg (L)")      # Svensson 2013 Table 3 'VP M3 = 2,680 L' (RSE 13.6%)

    # Bioavailability is implicitly F = 1 because CL and V are reported as
    # apparent F-relative values. Svensson 2013 Table 3 footnote b explicitly
    # states 'Estimated with typical value of F fixed to 1'.

    # Allometric scaling (theory-based, fixed): exponent 0.75 on apparent
    # clearances and 1 on apparent volumes around 70 kg.
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on apparent CL and Q (CL/F, Q1/F, Q2/F, CL_M2, Q_M2, CL_M3, Q_M3; unitless)")  # Svensson 2013 Results 'allometric scaling of disposition parameters with fixed coefficients ... 0.75 for clearances' (fixed)
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent on apparent Vc and Vp (V/F, VP1/F, VP2/F, V_M2, VP_M2, V_M3, VP_M3; unitless)") # Svensson 2013 Results 'allometric scaling of disposition parameters with fixed coefficients ... 1 for volumes' (fixed)

    # Drug-drug-interaction multiplicative factors on the bedaquiline, M2, and
    # M3 apparent clearances when on efavirenz at full induction (from 1 week
    # of co-administration onward). Svensson 2013 estimated a SHARED factor
    # for BDQ and M2 because 'the increase in clearance with induction was
    # not significantly different for BDQ and M2' (Results), and a SEPARATE
    # smaller factor for M3.
    e_efv_cl    <- 2.07  ; label("Multiplicative factor on bedaquiline and M2 apparent CL during efavirenz co-administration at full induction (unitless)") # Svensson 2013 Table 3 'EFV EFF BDQ and M2 = 2.07' (RSE 3.6%)
    e_efv_cl_m3 <- 1.12  ; label("Multiplicative factor on M3 apparent CL during efavirenz co-administration at full induction (unitless)")              # Svensson 2013 Table 3 'EFV EFF M3 = 1.12' (RSE 3.6%)

    # Inter-individual variability (between-subject, BSV). Final estimates
    # converted from the paper's CV% notation via omega^2 = log(1 + CV^2) for
    # log-normal IIV. Svensson 2013 Table 3 reports a NONMEM BLOCK(6) on the
    # bedaquiline-CL eta, the M2-CL eta, the M3-CL eta, and three EFV-effect
    # etas (one each for the subject-specific EFV induction effect on BDQ,
    # M2, and M3). The 6x6 lower-triangular correlation matrix is given by
    # the off-diagonal '%c' entries in Table 3 (footnote c: 'Correlation
    # estimated as a covariance; presented here structured as the omega or
    # sigma blocks used in the NONMEM code'). The nlmixr2lib implementation
    # here encodes the 3x3 BLOCK on (etalcl, etalcl_m2, etalcl_m3) and the
    # remaining diagonal IIVs on V, Q1, V_M2, VP_M2; but DROPS the three
    # EFV-effect etas (BSV EFV EFF-BDQ = 20.6%, EFF-M2 = 28.2%, EFF-M3 =
    # 32.7%) and all their off-diagonal correlations with the structural-CL
    # block because nlmixr2lib has no idiomatic encoding for an eta that is
    # gated by a binary covariate (the EFV-effect etas only apply when the
    # subject is in the post-induction window). The paper also reports
    # BOV F = 23.6%, BSV F = 24.3%, and BOV MTT = 55.4%; the between-occasion
    # variance components have no idiomatic encoding in nlmixr2lib and are
    # dropped here. See the vignette's Assumptions and deviations section
    # for the full list of dropped random effects with the original paper
    # values.
    #
    # BLOCK(3) on (CL, CL_M2, CL_M3): var_CL = log(1 + 0.237^2) = 0.054646;
    # var_CL_M2 = log(1 + 0.188^2) = 0.034733; var_CL_M3 = log(1 + 0.30^2) =
    # 0.086178. Off-diagonals from Svensson 2013 Table 3 correlations:
    # corr(CL, CL_M2)    = +0.296  -> cov = 0.296 * sqrt(0.054646*0.034733) = +0.012897
    # corr(CL, CL_M3)    = -0.129  -> cov = -0.129 * sqrt(0.054646*0.086178) = -0.008854
    # corr(CL_M2, CL_M3) = +0.713  -> cov = 0.713 * sqrt(0.034733*0.086178) = +0.039008
    # Lower-triangular order (NONMEM BLOCK): row 1; row 2; row 3.
    etalcl + etalcl_m2 + etalcl_m3 ~ c(0.054646,
                                       0.012897, 0.034733,
                                      -0.008854, 0.039008, 0.086178)  # Svensson 2013 Table 3 BSV CL = 23.7%, BSV CL_M2 = 18.8% (corr with CL = +29.6%), BSV CL_M3 = 30.0% (corr with CL = -12.9%, corr with CL_M2 = +71.3%)
    etalvc      ~ 0.113081  # var = log(1 + 0.346^2) = 0.113081 ; Svensson 2013 Table 3 BSV V = 34.6% CV (RSE 32%)
    etalq       ~ 0.034372  # var = log(1 + 0.187^2) = 0.034372 ; Svensson 2013 Table 3 BSV Q1 = 18.7% CV (RSE 15%)
    etalvc_m2   ~ 0.080214  # var = log(1 + 0.289^2) = 0.080214 ; Svensson 2013 Table 3 BSV V M2 = 28.9% CV (RSE 19%)
    etalvp_m2   ~ 0.064902  # var = log(1 + 0.259^2) = 0.064902 ; Svensson 2013 Table 3 BSV VP 1M2 = 25.9% CV (RSE 39%)

    # Residual error. The paper reports proportional residual errors on BDQ
    # (23.9% CV), M2 (17.7% CV), and M3 (15.0% CV) plus correlations across
    # residual errors for the three compounds (BDQ-M2 corr +14.9%, BDQ-M3
    # corr +7.5%, M2-M3 corr +11.7%; Svensson 2013 Table 3 'Prop error'
    # rows) and time-varying weighting factors of 1.87 for samples drawn in
    # the absorption phase (TAD < 6 h) and 3.28 for observations below the
    # limit of quantification (Svensson 2013 Table 3 'Error wt TAD < 6 h' /
    # 'Error wt < BLQ'). The cross-output residual correlations, the
    # absorption-phase weighting, and the BLQ weighting are all dropped here
    # because nlmixr2lib has no idiomatic encoding for any of them; see
    # vignette Assumptions and deviations.
    propSd     <- 0.239  ; label("Bedaquiline residual error (proportional SD, fraction)")  # Svensson 2013 Table 3 'Prop error BDQ = 23.9% CV' (RSE 5.3%)
    propSd_m2  <- 0.177  ; label("M2 metabolite residual error (proportional SD, fraction)") # Svensson 2013 Table 3 'Prop error M2 = 17.7% CV' (RSE 4.8%)
    propSd_m3  <- 0.150  ; label("M3 metabolite residual error (proportional SD, fraction)") # Svensson 2013 Table 3 'Prop error M3 = 15.0% CV' (RSE 10.2%)
  })

  model({
    # 1. Individual PK parameters with allometric scaling on apparent
    #    clearances (CL/F, Q1/F, Q2/F, CL_M2, Q_M2, CL_M3, Q_M3) and apparent
    #    volumes (V/F, VP1/F, VP2/F, V_M2, VP_M2, V_M3, VP_M3) around a 70 kg
    #    reference adult.
    allcl <- (WT / 70)^e_wt_cl_q
    allv  <- (WT / 70)^e_wt_vc_vp

    # 2. Drug-drug-interaction multiplicative factor on bedaquiline, M2, and
    #    M3 apparent CL. CONMED_EFV is a time- and subject-varying binary
    #    indicator that the subject has reached full efavirenz induction
    #    (from 1 week of co-administration onward). When CONMED_EFV = 0 the
    #    factor is 1 (no induction). The BDQ/M2 factor is the same (2.07);
    #    the M3 factor is smaller (1.12).
    ddi_bdq_m2 <- e_efv_cl^CONMED_EFV
    ddi_m3     <- e_efv_cl_m3^CONMED_EFV

    # 3. Individual PK parameters with IIV (etas applied multiplicatively on
    #    the log-scale around the typical-value log-transformed structural
    #    parameter).
    mtt   <- exp(lmtt)
    ka    <- exp(lka)
    nn    <- exp(lnn)

    cl    <- exp(lcl    + etalcl)    * allcl * ddi_bdq_m2
    vc    <- exp(lvc    + etalvc)    * allv
    q     <- exp(lq     + etalq)     * allcl
    vp    <- exp(lvp)                * allv
    q2    <- exp(lq2)                * allcl
    vp2   <- exp(lvp2)               * allv
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * allcl * ddi_bdq_m2
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * allv
    q_m2  <- exp(lq_m2)              * allcl
    vp_m2 <- exp(lvp_m2 + etalvp_m2) * allv
    cl_m3 <- exp(lcl_m3 + etalcl_m3) * allcl * ddi_m3
    vc_m3 <- exp(lvc_m3)             * allv
    q_m3  <- exp(lq_m3)              * allcl
    vp_m3 <- exp(lvp_m3)             * allv

    # 4. ODE system. Savic 2007 analytical transit-compartment input feeds
    #    the depot at rate (NN+1)/MTT (via the rxode2 transit() built-in)
    #    and the depot then absorbs into bedaquiline central at rate ka.
    #    Three-compartment bedaquiline disposition (central + peripheral1 +
    #    peripheral2) with linear elimination feeding the M2 central
    #    compartment; M2 has its own two-compartment disposition
    #    (central_m2 + peripheral1_m2) with linear elimination feeding the
    #    M3 central compartment; M3 has its own two-compartment disposition
    #    (central_m3 + peripheral1_m3) with linear elimination from
    #    central_m3. The transit() function reads dose amount from
    #    podo(depot); f(depot) <- 0 suppresses the bolus contribution so the
    #    chain delivers the full dose smoothly through the transit
    #    compartments.
    d/dt(depot)          <- transit(nn, mtt, 1) - ka * depot
    d/dt(central)        <-  ka * depot -
                              cl  * central / vc -
                              q   * central / vc + q  * peripheral1 / vp -
                              q2  * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)    <-  q   * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)    <-  q2  * central / vc - q2 * peripheral2 / vp2
    d/dt(central_m2)     <-  cl    * central    / vc -
                              cl_m2 * central_m2 / vc_m2 -
                              q_m2  * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2) <-  q_m2  * central_m2 / vc_m2 - q_m2 * peripheral1_m2 / vp_m2
    d/dt(central_m3)     <-  cl_m2 * central_m2 / vc_m2 -
                              cl_m3 * central_m3 / vc_m3 -
                              q_m3  * central_m3 / vc_m3 + q_m3 * peripheral1_m3 / vp_m3
    d/dt(peripheral1_m3) <-  q_m3  * central_m3 / vc_m3 - q_m3 * peripheral1_m3 / vp_m3

    # 5. Bioavailability. F is fixed at 1 (Svensson 2013 Table 3 footnote b).
    #    f(depot) <- 0 suppresses the bolus contribution; transit() drives
    #    the input rate over the absorption window.
    f(depot) <- 0

    # 6. Observed plasma concentrations. Compartmental amount (mg) / volume
    #    (L) = concentration in mg/L (= ug/mL).
    Cc    <- central    / vc
    Cc_m2 <- central_m2 / vc_m2
    Cc_m3 <- central_m3 / vc_m3

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
    Cc_m3 ~ prop(propSd_m3)
  })
}
