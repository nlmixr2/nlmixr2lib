Svensson_2014_bedaquiline <- function() {
  description <- "Three-compartment population PK model for bedaquiline (BDQ) and a two-compartment N-desmethyl metabolite M2 in healthy adult volunteers following single 400 mg oral doses, with Savic 2007 analytical transit-compartment absorption (non-integer NN feeding a first-order depot at rate ka), fixed allometric scaling on disposition (0.75 on CL/Q at 70 kg, 1 on Vc/Vp), and multiplicative rifampicin or rifapentine drug-drug-interaction factors of 4.78 and 3.96 on bedaquiline and M2 apparent clearance, applied at full induction from day 3 of rifamycin co-administration."
  reference <- paste(
    "Svensson E. M., Murray S., Karlsson M. O., Dooley K. E. (2014).",
    "Rifampicin and rifapentine significantly reduce concentrations of",
    "bedaquiline, a new anti-TB drug.",
    "Journal of Antimicrobial Chemotherapy 70(4):1106-1114.",
    "doi:10.1093/jac/dku504.",
    sep = " "
  )
  vignette <- "Svensson_2014_bedaquiline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (used for allometric scaling around 70 kg)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed body weight. Allometric power on apparent clearances (CL/F, Q1/F, Q2/F, CLM2/(F*fm), Q_M2/(F*fm)) at exponent 0.75 and on apparent volumes (V/F, VP1/F, VP2/F, V_M2/(F*fm), VP_M2/(F*fm)) at exponent 1 around a reference 70 kg adult (Svensson 2014 Methods 'Model development': 'Allometric scaling was applied to CL and V using body weight and fixed coefficients of 0.75 and 1, respectively.').",
      source_name        = "WT"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampicin co-administration at full CYP3A4 induction (1 = on daily 600 mg rifampicin for at least 3 days; 0 = not on rifampicin or pre-induction lag).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on rifampicin or first <3 days of co-administration)",
      notes              = "Subject- and time-varying indicator that the subject has reached the post-induction equilibrium with rifampicin co-administration. Svensson 2014 found that parameterising the bedaquiline and M2 apparent clearances as switching instantaneously to their induced values after 3 days of rifamycin administration provided the best objective-function fit ('Parameterizing the CLs to change after 3 days of rifamycin administration provided the best fit based on objective function value, and the magnitude of the estimated interaction effect remained similar over the evaluated range of timepoints for onset.'). Multiplicative factor on CL_BDQ and CL_M2: cl_eff = cl_base * 4.78^CONMED_RIF (Svensson 2014 Table 2 'Factor change BDQ/M2 CL with RIF = 4.78'). For simulation, set CONMED_RIF = 1 on observation rows that fall >= 3 days after the start of rifampicin co-administration and 0 otherwise. The same indicator drives the multiplicative effect on M2 CL because a separate factor for M2 did not improve the fit.",
      source_name        = "RIF"
    ),
    CONMED_RPT = list(
      description        = "Concomitant rifapentine co-administration at full CYP3A4 induction (1 = on daily 600 mg rifapentine for at least 3 days; 0 = not on rifapentine or pre-induction lag).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on rifapentine or first <3 days of co-administration)",
      notes              = "Subject- and time-varying indicator that the subject has reached the post-induction equilibrium with rifapentine co-administration. Same 3-day-lag-then-instantaneous-switch parameterisation as CONMED_RIF (Svensson 2014 Methods 'Model development'). Multiplicative factor on CL_BDQ and CL_M2: cl_eff = cl_base * 3.96^CONMED_RPT (Svensson 2014 Table 2 'Factor change BDQ/M2 CL with RPT = 3.96'). The same indicator drives the multiplicative effect on M2 CL because a separate factor for M2 did not improve the fit.",
      source_name        = "RPT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "19-55 years (study cohort) per Svensson 2014 Table 1",
    age_median     = "35.5 years (overall); 38 years (rifampicin arm), 34 years (rifapentine arm)",
    weight_range   = "57.3-122 kg (overall) per Svensson 2014 Table 1",
    weight_median  = "81.8 kg (overall); 80.5 kg (rifampicin arm), 83.3 kg (rifapentine arm)",
    sex_female_pct = 12.5,
    race_ethnicity = c(
      White                          = 87.5,
      `Black or African American`    = 6.2,
      `American Indian or Alaska Native` = 3.1,
      Asian                          = 3.1
    ),
    disease_state  = "Healthy adult volunteers in a Phase I two-arm open-label, two-period single-sequence drug-drug-interaction study (TMC207-CL002) assessing the effect of multiple-dose rifampicin (Arm 1, n = 13 completers) or rifapentine (Arm 2, n = 16 completers) at 600 mg daily on bedaquiline pharmacokinetics.",
    dose_range     = "Single 400 mg oral dose of bedaquiline on Day 1 (period 1, alone) and again on Day 29 (period 2, after 9 days of rifamycin pre-treatment with rifamycin co-administration continuing through the 14-day PK sampling). Bedaquiline given as one tablet with PK sampling at 0, 1, 2, 3, 4, 5, 6, 8, 12, 24 h, then every 24 h until 336 h (14 days) after each dose, plus a sample on Day 20 before rifamycin start.",
    regions        = "Single Phase I site; geographic detail not stated in publication.",
    studies        = "Combined fit of TMC207-CL002 (32 healthy volunteers, two arms, 1419 bedaquiline + 1419 M2 observations across two single 400 mg doses each) with historical bedaquiline DDI data from a similarly designed efavirenz interaction study in a comparable population (1083 bedaquiline + 1055 M2 observations after a single 400 mg dose alone or with efavirenz).",
    notes          = "Baseline demographics from Svensson 2014 Table 1. Three premature discontinuations all in the rifampicin arm (two failure-to-comply, one adverse event with unlikely drug relation). Bedaquiline and M2 concentrations were analysed in molar units in the original NONMEM fit (molecular weights 555.50 g/mol for bedaquiline and 541.47 g/mol for M2); the apparent parameters in Table 2 are mass-balance-equivalent to mg/L apparent disposition because volumes and clearances were estimated against ratios of molar concentrations. The bundled paper figures plot concentrations in nmol/L; multiply mg/L by 1000/MW to convert (BDQ: x1.800 nmol/L per mg/L; M2: x1.847 nmol/L per mg/L)."
  )

  ini({
    # Final parameter estimates from Svensson 2014 Table 2 'Final model parameter
    # estimates with precision obtained from non-parametric bootstrap (n=1000)'.
    # The model structure (3-compartment BDQ + 2-compartment M2 with NN transit
    # compartments + first-order ka absorption) was inherited from a previously
    # developed bedaquiline DDI model (Svensson 2014 reference 19), but all final
    # parameter values were re-estimated on the new TMC207-CL002 data merged
    # with the prior efavirenz-DDI dataset; values below are from this paper's
    # Table 2.

    # Structural absorption: Savic 2007 transit-compartment model (non-integer
    # number of transit compartments NN feeding a depot at rate KTR = (NN+1)/MTT,
    # which then absorbs into central at rate KA).
    lmtt    <- log(0.97)    ; label("Mean transit time MTT (h, on log scale)")                                   # Svensson 2014 Table 2 'MTT = 0.97 h' (RSE 11.5%)
    lka     <- log(0.12)    ; label("First-order absorption rate constant ka (1/h, on log scale)")               # Svensson 2014 Table 2 'KA = 0.12 1/h' (RSE 3.9%)
    lnn     <- log(8.41)    ; label("Number of transit compartments NN (unitless, on log scale)")                # Svensson 2014 Table 2 'NN = 8.41' (RSE 36.2%)

    # Structural disposition: 3-compartment bedaquiline (parent) +
    # 2-compartment M2 metabolite. Apparent volumes and clearances (i.e.,
    # CL/F, V/F for the parent and CL/(F*fm), V/(F*fm) for the metabolite,
    # where fm is the fraction of bedaquiline metabolised to M2 -- not
    # separately identifiable from the F-relative parameters).
    lcl     <- log(3.20)    ; label("Apparent bedaquiline clearance CL/F at 70 kg (L/h)")                         # Svensson 2014 Table 2 'CL/F = 3.20 L/h' (RSE 6.5%)
    lvc     <- log(16.2)    ; label("Apparent bedaquiline central volume V/F at 70 kg (L)")                       # Svensson 2014 Table 2 'V/F = 16.2 L' (RSE 12.9%)
    lq      <- log(4.71)    ; label("Apparent bedaquiline first-peripheral inter-compartmental clearance Q1/F at 70 kg (L/h)") # Svensson 2014 Table 2 'Q1/F = 4.71 L/h' (RSE 5.6%)
    lvp     <- log(2801)    ; label("Apparent bedaquiline first-peripheral volume VP1/F at 70 kg (L)")             # Svensson 2014 Table 2 'VP1/F = 2801 L' (RSE 10.1%)
    lq2     <- log(3.10)    ; label("Apparent bedaquiline second-peripheral inter-compartmental clearance Q2/F at 70 kg (L/h)") # Svensson 2014 Table 2 'Q2/F = 3.10 L/h' (RSE 6.0%)
    lvp2    <- log(137)     ; label("Apparent bedaquiline second-peripheral volume VP2/F at 70 kg (L)")            # Svensson 2014 Table 2 'VP2/F = 137 L' (RSE 10.4%)
    lcl_m2  <- log(13.1)    ; label("Apparent M2 metabolite clearance CLM2/(F*fm) at 70 kg (L/h)")                 # Svensson 2014 Table 2 'CLM2/F/fm = 13.1 L/h' (RSE 6.6%)
    lvc_m2  <- log(882)     ; label("Apparent M2 metabolite central volume V_M2/(F*fm) at 70 kg (L)")              # Svensson 2014 Table 2 'VM2/F/fm = 882 L' (RSE 4.9%)
    lq_m2   <- log(105)     ; label("Apparent M2 metabolite inter-compartmental clearance Q_M2/(F*fm) at 70 kg (L/h)") # Svensson 2014 Table 2 'Q1M2/F/fm = 105 L/h' (RSE 8.6%)
    lvp_m2  <- log(3349)    ; label("Apparent M2 metabolite peripheral volume VP_M2/(F*fm) at 70 kg (L)")          # Svensson 2014 Table 2 'VP1M2/F/fm = 3349 L' (RSE 3.8%)

    # Bioavailability is implicitly F = 1 because CL and V are reported as
    # apparent F-relative values (CL/F, V/F). The transit() function below
    # passes F = 1 directly; no `lfdepot` parameter is needed.

    # Allometric scaling (theory-based, fixed): exponent 0.75 on apparent
    # clearances and 1 on apparent volumes around 70 kg (Svensson 2014 Methods
    # 'Model development': 'Allometric scaling was applied to CL and V using
    # body weight and fixed coefficients of 0.75 and 1, respectively.').
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on CL and Q (CL/F, Q1/F, Q2/F, CL_M2, Q_M2; unitless)")  # Svensson 2014 Methods Model development (fixed)
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent on Vc and Vp (V/F, VP1/F, VP2/F, V_M2, VP_M2; unitless)") # Svensson 2014 Methods Model development (fixed)

    # Drug-drug-interaction multiplicative factors on the bedaquiline and M2
    # apparent clearances when on rifampicin (Factor RIF) or rifapentine
    # (Factor RPT) at full induction (from day 3 of co-administration onward).
    # The paper reports the SAME factor applies to bedaquiline CL and M2 CL
    # ('a separate fixed effect parameter for the change in M2 CL did not
    # improve the fit of the model significantly'), so a single factor per
    # rifamycin drives both apparent clearances.
    e_rif_cl <- 4.78  ; label("Multiplicative factor on bedaquiline and M2 apparent CL during rifampicin co-administration at full induction (unitless)")  # Svensson 2014 Table 2 'Factor change BDQ/M2 CL with RIF = 4.78' (RSE 9.1%)
    e_rpt_cl <- 3.96  ; label("Multiplicative factor on bedaquiline and M2 apparent CL during rifapentine co-administration at full induction (unitless)") # Svensson 2014 Table 2 'Factor change BDQ/M2 CL with RPT = 3.96' (RSE 5.0%)

    # Inter-individual variability. Final estimates converted from the paper's
    # %CV notation via omega^2 = log(1 + CV^2) for log-normal IIV. Svensson
    # 2014 Table 2 reports a NONMEM BLOCK(6) on the bedaquiline-CL eta, the
    # M2-CL eta, and four individual induction-effect etas (one each for the
    # subject-specific rifampicin effect on BDQ and M2, and the rifapentine
    # effect on BDQ and M2; 6x6 lower triangular correlation matrix). The
    # nlmixr2lib implementation here encodes the 2x2 BSV(CL, CL_M2) block and
    # the diagonal IIVs on V, Q1, V_M2, VP_M2, but DROPS the four
    # individual-induction etas (BSV RIF BDQ, BSV RIF M2, BSV RPT BDQ, BSV
    # RPT M2) and their inter-block correlations with BSV CL / BSV CLM2
    # because nlmixr2lib has no idiomatic encoding for an eta that is gated
    # by a binary covariate (the induction etas only apply when the subject
    # is in the respective arm). The paper also reports BOV F = 18.5%, BOV
    # MTT = 64.7%, and BSV F = 12.1% (between-occasion variability is the
    # paper's primary stochastic component on the two dosing occasions); BOV
    # has no idiomatic encoding in nlmixr2lib and is dropped here. See the
    # vignette's Assumptions and deviations section for the full list of
    # dropped random effects with the original paper values.
    etalcl + etalcl_m2 ~ c(0.07654, 0.04523, 0.09236)  # Svensson 2014 Table 2 BSV CL = 28.2% CV, BSV CLM2 = 31.1% CV, correlation 53.8%; var_cl = log(1 + 0.282^2) = 0.07654, cov = 0.538 * sqrt(0.07654 * 0.09236) = 0.04523, var_clm2 = log(1 + 0.311^2) = 0.09236
    etalvc      ~ 0.16769       # var = log(1 + 0.428^2) = 0.16769 ; Svensson 2014 Table 2 BSV V = 42.8% CV
    etalq       ~ 0.04372       # var = log(1 + 0.211^2) = 0.04372 ; Svensson 2014 Table 2 BSV Q1 = 21.1% CV
    etalvc_m2   ~ 0.12379       # var = log(1 + 0.364^2) = 0.12379 ; Svensson 2014 Table 2 BSV VM2 = 36.4% CV
    etalvp_m2   ~ 0.04330       # var = log(1 + 0.210^2) = 0.04330 ; Svensson 2014 Table 2 BSV VP1M2 = 21.0% CV

    # Residual error. The paper reports proportional residual errors on both
    # bedaquiline (15.7% CV) and M2 (12.2% CV) plus a 55% correlation between
    # the two residuals and a 2.19-fold weighting of the residual SD on
    # observations between 0 and 6 h after dose (Svensson 2014 Table 2
    # 'Weighting of samples 0-6 h = 2.19'). The cross-output residual
    # correlation and the early-sample weighting are dropped here because
    # nlmixr2lib has no idiomatic encoding for either; see vignette
    # Assumptions and deviations.
    propSd     <- 0.157  ; label("Bedaquiline residual error (proportional SD, fraction)")  # Svensson 2014 Table 2 'Prop err BDQ = 15.7% CV' (RSE 4.4%)
    propSd_m2  <- 0.122  ; label("M2 metabolite residual error (proportional SD, fraction)") # Svensson 2014 Table 2 'Prop err M2 = 12.2% CV' (RSE 5.0%)
  })

  model({
    # 1. Individual PK parameters with allometric scaling on apparent
    #    clearances (CL/F, Q1/F, Q2/F, CL_M2, Q_M2) and apparent volumes
    #    (V/F, VP1/F, VP2/F, V_M2, VP_M2) around a 70 kg reference adult.
    allcl <- (WT / 70)^e_wt_cl_q
    allv  <- (WT / 70)^e_wt_vc_vp

    # 2. Drug-drug-interaction multiplicative factor on bedaquiline and M2
    #    apparent CL. CONMED_RIF / CONMED_RPT are time- and subject-varying
    #    binary indicators that the subject has reached full rifamycin
    #    induction (from day 3 of co-administration onward). When neither
    #    indicator is on, the factor is 1.
    ddi <- e_rif_cl^CONMED_RIF * e_rpt_cl^CONMED_RPT

    # 3. Individual PK parameters with IIV (etas applied multiplicatively
    #    on the log-scale around the typical-value log-transformed structural
    #    parameter).
    mtt   <- exp(lmtt)
    ka    <- exp(lka)
    nn    <- exp(lnn)

    cl    <- exp(lcl    + etalcl)    * allcl * ddi
    vc    <- exp(lvc    + etalvc)    * allv
    q     <- exp(lq     + etalq)     * allcl
    vp    <- exp(lvp)                * allv
    q2    <- exp(lq2)                * allcl
    vp2   <- exp(lvp2)               * allv
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * allcl * ddi
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * allv
    q_m2  <- exp(lq_m2)              * allcl
    vp_m2 <- exp(lvp_m2 + etalvp_m2) * allv

    # 4. ODE system. Savic 2007 analytical transit-compartment input feeds
    #    the depot at rate (NN+1)/MTT (via the rxode2 transit() built-in)
    #    and the depot then absorbs into bedaquiline central at rate ka.
    #    Three-compartment bedaquiline disposition (central + peripheral1
    #    + peripheral2) with linear elimination feeding the metabolite
    #    central_m2 compartment; M2 has its own two-compartment disposition
    #    (central_m2 + peripheral1_m2) with linear elimination from
    #    central_m2. The transit() function reads dose amount from
    #    podo(depot); f(depot) <- 0 suppresses the bolus contribution so
    #    the chain delivers the full dose smoothly through the transit
    #    compartments.
    d/dt(depot)          <- transit(nn, mtt, 1) - ka * depot
    d/dt(central)        <-  ka * depot -
                              cl  * central / vc -
                              q   * central / vc + q  * peripheral1 / vp -
                              q2  * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)    <-  q   * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)    <-  q2  * central / vc - q2 * peripheral2 / vp2
    d/dt(central_m2)     <-  cl  * central / vc -
                              cl_m2 * central_m2 / vc_m2 -
                              q_m2  * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2) <-  q_m2  * central_m2 / vc_m2 - q_m2 * peripheral1_m2 / vp_m2

    # 5. Bioavailability. F is fixed at 1 (the F = 1 anchor is in ini()).
    #    f(depot) <- 0 suppresses bolus into depot; transit() drives the
    #    input rate over the absorption window.
    f(depot) <- 0

    # 6. Observed plasma concentrations. Compartmental amount (mg) / volume
    #    (L) = concentration in mg/L (= ug/mL). To convert to nmol/L for
    #    comparison with the paper figures, multiply by 1000/MW (BDQ: x1.800
    #    nmol/L per mg/L; M2: x1.847 nmol/L per mg/L).
    Cc    <- central    / vc
    Cc_m2 <- central_m2 / vc_m2

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
  })
}
