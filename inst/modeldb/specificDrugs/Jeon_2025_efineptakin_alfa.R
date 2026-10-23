Jeon_2025_efineptakin_alfa <- function() {
  description <- "Semi-mechanistic population PK model for efineptakin alfa (rhIL-7-hyFc, GX-I7 / NT-I7; a homodimer of two IL-7 moieties fused to a hybrid IgD/IgG4 Fc, 104 kDa) after a single SC or IM microdose in healthy adults. Three states: the dose enters an injection-site / interstitial absorption space (apparent volume V_I, route-specific) that exchanges bidirectionally with central (Q) and is taken up (k_uptake) into a lumped FcRn-recycling endosomal space (apparent volume V_D). In the endosome the drug binds FcRn at quasi-steady state and the complex leaves at k_recycle, a fraction FR = 0.715 returning to central and the rest to the injection site. In central the drug binds IL-7R (CD127) at quasi-steady state (QSS-TMDD) with internalisation k_int and is cleared linearly (CL). Amounts are pmol of IL-7 (two per drug molecule); the observation is the free central IL-7 concentration (pmol/L)."
  reference <- "Jeon HS, Lee SW, Jung W, Jung H, Choi D, Byun M-S, Yun H-y, Lee S, Kim A, Chae J-w, Lee H. Dose Optimization of rhIL-7-hyFc for Patients With Lymphopenia Using a Neonatal Fc Receptor-Mediated Recycling-Based and Target-Mediated Drug Disposition Pharmacokinetic Model. Clin Transl Sci. 2025;18(5):e70252. doi:10.1111/cts.70252. PMCID: PMC12094080."
  vignette <- "Jeon_2025_efineptakin_alfa"
  units <- list(
    time = "h",
    dosing = "pmol (of rhIL-7-hyFc; f(abs_site) = 2 converts to pmol IL-7)",
    concentration = "pmol/L (IL-7)"
  )

  compartmentData <- list(
    abs_site = list(
      analyte = "efineptakin alfa (IL-7 moiety)",
      units = "pmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "efineptakin alfa (IL-7 moiety)", units = "pmol", specimen = "serum", verified = TRUE),
    endosome = list(analyte = "efineptakin alfa (IL-7 moiety)", units = "pmol", specimen = "endosome", verified = TRUE)
  )

  covariateData <- list(
    ROUTE_SC = list(
      description = "Injection route: 1 = subcutaneous, 0 = intramuscular",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IM)",
      notes = "Assigned route (single-dose SC and IM cohorts). Switches the apparent injection-site volume V_I between its SC (1680 L) and IM (442 L) estimates, and through K12 = Q/V_I the rate of transfer from the injection site to central. The control stream encodes the same switch as 'VI = VISC; IF(COHORT.EQ.3) VI = VIIM' (cohort 3 = 60 ug/kg IM).",
      source_name = "COHORT (3 = IM)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling (p < 0.05 forward, p < 0.01 backward); not retained (Jeon 2025 Results 3.1). Control-stream $INPUT column AGE."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened by stepwise covariate modelling; not retained (Jeon 2025 Results 3.1). Control-stream $INPUT column SEX (coding not stated)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained. The authors attribute the absence of a weight effect to weight-based (ug/kg) dosing (Jeon 2025 Results 3.1). Control-stream $INPUT column WTKG. Weight enters only through the ug/kg dose."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained (Jeon 2025 Results 3.1). Control-stream $INPUT column HTCM."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained (Jeon 2025 Results 3.1). Control-stream $INPUT column BMI."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    n_observations = 410,
    study_id = "NCT02860715 (phase 1, randomized, double-blind, placebo-controlled; Seoul National University Hospital)",
    age_range = "19-45 years (eligibility); mean 27.1 +/- 6.0 years (Lee 2020)",
    weight_range = "mean 65.3 +/- 11.1 kg (Lee 2020)",
    sex_female_pct = NA,
    race_ethnicity = "Korean (single-centre study in the Republic of Korea; not tabulated)",
    disease_state = "Healthy adult volunteers",
    dose_range = "Single dose of rhIL-7-hyFc 20 ug/kg SC, 60 ug/kg SC or 60 ug/kg IM (8 active subjects per cohort; 30 randomized 8:2 to drug or placebo)",
    regions = "Republic of Korea",
    notes = "Demographics are not reported in Jeon 2025, which refers to the primary trial report (Lee SW et al., Clin Transl Sci 2020;13:1161-1169, doi:10.1111/cts.12800); mean age, height, weight and BMI above are from that report (men and women eligible; sex split not given). Concentrations are baseline-adjusted serum IL-7 (pre-dose mean subtracted), assumed by the authors to be unbound IL-7. LLOQ 0.031 ng/mL. The model was used to simulate multiple-dose IM regimens up to 10,800 ug/kg, well above the 60 ug/kg maximum studied."
  )

  ini({
    # Structural parameters. Every value is the final estimate in the Text S1
    # control stream ($THETA); each equals the 'Observed value - Median' column of
    # Jeon 2025 Table 1 and the final-model panel of Figure S1.
    lcl <- log(8.91); label("Clearance of free drug from central (L/h)") # Text S1 $THETA 1 (0, 8.91); Table 1 CL = 8.91 L/h (RSE 18.7%)
    lvc <- log(39); label("Central volume of distribution (L)") # Text S1 $THETA 2 (0, 39); Table 1 V_C = 39 L (RSE 20.2%)
    lq <- log(5.46); label("Inter-compartmental clearance between injection-site space and central (L/h)") # Text S1 $THETA 3 (0, 5.46); Table 1 Q = 5.46 L/h (RSE 34.4%)
    lvendosome <- log(0.138); label("Apparent volume of the FcRn-recycling endosomal space V_D (L)") # Text S1 $THETA 5 (0, 0.138); Table 1 V_D = 0.138 L (RSE 66.2%)
    lkrec <- log(0.00107); label("Exit rate constant of the drug-FcRn complex from the endosome, k_recycle (1/h)") # Text S1 $THETA 6 (0, 0.00107); Table 1 K_rec = 0.00107 1/h (RSE 21.6%)
    lkint <- log(0.375); label("Internalisation rate constant of the drug-IL-7R complex in central (1/h)") # Text S1 $THETA 7 (0, 0.375), estimated (not FIX); Table 1 K_int = 0.375 1/h (RSE 43.2%)
    lvabs_im <- log(442); label("Apparent injection-site volume V_I after IM injection (L)") # Text S1 $THETA 11 (0, 442) 'VIIM'; Table 1 V_I_IM = 442 L (RSE 29%)
    lvabs_sc <- log(1680); label("Apparent injection-site volume V_I after SC injection (L)") # Text S1 $THETA 12 (0, 1680) 'VISC'; Table 1 V_I_SC = 1680 L (RSE 17.5%)

    # Held fixed by the authors as $THETA ... FIX (provenance in Table S1).
    lkup <- fixed(log(0.00952)); label("Uptake rate constant from the injection-site space into the endosome, k_uptake (1/h)") # Text S1 $THETA 4 '(0, 0.00952) FIX'; Table S1 K_uptake = 0.00952 1/h (Waldmann 1990)
    lc_il7r_t <- fixed(log(1360)); label("Total IL-7R (CD127) concentration in central, R_tot (pmol/L)") # Text S1 $THETA 8 '(0, 1360) FIX'; Table S1 R_tot = 1360 pmol/L (sCD127, Monti 2013 Figure 3E)

    # Constants hard-coded in the Text S1 $PK block (not $THETAs); provenance in
    # Jeon 2025 Methods 2.1 and Table S1. Binding rate constants were measured by
    # Genexine by biolayer interferometry.
    c_fcrn_t <- fixed(91950000); label("Total endosomal FcRn concentration (pmol/L)") # Text S1 $PK 'FcRn = 91950000'; Table S1 FcRn = 91,950,000 pmol/L (Barber 2023)
    fr_recycle <- fixed(0.715); label("Fraction of the FcRn-bound drug recycled to central; the rest returns to the injection site (unitless)") # Text S1 $PK 'FR = 0.715'; Table S1 FR = 0.715 (Zhao 2013, Garg 2007)
    koff_fcrn <- fixed(12.132); label("Drug-FcRn dissociation rate constant (1/h)") # Text S1 $PK 'KOFF1 = 12.132'
    kon_fcrn <- fixed(0.0004068); label("Drug-FcRn association rate constant (L/pmol/h)") # Text S1 $PK 'KON1 = 4.068/10000'
    koff_il7r <- fixed(59.76); label("Drug-IL-7R dissociation rate constant (1/h)") # Text S1 $PK 'KOFF2 = 59.76'
    kon_il7r <- fixed(0.0020196); label("Drug-IL-7R association rate constant (L/pmol/h)") # Text S1 $PK 'KON2 = 20.196/10000'

    # IIV. Text S1 $OMEGA holds log-scale variances: sqrt(exp(0.347) - 1) = 64.4%
    # and sqrt(exp(0.659) - 1) = 96.6%, exactly the CV% in Table 1.
    etalcl ~ 0.347 # Text S1 $OMEGA 1 (ETA(1) on CL); Table 1 IIV CL 64.4% (shrinkage 5%)
    etalq ~ 0.659 # Text S1 $OMEGA 2 (ETA(2) on Q); Table 1 IIV Q 96.6% (shrinkage 3%)

    # Residual error. Text S1 $ERROR: W = SQRT(THETA(9)**2 * IPRED**2 + THETA(10)**2)
    # with $SIGMA 1 FIX, applied to the free central concentration.
    propSd <- 0.335; label("Proportional residual error (fraction)") # Text S1 $THETA 9 (0, 0.335); Table 1 proportional error 33.5% (RSE 6.2%)
    addSd <- 1.05; label("Additive residual error (pmol/L)") # Text S1 $THETA 10 (0, 1.05); Table 1 additive error 1.05 pmol/L (RSE 11.2%)
  })
  model({
    # IL-7 moieties per rhIL-7-hyFc homodimer. The control stream doses pmol of
    # drug into A(1) with F1 = F2 = 2 ('AMT: GX-I7 pmol x 2 = IL7 pmol'); every
    # state and concentration below is in pmol / pmol/L of IL-7.
    n_il7 <- 2

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq + etalq)
    vendosome <- exp(lvendosome)
    krec <- exp(lkrec)
    kint <- exp(lkint)
    kup <- exp(lkup)
    c_il7r_t <- exp(lc_il7r_t)

    # Route-specific injection-site volume (control stream: VI = VISC, IM cohort VIIM)
    vabs <- exp(lvabs_sc) * ROUTE_SC + exp(lvabs_im) * (1 - ROUTE_SC)

    kel <- cl / vc
    k12 <- q / vabs
    k21 <- q / vc

    # QSS constants, Jeon 2025 Eqs. (1) and (2)
    kss_fcrn <- (krec + koff_fcrn) / kon_fcrn
    kss_il7r <- (kint + koff_il7r) / kon_il7r

    # Free concentrations from the totals by the QSS quadratic. The control stream
    # writes free = 0.5 * (D + sqrt(D^2 + 4 * KSS * Ctot)); the form below is the
    # same root rationalised, which avoids the cancellation of two ~1e8 numbers in
    # the endosome (FcRn is 92 uM). D < 0 whenever Ctot = 0, so it never divides by 0.
    a_tot <- endosome / vendosome
    dif_fcrn <- a_tot - c_fcrn_t - kss_fcrn
    a_free <- 2 * kss_fcrn * a_tot / (sqrt(dif_fcrn^2 + 4 * kss_fcrn * a_tot) - dif_fcrn)
    a_cpx <- c_fcrn_t * a_free / (kss_fcrn + a_free)

    c_tot <- central / vc
    dif_il7r <- c_tot - c_il7r_t - kss_il7r
    c_free <- 2 * kss_il7r * c_tot / (sqrt(dif_il7r^2 + 4 * kss_il7r * c_tot) - dif_il7r)
    c_cpx <- c_il7r_t * c_free / (kss_il7r + c_free)

    # Text S1 $DES: A(1) injection site, A(2) central, A(3) FcRn-recycling endosome
    d/dt(abs_site) <- -kup * abs_site + krec * a_cpx * vendosome * (1 - fr_recycle) - k12 * abs_site + k21 * c_free * vc
    d/dt(central) <- k12 * abs_site - k21 * c_free * vc + krec * a_cpx * vendosome * fr_recycle - kel * c_free * vc - kint * c_cpx * vc
    d/dt(endosome) <- kup * abs_site - krec * a_cpx * vendosome

    f(abs_site) <- n_il7
    f(central) <- n_il7

    # Observation: free IL-7 concentration in central (Text S1 $ERROR IPRED)
    Cc <- c_free
    Cc ~ prop(propSd) + add(addSd)
  })
}
