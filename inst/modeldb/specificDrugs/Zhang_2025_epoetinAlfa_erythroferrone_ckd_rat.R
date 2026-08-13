Zhang_2025_epoetinAlfa_erythroferrone_ckd_rat <- function() {
  description <- paste(
    "Preclinical (rat). Mechanism-based PK/PD model of erythroferrone",
    "(ERFE), hemoglobin and red blood cells after intravenous recombinant",
    "human erythropoietin (rHuEPO / epoetin alfa) in adenine-induced",
    "chronic-kidney-disease anemic rats (Zhang 2025, model structure of",
    "Figure 1A). rHuEPO PK is a two-compartment model with parallel linear",
    "and Michaelis-Menten elimination, fixed from prior rat studies (Table",
    "S1). The erythroid response is a lifespan-based indirect-response",
    "chain: a zero-order progenitor input kIN0, linearly stimulated by",
    "rHuEPO concentration and damped by an (RBC0 / RBC)^gamma feedback,",
    "feeds BFU-E (precursor1) which amplifies AMPN-fold into CFU-E",
    "(precursor2), amplifies AMPC-fold into erythroblasts (precursor3),",
    "then reticulocytes (precursor4) and circulating erythrocytes;",
    "hemoglobin is derived algebraically from the RBC count and the mean",
    "corpuscular hemoglobin. Circulating ERFE receives two simultaneous",
    "inputs: a circadian baseline pool driven by a cosinor production rate",
    "kCINE (mesor RM, amplitude RA, acrophase tpeak), and an",
    "rHuEPO-stimulated pool whose production is proportional to the",
    "erythroblast mass (precursor3) and to an Emax / EC50 function of",
    "rHuEPO concentration, delayed by two transit compartments. Fitted in",
    "NONMEM 7.5 with interindividual variability fixed to zero, so the",
    "model is a typical-value mechanism with proportional residual error",
    "on each of the three observed outputs."
  )
  reference <- paste(
    "Zhang L, Xu P, Yan X.",
    "Mechanism-Based Pharmacokinetic/Pharmacodynamic Modeling of",
    "Erythroferrone in Anemic Rats with Chronic Kidney Disease and",
    "Chemotherapy-Induced Anemia: An Early Biomarker for Hemoglobin",
    "Response and rHuEPO Hyporesponsiveness.",
    "ACS Pharmacol Transl Sci. 2025;8(1):189-202.",
    "doi:10.1021/acsptsci.4c00575.",
    "rHuEPO PK parameters fixed from Fan X, Krzyzanski W, Wong RSM, Yan X.",
    "J Pharmacol Exp Ther. 2022;382(1):31-43, as tabulated in",
    "Zhang 2025 Supporting Information Table S1."
  )
  vignette <- "Zhang_2025_erythroferrone"

  units <- list(
    time                      = "h",
    dosing                    = "mIU/kg",
    concentration             = "mIU/mL",
    erythroferrone            = "ng/mL",
    hemoglobin                = "g/dL",
    red_blood_cells           = "1e12 cells/L"
  )

  # The three erythroferrone states are paper-mechanistic decompositions of
  # the measured circulating ERFE concentration (Zhang 2025 eqs 18-22 and
  # Figure 1A): erfe_base is the latent circadian baseline pool, erfe_induced
  # is the rHuEPO-driven marrow release pool feeding transit1 / transit2, and
  # erfe is the measured circulating concentration receiving both inputs.
  # They are declared paper-specific rather than registered as new canonical
  # compartments; see the vignette "Assumptions and deviations" section for
  # the naming proposal left for the register maintainer.
  paper_specific_compartments <- c("erfe", "erfe_base", "erfe_induced")

  compartmentData <- list(
    central      = list(analyte = "epoetin alfa",   units = "mIU/kg",       specimen = "plasma",         verified = TRUE),
    peripheral1  = list(analyte = "epoetin alfa",   units = "mIU/kg",       specimen = "plasma",         verified = TRUE),
    precursor1   = list(analyte = "BFU-E cells",    units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor2   = list(analyte = "CFU-E cells",    units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor3   = list(analyte = "erythroblasts",  units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor4   = list(analyte = "reticulocytes",  units = "1e12 cells/L", specimen = "blood cell",     verified = TRUE),
    erythrocytes = list(analyte = "red blood cells", units = "1e12 cells/L", specimen = "blood cell",    verified = TRUE),
    erfe_base    = list(analyte = "erythroferrone", units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    erfe_induced = list(analyte = "erythroferrone", units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    transit1     = list(analyte = "erythroferrone", units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    transit2     = list(analyte = "erythroferrone", units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    erfe         = list(analyte = "erythroferrone", units = "ng/mL",        specimen = "plasma",         verified = TRUE)
  )

  population <- list(
    species        = "rat (Sprague-Dawley, male)",
    n_subjects     = 17L,
    n_studies      = 1L,
    weight_range   = "250-300 g at study entry",
    sex_female_pct = 0,
    disease_state  = paste(
      "Adenine-induced chronic-kidney-disease anemia. Adenine 600 mg/kg",
      "once daily by gavage for 6 days, then 300 mg/kg once daily for 6",
      "days, then a 1-week stabilisation, on a standard diet supplemented",
      "with 1% (w/w) carbonyl iron. Anemia confirmed as a fall of more",
      "than 2 g/dL in hemoglobin versus healthy controls and persisting",
      "for at least 1 month (Zhang 2025 Figure S3A,B). Model baseline",
      "hemoglobin MCH * RBC0 / 10 = 12.4 g/dL."
    ),
    dose_range     = paste(
      "Intravenous rHuEPO (EPOGEN, epoetin alfa) 1350 IU/kg (n = 6) or",
      "450 IU/kg (n = 6) three times weekly for 2 weeks, or saline",
      "(n = 5). Doses enter the model in mIU/kg (1350 IU/kg = 1.35e6",
      "mIU/kg)."
    ),
    regions        = "Hong Kong SAR, China (The Chinese University of Hong Kong)",
    notes          = paste(
      "ERFE assayed at 0, 1, 2, 4, 6, 8, 10, 12 and 24 h after the first",
      "injection by validated ELISA (FineTest ER1573); hemoglobin and RBC",
      "on 'Week 0' days 0, 4, 10, 15, 20, 25, 30 and 34 on a BC2800VET",
      "analyser. A rotation sampling method limited blood loss. Fitted in",
      "NONMEM 7.5 by nonlinear regression; interindividual variability was",
      "fixed to zero (Zhang 2025 Methods, 'PK/PD Modeling'), so the",
      "packaged model is a typical-value mechanism. Parameter estimates",
      "from Zhang 2025 Table 1; rHuEPO PK from Table S1 (fixed)."
    )
  )

  ini({
    # =================================================================
    # rHuEPO PK -- ALL FIXED from prior rat studies, tabulated in
    # Zhang 2025 Supporting Information Table S1 ("PK parameters were
    # fixed based on previous studies"). Amounts are per kg body weight
    # because the volume is reported per kg.
    # =================================================================
    lvmax <- fixed(log(1993));   label("rHuEPO Michaelis-Menten capacity Vmax (mIU/h/kg)")            # Table S1: Vmax 1993 mIU h-1 kg-1
    lkm   <- fixed(log(67.28));  label("rHuEPO Michaelis-Menten affinity Km (mIU/mL)")                # Table S1: Km 67.28 mIU mL-1
    lvc   <- fixed(log(61.18));  label("rHuEPO central volume of distribution VEPO (mL/kg)")           # Table S1: VEPO 61.18 mL kg-1
    lkel  <- fixed(log(0.209));  label("rHuEPO linear elimination rate constant kel,EPO (1/h)")        # Table S1: Kel,EPO 0.209 h-1
    lk12  <- fixed(log(0.171));  label("rHuEPO central-to-peripheral rate constant kpt,EPO (1/h)")     # Table S1: Kpt,EPO 0.171 h-1
    lk21  <- fixed(log(0.148));  label("rHuEPO peripheral-to-central rate constant ktp,EPO (1/h)")     # Table S1: Ktp,EPO 0.148 h-1

    # =================================================================
    # ERYTHROID RESPONSE -- estimated (Zhang 2025 Table 1)
    # =================================================================
    ltrbc      <- log(1280);    label("Mean lifespan of circulating red blood cells TRBC (h)")         # Table 1: T RBC 1280 h (RSE 13.2%)
    ltp        <- log(20.7);    label("Mean lifespan of erythroid precursor cells T (h)")              # Table 1: T 20.7 h (RSE 8.5%)
    lrbase_rbc <- log(6.15);    label("Baseline red blood cell count RBC0 (1e12 cells/L)")             # Table 1: RBC0 6.15 (RSE 0.9%)
    lmch       <- log(20.2);    label("Mean corpuscular hemoglobin MCH (pg/cell)")                     # Table 1: MCH 20.2 pg cell-1 (RSE 0.9%)
    lgamma     <- log(3.9);     label("RBC feedback exponent gamma on (RBC0 / RBC) (unitless)")        # Table 1: gamma 3.9 (RSE 19.9%)
    lslope     <- log(0.0106);  label("Linear rHuEPO stimulation of progenitor input KEPO (mL/mIU)")   # Table 1: K EPO 0.0106 (RSE 19.2%)

    # Structural amplification factors: the number of CFU-E produced by one
    # BFU-E (AMPN) and the number of erythroblasts produced by one CFU-E
    # (AMPC). Both assumed equal to 2^5 = 32 by the authors, i.e. fixed.
    ampn <- fixed(32);          label("BFU-E to CFU-E amplification factor AMPN (cells/cell)")         # Methods after eq 11: AMPN = AMPC = 2^5
    ampc <- fixed(32);          label("CFU-E to erythroblast amplification factor AMPC (cells/cell)")  # Methods after eq 11: AMPN = AMPC = 2^5

    # =================================================================
    # ERYTHROFERRONE -- estimated (Zhang 2025 Table 1)
    # =================================================================
    lkout  <- log(2.82);   label("First-order elimination rate constant of circadian ERFE Kout (1/h)")     # Table 1: K out 2.82 h-1 (RSE 18.8%)
    lrm    <- log(7.25);   label("Mesor (mean baseline) of the ERFE circadian rhythm RM (ng/mL)")          # Table 1: RM 7.25 ng mL-1 (RSE 8.6%)
    lra    <- log(2.51);   label("Amplitude of the ERFE circadian rhythm RA (ng/mL)")                      # Table 1: RA 2.51 ng mL-1 (RSE 15.9%)
    ltacro <- log(14.9);   label("Acrophase (peak time) of the ERFE circadian rhythm tpeak (h)")           # Table 1: t peak 14.9 h (RSE 4%)
    lemax  <- log(47.2);   label("Maximum ERFE induction by rHuEPO per erythroblast Emax ((ng/mL)/(1e12 cells/L))")  # Table 1: E max 47.2 (RSE 37.5%)
    lec50  <- log(1360);   label("rHuEPO concentration giving 50% of maximum ERFE induction EC50 (mIU/mL)")  # Table 1: EC 50 1360 mIU mL-1 (RSE 69.9%)
    lktr   <- log(2.25);   label("Transit rate constant of the ERFE release chain Ktr (1/h)")              # Table 1: K tr 2.25 h-1 (RSE 20.5%)

    # =================================================================
    # RESIDUAL ERROR -- proportional on each observed output. Zhang 2025
    # Table 1 reports sigma for RBC / HGB / ERFE as 0.495 / 1.03 / 1.54
    # in a column headed "sigma (%)". Those numbers are the proportional
    # error VARIANCE expressed as a percentage, so the coefficient of
    # variation is sqrt(value / 100): 7.0% / 10.1% / 12.4%. This reading
    # is the only one of the three candidates that reproduces the width of
    # the published pcVPC prediction bands in Zhang 2025 Figure 6A-C
    # (taking the values as CV% directly gives bands ~20x too narrow;
    # taking them as unitless variances gives CV 70-124%, far too wide).
    # The discriminator is documented in the vignette Errata.
    # =================================================================
    propSd_erythrocytes <- 0.07036;  label("Proportional residual error on RBC (fraction; CV 7.0%)")    # Table 1: sigma RBC 0.495 (RSE 6.2%)
    propSd_hb           <- 0.1015;   label("Proportional residual error on HGB (fraction; CV 10.1%)")   # Table 1: sigma HGB 1.03 (RSE 6.2%)
    propSd_erfe         <- 0.1241;   label("Proportional residual error on ERFE (fraction; CV 12.4%)")  # Table 1: sigma ERFE 1.54 (RSE 8.1%)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Individual (typical-value) parameters. Interindividual
    #    variability was fixed to zero by the authors, so no etas.
    # ---------------------------------------------------------------
    vmax  <- exp(lvmax)
    km    <- exp(lkm)
    vc    <- exp(lvc)
    kel   <- exp(lkel)
    k12   <- exp(lk12)
    k21   <- exp(lk21)

    trbc  <- exp(ltrbc)
    tp    <- exp(ltp)
    rbc0  <- exp(lrbase_rbc)
    mch   <- exp(lmch)
    gamma <- exp(lgamma)
    slope <- exp(lslope)

    kout  <- exp(lkout)
    rm    <- exp(lrm)
    ra    <- exp(lra)
    tacro <- exp(ltacro)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)
    ktr   <- exp(lktr)

    # ---------------------------------------------------------------
    # 2. Derived structural quantities
    # ---------------------------------------------------------------
    # Progenitor input rate constant (Zhang 2025 Methods after eq 11):
    #   kIN0 = RBC0 / ((TRBC + TRET) * AMPN * AMPC)
    # TRET is never reported anywhere in Zhang 2025 or its Supporting
    # Information. In this model the reticulocyte pool is precursor4 and
    # its residence time is the precursor lifespan T, so TRET is taken
    # equal to T. That choice makes kIN0 exactly equal to precursor1(0) / T,
    # i.e. it closes the progenitor balance at the paper's own printed
    # initial conditions. Documented in the vignette Errata.
    kin0 <- rbc0 / ((trbc + tp) * ampn * ampc)

    # Baseline erythroblast pool. The paper writes precursor3(0) = RET0 but
    # never reports RET0; the zero-net-flux steady state of eq 8,
    # AMPC * precursor2(0) / T = precursor3(0) / T, fixes it. Documented in
    # the vignette Errata.
    ret0 <- rbc0 * tp / (trbc + tp)

    # Circadian angular frequency for a 24 h period (1/h).
    wcirc <- 2 * pi / 24

    # rHuEPO plasma concentration (mIU/mL).
    Cc <- central / vc

    # ---------------------------------------------------------------
    # 3. rHuEPO PK: two compartments, parallel linear and
    #    Michaelis-Menten elimination (Zhang 2025 eqs 1-2). Amounts are
    #    per kg, so km * vc converts Km to an amount (mIU/kg).
    # ---------------------------------------------------------------
    d/dt(central)     <- -vmax * central / (km * vc + central) -
                          (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---------------------------------------------------------------
    # 4. Erythroid response: lifespan-based indirect-response chain
    #    (Zhang 2025 eqs 6-11). Zero-order progenitor input kIN0 is
    #    stimulated linearly by rHuEPO and damped by the RBC feedback.
    # ---------------------------------------------------------------
    d/dt(precursor1)   <- kin0 * (1 + slope * Cc) *
                          (rbc0 / erythrocytes)^gamma - precursor1 / tp
    d/dt(precursor2)   <- ampn * precursor1 / tp - precursor2 / tp
    d/dt(precursor3)   <- ampc * precursor2 / tp - precursor3 / tp
    d/dt(precursor4)   <- precursor3 / tp - precursor4 / tp
    d/dt(erythrocytes) <- precursor4 / tp - erythrocytes / trbc

    # Hemoglobin (g/dL) from RBC count (1e12 cells/L) and MCH (pg/cell).
    # Zhang 2025 eq 11 is typeset as dAHGB/dt = MCH * ARBC / 10, which is
    # dimensionally a rate and would grow without bound; Figure 1A draws
    # HGB as derived from RBC, and MCH * RBC0 / 10 = 12.4 g/dL reproduces
    # the observed CKD baseline of Figure 2B exactly. Encoded as the
    # algebraic relation the paper intends; noted in the vignette Errata.
    hb <- mch * erythrocytes / 10

    # ---------------------------------------------------------------
    # 5. Circadian ERFE baseline pool (Zhang 2025 eq 18). kCINE is the
    #    time-dependent production rate that makes erfe_base track the
    #    cosinor RM + RA * cos(2*pi/24 * (t - tpeak)) exactly:
    #      kCINE = kout*RM + kout*RA*cos(theta) - (2*pi/24)*RA*sin(theta)
    #    The paper writes the second coefficient as kout,ERFE; only one
    #    kout is reported, and setting kout,ERFE = kout is what makes the
    #    cosinor an exact solution of eq 18.
    # ---------------------------------------------------------------
    kcine <- kout * rm + kout * ra * cos(wcirc * (t - tacro)) -
             wcirc * ra * sin(wcirc * (t - tacro))
    d/dt(erfe_base) <- kcine - kout * erfe_base

    # ---------------------------------------------------------------
    # 6. rHuEPO-induced ERFE release (Zhang 2025 eqs 19-21). Production
    #    is proportional to the erythroblast mass precursor3 and to an
    #    Emax / EC50 function of rHuEPO concentration; the marrow-to-blood
    #    delay is two transit compartments. The paper sets the ERFE
    #    production and elimination rate constants kin and kel equal to
    #    ktr to limit the number of estimated parameters.
    # ---------------------------------------------------------------
    d/dt(erfe_induced) <- ktr * precursor3 * emax * Cc / (ec50 + Cc) -
                          ktr * erfe_induced
    d/dt(transit1)     <- ktr * erfe_induced - ktr * transit1
    d/dt(transit2)     <- ktr * transit1     - ktr * transit2

    # Circulating (measured) ERFE: circadian input plus the delayed
    # rHuEPO-induced input (Zhang 2025 eq 22, with kel = ktr).
    d/dt(erfe) <- kout * erfe_base + ktr * transit2 - ktr * erfe

    # ---------------------------------------------------------------
    # 7. Initial conditions (Zhang 2025 Methods, inline expressions
    #    following eq 11, and following eq 22 for the ERFE states).
    #    Total circulating cells partition into reticulocytes and mature
    #    erythrocytes, so precursor4(0) + erythrocytes(0) = RBC0.
    # ---------------------------------------------------------------
    precursor1(0)   <- rbc0 * tp / ((trbc + tp) * ampc * ampn)
    precursor2(0)   <- rbc0 * tp / ((trbc + tp) * ampc)
    precursor3(0)   <- ret0
    precursor4(0)   <- rbc0 * tp / (trbc + tp)
    erythrocytes(0) <- rbc0 - rbc0 * tp / (trbc + tp)

    erfe_base(0)    <- rm + ra * cos(wcirc * (0 - tacro))
    erfe(0)         <- rm + ra * cos(wcirc * (0 - tacro))

    # ---------------------------------------------------------------
    # 8. Observations and residual error (Zhang 2025 Methods,
    #    "Interindividual variability (IIV) was fixed as zero. The
    #    residual variance was described by a proportional error model").
    # ---------------------------------------------------------------
    erfe         ~ prop(propSd_erfe)
    hb           ~ prop(propSd_hb)
    erythrocytes ~ prop(propSd_erythrocytes)
  })
}
