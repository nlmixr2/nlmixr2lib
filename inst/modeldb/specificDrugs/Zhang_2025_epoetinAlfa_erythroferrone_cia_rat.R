Zhang_2025_epoetinAlfa_erythroferrone_cia_rat <- function() {
  description <- paste(
    "Preclinical (rat). Mechanism-based PK/PD model of erythroferrone",
    "(ERFE), hemoglobin and red blood cells after intravenous recombinant",
    "human erythropoietin (rHuEPO / epoetin alfa) in carboplatin-induced",
    "chemotherapy-induced-anemia rats (Zhang 2025, model structure of",
    "Figure 1B). Two co-modelled drugs, both fixed from prior rat studies",
    "(Table S1): rHuEPO with two-compartment disposition and parallel",
    "linear and Michaelis-Menten elimination, and carboplatin with",
    "three-compartment linear disposition. The erythroid response is a",
    "Friberg-Karlsson myelosuppression chain: a self-renewing progenitor",
    "pool (precursor1) proliferating at kprol, inhibited linearly by",
    "carboplatin concentration and damped by an (RBC0 / RBC)^gamma",
    "feedback, feeding three transit compartments (precursor2 to",
    "precursor4) and a circulating erythrocyte pool, with kprol = kcirc =",
    "ktr = (N + 1) / MTT for N = 3 transits; hemoglobin is derived",
    "algebraically from the RBC count and the mean corpuscular hemoglobin.",
    "Circulating ERFE receives two simultaneous inputs: a circadian",
    "baseline pool driven by a cosinor production rate kCINE (mesor RM,",
    "amplitude RA, acrophase tpeak), and an rHuEPO-stimulated pool whose",
    "production is proportional to the erythroblast mass (precursor3) and",
    "to an Emax / EC50 function of rHuEPO concentration, delayed by two",
    "transit compartments. Because carboplatin drives precursor3 to a",
    "nadir by the time rHuEPO is given, the model reproduces the blunted",
    "ERFE induction that marks rHuEPO hyporesponsiveness. Fitted in NONMEM",
    "7.5 with interindividual variability fixed to zero, so the model is a",
    "typical-value mechanism with proportional residual error on each of",
    "the three observed outputs."
  )
  reference <- paste(
    "Zhang L, Xu P, Yan X.",
    "Mechanism-Based Pharmacokinetic/Pharmacodynamic Modeling of",
    "Erythroferrone in Anemic Rats with Chronic Kidney Disease and",
    "Chemotherapy-Induced Anemia: An Early Biomarker for Hemoglobin",
    "Response and rHuEPO Hyporesponsiveness.",
    "ACS Pharmacol Transl Sci. 2025;8(1):189-202.",
    "doi:10.1021/acsptsci.4c00575.",
    "rHuEPO and carboplatin PK parameters fixed from Fan X, Krzyzanski W,",
    "Wong RSM, Liu D, Yan X. ACS Pharmacol Transl Sci. 2023;6:1884-1897,",
    "doi:10.1021/acsptsci.3c00194, as tabulated in Zhang 2025 Supporting",
    "Information Table S1."
  )
  vignette <- "Zhang_2025_erythroferrone"

  units <- list(
    time                      = "h",
    dosing                    = "mIU/kg",
    concentration             = "mIU/mL",
    dosing_carboplatin        = "ug/kg",
    concentration_carboplatin = "ug/mL",
    erythroferrone            = "ng/mL",
    hemoglobin                = "g/dL",
    red_blood_cells           = "1e12 cells/L"
  )

  # Two groups of paper-specific states. (a) The three erythroferrone states
  # are paper-mechanistic decompositions of the measured circulating ERFE
  # concentration (Zhang 2025 eqs 18-22 and Figure 1B): erfe_base is the
  # latent circadian baseline pool, erfe_induced is the rHuEPO-driven marrow
  # release pool feeding transit1 / transit2, and erfe is the measured
  # circulating concentration receiving both inputs. (b) The carboplatin
  # disposition states carry a `_carb` sibling-drug suffix; carboplatin is
  # co-administered with rHuEPO and is not a metabolite of it, and `carb` is
  # not yet a registered sibling-drug suffix. Both groups are declared
  # paper-specific rather than registered as new canonicals; see the
  # vignette "Assumptions and deviations" section for the naming proposal
  # left for the register maintainer.
  paper_specific_compartments <- c(
    "erfe", "erfe_base", "erfe_induced",
    "central_carb", "peripheral1_carb", "peripheral2_carb"
  )

  compartmentData <- list(
    central          = list(analyte = "epoetin alfa",    units = "mIU/kg",       specimen = "plasma",         verified = TRUE),
    peripheral1      = list(analyte = "epoetin alfa",    units = "mIU/kg",       specimen = "plasma",         verified = TRUE),
    central_carb     = list(analyte = "carboplatin",     units = "ug/kg",        specimen = "plasma",         verified = TRUE),
    peripheral1_carb = list(analyte = "carboplatin",     units = "ug/kg",        specimen = "plasma",         verified = TRUE),
    peripheral2_carb = list(analyte = "carboplatin",     units = "ug/kg",        specimen = "plasma",         verified = TRUE),
    precursor1       = list(analyte = "BFU-E cells",     units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor2       = list(analyte = "CFU-E cells",     units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor3       = list(analyte = "erythroblasts",   units = "1e12 cells/L", specimen = "not applicable",    verified = TRUE),
    precursor4       = list(analyte = "reticulocytes",   units = "1e12 cells/L", specimen = "blood cell",     verified = TRUE),
    erythrocytes     = list(analyte = "red blood cells", units = "1e12 cells/L", specimen = "blood cell",     verified = TRUE),
    erfe_base        = list(analyte = "erythroferrone",  units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    erfe_induced     = list(analyte = "erythroferrone",  units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    transit1         = list(analyte = "erythroferrone",  units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    transit2         = list(analyte = "erythroferrone",  units = "ng/mL",        specimen = "not applicable", verified = TRUE),
    erfe             = list(analyte = "erythroferrone",  units = "ng/mL",        specimen = "plasma",         verified = TRUE)
  )

  population <- list(
    species        = "rat (Sprague-Dawley, male)",
    n_subjects     = 27L,
    n_studies      = 1L,
    weight_range   = "250-300 g at study entry",
    sex_female_pct = 0,
    disease_state  = paste(
      "Carboplatin-induced (chemotherapy-induced) anemia. A single",
      "intravenous carboplatin dose of 60 mg/kg via the tail vein",
      "produced a rapid fall in RBC and hemoglobin to a nadir of about",
      "8 g/dL on day 14, with recovery to normal by about 1 month",
      "(Zhang 2025 Figure S3C,D). Flow cytometry confirmed near-complete",
      "depletion of marrow erythroid precursors with less effect on",
      "splenic precursors (Figure 3C). Model baseline hemoglobin",
      "MCH * RBC0 / 10 = 13.6 g/dL."
    ),
    dose_range     = paste(
      "Single intravenous carboplatin 60 mg/kg at study time zero, then a",
      "single intravenous rHuEPO (EPOGEN, epoetin alfa) dose of",
      "1350 IU/kg (n = 9) or 450 IU/kg (n = 9) or saline (n = 9) one week",
      "later ('Week 0'). Doses enter the model as 6e4 ug/kg carboplatin",
      "and 1.35e6 or 4.5e5 mIU/kg rHuEPO."
    ),
    regions        = "Hong Kong SAR, China (The Chinese University of Hong Kong)",
    notes          = paste(
      "ERFE assayed at 0, 1, 2, 4, 6, 8, 10, 12 and 24 h after the rHuEPO",
      "or saline injection by validated ELISA (FineTest ER1573);",
      "hemoglobin and RBC on 'Week -1' days -8, -6, -4, -2, 0, 3, 6, 9,",
      "12, 13, 14, 16, 17, 18, 20, 21 and 24 on a BC2800VET analyser. A",
      "rotation sampling method limited blood loss. Additional rats",
      "(n = 3 per group per day) were sacrificed on days 0, 1, 4, 8, 11,",
      "14 and 15 after carboplatin for bone-marrow and spleen flow",
      "cytometry. Fitted in NONMEM 7.5 by nonlinear regression;",
      "interindividual variability was fixed to zero (Zhang 2025 Methods,",
      "'PK/PD Modeling'), so the packaged model is a typical-value",
      "mechanism. Parameter estimates from Zhang 2025 Table 2; rHuEPO and",
      "carboplatin PK from Table S1 (fixed). The hematological response to",
      "rHuEPO showed high between-animal variability in this cohort, with",
      "some rats hyporesponsive; that variability is not carried in the",
      "model because IIV was fixed to zero."
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
    # Carboplatin PK -- ALL FIXED (Table S1). Three-compartment linear
    # disposition; amounts per kg, concentration in ug/mL (= mg/L).
    # =================================================================
    lvc_carb  <- fixed(log(148.4));  label("Carboplatin central volume of distribution Vcarb (mL/kg)")            # Table S1: Vcarb 148.4 mL kg-1
    lkel_carb <- fixed(log(2.788));  label("Carboplatin linear elimination rate constant kel,carb (1/h)")          # Table S1: Kel,carb 2.788 h-1
    lk12_carb <- fixed(log(0.074));  label("Carboplatin central-to-first-peripheral rate constant k12 (1/h)")      # Table S1: K12 0.074 h-1
    lk21_carb <- fixed(log(0.425));  label("Carboplatin first-peripheral-to-central rate constant k21 (1/h)")      # Table S1: K21 0.425 h-1
    lk13_carb <- fixed(log(4.194));  label("Carboplatin central-to-second-peripheral rate constant k13 (1/h)")     # Table S1: K13 4.194 h-1
    lk31_carb <- fixed(log(5.611));  label("Carboplatin second-peripheral-to-central rate constant k31 (1/h)")     # Table S1: K31 5.611 h-1

    # =================================================================
    # ERYTHROID RESPONSE -- estimated (Zhang 2025 Table 2)
    # =================================================================
    lmtt       <- log(200);     label("Mean transit time of erythroid precursor maturation MTT (h)")       # Table 2: MTT 200 h (RSE 3.1%)
    lrbase_rbc <- log(6.46);    label("Baseline red blood cell count RBC0 (1e12 cells/L)")                 # Table 2: RBC0 6.46 (RSE 1.1%)
    lmch       <- log(21.1);    label("Mean corpuscular hemoglobin MCH (pg/cell)")                          # Table 2: MCH 21.1 pg cell-1 (RSE 1.1%)
    lgamma     <- log(0.189);   label("RBC feedback exponent gamma on (RBC0 / RBC) (unitless)")             # Table 2: gamma 0.189 (RSE 8%)
    lslope     <- log(0.172);   label("Linear carboplatin killing effect Kcarb on the precursor pool (mL/ug)")  # Table 2: K carb 0.172 (RSE 4.6%)

    # =================================================================
    # ERYTHROFERRONE -- estimated (Zhang 2025 Table 2)
    # =================================================================
    lkout  <- log(1.16);   label("First-order elimination rate constant of circadian ERFE Kout (1/h)")     # Table 2: K out 1.16 h-1 (RSE 19.9%)
    lrm    <- log(9.03);   label("Mesor (mean baseline) of the ERFE circadian rhythm RM (ng/mL)")          # Table 2: RM 9.03 ng mL-1 (RSE 18.2%)
    lra    <- log(3.69);   label("Amplitude of the ERFE circadian rhythm RA (ng/mL)")                      # Table 2: RA 3.69 ng mL-1 (RSE 16.5%)
    ltacro <- log(15.1);   label("Acrophase (peak time) of the ERFE circadian rhythm tpeak (h)")           # Table 2: t peak 15.1 h (RSE 7.4%)
    lemax  <- log(1.16);   label("Maximum ERFE induction by rHuEPO per erythroblast Emax ((ng/mL)/(1e12 cells/L))")  # Table 2: E max 1.16 (RSE 14.5%)
    lec50  <- log(126);    label("rHuEPO concentration giving 50% of maximum ERFE induction EC50 (mIU/mL)")  # Table 2: EC 50 126 mIU mL-1 (RSE 31.2%)
    lktr   <- log(2.26);   label("Transit rate constant of the ERFE release chain Ktr (1/h)")              # Table 2: K tr 2.26 h-1 (RSE 14.9%)

    # =================================================================
    # RESIDUAL ERROR -- proportional on each observed output. Zhang 2025
    # Table 2 reports sigma for RBC / HGB / ERFE as 0.748 / 1.51 / 2.2 in
    # a column headed "sigma (%)". Those numbers are the proportional
    # error VARIANCE expressed as a percentage, so the coefficient of
    # variation is sqrt(value / 100): 8.6% / 12.3% / 14.8%. This reading
    # is the only one of the three candidates that reproduces the width of
    # the published pcVPC prediction bands in Zhang 2025 Figure 6D-F
    # (taking the values as CV% directly gives bands ~20x too narrow;
    # taking them as unitless variances gives CV 86-148%, far too wide).
    # The discriminator is documented in the vignette Errata.
    # =================================================================
    propSd_erythrocytes <- 0.08649;  label("Proportional residual error on RBC (fraction; CV 8.6%)")    # Table 2: sigma RBC 0.748 (RSE 4.2%)
    propSd_hb           <- 0.1229;   label("Proportional residual error on HGB (fraction; CV 12.3%)")   # Table 2: sigma HGB 1.51 (RSE 4.2%)
    propSd_erfe         <- 0.1483;   label("Proportional residual error on ERFE (fraction; CV 14.8%)")  # Table 2: sigma ERFE 2.2 (RSE 7.9%)
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

    vc_carb  <- exp(lvc_carb)
    kel_carb <- exp(lkel_carb)
    k12_carb <- exp(lk12_carb)
    k21_carb <- exp(lk21_carb)
    k13_carb <- exp(lk13_carb)
    k31_carb <- exp(lk31_carb)

    mtt   <- exp(lmtt)
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
    # Zhang 2025 Methods after eq 17: "the time delay of cell maturation
    # was described by three transit compartments (ktr) and it is assumed
    # that kprol = kcirc = ktr = (N + 1) / MTT, where N is the number of
    # transit compartments". The chain precursor1 -> precursor2 ->
    # precursor3 -> precursor4 -> erythrocytes has N = 3 transits
    # (precursor2 to precursor4), so the shared rate is 4 / MTT = 0.02 1/h.
    # The erythroid chain rate is named ktr_rbc here because the paper
    # reuses the symbol ktr for the numerically different ERFE transit
    # rate constant (Table 2, Ktr = 2.26 1/h).
    ntr     <- 3
    ktr_rbc <- (ntr + 1) / mtt
    kprol   <- ktr_rbc
    kcirc   <- ktr_rbc

    # Circadian angular frequency for a 24 h period (1/h).
    wcirc <- 2 * pi / 24

    # rHuEPO plasma concentration (mIU/mL) and carboplatin plasma
    # concentration (ug/mL = mg/L).
    Cc      <- central / vc
    Cc_carb <- central_carb / vc_carb

    # ---------------------------------------------------------------
    # 3. rHuEPO PK: two compartments, parallel linear and
    #    Michaelis-Menten elimination (Zhang 2025 eqs 1-2). Amounts are
    #    per kg, so km * vc converts Km to an amount (mIU/kg).
    # ---------------------------------------------------------------
    d/dt(central)     <- -vmax * central / (km * vc + central) -
                          (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---------------------------------------------------------------
    # 4. Carboplatin PK: three compartments, linear elimination
    #    (Zhang 2025 eqs 3-5).
    # ---------------------------------------------------------------
    d/dt(central_carb)     <- -(kel_carb + k12_carb + k13_carb) * central_carb +
                               k21_carb * peripheral1_carb +
                               k31_carb * peripheral2_carb
    d/dt(peripheral1_carb) <-  k12_carb * central_carb - k21_carb * peripheral1_carb
    d/dt(peripheral2_carb) <-  k13_carb * central_carb - k31_carb * peripheral2_carb

    # ---------------------------------------------------------------
    # 5. Erythroid response: Friberg-Karlsson myelosuppression chain
    #    (Zhang 2025 eqs 12-17). The self-renewing progenitor pool is
    #    inhibited linearly by carboplatin and damped by the RBC
    #    feedback; rHuEPO has no effect on proliferation in this model
    #    (compare Figure 1B with Figure 1A).
    # ---------------------------------------------------------------
    d/dt(precursor1)   <- kprol * precursor1 * (1 - slope * Cc_carb) *
                          (rbc0 / erythrocytes)^gamma - ktr_rbc * precursor1
    d/dt(precursor2)   <- ktr_rbc * (precursor1 - precursor2)
    d/dt(precursor3)   <- ktr_rbc * (precursor2 - precursor3)
    d/dt(precursor4)   <- ktr_rbc * (precursor3 - precursor4)
    d/dt(erythrocytes) <- ktr_rbc * precursor4 - kcirc * erythrocytes

    # Hemoglobin (g/dL) from RBC count (1e12 cells/L) and MCH (pg/cell).
    # Zhang 2025 eq 17 is typeset as dAHGB/dt = MCH * ARBC / 10, which is
    # dimensionally a rate and would grow without bound; Figure 1B draws
    # HGB as derived from RBC, and MCH * RBC0 / 10 = 13.6 g/dL reproduces
    # the observed CIA baseline of Figure 2E. Encoded as the algebraic
    # relation the paper intends; noted in the vignette Errata.
    hb <- mch * erythrocytes / 10

    # ---------------------------------------------------------------
    # 6. Circadian ERFE baseline pool (Zhang 2025 eq 18). kCINE is the
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
    # 7. rHuEPO-induced ERFE release (Zhang 2025 eqs 19-21). Production
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
    # 8. Initial conditions. Zhang 2025 Methods after eq 17: "The initial
    #    conditions for AP1, AP2, AP3, AP4, ARBC are defined as RBC0."
    #    With kprol = kcirc = ktr_rbc that puts the whole chain at an
    #    exact steady state at t = 0. ERFE initial conditions from the
    #    text following eq 22.
    # ---------------------------------------------------------------
    precursor1(0)   <- rbc0
    precursor2(0)   <- rbc0
    precursor3(0)   <- rbc0
    precursor4(0)   <- rbc0
    erythrocytes(0) <- rbc0

    erfe_base(0)    <- rm + ra * cos(wcirc * (0 - tacro))
    erfe(0)         <- rm + ra * cos(wcirc * (0 - tacro))

    # ---------------------------------------------------------------
    # 9. Observations and residual error (Zhang 2025 Methods,
    #    "Interindividual variability (IIV) was fixed as zero. The
    #    residual variance was described by a proportional error model").
    # ---------------------------------------------------------------
    erfe         ~ prop(propSd_erfe)
    hb           ~ prop(propSd_hb)
    erythrocytes ~ prop(propSd_erythrocytes)
  })
}
