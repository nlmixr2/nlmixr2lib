Zhang_2024_nedosiran <- function() {
  description <- "Two-compartment population PK/PD model for nedosiran (GalNAc-conjugated siRNA targeting hepatic LDHA) in healthy volunteers and patients with primary hyperoxaluria, with dual parallel n-transit subcutaneous absorption (a 2-compartment slow chain and a 4-compartment fast chain), parallel linear and Michaelis-Menten elimination, and an effect-compartment-driven sigmoid-Imax indirect-response model inhibiting production of 24-h urinary oxalate; developed from 1978 plasma concentrations in 143 participants and 588 24-h urinary oxalate observations in 46 PH1 patients across five clinical studies (Zhang 2024)."
  reference   <- paste(
    "Zhang S, Amrite A, Tan B, Jamsen K, Pradhan S, Choy S, Plotkin H.",
    "Nedosiran population pharmacokinetic and pharmacodynamic modelling and",
    "simulation to guide clinical development and dose selection in patients",
    "with primary hyperoxaluria type 1. Br J Clin Pharmacol.",
    "2024;90(12):3176-3189. doi:10.1111/bcp.16194",
    sep = " "
  )
  vignette    <- "Zhang_2024_nedosiran"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Structure read from Zhang 2024 Figure 3 (page 3183),
  # which is the ONLY place the transit-chain lengths appear -- neither the
  # text nor Table 1 states how many transit compartments each pathway has.
  # Figure 3: the slow pathway is input -> Transit 1 -> central (2 states in
  # series, both leaving at ka1); the fast pathway is input -> Transit 1 ->
  # Transit 2 -> Transit 3 -> central (4 states in series, all leaving at ka2).
  compartmentData <- list(
    depot1         = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    transit1       = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    depot2         = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    transit2       = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    transit3       = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    transit4       = list(analyte = "nedosiran", units = "mg",        specimen = "administration site", verified = TRUE),
    central        = list(analyte = "nedosiran", units = "mg",        specimen = "plasma",              verified = TRUE),
    peripheral1    = list(analyte = "nedosiran", units = "mg",        specimen = "plasma",              verified = TRUE),
    effect         = list(analyte = "nedosiran", units = "ng/mL",     specimen = "not applicable",      verified = TRUE),
    oxalate_urine  = list(analyte = "oxalate",   units = "umol/24 h", specimen = "urine",               verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with a 70 kg reference. Exponents fixed to 0.75 on CL/F and Q/F and to 1 on Vc/F and Vp/F (Zhang 2024 Table 1, both rows marked FIX). A separately estimated exponent (0.83) is shared by ka1 and ka2.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, body-surface-area normalised",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhang 2024 reports this covariate as eGFR in mL/min/1.73 m^2; eGFR is a registered source alias of the canonical CRCL column. Power effects normalised to a reference of 90 mL/min/1.73 m^2 on CL/F (exponent 0.87) and Vc/F (exponent 0.19). The published equations write the covariate as eGFR*RMF, where RMF is a renal maturation function (Anderson & Holford 2011, cited as reference 31 and not on disk) applied ONLY when deriving eGFR for the virtual paediatric cohort of the simulation section; RMF = 1 in adults and no fitted parameter depends on its form, so this column carries the already-matured eGFR.",
      source_name        = "eGFR"
    ),
    DIS_PH1 = list(
      description        = "Primary hyperoxaluria type 1 patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer or PH2 patient -- the reference group pools both)",
      notes              = "Zhang 2024 Table 1 'Covariate effect of PH1 on ka1' = 1.54, a multiplicative factor on the slow-pathway absorption rate constant. The reference group is NOT the healthy volunteers alone: the pooled PK dataset is 85 non-PH volunteers, 46 PH1 and 12 PH2 patients, and only PH1 carries the effect, so DIS_PH1 = 0 covers volunteers AND PH2 patients. Discussion: the PH1 effect was placed on ka1 rather than Vc/F because it described the overall PK better, and it produced no clinically meaningful change in Cmax,ss or AUC0-tau,ss.",
      source_name        = "PH1"
    )
  )

  # Covariates the paper screened but did not retain in the final model. They
  # are documented here for provenance only and are deliberately absent from
  # model(); see checkModelConventions() treatment of covariatesDataExcluded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Zhang 2024 Discussion: 'After accounting for weight, other size metrics, such as age, BSA and BMI, were found not to be significant covariates to the PKs of nedosiran.' Also not significant on the PD model. Full screening list is Supporting Information Table S2 (not on disk)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened as an alternative size metric and not retained after body weight was included (Zhang 2024 Discussion)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as an alternative size metric and not retained after body weight was included (Zhang 2024 Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 143L,
    n_studies      = 5L,
    age_range      = "adults and children aged >= 6 years (14 participants aged 9 to <18 years contributed PK, 10 contributed PD); exact range is in Supporting Information Table S1, which is not on disk",
    weight_range   = "not reported in the main text (Supporting Information Table S1); the simulation reference subject is 70 kg and the weight-banded dosing cut-point is 50 kg",
    sex_female_pct = NA_real_,
    disease_state  = "85 healthy adult volunteers without a primary hyperoxaluria diagnosis, 46 patients with primary hyperoxaluria type 1 (PH1) and 12 with PH2; PHYOX5 additionally enrolled non-PH adults with severe renal impairment or kidney failure",
    dose_range     = "1.5, 3.0 or 6.0 mg/kg single s.c. dose (PHYOX1, PHYOX6); 170 mg s.c. single dose (PHYOX5); repeat flat s.c. doses of 170 mg (nedosiran sodium salt, 160 mg free-acid equivalent) once monthly, or 136 mg once monthly for body weight <50 kg (PHYOX2, PHYOX3)",
    regions        = "not reported in the main text; PHYOX6 was an ethno-bridging study in Japanese and Caucasian healthy adults",
    renal_function = "eGFR spanned normal through end-stage renal disease; all participants with eGFR <30 mL/min/1.73 m^2 were non-PH volunteers from PHYOX5",
    notes          = "PK: 1978 plasma nedosiran concentrations from 143 participants across PHYOX1, PHYOX2, PHYOX3, PHYOX5 and PHYOX6. PD: 588 24-h urinary oxalate observations from the 46 PH1 patients in PHYOX1, PHYOX2 and PHYOX3 (adjusted per 1.73 m^2 body surface area in participants aged <18 years). Placebo-arm PD data (n = 11) were excluded after showing no apparent change in 24-h Uox for up to 28 weeks. 13.1% of PK observations were below the 1.0 ng/mL LLOQ and were excluded. Baseline demographics are tabulated in Supporting Information Table S1, which is not on disk. Estimation: NONMEM VII level 4.0 with PsN 5.0.0; standard errors by importance sampling."
  )

  ini({
    # ==========================================================================
    # Disposition -- Zhang 2024 Table 1 (POP-PK model parameter estimates).
    #
    # SOURCE-CONFLICT NOTE: the in-text equations in Results section 3.2 print a
    # DIFFERENT parameter set from Table 1 (CL/F 5.73 vs 5.69; eGFR-on-CL
    # exponent 0.864 vs 0.87; Vc/F 145 vs 147; ka1 0.214 vs 0.22; WT-on-ka
    # 0.816 vs 0.83; PH1-on-ka1 1.56 vs 1.54; ka2 16.6 vs 16.7). These are not
    # roundings of one another, so the two sets come from different runs. This
    # file uses TABLE 1 throughout: it is the formal parameter-estimates table
    # (it carries %RSE and eta shrinkage), the text points to it for "the
    # parameter estimates for the POP-PK model", and it is the only source that
    # covers every parameter the model needs. See the vignette Errata.
    # ==========================================================================
    lcl   <- log(5.69)    ; label("Apparent clearance CL/F at 70 kg and eGFR 90 (L/h)")                      # Table 1: CL/F = 5.69 L/h (%RSE 11.0)
    lvc   <- log(147)     ; label("Apparent central volume Vc/F at 70 kg and eGFR 90 (L)")                   # Table 1: Vc/F = 147 L (%RSE 3.7)
    lq    <- log(2.45)    ; label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)")                # Table 1: Q/F = 2.45 L/h (%RSE 13.1)
    lvp   <- log(5300)    ; label("Apparent peripheral volume Vp/F at 70 kg (L)")                            # Table 1: Vp/F = 5300 L (%RSE 26.8)
    lvmax <- log(3.94)    ; label("Michaelis-Menten maximum metabolic rate Vmax (mg/h)")                     # Table 1: Vmax = 3.94 mg/h (%RSE 15.6)
    lkm   <- log(293)     ; label("Michaelis-Menten constant KM (ng/mL)")                                    # Table 1: KM = 293 ng/mL (%RSE 15.0)

    # ==========================================================================
    # Dual parallel absorption -- Zhang 2024 Table 1 and Figure 3.
    #
    # FR1 is the fraction of the dose absorbed via the SLOW pathway. It is
    # carried on the logit scale (the library's logitffo idiom) so individual
    # values stay inside (0, 1); see the etalogitffo comment for the evidence
    # that the paper estimated it on a bounded scale too.
    # ==========================================================================
    logitffo <- qlogis(0.698) ; label("Logit of FR1, the fraction of dose absorbed via the slow pathway (logit units)")  # Table 1: FR1 = 69.80% (%RSE 1.3)
    lka1     <- log(0.22)     ; label("Slow-pathway transit rate constant ka1 at 70 kg, non-PH1 (1/h)")                  # Table 1: ka1 = 0.22 1/h (%RSE 4.9)
    lka2     <- log(16.7)     ; label("Fast-pathway transit rate constant ka2 at 70 kg (1/h)")                           # Table 1: ka2 = 16.7 1/h (%RSE 2.5)

    # ==========================================================================
    # Covariate effects -- Zhang 2024 Table 1.
    # ==========================================================================
    e_wt_cl       <- fixed(0.75) ; label("Allometric exponent for WT on CL/F and Q/F (unitless)")            # Table 1: "Exponent for WT on CL/F and Q/F" = 0.75 FIX
    e_wt_vc       <- fixed(1)    ; label("Allometric exponent for WT on Vc/F and Vp/F (unitless)")           # Table 1: "Exponent for WT on Vc/F and Vp/F" = 1 FIX
    e_crcl_cl     <- 0.87        ; label("Exponent for eGFR on CL/F, normalised to 90 mL/min/1.73 m^2 (unitless)")  # Table 1: 0.87 (%RSE 17.5)
    e_crcl_vc     <- 0.19        ; label("Exponent for eGFR on Vc/F, normalised to 90 mL/min/1.73 m^2 (unitless)")  # Table 1: 0.19 (%RSE 18.0)
    e_wt_ka       <- 0.83        ; label("Allometric exponent for WT on ka1 and ka2 (unitless)")             # Table 1: "Exponent for WT on ka1 and ka2" = 0.83 (%RSE 11.5)
    e_dis_ph1_ka1 <- 1.54        ; label("Multiplicative effect of PH1 disease status on ka1 (unitless)")    # Table 1: "Covariate effect of PH1 on ka1" = 1.54 (%RSE 7.8)

    # ==========================================================================
    # Pharmacodynamics -- Zhang 2024 Table 2 (POP-PKPD model parameter
    # estimates) and Figure 3 (equations for Eff and ke0).
    #
    # Table 2 reports Kout in 1/week and the effect-compartment equilibration
    # half-life lambda in weeks; both are converted to the model's hour time
    # base by the /168 and *168 factors written out below.
    # ==========================================================================
    lrbase <- log(1420)                       ; label("Baseline 24-h urinary oxalate (umol/24 h)")                      # Table 2: baseline 24-h Uox = 1420 umol/24 h (%RSE 4.7)
    lkout  <- log(0.338 / 168)                ; label("First-order elimination rate of 24-h Uox, Kout (1/h)")           # Table 2: Kout = 0.338 1/week (%RSE 12.5); 1 week = 168 h
    lic50  <- log(1.04)                       ; label("Half-maximal inhibitory concentration in the effect compartment IC50 (ng/mL)")  # Table 2: IC50 = 1.04 ng/mL (%RSE 24.4)
    lhill  <- log(2.56)                       ; label("Hill coefficient gamma of the inhibitory concentration-effect relationship (unitless)")  # Table 2: gamma = 2.56 (%RSE 14.2)
    lke0   <- log(log(2) / (21.9 * 168))      ; label("Effect-compartment equilibration rate constant ke0 (1/h)")       # Table 2: lambda = 21.9 weeks (%RSE 17.3); Figure 3: ke0 = ln(2)/lambda

    # Imax is a fraction bounded in (0, 1) and is carried on the logit scale.
    # The paper prints it as a percentage (64.5%) with a "%CV" between-subject
    # variability; see the etalogitimax comment for why that variability cannot
    # be log-normal on Imax.
    logitimax <- qlogis(0.645)                ; label("Logit of the maximum inhibitory effect Imax (logit units)")      # Table 2: Imax = 64.5% (%RSE 3.3)

    # ==========================================================================
    # Between-subject variability.
    #
    # SCALE CONVENTION. Zhang 2024 Tables 1 and 2 report every variability term
    # on a STANDARD-DEVIATION scale, never as a variance:
    #   * the log-normal BSVs are labelled "%CV" and the table footnote defines
    #     %CV as the "approximate coefficient of variation", i.e. 100 * omega;
    #     omega^2 = CV^2 is used below;
    #   * the additive PD residual is labelled explicitly "(SD, umol/24 h) 205"
    #     -- read as a variance it would imply an SD of 14 umol/24 h against a
    #     1420 umol/24 h baseline, which is not credible;
    #   * therefore the lone unlabelled entry, "Between-subject variability for
    #     FR1 (additive) 0.35", is also an SD (on FR1's estimation scale).
    #
    # This reading is confirmed quantitatively against the paper's own
    # simulation, Figure 5C -- see the etalogitimax note below and the
    # vignette's "Recovering the variability scale" section.
    # ==========================================================================
    etalvc      ~ 0.271^2   # Table 1: BSV Vc/F  = 27.1 %CV (%RSE 11.1, shrinkage 20.6%)
    etalka1     ~ 0.421^2   # Table 1: BSV ka1   = 42.1 %CV (%RSE  9.5, shrinkage 16.3%)
    etalka2     ~ 0.276^2   # Table 1: BSV ka2   = 27.6 %CV (%RSE  7.5, shrinkage 13.3%)
    etalvmax    ~ 0.353^2   # Table 1: BSV Vmax  = 35.3 %CV (%RSE 14.8, shrinkage  9.2%)

    # Table 1: BSV FR1 = 0.35 "(additive)" (%RSE 12.6, shrinkage 30.2%). Read as
    # an SD on FR1's estimation scale, which must be bounded: an SD of 0.35
    # added directly to FR1 = 0.698 would put ~22% of subjects outside [0, 1]
    # and send the fast-pathway fraction (1 - FR1) negative. Carried here on the
    # logit scale, giving FR1 5th-95th percentiles of 0.565 to 0.804.
    etalogitffo ~ 0.35^2

    etalrbase     ~ 0.317^2  # Table 2: BSV baseline 24-h Uox = 31.7 %CV (%RSE  9.3, shrinkage  4.9%)
    etalic50      ~ 0.86^2   # Table 2: BSV IC50              = 86   %CV (%RSE 15.5, shrinkage 28.3%)

    # Table 2: BSV Imax = 29.4 %CV (%RSE 61.5, shrinkage 36.6%). Log-normal on
    # Imax is FALSIFIED by the paper's own Figure 5C: it would put Imax > 1 for
    # 6.8% of subjects, driving kin*(1 - Eff) negative and the simulated 2.5th
    # percentile of 24-h Uox below zero, whereas Figure 5C shows that percentile
    # plateauing at roughly +245 umol/24 h. Taking omega = 0.294 on the LOGIT
    # scale instead reproduces the published plateau interval almost exactly
    # (predicted median/2.5th/97.5th = 485/237/991 vs read 485/245/990). See the
    # vignette's "Recovering the variability scale" section.
    etalogitimax  ~ 0.294^2

    # ==========================================================================
    # Residual unexplained variability.
    # ==========================================================================
    propSd    <- 0.277 ; label("Proportional residual error on plasma nedosiran (fraction)")   # Table 1: proportional RUV = 27.70% (%RSE 5.2); Results 3.2: "an additive error model on the log-scale ... equates to a proportional RUV on the normal scale"
    addSd_uox <- 205   ; label("Additive residual error on 24-h urinary oxalate (umol/24 h)")  # Table 2: additive RUV = 205 umol/24 h SD (%RSE 8.7)
  })

  model({
    # ------------------------------------------------------------------------
    # 1. Individual disposition parameters.
    #    Zhang 2024 Results 3.2 equations, with Table 1 point estimates:
    #      CL/F = 5.69 * (WT/70)^0.75 * (eGFR/90)^0.87
    #      Vc/F = 147  * (WT/70)^1    * (eGFR/90)^0.19
    #    eGFR acts on CL/F and Vc/F only; Q/F and Vp/F carry weight alone.
    #    No BSV was estimated on CL/F, Q/F, Vp/F or KM (Table 1).
    # ------------------------------------------------------------------------
    cl   <- exp(lcl)             * (WT / 70)^e_wt_cl * (CRCL / 90)^e_crcl_cl
    vc   <- exp(lvc + etalvc)    * (WT / 70)^e_wt_vc * (CRCL / 90)^e_crcl_vc
    q    <- exp(lq)              * (WT / 70)^e_wt_cl
    vp   <- exp(lvp)             * (WT / 70)^e_wt_vc
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)

    # ------------------------------------------------------------------------
    # 2. Individual absorption parameters.
    #      ka1 = 0.22 * (WT/70)^0.83 * 1.54^PH1
    #      ka2 = 16.7 * (WT/70)^0.83
    #    ffo is FR1, the slow-pathway dose fraction; 1 - ffo goes to the fast
    #    pathway (Figure 3).
    # ------------------------------------------------------------------------
    ka1 <- exp(lka1 + etalka1) * (WT / 70)^e_wt_ka * e_dis_ph1_ka1^DIS_PH1
    ka2 <- exp(lka2 + etalka2) * (WT / 70)^e_wt_ka
    ffo <- expit(logitffo + etalogitffo)

    # ------------------------------------------------------------------------
    # 3. PD parameters. kin is not estimated directly: Figure 3 defines the
    #    24-h Uox pool by its baseline, so kin = baseline * kout holds the
    #    untreated system at rbase.
    # ------------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout)
    kin   <- rbase * kout
    imax  <- expit(logitimax + etalogitimax)
    ic50  <- exp(lic50 + etalic50)
    hill  <- exp(lhill)
    ke0   <- exp(lke0)

    # ------------------------------------------------------------------------
    # 4. Micro-constants.
    # ------------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Plasma concentration in ng/mL: amounts are mg and volumes are L, so
    # central/vc is mg/L and 1 mg/L = 1000 ng/mL. KM and IC50 are both reported
    # in ng/mL, so the Michaelis-Menten and Imax terms use this same scale.
    Cc <- 1000 * central / vc

    # Inhibitory effect driven by the effect-compartment concentration
    # (Figure 3): Eff = Imax * Ceff^gamma / (Ceff^gamma + IC50^gamma).
    eff <- imax * effect^hill / (effect^hill + ic50^hill)

    # ------------------------------------------------------------------------
    # 5. ODE system (Zhang 2024 Figure 3).
    #    Slow chain: depot1 -ka1-> transit1 -ka1-> central.
    #    Fast chain: depot2 -ka2-> transit2 -ka2-> transit3 -ka2-> transit4
    #                -ka2-> central.
    #    Central eliminates in parallel by linear CL/F and by a
    #    Michaelis-Menten pathway (liver ASGPR-mediated uptake, Discussion).
    # ------------------------------------------------------------------------
    d/dt(depot1)      <- -ka1 * depot1
    d/dt(transit1)    <-  ka1 * depot1 - ka1 * transit1

    d/dt(depot2)      <- -ka2 * depot2
    d/dt(transit2)    <-  ka2 * depot2   - ka2 * transit2
    d/dt(transit3)    <-  ka2 * transit2 - ka2 * transit3
    d/dt(transit4)    <-  ka2 * transit3 - ka2 * transit4

    d/dt(central)     <-  ka1 * transit1 + ka2 * transit4 -
                          kel * central - k12 * central + k21 * peripheral1 -
                          vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Sheiner effect compartment; the state holds a concentration (ng/mL).
    d/dt(effect)      <-  ke0 * (Cc - effect)

    # Indirect response: nedosiran inhibits production of 24-h urinary oxalate.
    d/dt(oxalate_urine) <- kin * (1 - eff) - kout * oxalate_urine
    oxalate_urine(0)    <- rbase

    # ------------------------------------------------------------------------
    # 6. Dose split between the two parallel absorption pathways. The dose is
    #    supplied as two records at the same time with the same amt, one to
    #    depot1 and one to depot2 (the Baverel 2015 tralokinumab idiom).
    # ------------------------------------------------------------------------
    f(depot1) <- ffo
    f(depot2) <- 1 - ffo

    # ------------------------------------------------------------------------
    # 7. Observations.
    # ------------------------------------------------------------------------
    uox <- oxalate_urine
    Cc  ~ prop(propSd)
    uox ~ add(addSd_uox)
  })
}
