Ashraf_2018_s_ketamine_ticlopidine <- function() {
  description <- paste(
    "Semi-PBPK. Joint semi-mechanistic population PK model of oral and",
    "intravenous S-ketamine, its primary metabolite norketamine, and the",
    "co-administered CYP2B6 mechanism-based inhibitor ticlopidine, with a",
    "mechanistic-static drug-drug-interaction layer, in healthy adult",
    "volunteers (Ashraf 2018). S-ketamine has three-compartment mammillary",
    "disposition with first-order absorption; norketamine has",
    "two-compartment disposition; ticlopidine has two-compartment",
    "disposition behind a four-transit absorption chain. All three",
    "compounds are eliminated by well-stirred gut-wall and hepatic",
    "clearance driven by weight-scaled physiological blood flows",
    "(QH = 3.75 * WT^0.75 L/h), with the S-ketamine gut wall additionally",
    "using the QGUT permeability model. For S-ketamine and norketamine the",
    "gut wall, portal vein and liver are solved by the simplified",
    "quasi-steady-state approximation (algebraic, not ODE states); for",
    "ticlopidine they are explicit ODE states because the portal-vein",
    "inhibitor concentration drives the interaction. All S-ketamine",
    "extracted at the gut wall and liver becomes norketamine (Fmet fixed",
    "to 1), so the norketamine disposition parameters are apparent.",
    "Ticlopidine inhibits hepatic CYP2B6 through the mechanistic static",
    "model: a reversible component AH = 1 / (1 + IPV / ki) and a",
    "time-dependent component BH = kdeg / (kdeg + kinact * IPV /",
    "(KI + IPV)) combine as fCLint = AH * BH * fm2B6 + (1 - fm2B6), which",
    "multiplies the S-ketamine intrinsic hepatic clearance. With no",
    "ticlopidine present the interaction term is exactly 1. Enzyme",
    "inactivation constants are in-vitro values fixed from Obach 2007; the",
    "CYP2B6-metabolised fraction 0.63 was estimated by log-likelihood",
    "profiling. No inhibition model is applied to norketamine or to",
    "ticlopidine's own (auto-inhibited) metabolism because neither could be",
    "supported by the data.",
    sep = " "
  )
  reference <- paste(
    "Ashraf MW, Peltoniemi MA, Olkkola KT, Neuvonen PJ, Saari TI.",
    "Semimechanistic Population Pharmacokinetic Model to Predict the",
    "Drug-Drug Interaction Between S-ketamine and Ticlopidine in Healthy",
    "Human Volunteers. CPT Pharmacometrics Syst Pharmacol.",
    "2018;7(10):687-697. doi:10.1002/psp4.12346",
    sep = " "
  )
  vignette <- "Ashraf_2018_s_ketamine_ticlopidine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Total body weight (kg).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives every physiological blood flow in the model through the",
        "hepatic blood flow QH = 3.75 * WT^0.75 L/h (Ashraf 2018 Table 1,",
        "from Brown 1997), from which the portal-vein, hepatic-artery,",
        "intestinal, mucosal and villous flows are derived as fixed",
        "fractions. It is the only subject-level covariate in the model.",
        "Observed weights across the five pooled studies span 50-88 kg",
        "(study means 59-70 kg; Supplementary Information S1 Table S2).",
        sep = " "
      ),
      source_name = "WTKG"
    ),
    DOSE_TICLOPIDINE_MG = list(
      description = "Ticlopidine dose amount currently being administered (mg); 0 when ticlopidine is not being given.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the mechanistic-static-model portal-vein inhibitor",
        "concentration as the Fahmi 2008 gut-input term",
        "fu_b,TIC * F * ka,TIC * Dose / QH. In the deposited NONMEM control",
        "stream (Supplementary Information S3, $DES) this is the reserved",
        "data item AMT: 'TCE = TFUB*(A(9)+((F7*TKA*AMT)/QH))'. rxode2 model",
        "code cannot read the amt of a dose record, so the administered",
        "ticlopidine dose is carried as an explicit column. Set it to the",
        "ticlopidine dose (250 mg in Ashraf 2018 Study III) on records",
        "inside the ticlopidine dosing period and to 0 during washout and",
        "in the placebo phase -- the washout zero is what lets CYP2B6",
        "activity recover over 4-5 days as in Ashraf 2018 Figure 3b.",
        sep = " "
      ),
      source_name = "AMT"
    ),
    CONMED_TICLOPIDINE = list(
      description = "Ticlopidine study-phase indicator (1 = ticlopidine phase, 0 = placebo phase).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo phase)",
      notes = paste(
        "Selects between the two proportional residual-error magnitudes",
        "that Ashraf 2018 estimated per study phase for S-ketamine and for",
        "norketamine (Table 2 rows RV SK,PLAC / RV SK,TICLO and",
        "RV NK,PLAC / RV NK,TICLO). In the deposited control stream this is",
        "the TICLO data item, which selects EPS(1)/EPS(2) in the placebo",
        "phase and EPS(3)/EPS(4)/EPS(5) in the ticlopidine phase. It has no",
        "structural effect: the control stream also gates the interaction",
        "term on TICLO == 1, but that gate is redundant because with no",
        "ticlopidine in the system the interaction term evaluates to",
        "exactly 1 (see the model block).",
        sep = " "
      ),
      source_name = "TICLO"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "S-ketamine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    central_snk = list(analyte = "norketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_snk = list(analyte = "norketamine", units = "mg", specimen = "plasma", verified = TRUE),
    depot_tic = list(analyte = "ticlopidine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1_tic = list(analyte = "ticlopidine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2_tic = list(analyte = "ticlopidine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3_tic = list(analyte = "ticlopidine", units = "mg", specimen = "administration site", verified = TRUE),
    transit4_tic = list(analyte = "ticlopidine", units = "mg", specimen = "administration site", verified = TRUE),
    gut_tic = list(analyte = "ticlopidine", units = "mg", specimen = "tissue", verified = TRUE),
    portal_tic = list(analyte = "ticlopidine", units = "mg", specimen = "whole blood", verified = TRUE),
    liver_tic = list(analyte = "ticlopidine", units = "mg", specimen = "tissue", verified = TRUE),
    central_tic = list(analyte = "ticlopidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_tic = list(analyte = "ticlopidine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 41L,
    n_studies = 5L,
    age_range = "19-35 years (per-study means 22.0-27.5 years)",
    weight_range = "50-88 kg (per-study means 59-70 kg)",
    sex_female_pct = 48.8,
    disease_state = "Healthy, non-smoking adult volunteers; no concomitant drug therapy for at least 14 days before study entry and no CYP-inhibiting or CYP-inducing drugs for 4 weeks before study entry.",
    dose_range = paste(
      "S-ketamine 0.2-0.3 mg/kg as an oral syrup (Studies I, II, III, IV",
      "oral part, V) and 0.1 mg/kg intravenously over 2 minutes (Study IV",
      "intravenous part). Ticlopidine 250 mg twice daily for 6 days before",
      "the S-ketamine dose (Study III phase 3).",
      sep = " "
    ),
    regions = "Single centre, Turku, Finland (clinical phases run June 2008 - March 2010).",
    co_medication = "None permitted. The ticlopidine phase of Study III is the only interaction arm carried into this model; the clarithromycin, St John's wort, itraconazole, rifampicin and grapefruit-juice arms of the same five studies were excluded and only their placebo phases pooled.",
    n_observations = paste(
      "13 timed plasma samples per individual (pre-dose and 20 min, 40 min,",
      "1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12 and 24 h after the S-ketamine dose;",
      "ticlopidine sampled 1.33-12 h). 67 placebo-phase and 11",
      "ticlopidine-phase S-ketamine / norketamine profiles (Ashraf 2018",
      "Figure S2 caption). LLOQ 0.025 ng/mL for S-ketamine and",
      "norketamine, 10 ng/mL for ticlopidine.",
      sep = " "
    ),
    notes = paste(
      "Data pooled from the placebo phases of five randomised,",
      "placebo-controlled, crossover healthy-volunteer studies (Hagelberg",
      "2010; Peltoniemi 2011, 2012a, 2012b, 2012c) plus the complete",
      "placebo-and-ticlopidine dataset of the ticlopidine interaction study",
      "(Peltoniemi 2011). 13 of the 41 volunteers took part in two studies",
      "and one in three; repeat participations were modelled as separate",
      "individuals to avoid the inter-occasion variability the model could",
      "not support (Supplementary Information S1 footnote). Demographics in",
      "Supplementary Information S1 Table S2; final parameter estimates in",
      "Table 2.",
      sep = " "
    )
  )

  ini({
    # ===================================================================
    # All estimated values are the final NONMEM mean theta estimates in
    # Ashraf 2018 Table 2. All fixed values are the FIX entries of the
    # deposited final control stream (Supplementary Information S3,
    # $THETA) and/or Table 1. The control stream's non-FIX $THETA entries
    # are INITIAL estimates (e.g. KV1 13.5 vs the final 14.4) and are NOT
    # used here; only its FIX entries are final.
    # ===================================================================

    # ---------------- S-ketamine structural parameters -----------------
    lka <- log(1.76)
    label("S-ketamine first-order absorption rate constant ka (1/h)") # Table 2, K a,SK = 1.76 /h (RSE 16.1%; SIR median 1.78 [1.47, 2.11])

    lvc <- log(14.4)
    label("S-ketamine central volume Vc (L)") # Table 2, V C,SK = 14.4 L (RSE 43.1%; SIR median 14.4 [5.90, 22.5])

    lq <- log(287)
    label("S-ketamine central-to-first-peripheral clearance Q1 (L/h)") # Table 2, Q PER1,SK = 287 L/h (RSE 11.3%)

    lvp <- log(102)
    label("S-ketamine first peripheral volume Vp1 (L)") # Table 2, V PER1,SK = 102 L (RSE 8.0%)

    lq2 <- log(22.4)
    label("S-ketamine central-to-second-peripheral clearance Q2 (L/h)") # Table 2, Q PER2,SK = 22.4 L/h (RSE 6.0%)

    lvp2 <- log(180)
    label("S-ketamine second peripheral volume Vp2 (L)") # Table 2, V PER2,SK = 180 L (RSE 3.8%)

    lclint_liver <- log(301)
    label("S-ketamine intrinsic hepatic clearance CLint,H at no inhibitor (L/h)") # Table 2, CL INT,H,SK = 301 L/h (RSE 14.2%; SIR median 304 [260, 352])

    lclint_gut <- log(1.19)
    label("S-ketamine intrinsic gut-wall clearance CLint,GW (L/h)") # Table 2, CL INT,GW,SK = 1.19 L/h (RSE 39.5%; the 1:253 gut-to-liver ratio quoted in Results)

    clperm <- fixed(4.1)
    label("S-ketamine gut-mucosal permeability clearance CLp for the QGUT model (L/h)") # $THETA 9 = 4.1 FIX; derived in Suppl. S2 from PAMPA Papp 1.49e-6 cm/s -> Peff,man 1.74e-4 cm/s x 0.66 m2 = 1.15 cm3/s = 4.1 L/h

    dur_iv <- fixed(0.03333)
    label("Duration of the intravenous S-ketamine injection (h)") # $PK D2 = 0.03333 h, i.e. the 2-minute i.v. injection of Study IV; only used when a dose record carries rate = -2

    # ---------------- Norketamine structural parameters ----------------
    # Apparent values: Fmet was set to 1 (all extracted S-ketamine becomes
    # norketamine) while the true metabolic fraction is about 0.80, so
    # these are apparent rather than absolute (Table 2 footnote b).
    lvc_snk <- log(88.4)
    label("Norketamine apparent central volume Vc (L)") # Table 2, V C,NK = 88.4 L (RSE 4.8%)

    lq_snk <- log(19.9)
    label("Norketamine apparent central-to-peripheral clearance Q1 (L/h)") # Table 2, Q PER1,NK = 19.9 L/h (RSE 12.3%)

    lvp_snk <- log(88.9)
    label("Norketamine apparent first peripheral volume Vp1 (L)") # Table 2, V PER1,NK = 88.9 L (RSE 5.6%)

    lclint_liver_snk <- log(73.5)
    label("Norketamine apparent intrinsic hepatic clearance CLint,H (L/h)") # Table 2, CL INT,H,NK = 73.5 L/h (RSE 7.8%)

    lclint_gut_snk <- log(44.4)
    label("Norketamine apparent intrinsic gut-wall clearance CLint,GW (L/h)") # Table 2, CL INT,GW,NK = 44.4 L/h (RSE 35.8%)

    # ---------------- Ticlopidine structural parameters ----------------
    lka_tic <- fixed(log(3.3))
    label("Ticlopidine absorption / transit rate constant ka (1/h)") # Table 2, K a,TIC = 3.3 /h (FIXED); $THETA 15 = 3.3 FIX, fixed to the published value of Palacharla 2018 because the sparse early absorption data could not support estimation

    lvc_tic <- log(50.3)
    label("Ticlopidine central volume Vc (L)") # Table 2, V C,TIC = 50.3 L (RSE 13.9%)

    lq_tic <- log(26.3)
    label("Ticlopidine central-to-peripheral clearance Q1 (L/h)") # Table 2, Q PER1,TIC = 26.3 L/h (RSE 24.4%)

    lvp_tic <- log(191)
    label("Ticlopidine first peripheral volume Vp1 (L)") # Table 2, V PER1,TIC = 191 L (RSE 67%; SIR median 208.7 [130, 393])

    lclint_liver_tic <- log(1505)
    label("Ticlopidine intrinsic hepatic clearance CLint,H (L/h)") # Table 2, CL INT,H,TIC = 1505 L/h (RSE 27.9%)

    clint_gut_tic <- fixed(0)
    label("Ticlopidine intrinsic gut-wall clearance CLint,GW (L/h)") # Table 2, CL INT,GW,TIC = 0 (FIXED); $THETA 20 = 0 FIX. Kept on the natural scale because log(0) is undefined. The gut wall therefore acts as a pure physiological transit for ticlopidine and metabolism happens only in the liver

    # ---------------- Drug binding constants (Table 1) ------------------
    fu <- fixed(0.70)
    label("S-ketamine fraction unbound in plasma (unitless)") # Table 1, f u,SK = 0.70 (Peltoniemi 2016)

    fu_snk <- fixed(0.50)
    label("Norketamine fraction unbound in plasma (unitless)") # Table 1, f u,NK = 0.50 (Hijazi 2002)

    fu_tic <- fixed(0.02)
    label("Ticlopidine fraction unbound in plasma (unitless)") # Table 1, f u,TIC = 0.02 (Ito 1992)

    bp <- fixed(0.50)
    label("S-ketamine blood-to-plasma concentration ratio (unitless)") # Table 1, BP RATIO = 0.50 for S-ketamine (Launiainen 2014)

    bp_snk <- fixed(1)
    label("Norketamine blood-to-plasma concentration ratio (unitless)") # Table 1, BP RATIO = 1 for norketamine (Yang 2007)

    bp_tic <- fixed(1)
    label("Ticlopidine blood-to-plasma concentration ratio (unitless)") # Table 1, BP RATIO = 1 for ticlopidine (Yang 2007)

    fu_gut <- fixed(1)
    label("S-ketamine fraction unbound in the gut wall (unitless)") # Table 1, f u,GW,SK = 1 (Yang 2007)

    fu_gut_snk <- fixed(1)
    label("Norketamine fraction unbound in the gut wall (unitless)") # Table 1, f u,GW,NK = 1 (Yang 2007)

    fu_gut_tic <- fixed(1)
    label("Ticlopidine fraction unbound in the gut wall (unitless)") # Table 1, f u,GW,TIC = 1 (Yang 2007)

    # ---------------- Physiological parameters (Table 1) ----------------
    qh_coef <- fixed(3.75)
    label("Hepatic blood flow coefficient in QH = qh_coef * WT^e_wt_qh (L/h/kg^0.75)") # Table 1, Q H = 3.75 * WTKG^0.75 (Brown 1997)

    e_wt_qh <- fixed(0.75)
    label("Allometric exponent of body weight on hepatic blood flow (unitless)") # Table 1, Q H = 3.75 * WTKG^0.75 (Brown 1997)

    fq_portal <- fixed(0.75)
    label("Portal-vein blood flow as a fraction of hepatic blood flow (unitless)") # Table 1, Q PV = 0.75 * Q H (Williams 1989)

    fq_hepatic_artery <- fixed(0.25)
    label("Hepatic-artery blood flow as a fraction of hepatic blood flow (unitless)") # Table 1, Q HA = 0.25 * Q H (Williams 1989)

    fq_intestinal <- fixed(0.40)
    label("Intestinal blood flow as a fraction of hepatic blood flow (unitless)") # Table 1, Q INT = 0.40 * Q H (Williams 1989)

    fq_mucosal <- fixed(0.80)
    label("Mucosal blood flow as a fraction of intestinal blood flow (unitless)") # Table 1, Q MU = 0.80 * Q INT (Yang 2007)

    fq_villous <- fixed(0.60)
    label("Villous blood flow as a fraction of mucosal blood flow (unitless)") # Table 1, Q VI = 0.60 * Q MU (Yang 2007)

    vgut <- fixed(1)
    label("Gut-wall volume (L)") # Table 1, V GW = 1 L

    vportal <- fixed(1)
    label("Portal-vein volume (L)") # Table 1, V PV = 1 L

    vliver <- fixed(1)
    label("Liver volume (L)") # Table 1, V H = 1 L

    # ---------------- CYP2B6 interaction constants ----------------------
    kirev_2b6 <- fixed(0.031)
    label("Ticlopidine equilibrium dissociation constant ki for reversible CYP2B6 inhibition (umol/L)") # Table 1, k i = 0.031 uM (Obach 2007); $THETA 25 = 0.031 FIX, used in the control stream as TKIN in TAH = 1/(1 + TCE/TKIN)

    ki_2b6 <- fixed(0.57)
    label("Ticlopidine concentration at half-maximal CYP2B6 inactivation rate KI (umol/L)") # Table 1, K I = 0.57 uM (Obach 2007); $THETA 21 = 0.57 FIX

    kinact_2b6 <- fixed(18)
    label("Maximum rate constant of CYP2B6 inactivation by ticlopidine kinact (1/h)") # Table 1, k inact = 0.30 /min (Obach 2007) = 18 /h; $THETA 22 = 18 FIX

    kdeg_2b6 <- fixed(0.017)
    label("Physiological degradation rate constant of hepatic CYP2B6 at zero inhibitor kdeg (1/h)") # $THETA 24 = 0.017 FIX (control-stream LAMBDA). Table 1 prints k deg = 0.00026 /min = 0.0156 /h; the deposited final model uses 0.017 /h and that value is used here (see vignette Errata)

    fm_cyp2b6 <- fixed(0.63)
    label("Fraction of S-ketamine hepatic metabolism mediated by CYP2B6 (unitless)") # Table 1 footnote and Results: estimated in the model by log-likelihood profiling, final value 0.63 (the in vitro value of Palacharla 2018 is 0.60); $THETA 23 = 0.63 FIX

    # ---------------- Interindividual variability -----------------------
    # NONMEM $OMEGA diagonal variances on the exponential IIV model
    # Phi_i = theta * exp(eta_i) (Methods, Stochastic model).
    etalclint_liver ~ 0.25 # Table 2, 'IIV on CL INT,H,SK' = 0.25 (RSE 11.9%; SIR median 0.26 [0.17, 0.37])
    etalclint_gut ~ 2.1 # Table 2, 'IIV on CL INT,GW,SK' = 2.1 (RSE 22.4%; SIR median 2.0 [1.11, 3.67])
    etalka ~ 0.42 # Table 2, 'IIV on K a,SK' = 0.42 (RSE 17.9%; SIR median 0.44 [0.29, 0.66])
    etalvp ~ 0.045 # Table 2, 'IIV on V PER1,SK' = 0.045 (RSE 22.3%; SIR median 0.046 [0.022, 0.076])
    etalclint_liver_snk ~ 0.10 # Table 2, 'IIV on CL INT,H,NK' = 0.10 (RSE 15.1%)
    etalclint_liver_tic ~ 0.12 # Table 2, 'IIV on CL INT,H,TIC' = 0.12 (RSE 76.2%; SIR median 0.13 [0.06, 0.28])

    # ---------------- Residual variability ------------------------------
    # Proportional error model C_obs = C_pred * (1 + eps). Table 2 reports
    # the NONMEM $SIGMA VARIANCES, so each nlmixr2 SD below is the square
    # root of the tabulated value. Ashraf 2018 estimated separate residual
    # magnitudes per study phase for S-ketamine and norketamine ($ERROR
    # EPS(1)-EPS(5), keyed on TICLO and CMT).
    propSd <- sqrt(0.086)
    label("S-ketamine proportional residual SD, placebo phase (fraction)") # Table 2, RV SK,PLAC = 0.086 variance (RSE 2.1%) -> SD 0.293

    propSd_ticlo <- sqrt(0.065)
    label("S-ketamine proportional residual SD, ticlopidine phase (fraction)") # Table 2, RV SK,TICLO = 0.065 variance (RSE 33%) -> SD 0.255

    propSd_snk <- sqrt(0.062)
    label("Norketamine proportional residual SD, placebo phase (fraction)") # Table 2, RV NK,PLAC = 0.062 variance (RSE 2.4%) -> SD 0.249

    propSd_snk_ticlo <- sqrt(0.064)
    label("Norketamine proportional residual SD, ticlopidine phase (fraction)") # Table 2, RV NK,TICLO = 0.064 variance (RSE 31%) -> SD 0.253

    propSd_tic <- sqrt(0.066)
    label("Ticlopidine proportional residual SD (fraction)") # Table 2, RV TICLO = 0.066 variance (RSE 9.5%) -> SD 0.257; ticlopidine is only ever observed in the ticlopidine phase so there is no phase split
  })

  model({
    # =================================================================
    # 1. Physiological blood flows and volumes (Ashraf 2018 Table 1;
    #    control stream $PK "Physiological parameters")
    # =================================================================
    qh <- qh_coef * WT^e_wt_qh
    qpv <- fq_portal * qh
    qha <- fq_hepatic_artery * qh
    qint <- fq_intestinal * qh
    qmu <- fq_mucosal * qint
    qvi <- fq_villous * qmu

    # Unbound fractions in BLOOD, fu,b = fu,plasma / (blood:plasma ratio)
    # (control stream KFUB / NKFUB / TFUB).
    fub <- fu / bp
    fub_snk <- fu_snk / bp_snk
    fub_tic <- fu_tic / bp_tic

    # =================================================================
    # 2. Individual parameters
    # =================================================================
    ka <- exp(lka + etalka)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    clint_liver <- exp(lclint_liver + etalclint_liver)
    clint_gut <- exp(lclint_gut + etalclint_gut)

    vc_snk <- exp(lvc_snk)
    q_snk <- exp(lq_snk)
    vp_snk <- exp(lvp_snk)
    clint_liver_snk <- exp(lclint_liver_snk + etalclint_liver_snk)
    clint_gut_snk <- exp(lclint_gut_snk)

    ka_tic <- exp(lka_tic)
    ktr_tic <- ka_tic # control stream TKTR = TKA: one rate constant drives both the depot and the four transits
    vc_tic <- exp(lvc_tic)
    q_tic <- exp(lq_tic)
    vp_tic <- exp(lvp_tic)
    clint_liver_tic <- exp(lclint_liver_tic + etalclint_liver_tic)

    # =================================================================
    # 3. Mechanistic static drug-drug-interaction model
    #    (Suppl. S2 "Drug-drug interaction model"; control stream $DES)
    #
    #    IPV is the ticlopidine concentration at the enzyme site: the
    #    unbound portal-vein concentration plus the Fahmi 2008 gut-input
    #    term F * ka * Dose / QH. Note that the in vitro constants are in
    #    umol/L while IPV is computed in mg/L; the model is transcribed
    #    exactly as deposited (see vignette Errata).
    # =================================================================
    ipv <- fub_tic * (portal_tic / vportal + ka_tic * DOSE_TICLOPIDINE_MG / qh)
    ah <- 1 / (1 + ipv / kirev_2b6) # reversible (competitive) component
    bh <- kdeg_2b6 / (kdeg_2b6 + kinact_2b6 * ipv / (ki_2b6 + ipv)) # time-dependent component
    zeta <- ah * bh * fm_cyp2b6 + (1 - fm_cyp2b6) # fCL'int,H
    ra_cyp2b6 <- 100 * ah * bh # per cent CYP2B6 activity remaining (Results Eq. for %RA with Activity[I]=0 = 100)
    clint_liver_inh <- clint_liver * zeta # CLint,H,[I] = CLint,H,[I]=0 * fCL'int,H

    # =================================================================
    # 4. Well-stirred extraction at the gut wall and the liver
    # =================================================================
    # S-ketamine gut wall, with the QGUT permeability model
    qgut <- qvi * clperm / (qvi + clperm)
    fgut <- qgut / (qgut + clint_gut * fu_gut)

    # S-ketamine liver, using the (possibly inhibited) intrinsic clearance
    eliver <- clint_liver_inh * fub / (qh + clint_liver_inh * fub)
    fliver <- 1 - eliver

    # Norketamine gut wall (villous flow, no QGUT model) and liver
    fgut_snk <- qvi / (qvi + clint_gut_snk * fu_gut_snk)
    fliver_snk <- qh / (qh + clint_liver_snk * fub_snk)

    # Ticlopidine gut wall and liver
    egut_tic <- clint_gut_tic * fu_gut_tic / (qvi + clint_gut_tic * fu_gut_tic)
    fgut_tic <- 1 - egut_tic
    fliver_tic <- qh / (qh + clint_liver_tic * fub_tic)
    eliver_tic <- 1 - fliver_tic

    # =================================================================
    # 5. Distribution micro-constants
    # =================================================================
    k23 <- q / vc
    k32 <- q / vp
    k24 <- q2 / vc
    k42 <- q2 / vp2
    k56 <- q_snk / vc_snk
    k65 <- q_snk / vp_snk
    k1112 <- q_tic / vc_tic
    k1211 <- q_tic / vp_tic

    # =================================================================
    # 6. S-ketamine: quasi-steady-state gut wall / portal vein / liver
    #    (Suppl. S2 S-ketamine equations; control stream $DES)
    # =================================================================
    a_gut <- ka * depot / (qvi / vgut)
    a_portal <- ((qvi / vgut) * fgut * a_gut + (qpv / vc) * central) / (qpv / vportal)
    a_liver <- ((qha / vc) * central + (qpv / vportal) * a_portal) / (qh / vliver)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- fliver * (qh / vliver) * a_liver -
      (qha / vc) * central - (qpv / vc) * central -
      k23 * central + k32 * peripheral1 -
      k24 * central + k42 * peripheral2
    d/dt(peripheral1) <- k23 * central - k32 * peripheral1
    d/dt(peripheral2) <- k24 * central - k42 * peripheral2
    dur(central) <- dur_iv

    # =================================================================
    # 7. Norketamine: all S-ketamine extracted at the gut wall and the
    #    liver enters the metabolite model directly (Fmet = 1)
    # =================================================================
    a_gut_snk <- (1 - fgut) * a_gut
    a_portal_snk <- ((qvi / vgut) * a_gut_snk * fgut_snk + (qpv / vc_snk) * central_snk) / (qpv / vportal)
    a_liver_snk <- ((qha / vc_snk) * central_snk + (qpv / vportal) * a_portal_snk +
      (qh / vliver) * eliver * a_liver) / (qh / vliver)

    d/dt(central_snk) <- fliver_snk * (qh / vliver) * a_liver_snk -
      (qpv / vc_snk) * central_snk - (qha / vc_snk) * central_snk -
      k56 * central_snk + k65 * peripheral1_snk
    d/dt(peripheral1_snk) <- k56 * central_snk - k65 * peripheral1_snk

    # =================================================================
    # 8. Ticlopidine: four-transit absorption into explicit gut wall,
    #    portal vein and liver compartments (Suppl. S2 ticlopidine
    #    equations; control stream $DES)
    # =================================================================
    d/dt(depot_tic) <- -ka_tic * depot_tic
    d/dt(transit1_tic) <- ka_tic * depot_tic - ktr_tic * transit1_tic
    d/dt(transit2_tic) <- ktr_tic * transit1_tic - ktr_tic * transit2_tic
    d/dt(transit3_tic) <- ktr_tic * transit2_tic - ktr_tic * transit3_tic
    d/dt(transit4_tic) <- ktr_tic * transit3_tic - ktr_tic * transit4_tic
    d/dt(gut_tic) <- ktr_tic * transit4_tic -
      (qvi / vgut) * egut_tic * gut_tic - (qvi / vgut) * fgut_tic * gut_tic
    d/dt(portal_tic) <- (qvi / vgut) * fgut_tic * gut_tic +
      (qpv / vc_tic) * central_tic - (qpv / vportal) * portal_tic
    d/dt(liver_tic) <- (qpv / vportal) * portal_tic + (qha / vc_tic) * central_tic -
      eliver_tic * (qh / vliver) * liver_tic - fliver_tic * (qh / vliver) * liver_tic
    d/dt(central_tic) <- fliver_tic * (qh / vliver) * liver_tic -
      (qpv / vc_tic) * central_tic - (qha / vc_tic) * central_tic -
      k1112 * central_tic + k1211 * peripheral1_tic
    d/dt(peripheral1_tic) <- k1112 * central_tic - k1211 * peripheral1_tic

    # =================================================================
    # 9. Observations. Amounts are mg and volumes L, so A/V is mg/L and
    #    the factor 1000 converts to ng/mL (control stream S2 = KV1/1000,
    #    S5 = NKV5/1000, S11 = TV1/1000).
    # =================================================================
    Cc <- 1000 * central / vc
    Cc_snk <- 1000 * central_snk / vc_snk
    Cc_tic <- 1000 * central_tic / vc_tic

    # Phase-specific proportional residual magnitudes ($ERROR EPS(1)-(4)
    # switched on the TICLO data item).
    sd_cc <- propSd * (1 - CONMED_TICLOPIDINE) + propSd_ticlo * CONMED_TICLOPIDINE
    sd_cc_snk <- propSd_snk * (1 - CONMED_TICLOPIDINE) + propSd_snk_ticlo * CONMED_TICLOPIDINE

    Cc ~ prop(sd_cc)
    Cc_snk ~ prop(sd_cc_snk)
    Cc_tic ~ prop(propSd_tic)
  })
}
