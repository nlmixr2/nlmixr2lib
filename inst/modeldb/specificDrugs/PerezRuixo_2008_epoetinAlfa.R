PerezRuixo_2008_epoetinAlfa <- function() {
  description <- paste(
    "Population PK/PD model for subcutaneous recombinant human",
    "erythropoietin (rHuEPO / epoetin alfa) in healthy adult male",
    "volunteers (Perez-Ruixo 2008). PK is the Olsson-Gisleskog 2007 prior",
    "(two-compartment with linear + Michaelis-Menten elimination and dual",
    "subcutaneous absorption: a fast sequential zero-order infusion into",
    "the depot of duration D1 feeding first-order absorption ka into",
    "central, plus a slower zero-order direct infusion into central of",
    "duration D2 after lag time tlag2; dose-dependent absolute",
    "bioavailability F = F0 + Emax(F)*Dose/(ED50(F)+Dose)). Endogenous EPO",
    "is maintained at the baseline BSL by a constant input rate kEPO",
    "derived from the steady-state balance against linear + MM elimination",
    "(equation 4). The PD layer is the maturation-structured cytokinetic",
    "model D: rHuEPO stimulates the progenitor production rate",
    "kin*C/(SC50+C) into a 10-stage bone-marrow precursor age chain",
    "(transfer rate Np/Tp), which feeds a 10-stage circulating",
    "reticulocyte age chain whose transfer rate (NR/TR)*(S0/SM) is",
    "inhibited by a 5-stage signal transduction (transit time tau)",
    "driven by C/(EC50+C). Output RET = sum of reticulocyte compartments",
    "reproduces the percentage of reticulocytes in % units. No",
    "demographic covariate effects were retained in either layer."
  )
  reference <- paste(
    "Perez-Ruixo JJ, Krzyzanski W, Hing J.",
    "Pharmacodynamic analysis of recombinant human erythropoietin",
    "effect on reticulocyte production rate and age distribution",
    "in healthy subjects.",
    "Clin Pharmacokinet. 2008;47(6):399-415.",
    "doi:10.2165/00003088-200847060-00004.",
    "PK structure and prior typical values + prior CV% adopted from",
    "Olsson-Gisleskog P, Jacqmin P, Perez-Ruixo JJ.",
    "Population pharmacokinetics meta-analysis of recombinant human",
    "erythropoietin in healthy subjects.",
    "Clin Pharmacokinet. 2007;46(2):159-73.",
    "doi:10.2165/00003088-200746020-00004",
    "(MAP-Bayesian POSTHOC PK inputs to the PD analysis; see Table I)."
  )
  vignette <- "PerezRuixo_2008_epoetinAlfa"
  units <- list(time = "h", dosing = "IU", concentration = "IU/L")

  # The reticulocyte age sub-compartments retic1..retic10 are paper-
  # mechanistic states of the maturation-structured cytokinetic model
  # (Perez-Ruixo 2008 Methods, "Pharmacodynamic Analysis"). They are not
  # general-purpose canonical compartments, so they are declared as
  # paper-specific so checkModelConventions() does not flag them.
  paper_specific_compartment_pattern <- "^retic[0-9]+$"

  covariateData <- list(
    DOSE = list(
      description        = paste(
        "Per-subject planned single subcutaneous rHuEPO dose in IU.",
        "Used to evaluate the dose-dependent absolute bioavailability",
        "F = F0 + Emax(F) * DOSE / (ED50(F) + DOSE) (Perez-Ruixo 2008",
        "eq 5; Olsson-Gisleskog 2007). Encoded as a per-record covariate",
        "column because rxode2 zero-order infusion absorption with a",
        "podo()-dependent bioavailability evaluates podo() to NA inside",
        "f() when the dose record uses rate = -2 (model-defined dur());",
        "this DOSE covariate is the supported workaround for the same",
        "math, and it gives the user a single explicit place to declare",
        "the dose-dependence at simulation time."
      ),
      units              = "IU",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in single-dose use. Should equal the",
        "amt value of the depot and central dose records the user",
        "constructs for the dual SC absorption (the same Dose value",
        "appears in both the f() math and the event-table amt). For",
        "multi-dose simulations not foreseen by Perez-Ruixo 2008,",
        "DOSE should be set to the most recent dose. Studied range",
        "20000-160000 IU (= 20-160 kIU). Reference value not required",
        "because F is hyperbolic in DOSE, not power-form."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 88L,
    n_studies      = 3L,
    age_range      = "18-45 years",
    weight_range   = "63.6-100 kg",
    sex_female_pct = 0,
    race_ethnicity = "Not separately tabulated in Perez-Ruixo 2008.",
    disease_state  = paste(
      "Healthy adult male volunteers; physical exam, ECG and standard",
      "laboratory tests confirmed healthy. Baseline serum erythropoietin",
      "< 30 IU/L; haemoglobin 13.8-16.4 g/dL; haematocrit 41-49%; baseline",
      "reticulocytes <= 3%; normal serum folate, vitamin B12 and iron",
      "parameters. All subjects received daily oral iron supplementation",
      "(equivalent to 210 mg elemental iron) through day 29."
    ),
    dose_range     = paste(
      "Single subcutaneous rHuEPO dose. Study A (1996, n=20):",
      "300/600/1200/2400 IU/kg. Study B (1996, n=20):",
      "450/900/1350/1800 IU/kg. Study C (2002, n=48; n=47 included in",
      "PD analysis): 20/40/60/90/120/160 kIU absolute. Studies A and B",
      "doses-per-kg were converted to absolute doses for the analyses."
    ),
    regions        = paste(
      "Study A and B: South Florida Bioavailability Clinic (Miami, FL,",
      "USA). Study C: Chiltern International (Buckinghamshire, UK)."
    ),
    n_pk_observations = 2018L,
    n_pd_observations = 1628L,
    notes          = paste(
      "Three open-label, randomized, placebo-controlled, parallel-group",
      "phase I studies pooled. Serum rHuEPO by modified DSL",
      "radioimmunoassay (LLOQ 7.8 IU/L; CV < 15% up to 5000 IU/L).",
      "Percentage of reticulocytes by Technicon H-1 flow cytometry.",
      "PK analysis used MAP-Bayesian POSTHOC against the Olsson-Gisleskog",
      "2007 prior; PD analysis used FOCE on the PD layer with individual",
      "PK predictions as input. Demographic summary from Perez-Ruixo",
      "2008 Methods (Study Design and Subject Eligibility Criteria)."
    )
  )

  ini({
    # =================================================================
    # PK STRUCTURAL PARAMETERS
    # Olsson-Gisleskog 2007 prior typical values reported in Perez-Ruixo
    # 2008 Table I (Prior mean column). The MAP-Bayesian POSTHOC PK
    # estimation in NONMEM used these as priors; the column 'Prior mean
    # (% CV)' lists the Olsson-Gisleskog typical value and the
    # log-normal IIV inherited as a prior from the upstream popPK.
    # =================================================================
    lka     <- log(0.034);           label("First-order SC absorption rate ka (1/h)")            # Perez-Ruixo 2008 Table I prior mean
    lcl     <- log(0.358);           label("Linear clearance CLI (L/h)")                          # Perez-Ruixo 2008 Table I prior mean
    lvc     <- log(3.89);            label("Central volume V2 (L)")                               # Perez-Ruixo 2008 Table I prior mean
    lvmax   <- log(211);             label("Saturable elimination Vmax (IU/h)")                   # Perez-Ruixo 2008 Table I prior mean
    ld1     <- log(0.725);           label("Faster-pathway zero-order infusion duration D1 (h)")  # Perez-Ruixo 2008 Table I prior mean
    ltlag2  <- log(2.72);            label("Slower-pathway lag time tlag2 (h)")                   # Perez-Ruixo 2008 Table I prior mean
    lfdepot <- log(0.62);            label("Minimum absolute bioavailability F0 (fraction)")      # Perez-Ruixo 2008 Table I prior mean
    lrbase  <- log(13.9);            label("Baseline endogenous EPO concentration BSL (IU/L)")    # Perez-Ruixo 2008 Table I prior mean
    logitfr <- logit(0.60);          label("Fast/slow absorption split fraction fr (logit scale)") # Perez-Ruixo 2008 Table I prior mean

    # PK structural parameters held fixed from Olsson-Gisleskog 2007
    # (Perez-Ruixo 2008 Table I footnote a). All six were fixed during
    # the POSTHOC PK estimation; they carry no IIV.
    lkm     <- fixed(log(394));      label("Michaelis-Menten constant Km (IU/L)")                       # Table I footnote a: Km fixed at 394 IU/L
    lq      <- fixed(log(0.044));    label("Intercompartmental clearance Q (L/h)")                       # Table I footnote a: Q fixed at 0.044 L/h
    lvp     <- fixed(log(1.64));     label("Peripheral volume V3 (L)")                                   # Table I footnote a: V3 fixed at 1.64 L
    ld2     <- fixed(log(37.8));     label("Slower-pathway zero-order infusion duration D2 (h)")         # Table I footnote a: D2 fixed at 37.8 h
    lemax_f <- fixed(log(0.649));    label("Maximum increase in absolute bioavailability Emax(F) (fraction)")  # Table I footnote a: Emax(F) fixed at 64.9%
    led50_f <- fixed(log(63200));    label("Dose giving 50% of the bioavailability increase ED50(F) (IU)")     # Table I footnote a: ED50(F) fixed at 63.2 kIU = 63200 IU

    # =================================================================
    # PK IIVs (Perez-Ruixo 2008 Table I prior %CV column; Olsson-Gisleskog
    # 2007 log-normal IIVs; eq 26). Variance on the log scale uses
    # omega^2 = log(CV^2 + 1).
    # IIV on fr is on the logit scale (paper Methods, "In order to keep
    # the individual values of fr between 0 and 1. IIV in this parameter
    # was modelled using a normal distribution in the logit domain");
    # the CV = 45% in Table I is back-transformed; the logit-scale
    # variance is approximated here as log(CV^2 + 1) and the
    # approximation is documented in the vignette Errata.
    # =================================================================
    etalka     ~ 0.1206  # log(0.36^2 + 1) = 0.1206; Table I prior CV(ka)    = 36%
    etalcl     ~ 0.1029  # log(0.33^2 + 1) = 0.1029; Table I prior CV(CLI)   = 33%
    etalvc     ~ 0.0913  # log(0.31^2 + 1) = 0.0913; Table I prior CV(V2)    = 31%
    etalvmax   ~ 0.0808  # log(0.29^2 + 1) = 0.0808; Table I prior CV(Vmax)  = 29%
    etald1     ~ 0.9398  # log(1.25^2 + 1) = 0.9398; Table I prior CV(D1)    = 125%
    etaltlag2  ~ 0.2487  # log(0.53^2 + 1) = 0.2487; Table I prior CV(tlag2) = 53%
    etalfdepot ~ 0.1147  # log(0.35^2 + 1) = 0.1147; Table I prior CV(F0)    = 35%
    etalrbase  ~ 0.0862  # log(0.30^2 + 1) = 0.0862; Table I prior CV(BSL)   = 30%
    etalogitfr ~ 0.1854  # log(0.45^2 + 1) = 0.1854; Table I prior CV(fr)    = 45% (logit-domain approximation)

    # =================================================================
    # PD STRUCTURAL PARAMETERS (Perez-Ruixo 2008 Table II Model D,
    # Original-dataset column). Model D = stimulation of progenitor
    # production via SC50 + 5-stage signal-transduction inhibition of
    # reticulocyte aging via EC50 / tau. The maximum stimulation factor
    # of the signal (Smax) was not identifiable in the presence of the
    # S0/SM ratio (Methods); the driving function therefore reduces to
    # C/(EC50+C) (Smax fixed at 1).
    # Np = NR = 10 age compartments and M = 5 signal-transduction
    # transit compartments are fixed structural choices (Methods).
    # =================================================================
    lret0 <- log(1.24);   label("Baseline reticulocyte percentage RET0 (% reticulocytes)")        # Table II Model D, Original dataset
    ltr   <- log(62.2);   label("Mean lifespan of circulating reticulocytes TR (h)")              # Table II Model D, Original dataset
    ltp   <- log(118);    label("Mean lifespan of precursor cells Tp (h)")                        # Table II Model D, Original dataset
    lsc50 <- log(7.61);   label("EPO concentration for 50% of max progenitor production SC50 (IU/L)")  # Table II Model D, Original dataset
    lec50 <- log(56.3);   label("EPO concentration for 50% of max signal EC50 (IU/L)")            # Table II Model D, Original dataset
    ltau  <- log(4.89);   label("Signal transduction transit time tau (h)")                       # Table II Model D, Original dataset

    # =================================================================
    # PD IIVs (Perez-Ruixo 2008 Table II Model D, IIV (omega) [%] block,
    # Original-dataset column). Log-normal IIV per eq 26;
    # omega^2 = log(CV^2 + 1). EC50 and tau have no IIV reported
    # (blank in the table) and are therefore typical-value only.
    # The SC50 IIV is large (107% CV; Table II footnote a references
    # SC50P, but Model D's stimulation term uses Model A's SC50 per
    # Methods 'The stimulation of kin was described as in model A as
    # well as equations for P1, ..., PNp.').
    # =================================================================
    etalret0 ~ 0.0386  # log(0.198^2 + 1) = 0.0386; Table II Model D, IIV(RET0) = 19.8%
    etaltr   ~ 0.1224  # log(0.362^2 + 1) = 0.1224; Table II Model D, IIV(TR)   = 36.2%
    etaltp   ~ 0.0723  # log(0.274^2 + 1) = 0.0723; Table II Model D, IIV(Tp)   = 27.4%
    etalsc50 ~ 0.7634  # log(1.07^2 + 1)  = 0.7634; Table II Model D, IIV(SC50) = 107%

    # =================================================================
    # RESIDUAL ERROR
    # PD residual: Perez-Ruixo 2008 eq 27 (Robs = Rpred + epsilon,
    # epsilon ~ N(0, sigma^2)). Table II Model D sigma = 63.1% (Original
    # dataset). The "[%]" column header is interpreted as a proportional
    # CV on the predicted reticulocyte percentage (the alternative
    # additive interpretation of 0.631 absolute % reticulocytes is also
    # consistent with the reported number; the proportional reading is
    # used here, documented in the vignette Errata).
    # PK residual: NOT reported in Perez-Ruixo 2008. The PK residual
    # lives in Olsson-Gisleskog 2007 (the upstream popPK paper, not on
    # disk for this extraction). Encoded as fixed(0) per the task hard
    # constraint ("If a needed value is unreported, encode fixed(0) +
    # an erratum rather than a class-typical placeholder"); documented
    # in the vignette Errata.
    # =================================================================
    propSd     <- fixed(0); label("Proportional residual error on Cc (placeholder; not reported in Perez-Ruixo 2008)")
    propSd_RET <- 0.631;    label("Proportional residual error on RET (fraction)")  # Table II Model D, sigma = 63.1%
  })

  model({
    # ---------------------------------------------------------------
    # 1. Individual structural parameters
    # ---------------------------------------------------------------
    ka     <- exp(lka     + etalka)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc     + etalvc)
    vmax   <- exp(lvmax   + etalvmax)
    d1_dur <- exp(ld1     + etald1)
    tlag2  <- exp(ltlag2  + etaltlag2)
    f0     <- exp(lfdepot + etalfdepot)
    bsl    <- exp(lrbase  + etalrbase)
    fr     <- expit(logitfr + etalogitfr)

    km     <- exp(lkm)
    q      <- exp(lq)
    vp     <- exp(lvp)
    d2_dur <- exp(ld2)
    emax_f <- exp(lemax_f)
    ed50_f <- exp(led50_f)

    ret0 <- exp(lret0 + etalret0)
    tr   <- exp(ltr   + etaltr)
    tp   <- exp(ltp   + etaltp)
    sc50 <- exp(lsc50 + etalsc50)
    ec50 <- exp(lec50)
    tau  <- exp(ltau)

    # ---------------------------------------------------------------
    # 2. Derived structural quantities
    # ---------------------------------------------------------------
    # Fixed integer structural choices (Methods).
    np <- 10
    nr <- 10

    # Endogenous EPO production rate, derived from the steady-state
    # balance against linear + Michaelis-Menten elimination at C = BSL
    # (Perez-Ruixo 2008 eq 4). kEPO is a constant input to central; the
    # accompanying central(0) = BSL * vc and peripheral1(0) = BSL * vp
    # initial conditions place the PK system at the endogenous SS at
    # t = 0 with no exogenous rHuEPO administered yet.
    kEPO <- cl * bsl + vmax * bsl / (km + bsl)

    # Maximal endogenous progenitor production rate, derived from the
    # steady-state baseline reticulocyte count RET0 (Perez-Ruixo 2008
    # eq 13). At C = BSL, production = kin * BSL / (SC50 + BSL) and the
    # zero-net-flux closure on the reticulocyte chain forces this to
    # equal RET0 / TR.
    kin <- (ret0 / tr) * (sc50 + bsl) / bsl

    # Baseline signal level: the steady-state of the signal-transduction
    # chain at C = BSL (Perez-Ruixo 2008 eq 23 implicit form). All five
    # transit compartments start at this value; the modulator S0/SM
    # equals 1 at baseline.
    sig0 <- bsl / (ec50 + bsl)

    # Per-compartment first-order transfer rates.
    k_p      <- np / tp           # precursor age transfer (h^-1)
    k_r      <- nr / tr           # baseline reticulocyte age transfer (h^-1)
    rate_sig <- 1 / tau           # signal transit rate (h^-1)

    # ---------------------------------------------------------------
    # 3. PK ODEs (Perez-Ruixo 2008 eqs 1-3). Two-compartment with linear
    # + Michaelis-Menten elimination, first-order ka absorption from
    # depot to central, and a constant endogenous production kEPO into
    # central. The dual-pathway SC absorption is realised at the dose
    # event level: dose into depot (cmt = depot, rate = -2, duration set
    # by dur(depot) = D1, bioavailability F * fr); simultaneous dose
    # into central (cmt = central, rate = -2, duration set by
    # dur(central) = D2, lag alag(central) = tlag2, bioavailability
    # F * (1 - fr)). Both pathways share the same dose-dependent
    # bioavailability F = F0 + Emax(F) * Dose / (ED50(F) + Dose).
    # ---------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + kEPO -
                          cl * central / vc -
                          vmax * central / (km * vc + central) -
                          q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <-  q  * central / vc - q * peripheral1 / vp

    # Dose-dependent absolute bioavailability (Perez-Ruixo 2008 eq 5).
    # DOSE comes from the event-table covariate column (declared in
    # covariateData$DOSE). The rxode2 combination of zero-order infusion
    # absorption (rate = -2) and podo()-dependent bioavailability inside
    # f() evaluates podo() to NA at the dose record, so the user-supplied
    # DOSE covariate is the supported encoding. The same DOSE value should
    # also be the amt of the depot and central dose records.
    #
    # Dual-pathway dose-record convention (see validation vignette): the
    # depot dose record uses time = t_dose; the central dose record uses
    # time = t_dose + dt (any positive offset, e.g. 1e-3 h) so the two
    # simultaneous rate = -2 records do not collide inside rxode2's
    # solver. alag(central) = tlag2 then shifts the actual administration
    # of the central dose to t_dose + dt + tlag2 ~ t_dose + tlag2 -- the
    # paper's intended slow-pathway delay -- without changing the per-
    # subject IIV semantics (tlag2 still carries etaltlag2).
    f_total <- f0 + emax_f * DOSE / (ed50_f + DOSE)
    f(depot)      <- f_total * fr
    f(central)    <- f_total * (1 - fr)
    dur(depot)    <- d1_dur
    dur(central)  <- d2_dur
    alag(central) <- tlag2

    # ---------------------------------------------------------------
    # 4. Observable: total EPO serum concentration (endogenous +
    # exogenous rHuEPO). With central holding the steady-state
    # endogenous amount BSL * vc at t = 0 and kEPO maintaining it,
    # central/vc reads the total observable concentration directly.
    # ---------------------------------------------------------------
    Cc <- central / vc

    # ---------------------------------------------------------------
    # 5. Signal transduction chain (Perez-Ruixo 2008 eqs 21-22; M = 5
    # transit compartments). The driving function C/(EC50+C) is bounded
    # in [0, 1]; Smax is implicitly fixed at 1 per Methods because the
    # ratio S0/SM in the reticulocyte equations renders Smax
    # unidentifiable.
    # ---------------------------------------------------------------
    d/dt(transit1) <- rate_sig * (Cc / (ec50 + Cc) - transit1)
    d/dt(transit2) <- rate_sig * (transit1 - transit2)
    d/dt(transit3) <- rate_sig * (transit2 - transit3)
    d/dt(transit4) <- rate_sig * (transit3 - transit4)
    d/dt(transit5) <- rate_sig * (transit4 - transit5)

    # Reticulocyte aging modulator. At baseline transit5 = sig0 so
    # aging_mod = 1 and the chain runs at its baseline rate NR/TR;
    # rising drug exposure pushes transit5 above sig0, shrinks
    # aging_mod, and slows reticulocyte aging (longer residence time).
    aging_mod <- sig0 / transit5

    # ---------------------------------------------------------------
    # 6. Bone-marrow precursor age chain (10 compartments;
    # Perez-Ruixo 2008 eqs 6-7). rHuEPO stimulates progenitor
    # production via the Emax-like SC50 term in d/dt(precursor1).
    # ---------------------------------------------------------------
    d/dt(precursor1)  <- kin * Cc / (sc50 + Cc) - k_p * precursor1
    d/dt(precursor2)  <- k_p * (precursor1  - precursor2)
    d/dt(precursor3)  <- k_p * (precursor2  - precursor3)
    d/dt(precursor4)  <- k_p * (precursor3  - precursor4)
    d/dt(precursor5)  <- k_p * (precursor4  - precursor5)
    d/dt(precursor6)  <- k_p * (precursor5  - precursor6)
    d/dt(precursor7)  <- k_p * (precursor6  - precursor7)
    d/dt(precursor8)  <- k_p * (precursor7  - precursor8)
    d/dt(precursor9)  <- k_p * (precursor8  - precursor9)
    d/dt(precursor10) <- k_p * (precursor9  - precursor10)

    # ---------------------------------------------------------------
    # 7. Circulating reticulocyte age chain (10 compartments;
    # Perez-Ruixo 2008 eqs 24-25 Model D). The first compartment is
    # fed by the terminal precursor at the standard rate k_p; the
    # transfer between reticulocyte compartments is modulated by the
    # signal ratio aging_mod = S0/SM.
    # ---------------------------------------------------------------
    d/dt(retic1)  <- k_p * precursor10 - k_r * aging_mod * retic1
    d/dt(retic2)  <- k_r * aging_mod * (retic1 - retic2)
    d/dt(retic3)  <- k_r * aging_mod * (retic2 - retic3)
    d/dt(retic4)  <- k_r * aging_mod * (retic3 - retic4)
    d/dt(retic5)  <- k_r * aging_mod * (retic4 - retic5)
    d/dt(retic6)  <- k_r * aging_mod * (retic5 - retic6)
    d/dt(retic7)  <- k_r * aging_mod * (retic6 - retic7)
    d/dt(retic8)  <- k_r * aging_mod * (retic7 - retic8)
    d/dt(retic9)  <- k_r * aging_mod * (retic8 - retic9)
    d/dt(retic10) <- k_r * aging_mod * (retic9 - retic10)

    # ---------------------------------------------------------------
    # 8. PD observable: total reticulocyte percentage (Perez-Ruixo
    # 2008 eq 10). Units carry through as %: each retic_j(0) = RET0/NR
    # in %; the chain preserves the unit; the sum is the reticulocyte
    # percentage in the same units as RET0 (e.g., 1.24%).
    # ---------------------------------------------------------------
    RET <- retic1 + retic2 + retic3 + retic4 + retic5 +
           retic6 + retic7 + retic8 + retic9 + retic10

    # ---------------------------------------------------------------
    # 9. Initial conditions (Perez-Ruixo 2008 eqs 11-12 and the implicit
    # endogenous PK steady state). The PD chains are at their pre-dose
    # steady state under the baseline endogenous EPO concentration BSL,
    # so kin = baseline-production satisfies the closure NR/TR * R_ss =
    # NR/TR * RET0/NR = RET0/TR = Np/Tp * P_ss.
    # ---------------------------------------------------------------
    central(0)     <- bsl * vc
    peripheral1(0) <- bsl * vp

    transit1(0) <- sig0
    transit2(0) <- sig0
    transit3(0) <- sig0
    transit4(0) <- sig0
    transit5(0) <- sig0

    precursor1(0)  <- (tp / np) * (ret0 / tr)
    precursor2(0)  <- (tp / np) * (ret0 / tr)
    precursor3(0)  <- (tp / np) * (ret0 / tr)
    precursor4(0)  <- (tp / np) * (ret0 / tr)
    precursor5(0)  <- (tp / np) * (ret0 / tr)
    precursor6(0)  <- (tp / np) * (ret0 / tr)
    precursor7(0)  <- (tp / np) * (ret0 / tr)
    precursor8(0)  <- (tp / np) * (ret0 / tr)
    precursor9(0)  <- (tp / np) * (ret0 / tr)
    precursor10(0) <- (tp / np) * (ret0 / tr)

    retic1(0)  <- ret0 / nr
    retic2(0)  <- ret0 / nr
    retic3(0)  <- ret0 / nr
    retic4(0)  <- ret0 / nr
    retic5(0)  <- ret0 / nr
    retic6(0)  <- ret0 / nr
    retic7(0)  <- ret0 / nr
    retic8(0)  <- ret0 / nr
    retic9(0)  <- ret0 / nr
    retic10(0) <- ret0 / nr

    # ---------------------------------------------------------------
    # 10. Per-output residual error.
    # ---------------------------------------------------------------
    Cc  ~ prop(propSd)
    RET ~ prop(propSd_RET)
  })
}
