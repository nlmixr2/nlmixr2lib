Chan_2008_maraviroc <- function() {
  description <- "Two-compartment population PK meta-analysis model for oral maraviroc (CCR5 antagonist) in healthy volunteers and asymptomatic HIV-infected adults, with hepatic-extraction-ratio parameterisation of clearance, dose-dependent absorption (sigmoid-Emax F_ABS and power-function ka), food effect on both, Asian-race covariates on hepatic extraction / peripheral volume / inter-compartmental clearance, an age effect on Q, and a TAD-dependent residual error (Chan 2008)"
  reference   <- "Chan PLS, Weatherley B, McFadyen L. A population pharmacokinetic meta-analysis of maraviroc in healthy volunteers and asymptomatic HIV-infected subjects. Br J Clin Pharmacol. 2008 Apr;65 Suppl 1:76-85. doi:10.1111/j.1365-2125.2008.03139.x"
  vignette    <- "Chan_2008_maraviroc"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL_ic, normalised to the median age of 30 years (Chan 2008 Methods, 'Covariate testing': 'Continuous covariates, age and weight, were modelled as multiplicative effects and normalized to their median values, 30 years and 71 kg, respectively.'). Implementation form: CL_ic = CL_ic_typ * (AGE / 30)^e_age_q. At age 60 the multiplier is (60/30)^0.349 = 1.274 (CL_ic is 27.4% greater than at age 30; Chan 2008 Results, 'The population estimate of CLic').",
      source_name        = "AGE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all non-Asian, including White, Black, and Other)",
      notes              = "Multiplicative-exponential covariate on E_H, V_p, and CL_ic: PARAM_Asian = PARAM_typ * exp(theta * RACE_ASIAN). Verified by back-calculation from paper Table 5 (typical Asian CL = 47.88 L/h; F_HEP_Asian = 0.398, matching CL_H = E_H * FQ = 0.662 * exp(-0.0948) * 59.59 = 35.89 L/h and CL_total = 35.89 + 12 = 47.89 L/h). Chan 2008 Methods, 'Covariate testing': 'In the final covariate search the race effect was tested as binary (Asians vs reference of all non-Asians).' Black subjects (3.4% of the cohort) were collapsed into the non-Asian reference because the preliminary Black-vs-Whites/Others test was not statistically significant.",
      source_name        = "ASIAN"
    ),
    FED = list(
      description        = "Fed-vs-fasted at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-dose-event indicator. 1 = oral tablet co-administered with a high-fat meal; 0 = fully fasted overnight with meals deferred at least 4 h post-dose. Drives a factorial multiplicative effect on ka (Chan 2008 Eq 6) and a log-multiplicative-exponential effect on ABSEmax / ED50 (Eqs 7-8). Light-meal and intermediate-timing food conditions were excluded from the analysis dataset.",
      source_name        = "FED"
    ),
    DOSE = list(
      description        = "Maraviroc tablet dose at the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-dose-record continuous covariate carrying the administered dose level in mg. Used by the dose-dependent sigmoid-Emax F_ABS (Chan 2008 Eq 4) and by the dose-dependent power-function ka (Eq 5). Supply DOSE at the same value as the rxode2 amt column on each dose record and carry it forward by last-observation-carried-forward (LOCF) so ka and F_ABS reflect the most recent dose during ongoing absorption.",
      source_name        = "Dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 413L,
    n_studies      = 17L,
    age_range      = "18-54 years (median 30)",
    age_median     = "30 years",
    weight_range   = "46-109 kg (median 71)",
    weight_median  = "71 kg",
    sex_female_pct = 23.2,
    race_ethnicity = c(White = 73.1, Asian = 23.0, Black = 3.4, Other = 0.5),
    disease_state  = "Pooled healthy volunteers (n = 365, 88.4%) and asymptomatic HIV-infected subjects (n = 48, 11.6%)",
    dose_range     = "Per-dose 25-1200 mg oral tablet; total daily doses 100-1800 mg/day; single- and multiple-dose regimens",
    regions        = "Multi-national; Pfizer Phase 1 and 2a clinical pharmacology studies",
    n_observations = 8951L,
    notes          = "Meta-analysis of 17 Pfizer studies (Chan 2008 Table 1: A4001003 through A4001043). 690 concentration-time profiles; 46 (6.7%) profiles obtained under fed conditions. HIV-infected subjects contributed 929 of the 8951 plasma observations. Demographic summaries from Chan 2008 Tables 2-3."
  )

  ini({
    # Structural parameters at the reference covariate set:
    #   AGE = 30 y, RACE_ASIAN = 0 (non-Asian), FED = 0 (fasted), DOSE = 100 mg.
    #
    # Hepatic extraction ratio (E_H) is bounded in [0, 1]. It is parameterised
    # on the logit scale so the intersubject-variability transformation in
    # Chan 2008 Eq 9 -- logit(E_H_i) = logit(E_H_typ) + eta -- reduces to a
    # standard additive eta on the logit scale. logit(0.662) = 0.6722.
    logiteh <- 0.6722       ; label("Logit of hepatic extraction ratio E_H at the reference covariate set (unitless)")   # Chan 2008 Table 4: E_H = 0.662 (RSE 1.5%)
    lvc     <- log(132)     ; label("Central volume of distribution (V_c, L)")                                          # Chan 2008 Table 4: V_c = 132 L (RSE 2.7%)
    lvp     <- log(277)     ; label("Peripheral volume of distribution (V_p, L) at the non-Asian reference")            # Chan 2008 Table 4: V_p = 277 L (RSE 4.2%)
    lq      <- log(16.4)    ; label("Inter-compartmental clearance (CL_ic, L/h) at the 30-y non-Asian reference")       # Chan 2008 Table 4: CL_ic = 16.4 L/h (RSE 3.9%)
    lka     <- log(0.277)   ; label("Absorption rate constant at 1-mg dose under fasted conditions (ka_1mg, 1/h)")      # Chan 2008 Table 4: ka_1mg = 0.277 /h (RSE 21.0%)
    ltlag   <- log(0.198)   ; label("Absorption lag time (Tlag, h)")                                                    # Chan 2008 Table 4: Tlag = 0.198 h (RSE 4.1%)
    led50   <- log(51.2)    ; label("Dose producing 50% of ABSEmax at the fasted reference (ED50, mg)")                 # Chan 2008 Table 4: ED50 = 51.2 mg (RSE 12.8%)
    lhill   <- log(1.39)    ; label("Sigmoidicity exponent g of the dose-vs-F_ABS sigmoid-Emax model (unitless)")       # Chan 2008 Table 4: g = 1.39 (RSE 15.3%)

    # Dose / food effect coefficients on absorption
    qka          <- 0.173   ; label("Dose-power exponent on ka in Chan 2008 Eq 5 (unitless)")                           # Chan 2008 Table 4: theta_ka (dose exp) = 0.173 (RSE 22.9%)
    thetaka1     <- 0.547   ; label("Multiplicative food-vs-fasted ratio applied to ka in Eq 6 (unitless)")             # Chan 2008 Table 4: theta_ka1 = 0.547 (RSE 18.6%); 45.3% reduction under high-fat meal
    thetaABSEmax <- -0.258  ; label("Log multiplicative high-fat-meal effect on ABSEmax in Eq 7 (unitless)")            # Chan 2008 Table 4: theta_ABSEmax = -0.258 (RSE 27.3%)
    thetaED50    <- 0.594   ; label("Log multiplicative high-fat-meal effect on ED50 in Eq 8 (unitless)")               # Chan 2008 Table 4: theta_ED50 = 0.594 (RSE 34.5%)

    # Race covariate effects (multiplicative-exponential: PARAM_i = PARAM_typ * exp(theta * RACE_ASIAN));
    # age covariate effect (power form: CL_ic = CL_ic_typ * (AGE / 30)^e_age_q).
    e_race_asian_eh <- -0.0948 ; label("Multiplicative-exponential Asian-vs-non-Asian effect on E_H (unitless)")        # Chan 2008 Table 4: Race effect on E_H = -0.0948 (RSE 15.0%); verified by back-calc from Table 5 F_HEP_Asian = 0.398
    e_race_asian_vp <- -0.637  ; label("Multiplicative-exponential Asian-vs-non-Asian effect on V_p (unitless)")        # Chan 2008 Table 4: Race effect on V_p = -0.637 (RSE 8.5%)
    e_race_asian_q  <- -0.298  ; label("Multiplicative-exponential Asian-vs-non-Asian effect on CL_ic (unitless)")      # Chan 2008 Table 4: Race effect on CL_ic = -0.298 (RSE 17.7%)
    e_age_q         <- 0.349   ; label("Power exponent of AGE on CL_ic (reference 30 years; unitless)")                 # Chan 2008 Table 4: Age effect on CL_ic = 0.349 (RSE 28.3%); verified by back-calc CL_ic(60)/CL_ic(30) = (60/30)^0.349 = 1.274

    # Fixed structural anchors (Chan 2008 Table 4 and Methods)
    absemax <- fixed(1)     ; label("Maximum F_ABS (Eq 4; fixed by the authors)")                                       # Chan 2008 Table 4: ABSEmax = 1 FIX
    fq      <- fixed(59.59) ; label("Hepatic plasma flow used in Eqs 2-3 (L/h; fixed)")                                 # Chan 2008 Table 4: FQ = 59.59 L/h FIX (from a prior mass-balance maraviroc analysis)
    clr     <- fixed(12)    ; label("Renal clearance (L/h; fixed)")                                                     # Chan 2008 Methods: CL_R fixed at 12 L/h (taken from a previous mass-balance maraviroc analysis)

    # TAD-dependent residual error (Chan 2008 Eq 12). The error is additive on
    # log-transformed concentration: ln(Y) = ln(Cc) + W(TAD) * eps with var(eps)
    # fixed at 1, so W(TAD) maps to nlmixr2's lnorm-residual SD on log scale.
    # Per the paper, Pmax and Base are 'expressed as CV%' but are the log-scale
    # SDs themselves (the standard small-omega NONMEM convention CV% ~ omega *
    # 100 on log-transformed data).
    pmax_res <- 0.742 ; label("Peak residual SD above baseline (Chan 2008 Eq 12; log-scale SD ~ CV)")                   # Chan 2008 Table 4: Pmax = 74.2% (RSE 3.9%)
    tmax_res <- 0.950 ; label("Time of peak residual SD after a dose (Chan 2008 Eq 12; h)")                             # Chan 2008 Table 4: Tmax = 0.950 h (RSE 8.8%)
    k_res    <- 0.403 ; label("Exponential decay-rate constant of the TAD-dependent residual SD (Chan 2008 Eq 12; 1/h)") # Chan 2008 Table 4: K = 0.403 1/h (RSE 6.4%)
    base_res <- 0.202 ; label("Baseline residual SD (Chan 2008 Eq 12; log-scale SD ~ CV)")                              # Chan 2008 Table 4: Base = 20.2% (RSE 3.3%)

    # Inter-subject variability (Chan 2008 Table 4). For the multiplicative-
    # exponential parameters the reported approximate %CV is omega * 100 on the
    # log scale, so omega^2 = (CV / 100)^2. For E_H the reported 8.3% is
    # (1 - E_H) * sqrt(omega^2) per Table 4 footnote -- the underlying logit-
    # scale variance is (0.083 / (1 - 0.662))^2 = 0.0603.
    etalogiteh ~ 0.0603    # Chan 2008 Table 4: w[E_H] = 8.3% with (1 - 0.662) * sqrt(omega^2) -> omega^2 = (0.083 / 0.338)^2 = 0.0603 on the logit(E_H) scale
    etaled50   ~ 0.347     # Chan 2008 Table 4: w[ED50] = 58.9% -> omega^2 = 0.589^2 = 0.347
    etalvc     ~ 0.0132    # Chan 2008 Table 4: w[V_c]  = 11.5% -> omega^2 = 0.115^2 = 0.0132
    etalq      ~ 0.0930    # Chan 2008 Table 4: w[CL_ic]= 30.5% -> omega^2 = 0.305^2 = 0.0930
    etalka     ~ 0.160     # Chan 2008 Table 4: w[ka]   = 40.0% -> omega^2 = 0.400^2 = 0.160
    etalvp     ~ 0.0773    # Chan 2008 Table 4: w[V_p]  = 27.8% -> omega^2 = 0.278^2 = 0.0773
  })

  model({
    # ---- Hepatic extraction ratio with race covariate and logit-scale IIV ----
    # Population typical E_H (no covariate, no eta):
    eh_pop_typ <- expit(logiteh)
    # Asian-race effect applied as multiplicative-exponential on linear E_H scale
    # (Chan 2008 Table 4 and Results: typical Asian E_H = 0.662 * exp(-0.0948) =
    # 0.6022; F_HEP_Asian = 1 - 0.6022 = 0.3978, matching Table 5 'F_HEP Asian =
    # 0.398' and CL_Asian = 35.89 + 12 = 47.89 L/h matching Table 5 47.88).
    eh_typ_cov <- eh_pop_typ * exp(e_race_asian_eh * RACE_ASIAN)
    # IIV via Eq 9: logit(E_H_i) = logit(eh_typ_cov) + etalogiteh.
    # Compact algebraic form using the odds-ratio R = eh_typ_cov / (1 - eh_typ_cov):
    eh_ratio <- eh_typ_cov / (1 - eh_typ_cov) * exp(etalogiteh)
    eh       <- eh_ratio / (1 + eh_ratio)

    # ---- Clearance decomposition (Chan 2008 Eqs 2-3 + text) ----
    clh  <- fq * eh        # Hepatic CL (Eq 2)
    cl   <- clh + clr      # Total CL = CL_H + CL_R (Methods text)
    fhep <- 1 - eh         # F_HEP = 1 - E_H (Eq 3)

    # ---- Volumes and inter-compartmental CL with race / age covariates ----
    # Race effects are multiplicative-exponential on linear PARAM scale; age
    # effect on CL_ic is a power form normalised to the median 30 y.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * exp(e_race_asian_vp * RACE_ASIAN)
    q  <- exp(lq  + etalq)  *
            exp(e_race_asian_q * RACE_ASIAN) *
            (AGE / 30)^e_age_q

    # ---- Absorption lag ----
    tlag_depot <- exp(ltlag)

    # ---- Dose- and food-dependent ka (Chan 2008 Eqs 5-6) ----
    # ka(dose, fed) = ka_1mg * DOSE^qka * (theta_ka1 if FED else 1)
    ka_1mg <- exp(lka + etalka)
    ka     <- ka_1mg * DOSE^qka * (FED * thetaka1 + (1 - FED))

    # ---- Dose-dependent F_ABS with food effects on ABSEmax and ED50 (Eqs 4, 7, 8) ----
    ed50_typ    <- exp(led50 + etaled50)
    hill        <- exp(lhill)
    absemax_eff <- absemax * exp(thetaABSEmax * FED)
    ed50_eff    <- ed50_typ * exp(thetaED50 * FED)
    fabs        <- absemax_eff * DOSE^hill / (ed50_eff^hill + DOSE^hill)

    # ---- Total bioavailability F = F_HEP * F_ABS (Chan 2008 Eq 1) ----
    fdepot <- fhep * fabs

    # ---- Micro-constants for the two-compartment disposition ODE ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: depot -> central -> peripheral1 ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Bioavailability and lag applied at dose event on the depot ----
    f(depot)    <- fdepot
    alag(depot) <- tlag_depot

    # ---- Concentration: dose mg / volume L = mg/L; multiply by 1000 for ng/mL
    # to match the paper's LLOQ of 0.5 ng/mL and Table 5 AUC units of ng h /mL.
    Cc <- central / vc * 1000

    # ---- TAD-dependent residual error (Chan 2008 Eq 12) ----
    # W(TAD) = Pmax * A * TAD^P * exp(-K * TAD) + Base, where P = K * Tmax and
    # A = exp(P) / Tmax^P. At TAD = Tmax the peak value is Base + Pmax.
    # Verification (paper Results): at TAD = 6 h W = 39.8% (paper: 'dropped to
    # 40%'); at TAD = 10 h W = 25.0% (paper: '25%'). Applied via lnorm() so
    # log(Cc) = log(Cc_pred) + W * eps, matching the paper's log-transformed
    # additive error structure.
    tad_safe <- max(tad(), 0)
    ppow     <- k_res * tmax_res
    a_res    <- exp(ppow) / (tmax_res^ppow)
    expSd    <- pmax_res * a_res * tad_safe^ppow * exp(-k_res * tad_safe) + base_res
    Cc      ~ lnorm(expSd)
  })
}
