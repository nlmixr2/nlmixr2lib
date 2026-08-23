Yin_2023_soticlestat <- function() {
  description <- paste(
    "Joint population PK / CH24H enzyme-occupancy (EO) / 24S-hydroxycholesterol (24HC)",
    "pharmacodynamic model for soticlestat (TAK-935), a first-in-class selective inhibitor of",
    "cholesterol 24-hydroxylase (CH24H / CYP46A1), in healthy adults (Yin 2023).",
    "Two-compartment PK with first-order absorption fed by an explicit two-compartment transit",
    "chain that is the dosing route for the tablet formulation, the oral solution dosing",
    "directly into the absorption depot; empirical dose-nonlinearity power terms on CL/F, Q/F",
    "and Vp/F; a plasma-to-brain effect compartment whose concentration drives both a sigmoid",
    "Emax brain CH24H enzyme-occupancy read-out and a semimechanistic sigmoid Imax inhibitory",
    "indirect-response turnover model for plasma 24HC.",
    "Weight-based allometry on CL/F, Vc/F, Q/F and Vp/F (reference 70 kg) was added by the",
    "authors after estimation to support the paediatric dosing simulations."
  )
  reference <- paste(
    "Yin W, Facius A, Wagner T, Tsai M, Asgharnejad M, Lahu G, Vakilynejad M.",
    "Population pharmacokinetics, enzyme occupancy, and 24S-hydroxycholesterol modeling of",
    "soticlestat, a novel cholesterol 24-hydroxylase inhibitor, in healthy adults.",
    "Clin Transl Sci. 2023;16(7):1149-1162. doi:10.1111/cts.13517.",
    "Structural equations from Appendix S1 (NONMEM control streams, three sections) and",
    "Appendix S2 (population PK differential equations); all shipped values from Table 2.",
    "The CH24H enzyme-occupancy sub-model estimated here is carried forward FIXED into the",
    "patient model of Yin W, Facius A, Asgharnejad M, Lahu G, Vakilynejad M.",
    "Population pharmacokinetics, enzyme occupancy, and pharmacodynamic modeling of",
    "soticlestat in patients with developmental and epileptic encephalopathies.",
    "Clin Transl Sci. 2024;17(3):e13722. doi:10.1111/cts.13722; see",
    "modellib('Yin_2024_soticlestat').",
    sep = " "
  )
  vignette <- "Yin_2023_soticlestat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  dosing <- c("depot", "transit1")

  covariateData <- list(
    DOSE = list(
      description        = "Soticlestat dose level of the current administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Empirical dose-nonlinearity covariate entering as (DOSE/300)^exponent on the TYPICAL",
        "VALUES of CL/F, Q/F and Vp/F (NONMEM $PK 'TVCL = THETA(4)*((DOSE/300.00)**THETA(9))'",
        "and its Q and Vp analogues; Table 2 'Dose effect (exponent)' rows). Reference 300 mg.",
        "DOSE is the PER-ADMINISTRATION dose, not the total daily dose: reading it as the daily",
        "dose overpredicts every published steady-state AUC answer key by 15-19%, whereas the",
        "per-administration reading reproduces all eleven Table S3 rows within 7.5%.",
        "All three exponents were ESTIMATED, not fixed (Table 2 reports an RSE and a bootstrap",
        "95% CI for each). Studied dose range 15-1350 mg (Table 1).",
        "No formal covariate screening was performed; dose was tested specifically to",
        "characterise the observed PK nonlinearity (paper, 'Covariate model development')."
      ),
      source_name        = "DOSE"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power (allometric) effect on CL/F, Vc/F, Q/F and Vp/F, reference 70 kg. Weight was",
        "NOT a covariate of the estimated model: 'Weight was not tested as a covariate because",
        "of the limited range in the analysis data set; however, weight was added to the model",
        "prior to the pediatric simulations using power functions with weight centered at",
        "70 kg' (paper, 'PopPK model'). The exponents are therefore fixed structural",
        "assumptions of the paediatric simulation, not estimates, and all four are wrapped in",
        "fixed(). Setting WT = 70 recovers the estimated adult model exactly.",
        "See the vignette Errata for the printed Vc/F sign and for the two exponents the paper",
        "did not print."
      ),
      source_name        = "WEIGHT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Collected at baseline across the four studies (paper, 'Study population and data",
        "collection'; Table S1 reports mean (SD) 34.7 (9.6) years) but not retained in any",
        "layer of the final model. The paediatric simulations explicitly assume no maturation",
        "effect: 'No differences in maturation stage were expected because the target",
        "population was at least 2 years of age.'"
      )
    ),
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Collected at baseline (Table S1, mean (SD) 25.6 (2.9) kg/m^2); not retained in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected at baseline (Table S1; 74/108 (69%) male); not retained in the final model."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected at baseline (Table S1; 78/108 (72%) White); not retained in the final model."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected at baseline (Table S1; 28/108 (26%) Black or African American); not retained in the final model."
    ),
    FORM_TABLET = list(
      description = "Tablet (vs oral solution) formulation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "The formulation difference IS in the model, but as a routing choice rather than a",
        "covariate: the tablet doses into transit1 and the oral solution doses into depot",
        "(Appendix S2: 'The transit compartment 1 is the dosing compartment for the tablet",
        "formulation'). Routing is expressed through the event table's cmt column, so no",
        "formulation covariate column is referenced in model()."
      )
    )
  )

  compartmentData <- list(
    transit1 = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte  = "soticlestat",
      units    = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    hc24 = list(
      analyte  = "24S-hydroxycholesterol",
      units    = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104,
    n_studies      = 4,
    age_range      = "19-55 years across the four studies; mean (SD) 34.7 (9.6) years",
    weight_range   = "not reported; mean (SD) body mass index 25.6 (2.9) kg/m^2",
    sex_female_pct = 31,
    race_ethnicity = c(White = 72, `Black or African American` = 26, Multiracial = 2),
    disease_state  = "healthy adults",
    dose_range     = "15-1350 mg single oral dose, or 100-600 mg q.d. and 300 mg b.i.d. for 10-14 days",
    notes          = paste(
      "Four phase I studies in healthy adults (Table 1): a single-rising-dose study",
      "(NCT02201056, n = 48), a multiple-rising-dose study (NCT02539134, n = 40), an",
      "open-label [18F]MNI-792 PET enzyme-occupancy study (NCT02497235, n = 11), and a",
      "relative-bioavailability and food-effect crossover study (NCT02906813, n = 9);",
      "108 participants enrolled in total (Table S1).",
      "The three model layers used three overlapping analysis sets: the population PK model",
      "1727 soticlestat concentrations from 104 individuals, the PK/EO model 20 brain CH24H",
      "occupancy observations from 11 individuals, and the PK/PD model 2270 plasma 24HC",
      "concentrations from 99 individuals. No data were excluded or imputed.",
      "n_subjects records the population PK analysis set.",
      "Soticlestat was given as an oral solution in every study except the food-effect study,",
      "which used both the oral solution and a 100 mg tablet. All dosing was oral; no",
      "intravenous data were available, so absolute bioavailability could not be estimated",
      "and every clearance and volume is an apparent (/F) parameter.",
      "Race percentages are from the Table S1 total column and are rounded."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Population PK model -- Table 2, 'Population PK model' panel,
    # 'Original / Estimate' column. The Appendix S1 $THETA / $OMEGA blocks
    # carry INITIAL estimates (CL 220, Vc 60, ka 2.0) and are used here only
    # for structure and for the FIX flags.
    # ---------------------------------------------------------------------
    lka  <- log(2.13)  ; label("First-order absorption rate constant from the depot (ka, 1/h)")           # Table 2, Absorption rate (ka), TV
    lktr <- log(2.47)  ; label("Transit rate constant of the tablet absorption chain (ktr, 1/h)")         # Table 2, Transit rate for tablets (ktr), TV
    lcl  <- log(204)   ; label("Apparent oral clearance at a 300 mg dose (CL/F, L/h)")                    # Table 2, Oral clearance (CL/F), TV
    lvc  <- log(65.5)  ; label("Apparent central volume of distribution (Vc/F, L)")                       # Table 2, Volume of distribution of the central compartment (Vc/F), TV
    lq   <- log(52.6)  ; label("Apparent intercompartmental clearance at a 300 mg dose (Q/F, L/h)")       # Table 2, Intercompartmental apparent clearance (Q/F), TV
    lvp  <- log(356)   ; label("Apparent peripheral volume of distribution at a 300 mg dose (Vp/F, L)")   # Table 2, Apparent volume of distribution for the peripheral compartment (Vp/F), TV

    # Empirical dose-nonlinearity exponents, all ESTIMATED (each row carries an
    # RSE and a bootstrap 95% CI). Reference dose 300 mg. Negative exponents mean
    # clearance and distribution volumes FALL as the dose rises, i.e. exposure is
    # more than dose proportional -- the paper's 'Observed nonlinearity in PK
    # parameters was addressed by including dose as a covariate on CL/F, Q/F and Vp/F'.
    e_dose_cl <- -0.278 ; label("Exponent on (DOSE/300) for CL/F (unitless)")  # Table 2, Oral clearance (CL/F), Dose effect (exponent)
    e_dose_q  <- -0.554 ; label("Exponent on (DOSE/300) for Q/F (unitless)")   # Table 2, Intercompartmental apparent clearance (Q/F), Dose effect (exponent)
    e_dose_vp <- -0.684 ; label("Exponent on (DOSE/300) for Vp/F (unitless)")  # Table 2, Apparent volume of distribution for the peripheral compartment (Vp/F), Dose effect (exponent)

    # Weight allometry, reference 70 kg. NOT estimated: 'weight was added to the
    # model prior to the pediatric simulations using power functions with weight
    # centered at 70 kg'. All four are fixed structural assumptions of the
    # paediatric simulation, so all four are fixed().
    #   - CL/F 0.75 is printed by the paper: '(Weight/70)^0.75'. Recovered
    #     independently as 0.763 (median over the eleven Table S3 weight/dose rows).
    #   - Vc/F is printed as '(Weight/70)^-1'; the sign is a typesetting error
    #     (a negative exponent on a volume predicts a LARGER central volume in a
    #     smaller subject) and is shipped as +1. Operator ruling, sidecar
    #     request-001 q2 = A.
    #   - Q/F and Vp/F exponents are NOT printed anywhere in the paper; the
    #     standard allometric values are assumed so that the paediatric
    #     distribution kinetics behind Figure 5 can be reproduced. Operator
    #     ruling, sidecar request-001 q2 = A, conditional on reproducing the
    #     paper's own answer keys -- see the vignette Errata.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on (WT/70) for CL/F (unitless)")  # Paper 'PopPK model', CL/F equation
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vc/F (unitless)")  # Paper 'PopPK model', Vc/F equation (printed -1; sign corrected)
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on (WT/70) for Q/F (unitless)")   # Not printed; standard allometric exponent assumed
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vp/F (unitless)")  # Not printed; standard allometric exponent assumed

    # Population PK between-subject variability. Table 2's BSV 'Estimate' column
    # holds omega STANDARD DEVIATIONS while its 95% CI column is on the omega
    # VARIANCE, so the ini() values below are the SQUARES of the printed numbers.
    # Two independent checks: (1) across the 8 BSV rows of Table 2 the SQUARE of
    # the printed estimate lies inside the row's own CI in all 8, while the
    # printed estimate itself lies outside its own CI in 6 of the 8 (the Vp/F
    # and ktr CIs are wide enough to admit either reading); (2) the squares
    # reproduce the prose CVs -- sqrt(exp(0.353^2)-1) = 36.4% matches 'a
    # coefficient of variation (CV) of 36%' and sqrt(exp(0.562^2)-1) = 60.9%
    # matches 'a CV of 61%', whereas reading the printed numbers as variances
    # would give 65% and 87%.
    etalcl  ~ 0.124609  # Table 2, CL/F BSV  0.353 (SD) -> 0.353^2
    etalvc  ~ 0.315844  # Table 2, Vc/F BSV  0.562 (SD) -> 0.562^2
    etalq   ~ 0.174724  # Table 2, Q/F BSV   0.418 (SD) -> 0.418^2
    etalvp  ~ 0.379456  # Table 2, Vp/F BSV  0.616 (SD) -> 0.616^2
    etalktr ~ 0.962361  # Table 2, ktr BSV   0.981 (SD) -> 0.981^2
    etalka  ~ fixed(0)  # Table 2, ka BSV row: zero, not estimated; Appendix S1 popPK $OMEGA 5 agrees

    # Soticlestat plasma residual error. Appendix S1 popPK $ERROR:
    # W = SQRT((THETA(1)/100*IPRED)**2 + THETA(2)**2), so THETA(1) is a PERCENT
    # (45.4% -> 0.454 as a fraction) and the additive term is in ng/mL.
    propSd <- 0.454        ; label("Proportional residual error for soticlestat plasma concentration (fraction)")  # Table 2, Residual variability Proportional (%) 45.4
    addSd  <- fixed(0.001) ; label("Additive residual error for soticlestat plasma concentration (ng/mL)")         # Table 2, Residual variability Additive (ng/mL) 0.001, Fixed

    # ---------------------------------------------------------------------
    # CH24H enzyme-occupancy (PK/EO) model -- Table 2, 'EO (PK/EO) model'
    # panel. A direct link between the effect-site concentration and brain
    # occupancy, with an effect compartment carrying the observed hysteresis.
    # ---------------------------------------------------------------------
    lke0     <- log(0.255)      ; label("Plasma-to-brain effect-site equilibration rate constant (kEO, 1/h)")     # Table 2, Delay rate (kEO), TV
    lemax    <- fixed(log(100)) ; label("Maximum brain CH24H enzyme occupancy (Emax, % occupancy)")               # Table 2, Maximum EO (Emax), TV, Fixed; Appendix S1 PK/EO $THETA 4 '(100) FIX'
    lec50    <- log(5.86)       ; label("Effect-site concentration for 50% maximum CH24H occupancy (EC50, ng/mL)")  # Table 2, Effect-site concentration for 50% maximum effect (EC50), TV
    lhill_eo <- log(0.769)      ; label("Sigmoidicity of the CH24H enzyme-occupancy relationship (unitless)")     # Table 2, Shape parameter (gamma), TV

    etalec50 ~ 0.478864  # Table 2, EO EC50 BSV 0.692 (SD) -> 0.692^2

    # Enzyme occupancy is a PERCENTAGE, so this additive residual is in percentage
    # points of occupancy. Table 2 labels the row '(ng/mL)', which is a carried-over
    # header from the PK panel -- see the vignette Errata. The proportional
    # component was fixed to zero, so the EO residual is purely additive.
    addSd_occ <- 2.88 ; label("Additive residual error for brain CH24H enzyme occupancy (% occupancy)")  # Table 2, EO Residual variability Additive 2.88; Proportional '0 Fixed'

    # ---------------------------------------------------------------------
    # 24HC (PK/PD) model -- Table 2, '24HC (PK/PD) model' panel. A
    # semimechanistic inhibitory indirect-response turnover model in which the
    # effect-site concentration inhibits 24HC production. The PK/EO parameters
    # above were FIXED at their PK/EO estimates for this step ('The PK/PD model
    # was developed by fixing EO model parameters to estimates obtained in the
    # PK/EO modeling step'; Appendix S1 PK/24HC $THETA 3-6 all carry FIX).
    # ---------------------------------------------------------------------
    lrbase <- log(45.9)   ; label("Typical baseline plasma 24HC concentration (BL24HC, ng/mL)")            # Table 2, Baseline 24HC (BL24HC), TV
    lkout  <- log(0.0182) ; label("First-order 24HC degradation rate constant (kout, 1/h)")                # Table 2, 24HC degradation rate (kout), TV
    limax  <- log(78.2)   ; label("Maximum inhibition of 24HC production (Imax, % of the production rate)")  # Table 2, Maximum 24HC inhibition (Imax), TV
    lic50  <- log(5.21)   ; label("Effect-site concentration for 50% of maximum 24HC inhibition (IC50, ng/mL)")  # Table 2, Effect-site concentration for 50% maximum effect (IC50), TV

    # The paper's printed 24HC equation carries a shape parameter on the
    # effect-site concentration, but Appendix S1 PK/24HC $THETA 11 sets it to
    # '1 FIX' and $DES writes the sigmoid without an exponent
    # ('EFF = 1 - IMAX*A(6)/(A(6)+IC50)'), so the exponent is inert. It is kept
    # explicit and fixed at 1 to preserve the published structure.
    lhill_hc24 <- fixed(log(1)) ; label("Sigmoidicity of the 24HC inhibition relationship (unitless)")  # Appendix S1 PK/24HC $THETA 11 IGAM '1 FIX'

    etalrbase ~ 0.224676  # Table 2, BL24HC BSV 0.474 (SD) -> 0.474^2
    etalkout  ~ 0.398161  # Table 2, kout BSV   0.631 (SD) -> 0.631^2

    # Plasma 24HC residual error. Appendix S1 PK/24HC $ERROR:
    # W = SQRT((THETA(1)*IPRED)**2 + THETA(2)**2). Unlike the popPK stream there
    # is no /100 here, so the fixed 0.001 is a FRACTION (0.1%), not 0.001%.
    # Table 2 labels the row 'Proportional (%)' -- see the vignette Errata. The
    # term is negligible either way; the 24HC residual is effectively additive.
    propSd_hc24 <- fixed(0.001) ; label("Proportional residual error for plasma 24HC (fraction)")  # Table 2, 24HC Residual variability Proportional (%) 0.001, Fixed
    addSd_hc24  <- 3.4          ; label("Additive residual error for plasma 24HC (ng/mL)")         # Table 2, 24HC Residual variability Additive (ng/mL) 3.4
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters (Appendix S1 popPK $PK)
    # ---------------------------------------------------------------------
    # The dose-nonlinearity terms act on the TYPICAL VALUE, inside the
    # exponential IIV, exactly as the control stream writes them
    # ('TVCL = THETA(4)*((DOSE/300)**THETA(9))' then 'CL = TVCL*EXP(ETA(1))').
    # Setting WT = 70 collapses every allometric factor to 1 and recovers the
    # estimated adult model.
    ka  <- exp(lka + etalka)
    ktr <- exp(lktr + etalktr)

    cl <- exp(lcl + etalcl) * (DOSE / 300)^e_dose_cl * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq + etalq) * (DOSE / 300)^e_dose_q * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (DOSE / 300)^e_dose_vp * (WT / 70)^e_wt_vp

    # ---------------------------------------------------------------------
    # 2. Effect-site and 24HC turnover parameters (Appendix S1 PK/24HC $PK)
    # ---------------------------------------------------------------------
    ke0     <- exp(lke0)
    emax    <- exp(lemax)
    ec50    <- exp(lec50 + etalec50)
    hill_eo <- exp(lhill_eo)

    rbase     <- exp(lrbase + etalrbase)
    kout      <- exp(lkout + etalkout)
    kin       <- rbase * kout
    imax      <- exp(limax)
    ic50      <- exp(lic50)
    hill_hc24 <- exp(lhill_hc24)

    # ---------------------------------------------------------------------
    # 3. ODE system (Appendix S2, confirmed against Appendix S1 $DES)
    # ---------------------------------------------------------------------
    # Dose in mg over a volume in L gives mg/L; the factor 1000 converts to
    # ng/mL and reproduces the NONMEM scale factor S2 = Vc/1000.
    Cc <- 1000 * central / vc

    # Tablet absorption route: transit1 -> transit2 -> depot, every transfer at
    # ktr ('The transit compartment 1 is the dosing compartment for the tablet
    # formulation', Appendix S2). The oral solution doses straight into depot,
    # which is the NONMEM DEFDOSE compartment. Both routes then share ka.
    # NONMEM slot map: A(4) = transit1, A(5) = transit2, A(1) = depot,
    # A(2) = central, A(3) = peripheral1, A(6) = effect, A(7) = hc24.
    d/dt(transit1)    <- -ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(depot)       <-  ktr * transit2 - ka * depot
    d/dt(central)     <-  ka * depot - cl / vc * central -
                          q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # Effect compartment. The state is a CONCENTRATION in ng/mL, not an amount:
    # NONMEM writes DADT(6) = KEO*(A(2)/S2 - A(6)), i.e. ke0*(Cc - effect).
    d/dt(effect) <- ke0 * (Cc - effect)

    # The occupancy sigmoid raises the effect-site concentration to a FRACTIONAL
    # power (hill_eo = 0.769), which is undefined for a negative argument. The
    # effect site is driven by a non-negative plasma concentration and so cannot
    # go negative mathematically, but a stiff-solver undershoot can return a tiny
    # negative value; max() keeps the fractional power defined.
    ce <- max(effect, 0)

    # Fractional inhibition of 24HC synthesis. Imax is a PERCENT in Table 2
    # (78.2), so it is divided by 100 to give the fraction the control stream's
    # IMAX carries directly.
    drugEff <- imax / 100 * ce^hill_hc24 / (ce^hill_hc24 + ic50^hill_hc24)

    # Indirect-response turnover: the pool starts at its own baseline and the
    # production rate is set so the drug-free system is at steady state
    # (NONMEM 'A_0(7) = BL' and 'KIN = BL*KOUT').
    d/dt(hc24) <- kin * (1 - drugEff) - kout * hc24
    hc24(0)    <- rbase

    # ---------------------------------------------------------------------
    # 4. Read-outs and observations (Appendix S1 $ERROR)
    # ---------------------------------------------------------------------
    # Brain CH24H enzyme occupancy (% occupancy). Emax is fixed at 100, so the
    # sigmoid returns a percentage directly.
    occ <- emax * ce^hill_eo / (ec50^hill_eo + ce^hill_eo)

    # Percent change from the individual's own baseline 24HC, the quantity the
    # paper plots in Figures 4c and 5c.
    hc24Chg <- 100 * (hc24 / rbase - 1)

    Cc   ~ prop(propSd) + add(addSd)
    occ  ~ add(addSd_occ)
    hc24 ~ prop(propSd_hc24) + add(addSd_hc24)
  })
}
