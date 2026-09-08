Qin_2025_ropeginterferon <- function() {
  description <- paste0(
    "Quasi-equilibrium target-mediated drug disposition (TMDD) ",
    "population pharmacokinetic model for subcutaneous ",
    "ropeginterferon alfa-2b (ropeg, a mono-PEGylated interferon ",
    "alfa-2b) in 126 pooled Chinese and Japanese participants: 48 ",
    "healthy volunteers given a single 90-300 ug dose (phase I ",
    "A17-101 and A17-102) and 78 patients with polycythaemia vera ",
    "given 50-500 ug every 2 weeks (phase II A19-201 slow titration ",
    "and A20-202 fast titration). First-order absorption with a lag ",
    "into a single serum compartment carrying linear clearance plus ",
    "saturable binding to a turnover target pool; the drug-target ",
    "complex is internalised at kint. Body mass index acts on ",
    "clearance as a power function. Under chronic dosing in patients ",
    "the target pool declines from its baseline binding capacity to ",
    "a lower steady-state capacity, which healthy single-dose ",
    "participants do not experience. THE PUBLISHED PARAMETER TABLE ",
    "LABELS EVERY RATE PER HOUR BUT THE VALUES ARE PER DAY; see the ",
    "vignette Errata for the three independent proofs. Companion ",
    "pharmacodynamic and exposure-response models in the ",
    "Qin_2025_ropeginterferon_* family."
  )
  reference <- paste(
    "Qin A, Shimoda K, Suo S, Fu R, Kirito K, Wu D, Liao J, Chen H, Wu L,",
    "Su X, Gao Y, Sato T, Li Y, Zhang J, Shen W, Wang W, Zhang L, Jin J,",
    "Komatsu N.",
    "Population pharmacokinetics-pharmacodynamics and exposure-response of",
    "ropeginterferon alfa-2b in Chinese and Japanese patients with",
    "polycythemia vera.",
    "Pharmacol Res Perspect. 2025;13(3):e70109.",
    "doi:10.1002/prp2.70109.",
    sep = " "
  )
  vignette <- "Qin_2025_ropeginterferon"
  units <- list(
    time          = "day (NOT hour: Qin 2025 Table 2 labels the rate constants h^-1, but the values are per day. The Methods text itself reports kdec in day^-1, and three independent anchors confirm the day reading -- see the vignette Errata)",
    dosing        = "ug (micrograms of ropeginterferon alfa-2b, subcutaneous)",
    concentration = "ug/L total serum ropeg (Cc), which is numerically identical to the ng/mL that Qin 2025 reports. Doses in ug with vc in L give ug/L directly, so no conversion constant appears in model()"
  )

  compartmentData <- list(
    depot        = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "administration site", verified = FALSE),
    central      = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "serum",                             verified = FALSE),
    total_target = list(analyte = "ropeg target (total, free plus drug-bound binding capacity)", units = "ug/L (= ng/mL)", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Baseline body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final population PK model",
        "(Qin 2025 Results 3.2.1: 'Body mass index (BMI) was included as",
        "a covariate of clearance'). Enters as the paper's Equation (1)",
        "power model, P_i = P_TV * (COV / COV_med)^theta, with COV_med the",
        "MEDIAN of the pooled 126-participant analysis set: 23.1 kg/m^2",
        "(Qin 2025 Table 1, 'Overall' column, range 17.4-32.2). The",
        "centring value is the pooled median rather than a per-study",
        "median because Equation (1) is written against the population",
        "median and Table 1 reports the Overall column for exactly this",
        "purpose. Qin 2025 Discussion confirms the direction and the",
        "functional form: 'BMI was identified as a significant covariate,",
        "with clearance increasing nonlinearly with BMI'. Note the paper",
        "screened weight, BMI and body surface area separately and",
        "retained BMI, so a downstream user must supply BMI itself rather",
        "than deriving clearance from weight."
      ),
      source_name        = "BMI (body mass index)"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator; 1 = healthy volunteer, 0 = patient with polycythaemia vera.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with polycythaemia vera; the 78 phase II participants of A19-201 and A20-202)",
      notes              = paste(
        "Gates the chronic decline of the target pool, NOT any structural",
        "PK parameter. Qin 2025 Methods 2.4.1 states it directly: 'In",
        "healthy subjects who received only a single dose, as well as in",
        "PV patients prior to the start time of target mediation (TSTART),",
        "Rtot was simply Rtot0.' So the decline from rtot0 to rtot_ss runs",
        "only when DIS_HEALTHY = 0 and only after t_start. Set",
        "DIS_HEALTHY = 1 to reproduce the 48 phase I single-dose healthy",
        "participants (A17-101, n = 18 Chinese; A17-102, n = 30 Japanese",
        "and Caucasian), for whom the binding capacity is held at rtot0",
        "for the whole profile. Because the paper reports no separate",
        "healthy-versus-patient effect on CL, Vc or Ka, this column",
        "changes nothing else in the model."
      ),
      source_name        = "healthy-volunteer study flag (A17-101 / A17-102 versus A19-201 / A20-202); the source control stream name is not published"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate search (Qin 2025 Table 1 lists",
        "it among the assessed demographics) but NOT retained in the final",
        "population PK model, which kept BMI alone on clearance. Pooled",
        "median 62.8 kg, range 43.6-91.0 (Table 1, Overall). Body weight",
        "IS retained in the companion week-24 JAK2 V617F",
        "exposure-response model",
        "(Qin_2025_ropeginterferon_jak2_week24), so the null result here is",
        "specific to the PK layer."
      )
    ),
    BSA = list(
      description = "Baseline body surface area.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but not retained. Pooled median 1.73 m^2, range 1.39-2.09 (Qin 2025 Table 1, Overall)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained. Pooled median 43.5 years, range 21.0-72.0 (Qin 2025 Table 1, Overall); the healthy phase I cohorts are much younger (median 27-30 years) than the PV cohorts (median 54-56 years)."
    ),
    SEXF = list(
      description = "Female-sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Qin 2025 Equation (2) is written with sex as the worked example of a categorical covariate, but no categorical covariate survived into the final model. Pooled 47/126 female (37.3%) (Table 1, Overall)."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened but not retained. Pooled median 111 mL/min, range 43.3-169 (Qin 2025 Table 1, Overall)."
    ),
    CREAT = list(
      description = "Baseline serum creatinine.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained. Pooled median 67.2 umol/L, range 41.0-106 (Qin 2025 Table 1, Overall)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained on any PK parameter. Pooled median 17.1 U/L, range 7.00-51.0, 8 records missing (Qin 2025 Table 1, Overall). ALT increase IS a modelled exposure-safety endpoint in the companion Qin_2025_ropeginterferon_alt_increase model."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained on any PK parameter. Pooled median 20.0 U/L, range 11.4-36.0, 8 records missing (Qin 2025 Table 1, Overall)."
    ),
    HGB = list(
      description = "Baseline haemoglobin concentration.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained. Pooled median 153 g/L, range 116-219 (Qin 2025 Table 1, Overall)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 126L,
    n_studies      = 4L,
    n_observations = "not reported as a record count; PK sampling was pre-dose and 1, 3, 6, 9, 12, 24, 36, 48, 72, 96, 120, 144, 168, 192, 240, 288, 336, 504 and 672 h post-dose in both phase I studies, and weeks 0 and 28 (A19-201) or weeks 0 and 12 (A20-202) pre-dose plus 48, 96 and 168 h post-dose with trough concentrations at every visit in the phase II studies",
    age_range      = "median 43.5 years, range 21.0-72.0 pooled (Qin 2025 Table 1, Overall); healthy phase I median 27.0-30.0 years, PV phase II median 54.0-56.0 years",
    weight_range   = "median 62.8 kg, range 43.6-91.0 (Qin 2025 Table 1, Overall)",
    bmi_range      = "median 23.1 kg/m^2, range 17.4-32.2 (Qin 2025 Table 1, Overall); this median is the centring value of the clearance covariate model",
    sex_female_pct = 37.3,
    race_ethnicity = "Chinese (A17-101 n = 18 healthy; A20-202 n = 49 PV) and Japanese (A19-201 n = 29 PV); A17-102 enrolled 30 healthy Japanese and Caucasian participants. Qin 2025 Methods 2.1: 'Only participants from Japan and China who were administered ropeg were included in the analyses.'",
    disease_state  = "48 healthy volunteers and 78 patients with polycythaemia vera. All A20-202 patients and all but two A19-201 patients carried the JAK2 V617F driver mutation; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline JAK2 V617F allele burden median 77.8% (A19-201) and 61.2% (A20-202)",
    dose_range     = "phase I single subcutaneous doses of 90-270 ug (A17-101) and 100-300 ug (A17-102); phase II subcutaneous doses every 2 weeks, A19-201 starting at 100 ug (or 50 ug on prior cytoreductive therapy) titrated in 50 ug steps to a 500 ug maximum (slow titration) and A20-202 starting at 250 ug with titration to 350 ug at week 2 and 500 ug from week 4 (fast titration)",
    regions        = "China (A17-101, A20-202) and Japan (A17-102, A19-201); A17-102 also enrolled Caucasian participants who were excluded from these analyses",
    trial_registration = "A17-102 NCT03546465; A19-201 NCT04182100; A20-202 NCT05485948; A17-101 CTR20190451 (chinadrugtrials.org.cn)",
    notes          = paste0(
      "Ropeginterferon alfa-2b is approved for polycythaemia vera by ",
      "the EMA (2019), the FDA (2021) and in Japan (2023) and China. ",
      "The pooled analysis set spans a 10-fold single-dose range in ",
      "healthy participants and up to 100 weeks of every-2-week ",
      "dosing in patients, which is what identifies both the linear ",
      "and the target-mediated elimination arms. All fixed and random ",
      "effects were estimated with relative standard errors below ",
      "15% (Qin 2025 Results 3.2.1)."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Table 2, "Parameter estimates of PopPK model". NONMEM
    # 7.5.0, FOCE with eta-epsilon interaction.
    #
    # TIME UNIT -- READ THIS BEFORE CHANGING ANY VALUE.
    # Table 2 prints "CL (L h-1)", "Ka (h-1)", "kint (h-1)",
    # "kdeg (h-1)", "kdec (h-1)" and "TSTART (h)". Those unit labels are
    # WRONG; every one of these values is per DAY. Four pieces of
    # evidence, three of them independent of each other:
    #
    #   (a) The paper contradicts its own table. Methods 2.4.1 states
    #       kdec is "in day-1" while Table 2 labels the same parameter
    #       h-1. Methods 2.4.1 likewise gives kout "in day-1" against
    #       Table 3's h-1.
    #   (b) Tmax. With ka = 0.18 and kel = CL/Vc = 0.753/3.29 = 0.2289
    #       per DAY, Tmax = ln(ka/kel)/(ka-kel) = 4.91 d = 118 h, and
    #       the terminal (flip-flop) half-life is ln(2)/ka = 92 h. The
    #       supplementary PK visual predictive check (Figure S2,
    #       stratified by dose) peaks at 100-150 h after dose in every
    #       single-dose panel. Read per HOUR the same algebra gives
    #       Tmax = 4.9 HOURS, which the VPC excludes outright.
    #   (c) Cmax. A single 300 ug dose gives 29.6 ng/mL under the day
    #       reading; the Figure S2 "DOSEUG == 300" panel peaks at
    #       25-30 ng/mL.
    #   (d) Cavg. Dose/(CL*tau) for 500 ug every 14 days is 47.4 ng/mL
    #       under the day reading, and for the ~210 ug average dose of
    #       the slow-titration arm is 19.9 ng/mL. Figure 3 puts the
    #       observed median Cavg,0-24W at about 37.5 ng/mL (A20-202,
    #       fast titration) and about 20 ng/mL (A19-201, slow
    #       titration). Read per HOUR, 500 ug every 336 h gives
    #       1.98 ng/mL -- a 20-fold miss.
    #
    # The values below are therefore entered UNCHANGED from Table 2 and
    # the model time unit is declared as days in units$time. Nothing is
    # rescaled; only the label is corrected.
    #
    # DOSE AND CONCENTRATION UNITS. Doses are in ug and vc is in L, so
    # central/vc is ug/L, which equals ng/mL exactly. KD and the target
    # capacities are reported in ng/mL and are used as-is.
    # ==================================================================

    # ----- Absorption -----
    lka   <- log(0.18)   ; label("First-order absorption rate constant from the subcutaneous depot (Ka, 1/day)")     # Qin 2025 Table 2: Ka 0.18 (RSE 1.04%), printed as h-1; per day (see the unit note above)
    ltlag <- log(0.62)   ; label("Absorption lag time on the subcutaneous depot (ALAG, day)")                        # Qin 2025 Table 2: absorption lag time 0.62 (RSE 0.421%), printed as h. Carried in DAYS for internal consistency with the rest of the table; a 0.62 h reading changes Cmax and AUC by under 1% because Tmax is ~118 h (see vignette Errata)

    # ----- Linear disposition -----
    lcl   <- log(0.753)  ; label("Linear clearance from the serum compartment at the median BMI of 23.1 kg/m^2 (CL, L/day)")  # Qin 2025 Table 2: CL 0.753 (RSE 0.846%), printed as L h-1; per day (see the unit note above)
    lvc   <- log(3.29)   ; label("Volume of the serum compartment (Vc, L)")                                          # Qin 2025 Table 2: Vc 3.29 (RSE 1.03%)

    # ----- Target-mediated disposition (quasi-equilibrium) -----
    lrtot0    <- log(0.317)  ; label("Baseline maximum target binding capacity, the initial condition of the total-target pool (Rtot0, ng/mL)")  # Qin 2025 Table 2: Rtot0 0.317 (RSE 0.904%)
    lrtot_ss  <- log(0.012)  ; label("Steady-state maximum target binding capacity reached under chronic dosing in patients (Rtot,SS, ng/mL)")   # Qin 2025 Table 2: Rtot,SS 0.012 (RSE 0.997%)
    lkint     <- log(0.0223) ; label("First-order elimination rate constant of the drug-target complex (kint, 1/day)")                           # Qin 2025 Table 2: kint 0.0223 (RSE 0.88%), printed as h-1; per day
    lkdeg     <- log(0.51)   ; label("First-order degradation rate constant of free target (kdeg, 1/day)")                                       # Qin 2025 Table 2: kdeg 0.51 (RSE 0.576%), printed as h-1; per day
    lkd       <- log(0.0662) ; label("Equilibrium dissociation constant of ropeg for its target (KD, ng/mL)")                                    # Qin 2025 Table 2: KD 0.0662 (RSE 0.981%)
    lkdecay   <- log(0.0255) ; label("First-order rate constant of the decline in binding capacity from Rtot0 to Rtot,SS (kdec, 1/day)")         # Qin 2025 Table 2: kdec 0.0255 (RSE 0.892%). Methods 2.4.1 states this one explicitly as "in day-1", which is the direct textual confirmation that Table 2's h-1 labels are wrong
    t_start   <- fixed(7)    ; label("Time after the first dose at which target mediation begins in patients (TSTART, day)")                     # Qin 2025 Table 2: TSTART 7 (FIX), printed as h; carried in DAYS with the rest of the table. Immaterial either way against the 27-day half-life of the kdec decline

    # ----- Covariate effect -----
    e_bmi_cl <- 0.813 ; label("Power exponent on (BMI / 23.1) for clearance (unitless)")  # Qin 2025 Table 2: "BMI effect on CL" 0.813 (RSE 9.39%). Applied through the paper's Equation (1) power form P_i = P_TV*(COV/COV_med)^theta with COV_med = 23.1 kg/m^2 (Table 1, Overall median)

    # ----- Inter-individual variability -----
    # Table 2 reports IIV as a coefficient of variation in per cent for
    # a log-normally distributed random effect (Methods 2.4.2: "IIV was
    # distributed log-normally"), so the internal variance is
    # omega^2 = log(CV^2 + 1). No correlations between PK random
    # effects are reported.
    etalka    ~ 0.397837  ; label("IIV variance on log Ka")     # Qin 2025 Table 2: 69.9% CV (RSE 8.2%, shrinkage 14.3%); log(0.699^2 + 1) = 0.397837
    etalcl    ~ 0.117435  ; label("IIV variance on log CL")     # Qin 2025 Table 2: 35.3% CV (RSE 7.64%, shrinkage 9.87%); log(0.353^2 + 1) = 0.117435
    etalvc    ~ 0.639175  ; label("IIV variance on log Vc")     # Qin 2025 Table 2: 94.6% CV (RSE 8.72%, shrinkage 16.6%); log(0.946^2 + 1) = 0.639175
    etalrtot0 ~ 2.128041  ; label("IIV variance on log Rtot0")  # Qin 2025 Table 2: 272% CV (RSE 6.63%, shrinkage 20.6%); log(2.72^2 + 1) = 2.128041
    etalkint  ~ 0.940983  ; label("IIV variance on log kint")   # Qin 2025 Table 2: 125% CV (RSE 10.3%, shrinkage 17.2%); log(1.25^2 + 1) = 0.940983

    # ----- Residual unexplained variability -----
    # Methods 2.4.2: "RUV was normally distributed and described as an
    # additive, proportional, or mixed RUV form"; Results 3.2.1 confirms
    # the final model kept BOTH terms.
    propSd <- 0.197 ; label("Proportional residual SD on total serum ropeg concentration (fraction)")  # Qin 2025 Table 2: proportional residual error 19.7% (RSE 1.01%, shrinkage 9.26%)
    addSd  <- 0.542 ; label("Additive residual SD on total serum ropeg concentration (ng/mL)")         # Qin 2025 Table 2: "additional residual error" 0.542 ng/mL (RSE 2.08%, shrinkage 9.26%)
  })

  model({
    # ----- Individual parameters -----
    ka     <- exp(lka + etalka)
    tlag   <- exp(ltlag)
    # Qin 2025 Equation (1), continuous-covariate power model:
    #   P_i = P_TV * (COV / COV_med)^theta,  COV_med = 23.1 kg/m^2.
    cl     <- exp(lcl + etalcl) * (BMI / 23.1)^e_bmi_cl
    vc     <- exp(lvc + etalvc)
    rtot0  <- exp(lrtot0 + etalrtot0)
    rtot_ss <- exp(lrtot_ss)
    kint   <- exp(lkint + etalkint)
    kdeg   <- exp(lkdeg)
    kd     <- exp(lkd)
    kdecay <- exp(lkdecay)

    # ----- Time-varying target binding capacity -----
    # Qin 2025 Methods 2.4.1: "The zero-order rate constant for target
    # synthesis (ksyn) was the product of Rtot and kdeg. Differences in
    # target turnover over time were described by allowing Rtot to
    # decrease from an initial starting value (Rtot0) to a new
    # steady-state value (Rtot,SS) with an estimated first-order rate
    # constant (kdec). In healthy subjects who received only a single
    # dose, as well as in PV patients prior to the start time of target
    # mediation (TSTART), Rtot was simply Rtot0."
    #
    # tdec is clamped at zero rather than written as an if/else so that
    # the expression is smooth for the ODE solver and never evaluates
    # exp() of a large positive argument before being multiplied by a
    # zero indicator.
    tdec      <- max(t - t_start, 0)
    rtot_cap  <- rtot0 + (1 - DIS_HEALTHY) * (rtot_ss - rtot0) * (1 - exp(-kdecay * tdec))
    ksyn      <- kdeg * rtot_cap

    # ----- Quasi-equilibrium binding -----
    # With total drug ctot = central/vc and total target total_target,
    # free drug solves ctot = cfree + total_target*cfree/(KD + cfree),
    # i.e. the positive root of the binding quadratic. This is the
    # Figure 1A structure: drug and target associate reversibly with
    # affinity KD, free target degrades at kdeg and the complex is
    # internalised at kint.
    ctot    <- central / vc
    disc    <- ctot - total_target - kd
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    complex <- total_target * cfree / (kd + cfree)

    total_target(0) <- rtot0

    # ----- ODEs -----
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - cl * cfree - kint * complex * vc
    d/dt(total_target) <-  ksyn - kdeg * (total_target - complex) - kint * complex
    alag(depot)        <- tlag

    # ----- Observation -----
    # The bioassay measures TOTAL serum ropeg, and Methods 2.4.1 names
    # the pharmacodynamic driver CS, "total serum ropeg concentrations".
    Cc <- ctot
    Cc ~ add(addSd) + prop(propSd)
  })
}
