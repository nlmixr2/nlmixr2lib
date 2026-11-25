SheikhGhadzi_2017_T2DM_IGI_OGTT <- function() {
  description <- "Methodology reference. Simulation (data-generating) model of Sheikh Ghadzi 2017, a study of the power to characterise an antihyperglycemic drug from oral glucose tolerance tests (OGTT) with and without insulin measurements. Glucose and insulin follow the integrated glucose-insulin (IGI) OGTT model of Jauslin 2007 for type 2 diabetes: two-compartment glucose disposition with insulin-independent (CLG) and effect-delayed insulin-dependent (CLGI) clearance, endogenous glucose production set by the fasting steady state, oral glucose absorbed through a depot and a transit compartment with bioavailability BIOG, and one-compartment insulin with secretion driven by an effect-delayed glucose power function (IPRG) and an Emax incretin effect of the transit glucose amount. A HYPOTHETICAL study compound (one-compartment oral PK, author-chosen constants) acts through an Emax function of its plasma concentration on one or more of seven mechanisms of action: incretin effect (moa_incr), basal insulin secretion (moa_bins), CLG (moa_clg), CLGI (moa_clgi), endogenous glucose production (moa_egp, inhibition), glucose absorption (moa_gabs, inhibition) and glucose sensitivity of insulin secretion (moa_gsen). Each Emax was titrated by the authors to a 10 percent reduction of the glucose AUC; the glucose-absorption effect is encoded as a first-order loss of transit glucose (rate Emax x Cc/(EC50 + Cc), 1/min), which reproduces that titration where the parenthesisation printed in the authors' code does not. The 0/1 switches select the mechanisms; the default is the incretin effect alone, the example of the authors' published NONMEM code. Placebo is no drug dose."
  reference <- paste(
    "Sheikh Ghadzi SM, Karlsson MO, Kjellsson MC. Implications for drug",
    "characterization in glucose tolerance tests without insulin: simulation",
    "study of power and predictions using model-based analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(10):686.",
    "doi:10.1002/psp4.12214. PMID 28575547; PMCID PMC5658280.",
    "Values from Table 1 (simulation values); equations from the NONMEM model",
    "file in Supplement 3. The IGI-OGTT model and its glucose/insulin values",
    "are from Jauslin PM, Silber HE, Frey N, et al. An integrated",
    "glucose-insulin model to describe oral glucose tolerance test data in",
    "type 2 diabetics. J Clin Pharmacol. 2007;47(10):1244-1255.",
    "doi:10.1177/0091270007302168.",
    sep = " "
  )
  vignette <- "SheikhGhadzi_2017_T2DM_IGI_OGTT"
  # Glucose amounts are in the NONMEM code's (mg/dL) x L units, i.e. 10 mg:
  # the glucose dose is in mg and f(glucose_depot) carries the 1/10 'L -> dL'
  # correction of Supplement 3 (F6 = BIOG/10).
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "mg/L (drug, Cc); glucose Gc in mg/dL; insulin Ic in mU/L"
  )

  paper_specific_compartments <- c(
    "glucose_depot",
    "glucose_transit",
    "glucose_peripheral",
    "glucose_effect",
    "insulin_effect"
  )

  compartmentData <- list(
    depot = list(
      analyte = "hypothetical antihyperglycemic study compound",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "hypothetical antihyperglycemic study compound",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    glucose_depot = list(
      analyte = "glucose",
      units = "(mg/dL) x L (10 mg)",
      specimen = "administration site",
      verified = TRUE
    ),
    glucose_transit = list(
      analyte = "glucose",
      units = "(mg/dL) x L (10 mg)",
      specimen = "administration site",
      verified = TRUE
    ),
    glucose = list(
      analyte = "glucose",
      units = "(mg/dL) x L (10 mg)",
      specimen = "plasma",
      verified = TRUE
    ),
    glucose_peripheral = list(
      analyte = "glucose",
      units = "(mg/dL) x L (10 mg)",
      specimen = "tissue",
      verified = TRUE
    ),
    insulin = list(
      analyte = "insulin",
      units = "mU",
      specimen = "plasma",
      verified = TRUE
    ),
    glucose_effect = list(
      analyte = "glucose effect on insulin secretion",
      units = "mg/dL",
      specimen = "not applicable",
      verified = TRUE
    ),
    insulin_effect = list(
      analyte = "insulin effect on glucose clearance",
      units = "mU/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 500L,
    n_studies = 0L,
    disease_state = "Virtual patients with type 2 diabetes mellitus undergoing a 75 g oral glucose tolerance test; the glucose-insulin parameters are the Jauslin 2007 IGI-OGTT estimates in type 2 diabetes.",
    dose_range = "Crossover of placebo and 50 mg oral study compound at time 0, followed by 75 g oral glucose at 30 min; glucose and insulin sampled at 0, 30, 60, 90, 120, 150, 180, 210 and 240 min.",
    notes = "Sheikh Ghadzi 2017 Methods 'Study design' and 'Simulation of data'. No real patients: relative power was assessed on 500 simulated subjects and accuracy/precision on 500 simulated studies of 15 subjects. The drug PK constants and the drug Emax/EC50 values are author-chosen (Cmax about 0.243 mg/L at 65 min, 4-h concentration 0.127 mg/L after 50 mg). No covariates."
  )

  ini({
    # -----------------------------------------------------------------------
    # Glucose (Table 1). Footnote a ('Fixed for analysis') -> fixed(); CLG,
    # CLGI and GSS were re-estimated in every analysis and are not fixed().
    # -----------------------------------------------------------------------
    lgss <- log(169); label("Baseline (fasting) glucose concentration GSS (mg/dL)") # Table 1 'Gss, mg/dL' 169
    lvg <- fixed(log(9.33)); label("Central glucose volume of distribution VG (L)") # Table 1 'VG, L' 9.33, footnote a
    lvpg <- fixed(log(8.56)); label("Peripheral glucose volume of distribution VP (L)") # Table 1 'Vp, L' 8.56, footnote a
    lqg <- fixed(log(0.442)); label("Glucose intercompartmental clearance Q (L/min)") # Table 1 'Q, L/min' 0.442, footnote a
    lclg <- log(0.0287); label("Insulin-independent glucose clearance CLG (L/min)") # Table 1 'CLG, L/min' 0.0287
    lclgi <- log(0.0059); label("Insulin-dependent glucose clearance CLGI (L/min per mU/L)") # Table 1 'CLGI, L/min' 0.0059
    lkge <- fixed(log(0.0289)); label("Effect-delay rate constant, glucose on insulin secretion KGE (1/min)") # Table 1 'KGE, min-1' 0.0289, footnote a
    liprg <- fixed(log(1.42)); label("Power of the glucose effect on insulin secretion IPRG (unitless)") # Table 1 'IPRG' 1.42, footnote a
    lbiog <- fixed(log(0.811)); label("Oral glucose bioavailability BIOG (fraction)") # Table 1 'BIOG' 0.811, footnote a
    lkabs <- fixed(log(0.0364)); label("Glucose absorption rate constant KABS (1/min)") # Table 1 'KABS, min-1' 0.0364, footnote a
    lincrmax <- fixed(log(1.47)); label("Maximal incretin effect INCRmax (unitless)") # Table 1 'INCRmax' 1.47, footnote a
    lincr50 <- fixed(log(14.8)); label("Transit glucose amount at half-maximal incretin effect INCR50 ((mg/dL) x L)") # Table 1 'INCR50' 14.8, footnote a

    # -----------------------------------------------------------------------
    # Insulin (Table 1). CLI and ISS carry footnote c (fixed only for the
    # analysis without insulin) and are not fixed() here.
    # -----------------------------------------------------------------------
    liss <- log(8.71); label("Baseline (fasting) insulin concentration ISS (mU/L)") # Table 1 'Iss, mU/L' 8.71, footnote c
    lvi <- fixed(log(6.09)); label("Insulin volume of distribution VI (L)") # Table 1 'VI, L' 6.09, footnote a
    lcli <- log(1.22); label("Insulin clearance CLI (L/min)") # Table 1 'CLI, L/min' 1.22, footnote c
    lkie <- fixed(log(0.0213)); label("Effect-delay rate constant, insulin on CLGI KIE (1/min)") # Table 1 'KIE, min-1' 0.0213, footnote a

    # -----------------------------------------------------------------------
    # Hypothetical study compound: one-compartment oral PK, author-chosen and
    # held fixed in the analysis (Table 1 footnote a; Methods: 'PK was kept
    # fixed').
    # -----------------------------------------------------------------------
    lka <- fixed(log(0.0365)); label("Drug absorption rate constant KaD (1/min)") # Table 1 'KaD, min-1' 0.0365, footnote a
    lcl <- fixed(log(0.682)); label("Drug clearance CLD (L/min)") # Table 1 'CLD, L/min' 0.682, footnote a
    lvc <- fixed(log(155)); label("Drug volume of distribution VD (L)") # Table 1 'VD, L' 155, footnote a

    # -----------------------------------------------------------------------
    # Drug effect: Emax functions of Cc sharing one EC50 (Supplement 3 $DES).
    # Each Emax was titrated to a 10 percent glucose-AUC reduction (Methods)
    # and re-estimated in the analysis, so none is fixed().
    # -----------------------------------------------------------------------
    lec50 <- log(0.1); label("Drug concentration at half-maximal effect EC50D (mg/L)") # Table 1 'EC50D' 0.1; Methods 'EC50 set to 0.1 mg/L'
    emax_incr <- 1.0; label("Maximal drug stimulation of the incretin effect (fraction)") # Table 1 'EmaxD(Incr)' 1.0
    emax_bins <- 0.5; label("Maximal drug stimulation of basal insulin secretion (fraction)") # Table 1 'EmaxD(BINS)' 0.5
    emax_clg <- 2.0; label("Maximal drug stimulation of CLG (fraction)") # Table 1 'EmaxD(CLG)' 2.0
    emax_clgi <- 0.5; label("Maximal drug stimulation of CLGI (fraction)") # Table 1 'EmaxD(CLGI)' 0.5
    emax_egp <- 0.7; label("Maximal drug inhibition of endogenous glucose production (fraction)") # Table 1 'EmaxD(EGP)' 0.7
    emax_gabs <- 0.02; label("Maximal drug-induced first-order loss rate of transit glucose (1/min)") # Table 1 'EmaxD(GABS)' 0.02
    emax_gsen <- 1.5; label("Maximal drug stimulation of the IPRG glucose sensitivity exponent (fraction)") # Table 1 'EmaxD(GSEN)' 1.5

    # Mechanism-of-action switches. Supplement 3 selects the simulated
    # mechanism by commenting the Emax term out of each drug-effect line;
    # 1 keeps the line, 0 comments it out. Default = the published example
    # ('7 MOA-EXAMPLE INCRETIN DRUG EFFECT'). Part 1/2 scenarios switch on
    # one mechanism, Part 3 two (Supplementary Table S1).
    moa_incr <- fixed(1); label("Switch: drug acts on the incretin effect (0/1)") # Supplement 3 'DINCR = 1 + EMAXD1*CONC/(EC50D+CONC)'
    moa_bins <- fixed(0); label("Switch: drug acts on basal insulin secretion (0/1)") # Supplement 3 'DBINS = 1; + EMAXD2*...'
    moa_clg <- fixed(0); label("Switch: drug acts on CLG (0/1)") # Supplement 3 'DCLG = 1; + EMAXD3*...'
    moa_clgi <- fixed(0); label("Switch: drug acts on CLGI (0/1)") # Supplement 3 'DCLGI = 1; + EMAXD4*...'
    moa_egp <- fixed(0); label("Switch: drug acts on endogenous glucose production (0/1)") # Supplement 3 'DEGP = 1; - EMAXD5*...'
    moa_gabs <- fixed(0); label("Switch: drug acts on glucose absorption (0/1)") # Supplement 3 'DGABS = 0; EMAXD6*...'
    moa_gsen <- fixed(0); label("Switch: drug acts on glucose sensitivity (0/1)") # Supplement 3 'DGSEN = 1; + EMAXD7*...'

    # -----------------------------------------------------------------------
    # Between-subject variability. Table 1 reports IIV in percent with
    # omega = (IIV/100)^2: Supplement 3 prints OM1 GSS 0.0196 for 14 and
    # OM16 KAD 0.462 for 68. Footnote b ('Fixed to zero for analysis') and
    # footnote a on the drug PK -> fixed(); CLG, CLGI, GSS, ISS, CLI and
    # EC50D IIVs were estimated in the analysis. BIOG and the Emax values
    # have no IIV. Table 1 'Parameter correlations' are correlation
    # coefficients (a covariance of 0.7225 between VG and VI would exceed
    # sqrt(0.09 x 0.1681)); covariance = r x SD1 x SD2.
    # -----------------------------------------------------------------------
    etalgss ~ 0.0196 # Table 1 Gss IIV 14; Supplement 3 '$OMEGA 0.0196 ;OM1 GSS'
    etalclg ~ 0.3481 # Table 1 CLG IIV 59 -> 0.59^2
    etalclgi ~ 0.2116 # Table 1 CLGI IIV 46 -> 0.46^2
    etalvpg ~ fixed(0.09) # Table 1 Vp IIV 30, footnote b -> 0.30^2
    etalkge ~ fixed(0.0064) # Table 1 KGE IIV 8, footnote b -> 0.08^2
    etaliprg ~ fixed(0.1225) # Table 1 IPRG IIV 35, footnote b -> 0.35^2
    etalkabs ~ fixed(0.0361) # Table 1 KABS IIV 19, footnote b -> 0.19^2
    etalincrmax ~ fixed(0.3025) # Table 1 INCRmax IIV 55, footnote b -> 0.55^2
    etalincr50 ~ fixed(1.2996) # Table 1 INCR50 IIV 114, footnote b -> 1.14^2
    # VG (30), Q (85), VI (41), footnote b; r(VG,Q) = -0.19, r(VG,VI) = 0.7225,
    # r(Q,VI) = -0.122 -> cov -0.04845, 0.08886750, -0.042517.
    etalvg + etalqg + etalvi ~ fixed(c(
      0.09,
      -0.04845, 0.7225,
      0.0888675, -0.042517, 0.1681
    )) # Table 1 VG/Q/VI IIV and 'Parameter correlations'
    etaliss ~ 0.2401 # Table 1 Iss IIV 49, footnote c -> 0.49^2
    etalcli ~ 0.0841 # Table 1 CLI IIV 29, footnote c -> 0.29^2
    etalkie ~ fixed(0.3364) # Table 1 KIE IIV 58, footnote b -> 0.58^2
    etalka ~ fixed(0.4624) # Table 1 KaD IIV 68, footnote a; Supplement 3 '$OMEGA 0.462 ;OM16 KAD'
    # CLD (15) and VD (29.6), footnote a; r(CLD,VD) = 0.0311 -> cov 0.00138084.
    etalcl + etalvc ~ fixed(c(
      0.0225,
      0.00138084, 0.087616
    )) # Table 1 CLD/VD IIV and 'CorrCLD-VD'
    etalec50 ~ 0.09 # Table 1 EC50D IIV 30 -> 0.30^2

    # -----------------------------------------------------------------------
    # Residual error: additive on log-transformed glucose and insulin
    # (Supplement 3 $ERROR: Y = LOG(IPRED) + EPS(1)*W with $SIGMA 1 FIX and
    # W = RESG or RESI), i.e. log-normal with SD RESG / RESI.
    # -----------------------------------------------------------------------
    expSd_Gc <- 0.0732; label("Log-scale residual SD on glucose (unitless)") # Table 1 'RESG' 0.0732
    expSd_Ic <- 0.252; label("Log-scale residual SD on insulin (unitless)") # Table 1 'RESI' 0.252
  })

  model({
    # Individual parameters (Supplement 3 $PK: THETA * EXP(ETA)).
    gss <- exp(lgss + etalgss)
    vg <- exp(lvg + etalvg)
    vpg <- exp(lvpg + etalvpg)
    qg <- exp(lqg + etalqg)
    clg <- exp(lclg + etalclg)
    clgi <- exp(lclgi + etalclgi)
    kge <- exp(lkge + etalkge)
    iprg <- exp(liprg + etaliprg)
    biog <- exp(lbiog)
    kabs <- exp(lkabs + etalkabs)
    incrmax <- exp(lincrmax + etalincrmax)
    incr50 <- exp(lincr50 + etalincr50)
    iss <- exp(liss + etaliss)
    vi <- exp(lvi + etalvi)
    cli <- exp(lcli + etalcli)
    kie <- exp(lkie + etalkie)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ec50 <- exp(lec50 + etalec50)

    # Micro-constants and fasting steady state (Supplement 3 $PK).
    k15 <- qg / vg
    k51 <- qg / vpg
    kg <- clg / vg
    kgi <- clgi / vg
    ki <- cli / vi
    egp0 <- gss * (clg + clgi * iss)
    eis0 <- cli * iss

    # Drug PK (Supplement 3 $DES, DADT(8) and DADT(9)).
    Cc <- central / vc

    # Drug effects (Supplement 3 $DES); the switches stand for the
    # commented-out Emax terms of the published code.
    edrug <- Cc / (ec50 + Cc)
    dincr <- 1 + moa_incr * emax_incr * edrug
    dbins <- 1 + moa_bins * emax_bins * edrug
    dclg <- 1 + moa_clg * emax_clg * edrug
    dclgi <- 1 + moa_clgi * emax_clgi * edrug
    degp <- 1 - moa_egp * emax_egp * edrug
    dgabs <- moa_gabs * emax_gabs * edrug
    dgsen <- 1 + moa_gsen * emax_gsen * edrug

    # Glucose fluxes.
    egp <- egp0 * degp
    aglu <- kabs * glucose_transit
    dglu <- k15 * glucose - k51 * glucose_peripheral
    eglu <- (kg * dclg + kgi * dclgi * insulin_effect) * glucose

    # Insulin secretion: glucose power function and incretin effect of the
    # transit glucose amount.
    relgss <- glucose_effect / gss + 1e-5
    isecg <- relgss^(iprg * dgsen)
    incr <- 1 + glucose_transit * incrmax * dincr / (glucose_transit + incr50)
    eis <- eis0 * dbins * isecg * incr
    eins <- ki * insulin

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    d/dt(glucose_depot) <- -kabs * glucose_depot
    # Glucose-absorption inhibition is a first-order loss of transit glucose
    # (rate dgabs, 1/min) alongside absorption. Supplement 3 prints
    # KABS*(A(6) - (1+DGABS)*A(7)), which with EmaxD(GABS) = 0.02 changes
    # the glucose AUC by under 1 percent, against the 10 percent titration
    # of Methods and the visible effect in Supplementary Figure S1.
    d/dt(glucose_transit) <- kabs * glucose_depot - (kabs + dgabs) * glucose_transit
    d/dt(glucose) <- egp - eglu - dglu + aglu
    d/dt(glucose_peripheral) <- dglu
    d/dt(insulin) <- eis - eins
    d/dt(glucose_effect) <- kge * (glucose / vg - glucose_effect)
    d/dt(insulin_effect) <- kie * (insulin / vi - insulin_effect)

    # Oral glucose: bioavailability with the 'L -> dL' correction (F6 =
    # BIOG/10) so a dose in mg lands in (mg/dL) x L units.
    f(glucose_depot) <- biog / 10

    # Fasting steady state (Supplement 3 A_0).
    glucose(0) <- gss * vg
    insulin(0) <- iss * vi
    glucose_effect(0) <- gss
    glucose_peripheral(0) <- k15 * gss * vg / k51
    insulin_effect(0) <- iss

    # Observations (Supplement 3 $ERROR, IPRED + 1E-5).
    Gc <- glucose / vg + 1e-5
    Ic <- insulin / vi + 1e-5
    Gc ~ lnorm(expSd_Gc)
    Ic ~ lnorm(expSd_Ic)
  })
}
