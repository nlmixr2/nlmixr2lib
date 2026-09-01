Chanu_2010_methoxyPolyethyleneGlycolEpoetinBeta <- function() {
  description <- paste(
    "One-compartment population PK with first-order subcutaneous absorption",
    "for methoxy polyethylene glycol-epoetin beta (C.E.R.A.), coupled to a",
    "red-blood-cell life-span hemoglobin PD model (Chanu 2010). Both IV and SC",
    "routes are supported. Body weight modifies clearance and volume and age",
    "modifies volume through normalized power models. The PD layer is the",
    "Krzyzanski life-span indirect-response model implemented as a delay",
    "differential equation: hemoglobin production is stimulated by an",
    "Emax function of serum C.E.R.A. concentration, and the same production",
    "term delayed by one red-blood-cell life span acts as the loss term. An",
    "additional SESA term, active only during the first life span after the",
    "switch, returns previously ESA-treated patients from their switch",
    "hemoglobin to their pre-ESA baseline. C-reactive protein and the prior",
    "weekly epoetin dose modify SC50 through normalized power models.",
    "Proportional-plus-additive residual error on concentration; additive",
    "residual error on hemoglobin."
  )
  reference <- paste(
    "Chanu P, Gieschke R, Charoin JE, Pannier A, Reigner B. Population",
    "pharmacokinetic/pharmacodynamic model for C.E.R.A. in both ESA-naive and",
    "ESA-treated chronic kidney disease patients with renal anemia.",
    "Journal of Clinical Pharmacology. 2010;50(5):507-520.",
    "doi:10.1177/0091270009343931"
  )
  vignette <- "Chanu_2010_methoxyPolyethyleneGlycolEpoetinBeta"
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL", hemoglobin = "g/dL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE for all three states: the analyte and
  # specimen were checked against Chanu 2010 (serum C.E.R.A. by ELISA,
  # Methods "PK and PD Assessments"; blood hemoglobin, Methods "Assumptions").
  compartmentData <- list(
    depot   = list(analyte = "methoxy polyethylene glycol-epoetin beta", units = "ug",   specimen = "administration site", verified = TRUE),
    central = list(analyte = "methoxy polyethylene glycol-epoetin beta", units = "ug",   specimen = "serum",               verified = TRUE),
    # hb is a CONCENTRATION state, not an amount state: the paper's PD model is
    # written directly on Hb in g/dL (Chanu 2010 Results eq for Hb'(t)), so
    # d/dt(hb) carries units (g/dL)/day. Recorded as g/dL rather than forced
    # into an amount unit that the source model never uses.
    hb      = list(analyte = "hemoglobin",                              units = "g/dL", specimen = "whole blood",         verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on clearance (exponent 0.571) and on volume",
        "of distribution (exponent 0.443); Chanu 2010 Table III 'Effect of BW",
        "on CL' / 'Effect of BW on V' with the log-linear parameterization of",
        "Methods 'Covariate Model (Stage 2)': TVP = theta_P * (BW/rBW)^theta_BW.",
        "The paper defines rBW as 'the median value of the population studied'",
        "but does not print the pooled median; the reference used here is 70 kg,",
        "which is the exact median in MAXIMA and PROTOS (265 of 400 patients)",
        "and rounds the size-weighted mean of the three study medians",
        "(135*66 + 122*70 + 143*70)/400 = 68.7 kg from Table I. Using 68.7 kg",
        "instead changes typical CL by 1.1% and typical V by 0.8%. See vignette",
        "Assumptions and deviations."
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on volume of distribution (exponent 0.267);",
        "Chanu 2010 Table III 'Effect of age on V'. Reference age 60 years is",
        "not printed in the paper; it is the size-weighted mean of the three",
        "study medians in Table I, (135*54 + 122*57 + 143*68)/400 = 59.9 years,",
        "rounded. The paper's own sensitivity statement -- 'a change in age (by",
        "a factor of 1.5) increased V by 11%' -- is reproduced exactly by the",
        "exponent alone (1.5^0.267 = 1.114) and is independent of the reference",
        "value. See vignette Assumptions and deviations."
      ),
      source_name        = "AGE"
    ),
    CRP = list(
      description        = "C-reactive protein (time-varying)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on SC50 (exponent 0.319); Chanu 2010 Table IV",
        "'Effect of CRP on SC50'. Chanu 2010 Table I lists CRP explicitly under",
        "'Continuous covariates (time varying)', so the column may be supplied",
        "per observation record. Reference 5 mg/L is not printed; it is the",
        "median in AMICUS and MAXIMA (257 of 400 patients) and rounds the",
        "size-weighted mean of the three study medians,",
        "(135*5 + 122*5 + 143*6.5)/400 = 5.5 mg/L (Table I). Using 5.5 mg/L",
        "instead raises typical SC50 by 3.1%. The paper's own sensitivity",
        "statement -- 'a 5-fold increase in CRP predicts a 67% increase in",
        "SC50' -- is reproduced by the exponent alone (5^0.319 = 1.671) and is",
        "independent of the reference value."
      ),
      source_name        = "CRP"
    ),
    ESAD = list(
      description        = "Previous weekly epoetin dose at the switch to C.E.R.A.",
      units              = "IU/week",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chanu 2010 calls this covariate DEPO, 'previous weekly EPO dose at the",
        "start of C.E.R.A. maintenance treatment'; normalized power effect on",
        "SC50 (exponent 0.303, Table IV 'Effect of DEPO on SC50'). Reference",
        "7000 IU/week is not printed; it rounds the size-weighted mean of the",
        "two maintenance-study medians, (122*9000 + 143*5000)/265 = 6842 IU/week",
        "(Table I; AMICUS is ESA-naive and reports NA). ESAD = 0 marks an",
        "ESA-naive patient: the covariate factor is then gated to 1 by the",
        "esadf indicator, following the same idiom as Naik_2013_peginesatide.R.",
        "The paper's own sensitivity statement -- 'a 5-fold increase in DEPO",
        "predicts a 63% increase in SC50' -- is reproduced by the exponent",
        "alone (5^0.303 = 1.629). ESAD ALSO acts as the ESA-naive / ESA-treated",
        "switch for the SESA term and for the hemoglobin initial condition:",
        "ESAD = 0 makes the model start at the individual's own pre-ESA",
        "baseline hb0 with SESA = 0 (the AMICUS correction setting), while",
        "ESAD > 0 makes it start at HGB_BL with SESA active (the MAXIMA /",
        "PROTOS maintenance setting)."
      ),
      source_name        = "DEPO"
    ),
    HGB_BL = list(
      description        = "Hemoglobin concentration at the switch from the previous ESA to C.E.R.A.",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chanu 2010 calls this Hbsw, 'the Hb concentration at the switch from",
        "the previous ESA treatment to C.E.R.A.' (Results, Structural PD Model",
        "narrative). Used as the initial condition of the hb state and inside",
        "the SESA term SESA = (Hb0 - Hbsw)/LS. Only read for previously",
        "ESA-treated patients (ESAD > 0); for ESA-naive patients the model uses",
        "the individual's own drawn pre-ESA baseline hb0 instead, so the column",
        "may be left at any value (conventionally 0 or the observed baseline)",
        "for those subjects. Distinct from the PD PARAMETER Hb0 (canonical",
        "lrbase), which is the pre-ESA baseline and is unobservable in the",
        "maintenance studies."
      ),
      source_name        = "Hbsw"
    ),
    OCC = list(
      description        = "Dosing / PK-sampling occasion index",
      units              = "(index)",
      type               = "categorical",
      reference_category = "1",
      notes              = paste(
        "Values 1, 2, 3 identify the interoccasion-variability occasion.",
        "Chanu 2010 Methods 'Interoccasion Variability Model' defines an",
        "occasion as 'each drug intake that was followed by serum C.E.R.A.",
        "evaluation'; Table II shows that PK sampling followed the doses given",
        "on study days 1, 57 and 127 (AMICUS) or 1, 57 and 141 (MAXIMA,",
        "PROTOS), i.e. three occasions per subject. Decomposed inside model()",
        "into binary indicators oc1..oc3 that multiplex the three IOV etas on",
        "log-CL, following the Jonsson_2011_ethambutol.R idiom (rxode2 parses",
        "but cannot simulate the eta ~ var | occ form). Records outside those",
        "three occasions should carry OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 400L,
    n_studies      = 3L,
    age_range      = "18-89 years (AMICUS 18-89, MAXIMA 23-87, PROTOS 21-85; Chanu 2010 Table I)",
    age_median     = "54 years (AMICUS), 57 years (MAXIMA), 68 years (PROTOS); pooled median not printed",
    weight_range   = "36-131 kg (AMICUS 43-106, MAXIMA 36-131, PROTOS 43-115; Chanu 2010 Table I)",
    weight_median  = "66 kg (AMICUS), 70 kg (MAXIMA), 70 kg (PROTOS); pooled median not printed",
    sex_female_pct = 40.5,
    race_ethnicity = c(White = 76.5, Black = 17.8, Asian = 4.5, Other = 1.3),
    disease_state  = paste(
      "Chronic kidney disease with renal anemia, on dialysis (hemodialysis 386",
      "of 400, peritoneal dialysis 14 of 400). AMICUS (n = 135) enrolled",
      "ESA-naive patients and studied correction of anemia; MAXIMA (n = 122)",
      "and PROTOS (n = 143) enrolled patients already on another",
      "erythropoiesis-stimulating agent and studied maintenance of hemoglobin.",
      "Dose was titrated in every study: to Hb >= 11 g/dL with a rise of >= 1",
      "g/dL from baseline in AMICUS, and to hold Hb within +/- 1 g/dL of",
      "baseline and between 10 and 13 g/dL in MAXIMA and PROTOS."
    ),
    dose_range     = paste(
      "C.E.R.A. intravenously or subcutaneously, once every 2 weeks (Q2W;",
      "n = 263) or once every 4 weeks (Q4W; n = 137), with protocol-driven",
      "titration. The publication does not tabulate the administered doses;",
      "the nine individual fits in Chanu 2010 Figure 4 report STARTING doses of",
      "19, 24 and 27 ug (the three AMICUS, ESA-naive panels) and 60, 100, 120,",
      "120, 120 and 200 ug (the six MAXIMA / PROTOS, ESA-treated panels)."
    ),
    regions        = "Not stated; three multicenter phase III studies (AMICUS, MAXIMA, PROTOS)",
    notes          = paste(
      "PK analysis population 400 patients with 4554 measurable serum",
      "concentrations; PD analysis population 400 patients with 10 089",
      "evaluable hemoglobin assessments (Chanu 2010 Table II). Serum C.E.R.A.",
      "was measured by a validated ELISA with a limit of quantification of 150",
      "pg/mL (0.15 ng/mL); the first below-LOQ value per subject was imputed as",
      "LOQ/2 and subsequent below-LOQ values were ignored. Estimation used",
      "NONMEM version 5 level 1 (FOCE INTER where practical, otherwise FO), in",
      "a sequential PK-then-PD fit: individual post hoc PK parameters from the",
      "PK model were used as the input to the PD model. Eight concentrations",
      "with |WRES| > 8 were dropped during PK model building. Demographics from",
      "Table I; study designs from Table II."
    )
  )

  ini({
    # === PK structural parameters (Chanu 2010 Table III) =====================
    # The paper reports back-transformed values on the natural scale, so each
    # entry is wrapped in log() here.
    lcl     <- log(0.749);        label("Typical clearance at WT 70 kg (L/day)")                                  # Chanu 2010 Table III: CL = 0.749 L/d (RSE 2%)
    lvc     <- log(4.72);         label("Typical volume of distribution at WT 70 kg and AGE 60 years (L)")        # Chanu 2010 Table III: V = 4.72 L (RSE 2%)
    lka     <- fixed(log(0.825)); label("First-order SC absorption rate constant (1/day)")                        # Chanu 2010 Table III: ka = 0.825 1/d, footnote a 'Fixed'
    lfdepot <- log(0.394);        label("Absolute bioavailability after SC dosing (fraction)")                    # Chanu 2010 Table III: F = 0.394 (RSE 4%)

    # PK covariate effects (Chanu 2010 Table III; log-linear / normalized power
    # parameterization of Methods 'Covariate Model (Stage 2)':
    # TVP = theta_P * (COV / rCOV)^theta_COV, combined multiplicatively).
    e_wt_cl  <- 0.571; label("Power exponent of (WT/70 kg) on clearance (unitless)")                              # Chanu 2010 Table III: Effect of BW on CL = 0.571 (RSE 13%)
    e_wt_vc  <- 0.443; label("Power exponent of (WT/70 kg) on volume of distribution (unitless)")                 # Chanu 2010 Table III: Effect of BW on V = 0.443 (RSE 17%)
    e_age_vc <- 0.267; label("Power exponent of (AGE/60 years) on volume of distribution (unitless)")             # Chanu 2010 Table III: Effect of age on V = 0.267 (RSE 19%)

    # === PK interindividual variability (Chanu 2010 Table III) ===============
    # Table III reports IIV as CV%. Those percentages are omega itself expressed
    # as a percentage (omega^2 = (CV%/100)^2), NOT the exact log-normal CV. The
    # discriminator is Chanu 2010 Figure 3, which plots every post hoc PD
    # parameter: SC50 spans the whole 10^-4 to 10^3 ng/mL axis with points at
    # both limits, i.e. at least 7 orders of magnitude. Reading CV% through the
    # exact relation omega^2 = log(CV^2 + 1) would give omega(SC50) = 1.86, a
    # prior 95% interval of only 3.2 decades, which post hoc estimates cannot
    # spread beyond; reading it as omega = 5.59 gives 9.5 decades, which the
    # figure is consistent with after shrinkage. The Smax panel agrees: its post
    # hoc points run from about 0.02 to about 15 (roughly 5.8 omega for n = 400),
    # implying omega near 1.2, against 1.42 on this convention and 1.05 on the
    # exact one. Both PD panels imply the same ~40-45% shrinkage on this reading
    # and mutually inconsistent shrinkage on the other. The same convention is
    # then applied to the PK rows for consistency (below CV ~ 30% the two
    # readings differ by under 4% anyway). See vignette 'Assumptions and
    # deviations'.
    etalcl ~ 0.0784        # Chanu 2010 Table III: IIV CL CV% = 28 (RSE 9%) -> 0.28^2
    etalvc ~ 0.0729        # Chanu 2010 Table III: IIV V  CV% = 27 (RSE 11%) -> 0.27^2
    etalka ~ fix(0.6724)   # Chanu 2010 Table III: IIV ka CV% = 82, per footnote a -> 0.82^2
    # Chanu 2010 Table III reports IIV on F as CV% = 0, fixed. A zero-variance
    # eta cannot be carried in ini(), so F has no eta in this implementation.

    # === PK interoccasion variability (Chanu 2010 Table III) =================
    # IOV was estimated on CL only in the final model (Table III lists a single
    # 'Random effects IOV' row, CL CV% = 9). Chanu 2010 Methods also mentions
    # preliminary IOV on F, which the final model does not retain. Encoded as
    # the occasion-indicator expansion rather than the eta ~ var | occ form,
    # which rxode2 parses but cannot simulate; occasions 2 and 3 are fixed equal
    # to occasion 1, the equivalent of a NONMEM $OMEGA BLOCK(1) SAME.
    etaiov_cl_1 ~ 0.0081        # Chanu 2010 Table III: IOV CL CV% = 9 (RSE 32%) -> 0.09^2
    etaiov_cl_2 ~ fix(0.0081)   # Chanu 2010 Table III: SAME-equivalent
    etaiov_cl_3 ~ fix(0.0081)   # Chanu 2010 Table III: SAME-equivalent

    # === PK residual error (Chanu 2010 Table III and Results) ================
    # The published final PK error model is the DOUBLE-EXPONENTIAL form of
    # Chanu 2010 Results: Cij = C*ij * exp(eps1ij) + m * exp(eps2ij), with
    # m = 0.150 ng/mL fixed (Table III), var(eps1) = sigma1^2 = 0.141 and
    # var(eps2) = sigma2^2 = 0.691. rxode2 5.1.7 parses lnorm() + add() but
    # rxSolve() cannot solve it, so the two exponential terms are encoded as the
    # moment-matched proportional and additive terms below:
    #   propSd = sqrt(exp(0.141) - 1)                     = 0.3891
    #   addSd  = 0.150 * sqrt(exp(0.691) * (exp(0.691) - 1)) = 0.2114
    # This preserves the SD of each term. It does NOT preserve two other
    # features of the published model: the log-normal shape of both terms, and
    # the positive mean 0.150 * exp(0.691/2) = 0.212 ng/mL that the published
    # additive term carries (the point of replacing the mean-zero additive term
    # was to remove a bias at the lowest concentrations, near the 0.15 ng/mL
    # LOQ). Both deviations matter only close to the LOQ. See vignette
    # 'Assumptions and deviations'.
    propSd <- 0.3891; label("Proportional residual error on serum C.E.R.A. (fraction)")                           # Chanu 2010 Table III: sigma1^2 (proportional) = 0.141 (RSE 6%)
    addSd  <- 0.2114; label("Additive residual error on serum C.E.R.A. (ng/mL)")                                  # Chanu 2010 Table III: m = 0.150 ng/mL fixed, sigma2^2 (additive) = 0.691 (RSE 12%)

    # === PD structural parameters (Chanu 2010 Table IV) ======================
    lemax    <- log(0.425);       label("Smax, maximum fractional increase in hemoglobin production rate (unitless)")  # Chanu 2010 Table IV: Smax = 0.425 (RSE 13%)
    lec50    <- log(0.898);       label("SC50 at CRP 5 mg/L and ESAD 7000 IU/week (ng/mL)")                            # Chanu 2010 Table IV: SC50 = 0.898 ng/mL (RSE 34%)
    lmtt_rbc <- fixed(log(61.3)); label("LS, apparent red-blood-cell life span (day)")                                 # Chanu 2010 Table IV: LS = 61.3 d, footnote a 'Fixed'
    lrbase   <- fixed(log(9.30)); label("Hb0, hemoglobin concentration before any ESA treatment (g/dL)")               # Chanu 2010 Table IV: Hb0 = 9.30 g/dL, footnote a 'Fixed'

    # PD covariate effects (Chanu 2010 Table IV; same normalized power
    # parameterization as the PK covariates, per Methods 'Covariate Model
    # (Stage 2)' for the PD analysis: 'The parameterization and the procedure of
    # covariates selection used were similar to that described for the
    # population PK analysis').
    e_crp_ec50  <- 0.319; label("Power exponent of (CRP/5 mg/L) on SC50 (unitless)")                              # Chanu 2010 Table IV: Effect of CRP on SC50 = 0.319 (RSE 52%)
    e_esad_ec50 <- 0.303; label("Power exponent of (ESAD/7000 IU/week) on SC50 (unitless)")                       # Chanu 2010 Table IV: Effect of DEPO on SC50 = 0.303 (RSE 36%)

    # === PD interindividual variability (Chanu 2010 Table IV) ================
    # Same CV% convention as the PK block above; see that comment for the
    # Figure 3 evidence. The Smax and SC50 variances are enormous by design --
    # the paper itself warns that 'given the very high interindividual
    # variability, particularly in SC50 ... results from the model should be
    # interpreted with some caution when used for simulations'.
    etalemax    ~ 2.0164   # Chanu 2010 Table IV: IIV Smax CV% = 142 (RSE 40%) -> 1.42^2
    etalec50    ~ 31.2481  # Chanu 2010 Table IV: IIV SC50 CV% = 559 (RSE 34%) -> 5.59^2
    etalmtt_rbc ~ 0.1024   # Chanu 2010 Table IV: IIV LS  CV% = 32 (RSE 45%) -> 0.32^2
    etalrbase   ~ 0.0625   # Chanu 2010 Table IV: IIV Hb0 CV% = 25 (RSE 53%) -> 0.25^2

    # === PD residual error (Chanu 2010 Table IV) =============================
    # Additive on hemoglobin: Hbij = Hb*ij + epsij, var(eps) = sigma^2 = 0.357.
    addSd_hb <- 0.5975; label("Additive residual error on hemoglobin (g/dL)")                                     # Chanu 2010 Table IV: sigma^2 = 0.357 (RSE 4%) -> sqrt(0.357)
  })

  model({
    # =====================================================================
    # PK: one compartment, first-order SC absorption, first-order
    # elimination (Chanu 2010 Methods 'Structural PK Model (Stage 1)').
    # The paper writes the closed-form solutions for SC and IV dosing;
    # the equivalent ODE system is used here so that both routes, repeat
    # dosing and the PD coupling are handled by one model.
    # =====================================================================

    # Interoccasion variability on CL, occasion-indicator expansion.
    oc1    <- (OCC == 1)
    oc2    <- (OCC == 2)
    oc3    <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual PK parameters. Covariate effects are multiplicative
    # normalized power terms (Chanu 2010 Methods: for two covariates on one
    # parameter, TVP = theta_P * (AGE/rAGE)^theta_AGE * (BW/rBW)^theta_BW).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (AGE / 60)^e_age_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # SC doses enter depot and carry the absolute bioavailability F; IV doses
    # enter central directly and use the default f(central) = 1.
    f(depot) <- exp(lfdepot)

    # Dose in ug, central in ug, vc in L -> ug/L = ng/mL.
    Cc <- central / vc

    # =====================================================================
    # PD: red-blood-cell life-span indirect-response model
    # (Chanu 2010 Methods 'Structural PD Model (Stage 1)'):
    #   Hb'(t) = S(t) - S(t - LS) + SESA
    #   S(t)   = Hb0/LS * (1 + E(C(t)))
    #   E(C)   = Smax * C / (SC50 + C)
    #   SESA   = (Hb0 - Hbsw)/LS   for 0 <= t <= LS, else 0
    # =====================================================================

    ls   <- exp(lmtt_rbc + etalmtt_rbc)
    hb0  <- exp(lrbase + etalrbase)
    emax <- exp(lemax + etalemax)

    # ESA-naive patients (AMICUS) have no prior epoetin dose. Chanu 2010 does
    # not print the DEPO value used for them; a normalized power term is
    # undefined at DEPO = 0, so the factor is gated to 1, the same idiom the
    # register records for ESAD in Naik_2013_peginesatide.R.
    esadf      <- (ESAD > 0)
    esad_ratio <- (ESAD * esadf + 7000 * (1 - esadf)) / 7000

    ec50 <- exp(lec50 + etalec50) * (CRP / 5)^e_crp_ec50 * esad_ratio^e_esad_ec50

    # Hemoglobin at the switch. For a previously ESA-treated patient this is the
    # observed Hbsw covariate; for an ESA-naive patient the switch value IS the
    # individual's pre-ESA baseline, so SESA vanishes identically.
    hbsw <- esadf * HGB_BL + (1 - esadf) * hb0

    # Delayed serum concentration. delay() interpolates the central state at
    # t - LS from the solver's dense output; before the start of integration the
    # constant initial condition central(0) = 0 is used, which is exactly the
    # paper's assumption that baseline production Hb0/LS 'was considered to be
    # constant during one LS before time of the Hb0 assessment'.
    cpast <- delay(central, ls) / vc

    snow  <- hb0 / ls * (1 + emax * Cc    / (ec50 + Cc))
    spast <- hb0 / ls * (1 + emax * cpast / (ec50 + cpast))

    # The 'virtual loss' of the previous ESA effect, valid only during the first
    # life span after the switch.
    sesa <- (t <= ls) * (hb0 - hbsw) / ls

    d/dt(hb) <- snow - spast + sesa

    hb(0) <- hbsw

    Cc ~ prop(propSd) + add(addSd)
    hb ~ add(addSd_hb)
  })
}
