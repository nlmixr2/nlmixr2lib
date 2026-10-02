Hornik_2019_methylprednisolone_il6 <- function() {
  description <- "Two-compartment population PK model of methylprednisolone with first-order formation from its sodium succinate prodrug, linked to an indirect-response model of interleukin-6 (IL-6) in which methylprednisolone inhibits and cardiopulmonary bypass (CPB) stimulates IL-6 production with partial drug-CPB interaction, in neonates undergoing cardiac surgery on CPB (Hornik 2019)"
  reference <- "Hornik CP, Gonzalez D, Dumond J, Wu H, Graham EM, Hill KD, Cohen-Wolkowiez M. Population Pharmacokinetic/Pharmacodynamic Modeling of Methylprednisolone in Neonates Undergoing Cardiopulmonary Bypass. CPT Pharmacometrics Syst Pharmacol. 2019;8(12):913-922. doi:10.1002/psp4.12470"
  vignette <- "Hornik_2019_methylprednisolone"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (methylprednisolone, Cc); pg/mL (IL-6, il6)"
  )

  # Two bookkeeping states carry the CPB time course that the NONMEM control
  # stream (Data S2) computed from record flags and the CPBSTART data column:
  # cpb_elapsed integrates time on CPB so the 0.5 h onset delay (CPBES) can be
  # read off the solve, and cpb_decay is the post-CPB exponential withdrawal
  # exp(-0.693 * TACPB / CPBH) of Eq. 6, written as a first-order decay that
  # starts at 1 and only runs while CPB_POST = 1.
  paper_specific_compartments <- c("cpb_elapsed", "cpb_decay")

  compartmentData <- list(
    depot = list(
      analyte = "methylprednisolone sodium succinate (prodrug)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "methylprednisolone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methylprednisolone", units = "mg", specimen = "tissue", verified = TRUE),
    il6 = list(analyte = "interleukin-6", units = "pg/mL", specimen = "plasma", verified = TRUE),
    cpb_elapsed = list(
      analyte = "time on cardiopulmonary bypass",
      units = "h",
      specimen = "not applicable",
      verified = TRUE
    ),
    cpb_decay = list(
      analyte = "post-bypass CPB-effect withdrawal fraction",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Weight at the time of the first sample (Table 1 median 3.2 kg, range 2.2-4.3). Power scaling normalized to the 3.2 kg cohort median: estimated exponent 1.24 on CL and Q, exponent 1 (not estimated) on Vc and Vp (Methods; Table 2; Data S1-S3 LSV = (WTKG/3.2)**1).",
      source_name = "WTKG"
    ),
    T_CPB = list(
      description = "Total cardiopulmonary bypass duration",
      units = "minutes",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Power effect on post-CPB clearance normalized to the 156.5 min cohort median (Table 1; Table 2 footnote). In Data S2 the column also sets the end of CPB (TACPB = TIME - CPBSTART - CPBTIME/60), which in this model is carried by the CPB_ON / CPB_POST indicators instead.",
      source_name = "CPBTIME"
    ),
    CPB_ON = list(
      description = "Cardiopulmonary bypass phase indicator (on the bypass circuit)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-CPB or post-CPB)",
      notes = "Time-varying. 1 from the start of CPB until separation from the circuit, 0 otherwise (the paper's ONCPB = STRT x ENDT, Eq. 6). Hornik 2019 does not resolve a rewarming sub-phase, so CPB_ON = 1 for the whole bypass run and CPB_REWARM is not used. Supply records at the CPB start and end times so the indicator switches at the right moments. The 0.5 h onset delay (CPBES) is applied inside the model from the elapsed time on CPB, so CPB_ON must switch to 1 at the actual start of CPB, not 0.5 h later.",
      source_name = "ONCPB (STRT, ENDT; derived from the POINT sample index in Data S2)"
    ),
    CPB_POST = list(
      description = "Post-cardiopulmonary-bypass period indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-CPB or on CPB)",
      notes = "Time-varying. 1 after separation from CPB, 0 before. Selects the post-CPB clearance (the paper's POSTCPB, Table 2 footnote) and switches the IL-6 CPB effect from the on-CPB plateau to its exponential withdrawal (the paper's 1 - ENDT, Eq. 6). The control streams code POSTCPB as STA4 = STA2 + STA3 (study-period flags in the dataset); the paper names the indicator POSTCPB and describes CL as pre-CPB versus post-CPB, which is the reading used here.",
      source_name = "POSTCPB (STA2 + STA3 in Data S1-S3); 1 - ENDT"
    ),
    RACHS1 = list(
      description = "Risk Adjustment for Congenital Heart Surgery 1 (RACHS-1) category",
      units = "(categorical; 1-6 integer)",
      type = "categorical",
      reference_category = "RACHS-1 < 4",
      notes = "Time-fixed. Decomposed inside model() into rachs1_high = (RACHS1 >= 4), which multiplies the IL-6 CPB effect by 2.59 (Results; Table 3; Data S2 RANKN = 1 if RANK >= 4). 42% of the cohort had RACHS-1 < 4 and 58% had RACHS-1 >= 4 (Table 1).",
      source_name = "RANK"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      notes = "Screened on CL (maturation, power and linear forms; Table S1) and on the IL-6 parameters (Figure S8; Table S2); not retained in the IL-6 model. It is retained on the IL-10 CPB effect (Hornik_2019_methylprednisolone_il10)."
    ),
    PNA = list(
      description = "Postnatal age",
      units = "days",
      type = "continuous",
      notes = "Screened on CL (Table S1) and on IL-6 baseline (Figure S8; Table S2); not retained."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened on CL (Table S1) and on IL-6 baseline (Figure S8; Table S2); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 64,
    n_subjects_pd = "62 neonates contributed 314 IL-6 observations (subjects 25 and 31 excluded from the IL-6 fit, Data S2)",
    n_studies = 1,
    age_range = "postnatal age 3-30 days at first sample (median 7); gestational age at birth 34.6-42 weeks (median 39); postmenstrual age 36-44 weeks (median 40)",
    weight_range = "2.2-4.3 kg (median 3.2)",
    sex_female_pct = 47,
    race_ethnicity = "White 58%, Black 25%, Latino 13%, Asian 2%, Latino/black 2%, Latino/white 2% (Table 1)",
    disease_state = "neonates with congenital heart disease undergoing cardiac surgery on cardiopulmonary bypass (median CPB time 156.5 min, range 64-251); RACHS-1 < 4 in 42%, >= 4 in 58%",
    dose_range = "methylprednisolone sodium succinate 30 mg/kg IV over 1 h, as one dose at CPB induction (29 neonates) or two doses (~10 h before CPB and at CPB induction; 35 neonates)",
    regions = "USA",
    notes = "Prospective randomized trial NCT00934843 (Methods). 290 methylprednisolone concentrations (1.07-12,700 ng/mL, none below the 1 ng/mL LLOQ); median IL-6 baseline 9.5 pg/mL (0.7-83.3). PK and PD were fit sequentially: individual PK estimates were fixed when the PD parameters were estimated."
  )

  ini({
    # --- Population PK (Table 2; Data S3 control stream) -------------------
    # Typical values for a 3.2 kg neonate. The paper labels them CL/F, Vc/F,
    # Q/F, Vp/F because the drug is dosed as the prodrug and F1 is fixed to 1
    # (Data S3 TVF1 = 1); they are encoded as CL, Vc, Q, Vp with the prodrug
    # dose entering the depot state.
    lcl <- log(3.88); label("Clearance pre-CPB for a 3.2 kg neonate (L/h)") # Table 2 'CL/F (L/hour, 3.2 kg)' = 3.88
    lvc <- log(8.92); label("Central volume of distribution for a 3.2 kg neonate (L)") # Table 2 'Vc/F (L, 3.2 kg)' = 8.92
    lq <- log(0.10); label("Intercompartmental clearance for a 3.2 kg neonate (L/h)") # Table 2 'Q/F (L/hour, 3.2 kg)' = 0.10
    lvp <- log(16.81); label("Peripheral volume of distribution for a 3.2 kg neonate (L)") # Table 2 'Vp/F (L, 3.2 kg)' = 16.81
    lka <- log(0.41); label("Formation rate constant of methylprednisolone from the sodium succinate prodrug, Kf (1/h)") # Table 2 'Kf (1/hour)' = 0.41

    e_wt_cl_q <- 1.24; label("Power exponent of (WT / 3.2 kg) on CL and Q (unitless)") # Table 2 'WT exponent on CL and Q' = 1.24
    e_wt_vc_vp <- fixed(1); label("Power exponent of (WT / 3.2 kg) on Vc and Vp (unitless)") # Methods; Data S3 LSV = (WTKG/3.2)**1, not estimated
    e_t_cpb_cl <- -0.47; label("Power exponent of (T_CPB / 156.5 min) on post-CPB CL (unitless)") # Table 2 'CPB time on CL' = -0.47

    # IIV reported as %CV (Table 2); omega^2 = log(1 + CV^2). No IIV on Vp or
    # Kf (Data S3 V3 = TVV3, KA = THETA(5)).
    etalcl ~ 0.201130 # Table 2 'IIV, CL' 47.2 %CV
    etalvc ~ 0.067374 # Table 2 'IIV, Vc' 26.4 %CV
    etalq ~ 0.100999 # Table 2 'IIV, Q' 32.6 %CV

    propSd <- 0.428; label("Proportional residual error, methylprednisolone (fraction)") # Table 2 'Proportional error,%' = 42.8; Data S3 Y = IPRED + IPRED*ERR(1)

    # --- IL-6 indirect response (Table 3; Data S2 control stream) ----------
    limax <- fixed(log(1)); label("Maximum fractional inhibition of IL-6 production by methylprednisolone, Imax (fraction)") # Table 3 'I max' = 1 FIX
    lic50 <- log(14); label("Methylprednisolone concentration giving half-maximal inhibition of IL-6 production, IC50 (ng/mL)") # Table 3 'IC 50, ng/mL' = 14
    lrbase <- log(7.9); label("Baseline IL-6 concentration before the first dose, IL-6base (pg/mL)") # Table 3 'IL-6 base, pg/mL' = 7.9
    lkout <- log(0.171); label("First-order IL-6 elimination rate constant, Rout (1/h)") # Table 3 'R out, 1/hour' = 0.171
    lhill <- log(2.53); label("Hill coefficient of the methylprednisolone inhibition (unitless)") # Table 3 'HILL' = 2.53
    lcpbe <- log(48.6); label("Fold increase in IL-6 production during CPB for RACHS-1 below 4, CPBE (unitless)") # Table 3 'CPBE' = 48.6
    pct_cpb_noint <- 21.4; label("Percentage of the CPB effect on IL-6 production that does not interact with methylprednisolone, PER (percent)") # Table 3 'Percent of CPB effect not interacting with MP' = 21.4
    lthalf_cpb <- log(9.08); label("Half-life of the post-CPB withdrawal of the CPB effect on IL-6 production, CPBH (h)") # Table 3 'CPB effect half-life, hour' = 9.08
    e_rachs1_cpbe <- 2.59; label("Multiplicative factor on CPBE for RACHS-1 of 4 or higher (unitless)") # Table 3 'RACHS-1 >= 4 on CPBE' = 2.59; Results CPBE = 48.6 x 2.59^(RACHS-1 >= 4)

    # IIV reported as %CV (Table 3); omega^2 = log(1 + CV^2). Data S2 carries
    # ETA(3) on BASE and ETA(6) on CPBEFFECT; the other PD etas are held at 0.
    etalrbase ~ 0.698147 # Table 3 'IIV, IL-6 base' 100.5 %CV
    etalcpbe ~ 0.529979 # Table 3 'IIV, CPBE' 83.6 %CV

    propSd_il6 <- 0.541; label("Proportional residual error, IL-6 (fraction)") # Table 3 'Proportional error,%' = 54.1; Data S2 Y = EFF*EXP(ERR(1)) under FOCE-I
  })
  model({
    # --- PK (Data S3 $PK; Table 2 footnote) ---------------------------------
    allom_clq <- (WT / 3.2)^e_wt_cl_q
    allom_v <- (WT / 3.2)^e_wt_vc_vp

    # CL = (3.88 x (1 - POSTCPB) + 3.88 x POSTCPB x (CPBtime/156.5)^-0.47) x (WT/3.2)^1.24
    cl <- exp(lcl + etalcl) *
      ((1 - CPB_POST) + CPB_POST * (T_CPB / 156.5)^e_t_cpb_cl) *
      allom_clq
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq + etalq) * allom_clq
    vp <- exp(lvp) * allom_v
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The prodrug dose (mg) enters depot and is converted to methylprednisolone
    # at rate Kf with F1 = 1 (Data S3; Discussion).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L -> ng/mL (Data S3 S2 = V2/1000)
    Cc <- 1000 * central / vc

    # --- CPB time course (Eq. 6; Data S2) -----------------------------------
    # Elapsed time on CPB; the CPB effect switches on 0.5 h after CPB starts
    # (CPBES = 1 if TACPBS > 0.5).
    d/dt(cpb_elapsed) <- CPB_ON
    cpbes <- 0
    if (cpb_elapsed > 0.5) {
      cpbes <- 1
    }

    # Post-CPB withdrawal exp(-0.693 * TACPB / CPBH): 1 until CPB ends, then
    # first-order decay with half-life CPBH (0.693 as written in Data S2).
    thalf_cpb <- exp(lthalf_cpb)
    cpb_decay(0) <- 1
    d/dt(cpb_decay) <- -(0.693 / thalf_cpb) * cpb_decay * CPB_POST

    # CPBV = 1 x ONCPB x CPBES + (1 - ENDT) x EXP(-0.693 x TACPB / CPBH)
    cpbv <- CPB_ON * cpbes + CPB_POST * cpb_decay

    rachs1_high <- 0
    if (RACHS1 >= 4) {
      rachs1_high <- 1
    }
    cpbe <- exp(lcpbe + etalcpbe) * e_rachs1_cpbe^rachs1_high
    cpbfx <- cpbe * cpbv # Eq. 6 CPBFX

    # --- IL-6 indirect response (Eqs. 3-5; Data S2 $DES) ---------------------
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    kin <- kout * rbase # Eq. 4

    dfx <- imax * Cc^hill / (ic50^hill + Cc^hill) # Eq. 5

    # Eq. 3 (partial interaction). The paper's STRT switch is not needed here:
    # before CPB starts cpbfx = 0 and both branches reduce to kin * (1 - DFX).
    fnoint <- pct_cpb_noint / 100
    il6(0) <- rbase
    d/dt(il6) <- kin * ((1 + cpbfx * (1 - fnoint)) * (1 - dfx) + cpbfx * fnoint) - kout * il6

    Cc ~ prop(propSd)
    il6 ~ prop(propSd_il6)
  })
}
