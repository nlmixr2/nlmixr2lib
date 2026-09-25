Hornik_2019_methylprednisolone_il10 <- function() {
  description <- "Two-compartment population PK model of methylprednisolone with first-order formation from its sodium succinate prodrug, linked to an indirect-response model of interleukin-10 (IL-10) in which methylprednisolone and cardiopulmonary bypass (CPB) both stimulate IL-10 production with complete (multiplicative) drug-CPB interaction, in neonates undergoing cardiac surgery on CPB (Hornik 2019)"
  reference <- "Hornik CP, Gonzalez D, Dumond J, Wu H, Graham EM, Hill KD, Cohen-Wolkowiez M. Population Pharmacokinetic/Pharmacodynamic Modeling of Methylprednisolone in Neonates Undergoing Cardiopulmonary Bypass. CPT Pharmacometrics Syst Pharmacol. 2019;8(12):913-922. doi:10.1002/psp4.12470"
  vignette <- "Hornik_2019_methylprednisolone"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (methylprednisolone, Cc); pg/mL (IL-10, il10)"
  )

  # il10 is the IL-10 indirect-response state (the IL-10 sibling of the
  # canonical il6 state). cpb_elapsed is a bookkeeping state that integrates
  # time on CPB so the 0.5 h onset delay (CPBES) of the NONMEM control stream
  # (Data S1) can be read off the solve. Unlike IL-6, the IL-10 model has no
  # post-CPB withdrawal term: Data S1 comments it out
  # (CPBV = 1*ONCPB*CPBES;+(1-ENDT)*EXP(-0.693*TACPB/CPBH)), so the CPB effect
  # stops when the patient comes off bypass.
  paper_specific_compartments <- c("il10", "cpb_elapsed")

  compartmentData <- list(
    depot = list(
      analyte = "methylprednisolone sodium succinate (prodrug)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "methylprednisolone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methylprednisolone", units = "mg", specimen = "tissue", verified = TRUE),
    il10 = list(analyte = "interleukin-10", units = "pg/mL", specimen = "plasma", verified = TRUE),
    cpb_elapsed = list(
      analyte = "time on cardiopulmonary bypass",
      units = "h",
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
      notes = "Time-fixed per subject. Power effect on post-CPB clearance normalized to the 156.5 min cohort median (Table 1; Table 2 footnote).",
      source_name = "CPBTIME"
    ),
    CPB_ON = list(
      description = "Cardiopulmonary bypass phase indicator (on the bypass circuit)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-CPB or post-CPB)",
      notes = "Time-varying. 1 from the start of CPB until separation from the circuit, 0 otherwise (the paper's ONCPB = STRT x ENDT, Eq. 6). Hornik 2019 does not resolve a rewarming sub-phase, so CPB_ON = 1 for the whole bypass run and CPB_REWARM is not used. Supply records at the CPB start and end times so the indicator switches at the right moments. The 0.5 h onset delay (CPBES) is applied inside the model from the elapsed time on CPB, so CPB_ON must switch to 1 at the actual start of CPB, not 0.5 h later.",
      source_name = "ONCPB (STRT, ENDT; derived from the POINT sample index in Data S1)"
    ),
    CPB_POST = list(
      description = "Post-cardiopulmonary-bypass period indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-CPB or on CPB)",
      notes = "Time-varying. 1 after separation from CPB, 0 before. Selects the post-CPB clearance only (the paper's POSTCPB, Table 2 footnote); the IL-10 CPB effect has no post-CPB term. The control streams code POSTCPB as STA4 = STA2 + STA3 (study-period flags in the dataset); the paper names the indicator POSTCPB and describes CL as pre-CPB versus post-CPB, which is the reading used here.",
      source_name = "POSTCPB (STA2 + STA3 in Data S1-S3)"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed: postmenstrual age (GA + PNA) at the time of the first sample (Supplementary Methods; Table 1 median 40 weeks, range 36-44). Power effect on the IL-10 CPB effect normalized to 40 weeks, CPBE = 45.7 x (PMA/40)^14.8 (Results; Table 4; Data S1). Weeks, not the register's default months.",
      source_name = "PMAW"
    )
  )

  covariatesDataExcluded <- list(
    RACHS1 = list(
      description = "Risk Adjustment for Congenital Heart Surgery 1 (RACHS-1) category",
      units = "(categorical; 1-6 integer)",
      type = "categorical",
      notes = "Screened on the IL-10 parameters (Table S3); not retained. It is retained on the IL-6 CPB effect (Hornik_2019_methylprednisolone_il6)."
    ),
    PNA = list(
      description = "Postnatal age",
      units = "days",
      type = "continuous",
      notes = "Screened on CL (Table S1) and on IL-10 baseline (Figure S9; Table S3); not retained."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened on CL (Table S1) and on IL-10 baseline (Figure S9; Table S3); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 64,
    n_subjects_pd = "64 neonates contributed 324 IL-10 observations",
    n_studies = 1,
    age_range = "postnatal age 3-30 days at first sample (median 7); gestational age at birth 34.6-42 weeks (median 39); postmenstrual age 36-44 weeks (median 40)",
    weight_range = "2.2-4.3 kg (median 3.2)",
    sex_female_pct = 47,
    race_ethnicity = "White 58%, Black 25%, Latino 13%, Asian 2%, Latino/black 2%, Latino/white 2% (Table 1)",
    disease_state = "neonates with congenital heart disease undergoing cardiac surgery on cardiopulmonary bypass (median CPB time 156.5 min, range 64-251); RACHS-1 < 4 in 42%, >= 4 in 58%",
    dose_range = "methylprednisolone sodium succinate 30 mg/kg IV over 1 h, as one dose at CPB induction (29 neonates) or two doses (~10 h before CPB and at CPB induction; 35 neonates)",
    regions = "USA",
    notes = "Prospective randomized trial NCT00934843 (Methods). 290 methylprednisolone concentrations (1.07-12,700 ng/mL, none below the 1 ng/mL LLOQ); median IL-10 baseline 1.3 pg/mL (0.1-9.5). PK and PD were fit sequentially: individual PK estimates were fixed when the PD parameters were estimated."
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

    # --- IL-10 indirect response (Table 4; Data S1 control stream) ---------
    lsmax <- log(2.28); label("Maximum fold stimulation of IL-10 production by methylprednisolone, Smax (unitless)") # Table 4 'S max' = 2.28
    lsc50 <- log(58.2); label("Methylprednisolone concentration giving half-maximal stimulation of IL-10 production, SC50 (ng/mL)") # Table 4 'SC 50, ng/mL' = 58.2
    lrbase <- log(1.52); label("Baseline IL-10 concentration before the first dose, IL-10base (pg/mL)") # Table 4 'IL-10 base, pg/mL' = 1.52
    lkout <- log(0.542); label("First-order IL-10 elimination rate constant, Rout (1/h)") # Table 4 'R out, 1/hour' = 0.542
    lhill <- log(3.58); label("Hill coefficient of the methylprednisolone stimulation (unitless)") # Table 4 'HILL' = 3.58
    lcpbe <- log(45.7); label("Fold increase in IL-10 production during CPB at a postmenstrual age of 40 weeks, CPBE (unitless)") # Table 4 'CPBE' = 45.7
    e_page_cpbe <- 14.8; label("Power exponent of (PMA / 40 weeks) on CPBE (unitless)") # Table 4 'PMA on CPBE' = 14.8; Results CPBE = 45.7 x (PMA/40)^14.8

    # IIV reported as %CV (Table 4); omega^2 = log(1 + CV^2). Data S1 carries
    # ETA(1) on EMAX (Smax), ETA(3) on BASE and ETA(6) on CPBEFFECT; the other
    # PD etas are held at 0.
    etalsmax ~ 0.792993 # Table 4 'IIV, S max' 110 %CV
    etalrbase ~ 0.349677 # Table 4 'IIV, IL-10 base' 64.7 %CV
    etalcpbe ~ 0.574454 # Table 4 'IIV, CPBE' 88.1 %CV

    propSd_il10 <- 0.538; label("Proportional residual error, IL-10 (fraction)") # Table 4 'Proportional error,%' = 53.8; Data S1 Y = EFF*EXP(ERR(1)) under FOCE-I
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

    # --- CPB time course (Eq. 6 without the withdrawal term; Data S1) -------
    # Elapsed time on CPB; the CPB effect switches on 0.5 h after CPB starts
    # (CPBES = 1 if TACPBS > 0.5) and off when CPB ends (CPBV = ONCPB x CPBES).
    d/dt(cpb_elapsed) <- CPB_ON
    cpbes <- 0
    if (cpb_elapsed > 0.5) {
      cpbes <- 1
    }
    cpbv <- CPB_ON * cpbes

    cpbe <- exp(lcpbe + etalcpbe) * (PAGE / 40)^e_page_cpbe
    cpbfx <- cpbe * cpbv

    # --- IL-10 indirect response (Eq. 2 analogue; Data S1 $DES) -------------
    smax <- exp(lsmax + etalsmax)
    sc50 <- exp(lsc50)
    hill <- exp(lhill)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    kin <- kout * rbase # Eq. 4 analogue

    sfx <- smax * Cc^hill / (sc50^hill + Cc^hill)

    # Complete interaction: CPB and drug effects multiply (Results; Data S1
    # DADT(4) = KIN*(1+CPBEFFECT*CPBV)*(1+EMAX*C^HILL/(EC50^HILL+C^HILL)) - KOUT*A(4)).
    il10(0) <- rbase
    d/dt(il10) <- kin * (1 + cpbfx) * (1 + sfx) - kout * il10

    Cc ~ prop(propSd)
    il10 ~ prop(propSd_il10)
  })
}
