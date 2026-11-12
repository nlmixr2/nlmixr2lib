Healy_2025_tramadol <- function() {
  description <- "Joint parent-and-metabolite population PK model for tramadol and its active CYP2D6-formed metabolite O-desmethyltramadol (M1) in neonates and young infants given intravenous tramadol, used to simulate oral dosing in children aged 3 months to 18 years with chronic pain (Healy 2025). Two-compartment tramadol and two-compartment M1 disposition, all clearances allometrically scaled (WT^0.75, volumes WT^1, 70 kg reference) with separate sigmoid postmenstrual-age maturation on the non-M1 tramadol clearance, the M1 formation clearance and the M1 elimination clearance. The tramadol-to-M1 formation clearance takes a separately estimated typical value and IIV for each of three CYP2D6 activity-score groups, and the M1 elimination clearance a separately estimated typical value for each of three OCT1 active-allele groups. First-order oral absorption (ka fixed from the literature) is added for the oral simulations."
  reference <- "Healy P, Allegaert K, Della Pasqua O. Evaluation of the effect of CYP2D6 and OCT1 polymorphisms on the pharmacokinetics of tramadol: Implications for clinical safety and dose rationale in paediatric chronic pain. Br J Clin Pharmacol. 2025;91(2):283-296. doi:10.1111/bcp.16201"
  vignette <- "Healy_2025_tramadol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (Eq. 1) to a 70 kg standard weight: exponent 0.75 on every clearance (CLPO, CLPM, CLMO, QP, QM) and 1 on every volume (VP1, VP2, VM1, VM2). Neonatal fitting cohort median 3.2 kg (0.93-6.14 kg, Table 1).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Postmenstrual age in WEEKS (the source maturation equation and its TM50 values are in weeks). Drives the sigmoid maturation Fmat = 1 / (1 + (PMA / TM50)^-Hill) (Eq. 2) separately on CLPO, CLPM and CLMO. Neonatal fitting cohort median 40 weeks (31-54 weeks, Table 1); simulation cohort 51-976 weeks (Table 3).",
      source_name = "PMA"
    ),
    CYP2D6_ASG2 = list(
      description = "CYP2D6 activity-score group 2 indicator (1 = activity score 1-1.5; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with CYP2D6_ASG3 = 0 (group G1, activity score 0-0.5: poor plus intermediate metabolizers)",
      notes = "Healy 2025 Table 2 grouping G2 (AS 1 and 1.5, extensive metabolizers). Selects the group-specific typical M1 formation clearance CLPM and its group-specific IIV. Paired with CYP2D6_ASG3; both 0 is the G1 reference. The source dataset column CYP (1/2/3) maps as CYP2D6_ASG2 = (CYP == 2).",
      source_name = "CYP"
    ),
    CYP2D6_ASG3 = list(
      description = "CYP2D6 activity-score group 3 indicator (1 = activity score 2-3; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with CYP2D6_ASG2 = 0 (group G1, activity score 0-0.5)",
      notes = "Healy 2025 Table 2 grouping G3 (AS 2 and 3, extensive with AS 2 plus ultrarapid metabolizers). Source dataset column CYP maps as CYP2D6_ASG3 = (CYP == 3).",
      source_name = "CYP"
    ),
    OCT1_ASG1 = list(
      description = "OCT1 (SLC22A1) activity group 1 indicator (1 = one active allele; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with OCT1_ASG2 = 0 (group G0, no active OCT1 allele)",
      notes = "Healy 2025 Table 2 OCT1 grouping G1 (1 active allele). Selects the group-specific typical M1 elimination clearance CLMO. Paired with OCT1_ASG2; both 0 is the G0 reference. The effect is non-monotonic across groups (CLMO 57.7 / 78.7 / 62.2 L/h/70 kg for G0 / G1 / G2), so the two indicators cannot be collapsed into an active-allele count. Source dataset column OCT (0/1/2) maps as OCT1_ASG1 = (OCT == 1).",
      source_name = "OCT"
    ),
    OCT1_ASG2 = list(
      description = "OCT1 (SLC22A1) activity group 2 indicator (1 = two active alleles; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with OCT1_ASG1 = 0 (group G0, no active OCT1 allele)",
      notes = "Healy 2025 Table 2 OCT1 grouping G2 (2 active alleles). Source dataset column OCT maps as OCT1_ASG2 = (OCT == 2).",
      source_name = "OCT"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Carried in the dataset ($INPUT SEX, Table 1: 32 male / 14 female) but no sex effect is in the final model."
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Carried in the dataset ($INPUT PNA); maturation is described on postmenstrual age only."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral dose site for the simulation scenarios (Figure 2 'Ka', 'F'). The model was fitted to intravenous data only; ka is fixed from an external capsule bioavailability study (ref. 28) and no bioavailability value is reported, so F = 1 is implicit. Dose into central to reproduce the intravenous fit."
    ),
    central = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Tramadol central compartment (control stream CENT, A(1); volume VP1). Loses drug by CLPO (other routes) and CLPM (formation of M1) in parallel, so total tramadol clearance CLPP = CLPO + CLPM."
    ),
    peripheral1 = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Tramadol peripheral compartment (control stream PER, A(2); volume VP2, exchange QP)."
    ),
    central_m1 = list(
      analyte = "O-desmethyltramadol (M1)",
      units = "mg (see notes)",
      specimen = "plasma",
      verified = TRUE,
      notes = "M1 central compartment (control stream MCENT, A(3); volume VM1). Fed by CLPM/VP1 * central in tramadol mass units; the outflows to peripheral1_m1 and by CLMO are multiplied by the 0.947 molar ratio (M1 249.3 / tramadol 263.4 g/mol) exactly as the control stream $DES writes them, and the output Cc_m1 = central_m1 / VM1 carries no 0.947 factor ($ERROR S3 = VM1)."
    ),
    peripheral1_m1 = list(
      analyte = "O-desmethyltramadol (M1)",
      units = "mg (see central_m1 notes)",
      specimen = "tissue",
      verified = TRUE,
      notes = "M1 peripheral compartment (control stream MPER, A(4); volume VM2, exchange QM)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 2L,
    age_range = "Postmenstrual age 31-54 weeks (median 40 weeks)",
    weight_range = "0.93-6.14 kg (median 3.2 kg)",
    sex_female_pct = 30.4,
    disease_state = "Neonates and young infants receiving intravenous tramadol for analgesia at University Hospital Leuven, Belgium",
    dose_range = "Intravenous loading dose averaging 2.1 mg/kg over 30 min, then continuous infusion 0.35 mg/kg/h",
    regions = "Belgium (single centre: University Hospital Leuven)",
    genotypes = "CYP2D6 activity-score group G3 (AS 2-3) 2, G2 (AS 1-1.5) 29, G1 (AS 0-0.5) 15; OCT1 group G0 (0 active alleles) 2, G1 (1 active) 20, G2 (2 active) 24 (Table 1)",
    notes = "Neonatal subset (n = 46 of 52 genotyped neonates; 122 patients in the original studies, refs 19-20) with 247 tramadol and 223 M1 concentrations in the final dataset; GC with nitrogen-selective detection, LLOQ 2.0 ng/mL tramadol and 2.5 ng/mL M1. Table 1 reports 32 male / 14 female. The paper's clinical-trial simulations extrapolate the model to oral dosing in a virtual cohort of 3 months to 18 years (PMA 51-976 weeks, 4.7-73.4 kg; Table 3)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Tramadol disposition. Healy 2025 Table 5 point estimates (final
    # estimates; e.g. CLPO 25.6 supersedes the $THETA initial value 30 in the
    # supplement control stream). All clearances in L/h/70 kg, volumes in
    # L/70 kg. QP, VP1 and VP2 carry FIX in the final joint control stream
    # ('Tramadol-specific parameters were then fixed in the subsequent step');
    # Table 5 lists the RSEs from the earlier tramadol-only step.
    # -----------------------------------------------------------------------
    lcl <- log(25.6); label("Tramadol clearance by routes other than M1 formation, CLPO (L/h/70 kg)") # Table 5 'CLPO' = 25.6 (RSE 10%)
    lq <- fixed(log(111)); label("Tramadol intercompartmental clearance QP (L/h/70 kg)") # Table 5 'QP' = 111 (RSE 58.3%); control stream THETA QP '(0, 111) FIX'
    lvc <- fixed(log(141)); label("Tramadol central volume VP1 (L/70 kg)") # Table 5 'VP1' = 141 (RSE 62.2%); control stream THETA VP1 '(0, 141) FIX'
    lvp <- fixed(log(102)); label("Tramadol peripheral volume VP2 (L/70 kg)") # Table 5 'VP2' = 102 (RSE 73.6%); control stream THETA VP2 '(0, 102) FIX'

    # -----------------------------------------------------------------------
    # M1 formation clearance CLPM, one typical value per CYP2D6 activity-score
    # group (Table 2: G1 AS 0-0.5, G2 AS 1-1.5, G3 AS 2-3).
    # -----------------------------------------------------------------------
    lcl_form_m1_asg1 <- log(5.08); label("Tramadol-to-M1 formation clearance CLPM, CYP2D6 group G1 (L/h/70 kg)") # Table 5 'CLPM (CYP2D6-G1)' = 5.08 (RSE 35%)
    lcl_form_m1_asg2 <- log(11.4); label("Tramadol-to-M1 formation clearance CLPM, CYP2D6 group G2 (L/h/70 kg)") # Table 5 'CLPM (CYP2D6-G2)' = 11.4 (RSE 35%)
    lcl_form_m1_asg3 <- log(14.7); label("Tramadol-to-M1 formation clearance CLPM, CYP2D6 group G3 (L/h/70 kg)") # Table 5 'CLPM (CYP2D6-G3)' = 14.7 (RSE 32%)

    # -----------------------------------------------------------------------
    # M1 disposition. CLMO has one typical value per OCT1 group (Table 2:
    # G0 / G1 / G2 = 0 / 1 / 2 active alleles). QM and VM1 are fixed to prior
    # estimates (Table 5 footnote a; VM1 scaled from an IV M1 study in dogs,
    # ref. 27).
    # -----------------------------------------------------------------------
    lcl_m1_asg0 <- log(57.7); label("M1 elimination clearance CLMO, OCT1 group G0 (L/h/70 kg)") # Table 5 'CLMO (OCT1-G0)' = 57.7 (RSE 35%)
    lcl_m1_asg1 <- log(78.7); label("M1 elimination clearance CLMO, OCT1 group G1 (L/h/70 kg)") # Table 5 'CLMO (OCT1-G1)' = 78.7 (RSE 34%)
    lcl_m1_asg2 <- log(62.2); label("M1 elimination clearance CLMO, OCT1 group G2 (L/h/70 kg)") # Table 5 'CLMO (OCT1-G2)' = 62.2 (RSE 32%)
    lq_m1 <- fixed(log(635)); label("M1 intercompartmental clearance QM (L/h/70 kg)") # Table 5 'QM' = 635 (footnote a, fixed)
    lvc_m1 <- fixed(log(78.9)); label("M1 central volume VM1 (L/70 kg)") # Table 5 'VM1' = 78.9 (footnote a, fixed); Methods 'fixed to a value of 78.9 L/70 kg'
    lvp_m1 <- log(101); label("M1 peripheral volume VM2 (L/70 kg)") # Table 5 'VM2' = 101 (RSE 101%)

    # -----------------------------------------------------------------------
    # Sigmoid maturation on postmenstrual age (Eq. 2), all fixed (Table 5
    # footnote a; the same constants are hard-coded in the control stream $PK).
    # -----------------------------------------------------------------------
    ltm50_cl <- fixed(log(39.1)); label("PMA at 50% maturation of CLPO (weeks)") # Table 5 'TM50 CLPO' = 39.1 (fixed)
    lhill_cl <- fixed(log(6.76)); label("Hill coefficient of CLPO maturation (unitless)") # Table 5 'Hill coefficient CLPO' = 6.76 (fixed)
    ltm50_cl_form_m1 <- fixed(log(39.8)); label("PMA at 50% maturation of CLPM (weeks)") # Table 5 'TM50 CLPM' = 39.8 (fixed)
    lhill_cl_form_m1 <- fixed(log(9)); label("Hill coefficient of CLPM maturation (unitless)") # Table 5 'Hill coefficient CLPM' = 9 (fixed)
    ltm50_cl_m1 <- fixed(log(47.7)); label("PMA at 50% maturation of CLMO (weeks)") # Table 5 'TM50 CLMO' = 47.7 (fixed)
    lhill_cl_m1 <- fixed(log(3.4)); label("Hill coefficient of CLMO maturation (unitless)") # Table 5 'Hill coefficient CLMO' = 3.4 (fixed)

    # Allometric exponents (Eq. 1 and the Methods text: 0.75 for clearance,
    # 1 for volume of distribution).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on clearances (unitless)") # Eq. 1 PWR = 0.75 for clearance
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on volumes (unitless)") # Eq. 1 PWR = 1 for volume of distribution

    # Oral absorption for the simulation scenarios (Methods: 'population value
    # of 1.28 h-1 obtained from a bioavailability study of tramadol capsules').
    lka <- fixed(log(1.28)); label("First-order oral absorption rate constant ka (1/h)") # Methods 'Pharmacokinetic modelling and extrapolation': Ka = 1.28 h-1 (ref. 28)

    # -----------------------------------------------------------------------
    # IIV: log-normal (Eq. 3). Table 5 heads the block 'Between-subject
    # variability (CV%)', but the values are the NONMEM $OMEGA variances (the
    # fixed ones match $OMEGA 0.211, 1.98, 1.45 verbatim), so they are used
    # as variances. CLPM carries a separate eta per CYP2D6 group
    # (control stream ETA(3) / ETA(4) / ETA(5)); CLMO has no eta (commented
    # out in the control stream).
    # -----------------------------------------------------------------------
    etalcl ~ 0.161 # Table 5 'CLPO (CYP1)' = 0.161 (RSE 28%)
    etalvc ~ fixed(0.211) # Table 5 'VP1' = 0.211; $OMEGA entry 0.211, not estimated in the joint run
    etalcl_form_m1_asg1 ~ 0.628 # Table 5 'CLPM (CYP1)' = 0.628 (RSE 42%)
    etalcl_form_m1_asg2 ~ 0.0915 # Table 5 'CLPM (CYP2)' = 0.0915 (RSE 47%)
    etalcl_form_m1_asg3 ~ 0.00878 # Table 5 'CLPM (CYP3)' = 0.00878 (RSE 25%)
    etalq_m1 ~ fixed(1.98) # Table 5 'QM' = 1.98 (footnote a: prior estimate from Allegaert et al.); $OMEGA entry 1.98, not estimated
    etalvc_m1 ~ fixed(1.45) # Table 5 'VM1' = 1.45 (footnote a: prior estimate from Allegaert et al.); $OMEGA entry 1.45, not estimated

    # -----------------------------------------------------------------------
    # Residual error: Y = F * (1 + EPS) + EPS_ADD (Eqs. 4-5), variances added
    # (combined2). Each SD is sqrt of the NONMEM $SIGMA variance. The mapping
    # follows the executed control stream ($ERROR: tramadol uses EPS(1),
    # EPS(2); M1 uses EPS(3), EPS(4); $SIGMA 0.0159 FIX, 0.00001 FIX, then two
    # estimated), not the Table 5 row labels, which place 0.00001 on M1 and
    # 0.01 on tramadol; see the vignette Assumptions and deviations.
    # -----------------------------------------------------------------------
    propSd <- fixed(0.126095); label("Proportional residual error for tramadol (fraction)") # $SIGMA(1) '0.0159 FIX' (Table 5 'Proportional error (tramadol)' 0.0159); sqrt(0.0159) = 0.126095
    addSd <- fixed(0.0031623); label("Additive residual error for tramadol (ng/mL)") # $SIGMA(2) '0.00001 FIX' -> EPS(2), the tramadol additive term in $ERROR (Table 5 prints 0.00001 on the M1 proportional row); sqrt(0.00001) = 0.0031623
    propSd_m1 <- 0.1; label("Proportional residual error for M1 (fraction)") # $SIGMA(3), estimated; Table 5 value 0.01 (printed on the 'Additive error (tramadol)' row, RSE 45%); sqrt(0.01) = 0.1
    addSd_m1 <- 0.0063246; label("Additive residual error for M1 (ng/mL)") # $SIGMA(4), estimated; Table 5 'Additive error (M1)' = 0.00004 (RSE 42%); sqrt(0.00004) = 0.0063246
  })

  model({
    # --- Group indicators (Table 2) -----------------------------------------
    cyp2d6_asg1 <- 1 - CYP2D6_ASG2 - CYP2D6_ASG3
    oct1_asg0 <- 1 - OCT1_ASG1 - OCT1_ASG2

    # --- Size (Eq. 1) and maturation (Eq. 2) ---------------------------------
    fsize_cl <- (WT / 70)^e_wt_cl
    fsize_v <- (WT / 70)^e_wt_vc
    fmat_cl <- 1 / (1 + (PAGE / exp(ltm50_cl))^(-exp(lhill_cl)))
    fmat_cl_form_m1 <- 1 / (1 + (PAGE / exp(ltm50_cl_form_m1))^(-exp(lhill_cl_form_m1)))
    fmat_cl_m1 <- 1 / (1 + (PAGE / exp(ltm50_cl_m1))^(-exp(lhill_cl_m1)))

    # --- Tramadol ----------------------------------------------------------
    cl <- exp(lcl + etalcl) * fsize_cl * fmat_cl
    q <- exp(lq) * fsize_cl
    vc <- exp(lvc + etalvc) * fsize_v
    vp <- exp(lvp) * fsize_v

    # --- M1 formation: group-specific typical value AND group-specific eta
    # (control stream: IF (CYP.EQ.k) CLPM = TVCLPM * EXP(ETA(2 + k))) --------
    cl_form_m1_asg1 <- exp(lcl_form_m1_asg1 + etalcl_form_m1_asg1)
    cl_form_m1_asg2 <- exp(lcl_form_m1_asg2 + etalcl_form_m1_asg2)
    cl_form_m1_asg3 <- exp(lcl_form_m1_asg3 + etalcl_form_m1_asg3)
    cl_form_m1 <- (cyp2d6_asg1 * cl_form_m1_asg1 +
      CYP2D6_ASG2 * cl_form_m1_asg2 +
      CYP2D6_ASG3 * cl_form_m1_asg3) *
      fsize_cl *
      fmat_cl_form_m1

    # --- M1 disposition ----------------------------------------------------
    cl_m1 <- (oct1_asg0 * exp(lcl_m1_asg0) +
      OCT1_ASG1 * exp(lcl_m1_asg1) +
      OCT1_ASG2 * exp(lcl_m1_asg2)) *
      fsize_cl *
      fmat_cl_m1
    q_m1 <- exp(lq_m1 + etalq_m1) * fsize_cl
    vc_m1 <- exp(lvc_m1 + etalvc_m1) * fsize_v
    vp_m1 <- exp(lvp_m1) * fsize_v

    # Total tramadol clearance (derived output; CLPP = CLPO + CLPM).
    cl_total <- cl + cl_form_m1

    ka <- exp(lka)

    # --- ODEs (control stream $DES verbatim; 0.947 = M1/tramadol molar
    # ratio applied to the M1 central outflows only) ------------------------
    mw_ratio <- 0.947
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + q / vp * peripheral1 - q / vc * central -
      cl / vc * central - cl_form_m1 / vc * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(central_m1) <- cl_form_m1 / vc * central -
      q_m1 / vc_m1 * central_m1 * mw_ratio + q_m1 / vp_m1 * peripheral1_m1 -
      cl_m1 / vc_m1 * central_m1 * mw_ratio
    d/dt(peripheral1_m1) <- q_m1 / vc_m1 * central_m1 * mw_ratio -
      q_m1 / vp_m1 * peripheral1_m1

    # --- Observations: mg / L x 1000 = ng/mL ($ERROR S1 = VP1, S3 = VM1) ---
    Cc <- 1000 * central / vc
    Cc_m1 <- 1000 * central_m1 / vc_m1
    Cc ~ add(addSd) + prop(propSd) + combined2()
    Cc_m1 ~ add(addSd_m1) + prop(propSd_m1) + combined2()
  })
}
