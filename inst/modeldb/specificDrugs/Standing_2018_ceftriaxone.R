Standing_2018_ceftriaxone <- function() {
  description <- "Three-compartment population PK model for intravenous ceftriaxone in hospitalised infants and young children (2-45 months) with complicated severe acute malnutrition in Kenya (Standing 2018). Disposition is parameterised on UNBOUND ceftriaxone: the central state holds total drug, the unbound concentration is derived from the total central concentration by instantaneous saturable (Michaelis-Menten) protein binding, and only unbound drug is cleared or distributes to the two peripheral (unbound) compartments. Allometric weight scaling (0.75 on clearances, 1 on volumes), a fixed Rhodin postmenstrual-age renal maturation and a power effect of serum creatinine relative to its age-expected value on clearance, a power effect of serum albumin on the binding capacity Bmax, and an oedema-corrected weight (any oedema lowers the weight used for volumes by 26%)."
  reference <- "Standing JF, Ongas MO, Ogwang C, Kagwanja N, Murunga S, Mwaringa S, Ali R, Mturi N, Timbwa M, Manyasi C, Mwalekwa L, Bandika VL, Ogutu B, Waichungo J, Kipper K, Berkley JA. Dosing of Ceftriaxone and Metronidazole for Children With Severe Acute Malnutrition. Clin Pharmacol Ther. 2018;104(6):1165-1174. doi:10.1002/cpt.1078. PMID 29574688. Structure, covariate equations and final estimates from the Supplementary Material section 2.1 (ceftriaxone NONMEM control stream)."
  vignette <- "Standing_2018_ceftriaxone_metronidazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg: exponent 0.75 on CL, Q2 and Q3 (actual weight) and 1 on V1, V2 and V3 (oedema-corrected weight, see DIS_EDEMA). Cohort median 5.88 kg (range 2.53-10.9; Table 1).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the fixed Rhodin 2009 sigmoidal renal maturation of unbound clearance, PMA^3.4 / (PMA^3.4 + 47.7^3.4) (stream $PK; Methods 'Biological prior information'). The maturation fraction is not normalised to a 70-kg adult, so a mature subject has a factor just below 1. The youngest study patient had a PMA of 46 weeks (Discussion). Supply postnatal age plus gestational age at birth (40 weeks for term).",
      source_name = "PMA"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to compute the age-expected mean serum creatinine (Ceriotti) that normalises CREAT: MSCR = -2.37330 - 12.91367 * log(AGE) + 23.93581 * AGE^0.5 umol/L (stream '$PK', AGEY = AGED / 365.25). Must be > 0. Cohort median 15 months (range 2-45; Table 1).",
      source_name = "AGED"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on unbound clearance, (CREAT / MSCR)^-0.261, where MSCR is the age-expected mean serum creatinine computed in model() from AGE (Results; stream THETA(10)). Cohort median 21 umol/L (range 4-178; Table 1).",
      source_name = "SECR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the protein-binding capacity Bmax, (ALB / 33.75)^0.555, centred on the ceftriaxone-analysis median 33.75 g/L (Table 2 footnote; stream THETA(9)). Lower albumin gives a lower Bmax and a higher unbound fraction. Range 7-57 g/L (Table 1).",
      source_name = "ALB"
    ),
    DIS_EDEMA = list(
      description = "Clinical (nutritional) oedema present, any grade",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no oedema)",
      notes = "1 = oedema grade 1, 2 or 3 (kwashiorkor-type bilateral pitting oedema), 0 = grade 0. Derived from the stream's ordinal OEDM column as 'IF (OEDM > 0) WTC = WT * (1 + THETA(11))'. The oedema-corrected weight WTC = WT * (1 - 0.256) replaces WT in the allometric scaling of V1, V2 and V3 only. 35 of 80 patients had oedema (grades 1/2/3 = 12/14/9; Table 1). The metronidazole model in the same paper uses a grade >= 2 cut instead (DIS_EDEMA_GE2).",
      source_name = "OEDM"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "TOTAL (bound + unbound) ceftriaxone amount; Cc = central / vc is the measured total plasma concentration (stream COMP(TOTAL), 'CTOT = A(1)/V1')."
    ),
    peripheral1 = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Unbound ceftriaxone amount in the first peripheral compartment (stream COMP(PREFREE1))."
    ),
    peripheral2 = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Unbound ceftriaxone amount in the second peripheral compartment (stream COMP(PERFREE2))."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80,
    n_studies = 1,
    age_range = "2-45 months",
    age_median = "15 months",
    weight_range = "2.53-10.9 kg",
    weight_median = "5.88 kg",
    sex_female_pct = 56.25,
    race_ethnicity = "Kenyan (East African) children",
    disease_state = "Complicated severe acute malnutrition (SAM) requiring hospital admission: kwashiorkor, mid-upper-arm circumference < 11.5 cm (< 11 cm at 1-5 months) or weight-for-height Z-score < -3. Oedema grade 0/1/2/3 = 45/12/14/9.",
    dose_range = "50 mg/kg (100 mg/kg if meningitis suspected) intravenous bolus every 24 h; sampling up to 72 h after the first dose",
    regions = "Kenya (three centres; NCT02746276, the PK stage of the FLACSAM trial)",
    renal_function = "Serum creatinine median 21 umol/L (range 4-178)",
    albumin = "Serum albumin median 33.75 g/L (range 7-57)",
    notes = "Standing 2018 Table 1 (ceftriaxone analysis column). 81 patients enrolled; one excluded from the ceftriaxone analysis for ceftriaxone in the pre-dose sample. Three post-dose samples per patient (randomised early 5 min / 30 min / 1 h, middle 2 / 4 / 8 h, late 24 / 48 / 72 h); total ceftriaxone in all and unbound ceftriaxone (ultrafiltration) in 190 of the 234 post-dose samples. NONMEM 7.3 FOCE-I."
  )

  ini({
    # Typical values for a 70-kg individual with serum creatinine equal to its
    # age-expected value (Table 2; Supplementary Material section 2.1 $THETA).
    # Table 2's parameter-name column is offset by one row from its unlabelled
    # third row down; the stream $THETA order (and its comments) is used.
    lcl <- log(4.61627); label("Unbound ceftriaxone clearance CL, 70 kg (L/h)") # stream THETA(1) TVCL = 4.61627; Table 2 'CL' 4.62 (95% CI 4.06-5.24)
    lvc <- log(14.4644); label("Total ceftriaxone central volume V1, 70 kg (L)") # stream THETA(2) TVV1 = 14.4644; Table 2 'V1' 14.46 (9.47-17.76)
    lq <- log(0.933031); label("Unbound intercompartmental clearance Q2, 70 kg (L/h)") # stream THETA(3) TVQ2 = 0.933031; Table 2 unlabelled row 0.93 (0.16-5.61)
    lvp <- log(6.05411); label("Unbound first peripheral volume V2, 70 kg (L)") # stream THETA(4) TVV2 = 6.05411; Table 2 value 6.05 (printed against 'Q2')
    lq2 <- log(29.9815); label("Unbound intercompartmental clearance Q3, 70 kg (L/h)") # stream THETA(5) TVQ3 = 29.9815; Table 2 value 29.98 (printed against 'V2')
    lvp2 <- log(19.2606); label("Unbound second peripheral volume V3, 70 kg (L)") # stream THETA(6) TVV3 = 19.2606; Table 2 value 19.26 (printed against 'Q3')

    # Saturable protein binding (Methods equations; stream $DES)
    lbmax_pb <- log(69.7743); label("Maximum protein-binding capacity Bmax at albumin 33.75 g/L (mg/L)") # stream THETA(7) TVBM = 69.7743; Table 2 value 69.77 (printed against 'V3')
    lkd_pb <- log(22.8853); label("Binding dissociation constant Kd (mg/L)") # stream THETA(8) TVKD = 22.8853; Table 2 value 22.89 (printed against 'BM'); Discussion 'Kd of 22.9 mg/L'

    # Covariate effects
    e_alb_bmax_pb <- 0.555443; label("Power exponent of albumin (ALB/33.75) on Bmax (unitless)") # stream THETA(9) 'ALB power cov' = 0.555443; Table 2 value 0.56 (printed against 'KD')
    e_creat_cl <- -0.261234; label("Power exponent of serum creatinine relative to age-expected value on CL (unitless)") # stream THETA(10) 'SECR power cov' = -0.261234; Table 2 -0.26
    e_dis_edema_wt <- -0.255959; label("Fractional change in the weight used for volume scaling with any oedema (unitless)") # stream THETA(11) 'WT correction for any oedema' = -0.255959; Results '26%'

    # Fixed a priori size and maturation priors (Methods; stream $PK literals)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL, Q2 and Q3 (unitless)") # stream '(WT/70)**0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of oedema-corrected WT on V1, V2 and V3 (unitless)") # stream '(WTC/70)'
    pma50_cl <- fixed(47.7); label("Postmenstrual age at half-maximal renal maturation of CL (weeks)") # stream 'PMA**3.4/(PMA**3.4 + 47.7**3.4)'; Rhodin 2009 (paper ref 10)
    h_pma_cl <- fixed(3.4); label("Hill coefficient of renal maturation of CL (unitless)") # stream 'PMA**3.4'

    # IIV (stream $OMEGA; variances on the log scale). CVs sqrt(omega):
    # CL 42%, V1 43%, Bmax 42% (Table 2 IIV column, Bmax's printed against 'V3').
    etalcl + etalvc ~ c(0.174903, 0.0807096, 0.18196) # stream $OMEGA BLOCK(2) 0.174903 / 0.0807096 / 0.18196
    etalbmax_pb ~ 0.172225 # stream $OMEGA 0.172225 (ETA(3) on BM)

    # Residual error (stream $SIGMA variances; Y = IPRED*(1+EPS) + EPS)
    propSd <- 0.19039; label("Proportional residual SD, total ceftriaxone (fraction)") # stream $SIGMA(3) 0.0362471 -> sqrt = 0.1904; Table 2 '19%'
    addSd <- 0.22361; label("Additive residual SD, total ceftriaxone (mg/L)") # stream $SIGMA(4) 0.05 (not FIX) -> sqrt = 0.2236
    propSd_Cu <- 0.26663; label("Proportional residual SD, unbound ceftriaxone (fraction)") # stream $SIGMA(1) 0.0710923 -> sqrt = 0.2666; Table 2 '27%'; additive $SIGMA(2) is 0 FIX
  })

  model({
    # Oedema-corrected weight for the volume terms (stream $PK)
    wtc <- WT * (1 + e_dis_edema_wt * DIS_EDEMA)

    # Age-expected mean serum creatinine, umol/L (Ceriotti; stream MSCR)
    creat_exp <- -2.37330 - 12.91367 * log(AGE) + 23.93581 * AGE^0.5

    # Fixed Rhodin renal maturation
    fmat_cl <- PAGE^h_pma_cl / (PAGE^h_pma_cl + pma50_cl^h_pma_cl)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat_cl * (CREAT / creat_exp)^e_creat_cl
    vc <- exp(lvc + etalvc) * (wtc / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (wtc / 70)^e_wt_vc
    q2 <- exp(lq2) * (WT / 70)^e_wt_cl
    vp2 <- exp(lvp2) * (wtc / 70)^e_wt_vc

    bmax_pb <- exp(lbmax_pb + etalbmax_pb) * (ALB / 33.75)^e_alb_bmax_pb
    kd_pb <- exp(lkd_pb)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Total central concentration and the unbound concentration from
    # instantaneous Michaelis-Menten binding (Methods equations; stream $DES)
    Cc <- central / vc
    Cu <- 0.5 * ((Cc - bmax_pb - kd_pb) + sqrt((Cc - bmax_pb - kd_pb)^2 + 4 * kd_pb * Cc))

    # Only unbound drug leaves the central compartment. The stream writes
    # K*A(1)*FU with FU = CFREE/CTOT; A(1)*FU = Cu*vc avoids 0/0 before a dose.
    d/dt(central) <- -(kel + k12 + k13) * Cu * vc + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * Cu * vc - k21 * peripheral1
    d/dt(peripheral2) <- k13 * Cu * vc - k31 * peripheral2

    Cc ~ add(addSd) + prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
