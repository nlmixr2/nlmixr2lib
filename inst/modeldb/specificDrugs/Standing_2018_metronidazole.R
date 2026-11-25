Standing_2018_metronidazole <- function() {
  description <- "Parent-metabolite population PK model for oral metronidazole and its hydroxy metabolite (hydroxymetronidazole) in hospitalised infants and young children (2-45 months) with complicated severe acute malnutrition in Kenya (Standing 2018). One-compartment disposition for each analyte with first-order absorption after a lag time. Hydroxymetronidazole is formed both during first pass and from circulating metronidazole through a well-stirred liver model: an estimated intrinsic clearance against a fixed allometric hepatic plasma flow (scaled by haematocrit) gives the extraction ratio, which sets the first-pass fraction converted and the systemic formation clearance. Allometric weight scaling (0.75 on clearances and hepatic flow, 1 on volumes, -0.25 on ka), a fixed CYP2A6 age maturation on intrinsic clearance, a fixed Rhodin postmenstrual-age maturation on metabolite clearance, and an oedema-corrected weight (oedema grade >= 2 lowers the weight used for volumes by 17%). Parameters are apparent (no intravenous data)."
  reference <- "Standing JF, Ongas MO, Ogwang C, Kagwanja N, Murunga S, Mwaringa S, Ali R, Mturi N, Timbwa M, Manyasi C, Mwalekwa L, Bandika VL, Ogutu B, Waichungo J, Kipper K, Berkley JA. Dosing of Ceftriaxone and Metronidazole for Children With Severe Acute Malnutrition. Clin Pharmacol Ther. 2018;104(6):1165-1174. doi:10.1002/cpt.1078. PMID 29574688. Structure, covariate equations and final estimates from the Supplementary Material section 2.2 (metronidazole NONMEM control stream) and section 3 (identical estimates in the initial-estimate perturbation table)."
  vignette <- "Standing_2018_ceftriaxone_metronidazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg: exponent 0.75 on intrinsic clearance, metronidazole clearance by other routes, hydroxymetronidazole clearance and hepatic blood flow; -0.25 on ka (actual weight); 1 on both volumes (oedema-corrected weight, see DIS_EDEMA_GE2). Cohort median 5.88 kg (range 2.53-10.9; Table 1).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the fixed Rhodin 2009 sigmoidal renal maturation of hydroxymetronidazole clearance, PMA^3.4 / (PMA^3.4 + 47.7^3.4) (stream $PK). Supply postnatal age plus gestational age at birth (40 weeks for term).",
      source_name = "PMA"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the fixed CYP2A6 maturation of intrinsic clearance (Upreti; paper ref 11), (1.1 - 0.15) * AGE^0.8 / (0.04^0.8 + AGE^0.8) + 0.15 (stream $PK, AGEY = AGED / 365.25). Not normalised to an adult, so a mature subject has a factor close to 1.1. Cohort median 17 months (range 2-45; Table 1).",
      source_name = "AGED"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Converts hepatic blood flow to hepatic plasma flow, QHP = QH * (1 - 0.91 * HCT/100) (Methods, after Price et al.; stream 'QHP = QH * (1 - HCT * 0.91)'). The stream column HCT is a volume FRACTION; the canonical column is percent, so model() divides by 100. Haematocrit is not tabulated in the paper (Table 1 reports haemoglobin, median 9 g/dL, range 4.4-12.8).",
      source_name = "HCT"
    ),
    DIS_EDEMA_GE2 = list(
      description = "Nutritional oedema of grade 2 or 3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oedema grade 0 or 1)",
      notes = "1 = oedema grade >= 2, 0 = grade 0 or 1. Derived from the stream's ordinal OEDM column as 'IF (OEDM > 1) WTC = WT * (1 + THETA(8))'. The oedema-corrected weight WTC = WT * (1 - 0.169) replaces WT in the allometric scaling of both volumes only (Results: 'an edema score of 2 or more leading to a decrease in apparent weight of 17%'). 21 of 69 patients in the metronidazole analysis (grades 2/3 = 13/8; Table 1). DIS_EDEMA_GE2 = 1 implies the any-grade indicator DIS_EDEMA = 1 used by the ceftriaxone model of the same paper.",
      source_name = "OEDM"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "metronidazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral suspension dose. Replaces the stream's two depots DEPOTPAR (F1 = 1 - ER) and DEPOTMET (F2 = ER), which share ka and the lag time; a single depot whose absorbed flux is split (1 - eh) / eh is mathematically identical and needs one dose record instead of two."
    ),
    central = list(
      analyte = "metronidazole",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Metronidazole amount (stream COMP(PARENT)); Cc = central / vc."
    ),
    central_2oh = list(
      analyte = "hydroxymetronidazole",
      units = "mg (metronidazole mass, no molecular-weight conversion)",
      specimen = "plasma",
      verified = TRUE,
      notes = "Hydroxymetronidazole amount (stream COMP(METAB)); Cc_2oh = central_2oh / vc_2oh. The stream moves mass from parent to metabolite without a molecular-weight ratio (metronidazole 171.15, hydroxymetronidazole 187.15 g/mol), so the 9% mass difference is absorbed into the apparent metabolite volume."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 69,
    n_studies = 1,
    age_range = "2-45 months",
    age_median = "17 months",
    weight_range = "2.53-10.9 kg",
    weight_median = "5.88 kg",
    sex_female_pct = 56.5,
    race_ethnicity = "Kenyan (East African) children",
    disease_state = "Complicated severe acute malnutrition (SAM) requiring hospital admission: kwashiorkor, mid-upper-arm circumference < 11.5 cm (< 11 cm at 1-5 months) or weight-for-height Z-score < -3. Oedema grade 0/1/2/3 = 37/11/13/8.",
    dose_range = "7.5 mg/kg oral suspension every 8 h; sampling up to 72 h after the first dose",
    regions = "Kenya (three centres; NCT02746276, the PK stage of the FLACSAM trial)",
    albumin = "Serum albumin median 30.7 g/L (range 7-57)",
    notes = "Standing 2018 Table 1 (metronidazole analysis column). 81 patients enrolled; 12 excluded from the metronidazole analysis for metronidazole and/or hydroxymetronidazole in the pre-dose sample. Three post-dose samples per patient (randomised early / middle / late times); metronidazole and hydroxymetronidazole assayed by HPLC-ESI-MS/MS (LLOQ 0.05 and 0.02 ug/mL). Samples after vomiting (1.6% of 572 doses) were retained. NONMEM 7.3 FOCE-I. Median extraction ratio 0.015 and median fraction metabolised 0.636 (Supplementary Material section 3)."
  )

  ini({
    # Typical values for a 70-kg individual (Table 2; Supplementary Material
    # section 2.2 $THETA, identical to the 'original' row of the section 3.2
    # perturbation table).
    lka <- log(0.393084); label("Absorption rate constant ka, 70 kg (1/h)") # stream THETA(1) TVKA = 0.393084; Table 2 'Ka' 0.39 (0.29-0.6)
    lclint <- log(0.743542); label("Apparent intrinsic clearance of metronidazole to hydroxymetronidazole CLI, 70 kg (L/h)") # stream THETA(2) TVCLI = 0.743542; Table 2 'CLIMET' 0.74 (0.51-1.08)
    lcl <- log(0.505611); label("Apparent metronidazole clearance by other routes, 70 kg (L/h)") # stream THETA(3) TVCL = 0.505611; Table 2 'CLMET' 0.51 (0.23-0.76)
    lvc <- log(56.3524); label("Apparent metronidazole volume, 70 kg (L)") # stream THETA(4) TVV1 = 56.3524; Table 2 'VMET' 56.35 (47.73-65.48)
    lcl_2oh <- log(1.45238); label("Apparent hydroxymetronidazole clearance, 70 kg (L/h)") # stream THETA(5) TVCLM = 1.45238; Table 2 'CLOH' 1.45 (1.08-1.95)
    lvc_2oh <- log(16.8914); label("Apparent hydroxymetronidazole volume, 70 kg (L)") # stream THETA(6) TVV2 = 16.8914; Table 2 'VOH' 16.89 (11.67-24.77)
    ltlag <- log(0.171638); label("Absorption lag time, both analytes (h)") # stream THETA(7) AL1 = 0.171638 (bounds 0.1-0.75); Table 2 'Alag' 0.17 (0.11-0.26)

    e_dis_edema_ge2_wt <- -0.16871; label("Fractional change in the weight used for volume scaling with oedema grade >= 2 (unitless)") # stream THETA(8) 'Oedema cov' = -0.16871; Table 2 'thetaODEM_MET' -0.17

    # Fixed a priori size and maturation priors (Methods; stream $PK literals)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on all clearances and hepatic blood flow (unitless)") # stream '(WT/70)**0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent of oedema-corrected WT on both volumes (unitless)") # stream '(WTC/70)'
    e_wt_ka <- fixed(-0.25); label("Allometric exponent of WT on ka (unitless)") # stream '(WT/70)**(-0.25)'
    fmax_age_clint <- fixed(1.1); label("Maximum CYP2A6 maturation fraction on CLI (unitless)") # stream '(1.1 - 0.15) * ... + 0.15'; Upreti (paper ref 11)
    fmin_age_clint <- fixed(0.15); label("CYP2A6 maturation fraction at birth on CLI (unitless)") # stream '+ 0.15'
    age50_clint <- fixed(0.04); label("Age at half-maximal CYP2A6 maturation of CLI (years)") # stream '0.04**0.8'
    h_age_clint <- fixed(0.8); label("Hill coefficient of CYP2A6 maturation of CLI (unitless)") # stream 'AGEY**0.8'
    pma50_cl_2oh <- fixed(47.7); label("Postmenstrual age at half-maximal renal maturation of hydroxymetronidazole clearance (weeks)") # stream 'PMA**3.4/(PMA**3.4 + 47.7**3.4)'; Rhodin 2009 (paper ref 10)
    h_pma_cl_2oh <- fixed(3.4); label("Hill coefficient of renal maturation of hydroxymetronidazole clearance (unitless)") # stream 'PMA**3.4'

    # IIV (stream $OMEGA, diagonal; CL, CLM and V2 are '0 FIX' and carry no eta).
    # CVs sqrt(omega): ka 79%, CLI 61%, V 38% (Table 2 IIV column).
    etalka ~ 0.626046 # stream $OMEGA(1) 0.626046
    etalclint ~ 0.376697 # stream $OMEGA(2) 0.376697
    etalvc ~ 0.141027 # stream $OMEGA(4) 0.141027

    # Residual error (stream $SIGMA variances; Y = IPRED*(1+EPS) + EPS)
    propSd <- 0.32335; label("Proportional residual SD, metronidazole (fraction)") # stream $SIGMA(1) 0.104557 -> sqrt = 0.3234; Table 2 '32%'; additive $SIGMA(2) is 0 FIX
    propSd_2oh <- 0.28805; label("Proportional residual SD, hydroxymetronidazole (fraction)") # stream $SIGMA(3) 0.0829706 -> sqrt = 0.2880; Table 2 '29%'
    addSd_2oh <- 0.12221; label("Additive residual SD, hydroxymetronidazole (mg/L)") # stream $SIGMA(4) 0.0149359 -> sqrt = 0.1222; Table 2 prints the variance as 0.01 (0.002-0.04)
  })

  model({
    # Oedema-corrected weight for the volume terms (stream $PK)
    wtc <- WT * (1 + e_dis_edema_ge2_wt * DIS_EDEMA_GE2)

    # Fixed maturation functions (stream $PK)
    fmat_clint <- (fmax_age_clint - fmin_age_clint) * AGE^h_age_clint /
      (age50_clint^h_age_clint + AGE^h_age_clint) + fmin_age_clint
    fmat_cl_2oh <- PAGE^h_pma_cl_2oh / (PAGE^h_pma_cl_2oh + pma50_cl_2oh^h_pma_cl_2oh)

    ka <- exp(lka + etalka) * (WT / 70)^e_wt_ka
    clint <- exp(lclint + etalclint) * (WT / 70)^e_wt_cl * fmat_clint
    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (wtc / 70)^e_wt_vc
    cl_2oh <- exp(lcl_2oh) * (WT / 70)^e_wt_cl * fmat_cl_2oh
    vc_2oh <- exp(lvc_2oh) * (wtc / 70)^e_wt_vc
    tlag <- exp(ltlag)

    # Well-stirred liver: 76 L/h adult hepatic blood flow scaled allometrically,
    # converted to plasma flow with 0.91 * haematocrit (stream QH, QHP, ER, CLTM)
    qh <- 76 * (WT / 70)^e_wt_cl
    fq <- qh * (1 - 0.91 * HCT / 100)
    eh <- clint / (fq + clint)
    clf_2oh <- fq * eh

    kel <- cl / vc
    kf_2oh <- clf_2oh / vc
    kel_2oh <- cl_2oh / vc_2oh

    # Absorbed drug: fraction (1 - eh) reaches the circulation as metronidazole
    # (stream F1, DEPOTPAR) and fraction eh is converted on first pass (F2, DEPOTMET)
    d/dt(depot) <- -ka * depot
    lag(depot) <- tlag
    d/dt(central) <- (1 - eh) * ka * depot - (kel + kf_2oh) * central
    d/dt(central_2oh) <- eh * ka * depot + kf_2oh * central - kel_2oh * central_2oh

    Cc <- central / vc
    Cc_2oh <- central_2oh / vc_2oh

    Cc ~ prop(propSd)
    Cc_2oh ~ add(addSd_2oh) + prop(propSd_2oh)
  })
}
