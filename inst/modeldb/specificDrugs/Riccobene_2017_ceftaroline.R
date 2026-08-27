Riccobene_2017_ceftaroline <- function() {
  description <- "Joint two-compartment ceftaroline fosamil (prodrug) plus two-compartment ceftaroline (active metabolite) population PK model for children and adults, fitted to 6633 plasma concentrations from 305 children aged 1 day to under 18 years pooled with healthy adults, adults with renal impairment and adult patients with ABSSSI or CABP (Riccobene 2017). Ceftaroline fosamil is assumed to be converted completely to ceftaroline, so the whole prodrug elimination clearance CLcf enters the ceftaroline central compartment. Body weight enters allometrically on all clearances (0.75) and volumes (1). Ceftaroline clearance additionally carries a BSA-normalised creatinine clearance term active only below 80 mL/min/1.73 m2, an age term active only above 50 years, a hemodialysis term, a healthy-versus-patient multiplier, and, for children aged 2 years or younger, a Rhodin-style postmenstrual-age renal maturation function that replaces the creatinine clearance term. Ceftaroline central volume carries the same healthy-versus-patient multiplier plus an exponentially decaying postmenstrual-age maturation excess in children aged 2 years or younger. An intramuscular depot (first-order ka, bioavailability fixed to 1) is carried from the upstream adult model; every study contributing to this analysis dosed intravenously. This model supported the FDA and EMA pediatric dose regimens for ceftaroline fosamil."
  reference   <- paste(
    "Riccobene TA, Khariton T, Knebel W, Das S, Li J, Jandourek A, Carrothers TJ, Bradley JS.",
    "Population PK modeling and target attainment simulations to support dosing of ceftaroline",
    "fosamil in pediatric patients with acute bacterial skin and skin structure infections and",
    "community-acquired bacterial pneumonia. J Clin Pharmacol. 2017;57(3):345-355.",
    "doi:10.1002/jcph.809.",
    "Parameter estimates are from Supplemental Table S2 and the covariate model from",
    "Supplemental Equation S1 of that article's supplemental material.",
    "The structural model was carried from the upstream pooled adult and pediatric",
    "ceftaroline fosamil / ceftaroline population PK analysis cited as reference 18",
    "(Riccobene T, Khariton T, Knebel W, O'Neal T, Ghahramani P. 23rd European Congress on",
    "Clinical Microbiology and Infectious Diseases, Berlin 2013, abstract P902);",
    "that abstract is a conference abstract and is not in nlmixr2lib.",
    sep = " "
  )
  vignette    <- "Riccobene_2017_ceftaroline"

  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of every clearance (exponent 0.75) and every volume (exponent 1) of both ceftaroline fosamil and ceftaroline, referenced to 70 kg. Supplemental Table S2 prints the scaling terms as (WT/70)^0.75 and (WT/70)^1 under each structural parameter and Supplemental Equation S1 writes them out explicitly. Weights ranged 1.5 to 100 kg among the 305 children (Results, Study Population).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at study entry.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters ceftaroline clearance twice. (1) Power term (AGE/50)^-0.807, active only above 50 years (Supplemental Equation S1: 'For AGE greater than 50, COV5 = (AGE/50)^-0.807'; COV5 is 1 otherwise). Encoded with max(AGE, 50) so the term is exactly 1 at and below the 50-year pivot. (2) Gate on the two maturation terms: Supplemental Equation S1 activates COV8 (renal maturation on CL) and COV9 (volume maturation on Vcc) only 'For age <= 2'. The gate is a step, so typical clearance is discontinuous at exactly 2 years; that is the published model and is reproduced here (see vignette Errata).",
      source_name        = "AGE"
    ),
    PAGE = list(
      description        = "Postmenstrual age (postnatal age plus gestational age at birth).",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplemental Equation S1 defines PMA_i = AGEW_i + GAGE_i, i.e. postnatal age in weeks plus gestational age at birth in weeks, and uses it in two maturation terms that are active only for AGE <= 2 years: the Rhodin-style renal maturation Hill function FPMA = PMA^1.6/(47.7^1.6 + PMA^1.6) on ceftaroline clearance, and the exponentially decaying volume excess 1 + 1.71*exp(-(PMA - 33)*log(2)/5.51) on ceftaroline central volume. Units here are WEEKS, matching the paper and the Germovsek_2018_meropenem precedent, not the months given as the register default; the reference values 47.7, 33 and 5.51 are all in weeks. Must be supplied for every subject (for an adult it is age in weeks plus 40); it is inert outside the AGE <= 2 gate.",
      source_name        = "PMA"
    ),
    CRCL = list(
      description        = "Creatinine clearance normalised by body surface area.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated by the Cockcroft-Gault equation on total body weight for adults and by the bedside Schwartz equation for children (Methods, Population PK Model). Enters ceftaroline clearance in two mutually exclusive ways. For AGE > 2 years: the power term (nCRCL/80)^0.472, active only below the 80 mL/min/1.73 m^2 reference and only in subjects not on a hemodialysis programme (Supplemental Equation S1 COV3). For AGE <= 2 years the creatinine clearance term is replaced by the postmenstrual-age maturation function (Methods, Population PK Model), and mild renal impairment is represented by scaling that maturation function linearly by CRCL/80 (Methods, Simulations: 'The FPMA was scaled from CrCL scaling of 0.625 (50/80) to 0.988 (79/80) for mild renal impairment'). Both forms are encoded with min(CRCL, 80) so they saturate at the reference and are exactly 1 for normal renal function. Source column nCRCL in the supplement.",
      source_name        = "nCRCL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator (1 = healthy adult subject, 0 = patient with an infection).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (healthy adult subject; the pooled phase 1 healthy-volunteer and renal-impairment cohort)",
      notes              = "Supplemental Table S2 prints the multiplicative terms theta16^PAT = 3.23 on ceftaroline clearance and theta15^PAT = 4.33 on ceftaroline central volume, and Supplemental Equation S1 introduces both with the line 'For patients'. The source column is PAT (1 = patient, 0 = healthy subject) and is re-expressed here as PAT = 1 - DIS_HEALTHY, so the printed typical values CLc = 3.28 L/h and Vcc = 3.6 L are the HEALTHY-subject reference and the multipliers move them to the patient state. That orientation is the one that reproduces the paper's own simulation tables: every age band of Tables 3 to 5, all of which are patients, lands within about 13 percent of Dose/CL under it and about 4-fold off under the opposite orientation (see vignette Errata). Note this is the opposite orientation from the sibling model Riccobene_2016_ceftaroline, whose all-healthy ELF cohort required the multipliers to be active.",
      source_name        = "PAT"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = subject is on an intermittent hemodialysis programme, 0 = not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on hemodialysis)",
      notes              = "Multiplicative power-form effect 0.331^RRT_HEMODIAL_STATUS on ceftaroline clearance between dialysis sessions, and it switches off the creatinine clearance term (Supplemental Equation S1 scopes COV3 to 'non-dialysis patients' and COV4 = 0.331 to 'dialysis patients during non-dialysis period'). Supplemental Table S2 labels the source column ESRD (end-stage renal disease, 1 = yes), so within this model the ESRD flag and the on-a-hemodialysis-programme flag are operationally the same column; the canonical RRT name is used because the term's scope is dialysis, not renal disease generally. No child in this analysis had end-stage renal disease, so the coefficient is informed by the pooled adult renal-impairment data.",
      source_name        = "ESRD"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 while a dialysis session is running, 0 in the interdialytic interval and in non-dialysed subjects).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no dialysis session running)",
      notes              = "Time-varying per-session gate. Supplemental Equation S1 REPLACES ceftaroline clearance entirely during a session (CLc_i = 10.9 L/h, with no covariate or eta terms) rather than adding a dialyser arm to the body clearance, so the model composes the two arms as cl_interdialytic*(1 - RRT_HEMODIAL_ACTIVE) + cl_hemodialysis*RRT_HEMODIAL_ACTIVE. No child in this analysis was dialysed.",
      source_name        = "(dialysis period flag, unnamed in the supplement)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_ceftaroline = list(
      analyte = "ceftaroline", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_ceftaroline = list(
      analyte = "ceftaroline", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 720L,
    n_studies      = 5L,
    age_range      = "Children 1 day to under 18 years (including full-term neonates under 28 days and preterm neonates of 32 to 37 weeks gestational age); adults from the pooled healthy-subject, renal-impairment and phase 2/3 patient studies of the upstream analysis",
    weight_range   = "1.5-100 kg among the 305 children",
    sex_female_pct = 100 * 132 / 305,
    race_ethnicity = "Not reported for the pooled analysis dataset",
    disease_state  = "Children with acute bacterial skin and skin structure infections, community-acquired bacterial pneumonia, complicated community-acquired bacterial pneumonia, or another suspected or confirmed infection requiring antibiotic therapy, plus healthy adult subjects, adults with varying degrees of renal impairment, and adult patients with ABSSSI or CABP. All children except one had normal renal function or mild renal impairment: of those aged 28 days or older, 241 had CrCL of 80 mL/(min 1.73 m^2) or more and 40 had CrCL of 50 to under 80 mL/(min 1.73 m^2).",
    dose_range     = "Children: single doses of 8, 10, 12 or 15 mg/kg (maximum 600 mg) as 1- to 1.5-h infusions, and multiple doses of 8, 10, 12 or 15 mg/kg q8h (maximum 400 or 600 mg) as 1- or 2-h infusions. Adults: 50 to 1000 mg ceftaroline fosamil, including the approved 600 mg q12h regimen.",
    regions        = "Multinational; the five pediatric studies are NCT00633126, NCT01298843, NCT01400867, NCT01530763 and NCT01669980",
    n_observations = "6633 measurable plasma concentrations (1799 ceftaroline fosamil, 4834 ceftaroline), of which 974 (234 ceftaroline fosamil, 740 ceftaroline) came from the 305 children",
    notes          = "n_subjects counts the 525 patients plus 195 healthy subjects or subjects with various degrees of renal impairment reported in Results, Study Population; 305 of the 525 patients were children. n_studies counts the five pediatric studies of Table 1 that this analysis added; the adult data come from a previously pooled dataset of healthy-subject and phase 2/3 patient studies. The children were 173 males and 132 females; the sex split above is therefore the pediatric split, not the whole-dataset split, which is not reported. Concentrations below the limit of quantification (0.01 or 0.05 mg/L) were excluded. Analysis was run in NONMEM 7.3."
  )

  ini({
    # ---------------------------------------------------------------------
    # Ceftaroline fosamil (parent prodrug) -- Supplemental Table S2, Fixed
    # Effects. The table prints back-transformed typical values at the
    # 70 kg allometric reference; Supplemental Equation S1 writes each as
    # exp(theta) * (WT/70)^exponent * exp(eta).
    # ---------------------------------------------------------------------
    lcl <- log(231);  label("Ceftaroline fosamil clearance CLcf at 70 kg (L/h)")                    # Suppl Table S2: CLcf (theta1) = 231 L/h (95% CI 207, 259); Suppl Eq S1 eq for exp(CLcf)_i
    lvc <- log(7.54); label("Ceftaroline fosamil central volume Vccf at 70 kg (L)")                 # Suppl Table S2: Vccf (theta2) = 7.54 L (95% CI 5.09, 11.1)
    lq  <- log(25);   label("Ceftaroline fosamil intercompartmental clearance Q1cf at 70 kg (L/h)") # Suppl Table S2: Q1cf (theta3) = 25 L/h (95% CI 13.7, 46.3)
    lvp <- log(6.75); label("Ceftaroline fosamil peripheral volume Vp1cf at 70 kg (L)")             # Suppl Table S2: Vp1cf (theta4) = 6.75 L (95% CI 5.08, 10.4)

    # Intramuscular route: carried from the upstream adult model. Every
    # study in this analysis dosed intravenously, so neither parameter is
    # informed by these data and Supplemental Table S2 marks both [FIXED].
    lka     <- fixed(log(0.544)); label("Ceftaroline fosamil intramuscular absorption rate ka1cf (1/h)")       # Suppl Table S2: ka1cf (theta5) = 0.544 1/h [FIXED]; Suppl Eq S1: ka1cf_i = 0.544*exp(eta)
    lfdepot <- fixed(log(1));     label("Ceftaroline fosamil intramuscular bioavailability FIMcf (fraction)")  # Suppl Table S2: FIMcf (theta6) = 1 [FIXED]; Suppl Eq S1: F2 = FIMcf

    # ---------------------------------------------------------------------
    # Ceftaroline (active metabolite) -- Supplemental Table S2.
    #
    # lcl_ceftaroline and lvc_ceftaroline carry the HEALTHY-SUBJECT typical
    # values: Supplemental Equation S1 introduces COV7 = 3.23 (on CLc) and
    # COV6 = 4.33 (on Vcc) with the line "For patients", and both default
    # to 1. See covariateData$DIS_HEALTHY and the vignette Errata for the
    # arithmetic that confirms this orientation against Tables 3 to 5.
    # ---------------------------------------------------------------------
    lcl_ceftaroline <- log(3.28); label("Ceftaroline clearance CLc at 70 kg in healthy subjects (L/h)")                     # Suppl Table S2: CLc (theta7) = 3.28 L/h (95% CI 3.14, 3.39)
    lvc_ceftaroline <- log(3.6);  label("Ceftaroline central volume Vcc at 70 kg in healthy subjects (L)")                  # Suppl Table S2: Vcc (theta8) = 3.6 L (95% CI 3.26, 3.91)
    lq_ceftaroline  <- log(8.47); label("Ceftaroline intercompartmental clearance Qc at 70 kg (L/h)")                       # Suppl Table S2: Qc (theta9) = 8.47 L/h (95% CI 7.17, 9.98)
    lvp_ceftaroline <- log(10);   label("Ceftaroline peripheral volume Vp1c at 70 kg (L)")                                  # Suppl Table S2: Vp1c (theta10) = 10 L (95% CI 9.37, 10.6)

    # Dialysis-active ceftaroline clearance. Supplemental Equation S1
    # REPLACES CLc with this value while a session runs, rather than adding
    # a dialyser arm on top of body clearance.
    lcl_hemodialysis_ceftaroline <- log(10.9); label("Ceftaroline clearance during a hemodialysis session, CLdial (L/h)")   # Suppl Table S2: CLdial (theta14) = 10.9 L/h (95% CI 10.2, 11.9)

    # ---------------------------------------------------------------------
    # Covariate effects -- Supplemental Table S2 and Supplemental Equation S1.
    # ---------------------------------------------------------------------
    # Shared allometric exponents, printed under every row of Supplemental
    # Table S2 as (WT/70)^0.75 on clearances and (WT/70)^1 on volumes. No
    # uncertainty is reported for either, so both are fixed.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent shared by all clearances (unitless)")  # Suppl Table S2: (WT/70)^0.75 under CLcf, Q1cf, CLc, Qc
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent shared by all volumes (unitless)")     # Suppl Table S2: (WT/70)^1 under Vccf, Vp1cf, Vcc, Vp1c

    e_crcl_cl_ceftaroline        <- 0.472;  label("Exponent of (CRCL/80) on ceftaroline CL below 80 mL/min/1.73 m^2, age over 2 y")  # Suppl Table S2: (nCRCL/80)^theta11, theta11 = 0.472 (95% CI 0.384, 0.553); Suppl Eq S1 COV3
    e_age_cl_ceftaroline         <- -0.807; label("Exponent of (AGE/50) on ceftaroline CL above 50 years")                           # Suppl Table S2: (AGE/50)^theta13, theta13 = -0.807 (95% CI -0.936, -0.653); Suppl Eq S1 COV5
    e_dis_healthy_cl_ceftaroline <- 3.23;   label("Patient multiplier on ceftaroline CL, applied as base^(1 - DIS_HEALTHY)")         # Suppl Table S2: theta16^PAT = 3.23 (95% CI 3.17, 3.46); Suppl Eq S1 COV7 = 3.23 for patients
    e_dis_healthy_vc_ceftaroline <- 4.33;   label("Patient multiplier on ceftaroline Vc, applied as base^(1 - DIS_HEALTHY)")         # Suppl Table S2: theta15^PAT = 4.33 (95% CI 3.94, 4.9); Suppl Eq S1 COV6 = 4.33 for patients
    e_rrt_hemodial_status_cl_ceftaroline <- 0.331; label("Hemodialysis multiplier on interdialytic ceftaroline CL (power-form base)") # Suppl Table S2: theta12^ESRD = 0.331 (95% CI 0.298, 0.392); Suppl Eq S1 COV4

    # Renal maturation on ceftaroline clearance, active only for AGE <= 2 y.
    # Supplemental Equation S1: FPMA_i = PMA_i^1.6 / (47.7^1.6 + PMA_i^1.6).
    tmat50   <- fixed(47.7); label("Postmenstrual age at 50 percent renal maturation, TM50MAT (weeks)")  # Suppl Table S2: TM50MAT = 47.7 weeks [FIXED]; cited to Rhodin 2009 (paper reference 24)
    hill_mat <- 1.6;         label("Hill coefficient for renal maturation, HillMAT (unitless)")          # Suppl Table S2: theta17 HillMAT = 1.6 (95% CI 1.31, 2.34)

    # Volume maturation on ceftaroline central volume, active only for
    # AGE <= 2 y. Supplemental Equation S1:
    #   COV9 = 1 + Bvol * exp(-(PMA - 33) * log(2) / Tvol)
    # i.e. a fractional volume excess of Bvol at PMA = 33 weeks that decays
    # with half-life Tvol. The 33-week pivot is a hardcoded constant of the
    # published equation, not a reported parameter.
    e_page_vc_ceftaroline <- 1.71; label("Fractional increase in ceftaroline Vc at 33 weeks PMA due to maturation, Bvol (unitless)")  # Suppl Table S2: theta18 Bvol = 1.71 (95% CI 0.67, 2.51)
    thalf_vc_ceftaroline  <- 5.51; label("Maturational half-life for the ceftaroline Vc excess, Tvol (weeks)")                        # Suppl Table S2: theta19 Tvol = 5.51 weeks (95% CI 3.87, 8.07)

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Supplemental Table S2,
    # "Inter-individual Variability". The listed VAR / COV entries form a
    # complete 4x4 lower triangle in the order CLcf, Vccf, CLc, Vcc, plus
    # two diagonal-only entries for ka1cf and Vp1c. Values are log-scale
    # variances: the printed %CV equals 100*sqrt(VAR) for every row (e.g.
    # VAR(CLcf) = 0.435 -> 66%), and every printed correlation equals
    # COV/sqrt(VAR*VAR), so the variances are carried unchanged.
    # ---------------------------------------------------------------------
    # Row-by-row source trace for the block below (comments cannot sit
    # inside the c(...) itself -- rxode2's label parser rejects them):
    #   row 1 (CLcf): VAR(CLcf)      = 0.435   (%CV 66;   95% CI 0.331, 0.894)
    #   row 2 (Vccf): COV(CLcf,Vccf) = -0.132  (r -0.167; 95% CI -0.171, 0.455)
    #                 VAR(Vccf)      = 1.43    (%CV 120;  95% CI 0.757, 2.52)
    #   row 3 (CLc):  COV(CLcf,CLc)  = 0.0459  (r 0.28;   95% CI 0.025, 0.0801)
    #                 COV(Vccf,CLc)  = -0.0684 (r -0.229; 95% CI -0.0822, 0.0158)
    #                 VAR(CLc)       = 0.062   (%CV 24.9; 95% CI 0.0541, 0.0802)
    #   row 4 (Vcc):  COV(CLcf,Vcc)  = 0.14    (r 0.536;  95% CI 0.0876, 0.2)
    #                 COV(Vccf,Vcc)  = 0.0797  (r 0.168;  95% CI 0.0269, 0.275)
    #                 COV(CLc,Vcc)   = 0.0639  (r 0.647;  95% CI 0.0469, 0.089)
    #                 VAR(Vcc)       = 0.158   (%CV 39.7; 95% CI 0.121, 0.21)
    etalcl + etalvc + etalcl_ceftaroline + etalvc_ceftaroline ~ c(
      0.435,
      -0.132, 1.43,
      0.0459, -0.0684, 0.062,
      0.14, 0.0797, 0.0639, 0.158
    )
    etalka             ~ fixed(0.114)   # VAR(ka1cf) = 0.114 (%CV 33.8); 95% CI printed as (0.114, 0.114), a zero-width interval
    etalvp_ceftaroline ~ 0.0148         # VAR(Vp1c)  = 0.0148 (%CV 12.2; 95% CI 0.0117, 0.0259)

    # ---------------------------------------------------------------------
    # Residual variability -- Supplemental Table S2, "Residual Variability".
    # The table reports VARIANCES; the parenthetical %CV / SD columns are
    # their square roots (e.g. propCF 0.137 -> %CV 37, addCF 0.00445 ->
    # SD 0.0667), so the SDs below are sqrt(variance). The reported
    # COV(propCF, propC) = 0.0305 (r = 0.407), a correlation between the two
    # analytes' proportional residuals, has no nlmixr2 representation and is
    # not carried (see vignette Errata).
    # ---------------------------------------------------------------------
    propSd             <- 0.3701; label("Ceftaroline fosamil proportional residual SD (fraction)") # Suppl Table S2: propCF = 0.137 variance (%CV 37) -> sqrt(0.137) = 0.3701
    addSd              <- 0.0667; label("Ceftaroline fosamil additive residual SD (mg/L)")         # Suppl Table S2: addCF  = 0.00445 variance (SD 0.0667)
    propSd_ceftaroline <- 0.2027; label("Ceftaroline proportional residual SD (fraction)")         # Suppl Table S2: propC  = 0.0411 variance (%CV 20.3) -> sqrt(0.0411) = 0.2027
    addSd_ceftaroline  <- 0.0523; label("Ceftaroline additive residual SD (mg/L)")                 # Suppl Table S2: addC   = 0.00274 variance (SD 0.0523)
  })

  model({
    # ----- Reference covariate values (Supplemental Equation S1)
    ref_wt   <- 70   # kg, allometric reference
    ref_age  <- 50   # years, pivot of the COV5 age term
    ref_crcl <- 80   # mL/min/1.73 m^2, pivot of the COV3 renal term
    ref_page <- 33   # weeks PMA, pivot of the COV9 volume-maturation term

    allom_cl <- (WT / ref_wt)^e_wt_cl_q
    allom_v  <- (WT / ref_wt)^e_wt_vc_vp

    # ----- Maturation gate. Supplemental Equation S1 activates COV8 and
    # COV9 only "For age <= 2" (years) and the Methods state that for these
    # children the postmenstrual-age maturation term "replaced the term
    # representing the effect of body surface area-normalized creatinine
    # clearance". The gate is therefore a step at exactly 2 years.
    is_mat <- (AGE <= 2)

    # ----- Renal function. Saturates at the 80 mL/min/1.73 m^2 reference,
    # so it is exactly 1 for normal renal function.
    renal_ratio <- min(CRCL, ref_crcl) / ref_crcl

    # COV3: power term on ceftaroline CL, only above 2 years of age and
    # only in subjects not on a hemodialysis programme. Zeroing the
    # exponent reproduces both gates without branching.
    crcl_cl <- renal_ratio^(e_crcl_cl_ceftaroline * (1 - is_mat) * (1 - RRT_HEMODIAL_STATUS))

    # COV5: age term on ceftaroline CL, only above 50 years. The max()
    # holds the ratio at 1 up to the pivot.
    age_cl <- (max(AGE, ref_age) / ref_age)^e_age_cl_ceftaroline

    # COV8: Rhodin-style renal maturation on ceftaroline CL for age <= 2 y,
    # scaled linearly by CRCL/80 to carry mild renal impairment in this age
    # band (Methods, Simulations). Inert (exactly 1) above 2 years.
    fmat   <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)
    mat_cl <- 1 + is_mat * (fmat * renal_ratio - 1)

    # COV9: exponentially decaying volume excess on ceftaroline Vc for
    # age <= 2 y. Inert (exactly 1) above 2 years.
    mat_vc <- 1 + is_mat * e_page_vc_ceftaroline *
      exp(-(PAGE - ref_page) * log(2) / thalf_vc_ceftaroline)

    # ----- Individual parameters: ceftaroline fosamil (parent prodrug)
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q  <- exp(lq)  * allom_cl
    vp <- exp(lvp) * allom_v
    ka <- exp(lka + etalka)

    # ----- Individual parameters: ceftaroline (active metabolite).
    # Clearance switches wholesale to the dialysis-session value while
    # RRT_HEMODIAL_ACTIVE is 1 (Supplemental Equation S1: "For dialysis
    # patients during dialysis, CLc_i = 10.9"). The patient multipliers are
    # applied as base^(1 - DIS_HEALTHY) because the source column is
    # PAT = 1 - DIS_HEALTHY.
    cl_ceftaroline <-
      exp(lcl_ceftaroline + etalcl_ceftaroline) * allom_cl * crcl_cl * age_cl * mat_cl *
        e_rrt_hemodial_status_cl_ceftaroline^RRT_HEMODIAL_STATUS *
        e_dis_healthy_cl_ceftaroline^(1 - DIS_HEALTHY) *
        (1 - RRT_HEMODIAL_ACTIVE) +
      exp(lcl_hemodialysis_ceftaroline) * RRT_HEMODIAL_ACTIVE
    vc_ceftaroline <- exp(lvc_ceftaroline + etalvc_ceftaroline) * allom_v * mat_vc *
      e_dis_healthy_vc_ceftaroline^(1 - DIS_HEALTHY)
    q_ceftaroline  <- exp(lq_ceftaroline)                       * allom_cl
    vp_ceftaroline <- exp(lvp_ceftaroline + etalvp_ceftaroline) * allom_v

    # ----- Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_ceftaroline <- cl_ceftaroline / vc_ceftaroline
    k12_ceftaroline <- q_ceftaroline  / vc_ceftaroline
    k21_ceftaroline <- q_ceftaroline  / vp_ceftaroline

    # ----- ODE system.
    # Ceftaroline fosamil is assumed to be converted completely to
    # ceftaroline, so its whole elimination clearance CLcf becomes the
    # formation flux into the ceftaroline central compartment: Supplemental
    # Equation S1 carries no formation fraction and no molar conversion, and
    # both analytes are reported in mg/L. Mass is transferred 1:1 in mg,
    # matching the sibling model Riccobene_2016_ceftaroline.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_ceftaroline) <-
      kel * central -
      kel_ceftaroline * central_ceftaroline -
      k12_ceftaroline * central_ceftaroline + k21_ceftaroline * peripheral1_ceftaroline
    d/dt(peripheral1_ceftaroline) <-
      k12_ceftaroline * central_ceftaroline - k21_ceftaroline * peripheral1_ceftaroline

    # ----- Intramuscular bioavailability (Supplemental Equation S1: F2 = FIMcf)
    f(depot) <- exp(lfdepot)

    # ----- Observations and residual error (Supplemental Equation S1:
    # C_ij = Chat_ij + Chat_ij*eps_p + eps_a for each analyte)
    Cc             <- central             / vc
    Cc_ceftaroline <- central_ceftaroline / vc_ceftaroline

    Cc             ~ add(addSd)             + prop(propSd)
    Cc_ceftaroline ~ add(addSd_ceftaroline) + prop(propSd_ceftaroline)
  })
}
