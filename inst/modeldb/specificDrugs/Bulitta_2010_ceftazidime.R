Bulitta_2010_ceftazidime <- function() {
  description <- "Three-compartment population PK model for ceftazidime after 5-min IV infusion in cystic fibrosis patients and healthy volunteers (Bulitta 2010), with allometric fat-free-mass scaling and a cystic-fibrosis-vs-healthy disease-group factor on total clearance."
  reference   <- "Bulitta JB, Landersdorfer CB, Huttner SJ, Drusano GL, Kinzig M, Holzgrabe U, Stephan U, Sorgel F. Population pharmacokinetic comparison and pharmacodynamic breakpoints of ceftazidime in cystic fibrosis patients and healthy volunteers. Antimicrob Agents Chemother. 2010;54(3):1275-1282. doi:10.1128/AAC.00936-09"
  vignette    <- "Bulitta_2010_ceftazidime"
  units       <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on clearances (exponent 0.75) and linear scaling on volumes (exponent 1.0) with reference 53 kg. Bulitta 2010 Methods 'Body size model' fixes the exponents; Table 3 footnote states that all CL and V estimates are group estimates for subjects of standard size, fat-free mass = 53 kg. FFM was derived per subject from total body weight, height, and sex via the Janmahasatian et al. (2005) formula (Methods, Subjects).",
      source_name        = "FFM"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cystic fibrosis patient)",
      notes              = "Multiplicative effects on CL and on V1 / V2 / V3 relative to the CF reference. Bulitta 2010 Table 3 parameterises the disease contrast through scale factors FCYFCL = 1.17 and FCYFVSS = 1.01 with the healthy volunteer cohort serving as the structural reference (e.g., CL_CF = 7.82 L/h = 6.68 L/h x 1.17 at FFM = 53 kg). The model file uses the canonical DIS_HEALTHY orientation (1 = healthy, 0 = patient) so the typical-value parameters represent the CF reference (DIS_HEALTHY = 0) and the e_healthy_* effects shift them toward the HV estimates when DIS_HEALTHY = 1. The reorientation preserves the published parameter values exactly: e_healthy_cl = log(1 / FCYFCL) and e_healthy_vc / vp / vp2 = log(1 / FCYFVSS). The intercompartmental clearances Q and Q2 are not stratified by disease group (Table 3 reports identical values 27.9 and 2.57 L/h for both cohorts). Source NONMEM dataset column name is not stated in the paper; the source orientation likely was a CF-indicator that was re-expressed as DIS_HEALTHY = 1 - CF in the package."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    n_cf           = 8L,
    n_hv           = 7L,
    age_range      = "10-45 years (CF 10-45; HV 19-33; Table 1)",
    age_median     = "CF 20 years; HV 22 years",
    weight_range   = "14.2-73.5 kg (CF 14.2-73.5; HV 56-71)",
    weight_median  = "CF 37.9 kg; HV 67 kg",
    ffm_range      = "13.9-60.5 kg (CF 13.9-57.7; HV 46.4-60.5; Janmahasatian formula)",
    ffm_median     = "CF 35.9 kg; HV 54.0 kg",
    height_range   = "105-191 cm (CF 105-176; HV 164-191)",
    height_median  = "CF 164 cm; HV 175 cm",
    bmi_range      = "12.5-23.7 kg/m^2 (CF 12.5-23.7; HV 19.5-23.7)",
    sex_female_pct = 47.0,
    race_ethnicity = "All Caucasian (Bulitta 2010 Methods, Subjects)",
    disease_state  = "8 cystic fibrosis patients with normal renal function and 7 healthy volunteers; CF patients smaller and leaner than HV (Tables 1-2). Study performed in 1983.",
    dose_range     = "Single 2 g IV over 5 min (one CF patient received 1 g, one 1.5 g, one 3 g per physician judgment)",
    regions        = "Germany (single-centre)",
    samples        = "21 plasma samples per subject between 0 and 12 h after end of infusion (0, 5, 10, 15, 20, 30, 45, 60, 90 min and 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12 h after end of infusion)",
    notes          = "Allometric scaling by fat-free mass reduced the unexplained between-subject variance by 32% on CL and by 18-26% on the peripheral volumes relative to linear scaling by total weight (Bulitta 2010 Table 5). The paper estimated the same model with three independent estimation algorithms (NONMEM, S-ADAPT, NPAG); this model file reproduces the NONMEM column of Table 3."
  )

  ini({
    # Structural parameters at the cystic-fibrosis reference (DIS_HEALTHY = 0)
    # and FFM = 53 kg, taken from Bulitta 2010 Table 3 (NONMEM column,
    # 'CF patients' rows). Cf. footnote b: 'All clearance and volume
    # estimates are group estimates of the respective PK parameter for
    # subjects of standard size (fat-free mass: 53 kg).'
    lcl  <- log(7.82);  label("Total clearance at CF reference, FFM = 53 kg (L/h)")                                    # Bulitta 2010 Table 3 NONMEM CF
    lvc  <- log(5.73);  label("Central volume at CF reference, FFM = 53 kg (V1, L)")                                   # Bulitta 2010 Table 3 NONMEM CF
    lvp  <- log(3.92);  label("Shallow peripheral volume at CF reference, FFM = 53 kg (V2, L)")                        # Bulitta 2010 Table 3 NONMEM CF
    lvp2 <- log(3.16);  label("Deep peripheral volume at CF reference, FFM = 53 kg (V3, L)")                           # Bulitta 2010 Table 3 NONMEM CF
    lq   <- log(27.9);  label("Intercompartmental clearance to shallow peripheral at FFM = 53 kg (Q, L/h)")            # Bulitta 2010 Table 3 NONMEM (same for both groups)
    lq2  <- log(2.57);  label("Intercompartmental clearance to deep peripheral at FFM = 53 kg (Q2, L/h)")              # Bulitta 2010 Table 3 NONMEM (same for both groups)

    # Allometric body-size exponents -- fixed per Bulitta 2010 Methods,
    # 'Body size model': clearance exponent fixed to 0.75, volume exponent
    # fixed to 1.0. The 'allometric scaling by FFM' final model is the
    # column reproduced here (Table 3 column NONMEM).
    e_ffm_cl <- fixed(0.75); label("Allometric exponent of FFM on all clearances CL / Q / Q2 (unitless)")              # Bulitta 2010 Methods, Body size model
    e_ffm_vc <- fixed(1.0);  label("Linear exponent of FFM on all volumes V1 / V2 / V3 (unitless)")                    # Bulitta 2010 Methods, Body size model

    # Healthy-volunteer effect (Bulitta 2010 Table 3 scale factors).
    # FCYFCL = CL_CF / CL_HV = 7.82 / 6.68 = 1.17 (Table 4 row 5).
    # FCYFVSS = V_CF / V_HV = 1.01, applied identically to V1, V2, V3
    # (Table 4 row 5; the three volumes in Table 3 share a single
    # estimated scale factor). On the DIS_HEALTHY orientation the
    # log-effects are the inverse: e_healthy_<param> = log(1 / FCYF).
    e_healthy_cl  <- log(1 / 1.17); label("Log effect of DIS_HEALTHY on total CL (= log(1 / FCYFCL))")                 # Bulitta 2010 Table 4 FFM allometric row, FCYFCL = 1.17
    e_healthy_vc  <- log(1 / 1.01); label("Log effect of DIS_HEALTHY on V1 (= log(1 / FCYFVSS))")                      # Bulitta 2010 Table 4 FFM allometric row, FCYFVSS = 1.01
    e_healthy_vp  <- log(1 / 1.01); label("Log effect of DIS_HEALTHY on V2 (= log(1 / FCYFVSS))")                      # Bulitta 2010 Table 4 FFM allometric row, FCYFVSS = 1.01
    e_healthy_vp2 <- log(1 / 1.01); label("Log effect of DIS_HEALTHY on V3 (= log(1 / FCYFVSS))")                      # Bulitta 2010 Table 4 FFM allometric row, FCYFVSS = 1.01

    # Between-subject variability (Bulitta 2010 Table 3 NONMEM column).
    # The published apparent CV%s describe an exponential-IIV model; the
    # log-scale variance is omega^2 = log(CV^2 + 1). Volumes are correlated
    # per Table 3 footnote e: r(V1,V2) = -0.67, r(V1,V3) = +0.56,
    # r(V2,V3) = -0.59. CL is independent. Table 3 footnote c states that
    # one joint variance was estimated across CF and HV cohorts per
    # parameter (small sample size).
    etalcl                          ~ 0.07555                                                                          # Bulitta 2010 Table 3: CV 28% on CL -> log(0.28^2 + 1)
    etalvc + etalvp + etalvp2       ~ c(0.18438,
                                        -0.07084, 0.06062,
                                         0.06834, -0.04128, 0.08075)                                                   # Bulitta 2010 Table 3 and footnote e: CV 45/25/29% on V1/V2/V3 with paired correlations

    # Residual error -- combined additive + proportional (Bulitta 2010
    # Methods, 'Between-subject variability and observation model';
    # Table 3 footer reports CV C = 0.122 and SD C = 0.059 mg/L).
    propSd <- 0.122; label("Proportional residual SD on Cc (fraction)")                                                # Bulitta 2010 Table 3 NONMEM CV C
    addSd  <- 0.059; label("Additive residual SD on Cc (mg/L)")                                                        # Bulitta 2010 Table 3 NONMEM SD C
  })

  model({
    # Allometric / linear size factors centred on the published
    # FFM = 53 kg reference (Bulitta 2010 Table 3 footnote b).
    size_cl <- (FFM / 53)^e_ffm_cl
    size_v  <- (FFM / 53)^e_ffm_vc

    # Individual parameters. The DIS_HEALTHY indicator shifts the CF
    # reference toward the healthy-volunteer estimates by log(1 / FCYF);
    # at DIS_HEALTHY = 0 the typical values equal the Table 3 'CF
    # patients' column, and at DIS_HEALTHY = 1 they recover the
    # 'healthy volunteers' column. Q and Q2 are not stratified by
    # disease group (Bulitta 2010 Table 3: identical values for both
    # cohorts).
    cl  <- exp(lcl  + etalcl  + e_healthy_cl  * DIS_HEALTHY) * size_cl
    vc  <- exp(lvc  + etalvc  + e_healthy_vc  * DIS_HEALTHY) * size_v
    vp  <- exp(lvp  + etalvp  + e_healthy_vp  * DIS_HEALTHY) * size_v
    vp2 <- exp(lvp2 + etalvp2 + e_healthy_vp2 * DIS_HEALTHY) * size_v
    q   <- exp(lq)                                            * size_cl
    q2  <- exp(lq2)                                           * size_cl

    # Three-compartment IV disposition micro-constants. The shallow
    # peripheral (peripheral1) carries the higher intercompartmental
    # clearance Q; the deep peripheral (peripheral2) carries the slower
    # Q2.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system. Doses are administered into the central compartment
    # as a short (~5 min) IV infusion via the data event table; the
    # 'duration of zero order input was fixed to 5 min and not
    # estimated' (Bulitta 2010 Table 3 footer).
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Ceftazidime plasma concentration in mg/L (dose in mg / vc in L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
