Aksenov_2018_uricAcid <- function() {
  description <- "Semi-mechanistic exposure-response model of uric acid disposition in adults with drug effects of xanthine-oxidase inhibitors (allopurinol via oxypurinol, febuxostat) and a uricosuric URAT1 inhibitor (lesinurad)"
  reference <- "Aksenov S, Peck CC, Eriksson UG, Stanski DR. Individualized treatment strategies for hyperuricemia informed by a semi-mechanistic exposure-response model of uric acid dynamics. Physiol Rep. 2018 Mar;6(5):e13614. doi:10.14814/phy2.13614"
  vignette <- "Aksenov_2018_uricAcid"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "mg/dL"
  )

  covariateData <- list(
    CRCL = list(
      description = "Renal function as glomerular filtration rate (Cockcroft-Gault creatinine clearance, actual body weight, NOT BSA-normalized) -- the source paper uses the symbol GFR but operationally computes a creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per simulation. Source paper Methods (p. 5) explicitly states 'GFR was approximated with creatinine clearance calculated using the Cockcroft-Gault formula and actual body weight without adjustment or normalization.' Reference value 130 mL/min for a typical 75-kg adult (Aksenov 2018 Appendix 2 Table 2). Range across study populations 22-164 mL/min (Methods p. 5).",
      source_name = "GFR"
    ),
    CP_OXY_NGML = list(
      description = "Plasma oxypurinol concentration (active metabolite of allopurinol; xanthine oxidase inhibitor)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input. Set to 0 to disable allopurinol/oxypurinol effect. Mean daily concentration on 300 mg/day allopurinol is approximately 10000 ng/mL (Aksenov 2018 Eq. 13 narrative; Schumacher 2008). The PK model for oxypurinol is not parameterized in this paper -- users must supply CP_OXY_NGML from an external PK source (e.g., Wright et al. 2013) or use a steady-state value.",
      source_name = "[P]_PIN (oxypurinol)"
    ),
    CP_FBX_NGML = list(
      description = "Plasma febuxostat concentration (xanthine oxidase inhibitor)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input. Set to 0 to disable febuxostat effect. The PK model for febuxostat is not parameterized in this paper -- users must supply CP_FBX_NGML from an external PK source (e.g., Bhattaram & Gobburu 2017 regulatory review) or use a steady-state value.",
      source_name = "[P]_PIN (febuxostat)"
    ),
    CP_LSN_NGML = list(
      description = "Plasma lesinurad concentration (uricosuric URAT1 reabsorption inhibitor)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input. Set to 0 to disable lesinurad effect. The PK model for lesinurad is not parameterized in this paper -- users must supply CP_LSN_NGML from an external PK source (e.g., Fleischmann et al. 2014, Shen et al. 2015) or use a steady-state value.",
      source_name = "[P]_RIN (lesinurad)"
    )
  )

  population <- list(
    n_subjects = 278,
    n_studies = 9,
    age_range = "adults",
    weight_range = "reference 75 kg (Appendix 2 Table 2)",
    sex_female_pct = NA_real_,
    disease_state = "Healthy adults and gout patients with hyperuricemia",
    dose_range = "Lesinurad 50-1600 mg single or 100-600 mg/day; allopurinol 300 mg/day; febuxostat 40 or 80 mg/day",
    regions = "United States and Japan (Study 125)",
    notes = "Estimation dataset D1: 278 subjects from 9 Phase I lesinurad studies (Studies 101, 102, 103, 105, 109, 110, 111, 117, 125) with 4455 serum-UA samples and 3058 urine-UA samples in 18 treatment groups (Aksenov 2018 Table 4). Mean baseline values per treatment group given in Aksenov 2018 Table 5; GFR range 108-151 mL/min in healthy cohorts and 22-164 mL/min when including renal-impaired qualification cohorts (Studies 104 and 120). Qualification datasets D3 (renal impairment, Studies 104 and 120, n = 39) and D4 (Phase III gout patients, Studies 301-304, n = 647 in active arms) used to qualify the final model under reduced GFR and in hyperuricemic patients (Aksenov 2018 Methods p. 4-6, Appendix 3). Reference subject characteristics for typical-value parameter estimates: 75 kg, BMI 22, male, GFR 130 mL/min (Appendix 2 Table 2)."
  )

  ini({
    # ----- UA disposition (mechanistic) -----
    cl_i <- 0.27;   label("Intestinal clearance of UA (L/h)")           # Aksenov 2018 Table 1
    v_ua <- 19;     label("Volume of distribution of UA (L)")           # Aksenov 2018 Table 1

    # Typical-value baseline production rate and fractional excretion (Appendix 2 Table 2;
    # used to derive the steady-state baseline serum UA via Eq. 5).
    kp0 <- 50.5;    label("Baseline UA production rate, typical 75-kg adult (mg/h)")  # Aksenov 2018 Appendix 2 Table 2
    fe0 <- 0.07;    label("Baseline fractional excretion of UA in urine (dimensionless)")  # Aksenov 2018 Appendix 2 Table 2

    # ----- Drug-effect parameters (Eq. 9 and 10) -----
    # XOI -- allopurinol (via oxypurinol).
    rmax_oxy <- 0.84;   label("Max fractional decrease in UA production by oxypurinol (dimensionless)")  # Aksenov 2018 Table 1
    p50_oxy  <- 14000;  label("Oxypurinol concentration at half-maximal production inhibition (ng/mL)")  # Aksenov 2018 Table 1

    # XOI -- febuxostat.  Default p50_fbx is the hyperuricemia/gout estimate (Table 1);
    # for healthy/normouricemic subjects override to p50_fbx = 87 ng/mL.
    rmax_fbx <- 1;      label("Max fractional decrease in UA production by febuxostat (fixed in source)")  # Aksenov 2018 Table 1; fixed at 1 per Bhattaram & Gobburu 2017
    p50_fbx  <- 120;    label("Febuxostat concentration at half-maximal production inhibition, hyperuricemic subjects (ng/mL)")  # Aksenov 2018 Table 1

    # Uricosuric -- lesinurad.  Default p50_lsn is the hyperuricemia/gout estimate (Table 1);
    # for healthy/normouricemic subjects override to p50_lsn = 11000 ng/mL.
    fmax_lsn <- 0.56;   label("Max increase in fractional excretion by lesinurad (dimensionless)")  # Aksenov 2018 Table 1; fixed during estimation
    p50_lsn  <- 23000;  label("Lesinurad concentration at half-maximal increase in fractional excretion, hyperuricemic subjects (ng/mL)")  # Aksenov 2018 Table 1

    # ----- Residual error (Eq. 4: y = f + sqrt(a^2 + (b*f)^2) * eps) -----
    addSd_sUA  <- 0.45;  label("Additive residual SD on serum UA (mg/dL)")            # Aksenov 2018 Table 1
    propSd_sUA <- 0.15;  label("Proportional residual SD on serum UA (fraction)")     # Aksenov 2018 Table 1
    addSd_uUA  <- 50;    label("Additive residual SD on cumulative urinary UA (mg)")  # Aksenov 2018 Table 1
    propSd_uUA <- 0.29;  label("Proportional residual SD on cumulative urinary UA (fraction)")  # Aksenov 2018 Table 1
  })

  model({
    # GFR is in mL/min in the source paper (Cockcroft-Gault); convert to L/h for ODE balance.
    gfr_lh <- CRCL * 60 / 1000

    # Drug effects (Aksenov 2018 Eq. 9 and Eq. 10).
    # Production inhibition by XOIs.  Each XOI follows the paper's Eq. 9 form; when both
    # XOI concentrations are non-zero (off-label combination, not in the paper's data),
    # they combine via a Bliss-independence form so the surviving fraction stays in [0, 1].
    # When only one XOI concentration is non-zero (the typical clinical case) the form
    # collapses exactly to Eq. 9 with that single drug.
    pin_oxy <- rmax_oxy * CP_OXY_NGML / (CP_OXY_NGML + p50_oxy)
    pin_fbx <- rmax_fbx * CP_FBX_NGML / (CP_FBX_NGML + p50_fbx)
    kp <- kp0 * (1 - pin_oxy) * (1 - pin_fbx)

    # Increase in fractional excretion by lesinurad (Eq. 10).
    fe <- fe0 + fmax_lsn * CP_LSN_NGML / (CP_LSN_NGML + p50_lsn)

    # Serum concentration of UA (Eq. 3).
    sua_mgL <- serum / v_ua
    sUA <- sua_mgL / 10

    # ODE system (Eq. 1 and Eq. 2).
    d/dt(serum) <- kp - cl_i * sua_mgL - gfr_lh * fe * sua_mgL
    d/dt(urine) <- gfr_lh * fe * sua_mgL

    # Drug-free steady-state baseline serum UA from Eq. 5 used as the initial condition.
    bl_sua_mgL <- kp0 / (cl_i + gfr_lh * fe0)
    serum(0) <- bl_sua_mgL * v_ua
    urine(0) <- 0

    # Output: cumulative urinary UA in mg.
    uUA <- urine

    # Residual error -- additive + proportional combination matches Eq. 4
    # variance form a^2 + (b*f)^2 (additive and proportional eps independent).
    sUA ~ add(addSd_sUA) + prop(propSd_sUA)
    uUA ~ add(addSd_uUA) + prop(propSd_uUA)
  })
}
