Braune_2018_meropenem <- function() {
  description <- "Two-compartment IV population PK model for meropenem in 19 septic critically ill adults with acute kidney injury receiving sustained low-efficiency dialysis (Braune 2018). Pmetrics NPAG non-parametric fit parameterised by a central volume and the rate constants Kcp / Kpc. Total clearance is the additive sum of three arms: a non-renal arm (2.6 L/h), a native renal arm proportional to 24-hour residual diuresis (1.5 L/h per 100 mL/24h), and a SLED arm (7.9 L/h) switched on only while a dialysis session is running. Blood/dialysate flow, ultrafiltration rate and body weight were screened but not retained."
  reference <- "Braune S, Konig C, Roberts JA, Nierhaus A, Steinmetz O, Baehr M, Kluge S, Langebrake C. Pharmacokinetics of meropenem in septic patients on sustained low-efficiency dialysis: a population pharmacokinetic study. Crit Care. 2018;22(1):25. doi:10.1186/s13054-018-1940-1"
  vignette <- "Braune_2018_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    RRT_CRRT_ACTIVE = list(
      description = "Sustained low-efficiency dialysis (SLED) session currently running",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SLED off; interdialytic interval)",
      notes = paste(
        "Braune 2018 Results, 'Pharmacokinetic model building': 'The value for term SLED is 1 when SLED is on,",
        "whereas it is 0 when SLED is off.' Per-time-point (within-subject time-varying) gate, NOT a subject-level",
        "treatment-status flag. SLED is a prolonged-intermittent modality and belongs to the RRT_CRRT_* family,",
        "which the register entry names explicitly ('sustained low-efficiency dialysis (SLED)'); it is not the",
        "intermittent-hemodialysis gate RRT_HEMODIAL_ACTIVE. SLED was delivered with the Genius batch system and a",
        "Fresenius FX 60 filter (1.4 m2); median session duration 315 min [IQR 275-354, range 80-470] (Table 1).",
        "The paper's own Monte Carlo simulations used a 5-hour session beginning 17 h after the first dose",
        "(Methods, 'Probability of target attainment'). Source-paper alias: the 'SLED' term of the TVCL equation."
      ),
      source_name = "SLED"
    ),
    URINE_VOL_24H = list(
      description = "Residual diuresis: total urine volume produced over 24 hours",
      units = "mL/24h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Braune 2018 Table 1 'Residual diuresis [mL/d]': median 0 [IQR 0-80], range 0-360 across the 19-subject",
        "cohort. Enters the printed equation TVCLNS = CLD + (CLN * RD / 100) as a purely PROPORTIONAL scaling of",
        "the native renal clearance arm normalised at 100 mL/24h, so cl_renal is exactly zero in an anuric subject",
        "and the non-SLED clearance reduces to the non-renal arm CLD. This is the multiplicative branch the",
        "register entry sanctions. Note the contrast with Ulldemolins_2015_meropenem.R, which uses the same",
        "RD/100 normalisation but as an ADDITIVE slope on a non-zero baseline clearance",
        "(CL = 3.68 + 0.22 * RD/100); Braune 2018 cites that paper (ref 50) for the underlying relationship but",
        "parameterises it as a separate clearance arm with its own NPAG distribution. Simulation values used by",
        "the paper: 0, 100 and 300 mL/24h (Figs. 3-5, Tables 3-5). Source-paper alias: 'RD' / 'residual diuresis'."
      ),
      source_name = "RD"
    )
  )

  covariatesDataExcluded <- list(
    BFR = list(
      description = "Blood flow rate through the SLED circuit",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened on CL; not retained. Braune 2018 Results: 'No other covariates could be identified as",
        "significant, e.g. blood/dialysate and ultrafiltration flow rate.' The Genius batch system runs blood and",
        "dialysate at a single common flow, which Table 1 reports as one combined 'Blood/dialysate flow [mL/min]'",
        "row: median 250 [IQR 208-278], range 170-350. BFR and DFR therefore carry the same value in this cohort."
      )
    ),
    DFR = list(
      description = "Dialysate flow rate through the SLED circuit",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened on CL; not retained. Braune 2018 Results. Reported jointly with BFR in the single Table 1",
        "'Blood/dialysate flow [mL/min]' row (median 250 [IQR 208-278], range 170-350) because the Genius batch",
        "system drives both at a common rate."
      )
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Ultrafiltration rate during SLED",
      units = "mL/h",
      type = "continuous",
      notes = paste(
        "Screened on CL; not retained. Braune 2018 Results: 'No other covariates could be identified as",
        "significant, e.g. blood/dialysate and ultrafiltration flow rate.' Table 1: median 500 mL/h",
        "[IQR 400-597], range 50-1000."
      )
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened as a covariate on meropenem PK; not retained. Braune 2018 Methods, 'Population pharmacokinetic",
        "modeling': 'Demographic and clinical characteristics, which were considered biologically plausible for",
        "affecting meropenem PK, such as residual diuresis, blood/dialysate flow and bodyweight, were tested for",
        "inclusion as covariates.' Table 1: median 81 kg [IQR 76-90], range 70-183. No allometric or linear weight",
        "term appears in the final model, so all volumes and clearances below are absolute, not per-kg."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "37-78 years",
    age_median = "66 years",
    weight_range = "70-183 kg",
    weight_median = "81 kg",
    sex_female_pct = 26.3,
    race_ethnicity = "Not reported (single-center German ICU cohort)",
    disease_state = paste(
      "Critically ill septic adults with acute kidney injury requiring renal replacement therapy. Median SOFA",
      "score 11 [IQR 9-13, range 5-16] on the first day of sampling and 10 [IQR 9-13] on ICU admission; median",
      "SAPS II 46 [IQR 34-50] on admission. Mechanical ventilation in 58% and vasopressors in 89% on sampling",
      "day 1. Median C-reactive protein 101 mg/L [IQR 51-162], procalcitonin 1.48 ug/L [IQR 0.72-3.02], serum",
      "albumin 15.5 g/L (SD 3.7). ICU mortality 47%; median ICU length of stay 36 days [IQR 23-79]."
    ),
    dose_range = paste(
      "Meropenem 0.5 g, 1 g or 2 g IV over 30 minutes 8-hourly, at the discretion of the treating physician",
      "(Methods, 'Dosing, administration and data collection')."
    ),
    regions = "Germany (Department of Intensive Care Medicine, University Medical Center Hamburg-Eppendorf). Enrolled July 2013 - November 2014. ClinicalTrials.gov NCT02287493.",
    renal_function = paste(
      "All patients had acute kidney injury and received SLED. Median residual diuresis 0 mL/24h [IQR 0-80,",
      "range 0-360], i.e. the median patient was anuric. Median serum creatinine 2.6 mg/dL [IQR 1.7-4.0,",
      "range 0.8-9.1] on sampling day 1, but the paper states this value is confounded by prior SLED or CRRT",
      "sessions and that endogenous renal function could therefore not be estimated by standard approaches",
      "(Limitations)."
    ),
    n_concentrations = 308L,
    notes = paste(
      "Demographics from Braune 2018 Table 1. Serial sampling on three consecutive days of SLED: a trough 1 h",
      "before infusion, then 10 min, 1 h, 2 h and 4 h after the start of SLED, and at the end of the session;",
      "SLED began no later than 3 h after the meropenem infusion. Post-SLED samples were NOT collected, so any",
      "rebound is unobserved (Limitations). Total meropenem measured by HPLC-UV, LLOQ 1 mg/L, intra-day",
      "precision CV 9.6% / 3.9% / 2.2% at 10 / 20 / 80 mg/L. Median observed trough 28.9 mg/L [IQR 21.6-36.9,",
      "range 10.2-95.8] across the pooled dose levels. Fitted with the Nonparametric Adaptive Grid (NPAG)",
      "algorithm in Pmetrics; Monte Carlo simulations (n = 1000) supported the PTA and FTA analyses.",
      "Single-vs-multicenter design: single-center prospective observational cohort."
    )
  )

  ini({
    # Structural parameters: Braune 2018 Table 2, 'Mean' column of the
    # Pmetrics NPAG non-parametric population distribution. Table 2 also
    # reports a 'Median' column; the mean is used as the typical value here
    # -- it is the column the Discussion quotes throughout ('CLSLED = 7.9
    # L/h', 'CLNS is composed of CLN (1.5 L/h) and CLD (2.6 L/h)', 'The
    # central volume of distribution (Vc) in our patients of 8.1 L') -- and
    # the median is noted per line. This follows the sibling Pmetrics NPAG
    # extractions Tsai_2023_ceftriaxone.R and Duke_2024_cefazolin.R.
    #
    # The clearance decomposition is the paper's own printed equation
    # (Results, 'Pharmacokinetic model building'; equation image
    # 13054_2018_1940_Article_Equa.gif in the publisher's file set, which
    # carries the multiplication signs the PDF text layer drops):
    #     TVCL   = CLSLED * (SLED) + CLNS
    #     TVCLNS = CLD + (CLN * RD / 100)
    # so the SLED arm is ADDITIVE on top of the non-SLED arms, not a
    # replacement for them. Total CL while a session runs and RD = 0 is
    # therefore 7.9 + 2.6 = 10.5 L/h. See the vignette 'Assumptions and
    # deviations' for the one Discussion sentence that reads loosely
    # against this and why the printed equation governs.
    lcl_nonren <- log(2.6)
    label("Non-SLED, non-renal clearance arm CLD (L/h)")
    # Braune 2018 Table 2: CLD mean 2.6, SD 1.2, CV 44.9%, median 2.3 L/h

    lcl_renal <- log(1.5)
    label("Native renal clearance arm CLN at URINE_VOL_24H = 100 mL/24h (L/h)")
    # Braune 2018 Table 2: CLN mean 1.5, SD 2.1, CV 134.7%, median 0.7 L/h

    lcl_crrt <- log(7.9)
    label("SLED-session clearance arm CLSLED (L/h)")
    # Braune 2018 Table 2: CLSLED mean 7.9, SD 4.2, CV 53.6%, median 6.8 L/h

    lvc <- log(8.1)
    label("Central volume of distribution Vc (L)")
    # Braune 2018 Table 2: Vc mean 8.1, SD 7.1, CV 87.9%, median 4.9 L

    lk12 <- log(10.3)
    label("Central-to-peripheral rate constant Kcp (1/h)")
    # Braune 2018 Table 2: KCP mean 10.3, SD 8.8, CV 85.4%, median 7.9 1/h

    lk21 <- log(1.8)
    label("Peripheral-to-central rate constant Kpc (1/h)")
    # Braune 2018 Table 2: KPC mean 1.8, SD 1.9, CV 104.4%, median 1.2 1/h

    # Interindividual variability. Pmetrics NPAG estimates a discrete
    # non-parametric distribution rather than a parametric omega matrix;
    # Table 2 summarises that distribution by its mean, SD and CV%. The
    # CV% is carried here into a log-normal random effect using the
    # standard omega^2 = log(CV^2 + 1) identity. This is a parametric
    # APPROXIMATION of a non-parametric distribution (see vignette
    # 'Assumptions and deviations'); it is required to reproduce the
    # paper's own Monte Carlo PTA / FTA simulations, which sample the
    # population distribution. Table 2 reports no correlations between
    # parameters, so the etas are left uncorrelated.
    #   CLD    :  44.9% CV -> omega^2 = log(0.449^2 + 1) = 0.183655
    #   CLN    : 134.7% CV -> omega^2 = log(1.347^2 + 1) = 1.034752
    #   CLSLED :  53.6% CV -> omega^2 = log(0.536^2 + 1) = 0.252544
    #   Vc     :  87.9% CV -> omega^2 = log(0.879^2 + 1) = 0.572471
    #   KCP    :  85.4% CV -> omega^2 = log(0.854^2 + 1) = 0.547726
    #   KPC    : 104.4% CV -> omega^2 = log(1.044^2 + 1) = 0.737133
    etalcl_nonren ~ 0.183655 # Braune 2018 Table 2 (CLD, CV 44.9%)
    etalcl_renal ~ 1.034752 # Braune 2018 Table 2 (CLN, CV 134.7%)
    etalcl_crrt ~ 0.252544 # Braune 2018 Table 2 (CLSLED, CV 53.6%)
    etalvc ~ 0.572471 # Braune 2018 Table 2 (Vc, CV 87.9%)
    etalk12 ~ 0.547726 # Braune 2018 Table 2 (KCP, CV 85.4%)
    etalk21 ~ 0.737133 # Braune 2018 Table 2 (KPC, CV 104.4%)

    # Residual error. Braune 2018 Methods states only that 'Additive
    # (lambda) and exponential (gamma) error models were both tested for
    # inclusion' and Results that 'A two-compartment linear model using an
    # additive error adequately described the serum concentrations'. The
    # fitted lambda is not reported anywhere in the paper, and the
    # publisher's supplementary file set for PMC5791175 contains only the
    # five figure images plus the equation image -- no Pmetrics model file
    # and no error-polynomial table. The additive term is therefore carried
    # as fixed(0) rather than substituted from the HPLC-UV assay precision
    # (which is a different quantity: an assay CV, not a fitted residual
    # SD). Same convention as the sibling Pmetrics extraction
    # Setiawan_2023_sulbactam.R. See vignette 'Assumptions and deviations'.
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; 0 -- magnitude not reported in the source)")
  })
  model({
    # Normalising residual diuresis for the native-renal clearance arm.
    # Hard-coded as '/100' in the paper's printed TVCLNS equation; 100
    # mL/24h is also the anuria cutoff the URINE_VOL_24H register entry
    # records.
    rd_ref <- 100

    # Individual clearance arms (L/h). Each carries its own log-normal
    # random effect because Table 2 reports a separate NPAG marginal
    # distribution for each.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (URINE_VOL_24H / rd_ref)
    cl_crrt <- exp(lcl_crrt + etalcl_crrt)

    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # Peripheral volume and inter-compartmental clearance, derived from the
    # rate constants the paper reports. Braune 2018 parameterises
    # distribution by Kcp / Kpc only, so neither of these is a published
    # value; both follow algebraically from q = Kcp * Vc and
    # Vp = Vc * Kcp / Kpc (typical values 83.4 L/h and 46.4 L).
    #
    # They are NOT decorative. Declaring `cl` and `vc` without a `q` / `vp`
    # pair makes rxode2 match this file against its ONE-compartment
    # analytic linear-compartment kernel: it then discards the explicit
    # d/dt() right-hand sides, drops `peripheral1` from the solve entirely,
    # and returns a mono-exponential profile decaying at exactly cl / vc.
    # That failure is silent -- `ui$linCmt` reads empty, the ODE block
    # parses and prints correctly, `checkModelConventions()` is clean, and
    # the AUC-recovery identity cl * AUCinf == Dose still holds because it
    # is blind to the number of compartments. Supplying q and vp restores a
    # consistent two-compartment specification; the vignette's
    # 'Two-compartment disposition is real' gate pins the solved profile to
    # the closed-form biexponential so a regression cannot pass silently.
    # This is a sharper form of the hazard documented in the
    # RRT_CRRT_ACTIVE / RRT_HEMODIAL_ACTIVE entries of
    # covariate-columns.md, which describes the gate going inert but not
    # the whole peripheral compartment disappearing.
    vp <- vc * k12 / k21
    q <- k12 * vc

    # Total clearance, written as a SINGLE assignment to `cl` including the
    # SLED gate. Routing the gated sum through a separate `cl_total`
    # variable can leave the dialysis arm silently inert, because rxode2
    # sometimes solves the linear system analytically from variables named
    # `cl` / `vc` and discards the explicit d/dt() -- see the encoding rule
    # in the RRT_CRRT_ACTIVE / RRT_HEMODIAL_ACTIVE entries of
    # covariate-columns.md and the mechanical gate
    # tests/testthat/test-modeldb-active-gate.R.
    cl <- cl_nonren + cl_renal + RRT_CRRT_ACTIVE * cl_crrt
    kel <- cl / vc

    # Two-compartment disposition parameterised by Vc and the rate
    # constants Kcp / Kpc, as Pmetrics / NPAG reports it (Table 2); the
    # peripheral volume is not separately identifiable and is not reported.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central/vc has units mg/L, matching the
    # paper's observed serum concentrations (HPLC-UV, LLOQ 1 mg/L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
