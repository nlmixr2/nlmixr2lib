Abouelhassan_2025_sulbactam <- function() {
  description <- paste(
    "Two-compartment population PK simulation model for sulbactam (given as",
    "sulbactam-durlobactam) in critically ill adults with ventilator-associated",
    "bacterial pneumonia, with and without continuous renal replacement therapy",
    "(CRRT) (Abouelhassan 2025). Without CRRT, total clearance is a non-renal",
    "arm plus a creatinine-clearance-driven renal arm, scaled by body weight and",
    "an infection-type shift (supplement Table S11). On CRRT the renal arm is",
    "replaced by the transmembrane clearance measured in an ex vivo CRRT",
    "circuit, a linear function of the effluent flow rate",
    "(CLTM = 0.114 + 0.973 x effluent, L/h), added to the non-renal clearance",
    "(supplement Table S12). Unbound concentration uses a fixed 38% protein",
    "binding."
  )
  reference <- paste(
    "Abouelhassan Y, Shen Y, Chen A, Ye X, Nicolau DP, Kuti JL. Ex vivo",
    "assessment of sulbactam-durlobactam clearance during continuous renal",
    "replacement therapy to guide dosing recommendations. Antimicrob Agents",
    "Chemother. 2025;69(1):e01674-23. doi:10.1128/aac.01674-23.",
    "CRRT clearance regression from the Results ('CRRT clearance studies');",
    "PK parameters and dispersion from supplement Tables S11 and S12.",
    sep = " "
  )
  vignette <- "Abouelhassan_2025_sulbactam_durlobactam_crrt"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function with reference 75 kg on central volume in both arms and",
        "on total clearance in the non-CRRT arm only (supplement Table S11 rows",
        "'Total Body CL' and 'Vc'; the Table S12 CRRT clearance carries no",
        "weight term). Simulated as 76.4 +/- 4.7 kg (mean +/- SD)."
      ),
      source_name = "Wt"
    ),
    CRCL = list(
      description = "Creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the renal clearance arm as (CRCL / 90)^0.932 in the non-CRRT arm",
        "(supplement Table S11). Ignored when RRT_CRRT_STATUS = 1, because the",
        "CRRT simulations replace the renal arm by the transmembrane clearance",
        "and do not account for residual renal function (Discussion). Units are",
        "printed as ml/min without BSA normalization; the simulated VABP",
        "population had a mean of 80 mL/min (range 30-129 mL/min)."
      ),
      source_name = "CLCR"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia indicator (1 = VABP)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference infection type of the underlying population PK model; not identified in the source)",
      notes = paste(
        "Proportional shifts on total clearance (CLINF = -0.294, non-CRRT arm",
        "only) and on central volume (VcINF = 1.38, both arms), supplement",
        "Table S11 footnotes a and b. The paper simulated only VABP (or",
        "HABP/VABP) patients, so set DIS_VABP = 1 to reproduce it. The Vc shift",
        "is described as applying to HABP/VABP patients, but the sulbactam",
        "clearance shift for HABP is not printed, so HABP is not encoded."
      ),
      source_name = "CLINF / VcINF (infection type)"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator (1 = on CRRT)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT; supplement Table S11 model)",
      notes = paste(
        "Selects between the two clearance models of the source. 0: Table S11",
        "VABP model, total CL = (4.79 + 0.648 x 13.6 x (CRCL/90)^0.932) x",
        "(WT/75)^1.03 x (1 + CLINF). 1: Table S12 CRRT model, total CL =",
        "CL_NR + CL_TM, with CL_TM from the effluent-flow regression and",
        "between-subject variability on CL_NR (30% CV) and CL_TM (12% CV)."
      ),
      source_name = "CRRT"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Total CRRT effluent flow rate",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper reports effluent flow in L/h; divided by 1000 inside model().",
        "Effluent is the total ultrafiltration rate in CVVH and dialysate plus",
        "fluid-removal rate in CVVHD (Methods, 'Statistical analyses'). The",
        "regression was fitted on 1-3 L/h ex vivo runs and simulated over",
        "1-5 L/h. Used only when RRT_CRRT_STATUS = 1."
      ),
      source_name = "Effluent flow rate (L/h)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 162L,
    n_studies = 2L,
    disease_state = paste(
      "Simulated adults with ventilator-associated bacterial pneumonia (VABP)",
      "caused by Acinetobacter baumannii-calcoaceticus complex, either with",
      "creatinine clearance 30-129 mL/min and no CRRT, or on CRRT (CVVH or",
      "CVVHD) with effluent flow rates of 1-5 L/h."
    ),
    weight_range = "76.4 +/- 4.7 kg (simulated mean +/- SD)",
    renal_function = "Non-CRRT: CLcr mean 80 mL/min (range 30-129); CRRT: no residual renal function",
    dose_range = "Sulbactam 1 g every 4, 6, 8 or 12 h as a 3-h intravenous infusion (co-administered with durlobactam 1 g)",
    regions = "Phase 2 and Phase 3 clinical trial programmes (underlying population PK model)",
    notes = paste(
      "The base two-compartment population PK model was fitted to 162 patients",
      "from the Phase 2 and Phase 3 sulbactam-durlobactam trials (supplement,",
      "'Monte Carlo simulations'); its demographics are not reported in this",
      "paper. The non-renal clearance comes from Phase 1 patients with chronic",
      "renal insufficiency (data on file). The CRRT transmembrane clearance was",
      "measured ex vivo with heparinized bovine blood circulated through a",
      "Prismaflex circuit (M100 and HF1400 hemofilters, CVVH and CVVHD modes,",
      "effluent 1-3 L/h, blood flow 200 mL/min, triplicate runs). No human",
      "subject was studied in this paper; the authors simulated 1,000 virtual",
      "patients per scenario with the Pmetrics non-parametric simulator."
    )
  )

  ini({
    # Non-renal and renal clearance arms (supplement Table S11, 'Total Body CL' row:
    # (4.79+(0.648*13.6*((CLCR/90)^0.932)))*((Wt/75)^1.03)*(1+CLINF)). 4.79 equals
    # (1 - 0.648) x 13.6 = 4.787 to rounding.
    lcl_nonren <- log(4.79); label("Non-renal clearance (L/h)") # Table S11 'Total Body CL' constant 4.79; Table S12 CL_NR = 4.79 +/- 1.44 L/h; supplement Methods: CL_NR 4.79 L/h derived from patients with chronic renal insufficiency
    lcl <- log(13.6); label("Clearance of the underlying population PK model whose renal fraction forms the renal arm, at CLcr 90 mL/min (L/h)") # Table S11 'Total Body CL' term 13.6
    f_renal <- 0.648; label("Fraction of the population PK clearance that is renal (unitless)") # Table S11 'Total Body CL' term 0.648
    e_crcl_cl_renal <- 0.932; label("Power exponent of (CLcr / 90) on the renal clearance arm (unitless)") # Table S11 'Total Body CL' exponent 0.932
    e_wt_cl <- 1.03; label("Power exponent of (WT / 75) on total clearance, non-CRRT arm (unitless)") # Table S11 'Total Body CL' exponent 1.03
    e_vabp_cl <- -0.294; label("Proportional shift in total clearance for VABP, non-CRRT arm (fraction)") # Table S11 footnote a: CLINF = -0.294 used to simulate VABP patients

    # Central volume (Table S11 / S12 'Vc' row: 12.1*(wt/75)^0.828*(1+VcINF))
    lvc <- log(12.1); label("Central volume of distribution at WT 75 kg (L)") # Table S11 'Vc' 12.1
    e_wt_vc <- 0.828; label("Power exponent of (WT / 75) on central volume (unitless)") # Table S11 'Vc' exponent 0.828
    e_vabp_vc <- 1.38; label("Proportional shift in central volume for VABP (printed as the HABP/VABP value) (fraction)") # Table S11 footnote b: VcINF = 1.38 for HABP/VABP; Table S12 'Vc' uses (1+1.38)

    lq <- log(7.82); label("Intercompartmental clearance (L/h)") # Table S11 'Q' = 7.82; Table S12 'Q' = 7.82
    lvp <- log(6.98); label("Peripheral volume of distribution (L)") # Table S11 'Vp' = 6.98 +/- 6.66; Table S12 'Vp' = 6.98 +/- 6.66

    # CRRT transmembrane clearance regression on effluent flow
    lcl_crrt <- log(0.114); label("Intercept of the CRRT transmembrane clearance vs effluent flow regression (L/h)") # Results 'CRRT clearance studies': SUL CLTM (L/h) = 0.114 + (0.973 x effluent flow rate, L/h), adjusted R2 0.831; Fig. S1a
    e_effl_cl_crrt <- 0.973; label("Slope of CRRT transmembrane clearance on effluent flow rate (L/h per L/h)") # Results 'CRRT clearance studies': slope 0.973; Table 2 effluent coefficient 0.973 (SE 0.038)

    # Between-subject variability. The authors used the Pmetrics simulator with
    # the mean +/- SD below; carried here as log-normal with
    # omega^2 = log(CV^2 + 1) -- see vignette Assumptions and deviations.
    etalvp ~ 0.647319 # Table S11/S12 'Vp' 6.98 +/- 6.66 L; CV 95.4%
    etalcl_nonren ~ 0.086178 # Table S12 'CL_NR' 4.79 +/- 1.44 L/h; supplement Methods: 30% CV added to CL_NR (CRRT simulations only)
    etalcl_crrt ~ 0.014297 # supplement Methods: 12% CV added to CL_TM

    fu <- fixed(0.62); label("Fraction of sulbactam unbound in plasma (unitless)") # supplement Methods: protein binding fixed at 38% for sulbactam during all simulations; 1 - 0.38

    # The source is a simulation study and reports no residual error model.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
  })

  model({
    effluent_lh <- RRT_CRRT_EFFLUENT_FLOW / 1000

    # Non-CRRT arm (Table S11): renal + non-renal clearance, weight and infection scaled
    cl_renal <- f_renal * exp(lcl) * (CRCL / 90)^e_crcl_cl_renal
    cl_nocrrt <- (exp(lcl_nonren) + cl_renal) * (WT / 75)^e_wt_cl * (1 + e_vabp_cl * DIS_VABP)

    # CRRT arm (Table S12): non-renal clearance plus transmembrane clearance, each
    # with its own variability; the non-CRRT arm above uses the typical non-renal
    # clearance because Table S11 gives it no dispersion
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_crrt <- (exp(lcl_crrt) + e_effl_cl_crrt * effluent_lh) * exp(etalcl_crrt)
    cl_oncrrt <- cl_nonren + cl_crrt

    cl <- RRT_CRRT_STATUS * cl_oncrrt + (1 - RRT_CRRT_STATUS) * cl_nocrrt
    vc <- exp(lvc) * (WT / 75)^e_wt_vc * (1 + e_vabp_vc * DIS_VABP)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
