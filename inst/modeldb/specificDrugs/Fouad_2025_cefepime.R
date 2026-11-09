Fouad_2025_cefepime <- function() {
  description <- paste(
    "Two-compartment population PK simulation model for cefepime (given as",
    "cefepime-taniborbactam) in adults, with and without continuous renal",
    "replacement therapy (CRRT) (Fouad 2025). Without CRRT, clearance is the",
    "sponsor population PK clearance, a power function of creatinine clearance",
    "with end-stage renal disease and complicated urinary tract infection",
    "shifts. On CRRT the renal clearance is replaced by the transmembrane",
    "clearance measured in an ex vivo CRRT circuit, a linear function of the",
    "effluent flow rate (CLTM = 0.41 + 0.86 x effluent, L/h), added to a",
    "non-renal clearance of 0.87 L/h. Unbound concentration uses a fixed 17.3%",
    "protein binding."
  )
  reference <- paste(
    "Fouad A, McGrath C, Buyukyanbolu E, Kuti JL. Ex vivo",
    "assessment and simulation to guide cefepime-taniborbactam dosing",
    "recommendations for patients receiving continuous renal replacement",
    "therapy. Antimicrob Agents Chemother. 2025;69(6):e00061-25.",
    "doi:10.1128/aac.00061-25.",
    "CRRT clearance regression from the Results ('Impact of CRRT variables on",
    "CLTM'); population PK parameters, non-renal clearance, dispersion and",
    "protein binding from the supplement ('Monte Carlo Simulations and Dose",
    "Recommendations').",
    sep = " "
  )
  vignette <- "Fouad_2025_cefepime_taniborbactam_crrt"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Baseline creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives clearance as (CRCL / 77)^0.696 in the non-CRRT arm (supplement,",
        "'CL (L/h) = 4.75*(BCLCR/77)^0.696*...'). The estimating equation is not",
        "stated. The paper simulated three groups at 60, 90 and 120 mL/min, each",
        "with a 15% CV. Ignored when RRT_CRRT_STATUS = 1, because the CRRT",
        "simulations replace the renal clearance by the transmembrane clearance",
        "and do not account for residual renal function (Discussion)."
      ),
      source_name = "BCLCR"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator (1 = ESRD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ESRD)",
      notes = paste(
        "Multiplies the non-CRRT clearance by (1 - 0.613) (supplement CL",
        "equation). The ESRD threshold of the sponsor model is not stated. The",
        "paper's simulations used ESRD = 0. Not used when RRT_CRRT_STATUS = 1."
      ),
      source_name = "ESRD"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection / acute pyelonephritis indicator (1 = cUTI/AP)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteers and other phase 1 participants)",
      notes = paste(
        "Multiplies the non-CRRT clearance by (1 + 0.281) (supplement CL",
        "equation, 'cUTI/AP'). The source pools cUTI and acute pyelonephritis",
        "into one indicator. The paper's control simulations applied the",
        "factor (the printed 5.1, 6.78 and 8.29 L/h include x 1.281), so set",
        "DIS_CUTI = 1 to reproduce them. Not used when RRT_CRRT_STATUS = 1."
      ),
      source_name = "cUTI/AP"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator (1 = on CRRT)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT; sponsor population PK clearance)",
      notes = paste(
        "Selects between the two clearance models of the source. 0: total CL =",
        "4.75 x (CRCL/77)^0.696 x (1 - 0.613)^ESRD x (1 + 0.281)^cUTI/AP. 1:",
        "total CL = CL_NR + CL_TM, with CL_NR 0.87 L/h (30% CV) and CL_TM from",
        "the effluent-flow regression. Volumes and Q are shared by both arms."
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
        "The regression was fitted on 2-4 L/h ex vivo runs (CVVH and CVVHD,",
        "ST150 and HF1400 hemofilters, blood flow 200 mL/min) and simulated over",
        "2-5 L/h. Used only when RRT_CRRT_STATUS = 1."
      ),
      source_name = "Effluent flow rate (L/h)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 544L,
    n_studies = 7L,
    disease_state = paste(
      "Simulated adults with hospital-acquired or ventilator-associated",
      "pneumonia, either with creatinine clearance 60-120 mL/min and no CRRT,",
      "or on CRRT (CVVH or CVVHD) at effluent flow rates of 2-5 L/h. The",
      "underlying population PK model was fitted to healthy participants,",
      "participants with renal impairment, and patients with complicated",
      "urinary tract infection."
    ),
    renal_function = "Non-CRRT: CLcr 60, 90 or 120 mL/min (15% CV); CRRT: no residual renal function",
    dose_range = "Cefepime 1 or 2 g every 8 or 12 h as a 4-h intravenous infusion (co-administered with taniborbactam 0.25 or 0.5 g)",
    regions = "Six phase 1 trials and one phase 3 trial of the cefepime-taniborbactam programme (underlying population PK model)",
    notes = paste(
      "The base population PK model was fitted to 607 adults (544 for",
      "cefepime) from six phase 1 trials and one phase 3 trial (sponsor report",
      "ICPD 00443-1, data on file); its demographics are not reported in this",
      "paper. The non-renal clearance comes from phase 1 participants with",
      "renal impairment (data on file). The CRRT transmembrane clearance was",
      "measured ex vivo with heparinized bovine blood circulated through a",
      "Prismaflex circuit (ST150 and HF1400 hemofilters, CVVH and CVVHD modes,",
      "effluent 2-4 L/h, blood flow 200 mL/min, duplicate runs). No human",
      "subject was studied in this paper; the authors simulated 1,000 virtual",
      "patients per scenario with the Pmetrics simulator."
    )
  )

  ini({
    # Non-CRRT clearance (supplement: CL (L/h) = 4.75*(BCLCR/77)^0.696*(1-0.613)^ESRD*(1+0.281)^cUTI/AP)
    lcl <- log(4.75); label("Clearance without CRRT at CLcr 77 mL/min (L/h)") # Supplement Monte Carlo section, cefepime CL equation, 4.75
    e_crcl_cl <- 0.696; label("Power exponent of (CRCL / 77) on clearance (unitless)") # Supplement cefepime CL equation, exponent 0.696
    e_esrd_cl <- -0.613; label("Fractional change in clearance with ESRD (fraction)") # Supplement cefepime CL equation, (1-0.613)^ESRD
    e_cuti_cl <- 0.281; label("Fractional change in clearance with cUTI/AP (fraction)") # Supplement cefepime CL equation, (1+0.281)^cUTI/AP

    lvc <- log(8.49); label("Central volume of distribution (L)") # Supplement Monte Carlo section, cefepime 'Vc (L), 8.49'
    lq <- log(5.94); label("Intercompartmental clearance (L/h)") # Supplement Monte Carlo section, cefepime 'Q (L/h), 5.94'
    lvp <- log(6.26); label("Peripheral volume of distribution (L)") # Supplement Monte Carlo section, cefepime 'Vp (L), 6.26'

    # CRRT arm: non-renal clearance plus ex vivo transmembrane clearance
    lcl_nonren <- log(0.87); label("Non-renal clearance (L/h)") # Supplement: cefepime CL_NR 0.87 L/h (patients with renal impairment, VENATORX-008)
    lcl_crrt <- log(0.41); label("Intercept of the CRRT transmembrane clearance vs effluent flow regression (L/h)") # Results 'Impact of CRRT variables on CLTM': cefepime CLTM (L/h) = 0.41 + (0.86 x effluent flow rate, L/h)
    e_effl_cl_crrt <- 0.86; label("Slope of CRRT transmembrane clearance on effluent flow rate (L/h per L/h)") # Results 'Impact of CRRT variables on CLTM': slope 0.86; Table 2 effluent coefficient 0.86 (SE 0.05)

    # Between-subject variability used in the Pmetrics simulations, carried as
    # log-normal with omega^2 = log(CV^2 + 1) (see vignette Assumptions)
    etalvc ~ 0.254211 # Supplement: CV 53.8% for Vc (ICPD 00443-1)
    etalvp ~ 0.0039111 # Supplement: CV '6.26%' for Vp (ICPD 00443-1), as printed
    etalq ~ 0.0861777 # Supplement: 30% CV assumed for Q for both drugs
    etalcl_nonren ~ 0.0861777 # Supplement: 30% CV added to CL_NR (CRRT simulations)

    fu <- fixed(0.827); label("Fraction of cefepime unbound in plasma (unitless)") # Supplement: protein binding fixed at 17.3% for cefepime during all simulations; 1 - 0.173

    # The source is a simulation study and reports no residual error model.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
  })

  model({
    effluent_lh <- RRT_CRRT_EFFLUENT_FLOW / 1000

    # Non-CRRT arm: sponsor population PK clearance
    cl_nocrrt <- exp(lcl) * (CRCL / 77)^e_crcl_cl *
      (1 + e_esrd_cl * RENALIMP_ESRD) * (1 + e_cuti_cl * DIS_CUTI)

    # CRRT arm: renal clearance replaced by the ex vivo transmembrane clearance
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_crrt <- exp(lcl_crrt) + e_effl_cl_crrt * effluent_lh
    cl_oncrrt <- cl_nonren + cl_crrt

    cl <- RRT_CRRT_STATUS * cl_oncrrt + (1 - RRT_CRRT_STATUS) * cl_nocrrt
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
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
