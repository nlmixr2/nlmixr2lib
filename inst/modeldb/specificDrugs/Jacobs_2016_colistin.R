Jacobs_2016_colistin <- function() {
  description <- "Two-state parent-metabolite population PK model for colistimethate sodium (CMS, prodrug) and colistin (active polymyxin) in critically ill ICU patients with acute renal failure requiring intermittent hemodialysis (n=8). One compartment each for CMS and colistin. CMS renal clearance is structurally fixed at 0 (anuric HD population); the estimated CMS clearance is therefore nonrenal (CL_NRCMS). Colistin disposition is parameterised in apparent units (V_col/f_m and CL_col/f_m) because the fraction f_m of nonrenally cleared CMS that becomes colistin is not separately identifiable from plasma data. Hemodialysis clearances of CMS (90 mL/min) and colistin (137 mL/min) are fixed experimental constants from Marchand 2010 (ref 7 of Jacobs 2016) and are gated on/off by the time-varying DIAL covariate; PK sampling in the source study was conducted between HD sessions, so DIAL = 0 over the model-fit data."
  reference <- paste(
    "Jacobs M, Gregoire N, Megarbane B, Gobin P, Balayn D, Marchand S, Mimoz O,",
    "Couet W. Population pharmacokinetics of colistin methanesulfonate and",
    "colistin in critically ill patients with acute renal failure requiring",
    "intermittent hemodialysis. Antimicrob Agents Chemother. 2016;60(3):1788-1793.",
    "doi:10.1128/AAC.01868-15"
  )
  vignette <- "Jacobs_2016_colistin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    DIAL = list(
      description        = "Hemodialysis-active indicator (1 during a dialysis session, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = "Time-varying within subject. Gates the fixed hemodialysis-clearance contributions for CMS (90 mL/min) and colistin (137 mL/min) -- those terms are zero when DIAL = 0. The Jacobs 2016 PK sampling was performed between HD sessions, so DIAL = 0 throughout the data the model was fit to (Methods, Sample collection); HD-active dynamics are exercised in the paper's HD simulation scenarios (Methods, Simulations; Figure 4) where 4-h sessions are imposed at user-chosen times. For non-HD patients leave DIAL = 0 throughout.",
      source_name        = "DIAL"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested in forward inclusion (P<0.05) and backward deletion (P<0.01) but not retained in the final model (Jacobs 2016 Methods, Population PK modeling; Results: no covariate was included due to nonsignificant decreases of OFV)."
    ),
    AGE = list(
      description = "Subject age",
      units       = "yr",
      type        = "continuous",
      notes       = "Tested but not retained (same screen)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested but not retained. All subjects had abolished renal function (anuric); the structural model fixes CMS renal clearance at 0."
    ),
    BUN = list(
      description = "Serum urea concentration",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Tested but not retained."
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units       = "C",
      type        = "continuous",
      notes       = "Tested but not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "36-82 years (median 65)",
    weight_range   = "52-100 kg (median 80)",
    sex_female_pct = 25,
    race_ethnicity = NULL,
    disease_state  = "Critically ill ICU patients with acute renal failure requiring intermittent hemodialysis, treated with colistimethate sodium (CMS / Colimycine) as part of standard care for severe Gram-negative infection (Pseudomonas aeruginosa, Acinetobacter baumannii, Klebsiella pneumoniae). All eight subjects had abolished CMS renal clearance, so the estimated CMS clearance is structurally nonrenal.",
    dose_range     = "First dose 0.4-9 MIU CMS (median 1.5 MIU); maintenance 0.4-2 MIU q8h (median 0.5 MIU q8h) administered as 1-h IV infusions of CMS reconstituted in 50 mL saline. 1 MIU CMS sodium is equivalent to 80 mg of CMS sodium (Jacobs 2016 Results). Four subjects also received CMS aerosol cotreatment; the aerosol pathway (about 9% of dose reaching systemic, 1.4% presystemically converted to colistin per Marchand 2014 ref 12) is not represented in this nlmixr2lib model.",
    regions        = "Two-site multicenter France: CHU Poitiers, Hopital Lariboisiere Paris.",
    notes          = "Demographics from Jacobs 2016 Table 1: 2 women, 6 men (25% female). Serum creatinine range 172-470 micromol/L (median 310). SAPS II 39-75 (median 58). Intermittent hemodialysis sessions: 4-h duration every 2 days using a 1.6 m^2 B3 polymethylmethacrylate membrane on a Gambro AK 200, blood flow 300 mL/min, dialysate flow 500 mL/min."
  )

  ini({
    # Structural population PK parameters from Jacobs 2016 Table 2, ICU-HD column.
    # The paper reports clearances in mL/min; converted to L/h via x60/1000.
    #   113 mL/min  = 6.78  L/h    (CL_NRCMS)
    #    33.3 mL/min = 1.998 L/h   (CL_col/f_m)
    lvc      <- log(21);    label("CMS systemic distribution volume V_CMS (L)")                                  # Table 2: 21 L (13% RSE)
    lcl      <- log(6.78);  label("CMS nonrenal clearance CL_NRCMS (L/h); CL_RCMS structurally fixed at 0")     # Table 2: 113 mL/min (14% RSE); CL_RCMS = 0 fixed
    lvc_col  <- log(28.3);  label("Colistin apparent distribution volume V_col/f_m (L)")                         # Table 2: 28.3 L (18% RSE)
    lcl_col  <- log(1.998); label("Colistin apparent clearance CL_col/f_m (L/h)")                                # Table 2: 33.3 mL/min (16% RSE)

    # Inter-individual variability. Jacobs 2016 Table 2 reports %CV; the
    # log-normal random-effect model implies omega^2 = log(CV^2 + 1).
    #   24% CV -> 0.0561
    #   31% CV -> 0.0918
    #   42% CV -> 0.1626
    etalvc      ~ 0.0561  # IIV V_CMS, 24% CV (RSE 45%) -- Table 2
    etalcl      ~ 0.0918  # IIV CL_NRCMS, 31% CV (RSE 35%) -- Table 2
    etalvc_col  ~ 0.1626  # IIV V_col/f_m, 42% CV (RSE 36%) -- Table 2
    etalcl_col  ~ 0.1626  # IIV CL_col/f_m, 42% CV (RSE 29%) -- Table 2

    # Residual error: combined additive + proportional, separate for CMS and colistin
    # (Jacobs 2016 Methods, Population PK modeling; Table 2). 0.04 mg/L is the LOQ
    # for both analytes and is handled via Beal M3 in the paper.
    propSd      <- 0.48;  label("CMS proportional residual SD (fraction)")           # Table 2: 48% (RSE 12%)
    addSd       <- 0.11;  label("CMS additive residual SD (mg/L)")                   # Table 2: 0.11 mg/L (RSE 28%)
    propSd_col  <- 0.15;  label("Colistin proportional residual SD (fraction)")      # Table 2: 15% (RSE 24%)
    addSd_col   <- 0.13;  label("Colistin additive residual SD (mg/L)")              # Table 2: 0.13 mg/L (RSE 27%)
  })

  model({
    # Hemodialysis-clearance constants. Fixed device/experimental values from
    # Marchand 2010 (ref 7 of Jacobs 2016, same CMS brand and HD apparatus),
    # not estimated in this paper. Encoded as numeric constants in model() to
    # match the convention used by other CRRT / HD models in the library for
    # fixed circuit-level constants (Leuppi-Taegtmeyer 2019, Liesenfeld 2013
    # mass-transfer coefficient). Units L/h: 90 mL/min = 5.40, 137 mL/min = 8.22.
    cl_hd_cms <- 5.40   # CL_HDCMS, L/h
    cl_hd_col <- 8.22   # CL_HDCOL, L/h

    # Individual PK parameters
    vc      <- exp(lvc + etalvc)
    cl      <- exp(lcl + etalcl)
    vc_col  <- exp(lvc_col + etalvc_col)
    cl_col  <- exp(lcl_col + etalcl_col)

    # Total clearances. The HD contribution is gated by DIAL (1 during an active
    # intermittent-HD session, 0 otherwise). Colistin formation from CMS is driven
    # by the NONRENAL CMS-clearance arm only (cl, not cl_tot): CMS removed by
    # dialysis leaves the body in the dialysate and does not become colistin.
    cl_tot      <- cl     + DIAL * cl_hd_cms
    cl_col_tot  <- cl_col + DIAL * cl_hd_col

    # Parent CMS and apparent colistin compartments. The colistin state is the
    # "apparent" mass = A_col / f_m, where f_m is the unknown CMS-to-colistin
    # conversion fraction (Jacobs 2016 Methods). Because both volume and clearance
    # for colistin scale by 1/f_m, the observed concentration
    # Ccol = central_col / vc_col is invariant to f_m and matches the published
    # plasma colistin reading.
    d/dt(central)     <- -cl_tot * central / vc
    d/dt(central_col) <-  cl * central / vc - cl_col_tot * central_col / vc_col

    Cc      <- central     / vc
    Cc_col  <- central_col / vc_col
    Cc      ~ add(addSd)     + prop(propSd)
    Cc_col  ~ add(addSd_col) + prop(propSd_col)
  })
}
