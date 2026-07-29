Mohamed_2012_colistin <- function() {
  description <- "Two-compartment population PK model for colistin methanesulfonate (CMS, prodrug) plus a one-compartment apparent model for colistin (formed metabolite) in critically ill patients, with concentration-dependent unbound fraction of colistin A and a semimechanistic Pseudomonas aeruginosa bacterial-kill PKPD (susceptible / resting compartments from Bulitta 2010)"
  reference <- paste(
    "Mohamed AF, Karaiskos I, Plachouras D, Karvanen M, Pontikis K, Jansson B,",
    "Papadomichelakis E, Antoniadou A, Giamarellou H, Armaganidis A, Cars O, Friberg LE.",
    "Application of a Loading Dose of Colistin Methanesulfonate in Critically Ill Patients:",
    "Population Pharmacokinetics, Protein Binding, and Prediction of Bacterial Kill.",
    "Antimicrob Agents Chemother. 2012 Aug;56(8):4241-9.",
    "doi:10.1128/AAC.06426-11.",
    "Semimechanistic bacterial-kill PKPD parameters are inherited from Mohamed",
    "AF, Cars O, Friberg LE (2011) PAGE abstract 2223 (the reference 29 in the",
    "main paper) and are restated in the Materials and Methods of the 2012",
    "publication.",
    sep = " "
  )
  vignette <- "Mohamed_2012_colistin"
  paper_specific_compartments <- c(
    "central_cms", "peripheral_cms", "central_col", "bact_s", "bact_r"
  )
  paper_specific_residual_sds <- c(
    "addSd_Ccms", "propSd_Ccms", "addSd_Cc_col", "propSd_Cc_col"
  )

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The published final model retained no covariates: a stepwise covariate
  # screen (gender, body weight, ideal body weight, age, serum creatinine,
  # creatinine clearance, serum albumin, APACHE II, septicemic state,
  # hemoglobin, hematocrit) failed to meet the dOFV > 10.83 inclusion
  # threshold (Mohamed 2012, Materials and Methods, "Population
  # pharmacokinetic modeling"; Results "PK model" and Table 3).
  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the stepwise covariate model (SCM); not retained. Body weight range 60-140 kg (Table 1).",
      source_name        = "WT"
    ),
    IBW = list(
      description        = "Ideal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in SCM; not retained. Ideal body weights 60-110 kg (Table 1).",
      source_name        = "IBW"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in SCM; not retained. Age range 32-88 years, mean 55.4 (Results).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened in SCM; not retained. 6 males / 4 females in the n=10 prospective cohort.",
      source_name        = "GENDER"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate measured on days 1, 7, 14, and 21 (Methods). Screened in SCM; not retained.",
      source_name        = "CREAT"
    ),
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cockcroft-Gault CrCL; capped at 130 mL/min (7.8 L/h) and modeled as time-varying when explored (Methods).",
        "Renal function correlation with CL_CMS was explored exhaustively (Eqs 1 and 2 in the paper).",
        "Reduction in OFV was only 3.6 (P > 0.05); CrCL did not meet the predefined dOFV > 10.83 inclusion criterion and was not retained (Results, Table 3)."
      ),
      source_name        = "CRCL"
    ),
    ALB = list(
      description        = "Serum albumin",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in SCM; not retained. Baseline albumin 1.9-3.8 g/dL (Table 1).",
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 2L,
    age_range      = "32-88 years",
    age_median     = "55.4 years (mean of the n=10 prospective cohort)",
    weight_range   = "60-140 kg (actual); 60-110 kg (ideal)",
    sex_female_pct = 40,
    disease_state  = "critically ill adults with infections caused by multidrug-resistant Gram-negative bacteria (most commonly ventilator-associated pneumonia); APACHE II scores 7-23",
    dose_range     = "Loading dose 480 mg (6 MU) CMS as 15-min IV infusion, then maintenance 80-240 mg (1-3 MU) q8h or 30-90 mg q12h (Materials and Methods)",
    regions        = "Greece (prospective cohort, Attikon University Hospital, Athens)",
    notes          = paste(
      "Combined analysis of the n=10 loading-dose cohort enrolled at Attikon",
      "University Hospital July 2009 - January 2010 (Table 1) plus the n=18",
      "patients from Plachouras et al. 2009 (Antimicrob Agents Chemother",
      "53:3430-3436, reference 34), giving N=28 patients with 267 CMS",
      "and colistin concentration-time samples taken after the first and",
      "eighth doses. Exclusion: continuous venovenous hemodiafiltration",
      "for renal replacement therapy."
    )
  )

  ini({
    # ============================================================
    # CMS pharmacokinetics (Table 2; combined n=28 final fit). The
    # paper reports IIV / IOV as the omega-on-log-scale magnitude
    # using the small-omega "CV%-equals-SD-on-log-scale" convention
    # (the 1.8x scaling relationship between CL_CMS IIV (42%) and
    # CL_col IIV (76%) only reconciles under this convention).
    # ============================================================
    lcl   <- log(13.1);  label("CMS clearance CL_CMS (L/h)")                            # Table 2 typical value (90% CI 11.6-15.1)
    lvc   <- log(11.8);  label("CMS central volume V1 (L)")                              # Table 2 typical value (90% CI 7.6-12.5)
    lq    <- log(206);   label("CMS intercompartmental clearance Q (L/h)")               # Table 2 typical value (90% CI 93-641)
    lvp   <- log(28.4);  label("CMS peripheral volume V2 (L)")                           # Table 2 typical value (90% CI 22.4-33.5)

    # ============================================================
    # Colistin pharmacokinetics (Table 2). Both reported parameters
    # are scaled to the unknown fraction f_m of CMS that converts to
    # colistin (Materials and Methods, "Population pharmacokinetic
    # modeling"). Using the apparent-amount parametrization in
    # model() recovers the observed (real) colistin concentration
    # without needing f_m.
    # ============================================================
    lcl_col <- log(8.2); label("Apparent colistin clearance CL_col / f_m (L/h)")         # Table 2 typical value (90% CI 6.0-11.2)
    lvc_col <- log(218); label("Apparent colistin volume V_col / f_m (L)")               # Table 2 typical value (90% CI 192-252)

    # ============================================================
    # Inter-individual variability (Table 2). The paper reports IIV
    # only for CL_CMS (42%), Q (111%), and CL_col (76% = 1.8 x 42).
    # The CL_CMS / CL_col correlation was estimated at 100% with
    # the colistin SD scaled from the CMS SD by a factor of 1.8;
    # implemented in model() as cl_col <- exp(lcl_col + 1.8 * etalcl)
    # so the two clearances share a single random effect.
    # ============================================================
    etalcl ~ 0.42^2     # IIV(CL_CMS) = 42 %; Table 2
    etalq  ~ 1.11^2     # IIV(Q)     = 111 %; Table 2

    # ============================================================
    # Inter-occasion variability (Table 2). For library simulation
    # without an explicit OCC column, IOV is approximated as added
    # log-scale random effects on the affected parameters (see
    # vignette "Assumptions and deviations" for the trade-off).
    # ============================================================
    etalvp     ~ 0.59^2  # IOV(V2_CMS) = 59 %; Table 2 (used as a per-subject random effect for simulation)
    etalcl_col ~ 0.48^2  # IOV(CL_col)  = 48 %; Table 2
    etalvc_col ~ 0.48^2  # IOV(V_col)   = 48 %; Table 2

    # ============================================================
    # Residual error (Table 2). Reported additive component is in
    # umol/L; converted to mg/L in model() using molar masses
    # 1743 g/mol (CMS) and 1163 g/mol (colistin) reported in the
    # Analytical method paragraph.
    # ============================================================
    propSd_Ccms   <- 0.23;   label("CMS proportional residual SD (fraction)")               # Table 2 (90% CI 0.12-0.21)
    addSd_Ccms    <- 0.071;  label("CMS additive residual SD (umol/L)")                     # Table 2 (90% CI 0.063-0.25)
    propSd_Cc_col <- 0.082;  label("Colistin proportional residual SD (fraction)")          # Table 2 (90% CI 0.065-0.12)
    addSd_Cc_col  <- 0.044;  label("Colistin additive residual SD (umol/L)")                # Table 2 (90% CI 0.018-0.053)

    # ============================================================
    # Plasma protein binding (Results "Plasma protein binding";
    # Eq 3). Colistin A unbound fraction is a saturable function
    # of total colistin A concentration; colistin B unbound
    # fraction is constant at 0.43. CA / CB ratio averaged 3.6 in
    # patients (range 2.5-4.8); the unbound colistin driving the
    # bacterial kill is the sum unb(CA) + unb(CB).
    # ============================================================
    fuA_max  <- fixed(0.312); label("Maximum unbound fraction of colistin A (Eq 3)")     # Eq 3 fitted parameter
    fuA_KCA  <- fixed(0.094); label("Colistin A concentration at 50% of fuA_max (mg/L)")  # Eq 3 fitted parameter
    fuB      <- fixed(0.43);  label("Constant unbound fraction of colistin B")            # Results "Plasma protein binding"
    ratio_AB <- fixed(3.6);   label("Average ratio of colistin A to colistin B")          # Results "Plasma protein binding"

    # ============================================================
    # Semimechanistic bacterial-kill PKPD against wild-type
    # P. aeruginosa ATCC 27853 (MIC 1 mg/L), Materials and Methods
    # "Predictions of bacterial kill" -- structure and parameters
    # carried in from Bulitta 2010 / Mohamed 2011 (reference 29 in
    # the paper, PAGE abstract 2223). All values are restated
    # explicitly in the current paper.
    # ============================================================
    Emax_kill <- fixed(35);     label("Maximum colistin kill rate (1/h)")                 # Methods, Predictions of bacterial kill
    EC50_kill <- fixed(2.9);    label("Unbound colistin for 50% Emax (mg/L)")             # Methods, Predictions of bacterial kill
    kgrow     <- fixed(0.99);   label("Bacterial growth rate constant (1/h)")             # Methods, Predictions of bacterial kill
    kdeath    <- fixed(0.18);   label("Natural bacterial death rate constant (1/h)")      # Methods, Predictions of bacterial kill
    konR      <- fixed(6.0e-5); label("Colistin-driven emergence of resistance (L/mg/h)") # Methods, Predictions of bacterial kill
    koffR     <- fixed(0.15);   label("Reversal of resistance rate constant (1/h)")       # Methods, Predictions of bacterial kill
    CFUmax    <- fixed(1.8e8);  label("Maximum bacterial count at stationary phase (CFU/mL)") # Methods, Predictions of bacterial kill
    CFU0      <- fixed(4.5e5);  label("Initial inoculum (CFU/mL)")                        # Methods, Predictions of bacterial kill
  })

  model({
    # ----- Molar-mass constants (Materials and Methods, Analytical method) -----
    # 1 umol/L = MW * 1e-3 mg/L. CMS 1743 g/mol, colistin 1163 g/mol.
    addSd_Ccms_mgL   <- addSd_Ccms   * 1.743e-3
    addSd_Cc_col_mgL <- addSd_Cc_col * 1.163e-3

    # ----- CMS individual PK parameters (apparent volumes; CL/F not relevant since IV) -----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # ----- Colistin individual apparent PK parameters -----
    # CL_col IIV (76%) = 1.8 x CL_CMS IIV (42%) with perfect correlation
    # (Table 2 footnote). Implemented as a shared eta scaled by 1.8.
    cl_col <- exp(lcl_col + 1.8 * etalcl + etalcl_col)
    vc_col <- exp(lvc_col + etalvc_col)

    # ----- 2-cmt CMS + 1-cmt apparent colistin -----
    # CMS->colistin flux into the apparent-amount compartment is the full
    # CL_CMS * C_CMS (no f_m factor needed; see vignette source trace).
    d/dt(central_cms)    <- -cl * (central_cms / vc) - q * (central_cms / vc) + q * (peripheral_cms / vp)
    d/dt(peripheral_cms) <-  q * (central_cms / vc) - q * (peripheral_cms / vp)
    d/dt(central_col)    <-  cl * (central_cms / vc) - cl_col * (central_col / vc_col)

    Ccms   <- central_cms / vc       # CMS plasma concentration (mg/L)
    Cc_col <- central_col / vc_col   # Observed (real) total colistin concentration (mg/L)

    # ----- Concentration-dependent protein binding (Eq 3) -----
    CA_tot  <- Cc_col * (ratio_AB / (1.0 + ratio_AB))  # total colistin A in plasma (mg/L)
    CB_tot  <- Cc_col * (1.0      / (1.0 + ratio_AB))  # total colistin B in plasma (mg/L)
    fuA_now <- fuA_max * CA_tot / (fuA_KCA + CA_tot)
    Cu_col  <- fuA_now * CA_tot + fuB * CB_tot       # unbound colistin (mg/L) driving bacterial kill

    # ----- Bulitta semimechanistic bacterial-kill model -----
    plat    <- 1.0 - (bact_s + bact_r) / CFUmax
    kill    <- Emax_kill * Cu_col / (EC50_kill + Cu_col)
    SR_rate <- konR * Cu_col                          # S -> R rate constant (1/h)

    d/dt(bact_s) <- kgrow * plat * bact_s - kdeath * bact_s - kill * bact_s - SR_rate * bact_s + koffR * bact_r
    d/dt(bact_r) <-                       - kdeath * bact_r                + SR_rate * bact_s - koffR * bact_r

    bact_s(0) <- CFU0
    bact_r(0) <- 0.0

    # log10 bacterial count for plotting (no residual error -- deterministic prediction).
    logCFU <- log(bact_s + bact_r) / log(10.0)

    # ----- Observation models -----
    Ccms   ~ add(addSd_Ccms_mgL)   + prop(propSd_Ccms)
    Cc_col ~ add(addSd_Cc_col_mgL) + prop(propSd_Cc_col)
  })
}
