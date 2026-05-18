Oualha_2014_epinephrine <- function() {
  description <- "Population PK/PD model for continuous IV epinephrine in critically ill children following cardiopulmonary bypass for repair of congenital heart defects (Oualha 2014). One-compartment open PK with first-order elimination plus an endogenous zero-order production rate q0 and circulating-volume-anchored Vc = 0.08*WT; allometric scaling of CL and q0 on body weight (exponents fixed to 3/4). Hemodynamic Emax sub-models for heart rate (HR) and the stroke-volume * systemic-vascular-resistance product (SV*SVR) with age power effects on basal HR and SV*SVR and a RACHS-1 categorical effect on SV*SVR_max. Glucose/lactate turnover sub-model: epinephrine stimulates the zero-order plasma glucose production rate via an Emax function; plasma lactate is produced at the rate of glucose elimination and itself follows first-order elimination. kGLY and kLAC are derived at steady state (Eq. 12-13)."
  reference <- paste(
    "Oualha M, Urien S, Spreux-Varoquaux O, Bordessoule A,",
    "D'Agostino I, Pouard P, Treluyer JM (2014).",
    "Pharmacokinetics, hemodynamic and metabolic effects of",
    "epinephrine to prevent post-operative low cardiac output",
    "syndrome in children.",
    "Critical Care 18(1):R23.",
    "doi:10.1186/cc13707.",
    sep = " "
  )
  vignette <- "Oualha_2014_epinephrine"
  units <- list(time = "hour", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL and q0 with paper-fixed exponent 0.75, and linear scaling on Vc (Vc = 0.08*WT L per Eq. 5 of Oualha 2014, citing Linderkamp et al. circulating-volume formula). Cohort WT range 2.5-58 kg (median 4.5 kg), Table 1.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Paper reports age in months and parameterises HR0 and SV*SVR0 as power-of-age effects (age in months). The canonical AGE column in nlmixr2lib carries years; model() converts AGE-years to AGE-months internally (age_mo = AGE * 12). Cohort age range 0.1-189 months (median 3.9 months) per Table 1.",
      source_name        = "AGE"
    ),
    RACHS1 = list(
      description        = "RACHS-1 (Risk Adjustment for Congenital Heart Surgery) category, 2-4 in this cohort.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "2 (low risk; baseline SV*SVR_max = 0.44)",
      notes              = "Oualha 2014 pools RACHS-1 categories 3 and 4 (higher risk) versus category 2 (lower risk) and estimates two separate typical values for SV*SVR_max: 0.44 (RACHS-1 = 2) and 0.26 (RACHS-1 = 3 or 4). Decomposed inside model() into a binary indicator rachs1_high (= (RACHS1 >= 3)) that selects the high-risk shift. Cohort distribution: 16 patients in category 2, 17 in category 3, 6 in category 4 (Table 1).",
      source_name        = "RACHS-1"
    ),
    CVP = list(
      description        = "Central venous pressure",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the MAP equation additively (Eq. 7: MAP = HR * SV*SVR + CVP). Cohort median 11 mmHg (range 8-15) per Results. The model supports a per-subject column; if absent the vignette uses the cohort median 11 mmHg as a constant.",
      source_name        = "CVP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "0.1-189 months (paediatric; including 5 preterm neonates at GA 33-36 weeks)",
    age_median     = "3.9 months",
    weight_range   = "2.5-58 kg",
    weight_median  = "4.5 kg",
    sex_female_pct = 33.3,
    race_ethnicity = "Not reported (single-centre French paediatric cardiovascular intensive care unit).",
    disease_state  = "Critically ill children requiring continuous IV epinephrine after open heart surgical repair for congenital heart defects under cardiopulmonary bypass, co-administered with milrinone (0.3-0.7 ug/kg/min) to prevent low cardiac output syndrome.",
    dose_range     = "Continuous IV infusion 0.01-0.23 ug/kg/min (median 0.07) via central venous catheter; infusion duration median 1.5 days (range 1-13).",
    regions        = "France (Hopital Necker Enfants-Malades, Paris).",
    notes          = "Per Table 1 of Oualha 2014. 26 of 39 patients were male (66.6%); 11 were preoperatively cyanotic; 19 were malnourished (< 2 SD); 9 met the LCOS definition. RACHS-1 distribution: category 2 n=16, category 3 n=17, category 4 n=6. Concentrations measured at baseline (C0, n=33), after at least 60 minutes of infusion (C1, n=39), and 40 minutes after a flow-rate change or before the 6 hour cutoff (C2, n=25); 97 PK observations, 434 HR, 464 MAP, 101 glucose, 140 lactate observations."
  )

  ini({
    # ---- PK structural parameters (Table 2 of Oualha 2014) -------------
    # CL and q0 are reported as per-kg typical values evaluated at the
    # allometric reference weight 1 kg, with BW exponents FIXED to 3/4
    # (rejected after the data could not distinguish 1 from 3/4 with
    # acceptable precision; see Oualha 2014 Results paragraph 1 of
    # "Epinephrine pharmacokinetics").
    lcl <- log(2)            ; label("Typical unit clearance (L/h/kg^0.75)")            # Table 2: theta_CL = 2 L/h/kg (per the table; algebraic units L/h/kg^0.75)
    lq0 <- log(0.15)         ; label("Typical unit endogenous Ep production rate (ug/h/kg^0.75)") # Table 2: theta_q0 = 0.15 ug/h/kg
    e_wt_cl <- fixed(0.75)   ; label("Body weight allometric exponent on CL (unitless, FIXED to 3/4)") # Table 2: theta_BW(CL) = 0.75 FIX
    e_wt_q0 <- fixed(0.75)   ; label("Body weight allometric exponent on q0 (unitless, FIXED to 3/4)") # Table 2: theta_BW(q0) = 0.75 FIX

    # Vc was not separately estimated; ascribed to the circulating volume
    # via Eq. 5 (Linderkamp 1977): V_circ = 0.08 * BW (L). Encoded as a
    # FIXED log(Vc/kg) with a FIXED linear weight exponent so the
    # parameter is provenance-traceable.
    lvc <- fixed(log(0.08))  ; label("Central volume per kg (Vc / WT, L/kg, FIXED to circulating volume)") # Eq. 5: V_circ = 0.08 * BW
    e_wt_vc <- fixed(1)      ; label("Body weight allometric exponent on Vc (unitless, FIXED to 1)")      # Eq. 5

    # IIV on log(CL) and log(q0) with correlation 0.88 (Table 2). Variances
    # are SD^2 where the source reports sqrt(omega^2). Off-diagonal:
    # cov = r * SD_CL * SD_q0 = 0.88 * 1 * 1.1 = 0.968. Source: Table 2
    # eta_CL = 1, eta_q0 = 1.1, corr(eta_CL, eta_q0) = 0.88.
    etalcl + etalq0 ~ c(1, 0.968, 1.21)

    # ---- Hemodynamic sub-model (Table 3, "Hemodynamic population parameters") ----
    lhr0      <- log(133)    ; label("Typical basal heart rate at age = 1 month (b/min)")            # Table 3: HR0 = 133 b/min
    lhrmax    <- log(180)    ; label("Maximal heart rate at infinite Ep concentration (b/min)")      # Table 3: HR_max = 180 b/min
    lc50hr    <- log(5.71)   ; label("Ep concentration producing 50% of HR_max effect (ug/L)")       # Table 3: C50_HR = 5.71 ug/L
    lsvsvr0   <- log(0.31)   ; label("Typical basal SV * SVR product at age = 1 month (mmHg.min/b)") # Table 3: SV*SVR0 = 0.31 (algebraic unit mmHg per beat per minute)
    lsvsvrmax <- log(0.44)   ; label("Typical maximal SV * SVR product when RACHS-1 = 2 (reference)")# Table 3: SV*SVR_max (RACHS-1=2) = 0.44
    lc50svsvr <- log(18)     ; label("Ep concentration producing 50% of SV*SVR_max effect (ug/L)")   # Table 3: C50_SV*SVR = 18 ug/L

    # Age power effects on basal HR0 and SV*SVR0. Paper formulae (Results,
    # paragraph 2 of "Epinephrine hemodynamics"):
    #   HR0_i      = HR0     * age_i ^ (-0.061)
    #   SV*SVR0_i  = SV*SVR0 * age_i ^ ( 0.094)
    # with age in months (paper convention, not years).
    e_age_hr0    <- -0.061   ; label("Age (months) power exponent on basal HR0 (unitless)")          # Table 3 / Results: theta_age (HR0) = -0.061
    e_age_svsvr0 <-  0.094   ; label("Age (months) power exponent on basal SV*SVR0 (unitless)")      # Table 3 / Results: theta_age (SV*SVR0) = 0.094

    # RACHS-1 categorical effect on SV*SVR_max: SV*SVR_max = 0.44 when
    # RACHS-1 = 2 (reference) and 0.26 when RACHS-1 in {3, 4}. Encoded as
    # an additive shift on log(SV*SVR_max): e = log(0.26 / 0.44) =
    # -0.5263, applied only when rachs1_high = (RACHS1 >= 3).
    e_rachs1_high_svsvrmax <- log(0.26 / 0.44) ; label("Additive log-shift on SV*SVR_max for RACHS-1 in {3,4} (unitless)") # Table 3: SV*SVR_max = 0.44 (RACHS-1=2) -> 0.26 (RACHS-1=3,4)

    # IIV on log(HR0), log(C50_HR), log(SV*SVR0), log(SV*SVR_max).
    etalhr0       ~ 0.0196     # SD = 0.14   ; Table 3 eta_HR0
    etalc50hr     ~ 1.4884     # SD = 1.22   ; Table 3 eta_C50HR
    etalsvsvr0    ~ 0.0169     # SD = 0.13   ; Table 3 eta_SV*SVR0
    etalsvsvrmax  ~ 0.0529     # SD = 0.23   ; Table 3 eta_SV*SVRmax

    # ---- Metabolic sub-model (Table 3, "Metabolic population parameters") ----
    # Paper reports R_GLY in mmol/L/min; the model time unit is hours so
    # R_GLY is stored as the per-hour equivalent (0.04 * 60 = 2.4) and
    # k_GLY / k_LAC are derived in model() from the steady-state
    # assumptions (Eq. 12-13).
    lgly0   <- log(5.46)     ; label("Basal plasma glucose level (mmol/L)")                          # Table 3: GLY0 = 5.46 mmol/L
    lgmax   <- log(1.69)     ; label("Maximal stimulation of glucose zero-order production rate (unitless, multiplier)") # Table 3: G_max = 1.69
    lrgly   <- log(0.04 * 60); label("Plasma glucose zero-order production rate (mmol/L/h)")         # Table 3: R_GLY = 0.04 mmol/L/min -> 2.4 mmol/L/h
    lc50gly <- log(0.52)     ; label("Ep concentration producing 50% of glucose-production stimulation (ug/L)") # Table 3: C50_GLY = 0.52 ug/L
    llac0   <- log(1.23)     ; label("Basal plasma lactate level (mmol/L)")                          # Table 3: LAC0 = 1.23 mmol/L

    etalgly0 ~ 0.0441         # SD = 0.21   ; Table 3 eta_GLY0
    etalgmax ~ 0.0454         # SD = 0.213  ; Table 3 eta_Gmax
    etalrgly ~ 1              # SD = 1      ; Table 3 eta_RGLY
    etallac0 ~ 0.1089         # SD = 0.33   ; Table 3 eta_LAC0

    # ---- Residual error ------------------------------------------------
    # Per-output residual errors (Table 2 and Table 3). Ep, HR, and MAP
    # use a proportional model; plasma glucose and lactate use an
    # additive (constant) model.
    propSd     <- 0.3        ; label("Proportional residual SD on Ep concentration (fraction)")      # Table 2: residual proportional = 0.3
    propSd_HR  <- 0.08       ; label("Proportional residual SD on heart rate (fraction)")             # Table 3: HR proportional residual = 0.08
    propSd_MAP <- 0.16       ; label("Proportional residual SD on mean arterial pressure (fraction)") # Table 3: MAP proportional residual = 0.16
    addSd_GLY  <- 2.23       ; label("Additive residual SD on plasma glucose (mmol/L)")               # Table 3: GLY additive residual = 2.23
    addSd_LAC  <- 0.5        ; label("Additive residual SD on plasma lactate (mmol/L)")               # Table 3: LAC additive residual = 0.5
  })

  model({
    # ---- Derived covariate terms -----------------------------------------
    # Paper age covariate uses months; the canonical AGE column is in years.
    age_mo <- AGE * 12
    # Binary high-risk RACHS-1 indicator (categories 3 or 4 vs reference 2).
    rachs1_high <- (RACHS1 >= 3)

    # ---- Individual PK parameters ----------------------------------------
    # Paper parameters are per L/h/kg^0.75 etc. with WT in kg and time in
    # hours. e_wt_cl and e_wt_q0 are FIXED to 3/4 in ini().
    cl <- exp(lcl + etalcl) * WT^e_wt_cl
    q0 <- exp(lq0 + etalq0) * WT^e_wt_q0
    vc <- exp(lvc) * WT^e_wt_vc            # Vc = 0.08 * WT (L), per Eq. 5

    # ---- Individual hemodynamic parameters -------------------------------
    hr0      <- exp(lhr0      + etalhr0     ) * age_mo^e_age_hr0
    hrmax    <- exp(lhrmax)
    c50hr    <- exp(lc50hr    + etalc50hr)
    svsvr0   <- exp(lsvsvr0   + etalsvsvr0  ) * age_mo^e_age_svsvr0
    svsvrmax <- exp(lsvsvrmax + etalsvsvrmax + e_rachs1_high_svsvrmax * rachs1_high)
    c50svsvr <- exp(lc50svsvr)

    # ---- Individual metabolic parameters ---------------------------------
    gly0   <- exp(lgly0  + etalgly0)
    gmax   <- exp(lgmax  + etalgmax)
    rgly   <- exp(lrgly  + etalrgly)
    c50gly <- exp(lc50gly)
    lac0   <- exp(llac0  + etallac0)

    # Derive k_GLY and k_LAC from the steady-state assumption (Eq. 12-13).
    # Units check: rgly [mmol/L/h] / gly0 [mmol/L] = kgly [1/h]. Similarly
    # klac = rgly / lac0 [1/h] using k_LAC = k_GLY * GLY0 / LAC0 = R_GLY /
    # LAC0 (substituting Eq. 12 into Eq. 13).
    kgly <- rgly / gly0
    klac <- rgly / lac0

    # ---- Initial conditions ----------------------------------------------
    # central holds total (exogenous + endogenous) Ep amount. At baseline,
    # the endogenous production q0 balances elimination CL * Cc, giving
    # Cc(0) = q0 / cl and amount A(0) = q0 * vc / cl (Eq. 3).
    central(0) <- q0 * vc / cl
    # Glucose and lactate compartments start at steady-state basal
    # levels (Eq. 12-13). Compartment names use the canonical lowercase
    # endogenous-species form (glucose, lactate); the published-name
    # observations GLY / LAC are derived algebraically below.
    glucose(0) <- gly0
    lactate(0) <- lac0

    # ---- ODE system ------------------------------------------------------
    # PK: dA/dt = q0 + Rate - cl * Cc (Eq. 1 augmented by the endogenous
    # zero-order source q0; the exogenous infusion Rate enters via the
    # event table).
    d/dt(central) <- q0 - cl * central / vc

    # Plasma Ep concentration (observation; ug/L). Paper Eq. 3 expresses
    # Cc as A/V + q0/CL; since A(0) was set to the SS amount q0*Vc/cl,
    # central / vc already represents the total observed concentration.
    Cc <- central / vc

    # Hemodynamic sub-model (Eq. 6, 7, 8). CVP enters MAP additively; the
    # cohort median 11 mmHg is used when CVP is not supplied as a
    # per-subject covariate (handled in the vignette).
    HR    <- hr0 + (hrmax - hr0) * Cc / (Cc + c50hr)
    SVSVR <- svsvr0 + (svsvrmax - svsvr0) * Cc / (Cc + c50svsvr)
    MAP   <- HR * SVSVR + CVP

    # Metabolic sub-model: Eq. 9-10 (glucose turnover stimulated by Ep);
    # Eq. 11 (lactate produced at the rate of glucose elimination and
    # itself first-order eliminated).
    Sglu          <- 1 + gmax * Cc / (Cc + c50gly)
    d/dt(glucose) <- rgly * Sglu - kgly * glucose
    d/dt(lactate) <- kgly * glucose - klac * lactate

    # Published-name observation aliases for the lowercase compartments
    # (paper notation GLY / LAC). Carry the residual-error model on the
    # uppercase aliases so they match the published symbols.
    GLY <- glucose
    LAC <- lactate

    # ---- Observation and residual error ----------------------------------
    Cc  ~ prop(propSd)
    HR  ~ prop(propSd_HR)
    MAP ~ prop(propSd_MAP)
    GLY ~ add(addSd_GLY)
    LAC ~ add(addSd_LAC)
  })
}
