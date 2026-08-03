Kleiber_2017_clonidine <- function() {
  description <- paste0(
    "One-compartment population PK model for intravenous clonidine in ",
    "critically ill neonates and children on venovenous or venoarterial ",
    "extracorporeal membrane oxygenation (ECMO) with concomitant continuous ",
    "venovenous hemofiltration (CVVH). Central compartment only with IV bolus ",
    "and IV infusion dosing (no absorption). Clearance carries the standard ",
    "allometric fixed exponent 0.75 on body weight (reference 70 kg), a ",
    "Hill-type maturation function of postnatal age (steep exponent 3.02 and ",
    "T50 = 1.13 weeks -- ~70% of mature clearance by 10 days of PNA), and a ",
    "multiplicative diuretic-use effect (CL x 0.659 when any diuretic is ",
    "active). Volume of distribution carries the standard fixed allometric ",
    "exponent 1 on body weight and a sigmoidal Emax effect of time-on-ECMO ",
    "(Emax = +55%, T50 = 51.7 h, Hill exponent 18.5 -- effectively a step ",
    "increase reached by 72 h on ECMO). Residual error is proportional only."
  )
  reference <- paste(
    "Kleiber N, Mathot RAA, Ahsman MJ, Wildschut ED, Tibboel D, de Wildt SN.",
    "Population pharmacokinetics of intravenous clonidine for sedation",
    "during paediatric extracorporeal membrane oxygenation and continuous",
    "venovenous hemofiltration.",
    "Br J Clin Pharmacol. 2017;83(6):1227-1239.",
    "doi:10.1111/bcp.13235.",
    sep = " "
  )
  vignette <- "Kleiber_2017_clonidine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "clonidine", units = "ug", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The most recently measured pre-ECMO body weight was used for both dose calculation and PK analysis (paper Methods, Clinical parameters). Time-fixed per subject in the source analysis. Standard fixed allometric exponents 0.75 on CL and 1 on V, centred at 70 kg (paper Methods, Eqs 4-5).",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Kleiber 2017 reports PNA and T50_PNA in weeks (Eq 8, Table 3); the canonical PNA carries months. model() converts PNA (months) to pna_weeks = PNA * 4.345 before entering the Hill maturation function so the paper-reported T50_PNA = 1.13 weeks and Hill exponent 3.02 stay in paper-natural units.",
      source_name        = "PNA"
    ),
    CONMED_DIURETIC = list(
      description        = "Binary indicator of any-diuretic coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diuretic exposure)",
      notes              = "Time-varying per-time-point indicator. Kleiber 2017 pools four diuretic classes into a single DIURETIC column: furosemide intermittent (63.6% of patients), furosemide infusion (18%), spironolactone (31.8%), and bumetanide (9.1%). Multiplicative effect on CL: CL x Theta2^DIURETIC with Theta2 = 0.659 (Table 3 / Table 4), i.e. 34.1% reduction in CL when any diuretic is active.",
      source_name        = "DIURETIC"
    ),
    T_ECMO = list(
      description        = "Time since ECMO cannulation start",
      units              = "hour",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; 0 before ECMO cannulation, monotonically increasing thereafter. Kleiber 2017 does not reverse the covariate at decannulation; the sigmoidal Emax effect saturates around 72 h so post-decannulation values sit in the saturated region (paper Discussion: 'Time dependent changes in clearance or volume after ECMO decannulation were not detected'). Enters as a sigmoidal Emax multiplier on V: V x (1 + Emax * t_ECMO^Hill / (T50_EC^Hill + t_ECMO^Hill)) with Emax = 0.55, T50_EC = 51.7 h, Hill = 18.5 (Table 3, Eq 9). At T_ECMO = 0 the multiplier evaluates to 1 (baseline V).",
      source_name        = "tec"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "3 days to 6 years postnatal at ECMO start (paper Discussion: youngest 3 days, oldest 6 years); gestational age at birth median 38.9 weeks (IQR 5.6)",
    age_median     = "1 month postnatal (IQR 6.4 months); postmenstrual age median 42.8 weeks (IQR 17.6)",
    weight_range   = "not reported explicitly",
    weight_median  = "4 kg (IQR 3.1)",
    sex_female_pct = 50.0,
    race_ethnicity = NULL,
    disease_state  = "Critically ill neonates and children on ECMO for respiratory or cardiac failure. Primary diagnoses (Table 2): pulmonary 31.8%, meconium aspiration syndrome 22.7%, cardiac 18.2%, congenital diaphragmatic hernia 13.6%, sepsis 9.1%, persistent pulmonary hypertension of the newborn 4.5%.",
    dose_range     = "IV clonidine infusion 0.1-1 ug/kg/h (median infusion dose 0.24 ug/kg/h, IQR 0.15); IV boluses 1-2 ug/kg per physician discretion (occasionally administered as a 10-min slow infusion when hemodynamically unstable). Simulation scenarios in the paper use 5 ug/kg boluses.",
    regions        = "Single centre, Erasmus MC-Sophia Pediatric Intensive Care Unit, Rotterdam, Netherlands.",
    notes          = "May 2007-July 2009. ECMO modality 68.2% VV / 31.8% VA (Table 2). Median ECMO duration 6.2 days (IQR 7.1); 90.9% of patients received concomitant CVVH at median flow 300 mL/min. Diuretic exposure common: 72.7% of patients received any diuretic. Concomitant sedatives were midazolam and morphine (all patients); some received pentobarbital, propofol, or ketamine. Survival 77.3%. 375 plasma clonidine samples (median 12.5 per patient, IQR 13); 4.2% below LOQ were ignored during model building. Clonidine measured by LC-MS/MS on serum; validated range 0.100-20.0 ug/L. Neonatal ECMO system used for children <10-12 kg (1.5 m^2 oxygenator, 350 mL prime); paediatric system used for larger children (2.5 m^2 oxygenator, 900 mL prime)."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Table 3 final-model column. Standardised to a
    # reference weight of 70 kg per the a priori allometric scaling
    # (paper Methods, Eqs 4-5). Doses in ug, volumes in L, clearances in L/h
    # so the mass-balance ug / L = ng/mL closes without unit conversion.
    # ----------------------------------------------------------------------
    lcl <- log(29.9); label("Clearance CL at reference WT = 70 kg, full PNA maturation, no diuretic exposure (L/h)")  # Table 3 final model
    lvc <- log(454);  label("Central volume of distribution V at reference WT = 70 kg, pre-ECMO baseline (L)")         # Table 3 final model

    # ----------------------------------------------------------------------
    # Allometric exponents -- paper Methods, Eqs 4-5: parameter values were
    # standardised a priori to a bodyweight of 70 kg using fixed exponents
    # 0.75 on CL and 1 on V. The a priori allometric scaling produced
    # dOFV = -35.3 (P < 0.001, Table 4). Not estimated in NONMEM.
    # ----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric (WT) exponent on CL (unitless)")  # paper Methods Eq 4
    e_wt_vc  <- fixed(1);    label("Allometric (WT) exponent on V (unitless)")   # paper Methods Eq 5

    # ----------------------------------------------------------------------
    # Postnatal-age Hill maturation of CL -- Eq 8 / Table 3 final model.
    # Cl_i = Cl_pop * (WT/70)^0.75 * PNA^Theta1 / (T50_PNA^Theta1 + PNA^Theta1) * Theta2^DIURETIC.
    # Paper reports PNA and T50_PNA in weeks; the canonical PNA column carries
    # months, so model() converts to weeks internally via pna_weeks = PNA * 4.345
    # and the coefficients below stay in paper-natural weeks units.
    # ----------------------------------------------------------------------
    hill_pna_cl <- 3.02; label("Hill exponent for PNA maturation of CL (Theta1, unitless)")            # Table 3 final model
    t50_pna_cl  <- 1.13; label("Postnatal age at 50% of mature CL (T50_PNA, weeks)")                    # Table 3 final model

    # ----------------------------------------------------------------------
    # Diuretic effect on CL -- Eq 6 form (Table 3 / Table 4 multivariate step 4).
    # Multiplicative on CL: Theta2^DIURETIC with Theta2 = 0.659 (Table 3 RSE 7.5%),
    # i.e. 34.1% reduction in CL when any diuretic is active (paper Table 4
    # remark: "Cl decreased by 34% with use of diuretics").
    # ----------------------------------------------------------------------
    e_diur_cl <- 0.659; label("Multiplicative CL ratio when any diuretic is active (Theta2, unitless)")  # Table 3 final model

    # ----------------------------------------------------------------------
    # On-ECMO sigmoidal Emax effect on V -- Eq 9 (Table 3 final model).
    # V_i = V_pop * (WT/70) * (1 + Emax * t_ECMO^Theta3 / (T50_EC^Theta3 + t_ECMO^Theta3)).
    # Steep Hill exponent (18.5) makes the transition essentially a step: below
    # ~40 h the multiplier is near 1, above ~65 h it saturates near 1.55; paper
    # Discussion: "Maximal VD was reached 72 h after initiation of ECMO".
    # ----------------------------------------------------------------------
    emax_tec_v <- 0.55; label("Maximal fractional increase in V during ECMO (Emax, unitless)")       # Table 3 final model
    t50_tec_v  <- 51.7; label("Time on ECMO at 50% of maximal V effect (T50_EC, hour)")              # Table 3 final model
    hill_tec_v <- 18.5; label("Hill exponent for t_ECMO effect on V (Theta3, unitless)")             # Table 3 final model

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Table 3 final model.
    # Reported as CV% (40% on CL, 44% on V). Variances on the log scale are
    # computed from CV via omega^2 = log(CV^2 + 1):
    #   log(0.40^2 + 1) = log(1.16)   = 0.148420
    #   log(0.44^2 + 1) = log(1.1936) = 0.176974
    # IIV encoded as diagonal (independent etas); the source paper did not
    # report an OMEGA-block correlation between CL and V.
    # ----------------------------------------------------------------------
    etalcl ~ 0.148420  # IIV variance on log CL (40% CV, Table 3 final model)
    etalvc ~ 0.176974  # IIV variance on log V (44% CV, Table 3 final model)

    # ----------------------------------------------------------------------
    # Residual error -- Table 3 final model. Proportional only (Table 3
    # reports "Proportional 0.208 (RSE 6.5%)" for the final model; the
    # structural model reported 0.281). The paper's Discussion says "the
    # residual error decreased from 28% to 20%", matching these two values.
    # The paper's initial-model description (Methods, Eq 3) mentions both
    # proportional and additive terms, but only the proportional was retained
    # in the final fit (the additive was dropped during model building --
    # the paper text says "The proportional error in the residual error
    # model was very small and could be eliminated", which is inconsistent
    # with Table 3 showing the proportional term present in both structural
    # and final models; the Discussion percentages match the proportional
    # row, so the paper text is treated as a typo and the additive term is
    # taken as the one that was eliminated).
    # ----------------------------------------------------------------------
    propSd <- 0.208; label("Proportional residual error on Cc (fraction)")  # Table 3 final model
  })

  model({
    # ------------------------------------------------------------------
    # Reference (centring) values.
    # ------------------------------------------------------------------
    ref_wt <- 70  # kg -- allometric a priori standardisation weight (paper Eqs 4-5)

    # ------------------------------------------------------------------
    # Postnatal-age Hill maturation of CL. Canonical PNA is in months;
    # convert to weeks (paper's units) so T50_PNA and Hill_PNA stay in
    # the paper-natural values.
    # ------------------------------------------------------------------
    pna_weeks  <- PNA * 4.345
    mat_cl_pna <- pna_weeks^hill_pna_cl / (t50_pna_cl^hill_pna_cl + pna_weeks^hill_pna_cl)

    # ------------------------------------------------------------------
    # Diuretic effect on CL (Eq 6). Theta2^DIURETIC evaluates to 1 when
    # CONMED_DIURETIC = 0 and to Theta2 = 0.659 when CONMED_DIURETIC = 1.
    # ------------------------------------------------------------------
    diur_cl <- e_diur_cl^CONMED_DIURETIC

    # ------------------------------------------------------------------
    # On-ECMO sigmoidal Emax multiplier on V (Eq 9). At T_ECMO = 0 the
    # ratio T_ECMO^18.5 / (T50^18.5 + T_ECMO^18.5) = 0 so the multiplier
    # collapses to 1 (baseline V); for T_ECMO >> T50 it saturates at
    # 1 + Emax = 1.55. Post-decannulation samples keep the ECMO effect
    # because the Hill has saturated and the paper does not describe
    # reversal.
    # ------------------------------------------------------------------
    tec_v <- 1 + emax_tec_v * T_ECMO^hill_tec_v / (t50_tec_v^hill_tec_v + T_ECMO^hill_tec_v)

    # ------------------------------------------------------------------
    # Individual parameters.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl * mat_cl_pna * diur_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc  * tec_v

    # ------------------------------------------------------------------
    # Micro-constant and ODE. Central compartment only (IV bolus + IV
    # infusion; no absorption).
    # ------------------------------------------------------------------
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # Observation and residual error. Dose in ug, Vc in L => Cc in
    # ug/L = ng/mL, matching the paper's plasma-concentration units.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
