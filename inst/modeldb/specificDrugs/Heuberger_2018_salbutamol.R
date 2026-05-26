Heuberger_2018_salbutamol <- function() {
  description <- paste(
    "Semi-physiological PK simulation model for inhaled and oral salbutamol",
    "with its sulphate metabolite (S-SAL) in adult elite athletes. Eight",
    "compartments (gut, two-compartment parent disposition, parent plasma",
    "metabolite arm, cumulative parent urine, cumulative S-SAL urine,",
    "cumulative urine volume) with allometric scaling on disposition and",
    "physiological scaling on the cardiac-output-driven urine production",
    "rate, synthesised from literature (Auclair 2000 dog model, Morgan 1986",
    "renal CL, Holt 1968 cardiac output, Moerkeberg 2009 haematocrit) and",
    "calibrated to Haase 2009 inhaled-salbutamol data (Heuberger 2018)."
  )
  reference <- paste(
    "Heuberger JAAC, van Dijkman SC, Cohen AF. Futility of current urine",
    "salbutamol doping control. Br J Clin Pharmacol. 2018;84(8):1830-1838.",
    "doi:10.1111/bcp.13619. Final NONMEM code in Data S1 supplement."
  )
  vignette <- "Heuberger_2018_salbutamol"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg. Drives allometric scaling on renal",
        "clearances (exponent 0.75), distribution volumes (exponent 1.0),",
        "intercompartmental clearance (exponent 0.75), and cardiac output",
        "(exponent 0.79). Heuberger 2018 simulation cohort was generated",
        "as Normal(mean = 84 kg, SD = 17 kg) for elite athletes."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000,
    n_studies      = 0,
    age_range      = "Adult (not explicit in source)",
    weight_range   = "Simulated cohort: Normal(mean 84 kg, SD 17 kg)",
    weight_median  = "84 kg (simulated typical elite athlete)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy adult elite athletes (cyclists)",
    dose_range     = paste(
      "Inhalation: single 1600 ug (Haase validation) or 800 ug BID",
      "steady-state (WADA allowed); Oral: 8 mg QD x 14 days (anabolic",
      "scenario). The model accepts any inhalation or oral regimen."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Virtual cohort. Parameters synthesised from literature, not fit to",
      "individual data. Calibration target was the Haase 2009 [BJSM",
      "doi:10.1136/bjsm.2008.052522] mean concentration-time profile after",
      "a single 1600 ug salbutamol inhalation in 13 exercised, dehydrated",
      "cyclists. See Table 1 of Heuberger 2018 for the final parameter set",
      "(source labels: Holt, Moerkeberg, Auclair, Morgan)."
    )
  )

  ini({
    # All values reproduce Heuberger 2018 Table 1 (final, calibrated). The
    # supplement (Data S1) ships the NONMEM code that defines the
    # parameter structure; see comments below for the exact lineage.

    # --- Structural absorption parameters (Auclair 2000 dog model,
    # calibrated to Haase 2009 by the authors). All wrapped in fixed():
    # this is a synthesis / simulation model with no individual-data fit,
    # so the source authors held every value constant.
    lka      <- fixed(log(0.5));        label("Gut absorption rate constant ka (1/h)")            # Table 1: 0.5 1/h (Auclair, footnote a: adjusted from dog 1.5 1/h)
    lalag    <- fixed(log(1.5));        label("Gut absorption lag time ALAG (h)")                  # Table 1: 1.5 h (Auclair)
    lfdepot  <- fixed(log(0.80));       label("Gut bioavailability of inhaled dose (fraction)")    # Table 1: 80 % (Auclair). 80 % of inhaled dose is swallowed and absorbed via gut; the residual 20 % is delivered directly to central via a separate dose record (the "lung-direct" component).

    # --- Renal clearances (Morgan 1986 human IV/oral data, kept as-is).
    lcl      <- fixed(log(17.5));       label("Renal clearance of salbutamol at 70 kg (L/h)")     # Table 1: 17.5 L/h (Morgan)
    lcl_sulf <- fixed(log(5.91));       label("Renal clearance of S-SAL at 70 kg (L/h)")          # Table 1: 5.91 L/h (Morgan)

    # --- Distribution (Auclair 2000 dog model, calibrated; Table 1 values
    # are L/kg, multiplied here by 70 kg for the model's reference value).
    lvc      <- fixed(log(1.12 * 70));  label("Central volume Vc at 70 kg (L)")                    # Table 1: 1.12 L/kg (Auclair, footnote a: adjusted from dog 1.4 L/kg). V_metabolite = V_central (supplement: S3 = VC and KE_REN_SM = CL_REN_SM / VC).
    lvp      <- fixed(log(1.92 * 70));  label("Peripheral volume Vp at 70 kg (L)")                 # Table 1: 1.92 L/kg (Auclair, footnote a: adjusted from dog 4.8 L/kg)
    lq       <- fixed(log(0.56));       label("Intercompartmental clearance Q at 70 kg (L/h)")    # Table 1: 0.56 L/h (Auclair, footnote a: adjusted from dog 1.4 L/h)

    # --- Cardiac output (Holt 1968) and haematocrit (Moerkeberg 2009).
    # CO (L/min) = lco_coef * WT^0.79; reference value 0.166 at 1 kg gives
    # 0.166 * 70^0.79 = 4.78 L/min for a 70 kg adult (Methods equation 1).
    lco_coef <- fixed(log(0.166));      label("Cardiac output coefficient (L/min/kg^0.79)")        # Table 1: 0.166 (Holt). CO = 0.166 * WT^0.79.
    hct      <- fixed(0.409);           label("Typical haematocrit in trained cyclists (fraction)") # Table 1: 0.409 (Moerkeberg)

    # --- Fixed physiological constants (Methods equation 1, supplement
    # "Parameter definitions"). Not estimated. Kept inline in ini() so a
    # downstream user can refit by relaxing fixed().
    fr_kid    <- fixed(0.21);   label("Fraction of cardiac output to kidneys (unitless)")          # Methods: 21 % of CO flows through the kidneys.
    fr_bow    <- fixed(0.19);   label("Fraction of plasma entering Bowman's capsule (unitless)")   # Methods: 19 % of plasma enters renal capsule.
    fr_elim   <- fixed(0.008);  label("Fraction of glomerular filtrate becoming urine (unitless)") # Methods: 99.2 % is reabsorbed, 0.8 % leaves as urine.
    cl_hep_r  <- fixed(5.35);   label("Ratio CL_renal_salbutamol / CL_hepatic_salbutamol (unitless)") # Supplement: CL_HEP_SA = CL_REN_SA / 5.35, where 5.35 = 64.2 / 12 is the ratio of total exposure (SAL / S-SAL) after IV admin per Morgan 1986.

    # --- Allometric exponents (canonical, fixed).
    e_wt_cl  <- fixed(0.75);    label("Allometric exponent on CL and Q (unitless)")                # Supplement: (WT/70)^0.75 on CL_REN_SA, CL_REN_SM, Q3.
    e_wt_vc  <- fixed(1.0);     label("Allometric exponent on Vc and Vp (unitless)")               # Supplement: VC = THETA * WT (linear).
    e_wt_co  <- fixed(0.79);    label("Allometric exponent on cardiac output (unitless)")          # Holt 1968: CO scales as WT^0.79.

    # --- IIV (log-normal). Heuberger 2018 reports the variability as the
    # %CV percent that the simulated parameter distribution exhibits;
    # convert to the lognormal variance scale via omega^2 = log(CV^2 + 1).
    # The two parameters explicitly fixed by the source (omega = 0.05 in
    # NONMEM, per Table 1 footnote b) are encoded as fixed(0.0515) to
    # mirror that small-CV approximation: log(1 + 0.23^2) = 0.0515.
    etalco_coef ~ fixed(0.0515)  # Table 1 CV 23 % (Holt), footnote b fixed
    etahct      ~ fixed(0.00489) # Table 1 CV 7  % (Moerkeberg)
    etalfdepot  ~ fixed(0.0515)  # Table 1 CV 23 % (Auclair), footnote b fixed
    etalka      ~ fixed(0.281)   # Table 1 CV 57 % (Auclair)
    etalalag    ~ fixed(0.524)   # Table 1 CV 83 % (Auclair)
    etalcl      ~ fixed(0.0606)  # Table 1 CV 25 % (Morgan)
    etalcl_sulf ~ fixed(0.0606)  # Table 1 CV 25 % (Morgan)
    etalvc      ~ fixed(0.334)   # Table 1 CV 63 % (Auclair)
    etalvp      ~ fixed(0.223)   # Table 1 CV 50 % (Auclair)
    etalq       ~ fixed(0.128)   # Table 1 CV 37 % (Auclair)

    # --- Residual error (Table 1: proportional 23 % for both plasma
    # outputs). Heuberger 2018 does not report a residual-error model for
    # the urine outputs; we share the plasma proportional error magnitude
    # for the simulated urine outputs so a downstream user can produce a
    # VPC for the urine compartments.
    propSd       <- fixed(0.23);  label("Proportional residual SD on plasma salbutamol (fraction)")    # Table 1: 23 %
    propSd_sulf  <- fixed(0.23);  label("Proportional residual SD on plasma S-SAL (fraction)")         # Table 1: 23 %
    propSd_urineSal  <- fixed(0.23); label("Proportional residual SD on urine salbutamol (fraction)")  # not reported; shared with plasma (see vignette Assumptions)
    propSd_urineSulf <- fixed(0.23); label("Proportional residual SD on urine S-SAL (fraction)")       # not reported; shared with plasma (see vignette Assumptions)
  })

  model({
    # Physiological constants and derived urine production rate.
    # Equation 1 of the paper: UR_PROD = CO * fr_kid * (1 - HCT) * fr_bow * fr_elim (L/min).
    # Multiplied by 60 to get L/h. For a 70 kg adult with typical values
    # the equation gives 0.054 L/h ~ 1.30 L/day; the paper reports 1.2
    # L/day, in agreement.
    co_coef <- exp(lco_coef + etalco_coef)
    co      <- co_coef * WT^e_wt_co                                # L/min
    ht      <- hct * exp(etahct)
    ur_prod_h <- co * fr_kid * (1 - ht) * fr_bow * fr_elim * 60    # L/h

    # Individual PK parameters with allometric scaling (reference 70 kg).
    ka       <- exp(lka      + etalka)
    alag_t   <- exp(lalag    + etalalag)
    cl       <- exp(lcl      + etalcl)      * (WT / 70)^e_wt_cl
    cl_sulf  <- exp(lcl_sulf + etalcl_sulf) * (WT / 70)^e_wt_cl
    vc       <- exp(lvc      + etalvc)      * (WT / 70)^e_wt_vc
    vp       <- exp(lvp      + etalvp)      * (WT / 70)^e_wt_vc
    q        <- exp(lq       + etalq)       * (WT / 70)^e_wt_cl

    # Hepatic clearance of salbutamol (parent -> S-SAL in circulation),
    # derived from the renal clearance by the AUC-ratio constant 5.35
    # (supplement Parameter definitions).
    cl_hep   <- cl / cl_hep_r

    # Micro-constants.
    kel_sal       <- cl       / vc
    kel_hep_sal   <- cl_hep   / vc
    kel_sulf      <- cl_sulf  / vc        # V_metabolite = Vc per supplement
    k12           <- q        / vc
    k21           <- q        / vp

    # ODE system (supplement). The gut compartment empties at rate ka into
    # central (50 % parent) and at rate ka into central_sulf (50 % via
    # first-pass sulphation), so the effective half-life of gut emptying
    # is ln(2) / (2 * ka). The lung-absorbed 20 % of inhaled dose is
    # introduced directly to central via a separate dose record (see
    # vignette for the dosing recipe).
    d/dt(depot)        <- -2 * ka * depot
    d/dt(central)      <-      ka * depot - (kel_sal + kel_hep_sal) * central -
                                k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-                      k12 * central - k21 * peripheral1
    d/dt(central_sulf) <-      ka * depot + kel_hep_sal * central - kel_sulf * central_sulf
    d/dt(urine)        <-      kel_sal     * central
    d/dt(urine_sulf)   <-      kel_sulf    * central_sulf
    d/dt(urine_vol)    <-      ur_prod_h

    # Bioavailability and lag on the gut depot.
    f(depot)    <- exp(lfdepot + etalfdepot)
    alag(depot) <- alag_t

    # Observation variables.
    # Cc / Cc_sulf: plasma concentrations (dose in ug, vc in L ->
    # ug/L = ng/mL).
    # urineSal / urineSulf: urinary concentrations (cumulative amount in
    # ug divided by cumulative urine volume in L -> ug/L = ng/mL).
    Cc        <- central       / vc
    Cc_sulf   <- central_sulf  / vc
    urineSal  <- urine         / urine_vol
    urineSulf <- urine_sulf    / urine_vol

    Cc        ~ prop(propSd)
    Cc_sulf   ~ prop(propSd_sulf)
    urineSal  ~ prop(propSd_urineSal)
    urineSulf ~ prop(propSd_urineSulf)
  })
}
