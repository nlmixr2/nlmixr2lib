Flint_2017_s_ketamine <- function() {
  description <- "Joint two-compartment S-ketamine + one-compartment S-norketamine population PK model for continuous intravenous S-ketamine infusion during prolonged sedation in pediatric intensive care patients aged 0.02-12.5 years (Flint 2017). The parent S-ketamine has two-compartment disposition (CL = 112 L/h, V1 = 7.73 L, Q = 196 L/h, V2 = 545 L at 70 kg) and feeds the active metabolite S-norketamine, modelled as one apparent central compartment with Clsnk/Fm = 53.2 L/h and Vsnk/Fm = 1 L (fixed; Fm is not identifiable). Body weight is allometrically scaled with fixed exponents 0.75 for clearances and 1.0 for volumes referenced to 70 kg; time after the first S-ketamine dose acts as a linear positive multiplier on Clsnk (0.870 percent per hour), the only retained covariate at backward elimination."
  reference   <- paste(
    "Flint RB, Brouwer CNM, Kranzlin ASC, Lie-A-Huen L, Bos AP,",
    "Mathot RAA. Pharmacokinetics of S-ketamine during prolonged",
    "sedation at the pediatric intensive care unit.",
    "Pediatr Anesth. 2017;27(11):1098-1107. doi:10.1111/pan.13239",
    sep = " "
  )
  vignette    <- "Flint_2017_s_ketamine"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at study entry (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling with reference weight 70 kg (Flint 2017 Methods Section 2.4): clearances scale as (WT/70)^0.75 and volumes scale as (WT/70)^1.0, with the exponents fixed at the canonical Holford allometric values (not estimated). The same allometric structure applies to both S-ketamine and S-norketamine apparent clearances and volumes.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 25L,
    n_studies      = 1L,
    age_range      = "0.02-12.5 years (median 0.42 years; 19/25 below 2 years)",
    age_median     = "0.42 years",
    weight_range   = "3.4-35 kg (median 7.0 kg)",
    weight_median  = "7.0 kg",
    sex_female_pct = 52.0,
    disease_state  = "Mechanically ventilated children admitted to the pediatric intensive care unit (Emma Children's Hospital PICU) requiring prolonged sedation. Predominant indication was respiratory insufficiency from lower respiratory tract infection (21/25; predominantly bronchiolitis and pneumonia); other indications were encephalopathy, seizures, post-surgical recovery, and cognitive impairment (Flint 2017 Table 2).",
    dose_range     = "Continuous intravenous infusion 0.3-3.6 mg/kg/h (neonates 0.1-1.8 mg/kg/h, infants 0.3-3.6 mg/kg/h) per the PICU sedation protocol; cumulative dose 12.2-250.1 mg/kg (median 66.1); duration of infusion 9.6-140.7 h (median 53.5). Loading dose 0.5-2 mg/kg IV push and bolus loading doses equal to one hour of the current infusion rate were administered when sedation depth was inadequate (COMFORT-B above 15) -- the model is therefore appropriate for continuous-infusion regimens with occasional intercurrent bolus loading.",
    regions        = "Single-center study at Emma Children's Hospital, Academic Medical Center, Amsterdam, The Netherlands (August 2011 - May 2012).",
    n_observations = "86 plasma concentrations of S-ketamine and S-norketamine across 25 patients (one daily random-time sample during infusion plus two post-infusion samples at 1 h and 4 h after discontinuation). 4 of 86 S-ketamine and 2 of 86 S-norketamine samples were below quantification (LLOQ 20 ng/mL and 10 ng/mL respectively) and were excluded from the final fit after a sensitivity check showed M3-method inclusion did not change the parameter estimates.",
    co_medication  = "All subjects received lorazepam as standard sedation co-medication (0.05-0.1 mg/kg q4h IV/oral plus 0.003-0.01 mg/kg/h IV). Other commonly used concomitant medications (Flint 2017 Table 2) were midazolam (16/25), chloral hydrate (14/25), rocuronium (14/25), propofol (14/25), and morphine (13/25). Concomitant medication was not retained as a covariate in the final model.",
    notes          = "Median COMFORT-B target was 10-15 (adequate sedation) and Visual Analog Scale target was below 4 (adequate analgesia); infusion rate was titrated to these clinical endpoints rather than to a target plasma concentration. Demographics in Table 2; final parameter estimates in Table 3."
  )

  ini({
    # ---------------- Parent (S-ketamine) structural parameters ----------------
    # Reference body weight 70 kg per the standard allometric anchor used
    # throughout Flint 2017 Methods Section 2.4 and Table 3. Point estimates
    # are the final-model means from Table 3, allometrically standardised to
    # 70 kg.
    lcl    <- log(112);  label("S-ketamine clearance at 70 kg (L/h)")                           # Flint 2017 Table 3: Clsk = 112 L/h/70 kg, %SE 9
    lvc    <- log(7.73); label("S-ketamine central volume of distribution at 70 kg (L)")        # Flint 2017 Table 3: V1sk = 7.73 L/70 kg, %SE 31
    lq     <- log(196);  label("S-ketamine intercompartmental clearance at 70 kg (L/h)")        # Flint 2017 Table 3: Qsk = 196 L/h/70 kg, %SE 22
    lvp    <- log(545);  label("S-ketamine peripheral volume of distribution at 70 kg (L)")     # Flint 2017 Table 3: V2sk = 545 L/70 kg, %SE 19

    # ---------------- Metabolite (S-norketamine) apparent parameters -----------
    # Apparent parameters because the fraction of S-ketamine metabolised to
    # S-norketamine, Fm, is not identifiable from the design; Flint 2017
    # estimates and reports clearance and volume as Clsnk/Fm and Vsnk/Fm
    # (Methods Section 2.4 and Table 3 footnote a). Vsnk/Fm is fixed at unity
    # because varying it from 0.1 to 10 did not change the goodness of fit or
    # the other parameter values (Results Section 3.3).
    lcl_snk <- log(53.2);       label("S-norketamine apparent clearance Clsnk/Fm at 70 kg (L/h)") # Flint 2017 Table 3: Clsnk/Fm = 53.2 L/h/70 kg, %SE 27
    lvc_snk <- fixed(log(1));   label("S-norketamine apparent central volume Vsnk/Fm at 70 kg, fixed (L)")  # Flint 2017 Table 3 footnote a: Vsnk/Fm fixed at 1

    # ---------------- Allometric covariate effects ----------------------------
    # Fixed at the canonical Holford allometric exponents (Flint 2017 Methods
    # Section 2.4): 0.75 for clearances and 1.0 for volumes; reference weight
    # 70 kg. The exponents are declared as fixed structural assumptions, not
    # estimated.
    e_wt_cl     <- fixed(0.75); label("Allometric exponent of (WT/70) on S-ketamine CL (unitless)")           # Flint 2017 Methods Section 2.4 (fixed)
    e_wt_vc     <- fixed(1.0);  label("Allometric exponent of (WT/70) on S-ketamine V1 (unitless)")           # Flint 2017 Methods Section 2.4 (fixed)
    e_wt_q      <- fixed(0.75); label("Allometric exponent of (WT/70) on S-ketamine Q (unitless)")            # Flint 2017 Methods Section 2.4 (fixed)
    e_wt_vp     <- fixed(1.0);  label("Allometric exponent of (WT/70) on S-ketamine V2 (unitless)")           # Flint 2017 Methods Section 2.4 (fixed)
    e_wt_cl_snk <- fixed(0.75); label("Allometric exponent of (WT/70) on S-norketamine Clsnk/Fm (unitless)")  # Flint 2017 Methods Section 2.4 (fixed)
    e_wt_vc_snk <- fixed(1.0);  label("Allometric exponent of (WT/70) on S-norketamine Vsnk/Fm (unitless)")   # Flint 2017 Methods Section 2.4 (fixed)

    # ---------------- Time-after-first-dose effect on Clsnk -------------------
    # Linear positive slope on S-norketamine apparent clearance: Clsnk/Fm
    # increases by 0.870 percent per hour relative to its baseline value
    # (Flint 2017 Table 3 and Results Section 3.2/3.3). Implemented inside
    # model() as a multiplicative factor (1 + e_t_cl_snk * t), where t is
    # simulation time in hours. Time is treated as the covariate because the
    # paper's parameterisation is time-after-first-dose and the model
    # convention is that t = 0 corresponds to the first S-ketamine dose.
    e_t_cl_snk <- 0.00870; label("Linear slope of time on S-norketamine Clsnk/Fm (1/h)")  # Flint 2017 Table 3: 0.870%/h, %SE 70

    # ---------------- IIV (variances on log scale) ----------------------------
    # Flint 2017 Table 3 reports interpatient variability as percentages
    # (CV%), the standard NONMEM exponential parameterisation: Clsk 40.2%,
    # Clsnk 104%, correlation r = 0.506. The log-scale variances follow the
    # log-normal-from-CV relationship omega^2 = log(CV^2 + 1). The
    # cross-covariance entry is r * sd_lcl * sd_lcl_snk where sd_X = sqrt(omega^2_X).
    #   omega^2_lcl     = log(0.402^2 + 1) = 0.14965
    #   omega^2_lcl_snk = log(1.04^2  + 1) = 0.73326
    #   sd_lcl          = sqrt(0.14965)    = 0.38685
    #   sd_lcl_snk      = sqrt(0.73326)    = 0.85631
    #   cov             = 0.506 * 0.38685 * 0.85631 = 0.16762
    # Only Clsk and Clsnk had identifiable IIV; the other parameters had no
    # IIV retained in the final model (Results Section 3.1: "the data did not
    # contain enough information to estimate interpatient variability for the
    # other parameters").
    etalcl + etalcl_snk ~ c(0.14965, 0.16762, 0.73326)  # Flint 2017 Table 3: CV(Clsk) = 40.2%, CV(Clsnk) = 104%, r = 0.506

    # ---------------- Residual error (proportional, by output) ----------------
    # Flint 2017 Section 3.1 reports a proportional error model for both
    # compounds; Table 3 gives 41.9% (S-ketamine) and 47.5% (S-norketamine).
    propSd     <- 0.419; label("S-ketamine proportional residual SD (fraction)")     # Flint 2017 Table 3: 41.9%, %SE 13
    propSd_snk <- 0.475; label("S-norketamine proportional residual SD (fraction)")  # Flint 2017 Table 3: 47.5%, %SE 11
  })

  model({
    # Reference body weight for allometric scaling (Flint 2017 Methods
    # Section 2.4 and Eq. on page 3: "P_70kg is the parameter in patient
    # with a standardized bodyweight of 70 kg").
    ref_wt <- 70

    # ---------- Individual parameters: S-ketamine (parent) ---------------------
    # Allometric scaling with fixed canonical exponents (0.75 for clearances,
    # 1.0 for volumes). IIV applied to lcl only -- Flint 2017 retained IIV
    # only on the two clearances.
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc)          * (WT / ref_wt)^e_wt_vc
    q  <- exp(lq)           * (WT / ref_wt)^e_wt_q
    vp <- exp(lvp)          * (WT / ref_wt)^e_wt_vp

    # ---------- Individual parameters: S-norketamine (apparent) ----------------
    # Time-after-first-dose multiplier on Clsnk: (1 + e_t_cl_snk * t) where t
    # is the simulation-time variable in hours. The convention is that t = 0
    # corresponds to the first S-ketamine dose, which matches Flint 2017's
    # parameterisation of "time after first dose" exactly when the event
    # table starts the infusion at time zero.
    cl_snk <- exp(lcl_snk + etalcl_snk) *
              (WT / ref_wt)^e_wt_cl_snk *
              (1 + e_t_cl_snk * t)
    vc_snk <- exp(lvc_snk) * (WT / ref_wt)^e_wt_vc_snk

    # ---------- Micro-constants -----------------------------------------------
    # Standard 2-compartment + 1-compartment micro-constants derived from
    # macro parameters (CL, V1, Q, V2 for the parent; Clsnk_app, Vsnk_app for
    # the metabolite).
    kel    <- cl / vc
    k12    <- q  / vc
    k21    <- q  / vp
    kelsnk <- cl_snk / vc_snk

    # ---------- ODE system ----------------------------------------------------
    # Two-compartment S-ketamine fed by a continuous IV infusion into
    # `central`. The flux from `central` to `central_snk` carries the
    # parent's elimination clearance times the parent concentration (mass /
    # time of S-ketamine cleared per hour). On the apparent-metabolite scale,
    # because Vsnk and Clsnk both fold the unknown fraction Fm into their
    # apparent values (Vsnk/Fm and Clsnk/Fm), the apparent flux equals the
    # parent elimination flux directly -- the Fm cancels out exactly. The
    # observed apparent S-norketamine concentration is therefore
    # central_snk / vc_snk in the same mass-per-volume units as Cc, scaled
    # implicitly by Fm.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_snk) <-  kel * central - kelsnk * central_snk

    # ---------- Observations ---------------------------------------------------
    # S-ketamine plasma concentration (mg/L = ug/mL); multiply by 1000 to get
    # ng/mL when comparing to Flint 2017 figures and assay-range statements.
    # S-norketamine concentration is reported in S-norketamine-equivalent
    # mg/L absorbing the unidentifiable Fm factor (so the observation is
    # numerically equal to what the LC-MS/MS assay returns, expressed in
    # the same units as Cc).
    Cc     <- central     / vc
    Cc_snk <- central_snk / vc_snk

    Cc     ~ prop(propSd)
    Cc_snk ~ prop(propSd_snk)
  })
}
