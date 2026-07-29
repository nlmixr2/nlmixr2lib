Abduljalil_2009_clarithromycin <- function() {
  description <- "Semimechanistic population pharmacokinetic model for oral clarithromycin and its 14-(R)-hydroxy metabolite during repeated b.i.d. administration (Abduljalil 2009): a single-phase Weibull absorption (kw, lambda) into a one-compartment parent disposition with linear distribution and a parent clearance that is partly inhibited by the parent's own concentration in a hypothetical effect-style inhibition compartment (Imax form with FCLp = fraction of CLp not subject to inhibition and IC50 = inhibition-compartment concentration giving 50% of maximum inhibition); all parent metabolic clearance feeds a parallel one-compartment metabolite disposition (14-OH-clarithromycin). Body weight enters allometrically with fixed exponents 0.75 on CL and 1.0 on V (parent and metabolite), reference 70 kg."
  reference <- paste(
    "Abduljalil K, Kinzig M, Bulitta J, Horkovics-Kovats S, Sorgel F,",
    "Rodamer M, Fuhr U. Modeling the autoinhibition of clarithromycin",
    "metabolism during repeated oral administration. Antimicrob Agents",
    "Chemother. 2009 Jul;53(7):2892-2901. doi:10.1128/AAC.01193-08.",
    "PMID 19414575.",
    sep = " "
  )
  vignette <- "Abduljalil_2009_clarithromycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight; allometrically scales CL and V for parent and metabolite",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Reference 70 kg with fixed allometric exponents 0.75 on CL (parent and metabolite) and 1.00 on V (parent and metabolite); see Abduljalil 2009 Methods (Population pharmacokinetic analysis, Model development paragraph 4): CL_i = CL * (BW/70)^0.75 and V_i = V * (BW/70).",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "19-41 years (mean 28, SD 8)",
    age_median     = "approximately 28 years (mean)",
    weight_range   = "45.1-86.1 kg (mean 66.5, SD 11.8)",
    weight_median  = "66.5 kg (mean)",
    height_range   = "150.0-186.0 cm (mean 168.4, SD 9.7)",
    sex_female_pct = 41.7,
    race_ethnicity = c(White = 100),
    disease_state  = "Healthy adult Caucasian volunteers (7 men, 5 women), nonsmokers or former smokers, judged healthy on medical history, vital signs, physical examination, neurology, 12-lead ECG, and clinical chemistry / hematology / urine / virology screens.",
    dose_range     = "Clarithromycin 500 mg oral suspension (Klacid, Abbott) every 12 h for 7 doses (4 consecutive days) in the fasting state with 240 mL low-carbonation calcium-poor mineral water.",
    regions        = "Republic of Moldova (single-centre study, approved by the Ethics Committee of the Ministry of Health Clinic Hospital of the Republic of Moldavia, Chisinau).",
    notes          = "Demographics from Abduljalil 2009 Methods, Subjects and treatments paragraph. A total of 624 clarithromycin and 624 14-(R)-hydroxy-clarithromycin plasma samples were comodeled by NONMEM V FOCE-I, ADVAN6 differential-equation system. Grapefruit products, methylxanthine-containing foods and beverages, alcohol, and fatty foods were restricted around the study."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification of each item):
  # * Absorption: Abduljalil 2009 Methods (Model development paragraph 2)
  #   describes a single-phase Weibull cumulative absorbed fraction
  #     WB(Tw) = 1 - exp[-(kw*Tw)^lambda]
  #   where Tw is time since the previous dose. The corresponding
  #   instantaneous absorption rate constant is
  #     k_abs(Tw) = lambda * kw * (kw*Tw)^(lambda - 1)
  #   so the amount remaining to be absorbed obeys
  #     d/dt(depot) = -k_abs(tad) * depot
  #   with tad = rxode2 tad(), the time since the most recent dose event.
  #   With lambda = 2.23 > 1 the rate starts at zero at tad = 0 and
  #   builds up to a maximum around the inflection of the Weibull
  #   cumulative curve, capturing the delayed onset of absorption that
  #   the paper attributes to intestinal CYP inactivation during
  #   first-pass metabolism. For multi-dose simulations tad resets at
  #   each dose, so each dose follows its own Weibull profile starting
  #   from tad = 0; this is the same approximation Abduljalil 2009 used
  #   in NONMEM ADVAN6.
  # * Parent disposition: one-compartment with apparent volume Vp/F and
  #   nonlinear apparent clearance CLp(t)/F. The fraction not subject
  #   to inhibition (FCLp) is bounded in [0, 1] and is encoded with a
  #   logit transform (logitfclp) so the parameter cannot wander
  #   outside the published feasible range during simulation.
  # * Inhibition compartment: modelled as an effect-compartment-style
  #   state holding the inhibition-driving concentration in the same
  #   units as the parent plasma concentration (ug/mL). The first-order
  #   exchange constant ki = 2.01 1/h corresponds to a t_1/2 of
  #   ln(2)/2.01 = 0.345 h (Abduljalil 2009 Discussion paragraph 3
  #   reports "approximately 0.35 h"). The driving concentration is the
  #   parent plasma concentration Cp = central / vc. The overall
  #   inhibition factor is
  #     INH = FCLp + (1 - FCLp) / (1 + effect / IC50)
  #   so that INH = 1 (no inhibition) when effect = 0 and INH = FCLp
  #   (maximum inhibition) when effect >> IC50. The apparent clearance
  #   acting on the parent compartment at any instant is CLp * INH.
  # * Metabolite formation: Abduljalil 2009 Discussion paragraph 4
  #   states "the final model assumption was that all clarithromycin
  #   molecules are converted to the metabolite", so the full parent
  #   apparent clearance CLp * INH * Cp is also the metabolite formation
  #   rate (mass per time) entering the metabolite central compartment.
  #   Metabolite elimination is linear (CL_met).
  # * Allometric WT scaling: fixed exponents 0.75 on CL (parent and
  #   metabolite) and 1.00 on V (parent and metabolite) with reference
  #   weight 70 kg, encoded with `fixed()` per Abduljalil 2009
  #   Methods (Model development paragraph 4): "Individual body weights
  #   were related to a standard weight, i.e., 70 kg". The paper
  #   reports no fitted exponents for these terms.
  # * IIV scale: NONMEM reports OMEGA on the internal log scale and
  #   Abduljalil 2009 Table 1 lists between-subject variability as
  #   percent values. The packaged model uses the exact log-normal
  #   conversion omega^2 = log(CV^2 + 1) so the simulated cross-subject
  #   geometric coefficient of variation matches the reported value.
  # * Residual error: Abduljalil 2009 Results paragraph 4 reports
  #   "residual additive error models resulted in the best model
  #   convergence" for both analytes; the Table 1 additive error of
  #   0.12 ug/mL (parent) and 0.01 ug/mL (metabolite) is encoded as
  #   `addSd` and `addSd_ohcla` respectively.
  ini({
    # ----- Structural parameters (Abduljalil 2009 Table 1) -----
    # Weibull absorption. lka here is the canonical absorption-rate
    # log parameter; the Weibull form differs from a first-order ka
    # but the same canonical name is reused because the role
    # ("absorption rate constant") matches.
    lka      <- log(0.56);  label("Weibull absorption rate constant kw (1/h)")         # Abduljalil 2009 Table 1: kw = 0.56 (95% CI 0.42-0.69)
    llambda  <- log(2.23);  label("Weibull shape parameter lambda (unitless)")         # Abduljalil 2009 Table 1: lambda = 2.23 (95% CI 1.67-2.77)

    # Parent disposition. CLp here is the basic (uninhibited) apparent
    # total clearance; the effective clearance at any instant is
    # cl * INH where INH is the inhibition factor (see model()).
    lvc      <- log(172);   label("Apparent parent central volume Vp/F at 70 kg (L)") # Abduljalil 2009 Table 1: Vp = 172 (95% CI 145-198)
    lcl      <- log(60);    label("Apparent parent basic (uninhibited) clearance CLp/F at 70 kg (L/h)")  # Abduljalil 2009 Table 1: CLp = 60 (95% CI 40-80)

    # Autoinhibition mechanism.
    logitfclp <- log(0.10 / (1 - 0.10));  label("Logit of fraction of CLp not subject to inhibition (FCLp; unitless)")  # Abduljalil 2009 Table 1: FCLp = 0.10 (95% CI 0.02-0.17)
    lki      <- log(2.01);  label("Inhibition-compartment first-order exchange rate ki (1/h)")  # Abduljalil 2009 Table 1: ki = 2.01 (95% CI 0.09-3.93); corresponds to t1/2 = ln(2)/ki ~ 0.35 h
    lic50    <- log(0.77);  label("Inhibition-compartment concentration giving 50% of maximum inhibition IC50 (ug/mL)")  # Abduljalil 2009 Table 1: IC50 = 0.77 (95% CI 0.23-1.28)

    # Metabolite disposition (14-(R)-hydroxy-clarithromycin).
    lcl_ohcla <- log(50.2); label("Apparent 14-OH-clarithromycin clearance CL_met/F at 70 kg (L/h)")  # Abduljalil 2009 Table 1: CL_met = 50.2 (95% CI 42.3-58.1)
    lvc_ohcla <- log(34);   label("Apparent 14-OH-clarithromycin central volume V_met/F at 70 kg (L)")  # Abduljalil 2009 Table 1: V_met = 34 (95% CI 12-56)

    # ----- Allometric exponents (fixed) -----
    e_wt_cl   <- fixed(0.75);  label("Allometric WT exponent on parent CL (unitless)")  # Abduljalil 2009 Methods (Model development paragraph 4); standard adult allometric exponent, not estimated
    e_wt_vc   <- fixed(1.00);  label("Allometric WT exponent on parent Vc (unitless)")  # Abduljalil 2009 Methods (Model development paragraph 4); standard adult allometric exponent, not estimated
    e_wt_cl_ohcla <- fixed(0.75); label("Allometric WT exponent on metabolite CL (unitless)")  # Abduljalil 2009 Methods (Model development paragraph 4)
    e_wt_vc_ohcla <- fixed(1.00); label("Allometric WT exponent on metabolite Vc (unitless)")  # Abduljalil 2009 Methods (Model development paragraph 4)

    # ----- Between-subject variability (Abduljalil 2009 Table 1) -----
    # Exact log-normal conversion omega^2 = log(CV^2 + 1).
    etalka         ~ 0.18664   # Abduljalil 2009 Table 1: CV = 45.3%; log(1 + 0.453^2) = 0.18664
    etalvc         ~ 0.06203   # Abduljalil 2009 Table 1: CV = 25.3%; log(1 + 0.253^2) = 0.06203
    etalcl         ~ 0.02981   # Abduljalil 2009 Table 1: CV = 17.4%; log(1 + 0.174^2) = 0.02981
    etalcl_ohcla   ~ 0.07502   # Abduljalil 2009 Table 1: CV = 27.9%; log(1 + 0.279^2) = 0.07502

    # ----- Residual error (Abduljalil 2009 Table 1) -----
    addSd        <- 0.12;     label("Additive residual SD on clarithromycin plasma concentration (ug/mL)")  # Abduljalil 2009 Table 1: additive error of clarithromycin = 0.12 ug/mL
    addSd_ohcla  <- 0.01;     label("Additive residual SD on 14-OH-clarithromycin plasma concentration (ug/mL)")  # Abduljalil 2009 Table 1: additive error of 14-OH-clarithromycin = 0.01 ug/mL
  })
  model({
    # ----- Reference weight -----
    ref_wt <- 70  # Abduljalil 2009 Methods (Model development paragraph 4)

    # ----- Individual parameters (allometric WT scaling, exponents fixed) -----
    kw      <- exp(lka     + etalka)
    lambda  <- exp(llambda)
    vc      <- exp(lvc     + etalvc) * (WT / ref_wt)^e_wt_vc
    cl      <- exp(lcl     + etalcl) * (WT / ref_wt)^e_wt_cl
    ki      <- exp(lki)
    ic50    <- exp(lic50)
    fclp    <- exp(logitfclp) / (1 + exp(logitfclp))
    cl_ohcla <- exp(lcl_ohcla + etalcl_ohcla) * (WT / ref_wt)^e_wt_cl_ohcla
    vc_ohcla <- exp(lvc_ohcla)               * (WT / ref_wt)^e_wt_vc_ohcla

    # ----- Driving concentrations -----
    Cp        <- central / vc
    Cc_ohcla  <- central_ohcla / vc_ohcla

    # ----- Autoinhibition factor -----
    # INH = FCLp + (1 - FCLp) / (1 + effect / IC50): equals 1 (no
    # inhibition) when effect = 0, approaches FCLp (the non-inhibitable
    # residual fraction) as effect grows. `effect` is a concentration
    # (ug/mL) driven by Cp via a first-order exchange with rate ki, in
    # the Sheiner-Holford effect-compartment form.
    inh <- fclp + (1 - fclp) / (1 + effect / ic50)

    # ----- Weibull absorption rate constant (1/h) -----
    # k_abs(tad) = lambda * kw * (kw*tad)^(lambda - 1). With lambda > 1
    # the rate is 0 at tad = 0 and rises through the absorption window.
    kw_tad      <- kw * tad()
    k_abs       <- lambda * kw * kw_tad^(lambda - 1)

    # ----- ODE system -----
    # Parent absorption + linear distribution + inhibited elimination.
    d/dt(depot)          <- -k_abs * depot
    d/dt(central)        <-  k_abs * depot - cl * inh * Cp

    # Effect-compartment-style inhibition driver. Holds the inhibition
    # concentration (ug/mL) that lags Cp with first-order rate ki.
    d/dt(effect)         <-  ki * (Cp - effect)

    # Metabolite formation = parent's inhibited apparent clearance flux.
    d/dt(central_ohcla)  <-  cl * inh * Cp - cl_ohcla * Cc_ohcla

    # ----- Observations -----
    # Both outputs in ug/mL (dose mg, vc L, central mg gives mg/L = ug/mL).
    Cc  <- Cp
    Cc           ~ add(addSd)
    Cc_ohcla     ~ add(addSd_ohcla)
  })
}
