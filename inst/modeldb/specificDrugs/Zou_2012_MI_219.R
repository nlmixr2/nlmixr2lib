Zou_2012_MI_219 <- function() {
  description <- "Predicted-human two-compartment IV PK model for MI-219 (a small-molecule HDM2/p53 inhibitor) in adults, with parameters projected from NONMEM-based interspecies allometric scaling of single-dose IV plasma profiles in rats (5 mg/kg), beagle dogs (2 mg/kg), and cynomolgus monkeys (10 mg/kg). Linear elimination from the central compartment; mouse data were excluded from the joint NONMEM fit because the mouse profile was not superimposable on the other species under Wajima / Dedrick normalisation. The model file encodes the predicted human typical values at a 70 kg reference body weight (Zou 2012 Table 5, NONMEM column)."
  reference <- paste(
    "Zou P, Zheng N, Yu Y, Yu S, Sun W, McEachern D, Yang Y, Yu LX,",
    "Wang S, Sun D. (2012). Preclinical pharmacokinetics of MI-219, a",
    "novel human double minute 2 (HDM2) inhibitor and prediction of",
    "human pharmacokinetics. J Pharm Pharm Sci 15(2):265-280.",
    "doi:10.18433/j34s4n.",
    sep = " "
  )
  vignette <- "Zou_2012_MI_219"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species          = "human (predicted from preclinical interspecies allometric scaling using rat, dog, and monkey IV PK)",
    n_subjects       = NA_integer_,
    n_studies        = 0L,
    age_range        = NA_character_,
    weight_range     = "70 kg reference body weight",
    weight_median    = "70 kg",
    sex_female_pct   = NA_real_,
    disease_state    = "Predicted first-in-human PK for the HDM2 / p53 small-molecule inhibitor MI-219; no human subjects were dosed in Zou 2012.",
    dose_range       = "5 mg/kg IV bolus (Zou 2012 Figure 5 and Figure 7 simulation dose used for the human-PK projection).",
    regions          = NA_character_,
    notes            = "Predicted human PK parameters were obtained by joint NONMEM (Version VII) allometric fit of single-dose IV plasma profiles in rat (Sprague-Dawley, n = 3, 5 mg/kg, 195 +/- 15 g), beagle dog (n = 3, 2 mg/kg, 8.0 +/- 1.1 kg), and cynomolgus monkey (n = 6, 10 mg/kg, 2.5 +/- 0.26 kg). The structural model is two-compartment with linear elimination from the central compartment and body-weight allometry P = PTV * WT^b applied to every PK parameter. Per-species PTV and exponent b are reported in Zou 2012 Table S1 of the supplement; the supplement is not on disk for this extraction. Mouse data (CD-1, n = 33, 10 mg/kg, 18 +/- 1.7 g) were collected but excluded from the joint NONMEM fit because the mouse plasma profile was not superimposable on rat / dog / monkey curves under Wajima or Dedrick normalisation (Zou 2012 Results, page 274). The parameter encoding below is the predicted human PK at 70 kg from Zou 2012 Table 5, NONMEM column."
  )

  ini({
    # Structural PK parameters - Zou 2012 Table 5 (NONMEM column) predicted
    # human values for a 70 kg adult. The paper reports the two-compartment
    # macroconstants alpha, beta, A, B together with the marginal point
    # estimates CL (0.186 L/h/kg) and Vdss (0.843 L/kg). The microconstant
    # form encoded here (CL, Vc, Q, Vp) is derived deterministically from
    # the four reported macroconstants and the simulated 5 mg/kg IV bolus
    # dose (D = 350 mg in a 70 kg subject) using the standard 2-compartment
    # IV-bolus closed form:
    #
    #   Vc  = D / (A + B)
    #   k21 = (A*beta + B*alpha) / (A + B)
    #   k10 = alpha * beta / k21
    #   k12 = alpha + beta - k10 - k21
    #   CL  = k10 * Vc       Q  = k12 * Vc       Vp = Q / k21
    #
    # With A = 21714 ng/mL, B = 3109 ng/mL, alpha = 2.839 /h, beta = 0.161 /h,
    # D = 350 mg:
    #   Vc = 14.10 L (0.2014 L/kg)
    #   k21 = 0.4964 /h    k10 = 0.9208 /h    k12 = 1.583 /h
    #   CL = 12.98 L/h (0.1854 L/h/kg; matches Table 5 NONMEM CL = 0.186 L/h/kg
    #                                  within rounding)
    #   Q  = 22.32 L/h (0.3189 L/h/kg)
    #   Vp = 44.96 L  (0.6422 L/kg)
    #   Vss = Vc + Vp = 59.06 L (0.8437 L/kg; matches Table 5 NONMEM
    #                            Vdss = 0.843 L/kg within rounding)
    #
    # No IIV or residual variability is encoded because the NONMEM Omega and
    # Sigma estimates are reported only in Zou 2012 Table S1 of the
    # supplement, which is not on disk for this extraction. The model is a
    # typical-value forward-prediction; see vignette Assumptions and deviations.
    lcl <- log(12.98)   ; label("Predicted human clearance (L/h) at 70 kg")                       # derived from Zou 2012 Table 5 (NONMEM column): alpha = 2.839, beta = 0.161, A = 21714, B = 3109, dose 5 mg/kg
    lvc <- log(14.10)   ; label("Predicted human central volume (L) at 70 kg")                    # derived from Zou 2012 Table 5: Vc = Dose / (A + B), Dose = 5 mg/kg * 70 kg = 350 mg
    lq  <- log(22.32)   ; label("Predicted human intercompartmental clearance Q (L/h) at 70 kg")  # derived from Zou 2012 Table 5: Q  = k12 * Vc, k12 = 1.583 /h
    lvp <- log(44.96)   ; label("Predicted human peripheral volume Vp (L) at 70 kg")              # derived from Zou 2012 Table 5: Vp = Q / k21,  k21 = 0.4964 /h
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment IV with linear elimination from the central
    # compartment (Zou 2012 Methods: "a two-compartment pharmacokinetic
    # (PK) model with linear elimination from the central compartment was
    # developed to describe the concentration-time profile of MI-219 in
    # different animal species"). Doses enter `central` directly.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Plasma concentration; dose units mg, Vc units L, so Cc units mg/L
    # (= ug/mL = 1000 ng/mL). To compare against Zou 2012 Figure 7 (which
    # plots ng/mL), multiply Cc by 1000.
    Cc <- central / vc
  })
}
