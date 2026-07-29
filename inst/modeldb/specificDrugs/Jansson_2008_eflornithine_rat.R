Jansson_2008_eflornithine_rat <- function() {
  description <- "Preclinical (rat, Sprague-Dawley). Stereoselective two-enantiomer population PK model of racemic eflornithine after single oral or IV doses in male Sprague-Dawley rats (Jansson 2008). Each enantiomer (L = active, D) carries its own 2-compartment disposition (CL, Vc) with shared Q and Vp from the IV fit; oral absorption is modeled with a shared Savic 2007 transit-compartment chain (continuous number of compartments via Stirling approximation) feeding per-enantiomer depots that drain to central via saturable Michaelis-Menten kinetics (Tmax, Kt). Bioavailability differs between enantiomers and shifts upward at the highest oral dose level (3000 mg/kg) via a categorical indicator. Racemic plasma concentration is the algebraic sum Cc_rac = Cc_l + Cc_d."
  reference <- "Jansson R, Malm M, Roth C, Ashton M. Enantioselective and nonlinear intestinal absorption of eflornithine in the rat. Antimicrob Agents Chemother. 2008;52(8):2842-2848. doi:10.1128/aac.00050-08"
  vignette <- "Jansson_2008_eflornithine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_HIGH_EFL = list(
      description       = "High-dose oral eflornithine indicator (1 = 3000 mg/kg oral dose record; 0 = 750-2000 mg/kg oral or any IV dose).",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (750-2000 mg/kg oral, or IV at 375 / 1000 mg/kg)",
      notes             = "Per-dose-record indicator that scales bioavailability F upward for both enantiomers at the highest dose level. Jansson 2008 Results paragraph reports OFV drop -11.4 (P < 0.01) for the categorical encoding vs linear or power dose-F relationships. Apply at the dose record; observation rows inherit the indicator from the preceding dose.",
      source_name       = "derived from administered dose"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 69L,
    n_studies      = 1L,
    sex_female_pct = 0,
    weight_range   = "260-320 g",
    weight_median  = "approximately 290 g (study midpoint)",
    disease_state  = "healthy male Sprague-Dawley rats",
    dose_range     = "Oral 750, 1500, 2000, or 3000 mg/kg of racemic eflornithine hydrochloride solution by gavage (10 mL/kg); IV 375 or 1000 mg/kg as 3-min infusion modeled as bolus via the jugular vein catheter (3.3 mL/kg).",
    regions        = "Sweden (Goteborg University)",
    notes          = "n = 69 male Sprague-Dawley rats split across 6 dose-route arms (Jansson 2008 Table 1): 10 oral 750 mg/kg + 9 oral 1500 mg/kg + 7 oral 2000 mg/kg + 10 oral 3000 mg/kg + 5 IV 375 mg/kg + 5 IV 1000 mg/kg for chiral analysis (1-3 samples / rat), plus 4 oral 750 + 4 oral 1500 + 4 oral 2000 + 4 oral 3000 + 4 IV 375 + 3 IV 1000 for racemic analysis (8-16 samples / rat). Animals acclimatized at least 5 days; jugular vein catheters tunneled subcutaneously; allowed 16 h recovery after surgery. Chiral assay LLOQ 83 uM (300 uL plasma); racemic assay LLOQ 5 uM (50 uL plasma). Eflornithine MW 182.17 g/mol."
  )

  ini({
    # Two-compartment disposition (IV fit, Jansson 2008 Table 2).
    # Per-enantiomer CL and Vc; Q and Vp shared between L and D.
    # Values for a typical 290 g rat: CL = paper_value_mL_per_min_per_kg * 60 / 1000 * 0.290 (L/h);
    # V = paper_value_L_per_kg * 0.290 (L); Q analogous to CL.
    lcl_l <- log(0.2523); label("L-eflornithine clearance CL_L at 290 g rat (L/h)")        # Jansson 2008 Table 2: CL_L = 14.5 mL/min/kg; 14.5 * 60 / 1000 * 0.290 = 0.2523 L/h
    lcl_d <- log(0.2192); label("D-eflornithine clearance CL_D at 290 g rat (L/h)")        # Jansson 2008 Table 2: CL_D = 12.6 mL/min/kg; 12.6 * 60 / 1000 * 0.290 = 0.2192 L/h
    lvc_l <- log(0.1172); label("L-eflornithine central volume Vc_L at 290 g rat (L)")     # Jansson 2008 Table 2: Vc_L = 0.404 L/kg; 0.404 * 0.290 = 0.1172 L
    lvc_d <- log(0.1163); label("D-eflornithine central volume Vc_D at 290 g rat (L)")     # Jansson 2008 Table 2: Vc_D = 0.401 L/kg; 0.401 * 0.290 = 0.1163 L
    lq    <- log(0.02088); label("Shared inter-compartmental clearance Q at 290 g rat (L/h)")  # Jansson 2008 Table 2: Q = 1.2 mL/min/kg (shared); 1.2 * 60 / 1000 * 0.290 = 0.02088 L/h
    lvp   <- log(0.1293); label("Shared peripheral volume Vp at 290 g rat (L)")            # Jansson 2008 Table 2: Vp = 0.446 L/kg (shared); 0.446 * 0.290 = 0.1293 L

    # Oral absorption (Jansson 2008 Table 3). Shared Savic 2007 transit
    # chain delivers drug to per-enantiomer depots; the depot-to-central
    # transfer is Michaelis-Menten with Tmax (max absorption rate) and Kt
    # (amount-in-depot giving half-max rate). MW conversion uses eflornithine
    # MW 182.17 g/mol; Tmax (umol/min/kg) -> mg/h is value * MW * 60 / 1000;
    # then * 0.290 kg for a 290 g rat. Kt (umol/kg) -> mg is value * MW / 1000;
    # then * 0.290 kg for a 290 g rat.
    lmtt  <- log(1.467);  label("Shared mean transit time MTT (h)")                        # Jansson 2008 Table 3: MTT = 88 min = 88/60 = 1.467 h
    lnn   <- log(1.424);  label("Shared transit-compartment count n (continuous, dimensionless)")  # Jansson 2008 Table 3: n = 1.424 (Stirling-approximation Savic 2007 chain)

    ltmax_abs_l <- log(35.18); label("L-eflornithine saturable-absorption Vmax Tmax_L at 290 g rat (mg/h)")  # Jansson 2008 Table 3: Tmax_L = 11.1 umol/min/kg; 11.1 * 182.17e-3 * 60 * 0.290 = 35.18 mg/h
    ltmax_abs_d <- log(45.95); label("D-eflornithine saturable-absorption Vmax Tmax_D at 290 g rat (mg/h)")  # Jansson 2008 Table 3: Tmax_D = 14.5 umol/min/kg; 14.5 * 182.17e-3 * 60 * 0.290 = 45.95 mg/h
    lkt_abs_l   <- log(82.41); label("L-eflornithine saturable-absorption half-saturation amount Kt_L at 290 g rat (mg)")  # Jansson 2008 Table 3: Kt_L = 1560 umol/kg; 1560 * 182.17e-3 * 0.290 = 82.41 mg
    lkt_abs_d   <- log(41.42); label("D-eflornithine saturable-absorption half-saturation amount Kt_D at 290 g rat (mg)")  # Jansson 2008 Table 3: Kt_D = 784 umol/kg; 784 * 182.17e-3 * 0.290 = 41.42 mg

    # Bioavailability (Jansson 2008 Table 3). Per-enantiomer fdepot at the
    # 750-2000 mg/kg oral dose range; the highest-dose categorical effect
    # adds a multiplicative bump on F for both enantiomers.
    lfdepot_l <- log(0.41);  label("L-eflornithine oral bioavailability F_L at 750-2000 mg/kg (fraction)")  # Jansson 2008 Table 3: F_L = 41% at 750-2000 mg/kg
    lfdepot_d <- log(0.623); label("D-eflornithine oral bioavailability F_D at 750-2000 mg/kg (fraction)")  # Jansson 2008 Table 3: F_D = 62.3% at 750-2000 mg/kg
    e_dose_high_efl_fdepot_l <- 0.1463; label("Relative change in F_L at 3000 mg/kg vs 750-2000 mg/kg reference (fraction)")  # Jansson 2008 Table 3: F_L (3000) = 47% vs 41%; (47/41) - 1 = 0.1463
    e_dose_high_efl_fdepot_d <- 0.3275; label("Relative change in F_D at 3000 mg/kg vs 750-2000 mg/kg reference (fraction)")  # Jansson 2008 Table 3: F_D (3000) = 82.7% vs 62.3%; (82.7/62.3) - 1 = 0.3275

    # Inter-individual variability (omega^2 = log(1 + CV^2) for log-normal).
    # Per-enantiomer eta on CL (paper reports separately for L and D).
    # Shared eta on MTT and on F (paper: "Interindividual variability values
    # for mean transit time and bioavailability could not be estimated
    # separately for D- and L-eflornithine and were therefore assumed
    # to be identical"; Discussion paragraph + Table 3 footnote).
    etalcl_l   ~ 0.02197  # Jansson 2008 Table 2: IIV CL_L = 14.9% CV; log(1 + 0.149^2) = 0.02197
    etalcl_d   ~ 0.03263  # Jansson 2008 Table 2: IIV CL_D = 18.2% CV; log(1 + 0.182^2) = 0.03263
    etalmtt    ~ 0.1349   # Jansson 2008 Table 3: IIV MTT = 38% CV (shared L/D); log(1 + 0.38^2) = 0.1349
    etalfdepot ~ 0.00995  # Jansson 2008 Table 3: IIV F = 10% CV (shared L/D); log(1 + 0.10^2) = 0.00995

    # Residual error (proportional only, Jansson 2008 Methods + Tables 2-3).
    # Paper reports separate proportional residuals for the IV fit (sigma =
    # 19.9% CV) and the oral fit (sigma = 27.7% CV). The multi-output
    # encoding here applies the larger oral residual to all three outputs
    # (L, D, racemic) for simulation simplicity; the IV residual is
    # documented in the vignette Assumptions section.
    propSd_l   <- 0.277; label("Proportional residual error on L-eflornithine plasma Cc_l (fraction)")  # Jansson 2008 Table 3: oral sigma = 27.7% proportional
    propSd_d   <- 0.277; label("Proportional residual error on D-eflornithine plasma Cc_d (fraction)")  # Jansson 2008 Table 3: oral sigma = 27.7% (applied to D)
    propSd_rac <- 0.277; label("Proportional residual error on racemic plasma Cc_rac (fraction)")       # Jansson 2008 Table 3: oral sigma = 27.7% (applied to racemic sum)
  })

  model({
    # Shared Savic 2007 transit-compartment chain parameters.
    n   <- exp(lnn)
    mtt <- exp(lmtt + etalmtt)

    # Enantiomer-specific disposition (IV-derived) and saturable-absorption
    # (oral-derived) parameters.
    cl_l <- exp(lcl_l + etalcl_l)
    cl_d <- exp(lcl_d + etalcl_d)
    vc_l <- exp(lvc_l)
    vc_d <- exp(lvc_d)
    q    <- exp(lq)
    vp   <- exp(lvp)
    tmax_abs_l <- exp(ltmax_abs_l)
    tmax_abs_d <- exp(ltmax_abs_d)
    kt_abs_l   <- exp(lkt_abs_l)
    kt_abs_d   <- exp(lkt_abs_d)

    # Per-enantiomer bioavailability with categorical high-dose adjustment.
    # F enters via the bio argument of rxode2::transit(); the dose record
    # itself carries the 50/50 racemic split as separate cmt = depot_l and
    # cmt = depot_d dose rows so each transit chain reads its own podo() /
    # tad() (rxode2 transit() requires that the dose target compartment be
    # the same compartment whose d/dt() carries the transit() call).
    fdepot_l <- exp(lfdepot_l + etalfdepot) * (1 + e_dose_high_efl_fdepot_l * DOSE_HIGH_EFL)
    fdepot_d <- exp(lfdepot_d + etalfdepot) * (1 + e_dose_high_efl_fdepot_d * DOSE_HIGH_EFL)

    # Micro-rate constants for the two-compartment disposition; Q and Vp
    # are shared between L and D so k21_l = k21_d = q / vp.
    kel_l <- cl_l / vc_l
    kel_d <- cl_d / vc_d
    k12_l <- q / vc_l
    k12_d <- q / vc_d
    k21_l <- q / vp
    k21_d <- q / vp

    # ODE system. Transit chain feeds per-enantiomer depots via the rxode2
    # transit() built-in; depots drain to central via Michaelis-Menten
    # absorption; central + peripheral1 form the two-compartment
    # disposition. f(depot_*) <- 0 suppresses the bolus contribution so
    # the entire dose enters via the delayed Savic input function.
    d/dt(depot_l)       <- transit(n, mtt, fdepot_l) - tmax_abs_l * depot_l / (kt_abs_l + depot_l)
    d/dt(depot_d)       <- transit(n, mtt, fdepot_d) - tmax_abs_d * depot_d / (kt_abs_d + depot_d)
    d/dt(central_l)     <- tmax_abs_l * depot_l / (kt_abs_l + depot_l) - kel_l * central_l - k12_l * central_l + k21_l * peripheral1_l
    d/dt(central_d)     <- tmax_abs_d * depot_d / (kt_abs_d + depot_d) - kel_d * central_d - k12_d * central_d + k21_d * peripheral1_d
    d/dt(peripheral1_l) <- k12_l * central_l - k21_l * peripheral1_l
    d/dt(peripheral1_d) <- k12_d * central_d - k21_d * peripheral1_d

    f(depot_l) <- 0
    f(depot_d) <- 0

    # Plasma concentrations (dose mg / volume L = mg/L = ug/mL).
    Cc_l   <- central_l / vc_l
    Cc_d   <- central_d / vc_d
    Cc_rac <- Cc_l + Cc_d

    Cc_l   ~ prop(propSd_l)
    Cc_d   ~ prop(propSd_d)
    Cc_rac ~ prop(propSd_rac)
  })
}
