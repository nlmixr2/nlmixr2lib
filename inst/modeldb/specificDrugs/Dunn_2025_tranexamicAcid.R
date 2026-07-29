Dunn_2025_tranexamicAcid <- function() {
  description <- "Two-compartment population PK model for tranexamic acid (TXA) with parallel first-order intramuscular and first-order oral absorption (oral lag time) and first-order elimination, in pregnant individuals receiving IV, IM, or oral TXA for prevention or treatment of postpartum hemorrhage (Dunn 2025)."
  reference   <- "Dunn A, Felfeli M, Seifert SM, Gilliot S, Ducloy-Bouthors A-S, Shakur-Still H, Geer A, Grassin-Delyle S, Luban NL, van den Anker JN, Gobburu JVS, Roberts I, Ahmadzia HK. Evaluating tranexamic acid dosing strategies for postpartum hemorrhage: a population pharmacokinetic approach in pregnant individuals. J Clin Pharmacol. 2025;65(10):1262-1272. doi:10.1002/jcph.70031"
  vignette    <- "Dunn_2025_tranexamicAcid"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling with reference 80 kg (median pooled-cohort body weight, Dunn 2025 Table 1). Shared exponent 0.89 applied to CL and Q; shared exponent 0.44 applied to Vc and Vp (Dunn 2025 Results / Table 2). All pregnant participants; pregnancy status is a population-level fact (no PREG covariate effect estimated in the source).",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 221L,
    n_studies       = 4L,
    age_range       = "22-47 years (mean 33)",
    age_median      = "33 years",
    weight_range    = "47-156 kg (mean 81.7)",
    weight_median   = "78 kg",
    bmi_range       = "18.0-55.8 kg/m^2 (mean 31.0)",
    bmi_median      = "30.0 kg/m^2",
    sex_female_pct  = 100,
    disease_state   = "Pregnant or immediately postpartum individuals undergoing or at risk of postpartum hemorrhage during caesarean delivery",
    dose_range      = "IV: fixed 0.5 g (n=34), 1 g (n=97), 1.5 g (n=1), 2 g (n=3), or weight-based 5/10/15 mg/kg (n=10 each); IM: fixed 1 g (n=26); oral: fixed 4 g (n=30). Dunn 2025 Table 1.",
    regions         = "Pooled across four trials: NCT03863964 (USA), NCT02797119 / TRACES (France), NCT04274335 / WOMAN-PharmacoTXA (multi-national), NCT03287336 (USA). Dunn 2025 Methods.",
    n_observations  = 1303L,
    notes           = "Pooled population PK dataset from four trials in pregnant participants receiving TXA at the time of caesarean delivery. Pregnancy status is fixed at 100% (all subjects pregnant); see Dunn 2025 Table 1 for the per-study demographic breakdown. Routes of administration carried in the dosing dataset via the rxode2 cmt column: 'central' for IV (bolus or infusion), 'depot_im' for IM, 'depot_oral' for oral. IM bioavailability was held at 1.0 because the typical estimate approached unity (Dunn 2025 Results, Structural PK Model)."
  )

  ini({
    # Structural parameters from Dunn 2025 Table 2 Final PK Model. All clearances
    # and volumes are reported for a typical participant of 80 kg actual body weight.
    lcl          <- log(8.59);  label("Clearance for the reference 80 kg subject (CL, L/h)")            # Dunn 2025 Table 2
    lvc          <- log(10.7);  label("Central volume of distribution for the reference 80 kg subject (Vc, L)")  # Dunn 2025 Table 2
    lq           <- log(28.9);  label("Inter-compartmental clearance for the reference 80 kg subject (Q, L/h)")  # Dunn 2025 Table 2
    lvp          <- log(15.0);  label("Peripheral volume of distribution for the reference 80 kg subject (Vp, L)") # Dunn 2025 Table 2
    lka_oral     <- log(0.18);  label("First-order oral absorption rate constant (Ka_oral, 1/h)")       # Dunn 2025 Table 2
    lka_im       <- log(2.31);  label("First-order intramuscular absorption rate constant (Ka_IM, 1/h)") # Dunn 2025 Table 2
    lfdepot_oral <- log(0.56);  label("Oral bioavailability (F_oral, fraction)")                        # Dunn 2025 Table 2
    llag_oral    <- log(0.16);  label("Oral absorption lag time (Tlag_oral, h)")                        # Dunn 2025 Table 2

    # Allometric body-weight scaling. Exponents were estimated (not fixed); pooled
    # across the clearance and volume categories per Dunn 2025 Results: "an exponent
    # of 0.89 was estimated [for CL and Q]" and "a final pooled estimate of 0.44 was
    # retained [for Vc and Vp]." Form: P_i = P_pop * (WT/80)^exp.
    e_wt_cl_q  <- 0.89; label("Shared allometric exponent of (WT/80) on CL and Q (unitless)")  # Dunn 2025 Results / Table 2
    e_wt_vc_vp <- 0.44; label("Shared allometric exponent of (WT/80) on Vc and Vp (unitless)") # Dunn 2025 Results / Table 2

    # Inter-individual variability. Dunn 2025 Table 2 reports BSV as %CV; the
    # log-normal variance is omega^2 = log(1 + CV^2). The paper states (Methods)
    # that BSV was NOT imposed on Q, Ka_oral, or Ka_IM due to high shrinkage,
    # so etas exist only for CL, Vc, Vp, F_oral, and Tlag_oral.
    etalcl          ~ log(1 + 0.324^2)  # Dunn 2025 Table 2 BSV CL = 32.4 %CV
    etalvc          ~ log(1 + 0.595^2)  # Dunn 2025 Table 2 BSV Vc = 59.5 %CV
    etalvp          ~ log(1 + 0.464^2)  # Dunn 2025 Table 2 BSV Vp = 46.4 %CV
    etalfdepot_oral ~ log(1 + 0.448^2)  # Dunn 2025 Table 2 BSV F_oral = 44.8 %CV (shrinkage 67%; retained per paper)
    etallag_oral    ~ log(1 + 0.646^2)  # Dunn 2025 Table 2 BSV Tlag_oral = 64.6 %CV (shrinkage 69%; retained per paper)

    # Residual error: combined additive + proportional. Dunn 2025 Table 2:
    #   sigma_additive     = 0.67 mg/L
    #   sigma_proportional = 27.2 % (i.e., fractional SD 0.272)
    addSd  <- 0.67;  label("Additive residual error (mg/L)")             # Dunn 2025 Table 2
    propSd <- 0.272; label("Proportional residual error (fraction)")     # Dunn 2025 Table 2
  })

  model({
    # Individual PK parameters at reference WT = 80 kg, with allometric
    # body-weight scaling per Dunn 2025 (shared exponents within
    # clearance and volume categories).
    cl       <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl_q
    vc       <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc_vp
    q        <- exp(lq)           * (WT / 80)^e_wt_cl_q
    vp       <- exp(lvp + etalvp) * (WT / 80)^e_wt_vc_vp
    ka_oral  <- exp(lka_oral)
    ka_im    <- exp(lka_im)
    foral    <- exp(lfdepot_oral + etalfdepot_oral)
    tlagoral <- exp(llag_oral    + etallag_oral)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Parallel absorption: route is selected by the cmt column on each dose row.
    # IV doses go directly to central (cmt = "central"); IM doses go to depot_im
    # (cmt = "depot_im"); oral doses go to depot_oral (cmt = "depot_oral").
    d/dt(depot_im)    <- -ka_im   * depot_im
    d/dt(depot_oral)  <- -ka_oral * depot_oral
    d/dt(central)     <-  ka_im   * depot_im +
                          ka_oral * depot_oral -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Route-specific bioavailability and lag time. IM bioavailability was fixed
    # at 1.0 in Dunn 2025 (typical estimate approached unity; see Results,
    # Structural PK Model), so no f(depot_im) adjustment is needed.
    f(depot_oral)   <- foral
    lag(depot_oral) <- tlagoral

    # Concentration: dose in mg, vc in L -> mg/L (== ug/mL), matching the
    # Dunn 2025 Table 2 concentration units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
