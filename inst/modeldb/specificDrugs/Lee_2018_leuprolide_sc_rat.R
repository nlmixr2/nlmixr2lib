Lee_2018_leuprolide_sc_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "One-compartment PK with first-order absorption for leuprolide solution",
    "after a single subcutaneous dose in normal male Wistar rats, linked to a",
    "testosterone-suppression PD model. Testosterone is a turnover pool whose",
    "zero-order production is STIMULATED by a sigmoid Emax function of plasma",
    "leuprolide (the flare-up effect) and simultaneously gated by a",
    "Gabrielsson-Hjorth moderator pool whose own production is inversely",
    "proportional to the prevailing testosterone level. Clearance and volume",
    "are apparent (CL/F, Vd/F) because only subcutaneous data informed them.",
    "Group 3 of Lee 2018; parameter values from Tables 1 and 4 and",
    "Equations 1-3 and 11-16."
  )
  reference <- paste(
    "Lee DS, Kim SJ, Choi GW, Lee YB, Cho HY.",
    "Pharmacokinetic-Pharmacodynamic Model for the Testosterone-Suppressive",
    "Effect of Leuprolide in Normal and Prostate Cancer Rats.",
    "Molecules. 2018;23(4):909. doi:10.3390/molecules23040909",
    sep = " "
  )
  vignette <- "Lee_2018_leuprolide"

  units <- list(time = "h", dosing = "ng", concentration = "ng/mL")

  dosing <- c("depot")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    ),
    TT = list(
      analyte = "testosterone",
      units = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    moderator1 = list(
      analyte = "hypothalamo-pituitary-gonadal feedback signal",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 5,
    n_studies = 1,
    disease_state = "normal (no prostate cancer)",
    dose_range = "single 0.1 mg/kg leuprolide acetate solution subcutaneously",
    regions = "Korea",
    notes = paste(
      "Group 3 of the six-group design in Table 5: adult male Wistar rats,",
      "n = 5, single subcutaneous 0.1 mg/kg leuprolide acetate solution.",
      "Sampling as for Group 1 (Methods 4.3). Absolute bioavailability against",
      "the intravenous Group 1 was 50.6% (Table 1), so CL and Vd from this",
      "group are apparent values; they are used unchanged here because",
      "Equations (2) and (3) divide by Vd with no F term. Body weight is NOT",
      "reported; the paper's own CL = Dose / AUC0-inf (Methods 4.6)",
      "back-solves the administered amount as 53.33 h*ng/mL x 514.46 mL/h =",
      "27,436 ng, implying a body weight of about 274 g and agreeing with the",
      "26,228 ng the intravenous group implies. The concurrent vehicle-only",
      "group (Group 2, n = 5) supplied the baseline testosterone turnover",
      "parameters kin and kout."
    )
  )

  ini({
    # ---------------- Leuprolide disposition (Table 1, SC / Group 3) --------
    lka <- log(16.67)
    label("First-order absorption rate constant ka (1/h)")
    # Table 1, column 'SC (Group 3)', row 'k a (h - 1)' = 16.67 +/- 2.55 (SE, n = 5)

    lcl <- log(514.46)
    label("Apparent clearance CL/F (mL/h)")
    # Table 1, column 'SC (Group 3)', row 'CL (mL/h)' = 514.46 +/- 40.10 (SE, n = 5)

    lvc <- log(487.40)
    label("Apparent central volume of distribution Vd/F (mL)")
    # Table 1, column 'SC (Group 3)', row 'V d (mL)' = 487.40 +/- 29.02 (SE, n = 5)

    # ---------------- Testosterone turnover (Table 4, Solution-SC) ----------
    lkin_tt <- log(0.68)
    label("Testosterone zero-order production rate constant kin (ng/mL/h)")
    # Table 4, column 'Solution-SC (Group 3)', row 'k in (h - 1)' = 0.68 (no SE; fixed from the Group 2 baseline fit)

    lkout_tt <- log(0.16)
    label("Testosterone first-order loss rate constant kout (1/h)")
    # Table 4, column 'Solution-SC (Group 3)', row 'k out (h - 1)' = 0.16 (no SE; fixed from the Group 2 baseline fit)

    lrbase_tt <- log(4.353)
    label("Baseline plasma testosterone R0 (ng/mL)")
    # Equation (12): R 0, wistar = 4.353. Results 2.3 reports the same quantity
    # as an observed mean basal concentration of 4.35 +/- 1.45 ng/mL (SE, n = 5).

    # ---------------- Drug effect on testosterone production (Table 4) ------
    lemax <- log(183.50)
    label("Maximum fractional increase in testosterone production Emax (unitless)")
    # Table 4, column 'Solution-SC (Group 3)', row 'E max' = 183.50 +/- 12.87 (SE, n = 5)

    lec50 <- log(6.17)
    label("Plasma leuprolide concentration at half-maximal effect EC50 (ng/mL)")
    # Table 4, column 'Solution-SC (Group 3)', row 'EC 50 (ng/mL)' = 6.17 +/- 2.49 (SE, n = 5)

    lhill <- log(2.02)
    label("Hill coefficient h of the leuprolide effect (unitless)")
    # Table 4, column 'Solution-SC (Group 3)', row 'h' = 2.02 +/- 0.05 (SE, n = 5)

    # ---------------- Feedback (moderator) pool (Table 4) -------------------
    # Equation (15): dF/dt = kf,on / R - kf,off * F. Table 4 heads both rows
    # 'h - 1'; that is right for kf,off but not for kf,on, which by Equation
    # (15) must carry ng/mL/h. See the vignette Errata.
    lkin_moderator1 <- log(0.14)
    label("Feedback-pool onset rate constant kf,on (ng/mL/h; Table 4 prints 1/h)")
    # Table 4, column 'Solution-SC (Group 3)', row 'k f,on (h - 1)' = 0.14 +/- 0.03 (SE, n = 5)

    lkout_moderator1 <- log(0.02)
    label("Feedback-pool offset rate constant kf,off (1/h)")
    # Table 4, column 'Solution-SC (Group 3)', row 'k f,off (h - 1)' = 0.02 +/- 0.014 (SE, n = 5)

    # ---------------- Residual unexplained variability ----------------------
    # Not reported: the paper fits each group's mean profile in Berkeley
    # Madonna by Runge-Kutta least squares (Methods 4.7) and every '+/-' in
    # its tables is the standard error of the estimate across n = 5 animals,
    # not a variance component. Encoded as fixed(0) per the standing policy
    # for unreported RUV; see vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual SD on plasma leuprolide (fraction; 0 -- not reported in the source)")

    propSd_TT <- fixed(0)
    label("Proportional residual SD on plasma testosterone (fraction; 0 -- not reported in the source)")
  })

  model({
    # ---------------- Individual parameters ---------------------------------
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kin_tt <- exp(lkin_tt)
    kout_tt <- exp(lkout_tt)
    rbase_tt <- exp(lrbase_tt)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kin_moderator1 <- exp(lkin_moderator1)
    kout_moderator1 <- exp(lkout_moderator1)

    kel <- cl / vc

    # ---------------- Leuprolide PK -----------------------------------------
    # Equation (1): d(Drug)/dt = dose - ka * Drug (the `dose` term is the
    # instantaneous input, carried by the event table rather than written into
    # the right-hand side).
    # Equation (2): d(Ap)/dt = ka * Drug - CL / Vd * Ap.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Equation (3): Cp = Ap / Vd. No bioavailability term appears anywhere in
    # Equations (1)-(3); the 50.6% absolute bioavailability of Table 1 is
    # absorbed into the apparent CL/F and Vd/F above, and applying it again
    # here would halve every predicted concentration.
    Cc <- central / vc

    # ---------------- Testosterone PD ---------------------------------------
    # Equation (16). max(0, Cc) guards the fractional power against a solver
    # undershoot into slightly negative concentrations.
    ce <- emax * max(0, Cc)^hill / (max(0, Cc)^hill + ec50^hill)

    # Equations (14) and (15).
    d/dt(TT) <- (1 + ce) * kin_tt * moderator1 - kout_tt * TT
    d/dt(moderator1) <- kin_moderator1 / TT - kout_moderator1 * moderator1

    # Equation (12) supplies TT(0). F(0) = 1 is the value at which Equation
    # (14) collapses onto the baseline model of Equation (11); see the
    # Lee_2018_leuprolide_iv_rat file for the full argument and the
    # sensitivity check.
    TT(0) <- rbase_tt
    moderator1(0) <- 1

    # ---------------- Observations ------------------------------------------
    Cc ~ prop(propSd)
    TT ~ prop(propSd_TT)
  })
}
