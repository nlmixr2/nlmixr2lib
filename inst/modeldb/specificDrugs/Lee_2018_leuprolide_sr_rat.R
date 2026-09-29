Lee_2018_leuprolide_sr_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Three-section sustained-release depot PK for the leuprolide acetate",
    "microsphere depot (Lucrin depot) after a single subcutaneous dose in",
    "normal male Wistar rats, linked to a testosterone-suppression PD model.",
    "The dose is split across a non-capsuled section released at the solution",
    "absorption rate ka, a diffusive section released at kd after a lag, and",
    "an erosive section released through a two-compartment transit chain at",
    "kt after a longer lag. Testosterone is a turnover pool whose zero-order",
    "production is STIMULATED by a sigmoid Emax function of plasma leuprolide",
    "(the flare-up effect) and simultaneously gated by a Gabrielsson-Hjorth",
    "moderator pool whose own production is inversely proportional to the",
    "prevailing testosterone level. Group 4 of Lee 2018; parameter values",
    "from Tables 1, 2 and 4 and Equations 4-10 and 11-16."
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

  # Equations (4), (5) and (6) each receive their own share of the dose, so
  # the event table carries three simultaneous dose records -- one per release
  # section -- and the f() fractions below split them.
  dosing <- c("depot", "depot2", "transit1")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
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
    dose_range = "single Lucrin depot dispersion equivalent to 0.1 mg/kg leuprolide acetate subcutaneously",
    regions = "Korea",
    notes = paste(
      "Group 4 of the six-group design in Table 5: adult male Wistar rats,",
      "n = 5, single subcutaneous Lucrin depot (leuprolide acetate 3.75 mg",
      "microspheres) dispersion equivalent to 0.1 mg/kg leuprolide acetate.",
      "Sampling as for Group 1 (Methods 4.3). The disposition parameters ka,",
      "CL/F and Vd/F are NOT re-estimated here: Methods 4.7 states that the",
      "depot model was built 'based on the PK parameters such as ka, CL, and",
      "Vd in solution-administered groups', so they are carried unchanged",
      "from the subcutaneous solution group (Group 3, Table 1). Only the",
      "release parameters of Table 2 and the PD parameters of Table 4 are",
      "specific to this group. Body weight is NOT reported; the vignette uses",
      "the 27,436 ng the subcutaneous solution group back-solves from",
      "CL = Dose / AUC0-inf. The concurrent vehicle-only group (Group 2,",
      "n = 5) supplied the baseline testosterone turnover parameters."
    )
  )

  ini({
    # ---------------- Leuprolide disposition (Table 1, SC / Group 3) --------
    # Carried unchanged from the solution-administered group per Methods 4.7;
    # apparent values, since only subcutaneous data informed them.
    lka <- log(16.67)
    label("First-order absorption rate constant ka out of the non-capsuled section (1/h)")
    # Table 1, column 'SC (Group 3)', row 'k a (h - 1)' = 16.67 +/- 2.55 (SE, n = 5)

    lcl <- log(514.46)
    label("Apparent clearance CL/F (mL/h)")
    # Table 1, column 'SC (Group 3)', row 'CL (mL/h)' = 514.46 +/- 40.10 (SE, n = 5)

    lvc <- log(487.40)
    label("Apparent central volume of distribution Vd/F (mL)")
    # Table 1, column 'SC (Group 3)', row 'V d (mL)' = 487.40 +/- 29.02 (SE, n = 5)

    # ---------------- Depot release (Table 2, Wistar / Group 4) -------------
    # Results 2.4: 'the kd and kt used in this model are rate constants that
    # include both the release process in the biological systems and the
    # absorption process of the released drug, unlike the ka which contains
    # only the absorption constant of the molecular state.' They are therefore
    # apparent first-order input rate constants, which is exactly what the
    # canonical ka2 (second parallel-absorption depot) and ktr (transit-chain
    # rate) denote.
    lka2 <- log(0.08)
    label("Diffusive release rate constant kd out of the diffusion section (1/h)")
    # Table 2, column 'Wistar Rats (Group 4)', row 'k d (h - 1)' = 0.08 +/- 0.01 (SE, n = 5)

    lktr <- log(0.0078)
    label("Erosive release rate constant kt of the two-stage erosion transit chain (1/h)")
    # Table 2, column 'Wistar Rats (Group 4)', row 'k t (h - 1)' = 0.0078 +/- 0.0002 (SE, n = 5)

    # Equation (7) constrains NR + DR + ER = 1, so the three fractions live on
    # a 2-simplex and only two quantities are free. They are encoded by
    # stick-breaking -- the first logit is NR itself, the second is DR
    # CONDITIONAL on the dose not being in the non-capsuled section -- which
    # reproduces the printed 0.18 / 0.28 / 0.54 exactly and cannot leak outside
    # the simplex under any downstream perturbation. The register's `logitfrel`
    # entry anticipates this ('Three or more processes would take `logitfrel2`
    # etc., but a softmax / stick-breaking encoding should be considered at
    # that point').
    logitfrel <- log(0.18 / (1 - 0.18))
    label("Non-capsuled release fraction NR, logit scale (unitless)")
    # Table 2, column 'Wistar Rats (Group 4)', row 'N R' = 0.18 +/- 0.04 (SE, n = 5)

    logitfrel2 <- log((0.28 / (1 - 0.18)) / (1 - 0.28 / (1 - 0.18)))
    label("Diffusive release fraction DR conditional on not being non-capsuled, logit scale (unitless)")
    # Table 2, column 'Wistar Rats (Group 4)', row 'D R' = 0.28 +/- 0.07 (SE, n = 5);
    # conditional value 0.28 / (1 - 0.18) = 0.34146, so that the erosive
    # remainder is (1 - 0.18) * (1 - 0.34146) = 0.54 = Table 2 row 'E R' = 0.54 +/- 0.08

    ltlag <- log(0.47)
    label("Lag time before diffusive release begins t_lag,d (h)")
    # Table 2, column 'Wistar Rats (Group 4)', row 't lag,d (h)' = 0.47 +/- 0.09 (SE, n = 5)

    ltlag2 <- log(3.61)
    label("Lag time before erosive release begins t_lag,e (h)")
    # Table 2, column 'Wistar Rats (Group 4)', row 't lag ,e (h)' = 3.61 +/- 1.53 (SE, n = 5)

    # ---------------- Testosterone turnover (Table 4, SR-SC Group 4) --------
    lkin_tt <- log(0.68)
    label("Testosterone zero-order production rate constant kin (ng/mL/h)")
    # Table 4, column 'SR-SC (Group 4)', row 'k in (h - 1)' = 0.68 (no SE; fixed from the Group 2 baseline fit)

    lkout_tt <- log(0.16)
    label("Testosterone first-order loss rate constant kout (1/h)")
    # Table 4, column 'SR-SC (Group 4)', row 'k out (h - 1)' = 0.16 (no SE; fixed from the Group 2 baseline fit)

    lrbase_tt <- log(4.353)
    label("Baseline plasma testosterone R0 (ng/mL)")
    # Equation (12): R 0, wistar = 4.353. Results 2.3 reports the same quantity
    # as an observed mean basal concentration of 4.35 +/- 1.45 ng/mL (SE, n = 5).

    # ---------------- Drug effect on testosterone production (Table 4) ------
    lemax <- log(380.00)
    label("Maximum fractional increase in testosterone production Emax (unitless)")
    # Table 4, column 'SR-SC (Group 4)', row 'E max' = 380.00 +/- 87.54 (SE, n = 5)

    lec50 <- log(1.80)
    label("Plasma leuprolide concentration at half-maximal effect EC50 (ng/mL)")
    # Table 4, column 'SR-SC (Group 4)', row 'EC 50 (ng/mL)' = 1.80 +/- 0.57 (SE, n = 5)

    lhill <- log(2.00)
    label("Hill coefficient h of the leuprolide effect (unitless)")
    # Table 4, column 'SR-SC (Group 4)', row 'h' = 2.00 +/- 0.01 (SE, n = 5)

    # ---------------- Feedback (moderator) pool (Table 4) -------------------
    # Equation (15): dF/dt = kf,on / R - kf,off * F. Table 4 heads both rows
    # 'h - 1'; that is right for kf,off but not for kf,on, which by Equation
    # (15) must carry ng/mL/h. See the vignette Errata.
    lkin_moderator1 <- log(0.40)
    label("Feedback-pool onset rate constant kf,on (ng/mL/h; Table 4 prints 1/h)")
    # Table 4, column 'SR-SC (Group 4)', row 'k f,on (h - 1)' = 0.40 +/- 0.25 (SE, n = 5)

    lkout_moderator1 <- log(0.04)
    label("Feedback-pool offset rate constant kf,off (1/h)")
    # Table 4, column 'SR-SC (Group 4)', row 'k f,off (h - 1)' = 0.04 +/- 0.032 (SE, n = 5)

    # ---------------- Residual unexplained variability ----------------------
    # Not reported; see vignette Errata and the standing policy for unreported
    # RUV. Every '+/-' in the source tables is the standard error of the
    # estimate across n = 5 animals, not a variance component.
    propSd <- fixed(0)
    label("Proportional residual SD on plasma leuprolide (fraction; 0 -- not reported in the source)")

    propSd_TT <- fixed(0)
    label("Proportional residual SD on plasma testosterone (fraction; 0 -- not reported in the source)")
  })

  model({
    # ---------------- Individual parameters ---------------------------------
    ka <- exp(lka)
    ka2 <- exp(lka2)
    ktr <- exp(lktr)
    cl <- exp(lcl)
    vc <- exp(lvc)
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)
    kin_tt <- exp(lkin_tt)
    kout_tt <- exp(lkout_tt)
    rbase_tt <- exp(lrbase_tt)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kin_moderator1 <- exp(lkin_moderator1)
    kout_moderator1 <- exp(lkout_moderator1)

    kel <- cl / vc

    # ---------------- Release-section fractions, Equation (7) ---------------
    # Stick-breaking back-transform: frel is NR, frel2 is DR and frel3 is ER,
    # and frel + frel2 + frel3 == 1 identically.
    frel <- expit(logitfrel)
    frel2 <- (1 - frel) * expit(logitfrel2)
    frel3 <- (1 - frel) * (1 - expit(logitfrel2))

    # ---------------- Sustained-release depot, Equations (4)-(9) ------------
    # Equation (4): d(NS)/dt = dose * NR - ka * NS        -> depot
    # Equation (5): d(DS)/dt = dose * DR - kd * DS        -> depot2
    # Equation (6): d(ES1)/dt = dose * ER - kt * ES1      -> transit1
    # Equation (8): d(ES_2..n)/dt = kt * ES_(n-1) - kt * ES_n -> transit2
    # Table 2 row 'ES n' = 2 and Results 2.4: 'We have attempted to model 2 to
    # 10 transit compartments (ES2-ES10) and obtained the best fit in 2
    # transit compartments; therefore, the number of transit models was fixed
    # at 2.' The chain length is therefore structural, not a parameter.
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2

    # Equation (9): d(Ap)/dt = ka*NS + kd*DS + kt*ES_n - CL/Vd * Ap
    d/dt(central) <- ka * depot + ka2 * depot2 + ktr * transit2 - kel * central

    # The `dose * NR`, `dose * DR` and `dose * ER` inputs of Equations (4)-(6)
    # are the instantaneous dose split; each section receives its share of the
    # same administered amount through its own dose record.
    f(depot) <- frel
    f(depot2) <- frel2
    f(transit1) <- frel3

    # Results 2.4: 'The Ap is absorbed from each section, NS, DS, and ES, at
    # various release times (the lag time of drug release in DS, t lag,d; the
    # lag time of drug release in ES, t lag,e).' The non-capsuled section is
    # the free, un-entrapped drug and carries no lag.
    lag(depot2) <- tlag
    lag(transit1) <- tlag2

    # Equation (10): C = Ap / Vd.
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
    # Lee_2018_leuprolide_iv_rat file for the full argument.
    TT(0) <- rbase_tt
    moderator1(0) <- 1

    # ---------------- Observations ------------------------------------------
    Cc ~ prop(propSd)
    TT ~ prop(propSd_TT)
  })
}
