Lee_2018_leuprolide_iv_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "One-compartment PK of leuprolide solution after a single intravenous",
    "dose in normal male Wistar rats, linked to a testosterone-suppression",
    "PD model. Testosterone is a turnover pool whose zero-order production",
    "is STIMULATED by a sigmoid Emax function of plasma leuprolide (the",
    "flare-up effect) and simultaneously gated by a Gabrielsson-Hjorth",
    "moderator pool whose own production is inversely proportional to the",
    "prevailing testosterone level. The moderator is what turns the initial",
    "flare into the sustained suppression below baseline that follows it.",
    "Group 1 of Lee 2018; parameter values from Tables 1 and 4 and",
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

  # Amounts are carried in nanograms and volumes in millilitres, so that
  # `central / vc` is directly in ng/mL -- the unit of both the UPLC-MS/MS
  # leuprolide assay and the testosterone ELISA (Methods 4.4 and 4.5), and
  # the unit of every concentration in Tables 1 and 4.
  units <- list(time = "h", dosing = "ng", concentration = "ng/mL")

  dosing <- c("central")

  # The paper fits each experimental group separately and reports no
  # covariate effects; group membership (strain, route, formulation) is
  # carried by the choice of model file, not by a data column.
  covariateData <- list()

  compartmentData <- list(
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
    dose_range = "single 0.1 mg/kg leuprolide acetate solution intravenously",
    regions = "Korea",
    notes = paste(
      "Group 1 of the six-group design in Table 5: adult male Wistar rats,",
      "n = 5, single intravenous 0.1 mg/kg leuprolide acetate solution.",
      "Blood was drawn pre-dose and at 0.25, 1, 2, 4 and 8 h and 1, 2, 3, 4,",
      "5, 6, 7, 11 and 14 days (Methods 4.3). Body weight is NOT reported, so",
      "the administered amount cannot be computed from the mg/kg dose; the",
      "paper's own definition CL = Dose / AUC0-inf (Methods 4.6) back-solves",
      "it as 105.50 h*ng/mL x 248.61 mL/h = 26,228 ng, implying a body weight",
      "of about 262 g. The vignette uses that amount and shows that the",
      "subcutaneous group gives an independent and consistent estimate.",
      "The concurrent vehicle-only group (Group 2, n = 5) supplied the",
      "baseline testosterone turnover parameters kin and kout."
    )
  )

  ini({
    # ---------------- Leuprolide disposition (Table 1, IV / Group 1) --------
    # Table 1 reports k, CL and Vd as separate per-animal means, and the three
    # are not mutually consistent at the 10% level (CL / Vd = 1.288 1/h against
    # a reported k of 1.42 1/h and a t1/2 of 0.52 h, i.e. 1.33 1/h). Equation
    # (2) writes the elimination term as CL / Vd, so CL and Vd are the two
    # values transcribed here and k is left derived, as the paper's own
    # equation requires.
    lcl <- log(248.61)
    label("Clearance CL (mL/h)")
    # Table 1, column 'IV (Group 1)', row 'CL (mL/h)' = 248.61 +/- 35.07 (SE, n = 5)

    lvc <- log(192.95)
    label("Central volume of distribution Vd (mL)")
    # Table 1, column 'IV (Group 1)', row 'V d (mL)' = 192.95 +/- 54.24 (SE, n = 5)

    # ---------------- Testosterone turnover (Table 4, Solution-IV) ----------
    # Equation (11), the baseline model fitted to the vehicle-only Wistar
    # group (Group 2): dR/dt = kin - kout * R. Table 4 prints kin and kout
    # without an SE because they were fixed from that baseline fit and reused
    # for every Wistar drug-effect fit ("These were estimated by curve fitting
    # with the vehicle-administered group and used as the basic parameters of
    # the drug-effect model described below (Group 2)").
    lkin_tt <- log(0.68)
    label("Testosterone zero-order production rate constant kin (ng/mL/h)")
    # Table 4, column 'Solution-IV (Group 1)', row 'k in (h - 1)' = 0.68 (no SE; fixed from the Group 2 baseline fit)

    lkout_tt <- log(0.16)
    label("Testosterone first-order loss rate constant kout (1/h)")
    # Table 4, column 'Solution-IV (Group 1)', row 'k out (h - 1)' = 0.16 (no SE; fixed from the Group 2 baseline fit)

    lrbase_tt <- log(4.353)
    label("Baseline plasma testosterone R0 (ng/mL)")
    # Equation (12): R 0, wistar = 4.353. Results 2.3 reports the same quantity
    # as an observed mean basal concentration of 4.35 +/- 1.45 ng/mL (SE, n = 5).

    # ---------------- Drug effect on testosterone production (Table 4) ------
    # Equation (16): CE = Emax * C^h / (C^h + EC50^h). Equation (14) enters it
    # as (1 + CE), so Emax is a dimensionless fractional INCREMENT on the
    # production rate, not an attained level -- at full effect production is
    # (1 + Emax) times basal. It is large because the flare-up raises plasma
    # testosterone roughly two orders of magnitude above baseline within hours.
    lemax <- log(303.77)
    label("Maximum fractional increase in testosterone production Emax (unitless)")
    # Table 4, column 'Solution-IV (Group 1)', row 'E max' = 303.77 +/- 12.90 (SE, n = 5)

    lec50 <- log(3.48)
    label("Plasma leuprolide concentration at half-maximal effect EC50 (ng/mL)")
    # Table 4, column 'Solution-IV (Group 1)', row 'EC 50 (ng/mL)' = 3.48 +/- 1.74 (SE, n = 5)

    lhill <- log(2.00)
    label("Hill coefficient h of the leuprolide effect (unitless)")
    # Table 4, column 'Solution-IV (Group 1)', row 'h' = 2.00 +/- 0.61 (SE, n = 5)

    # ---------------- Feedback (moderator) pool (Table 4) -------------------
    # Equation (15): dF/dt = kf,on / R - kf,off * F. The pool is a
    # Gabrielsson-Hjorth moderator: it carries no mass, it multiplies the
    # testosterone production rate in Equation (14), and its own production is
    # INVERSELY proportional to the prevailing testosterone level, which is
    # what makes the loop negative. Table 4 heads both rows 'h - 1'; that is
    # right for kf,off but not for kf,on, which by Equation (15) must carry
    # ng/mL/h for dF/dt to be a pure 1/h rate. The label states the units the
    # equation requires. See the vignette Errata.
    lkin_moderator1 <- log(0.29)
    label("Feedback-pool onset rate constant kf,on (ng/mL/h; Table 4 prints 1/h)")
    # Table 4, column 'Solution-IV (Group 1)', row 'k f,on (h - 1)' = 0.29 +/- 0.19 (SE, n = 5)

    lkout_moderator1 <- log(0.059)
    label("Feedback-pool offset rate constant kf,off (1/h)")
    # Table 4, column 'Solution-IV (Group 1)', row 'k f,off (h - 1)' = 0.059 +/- 0.017 (SE, n = 5)

    # ---------------- Residual unexplained variability ----------------------
    # The paper fits each group's mean profile in Berkeley Madonna by
    # Runge-Kutta least squares (Methods 4.7) and reports neither a residual
    # error model nor between-animal variability -- the '+/-' in every table is
    # the standard error of the parameter estimate across n = 5 animals, not a
    # variance component. Encoded as fixed(0) per the standing policy for
    # unreported RUV rather than inventing a magnitude; see vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual SD on plasma leuprolide (fraction; 0 -- not reported in the source)")

    propSd_TT <- fixed(0)
    label("Proportional residual SD on plasma testosterone (fraction; 0 -- not reported in the source)")
  })

  model({
    # ---------------- Individual parameters ---------------------------------
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
    # Equation (2) with the absorption term dropped for the intravenous route
    # (Equation (1), the `Drug` depot, applies only to the SC solution group):
    #   d(Ap)/dt = -CL / Vd * Ap
    d/dt(central) <- -kel * central

    # Equation (3): Cp = Ap / Vd. The paper writes `C` in Equation (16).
    Cc <- central / vc

    # ---------------- Testosterone PD ---------------------------------------
    # Equation (16). max(0, Cc) guards the fractional power against a solver
    # undershoot into slightly negative concentrations; Cc is a concentration
    # and the guard is inactive for every physically meaningful value.
    ce <- emax * max(0, Cc)^hill / (max(0, Cc)^hill + ec50^hill)

    # Equation (14): dR/dt = (1 + CE) * kin * F - kout * R, and
    # Equation (15): dF/dt = kf,on / R - kf,off * F.
    d/dt(TT) <- (1 + ce) * kin_tt * moderator1 - kout_tt * TT
    d/dt(moderator1) <- kin_moderator1 / TT - kout_moderator1 * moderator1

    # Equation (12) supplies the testosterone initial condition. The paper does
    # not state F(0); it is set to 1 because that is the value at which
    # Equation (14) collapses exactly onto the baseline model of Equation (11),
    # which is how the paper describes the drug-effect model being built ("the
    # drug-effect model was established based on kin and kout on the baseline
    # model"). The choice is not load-bearing: starting F instead at its own
    # quasi-steady state kf,on / (kf,off * R0) moves the flare peak by about
    # 13% and leaves the nadir and the recovered plateau unchanged, because the
    # coupled system is attracted to the same joint steady state either way.
    TT(0) <- rbase_tt
    moderator1(0) <- 1

    # ---------------- Observations ------------------------------------------
    Cc ~ prop(propSd)
    TT ~ prop(propSd_TT)
  })
}
