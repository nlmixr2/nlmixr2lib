Lee_2018_leuprolide_sr_cancer_rat <- function() {
  description <- paste(
    "Preclinical (rat, prostate cancer).",
    "Three-section sustained-release depot PK for the leuprolide acetate",
    "microsphere depot (Lucrin depot) after a single subcutaneous dose in male",
    "Iar:COP rats bearing spontaneous Dunning R-3327 prostate adenocarcinoma,",
    "linked to a testosterone-suppression PD model. The dose is split across a",
    "non-capsuled section released at the solution absorption rate ka, a",
    "diffusive section released at kd after a lag, and an erosive section",
    "released through a two-compartment transit chain at kt after a longer",
    "lag. Testosterone is a turnover pool whose zero-order production is",
    "STIMULATED by a sigmoid Emax function of plasma leuprolide (the flare-up",
    "effect) and simultaneously gated by a Gabrielsson-Hjorth moderator pool",
    "whose own production is inversely proportional to the prevailing",
    "testosterone level. Group 6 of Lee 2018; parameter values from Tables 1,",
    "2 and 4 and Equations 4-10 and 11-16. NOTE: this group's published",
    "kf,on / kf,off pair does not reproduce the paper's own Figure 2D plateau",
    "-- see the vignette Errata before using it."
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
    species = "rat (Iar:COP, Copenhagen)",
    n_subjects = 5,
    n_studies = 1,
    disease_state = "spontaneous Dunning R-3327 prostate adenocarcinoma",
    dose_range = "single Lucrin depot dispersion equivalent to 0.1 mg/kg leuprolide acetate subcutaneously",
    regions = "Korea",
    notes = paste(
      "Group 6 of the six-group design in Table 5: adult male Iar:COP",
      "(Copenhagen) rats, n = 5, single subcutaneous Lucrin depot dispersion",
      "equivalent to 0.1 mg/kg leuprolide acetate. Methods 4.3: 'Dunning",
      "observed the appearance of prostate cancer in male Iar:COP rats, and",
      "Dunning R-3327 adenocarcinoma is a spontaneously occurred prostate",
      "tumor found in male Iar:COP rats.' Sampling as for Group 1. The",
      "disposition parameters ka, CL/F and Vd/F are carried unchanged from the",
      "Wistar subcutaneous solution group (Group 3, Table 1) per Methods 4.7;",
      "only the release parameters of Table 2 and the PD parameters of Table 4",
      "are specific to this group. The concurrent vehicle-only Iar:COP group",
      "(Group 5, n = 5) supplied the baseline testosterone turnover",
      "parameters. Results 2.4 notes that ka, CL and Vd did not differ",
      "significantly between strains, whereas Emax and EC50 did (p < 0.05)."
    )
  )

  ini({
    # ---------------- Leuprolide disposition (Table 1, SC / Group 3) --------
    # Carried unchanged from the Wistar solution group per Methods 4.7; the
    # paper does not re-estimate disposition in the Iar:COP strain.
    lka <- log(16.67)
    label("First-order absorption rate constant ka out of the non-capsuled section (1/h)")
    # Table 1, column 'SC (Group 3)', row 'k a (h - 1)' = 16.67 +/- 2.55 (SE, n = 5)

    lcl <- log(514.46)
    label("Apparent clearance CL/F (mL/h)")
    # Table 1, column 'SC (Group 3)', row 'CL (mL/h)' = 514.46 +/- 40.10 (SE, n = 5)

    lvc <- log(487.40)
    label("Apparent central volume of distribution Vd/F (mL)")
    # Table 1, column 'SC (Group 3)', row 'V d (mL)' = 487.40 +/- 29.02 (SE, n = 5)

    # ---------------- Depot release (Table 2, Iar:COP / Group 6) ------------
    lka2 <- log(0.08)
    label("Diffusive release rate constant kd out of the diffusion section (1/h)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 'k d (h - 1)' = 0.08 +/- 0.01 (SE, n = 5)

    lktr <- log(0.0193)
    label("Erosive release rate constant kt of the two-stage erosion transit chain (1/h)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 'k t (h - 1)' = 0.0193 +/- 0.0083 (SE, n = 5).
    # Results 2.4 explains why this one differs between strains while kd does
    # not: 'since the release in the erosive section proceeds very slowly, the
    # formulations are greatly influenced by physiological properties such as
    # pH during release over several weeks'.

    # Equation (7) constrains NR + DR + ER = 1; stick-breaking encoding, as in
    # the Wistar sibling file.
    logitfrel <- log(0.08 / (1 - 0.08))
    label("Non-capsuled release fraction NR, logit scale (unitless)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 'N R' = 0.08 +/- 0.02 (SE, n = 5)

    logitfrel2 <- log((0.43 / (1 - 0.08)) / (1 - 0.43 / (1 - 0.08)))
    label("Diffusive release fraction DR conditional on not being non-capsuled, logit scale (unitless)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 'D R' = 0.43 +/- 0.08 (SE, n = 5);
    # conditional value 0.43 / (1 - 0.08) = 0.46739, so that the erosive
    # remainder is (1 - 0.08) * (1 - 0.46739) = 0.49 = Table 2 row 'E R' = 0.49 +/- 0.08

    ltlag <- log(0.35)
    label("Lag time before diffusive release begins t_lag,d (h)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 't lag,d (h)' = 0.35 +/- 0.06 (SE, n = 5)

    ltlag2 <- log(2.58)
    label("Lag time before erosive release begins t_lag,e (h)")
    # Table 2, column 'Iar:COP Rats (Group 6)', row 't lag ,e (h)' = 2.58 +/- 0.22 (SE, n = 5)

    # ---------------- Testosterone turnover (Table 4, SR-SC Group 6) --------
    lkin_tt <- log(0.35)
    label("Testosterone zero-order production rate constant kin (ng/mL/h)")
    # Table 4, column 'SR-SC (Group 6)', row 'k in (h - 1)' = 0.35 (no SE; fixed from the Group 5 baseline fit)

    lkout_tt <- log(0.06)
    label("Testosterone first-order loss rate constant kout (1/h)")
    # Table 4, column 'SR-SC (Group 6)', row 'k out (h - 1)' = 0.06 (no SE; fixed from the Group 5 baseline fit)

    lrbase_tt <- log(4.094)
    label("Baseline plasma testosterone R0 (ng/mL)")
    # Equation (13): R 0, Iar:COP = 4.094. Results 2.3 reports the same quantity
    # as an observed mean basal concentration of 4.09 +/- 0.97 ng/mL (SE, n = 5).
    # Note that kin / kout = 0.35 / 0.06 = 5.83 ng/mL is NOT equal to this
    # value, so the baseline model of Equation (11) is not at rest at R0 for
    # this group; see the vignette Errata.

    # ---------------- Drug effect on testosterone production (Table 4) ------
    lemax <- log(634.50)
    label("Maximum fractional increase in testosterone production Emax (unitless)")
    # Table 4, column 'SR-SC (Group 6)', row 'E max' = 634.50 +/- 144.73 (SE, n = 5); * p < 0.05 versus Wistar SR-SC

    lec50 <- log(3.34)
    label("Plasma leuprolide concentration at half-maximal effect EC50 (ng/mL)")
    # Table 4, column 'SR-SC (Group 6)', row 'EC 50 (ng/mL)' = 3.34 +/- 0.56 (SE, n = 5); * p < 0.05 versus Wistar SR-SC

    lhill <- log(3.18)
    label("Hill coefficient h of the leuprolide effect (unitless)")
    # Table 4, column 'SR-SC (Group 6)', row 'h' = 3.18 +/- 1.27 (SE, n = 5)

    # ---------------- Feedback (moderator) pool (Table 4) -------------------
    # Equation (15): dF/dt = kf,on / R - kf,off * F. This group is the one
    # place where the published numbers do not reproduce the paper's own
    # figure: the joint steady state of Equations (14) and (15) is
    # sqrt(kin * kf,on / (kout * kf,off)), which these values put at 1.04
    # ng/mL against the roughly 4.2 ng/mL plateau drawn in Figure 2D and the
    # 4.39 ng/mL observed mean of Results 2.3. The same closed form reproduces
    # all three Wistar groups. The printed values are transcribed unchanged --
    # nothing is tuned -- and the discrepancy is documented in the vignette
    # Errata. Note also that this is the only group whose kf,on is smaller
    # than its kf,off and whose kf,on SE exceeds its point estimate.
    lkin_moderator1 <- log(0.083)
    label("Feedback-pool onset rate constant kf,on (ng/mL/h; Table 4 prints 1/h)")
    # Table 4, column 'SR-SC (Group 6)', row 'k f,on (h - 1)' = 0.083 +/- 0.14 (SE, n = 5)

    lkout_moderator1 <- log(0.45)
    label("Feedback-pool offset rate constant kf,off (1/h)")
    # Table 4, column 'SR-SC (Group 6)', row 'k f,off (h - 1)' = 0.45 +/- 0.04 (SE, n = 5)

    # ---------------- Residual unexplained variability ----------------------
    # Not reported; see vignette Errata and the standing policy for unreported
    # RUV.
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
    frel <- expit(logitfrel)
    frel2 <- (1 - frel) * expit(logitfrel2)
    frel3 <- (1 - frel) * (1 - expit(logitfrel2))

    # ---------------- Sustained-release depot, Equations (4)-(9) ------------
    # Equations (4), (5), (6) and (8); Table 2 row 'ES n' = 2 fixes the
    # erosion chain at two stages.
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2

    # Equation (9).
    d/dt(central) <- ka * depot + ka2 * depot2 + ktr * transit2 - kel * central

    f(depot) <- frel
    f(depot2) <- frel2
    f(transit1) <- frel3

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

    # Equation (13) supplies TT(0). F(0) = 1 is the value at which Equation
    # (14) collapses onto the baseline model of Equation (11); see the
    # Lee_2018_leuprolide_iv_rat file for the full argument.
    TT(0) <- rbase_tt
    moderator1(0) <- 1

    # ---------------- Observations ------------------------------------------
    Cc ~ prop(propSd)
    TT ~ prop(propSd_TT)
  })
}
