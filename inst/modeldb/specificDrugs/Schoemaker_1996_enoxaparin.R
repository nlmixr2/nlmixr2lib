Schoemaker_1996_enoxaparin <- function() {
  description <- paste(
    "One-compartment population PK model with intravenous bolus input and",
    "an estimated constant basal anti-Xa activity for the low molecular",
    "weight heparin enoxaparine (trade name Clexane) in healthy volunteers",
    "(Schoemaker & Cohen 1996, Example 2 / Table 3, Solution 2). Enoxaparin",
    "amount in the central compartment plus an additive endogenous baseline",
    "reproduces the lingering low post-dose anti-Xa activity that would",
    "otherwise force a second compartment if pre-value subtraction were",
    "applied; the authors recommend the basal-activity formulation over the",
    "competing two-compartment model (Solution 1, Table 2) because it",
    "matches the dose / AUC clearance estimate from the upstream",
    "Stiekema 1993 paper. Anti-Xa activity is the surrogate concentration",
    "measure; doses are in anti-Xa IU and concentration is in IU/mL.",
    "Validation of this model and the companion dalteparin PK/PD model",
    "share a single vignette."
  )
  reference <- paste(
    "Schoemaker RC, Cohen AF.",
    "Estimating impossible curves using NONMEM.",
    "Br J Clin Pharmacol. 1996 Sep;42(3):283-290.",
    "doi:10.1046/j.1365-2125.1996.04231.x.",
    "Anti-Xa observation data are from the upstream cross-over study",
    "Stiekema JCJ, van Griensven JMT, van Dinther TG, Cohen AF.",
    "A cross-over comparison of the anti-clotting effects of three low",
    "molecular weight heparins and glycosaminoglycuronan.",
    "Br J Clin Pharmacol. 1993;35:51-56.",
    "The popPK fits and parameter estimates are original to Schoemaker 1996.",
    sep = " "
  )
  vignette <- "Schoemaker_1996_low_molecular_weight_heparin_modeling"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "Not reported in Schoemaker 1996; see the upstream Stiekema 1993 paper.",
    weight_range   = "Not reported in Schoemaker 1996; see the upstream Stiekema 1993 paper.",
    sex_female_pct = NA,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Healthy volunteers in an open randomised cross-over study comparing",
      "the anti-clotting effects of three low molecular weight heparins",
      "(Fragmin / Clexane / Fraxiparine) and the glycosaminoglycuronan",
      "danaparoid (Orgaran). 12 subjects received intravenous administrations",
      "of all four drugs on separate occasions. Only the Clexane (enoxaparine)",
      "arm is modelled here (Schoemaker 1996 Example 2)."
    ),
    dose_range     = paste(
      "Intravenous bolus. The Schoemaker 1996 paper does not list the absolute",
      "dose amount in either anti-Xa IU or mg; see the upstream Stiekema 1993",
      "report for the dose specification. Figure 2 (subjects 4 and 5) shows",
      "peak anti-Xa activity in the 0.4-1.0 IU/mL range, consistent with a",
      "single IV bolus of approximately 5000 anti-Xa IU (~50 mg of enoxaparine)."
    ),
    regions        = "Single-centre, Centre for Human Drug Research, Leiden, The Netherlands.",
    notes          = paste(
      "Baseline demographics are reported in the upstream Stiekema 1993 paper",
      "(reference [10] of Schoemaker 1996) which was NOT on disk at extraction",
      "time. The Schoemaker 1996 paper itself reports only that 12 subjects",
      "contributed to the cross-over study. With ordinary nonlinear regression,",
      "five subjects were adequately fit by a mono-exponential function, four",
      "required a bi-exponential function, and three provided no adequate fit",
      "at all; NONMEM with the one-compartment + basal activity model fitted",
      "all 12 subjects simultaneously by borrowing information across the",
      "cohort. Estimation proceeded via the first-order method with a final",
      "conditional-method run (Schoemaker 1996 Methods); NONMEM Version IV /",
      "NMTRAN II / PREDPP III, FORTRAN PowerStation 1.0 on MS-DOS."
    )
  )

  ini({
    # Final NONMEM estimates from Schoemaker 1996 Table 3 (Solution 2: one-
    # compartment IV bolus with an estimated constant basal anti-Xa activity).
    # Table 3 reports t1/2 (min), CL (mL/min) and Base (IU/mL). Values are
    # converted to h-based units to match the package convention; the typical
    # central volume is derived from the typical CL and t1/2 via the standard
    # one-compartment relation Vc = CL * t1/2 / ln(2).
    lcl    <- log(1.974)
    label("Clearance (CL, L/h)")
    # Schoemaker 1996 Table 3: CL = 32.9 mL/min; converted: 32.9 mL/min * 60 / 1000 = 1.974 L/h
    lvc    <- log(6.170)
    label("Central volume of distribution (Vc, L)")
    # Schoemaker 1996 Table 3: derived from CL = 32.9 mL/min and t1/2 = 130 min via
    # Vc = CL * t1/2 / ln(2) = 1.974 L/h * (130/60) h / 0.6931 = 6.170 L
    lrbase <- log(0.0254)
    label("Basal anti-Xa activity baseline (IU/mL)")
    # Schoemaker 1996 Table 3: Base = 0.0254 IU/mL (s.e. 0.00333; CV 37.7%)

    # IIV. Schoemaker 1996 used the constant-coefficient-of-variation log-normal
    # IIV model (their equation 6): theta_i = theta * exp(eta_i). The variance
    # of each eta on the log-scale equals log(1 + CV^2). The paper reports IIV
    # on CL (15.1%) and t1/2 (8.6%) as independent etas in the (CL, t1/2)
    # parameterisation. Re-expressed in the (CL, Vc) parameterisation used
    # here (Vc = CL * t1/2 / ln(2), so eta_lvc = eta_lcl + eta_lt12 with
    # eta_lt12 independent of eta_lcl), the implied joint distribution is:
    #   Var(eta_lcl) = log(1 + 0.151^2) = 0.02256
    #   Var(eta_lvc) = Var(eta_lcl) + Var(eta_lt12)
    #                = 0.02256 + log(1 + 0.086^2) = 0.02992
    #   Cov(eta_lcl, eta_lvc) = Var(eta_lcl) = 0.02256
    # (Correlation = 0.868). This block-diagonal encoding preserves the paper's
    # reported IIV on CL (15.1%) and t1/2 (8.6%) exactly: when sampled, the
    # implied IIV on t1/2 = sqrt(Var(eta_lvc - eta_lcl)) = sqrt(0.00736) ~ 8.6%
    # and IIV on CL = sqrt(0.02256) ~ 15.1%, matching Table 3.
    etalcl + etalvc ~ c(0.02256, 0.02256, 0.02992)
    # Schoemaker 1996 Table 3: CV(CL) = 15.1%, CV(t1/2) = 8.6%
    etalrbase ~ 0.1329
    # Schoemaker 1996 Table 3: CV(Base) = 37.7% -> log(1 + 0.377^2) = 0.1329

    # Residual error. Schoemaker 1996 Table 3 reports an additive (constant-
    # variance, not constant-CV) residual error on anti-Xa activity, applied
    # in the original linear scale. The paper notes that an additive error is
    # appropriate here because the anti-Xa assay is not very accurate in the
    # lower concentration range (variance does not decrease at low Cc as it
    # would for typical PK data).
    addSd <- 0.0274
    label("Additive residual SD on anti-Xa activity (IU/mL)")
    # Schoemaker 1996 Table 3: s.d. residual error = 0.0274 IU/mL
  })
  model({
    # Individual structural parameters.
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    rbase <- exp(lrbase + etalrbase)

    # Micro-constant. One-compartment IV bolus elimination.
    kel <- cl / vc

    # ODE. Dose is administered as an IV bolus directly into the central
    # compartment (cmt = "central" in the event table). The state holds the
    # drug amount in anti-Xa IU.
    d/dt(central) <- -kel * central

    # Observation. central / vc gives anti-Xa activity in IU/L; divide by
    # 1000 to convert to IU/mL (matching the Schoemaker 1996 Figure 2 axis
    # and Table 3 baseline units). The additive basal activity reflects
    # endogenous anti-Xa-like activity that lingers below the assay's
    # post-dose detection capability (Schoemaker 1996 Solution 2 prose,
    # equation 9: C_j = base + (D/V) * exp(-k * t_j)).
    Cc <- central / vc / 1000 + rbase
    Cc ~ add(addSd)
  })
}
