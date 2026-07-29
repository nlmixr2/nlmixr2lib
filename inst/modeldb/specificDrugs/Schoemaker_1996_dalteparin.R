Schoemaker_1996_dalteparin <- function() {
  description <- paste(
    "One-compartment population PK/PD model for the low molecular weight",
    "heparin dalteparine (trade name Fragmin) in healthy volunteers, fitted",
    "simultaneously to intravenous and subcutaneous administration data",
    "(Schoemaker & Cohen 1996, Example 3 / Table 4). The kinetic sub-model",
    "is a one-compartment IV bolus / first-order SC absorption disposition",
    "with an estimated constant basal anti-Xa activity (extending the",
    "Schoemaker 1996 Example 2 enoxaparine model with a depot compartment",
    "and bioavailability). The pharmacodynamic sub-model links anti-Xa",
    "activity (Cc) to activated partial thromboplastin time (APTT) through",
    "an exponential concentration-effect relationship parameterised by I10",
    "(the anti-Xa activity increment required to produce a 10% increase in",
    "APTT). Common kinetic parameters are shared between IV and SC routes",
    "within each subject; only F (bioavailability) and ka (absorption rate)",
    "differ between routes. Validation of this model and the companion",
    "enoxaparine PK model share a single vignette."
  )
  reference <- paste(
    "Schoemaker RC, Cohen AF.",
    "Estimating impossible curves using NONMEM.",
    "Br J Clin Pharmacol. 1996 Sep;42(3):283-290.",
    "doi:10.1046/j.1365-2125.1996.04231.x.",
    "Anti-Xa and APTT observation data are from the upstream cross-over study",
    "Kroon C, de Boer A, Kroon JM, Schoemaker RC, Briet E, Cohen AF.",
    "Comparison of the bioavailability of heparin, Fragmin, Fraxiparine and",
    "Orgaran after subcutaneous administration in healthy volunteers.",
    "Br J Clin Pharmacol. 1993;35:548P.",
    "The popPK/PD fits and parameter estimates are original to Schoemaker 1996.",
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
    age_range      = "Not reported in Schoemaker 1996; see the upstream Kroon 1993 abstract.",
    weight_range   = "Not reported in Schoemaker 1996; see the upstream Kroon 1993 abstract.",
    sex_female_pct = NA,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Healthy volunteers in an open randomised cross-over study comparing",
      "the bioavailability of heparin, Fragmin (dalteparine), Fraxiparine",
      "(nadroparine) and Orgaran (danaparoid) after subcutaneous administration.",
      "12 subjects were randomised to eight different occasions each, on which",
      "one of the four drugs was administered either subcutaneously or",
      "intravenously. Only the Fragmin (dalteparine) arm is modelled here",
      "(Schoemaker 1996 Example 3); the same 12 subjects received both IV and",
      "SC dalteparine, enabling simultaneous PK estimation across the two",
      "routes within each subject."
    ),
    dose_range     = paste(
      "Intravenous bolus and subcutaneous administration. The Schoemaker 1996",
      "paper does not list the absolute dose amount in either anti-Xa IU or mg;",
      "see the upstream Kroon 1993 abstract for the dose specification. Figure 4",
      "(subject 3) shows peak anti-Xa activity in the 0.3-0.4 IU/mL range and",
      "peak APTT response of 55-60 seconds (versus baseline ~30 s)."
    ),
    regions        = "Single-centre, Centre for Human Drug Research, Leiden, The Netherlands.",
    notes          = paste(
      "Baseline demographics are reported in the upstream Kroon 1993 abstract",
      "(reference [11] of Schoemaker 1996; Br J Clin Pharmacol 35:548P) which",
      "was NOT on disk at extraction time. The Schoemaker 1996 paper itself",
      "reports only that 12 subjects contributed to the cross-over study.",
      "The complex multi-route PK/PD fit was advanced in steps to ensure",
      "successful NONMEM convergence: first the IV data were analysed alone;",
      "then all data (IV + SC) were analysed together with the IV-derived PK",
      "parameter estimates fixed to recover the bioavailability and absorption",
      "estimates; a final analysis was run with no parameters restricted.",
      "Initial absorption rate IIV (CV on t1/2A) was estimated at 0% and",
      "encoded here as fixed(0). Estimation proceeded via the first-order",
      "method with a final conditional-method run (Schoemaker 1996 Methods);",
      "NONMEM Version IV / NMTRAN II / PREDPP III, FORTRAN PowerStation 1.0",
      "on MS-DOS."
    )
  )

  ini({
    # Final NONMEM estimates from Schoemaker 1996 Table 4 (one-compartment IV
    # bolus + first-order SC absorption + constant basal anti-Xa activity;
    # exponential anti-Xa -> APTT concentration-effect model). Table 4
    # reports times (min), CL (mL/min), F (%), Base (IU/mL), APTT0 (s),
    # and I10 (IU/mL). Values are converted to h-based units; the typical
    # central volume is derived from the typical CL and elimination t1/2E
    # via Vc = CL * t1/2E / ln(2).
    lka     <- log(0.5429)
    label("Absorption rate constant (ka, 1/h)")
    # Schoemaker 1996 Table 4: t1/2A = 76.6 min; ka = ln(2)/(76.6/60) h = 0.5429 1/h
    lcl     <- log(2.766)
    label("Clearance (CL, L/h)")
    # Schoemaker 1996 Table 4: CL = 46.1 mL/min; converted: 46.1 * 60 / 1000 = 2.766 L/h
    lvc     <- log(6.166)
    label("Central volume of distribution (Vc, L)")
    # Schoemaker 1996 Table 4: derived from CL = 46.1 mL/min and t1/2E = 92.7 min via
    # Vc = CL * t1/2E / ln(2) = 2.766 L/h * (92.7/60) h / 0.6931 = 6.166 L
    lfdepot <- log(0.705)
    label("Bioavailability for subcutaneous administration (fraction)")
    # Schoemaker 1996 Table 4: F = 70.5% (s.e. 7.81%; CV 34.8%)
    lrbase  <- log(0.0192)
    label("Basal anti-Xa activity baseline (IU/mL)")
    # Schoemaker 1996 Table 4: Base = 0.0192 IU/mL (s.e. 0.00271; CV 16.6%)

    # PD parameters. The Schoemaker 1996 exponential anti-Xa -> APTT
    # concentration-effect model (their equation 10):
    #   APTT_j = APTT0 * exp((ln(1.1)/I10) * Xa_j) + e2_j
    # The I10 parameter is the anti-Xa activity (IU/mL) that produces a 10%
    # increase in APTT relative to the basal APTT level. The authors chose
    # this somewhat unconventional reparameterisation because (i) the
    # individual concentration-effect plots showed a slightly curvilinear
    # relationship with no indication of a maximum effect, and (ii) the
    # alternative exponential-slope form e^(slope * Xa) hides the natural
    # interpretation that I10 provides.
    lrbase_APTT <- log(30.2)
    label("Basal APTT level (APTT0, s)")
    # Schoemaker 1996 Table 4: APTT0 = 30.2 s (s.e. 0.910; CV 10.6%)
    lI10        <- log(0.0665)
    label("Anti-Xa activity increment producing a 10% APTT increase (I10, IU/mL)")
    # Schoemaker 1996 Table 4: I10 = 0.0665 IU/mL (s.e. 0.00439; CV 19.6%)

    # IIV. Schoemaker 1996 used the constant-coefficient-of-variation log-normal
    # IIV model (their equation 6): theta_i = theta * exp(eta_i). The variance
    # of each eta on the log-scale equals log(1 + CV^2). The paper reports IIV
    # on CL (14.1%) and elimination t1/2E (8.3%) as independent etas in the
    # (CL, t1/2E) parameterisation. Re-expressed in the (CL, Vc) parameterisation
    # used here (Vc = CL * t1/2E / ln(2), so eta_lvc = eta_lcl + eta_lt12E with
    # eta_lt12E independent of eta_lcl), the implied joint distribution is:
    #   Var(eta_lcl) = log(1 + 0.141^2) = 0.01970
    #   Var(eta_lvc) = Var(eta_lcl) + Var(eta_lt12E)
    #                = 0.01970 + log(1 + 0.083^2) = 0.02657
    #   Cov(eta_lcl, eta_lvc) = Var(eta_lcl) = 0.01970
    # (Correlation = 0.861). This block-diagonal encoding preserves the paper's
    # reported IIV on CL (14.1%) and t1/2E (8.3%) exactly: when sampled, the
    # implied IIV on t1/2E = sqrt(Var(eta_lvc - eta_lcl)) = sqrt(0.00687) ~ 8.3%
    # and IIV on CL = sqrt(0.01970) ~ 14.1%, matching Table 4.
    etalcl + etalvc ~ c(0.01970, 0.01970, 0.02657)
    # Schoemaker 1996 Table 4: CV(CL) = 14.1%, CV(t1/2E) = 8.3%

    # Absorption-rate IIV: Schoemaker 1996 Table 4 reports CV(t1/2A) = 0.0%
    # (no inter-individual variability estimated on the absorption half-life;
    # the authors note that absorption and elimination half-lives cannot be
    # separated very well without independent IV information, and in the
    # final joint IV+SC fit the absorption-half-life IIV collapsed to 0).
    # No etalka entry is declared so the simulated population has a single
    # shared absorption-rate trajectory matching the paper's fixed-ka fit.

    etalfdepot     ~ 0.1143
    # Schoemaker 1996 Table 4: CV(F) = 34.8% -> log(1 + 0.348^2) = 0.1143
    etalrbase      ~ 0.02719
    # Schoemaker 1996 Table 4: CV(Base) = 16.6% -> log(1 + 0.166^2) = 0.02719
    etalrbase_APTT ~ 0.01117
    # Schoemaker 1996 Table 4: CV(APTT0) = 10.6% -> log(1 + 0.106^2) = 0.01117
    etalI10        ~ 0.03770
    # Schoemaker 1996 Table 4: CV(I10) = 19.6% -> log(1 + 0.196^2) = 0.03770

    # Residual error. Schoemaker 1996 Table 4 reports two independent additive
    # (constant-variance, not constant-CV) residuals: one on anti-Xa activity
    # (Cc, IU/mL) in the kinetic sub-model, and one on APTT (s) in the
    # exponential PD sub-model. Both are applied in the original linear scale
    # because anti-Xa assay precision does not improve at low concentrations
    # (Schoemaker 1996 Example 2 Solution 2 discussion).
    addSd      <- 0.0207
    label("Additive residual SD on anti-Xa activity (IU/mL)")
    # Schoemaker 1996 Table 4: s.d. residual error anti-Xa vs time = 0.0207 IU/mL
    addSd_APTT <- 1.70
    label("Additive residual SD on APTT (s)")
    # Schoemaker 1996 Table 4: s.d. residual error APTT vs anti-Xa = 1.70 s
  })
  model({
    # Individual structural parameters. No etalka -- Schoemaker 1996 Table 4
    # reports CV(t1/2A) = 0%, so ka is a typical-value-only parameter.
    ka         <- exp(lka)
    cl         <- exp(lcl + etalcl)
    vc         <- exp(lvc + etalvc)
    fdepot     <- exp(lfdepot + etalfdepot)
    rbase      <- exp(lrbase + etalrbase)
    rbase_APTT <- exp(lrbase_APTT + etalrbase_APTT)
    I10        <- exp(lI10 + etalI10)

    # Micro-constant.
    kel <- cl / vc

    # ODE. The depot compartment receives subcutaneous doses (cmt = "depot");
    # the central compartment receives intravenous bolus doses
    # (cmt = "central"). For SC doses, bioavailability fdepot multiplies the
    # absorbed amount via f(depot) below; for IV doses, the bolus goes
    # directly into central with F = 1 (rxode2 default).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability of the depot compartment (subcutaneous route).
    f(depot) <- fdepot

    # Observation 1: anti-Xa activity (IU/mL). central is in anti-Xa IU and
    # vc is in L, so central/vc is in IU/L; divide by 1000 to convert to
    # IU/mL. The additive basal activity reflects lingering low / endogenous
    # anti-Xa-like activity that survives below the assay's detection
    # capability (Schoemaker 1996 equation 9).
    Cc <- central / vc / 1000 + rbase
    Cc ~ add(addSd)

    # Observation 2: APTT (s). Exponential anti-Xa -> APTT concentration-
    # effect model (Schoemaker 1996 equation 10):
    #   APTT_j = APTT0 * exp((ln(1.1) / I10) * Xa_j) + e2_j
    # where Xa_j is the model-predicted anti-Xa activity at time j (above)
    # and e2_j is the additive APTT residual (s.d. addSd_APTT). The
    # ln(1.1)/I10 form is the natural-log slope of the exponential effect
    # curve, parameterised so that I10 is the anti-Xa activity (IU/mL)
    # required to produce a 10% APTT increase above APTT0.
    APTT <- rbase_APTT * exp((log(1.1) / I10) * Cc)
    APTT ~ add(addSd_APTT)
  })
}
