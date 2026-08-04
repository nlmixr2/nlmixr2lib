Tsai_2023_ceftriaxone <- function() {
  description <- "Two-compartment population PK model for intravenous ceftriaxone in Indigenous Australian adults with end-stage renal disease on three-times-weekly intermittent high-flux hemodialysis, receiving a novel 2 g three-times-weekly post-dialysis regimen. PK is parameterised on unbound drug: the central state carries unbound ceftriaxone and an explicit second-order albumin-binding exchange (k1 on / k2 off) against a capacity bmax derived from serum albumin carries the bound drug, so total and unbound plasma concentrations are both model outputs. Clearance is replaced (not augmented) by a > 10-fold higher dialytic clearance while a session is running, gated by the time-varying RRT_HEMODIAL_ACTIVE covariate; interdialytic clearance falls with serum bilirubin through an inverse-power relationship. Estimated with the Pmetrics non-parametric adaptive grid (NPAG). Tsai 2023, n = 16 subjects, 122 total-and-unbound plasma samples."
  reference <- "Tsai D, Zam BB, Tongs C, Chiong F, Sajiv C, Pawar B, Ashok A, Cooper BP, Tong SYC, Janson S, Wallis SC, Roberts JA, Parker SL. Validating a novel three-times-weekly post-hemodialysis ceftriaxone regimen in infected Indigenous Australian patients - a population pharmacokinetic study. J Antimicrob Chemother. 2023;78(8):2032-2038. doi:10.1093/jac/dkad190"
  vignette <- "Tsai_2023_ceftriaxone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "unbound ceftriaxone", units = "mg", specimen = "plasma", verified = FALSE),
    complex     = list(analyte = "bound ceftriaxone", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "ceftriaxone", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    TBILI = list(
      description        = "Total serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. The only covariate retained in the final model (Tsai 2023 Results, 'Population pharmacokinetic model and model diagnostics': serum bilirubin was 'the only covariate retained in the final model'). Enters as an inverse-power effect on the interdialytic clearance arm, CL = CLnHD * (14.1 / bili)^0.5 (Tsai 2023 Results equation 'When dialysis is off', reproduced verbatim as 'CL=CL_nHD*(14.1/Bili)**0.5' in the Table S2 Pmetrics model file). Clearance and bilirubin followed an inverse-power relationship with r^2 = 0.74 (Figure 2); with the single highest-bilirubin subject (72 umol/L) removed, r^2 = 0.70. Cohort median 10 umol/L (IQR 6-14); 7 of 16 subjects had bilirubin > 10 umol/L, all acute (2 cholecystitis, 5 severe infection). Must be strictly positive: the covariate enters as a denominator, so TBILI = 0 is undefined. The 14.1 umol/L reference value is hard-coded in the Table S2 model file and the paper does not state its derivation (it is not the cohort median of 10 umol/L); it is carried here exactly as printed.",
      source_name        = "Bili"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. Not a covariate on any structural PK parameter; instead it sets the albumin-binding capacity of the central compartment through the Table S2 secondary variable Bmax1 = Alb * Vc * 16.7 (mg). The 16.7 mg ceftriaxone per g albumin constant encodes 2 binding sites per albumin molecule: 2 * 554.58 (ceftriaxone g/mol) / 66437 (albumin g/mol) * 1000 mg/g = 16.7, consistent with the paper's complex-binding definition (Methods: 'N is the number of ceftriaxone-binding sites per albumin molecule'). Cohort median 36 g/L (IQR 33-39); only one subject had albumin < 30 g/L (at 28 g/L), which the Discussion cites as the reason serum albumin had limited independent influence on the PK model.",
      source_name        = "Alb"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 while an intermittent high-flux hemodialysis session is running, 0 in the interdialytic interval)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no dialysis session running)",
      notes              = "Time-varying within subject. Implemented in the source as the Pmetrics conditional 'IF (HDx.EQ.1) CL = CL_HD' (Table S2 secondary variables), i.e. the dialytic clearance REPLACES the interdialytic clearance arm for the duration of the session rather than being added to it. This is the opposite composition rule from the additive dialysis-arm precedents (Veinstein 2013 gentamicin, Eyler 2014 ertapenem, Jacobs 2016 colistin); it is encoded here as the paper wrote it. Because CL_HD replaces CL entirely, the bilirubin covariate does not act during a dialysis session. Ceftriaxone clearance was > 10-fold higher during dialysis (8.76 vs 0.83 L/h), which the authors attribute to the high-flux membranes used (FX80 / FX100 / FX120, Fresenius); prior studies using low-flux dialyzers reported no dialytic enhancement. Doses in this study were given post-dialysis (within 5 min of the end of the session), so RRT_HEMODIAL_ACTIVE = 0 at the dosing times of the observed data; it is set to 1 only to reproduce the paper's counterfactual 'administered during dialysis' simulation.",
      source_name        = "HDx"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 16L,
    n_studies        = 1L,
    n_samples        = 122L,
    age_median       = "57 years (IQR 51-64); full range not reported",
    weight_median    = "71 kg (IQR 59-83); full range not reported",
    sex_female_pct   = 81.25,
    race_ethnicity   = "100% Indigenous Australian (an explicit inclusion criterion, identified by electronic health record). The Discussion notes that inter-ethnic PK differences are considered unlikely for ceftriaxone, so the authors regard the estimates as transferable to other ethnic origins.",
    disease_state    = "Adults with end-stage renal disease established on three-times-weekly intermittent hemodialysis, treated with ceftriaxone for an active infection. Sources of infection: respiratory 11, urinary 1, bacteremia 1, skin and soft tissue 1, intra-abdominal 1. Baseline laboratory values (median, IQR): albumin 36 g/L (33-39), urea 15.1 mmol/L (12.2-19.2), total bilirubin 10 umol/L (6-14), ALP 216 U/L (171-338), GGT 74 U/L (42-133), ALT 23 U/L (13-29). Seven subjects had bilirubin > 10 umol/L, all of acute cause.",
    renal_function   = "Anuric / end-stage renal disease requiring intermittent hemodialysis three times weekly. Exact renal function was not quantifiable because serum creatinine depended on time since the last dialysis session (stated study limitation). Dialyzers were high-flux: FX80 in 3 subjects (19%), FX100 in 10 (63%), FX120 in 3 (19%) (Fresenius Medical Care).",
    dose_range       = "2 g ceftriaxone (Ceftriaxone-AFT) dissolved in 10 mL water-for-injection, slowly injected into the fistula or central venous cannula within 5 min after the conclusion of each dialysis session, three times weekly.",
    regions          = "Australia (renal dialysis unit of a remote Northern Territory primary referral centre, Alice Springs Hospital)",
    protein_binding  = "Measured directly rather than assumed: median unbound fraction 0.29 (IQR 0.20-0.40), substantially higher than the 0.04-0.17 reported for healthy individuals. Median pre-dialysis unbound trough was 18.2 mg/L (IQR 9.7-25.9) over a 2-day interval and 8.8 mg/L (IQR 7.1-17.7) over a 3-day interval; unbound concentrations fell by a median 70% (IQR 64-74%) across each dialysis session.",
    notes            = "Prospective population PK study. Plasma sampled over two dosing intervals: immediately before and after dialysis, then 5, 15, 60 and 1440 min after the dose, then at 2880 min (or immediately before the next dialysis session for a 48 h interval), and immediately before the next session for a 72 h interval. Total and unbound ceftriaxone assayed by validated UPLC-MS/MS (Table S1); unbound fraction isolated by ultracentrifugation at 37 C with Centrifree devices. One sample was below LLOQ (drawn before any ceftriaxone was given) and was carried into the analysis as 0.0 mg/L. Exclusion criterion: pregnancy."
  )

  ini({
    # Structural parameters: Tsai 2023 Table 2, 'Mean' column of the
    # Pmetrics NPAG non-parametric population distribution. Table 2 also
    # reports a 'Median' column; the mean is used as the typical value
    # here and the median is noted per line. Only the seven primary
    # variables of the Table S2 model file are estimated; Vd and the two
    # half-lives in Table 2 are flagged there as manually calculated
    # derived quantities and are not model parameters.
    lcl <- log(0.83)
    label("Interdialytic (dialysis-off) clearance CLnHD (L/h)")
    # Tsai 2023 Table 2: CLnHD mean 0.83, SD 0.28, CV 33.33%, median 0.77 L/h

    lcl_hemodialysis <- log(8.76)
    label("Intradialytic (dialysis-on) clearance CLHD (L/h)")
    # Tsai 2023 Table 2: CLHD mean 8.76, SD 1.85, CV 21.10%, median 9.13 L/h

    lvc <- log(2.25)
    label("Central volume of distribution Vc (L)")
    # Tsai 2023 Table 2: Vc mean 2.25, SD 1.71, CV 75.81%, median 1.91 L

    lk1 <- log(1.18)
    label("Second-order ceftriaxone-albumin association rate constant Kon (L/mg/h)")
    # Tsai 2023 Table 2: Kon mean 1.18, SD 0.33, CV 27.76%, median 1.18 L/mg/h

    lk2 <- log(206.37)
    label("First-order ceftriaxone-albumin dissociation rate constant Koff (1/h)")
    # Tsai 2023 Table 2: Koff mean 206.37, SD 48.38, CV 23.44%, median 202.36 1/h

    lk12 <- log(15.37)
    label("Central-to-peripheral rate constant Kcp (1/h)")
    # Tsai 2023 Table 2: Kcp mean 15.37, SD 7.24, CV 47.10%, median 17.01 1/h

    lk21 <- log(0.78)
    label("Peripheral-to-central rate constant Kpc (1/h)")
    # Tsai 2023 Table 2: Kpc mean 0.78, SD 0.50, CV 64.05%, median 0.55 1/h

    # Bilirubin effect on the interdialytic clearance arm. The exponent is
    # hard-coded in the Table S2 model file rather than reported as an
    # estimated THETA in Table 2, so it is encoded as fixed().
    e_tbili_cl <- fixed(0.5)
    label("Inverse-power exponent of total bilirubin on interdialytic CL (unitless)")
    # Tsai 2023 Table S2 secondary variables: CL=CL_nHD*(14.1/Bili)**0.5;
    # same form printed in Results ('When dialysis is off'). Equivalent to
    # (TBILI / 14.1)^-0.5. Supported by the Figure 2 inverse-power fit
    # (r^2 = 0.74; r^2 = 0.70 excluding the 72 umol/L subject).

    # Interindividual variability. Pmetrics NPAG estimates a discrete
    # non-parametric distribution rather than a parametric omega matrix;
    # Table 2 summarises that distribution by its mean, SD and CV%. The
    # CV% is carried here into a log-normal random effect using the
    # standard omega^2 = log(CV^2 + 1) identity. This is a parametric
    # APPROXIMATION of a non-parametric distribution (see vignette
    # 'Assumptions and deviations'); it is required to reproduce the
    # paper's own Monte Carlo PTA simulations, which sample the
    # population distribution.
    #   CLnHD : 33.33% CV -> omega^2 = log(0.3333^2 + 1) = 0.105341
    #   CLHD  : 21.10% CV -> omega^2 = log(0.2110^2 + 1) = 0.043558
    #   Vc    : 75.81% CV -> omega^2 = log(0.7581^2 + 1) = 0.454075
    #   Kon   : 27.76% CV -> omega^2 = log(0.2776^2 + 1) = 0.074237
    #   Koff  : 23.44% CV -> omega^2 = log(0.2344^2 + 1) = 0.053487
    #   Kcp   : 47.10% CV -> omega^2 = log(0.4710^2 + 1) = 0.200359
    #   Kpc   : 64.05% CV -> omega^2 = log(0.6405^2 + 1) = 0.343760
    etalcl              ~ 0.105341  # Tsai 2023 Table 2 (CLnHD, CV 33.33%)
    etalcl_hemodialysis ~ 0.043558  # Tsai 2023 Table 2 (CLHD,  CV 21.10%)
    etalvc              ~ 0.454075  # Tsai 2023 Table 2 (Vc,    CV 75.81%)
    etalk1              ~ 0.074237  # Tsai 2023 Table 2 (Kon,   CV 27.76%)
    etalk2              ~ 0.053487  # Tsai 2023 Table 2 (Koff,  CV 23.44%)
    etalk12             ~ 0.200359  # Tsai 2023 Table 2 (Kcp,   CV 47.10%)
    etalk21             ~ 0.343760  # Tsai 2023 Table 2 (Kpc,   CV 64.05%)

    # Residual error. Table S2 '#Error' block gives one assay-error
    # polynomial per output equation, identical for both:
    #   0.3, 0.1, 0, 0   ->  SD = 0.3 + 0.1 * conc  (C2 = C3 = 0)
    # so each output carries a 0.3 mg/L additive plus 10% proportional
    # term. The C1 = 0.1 slope is consistent with the Table S1 assay
    # validation (precision within 8.4% total, 7.5% unbound).
    # NOTE: Pmetrics multiplies this assay polynomial by an estimated
    # noise-inflation factor gamma; the Table S2 file sets the gamma
    # STARTING value 'G=2', and the paper does not report the final
    # estimated gamma anywhere. The assay polynomial is therefore carried
    # here unscaled (equivalent to gamma = 1), which is the minimum-
    # assumption reading of the on-disk file. See vignette 'Assumptions
    # and deviations'.
    addSd <- 0.3
    label("Additive residual error on total Cc (mg/L)")
    # Tsai 2023 Table S2 #Error, output 2 (total): C0 = 0.3
    propSd <- 0.1
    label("Proportional residual error on total Cc (fraction)")
    # Tsai 2023 Table S2 #Error, output 2 (total): C1 = 0.1
    addSd_Cunbound <- 0.3
    label("Additive residual error on unbound Cunbound (mg/L)")
    # Tsai 2023 Table S2 #Error, output 1 (unbound): C0 = 0.3
    propSd_Cunbound <- 0.1
    label("Proportional residual error on unbound Cunbound (fraction)")
    # Tsai 2023 Table S2 #Error, output 1 (unbound): C1 = 0.1
  })

  model({
    # Stoichiometric constant for the albumin-binding capacity, carried
    # exactly as hard-coded in the Table S2 secondary variable
    # Bmax1 = Alb * Vc * 16.7. Units: mg ceftriaxone bound per g albumin.
    # It encodes 2 binding sites per albumin molecule --
    #   2 * 554.58 (ceftriaxone g/mol) / 66437 (albumin g/mol) * 1000 = 16.7
    # -- matching the paper's complex-binding definition (Methods: Bmax is
    # a function of N binding sites, MCTX and MAlb).
    bmax_per_g_alb <- 16.7

    # Reference bilirubin for the inverse-power clearance covariate
    # (umol/L). Hard-coded in the Table S2 model file; the paper does not
    # state its derivation.
    tbili_ref <- 14.1

    # Individual parameters.
    cl              <- exp(lcl + etalcl) * (tbili_ref / TBILI)^e_tbili_cl
    cl_hemodialysis <- exp(lcl_hemodialysis + etalcl_hemodialysis)
    vc              <- exp(lvc + etalvc)
    k1              <- exp(lk1 + etalk1)
    k2              <- exp(lk2 + etalk2)
    k12             <- exp(lk12 + etalk12)
    k21             <- exp(lk21 + etalk21)

    # Dialysis REPLACES the interdialytic clearance arm rather than adding
    # to it (Table S2: 'IF (HDx.EQ.1) CL = CL_HD'). Note this differs from
    # the additive dialysis-arm convention used by Veinstein 2013 /
    # Eyler 2014 / Jacobs 2016; it is encoded as Tsai 2023 wrote it.
    cl_total <- (1 - RRT_HEMODIAL_ACTIVE) * cl + RRT_HEMODIAL_ACTIVE * cl_hemodialysis
    kel      <- cl_total / vc

    # Albumin-binding capacity of the central compartment, as a MASS (mg)
    # rather than a concentration -- so it is directly comparable with the
    # bound-drug amount held in the 'complex' state (Table S2).
    bmax <- ALB * vc * bmax_per_g_alb

    # ODE system, transcribed from the Table S2 '#Differential equations'
    # block. X(1) -> central (unbound drug), X(2) -> complex (albumin-bound
    # drug), X(3) -> peripheral1. Elimination and inter-compartmental
    # distribution act on unbound drug only; the bound state exchanges
    # solely with central.
    #   XP(1) = RATEIV(1) - (Ke + Kcp)*X(1) - (Kon/Vc)*(Bmax1-X(2))*X(1)
    #                     + Koff*X(2) + Kpc*X(3)
    #   XP(2) =             (Kon/Vc)*(Bmax1-X(2))*X(1) - Koff*X(2)
    #   XP(3) =  Kcp*X(1) - Kpc*X(3)
    # The dose enters 'central' (Pmetrics RATEIV(1)) via the event table.
    d/dt(central) <- -(kel + k12) * central -
      (k1 / vc) * (bmax - complex) * central + k2 * complex + k21 * peripheral1
    d/dt(complex) <-
      (k1 / vc) * (bmax - complex) * central - k2 * complex
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Output equations (Table S2 '#Output equations').
    #   Y(1) = X(1)/Vc          -> unbound plasma concentration
    #   Y(2) = (X(2)+X(1))/Vc   -> total plasma concentration
    Cunbound <- central / vc
    Cc       <- (complex + central) / vc

    Cc       ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
  })
}
