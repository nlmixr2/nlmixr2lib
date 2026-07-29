Zvada_2010_rifapentine <- function() {
  paper_specific_etas <- c("etalcl", "etalcl_desrpt")
  # etalcl is shared by both CL1/F and CL2/F (the pre- and post-MTIME apparent
  # clearances) per Zvada 2010 Table 2, which reports identical IIV 19.2% and
  # IOV 12.2% for the two rows -- consistent with a single random effect carrying
  # through the autoinduction step. Similarly etalcl_desrpt is shared between
  # CLM1/F and CLM2/F (Zvada 2010 Table 4 identical IIV 23.9%). Neither eta
  # therefore pairs 1-to-1 with a single lX parameter; both are declared as
  # paper_specific_etas so checkModelConventions() does not flag the
  # one-to-one-matching convention.
  description <- "Parent-metabolite population pharmacokinetic model for single-dose 900 mg oral rifapentine (RFP) and its primary active metabolite 25-O-desacetyl rifapentine (25-DRFP) in 34 healthy adult male volunteers, with characterization of food effect on bioavailability for four meal types (high-fat English breakfast A, low-fat bulky maize porridge B, high-fat bulky maize porridge with lard C, and low-fat high-fluid chicken noodle soup D) relative to fasted reference (meal E). Parent RFP is described by a one-compartment model with Savic transit absorption (NN = 10.9, MTT = 1.45 h) and a step-function autoinduction of apparent clearance at MTIME = 43 h (CL1/F = 2.14 to CL2/F = 3.22 L/h). All RFP is assumed to convert to 25-DRFP; metabolite disposition is two-compartment with its own step-function clearance switch at MTIME_M = 46.8 h (CLM1/F = 1.81 to CLM2/F = 4.63 L/h). Meal effects multiply the typical bioavailability via TVF = 1 * (1 + sum of per-meal fractional changes). Inter-individual variability is a 3x3 block on CL/F, MTT, and F (correlations rho_F_MTT = 0.65 and rho_CL_MTT = -0.56; cov(CL, F) assumed 0 since not reported); a single shared eta on CL/F applies to both CL1 and CL2. Inter-occasion variability (Table 2 IOV columns) is omitted because nlmixr2lib does not standardize an OCC encoding; the cross-over IOV magnitudes are noted in the vignette Errata. Residual variability is combined additive + proportional on plasma RFP (additive 0.206 mg/L, proportional 10.6%) and on plasma 25-DRFP (additive 0.211 mg/L, proportional 19.1%)."
  reference <- paste(
    "Zvada SP, Van Der Walt JS, Smith PJ, Fourie PB, Roscigno G,",
    "Mitchison D, Simonsson USH, McIlleron HM. (2010).",
    "Effects of four different meal types on the population",
    "pharmacokinetics of single-dose rifapentine in healthy male volunteers.",
    "Antimicrob Agents Chemother 54(8):3390-3394.",
    "doi:10.1128/AAC.00345-10.",
    sep = " "
  )
  vignette <- "Zvada_2010_rifapentine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    MEAL_A = list(
      description        = "Binary indicator for meal A administered 30 min before the oral rifapentine dose (1 = meal A, 0 = otherwise). Meal A is a high-fat English breakfast (Zvada 2010 Table 1: 2 rashers of bacon (20 g), 1 fried egg (50 g), 1 slice white toast (30 g) with butter (7 g) and marmalade (10 g), 2 cups decaffeinated coffee (400 ml) with full-cream milk (100 ml) and sugar (10 g); 18.9 g protein, 27 g fat, 38 g carbohydrate, 1,966 kJ, 627 g total weight).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted = meal E, 200 ml of water)",
      notes              = "Mutually exclusive with MEAL_B, MEAL_C, and MEAL_D - at most one of the four meal indicators can be 1 on any single dosing record; the fasted reference is encoded as all four = 0. Per Zvada 2010 Table 3 meal A produced the largest food effect, increasing RFP oral bioavailability by 85.7% (RSE 20.1%) relative to fasting.",
      source_name        = "meal A (Zvada 2010 Table 1 / Table 3)"
    ),
    MEAL_B = list(
      description        = "Binary indicator for meal B administered 30 min before the oral rifapentine dose (1 = meal B, 0 = otherwise). Meal B is a low-fat bulky maize-meal porridge breakfast (Zvada 2010 Table 1: 1.5 cups soft maize meal porridge (375 g cooked) with 3 teaspoons sugar (15 g), 1 cup decaffeinated coffee (200 ml) with full-cream milk (50 ml) and 1 teaspoon sugar (5 g); 6 g protein, 3 g fat, 66 g carbohydrate, 1,285 kJ, 645 g total weight).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted = meal E)",
      notes              = "Mutually exclusive with MEAL_A, MEAL_C, MEAL_D. Zvada 2010 Table 3 reports a +32.7% (RSE 40.4%) increase in RFP oral bioavailability under meal B relative to fasting.",
      source_name        = "meal B (Zvada 2010 Table 1 / Table 3)"
    ),
    MEAL_C = list(
      description        = "Binary indicator for meal C administered 30 min before the oral rifapentine dose (1 = meal C, 0 = otherwise). Meal C is a high-fat bulky maize-meal porridge breakfast (Zvada 2010 Table 1: 1.5 cups soft maize meal porridge (375 g cooked) with 3 teaspoons sugar (15 g) and 5 teaspoons of lard (25 g), 1 cup decaffeinated coffee (200 ml) with full-cream milk (50 ml) and 1 teaspoon sugar (5 g); 6 g protein, 28 g fat, 66 g carbohydrate, 2,229 kJ, 670 g total weight).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted = meal E)",
      notes              = "Mutually exclusive with MEAL_A, MEAL_B, MEAL_D. Zvada 2010 Table 3 reports a +45.7% (RSE 29.3%) increase in RFP oral bioavailability under meal C relative to fasting; the paper Discussion notes this is smaller than meal A's +85.7% despite similar total fat content, supporting the authors' hypothesis that eggs in meal A are independently relevant.",
      source_name        = "meal C (Zvada 2010 Table 1 / Table 3)"
    ),
    MEAL_D = list(
      description        = "Binary indicator for meal D administered 30 min before the oral rifapentine dose (1 = meal D, 0 = otherwise). Meal D is a low-fat high-fluid chicken noodle soup breakfast (Zvada 2010 Table 1: 2 cups reconstituted powdered chicken noodle soup (400 ml), 1 cup decaffeinated coffee (200 ml) with skim milk (50 ml) and 1 teaspoon sugar (5 g); 9 g protein, 4 g fat, 28 g carbohydrate, 774 kJ, 660 g total weight).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted = meal E)",
      notes              = "Mutually exclusive with MEAL_A, MEAL_B, MEAL_C. Zvada 2010 Table 3 reports a +48.9% (RSE 30.7%) increase in RFP oral bioavailability under meal D relative to fasting; the paper Discussion attributes part of this effect to monosodium glutamate (MSG) in the soup accelerating gastric emptying.",
      source_name        = "meal D (Zvada 2010 Table 1 / Table 3)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "23.9 (4.82) years (mean (SD))",
    age_median     = NULL,
    weight_range   = "74.4 (12.3) kg (mean (SD))",
    weight_median  = NULL,
    height_mean    = "177.2 (7.33) cm",
    sex_female_pct = 0,
    race_ethnicity = "not specified; the study was conducted in Cape Town, South Africa.",
    disease_state  = "Healthy adult male volunteers (no TB; HIV- and HBV-negative; no clinically relevant cardiovascular, hepatic, neurologic, endocrine, or other major systemic disease).",
    dose_range     = "Single 900 mg oral rifapentine (6 x 150 mg Priftin tablets; Hoechst Marion Roussel, Italy) with 200 ml of water, given 30 min after each of five test meals (A, B, C, D) or fasted (E) in a five-way crossover with a 14-day washout between occasions.",
    regions        = "Cape Town, South Africa (Groote Schuur Hospital).",
    sampling_design = "Open-label, randomized, sequential, five-way crossover. At each visit a 20-ml predose sample was collected and 10-ml samples were drawn at 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 24, 36, 48, and 72 h after the dose. Each dose was separated by a 14-day washout period. Seven participants were randomized to each of five different meal sequences. 2,272 RFP and 25-DRFP plasma concentration samples available from 34 participants (one withdrew after a single visit). Less than 1% of concentration-time data were below the limit of quantification.",
    assay          = "Validated tandem HPLC method (Division of Clinical Pharmacology, Cape Town) over the concentration range 0.5 to 30 ug/ml; calibration curve linearity r^2 = 0.9975 for RFP and r^2 = 0.9946 for 25-DRFP.",
    notes          = paste(
      "Inclusion criteria: weight >= 50 kg, normal physical examination and",
      "baseline laboratory evaluation, nonsmoker. Exclusion criteria: history of TB,",
      "active allergies, excessive coffee or alcohol consumption, recent blood",
      "donation > 500 ml, clinically relevant systemic disease, women (to limit",
      "fetal exposure risk), HIV+, HBV+. The study protocol (M000473/1LO1) was",
      "approved by the Research Ethics Committee of the University of Cape Town",
      "and the Medicines Control Council of South Africa.",
      "Three volunteers were withdrawn for elevated ALT (one after first dose at",
      "2x ULN; two after second dose at 2.2x and 6.9x ULN); levels returned to",
      "normal after withdrawal and there was no relationship between toxicity and",
      "drug concentration or meal type."
    )
  )

  ini({
    # =========================================================================
    # Rifapentine (parent) structural parameters - one-compartment model with
    # Savic transit absorption (NN, MTT, ka) and step-function autoinduction of
    # apparent clearance at MTIME = 43 h. All values from Zvada 2010 Table 2.
    # Time unit is hours throughout (matches the source).
    # =========================================================================
    lcl1   <- log(2.14)
    label("Rifapentine apparent oral clearance CL1/F before MTIME (L/h)")
    # Zvada 2010 Table 2 'CL1/F (liters/h)' = 2.14 (RSE 13.6%). Baseline clearance
    # before the autoinduction step at MTIME = 43 h. The same eta_CL applies to
    # CL1 and CL2 (Zvada Table 2 reports identical IIV 19.2% and IOV 12.2% for the
    # two rows, consistent with a single shared random effect).

    lcl2   <- log(3.22)
    label("Rifapentine apparent oral clearance CL2/F after MTIME (L/h)")
    # Zvada 2010 Table 2 'CL2/F (liters/h)' = 3.22 (RSE 11.9%). Autoinduced
    # clearance after MTIME; CL2/CL1 ratio = 1.50.

    lt_switch <- fixed(log(43))
    label("Rifapentine clearance autoinduction switch time (h)")
    # Zvada 2010 Table 2 'MTIME (h)' = 43 (RSE 2.6%). Reported without an IIV /
    # IOV column in Table 2 and described as the population estimate of the
    # event time at which CL/F changes; held fixed in this implementation since
    # no random effect on MTIME was reported. Symbol renamed from MTIME to
    # t_switch in this nlmixr2 implementation because mtime() is a reserved
    # rxode2 modeling function (used to declare model-event times); the
    # autoinduction step is encoded inside model() via tad(depot).

    lvc    <- log(60.6)
    label("Rifapentine apparent central volume V/F (L)")
    # Zvada 2010 Table 2 'V/F (liters)' = 60.6 (RSE 9.2%). No IIV reported.

    lka    <- log(1.66)
    label("Rifapentine first-order absorption rate constant ka (1/h)")
    # Zvada 2010 Table 2 'ka (h^-1)' = 1.66 (RSE 13.1%). No IIV reported.

    lmtt   <- log(1.45)
    label("Rifapentine mean transit time MTT (h)")
    # Zvada 2010 Table 2 'MTT (h)' = 1.45 (RSE 10.8%). Savic transit-absorption
    # chain mean transit time.

    lnn    <- log(10.9)
    label("Rifapentine Savic transit-chain shape parameter NN (unitless)")
    # Zvada 2010 Table 2 'NN' = 10.9 (RSE 9.6%). Non-integer; rxode2's
    # transit(n, mtt, bio) closed form accepts non-integer n via the gamma-density
    # evaluation.

    lfdepot <- fixed(log(1))
    label("Rifapentine reference oral bioavailability F (FIXED to 1 as fasting anchor)")
    # Zvada 2010 Table 2 'F = 1 (fixed)' and Results paragraph 1:
    # 'TVF = theta_F * (1 + RXF); theta_F is the value of F under fasting
    # conditions (fixed to 1)'. The meal effects are encoded below as fractional
    # increments multiplied onto this 1.0 anchor.

    # =========================================================================
    # Meal effects on rifapentine bioavailability (Zvada 2010 Table 3). Each
    # meal increases F relative to the fasted (meal E) reference; the encoding
    # is TVF = 1 * (1 + sum of per-meal fractional changes), so the four
    # coefficients below are the per-meal RXF terms.
    # =========================================================================
    e_meal_a_fdepot <- 0.857
    label("Fractional increase in rifapentine F for meal A vs fasted (unitless)")
    # Zvada 2010 Table 3 'Meal A: +85.7%' (RSE 20.1%). High-fat English breakfast.

    e_meal_b_fdepot <- 0.327
    label("Fractional increase in rifapentine F for meal B vs fasted (unitless)")
    # Zvada 2010 Table 3 'Meal B: +32.7%' (RSE 40.4%). Low-fat bulky maize porridge.

    e_meal_c_fdepot <- 0.457
    label("Fractional increase in rifapentine F for meal C vs fasted (unitless)")
    # Zvada 2010 Table 3 'Meal C: +45.7%' (RSE 29.3%). High-fat bulky maize porridge with lard.

    e_meal_d_fdepot <- 0.489
    label("Fractional increase in rifapentine F for meal D vs fasted (unitless)")
    # Zvada 2010 Table 3 'Meal D: +48.9%' (RSE 30.7%). Low-fat high-fluid chicken noodle soup.

    # =========================================================================
    # 25-O-desacetyl rifapentine (active metabolite, suffix '_desrpt')
    # Two-compartment disposition with step-function autoinduction of apparent
    # clearance at MTIME_M = 46.8 h. The parent-to-metabolite formation flux is
    # the full RFP elimination rate (CL_RFP / V_RFP * central_RFP) since the
    # paper assumes all RFP is converted to 25-DRFP (Figure 1 caption);
    # metabolite parameters are apparent (with respect to RFP F and the
    # implicit unit conversion between RFP mass and metabolite mass). All
    # values from Zvada 2010 Table 4.
    # =========================================================================
    lcl1_desrpt <- log(1.81)
    label("25-DRFP apparent oral clearance CLM1/F before MTIME_M (L/h)")
    # Zvada 2010 Table 4 'CLM1/F (liters/h)' = 1.81 (RSE 8.7%). Same single
    # eta_CLM shared between CLM1 and CLM2 (matching the identical IIV 23.9%
    # and IOV 30.2% reported for both rows).

    lcl2_desrpt <- log(4.63)
    label("25-DRFP apparent oral clearance CLM2/F after MTIME_M (L/h)")
    # Zvada 2010 Table 4 'CLM2/F (liters/h)' = 4.63 (RSE 8.1%); CLM2/CLM1 = 2.56.

    lt_switch_desrpt <- fixed(log(46.8))
    label("25-DRFP clearance autoinduction switch time (h)")
    # Zvada 2010 Table 4 'MTIME (h)' = 46.8 (RSE 0.3%). Held fixed since no
    # random effect on MTIME was reported. Renamed from MTIME to t_switch_desrpt
    # in this implementation; see lt_switch label note for rationale.

    lvc_desrpt <- log(6.36)
    label("25-DRFP apparent central volume VMC/F (L)")
    # Zvada 2010 Table 4 'V_MC/F (liters)' = 6.36 (RSE 5.8%). No IIV reported.

    lvp_desrpt <- log(22.1)
    label("25-DRFP apparent peripheral volume VMP/F (L)")
    # Zvada 2010 Table 4 'V_MP/F (liters)' = 22.1 (RSE 4.8%). No IIV reported.

    lq_desrpt <- log(4.4)
    label("25-DRFP apparent inter-compartmental clearance QM/F (L/h)")
    # Zvada 2010 Table 4 'QM (liters/h)' = 4.4 (RSE 6.9%). No IIV reported.

    # =========================================================================
    # Inter-individual variability (Zvada 2010 Table 2 / Table 4 'IIV' columns,
    # reported as CV percent). Internal-scale variance is omega^2 = log(1 + CV^2).
    # Parent block: 3x3 correlated etas on CL/F, MTT, F per Discussion paragraph 4:
    #   'correlation between F and MTT (correlation coefficient = 0.65) and
    #    between CL/F and MTT (correlation coefficient = -0.56) ... in the RFP model'.
    # The cov(CL, F) entry is set to 0 because the paper does not report it.
    # The block ordering is etalcl + etalmtt + etalfdepot, with the lower-
    # triangle covariance entries c(var_cl, cov_cl_mtt, var_mtt, cov_cl_f,
    # cov_mtt_f, var_f).
    # =========================================================================
    # Lower-triangle order for the 3x3 block:
    #   var(cl), cov(cl,mtt), var(mtt), cov(cl,f), cov(mtt,f), var(f)
    # Diagonal variances are log(1 + CV^2):
    #   var(cl)     = log(1 + 0.192^2) = 0.0362008 (Zvada 2010 Table 2 IIV CL/F = 19.2% CV, RSE 24.3%)
    #   var(mtt)    = log(1 + 0.095^2) = 0.0089845 (Zvada 2010 Table 2 IIV MTT  =  9.5% CV, RSE 87.0%)
    #   var(f)      = log(1 + 0.214^2) = 0.0447783 (Zvada 2010 Table 2 IIV F    = 21.4% CV, RSE 35.3%)
    # Off-diagonal covariances (Zvada 2010 Discussion paragraph 4):
    #   cov(cl, mtt) = -0.56 * sqrt(0.0362008 * 0.0089845) = -0.0100994 (rho = -0.56)
    #   cov(f, mtt)  =  0.65 * sqrt(0.0447783 * 0.0089845) =  0.0130375 (rho =  0.65)
    #   cov(cl, f)   = 0 (rho(cl,f) is not reported by Zvada 2010; assumed 0; see vignette Assumptions and deviations).
    # Block is positive definite (eigenvalues 0.0499 / 0.0381 / 0.00202).
    etalcl + etalmtt + etalfdepot ~ c(0.0362008, -0.0100994, 0.0089845, 0, 0.0130375, 0.0447783)

    # Zvada 2010 Table 4 IIV CLM/F = 23.9% (RSE 23.9%); omega^2 = log(1 + 0.239^2) = 0.0555492.
    # Single shared eta across CLM1 and CLM2 (same as the parent CL).
    etalcl_desrpt ~ 0.0555492

    # =========================================================================
    # Residual error (Zvada 2010 Table 2 / Table 4 'Residual variability' block).
    # Combined additive + proportional model on each observed concentration.
    # =========================================================================
    addSd <- 0.206
    label("Rifapentine additive residual SD (mg/L)")
    # Zvada 2010 Table 2 'Additive error (mg/liter) = 0.206' (RSE 12.3%).

    propSd <- 0.106
    label("Rifapentine proportional residual SD (fraction)")
    # Zvada 2010 Table 2 'Proportional error (%) = 10.6' (RSE 1.81%).

    addSd_desrpt <- 0.211
    label("25-DRFP additive residual SD (mg/L)")
    # Zvada 2010 Table 4 'Additive error (mg/liter) = 0.211' (RSE 9.7%).

    propSd_desrpt <- 0.191
    label("25-DRFP proportional residual SD (fraction)")
    # Zvada 2010 Table 4 'Proportional error (%) = 19.1' (RSE 0.7%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Time-after-dose used to switch CL between pre- and post-MTIME values.
    #    The paper's mpast(1) indicator is 0 before MTIME and 1 after; the
    #    rxode2 equivalent is tad(depot) (time after the most recent dose into
    #    the depot compartment), which resets at each crossover-occasion dose
    #    so the autoinduction step correctly fires at +43 h post each dose.
    # -----------------------------------------------------------------------
    tsd  <- tad(depot)
    t_switch  <- exp(lt_switch)
    t_switch_desrpt <- exp(lt_switch_desrpt)

    # -----------------------------------------------------------------------
    # 2. Individual rifapentine parameters. The shared etalcl applies to both
    #    pre- and post-MTIME clearances (Zvada Table 2 reports identical 19.2%
    #    IIV for the two rows, consistent with one shared eta).
    # -----------------------------------------------------------------------
    cl1 <- exp(lcl1 + etalcl)
    cl2 <- exp(lcl2 + etalcl)
    cl  <- cl1 * (tsd < t_switch) + cl2 * (tsd >= t_switch)
    vc  <- exp(lvc)
    ka  <- exp(lka)
    mtt <- exp(lmtt + etalmtt)
    nn  <- exp(lnn)

    # -----------------------------------------------------------------------
    # 3. Bioavailability with meal-effect covariates. TVF = theta_F * (1 + RXF)
    #    per Zvada 2010 Methods 'TVF = theta_F * (1 + RXF)' with theta_F fixed
    #    at 1 under fasting. The individual fdepot multiplies the typical-F
    #    expression by exp(etalfdepot).
    # -----------------------------------------------------------------------
    rxf <- e_meal_a_fdepot * MEAL_A +
           e_meal_b_fdepot * MEAL_B +
           e_meal_c_fdepot * MEAL_C +
           e_meal_d_fdepot * MEAL_D
    fdepot <- exp(lfdepot + etalfdepot) * (1 + rxf)

    # -----------------------------------------------------------------------
    # 4. Individual 25-DRFP parameters. Two-compartment metabolite disposition
    #    with the same step-function autoinduction structure as the parent.
    #    A single shared etalcl_desrpt applies to both CLM1 and CLM2 (Zvada
    #    Table 4 reports identical 23.9% IIV for the two rows).
    # -----------------------------------------------------------------------
    cl1_desrpt <- exp(lcl1_desrpt + etalcl_desrpt)
    cl2_desrpt <- exp(lcl2_desrpt + etalcl_desrpt)
    cl_desrpt  <- cl1_desrpt * (tsd < t_switch_desrpt) + cl2_desrpt * (tsd >= t_switch_desrpt)
    vc_desrpt  <- exp(lvc_desrpt)
    vp_desrpt  <- exp(lvp_desrpt)
    q_desrpt   <- exp(lq_desrpt)

    # -----------------------------------------------------------------------
    # 5. Micro-constants. The RFP elimination rate constant kel feeds the
    #    metabolite formation: rate of mass leaving RFP central = rate of mass
    #    entering 25-DRFP central. The paper does not apply an MW correction
    #    (RFP 877.03 g/mol vs 25-DRFP 835.0 g/mol, ~5% difference), so the
    #    metabolite apparent volume VMC/F implicitly absorbs that factor.
    # -----------------------------------------------------------------------
    kel       <- cl / vc
    kelm      <- cl_desrpt / vc_desrpt
    k_mc_mp   <- q_desrpt / vc_desrpt
    k_mp_mc   <- q_desrpt / vp_desrpt

    # -----------------------------------------------------------------------
    # 6. ODE system. Compartments: depot (transit-absorption sink); central
    #    (RFP plasma pool); central_desrpt + peripheral1_desrpt (25-DRFP). The
    #    rxode2 closed-form transit(nn, mtt, fdepot) returns the gamma-density
    #    input rate from the most recent dose; f(depot) <- 0 disables normal
    #    depot accumulation so transit() is the only depot input.
    # -----------------------------------------------------------------------
    d/dt(depot)              <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)            <-  ka * depot - kel * central
    d/dt(central_desrpt)     <-  kel * central - kelm * central_desrpt -
                                 k_mc_mp * central_desrpt + k_mp_mc * peripheral1_desrpt
    d/dt(peripheral1_desrpt) <-  k_mc_mp * central_desrpt - k_mp_mc * peripheral1_desrpt

    f(depot) <- 0

    # -----------------------------------------------------------------------
    # 7. Observations. Plasma concentrations are reported in mg/L; with dose
    #    in mg and volumes in L, central / vc has units of mg/L directly.
    # -----------------------------------------------------------------------
    Cc        <- central        / vc
    Cc_desrpt <- central_desrpt / vc_desrpt

    Cc        ~ add(addSd)            + prop(propSd)
    Cc_desrpt ~ add(addSd_desrpt)  + prop(propSd_desrpt)
  })
}
