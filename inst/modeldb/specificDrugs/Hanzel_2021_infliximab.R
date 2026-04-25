Hanzel_2021_infliximab <- function() {
  description <- "Two-compartment population PK model of subcutaneous and intravenous infliximab CT-P13 (biosimilar) in adults with Crohn's disease and ulcerative colitis (Hanzel 2021)"
  reference <- "Hanzel J, Bukkems LH, Gecse KB, D'Haens GR, Mathot RAA. Population pharmacokinetics of subcutaneous infliximab CT-P13 in Crohn's disease and ulcerative colitis. Aliment Pharmacol Ther. 2021;54(10):1309-1319. doi:10.1111/apt.16609"
  vignette <- "Hanzel_2021_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, Vc, Vp, Q with reference 70 kg per Hanzel 2021 Table 3 (final model). Time-varying weight per Table 1 (covariate-parameter relations).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin (time-varying)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL with reference 44 g/L per Hanzel 2021 Table 3 (final model). Note unit is g/L (SI convention) -- distinct from g/dL used by some other infliximab popPK papers (e.g., Fasanmade 2009).",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (time-varying)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Multiplicative power-of-coefficient effect on CL: CL = CL_typ * theta_ATI^ADA_POS, with theta_ATI = 1.39 per Hanzel 2021 Table 3 (final model). Source paper labels this covariate 'ATI' (antibodies to infliximab); renamed to canonical ADA_POS per covariate-columns.md. Time-varying per Table 1 (covariate-parameter relations).",
      source_name        = "ATI"
    )
  )

  population <- list(
    n_subjects     = 175L,
    n_studies      = 1L,
    age_range      = "18-70 years (adults)",
    age_median     = "36 years",
    weight_range   = "43-118 kg",
    weight_median  = "69 kg",
    sex_female_pct = 54,
    race_ethnicity = "Not reported.",
    disease_state  = "Crohn's disease (55%) or ulcerative colitis (45%); moderate-to-severe activity at baseline.",
    dose_range     = "5 mg/kg IV induction at weeks 0 and 2 (and additional IV doses through week 22 for the part-2 IV arm), followed by subcutaneous maintenance doses of 120, 180 or 240 mg every 2 weeks (weight-based 120 mg if <=80 kg, 240 mg if >80 kg in part 2).",
    regions        = "Single multinational CT-P13 1.6 study (NCT02883452).",
    immunomodulators = "84% on immunomodulators pre-treatment, 46% during treatment.",
    corticosteroids  = "61% on corticosteroids pre-treatment, 41% during treatment.",
    ada_incidence    = "33% developed neutralising antibodies during treatment.",
    notes          = "Baseline demographics from Hanzel 2021 Table 2 (n = 175). 2772 PK samples used after exclusions of below-LLOQ and ADA-driven-low samples (4.4% below LLOQ before exclusions). Reference covariate values for the typical patient: WT = 70 kg, ALB = 44 g/L, ADA-negative."
  )

  ini({
    # Structural parameters -- typical values for the reference patient
    # (70 kg, ALB 44 g/L, ADA-negative) per Hanzel 2021 Table 3 (Final model column).
    lka     <- log(0.273); label("First-order SC absorption rate constant (Ka, 1/day)")            # Hanzel 2021 Table 3: Ka = 0.273 /d
    lcl     <- log(0.355); label("Clearance for the reference patient (CL, L/day)")                # Hanzel 2021 Table 3: CL = 0.355 L/d
    lvc     <- log(3.10);  label("Central volume of distribution for the reference patient (Vc, L)")     # Hanzel 2021 Table 3: Vc = 3.10 L
    lvp     <- log(1.93);  label("Peripheral volume of distribution for the reference patient (Vp, L)")  # Hanzel 2021 Table 3: Vp = 1.93 L
    lq      <- log(0.598); label("Inter-compartmental clearance for the reference patient (Q, L/day)")   # Hanzel 2021 Table 3: Q = 0.598 L/d
    lfdepot <- log(0.791); label("Subcutaneous bioavailability (F1, fraction)")                          # Hanzel 2021 Table 3: F1 = 79.1%

    # Covariate effect parameters from Hanzel 2021 Table 3 (Final model column).
    # Continuous covariates use a power model relative to the reference value;
    # the categorical ATI effect uses theta_ATI^ADA_POS (a power function with the
    # 0/1 indicator as the on-off switch, per the Methods section).
    e_wt_cl  <-  0.666; label("Power exponent of body weight on CL ((WT/70)^e_wt_cl)")                # Hanzel 2021 Table 3
    e_wt_vc  <-  0.385; label("Power exponent of body weight on Vc ((WT/70)^e_wt_vc)")                # Hanzel 2021 Table 3
    e_wt_vp  <-  1.08;  label("Power exponent of body weight on Vp ((WT/70)^e_wt_vp)")                # Hanzel 2021 Table 3
    e_wt_q   <-  1.26;  label("Power exponent of body weight on Q ((WT/70)^e_wt_q)")                  # Hanzel 2021 Table 3
    e_alb_cl <- -0.826; label("Power exponent of serum albumin on CL ((ALB/44)^e_alb_cl)")            # Hanzel 2021 Table 3
    e_ada_cl <-  1.39;  label("ATI multiplicative effect on CL (e_ada_cl^ADA_POS; +39% when ADA-positive)") # Hanzel 2021 Table 3

    # Inter-individual variability on CL, F1, Vc, Ka -- modelled as a 4x4 block.
    # Variances are omega^2 on log-scale: omega^2 = log(1 + CV^2). Reported %CV
    # values are: CL 27.7%, F1 16.4%, Vc 21.4%, Ka 48.5% (Hanzel 2021 Table 3).
    # Off-diagonals are correlations from Table 3 ("Correlation between ..."),
    # converted to covariances via cov_ij = corr_ij * sqrt(var_i * var_j).
    #   var_cl     = log(1 + 0.277^2) = 0.073928
    #   var_f      = log(1 + 0.164^2) = 0.026541
    #   var_vc     = log(1 + 0.214^2) = 0.044778
    #   var_ka     = log(1 + 0.485^2) = 0.211253
    #   cov(CL,F)  = -0.013    * sqrt(var_cl * var_f)  = -0.000576
    #   cov(CL,Vc) =  0.028    * sqrt(var_cl * var_vc) =  0.001611
    #   cov(CL,Ka) = -0.046    * sqrt(var_cl * var_ka) = -0.005749
    #   cov(F,Vc)  =  0.00008  * sqrt(var_f  * var_vc) =  2.76e-6
    #   cov(F,Ka)  =  0.003    * sqrt(var_f  * var_ka) =  0.000225
    #   cov(Vc,Ka) = -0.069    * sqrt(var_vc * var_ka) = -0.006711
    etalcl + etalfdepot + etalvc + etalka ~
      c(0.073928,
        -0.000576, 0.026541,
         0.001611, 2.76e-6,   0.044778,
        -0.005749, 0.000225, -0.006711, 0.211253)   # Hanzel 2021 Table 3

    # Residual error model -- combined additive (mg/L) + proportional from
    # Hanzel 2021 Table 3 (Final model column).
    addSd  <- 1.66;  label("Additive residual error (mg/L)")              # Hanzel 2021 Table 3
    propSd <- 0.102; label("Proportional residual error (fraction)")      # Hanzel 2021 Table 3
  })
  model({
    # Individual PK parameters. Reference patient: 70 kg, albumin 44 g/L, ADA-negative.
    # Covariate forms per Hanzel 2021 Table 3 / Methods (full fixed-effect modelling):
    #   CL = CL_typ * (WT/70)^e_wt_cl  * (ALB/44)^e_alb_cl * e_ada_cl^ADA_POS
    #   Vc = Vc_typ * (WT/70)^e_wt_vc
    #   Vp = Vp_typ * (WT/70)^e_wt_vp
    #   Q  = Q_typ  * (WT/70)^e_wt_q
    #   Ka = Ka_typ
    #   F1 = F1_typ                                          (applied at depot for SC dosing)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (ALB / 44)^e_alb_cl *
      e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp)          * (WT / 70)^e_wt_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Bioavailability is applied to depot dosing (SC); IV doses go directly to
    # central with implied F = 1, consistent with the source's joint IV+SC analysis.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Concentration: dose in mg, volume in L -> mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
