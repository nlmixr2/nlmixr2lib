Jorga_2000_tolcapone_fluctuators <- function() {
  description <- "Two-compartment population PK model with first-order absorption (no lag) for tolcapone in parkinsonian patients with fluctuating levodopa response, with effects of lean body weight and serum protein on clearance, lean body weight and dose group on central volume, serum albumin and dose group on peripheral volume, and concomitant food on bioavailability (Jorga 2000, fluctuator dataset, n=215)"
  reference <- paste(
    "Jorga K, Fotteler B, Banken L, Snell P, Steimer JL.",
    "Population pharmacokinetics of tolcapone in parkinsonian patients in dose finding studies.",
    "Br J Clin Pharmacol. 2000;49(1):39-48.",
    "doi:10.1046/j.1365-2125.2000.00113.x"
  )
  vignette <- "Jorga_2000_tolcapone"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass (paper alias 'LBW', lean body weight; computed via James 1976 formula from total body weight and height)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on CL and Vc with reference 55 kg (population median, Table 1). Computed per James 1976: LBW(kg) = (1.10*BW - 128*BW^2/HT^2) for males, (1.07*BW - 148*BW^2/HT^2) for females, with BW in kg and HT in cm (Jorga 2000 Methods).",
      source_name        = "LBW"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on CL with reference 72 g/L (population median, Table 1). Source column 'Protein' in the paper.",
      source_name        = "Protein"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on Vp with reference 44 g/L (population median, Table 1). Source column 'Albumin' in the paper.",
      source_name        = "Albumin"
    ),
    DOSE_50MG = list(
      description        = "Subject-level indicator for the 50 mg three-times-daily tolcapone dose group (1 = 50 mg t.i.d. arm, 0 = 200 mg or 400 mg t.i.d. arm)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (200 mg or 400 mg t.i.d. arm; combined with DOSE_400MG = 0 selects the 200 mg reference)",
      notes              = "Paper indicator I_Dose50mg; subject-level (each fluctuator patient was randomized to one of 50/200/400 mg t.i.d.). Multiplicative effect (1 + e_dose_50mg_vc_vp * DOSE_50MG) applied to Vc and to Vp; e_dose_50mg_vc_vp = -0.45 means V is 55% of the 200 mg reference at the 50 mg dose. The dose-dependent V is an empirical finding the authors note is plausibly driven by a few high-V outliers in the small-volume cohort and was NOT confirmed in the nonfluctuator dataset (Jorga 2000 Discussion).",
      source_name        = "(derived from study-arm dose level)"
    ),
    DOSE_400MG = list(
      description        = "Subject-level indicator for the 400 mg three-times-daily tolcapone dose group (1 = 400 mg t.i.d. arm, 0 = 50 mg or 200 mg t.i.d. arm)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (50 mg or 200 mg t.i.d. arm; combined with DOSE_50MG = 0 selects the 200 mg reference)",
      notes              = "Paper indicator I_Dose400mg; subject-level. Multiplicative effect (1 + e_dose_400mg_vc_vp * DOSE_400MG) applied to Vc and to Vp; e_dose_400mg_vc_vp = +0.40 means V is 140% of the 200 mg reference at the 400 mg dose. See DOSE_50MG notes for the empirical caveat.",
      source_name        = "(derived from study-arm dose level)"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (1 = dose taken with concomitant food, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Paper indicator I_Food. Multiplicative effect (1 + e_food_f * FED) applied on F1; e_food_f = -0.12 corresponds to a ~12% reduction in relative bioavailability in the fed state for the fluctuator cohort (Jorga 2000 Table 3 theta_Food = 0.88; Discussion: 10-15% reduction in fluctuators). Reference fasted F1 fixed at 0.6.",
      source_name        = "Food"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 215L,
    n_studies      = 2L,
    study_names    = c("Fluctuator 50 mg t.i.d. (n = 75)",
                       "Fluctuator 200 mg t.i.d. (n = 74)",
                       "Fluctuator 400 mg t.i.d. (n = 66)"),
    n_observations = 981L,
    age_range      = "34-82 years",
    age_median     = "65 years",
    weight_range   = "36-153 kg",
    weight_median  = "71 kg",
    lbm_range      = "25-83 kg",
    lbm_median     = "55 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Caucasian = 98.0, Black = 0.24, Asian = 0.73, Other = 1.46),
    disease_state  = "Parkinson's disease with fluctuating motor response to levodopa/AADC inhibitor therapy ('fluctuators')",
    dose_range     = "50, 200, or 400 mg tolcapone three times daily for 6 weeks, in addition to ongoing levodopa-carbidopa (Sinemet) or levodopa-benserazide (Madopar) therapy",
    regions        = "Three multicentre Phase II dose-finding studies across 49 centres worldwide",
    sampling       = "Sparse: 5-8 plasma samples per patient on 2-5 occasions across study days 14, 21/28, and 42; samples taken pre-dose, near Cmax, and during the decline phase",
    notes          = "Demographics aggregated for the combined fluctuator + nonfluctuator population in Table 1; n=412 enrolled, n=275 with PK data (215 fluctuators + 60 nonfluctuators), and the demographics here are taken from the combined-cohort 'Fluctuators' column of Table 1 where reported (n=315 for the broader fluctuator-arm denominator)."
  )

  ini({
    # Structural PK parameters -- typical values for a 200 mg t.i.d. fluctuator
    # at the population-median LBM (55 kg), serum protein (72 g/L), serum
    # albumin (44 g/L), under fasted conditions. From Jorga 2000 Table 3
    # (Fluctuator model, Final estimate column).
    lka      <- log(1.7);  label("Absorption rate constant (ka, 1/h)")                        # Table 3 (ka = 1.7 /h)
    lcl      <- log(4.8);  label("Apparent clearance for the 200 mg reference (CL/F1, L/h)")  # Table 3 (CL = 4.8 L/h)
    lvc      <- log(16);   label("Central volume of distribution for the 200 mg reference (Vc, L)")    # Table 3 (Vc = 16 L)
    lvp      <- log(12);   label("Peripheral volume of distribution for the 200 mg reference (Vp, L)") # Table 3 (Vp = 12 L)
    lq       <- log(5.2);  label("Inter-compartmental clearance (Q, L/h)")                    # Table 3 (Q = 5.2 L/h)
    # Fasted absolute bioavailability fixed at 0.6 from upstream IV/PO single-
    # dose study (Jorga et al. 1998, Eur J Clin Pharmacol 54:443-447; ref [22]).
    lfdepot  <- fixed(log(0.6)); label("Fasted absolute bioavailability (F1, fraction)")     # Methods: F1 fixed to 0.6 [ref 22]

    # Covariate effects on CL: power-form (LBM/55)^e_lbm_cl and (TPRO/72)^e_tpro_cl
    e_lbm_cl    <-  0.73; label("LBM power exponent on CL")              # Table 3 (theta_LBW(CL) = 0.73)
    e_tpro_cl   <- -0.81; label("Total serum protein power exponent on CL")   # Table 3 (theta_Protein(CL) = -0.81)

    # Covariate effects on Vc: power-form (LBM/55)^e_lbm_vc, plus dose-group multipliers
    e_lbm_vc    <-  0.65; label("LBM power exponent on Vc")              # Table 3 (theta_LBW(Vc) = 0.65)

    # Covariate effects on Vp: power-form (ALB/44)^e_alb_vp, plus dose-group multipliers
    e_alb_vp    <-  2.82; label("Albumin power exponent on Vp")          # Table 3 (theta_Albumin(Vp) = 2.82)

    # Dose-group multiplicative effects on Vc and Vp (shared; paper estimates a
    # single theta_Dose50mg(V) and theta_Dose400mg(V) applied identically to Vc
    # and Vp). Encoding (1 + e_dose_NNN_vc_vp * DOSE_NNN_MG): paper theta = 0.55
    # for 50 mg group -> e = -0.45 (V is 55% of 200 mg reference); paper theta
    # = 1.40 for 400 mg group -> e = +0.40 (V is 140% of 200 mg reference).
    e_dose_50mg_vc_vp  <- -0.45; label("Relative change in Vc and Vp for the 50 mg dose group vs 200 mg reference (fraction)")   # Table 3 (theta_Dose50mg(V) = 0.55 -> 0.55 - 1)
    e_dose_400mg_vc_vp <-  0.40; label("Relative change in Vc and Vp for the 400 mg dose group vs 200 mg reference (fraction)")  # Table 3 (theta_Dose400mg(V) = 1.40 -> 1.40 - 1)

    # Food effect on F1: (1 + e_food_f * FED) multiplier on fasted F1 = 0.6.
    # Paper theta_Food = 0.88 -> e_food_f = -0.12.
    e_food_f    <- -0.12; label("Relative change in F1 for fed vs fasted dosing (fraction)")  # Table 3 (theta_Food(F) = 0.88 -> 0.88 - 1)

    # Inter-individual variability (paper reports omega^2 on the log-normal
    # internal scale; sqrt(exp(omega^2) - 1) reproduces the parenthetical CV%).
    # omega^2 = 0.08 -> CV ~= 29% (Table 3)
    # omega^2 = 0.42 -> CV ~= 72% (Table 3)
    # omega^2 = 0.28 -> CV ~= 57% (Table 3)
    etalcl ~ 0.08
    etalvc ~ 0.42
    etalvp ~ 0.28

    # Residual error -- log-normal multiplicative (paper equation
    # DV = CP * exp(eps_mult) + eps_add with eps_add absent for the fluctuator
    # model; Table 3 lists only the multiplicative variance). For
    # |propSd| ~ 0.47 the log-normal multiplicative and linear-proportional
    # forms differ but the linear-proportional encoding is the convention used
    # across nlmixr2lib. omega^2(eps_mult) = 0.22 -> SD = sqrt(0.22) ~= 0.469.
    propSd <- 0.469; label("Proportional residual error (fraction)")     # Table 3 (eps_mult variance = 0.22 -> SD = sqrt(0.22))
  })

  model({
    # Individual PK parameters with covariate effects.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
          (LBM  / 55)^e_lbm_cl *
          (TPRO / 72)^e_tpro_cl
    vc <- exp(lvc + etalvc) *
          (LBM  / 55)^e_lbm_vc *
          (1 + e_dose_50mg_vc_vp  * DOSE_50MG) *
          (1 + e_dose_400mg_vc_vp * DOSE_400MG)
    vp <- exp(lvp + etalvp) *
          (ALB  / 44)^e_alb_vp *
          (1 + e_dose_50mg_vc_vp  * DOSE_50MG) *
          (1 + e_dose_400mg_vc_vp * DOSE_400MG)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment first-order absorption -- no lag time in the fluctuator
    # model per Jorga 2000 Results "Final estimates" paragraph
    # ("in contrast to the nonfluctuator model, the t_lag was excluded from
    # the fluctuator model").
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Fasted bioavailability with multiplicative food effect; applied to the
    # depot compartment (oral dosing).
    f(depot) <- exp(lfdepot) * (1 + e_food_f * FED)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
