Padari_2018_penicillin_G <- function() {
  description <- "Two-compartment IV population PK model for penicillin G (benzylpenicillin) in preterm and term neonates (Padari 2018; pooled with Metsvaht 2007 GA <=28 wk cohort). CL and Q are allometrically scaled to body weight (fixed exponent 0.75) with a fixed Rhodin-style postmenstrual-age (PMA) sigmoidal renal-maturation function on CL; Vc and Vp are allometrically scaled (fixed exponent 1.0)."
  reference   <- "Padari H, Metsvaht T, Germovsek E, Barker CI, Kipper K, Herodes K, Standing JF, Oselin K, Tasa T, Soeorg H, Lutsar I. Pharmacokinetics of penicillin G in preterm and term neonates. Antimicrob Agents Chemother. 2018;62(5):e02238-17. doi:10.1128/AAC.02238-17"
  vignette    <- "Padari_2018_penicillin_G"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight on the PK sampling day",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used for allometric scaling on CL and Q (fixed exponent 0.75) and on Vc and Vp (fixed exponent 1) with reference 70 kg per Padari 2018 Materials and Methods 'PK analyses' (Germovsek 2017 recommendation).",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the fixed Rhodin-style GFR-maturation Hill function on CL (Tmat50 = 47.7 weeks PMA, Hill = 3.4) per Padari 2018 Materials and Methods reference 48 (Rhodin et al. 2009). The canonical PAGE unit in inst/references/covariate-columns.md is months; this model uses weeks because the Rhodin function is defined in weeks of PMA. Convert from months as PAGE_weeks = PAGE_months * 4.35 if the upstream cohort is in months.",
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_studies      = 2L,
    age_range      = "PNA 0-3 days at PK sampling; GA 24-42 weeks (pooled across both source studies)",
    age_median     = "PMA 32.3 weeks (popPK cohort median; Padari 2018 Results 'popPK analysis')",
    weight_range   = "Pooled cohort spans ~0.5 to ~4 kg (Padari 2018 current cohort 2.0-3.8 kg; Metsvaht 2007 GA <=28 wk cohort included for the popPK fit)",
    weight_median  = "1.28 kg (popPK cohort median; Padari 2018 Results 'popPK analysis')",
    sex_female_pct = 35.3,
    race_ethnicity = "Not reported (single-centre Tartu University Hospital cohort plus Metsvaht 2007 historical cohort)",
    disease_state  = "Neonates of GA >=32 weeks (current study, n = 17) pooled with neonates of GA <=28 weeks from Metsvaht 2007 (n = 18), all treated for confirmed or suspected early-onset sepsis (EOS), congenital pneumonia, suspected congenital infection, or meconium aspiration syndrome",
    dose_range     = "Penicillin G 25,000 IU/kg or 50,000 IU/kg every 12 h as a 3-minute IV infusion (1 IU = 0.6 mg; equivalent to 15 mg/kg or 30 mg/kg)",
    regions        = "Estonia (Tartu University Hospital)",
    gestational_age_range = "24-42 weeks GA at birth (pooled across cohorts)",
    postmenstrual_age_range = "Equivalent to GA at the PK sampling day given PNA < 4 days in all subjects",
    samples_plasma = "Sparse sampling at trough, 5 min, 1 h, 3 h, 8 h, and 12 h after the steady-state dose (>= 36 h of therapy, typically after the 5th dose)",
    notes          = "Sex split (12 of 35 female = 35.3%) inferred from Padari 2018 Table 1 demographics for the current cohort (4 + 7 = 11 male of 17 enrolled); the Metsvaht 2007 cohort sex split is not retabulated. Concomitant gentamicin 4 mg/kg q24h was administered to all current-cohort subjects; no other potentially nephrotoxic drugs on the PK sampling day. Baseline median serum creatinine 52-61 umol/L, albumin 31-32 g/L, bilirubin 131-156 umol/L (Padari 2018 Table 1). Final popPK model retained no other covariates (BW, GA, SCR, ventilation, CPAP all non-significant at p < 0.01)."
  )

  ini({
    # Structural parameters (Padari 2018 Table 3 final-model column;
    # population estimates standardised to 70 kg using fixed allometric
    # exponents 0.75 on CL/Q and 1 on Vc/Vp per the Germovsek 2017
    # recommendation cited in Materials and Methods 'PK analyses').
    lcl <- log(13.2); label("Clearance standardised to 70 kg (L/h/70 kg)")                       # Padari 2018 Table 3: CL = 13.2 L/h/70 kg
    lvc <- log(10.3); label("Central volume of distribution standardised to 70 kg (L/70 kg)")     # Padari 2018 Table 3: V1 = 10.3 L/70 kg
    lq  <- log(55.6); label("Intercompartmental clearance standardised to 70 kg (L/h/70 kg)")    # Padari 2018 Table 3: Q  = 55.6 L/h/70 kg
    lvp <- log(29.8); label("Peripheral volume of distribution standardised to 70 kg (L/70 kg)") # Padari 2018 Table 3: V2 = 29.8 L/70 kg

    # Allometric exponents (fixed at the canonical Germovsek 2017 /
    # Anderson-Holford small-molecule renally-eliminated-drug pattern;
    # Padari 2018 Materials and Methods 'PK analyses').
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless)")          # Padari 2018 Methods (Germovsek 2017 recommendation)
    e_wt_vc_vp <- fixed(1.00); label("Shared allometric exponent on Vc and Vp (unitless)")         # Padari 2018 Methods (Germovsek 2017 recommendation)

    # Rhodin renal-maturation Hill function. Padari 2018 cites reference
    # 48 (Rhodin et al. 2009 Pediatr Nephrol 24:67-76) for the sigmoid
    # PMA-dependent renal-maturation function; the Rhodin 2009 GFR
    # parameters Tmat50 = 47.7 weeks and Hill = 3.4 are the de facto
    # standard and are also used in Germovsek 2018 meropenem in this
    # registry. Verification: the typical-individual prediction at the
    # popPK cohort median (WT = 1.28 kg, PMA = 32.3 weeks) with these
    # fixed maturation parameters reproduces Padari 2018 Results
    # 'popPK analysis' (CL = 0.15 L/h, Vc = 0.19 L, Q = 2.76 L/h, Vp = 0.54 L).
    tmat50   <- fixed(47.7); label("PMA at 50% renal maturation (weeks)")                        # Padari 2018 Methods reference 48 (Rhodin 2009)
    hill_mat <- fixed(3.4);  label("Hill coefficient for renal maturation (unitless)")           # Padari 2018 Methods reference 48 (Rhodin 2009)

    # Inter-individual variability. Padari 2018 Table 3 reports CV (%)
    # without an OMEGA covariance block, so the IIVs are encoded as
    # independent log-normal variances using omega^2 = log(CV^2 + 1):
    #   CL: CV 39% -> omega^2 = log(1 + 0.39^2) = 0.14164
    #   V1: CV 23% -> omega^2 = log(1 + 0.23^2) = 0.05154
    #   V2: CV 35% -> omega^2 = log(1 + 0.35^2) = 0.11556
    # Q has no CV / shrinkage entry in Table 3 and therefore no IIV.
    etalcl ~ 0.14164  # Padari 2018 Table 3: CV 39%; eta shrinkage 2.00%
    etalvc ~ 0.05154  # Padari 2018 Table 3: CV 23%; eta shrinkage 55.1% (V1)
    etalvp ~ 0.11556  # Padari 2018 Table 3: CV 35%; eta shrinkage 23.8% (V2)

    # Residual error. Padari 2018 Table 3 footnote reports a combined
    # proportional + additive model with proportional residual error
    # 13% and additive residual error 0.278 mg/L on the linear scale.
    propSd <- 0.13;  label("Proportional residual error (fraction)")                              # Padari 2018 Table 3 footnote
    addSd  <- 0.278; label("Additive residual error (mg/L)")                                      # Padari 2018 Table 3 footnote
  })

  model({
    # Rhodin renal-function maturation Hill function on PMA.
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # Individual PK parameters with allometric size scaling (reference
    # 70 kg) and the PMA-driven renal-maturation factor on CL.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants for the two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation. Dose in mg, V in L -> concentrations in mg/L
    # (equivalent to ug/mL). Penicillin G dose conversion: 1 IU = 0.6 mg
    # (Padari 2018 Materials and Methods 'Study drug administration').
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
