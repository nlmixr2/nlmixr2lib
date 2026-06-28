Padari_2012_meropenem <- function() {
  description <- "One-compartment IV population PK model for meropenem in very-low-birth-weight neonates (gestational age <=32 weeks, birth weight <1,500 g; n=19; Padari 2012). Vss scales linearly with current body weight; CL follows the Rhodin (2009) fixed renal-maturation function (allometric exponent 0.75 on CL, Hill-type postmenstrual-age maturation with TM50 = 47.7 weeks and Hill = 3.4); serum creatinine, postnatal age, and gestational age were screened and did not improve fit and are not retained."
  reference <- "Padari H, Metsvaht T, Korgvee LT, Germovsek E, Ilmoja ML, Kipper K, Herodes K, Standing JF, Oselin K, Lutsar I. Short versus long infusion of meropenem in very-low-birth-weight neonates. Antimicrob Agents Chemother. 2012;56(9):4760-4764. doi:10.1128/AAC.00655-12"
  vignette <- "Padari_2012_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight at PK sampling",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Linear scaling on Vss (V = TVV/70 * WT) and allometric scaling on CL (Rhodin 2009 fixed exponent 0.75 with reference 70 kg) per Padari 2012 Methods 'Statistical and PK analyses' / Results 'Population modeling'.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the fixed Rhodin (2009) renal-maturation Hill function on CL (TM50 = 47.7 weeks, Hill = 3.4) per Padari 2012 Results 'Population modeling' ('the fixed Rhodin model (25) was used'). Padari 2012 reports neonates with PMA approximately equal to GA + PNA; pooled cohort mean PMA approximately 29 weeks (Table 1: GA 26.9 / 25.8 weeks + PNA 15.6 / 20.5 days at PK sampling).",
      source_name        = "PMA"
    )
  )

  covariatesDataExcluded <- list(
    SCR = list(
      description = "Serum creatinine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened in covariate analysis (corrected for postnatal age per Ceriotti 2008 reference SCR; Padari 2012 Methods 'Statistical and PK analyses') but did not significantly improve model fit and was not retained in the final model (Padari 2012 Results 'Population modeling')."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = "Screened as a covariate on CL but no significant correlation with meropenem clearance was observed; not retained in the final model (Padari 2012 Results 'Noncompartmental PK analysis' and 'Population modeling')."
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened as a covariate but did not improve model fit and was not retained in the final model (Padari 2012 Results 'Population modeling'). GA enters the retained model only indirectly via PAGE = GA + PNA."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19L,
    n_studies      = 1L,
    age_range      = "PNA <=56 days at enrolment; PMA approximately 28-32 weeks (gestational + postnatal at PK sampling)",
    age_median     = "PNA 15.6 days (group 1, short infusion) / 20.5 days (group 2, prolonged infusion) at PK sampling",
    weight_range   = "Birth weight 0.84-0.90 kg (group means; Padari 2012 Table 1); current body weight 0.97-0.98 kg at PK sampling",
    weight_median  = "0.977 kg (pooled current body weight at PK sampling, Padari 2012 Table 1)",
    sex_female_pct = 36.8,
    race_ethnicity = "Not reported (single-centre Estonian neonatal ICU cohort)",
    disease_state  = "Very-low-birth-weight preterm neonates (GA <=32 weeks at birth, BW <1,500 g) receiving meropenem treatment for sepsis (n=16), pneumonia (n=3), or necrotizing enterocolitis due to proven or highly suspected resistant pathogens, or for clinical deterioration on empirical antibiotics. 95% required respiratory support; 42% required vasoactive treatment. Excluded major uncorrected congenital malformations and infants expected to die within 24 h.",
    dose_range     = "Meropenem 20 mg/kg q12h IV; group 1 (n=9) administered over 30 min; group 2 (n=10) administered over 4 h (after at least two prolonged infusions to reach steady state). After PK sampling, all subjects returned to 30-min infusions.",
    regions        = "Estonia (Tartu University Hospital and Tallinn Children's Hospital NICUs)",
    gestational_age_range = "GA at birth 25.8-26.9 weeks (group means; Padari 2012 Table 1)",
    samples_plasma = "Approximately 6 plasma samples per subject (pre-dose and 0.5, 1.5, 4, 8, 12 h after the 4th-7th doses of meropenem)",
    notes          = "Demographics from Padari 2012 Table 1 (n=9 short / n=10 prolonged). Mean current body weight at PK sampling 0.984 kg (group 1) and 0.969 kg (group 2). Mean serum creatinine at enrolment 51.4 umol/L (group 1) and 44.8 umol/L (group 2). Concomitant vancomycin or ibuprofen in 4/3 subjects. 11/19 had positive blood cultures (coagulase-negative staphylococci most common). Trial registration EU CTR 2009-017823-24. Sex_female_pct = 7/19 = 36.8% (groups had 6 males each of 9 and 10 subjects, so 12 males / 7 females total)."
  )

  ini({
    # Structural PK (Padari 2012 Abstract / Results 'Population modeling').
    #
    # Padari 2012 reports the final estimates as Vss = 0.301 L/kg and
    # CL = 0.061 L/h/kg with the maturation function FIXED at Rhodin (2009).
    # Rhodin's published form is
    #     CL_i = CL_std * (WT_i/70)^0.75 * F_mat(PMA_i)
    # with F_mat = PMA^Hill / (PMA^Hill + TM50^Hill), TM50 = 47.7 wk, Hill = 3.4.
    # The paper does not publish the underlying NONMEM THETA(CL_std); we
    # back-calculate it from the reported 0.061 L/h/kg using WT = 1 kg
    # (rounded VLBW-neonate reference per the skill's "undefined
    # reference value -> rounded standard" rule) and PMA = 29 weeks
    # (pooled cohort mean GA + PNA at PK sampling from Table 1):
    #     F_mat(29) = 29^3.4 / (29^3.4 + 47.7^3.4) = 0.1555
    #     size(1)   = (1/70)^0.75                    = 0.04132
    #     TVCL_70   = 0.061 / (0.04132 * 0.1555)     = 9.49 L/h
    # so CL_i (1 kg, 29 wk PMA) = 9.49 * 0.04132 * 0.1555 = 0.061 L/h.
    # The volume back-calculation is direct: TVV_70 = 0.301 L/kg * 70 kg = 21.07 L.
    lcl <- log(9.49);  label("Plasma clearance standardised to 70 kg fully mature (L/h)")  # Padari 2012 Abstract / Results 'Population modeling': back-calc from CL = 0.061 L/h/kg with Rhodin fixed scaling
    lvc <- log(21.07); label("Central volume of distribution standardised to 70 kg (L)")    # Padari 2012 Abstract: Vss = 0.301 L/kg, scaled to 70 kg = 21.07 L

    # Rhodin (2009) fixed renal-maturation parameters per Padari 2012 Results
    # 'Population modeling' ("the fixed Rhodin model (25) was used"). Rhodin et al.
    # 2009 Pediatr Nephrol 24:67-76 reports TM50 = 47.7 wk PMA and Hill = 3.4 for
    # human renal clearance maturation; the allometric exponent on CL is fixed
    # at 0.75 (Rhodin's published convention).
    tmat50   <- fixed(47.7); label("PMA at 50% renal maturation (weeks; Rhodin 2009 fixed)")   # Padari 2012 ref (25): Rhodin 2009
    hill_mat <- fixed(3.4);  label("Hill coefficient for renal maturation (unitless; Rhodin 2009 fixed)")  # Padari 2012 ref (25): Rhodin 2009
    e_wt_cl  <- fixed(0.75); label("Allometric exponent on CL (unitless; Rhodin 2009 fixed)")  # Padari 2012 ref (25): Rhodin 2009 fixed allometric

    # Inter-individual variability. Padari 2012 reports a final model with a
    # prediction-corrected VPC (Figure 3) and goodness-of-fit plots (Figure 2)
    # but does not publish OMEGA estimates. Per the nlmixr2lib standing policy
    # for unreported variance components, IIVs are set to fixed(0). The
    # vignette therefore reproduces typical-individual predictions (deterministic
    # in covariates) and compares those to the published NCA group means in
    # Table 2 rather than to a stochastic VPC.
    etalcl ~ fixed(0)   # Padari 2012 does not report IIV on CL
    etalvc ~ fixed(0)   # Padari 2012 does not report IIV on V

    # Residual error. Padari 2012 does not report the residual-error model or
    # magnitude. Per the standing policy for unreported RUV, propSd is encoded
    # as a small placeholder (0.001 = 0.1% CV) so the model parses and simulates
    # cleanly while staying effectively deterministic. The vignette uses
    # rxode2::zeroRe() for typical-individual replication and so does not
    # depend on the placeholder magnitude.
    propSd <- fixed(0.001); label("Proportional residual SD placeholder (unreported in Padari 2012; see Errata)")  # Padari 2012 does not report residual error
  })

  model({
    # Rhodin (2009) renal-maturation Hill function on PMA.
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # Individual PK parameters. Vss scales linearly with current body weight;
    # CL follows Rhodin (2009): allometric 0.75 on WT plus Hill maturation on PMA.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
