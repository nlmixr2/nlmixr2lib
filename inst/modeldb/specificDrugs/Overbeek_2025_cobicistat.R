Overbeek_2025_cobicistat <- function() {
  description <- "Well-stirred liver model for oral cobicistat (CYP3A pharmacokinetic booster) pooling healthy volunteers, postpartum women with HIV, patients with rheumatoid arthritis and patients with solid tumours, with Erlang-type absorption through three transit compartments, a mechanistic hepatic-extraction central/liver disposition driven by unbound intrinsic clearance per litre of liver, a priori allometric scaling to 70 kg, and a higher intrinsic clearance in the PROACTIVE (olaparib-boosting) cohort (Overbeek 2025)"
  reference   <- "Overbeek JK, van Erp NP, Burger DM, den Broeder AA, Koolen SLW, Huitema ADR, ter Heine R. Population Pharmacokinetics of Cobicistat and its Effect on the Pharmacokinetics of the Anticancer Drug Olaparib. Clin Pharmacokinet. 2025;64(3):425-435. doi:10.1007/s40262-025-01480-w"
  vignette    <- "Overbeek_2025_cobicistat_olaparib_boosting"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters in three distinct ways. (1) A priori allometric scaling to a standardised total body weight of 70 kg with exponents 0.75 for flow (hepatic plasma flow), 1 for volume (Vc) and -0.25 for the absorption rate constant ktr (Overbeek 2025 Methods 2.3, final paragraph before Eq 5). (2) The liver volume that scales intrinsic clearance, VL = 0.10 * TBW^0.59 (Eq 5), which is a function of the raw weight in kg and is deliberately NOT re-normalised to 70 kg. (3) Because CLint = theta_CLint * VL (Eq 6), the weight dependence of clearance is carried entirely by VL and no separate allometric exponent is applied to CLint. Cohort median 74 kg, range 52-124 kg (Table 1).",
      source_name        = "WEIGHT"
    ),
    STUDY_PROACTIVE = list(
      description        = "PROACTIVE study (olaparib-boosting cohort) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled DATE-4, PANNA and PRACTICAL cohorts)",
      notes              = "1 = subject enrolled in the PROACTIVE trial (NCT05078671; patients with solid tumours receiving cobicistat 150 mg twice daily to boost olaparib); 0 = subject enrolled in DATE-4, PANNA or PRACTICAL (all cobicistat 150 mg once daily). Multiplicative power-form effect on cobicistat intrinsic clearance per Overbeek 2025 Eq 7 (P = P0 * theta^STUDY), the only one of the four candidate study covariates retained in the final model (Results 3.1). The paper attributes the 1.21-fold higher intrinsic clearance to selection of a cohort with relatively high CYP3A activity rather than to olaparib itself, and the indicator is confounded with both the twice-daily regimen and the solid-tumour population, so it should be read as a cohort effect and not as an olaparib drug-drug interaction (Discussion, paragraphs 3 and 4). The NONMEM data item is named OLAP in the control stream (Online Resource Material 1, $INPUT), which is why Table 2 labels the coefficient CLint-olaparib.",
      source_name        = "OLAP"
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "cobicistat", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "cobicistat", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "cobicistat", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "cobicistat", units = "mg", specimen = "administration site", verified = TRUE),
    liver    = list(analyte = "cobicistat", units = "mg", specimen = "tissue", verified = TRUE),
    central  = list(analyte = "cobicistat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 66L,
    n_studies      = 4L,
    n_observations = 683L,
    age_range      = "21-78 years (median 51.5)",
    age_median     = "51.5 years",
    weight_range   = "52-124 kg (median 74)",
    weight_median  = "74 kg",
    sex_female_pct = 63.6,
    disease_state  = "Pooled healthy volunteers (DATE-4, n = 16), postpartum women living with HIV (PANNA, n = 12), patients with rheumatoid arthritis (PRACTICAL, n = 26) and patients with solid tumours (PROACTIVE, n = 12)",
    dose_range     = "Cobicistat 150 mg orally once daily (DATE-4, PANNA, PRACTICAL) or 150 mg twice daily (PROACTIVE), all at steady state after at least 7 days of treatment",
    regions        = "The Netherlands (Radboud university medical center and collaborating sites)",
    co_medication  = "Cobicistat was co-administered as a booster with atazanavir (DATE-4), elvitegravir (PANNA), tofacitinib (PRACTICAL) or olaparib (PROACTIVE)",
    notes          = "Demographics from Overbeek 2025 Table 1. Sex is reported as 24 of 66 male (36%), so 63.6% female. All studies used dense PK sampling over one 12 h or 24 h dosing interval at steady state. Only the postpartum (non-pregnant) PANNA occasion was included; the third-trimester data were excluded because cobicistat PK was markedly different in pregnancy (Methods 2.1). 14 of 683 samples (2.1%) were below the 0.03 mg/L limit of quantification and were retained using the 'All data' (M1) method."
  )

  ini({
    # ---- Structural parameters, at the 70 kg reference weight -------------
    # (Overbeek 2025 Table 2; NONMEM $THETA of Online Resource Material 1.)
    lktr   <- log(3.92)  ; label("Erlang absorption transition rate constant (ktr, 1/h)")                              # Overbeek 2025 Table 2: k_tr = 3.92 /h (RSE 12.3%); Online Resource Material 1 $THETA 1
    lvc    <- log(69.7)  ; label("Apparent central volume of distribution (Vc/F, L)")                                  # Overbeek 2025 Table 2: V_c = 69.7 L (RSE 5.6%); Online Resource Material 1 $THETA 2
    lclint <- log(322)   ; label("Apparent unbound intrinsic clearance per litre of liver (CLint/F, L/h/L liver)")     # Overbeek 2025 Table 2: CL_int = 322 L/h/L_liver (RSE 6.7%); Online Resource Material 1 $THETA 3

    # Prehepatic bioavailability. The paper reports all volume and clearance
    # parameters relative to the unknown absolute bioavailability (Methods 2.3),
    # so F1 is an anchor rather than an estimate; hepatic first pass is carried
    # mechanistically by routing the absorbed dose through the liver compartment.
    lfdepot <- fixed(log(1)) ; label("Prehepatic bioavailability (F1, unitless)")                                      # Online Resource Material 1 $PK: F1=1*EXP(ETA(1)) with $OMEGA 1 = 0 FIX

    # ---- Study covariate on intrinsic clearance (Eq 7, power form) --------
    e_study_proactive_clint <- 1.21 ; label("Ratio of intrinsic clearance in the PROACTIVE cohort to the other cohorts (unitless)")  # Overbeek 2025 Table 2: CL_int-olaparib = 1.21 (RSE 16.3%), p < 0.001; Online Resource Material 1 $THETA 4

    # ---- Physiological constants of the well-stirred liver model ----------
    # Assumed rather than estimated (Methods 2.3). hct was fixed because
    # haematocrit was not measured in all four studies; the paper's sensitivity
    # analysis over 0.30-0.50 changed the estimates negligibly (Results 3.1).
    q_liver      <- fixed(90)    ; label("Hepatic blood flow at the 70 kg reference weight (QH, L/h)")                 # Overbeek 2025 Methods 2.3: "assuming a hepatic blood flood (QH) of 90 L/h"
    hct          <- fixed(0.44)  ; label("Haematocrit, assumed (Ht, unitless)")                                        # Overbeek 2025 Methods 2.3: "hematocrit (Ht) of 0.44"; sensitivity analysis over 0.30-0.50 in Results 3.1
    fu           <- fixed(0.025) ; label("Fraction unbound in plasma, literature value (fu, unitless)")                # Overbeek 2025 Methods 2.3: "an unbound fraction in plasma (fu) of 0.025 for cobicistat"
    v_liver_coef <- fixed(0.10)  ; label("Coefficient of the liver-volume relation, from Small 2017 (L/kg^0.59)")    # Overbeek 2025 Eq 5: VL = 0.10 * TBW^0.59, citing reference [31] (Small BG et al., Biopharm Drug Dispos 2017;38(4):290-300)
    e_wt_v_liver <- fixed(0.59)  ; label("Body-weight exponent of the liver-volume relation, from Small 2017 (unitless)")  # Overbeek 2025 Eq 5: VL = 0.10 * TBW^0.59, citing reference [31] (Small BG et al., Biopharm Drug Dispos 2017;38(4):290-300)

    # ---- A priori allometric exponents -----------------------------------
    e_wt_fq  <- fixed(0.75)  ; label("Allometric exponent on hepatic plasma flow (unitless)")                          # Overbeek 2025 Methods 2.3: allometric components of 0.75 for flow parameters; Online Resource Material 1 $PK ALLOQHP
    e_wt_vc  <- fixed(1)     ; label("Allometric exponent on central volume (unitless)")                               # Overbeek 2025 Methods 2.3: allometric component of 1 for volume parameters; Online Resource Material 1 $PK ALLOV
    e_wt_ktr <- fixed(-0.25) ; label("Allometric exponent on the transition rate constant (unitless)")                 # Overbeek 2025 Methods 2.3: allometric component of -0.25 for absorption parameters; Online Resource Material 1 $PK ALLOKTR

    # ---- Inter-individual variability ------------------------------------
    # Log-normal (Eq 4). Table 2 reports these as %CV via sqrt(exp(omega^2) - 1)
    # (Table 2 footnote), so the variances below are the NONMEM $OMEGA values
    # and reproduce the printed %CV exactly.
    etalktr   ~ 0.488   # Overbeek 2025 Online Resource Material 1 $OMEGA 2; Table 2 k_tr BSV 79.3% = sqrt(exp(0.488) - 1) * 100
    etalvc    ~ 0.129   # Overbeek 2025 Online Resource Material 1 $OMEGA 3; Table 2 V_c BSV 37.1% = sqrt(exp(0.129) - 1) * 100
    etalclint ~ 0.208   # Overbeek 2025 Online Resource Material 1 $OMEGA 4; Table 2 CL_int BSV 48.1% = sqrt(exp(0.208) - 1) * 100

    # ---- Residual unexplained variability ---------------------------------
    # NONMEM $ERROR: Y = IPRED + IPRED*ERR(1) + ERR(2), i.e. independent
    # proportional and additive terms combining in variance. This is nlmixr2's
    # default combined2 form, so propSd / addSd are the square roots of the
    # $SIGMA variances. Table 2's 17.3% is the same variance presented through
    # the sqrt(exp(sigma^2) - 1) transformation used for the BSV terms.
    propSd <- 0.172   ; label("Proportional residual standard deviation (unitless)")   # Overbeek 2025 Online Resource Material 1 $SIGMA 1 = 0.0296; sqrt(0.0296) = 0.1720; Table 2 reports 17.3% (RSE 3%)
    addSd  <- 0.0332  ; label("Additive residual standard deviation (mg/L)")           # Overbeek 2025 Online Resource Material 1 $SIGMA 2 = 0.0011; sqrt(0.0011) = 0.03317; Table 2 reports 0.033 mg/L (RSE 13.2%)
  })

  model({
    # ---- Individual parameters with allometric scaling to 70 kg ----------
    ktr <- exp(lktr + etalktr) * (WT / 70)^e_wt_ktr
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc

    # Liver volume (Eq 5). A function of the raw body weight in kg, not of
    # WT/70: at 70 kg it gives 0.10 * 70^0.59 = 1.226 L.
    v_liver <- v_liver_coef * WT^e_wt_v_liver

    # Hepatic plasma flow (Eq 1), allometrically scaled as a flow parameter:
    # QHP = QH * (1 - Ht). At 70 kg this is 90 * 0.56 = 50.4 L/h.
    fq <- q_liver * (1 - hct) * (WT / 70)^e_wt_fq

    # Intrinsic clearance (Eq 6) with the PROACTIVE study covariate (Eq 7).
    clint <- exp(lclint + etalclint) * v_liver *
      e_study_proactive_clint^STUDY_PROACTIVE

    # ---- Well-stirred liver model (Eqs 2-3) ------------------------------
    # At 70 kg in the non-PROACTIVE cohorts: CLint = 322 * 1.226 = 394.9 L/h,
    # EH = 394.9 * 0.025 / (50.4 + 394.9 * 0.025) = 0.164 and
    # CLH = 0.164 * 50.4 = 8.26 L/h, matching the apparent clearance of
    # 8.3 L/h quoted in the Discussion (limitations paragraph).
    eh  <- (clint * fu) / (fq + clint * fu)
    clh <- eh * fq

    # ---- Concentrations ---------------------------------------------------
    Cc     <- central / vc
    Cliver <- liver / v_liver

    # ---- ODE system -------------------------------------------------------
    # Erlang absorption: depot -> transit1 -> transit2 -> transit3 -> liver,
    # every transfer sharing the single rate constant ktr (Results 3.1;
    # Online Resource Material 1 K14 = K45 = K56 = K62 = KTR). The absorbed
    # dose enters the liver, so hepatic first pass is produced by the model
    # rather than by a separate bioavailability term.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3

    # Liver: absorption input, hepatic plasma flow in from central (QHP * Cc),
    # flow out to central carrying the non-extracted fraction
    # (QHP * (1 - EH) * Cliver) and hepatic elimination (CLH * Cliver).
    # Online Resource Material 1 writes these as the rate constants
    # K23 = QHP*(1-EH)/VL, K32 = QHP/V and K20 = CLH/VL.
    d/dt(liver)    <-  ktr * transit3 + fq * Cc -
      fq * (1 - eh) * Cliver - clh * Cliver
    d/dt(central)  <-  fq * (1 - eh) * Cliver - fq * Cc

    # ---- Prehepatic bioavailability ---------------------------------------
    f(depot) <- exp(lfdepot)

    # ---- Residual error ---------------------------------------------------
    Cc ~ prop(propSd) + add(addSd)
  })
}
