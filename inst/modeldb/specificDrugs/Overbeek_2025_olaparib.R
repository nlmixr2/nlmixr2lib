Overbeek_2025_olaparib <- function() {
  description <- "Well-stirred liver model for oral olaparib in patients with solid tumours, with Erlang-type absorption through one transit compartment, a mechanistic hepatic-extraction central/liver disposition driven by unbound intrinsic clearance per litre of liver, a priori allometric scaling to 70 kg, and concomitant cobicistat raising prehepatic bioavailability 1.65-fold while lowering intrinsic clearance to 0.37-fold with its own reduced between-subject variability (Overbeek 2025)"
  reference   <- "Overbeek JK, van Erp NP, Burger DM, den Broeder AA, Koolen SLW, Huitema ADR, ter Heine R. Population Pharmacokinetics of Cobicistat and its Effect on the Pharmacokinetics of the Anticancer Drug Olaparib. Clin Pharmacokinet. 2025;64(3):425-435. doi:10.1007/s40262-025-01480-w"
  vignette    <- "Overbeek_2025_cobicistat_olaparib_boosting"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters in three distinct ways. (1) A priori allometric scaling to a standardised total body weight of 70 kg with exponents 0.75 for flow (hepatic plasma flow), 1 for volume (Vc) and -0.25 for the absorption rate constant ktr (Overbeek 2025 Methods 2.3, final paragraph before Eq 5). (2) The liver volume that scales intrinsic clearance, VL = 0.10 * TBW^0.59 (Eq 5), which is a function of the raw weight in kg and is deliberately NOT re-normalised to 70 kg. (3) Because CLint = theta_CLint * VL (Eq 6), the weight dependence of clearance is carried entirely by VL and no separate allometric exponent is applied to CLint. The olaparib cohort (PROACTIVE) had a median weight of 67 kg, range 54-104 kg (Table 1).",
      source_name        = "WEIGHT"
    ),
    CONMED_COBICISTAT = list(
      description        = "Concomitant cobicistat (PK boosting) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (olaparib 300 mg twice daily monotherapy)",
      notes              = "1 = the boosted arm, olaparib 100 mg twice daily co-administered with cobicistat 150 mg twice daily; 0 = the monotherapy arm, olaparib 300 mg twice daily. Drives two simultaneous multiplicative power-form effects (Overbeek 2025 Results 3.2): a 1.65-fold increase in prehepatic bioavailability, attributed to inhibition of intestinal CYP3A and P-glycoprotein, and a reduction of intrinsic clearance to 0.37-fold, attributed to hepatic CYP3A inhibition. It also selects which of two between-subject variance terms applies to intrinsic clearance, because the paper estimated a separate omega per arm and found variability lower under boosting. Time-varying within a subject: PROACTIVE was a randomised cross-over trial in which every patient contributed one week on each arm, so the indicator switches when the arm changes. The indicator is confounded with the olaparib dose level by design (100 mg boosted vs 300 mg unboosted), but the dose difference is carried by the amt column rather than by this covariate, so the coefficients are dose-independent. Named BOOST in the NONMEM control stream (Online Resource Material 2, $INPUT).",
      source_name        = "BOOST"
    )
  )

  covariatesDataExcluded <- list(
    AUC_COBICISTAT = list(
      description = "Cobicistat area under the plasma concentration-time curve over one dosing interval (AUCtau)",
      units       = "mg*h/L",
      type        = "continuous",
      notes       = "Screened but not retained. Overbeek 2025 Methods 2.5 pre-specified that, if a physiologically plausible relationship were seen, cobicistat AUCtau (computed as Dose / (CLH/F), Eq 8) would be tested as a covariate on olaparib intrinsic clearance. Results 3.2 and Online Resource Fig. 3 report no association between cobicistat exposure and the boosted-to-monotherapy intrinsic-clearance ratio, which the Discussion interprets as saturation of the CYP3A-inhibiting effect at 150 mg twice daily. The column is carried in the analysis dataset as AUCCOBI (Online Resource Material 2, $INPUT) but appears in no $PK expression of the final model, so no coefficient exists to encode."
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    liver    = list(analyte = "olaparib", units = "mg", specimen = "tissue", verified = TRUE),
    central  = list(analyte = "olaparib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    n_observations = 261L,
    age_range      = "55-78 years (median 63)",
    age_median     = "63 years",
    weight_range   = "54-104 kg (median 67)",
    weight_median  = "67 kg",
    sex_female_pct = 58.3,
    disease_state  = "Patients with solid tumours receiving olaparib",
    dose_range     = "Olaparib 300 mg orally twice daily as monotherapy for 1 week, and olaparib 100 mg orally twice daily with cobicistat 150 mg twice daily for 1 week",
    regions        = "The Netherlands",
    notes          = "The PROACTIVE trial (NCT05078671), a randomised cross-over comparison of reduced-dose boosted olaparib against standard olaparib monotherapy; demographics from Overbeek 2025 Table 1. Sex is reported as 5 of 12 male (42%), so 58.3% female. All 261 olaparib samples were above the limit of quantification. Because only 12 patients contributed, the paper's predictive checks used 12.5th / 50th / 87.5th percentiles rather than the usual 2.5th / 97.5th (Methods 2.6)."
  )

  ini({
    # ---- Structural parameters, at the 70 kg reference weight -------------
    # (Overbeek 2025 Table 3; NONMEM $THETA of Online Resource Material 2.)
    lktr   <- log(3.51)  ; label("Erlang absorption transition rate constant (ktr, 1/h)")                              # Overbeek 2025 Table 3: k_tr = 3.51 /h (RSE 15.2%); Online Resource Material 2 $THETA 2
    lvc    <- log(31.6)  ; label("Apparent central volume of distribution (Vc/F, L)")                                  # Overbeek 2025 Table 3: V_c = 31.6 L (RSE 8.7%); Online Resource Material 2 $THETA 3
    lclint <- log(45.6)  ; label("Apparent unbound intrinsic clearance per litre of liver, unboosted (CLint/F, L/h/L liver)")  # Overbeek 2025 Table 3: Cl_int = 45.6 L/h/L_liver (RSE 14.4%); Online Resource Material 2 $THETA 4

    # Prehepatic bioavailability in the monotherapy arm. The paper reports all
    # volume and clearance parameters relative to the unknown absolute
    # bioavailability (Methods 2.3), so F1 is an anchor rather than an estimate;
    # hepatic first pass is carried mechanistically by routing the absorbed dose
    # through the liver compartment.
    lfdepot <- fixed(log(1)) ; label("Prehepatic bioavailability without cobicistat (F1, unitless)")                   # Online Resource Material 2 $THETA 1 = 1 FIX, with $OMEGA 1 = 0 FIX

    # ---- Cobicistat covariate effects (Eq 7, power form) -----------------
    e_conmed_cobicistat_fdepot <- 1.65 ; label("Ratio of prehepatic bioavailability with cobicistat to without (unitless)")     # Overbeek 2025 Table 3: F1 cobicistat = 1.65 (RSE 6%); Results 3.2 "65% increase in prehepatic bioavailability"; Online Resource Material 2 $THETA 5
    e_conmed_cobicistat_clint  <- 0.37 ; label("Ratio of intrinsic clearance with cobicistat to without (unitless)")            # Overbeek 2025 Table 3: CL_int-cobicistat = 0.37 (RSE 6.5%); Results 3.2 "63% decrease in intrinsic clearance"; Online Resource Material 2 $THETA 6

    # ---- Physiological constants of the well-stirred liver model ----------
    # Assumed rather than estimated (Methods 2.3). hct was fixed because
    # haematocrit was not measured in all studies; the paper's sensitivity
    # analysis over 0.30-0.50 changed the estimates negligibly (Results 3.2).
    q_liver      <- fixed(90)    ; label("Hepatic blood flow at the 70 kg reference weight (QH, L/h)")                 # Overbeek 2025 Methods 2.3: "assuming a hepatic blood flood (QH) of 90 L/h"
    hct          <- fixed(0.44)  ; label("Haematocrit, assumed (Ht, unitless)")                                        # Overbeek 2025 Methods 2.3: "hematocrit (Ht) of 0.44"; sensitivity analysis over 0.30-0.50 in Results 3.2
    fu           <- fixed(0.181) ; label("Fraction unbound in plasma, literature value (fu, unitless)")                # Overbeek 2025 Methods 2.3: "an unbound fraction in plasma (fu) of ... 0.181 for olaparib"
    v_liver_coef <- fixed(0.10)  ; label("Coefficient of the liver-volume relation, from Small 2017 (L/kg^0.59)")      # Overbeek 2025 Eq 5: VL = 0.10 * TBW^0.59, citing reference [31] (Small BG et al., Biopharm Drug Dispos 2017;38(4):290-300)
    e_wt_v_liver <- fixed(0.59)  ; label("Body-weight exponent of the liver-volume relation, from Small 2017 (unitless)")  # Overbeek 2025 Eq 5: VL = 0.10 * TBW^0.59, citing reference [31] (Small BG et al., Biopharm Drug Dispos 2017;38(4):290-300)

    # ---- A priori allometric exponents -----------------------------------
    e_wt_fq  <- fixed(0.75)  ; label("Allometric exponent on hepatic plasma flow (unitless)")                          # Overbeek 2025 Methods 2.3: allometric components of 0.75 for flow parameters; Online Resource Material 2 $PK ALLOCL
    e_wt_vc  <- fixed(1)     ; label("Allometric exponent on central volume (unitless)")                               # Overbeek 2025 Methods 2.3: allometric component of 1 for volume parameters; Online Resource Material 2 $PK ALLOV
    e_wt_ktr <- fixed(-0.25) ; label("Allometric exponent on the transition rate constant (unitless)")                 # Overbeek 2025 Methods 2.3: allometric component of -0.25 for absorption parameters; Online Resource Material 2 $PK ALLOKA

    # ---- Inter-individual variability ------------------------------------
    # Log-normal (Eq 4). Table 3 reports these as %CV via sqrt(exp(omega^2) - 1)
    # (Table 3 footnote), so the variances below are the NONMEM $OMEGA values
    # and reproduce the printed %CV exactly. Intrinsic clearance carries two
    # variances, one per treatment arm: the model improved by 51 OFV points
    # when they were separated, and variability was lower under boosting
    # (Results 3.2).
    etalktr               ~ 0.249   # Overbeek 2025 Online Resource Material 2 $OMEGA 2; Table 3 k_tr BSV 53.2% = sqrt(exp(0.249) - 1) * 100
    etalvc                ~ 0.0477  # Overbeek 2025 Online Resource Material 2 $OMEGA 3; Table 3 V_c BSV 22.1% = sqrt(exp(0.0477) - 1) * 100
    etalclint_nocobicistat ~ 0.244  # Overbeek 2025 Online Resource Material 2 $OMEGA 4; Table 3 CL_int-without cobicistat BSV 52.6% = sqrt(exp(0.244) - 1) * 100
    etalclint_cobicistat   ~ 0.171  # Overbeek 2025 Online Resource Material 2 $OMEGA 5; Table 3 CL_int-with cobicistat BSV 43.2% = sqrt(exp(0.171) - 1) * 100

    # ---- Residual unexplained variability ---------------------------------
    # NONMEM $ERROR: Y = F + F*ERR(1) + ERR(2), i.e. independent proportional
    # and additive terms combining in variance. This is nlmixr2's default
    # combined2 form, so propSd / addSd are the square roots of the $SIGMA
    # variances. Table 3's 12.0% is the same variance presented through the
    # sqrt(exp(sigma^2) - 1) transformation used for the BSV terms.
    propSd <- 0.1192 ; label("Proportional residual standard deviation (unitless)")   # Overbeek 2025 Online Resource Material 2 $SIGMA 1 = 0.0142; sqrt(0.0142) = 0.11916; Table 3 reports 12.0% (RSE 33.1%)
    addSd  <- 0.3688 ; label("Additive residual standard deviation (mg/L)")           # Overbeek 2025 Online Resource Material 2 $SIGMA 2 = 0.136; sqrt(0.136) = 0.36878; Table 3 reports 0.369 mg/L (RSE 27.9%)
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

    # Intrinsic clearance (Eq 6) with the cobicistat covariate (Eq 7). The
    # control stream selects between two etas with IF (BOOST.EQ.0) / IF
    # (BOOST.EQ.1); the indicator-weighted sum below is the equivalent form,
    # since exactly one of the two weights is 1 for any given record.
    clint <- exp(lclint) * v_liver *
      e_conmed_cobicistat_clint^CONMED_COBICISTAT *
      exp((1 - CONMED_COBICISTAT) * etalclint_nocobicistat +
            CONMED_COBICISTAT * etalclint_cobicistat)

    # ---- Well-stirred liver model (Eqs 2-3) ------------------------------
    # At 70 kg without cobicistat: CLint = 45.6 * 1.226 = 55.9 L/h,
    # EH = 55.9 * 0.181 / (50.4 + 55.9 * 0.181) = 0.167 and
    # CLH = 0.167 * 50.4 = 8.43 L/h, matching the olaparib clearance of
    # 8.4 L/h quoted in the Discussion. With cobicistat CLint falls to
    # 20.7 L/h, EH to 0.069 and CLH to 3.49 L/h.
    eh  <- (clint * fu) / (fq + clint * fu)
    clh <- eh * fq

    # ---- Concentrations ---------------------------------------------------
    Cc     <- central / vc
    Cliver <- liver / v_liver

    # ---- ODE system -------------------------------------------------------
    # Erlang absorption through a single transit compartment:
    # depot -> transit1 -> liver, both transfers sharing the rate constant ktr
    # (Results 3.2; Online Resource Material 2 K12 = K23 = KTR). The absorbed
    # dose enters the liver, so hepatic first pass is produced by the model
    # rather than by a separate bioavailability term.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ktr * transit1

    # Liver: absorption input, hepatic plasma flow in from central (QHP * Cc),
    # flow out to central carrying the non-extracted fraction
    # (QHP * (1 - EH) * Cliver) and hepatic elimination (CLH * Cliver).
    # Online Resource Material 2 writes these as the rate constants
    # K34 = QHP*(1-EH)/VL, K43 = QHP/V and K30 = CLH/VL.
    d/dt(liver)    <-  ktr * transit1 + fq * Cc -
      fq * (1 - eh) * Cliver - clh * Cliver
    d/dt(central)  <-  fq * (1 - eh) * Cliver - fq * Cc

    # ---- Prehepatic bioavailability, raised by cobicistat -----------------
    f(depot) <- exp(lfdepot) *
      e_conmed_cobicistat_fdepot^CONMED_COBICISTAT

    # ---- Residual error ---------------------------------------------------
    Cc ~ prop(propSd) + add(addSd)
  })
}
