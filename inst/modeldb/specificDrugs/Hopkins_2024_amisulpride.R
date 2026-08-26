Hopkins_2024_amisulpride <- function() {
  description <- paste(
    "Joint PK/PD 'Distribution Model' for nonracemic amisulpride (SEP-4199, a fixed 85:15 ratio of",
    "aramisulpride to esamisulpride) and its brain dopamine D2 receptor occupancy in healthy adults",
    "(Hopkins 2024). Absorption is a three-parallel-depot construct with staggered lag times,",
    "representing regiospecific uptake along the gastrointestinal tract; each depot's first-order rate",
    "constant is scaled by a shared Weibull in vitro dissolution term fitted separately per formulation",
    "prototype, so a single model spans oral solution, immediate-release tablet, five controlled-release",
    "tablet strengths (10/15/20/25/40 percent rate-controlling polymer, one of them fed) and two",
    "multiparticulate (MUPS) capsules. Plasma disposition of total amisulpride enantiomers is two",
    "compartmental with linear clearance. The pharmacodynamic layer is the paper's named contribution:",
    "unbound esamisulpride distributes from plasma into a first brain compartment (plausibly the",
    "blood-brain / blood-CSF barrier), transits to a second brain compartment (plausibly parenchyma),",
    "and D2 receptor occupancy is an Emax-type function of the second compartment with a fixed binding",
    "affinity. This two-step distribution is what simultaneously reproduces the three observations that",
    "defeated conventional PK/PD forms: plasma-brain hysteresis (24 h plasma washout vs 5 day occupancy",
    "washout), a dose-response that does not saturate over the tested range, and an absence of occupancy",
    "accumulation over 7 daily doses. The model supported the model-informed decision to carry a",
    "controlled-release formulation into Phase III as dose-equivalent to the immediate-release form",
    "despite not being bioequivalent in plasma, because the reduced peak plasma concentrations lower QT",
    "prolongation while brain occupancy is preserved."
  )
  reference <- paste(
    "Hopkins SC, Toongsuwan S, Corriveau TJ, Watanabe T, Tsushima Y, Asada T, Lew R, Shi L, Zann V,",
    "Snowden TJ, van der Graaf PH, Darpo B, Searle GE, Rabiner EA, Wilding I, Szabo ST, Galluppi GR,",
    "Koblan KS. Discovery and Model-Informed Drug Development of a Controlled-Release Formulation of",
    "Nonracemic Amisulpride that Reduces Plasma Exposure but Achieves Pharmacodynamic Bioequivalence in",
    "the Brain. Clin Pharmacol Ther. 2024;116(6):1553-1562. doi:10.1002/cpt.3311."
  )
  vignette <- "Hopkins_2024_amisulpride"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Doses enter all three parallel absorption compartments simultaneously; the
  # NONMEM dataset carries three dose records per administration (CMT 1, 2, 3),
  # each scaled by its own bioavailability fraction F1 / F2 / F3 and, for
  # compartments 2 and 3, delayed by ALAG2 / ALAG3.
  dosing <- c("depot1", "depot2", "depot3")

  covariateData <- list(
    OCC = list(
      description        = "Integer occasion index for the inter-occasion variability terms.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Study 5 dosed each subject on several single-dose periods separated by a minimum 7 day",
        "washout, plus a multiple-dose cohort, so occasion is the dosing period. The NONMEM PK",
        "control stream declares OCC in $INPUT and maps it through seven $ABBR REPLACE blocks",
        "(ETA(13,20,27,34,41,48,55) and siblings), i.e. seven occasions carrying IOV on CL, Q,",
        "ALAG2, ALAG3, F1, F2 and F3. Records outside the seven modelled occasions take OCC = 0,",
        "which zeroes every indicator and leaves only the IIV terms."
      ),
      source_name        = "OCC"
    ),
    FORM_SEP4199_CR10 = list(
      description        = "1 = the dose was given as the 10 percent rate-controlling-polymer controlled-release SEP-4199 tablet; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN1 of the NONMEM control stream. Selects F1 / F2 / F3 and the Weibull dissolution pair (TKA1, TGAMA1) = (1.6536, 1.3732).",
      source_name        = "RGN"
    ),
    FORM_SEP4199_CR25 = list(
      description        = "1 = the dose was given as the 25 percent rate-controlling-polymer controlled-release SEP-4199 tablet under fasted conditions; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN2. Fasted counterpart of FORM_SEP4199_CR25_FED. Weibull pair (0.1824, 1.1305).",
      source_name        = "RGN"
    ),
    FORM_SEP4199_IR = list(
      description        = "1 = the dose was given as the immediate-release SEP-4199 tablet; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN3. The clinical reference formulation: the 200 and 400 mg IR tablet carried through the Phase II bipolar depression study (NCT03543410) against which the controlled-release prototypes were compared. Weibull pair (5.7466, 1.4364).",
      source_name        = "RGN"
    ),
    FORM_SEP4199_SOLUTION = list(
      description        = "1 = the dose was given as an oral solution of SEP-4199; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Regimen RGN4, used by the Study 3 single-dose and Study 4 timed-aliquot PET cohorts. The D2",
        "receptor occupancy control stream sets the Weibull dissolution term WB = 1 for this regimen",
        "(IF (RGN.EQ.4) WB = 1) because a solution has no solid-state dissolution step; the PK control",
        "stream instead carries a nominal fast Weibull pair (7.98, 0.905) that reaches WB > 0.95 within",
        "about 30 min. This model follows the D2 occupancy script and switches the term off for solution."
      ),
      source_name        = "RGN"
    ),
    FORM_SEP4199_CR15 = list(
      description        = "1 = the dose was given as the 15 percent rate-controlling-polymer controlled-release SEP-4199 tablet; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN5. Weibull pair (0.3276, 1.0949).",
      source_name        = "RGN"
    ),
    FORM_SEP4199_CR25_FED = list(
      description        = "1 = the dose was given as the 25 percent rate-controlling-polymer controlled-release SEP-4199 tablet with food; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Regimen RGN6, the fed arm of the 25 percent CR tablet. Its Weibull pair (1.1305, 5.7466) was",
        "not fitted to an in vitro dissolution curve: the supplement states that the parameters for the",
        "25 percent CR tablet with food and for the solution were interpolated from the other fits and",
        "the observed in vivo absorption curves. The two numbers are numerically the TKA1 of RGN2 and",
        "the TGAMA1 of RGN3, so this pair should be read as an operational interpolation rather than an",
        "independent dissolution measurement."
      ),
      source_name        = "RGN"
    ),
    FORM_SEP4199_MUPS30 = list(
      description        = "1 = the dose was given as the multiple-unit pellet system (MUPS) SEP-4199 capsule containing 30 percent polymer; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN7. Weibull pair (0.1800, 0.6988). Belongs to the capsule dosage-form class, which carries its own KA1 / KA2 / KA3 set.",
      source_name        = "RGN"
    ),
    FORM_SEP4199_MUPS225 = list(
      description        = "1 = the dose was given as the multiple-unit pellet system (MUPS) SEP-4199 capsule containing 22.5 percent polymer; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN8. Weibull pair (0.6870, 1.0411). Capsule dosage-form class.",
      source_name        = "RGN"
    ),
    FORM_SEP4199_CR20 = list(
      description        = "1 = the dose was given as the 20 percent rate-controlling-polymer controlled-release SEP-4199 tablet; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN9. Weibull pair (0.8457, 0.5496); Table 5 of the supplement rounds TKA1 to 0.80 while the executed control stream carries 0.8457.",
      source_name        = "RGN"
    ),
    FORM_SEP4199_CR40 = list(
      description        = "1 = the dose was given as the 40 percent rate-controlling-polymer controlled-release SEP-4199 tablet; 0 = any of the other nine study formulations.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Regimen RGN10. Weibull pair (0.1411, 0.7764); Table 5 of the supplement rounds TGAMA1 to 0.70 while the executed control stream carries 0.7764.",
      source_name        = "RGN"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, female indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "SEX is carried in the NONMEM $INPUT record of both control streams and is tabulated per study in Table 1 of the paper, but it is not referenced in $PK or $DES of either script and no sex effect is reported."
    ),
    FED = list(
      description = "Fed-versus-fasted indicator for the dosing occasion.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "FED is carried in the NONMEM $INPUT record but is not referenced in $PK or $DES. The single fed",
        "arm of the study (the 25 percent CR tablet with food) is instead carried structurally as its own",
        "regimen, FORM_SEP4199_CR25_FED, with a distinct F1 / F2 / F3 triple and Weibull pair, so the",
        "food effect is absorbed into the formulation indicator rather than modelled as a separate",
        "covariate coefficient."
      )
    )
  )

  compartmentData <- list(
    depot1 = list(analyte = "amisulpride enantiomers (total)", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "amisulpride enantiomers (total)", units = "mg", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "amisulpride enantiomers (total)", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amisulpride enantiomers (total)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amisulpride enantiomers (total)", units = "mg", specimen = "plasma", verified = TRUE),
    # The two brain states of the Distribution Model are concentration-like, not
    # amounts: the influx term is a rate constant multiplying a plasma
    # concentration, so effect1 and effect2 carry ng/mL and are directly
    # comparable with the fixed binding affinity kd = 9 ng/mL. The paper is
    # explicit that the model does not include a mechanistic description of
    # brain physiology, so the specimen is recorded as not applicable.
    effect1 = list(analyte = "esamisulpride (unbound)", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect2 = list(analyte = "esamisulpride (unbound)", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 181L,
    n_studies      = 6L,
    age_range      = "32-39 years (per-study means; Table 1)",
    disease_state  = "healthy volunteers",
    sex_female_pct = 0,
    race_ethnicity = c(White = 68, Black = 19, Other = 14),
    dose_range     = "25-700 mg single oral doses, plus 200 and 400 mg once daily for 7 days",
    regions        = "UK (Studies 1, 3, 4, 5 conducted under MHRA Clinical Trial Authorization)",
    notes          = paste(
      "Table 1 of the paper gives baseline demographics for the four PET / translational-pharmaceutics",
      "studies (Study 1 N = 6, Study 3 N = 11, Study 4 N = 11, Study 5 Part 1 N = 17 and Part 2 N = 18,",
      "imaging arm N = 37). Supplement Table 1 lists all six studies pooled for the PK/PD analysis,",
      "adding Study 2 (N = 33, aramisulpride polysomnography) and an ascending single oral dose study",
      "(N = 48) in Caucasian and Japanese subjects, for the stated total of 181 subjects, 11,626 plasma",
      "samples and 157 D2 receptor occupancy measurements. Sex percentages are given per study rather",
      "than pooled; Studies 1 and 3 were all male, Study 4 was 64 percent male and Study 5 was 53 and 44",
      "percent male in Parts 1 and 2. The race percentages recorded here are those of the N = 37 imaging",
      "arm, the only pooled column in Table 1. sex_female_pct is set to 0 because the only pooled",
      "column reported (the imaging arm) is 78.4 percent male and no pooled female percentage is given;",
      "see the vignette Assumptions section."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Plasma disposition. Values are the $THETA ... FIX entries of the PK
    # control stream (supplement, "PK Fitting NONMEM Script"), which is the
    # parameter set that generated the individual PK predictions consumed by
    # the D2 receptor occupancy script. Supplement Table 2 reports a different
    # set from a separate, explicitly non-converged fitting run; see the
    # vignette Assumptions and deviations section.
    # ---------------------------------------------------------------------
    lcl <- fixed(log(106));  label("Apparent clearance from the central compartment CL/F (L/h)")           # PK $THETA TH1 = 106 FIX
    lq  <- fixed(log(62));   label("Apparent intercompartmental clearance Q/F (L/h)")                      # PK $THETA TH2 = 62 FIX
    lvc <- fixed(log(289));  label("Apparent central volume of distribution Vc/F (L)")                     # PK $THETA TH3 = 289 FIX
    lvp <- fixed(log(981));  label("Apparent peripheral volume of distribution Vp/F (L)")                  # PK $THETA TH4 = 981 FIX

    # ALAG3 is coded in $PK as TAL3 = TAL2 + TAL3T, so THETA(6) is the
    # INCREMENT of the third absorption compartment's lag over the second's,
    # not an absolute lag. The absolute lag is 0.241 + 2.4 = 2.641 h.
    ltlag2 <- fixed(log(0.241)); label("Lag time before absorption compartment 2 opens (h)")               # PK $THETA TH5 = 0.241 FIX
    ltlag3 <- fixed(log(2.4));   label("Additional lag of absorption compartment 3 beyond compartment 2 (h)") # PK $THETA TH6 = 2.4 FIX

    # ---------------------------------------------------------------------
    # Asymptotic (Weibull-scaled) absorption rate constants, one triple per
    # dosage-form class. FORM = 1 solution, 2 tablet, 3 MUPS capsule.
    # ---------------------------------------------------------------------
    lkamax1_sol <- fixed(log(0.657)); label("Absorption rate constant, compartment 1, oral solution (1/h)")  # PK $THETA TH7  = 0.657 FIX
    lkamax1_tab <- fixed(log(0.214)); label("Absorption rate constant, compartment 1, tablet (1/h)")         # PK $THETA TH8  = 0.214 FIX
    lkamax1_cap <- fixed(log(0.197)); label("Absorption rate constant, compartment 1, MUPS capsule (1/h)")   # PK $THETA TH9  = 0.197 FIX
    lkamax2_sol <- fixed(log(0.283)); label("Absorption rate constant, compartment 2, oral solution (1/h)")  # PK $THETA TH10 = 0.283 FIX
    lkamax2_tab <- fixed(log(1.33));  label("Absorption rate constant, compartment 2, tablet (1/h)")         # PK $THETA TH11 = 1.33 FIX
    lkamax2_cap <- fixed(log(0.433)); label("Absorption rate constant, compartment 2, MUPS capsule (1/h)")   # PK $THETA TH12 = 0.433 FIX
    lkamax3_sol <- fixed(log(8.05));  label("Absorption rate constant, compartment 3, oral solution (1/h)")  # PK $THETA TH13 = 8.05 FIX
    lkamax3_tab <- fixed(log(1.36));  label("Absorption rate constant, compartment 3, tablet (1/h)")         # PK $THETA TH14 = 1.36 FIX
    lkamax3_cap <- fixed(log(1.3));   label("Absorption rate constant, compartment 3, MUPS capsule (1/h)")   # PK $THETA TH15 = 1.3 FIX

    # ---------------------------------------------------------------------
    # Per-regimen fractions absorbed through each of the three parallel
    # absorption processes. These are relative-bioavailability multipliers on
    # an apparent-clearance model, so they are not constrained to sum to 1.
    # ---------------------------------------------------------------------
    lf1_cr10     <- fixed(log(1.39));    label("Fraction absorbed via process 1, 10% CR tablet (unitless)")            # PK $THETA TH16 = 1.39 FIX
    lf1_cr25     <- fixed(log(0.759));   label("Fraction absorbed via process 1, 25% CR tablet (unitless)")            # PK $THETA TH17 = 0.759 FIX
    lf1_ir       <- fixed(log(1.39));    label("Fraction absorbed via process 1, IR tablet (unitless)")                # PK $THETA TH18 = 1.39 FIX
    lf1_solution <- fixed(log(1.99));    label("Fraction absorbed via process 1, oral solution (unitless)")            # PK $THETA TH19 = 1.99 FIX
    lf1_cr15     <- fixed(log(0.807));   label("Fraction absorbed via process 1, 15% CR tablet (unitless)")            # PK $THETA TH20 = 0.807 FIX
    lf1_cr25fed  <- fixed(log(1.04));    label("Fraction absorbed via process 1, 25% CR tablet fed (unitless)")        # PK $THETA TH21 = 1.04 FIX
    lf1_mups30   <- fixed(log(0.0092));  label("Fraction absorbed via process 1, MUPS capsule 30% polymer (unitless)") # PK $THETA TH22 = 0.0092 FIX
    lf1_mups225  <- fixed(log(1.07));    label("Fraction absorbed via process 1, MUPS capsule 22.5% polymer (unitless)") # PK $THETA TH23 = 1.07 FIX
    lf1_cr20     <- fixed(log(0.942));   label("Fraction absorbed via process 1, 20% CR tablet (unitless)")            # PK $THETA TH24 = 0.942 FIX
    lf1_cr40     <- fixed(log(0.638));   label("Fraction absorbed via process 1, 40% CR tablet (unitless)")            # PK $THETA TH25 = 0.638 FIX

    lf2_cr10     <- fixed(log(0.0599));  label("Fraction absorbed via process 2, 10% CR tablet (unitless)")            # PK $THETA TH26 = 0.0599 FIX
    lf2_cr25     <- fixed(log(0.36));    label("Fraction absorbed via process 2, 25% CR tablet (unitless)")            # PK $THETA TH27 = 0.36 FIX
    lf2_ir       <- fixed(log(0.224));   label("Fraction absorbed via process 2, IR tablet (unitless)")                # PK $THETA TH28 = 0.224 FIX
    lf2_solution <- fixed(log(0.198));   label("Fraction absorbed via process 2, oral solution (unitless)")            # PK $THETA TH29 = 0.198 FIX
    lf2_cr15     <- fixed(log(0.201));   label("Fraction absorbed via process 2, 15% CR tablet (unitless)")            # PK $THETA TH30 = 0.201 FIX
    lf2_cr25fed  <- fixed(log(2.21e-7)); label("Fraction absorbed via process 2, 25% CR tablet fed (unitless)")        # PK $THETA TH31 = 0.000000221 FIX
    lf2_mups30   <- fixed(log(0.98));    label("Fraction absorbed via process 2, MUPS capsule 30% polymer (unitless)") # PK $THETA TH32 = 0.98 FIX
    lf2_mups225  <- fixed(log(0.299));   label("Fraction absorbed via process 2, MUPS capsule 22.5% polymer (unitless)") # PK $THETA TH33 = 0.299 FIX
    lf2_cr20     <- fixed(log(0.251));   label("Fraction absorbed via process 2, 20% CR tablet (unitless)")            # PK $THETA TH34 = 0.251 FIX
    lf2_cr40     <- fixed(log(0.439));   label("Fraction absorbed via process 2, 40% CR tablet (unitless)")            # PK $THETA TH35 = 0.439 FIX

    lf3_cr10     <- fixed(log(0.542));   label("Fraction absorbed via process 3, 10% CR tablet (unitless)")            # PK $THETA TH36 = 0.542 FIX
    lf3_cr25     <- fixed(log(0.287));   label("Fraction absorbed via process 3, 25% CR tablet (unitless)")            # PK $THETA TH37 = 0.287 FIX
    lf3_ir       <- fixed(log(0.392));   label("Fraction absorbed via process 3, IR tablet (unitless)")                # PK $THETA TH38 = 0.392 FIX
    lf3_solution <- fixed(log(0.00842)); label("Fraction absorbed via process 3, oral solution (unitless)")            # PK $THETA TH39 = 0.00842 FIX
    lf3_cr15     <- fixed(log(0.521));   label("Fraction absorbed via process 3, 15% CR tablet (unitless)")            # PK $THETA TH40 = 0.521 FIX
    lf3_cr25fed  <- fixed(log(0.234));   label("Fraction absorbed via process 3, 25% CR tablet fed (unitless)")        # PK $THETA TH41 = 0.234 FIX
    lf3_mups30   <- fixed(log(1.15));    label("Fraction absorbed via process 3, MUPS capsule 30% polymer (unitless)") # PK $THETA TH42 = 1.15 FIX
    lf3_mups225  <- fixed(log(0.885));   label("Fraction absorbed via process 3, MUPS capsule 22.5% polymer (unitless)") # PK $THETA TH43 = 0.885 FIX
    lf3_cr20     <- fixed(log(0.491));   label("Fraction absorbed via process 3, 20% CR tablet (unitless)")            # PK $THETA TH44 = 0.491 FIX
    lf3_cr40     <- fixed(log(0.393));   label("Fraction absorbed via process 3, 40% CR tablet (unitless)")            # PK $THETA TH45 = 0.393 FIX

    # ---------------------------------------------------------------------
    # Weibull in vitro dissolution pairs, fitted outside NONMEM to the
    # measured dissolution curves (supplement Figure 1 / Table 5) and
    # hardcoded per regimen in $PK. `ra` is the rate-scaling parameter and
    # `gam1` the shape, in ka(t) = kamax * (1 - exp(-(ra * tad)^gam1)).
    # Values are taken from the executed $PK block, which carries four
    # decimal places; Table 5 rounds two of them differently (RGN9 TKA1 to
    # 0.80 and RGN10 TGAMA1 to 0.70).
    # ---------------------------------------------------------------------
    lra_cr10     <- fixed(log(1.6536)); label("Weibull dissolution rate scaling, 10% CR tablet (1/h)")      # PK $PK IF (RGN.EQ.1) TKA1 = 1.6536
    lra_cr25     <- fixed(log(0.1824)); label("Weibull dissolution rate scaling, 25% CR tablet (1/h)")      # PK $PK IF (RGN.EQ.2) TKA1 = 0.1824
    lra_ir       <- fixed(log(5.7466)); label("Weibull dissolution rate scaling, IR tablet (1/h)")          # PK $PK IF (RGN.EQ.3) TKA1 = 5.7466
    lra_solution <- fixed(log(7.98));   label("Weibull dissolution rate scaling, oral solution (1/h)")      # PK $PK IF (RGN.EQ.4) TKA1 = 7.98 (unused; solution bypasses the term)
    lra_cr15     <- fixed(log(0.3276)); label("Weibull dissolution rate scaling, 15% CR tablet (1/h)")      # PK $PK IF (RGN.EQ.5) TKA1 = 0.3276
    lra_cr25fed  <- fixed(log(1.1305)); label("Weibull dissolution rate scaling, 25% CR tablet fed (1/h)")  # PK $PK IF (RGN.EQ.6) TKA1 = 1.1305 (interpolated)
    lra_mups30   <- fixed(log(0.18));   label("Weibull dissolution rate scaling, MUPS capsule 30% polymer (1/h)")   # PK $PK IF (RGN.EQ.7) TKA1 = 0.1800
    lra_mups225  <- fixed(log(0.687));  label("Weibull dissolution rate scaling, MUPS capsule 22.5% polymer (1/h)") # PK $PK IF (RGN.EQ.8) TKA1 = 0.6870
    lra_cr20     <- fixed(log(0.8457)); label("Weibull dissolution rate scaling, 20% CR tablet (1/h)")      # PK $PK IF (RGN.EQ.9) TKA1 = 0.8457
    lra_cr40     <- fixed(log(0.1411)); label("Weibull dissolution rate scaling, 40% CR tablet (1/h)")      # PK $PK IF (RGN.EQ.10) TKA1 = 0.1411

    lgam1_cr10     <- fixed(log(1.3732)); label("Weibull dissolution shape, 10% CR tablet (unitless)")      # PK $PK IF (RGN.EQ.1) TGAMA1 = 1.3732
    lgam1_cr25     <- fixed(log(1.1305)); label("Weibull dissolution shape, 25% CR tablet (unitless)")      # PK $PK IF (RGN.EQ.2) TGAMA1 = 1.1305
    lgam1_ir       <- fixed(log(1.4364)); label("Weibull dissolution shape, IR tablet (unitless)")          # PK $PK IF (RGN.EQ.3) TGAMA1 = 1.4364
    lgam1_solution <- fixed(log(0.905));  label("Weibull dissolution shape, oral solution (unitless)")      # PK $PK IF (RGN.EQ.4) TGAMA1 = 0.905 (unused; solution bypasses the term)
    lgam1_cr15     <- fixed(log(1.0949)); label("Weibull dissolution shape, 15% CR tablet (unitless)")      # PK $PK IF (RGN.EQ.5) TGAMA1 = 1.0949
    lgam1_cr25fed  <- fixed(log(5.7466)); label("Weibull dissolution shape, 25% CR tablet fed (unitless)")  # PK $PK IF (RGN.EQ.6) TGAMA1 = 5.7466 (interpolated)
    lgam1_mups30   <- fixed(log(0.6988)); label("Weibull dissolution shape, MUPS capsule 30% polymer (unitless)")   # PK $PK IF (RGN.EQ.7) TGAMA1 = 0.6988
    lgam1_mups225  <- fixed(log(1.0411)); label("Weibull dissolution shape, MUPS capsule 22.5% polymer (unitless)") # PK $PK IF (RGN.EQ.8) TGAMA1 = 1.0411
    lgam1_cr20     <- fixed(log(0.5496)); label("Weibull dissolution shape, 20% CR tablet (unitless)")      # PK $PK IF (RGN.EQ.9) TGAMA1 = 0.5496
    lgam1_cr40     <- fixed(log(0.7764)); label("Weibull dissolution shape, 40% CR tablet (unitless)")      # PK $PK IF (RGN.EQ.10) TGAMA1 = 0.7764

    # ---------------------------------------------------------------------
    # Distribution Model (brain) layer. Final estimates from supplement
    # Table 3, "NONMEM D2 Receptor Occupancy Fitting Results"; the $THETA
    # entries of the occupancy control stream are the initial estimates for
    # that run, not the final values.
    # ---------------------------------------------------------------------
    lkinf_effect1 <- log(0.0221); label("Influx rate constant from plasma into brain compartment 1, KE1 (1/h)")        # Supplement Table 3: KE1 = 0.0221
    lkeff_effect1 <- log(0.0685); label("Efflux rate constant out of brain compartment 1 back to plasma, KE2 (1/h)")   # Supplement Table 3: KE2 = 0.0685
    lkinf_effect2 <- log(0.42);   label("Distribution rate constant from brain compartment 1 to compartment 2, KT (1/h)") # Supplement Table 3: KT = 0.42
    lkeff_effect2 <- log(0.0602); label("Efflux rate constant out of brain compartment 2, KE3 (1/h)")                  # Supplement Table 3: KE3 = 0.0602

    kd     <- fixed(9);    label("D2 receptor binding affinity in brain compartment 2 units, KD (ng/mL)")  # Supplement Table 3: KD = 9 (FIX); occupancy control stream $THETA TH5 9 FIX
    fu     <- fixed(0.83); label("Fraction of amisulpride unbound in plasma, FU (unitless)")               # Supplement Table 3: FU = 0.83 (FIX); occupancy control stream $THETA TH6 0.83 FIX
    fenant <- fixed(0.15); label("Fraction of total measured amisulpride that is the D2-active enantiomer esamisulpride, SP (unitless)") # Supplement Table 3: SP = 0.15 (FIX); the fixed 85:15 aramisulpride:esamisulpride ratio of SEP-4199 (Methods, Drugs)

    # ---------------------------------------------------------------------
    # Inter-individual variability. PK values are the $OMEGA ... FIX diagonal
    # of the PK control stream. The stream fixes the absorption-rate and
    # absorbed-fraction IIVs to exactly zero, so no eta is carried for
    # kamax1-3 or f1-3; a zero-variance eta would make OMEGA singular.
    # ---------------------------------------------------------------------
    etalcl    ~ 0.105  # PK $OMEGA ETA1 - CL IIV  = 0.105 FIX
    etalq     ~ 0.719  # PK $OMEGA ETA2 - Q IIV   = 0.719 FIX
    etalvc    ~ 0.431  # PK $OMEGA ETA3 - Vcen IIV  = 0.431 FIX
    etalvp    ~ 0.389  # PK $OMEGA ETA4 - Vperi IIV = 0.389 FIX
    etaltlag2 ~ 0.01   # PK $OMEGA ETA5 - ALAG2 IIV = 0.01 FIX
    etaltlag3 ~ 0.01   # PK $OMEGA ETA6 - ALAG3 IIV = 0.01 FIX

    # D2 receptor occupancy IIV, supplement Table 3 "IIV / Estimate" column.
    etalkinf_effect1 ~ 0.0314  # Supplement Table 3: KE1 - IIV = 0.0314
    etalkeff_effect1 ~ 1e-06   # Supplement Table 3: KE2 - IIV = 1.00E-06 (effectively no IIV, but not fixed to zero)
    etalkinf_effect2 ~ 0.706   # Supplement Table 3: KT  - IIV = 0.706
    etalkeff_effect2 ~ 0.0341  # Supplement Table 3: KE3 - IIV = 0.0341

    # ---------------------------------------------------------------------
    # Inter-occasion variability, PK layer only. NONMEM encodes this as
    # $OMEGA BLOCK(7) followed by six BLOCK(7) SAME, i.e. seven occasions
    # sharing one 7 x 7 variance block whose off-diagonals are all zero.
    # Occasion 1 carries the estimated variance and occasions 2-7 repeat it
    # via fixed().
    # ---------------------------------------------------------------------
    etaiov_lcl_1 ~ 0.048         # PK $OMEGA BLOCK(7) element (1,1) = 0.048, CL IOV
    etaiov_lcl_2 ~ fixed(0.048)  # BLOCK(7) SAME
    etaiov_lcl_3 ~ fixed(0.048)  # BLOCK(7) SAME
    etaiov_lcl_4 ~ fixed(0.048)  # BLOCK(7) SAME
    etaiov_lcl_5 ~ fixed(0.048)  # BLOCK(7) SAME
    etaiov_lcl_6 ~ fixed(0.048)  # BLOCK(7) SAME
    etaiov_lcl_7 ~ fixed(0.048)  # BLOCK(7) SAME

    etaiov_lq_1 ~ 0.293          # PK $OMEGA BLOCK(7) element (2,2) = 0.293, Q IOV
    etaiov_lq_2 ~ fixed(0.293)   # BLOCK(7) SAME
    etaiov_lq_3 ~ fixed(0.293)   # BLOCK(7) SAME
    etaiov_lq_4 ~ fixed(0.293)   # BLOCK(7) SAME
    etaiov_lq_5 ~ fixed(0.293)   # BLOCK(7) SAME
    etaiov_lq_6 ~ fixed(0.293)   # BLOCK(7) SAME
    etaiov_lq_7 ~ fixed(0.293)   # BLOCK(7) SAME

    etaiov_ltlag2_1 ~ 0.0103         # PK $OMEGA BLOCK(7) element (3,3) = 0.0103, ALAG2 IOV
    etaiov_ltlag2_2 ~ fixed(0.0103)  # BLOCK(7) SAME
    etaiov_ltlag2_3 ~ fixed(0.0103)  # BLOCK(7) SAME
    etaiov_ltlag2_4 ~ fixed(0.0103)  # BLOCK(7) SAME
    etaiov_ltlag2_5 ~ fixed(0.0103)  # BLOCK(7) SAME
    etaiov_ltlag2_6 ~ fixed(0.0103)  # BLOCK(7) SAME
    etaiov_ltlag2_7 ~ fixed(0.0103)  # BLOCK(7) SAME

    etaiov_ltlag3_1 ~ 0.00964         # PK $OMEGA BLOCK(7) element (4,4) = 0.00964, ALAG3 IOV
    etaiov_ltlag3_2 ~ fixed(0.00964)  # BLOCK(7) SAME
    etaiov_ltlag3_3 ~ fixed(0.00964)  # BLOCK(7) SAME
    etaiov_ltlag3_4 ~ fixed(0.00964)  # BLOCK(7) SAME
    etaiov_ltlag3_5 ~ fixed(0.00964)  # BLOCK(7) SAME
    etaiov_ltlag3_6 ~ fixed(0.00964)  # BLOCK(7) SAME
    etaiov_ltlag3_7 ~ fixed(0.00964)  # BLOCK(7) SAME

    etaiov_lf1_1 ~ 0.165         # PK $OMEGA BLOCK(7) element (5,5) = 0.165, F1 IOV
    etaiov_lf1_2 ~ fixed(0.165)  # BLOCK(7) SAME
    etaiov_lf1_3 ~ fixed(0.165)  # BLOCK(7) SAME
    etaiov_lf1_4 ~ fixed(0.165)  # BLOCK(7) SAME
    etaiov_lf1_5 ~ fixed(0.165)  # BLOCK(7) SAME
    etaiov_lf1_6 ~ fixed(0.165)  # BLOCK(7) SAME
    etaiov_lf1_7 ~ fixed(0.165)  # BLOCK(7) SAME

    etaiov_lf2_1 ~ 0.05         # PK $OMEGA BLOCK(7) element (6,6) = 0.05, F2 IOV
    etaiov_lf2_2 ~ fixed(0.05)  # BLOCK(7) SAME
    etaiov_lf2_3 ~ fixed(0.05)  # BLOCK(7) SAME
    etaiov_lf2_4 ~ fixed(0.05)  # BLOCK(7) SAME
    etaiov_lf2_5 ~ fixed(0.05)  # BLOCK(7) SAME
    etaiov_lf2_6 ~ fixed(0.05)  # BLOCK(7) SAME
    etaiov_lf2_7 ~ fixed(0.05)  # BLOCK(7) SAME

    etaiov_lf3_1 ~ 0.2         # PK $OMEGA BLOCK(7) element (7,7) = 0.2, F3 IOV
    etaiov_lf3_2 ~ fixed(0.2)  # BLOCK(7) SAME
    etaiov_lf3_3 ~ fixed(0.2)  # BLOCK(7) SAME
    etaiov_lf3_4 ~ fixed(0.2)  # BLOCK(7) SAME
    etaiov_lf3_5 ~ fixed(0.2)  # BLOCK(7) SAME
    etaiov_lf3_6 ~ fixed(0.2)  # BLOCK(7) SAME
    etaiov_lf3_7 ~ fixed(0.2)  # BLOCK(7) SAME

    # ---------------------------------------------------------------------
    # Residual error. NONMEM reports EPS variances; nlmixr2 wants standard
    # deviations, so each is the square root of the reported $SIGMA entry.
    # ---------------------------------------------------------------------
    propSd     <- sqrt(0.15); label("Proportional residual error on plasma concentration (fraction)")  # PK / occupancy $SIGMA EPS(1) variance 0.15 -> SD 0.3873
    addSd      <- sqrt(10);   label("Additive residual error on plasma concentration (ng/mL)")         # PK / occupancy $SIGMA EPS(2) variance 10 -> SD 3.1623
    addSd_rod2 <- sqrt(52.2); label("Additive residual error on D2 receptor occupancy (percent)")      # Occupancy $SIGMA EPS(3) variance 52.2 -> SD 7.2250
  })

  model({
    # -------------------------------------------------------------------
    # 1. Formulation-derived terms.
    #
    # The ten regimen indicators are mutually exclusive and exhaustive, so
    # each sum below selects exactly one value.
    # -------------------------------------------------------------------
    is_solution <- FORM_SEP4199_SOLUTION
    is_capsule  <- FORM_SEP4199_MUPS30 + FORM_SEP4199_MUPS225
    is_tablet   <- 1 - is_solution - is_capsule

    ra <-
      exp(lra_cr10)     * FORM_SEP4199_CR10 +
      exp(lra_cr25)     * FORM_SEP4199_CR25 +
      exp(lra_ir)       * FORM_SEP4199_IR +
      exp(lra_solution) * FORM_SEP4199_SOLUTION +
      exp(lra_cr15)     * FORM_SEP4199_CR15 +
      exp(lra_cr25fed)  * FORM_SEP4199_CR25_FED +
      exp(lra_mups30)   * FORM_SEP4199_MUPS30 +
      exp(lra_mups225)  * FORM_SEP4199_MUPS225 +
      exp(lra_cr20)     * FORM_SEP4199_CR20 +
      exp(lra_cr40)     * FORM_SEP4199_CR40

    gam1 <-
      exp(lgam1_cr10)     * FORM_SEP4199_CR10 +
      exp(lgam1_cr25)     * FORM_SEP4199_CR25 +
      exp(lgam1_ir)       * FORM_SEP4199_IR +
      exp(lgam1_solution) * FORM_SEP4199_SOLUTION +
      exp(lgam1_cr15)     * FORM_SEP4199_CR15 +
      exp(lgam1_cr25fed)  * FORM_SEP4199_CR25_FED +
      exp(lgam1_mups30)   * FORM_SEP4199_MUPS30 +
      exp(lgam1_mups225)  * FORM_SEP4199_MUPS225 +
      exp(lgam1_cr20)     * FORM_SEP4199_CR20 +
      exp(lgam1_cr40)     * FORM_SEP4199_CR40

    # -------------------------------------------------------------------
    # 2. Occasion indicators and the inter-occasion variability terms.
    #    OCC = 0 (or any value outside 1..7) zeroes every indicator and
    #    leaves only the IIV contribution, matching the NONMEM $ABBR
    #    REPLACE mapping.
    # -------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)

    iov_lcl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 + oc3 * etaiov_lcl_3 +
      oc4 * etaiov_lcl_4 + oc5 * etaiov_lcl_5 + oc6 * etaiov_lcl_6 + oc7 * etaiov_lcl_7
    iov_lq <- oc1 * etaiov_lq_1 + oc2 * etaiov_lq_2 + oc3 * etaiov_lq_3 +
      oc4 * etaiov_lq_4 + oc5 * etaiov_lq_5 + oc6 * etaiov_lq_6 + oc7 * etaiov_lq_7
    iov_ltlag2 <- oc1 * etaiov_ltlag2_1 + oc2 * etaiov_ltlag2_2 + oc3 * etaiov_ltlag2_3 +
      oc4 * etaiov_ltlag2_4 + oc5 * etaiov_ltlag2_5 + oc6 * etaiov_ltlag2_6 + oc7 * etaiov_ltlag2_7
    iov_ltlag3 <- oc1 * etaiov_ltlag3_1 + oc2 * etaiov_ltlag3_2 + oc3 * etaiov_ltlag3_3 +
      oc4 * etaiov_ltlag3_4 + oc5 * etaiov_ltlag3_5 + oc6 * etaiov_ltlag3_6 + oc7 * etaiov_ltlag3_7
    iov_lf1 <- oc1 * etaiov_lf1_1 + oc2 * etaiov_lf1_2 + oc3 * etaiov_lf1_3 +
      oc4 * etaiov_lf1_4 + oc5 * etaiov_lf1_5 + oc6 * etaiov_lf1_6 + oc7 * etaiov_lf1_7
    iov_lf2 <- oc1 * etaiov_lf2_1 + oc2 * etaiov_lf2_2 + oc3 * etaiov_lf2_3 +
      oc4 * etaiov_lf2_4 + oc5 * etaiov_lf2_5 + oc6 * etaiov_lf2_6 + oc7 * etaiov_lf2_7
    iov_lf3 <- oc1 * etaiov_lf3_1 + oc2 * etaiov_lf3_2 + oc3 * etaiov_lf3_3 +
      oc4 * etaiov_lf3_4 + oc5 * etaiov_lf3_5 + oc6 * etaiov_lf3_6 + oc7 * etaiov_lf3_7

    # -------------------------------------------------------------------
    # 3. Individual parameters.
    # -------------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_lcl)
    q  <- exp(lq + etalq + iov_lq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)

    # $PK: TAL2 = THETA(5)*EXP(...); TAL3 = TAL2 + THETA(6)*EXP(...)
    tlag2 <- exp(ltlag2 + etaltlag2 + iov_ltlag2)
    tlag3 <- tlag2 + exp(ltlag3 + etaltlag3 + iov_ltlag3)

    kamax1 <- exp(lkamax1_sol) * is_solution + exp(lkamax1_tab) * is_tablet + exp(lkamax1_cap) * is_capsule
    kamax2 <- exp(lkamax2_sol) * is_solution + exp(lkamax2_tab) * is_tablet + exp(lkamax2_cap) * is_capsule
    kamax3 <- exp(lkamax3_sol) * is_solution + exp(lkamax3_tab) * is_tablet + exp(lkamax3_cap) * is_capsule

    f1 <- exp(iov_lf1) * (
      exp(lf1_cr10)     * FORM_SEP4199_CR10 +
        exp(lf1_cr25)     * FORM_SEP4199_CR25 +
        exp(lf1_ir)       * FORM_SEP4199_IR +
        exp(lf1_solution) * FORM_SEP4199_SOLUTION +
        exp(lf1_cr15)     * FORM_SEP4199_CR15 +
        exp(lf1_cr25fed)  * FORM_SEP4199_CR25_FED +
        exp(lf1_mups30)   * FORM_SEP4199_MUPS30 +
        exp(lf1_mups225)  * FORM_SEP4199_MUPS225 +
        exp(lf1_cr20)     * FORM_SEP4199_CR20 +
        exp(lf1_cr40)     * FORM_SEP4199_CR40)

    f2 <- exp(iov_lf2) * (
      exp(lf2_cr10)     * FORM_SEP4199_CR10 +
        exp(lf2_cr25)     * FORM_SEP4199_CR25 +
        exp(lf2_ir)       * FORM_SEP4199_IR +
        exp(lf2_solution) * FORM_SEP4199_SOLUTION +
        exp(lf2_cr15)     * FORM_SEP4199_CR15 +
        exp(lf2_cr25fed)  * FORM_SEP4199_CR25_FED +
        exp(lf2_mups30)   * FORM_SEP4199_MUPS30 +
        exp(lf2_mups225)  * FORM_SEP4199_MUPS225 +
        exp(lf2_cr20)     * FORM_SEP4199_CR20 +
        exp(lf2_cr40)     * FORM_SEP4199_CR40)

    f3 <- exp(iov_lf3) * (
      exp(lf3_cr10)     * FORM_SEP4199_CR10 +
        exp(lf3_cr25)     * FORM_SEP4199_CR25 +
        exp(lf3_ir)       * FORM_SEP4199_IR +
        exp(lf3_solution) * FORM_SEP4199_SOLUTION +
        exp(lf3_cr15)     * FORM_SEP4199_CR15 +
        exp(lf3_cr25fed)  * FORM_SEP4199_CR25_FED +
        exp(lf3_mups30)   * FORM_SEP4199_MUPS30 +
        exp(lf3_mups225)  * FORM_SEP4199_MUPS225 +
        exp(lf3_cr20)     * FORM_SEP4199_CR20 +
        exp(lf3_cr40)     * FORM_SEP4199_CR40)

    kinf_effect1 <- exp(lkinf_effect1 + etalkinf_effect1)
    keff_effect1 <- exp(lkeff_effect1 + etalkeff_effect1)
    kinf_effect2 <- exp(lkinf_effect2 + etalkinf_effect2)
    keff_effect2 <- exp(lkeff_effect2 + etalkeff_effect2)

    # -------------------------------------------------------------------
    # 4. Micro-constants.
    # -------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -------------------------------------------------------------------
    # 5. Weibull dissolution scaling of the three absorption rate
    #    constants, WB = 1 - exp(-(KAW * TAD)^GAMAW) in $PK. TAD is time
    #    after the most recent dose; depot1 takes the dose with no lag, so
    #    tad(depot1) is the NONMEM TAD column. The occupancy control stream
    #    sets WB = 1 for the oral solution, which has no dissolution step.
    # -------------------------------------------------------------------
    wb <- is_solution + (1 - is_solution) * (1 - exp(-(ra * tad(depot1))^gam1))

    # $ERROR computes CT = A(4)/V4 directly, which means the NONMEM dataset
    # carried dose amounts in ug: the assay curve range and every reported
    # concentration are ng/mL (equivalently ug/L), and 289 L is an apparent
    # central volume in litres. Doses here are declared in mg to match the
    # published dose levels, so the quotient central/vc is in mg/L and needs
    # the 1000-fold conversion to ng/mL. The factor is a pure unit
    # conversion, not a fitted scaling parameter. It is applied before the
    # effect compartments because their influx is driven by the plasma
    # concentration on the ng/mL scale, the same scale as kd = 9 ng/mL.
    Cc <- 1000 * central / vc

    # -------------------------------------------------------------------
    # 6. ODE system. Compartments 1-5 are $DES of the PK script; 6 and 7
    #    are the Distribution Model states added by the occupancy script.
    # -------------------------------------------------------------------
    d/dt(depot1) <- -kamax1 * wb * depot1
    d/dt(depot2) <- -kamax2 * wb * depot2
    d/dt(depot3) <- -kamax3 * wb * depot3
    d/dt(central) <- kamax1 * wb * depot1 + kamax2 * wb * depot2 + kamax3 * wb * depot3 +
      k21 * peripheral1 - (kel + k12) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # DADT(6) = KE1*FU*SP*(A(4)/V4) - KE2*A(6) - KT*A(6)
    d/dt(effect1) <- kinf_effect1 * fu * fenant * Cc - keff_effect1 * effect1 - kinf_effect2 * effect1
    # DADT(7) = KT*A(6) - KE3*A(7)
    d/dt(effect2) <- kinf_effect2 * effect1 - keff_effect2 * effect2

    # -------------------------------------------------------------------
    # 7. Bioavailability and lag times.
    # -------------------------------------------------------------------
    f(depot1) <- f1
    f(depot2) <- f2
    f(depot3) <- f3
    alag(depot2) <- tlag2
    alag(depot3) <- tlag3

    # -------------------------------------------------------------------
    # 8. Observations. $ERROR: CT = A(4)/V4 with a combined error model,
    #    and D2P = (A(7)/(KD + A(7)))*100 with an additive error model.
    # -------------------------------------------------------------------
    rod2 <- 100 * effect2 / (kd + effect2)

    Cc ~ prop(propSd) + add(addSd)
    rod2 ~ add(addSd_rod2)
  })
}
