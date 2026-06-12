Billard_1995_methoxsalen_blood <- function() {
  description <- "Three-compartment intravenous population PK model for 8-methoxypsoralen (8-MOP, methoxsalen) whole blood concentrations in healthy adult volunteers receiving 5/10/15 mg over 60 min (Billard 1995)"
  reference <- "Billard V, Gambus PL, Barr J, Minto CF, Corash L, Tessman JW, Stickney JL, Shafer SL. The pharmacokinetics of 8-methoxypsoralen following i.v. administration in humans. Br J Clin Pharmacol. 1995;40(4):347-360. doi:10.1111/j.1365-2125.1995.tb04557.x"
  vignette <- "Billard_1995_methoxsalen"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight (single-occasion 60-min IV infusion study). Power scaling with allometric exponent 1.0 (fixed; the paper applies weight 'in simple proportion' to volumes and clearances and did not estimate the exponent); reference 70 kg. Paper reports the structural parameters as L/kg and L/kg/min (Table 5); the 70 kg encoding is a presentation choice that leaves the per-kg structural form unchanged. The weight-proportional model decreased the NONMEM -2LL by 10 vs the non-proportional whole-blood model (Results: Compartmental analysis); the improvement was substantially smaller than for the plasma fit (where -2LL dropped by 77).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in Methods (Compartmental analysis paragraph) but not retained; weight was the only covariate that significantly improved the model. Demographics in Table 1: 32.6 +/- 8.4 years (mean +/- s.d.) across 18 subjects."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened (formula BSA = WT^0.425 * HT^0.725 * 0.007184 footnoted in Methods). 'Models in which the volumes and clearances were proportional to body surface area or lean body mass resulted in similar, but not better, log likelihood values to the weight-proportional model' (Discussion: Compartmental analysis). Not retained."
    ),
    LBM = list(
      description = "Lean body mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (sex-specific Boer formulae footnoted in Methods: LBM_men = 1.1*WT - 128*(WT/HT)^2; LBM_women = 1.07*WT - 148*(WT/HT)^2). Not retained for the same reason as BSA."
    ),
    HCT = list(
      description = "Haematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened (Methods: Compartmental analysis paragraph) but not retained. Note: haematocrit was used downstream in the partition-coefficient analysis (Methods: Partition coefficient; Discussion: Partition coefficient analysis) to back-calculate red-cell concentrations from measured whole-blood and plasma, but did not improve the structural PK fit when tested as a covariate on the volumes and clearances."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened (Methods: Compartmental analysis paragraph) but not retained."
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 18L,
    n_studies            = 1L,
    age_range            = "18-40 years (inclusion criterion)",
    age_mean_sd          = "32.6 +/- 8.4 years",
    weight_mean_sd       = "79.8 +/- 13.7 kg",
    height_mean_sd       = "180.6 +/- 9.8 cm",
    bsa_mean_sd          = "2.00 +/- 0.21 m^2",
    lbm_mean_sd          = "61.9 +/- 9.3 kg",
    haematocrit_mean_sd  = "43.5 +/- 3.2 %",
    albumin_mean_sd      = "4.96 +/- 0.37 g/dL",
    sex_female_pct       = NA_real_,
    disease_state        = "Healthy adult volunteers (Stanford University School of Medicine / Palo Alto VA Medical Center). No history of significant medical illness; no chronic tobacco / alcohol / medication / illicit drug use; normal laboratory blood and urine tests including HBsAg and HIV antibody.",
    dose_range           = "Single IV infusion of 5, 10, or 15 mg 8-MOP over 60 min (n = 6 per dose group). Drug diluted to 600 mL total volume in sterile saline; infusion rate 10 mL/min via volumetric pump.",
    regions              = "Single-center United States study (Palo Alto, CA).",
    notes                = "Same cohort as the plasma fit (Billard_1995_methoxsalen_plasma); the paper ran two parallel NONMEM analyses on the same 18 subjects, one on plasma concentrations and one on whole-blood concentrations. 'Whole blood 8-MOP concentration decreased below the limit of detection earlier than plasma' (Results: Moment analysis), so the whole-blood fit covers a shorter post-infusion window than the plasma fit. Both men and women enrolled (the LBM formula in Methods is sex-specific), but the paper does not report the sex breakdown; sex_female_pct left as NA."
  )

  ini({
    # Structural parameters -- 3-compartment mammillary IV model fit to
    # whole-blood concentrations. Table 5 column 'Whole blood Estimate'.
    # Paper reports L/kg and L/kg/min; encoded here for a reference 70 kg
    # subject, with a power-1 weight effect applied in model() so the
    # per-kg structural form is preserved:
    #   V1  = 0.061 L/kg     * 70 kg =  4.270 L
    #   V2  = 1.15  L/kg     * 70 kg = 80.50  L
    #   V3  = 0.21  L/kg     * 70 kg = 14.70  L
    #   CL1 = 0.015 L/kg/min * 70 kg =  1.050 L/min
    #   CL2 = 0.011 L/kg/min * 70 kg =  0.770 L/min
    #   CL3 = 0.015 L/kg/min * 70 kg =  1.050 L/min
    lvc  <- log(4.270);  label("Central volume V1 for the reference 70 kg subject (L)")                              # Billard 1995 Table 5 whole blood V1 = 0.061 L/kg * 70 kg
    lvp  <- log(80.50);  label("Rapidly equilibrating peripheral volume V2 for the reference 70 kg subject (L)")     # Billard 1995 Table 5 whole blood V2 = 1.15 L/kg * 70 kg
    lvp2 <- log(14.70);  label("Slowly equilibrating peripheral volume V3 for the reference 70 kg subject (L)")      # Billard 1995 Table 5 whole blood V3 = 0.21 L/kg * 70 kg
    lcl  <- log(1.050);  label("Systemic clearance CL1 for the reference 70 kg subject (L/min)")                     # Billard 1995 Table 5 whole blood CL1 = 0.015 L/kg/min * 70 kg
    lq   <- log(0.770);  label("Intercompartmental clearance CL2 to V2 for the reference 70 kg subject (L/min)")     # Billard 1995 Table 5 whole blood CL2 = 0.011 L/kg/min * 70 kg
    lq2  <- log(1.050);  label("Intercompartmental clearance CL3 to V3 for the reference 70 kg subject (L/min)")     # Billard 1995 Table 5 whole blood CL3 = 0.015 L/kg/min * 70 kg

    # Weight effect -- 'volumes and clearances were proportional to weight'
    # (Abstract point 4; Results: Compartmental analysis). The exponent is
    # an explicit structural anchor (=1.0), not an estimated quantity, so
    # it is wrapped in fixed(). Six separate fixed-1.0 parameters carry
    # one source-trace comment per affected structural parameter.
    e_wt_cl  <- fixed(1.0); label("Allometric power exponent of body weight on systemic clearance CL1 (unitless)")              # Billard 1995 Results: Compartmental analysis -- 'weight applied to volumes and clearances in simple proportion'
    e_wt_q   <- fixed(1.0); label("Allometric power exponent of body weight on intercompartmental clearance CL2 (unitless)")    # Billard 1995 Results: Compartmental analysis
    e_wt_q2  <- fixed(1.0); label("Allometric power exponent of body weight on intercompartmental clearance CL3 (unitless)")    # Billard 1995 Results: Compartmental analysis
    e_wt_vc  <- fixed(1.0); label("Allometric power exponent of body weight on central volume V1 (unitless)")                   # Billard 1995 Results: Compartmental analysis
    e_wt_vp  <- fixed(1.0); label("Allometric power exponent of body weight on rapidly equilibrating peripheral volume V2 (unitless)") # Billard 1995 Results: Compartmental analysis
    e_wt_vp2 <- fixed(1.0); label("Allometric power exponent of body weight on slowly equilibrating peripheral volume V3 (unitless)")  # Billard 1995 Results: Compartmental analysis

    # Inter-individual variability. Methods (NONMEM parameter model
    # paragraph) defines the IIV as P_i = P_TV * exp(eta_i) with
    # eta ~ N(0, omega^2). The Table 5 'CV' column is omega itself (the
    # log-domain SD; the paper explicitly notes 'omega is thus the standard
    # deviation of the parameter in the log domain'); omega^2 is the
    # variance and is the value nlmixr2 expects on the diagonal. Off-
    # diagonals are not reported -- IIVs are modelled as independent.
    etalvc  ~ 0.0049  # Billard 1995 Table 5 whole blood V1  omega = 0.07; omega^2 = 0.07^2 = 0.0049
    etalvp  ~ 0.4761  # Billard 1995 Table 5 whole blood V2  omega = 0.69; omega^2 = 0.69^2 = 0.4761
    etalvp2 ~ 0.0784  # Billard 1995 Table 5 whole blood V3  omega = 0.28; omega^2 = 0.28^2 = 0.0784
    etalcl  ~ 0.1764  # Billard 1995 Table 5 whole blood CL1 omega = 0.42; omega^2 = 0.42^2 = 0.1764
    etalq   ~ 0.0009  # Billard 1995 Table 5 whole blood CL2 omega = 0.03; omega^2 = 0.03^2 = 0.0009
    etalq2  ~ 0.2601  # Billard 1995 Table 5 whole blood CL3 omega = 0.51; omega^2 = 0.51^2 = 0.2601

    # Residual error -- 'A log normal residual error model was assumed'
    # (Methods: NONMEM residual model paragraph). Table 5 reports the
    # log-domain SD ('Residual CV (sigma)') as 8% for the whole blood fit.
    expSd <- 0.08; label("Log-normal residual error (log-scale SD)")  # Billard 1995 Table 5 whole blood Residual CV (sigma) = 8%
  })

  model({
    # Individual PK parameters with weight-proportional scaling (Results:
    # Compartmental analysis). Reference 70 kg; allometric exponents fixed
    # at 1.0 for every parameter (paper's weight-scaled model).
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl
    q   <- exp(lq   + etalq)   * (WT / 70)^e_wt_q
    q2  <- exp(lq2  + etalq2)  * (WT / 70)^e_wt_q2
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    vp  <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vp
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vp2

    # Micro-rate constants for the 3-compartment mammillary form.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3-compartment mammillary IV model (Methods: Compartmental analysis;
    # Discussion: Compartmental analysis). All drug input is to the
    # central compartment via the IV catheter; no oral / depot absorption.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                       k13 * central - k31 * peripheral2

    # Cc has units mg/L given dose in mg and vc in L; convert to ng/mL by
    # multiplying by 1000. Paper reports whole blood concentrations in
    # ng/mL throughout Tables 4 and Figures 2, 4.
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
  })
}
