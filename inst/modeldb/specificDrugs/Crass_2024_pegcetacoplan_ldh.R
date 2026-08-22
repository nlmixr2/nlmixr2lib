Crass_2024_pegcetacoplan_ldh <- function() {
  description <- "Direct sigmoidal Emax PK/PD model relating pegcetacoplan serum concentration to serum lactate dehydrogenase in adults with paroxysmal nocturnal hemoglobinuria (Crass 2024). Lactate dehydrogenase is a direct (no-delay) function of the concurrent pegcetacoplan concentration, ldh = rbase * (1 - emax * Cc^hill / (ec50^hill + Cc^hill)), with the pegcetacoplan concentration supplied by the companion one-compartment-plus-single-transit population PK model of the same paper, whose parameters are carried here as fixed values because the PK/PD model was fitted sequentially on individual empirical Bayes estimates of the PK parameters. Baseline lactate dehydrogenase and the maximal fractional suppression are each estimated separately for patients who were eculizumab-naive versus eculizumab-treated at baseline, and concurrent eculizumab co-treatment adds a further shift to the maximal suppression in the baseline-eculizumab stratum. The baseline is estimated in the log domain and the maximal effect in the logit domain. Inter-individual variability is log-normal on the baseline, the logit maximal effect, and EC50 as a correlated 3x3 block; residual error is additive on the log scale (log-normal)."
  reference <- "Crass RL, Smith B, Adriaens S, Chapel S, Langdon G. Population Pharmacokinetic and Pharmacokinetic/Pharmacodynamic Analyses of Pegcetacoplan in Patients with Paroxysmal Nocturnal Hemoglobinuria. Drugs R D. 2024;24(4):565-577. doi:10.1007/s40268-024-00500-7"
  vignette <- "Crass_2024_pegcetacoplan"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(
      analyte = "pegcetacoplan", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "pegcetacoplan", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "pegcetacoplan", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline body weight. Enters only the carried-forward PK layer (power exponents on CL and Vc referenced to 70 kg). Body weight was screened as a covariate on the lactate-dehydrogenase PK/PD parameters and was not retained (Crass 2024 Sect. 3.4.2: 'All other factors (i.e., sex, baseline CrCl, race, age, body weight, and baseline C3 level) did not have meaningful effects on the LDH response').",
      source_name        = "BWT"
    ),
    DIS_PNH = list(
      description        = "Paroxysmal nocturnal hemoglobinuria indicator: 1 = patient with PNH, 0 = non-PNH participant.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-PNH participant)",
      notes              = "Time-fixed per subject. Enters only the carried-forward PK layer, where it applies the fractional clearance increase e_dis_pnh_cl = 0.257. Every subject in the lactate-dehydrogenase PK/PD analysis set is a patient with PNH (DIS_PNH = 1), so the term is inert in this cohort; it is retained so the PK layer is the same structure as the parent model `Crass_2024_pegcetacoplan`.",
      source_name        = "PNH"
    ),
    FORM_PEGCET_LYOPHILIZED = list(
      description        = "Pegcetacoplan lyophilized-powder formulation indicator: 1 = lyophilized powder reconstituted before subcutaneous administration, 0 = a ready-to-use solution formulation.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ready-to-use solution formulations)",
      notes              = "Per-regimen categorical indicator. Enters only the carried-forward PK layer, scaling subcutaneous bioavailability by (1 + 0.220 x FORM_PEGCET_LYOPHILIZED). Retained so the PK layer matches the parent model `Crass_2024_pegcetacoplan`.",
      source_name        = "FORM (level 4 = POWDER)"
    ),
    CONMED_ECULIZUMAB_BL = list(
      description        = "Eculizumab (complement C5 inhibitor) treatment status at baseline: 1 = receiving eculizumab at study entry, 0 = eculizumab-naive at study entry.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (eculizumab-naive at baseline; the PADDOCK, PALOMINO, and PRINCE cohorts)",
      notes              = "Time-fixed per subject. Selects BOTH the baseline lactate-dehydrogenase stratum and the maximal-suppression stratum: `IF(BECU.EQ.1) TVBLDH = THETA(2)` and `IF (BECU.EQ.1) TVEMAX = THETA(4) + THETA(8)*ECU` in the ESM Table 1 LDH control stream. Eculizumab-naive patients have a typical baseline of exp(7.56) = 1920 U/L and a maximal 91.7% suppression; patients already on eculizumab at baseline have a typical baseline of exp(5.52) = 249 U/L and a maximal 20.0% suppression, because terminal-complement blockade has already suppressed intravascular hemolysis. The PK/PD analysis set is 54% eculizumab-treated / 46% eculizumab-naive at baseline (Crass 2024 Sect. 3.1.2).",
      source_name        = "BECU"
    ),
    CONMED_ECULIZUMAB = list(
      description        = "Concurrent eculizumab (complement C5 inhibitor) co-administration indicator: 1 = receiving eculizumab together with pegcetacoplan at the current time, 0 = pegcetacoplan monotherapy at the current time.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pegcetacoplan monotherapy at the current time)",
      notes              = "TIME-VARYING within a subject. Distinct from CONMED_ECULIZUMAB_BL, which is the fixed baseline status. It is non-zero during the PEGASUS dual-therapy periods (the 4-week run-in on eculizumab + pegcetacoplan, and the 4-week transition for patients randomised to eculizumab who later crossed to pegcetacoplan; Crass 2024 Table 1 footnote a). It enters the logit maximal suppression additively but only within the baseline-eculizumab stratum: `IF (BECU.EQ.1) TVEMAX = THETA(4) + THETA(8)*ECU` with THETA(8) = 0.783, which lifts the maximal suppression from inverse-logit(-1.39) = 20.0% during monotherapy to inverse-logit(-1.39 + 0.783) = 35.3% during co-treatment (Crass 2024 Sect. 3.4.2 reports 35.0%). For an eculizumab-naive patient the term is not applied at all, matching the control stream's guard.",
      source_name        = "ECU"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the full covariate model on the lactate-dehydrogenase PK/PD parameters and eliminated by the backward-elimination step (alpha = 0.001); no effect is reported in ESM Table 5 or Crass 2024 Table 2. Analysis-set median 45 years (range 19-81 years).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened and not retained on the lactate-dehydrogenase response (Crass 2024 Sect. 3.4.2). Retained on the hemoglobin response in the sibling model `Crass_2024_pegcetacoplan_hemoglobin`, so the contrast between the two endpoints is a genuine model-selection outcome and not an omission here.",
      source_name        = "SEXF"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance estimated with the Cockcroft-Gault equation, in mL/min (NOT body-surface-area normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained on the lactate-dehydrogenase response (Crass 2024 Sect. 3.4.2). Retained on the hemoglobin maximal response in the sibling model `Crass_2024_pegcetacoplan_hemoglobin`.",
      source_name        = "BCRCL"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Race was screened and not retained on the lactate-dehydrogenase response (Crass 2024 Sect. 3.4.2). The PK/PD analysis set is 39% White and 38% Asian (Sect. 3.1.2). Recorded here for provenance only.",
      source_name        = "RACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 165,
    n_studies      = 5,
    age_median     = "45 years",
    age_range      = "19-81 years",
    sex_female_pct = 56,
    race_ethnicity = c(White = 39, Asian = 38),
    disease_state  = "adults with paroxysmal nocturnal hemoglobinuria; 54% receiving eculizumab at baseline (PHAROAH, PEGASUS) and 46% complement C5-inhibitor naive (PADDOCK, PALOMINO, PRINCE)",
    dose_range     = "Subcutaneous 25-1080 mg, once daily or twice weekly (including the approved 1080 mg twice-weekly regimen)",
    regions        = "Multinational; the PRINCE cohort is predominantly Asian (36/50, 72%)",
    notes          = "Crass 2024 Sect. 3.1.2 and Table 1. 3423 lactate-dehydrogenase samples were analysed. The upper limit of normal used by the phase 3 reference laboratory was 226 U/L for the 1.5 x ULN threshold quoted in Sect. 3.4.3. The PK/PD models were fitted sequentially, conditioned on individual empirical Bayes estimates of the PK parameters from the population PK model (Sect. 2.2)."
  )

  ini({
    # -----------------------------------------------------------------
    # PK layer carried forward from the companion population PK model
    # `Crass_2024_pegcetacoplan` (Crass 2024 ESM Table 3). The LDH
    # PK/PD control stream (ESM Table 1, "LDH Final Updated Model")
    # reads the individual PK parameters straight from the $INPUT
    # columns ICL / IV2 / IKA / IF1 -- the post-hoc empirical Bayes
    # estimates produced by the population PK run -- so none of these
    # values were re-estimated here. They are wrapped in fixed() to
    # record that. The IIV block is carried forward as well so the file
    # can be simulated stand-alone; see the vignette's Assumptions
    # section.
    lcl <- fixed(log(0.0117)); label("Clearance in a 70 kg non-PNH participant (L/h)")     # ESM Table 3, theta 1 (carried from the population PK model)
    lvc <- fixed(log(3.96)); label("Central volume of distribution in a 70 kg participant (L)") # ESM Table 3, theta 2 (carried from the population PK model)
    lka <- fixed(log(0.0370)); label("First-order transit / absorption rate constant (1/h)") # ESM Table 3, theta 3 (carried from the population PK model)
    lfdepot <- fixed(log(0.758)); label("Subcutaneous bioavailability, solution formulations (fraction)") # ESM Table 3, theta 4 (carried from the population PK model)
    e_wt_cl <- fixed(0.646); label("Body-weight power exponent on clearance (unitless)")   # ESM Table 3, theta 10 (carried from the population PK model)
    e_wt_vc <- fixed(0.809); label("Body-weight power exponent on central volume (unitless)") # ESM Table 3, theta 11 (carried from the population PK model)
    e_dis_pnh_cl <- fixed(0.257); label("Fractional change in clearance in PNH patients (unitless)") # ESM Table 3, theta 9 (carried from the population PK model)
    e_form_pegcet_lyophilized_fdepot <- fixed(0.220); label("Fractional change in SC bioavailability, lyophilized powder (unitless)") # ESM Table 3, theta 8 (carried from the population PK model)

    # -----------------------------------------------------------------
    # Lactate-dehydrogenase PK/PD structural parameters. ESM Table 5
    # "PK/PD Parameter Estimates for the LDH PK/PD Model" reports the
    # estimates on their transformed (log / logit) scales together with
    # the back-transformed values; the back-transforms are entered here
    # in their transformed form so the model reproduces the source
    # control stream verbatim.
    #   exp(7.56)      = 1920 U/L   (ESM Table 5 "Transformed Estimate")
    #   exp(5.52)      =  249 U/L
    #   ilogit(2.40)   = 0.917
    #   ilogit(-1.39)  = 0.200
    #   exp(5.23)      =  187 ug/mL
    lrbase <- 7.56; label("Log baseline LDH, eculizumab-naive (log U/L)")                  # ESM Table 5, theta 1 "Log BaseLDH" 7.56 (95% CI 7.46-7.67); transformed 1920 U/L (1740-2140); Table 2
    lrbase_ecu_bl <- 5.52; label("Log baseline LDH, eculizumab-treated at baseline (log U/L)") # ESM Table 5, theta 2 "Log BaseLDH BECU" 5.52 (95% CI 5.42-5.61); transformed 249 U/L (226-273); Table 2
    logitemax <- 2.40; label("Logit maximal fractional LDH suppression, eculizumab-naive (logit units)") # ESM Table 5, theta 3 "Logit Emax" 2.40 (95% CI 2.17-2.63); transformed 0.917 (0.898-0.933); Table 2 "-92 (-93 to -90)"
    logitemax_ecu_bl <- -1.39; label("Logit maximal fractional LDH suppression, eculizumab-treated at baseline (logit units)") # ESM Table 5, theta 4 "Logit Emax BECU" -1.39 (95% CI -1.65 to -1.12); transformed 0.200 (0.161-0.246); Table 2 "-20 (-25 to -16)"
    lec50 <- 5.23; label("Log pegcetacoplan concentration producing half the maximal LDH response (log ug/mL)") # ESM Table 5, theta 5 "Log EC50" 5.23 (95% CI 5.11-5.35); transformed 187 ug/mL (166-211); Table 2
    lhill <- log(3.42); label("Hill sigmoidicity coefficient for the LDH response (unitless)") # ESM Table 5, theta 6 Hill 3.42 (95% CI 2.84-4.00), reported untransformed; Table 2

    # -----------------------------------------------------------------
    # Covariate effect: concurrent eculizumab co-treatment, applied
    # additively on the logit scale and only within the
    # baseline-eculizumab stratum.
    e_conmed_eculizumab_emax <- 0.783; label("Additive shift in maximal LDH suppression during eculizumab co-treatment (logit units)") # ESM Table 5, theta 8 "Eculizumab Co-treatment on Logit Emax" 0.783 (95% CI 0.605-0.961)

    # -----------------------------------------------------------------
    # Inter-individual variability, PK layer: carried forward from the
    # population PK model exactly as in `Crass_2024_pegcetacoplan`
    # (variances recovered from the ESM Table 3 CV% rows via
    # omega^2 = log(CV^2 + 1); the CL-Vc covariance from rho = 0.571).
    etalcl + etalvc ~ fixed(c(
      0.042754,
      0.024984, 0.044778
    )) # ESM Table 3 IIV rows ETA1-CL 20.9%, rho 0.571, ETA2-V2 21.4% (carried from the population PK model)
    etalka ~ fixed(0.243436) # ESM Table 3 IIV row ETA3-KA 52.5% (carried from the population PK model)

    # Inter-individual variability, PD layer: a single $OMEGA BLOCK(3)
    # over (BaseLDH, logit Emax, EC50). ESM Table 5 reports variances
    # and covariances directly; the reported correlations reproduce as
    # 0.265 / sqrt(0.182 * 0.745) = 0.720 and
    # 0.167 / sqrt(0.745 * 0.121) = 0.556. The remaining block element
    # omega3,1 (EC50:BaseLDH) is absent from ESM Table 5 and is entered
    # as 0 in the ESM Table 1 control stream's $OMEGA BLOCK(3), so it is
    # encoded as 0 here; the resulting 3x3 matrix is positive definite
    # (eigenvalues 0.882, 0.144, 0.022).
    etalrbase + etalogitemax + etalec50 ~ c(
      0.182,
      0.265, 0.745,
      0, 0.167, 0.121
    ) # ESM Table 5: omega1,1^2 BaseLDH 0.182 (44.7 CV%), omega2,1 0.265 (rho 0.719), omega2,2^2 Emax 0.745, omega3,1 not reported (0 in the control stream), omega3,2 0.167 (rho 0.556), omega3,3^2 EC50 0.121 (35.9 CV%)

    # -----------------------------------------------------------------
    # Residual error. The source $ERROR block for LDH is
    # log-transformed-both-sides (IPRED = LOG(EFF); Y = IPRED + W*ERR(1)
    # with $SIGMA 1 FIXED), so THETA(7) is a standard deviation on the
    # natural-log scale and maps to nlmixr2's lnorm() error model.
    expSd <- 0.305; label("Log-scale residual SD on LDH")                                  # ESM Table 5, theta 7 "Log additive" 0.305 (95% CI 0.297-0.313), transformed 30.5% (29.7-31.3%)
  })

  model({
    # -----------------------------------------------------------------
    # PK layer (identical to `Crass_2024_pegcetacoplan`). ESM Table 1
    # LDH control stream $DES reproduces the population PK $DES exactly
    # and adds an inert fourth compartment, DADT(4) = 0, which never
    # enters the prediction; the LDH response is computed algebraically
    # in $ERROR. That inert state is therefore not carried here.
    cl <- exp(lcl + etalcl) * (1 + e_dis_pnh_cl * DIS_PNH) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central

    f(depot) <- exp(lfdepot) *
      (1 + e_form_pegcet_lyophilized_fdepot * FORM_PEGCET_LYOPHILIZED)

    # Pegcetacoplan serum concentration in ug/mL; the driver of the LDH
    # response. Reported as an output column but not fitted in this run
    # (the LDH control stream sets IGNORE=(TYPE.EQ.2) so pegcetacoplan
    # concentration records are excluded).
    Cc <- central / vc

    # -----------------------------------------------------------------
    # LDH PK/PD layer. ESM Table 1 LDH control stream:
    #   TVBLDH = THETA(1)
    #   IF(BECU.EQ.1) TVBLDH = THETA(2)
    #   TVEMAX = THETA(3)
    #   IF (BECU.EQ.1) TVEMAX = THETA(4) + THETA(8)*ECU
    #   BASELDH = EXP(TVBLDH + ETA(1))
    #   EMAX    = EXP(TVEMAX+ETA(2))/(EXP(TVEMAX+ETA(2))+1)
    #   EC50    = EXP(TVEC50 + ETA(3))
    #   EDRUG   = (EMAX*CP**HILL)/(EC50**HILL+CP**HILL)
    #   EFF     = BASELDH*(1-EDRUG)
    # which is the paper's Sect. 3.2.2 equation with a negative Emax:
    #   Biomarker = Baseline * (1 + Emax*C^Hill/(EC50^Hill + C^Hill)).
    # The two IF-guarded lines are switches (a full replacement of the
    # typical value), not additive shifts, so they are encoded as an
    # indicator-weighted mixture of the two stratum parameters.
    lrbase_typ <- (1 - CONMED_ECULIZUMAB_BL) * lrbase +
      CONMED_ECULIZUMAB_BL * lrbase_ecu_bl
    logitemax_typ <- (1 - CONMED_ECULIZUMAB_BL) * logitemax +
      CONMED_ECULIZUMAB_BL *
        (logitemax_ecu_bl + e_conmed_eculizumab_emax * CONMED_ECULIZUMAB)

    rbase <- exp(lrbase_typ + etalrbase)
    emax <- exp(logitemax_typ + etalogitemax) /
      (exp(logitemax_typ + etalogitemax) + 1)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    edrug <- emax * Cc^hill / (ec50^hill + Cc^hill)
    ldh <- rbase * (1 - edrug)

    ldh ~ lnorm(expSd)
  })
}
