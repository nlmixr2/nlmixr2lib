Crass_2024_pegcetacoplan_hemoglobin <- function() {
  description <- "Direct sigmoidal Emax PK/PD model relating pegcetacoplan serum concentration to hemoglobin in adults with paroxysmal nocturnal hemoglobinuria (Crass 2024). Hemoglobin is a direct (no-delay) function of the concurrent pegcetacoplan concentration, hb = rbase * (1 + emax * Cc^hill / (ec50^hill + Cc^hill)), with the pegcetacoplan concentration supplied by the companion one-compartment-plus-single-transit population PK model of the same paper, whose parameters are carried here as fixed values because the PK/PD model was fitted sequentially on individual empirical Bayes estimates of the PK parameters. The maximal proportional increase in hemoglobin carries a power effect of baseline creatinine clearance and a fractional effect of female sex; the baseline-eculizumab effect on the maximal response was fixed to zero in the source control stream. Inter-individual variability is log-normal on the hemoglobin baseline and the maximal effect (correlated block) and on EC50 (diagonal); residual error on hemoglobin is additive."
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
      notes              = "Time-fixed baseline body weight. Enters only the carried-forward PK layer (power exponents on CL and Vc referenced to 70 kg). Body weight was screened as a covariate on the hemoglobin PK/PD parameters and was not retained (Crass 2024 Sect. 3.4.2: 'There were no meaningful effects on Hb with all other covariates (i.e., prior eculizumab treatment, race, age, body weight, and baseline C3 level)').",
      source_name        = "BWT"
    ),
    DIS_PNH = list(
      description        = "Paroxysmal nocturnal hemoglobinuria indicator: 1 = patient with PNH, 0 = non-PNH participant.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-PNH participant)",
      notes              = "Time-fixed per subject. Enters only the carried-forward PK layer, where it applies the fractional clearance increase e_dis_pnh_cl = 0.257. Every subject in the hemoglobin PK/PD analysis set is a patient with PNH (DIS_PNH = 1), so the term is inert in this cohort; it is retained so the PK layer is byte-for-byte the same structure as the parent model `Crass_2024_pegcetacoplan`.",
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
    CRCL = list(
      description        = "Baseline creatinine clearance estimated with the Cockcroft-Gault equation, in mL/min (NOT body-surface-area normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Enters the maximal hemoglobin response as a power term referenced to 120 mL/min (ESM Table 1 Hb control stream: MCRCL=120; TVEMAX = THETA(2)*((CRCL/MCRCL)**THETA(6))), with exponent 0.641. The control stream substitutes the 120 mL/min reference whenever BCRCL is missing or non-positive. Crass 2024 ESM Table 4 footnote defines CrCl as 'creatinine clearance calculated using the Cockcroft-Gault equation', which returns raw mL/min rather than mL/min/1.73 m^2; the canonical CRCL column accommodates both forms (see the register entry's Delattre 2010 / Chen 2023 / Wada 2023 precedents for raw Cockcroft-Gault use). The paper cites the 5th percentile as 40 mL/min and the 95th percentile as 191 mL/min in this cohort (Sect. 3.4.2).",
      source_name        = "BCRCL"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Applies a fractional change to the maximal hemoglobin response, `emax * (1 + e_sexf_emax * SEXF)` with e_sexf_emax = -0.337 (ESM Table 1 Hb control stream TVEMAX line; ESM Table 4 theta 7). The source column is already coded SEXF (1 = female), matching the canonical orientation directly, so no value transformation is needed. At the reference CrCl of 120 mL/min this gives a 33.8% maximal increase in females versus 51.0% in males (Crass 2024 Sect. 3.4.2). The PK/PD analysis set is 56% female (Sect. 3.1.2).",
      source_name        = "SEXF"
    ),
    CONMED_ECULIZUMAB_BL = list(
      description        = "Eculizumab (complement C5 inhibitor) treatment status at baseline: 1 = receiving eculizumab at study entry, 0 = eculizumab-naive at study entry.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (eculizumab-naive at baseline; the PADDOCK, PALOMINO, and PRINCE cohorts)",
      notes              = "Time-fixed per subject. Present in the source TVEMAX line as `TVEMAX = TVEMAX*(1+THETA(7)*SEXF)*(1+THETA(8)*BECU)` with THETA(8) entered as `(0 FIX)` in ESM Table 1, so the baseline-eculizumab effect on the maximal hemoglobin response is fixed at zero and is not reported in ESM Table 4. The term is retained in the model with `e_conmed_eculizumab_bl_emax <- fixed(0)` so the structural form matches the published control stream and the fixed status is visible; it is what makes Crass 2024 Table 2 report the same 51% Emax for eculizumab-naive and eculizumab-treated patients. The PK/PD analysis set is 54% eculizumab-treated / 46% eculizumab-naive at baseline (Sect. 3.1.2).",
      source_name        = "BECU"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the full covariate model on the hemoglobin PK/PD parameters and eliminated by the backward-elimination step (alpha = 0.001); no effect is reported in ESM Table 4 or Crass 2024 Table 2. Crass 2024 Sect. 3.4.2 states there were no meaningful effects on Hb from age. Analysis-set median 45 years (range 19-81 years).",
      source_name        = "AGE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "The source Hb control stream derives RASN (RACE 3 or 4) and JPN (RACE 4) in its $PK block, but neither enters the final TVEMAX / TVEC50 / TVBASEHGB expressions; race was screened and not retained (Crass 2024 Sect. 3.4.2). The PK/PD analysis set is 39% White and 38% Asian (Sect. 3.1.2). Recorded here for provenance only.",
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
    notes          = "Crass 2024 Sect. 3.1.2 and Table 1. 4498 hemoglobin samples were analysed. Median baseline complement C3 across the PNH studies was 0.910 g/L (range 0.47-1.64 g/L). The PK/PD models were fitted sequentially, conditioned on individual empirical Bayes estimates of the PK parameters from the population PK model (Sect. 2.2)."
  )

  ini({
    # -----------------------------------------------------------------
    # PK layer carried forward from the companion population PK model
    # `Crass_2024_pegcetacoplan` (Crass 2024 ESM Table 3). The Hb PK/PD
    # control stream (ESM Table 1, "Hb Final Updated Model") reads the
    # individual PK parameters straight from the $INPUT columns
    # ICL / IV2 / IKA / IF1 -- the post-hoc empirical Bayes estimates
    # produced by the population PK run -- so none of these values were
    # re-estimated here. They are wrapped in fixed() to record that.
    # The IIV block is carried forward as well so the file can be
    # simulated stand-alone; see the vignette's Assumptions section.
    lcl <- fixed(log(0.0117)); label("Clearance in a 70 kg non-PNH participant (L/h)")     # ESM Table 3, theta 1 (carried from the population PK model)
    lvc <- fixed(log(3.96)); label("Central volume of distribution in a 70 kg participant (L)") # ESM Table 3, theta 2 (carried from the population PK model)
    lka <- fixed(log(0.0370)); label("First-order transit / absorption rate constant (1/h)") # ESM Table 3, theta 3 (carried from the population PK model)
    lfdepot <- fixed(log(0.758)); label("Subcutaneous bioavailability, solution formulations (fraction)") # ESM Table 3, theta 4 (carried from the population PK model)
    e_wt_cl <- fixed(0.646); label("Body-weight power exponent on clearance (unitless)")   # ESM Table 3, theta 10 (carried from the population PK model)
    e_wt_vc <- fixed(0.809); label("Body-weight power exponent on central volume (unitless)") # ESM Table 3, theta 11 (carried from the population PK model)
    e_dis_pnh_cl <- fixed(0.257); label("Fractional change in clearance in PNH patients (unitless)") # ESM Table 3, theta 9 (carried from the population PK model)
    e_form_pegcet_lyophilized_fdepot <- fixed(0.220); label("Fractional change in SC bioavailability, lyophilized powder (unitless)") # ESM Table 3, theta 8 (carried from the population PK model)

    # -----------------------------------------------------------------
    # Hemoglobin PK/PD structural parameters. ESM Table 4 "Parameter
    # Estimates for the Hb PK/PD Model"; the same four values appear in
    # Crass 2024 Table 2 (Hb column).
    lrbase <- log(8.74); label("Baseline hemoglobin (g/dL)")                               # ESM Table 4, theta 1 BaseHb 8.74 g/dL (95% CI 8.55-8.93); Table 2 "Baseline concentration"
    lemax <- log(0.510); label("Maximal proportional increase in hemoglobin at CrCl 120 mL/min in males (fraction)") # ESM Table 4, theta 2 Emax 0.510 (95% CI 0.431-0.589); Table 2 "Emax, % change from baseline" 51 (43-59)
    lec50 <- log(337); label("Pegcetacoplan concentration producing half the maximal hemoglobin response (ug/mL)") # ESM Table 4, theta 3 EC50 337 ug/mL (95% CI 302-372); Table 2
    lhill <- log(4.66); label("Hill sigmoidicity coefficient for the hemoglobin response (unitless)") # ESM Table 4, theta 4 Hill 4.66 (95% CI 3.81-5.52); Table 2

    # -----------------------------------------------------------------
    # Covariate effects on the maximal hemoglobin response.
    e_crcl_emax <- 0.641; label("Power exponent of baseline creatinine clearance on Emax (unitless)") # ESM Table 4, theta 6 "CrCl on Emax" 0.641 (95% CI 0.433-0.850)
    e_sexf_emax <- -0.337; label("Fractional change in Emax for female sex (unitless)")    # ESM Table 4, theta 7 "Female sex on Emax" -0.337 (95% CI -0.457 to -0.218)
    # THETA(8) is entered as `(0 FIX)` in the ESM Table 1 Hb control
    # stream and is absent from ESM Table 4, which is why Crass 2024
    # Table 2 reports an identical 51% Emax for eculizumab-naive and
    # eculizumab-treated patients.
    e_conmed_eculizumab_bl_emax <- fixed(0); label("Fractional change in Emax for baseline eculizumab treatment (unitless)") # ESM Table 1 Hb control stream $THETA line 8 "(0 FIX) ;8 BECU ON EMAX"

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

    # Inter-individual variability, PD layer. ESM Table 4 reports these
    # as variances directly (the parenthesised CV% values are the
    # back-transforms, and the reported rho = -0.580 reproduces as
    # -0.0437 / sqrt(0.0162 * 0.351)). The source $OMEGA is a BLOCK(2)
    # over (BaseHb, Emax) plus a separate diagonal $OMEGA for EC50.
    etalrbase + etalemax ~ c(
      0.0162,
      -0.0437, 0.351
    ) # ESM Table 4: omega1,1^2 BaseHb 0.0162 (12.8 CV%), omega2,1 Emax:BaseHb -0.0437 (rho -0.580), omega2,2^2 Emax 0.351 (64.8 CV%)
    etalec50 ~ 0.291 # ESM Table 4: omega3,3^2 EC50 0.291 (58.1 CV%); separate diagonal $OMEGA

    # -----------------------------------------------------------------
    # Residual error. The source $ERROR block for Hb is on the natural
    # scale (IPRED = EFF, not LOG(EFF)) with $SIGMA 1 FIXED, so THETA(5)
    # is an additive standard deviation in g/dL.
    addSd <- 0.913; label("Additive residual SD on hemoglobin (g/dL)")                     # ESM Table 4, theta 5 "RE Additive (g/dL)" 0.913 (95% CI 0.888-0.938)
  })

  model({
    # -----------------------------------------------------------------
    # PK layer (identical to `Crass_2024_pegcetacoplan`). ESM Table 1
    # Hb control stream $DES reproduces the population PK $DES exactly
    # and adds an inert fourth compartment, DADT(4) = 0, which never
    # enters the prediction; the hemoglobin response is computed
    # algebraically in $ERROR. That inert state is therefore not
    # carried here.
    cl <- exp(lcl + etalcl) * (1 + e_dis_pnh_cl * DIS_PNH) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central

    f(depot) <- exp(lfdepot) *
      (1 + e_form_pegcet_lyophilized_fdepot * FORM_PEGCET_LYOPHILIZED)

    # Pegcetacoplan serum concentration in ug/mL; the driver of the
    # hemoglobin response. Reported as an output column but not fitted
    # in this run (the Hb control stream sets IGNORE=(TYPE.EQ.2) so
    # pegcetacoplan concentration records are excluded).
    Cc <- central / vc

    # -----------------------------------------------------------------
    # Hemoglobin PK/PD layer. ESM Table 1 Hb control stream:
    #   TVEMAX = THETA(2)*((CRCL/MCRCL)**THETA(6))
    #   TVEMAX = TVEMAX*(1+THETA(7)*SEXF)*(1+THETA(8)*BECU)
    #   BASEHGB = TVBASEHGB*EXP(ETA(1))
    #   EMAX    = TVEMAX*EXP(ETA(2))
    #   EC50    = TVEC50*EXP(ETA(3))
    #   EDRUG   = (EMAX*CP**HILL)/(EC50**HILL+CP**HILL)
    #   EFF     = BASEHGB*(1+EDRUG)
    # which is the paper's Sect. 3.2.2 equation
    #   Biomarker = Baseline Biomarker * (1 + Emax*C^Hill/(EC50^Hill + C^Hill)).
    # The creatinine-clearance reference is MCRCL = 120 mL/min.
    rbase <- exp(lrbase + etalrbase)
    emax <- exp(lemax + etalemax) *
      (CRCL / 120)^e_crcl_emax *
      (1 + e_sexf_emax * SEXF) *
      (1 + e_conmed_eculizumab_bl_emax * CONMED_ECULIZUMAB_BL)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    edrug <- emax * Cc^hill / (ec50^hill + Cc^hill)
    hb <- rbase * (1 + edrug)

    hb ~ add(addSd)
  })
}
