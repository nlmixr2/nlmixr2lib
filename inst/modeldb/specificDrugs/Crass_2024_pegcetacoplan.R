Crass_2024_pegcetacoplan <- function() {
  description <- "One-compartment population PK model for pegcetacoplan with a single-transit-compartment subcutaneous absorption chain and direct intravenous input, pooled across 11 studies in healthy adults, adults with renal impairment, and adults with paroxysmal nocturnal hemoglobinuria (Crass 2024). Subcutaneous dose enters a depot that transfers to one transit compartment and then to the central compartment, all at the same first-order rate ka; intravenous dose enters the central compartment directly. Clearance carries a fractional increase in patients with paroxysmal nocturnal hemoglobinuria, and both clearance and central volume carry estimated body-weight power exponents referenced to 70 kg. Subcutaneous bioavailability carries a fractional increase for the lyophilized-powder formulation relative to the ready-to-use solution formulations. Inter-individual variability is log-normal on clearance and central volume (correlated block) and on ka (diagonal). Residual error is additive on the log scale (log-normal) and is stratified across three paper-defined sub-populations (healthy participants, phase 1/2 paroxysmal nocturnal hemoglobinuria studies, phase 3 paroxysmal nocturnal hemoglobinuria studies), switched at runtime via the canonical DIS_PNH and STUDY_PEGCET_PHASE3 covariates."
  reference <- "Crass RL, Smith B, Adriaens S, Chapel S, Langdon G. Population Pharmacokinetic and Pharmacokinetic/Pharmacodynamic Analyses of Pegcetacoplan in Patients with Paroxysmal Nocturnal Hemoglobinuria. Drugs R D. 2024;24(4):565-577. doi:10.1007/s40268-024-00500-7"
  vignette <- "Crass_2024_pegcetacoplan"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  paper_specific_residual_sds <- c(
    "expSd_healthy", "expSd_pnh_ph12", "expSd_pnh_ph3"
  )

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
      notes              = "Time-fixed baseline body weight. Enters CL and Vc as estimated power exponents referenced to 70 kg (ESM Table 1 PK control stream: MWT=70; TVCL = THETA(1)*(1+PNH*THETA(9))*((WT/MWT)**THETA(10)); TVV2 = THETA(2)*((WT/MWT)**THETA(11))). The control stream substitutes the 70 kg reference whenever BWT is missing or non-positive. Analysis-set median 70.0 kg, range 41-156 kg (ESM Table 2, Total column); the 5th-95th percentile band cited in the paper's forest-plot discussion is 54-95 kg.",
      source_name        = "BWT"
    ),
    DIS_PNH = list(
      description        = "Paroxysmal nocturnal hemoglobinuria indicator: 1 = patient with PNH, 0 = healthy participant or renal-impairment participant.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-PNH participant: healthy adults from CP0713-1, CP1014, 101, 102, 401 and the healthy plus renal-impaired participants of AIRIS)",
      notes              = "Time-fixed per subject. Carries a fractional (linear) increase in clearance, `cl * (1 + e_dis_pnh_cl * DIS_PNH)` with e_dis_pnh_cl = 0.257 (ESM Table 1 PK control stream TVCL line; ESM Table 3 theta 9). Also selects which of the three residual-error strata applies. Analysis-set split 124 healthy (44%) / 160 PNH (56%) (ESM Table 2, Total column).",
      source_name        = "PNH"
    ),
    STUDY_PEGCET_PHASE3 = list(
      description        = "Phase 3 pegcetacoplan study cohort indicator: 1 = participant enrolled in PEGASUS (NCT03500549) or PRINCE (NCT04085601), 0 = any other study in the pooled analysis.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1 / phase 1b / phase 2a studies: CP0713-1, CP1014, 101, 102, 401, AIRIS, PHAROAH, PADDOCK, PALOMINO)",
      notes              = "Time-fixed per subject. Used only in the residual-error switch, never in the structural model. The source control stream builds two study-group flags, PNH3 (STUD 302 or 308, i.e. the phase 3 studies PEGASUS and PRINCE) and PNH2 (STUD 202, 204, or 514, i.e. the phase 1/2 PNH studies PHAROAH, PADDOCK, PALOMINO), and selects W from THETA(5) / THETA(6) / THETA(7) accordingly. Because PNH2 and PNH3 partition the PNH cohort exactly, the two source flags are re-expressed here as DIS_PNH x STUDY_PEGCET_PHASE3.",
      source_name        = "(derived from STUD; STUD 302 = PEGASUS and STUD 308 = PRINCE map to 1)"
    ),
    FORM_PEGCET_LYOPHILIZED = list(
      description        = "Pegcetacoplan lyophilized-powder formulation indicator: 1 = lyophilized powder reconstituted before subcutaneous administration, 0 = a ready-to-use solution formulation (sorbitol, dextrose, or mannitol vehicle).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ready-to-use solution formulations; source FORM levels 1 = sorbitol, 2 = dextrose, 3 = mannitol)",
      notes              = "Per-regimen categorical indicator. Carries a fractional (linear) increase in subcutaneous bioavailability, `fdepot * (1 + e_form_pegcet_lyophilized_fdepot * FORM_PEGCET_LYOPHILIZED)` with coefficient 0.220 (ESM Table 1 PK control stream `F1 = THETA(4)*(1+FORM4*THETA(8))`; ESM Table 3 theta 8, labelled 'Lyophilized formulation on F1'). Applies only to the subcutaneous depot; intravenous doses enter the central compartment directly and are not scaled.",
      source_name        = "FORM (level 4 = POWDER)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 284,
    n_studies      = 11,
    age_median     = "36 years",
    age_range      = "19-81 years",
    weight_median  = "70.0 kg",
    weight_range   = "41-156 kg",
    sex_female_pct = 44.4,
    race_ethnicity = c(White = 54.2, Black = 2.5, Asian = 29.9, Other = 7.7, Missing = 5.6),
    disease_state  = "44% healthy adults (including one renal-impairment study); 56% adults with paroxysmal nocturnal hemoglobinuria, both complement C5-inhibitor naive and C5-inhibitor experienced",
    dose_range     = "Subcutaneous 25-2600 mg as single doses, once daily, twice weekly, or once weekly (including the approved 1080 mg twice-weekly regimen); intravenous 200-2300 mg single doses (study 401)",
    regions        = "Multinational; the PRINCE phase 3 cohort is predominantly Asian (32/45, 71%)",
    notes          = "Baseline demographics from Crass 2024 ESM Table 2 (categorical and continuous covariates by study) and Table 1 (study designs and dosing regimens). Median baseline complement C3 across all participants 1.00 g/L (range 0.47-1.64 g/L). 5195 PK samples were analysed, of which 4737 (91%) were quantifiable and 458 (9%) were below the lower limit of quantification."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters. Reference subject: 70 kg, non-PNH
    # (healthy), ready-to-use solution formulation.
    # ESM Table 3 "Parameter Estimates for the PK Model".
    lcl <- log(0.0117); label("Clearance in a 70 kg non-PNH participant (L/h)")            # ESM Table 3, theta 1 TVCL 0.0117 L/h (95% CI 0.0106-0.0128)
    lvc <- log(3.96); label("Central volume of distribution in a 70 kg participant (L)")   # ESM Table 3, theta 2 TVV2 3.96 L (95% CI 3.57-4.34)
    lka <- log(0.0370); label("First-order transit / absorption rate constant (1/h)")      # ESM Table 3, theta 3 TVKA 0.0370 1/h (95% CI 0.0343-0.0397)
    lfdepot <- log(0.758); label("Subcutaneous bioavailability, solution formulations (fraction)") # ESM Table 3, theta 4 F1 0.758 (95% CI 0.681-0.834)

    # -----------------------------------------------------------------
    # Covariate effects. All four are estimated (no FIX flag in the ESM
    # Table 1 control stream; all four carry an ASE and a 95% CI in
    # ESM Table 3).
    e_wt_cl <- 0.646; label("Body-weight power exponent on clearance (unitless)")          # ESM Table 3, theta 10 "WT on CL" 0.646 (95% CI 0.515-0.776)
    e_wt_vc <- 0.809; label("Body-weight power exponent on central volume (unitless)")     # ESM Table 3, theta 11 "WT on V2" 0.809 (95% CI 0.643-0.975)
    e_dis_pnh_cl <- 0.257; label("Fractional change in clearance in PNH patients (unitless)") # ESM Table 3, theta 9 "PNH on CL" 0.257 (95% CI 0.197-0.316)
    e_form_pegcet_lyophilized_fdepot <- 0.220; label("Fractional change in SC bioavailability, lyophilized powder (unitless)") # ESM Table 3, theta 8 "Lyophilized formulation on F1" 0.220 (95% CI 0.145-0.296)

    # -----------------------------------------------------------------
    # Inter-individual variability. ESM Table 3 reports IIV as CV%;
    # variances are recovered as omega^2 = log(CV^2 + 1) for the
    # log-normal etas, and the CL-Vc covariance from the reported
    # correlation rho = 0.571 as rho * sqrt(var_cl * var_vc).
    #   CL:  CV 20.9% -> log(0.209^2 + 1) = 0.042754
    #   Vc:  CV 21.4% -> log(0.214^2 + 1) = 0.044778
    #   cov: 0.571 * sqrt(0.042754 * 0.044778) = 0.024984
    #   ka:  CV 52.5% -> log(0.525^2 + 1) = 0.243436
    etalcl + etalvc ~ c(
      0.042754,
      0.024984, 0.044778
    ) # ESM Table 3 IIV rows: ETA1-CL 20.9% (18.8-22.9), rho(ETA1,ETA2) 0.571, ETA2-V2 21.4% (18.6-23.8); $OMEGA BLOCK(2)
    etalka ~ 0.243436 # ESM Table 3 IIV row ETA3-KA 52.5% (46.0-58.5); separate diagonal $OMEGA

    # -----------------------------------------------------------------
    # Residual error. The source $ERROR block is log-transformed-both-
    # sides: IPRED = LOG(F); Y = IPRED + W*ERR(1) with $SIGMA 1 FIXED,
    # so each W is a standard deviation on the natural-log scale and
    # maps to nlmixr2's lnorm() error model. W is switched across three
    # paper-defined sub-populations by the source PNH / PNH2 / PNH3
    # flags. All three names are declared in
    # `paper_specific_residual_sds` because the canonical residual-error
    # matcher recognises only bare `expSd` or per-output `expSd_<obs>`
    # suffixes, not per-cohort suffixes.
    expSd_healthy <- 0.200; label("Log-scale residual SD, healthy participants")           # ESM Table 3, theta 5 "RE healthy subjects" 20.0% (19.4-20.6)
    expSd_pnh_ph12 <- 0.326; label("Log-scale residual SD, phase 1/2 PNH studies")         # ESM Table 3, theta 6 "RE PNH Phase 1 and 2" 32.6% (30.7-34.6)
    expSd_pnh_ph3 <- 0.163; label("Log-scale residual SD, phase 3 PNH studies")            # ESM Table 3, theta 7 "RE PNH Phase 3" 16.3% (15.7-16.9)
  })

  model({
    # -----------------------------------------------------------------
    # Individual parameters. The body-weight terms use the source
    # control stream's 70 kg reference (MWT = 70). The PNH effect on
    # clearance is fractional-linear, exactly as written in the source
    # TVCL line, not a power or exponential form.
    cl <- exp(lcl + etalcl) * (1 + e_dis_pnh_cl * DIS_PNH) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    kel <- cl / vc

    # -----------------------------------------------------------------
    # ODE system. ESM Table 1 PK control stream $DES, with the source
    # compartment numbering DOSESC1 = 1, CENTRAL = 2, TRANSIT1 = 3:
    #   DADT(1) = -KA*A(1)
    #   DADT(2) =  KA*A(3) - (CL/V2)*A(2)
    #   DADT(3) =  KA*A(1) - KA*A(3)
    # i.e. a subcutaneous depot feeding one transit compartment which
    # feeds the central compartment, every transfer governed by the
    # same rate constant ka. Intravenous doses (study 401) are given
    # directly into `central`.
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central

    # Bioavailability applies to the subcutaneous depot only; the
    # intravenous route enters `central` and is not scaled.
    f(depot) <- exp(lfdepot) *
      (1 + e_form_pegcet_lyophilized_fdepot * FORM_PEGCET_LYOPHILIZED)

    # Serum concentration. The source sets S2 = V2 with dose in mg and
    # concentration in mcg/mL = mg/L, so central (mg) / vc (L) is
    # already in ug/mL and no scaling factor is needed.
    Cc <- central / vc

    # -----------------------------------------------------------------
    # Residual error switch. Exactly one of the three indicator
    # products is 1 for any given subject: non-PNH participants take
    # the healthy stratum, and PNH patients take the phase 3 stratum
    # when enrolled in PEGASUS / PRINCE and the phase 1/2 stratum
    # otherwise.
    expSd <- expSd_healthy * (1 - DIS_PNH) +
      DIS_PNH * (STUDY_PEGCET_PHASE3 * expSd_pnh_ph3 +
        (1 - STUDY_PEGCET_PHASE3) * expSd_pnh_ph12)

    Cc ~ lnorm(expSd)
  })
}
