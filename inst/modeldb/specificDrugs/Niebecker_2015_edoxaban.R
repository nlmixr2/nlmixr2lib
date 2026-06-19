Niebecker_2015_edoxaban <- function() {
  description <- "Two-compartment population PK model with first-order absorption and a lag time for edoxaban in adults; pooled phase 1 healthy volunteers (13 studies) and Hokusai-VTE phase 3 patients with deep-vein thrombosis or pulmonary embolism (Niebecker 2015). Apparent clearance is split into a non-renal component and a piecewise-linear renal component driven by creatinine clearance, with a phase-3 patient effect on the upper-CLcr slope and on Q/F. Asian race increases Vc/F; concomitant P-glycoprotein inhibitors increase phase-1 CL/F and F. The fed-state study 6 has a slower ka and higher non-renal CL/F (FED covariate)."
  reference <- paste(
    "Niebecker R, Jonsson S, Karlsson MO, Miller R, Nyberg J, Krekels EHJ, Simonsson USH.",
    "Population pharmacokinetics of edoxaban in patients with symptomatic deep-vein",
    "thrombosis and/or pulmonary embolism - the Hokusai-VTE phase 3 study.",
    "Br J Clin Pharmacol. 2015;80(6):1374-1387.",
    "doi:10.1111/bcp.12727.",
    sep = " "
  )
  vignette <- "Niebecker_2015_edoxaban"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL/F, Vc/F, Vp/F, Q/F with reference weight 70 kg. Allometric exponents are fixed per Niebecker 2015 Table 3 footnote (paragraph mark mark): CL/F = (WT/70)^(3/4), Vc/F = (WT/70)^1, Vp/F = (WT/70)^(3/4), Q/F = (WT/70)^1. The Vp/F-at-3/4 and Q/F-at-1 assignment is published as printed; note that this is the opposite of the more common volume-at-1 / clearance-at-3/4 grouping. Hokusai-VTE patient weight: 60-108 kg (10th-90th percentile), median 80.5 kg.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance, raw Cockcroft-Gault (mL/min), NOT BSA-normalized.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed via Cockcroft-Gault per Niebecker 2015 Methods reference [21] (Cockcroft & Gault 1976). Truncated at 150 mL/min in the source analysis: values above 150 are clamped to 150 before entering the piecewise-linear effect on the renal component of CL/F. The piecewise breakpoint is at CRCL = 90 mL/min: slope1 = 0.202 L/h/(mL/min) below 90; slope2 = 0.0321 L/h/(mL/min) above 90 in the phase 1 healthy-volunteer pool, and slope2 = 0.0321 * (1 + 2.74) = 0.120 L/h/(mL/min) above 90 in the Hokusai-VTE phase 3 patient cohort. Hokusai-VTE patient CRCL: 57.5-151 mL/min (10th-90th percentile), median 99 mL/min, range 14.1-247.",
      source_name        = "CLcr"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = non-Asian).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Multiplicative effect on Vc/F: Vc/F is 22.6% higher in Asian patients than non-Asians (Niebecker 2015 Table 3 final model: theta_Asian = 0.226). The source paper dichotomized race after finding that the clinically significant difference was between Asian and non-Asian subjects only (Results, 'Concerning the impact of race on Vc/F, a clinically significant difference was only found for Asians vs. non-Asians'). Asian fraction in Hokusai-VTE = 20.1% (740 of 3707 patients).",
      source_name        = "Asian race indicator"
    ),
    PGP_INH = list(
      description        = "Concomitant P-glycoprotein inhibitor coadministration indicator (1 = on a P-gp inhibitor, 0 = no concomitant P-gp inhibitor).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant P-gp inhibitor)",
      notes              = "Niebecker 2015 pools the following P-gp inhibitors into the indicator: verapamil, quinidine, dronedarone, erythromycin, azithromycin, clarithromycin, ketoconazole, itraconazole (the 'selected strong P-gp inhibitors' that mandated 50% dose reduction in Hokusai-VTE). Additional P-gp inhibitors recorded but not triggering dose reduction (amiodarone, captopril, carvedilol, conivaptan, ciclosporin, diltiazem, felodipine, quercetin, ranolazin, ticagrelor) were tested in covariate analysis but not retained. Effects on CL/F (+33.4%) and F (+125%) are applied ONLY to phase 1 healthy-volunteer subjects (STUDY_HOKVTE = 0) per Table 3 final-model column. The phase 3 effect on CL/F (~12% decrease) did not meet clinical significance and is excluded from the final model.",
      source_name        = "P-gp inhibitor co-administration"
    ),
    STUDY_HOKVTE = list(
      description        = "Hokusai-VTE phase 3 cohort indicator (1 = subject in Hokusai-VTE phase 3 VTE-patient cohort, 0 = subject in the pooled phase 1 healthy-volunteer studies).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1 healthy-volunteer pool: 13 studies, 443 subjects, 8652 PK observations)",
      notes              = "Hokusai-VTE: NCT00986154 phase 3 randomized double-blind double-dummy multicenter VTE study; 3707 patients with symptomatic DVT and/or PE receiving 60 mg edoxaban once daily (or 30 mg if dose-reduced for WT <= 60 kg, CLcr 30-50 mL/min, or concomitant P-gp inhibitor); 9531 PK observations sampled month 3 / month 12 / on-event. Used to switch (a) the upper-CLcr slope on the renal component of CL/F (+274% scaling relative to phase 1, so slope2_HV = 0.0321 -> slope2_phase3 = 0.120 L/h/(mL/min)) and (b) Q/F (+64.6% in patients vs healthy volunteers). Subject-level / time-fixed.",
      source_name        = "Study phase (phase 1 vs phase 3 Hokusai-VTE)"
    ),
    FED = list(
      description        = "Fed-state-at-dosing indicator (1 = administered with food, 0 = administered fasted).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "In Niebecker 2015 only study 6 (the dronedarone DDI crossover; Table 1 row 6) administered edoxaban under fed conditions; all other 12 phase 1 studies and the Hokusai-VTE phase 3 study used overnight-fast dosing. Two structural effects: ka is reduced 69% in the fed state (e_fed_ka = -0.690 in Niebecker's 'fractional change in ka study 6' notation) and apparent non-renal CL/F increases from 15.2 to 18.3 L/h (e_fed_cl_nonren = 0.204 = 18.3/15.2 - 1, in Niebecker's 'CLnr/F study 6' notation). The food-effect interpretation matches the Table 1 study-design column (Fed vs Overnight fast); the source paper used 'study 6' as the covariate name and did not separately mechanistically attribute the effect. For the targeted Hokusai-VTE phase 3 simulation population, FED = 0 always, so both effects contribute nothing.",
      source_name        = "Study 6 indicator (the only Fed-state phase 1 study)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 4130L,
    n_observations   = 17406L,
    n_studies        = 14L,
    age_range        = "Phase 1: 22.0-43.8 years (10th-90th percentile); Hokusai-VTE: 32.6-76.0 years (10th-90th); pooled mean 32 (range 18-67 across the analysis set; Hokusai-VTE mean 55.6).",
    age_median       = "Phase 1 median 30.0 years; Hokusai-VTE median 57.0 years",
    weight_range     = "Phase 1: 63.0-94.4 kg (10th-90th percentile); Hokusai-VTE: 60-108 kg (10th-90th); pooled mean 79 kg (range 50-111).",
    weight_median    = "Phase 1 median 79.3 kg; Hokusai-VTE median 80.5 kg",
    sex_female_pct   = 38,
    race_ethnicity   = c(White = 60.6, Black = 8.6, Asian = 18.1, Other = 12.4),
    disease_state    = "Hokusai-VTE: symptomatic deep-vein thrombosis (DVT) with or without pulmonary embolism (PE), or PE alone. Phase 1: healthy adult volunteers across 13 studies including a renal-impairment study and 5 P-gp DDI crossover studies.",
    dose_range       = "Hokusai-VTE: 60 mg orally once daily (or 30 mg if dose-reduced for WT <= 60 kg, CLcr 30-50 mL/min, or concomitant P-gp inhibitor); Phase 1: 15-60 mg single or multiple oral doses across the 13 studies.",
    regions          = "Hokusai-VTE: 439 centers across 37 countries; phase 1 studies in healthy-volunteer units in the US and Europe.",
    notes            = "Pooled population PK analysis. Demographics from Niebecker 2015 Table 2; race distribution computed from the pooled cohort percentages (Hokusai-VTE: 71.0% White, 3.23% Black, 20.1% Asian, 5.48% Other, 0.244% Missing; phase 1: 43.6% White, 51.7% Black, 2.03% Asian-non-Japanese, 2.48% Other; pooled fractions weighted by subject count). Missing race in 9 Hokusai-VTE patients imputed as White. Final analysis dataset: 17,406 plasma concentrations from 4,130 individuals after exclusion of 633 phase 1 and 144 Hokusai-VTE below-LLOQ samples (LLOQ 0.764 ng/mL by the Advion LC-MS/MS assay; 1 ng/mL by the BioDynamics LC-MS/MS assay used for the renal-impairment phase 1 study)."
  )

  ini({
    # Structural PK parameters -- Niebecker 2015 Table 3 final-model column.
    # All "/F" quantities are apparent values (oral dosing; bioavailability not
    # separately identifiable from CL and V on oral PK alone).
    lka         <- log(3.36);  label("Absorption rate constant ka (1/h)")                              # Table 3 final model: ka = 3.36 1/h (RSE 4.74%)
    lcl_nonren  <- log(15.2);  label("Apparent non-renal clearance CLnr/F at WT = 70 kg, fasted (L/h)") # Table 3 final model: CLnr/F = 15.2 L/h (RSE 2.20%)
    lvc         <- log(209);   label("Apparent central volume Vc/F at WT = 70 kg, non-Asian (L)")      # Table 3 final model: Vc/F = 209 L (RSE 1.61%)
    lvp         <- log(92.3);  label("Apparent peripheral volume Vp/F at WT = 70 kg (L)")              # Table 3 final model: Vp/F = 92.3 L (RSE 2.66%)
    lq          <- log(5.91);  label("Apparent intercompartmental clearance Q/F at WT = 70 kg, healthy volunteer (L/h)") # Table 3 final model: Q/F = 5.91 L/h (RSE 3.44%)
    ltlag       <- fixed(log(0.250)); label("Absorption lag time tlag (h)")                            # Table 3 final model: tlag = 0.250 h FIXED (Table 3 footnote sect: tlag fixed to phase 1 estimate)

    # Fixed bioavailability anchor (F = 1 typical-value; covariate effects from
    # PGP_INH are applied multiplicatively to F via f(depot) in model()).
    lfdepot     <- fixed(log(1)); label("Apparent bioavailability anchor F (fraction)")                # Implicit: typical F not separately identifiable on oral PK

    # Piecewise-linear creatinine-clearance effect on the renal component of CL/F.
    # Niebecker 2015 Table 3 footnote paragraph mark: Typical CL/F = CLnr/F + theta_Slope1 * CLcr
    # for CLcr <= 90 mL/min; for CLcr > 90 mL/min the second slope kicks in. In phase 3
    # (Hokusai-VTE patients) the second slope is scaled by (1 + theta_Scale) = (1 + 2.74) = 3.74.
    # CLcr is truncated at 150 mL/min before entering this formula.
    e_crcl_cl_renal_slope1 <- 0.202;  label("Renal CL slope below CLcr = 90 mL/min (L/h per mL/min, both cohorts)") # Table 3 final model: theta_Slope1 = 0.202 (RSE 2.22%)
    e_crcl_cl_renal_slope2 <- 0.0321; label("Renal CL slope above CLcr = 90 mL/min in phase 1 healthy volunteers (L/h per mL/min)") # Table 3 final model: theta_Slope2 = 0.0321 (RSE 4.74%)
    e_study_hokvte_cl_renal_slope2 <- 2.74; label("Fractional increase in renal CL slope2 for phase 3 patients vs phase 1 (unitless; +274%)") # Table 3 final model: 'Scaling parameter for slope 2 in phase 3, %' = 274 (RSE 8.02%)

    # Asian race effect on Vc/F (multiplicative): Vc/F_Asian = Vc/F_non-Asian * (1 + 0.226).
    e_race_asian_vc <- 0.226; label("Asian race fractional effect on Vc/F (unitless; +22.6% vs non-Asian)") # Table 3 final model: theta_Asian on Vc/F = 0.226 (RSE 13.6%; reported as 22.6%)

    # Phase 3 patient effect on Q/F (multiplicative): Q/F_phase3 = Q/F_phase1 * (1 + 0.646).
    e_study_hokvte_q <- 0.646; label("Phase 3 patient fractional effect on Q/F (unitless; +64.6% vs healthy volunteer)") # Table 3 final model: theta_Patients on Q/F = 0.646 (RSE 19.5%; reported as 64.6%)

    # Concomitant P-gp inhibitor effects -- estimated on the phase 1 subset only
    # (Table 3 final model row 'P-gp inhibitors on CL, phase 1' and 'P-gp inhibitors on F, phase 1').
    # Applied multiplicatively: CL/F_with_PGP = CL/F_without * (1 + 0.334) and F_with_PGP = F_without * (1 + 1.25).
    e_pgp_inh_cl <- 0.334; label("P-gp-inhibitor fractional effect on CL/F (phase 1 only; +33.4%)") # Table 3 final model: theta_P-gp on CL = 0.334 (RSE 9.39%; reported as 33.4%)
    e_pgp_inh_f  <- 1.25;  label("P-gp-inhibitor fractional effect on F (phase 1 only; +125%)")    # Table 3 final model: theta_P-gp on F = 1.25 (RSE 5.19%; reported as 125%)

    # Fed-state ("study 6") effects on ka and the non-renal component of CL/F.
    # Table 3 final model: ka in study 6 = ka * (1 + (-0.690)) = ka * 0.310; CLnr/F in study 6 = 18.3 L/h.
    # Encoded as multiplicative fractional changes on the typical values for fasted state, switched by FED.
    e_fed_ka         <- -0.690; label("Fed-state fractional effect on ka (unitless; -69.0% absorption slowdown with food)") # Table 3 final model: 'Fractional change in ka study 6' = -0.690 (RSE 1.19%)
    e_fed_cl_nonren  <-  0.204; label("Fed-state fractional effect on CLnr/F (unitless; +20.4% non-renal CL with food)")    # Table 3 final model: CLnr/F_study6 = 18.3 L/h vs CLnr/F = 15.2 L/h; (18.3/15.2 - 1) = 0.204

    # IIV scaling parameters -- the Niebecker 2015 random-effects structure uses
    # ONE eta per pair, with the second parameter of each pair receiving the eta
    # scaled by a multiplicative factor. theta_scale_cl_vc transfers etalcl from
    # CL/F to Vc/F (estimated at 1.56); theta_scale_vp_q transfers etalvp from
    # Vp/F to Q/F (fixed at 1.0, so Vp/F and Q/F share the same eta magnitude).
    # Table 3 final model footnote paragraph mark mark.
    theta_scale_cl_vc <- 1.56;        label("IIV scaling factor for etalcl applied to Vc/F (unitless)")  # Table 3 final model: theta_Scale1 = 1.56 (RSE 2.47%)
    theta_scale_vp_q  <- fixed(1.00); label("IIV scaling factor for etalvp applied to Q/F (unitless, FIXED)") # Table 3 final model: theta_Scale2 = 1.00 FIXED

    # Allometric exponents -- fixed per Table 3 footnote paragraph mark mark
    # (paper-as-printed): CL/F gets 3/4, Vc/F gets 1, Vp/F gets 3/4, Q/F gets 1.
    # The Vp-at-3/4 and Q-at-1 assignment is the opposite of the more common
    # volume-at-1 / clearance-at-3/4 grouping; the values are reproduced as
    # printed in the source paper.
    allo_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless, fixed at 3/4)")  # Table 3 footnote paragraph mark mark: CL/F (WT/70)^(3/4)
    allo_vc <- fixed(1.00); label("Allometric exponent on Vc/F (unitless, fixed at 1)")     # Table 3 footnote paragraph mark mark: Vc/F (WT/70)^1
    allo_vp <- fixed(0.75); label("Allometric exponent on Vp/F (unitless, fixed at 3/4; paper-as-printed)") # Table 3 footnote paragraph mark mark: Vp/F (WT/70)^(3/4)
    allo_q  <- fixed(1.00); label("Allometric exponent on Q/F (unitless, fixed at 1; paper-as-printed)")    # Table 3 footnote paragraph mark mark: Q/F (WT/70)^1

    # Inter-individual variability (log-normal). The paper reports CV%; the
    # log-scale variance is omega^2 = log(1 + CV^2). etalcl is shared between
    # CL/F (full magnitude) and Vc/F (scaled by theta_scale_cl_vc = 1.56);
    # etalvp is shared between Vp/F (full magnitude) and Q/F (scaled by
    # theta_scale_vp_q = 1.00 fixed). Correlation 42.7% between etalcl and etalvp.
    #   IIV CL/F = 14.9% CV -> omega^2 = log(1 + 0.149^2) = 0.02196
    #   IIV Vp/F = 52.7% CV -> omega^2 = log(1 + 0.527^2) = 0.24505
    #   cov(etalcl, etalvp) = 0.427 * sqrt(0.02196 * 0.24505) = 0.03133
    #   IIV tlag  = 58.5% CV -> omega^2 = log(1 + 0.585^2) = 0.29442
    etalcl + etalvp ~ c(0.02196,
                        0.03133, 0.24505)                                                              # Table 3 final model: IIV CL/F 14.9% CV (RSE 7.10%); IIV Vp/F 52.7% CV (RSE 8.57%); correlation 42.7% (RSE 13.7%)
    etaltlag ~ 0.29442                                                                                 # Table 3 final model: IIV tlag 58.5% CV (RSE 7.83%)

    # IIV on the proportional residual SD itself: each individual has their own
    # residual-error magnitude, sampled log-normally with 33.3% CV. Encoded with
    # the canonical anchor idiom (Muller_2010_clindamycin, Dogterom_2018_asenapine):
    # a fixed log-anchor lrv pairs the eta with a typical-value fixed effect, and
    # etalrv carries the reported variance.
    #   omega^2 = log(1 + 0.333^2) = 0.10516
    lrv    <- fixed(log(1)); label("Residual-variability scaling anchor (fixed log(1))")               # structural anchor; pairs etalrv with a typical-value fixed effect
    etalrv ~ 0.10516                                                                                   # Table 3 final model: 'IIV on Residual unexplained variability' 33.3% CV (RSE 6.89%)

    # Residual error -- proportional in linear space (Niebecker's "additive on
    # log-scale" formulation is equivalent to proportional in linear space).
    # Two magnitudes: a base phase-1 SD and an incremental phase-3 SD that
    # combine in quadrature for phase-3 observations.
    propSdPhase1    <- 0.142; label("Proportional residual SD in phase 1 healthy volunteers (fraction)")              # Table 3 final model: phase 1 RUV = 14.2% CV (RSE 2.80%)
    propSdPhase3Inc <- 0.544; label("Incremental proportional residual SD added in phase 3 patients (fraction)")      # Table 3 final model: incremental phase 3 RUV = 54.4% CV (RSE 1.80%)
  })

  model({
    # Allometric body-weight scaling. Reference 70 kg per Table 3 footnote.
    wt_ratio <- WT / 70

    # Phase-1-only P-gp-inhibitor effect indicator
    # (P-gp effects on CL/F and F are applied only when STUDY_HOKVTE = 0).
    pgp_phase1 <- PGP_INH * (1 - STUDY_HOKVTE)

    # Renal-CL slope above CLcr = 90 mL/min, scaled by (1 + theta_Scale) in phase 3 patients.
    # In phase 1: slope2_eff = 0.0321 L/h/(mL/min);
    # in phase 3: slope2_eff = 0.0321 * (1 + 2.74) = 0.120 L/h/(mL/min).
    slope2_eff <- e_crcl_cl_renal_slope2 * (1 + e_study_hokvte_cl_renal_slope2 * STUDY_HOKVTE)

    # Truncate CLcr at 150 mL/min per Niebecker 2015 Methods, base model development.
    # Implemented via indicator multiplication (avoids rxode2 parser ambiguity for min/max).
    crcl_eff <- CRCL * (CRCL <= 150) + 150 * (CRCL > 150)

    # Piecewise-linear renal CL component (apparent units L/h, before WT and IIV scaling).
    # For CLcr <= 90: cl_renal_typ = slope1 * CLcr
    # For CLcr >  90: cl_renal_typ = slope1 * 90 + slope2_eff * (CLcr - 90)
    # Implemented via indicator multiplication to keep the parse explicit.
    crcl_below_90 <- crcl_eff * (crcl_eff <= 90) + 90 * (crcl_eff > 90)
    crcl_above_90 <- (crcl_eff - 90) * (crcl_eff > 90)
    cl_renal_typ <- e_crcl_cl_renal_slope1 * crcl_below_90 +
                    slope2_eff * crcl_above_90

    # Non-renal CL component with fed-state effect (FED = 1 in study 6 only).
    cl_nonren_typ <- exp(lcl_nonren) * (1 + e_fed_cl_nonren * FED)

    # Total apparent CL/F = (CLnr + CLr) * (WT/70)^(3/4) * exp(etalcl) * P-gp factor.
    # The P-gp +33.4% factor is applied to phase 1 subjects only via pgp_phase1.
    cl <- (cl_nonren_typ + cl_renal_typ) *
          wt_ratio^allo_cl *
          (1 + e_pgp_inh_cl * pgp_phase1) *
          exp(etalcl)

    # Apparent central volume Vc/F = theta_Vc * (WT/70)^1 * exp(etalcl * 1.56) * Asian factor.
    # etalcl is shared with CL/F at full magnitude; for Vc/F it is scaled by theta_scale_cl_vc = 1.56.
    vc <- exp(lvc + etalcl * theta_scale_cl_vc) *
          wt_ratio^allo_vc *
          (1 + e_race_asian_vc * RACE_ASIAN)

    # Apparent peripheral volume Vp/F = theta_Vp * (WT/70)^(3/4) * exp(etalvp).
    # Paper-as-printed: allometric exponent 3/4 on Vp/F.
    vp <- exp(lvp + etalvp) * wt_ratio^allo_vp

    # Apparent intercompartmental clearance Q/F = theta_Q * (WT/70)^1 * exp(etalvp * 1.00) * phase 3 factor.
    # etalvp is shared with Vp/F at full magnitude; for Q/F it is scaled by theta_scale_vp_q = 1.00 (FIXED).
    # Paper-as-printed: allometric exponent 1 on Q/F.
    q <- exp(lq + etalvp * theta_scale_vp_q) *
         wt_ratio^allo_q *
         (1 + e_study_hokvte_q * STUDY_HOKVTE)

    # Absorption rate with fed-state effect (FED = 1 in study 6 only).
    ka <- exp(lka) * (1 + e_fed_ka * FED)

    # Absorption lag time (fixed structural parameter, with IIV).
    tlag <- exp(ltlag + etaltlag)

    # Micro-constants for the two-compartment ODE.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order oral absorption from depot to central, with lag.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1

    # Absorption lag (alag on depot).
    alag(depot) <- tlag

    # Bioavailability with phase-1-only P-gp inhibitor effect (+125%).
    f(depot) <- exp(lfdepot) * (1 + e_pgp_inh_f * pgp_phase1)

    # Plasma concentration. Dose in mg, volumes in L, so central/vc is in mg/L = ug/mL;
    # multiply by 1000 to convert to ng/mL (bioanalytical units of the source assays).
    Cc <- (central / vc) * 1000

    # Residual error: combined proportional SD with phase-3 increment in quadrature
    # and individual log-normal scaling via etalrv (33.3% CV on the residual SD).
    # Anchor lrv is fixed at log(1) so exp(lrv + etalrv) = exp(etalrv) when etalrv = 0.
    propSdEff <- sqrt(propSdPhase1^2 + (STUDY_HOKVTE * propSdPhase3Inc)^2) *
                 exp(lrv + etalrv)
    Cc ~ prop(propSdEff)
  })
}
