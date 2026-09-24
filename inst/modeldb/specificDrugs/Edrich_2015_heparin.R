Edrich_2015_heparin <- function() {
  description <- paste(
    "One-compartment population PK + linear direct-effect PD model for",
    "intravenous heparin in adults undergoing catheter-based ablation of",
    "atrial fibrillation (Edrich 2015). Because only heparin doses and the",
    "resulting activated clotting times (ACT) were recorded, the PK and PD",
    "layers are not separately identifiable: the central compartment holds",
    "heparin units in the estimated blood volume and a single multiplicative",
    "sensitivity coefficient k_ACT maps that scaled concentration onto the",
    "ACT in seconds, ACT = ACT_BASE + k_ACT * Cc. The central volume is not",
    "estimated but set to the weight-based estimated blood volume, which is",
    "sex-dependent. k_ACT carries a four-level multiplicative factor for the",
    "patient's chronic oral anticoagulant at presentation (none = reference;",
    "warfarin with INR < 2; warfarin with INR >= 2; dabigatran stopped about",
    "27 h earlier), which is the paper's finding: warfarin patients are about",
    "twice as heparin-sensitive as dabigatran or unanticoagulated patients.",
    "Clearance carried no group effect in the final model. The point-of-care",
    "instrument ceiling of 400 s is NOT applied to the prediction here (it is",
    "an estimation device in the source control stream, not pharmacology);",
    "see the vignette for its consequences.",
    sep = " "
  )
  reference <- paste(
    "Edrich T, Frendl G, Michaud G, Paschalidis ICh.",
    "Heparin requirements for full anticoagulation are higher for patients on",
    "dabigatran than for those on warfarin - a model-based study.",
    "Clin Pharmacol Adv Appl. 2015;7:19-25. doi:10.2147/CPAA.S72185.",
    sep = " "
  )
  vignette <- "Edrich_2015_heparin"
  units <- list(time = "min", dosing = "IU", concentration = "IU/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The source's own comment on the $ERROR block names the
  # matrix explicitly -- F is a 'scaled drug concentration' in 'Units
  # heparin/L blood volume' -- and V1 is set to the estimated BLOOD volume,
  # so whole blood (not plasma) is the correct specimen here.
  compartmentData <- list(
    central = list(analyte = "heparin", units = "IU", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the central volume of distribution directly: the source fixes",
        "V1 to the weight-based estimated blood volume rather than estimating",
        "it (supplement $PK block). Linear in WT (exponent 1, structural --",
        "not an estimated allometric exponent). Edrich 2015 does NOT tabulate",
        "body weight anywhere; Table 1 reports BMI only. See the vignette for",
        "the back-solve that recovers a cohort-typical weight of about 77 kg",
        "from the paper's own printed half-life range and group clearances."
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the per-kg blood-volume constant used for the central",
        "volume. The source codes sex as M1F2 (1 = male, 2 = female) and",
        "writes V1 = (0.075 - 0.005 * (M1F2 - 1)) * WT, so M1F2 - 1 is",
        "exactly SEXF and no value transformation beyond the offset is",
        "needed. Sex does NOT act on k_ACT in the final model: the paper's",
        "male-vs-female k_ACT difference in group D (0.14 vs 0.12) is a",
        "post-hoc subgroup comparison (Table 3, Results), not a term in the",
        "GROUPMOD expression."
      ),
      source_name = "M1F2"
    ),
    ACT_BASE = list(
      description = "Pre-heparin baseline activated clotting time",
      units = "s",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Read per subject from the data set and added to every prediction",
        "(source $ERROR block: Predicted_ACT = F * KACT + BASEACT). It is a",
        "covariate, not an estimated parameter -- the source fits no",
        "population baseline. Baseline ACT differs by group and is one of the",
        "paper's findings: medians 144 s (no anticoagulant), 155 s",
        "(dabigatran), 169 s (warfarin INR < 2), 182 s (warfarin INR >= 2)",
        "per Table 1. Measured on a Hemochron Signature Elite whole-blood",
        "microcoagulation system with a Hemochron Jr cartridge, range",
        "0-400 s. Pre-heparin ACT was missing for 21-33% of patients; the",
        "source imputed the group median in those cases (Methods)."
      ),
      source_name = "BASEACT"
    ),
    CONMED_WARFARIN = list(
      description = "Chronic warfarin at presentation (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chronic oral anticoagulant)",
      notes = paste(
        "Time-fixed at presentation. Combined with INR_BASE to reproduce the",
        "source's two warfarin indicators: COUMLO is warfarin with an",
        "unintentionally low INR < 2 and COUMHI is warfarin with INR >= 2",
        "within the last 3 days (Methods). Mutually exclusive with",
        "CONMED_DABIGATRAN; both zero denotes the source's NEITHER reference",
        "group. Warfarin patients took their evening dose a median 15 h",
        "before the procedure."
      ),
      source_name = "COUMLO / COUMHI"
    ),
    CONMED_DABIGATRAN = list(
      description = "Chronic dabigatran at presentation (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chronic oral anticoagulant)",
      notes = paste(
        "Time-fixed at presentation. The source's PRADAX indicator. Mutually",
        "exclusive with CONMED_WARFARIN. The last dabigatran dose was taken a",
        "median 27 h (IQR 24-31) before the procedure, roughly two",
        "elimination half-lives, so this indicator denotes recent-and-stopped",
        "rather than ongoing dabigatran exposure -- the paper's conclusion is",
        "that at that withholding interval heparin sensitivity is",
        "indistinguishable from no anticoagulation at all."
      ),
      source_name = "PRADAX"
    ),
    INR_BASE = list(
      description = "Last pre-procedural international normalized ratio",
      units = "(unitless ratio; INR has no units)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only to split CONMED_WARFARIN into the source's COUMLO",
        "(INR < 2) and COUMHI (INR >= 2) strata at the 2.0 cut point given in",
        "Methods; it carries no continuous effect of its own. Group means",
        "1.2 (dabigatran), 1.8 (warfarin INR < 2), 2.3 (warfarin INR >= 2),",
        "1.0 (no anticoagulant) per Results and Table 1."
      ),
      source_name = "INR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 188L,
    n_studies = 1L,
    age_median = "64 y (dabigatran), 58 y (warfarin INR < 2), 62 y (warfarin INR >= 2), 58 y (no anticoagulant); Table 1 medians",
    weight_range = "Not reported. Table 1 gives BMI only (medians 29.3, 28.7, 29.3, 26.3 kg/m^2 by group); see covariateData[['WT']]$notes",
    bmi_range = "Group medians 26.3-29.3 kg/m^2 (Table 1)",
    sex_female_pct = 31.5,
    race_ethnicity = "Caucasian 94% (dabigatran), 95% (warfarin INR < 2), 94% (warfarin INR >= 2), 81% (no anticoagulant) per Table 1; about 92% of the pooled cohort. No other race categories reported, and the paper found no significant difference in any model parameter by race.",
    disease_state = "Adults with atrial fibrillation presenting for catheter-based atrial ablation requiring full intraprocedural anticoagulation with heparin.",
    dose_range = "Not reported. Heparin dosing was not standardized in this retrospective cohort; each patient received a median of three intravenous heparin boluses and six ACT measurements per case, with no heparin infusions (Results). The administered amounts are not tabulated -- Table 2 reports only the ACT response normalized per unit of heparin per kg.",
    regions = "Single-centre: Brigham and Women's Hospital, Boston, MA, USA (IRB-approved retrospective chart review, January 2011 to June 2012)",
    n_group_dabigatran = 66L,
    n_group_warfarin_inr_low = 42L,
    n_group_warfarin_inr_high = 53L,
    n_group_no_anticoagulant = 27L,
    notes = paste(
      "Baseline demographics from Edrich 2015 Table 1. Group W (warfarin,",
      "n = 95) was subdivided by the last pre-procedural INR into W_low",
      "(INR < 2, n = 42) and W_high (INR >= 2, n = 53). Estimation: NONMEM 7",
      "with PLTTools, FOCE with interaction, ADVAN1. A two-compartment",
      "alternative was rejected on objective function (10,424 vs 10,372).",
      "Adding the group factor to k_ACT improved the objective function from",
      "10,574 to 10,372 and the RMS error from 68.6 to 51.1 ACT-seconds",
      "(20.5% to 15.3% of the average ACT). All ACT values were measured on",
      "a Hemochron Signature Elite whole-blood system with a 0-400 s range;",
      "the first post-heparin ACT exceeded 400 s in 69%, 77%, 24% and 7% of",
      "the W_low, W_high, D and N groups respectively, so the group-W entries",
      "of Table 2 are right-censored and understate the true response (the",
      "paper says so explicitly in the Discussion)."
    )
  )

  ini({
    # =================================================================
    # PHARMACOKINETICS
    #
    # The source's final THETA / OMEGA / SIGMA estimates are never
    # tabulated; the supplement prints only the NONMEM starting values
    # in 'LOWER STARTING UPPER' form. Every value below is therefore
    # recovered from the paper's reported RESULTS, and each recovered
    # value is checked against the corresponding $THETA bound as an
    # independent consistency test (all pass; see the trailing
    # comments and the vignette source-trace table).
    #
    # Clearance: the final model carries no group effect on CL ('the
    # clearances ... did not differ significantly among groups'), so one
    # pooled typical value is used -- the n-weighted mean of the four
    # group medians of individual CL reported in Results:
    #   (22.4*42 + 22.1*53 + 26.5*66 + 23.4*27) / 188 = 23.9 mL/min
    # = 0.0239 L/min. For a log-normally distributed individual CL the
    # median equals the typical value, so a weighted mean of group
    # medians estimates THETA(1).
    # =================================================================
    lcl <- log(23.9 / 1000)
    label("Heparin clearance from the blood volume (L/min)")                                  # Results: group medians of individual CL 22.4 / 22.1 / 26.5 / 23.4 mL/min; n-weighted mean 23.9 mL/min. Inside the supplement bound (0, 10) L/min for THETA(1)

    # =================================================================
    # Central volume = estimated BLOOD volume, not an estimated
    # parameter. Supplement $PK:
    #   V1 = (0.075 - 0.005*(M1F2-1)) * WT * (1 + THETA(6))   [Liters]
    # with M1F2 = 1 for male and 2 for female, so the numeric literal
    # is 0.075 L/kg for men and 0.070 L/kg for women. The comment on
    # that same code line and the Methods prose both instead say
    # 70 mL/kg for men and 65 mL/kg for women. The code literal is
    # adopted here per the standing 'trust the printed equation' rule,
    # because it is what actually ran and k_ACT was estimated
    # conditional on it; the discrepancy (a 7.1% scale factor on every
    # predicted ACT increment) is recorded in the vignette Errata.
    #
    # Sex-stratified because the paper reports the same quantity once
    # per sex rather than a reference value plus an offset; both strata
    # therefore carry an explicit suffix.
    #
    # THETA(6), an estimated multiplicative adjustment to this
    # weight-based blood volume (bounds -0.5 to 0.5, start 0.1), is
    # never reported. It is taken as 0 here, i.e. V1 is exactly the
    # per-kg blood volume times weight. Recorded in vignette Errata.
    # =================================================================
    lvc_male <- fixed(log(0.075))
    label("Blood volume per kg body weight, males (L/kg)")                                    # Supplement $PK: V1 = (0.075 - 0.005*(M1F2-1))*WT; M1F2 = 1 for male
    lvc_female <- fixed(log(0.070))
    label("Blood volume per kg body weight, females (L/kg)")                                  # Supplement $PK: same expression with M1F2 = 2 for female

    # =================================================================
    # PHARMACODYNAMICS -- linear direct effect
    #
    # Supplement $ERROR: Predicted_ACT = F*KACT + BASEACT, where F is
    # the scaled concentration in heparin units per litre of blood
    # volume. KACT = THETA(2) * GROUPMOD and the control stream's own
    # comment states that group 'neither' is the reference with
    # multiplier 1, so THETA(2) is the no-anticoagulant value.
    #
    # Table 3 reports medians of INDIVIDUAL k_ACT by group; those
    # medians estimate the group typical values (median = typical for a
    # log-normal eta, and KACT carries EXP(ETA(2)) in $PK).
    # =================================================================
    lslope <- log(0.12)
    label("Sensitivity of ACT to scaled heparin concentration, no-anticoagulant reference group (s*L/IU)")      # Table 3: median individual k_ACT 0.12 s*L/units in group n (the GROUPMOD reference). Inside the supplement bound (0, 20) for THETA(2)

    # =================================================================
    # Group multipliers on k_ACT, recovered as the ratio of each
    # group's Table 3 median to the reference group's. Written as the
    # explicit ratio so the provenance is readable. Each recovered
    # value falls inside the corresponding supplement $THETA bound,
    # which is an independent confirmation that group n is the
    # reference -- THETA(3) is bounded above by 2 and 0.23/0.12 = 1.92
    # only just fits.
    # =================================================================
    e_warflo_slope <- 0.23 / 0.12
    label("Multiplier on k_ACT for chronic warfarin with INR < 2 (unitless)")                 # Table 3: median k_ACT 0.23 (group W low) / 0.12 (group n) = 1.92. Inside the supplement bound (1, 2) for THETA(3)
    e_warfhi_slope <- 0.30 / 0.12
    label("Multiplier on k_ACT for chronic warfarin with INR >= 2 (unitless)")                # Table 3: median k_ACT 0.30 (group W high) / 0.12 (group n) = 2.50. Inside the supplement bound (1, 4) for THETA(4)
    e_dabigatran_slope <- 0.13 / 0.12
    label("Multiplier on k_ACT for chronic dabigatran stopped about 27 h earlier (unitless)")  # Table 3: median k_ACT 0.13 (group D) / 0.12 (group n) = 1.08. Inside the supplement bound (0.6, 1.5) for THETA(5)

    # =================================================================
    # RESIDUAL ERROR
    #
    # The source declares both components -- $SIGMA(1) proportional and
    # $SIGMA(2) additive, with Y = IPRED*(1 + EPS(1)) + EPS(2) -- but
    # reports only the final model's combined root-mean-square error:
    # 51.1 ACT-seconds, equal to 15.3% of the average ACT (Results).
    # Those are one quantity written two ways (they agree at an average
    # ACT of 51.1/0.153 = 334 s), so they cannot be split between the
    # two components. The whole RMS is assigned to the proportional
    # term and the unrecoverable additive term is held at zero rather
    # than invented. Recorded in vignette Errata.
    # =================================================================
    propSd <- 0.153
    label("Proportional residual SD on ACT (fraction)")                                       # Results: final-model RMS error 51.1 ACT-seconds = 15.3% of the average ACT
    addSd <- fixed(0)
    label("Additive residual SD on ACT (s; 0 -- the split between the two declared residual components is not reported)")   # Supplement $SIGMA(2) declares an additive component; its final value is not reported

    # =================================================================
    # INTER-INDIVIDUAL VARIABILITY
    #
    # The source fits $OMEGA BLOCK(2) with ETA(1) on CL and ETA(2) on
    # KACT (per $PK and the Results sentence naming 'EXP(ETA(2))' as
    # the k_ACT inter-individual variability; the $OMEGA block's own
    # comment mis-labels ETA2 as applying to V1, which has no eta at
    # all). Only the starting values (0.3; 0.1 0.3) are printed, never
    # the final estimates, so this is a typical-value model with no
    # etas. The Table 3 and Results interquartile ranges of individual
    # CL and k_ACT are post-hoc empirical-Bayes spreads, shrunken
    # toward the typical value, and are NOT valid OMEGA estimates.
    # Recorded in vignette Errata.
    # =================================================================
  })

  model({
    # --- Group indicators, reconstructing the source's GROUPMOD -------
    # Supplement $PK:
    #   GROUPMOD = NEITHER*1 + COUMLO*THETA(3)
    #              + COUMHI*THETA(4) + PRADAX*THETA(5)
    # The four indicators are mutually exclusive and exhaustive, so the
    # reference group is recovered as 1 - warfarin - dabigatran. The
    # warfarin arm is split at the INR = 2.0 cut point given in Methods.
    warflo <- CONMED_WARFARIN * (INR_BASE < 2)
    warfhi <- CONMED_WARFARIN * (INR_BASE >= 2)
    neither <- 1 - CONMED_WARFARIN - CONMED_DABIGATRAN

    groupmod <-
      neither * 1 +
      warflo * e_warflo_slope +
      warfhi * e_warfhi_slope +
      CONMED_DABIGATRAN * e_dabigatran_slope

    # --- Individual parameters ---------------------------------------
    # vc is the estimated blood volume: a per-kg constant (L/kg) chosen
    # by sex, times body weight. Linear in WT by construction, not via
    # an estimated allometric exponent.
    cl <- exp(lcl)
    vc <- (exp(lvc_male) * (1 - SEXF) + exp(lvc_female) * SEXF) * WT
    slope <- exp(lslope) * groupmod

    kel <- cl / vc

    # --- ODE: single central compartment, IV bolus dosing -------------
    d/dt(central) <- -kel * central

    # --- Outputs ------------------------------------------------------
    # Doses are in heparin units and vc is in L, so central / vc is in
    # units per litre of blood volume -- the source's F. Multiplying by
    # k_ACT (s*L/units) gives an ACT increment in seconds, which adds
    # to the subject's own pre-heparin baseline. Heparin concentration
    # itself was never assayed in this study, so Cc is a structural
    # intermediate and carries no residual error; ACT is the only
    # observed endpoint.
    #
    # The source's IF (IPRED.GT.400) IPRED = 400 is deliberately NOT
    # reproduced: its stated purpose is 'do not penalize if model
    # guesses over 400', i.e. it is a likelihood device for the
    # instrument's measurement ceiling, not part of the pharmacology.
    # Applying it here would silently censor every simulated ACT.
    Cc <- central / vc
    ACT <- ACT_BASE + slope * Cc

    ACT ~ add(addSd) + prop(propSd)
  })
}
