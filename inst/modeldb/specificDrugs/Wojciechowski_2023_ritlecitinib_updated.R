Wojciechowski_2023_ritlecitinib_updated <- function() {
  description <- "Updated model (iteration 2 of 3) for oral ritlecitinib: two-compartment disposition with first-order absorption and a direct-response non-stationary (autoinhibitory) Imax effect of the peripheral-compartment concentration on both apparent clearance and bioavailability. Base-model parameters were re-estimated on the pooled healthy participant, rheumatoid arthritis, ulcerative colitis, alopecia areata, vitiligo and moderate hepatic impairment data (668 individuals, 5187 concentrations) and a stepwise covariate analysis was run. Carries rheumatoid arthritis, ulcerative colitis, alopecia areata and vitiligo effects on CL/F; ulcerative colitis and moderate hepatic impairment effects on F; high-fat-meal, 800 mg dose and capsule effects on ka; an over-encapsulated-capsule effect on loss from the depot; and inflammatory-disease scaling of the IIV and proportional residual error magnitudes. Allometric weight scaling with fixed exponents 0.75 on clearances and 1 on volumes referenced to 70 kg."
  reference <- paste(
    "Wojciechowski J, Purohit VS, Huh Y, Banfield C, Nicholas T.",
    "Evolution of Ritlecitinib Population Pharmacokinetic Models During",
    "Clinical Drug Development. Clin Pharmacokinet. 2023;62(12):1765-1779.",
    "doi:10.1007/s40262-023-01318-3. PMCID PMC10684409.",
    "Updated-model column of Table 2; univariable covariate screen in",
    "Table S3 and the complete NONMEM control stream in the Electronic",
    "Supplementary Material.",
    sep = " "
  )
  vignette <- "Wojciechowski_2023_ritlecitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "ritlecitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ritlecitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ritlecitinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight; allometrically scales CL/F, Q/F, Vc/F and Vp/F",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight (NONMEM column BWT). Reference 70 kg with fixed exponents 0.75 on CL/F and Q/F and 1.00 on Vc/F and Vp/F; Wojciechowski 2023 Table 2 footnote. Updated-model analysis population median 75.0 kg (range 35.1-164). Table S3 run 5 tested an estimated weight exponent on Vc/F but the fixed allometric form was retained.",
      source_name        = "BWT"
    ),
    DIS_RA = list(
      description        = "Rheumatoid arthritis patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant, moderate hepatic impairment participant, or a patient with one of the other pooled indications)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 1). Carries a multiplicative effect on CL/F and contributes to the inflammatory-disease group (PTST in {1, 2, 3, 5}) that scales the IIV and residual-error magnitudes. The CL/F reference category is healthy participants.",
      source_name        = "PTST"
    ),
    DIS_UC = list(
      description        = "Ulcerative colitis patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant, moderate hepatic impairment participant, or a patient with one of the other pooled indications)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 2). Carries multiplicative effects on both CL/F and F, and contributes to the inflammatory-disease group that scales the IIV and residual-error magnitudes. The F effect's reference group is healthy participants plus RA, alopecia areata and vitiligo patients (ESM $PK: reference subjects for COVPTSTF are healthy participants, RA, AA, vitiligo). Moderate-to-severe UC per the phase IIb study.",
      source_name        = "PTST"
    ),
    DIS_ALOPECIA_AREATA = list(
      description        = "Alopecia areata patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant, moderate hepatic impairment participant, or a patient with one of the other pooled indications)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 3). Carries a multiplicative effect on CL/F and contributes to the inflammatory-disease group that scales the IIV and residual-error magnitudes.",
      source_name        = "PTST"
    ),
    DIS_VITILIGO = list(
      description        = "Vitiligo patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant, moderate hepatic impairment participant, or a patient with one of the other pooled indications)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 5). Carries a multiplicative effect on CL/F and contributes to the inflammatory-disease group that scales the IIV and residual-error magnitudes. Non-segmental vitiligo per the phase IIb study.",
      source_name        = "PTST"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; in this analysis population every non-flagged participant had normal hepatic function)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 6). Classified by Child-Pugh score (Wojciechowski 2023 Sect. 2.4.2), from the dedicated phase I hepatic impairment study NCT04016077 (n = 10). Carries a multiplicative effect on F only. Note that moderate hepatic impairment participants are EXCLUDED from the inflammatory-disease group that scales the IIV and residual-error magnitudes (the control stream tests PTST > 0 AND PTST < 6).",
      source_name        = "PTST"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted, or food intake not controlled)",
      notes              = "Per-dose-record indicator derived from the NONMEM FOOD column (0 = fasted, 1 = high-fat meal, 2 = not controlled for food); the control stream applies the effect only when FOOD = 1, so both fasted and not-controlled records are the reference.",
      source_name        = "FOOD"
    ),
    DOSE = list(
      description        = "Administered ritlecitinib dose level for the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only as a step-function switch: the control stream applies a separate ka effect when DOSE = 800 exactly (IF (DOSE.EQ.800)). Every dose level below 800 mg is the reference. Updated-model dose range 5-800 mg.",
      source_name        = "DOSE"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (any capsule)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet or oral solution)",
      notes              = "Per-dose-record indicator derived from the NONMEM FORM column (2 = tablet, 11 = capsule); the control stream comment names the reference as tablet and solution. Carries a multiplicative effect on ka only. The capsule arm is the pilot capsule formulation evaluated in the phase I study NCT04004663; Wojciechowski 2023 Discussion reports it gave a 14% lower Cmax and approximately no change in AUC(tau) relative to the tablet.",
      source_name        = "FORM"
    ),
    FORM_RIT_OVERENCAP = list(
      description        = "Over-encapsulated ritlecitinib capsule indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet, oral solution, pilot capsule, or large-API-particle-size capsule)",
      notes              = "Per-dose-record indicator for the over-encapsulated capsule arm only. The control stream selects it with IF (FORM.EQ.11 .AND. TRT.EQ.4), i.e. a capsule record within treatment arm 4; TRT is not listed in the deposited $INPUT record, so the deposited stream appears to be a lightly edited copy of the executed one and the condition is represented here by a single dedicated indicator column. Scales the rate of loss from the depot compartment, NOT the amount reaching the systemic circulation: Table S3 tested both parameterisations (run 17 loss from depot, dOFV -40.3; run 16 amount reaching systemic circulation, dOFV -39.4) and the loss-from-depot form was carried forward.",
      source_name        = "FORM (with TRT)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 668L,
    n_studies      = 10L,
    n_observations = 5187L,
    age_range      = "18.0-74.0 years",
    age_median     = "43.0 years",
    weight_range   = "35.1-164 kg",
    weight_median  = "75.0 kg",
    sex_female_pct = 45.2,
    race_ethnicity = c(White = 79.5, Asian = 12.0, Black = 2.4, Other = 1.3, Missing = 4.8),
    disease_state  = "Healthy participants (98; 14.7%), rheumatoid arthritis (42; 6.3%), ulcerative colitis (150; 22.5%), alopecia areata (70; 10.5%), vitiligo (298; 44.6%) and moderate hepatic impairment (10; 1.5%).",
    dose_range     = "5-800 mg/day oral ritlecitinib (single and multiple dose; tablet, capsule and oral-solution formulations, fasted and high-fat-meal conditions)",
    regions        = "Global (12-trial ritlecitinib clinical development programme; see Table S1 of the Electronic Supplementary Material)",
    renal_function = "Normal (severe renal impairment participants entered only at the final-model iteration)",
    hepatic_function = "Normal, plus 10 participants with Child-Pugh moderate hepatic impairment (NCT04016077)",
    notes          = "Demographics from Wojciechowski 2023 Table 1, updated-model column. Median baseline creatinine clearance 113 mL/min (range 45.4-271) and median albumin 4.50 g/dL (range 2.80-5.40). Below-limit-of-quantification observations were excluded during estimation. Age, sex, race, baseline creatinine clearance and moderate hepatic impairment on CL/F were screened (Table S3) but not retained in the final covariate model."
  )

  # Implementation notes (see the vignette 'Assumptions and deviations'
  # section for the full justification of each item):
  # * Structure. Identical to the base model
  #   (Wojciechowski_2023_ritlecitinib_base); Sect. 2.4 states the base
  #   model was re-estimated on the pooled dataset to provide a reference
  #   structural model for the stepwise covariate analysis, and Sect. 4
  #   confirms the structural elements were retained across iterations.
  #   The deposited $DES block for this iteration is
  #     AUTOINH = 1 + IMAXP*PERIC/(IC50P + PERIC)
  #     INHF    = 1 + (1 - AUTOINH)
  #     KAD     = KA*COVFOODKA*COVDOSEKA*COVFORMKA
  #     FD      = COVFORMF   ; fraction into the depot
  #     FA      = COVPTSTF   ; fraction absorbed
  #     DADT(1) = -KAD*A(1)*FD
  #     DADT(2) = KAD*A(1)*INHF*FA - Q*CENTC + Q*PERIC - CL*AUTOINH*CENTC
  #     DADT(3) = Q*CENTC - Q*PERIC
  # * CL/F point estimate. Table 2 prints 113 (95% CI 105-125) for the
  #   updated model, but 105-125 is symmetric about 115, not 113, and
  #   every other estimate in Table 2 sits exactly at its own interval
  #   midpoint (the footnote states the intervals are asymptotic). The
  #   final model's frequentist-prior block in the Electronic
  #   Supplementary Material - which by construction carries the updated
  #   model's final estimates - lists "115 FIX ; PTVCL". Both lines of
  #   evidence give 115, so the printed 113 is a typographical error and
  #   115 is used here. Every other $THETAP prior value matches its
  #   Table 2 updated-model entry exactly. See the vignette Errata.
  # * Over-encapsulated capsule. The FD term multiplies the depot LOSS
  #   RATE and is deliberately not applied to the amount entering the
  #   central compartment: Table S3 evaluated both forms and selected
  #   "effect of over-encapsulated capsule on loss from depot" (run 17)
  #   over "effect on amount reaching systemic circulation" (run 16).
  #   Encoded exactly as deposited. Because the central-compartment input
  #   is not scaled by FD, total absorbed drug is amplified by 1/FD when
  #   FORM_RIT_OVERENCAP = 1; this is a property of the published
  #   parameterisation, not of the transcription. The term is inert at
  #   the reference formulation.
  # * Inlined state expressions, units, IIV scale, residual error and the
  #   inflammatory-disease SD scaling all follow the base model; see the
  #   implementation notes in Wojciechowski_2023_ritlecitinib_base.R.
  # * Inflammatory-disease group. PTST in {1, 2, 3, 5} = RA, UC, alopecia
  #   areata and vitiligo. Healthy participants (PTST = 0) and moderate
  #   hepatic impairment (PTST = 6) are the reference.
  ini({
    # ----- Structural parameters (Wojciechowski 2023 Table 2, updated-model column) -----
    lcl   <- log(115);   label("Apparent clearance CL/F at 70 kg, uninhibited (L/h)")                          # ESM final-model $THETAP: "115 FIX ; PTVCL"; Table 2 updated prints 113 (95% CI 105-125) whose midpoint is 115 - see implementation notes and vignette Errata
    lvc   <- log(149);   label("Apparent central volume Vc/F at 70 kg (L)")                                    # Table 2 updated: 149 (95% CI 143-155); ESM $THETAP PTVV2 149
    lq    <- log(0.304); label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")                    # Table 2 updated: 0.304 (95% CI 0.264-0.344); ESM $THETAP PTVQ 0.304
    lvp   <- log(4.67);  label("Apparent peripheral volume Vp/F at 70 kg (L)")                                 # Table 2 updated: 4.67 (95% CI 4.29-5.05); ESM $THETAP PTVV3 4.67
    lka   <- log(8.51);  label("First-order absorption rate constant ka (1/h)")                                # Table 2 updated: 8.51 (95% CI 6.71-10.3); ESM $THETAP PTVKA 8.51

    # ----- Non-stationary (autoinhibition) parameters -----
    imax  <- -0.488;     label("Maximum non-stationary fractional effect Imax,P on CL/F and F (unitless)")     # Table 2 updated: -0.488 (95% CI -0.526 to -0.450); ESM $THETAP PTVIMAXP -0.488
    lic50 <- log(15.1);  label("Peripheral concentration giving half the maximum non-stationary effect IC50,P (ng/mL)")  # Table 2 updated: 15.1 (95% CI 11.6-18.6); ESM $THETAP PTVIC50P 15.1

    # ----- Allometric exponents (fixed, Table 2 footnote) -----
    e_wt_cl_q  <- fixed(0.75); label("Allometric WT exponent shared by CL/F and Q/F (unitless)")               # Table 2 footnote: exponent 0.75 on clearance parameters, reference 70 kg
    e_wt_vc_vp <- fixed(1.00); label("Allometric WT exponent shared by Vc/F and Vp/F (unitless)")              # Table 2 footnote: exponent 1 on volume parameters, reference 70 kg

    # ----- Categorical covariate effects on CL/F (reference = healthy participants) -----
    e_ra_cl              <- -0.496; label("Fractional effect of rheumatoid arthritis on CL/F (unitless)")      # Table 2 updated: -0.496 (95% CI -0.587 to -0.405); ESM $PK THETA(12) PTSTCL1
    e_uc_cl              <- -0.560; label("Fractional effect of ulcerative colitis on CL/F (unitless)")        # Table 2 updated: -0.560 (95% CI -0.619 to -0.501); ESM $PK THETA(13) PTSTCL2
    e_alopecia_areata_cl <- -0.322; label("Fractional effect of alopecia areata on CL/F (unitless)")           # Table 2 updated: -0.322 (95% CI -0.448 to -0.196); ESM $PK THETA(14) PTSTCL3; ESM $THETAP PTVPTSTCL3 -0.322
    e_vitiligo_cl        <- -0.214; label("Fractional effect of vitiligo on CL/F (unitless)")                  # Table 2 updated: -0.214 (95% CI -0.293 to -0.135); ESM $PK THETA(15) PTSTCL5

    # ----- Categorical covariate effects on F (reference = healthy participants, RA, AA, vitiligo) -----
    e_uc_fdepot         <- -0.224; label("Fractional effect of ulcerative colitis on bioavailability F (unitless)")          # Table 2 updated: -0.224 (95% CI -0.305 to -0.143; the printed upper bound drops its minus sign); ESM $PK THETA(16) PTSTF2
    e_hepimp_mod_fdepot <-  0.255; label("Fractional effect of moderate hepatic impairment on bioavailability F (unitless)") # Table 2 updated: 0.255 (95% CI 0.135-0.375); ESM $PK THETA(17) PTSTF6

    # ----- Categorical covariate effects on ka -----
    e_fed_highfat_ka   <- -0.750; label("Fractional effect of a high-fat meal on ka (unitless)")               # Table 2 updated: -0.750 (95% CI -0.800 to -0.700); ESM $PK THETA(18) FOODKA1
    e_dose800_ka       <- -0.833; label("Fractional effect of the 800 mg dose level on ka (unitless)")         # Table 2 updated: -0.833 (95% CI -0.876 to -0.790); ESM $PK THETA(19) DOSEKA800
    e_form_capsule_ka  <- -0.598; label("Fractional effect of the capsule formulation on ka (unitless)")       # Table 2 updated: -0.598 (95% CI -0.701 to -0.495); ESM $PK THETA(21) FORMKA11

    # ----- Formulation effect on the rate of loss from the depot -----
    e_form_rit_overencap_depotloss <- -0.134; label("Fractional effect of over-encapsulated capsules on the rate of loss from the depot (unitless)")  # Table 2 updated: -0.134 (95% CI -0.163 to -0.105); ESM $PK THETA(20) FORMF114 entering $DES as FD

    # ----- Inflammatory-disease scaling of the variability magnitudes -----
    # These multiply the SD (eta and residual SD), not the variance.
    e_inflam_iiv_cl <- 1.69;  label("Fractional effect of inflammatory disease on the SD of IIV on CL/F (unitless)")  # Table 2 updated: 1.69 (95% CI 1.18-2.20); ESM $PK THETA(10) VARCLPTST; ESM $THETAP PTVVARCLPTST 1.69
    e_inflam_iiv_vc <- 2.39;  label("Fractional effect of inflammatory disease on the SD of IIV on Vc/F (unitless)")  # Table 2 updated row 2 of the duplicated omega^2 label: 2.39 (95% CI 1.40-3.38); ESM $PK THETA(11) VARV2PTST; ESM $THETAP PTVVARV2PTST 2.39
    e_inflam_propsd <- 0.306; label("Fractional effect of inflammatory disease on the proportional residual SD (unitless)")  # Table 2 updated: 0.306 (95% CI 0.265-0.347); ESM $ERROR THETA(9) RUVPROPTST; ESM $THETAP PTVRUVPROPTST 0.306

    # ----- Between-subject variability (diagonal) -----
    etalcl ~ 0.039204  # Table 2 updated: 19.8 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.198^2; confirmed by ESM $OMEGAP PPPVCL 0.039204
    etalvc ~ 0.015625  # Table 2 updated: 12.5 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.125^2; confirmed by ESM $OMEGAP PPPVV2 0.015625

    # ----- Residual error -----
    propSd <- 0.359; label("Proportional residual SD for the reference (healthy participant) group (fraction)")  # Table 2 updated: 35.9 "% CV"; ESM $ERROR THETA(8) RUVPRO with $SIGMA 1 FIX; ESM $THETAP PTVRUVPRO 0.359
  })
  model({
    ref_wt <- 70  # Table 2 footnote: allometric reference weight

    # ----- Derived covariate terms -----
    # Inflammatory-disease group = PTST in {1, 2, 3, 5}; healthy
    # participants (PTST = 0) and moderate hepatic impairment (PTST = 6)
    # are the reference.
    inflam <- DIS_RA + DIS_UC + DIS_ALOPECIA_AREATA + DIS_VITILIGO

    covptstcl <- 1 + e_ra_cl * DIS_RA + e_uc_cl * DIS_UC +
      e_alopecia_areata_cl * DIS_ALOPECIA_AREATA + e_vitiligo_cl * DIS_VITILIGO
    covptstf  <- 1 + e_uc_fdepot * DIS_UC + e_hepimp_mod_fdepot * HEPIMP_MOD
    covformf  <- 1 + e_form_rit_overencap_depotloss * FORM_RIT_OVERENCAP
    covvarcl  <- 1 + e_inflam_iiv_cl * inflam
    covvarvc  <- 1 + e_inflam_iiv_vc * inflam
    covfoodka <- 1 + e_fed_highfat_ka * FED_HIGHFAT
    covdoseka <- 1 + e_dose800_ka * (DOSE == 800)
    covformka <- 1 + e_form_capsule_ka * FORM_CAPSULE

    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl * covvarcl) * (WT / ref_wt)^e_wt_cl_q * covptstcl
    vc <- exp(lvc + etalvc * covvarvc) * (WT / ref_wt)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / ref_wt)^e_wt_cl_q
    vp <- exp(lvp) * (WT / ref_wt)^e_wt_vc_vp

    kad <- exp(lka) * covfoodka * covdoseka * covformka

    # IC50 is reported in ng/mL; the peripheral concentration below is in
    # mg/L, so convert exactly as the control stream does (THETA(6)/1000).
    ic50 <- exp(lic50) / 1000

    # ----- ODE system -----
    # covformf scales the depot LOSS rate only (the published
    # "loss from depot" parameterisation); covptstf scales the fraction
    # absorbed into the central compartment. Concentrations are written
    # inline; see implementation notes.
    d/dt(depot) <- -kad * depot * covformf
    d/dt(central) <-
      kad * depot *
        (2 - (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp))) *
        covptstf -
      q * (central / vc) + q * (peripheral1 / vp) -
      cl * (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp)) *
        (central / vc)
    d/dt(peripheral1) <- q * (central / vc) - q * (peripheral1 / vp)

    # ----- Observation and residual error -----
    Cc <- central / vc * 1000
    propSdEff <- propSd * (1 + e_inflam_propsd * inflam)
    Cc ~ prop(propSdEff)
  })
}
