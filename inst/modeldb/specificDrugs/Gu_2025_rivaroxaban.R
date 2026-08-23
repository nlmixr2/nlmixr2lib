Gu_2025_rivaroxaban <- function() {
  description <- "Two-compartment oral population PK model with first-order absorption, an absorption lag time, and dose-level-dependent relative bioavailability for rivaroxaban in Chinese healthy volunteers and patients treated with radiofrequency ablation for non-valvular atrial fibrillation (Gu 2025)"
  reference <- "Gu F, Tang K, Zhang C, Hu M, Sun J, Yu X, Tian M, Zhang C, Chen Y. Population pharmacokinetic analysis of rivaroxaban in healthy volunteers and patients with radiofrequency ablation of non-valvular atrial fibrillation in China. Front Pharmacol. 2025;16:1562259. doi:10.3389/fphar.2025.1562259"
  vignette <- "Gu_2025_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalized:",
        "CRCL = [140 - AGE (years)] x WT (kg) x 0.85 (if female) / [72 x CREAT (mg/dL)]",
        "(Gu 2025 Table 1 footnote a). Enters CL/F as the power term (CRCL / 97.7)^1.53,",
        "with 97.7 mL/min the reference value printed in the Gu 2025 Results section 3.2",
        "final-model equation. Cohort values: healthy volunteers 125 +/- 17.8 (range 91.3-161)",
        "mL/min, NVAF patients 86.7 +/- 24.3 (range 33.6-158) mL/min (Gu 2025 Table 1).",
        "Same raw-Cockcroft-Gault usage as Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R,",
        "Wada_2023_sparsentan.R and Shu_2024_posaconazole.R."
      ),
      source_name        = "CRCL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with non-valvular atrial fibrillation treated by radiofrequency ablation, Gu 2025 Study 2)",
      notes              = paste(
        "1 = healthy Chinese volunteer from the bioequivalence study (Gu 2025 Study 1, n = 36);",
        "0 = NVAF patient after radiofrequency ablation (Gu 2025 Study 2, n = 105).",
        "Gu 2025 calls this the 'morbid state' covariate and reports the two typical clearances",
        "directly (TVCL = 6.48 L/h for healthy volunteers, 8.35 L/h for patients; Results section 3.2),",
        "so the structural lcl is anchored on the PATIENT state and e_dis_healthy_cl =",
        "log(6.48 / 8.35) restores the healthy-volunteer typical value, matching the canonical",
        "DIS_HEALTHY orientation used by Chen_2023_nemonoxacin.R and Galluppi_2021_ulotaront.R."
      ),
      source_name        = "morbid state"
    ),
    SNP_ABCB1_RS1045642_HOM = list(
      description        = "ABCB1 rs1045642 (c.3435C>T) homozygous-variant indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (GG or AG/CG, i.e. the pooled c.3435CC wild-type homozygotes and CT heterozygotes)",
      notes              = paste(
        "1 = AA genotype at rs1045642 (equivalently c.3435TT, the homozygous-variant stratum);",
        "0 = the pooled GG and AG/CG genotypes. Gu 2025 first screened rs1045642 as a four-category",
        "variable, found that only the AA stratum separated from the other three, and reclassified it",
        "as the two-category AA-vs-non-AA covariate used in the final model (Gu 2025 Results section 3.2).",
        "Genotype counts: healthy volunteers GG 18 / AG-CG 11 / AA 7; patients GG 42 / AG-CG 43 / AA 20",
        "(Gu 2025 Table 2). Note that the Gu 2025 Discussion labels the AA stratum the 'wild genotype',",
        "which conflicts with the standard dbSNP orientation for rs1045642 (A = the c.3435T variant allele,",
        "G = the c.3435C reference allele) and with the observed allele frequencies in this cohort",
        "(G = 0.60, A = 0.40 in Study 2, matching the Han-Chinese c.3435C frequency); the model file",
        "follows the unambiguous AA-vs-non-AA coding of the printed final-model equation, not the prose label."
      ),
      source_name        = "AA (ABCB1 rs1045642)"
    ),
    DOSE_10MG = list(
      description        = "10 mg rivaroxaban dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (15 mg reference dose level, when DOSE_20MG is also 0)",
      notes              = paste(
        "1 = the dose record is a 10 mg rivaroxaban dose. Relative bioavailability of the 10 mg dose",
        "group was estimated as 1.363 relative to the 15 mg reference and then FIXED because of the",
        "small 10 mg sample size (2 of 105 patients; Gu 2025 Results section 3.2 and Table 3).",
        "Paired with DOSE_20MG so that both indicators = 0 selects the 15 mg reference (F1 = 1)."
      ),
      source_name        = "F 10mg"
    ),
    DOSE_20MG = list(
      description        = "20 mg rivaroxaban dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (15 mg reference dose level, when DOSE_10MG is also 0)",
      notes              = paste(
        "1 = the dose record is a 20 mg rivaroxaban dose (all 36 healthy volunteers in Study 1 and",
        "2 of 105 patients in Study 2). Relative bioavailability 0.537 versus the 15 mg reference,",
        "estimated (%RSE 15.6; Gu 2025 Table 3). Paired with DOSE_10MG so that both indicators = 0",
        "selects the 15 mg reference (F1 = 1)."
      ),
      source_name        = "F 20mg"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the forward-selection covariate step; not significant (Gu 2025 Results section 3.2). Note that sex is nonetheless an input to the Cockcroft-Gault CRCL that IS retained."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2). Age is an input to the retained Cockcroft-Gault CRCL. Cohort 31.6 +/- 8.9 years (healthy) and 64.1 +/- 8.4 years (patients)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Weight entered the full model as an effect on the peripheral volume V3/F but did NOT meet the backward-elimination retention criterion and was removed (Gu 2025 Results section 3.2). No point estimate is published for the removed effect."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^12 cells/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2). No canonical register entry exists for platelet count; documented here only, since the effect was not retained."
    ),
    GLU = list(
      description = "Plasma glucose",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a covariate in its own right; not significant (Gu 2025 Results section 3.2). Creatinine is an input to the retained Cockcroft-Gault CRCL."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    CA = list(
      description = "Serum calcium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2). No canonical register entry exists for serum calcium; documented here only, since the effect was not retained."
    ),
    POT = list(
      description = "Serum potassium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    CPK = list(
      description = "Creatine kinase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "56 of 105 patients (53.3%) received a P-gp inhibitor (amiodarone, propafenone, or felodipine; Gu 2025 Table 1 and Discussion). Post hoc individual clearances were significantly lower in this group (t-test p = 0.046; Gu 2025 Figure 4), but the indicator was NOT retained as a covariate in the final model and no point estimate is published."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "6 of 105 patients (5.7%) received a CYP3A4 inhibitor (Gu 2025 Table 1). Post hoc individual clearances did not differ (t-test p = 0.877; Gu 2025 Figure 4) and the indicator was not retained in the final model."
    ),
    SNP_CYP3A4_RS2242480 = list(
      description = "CYP3A4 rs2242480 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "One of the ten genotyped SNPs screened as covariates; not significant (Gu 2025 Results section 3.2). Genotype counts in Gu 2025 Table 2."
    ),
    SNP_CYP3A4_RS2246709 = list(
      description = "CYP3A4 rs2246709 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_CYP3A4_RS3735451 = list(
      description = "CYP3A4 rs3735451 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_CYP3A5_RS776746 = list(
      description = "CYP3A5 rs776746 (CYP3A5*3) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_ABCB1_RS1128503 = list(
      description = "ABCB1 rs1128503 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_ABCB1_RS2032582 = list(
      description = "ABCB1 rs2032582 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2). Genotype distribution departed from Hardy-Weinberg equilibrium in both studies (Gu 2025 Table 2)."
    ),
    SNP_ABCB1_RS4148738 = list(
      description = "ABCB1 rs4148738 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_ABCB1_RS4728709 = list(
      description = "ABCB1 rs4728709 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    ),
    SNP_ABCG2_RS3114018 = list(
      description = "ABCG2 rs3114018 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not significant (Gu 2025 Results section 3.2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 141L,
    n_studies      = 2L,
    n_observations = 1506L,
    age_range      = "18-82 years (healthy volunteers 18-48; NVAF patients 34-82)",
    age_median     = "31.6 +/- 8.9 years (healthy volunteers); 64.1 +/- 8.4 years (patients)",
    weight_range   = "42.0-107.0 kg",
    weight_median  = "64.2 +/- 5.8 kg (healthy volunteers); 67.4 +/- 11.4 kg (patients)",
    sex_female_pct = 34.8,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy volunteers (bioequivalence study) pooled with patients treated by radiofrequency ablation for non-valvular atrial fibrillation",
    dose_range     = "Single oral 20 mg (healthy volunteers, Study 1); 10, 15, or 20 mg once daily orally (patients, Study 2)",
    regions        = "China (Huzhou Central Hospital)",
    renal_function = "Creatinine clearance 33.6-161 mL/min (Cockcroft-Gault); 42 patients with normal renal function, 50 with mild and 13 with moderate impairment (Gu 2025 Figure 3)",
    co_medication  = "P-glycoprotein inhibitors in 56/105 patients (53.3%); CYP3A4 inhibitors in 6/105 (5.7%)",
    notes          = paste(
      "Study 1 (CTR20202135): single-centre, single-dose, randomized, two-formulation, four-period",
      "crossover bioequivalence study in 36 healthy Chinese volunteers dosed 20 mg fasted, 18 samples",
      "per period over 48 h; only the reference-formulation (Xarelto) periods contributed the 1,296",
      "modelled concentrations. Study 2 (ChiCTR2500095918): real-world single-centre study of 105",
      "hospitalized NVAF patients after radiofrequency ablation, contributing 105 trough (30 min",
      "pre-dose) and 105 peak (2-4 h post-dose) concentrations. Baseline demographics per Gu 2025",
      "Table 1; sex_female_pct is the pooled 49/141 across both studies."
    )
  )

  ini({
    # Structural parameters -- reference subject is an NVAF PATIENT
    # (DIS_HEALTHY = 0) with CRCL = 97.7 mL/min, non-AA rs1045642 genotype
    # (SNP_ABCB1_RS1045642_HOM = 0), receiving the 15 mg reference dose
    # (DOSE_10MG = DOSE_20MG = 0). Values from Gu 2025 Table 3 and the
    # final-model equations printed in Results section 3.2.
    lcl     <- log(8.35);     label("Apparent clearance CL/F (L/h)")                          # Gu 2025 Table 3: CL/F = 8.35 L/h for patients (%RSE 15.8; bootstrap 8.25 [5.47-10.9])
    lvc     <- log(19.7);     label("Apparent central volume V2/F (L)")                        # Gu 2025 Table 3: V2/F = 19.7 L (%RSE 15.8; bootstrap 19.2 [12.3-24.3])
    lka     <- log(0.46);     label("Absorption rate constant ka (1/h)")                       # Gu 2025 Table 3: Ka = 0.46 1/h (%RSE 11.6; bootstrap 0.462 [0.373-0.525])
    lq      <- log(7.64);     label("Apparent intercompartmental clearance Q/F (L/h)")         # Gu 2025 Table 3: Q/F = 7.64 L/h (%RSE 10.8; bootstrap 7.51 [5.55-9.19])
    lvp     <- log(71.8);     label("Apparent peripheral volume V3/F (L)")                     # Gu 2025 Table 3: V3/F = 71.8 L (%RSE 18.9; bootstrap 71.8 [51.7-88.6])
    ltlag   <- log(0.168);    label("Absorption lag time (h)")                                 # Gu 2025 Table 3: ALAG = 0.168 h (%RSE 14.8; bootstrap 0.169 [0.123-0.196])
    lfdepot <- fixed(log(1)); label("Reference relative bioavailability F1 at the 15 mg dose level (unitless)")  # Gu 2025 Table 3: F 15mg = 1 FIX

    # Covariate effects on CL/F. The Gu 2025 Results section 3.2 final-model
    # equation is
    #   CL_i = TVCL * exp(1.53 * ln(CRCL_i / 97.7) - AA * 0.204 + eta_CL)
    # i.e. a power term on CRCL and a log-additive shift for the AA genotype.
    # Table 3 reports the genotype effect on the multiplicative scale as
    # CL_A642 = 0.815 = exp(-0.204), which is the cross-check between the two.
    # NOTE: this exponent is encoded exactly as published, but it does not
    # reproduce Gu 2025 Figures 2, 3 or 4 -- it predicts a much steeper renal
    # gradient than the paper's own results show. See the vignette
    # "Assumptions and deviations" section for the three-way check.
    e_crcl_cl                    <-  1.53;   label("Power exponent of creatinine clearance on CL/F, referenced to 97.7 mL/min (unitless)")  # Gu 2025 Table 3 CL_crcl = 1.53 (%RSE 11.9; bootstrap 1.53 [1-1.97]); reference 97.7 mL/min from the Results section 3.2 equation
    e_snp_abcb1_rs1045642_hom_cl <- -0.204;  label("log-scale multiplier of the ABCB1 rs1045642 AA genotype on CL/F (unitless; exp(-0.204) = 0.815)")  # Gu 2025 Results section 3.2 equation; Table 3 CL_A642 = 0.815 (%RSE 7.29; bootstrap 0.821 [0.615-0.96])

    # Morbid state ("healthy volunteer" vs "patient"). Gu 2025 reports the two
    # typical clearances rather than a coefficient: TVCL = 6.48 L/h for healthy
    # volunteers and 8.35 L/h for patients (Results section 3.2, Table 3). lcl is
    # anchored on the patient state, so the healthy shift is log(6.48 / 8.35).
    e_dis_healthy_cl <- log(6.48 / 8.35); label("log-scale multiplier of healthy-participant status on CL/F (unitless; exp(e_dis_healthy_cl) = 0.776)")  # Gu 2025 Table 3: CL/F = 6.48 L/h for healthy volunteers (%RSE 15.7; bootstrap 6.34 [4.89-7.38]) vs 8.35 L/h for patients

    # Dose-level relative bioavailability, referenced to the 15 mg dose group.
    e_dose_10mg_fdepot <- fixed(log(1.363)); label("log-scale multiplier of the 10 mg dose level on F1 (unitless; exp(e_dose_10mg_fdepot) = 1.363)")  # Gu 2025 Table 3: F 10mg = 1.363 FIX (estimated at the base-model stage, then held constant because only 2 of 105 patients received 10 mg)
    e_dose_20mg_fdepot <- log(0.537);        label("log-scale multiplier of the 20 mg dose level on F1 (unitless; exp(e_dose_20mg_fdepot) = 0.537)")  # Gu 2025 Table 3: F 20mg = 0.537 (%RSE 15.6; bootstrap 0.524 [0.413-0.609])

    # IIV. Gu 2025 Table 3 reports each IIV as a bare percentage in a column
    # shared with the proportional residual error, whose printed 26.3% is
    # unambiguously sqrt(sigma^2) x 100; the IIV percentages are therefore read
    # on the same scale, omega = pct / 100, and the variances below are pct^2.
    # See the vignette "Assumptions and deviations" section for the full
    # convention check (the alternative lognormal-CV reading,
    # omega^2 = log(1 + CV^2), changes omega by at most 12%).
    etalcl   ~ 0.126025   # Gu 2025 Table 3: IIV_CL   = 35.5% (%RSE 6.69; bootstrap 34.8 [25-41]);    omega^2 = 0.355^2
    etalvc   ~ 0.350464   # Gu 2025 Table 3: IIV_V2   = 59.2% (%RSE 13.2; bootstrap 57.7 [19.7-72.5]); omega^2 = 0.592^2
    etalka   ~ 0.077284   # Gu 2025 Table 3: IIV_Ka   = 27.8% (%RSE 40.9; bootstrap 24.2 [0.313-40.2]); omega^2 = 0.278^2
    etalq    ~ 0.409600   # Gu 2025 Table 3: IIV_Q    = 64.0% (%RSE 14.0; bootstrap 63.1 [48.2-73.5]); omega^2 = 0.640^2
    etalvp   ~ 0.434281   # Gu 2025 Table 3: IIV_V3   = 65.9% (%RSE 11.4; bootstrap 64.6 [41.6-79.6]); omega^2 = 0.659^2
    etaltlag ~ 0.602176   # Gu 2025 Table 3: IIV_ALAG = 77.6% (%RSE 13.7; bootstrap 76.2 [52.2-95.4]); omega^2 = 0.776^2

    # Residual error -- proportional only. The Gu 2025 Methods describe a
    # combined proportional + additive residual model, but Table 3 of the final
    # model reports the proportional term alone, so the additive term was not
    # retained.
    propSd <- 0.263; label("Proportional residual error (fraction)")  # Gu 2025 Table 3: delta = 26.3% (%RSE 1.32; bootstrap 26.1 [20.4-29.7])
  })

  model({
    # Individual parameters. Reference: NVAF patient (DIS_HEALTHY = 0),
    # CRCL = 97.7 mL/min, SNP_ABCB1_RS1045642_HOM = 0, 15 mg dose.
    cl <- exp(lcl + etalcl) *
      (CRCL / 97.7)^e_crcl_cl *
      exp(e_snp_abcb1_rs1045642_hom_cl * SNP_ABCB1_RS1045642_HOM) *
      exp(e_dis_healthy_cl * DIS_HEALTHY)
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)
    tlag <- exp(ltlag + etaltlag)

    # Relative bioavailability: 1 at the 15 mg reference dose level, 1.363 at
    # 10 mg and 0.537 at 20 mg (Gu 2025 Table 3 / Results section 3.2).
    fdepot <- exp(lfdepot) *
      exp(e_dose_10mg_fdepot * DOSE_10MG) *
      exp(e_dose_20mg_fdepot * DOSE_20MG)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Dose is in mg and volumes are in L, so central / vc is mg/L; x 1000
    # converts to the ng/mL used throughout Gu 2025.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
