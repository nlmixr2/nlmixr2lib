Fediuk_2021_varenicline <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "elimination for oral varenicline in adolescent smokers aged 12-20 years",
    "(Fediuk 2021), pooled from two phase 1 studies and one phase 4 study.",
    "Apparent clearance and apparent volume scale with body weight (power",
    "functions on WT/70) and carry multiplicative race (Black, Other) and",
    "female-sex factors; the phase 1 and phase 4 studies have separate",
    "combined additive + proportional residual errors."
  )
  reference <- paste(
    "Fediuk DJ, Sweeney K, Sahasrabudhe V, McRae T, Byon W.",
    "Population pharmacokinetics and exposure-response analyses of",
    "varenicline in adolescent smokers.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(7):769-781.",
    "doi:10.1002/psp4.12645"
  )
  vignette <- "Fediuk_2021_varenicline"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function (WT/70)^e on both CL/F and V/F (Fediuk 2021 Methods covariate equations; Text S1 'NWT = BWT/70'). Reference subject 70 kg. Cohort median 62.1 kg, range 35.0-121 kg (Table 2).",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male (the Fediuk 2021 reference subject is a 70 kg White male)",
      notes = "Text S1 derives 'NSEX = SEX - 1' from a SEX column coded 1 = male, 2 = female, so NSEX is identically the canonical SEXF (Methods: 'NSEX describes male = 0 and female = 1'). Power-of-indicator factor on CL/F and V/F.",
      source_name = "NSEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White (the reference race)",
      notes = "Text S1 'IF (RACE.EQ.2) RACE2=1'. Power-of-indicator factor on CL/F and V/F. 31 of 218 subjects (14.2%; Table 2), 19 of them from phase 1 study 1.",
      source_name = "RACE2"
    ),
    RACE_OTHER = list(
      description = "Composite 'Other' race indicator pooling Asian, Hispanic, American Indian and mixed race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White (the reference race)",
      notes = "Text S1 'IF (RACE.EQ.3) RACE3=1' and 'IF (RACE.EQ.4) RACE3=1'. Methods: 'other (Asian, Hispanic, American Indian, and mixed race [RACE3])'; Results: Asian and 'other' were pooled because the phase 1 studies had only one Asian subject. 45 of 218 subjects (Table 2: Asian 26 + Other 19). Same composite as the Ravva 2009 adult varenicline model.",
      source_name = "RACE3"
    ),
    STUDY_PHASE4 = list(
      description = "Phase 4 study-stratum indicator: 1 = record from the phase 4 efficacy/safety study NCT01312909, 0 = record from one of the two phase 1 PK studies",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = phase 1 studies (A3051029 single dose, NCT00463918 multiple dose)",
      notes = "Selects the residual-error magnitudes ONLY (Text S1 $ERROR: 'IND = 0; IF(PROT.EQ.1073) IND = 1', with EPS(1)/EPS(2) for phase 1 and EPS(3)/EPS(4) for phase 4). It touches no structural or covariate parameter, so the typical-value prediction is identical for either value. Use 0 to reproduce the intensively sampled phase 1 profiles and 1 to reproduce the sparse, outpatient phase 4 scatter.",
      source_name = "PROT (IND)"
    )
  )

  # Screened or collected by the source analysis but not in the final model;
  # documentation only, never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Not tested: 'since BWT and age were correlated, age was not included as a covariate on CL/F or V/F' (Fediuk 2021 Methods). Text S1 computes NAGE = AGE/16 but never uses it. Cohort median 16 years, range 12-20 (Table 2).",
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Not tested: 'CrCl ... was not included as a covariate on CL/F since BWT and CrCl were correlated and most adolescents were expected to have normal renal function' (Fediuk 2021 Methods). Text S1 computes PBCCL = BCCL/125 (capped at 150 mL/min) but never uses it. Cohort median 124 mL/min, range 51.7-257 (Table 2).",
      source_name = "BCCL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "varenicline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "varenicline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 218L,
    n_studies = 3L,
    n_observations = 1097L,
    age_range = "12-20 years",
    age_median = "16 years",
    weight_range = "35.0-121 kg",
    weight_median = "62.1 kg",
    sex_female_pct = 37.6,
    race_ethnicity = c(White = 65.1, Black = 14.2, Asian = 11.9, Other = 8.72),
    renal_function = "Creatinine clearance median 124 mL/min (range 51.7-257); about 99% normal renal function (CrCl > 80 mL/min), 1% mild impairment.",
    disease_state = "Adolescent cigarette smokers: healthy adolescent smokers in the two phase 1 studies; nicotine-dependent adolescent smokers (Fagerstrom Test for Nicotine Dependence score >= 4) motivated to quit in the phase 4 study.",
    dose_range = "Oral varenicline 0.5 or 1 mg single dose (phase 1 study 1); 0.5 mg q.d. or 0.5 mg b.i.d. for body weight <= 55 kg and 0.5 mg b.i.d. or 1 mg b.i.d. for body weight > 55 kg, for 14 days (phase 1 study 2) or 12 weeks (phase 4), with a 1- or 2-week up-titration.",
    regions = "Multicenter (Pfizer-sponsored studies A3051029, NCT00463918, NCT01312909).",
    notes = "Baseline demographics from Fediuk 2021 Table 2 ('popPK analysis', 'Total' column): 22 subjects from phase 1 study 1 (single dose, intensive sampling to 48 h), 57 from phase 1 study 2 (multiple dose, intensive sampling on days 1, 8 and 14), and 139 from the phase 4 study (one random sample at weeks 3, 6 and 12). Fewer than 10% of observations were below the 0.1 ng/mL LLOQ. Estimation by NONMEM 7.3 FOCE-I."
  )

  # nlmixr2 takes one residual-SD symbol per output, so the study-stratum SDs
  # are separate ini() parameters combined into one symbol inside model() with
  # the STUDY_PHASE4 indicator (the Rich_2026_momelotinib.R pattern), declared
  # here so checkModelConventions() accepts the stratum suffixes.
  paper_specific_residual_sds <- c(
    "propSdPh1",
    "addSdPh1",
    "propSdPh4",
    "addSdPh4"
  )

  ini({
    # Structural parameters -- reference subject: 70 kg White male
    # (Fediuk 2021 Methods: 'The reference subject was defined as a 70 kg, white male').
    lka <- log(0.860); label("First-order absorption rate constant ka (1/h)") # Table 3 'ka (h-1)' = 0.860 (RSE 12.3%); Text S1 THETA(3)
    lcl <- log(12.5); label("Apparent clearance CL/F for the reference subject (L/h)") # Table 3 'CL/F (L/h)' = 12.5 (RSE 5.06%); Text S1 THETA(1)
    lvc <- log(231); label("Apparent volume of distribution V/F for the reference subject (L)") # Table 3 'V/F (L)' = 231 (RSE 5.02%); Text S1 THETA(2)

    # Body-weight power exponents on WT/70
    e_wt_cl <- 0.567; label("Power exponent of WT/70 on CL/F (unitless)") # Table 3 CL/F 'Body weight' = 0.567 (RSE 23.3%); Text S1 THETA(4)
    e_wt_vc <- 0.872; label("Power exponent of WT/70 on V/F (unitless)") # Table 3 V/F 'Body weight' = 0.872 (RSE 10.3%); Text S1 THETA(5)

    # Categorical covariate effects: power-of-indicator multipliers, theta^indicator
    e_race_black_cl <- 1.01; label("Multiplicative factor on CL/F for Black vs White race (unitless)") # Table 3 CL/F 'Black race' = 1.01 (RSE 8.59%); Text S1 THETA(6)
    e_race_other_cl <- 1.12; label("Multiplicative factor on CL/F for Other vs White race (unitless)") # Table 3 CL/F 'Other race' = 1.12 (RSE 6.57%); Text S1 THETA(7)
    e_sexf_cl <- 0.850; label("Multiplicative factor on CL/F for female vs male (unitless)") # Table 3 CL/F 'Female sex' = 0.850 (RSE 5.87%); Text S1 THETA(8)
    e_race_black_vc <- 0.757; label("Multiplicative factor on V/F for Black vs White race (unitless)") # Table 3 V/F 'Black race' = 0.757 (RSE 5.46%); Text S1 THETA(9)
    e_race_other_vc <- 0.854; label("Multiplicative factor on V/F for Other vs White race (unitless)") # Table 3 V/F 'Other race' = 0.854 (RSE 4.78%); Text S1 THETA(10)
    e_sexf_vc <- 0.861; label("Multiplicative factor on V/F for female vs male (unitless)") # Table 3 V/F 'Female sex' = 0.861 (RSE 4.02%); Text S1 THETA(11)

    # IIV: full 3x3 block on CL/F, V/F, ka (exponential etas), NONMEM
    # ETA(1)/ETA(2)/ETA(3) order. Table 3 final estimates; the Text S1
    # $OMEGA initials print the two covariances as -0.00161 and -0.0581.
    etalcl + etalvc + etalka ~ c(
      0.102,
      -0.00162, 0.0182,
      -0.0582, 0.0307, 0.174
    ) # Table 3 omega2 CL/F 0.102, V/F 0.0182, ka 0.174; COV CL/F,V/F -0.00162, CL/F,ka -0.0582, V/F,ka 0.0307

    # Residual error: Text S1 $ERROR Y = F*(1+EPS1)+EPS2 (phase 1) or
    # F*(1+EPS3)+EPS4 (phase 4), diagonal $SIGMA -> combined additive +
    # proportional with summed variances (nlmixr2 default combined2).
    # Table 3 prints sigma^2; the SDs below are its square roots, as printed in
    # the Table 3 footnote (additive) and the Results text (proportional).
    propSdPh1 <- 0.291; label("Proportional residual SD, phase 1 studies (fraction)") # Table 3 'Phase 1 proportional' sigma2 = 0.0847; Results '29.1%'
    addSdPh1 <- 0.240; label("Additive residual SD, phase 1 studies (ng/mL)") # Table 3 'Phase 1 additive' sigma2 = 0.0577 (0.240 ng/mL)
    propSdPh4 <- 0.432; label("Proportional residual SD, phase 4 study (fraction)") # Table 3 'Phase 4 proportional' sigma2 = 0.187; Results '43.2%'
    addSdPh4 <- 0.580; label("Additive residual SD, phase 4 study (ng/mL)") # Table 3 'Phase 4 additive' sigma2 = 0.336 (0.580 ng/mL)
  })

  model({
    # Individual parameters (Fediuk 2021 Methods covariate equations; Text S1 $PK)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      e_race_black_cl^RACE_BLACK *
      e_race_other_cl^RACE_OTHER *
      e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc *
      e_race_black_vc^RACE_BLACK *
      e_race_other_vc^RACE_OTHER *
      e_sexf_vc^SEXF

    kel <- cl / vc

    # One-compartment model with first-order absorption (ADVAN2 TRANS2)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L; x1000 -> ng/mL (Text S1 'S2 = V/1000')
    Cc <- 1000 * central / vc

    # Study-stratum residual error (Text S1 $ERROR IND switch)
    propSdCc <- propSdPh4 * STUDY_PHASE4 + propSdPh1 * (1 - STUDY_PHASE4)
    addSdCc <- addSdPh4 * STUDY_PHASE4 + addSdPh1 * (1 - STUDY_PHASE4)
    Cc ~ add(addSdCc) + prop(propSdCc)
  })
}
