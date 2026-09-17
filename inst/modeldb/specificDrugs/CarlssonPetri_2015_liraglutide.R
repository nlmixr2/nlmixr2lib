CarlssonPetri_2015_liraglutide <- function() {
  description <- "Liraglutide population PK model in pediatric (10-17 y) and adult subjects with type 2 diabetes (Carlsson Petri 2015 pooled three-trial analysis)"
  reference <- "Carlsson Petri KC, Jacobsen LV, Klein DJ. Comparable liraglutide pharmacokinetics in pediatric and adult populations with type 2 diabetes: a population pharmacokinetic analysis. Clin Pharmacokinet. 2015;54(6):663-670. doi:10.1007/s40262-014-0229-z"
  vignette <- "CarlssonPetri_2015_liraglutide"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")
  # Time in hours; subcutaneous doses in nmol; concentrations in nmol/L (nM).
  # The paper works natively in nmol: Methods Sect. 2.2 derives exposure as
  # AUC24 = dose (in nmol) / (CL/F), and Figs. 1 and 4 plot liraglutide
  # concentration in nM. Liraglutide MW = 3751.2 g/mol, so the clinical doses
  # convert as 0.3 mg = 79.97 nmol, 0.6 mg = 159.9 nmol, 1.2 mg = 319.9 nmol,
  # and 1.8 mg = 479.8 nmol. The same MW is used by the sibling
  # Overgaard_2016_liraglutide model.

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "liraglutide", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "liraglutide", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F only; reference weight 90 kg per Carlsson Petri 2015 Methods Sect. 2.2 ('approximate average body weight for patients with T2D in the trial populations analyzed here') and the Online Resource Table 1 CL/F equation, which centres on 90 kg. Observed range across the three pooled trials 57-214 kg (main Table 1); the covariate forest plot (Fig. 3) is drawn at 53 kg and 216 kg. Body weight was measured at randomization for Trial 1 and at sampling time for Trials 2 and 3.",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) is the paper's reference category",
      notes = "Carlsson Petri 2015 Online Resource Table 1 reports the log-scale contrast 'Cov CL-Gen' = +0.365, applied in the CL/F equation as exp(CovCLGen) when the subject is male. Implemented as exp(e_male_cl * (1 - SEXF)) so SEXF = 1 (female) evaluates to 1.00 (reference) and SEXF = 0 (male) evaluates to exp(0.365) = 1.441. The higher male CL/F reproduces the main text's 31 % lower male AUC24 (1/1.441 = 0.694). The paper's reference subject is explicitly 'an adult female subject weighing 90 kg' (Online Resource Table 1 footnote).",
      source_name = "Gen"
    ),
    CHILD = list(
      description = "Indicator for the pediatric age category (10-17 years), 0 = adult",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult) is the paper's reference category",
      notes = "Carlsson Petri 2015 treats age as a two-level category, pediatric versus adult, rather than as a continuous covariate: Methods Sect. 2.2 states that the narrow pediatric age range (10-17 y) and the gap to the youngest adult (33 y) made a continuous age covariate infeasible. Pediatric = Trial 1 (10-17 y); adult = Trials 2 and 3 (33-73 y). Online Resource Table 1 reports 'Cov CL-AGEgr' = +0.107, applied as exp(CovCLPaed) when the subject is pediatric, i.e. pediatric CL/F is 11 % higher and pediatric AUC24 is 10 % lower than an adult of the same weight and sex. The paper judges this effect NOT pharmacokinetically relevant (RSE 91 %; 90 % CI on the AUC24 ratio 0.78-1.03, only just outside the 0.80-1.25 bioequivalence window) but retains it in the full model, so it is encoded here as published.",
      source_name = "AGEgr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 57L, # 13 pediatric (Trial 1) + 12 adult (Trial 2) + 32 adult (Trial 3); Carlsson Petri 2015 Table 1
    n_studies = 3L, # Trial 1 NCT00943501 (pediatric), Trial 2 NCT00993304 (adult), Trial 3 NCT00873223 (adult)
    age_range = "10-73 years (Trial 1 pediatric 10-17; Trial 2 adult 54-73; Trial 3 adult 33-68)",
    weight_range = "57-214 kg overall (Trial 1 median 106, range 57-214; Trial 2 median 83, range 72-104; Trial 3 median 96, range 58-140)",
    weight_median = "Trial 1 106 kg, Trial 2 83 kg, Trial 3 96 kg",
    sex_female_pct = 40.4, # 23/57 female: 8/13 (62 %) Trial 1, 6/12 (50 %) Trial 2, 9/32 (28 %) Trial 3; Carlsson Petri 2015 Table 1
    disease_state = "Type 2 diabetes. Pediatric subjects (Trial 1) had HbA1c 6.5-11 %, fasting plasma glucose 6.1-13.3 mmol/L, and BMI above the 85th percentile for age and sex, treated with diet and exercise alone or with a stable metformin dose; the majority were post-pubertal by Tanner stage. Subjects with impaired renal function were excluded from all three trials.",
    dose_range = "0.3-1.8 mg once daily subcutaneously. Trial 1 escalated weekly 0.3 -> 0.6 -> 0.9 -> 1.2 -> 1.8 mg with PK sampling at 0.3, 0.6, 1.2 and 1.8 mg; Trials 2 and 3 sampled at steady state on 1.8 mg only. Protocol deviations contributed some 0.9 mg and 1.5 mg observations (Table 1 footnote a).",
    regions = "Multi-national; pediatric subjects were treated in Europe and the USA (Trial 1, NCT00943501).",
    trials = c("NCT00943501", "NCT00993304", "NCT00873223"),
    notes = "Baseline demographics and PK sampling schedules are in Carlsson Petri 2015 Table 1. The pediatric median body weight (106 kg) exceeds both adult trials, driven partly by one 214 kg subject. Trial 1 sampled -0.25, 2, 4, 8, 10, 11, 12, 14, 24 h for the first four subjects and -0.25, 2, 5, 8, 10, 13 h thereafter (protocol amendment shortening the clinic stay), plus 24, 48 and 72 h after the final week-5 dose."
  )

  ini({
    # Structural parameters. Carlsson Petri 2015 Online Resource (ESM) Table 1,
    # 'Full model' column. The structural model is one compartment with
    # first-order absorption and first-order elimination, carried over from the
    # adult liraglutide popPK model of ref [8] and re-fitted to the pooled
    # three-trial dataset (Methods Sect. 2.2).
    lka <- log(0.0657); label("Absorption rate constant (1/h)") # ESM Table 1 'KA', full model 0.0657, RSE 6 % (base model 0.0657, RSE 7 %). The table's unit column reads 'L/h', which is a typographical error: a first-order absorption rate constant is 1/h. Confirmed by magnitude, ln(2)/0.0657 = 10.6 h absorption half-life, consistent with liraglutide's ~8-12 h tmax and with the sibling Overgaard 2016 value of 0.0806 1/h.
    lcl <- log(1.06); label("Apparent clearance CL/F (L/h) for the reference subject") # ESM Table 1 'CL/F', full model 1.06, RSE 10 % (base model 1.49, RSE 9 %). Table footnote: 'Reference value corresponding to an adult female subject weighing 90 kg.'
    lvc <- log(15.3); label("Apparent volume of distribution V/F (L)") # ESM Table 1 'V/F', full model 15.3, RSE 19 % (base model 18.1, RSE 15 %). No covariate effects were estimated on V/F.
    lfdepot <- fixed(log(1)); label("Bioavailability of the subcutaneous depot (fraction)") # ESM Table 1 'F', reported as '1 (fixed)' in both the full and base models. CL/F and V/F are therefore apparent parameters and F carries no additional information; it is encoded explicitly to preserve the published fixed status.

    # Covariate effects on CL/F. ESM Table 1 footnote ** gives the equation
    #   CL_i/F = (CL/F) * (BWT/90 kg)^CovCLBWT
    #            * exp(CovCLGen if male) * exp(CovCLPaed if paediatric) * exp(eta_i)
    # Every coefficient below is independently falsified by an effect size the
    # main text prints in Results (see the trailing comments and the vignette's
    # source-trace table).
    e_wt_cl <- 0.929; label("Body weight power exponent on CL/F (unitless)") # ESM Table 1 'Cov CL-BWT', 0.929, RSE 15 %; applied as (WT/90)^e_wt_cl. Falsified against Results: AUC24 ratio at 53 kg = (53/90)^-0.929 = 1.635 versus the printed '64 % higher', and at 216 kg = (216/90)^-0.929 = 0.443 versus the printed '56 % lower'.
    e_male_cl <- 0.365; label("Log-scale male sex effect on CL/F") # ESM Table 1 'Cov CL-Gen', 0.365, RSE 29 %; applied as exp(e_male_cl * (1 - SEXF)) because female is the reference. Falsified against Results: male AUC24 ratio = exp(-0.365) = 0.694 versus the printed 'a 31 % lower drug exposure compared with the reference female subject of the same body weight'.
    e_child_cl <- 0.107; label("Log-scale pediatric age-category effect on CL/F") # ESM Table 1 'Cov CL-AGEgr', 0.107, RSE 91 %; applied as exp(e_child_cl * CHILD) because adult is the reference. Falsified against Results: pediatric AUC24 ratio = exp(-0.107) = 0.899 versus the printed '10 % decrease in AUC24 compared with an adult subject of the same weight and gender'.

    # Between-subject variability. ESM Table 1 reports BSV as a percent CV on
    # the NONMEM variance scale, i.e. omega^2 = (CV/100)^2, NOT the exact
    # log-normal omega^2 = log(1 + CV^2). The scale is settled by the Results
    # statement that 'the inclusion of covariates reduced between-subject
    # variability in CL/F by 50 % compared with the base model without
    # covariates': on the variance scale 1 - 0.34^2/0.48^2 = 49.8 %, which
    # rounds to the printed 50 %, whereas the exact log-normal reading gives
    # 47.2 % and the CV-scale reading gives 29.2 %.
    etalcl ~ 0.1156 # ESM Table 1 'BSV in CL/F' = 34 %CV in the full model (shrinkage 1 %); base model 48 %CV (shrinkage 11 %). omega^2 = 0.34^2 = 0.1156.
    etalvc ~ fixed(0) # Methods Sect. 2.2 states 'Between-subject variability parameters were included for CL/F and Vd/F (data not shown)', but ESM Table 1 reports a BSV estimate for CL/F only. The V/F variance is genuinely unpublished, so it is encoded as a structural zero rather than invented. See the vignette 'Assumptions and deviations' section.

    propSd <- 0.23; label("Proportional residual error (fraction)") # ESM Table 1 'Sigma', residual error (proportional), 23 %CV in both the full model (shrinkage 7 %) and the base model (shrinkage 8 %).
  })

  model({
    ka <- exp(lka)
    # CL/F equation, ESM Table 1 footnote **. Female and adult are the
    # reference categories, so the male and pediatric terms switch on via
    # (1 - SEXF) and CHILD respectively.
    cl <- exp(lcl + etalcl) *
      (WT / 90)^e_wt_cl *
      exp(e_male_cl * (1 - SEXF) + e_child_cl * CHILD)
    vc <- exp(lvc + etalvc)

    f(depot) <- exp(lfdepot)
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
