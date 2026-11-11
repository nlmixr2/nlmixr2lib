Saqr_2024_tacrolimus <- function() {
  description <- "Steady-state apparent oral clearance (CL/F) population model for immediate-release tacrolimus whole-blood troughs in adult kidney transplant recipients in the first 6 months post-transplant (Saqr 2024, DeKAF Genomics development cohort). The model has no compartments and no dose events: the trough is approximated by the average steady-state concentration, Css = total daily dose / (24 * CL/F) (Saqr 2024 Equations 1-2), with the current total daily dose supplied as the covariate DOSE_TAC_MGD. CL/F = 27.2 L/h in the reference patient (2 CYP3A loss-of-function alleles, no steroids, age 51.41 years, day 9 or later post-transplant, no diabetes, no calcium channel blocker, antiviral or ACE inhibitor), multiplied by eight composite CYP3A4/5 loss-of-function-allele-count x steroid-use group factors (0.68 to 2.65), diabetes at transplant (0.90), the first 9 days post-transplant (1.21), calcium channel blocker (0.93), prophylactic antiviral (0.96) and ACE inhibitor (1.04) co-medication, and a power effect of age ((AGE / 51.41)^-0.35). Two-level log-normal random effects on CL/F (between-subject and nested between-site) and combined proportional + additive residual error."
  reference <- "Saqr A, Al-Kofahi M, Mohamed M, Dorr C, Remmel RP, Onyeaghala G, Oetting WS, Guan W, Mannon RB, Matas AJ, Israni A, Jacobson PA. Steroid-tacrolimus drug-drug interaction and the effect of CYP3A genotypes. Br J Clin Pharmacol. 2024;90(11):2837-2848. doi:10.1111/bcp.16172"
  vignette <- "Saqr_2024_tacrolimus"
  units <- list(
    time = "h",
    dosing = "mg/day (no dose events: the current total daily dose is supplied as the covariate DOSE_TAC_MGD)",
    concentration = "ng/mL (whole-blood trough)"
  )

  covariateData <- list(
    DOSE_TAC_MGD = list(
      description = "Current total daily oral tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying; the total daily dose the patient was taking at the trough (Saqr 2024 Equation 2 'total daily dose tacrolimus'). Structural input, not a covariate effect: the predicted trough is DOSE_TAC_MGD / 24 / CL/F (mg/h divided by L/h, times 1000 for ng/mL). The NONMEM code supplies the same quantity as a dose RATE column (CSS = RATE / CL). Saqr 2024 Results: median 6 mg/day (range 0.1-45) across both cohorts.",
      source_name = "total daily dose (NONMEM RATE)"
    ),
    CYP3A5_STAR1_HET = list(
      description = "CYP3A5 genotype with exactly one functional *1 allele (one loss-of-function allele: *1/*3, *1/*6 or *1/*7)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (paired with CYP3A5_STAR1_HOM = 0: two CYP3A5 loss-of-function alleles)",
      notes = "Time-fixed germline genotype. Saqr 2024 counts CYP3A5*3 (rs776746), *6 (rs10264272) and *7 (rs41303343) as loss-of-function (LOF) alleles, so a *1 allele here is one carrying none of the three. Enters only through the composite LOF count reconstructed in model(): CYP3A5 LOF alleles = 2 - CYP3A5_STAR1_HET - 2 * CYP3A5_STAR1_HOM.",
      source_name = "LOF3A (composite LOF group x steroid code)"
    ),
    CYP3A5_STAR1_HOM = list(
      description = "CYP3A5*1/*1 genotype (no CYP3A5 loss-of-function allele)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (paired with CYP3A5_STAR1_HET = 0: two CYP3A5 loss-of-function alleles)",
      notes = "Time-fixed germline genotype; see CYP3A5_STAR1_HET. Enters only through the composite LOF count reconstructed in model().",
      source_name = "LOF3A (composite LOF group x steroid code)"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A4*1/*1)",
      notes = "Time-fixed germline genotype; counted as one CYP3A LOF allele in the composite LOF count. A carrier indicator cannot distinguish *22/*22 from *1/*22; this only matters for the 4 DeKAF patients with CYP3A5*3/*3 + CYP3A4*22/*22 (4 LOF alleles), whom Saqr 2024 pooled with the 3-LOF group anyway, so the reconstruction (LOF >= 3 -> 3/4-LOF group) is exact for every diplotype in Saqr 2024 Table S1.",
      source_name = "CYP3A4*22"
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid (prednisone or methylprednisolone) use at the time of the trough",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no steroid at the trough)",
      notes = "Time-varying per trough record (Saqr 2024 Methods 'Covariate analysis': 'At each trough measurement, the presence of steroids (yes/no) was determined'). Dose and duration were not recorded. Enters only through the LOF x steroid group factor, so the steroid effect differs by LOF group (Table 4: +11.8%, +5%, +3% and +10.3% CL/F for 0, 1, 2 and 3/4 LOF alleles).",
      source_name = "LOF3A (steroid component)"
    ),
    DIS_DIAB = list(
      description = "Diabetes at the time of transplant",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes)",
      notes = "Time-fixed. Multiplies CL/F by 0.90 (Saqr 2024 Table 3). 620/1608 (38.6%) of the development cohort (Table 1). Missing values were coded as no diabetes in the NONMEM code.",
      source_name = "DIAB"
    ),
    POD = list(
      description = "Days since kidney transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Enters only as the indicator POD < 9 ('Before Day 9 post-transplant', Saqr 2024 Table 3), which multiplies CL/F by 1.21; the reference category is day 9 or later (NONMEM DAYSPL = 1).",
      source_name = "DAYSPL"
    ),
    AGE = list(
      description = "Recipient age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE / 51.41)^-0.35 on CL/F, centred at the cohort median age 51.41 years (Saqr 2024 Table 3 and NONMEM code CLAGE = (AGE/51.41)**THETA(13)).",
      source_name = "AGE"
    ),
    CONMED_CCB = list(
      description = "Calcium channel blocker co-medication at the time of the trough",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CCB)",
      notes = "Time-varying per trough. Multiplies CL/F by 0.93 (Saqr 2024 Table 3). Agents were not recorded individually (pooled class). 37.8% of development-cohort troughs (Table 1).",
      source_name = "CCB"
    ),
    CONMED_ANTIVIRAL = list(
      description = "Prophylactic antiviral co-medication at the time of the trough",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no antiviral)",
      notes = "Time-varying per trough. Multiplies CL/F by 0.96 (Saqr 2024 Table 3). Saqr 2024 Methods: 'prophylactic antiviral agents'; the specific agents were not recorded (Discussion). 56.7% of development-cohort troughs (Table 1). Missing values were coded as no antiviral in the NONMEM code.",
      source_name = "ANTIVIRAL"
    ),
    CONMED_ACEI = list(
      description = "Angiotensin-converting enzyme inhibitor co-medication at the time of the trough",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ACE inhibitor)",
      notes = "Time-varying per trough. Multiplies CL/F by 1.04 (Saqr 2024 Table 3). 13.5% of development-cohort troughs (Table 1).",
      source_name = "ACE"
    ),
    SITEN = list(
      description = "Enrolling transplant centre (nesting level for the between-site random effect)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Integer centre index. Saqr 2024 Equation 3 adds an inter-site random effect on CL/F nested above the subject level (NONMEM $LEVEL UNIT=(2[1])); the DeKAF development cohort had seven centres. Consumed by the ini() nesting syntax (etalcl_site ~ ... | SITEN), not by model(); a simulation needs at least two distinct values.",
      source_name = "UNIT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1608L,
    n_studies = 1L,
    n_observations = 28718L,
    age_range = ">= 18 years",
    age_median = "51.41 years (reference age of the power effect)",
    sex_female_pct = 36.7,
    race_ethnicity = c(Caucasian = 78.1, `African American` = 16.6, Asian = 3.2, Hispanic = 2.1),
    disease_state = "Adult kidney transplant recipients on immediate-release tacrolimus maintenance immunosuppression in the first 6 months post-transplant; 59.4% living donor, 38.6% diabetic at transplant, 15.6% prior kidney transplant. Immunosuppression per centre protocol (early steroid withdrawal, steroid avoidance, or late steroid withdrawal with prednisone 5-10 mg/day).",
    dose_range = "Oral immediate-release tacrolimus twice daily, titrated to troughs of 8-12 ng/mL (months 1-3) and 6-10 ng/mL (months 3-6); total daily dose median 6 mg/day (range 0.1-45) across both cohorts.",
    regions = "United States (seven DeKAF transplant centres)",
    cyp3a_lof_groups = c(`0-LOF` = 4.4, `1-LOF` = 19.7, `2-LOF` = 67.1, `3/4-LOF` = 8.0, missing = 0.8),
    notes = "Development cohort: DeKAF Genomics (NCT00270712, enrolled 2005-2011), 1608 patients and 28718 troughs (Saqr 2024 Table 1); 63.6% of troughs were taken on steroids. External validation in GEN03 (NCT01714440, 854 patients, 15828 troughs) gave MPE -1.42 ng/mL, MPPE -17.8% and RMSE 2.6 ng/mL (Table 5). Troughs were sampled twice weekly for 8 weeks and twice monthly in months 3-6 (up to 24 per patient). Weight, serum creatinine, MDRD GFR, prior transplant and donor type were tested and not retained."
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a time-varying covariate on CL/F in stepwise covariate modelling and not retained (Saqr 2024 Methods 'Covariate analysis').",
      source_name = "body weight"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a time-varying covariate on CL/F and not retained (Saqr 2024 Methods).",
      source_name = "serum creatinine"
    )
  )

  ini({
    # Saqr 2024 Table 3. NONMEM: TVCL = THETA(1)/1000 with CSS = RATE/CL; the
    # /1000 converts the mg/h dose rate to a ng/mL trough, so THETA(1) is CL/F
    # in L/h.
    lcl <- log(27.2); label("Apparent oral clearance CL/F, reference patient (L/h)") # Table 3 'typical value of CL/F' 27.2 L/h (RSE 2%)

    # LOF x steroid group factors, relative to 2-LOF without steroids (Table 3)
    e_lof0_cl <- 2.37; label("CL/F factor, 0 CYP3A LOF alleles without steroids (unitless)") # Table 3 '0-LOF without steroids' 2.37
    e_lof0_steroid_cl <- 2.65; label("CL/F factor, 0 CYP3A LOF alleles with steroids (unitless)") # Table 3 '0-LOF with steroids' 2.65
    e_lof1_cl <- 1.78; label("CL/F factor, 1 CYP3A LOF allele without steroids (unitless)") # Table 3 '1-LOF without steroids' 1.78
    e_lof1_steroid_cl <- 1.87; label("CL/F factor, 1 CYP3A LOF allele with steroids (unitless)") # Table 3 '1-LOF with steroids' 1.87
    e_lof2_steroid_cl <- 1.03; label("CL/F factor, 2 CYP3A LOF alleles with steroids (unitless)") # Table 3 '2-LOF with steroids' 1.03
    e_lof34_cl <- 0.68; label("CL/F factor, 3 or 4 CYP3A LOF alleles without steroids (unitless)") # Table 3 '3/4-LOF without steroids' 0.68
    e_lof34_steroid_cl <- 0.75; label("CL/F factor, 3 or 4 CYP3A LOF alleles with steroids (unitless)") # Table 3 '3/4-LOF with steroids' 0.75

    # Clinical covariate factors (Table 3)
    e_diab_cl <- 0.90; label("CL/F factor, diabetes at transplant (unitless)") # Table 3 'Diabetes at time of transplant' 0.90
    e_pod_lt9_cl <- 1.21; label("CL/F factor, before day 9 post-transplant (unitless)") # Table 3 'Before Day 9 post-transplant' 1.21
    e_ccb_cl <- 0.93; label("CL/F factor, calcium channel blocker co-medication (unitless)") # Table 3 'CCB use' 0.93
    e_antiviral_cl <- 0.96; label("CL/F factor, prophylactic antiviral co-medication (unitless)") # Table 3 'Antiviral use' 0.96
    e_acei_cl <- 1.04; label("CL/F factor, ACE inhibitor co-medication (unitless)") # Table 3 'ACE inhibitor use' 1.04
    e_age_cl <- -0.35; label("Power exponent of (AGE / 51.41) on CL/F (unitless)") # Table 3 'Age' (Age/51.41)^-0.35

    # Two-level random effects on CL/F (Equation 3, CL/F = TVCL * exp(eta_IIV + eta_ISV)).
    # Table 3 reports CV%; converted with omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.449^2) # Table 3 'Interindividual variability' 44.9 CV%
    # Numeric literal: rxode2 evaluates an expression before '| SITEN' without
    # its outer log(), so log(1 + 0.189^2) would become 1.0357.
    etalcl_site ~ 0.0350978 | SITEN # Table 3 'Inter-site variability' 18.9 CV%; log(1 + 0.189^2)

    # Residual error, Equation 4 (NONMEM Y = F*(1+ERR(1)) + ERR(2)).
    propSd <- 0.272; label("Proportional residual error (fraction)") # Table 3 'Proportional residual unexplained variability' 27.2%
    addSd <- sqrt(8.5); label("Additive residual error (ng/mL)") # Table 3 'Additive residual unexplained variability' 8.5, read as a variance; see vignette
  })

  model({
    # Composite CYP3A loss-of-function allele count (Saqr 2024 Methods
    # 'Covariate analysis'): CYP3A5 *3/*6/*7 alleles plus CYP3A4*22.
    lof <- (2 - CYP3A5_STAR1_HET - 2 * CYP3A5_STAR1_HOM) + SNP_CYP3A4_RS35599367
    lof0 <- 0
    lof1 <- 0
    lof34 <- 0
    if (lof < 0.5) lof0 <- 1
    if (lof > 0.5 & lof < 1.5) lof1 <- 1
    if (lof > 2.5) lof34 <- 1
    lof2 <- 1 - lof0 - lof1 - lof34

    # LOF x steroid group factor; 2-LOF without steroids is the reference.
    f_lof_steroid <- lof0 * (e_lof0_cl * (1 - CONMED_STEROID) + e_lof0_steroid_cl * CONMED_STEROID) +
      lof1 * (e_lof1_cl * (1 - CONMED_STEROID) + e_lof1_steroid_cl * CONMED_STEROID) +
      lof2 * ((1 - CONMED_STEROID) + e_lof2_steroid_cl * CONMED_STEROID) +
      lof34 * (e_lof34_cl * (1 - CONMED_STEROID) + e_lof34_steroid_cl * CONMED_STEROID)

    pod_lt9 <- 0
    if (POD < 9) pod_lt9 <- 1

    cl <- exp(lcl + etalcl + etalcl_site) *
      f_lof_steroid *
      e_diab_cl^DIS_DIAB *
      e_pod_lt9_cl^pod_lt9 *
      e_ccb_cl^CONMED_CCB *
      e_antiviral_cl^CONMED_ANTIVIRAL *
      e_acei_cl^CONMED_ACEI *
      (AGE / 51.41)^e_age_cl

    # Equations 1-2: trough ~ Css = daily dose / (24 * CL/F); x 1000 converts
    # mg/L to ng/mL.
    Cc <- 1000 * DOSE_TAC_MGD / (24 * cl)
    Cc ~ prop(propSd) + add(addSd)
  })
}
