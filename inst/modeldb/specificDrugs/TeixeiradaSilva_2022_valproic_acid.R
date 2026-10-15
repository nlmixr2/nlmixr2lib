TeixeiradaSilva_2022_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption for total serum valproic acid in Caucasian (Spanish) paediatric and adult outpatients followed by therapeutic drug monitoring (Teixeira-da-Silva 2022 final model). Apparent clearance scales allometrically with total body weight (fixed exponent 0.75), carries a small power effect of age centred on 15 years, and is increased by concomitant phenytoin (x1.640), phenobarbital (x1.386) and carbamazepine (x1.521). Apparent volume is fixed at 14 L per 70 kg (linear in weight) and the formulation-specific absorption rate constants are fixed from the literature (oral solution 2.64, gastro-resistant tablet 0.78, modified-release coated tablet 0.38 1/h)."
  reference <- "Teixeira-da-Silva P, Perez-Blanco JS, Santos-Buelga D, Otero MJ, Garcia MJ. Population Pharmacokinetics of Valproic Acid in Pediatric and Adult Caucasian Patients. Pharmaceutics. 2022;14(4):811. doi:10.3390/pharmaceutics14040811. PMCID PMC9031051. Final-model estimates from Table 2 (development dataset) and Equations (1)-(2); model structure from the NONMEM control stream in Supplementary Table S2."
  vignette <- "TeixeiradaSilva_2022_valproic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = "70 kg (Equations (1) and (2))",
      notes = "Allometric scaling fixed a priori: exponent 0.75 on CL/F and 1 on V/F, both normalised to 70 kg (Methods 2.3; Supplementary Table S2 'TVCL = THETA(1) * (TBW/70)**0.75', 'TVV = THETA(2) * (TBW/70)**1'). Development cohort 6.70-125.00 kg, median 60.00 kg (Table 1).",
      source_name = "TBW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = "15 years (Equation (1))",
      notes = "Power effect on CL/F centred on 15 years, (AGE/15)^-0.0154 (Equation (1); Supplementary Table S2 'CLAGE = ((AGE/15)**THETA(3))'). The centring value is the cut-off the authors saw on a plot of CL/F against age (Results 3.1). The exponent is imprecise (64.3% RSE; bootstrap 95% CI -0.034 to 0.004 crosses zero), and its effect is small: about +4% at 1 year and -1% at 35 years. Development cohort 0.11-89.42 years (Table 1).",
      source_name = "AGE"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Proportional effect on CL/F, (1 + 0.640)^PHT = 1.640 (Equation (1); Supplementary Table S2 'CLPHT = ( 1 + THETA(6))'). 18 of 836 development patients (39 samples) received phenytoin (Supplementary Table S1). Only mono- or dual therapy was admitted, so the three inducer indicators are mutually exclusive in the source data.",
      source_name = "PHT"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes = "Proportional effect on CL/F, (1 + 0.386)^PB = 1.386 (Equation (1); Supplementary Table S2 'CLPB = ( 1 + THETA(5))'). 19 of 836 development patients (37 samples) received phenobarbital (Supplementary Table S1).",
      source_name = "PB"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes = "Proportional effect on CL/F, (1 + 0.521)^CBZ = 1.521 (Equation (1); Supplementary Table S2 'CLCBZ = ( 1 + THETA(4))'). Table 2 prints the coefficient as 0.512, a transposition of the 0.521 that Equation (1) and the Discussion ('an increase CL/F of about 52%', and CL/F ratios of 0.818 / 0.538 = 1.520 and 2.258 / 1.484 = 1.522) both use; the equation value is encoded. 58 of 836 development patients (139 samples) received carbamazepine (Supplementary Table S1).",
      source_name = "CBZ"
    ),
    FORM_TABLET = list(
      description = "Gastro-resistant (enteric-coated) valproic acid tablet indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral solution / syrup, Ka = 2.64 1/h)",
      notes = "Selects the FIXED gastro-resistant-tablet absorption rate constant Ka = 0.78 1/h (Methods 2.3 and Table 2 footnote; Supplementary Table S2 'IF (FFS.EQ.2) KA=0.78'). The source FFS column is four-level {0 = not available, 1 = oral solution, 2 = gastro-resistant tablet, 3 = modified-release coated tablet}; FFS = 0 was also given Ka = 0.78, so a record with unknown formulation maps to FORM_TABLET = 1. Methods 2.3 states that missing formulation information was imputed as oral solution below 12 years of age and gastro-resistant tablet from 12 years. Pairs with FORM_VPA_SR; both 0 = oral solution.",
      source_name = "FFS"
    ),
    FORM_VPA_SR = list(
      description = "Modified-release (prolonged-release) coated valproic acid tablet indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral solution / syrup, Ka = 2.64 1/h)",
      notes = "Selects the FIXED modified-release-coated-tablet absorption rate constant Ka = 0.38 1/h (Methods 2.3 and Table 2 footnote; Supplementary Table S2 'IF (FFS.EQ.3) KA=0.38'). Available presentations were 300 mg and 500 mg coated prolonged-release tablets (Methods 2.2). Set FORM_TABLET = 0 when FORM_VPA_SR = 1.",
      source_name = "FFS"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Statistically significant on CL/F in the stepwise forward step (p < 0.005) but not retained in the final model because it changed CL/F by less than 10% (Results 3.1); no coefficient is reported. 385 of 836 development patients (46.1%) were female (Table 1).",
      source_name = "SEX"
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Statistically significant on CL/F in the stepwise forward step (p < 0.005) but not retained because it changed CL/F by less than 10% (Results 3.1 and Discussion); no coefficient is reported. 50 of 836 development patients received lamotrigine (Supplementary Table S1).",
      source_name = "LTG"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 836,
    n_studies = 1,
    n_observations = 1751,
    age_range = "0.11-89.42 years",
    age_median = "32.42 years",
    weight_range = "6.70-125.00 kg",
    weight_median = "60.00 kg",
    sex_female_pct = 46.1,
    race_ethnicity = "Caucasian (Spanish)",
    disease_state = "Outpatients receiving valproic acid (epilepsy and other neurological or psychiatric indications) in mono- or dual antiepileptic therapy",
    dose_range = "150-4500 mg/day (median 1000 mg/day); oral solution 200 mg/mL, gastro-resistant tablets 200 and 500 mg, prolonged-release coated tablets 300 and 500 mg",
    regions = "Spain (therapeutic drug monitoring programme, University Hospital of Salamanca)",
    co_medication = "Development dataset: carbamazepine 58, phenytoin 18, phenobarbital 19, lamotrigine 50, topiramate 13, ethosuximide 3, clobazam 8, primidone 1, other dual therapies 280 patients (Supplementary Table S1). At most one antiepileptic in addition to valproic acid.",
    age_groups = "Development dataset by age: 28 days-2 years 33 patients (47 samples), 2-11 years 208 (419), 12-18 years 70 (157), > 18 years 525 (1128) (Supplementary Table S1)",
    notes = "Development dataset of 836 patients and 1751 serum samples (Table 1 and Supplementary Table S1; the Abstract says 776 patients, which is the external-evaluation sample count). External evaluation dataset of 368 patients and 776 samples. 97.5% of samples were steady-state troughs from routine therapeutic drug monitoring, measured by fluorescence polarisation immunoassay (Abbott TDx). Fitted in NONMEM 7.5.1 with FOCEI and log-transform-both-sides. Patients with BMI outside 16.0-39.9 kg/m2 were excluded."
  )

  ini({
    # ----------------------------------------------------------------
    # Absorption - FIXED from the literature (Methods 2.3, Table 2
    # footnote, Supplementary Table S2). Oral solution is the reference
    # formulation, so lka is the solution value and each tablet
    # indicator carries the log-ratio shift -- the same encoding as the
    # sibling Zhang_2023_valproic_acid_base.R.
    # ----------------------------------------------------------------
    lka <- fixed(log(2.64)); label("Absorption rate constant, oral solution reference (1/h)") # Methods 2.3 and Table 2 footnote (Ka = 2.64 1/h for oral solution, fixed)
    e_form_tablet_ka <- fixed(log(0.78 / 2.64)); label("Log-ratio shift on Ka for gastro-resistant tablet vs oral solution") # Methods 2.3 and Table 2 footnote (Ka = 0.78 1/h for gastro-resistant tablets, fixed)
    e_form_vpa_sr_ka <- fixed(log(0.38 / 2.64)); label("Log-ratio shift on Ka for modified-release coated tablet vs oral solution") # Methods 2.3 and Table 2 footnote (Ka = 0.38 1/h for modified-release coated tablets, fixed)

    # ----------------------------------------------------------------
    # Structural parameters (apparent, total serum valproic acid).
    # ----------------------------------------------------------------
    lcl <- log(0.646); label("Apparent clearance CL/F for a 70 kg, 15-year-old patient on monotherapy (L/h)") # Table 2 final model, development dataset (CL/F = 0.646, RSE 1.20%); Equation (1)
    lvc <- fixed(log(14)); label("Apparent volume of distribution V/F for a 70 kg patient (L)") # Table 2 footnote and Equation (2) (V/F fixed at 14 L for a 70 kg patient)

    # ----------------------------------------------------------------
    # Allometric exponents, fixed a priori (Methods 2.3).
    # ----------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL/F (unitless)") # Methods 2.3 and Equation (1) (exponent 0.75, fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent on (WT/70) for V/F (unitless)") # Methods 2.3 and Equation (2) (exponent 1, fixed)

    # ----------------------------------------------------------------
    # Covariate effects on CL/F. Supplementary Table S2 codes each
    # comedication as (1 + THETA) when present and 1 when absent, and
    # age as (AGE/15)**THETA(3).
    # ----------------------------------------------------------------
    e_age_cl <- -0.0154; label("Power exponent on (AGE/15) for CL/F (unitless)") # Table 2 row 'AGE' = -0.0154 (RSE 64.3%); Equation (1)
    e_conmed_pht_cl <- 0.640; label("Fractional increase in CL/F with concomitant phenytoin (unitless)") # Table 2 row 'PHT' = 0.640 (RSE 24.2%); Equation (1) factor 1.640
    e_conmed_pb_cl <- 0.386; label("Fractional increase in CL/F with concomitant phenobarbital (unitless)") # Table 2 row 'PB' = 0.386 (RSE 23.2%); Equation (1) factor 1.386
    e_conmed_cbz_cl <- 0.521; label("Fractional increase in CL/F with concomitant carbamazepine (unitless)") # Equation (1) factor 1.521 and Discussion 'about 52%'; Table 2 row 'CBZ' prints the transposed 0.512 (RSE 13.4%)

    # ----------------------------------------------------------------
    # IIV. Exponential IIV on CL/F only; Supplementary Table S2 fixes
    # the V/F eta variance to 0 ('$OMEGA ... 0 FIX'). That zero-variance
    # term is omitted rather than written as a fixed(0) eta (same
    # treatment as Zhang_2024_valproic_acid.R). Table 2 reports the CL/F
    # IIV as a CV%, converted as omega^2 = log(CV^2 + 1):
    #   log(0.268^2 + 1) = 0.06938
    # ----------------------------------------------------------------
    etalcl ~ 0.06938 # Table 2 row 'IIV_CL/F (%)' = 26.8 (RSE 5.50%, shrinkage 19.0%); bootstrap median 26.6

    # ----------------------------------------------------------------
    # Residual error. Log-transform-both-sides with an additive error on
    # log(concentration) (Methods 2.3; Supplementary Table S2
    # 'Y = IPRED + EPS(1)' with IPRED = LOG(F)), encoded as lnorm().
    # Table 2 prints RUV = 57.7% for the development fit, but that value
    # lies far outside its own bootstrap 95% CI (25.8-30.4%), and the
    # individual-prediction RMSE of about 20% and population-prediction
    # RMSE of about 37% in Table 3 are consistent only with a residual CV
    # near 28% (see the vignette). The bootstrap median of 28.1% is used,
    # converted as sd = sqrt(log(CV^2 + 1)):
    #   sqrt(log(0.281^2 + 1)) = 0.2757
    # ----------------------------------------------------------------
    expSd <- 0.2757; label("Residual SD on the log scale (log-transform-both-sides)") # Table 2 row 'RUV (%)' bootstrap median 28.1 (95% CI 25.8-30.4); final-model column prints 57.7
  })

  model({
    # 1. Formulation-specific absorption rate constant. Oral solution is
    #    the reference (FORM_TABLET = 0 and FORM_VPA_SR = 0).
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Apparent clearance, Equation (1):
    #    CL/F = 0.646 * (TBW/70)^0.75 * (1 + PHT effect)^PHT
    #           * (1 + PB effect)^PB * (1 + CBZ effect)^CBZ * (AGE/15)^-0.0154
    cl <- exp(lcl + e_wt_cl * log(WT / 70) + e_age_cl * log(AGE / 15) + etalcl) *
      (1 + e_conmed_pht_cl * CONMED_PHT) *
      (1 + e_conmed_pb_cl * CONMED_PB) *
      (1 + e_conmed_cbz_cl * CONMED_CBZ)

    # 3. Apparent volume of distribution, Equation (2): V/F = 14 * (TBW/70)^1
    vc <- exp(lvc + e_wt_vc * log(WT / 70))

    # 4. Micro-constant
    kel <- cl / vc

    # 5. One-compartment ODE system with first-order oral absorption
    #    (NONMEM ADVAN2 TRANS2)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 6. Observation (total serum valproic acid, mg/L) and residual error
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
