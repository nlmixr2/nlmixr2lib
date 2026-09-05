Tang_2023_tenecteplase <- function() {
  description <- paste(
    "Two-compartment linear population PK model for intravenous bolus",
    "tenecteplase in adults with acute myocardial infarction (Tang 2023;",
    "phase II TIMI 10B, 785 plasma concentrations from 103 patients).",
    "Clearance is the sum of a renal-function-independent intercept and an",
    "additive linear effect of Cockcroft-Gault creatinine clearance",
    "normalized to 70 kg body weight (CRCL-NM), centered on the cohort",
    "median 82.6 mL/min/70 kg; allometric body-weight scaling is fixed at",
    "0.75 on CL and Q and 1 on Vc and Vp with a 70 kg reference. Because",
    "the sandwich ELISA cannot distinguish tenecteplase from endogenous",
    "tissue plasminogen activator, the predicted observation is the sum of",
    "the model-predicted tenecteplase concentration and a per-subject",
    "background signal (bl_tpa, typical 18.6 ng/mL). Inter-individual",
    "variability on Vc is perfectly correlated with that on CL and is",
    "constructed as vc_eta_scale * etalcl; the log-additive residual error",
    "itself carries inter-individual variability on its magnitude. The",
    "authors caution that the statistically selected renal contribution to",
    "clearance is biologically implausible for a 65 kDa protein cleared",
    "mainly by hepatic metabolism, and that it may reflect collinearity",
    "between CRCL-NM, age and serum creatinine in the stepwise covariate",
    "search."
  )
  reference <- paste(
    "Tang F, Langenhorst J, Dang S, Kassir N, Owen R, Purdon B, Magnusson",
    "MO, Deng R (2023). Population Pharmacokinetics of Tenecteplase in",
    "Patients With Acute Myocardial Infarction and Application to Patients",
    "With Acute Ischemic Stroke. The Journal of Clinical Pharmacology",
    "63(2):197-209. doi:10.1002/jcph.2164.",
    sep = " "
  )
  vignette <- "Tang_2023_tenecteplase"
  units    <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "tenecteplase", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "tenecteplase", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives allometric scaling of CL and Q (exponent 0.750) and of Vc",
        "and Vp (exponent 1.00) on a 70 kg reference weight, per Eq. 6",
        "'CovEff_WT = (WT/70)^theta_WT' and Eq. 9. Both exponents are fixed",
        "rather than estimated: Results 'Model Development' states 'Fixing",
        "the exponents for the allometric scaling to 1.0 and 0.75 for Vc",
        "and Vp and for CL and Q, respectively, did not significantly",
        "worsen the model fit (run 13)', and Table 1 reports both exponent",
        "rows with no RSE. Weight was carried as a mechanistic covariate in",
        "the base model rather than entering the stepwise search: Methods",
        "'Covariate Model' -- 'WT (on CL, Vc, Vp, and Q) was considered as",
        "a mechanistic covariate and was included in the base model.'",
        "Cohort weight 84.1 +/- 20.0 kg (Table S1, all doses); the",
        "reference patient used for the paper's own simulations weighs",
        "81.8 kg, the cohort median (Results 'Evaluation of Covariate",
        "Effect by Simulation')."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault, normalized to 70 kg body weight",
      units              = "mL/min/70 kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CRCL-NM. Methods 'Covariate Model': 'creatinine",
        "clearance (estimated using the Cockroft-Gault equation and",
        "normalized to 70 kg [CRCL-NM])'. NOTE the normalization basis:",
        "this paper normalizes to 70 kg BODY WEIGHT, not to 1.73 m^2 BSA,",
        "so values are in mL/min/70 kg and are NOT interchangeable with a",
        "BSA-normalized eGFR. The canonical CRCL column already carries",
        "non-BSA-normalized Cockcroft-Gault variants (see the CRCL entry in",
        "inst/references/covariate-columns.md and the",
        "Chen_2023_nemonoxacin.R / Delattre_2010_amikacin.R /",
        "Georges_2009_ceftazidime.R precedents); the per-model",
        "normalization basis is recorded here as that entry requires.",
        "The effect on CL is additive linear, centered on the cohort median",
        "82.6 mL/min/70 kg, per Eq. 9: 'CL_total = (CL_population + theta x",
        "(CRCL-NM_i - 82.6)) x (WT_i/70)^0.75'. Cohort CRCL-NM 83.3 +/-",
        "25.5 mL/min/70 kg (Table S1, all doses). Renal-function strata",
        "used in the paper's forest plots (Methods 'Evaluation of Covariate",
        "Effect by Simulation'): moderate to severe impairment < 60, mild",
        "impairment 60-90, normal 90-120, supranormal >= 120 mL/min/70 kg.",
        "The authors explicitly caution against a mechanistic reading of",
        "this covariate -- Discussion: 'given the biological unlikelihood of",
        "a renal contribution to CL, statistical significance of CRCL-NM",
        "needs to be interpreted with caution.'"
      ),
      source_name        = "CRCL-NM"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  # Documentation only: checkModelConventions() does not require these to be
  # referenced in model(). Age and race WERE selected on CL in the paper's
  # alternative covariate analysis (run 50, Table S4), which is reported as a
  # sensitivity analysis rather than the final model and is therefore not
  # encoded here; see the vignette's Assumptions and deviations section.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = paste(
        "Cohort 56.4 +/- 11.4 years (Table S1); median 56 years (Table S4",
        "footnote a). Tested on CL and Vc as a prespecified structural",
        "covariate and not selected in the final model's SCM. In the",
        "alternative covariate analysis that removed CRCL-NM and CREA from",
        "the structural set, age WAS selected on CL with an effect of",
        "-0.00539 per year centered at 56 years (run 50, Table S4), but",
        "that run has a worse OFV than the final model (-441.92 vs",
        "-449.46). Results 'Model Development' notes age, CREA and CRCL-NM",
        "were strongly correlated (Figure S2)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator", units = "unitless", type = "categorical",
      notes = paste(
        "75 of 103 patients (72.8%) were male, i.e. 27.2% female (Table",
        "S1). Tested on CL and Vc as a prespecified structural covariate;",
        "not selected."
      )
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator", units = "unitless", type = "categorical",
      notes = paste(
        "14 of 103 patients (13.6%) (Table S1). Tested on CL and Vc; not",
        "selected in the final model. Selected in the alternative covariate",
        "analysis with an effect of -0.130 on CL relative to the reference",
        "category (run 50, Table S4). Results 'Model Development': the",
        "'other' and 'Asian' subcategories had limited patient numbers and",
        "were pooled into the 'White' reference category for the SCM, so",
        "the alternative model's reference is White + Asian + other (Table",
        "S4 footnote b)."
      )
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator", units = "unitless", type = "categorical",
      notes = paste(
        "14 of 103 patients (13.6%) (Table S1). Not selected in the final",
        "model; effect -0.0542 on CL in the alternative covariate analysis",
        "(run 50, Table S4, RSE 58.6%)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator", units = "unitless", type = "categorical",
      notes = "1 of 103 patients (1.0%) (Table S1); pooled into the White reference category for the SCM."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Cohort 38.9 +/- 42.1 U/L (Table S1). Tested on CL as a prespecified structural covariate; not selected."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = paste(
        "Cohort 53.7 +/- 33.5 U/L (Table S1). Tested on CL; not selected.",
        "Uniquely among the covariates, ALT was entered as the per-patient",
        "MEAN from baseline to discharge rather than the baseline value,",
        "because 'a large proportion of patients had missing baseline",
        "values, and there appeared to be no relevant consistent trends in",
        "ALT values at different sampling times' (Methods 'Covariate",
        "Model')."
      )
    ),
    CREAT = list(
      description = "Serum creatinine", units = "mg/dL", type = "continuous",
      notes = paste(
        "Cohort 1.08 +/- 0.956 mg/dL (Table S1). Tested on CL; not",
        "selected, and removed from the structural covariate set entirely",
        "in the alternative covariate analysis. Strongly correlated with",
        "age and CRCL-NM (Figure S2)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 103L,
    n_studies      = 1L,
    age_median     = "56 years (mean 56.4 +/- 11.4)",
    weight_median  = "81.8 kg (mean 84.1 +/- 20.0)",
    sex_female_pct = 100 * (103 - 75) / 103,
    race_ethnicity = c(White = 100 * 72 / 103, Black = 100 * 14 / 103,
                       Hispanic = 100 * 14 / 103, Other = 100 * 2 / 103,
                       Asian = 100 * 1 / 103),
    disease_state  = paste(
      "Adults with acute myocardial infarction enrolled in the phase II",
      "TIMI 10B study. All patients received 150-325 mg of oral or",
      "intravenous aspirin daily, and heparin before or as soon as",
      "possible after study drug. Concomitant medications recorded from 48",
      "hours before treatment through hospital discharge were near",
      "universal for anticoagulants (100%), platelet aggregation",
      "inhibitors (100%) and antihypertensives (99.0%), and common for",
      "benzodiazepines (89.3%), acetaminophen (75.7%), antihistamines",
      "(50.5%) and H2 inhibitors (45.6%) (Table S1)."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance normalized to 70 kg (CRCL-NM)",
      "83.3 +/- 25.5 mL/min/70 kg overall (Table S1); the reference",
      "patient used in the paper's simulations has the cohort median 82.6",
      "mL/min/70 kg. Serum creatinine 1.08 +/- 0.956 mg/dL. The cohort",
      "spans moderate-to-severe impairment (< 60) through supranormal",
      "(>= 120) mL/min/70 kg; Figure S2 footnote notes four patients with",
      "extreme creatinine values."
    ),
    dose_range     = paste(
      "Single intravenous bolus tenecteplase 30 mg (n = 52), 40 mg",
      "(n = 31) or 50 mg (n = 20). The 50 mg dose was suspended on 22",
      "August 1996 and replaced by 40 mg as a conservative safety measure",
      "after 3 intracranial hemorrhages among 78 patients dosed at 50 mg.",
      "The dose approved for acute myocardial infarction is 45 mg at the",
      "reference weight of 81.8 kg; the dose anticipated for acute",
      "ischemic stroke is 0.25 mg/kg."
    ),
    notes          = paste(
      "785 PK observations from 103 tenecteplase-treated patients in TIMI",
      "10B. Samples at baseline and 2, 30, 60, 90, 120, 180 and 360",
      "minutes after the start of study drug administration; all samples",
      "were above the 8 ng/mL lower limit of quantification. The sandwich",
      "ELISA may measure bound as well as free tenecteplase and cannot",
      "distinguish tenecteplase from endogenous tissue plasminogen",
      "activator, which is why the model carries an additive per-subject",
      "background signal. NONMEM 7.4.4 with FOCE+interaction;",
      "PsN 4.9.0; prediction-corrected VPC with 500 simulations. Final",
      "model is run 34 (OFV -449.46, condition number 38.97). Baseline",
      "demographics from Table S1; final parameter estimates from Table 1",
      "and Table S4. The model was partially externally validated against",
      "summary statistics from 75 patients with acute ischemic stroke in",
      "study N1811s (doses 0.1, 0.2, 0.4 and 0.5 mg/kg), where it",
      "overpredicted the observed mean concentrations by a factor of 1.39",
      "while 71.7% of the derived observation distribution fell inside the",
      "model's 90% prediction interval."
    )
  )

  ini({
    # Structural parameters -- Tang 2023 Table 1 (final model, run 34),
    # reproduced identically in Table S4 'Final model' column. All typical
    # values are for a patient of WT = 70 kg and CRCL-NM = 82.6 mL/min/70 kg.
    #
    # Final covariate model, Eq. 9 verbatim:
    #   CL_total = (CL_population + theta x (CRCL-NM_i - 82.6)) x (WT_i/70)^0.75
    # with the volumes and Q scaled by Eq. 6 'CovEff_WT = (WT/70)^theta_WT'.
    # Eq. 1 gives exponential IIV 'theta_pi = theta_p x exp(eta_pi)'.
    #
    # CL intercept: clearance at the cohort-median CRCL-NM of 82.6 mL/min/70 kg
    # for a 70 kg patient. Eq. 8 shows the same model written as the sum of a
    # renal and a nonrenal clearance; the nonrenal component implied here is
    # 91.3 - 0.322 x 82.6 = 64.7 mL/min. This intercept-plus-slope encoding of a
    # linear renal-function effect on CL follows Chen_2023_nemonoxacin.R.
    lcl <- log(91.3);  label("CL at WT = 70 kg and CRCL-NM = 82.6 mL/min/70 kg (mL/min)") # Table 1 row 'CL, mL/min' = 91.3 (RSE 2.52%); Table S4 final-model column
    lvc <- log(3610);  label("Central volume of distribution Vc at WT = 70 kg (mL)")      # Table 1 row 'Vc, mL' = 3610 (RSE 2.94%)
    lq  <- log(10.0);  label("Intercompartmental clearance Q at WT = 70 kg (mL/min)")     # Table 1 row 'Q, mL/min' = 10.0 (RSE 5.89%)
    lvp <- log(914);   label("Peripheral volume of distribution Vp at WT = 70 kg (mL)")   # Table 1 row 'Vp, mL' = 914 (RSE 5.70%)

    # Allometric exponents, fixed rather than estimated. Results 'Model
    # Development': 'Fixing the exponents for the allometric scaling to 1.0 and
    # 0.75 for Vc and Vp and for CL and Q, respectively, did not significantly
    # worsen the model fit (run 13).' Table 1 reports both rows with no RSE.
    e_wt_cl_q  <- fixed(0.750); label("Allometric exponent of body weight on CL and Q (unitless)")   # Table 1 row 'Exponent for WT on CL and Q' = 0.750, no RSE
    e_wt_vc_vp <- fixed(1.00);  label("Allometric exponent of body weight on Vc and Vp (unitless)")  # Table 1 row 'Exponent for WT on volume' = 1.00, no RSE

    # Additive linear CRCL-NM effect on CL, centered on the cohort median 82.6
    # mL/min/70 kg (Eq. 9). Interpreted by the authors as the renal clearance of
    # tenecteplase relative to that of creatinine: Results 'Model Development' --
    # 'the CRCL-NM covariate effect (theta) in Equation 8 was estimated to be
    # 0.322, suggesting that the renal CL of tenecteplase was 0.322 compared with
    # the renal CL of creatinine, which was higher than expected.' The implied
    # renal fraction of total CL at the median covariate values is
    # 0.322 x 82.6 / 91.3 = 29.1%, matching the paper's stated 29.2%.
    e_crcl_cl <- 0.322; label("Slope of the additive linear CRCL-NM effect on CL (mL/min per mL/min/70 kg)") # Table 1 row 'CRCL-NM effect on CL' = 0.322 (RSE 10.6%); Eq. 9

    # Background signal. The assay cannot distinguish tenecteplase from
    # endogenous tissue plasminogen activator, so Eq. 2 predicts the observation
    # as 'y_ij = C_ij + BS_i' -- the model-predicted tenecteplase concentration
    # plus a per-subject background estimated from the baseline samples. Encoded
    # on the log scale because Eq. 1 gives it exponential IIV.
    lbl_tpa <- log(18.6); label("Baseline endogenous t-PA assay background signal BS (ng/mL)") # Table 1 row 'BS, ng/mL' = 18.6 (RSE 5.48%); Eq. 2

    # IIV on Vc is perfectly correlated with IIV on CL. Results 'Model
    # Development': 'The correlation between IIV terms of Vc and CL was
    # significant but resulted in a high condition number (>1000) and was
    # estimated close to 1.0/100% (run 4); fixing the correlation to 100% (e.g.,
    # a single interindividual variability term (ETA) was sampled for CL and Vc
    # with a population-wide scalar from individual CL to individual Vc)
    # stabilized the model fit (run 7).' Encoded as eta_Vc = vc_eta_scale x
    # etalcl, the Hirt_2009_efavirenz.R / Prytula_2016_tacrolimus.R
    # shared-eta-scaler pattern. Implied omega for Vc is 0.913 x 0.199 = 0.182.
    vc_eta_scale <- 0.913; label("Scaling factor relating eta_Vc to eta_CL (correlation fixed to 1; eta_Vc = vc_eta_scale * etalcl)") # Table 1 row 'IIV CL - IIV Vc scaler' = 0.913 (RSE 5.87%)

    # Inter-individual variability. Table 1 reports these in a column headed
    # 'CV', but the values are omega standard deviations on the exponential-IIV
    # scale of Eq. 1, not variances and not sqrt(exp(omega^2)-1) coefficients of
    # variation. Two independent signals fix the scale: (a) the Table S3 and
    # Table S4 captions state 'The RSE for IIV and RUV parameters is reported on
    # the approximate SD scale'; (b) the last row is the log-additive residual
    # error, whose sigma on the log scale IS its approximate CV, so the column
    # cannot be holding variances. Methods Eq. 1 defines eta_pi as 'a normally
    # distributed random variable with mean 0 and standard deviation (SD)
    # omega_p'. Recorded here as omega^2 = SD^2.
    etalcl     ~ 0.199^2  # Table 1 row 'IIV CL, CV'  = 0.199 (RSE 11.0%, shrinkage 12.0%) -> variance 0.0396
    etalbl_tpa ~ 0.285^2  # Table 1 row 'IIV BS, CV'  = 0.285 (RSE 18.5%, shrinkage 25.5%) -> variance 0.0812
    etaexpSd   ~ 0.789^2  # Table 1 row 'IIV RUV, CV' = 0.789 (RSE 10.7%, shrinkage 0%)    -> variance 0.6225

    # Residual error: log-additive (Eq. 3), i.e. exponential / log-normal on the
    # linear scale, which nlmixr2 writes as `~ lnorm(expSd)`. Eq. 3 verbatim:
    # 'log y_ij = log yhat_ij + eps_add,ij x exp(eta_pi)' -- note the residual
    # magnitude itself carries exponential IIV (etaexpSd above), which Results
    # 'Model Development' introduces at run 14. Methods: 'Implementing IIV on the
    # residual model reduces the influence of outliers with less need to exclude
    # individual observations.'
    expSd <- 0.311; label("Log-additive residual error SD, typical subject (log-scale SD)") # Table 1 row 'LogAdd RUV tenecteplase, CV' = 0.311 (RSE 7.74%, shrinkage 0.907%); Eq. 3
  })

  model({
    # 1 mg/mL = 1e6 ng/mL. Doses enter in mg and volumes are in mL, so the
    # amount/volume ratio is mg/mL and must be converted to the ng/mL reporting
    # units of the assay and of Table 1's BS row.
    mgPerMlToNgPerMl <- 1e6

    # Allometric scaling on a 70 kg reference (Eq. 6). The same exponent is
    # shared by CL and Q, and by Vc and Vp, exactly as the paper constrains it:
    # Methods 'Base Model' -- 'the WT covariate effect was identical between
    # central volume (Vc) and peripheral volume (Vp) and between clearance (CL)
    # and intercompartmental clearance (Q)'.
    wtCl  <- (WT / 70)^e_wt_cl_q
    wtVol <- (WT / 70)^e_wt_vc_vp

    # Eq. 9: CL_total = (CL_population + theta x (CRCL-NM - 82.6)) x (WT/70)^0.75,
    # then exponential IIV per Eq. 1. The intercept stays positive across the
    # whole physiological CRCL-NM range -- at CRCL-NM = 0 it equals the Eq. 8
    # nonrenal clearance of 91.3 - 0.322 x 82.6 = 64.7 mL/min -- so no guard
    # against a negative clearance is required.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 82.6)) * wtCl * exp(etalcl)

    # Single shared eta: eta_Vc = vc_eta_scale x etalcl (correlation fixed to 1).
    vc <- exp(lvc) * wtVol * exp(vc_eta_scale * etalcl)

    # No IIV was supported on Q or Vp. Results 'Model Development': 'The
    # inclusion of IIV terms for Q and Vp was not supported by the data (runs 2
    # and 3, respectively).'
    q  <- exp(lq) * wtCl
    vp <- exp(lvp) * wtVol

    # Per-subject assay background signal (Eq. 2), exponential IIV per Eq. 1.
    bl_tpa <- exp(lbl_tpa + etalbl_tpa)

    # Per-subject residual error magnitude (Eq. 3): eps is scaled by exp(eta),
    # so the individual log-scale residual SD is expSd x exp(etaexpSd).
    expSdi <- expSd * exp(etaexpSd)

    d/dt(central)     <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-                        (q / vc) * central - (q / vp) * peripheral1

    # Cc is what the assay reports: model-predicted tenecteplase plus the
    # endogenous t-PA background (Eq. 2). This folds the baseline into Cc
    # following the Bauer_2023_vonicogAlfa.R / Visscher_2025_parathyroidHormone.R
    # convention. Set bl_tpa's typical value to 0 to recover the
    # tenecteplase-only concentration.
    Cc <- central / vc * mgPerMlToNgPerMl + bl_tpa

    Cc ~ lnorm(expSdi)
  })
}
