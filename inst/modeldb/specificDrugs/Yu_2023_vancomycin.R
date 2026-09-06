Yu_2023_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in 71 critically ill adults receiving intermittent intravenous vancomycin during continuous renal replacement therapy (Yu 2023). Total clearance is a power function of 24-hour residual diuresis, CL = 1.05 * 1.90^(log10(URINE_VOL_24H + 10) / 2.3) L/h, so an anuric subject retains only the CRRT-mediated clearance (1.39 L/h) while a subject producing 3 L/day reaches 2.77 L/h. Central volume is 69.0 L with no interindividual variability, because only trough concentrations were available. Age, sex, body weight, BMI, daily dose, serum creatinine, blood urea nitrogen, and CRRT modality were screened but not retained."
  reference <- "Yu Z, Liu J, Yu H, Zhou L, Zhu J, Liang G, Yang Y, Zheng Y, Han Y, Xu J, Han G, Yu L, Zhao Y. Population pharmacokinetics and individualized dosing of vancomycin for critically ill patients receiving continuous renal replacement therapy: the role of residual diuresis. Front Pharmacol. 2023;14:1298397. doi:10.3389/fphar.2023.1298397"
  vignette <- "Yu_2023_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    URINE_VOL_24H = list(
      description        = "Residual diuresis: total urine volume produced over 24 hours, measured on the day of vancomycin concentration monitoring",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2023 Methods ('24-h urine volume (UV)') and Table 1: median 160 mL/24h (IQR 7.00-780, range 0.00-6220); 33/71 (46.5%) at or below 100 mL, 12/71 (16.9%) 100-500 mL, 6/71 (8.45%) at or above 2500 mL. The only covariate retained in the final model (Supplementary Table S1 forward-inclusion / backward-elimination sequence: adding UV on CL drops OFV by 16.14, p < 0.001; removing it raises OFV by 16.14). Enters CL as the multiplicative power term e_urine_vol_24h_cl^(log10(URINE_VOL_24H + 10) / 2.3), printed in the Table 2 footnote as 'CL (L/h)=1.05*1.90^(LOG(UV+10)/2.3)'. The +10 mL offset keeps the term finite for the anuric subjects (UV = 0) that make up nearly half the cohort. The divisor 2.3 is the cohort median-normalising constant log10(160 + 10) = 2.230 rounded to two significant figures, consistent with Yu 2023 Methods ('the effects of continuous covariates were modeled using a median-normalized model'), so LOG in the footnote is base-10 and not the NONMEM natural log. The natural-log reading is falsified by the paper's own Figure 3: it would put CL at 9.79 L/h for a 3000 mL/24h subject, giving 1500 mg q12h a probability of AUC >= 400 mg*h/L of about 22%, whereas Figure 3A places 1500 mg q12h at 100% across every simulated urine volume. See the vignette section 'Which logarithm?' for the full adjudication.",
      source_name        = "UV"
    )
  )

  # Screened during covariate model building but NOT retained in the final
  # model, so they are documentation only and are not referenced in model().
  # Sources: Yu 2023 Methods (covariate screening list), Results, and
  # Supplementary Table S1 (key process of covariate screening).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = "Screened as a demographic covariate; not retained. Yu 2023 Table 1: mean 61.6 (SD 14.6) years, range 24-87."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on V (Supplementary Table S1 step 8: dOFV -0.023, p > 0.05); not retained. Yu 2023 Table 1: 27/71 female (38.0%), 44/71 male (62.0%)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on V as a power term (Supplementary Table S1 step 6: dOFV -1.435, p > 0.05); not retained. Yu 2023 Table 1: mean 64.2 (SD 14.7) kg, range 32.1-100. Note that this model does NOT scale either CL or V with body weight."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on CL as a linear term (Supplementary Table S1 step 4: dOFV -9.4, p < 0.01) and on V as a power term (step 7: dOFV -2.679, p > 0.05); the CL effect entered on forward inclusion but was removed on backward elimination (step 10: dOFV +9.405, p > 0.001). Yu 2023 Table 1: mean 23.3 (SD 3.97) kg/m^2, range 12.7-33.8."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened on CL as a power term (Supplementary Table S1 step 3: dOFV -6.136, p < 0.05); entered on forward inclusion, removed on backward elimination (step 11: dOFV +6.136, p > 0.001). Yu 2023 Table 1: median 104 umol/L (IQR 62.0-202, range 5.97-661). The Discussion notes that serum creatinine may not reflect renal function in patients with AKI receiving CRRT, which is the paper's motivation for using residual diuresis instead."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as a renal-function covariate; never reached the model. Yu 2023 Table 1: median 7.80 mmol/L (IQR 5.20-13.4, range 0.93-45.2)."
    ),
    DOSE_VANCOMYCIN_MGD = list(
      description = "Total daily vancomycin dose",
      units       = "mg/day",
      type        = "continuous",
      notes       = "Screened as a dosing-regimen covariate; never reached the model. Yu 2023 Table 1: median 1000 mg/day (IQR 500-2000, range 500-3000), equal to a median 15.4 mg/kg/day (range 5.26-50.0)."
    ),
    RRT_CVVHDF_STATUS = list(
      description = "CRRT modality indicator: 1 = continuous venovenous hemodiafiltration (CVVHDF), 0 = continuous venovenous hemofiltration (CVVH)",
      units       = "binary",
      type        = "categorical",
      notes       = "Screened on CL as a power term (Supplementary Table S1 step 5: dOFV -4.982, p < 0.05); entered on forward inclusion, removed on backward elimination (step 9: dOFV +4.977, p > 0.001). Yu 2023 coded CVVH as 1 and CVVHDF as 2; recorded here as a 0/1 indicator with CVVH as the reference. Table 1: CVVH 40/71 (56.3%), CVVHDF 31/71 (43.7%). Documentation only, so this name is NOT registered in inst/references/covariate-columns.md; the RRT_CRRT_STATUS entry there anticipates modality-specific siblings for a future model that actually retains one. Every subject in this cohort was on CRRT, so an on/off RRT_CRRT_STATUS column would be constant at 1 and carries no information here."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 71L,
    n_studies        = 1L,
    n_centers        = 4L,
    n_concentrations = 113L,
    age_range        = "24-87 years",
    age_mean         = "61.6 years (SD 14.6)",
    weight_range     = "32.1-100 kg",
    weight_mean      = "64.2 kg (SD 14.7)",
    sex_female_pct   = 38.0,
    race_ethnicity   = "Not reported (four-center Chinese ICU cohort)",
    disease_state    = "Critically ill ICU adults receiving intermittent intravenous vancomycin while on continuous renal replacement therapy. BMI mean 23.3 kg/m^2 (SD 3.97, range 12.7-33.8). Infection site: pulmonary 31/71 (43.7%), bloodstream 18/71 (25.4%), intra-abdominal 17/71 (23.9%), skin and soft tissue 3/71 (4.23%), other or undefined 17/71 (23.9%). Pathogen: undefined 35/71 (49.3%), other 16/71 (22.5%), Enterococcus 10/71 (14.1%), coagulase-negative staphylococci 7/71 (9.86%), MRSA 6/71 (8.45%), MSSA 2/71 (2.82%).",
    renal_function   = "All subjects on CRRT: CVVH 40/71 (56.3%), CVVHDF 31/71 (43.7%). 24-h urine volume median 160 mL (IQR 7.00-780, range 0.00-6220); 33/71 (46.5%) at or below 100 mL, 12/71 (16.9%) 100-500 mL, 6/71 (8.45%) at or above 2500 mL. Serum creatinine median 104 umol/L (IQR 62.0-202, range 5.97-661). Blood urea nitrogen median 7.80 mmol/L (IQR 5.20-13.4, range 0.93-45.2). CRRT settings (effluent / blood / dialysate flow rates) were not collected, which the paper lists as a study limitation.",
    dose_range       = "Intermittent intravenous vancomycin; total daily dose median 1000 mg/day (IQR 500-2000, range 500-3000), equal to a median 15.4 mg/kg/day (range 5.26-50.0).",
    regions          = "China (Sir Run Run Shaw Hospital and the Second Affiliated Hospital, Zhejiang University School of Medicine; Affiliated Xiaoshan Hospital, Hangzhou Normal University; Zhejiang Zhoushan Hospital). Retrospective, January 2019 to October 2022.",
    notes            = "Demographics from Yu 2023 Table 1. 191 trough concentrations from 101 patients were measured; after exclusions (age < 18, missing data, non-continuous CRRT, sampling more than 48 h after the last dose) 113 trough concentrations from 71 patients entered the model. Vancomycin measured by LC-MS/MS at all four centers. NONMEM 7.5.0 with PDxPop 5.3.1, FOCE with interaction. Model selection compared one- and two-compartment structures (AIC 800.588 vs 803.682, favouring one compartment). Bootstrap: 1000 replicates, 100% success rate (Table 2). Sampling was almost entirely trough-only, which is why interindividual variability on V could not be estimated."
  )

  ini({
    # Structural parameters. Yu 2023 Table 2 "Final model / Estimate [RSE (%)]"
    # column, and the Table 2 footnote which prints the final model verbatim:
    #   CL (L/h) = 1.05 * 1.90^(LOG(UV + 10) / 2.3)
    #   V  (L)   = 69.0
    #
    # lcl is therefore NOT the typical clearance of a typical subject: it is
    # the coefficient of the covariate power term, i.e. the clearance that
    # would apply if the power term equalled 1 (UV + 10 = 1 mL, outside the
    # observed range). The typical anuric subject (UV = 0) has
    #   CL = 1.05 * 1.90^(log10(10) / 2.3) = 1.39 L/h,
    # and the typical subject at the cohort median UV = 160 mL/24h has
    #   CL = 1.05 * 1.90^(log10(170) / 2.3) = 1.96 L/h.
    lcl <- log(1.05); label("Clearance coefficient of the urine-volume power term (L/h)")  # Yu 2023 Table 2: CL = 1.05 (RSE 17.2%); bootstrap median 1.07, 95% CI 0.721-1.53
    lvc <- log(69.0); label("Central volume of distribution (L)")                          # Yu 2023 Table 2: V  = 69.0 (RSE 6.61%); bootstrap median 68.6, 95% CI 59.9-77.7

    # Covariate effect: 24-hour residual diuresis on CL, entered as the BASE of
    # a power whose exponent is the median-normalised log10 urine volume
    # (Yu 2023 Table 2 footnote). See covariateData$URINE_VOL_24H$notes and the
    # vignette section "Which logarithm?" for why LOG is read as base-10.
    e_urine_vol_24h_cl <- 1.90; label("Base of the (URINE_VOL_24H + 10) power term on CL (unitless)")  # Yu 2023 Table 2: theta_UV-CL = 1.90 (RSE 16.5%); bootstrap median 1.91, 95% CI 1.31-2.71

    # Interindividual variability. Yu 2023 Table 2 row "omega (%) for CL": 12.1
    # (RSE 47.9%), bootstrap median 11.5, 95% CI 2.45-22.0. The Table 2
    # abbreviation footnote defines omega as the "interindividual VARIANCE for
    # CL", so 12.1% is omega^2 = 0.121 (omega = 0.348, i.e. 35.9% CV), not an
    # omega of 0.121. Two independent checks in the paper confirm the variance
    # reading, both worked through in the vignette:
    #   - Figure 2I (the ETA histogram) spans roughly -0.9 to +0.5 with a peak
    #     density near 1.5; an omega of 0.121 would confine it to about +/-0.36
    #     with a peak density near 3.3.
    #   - Regressing qnorm(PTA) from the Figure 3A 500 mg q12h curve on the
    #     model-implied log AUC recovers a slope of 1/omega = 2.91, i.e.
    #     omega = 0.343, within 1.3% of sqrt(0.121) = 0.3479.
    # No IIV was estimable on V ("Intre-individual variability for V was not
    # estimated", Table 2 footnote), so V carries no eta.
    etalcl ~ 0.121  # Yu 2023 Table 2, row omega (%) for CL = 12.1, i.e. omega^2 = 0.121 (omega = 0.348, 35.9% CV)

    # Proportional residual error. Yu 2023 Table 2 row "sigma (%)": 9.78
    # (RSE 29.1%), bootstrap median 9.60, 95% CI 4.70-15.2. Reported on the
    # same variance-as-percent scale as omega, so sigma^2 = 0.0978 and the
    # residual SD is sqrt(0.0978) = 0.3127 (31.3% proportional error). Yu 2023
    # Results: "An exponential model and a proportional model were used to
    # describe the interindividual variability and residual variability".
    propSd <- sqrt(0.0978); label("Proportional residual error (fraction)")  # Yu 2023 Table 2: sigma^2 = 0.0978 (reported as "sigma (%) = 9.78")
  })

  model({
    # Total vancomycin clearance: CRRT-mediated clearance plus the residual
    # native renal contribution, both absorbed into one power term in the
    # 24-hour urine volume (Yu 2023 Table 2 footnote, Discussion: "when the
    # 24-h urine volume was close to 0, the CL was the intrinsic clearance of
    # CRRT. When the 24-h urine volume increased, residual renal function
    # contributed to the total clearance").
    #
    # URINE_VOL_24H is in mL/24h; the +10 offset is inside the log so that
    # anuric subjects (UV = 0) give a finite exponent of log10(10)/2.3.
    cl <- exp(lcl + etalcl) * e_urine_vol_24h_cl^(log10(URINE_VOL_24H + 10) / 2.3)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
