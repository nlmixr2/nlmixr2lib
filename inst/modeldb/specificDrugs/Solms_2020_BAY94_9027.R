Solms_2020_BAY94_9027 <- function() {
  description <- "One-compartment population PK model for BAY 94-9027 (damoctocog alfa pegol, Jivi, an extended-half-life site-specifically PEGylated B-domain-deleted recombinant factor VIII) in 198 male patients aged 2-62 years with severe haemophilia A pooled from the BAY 94-9027 phase I (NCT01184820), PROTECT VIII (NCT01580293), and PROTECT VIII Kids (NCT01775618) trials (Solms 2020). Final chromogenic-assay model has lean body weight (LBW) as a power-form covariate on CL and Vc and von Willebrand factor antigen (VWF) as a power-form covariate on CL; between-subject variability is a BLOCK(2) on CL and Vc with correlation 0.449; residual error is combined additive plus proportional. NONMEM M3 likelihood was used for samples below the chromogenic-assay lower limit of quantitation (1.5-3 IU/dL)."
  reference <- "Solms A, Iorio A, Ahsman MJ, Vis P, Shah A, Berntorp E, Garmann D. Favorable Pharmacokinetic Characteristics of Extended-Half-Life Recombinant Factor VIII BAY 94-9027 Enable Robust Individual Profiling Using a Population Pharmacokinetic Approach. Clin Pharmacokinet. 2020 May;59(5):605-616. doi:10.1007/s40262-019-00832-7. PMID:31749076."
  vignette <- "Solms_2020_BAY94_9027"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight (canonical column LBM; source paper uses LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference LBW = 49.1 kg, per the Solms 2020 Table 2 footnotes b and c (CL = 1.09 * (LBW/49.1)^0.707 * (VWF/110)^-0.604; Vc = 26.2 * (LBW/49.1)^0.887). Solms 2020 Table 1 reports the cohort median LBW as 49.4 kg (range 10-75 kg); the 49.1 kg used in the centering footnote is the value the authors used in the fitted model and is preserved here exactly. Body-composition formula is not specified in the paper; in haemophilia popPK literature LBW is most commonly computed via the Hume (1966) or James (1976) formula. Stored under canonical LBM (lean body mass); LBW and LBM refer to the same quantity (Garmann_2017_BAY81_8973 follows the same convention). Body weight, height, BMI, and age were screened during covariate analysis but only LBW was retained in the final model (LBW gave the largest OFV decrease, -509.92, of all size-related covariates).",
      source_name        = "LBW"
    ),
    VWF = list(
      description        = "Plasma von Willebrand factor (VWF) antigen concentration; FVIII-protective carrier protein.",
      units              = "IU/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL only with reference 110 IU/dL (Solms 2020 Table 2 footnote b: CL = 1.09 * (LBW/49.1)^0.707 * (VWF/110)^-0.604). The centering value 110 IU/dL matches the cohort median (Solms 2020 Table 1, n=145; range 47-366 IU/dL). The exponent is negative: higher VWF protects FVIII from clearance. VWF was measured at baseline in the phase I and PROTECT VIII studies; no VWF measurements were available in PROTECT VIII Kids, and the CL-VWF relationship is informed by adult and adolescent data only. The authors note (Discussion) that the adult/adolescent VWF range (47-366 IU/dL) covers the paediatric range typically observed (Yee et al. ref 17), so the relationship is expected to extend across the modelled age range. The paper does not characterise within-subject VWF time course; simulations should treat VWF as a time-fixed baseline covariate. VWF assay was antigen-based (VWF:Ag); record per-model in covariateData[[VWF]]$notes when reusing this covariate for other FVIII products."
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened for effects on CL and Vc in the univariate analysis; LBW had the largest OFV decrease among the size-related covariates (BW, height, BMI, LBW) and was retained in preference. After LBW was added to CL and Vc, no visible residual effect of weight remained on either parameter (Solms 2020 Results, base-model section)."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a size-related covariate; LBW preferred. No residual effect of height after LBW inclusion (Solms 2020 Results, base-model section)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a size-related covariate; LBW preferred. No residual effect of BMI after LBW inclusion (Solms 2020 Results, base-model section)."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Univariate analysis identified age as eligible for the multivariate step alongside LBW and VWF, but in the stepwise deletion process the best covariate model retained only LBW on CL and Vc and VWF on CL. After LBW and VWF were added, no visible effect of age remained on either parameter (Solms 2020 Results, covariate-analysis section)."
    ),
    RACE = list(
      description = "Race",
      units       = NA_character_,
      type        = "categorical",
      notes       = "Screened (white 144/198; Asian 35/198; black 9/198; Native American 1/198; not reported 9/198). 'There was no significant relationship between CL or Vc and race or geographic region (Asia vs. other)' (Solms 2020 Results, covariate-analysis section)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 198L,
    n_studies      = 3L,
    age_range      = "2-62 years (final-model n = 198 including 53 patients aged < 12 years; phase I 21-58 years, PROTECT VIII 12-62 years, PROTECT VIII Kids 2-11 years)",
    age_median     = "28.5 years (mean 28.2, SD 17.6)",
    weight_range   = "12-126 kg",
    weight_median  = "67.0 kg (mean 62.5, SD 27.3)",
    height_range   = "87-192 cm",
    height_median  = "171 cm (mean 160, SD 26.8)",
    bmi_range      = "13-42 kg/m^2",
    bmi_median     = "22.0 kg/m^2 (mean 22.7, SD 5.6)",
    lbw_range      = "10-75 kg",
    lbw_median     = "49.4 kg (mean 44.5, SD 16.1); the model-centering footnote uses 49.1 kg (Solms 2020 Table 2)",
    vwf_range      = "47-366 IU/dL (n = 145; not measured in 53 patients from PROTECT VIII Kids)",
    vwf_median     = "110 IU/dL (mean 122, SD 55.3); used as the VWF covariate-centering reference",
    sex_female_pct = 0,
    race_ethnicity = c(White = 144L, Asian = 35L, Black = 9L, NativeAmerican = 1L, NotReported = 9L),
    disease_state  = "Severe haemophilia A (FVIII activity < 1 IU/dL) in previously treated patients with no history or current evidence of FVIII inhibitors. >= 150 prior FVIII exposure days for phase I and PROTECT VIII; > 50 prior exposure days for PROTECT VIII Kids.",
    dose_range     = "Intravenous BAY 94-9027 across phase I single-dose pharmacokinetic profiling, PROTECT VIII prophylaxis (twice weekly, every-5-day, weekly schedules with 25-60 IU/kg), and PROTECT VIII Kids (paediatric dosing). Simulation analyses in the paper used single 25-60 IU/kg doses and steady-state prophylaxis with 25-60 IU/kg.",
    regions        = "Multinational (phase I, PROTECT VIII, PROTECT VIII Kids)",
    notes          = "Final model dataset: 2196 chromogenic-assay observations (455 / 21% BLQ) from 198 male patients across the three BAY 94-9027 trials. Eight patients (seven with anti-PEG antibodies and/or perceived loss of efficacy; one with drug hypersensitivity) were excluded from the final model. A separate fit to 1648 one-stage assay observations from 146 phase I and PROTECT VIII subjects gave parameter estimates similar but not identical to the chromogenic fit (Solms 2020 Table 2, one-stage assay column); only the chromogenic fit is encoded here. Haemophilia A is X-linked recessive and the BAY 94-9027 trials enrolled males; sex_female_pct = 0 reflects this. BAY 94-9027 is approved for prophylaxis and on-demand treatment in patients aged >= 12 years, so inferences from this model are limited to the >= 12-year subpopulation even though the model was fit using all 2-62 year data."
  )

  ini({
    # Structural parameters at the paper's reference covariates (LBW = 49.1 kg,
    # VWF = 110 IU/dL). Volumes in dL and clearances in dL/h match the paper's
    # units. Source: Solms 2020 Table 2 ("Chromogenic assay" column, final
    # model).
    lcl <- log(1.09); label("Clearance for the reference 49.1 kg LBW and 110 IU/dL VWF patient (CL, dL/h)") # Solms 2020 Table 2: CL = 1.09 dL/h (2.19% RSE; 95% CI 1.04-1.14)
    lvc <- log(26.2); label("Central volume of distribution for the reference 49.1 kg LBW patient (Vc, dL)") # Solms 2020 Table 2: Vc = 26.2 dL (1.18% RSE; 95% CI 25.6-26.8)

    # Covariate effects: power-form scaling for LBW on CL and Vc, and VWF on CL.
    # Form (Solms 2020 Table 2 footnotes b, c):
    #   CL = CL_pop * (LBW/49.1)^e_lbw_cl * (VWF/110)^e_vwf_cl
    #   Vc = Vc_pop * (LBW/49.1)^e_lbw_vc
    e_lbw_cl <-  0.707; label("Power exponent of LBW on CL (unitless)")  # Solms 2020 Table 2: 0.707 (7.93% RSE; 95% CI 0.597-0.817)
    e_vwf_cl <- -0.604; label("Power exponent of VWF on CL (unitless)")  # Solms 2020 Table 2: -0.604 (10.2% RSE; 95% CI -0.725 to -0.483)
    e_lbw_vc <-  0.887; label("Power exponent of LBW on Vc (unitless)")  # Solms 2020 Table 2: 0.887 (2.78% RSE; 95% CI 0.839-0.935)

    # Inter-individual variability: BLOCK(2) on CL and Vc (Solms 2020 Table 2).
    # The paper reports BSV as % CV. Following the convention used by the same
    # author group in Garmann_2017_BAY81_8973 (log-normal BSV model
    # P_j = P_pop * exp(eta_j), so the eta-scale variance is the exact
    # log-normal map omega^2 = log((CV%/100)^2 + 1)):
    #   CV(CL) = 24.0% -> omega^2 = log(1 + 0.240^2) = 0.0560090
    #   CV(Vc) = 12.8% -> omega^2 = log(1 + 0.128^2) = 0.0162519
    # The reported correlation is 0.449; corresponding covariance is
    #   cov = 0.449 * sqrt(0.0560090 * 0.0162519) = 0.0134738.
    # Inter-occasion variability (IOV) was estimated during model development
    # (IOV variance 0.0195 for CL and 0.00475 for Vc; base model section) but
    # was not retained in the final model because it had no effect on the
    # population estimates of CL or Vc; see vignette deviations.
    etalcl + etalvc ~ c(0.0560090,
                        0.0134738,  0.0162519)  # Solms 2020 Table 2: BSV CL = 24.0% CV (16.4% RSE; 95% CI 19.7-27.8); BSV Vc = 12.8% CV (17.5% RSE; 95% CI 10.4-14.9); Correlation(BSV CL, BSV Vc) = 0.449 (19.2% RSE; 95% CI 0.28-0.618)

    # Combined proportional + additive residual unexplained variability
    # (Solms 2020 Methods 2.2; Table 2). Both components are reported as the
    # estimated variance (Table 2 footnote a explicitly distinguishes the SE
    # of the variance/covariance estimate from the derived CV); the additive
    # variance is 1.78 IU^2/dL^2 (-> SD = sqrt(1.78) = 1.3342 IU/dL) and the
    # proportional variance is 0.175 (-> CV = sqrt(0.175) = 0.4183, the 41.8%
    # CV reported in Table 2 footnote d). The combined model is
    #   obs = pred * (1 + eps_prop) + eps_add.
    addSd  <- sqrt(1.78);  label("Additive residual error SD (IU/dL)")     # Solms 2020 Table 2: additive RUV variance = 1.78 (19.2% RSE; 95% CI 1.11-2.45); sqrt(1.78) = 1.3342 IU/dL
    propSd <- sqrt(0.175); label("Proportional residual error (fraction)") # Solms 2020 Table 2: proportional RUV variance = 0.175 (5.66% RSE; 95% CI 0.156-0.194); footnote d reports the corresponding CV as 41.8%, i.e. sqrt(0.175) = 0.4183
  })
  model({
    # Individual PK parameters with covariate effects (Solms 2020 Table 2
    # footnotes b, c).
    cl <- exp(lcl + etalcl) * (LBM / 49.1)^e_lbw_cl * (VWF / 110)^e_vwf_cl
    vc <- exp(lvc + etalvc) * (LBM / 49.1)^e_lbw_vc

    kel <- cl / vc

    # One-compartment IV model: BAY 94-9027 is administered as an intravenous
    # infusion; doses enter the central compartment directly. The infusion
    # duration is supplied through the user's event table (rate / dur), not
    # the structural model. A two-compartment model was tested during base-
    # model development but the peripheral compartment could not be identified
    # (Q RSE > 70%) and goodness-of-fit showed no trend, so the authors
    # retained the one-compartment model (Solms 2020 Results, base-model
    # section).
    d/dt(central) <- -kel * central

    # FVIII activity: dose in IU and central volume in dL -> central / vc has
    # units IU/dL (chromogenic FVIII activity, equivalent to % of normal FVIII).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
