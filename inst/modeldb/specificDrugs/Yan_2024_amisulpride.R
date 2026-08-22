Yan_2024_amisulpride <- function() {
  description <- "Population PK model for oral amisulpride in Chinese adult inpatients with schizophrenia (Yan 2024). One-compartment disposition with first-order absorption and first-order elimination, parameterised on the apparent (oral) scale as CL/F and V/F because the therapeutic-drug-monitoring dataset was oral-only and bioavailability was not identifiable. The absorption rate constant Ka was fixed to 0.18 1/h, carried from earlier amisulpride work in Chinese patients, because the data are almost entirely steady-state pre-dose troughs and the absorption phase could not support an estimate. Estimated creatinine clearance (Cockcroft-Gault, raw mL/min) enters apparent clearance as a power function centred on the cohort median 114.42 mL/min; it was the only covariate retained after forward inclusion and backward elimination. Inter-individual variability is exponential on CL/F only, and residual variability is proportional. The paper's primary purpose was an external evaluation of five previously published amisulpride popPK models against this cohort; all five showed unacceptable simulation-based bias, and this model was then developed on the same independent dataset and used in Monte Carlo simulations of remedial dosing after a delayed or missed dose."
  reference   <- "Yan D, Ju G, Liu X, Shao Q, Zhang Y, Wang N, Yan K. External Validation of the Population Pharmacokinetic Models of Amisulpride and Remedial Strategies for Delayed or Missed Doses. Drug Des Devel Ther. 2024;18:6345-6358. doi:10.2147/DDDT.S469149."
  vignette    <- "Yan_2024_amisulpride"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "amisulpride", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amisulpride", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance calculated with the Cockcroft-Gault equation. NOT body-surface-area normalised: the source reports raw Cockcroft-Gault mL/min (Yan 2024 Methods, 'Data Accumulation and Blood Sampling': 'the eCLcr calculated using the Cockcroft-Gault formula', citing Cockcroft & Gault 1976). The paper separately tabulates a CKD-EPI eGFR column computed without the race component, which was screened but not retained in the final model.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in the source analysis. The only covariate retained by the forward-inclusion /",
        "backward-elimination search (Yan 2024 Results, 'PopPK Model Development and Validation':",
        "dOFV = -27.9, p < 0.001, dAIC = 28). Applied to apparent clearance as a power function centred on",
        "the cohort median: CL/F = 45.1 * (eCLcr / 114.42)^0.364.",
        "The centring value is not printed in Table 5; it is fixed by the Yan 2024 Discussion sentence",
        "'the typical clearance rate of amisulpride is 45.1 L/h (average eCLcr 114.42 mL/min)' read together",
        "with Table 1, whose eCLcr row is 114.42 (23.29, 239.51) mL/min. Cohort mean 116.7 +/- 33.53 mL/min.",
        "Do not extrapolate outside the observed 23.29-239.51 mL/min range.",
        "The power-on-median form matches the two amisulpride models this one was built from (Yan 2024 Table 3:",
        "Liu 2023 CL/F = 32.6 * (eCLcr/114.3)^0.485 and Li 2023 CL/F = 60.5 * (eGFR/113.87)^0.817).",
        "Mechanistically well supported: within 24 h about 90 percent of amisulpride is excreted renally and",
        "active tubular secretion (SLC22 organic ion transporters) is the predominant elimination pathway",
        "(Yan 2024 Introduction and Discussion).",
        "The raw (non-BSA-normalised) Cockcroft-Gault form matches the CRCL register's Delattre 2010 amikacin",
        "and Shu 2024 posaconazole precedents; the register requires the assay form to be documented per model,",
        "which this note does.",
        sep = " "
      ),
      source_name        = "eCLcr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate but not retained in the final model (Yan 2024 Results, 'PopPK Model Development and Validation'; Supplementary Table S1). Enters the retained covariate indirectly, as an input to the Cockcroft-Gault eCLcr. Yan 2024 explicitly compared a model using eCLcr against one using age, sex, weight and creatinine as separate covariates and 'found no significant differences in the parameter fitting result' (Supplementary Table S2). Cohort median 32 years, range 18-67, mean 34.42 +/- 10.56 (Yan 2024 Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained (Yan 2024 Results; Supplementary Tables S1 and S2). Enters the retained covariate indirectly, as an input to the Cockcroft-Gault eCLcr. Cohort median 62 kg, range 40-109, mean 62.93 +/- 11.71 (Yan 2024 Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Yan 2024 Results; Supplementary Tables S1 and S2). Enters the retained covariate indirectly, as an input to the Cockcroft-Gault eCLcr. 211 of 361 patients were female (Yan 2024 Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Yan 2024 Results; Supplementary Tables S1 and S2). Enters the retained covariate indirectly, as an input to the Cockcroft-Gault eCLcr. The Discussion states that eCLcr predicted amisulpride clearance better than serum creatinine alone. Cohort median 60 umol/L, range 30-148.9, mean 66.34 +/- 15.7 (Yan 2024 Table 1)."
    ),
    UA = list(
      description = "Serum uric acid",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Yan 2024 Methods, 'Data Accumulation and Blood Sampling', lists uric acid among the renal-function laboratory indicators collected; Supplementary Table S1). Cohort median 297 umol/L, range 81-698.9, mean 305.03 +/- 89.86 (Yan 2024 Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Yan 2024 Methods and Supplementary Table S1). Cohort median 3.5 mmol/L, range 0.6-10.1, mean 3.72 +/- 1.37 (Yan 2024 Table 1)."
    ),
    DOSE_AMISULPRIDE_MGD = list(
      description = "Daily amisulpride dose",
      units       = "mg/day",
      type        = "continuous",
      notes       = "Daily dosage was included among the screened covariates (Yan 2024 Results, 'PopPK Model Development and Validation': 'the demographic factors, laboratory indicators, and daily dosage were considered as covariates') but was not retained, i.e. the final model is dose-linear. Cohort median 600 mg/day, Table 1 range 200-1200 mg/day, mean 555.90 +/- 192.56."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 361L,
    n_studies      = 1L,
    age_range      = "18-67 years (median 32, mean 34.42 +/- 10.56; Yan 2024 Table 1)",
    age_median     = "32 years",
    weight_range   = "40-109 kg (median 62, mean 62.93 +/- 11.71; Yan 2024 Table 1)",
    weight_median  = "62 kg",
    sex_female_pct = 58.4,
    race_ethnicity = "Chinese (single-centre cohort at the Xi'an Mental Health Center, Xi'an, Shaanxi)",
    disease_state  = "Adult inpatients diagnosed with schizophrenia, treated with oral amisulpride for at least 72 h (at least five half-lives) before therapeutic drug monitoring.",
    dose_range     = "200-1200 mg/day oral amisulpride (median 600 mg/day, mean 555.90 +/- 192.56; Yan 2024 Table 1). 302 of 361 patients received twice-daily dosing and 59 received once-daily dosing. The Yan 2024 Results text states a daily dosage range of 200-2000 mg, which disagrees with the 200-1200 mg range printed in its own Table 1; see the vignette Errata.",
    regions        = "China (Xi'an, Shaanxi), 2017 to 2021",
    n_observations = "390 steady-state serum concentrations from 361 patients. Samples were drawn at 06:00 immediately before the next dose, i.e. trough concentrations only. Observed concentrations 54.7-1955.7 ng/mL (median 471.8, mean 546.03 +/- 341.76); dose-corrected concentrations 0.17-3.26 ng/mL per mg (median 0.88, mean 0.970 +/- 0.474) (Yan 2024 Table 1).",
    renal_function = "Cockcroft-Gault eCLcr median 114.42 mL/min (range 23.29-239.51, mean 116.7 +/- 33.53); CKD-EPI eGFR (no race component) median 118.5 (range 79.9-146.5); serum creatinine median 60 umol/L (range 30-148.9) (Yan 2024 Table 1).",
    co_medication  = "Concomitant medication was recorded but could not be entered as a covariate: the group sizes were too small and 'significant differences between the groups' prevented its inclusion (Yan 2024 Results and Discussion). The Discussion notes published reports of elevated amisulpride concentrations under co-treatment with clozapine or lithium, attributed to competition in the renal clearance pathway.",
    notes          = paste(
      "Retrospective single-centre therapeutic-drug-monitoring study. Amisulpride was quantified by a validated",
      "LC-MS/MS assay (Shimadzu 8050) with a linear range of 20-2000 ng/mL and intra- and inter-day RSD within",
      "5 percent (Yan 2024 Methods).",
      "Model estimated in NONMEM 7.4.0 with FOCE-I; internal validation by goodness-of-fit plots, a 1000-sample",
      "bootstrap and prediction-corrected VPC (Yan 2024 Table 5, Figures 1 and 2).",
      "This cohort was first used as the EXTERNAL evaluation dataset for five previously published amisulpride",
      "popPK models (Yan 2024 Tables 2-4). None passed the simulation-based diagnostics (NPDE global test",
      "p < 0.001 for all five), which is what motivated developing the present model. The five evaluated models",
      "are Liu 2023 (doi:10.5414/CP204334), Reeves 2016 (doi:10.1007/s00213-016-4379-6),",
      "Glatard 2020 (doi:10.1007/s40262-019-00821-w), Huang 2021 (doi:10.2147/DDDT.S327506) and",
      "Li 2023 (doi:10.3389/fphar.2023.1215065); they are NOT extracted from this paper's summary table --",
      "see the vignette Errata for why.",
      "The observed mean concentration (546 ng/mL) sits well above the AGNP therapeutic reference range of",
      "100-320 ng/mL, which the authors argue reflects genuinely higher effective exposures in Chinese patients",
      "rather than overdosing (Yan 2024 Discussion).",
      "Baseline demographics are in Yan 2024 Table 1.",
      sep = " "
    )
  )

  ini({
    # ====================================================================
    # Structural PK parameters (Yan 2024 Table 5).
    #
    # One-compartment model on the APPARENT (oral) scale. Amisulpride was
    # given only orally and the dataset is trough-only TDM, so
    # bioavailability F is not identifiable and the reported disposition
    # parameters are CL/F and V/F. No separate bioavailability term is
    # encoded: the nominal dose enters the depot compartment and F is
    # absorbed into lcl and lvc. (Amisulpride's absolute bioavailability is
    # about 48 percent -- Yan 2024 Introduction -- but that value is not
    # part of this model.)
    # ====================================================================

    # Ka was NOT estimated. Yan 2024 Methods, 'PopPK Model Development and
    # Validation': 'The dataset predominantly consists of trough
    # concentration data following long-term administration, making it
    # challenging to accurately estimate the drug absorption process. Based
    # on articles with favorable predictive performance, taking into
    # account the racial similarity among the Chinese population, which the
    # absorption rate constant (Ka) is fixed at 0.18.' Table 5 accordingly
    # prints '0.18 (fixed)' with no RSE. The value is Huang 2021's estimate
    # (Yan 2024 Table 3, model M4: Ka = 0.18).
    lka <- fixed(log(0.18)); label("First-order absorption rate constant Ka (1/h; taken from Huang 2021 amisulpride popPK in Chinese patients)")  # Yan 2024 Table 5: ka = 0.18 (fixed); bootstrap 0.18

    lcl <- log(45.1);        label("Apparent clearance CL/F at the cohort-median eCLcr of 114.42 mL/min (L/h)")  # Yan 2024 Table 5: CL/F = 45.1 L/h, RSE 4 percent, bootstrap median 44.9 (95% CI 41.5-48.16)

    lvc <- log(466);         label("Apparent volume of distribution V/F (L)")  # Yan 2024 Table 5: V/F = 466 L, RSE 20 percent, bootstrap median 461.6 (95% CI 319.5-742.0)

    # ====================================================================
    # Covariate effect (Yan 2024 Table 5, row 'eCLcr on CL/F').
    #
    # Table 5 prints only the exponent, 0.364; the centring value comes
    # from the Discussion ('the typical clearance rate of amisulpride is
    # 45.1 L/h (average eCLcr 114.42 mL/min)') cross-read with the Table 1
    # eCLcr median of 114.42 mL/min. The functional form is the power-on-
    # median model used by both amisulpride studies this model was built
    # from (Yan 2024 Table 3, models M1 and M5).
    # ====================================================================
    e_crcl_cl <- 0.364; label("Power exponent of estimated creatinine clearance on apparent clearance (unitless; CL/F ~ (CRCL/114.42)^e_crcl_cl)")  # Yan 2024 Table 5: eCLcr on CL/F = 0.364, RSE 19 percent, bootstrap median 0.36 (95% CI 0.223-0.495)

    # ====================================================================
    # Inter-individual variability (Yan 2024 Table 5, row 'Interindividual
    # variability CL/F' = 0.043).
    #
    # SCALE. Reported as a bare decimal with no unit and no percent sign,
    # so it could in principle be omega^2 (variance, CV 20.7 percent) or
    # omega (CV 4.3 percent). It is the VARIANCE, on two independent
    # grounds:
    #
    #   (1) The Yan 2024 Methods define the IIV model as Pi = PTV *
    #       exp(eta_i) where 'eta_i was a random variable with a zero
    #       average and omega^2 variance' -- i.e. the reported IIV
    #       parameter IS omega^2.
    #   (2) Digitising Yan 2024 Figure 1A (observations vs IPRED) and 1B
    #       (observations vs PRED) and pairing records by their shared
    #       observation value gives SD(log(IPRED/PRED)) ~ 0.30-0.41
    #       (interquartile-based estimate 0.33). Because MAP-Bayes IPREDs
    #       shrink toward PRED, that is a LOWER bound on omega. It is
    #       compatible with sqrt(0.043) = 0.207 and excludes 0.043 by
    #       roughly an order of magnitude.
    #
    # A third, coarser check agrees: the Table 1 dose-corrected
    # concentration spread (mean 0.970 +/- 0.474, i.e. log-SD ~ 0.46) minus
    # the once-daily/twice-daily regimen split and the eCLcr-explained part
    # leaves ~0.26 for IIV + nothing else, again near 0.207 and nowhere
    # near 0.043.
    #
    # IIV was estimated on CL/F only; Table 5 reports no IIV on V/F.
    # ====================================================================
    etalcl ~ 0.043  # Yan 2024 Table 5: IIV on CL/F = 0.043 (variance), RSE 32 percent, bootstrap median 0.041 (95% CI 0.02-0.08); omega = 0.207, CV 20.9 percent

    # ====================================================================
    # Residual error (Yan 2024 Table 5, row 'Proportional residual' =
    # 0.314). Yan 2024 Methods tested additive, proportional and combined
    # structures and Table 5 reports a single proportional term.
    #
    # SCALE. 0.314 is the residual STANDARD DEVIATION (a 31.4 percent
    # proportional error), not a variance. This is measured directly from
    # the paper's own figures: for a proportional error, Yan 2024 Figure 1E
    # plots |IWRES| = |DV - IPRED| / (sigma * IPRED), so the ratio of the
    # relative residual read off Figure 1A to the |IWRES| read off Figure
    # 1E IS sigma. Digitising both panels and matching the nine most
    # extreme records by their IPRED gives that ratio as 0.312, 0.312,
    # 0.314, 0.313, 0.312, 0.312, 0.313, 0.313 and 0.312 -- i.e. sigma =
    # 0.313, reproducing the printed 0.314 to three digits and excluding
    # the variance reading (sqrt(0.314) = 0.56) outright.
    #
    # Note that Table 5 therefore mixes conventions: the IIV row is a
    # variance and the residual row is an SD. The Methods wording ('sigma
    # variance') would suggest otherwise, but the figures are unambiguous
    # and are the model that was actually fitted. See the vignette Errata.
    # ====================================================================
    propSd <- 0.314; label("Proportional residual error (fraction)")  # Yan 2024 Table 5: proportional residual = 0.314, RSE 6 percent, bootstrap median 0.314 (95% CI 0.261-0.346)
  })

  model({
    # ====================================================================
    # 1. Individual PK parameters.
    #
    #     CL/F = 45.1 * (eCLcr / 114.42)^0.364 * exp(eta_CL)
    #     V/F  = 466
    #
    # centred on the Table 1 cohort-median eCLcr (see ini()).
    # ====================================================================
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 114.42)^e_crcl_cl
    vc <- exp(lvc)

    # ====================================================================
    # 2. Micro-constants
    # ====================================================================
    kel <- cl / vc

    # ====================================================================
    # 3. ODE system -- one compartment with first-order absorption.
    # Oral dose records target cmt = depot with the nominal dose in mg;
    # bioavailability is absorbed into the apparent parameters CL/F and
    # V/F, so no f(depot) term is applied.
    # ====================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ====================================================================
    # 4. Observation and residual error.
    # central is in mg and vc in L, so central/vc is mg/L == ug/mL; the
    # factor 1000 converts to the ng/mL used throughout Yan 2024.
    # ====================================================================
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
