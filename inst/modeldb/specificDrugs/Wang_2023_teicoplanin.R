Wang_2023_teicoplanin <- function() {
  description <- paste(
    "Two-compartment IV-infusion population PK model for teicoplanin in 108",
    "critically ill adults in a Chinese medical-surgical intensive care unit",
    "(Wang 2023). CKD-EPI estimated glomerular filtration rate is the only",
    "retained covariate and enters clearance as a LINEAR term centred on 50",
    "mL/min/1.73 m^2: CL = 0.838 * (1 + (eGFR - 50) * 0.00823) L/h.",
    "Between-subject variability is exponential on CL, Vc and Vp; the",
    "inter-compartmental clearance Q carries no random effect because its",
    "omega was fixed to zero in the published fit. Residual error is",
    "proportional. The paper's Monte Carlo simulations use the model to",
    "recommend teicoplanin loading and maintenance regimens stratified by",
    "renal function and by whether the infection is deep-seated",
    "(trough target 15 mg/L) or not (trough target 10 mg/L)."
  )
  reference <- paste(
    "Wang Y, Yao F, Chen S, Ouyang X, Lan J, Wu Z, Wang Y, Chen J, Wang X,",
    "Chen C. Optimal teicoplanin dosage regimens in critically ill patients:",
    "population pharmacokinetics and dosing simulations based on renal",
    "function and infection type.",
    "Drug Des Devel Ther. 2023;17:2259-2271. doi:10.2147/DDDT.S413662"
  )
  vignette <- "Wang_2023_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Body-surface-area-normalized glomerular filtration rate estimated with the CKD-EPI equation, the only covariate retained in the final model (linear effect on clearance)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. Creatinine-based CKD-EPI estimate (Wang 2023 Methods,",
        "'Patients and Data Collection'; the paper cites the CKD-EPI equation as its",
        "reference 20), already normalized to 1.73 m^2 body surface area, so it is",
        "stored under the canonical CRCL column per inst/references/covariate-columns.md.",
        "Enters the final model as the LINEAR centred term printed in the Wang 2023",
        "Table 2 notes: CL (L/h) = 0.838 * (1 + (eGFR - 50) * 0.0082). The generic",
        "form is given in the Supplementary Table S1 notes as",
        "P_ij = P_tv,j * [1 + theta_j * (COV - COV_ave)] * e^eta_i, where COV_ave is",
        "described as the average of the continuous covariate; the numeric centring",
        "constant 50 mL/min/1.73 m^2 is printed explicitly in the Table 2 notes and is",
        "used here. The cohort MEDIAN eGFR is lower (29.2 mL/min/1.73 m^2, range",
        "4.89-170, Wang 2023 Table 1), consistent with 50 being the mean of a",
        "right-skewed renal-function distribution rather than the median.",
        "Because the effect is linear rather than a power term, the multiplier is",
        "0.589 at eGFR = 0 and 1.658 at eGFR = 130; it stays positive across and well",
        "beyond the observed eGFR range, so no lower guard is required. The paper",
        "warns against extrapolating outside the studied cohort, in particular to",
        "highly augmented renal clearance (Wang 2023 Discussion, limitations).",
        "43 of 108 model-building subjects (39.8%) were receiving continuous renal",
        "replacement therapy; CRRT status was tested and NOT retained (Supplementary",
        "Table S1: CL-CRRT dOFV -6.685, below the -10.828 backward-elimination",
        "threshold), so eGFR is the sole renal descriptor in the model and is applied",
        "to CRRT and non-CRRT records alike."
      ),
      source_name        = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 61 kg (range 40-115) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -0.930) and on Vc (dOFV -1.448) and not retained (Supplementary Table S1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Median 33.6 g/L (range 17.5-49.2) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -3.700) and on Vc (dOFV -2.688) and not retained (Supplementary Table S1). The Discussion attributes the null result to the cohort's albumin being within the normal range, in contrast with earlier teicoplanin models in which albumin acted on the volume of distribution."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Median 161 umol/L (range 24.0-696) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -23.839, the largest single drop in the forward step) and on Vc (dOFV -7.590). Not retained: the authors state they deliberately selected only mutually independent covariates (Wang 2023 Methods, 'PK Modeling'), and eGFR - which is itself derived from serum creatinine - was carried forward instead."
    ),
    CRCL_COCKCROFT = list(
      description = "Cockcroft-Gault creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Median 34.8 mL/min (range 8.82-296) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -12.856) and on Vc (dOFV -1.745) and not retained in favour of the BSA-normalized CKD-EPI eGFR (Supplementary Table S1). Documented here only as a screened covariate; the canonical CRCL column in this model carries the CKD-EPI eGFR, not this Cockcroft-Gault value."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 56.5 years (range 21-91) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -2.369) and on Vc (dOFV -0.210) and not retained (Supplementary Table S1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "38 of 108 model-building subjects were female (Wang 2023 Table 1). Screened on CL (dOFV -1.356) and on Vc (dOFV +38.137, i.e. a worse fit) and not retained (Supplementary Table S1). The supplement's categorical-covariate form is P_ij = P_tv,j * theta_j^COV_gender * e^eta_i with COV_gender = 1 for male and 0 for female."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "43 of 108 model-building subjects (39.8%) were on CRRT (Wang 2023 Table 1). Screened on CL (dOFV -6.685) and on Vc (dOFV -1.770) and not retained (Supplementary Table S1). The Discussion notes this contrasts with reports that CRRT modality influences teicoplanin disposition, and recommends therapeutic drug monitoring for these patients; the influence of CRRT settings (blood-flow rate, dialysate flow rate, filter) is listed as a study limitation."
    ),
    APACHE2 = list(
      description = "Baseline APACHE II severity-of-illness score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Median 25.5 (range 10-37) in the model-building group (Wang 2023 Table 1). Screened on CL (dOFV -2.619) and on Vc (dOFV -4.827) and not retained (Supplementary Table S1)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 108L,
    n_studies        = 1L,
    n_concentrations = 304L,
    age_range        = "21-91 years (Wang 2023 Table 1, model-establishment group)",
    age_median       = "56.5 years",
    weight_range     = "40-115 kg (Wang 2023 Table 1, model-establishment group)",
    weight_median    = "61 kg",
    sex_female_pct   = 35.2,
    race_ethnicity   = "Not reported by category; single-centre Chinese ICU cohort (Guangdong Provincial People's Hospital, Guangzhou)",
    disease_state    = paste(
      "Critically ill adults (age >= 18 years) admitted to a medical-surgical",
      "intensive care unit and treated intravenously with teicoplanin for a",
      "diagnosed or clinically suspected Gram-positive infection. Infection sites",
      "recorded in the model-establishment group (Wang 2023 Table 1; the counts sum",
      "to more than the 108 subjects, so patients could carry more than one site):",
      "pneumonia 73, intra-abdominal infection 13, bacteraemia 10, urinary tract",
      "infection 5, infective endocarditis 4, skin and soft tissue infection 1,",
      "intracranial infection 1, unknown 20. Median baseline APACHE II score 25.5",
      "(range 10-37); median duration of teicoplanin therapy 11 days (range 3-33).",
      "The cohort was older and substantially renally impaired: median eGFR 29.2",
      "mL/min/1.73 m^2 and 39.8% receiving CRRT. Patients missing any of the",
      "required demographic or laboratory data were excluded."
    ),
    dose_range       = paste(
      "Teicoplanin (Targocid, Sanofi) by intravenous infusion over 1 h. Across the",
      "full 151-patient cohort the observed regimens were (Wang 2023 Methods,",
      "'Dosing Regimen'): group 1, 400 mg q12h x 3 loading then 400 mg q24h (n = 35);",
      "group 2, 600 mg q12h for 3-5 doses then 400-600 mg q24h (n = 17); group 3,",
      "800 mg q12h for 3-5 doses then 400-800 mg q24h (n = 8); group 4, irregular",
      "doses at irregular intervals (n = 91). The published Monte Carlo simulations",
      "(Wang 2023 Tables 3, S2-1 and S2-2) evaluated 400/600 mg q12h x 3 or x 5",
      "loading with 400/600 mg q24h maintenance for ordinary Gram-positive",
      "infections, and 800 mg q12h x 3 or x 5 loading with 400/600/800 mg q24h",
      "maintenance for deep-seated infection."
    ),
    regions          = "China (Guangdong Provincial People's Hospital ICU, Guangzhou; July 2018 to January 2020)",
    renal_function   = paste(
      "Model-establishment group (Wang 2023 Table 1): CKD-EPI eGFR median 29.2",
      "mL/min/1.73 m^2 (range 4.89-170); Cockcroft-Gault creatinine clearance median",
      "34.8 mL/min (range 8.82-296); serum creatinine median 161 umol/L (range",
      "24.0-696); 43/108 (39.8%) on continuous renal replacement therapy."
    ),
    observed_exposure = paste(
      "Model-establishment group (Wang 2023 Table 1): observed plasma Cmax median",
      "39.9 mg/L (range 8.35-120) and observed plasma Cmin median 10.5 mg/L (range",
      "0.86-62.7)."
    ),
    screened_covariates = paste(
      "Tested-but-not-retained (Wang 2023 Supplementary Table S1, dOFV against the",
      "two-compartment base model at OFV 1564.335): on CL - sex -1.356, age -2.369,",
      "body weight -0.930, albumin -3.700, serum creatinine -23.839, Cockcroft-Gault",
      "CrCL -12.856, CRRT -6.685, APACHE II -2.619; on Vc - sex +38.137, age -0.210,",
      "body weight -1.448, albumin -2.688, serum creatinine -7.590, Cockcroft-Gault",
      "CrCL -1.745, eGFR -5.782, CRRT -1.770, APACHE II -4.827. Only eGFR on CL",
      "(dOFV -20.374, final-model OFV 1543.961) survived backward elimination at",
      "dOFV < -10.828. Serum creatinine gave a larger drop than eGFR but was not",
      "carried forward because the authors restricted themselves to mutually",
      "independent covariates and eGFR is derived from serum creatinine."
    ),
    external_validation = paste(
      "An independent group of 43 patients contributing one concentration each",
      "(Wang 2023 Results). Mean prediction error -2.22 mg/L and RMSE 8.88 mg/L",
      "against the final model."
    ),
    notes            = paste(
      "Single-centre prospective observational study. 347 plasma teicoplanin",
      "concentrations from 151 ICU patients in total; the population PK model was",
      "built on the 304 samples from the 108 patients who contributed more than two",
      "samples, and the remaining 43 single-sample patients formed the external",
      "validation group. Sampling: before one loading dose, at the end of a",
      "maintenance-dose infusion, and before the next maintenance dose. Assay:",
      "HPLC-MS/MS quantifying the teicoplanin A2-2 component (about 50% of the",
      "five-component mixture and the most active) as the surrogate for teicoplanin,",
      "with daptomycin as internal standard; linear 1.0-100.0 mg/L, LLOQ 1.0 mg/L,",
      "intra- and inter-batch precision RSD <= 10.9%. Estimation: NONMEM 7.3, FOCE",
      "with interaction, run through Pirana 2.9.0. Qualification: goodness-of-fit",
      "plots, 1000-sample bootstrap (100% success rate), VPC (Supplementary Figure",
      "S1), NPDE, and the external validation above."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects. Wang 2023 Table 2, 'Final Model / Estimate'
    # column; the same four values are restated in Results ('The estimates of
    # CL, Vc, Q and Vp were 0.838 L/h, 14.4 L, 3.08 L/h and 51.6 L') and in the
    # typical-value equations printed in the Table 2 notes.
    #
    # lcl is the clearance of the REFERENCE subject at eGFR = 50 mL/min/1.73 m^2,
    # which is the centring constant printed in the Table 2 note - not the cohort
    # median eGFR of 29.2. At the cohort median the model gives
    # 0.838 * (1 + (29.2 - 50) * 0.00823) = 0.694 L/h.
    # ------------------------------------------------------------------
    lcl <- log(0.838); label("Clearance at eGFR = 50 mL/min/1.73 m^2 (L/h)")   # Wang 2023 Table 2: CL = 0.838 L/h (RSE 8.4%; bootstrap median 0.86, 95% CI 0.68-0.99)
    lvc <- log(14.4);  label("Central volume Vc (L)")                          # Wang 2023 Table 2: Vc = 14.4 L (RSE 7.2%; bootstrap median 14.3, 95% CI 10.9-16.6)
    lq  <- log(3.08);  label("Inter-compartmental clearance Q (L/h)")          # Wang 2023 Table 2: Q = 3.08 L/h (RSE 14.3%; bootstrap median 3.25, 95% CI 2.49-9.62)
    lvp <- log(51.6);  label("Peripheral volume Vp (L)")                       # Wang 2023 Table 2: Vp = 51.6 L (RSE 10.7%; bootstrap median 52.1, 95% CI 41.8-66.7)

    # ------------------------------------------------------------------
    # Covariate effect: LINEAR (not power) slope of eGFR on CL, centred at
    # 50 mL/min/1.73 m^2. Wang 2023 Table 2 note prints
    #   CL (L/h) = 0.838 * (1 + (eGFR - 50) * 0.0082)
    # with the coefficient rounded to two significant figures; the Table 2 row
    # and the Abstract both give the unrounded 0.00823, which is used here.
    # The minus sign inside (eGFR - 50) is not recoverable from a text extraction
    # of the PDF (the notes are typeset in a symbol font whose minus glyph drops
    # out); it was confirmed by rendering page 6 of the PDF to an image, and the
    # subtraction is corroborated by the generic form in the Supplementary
    # Table S1 notes, P_ij = P_tv,j * [1 + theta_j * (COV - COV_ave)] * e^eta_i.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.00823; label("Linear slope of eGFR on CL, per mL/min/1.73 m^2 above 50")  # Wang 2023 Table 2: theta_eGFR-CL = 0.00823 (RSE 28.3%; bootstrap median 0.0083, 95% CI 0.0036-0.016); Abstract 'eGFR adjustment factor of 0.00823'

    # ------------------------------------------------------------------
    # Between-subject variability. Wang 2023 Table 2 reports omega as a
    # percentage. Those percentages are omega (the SD on the log scale) x 100,
    # the usual NONMEM approximate-CV convention, NOT sqrt(exp(omega^2) - 1):
    # squaring the Vp percentage reproduces the exponent printed in the Table 2
    # note exactly (0.622^2 = 0.3869, printed as e^0.387), which pins the
    # convention. The variances below are therefore the squared percentages.
    # ------------------------------------------------------------------
    etalcl ~ 0.142   # Wang 2023 Table 2: omega_CL = 37.7% (RSE 20%)   -> omega^2 = 0.377^2 = 0.1421
    etalvc ~ 0.0708  # Wang 2023 Table 2: omega_Vc = 26.6% (RSE 43.4%) -> omega^2 = 0.266^2 = 0.0708
    etalvp ~ 0.387   # Wang 2023 Table 2: omega_Vp = 62.2% (RSE 18%)   -> omega^2 = 0.622^2 = 0.3869; the Table 2 note prints this as e^0.387

    # Wang 2023 Table 2 reports omega_Q as '0 FIX': the inter-compartmental
    # clearance carries no between-subject variability in the published fit
    # (Results also states the IIV of Q was 0%). No eta is declared on Q. It is
    # omitted rather than written as `etalq ~ fixed(0)` because a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve (same handling as Kim_2026_midazolam_postecmo.R and
    # Wattanakul_2024_primaquine_motherinfant.R). Mathematically identical: a
    # zero-variance eta contributes nothing.
    #
    # No off-diagonal covariances were published.

    # ------------------------------------------------------------------
    # Residual variability: proportional only (Wang 2023 Results, 'Exponential
    # and proportional error models were selected to describe the
    # inter-individual and residual variability, respectively').
    #
    # The magnitude is printed two ways, with the digits transposed: Table 2
    # gives 31.7% and the Results text gives 'the intra-individual variability
    # was 37.1%'. Table 2 is used, because its own bootstrap 95% CI for this
    # parameter (21.4-36.4%) contains 31.7% and excludes 37.1%.
    # ------------------------------------------------------------------
    propSd <- 0.317; label("Proportional residual error (fraction)")  # Wang 2023 Table 2: sigma proportional error = 31.7% (RSE 10.6%; bootstrap median 30.3%, 95% CI 21.4-36.4%), encoded as a fraction
  })

  model({
    # Individual PK parameters. eGFR (canonical column CRCL) enters CL only, as
    # the linear centred term printed in the Wang 2023 Table 2 notes. Vc and Vp
    # carry log-normal IIV with no covariate; Q carries neither.
    cl <- exp(lcl + etalcl) * (1 + (CRCL - 50) * e_crcl_cl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Intravenous infusion into the central compartment (1 h infusions in the
    # study; Wang 2023 Methods, 'Dosing Regimen'), first-order distribution to
    # peripheral1 and first-order elimination from central. Dose in mg and
    # volumes in L give central / vc in mg/L.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
