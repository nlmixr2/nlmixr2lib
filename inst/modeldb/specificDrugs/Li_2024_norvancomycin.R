Li_2024_norvancomycin <- function() {
  description <- "Two-compartment intravenous population PK model for norvancomycin (demethylvancomycin) in Chinese adults hospitalised with community-acquired pneumonia caused by gram-positive cocci, developed from prospectively collected peak and trough serum concentrations at a single centre in Shijiazhuang. Clearance carries two covariates -- a median-centered power function of age and a median-centered exponential function of serum creatinine; the volumes and the intercompartmental clearance carry none. The model underpins the paper's Monte Carlo dosing recommendations against the AUC24h/MIC >= 361 PK/PD breakpoint."
  reference   <- "Li Y, Jiao X, Sun G, Wang F, Wu X, Dong W, Lu W, Zhang Z, Yuan Y, Zhang Z. Population Pharmacokinetics and Dosing Optimization of Norvancomycin for Chinese Patients with Community-Acquired Pneumonia. Infect Drug Resist. 2024;17:5881-5893. doi:10.2147/IDR.S496776"
  vignette    <- "Li_2024_norvancomycin"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the source: norvancomycin was given as
  # a 1-h intravenous infusion and quantified in SERUM by HPLC-UV at 236 nm,
  # calibration range 0.25-100 ug/mL, LOD 0.05 ug/mL (Methods, "Assay of Serum
  # NVCM and Creatinine"). The paper reports concentrations in ug/mL, which is
  # numerically identical to the mg/L used here and in its own AUC units
  # (mg.h/L). The peripheral compartment is a mathematical disposition
  # compartment with no named tissue in the source; it is recorded as serum
  # because the model was fitted only to serum data.
  compartmentData <- list(
    central     = list(analyte = "norvancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "norvancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Age at enrolment",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as the median-centered power term (AGE / 57.5)^-0.426. The centering constant 57.5 years is the",
        "cohort MEDIAN age of Table 1 (mean 54.91, SD 15.66, median 57.50, range 27-80), not the mean; the equation printed",
        "in Results, 'PPK Modeling', divides by 57.5 explicitly. Retained in the final model after backward elimination",
        "together with Scr (joint dOFV = -23.95 for age and Scr on CL). Cohort range 27-80 years, so the covariate is",
        "uninformative outside that band; the paper's dosing tables (Table 3) stratify age only as 18-64 and 65-80 years.",
        sep = " "
      ),
      source_name        = "Age"
    ),
    CREAT = list(
      description        = "Serum creatinine, measured by an enzymatic method (Creatinine plus ver.2, Roche Diagnostics)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNITS ARE LOAD-BEARING: this column is umol/L, not mg/dL. The effect is an EXPONENTIAL-LINEAR term centered on the",
        "cohort median, exp(-0.00886 * (CREAT - 59)), so the coefficient -0.00886 has units of 1/(umol/L). Supplying mg/dL",
        "would silently flatten the covariate almost completely -- 1 mg/dL = 88.4 umol/L, so a mg/dL-scaled column spans",
        "about 0.3-1.5 where the fitted term expects 26-130, and the whole renal-function effect would collapse to a",
        "1-2% modulation. The centering constant 59 umol/L is the cohort MEDIAN of Table 1 (mean 64.09, SD 21.77,",
        "median 59.00, range 26-130) and is printed explicitly in the Results equation. Cohort range 26-130 umol/L;",
        "patients on renal replacement therapy were excluded and the Discussion states that extrapolation to renal",
        "insufficiency or failure is not supported (the paper's own Table 3 nonetheless tabulates a 133-<178 umol/L",
        "stratum, which is an extrapolation beyond the observed maximum of 130).",
        sep = " "
      ),
      source_name        = "Scr"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Entered the model at forward inclusion together with age and Scr (dOFV < -3.84, p < 0.05) but was removed at",
        "backward elimination; only age and Scr on CL survived (Results, 'PPK Modeling'). No coefficient is published.",
        "Cohort mean 165.94 cm, SD 6.76, median 165, range 153-185 (Table 1).",
        sep = " "
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate and not selected at forward inclusion; the final model carries NO body-size term",
        "at all, so the parameters below are whole-body values for a cohort of mean weight 64.75 kg (SD 11.24, median 64,",
        "range 46-90; Table 1) and must not be allometrically rescaled. Cohort BMI mean 23.50 kg/m2 (range 16.07-30.42).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = NULL,
      type        = "categorical",
      notes       = "Screened as a candidate covariate and not selected (Results, 'PPK Modeling'). Cohort 17 male / 17 female of 34."
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate and NOT selected -- the raw serum creatinine gave the better fit, which the",
        "Discussion flags as a departure from the earlier norvancomycin and vancomycin literature where CLcr is the usual",
        "renal covariate. Cohort mean 107.57 mL/min, SD 41.06, median 102.93, range 35.61-195.17 (Table 1, footnote a).",
        "Not BSA-normalized.",
        sep = " "
      )
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate and not selected. Cohort mean 99.97 mg/L, SD 77.02, range 1.59-235.20 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate covariate and not selected. Cohort mean 28.43 U/L, median 20.40, range 3.40-88.00 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate covariate and not selected. Cohort mean 28.52 U/L, median 21.55, range 6.80-87.90 (Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate covariate and not selected. Cohort mean 90.09 U/L, median 81.50, range 50.00-187.00 (Table 1)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate covariate and not selected. Cohort mean 12.47 umol/L, median 11.05, range 2.63-73.30 (Table 1)."
    ),
    BILI_DIRECT = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function candidate covariate and not selected. Cohort mean 7.15 umol/L, median 5.00, range 1.17-68.70 (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_observations = 231L,
    n_studies      = 1L,
    age_range      = "27-80 years; mean 54.91, SD 15.66, median 57.50 (Table 1)",
    age_median     = "57.5 years",
    weight_range   = "46-90 kg; mean 64.75, SD 11.24, median 64.00 (Table 1)",
    weight_median  = "64 kg",
    height_range   = "153-185 cm; mean 165.94, SD 6.76, median 165.00 (Table 1)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported beyond nationality; single-country Chinese cohort recruited in Shijiazhuang, Hebei.",
    disease_state  = "Hospitalised adults (>= 18 years) with community-acquired pneumonia, with suspected or confirmed pulmonary infection attributable to gram-positive bacteria. Recruited from Respiratory and Critical Care Medicine (23), Neurosurgery (6) and Haematology (5). Pregnant patients and patients on renal replacement therapy were excluded. Ten Staphylococcus isolates were recovered from 34 patients (4 methicillin-sensitive S. aureus, 4 MRSA, 2 S. epidermidis), all susceptible to norvancomycin and vancomycin with MIC <= 1 mg/L.",
    renal_function = "Serum creatinine mean 64.09 umol/L, SD 21.77, median 59.00, range 26-130 (Table 1); Cockcroft-Gault creatinine clearance mean 107.57 mL/min, SD 41.06, median 102.93, range 35.61-195.17. The cohort is essentially renally intact -- the observed maximum serum creatinine of 130 umol/L is below the 133 umol/L lower edge of the paper's own second dosing stratum, and the Discussion states that the number of patients with mild-to-moderate renal insufficiency was too small to support extrapolation.",
    hepatic_function = "ALT mean 28.43 U/L, AST mean 28.52 U/L, ALP mean 90.09 U/L, total bilirubin mean 12.47 umol/L, direct bilirubin mean 7.15 umol/L (Table 1). No hepatic marker was retained as a covariate.",
    dose_range     = "800 mg intravenously every 12 h, each dose infused over 1 h by syringe pump (Methods, 'Dosage Regimen and Sampling'). This single regimen generated all of the modelled data; the 200-800 mg q8h/q12h/q24h grid in Table 3 is Monte Carlo simulation output, not observed dosing.",
    regions        = "China (The Second Hospital of Hebei Medical University, Shijiazhuang, Hebei)",
    notes          = paste(
      "Prospective, single-centre, open-label observational study run between November 2020 and March 2022; Chinese Clinical",
      "Trial Registry ChiCTR2000039794; ethics approval 2020EC07-05-2. 231 serum norvancomycin concentrations from 34",
      "patients (3-8 samples each), of which 115 (49.8%) were peaks and 116 (50.2%) troughs. Troughs were drawn 0.5 h",
      "before a dose and peaks 0.5-1.5 h after the end of the 1-h infusion, repeated on the second, third, fourth and last",
      "day of therapy. Norvancomycin was quantified by HPLC-UV (Waters e2695 / 2489, Diamonsil C18 150 x 4.6 mm 5 um, 236",
      "nm), calibration range 0.25-100 ug/mL, LOD 0.05 ug/mL, vancomycin as internal standard. Fitted in NONMEM 7.5.0 by",
      "FOCE-I with interaction; evaluated by 1000-sample bootstrap in PsN 5.2.6, VPC, and NPDE. Baseline demographics are",
      "Table 1 and the parameter estimates Table 2 of Li 2024.",
      sep = " "
    )
  )

  ini({
    # All point estimates are the "Final model estimates" column of Li 2024
    # Table 2, and every one of them is reproduced in the four typical-value
    # equations printed in Results, "PPK Modeling" (recovered from the PDF with
    # `pdftotext -layout`; the docling-derived markdown renders that display
    # block as `<!-- formula-not-decoded -->`):
    #
    #   CL (L/h) = 3.15 x (Age/57.5)^-0.425 x e^((Scr-59) x (-0.00886)) x e^etaCL
    #   V1 (L)   = 12.3 x e^etaV1
    #   V2 (L)   = 115  x e^etaV2
    #   Q (L/h)  = 5.21 x e^etaQ
    #
    # Interindividual variability was exponential (Methods, "Basic Model":
    # Pi = TV(P) x e^eta_i), so the typical values are entered on the log scale.

    lcl <- log(3.15); label("Clearance at the reference age of 57.5 y and serum creatinine of 59 umol/L (L/h)")  # Li 2024 Table 2, CL = 3.15 L/h (RSE 8.2%; bootstrap median 3.13, 95% CI 2.69-3.58)
    lvc <- log(12.3); label("Central volume of distribution (L)")                                                # Li 2024 Table 2, V1 = 12.3 L (RSE 7.4%; bootstrap median 12.3, 95% CI 10.3-14.5)
    lvp <- log(115);  label("Peripheral volume of distribution (L)")                                             # Li 2024 Table 2, V2 = 115 L (RSE 18.6%; bootstrap median 113, 95% CI 74.3-180.5)
    lq  <- log(5.21); label("Intercompartmental clearance (L/h)")                                                # Li 2024 Table 2, Q = 5.21 L/h (RSE 7.2%; bootstrap median 5.17, 95% CI 4.17-6.03)

    # Covariate effects on clearance. Both are centered on the Table 1 cohort
    # MEDIAN of the covariate, and both centering constants are printed inside
    # the Results equation rather than being inferred: Age/57.5 against a Table
    # 1 median age of 57.50 y, and (Scr - 59) against a Table 1 median serum
    # creatinine of 59.00 umol/L. The two are therefore self-consistent and no
    # back-solving was needed.
    #
    # The age term is a POWER of the centered ratio; the creatinine term is an
    # EXPONENTIAL of the centered DIFFERENCE, so e_creat_cl is a slope in
    # 1/(umol/L) and NOT a power exponent (contrast Blackman_2026_methotrexate,
    # where the identically-named parameter is an exponent). See
    # covariateData$CREAT for why feeding this model a mg/dL creatinine column
    # would silently destroy the covariate effect.
    #
    # DISCREPANCY, 3rd significant figure of the age exponent: the Results
    # equation prints -0.425 while Table 2 prints -0.426 in BOTH the final-model
    # column AND the bootstrap-median column. The two independent Table 2 prints
    # are taken as authoritative and -0.426 is used; the difference moves typical
    # CL by under 0.03% anywhere in the observed 27-80 y range, so nothing
    # downstream depends on the choice.
    e_age_cl   <- -0.426;   label("Power exponent on (AGE / 57.5 y) for CL (unitless)")                                # Li 2024 Table 2, Age effect on CL = -0.426 (RSE 39.9%; bootstrap median -0.426, 95% CI -0.765 to -0.003)
    e_creat_cl <- -0.00886; label("Exponential slope on (CREAT - 59 umol/L) for CL (per umol/L)")                      # Li 2024 Table 2, Scr effect on CL = -0.00886 (RSE 22.2%; bootstrap median -0.00901, 95% CI -0.016 to -0.004)

    # Interindividual variability. The Table 2 abbreviation footnote states
    # "IIV, interindividual variability (variance value)", so the tabulated
    # numbers are omega^2 and are entered here unchanged as variances.
    etalcl ~ 0.0803  # Li 2024 Table 2, IIV of CL = 0.0803 (variance; RSE 17.1%; bootstrap median 0.0671, 95% CI 0.033-0.115) -> 28.9% CV
    etalvc ~ 0.0424  # Li 2024 Table 2, IIV of V1 = 0.0424 (variance; RSE 46.4%; bootstrap median 0.0381, 95% CI 0.005-0.106) -> 20.8% CV
    etalvp ~ 0.996   # Li 2024 Table 2, IIV of V2 = 0.996  (variance; RSE 22.6%; bootstrap median 0.898, 95% CI 0.271-1.88)  -> 128% CV
    etalq  ~ 0.0968  # Li 2024 Table 2, IIV of Q  = 0.0968 (variance; RSE 25.2%; bootstrap median 0.0952, 95% CI 0.019-0.239) -> 31.9% CV

    # Residual error. REPORTING GAP, encoded per the standing policy for an
    # unreported residual component. Results, "PPK Modeling" states that "the
    # residual variability was fitted to a combined error model", and Methods,
    # "Basic Model" defines that model as Cobs = Cpred x (1 + eps1) + eps2. But
    # Table 2 tabulates only one residual row, "RV (proportional) = 0.0401",
    # and the additive term eps2 is never given a value anywhere in the paper.
    # The proportional term is therefore taken from Table 2 and the additive
    # term is fixed at zero rather than invented; the practical consequence is
    # that this model under-disperses the very lowest concentrations, which for
    # a trough-and-peak dataset with an LOD of 0.05 ug/mL and observed troughs
    # around 10 mg/L is a small effect. See the vignette Errata.
    #
    # The Table 2 footnote states "RV, residual variability (variance value)",
    # so 0.0401 is sigma^2 and propSd is its square root, sqrt(0.0401) = 0.2002.
    propSd <- 0.2002498; label("Proportional residual error (fraction)")                     # Li 2024 Table 2, RV (proportional) = 0.0401 (variance; RSE 9.5%; bootstrap median 0.0390, 95% CI 0.028-0.053) -> SD 0.2002
    addSd  <- fixed(0);  label("Additive residual SD (mg/L; 0 -- not reported in the source)")  # Li 2024: Methods declare a combined error model but Table 2 publishes no additive term
  })

  model({
    # Covariate model on clearance, exactly as printed in Li 2024 Results,
    # "PPK Modeling":
    #
    #   CL_i = 3.15 * (AGE / 57.5)^-0.426 * exp(-0.00886 * (CREAT - 59)) * exp(eta_CL)
    #
    # At the reference covariate values (AGE = 57.5 y, CREAT = 59 umol/L) both
    # covariate terms collapse to 1 and typical CL is exactly the tabulated
    # 3.15 L/h. That reference point is corroborated arithmetically by the
    # paper's own dosing analysis, which defines AUCss,24h = daily dose / CL:
    # the observed 800 mg q12h regimen gives 1600 / 3.15 = 507.9 mg.h/L against
    # the reported cohort mean AUCss,24h of 505.63 mg.h/L (Results, "Dosing
    # Optimization"). The vignette runs this as a quantitative gate.
    cl <- exp(lcl + etalcl) * (AGE / 57.5)^e_age_cl * exp(e_creat_cl * (CREAT - 59))
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Two-compartment disposition with first-order elimination. Results, "PPK
    # Modeling": the data "conformed better to a two-compartment model
    # (OFV = 513.504) with first-order elimination, as opposed to a
    # one-compartment model (OFV = 585.283)".
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Doses enter the central compartment directly as intravenous infusions.
    # Every modelled dose was infused over 1 h by syringe pump (Methods,
    # "Dosage Regimen and Sampling"); the infusion duration is a property of
    # the event table (rate / dur), not of the model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Serum norvancomycin concentration in mg/L (doses in mg, volumes in L).
    # The source reports concentrations as ug/mL, which is numerically the
    # same unit, and its AUC targets as mg.h/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
