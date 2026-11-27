Gomes_2017_gentamicin_endocarditis <- function() {
  description <- "One-compartment intravenous population PK model for gentamicin in adults with infective endocarditis (the 'endocarditis model' of Gomes 2017), fitted by iterative two-stage Bayesian analysis in MWPharm 3.82 for therapeutic drug monitoring. Clearance is the sum of a metabolic arm scaled linearly by total body weight (fixed at the MWPharm ICU-model value) and a renal arm proportional to Cockcroft-Gault creatinine clearance; volume is a per-kg coefficient on the MWPharm fat-corrected lean body mass. Log-normal between-subject variability on the renal fraction and the volume coefficient; residual error is the fixed EMIT assay-error polynomial."
  reference <- paste(
    "Gomes A, van der Wijk L, Proost JH, Sinha B, Touw DJ.",
    "Pharmacokinetic modeling of gentamicin in treatment of infective",
    "endocarditis: Model development and validation of existing models.",
    "PLoS One. 2017;12(5):e0177324. doi:10.1371/journal.pone.0177324.",
    "PMCID: PMC5419648.",
    sep = " "
  )
  vignette <- "Gomes_2017_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters twice. (1) Linearly on the metabolic clearance arm as",
        "cl_nonren = CLm * (WT / 70) (Methods, 'CL = CLm*(BW/70) + fr*CLcr');",
        "70 kg is the unit denominator of the MWPharm CLm (L/h/70kgBW), not",
        "a cohort statistic. (2) Through the fat-corrected lean body mass",
        "LBMc = LBM + 0.4 * (WT - LBM) that scales the volume (see HT).",
        "Modeling cohort mean 76.2 kg, range 46-121 kg (Table 1)."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Input to the lean-body-mass formula. The paper cites Chennavasin",
        "and Brater (reference 27) for LBM and Schwartz et al. (reference 28)",
        "for the fat correction but prints neither formula; both were",
        "back-solved from the per-patient LBM and LBMc columns of the S1 Data",
        "workbook (sheet 'The Hague', all 65 modeling-cohort patients) and",
        "reproduce every tabulated value to its 0.1 kg rounding:",
        "LBM = min(WT, 50 + 2.3 * (HT / 2.54 - 60)) for men and",
        "LBM = min(WT, 45.5 + 2.3 * (HT / 2.54 - 60)) for women (the Devine",
        "ideal body weight with height in inches, capped at total body",
        "weight), then LBMc = LBM + 0.4 * (WT - LBM). Modeling cohort mean",
        "173.9 cm, range 149-193 cm (Table 1)."
      ),
      source_name = "Height"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the intercept of the lean-body-mass formula only (50 kg for",
        "men, 45.5 kg for women). Sex has no separate effect on clearance or",
        "on the volume coefficient (its effect on creatinine clearance is",
        "already inside the supplied CRCL). S1 Data codes sex as 'M' / 'V'",
        "(Dutch 'vrouw' = female). Modeling cohort 21 F / 44 M (Table 1)."
      ),
      source_name = "Gender"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance by Cockcroft-Gault in ABSOLUTE mL/min, NOT",
        "BSA-normalized. Stored under the canonical CRCL column, which",
        "accepts raw mL/min when the source applies no BSA normalization."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The renal arm is fr * CLcr with CLcr in L/h; model() converts the",
        "mL/min column by 0.06 (60 min/h / 1000 mL/L). The absolute (not",
        "/1.73 m^2) scale is fixed by the paper's own Figure 1: its ICU-model",
        "typical curve for the 'standard patient' (male, 70 y, 76.2 kg,",
        "174 cm, CLcr 58 mL/min/1.73 m^2) decays at 0.156 1/h, which the",
        "ICU-model parameters reproduce only with the de-normalized",
        "58 * BSA / 1.73 = 64.0 mL/min (Du Bois BSA 1.91 m^2); the",
        "normalized 58 mL/min would give 0.142 1/h. The S1 Data workbook",
        "carries both CLcr (mL/min) and CLcr_n (/1.73 m^2) per patient.",
        "Modeling cohort mean 64.3 mL/min, range 8.7-157.5 mL/min (Table 1).",
        "Dialysis patients were excluded, and fewer than 10% of the cohort had",
        "CLcr below 15 mL/min (Discussion)."
      ),
      source_name = "CLcr"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 65,
    n_studies = 1,
    n_observations = 221,
    age_range = "32-92 years (mean 69.3)",
    weight_range = "46-121 kg (mean 76.2)",
    height_range = "149-193 cm (mean 173.9)",
    sex_female_pct = 32.3,
    race_ethnicity = "Not reported.",
    disease_state = "Infective endocarditis in adults (>= 18 years) treated with intravenous gentamicin; dialysis patients excluded",
    renal_function = "Cockcroft-Gault creatinine clearance mean 64.3 mL/min (range 8.7-157.5); mean 58.4 mL/min/1.73 m^2 (range 7.8-153.2)",
    dose_range = "Intravenous gentamicin, standard 3 mg/kg once daily as a 30-minute infusion (Figures 1 and 2); individual dosing records were not published",
    regions = "Netherlands (modeling cohort: The Hague Hospitals, 2011-2013)",
    notes = paste(
      "Retrospective observational two-cohort study of therapeutic-drug-",
      "monitoring data. Modeling cohort: 65 patients / 221 serum samples from",
      "The Hague Hospitals. A separate validation cohort of 14 patients /",
      "30 samples from the University Medical Center Groningen (2012-2014)",
      "was used only to compare the endocarditis, ICU and standard models",
      "(Table 4) and did not contribute to these estimates. Separate",
      "endocarditis episodes of one patient were entered as separate",
      "patients. Serum gentamicin by EMIT on an Architect analyser",
      "(lower limit of quantification 0.2 mg/L) for the modeling cohort."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Table 3, row 'Endocarditis model' (rounded) and the S2 Data workbook,
    # sheet 'Table 2 and 3', row 'Fixed/bayesian/bayesian' under 'Start
    # model ICU' (unrounded): the final model of Table 2 (AIC 1195,
    # SumWSS/df 5.25). Parameters are reported as a population mean and an
    # interindividual SD on the NATURAL scale; inter-individual variability
    # was assumed log-normal (Methods, Phase 1), so the mean is carried as
    # the typical value and the SD is converted to a log-scale variance by
    # omega^2 = log(1 + (SD / mean)^2).
    # ---------------------------------------------------------------------

    # Metabolic (non-renal) clearance per 70 kg total body weight. FIXED at
    # the ICU-model value because the cohort held too few patients with
    # poor renal function to estimate it (Discussion); Table 3 prints its
    # SD as 0, so it carries no between-subject variability.
    lcl_nonren <- fixed(log(0.2765))
    label("Metabolic clearance CLm at 70 kg total body weight (L/h/70 kg)") # Table 3 'Endocarditis model' CLm mean 0.277, SD 0; S2 Data 'Fixed/bayesian/bayesian' CLm 0.2765

    # fr, renal gentamicin clearance as a fraction of creatinine clearance:
    # the slope of the renal clearance arm on CLcr (both in L/h). Carried on
    # the log scale so its log-normal between-subject variability can be
    # attached.
    le_crcl_cl_renal <- log(0.6977)
    label("Renal gentamicin clearance as a fraction fr of creatinine clearance (unitless)") # Table 3 'Endocarditis model' fr mean 0.698 (95% CI 0.610-0.794); S2 Data fr 0.6977

    # Volume coefficient per kg of fat-corrected lean body mass; the
    # absolute volume vc = vd * LBMc is derived in model().
    lvd <- log(0.3115)
    label("Volume of distribution per kg fat-corrected lean body mass (L/kg LBMc)") # Table 3 'Endocarditis model' Vd mean 0.312 (95% CI 0.292-0.331) L/kgLBMc; S2 Data V 0.3115

    # MWPharm fat-distribution fraction in LBMc = LBM + fd * (WT - LBM).
    # Not printed in the paper; back-solved from the S1 Data LBM and LBMc
    # columns, where (LBMc - LBM) / (WT - LBM) has median 0.4000 across all
    # patients with WT > LBM (range 0.385-0.419, i.e. the 0.1 kg rounding
    # of the tabulated values). It is the 40% excess-weight correction of
    # Schwartz et al., the paper's reference 28.
    e_wt_vd <- fixed(0.4)
    label("Fat-distribution fraction fd weighting (WT - LBM) inside LBMc (unitless)") # back-solved from S1 Data sheet 'The Hague' LBM / LBMc / Weight columns; Methods cite reference 28 for the fat correction

    # IIV, log-scale variances from the S2 Data natural-scale SDs:
    #   fr SD 0.3578 / 0.6977 -> CV 51.28 % -> log(1 + 0.51283^2) = 0.23348
    #   Vd SD 0.0764 / 0.3115 -> CV 24.53 % -> log(1 + 0.24526^2) = 0.05841
    # The fr-Vd correlation is not in Table 3; the S2 Data bootstrap sheet
    # ('CI table 3 by Bootstrap', column 'r', 1000 resamples) gives mean
    # 0.310 (95% CI 0.072-0.511), carried as the correlation of the two
    # log-scale effects: covariance 0.310 * sqrt(0.23348 * 0.05841).
    etale_crcl_cl_renal + etalvd ~ c(
      0.23348,
      0.036204,
      0.058414
    ) # Table 3 'Endocarditis model' SDs fr 0.358 and Vd 0.076; S2 Data SDs 0.3578 and 0.0764; correlation = S2 bootstrap mean r 0.310

    # Residual error: the EMIT assay-error polynomial used for modeling
    # (Methods, Gentamicin assay): SD = 0.0766 + 0.0006 * C + 0.0064 * C^2
    # mg/L. Assay characteristics, not fitted, hence fixed (the S2 Data
    # bootstrap sheet holds S0, S1, S2 constant at these values).
    addSd <- fixed(0.0766)
    label("Assay SD polynomial intercept S0 (mg/L)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
    propSd <- fixed(0.0006)
    label("Assay SD polynomial linear coefficient S1 (fraction)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
    quadSd <- fixed(0.0064)
    label("Assay SD polynomial quadratic coefficient S2 (L/mg)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
  })

  model({
    # Lean body mass (kg): Devine ideal body weight with height in inches,
    # capped at total body weight; then the MWPharm fat correction. Both
    # back-solved from the S1 Data per-patient LBM / LBMc columns.
    ibw <- (50 - 4.5 * SEXF) + 2.3 * (HT / 2.54 - 60)
    lbm <- ibw
    if (WT < ibw) lbm <- WT
    lbmc <- lbm + e_wt_vd * (WT - lbm)

    # Cockcroft-Gault creatinine clearance from mL/min to L/h.
    clcr <- 0.06 * CRCL

    # Methods: CL = CLm * (BW / 70) + fr * CLcr.
    cl_nonren <- exp(lcl_nonren) * (WT / 70)
    fr <- exp(le_crcl_cl_renal + etale_crcl_cl_renal)
    cl_renal <- fr * clcr
    cl <- cl_nonren + cl_renal

    # Volume = Vd (L/kg LBMc) * LBMc.
    vd <- exp(lvd + etalvd)
    vc <- vd * lbmc

    kel <- cl / vc

    # One compartment, first-order elimination (Methods, Phase 1).
    d/dt(central) <- -kel * central

    Cc <- central / vc
    sdCc <- addSd + propSd * Cc + quadSd * Cc^2
    Cc ~ add(sdCc)
  })
}
