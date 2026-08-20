Gaspar_2025_fexofenadine <- function() {
  description <- paste0(
    "Population pharmacokinetic model for low-dose (25 mg) oral fexofenadine ",
    "used as a P-glycoprotein (Pgp) phenotyping probe in the Geneva cocktail, ",
    "fit to 1089 dried-blood-spot concentrations from 449 hospitalized older ",
    "adult polymorbid patients (OptimAT, NCT03477331; median age 71 years) ",
    "pooled with 10 healthy volunteers (Geneva cocktail study, NCT01731067). ",
    "Two-compartment model with sequential zero-order (D1 = 1.59 h) plus ",
    "first-order (ka2 = 0.282/h) absorption -- both absorption parameters ",
    "carried unchanged from Piscitelli 2023 -- and linear elimination from ",
    "the central compartment. Apparent clearance CL/F = 116 L/h carries three ",
    "log-additive covariate effects: creatinine clearance (power 0.33 on ",
    "CRCL/77.2), age (power -0.59 on AGE/71), and concomitant Pgp-inhibitor ",
    "use (exp(-0.38), a 32 percent reduction in CL/F). Interindividual ",
    "variability is exponential on CL/F only (82 percent CV); residual ",
    "variability is proportional (23 percent CV). The paper's headline output ",
    "is a grid of typical AUC(0-6 h) values stratified by CKD stage, age ",
    "decade, and Pgp-inhibitor status, intended to reinterpret the ",
    "healthy-volunteer Pgp-activity thresholds of Bosilkovska 2014 in an ",
    "older adult inpatient population."
  )
  reference <- paste(
    "Gaspar F, Jacost-Descombes C, Gosselin P, Reny JL, Guidi M, Csajka C,",
    "Samer C, Daali Y, Terrier J. (2025). Improving Understanding of",
    "Fexofenadine Pharmacokinetics to Assess Pgp Phenotypic Activity in Older",
    "Adult Patients Using Population Pharmacokinetic Modeling.",
    "Clin Pharmacokinet 64:275-284. doi:10.1007/s40262-024-01470-4.",
    "The sequential zero-order / first-order absorption parameters D1 and ka2",
    "were held at the literature values of Piscitelli J, Nikanjam M,",
    "Capparelli EV, Blaquera CL, Penzak SR, Nolin TD, et al. (2023)",
    "Ther Drug Monit 45(4):539-545 (Gaspar 2025 reference 12).",
    sep = " "
  )
  vignette <- "Gaspar_2025_fexofenadine"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    CRCL = list(
      description        = paste0(
        "Creatinine clearance estimated with the Cockcroft-Gault equation. ",
        "Enters CL/F as a power function normalized to the population median ",
        "of 77.2."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Gaspar 2025 Methods 2.3.2 states CLcr was 'calculated using the ",
        "Cockcroft-Gault equation'; Table 2 labels the column 'Estimated ",
        "creatinine clearance, mL/min/1.73 m^2' with a median of 77.2 ",
        "(range 5.7-147.3). Note that Cockcroft-Gault returns mL/min and is ",
        "not intrinsically BSA-normalized, so the paper's own unit label and ",
        "its stated estimating equation are inconsistent; the values are used ",
        "here exactly as tabulated. Time-fixed (a single baseline value per ",
        "subject). Missing values in the hospitalized cohort were imputed at ",
        "the patient-population median and missing laboratory values in the ",
        "healthy volunteers were set to healthy-volunteer typical values ",
        "(Methods 2.3.2, final paragraph); Table 2 reports 0 percent missing ",
        "for CLcr in OptimAT and no CLcr at all for the 10 volunteers. ",
        "NORMALIZATION CONFLICT: the Table 3 footnote equation divides by ",
        "'CLcr_m, the median value of CLcr in the population' (77.2), while ",
        "Methods Eq. 1 defines Cov_weight as a 'typical value' and names ",
        "'CLcr (100 mL/min)'. This model uses the Table 3 footnote (77.2) -- ",
        "see the model file comment on `e_crcl_cl` for the full rationale."
      ),
      source_name        = "CLcr"
    ),
    AGE = list(
      description        = paste0(
        "Age at inclusion. Enters CL/F as a power function normalized to the ",
        "population median of 71 years; older age lowers CL/F."
      ),
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Gaspar 2025 Table 2: median 71 years (range 25-97) in the 449 ",
        "OptimAT hospitalized patients and 23 years (range 20-36) in the 10 ",
        "Geneva cocktail healthy volunteers; 0 percent missing in both ",
        "cohorts. Time-fixed. Age was retained as an effect on CL/F ",
        "independent of renal function (Results 3.2: 'Age accounted for 20 ",
        "percent of the initial variance in fexofenadine clearance, and CLcr ",
        "contributed an additional 13 percent'). The normalizing median is ",
        "AGE_m in the Table 3 footnote equation; 71 years is the OptimAT ",
        "median from Table 2, and pooling the 10 volunteers shifts the median ",
        "of the 459 subjects negligibly."
      ),
      source_name        = "AGE"
    ),
    CONMED_PGP_INH = list(
      description        = paste0(
        "1 = subject concomitantly treated with a P-glycoprotein inhibitor, ",
        "0 = no Pgp-inhibitor coadministration. Lowers CL/F by 32 percent."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant Pgp inhibitor)",
      notes              = paste0(
        "Gaspar 2025 Table 2: 14.0 percent of the 449 OptimAT patients and 0 ",
        "percent of the 10 healthy volunteers (Results 3.1 rounds this to ",
        "'approximately 15 percent'). The Table 2 footnote a lists the pooled ",
        "inhibitor set with within-group frequencies: amiodarone (n = 24, ",
        "38.1 percent), diltiazem (n = 13, 20.6 percent), clarithromycin ",
        "(n = 7, 11.1 percent), quetiapine (n = 7, 11.1 percent), duloxetine ",
        "(n = 4, 6.3 percent), fluoxetine (n = 3, 4.8 percent), ketoconazole ",
        "(n = 2, 3.2 percent), paroxetine (n = 2, 3.2 percent), verapamil ",
        "(n = 1, 1.6 percent). Coded 'presence = 1, absence = 0' per the same ",
        "footnote. Time-fixed over the 6 h sampling window. Enters as a ",
        "log-additive shift on log(CL/F) (Table 3 footnote equation term ",
        "beta_PGP * PGP_i), i.e. a multiplicative exp(-0.38) = 0.684 factor."
      ),
      source_name        = "PGP"
    )
  )

  # Covariates Gaspar 2025 screened (Methods 2.3.2) but did NOT retain in the
  # final model. Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = paste0(
        "Screened per Methods 2.3.2 ('demographic factors (e.g., body ",
        "weight, age, gender)'). Methods Eq. 1 names 70 kg as the weight ",
        "normalizing value, but weight was not retained on any parameter: ",
        "Results 3.2 lists only age, CLcr, co-administration of Pgp ",
        "inhibitors, and population type as univariately significant, and ",
        "only the first three survived the multivariate step. Table 2 median ",
        "77.2 kg (range 37-130) in OptimAT, 70.1 kg (69-86) in volunteers."
      )
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = paste0(
        "Screened as 'gender' (Methods 2.3.2); not retained. Table 2 reports ",
        "69 percent male in OptimAT and 100 percent male among the 10 ",
        "healthy volunteers."
      )
    ),
    HT = list(
      description = "Height", units = "cm", type = "continuous",
      notes = "Table 2 demographic (median 172 cm OptimAT, 170 cm volunteers); not retained."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Table 2 demographic (median 26.5 OptimAT, 22 volunteers); not retained."
    ),
    SMOKE_CURRENT = list(
      description = "Current-smoker indicator", units = "(binary)", type = "binary",
      notes = paste0(
        "Screened as the 'environmental influences (e.g., smoking status)' ",
        "term in Methods 2.3.2; not retained. Table 2: 24.3 percent smokers ",
        "in OptimAT, 0 percent among volunteers."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened laboratory covariate (Methods 2.3.2); not retained. Table 2 median 27 U/L (5-424)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened laboratory covariate (Methods 2.3.2); not retained. Table 2 median 26 U/L (6-286)."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = "Screened laboratory covariate (Methods 2.3.2); not retained. Table 2 median 8 umol/L (3-41)."
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous",
      notes = paste0(
        "Listed among the screened laboratory measurements in Methods 2.3.2; ",
        "not retained, and no albumin summary appears in Table 2."
      )
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "U/L", type = "continuous",
      notes = paste0(
        "Listed among the screened laboratory measurements in Methods 2.3.2; ",
        "not retained, and no GGT summary appears in Table 2."
      )
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "U/L", type = "continuous",
      notes = "Screened laboratory covariate (Methods 2.3.2); not retained. Table 2 median 67 U/L (14-1406)."
    ),
    BUN = list(
      description = "Blood urea", units = "(not reported)", type = "continuous",
      notes = paste0(
        "Screened as 'urea' (Methods 2.3.2); not retained, and no urea ",
        "summary or unit appears in Table 2."
      )
    ),
    STUDY_OPTIMAT = list(
      description = "1 = OptimAT hospitalized patient, 0 = Geneva cocktail healthy volunteer",
      units = "(binary)", type = "binary",
      notes = paste0(
        "Methods 2.3.2: 'Population type (pooled healthy volunteers vs ",
        "hospitalized patients) was treated as a unique covariate.' It was ",
        "univariately significant on CL/F (dOFV < -12.66, Results 3.2) but ",
        "did NOT survive backward deletion -- 'Subsequent multivariate ",
        "analyses retained age, CLcr, and co-administration of Pgp ",
        "inhibitors' -- so the population difference is fully absorbed by ",
        "the age and CLcr terms."
      )
    ),
    CONMED_PGP_IND = list(
      description = "1 = subject concomitantly treated with a P-glycoprotein inducer, 0 = otherwise",
      units = "(binary)", type = "binary",
      notes = paste0(
        "Collected and categorized (Methods 2.3.2: Pgp-modulating drug use ",
        "was 'categorized as no interaction, Pgp inhibitor, or Pgp inducer'; ",
        "Table 2 reports 2.0 percent inducer use) but never entered into the ",
        "model: the Discussion states 'as only a few patients were treated ",
        "with a Pgp inducer, we have not tested the effect of this covariate ",
        "into the model.' NOTE: CONMED_PGP_IND is NOT a registered canonical ",
        "in inst/references/covariate-columns.md -- it is named here only to ",
        "preserve the paper's covariate screen. A future model that actually ",
        "fits a Pgp-inducer effect should register the canonical then, as the ",
        "sibling of CONMED_PGP_INH."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "fexofenadine", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "fexofenadine", units = "mg",
      specimen = "whole blood", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "fexofenadine", units = "mg",
      specimen = "whole blood", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 459L,
    n_studies      = 2L,
    n_observations = paste0(
      "1089 fexofenadine concentrations (1059 from 449 OptimAT patients, 30 ",
      "from 10 Geneva cocktail healthy volunteers), sampled at baseline and ",
      "2, 3 and 6 h post-dose (Methods 2.2, Results 3.1)"
    ),
    age_range      = "25-97 years (OptimAT patients); 20-36 years (healthy volunteers)",
    age_median     = "71 years (OptimAT patients); 23 years (healthy volunteers)",
    weight_range   = "37-130 kg (OptimAT patients); 69-86 kg (healthy volunteers)",
    weight_median  = "77.2 kg (OptimAT patients); 70.1 kg (healthy volunteers)",
    sex_female_pct = 31,
    race_ethnicity = c(Caucasian = 94.8, Hispanic = 2.0, Asian = 1.3, African = 1.1, Other = 0.9),
    disease_state  = paste0(
      "Hospitalized polymorbid older adults enrolled in the OptimAT ",
      "antithrombotic-optimization cohort (Geneva University Hospitals, ",
      "January 2018 - July 2021), pooled with healthy volunteers from the ",
      "Geneva cocktail study as a normal-Pgp-activity control"
    ),
    renal_function = paste0(
      "Estimated creatinine clearance (Cockcroft-Gault) median 77.2 ",
      "mL/min/1.73 m^2 (range 5.7-147.3). CKD stage distribution in OptimAT ",
      "(Table 2): stage 1 (CLcr > 90) 31.2 percent, stage 2 (60-89) 39.9 ",
      "percent, stage 3 (30-59) 20.1 percent, stage 4 (15-29) 3.7 percent, ",
      "stage 5 (< 15) 5.1 percent"
    ),
    co_medication  = paste0(
      "Pgp inhibitors 14.0 percent, Pgp inducers 2.0 percent (OptimAT, ",
      "Table 2); none in the healthy volunteers. Fexofenadine was dosed as ",
      "part of a six-drug phenotyping cocktail also containing bupropion 25 ",
      "mg, flurbiprofen 25 mg, dextromethorphan 10 mg, omeprazole 5 mg and ",
      "midazolam 1 mg (Methods 2.2)"
    ),
    dose_range     = "single 25 mg oral fexofenadine dose (Geneva cocktail capsule)",
    regions        = "Switzerland (single centre, Geneva University Hospitals)",
    notes          = paste0(
      "Concentrations were measured in 10 uL capillary whole-blood dried ",
      "blood spots collected by finger prick onto Whatman 903 cards and ",
      "quantified by LC-MS/MS (Methods 2.2; assay of Bosilkovska 2014). The ",
      "paper's text, tables and title nonetheless refer throughout to ",
      "'plasma concentrations' while the Conclusion says 'blood exposure'; ",
      "the measured matrix is whole blood, which is what compartmentData ",
      "records. The model was fit by SAEM in MONOLIX 2021.1.0 and evaluated ",
      "by prediction-corrected VPC (n = 500 simulations) and nonparametric ",
      "bootstrap (n = 1000, Rsmlx); Results 3.3 reports every parameter ",
      "within its bootstrap 95 percent CI with weighted difference < 20 ",
      "percent. The 6 h sampling window is a validated Pgp-phenotyping ",
      "design and does not span the full absorption and elimination phases ",
      "(Discussion, limitations), so the absorption parameters are not ",
      "identifiable from these data and were carried from Piscitelli 2023."
    )
  )

  ini({
    # ========================================================================
    # Absorption -- both parameters carried unchanged from Piscitelli 2023
    # (Gaspar 2025 reference 12) and not re-estimated here. Gaspar 2025
    # Methods 2.3.1: "Literature-sourced values for D1 and ka2 were 1.59
    # hours and 0.282 h-1, respectively [12]." Table 3 lists both with an
    # asterisk, whose footnote reads "Parameters marked with an asterisk (*)
    # were fixed during the modeling process".
    # ========================================================================
    ld1 <- fixed(log(1.59))
    label("Zero-order absorption duration D1 (h); the sequential zero-order input precedes the first-order phase; literature value from Piscitelli 2023")
    # Gaspar 2025 Methods 2.3.1 and Table 3 row "D1". Table 3's row label
    # reads "D1 (h -1)", which is a typographical slip -- a duration cannot
    # carry reciprocal-time units, and Methods 2.3.1 states "1.59 hours".

    lka <- fixed(log(0.282))
    label("First-order absorption rate constant ka2 (1/h); literature value from Piscitelli 2023")
    # Gaspar 2025 Methods 2.3.1 gives 0.282 1/h; Table 3 row "k a (h -1)"
    # rounds the same value to 0.28. The 3-significant-figure Methods value
    # is used here because it is the one attributed to the source
    # publication. Mean absorption time = D1/2 + 1/ka = 0.795 + 3.55 = 4.34 h.

    # ========================================================================
    # Disposition -- Gaspar 2025 Table 3, "Final model Estimate (RSE %)"
    # column. Bootstrap medians and 95% CIs from the adjacent column are
    # quoted for reference; they agree closely (Results 3.3).
    # ========================================================================
    lcl <- log(116)
    label("Apparent clearance CL/F (L/h) for the reference subject (CLcr 77.2 mL/min/1.73 m^2, age 71 years, no Pgp inhibitor)")
    # Gaspar 2025 Table 3: CL = 116 L/h (RSE 4.1%); bootstrap median 116
    # (95% CI 110-122). The covariate-free base model gave CL/F = 114.64 L/h
    # (RSE 3.99%, Results 3.2), which the centred covariate model reproduces
    # to 1.2% -- see the `e_crcl_cl` comment on the centring choice.

    lvc <- log(12)
    label("Apparent central volume of distribution V1/F (L)")
    # Gaspar 2025 Table 3: V1 = 12 L (RSE 19.9%); bootstrap median 12
    # (95% CI 9.0-15.0). Base model V1/F was 14.9 L (RSE 19.5%).

    lq <- log(44)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Gaspar 2025 Table 3: Q = 44 L/h (RSE 7.2%); bootstrap median 44.13
    # (95% CI 40.5-47.8). Base model Q/F was 52.02 L/h (RSE 8.71%).

    lvp <- log(157)
    label("Apparent peripheral volume of distribution V2/F (L)")
    # Gaspar 2025 Table 3: V2 = 157 L (RSE 8.9%); bootstrap median 156
    # (95% CI 140-174). Base model V2/F was 145.7 L (RSE 8.61%).

    # ========================================================================
    # Covariate effects on log(CL/F). Gaspar 2025 Table 3 footnote states the
    # final model equation verbatim:
    #
    #   log(CLi) = log(CLpop) + beta_GFR * log(CLcr_i / CLcr_m)
    #                         + beta_AGE * log(AGE_i / AGE_m)
    #                         + beta_PGP * PGP_i + eta_CL
    #
    # "where ... CLcr_m and AGE_m the median values of CLcr and AGE in the
    # population, and eta_CL the interindividual variability of CL. PGP_i = 1
    # if individual i is taking a Pgp inhibitor and 0 otherwise."
    # ========================================================================
    e_crcl_cl <- 0.33
    label("Power exponent on (CRCL / 77.2) for CL/F (unitless)")
    # Gaspar 2025 Table 3: beta CLcr = 0.33 (RSE 25.9%); bootstrap median 0.32
    # (95% CI 0.22-0.42). The Table 3 footnote quotes 0.32 (the bootstrap
    # median) rather than the 0.33 point estimate; the Estimate column is
    # used here.
    #
    # CENTRING VALUE (77.2): the Table 3 footnote defines the divisor as
    # "CLcr_m ... the median value of CLcr in the population", which Table 2
    # gives as 77.2 mL/min/1.73 m^2. Methods Eq. 1 instead defines the
    # divisor as a "typical value" and names "CLcr (100 mL/min)". The
    # footnote equation is the paper's own statement of the FINAL model and
    # is therefore preferred, and it is corroborated by intercept continuity:
    # centring at the median makes CLpop = 116 L/h essentially the base
    # model's covariate-free CL/F of 114.64 L/h (1.2% apart), whereas
    # centring at 100 would place the typical-subject CL/F at 106.5 L/h, 7%
    # below the base model. The choice is also low-impact: reproducing the
    # Table 4 AUC(0-6) ratio grid gives a mean absolute error of 6.0% with
    # 77.2 and 4.8% with 100 (see the validation vignette), because the
    # divisor largely cancels in a between-subject ratio.

    e_age_cl <- -0.59
    label("Power exponent on (AGE / 71) for CL/F (unitless)")
    # Gaspar 2025 Table 3 footnote: "The coefficients beta CLcr, beta AGE and
    # beta PGP quantify the impact of CLcr, AGE and Pgp on clearance, with
    # values of 0.32, - 0.59, and - 0.39, respectively."
    #
    # VALUE CONFLICT WITHIN TABLE 3, resolved against the paper's own
    # simulations. The Table 3 body reports beta Age = -0.17 (RSE 23.5%) with
    # a bootstrap median of -0.17 (95% CI -0.23 to -0.11), i.e. the Estimate,
    # RSE and both bootstrap columns are mutually consistent at -0.17. The
    # Table 3 footnote equation, which is the paper's authoritative statement
    # of the fitted final model, gives -0.59. Three independent lines of
    # evidence select -0.59:
    #   (1) Results 3.2: "For a 90-year-old patient with normal renal
    #       function, CL/F was reduced by approximately 59% compared with a
    #       20-year-old patient with similar renal function." Under the
    #       footnote equation (90/20)^-0.59 = 0.412, a 58.8% reduction --
    #       an exact match. Under -0.17 it is (90/20)^-0.17 = 0.774, a 22.6%
    #       reduction, which contradicts the prose.
    #   (2) The Table 4 AUC(0-6) ratio grid, generated by simulating the
    #       fitted model, is reproduced to a mean absolute error of 4.8-6.0%
    #       with -0.59 but is missed by 25-26% on average and by up to 40% in
    #       the oldest age decades with -0.17, which flattens the age trend
    #       that is the paper's central finding.
    #   (3) The Section 3.4 absolute AUC(0-6) values show the same pattern.
    # The vignette carries the falsification chunk that produces (2) and (3)
    # so the choice is auditable and reversible. -0.17 is treated as a
    # transcription error in the Table 3 body.

    e_conmed_pgp_inh_cl <- -0.38
    label("Log-additive effect of concomitant Pgp-inhibitor use on CL/F (unitless)")
    # Gaspar 2025 Table 3: beta Pgp = -0.38 (RSE 24.2%); bootstrap median
    # -0.39 (95% CI -0.50 to -0.28). The Table 3 footnote quotes -0.39 (the
    # bootstrap median); the Estimate column is used here. exp(-0.38) = 0.684,
    # a 31.6% reduction in CL/F, which reproduces the Abstract's
    # "Co-administration of Pgp inhibitors led to a 35% increase in AUC 0-6"
    # (the Table 4 stage-1 / age-20-30 pair gives 1.55/1.15 = 1.35).
    # Results 3.2 additionally claims this effect "reduced CL/F by 48%",
    # which would require beta Pgp = -0.65; that prose figure is inconsistent
    # with both Table 3 columns, with the Table 3 footnote, and with the
    # Abstract's 35% AUC increase, and is not used.

    # ========================================================================
    # Interindividual variability -- retained on CL/F only. Results 3.2: "IIV
    # was retained only for CL/F, and a proportional error model was used to
    # characterize the residual variability."
    # ========================================================================
    etalcl ~ 0.514
    # Gaspar 2025 Table 3: IIV CL = 82 CV% (RSE 3.9%); bootstrap median 83
    # (95% CI 78.5-87.5). Converted to the log-scale variance nlmixr2 stores
    # with omega^2 = log(CV^2 + 1) = log(0.82^2 + 1) = log(1.6724) = 0.5143.

    # ========================================================================
    # Residual variability
    # ========================================================================
    propSd <- 0.23
    label("Proportional residual error (fraction)")
    # Gaspar 2025 Table 3: proportional error model 23 CV% (RSE 2.9%);
    # bootstrap median 23 with a 95% CI printed as "0.21 to 0.25", i.e. the
    # CI is on the fraction scale while the point estimate is in CV% -- the
    # two agree once the scales are reconciled (0.23 = 23%).
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual apparent clearance. Log-additive covariate model exactly
    #    as printed in the Gaspar 2025 Table 3 footnote equation. The
    #    divisors are the population medians from Table 2: CLcr_m = 77.2
    #    mL/min/1.73 m^2 and AGE_m = 71 years.
    # ----------------------------------------------------------------------
    cl <- exp(lcl + etalcl
              + e_crcl_cl * log(CRCL / 77.2)
              + e_age_cl * log(AGE / 71)
              + e_conmed_pgp_inh_cl * CONMED_PGP_INH)

    # ----------------------------------------------------------------------
    # 2. Remaining individual parameters. No IIV was retained on any of
    #    these (Results 3.2), so each is the typical value.
    # ----------------------------------------------------------------------
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)
    d1 <- exp(ld1)

    # ----------------------------------------------------------------------
    # 3. Micro-constants.
    # ----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 4. Two-compartment ODE system with linear elimination from central
    #    (Gaspar 2025 Fig. 1).
    # ----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------------
    # 5. Sequential zero-order plus first-order absorption (Gaspar 2025
    #    Methods 2.3.1 and Fig. 1): the dose is delivered into the depot
    #    uniformly over [0, D1] (zero-order phase) and the depot then drains
    #    into central at ka2 (first-order phase). Dose records must carry
    #    rate = -2 for rxode2 to apply the modelled duration. This is the
    #    only reading of "sequential zero- and first-order absorption" that
    #    is identifiable from the two parameters the paper reports: a split
    #    zero-order / first-order pathway would additionally require an
    #    unreported fraction, and a zero-order input of the whole dose
    #    directly into central would leave ka2 with nothing to act on.
    #    Same encoding as Heathman_2024_efavirenz.
    # ----------------------------------------------------------------------
    dur(depot) <- d1

    # ----------------------------------------------------------------------
    # 6. Observation. Amounts are in mg and volumes in L, so central / vc is
    #    mg/L; the factor 1000 converts to the ng/mL scale in which Gaspar
    #    2025 reports concentrations and AUC(0-6) (Table 1, Section 3.4).
    # ----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
