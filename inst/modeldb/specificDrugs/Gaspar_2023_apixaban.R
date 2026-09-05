Gaspar_2023_apixaban <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model with first-order ",
    "absorption and an absorption lag time for oral apixaban in a real-life ",
    "cohort of hospitalized, mostly elderly patients from the OptimAT study ",
    "(Gaspar 2023). Apparent oral clearance CL/F = 3.2 L/h at a ",
    "Cockcroft-Gault creatinine clearance of 100 mL/min and a P-glycoprotein ",
    "phenotyping index (fexofenadine AUC from 2 to 6 h) of 205 mg*h/L, ",
    "scaled by the power terms (CRCL/100)^0.47 and (AUC_FEXO/205)^-0.17. ",
    "Renal function and P-gp phenotypic activity together explained 58% of ",
    "the between-patient variance in clearance (41% and 17%, ",
    "respectively). Apparent central and peripheral volumes are 27 L and 51 ",
    "L, apparent intercompartmental clearance 16 L/h, absorption rate ",
    "constant 0.82 1/h and lag time 0.16 h. Variability is carried as ",
    "interoccasion variability on clearance (30% CV, one occasion per ",
    "hospitalization) rather than interindividual variability, with ",
    "interindividual variability on lag time (108% CV), absorption rate ",
    "(58% CV), intercompartmental clearance (65% CV) and peripheral volume ",
    "(88% CV); interindividual variability on the central volume was tested ",
    "and not retained. Residual variability is proportional (8.4% CV). ",
    "Neither CYP3A4/5 phenotypic activity nor any demographic or ",
    "laboratory covariate was retained. Parameter values are taken from ",
    "the publication's Table 2 (Final model column) and the final-model ",
    "equation printed in the Table 2 note."
  )
  reference <- paste0(
    "Gaspar F, Terrier J, Favre S, Gosselin P, Fontana P, Daali Y, ",
    "Lenoir C, Samer CF, Rollason V, Reny JL, Csajka C, Guidi M. ",
    "Population pharmacokinetics of apixaban in a real-life hospitalized ",
    "population from the OptimAT study. ",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1541-1552. ",
    "doi:10.1002/psp4.13032."
  )
  vignette <- "Gaspar_2023_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Apixaban was measured in whole blood collected as
  # dried blood spots and converted to a PLASMA concentration by the
  # Foerster relationship Conc_plasma = Conc_DBS * 1.46 + 13.84 (Gaspar 2023
  # Methods "Analytical method"), so the model was fitted to, and predicts,
  # plasma concentrations.
  compartmentData <- list(
    depot       = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "RAW Cockcroft-Gault creatinine clearance in mL/min, NOT ",
        "BSA-normalized (Gaspar 2023 Methods 'Covariate model' lists 'creatinine ",
        "clearance [CLcr; Cockcroft and Gault equation]' among the laboratory ",
        "factors, and Table 1 reports it in mL/min). Normalized to the typical ",
        "value of 100 mL/min stated in Methods 'Covariate model' ('Cov_weight is ",
        "Cov weighted mean or typical value in case of body weight (70 kg) and ",
        "CLcr (100 mL/min)') and printed in the Table 2 note as the denominator ",
        "of the final-model equation. Cohort median 57.2 mL/min, range 23-136 ",
        "(Table 1), i.e. a renal-impairment-rich hospitalized population; the ",
        "paper's own simulations nevertheless extrapolate the term down to ",
        "chronic-kidney-disease stage 5 (CLcr < 15 mL/min). Enters clearance ",
        "only, as the power term (CRCL / 100)^0.47. Note that body weight and ",
        "age were each significant on CL in univariable screening but were ",
        "dropped in the multivariable step; the Discussion attributes this to ",
        "both being inputs to the Cockcroft-Gault equation, so body size and age ",
        "reach clearance THROUGH this column and must not be added again as a ",
        "separate allometric term."
      ),
      source_name        = "CLcr"
    ),
    AUC_FEXO = list(
      description        = "Fexofenadine area under the curve from 2 to 6 h, the Geneva cocktail P-glycoprotein phenotyping index",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Individual P-glycoprotein phenotypic activity measured with the Geneva ",
        "cocktail (Bosilkovska 2014) as the area under the fexofenadine ",
        "concentration-time curve from 2 to 6 h after the probe dose (Gaspar 2023 ",
        "Methods 'Phenotyping'). ORIENTATION IS INVERSE: fexofenadine is a P-gp ",
        "substrate, so a HIGHER AUC_FEXO means LOWER P-gp activity. The Methods ",
        "give the phenotype-category means (mean +/- SD): 285.5 +/- 67.1 mg*h/L ",
        "for reduced-activity poor metabolizers, 100.1 +/- 47.5 for normal ",
        "metabolizers and 50.4 +/- 15.3 for increased-activity ultra-rapid ",
        "metabolizers; the OptimAT cohort was 47% PM, 51% NM and 2% UM (Table 1). ",
        "Normalized to 205 mg*h/L, the constant printed in the Table 2 ",
        "final-model equation. Enters clearance only, as the power term ",
        "(AUC_FEXO / 205)^-0.17, so a poor metabolizer has 16% lower and an ",
        "ultra-rapid metabolizer 12% higher clearance than a normal metabolizer ",
        "(Gaspar 2023 Results 'Covariate analyses'). The value is probe- and ",
        "window-specific: a P-gp activity measured with a different probe ",
        "(digoxin, talinolol) or over a different window is NOT numerically ",
        "interchangeable with this column."
      ),
      source_name        = "P-gp AUC2-6"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for inter-occasion-variability multiplexing",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste0(
        "One occasion is one hospitalization: 'One occasion refers to one ",
        "hospitalization, with six patients being hospitalized twice' (Gaspar ",
        "2023 Methods 'Population pharmacokinetic analysis'), so the 100 ",
        "patients contributed 106 occasions and the maximum occasion count in ",
        "the source data is TWO. Two occasions are therefore encoded here, ",
        "against the four of Blackman_2026_methotrexate.R. The paper fits a ",
        "single inter-occasion variance on clearance shared across occasions ",
        "(Table 2 row 'IOV_CL (CV%)'), so the second variance is fix()ed equal ",
        "to the first and extending the chain to more occasions would be a ",
        "mechanical copy of that line. Decomposed inside model() into binary ",
        "indicators oc1 and oc2 that multiplex the two IOV etas on log-CL. For ",
        "single-occasion records pass OCC = 1 so the first IOV eta applies."
      ),
      source_name        = "occasion"
    )
  )

  # Covariates that Gaspar 2023 screened but did NOT retain in the final
  # model. Documentation only: they are deliberately absent from model(),
  # because the paper publishes no point estimate for any of them.
  covariatesDataExcluded <- list(
    CYP3A4 = list(
      description        = "CYP3A4/5 phenotypic activity, Geneva cocktail midazolam metabolic ratio at 2 h",
      units              = "(ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Measured as the 1'-hydroxymidazolam / midazolam concentration ratio at ",
        "2 h (Gaspar 2023 Methods 'Phenotyping'); category means (mean +/- SD) ",
        "75.3 +/- 17.5 for poor metabolizers, 27.0 +/- 14.8 for normal and 13.8 ",
        "+/- 4.0 for ultra-rapid metabolizers, with the cohort 19% PM, 78% NM, ",
        "3% UM (Table 1). NOT retained: 'Neither CYP3A4/5 phenotype (dOFV = 2.8, ",
        "p = 0.08) ... showed any influence on apixaban PK' (Results 'Covariate ",
        "analyses'). The Discussion attributes the null result to the small ",
        "number of PM and UM carriers rather than to a genuine absence of a ",
        "CYP3A4/5 effect, so this is a power limitation and not a negative ",
        "finding. Sibling to the retained AUC_FEXO P-gp phenotyping index, ",
        "measured in the same cocktail administration."
      ),
      source_name        = "CYP3A4/5 ratio"
    ),
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL (dOFV = -4.6, p < 0.05) and on Q (dOFV = -6.0, p < 0.05) in univariable screening but not retained after the multivariable and backward-deletion steps; no point estimate published. Median 75 kg, range 44-126 (Table 1). Body weight reaches clearance indirectly through CRCL, which is Cockcroft-Gault-derived.",
      source_name = "Body weight"
    ),
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL (dOFV = -20, p < 0.01), Vp (dOFV < -4.0, p < 0.05) and Q (dOFV < -5.5, p < 0.02) in univariable screening but not retained; no point estimate published. Median 77 years, range 51-94 (Table 1). Age reaches clearance indirectly through CRCL.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)", units = "(binary)", type = "binary",
      reference_category = "0 (male)",
      notes = "Reported as 'gender'; significant on Vp (dOFV < -4.0, p < 0.05) and Q (dOFV < -5.5, p < 0.02) in univariable screening but not retained; no point estimate published. 42% female (Table 1).",
      source_name = "Gender"
    ),
    SMOKE = list(
      description = "Current-smoker indicator", units = "(binary)", type = "binary",
      reference_category = "0 (non-smoker)",
      notes = "Significant on CL in univariable screening (dOFV = -10, p < 0.01) but not retained; no point estimate published and the cohort smoking prevalence is not tabulated.",
      source_name = "Smoking status"
    ),
    BUN = list(
      description = "Blood urea", units = "mmol/L", type = "continuous",
      reference_category = NULL,
      notes = "Reported by the paper as 'urea'. Significant on CL in univariable screening (dOFV = -12, p < 0.01) but not retained; no point estimate published. Median 7.3 mmol/L, range 1.1-37.7 (Table 1). Collinear with CRCL.",
      source_name = "Urea"
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL in univariable screening (dOFV = -8.8, p < 0.01) but not retained; no point estimate published. Median 37 g/L, range 24-62 (Table 1).",
      source_name = "Albumin"
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL in univariable screening (dOFV = -5.8, p < 0.05) but not retained; no point estimate published. Median 21 U/L, range 9-286 (Table 1).",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Screened and not significant ('any of the other tested covariates (dOFV > -3.7, p > 0.05)'); no point estimate published. Median 22 U/L, range 5-424 (Table 1).",
      source_name = "ALT"
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Screened and not significant; no point estimate published. Median 69 U/L, range 14-295 (Table 1).",
      source_name = "Alkaline phosphatase"
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "U/L", type = "continuous",
      reference_category = NULL,
      notes = "Screened and not significant; no point estimate published. Median 41 U/L, range 6-662 (Table 1).",
      source_name = "Gamma-glutamyl transferase"
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      reference_category = NULL,
      notes = "Screened and not significant; no point estimate published. Median 8 umol/L, range 2.9-37 (Table 1).",
      source_name = "Bilirubin"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 1L,
    n_observations = 825L,
    n_occasions    = 106L,
    age_range      = "51-94 years",
    age_median     = "77 years",
    weight_range   = "44-126 kg",
    weight_median  = "75 kg",
    sex_female_pct = 42,
    disease_state  = "hospitalized patients receiving apixaban; 89% for atrial fibrillation, 11% for venous thromboembolism",
    renal_function = "Cockcroft-Gault creatinine clearance median 57.2 mL/min, range 23-136 mL/min",
    dose_range     = "2.5 mg b.i.d. (n = 40), 5 mg b.i.d. (n = 56), 10 mg b.i.d. (n = 4)",
    regions        = "Switzerland (Geneva University Hospitals)",
    notes          = paste0(
      "OptimAT study (NCT03477331), a single-centre prospective observational ",
      "study of patients admitted to the Geneva University Hospitals between ",
      "January 2018 and November 2019 with a prescription of apixaban. Baseline ",
      "demographics and laboratory values are in Gaspar 2023 Table 1. Each ",
      "participant underwent one day of multiple finger-prick dried-blood-spot ",
      "sampling after the morning dose, scheduled predose and at 0.5, 1, 2, 3, ",
      "4, 6 and 8 h; a median of 8 samples per individual (range 7-8) was ",
      "available. Phenotype distribution: CYP3A4/5 19% PM / 78% NM / 3% UM and ",
      "P-gp 47% PM / 51% NM / 2% UM. Model building used OptimAT only; an ",
      "independent cohort of 63 patients from the DAPHNE study (NCT03112525) ",
      "was used for external validation and gave a mean prediction error of 0% ",
      "(95% CI -1% to 1%) with a precision of 11%."
    )
  )

  ini({
    # Structural parameters, Gaspar 2023 Table 2 "Final model Estimate (RSE %)".
    # CL/F, Vc/F, Q/F and Vp/F are APPARENT (oral-only data; bioavailability is
    # not identifiable and is absorbed into each of them), which is how the
    # paper itself labels them in Results "Base model".
    ltlag <- log(0.16); label("Absorption lag time (h)")                        # Table 2, Tlag = 0.16 h (RSE 19%; bootstrap 0.18 [0.07-0.28])
    lka   <- log(0.82); label("First-order absorption rate constant (1/h)")     # Table 2, ka = 0.82 1/h (RSE 11%; bootstrap 0.92 [0.54-1.87])
    lcl   <- log(3.2);  label("Apparent oral clearance CL/F (L/h)")             # Table 2, CL = 3.2 L/h (RSE 6%; bootstrap 3.1 [2.4-3.6]); typical value at CRCL = 100 mL/min and AUC_FEXO = 205 mg*h/L per the Table 2 note equation
    lvc   <- log(27);   label("Apparent central volume of distribution Vc/F (L)")     # Table 2, Vc = 27 L (RSE 8%; bootstrap 29 [18-44])
    lq    <- log(16);   label("Apparent intercompartmental clearance Q/F (L/h)")      # Table 2, Q = 16 L/h (RSE 12%; bootstrap 14 [8.0-27])
    lvp   <- log(51);   label("Apparent peripheral volume of distribution Vp/F (L)")  # Table 2, Vp = 51 L (RSE 12%; bootstrap 56 [25-118])

    # Covariate effects on clearance. Both enter as centred power terms in the
    # final-model equation printed in the Table 2 note:
    #   CL_occ,i = CL * (CLcr/100)^beta_CLcr * (AUC2-6/205)^beta_PgP * exp(IOV_CL,occ,i)
    e_crcl_cl     <- 0.47;  label("Power exponent of creatinine clearance on CL/F (unitless)")            # Table 2, beta_CLcr = 0.47 (RSE 19%; bootstrap 0.47 [0.22-0.74]); centred at 100 mL/min
    e_auc_fexo_cl <- -0.17; label("Power exponent of the fexofenadine AUC2-6 P-gp index on CL/F (unitless)")  # Table 2, beta_P-gp = -0.17 (RSE 32%; bootstrap -0.21 [-0.44 to -0.05]); centred at 205 mg*h/L

    # Interindividual variability. Table 2 reports each as a CV%; for the
    # log-normally distributed parameters assumed throughout (Methods "Base
    # model": "A log-normal distribution was assumed for all PK parameters"),
    # the variance is omega^2 = log(CV^2 + 1).
    # NOTE: there is deliberately NO eta on CL and NO eta on Vc. Results "Base
    # model" states that IOV on CL replaced IIV on CL ("The use of IOV instead
    # of IIV on apparent clearance (CL) ... markedly improved data description")
    # and that IIV on Vc was discarded once the peripheral compartment was added
    # ("allowed discarding IIV on Vc (dOFV = 1.7, p > 0.05)").
    etaltlag ~ 0.77307   # Table 2, IIV Tlag = 108% CV; log(1.08^2 + 1) = 0.77307
    etalka   ~ 0.28998   # Table 2, IIV ka   =  58% CV; log(0.58^2 + 1) = 0.28998
    etalq    ~ 0.35242   # Table 2, IIV Q    =  65% CV; log(0.65^2 + 1) = 0.35242
    etalvp   ~ 0.57346   # Table 2, IIV Vp   =  88% CV; log(0.88^2 + 1) = 0.57346

    # Interoccasion variability on clearance, one occasion per hospitalization.
    etaiov_cl_1 ~ 0.08618       # Table 2, IOV_CL = 30% CV; log(0.30^2 + 1) = 0.08618 (estimated, occasion 1)
    etaiov_cl_2 ~ fix(0.08618)  # occasion-2 variance equal to occasion-1: the paper reports ONE shared IOV variance

    # Residual unexplained variability. Results "Base model": "The proportional
    # error model adequately captured the RUV."
    propSd <- 0.084; label("Proportional residual error (fraction)")  # Table 2, proportional error model = 8.4% CV (RSE 3%; bootstrap 8.4 [7.5-9.4])
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators to
    #    multiplex the two inter-occasion-variability etas on log-CL. For
    #    single-occasion data pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # 2. Individual parameters. Clearance carries the two retained covariate
    #    power terms and the occasion-specific IOV eta but NO interindividual
    #    eta; the central volume carries neither a covariate nor an eta.
    #    Reference covariates: CRCL 100 mL/min, AUC_FEXO 205 mg*h/L.
    cl   <- exp(lcl + iov_cl) * (CRCL / 100)^e_crcl_cl * (AUC_FEXO / 205)^e_auc_fexo_cl
    vc   <- exp(lvc)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment ODE system with first-order oral absorption.
    #    Bioavailability F is not identifiable from oral-only data and is
    #    absorbed into the apparent parameters CL/F, Vc/F, Q/F and Vp/F, so no
    #    f(depot) term is applied.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Absorption lag time. Results "Base model": "The inclusion of lag-time
    #    and associated IIV in addition to ka (Tlag, dOFV = -104, p < 0.01)
    #    captured the apixaban absorption process more adequately than the
    #    sequential and simultaneous first- and zero-order absorption models."
    alag(depot) <- tlag

    # 6. Observation and error.
    #    central is in mg and vc in L, so central / vc is mg/L; multiply by 1000
    #    to report ng/mL (1 mg/L = 1 ug/mL = 1000 ng/mL), matching the paper's
    #    reporting units for concentrations, the 50 ng/mL pre-procedural target
    #    and the 230 ng/mL bleeding-risk threshold.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
