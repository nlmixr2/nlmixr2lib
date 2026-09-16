Tao_2026_levofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for intravenous levofloxacin",
    "in Chinese pediatric patients (0.16-16 years) with severe refractory Mycoplasma pneumoniae",
    "pneumonia (Tao 2026). Body weight enters both clearance and central volume as an estimated",
    "power term referenced to the cohort mean of 24.05 kg, with exponents of 1.11 on clearance and",
    "1.83 on central volume, i.e. steeper than the conventional allometric 0.75 / 1. Serum",
    "creatinine enters clearance as a second power term referenced to the cohort mean of 31.77",
    "umol/L with an exponent of -0.20; estimated glomerular filtration rate dropped the objective",
    "function almost identically (-215.89 versus -216.15) but was not retained. Age was screened",
    "and rejected, and a maturation function on clearance did not improve the fit. Interindividual",
    "variability is carried on clearance only; the authors fixed it to zero on central volume,",
    "peripheral volume and intercompartmental clearance because of high shrinkage, so this model",
    "carries a single eta. Residual error is proportional. The companion exposure-response analysis",
    "(a Cox proportional hazards model, not reproducible as an ODE and therefore not encoded here)",
    "identified a steady-state AUC(0-24 h) of 30.74 mg*h/L as the cutoff above which levofloxacin",
    "treatment duration was shorter, and 33.72 mg*h/L for faster cough resolution; the authors used",
    "the former as the target for Monte Carlo dose finding and recommend 11 mg/kg q24h or 5.5 mg/kg",
    "q12h over the guideline 8-10 mg/kg regimens.",
    sep = " "
  )
  reference <- paste(
    "Tao X, Zhou Y, Xu S, Su Y, Tang Y, Zhang S, Chang Z, Xie F, Lv M.",
    "Population pharmacokinetics and exposure-response analysis of levofloxacin in Chinese",
    "pediatric patients with severe refractory Mycoplasma pneumoniae pneumonia.",
    "Antimicrob Agents Chemother. 2026;70(6):e01853-25. doi:10.1128/aac.01853-25",
    sep = " "
  )
  vignette <- "Tao_2026_levofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "levofloxacin", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    # Empirical distribution compartment. Only plasma was assayed (Methods
    # 'Levofloxacin dosing and concentration analysis': scavenged plasma,
    # HPLC-MS/MS), so the peripheral state is a lumped distribution volume
    # rather than a named matrix; it is labelled plasma to match the sibling
    # levofloxacin two-compartment entries.
    peripheral1 = list(
      analyte = "levofloxacin", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mean 24.05 kg (SD 9.00), median 22.50 kg (IQR 18.00-28.00), Tao 2026 Table 1.",
        "The covariate equations printed beneath Table 2 normalize to the MEAN, 24.05 kg,",
        "not the median. Body weight was the single body-size descriptor carried into the",
        "stepwise screen because it correlates strongly with age (r = 0.84) and height",
        "(r = 0.88) (Results 'PopPK analysis of levofloxacin', Fig. S2), so age and height",
        "were deliberately excluded to avoid collinearity. Retained on both clearance",
        "(exponent 1.11) and central volume (exponent 1.83). The Monte Carlo scenarios in",
        "Table 3 span 9-45 kg."
      ),
      source_name        = "BW"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mean 31.77 umol/L (SD 8.03), median 30.70 umol/L (IQR 26.80-37.00), Tao 2026 Table 1.",
        "The Table 2 covariate equation normalizes to the MEAN, 31.77 umol/L. Reported and",
        "modelled in umol/L, NOT mg/dL: divide by 88.4 to convert (31.77 umol/L = 0.36 mg/dL,",
        "the expected order of magnitude for a healthy 7-year-old). Retained on clearance only,",
        "as a power term with the estimated exponent -0.20. Levofloxacin is almost entirely",
        "eliminated unchanged by the kidneys (Discussion), so this is the model's renal-function",
        "term. The Monte Carlo scenarios in Table 3 span 19-45 umol/L; the cohort contained no",
        "renally impaired children, so do not extrapolate the power term to elevated creatinine.",
        "eGFR by the modified Schwartz equation was the competing renal descriptor and dropped",
        "the objective function almost as far (-215.89 versus -216.15) but was not retained."
      ),
      source_name        = "SCr"
    )
  )

  # Covariates screened and NOT retained in the final model (Results 'PopPK
  # analysis of levofloxacin': 'Other tested covariates, including AST, ALT,
  # TBIL, DBIL, ALB, and eGFR, were not significant'). Age and height were
  # excluded before screening rather than screened and rejected, because of
  # their collinearity with body weight (Fig. S2), and are recorded here for
  # the same reason: a user reading this model should be able to see that the
  # absence of an age term is a modelling decision, not an oversight.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years", type = "continuous",
      notes = paste(
        "Mean 7.07 years (SD 2.77), median 6.92 (IQR 5.25-8.58), full range 0.16-16 years,",
        "Tao 2026 Table 1. Not entered into the stepwise screen because of its correlation with",
        "body weight (r = 0.84, Fig. S2). The Discussion records that age was separately found",
        "not to predict clearance and that 'the inclusion of maturation of CL with age did not",
        "further improve the data fit', which the authors attribute to only 7 patients being",
        "under 2 years old (11 concentrations). Several other pediatric levofloxacin models do",
        "retain an age effect (Table S1: Denti 2018, Garcia-Prats 2019, van der Laan 2021,",
        "White 2024), so the absence of one here is a property of this cohort."
      )
    ),
    HT = list(
      description = "Height.",
      units = "cm", type = "continuous",
      notes = paste(
        "Collected from medical records (Methods) and plotted in the Fig. S2 covariate",
        "correlation matrix, but no summary statistic is printed in Table 1. Not entered into",
        "the stepwise screen because of its correlation with body weight (r = 0.88, Fig. S2)."
      )
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the modified Schwartz equation.",
      units = "mL/min/1.73 m^2", type = "continuous",
      notes = paste(
        "Mean 149.61 (SD 34.48), median 144.30 (IQR 125.58-167.43) mL/min/1.73 m^2,",
        "Tao 2026 Table 1. Screened against clearance and, on its own, significant: it dropped",
        "the objective function by 215.89 points against the base model. Serum creatinine",
        "dropped it by 216.15 and was chosen instead 'taking into account both statistical",
        "performance and clinical practicality' (Discussion). The two are alternative encodings",
        "of the same renal signal, so only one is carried; use CREAT with this model."
      )
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L", type = "continuous",
      notes = "Mean 38.72 g/L (SD 4.98), Tao 2026 Table 1. Screened and not significant (Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L", type = "continuous",
      notes = "Mean 34.87 U/L (SD 26.14), Tao 2026 Table 1. Screened and not significant (Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L", type = "continuous",
      notes = "Mean 36.39 U/L (SD 54.24), Tao 2026 Table 1. Screened and not significant (Results)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "umol/L", type = "continuous",
      notes = "Mean 5.71 umol/L (SD 2.22), Tao 2026 Table 1. Screened and not significant (Results)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin.",
      units = "umol/L", type = "continuous",
      notes = "Mean 1.86 umol/L (SD 0.68), Tao 2026 Table 1. Screened and not significant (Results)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)", type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "91 of 191 patients female (47.64%), Tao 2026 Table 1. Sex is the one Table 1 variable",
        "explicitly excluded from the Pearson correlation screen ('Correlations between every two",
        "potential covariates listed in Table 1, with the exception of sex') and does not appear",
        "in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 191,
    n_studies      = 1,
    age_range      = "0.16-16 years",
    age_median     = "6.92 years (IQR 5.25-8.58); mean 7.07 (SD 2.77)",
    weight_median  = "22.50 kg (IQR 18.00-28.00); mean 24.05 (SD 9.00)",
    weight_range   = "IQR 18.00-28.00 kg; the full range is not reported, but the Monte Carlo scenarios span 9-45 kg",
    sex_female_pct = 47.64,
    race_ethnicity = c(Asian = 100),
    disease_state  = "severe refractory Mycoplasma pneumoniae pneumonia (SRMPP), defined per the 2023 Chinese guidelines",
    renal_function = paste(
      "serum creatinine mean 31.77 umol/L (SD 8.03); estimated glomerular filtration rate by the",
      "modified Schwartz equation mean 149.61 mL/min/1.73 m^2 (SD 34.48). No renally impaired",
      "children were enrolled and the authors note the eGFR range was narrow."
    ),
    co_medication  = "all patients received concomitant glucocorticoids as part of standard care; the authors caution against extrapolating the exposure-response relationship to children not receiving them",
    dose_range     = "intravenous levofloxacin 8-10 mg/kg per dose, q12h if under 5 years and q24h if 5 years or older, not exceeding 750 mg/day, infused over 0.5-2 h; mean daily dose 11.32 mg/kg (SD 3.06)",
    regions        = "China (single center, Children's Hospital Affiliated to Zhengzhou University)",
    notes          = paste(
      "Single-center, prospective, open-label PK/PD study run between April 2023 and April 2024",
      "(ethics approval 2022-K-L061). Baseline demographics are Tao 2026 Table 1. The PK data set",
      "is sparse and opportunistic: 293 scavenged plasma samples from 191 patients, roughly 1.5",
      "samples per child, which is why interindividual variability on the distribution parameters",
      "could not be estimated. Only 7 patients were under 2 years of age (11 concentrations), so",
      "the authors caution against applying the model there. A separate exposure-response cohort",
      "of 161 patients (Table S2, after excluding 4 who discontinued for adverse effects) had a",
      "median steady-state AUC(0-24 h) of 35.59 mg*h/L (IQR 33.11-41.54) on a median daily dose of",
      "10.00 mg/kg, all of whom achieved clinical cure or improvement. Model estimation used",
      "Phoenix NLME 8.5 with first-order conditional estimation-extended least squares."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural disposition. Tao 2026 Table 2, 'Final model estimate'
    # column. The typical values are referenced to the cohort MEAN body
    # weight (24.05 kg) and MEAN serum creatinine (31.77 umol/L), which
    # are the normalizing constants printed in the Table 2 footnote
    # equations:
    #   CL = theta_CL * (BW/24.05)^theta_BW_CL * (SCr/31.77)^theta_SCr_CL * exp(eta_CL)
    #   V1 = theta_V1 * (BW/24.05)^theta_BW_V1
    # A reader who reached for the Table 1 medians (22.50 kg, 30.70
    # umol/L) would shift every typical value; the equations are explicit
    # that the means are the reference. The bracketed bootstrap medians
    # below are the Table 2 right-hand column and agree with the point
    # estimates throughout, which is the paper's robustness evidence.
    # ------------------------------------------------------------------
    lcl <- log(6.85); label("Typical clearance at 24.05 kg body weight and 31.77 umol/L serum creatinine (L/h)")  # Table 2 theta_CL = 6.85 L/h (RSE 4.33%); bootstrap median 6.88, 95% CI 6.34-7.48. Weight-normalized this is 6.85/24.05 = 0.285 L/kg/h, the 0.28 L/kg/h quoted in the Discussion.
    lvc <- log(25.04); label("Typical central volume of distribution at 24.05 kg body weight (L)")                # Table 2 theta_V1 = 25.04 L (RSE 8.78%); bootstrap median 25.15, 95% CI 20.98-29.94
    lq <- log(1.20); label("Intercompartmental clearance (L/h)")                                                  # Table 2 theta_Q = 1.20 L/h (RSE 29.54%); bootstrap median 1.19, 95% CI 0.64-2.72. No covariate enters Q.
    lvp <- log(5.28); label("Peripheral volume of distribution (L)")                                              # Table 2 theta_V2 = 5.28 L (RSE 18.84%); bootstrap median 5.32, 95% CI 3.51-8.38. No covariate enters V2.

    # ------------------------------------------------------------------
    # Covariate exponents, Tao 2026 Table 2. Both covariates enter as
    # power functions (Results: 'the relationships best described by
    # power models'; 'SCr was also identified as a significant covariate
    # for CL, likewise modeled using a power function').
    #
    # The body-weight exponents are ESTIMATED, not fixed at the
    # conventional allometric 0.75 / 1, and both land well above them
    # (1.11 on clearance, 1.83 on central volume). The confidence
    # intervals exclude the allometric values for V1 (bootstrap 95% CI
    # 1.24-2.18) but not for CL (0.76-1.29). A consequence worth knowing
    # before simulating outside the observed 9-45 kg window: with an
    # exponent of 1.83 the central volume grows faster than body weight,
    # so weight-normalized V1 rises with size rather than staying flat.
    # ------------------------------------------------------------------
    e_wt_cl <- 1.11; label("Power exponent on (WT / 24.05 kg) for CL (unitless)")           # Table 2 theta_BW_CL = 1.11 (RSE 12.88%); bootstrap median 1.09, 95% CI 0.76-1.29
    e_creat_cl <- -0.20; label("Power exponent on (CREAT / 31.77 umol/L) for CL (unitless)") # Table 2 theta_SCr_CL = -0.20 (RSE 24.86%); bootstrap median -0.20, 95% CI -0.31 to -0.10
    e_wt_vc <- 1.83; label("Power exponent on (WT / 24.05 kg) for Vc (unitless)")            # Table 2 theta_BW_V1 = 1.83 (RSE 14.66%); bootstrap median 1.77, 95% CI 1.24-2.18

    # ------------------------------------------------------------------
    # Interindividual variability, exponential (Methods: 'Interindividual
    # variability of PK parameters was modeled exponentially').
    #
    # VARIANCE, NOT SD. The Table 2 row is headed 'CL (omega 2)', i.e. the
    # table prints omega-squared directly, so 0.013 is a variance on the
    # log scale and needs no conversion. The arithmetic corroborates it:
    # sqrt(0.013) = 0.114, an 11.4 percent coefficient of variation, which
    # is small but consistent with a model whose two covariates absorb most
    # of the between-subject signal. The alternative reading, 0.013 as a
    # standard deviation, would give a 1.3 percent CV, i.e. essentially no
    # between-subject variability at all, and would make the Table 3 Monte
    # Carlo target attainment step almost vertically from 0 to 100 percent
    # rather than spanning 14.6-88.8 percent across the dosing regimens.
    #
    # Only one eta exists. Results: 'Interindividual variability for the
    # central (V1) and peripheral (V2) volumes of distribution and
    # intercompartmental clearance was fixed to zero due to high
    # shrinkage.' The sparse opportunistic design (293 samples from 191
    # patients) is the reason. Those three etas are omitted rather than
    # written as `~ fixed(0)` because a zero-variance diagonal makes OMEGA
    # singular and breaks the Cholesky sampler rxSolve uses.
    # ------------------------------------------------------------------
    etalcl ~ 0.013  # Table 2 interindividual variability, row 'CL (omega 2)' = 0.013 (RSE 21.08%); bootstrap median 0.012, 95% CI 0.0073-0.018. sqrt(0.013) = 0.114 on the log scale.

    # ------------------------------------------------------------------
    # Residual unexplained variability. Additive, proportional and
    # combined error models were all evaluated and 'residual variability
    # was best captured by a proportional model' (Results). Table 2 prints
    # it as a percentage, 35.0, which is the proportional SD as a fraction
    # of the prediction.
    # ------------------------------------------------------------------
    propSd <- 0.350; label("Proportional residual error (fraction)")  # Table 2 residual variability, 'Proportional error (%)' = 35.0 (RSE 5.90%); bootstrap median 35.0, 95% CI 31.0-39.0
  })

  model({
    # 1. Individual parameters. Covariate model exactly as printed in the
    #    Tao 2026 Table 2 footnote:
    #      CL = theta_CL * (BW/24.05)^theta_BW_CL * (SCr/31.77)^theta_SCr_CL * exp(eta_CL)
    #      V1 = theta_V1 * (BW/24.05)^theta_BW_V1
    #    Q and V2 carry no covariate and no eta.
    cl <- exp(lcl + etalcl) * (WT / 24.05)^e_wt_cl * (CREAT / 31.77)^e_creat_cl
    vc <- exp(lvc) * (WT / 24.05)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Two compartments with first-order elimination from
    #    the central compartment (Results: 'A two-compartment model with
    #    first-order elimination best described the data'). Levofloxacin
    #    was given intravenously by syringe-pump infusion over 0.5-2 h, so
    #    doses go to `central` with a rate or duration; there is no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation and error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
